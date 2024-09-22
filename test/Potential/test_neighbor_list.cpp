#include "copyright.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/formulas.h"
#include "../../src/Math/log_scale_spline.h"
#include "../../src/Math/math_enumerators.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Numerics/host_popc.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/cellgrid.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/forward_exclusionmask.h"
#include "../../src/Potential/local_exclusionmask.h"
#include "../../src/Potential/pme_util.h"
#include "../../src/Potential/pmigrid.h"
#include "../../src/Potential/ppitable.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/tile_manager.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/local_arrangement.h"
#include "../../src/Structure/structure_enumerators.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/brickwork.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/synthesis_abstracts.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

#ifndef STORMM_USE_HPC
using stormm::int2;
using stormm::int3;
using stormm::int4;
using stormm::double2;
using stormm::double3;
using stormm::double4;
using stormm::float4;
using stormm::short4;
using stormm::uint2;
#endif
using stormm::llint;
using stormm::ullint;
using stormm::Ecumenical4;
using stormm::Ecumenical8;
using stormm::double_type_index;
using stormm::double4_type_index;
using stormm::float_type_index;
using stormm::int_type_index;
using stormm::llint_type_index;
using stormm::constants::ExceptionResponse;
using stormm::constants::warp_size_int;
using stormm::constants::warp_bits_mask_int;
using stormm::data_types::int95_t;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getBaseName;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::errors::rtWarn;
using stormm::numerics::hostPopcll;
using stormm::parse::addLeadingWhiteSpace;
using stormm::parse::char4ToString;
using stormm::parse::findAlignmentWidth;
using stormm::parse::lowercase;
using stormm::parse::NumberFormat;
using stormm::parse::polyNumericVector;
using stormm::parse::uppercase;
using stormm::random::Xoshiro256ppGenerator;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using stormm::stmath::BasisFunctions;
using stormm::stmath::ipow;
using stormm::stmath::LogScaleSpline;
using stormm::stmath::LogSplineForm;
using stormm::stmath::reduceUniqueValues;
using stormm::stmath::tileVector;
using stormm::synthesis::AtomGraphSynthesis;
using stormm::synthesis::Brickwork;
using stormm::synthesis::PhaseSpaceSynthesis;
using stormm::synthesis::SyNonbondedKit;
using stormm::topology::AtomGraph;
using stormm::topology::ChemicalDetailsKit;
using stormm::topology::NonbondedKit;
using stormm::trajectory::CoordinateFileKind;
using stormm::trajectory::CoordinateFrame;
using stormm::trajectory::CoordinateFrameReader;
using stormm::trajectory::PhaseSpace;
using stormm::trajectory::TrajectoryKind;
using namespace stormm::card;
using namespace stormm::energy;
using namespace stormm::structure;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Loop over all exclusions in a given array based on the supplied bounds.  Check for non-reflexive
// exclusions and also for self-exclusions (neither of which should exist).  Modify formal
// arguments if either condition is found.  Return a vector of booleans to indicate that exclusions
// listed in the topology were found within a compact "forward" exclusion mask.
//
// Arguments
//   natom:                     Number of atoms in the system
//   bounds:                    Bounds array for exclusions in excl_list
//   excl_list:                 Exclusions list, rather like EXCLUDED_ATOMS_LIST in an Amber prmtop
//                              but with double-counting for all exclusions
//   femask:                    A very compact expression of the system's masks
//   self_exclusion_detected:   Indicator of self exclusions detected in error (returned)
//   non_reflexive_exclusions:  Indicator of non-reflexive elements of excl_list (returned)
//-------------------------------------------------------------------------------------------------
std::vector<bool> checkForwardMaskExcl(const int natom, const int* bounds, const int* excl_list,
                                       const ForwardExclusionMask &femask,
                                       bool *self_exclusion_detected,
                                       std::vector<int2> *non_reflexive_exclusions) {
  std::vector<bool> result(bounds[natom], false);
  for (int i = 0; i < natom; i++) {
    for (int j = bounds[i]; j < bounds[i + 1]; j++) {
      const int j_atom_idx = excl_list[j];
      *self_exclusion_detected = (*self_exclusion_detected || j_atom_idx == i);
      bool mirror_found = false;
      for (int k = bounds[j_atom_idx]; k < bounds[j_atom_idx + 1]; k++) {
        mirror_found = (mirror_found || excl_list[k] == i);
      }
      if (mirror_found == false) {
        non_reflexive_exclusions->push_back({i, j_atom_idx});
      }
      result[j] = femask.testExclusion(i, j_atom_idx);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Check the placement of atoms within a cell grid.
//
// Arguments:
//   poly_ps:   A synthesis of coordinates, containing multiple systems which are all represented
//              in the cell grid cg
//   sysno:     Index of a system to examine within the coordinate synthesis and thus the cell grid
//   cg:        The cell grid holding all atoms of the larger synthesis
//   tol:       Tolerance by which coordinates extracted from the cell grid should match their
//              original counterparts after double-precision re-imaging
//   do_tests:  Indicate whether tests can be run
//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void checkCellGridPlacements(const PhaseSpaceSynthesis &poly_ps, const int sysno,
                             const CellGrid<T, Tacc, Tcalc, T4> &cg, const double tol,
                             const TestPriority do_tests) {
  const CoordinateFrame cf = poly_ps.exportCoordinates(sysno);
  const int sys_atom_offset = poly_ps.getAtomOffset(sysno);
  const CoordinateFrameReader cfr = cf.data();
  const AtomGraph *ag_ptr = poly_ps.getSystemTopologyPointer(sysno);
  const NonbondedKit<double> nbk = ag_ptr->getDoublePrecisionNonbondedKit();
  const CellGridReader<T, Tacc, Tcalc, T4> cgr = cg.data();
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool fpconv = (isFloatingPointScalarType<T>() == false);
  const T value_zero = 0.0;

  // Extract the nonbonded parameters from the attached topology synthesis.
  const SyNonbondedKit<double, double2> synbk =
    cg.getTopologySynthesisPointer()->getDoublePrecisionNonbondedKit();
  
  // Loop through the system's chains, visiting all cells within each chain and checking the atom
  // placements within it one by one.
  const ullint cdims = cgr.system_cell_grids[sysno];
  const int na = ((cdims >> 28) & 0xfffLLU);
  const int nb = ((cdims >> 40) & 0xfffLLU);
  const int nc = ((cdims >> 52) & 0xfffLLU);
  const int sysc_offset = (cdims & 0xfffffffLLU);
  const double d_na = na;
  const double d_nb = nb;
  const double d_nc = nc;

  // Check that the atoms are arranged as expected
  bool back_to_back = true;
  bool big_space = true;
  for (int k = 0; k < nc; k++) {
    for (int j = 0; j < nb; j++) {
      for (int i = 0; i < na; i++) {
        const uint2 cijk_bounds = cgr.cell_limits[sysc_offset + (((k * nb) + j) * na) + i];
        if (i < na - 1) {
          const uint2 cijk_next = cgr.cell_limits[sysc_offset + (((k * nb) + j) * na) + i + 1];
          back_to_back = (back_to_back &&
                          (cijk_bounds.x + ((cijk_bounds.y >> 16) & 0xffff) == cijk_next.x));
        }
        else {
          if (k < nc - 1 || j < nb - 1) {
            const uint2 cijk_next = cgr.cell_limits[sysc_offset + (((k * nb) + j) * na) + i + 1];
            const size_t pooled_spc = cijk_next.x -
                                      (cijk_bounds.x + ((cijk_bounds.y >> 16) & 0xffff));

            // While it is not required that cell grids place a full cell's maximum capacity as
            // space after the end of each chain, it is very likely to happen and in all of these
            // cases is is definitely expected.
            big_space = (big_space && pooled_spc > cgr.cell_base_capacity);
          }
        }
      }
    }
  }
  check(back_to_back, "Cells in system " + getBaseName(ag_ptr->getFileName()) + " do not obey the "
        "'back to back' rule in terms of their atom placement in the chains.\n", do_tests);
  check(big_space, "Cells in system " + getBaseName(ag_ptr->getFileName()) + " do not obey the "
        "'big space' rule in terms of their atom placement in the chains.\n", do_tests);

  // Check individual atom properties and coordinates
  std::vector<double> x_img(cfr.natom), y_img(cfr.natom), z_img(cfr.natom);
  std::vector<double> x_chk(cfr.natom), y_chk(cfr.natom), z_chk(cfr.natom);
  std::vector<int> presence(cfr.natom, 0);
  std::vector<int> img_lj_type_indices(cfr.natom), chk_lj_type_indices(cfr.natom);
  std::vector<double> img_charges(cfr.natom), chk_charges(cfr.natom);
  int n_chain_cell_failures = 0;
  int pos = 0;
  for (int i = 0; i < na; i++) {
    const double di = static_cast<double>(i) / d_na;
    for (int j = 0; j < nb; j++) {
      const double dj = static_cast<double>(j) / d_nb;
      for (int k = 0; k < nc; k++) {
        const double dk = static_cast<double>(k) / d_nc;
        const double orig_x = (cfr.invu[0] * di) + (cfr.invu[3] * dj) + (cfr.invu[6] * dk);
        const double orig_y =                      (cfr.invu[4] * dj) + (cfr.invu[7] * dk);
        const double orig_z =                                           (cfr.invu[8] * dk);
        const uint2 cijk_bounds = cgr.cell_limits[sysc_offset + (((k * nb) + j) * na) + i];
        const uint mhlim = cijk_bounds.x + ((cijk_bounds.y >> 16) & 0xffff);
        for (uint m = cijk_bounds.x; m < mhlim; m++) {

          // Check that the atom's chain cell identification is consistent with the cell that
          // it is actually in.
          if (cgr.img_atom_chn_cell[m] != i) {
            n_chain_cell_failures++;
          }

          // Extract the atom's position / property tuple and return it to the broader system
          // context.
          const T4 atom_tuple = cgr.image[m];
          double img_atom_x, img_atom_y, img_atom_z;
          if (fpconv) {
            img_atom_x = llround(atom_tuple.x * cgr.lpos_scale);
            img_atom_y = llround(atom_tuple.y * cgr.lpos_scale);
            img_atom_z = llround(atom_tuple.z * cgr.lpos_scale);
          }
          else {
            img_atom_x = atom_tuple.x;
            img_atom_y = atom_tuple.y;
            img_atom_z = atom_tuple.z;
          }
          x_img[pos] = img_atom_x + orig_x;
          y_img[pos] = img_atom_y + orig_y;
          z_img[pos] = img_atom_z + orig_z;

          // Pull out the atom properties
          switch (cg.getTheme()) {
          case NonbondedTheme::ELECTROSTATIC:

            // The charge is represented as a real number, if possible, or reinterpreted bitwise
            // as a signed integer of the appropriate size.
            if (fpconv) {
              if (ct == llint_type_index) {
                const Ecumenical8 qconv = { .lli = static_cast<llint>(atom_tuple.w) };
                img_charges[pos] = qconv.d;
              }
              else if (ct == int_type_index) {
                const Ecumenical4 qconv = { .i = static_cast<int>(atom_tuple.w) };
                img_charges[pos] = qconv.f;
              }
              else {

                // This case should never be reached, but checks on the validity of the image data
                // type occur elsewhere.
                img_charges[pos] = value_zero;
              }
            }
            else {
              img_charges[pos] = atom_tuple.w;
            }
            break;
          case NonbondedTheme::VAN_DER_WAALS:

            // The Lennard-Jones index is represented as an integer, or reinterpreted bitwise as a
            // floating point number of the appropriate size.
            if (fpconv) {
              img_lj_type_indices[pos] = atom_tuple.w;
            }
            else {
              if (ct == double_type_index) {
                const Ecumenical8 tconv = { .d = static_cast<double>(atom_tuple.w) };
                img_lj_type_indices[pos] = tconv.lli;
              }
              else if (ct == float_type_index) {
                const Ecumenical4 tconv = { .f = static_cast<float>(atom_tuple.w) };
                img_lj_type_indices[pos] = tconv.i;
              }
              else {
                img_lj_type_indices[pos] = value_zero;
              }
            }
            break;
          case NonbondedTheme::ALL:
            
            // The combined charge and Lennard-Jones indices are packed into the low and high bits
            // of an integer, respectively.  The demarcation varies by how many bits are available
            // to work with in the image data type.
            if (ct == int_type_index) {
              const Ecumenical4 tconv = { .i = static_cast<int>(atom_tuple.w) };
              img_lj_type_indices[pos] = ((tconv.ui >> sp_charge_index_bits) & 0x1ff);
              const int img_q_idx = (tconv.ui & 0x7fffff);
              img_charges[pos] = synbk.q_params[img_q_idx];
            }
            else if (ct == llint_type_index) {
              img_lj_type_indices[pos] = ((static_cast<llint>(atom_tuple.w) >>
                                           dp_charge_index_bits) & 0xffffffffLL);
              const int img_q_idx = (static_cast<llint>(atom_tuple.w) & 0xffffffffLL);
              img_charges[pos] = synbk.q_params[img_q_idx];
            }
            else if (ct == float_type_index) {
              const Ecumenical4 tconv = { .f = static_cast<float>(atom_tuple.w) };
              img_lj_type_indices[pos] = ((tconv.ui >> sp_charge_index_bits) & 0x1ff);
              const int img_q_idx = (tconv.ui & 0x7fffff);
              img_charges[pos] = synbk.q_params[img_q_idx];
            }
            else if (ct == double_type_index) {
              
              // Because the bit string had to be converted anyway, it was taken as an unsigned
              // long long int in its inception.  Reading as a signed long long int would probably
              // be safe as well: there are not enough Lennard-Jones types to overflow the format.
              const Ecumenical8 tconv = { .d = static_cast<double>(atom_tuple.w) };
              img_lj_type_indices[pos] = ((tconv.ulli >> dp_charge_index_bits) & 0xffffffffLLU);
              const int img_q_idx = (tconv.ulli & 0xffffffffLLU);
              img_charges[pos] = synbk.q_params[img_q_idx];
            }
            else {
              img_charges[pos] = 0.0;
              img_lj_type_indices[pos] = 0;
            }
            break;
          }
          
          // Find the atom in the original structure and re-image it.  The two sets of coordinates
          // should match.
          const int topological_idx = cgr.nonimg_atom_idx[m] - sys_atom_offset;
          double atom_x = cfr.xcrd[topological_idx];
          double atom_y = cfr.ycrd[topological_idx];
          double atom_z = cfr.zcrd[topological_idx];
          imageCoordinates<double, double>(&atom_x, &atom_y, &atom_z, cfr.umat, cfr.invu,
                                           cfr.unit_cell, ImagingMethod::PRIMARY_UNIT_CELL);
          presence[topological_idx] += 1;
          x_chk[pos] = atom_x;
          y_chk[pos] = atom_y;
          z_chk[pos] = atom_z;
          chk_charges[pos] = nbk.charge[topological_idx];
          chk_lj_type_indices[pos] = nbk.lj_idx[topological_idx];          
          pos++;
        }
      }
    }
  }
  check(x_img, RelationalOperator::EQUAL, Approx(x_chk).margin(tol), "Cartesian X coordinates of "
        "particles in the cell grid image do not match their counterparts imaged from the "
        "coordinate frame.  Synthesis structure index: " + std::to_string(sysno) + ".", do_tests);
  check(y_img, RelationalOperator::EQUAL, Approx(y_chk).margin(tol), "Cartesian X coordinates of "
        "particles in the cell grid image do not match their counterparts imaged from the "
        "coordinate frame.  Synthesis structure index: " + std::to_string(sysno) + ".", do_tests);
  check(z_img, RelationalOperator::EQUAL, Approx(z_chk).margin(tol), "Cartesian X coordinates of "
        "particles in the cell grid image do not match their counterparts imaged from the "
        "coordinate frame.  Synthesis structure index: " + std::to_string(sysno) + ".", do_tests);
  check(presence, RelationalOperator::EQUAL, std::vector<int>(cfr.natom, 1), "Various atoms were "
        "over- or under-represented in the cell grid for system " + std::to_string(sysno) + ".",
        do_tests);
  check(n_chain_cell_failures, RelationalOperator::EQUAL, 0, "The chain cells of some atoms were "
        "not properly recorded.", do_tests);
  bool perform_q_check = false;
  bool perform_lj_check = false;
  switch (cg.getTheme()) {
  case NonbondedTheme::ELECTROSTATIC:
    perform_q_check = true;
    break;
  case NonbondedTheme::VAN_DER_WAALS:
    perform_lj_check = true;
    break;
  case NonbondedTheme::ALL:
    perform_q_check = true;
    perform_lj_check = true;
    break;
  }
  if (perform_q_check) {
    check(img_charges, RelationalOperator::EQUAL, Approx(chk_charges).margin(tol), "Charges "
          "presented in the image do not correspond to the topological values.  Non-bonded "
          "theme: " + getEnumerationName(cg.getTheme()) + ".", do_tests);
  }
  if (perform_lj_check) {
    check(img_lj_type_indices, RelationalOperator::EQUAL, Approx(chk_lj_type_indices).margin(tol),
          "Lennard-Jones type indices presented in the image do not correspond to the topological "
          "values.  Non-bonded theme: " + getEnumerationName(cg.getTheme()) + ".", do_tests);
  }

  // Check that the functions for finding the chain and cell locations of particular atoms are
  // working properly.
  int n_image_func_fail = 0;
  int n_chain_func_fail = 0;
  int n_cell_func_fail  = 0;
  int n_cell_xyz_fail = 0;
  int total_atoms = 0;
  for (int i = 0; i < synbk.nsys; i++) {
    const int j_offset = synbk.atom_offsets[i];
    const ullint i_gdims = cgr.system_cell_grids[i];
    const int cell_offset = (i_gdims & 0xfffffff);
    const int cell_na = ((i_gdims >> 28) & 0xfff);
    const int cell_nb = ((i_gdims >> 40) & 0xfff);
    for (int j = 0; j < synbk.atom_counts[i]; j++) {
      const uint image_loc = cg.getImageIndex(i, j);
      n_image_func_fail += (cgr.nonimg_atom_idx[image_loc] != j_offset + j);
      const int2 chain_loc = cg.getChainLocation(i, j);
      n_chain_func_fail += (cgr.chain_limits[chain_loc.x] + chain_loc.y != image_loc);
      const int4 cell_loc = cg.getCellLocation(i, j);
      const uint2 ij_cell_bounds = cgr.cell_limits[cell_loc.w];
      n_cell_func_fail += (image_loc < ij_cell_bounds.x ||
                           image_loc >= ij_cell_bounds.x + (ij_cell_bounds.y >> 16));
      n_cell_xyz_fail += ((((cell_loc.z * cell_nb) + cell_loc.y) * cell_na) + cell_loc.x +
                          cell_offset != cell_loc.w);
    }
    total_atoms += synbk.atom_counts[i];
  }
  check(n_image_func_fail == 0, "The CellGrid image indexing member function fails to identify " +
        std::to_string(n_image_func_fail) + " out of " + std::to_string(total_atoms) + " atoms.",
        do_tests);
  check(n_chain_func_fail == 0, "The CellGrid chain indexing member function fails to pinpoint " +
        std::to_string(n_chain_func_fail) + " out of " + std::to_string(total_atoms) + " atoms.",
        do_tests);
  check(n_cell_func_fail == 0, "The CellGrid cell indexing member function fails to place " +
        std::to_string(n_cell_func_fail) + " out of " + std::to_string(total_atoms) + " atoms.",
        do_tests);
  check(n_cell_xyz_fail == 0, "The CellGrid cell indexing member function mismatches the A, B, "
        "and C system-specific cell indices of " + std::to_string(n_cell_xyz_fail) + " out of " +
        std::to_string(total_atoms) + " atoms with the overall cell indices.", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Test the various components of a profile for collisions.
//
// Arguments:
//   mode_name:    The atom profile mode (see the LocalExclusionMask class)
//   components:   Array of different bitmasks delimiting the extent of each profile component
//   descriptors:  Array of names for each component
//-------------------------------------------------------------------------------------------------
void screenProfileCollisions(const std::string &mode_name, const std::vector<ullint> &components,
                             const std::vector<std::string> &descriptors,
                             const std::vector<int> &component_lengths) {
  const size_t ncomp = components.size();
  for (size_t i = 1; i < ncomp; i++) {
    for (size_t j = 0; j < i; j++) {
      check((components[i] & components[j]), RelationalOperator::EQUAL, 0LLU, "Components '" +
            descriptors[i] + "' and '" + descriptors[j] + "' of mask mode " + mode_name +
            " overlap in " + std::to_string(hostPopcll(components[i] & components[j])) + " bits.");
    }

    // The mode mask is expected to be the first component.  Its extent is checked in the
    // testLocalExclusionMaskSizing() function below.
    check(hostPopcll(components[i]), RelationalOperator::EQUAL, component_lengths[i],
          "Component " + descriptors[i] + " does not have the correct nuber of bits.");
  }
}

//-------------------------------------------------------------------------------------------------
// Test bit packing in the LocalExclusionMask object.  These tests will ensure that edits and
// adjustments to the limits always fit the 64-bit format, do not have different parts of the
// profiles colliding, and can deliver the proper information when implemented.
//-------------------------------------------------------------------------------------------------
void testLocalExclusionMaskSizing() {

  // General dimensions
  check(hostPopcll(lmask_mode_bitmask), RelationalOperator::EQUAL, lmask_mode_bits, "The mode "
        "selector bitmask does not have the correct number of bits.");
  const std::vector<ullint> all_modes = { lmask_mode_a, lmask_mode_b, lmask_mode_c, lmask_mode_d,
                                          lmask_mode_e, lmask_mode_f };
  const std::vector<std::string> mode_names = { "A", "B", "C", "D", "E", "G" };
  for (size_t i = 0; i < all_modes.size(); i++) {
    check((all_modes[i] & lmask_mode_bitmask), RelationalOperator::EQUAL, all_modes[i],
          "The mode " + mode_names[i] + " selector is invalid by the allotted mode bits.");
  }

  // Mode A  
  check((2 * lmask_long_local_span) + 1 + lmask_mode_bits, RelationalOperator::LESS_THAN_OR_EQUAL,
        64, "The number of bits in the A mode mask exceeds the 64-bit integer format.");
  const std::vector<ullint> a_masks = { lmask_mode_bitmask, lmask_a_excl };
  const std::vector<std::string> a_mask_parts = { "mode mask", "local exclusion mask" };
  const std::vector<int> a_mask_lengths = { lmask_mode_bits, (2 * lmask_long_local_span) + 1 };
  screenProfileCollisions(mode_names[0], a_masks, a_mask_parts, a_mask_lengths);

  // Mode B
  check((2 * (lmask_short_local_span + lmask_short_extra_span + lmask_b_shft_bits)) + 1 +
        lmask_mode_bits, RelationalOperator::LESS_THAN_OR_EQUAL, 64, "The number of bits in the "
        "B mode mask exceeds the 64-bit integer format.");
  const std::vector<ullint> b_masks = { lmask_mode_bitmask, lmask_b_excl, lmask_b_lower_mask,
                                        lmask_b_upper_mask, lmask_b_lower_shft,
                                        lmask_b_upper_shft };
  const std::vector<std::string> b_mask_parts = { "mode mask", "local exclusion mask",
                                                  "lower extension mask", "upper extension mask",
                                                  "lower extension shift",
                                                  "upper extension shift" };
  const std::vector<int> b_mask_lengths = { lmask_mode_bits, (2 * lmask_short_local_span) + 1,
                                            lmask_short_extra_span, lmask_short_extra_span,
                                            lmask_b_shft_bits, lmask_b_shft_bits };
  screenProfileCollisions(mode_names[1], b_masks, b_mask_parts, b_mask_lengths);

  // Mode C
  check((2 * lmask_short_local_span) + 1 + lmask_c_alt_mask_bits + lmask_c_shft_bits +
        lmask_mode_bits, RelationalOperator::LESS_THAN_OR_EQUAL, 64, "The number of bits in the "
        "C mode mask exceeds the 64-bit integer format.");
  const std::vector<ullint> c_masks = { lmask_mode_bitmask, lmask_c_excl, lmask_c_alt_mask,
                                        lmask_c_shft };
  const std::vector<std::string> c_mask_parts = { "mode mask", "local exclusion mask",
                                                  "alternate mask", "large index shift" };
  const std::vector<int> c_mask_lengths = { lmask_mode_bits, (2 * lmask_short_local_span) + 1,
                                            lmask_long_extra_span, lmask_c_shft_bits };
  screenProfileCollisions(mode_names[2], c_masks, c_mask_parts, c_mask_lengths);

  // Mode D 
  check((2 * lmask_short_local_span) + 1 + lmask_d_array_idx_bits + lmask_d_array_cnt_bits +
        lmask_mode_bits, RelationalOperator::LESS_THAN_OR_EQUAL, 64, "The number of bits in the "
        "E mode mask exceeds the 64-bit integer format.");
  const std::vector<ullint> d_masks = { lmask_mode_bitmask, lmask_d_excl, lmask_d_array_idx,
                                        lmask_d_array_cnt };
  const std::vector<std::string> d_mask_parts = { "mode mask", "local exclusion mask",
                                                  "extra mask lower index bound",
                                                  "extra mask count" };
  const std::vector<int> d_mask_lengths = { lmask_mode_bits, (2 * lmask_short_local_span) + 1,
                                            lmask_d_array_idx_bits, lmask_d_array_cnt_bits };
    
  screenProfileCollisions(mode_names[4], d_masks, d_mask_parts, d_mask_lengths);

  // Mode E 
  check(lmask_e_array_idx_bits + lmask_e_array_cnt_bits + lmask_mode_bits,
        RelationalOperator::LESS_THAN_OR_EQUAL, 64, "The number of bits in the F mode mask "
        "exceeds the 64-bit integer format.");
  const std::vector<ullint> e_masks = { lmask_mode_bitmask, lmask_e_array_idx,
                                        lmask_e_array_cnt };
  const std::vector<std::string> e_mask_parts = { "mode mask", "extra mask lower index bound",
                                                  "extra mask count" };
  const std::vector<int> e_mask_lengths = { lmask_mode_bits, lmask_e_array_idx_bits,
                                            lmask_e_array_cnt_bits };
  screenProfileCollisions(mode_names[5], e_masks, e_mask_parts, e_mask_lengths);
}

//-------------------------------------------------------------------------------------------------
// Print a human-readable list of exclusions and the index atom for the purposes of investigating
// a failed test.
//
// Arguments:
//   index_atom:  The index atom
//   exclusions:  List of excluded atoms
//-------------------------------------------------------------------------------------------------
std::string printIndexAndExclusions(const int index_atom, const std::vector<int> &exclusions) {
  std::string result = "Index atom: " + std::to_string(index_atom) + "\n";
  result += "Exclusions: ";
  const int nexcl = exclusions.size();
  int lpos = 0;
  const int nalign = findAlignmentWidth(exclusions, 0);
  int n_printed = 0;
  bool snooze_reports = false;
  for (int i = 0; i < nexcl; i++) {
    if (n_printed >= 64 && i < nexcl - 4) {
      if (snooze_reports == false) {
        result += " ...\n";
        lpos = 0;
        snooze_reports = true;
      }
      continue;
    }
    result += addLeadingWhiteSpace(std::to_string(exclusions[i]), nalign);
    lpos++;
    if (lpos == 8) {
      lpos = 0;
      result += "\n            ";
    }
    else if (i < nexcl - 1) {
      result += ", ";
    }
  }
  if (lpos > 0) {
    result += "\n";
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Run repeated trials, forming a list of exclusions within a set list of boundaries, to see
// whether each exclusion list will fit within the specified mode.  Return a string containing an
// error report if any tests failed.  A zero length in the resulting string will be taken as an
// indicator of success.
//
// Arguments:
//   xrs:      Random number generator to differentiate sequences and permit repeated tests
//   limits:   A series of limits as to where the exclusions can occur relative to the index atom
//             (the topological position of the index atom will be randomly selected).  The lower
//             limit is placed in the "x" member of the tuple and the upper limit is placed in the
//             "y" member.
//   mode:     Character code for the mask mode
//   ntrials:  The number of trials to perform
//-------------------------------------------------------------------------------------------------
void testExclusionZones(Xoshiro256ppGenerator *xrs, const std::vector<int2> &limits,
                        const char mode, const int ntrials = 10) {
  std::string failed_results;
  std::string mistaken_passes;
  const int n_limits = limits.size();
  int min_index = limits[0].y;
  for (int i = 0; i < n_limits; i++) {
    min_index = std::min(min_index, limits[i].x);
  }
  min_index = abs(min_index);
  for (int trial = 0; trial < ntrials; trial++) {
    const int index_atom = static_cast<int>(1000.0 * xrs->uniformRandomNumber()) + min_index;
    std::vector<uint2> tmp_secondary_masks;

    // Test whether examples of exclusions that should pass fit the profile format.
    std::vector<int> exclusions;
    for (int i = 0; i < n_limits; i++) {
      for (int j = index_atom + limits[i].x; j <= index_atom + limits[i].y;
           j += static_cast<int>(1.0 + (3.0 * xrs->uniformRandomNumber()))) {
        exclusions.push_back(j);
      }
    }
    switch (mode) {
    case 'A':
      if (lMaskFitsModeA(index_atom, exclusions) == false) {
        failed_results += "\n" + printIndexAndExclusions(index_atom, exclusions);        
      }
      break;
    case 'B':
    case 'C':
    case 'D':
    case 'E':
      break;
    default:
      rtErr("Unrecognized mode " + std::string(1, mode) + ".", "testExclusionZones");
    }

    // Some cases involve variable shifts, which are implicit in the supplied limits.  However,
    // the format might accept a new "disruptive" atom if it is near enough to the stated limits
    // but there are blanks on the other side the span such that the bounds could be moved, within
    // the limits of the format, to yet conform to the format.  In order to ensure that added atoms
    // outside the stated limits will break the format, place atoms at the limits.
    switch (mode) {
    case 'A':
    case 'E':

      // Mode A operates under strict limits, with no variable extensions.  Mode E is designed to
      // be unbreakable, at least within the scope of a molecular system that fits in computer
      // memory).
      break;
    case 'B':
    case 'C':
    case 'D':
      for (int i = 0; i < n_limits; i++) {
        exclusions.push_back(index_atom + limits[i].x);
        exclusions.push_back(index_atom + limits[i].y);
      }
      reduceUniqueValues(&exclusions);
      if (mode == 'B' && lMaskFitsModeB(index_atom, exclusions) == false) {
        failed_results += "\n" + printIndexAndExclusions(index_atom, exclusions);        
      }
      else if (mode == 'C' && lMaskFitsModeC(index_atom, exclusions) == false) {
        failed_results += "\n" + printIndexAndExclusions(index_atom, exclusions);        
      }
      else if (mode == 'D' &&
               lMaskFitsModeD(index_atom, exclusions, tmp_secondary_masks) == false) {
        failed_results += "\n" + printIndexAndExclusions(index_atom, exclusions);        
      }
      break;
    default:
      break;
    }
    
    // Test the inverse case, creating lists of exclusions that should violate the format in one
    // way or another, to ensure that they do NOT pass.
    int disruptor;
    bool non_disruptive;
    do {
      disruptor = index_atom + static_cast<int>(500.0 - (1000.0 * xrs->uniformRandomNumber()));
      non_disruptive = false;
      for (int i = 0; i < n_limits; i++) {
        non_disruptive = (non_disruptive ||
                          (disruptor >= index_atom + limits[i].x &&
                           disruptor <= index_atom + limits[i].y));
      }

      // Do not admit negative numbers to the exclusions list
      non_disruptive = (non_disruptive || disruptor < 0);
    } while (non_disruptive);
    exclusions.push_back(disruptor);
    reduceUniqueValues(&exclusions);
    switch (mode) {
    case 'A':
      if (lMaskFitsModeA(index_atom, exclusions)) {
        mistaken_passes += "\n" + printIndexAndExclusions(index_atom, exclusions);        
      }
      break;
    case 'B':
      if (lMaskFitsModeB(index_atom, exclusions)) {
        mistaken_passes += "\n" + printIndexAndExclusions(index_atom, exclusions);
      }
      break;
    case 'C':
      if (lMaskFitsModeC(index_atom, exclusions)) {
        mistaken_passes += "\n" + printIndexAndExclusions(index_atom, exclusions);        
      }
      break;
    case 'D':
      if (lMaskFitsModeD(index_atom, exclusions, tmp_secondary_masks)) {
        mistaken_passes += "\n" + printIndexAndExclusions(index_atom, exclusions);        
      }
      break;
    case 'E':
      break;
    }
  }
  check(failed_results.size() == 0, "One or more lists of exclusions did not fit Mode " +
        std::string(1, mode) + " as expected.\n" + failed_results);

  // Only check for mistaken list admissions in limited formats
  switch (mode) {
  case 'A':
  case 'B':
  case 'C':
    check(mistaken_passes.size() == 0, "One or more lists of exclusions was found to fit Mode " +
          std::string(1, mode) + ", but should not have.\n" + mistaken_passes);
    break;
  case 'D':
  case 'E':
    break;
  }
}

//-------------------------------------------------------------------------------------------------
// Test the individual gears of the LocalExclusionMask.  Many of the functions supporting the
// LocalExclusionMask are free functions which this routine can drive.
//
// Arguments:
//   xrs:  Random number generator to differentiate sequences and permit repeated tests
//-------------------------------------------------------------------------------------------------
void testLocalExclusionMaskMechanics(Xoshiro256ppGenerator *xrs) {

  // Test some exclusion lists which should fit mode A.
  const std::vector<int2> a_limits(1, { -lmask_long_local_span, lmask_long_local_span });
  testExclusionZones(xrs, a_limits, 'A');

  // Test some exclusion lists which should fit mode B.
  std::vector<int2> b_limits(3);
  const int b_max_shift = ipow(2, lmask_b_shft_bits) - 1;
  b_limits[0] = { -lmask_short_local_span, lmask_short_local_span };
  b_limits[1] = { -lmask_short_local_span - b_max_shift - lmask_short_extra_span,
                  -lmask_short_local_span - b_max_shift - 1 };
  b_limits[2] = { lmask_short_local_span + b_max_shift,
                  lmask_short_local_span + b_max_shift + lmask_short_extra_span -1 };
  testExclusionZones(xrs, b_limits, 'B');
  const int min_diff = lmask_long_local_span - lmask_short_local_span + 1;
  b_limits[1] = { -lmask_short_local_span - min_diff - lmask_short_extra_span,
                  -lmask_short_local_span - min_diff - 1};
  b_limits[2] = { lmask_short_local_span + min_diff,
                  lmask_short_local_span + min_diff + lmask_short_extra_span - 1 };
  testExclusionZones(xrs, b_limits, 'B');
  const int spacer_low = std::max(lmask_short_extra_span - 6, 0);
  const int spacer_hgh = std::max(lmask_short_extra_span - 3, 0);
  b_limits[1] = { -lmask_short_local_span - min_diff - spacer_low - lmask_short_extra_span,
                  -lmask_short_local_span - min_diff - spacer_low - 1};
  b_limits[2] = { lmask_short_local_span + min_diff + spacer_hgh,
                  lmask_short_local_span + min_diff + spacer_hgh + lmask_short_extra_span - 1 };
  testExclusionZones(xrs, b_limits, 'B');

  // Test mode C.
  std::vector<int2> c_limits(2);
  c_limits[0] = { -lmask_short_local_span, lmask_short_local_span };
  const int lowest_shift  = -ipow(2, lmask_c_shft_bits - 1);
  c_limits[1] = { lowest_shift, lowest_shift + lmask_long_extra_span - 1 };
  testExclusionZones(xrs, c_limits, 'C');
  const int highest_shift =  ipow(2, lmask_c_shft_bits - 1) - 1;
  c_limits[1] = { highest_shift, highest_shift + lmask_long_extra_span - 1 };
  testExclusionZones(xrs, c_limits, 'C');
  const int middle_shift_plus = highest_shift / 3;
  c_limits[1] = { middle_shift_plus, middle_shift_plus + lmask_long_extra_span - 1 };
  testExclusionZones(xrs, c_limits, 'C');
  const int middle_shift_minus = lowest_shift / 5;
  c_limits[1] = { middle_shift_minus, middle_shift_minus + lmask_long_extra_span - 1 };
  testExclusionZones(xrs, c_limits, 'C');

  // Test mode D.  The D mode profile can accommodate local exclusions plus up to
  // 2^(lmask_e_cnt_bits) - 1 additional masks.
  const int n_dsegments = ipow(2, lmask_d_array_cnt_bits);
  std::vector<int2> d_limits(n_dsegments);
  d_limits[0] = { -lmask_short_local_span, lmask_short_local_span };
  for (int i = 1; i < n_dsegments; i++) {
    bool collision;
    int new_start;
    do {

      // Select a stretch of atoms that is nowhere colliding with a preceding stretch of atoms.
      new_start = 65535 * xrs->uniformRandomNumber();
      collision = false;
      for (int j = 0; j < i; j++) {
        collision = (collision || (new_start >= d_limits[j].x && new_start <= d_limits[j].y));
      }
    } while (collision);
    d_limits[i] = { new_start, new_start + 31 };
  }
  testExclusionZones(xrs, d_limits, 'D');
}

//-------------------------------------------------------------------------------------------------
// Compute the exclusions map of all atoms to all others in a topology.
//
// Arguments:
//   ag:  The topology of interest
//-------------------------------------------------------------------------------------------------
std::vector<bool> makeFullExclusionMap(const AtomGraph *ag) {
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  std::vector<bool> result;
  const int max_atoms = 8192;
  if (nbk.natom < max_atoms) {
    result.resize(nbk.natom * nbk.natom, false);
    for (int i = 0; i < nbk.natom; i++) {
      result[(i * nbk.natom) + i] = true;
      for (int j = nbk.nb11_bounds[i]; j < nbk.nb11_bounds[i + 1]; j++) {
        const int k = nbk.nb11x[j];
        result[(i * nbk.natom) + k] = true;
      }
      for (int j = nbk.nb12_bounds[i]; j < nbk.nb12_bounds[i + 1]; j++) {
        const int k = nbk.nb12x[j];
        result[(i * nbk.natom) + k] = true;
      }
      for (int j = nbk.nb13_bounds[i]; j < nbk.nb13_bounds[i + 1]; j++) {
        const int k = nbk.nb13x[j];
        result[(i * nbk.natom) + k] = true;
      }
      for (int j = nbk.nb14_bounds[i]; j < nbk.nb14_bounds[i + 1]; j++) {
        const int k = nbk.nb14x[j];
        result[(i * nbk.natom) + k] = true;
      }
    }
  }
  else {
    result.resize(max_atoms * max_atoms, false);
    for (int i = 0; i < max_atoms; i++) {
      result[(i * max_atoms) + i] = true;
      for (int j = nbk.nb11_bounds[i]; j < nbk.nb11_bounds[i + 1]; j++) {
        const int k = nbk.nb11x[j];
        if (k < max_atoms) {
          result[(i * max_atoms) + k] = true;
        }
      }
      for (int j = nbk.nb12_bounds[i]; j < nbk.nb12_bounds[i + 1]; j++) {
        const int k = nbk.nb12x[j];
        if (k < max_atoms) {
          result[(i * max_atoms) + k] = true;
        }
      }
      for (int j = nbk.nb13_bounds[i]; j < nbk.nb13_bounds[i + 1]; j++) {
        const int k = nbk.nb13x[j];
        if (k < max_atoms) {
          result[(i * max_atoms) + k] = true;
        }
      }
      for (int j = nbk.nb14_bounds[i]; j < nbk.nb14_bounds[i + 1]; j++) {
        const int k = nbk.nb14x[j];
        if (k < max_atoms) {
          result[(i * max_atoms) + k] = true;
        }
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Compile a list of missing or erroneous exclusions into a human-readable string.
//
// Arguments:
//   bad_exclusions:  List of bad exclusions, with the ith atom in the "x" member, the jth atom in
//                    the "y" member, and the profile mode of the exclusion in the "z" member
//   descriptor:      Single word to descirbe the exclusions (will be capitalized as needed)
//   cdk:             Contains atom names for more precise reporting of errors
//-------------------------------------------------------------------------------------------------
std::string compileExclusionReport(const std::vector<int3> &bad_exclusions,
                                   const std::string &descriptor, const ChemicalDetailsKit &cdk) {
  const int n_bad = std::min(bad_exclusions.size(), static_cast<size_t>(8));
  std::string result;
  std::string actual_desc = lowercase(descriptor);
  if (n_bad <= 8) {
    if (actual_desc.size() > 0) {
      actual_desc[0] = uppercase(actual_desc[0]);
    }
    result += actual_desc + " exclusions include: ";      
  }
  else {
    result += "Examples of " + actual_desc + " exclusions include: ";      
  }
  for (int i = 0; i < n_bad; i++) {
    result += "[ " + std::to_string(bad_exclusions[i].x) + " (" +
              char4ToString(cdk.atom_names[bad_exclusions[i].x]) + ") -> " +
              std::to_string(bad_exclusions[i].y) + " (" +
              char4ToString(cdk.atom_names[bad_exclusions[i].y]) + "), Mode " +
              std::string(1, static_cast<char>('A' + bad_exclusions[i].z)) + " ]";
    if (i < n_bad - 1) {
      result += " ";
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Test the LocalExclusionMask on typical molecular systems.
//
// Arguments:
//   ag:        The topology of interest
//   timer:     Tracking object for execution wall time
//   do_tests:  Indicate whether tests are possible
//-------------------------------------------------------------------------------------------------
void testLocalExclusionMask(const AtomGraph *ag, StopWatch *timer, const TestPriority do_tests) {
  const int build_time = timer->addCategory("Local Excl. Mask (AtomGraph)");
  timer->assignTime(0);
  const LocalExclusionMask lemask(ag);
  timer->assignTime(build_time);

  // Loop over all exclusions in the topology
  const std::vector<bool> chk_excl = makeFullExclusionMap(ag);
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  std::vector<int3> missed_pairs;
  std::vector<int3> wrong_pairs;
  const int testable_atoms = std::min(nbk.natom, 8192);
  for (int i = 0; i < testable_atoms; i++) {
    for (int j = 0; j < testable_atoms; j++) {
      const bool lmask_excludes = lemask.testExclusion(i, j);
      if (lmask_excludes && chk_excl[(i * testable_atoms) + j] == false) {
        wrong_pairs.push_back({ i, j, lemask.getMode(i) });
      }
      if (lmask_excludes == false && chk_excl[(i * testable_atoms) + j]) {
        missed_pairs.push_back({ i, j, lemask.getMode(i) });
      }
    }
  }
  std::string missed_pair_report = compileExclusionReport(missed_pairs, "missing", cdk);
  std::string wrong_pair_report  = compileExclusionReport(wrong_pairs, "wrong", cdk);
  check(missed_pairs.size() == 0, "A total of " + std::to_string(missed_pairs.size()) +
        " exclusions were missed by the LocalExclusionMask applied to topology " +
        getBaseName(ag->getFileName()) + ".  " + missed_pair_report);
  check(wrong_pairs.size() == 0, "A total of " + std::to_string(wrong_pairs.size()) +
        " exclusions were reported by the LocalExclusionMask applied to topology " +
        getBaseName(ag->getFileName()) + ", but have no basis in the actual system.  " +
        wrong_pair_report);
}

//-------------------------------------------------------------------------------------------------
// Test the LocalExclusionMask on syntheses of typical molecular systems.
//
// Arguments:
//   poly_ag:   The topology synthesis to analyze
//   timer:     Tracking object for execution wall time
//   do_tests:  Indicate whether tests are possible
//-------------------------------------------------------------------------------------------------
void testLocalExclusionMask(const AtomGraphSynthesis &poly_ag, StopWatch *timer,
                            const TestPriority do_tests) {
  const int build_time = timer->addCategory("Local Excl. Mask (Synthesis)");
  timer->assignTime(0);
  const LocalExclusionMask lemask(poly_ag);
  timer->assignTime(build_time);
  
  // Loop over all systems and perform the exclusion check
  for (int pos = 0; pos < poly_ag.getSystemCount(); pos++) {
    const AtomGraph *ag_ptr = poly_ag.getSystemTopologyPointer(pos);
    const std::vector<bool> chk_excl = makeFullExclusionMap(ag_ptr);
    const NonbondedKit<double> nbk = ag_ptr->getDoublePrecisionNonbondedKit();
    const ChemicalDetailsKit cdk = ag_ptr->getChemicalDetailsKit();
    std::vector<int3> missed_pairs;
    std::vector<int3> wrong_pairs;
    const int pos_offset = poly_ag.getAtomOffset(pos);
    for (int i = 0; i < nbk.natom; i++) {
      const int pos_i = pos_offset + i;
      for (int j = 0; j < nbk.natom; j++) {
        const int pos_j = pos_offset + j;
        const bool lmask_excludes = lemask.testExclusion(pos_i, pos_j);
        if (lmask_excludes && chk_excl[(i * nbk.natom) + j] == false) {
          wrong_pairs.push_back({ i, j, lemask.getMode(pos_i) });
        }
        if (lmask_excludes == false && chk_excl[(i * nbk.natom) + j]) {
          missed_pairs.push_back({ i, j, lemask.getMode(pos_i) });
        }
      }
    }
    std::string missed_pair_report = compileExclusionReport(missed_pairs, "missing", cdk);
    std::string wrong_pair_report  = compileExclusionReport(wrong_pairs, "wrong", cdk);
    check(missed_pairs.size() == 0, "A total of " + std::to_string(missed_pairs.size()) +
          " exclusions were missed by the LocalExclusionMask applied to topology " +
          getBaseName(ag_ptr->getFileName()) + ", when it was compiled into a topology "
          "synthesis.  " + missed_pair_report);
    check(wrong_pairs.size() == 0, "A total of " + std::to_string(wrong_pairs.size()) +
          " exclusions were reported by the LocalExclusionMask applied to topology " +
          getBaseName(ag_ptr->getFileName()) + ", but have no basis in the actual system.  This "
          "occurred when the topology was compiled into a synthesis.  " + wrong_pair_report);
  }
}

//-------------------------------------------------------------------------------------------------
// Check the fidelity of a compilation of logarithmic spline tables, in both precision modes and
// over all energy and force functions.
//
// Arguments:
//   nrg_tab:           The compilation of spline tables to test
//   xrs:               Source of pseudo-random numbers to generate displacements for energy and
//                      force calculations
//   min_displacement:  The minimum displacement between particles
//-------------------------------------------------------------------------------------------------
void checkParticleInteractionTable(const PPITable &nrg_tab,
                                   Xoshiro256ppGenerator *xrs,
                                   const double min_displacement = 2.0) {
  const int npts = 100;
  std::vector<double> nrg_rf(npts), frc_rf(npts), nrg_chk(npts), frc_chk(npts);
  std::vector<double> nrg_rd(npts), frc_rd(npts);
  std::vector<double> nrgx_rf(npts), frcx_rf(npts), nrgx_chk(npts), frcx_chk(npts);
  std::vector<double> nrgx_rd(npts), frcx_rd(npts);
  const double cutoff = nrg_tab.getCutoff();
  const double ew_coeff = nrg_tab.getEwaldCoefficient();
  for (int i = 0; i < npts; i++) {
    const double r = (xrs->uniformRandomNumber() * (cutoff - min_displacement)) + min_displacement;
    nrg_rf[i] = nrg_tab.evaluateByRealArg(r, LogSplineForm::ELEC_PME_DIRECT);
    frc_rf[i] = nrg_tab.evaluateByRealArg(r, LogSplineForm::DELEC_PME_DIRECT);
    nrg_chk[i] = elecPMEDirectSpace(ew_coeff, amber_ancient_bioq, r, r * r, 0);
    frc_chk[i] = elecPMEDirectSpace(ew_coeff, amber_ancient_bioq, r, r * r, 1) / r;
    nrg_rd[i] = nrg_tab.evaluateByRealArg(r, LogSplineForm::ELEC_PME_DIRECT,
                                          PrecisionModel::DOUBLE);
    frc_rd[i] = nrg_tab.evaluateByRealArg(r, LogSplineForm::DELEC_PME_DIRECT,
                                          PrecisionModel::DOUBLE);
    nrgx_rf[i] = nrg_tab.evaluateByRealArg(r, LogSplineForm::ELEC_PME_DIRECT_EXCL);
    frcx_rf[i] = nrg_tab.evaluateByRealArg(r, LogSplineForm::DELEC_PME_DIRECT_EXCL);
    nrgx_chk[i] = elecPMEDirectSpace(ew_coeff, amber_ancient_bioq, r, r * r, 0) -
                  (amber_ancient_bioq / r);
    frcx_chk[i] = (elecPMEDirectSpace(ew_coeff, amber_ancient_bioq, r, r * r, 1) / r) +
                  (amber_ancient_bioq / (r * r * r));
    nrgx_rd[i] = nrg_tab.evaluateByRealArg(r, LogSplineForm::ELEC_PME_DIRECT_EXCL,
                                           PrecisionModel::DOUBLE);
    frcx_rd[i] = nrg_tab.evaluateByRealArg(r, LogSplineForm::DELEC_PME_DIRECT_EXCL,
                                           PrecisionModel::DOUBLE);
  }
  check(nrg_rf, RelationalOperator::EQUAL, Approx(nrg_chk).margin(7.5e-5), "Energies computed "
        "with splines in " + getEnumerationName(PrecisionModel::SINGLE) + " precision do not "
        "meet expectations.");
  check(nrg_rd, RelationalOperator::EQUAL, Approx(nrg_chk).margin(7.5e-5), "Energies computed "
        "with splines in " + getEnumerationName(PrecisionModel::DOUBLE) + " precision do not "
        "meet expectations.");
  check(frc_rf, RelationalOperator::EQUAL, Approx(frc_chk).margin(1.5e-4), "Forces computed "
        "with splines in " + getEnumerationName(PrecisionModel::SINGLE) + " precision do not "
        "meet expectations.");
  check(frc_rd, RelationalOperator::EQUAL, Approx(frc_chk).margin(1.5e-4), "Forces computed "
        "with splines in " + getEnumerationName(PrecisionModel::DOUBLE) + " precision do not "
        "meet expectations.");
  check(nrgx_rf, RelationalOperator::EQUAL, Approx(nrgx_chk).margin(1.5e-5), "Excluded energies "
        "computed with splines in " + getEnumerationName(PrecisionModel::SINGLE) + " precision "
        "do not meet expectations.");
  check(nrgx_rd, RelationalOperator::EQUAL, Approx(nrgx_chk).margin(1.5e-5), "Excluded energies "
        "computed with splines in " + getEnumerationName(PrecisionModel::DOUBLE) + " precision "
        "do not meet expectations.");
  check(frcx_rf, RelationalOperator::EQUAL, Approx(frcx_chk).margin(1.5e-5), "Excluded forces "
        "computed with splines in " + getEnumerationName(PrecisionModel::SINGLE) + " precision "
        "do not meet expectations.");
  check(frcx_rd, RelationalOperator::EQUAL, Approx(frcx_chk).margin(1.5e-5), "Excluded forces "
        "computed with splines in " + getEnumerationName(PrecisionModel::DOUBLE) + " precision "
        "do not meet expectations.");
}

//-------------------------------------------------------------------------------------------------
// Test the particle-particle interaction potentials as they are collected for use in simulations.
//
// Arguments:
//   cutoff:            The extent of particle-particle interactions
//   mantissa_bits:     The number of bits out of the floating point argument's mantissa to use
//                      when indexing into the tables
//   bss:               The basis functions to be used in the spline approximation
//   xrs:               Source of pseudo-random numbers to generate displacements for energy and
//                      force calculations
//   min_displacement:  The minimum displacement between particles
//-------------------------------------------------------------------------------------------------
template <typename T4>
void testPPSplineTables(const double cutoff, const int mantissa_bits, const BasisFunctions bss,
                        Xoshiro256ppGenerator *xrs, const double min_displacement = 2.0) {

  const std::vector<LogSplineForm> components = { LogSplineForm::ELEC_PME_DIRECT,
                                                  LogSplineForm::DELEC_PME_DIRECT,
                                                  LogSplineForm::ELEC_PME_DIRECT_EXCL,
                                                  LogSplineForm::DELEC_PME_DIRECT_EXCL };
  const double padded_cutoff = exp2(ceil(log2(cutoff)));
  if (std::type_index(typeid(T4)).hash_code() == double4_type_index && mantissa_bits == 4) {
    
    // Test the various constructions for the table with otherwise cheap parameters.
    for (int i = 0; i < 4; i++) {
      const LogScaleSpline<T4> lfunc(components[i], ewaldCoefficient(cutoff, default_dsum_tol),
                                     amber_ancient_bioq, mantissa_bits, padded_cutoff, 0.015625,
                                     TableIndexing::ARG, bss, 0, 0.0);
      const PPITable nrg_tab(lfunc);
      checkParticleInteractionTable(nrg_tab, xrs, min_displacement);
    }
    for (int i = 0; i < 3; i++) {
      const LogScaleSpline<T4> lfunc(components[i], ewaldCoefficient(cutoff, default_dsum_tol),
                                     amber_ancient_bioq, mantissa_bits, padded_cutoff, 0.015625,
                                     TableIndexing::ARG, bss, 0, 0.0);
      for (int j = i + 1; j < 4; j++) {
        const LogScaleSpline<T4> kfunc(components[j], ewaldCoefficient(cutoff, default_dsum_tol),
                                       amber_ancient_bioq, mantissa_bits, padded_cutoff, 0.015625,
                                       TableIndexing::ARG, bss, 0, 0.0);
        const PPITable nrg_tab(lfunc, kfunc);
        checkParticleInteractionTable(nrg_tab, xrs, min_displacement);
      }
    }
    for (int i = 0; i < 2; i++) {
      const LogScaleSpline<T4> lfunc(components[i], ewaldCoefficient(cutoff, default_dsum_tol),
                                     amber_ancient_bioq, mantissa_bits, padded_cutoff, 0.015625,
                                     TableIndexing::ARG, bss, 0, 0.0);
      for (int j = i + 1; j < 3; j++) {
        const LogScaleSpline<T4> jfunc(components[j], ewaldCoefficient(cutoff, default_dsum_tol),
                                       amber_ancient_bioq, mantissa_bits, padded_cutoff, 0.015625,
                                       TableIndexing::ARG, bss, 0, 0.0);
        for (int k = j + 1; k < 4; k++) {
          const LogScaleSpline<T4> kfunc(components[k], ewaldCoefficient(cutoff, default_dsum_tol),
                                         amber_ancient_bioq, mantissa_bits, padded_cutoff,
                                         0.015625, TableIndexing::ARG, bss, 0, 0.0);
          const PPITable nrg_tab(lfunc, jfunc, kfunc);
          checkParticleInteractionTable(nrg_tab, xrs, min_displacement);
        }
      }
    }

    // Try constructing a table without templated parameters.  This will only be executed once for
    // a particular templated form of this function.
    const PPITable fixed_nrg_tab(NonbondedTheme::ELECTROSTATIC, BasisFunctions::POLYNOMIAL,
                                 TableIndexing::SQUARED_ARG, 10.5, 0.0, 1.0e-6, 5,
                                 amber_ancient_bioq, 0.015625);
    checkParticleInteractionTable(fixed_nrg_tab, xrs, min_displacement);
  }

  // Test the energy and force values coming out of the table
  const double ew_coeff = ewaldCoefficient(cutoff, default_dsum_tol);
  const LogScaleSpline<T4> mfunc(components[0], ew_coeff, amber_ancient_bioq, mantissa_bits,
                                 padded_cutoff, 0.015625, TableIndexing::SQUARED_ARG, bss, 0, 0.0);
  const PPITable nrg_tab_ii(mfunc);
  checkParticleInteractionTable(nrg_tab_ii, xrs, min_displacement);
}

//-------------------------------------------------------------------------------------------------
// Test the particle migration features within a CellGrid object.
//
// Arguments:
//   cg:        The cell grid neighbor list ready for testing
//   poly_ps:   The associated coordinate synthesis (provided as a formal argument to avoid having
//              to const-cast the pointer if it were harvested from cg itself)
//   xrs:       Source of random numbers
//   do_tests:  Indicate whether testing is possible
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void testCpuCellMigration(CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> *cg, PhaseSpaceSynthesis *poly_ps,
                          Xoshiro256ppGenerator *xrs, const TestPriority do_tests) {

  // Perturb particle positions in the associated coordinate synthesis.  Place the perturbed
  // positions in the next move cycle.
  PsSynthesisWriter poly_psw = poly_ps->data();
  std::vector<std::vector<double>> x_moves(poly_psw.system_count);
  std::vector<std::vector<double>> y_moves(poly_psw.system_count);
  std::vector<std::vector<double>> z_moves(poly_psw.system_count);  
  for (int i = 0; i < poly_psw.system_count; i++) {
    const size_t illim = poly_psw.atom_starts[i];
    const size_t ihlim = illim + poly_psw.atom_counts[i];
    x_moves[i].resize(ihlim - illim);
    y_moves[i].resize(ihlim - illim);
    z_moves[i].resize(ihlim - illim);
    for (size_t j = illim; j < ihlim; j++) {
      const double xpos = hostInt95ToDouble(poly_psw.xcrd[j], poly_psw.xcrd_ovrf[j]);
      const double ypos = hostInt95ToDouble(poly_psw.ycrd[j], poly_psw.ycrd_ovrf[j]);
      const double zpos = hostInt95ToDouble(poly_psw.zcrd[j], poly_psw.zcrd_ovrf[j]);
      x_moves[i][j - illim] = xrs->gaussianRandomNumber();
      y_moves[i][j - illim] = xrs->gaussianRandomNumber();
      z_moves[i][j - illim] = xrs->gaussianRandomNumber();
      const double nxpos = xpos + (x_moves[i][j - illim] * poly_psw.gpos_scale);
      const double nypos = ypos + (y_moves[i][j - illim] * poly_psw.gpos_scale);
      const double nzpos = zpos + (z_moves[i][j - illim] * poly_psw.gpos_scale);
      const int95_t inxpos = hostDoubleToInt95(nxpos);
      const int95_t inypos = hostDoubleToInt95(nypos);
      const int95_t inzpos = hostDoubleToInt95(nzpos);
      poly_psw.xalt[j] = inxpos.x;
      poly_psw.yalt[j] = inypos.x;
      poly_psw.zalt[j] = inzpos.x;
      poly_psw.xalt_ovrf[j] = inxpos.y;
      poly_psw.yalt_ovrf[j] = inypos.y;
      poly_psw.zalt_ovrf[j] = inzpos.y;
    }
  }
  cg->updatePositions();
  cg->updateCyclePosition();
  poly_ps->updateCyclePosition();
  CellGridWriter<Tcoord, Tacc, Tcalc, Tcoord4> cgw = cg->data();
  
  std::vector<int> base_misalignments(poly_psw.system_count, 0);
  bool particles_unplaceable = false;
  const int xfrm_stride = roundUp(9, warp_size_int);
  for (int i = 0; i < poly_psw.system_count; i++) {
    const PhaseSpace ps = poly_ps->exportSystem(i);
    const PhaseSpaceReader psr = ps.data();

    // Check that the structural evolution is reflected in the PhaseSpace object
    for (int j = 0; j < psr.natom; j++) {
      base_misalignments[i] += (fabs(psr.xcrd[j] - x_moves[i][j] - psr.xalt[j]) > 1.0e-6 ||
                                fabs(psr.ycrd[j] - y_moves[i][j] - psr.yalt[j]) > 1.0e-6 ||
                                fabs(psr.zcrd[j] - z_moves[i][j] - psr.zalt[j]) > 1.0e-6);
    }
    std::vector<double> cg_xpos(psr.natom), cg_ypos(psr.natom), cg_zpos(psr.natom);
    std::vector<double> ps_xpos(psr.natom), ps_ypos(psr.natom), ps_zpos(psr.natom);

    // Check that the cell grid accounts for each new particle position correctly.
    const ullint grid_bounds = cgw.system_cell_grids[i];
    const int cell_start = (grid_bounds & 0xfffffff);
    const int ncell_a = ((grid_bounds >> 28) & 0xfff);
    const int ncell_b = ((grid_bounds >> 40) & 0xfff);
    const int ncell_c = (grid_bounds >> 52);
    const uint base_img_idx = cgw.cell_limits[cell_start].x;
    const Tcoord* invu_p = &cgw.system_cell_invu[xfrm_stride * i];
    for (int j = 0; j < psr.natom; j++) {
      const int synth_idx = poly_psw.atom_starts[i] + j;
      const uint img_idx = cgw.img_atom_idx[synth_idx];
      const uint img_del_idx = img_idx - base_img_idx;
      const int chain_idx = img_del_idx / static_cast<uint>(cgw.cell_base_capacity * ncell_a);
      const int cell_cidx = chain_idx / ncell_b;
      const int cell_bidx = chain_idx - (cell_cidx * ncell_b);
      int cell_aidx = -1;
      const int base_chain_cell = cell_start + (((cell_cidx * ncell_b) + cell_bidx) * ncell_a);
      for (int k = 0; k < ncell_a; k++) {
        const uint2 kcell_lims = cgw.cell_limits[base_chain_cell + k];
        const uint kcell_max = kcell_lims.x + (kcell_lims.y >> 16);
        if (img_idx >= kcell_lims.x && img_idx < kcell_max) {
          cell_aidx = k;
        }
      }
      particles_unplaceable = ((cell_aidx == -1) || particles_unplaceable);
      const Tcoord4 xyz_prop = cgw.image[img_idx];
      cg_xpos[j] = xyz_prop.x * cgw.lpos_scale;
      cg_ypos[j] = xyz_prop.y * cgw.lpos_scale;
      cg_zpos[j] = xyz_prop.z * cgw.lpos_scale;
      const double dcell_aidx = cell_aidx;
      const double dcell_bidx = cell_bidx;
      const double dcell_cidx = cell_cidx;
      cg_xpos[j] += (invu_p[0] * dcell_aidx) + (invu_p[3] * dcell_bidx) + (invu_p[6] * dcell_cidx);
      cg_ypos[j] +=                            (invu_p[4] * dcell_bidx) + (invu_p[7] * dcell_cidx);
      cg_zpos[j] +=                                                       (invu_p[8] * dcell_cidx);
      ps_xpos[j] = psr.xcrd[j];
      ps_ypos[j] = psr.ycrd[j];
      ps_zpos[j] = psr.zcrd[j];
    }
    imageCoordinates<double, double>(&ps_xpos, &ps_ypos, &ps_zpos, psr.umat, psr.invu,
                                     psr.unit_cell, ImagingMethod::PRIMARY_UNIT_CELL);
    const AtomGraph *ag = poly_ps->getSystemTopologyPointer(i);
    const std::string err_append = " (system " + getBaseName(ag->getFileName()) + ").";
    check(cg_xpos, RelationalOperator::EQUAL, ps_xpos, "CellGrid X coordinates do not match those "
          "of the underlying coordinate synthesis after particle migration.", do_tests);
    check(cg_ypos, RelationalOperator::EQUAL, ps_ypos, "CellGrid Y coordinates do not match those "
          "of the underlying coordinate synthesis after particle migration.", do_tests);
    check(cg_zpos, RelationalOperator::EQUAL, ps_zpos, "CellGrid Z coordinates do not match those "
          "of the underlying coordinate synthesis after particle migration.", do_tests);
  }
  check(base_misalignments, RelationalOperator::EQUAL, std::vector<int>(poly_psw.system_count, 0),
        "Some particles were not properly updated in the coordinate synthesis, as judged by "
        "comparing the two time cycles of exported PhaseSpace objects.", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Check the mechanics of cell grid placement using the pre-determined neighbor list cell origins.
//
// Arguments:
//   tsm:         A collection of test systems
//   system_idx:  The index of the system within the collection to test
//   fp_bits:     The number of fixed-precision bits after the point in the coordinate synthesis
//                that will be constructed from the testing apparatus (ions)
//   xrs:         Source of random numbers for selecting a subset of neighbor list cell origins
//                to test
//-------------------------------------------------------------------------------------------------
void testCellOrigins(const TestSystemManager &tsm, const std::vector<int> &system_idx,
                     const int fp_bits, Xoshiro256ppGenerator *xrs) {
  PhaseSpaceSynthesis poly_ps = tsm.exportPhaseSpaceSynthesis(system_idx, 200.0, 9286943, fp_bits);
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(system_idx);
  PhaseSpaceSynthesis poly_ps_copy = poly_ps;
  PsSynthesisWriter poly_psw = poly_ps.data();
  CellGrid<double, llint, double, double4> cg(&poly_ps, poly_ag, 4.8, 0.05, 4,
                                              NonbondedTheme::ALL);
  CellGridWriter cgw = cg.data();
  const CellOriginsReader corg = cg.getRulers();
  const int xfrm_stride = roundUp(9, warp_size_int);
  std::vector<double3> xfrm_orig, tabl_orig;
  std::vector<int4> sys_membership;
  std::vector<UnitCellType> sys_unit_cell;
  for (int i = 0; i < poly_psw.system_count; i++) {

    // Check the accuracy of cell grid origins
    const ullint sys_grids = cgw.system_cell_grids[i];
    const int cell_na = ((sys_grids >> 28) & 0xfff);
    const int cell_nb = ((sys_grids >> 40) & 0xfff);
    const int cell_nc = (sys_grids >> 52);
    const double cell_dna = cell_na;
    const double cell_dnb = cell_nb;
    const double cell_dnc = cell_nc;
    const double* umat = &poly_psw.umat[(i * xfrm_stride)];
    const double* invu = &poly_psw.invu[(i * xfrm_stride)];
    for (int j = 0; j < cell_na; j++) {
      const double dj = j;
      for (int k = 0; k < cell_nb; k++) {
        const double dk = k;
        for (int m = 0; m < cell_nc; m++) {

          // Test only a portion of all neighbor list cell origins.
          if (xrs->uniformRandomNumber() < 0.92) {
            continue;
          }
          const double dm = m;

          // The fractional coordinates of the neighbor list cell origin are alreayd at hand.
          // Compute the Cartesian position.
          const double da = dj / cell_dna;
          const double db = dk / cell_dnb;
          const double dc = dm / cell_dnc;
          const double dx = (invu[0] * da) + (invu[3] * db) + (invu[6] * dc);
          const double dy =                  (invu[4] * db) + (invu[7] * dc);
          const double dz =                                   (invu[8] * dc);

          // Compute the Cartesian position of the neighbor list cell origin based on the cell
          // grid's rulers.
          const int ruler_j = (i * corg.stride) + j;
          const int ruler_k = (i * corg.stride) + k;
          const int ruler_m = (i * corg.stride) + m;
          const int95_t pt_ax = { corg.ax[ruler_j], corg.ax_ovrf[ruler_j] };
          const int95_t pt_bx = { corg.bx[ruler_k], corg.bx_ovrf[ruler_k] };
          const int95_t pt_by = { corg.by[ruler_k], corg.by_ovrf[ruler_k] };
          const int95_t pt_cx = { corg.cx[ruler_m], corg.cx_ovrf[ruler_m] };
          const int95_t pt_cy = { corg.cy[ruler_m], corg.cy_ovrf[ruler_m] };
          const int95_t pt_cz = { corg.cz[ruler_m], corg.cz_ovrf[ruler_m] };
          const double pt_x = hostInt95ToDouble(hostSplitFPSum(hostSplitFPSum(pt_ax, pt_bx),
                                                               pt_cx)) * poly_psw.inv_gpos_scale_f;
          const double pt_y = hostInt95ToDouble(hostSplitFPSum(pt_by, pt_cy)) *
                              poly_psw.inv_gpos_scale_f;
          const double pt_z = hostInt95ToDouble(pt_cz) * poly_psw.inv_gpos_scale_f;
          if (fabs(dx - pt_x) > 1.0e-5 || fabs(dy - pt_y) > 1.0e-5 || fabs(dz - pt_z) > 1.0e-5) {
            xfrm_orig.push_back({ dx, dy, dz });
            tabl_orig.push_back({ pt_x, pt_y, pt_z });
            sys_membership.push_back({ j, k, m, i });
            sys_unit_cell.push_back(poly_ps.getSystemTopologyPointer(i)->getUnitCellType());
          }
        }
      }
    }
  }
  std::string err_msg("Examples of errors: ");
  if (sys_membership.size() > 0) {
    const int n_report = std::min(8, static_cast<int>(sys_membership.size()));
    for (int i = 0; i < n_report; i++) {
      err_msg += "Cell [ " + std::to_string(sys_membership[i].x) + " " +
                 std::to_string(sys_membership[i].y) + " " + std::to_string(sys_membership[i].z) +
                 " ] Transform gets [ " +
                 realToString(xfrm_orig[i].x, 9, 4, NumberFormat::STANDARD_REAL) + " " +
                 realToString(xfrm_orig[i].x, 9, 4, NumberFormat::STANDARD_REAL) + " " +
                 realToString(xfrm_orig[i].x, 9, 4, NumberFormat::STANDARD_REAL) +
                 " ], rulers get [ " +
                 realToString(tabl_orig[i].x, 9, 4, NumberFormat::STANDARD_REAL) + " " +
                 realToString(tabl_orig[i].x, 9, 4, NumberFormat::STANDARD_REAL) + " " +
                 realToString(tabl_orig[i].x, 9, 4, NumberFormat::STANDARD_REAL) +
                 " ], system type " + getEnumerationName(sys_unit_cell[i]) + ".";
      if (i < n_report - 1) {
        err_msg += "  ";
      }
    }
  }
  check(sys_membership.size(), RelationalOperator::EQUAL, 0, "Some of the tested neighbor list "
        "cells were not placed at the same locations when calculated with the CellGrid's rulers "
        "as opposed to direct matrix transformations.  " + err_msg, tsm.getTestingStatus());
  
}

//-------------------------------------------------------------------------------------------------
// Run tests of the GPU migration kernels.  These complement tests for the CPU updates.
//
// Arguments:
//   ions:     Group of (ionic) systems for testing
//   ncopy:    The number of copies of the system in the ions test system manager to collect into
//             a synthesis
//   crd_tol:  Tolerance for coordinate placements in the new neighbor list
//   stage:    Indicate, by number, one of several tests to perform.  1: add random noise in a way
//             that will generate a large number of migrators and perhaps a handful of "wanderers".
//             2: arrange ions around a particular chain and then move a large number of them into
//             that chain, without making them wanderers but in such a way as to overload the
//             notes buffer.
//   xrs:      Source of random numbers for reposition particles
//   gpu:      Specifications of the GPU to use in calculations
//-------------------------------------------------------------------------------------------------
#ifdef STORMM_USE_HPC
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void testGpuCellMigration(const TestSystemManager &ions, const int ncopy, const double crd_tol,
                          const int stage, Xoshiro256ppGenerator *xrs, const GpuDetails &gpu) {
  const size_t ct_coord = std::type_index(typeid(Tcoord)).hash_code();
  const bool tcoord_is_long = (ct_coord == double_type_index || ct_coord == llint_type_index);
  const bool tcoord_is_integral = (ct_coord == int_type_index || ct_coord == llint_type_index);
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  PhaseSpaceSynthesis poly_ps = ions.exportPhaseSpaceSynthesis(std::vector<int>(ncopy, 0), 0.0,
                                                               71583329, 35);
  PsSynthesisWriter poly_psw = poly_ps.data();
  const int xfrm_stride = roundUp(9, warp_size_int);
  switch (stage) {
  case 1:
    break;
  case 2:
    for (int i = 0; i < poly_psw.system_count; i++) {
      const int llim = poly_psw.atom_starts[i];
      const int hlim = llim + poly_psw.atom_counts[i];
      const double xspan = poly_psw.invu[(i * xfrm_stride)    ];
      const double yspan = poly_psw.invu[(i * xfrm_stride) + 4];
      const double zspan = poly_psw.invu[(i * xfrm_stride) + 8];
      for (int j = llim; j < hlim; j++) {
        const double orig_x = xspan * xrs->uniformRandomNumber();
        const double orig_y = yspan * (0.17 + (0.49 * xrs->uniformRandomNumber()));
        const double orig_z = zspan * (0.17 + (0.49 * xrs->uniformRandomNumber()));
        const int95_t orig_xi = hostDoubleToInt95(orig_x * poly_psw.gpos_scale);
        const int95_t orig_yi = hostDoubleToInt95(orig_y * poly_psw.gpos_scale);
        const int95_t orig_zi = hostDoubleToInt95(orig_z * poly_psw.gpos_scale);
        poly_psw.xcrd[j] = orig_xi.x;
        poly_psw.ycrd[j] = orig_yi.x;
        poly_psw.zcrd[j] = orig_zi.x;
        poly_psw.xcrd_ovrf[j] = orig_xi.y;
        poly_psw.ycrd_ovrf[j] = orig_yi.y;
        poly_psw.zcrd_ovrf[j] = orig_zi.y;
        poly_psw.xalt[j] = orig_xi.x;
        poly_psw.yalt[j] = orig_yi.x;
        poly_psw.zalt[j] = orig_zi.x;
        poly_psw.xalt_ovrf[j] = orig_xi.y;
        poly_psw.yalt_ovrf[j] = orig_yi.y;
        poly_psw.zalt_ovrf[j] = orig_zi.y;
      }
    }
    break;
  default:
    break;
  }
  AtomGraphSynthesis poly_ag = ions.exportAtomGraphSynthesis(std::vector<int>(ncopy, 0));
  poly_ag.upload();
  CellGrid<Tcoord, Tacc,
           Tcalc, Tcoord4> cg(&poly_ps, &poly_ag, 8.0, 0.05, 4, NonbondedTheme::ALL, 1000);
  cg.upload();
  
  // Introduce a perturbation to the official coordinates for the neighbor list to adopt
  switch (stage) {
  case 1:
    addRandomNoise(xrs, poly_psw.xalt, poly_psw.xalt_ovrf, poly_psw.yalt, poly_psw.yalt_ovrf,
                   poly_psw.zalt, poly_psw.zalt_ovrf, poly_ps.getPaddedAtomCount(), 2.5,
                   poly_psw.gpos_scale_f);
    break;
  case 2:
    for (int i = 0; i < poly_psw.system_count; i++) {
      const int llim = poly_psw.atom_starts[i];
      const int hlim = llim + poly_psw.atom_counts[i];
      const double yspan = poly_psw.invu[(i * xfrm_stride) + 4];
      const double zspan = poly_psw.invu[(i * xfrm_stride) + 8];
      for (int j = llim; j < hlim; j++) {
        const double orig_x = hostInt95ToDouble(poly_psw.xcrd[j], poly_psw.xcrd_ovrf[j]) *
                              poly_psw.inv_gpos_scale_f;
        const double next_x = orig_x + xrs->gaussianRandomNumber();
        const double next_y = yspan * (0.34 + (0.16 * xrs->uniformRandomNumber()));
        const double next_z = zspan * (0.34 + (0.16 * xrs->uniformRandomNumber()));
        const int95_t next_xi = hostDoubleToInt95(next_x * poly_psw.gpos_scale);
        const int95_t next_yi = hostDoubleToInt95(next_y * poly_psw.gpos_scale);
        const int95_t next_zi = hostDoubleToInt95(next_z * poly_psw.gpos_scale);
        poly_psw.xalt[j] = next_xi.x;
        poly_psw.yalt[j] = next_yi.x;
        poly_psw.zalt[j] = next_zi.x;
        poly_psw.xalt_ovrf[j] = next_xi.y;
        poly_psw.yalt_ovrf[j] = next_yi.y;
        poly_psw.zalt_ovrf[j] = next_zi.y;
      }
    }
    break;
  default:
    break;
  }
  poly_ps.upload();

  // Update the GPU positions.  The cell grid's BLACK image on the GPU device should then match
  // the BLACK image in the coordinate synthesis, which is consistent between the CPU and GPU.
  // Download the cell grid and compare its BLACK image to that of the coordinate synthesis on the
  // CPU.  After the migration is complete, the CellGrid's coordinate cycle will have ticked
  // forward, making BLACK its current image.
  const CoreKlManager launcher(gpu, poly_ag);
  launchMigration(&cg, poly_ps, launcher);
  CellGridWriter<Tcoord, Tacc, Tcalc, Tcoord4> host_cgw = cg.data();
  cg.download();
  std::vector<uint2> image_indexing_failures;
  const size_t padded_natom = poly_ps.getPaddedAtomCount();
  std::vector<double> nbl_x(padded_natom), nbl_y(padded_natom), nbl_z(padded_natom);
  std::vector<double> syn_x(padded_natom), syn_y(padded_natom), syn_z(padded_natom);
  for (int i = 0; i < host_cgw.total_cell_count; i++) {
    const uint2 ic_bounds = host_cgw.cell_limits_alt[i];
    const uint llim = ic_bounds.x;
    const uint hlim = ic_bounds.x + (ic_bounds.y >> 16);
    const int sys_idx = (ic_bounds.y & 0xffffU);
    const double* umat = &poly_psw.umat_alt[sys_idx * xfrm_stride];
    const double* invu = &poly_psw.invu_alt[sys_idx * xfrm_stride];
    const ullint sys_gdims = host_cgw.system_cell_grids[sys_idx];
    const int cell_na = ((sys_gdims >> 28) & 0xfff);
    const int cell_nb = ((sys_gdims >> 40) & 0xfff);
    const int cell_nc = (sys_gdims >> 52);
    const double dcell_na = cell_na;
    const double dcell_nb = cell_nb;
    const double dcell_nc = cell_nc;
    const int cell_abcs = (sys_gdims & 0xfffffff);
    const int cell_del = i - cell_abcs;
    const int cell_cidx = cell_del / (cell_na * cell_nb);
    const int cell_bidx = (cell_del - (cell_cidx * cell_na * cell_nb)) / cell_na;
    const int cell_aidx = cell_del - (((cell_cidx * cell_nb) + cell_bidx) * cell_na);
    const double da_origin = static_cast<double>(cell_aidx) / static_cast<double>(cell_na);
    const double db_origin = static_cast<double>(cell_bidx) / static_cast<double>(cell_nb);
    const double dc_origin = static_cast<double>(cell_cidx) / static_cast<double>(cell_nc);
    for (uint j = llim; j < hlim; j++) {

      // Place the atom into the primary unit cell image, based on its neighbor list cell and
      // location within that cell.
      const Tcoord4 crdq = host_cgw.image_alt[j];
      const int topl_idx = host_cgw.nonimg_atom_idx_alt[j];
      const uint cg_img_idx = host_cgw.img_atom_idx_alt[topl_idx];
      if (cg_img_idx != j) {
        image_indexing_failures.push_back({ j, cg_img_idx }); 
      }
      double da = crdq.x;
      double db = crdq.y;
      double dc = crdq.z;
      if (tcoord_is_integral) {
        da *= host_cgw.inv_lpos_scale;
        db *= host_cgw.inv_lpos_scale;
        dc *= host_cgw.inv_lpos_scale;
      }
      const double frac_da = (umat[0] * da) + (umat[3] * db) + (umat[6] * dc) + da_origin;
      const double frac_db =                  (umat[4] * db) + (umat[7] * dc) + db_origin;
      const double frac_dc =                                   (umat[8] * dc) + dc_origin;
      const double atom_x = (invu[0] * frac_da) + (invu[3] * frac_db) + (invu[6] * frac_dc);
      const double atom_y =                       (invu[4] * frac_db) + (invu[7] * frac_dc);
      const double atom_z =                                             (invu[8] * frac_dc);
      double sy_da = hostInt95ToDouble(poly_psw.xalt[topl_idx],
                                       poly_psw.xalt_ovrf[topl_idx]) * poly_psw.inv_gpos_scale_f;
      double sy_db = hostInt95ToDouble(poly_psw.yalt[topl_idx],
                                       poly_psw.yalt_ovrf[topl_idx]) * poly_psw.inv_gpos_scale_f;
      double sy_dc = hostInt95ToDouble(poly_psw.zalt[topl_idx],
                                       poly_psw.zalt_ovrf[topl_idx]) * poly_psw.inv_gpos_scale_f;
      imageCoordinates(&sy_da, &sy_db, &sy_dc, umat, invu, poly_psw.unit_cell,
                       ImagingMethod::PRIMARY_UNIT_CELL, poly_psw.gpos_scale_f);
      syn_x[topl_idx] = sy_da;
      syn_y[topl_idx] = sy_db;
      syn_z[topl_idx] = sy_dc;
      nbl_x[topl_idx] = atom_x;
      nbl_y[topl_idx] = atom_y;
      nbl_z[topl_idx] = atom_z;
    }
  }
  check(nbl_x, RelationalOperator::EQUAL, Approx(syn_x).margin(crd_tol), "Cartesian X coordinates "
        "in a neighbor list were not updated in accord with the underlying coordinate synthesis.");
  check(nbl_y, RelationalOperator::EQUAL, Approx(syn_y).margin(crd_tol), "Cartesian Y coordinates "
        "in a neighbor list were not updated in accord with the underlying coordinate synthesis.");
  check(nbl_z, RelationalOperator::EQUAL, Approx(syn_z).margin(crd_tol), "Cartesian Z coordinates "
        "in a neighbor list were not updated in accord with the underlying coordinate synthesis.");
}
#endif

//-------------------------------------------------------------------------------------------------
// Encapsulate the checks on tile interaction coverage, to be applied to interactions between
// atoms in distinct tiles as well as "self" interactions between atoms in the same tile.
//
// Arguments:
//   off_tile_interactions:  Flag to indicate that the interactions with fictitious placeholder
//                           atoms should be included in the tally
//   snd_natom:              The number of sending atoms
//   rcv_natom:              The number of receiving atoms
//   intr:                   A mask of confirmed interactions between the sending and receiving
//                           interactions, computed in terms of the ranged cutoff
//   test_context:           Description of the context of the test
//-------------------------------------------------------------------------------------------------
void tallyTileInteractions(const bool off_tile_interactions, const int snd_natom,
                           const int rcv_natom, const std::vector<std::vector<int>> &intr,
                           const std::string &test_context) {
  check(off_tile_interactions == false, "Off-tile interactions were included in evaluating " +
        std::to_string(snd_natom) + " x " + std::to_string(rcv_natom) + " atoms with warp "
        "size " + std::to_string(warp_size_int) + ".");
  int missing_interactions = 0;
  int double_counting      = 0;
  const bool ignore_diag = (test_context == "tile self interactions");
  for (int i = 0; i < snd_natom; i++) {
    for (int j = 0; j < rcv_natom; j++) {
      missing_interactions += (intr[i][j] == 0 && (ignore_diag == false || i != j));
      double_counting += (intr[i][j] > 1);
    }
  }
  check(missing_interactions, RelationalOperator::EQUAL, 0, "Some interactions between atoms "
        "were not counted in evaluating " + std::to_string(snd_natom) + " x " +
        std::to_string(rcv_natom) + " atoms with warp size " + std::to_string(warp_size_int) +
        " (" + test_context + ").");
  check(double_counting, RelationalOperator::EQUAL, 0, "Some interactions between atoms "
        "were counted multiple times in evaluating " + std::to_string(snd_natom) + " x " +
        std::to_string(rcv_natom) + " atoms with warp size " + std::to_string(warp_size_int) +
        " (" + test_context + ").");
}

//-------------------------------------------------------------------------------------------------
// Test the tile plans used to manage interacting particle sets in the neutral-territory
// decomposition.
//-------------------------------------------------------------------------------------------------
void testTilePlans() {
  const int2 lp = { 320, 256 };
  const TileManager tlmn(lp);
  for (int i = 0; i <= tlmn.getMaximumBatchDegeneracy(); i++) {
    const int snd_deg = ipow(2, i);
    const int snd_natom = warp_size_int / snd_deg;
    const std::vector<int> send_atoms = tlmn.getSendingAtomLayout(snd_deg);

    // Check the sending atom layout
    std::vector<int> chk_send_atoms(warp_size_int);
    int popcon = 0;
    for (int j = 0; j < snd_deg; j++) {
      for (int k = 0; k < snd_natom; k++) {
        chk_send_atoms[popcon] = k;
        popcon++;
      }
    }
    check(send_atoms, RelationalOperator::EQUAL, chk_send_atoms, "The layout of sending atom "
          "relative indices produced by a TileManager does not meet expectations.");

    // Loop over various receiving atom configurations for interactions among distinct cells
    for (int j = 0; j <= tlmn.getMaximumBatchDegeneracy(); j++) {
      const int rcv_deg = ipow(2, j);
      const int rcv_natom = warp_size_int / rcv_deg;
      std::vector<std::vector<int>> intr(snd_natom, std::vector<int>(rcv_natom, 0));
      std::vector<int> recv_atoms = tlmn.getReadAssignments(snd_deg, rcv_deg);
      const int tile_depth = (snd_natom * rcv_natom) / warp_size_int;
      bool off_tile_interactions = false;
      for (int k = 0; k < tile_depth; k++) {

        // Loop over all "threads" of the warp and record the interaction
        for (int m = 0; m < warp_size_int; m++) {
          if (send_atoms[m] < 0 || send_atoms[m] >= snd_natom ||
              recv_atoms[m] < 0 || recv_atoms[m] >= rcv_natom) {
            off_tile_interactions = true;
          }
          else {
            intr[send_atoms[m]][recv_atoms[m]] += 1;
          }
        }

        // Advance the receiving atoms, as if by shuffle instructions in the warp.  This is not
        // done on the last iteration, as the reduction shuffles anticipate the data configuration
        // after the final iteration is calculated, not calculated and moved again.
        if (k < tile_depth - 1) {
          const int zero_idx = recv_atoms[0];
          for (int m = 0; m < warp_bits_mask_int; m++) {
            recv_atoms[m] = recv_atoms[m + 1];
          }
          recv_atoms[warp_bits_mask_int] = zero_idx;
        }
      }
      tallyTileInteractions(off_tile_interactions, snd_natom, rcv_natom, intr, "distinct tiles");
      
      // Apply the reorganization that would prepare for reduction of the accumulated forces.
      // The receiving atom layout, after applying the reduction map, should look like a sending
      // atom layout at the same level of degeneracy.
      std::vector<int> reduce_deck(warp_size_int);
      const std::vector<int> reduce_insr = tlmn.getReductionPreparations(snd_deg, rcv_deg);
      for (int k = 0; k < warp_size_int; k++) {
        reduce_deck[k] = recv_atoms[reduce_insr[k]];
      }
      const std::vector<int> reduce_ansr = tlmn.getSendingAtomLayout(rcv_deg);
      check(reduce_deck, RelationalOperator::EQUAL, reduce_ansr, "The reduction-prepared "
            "receiving atom layout is not suitable for collecting forces.");
    }

    // Loop over various receiving atom configurations for interactions within one cell
    for (int j = 0; j <= tlmn.getMaximumBatchDegeneracy(); j++) {
      const int rcv_deg = ipow(2, j);
      const int rcv_natom = warp_size_int / rcv_deg;
      std::vector<std::vector<int>> intr(snd_natom, std::vector<int>(rcv_natom, 0));
      std::vector<int> recv_atoms = tlmn.getSelfAssignments(snd_deg, rcv_deg);
      int tile_depth = (snd_natom * rcv_natom) / warp_size_int;
      if (i == j) {
        tile_depth /= 2;
      }
      bool off_tile_interactions = false;
      for (int k = 0; k < tile_depth; k++) {

        // Loop over all "threads" of the warp and record the interaction
        for (int m = 0; m < warp_size_int; m++) {
          if (send_atoms[m] < 0 || send_atoms[m] >= snd_natom ||
              recv_atoms[m] < 0 || recv_atoms[m] >= rcv_natom) {
            off_tile_interactions = true;
          }
          else {
            intr[send_atoms[m]][recv_atoms[m]] += 1;
            if (i == j) {

              // Count the pair both ways, so that the basic coverage checking can be applied
              intr[recv_atoms[m]][send_atoms[m]] += 1;
            }
          }
        }

        // If this is the first iteration of the tile, apply the scaling factors, again so that
        // the basic coverage checking can be applied.
        if (i == j && k == 0) {
          const std::vector<float> tscl = tlmn.getThreadScalings(snd_deg);
          for (int m = 0; m < warp_size_int; m++) {
            const int scl_count = round(intr[send_atoms[m]][recv_atoms[m]] * tscl[m]);
            intr[send_atoms[m]][recv_atoms[m]] = scl_count;
          }
        }

        // Advance the receiving atoms, as if by shuffle instructions in the warp
        const int zero_idx = recv_atoms[0];
        for (int m = 0; m < warp_bits_mask_int; m++) {
          recv_atoms[m] = recv_atoms[m + 1];
        }
        recv_atoms[warp_bits_mask_int] = zero_idx;
      }
      tallyTileInteractions(off_tile_interactions, snd_natom, rcv_natom, intr,
                            "tile self interactions");
    }
  }
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {
  
  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  StopWatch timer;
  Xoshiro256ppGenerator xrs(oe.getRandomSeed());
#ifdef STORMM_USE_HPC
  HpcConfig gpu_config(ExceptionResponse::WARN);
  std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> create_to_engage_gpu(1);
#endif

  // Section 1
  section("Forward exclusion list analysis");

  // Section 2
  section("Cell grid creation and mechanics");

  // Section 3
  section("Local exclusion mask analysis");

  // Section 4
  section("Particle-particle interaction spline tables");

  // Section 5
  section("Non-bonded tile mechanics");
  
  // Locate topologies and coordinate files
  const char osc = osSeparator();
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string base_ptl_name = oe.getStormmSourcePath() + osc + "test" + osc + "Potential";
  const std::vector<std::string> system_names = { "trpcage", "trpcage_in_water", "dhfr_cmap",
                                                  "ala_dipeptide", "ahector", "tamavidin",
                                                  "ubiquitin", "drug_example",
                                                  "symmetry_C1_in_water", "symmetry_L1_vs",
                                                  "bromobenzene_vs", "drug_example_vs" };
  TestSystemManager tsm(base_top_name, "top", system_names, base_crd_name, "inpcrd", system_names,
                        ExceptionResponse::SILENT);
  AtomGraph trpi_ag = tsm.exportAtomGraph(0);
  AtomGraph trpw_ag = tsm.exportAtomGraph(1);
  AtomGraph dhfr_ag = tsm.exportAtomGraph(2);
  AtomGraph alad_ag = tsm.exportAtomGraph(3);
  AtomGraph ahec_ag = tsm.exportAtomGraph(4);
  PhaseSpace trpi_ps = tsm.exportPhaseSpace(0);
  PhaseSpace trpw_ps = tsm.exportPhaseSpace(1);
  PhaseSpace dhfr_ps = tsm.exportPhaseSpace(2);
  PhaseSpace alad_ps = tsm.exportPhaseSpace(3);
  PhaseSpace ahec_ps = tsm.exportPhaseSpace(4);

  // Set up the forward exclusion lists and check each topology's exclusions
  section(1);
  const std::vector<AtomGraph*> all_topologies = { &trpi_ag, &trpw_ag, &dhfr_ag, &alad_ag,
                                                   &ahec_ag };
  const std::vector<ForwardExclusionMask> all_femasks = { ForwardExclusionMask(&trpi_ag),
                                                          ForwardExclusionMask(&trpw_ag),
                                                          ForwardExclusionMask(&dhfr_ag),
                                                          ForwardExclusionMask(&alad_ag),
                                                          ForwardExclusionMask(&ahec_ag) };
  const std::vector<int> primary_mask_count_answer = { 99, 99, 140, 11, 145 };
  const std::vector<int> extended_mask_count_answer = { 0, 0, 0, 0, 20 };
  const std::vector<int> secondary_mask_count_answer = { 0, 0, 0, 0, 0 };
  std::vector<int> primary_mask_counts(all_topologies.size());
  std::vector<int> extended_mask_counts(all_topologies.size());
  std::vector<int> secondary_mask_counts(all_topologies.size());

  // Scan the exclusion lists in each topology, making sure that the forward masks cover all
  // of them.  Also check the topology exclusions themselves, in one additional context, to
  // make sure that they contain no self exclusions and that every exclusion is reflexive.
  for (size_t i = 0; i < all_topologies.size(); i++) {
    bool self_exclusion_detected = false;
    std::vector<int2> non_reflexive_exclusions;
    const NonbondedKit<double> nbk = all_topologies[i]->getDoublePrecisionNonbondedKit();
    const std::vector<bool> nb11x_satisfied =
      checkForwardMaskExcl(nbk.natom, nbk.nb11_bounds, nbk.nb11x, all_femasks[i],
                           &self_exclusion_detected, &non_reflexive_exclusions);
    const std::vector<bool> nb12x_satisfied =
      checkForwardMaskExcl(nbk.natom, nbk.nb12_bounds, nbk.nb12x, all_femasks[i],
                           &self_exclusion_detected, &non_reflexive_exclusions);
    const std::vector<bool> nb13x_satisfied =
      checkForwardMaskExcl(nbk.natom, nbk.nb13_bounds, nbk.nb13x, all_femasks[i],
                           &self_exclusion_detected, &non_reflexive_exclusions);
    const std::vector<bool> nb14x_satisfied =
      checkForwardMaskExcl(nbk.natom, nbk.nb14_bounds, nbk.nb14x, all_femasks[i],
                           &self_exclusion_detected, &non_reflexive_exclusions);
    check(self_exclusion_detected == false, "Self exclusions were detected in topology " +
          all_topologies[i]->getFileName() + ".", tsm.getTestingStatus());
    check(non_reflexive_exclusions.size() == 0LLU, "Non-reflexive exclusions were detected in "
          "topology " + all_topologies[i]->getFileName() + ".", tsm.getTestingStatus());
    if (nb11x_satisfied.size() > 0LLU) {
      check(nb11x_satisfied, RelationalOperator::EQUAL,
            std::vector<bool>(nbk.nb11_bounds[nbk.natom], true), "Some 1:1 exclusions in "
            "topology " + all_topologies[i]->getFileName() + " were not satisfied by the forward "
            "mask.", tsm.getTestingStatus());
    }
    if (nb12x_satisfied.size() > 0LLU) {
      check(nb12x_satisfied, RelationalOperator::EQUAL,
            std::vector<bool>(nbk.nb12_bounds[nbk.natom], true), "Some 1:2 exclusions in "
            "topology " + all_topologies[i]->getFileName() + " were not satisfied by the forward "
            "mask.", tsm.getTestingStatus());
    }
    if (nb13x_satisfied.size() > 0LLU) {
      check(nb13x_satisfied, RelationalOperator::EQUAL,
            std::vector<bool>(nbk.nb13_bounds[nbk.natom], true), "Some 1:3 exclusions in "
            "topology " + all_topologies[i]->getFileName() + " were not satisfied by the forward "
            "mask.", tsm.getTestingStatus());
    }
    if (nb14x_satisfied.size() > 0LLU) {
      check(nb14x_satisfied, RelationalOperator::EQUAL,
            std::vector<bool>(nbk.nb14_bounds[nbk.natom], true), "Some 1:4 exclusions in "
            "topology " + all_topologies[i]->getFileName() + " were not satisfied by the forward "
            "mask.", tsm.getTestingStatus());
    }
    primary_mask_counts[i]   = all_femasks[i].getPrimaryMaskCount();
    extended_mask_counts[i]  = all_femasks[i].getExtendedMaskCount();
    secondary_mask_counts[i] = all_femasks[i].getSecondaryMaskCount();
    check((nbk.nb11_bounds[nbk.natom] + nbk.nb12_bounds[nbk.natom] + nbk.nb13_bounds[nbk.natom] +
           nbk.nb14_bounds[nbk.natom]) / 2, RelationalOperator::EQUAL,
          all_femasks[i].getTotalExclusionsCount(), "The total number of exclusions expressed in "
          "topology " + getBaseName(all_topologies[i]->getFileName()) + " is not refelcted in the "
          "corresponding ForwardExclusionsMask.", tsm.getTestingStatus());
  }
  check(primary_mask_counts, RelationalOperator::EQUAL, primary_mask_count_answer, "Unique "
        "primary exclusion mask counts in each system do not match.", tsm.getTestingStatus());
  check(extended_mask_counts, RelationalOperator::EQUAL, extended_mask_count_answer, "Unique "
        "extended exclusion mask counts in each system do not match.", tsm.getTestingStatus());
  check(secondary_mask_counts, RelationalOperator::EQUAL, secondary_mask_count_answer, "Unique "
        "secondary exclusion mask counts in each system do not match.", tsm.getTestingStatus());
  
  // Form a cell grid using some of the periodic systems.  Begin by testing the cell decomposition
  // methods.
  section(2);
  int ncoarsened = 0;
  double most_coarsened = 1.0; 
  for (int i = 4; i <= 25; i++) {
    for (int j = 4; j <= 25; j++) {
      for (int k = 4; k <= 25; k++) {
        const int3 opt_ijk = optimizeCellConfiguration(i, j, k, 3);
        const double orig_cc = i * j * k;
        const double opt_cc  = opt_ijk.x * opt_ijk.y * opt_ijk.z;
        ncoarsened += (opt_cc / orig_cc < 0.90);
        most_coarsened = std::min(opt_cc / orig_cc, most_coarsened);
      }
    }
  }
  check(ncoarsened, RelationalOperator::EQUAL, 658, "The number of cell grid dimensional "
        "combinations resulting in coarser grids (to satisfy FFT dimension requirements) did not "
        "meet expectations.");
  check(most_coarsened, RelationalOperator::EQUAL, Approx(0.786527082).margin(1.0e-5), "The most "
        "significant shift in the cell grid dimensions did not meet expectations.");
  check(computeMigrationRate(5.5, 0.1), RelationalOperator::EQUAL,
        Approx(0.10500884).margin(1.0e-6), "The migratory rate of particles within a 5.5A cell "
        "moving an average of 0.1A per step does not meet expectations.");
  check(computeMigrationRate(6.5, 0.1), RelationalOperator::EQUAL,
        Approx(0.09174435).margin(1.0e-6), "The migratory rate of particles within a 5.5A cell "
        "moving an average of 0.1A per step does not meet expectations.");
  const std::vector<int> pbc_systems = tsm.getQualifyingSystems({ UnitCellType::ORTHORHOMBIC,
                                                                  UnitCellType::TRICLINIC });
  PhaseSpaceSynthesis poly_ps = tsm.exportPhaseSpaceSynthesis(pbc_systems);
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(pbc_systems);
  CellGrid<double, llint, double, double4> cg(poly_ps, poly_ag, 5.0, 0.25, 4, NonbondedTheme::ALL);
  for (size_t i = 0; i < poly_ps.getSystemCount(); i++) {
    checkCellGridPlacements<double, llint, double, double4>(poly_ps, i, cg, 1.0e-8,
                                                            tsm.getTestingStatus());
  }
  PMIGrid pmig(cg, NonbondedTheme::ELECTROSTATIC, 5, PrecisionModel::DOUBLE);
#ifdef STORMM_USE_HPC  
  pmig.prepareWorkUnits(gpu);
#endif

  // Test the CellGrid object's particle migration
  section(2);
  const std::vector<std::string> ion_name = { "fivek_ions" };
  TestSystemManager ions(base_top_name, "top", ion_name, base_crd_name, "inpcrd", ion_name);
  testCellOrigins(ions, std::vector<int>(8, 0), 36, &xrs);
  testCellOrigins(tsm, pbc_systems, 64, &xrs);
  testCpuCellMigration(&cg, &poly_ps, &xrs, tsm.getTestingStatus());
#ifdef STORMM_USE_HPC
  for (int i = 1; i <= 3; i++) {
    for (int j = 1; j <= 2; j++) {
      testGpuCellMigration<double, llint, double, double4>(ions, i, 4.5e-6, j, &xrs, gpu);
      testGpuCellMigration<float, int, float, float4>(ions, i, 4.5e-6, j, &xrs, gpu);
    }
  }
#endif

  // Check the mechanics of the LocalExclusionMask object
  section(3);
  testLocalExclusionMaskSizing();
  testLocalExclusionMaskMechanics(&xrs);
  timer.assignTime(0);
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    testLocalExclusionMask(tsm.getTopologyPointer(i), &timer, tsm.getTestingStatus());
  }
  testLocalExclusionMask(poly_ag, &timer, tsm.getTestingStatus());

  // Check the larger kinase system with local exclusion masks
  const std::vector<std::string> kinase_str(1, "kinase");
  TestSystemManager kinase_tsm(base_top_name, "top", kinase_str, base_crd_name, "inpcrd",
                               kinase_str, ExceptionResponse::SILENT);
  testLocalExclusionMask(kinase_tsm.getTopologyPointer(0), &timer, kinase_tsm.getTestingStatus());
  
  // Test the particle-particle interaction table
  section(4);
  testPPSplineTables<float4>(10.0, 5, BasisFunctions::POLYNOMIAL, &xrs);
  testPPSplineTables<double4>(10.0, 4, BasisFunctions::POLYNOMIAL, &xrs);

  // Test the tile plans
  section(5);
  testTilePlans();
  
  // Print results
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
