// -*-c++-*-
#include "copyright.h"
#include "../../src/Accelerator/core_kernel_manager.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/MolecularMechanics/dynamics.h"
#include "../../src/MolecularMechanics/hpc_dynamics.h"
#include "../../src/MolecularMechanics/hpc_minimization.h"
#include "../../src/MolecularMechanics/kinetic.h"
#include "../../src/MolecularMechanics/minimization.h"
#include "../../src/MolecularMechanics/mm_controls.h"
#include "../../src/Namelists/nml_dynamics.h"
#include "../../src/Namelists/nml_minimize.h"
#include "../../src/Namelists/nml_pppm.h"
#include "../../src/Namelists/nml_precision.h"
#include "../../src/Namelists/nml_random.h"
#include "../../src/Namelists/nml_solvent.h"
#include "../../src/Potential/cacheresource.h"
#include "../../src/Potential/cellgrid.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/hpc_pme_potential.h"
#include "../../src/Potential/hpc_nonbonded_potential.h"
#include "../../src/Potential/hpc_valence_potential.h"
#include "../../src/Potential/ppitable.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/tile_manager.h"
#include "../../src/Random/random.h"
#include "../../src/Structure/hpc_virtual_site_handling.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/implicit_solvent_workspace.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/nonbonded_workunit.h"
#include "../../src/Synthesis/static_mask_synthesis.h"
#include "../../src/Synthesis/valence_workunit.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_intake.h"
#include "../../src/Trajectory/coordinate_copy.h"
#include "../../src/Trajectory/coordinate_intake.h"
#include "../../src/Trajectory/coordinate_series.h"
#include "../../src/Trajectory/hpc_integration.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/thermostat.h"
#include "../../src/Trajectory/trajectory_enumerators.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::energy;
using namespace stormm::mm;
using namespace stormm::namelist;
using namespace stormm::random;
using namespace stormm::stmath;
using namespace stormm::structure;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;
#if (CUDART_VERSION < 13000)
using stormm::data_types::double4_16a;
using stormm::data_types::ulonglong4_16a;
#endif

//-------------------------------------------------------------------------------------------------
// Download an array of values based on an arbitrary starting index and stated length.
//
// Arguments:
//   buffer:     An array to hold the result, allocated by this routine to reside on the CPU host
//   source:     The source array, pointing to memory on the GPU device
//   start_idx:  Stating index within the source array at which to begin the download.  It is
//               trusted that the GPU array has enough memory to accommodate the start (and the
//               stated length n).
//   n:          Trusted length of the arrays on the GPU device
//   desc:       Description of the quantity represented by the integers
//   caller:     Name of the calling function
//-------------------------------------------------------------------------------------------------
template <typename T>
void downloadValues(std::vector<T> *buffer, const T* source, const size_t start_idx,
                    const size_t n, const std::string &desc, const char* caller) {
  buffer->resize(n);
  if (cudaMemcpy(buffer->data(), &source[start_idx], n * sizeof(T), cudaMemcpyDeviceToHost) !=
      cudaSuccess) {
    rtErr("Error in cudaMemcpy (downloading device " + desc + ").", caller);
  }
}

//-------------------------------------------------------------------------------------------------
// Download an array of integer data type accumulators along with an associated array of 32-bit
// integer overflow accumulators.
//
// Arguments:
//   base_acc:   Array to hold the downloaded base accumulators.  This array will be allocated,
//               filled, and returned.
//   ovrf_acc:   Array to hold the downloaded overflow accumulators.  This array will be allocated,
//               filled, and returned.
//   base_devc:  The array of base accumulators, stored on the GPU device
//   ovrf_devc:  The array of base accumulators, stored on the GPU device
//   n:          Trusted length of the arrays on the GPU device
//   desc:       Description of the quantity represented by the integers
//   caller:     Name of the calling function
//-------------------------------------------------------------------------------------------------
template <typename Tbase>
void downloadSplitValues(std::vector<Tbase> *base_acc, std::vector<int> *ovrf_acc,
                         const Tbase* base_devc, const int* ovrf_devc, const size_t n,
                         const std::string &desc, const char* caller) {
  base_acc->resize(n);
  ovrf_acc->resize(n);
  if (cudaMemcpy(base_acc->data(), base_devc, n * sizeof(Tbase), cudaMemcpyDeviceToHost) !=
      cudaSuccess) {
    rtErr("Error in cudaMemcpy (downloading device " + desc + ").", caller);
  }
  if (cudaMemcpy(ovrf_acc->data(), ovrf_devc, n * sizeof(int), cudaMemcpyDeviceToHost) !=
      cudaSuccess) {
    rtErr("Error in cudaMemcpy (downloading device " + desc + ").", caller);
  }
}

//-------------------------------------------------------------------------------------------------
// Scramble the coordinates of a synthesis of systems by applying a slight random perturbation to
// all positions, shifting all positions by a random (and much more significant) value such that
// each system will still tesselate in three dimensions but pack into the neighbor list cells in a
// different way, and swap the coordinates of molecules with identical stoichiometry (e.g. water),
// again at random.
//
// Arguments:
//   poly_ps:  The coordinate synthesis to modify
//   noise:    Multiplier for the amount to perturb all coordinates
//   shift:    Multiplier for the amount to shift all coordinates
//-------------------------------------------------------------------------------------------------
void scrambleSynthesis(PhaseSpaceSynthesis *poly_ps, Xoroshiro128pGenerator *xrs,
                       const double noise = 0.01, const double shift = 11.0) {
  PsSynthesisWriter poly_psw = poly_ps->data();
  std::vector<bool> scrambled(poly_psw.system_count, false);
  for (int i = 0; i < poly_psw.system_count; i++) {
    if (scrambled[i]) {
      continue;
    }
    const AtomGraph* agi = poly_ps->getSystemTopologyPointer(i);
    for (int j = i + 1; j < poly_psw.system_count; j++) {
      if (scrambled[j]) {
        continue;
      }
      const AtomGraph* agj = poly_ps->getSystemTopologyPointer(j);

      // This filter is not absolute, but sufficient to identify identical topologies within a
      // controlled test case.
      if (agi->getAtomCount() != agj->getAtomCount() ||
          agi->getResidueCount() != agj->getResidueCount() ||
          agi->getMoleculeCount() != agj->getMoleculeCount()) {
        continue;
      }
      
      // Perturb and translate any system after the first in a group.
      const double xtrans = xrs->gaussianRandomNumber() * shift;
      const double ytrans = xrs->gaussianRandomNumber() * shift;
      const double ztrans = xrs->gaussianRandomNumber() * shift;
      const int j_llim = poly_psw.atom_starts[j];
      const int j_hlim = j_llim + poly_psw.atom_counts[j];
      for (int k = j_llim; k < j_hlim; k++) {
        const double xpert = xrs->gaussianRandomNumber() * noise;
        const double ypert = xrs->gaussianRandomNumber() * noise;
        const double zpert = xrs->gaussianRandomNumber() * noise;
        poly_psw.xcrd[k] += xpert + xtrans;
        poly_psw.ycrd[k] += ypert + ytrans;
        poly_psw.zcrd[k] += zpert + ztrans;
        poly_psw.xalt[k] += xpert + xtrans;
        poly_psw.yalt[k] += ypert + ytrans;
        poly_psw.zalt[k] += zpert + ztrans;
      }

      // The two systems share the same topology and can therefore be considered identical.
      // Scramble the positions of equivalent molecules in this (j^th) system, and any subsequent
      // systems that match the i^th system.
      scrambled[j] = true;
      const ChemicalDetailsKit cdk = agi->getChemicalDetailsKit();
      std::vector<int> mol_list;
      std::vector<bool> swapped(cdk.nmol, false);
      for (int k = 0; k < cdk.nmol; k++) {
        if (swapped[k]) {
          continue;
        }
        mol_list.resize(1);
        mol_list[0] = k;
        const int nk_atom = cdk.mol_limits[k + 1] - cdk.mol_limits[k];
        for (int m = k + 1; m < cdk.nmol; m++) {
          if (cdk.mol_limits[m + 1] - cdk.mol_limits[m] == nk_atom) {
            bool match = true;
            for (int p = 0; p < nk_atom; p++) {
              const int katom = cdk.mol_contents[cdk.mol_limits[k] + p];
              const int matom = cdk.mol_contents[cdk.mol_limits[m] + p];
              match = (match && (cdk.z_numbers[katom] == cdk.z_numbers[matom]));
            }
            if (match) {
              mol_list.push_back(m);
              swapped[m] = true;
            }
          }
        }
        const size_t mol_list_size = mol_list.size();
        const double nswap_mol = mol_list_size;
        for (size_t m = 0; m < mol_list_size; m++) {

          // Choose two molecules from the list at random and swap their positions.
          size_t mol_a = mol_list[static_cast<int>(xrs->uniformRandomNumber() * nswap_mol)];
          size_t mol_b = mol_list[static_cast<int>(xrs->uniformRandomNumber() * nswap_mol)];
          mol_a -= (mol_a == mol_list_size);
          mol_b -= (mol_b == mol_list_size);
          for (int p = 0; p < nk_atom; p++) {
            const int atom_a = cdk.mol_contents[cdk.mol_limits[mol_list[mol_a]] + p];
            const int atom_b = cdk.mol_contents[cdk.mol_limits[mol_list[mol_b]] + p];
            const int synth_a = j_llim + atom_a;
            const int synth_b = j_llim + atom_b;
            std::swap(poly_psw.xcrd[synth_a], poly_psw.xcrd[synth_b]);
            std::swap(poly_psw.ycrd[synth_a], poly_psw.ycrd[synth_b]);
            std::swap(poly_psw.zcrd[synth_a], poly_psw.zcrd[synth_b]);
            std::swap(poly_psw.xalt[synth_a], poly_psw.xalt[synth_b]);
            std::swap(poly_psw.yalt[synth_a], poly_psw.yalt[synth_b]);
            std::swap(poly_psw.zalt[synth_a], poly_psw.zalt[synth_b]);
          }
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Download forces from a cell grid for a specific system, packaged into a convenient, real-valued
// coordinate object for a single system with memory stored on the CPU host.  This routine
// encapsualtes the procedure for reuse.
//
// Arguments:
//   gpu_cgr:       A pointer to the cell grid's abstract, oriented to the proper point in the
//                  white black time cycle.  The abstract contains additional pointers which
//                  address data on the GPU device.
//   cpu_cgr:       A pointer to the cell grid's abstract, oriented to the proper point in the
//                  white black time cycle.  The abstract contains additional pointers which
//                  address data on the CPU host.
//   poly_ps:       The coordinate synthesis from which the cell grid is derived
//   system_index:  Index of the system to package, as numbered within the synthesis
//   data_kind:     The type of data that is to be extracted from the GPU memory
//   caller:        Name of the calling function
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4> CoordinateFrame
packageCellGridDataFromDevice(const CellGridReader<Tcoord, Tacc, Tcalc, Tcoord4> *cpu_cgr,
                              const CellGridReader<Tcoord, Tacc, Tcalc, Tcoord4> *gpu_cgr,
                              const PhaseSpaceSynthesis &poly_ps, const int system_index,
                              const char* caller) {
  const bool tacc_is_long = (std::type_index(typeid(Tacc)).hash_code() == llint_type_index);

  // The GPU requires a download using HPC utilities
  const int natom = poly_ps.getAtomCount(system_index);
  const int synth_start_idx = poly_ps.getAtomOffset(system_index);
  CoordinateFrame result(natom, poly_ps.getUnitCellType(), HybridFormat::HOST_ONLY);
  const int imin_chain = cpu_cgr->system_chain_bounds[system_index];
  const int imax_chain = cpu_cgr->system_chain_bounds[system_index + 1];
  const int img_start_idx = cpu_cgr->chain_limits[imin_chain];
  const int img_final_idx = cpu_cgr->chain_limits[imax_chain];
  const size_t img_len = img_final_idx - img_start_idx;
  std::vector<Tacc> img_xdata, img_ydata, img_zdata;
  std::vector<int>  img_xdata_ovrf, img_ydata_ovrf, img_zdata_ovrf;
  std::vector<uint> img_atom_indices;
  downloadSplitValues(&img_xdata, &img_xdata_ovrf, &gpu_cgr->xfrc[img_start_idx],
                      &gpu_cgr->xfrc_ovrf[img_start_idx], img_len, "X forces", caller);
  downloadSplitValues(&img_ydata, &img_ydata_ovrf, &gpu_cgr->yfrc[img_start_idx],
                      &gpu_cgr->yfrc_ovrf[img_start_idx], img_len, "Y forces", caller);
  downloadSplitValues(&img_zdata, &img_zdata_ovrf, &gpu_cgr->zfrc[img_start_idx],
                      &gpu_cgr->zfrc_ovrf[img_start_idx], img_len, "Z forces", caller);
  downloadValues(&img_atom_indices, gpu_cgr->img_atom_idx, synth_start_idx, natom,
                 "Image atom indices", "checkCellGridPositions");
  CoordinateFrameWriter resultw = result.data();
  const double inv_scl = cpu_cgr->inv_frc_scale;
  for (int j = 0; j < resultw.natom; j++) {
    const uint img_idx = img_atom_indices[j] - img_start_idx;
    if (tacc_is_long) {
      resultw.xcrd[j] = hostInt95ToDouble(img_xdata[img_idx], img_xdata_ovrf[img_idx]) * inv_scl;
      resultw.ycrd[j] = hostInt95ToDouble(img_ydata[img_idx], img_ydata_ovrf[img_idx]) * inv_scl;
      resultw.zcrd[j] = hostInt95ToDouble(img_zdata[img_idx], img_zdata_ovrf[img_idx]) * inv_scl;
    }
    else {
      resultw.xcrd[j] = hostInt63ToDouble(img_xdata[img_idx], img_xdata_ovrf[img_idx]) * inv_scl;
      resultw.ycrd[j] = hostInt63ToDouble(img_ydata[img_idx], img_ydata_ovrf[img_idx]) * inv_scl;
      resultw.zcrd[j] = hostInt63ToDouble(img_zdata[img_idx], img_zdata_ovrf[img_idx]) * inv_scl;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Check the particle properties resident in the neighbor list against those found in the
// non-bonded topology synthesis abstract.
//
// Arguments:
//   cg:        The neighbor list, assumed to be of theme "ALL" (both Lennard-Jones and
//              electrostatic properties)
//   poly_nbk:  Non-bonded parameters from the original (reference) topology synthesis
//   do_tests:  Indicate whether testing is possible
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcalc2, typename Tcoord4>
void checkNeighborListProperties(const CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> &cg,
                                 const SyNonbondedKit<Tcalc, Tcalc2> &poly_nbk,
                                 const TestPriority do_tests) {

  // Verify the neighbor list property contents
  switch (cg.getTheme()) {
  case NonbondedTheme::ELECTROSTATIC:
  case NonbondedTheme::VAN_DER_WAALS:
    rtErr("The neighbor list must contain both electrostatic and van-der Waals properties.  "
          "Instead, it contains " + getEnumerationName(cg.getTheme()) + " particle properties.",
          "checkNeighborListProperties");
    break;
  case NonbondedTheme::ALL:
    break;
  }
  std::vector<CoordinateCycle> all_cyc = { CoordinateCycle::WHITE, CoordinateCycle::BLACK };
  const bool tcoord_is_real = isFloatingPointScalarType<Tcoord>();
  for (size_t pos = 0; pos < all_cyc.size(); pos++) {
    const CellGridReader<Tcoord, Tacc, Tcalc, Tcoord4> cgr = cg.data(all_cyc[pos]);
    for (int sysidx = 0; sysidx < poly_nbk.nsys; sysidx++) {
      const ullint gdims = cgr.system_cell_grids[sysidx];
      const int cell_start = (gdims & 0xfffffff);
      const int cell_na = ((gdims >> 28) & 0xfff);
      const int cell_nb = ((gdims >> 40) & 0xfff);
      const int cell_nc = ((gdims >> 52) & 0xfff);
      const int total_cells = cell_na * cell_nb * cell_nc;
      const int natom = poly_nbk.atom_counts[sysidx];
      std::vector<double> q_from_cg(natom, 0.0);
      std::vector<int> lj_from_cg(natom, -1);
      std::vector<bool> q_covered(natom, false), lj_covered(natom, false);
      const int ag_llim = poly_nbk.atom_offsets[sysidx];
      const int ag_hlim = ag_llim + natom;
      for (int i = 0; i < total_cells; i++) {
        const uint2 clims = cgr.cell_limits[cell_start + i];
        const uint c_llim = clims.x;
        const uint c_hlim = clims.x + (clims.y >> 16);
        for (int j = c_llim; j < c_hlim; j++) {
          const int topl_idx = cgr.nonimg_atom_idx[j];
          const Tcoord prop = cgr.image[j].w;
          const int q_idx = sourceIndex<Tcoord>(NonbondedTheme::ELECTROSTATIC, NonbondedTheme::ALL,
                                                cgr.image[j].w, tcoord_is_real);
          q_from_cg[topl_idx - ag_llim] = poly_nbk.q_params[q_idx];
          lj_from_cg[topl_idx - ag_llim] = sourceIndex<Tcoord>(NonbondedTheme::VAN_DER_WAALS,
                                                               NonbondedTheme::ALL,
                                                               cgr.image[j].w, tcoord_is_real);
          q_covered[topl_idx - ag_llim] = true;
          lj_covered[topl_idx - ag_llim] = true;
        }
      }
      std::vector<double> q_from_ag(natom, 0.0);
      std::vector<int> lj_from_ag(natom, -1);
      bool all_q_covered = true;
      bool all_lj_covered = true;
      for (int i = ag_llim; i < ag_hlim; i++) {
        q_from_ag[i - ag_llim] = poly_nbk.charge[i];
        lj_from_ag[i - ag_llim] = poly_nbk.lj_idx[i];
        all_q_covered = (all_q_covered && q_covered[i - ag_llim]);
        all_lj_covered = (all_lj_covered && lj_covered[i - ag_llim]);
      }
      check(all_q_covered, "A total of " + std::to_string(natom - sum<int>(q_covered)) +
            " particles' charges (of " + std::to_string(natom) + ") were not found in the "
            "neighbor list.", do_tests);
      check(all_lj_covered, "A total of " + std::to_string(natom - sum<int>(lj_covered)) +
            " particles' Lennard-Jones indices (of " + std::to_string(natom) + ") were not found "
            "in the neighbor list.", do_tests);
      check(q_from_cg, RelationalOperator::EQUAL, q_from_ag, "Charges obtained by decoding table "
            "indices in a cell grid do not match those in the underlying topology synthesis.  "
            "System index: " + std::to_string(sysidx) + ".  Coordinate cycle: " +
            getEnumerationName(all_cyc[pos]) + ".", do_tests);
      check(lj_from_cg, RelationalOperator::EQUAL, lj_from_ag, "Lennard-Jones table indices "
            "decoded from a cell grid do not match those in the underlying topology synthesis.  "
            "System index: " + std::to_string(sysidx) + ".  Coordinate cycle: " +
            getEnumerationName(all_cyc[pos]) + ".", do_tests);
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Compare the forces stored on the CPU host and on the GPU device levels of a single cell grid,
// system by system.  Arguments follow from packageCellGridForcesFromDevice(), above, in addition
// to:
//
// Arguments:
//   tol:       Tolerance for a successful comparison
//   do_tests:  Indicate whether testing is feasible
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void checkCellGridForces(const CellGridReader<Tcoord, Tacc, Tcalc, Tcoord4> *gpu_cgr,
                         const CellGridReader<Tcoord, Tacc, Tcalc, Tcoord4> *cpu_cgr,
                         const PhaseSpaceSynthesis &poly_ps, const double tol,
                         const std::string &desc, const TestPriority do_tests) {
  const bool tacc_is_long = (std::type_index(typeid(Tacc)).hash_code() == llint_type_index);

  // Check the results for each system
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    CoordinateFrame cpu_frc_i(poly_ps.getAtomCount(i), poly_ps.getUnitCellType(),
                              HybridFormat::HOST_ONLY);
    CoordinateFrameWriter cpu_frc_iw = cpu_frc_i.data();
    for (int j = 0; j < cpu_frc_iw.natom; j++) {
      const int syn_idx = poly_ps.getAtomOffset(i) + j;
      const uint img_idx = cpu_cgr->img_atom_idx[syn_idx];
      if (tacc_is_long) {
        cpu_frc_iw.xcrd[j] = hostInt95ToDouble(cpu_cgr->xfrc[img_idx],
                                               cpu_cgr->xfrc_ovrf[img_idx]) *
                             cpu_cgr->inv_frc_scale;
        cpu_frc_iw.ycrd[j] = hostInt95ToDouble(cpu_cgr->yfrc[img_idx],
                                               cpu_cgr->yfrc_ovrf[img_idx]) *
                             cpu_cgr->inv_frc_scale;
        cpu_frc_iw.zcrd[j] = hostInt95ToDouble(cpu_cgr->zfrc[img_idx],
                                               cpu_cgr->zfrc_ovrf[img_idx]) *
                             cpu_cgr->inv_frc_scale;
      }
      else {
        cpu_frc_iw.xcrd[j] = hostInt63ToDouble(cpu_cgr->xfrc[img_idx],
                                               cpu_cgr->xfrc_ovrf[img_idx]) *
                             cpu_cgr->inv_frc_scale;
        cpu_frc_iw.ycrd[j] = hostInt63ToDouble(cpu_cgr->yfrc[img_idx],
                                               cpu_cgr->yfrc_ovrf[img_idx]) *
                             cpu_cgr->inv_frc_scale;
        cpu_frc_iw.zcrd[j] = hostInt63ToDouble(cpu_cgr->zfrc[img_idx],
                                               cpu_cgr->zfrc_ovrf[img_idx]) *
                             cpu_cgr->inv_frc_scale;
      }
    }

    // The GPU requires a download using HPC utilities
    const CoordinateFrame gpu_frc_i = packageCellGridDataFromDevice(cpu_cgr, gpu_cgr, poly_ps, i,
                                                                    "checkCellGridForces");
    const CoordinateFrameReader gpu_frc_ir = gpu_frc_i.data();

    // Compare results along each Cartesian axis.
    const std::vector<CartesianDimension> dims = { CartesianDimension::X, CartesianDimension::Y,
                                                   CartesianDimension::Z };
    for (size_t j = 0; j < 3; j++) {
      const std::vector<double> cpu_val = cpu_frc_i.getCoordinateHandle(dims[j])->readHost();
      const std::vector<double> gpu_val = gpu_frc_i.getCoordinateHandle(dims[j])->readHost();
      check(gpu_val, RelationalOperator::EQUAL, Approx(cpu_val).margin(tol), "Non-bonded forces "
            "computed for system " + std::to_string(i) + ", based on topology " +
            getBaseName(poly_ps.getSystemTopologyPointer(i)->getFileName()) + ", do not match "
            "along the " + getEnumerationName(dims[j]) + " axis.  Precision model: " +
            getStormmScalarTypeName<Tcalc>() + ".  Situation: " + desc + ".",
            do_tests);
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Check the cell grid's particle positions, comparing locations on the GPU device to those on the
// CPU host.  Descriptions of input parameters follow from checkCellGridForces(), above.
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void checkCellGridPositions(const CellGridReader<Tcoord, Tacc, Tcalc, Tcoord4> *gpu_cgr,
                            const CellGridReader<Tcoord, Tacc, Tcalc, Tcoord4> *cpu_cgr,
                            const PhaseSpaceSynthesis &poly_ps, const double tol,
                            const TestPriority do_tests) {
  const bool tcoord_is_long = (std::type_index(typeid(Tcoord)).hash_code() == llint_type_index ||
                               std::type_index(typeid(Tcoord)).hash_code() == double_type_index);
  const bool tcoord_is_int = isSignedIntegralScalarType<Tcoord>();
  for (int system_index = 0; system_index < poly_ps.getSystemCount(); system_index++) {
    const int natom = poly_ps.getAtomCount(system_index);
    const int synth_start_idx = poly_ps.getAtomOffset(system_index);
    const int imin_chain = cpu_cgr->system_chain_bounds[system_index];
    const int imax_chain = cpu_cgr->system_chain_bounds[system_index + 1];
    const int img_start_idx = cpu_cgr->chain_limits[imin_chain];
    const int img_final_idx = cpu_cgr->chain_limits[imax_chain];
    const int img_len = img_final_idx - img_start_idx;
    std::vector<uint> img_atom_indices(natom);
    std::vector<Tcoord4> img_data;
    std::vector<double> gpu_xdata(img_len), gpu_ydata(img_len), gpu_zdata(img_len);
    std::vector<double> cpu_xdata(img_len), cpu_ydata(img_len), cpu_zdata(img_len);
    downloadValues(&img_data, gpu_cgr->image, img_start_idx, img_len, "XYZ-P image coordinates",
                   "checkCellGridPositions");
    downloadValues(&img_atom_indices, gpu_cgr->img_atom_idx, synth_start_idx, natom,
                   "Image atom indices", "checkCellGridPositions");
    if (tcoord_is_int) {
      const double inv_scale = cpu_cgr->inv_lpos_scale;
      for (int i = 0; i < img_len; i++) {
        gpu_xdata[i] = static_cast<double>(img_data[i].x) * inv_scale;
        gpu_ydata[i] = static_cast<double>(img_data[i].y) * inv_scale;
        gpu_zdata[i] = static_cast<double>(img_data[i].z) * inv_scale;
        const Tcoord4 cpu_quad = cpu_cgr->image[img_start_idx + i];
        cpu_xdata[i] = static_cast<double>(cpu_quad.x) * inv_scale;
        cpu_ydata[i] = static_cast<double>(cpu_quad.y) * inv_scale;
        cpu_zdata[i] = static_cast<double>(cpu_quad.z) * inv_scale;
      }
    }
    else {
      for (int i = 0; i < img_len; i++) {
        gpu_xdata[i] = img_data[i].x;
        gpu_ydata[i] = img_data[i].y;
        gpu_zdata[i] = img_data[i].z;
        const Tcoord4 cpu_quad = cpu_cgr->image[img_start_idx + i];
        cpu_xdata[i] = cpu_quad.x;
        cpu_ydata[i] = cpu_quad.y;
        cpu_zdata[i] = cpu_quad.z;
      }
    }
  
    // Compare results along each Cartesian axis.
    const std::vector<CartesianDimension> dims = { CartesianDimension::X, CartesianDimension::Y,
                                                   CartesianDimension::Z };
    for (int i = 0; i < 3; i++) {
      std::vector<double> cpu_val(natom), gpu_val(natom);
      for (int j = 0; j < natom; j++) {
        const int cpu_idx = cpu_cgr->img_atom_idx[synth_start_idx + j] - img_start_idx;
        const int gpu_idx = img_atom_indices[j] - img_start_idx;
        if (i == 0) {
          cpu_val[j] = cpu_xdata[cpu_idx];
          gpu_val[j] = gpu_xdata[gpu_idx];
        }
        else if (i == 1) {
          cpu_val[j] = cpu_ydata[cpu_idx];
          gpu_val[j] = gpu_ydata[gpu_idx];
        }
        else if (i == 2) {
          cpu_val[j] = cpu_zdata[cpu_idx];
          gpu_val[j] = gpu_zdata[gpu_idx];
        }
      }
      check(gpu_val, RelationalOperator::EQUAL, Approx(cpu_val).margin(tol), "Local neighbor list "
            "positions computed for system " + std::to_string(system_index) + ", based on "
            "topology " +
            getBaseName(poly_ps.getSystemTopologyPointer(system_index)->getFileName()) +
            ", do not match along the " + getEnumerationName(dims[i]) + " axis.", do_tests);
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Compare data accumulated in the coordinate synthesis between the CPU host and GPU device.
//
// Arguments:
//   poly_ps:          The coordinates, including forces on both the CPU host and GPU device
//   cg:               The cell grid, which may be needed to augment the set of forces on the GPU
//   gpu:              Details of the GPU on which the data is stored (and which will perform the
//                     transfer)
//   tol:              Tolerance for a successful comparison
//   data_kind:        Type of coordinates to compare: particle positions, velocities, or forces
//   take_next_image:  Indicate that the developing image, as opposed to the current image, should
//                     be compared between he GPU and CPU data sets 
//   do_tests:         Indicate whether testing is feasible
//   valence_prec:     Precision with which valence interactions and phase space adjustments were
//                     performed
//   nonbonded_prec:   Precision with which nonbonded neighbor list interactions were performed
//   situation:        Describe the point in the time step cycle at which comparsions were made
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void checkSynthesisMatch(const PhaseSpaceSynthesis &poly_ps,
                         const CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> &cg, const GpuDetails &gpu,
                         const double tol, const TrajectoryKind data_kind,
                         const bool take_next_image, const TestPriority do_tests,
                         const PrecisionModel valence_prec, const PrecisionModel nonbonded_prec,
                         const std::string &situation) {
  const CoordinateCycle orientation = (take_next_image) ?
                                      getNextCyclePosition(poly_ps.getCyclePosition()) :
                                      poly_ps.getCyclePosition();
  const PsSynthesisReader cpu_syn = poly_ps.data(orientation);
  const PsSynthesisReader gpu_syn = poly_ps.data(orientation, HybridTargetLevel::DEVICE);
  const CellGridReader<Tcoord, Tacc, Tcalc, Tcoord4> cpu_cgr = cg.data();
  const CellGridReader<Tcoord, Tacc, Tcalc, Tcoord4> gpu_cgr = cg.data(HybridTargetLevel::DEVICE);
  for (int i = 0; i < cpu_syn.system_count; i++) {
    CoordinateFrame cpu_rslt_i(cpu_syn.atom_counts[i], cpu_syn.unit_cell,
                               HybridFormat::HOST_MOUNTED);
    CoordinateFrame gpu_rslt_i(cpu_syn.atom_counts[i], cpu_syn.unit_cell,
                               HybridFormat::HOST_MOUNTED);
    coordCopy(&cpu_rslt_i, poly_ps, i, data_kind, orientation);
    coordCopy(&gpu_rslt_i, poly_ps, i, data_kind, orientation, HybridTargetLevel::HOST,
              HybridTargetLevel::DEVICE, gpu);
    CoordinateFrame gpu_cg_i;
    switch (data_kind) {
    case TrajectoryKind::POSITIONS:
    case TrajectoryKind::VELOCITIES:
      break;
    case TrajectoryKind::FORCES:
      gpu_cg_i = packageCellGridDataFromDevice(&cpu_cgr, &gpu_cgr, poly_ps, i,
                                               "checkSynthesisForces");
      break;
    }

    // Compare results along each Cartesian axis.
    const std::vector<CartesianDimension> dims = { CartesianDimension::X, CartesianDimension::Y,
                                                   CartesianDimension::Z };
    const size_t natom = cpu_rslt_i.getAtomCount();
    for (size_t j = 0; j < 3; j++) {
      const std::vector<double> cpu_sy_val = cpu_rslt_i.getCoordinateHandle(dims[j])->readHost();
      const std::vector<double> gpu_sy_val = gpu_rslt_i.getCoordinateHandle(dims[j])->readHost();
      std::vector<double> gpu_total_val(natom);
      std::string desc;
      switch (data_kind) {
      case TrajectoryKind::POSITIONS:
        desc = "Positions";
        for (size_t k = 0; k < natom; k++) {
          gpu_total_val[k] = gpu_sy_val[k];
        }
        break;
      case TrajectoryKind::VELOCITIES:
        desc = "Velocities";
        for (size_t k = 0; k < natom; k++) {
          gpu_total_val[k] = gpu_sy_val[k];
        }
        break;
      case TrajectoryKind::FORCES:
        {
          desc = "Combined forces";
          const std::vector<double> gpu_cg_val = gpu_cg_i.getCoordinateHandle(dims[j])->readHost();
          for (size_t k = 0; k < natom; k++) {
            gpu_total_val[k] = gpu_sy_val[k] + gpu_cg_val[k];
          }
        }
        break;
      }
      check(gpu_total_val, RelationalOperator::EQUAL, Approx(cpu_sy_val).margin(tol),
            desc + " computed for system " + std::to_string(i) + ", based on topology " +
            getBaseName(poly_ps.getSystemTopologyPointer(i)->getFileName()) + ", do not match "
            "along the " + getEnumerationName(dims[j]) + " axis.  Valence interactions were "
            "computed in " + getEnumerationName(valence_prec) + ", neighbor list interactions "
            "in " + getEnumerationName(nonbonded_prec) + ".  Situation: " + situation + ".",
            do_tests);
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Check the charge density mapped to the particle-mesh interaction grid, or the electrostatic
// potential computed after the convolution, comparing results on the GPU device to those on the
// CPU host.
//
// Arguments:
//   pmigr_gpu:  Reader for results on the GPU device
//   pmigr_cpu:  Reader for results on the CPU hostm the reference calculation
//   tol:        Tolerance for a successful comparison
//   desc:       Description of the comparison to be made
//   do_tests:   Indicate whether testing is possible
//-------------------------------------------------------------------------------------------------
void checkMeshgridMatch(const PMIGridReader &pmigr_gpu, const PMIGridReader &pmigr_cpu,
                        const double tol, const std::string &desc, const TestPriority do_tests) {
  std::vector<uint4> gpu_dims;
  downloadValues<uint4>(&gpu_dims, pmigr_gpu.dims, 0, pmigr_cpu.nsys, "particle-mesh interaction "
                        "grid dimensions", "checkMeshgridMatch");
  for (int sys_idx = 0; sys_idx < pmigr_cpu.nsys; sys_idx++) {
    const uint4 t_dims = pmigr_cpu.dims[sys_idx];
    const std::vector<uint> cpu_dimv = { t_dims.x, t_dims.y, t_dims.z, t_dims.w };
    const std::vector<uint> gpu_dimv = { gpu_dims[sys_idx].x, gpu_dims[sys_idx].y,
                                         gpu_dims[sys_idx].z, gpu_dims[sys_idx].w };
    check(gpu_dimv, RelationalOperator::EQUAL, cpu_dimv, "Grid dimensions and array placement "
          "on the GPU do not match those on the CPU.");
    const uint npts = t_dims.x * t_dims.y * t_dims.z;
    switch (pmigr_cpu.mode) {
    case PrecisionModel::DOUBLE:
      {
        std::vector<double> cpu_result(t_dims.x * t_dims.y * t_dims.z);
        std::vector<double> gpu_result(t_dims.x * t_dims.y * t_dims.z);
        downloadValues(&gpu_result, pmigr_gpu.ddata, t_dims.w, npts,
                       "mesh data in " + getEnumerationName(pmigr_cpu.mode), "checkMeshgridMatch");
        for (uint i = 0; i < npts; i++) {
          cpu_result[i] = pmigr_cpu.ddata[t_dims.w + i];
        }
        check(gpu_result, RelationalOperator::EQUAL, Approx(cpu_result).margin(tol), "Mesh grid "
              "data for system " + std::to_string(sys_idx) + " obtained by the GPU does not match "
              "that obtained by the CPU.  Precision of the mesh grid: " +
              getEnumerationName(pmigr_cpu.mode) + ".  Process examined: " + desc + ".", do_tests);
      }
      break;
    case PrecisionModel::SINGLE:
      {
        std::vector<float> cpu_result(t_dims.x * t_dims.y * t_dims.z);
        std::vector<float> gpu_result(t_dims.x * t_dims.y * t_dims.z);
        downloadValues(&gpu_result, pmigr_gpu.fdata, t_dims.w, npts,
                       "mesh data in " + getEnumerationName(pmigr_cpu.mode), "checkMeshgridMatch");
        for (uint i = 0; i < npts; i++) {
          cpu_result[i] = pmigr_cpu.fdata[t_dims.w + i];
        }
        check(gpu_result, RelationalOperator::EQUAL, Approx(cpu_result).margin(tol), "Mesh grid "
              "data for system " + std::to_string(sys_idx) + " obtained by the GPU does not match "
              "that obtained by the CPU.  Precision of the mesh grid: " +
              getEnumerationName(pmigr_cpu.mode) + ".  Process examined: " + desc + ".", do_tests);
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Check the forward FFT and convolution results obtained on the GPU device by comparing them to
// results obtained on the CPU host.
//
// Arguments:
//   conv_gpu:  Framework for the GPU-based convolutions
//   conv_cpu:  Framework for the CPU-based calculations
//   pmig:      Particle-mesh interaction grids holding the actual data upon which each convolution
//              rests.  This also indicates the precision model for such calculations.
//   tol:       Tolerance for a successful comparison
//   desc:      Description of the precise moment at which the frequency-space data is tested
//   do_tests:  Indicate whether testing is possible
//-------------------------------------------------------------------------------------------------
void checkConvolutionMatch(const ConvolutionManager &conv_gpu, const ConvolutionManager &conv_cpu,
                           const PMIGrid &pmig, const double tol, const std::string &desc,
                           const TestPriority do_tests) {
  const ConvolutionReader<double, double2> dconvr_gpu = conv_gpu.dpData();
  const ConvolutionReader<float, float2> fconvr_gpu = conv_gpu.spData();
  const ConvolutionReader<double, double2> dconvr_cpu = conv_cpu.dpData();
  const ConvolutionReader<float, float2> fconvr_cpu = conv_cpu.spData();
  const PMIGridReader pmigr_cpu = pmig.data();
  const PMIGridReader pmigr_gpu = pmig.data(HybridTargetLevel::DEVICE);
  for (int batch_idx = 0; batch_idx < dconvr_gpu.fft_ops->size(); batch_idx++) {
    const FFTStage& tgrp = dconvr_gpu.fft_ops->at(batch_idx);

    // The list of systems assigned to each GPU batch will be available on the host in all cases.
    // In the convolution manager used by the CPU host, there may not be a valid list or it may
    // differ from the list used by the GPU.
    const std::vector<int> system_list = conv_gpu.getFFTGroupSystemList(batch_idx);
    const size_t npts = tgrp.getFrequencyCount(0) * tgrp.getFrequencyCount(1) *
                        tgrp.getFrequencyCount(2);
    std::vector<double> cpu_real(npts), gpu_real(npts), cpu_imag(npts), gpu_imag(npts);
    for (size_t i = 0; i < tgrp.getBatchCount(); i++) {

      // The dimensions of each mesh grid will be consistent between the CPU and GPU calculations.
      const uint4 i_dims = pmigr_cpu.dims[system_list[i]];
      switch (pmigr_cpu.mode) {
      case PrecisionModel::DOUBLE:
        {
          std::vector<double2> cpu_fft(npts), gpu_fft;
          switch (pmigr_cpu.fftm) {
          case FFTMode::IN_PLACE:
            {
              std::string ext_rep;
              if (tgrp.getBatchCount() > 1) {
                ext_rep = ", member " + std::to_string(i) + " of batch group " +
                          std::to_string(batch_idx);
              }
              else {
                ext_rep = "";
              }
              downloadValues(&gpu_fft, reinterpret_cast<const double2*>(pmigr_gpu.ddata),
                             i_dims.w / 2, npts, "FFT data for system " +
                             std::to_string(system_list[i]) + ext_rep, "checkConvolutionMatch");
            }
            break;
          case FFTMode::OUT_OF_PLACE:
            downloadValues(&gpu_fft, tgrp.getFrequencyData<double2>(),
                           i * tgrp.getFrequencyBatchStride(), npts, "FFT data for system " +
                           std::to_string(system_list[i]), "checkConvolutionMatch");
            break;
          }
          for (size_t i = 0; i < npts; i++) {
            cpu_real[i] = cpu_fft[i].x;
            cpu_imag[i] = cpu_fft[i].y;
            gpu_real[i] = gpu_fft[i].x;
            gpu_imag[i] = gpu_fft[i].y;
          }
        }
        break;
      case PrecisionModel::SINGLE:
        {
          std::vector<float2> cpu_fft(npts), gpu_fft;
          switch (pmigr_cpu.fftm) {
          case FFTMode::IN_PLACE:
            {
              std::string ext_rep;
              if (tgrp.getBatchCount() > 1) {
                ext_rep = ", member " + std::to_string(i) + " of batch group " +
                          std::to_string(batch_idx);
              }
              else {
                ext_rep = "";
              }
              downloadValues(&gpu_fft, reinterpret_cast<const float2*>(pmigr_gpu.fdata),
                             i_dims.w / 2, npts, "FFT data for system " +
                             std::to_string(system_list[i]) + ext_rep, "checkConvolutionMatch");
            }
            break;
          case FFTMode::OUT_OF_PLACE:
            downloadValues(&gpu_fft, tgrp.getFrequencyData<float2>(),
                           i * tgrp.getFrequencyBatchStride(), npts, "FFT data for system " +
                           std::to_string(system_list[i]), "checkConvolutionMatch");
            break;
          }
          std::vector<double> cpu_real(npts), gpu_real(npts), cpu_imag(npts), gpu_imag(npts);
          for (size_t i = 0; i < npts; i++) {
            cpu_real[i] = cpu_fft[i].x;
            cpu_imag[i] = cpu_fft[i].y;
            gpu_real[i] = gpu_fft[i].x;
            gpu_imag[i] = gpu_fft[i].y;
          }
        }
        break;
      }
      check(gpu_real, RelationalOperator::EQUAL, Approx(cpu_real).margin(tol), "The real part of "
            "frequency-space data generated on the GPU does not match that on the CPU.  Precision "
            "of the mesh grid: " + getEnumerationName(pmigr_cpu.mode) + ".  Particular "
            "situation: " + desc + ".", do_tests);
      check(gpu_imag, RelationalOperator::EQUAL, Approx(cpu_imag).margin(tol), "The imaginary "
            "part of frequency-space data generated on the GPU does not match that on the CPU.  "
            "Precision of the mesh grid: " + getEnumerationName(pmigr_cpu.mode) + ".  Particular "
            "situation: " + desc + ".", do_tests);
    }
  }
}
                           
//-------------------------------------------------------------------------------------------------
// Carry out multiple steps of PME dynamics in a synthesis of systems using the GPU while
// replicating the process on the CPU.  Initial velocities are set to zero or initialized according
// to a Maxwell distribution on the CPU, then ported to the GPU.
//
// Arguments:
//   tsm:            A collection of test systems, ready to unpack into new coordinate and topology
//                   syntheses
//   dyncon:         Mock user input from the &dynamics namelist control block
//   preccon:        Mock user input from the &precision namelist control block
//   rngcon:         Mock user input from the &random namelist control block
//   pmecon:         Mock user input from the &pppm namelist control block
//   gpu:            Details of the GPU that will carry out the HPC calculations
//   ngbr_frc_tol:   Tolerance for a successful comparison of forces derived from the neighbor list
//                   on the CPU host versus those calculated on the GPU device
//   qgrid_tol:      Tolerance for a successful comparison of charge density grid values on the CPU
//                   host versus those calculated on the GPU device
//   ugrid_tol:      Tolerance for a successful comparison of electrostatic potential grid values
//                   on the CPU host versus those calculated on the GPU device
//   total_frc_tol:  Tolerance for a successful comparison of total forces calculated by the CPU
//                   host versus those calculated by the GPU device
//   vel_tol:        Tolerance for a successful comparison of particle velocities calculated on the
//                   CPU host and those calcualted on the GPU device
//   pos_tol:        Tolerance for a successful comparison of particle positions calculated on the
//                   CPU host and those calcualted on the GPU device
//-------------------------------------------------------------------------------------------------
void stepwisePMELaboratory(const TestSystemManager &tsm, const DynamicsControls &dyncon,
                           const PrecisionControls &preccon, const RandomControls &rngcon,
                           const PPPMControls &pmecon, const GpuDetails &gpu,
                           const double ngbr_frc_tol, const double qgrid_tol,
                           const double ugrid_tol, const double total_frc_tol,
                           const double vel_tol, const double pos_tol) {
  const std::vector<UnitCellType> pbcs = { UnitCellType::ORTHORHOMBIC, UnitCellType::TRICLINIC };
  const std::vector<int> pbcs_idx = tsm.getQualifyingSystems(pbcs);
  
  // Prepare basic resources.  The cutoffs specified in &dynamics and &pppm namelist control blocks
  // may conflict, so resolve the conflicts based on namelist priority and whether the "user" has
  // specified the relevant keywords in one namelist or another.
  Xoroshiro128pGenerator xrs;
  PhaseSpaceSynthesis poly_ps = tsm.exportPhaseSpaceSynthesis(pbcs_idx, 0.00,
                                                              rngcon.getRandomSeed(),
                                                              preccon.getGlobalPosScalingBits(),
                                                              preccon.getVelocityScalingBits(),
                                                              preccon.getForceScalingBits());
  scrambleSynthesis(&poly_ps, &xrs);
  ScoreCard sc(poly_ps.getSystemCount(), dyncon.getStepCount(), 32);
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(pbcs_idx);
  Thermostat tst(poly_ag, ThermostatKind::NONE);
  tst.setGeometryConstraints(dyncon.constrainGeometry());
  tst.setRattleTolerance(dyncon.getRattleTolerance());
  tst.setTimeStep(dyncon.getTimeStep());
  LocalExclusionMask lem(poly_ag);
  MolecularMechanicsControls mmctrl_fe(dyncon, pmecon);
  MolecularMechanicsControls mmctrl_fx(dyncon, pmecon);
  
  // Upload the main topological and coordinate data, along with basic resources.
  poly_ag.upload();
  poly_ps.upload();
  tst.upload();
  lem.upload();

  // Create abstracts to basic resources on the GPU device
  const PrecisionModel val_prec = preccon.getValenceMethod();
  const PrecisionModel nb_prec = preccon.getNonbondedMethod();
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  const CoordinateCycle poly_ps_next_stage = getNextCyclePosition(poly_ps.getCyclePosition());
  PsSynthesisWriter poly_psw = poly_ps.data(devc);
  PsSynthesisWriter poly_psw_alt = poly_ps.data(poly_ps_next_stage, devc);
  PsSynthesisReader poly_psr(poly_psw);
  PsSynthesisReader poly_psr_alt(poly_psw_alt);
  PsSynthesisBorders pssb = poly_ps.borders(devc);
  PsSynthesisBorders pssb_alt = poly_ps.borders(poly_ps_next_stage, devc);
  const LocalExclusionMaskReader lemr = lem.data(devc);
  ScoreCardWriter scw = sc.data(devc);
  ThermostatWriter<double> d_tstw = tst.dpData(devc);
  ThermostatWriter<float> f_tstw = tst.spData(devc);
  ThermostatWriter<double> host_dtstw = tst.dpData();
  ThermostatWriter<float> host_ftstw = tst.spData();
  ThermostatReader<double> d_tstr(d_tstw);
  ThermostatReader<float> f_tstr(f_tstw);
  ThermostatReader<double> host_dtstr(host_dtstw);
  ThermostatReader<float> host_ftstr(host_ftstw);

  // Create host-oriented abstracts for the tandem process
  PsSynthesisWriter host_poly_psw = poly_ps.data();
  PsSynthesisWriter host_poly_psw_alt = poly_ps.data(poly_ps_next_stage);
  PsSynthesisReader host_poly_psr(host_poly_psw);
  PsSynthesisReader host_poly_psr_alt(host_poly_psw_alt);
  PsSynthesisBorders host_pssb = poly_ps.borders();
  PsSynthesisBorders host_pssb_alt = poly_ps.borders(poly_ps_next_stage);
  const LocalExclusionMaskReader host_lemr = lem.data();
  ScoreCardWriter host_scw = sc.data();

  // Obtain topology abstracts for both precision models on the GPU.
  const SyNonbondedKit<double,
                       double2> dpoly_nbk = poly_ag.getDoublePrecisionNonbondedKit(devc);
  const SyNonbondedKit<float,
                       float2> fpoly_nbk = poly_ag.getSinglePrecisionNonbondedKit(devc);
  const SyValenceKit<double> dpoly_vk = poly_ag.getDoublePrecisionValenceKit(devc);
  const SyValenceKit<float> fpoly_vk = poly_ag.getSinglePrecisionValenceKit(devc);
  const SyAtomUpdateKit<double,
                        double2,
                        double4_16a> dpoly_auk = poly_ag.getDoublePrecisionAtomUpdateKit(devc);
  const SyAtomUpdateKit<float,
                        float2,
                        float4> fpoly_auk = poly_ag.getSinglePrecisionAtomUpdateKit(devc);
  const SyRestraintKit<double,
                       double2,
                       double4_16a> dpoly_rk = poly_ag.getDoublePrecisionRestraintKit(devc);
  const SyRestraintKit<float,
                       float2,
                       float4> fpoly_rk = poly_ag.getSinglePrecisionRestraintKit(devc);

  // Obtain topology abstracts for both precision models on the CPU host.
  const SyNonbondedKit<double, double2> host_dpoly_nbk = poly_ag.getDoublePrecisionNonbondedKit();
  const SyNonbondedKit<float, float2> host_fpoly_nbk = poly_ag.getSinglePrecisionNonbondedKit();
  const SyValenceKit<double> host_dpoly_vk = poly_ag.getDoublePrecisionValenceKit();
  const SyValenceKit<float> host_fpoly_vk = poly_ag.getSinglePrecisionValenceKit();
  const SyAtomUpdateKit<double,
                        double2,
                        double4_16a> host_dpoly_auk = poly_ag.getDoublePrecisionAtomUpdateKit();
  const SyAtomUpdateKit<float,
                        float2, float4> host_fpoly_auk = poly_ag.getSinglePrecisionAtomUpdateKit();
  const SyRestraintKit<double,
                       double2,
                       double4_16a> host_dpoly_rk = poly_ag.getDoublePrecisionRestraintKit();
  const SyRestraintKit<float,
                       float2, float4> host_fpoly_rk = poly_ag.getSinglePrecisionRestraintKit();
  
  // Create launch parameters for kernel execution.
  const CoreKlManager launcher(gpu, poly_ag);

  // Read additional information from the &dynamics namelist.
  const int nscm = dyncon.getCenterOfMassMotionPurgeFrequency();

  // Set up both cell grids, even though only one will be used in testing.  This will simplify the
  // mock dynamics loop.
  CellGrid<double, llint, double, double4_16a> d_cg(poly_ps, poly_ag, mmctrl_fe.getLongestCutoff(),
                                                    0.02, 4, NonbondedTheme::ALL);
  const int2 pair_fe_lp = launcher.getPMEPairsKernelDims(nb_prec, nb_prec, NeighborListKind::MONO,
                                                         d_cg.getTinyBoxPresence(),
                                                         EvaluateForce::YES, EvaluateEnergy::YES,
                                                         ClashResponse::NONE);
  const int2 pair_fx_lp = launcher.getPMEPairsKernelDims(nb_prec, nb_prec, NeighborListKind::MONO,
                                                         d_cg.getTinyBoxPresence(),
                                                         EvaluateForce::YES, EvaluateEnergy::NO,
                                                         ClashResponse::NONE);
  const int2 vale_fx_lp = launcher.getValenceKernelDims(val_prec, EvaluateForce::YES,
                                                        EvaluateEnergy::NO,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::MOVE_PARTICLES,
                                                        ClashResponse::NONE);
  const int2 vale_fe_lp = launcher.getValenceKernelDims(val_prec, EvaluateForce::YES,
                                                        EvaluateEnergy::YES,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::MOVE_PARTICLES,
                                                        ClashResponse::NONE);
  const int2 intg_vvi_lp = launcher.getIntegrationKernelDims(val_prec, AccumulationMethod::SPLIT,
                                                             IntegrationStage::VELOCITY_ADVANCE);
  const int2 intg_vc_lp = launcher.getIntegrationKernelDims(val_prec, AccumulationMethod::SPLIT,
                                                            IntegrationStage::VELOCITY_CONSTRAINT);
  const int2 intg_ke_lp = launcher.getIntegrationKernelDims(val_prec, AccumulationMethod::SPLIT,
                                                            IntegrationStage::CALC_KINETIC);
  const int2 intg_vvii_lp = launcher.getIntegrationKernelDims(val_prec, AccumulationMethod::SPLIT,
                                                              IntegrationStage::POSITION_ADVANCE);
  const int2 intg_gc_lp = launcher.getIntegrationKernelDims(val_prec, AccumulationMethod::SPLIT,
                                                            IntegrationStage::GEOMETRY_CONSTRAINT);
  const int2 migr_one_lp = launcher.getMigrationKernelDims(nb_prec, NeighborListKind::MONO, 1,
                                                           poly_ps.getGlobalPositionBits(),
                                                           d_cg.getTotalChainCount());
  const int2 migr_two_lp = launcher.getMigrationKernelDims(nb_prec, NeighborListKind::MONO, 2,
                                                           poly_ps.getGlobalPositionBits(),
                                                           d_cg.getTotalChainCount());
  const int2 vs_xfer_lp = launcher.getVirtualSiteKernelDims(val_prec,
                                                            VirtualSiteActivity::TRANSMIT_FORCES);
  const int2 vs_place_lp = launcher.getVirtualSiteKernelDims(val_prec,
                                                             VirtualSiteActivity::PLACEMENT);
  CellGrid<float, int, float, float4> f_cg(poly_ps, poly_ag, mmctrl_fe.getLongestCutoff(), 0.02, 4,
                                           NonbondedTheme::ALL);

  // Check the neighbor list to ensure that all atoms have the correct properties
  switch (nb_prec) {
  case PrecisionModel::DOUBLE:
    checkNeighborListProperties(d_cg, host_dpoly_nbk, tsm.getTestingStatus());
    break;
  case PrecisionModel::SINGLE:
    checkNeighborListProperties(f_cg, host_fpoly_nbk, tsm.getTestingStatus());
    break;
  }

  // Detect a tiny unit cell, one with only four neighbor list cells along any one dimension.  This
  // will determine which call to the GPU-enabled PME pairs kernel launcher to make, so that
  // kernels designed to handle tiny boxes can be put to work (making use of the unit cell borders
  // abstract from the coordinate synthesis).
  const bool has_tiny_box = (d_cg.getTinyBoxPresence() == TinyBoxPresence::YES);

  // Detect the presence of virtual sites in any of the systems, which will necessitate additional
  // kernel launches or function calls in the staged time step.
  const bool has_virtual_sites = (poly_ag.getVirtualSiteCount() > 0);
  const CoordinateCycle cg_next_stage = getNextCyclePosition(d_cg.getCyclePosition());
  d_cg.checkViability();
  f_cg.checkViability();

  // Create supporting structures for the particle-mesh and convolution operations
  PMIGrid d_pmig(&d_cg, NonbondedTheme::ELECTROSTATIC, 4, PrecisionModel::DOUBLE,
                 FFTMode::OUT_OF_PLACE, 60);
  PMIGrid f_pmig(&f_cg, NonbondedTheme::ELECTROSTATIC, 4, PrecisionModel::SINGLE,
                 FFTMode::OUT_OF_PLACE, 28);
  PMIGridAccumulator host_dpmacc = d_pmig.fpData();
  PMIGridAccumulator d_pmacc = d_pmig.fpData(devc);
  PMIGridWriter host_dpmwrt = d_pmig.data();
  PMIGridWriter d_pmwrt = d_pmig.data(devc);
  PMIGridReader host_dpmrdr(host_dpmwrt);
  PMIGridReader d_pmrdr(d_pmwrt);
  PMIGridAccumulator host_fpmacc = f_pmig.fpData();
  PMIGridAccumulator f_pmacc = f_pmig.fpData(devc);
  PMIGridWriter host_fpmwrt = f_pmig.data();
  PMIGridWriter f_pmwrt = f_pmig.data(devc);
  PMIGridReader host_fpmrdr(host_fpmwrt);
  PMIGridReader f_pmrdr(f_pmwrt);
  ConvolutionManager d_cvol(&d_pmig, pmecon.getEwaldCoefficient(), gpu);
  ConvolutionManager host_dcvol(&d_pmig, pmecon.getEwaldCoefficient());
  ConvolutionManager f_cvol(&f_pmig, pmecon.getEwaldCoefficient(), gpu);
  ConvolutionManager host_fcvol(&f_pmig, pmecon.getEwaldCoefficient());
  ConvolutionWriter host_dcvolw = host_dcvol.dpData();
  ConvolutionWriter host_fcvolw = host_fcvol.spData();
  ConvolutionWriter d_cvolw = d_cvol.dpData();
  ConvolutionWriter f_cvolw = f_cvol.spData();
  
  // Calculate additional launch parameters for kernels involving particle-mesh interactions.
  size_t ct_tmat;
  bool use_overflow;
  switch (nb_prec) {
  case PrecisionModel::DOUBLE:
    ct_tmat = double_type_index;
    use_overflow = d_pmig.useOverflowAccumulation();
    break;
  case PrecisionModel::SINGLE:
    ct_tmat = float_type_index;
    use_overflow = f_pmig.useOverflowAccumulation();
    break;
  }
  const int2 qmap_lp = launcher.getDensityMappingKernelDims(pmecon.getDensityMappingMethod(),
                                                            nb_prec, nb_prec, use_overflow,
                                                            ct_tmat,
                                                            pmecon.getInterpolationOrder());
  const int2 fintrp_lp = launcher.getForceGatheringKernelDims(QMapMethod::GENERAL_PURPOSE, nb_prec,
                                                              ct_tmat,
                                                              pmecon.getInterpolationOrder());
  
  // Upload the neighbor list data and other components
  d_cg.upload();
  f_cg.upload();
  d_pmig.upload();
  f_pmig.upload();
  d_cvol.upload();
  f_cvol.upload();

  // Create abstracts to data on the GPU device.
  CellGridWriter<double, llint, double, double4_16a> dd_cgw = d_cg.data(devc);
  CellGridWriter<double, llint, double, double4_16a> dd_cgw_alt = d_cg.data(cg_next_stage, devc);
  CellGridWriter<float, int, float, float4> ff_cgw = f_cg.data(devc);
  CellGridWriter<float, int, float, float4> ff_cgw_alt = f_cg.data(cg_next_stage, devc);
  CellGridReader<double, llint, double, double4_16a> dd_cgr = d_cg.data(devc);
  CellGridReader<double, llint, double, double4_16a> dd_cgr_alt = d_cg.data(cg_next_stage, devc);
  CellGridReader<float, int, float, float4> ff_cgr = f_cg.data(devc);
  CellGridReader<float, int, float, float4> ff_cgr_alt = f_cg.data(cg_next_stage, devc);
  CellGridWriter<void, void, void, void> dv_cgw = d_cg.templateFreeData(devc);
  CellGridWriter<void, void, void, void> dv_cgw_alt = d_cg.templateFreeData(cg_next_stage, devc);
  CellGridWriter<void, void, void, void> fv_cgw = f_cg.templateFreeData(devc);
  CellGridWriter<void, void, void, void> fv_cgw_alt = f_cg.templateFreeData(cg_next_stage, devc);
  CellGridReader<void, void, void, void> dv_cgr(dv_cgw);
  CellGridReader<void, void, void, void> dv_cgr_alt(dv_cgw_alt);
  CellGridReader<void, void, void, void> fv_cgr(fv_cgw);
  CellGridReader<void, void, void, void> fv_cgr_alt(fv_cgw_alt);
  CellOriginsReader d_corg = d_cg.getRulers(devc);
  CellOriginsReader d_corg_alt = d_cg.getRulers(cg_next_stage, devc);
  CellOriginsReader f_corg = f_cg.getRulers(devc);
  CellOriginsReader f_corg_alt = f_cg.getRulers(cg_next_stage, devc);
  CellOriginsReader host_dcorg = d_cg.getRulers();
  CellOriginsReader host_dcorg_alt = d_cg.getRulers(cg_next_stage);
  CellOriginsReader host_fcorg = f_cg.getRulers();
  CellOriginsReader host_fcorg_alt = f_cg.getRulers(cg_next_stage);
  
  // Create abstracts to data on the CPU host for calculations that will track and validate
  // dynamics on the GPU.
  CellGridWriter<double, llint, double, double4_16a> host_dcgw = d_cg.data();
  CellGridWriter<float, int, float, float4> host_fcgw = f_cg.data();
  CellGridWriter<double, llint, double, double4_16a> host_dcgw_alt = d_cg.data(cg_next_stage);
  CellGridWriter<float, int, float, float4> host_fcgw_alt = f_cg.data(cg_next_stage);
  CellGridReader<double, llint, double, double4_16a> host_dcgr(host_dcgw);
  CellGridReader<float, int, float, float4> host_fcgr(host_fcgw);
  CellGridReader<double, llint, double, double4_16a> host_dcgr_alt(host_dcgw_alt);
  CellGridReader<float, int, float, float4> host_fcgr_alt(host_fcgw_alt);
  CellGridWriter<void, void, void, void> host_dcgw_v = d_cg.templateFreeData();
  CellGridWriter<void, void, void, void> host_fcgw_v = f_cg.templateFreeData();
  CellGridWriter<void, void, void, void> host_dcgw_v_alt = d_cg.templateFreeData(cg_next_stage);
  CellGridWriter<void, void, void, void> host_fcgw_v_alt = f_cg.templateFreeData(cg_next_stage);
  CellGridReader<void, void, void, void> host_dcgr_v(host_dcgw_v);
  CellGridReader<void, void, void, void> host_fcgr_v(host_fcgw_v);
  CellGridReader<void, void, void, void> host_dcgr_v_alt(host_dcgw_v_alt);
  CellGridReader<void, void, void, void> host_fcgr_v_alt(host_fcgw_v_alt);

  // Prepare additional resources for CPU and GPU dynamics.
  int log_tab_bits;
  switch (preccon.getNonbondedMethod()) {
  case PrecisionModel::DOUBLE:
    log_tab_bits = 7;
    break;
  case PrecisionModel::SINGLE:
    log_tab_bits = 5;
    break;
  }
  PPITable nrg_tab(NonbondedTheme::ELECTROSTATIC, BasisFunctions::MIXED_FRACTIONS,
                   TableIndexing::SQUARED_ARG, mmctrl_fe.getElectrostaticCutoff(), 0.0,
                   pmecon.getDirectSumTolerance(), log_tab_bits);
  const size_t ct_ngbr_crd = (nb_prec == PrecisionModel::DOUBLE) ? double_type_index :
                                                                   float_type_index;
  mmctrl_fe.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                                  ClashResponse::NONE, VwuGoal::ACCUMULATE, val_prec, nb_prec,
                                  pmecon.getDensityMappingMethod(), nb_prec, ct_ngbr_crd, 5,
                                  NeighborListKind::MONO, d_cg.getTinyBoxPresence(), poly_ag);
  mmctrl_fx.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::NO,
                                  ClashResponse::NONE, VwuGoal::ACCUMULATE, val_prec, nb_prec,
                                  pmecon.getDensityMappingMethod(), nb_prec, ct_ngbr_crd, 5,
                                  NeighborListKind::MONO, d_cg.getTinyBoxPresence(), poly_ag);
  TileManager tlmn_fe(pair_fe_lp);
  TileManager tlmn_fx(pair_fx_lp);
  CacheResource vale_fe_cache(vale_fe_lp.x, maximum_valence_work_unit_atoms);
  CacheResource vale_fx_cache(vale_fx_lp.x, maximum_valence_work_unit_atoms);
  CacheResource nonb_fe_cache(pair_fe_lp.x, small_block_max_atoms);
  CacheResource nonb_fx_cache(pair_fx_lp.x, small_block_max_atoms);
  CacheResource vsite_cache(std::max(vs_xfer_lp.x, vs_place_lp.x),
                            std::max(vs_xfer_lp.y, vs_place_lp.y));

  // Upload the additional resources to the GPU
  nrg_tab.upload();
  mmctrl_fe.upload();
  mmctrl_fx.upload();
  tlmn_fe.upload();
  tlmn_fx.upload();
  
  // Obtain abstracts to the additional resources, as necessary, on the CPU host or GPU device.
  // The particle-particle interaction table (PPITable) is not used to accelerate calculations on
  // the host, nor are the tile management table (TileManager) or molecular mechanics control
  // object (MolecularMechanicsControls, used to synchronize GPU work units).
  const PPIKit<double, double4_16a> dppi_direct = nrg_tab.dpData(devc);
  const PPIKit<float, float4> fppi_direct = nrg_tab.spData(devc);
  MMControlKit<double> d_ctrl_fx = mmctrl_fx.dpData(devc);
  MMControlKit<double> d_ctrl_fe = mmctrl_fe.dpData(devc);
  MMControlKit<float> f_ctrl_fx = mmctrl_fx.spData(devc);
  MMControlKit<float> f_ctrl_fe = mmctrl_fe.spData(devc);
  TilePlan tlpn_fe = tlmn_fe.data(devc);
  TilePlan tlpn_fx = tlmn_fx.data(devc);
  CacheResourceKit<double> d_vale_fe_res = vale_fe_cache.dpData(devc);
  CacheResourceKit<double> d_vale_fx_res = vale_fx_cache.dpData(devc);
  CacheResourceKit<float> f_vale_fe_res = vale_fe_cache.spData(devc);
  CacheResourceKit<float> f_vale_fx_res = vale_fx_cache.spData(devc);
  CacheResourceKit<double> dvs_work_res = vsite_cache.dpData(devc);
  CacheResourceKit<float> fvs_work_res = vsite_cache.spData(devc);
  
  // To check GPU dynamics, break down the dynamics cycle and check stage by stage.
  const int ntpr = dyncon.getDiagnosticPrintFrequency();
  for (int step_idx = 0; step_idx < dyncon.getStepCount(); step_idx++) {
    const bool on_energy_step = (ntpr > 0) ? (step_idx % ntpr == 0) : false;
    MMControlKit<double> *d_mmctrl_ptr;
    MMControlKit<float> *f_mmctrl_ptr;
    CacheResourceKit<double> *d_vale_res_ptr;
    CacheResourceKit<float> *f_vale_res_ptr;
    ScoreCardWriter *scw_ptr;
    TilePlan *tile_plan_ptr;
    EvaluateEnergy eval_nrg;
    int2 pair_lp, vale_lp;
    if (on_energy_step) {

      // Initialize the energy tracking on both the CPU host and GPU device
      sc.initialize(devc, gpu);
      sc.initialize();

      // Set pointers to resources for force and energy computations
      d_mmctrl_ptr = &d_ctrl_fe;
      f_mmctrl_ptr = &f_ctrl_fe;
      d_vale_res_ptr = &d_vale_fe_res;
      f_vale_res_ptr = &f_vale_fe_res;
      scw_ptr = &scw;
      eval_nrg = EvaluateEnergy::YES;
      pair_lp = pair_fe_lp;
      vale_lp = vale_fe_lp;
      tile_plan_ptr = &tlpn_fe;
    }
    else {

      // Set pointers to resources for force computations
      d_mmctrl_ptr = &d_ctrl_fx;
      f_mmctrl_ptr = &f_ctrl_fx;
      d_vale_res_ptr = &d_vale_fx_res;
      f_vale_res_ptr = &f_vale_fx_res;
      scw_ptr = nullptr;
      eval_nrg = EvaluateEnergy::NO;
      pair_lp = pair_fx_lp;
      vale_lp = vale_fx_lp;
      tile_plan_ptr = &tlpn_fx;
    }

    // Lay out pointers which can address the proper object, in the proper point of the black and
    // white time cycle.
    PsSynthesisWriter *crdw_ptr, *host_crdw_ptr;
    PsSynthesisReader *host_crdr_ptr;
    PsSynthesisBorders *borders_ptr, *host_borders_ptr;
    CellGridWriter<double, llint, double, double4_16a> *dd_cgw_ptr, *host_dcgw_ptr;
    CellGridWriter<float, int, float, float4> *ff_cgw_ptr, *host_fcgw_ptr;
    CellGridReader<double, llint, double, double4_16a> *dd_cgr_ptr;
    CellGridReader<void, void, void, void> *dv_cgr_ptr, *fv_cgr_ptr;
    CellGridReader<float, int, float, float4> *ff_cgr_ptr;
    CellGridWriter<void, void, void, void> *dv_cgw_ptr, *fv_cgw_ptr;
    CellGridWriter<void, void, void, void> *host_dv_cgw_ptr, *host_fv_cgw_ptr;
    CellGridReader<double, llint, double, double4_16a> *host_dcgr_ptr;
    CellGridReader<float, int, float, float4> *host_fcgr_ptr;
    CellOriginsReader *dcorg_ptr, *fcorg_ptr;
    if (step_idx & 0x1) {
      crdw_ptr = &poly_psw_alt;
      borders_ptr = &pssb_alt;
      host_borders_ptr = &host_pssb_alt;
      host_crdw_ptr = &host_poly_psw_alt;
      host_crdr_ptr = &host_poly_psr_alt;
      dd_cgw_ptr = &dd_cgw_alt;
      ff_cgw_ptr = &ff_cgw_alt;
      dd_cgr_ptr = &dd_cgr_alt;
      ff_cgr_ptr = &ff_cgr_alt;
      dv_cgw_ptr = &dv_cgw_alt;
      fv_cgw_ptr = &fv_cgw_alt;
      dv_cgr_ptr = &dv_cgr_alt;
      fv_cgr_ptr = &fv_cgr_alt;
      host_dv_cgw_ptr = &host_dcgw_v_alt;
      host_fv_cgw_ptr = &host_fcgw_v_alt;
      host_dcgw_ptr = &host_dcgw_alt;
      host_fcgw_ptr = &host_fcgw_alt;
      host_dcgr_ptr = &host_dcgr_alt;
      host_fcgr_ptr = &host_fcgr_alt;
      dcorg_ptr = &d_corg_alt;
      fcorg_ptr = &f_corg_alt;
    }
    else {
      crdw_ptr = &poly_psw;
      borders_ptr = &pssb;
      host_borders_ptr = &host_pssb;
      host_crdw_ptr = &host_poly_psw;
      host_crdr_ptr = &host_poly_psr;
      dd_cgw_ptr = &dd_cgw;
      ff_cgw_ptr = &ff_cgw;
      dd_cgr_ptr = &dd_cgr;
      ff_cgr_ptr = &ff_cgr;
      dv_cgw_ptr = &dv_cgw;
      fv_cgw_ptr = &fv_cgw;
      dv_cgr_ptr = &dv_cgr;
      fv_cgr_ptr = &fv_cgr;
      host_dv_cgw_ptr = &host_dcgw_v;
      host_fv_cgw_ptr = &host_fcgw_v;
      host_dcgw_ptr = &host_dcgw;
      host_fcgw_ptr = &host_fcgw;
      host_dcgr_ptr = &host_dcgr;
      host_fcgr_ptr = &host_fcgr;
      dcorg_ptr = &d_corg;
      fcorg_ptr = &f_corg;
    }

    // Initialize all forces on both the CPU host and GPU device
    poly_ps.initializeForces();
    d_cg.initializeForces();
    f_cg.initializeForces();
    poly_ps.initializeForces(gpu, devc);
    d_cg.initializeForces(devc, gpu);
    f_cg.initializeForces(devc, gpu);
    
    // Compute the non-bonded particle-particle interactions.  
    switch (nb_prec) {
    case PrecisionModel::DOUBLE:
      {
        // The molecular mechanics controls should have the block counter for PME neighbor list
        // interactions set to the number of blocks in the launch grid.
        //checkMMControls(d_mm_ctrl_ptr, 
        
        // GPU process: particle-particle pairwise interactions
        if (has_tiny_box) {
          launchPMEPairs(dpoly_nbk, lemr, dppi_direct, *borders_ptr, dd_cgw_ptr, tile_plan_ptr,
                         &scw, d_mmctrl_ptr, EvaluateForce::YES, eval_nrg, pair_lp, 0.0, 0.0);
        }
        else {
          launchPMEPairs(dpoly_nbk, lemr, dppi_direct, dd_cgw_ptr, tile_plan_ptr,
                         &scw, d_mmctrl_ptr, EvaluateForce::YES, eval_nrg, pair_lp, 0.0, 0.0);
        }

        // CPU process: particle-particle pairwise interactions
        evaluateParticleParticleEnergy<double, llint,
                                       double, double2, double4_16a>(host_dv_cgw_ptr, &host_scw,
                                                                     host_poly_psr, host_dpoly_nbk,
                                                                     host_lemr,
                                                                     mmctrl_fe.getLongestCutoff(),
                                                                     pmecon.getEwaldCoefficient(),
                                                                     dyncon.getVdwSummation(),
                                                                     EvaluateForce::YES,
                                                                     NonbondedTheme::ALL);

        // Check progress
        checkCellGridForces<double, llint, double, double4_16a>(dd_cgr_ptr, host_dcgr_ptr, poly_ps,
                                                                ngbr_frc_tol, "particle-particle "
                                                                "forces evaluated",
                                                                tsm.getTestingStatus());
        
        // GPU process: charge mapping for particle-mesh interactions
        switch (pmecon.getDensityMappingMethod()) {
        case QMapMethod::ACC_SHARED:
        case QMapMethod::AUTOMATIC:
          launchShrAccDensityKernel(&d_pmwrt, use_overflow, d_mmctrl_ptr, *dv_cgr_ptr,
                                    double_type_index, dpoly_nbk, qmap_lp);
          break;
        case QMapMethod::GENERAL_PURPOSE:
          launchPMIGridInitialization(&d_pmacc, gpu);
          launchGenPrpDensityKernel(&d_pmacc, *dv_cgr_ptr, double_type_index, dpoly_nbk, qmap_lp);
          launchPMIGridRealConversion(&d_pmwrt, d_pmacc, gpu);
          break;
        }
        
        // CPU process: charge mapping for particle-mesh interactions
        d_pmig.initialize();
        mapDensity<double, llint,
                   double, double2, double4_16a>(&host_dpmacc, &host_dpmwrt, *host_dcgr_ptr,
                                                 host_dpoly_nbk);

        // Check progress
        checkMeshgridMatch(d_pmrdr, host_dpmrdr, qgrid_tol, "charge mapping",
                           tsm.getTestingStatus());

        // GPU process: convolution
        applyConvolution(&d_cvolw, *borders_ptr, d_pmrdr, gpu, &scw);

        // CPU process: convolution
        applyConvolution(&host_dcvolw, *host_borders_ptr, host_dpmrdr, &host_scw);

        // Check progress
        checkMeshgridMatch(d_pmrdr, host_dpmrdr, ugrid_tol, "convolution applied",
                           tsm.getTestingStatus());

        // GPU process: gather forces due to the mesh-based potential
        launchGenForceGatheringKernel(dv_cgw_ptr, d_pmrdr, ct_tmat, *borders_ptr, dpoly_nbk,
                                      fintrp_lp);

        // CPU process: gather forces due to the mesh-based potential
        gatherForces<double, llint, double, double2, double4>(host_dcgw_ptr, host_dpmrdr,
                                                              *host_borders_ptr, host_dpoly_nbk);
        
        // Check progress
        checkCellGridForces<double, llint, double, double4_16a>(dd_cgr_ptr, host_dcgr_ptr, poly_ps,
                                                                ngbr_frc_tol, "mesh-based forces "
                                                                "interpolated back to particles",
                                                                tsm.getTestingStatus());
      }
      break;
    case PrecisionModel::SINGLE:
      {
        // GPU process: particle-particle pairwise interactions
        if (has_tiny_box) {
          launchPMEPairs(fpoly_nbk, lemr, fppi_direct, *borders_ptr, ff_cgw_ptr, tile_plan_ptr,
                         &scw, f_mmctrl_ptr, EvaluateForce::YES, eval_nrg, pair_lp, 0.0, 0.0);
        }
        else {
          launchPMEPairs(fpoly_nbk, lemr, fppi_direct, ff_cgw_ptr, tile_plan_ptr, &scw,
                         f_mmctrl_ptr, EvaluateForce::YES, eval_nrg, pair_lp, 0.0, 0.0);
        }

        // CPU process: particle-particle pairwise interactions
        evaluateParticleParticleEnergy<float, int,
                                       float, float2, float4>(host_fv_cgw_ptr, &host_scw,
                                                              host_poly_psr, host_fpoly_nbk,
                                                              host_lemr,
                                                              mmctrl_fe.getLongestCutoff(),
                                                              pmecon.getEwaldCoefficient(),
                                                              dyncon.getVdwSummation(),
                                                              EvaluateForce::YES,
                                                              NonbondedTheme::ALL);
        
        // Check progress
        checkCellGridForces<float, int, float, float4>(ff_cgr_ptr, host_fcgr_ptr, poly_ps,
                                                       ngbr_frc_tol, "particle-particle forces "
                                                       "evaluated", tsm.getTestingStatus());
        
        // GPU process: charge mapping for particle-mesh interactions
        switch (pmecon.getDensityMappingMethod()) {
        case QMapMethod::ACC_SHARED:
        case QMapMethod::AUTOMATIC:
          launchShrAccDensityKernel(&f_pmwrt, use_overflow, f_mmctrl_ptr, *fv_cgr_ptr,
                                    float_type_index, fpoly_nbk, qmap_lp);
          break;
        case QMapMethod::GENERAL_PURPOSE:
          launchPMIGridInitialization(&f_pmacc, gpu);
          launchGenPrpDensityKernel(&f_pmacc, *fv_cgr_ptr, float_type_index, fpoly_nbk, qmap_lp);
          launchPMIGridRealConversion(&f_pmwrt, f_pmacc, gpu);
          break;
        }
        
        // CPU process: charge mapping for particle-mesh interactions
        f_pmig.initialize();
        mapDensity<float, int,
                   float, float2, float4>(&host_fpmacc, &host_fpmwrt, *host_fcgr_ptr,
                                          host_fpoly_nbk);

        // Check progress
        checkMeshgridMatch(f_pmrdr, host_fpmrdr, qgrid_tol, "charge mapping",
                           tsm.getTestingStatus());

        // GPU process: convolution
        applyConvolution(&f_cvolw, *borders_ptr, f_pmrdr, gpu, &scw);

        // CPU process: convolution
        applyConvolution(&host_fcvolw, *host_borders_ptr, host_fpmrdr, &host_scw);

        // Check progress
        checkMeshgridMatch(f_pmrdr, host_fpmrdr, ugrid_tol, "convolution applied",
                           tsm.getTestingStatus());

        // GPU process: gather forces due to the mesh-based potential
        launchGenForceGatheringKernel(fv_cgw_ptr, f_pmrdr, ct_tmat, *borders_ptr, fpoly_nbk,
                                      fintrp_lp);

        // CPU process: gather forces due to the mesh-based potential
        gatherForces<float, int, float, float2, float4>(host_fcgw_ptr, host_fpmrdr,
                                                        *host_borders_ptr, host_fpoly_nbk);
        
        // Check progress
        checkCellGridForces<float, int, float, float4>(ff_cgr_ptr, host_fcgr_ptr, poly_ps,
                                                       ngbr_frc_tol, "mesh-based forces "
                                                       "interpolated back to particles",
                                                       tsm.getTestingStatus());
      }
      break;
    }

    // Compute the valence interactions.
    switch (val_prec) {
    case PrecisionModel::DOUBLE:
      {
        // GPU process: valence interactions
        launchValence(dpoly_vk, dpoly_rk, *dd_cgr_ptr, d_mmctrl_ptr, crdw_ptr, dpoly_auk, &d_tstw,
                      &scw, d_vale_res_ptr, EvaluateForce::YES, eval_nrg,
                      VwuGoal::ACCUMULATE, vale_lp, 0.0, 0.0);
        
        // CPU process: valence interactions, with contributions from the neighbor list.  A
        // GPU device-oriented abstract is used to supply the time step, but this constant is
        // valid on the CPU as well as the GPU, and one such particular abstract is used to
        // keep this parameter throughout this exercise.
        contributeCellGridForces<double, llint, double, double4_16a>(host_crdw_ptr,
                                                                     *host_dcgr_ptr);
        evalValeRestMM<double, double2, double4_16a>(host_crdw_ptr, &sc, host_dpoly_vk,
                                                     host_dpoly_rk, host_dpoly_auk,
                                                     EvaluateForce::YES, VwuTask::ALL_TASKS,
                                                     d_tstw.step);
        
        // Check progress
        checkSynthesisMatch(poly_ps, d_cg, gpu, total_frc_tol, TrajectoryKind::FORCES, false,
                            tsm.getTestingStatus(), val_prec, nb_prec, "complete forces "
                            "calculated");

        // Transmit forces from virtual sites to particles with mass, if needed.
        if (has_virtual_sites) {

          // CPU process
          transmitVirtualSiteForces<double, double2, double4_16a>(host_crdw_ptr, host_dpoly_vk,
                                                                  host_dpoly_auk);
        }
        
        // GPU process: advance the velocities in the first velocity-Verlet update.  This will
        // combine the forces in the cell grid with those in the coordinate synthesis, which the
        // checking routine above did by downloading but sets of forces and combining them in
        // temporary arrays on the CPU host. 
        launchIntegrationProcess(crdw_ptr, d_vale_res_ptr, d_mmctrl_ptr, scw_ptr, *dd_cgr_ptr,
                                 dpoly_vk, dpoly_auk, d_tstw, intg_vvi_lp,
                                 IntegrationStage::VELOCITY_ADVANCE);
        
        // CPU process: advance the velocities in the first velocity-Verlet update.
        velocityVerletVelocityUpdate<double, double2, double4_16a>(host_crdw_ptr, host_dpoly_auk,
                                                                   host_dtstr);

        // Check progress
        checkSynthesisMatch(poly_ps, d_cg, gpu, vel_tol, TrajectoryKind::VELOCITIES, true,
                            tsm.getTestingStatus(), val_prec, nb_prec, "initial Velocity-Verlet "
                            "velocity update");

        // Constrain velocities by RATTLE, if requested in the dynamics namelist control block
        if (d_tstw.cnst_geom) {
          
          // GPU process
          launchIntegrationProcess(crdw_ptr, d_vale_res_ptr, d_mmctrl_ptr, scw_ptr, *dd_cgr_ptr,
                                   dpoly_vk, dpoly_auk, d_tstw, intg_vvi_lp,
                                   IntegrationStage::VELOCITY_CONSTRAINT);

          // CPU process
          rattleVelocities<double, double2, double4_16a>(host_crdw_ptr, host_dpoly_vk,
                                                         host_dpoly_auk, d_tstr.dt,
                                                         d_tstr.rattle_tol, d_tstr.rattle_iter);
          settleVelocities<double, double2, double4_16a>(host_crdw_ptr, host_dpoly_vk,
                                                         host_dpoly_auk);

          // Check progress
          checkSynthesisMatch(poly_ps, d_cg, gpu, vel_tol, TrajectoryKind::VELOCITIES, true,
                              tsm.getTestingStatus(), val_prec, nb_prec, "velocity constraints "
                              "applied");
        }

        // GPU process: compute the kinetic energy
        if (on_energy_step) {
          launchIntegrationProcess(crdw_ptr, d_vale_res_ptr, d_mmctrl_ptr, scw_ptr, *dd_cgr_ptr,
                                   dpoly_vk, dpoly_auk, d_tstw, intg_ke_lp,
                                   IntegrationStage::CALC_KINETIC);
        }

        // GPU process: advance the velocities and then positions in the second velocity-verlet
        // update
        launchIntegrationProcess(crdw_ptr, d_vale_res_ptr, d_mmctrl_ptr, scw_ptr, *dd_cgr_ptr,
                                 dpoly_vk, dpoly_auk, d_tstw, intg_vvii_lp,
                                 IntegrationStage::POSITION_ADVANCE);

        // CPU process: advance the velocities, then the positions
        velocityVerletCoordinateUpdate<double, double2, double4_16a>(host_crdw_ptr, host_dpoly_auk,
                                                                     host_dtstr);

        // Check progress
        checkSynthesisMatch(poly_ps, d_cg, gpu, vel_tol, TrajectoryKind::VELOCITIES, true,
                            tsm.getTestingStatus(), val_prec, nb_prec, "Velocity-Verlet position "
                            "update");
        if (has_virtual_sites == false) {
          checkSynthesisMatch(poly_ps, d_cg, gpu, pos_tol, TrajectoryKind::POSITIONS, true,
                              tsm.getTestingStatus(), val_prec, nb_prec, "Velocity-Verlet "
                              "position update");
        }
        
        // Constrain geometries by SHAKE, if requested in the dynamics namelist control block
        if (d_tstw.cnst_geom) {

          // GPU process
          launchIntegrationProcess(crdw_ptr, d_vale_res_ptr, d_mmctrl_ptr, scw_ptr, *dd_cgr_ptr,
                                   dpoly_vk, dpoly_auk, d_tstw, intg_gc_lp,
                                   IntegrationStage::GEOMETRY_CONSTRAINT);

          // CPU process
          shakePositions(host_crdw_ptr, host_dpoly_vk, host_dpoly_auk, d_tstw.dt,
                         dyncon.getRattleTolerance(), dyncon.getRattleIterations());
          settlePositions<double, double2, double3, double4_16a>(host_crdw_ptr, host_dpoly_vk,
                                                                 host_dpoly_auk, d_tstw.dt);

          // Check progress
          checkSynthesisMatch(poly_ps, d_cg, gpu, vel_tol, TrajectoryKind::VELOCITIES, true,
                              tsm.getTestingStatus(), val_prec, nb_prec, "position constraints "
                              "applied");
          if (has_virtual_sites == false) {
            checkSynthesisMatch(poly_ps, d_cg, gpu, pos_tol, TrajectoryKind::POSITIONS, true,
                                tsm.getTestingStatus(), val_prec, nb_prec, "position constraints "
                                "applied");
          }
        }

        // Place virtual sites based on the updated coordinates, if needed
        if (has_virtual_sites) {

          // CPU process
          placeVirtualSites<double, double2, double4_16a>(host_crdw_ptr, host_dpoly_vk,
                                                          host_dpoly_auk);

          // Check progress
          checkSynthesisMatch(poly_ps, d_cg, gpu, pos_tol, TrajectoryKind::POSITIONS, true,
                              tsm.getTestingStatus(), val_prec, nb_prec, "virtual sites "
                              "repositioned");
        }
      }
      break;
    case PrecisionModel::SINGLE:
      {
        // GPU process: valence interactions, with contributions from the neighbor list
        launchValence(fpoly_vk, fpoly_rk, *ff_cgr_ptr, f_mmctrl_ptr, crdw_ptr, fpoly_auk, &f_tstw,
                      &scw, f_vale_res_ptr, EvaluateForce::YES, eval_nrg, VwuGoal::ACCUMULATE,
                      AccumulationMethod::SPLIT, vale_lp, 0.0, 0.0);

        // CPU process: valence interactions, with contributions from the neighbor list.  A
        // GPU device-oriented abstract is used to supply the time step, but this constant is
        // valid on the CPU as well as the GPU, and one such particular abstract is used to
        // keep this parameter throughout this exercise.
        contributeCellGridForces<float, int, float, float4>(host_crdw_ptr, *host_fcgr_ptr);
        evalValeRestMM<float, float2, float4>(host_crdw_ptr, &sc, host_fpoly_vk, host_fpoly_rk,
                                              host_fpoly_auk, EvaluateForce::YES,
                                              VwuTask::ALL_TASKS, d_tstw.step);
        
        // Check progress
        checkSynthesisMatch(poly_ps, f_cg, gpu, total_frc_tol, TrajectoryKind::FORCES, false,
                            tsm.getTestingStatus(), val_prec, nb_prec, "complete forces "
                            "calculated");

        // Transmit forces from virtual sites to particles with mass, if needed, on the CPU.  The
        // following GPU kernel, launched by launchIntegrationProcess() below in the
        // Velocity-Verlet initial velocity update, will do the transmission on the device data.
        if (has_virtual_sites) {
          transmitVirtualSiteForces<float, float2, float4>(host_crdw_ptr, host_fpoly_vk,
                                                           host_fpoly_auk);
        }
        
        // GPU process: advance the velocities in the first velocity-Verlet update.  This will
        // combine the forces in the cell grid with those in the coordinate synthesis, which the
        // checking routine above did by downloading but sets of forces and combining them in
        // temporary arrays on the CPU host. 
        launchIntegrationProcess(crdw_ptr, f_vale_res_ptr, f_mmctrl_ptr, scw_ptr, *ff_cgr_ptr,
                                 fpoly_vk, fpoly_auk, f_tstw, intg_vvi_lp,
                                 AccumulationMethod::SPLIT, IntegrationStage::VELOCITY_ADVANCE);
        
        // CPU process: advance the velocities in the first velocity-Verlet update.
        velocityVerletVelocityUpdate<float, float2, float4>(host_crdw_ptr, host_fpoly_auk,
                                                            host_ftstr);

        // Check progress
        checkSynthesisMatch(poly_ps, f_cg, gpu, vel_tol, TrajectoryKind::VELOCITIES, true,
                            tsm.getTestingStatus(), val_prec, nb_prec, "initial Velocity-Verlet "
                            "velocity update");

        // Constrain velocities by RATTLE, if requested in the dynamics namelist control block
        if (d_tstw.cnst_geom) {
          
          // GPU process
          launchIntegrationProcess(crdw_ptr, f_vale_res_ptr, f_mmctrl_ptr, scw_ptr, *ff_cgr_ptr,
                                   fpoly_vk, fpoly_auk, f_tstw, intg_vvi_lp,
                                   AccumulationMethod::SPLIT,
                                   IntegrationStage::VELOCITY_CONSTRAINT);

          // CPU process
          rattleVelocities<float, float2, float4>(host_crdw_ptr, host_fpoly_vk, host_fpoly_auk,
                                                  f_tstr.dt, f_tstr.rattle_tol,
                                                  f_tstr.rattle_iter);
          settleVelocities<float, float2, float4>(host_crdw_ptr, host_fpoly_vk, host_fpoly_auk);

          // Check progress
          checkSynthesisMatch(poly_ps, f_cg, gpu, vel_tol, TrajectoryKind::VELOCITIES, true,
                              tsm.getTestingStatus(), val_prec, nb_prec, "velocity constraints "
                              "applied");
        }

        // GPU process: compute the kinetic energy
        if (on_energy_step) {
          launchIntegrationProcess(crdw_ptr, f_vale_res_ptr, f_mmctrl_ptr, scw_ptr, *ff_cgr_ptr,
                                   fpoly_vk, fpoly_auk, f_tstw, intg_ke_lp,
                                   AccumulationMethod::SPLIT, IntegrationStage::CALC_KINETIC);
        }

        // GPU process: advance the velocities and then positions in the second velocity-verlet
        // update
        launchIntegrationProcess(crdw_ptr, f_vale_res_ptr, f_mmctrl_ptr, scw_ptr, *ff_cgr_ptr,
                                 fpoly_vk, fpoly_auk, f_tstw, intg_vvii_lp,
                                 AccumulationMethod::SPLIT, IntegrationStage::POSITION_ADVANCE);

        // CPU process: advance the velocities, then the positions
        velocityVerletCoordinateUpdate<float, float2, float4>(host_crdw_ptr, host_fpoly_auk,
                                                              host_ftstr);

        // Check progress
        checkSynthesisMatch(poly_ps, f_cg, gpu, vel_tol, TrajectoryKind::VELOCITIES, true,
                            tsm.getTestingStatus(), val_prec, nb_prec, "Velocity-Verlet position "
                            "update");
        if (has_virtual_sites == false) {
          checkSynthesisMatch(poly_ps, f_cg, gpu, pos_tol, TrajectoryKind::POSITIONS, true,
                              tsm.getTestingStatus(), val_prec, nb_prec, "Velocity-Verlet "
                              "position update");
        }

        // Constrain geometries by SHAKE, if requested in the dynamics namelist control block
        if (f_tstw.cnst_geom) {

          // GPU process
          launchIntegrationProcess(crdw_ptr, f_vale_res_ptr, f_mmctrl_ptr, scw_ptr, *ff_cgr_ptr,
                                   fpoly_vk, fpoly_auk, f_tstw, intg_gc_lp,
                                   AccumulationMethod::SPLIT,
                                   IntegrationStage::GEOMETRY_CONSTRAINT);

          // CPU process
          shakePositions<float, float2, float4>(host_crdw_ptr, host_fpoly_vk, host_fpoly_auk,
                                                f_tstw.dt, dyncon.getRattleTolerance(),
                                                dyncon.getRattleIterations());
          settlePositions<float, float2, float3, float4>(host_crdw_ptr, host_fpoly_vk,
                                                         host_fpoly_auk, f_tstw.dt);

          // Check progress
          checkSynthesisMatch(poly_ps, f_cg, gpu, vel_tol, TrajectoryKind::VELOCITIES, true,
                              tsm.getTestingStatus(), val_prec, nb_prec, "position constraints "
                              "applied");
          if (has_virtual_sites == false) {
            checkSynthesisMatch(poly_ps, f_cg, gpu, pos_tol, TrajectoryKind::POSITIONS, true,
                                tsm.getTestingStatus(), val_prec, nb_prec, "position constraints "
                                "applied");
          }
        }

        // Place virutal sites based on the updated coordinates
        if (has_virtual_sites) {

          // CPU process
          placeVirtualSites<float, float2, float4>(host_crdw_ptr, host_fpoly_vk, host_fpoly_auk);

          // Check progress
          checkSynthesisMatch(poly_ps, f_cg, gpu, pos_tol, TrajectoryKind::POSITIONS, true,
                              tsm.getTestingStatus(), val_prec, nb_prec, "virtual sites "
                              "repositioned");
        }
      }
      break;
    }

    // Migrate particles within the neighbor list
    switch (nb_prec) {
    case PrecisionModel::DOUBLE:
      
      // GPU process: migrate particles within the neighbor list
      launchMigration(dd_cgw_ptr, *dcorg_ptr, *crdw_ptr, migr_one_lp, migr_two_lp);

      // CPU process: migrate particles within the neighbor list
      migrate(host_dcgw_ptr, *host_crdr_ptr);
      break;
    case PrecisionModel::SINGLE:
      
      // GPU process: migrate particles within the neighbor list
      launchMigration(ff_cgw_ptr, *fcorg_ptr, *crdw_ptr, migr_one_lp, migr_two_lp);

      // CPU process: migrate particles within the neighbor list
      migrate(host_fcgw_ptr, *host_crdr_ptr);
      break;
    }

    // Advance the two cell grids (only one is in use) and coordinate synthesis.
    d_cg.updateCyclePosition();
    f_cg.updateCyclePosition();
    poly_ps.updateCyclePosition();

    // Advance step counters.
    d_tstw.step += 1;
    f_tstw.step += 1;
    d_ctrl_fx.step += 1;
    d_ctrl_fe.step += 1;
    f_ctrl_fx.step += 1;
    f_ctrl_fe.step += 1;
    tst.incrementStep();

    // Check progress in the cell grid.  The cell grid abstract must be reset based on the updated
    // point in the coordinate time cycle, so that the updated positions may be checked rather than
    // the ones that were used for force computations earlier in the cycle.
    switch (nb_prec) {
    case PrecisionModel::DOUBLE:
      if (step_idx & 0x1) {
        dd_cgr_ptr = &dd_cgr_alt;
        host_dcgr_ptr = &host_dcgr_alt;
      }
      else {
        dd_cgr_ptr = &dd_cgr;
        host_dcgr_ptr = &host_dcgr;
      } 
      checkCellGridPositions<double, llint, double, double4_16a>(dd_cgr_ptr, host_dcgr_ptr,
                                                                 poly_ps, pos_tol,
                                                                 tsm.getTestingStatus());
      break;
    case PrecisionModel::SINGLE:
      if (step_idx & 0x1) {
        ff_cgr_ptr = &ff_cgr_alt;
        host_fcgr_ptr = &host_fcgr_alt;
      }
      else {
        ff_cgr_ptr = &ff_cgr;
        host_fcgr_ptr = &host_fcgr;
      }
      checkCellGridPositions<float, int, float, float4>(ff_cgr_ptr, host_fcgr_ptr, poly_ps,
                                                        pos_tol, tsm.getTestingStatus());
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Carry out one step for dynamics in a particular system using the GPU, then replicate the process
// on the CPU.  Initial velocities are set to zero or initialized according to a Maxwell
// distribution on the CPU, then ported to the GPU.
//
// Arguments:
//   tsm:         Collection of available test systems
//   test_index:  Index of the system to pull, replicate, and test
//   gpu:         Specifications of the GPU to use in calculations
//   prec:        The precision in which to perform GPU calculations
//   gb_model:    The type of implicit solvent model to apply to the topology
//   use_rattle:  Whether to apply geometry constraints to the system
//   nrep:        The number of replicas to make in the GPU-capable synthesis
//   nstep:       The number of steps over which to test for CPU and GPU tracking
//   pert_sigma:  Sigma width of the Gaussian perturbation to apply to initial positions (ensures a
//                nonzero force on particles)
//   init_temp:   The initial temperature at which to start the system
//   time_step:   The time step to apply
//   rng_seed:    Random number seed for the Andersen velocity kick-start
//-------------------------------------------------------------------------------------------------
void singleStepLaboratory(const TestSystemManager &tsm, const int test_index,
                          const GpuDetails &gpu,
                          const PrecisionModel prec = PrecisionModel::DOUBLE,
                          const ImplicitSolventModel gb_model = ImplicitSolventModel::NONE,
                          const ApplyConstraints use_rattle = ApplyConstraints::NO,
                          const int nrep = 1, const int nstep = 1, const double pert_sigma = 0.0,
                          const double init_temp = 0.0, const double time_step = 1.0,
                          const int rng_seed = 601847294) {

  // Create a random number generator with a unique seed from the thermostat.  This generator will
  // manage perturbations.
  Xoroshiro128pGenerator xrs(rng_seed + 5);

  // Obtain the correct topology.  Set the implicit solvent model.
  AtomGraph ag = tsm.exportAtomGraph(test_index);
  AtomicRadiusSet rads;
  switch (gb_model) {
  case ImplicitSolventModel::NONE:
    rads = AtomicRadiusSet::NONE;
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
    rads = AtomicRadiusSet::MBONDI2;
    break;
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    rads = AtomicRadiusSet::MBONDI3;
    break;
  }
  ag.setImplicitSolventModel(gb_model, 80.0, 0.0, rads, ExceptionResponse::WARN);
  
  // Create other components for propagating the dynamics.
  StaticExclusionMask se(ag);
  RestraintApparatus ra(&ag);
  Thermostat tst(ag, ThermostatKind::NONE, init_temp);
  tst.setGeometryConstraints(use_rattle);
  NeckGeneralizedBornTable ngb_tab;
  
  // The dynamics controls guide operations in the CPU routine
  DynamicsControls dyncon;
  dyncon.setTimeStep(time_step);
  dyncon.setStepCount(1);
  dyncon.setDiagnosticPrintFrequency(1);
  dyncon.setGeometricConstraints(use_rattle);

  // Create replicas of the coordinate system
  std::vector<PhaseSpace> ps_cpu;
  PhaseSpace ps = tsm.exportPhaseSpace(test_index);
  MinimizeControls mincon;
  mincon.setTotalCycles(50);
  mincon.setClashDampingCycles(25);
  mincon.setDiagnosticPrintFrequency(5);
  ScoreCard scmin = minimize(&ps, ag, ngb_tab, ra, se, mincon);
  for (int i = 0; i < nrep; i++) {
    ps_cpu.push_back(ps);
    PhaseSpaceWriter psw = ps_cpu.back().data();
    addRandomNoise(&xrs, psw.xcrd, psw.ycrd, psw.zcrd, psw.natom, pert_sigma);
    velocityKickStart(&ps_cpu[i], ag, &tst, dyncon, EnforceExactTemperature::YES);
  }
  std::vector<PhaseSpace> psv_ref = ps_cpu;

  // Create components for GPU dynamics of the same replicated system
  const std::vector<AtomGraph*> agv(1, &ag);
  const std::vector<StaticExclusionMask*> sev(1, &se);
  int gpos_bits, vel_bits, frc_bits;
  const int lpos_bits = 26;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    gpos_bits = 48;
    vel_bits  = 48;
    frc_bits  = 48;
    break;
  case PrecisionModel::SINGLE:
    gpos_bits = 32;
    vel_bits  = 36;
    frc_bits  = 40;    
    break;
  }
  PhaseSpaceSynthesis poly_ps(psv_ref, incrementingSeries(0, nrep), agv, std::vector<int>(nrep, 0),
                              gpos_bits, lpos_bits, vel_bits, frc_bits);
  AtomGraphSynthesis poly_ag(agv, std::vector<int>(nrep, 0));
  StaticExclusionMaskSynthesis poly_se(sev, std::vector<int>(nrep, 0));
  InitializationTask init_order;
  switch (gb_model) {
  case ImplicitSolventModel::NONE:
    init_order = InitializationTask::GENERAL_DYNAMICS;
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    init_order = InitializationTask::GB_DYNAMICS;
    break;
  }
  Thermostat poly_tst(poly_ag, ThermostatKind::NONE, init_temp);
  poly_tst.setGeometryConstraints(use_rattle);
  poly_tst.setRandomCacheDepth(1);
  poly_tst.initializeRandomStates(rng_seed, 25,  HybridTargetLevel::DEVICE, gpu);
  poly_ag.loadNonbondedWorkUnits(poly_se, init_order, poly_tst.getRandomCacheDepth(), gpu);
  
  // Upload critical data to the GPU  
  poly_ps.upload();
  poly_ag.upload();
  poly_se.upload();  
  const CoreKlManager launcher(gpu, poly_ag);  
  MolecularMechanicsControls mmctrl;
  mmctrl.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                               VwuGoal::MOVE_PARTICLES, prec, prec, poly_ag);
  const int2 vale_lp = launcher.getValenceKernelDims(prec, EvaluateForce::YES, EvaluateEnergy::YES,
                                                     AccumulationMethod::SPLIT,
                                                     VwuGoal::MOVE_PARTICLES, ClashResponse::NONE);
  const int2 nonb_lp = launcher.getNonbondedKernelDims(prec, poly_ag.getNonbondedWorkType(),
                                                       EvaluateForce::YES, EvaluateEnergy::YES,
                                                       AccumulationMethod::SPLIT, gb_model,
                                                       ClashResponse::NONE);
  CacheResource vale_tb_space(vale_lp.x, maximum_valence_work_unit_atoms);
  CacheResource nonb_tb_space(nonb_lp.x, small_block_max_atoms);
  ImplicitSolventWorkspace ism_space(poly_ag.getSystemAtomOffsets(), poly_ag.getSystemAtomCounts(),
                                     prec);
  
  // Run dynamics one step at a time for the requested number of cycles.
  for (int step_idx = 0; step_idx < nstep; step_idx++) {
    
    // Run CPU dynamics, one step.
    ScoreCard scmd_cpu(nrep, 1, 32);
    for (int i = 0; i < nrep; i++) {
      dynamics(&ps_cpu[i], &tst, &scmd_cpu, ag, ngb_tab, se, ra, dyncon, i);
    }
    
    // Compute forces on the GPU.  Move particles.
    ScoreCard scmd_gpu(nrep, 1, 32);
    launchNonbonded(prec, poly_ag, poly_se, &mmctrl, &poly_ps, &poly_tst, &scmd_gpu,
                    &nonb_tb_space, &ism_space, EvaluateForce::YES, EvaluateEnergy::YES, launcher);
    launchValence(prec, poly_ag, &mmctrl, &poly_ps, &poly_tst, &scmd_gpu, &vale_tb_space,
                  EvaluateForce::YES, EvaluateEnergy::YES, VwuGoal::MOVE_PARTICLES, launcher);
    mmctrl.incrementStep();
    poly_ps.updateCyclePosition();
    
    // Check that the original states are identical, to within the precision of the fixed-point
    // representation.
    poly_ps.download();
    std::vector<PhaseSpace> ps_gpu;
    for (int i = 0; i < nrep; i++) {
      ps_gpu.emplace_back(poly_ps.exportSystem(i));
    }
    const TrajectoryKind tpos = TrajectoryKind::POSITIONS;
    const TrajectoryKind tvel = TrajectoryKind::VELOCITIES;
    for (int i = 0; i < nrep; i++) {

      // Analyze the first three steps and the final step.
      if (step_idx < 3 || step_idx == nstep - 1) {

        // Check the current coordinates (positions and velocities)
        const std::vector<double> cpu_pos_curr = ps_cpu[i].getInterlacedCoordinates(tpos);
        const std::vector<double> gpu_pos_curr = ps_gpu[i].getInterlacedCoordinates(tpos);
        check(gpu_pos_curr, RelationalOperator::EQUAL, cpu_pos_curr, "The GPU and CPU do not "
              "indicate the same particle positions in replica " + std::to_string(i) + " after " +
              std::to_string(step_idx + 1) + " steps.  Precision model: " +
              getEnumerationName(prec) + ", GB model: " + getEnumerationName(gb_model) +
              ", constraints: " + getEnumerationName(use_rattle) + ".", tsm.getTestingStatus());
        const std::vector<double> cpu_vel_curr = ps_cpu[i].getInterlacedCoordinates(tvel);
        const std::vector<double> gpu_vel_curr = ps_gpu[i].getInterlacedCoordinates(tvel);
        check(gpu_vel_curr, RelationalOperator::EQUAL, cpu_vel_curr, "The GPU and CPU do not "
              "indicate the same velocities in replica " + std::to_string(i) + " after " +
              std::to_string(step_idx + 1) + " steps.  Precision model: " +
              getEnumerationName(prec) + ", GB model: " + getEnumerationName(gb_model) +
              ", constraints: " + getEnumerationName(use_rattle) + ".", tsm.getTestingStatus());

        // Check the previous coordinates
        ps_cpu[i].updateCyclePosition();
        ps_gpu[i].updateCyclePosition();
        const std::vector<double> cpu_pos_orig = ps_cpu[i].getInterlacedCoordinates(tpos);
        const std::vector<double> gpu_pos_orig = ps_gpu[i].getInterlacedCoordinates(tpos);
        check(gpu_pos_orig, RelationalOperator::EQUAL, cpu_pos_orig, "The GPU and CPU do not "
              "indicate the same reference particle positions in replica " + std::to_string(i) +
              " after " + std::to_string(step_idx + 1) + " steps.  Precision model: " +
              getEnumerationName(prec) + ", GB model: " + getEnumerationName(gb_model) +
              ", constraints: " + getEnumerationName(use_rattle) + ".", tsm.getTestingStatus());
        const std::vector<double> cpu_vel_orig = ps_cpu[i].getInterlacedCoordinates(tvel);
        const std::vector<double> gpu_vel_orig = ps_gpu[i].getInterlacedCoordinates(tvel);
        check(gpu_vel_orig, RelationalOperator::EQUAL, cpu_vel_orig, "The GPU and CPU do not "
              "indicate the same reference velocities in replica " + std::to_string(i) +
              " after " + std::to_string(step_idx + 1) + " steps.  Precision model: " +
              getEnumerationName(prec) + ", GB model: " + getEnumerationName(gb_model) +
              ", constraints: " + getEnumerationName(use_rattle) + ".", tsm.getTestingStatus());
        ps_cpu[i].updateCyclePosition();
        ps_gpu[i].updateCyclePosition();
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Test a rigid water system to ensure that SETTLE constraints are operating in a way that
// conserves energy.
//
// Arguments:
//   base_crd_name:  Root path to testing input coordinates (${STORMM_SOURCE}/test/Trajectory)
//   base_top_name:  Root path to testing topologies (${STORMM_SOURCE}/test/Topology)
//   prec:  The precision in which to perform calculations
//-------------------------------------------------------------------------------------------------
void checkRigidWater(const std::string &base_crd_name, const std::string &base_top_name,
                     const PrecisionModel prec) {

  // Construct the simple TIP3P system.
  const std::vector<std::string> waters(1, std::string("tip3p"));
  TestSystemManager tsm(base_top_name, "top", waters, base_crd_name, "inpcrd", waters);
  PhaseSpace tip3p_ps = tsm.exportPhaseSpace(0);
  AtomGraph tip3p_ag = tsm.exportAtomGraph(0);
  
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
    }
    break;
  case PrecisionModel::SINGLE:
    {
      
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
// Test each stage of the integration using the individual kernels.  This process involves many
// more memory "round trips" than a single, fused kernel, but permits greater flexibility to reach
// in at each stage and make modifications or collect metrics like the kinetic energy.
//
// Arguments:
//   tsm:       A collection of test systems, also offering an indication of whether testing is
//              possible
//   bcs:       Boundary conditions of systems to pluck from the test collection
//   mincon:    Energy minimization control input, for structural relaxation to ensure smooth
//              dynamics
//   dyncon:    Molecular dynamics control input
//   watcon:    Control input containing solvent details
//   prec:      The precision in which to calculate particle movement and constraints
//   gpu:       Details of the GPU that will perform the calculations
//-------------------------------------------------------------------------------------------------
void piecewiseIntegrationTests(const TestSystemManager &tsm, const std::vector<UnitCellType> &bcs,
                               const MinimizeControls &mincon, const DynamicsControls &dyncon,
                               const SolventControls &watcon, const PrecisionModel prec,
                               const GpuDetails &gpu) {

  // Find all systems in isolated boundary conditions, then in periodic boundary conditions.
  PhaseSpaceSynthesis poly_ps_direct = tsm.exportPhaseSpaceSynthesis(bcs);
  poly_ps_direct.upload();
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(bcs);
  const int nsys = poly_ag.getSystemCount();
  const std::vector<int> all_members = incrementingSeries<int>(0, poly_ag.getSystemCount());
  ScoreCard sc_direct(poly_ps_direct.getSystemCount(), 1, 36);
  ScoreCard sc_pieces(poly_ps_direct.getSystemCount(), 1, 36);
  Xoshiro256ppGenerator xrs(39945273);

  // Create general resources
  Thermostat tst_direct(poly_ag, ThermostatKind::NONE, 100.0);
  tst_direct.setTimeStep(dyncon.getTimeStep());
  tst_direct.setGeometryConstraints(dyncon.constrainGeometry());
  tst_direct.setRattleTolerance(dyncon.getRattleTolerance());
  tst_direct.setRattleIterations(dyncon.getRattleIterations());
  Thermostat tst_pieces = tst_direct;
  MolecularMechanicsControls mmctrl_direct;
  MolecularMechanicsControls mmctrl_pieces;
  std::vector<IntegrationStage> intg_stages(1, IntegrationStage::VELOCITY_ADVANCE);
  switch (dyncon.constrainGeometry()) {
  case ApplyConstraints::YES:
    intg_stages.push_back(IntegrationStage::VELOCITY_CONSTRAINT);
    intg_stages.push_back(IntegrationStage::POSITION_ADVANCE);
    intg_stages.push_back(IntegrationStage::GEOMETRY_CONSTRAINT);
    break;
  case ApplyConstraints::NO:
    intg_stages.push_back(IntegrationStage::POSITION_ADVANCE);
    break;
  }

  // Lay out arrays of PhaseSpace objects to record the positions, velocities, and forces of each
  // particle at each step of the dynamics.
  std::vector<PhaseSpace> calc_forces, velocity_adv, velocity_cnst;
  std::vector<PhaseSpace> calc_kinetic, position_adv, geom_cnst, full_step;
  const int nstep = 10;
  calc_forces.reserve(nstep);
  velocity_adv.reserve(nstep);
  velocity_cnst.reserve(nstep);
  calc_kinetic.reserve(nstep);
  position_adv.reserve(nstep);
  geom_cnst.reserve(nstep);
  full_step.reserve(nstep);

  // Common descriptions to help developers understand the meaning of test failures
  const std::string gen_desc(" recorded by a piecewise dynamics workflow do not match those in "
                             "which the valence force calculation is fused to integration of the "
                             "equations of motion.  The comparison involves root-mean squared "
                             "deviations of each system's atomic positions at each step, listing "
                             "all systems 0-" + std::to_string(poly_ps_direct.getSystemCount()) +
                             " for step 0, then for step 1, step 2, ...");
  const std::string curr_stage_desc(" (in the current stage of the coordinate time cycle)");
  const std::string alt_stage_desc(" (in the alternate stage of the coordinate time cycle)");

  // Bifurcate the integration along boundary condition types, implicit versus explicit solvent
  // dynamics.  Each uses a different set of force calculations.
  switch (poly_ag.getUnitCellType()) {
  case UnitCellType::NONE:
    {
      // Set the implicit solvent model
      NeckGeneralizedBornTable ngb_tab;
      const ImplicitSolventModel gb_model = watcon.getImplicitSolventModel();
      poly_ag.setImplicitSolventModel(gb_model, ngb_tab.dpData(), watcon.getPBRadiiSet());

      // Create resources for both paths
      StaticExclusionMaskSynthesis poly_se(poly_ag.getSystemTopologyPointer(), all_members);
      poly_ag.loadNonbondedWorkUnits(poly_se, InitializationTask::NONE, 0, gpu);
      const CoreKlManager launcher(gpu, poly_ag);
      poly_ag.upload();
      poly_se.upload();

      // Perform a brief energy minimization, then clone the coordinates.
      ScoreCard emin = launchMinimization(poly_ag, poly_se, &poly_ps_direct, mincon, gpu, prec);

      DynamicsControls mod_dyncon = dyncon;
      mod_dyncon.setThermostatSeed(dyncon.getThermostatSeed() + 715829320);
      mod_dyncon.setThermostatKind("langevin");
      const double temperature_start = dyncon.getInitialTemperatureTargets()[0];
      Thermostat kickstarter(poly_ag, ThermostatKind::LANGEVIN, temperature_start,
                             temperature_start, 0, 0);
      velocityKickStart(&poly_ps_direct, poly_ag, &kickstarter, mod_dyncon, prec,
                        EnforceExactTemperature::YES);
      poly_ps_direct.upload(TrajectoryKind::VELOCITIES);
      
      // Copy the state of each system so that two branches of dynamics can start from the same
      // set of states.
      PhaseSpaceSynthesis poly_ps_pieces = poly_ps_direct;
      
      // Create specific resources and carry out the "direct" approach, wherein the complete time
      // step is fused into a single kernel.
      mmctrl_direct.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                                          VwuGoal::MOVE_PARTICLES, prec, prec, poly_ag);
      const int2 vale_direct_lp = launcher.getValenceKernelDims(prec, EvaluateForce::YES,
                                                                EvaluateEnergy::YES,
                                                                AccumulationMethod::SPLIT,
                                                                VwuGoal::MOVE_PARTICLES,
                                                                ClashResponse::NONE);
      const int2 nonb_direct_lp = launcher.getNonbondedKernelDims(prec,
                                                                  poly_ag.getNonbondedWorkType(),
                                                                  EvaluateForce::YES,
                                                                  EvaluateEnergy::YES,
                                                                  AccumulationMethod::SPLIT,
                                                                  gb_model, ClashResponse::NONE);
      CacheResource vale_dir_tb_space(vale_direct_lp.x, maximum_valence_work_unit_atoms);
      CacheResource nonb_dir_tb_space(nonb_direct_lp.x, small_block_max_atoms);
      ImplicitSolventWorkspace ism_dir_space(poly_ag.getSystemAtomOffsets(),
                                             poly_ag.getSystemAtomCounts(), prec);
      poly_ps_direct.initializeForces(gpu, HybridTargetLevel::DEVICE);
      for (int i = 0; i < nstep; i++) {
        sc_direct.initialize(HybridTargetLevel::DEVICE, gpu);
        launchNonbonded(prec, poly_ag, poly_se, &mmctrl_direct, &poly_ps_direct, &tst_direct,
                        &sc_direct, &nonb_dir_tb_space, &ism_dir_space, EvaluateForce::YES,
                        EvaluateEnergy::YES, launcher);
        launchValence(prec, poly_ag, &mmctrl_direct, &poly_ps_direct, &tst_direct, &sc_direct,
                      &vale_dir_tb_space, EvaluateForce::YES, EvaluateEnergy::YES,
                      VwuGoal::MOVE_PARTICLES, launcher);
        for (int j = 0; j < nsys; j++) {
          full_step.push_back(poly_ps_direct.exportSystem(j, HybridTargetLevel::DEVICE));
        }
        poly_ps_direct.updateCyclePosition();
        tst_direct.incrementStep();
        mmctrl_direct.incrementStep();
      }
      
      // Create specific resources and carry out the "piecewise" approach, wherein the time step
      // is broken up into individual kernels.
      mmctrl_pieces.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                                          VwuGoal::ACCUMULATE, prec, prec, poly_ag);
      const int2 vale_pieces_lp = launcher.getValenceKernelDims(prec, EvaluateForce::YES,
                                                                EvaluateEnergy::YES,
                                                                AccumulationMethod::SPLIT,
                                                                VwuGoal::ACCUMULATE,
                                                                ClashResponse::NONE);
      const int2 nonb_pieces_lp = launcher.getNonbondedKernelDims(prec,
                                                                  poly_ag.getNonbondedWorkType(),
                                                                  EvaluateForce::YES,
                                                                  EvaluateEnergy::YES,
                                                                  AccumulationMethod::SPLIT,
                                                                  gb_model, ClashResponse::NONE);
      CacheResource vale_pcs_tb_space(vale_pieces_lp.x, maximum_valence_work_unit_atoms);
      CacheResource nonb_pcs_tb_space(nonb_pieces_lp.x, small_block_max_atoms);
      ImplicitSolventWorkspace ism_pcs_space(poly_ag.getSystemAtomOffsets(),
                                             poly_ag.getSystemAtomCounts(), prec);
      for (int i = 0; i < nstep; i++) {
        poly_ps_pieces.initializeForces(gpu, HybridTargetLevel::DEVICE);
        launchNonbonded(prec, poly_ag, poly_se, &mmctrl_pieces, &poly_ps_pieces, &tst_pieces,
                        &sc_pieces, &nonb_pcs_tb_space, &ism_pcs_space, EvaluateForce::YES,
                        EvaluateEnergy::YES, launcher);
        launchValence(prec, poly_ag, &mmctrl_pieces, &poly_ps_pieces, &tst_pieces, &sc_pieces,
                      &vale_pcs_tb_space, EvaluateForce::YES, EvaluateEnergy::YES,
                      VwuGoal::ACCUMULATE, launcher);
        for (int j = 0; j < nsys; j++) {
          calc_forces.push_back(poly_ps_pieces.exportSystem(j, HybridTargetLevel::DEVICE));
        }
        for (size_t j = 0; j < intg_stages.size(); j++) {
          const int2 klaunch_params = launcher.getIntegrationKernelDims(prec,
                                                                        AccumulationMethod::SPLIT,
                                                                        intg_stages[j]);
          launchIntegrationProcess(&poly_ps_pieces, &vale_pcs_tb_space, &tst_pieces,
                                   &mmctrl_pieces, &sc_pieces, poly_ag, launcher, prec,
                                   AccumulationMethod::SPLIT, intg_stages[j]);
          switch (intg_stages[j]) {
          case IntegrationStage::CALC_FORCES:
            break;
          case IntegrationStage::VELOCITY_ADVANCE:
            for (int k = 0; k < nsys; k++) {
              velocity_adv.push_back(poly_ps_pieces.exportSystem(k, HybridTargetLevel::DEVICE));
            }
            break;
          case IntegrationStage::VELOCITY_CONSTRAINT:
            for (int k = 0; k < nsys; k++) {
              velocity_cnst.push_back(poly_ps_pieces.exportSystem(k, HybridTargetLevel::DEVICE));
            }
            break;
          case IntegrationStage::CALC_KINETIC:
            for (int k = 0; k < nsys; k++) {
              calc_kinetic.push_back(poly_ps_pieces.exportSystem(k, HybridTargetLevel::DEVICE));
            }
            break;
          case IntegrationStage::POSITION_ADVANCE:
            for (int k = 0; k < nsys; k++) {
              position_adv.push_back(poly_ps_pieces.exportSystem(k, HybridTargetLevel::DEVICE));
            }
            break;
          case IntegrationStage::GEOMETRY_CONSTRAINT:
            for (int k = 0; k < nsys; k++) {
              geom_cnst.push_back(poly_ps_pieces.exportSystem(k, HybridTargetLevel::DEVICE));
            }
            break;
          }
        }
        poly_ps_pieces.updateCyclePosition();
        tst_pieces.incrementStep();
        mmctrl_pieces.incrementStep();
      }

      // Forces will not be consistent between the two simulation lines--the piecewise simulation
      // records the complete forces but the fused kernels compute the last force contributions
      // and use them immediately, never to write them back to the coordinate synthesis.
      std::vector<double> position_deviations(nstep * nsys, 0.0);
      std::vector<double> velocity_deviations(nstep * nsys, 0.0);
      std::vector<double> pos_alt_deviations(nstep * nsys, 0.0);
      std::vector<double> vel_alt_deviations(nstep * nsys, 0.0);
      for (int i = 0; i < nstep; i++) {
        for (int j = 0; j < nsys; j++) {
          const size_t ij_idx = (i * nsys) + j;
          PhaseSpaceReader psij_direct = full_step[ij_idx].data();
          PhaseSpaceReader psij_pieces = full_step[ij_idx].data();
          for (int k = 0; k < psij_direct.natom; k++) {
            const double dx = psij_pieces.xcrd[k] - psij_direct.xcrd[k];
            const double dy = psij_pieces.ycrd[k] - psij_direct.ycrd[k];
            const double dz = psij_pieces.zcrd[k] - psij_direct.zcrd[k];
            position_deviations[ij_idx] += (dx * dx) + (dy * dy) + (dz * dz);
            const double dx_a = psij_pieces.xalt[k] - psij_direct.xalt[k];
            const double dy_a = psij_pieces.yalt[k] - psij_direct.yalt[k];
            const double dz_a = psij_pieces.zalt[k] - psij_direct.zalt[k];
            pos_alt_deviations[ij_idx] += (dx_a * dx_a) + (dy_a * dy_a) + (dz_a * dz_a);
            const double vx = psij_pieces.xvel[k] - psij_direct.xvel[k];
            const double vy = psij_pieces.yvel[k] - psij_direct.yvel[k];
            const double vz = psij_pieces.zvel[k] - psij_direct.zvel[k];
            velocity_deviations[ij_idx] += (vx * vx) + (vy * vy) + (vz * vz);
            const double vx_a = psij_pieces.vxalt[k] - psij_direct.vxalt[k];
            const double vy_a = psij_pieces.vyalt[k] - psij_direct.vyalt[k];
            const double vz_a = psij_pieces.vzalt[k] - psij_direct.vzalt[k];
            vel_alt_deviations[ij_idx] += (vx_a * vx_a) + (vy_a * vy_a) + (vz_a * vz_a);
          }
          position_deviations[ij_idx] = sqrt(position_deviations[ij_idx] / psij_direct.natom);
          pos_alt_deviations[ij_idx] = sqrt(pos_alt_deviations[ij_idx] / psij_direct.natom);
          velocity_deviations[ij_idx] = sqrt(velocity_deviations[ij_idx] / psij_direct.natom);
          vel_alt_deviations[ij_idx] = sqrt(vel_alt_deviations[ij_idx] / psij_direct.natom);
        }
      }
      check(position_deviations, RelationalOperator::EQUAL, std::vector<double>(nstep * nsys, 0.0),
            "Particle positions" + curr_stage_desc + gen_desc, tsm.getTestingStatus());
      check(pos_alt_deviations, RelationalOperator::EQUAL, std::vector<double>(nstep * nsys, 0.0),
            "Particle positions" + alt_stage_desc + gen_desc, tsm.getTestingStatus());
      check(velocity_deviations, RelationalOperator::EQUAL, std::vector<double>(nstep * nsys, 0.0),
            "Particle velocities" + curr_stage_desc + gen_desc, tsm.getTestingStatus());
      check(vel_alt_deviations, RelationalOperator::EQUAL, std::vector<double>(nstep * nsys, 0.0),
            "Particle velocities" + alt_stage_desc + gen_desc, tsm.getTestingStatus());
    }
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    {
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Initialize the test environment
  const TestEnvironment oe(argc, argv);
  StopWatch timer;

  // Prep the GPU
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  const Hybrid<int> array_to_trigger_gpu_mapping(1);

  // Section 1: test the thermostat mechanism
  section("Thermostat diagnostics");

  // Section 2: test the propagation of atoms against CPU results
  section("Compare GPU and CPU dynamics");

  // Section 3: test consistency of system propagation
  section("Self-consistency of GPU dynamics");
  
  // Read topology and starting coordinate files
  const char osc = osSeparator();
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";

  // Check the GPU PME dynamics process.
  const std::vector<std::string> pbc_mols_ph = { "bromobenzene", "tip3p", "tip4p",
                                                 "trpcage_in_water", "trpcage_in_water",
                                                 "ubiquitin", "ubiquitin", "drug_example" };
  TestSystemManager pbc_tsm_ph(base_top_name, "top", pbc_mols_ph, base_crd_name, "inpcrd",
                               pbc_mols_ph);
  DynamicsControls dyncon;
  dyncon.setDiagnosticPrintFrequency(2);
  dyncon.setStepCount(5);
  dyncon.setCutoff(7.0);
  dyncon.setGeometricConstraints("yes");
  dyncon.setRattleTolerance(1.0e-6);
  dyncon.setTimeStep(1.0);
  PrecisionControls preccon;
  preccon.setGlobalPosScalingBits(38);
  preccon.setVelocityScalingBits(48);
  preccon.setForceScalingBits(40);
  const std::vector<PrecisionModel> all_prec = { PrecisionModel::DOUBLE, PrecisionModel::SINGLE };
  RandomControls rngcon;
  PPPMControls pmecon;
  const std::vector<double> all_cutoffs = { 7.0, 8.8 };
  for (size_t i = 0; i < all_prec.size(); i++) {
    for (size_t j = 0; j < all_prec.size(); j++) {

      // Different cutoffs are tested to force some of the simulations to produce neighbor lists
      // with four cells along each axis.  This activates the "small box" contingency.
      for (size_t k = 0; k < all_cutoffs.size(); k++) {
        pmecon.setCutoff(all_cutoffs[k]);
        preccon.setNonbondedMethod(all_prec[i]);
        preccon.setValenceMethod(all_prec[j]);
      
        // For each combination of precision settings, there are appropriate tolerances for force,
        // velocity, and position calculations.
        double qgrid_tol, ugrid_tol, ngbr_frc_tol, total_frc_tol, vel_tol, pos_tol;
        switch (preccon.getNonbondedMethod()) {
        case PrecisionModel::DOUBLE:
          switch (preccon.getValenceMethod()) {
          case PrecisionModel::DOUBLE:
            ngbr_frc_tol = 1.6e-5;
            qgrid_tol = 1.4e-7;
            ugrid_tol = 5.0e-5;
            total_frc_tol = 6.0e-4;
            vel_tol = 6.8e-7;
            pos_tol = 6.8e-7;
            break;
          case PrecisionModel::SINGLE:
            ngbr_frc_tol = 1.3e-3;
            qgrid_tol = 3.7e-7;
            ugrid_tol = 7.8e-5;
            total_frc_tol = 1.9e-3;
            vel_tol = 6.0e-6;
            pos_tol = 5.8e-6;
            break;
          }
          break;
        case PrecisionModel::SINGLE:
          switch (preccon.getValenceMethod()) {
          case PrecisionModel::DOUBLE:
            ngbr_frc_tol = 5.5e-3;
            qgrid_tol = 6.0e-7;
            ugrid_tol = 4.0e-4;
            total_frc_tol = 5.0e-4;
            vel_tol = 2.5e-6;
            pos_tol = 1.5e-6;
            break;
          case PrecisionModel::SINGLE:
            ngbr_frc_tol = 2.5e-3;
            qgrid_tol = 6.0e-7;
            ugrid_tol = 3.1e-4;
            total_frc_tol = 2.5e-3;
            vel_tol = 5.0e-6;
            pos_tol = 1.6e-6;
            break;
          }
          break;
        }
        stepwisePMELaboratory(pbc_tsm_ph, dyncon, preccon, rngcon, pmecon, gpu, ngbr_frc_tol,
                              qgrid_tol, ugrid_tol, total_frc_tol, vel_tol, pos_tol);
      }
    }
  }

  // Check the GPU implicit solvent dynamics process
  const std::vector<std::string> system_names = { "trpcage", "ala_dipeptide", "med_1", "med_5" };
  TestSystemManager tsm(base_top_name, "top", system_names, base_crd_name, "inpcrd", system_names);
  
  // Try various systems with a single step of dynamics, no thermostat
  singleStepLaboratory(tsm, 0, gpu, PrecisionModel::DOUBLE, ImplicitSolventModel::HCT_GB,
                       ApplyConstraints::NO, 1, 3, 0.02, 300.0, 0.5);
  singleStepLaboratory(tsm, 1, gpu, PrecisionModel::DOUBLE, ImplicitSolventModel::NONE,
                       ApplyConstraints::NO, 1, 3, 0.02, 300.0, 0.5);
  singleStepLaboratory(tsm, 2, gpu, PrecisionModel::DOUBLE, ImplicitSolventModel::NECK_GB,
                       ApplyConstraints::NO, 3, 8, 0.02, 300.0, 0.5);
  singleStepLaboratory(tsm, 3, gpu, PrecisionModel::DOUBLE, ImplicitSolventModel::OBC_GB,
                       ApplyConstraints::NO, 1, 3, 0.02, 300.0, 0.5);
  singleStepLaboratory(tsm, 1, gpu, PrecisionModel::DOUBLE, ImplicitSolventModel::OBC_GB_II,
                       ApplyConstraints::NO, 3, 8, 0.02, 300.0, 0.5);
  singleStepLaboratory(tsm, 0, gpu, PrecisionModel::DOUBLE, ImplicitSolventModel::NECK_GB_II,
                       ApplyConstraints::NO, 3, 1, 0.02, 300.0, 0.5);
  singleStepLaboratory(tsm, 1, gpu, PrecisionModel::DOUBLE, ImplicitSolventModel::NECK_GB_II,
                       ApplyConstraints::YES, 3, 1, 0.02, 300.0, 0.5);

  // Isolate the alanine dipeptide system
  AtomGraph alad_ag = tsm.exportAtomGraph(0);
  const ImplicitSolventModel born_model = ImplicitSolventModel::HCT_GB;
  alad_ag.setImplicitSolventModel(born_model);
  
  // Try a rather weird means of initializing a vector of PhaseSpace objects, to test the
  // first-classness of the PhaseSpace object itself.
  std::vector<PhaseSpace> alad_ps_vec(1, tsm.exportPhaseSpace(0));
  const std::vector<AtomGraph*> alad_ag_vec(1, &alad_ag);
  const std::vector<int> alad_tiling(1, 0);
  AtomGraphSynthesis alad_poly_ag(alad_ag_vec, alad_tiling, ExceptionResponse::WARN, gpu, &timer);
  StaticExclusionMaskSynthesis alad_poly_se(alad_poly_ag.getUniqueTopologies(),
                                            alad_poly_ag.getTopologyIndices());
  alad_poly_ag.loadNonbondedWorkUnits(alad_poly_se, InitializationTask::GB_LANGEVIN_DYNAMICS, 15,
                                      gpu);
  PhaseSpaceSynthesis alad_poly_ps(alad_ps_vec, alad_ag_vec, alad_tiling);
  alad_poly_ag.upload();
  alad_poly_se.upload();
  alad_poly_ps.upload();
  
  // Test the thermostat construction
  Thermostat alad_heat_bath(alad_poly_ag, ThermostatKind::LANGEVIN, 300.0);
  alad_heat_bath.setRandomCacheDepth(alad_poly_ag.getRandomCacheDepth());
  alad_heat_bath.initializeRandomStates(9815734, 25, HybridTargetLevel::DEVICE, gpu);
  alad_heat_bath.setTimeStep(1.0);
  std::vector<double> gpu_result, cpu_result;
  std::vector<ullint4> gpu_gstate, cpu_gstate;
  Xoshiro256ppGenerator xrs(9815734, 25);
  for (int i = 0; i < alad_heat_bath.getAtomCount(); i+= 1024) {
    xrs.setState(alad_heat_bath.getGeneratorState(i, HybridTargetLevel::HOST));
    for (int j = 0; j < alad_heat_bath.getRandomCacheDepth() * 3; j++) {
      gpu_result.push_back(alad_heat_bath.getCachedRandomResult(i, j, HybridTargetLevel::DEVICE));
      cpu_result.push_back(xrs.spGaussianRandomNumber());
    }
    gpu_gstate.push_back(alad_heat_bath.getGeneratorState(i, HybridTargetLevel::DEVICE));
    cpu_gstate.push_back(xrs.revealState());
  }

  // Check dynamics of a small system involving periodic boundary conditions and constraints.
  checkRigidWater(base_crd_name, base_top_name, PrecisionModel::SINGLE);
  checkRigidWater(base_crd_name, base_top_name, PrecisionModel::DOUBLE);
  
  // Check the integration of equations of motion, one kernel at a time and in fused segments
  const std::vector<std::string> pbc_mols = { "tip3p", "tip4p", "tamavidin", "trpcage_in_water" };
  const std::vector<std::string> iso_mols = { "bromobenzene_iso", "drug_example_iso", "trpcage" };
  TestSystemManager pbc_tsm(base_top_name, "top", pbc_mols, base_crd_name, "inpcrd", pbc_mols);
  TestSystemManager iso_tsm(base_top_name, "top", iso_mols, base_crd_name, "inpcrd", iso_mols);
  const std::vector<UnitCellType> iso_bcs = { UnitCellType::NONE };
  const std::vector<UnitCellType> img_bcs = { UnitCellType::ORTHORHOMBIC,
                                              UnitCellType::TRICLINIC };

  // Establish some simple and inexpensive minimization to do to the structures
  MinimizeControls mincon;
  mincon.setClashDampingCycles(50);
  mincon.setSteepestDescentCycles(100);
  mincon.setTotalCycles(200);
  mincon.setCheckpointProduction(false);
  mincon.setDiagnosticPrintFrequency(5);
  
  // Establish some very short dynamics simulations
  DynamicsControls dyncon_pbc;
  dyncon_pbc.setTimeStep(1.0);
  dyncon_pbc.setCutoff(9.0);
  dyncon_pbc.setCenterOfMassMotionPurgeFrequency(0);
  dyncon_pbc.setGeometricConstraints(ApplyConstraints::YES);
  DynamicsControls dyncon_iso = dyncon_pbc;

  // Establish an implicit solvent system for isolate boundary conditions
  SolventControls watcon;
  watcon.setImplicitSolventModel(ImplicitSolventModel::NONE);
  watcon.choosePBRadiiSet(AtomicRadiusSet::MBONDI3);
  piecewiseIntegrationTests(iso_tsm, iso_bcs, mincon, dyncon_iso, watcon, PrecisionModel::DOUBLE,
                            gpu);
  piecewiseIntegrationTests(iso_tsm, iso_bcs, mincon, dyncon_iso, watcon, PrecisionModel::SINGLE,
                            gpu);
#if 0
  piecewiseIntegrationTests(pbc_tsm, img_bcs, mincon, dyncon_pbc, watcon, PrecisionModel::DOUBLE,
                            gpu);
  piecewiseIntegrationTests(pbc_tsm, img_bcs, mincon, dyncon_pbc, watcon, PrecisionModel::SINGLE,
                            gpu);
#endif
  
  // Display timings and test results
  if (oe.getDisplayTimingsOrder()) {
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
  return countGlobalTestFailures();
}
