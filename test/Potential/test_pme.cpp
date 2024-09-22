#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
#    include <cuda.h>
#    include <cuda_runtime.h>
#  endif
#endif
#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Accelerator/core_kernel_manager.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/math_enumerators.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/MolecularMechanics/mm_controls.h"
#include "../../src/Parsing/parsing_enumerators.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Potential/cellgrid.h"
#include "../../src/Potential/energy_enumerators.h"
#ifdef STORMM_USE_HPC
#  include "../../src/Potential/hpc_pme_potential.h"
#endif
#include "../../src/Potential/map_density.h"
#include "../../src/Potential/pme_potential.h"
#include "../../src/Potential/pme_util.h"
#include "../../src/Potential/pmigrid.h"
#include "../../src/Potential/ppitable.h"
#include "../../src/Potential/tile_manager.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/hpc_phasespace_synthesis.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Topology/atomgraph_analysis.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::mm;
using namespace stormm::parse;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::symbols;
using namespace stormm::topology;
using namespace stormm::synthesis;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Check the properties of atoms in the cell grid.
//
// Arguments:
//   cg:       The cell grid, containing an image of all systems with coordinates and properties.
//   poly_ag:  The topology synthesis, containing topology data which can be compared to the
//             cell grid's atom properties.
//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void inspectCellGrids(const CellGrid<T, Tacc, Tcalc, T4> &cg, const AtomGraphSynthesis &poly_ag,
                      const TestPriority do_tests) {
  const CellGridReader<T, Tacc, Tcalc, T4> cgr = cg.data();
  const SyNonbondedKit<double, double2> synbk = poly_ag.getDoublePrecisionNonbondedKit();
  const bool t_is_real = isFloatingPointScalarType<T>();
  const bool t_is_double = (sizeof(T) == 8);
  
  // Count the number of atoms in cells.
  check(cgr.system_count, RelationalOperator::EQUAL, poly_ag.getSystemCount(), "The number of "
        "systems in the cell grid does not equal the number of systems in the associated topology "
        "synthesis", do_tests);
  std::vector<int> atom_counts_ag(cgr.system_count), atom_counts_cg(cgr.system_count);
  std::vector<double> accounting_devs(cgr.system_count), q_presence_devs(cgr.system_count);
  for (int pos = 0; pos < cgr.system_count; pos++) {
    const ullint pcell_dims = cgr.system_cell_grids[pos];
    const int pcell_init = (pcell_dims & 0xfffffffLLU);
    const int pcell_na = ((pcell_dims >> 28) & 0xfffLLU);
    const int pcell_nb = ((pcell_dims >> 40) & 0xfffLLU);
    const int pcell_nc = ((pcell_dims >> 52) & 0xfffLLU);
    const int pncell= pcell_na * pcell_nb * pcell_nc;
    const int imax = pcell_init + pncell;
    int natom = 0;
    const int patom_offset = poly_ag.getAtomOffset(pos);
    const int patom_count  = poly_ag.getAtomCount(pos);
    std::vector<int> accounted(patom_count, 0);
    std::vector<double> q_present(patom_count, 0.0);
    std::vector<double> q_expected(patom_count, 0.0);
    if (cgr.theme != NonbondedTheme::VAN_DER_WAALS) {
      if (t_is_double) {
        q_expected = poly_ag.getSystemTopologyPointer(pos)->getPartialCharge<double>();
      }
      else {
        for (int i = 0; i < patom_count; i++) {
          q_expected[i] = poly_ag.getSystemTopologyPointer(pos)->getPartialCharge<float>(i);
        }
      }
    }
    double tchrg = 0.0;
    for (int i = pcell_init; i < imax; i++) {
      const uint2 icell_limits = cgr.cell_limits[i];
      const uint jlim = icell_limits.x + (icell_limits.y >> 16);
      natom += (icell_limits.y >> 16);
      for (uint j = icell_limits.x; j < jlim; j++) {
        accounted[cgr.nonimg_atom_idx[j] - patom_offset] += 1;
        switch (cgr.theme) {
        case NonbondedTheme::ELECTROSTATIC:
          if (t_is_real) {
            q_present[cgr.nonimg_atom_idx[j] - patom_offset] = cgr.image[j].w;
          }
          else {
            if (t_is_double) {
              const Ecumenical8 conv = { .lli = static_cast<llint>(cgr.image[j].w) };
              q_present[cgr.nonimg_atom_idx[j] - patom_offset] = conv.d;
            }
            else {
              const Ecumenical4 conv = { .i = static_cast<int>(cgr.image[j].w) };
              q_present[cgr.nonimg_atom_idx[j] - patom_offset] = conv.f;
            }
          }
          break;
        case NonbondedTheme::VAN_DER_WAALS:
          q_present[cgr.nonimg_atom_idx[j] - patom_offset] = 0.0;
          break;
        case NonbondedTheme::ALL:
          {
            int q_idx;
            if (t_is_double) {
              q_idx = (static_cast<llint>(cgr.image[j].w) & dp_charge_index_mask);
            }
            else {
              q_idx = (static_cast<int>(cgr.image[j].w) & sp_charge_index_mask);
            }
            q_present[cgr.nonimg_atom_idx[j] - patom_offset] = synbk.q_params[q_idx];
          }
          break;
        }
      }
    }
    atom_counts_cg[pos] = natom;
    std::vector<int> accounted_ans(patom_count, 1);
    switch (cg.getTheme()) {
    case NonbondedTheme::ELECTROSTATIC:
      for (int i = patom_offset; i < patom_offset + patom_count; i++) {
        if (fabs(synbk.charge[i]) < 1.0e-6) {
          accounted_ans[i - patom_offset] = 0;
        }
      }
      break;
    case NonbondedTheme::VAN_DER_WAALS:
      for (int i = patom_offset; i < patom_offset + patom_count; i++) {

        // This will fail if the topology contains NBFix pair-specific Lennard-Jones interactions
        // such that the self-interaction is zero but an off-diagonal term is not.  Such a case is
        // very rare and not what this function is intended to test for.
        const int tlj_idx = synbk.lj_idx[i];
        const int nlj_typ = synbk.n_lj_types[pos];
        const int tsys_offset = synbk.ljabc_offsets[pos];
        if (fabs(synbk.ljab_coeff[tsys_offset + ((nlj_typ + 1) * tlj_idx)].y) < 1.0e-6) {
          accounted_ans[i - patom_offset] = 0;
        }
      }
      break;
    case NonbondedTheme::ALL:
      break;
    }
    atom_counts_ag[pos] = sum<int>(accounted_ans);
    accounting_devs[pos] = meanUnsignedError(accounted, accounted_ans);
    q_presence_devs[pos] =  meanUnsignedError(q_present, q_expected);
  }
  check(atom_counts_cg, RelationalOperator::EQUAL, atom_counts_ag, "The numbers of atoms in each "
        "system of the cell grid do not match the corresponding systems of the topology "
        "synthesis.  Cell grid theme: " + getEnumerationName(cgr.theme) + ".", do_tests);
  check(accounting_devs, RelationalOperator::EQUAL, std::vector<double>(cgr.system_count, 0.0),
        "The accounting of atoms in the cell grid is incomplete or over-stocked.  Cell grid "
        "theme: " + getEnumerationName(cgr.theme) + ".", do_tests);
  check(q_presence_devs, RelationalOperator::EQUAL, std::vector<double>(cgr.system_count, 0.0),
        "The accounting of charge in the cell grid is incorrect.  Cell grid theme: " +
        getEnumerationName(cgr.theme) + ".", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Compare the grid computed within one system of a PMIGrid object, whether using CPU methods or a
// GPU kernel, to that produced by an independent method (the calling function uses a CPU-based
// approach involving only the one system's standard, double-precision coordinates and its one
// topology object).
//
// Arguments:
//   cg:                The cell grid composed for the synthesis of systems.  The necessary
//                      topology synthesis can be extracted as needed.
//   pm:                The set of particle-mesh interaction grids for the same synthesis as was
//                      used to construct cg
//   sysid:             Index of the system of interest within the synthesis
//   q_compare:         Linearized form of the particle-mesh interaction grid for the one system
//   hardware:          Keyword for the hardware responsible for the calculation being checked,
//                      e.g. "CPU" or "GPU"
//   pm_fp_msg:         Message concerning the use of fixed-precision bits in the calculation
//   chrg_value_tol:    Tolerance for the density to deviate at any particular grid point
//   do_tests:          Indicate whether tests are viable
//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void compareSystemDensityField(const CellGrid<T, Tacc, Tcalc, T4> &cg, const PMIGrid &pm,
                               const int sysid, const std::vector<double> &q_compare,
                               const std::string &hardware, const std::string &pm_fp_msg,
                               const double chrg_value_tol, const TestPriority do_tests) {
  const std::vector<double> derived_q_grid = pm.getGrid(sysid);
  const AtomGraph* ag_ptr = pm.getCoordinateSynthesisPointer()->getSystemTopologyPointer(sysid);
  check(derived_q_grid, RelationalOperator::EQUAL, Approx(q_compare).margin(chrg_value_tol),
        "The density grid computed based on a CellGrid and PMIGrid (" + hardware + " computation) "
        "did not agree with a simpler calculation performed in uniform, double precision.  System "
        "topology file: " + getBaseName(ag_ptr->getFileName()) + ".  Cell grid coordinate type: " +
        getStormmScalarTypeName<T>() + ".  Cell grid non-bonded content: " +
        getEnumerationName(cg.getTheme()) + ".  Interpolation order: " +
        std::to_string(pm.getInterpolationOrder()) + ".  Particle-mesh interaction grid precision "
        "mode: " + getEnumerationName(pm.getMode()) + "." + pm_fp_msg, do_tests);
}

//-------------------------------------------------------------------------------------------------
// Finish setting up the particle-mesh interaction grids begun in spatialDecompositionOuter()
// below and perform density spreading.
//
// Arguments:
//   cg:                The cell grid composed for the synthesis of systems.  The necessary
//                      topology synthesis can be extracted as needed.
//   pmi_grid_prec:     Precision with which to represent the particle-mesh interaction grids.
//                      Like cell_grid_prec, this carries over into a fixed-precision
//                      representation given a nonzero setting in pm_acc_bits.
//   pm_acc_bits:       Number of bits after the decimal in which to perform the fixed-precision
//                      accumulation on the particle-mesh interaction grids
//   pm_order:          The order of particle-mesh density interpolation
//   pm_theme:          Indicate the type of non-bonded interaction for which to spread density
//   do_tests:          Indicate whether tests are viable
//   chrg_conserv_tol:  Tolerance for the total density conservation test
//   chrg_value_tol:    Tolerance for the density to deviate at any particular grid point
//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void spatialDecompositionInner(const CellGrid<T, Tacc, Tcalc, T4> &cg,
                               const PrecisionModel pmi_grid_prec, const int pm_acc_bits,
                               const int pm_order, const NonbondedTheme pm_theme,
                               const TestPriority do_tests, const double chrg_conserv_tol,
                               const double chrg_value_tol, const CoreKlManager &launcher) {
  PMIGrid pm(cg, pm_theme, pm_order, pmi_grid_prec, FFTMode::OUT_OF_PLACE, pm_acc_bits);
  const AtomGraphSynthesis *poly_ag = cg.getTopologySynthesisPointer();
  const PhaseSpaceSynthesis *poly_ps = cg.getCoordinateSynthesisPointer();
  mapDensity(&pm, poly_ag);
  std::vector<double> system_total_q(pm.getSystemCount());
  std::vector<double> chrg_conservation(pm.getSystemCount());
  const std::vector<AtomGraph*> system_topologies = poly_ag->getSystemTopologyPointer();
  for (int i = 0; i < pm.getSystemCount(); i++) {
    switch (pm_theme) {
    case NonbondedTheme::ELECTROSTATIC:
      system_total_q[i] = sum<double>(system_topologies[i]->getPartialCharge<double>());
      break;
    case NonbondedTheme::VAN_DER_WAALS:
      {
        double disp_src = 0.0;
        const NonbondedKit<double> i_nbk = system_topologies[i]->getDoublePrecisionNonbondedKit();
        for (int j = 0; j < i_nbk.natom; j++) {
          const double ljb = i_nbk.ljb_coeff[(i_nbk.n_lj_types + 1) * i_nbk.lj_idx[j]];
          if (ljb > 1.0e-6) {
            disp_src += 0.5 * sqrt(ljb);
          }
        }
        system_total_q[i] = disp_src;
      }
      break;
    case NonbondedTheme::ALL:
      break;
    }
    chrg_conservation[i] = pm.getTotalOnGrid(i);
  }

  // Check some features of the CellGrid object
  const CellGridReader<T, Tacc, Tcalc, T4> cgr = cg.data();
  int chain_count = 0;
  int chncon = 0;
  bool cell_starts_match_chains = true;
  for (int i = 0; i < cgr.system_count; i++) {
    const ullint igrid = cgr.system_cell_grids[i];
    const int cell_init = (igrid & 0xfffffff);
    const int cell_na = ((igrid >> 28) & 0x3ff);
    const int cell_nb = ((igrid >> 40) & 0x3ff);
    const int cell_nc = (igrid >> 52);
    chain_count += cell_nb * cell_nc;
    for (int k = 0; k < cell_nc; k++) {
      for (int j = 0; j < cell_nb; j++) {
        const int first_cell = cell_init + (((k * cell_nb) + j) * cell_na);
        cell_starts_match_chains = (cell_starts_match_chains &&
                                    (cgr.cell_limits[first_cell].x == cgr.chain_limits[chncon]));
        chncon++;
      }
    }
  }
  check(chain_count, RelationalOperator::EQUAL, cgr.total_chain_count, "The total number of "
        "chains in a cell grid (theme " + getEnumerationName(cgr.theme) + ") does not meet "
        "expectations.", do_tests);
  check(cell_starts_match_chains, "The starting atom indices of chains in the image do not "
        "correspond to those of their first cells.", do_tests);
  
  // Check the total charge conservation on the mesh.
  const std::string pm_fp_msg = "  Particle-mesh interaction grid fixed-precision bits: " +
                                std::to_string(pm_acc_bits) + ".";
  check(chrg_conservation, RelationalOperator::EQUAL,
        Approx(system_total_q).margin(chrg_conserv_tol), "The total charge on the grids does not "
        "match that in the underlying topologies.  Cell grid coordinate type: " +
        getStormmScalarTypeName<T>() + ".  Cell grid non-bonded content: " +
        getEnumerationName(cg.getTheme()) + ".  Particle-mesh interaction grid precision mode: " +
        getEnumerationName(pm.getMode()) + "." + pm_fp_msg, do_tests);

  // Perform a naive loop over all atoms in each system, setting up a new grid in the process to
  // check the values produced by the more complex CellGrid / PMIGrid approach.  Begin by
  // converting the PMIGrid to its real representation (after accumulation is complete).
  if (pm.dataIsReal() == false) {
    pm.convertToReal();
  }
  const PMIGridReader pm_rdr = pm.data();
  std::vector<std::vector<double>> simple_q_grids(pm.getSystemCount());
  for (int sysid = 0; sysid < pm.getSystemCount(); sysid++) {
    const CoordinateFrame cf = poly_ps->exportCoordinates(sysid);
    const int nx = cg.getMeshSubdivisions();
    const int na = nx * cg.getCellCount(sysid, UnitCellAxis::A);
    const int nb = nx * cg.getCellCount(sysid, UnitCellAxis::B);
    const int nc = nx * cg.getCellCount(sysid, UnitCellAxis::C);
    simple_q_grids[sysid] = mapDensity(&cf, system_topologies[sysid], pm_theme, pm_rdr.fftm, na,
                                       nb, nc, pm_order);
    compareSystemDensityField(cg, pm, sysid, simple_q_grids[sysid], "CPU", pm_fp_msg,
                              chrg_value_tol, do_tests);
  }
#ifdef STORMM_USE_HPC
  // Make a copy of the CellGrid so that its contents can be uploaded to the device memory.  Repeat
  // the tests with HPC kernels being responsible for the mapping.
  if (pm_acc_bits > 0) {
    CellGrid cg_copy = cg;
    cg_copy.upload();
    PMIGrid pm_copy(cg_copy, pm_theme, pm_order, pmi_grid_prec, FFTMode::OUT_OF_PLACE,
                    pm_acc_bits);
    pm_copy.upload();
    
    // The general-purpose density mapping kernel does not require any work units tailored to
    // either the cell grid or the particle-mesh interaction grids, and therefore does not need
    // either object to have been created with the particular GPU specifications as part of the
    // input.
    mapDensity(&pm_copy, nullptr, &cg_copy, poly_ag, launcher, QMapMethod::GENERAL_PURPOSE);
    pm_copy.download();
    for (int sysid = 0; sysid < pm_copy.getSystemCount(); sysid++) {
      chrg_conservation[sysid] = pm_copy.getTotalOnGrid(sysid);
      compareSystemDensityField(cg, pm_copy, sysid, simple_q_grids[sysid], "GPU", pm_fp_msg,
                                chrg_value_tol, do_tests);
    }
    check(chrg_conservation, RelationalOperator::EQUAL,
          Approx(system_total_q).margin(chrg_conserv_tol), "The total charge on grids computed by "
          "a GPU kernel (" + getEnumerationName(QMapMethod::GENERAL_PURPOSE) + ") does not match "
          "that in the underlying topologies.  Cell grid coordinate type: " +
          getStormmScalarTypeName<T>() + ".  Cell grid non-bonded content: " +
          getEnumerationName(cg.getTheme()) + ".  Interpolation order: " +
          std::to_string(pm_order) + ".  Particle-mesh interaction grid precision mode: " +
          getEnumerationName(pm.getMode()) + "." + pm_fp_msg, do_tests);
  }

  // In all cases, try the __shared__ accumulation kernel.
  CellGrid cg_copy_ii = cg;
  cg_copy_ii.upload();
  PMIGrid pm_copy_ii(cg_copy_ii, pm_theme, pm_order, pmi_grid_prec, FFTMode::OUT_OF_PLACE,
                     pm_acc_bits);
  pm_copy_ii.prepareWorkUnits(QMapMethod::ACC_SHARED, launcher.getGpu());
  pm_copy_ii.setRealDataFormat();
  pm_copy_ii.upload();
  MolecularMechanicsControls mm_ctrl;
  const PrecisionModel tcalc_prec = (std::type_index(typeid(Tcalc)).hash_code() ==
                                     double_type_index) ?
                                    PrecisionModel::DOUBLE : PrecisionModel::SINGLE;
  const size_t cg_tmat = std::type_index(typeid(T)).hash_code();
  mm_ctrl.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                                ClashResponse::NONE, VwuGoal::ACCUMULATE, tcalc_prec, tcalc_prec,
                                pm_copy_ii.getWorkUnitConfiguration(), pm_copy_ii.getMode(),
                                cg_tmat, pm_order, NeighborListKind::MONO, TinyBoxPresence::NO,
                                *poly_ag);
  mapDensity(&pm_copy_ii, &mm_ctrl, &cg_copy_ii, poly_ag, launcher, QMapMethod::ACC_SHARED);
  pm_copy_ii.download();
  for (int sysid = 0; sysid < pm_copy_ii.getSystemCount(); sysid++) {
    chrg_conservation[sysid] = pm_copy_ii.getTotalOnGrid(sysid);
    compareSystemDensityField(cg, pm_copy_ii, sysid, simple_q_grids[sysid],
                              "GPU, __shared__ accumulation kernel", pm_fp_msg,
                              chrg_value_tol, do_tests);
  }
  check(chrg_conservation, RelationalOperator::EQUAL,
        Approx(system_total_q).margin(chrg_conserv_tol), "The total charge on grids computed by "
        "a GPU kernel (" + getEnumerationName(QMapMethod::ACC_SHARED) + ") does not match "
        "that in the underlying topologies.  Cell grid coordinate type: " +
        getStormmScalarTypeName<T>() + ".  Cell grid non-bonded content: " +
        getEnumerationName(cg.getTheme()) + ".  Interpolation order: " +
        std::to_string(pm_order) + ".  Particle-mesh interaction grid precision mode: " +
        getEnumerationName(pm.getMode()) + "." + pm_fp_msg, do_tests);
#endif
}

//-------------------------------------------------------------------------------------------------
// Set up a series of spatial decompositions and particle-mesh interaction grids for the collection
// of systems at hand, given some directives for the precision of each component.  This outer
// function will set up the cell grid, then call the inner function above to set up the
// particle-mesh interaction grid.
//
// Arguments:
//   poly_ag:           The topology synthesis
//   poly_ps:           The coordinate synthesis matching poly_ag
//   cell_grid_prec:    The precision model in which to create the cell grid.  Given a nonzero
//                      value of cg_crd_bits, this can give rise to fixed-precision representations
//                      of the local coordinate and decomposition cell dimensions.
//   cast_cg_in_fp:     Flag to have the cell grid coordinates represented in fixed precision.  The
//                      number of bits after the decimal will be determined automatically by the
//                      dimensions of each spatial decomposition cell and the size of the data type
//                      indicated by cell_grid_prec.
//   mesh_ticks:        The number of particle-mesh interaction elements spanning the spatial
//                      decomposition cells in all directions (i.e. 4 will subdivide each spatial
//                      decomposition cell into 64 grid elements for the charge spreading)
//   cg_theme:          The nature of non-bonded interactions for which the spatial decomposition
//                      is prepared
//   pmi_grid_prec:     Precision with which to represent the particle-mesh interaction grids.
//                      Like cell_grid_prec, this carries over into a fixed-precision
//                      representation given a nonzero setting in pm_acc_bits.
//   pm_acc_bits:       Number of bits after the decimal in which to perform the fixed-precision
//                      accumulation on the particle-mesh interaction grids
//   pm_order:          The order of particle-mesh density interpolation
//   pm_theme:          Indicate the type of non-bonded interaction for which to spread density
//   do_tests:          Indicate whether tests are viable
//   chrg_conserv_tol:  Tolerance for the total charge conservation test
//   chrg_value_tol:    Tolerance for the density to deviate at any particular grid point
//-------------------------------------------------------------------------------------------------
void spatialDecompositionOuter(const AtomGraphSynthesis &poly_ag,
                               const PhaseSpaceSynthesis &poly_ps,
                               const PrecisionModel cell_grid_prec, const bool cast_cg_in_fp,
                               const int mesh_ticks, const NonbondedTheme cg_theme,
                               const PrecisionModel pmi_grid_prec, const int pm_acc_bits,
                               const int pm_order, const NonbondedTheme pm_theme,
                               const TestPriority do_tests, const double chrg_conserv_tol,
                               const double chrg_value_tol, const CoreKlManager &launcher) {
  switch (cell_grid_prec) {
  case PrecisionModel::DOUBLE:
    if (cast_cg_in_fp) {
      section(1);
      CellGrid<llint, llint, double, llint4> cg(poly_ps, poly_ag, 4.5, 0.25, mesh_ticks, cg_theme);
      inspectCellGrids(cg, poly_ag, do_tests);
      section(2);
      spatialDecompositionInner<llint, llint, double, llint4>(cg, pmi_grid_prec, pm_acc_bits,
                                                              pm_order, pm_theme, do_tests,
                                                              chrg_conserv_tol, chrg_value_tol,
                                                              launcher);
    }
    else {
      section(1);
      CellGrid<double, llint, double, double4> cg(poly_ps, poly_ag, 4.5, 0.25, mesh_ticks,
                                                  cg_theme);
      inspectCellGrids(cg, poly_ag, do_tests);
      section(2);
      spatialDecompositionInner<double, llint, double, double4>(cg, pmi_grid_prec, pm_acc_bits,
                                                                pm_order, pm_theme, do_tests,
                                                                chrg_conserv_tol, chrg_value_tol,
                                                                launcher);
    }
    break;
  case PrecisionModel::SINGLE:
    if (cast_cg_in_fp) {
      section(1);
      CellGrid<int, int, float, int4> cg(poly_ps, poly_ag, 4.5, 0.25, mesh_ticks, cg_theme);
      inspectCellGrids(cg, poly_ag, do_tests);
      section(2);
      spatialDecompositionInner<int, int, float, int4>(cg, pmi_grid_prec, pm_acc_bits, pm_order,
                                                       pm_theme, do_tests, chrg_conserv_tol,
                                                       chrg_value_tol, launcher);
    }
    else {
      section(1);
      CellGrid<float, int, float, float4> cg(poly_ps, poly_ag, 4.5, 0.25, mesh_ticks, cg_theme);
      inspectCellGrids(cg, poly_ag, do_tests);
      section(2);
      spatialDecompositionInner<float, int, float, float4>(cg, pmi_grid_prec, pm_acc_bits,
                                                           pm_order, pm_theme, do_tests,
                                                           chrg_conserv_tol, chrg_value_tol,
                                                           launcher);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
// Compute non-bonded forces in a PhaseSpace object broken down into a neighbor list decomposition
// using a simple all-to-all method.  The non-bonded accumulated energy, with the electrostatic
// sum in the "x" member and the van-der Waals sum in the "y" member, is returned.
//
// Arguments:
//   ps_chk:       Pre-allocated copy of the coordinate set in question (forces will be initialized
//                 and then accumulated)
//   ag:           Topology for the system in question
//   elec_cutoff:  Cutoff for electrostatic interactions
//   vdw_cutoff:   Cutoff for Lennard-Jones interactions
//   qqew_coeff:   Splitting coefficient for electrostatic interactions
//   ljew_coeff:   Splitting coefficient for Lennard-Jones interactions
//   vdw_sum:      The method for computing the particle-particle Lennard-Jones sum
//   lema:         Exclusions within the synthesis of systems
//   offset:       The offset of atom indices by which to look up exclusions
//   theme:        The type of non-bonded calculations to evaluate
//   do_tests:     Indicate whether tests are possible
//-------------------------------------------------------------------------------------------------
double2 cutoffAllToAll(PhaseSpace *ps_chk, const AtomGraph &ag, const double elec_cutoff,
                       const double vdw_cutoff, const double qqew_coeff, const double ljew_coeff,
                       const VdwSumMethod vdw_sum, const LocalExclusionMask &lema,
                       const int offset, const NonbondedTheme theme, const TestPriority do_tests) {
  ps_chk->initializeForces();
  PhaseSpaceWriter chkw = ps_chk->data();
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  const double elec_cutsq = elec_cutoff * elec_cutoff;
  const double vdw_cutsq = vdw_cutoff * vdw_cutoff;
  const double qq_bfac = 2.0 * qqew_coeff / sqrt(pi);
  double2 chksum = { 0.0, 0.0 };
  for (int i = 1; i < chkw.natom; i++) {
    const double atmi_x = chkw.xcrd[i];
    const double atmi_y = chkw.ycrd[i];
    const double atmi_z = chkw.zcrd[i];
    const double atmi_q = nbk.charge[i] * nbk.coulomb_constant;
    const int atmi_ljidx = nbk.lj_idx[i] * nbk.n_lj_types;
    for (int j = 0; j < i; j++) {
      double dx = chkw.xcrd[j] - atmi_x;
      double dy = chkw.ycrd[j] - atmi_y;
      double dz = chkw.zcrd[j] - atmi_z;
      imageCoordinates<double, double>(&dx, &dy, &dz, chkw.umat, chkw.invu, chkw.unit_cell,
                                       ImagingMethod::MINIMUM_IMAGE, 1.0);
      const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
      const double invr2 = 1.0 / r2;
      const double qij = atmi_q * nbk.charge[j];
      const int ljidx_ij = atmi_ljidx + nbk.lj_idx[j];
      double fmag = 0.0;
      if (lema.testExclusion(offset + i, offset + j)) {
        switch (theme) {
        case NonbondedTheme::ELECTROSTATIC:
        case NonbondedTheme::ALL:
          {
            const double invr = sqrt(invr2);
            const double u_quant = erfc(qqew_coeff / invr) * invr;
            const double exp_quant = qq_bfac * exp(-qqew_coeff * qqew_coeff * r2);
            fmag += -qij * (((exp_quant + u_quant) * invr) - invr2) * invr;
            chksum.x += qij * (u_quant - invr);
          }
          break;
        case NonbondedTheme::VAN_DER_WAALS:
          break;
        }
        switch (theme) {
        case NonbondedTheme::VAN_DER_WAALS:
        case NonbondedTheme::ALL:
          switch (vdw_sum) {
          case VdwSumMethod::CUTOFF:
          case VdwSumMethod::SMOOTH:
            break;
          case VdwSumMethod::PME:

            // If an infinite sum is in effect for PME van-der Waals interactions, then the primary
            // image of the inverse r^6 interaction must be subtracted.  Otherwise, the exclusion
            // implies that no interaction will be calculated.
            chksum.y += nbk.ljb_coeff[ljidx_ij] * invr2 * invr2 * invr2;
            break;
          }
          break;
        case NonbondedTheme::ELECTROSTATIC:
          break;
        }
      }
      else {
        switch (theme) {
        case NonbondedTheme::ELECTROSTATIC:
        case NonbondedTheme::ALL:
          if (r2 < elec_cutsq) {
            const double invr = sqrt(invr2);
            const double u_quant = erfc(qqew_coeff / invr) * invr;
            const double exp_quant = qq_bfac * exp(-qqew_coeff * qqew_coeff * r2);
            fmag += -qij * (exp_quant + u_quant) * invr2;
            chksum.x += qij * u_quant;
          }
          break;
        case NonbondedTheme::VAN_DER_WAALS:
          break;
        }
        switch (theme) {
        case NonbondedTheme::VAN_DER_WAALS:
        case NonbondedTheme::ALL:
          if (r2 < vdw_cutsq) {
            switch (vdw_sum) {
            case VdwSumMethod::CUTOFF:
              {
                const double invr4 = invr2 * invr2;
                const double invr6 = invr4 * invr2;
                const double lja = nbk.lja_coeff[ljidx_ij];
                const double ljb = nbk.ljb_coeff[ljidx_ij];
                fmag += ((6.0 * ljb) - (12.0 * lja * invr6)) * invr4 * invr4;
                chksum.y += ((lja * invr6) - ljb) * invr6;
              }
              break;
            case VdwSumMethod::PME:
            case VdwSumMethod::SMOOTH:
              break;
            }
          }
          break;
        case NonbondedTheme::ELECTROSTATIC:
          break;
        }
      }
      const double fmag_dx = dx * fmag;
      const double fmag_dy = dy * fmag;
      const double fmag_dz = dz * fmag;
      chkw.xfrc[i] += fmag_dx;
      chkw.yfrc[i] += fmag_dy;
      chkw.zfrc[i] += fmag_dz;
      chkw.xfrc[j] -= fmag_dx;
      chkw.yfrc[j] -= fmag_dy;
      chkw.zfrc[j] -= fmag_dz;
    }
  }
  return chksum;
}

//-------------------------------------------------------------------------------------------------
// Verify the accuracy of particle-particle interactions computed using a basic neighbor list with
// the simplest possible function, one which uses no neighbor list.  Descriptions of input
// parameters follow from cutoffAllToAll(), above, in addition to:
//
// Arguments:
//   ps:           The coordinates set in question, containing positions and accumulated forces for
//                 all particles
//   nrg:          Non-bonded energy sum associated with ps, with the electrostatic sum in the "x"
//                 member of the tuple and the Lennard-Jones sum in the "y" member of the tuple
//-------------------------------------------------------------------------------------------------
void checkParticlePairInteractions(PhaseSpace *ps_chk, const PhaseSpace &ps, const AtomGraph &ag,
                                   const double elec_cutoff, const double vdw_cutoff,
                                   const double qqew_coeff, const double ljew_coeff,
                                   const VdwSumMethod vdw_sum, const LocalExclusionMask &lema,
                                   const int offset, const double2 nrg, const NonbondedTheme theme,
                                   const TestPriority do_tests) {
  const double2 chksum = cutoffAllToAll(ps_chk, ag, elec_cutoff, vdw_cutoff, qqew_coeff,
                                        ljew_coeff, vdw_sum, lema, offset, theme, do_tests);

  // Check that the cell-based sum agrees with the simple procedure for these individual systems.
  bool check_elec = true;
  bool check_vdw = true;
  switch (theme) {
  case NonbondedTheme::ELECTROSTATIC:
    check_vdw = false;
    break;
  case NonbondedTheme::VAN_DER_WAALS:
    check_elec = false;
    break;
  case NonbondedTheme::ALL:
    break;
  }
  PhaseSpaceWriter chkw = ps_chk->data();
  std::vector<double> ref_fx(chkw.natom), ref_fy(chkw.natom), ref_fz(chkw.natom);
  const PhaseSpaceReader psr = ps.data();
  std::vector<double> tst_fx(psr.natom), tst_fy(psr.natom), tst_fz(psr.natom);
  for (int i = 0; i < psr.natom; i++) {
    ref_fx[i] = chkw.xfrc[i];
    ref_fy[i] = chkw.yfrc[i];
    ref_fz[i] = chkw.zfrc[i];
    tst_fx[i] = psr.xfrc[i];
    tst_fy[i] = psr.yfrc[i];
    tst_fz[i] = psr.zfrc[i];
  }
  if (check_elec) {
    check(chksum.x, RelationalOperator::EQUAL, nrg.x, "The electrostatic energy sums for a system "
          "based on topology " + getBaseName(ag.getFileName()) + " do not agree when computed "
          "with a simple nested loop versus a basic neighbor list and cell decomposition.  Energy "
          "computation protocol: " + getEnumerationName(theme) + ".", do_tests);
  }
  if (check_vdw) {
    check(chksum.y, RelationalOperator::EQUAL, nrg.y, "The van-der Waals energy sums for a system "
          "based on topology " + getBaseName(ag.getFileName()) + " do not agree when computed "
          "with a simple nested loop versus a basic neighbor list and cell decomposition.  Energy "
          "computation protocol: " + getEnumerationName(theme) + ".", do_tests);
  }
  check(tst_fx, RelationalOperator::EQUAL, ref_fx, "Cartesian X forces for a system based on "
        "topology " + getBaseName(ag.getFileName()) + " do not agree when computed with a simple "
        "nested loop versus a basic neighbor list and cell decomposition.  Energy computation "
        "protocol: " + getEnumerationName(theme) + ".", do_tests);
  check(tst_fy, RelationalOperator::EQUAL, ref_fy, "Cartesian Y forces for a system based on "
        "topology " + getBaseName(ag.getFileName()) + " do not agree when computed with a simple "
        "nested loop versus a basic neighbor list and cell decomposition.  Energy computation "
        "protocol: " + getEnumerationName(theme) + ".", do_tests);
  check(tst_fz, RelationalOperator::EQUAL, ref_fz, "Cartesian Z forces for a system based on "
        "topology " + getBaseName(ag.getFileName()) + " do not agree when computed with a simple "
        "nested loop versus a basic neighbor list and cell decomposition.  Energy computation "
        "protocol: " + getEnumerationName(theme) + ".", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Examine cell grid contents in order to make an additional checkon the object's performance and
// presentation of information.
//
// Arguments:
//   cg:        The cell grid in question, containing an indication of its non-bonded information
//   do_tests:  Indicate whether tests are possible
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void checkCellGridContents(const CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> &cg,
                           const TestPriority do_tests) {
  const PhaseSpaceSynthesis *poly_ps = cg.getCoordinateSynthesisPointer();
  const AtomGraphSynthesis *poly_ag = cg.getTopologySynthesisPointer();
  const CellGridReader<Tcoord, Tacc, Tcalc, Tcoord4> cgr = cg.data();
  const SyNonbondedKit<double, double2> poly_nbk = poly_ag->getDoublePrecisionNonbondedKit();
  const bool tcoord_is_real = isFloatingPointScalarType<Tcoord>();
  
  // Step through each system and determine, based on the type of non-bonded interactions, what
  // should be in the cell grid.
  TinyBoxPresence chk_tiny_box = TinyBoxPresence::NO;
  for (int i = 0; i < poly_ps->getSystemCount(); i++) {
    
    // Trust the cell grid for the spatial decompositon dimensions
    const ullint cell_layout = cgr.system_cell_grids[i];
    const int cell_start = (cell_layout & 0xfffffff);
    const int ncell_a = ((cell_layout >> 28) & 0xfff);
    const int ncell_b = ((cell_layout >> 40) & 0xfff);
    const int ncell_c = (cell_layout >> 52);
    if (ncell_a == 4 || ncell_b == 4 || ncell_c == 4) {
      chk_tiny_box = TinyBoxPresence::YES;
    }
    const double dnc_a = ncell_a;
    const double dnc_b = ncell_b;
    const double dnc_c = ncell_c;
    const int total_cells = ncell_a * ncell_b * ncell_c;
    const PhaseSpace ps = poly_ps->exportSystem(i);
    const AtomGraph *ag = poly_ps->getSystemTopologyPointer(i);
    const PhaseSpaceReader psr = ps.data();
    const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
    std::vector<int> cell_counts(total_cells, 0);
    std::vector<int3> cell_locations(psr.natom);
    for (int j = 0; j < psr.natom; j++) {
      bool incl_atom;
      switch (cg.getTheme()) {
      case NonbondedTheme::ELECTROSTATIC:
        incl_atom = fabs(nbk.charge[j]) > stormm::constants::small;
        break;
      case NonbondedTheme::VAN_DER_WAALS:
        {
          const VdwCombiningRule lj_rule = inferCombiningRule(ag);
          incl_atom = hasVdwProperties<double>(nbk, j, lj_rule);
        }
        break;
      case NonbondedTheme::ALL:
        incl_atom = true;
        break;
      }
      if (incl_atom) {
        double frac_x = (psr.umat[0] * psr.xcrd[j]) + (psr.umat[3] * psr.ycrd[j]) +
                        (psr.umat[6] * psr.zcrd[j]);
        double frac_y = (psr.umat[1] * psr.xcrd[j]) + (psr.umat[4] * psr.ycrd[j]) +
                        (psr.umat[7] * psr.zcrd[j]);
        double frac_z = (psr.umat[2] * psr.xcrd[j]) + (psr.umat[5] * psr.ycrd[j]) +
                        (psr.umat[8] * psr.zcrd[j]);
        frac_x -= floor(frac_x);
        frac_y -= floor(frac_y);
        frac_z -= floor(frac_z);
        frac_x *= dnc_a;
        frac_y *= dnc_b;
        frac_z *= dnc_c;
        const int jcell_a = frac_x;
        const int jcell_b = frac_y;
        const int jcell_c = frac_z;
        cell_counts[(((jcell_c * ncell_b) + jcell_b) * ncell_a) + jcell_a] += 1;
        cell_locations[j] = { jcell_a, jcell_b, jcell_c };
      }
      else {
        cell_locations[j] = { -1, -1, -1 };
      }
    }
    std::vector<std::vector<int>> chk_contents(total_cells);
    std::vector<std::vector<double>> chk_charges(total_cells);
    std::vector<std::vector<int>> chk_ljtypes(total_cells);
    for (int j = 0; j < total_cells; j++) {
      chk_contents[j].reserve(cell_counts[j]);
      chk_charges[j].reserve(cell_counts[j]);
      chk_ljtypes[j].reserve(cell_counts[j]);
    }
    for (int j = 0; j < psr.natom; j++) {
      if (cell_locations[j].x >= 0) {
        const int jcell_loc = (((cell_locations[j].z * ncell_b) + cell_locations[j].y) * ncell_a) +
                              cell_locations[j].x;
        chk_contents[jcell_loc].push_back(j);
        chk_charges[jcell_loc].push_back(nbk.charge[j]);
        chk_ljtypes[jcell_loc].push_back(nbk.lj_idx[j]);
      }
    }
    
    // Step through each cell and confirm the identity and location of each atom.  Log any missing
    // atoms in an array of tuples with the topological atom ID in the "x" member, the cell in
    // which it should be found in the "y" member, and the cell in which it is found in the "z"
    // member (or -1 if the atom is not present at all).  Also check that the total number of atoms
    // in each cell matches what is reported by the cell grid, using another list of tuples with
    // the cell index in the "x" member, the number of atoms displayed by the cell grid in the "y"
    // member, and the number of atoms counted from the topology in the "z" member.
    std::vector<int3> cell_miscounts;
    std::vector<int3> missing_atoms;
    std::vector<double4> erroneous_charges;
    std::vector<int4> erroneous_ljtypes;
    const int synth_atom_offset = poly_ps->getAtomOffset(i);
    for (int j = 0; j < total_cells; j++) {
      const int cg_jcount = (cgr.cell_limits[cell_start + j].y >> 16);
      const uint cg_llim = cgr.cell_limits[cell_start + j].x;
      const uint cg_hlim = cg_llim + cg_jcount;
      if (cg_jcount != cell_counts[j]) {
        cell_miscounts.push_back({ j, cg_jcount, cell_counts[j] });
      }
      for (int k = 0; k < cell_counts[j]; k++) {
        const int atom_idx = chk_contents[j][k];
        const uint img_atom_idx = cgr.img_atom_idx[synth_atom_offset + atom_idx];
        bool found = false;
        for (uint m = cg_llim; m < cg_hlim; m++) {
          found = (found || cgr.nonimg_atom_idx[m] - synth_atom_offset == atom_idx);
        }
        if (found == false) {
          const uint min_sys_atom_idx = cgr.cell_limits[cell_start].x;
          const uint img_idx_delta = img_atom_idx - min_sys_atom_idx;
          const uint chain_length = ncell_a * cgr.cell_base_capacity;
          const uint sys_chain_loc = (img_idx_delta / chain_length);
          const uint pcell_c = sys_chain_loc / ncell_b;
          const uint pcell_b = sys_chain_loc -  (pcell_c * ncell_b);
          int actual_cell = -1;
          for (int m = 0; m < ncell_a; m++) {
            const int test_cell = cell_start + (((pcell_c * ncell_b) + pcell_b) * ncell_a);
            const uint2 tclim = cgr.cell_limits[test_cell];
            if (img_atom_idx >= tclim.x && img_atom_idx < tclim.x + (tclim.y >> 16)) {
              actual_cell = test_cell - cell_start;
            }
          }
          missing_atoms.push_back({ atom_idx, j, actual_cell });
        }
        else {

          // If the atom was found as expected, compare its non-bonded properties.
          double charge;
          int lj_type_idx;
          bool read_charge = true;
          bool read_ljtype = true;
          switch (cgr.theme) {
          case NonbondedTheme::ELECTROSTATIC:
            read_ljtype = false;
            break;
          case NonbondedTheme::VAN_DER_WAALS:
            read_charge = false;
            break;
          case NonbondedTheme::ALL:
            break;
          }
          if (read_charge) {
            charge = sourceMagnitude(NonbondedTheme::ELECTROSTATIC, cgr.theme,
                                     cgr.image[img_atom_idx].w, tcoord_is_real, i, poly_nbk);
            if (fabs(charge - chk_charges[j][k]) > 1.0e-6) {
              erroneous_charges.push_back({ charge, chk_charges[j][k],
                                            static_cast<double>(j), static_cast<double>(k) });
            }
          }
          if (read_ljtype) {
            lj_type_idx = sourceIndex(NonbondedTheme::VAN_DER_WAALS, cgr.theme,
                                      cgr.image[img_atom_idx].w, tcoord_is_real);
            if (lj_type_idx != chk_ljtypes[j][k]) {
              erroneous_ljtypes.push_back({ lj_type_idx, chk_ljtypes[j][k], j, k });
            }
          }
        }
      }
    }
    std::string xmpl_str, chrg_str, lj_str;
    const int nmisc_rep = std::min(static_cast<int>(cell_miscounts.size()), 8); 
    for (int j = 0; j < nmisc_rep; j++) {
      const int xcl_c = cell_miscounts[j].x / (ncell_a * ncell_b);
      const int xcl_b = (cell_miscounts[j].x - (xcl_c * ncell_a * ncell_b)) / ncell_a;
      const int xcl_a = cell_miscounts[j].x - (((xcl_c * ncell_b) + xcl_b) * ncell_a);
      xmpl_str += "[ Cell " + std::to_string(xcl_a) + " " + std::to_string(xcl_b) + " " +
                  std::to_string(xcl_c) + ", grid has " + std::to_string(cell_miscounts[j].y) +
                  ", counted " + std::to_string(cell_miscounts[j].z) + " ]";
      xmpl_str += listSeparator(j, nmisc_rep);
    }
    const int nchrg_rep = std::min(static_cast<int>(erroneous_charges.size()), 8);
    for (int j = 0; j < nchrg_rep; j++) {
      const int cnum = erroneous_charges[j].z;
      const int xcl_c = cnum / (ncell_a * ncell_b);
      const int xcl_b = (cnum - (xcl_c * ncell_a * ncell_b)) / ncell_a;
      const int xcl_a = cnum - (((xcl_c * ncell_b) + xcl_b) * ncell_a);
      chrg_str += "[ Cell " + std::to_string(xcl_a) + " " + std::to_string(xcl_b) + " " +
                  std::to_string(xcl_c) + ", grid " + std::to_string(erroneous_charges[j].x) +
                  ", topology " + std::to_string(erroneous_charges[j].y) + " ]";
      chrg_str += listSeparator(j, nchrg_rep);
    }
    const int nlj_rep = std::min(static_cast<int>(erroneous_ljtypes.size()), 8);
    for (int j = 0; j < nlj_rep; j++) {
      const int cnum = erroneous_ljtypes[j].z;
      const int xcl_c = cnum / (ncell_a * ncell_b);
      const int xcl_b = (cnum - (xcl_c * ncell_a * ncell_b)) / ncell_a;
      const int xcl_a = cnum - (((xcl_c * ncell_b) + xcl_b) * ncell_a);
      lj_str += "[ Cell " + std::to_string(xcl_a) + " " + std::to_string(xcl_b) + " " +
                std::to_string(xcl_c) + ", grid " + std::to_string(erroneous_ljtypes[j].x) +
                ", topology " + std::to_string(erroneous_ljtypes[j].y) + " ]";
      lj_str += listSeparator(j, nchrg_rep);
    }
    check(cell_miscounts.size() == 0, "Some cells of the grid were found to have different "
          "numbers of atoms than were found by a simple re-analysis.  System: " +
          getBaseName(ag->getFileName()) + ".  Non-bonded theme: " +
          getEnumerationName(cg.getTheme()) + ".  Examples of cells include " + xmpl_str + ".",
          do_tests);
    check(missing_atoms.size(), RelationalOperator::EQUAL, 0, "Some atoms of the topology were "
          "not found on the grid as expected.  System: " + getBaseName(ag->getFileName()) +
          ".  Non-bonded theme: " + getEnumerationName(cg.getTheme()) + ".", do_tests);
    check(erroneous_charges.size(), RelationalOperator::EQUAL, 0, "Some charges of the topology "
          "were not properly conveyed by the cell grid, based on a simple re-analysis.  "
          "Examples of discrepancies include " + chrg_str + ".", do_tests);
    check(erroneous_ljtypes.size(), RelationalOperator::EQUAL, 0, "Some Lennard-Jones atom types "
          "of the topology were not properly conveyed by the cell grid, based on a simple "
          "re-analysis.  Examples of discrepancies include " + lj_str + ".", do_tests);
  }

  // Check the representation of the synthesis
  check(chk_tiny_box == cg.getTinyBoxPresence(), "The mannual detection of a tiny box based on "
        "the lengths of systems in the neighbor list decomposition does not match what the "
        "CellGrid object reports.", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Compare two force arrays.  This encapsulates the error message to elimiate repetitive code.
//
// Arguments:
//   test_frc:  The array of forces to test
//   ref_frc:   The array of forces to be trusted
//   dim:       The Cartesian dimension along which forces apply to atoms
//   nrg_kind:  The type of potential giving rise to the forces
//   do_tests:  Indicate whether testing is possible
//-------------------------------------------------------------------------------------------------
void compareForceArrays(const std::vector<double> &test_frc, const std::vector<double> &ref_frc,
                        const CartesianDimension dim, const NonbondedTheme nrg_kind,
                        const TestPriority do_tests) {
  check(test_frc, RelationalOperator::EQUAL, Approx(ref_frc).margin(7.0e-6), "Cartesian " +
        getEnumerationName(dim) + " forces emerging from " + getEnumerationName(nrg_kind) +
        " interactions computed with a CellGrid neighbor list and later transferred to a "
        "PhaseSpaceSynthesis do not agree with those obtained from a simple nested loop "
        "calculation.", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Perform finite difference tests on the forces for small systems.
//
// Arguments:
//   tsm:                  Collection of test systems
//   nbkinds:              Various types of neighbor lists to test
//   critical_atom_count:  The upper limit of system size to test
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void finiteDifferenceNeighborListTest(const TestSystemManager &tsm,
                                      const std::vector<NonbondedTheme> &nbkinds,
                                      const int critical_atom_count = 1024,
                                      const int fd_perturbation_bits = 16) {
  const double ew_coeff = ewaldCoefficient(default_pme_cutoff, default_dsum_tol);
  const std::vector<int> little_indices = tsm.getQualifyingSystems(1024, RelationalOperator::LE);
  PhaseSpaceSynthesis little_systems = tsm.exportPhaseSpaceSynthesis(little_indices);
  PsSynthesisWriter little_psw = little_systems.data();
  AtomGraphSynthesis little_topologies = tsm.exportAtomGraphSynthesis(little_indices);
  const double pos_incr = pow(2.0, -fd_perturbation_bits);
  const int95_t p_pert = hostDoubleToInt95(0.5 * pos_incr * little_psw.gpos_scale_f);
  const int95_t n_pert = hostDoubleToInt95(-pos_incr * little_psw.gpos_scale_f);
  const int95_t z_pert = { 0LL, 0 };
  const std::vector<int95_t> perturbations = { z_pert, z_pert, z_pert, p_pert, z_pert, z_pert,
                                               n_pert, z_pert, z_pert, p_pert, p_pert, z_pert,
                                               z_pert, n_pert, z_pert, z_pert, p_pert, p_pert,
                                               z_pert, z_pert, n_pert, z_pert, z_pert, p_pert };
  std::vector<int> target_atoms(little_indices.size());
  const SyNonbondedKit<double,
                       double2> little_nbk = little_topologies.getDoublePrecisionNonbondedKit();
  for (size_t i = 0; i < nbkinds.size(); i++) {
    const int little_nsys = little_systems.getSystemCount();
    for (int j = 0; j < little_nsys; j++) {
      
      // Obtain the basic topology and limits of the system in the coordinate synthesis
      int n_qual = 0;
      int k = little_nbk.atom_offsets[j];
      const int k_hlim = little_nbk.atom_offsets[j] + little_nbk.atom_counts[j];
      const int ljt_offset = little_nbk.ljabc_offsets[j];
      const int nlj_types = little_nbk.n_lj_types[j];
      while (k < k_hlim && n_qual < 10) {
        const int ljidx = ljt_offset + (little_nbk.lj_idx[k] * (nlj_types + 1));

        // Slightly different qualifications on whether to choose this atom than are used to
        // determine whether to place it in the cell grid.
        switch (nbkinds[i]) {
        case NonbondedTheme::ELECTROSTATIC:
          n_qual += (fabs(little_nbk.charge[k]) > 1.0e-6);
          break;
        case NonbondedTheme::VAN_DER_WAALS:
          n_qual += (little_nbk.ljab_coeff[ljidx].y > 1.0e-2);
          break;
        case NonbondedTheme::ALL:
          n_qual += (fabs(little_nbk.charge[k]) > 1.0e-6 &&
                     little_nbk.ljab_coeff[ljidx].y > 1.0e-2);
          break;
        }
        k++;
      }
      target_atoms[j] = (n_qual == 10) ? k - 1 : little_nbk.atom_offsets[j] +
                                                 (little_nbk.atom_counts[j] / 2);
    }
    std::vector<double> chosen_particle_forces(3 * little_nsys);
    std::vector<double> chosen_particle_fdfrc(3 * little_nsys);
    std::vector<std::vector<double>> chosen_particle_fdnrg(little_nsys,
                                                           std::vector<double>(7, 0.0));
    for (int np = 0; np < 7; np++) {

      // Perturb the chosen particle and calculate energies
      for (int j = 0; j < little_nsys; j++) {
        const int95_t adj_x = hostSplitFPSum(perturbations[ 3 * np     ],
                                             little_psw.xcrd[target_atoms[j]],
                                             little_psw.xcrd_ovrf[target_atoms[j]]);
        const int95_t adj_y = hostSplitFPSum(perturbations[(3 * np) + 1],
                                             little_psw.ycrd[target_atoms[j]],
                                             little_psw.ycrd_ovrf[target_atoms[j]]);
        const int95_t adj_z = hostSplitFPSum(perturbations[(3 * np) + 2],
                                             little_psw.zcrd[target_atoms[j]],
                                             little_psw.zcrd_ovrf[target_atoms[j]]);
        little_psw.xcrd[target_atoms[j]]      = adj_x.x;
        little_psw.xcrd_ovrf[target_atoms[j]] = adj_x.y;
        little_psw.ycrd[target_atoms[j]]      = adj_y.x;
        little_psw.ycrd_ovrf[target_atoms[j]] = adj_y.y;
        little_psw.zcrd[target_atoms[j]]      = adj_z.x;
        little_psw.zcrd_ovrf[target_atoms[j]] = adj_z.y;
      }
      CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> cg(little_systems, little_topologies,
                                                0.5 * default_pme_cutoff, 0.1, 4, nbkinds[i]);
      const CellGridReader cgr = cg.data();
      const LocalExclusionMask little_lema(little_topologies, nbkinds[i]);
      ScoreCard sc(little_nsys);
      evaluateParticleParticleEnergy<Tcoord, Tacc,
                                     Tcalc, Tcoord4>(&cg, &sc, little_lema, default_pme_cutoff,
                                                     default_pme_cutoff, ew_coeff, ew_coeff,
                                                     VdwSumMethod::CUTOFF, EvaluateForce::YES,
                                                     nbkinds[i]);
      for (int j = 0; j < little_nsys; j++) {
        
        // Compute the base energy. Record the force on the chosen particle.
        const double total_nrg = sc.reportInstantaneousStates(StateVariable::ELECTROSTATIC, j) +
                                 sc.reportInstantaneousStates(StateVariable::VDW, j);
        const uint target_img_idx = cgr.img_atom_idx[target_atoms[j]];
        if (np == 0) {
          chosen_particle_forces[(3 * j)    ] = hostInt95ToDouble(cgr.xfrc[target_img_idx],
                                                                  cgr.xfrc_ovrf[target_img_idx]) *
                                                cgr.inv_frc_scale;
          chosen_particle_forces[(3 * j) + 1] = hostInt95ToDouble(cgr.yfrc[target_img_idx],
                                                                  cgr.yfrc_ovrf[target_img_idx]) *
                                                cgr.inv_frc_scale;
          chosen_particle_forces[(3 * j) + 2] = hostInt95ToDouble(cgr.zfrc[target_img_idx],
                                                                  cgr.zfrc_ovrf[target_img_idx]) *
                                                cgr.inv_frc_scale;
        }
        chosen_particle_fdnrg[j][np] = total_nrg;
      }
    }

    // Make the final perturbation to return the chosen particles to their original positions.
    for (int j = 0; j < little_nsys; j++) {
      const int95_t adj_x = hostSplitFPSum(perturbations[21], little_psw.xcrd[target_atoms[j]],
                                           little_psw.xcrd_ovrf[target_atoms[j]]);
      const int95_t adj_y = hostSplitFPSum(perturbations[22], little_psw.ycrd[target_atoms[j]],
                                           little_psw.ycrd_ovrf[target_atoms[j]]);
      const int95_t adj_z = hostSplitFPSum(perturbations[23], little_psw.zcrd[target_atoms[j]],
                                           little_psw.zcrd_ovrf[target_atoms[j]]);
      little_psw.xcrd[target_atoms[j]]      = adj_x.x;
      little_psw.xcrd_ovrf[target_atoms[j]] = adj_x.y;
      little_psw.ycrd[target_atoms[j]]      = adj_y.x;
      little_psw.ycrd_ovrf[target_atoms[j]] = adj_y.y;
      little_psw.zcrd[target_atoms[j]]      = adj_z.x;
      little_psw.zcrd_ovrf[target_atoms[j]] = adj_z.y;
    }

    // Measure the finite difference forces on each chosen atom.
    for (int j = 0; j < little_nsys; j++) {
      chosen_particle_fdfrc[(3 * j)    ] = (chosen_particle_fdnrg[j][2] -
                                            chosen_particle_fdnrg[j][1]) / pos_incr;
      chosen_particle_fdfrc[(3 * j) + 1] = (chosen_particle_fdnrg[j][4] -
                                            chosen_particle_fdnrg[j][3]) / pos_incr;
      chosen_particle_fdfrc[(3 * j) + 2] = (chosen_particle_fdnrg[j][6] -
                                            chosen_particle_fdnrg[j][5]) / pos_incr;
    }
    const double test_margin = (nbkinds[i] == NonbondedTheme::ALL) ? 1.0e-4 : 7.5e-5;
    check(chosen_particle_forces, RelationalOperator::EQUAL,
          Approx(chosen_particle_fdfrc).margin(test_margin), "The non-bonded forces computed for "
          "selected particles from various small systems in periodic systems do not match "
          "finite-difference calculations when the CellGrid is of theme " +
          getEnumerationName(nbkinds[i]) + ".", tsm.getTestingStatus());
  }
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
// Run a single test of the interaction kernel suited for particular conditions.
//
// Arguments:
//   poly_ps:   The collection of coordinate sets, used to store forces after neighbor list-based
//              GPU calculations and to run CPU-based calculations
//   mmctrl:    Molecular mechanics control information and task counters
//   sc:        Energy tracking object
//   poly_ag:   The synthesis of collated topologies, used to create neighbor list cell grids
//   lem:       An expression of each particle's excluded interactions based on relative order in
//              the topology
//   nrg_tab:   Tablulated cubic splines for elecrostatic energy and force interpolation
//   launcher:  Holds launch parameters for various pair interaction kernels
//   elec_cut:  The cutoff for electrostatic particle-particle interactions
//   vdw_cut:   The cutoff for van-der Waals particle-particle interactions
//   coord_v:   Precision in which coordinates are stored and force accumulated
//   calc_v:    Precision in which arithmetic calculations are performed
//   ngbr_v:    Indicate separate neighbor lists for electrostatics and van-der Waals interactions
//              should be constructed, or if a single merged neighbor list is requested
//   force_v:   Indicate whether to compute forces between atoms
//   energy_v:  Indicate whether to accumulate the non-bonded energy due to pairwise interactions
//   clash_v:   Indicate whether to employ clash mitigation
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcoord4>
std::vector<double> pairIKT(PhaseSpaceSynthesis *poly_ps, MolecularMechanicsControls *mmctrl,
                            ScoreCard *sc, const AtomGraphSynthesis &poly_ag,
                            const LocalExclusionMask &lem, const PPITable &nrg_tab,
                            const CoreKlManager &launcher, const double elec_cut,
                            const double vdw_cut, const PrecisionModel coord_v,
                            const PrecisionModel calc_v, const NeighborListKind ngbr_v,
                            const EvaluateForce force_v, const EvaluateEnergy energy_v,
                            const ClashResponse clash_v) {

  // Prime the work unit counters.  The essential action is to set the non-bondedd work unit
  // counter initially to the number of warps in the entire launch grid, from which the tower-plate
  // kernel will increment it all the way to the total number of neighbor list cells (in the
  // combined grid, or in both electrostatic and van-der Waals grids) during its asynchronous work.
  size_t coord_type_index;
  switch (coord_v) {
  case PrecisionModel::DOUBLE:
    coord_type_index = double_type_index;
    break;
  case PrecisionModel::SINGLE:
    coord_type_index = float_type_index;
    break;
  }

  // Pull the GPU specifications out of the launcher (reduces the list of input arguments)
  const GpuDetails &gpu = launcher.getGpu();
  
  // Prepare the coordinate synthesis to receive computed forces
  PsSynthesisWriter poly_psw = poly_ps->data(HybridTargetLevel::DEVICE);
  psyInitializeForces(&poly_psw, -1, gpu);
  std::vector<double> result;
  sc->initialize(HybridTargetLevel::DEVICE);
  switch (ngbr_v) {
  case NeighborListKind::DUAL:
    {
      CellGrid<Tcoord, Tacc, Tcoord, Tcoord4> cg_qq(poly_ps, poly_ag, 0.5 * elec_cut, 0.1, 4,
                                                    NonbondedTheme::ELECTROSTATIC);
      CellGrid<Tcoord, Tacc, Tcoord, Tcoord4> cg_lj(poly_ps, poly_ag, 0.5 * vdw_cut, 0.1, 4,
                                                    NonbondedTheme::VAN_DER_WAALS);
      const TinyBoxPresence has_tiny_box = (cg_qq.getTinyBoxPresence() == TinyBoxPresence::YES ||
                                            cg_lj.getTinyBoxPresence() == TinyBoxPresence::YES) ?
                                           TinyBoxPresence::YES : TinyBoxPresence::NO;
      TileManager tlmn(launcher.getPMEPairsKernelDims(coord_v, calc_v, ngbr_v, has_tiny_box,
                                                      force_v, energy_v, clash_v));
      mmctrl->primeWorkUnitCounters(launcher, force_v, energy_v, clash_v, VwuGoal::MOVE_PARTICLES,
                                    PrecisionModel::SINGLE, calc_v, QMapMethod::ACC_SHARED,
                                    PrecisionModel::SINGLE, coord_type_index, 5, ngbr_v,
                                    has_tiny_box, poly_ag);
      mmctrl->upload();

      // Compute the pair interactions on the GPU
      cg_qq.upload();
      cg_lj.upload();
      launchPMEPairs(calc_v, lem, nrg_tab, &cg_qq, &cg_lj, &tlmn, sc, mmctrl, force_v, energy_v,
                     launcher);

      // Transmit the forces from the cell grid to the coordinate synthesis
      cg_qq.contributeForces(HybridTargetLevel::DEVICE, gpu);
      cg_lj.contributeForces(HybridTargetLevel::DEVICE, gpu);
    }
    break;
  case NeighborListKind::MONO:
    {
      CellGrid<Tcoord, Tacc, Tcoord, Tcoord4> cg(poly_ps, poly_ag, 0.5 * vdw_cut, 0.1, 4,
                                                 NonbondedTheme::ALL);
      TileManager tlmn(launcher.getPMEPairsKernelDims(coord_v, calc_v, ngbr_v,
                                                      cg.getTinyBoxPresence(), force_v, energy_v,
                                                      clash_v));
      mmctrl->primeWorkUnitCounters(launcher, force_v, energy_v, clash_v, VwuGoal::MOVE_PARTICLES,
                                    PrecisionModel::SINGLE, calc_v, QMapMethod::ACC_SHARED,
                                    PrecisionModel::SINGLE, coord_type_index, 5, ngbr_v,
                                    cg.getTinyBoxPresence(), poly_ag);
      mmctrl->upload();
      
      // Compute the pair interactions on the GPU
      cg.upload();
      launchPMEPairs(calc_v, lem, nrg_tab, &cg, &tlmn, sc, mmctrl, force_v, energy_v, launcher);

      // Transmit the forces from the cell grid to the coordinate synthesis
      cg.contributeForces(HybridTargetLevel::DEVICE, gpu);
    }
    break;
  }
  
  // Extract the forces for all systems and return the result
  for (int i = 0; i < poly_ps->getSystemCount(); i++) {
    const std::vector<double> xyz_i = poly_ps->getInterlacedCoordinates(i, TrajectoryKind::FORCES,
                                                                        HybridTargetLevel::DEVICE);
    result.insert(result.end(), xyz_i.begin(), xyz_i.end());
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Execute one test of a specific pair interaction kernel.
//
// Arguments:
//   gpu:         Details of the GPU to use
//   nrg_tab:     Table of interactions for electrostatic energy and forces.  Also contains the
//                cutoff to be used by the calculations.
//   dyncon:      Dynamics control input, containing numbers of steps as well as cutoffs.  The
//                dyncon object will be checked for consistency with the electrostatic cutoff in
//                the nrg_tab object.
//   atom_limit:  The extreme number of atoms for systems to select from tsm. 
//   filter:      Instruction as to the limiting behavior of atom_limit.  Default "equal" (systems
//                will be selected from tsm if they have exactly atom_limit particles).
//   calc_v:      The calculation precision model to employ
//   coord_v:     The coordinate precision model to employ
//   ngbr_v:      The neighbor list layout to employ
//   clash_v:     Clash mitigation in effect
//   force_v:     Indicate whether force calculations are requested
//   energy_v:    Indicate whether energy calculations are requested
//   force_err:   Error tolerance for force calculations
//   energy_err:  Error tolerance for energy calculations
//   do_tests:    Flag to ensure that testing is possible
//-------------------------------------------------------------------------------------------------
void pairInteractionKernelTest(PhaseSpaceSynthesis *poly_ps, MolecularMechanicsControls *mmctrl,
                               ScoreCard *sc, const AtomGraphSynthesis &poly_ag,
                               const LocalExclusionMask &lem, const CoreKlManager &launcher,
                               const PPITable &nrg_tab, const PrecisionModel calc_v,
                               const PrecisionModel coord_v, const NeighborListKind ngbr_v,
                               const ClashResponse clash_v, const EvaluateForce force_v,
                               const EvaluateEnergy energy_v, const double force_err,
                               const double energy_err, const TestPriority do_tests) {
  const double elec_cut = mmctrl->getElectrostaticCutoff();
  const double vdw_cut = mmctrl->getVanDerWaalsCutoff();
  if (fabs(elec_cut - nrg_tab.getCutoff()) > 1.0e-6) {
    rtErr("The electrostatic cutoffs obtained from the energy table (" +
          realToString(nrg_tab.getCutoff(), 9, 4, NumberFormat::STANDARD_REAL) + ") and the "
          "dynamics controls (" + realToString(elec_cut, 9, 4, NumberFormat::STANDARD_REAL) +
          ") do not agree.", "pairInteractionKernelTest");
  }
  
  // Branch over the neighbor list composition
  std::vector<double> gpu_frc;
  switch (calc_v) {
  case PrecisionModel::DOUBLE:
    switch (coord_v) {
    case PrecisionModel::DOUBLE:
      gpu_frc = pairIKT<double, llint, double4>(poly_ps, mmctrl, sc, poly_ag, lem, nrg_tab,
                                                launcher, elec_cut, vdw_cut, coord_v, calc_v,
                                                ngbr_v, force_v, energy_v, clash_v);
      break;
    case PrecisionModel::SINGLE:
      gpu_frc = pairIKT<float, int, float4>(poly_ps, mmctrl, sc, poly_ag, lem, nrg_tab, launcher,
                                            elec_cut, vdw_cut, coord_v, calc_v, ngbr_v, force_v,
                                            energy_v, clash_v);
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (coord_v) {
    case PrecisionModel::DOUBLE:
      gpu_frc = pairIKT<double, llint, double4>(poly_ps, mmctrl, sc, poly_ag, lem, nrg_tab,
                                                launcher, elec_cut, vdw_cut, coord_v, calc_v,
                                                ngbr_v, force_v, energy_v, clash_v);
      break;
    case PrecisionModel::SINGLE:
      gpu_frc = pairIKT<float, int, float4>(poly_ps, mmctrl, sc, poly_ag, lem, nrg_tab, launcher,
                                            elec_cut, vdw_cut, coord_v, calc_v, ngbr_v, force_v,
                                            energy_v, clash_v);
      break;
    }
    break;
  }

  // Assemble the vectors of energies from the GPU
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  std::vector<double> gpu_qq_nrg = sc->reportInstantaneousStates(StateVariable::ELECTROSTATIC,
                                                                 devc);
  std::vector<double> gpu_lj_nrg = sc->reportInstantaneousStates(StateVariable::VDW, devc);
  
  // Test the force and energy computation kernels
  std::vector<double> cpu_frc;
  const PsSynthesisWriter poly_psw = poly_ps->data();
  std::vector<double2> nb_nrg(poly_psw.system_count);
  std::vector<double> cpu_qq_nrg(poly_psw.system_count), cpu_lj_nrg(poly_psw.system_count);
  for (int i = 0; i < poly_psw.system_count; i++) {
    PhaseSpace ps_i = poly_ps->exportSystem(i);
    const AtomGraph *ag_i = poly_ps->getSystemTopologyPointer(i);
    const LocalExclusionMask lema_i(ag_i);
    nb_nrg[i] = evaluateParticleParticleEnergy(&ps_i, ag_i, lema_i, PrecisionModel::DOUBLE,
                                               elec_cut, vdw_cut, nrg_tab.getEwaldCoefficient());
    const std::vector<double> xyz_i = ps_i.getInterlacedCoordinates(TrajectoryKind::FORCES);
    cpu_frc.insert(cpu_frc.end(), xyz_i.begin(), xyz_i.end());
    cpu_qq_nrg[i] = nb_nrg[i].x;
    cpu_lj_nrg[i] = nb_nrg[i].y;
  }
  
  // Find errant atoms and their system indices.  The vectors of forces do not present a simple
  // way to understand exactly where the erroneous forces lie.
  int rslt_llim = 0;
  std::vector<double3> errant_atom_list;
  for (int i = 0; i < poly_psw.system_count; i++) {
    const int rslt_hlim = rslt_llim + poly_psw.atom_counts[i];
    for (int j = rslt_llim; j < rslt_hlim; j++) {
      const double dfx = cpu_frc[(3 * j)    ] - gpu_frc[(3 * j)    ];
      const double dfy = cpu_frc[(3 * j) + 1] - gpu_frc[(3 * j) + 1];
      const double dfz = cpu_frc[(3 * j) + 2] - gpu_frc[(3 * j) + 2];
      if (fabs(dfx) > force_err || fabs(dfy) > force_err || fabs(dfz) > force_err) {
        const double d_offset = poly_psw.atom_starts[i];
        errant_atom_list.push_back({ static_cast<double>(i),
                                     d_offset + static_cast<double>(j - rslt_llim),
                                     sqrt((dfx * dfx) + (dfy * dfy) + (dfz * dfz)) });
      }
    }
    rslt_llim = rslt_hlim;
  }
  std::sort(errant_atom_list.begin(), errant_atom_list.end(),
            []( double3 a, double3 b ) { return a.z > b.z; });
  const int nreport = std::min(static_cast<int>(errant_atom_list.size()), 8);
  std::string errant_atoms("Examples of large force discrepancies include: ");
  int curr_sys_idx = -1;
  for (int i = 0; i < nreport; i++) {
    if (static_cast<int>(errant_atom_list[i].x) != curr_sys_idx) {
      curr_sys_idx = errant_atom_list[i].x;
      if (i > 0) {
        errant_atoms.append(", ");
      }
      errant_atoms.append("(System " + std::to_string(static_cast<int>(curr_sys_idx)) + ") ");
    }
    errant_atoms += std::to_string(static_cast<int>(errant_atom_list[i].y)) + " (";
    errant_atoms += std::to_string(static_cast<int>(errant_atom_list[i].y) -
                                   poly_psw.atom_starts[static_cast<int>(errant_atom_list[i].x)]) +
                    ") ";
  }
  errant_atoms.append(". (Atom indices are indicated by the system index of the synthesis, "
                      "the index of the atom in the entire synthesis, followed by the index of "
                      "the atom within its own topology in parentheses.)");
  std::string sys_desc("System size");
  if (poly_psw.system_count == 1) {
    sys_desc += ": " + std::to_string(poly_psw.atom_counts[0]) + " particles.  ";
  }
  else {
    const int highest_natom = maxValue(poly_psw.atom_counts, poly_psw.system_count);
    const int lowest_natom  = minValue(poly_psw.atom_counts, poly_psw.system_count);
    sys_desc += "s: (" + std::to_string(poly_psw.system_count) + " total) " +
    std::to_string(lowest_natom) + " - " + std::to_string(highest_natom) + " particles.  ";
  }
  switch (ngbr_v) {
  case NeighborListKind::DUAL:
    sys_desc += "Cutoffs " + realToString(elec_cut, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
      realToString(vdw_cut, 7, 4, NumberFormat::STANDARD_REAL);
    break;
  case NeighborListKind::MONO:
    sys_desc += "Cutoff " + realToString(vdw_cut, 7, 4, NumberFormat::STANDARD_REAL);
    break;
  }
  const std::string intr_desc = "Interactions were computed using " + getEnumerationName(coord_v) +
                                "-precision coordinates, " + getEnumerationName(calc_v) +
                                "-precision arithmetic, a " + getEnumerationName(ngbr_v) +
                                "-type neighbor list, clash mitigation " +
                                getEnumerationName(clash_v) + ", and " +
                                (energy_v == EvaluateEnergy::YES ? "energy evaluation" :
                                                                   "no energy evaluation") + ".  ";
  switch (force_v) {
  case EvaluateForce::YES:
    check(gpu_frc, RelationalOperator::EQUAL, Approx(cpu_frc).margin(force_err), "Forces do "
          "not agree with a CPU-based reference calculation.  " + intr_desc + sys_desc + ".  " +
          errant_atoms, do_tests);
    break;
  case EvaluateForce::NO:
    break;
  }
  switch (energy_v) {
  case EvaluateEnergy::YES:
    check(gpu_qq_nrg, RelationalOperator::EQUAL, Approx(cpu_qq_nrg).margin(energy_err),
          "Electrostatic energies evaluated by the GPU do not agree with those evaluated by "
          "independent CPU methods.  " + intr_desc + sys_desc + ".", do_tests);
    check(gpu_lj_nrg, RelationalOperator::EQUAL, Approx(cpu_lj_nrg).margin(energy_err),
          "Lennard-Jones energies evaluated by the GPU do not agree with those evaluated by "
          "independent CPU methods.  " + intr_desc + sys_desc + ".", do_tests);
    break;
  case EvaluateEnergy::NO:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
// Execute combinatorial tests of pair interaction kernels to cover different situations.
// Descriptions of input parameters follow from pairInteractionKernelTest(), above, with these
// modifications:
//
// Overloaded:
//   - Provide a collection of test systems in order to construct a synthesis
//   - Provide the synthesis directly
//
// Arguments:
//   tsm:           Collection of test systems (this is used to create the synthesis for testing)
//   error_tol:     The error tolerance for force ("x" member of the tuple) and energy
//                  calculations ("y" member of the tuple)
//   nt_warp_mult:  The multiplicity for evaluate neutral-territory tower / plate interactions
//                  (this does not apply to tower / tower interactions, for which each cell's
//                  assignment in any neighbor list is completed by one and only one warp)
//   calc_x:        Array of calculation precision models to employ.  In this and other input
//                  lists, an empty array triggers default behavior as can be inspected inside the
//                  function.  Here, the default is a solitary SINGLE precision iteraction of the
//                  kernel.
//   coord_x:       Array of coordinate precision models to employ.  Default SINGLE precision.
//   ngbr_x:        Array of neighbor list layouts.  Default MONO.
//   clash_x:       Array of clash handling scenarios.  Default NONE (no clash mitigation).
//   force_x:       Array of flags to trigger force calculation.  Default YES.
//   energy_x:      Array of flags to trigger energy calculation.  Default YES.
//-------------------------------------------------------------------------------------------------
void pairInteractionKernelLoop(PhaseSpaceSynthesis *poly_ps, AtomGraphSynthesis *poly_ag,
                               const GpuDetails &gpu, const PPITable &nrg_tab,
                               const DynamicsControls &dyncon,
                               const std::vector<double2> &error_tol, const TestPriority do_tests,
                               const int nt_warp_mult = 1,
                               const std::vector<PrecisionModel> calc_x = {},
                               const std::vector<PrecisionModel> coord_x = {},
                               const std::vector<NeighborListKind> ngbr_x = {},
                               const std::vector<ClashResponse> clash_x = {},
                               const std::vector<EvaluateForce> force_x = {},
                               const std::vector<EvaluateEnergy> energy_x = {}) {

  // Construct other essential resources for the test
  LocalExclusionMask lem(poly_ag, NonbondedTheme::ALL);
  const CoreKlManager launcher(gpu, poly_ag);
  poly_ps->upload();
  poly_ag->upload();
  lem.upload();
  ScoreCard sc(poly_ps->getSystemCount(), 1, 32);
  int largest_system = 0;
  for (int i = 0; i < poly_ps->getSystemCount(); i++) {
    largest_system = std::max(poly_ps->getAtomCount(i), largest_system);
  }
  MolecularMechanicsControls mmctrl(dyncon);
  mmctrl.setNTWarpMultiplicity(nt_warp_mult);
  
  // Set the arrays of testing conditions
  const std::vector<PrecisionModel> calc_exec = (calc_x.size() == 0) ?
    std::vector<PrecisionModel>(1, PrecisionModel::SINGLE) : calc_x;
  const std::vector<PrecisionModel> coord_exec = (coord_x.size() == 0) ?
    std::vector<PrecisionModel>(1, PrecisionModel::SINGLE) : coord_x;
  const std::vector<NeighborListKind> ngbr_exec = (ngbr_x.size() == 0) ?
    std::vector<NeighborListKind>(1, NeighborListKind::MONO) : ngbr_x;
  const std::vector<ClashResponse> clash_exec = (clash_x.size() == 0) ?
    std::vector<ClashResponse>(1, ClashResponse::NONE) : clash_x;
  const std::vector<EvaluateForce> force_exec = (force_x.size() == 0) ?
    std::vector<EvaluateForce>(1, EvaluateForce::YES) : force_x;
  const std::vector<EvaluateEnergy> energy_exec = (energy_x.size() == 0) ?
    std::vector<EvaluateEnergy>(1, EvaluateEnergy::YES) : energy_x;

  // Loop over all enumerated possibilities
  int test_no = 0;
  for (size_t calc_idx = 0; calc_idx < calc_exec.size(); calc_idx++) {
    for (size_t coord_idx = 0; coord_idx < coord_exec.size(); coord_idx++) {
      for (size_t ngbr_idx = 0; ngbr_idx < ngbr_exec.size(); ngbr_idx++) {
        for (size_t clash_idx = 0; clash_idx < clash_exec.size(); clash_idx++) {
          for (size_t force_idx = 0; force_idx < force_exec.size(); force_idx++) {
            for (size_t energy_idx = 0; energy_idx < energy_exec.size(); energy_idx++) {

              // Skip cases where neither force nor energy are to be evaluated.  This case is not
              // filtered out by the production routines.  Also skip testing of system where the
              // energy alone is evaluated and the largest system is very large.  
              if (force_exec[force_idx] == EvaluateForce::NO) {
                if (energy_exec[energy_idx] == EvaluateEnergy::NO || largest_system > 3105) {
                  continue;
                }
              }

              // Skip large systems in redundant contexts
              if (poly_ps->getPaddedAtomCount() >= 10000 &&
                  (calc_exec[calc_idx] == PrecisionModel::DOUBLE ||
                   coord_exec[coord_idx] == PrecisionModel::DOUBLE) &&
                  energy_exec[energy_idx] == EvaluateEnergy::NO) {
                continue;
              }

              // Perform the test
              pairInteractionKernelTest(poly_ps, &mmctrl, &sc, *poly_ag, lem, launcher, nrg_tab,
                                        calc_exec[calc_idx], coord_exec[coord_idx],
                                        ngbr_exec[ngbr_idx], clash_exec[clash_idx],
                                        force_exec[force_idx], energy_exec[energy_idx],
                                        error_tol[test_no].x, error_tol[test_no].y, do_tests);
            }
          }
        }
      }
      test_no++;
    }
  }
}

void pairInteractionKernelLoop(const TestSystemManager &tsm, const GpuDetails &gpu,
                               const PPITable &nrg_tab, const DynamicsControls &dyncon,
                               const int atom_limit,
                               const std::vector<double2> &error_tol,
                               const RelationalOperator filter = RelationalOperator::EQ,
                               const int nt_warp_mult = 1,
                               const std::vector<PrecisionModel> calc_x = {},
                               const std::vector<PrecisionModel> coord_x = {},
                               const std::vector<NeighborListKind> ngbr_x = {},
                               const std::vector<ClashResponse> clash_x = {},
                               const std::vector<EvaluateForce> force_x = {},
                               const std::vector<EvaluateEnergy> energy_x = {}) {

  // Build the synthesis of systems outside the loop over combinations of run conditions
  const std::vector<int> systems_idx = tsm.getQualifyingSystems(atom_limit, filter);
  PhaseSpaceSynthesis poly_ps = tsm.exportPhaseSpaceSynthesis(systems_idx);
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(systems_idx);
  pairInteractionKernelLoop(&poly_ps, &poly_ag, gpu, nrg_tab, dyncon, error_tol,
                            tsm.getTestingStatus(), nt_warp_mult, calc_x, coord_x, ngbr_x,
                            clash_x, force_x, energy_x);
}

//-------------------------------------------------------------------------------------------------
// Place an ion particle in the designated neighbor list cell and cubelet.
//
// Arguments:
//   poly_psw:   The synthesis of coordinates
//   sys_idx:    Index of the system in which to reposition the particle
//   invu:       Vectorized copy of the transformation matrix taking fractional coordinates into
//               real space
//   cell_aidx:  Index of the neighbor list cell of interest, along the unit cell A axis
//   cell_bidx:  Index of the neighbor list cell of interest, along the unit cell B axis
//   cell_cidx:  Index of the neighbor list cell of interest, along the unit cell C axis
//   cblt_aidx:  Index of the cubelet within the subcell, along the unit cell A axis
//   cblt_bidx:  Index of the cubelet within the subcell, along the unit cell B axis
//   cblt_cidx:  Index of the cubelet within the subcell, along the unit cell C axis
//   nplaced:    The number of atoms placed in the system thus far
//   xrs:        Source of random numbers for placing a particle in the selected bin
//-------------------------------------------------------------------------------------------------
void placeIon(PsSynthesisWriter *poly_psw, const int sys_idx, const std::vector<double> &invu,
              const int cell_aidx, const int cell_bidx, const int cell_cidx, const int cblt_aidx,
              const int cblt_bidx, const int cblt_cidx, const int nplaced,
              Xoshiro256ppGenerator *xrs) {
  const double xupt = 0.2 * (static_cast<double>(cell_aidx) +
                             (0.25 * static_cast<double>(cblt_aidx)) +
                             ((xrs->uniformRandomNumber() - 0.5) * 0.0625));
  const double yupt = 0.2 * (static_cast<double>(cell_bidx) +
                             (0.25 * static_cast<double>(cblt_bidx)) +
                             ((xrs->uniformRandomNumber() - 0.5) * 0.0625));
  const double zupt = 0.2 * (static_cast<double>(cell_cidx) +
                             (0.25 * static_cast<double>(cblt_cidx)) +
                             ((xrs->uniformRandomNumber() - 0.5) * 0.0625));
  const double xpt = (invu[0] * xupt) + (invu[3] * yupt) + (invu[6] * zupt);
  const double ypt = (invu[1] * xupt) + (invu[4] * yupt) + (invu[7] * zupt);
  const double zpt = (invu[2] * xupt) + (invu[5] * yupt) + (invu[8] * zupt);
  const int95_t ixpt = hostDoubleToInt95(xpt * poly_psw->gpos_scale_f);
  const int95_t iypt = hostDoubleToInt95(ypt * poly_psw->gpos_scale_f);
  const int95_t izpt = hostDoubleToInt95(zpt * poly_psw->gpos_scale_f);
  const int atom_offset = poly_psw->atom_starts[sys_idx];
  poly_psw->xcrd[atom_offset + nplaced] = ixpt.x;
  poly_psw->ycrd[atom_offset + nplaced] = iypt.x;
  poly_psw->zcrd[atom_offset + nplaced] = izpt.x;
  poly_psw->xcrd_ovrf[atom_offset + nplaced] = ixpt.y;
  poly_psw->ycrd_ovrf[atom_offset + nplaced] = iypt.y;
  poly_psw->zcrd_ovrf[atom_offset + nplaced] = izpt.y;
}

//-------------------------------------------------------------------------------------------------
// Arrange the ions in a small box, given a series of instructions as to how many ions to place in
// any given segment of the box and placing others at random.  Each ion will be placed in the
// center of some cubelet in its neighbor list cell, up to 64 ions per cell.  In this way, a
// minimal distance will be guaranteed between any given pair of ions, even though there will still
// be some degree of randomness to prevent a regular formation from masking a mistake in the
// energetics.
//
// Arguments:
//   poly_ps:  The synthesis of coordinates
//   sys_idx:  Index of the system in which to reposition all particles
//   fill:     Series of instructions on how to fill cells with the available particles.  The
//             neighbor list cell a, b, and c positions are given in the "x", "y", and "z" members
//             of the tuple.  The number of particles to place in it is given by the "w" member.
//   xrs:      Source of random numbers for selecting bins and placing particles
//-------------------------------------------------------------------------------------------------
void arrangeIons(PhaseSpaceSynthesis *poly_ps, const int sys_idx, const std::vector<int4> &fill,
                 Xoshiro256ppGenerator *xrs) {
  PsSynthesisWriter poly_psw = poly_ps->data();
  if (sys_idx < 0 || sys_idx >= poly_psw.system_count) {
    rtErr("System index " + std::to_string(sys_idx) + " is invalid for a synthesis of " +
          std::to_string(poly_psw.system_count) + " systems.", "arrangeIons");
  }
  const int natom = poly_psw.atom_counts[sys_idx];
  const int nsubdiv = std::max(5, static_cast<int>(ceil(cbrt(static_cast<double>(natom) / 64.0))));
  const int atom_offset = poly_psw.atom_starts[sys_idx];
  const size_t ninsr = fill.size();
  std::vector<bool> occupied(nsubdiv * nsubdiv * nsubdiv * 64, false); 
  int nplaced = 0;
  std::vector<double> invu(9);
  const int xfrm_stride = roundUp(9, warp_size_int);
  for (int i = 0; i < 9; i++) {
    invu[i] = poly_psw.invu[(sys_idx * xfrm_stride) + i];
  }
  for (int i = 0; i < ninsr; i++) {
    
    // Identify the cell into which the particle must be placed.
    const int cell_aidx = (fill[i].x >= 0) ?
                          fill[i].x % nsubdiv : fill[i].x + (roundUp(-fill[i].x, nsubdiv));
    const int cell_bidx = (fill[i].y >= 0) ?
                          fill[i].y % nsubdiv : fill[i].y + (roundUp(-fill[i].y, nsubdiv));
    const int cell_cidx = (fill[i].z >= 0) ?
                          fill[i].z % nsubdiv : fill[i].z + (roundUp(-fill[i].z, nsubdiv));
    const int cell_idx = (((cell_cidx * nsubdiv) + cell_bidx) * nsubdiv) + cell_aidx;

    const int ni = std::min(std::min(fill[i].w, 64), natom - nplaced);
    for (int j = 0; j < ni; j++) {

      // Try placing the particle in a bin with random selection.  If that fails, place it in the
      // first available bin found during a linear search.
      int k = 0;
      bool placed = false;
      while (k < 8 && placed == false) {
        const int bin_a = static_cast<int>(xrs->uniformRandomNumber() * 4.0);
        const int bin_b = static_cast<int>(xrs->uniformRandomNumber() * 4.0);
        const int bin_c = static_cast<int>(xrs->uniformRandomNumber() * 4.0);
        const int occ_idx = (64 * cell_idx) + (((bin_c * 4) + bin_b) * 4) + bin_a;
        if (occupied[occ_idx] == false) {
          occupied[occ_idx] = true;
          placeIon(&poly_psw, sys_idx, invu, cell_aidx, cell_bidx, cell_cidx, bin_a, bin_b, bin_c,
                   nplaced, xrs);
          placed = true;
          nplaced++;
        }
        k++;
      }
      k = 0;
      while (k < 64 && placed == false) {
        const int occ_idx = (64 * cell_idx) + k;
        if (occupied[occ_idx] == false) {
          occupied[occ_idx] = true;
          const int bin_c = (k >> 4);
          const int bin_b = ((k - (bin_c * 16)) >> 2);
          const int bin_a = k - (((bin_c * 4) + bin_b) * 4);
          placeIon(&poly_psw, sys_idx, invu, cell_aidx, cell_bidx, cell_cidx, bin_a, bin_b,
                   bin_c, nplaced, xrs);
          placed = true;
          nplaced++;
        }
        k++;
      }
    }
  }
  int jcounter = 0;
  for (int i = nplaced; i < natom; i++) {
    int j = 0;
    bool placed = false;
    while (j < 16 && placed == false) {
      const int cell_idx = xrs->uniformRandomNumber() * 125.0;
      const int cell_cidx = (cell_idx / 25);
      const int cell_bidx = (cell_idx - (cell_cidx * 25)) / 5;
      const int cell_aidx = cell_idx - (((cell_cidx * 5) + cell_bidx) * 5);
      const int bin_a = static_cast<int>(xrs->uniformRandomNumber() * 4.0);
      const int bin_b = static_cast<int>(xrs->uniformRandomNumber() * 4.0);
      const int bin_c = static_cast<int>(xrs->uniformRandomNumber() * 4.0);
      const int occ_idx = (64 * cell_idx) + (((bin_c * 4) + bin_b) * 4) + bin_a;
      if (occupied[occ_idx] == false) {
        occupied[occ_idx] = true;
        placeIon(&poly_psw, sys_idx, invu, cell_aidx, cell_bidx, cell_cidx, bin_a, bin_b, bin_c,
                 nplaced, xrs);
        placed = true;
        nplaced++;
      }
      j++;
    }
    j = 0;
    while (j < 125 && placed == false) {
      const int cell_idx = (j + jcounter) % 125; 
      const int cell_cidx = (cell_idx / 25);
      const int cell_bidx = (cell_idx - (cell_cidx * 25)) / 5;
      const int cell_aidx = cell_idx - (((cell_cidx * 5) + cell_bidx) * 5);
      int k = 0;
      while (k < 64 && placed == false) {
        const int occ_idx = (64 * j) + k;
        if (occupied[occ_idx] == false) {
          occupied[occ_idx] = true;
          const int bin_c = (k >> 4);
          const int bin_b = ((k - (bin_c * 16)) >> 2);
          const int bin_a = k - (((bin_c * 4) + bin_b) * 4);
          placeIon(&poly_psw, sys_idx, invu, cell_aidx, cell_bidx, cell_cidx, bin_a, bin_b,
                   bin_c, nplaced, xrs);
          placed = true;
          jcounter++;
          nplaced++;
        }
        k++;
      }
      j++;
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Test the GPU functionality related to neighbor list force and energy calculations.
//
// Arguments:
//   tsm:      The collection of test systems
//   testdir:  Location of STORMM test materials, contains the path to STORMM's source code
//   xrs:      Source of random numbers for testing
//   gpu:      Details of the GPU to use
//-------------------------------------------------------------------------------------------------
void runGpuTests(const TestSystemManager &tsm, const std::string &testdir,
                 Xoshiro256ppGenerator *xrs, const GpuDetails &gpu) {

  // Create an exemplary namelist
  const std::string dyn_str("&dynamics\n  cut = " + realToString(default_pme_cutoff) + "\n&end\n");
  const TextFile dyn_tf(dyn_str, TextOrigin::RAM);
  int start_line = 0;
  bool found_nml = false;
  DynamicsControls dyncon(dyn_tf, &start_line, &found_nml);

  // Create the energy and force spline table
  PPITable nrg_tab(NonbondedTheme::ELECTROSTATIC, BasisFunctions::MIXED_FRACTIONS,
                   TableIndexing::SQUARED_ARG);
  const std::vector<PrecisionModel> all_prec = { PrecisionModel::DOUBLE, PrecisionModel::SINGLE };
  const std::vector<NeighborListKind> all_ngbr = { NeighborListKind::DUAL,
                                                   NeighborListKind::MONO };
  const std::vector<ClashResponse> all_clash = { ClashResponse::NONE };
  const std::vector<EvaluateEnergy> all_energy = { EvaluateEnergy::YES, EvaluateEnergy::NO };
  const std::vector<EvaluateForce> all_force = { EvaluateForce::YES, EvaluateForce::NO };

  // The tests will proceed with calculation precision as the outer loop variable in
  // pairInteractionKernelLoop(), followed by coordinate precision.  Arrange vectors of tolerances
  // accordingly.
  const std::vector<double2> brbz_tol = { { 3.0e-7, 1.0e-7 }, { 4.0e-7, 5.0e-7 },
                                          { 3.6e-7, 5.0e-7 }, { 5.0e-5, 1.0e-6 } };
  pairInteractionKernelLoop(tsm, gpu, nrg_tab, dyncon, 12, brbz_tol, RelationalOperator::EQ, 1,
                            all_prec, all_prec, all_ngbr, all_clash, all_force, all_energy);
  const std::vector<double2> drug_tol = { { 8.4e-7, 1.0e-6 }, { 6.0e-6, 5.0e-6 },
                                          { 6.0e-6, 1.1e-5 }, { 5.0e-5, 1.2e-5 } };
  pairInteractionKernelLoop(tsm, gpu, nrg_tab, dyncon, 53, drug_tol, RelationalOperator::EQ, 1,
                            all_prec, all_prec, all_ngbr, all_clash, all_force, all_energy);
  const std::vector<double2> ubiq_tol = { { 2.0e-6, 1.0e-6 }, { 4.3e-4, 3.0e-4 },
                                          { 3.0e-4, 1.6e-3 }, { 4.0e-4, 1.3e-3 } };
  pairInteractionKernelLoop(tsm, gpu, nrg_tab, dyncon, 3105, ubiq_tol, RelationalOperator::EQ, 1,
                            all_prec, all_prec, all_ngbr, all_clash, all_force, all_energy);
  pairInteractionKernelLoop(tsm, gpu, nrg_tab, dyncon, 12, brbz_tol, RelationalOperator::EQ, 2,
                            all_prec, all_prec, all_ngbr, all_clash, all_force, all_energy);
  pairInteractionKernelLoop(tsm, gpu, nrg_tab, dyncon, 53, drug_tol, RelationalOperator::EQ, 2,
                            all_prec, all_prec, all_ngbr, all_clash, all_force, all_energy);
  pairInteractionKernelLoop(tsm, gpu, nrg_tab, dyncon, 3105, ubiq_tol, RelationalOperator::EQ, 2,
                            all_prec, all_prec, all_ngbr, all_clash, all_force, all_energy);
  const std::vector<double2> tp4p_tol = { { 1.5e-6, 1.0e-6 }, { 4.5e-4, 3.5e-4 },
                                          { 3.9e-4, 2.7e-3 }, { 6.0e-4, 2.7e-3 } };
  pairInteractionKernelLoop(tsm, gpu, nrg_tab, dyncon, 1024, tp4p_tol, RelationalOperator::EQ, 2,
                            all_prec, all_prec, all_ngbr, all_clash, all_force, all_energy);

  // Check syntheses of multiple systems, including "tiny" and non-tiny boxes.
  const std::vector<double2> many_tol = { { 1.0e-6, 1.0e-6 }, { 8.0e-4, 3.8e-4 },
                                          { 8.0e-4, 2.2e-3 }, { 8.0e-3, 2.1e-3 } };
  pairInteractionKernelLoop(tsm, gpu, nrg_tab, dyncon, 768, many_tol,
                            RelationalOperator::LE, 1, all_prec, all_prec, all_ngbr, all_clash,
                            all_force, all_energy);
  const std::vector<double2> rnbw_tol = { { 1.6e-6, 2.1e-6 }, { 8.0e-4, 3.8e-4 },
                                          { 1.5e-3, 1.2e-2 }, { 8.0e-3, 1.2e-2 } };
  pairInteractionKernelLoop(tsm, gpu, nrg_tab, dyncon, 768, rnbw_tol, RelationalOperator::GE, 1,
                            all_prec, all_prec, all_ngbr, all_clash, all_force, all_energy);

  // Open the large kinase system and perform additional tests
  const std::vector<std::string> kinase_str(1, "kinase");
  const TestSystemManager kinase_tsm(testdir + "Topology", "top", kinase_str,
                                     testdir + "Trajectory", "inpcrd", kinase_str);
  const std::vector<double2> kins_tol = { { 2.5e-6, 7.0e-6 }, { 1.1e-3, 2.5e-2 },
                                          { 4.5e-3, 9.6e-2 }, { 4.5e-3, 9.0e-2 } };
  pairInteractionKernelLoop(kinase_tsm, gpu, nrg_tab, dyncon, 52889, kins_tol,
                            RelationalOperator::EQ, 2, all_prec, all_prec, all_ngbr, all_clash,
                            all_force, all_energy);
  dyncon.setCutoff(12.0);
  PPITable long_reach_nrg_tab(NonbondedTheme::ELECTROSTATIC, BasisFunctions::MIXED_FRACTIONS,
                              TableIndexing::SQUARED_ARG, dyncon.getElectrostaticCutoff());
  pairInteractionKernelLoop(kinase_tsm, gpu, long_reach_nrg_tab, dyncon, 52889, kins_tol,
                            RelationalOperator::EQ, 2, all_prec, all_prec, all_ngbr, all_clash,
                            all_force, all_energy);

  // Open the large DHFR (JAC benchmark) system and perform additional tests
  const std::vector<std::string> jac_str(1, "jac");
  const TestSystemManager jac_tsm(testdir + "Topology", "top", jac_str, testdir + "Trajectory",
                                  "inpcrd", jac_str);
  pairInteractionKernelLoop(jac_tsm, gpu, long_reach_nrg_tab, dyncon, 23558, kins_tol,
                            RelationalOperator::EQ, 2, all_prec, all_prec, all_ngbr, all_clash,
                            all_force, all_energy);
  dyncon.setCutoff(default_pme_cutoff);

  // Test a group of the protein-in-water systems
  const std::vector<int> ubiq_idx = tsm.getQualifyingSystems(3105, RelationalOperator::EQ);
  const std::vector<int> systems_vec(8, ubiq_idx[0]);
  PhaseSpaceSynthesis ubiq_poly_ps = tsm.exportPhaseSpaceSynthesis(systems_vec);
  PsSynthesisWriter ubiq_poly_psw = ubiq_poly_ps.data();
  for (int i = 0; i < ubiq_poly_ps.getSystemCount(); i++) {
    addRandomNoise(xrs, ubiq_poly_psw.xcrd, ubiq_poly_psw.xcrd_ovrf, ubiq_poly_psw.ycrd,
                   ubiq_poly_psw.ycrd_ovrf, ubiq_poly_psw.zcrd, ubiq_poly_psw.zcrd_ovrf,
                   ubiq_poly_ps.getAtomCount(i), 0.01, ubiq_poly_psw.gpos_scale_f);
  }
  AtomGraphSynthesis ubiq_poly_ag = tsm.exportAtomGraphSynthesis(systems_vec);
  pairInteractionKernelLoop(&ubiq_poly_ps, &ubiq_poly_ag, gpu, nrg_tab, dyncon, ubiq_tol,
                            tsm.getTestingStatus(), 2, all_prec, all_prec, all_ngbr, all_clash,
                            all_force, all_energy);

  // Iterative tests with 1000 ions in random arrangements
  const std::vector<std::string> thsnd(1, "thousand_ions");
  const TestSystemManager thsnd_tsm(testdir + "Topology", "top", thsnd, testdir + "Trajectory",
                                    "inpcrd", thsnd);
  const std::vector<int> zeros(50, 0);
  PhaseSpaceSynthesis poly_thousand_ps = thsnd_tsm.exportPhaseSpaceSynthesis(zeros);
  AtomGraphSynthesis poly_thousand_ag = thsnd_tsm.exportAtomGraphSynthesis(zeros);
  std::vector<int4> fill(17);
  const std::vector<int3> tower_plate = { {  0,  0,  0 }, {  0,  0, -2 }, {  0,  0, -1 },
                                          {  0,  0,  1 }, {  0,  0,  2 }, { -2, -2,  0 },
                                          { -1, -2,  0 }, {  0, -2,  0 }, {  1, -2,  0 },
                                          {  2, -2,  0 }, { -2, -1,  0 }, { -1, -1,  0 },
                                          {  0, -1,  0 }, {  1, -1,  0 }, {  2, -1,  0 },
                                          { -2,  0,  0 }, { -1,  0,  0 } };
  for (int i = 0; i < 5; i++) {

    // Try placing various numbers of atoms in the tower and plate cells.  Then, try placing the
    // entire arrangement to center at various positions along each of the principal unit cell
    // axes.  Next, try selecting neighbor list cells at random and placing various numbers of
    // particles in them.  This will generate some very uneven distributions that may detect less
    // obvious errors in the neighbor list evaluation.
    for (int j = 0; j < 17; j++) {
      fill[j].x = tower_plate[j].x + i;
      fill[j].y = tower_plate[j].y;
      fill[j].z = tower_plate[j].z;
      fill[j].w = xrs->uniformRandomNumber() * 16.0;
    }
    arrangeIons(&poly_thousand_ps, 3 * i, fill, xrs);
    for (int j = 0; j < 17; j++) {
      fill[j].x = tower_plate[j].x;
      fill[j].y = tower_plate[j].y + i;
      fill[j].z = tower_plate[j].z;
      fill[j].w = xrs->uniformRandomNumber() * 16.0;
    }
    arrangeIons(&poly_thousand_ps, (3 * i) + 1, fill, xrs);
    for (int j = 0; j < 17; j++) {
      fill[j].x = tower_plate[j].x;
      fill[j].y = tower_plate[j].y;
      fill[j].z = tower_plate[j].z + i;
      fill[j].w = xrs->uniformRandomNumber() * 16.0;
    }
    arrangeIons(&poly_thousand_ps, (3 * i) + 2, fill, xrs);
  }
  for (int i = 15; i < 50; i++) {
    for (int j = 0; j < 17; j++) {
      fill[j].x = xrs->uniformRandomNumber() * 5.0;
      fill[j].y = xrs->uniformRandomNumber() * 5.0;
      fill[j].z = xrs->uniformRandomNumber() * 5.0;
      fill[j].w = xrs->uniformRandomNumber() * 16.0;
    }
    arrangeIons(&poly_thousand_ps, i, fill, xrs);
  }
  const std::vector<double2> ions_tol = { { 7.2e-6, 2.9e-6 }, { 9.0e-3, 8.4e-3 },
                                          { 2.4e-3, 3.3e-3 }, { 7.2e-3, 8.5e-3 } };
  pairInteractionKernelLoop(&poly_thousand_ps, &poly_thousand_ag, gpu, nrg_tab, dyncon, ions_tol,
                            TestPriority::CRITICAL, 1, all_prec, all_prec, all_ngbr, all_clash,
                            all_force, all_energy);

  // Try a different cutoff for the van-der Waals interactions on selected systems
  dyncon.setVanDerWaalsCutoff(10.0);
  const std::vector<NeighborListKind> dual_ngbr = { NeighborListKind::DUAL };
  pairInteractionKernelLoop(tsm, gpu, nrg_tab, dyncon, 12, brbz_tol, RelationalOperator::EQ, 1,
                            all_prec, all_prec, dual_ngbr, all_clash, all_force, all_energy);
  pairInteractionKernelLoop(tsm, gpu, nrg_tab, dyncon, 53, drug_tol, RelationalOperator::EQ, 1,
                            all_prec, all_prec, dual_ngbr, all_clash, all_force, all_energy);
  pairInteractionKernelLoop(tsm, gpu, nrg_tab, dyncon, 3105, ubiq_tol, RelationalOperator::EQ, 1,
                            all_prec, all_prec, dual_ngbr, all_clash, all_force, all_energy);
  pairInteractionKernelLoop(&ubiq_poly_ps, &ubiq_poly_ag, gpu, nrg_tab, dyncon, ubiq_tol,
                            tsm.getTestingStatus(), 2, all_prec, all_prec, dual_ngbr, all_clash,
                            all_force, all_energy);
}
#endif

//-------------------------------------------------------------------------------------------------
// Perform additional tests of PME-related functions.
//-------------------------------------------------------------------------------------------------
void runMiscellaneousTests() {
  const double test_cutoff = 9.0;
  const double test_dsum_tol = 5.0e-6;
  const double ew_coeff = ewaldCoefficient(test_cutoff, test_dsum_tol);
  const double dsum_tol = recoverDirectSumTolerance(test_cutoff, ew_coeff);
  check(dsum_tol, RelationalOperator::EQUAL, Approx(test_dsum_tol).margin(1.0e-6), "The direct "
        "sum tolerance recovered from a computed Ewald coefficient and a known cutoff does not "
        "match the known tolerance used in computation of that same Ewald coefficient.");
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
  const int setup_walltime = timer.addCategory("System setup");
  const int qspr_walltime = timer.addCategory("Charge spreading tests");
  const int ppdir_walltime = timer.addCategory("Particle-particle tests");
  Xoshiro256ppGenerator xrs(oe.getRandomSeed());

#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> create_this_to_engage_gpu(1);
  const int gpu_walltime = timer.addCategory("GPU kernel testing");
#endif
  
  // Section 1                                                    
  section("Test mechanics in the PMIGrid object");

  // Section 2
  section("Test charge spreading in various modes");

  // Section 3
  section("Test simple particle-particle interactions");

  // Section 4
  section("Miscellaneous tests of PME-related utilities");
  
  // Create a synthesis of systems and the associated particle-mesh interaction grid
  section(1);
  const std::vector<std::string> pbc_systems = { "bromobenzene", "bromobenzene_vs", "tip3p",
                                                 "tip4p", "ubiquitin", "tamavidin", "drug_example",
                                                 "drug_example_dry", "trpcage_in_water" };
  const char osc = osSeparator();
  const std::string testdir = oe.getStormmSourcePath() + osc + "test" + osc;
  TestSystemManager tsm(testdir + "Topology", "top", pbc_systems, testdir + "Trajectory", "inpcrd",
                        pbc_systems);
  timer.assignTime(setup_walltime);
  
  // Create a synthesis of systems and the associated particle-mesh interaction grids
  const std::vector<int> psys = tsm.getQualifyingSystems({ UnitCellType::ORTHORHOMBIC,
                                                           UnitCellType::TRICLINIC });
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(psys);
  PhaseSpaceSynthesis poly_ps = tsm.exportPhaseSpaceSynthesis(psys);
#ifdef STORMM_USE_HPC
  const CoreKlManager launcher(gpu, poly_ag);
  poly_ag.upload();
  poly_ps.upload();
#else
  const CoreKlManager launcher(null_gpu, poly_ag);
#endif
  // Test the baseline, double / double case.  No fixed-precision accumulation.  The following
  // funciton toggles between test sections 1 and 2 internally.
  spatialDecompositionOuter(poly_ag, poly_ps, PrecisionModel::DOUBLE, false, 4,
                            NonbondedTheme::ELECTROSTATIC, PrecisionModel::DOUBLE, 0,
                            5, NonbondedTheme::ELECTROSTATIC, tsm.getTestingStatus(), 1.0e-8,
                            1.0e-7, launcher);

  // Fixed-precision accumulation, even in what would appear to be pretty fine increments, can
  // introduce a great deal of error.
  spatialDecompositionOuter(poly_ag, poly_ps, PrecisionModel::DOUBLE, false, 4,
                            NonbondedTheme::ELECTROSTATIC, PrecisionModel::DOUBLE, 32,
                            5, NonbondedTheme::ELECTROSTATIC, tsm.getTestingStatus(), 3.5e-5,
                            1.0e-6, launcher);
  spatialDecompositionOuter(poly_ag, poly_ps, PrecisionModel::DOUBLE, false, 4,
                            NonbondedTheme::ELECTROSTATIC, PrecisionModel::DOUBLE, 48,
                            5, NonbondedTheme::ELECTROSTATIC, tsm.getTestingStatus(), 1.0e-8,
                            1.0e-7, launcher);

  // Single-precision calculation of the B-splines introduces error of its own, but not as much as
  // a 32-bit fixed-precision model.  Restricting the precision to something that can be guaranteed
  // to fit in a single 32-bit accumulator is worse yet.  Spreading with lower interpolation order
  // is less damaging when combined with fixed-precision, as there are fewer small numbers to deal
  // with.
  spatialDecompositionOuter(poly_ag, poly_ps, PrecisionModel::SINGLE, false, 4,
                            NonbondedTheme::ELECTROSTATIC, PrecisionModel::SINGLE, 0,
                            5, NonbondedTheme::ELECTROSTATIC, tsm.getTestingStatus(), 9.0e-6,
                            1.0e-6, launcher);
  spatialDecompositionOuter(poly_ag, poly_ps, PrecisionModel::SINGLE, false, 4,
                            NonbondedTheme::ELECTROSTATIC, PrecisionModel::SINGLE, 28,
                            5, NonbondedTheme::ELECTROSTATIC, tsm.getTestingStatus(), 5.0e-4,
                            5.0e-5, launcher);
  spatialDecompositionOuter(poly_ag, poly_ps, PrecisionModel::SINGLE, false, 4,
                            NonbondedTheme::ELECTROSTATIC, PrecisionModel::SINGLE, 24,
                            4, NonbondedTheme::ELECTROSTATIC, tsm.getTestingStatus(), 4.4e-3,
                            5.0e-4, launcher);
  spatialDecompositionOuter(poly_ag, poly_ps, PrecisionModel::SINGLE, false, 4,
                            NonbondedTheme::ELECTROSTATIC, PrecisionModel::SINGLE, 28,
                            6, NonbondedTheme::ELECTROSTATIC, tsm.getTestingStatus(), 7.8e-4,
                            8.0e-5, launcher);

  // The precision model of the coordinate representation does not affect the total charge
  // conservation, but it can affect the representation of individual charges on the mesh.
  spatialDecompositionOuter(poly_ag, poly_ps, PrecisionModel::DOUBLE, false, 4,
                            NonbondedTheme::ELECTROSTATIC, PrecisionModel::SINGLE, 48,
                            5, NonbondedTheme::ELECTROSTATIC, tsm.getTestingStatus(), 4.0e-6,
                            4.0e-7, launcher);
  spatialDecompositionOuter(poly_ag, poly_ps, PrecisionModel::SINGLE, false, 4,
                            NonbondedTheme::ELECTROSTATIC, PrecisionModel::DOUBLE, 0,
                            5, NonbondedTheme::ELECTROSTATIC, tsm.getTestingStatus(), 1.0e-5,
                            1.0e-6, launcher);

  // Try a fixed-precision representation of the coordinates
  spatialDecompositionOuter(poly_ag, poly_ps, PrecisionModel::SINGLE, true, 4,
                            NonbondedTheme::ELECTROSTATIC, PrecisionModel::SINGLE, 0,
                            5, NonbondedTheme::ELECTROSTATIC, tsm.getTestingStatus(), 7.0e-6,
                            4.0e-7, launcher);

  // Try mapping dispersion sources
  spatialDecompositionOuter(poly_ag, poly_ps, PrecisionModel::DOUBLE, false, 4,
                            NonbondedTheme::VAN_DER_WAALS, PrecisionModel::DOUBLE, 0,
                            5, NonbondedTheme::VAN_DER_WAALS, tsm.getTestingStatus(), 1.0e-3,
                            1.0e-6, launcher);
  spatialDecompositionOuter(poly_ag, poly_ps, PrecisionModel::SINGLE, false, 4,
                            NonbondedTheme::VAN_DER_WAALS, PrecisionModel::SINGLE, 32,
                            5, NonbondedTheme::VAN_DER_WAALS, tsm.getTestingStatus(), 1.0e-3,
                            2.4e-6, launcher);

  // Check some input traps
  CellGrid<double, llint, double, double4> cg_test(poly_ps, poly_ag, 4.5, 0.25, 4,
                                                   NonbondedTheme::ELECTROSTATIC);
  CHECK_THROWS_SOFT(PMIGrid pm_bad(cg_test, NonbondedTheme::ELECTROSTATIC, 5,
                                   PrecisionModel::DOUBLE, FFTMode::OUT_OF_PLACE, 16), "A set of "
                    "particle-mesh interaction grids was created with an inaccurate "
                    "fixed-precision representation.", tsm.getTestingStatus());
  CHECK_THROWS_SOFT(PMIGrid pm_bad(cg_test, NonbondedTheme::ELECTROSTATIC, 5,
                                   PrecisionModel::SINGLE, FFTMode::OUT_OF_PLACE, 58), "A set "
                    "of " + getEnumerationName(NonbondedTheme::ELECTROSTATIC) + " particle-mesh "
                    "interaction grids was created with a risky fixed-precision representation.",
                    tsm.getTestingStatus());
  CHECK_THROWS_SOFT(PMIGrid pm_bad(cg_test, NonbondedTheme::VAN_DER_WAALS, 5,
                                   PrecisionModel::SINGLE, FFTMode::IN_PLACE, 54), "A set of " +
                    getEnumerationName(NonbondedTheme::VAN_DER_WAALS) + " particle-mesh "
                    "interaction grids was created with a risky fixed-precision representation.",
                    tsm.getTestingStatus());
  timer.assignTime(qspr_walltime);

  // Test particle-particle interactions using a simplified neighbor list
  section(3);
  const double ew_coeff = ewaldCoefficient(default_pme_cutoff, default_dsum_tol);
  std::vector<double2> nrg_results;
  const int nsys = tsm.getSystemCount();
  std::vector<std::vector<double>> qqx_ref_frc(nsys), qqy_ref_frc(nsys), qqz_ref_frc(nsys);
  std::vector<std::vector<double>> qqx_twr_frc(nsys), qqy_twr_frc(nsys), qqz_twr_frc(nsys);
  std::vector<std::vector<double>> ljx_ref_frc(nsys), ljy_ref_frc(nsys), ljz_ref_frc(nsys);
  std::vector<std::vector<double>> ljx_twr_frc(nsys), ljy_twr_frc(nsys), ljz_twr_frc(nsys);
  std::vector<std::vector<double>> nbx_ref_frc(nsys), nby_ref_frc(nsys), nbz_ref_frc(nsys);
  std::vector<std::vector<double>> nbx_twr_frc(nsys), nby_twr_frc(nsys), nbz_twr_frc(nsys);
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    bool skip = false;
    switch (tsm.getTopologyReference(i).getUnitCellType()) {
    case UnitCellType::NONE:
      skip = true;
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      break;
    }
    if (skip) {
      continue;
    }
    PhaseSpace ps = tsm.exportPhaseSpace(i);
    PhaseSpace ps_chk = ps;
    PhaseSpaceWriter psw = ps.data();
    const AtomGraph ag = tsm.getTopologyReference(i);
    const LocalExclusionMask lema_qq(ag, NonbondedTheme::ELECTROSTATIC);
    const LocalExclusionMask lema_lj(ag, NonbondedTheme::VAN_DER_WAALS);
    const LocalExclusionMask lema_all(ag, NonbondedTheme::ALL);
    const double2 pme_qqe  = evaluateParticleParticleEnergy(&ps, ag, lema_qq,
                                                            PrecisionModel::DOUBLE,
                                                            default_pme_cutoff, default_pme_cutoff,
                                                            ew_coeff, ew_coeff,
                                                            VdwSumMethod::CUTOFF,
                                                            EvaluateForce::YES,
                                                            NonbondedTheme::ELECTROSTATIC);
    checkParticlePairInteractions(&ps_chk, ps, ag, default_pme_cutoff, default_pme_cutoff,
                                  ew_coeff, ew_coeff, VdwSumMethod::CUTOFF, lema_qq, 0, pme_qqe,
                                  NonbondedTheme::ELECTROSTATIC, tsm.getTestingStatus());
    qqx_ref_frc[i].resize(psw.natom);
    qqy_ref_frc[i].resize(psw.natom);
    qqz_ref_frc[i].resize(psw.natom);
    for (int j = 0; j < psw.natom; j++) {
      qqx_ref_frc[i][j] = psw.xfrc[j];
      qqy_ref_frc[i][j] = psw.yfrc[j];
      qqz_ref_frc[i][j] = psw.zfrc[j];
    }
    ps.initializeForces();
    const double2 pme_lje  = evaluateParticleParticleEnergy(&ps, ag, lema_lj,
                                                            PrecisionModel::DOUBLE,
                                                            default_pme_cutoff, default_pme_cutoff,
                                                            ew_coeff, ew_coeff,
                                                            VdwSumMethod::CUTOFF,
                                                            EvaluateForce::YES,
                                                            NonbondedTheme::VAN_DER_WAALS);
    checkParticlePairInteractions(&ps_chk, ps, ag, default_pme_cutoff, default_pme_cutoff,
                                  ew_coeff, ew_coeff, VdwSumMethod::CUTOFF, lema_lj, 0, pme_lje,
                                  NonbondedTheme::VAN_DER_WAALS, tsm.getTestingStatus());
    ljx_ref_frc[i].resize(psw.natom);
    ljy_ref_frc[i].resize(psw.natom);
    ljz_ref_frc[i].resize(psw.natom);
    for (int j = 0; j < psw.natom; j++) {
      ljx_ref_frc[i][j] = psw.xfrc[j];
      ljy_ref_frc[i][j] = psw.yfrc[j];
      ljz_ref_frc[i][j] = psw.zfrc[j];
    }
    ps.initializeForces();
    const double2 pme_alle = evaluateParticleParticleEnergy(&ps, ag, lema_all,
                                                            PrecisionModel::DOUBLE,
                                                            default_pme_cutoff, default_pme_cutoff,
                                                            ew_coeff, ew_coeff,
                                                            VdwSumMethod::CUTOFF,
                                                            EvaluateForce::YES,
                                                            NonbondedTheme::ALL);
    checkParticlePairInteractions(&ps_chk, ps, ag, default_pme_cutoff, default_pme_cutoff,
                                  ew_coeff, ew_coeff, VdwSumMethod::CUTOFF, lema_all, 0, pme_alle,
                                  NonbondedTheme::ALL, tsm.getTestingStatus());
    nbx_ref_frc[i].resize(psw.natom);
    nby_ref_frc[i].resize(psw.natom);
    nbz_ref_frc[i].resize(psw.natom);
    for (int j = 0; j < psw.natom; j++) {
      nbx_ref_frc[i][j] = psw.xfrc[j];
      nby_ref_frc[i][j] = psw.yfrc[j];
      nbz_ref_frc[i][j] = psw.zfrc[j];
    }
    nrg_results.push_back(pme_alle);
  }
  
  // Run calculations based on the synthesis of all systems, then compare the results to those
  // obtained with the basic neighbor list method.
  const std::vector<NonbondedTheme> nbkinds = { NonbondedTheme::ELECTROSTATIC,
                                                NonbondedTheme::VAN_DER_WAALS,
                                                NonbondedTheme::ALL };
  std::vector<double> tower_plate_qq, tower_plate_lj, tower_plate_allqq, tower_plate_alllj;
  std::vector<std::vector<double>> tower_plate_qq_xfrc, tower_plate_qq_yfrc, tower_plate_qq_zfrc;
  std::vector<std::vector<double>> tower_plate_lj_xfrc, tower_plate_lj_yfrc, tower_plate_lj_zfrc;
  std::vector<std::vector<double>> tower_plate_xfrc, tower_plate_yfrc, tower_plate_zfrc;
  tower_plate_qq_xfrc.resize(nsys);
  tower_plate_qq_yfrc.resize(nsys);
  tower_plate_qq_zfrc.resize(nsys);
  tower_plate_lj_xfrc.resize(nsys);
  tower_plate_lj_yfrc.resize(nsys);
  tower_plate_lj_zfrc.resize(nsys);
  tower_plate_xfrc.resize(nsys);
  tower_plate_yfrc.resize(nsys);
  tower_plate_zfrc.resize(nsys);
  PsSynthesisWriter poly_psw = poly_ps.data();
  for (size_t i = 0; i < nbkinds.size(); i++) {
    CellGrid<double, llint, double, double4> cg(poly_ps, poly_ag, 0.5 * default_pme_cutoff, 0.1,
                                                4, nbkinds[i]);
    checkCellGridContents(cg, tsm.getTestingStatus());
    const LocalExclusionMask poly_lema(poly_ag, nbkinds[i]);
    ScoreCard sc(poly_ps.getSystemCount());
    evaluateParticleParticleEnergy<double, llint,
                                   double, double4>(&cg, &sc, poly_lema, default_pme_cutoff,
                                                    default_pme_cutoff, ew_coeff, ew_coeff,
                                                    VdwSumMethod::CUTOFF, EvaluateForce::YES,
                                                    nbkinds[i]);
    poly_ps.initializeForces();
    cg.contributeForces();
    const std::vector<double> qq_nrg = sc.reportInstantaneousStates(StateVariable::ELECTROSTATIC);
    const std::vector<double> lj_nrg = sc.reportInstantaneousStates(StateVariable::VDW);
    for (int j = 0; j < poly_ps.getSystemCount(); j++) {

      // Log the energies and forces
      const int k_llim = poly_psw.atom_starts[j];
      const int k_hlim = k_llim + poly_psw.atom_counts[j];
      switch (nbkinds[i]) {
      case NonbondedTheme::ELECTROSTATIC:
        tower_plate_qq.push_back(qq_nrg[j]);
        tower_plate_qq_xfrc[j].resize(poly_psw.atom_counts[j]);
        tower_plate_qq_yfrc[j].resize(poly_psw.atom_counts[j]);
        tower_plate_qq_zfrc[j].resize(poly_psw.atom_counts[j]);
        for (int k = k_llim; k < k_hlim; k++) {
          tower_plate_qq_xfrc[j][k - k_llim] = hostInt95ToDouble(poly_psw.xfrc[k],
                                                                 poly_psw.xfrc_ovrf[k]) *
                                               poly_psw.inv_frc_scale;
          tower_plate_qq_yfrc[j][k - k_llim] = hostInt95ToDouble(poly_psw.yfrc[k],
                                                                 poly_psw.yfrc_ovrf[k]) *
                                               poly_psw.inv_frc_scale;
          tower_plate_qq_zfrc[j][k - k_llim] = hostInt95ToDouble(poly_psw.zfrc[k],
                                                                 poly_psw.zfrc_ovrf[k]) *
                                               poly_psw.inv_frc_scale;
        }
        break;
      case NonbondedTheme::VAN_DER_WAALS:
        tower_plate_lj.push_back(lj_nrg[j]);
        tower_plate_lj_xfrc[j].resize(poly_psw.atom_counts[j]);
        tower_plate_lj_yfrc[j].resize(poly_psw.atom_counts[j]);
        tower_plate_lj_zfrc[j].resize(poly_psw.atom_counts[j]);
        for (int k = k_llim; k < k_hlim; k++) {
          tower_plate_lj_xfrc[j][k - k_llim] = hostInt95ToDouble(poly_psw.xfrc[k],
                                                                 poly_psw.xfrc_ovrf[k]) *
                                               poly_psw.inv_frc_scale;
          tower_plate_lj_yfrc[j][k - k_llim] = hostInt95ToDouble(poly_psw.yfrc[k],
                                                                 poly_psw.yfrc_ovrf[k]) *
                                               poly_psw.inv_frc_scale;
          tower_plate_lj_zfrc[j][k - k_llim] = hostInt95ToDouble(poly_psw.zfrc[k],
                                                                 poly_psw.zfrc_ovrf[k]) *
                                               poly_psw.inv_frc_scale;
        }
        break;
      case NonbondedTheme::ALL:
        tower_plate_allqq.push_back(qq_nrg[j]);
        tower_plate_alllj.push_back(lj_nrg[j]);
        tower_plate_xfrc[j].resize(poly_psw.atom_counts[j]);
        tower_plate_yfrc[j].resize(poly_psw.atom_counts[j]);
        tower_plate_zfrc[j].resize(poly_psw.atom_counts[j]);
        for (int k = k_llim; k < k_hlim; k++) {
          tower_plate_xfrc[j][k - k_llim] = hostInt95ToDouble(poly_psw.xfrc[k],
                                                              poly_psw.xfrc_ovrf[k]) *
                                            poly_psw.inv_frc_scale;
          tower_plate_yfrc[j][k - k_llim] = hostInt95ToDouble(poly_psw.yfrc[k],
                                                              poly_psw.yfrc_ovrf[k]) *
                                            poly_psw.inv_frc_scale;
          tower_plate_zfrc[j][k - k_llim] = hostInt95ToDouble(poly_psw.zfrc[k],
                                                              poly_psw.zfrc_ovrf[k]) *
                                            poly_psw.inv_frc_scale;
        }
        break;
      }
    }
  }

  // Check that the synthesis forces agree with the reference calculations
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    compareForceArrays(tower_plate_qq_xfrc[i], qqx_ref_frc[i], CartesianDimension::X,
                       NonbondedTheme::ELECTROSTATIC, tsm.getTestingStatus());
    compareForceArrays(tower_plate_qq_yfrc[i], qqy_ref_frc[i], CartesianDimension::Y,
                       NonbondedTheme::ELECTROSTATIC, tsm.getTestingStatus());
    compareForceArrays(tower_plate_qq_zfrc[i], qqz_ref_frc[i], CartesianDimension::Z,
                       NonbondedTheme::ELECTROSTATIC, tsm.getTestingStatus());
    compareForceArrays(tower_plate_lj_xfrc[i], ljx_ref_frc[i], CartesianDimension::X,
                       NonbondedTheme::VAN_DER_WAALS, tsm.getTestingStatus());
    compareForceArrays(tower_plate_lj_yfrc[i], ljy_ref_frc[i], CartesianDimension::Y,
                       NonbondedTheme::VAN_DER_WAALS, tsm.getTestingStatus());
    compareForceArrays(tower_plate_lj_zfrc[i], ljz_ref_frc[i], CartesianDimension::Z,
                       NonbondedTheme::VAN_DER_WAALS, tsm.getTestingStatus());
    compareForceArrays(tower_plate_xfrc[i], nbx_ref_frc[i], CartesianDimension::X,
                       NonbondedTheme::ALL, tsm.getTestingStatus());
    compareForceArrays(tower_plate_yfrc[i], nby_ref_frc[i], CartesianDimension::Y,
                       NonbondedTheme::ALL, tsm.getTestingStatus());
    compareForceArrays(tower_plate_zfrc[i], nbz_ref_frc[i], CartesianDimension::Z,
                       NonbondedTheme::ALL, tsm.getTestingStatus());
  }

  // Perform finite difference tests
  finiteDifferenceNeighborListTest<double, llint, double, double4>(tsm, nbkinds, 1024, 9);
  timer.assignTime(ppdir_walltime);
  
  // Check the HPC kernels
#ifdef STORMM_USE_HPC
  runGpuTests(tsm, testdir, &xrs, gpu);
  timer.assignTime(gpu_walltime);
#endif

  // Test some additional PME-related functions
  section(4);
  runMiscellaneousTests();
  
  // Summary evaluation
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
