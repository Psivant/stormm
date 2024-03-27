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
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/math_enumerators.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/MolecularMechanics/mm_controls.h"
#include "../../src/Potential/cellgrid.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/map_density.h"
#include "../../src/Potential/pmigrid.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::card;
using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::energy;
using namespace stormm::mm;
using namespace stormm::review;
using namespace stormm::stmath;
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
        if (fabs(synbk.ljb_coeff[tsys_offset + ((nlj_typ + 1) * tlj_idx)]) < 1.0e-6) {
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
                                cg_tmat, pm_order, *poly_ag);
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
      CellGrid<llint, llint, double, llint4> cg(poly_ps, poly_ag, 4.5, 0.25, mesh_ticks, cg_theme);
      inspectCellGrids(cg, poly_ag, do_tests);
      spatialDecompositionInner<llint, llint, double, llint4>(cg, pmi_grid_prec, pm_acc_bits,
                                                              pm_order, pm_theme, do_tests,
                                                              chrg_conserv_tol, chrg_value_tol,
                                                              launcher);
    }
    else {
      CellGrid<double, llint, double, double4> cg(poly_ps, poly_ag, 4.5, 0.25, mesh_ticks,
                                                  cg_theme);
      inspectCellGrids(cg, poly_ag, do_tests);
      spatialDecompositionInner<double, llint, double, double4>(cg, pmi_grid_prec, pm_acc_bits,
                                                                pm_order, pm_theme, do_tests,
                                                                chrg_conserv_tol, chrg_value_tol,
                                                                launcher);
    }
    break;
  case PrecisionModel::SINGLE:
    if (cast_cg_in_fp) {
      CellGrid<int, int, float, int4> cg(poly_ps, poly_ag, 4.5, 0.25, mesh_ticks, cg_theme);
      inspectCellGrids(cg, poly_ag, do_tests);
      spatialDecompositionInner<int, int, float, int4>(cg, pmi_grid_prec, pm_acc_bits, pm_order,
                                                       pm_theme, do_tests, chrg_conserv_tol,
                                                       chrg_value_tol, launcher);
    }
    else {
      CellGrid<float, int, float, float4> cg(poly_ps, poly_ag, 4.5, 0.25, mesh_ticks, cg_theme);
      inspectCellGrids(cg, poly_ag, do_tests);
      spatialDecompositionInner<float, int, float, float4>(cg, pmi_grid_prec, pm_acc_bits,
                                                           pm_order, pm_theme, do_tests,
                                                           chrg_conserv_tol, chrg_value_tol,
                                                           launcher);
    }
    break;
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
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> create_this_to_engage_gpu(1);
#endif
  
  // Section 1                                                    
  section("Test mechanics in the PMIGrid object");

  // Section 2
  section("Test charge spreading in various modes");

  // Create a synthesis of systems and the associated particle-mesh interaction grid
  section(1);
  const std::vector<std::string> pbc_systems = { "bromobenzene", "bromobenzene_vs", "tip3p",
                                                 "tip4p", "ubiquitin", "tamavidin", "drug_example",
                                                 "trpcage_in_water" };
  const char osc = osSeparator();
  const std::string testdir = oe.getStormmSourcePath() + osc + "test" + osc;
  TestSystemManager tsm(testdir + "Topology", "top", pbc_systems, testdir + "Trajectory", "inpcrd",
                        pbc_systems);

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
  // Test the baseline, double / double case.  No fixed-precision accumulation.
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
