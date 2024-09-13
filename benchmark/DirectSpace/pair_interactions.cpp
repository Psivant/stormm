#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Accelerator/core_kernel_manager.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/Math/math_enumerators.h"
#include "../../src/MolecularMechanics/mm_controls.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/cellgrid.h"
#include "../../src/Potential/energy_enumerators.h"
#ifdef STORMM_USE_HPC
#  include "../../src/Potential/hpc_pme_potential.h"
#endif
#include "../../src/Potential/local_exclusionmask.h"
#include "../../src/Potential/pme_potential.h"
#include "../../src/Potential/pme_util.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/hpc_phasespace_synthesis.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinate_copy.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/unit_test_enumerators.h"

using namespace stormm::card;
using namespace stormm::data_types;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::parse;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;

//-------------------------------------------------------------------------------------------------
// Check the forces computed by the GPU kernel against CPU results.
//
// Arguments:
//   poly_ps:    The synthesis of coordinates, including forces on all particles
//   eval_frc:   Indicate whether forces were evaluated as part of the run
//   timer:      Time tracking object, to absorb the time running the prior CPU calculation as well
//               as the test to compare the results
//   time_test:  The integer code for the test time tracking
//-------------------------------------------------------------------------------------------------
void checkPairBenchmarkForces(const PhaseSpaceSynthesis &poly_ps, const EvaluateForce eval_frc,
                              StopWatch *timer, const int time_test) {
#ifdef STORMM_USE_HPC
  const UnitCellType system_uc = poly_ps.getUnitCellType();
  CoordinateFrame cpu_frc(poly_ps.getAtomCount(0), system_uc, HybridFormat::HOST_MOUNTED);
  CoordinateFrame gpu_frc(poly_ps.getAtomCount(0), system_uc, HybridFormat::HOST_MOUNTED);
  const CoordinateFrameReader cpu_cfr = cpu_frc.data();
  const CoordinateFrameReader gpu_cfr = gpu_frc.data();
  std::vector<double> cpu_chk(cpu_cfr.natom, 0.0);
  std::vector<double> gpu_chk(cpu_cfr.natom, 0.0);
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      coordCopy(&cpu_frc, poly_ps, i, TrajectoryKind::FORCES, HybridTargetLevel::HOST,
                HybridTargetLevel::HOST);
      coordCopy(&gpu_frc, poly_ps, i, TrajectoryKind::FORCES, HybridTargetLevel::HOST,
                HybridTargetLevel::DEVICE);
      for (int j = 0; j < cpu_cfr.natom; j++) {
        cpu_chk[j] = cpu_cfr.xcrd[j];
        gpu_chk[j] = gpu_cfr.xcrd[j];
      }
      check(gpu_chk, RelationalOperator::EQUAL, Approx(cpu_chk).margin(1.0e-2), "Forces on "
            "particles along the Cartesian X axis calculated using the GPU kernel running in " +
            getEnumerationName(PrecisionModel::SINGLE) + "-precision do not agree with CPU "
            "forces calculated in " + getEnumerationName(PrecisionModel::DOUBLE) +
            "-precision.");
      for (int j = 0; j < cpu_cfr.natom; j++) {
        cpu_chk[j] = cpu_cfr.ycrd[j];
        gpu_chk[j] = gpu_cfr.ycrd[j];
      }
      check(gpu_chk, RelationalOperator::EQUAL, Approx(cpu_chk).margin(1.0e-2), "Forces on "
            "particles along the Cartesian Y axis calculated using the GPU kernel running in " +
            getEnumerationName(PrecisionModel::SINGLE) + "-precision do not agree with CPU "
            "forces calculated in " + getEnumerationName(PrecisionModel::DOUBLE) +
            "-precision.");
      for (int j = 0; j < cpu_cfr.natom; j++) {
        cpu_chk[j] = cpu_cfr.zcrd[j];
        gpu_chk[j] = gpu_cfr.zcrd[j];
      }
      check(gpu_chk, RelationalOperator::EQUAL, Approx(cpu_chk).margin(1.0e-2), "Forces on "
            "particles along the Cartesian Z axis calculated using the GPU kernel running in " +
            getEnumerationName(PrecisionModel::SINGLE) + "-precision do not agree with CPU "
            "forces calculated in " + getEnumerationName(PrecisionModel::DOUBLE) +
            "-precision.");
      break;
    case EvaluateForce::NO:
      break;
    }
  }
#endif
  timer->assignTime(time_test);
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main (const int argc, const char* argv[]) {

  // Baseline variables
  TestEnvironment oe(argc, argv, ExceptionResponse::SILENT);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  StopWatch timer;
#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> force_gpu_to_engage(1);

  // Lay out time categories for program profiling
  const int time_load  = timer.addCategory("Load input files");
  const int time_build = timer.addCategory("Class object setup");
  const int time_test  = timer.addCategory("CPU testing");
  const int time_init  = timer.addCategory("Initialization");
  const int time_pairs = timer.addCategory("Pairs kernel with initialization");
  
  // Take in additional command line inputs
  int n_trials = 4;
  int n_repeats = 100;
  int n_replica = 1;
  int fp_bits = 24;
  int warp_mult = 1;
  int ig_seed = 616804732;
  double elec_cutoff = default_pme_cutoff;
  double vdw_cutoff = default_pme_cutoff;
  double cutoff_pad = 0.15;
  EvaluateForce eval_frc = EvaluateForce::YES;
  EvaluateEnergy eval_nrg = EvaluateEnergy::NO;
  bool use_dual_cg = false;
  std::string inpcrd_file, topology_file;
  for (int i = 0; i < argc; i++) {

    // Keywords that alter behavior
    if (strcmpCased(argv[i], "-dual_cg", CaseSensitivity::NO)) {
      use_dual_cg = true;
    }
    
    // Keywords that are prompts for a value input
    if (i < argc - 1) {
      if (strcmpCased(argv[i], "-trials", CaseSensitivity::NO) &&
          verifyNumberFormat(argv[i + 1], NumberFormat::INTEGER)) {
        n_trials = atoi(argv[i + 1]);
      }
      else if (strcmpCased(argv[i], "-iter", CaseSensitivity::NO) &&
               verifyNumberFormat(argv[i + 1], NumberFormat::INTEGER)) {
        n_repeats = atoi(argv[i + 1]);
      }
      else if (strcmpCased(argv[i], "-elec_cutoff", CaseSensitivity::NO) &&
               (verifyNumberFormat(argv[i + 1], NumberFormat::STANDARD_REAL) ||
                verifyNumberFormat(argv[i + 1], NumberFormat::SCIENTIFIC))) {
        elec_cutoff = atof(argv[i + 1]);
      }
      else if (strcmpCased(argv[i], "-vdw_cutoff", CaseSensitivity::NO) &&
               (verifyNumberFormat(argv[i + 1], NumberFormat::STANDARD_REAL) ||
                verifyNumberFormat(argv[i + 1], NumberFormat::SCIENTIFIC))) {
        vdw_cutoff = atof(argv[i + 1]);
      }
      else if (strcmpCased(argv[i], "-pad", CaseSensitivity::NO) &&
               (verifyNumberFormat(argv[i + 1], NumberFormat::STANDARD_REAL) ||
                verifyNumberFormat(argv[i + 1], NumberFormat::SCIENTIFIC))) {
        cutoff_pad = atof(argv[i + 1]);
      }
      else if (strcmpCased(argv[i], "-replicas", CaseSensitivity::NO)) {
        n_replica = atof(argv[i + 1]);
      }
      else if (strcmpCased(argv[i], "-fp_bits", CaseSensitivity::NO)) {
        fp_bits = atof(argv[i + 1]);
      }
      else if (strcmpCased(argv[i], "-nt_warps", CaseSensitivity::NO)) {
        warp_mult = atof(argv[i + 1]);
      }
      else if (strcmpCased(argv[i], "-ig_seed", CaseSensitivity::NO)) {
        ig_seed = atoi(argv[i + 1]);
      }
      else if (strcmpCased(argv[i], "-c", CaseSensitivity::NO)) {
        inpcrd_file = argv[i + 1];
      }
      else if (strcmpCased(argv[i], "-p", CaseSensitivity::NO)) {
        topology_file = argv[i + 1];
      }
      else if (strcmpCased(argv[i], "-frc", CaseSensitivity::NO)) {
        if (strcmpCased(argv[i + 1], "no", CaseSensitivity::NO)) {
          eval_frc = EvaluateForce::NO;
        }
        else if (strcmpCased(argv[i + 1], "yes", CaseSensitivity::NO)) {
          eval_frc = EvaluateForce::YES;
        }
        else {
          rtErr("Bad argument to -frc '" + std::string(argv[i + 1]) + "'.", "pair_interactions");
        }
      }
      else if (strcmpCased(argv[i], "-nrg", CaseSensitivity::NO)) {
        if (strcmpCased(argv[i + 1], "no", CaseSensitivity::NO)) {
          eval_nrg = EvaluateEnergy::NO;
        }
        else if (strcmpCased(argv[i + 1], "yes", CaseSensitivity::NO)) {
          eval_nrg = EvaluateEnergy::YES;
        }
        else {
          rtErr("Bad argument to -nrg '" + std::string(argv[i + 1]) + "'.", "pair_interactions");
        }
      }
    }
  }

  // Input checks
  if (n_replica <= 0) {
    rtErr("A replica count of " + std::to_string(n_replica) + " is invalid.\n",
          "pair_interactions");
  }
  
  // A Hybrid object was created to engage the GPU.  Absorb any bootup time into "miscellaneous."
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
  }
  
  // Check that a system has been provided.  Load the basic topology and input coordinates.
  TestSystemManager tsm(std::vector<std::string>(1, topology_file),
                        std::vector<std::string>(1, inpcrd_file), ExceptionResponse::DIE);
  timer.assignTime(time_load);
  
  // Stage the input parameters: cutoffs for each interaction.
  MolecularMechanicsControls mmctrl(0.01, 1, 1, warp_mult, elec_cutoff, vdw_cutoff);

  // Create the basic objects
  PhaseSpaceSynthesis poly_ps = tsm.exportPhaseSpaceSynthesis(std::vector<int>(n_replica, 0));
  PsSynthesisWriter host_psw = poly_ps.data();
  Xoshiro256ppGenerator xrs(ig_seed);
  for (int i = 1; i < host_psw.system_count; i++) {
    const size_t aoff = host_psw.atom_starts[i];
    addRandomNoise(&xrs, &host_psw.xcrd[aoff], &host_psw.xcrd_ovrf[aoff], &host_psw.ycrd[aoff],
                   &host_psw.ycrd_ovrf[aoff], &host_psw.zcrd[aoff], &host_psw.zcrd_ovrf[aoff],
                   host_psw.atom_counts[i], 0.01, host_psw.gpos_scale_f);
  }
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(std::vector<int>(n_replica, 0));
  ScoreCard sc(poly_ps.getSystemCount(), 1, 32);
  use_dual_cg = (use_dual_cg || (fabs(elec_cutoff - vdw_cutoff) > 1.0e-6));
  if (use_dual_cg == false) {
    elec_cutoff = vdw_cutoff;
  }
  const NeighborListKind layout = (use_dual_cg) ? NeighborListKind::DUAL : NeighborListKind::MONO;
  CoreKlManager launcher(gpu, poly_ag);
  PPITable direct_space_table(NonbondedTheme::ELECTROSTATIC, BasisFunctions::MIXED_FRACTIONS,
                              TableIndexing::SQUARED_ARG, elec_cutoff);
  const double ew_coeff = direct_space_table.getEwaldCoefficient();
  timer.assignTime(time_build);
  
  // Upload the basic objects
  poly_ps.upload();
  poly_ag.upload();
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps.data(devc);
  SyNonbondedKit<float, float2> poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(devc);
  MMControlKit<float> ctrl = mmctrl.spData(devc);
  PPIKit<float, float4> nrg_tab = direct_space_table.spData(devc);
  ScoreCardWriter scw = sc.data(devc);
  timer.assignTime(time_build);

  // Create the neighbor list cell grids and run calculations
  if (use_dual_cg) {
    CellGrid<float, int, float, float4> cg_qq(&poly_ps, &poly_ag, 0.5 * elec_cutoff, cutoff_pad, 4,
                                              NonbondedTheme::ELECTROSTATIC, gpu);
    CellGrid<float, int, float, float4> cg_lj(&poly_ps, &poly_ag, 0.5 * vdw_cutoff, cutoff_pad, 4,
                                              NonbondedTheme::VAN_DER_WAALS, gpu);
    const TinyBoxPresence has_tiny_box = (cg_qq.getTinyBoxPresence() == TinyBoxPresence::YES ||
                                          cg_lj.getTinyBoxPresence() == TinyBoxPresence::YES) ?
                                         TinyBoxPresence::YES : TinyBoxPresence::NO;
    const int2 tp_bt = launcher.getPMEPairsKernelDims(PrecisionModel::SINGLE,
                                                      PrecisionModel::SINGLE,
                                                      NeighborListKind::DUAL, has_tiny_box,
                                                      eval_frc, eval_nrg, ClashResponse::NONE,
                                                      PairStance::TOWER_PLATE);
    const int2 tt_bt = launcher.getPMEPairsKernelDims(PrecisionModel::SINGLE,
                                                      PrecisionModel::SINGLE,
                                                      NeighborListKind::DUAL, has_tiny_box,
                                                      eval_frc, eval_nrg, ClashResponse::NONE,
                                                      PairStance::TOWER_TOWER);
    TileManager tlmn(tp_bt);
    TilePlan tlpn = tlmn.data(devc);
    mmctrl.primeWorkUnitCounters(launcher, eval_frc, eval_nrg, ClashResponse::NONE,
                                 VwuGoal::MOVE_PARTICLES, PrecisionModel::SINGLE,
                                 PrecisionModel::SINGLE, QMapMethod::ACC_SHARED,
                                 PrecisionModel::SINGLE, float_type_index, 5, layout, has_tiny_box,
                                 poly_ag);
    mmctrl.upload();
    cg_qq.upload();
    cg_lj.upload();

    // Make the appropriate exclusion masks
    LocalExclusionMask lem(poly_ag, NonbondedTheme::ALL);
    lem.upload();
    const LocalExclusionMaskReader lemr = lem.data(devc);
    
    // Take abstracts and run the kernel repeatedly.
    CellGridWriter<float, int, float, float4> cg_qqw = cg_qq.data(devc);
    CellGridWriter<float, int, float, float4> cg_ljw = cg_lj.data(devc);
    timer.assignTime(time_build);

    // Run the initialization of forces multiple times to gain a sense of the overhead paid in the
    // test.
    for (int i = 0; i < n_trials; i++) {
      for (int j = 0; j < n_repeats; j++) {
        cg_qq.initializeForces(devc, gpu);
        cg_lj.initializeForces(devc, gpu);
      }
      cudaDeviceSynchronize();
      timer.assignTime(time_init);
    }

    // Run the pair interactions evaluation and continue clearing the buffers
    for (int i = 0; i < n_trials; i++) {
      for (int j = 0; j < n_repeats; j++) {
        cg_qq.initializeForces(devc, gpu);
        cg_lj.initializeForces(devc, gpu);
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cg_qqw, &cg_ljw, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, has_tiny_box, tp_bt, tt_bt, 0.0, 0.0);
        ctrl.step +=1;
      }
      cudaDeviceSynchronize();
      timer.assignTime(time_pairs);
    }
    poly_ps.initializeForces(gpu, devc);
    cg_qq.contributeForces(devc, gpu);
    cg_lj.contributeForces(devc, gpu);
  }
  else {
    CellGrid<float, int, float, float4> cg(&poly_ps, &poly_ag, 0.5 * vdw_cutoff, cutoff_pad, 4,
                                           NonbondedTheme::ALL, gpu);
    const TinyBoxPresence has_tiny_box = cg.getTinyBoxPresence();
    const int2 tp_bt = launcher.getPMEPairsKernelDims(PrecisionModel::SINGLE,
                                                      PrecisionModel::SINGLE,
                                                      NeighborListKind::DUAL,
                                                      cg.getTinyBoxPresence(), eval_frc, eval_nrg,
                                                      ClashResponse::NONE,
                                                      PairStance::TOWER_PLATE);
    const int2 tt_bt = launcher.getPMEPairsKernelDims(PrecisionModel::SINGLE,
                                                      PrecisionModel::SINGLE,
                                                      NeighborListKind::DUAL,
                                                      cg.getTinyBoxPresence(), eval_frc, eval_nrg,
                                                      ClashResponse::NONE,
                                                      PairStance::TOWER_TOWER);
    TileManager tlmn(tp_bt);
    TilePlan tlpn = tlmn.data(devc);
    mmctrl.primeWorkUnitCounters(launcher, eval_frc, eval_nrg, ClashResponse::NONE,
                                 VwuGoal::MOVE_PARTICLES, PrecisionModel::SINGLE,
                                 PrecisionModel::SINGLE, QMapMethod::ACC_SHARED,
                                 PrecisionModel::SINGLE, float_type_index, 5, layout,
                                 cg.getTinyBoxPresence(), poly_ag);
    mmctrl.upload();
    cg.upload();

    // Make the appropriate exclusion masks
    LocalExclusionMask lem(poly_ag, NonbondedTheme::ALL);
    lem.upload();
    const LocalExclusionMaskReader lemr = lem.data(devc);
    
    // Take abstracts and run the kernel repeatedly.
    CellGridWriter<float, int, float, float4> cgw = cg.data(devc);
    timer.assignTime(time_build);

    // Test the basic force initialization
    for (int i = 0; i < n_trials; i++) {
      for (int j = 0; j < n_repeats; j++) {
        cg.initializeForces(devc, gpu);
      }
      cudaDeviceSynchronize();
      timer.assignTime(time_init);
    }

    // Run the pair interactions evaluation and continue clearing the buffers
    for (int i = 0; i < n_trials; i++) {
      for (int j = 0; j < n_repeats; j++) {
        cg.initializeForces(devc, gpu);
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw, &tlpn, &scw, &ctrl, eval_frc, eval_nrg,
                       cg.getTinyBoxPresence(), tp_bt, tt_bt, 0.0, 0.0);
        ctrl.step +=1;
      }
      cudaDeviceSynchronize();
      timer.assignTime(time_pairs);
    }
    poly_ps.initializeForces(gpu, devc);
    cg.contributeForces(devc, gpu);
  }

  // Run a separate check
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    PhaseSpace ps_i(poly_ps.getAtomCount(i), poly_ps.getUnitCellType(),
                    HybridFormat::HOST_ONLY);
    coordCopy(&ps_i, poly_ps, i);
    ps_i.initializeForces();
    const AtomGraph *ag_i = poly_ps.getSystemTopologyPointer(i);
    const LocalExclusionMask lem_i(ag_i);
    evaluateParticleParticleEnergy(&ps_i, ag_i, lem_i, PrecisionModel::SINGLE, elec_cutoff,
                                   vdw_cutoff, ew_coeff);
    coordCopy(&poly_ps, i, ps_i);
  }
  checkPairBenchmarkForces(poly_ps, eval_frc, &timer, time_test);
#else // STORMM_USE_HPC
  rtWarn("This benchmarking program requires GPU support to run.", "pair_interactions");
#endif // STORMM_USE_HPC

  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return 0;
}
