// -*-c++-*-
#include "../../src/Accelerator/core_kernel_manager.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Constants/scaling.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/MolecularMechanics/dynamics.h"
#include "../../src/MolecularMechanics/kinetic.h"
#include "../../src/MolecularMechanics/minimization.h"
#include "../../src/MolecularMechanics/mm_controls.h"
#include "../../src/Potential/cacheresource.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/hpc_valence_potential.h"
#include "../../src/Potential/hpc_nonbonded_potential.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Random/random.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/implicit_solvent_workspace.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/nonbonded_workunit.h"
#include "../../src/Synthesis/static_mask_synthesis.h"
#include "../../src/Synthesis/valence_workunit.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_intake.h"
#include "../../src/Trajectory/coordinate_intake.h"
#include "../../src/Trajectory/coordinate_series.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/thermostat.h"
#include "../../src/Trajectory/trajectory_enumerators.h"
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
using namespace stormm::random;
using namespace stormm::stmath;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Carry out one step for dynamics in a particular system usign the GPU, then replicate the process
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
//   rng_seed:    Random number seed for the Anderser velocity kick-start
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
              std::to_string(step_idx + 1) + " steps.", tsm.getTestingStatus());
        const std::vector<double> cpu_vel_curr = ps_cpu[i].getInterlacedCoordinates(tvel);
        const std::vector<double> gpu_vel_curr = ps_gpu[i].getInterlacedCoordinates(tvel);
        check(gpu_vel_curr, RelationalOperator::EQUAL, cpu_vel_curr, "The GPU and CPU do not "
              "indicate the same velocities in replica " + std::to_string(i) + " after " +
              std::to_string(step_idx + 1) + " steps.", tsm.getTestingStatus());

        // Check the previous coordinates
        ps_cpu[i].updateCyclePosition();
        ps_gpu[i].updateCyclePosition();
        const std::vector<double> cpu_pos_orig = ps_cpu[i].getInterlacedCoordinates(tpos);
        const std::vector<double> gpu_pos_orig = ps_gpu[i].getInterlacedCoordinates(tpos);
        check(gpu_pos_orig, RelationalOperator::EQUAL, cpu_pos_orig, "The GPU and CPU do not "
              "indicate the same reference particle positions in replica " + std::to_string(i) +
              " after " + std::to_string(step_idx + 1) + " steps.", tsm.getTestingStatus());
        const std::vector<double> cpu_vel_orig = ps_cpu[i].getInterlacedCoordinates(tvel);
        const std::vector<double> gpu_vel_orig = ps_gpu[i].getInterlacedCoordinates(tvel);
        check(gpu_vel_orig, RelationalOperator::EQUAL, cpu_vel_orig, "The GPU and CPU do not "
              "indicate the same reference velocities in replica " + std::to_string(i) +
              " after " + std::to_string(step_idx + 1) + " steps.", tsm.getTestingStatus());
        ps_cpu[i].updateCyclePosition();
        ps_gpu[i].updateCyclePosition();
      }
    }
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
#if 0
  check(gpu_result, RelationalOperator::EQUAL, Approx(cpu_result).margin(1.0e-6),
        "Random numbers generated by the CPU and GPU do not agree for a Langevin thermostat.");
#endif
  // Display timings and test results
  if (oe.getDisplayTimingsOrder()) {
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
  return countGlobalTestFailures();
}
