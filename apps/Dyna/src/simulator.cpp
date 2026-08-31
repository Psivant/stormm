#include <vector>
#include "copyright.h"
#include "../../../src/Accelerator/gpu_details.h"
#include "../../../src/Analysis/hydrogen_bond_analysis.h"
#ifdef STORMM_USE_HPC
#  include "../../../src/Accelerator/core_kernel_manager.h"
#  include "../../../src/Accelerator/hpc_config.h"
#  include "../../../src/MolecularMechanics/hpc_dynamics.h"
#  include "../../../src/MolecularMechanics/hpc_minimization.h"
#else
#  include "../../../src/MolecularMechanics/dynamics.h"
#  include "../../../src/MolecularMechanics/minimization.h"
#endif
#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/Debug/hpc_debug.h"
#include "../../../src/DataTypes/common_types.h"
#include "../../../src/DataTypes/stormm_vector_types.h"
#include "../../../src/Debug/pruning.h"
#include "../../../src/Math/vector_ops.h"
#include "../../../src/MolecularMechanics/dynamics_intervention.h"
#include "../../../src/MolecularMechanics/minimization.h"
#include "../../../src/Namelists/input_transcript.h"
#include "../../../src/Namelists/nml_dynamics.h"
#include "../../../src/Namelists/nml_minimize.h"
#include "../../../src/Namelists/nml_pppm.h"
#include "../../../src/Namelists/nml_precision.h"
#include "../../../src/Namelists/nml_random.h"
#include "../../../src/Namelists/nml_remd.h"
#include "../../../src/Namelists/user_settings.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Potential/energy_enumerators.h"
#include "../../../src/Potential/local_exclusionmask.h"
#include "../../../src/Potential/static_exclusionmask.h"
#include "../../../src/Reporting/error_format.h"
#include "../../../src/Reporting/help_messages.h"
#include "../../../src/Reporting/present_analysis.h"
#include "../../../src/Reporting/present_debug.h"
#include "../../../src/Reporting/present_energy.h"
#include "../../../src/Reporting/progress_bar.h"
#include "../../../src/Debug/watcher.h"
#include "../../../src/Sampling/exchange_nexus.h"
#include "../../../src/Synthesis/atomgraph_synthesis.h"
#include "../../../src/Synthesis/cull_synthesis.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/static_mask_synthesis.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/Topology/atomgraph_enumerators.h"
#include "../../../src/Trajectory/phasespace.h"
#include "../../../src/Trajectory/thermostat.h"
#include "../../../src/Trajectory/trajectory_enumerators.h"
#include "../../../src/UnitTesting/stopwatch.h"
#include "../../../src/UnitTesting/unit_test.h"
#include "setup.h"

using namespace stormm::analysis;
using namespace stormm::card;
using namespace stormm::chemistry;
using namespace stormm::data_types;
using namespace stormm::debug;
using namespace stormm::display;
using namespace stormm::energy;
using namespace stormm::mm;
using namespace stormm::namelist;
using namespace stormm::parse;
using namespace stormm::random;
using namespace stormm::reporting;
using namespace stormm::restraints;
using namespace stormm::review;
using namespace stormm::sampling;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;
using namespace dyna_app::setup;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Wall time tracking
  StopWatch timer("Timings for dynamics.stormm");
  const int file_parse_tm = timer.addCategory("File Parsing");
  const int gen_setup_tm  = timer.addCategory("Setup, General");
  const int min_setup_tm  = timer.addCategory("Setup, Minimization");
  const int min_run_tm    = timer.addCategory("Run, Minimization");
  const int dyn_setup_tm  = timer.addCategory("Setup, Dynamics");
  const int dyn_run_tm    = timer.addCategory("Run, Dynamics");
  const int download_tm   = timer.addCategory("GPU Data Download");
  const int output_tm     = timer.addCategory("Trajectory Output");
  
  // Engage the GPU
#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  const Hybrid<int> array_to_trigger_gpu_mapping(1);
  dyna_tk.setGpu(gpu);
  const HybridTargetLevel operating_tier = HybridTargetLevel::DEVICE;
#else
  const GpuDetails gpu = null_gpu;
  const HybridTargetLevel operating_tier = HybridTargetLevel::HOST;
#endif
  timer.assignTime(gen_setup_tm);
  
  // Parse the command line
  CommandLineParser clip("dynamics.stormm", "The principal molecular dynamics engine in STORMM.");
  clip.addStandardApplicationInputs();
  const std::vector<std::string> my_namelists = { "&files", "&minimize", "&dynamics", "&remd",
                                                  "&restraint", "&solvent", "&random", "&report",
                                                  "&precision", "&debug", "&analysis" };
  clip.addControlBlocks(my_namelists);
  if (displayNamelistHelp(argc, argv, my_namelists) && clip.doesProgramExitOnHelp()) {
    return 0;
  }
  clip.parseUserInput(argc, argv);
  
  // Read information from the command line and initialize the UserSettings object
  UserSettings ui(clip, { "-pe", "-ce" });
  
  // Read topologies and coordinate files.  Assemble critical details about each system.
  SystemCache sc(ui.getFilesNamelistInfo(), ui.getRestraintNamelistInfo(),
                 ui.getDynamicsNamelistInfo(), ui.getExceptionBehavior(), MapRotatableGroups::NO,
                 ui.getPrintingPolicy());
  timer.assignTime(file_parse_tm);
  
  // Prepare a synthesis of systems from the user input.
  const std::vector<AtomGraph*> agv = sc.getTopologyPointer();

  // Preview the implicit solvent model.
  NeckGeneralizedBornTable ngb_tab;
  if (ui.getSolventPresence()) {
    const SolventControls& isvcon = ui.getSolventNamelistInfo();
    for (int i = 0; i < sc.getTopologyCount(); i++) {
      AtomGraph *ag = sc.getTopologyPointer(i);
      ag->setImplicitSolventModel(isvcon.getImplicitSolventModel(), isvcon.getExternalDielectric(),
                                  isvcon.getSaltConcentration(), isvcon.getPBRadiiSet(),
                                  ui.getExceptionBehavior());
    }
  }

  // Create the synthesis of systems, including exclusion masks and non-bonded work units as
  // necessary.
  int system_count = sc.getSystemCount();
  const PrecisionControls& preccon = ui.getPrecisionNamelistInfo();
  const DynamicsControls& dyncon = ui.getDynamicsNamelistInfo();
  const PPPMControls& pmecon = ui.getPPPMNamelistInfo();
  const ThermostatKind tstat_choice = dyncon.getThermostatKind();
  PhaseSpaceSynthesis poly_ps = sc.exportCoordinateSynthesis(preccon.getGlobalPosScalingBits(),
                                                             preccon.getVelocityScalingBits(),
                                                             preccon.getForceScalingBits());
  AtomGraphSynthesis poly_ag = sc.exportTopologySynthesis(gpu, ui.getExceptionBehavior());
  if (ui.getSolventPresence()) {
    const SolventControls& isvcon = ui.getSolventNamelistInfo();
    poly_ag.setImplicitSolventModel(isvcon.getImplicitSolventModel(), ngb_tab,
                                    isvcon.getPBRadiiSet(), isvcon.getExternalDielectric(),
                                    isvcon.getSaltConcentration(), ui.getExceptionBehavior());
  }

  // Implicit and explicit solvent simulations utilize different exclusion masks.  Declare each as
  // a blank object and then populate the relevant one.
  StaticExclusionMaskSynthesis poly_se;
  LocalExclusionMask lem;
  switch (poly_ag.getUnitCellType()) {
  case UnitCellType::NONE:
    poly_se = createMaskSynthesis(sc, poly_ag);
    dyna_tk.setExclusionMasks(poly_se);
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    lem = LocalExclusionMask(poly_ag);
    dyna_tk.setExclusionMasks(lem);
    break;
  }
  
  // The synthesis cache map can be laid out as an incrementing series (0, 1, 2, ..., k) over all
  // systems in the synthesis, owing to the way that the SystemCache was used to produce the
  // synthesis.  Any directives from user input to replicate certain systems, e.g.
  // "-sys { -p <topology> -c <input coordinates> -n 50 }", was performed in the construction of
  // the SystemCache itself.  Calls to the ith system of the synthesis will go to the ith member
  // of the SystemCache in a 1:1 mapping, in this case.
  SynthesisCacheMap scmap(incrementingSeries(0, poly_ag.getSystemCount()), &sc, &poly_ag,
                          &poly_ps);

  // With complete topology and coordinate syntheses, the dynamics intervention can be loaded with
  // the basic resources.
  dyna_tk.setCoordinateSynthesis(&poly_ps);
  dyna_tk.setTopologySynthesis(&poly_ag);
  
  // Prepare a debugging ledger in the form of a Watcher class object.  Attach the object
  // immediately to the global class object used in interventions.  An empty array is initialized
  // to conserve memory if no Watcher is needed, but also to keep the Watcher object from falling
  // out of scope if it is created in the conditional branch below.
  const DebugControls& dbgcon = ui.getDebugNamelistInfo();
  std::vector<Watcher> monitor;
  std::vector<PhaseSpaceSynthesis> coordinate_workspaces;
  if (ui.getDebugPresence()) {
    monitor.reserve(1);
    monitor.emplace_back(poly_ps, poly_ag, dbgcon);
    coordinate_workspaces.reserve(1);
#ifdef STORMM_USE_HPC
    coordinate_workspaces.emplace_back(poly_ps, HybridFormat::HOST_MOUNTED);
#else
    coordinate_workspaces.emplace_back(poly_ps, HybridFormat::HOST_ONLY);
#endif
    dyna_tk.setAnomalyReporting(&monitor[0]);
    dyna_tk.setDebugging(&monitor[0], dbgcon, dyncon, &coordinate_workspaces[0], operating_tier);
  }
  
  // Check for analyses and create the objects as appropriate.
  const std::vector<AnalysisControls>& trkcon_v = ui.getAnalysisNamelistInfo();
  std::vector<HydrogenBondAnalysis> hbtrack_v;
  int n_hbtrack = 0;
  for (size_t i = 0; i < trkcon_v.size(); i++) {
    n_hbtrack += (trkcon_v[i].getHBondMaskCount() > 0);
  }
  hbtrack_v.reserve(n_hbtrack);
  for (size_t i = 0; i < trkcon_v.size(); i++) {
    hbtrack_v.emplace_back(poly_ps, poly_ag, scmap, trkcon_v[i], dyncon);
  }
  for (size_t i = 0; i < hbtrack_v.size(); i++) {
    if (hbtrack_v[i].getTotalCandidates() > 0) {
#ifdef STORMM_USE_HPC
      hbtrack_v[i].uploadDefinitions();
      dyna_tk.setAnalysis(&hbtrack_v[i], dyncon, HybridTargetLevel::DEVICE);
#else
      dyna_tk.setAnalysis(&hbtrack_v[i], dyncon);
#endif
    }
  }
  
  // Initialization of the REMD control apparatus
  if (ui.getRemdPresence()) {
    RemdControls remdcon = ui.getRemdNamelistInfo();
    int total_swap_count = remdcon.getTotalSwapCount();
    std::string remd_type = remdcon.getRemdType();
    int frequency_swaps_count = remdcon.getFrequencyOfSwaps();
    std::string swap_storage = remdcon.getSwapStore();
    std::string temp_distribution = remdcon.getTemperatureDistributionMethod();
    double exchange_probability = remdcon.getExchangeProbability();
    double tolerance = remdcon.getTolerance();
    int max_replicas = remdcon.getMaxReplicas();
    double low_temperature = remdcon.getLowTemperature();
    double high_temperature = remdcon.getHighTemperature();
    ExchangeNexus remd_a(system_count, total_swap_count, remd_type, frequency_swaps_count,
                         swap_storage, temp_distribution, exchange_probability, tolerance,
                         max_replicas, low_temperature, high_temperature);
    remd_a.setAtomGraphSynthesis(&poly_ag);
    std::vector<double> t_replicas = remd_a.getTempDistribution();
    std::vector<int> remd_top(t_replicas.size());
    for(int i = 0; i < remd_top.size(); i++) {
      remd_top[i] = 0;
    }
    scmap = SynthesisCacheMap(remd_top, &sc, &poly_ag, &poly_ps);
  }

  // Set the progress bar to system_count
  ProgressBar progress_bar;
  const ReportControls repcon = ui.getReportNamelistInfo();
#ifndef STORMM_USE_HPC
  progress_bar.setCycleCount(system_count);
#else
  progress_bar.setCycleCount(1);
#endif
  progress_bar.setStyle(repcon.getProgressBarStyle());
  timer.assignTime(gen_setup_tm);

  // Perform minimizations as requested.
  if (ui.getMinimizePresence()) {
    const MinimizeControls mincon = ui.getMinimizeNamelistInfo();
#ifdef STORMM_USE_HPC
    switch (poly_ag.getUnitCellType()) {
    case UnitCellType::NONE:
      {
        // Isolated boundary conditions involve all-to-all interactions and open the door to
        // implicit solvent models.  However, the non-bonded work units for minimizations in such
        // a case are not the same as those for dynamics.
        InitializationTask ism_prep;
        const SolventControls& isvcon = ui.getSolventNamelistInfo();
        switch (isvcon.getImplicitSolventModel()) {
        case ImplicitSolventModel::NONE:
          ism_prep = InitializationTask::GENERAL_MINIMIZATION;
          break;
        case ImplicitSolventModel::HCT_GB:
        case ImplicitSolventModel::OBC_GB:
        case ImplicitSolventModel::OBC_GB_II:
        case ImplicitSolventModel::NECK_GB:
        case ImplicitSolventModel::NECK_GB_II:
          ism_prep = InitializationTask::GB_MINIMIZATION;
          break;
        }
        poly_ag.loadNonbondedWorkUnits(poly_se, ism_prep, 0, gpu);
      }
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      rtErr("Minimization is not yet operational for periodic boundary conditions.", "main");
    }

    // Upload data to prepare for energy minimizations
    poly_ps.upload();
    poly_ag.upload();
    switch (poly_ag.getUnitCellType()) {
    case UnitCellType::NONE:
      poly_se.upload();
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      break;
    }
    timer.assignTime(min_setup_tm);
    
    // Perform energy minimization, with additional branching over the unit cell type.
    ScoreCard emin;
    switch (poly_ag.getUnitCellType()) {
    case UnitCellType::NONE:
      emin = launchMinimization(poly_ag, poly_se, &poly_ps, mincon, gpu,
                                preccon.getValenceMethod(), preccon.getEnergyScalingBits(),
                                &timer, &progress_bar, "Run, minimization");
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      rtErr("Minimization is not yet operational for periodic boundary conditions.", "main");
    }
    emin.computeTotalEnergy(HybridTargetLevel::DEVICE, gpu);
    
    // Download the energies and also coordinates, to prepare for a CPU-based velocity seeding.
    emin.download();
    poly_ps.download();
    timer.assignTime(download_tm);
#else
    std::vector<ScoreCard> all_mme;
    all_mme.reserve(system_count);

    // Print out the stage for progress bar, reset bar
    progress_bar.setTitle("Run, minimization");
    progress_bar.setCycleCount(system_count);
    progress_bar.reset();
    switch (poly_ps.getUnitCellType()) {
    case UnitCellType::NONE:
    
      // Loop over all systems 
      for (int i = 0; i < system_count; i++) {
        const int icache_top = scmap.getTopologyCacheIndex(i);
        PhaseSpace ps = poly_ps.exportSystem(i);
        AtomGraph *ag = sc.getTopologyPointer(icache_top);
        const RestraintApparatus& ra = sc.getRestraints(i);
        all_mme.emplace_back(minimize(&ps, *ag, ra, sc.getSystemStaticMask(i), mincon));
        poly_ps.importSystem(ps, i);
        progress_bar.update();
      }
      timer.assignTime(min_run_tm);
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      rtErr("Minimization is not yet operational for periodic boundary conditions.", "main");
    }
    progress_bar.finalizeTerminalOutput();
#endif

    // Check the results at this stage: alert the user if any of them have high bond strain,
    // and prune the synthesis of problematic systems.
#if 0
    if (ui.getDebugPresence() && dbgcon.runSanityChecks()) {
      const std::vector<int> survivors = checkBondAndAngleSanity(poly_ps, dbgcon);
      if (survivors.size() < poly_ps.getSystemCount()) {
        AtomGraphSynthesis new_poly_ag = cullSynthesis(poly_ag, survivors, gpu);
        PhaseSpaceSynthesis new_poly_ps = cullSynthesis(poly_ps, survivors, gpu);
        SynthesisCacheMap new_scmap = remapSynthesis(scmap, survivors);
      }
    }
#endif    
    // Print restart files from energy minimization
    if (mincon.getCheckpointProduction()) {
      progress_bar.reset();
      progress_bar.setTitle("Checkpoint files ");
      progress_bar.setCycleCount(system_count);
      for (int i = 0; i < system_count; i++) {
        const PhaseSpace ps = poly_ps.exportSystem(i);
        ps.exportToFile(sc.getCheckpointName(i), 0.0, TrajectoryKind::POSITIONS,
                        CoordinateFileKind::AMBER_ASCII_RST, ui.getPrintingPolicy(),
                        repcon.getAsciiSalvageStyle());
        progress_bar.update();
      }
      progress_bar.finalizeTerminalOutput();
      timer.assignTime(output_tm);
    }
  }
  
  // Initialize trajectories, or stop if the trajectories already exist and overwriting is not
  // enabled.
  if (dyncon.getTrajectoryPrintFrequency() > 0) {
    std::vector<int> system_indices(1);
    switch (ui.getPrintingPolicy()) {
    case PrintSituation::OVERWRITE:
    case PrintSituation::OPEN_NEW:
      for (int i = 0; i < poly_ps.getSystemCount(); i++) {
        system_indices[0] = i;
        const int sysc_idx = scmap.getSystemCacheIndex(i);
        const std::string& traj_name = sc.getTrajectoryName(sysc_idx);
        initializeTrajectory(traj_name, ui.getPrintingPolicy(), sc.getTrajectoryKind(sysc_idx),
                             poly_ps.getAtomCount(i), 0.0);
      }
      break;
    case PrintSituation::APPEND:
    case PrintSituation::UNKNOWN:
      break;
    }
  }
  Thermostat tst(poly_ag, dyncon, sc, incrementingSeries(0, sc.getSystemCount()), gpu);
  
  // Kick-start dynamics if necessary.  A CPU-based routine is used for this, as it will involve a
  // a great deal of code to get it working on the GPU.  This will modify the thermostat's random
  // state.  Create a dummy thermostat to ensure that the same random numbers are not used to seed
  // velocities and then perform a first stochastic velocity modification.
  DynamicsControls mod_dyncon = dyncon;
  mod_dyncon.setThermostatSeed(dyncon.getThermostatSeed() + 715829320);
  mod_dyncon.setThermostatKind("langevin");
  Thermostat kickstarter(poly_ag, mod_dyncon, sc, incrementingSeries(0, sc.getSystemCount()), gpu);
  velocityKickStart(&poly_ps, poly_ag, &kickstarter, mod_dyncon, preccon.getValenceMethod(),
                    EnforceExactTemperature::YES);
#ifndef STORMM_USE_HPC
  // In order to perform CPU-based dynamics on systems in implicit solvent, thermostats
  // must be created for each system and persist throughout the entire simulation.  If
  // thermostats are created anew for each epoch of the replica exchange, they would need
  // to be initiated with different random seeds each time.  The most efficient way is to
  // create these thermostats once and let them keep charging forward.
  std::vector<Thermostat> tst_vec;
  switch (poly_ag.getUnitCellType()) {
  case UnitCellType::NONE:
    tst_vec.reserve(system_count);
    for (int i = 0; i < system_count; i++) {
      const AtomGraph* ag = poly_ps.getSystemTopologyPointer(i);
      const int sys_cache_idx = scmap.getSystemCacheIndex(i);
      const std::string &sys_lbl = sc.getSystemLabel(sys_cache_idx);
      int tstat_idx = findStringInVector(dyncon.getThermostatLabels(), sys_lbl);
      if (tstat_idx == dyncon.getThermostatLabels().size()) {
        tstat_idx = findStringInVector(dyncon.getThermostatLabels(), "all");
      }
      const std::vector<double> &t_init_targets = dyncon.getInitialTemperatureTargets();
      const std::vector<double> &t_finl_targets = dyncon.getFinalTemperatureTargets();
      tst_vec.emplace_back(ag->getAtomCount(), dyncon.getThermostatKind(),
                           t_init_targets[tstat_idx], t_finl_targets[tstat_idx],
                           dyncon.getThermostatEvolutionStart(),
                           dyncon.getThermostatEvolutionEnd(), PrecisionModel::SINGLE,
                           dyncon.getThermostatSeed() + i);
      tst_vec.back().setGeometryConstraints(dyncon.constrainGeometry());
      tst_vec.back().setRattleTolerance(dyncon.getRattleTolerance());
      tst_vec.back().setRattleIterations(dyncon.getRattleIterations());
    }
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    break;
  }
#endif
  
  // Run dynamics
  if (ui.getDynamicsPresence()) {
    const DynamicsControls dyncon = ui.getDynamicsNamelistInfo();
#ifdef STORMM_USE_HPC
    switch (poly_ag.getUnitCellType()) {
    case UnitCellType::NONE:
      {
        // Isolated boundary conditions involve all-to-all interactions and open the door to
        // implicit solvent models.  However, the non-bonded work units for minimizations in such
        // a case are not the same as those for dynamics.
        InitializationTask ism_prep;
        const SolventControls& isvcon = ui.getSolventNamelistInfo();
        switch (isvcon.getImplicitSolventModel()) {
        case ImplicitSolventModel::NONE:
          switch (dyncon.getThermostatKind()) {
          case ThermostatKind::NONE:
          case ThermostatKind::BERENDSEN:
            ism_prep = InitializationTask::GENERAL_DYNAMICS;
            break;
          case ThermostatKind::ANDERSEN:
          case ThermostatKind::LANGEVIN:
            ism_prep = InitializationTask::LANGEVIN_DYNAMICS;
            break;
          }
          break;
        case ImplicitSolventModel::HCT_GB:
        case ImplicitSolventModel::OBC_GB:
        case ImplicitSolventModel::OBC_GB_II:
        case ImplicitSolventModel::NECK_GB:
        case ImplicitSolventModel::NECK_GB_II:
          switch (dyncon.getThermostatKind()) {
          case ThermostatKind::NONE:
          case ThermostatKind::BERENDSEN:
            ism_prep = InitializationTask::GB_DYNAMICS;
            break;
          case ThermostatKind::ANDERSEN:
          case ThermostatKind::LANGEVIN:
            ism_prep = InitializationTask::GB_LANGEVIN_DYNAMICS;
            break;
          }
          break;
        }

        // Build the static exclusion mask synthesis, if it has not been built already.  Otherwise,
        // load non-bonded work units appropriate for dynamics as opposed to energy minimization.
        // The difference lies in array initialization instructions assigned to each work unit.
        poly_ag.loadNonbondedWorkUnits(poly_se, ism_prep, dyncon.getThermostatCacheDepth(), gpu);
      }
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      break;
    }
    
    // Upload data to prepare for dynamics.  If energy minimizations were performed on the GPU, the
    // coordinates were downloaded afterward in order to do the velocity kick-start on the host.
    // The most current coordinates must therefore be uploaded to the GPU.
    tst.uploadPartitions();
    poly_ps.upload();
    poly_ag.upload();
    poly_se.upload();
    timer.assignTime(dyn_setup_tm);
    
    // Perform molecular dynamics
    ScoreCard edyn;
    switch (poly_ag.getUnitCellType()) {
    case UnitCellType::NONE:
      edyn = launchDynamics(poly_ag, poly_se, &tst, &poly_ps, dyncon, repcon, sc, scmap, gpu,
                            preccon.getValenceMethod(), preccon.getNonbondedMethod(),
                            preccon.getEnergyScalingBits(), &timer, &progress_bar,
                            "Run, dynamics    ");
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      edyn = launchDynamics(poly_ag, &poly_ps, dyncon, pmecon, preccon, repcon, sc, scmap, gpu,
                            &timer, &progress_bar, "Run, dynamics    ");
      break;
    }
    cudaDeviceSynchronize();
    timer.assignTime(dyn_run_tm);
    edyn.download();
    timer.assignTime(download_tm);
#else // STORMM_USE_HPC
    const int nstep = dyncon.getStepCount();
    const int ntpr  = dyncon.getDiagnosticPrintFrequency();
    ScoreCard edyn(system_count, ((nstep + ntpr - 1) / ntpr) + 1, preccon.getEnergyScalingBits());
    if (ui.getRemdPresence()) {
      const RemdControls& remdcon = ui.getRemdNamelistInfo();
      const int total_swap_count = remdcon.getTotalSwapCount();
      const std::string remd_type = remdcon.getRemdType();
      const int nstep = remdcon.getFrequencyOfSwaps();
      const std::string swap_storage = remdcon.getSwapStore();
      const std::string temp_distribution = remdcon.getTemperatureDistributionMethod();
      const double exchange_probability = remdcon.getExchangeProbability();
      const double tolerance = remdcon.getTolerance();
      const int max_replicas = remdcon.getMaxReplicas();
      const double low_temperature = remdcon.getLowTemperature();
      const double high_temperature = remdcon.getHighTemperature();
      const int ntpr = dyncon.getDiagnosticPrintFrequency();
      ExchangeNexus remd_a(system_count, total_swap_count, remd_type, nstep, swap_storage,
                           temp_distribution, exchange_probability, tolerance, max_replicas,
                           low_temperature, high_temperature);
      
      // Initialize ScoreCard for REMD
      edyn = ScoreCard(system_count, ((nstep + ntpr - 1) / ntpr) + 1,
                       preccon.getEnergyScalingBits());

      // Initialize progress bars for REMD
      ProgressBar remd_progress("Run, REMD        ", remdcon.getTotalSwapCount());
      
      // As elsewhere, branch based on the type of boundary conditions
      switch (poly_ps.getUnitCellType()) {
      case UnitCellType::NONE:
        {
          // REMD outer loop for swap attempts
          for (int epoch = 0; epoch < remdcon.getTotalSwapCount(); epoch++) {

            // Update REMD progress bar
            remd_progress.update();

            // Inner loop over each system replica
            for (int i = 0; i < system_count; i++) {

              // The topology index which the system is based on must be obtained from the
              // synthesis cache map.
              const int icache_top = scmap.getTopologyCacheIndex(i);

              // Export system state
              PhaseSpace ps = poly_ps.exportSystem(i);
              const AtomGraph* ag = poly_ag.getSystemTopologyPointer(i);
              const RestraintApparatus *ra = poly_ag.getSystemRestraintPointer(i);
              timer.assignTime(dyn_setup_tm);
              ScoreCard iedyn(1, ((nstep + ntpr - 1) / ntpr) + 1, preccon.getEnergyScalingBits());
              dynamics(&ps, &tst_vec[i], &iedyn, *ag, ngb_tab, sc.getSystemStaticMask(icache_top),
                       *ra, dyncon, 0);
              timer.assignTime(dyn_run_tm);

              // Import dynamics data into the main scorecard
              edyn.importCard(iedyn, i, 0);
              
              // Update progress bar for dynamics within each replica
              progress_bar.update();
            }
          }
        }
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        rtErr("Replica Exchange molecular dynamics is not yet operational for periodic boundary "
              "conditions.", "main");
        break;
      }
      remd_progress.finalizeTerminalOutput();
    }
    else {

      // Reset Progress bar before loop, print the purpose of the loop.  Technically, we do not
      // need to do a setCycleCount() here, since system_count is the same throughout.
      progress_bar.setTitle("Run, dynamics    ");
      progress_bar.setCycleCount(system_count);
      progress_bar.reset();
      switch (poly_ps.getUnitCellType()) {
      case UnitCellType::NONE:
        for (int i = 0; i < system_count; i++) {

          // The topology index which the system is based on must be obtained from the
          // synthesis cache map.
          const int icache_top = scmap.getTopologyCacheIndex(i);
          const int icache_sys = scmap.getSystemCacheIndex(i);
          const std::string itraj_name = sc.getTrajectoryName(icache_sys);

          // Unpack the individual system for CPU calculations.
          PhaseSpace ps = poly_ps.exportSystem(i);
          const AtomGraph *ag = poly_ps.getSystemTopologyPointer(i);
          const RestraintApparatus& ra = sc.getRestraints(i);
          timer.assignTime(dyn_setup_tm);
          ScoreCard iedyn(1, ((nstep + ntpr - 1) / ntpr) + 1, preccon.getEnergyScalingBits());
          dynamics(&ps, &tst_vec[i], &iedyn, *ag, ngb_tab, sc.getSystemStaticMask(icache_top), ra,
                   dyncon, 0, itraj_name);
          timer.assignTime(dyn_run_tm);
          edyn.importCard(iedyn, i, 0);
        
          // Update progress bar at the beginning of the loop
          progress_bar.update();
        }
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        switch (preccon.getNonbondedMethod()) {
        case PrecisionModel::DOUBLE:
          {
            CellGrid<double, llint, double, double4_16a> cg(poly_ps, poly_ag,
                                                            dyncon.getVanDerWaalsCutoff(), 0.02,
                                                            pmecon.getMeshSubdivisions(),
                                                            NonbondedTheme::ALL);
            cg.checkViability(scmap);
            dyna_tk.setNeighborList<double, llint, double, double4_16a>(&cg);
            PMIGrid pmig(&cg, NonbondedTheme::ELECTROSTATIC, pmecon.getInterpolationOrder(),
                         PrecisionModel::DOUBLE, FFTMode::OUT_OF_PLACE,
                         preccon.getChargeMeshScalingBits(), preccon.getChargeMeshScalingBits());
            ConvolutionManager cvol(&pmig, pmecon.getEwaldCoefficient());
            dynamics<double, llint, double, double4_16a>(&poly_ps, &cg, &pmig, &cvol, &edyn, &tst,
                                                         poly_ag, lem, dyncon, preccon, pmecon);
          }
          break;
        case PrecisionModel::SINGLE:
          {
            CellGrid<float, int, float, float4> cg(poly_ps, poly_ag,
                                                   dyncon.getVanDerWaalsCutoff(), 0.02,
                                                   pmecon.getMeshSubdivisions(),
                                                   NonbondedTheme::ALL);
            cg.checkViability(scmap);
            dyna_tk.setNeighborList<float, int, float, float4>(&cg);
            PMIGrid pmig(&cg, NonbondedTheme::ELECTROSTATIC, pmecon.getInterpolationOrder(),
                         PrecisionModel::SINGLE, FFTMode::OUT_OF_PLACE,
                         preccon.getChargeMeshScalingBits(), preccon.getChargeMeshScalingBits());
            ConvolutionManager cvol(&pmig, pmecon.getEwaldCoefficient());
            dynamics<float, int, float, float4>(&poly_ps, &cg, &pmig, &cvol, &edyn, &tst, poly_ag,
                                                lem, dyncon, preccon, pmecon);
          }
          break;
        }
      }
      progress_bar.finalizeTerminalOutput();
    }
#endif // STORMM_USE_HPC
    
    // Turn the energy tracking data into an output report
    createDiagnosticReport(edyn, scmap, ui);
    
    // Turn the analyses and debugging operations into their own output reports.
#ifdef STORMM_USE_HPC
    dyna_tk.downloadDebugging();
    dyna_tk.downloadAnalysisData();
#endif
    createDebugReport(ui);
    createAnalysisReport(ui);
    
    // Print restart files from dynamics
#ifdef STORMM_USE_HPC
    poly_ps.download();
#endif
    // Reset the progress bar before writing all checkpoint files.
    progress_bar.reset();
    progress_bar.setTitle("Checkpoint files ");
    progress_bar.setCycleCount(system_count);
    for (int i = 0; i < system_count; i++) {
      const PhaseSpace ps = poly_ps.exportSystem(i);
      ps.exportToFile(sc.getCheckpointName(i), 0.0, TrajectoryKind::POSITIONS,
                      CoordinateFileKind::AMBER_ASCII_RST, ui.getPrintingPolicy(),
                      repcon.getAsciiSalvageStyle());

      // Update progress bar at the beginning of the loop
      progress_bar.update();
    }
    progress_bar.finalizeTerminalOutput();
    
    // At the end of the progress bar, endl for the rest of the program
    std::cout << std::endl;
    timer.assignTime(output_tm);
  }
  
  // Summarize the results
  if (repcon.printWallTimeData()) {
    timer.assignTime(output_tm);
    timer.printResults();
  }

  return 0;
}
