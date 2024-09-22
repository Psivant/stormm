#include <vector>
#  include "../../../src/Accelerator/gpu_details.h"
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
#include "../../../src/Math/vector_ops.h"
#include "../../../src/MolecularMechanics/minimization.h"
#include "../../../src/Namelists/input_transcript.h"
#include "../../../src/Namelists/nml_dynamics.h"
#include "../../../src/Namelists/nml_minimize.h"
#include "../../../src/Namelists/nml_precision.h"
#include "../../../src/Namelists/nml_random.h"
#include "../../../src/Namelists/user_settings.h"
#include "../../../src/Potential/energy_enumerators.h"
#include "../../../src/Reporting/error_format.h"
#include "../../../src/Reporting/help_messages.h"
#include "../../../src/Reporting/present_energy.h"
#include "../../../src/Reporting/progress_bar.h"
#include "../../../src/Synthesis/atomgraph_synthesis.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/static_mask_synthesis.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/Topology/atomgraph_enumerators.h"
#include "../../../src/Trajectory/phasespace.h"
#include "../../../src/Trajectory/thermostat.h"
#include "../../../src/UnitTesting/stopwatch.h"
#include "../../../src/UnitTesting/unit_test.h"
#include "setup.h"

using namespace stormm::card;
using namespace stormm::chemistry;
using namespace stormm::display;
using namespace stormm::energy;
using namespace stormm::mm;
using namespace stormm::namelist;
using namespace stormm::random;
using namespace stormm::reporting;
using namespace stormm::restraints;
using namespace stormm::review;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;
using namespace dyna_app::setup;

//-------------------------------------------------------------------------------------------------
// Display a general help message for this program.
//-------------------------------------------------------------------------------------------------
void displayGeneralHelpMessage(const std::vector<std::string> &all_namelists) {
  const std::string base_msg =
    terminalFormat("A program for carrying out dynamics on molecular systems using standard "
                   "coordinates and topologies.\n\n", "dynamics");
  std::string list_of_nml;
  for (size_t i = 0; i < all_namelists.size(); i++) {
    list_of_nml += "  - " + all_namelists[i] + "\n";
  }
  const std::string nml_msg =
    terminalFormat("Applicable namelists (re-run with one of these terms as the command-line "
                   "argument, IN QUOTES, i.e. \"&files\" or '&files', for further details):\n" +
                   list_of_nml, nullptr, nullptr, 0, 2, 2);
  printf("%s", base_msg.c_str());
  printf("%s", nml_msg.c_str());
  printf("\n");
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Check for a help message
  const std::vector<std::string> my_namelists = {
    "&files", "&dynamics", "&restraint", "&solvent", "&random", "&minimize", "&report",
    "&precision"
  };
  if (detectHelpSignal(argc, argv)) {
    displayGeneralHelpMessage(my_namelists);
    return 0;
  }
  if (displayNamelistHelp(argc, argv, my_namelists)) {
    return 0;
  }
  
  // Wall time tracking
  StopWatch timer("Timings for dynamics.stormm");
  const int file_parse_tm = timer.addCategory("File Parsing");
  const int gen_setup_tm  = timer.addCategory("Setup, General");
  const int min_setup_tm  = timer.addCategory("Setup, Minimization");
  const int dyn_setup_tm  = timer.addCategory("Setup, Dynamics");
  const int min_run_tm    = timer.addCategory("Run, Minimization");
  const int dyn_run_tm    = timer.addCategory("Run, Dynamics");
  const int download_tm   = timer.addCategory("GPU Data Download");
  const int output_tm     = timer.addCategory("Trajectory Output");
  
  // Engage the GPU
#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  const Hybrid<int> array_to_trigger_gpu_mapping(1);
#else
  const GpuDetails gpu = null_gpu;
#endif
  timer.assignTime(gen_setup_tm);

  // Read information from the command line and initialize the UserSettings object
  UserSettings ui(argc, argv, AppName::DYNAMICS);
  
  // Read topologies and coordinate files.  Assemble critical deatils about each system.
  SystemCache sc(ui.getFilesNamelistInfo(), ui.getExceptionBehavior(), MapRotatableGroups::NO,
                 ui.getPrintingPolicy());
  timer.assignTime(file_parse_tm);

  // Prepare a synthesis of systems from the user input.
  const std::vector<AtomGraph*> agv = sc.getTopologyPointer();

  // Create the synthesis of systems, including exclusion masks and non-bonded work units as
  // necessary.
  const int system_count = sc.getSystemCount();
  const PrecisionControls& preccon = ui.getPrecisionNamelistInfo();
  const DynamicsControls& dyncon = ui.getDynamicsNamelistInfo();
  const ThermostatKind tstat_choice = dyncon.getThermostatKind();
  PhaseSpaceSynthesis poly_ps = sc.exportCoordinateSynthesis(preccon.getGlobalPosScalingBits(),
                                                             preccon.getVelocityScalingBits(),
                                                             preccon.getForceScalingBits());
  AtomGraphSynthesis poly_ag = sc.exportTopologySynthesis(gpu, ui.getExceptionBehavior());
  StaticExclusionMaskSynthesis poly_se = createMaskSynthesis(sc);
  SynthesisCacheMap scmap(incrementingSeries(0, poly_ps.getSystemCount()), &sc, &poly_ag,
                          &poly_ps);

  // Set the progress bar to system_count
  ProgressBar progress_bar;
  progress_bar.initialize(system_count);
  // Preview the implicit solvent model.
  const SolventControls& isvcon = ui.getSolventNamelistInfo();
  NeckGeneralizedBornTable ngb_tab;
  poly_ag.setImplicitSolventModel(isvcon.getImplicitSolventModel(), ngb_tab,
                                  isvcon.getPBRadiiSet(), isvcon.getExternalDielectric(),
                                  isvcon.getSaltConcentration(), ui.getExceptionBehavior());
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
      rtErr("Minimization and dynamics are not yet operational for periodic boundary conditions.",
            "main");
    }

    // Upload data to prepare for energy minimizations
    poly_ps.upload();
    poly_ag.upload();
    poly_se.upload();
    timer.assignTime(min_setup_tm);
    
    // Perform energy minimization
    ScoreCard emin = launchMinimization(poly_ag, poly_se, &poly_ps, mincon, gpu,
                                        preccon.getValenceMethod(),
                                        preccon.getEnergyScalingBits());
    emin.computeTotalEnergy(HybridTargetLevel::DEVICE, gpu);
    cudaDeviceSynchronize();
    timer.assignTime(min_run_tm);

    // Download the energies and also coordinates, to prepare for a CPU-based velocity seeding.
    emin.download();
    poly_ps.download();
    timer.assignTime(download_tm);
#else
    std::vector<ScoreCard> all_mme;
    all_mme.reserve(system_count);
    // Print out the stage for progress bar, reset bar
    std::cout << "Minimization" << std::endl;
    progress_bar.reset();
    for (int i = 0; i < system_count; i++) {
      progress_bar.update(); // Update progress bar.
      PhaseSpace ps = poly_ps.exportSystem(i);
      AtomGraph *ag = sc.getSystemTopologyPointer(i);
      ag->setImplicitSolventModel(isvcon.getImplicitSolventModel());
      const RestraintApparatus& ra = sc.getRestraints(i);
      switch(ps.getUnitCellType()) {
      case UnitCellType::NONE:
        all_mme.emplace_back(minimize(&ps, *ag, ra, sc.getSystemStaticMask(i), mincon));
        timer.assignTime(min_run_tm);
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        rtErr("Minimization and dynamics are not yet operational for periodic boundary "
              "conditions.", "main");
      }
      poly_ps.import(ps, i);
    }
    // After end of loop, endl for the rest of the program to print correctly
    std::cout << std::endl;
#endif
    // Print restart files from energy minimization
    if (mincon.getCheckpointProduction()) {
      for (int i = 0; i < system_count; i++) {
        const PhaseSpace ps = poly_ps.exportSystem(i);
        ps.exportToFile(sc.getCheckpointName(i), 0.0, TrajectoryKind::POSITIONS,
                        CoordinateFileKind::AMBER_ASCII_RST, ui.getPrintingPolicy());
      }
      timer.assignTime(output_tm);
    }
  }

  // Kick-start dynamics if necessary.  A CPU-based routine is used for this, as it will involve a
  // a great deal of code to get it working on the GPU.  This will modify the thermostat's random
  // state.  Create a dummy thermostat to ensure that the same random numbers are not used to seed
  // velocities and then perform a first stochastic velocity modification.
  DynamicsControls mod_dyncon = dyncon;
  mod_dyncon.setThermostatSeed(dyncon.getThermostatSeed() + 715829320);
  mod_dyncon.setThermostatKind("langevin");
  Thermostat tst(poly_ag, mod_dyncon, sc, incrementingSeries(0, sc.getSystemCount()), gpu);
  velocityKickStart(&poly_ps, poly_ag, &tst, mod_dyncon, preccon.getValenceMethod(),
                    EnforceExactTemperature::YES);
  
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
      rtErr("Minimization and dynamics are not yet operational for periodic boundary conditions.",
            "main");
    }

    // Upload data to prepare for dynamics.  If energy minimizations were performed on the GPU, the
    // coordinates were downloaded afterward in order to do the velocity kick-start on the host.
    // The most current coordinates must therefore be uploaded to the GPU.
    tst = Thermostat(poly_ag, dyncon, sc, incrementingSeries(0, sc.getSystemCount()), gpu);
    tst.uploadPartitions();
    poly_ps.upload();
    poly_ag.upload();
    poly_se.upload();
    timer.assignTime(dyn_setup_tm);

    // Perform molecular dynamics
    ScoreCard edyn = launchDynamics(poly_ag, poly_se, &tst, &poly_ps, dyncon, sc, scmap, gpu,
                                    preccon.getValenceMethod(), preccon.getNonbondedMethod(),
                                    preccon.getEnergyScalingBits(), &timer, "main_dynamics");
    cudaDeviceSynchronize();
    timer.assignTime(dyn_run_tm);
    edyn.download();
    timer.assignTime(download_tm);
#else
    const int nstep = dyncon.getStepCount();
    const int ntpr  = dyncon.getDiagnosticPrintFrequency();
    ScoreCard edyn(system_count, ((nstep + ntpr - 1) / ntpr) + 1, preccon.getEnergyScalingBits());

    // Reset Progress bar before loop, print the purpose of the loop.  Technically, we do not need
    // to do a setIterations() here, since system_count is the same throughout.
    progress_bar.reset();
    progress_bar.setIterations(system_count);
    std::cout << "Dynamics" << std::endl;

    for (int i = 0; i < system_count; i++) {
      
      // Update progress bar at the beginning of the loop
      progress_bar.update();
      PhaseSpace ps = poly_ps.exportSystem(i);
      AtomGraph *ag = sc.getSystemTopologyPointer(i);
      Thermostat itst(ag->getAtomCount(), dyncon.getThermostatKind(), 298.15, 298.15,
                      dyncon.getThermostatEvolutionStart(), dyncon.getThermostatEvolutionEnd());
      itst.setGeometryConstraints(dyncon.constrainGeometry());
      itst.setRattleTolerance(dyncon.getRattleTolerance());
      itst.setRattleIterations(dyncon.getRattleIterations());

      // Set the implicit solvent model here if it has not been done already.
      if (ui.getMinimizePresence() == false) {
        ag->setImplicitSolventModel(isvcon.getImplicitSolventModel());
      }
      const RestraintApparatus& ra = sc.getRestraints(i);
      timer.assignTime(dyn_setup_tm);
      ScoreCard iedyn(1, ((nstep + ntpr - 1) / ntpr) + 1, preccon.getEnergyScalingBits());
      switch(ps.getUnitCellType()) {
      case UnitCellType::NONE:
        dynamics(&ps, &itst, &iedyn, *ag, ngb_tab, sc.getSystemStaticMask(i), ra, dyncon, 0);
        timer.assignTime(dyn_run_tm);
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        rtErr("Minimization and dynamics are not yet operational for periodic boundary "
              "conditions.", "main");
      }
      edyn.import(iedyn, i, 0);
    }
    // At the end of the progress bar, endl for the rest of the program
    std::cout << std::endl;
#endif

    // Turn the energy tracking data into an output report
    createDiagnosticReport(edyn, scmap, ui);

    // Print restart files from dynamics
#ifdef STORMM_USE_HPC
    poly_ps.download();
#endif
    // Progress Bar reset before loop, and purpose
    std::cout << "Exporting data" << std::endl;
    progress_bar.reset();
    for (int i = 0; i < system_count; i++) {
      progress_bar.update(); // Update progress bar at the beginning of the loop
      const PhaseSpace ps = poly_ps.exportSystem(i);
      ps.exportToFile(sc.getCheckpointName(i), 0.0, TrajectoryKind::POSITIONS,
                      CoordinateFileKind::AMBER_ASCII_RST, ui.getPrintingPolicy());
    }
    // At the end of the progress bar, endl for the rest of the program
    std::cout << std::endl;
    timer.assignTime(output_tm);
  }
  // Summarize the results
  timer.assignTime(output_tm);
  timer.printResults();
  
  return 0;
}
