#include <vector>
#include "copyright.h"
#include "../../../src/Accelerator/gpu_details.h"
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
#include "../../../src/DataTypes/common_types.h"
#include "../../../src/DataTypes/stormm_vector_types.h"
#include "../../../src/Math/vector_ops.h"
#include "../../../src/MolecularMechanics/minimization.h"
#include "../../../src/Namelists/input_transcript.h"
#include "../../../src/Namelists/nml_dynamics.h"
#include "../../../src/Namelists/nml_minimize.h"
#include "../../../src/Namelists/nml_precision.h"
#include "../../../src/Namelists/nml_random.h"
#include "../../../src/Namelists/user_settings.h"
#include "../../../src/Potential/energy_enumerators.h"
#include "../../../src/Potential/local_exclusionmask.h"
#include "../../../src/Potential/static_exclusionmask.h"
#include "../../../src/Reporting/error_format.h"
#include "../../../src/Reporting/help_messages.h"

using namespace stormm::card;
using namespace stormm::chemistry;
using namespace stormm::data_types;
using namespace stormm::display;
using namespace stormm::energy;
using namespace stormm::mm;
using namespace stormm::namelist;
using namespace stormm::random;
using namespace stormm::reporting;
using namespace stormm::review;
using namespace stormm::synthesis;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Wall time tracking
  Stopwatch timer("Timings for permability.stormm.cuda");
  const int file_parse_tm = timer.addCategory("File Parsing");
  const int gen_setup_tm  = timer.addCategory("Setup, General");
  const int min_setup_tm  = timer.addCategory("Setup, Minimization");
  const int dyn_setup_tm  = timer.addCategory("Setup, Dynamics");
  const int min_run_tm    = timer.addCategory("Run, Minimization");
  const int dyn_run_tm    = timer.addCategory("Run, Dynamics");
  const int analysis_tm   = timer.addCategory("Trajectory Analysis");
  const int regression_tm = timer.addCategory("Model Training");
  const int evaluation_tm = timer.addCategory("Model Evaluation");

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

  // Parse the command line
  CommandLineParser clip("permeability.stormm", "Application of molecular dynamics and molecular "
                         "property evaluation for design and evaluation of QSAR models of "
                         "membrane permability");
  const std::vector<std::string> lib_namelists = { "&files", "&minimize", "&dynamics", "&solvent",
                                                   "&random", "&report", "&precision" };
  clip.addControlBlocks(my_namelists);
  if (displayNamelistHelp(argc, argv, my_namelists) && clip.doesProgramExitOnHelp()) {
    return 0;
  }
  const std::vector<NamelistToken> perm_specific_namelists = {
    NamelistToken(std::string("&permeability"), permeabilityInput)
  };
  clip.addCustomNamelists(perm_specific_namelists);
  if (displayNamelistHelp(argc, argv, my_namelists, perm_specific_namelists) &&
      clip.doesProgramExitOnHelp()) {
    return 0;
  }
  clip.parseUserInput(argc, argv);

  // Read information from the command line and initialize the UserSettings object
  UserSettings ui(clip, { "-pe", "-ce" });

  // Read topologies and coordinate files.  Assemble critical deatils about each system.
  SystemCache sc(ui.getFilesNamelistInfo(), ui.getExceptionBehavior(), MapRotatableGroups::NO,
                 ui.getPrintingPolicy());
  timer.assignTime(file_parse_tm);

  // Prepare a synthesis of systems from the user input.
  const std::vector<AtomGraph*> agv = sc.getTopologyPointer();

  // Preview the implicit solvent model.  Set all systems to work with the high external
  // dielectric.
  NeckGeneralizedBornTable ngb_tab;
  if (ui.getSolventPresence()) {
    const SolventControls& isvcon = ui.getSolventNamelistInfo();
    for (int i = 0; i < sc.getTopologyCount(); i++) {
      AtomGraph *ag = sc.getTopologyPointer(i);
      ag->setImplicitSolventModel(isvcon.getImplicitSolventModel());
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
                                    isvcon.getPBRadiiSet(), permcon.getHighDielectric(),
                                    isvcon.getSaltConcentration(), ui.getExceptionBehavior());
  }
  StaticExclusionMaskSynthesis poly_se = createMaskSynthesis(sc, poly_ag);
  LocalExclusionMask lem(poly_ag);
  SynthesisCacheMap scmap(incrementingSeries(0, poly_ag.getSystemCount()), &sc, &poly_ag,
                          &poly_ps);

  // Conduct energy minimization, if requested.
  if (ui.getMinimizePresence()) {
    const MinimizeControls mincon = ui.getMinimizeNamelistInfo();

    // The unit cell type for each molecule is expected to be isolated boundary conditions.  An
    // implicit solvent model must be present.
    InitializationTask ism_prep;
    const SolventControls& isvcon = ui.getSolventNamelistInfo();
    poly_ag.loadNonbondedWorkUnits(poly_se, InitializationTask::GB_MINIMIZATION, 0, gpu);
  }
}
