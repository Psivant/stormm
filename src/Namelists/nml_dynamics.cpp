#include "copyright.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "namelist_element.h"
#include "nml_dynamics.h"

namespace stormm {
namespace namelist {

using constants::getEnumerationName;
using constants::translatePrecisionModel;
using parse::realToString;
using parse::NumberFormat;
using parse::strcmpCased;
using structure::translateApplyConstraints;
using structure::translateRattleMethod;
using trajectory::translateThermostatKind;
  
//-------------------------------------------------------------------------------------------------
DynamicsControls::DynamicsControls(const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    policy{policy_in},
    total_step_count{default_dynamics_nstlim},
    diagnostic_frequency{default_dynamics_ntpr},
    trajectory_frequency{default_dynamics_ntwx},
    time_step{default_dynamics_time_step},
    electrostatic_cutoff{default_electrostatic_cutoff},
    van_der_waals_cutoff{default_van_der_waals_cutoff},
    constrain_geometry{std::string(default_geometry_constraint_behavior)},
    rattle_tolerance{default_rattle_tolerance},
    rattle_iterations{default_rattle_max_iter},
    rattle_protocol{std::string(default_rattle_protocol)},
    thermostat_kind{std::string(default_thermostat_kind)},
    thermo_evolution_start{default_tstat_evo_window_start},
    thermo_evolution_end{default_tstat_evo_window_end},
    thermostat_cache_depth{default_thermostat_cache_depth},
    thermostat_seed{default_thermostat_random_seed},
    thermostat_cache_config{std::string(default_thermostat_cache_config)},
    andersen_frequency{default_andersen_frequency},
    langevin_frequency{default_langevin_frequency},
    nt_warp_multiplicity{default_nt_warp_multiplicity},
    initial_temperature_targets{}, final_temperature_targets{}, thermostat_labels{},
    thermostat_label_indices{}, thermostat_masks{},
    nml_transcript{"dynamics"}
{
  // Always initialize the thermostat groups with at least one entry comprising all possible
  // systems.  This will be overwritten if namelist input is taken.
  setThermostatGroup();
}

//-------------------------------------------------------------------------------------------------
DynamicsControls::DynamicsControls(const TextFile &tf, int *start_line, bool *found_nml,
                                   const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    DynamicsControls(policy_in)
{
  NamelistEmulator t_nml = dynamicsInput(tf, start_line, found_nml, policy, wrap);
  nml_transcript = t_nml;
  t_nml.assignVariable(&total_step_count, "nstlim");
  t_nml.assignVariable(&diagnostic_frequency, "ntpr");
  t_nml.assignVariable(&trajectory_frequency, "ntwx");
  t_nml.assignVariable(&com_motion_purge_ferquency, "nscm");
  t_nml.assignVariable(&time_step, "dt");
  t_nml.assignVariable(&electrostatic_cutoff, "elec_cut");
  t_nml.assignVariable(&van_der_waals_cutoff, "vdw_cut");
  t_nml.assignVariable(&electrostatic_cutoff, "cut");
  t_nml.assignVariable(&van_der_waals_cutoff, "cut");
  t_nml.assignVariable(&rattle_tolerance, "tol");
  t_nml.assignVariable(&rattle_iterations, "rattle_iter");
  setCpuRattleMethod(t_nml.getStringValue("rattle_style"));
  
  // Detect whether RATTLE and SETTLE should be activated
  if (t_nml.getKeywordStatus("rigid_geom") != InputStatus::MISSING) {
    setGeometricConstraints(t_nml.getStringValue("rigid_geom"));
  }
  
  // Detect the thermostat kind and assign global properties.  Then, fill out arrays of temperature
  // targets.  Detect label groups, system indices, and atom masks marking subgroups within each
  // named system.
  t_nml.assignVariable(&thermostat_kind, "ntt");
  t_nml.assignVariable(&thermo_evolution_start, "tevo_start");
  t_nml.assignVariable(&thermo_evolution_end, "tevo_end");
  t_nml.assignVariable(&thermostat_cache_depth, "tcache_depth");
  t_nml.assignVariable(&thermostat_seed, "thermostat_seed");
  t_nml.assignVariable(&thermostat_cache_config, "tcache_config");
  t_nml.assignVariable(&andersen_frequency, "vrand");
  t_nml.assignVariable(&langevin_frequency, "gamma_ln");
  t_nml.assignVariable(&nt_warp_multiplicity, "nt_mult");
  const int ntstat = t_nml.getKeywordEntries("temperature");
  initial_temperature_targets.resize(ntstat);
  final_temperature_targets.resize(ntstat);
  thermostat_labels.resize(ntstat);
  thermostat_label_indices.resize(ntstat);
  thermostat_masks.resize(ntstat);
  for (int i = 0; i < ntstat; i++) {
    t_nml.assignVariable(&initial_temperature_targets[i], "temperature", "tempi", i);
    t_nml.assignVariable(&final_temperature_targets[i], "temperature", "temp0", i);
    t_nml.assignVariable(&thermostat_labels[i], "temperature", "-label", i);
    t_nml.assignVariable(&thermostat_label_indices[i], "temperature", "-n", i);
    t_nml.assignVariable(&thermostat_masks[i], "temperature", "-mask", i);
  }
  
  // Validate input
  validateStepCount();
  validateDiagnosticPrintFrequency();
  validateTrajectoryPrintFrequency();
  validateCenterOfMassMotionPurgeFrequency();
  validateTimeStep();
  validateCutoffs();
  validateRattleTolerance();
  validateRattleIterations();
  validateThermostatKind();
  validateCacheConfiguration();
  validateNTWarpMultiplicity();
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getStepCount() const {
  return total_step_count;
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getDiagnosticPrintFrequency() const {
  return diagnostic_frequency;
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getTrajectoryPrintFrequency() const {
  return trajectory_frequency;
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getCenterOfMassMotionPurgeFrequency() const {
  return com_motion_purge_ferquency;
}

//-------------------------------------------------------------------------------------------------
double DynamicsControls::getTimeStep() const {
  return time_step;
}

//-------------------------------------------------------------------------------------------------
double DynamicsControls::getElectrostaticCutoff() const {
  return electrostatic_cutoff;
}

//-------------------------------------------------------------------------------------------------
double DynamicsControls::getVanDerWaalsCutoff() const {
  return van_der_waals_cutoff;
}

//-------------------------------------------------------------------------------------------------
ApplyConstraints DynamicsControls::constrainGeometry() const {
  return translateApplyConstraints(constrain_geometry);
}

//-------------------------------------------------------------------------------------------------
double DynamicsControls::getRattleTolerance() const {
  return rattle_tolerance;
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getRattleIterations() const {
  return rattle_iterations;
}

//-------------------------------------------------------------------------------------------------
RattleMethod DynamicsControls::getCpuRattleMethod() const {
  return translateRattleMethod(rattle_protocol);
}

//-------------------------------------------------------------------------------------------------
ThermostatKind DynamicsControls::getThermostatKind() const {
  return translateThermostatKind(thermostat_kind);
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getThermostatEvolutionStart() const {
  return thermo_evolution_start;
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getThermostatEvolutionEnd() const {
  return thermo_evolution_end;
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getThermostatCacheDepth() const {
  return thermostat_cache_depth;
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getThermostatSeed() const {
  return thermostat_seed;
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getThermostatLayerCount() const {
  const int nt = initial_temperature_targets.size();
  if (nt != final_temperature_targets.size() || nt != thermostat_labels.size() ||
      nt != thermostat_label_indices.size() || nt != thermostat_masks.size()) {
    rtErr("A consistent number of initial temperatures (" + std::to_string(nt) + "), final "
          "temperatures (" + std::to_string(final_temperature_targets.size()) + "), thermostat "
          "labels (" + std::to_string(thermostat_labels.size()) + "), label indices (" +
          std::to_string(thermostat_label_indices.size()) + "), and atom masks (" +
          std::to_string(thermostat_masks.size()) + ").", "DynamicsControls",
          "getThermosatLayerCount");
  }
  return nt;
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getAndersenFrequency() const {
  return andersen_frequency;
}

//-------------------------------------------------------------------------------------------------
double DynamicsControls::getLangevinFrequency() const {
  return langevin_frequency;
}

//-------------------------------------------------------------------------------------------------
PrecisionModel DynamicsControls::getThermostatCacheConfig() const {
  return translatePrecisionModel(thermostat_cache_config);
}

//-------------------------------------------------------------------------------------------------
const std::vector<double>& DynamicsControls::getInitialTemperatureTargets() const {
  return initial_temperature_targets;
}

//-------------------------------------------------------------------------------------------------
const std::vector<double>& DynamicsControls::getFinalTemperatureTargets() const {
  return final_temperature_targets;
}

//-------------------------------------------------------------------------------------------------
const std::vector<std::string>& DynamicsControls::getThermostatLabels() const {
  return thermostat_labels;
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& DynamicsControls::getThermostatLabelIndices() const {
  return thermostat_label_indices;
}

//-------------------------------------------------------------------------------------------------
const std::vector<std::string>& DynamicsControls::getThermostatMasks() const {
  return thermostat_masks;
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getNTWarpMultiplicity() const {
  return nt_warp_multiplicity;
}

//-------------------------------------------------------------------------------------------------
const NamelistEmulator& DynamicsControls::getTranscript() const {
  return nml_transcript;
}
  
//-------------------------------------------------------------------------------------------------
void DynamicsControls::setStepCount(const int total_step_count_in) {
  total_step_count = total_step_count_in;
  validateStepCount();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setDiagnosticPrintFrequency(const int diagnostic_frequency_in) {
  diagnostic_frequency = diagnostic_frequency_in;
  validateDiagnosticPrintFrequency();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setTrajectoryPrintFrequency(const int trajectory_frequency_in) {
  trajectory_frequency = trajectory_frequency_in;
  validateTrajectoryPrintFrequency();
}

//-------------------------------------------------------------------------------------------------
void
DynamicsControls::setCenterOfMassMotionPurgeFrequency(const int com_motion_purge_ferquency_in) {
  com_motion_purge_ferquency = com_motion_purge_ferquency_in;
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setTimeStep(const double time_step_in) {
  time_step = time_step_in;
  validateTimeStep();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setElectrostaticCutoff(const double cutoff_in) {
  electrostatic_cutoff = cutoff_in;
  validateCutoffs();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setVanDerWaalsCutoff(const double cutoff_in) {
  van_der_waals_cutoff = cutoff_in;
  validateCutoffs();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setCutoff(const double cutoff_in) {
  electrostatic_cutoff = cutoff_in;
  van_der_waals_cutoff = cutoff_in;
  validateCutoffs();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setGeometricConstraints(const std::string &constrain_geometry_in) {
  constrain_geometry = constrain_geometry_in;
  try {
    const ApplyConstraints trial = translateApplyConstraints(constrain_geometry_in);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Invalid constraint activation \"" + constrain_geometry_in + "\" provided.",
            "DynamicsControls", "setGeometricConstraints");
    case ExceptionResponse::WARN:
      rtWarn("Invalid constraint activation \"" + constrain_geometry_in + "\" provided.  The "
             "default of " + std::string(default_geometry_constraint_behavior) + " will be "
             "reinstated.", "DynamicsControls", "setGeometricConstraints");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    constrain_geometry = std::string(default_geometry_constraint_behavior);
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setGeometricConstraints(const ApplyConstraints constrain_geometry_in) {
  constrain_geometry = getEnumerationName(constrain_geometry_in);
}
  
//-------------------------------------------------------------------------------------------------
void DynamicsControls::setCpuRattleMethod(const std::string &rattle_protocol_in) {
  rattle_protocol = rattle_protocol_in;
  try {
    const RattleMethod interp = translateRattleMethod(rattle_protocol);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("An invalid RATTLE approach for CPU operations (" + rattle_protocol + ") was "
            "specified.", "DynamicsControls", "setCpuRattleMethod");
    case ExceptionResponse::WARN:
      rtWarn("An invalid RATTLE approach for CPU operations (" + rattle_protocol + ") was "
             "specified.  The default of " + std::string(default_rattle_protocol) + " will be "
             "restored.", "DynamicsControls", "setCpuRattleMethod");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    rattle_protocol = std::string(default_rattle_protocol);
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setRattleTolerance(const double rattle_tolerance_in) {
  rattle_tolerance = rattle_tolerance_in;
  validateRattleTolerance();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setThermostatKind(const std::string &thermostat_kind_in) {
  thermostat_kind = thermostat_kind_in;
  validateThermostatKind();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setThermostatKind(const ThermostatKind thermostat_kind_in) {
  thermostat_kind = getEnumerationName(thermostat_kind_in);
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setThermostatEvolutionStart(const int thermo_evolution_start_in) {
  thermo_evolution_start = thermo_evolution_start_in;
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setThermostatEvolutionEnd(const int thermo_evolution_end_in) {
  thermo_evolution_end = thermo_evolution_end_in;
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setThermostatCacheDepth(const int depth_in) {
  thermostat_cache_depth = depth_in;
  validateCacheDepth();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setThermostatSeed(const int igseed) {
  thermostat_seed = igseed;
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setAndersenFrequency(const int frequency_in) {
  andersen_frequency = frequency_in;
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setLangevinFrequency(const double frequency_in) {
  langevin_frequency = frequency_in;
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setThermostatCacheConfig(const std::string &cache_config_in) {
  thermostat_cache_config = cache_config_in;
  validateCacheConfiguration();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setThermostatCacheConfig(const PrecisionModel cache_config_in) {
  thermostat_cache_config = getEnumerationName(cache_config_in);
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setThermostatGroup(const double initial_target, const double final_target,
                                          const std::string &label, const int label_index,
                                          const std::string &mask) {
  validateTemperature(initial_target);
  validateTemperature(final_target);
  initial_temperature_targets.push_back(initial_target);
  final_temperature_targets.push_back(final_target);
  thermostat_labels.push_back(label);
  thermostat_label_indices.push_back(label_index);
  thermostat_masks.push_back(mask);
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setNTWarpMultiplicity(const int mult_in) {
  nt_warp_multiplicity = mult_in;
  validateNTWarpMultiplicity();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateStepCount() {
  if (total_step_count < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A negative value for the number of dynamics steps is invalid.  This error may be the "
            "result of trying to supply too large a number of steps (greater than 2.1 billion, "
            "2^31, which overflows the signed integer format.  Use checkpoint files to carry out "
            "runs with very large numbers of total steps in multiple segments.",
            "DynamicsControls", "validateStepCount");
    case ExceptionResponse::WARN:
      rtWarn("A negative value for the number of dynamics steps is invalid.  This error may be "
             "the result of trying to supply too large a number of steps (greater than 2.1 "
             "billion, 2^31, which overflows the signed integer format.  Use checkpoint files to "
             "carry out runs with very large numbers of total steps in multiple segments.  The "
             "default of " + std::to_string(default_dynamics_nstlim) + " will be restored.",
             "DynamicsControls", "validateStepCount");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    total_step_count = default_dynamics_nstlim;
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateDiagnosticPrintFrequency() {
  if (diagnostic_frequency < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A negative value of the diagnostic print frequency (" +
            std::to_string(diagnostic_frequency) + ") is invalid.", "DynamicsControls",
            "validateDiagnosticPrintFrequency");
    case ExceptionResponse::WARN:
      rtWarn("A negative value of the diagnostic print frequency (" +
             std::to_string(diagnostic_frequency) + ") is invalid.  The default of " +
             std::to_string(default_dynamics_ntpr) + " will be restored.", "DynamicsControls",
             "validateDiagnosticPrintFrequency");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    diagnostic_frequency = default_dynamics_ntpr;
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateTrajectoryPrintFrequency() {
  if (trajectory_frequency < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A negative value of the trajectory print frequency (" +
            std::to_string(trajectory_frequency) + ") is invalid.", "DynamicsControls",
            "validateTrajectoryPrintFrequency");
    case ExceptionResponse::WARN:
      rtWarn("A negative value of the trajectory print frequency (" +
             std::to_string(trajectory_frequency) + ") is invalid.  The default of " +
             std::to_string(default_dynamics_ntpr) + " will be restored.", "DynamicsControls",
             "validateTrajectoryPrintFrequency");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    trajectory_frequency = default_dynamics_ntwx;
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateCenterOfMassMotionPurgeFrequency() {
  if (com_motion_purge_ferquency < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A negative value of the center of mass motion purge frequency (" +
            std::to_string(com_motion_purge_ferquency) + ") is invalid.", "DynamicsControls",
            "validateCenterOfMassMotionPurgeFrequency");
    case ExceptionResponse::WARN:
      rtWarn("A negative value of the center of mass motion purge frequency (" +
             std::to_string(com_motion_purge_ferquency) + ") is invalid and will be adjusted to "
             "zero (no purging).", "DynamicsControls", "validateCenterOfMassMotionPurgeFrequency");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    com_motion_purge_ferquency = 0;
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateTimeStep() {
  if (time_step < minimum_dynamics_time_step) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A time step of " + realToString(time_step, 11, 4, NumberFormat::SCIENTIFIC) +
            "fs is too small for the dynamics to be accurate.  Use a time step larger than " +
            realToString( minimum_dynamics_time_step, 11, 4, NumberFormat::SCIENTIFIC) +
            "fs or a special command line option to ignore this input trap.", "DynamicsControls",
            "validateTimeStep");
    case ExceptionResponse::WARN:
      rtWarn("A time step of " + realToString(time_step, 11, 4, NumberFormat::SCIENTIFIC) +
             "fs is probably too small for the dynamics to be accurate.  The minimum step of " +
             realToString(minimum_dynamics_time_step, 11, 4, NumberFormat::SCIENTIFIC) + " will "
             "be taken instead.", "DynamicsControls", "validateTimeStep");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    time_step = minimum_dynamics_time_step;
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateCutoffs() {
  if (electrostatic_cutoff < minimum_elec_cutoff) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The electrostatic cutoff (" +
            realToString(electrostatic_cutoff, 9, 4, NumberFormat::STANDARD_REAL) + ") cannot be "
            "set below a minimum value of " +
            realToString(minimum_elec_cutoff, 9, 4, NumberFormat::STANDARD_REAL) + ".",
            "DynamicsControls", "validateCutoffs");
    case ExceptionResponse::WARN:
      rtErr("The electrostatic cutoff (" +
            realToString(electrostatic_cutoff, 9, 4, NumberFormat::STANDARD_REAL) + ") cannot be "
            "set below a minimum value of " +
            realToString(minimum_elec_cutoff, 9, 4, NumberFormat::STANDARD_REAL) + ".  This "
            "minimum value will be applied.", "DynamicsControls", "validateCutoffs");
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (van_der_waals_cutoff < minimum_vdw_cutoff) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A van-der Waals cutoff of " +
            realToString(van_der_waals_cutoff, 9, 4, NumberFormat::STANDARD_REAL) + " is too "
            "short to permit accurate force and energy computations.", "DynamicsControls",
            "validateCutoffs");
    case ExceptionResponse::WARN:
      rtWarn("A van-der Waals cutoff of " +
             realToString(van_der_waals_cutoff, 9, 4, NumberFormat::STANDARD_REAL) + " is too "
             "short to permit accurate force and energy computations.  The minimum value of " +
             realToString(minimum_vdw_cutoff, 9, 4, NumberFormat::STANDARD_REAL) + " will be "
             "applied.", "DynamicsControls", "validateCutoffs");
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateRattleTolerance() {
  if (rattle_tolerance < minimum_rattle_tolerance) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A tolerance of " + realToString(rattle_tolerance, 11, 4, NumberFormat::SCIENTIFIC) +
            " is likely to be unattainable.  Tolerances of less than " +
            realToString(minimum_rattle_tolerance, 11, 4, NumberFormat::SCIENTIFIC) + " are not "
            "likely to improve energy conservation any further.", "DynamicsControls",
            "validateRattleTolerance");
      break;
    case ExceptionResponse::WARN:
      rtWarn("A tolerance of " + realToString(rattle_tolerance, 11, 4, NumberFormat::SCIENTIFIC) +
            realToString(minimum_rattle_tolerance, 11, 4, NumberFormat::SCIENTIFIC) + " are not "
            "likely to improve energy conservation any further.  This tolerance will be taken "
             "instead.", "DynamicsControls", "validateRattleTolerance");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    rattle_tolerance = minimum_rattle_tolerance;    
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateRattleIterations() {
  if (rattle_iterations > maximum_rattle_iterations) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("At most " + std::to_string(maximum_rattle_iterations) + " may be attempted when "
            "constraining bond lengths (" + std::to_string(rattle_iterations) + " requested).",
            "DynamicsControls", "validateRattleIterations");
      break;
    case ExceptionResponse::WARN:
      rtErr("At most " + std::to_string(maximum_rattle_iterations) + " may be attempted when "
            "constraining bond lengths (" + std::to_string(rattle_iterations) + " requested).  "
            "The maximum iteration count will be taken.", "DynamicsControls",
            "validateRattleIterations");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    rattle_iterations = maximum_rattle_iterations;
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateTemperature(const double t) const {
  if (t < 0.0 || t > 10000.0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A temperature of " + realToString(t, 9, 4, NumberFormat::STANDARD_REAL) + " is "
            "unrealistic for molecular dynamics.", "DynamicsControls", "validateTemperature");
    case ExceptionResponse::WARN:
      rtWarn("A temperature of " + realToString(t, 9, 4, NumberFormat::STANDARD_REAL) + " is "
             "unrealistic for molecular dynamics.", "DynamicsControls", "validateTemperature");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateThermostatKind() {
  try {
    const ThermostatKind trial = translateThermostatKind(thermostat_kind);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("An invalid thermostat type " + thermostat_kind + " was detected in the input.",
            "DynamicsControls", "validateThermostatKind");
    case ExceptionResponse::WARN:
      rtWarn("An invalid thermostat type " + thermostat_kind + " was detected in the input.  It "
             "will be replaced with the default of " + std::string(default_thermostat_kind) + ".",
             "DynamicsControls", "validateThermostatKind");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    thermostat_kind = std::string(default_thermostat_kind);
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateCacheDepth() {
  switch (this->getThermostatKind()) {
  case ThermostatKind::NONE:
  case ThermostatKind::BERENDSEN:
    break;
  case ThermostatKind::ANDERSEN:
  case ThermostatKind::LANGEVIN:
    if (thermostat_cache_depth < 1) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("At least one layer of random numbers must be held in reserve for each atom when "
              "implementing a " + getEnumerationName(this->getThermostatKind()) + " thermostat.",
              "DynamicsControls", "validateCacheDepth");
      case ExceptionResponse::WARN:
        rtWarn("At least one layer of random numbers must be held in reserve for each atom when "
               "implementing a " + getEnumerationName(this->getThermostatKind()) + " thermostat.  "
               "The default of " + std::to_string(default_thermostat_cache_depth) + " will be "
               "reinstated.", "DynamicsControls", "validateCacheDepth");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
      thermostat_cache_depth = default_thermostat_cache_depth;
    }
    else if (thermostat_cache_depth > maximum_thermostat_cache_depth) {
      int storage_req;
      const PrecisionModel tc_prec = this->getThermostatCacheConfig();
      switch (tc_prec) {
      case PrecisionModel::DOUBLE:
        break;
        storage_req = 24 * thermostat_cache_depth;
      case PrecisionModel::SINGLE:
        storage_req = 12 * thermostat_cache_depth;
        break;
      }
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr(std::to_string(thermostat_cache_depth) + "layers of random numbers stored in " +
              getEnumerationName(tc_prec) + " require " + std::to_string(storage_req) + " bytes "
              "per atom.  This is unreasonably large and will not confer much optimization.",
              "DynamicsControls", "validateCacheDepth");
      case ExceptionResponse::WARN:
        rtWarn(std::to_string(thermostat_cache_depth) + "layers of random numbers stored in " +
               getEnumerationName(tc_prec) + " require " + std::to_string(storage_req) + " bytes "
               "per atom.  This is unreasonably large and will not confer much optimization "
               "beyond the maximum recommended value of " +
               std::to_string(maximum_thermostat_cache_depth) + ", which will be reinstated.",
               "DynamicsControls", "validateCacheDepth");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
      thermostat_cache_depth = maximum_thermostat_cache_depth;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateCacheConfiguration() {
  try {
    const PrecisionModel trial = translatePrecisionModel(thermostat_cache_config);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Invalid cache configuration " + thermostat_cache_config + ".",
            "DynamicsControls", "validateCacheConfiguration");
    case ExceptionResponse::WARN:
      rtWarn("Invalid cache configuration " + thermostat_cache_config + ".  The default "
             "of " + getEnumerationName(translatePrecisionModel(default_thermostat_cache_config)) +
             " will be reinstated.", "DynamicsControls", "validateCacheConfiguration");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    thermostat_cache_config = std::string(default_thermostat_cache_config);
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateNTWarpMultiplicity() {
  if (nt_warp_multiplicity < 0 || nt_warp_multiplicity > maximum_nt_warp_multiplicity) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The maximum number of warps that can be devoted to any one neutral territory "
            "decomposition is " + std::to_string(maximum_nt_warp_multiplicity) + ".",
            "DynamicsControls", "validateNTWarpMultiplicity");
    case ExceptionResponse::WARN:
      rtWarn("The maximum number of warps that can be devoted to any one neutral territory "
             "decomposition is " + std::to_string(maximum_nt_warp_multiplicity) + ".  The "
             "minimum value (1) will be applied.", "DynamicsControls",
             "validateNTWarpMultiplicity");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    nt_warp_multiplicity = 1;
  }
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator dynamicsInput(const TextFile &tf, int *start_line, bool *found,
                               const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("dynamics", CaseSensitivity::AUTOMATIC, policy, "Wraps directives needed "
                         "to propagate dynamics of a molecular system.");

  // Trajectory length and outut keywords
  t_nml.addKeyword("nstlim", NamelistType::INTEGER, std::to_string(default_dynamics_nstlim));
  t_nml.addKeyword("ntpr", NamelistType::INTEGER, std::to_string(default_dynamics_ntpr));
  t_nml.addKeyword("ntwx", NamelistType::INTEGER, std::to_string(default_dynamics_ntwx));
  t_nml.addKeyword("nscm", NamelistType::INTEGER, std::to_string(default_dynamics_nscm));
  t_nml.addKeyword("dt", NamelistType::REAL, std::to_string(default_dynamics_time_step));
  t_nml.addKeyword("elec_cut", NamelistType::REAL, std::to_string(default_electrostatic_cutoff));
  t_nml.addKeyword("vdw_cut", NamelistType::REAL, std::to_string(default_van_der_waals_cutoff));
  t_nml.addKeyword("cut", NamelistType::REAL, std::to_string(default_van_der_waals_cutoff));
  
  // Constraint keywords
  t_nml.addKeyword("rigid_geom", NamelistType::STRING);
  t_nml.addKeyword("tol", NamelistType::REAL,
                   realToString(default_rattle_tolerance, 9, 2, NumberFormat::SCIENTIFIC));
  t_nml.addKeyword("rattle_iter", NamelistType::INTEGER, std::to_string(default_rattle_max_iter));
  t_nml.addKeyword("rattle_style", NamelistType::STRING, std::string(default_rattle_protocol));

  // Thermostating keywords
  const NumberFormat std_real = NumberFormat::STANDARD_REAL;
  t_nml.addKeyword("ntt", NamelistType::STRING, std::string(default_thermostat_kind));
  t_nml.addKeyword("tevo_start", NamelistType::INTEGER,
                   std::to_string(default_tstat_evo_window_start));
  t_nml.addKeyword("tevo_end", NamelistType::INTEGER,
                   std::to_string(default_tstat_evo_window_end));
  t_nml.addKeyword("tcache_depth", NamelistType::INTEGER,
                   std::to_string(default_thermostat_cache_depth));
  t_nml.addKeyword("thermostat_seed", NamelistType::INTEGER,
                   std::to_string(default_thermostat_random_seed));
  t_nml.addKeyword("tcache_config", NamelistType::STRING,
                   std::string(default_thermostat_cache_config));
  t_nml.addKeyword("vrand", NamelistType::INTEGER, std::to_string(default_andersen_frequency));
  t_nml.addKeyword("gamma_ln", NamelistType::REAL, realToString(default_langevin_frequency, 9, 6,
                                                                NumberFormat::STANDARD_REAL));
  t_nml.addKeyword("nt_mult", NamelistType::INTEGER, std::to_string(default_nt_warp_multiplicity));
  const std::string tempr_help("Specify the temperature, or the temperature profile, at which "
                               "to maintain a system or group of particles within a system.");
  const std::vector<std::string> tempr_keys_help = {
    "Initial temperature to set at the outset of dynamics, and to maintain until the start of the "
    "evolution window (if the window is specified)",
    "Equilibrium temperature to maintain, after the completion of the evolution window (if such a "
    "window is specified)", "The system label to which this temperature profile applies",
    "The index within the label group to which this temperature profile applies.  The default "
    "value of " + std::to_string(-1) + " implies that all members of the label group will be "
    "affected.", "Atom mask for atoms in the system affected by this temperature profile"
  };
  t_nml.addKeyword("temperature", { "tempi", "temp0", "-label", "-n", "-mask" },
                   { NamelistType::REAL, NamelistType::REAL, NamelistType::STRING,
                     NamelistType::INTEGER, NamelistType::STRING },
                   { realToString(default_simulation_temperature, 9, 4, std_real),
                     realToString(default_simulation_temperature, 9, 4, std_real),
                     std::string("all"), std::to_string(-1), std::string("@=") },
                   DefaultIsObligatory::YES, InputRepeats::YES, tempr_help, tempr_keys_help,
                   { KeyRequirement::REQUIRED, KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                     KeyRequirement::REQUIRED, KeyRequirement::REQUIRED });
  
  // Help messages for each trajectory keyword
  t_nml.addHelp("nstlim", "Number of dynamics steps to carry out");
  t_nml.addHelp("ntpr", "The frequency with which to calculate and print energy diagnostics");
  t_nml.addHelp("ntwx", "The frequency with which to print trajectory snapshots");
  t_nml.addHelp("nscm", "The frequency with which to purge motion of the center of mass");
  t_nml.addHelp("dt", "The time step, in units of femtoseconds");

  // Help messages for interaction cutoffs
  t_nml.addHelp("elec_cut", "The inter-particle distance at which to begin discounting "
                "electrostatic interactions (this applies to all methods for evaluating the "
                "electrostatic potential.");
  t_nml.addHelp("vdw_cut", "The inter-particle distance at which to begin discounting "
                "van-der Waals interactions (this applies to all methods for evaluating the "
                "van-der Waals potential).");
  t_nml.addHelp("cut", "The inter-particle distance at which to begin neglecting pairwise, "
                "particle-particle interactions");
  
  // Help messages for geometry constraints keywords
  t_nml.addHelp("rigid_geom", "Indicate whether to enforce rigid geometries, namely bond length "
                "constraints and rigid water molecules");
  t_nml.addHelp("tol", "Tolerance by which to constrain rigid bonds involving hydrogen atoms");
  t_nml.addHelp("rattle_iter", "Maximum number of iterations to use when attempting to converge "
                "constrained bond lengths");
  t_nml.addHelp("rattle_style", "The manner in which to converge 'hub and spoke' constrained "
                "bonds connected to a single central atom");

  // Help messages for non-struct thermostat keywords
  t_nml.addHelp("ntt", "Indicate the type of thermostat, whether by a human-readable string or "
                "the Amber numeric cognate.  Options are case-insensitive and include \"none\" "
                "(0), \"berendsen\" (1), \"mass_andersen\" (2), or \"langevin\" (3).");
  t_nml.addHelp("tevo_start", "The step number at which to begin a linear evolution from the "
                "initial temperature to the equilibrium temperature");
  t_nml.addHelp("tevo_end", "The step number at which to conclude a linear evolution from the "
                "initial temperature to the equilibrium temperature");
  t_nml.addHelp("vrand", "Step count between complete velocity reassignments in the \"massive\" "
                "Andersen thermostat scheme.  This takes its name from the Amber keyword and "
                "pertains to a similar thermostat.");
  t_nml.addHelp("gamma_ln", "Collision frequency for a Langevin thermostat, in units of inverse "
                "femtoseconds.  These values should be 1/1000th of those fed to Amber, to comport "
                "with STORMM's internal time units of femtoseconds, not picoseconds.");
  t_nml.addHelp("nt_mult", "The number of warps that will cooperate to solve any given neutral "
                "territory decomposition subdomain.  By default, STORMM will look at the workload "
                "and the available GPU, then try to guess the best value between 1 and 8.  "
                "Maximum 8 warps per NT decomposition subdomain."); 
  t_nml.addHelp("tcache_depth", "Quantity of random numbers to produce for each atom's x, y, and "
                "z moves each time the random number generators are taken out of global memory "
                "and used.  Each use of a random generator increments it, and the result must be "
                "written back to global memory.  Each atom has its own generator state.  Creating "
                "and storing more sets of random numbers with each state checkout can thereby "
                "optimize memory traffic.  The default of 1 is good in most situations, but "
                "perhaps not all.");
  t_nml.addHelp("thermostat_seed", "Random number seed for the first random number state vector "
                "in the simulation.  In most situations, this is the state vector for the first "
                "atom.  Subsequent random number state vectors (guiding other variables, such as "
                "other atoms) will be initialized based on the long jump function for the "
                "XOR-shift generator.");
  t_nml.addHelp("tcache_config", "Configures the random number cache to hold SINGLE or DOUBLE "
                "precision random number results.  The cache configuration is independent of the "
                "precision with which the random numbers are formed.");
  
  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.  All calls to this function should
  // proceed in consecutive calls, to make use of the updates to start_line and avoid reading any
  // instance of this namelist twice or skipping instances of it in the search for some other
  // namelist.  An alternative is to keep an independent counter to track progress through the
  // input file in search for &rst namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  return t_nml;
}

} // namespace namelist
} // namespace stormm
