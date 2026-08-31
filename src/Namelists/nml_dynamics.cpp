#include "copyright.h"
#include "Constants/symbol_values.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Reporting/error_format.h"
#include "Topology/atomgraph_constants.h"
#include "namelist_common.h"
#include "namelist_element.h"
#include "nml_dynamics.h"

namespace stormm {
namespace namelist {

using constants::getEnumerationName;
using constants::translatePrecisionModel;
using energy::translateVdwSumMethod;
using parse::realToString;
using parse::NumberFormat;
using parse::strcmpCased;
using parse::TextOrigin;
using structure::translateApplyConstraints;
using structure::translateRattleMethod;
using symbols::amber_ancient_bioq;
using topology::amber_default_elec14_screen;
using topology::amber_default_vdw14_screen;
using trajectory::translateThermostatKind;
using trajectory::translateBarostatKind;
  
//-------------------------------------------------------------------------------------------------
DynamicsControls::DynamicsControls(const ExceptionResponse policy_in) :
    policy{policy_in},
    total_step_count{default_dynamics_nstlim},
    diagnostic_frequency{default_dynamics_ntpr},
    trajectory_frequency{default_dynamics_ntwx},
    com_motion_purge_frequency{default_dynamics_nscm},
    time_step{default_dynamics_time_step},
    electrostatic_cutoff{default_electrostatic_cutoff},
    van_der_waals_cutoff{default_van_der_waals_cutoff},
    coulomb{amber_ancient_bioq},
    vdw_style{VdwSumMethod::CUTOFF},
    elec_14_screening{amber_default_elec14_screen},
    vdw_14_screening{amber_default_vdw14_screen},
    coulomb_set_by_user{false},
    elec_14_set_by_user{false},
    vdw_14_set_by_user{false},
    constrain_geometry{std::string(default_geometry_constraint_behavior)},
    use_shake{std::string(default_geometry_constraint_behavior)},
    use_settle{std::string(default_geometry_constraint_behavior)},
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
    barostat_kind{std::string(default_barostat_kind)},
    mcbarostat_frequency_a{default_mcbarostat_frequency},
    mcbarostat_frequency_b{default_mcbarostat_frequency},
    mcbarostat_frequency_c{default_mcbarostat_frequency},
    mcbarostat_rescale_a{default_mcbarostat_rescale},
    mcbarostat_rescale_b{default_mcbarostat_rescale},
    mcbarostat_rescale_c{default_mcbarostat_rescale},
    nt_warp_multiplicity{default_nt_warp_multiplicity},
    initial_temperature_targets{}, final_temperature_targets{}, thermostat_labels{},
    thermostat_label_indices{}, thermostat_masks{}, external_pressures{}, barostat_labels{},
    barostat_label_indices{},
    nml_transcript{"dynamics"}
{
  // Always initialize the thermostat groups with at least one entry comprising all possible
  // systems.  This will be overwritten if namelist input is taken.
  setThermostatGroup();

  // Always initialize the barostat groups with at least on entry, comprising all possible systems
  // and giving them no pressure-based rescaling. (This would be valid for systems without
  // periodic boundary conditions, as well.)
  setBarostatGroup();
  
  // Load in a blank namelist so that certain keywords will be present, as if this were the means
  // by which the data was loaded.
  std::string tfs("&dynamics\n&end\n");
  TextFile tf(tfs, TextOrigin::RAM);
  int start_line = 0;
  bool found;
  nml_transcript = dynamicsInput(tf, &start_line, &found, ExceptionResponse::SILENT);
}

//-------------------------------------------------------------------------------------------------
DynamicsControls::DynamicsControls(const TextFile &tf, int *start_line, bool *found_nml,
                                   const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    DynamicsControls(policy_in)
{
  NamelistEmulator t_nml = dynamicsInput(tf, start_line, found_nml, policy, wrap);
  nml_transcript = t_nml;
  
  // Interpret common keywords
  addRangedInteractionInterpretation(&electrostatic_cutoff, &van_der_waals_cutoff, &vdw_style,
                                     t_nml, policy);
  
  // Interpret namelist-specific keywords
  t_nml.assignVariable(&total_step_count, "nstlim");
  t_nml.assignVariable(&diagnostic_frequency, "ntpr");
  t_nml.assignVariable(&trajectory_frequency, "ntwx");
  t_nml.assignVariable(&com_motion_purge_frequency, "nscm");
  t_nml.assignVariable(&time_step, "dt");
  t_nml.assignVariable(&coulomb, "coulomb");
  t_nml.assignVariable(&elec_14_screening, "scee");
  t_nml.assignVariable(&vdw_14_screening, "scnb");
  t_nml.assignVariable(&rattle_tolerance, "tol");
  t_nml.assignVariable(&rattle_iterations, "rattle_iter");
  setCpuRattleMethod(t_nml.getStringValue("rattle_style"));
  
  // Detect whether SHAKE, RATTLE, and SETTLE should be activated
  if (t_nml.getKeywordStatus("rigid_geom") != InputStatus::MISSING) {
    setGeometricConstraints(t_nml.getStringValue("rigid_geom"));
  }
  if (t_nml.getKeywordStatus("rigid_h") != InputStatus::MISSING) {
    setUseShake(t_nml.getStringValue("rigid_h"));
  }
  if (t_nml.getKeywordStatus("rigid_wat") != InputStatus::MISSING) {
    setUseSettle(t_nml.getStringValue("rigid_wat"));
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
  t_nml.assignVariable(&barostat_kind, "ntp");
  t_nml.assignVariable(&mcbarostat_frequency_a, "mcb_freq_a");
  t_nml.assignVariable(&mcbarostat_frequency_b, "mcb_freq_b");
  t_nml.assignVariable(&mcbarostat_frequency_c, "mcb_freq_c");
  t_nml.assignVariable(&mcbarostat_frequency_a, "mcb_freq");
  t_nml.assignVariable(&mcbarostat_frequency_b, "mcb_freq");
  t_nml.assignVariable(&mcbarostat_frequency_c, "mcb_freq");
  t_nml.assignVariable(&mcbarostat_rescale_a, "mcb_factor_a");
  t_nml.assignVariable(&mcbarostat_rescale_b, "mcb_factor_b");
  t_nml.assignVariable(&mcbarostat_rescale_c, "mcb_factor_c");
  t_nml.assignVariable(&mcbarostat_rescale_a, "mcb_factor");
  t_nml.assignVariable(&mcbarostat_rescale_b, "mcb_factor");
  t_nml.assignVariable(&mcbarostat_rescale_c, "mcb_factor");
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
  const int nbstat = t_nml.getKeywordEntries("pressure");
  external_pressures.resize(nbstat);
  barostat_labels.resize(nbstat);
  barostat_label_indices.resize(nbstat);
  for (int i = 0; i < nbstat; i++) {
    t_nml.assignVariable(&external_pressures[i], "pressure", "pres0", i);
    t_nml.assignVariable(&barostat_labels[i], "pressure", "-label", i);
    t_nml.assignVariable(&barostat_label_indices[i], "pressure", "-n", i);
  }
  
  // If there is an active barostat, check to ensure that there is an active thermostat.
  switch (this->getBarostatKind()) {
  case BarostatKind::NONE:
    break;
  case BarostatKind::MONTE_CARLO:
    switch (this->getThermostatKind()) {
    case ThermostatKind::NONE:
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Isobaric simulations require a thermostat to run properly.  Specify some valid "
              "thermostat with the ntt keyword.", "DynamicsControls");
      case ExceptionResponse::WARN:
        rtWarn("Isobaric simulations require a thermostat to run properly.  A valid "
               "thermostat may be specified with the ntt keyword.  This set of simulations will "
               "use a " + getEnumerationName(ThermostatKind::LANGEVIN) + " thermostat.",
               "DynamicsControls");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
      setThermostatKind(ThermostatKind::LANGEVIN);
      break;
    case ThermostatKind::ANDERSEN:
    case ThermostatKind::LANGEVIN:
    case ThermostatKind::BERENDSEN:
      break;
    }
    break;
  }

  
  // Specifying temperatures without a functioning thermostat is an error.  Specifying pressures
  // without a functioning barostat is an error.  This does not apply if the default settings for
  // all systems' temperatures or pressures are in effect.  Check the input status of the
  // respective keywords.
  InputStatus ntt_status = (ntstat > 1) ? InputStatus::USER_SPECIFIED : InputStatus::DEFAULT;
  int tcon = 0;
  while (ntt_status == InputStatus::DEFAULT && tcon < ntstat) {
    if (t_nml.getKeywordStatus("temperature", "tempi", tcon) == InputStatus::USER_SPECIFIED ||
        t_nml.getKeywordStatus("temperature", "temp0", tcon) == InputStatus::USER_SPECIFIED ||
        t_nml.getKeywordStatus("temperature", "-label", tcon) == InputStatus::USER_SPECIFIED ||
        t_nml.getKeywordStatus("temperature", "-n", tcon) == InputStatus::USER_SPECIFIED ||
        t_nml.getKeywordStatus("temperature", "-mask", tcon) == InputStatus::USER_SPECIFIED) {
      ntt_status = InputStatus::USER_SPECIFIED;
    }
    tcon++;
  }
  switch (ntt_status) {
  case InputStatus::MISSING:
  case InputStatus::DEFAULT:
    break;
  case InputStatus::USER_SPECIFIED:
    switch (this->getThermostatKind()) {
    case ThermostatKind::NONE:
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Thermostating was requested without specifying a means of implementation.",
              "DynamicsControls");
      case ExceptionResponse::WARN:
        rtWarn("Thermostating was requested without specifying a means of implementation.  A " +
               getEnumerationName(ThermostatKind::LANGEVIN) + " thermostat will be used.",
               "DynamicsControls");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
      setThermostatKind(ThermostatKind::LANGEVIN);
      break;
    case ThermostatKind::LANGEVIN:
    case ThermostatKind::ANDERSEN:
    case ThermostatKind::BERENDSEN:
      break;
    }
    break;
  }
  InputStatus ntp_status = (nbstat > 1) ? InputStatus::USER_SPECIFIED : InputStatus::DEFAULT;
  int bcon = 0;
  while (ntp_status == InputStatus::DEFAULT && bcon < nbstat) {
    if (t_nml.getKeywordStatus("pressure", "pres0", bcon) == InputStatus::USER_SPECIFIED ||
        t_nml.getKeywordStatus("pressure", "-label", bcon) == InputStatus::USER_SPECIFIED ||
        t_nml.getKeywordStatus("pressure", "-n", bcon) == InputStatus::USER_SPECIFIED) {
      ntp_status = InputStatus::USER_SPECIFIED;
    }
    bcon++;
  }
  switch (ntp_status) {
  case InputStatus::MISSING:
  case InputStatus::DEFAULT:
    break;
  case InputStatus::USER_SPECIFIED:
    switch (this->getThermostatKind()) {
    case ThermostatKind::NONE:
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Isobaric simulations were requested without specifying a means of implementation.",
              "DynamicsControls");
      case ExceptionResponse::WARN:
        rtWarn("Isobaric simulations were requested without specifying a means of "
               "implementation.  A " + getEnumerationName(BarostatKind::MONTE_CARLO) +
               " barostat will be used.", "DynamicsControls");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
      setBarostatKind(BarostatKind::MONTE_CARLO);
      break;
    case ThermostatKind::LANGEVIN:
    case ThermostatKind::ANDERSEN:
    case ThermostatKind::BERENDSEN:
      break;
    }
  }
  
  // Validate input
  validateStepCount();
  validateDiagnosticPrintFrequency();
  validateTrajectoryPrintFrequency();
  validateCenterOfMassMotionPurgeFrequency();
  validateTimeStep();
  validateRattleTolerance();
  validateRattleIterations();
  validateThermostatKind();
  validateCacheConfiguration();
  validateNTWarpMultiplicity();
  for (size_t i = 0; i < external_pressures.size(); i++) {
    validatePressure(external_pressures[i]);
  }
  validateMCBarostatFrequency();
  validateMCBarostatRescaling();
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
  return com_motion_purge_frequency;
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
double DynamicsControls::getLennardJonesCutoff() const {
  return van_der_waals_cutoff;
}

//-------------------------------------------------------------------------------------------------
double DynamicsControls::getCoulombConstant() const {
  return coulomb;
}

//-------------------------------------------------------------------------------------------------
VdwSumMethod DynamicsControls::getVdwSummation() const {
  return vdw_style;
}

//-------------------------------------------------------------------------------------------------
double DynamicsControls::getElec14Screening() const {
  return elec_14_screening;
}

//-------------------------------------------------------------------------------------------------
double DynamicsControls::getVdw14Screening() const {
  return vdw_14_screening;
}

//-------------------------------------------------------------------------------------------------
bool DynamicsControls::coulombSetByUser() const {
  return coulomb_set_by_user;
}

//-------------------------------------------------------------------------------------------------
bool DynamicsControls::elec14SetByUser() const {
  return elec_14_set_by_user;
}

//-------------------------------------------------------------------------------------------------
bool DynamicsControls::vdw14SetByUser() const {
  return vdw_14_set_by_user;
}

//-------------------------------------------------------------------------------------------------
ApplyConstraints DynamicsControls::constrainGeometry() const {
  return translateApplyConstraints(constrain_geometry);
}

//-------------------------------------------------------------------------------------------------
ApplyConstraints DynamicsControls::useShake() const {
  return translateApplyConstraints(use_shake);
}

//-------------------------------------------------------------------------------------------------
ApplyConstraints DynamicsControls::useSettle() const {
  return translateApplyConstraints(use_settle);
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
BarostatKind DynamicsControls::getBarostatKind() const {
  return translateBarostatKind(barostat_kind);
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getMCBarostatFrequency() const {
  if (mcbarostat_frequency_b != mcbarostat_frequency_a ||
      mcbarostat_frequency_c != mcbarostat_frequency_a) {
    rtErr("Anisotropic rescaling is in effect, with frequencies { " +
          std::to_string(mcbarostat_frequency_a) + ", " + std::to_string(mcbarostat_frequency_b) +
          ", " + std::to_string(mcbarostat_frequency_c) + " }.  A meaningful move rate can only "
          "be given for a specific dimension.", "DynamicsControls", "getMCBarostatFrequency");
  }
  return mcbarostat_frequency_a;
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getMCBarostatFrequency(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return mcbarostat_frequency_a;
  case CartesianDimension::Y:
    return mcbarostat_frequency_b;
  case CartesianDimension::Z:
    return mcbarostat_frequency_c;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int DynamicsControls::getMCBarostatFrequency(const UnitCellAxis dim) const {
  switch (dim) {
  case UnitCellAxis::A:
    return mcbarostat_frequency_a;
  case UnitCellAxis::B:
    return mcbarostat_frequency_b;
  case UnitCellAxis::C:
    return mcbarostat_frequency_c;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double DynamicsControls::getMCBarostatRescaling() const {
  if (mcbarostat_rescale_b != mcbarostat_rescale_a ||
      mcbarostat_rescale_c != mcbarostat_rescale_a) {
    rtErr("Anisotropic rescaling is in effect, with rescaling factors { " +
          std::to_string(mcbarostat_rescale_a) + ", " + std::to_string(mcbarostat_rescale_b) +
          ", " + std::to_string(mcbarostat_rescale_c) + " }.  A meaningful move rate can only "
          "be given for a specific dimension.", "DynamicsControls", "getMCBarostatRescale");
  }
  return mcbarostat_rescale_a;
}

//-------------------------------------------------------------------------------------------------
double DynamicsControls::getMCBarostatRescaling(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return mcbarostat_rescale_a;
  case CartesianDimension::Y:
    return mcbarostat_rescale_b;
  case CartesianDimension::Z:
    return mcbarostat_rescale_c;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double DynamicsControls::getMCBarostatRescaling(const UnitCellAxis dim) const {
  switch (dim) {
  case UnitCellAxis::A:
    return mcbarostat_rescale_a;
  case UnitCellAxis::B:
    return mcbarostat_rescale_b;
  case UnitCellAxis::C:
    return mcbarostat_rescale_c;
  }
  __builtin_unreachable();
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
const std::vector<double>& DynamicsControls::getExternalPressures() const {
  return external_pressures;
}
  
//-------------------------------------------------------------------------------------------------
const std::vector<std::string>& DynamicsControls::getBarostatLabels() const {
  return barostat_labels;
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& DynamicsControls::getBarostatLabelIndices() const {
  return barostat_label_indices;
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
DynamicsControls::setCenterOfMassMotionPurgeFrequency(const int com_motion_purge_frequency_in) {
  com_motion_purge_frequency = com_motion_purge_frequency_in;
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setTimeStep(const double time_step_in) {
  time_step = time_step_in;
  validateTimeStep();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setElectrostaticCutoff(const double cutoff_in) {
  electrostatic_cutoff = cutoff_in;
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setVanDerWaalsCutoff(const double cutoff_in) {
  van_der_waals_cutoff = cutoff_in;
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setCutoff(const double cutoff_in) {
  electrostatic_cutoff = cutoff_in;
  van_der_waals_cutoff = cutoff_in;
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setVdwSummation(const std::string &vdw_method_in) {
  vdw_style = translateVdwSumMethod(vdw_method_in);
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setVdwSummation(const VdwSumMethod vdw_method_in) {
  vdw_style = vdw_method_in;
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setCoulombConstant(const double coulomb_in) {
  coulomb = coulomb_in;
  coulomb_set_by_user = true;
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setElec14Screening(const double screening_in) {
  elec_14_screening = screening_in;
  elec_14_set_by_user = true;
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setVdw14Screening(const double screening_in) {
  vdw_14_screening = screening_in;
  vdw_14_set_by_user = true;
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
void DynamicsControls::setUseShake(const std::string &use_shake_in) {
  use_shake = use_shake_in;
  try {
    const ApplyConstraints trial = translateApplyConstraints(use_shake_in);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Invalid constraint activation \"" + use_shake_in + "\" provided.",
            "DynamicsControls", "setUseShake");
    case ExceptionResponse::WARN:
      rtWarn("Invalid constraint activation \"" + use_shake_in + "\" provided.  The default of " +
             std::string(default_geometry_constraint_behavior) + " will be reinstated.",
             "DynamicsControls", "setUseShake");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    use_shake = std::string(default_geometry_constraint_behavior);
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setUseShake(const ApplyConstraints use_shake_in) {
  use_shake = getEnumerationName(use_shake_in);
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setUseSettle(const std::string &use_settle_in) {
  use_settle = use_settle_in;
  try {
    const ApplyConstraints trial = translateApplyConstraints(use_settle_in);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Invalid constraint activation \"" + use_settle_in + "\" provided.",
            "DynamicsControls", "setUseSettle");
    case ExceptionResponse::WARN:
      rtWarn("Invalid constraint activation \"" + use_settle_in + "\" provided.  The default of " +
             std::string(default_geometry_constraint_behavior) + " will be reinstated.",
             "DynamicsControls", "setUseSettle");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    use_settle = std::string(default_geometry_constraint_behavior);
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setUseSettle(const ApplyConstraints use_settle_in) {
  use_settle = getEnumerationName(use_settle_in);
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
void DynamicsControls::setBarostatKind(const std::string &barostat_kind_in) {
  barostat_kind = barostat_kind_in;
  validateBarostatKind();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setBarostatKind(const BarostatKind barostat_kind_in) {
  barostat_kind = getEnumerationName(barostat_kind_in);
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setMCBarostatFrequency(const int frequency_in) {
  mcbarostat_frequency_a = frequency_in;
  mcbarostat_frequency_b = frequency_in;
  mcbarostat_frequency_c = frequency_in;
  validateMCBarostatFrequency();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setMCBarostatFrequency(const int frequency_in,
                                              const CartesianDimension dim) {
  switch (dim) {
  case CartesianDimension::X:
    mcbarostat_frequency_a = frequency_in;
    break;
  case CartesianDimension::Y:
    mcbarostat_frequency_b = frequency_in;
    break;
  case CartesianDimension::Z:
    mcbarostat_frequency_c = frequency_in;
    break;
  }
  validateMCBarostatFrequency();
}
  
//-------------------------------------------------------------------------------------------------
void DynamicsControls::setMCBarostatFrequency(const int frequency_in, const UnitCellAxis dim) {
  switch (dim) {
  case UnitCellAxis::A:
    mcbarostat_frequency_a = frequency_in;
    break;
  case UnitCellAxis::B:
    mcbarostat_frequency_b = frequency_in;
    break;
  case UnitCellAxis::C:
    mcbarostat_frequency_c = frequency_in;
    break;
  }
  validateMCBarostatFrequency();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setMCBarostatRescaling(const double rescale_in) {
  mcbarostat_rescale_a = rescale_in;
  mcbarostat_rescale_b = rescale_in;
  mcbarostat_rescale_c = rescale_in;
  validateMCBarostatRescaling();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setMCBarostatRescaling(const double rescale_in,
                                              const CartesianDimension dim) {
  switch (dim) {
  case CartesianDimension::X:
    mcbarostat_rescale_a = rescale_in;
    break;
  case CartesianDimension::Y:
    mcbarostat_rescale_b = rescale_in;
    break;
  case CartesianDimension::Z:
    mcbarostat_rescale_c = rescale_in;
    break;
  }
  validateMCBarostatRescaling();
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::setMCBarostatRescaling(const double rescale_in, const UnitCellAxis dim) {
  switch (dim) {
  case UnitCellAxis::A:
    mcbarostat_rescale_a = rescale_in;
    break;
  case UnitCellAxis::B:
    mcbarostat_rescale_b = rescale_in;
    break;
  case UnitCellAxis::C:
    mcbarostat_rescale_c = rescale_in;
    break;
  }
  validateMCBarostatRescaling();
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
void DynamicsControls::setBarostatGroup(const double pressure_target, const std::string &label,
                                        const int label_index) {
  validatePressure(pressure_target);
  external_pressures.push_back(pressure_target);
  barostat_labels.push_back(label);
  barostat_label_indices.push_back(label_index);
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
  if (com_motion_purge_frequency < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A negative value of the center of mass motion purge frequency (" +
            std::to_string(com_motion_purge_frequency) + ") is invalid.", "DynamicsControls",
            "validateCenterOfMassMotionPurgeFrequency");
    case ExceptionResponse::WARN:
      rtWarn("A negative value of the center of mass motion purge frequency (" +
             std::to_string(com_motion_purge_frequency) + ") is invalid and will be adjusted to "
             "zero (no purging).", "DynamicsControls", "validateCenterOfMassMotionPurgeFrequency");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    com_motion_purge_frequency = 0;
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
void DynamicsControls::validateBarostatKind() {
  try {
    const BarostatKind trial = translateBarostatKind(barostat_kind);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("An invalid barostat type " + barostat_kind + " was detected in the input.",
            "DynamicsControls", "validateBarostatKind");
    case ExceptionResponse::WARN:
      rtWarn("An invalid barostat type " + barostat_kind + " was detected in the input.  It "
             "will be replaced with the default of " + std::string(default_barostat_kind) + ".",
             "DynamicsControls", "validateBarostatKind");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    barostat_kind = std::string(default_barostat_kind);
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validatePressure(const double p) const {
  if (p < 0.0 || p > 10000.0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A pressure of " + realToString(p, 9, 4, NumberFormat::STANDARD_REAL) + " is "
            "unrealistic for molecular dynamics.", "DynamicsControls", "validatePressure");
    case ExceptionResponse::WARN:
      rtWarn("A pressure of " + realToString(p, 9, 4, NumberFormat::STANDARD_REAL) + " is "
             "unrealistic for molecular dynamics.", "DynamicsControls", "validatePressure");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateMCBarostatFrequency() {
  if (mcbarostat_frequency_a < 0 || mcbarostat_frequency_b < 0 || mcbarostat_frequency_c < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("It is invalid to attempt Monte-Carlo moves at every { " +
            std::to_string(mcbarostat_frequency_a) + ", " +
            std::to_string(mcbarostat_frequency_b) + ", " +
            std::to_string(mcbarostat_frequency_c) + " } steps along each axis.",
            "DynamicsControls", "validateMCBarostatFrequency");
    case ExceptionResponse::WARN:
      rtErr("It is invalid to attempt Monte-Carlo moves at every { " +
            std::to_string(mcbarostat_frequency_a) + ", " +
            std::to_string(mcbarostat_frequency_b) + ", " +
            std::to_string(mcbarostat_frequency_c) + " } steps along each axis.  The default move "
            "frequency of " + std::to_string(default_mcbarostat_frequency) + " will be applied "
            "throughout.", "DynamicsControls", "validateMCBarostatFrequency");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    mcbarostat_frequency_a = default_mcbarostat_frequency;
    mcbarostat_frequency_b = default_mcbarostat_frequency;
    mcbarostat_frequency_c = default_mcbarostat_frequency;
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsControls::validateMCBarostatRescaling() {
  if (mcbarostat_rescale_a < 0.0 || mcbarostat_rescale_b < 0.0 || mcbarostat_rescale_c < 0.0 ||
      mcbarostat_rescale_a > 0.1 || mcbarostat_rescale_b > 0.1 || mcbarostat_rescale_c > 0.1) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Monte-Carlo barostat rescaling factors of { " +
            realToString(mcbarostat_rescale_a, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
            realToString(mcbarostat_rescale_b, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
            realToString(mcbarostat_rescale_c, 7, 4, NumberFormat::STANDARD_REAL) + " } are "
            "invalid.", "DynamicsControls", "validateMCBarostatRescaling");
    case ExceptionResponse::WARN:
      rtWarn("Monte-Carlo barostat rescaling factors of { " +
             realToString(mcbarostat_rescale_a, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
             realToString(mcbarostat_rescale_b, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
             realToString(mcbarostat_rescale_c, 7, 4, NumberFormat::STANDARD_REAL) + " } are "
             "invalid.  The default value of " +
             realToString(default_mcbarostat_rescale, 7, 4, NumberFormat::STANDARD_REAL) +
             " will be applied throughout.", "DynamicsControls", "validateMCBarostatRescaling");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    mcbarostat_rescale_a = default_mcbarostat_rescale;
    mcbarostat_rescale_b = default_mcbarostat_rescale;
    mcbarostat_rescale_c = default_mcbarostat_rescale;
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

  // Common keyword handling
  addRangedInteractionControls(&t_nml);
  
  // Trajectory length and outut keywords
  t_nml.addKeyword("nstlim", NamelistType::INTEGER, std::to_string(default_dynamics_nstlim));
  t_nml.addKeyword("ntpr", NamelistType::INTEGER, std::to_string(default_dynamics_ntpr));
  t_nml.addKeyword("ntwx", NamelistType::INTEGER, std::to_string(default_dynamics_ntwx));
  t_nml.addKeyword("nscm", NamelistType::INTEGER, std::to_string(default_dynamics_nscm));
  t_nml.addKeyword("dt", NamelistType::REAL, std::to_string(default_dynamics_time_step));
  t_nml.addKeyword("coulomb", NamelistType::REAL, std::to_string(amber_ancient_bioq));
  t_nml.addKeyword("scee", NamelistType::REAL, std::to_string(amber_default_elec14_screen));
  t_nml.addKeyword("scnb", NamelistType::REAL, std::to_string(amber_default_vdw14_screen));
  
  // Constraint keywords
  t_nml.addKeyword("rigid_geom", NamelistType::STRING);
  t_nml.addKeyword("rigid_h", NamelistType::STRING);
  t_nml.addKeyword("rigid_wat", NamelistType::STRING);
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
    "A particular system index within the label group to which this temperature profile applies.  "
    "The default value of " + std::to_string(-1) + " implies that all members of the label group "
    "will be affected.  Each declaration of the temperature keyword struct may single out one "
    "such system for special temperature regulation.",
    "Atom mask for atoms in the system affected by this temperature profile"
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

  // Barostating keywords
  t_nml.addKeyword("ntp", NamelistType::STRING, std::string(default_barostat_kind));
  t_nml.addKeyword("mcb_freq", NamelistType::INTEGER,
                   std::to_string(default_mcbarostat_frequency));
  t_nml.addKeyword("mcb_freq_a", NamelistType::INTEGER,
                   std::to_string(default_mcbarostat_frequency));
  t_nml.addKeyword("mcb_freq_b", NamelistType::INTEGER,
                   std::to_string(default_mcbarostat_frequency));
  t_nml.addKeyword("mcb_freq_c", NamelistType::INTEGER,
                   std::to_string(default_mcbarostat_frequency));
  t_nml.addKeyword("mcb_factor", NamelistType::REAL,
                   realToString(default_mcbarostat_rescale, 10, 7, NumberFormat::STANDARD_REAL));
  t_nml.addKeyword("mcb_factor_a", NamelistType::REAL,
                   realToString(default_mcbarostat_rescale, 10, 7, NumberFormat::STANDARD_REAL));
  t_nml.addKeyword("mcb_factor_b", NamelistType::REAL,
                   realToString(default_mcbarostat_rescale, 10, 7, NumberFormat::STANDARD_REAL));
  t_nml.addKeyword("mcb_factor_c", NamelistType::REAL,
                   realToString(default_mcbarostat_rescale, 10, 7, NumberFormat::STANDARD_REAL));
  const std::string pressure_help("Specify the external pressure at which to maintain a system or "
                                  "group of systems.");
  const std::vector<std::string> pressure_keys_help = {
    "The target external pressure to be applied by this barostat",
    "A label from the systems cache specifying a group of systems to which this barostat applies",
    "A particular system index within the label group to which this barostat applies"
  };
  t_nml.addKeyword("pressure", { "pres0", "-label", "-n" },
                   { NamelistType::REAL, NamelistType::STRING, NamelistType::INTEGER },
                   { realToString(default_external_pressure, 9, 4, std_real), std::string("all"),
                     std::to_string(-1) }, DefaultIsObligatory::YES, InputRepeats::YES,
                   pressure_help, pressure_keys_help,
                   { KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                     KeyRequirement::REQUIRED });
  
  // Help messages for each trajectory keyword
  t_nml.addHelp("nstlim", "Number of dynamics steps to carry out");
  t_nml.addHelp("ntpr", "The frequency with which to calculate and print energy diagnostics");
  t_nml.addHelp("ntwx", "The frequency with which to print trajectory snapshots");
  t_nml.addHelp("nscm", "The frequency with which to purge motion of the center of mass");
  t_nml.addHelp("dt", "The time step, in units of femtoseconds");

  // Help messages for interaction cutoffs
  t_nml.addHelp("coulomb", "Define Coulomb's constant for the simulation, in units of "
                "kcal/mol-e^2, where e is the charge of a proton (atomic unit of charge).");
  t_nml.addHelp("scee", "The unitless screening factor on electrostatic 1:4 interactions.  To "
                "attenuate such interactions between near neighbors to 83.333% of their nominal "
                "value, specify 6/5, i.e. 1.2.");
  t_nml.addHelp("scnb", "The unitless screening factor on van-der Waals 1:4 interactions.  To "
                "attenuate such interactions between near neighbors to 50% of their nominal "
                "value, specify 1/2, i.e. 0.5.");
  
  // Help messages for geometry constraints keywords
  t_nml.addHelp("rigid_geom", "Indicate whether to enforce all rigid geometries, both bond length "
                "constraints and rigid water molecules");
  t_nml.addHelp("rigid_h", "Indicate whether to enforce rigid bond lengths to hydrogen atoms");
  t_nml.addHelp("rigid_wat", "Indicate whether to enforce the geometry of trigonal water "
                "molecules bonded in a ring.  Rigid water models such as TIP3P, SPC/E, TIP4P-Ew, "
                "and others make use of these analytic constraints.");
  t_nml.addHelp("tol", "Tolerance by which to constrain rigid bonds involving hydrogen atoms.  "
                "The units of this tolerance are squared Angstroms, implying that the rigid "
                "geometry will be correct to within the square root of this tolerance after all "
                "iterations are complete.");
  t_nml.addHelp("rattle_iter", "Maximum number of iterations to use when attempting to converge "
                "constrained bond lengths");
  t_nml.addHelp("rattle_style", "The manner in which to converge 'hub and spoke' constrained "
                "bonds connected to a single central atom");

  // Help messages for non-struct thermostat keywords
  t_nml.addHelp("ntt", "Indicate the type of thermostat, whether by a human-readable string or "
                "the AMBER numeric cognate.  Options are case-insensitive and include \"none\" "
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
  t_nml.addHelp("ntp", "Indicate the type of barostat, either by a human-readable string or the "
                "AMBER numeric cognate.   Options are case-insensitive and include \"none\" "
                "(0) or \"monte-carlo\" (2).");
  t_nml.addHelp("mcb_freq", "The number of time steps between Monte-Carlo barostat volume "
                "adjustments.  This will apply to all axes of each system's unit cells, if "
                "specified.");
  t_nml.addHelp("mcb_freq_a", "The number of time steps between Monte-Carlo barostat volume "
                "adjustments along each system's unit cell A axis");
  t_nml.addHelp("mcb_freq_b", "The number of time steps between Monte-Carlo barostat volume "
                "adjustments along each system's unit cell B axis");
  t_nml.addHelp("mcb_freq_c", "The number of time steps between Monte-Carlo barostat volume "
                "adjustments along each system's unit cell C axis");
  t_nml.addHelp("mcb_factor", "The proportion by which the unit cell can be rescaled along all "
                "axes in Monte-Carlo barostating.  Specifying this parameter will supercede "
                "rescaling inputs for individual axes.");
  t_nml.addHelp("mcb_factor_a", "The proportion by which the unit cell can be rescaled along its "
                "A axis in Monte-Carlo barostating");
  t_nml.addHelp("mcb_factor_b", "The proportion by which the unit cell can be rescaled along its "
                "B axis in Monte-Carlo barostating");
  t_nml.addHelp("mcb_factor_c", "The proportion by which the unit cell can be rescaled along its "
                "C axis in Monte-Carlo barostating");
  
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
