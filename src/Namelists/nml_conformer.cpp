#include "copyright.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Reporting/error_format.h"
#include "Structure/local_arrangement.h"
#include "Structure/structure_enumerators.h"
#include "input.h"
#include "nml_conformer.h"

namespace stormm {
namespace namelist {

using constants::CaseSensitivity;
using errors::rtErr;
using errors::rtWarn;
using parse::CaseSensitivity;
using parse::NumberFormat;
using parse::realToString;
using parse::strcmpCased;
using parse::WrapTextSearch;
using structure::imageValue;
using structure::ImagingMethod;
using structure::translateSamplingIntensity;
using synthesis::translateSystemGrouping;
using synthesis::translateVariableTorsionAdjustment;
  
//-------------------------------------------------------------------------------------------------
ConformerControls::ConformerControls(const ExceptionResponse policy_in) :
    policy{policy_in}, core_atom_mask{std::string("")}, core_data_item{std::string("")},
    core_rk2{0.0}, core_rk3{0.0}, core_r2{0.0}, core_r3{0.0}, anchor_conformation{std::string("")},
    sample_chirality{false}, sample_cis_trans{false}, prevent_hbonds{false},
    running_states{default_conf_running_states},
    final_states{default_conf_final_states},
    rotation_sample_count{default_conf_rotation_samples},
    rotatable_bond_limit{default_conf_max_rotatable_bonds},
    cis_trans_sample_count{default_conf_cis_trans_samples},
    max_seeding_attempts{default_conf_max_seeding_attempts},
    clash_pair_tolerance{default_conf_clash_pairs},
    sampling_trial_limit{default_conf_sample_trials},
    maximum_trial_limit{default_conf_max_system_trials},
    rmsd_tolerance{default_conf_rmsd_tolerance},
    group_method{default_conf_output_grouping},
    sample_effort{default_conf_sampling_effort},
    rotation_snap_threshold{stod(std::string(default_conf_rotation_snap)) * pi / 180.0},
    cis_trans_snap_threshold{stod(std::string(default_conf_cis_trans_snap)) * pi / 180.0},
    adjustment_method{default_conf_adjustment_method},
    rotation_sample_values{}, cis_trans_sample_values{},
    nml_transcript{"conformer"}
{}

//-------------------------------------------------------------------------------------------------
ConformerControls::ConformerControls(const TextFile &tf, int *start_line, bool *found_nml,
                                     const ExceptionResponse policy_in,
                                     const WrapTextSearch wrap) :
    ConformerControls(policy_in)
{
  const NamelistEmulator t_nml = conformerInput(tf, start_line, found_nml, policy, wrap);
  nml_transcript = t_nml;
  if (t_nml.getKeywordStatus("core_mask") != InputStatus::MISSING) {
    core_atom_mask = t_nml.getStringValue("core_mask", "atoms");
    core_data_item = t_nml.getStringValue("core_mask", "data_item");

    // Get the rk2 restraint value
    if (t_nml.getKeywordStatus("core_mask", "rk2") != InputStatus::MISSING) {
      core_rk2 = t_nml.getRealValue("core_mask", "rk2");
    }
    else if (t_nml.getKeywordStatus("core_mask", "repulsion") != InputStatus::MISSING) {
      core_rk2 = t_nml.getRealValue("core_mask", "repulsion");
    }
    else if (t_nml.getKeywordStatus("core_mask", "stiffness") != InputStatus::MISSING) {
      core_rk2 = t_nml.getRealValue("core_mask", "stiffness");
    }

    // Get the rk3 restraint value
    if (t_nml.getKeywordStatus("core_mask", "rk3") != InputStatus::MISSING) {
      core_rk3 = t_nml.getRealValue("core_mask", "rk3");
    }
    else if (t_nml.getKeywordStatus("core_mask", "attraction") != InputStatus::MISSING) {
      core_rk3 = t_nml.getRealValue("core_mask", "attraction");
    }
    else if (t_nml.getKeywordStatus("core_mask", "stiffness") != InputStatus::MISSING) {
      core_rk3 = t_nml.getRealValue("core_mask", "stiffness");
    }

    // Get the r2 restraint value
    if (t_nml.getKeywordStatus("core_mask", "r2") != InputStatus::MISSING) {
      core_r2 = t_nml.getRealValue("core_mask", "r2");
    }
    else if (t_nml.getKeywordStatus("core_mask", "demand") != InputStatus::MISSING) {
      core_r2 = t_nml.getRealValue("core_mask", "demand");
    }

    // Get the r3 restraint value
    if (t_nml.getKeywordStatus("core_mask", "r3") != InputStatus::MISSING) {
      core_r3 = t_nml.getRealValue("core_mask", "r3");
    }
    else if (t_nml.getKeywordStatus("core_mask", "grace") != InputStatus::MISSING) {
      core_r3 = t_nml.getRealValue("core_mask", "grace");
    }
  }
  if (t_nml.getKeywordStatus("anchor_conf") != InputStatus::MISSING) {
    anchor_conformation = t_nml.getStringValue("anchor_conf");
  }
  sample_chirality = strcmpCased(t_nml.getStringValue("sample_chirality"), "true",
                                 CaseSensitivity::NO);
  sample_cis_trans = strcmpCased(t_nml.getStringValue("sample_cis_trans"), "true",
                                 CaseSensitivity::NO);
  prevent_hbonds = strcmpCased(t_nml.getStringValue("prevent_hbonds"), "true",
                               CaseSensitivity::NO);
  running_states = t_nml.getIntValue("running_states");
  final_states = t_nml.getIntValue("final_states");
  rotatable_bond_limit = t_nml.getIntValue("max_rotatable_bonds");
  max_seeding_attempts = t_nml.getIntValue("max_seeding_attempts");
  clash_pair_tolerance = t_nml.getIntValue("clash_pair_tol");
  maximum_trial_limit = t_nml.getIntValue("trial_limit");
  sampling_trial_limit = t_nml.getIntValue("local_trial_limit");
  rmsd_tolerance = t_nml.getRealValue("rmsd_tol");
  group_method = t_nml.getStringValue("grouping");
  sample_effort = t_nml.getStringValue("effort");
  adjustment_method = t_nml.getStringValue("rotamer_adjustment");

  // Pull in the rotatable bond and cis-trans isomeric bond samples
  processSamplingValues(&rotation_sample_count, "rotation_sample_count", &rotation_sample_values,
                        "rotation_sample", 120.0, 10.0, t_nml);
  processSamplingValues(&cis_trans_sample_count, "cis_trans_sample_count",
                        &cis_trans_sample_values, "cis_trans_sample", 180.0, 5.0, t_nml);
  
  // Validate input
  validateSampleChirality(t_nml.getStringValue("sample_chirality"));
  validateSampleCisTrans(t_nml.getStringValue("sample_cis_trans"));
  validatePreventHBonds(t_nml.getStringValue("prevent_hbonds"));
  validateStateCounts();
  validateGroupingMethod();
  validateSeedingAttempts();
  validateClashCounts();
  validateSamplingIntensity();
  validateTorsionAdjustmentProtocol();
  validateTorsionSnapping(&rotation_snap_threshold, "rotatable bond");
  validateTorsionSnapping(&cis_trans_snap_threshold, "cis-trans isomeric bond");
}

//-------------------------------------------------------------------------------------------------
const std::string& ConformerControls::getCoreAtomMask() const {
  return core_atom_mask;
}

//-------------------------------------------------------------------------------------------------
const std::string& ConformerControls::getCoreDataItemName() const {
  return core_data_item;
}

//-------------------------------------------------------------------------------------------------
double ConformerControls::getCoreRK2Value() const {
  return core_rk2;
}

//-------------------------------------------------------------------------------------------------
double ConformerControls::getCoreRK3Value() const {
  return core_rk3;
}

//-------------------------------------------------------------------------------------------------
double ConformerControls::getCoreR2Value() const {
  return core_r2;
}

//-------------------------------------------------------------------------------------------------
double ConformerControls::getCoreR3Value() const {
  return core_r3;
}

//-------------------------------------------------------------------------------------------------
bool ConformerControls::sampleChirality() const {
  return sample_chirality;
}

//-------------------------------------------------------------------------------------------------
bool ConformerControls::sampleCisTrans() const {
  return sample_cis_trans;
}

//-------------------------------------------------------------------------------------------------
bool ConformerControls::preventHydrogenBonding() const {
  return prevent_hbonds;
}

//-------------------------------------------------------------------------------------------------
int ConformerControls::getRunningStateCount() const {
  return running_states;
}

//-------------------------------------------------------------------------------------------------
int ConformerControls::getFinalStateCount() const {
  return final_states;
}

//-------------------------------------------------------------------------------------------------
int ConformerControls::getRotationSampleCount() const {
  return rotation_sample_count;
}

//-------------------------------------------------------------------------------------------------
int ConformerControls::getCisTransSampleCount() const {
  return cis_trans_sample_count;
}

//-------------------------------------------------------------------------------------------------
const std::vector<double>& ConformerControls::getRotationSampleValues() const {
  return rotation_sample_values;
}

//-------------------------------------------------------------------------------------------------
const std::vector<double>& ConformerControls::getCisTransSampleValues() const {
  return cis_trans_sample_values;
}

//-------------------------------------------------------------------------------------------------
int ConformerControls::getRotatableBondLimit() const {
  return rotatable_bond_limit;
}

//-------------------------------------------------------------------------------------------------
int ConformerControls::getMaxSeedingAttempts() const {
  return max_seeding_attempts;
}

//-------------------------------------------------------------------------------------------------
int ConformerControls::getClashPairTolerance() const {
  return clash_pair_tolerance;
}

//-------------------------------------------------------------------------------------------------
int ConformerControls::getSamplingTrialLimit() const {
  return sampling_trial_limit;
}

//-------------------------------------------------------------------------------------------------
int ConformerControls::getMaximumTrialLimit() const {
  return maximum_trial_limit;
}

//-------------------------------------------------------------------------------------------------
double ConformerControls::getRMSDTolerance() const {
  return rmsd_tolerance;
}

//-------------------------------------------------------------------------------------------------
SystemGrouping ConformerControls::getGroupingMethod() const {
  return translateSystemGrouping(group_method);
}

//-------------------------------------------------------------------------------------------------
SamplingIntensity ConformerControls::getSamplingIntensity() const {
  return translateSamplingIntensity(sample_effort);
}

//-------------------------------------------------------------------------------------------------
double ConformerControls::getRotatableBondSnapThreshold() const {
  return rotation_snap_threshold;
}

//-------------------------------------------------------------------------------------------------
double ConformerControls::getCisTransBondSnapThreshold() const {
  return cis_trans_snap_threshold;
}

//-------------------------------------------------------------------------------------------------
VariableTorsionAdjustment ConformerControls::getTorsionAdjustmentProtocol() const {
  return translateVariableTorsionAdjustment(adjustment_method);
}

//-------------------------------------------------------------------------------------------------
const NamelistEmulator& ConformerControls::getTranscript() const {
  return nml_transcript;
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validateCoreRestraint() const {
  if (core_r2 < 0.0 || core_r2 > core_r3) {
    rtErr("Displacement limits of " + realToString(core_r2, 7, 4, NumberFormat::STANDARD_REAL) +
          " to " + realToString(core_r3, 7, 4, NumberFormat::STANDARD_REAL) + " are invalid.",
          "ConformerControls", "validateCoreRestraint");
  }
  if (core_rk2 < 0.0) {
    rtErr("A repulsive restraint penalty of " +
          realToString(core_rk2, 7, 4, NumberFormat::STANDARD_REAL) + " kcal/mol-A^2 is "
          "unacceptable.", "ConformerControls", "validateCoreRestraint");
  }
  if (core_rk3 < 0.0) {
    rtErr("A restraint penalty of " + realToString(core_rk3, 7, 4, NumberFormat::STANDARD_REAL) +
          " kcal/mol-A^2 is unacceptable.", "ConformerControls", "validateCoreRestraint");
  }
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validateSampleChirality(const std::string &directive) const {
  if (strcmpCased(directive, "true", CaseSensitivity::NO) == false &&
      strcmpCased(directive, "false", CaseSensitivity::NO) == false) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The sample_chirality keyword accepts input of 'true' or 'false'.  Input " +
            directive + " is invalid.", "ConformerControls", "validateSampleChirality");
    case ExceptionResponse::WARN:
      rtWarn("The sample_chirality keyword accepts input of 'true' or 'false'.  Input " +
             directive + " is invalid.  Chirality will not be sampled.", "ConformerControls",
             "validateSampleChirality");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validateSampleCisTrans(const std::string &directive) const {
  if (strcmpCased(directive, "true", CaseSensitivity::NO) == false &&
      strcmpCased(directive, "false", CaseSensitivity::NO) == false) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The sample_cis_trans keyword accepts input of 'true' or 'false'.  Input " +
            directive + " is invalid.", "ConformerControls", "validateSampleCisTrans");
    case ExceptionResponse::WARN:
      rtWarn("The sample_cis_trans keyword accepts input of 'true' or 'false'.  Input " +
             directive + " is invalid.  Cis- and trans- states of molecules will not be sampled.",
             "ConformerControls", "validateSampleCisTrans");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validatePreventHBonds(const std::string &directive) const {
  if (strcmpCased(directive, "true", CaseSensitivity::NO) == false &&
      strcmpCased(directive, "false", CaseSensitivity::NO) == false) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The prevent_hbonds keyword accepts input of 'true' or 'false' to apply restraints "
            "that will inhibit hydrogen bond formation within the molecules.  Input " +
            directive + " is invalid.", "ConformerControls", "validateSampleCisTrans");
    case ExceptionResponse::WARN:
      rtWarn("The prevent_hbonds keyword accepts input of 'true' or 'false' to apply restraints "
            "that will inhibit hydrogen bond formation within the molecules.  Input " +
             directive + " is invalid.  Intramolecular hydrogen bonds will be free to form.",
             "ConformerControls", "validatePreventHBonds");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validateStateCounts() {
  if (running_states > active_states_limit) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The number of running states (" + std::to_string(running_states) + ") cannot exceed "
            "the maximum allowed number of " + std::to_string(active_states_limit) + ".",
            "ConformerControls", "validateStateCounts");
    case ExceptionResponse::WARN:
      rtWarn("The number of running states (" + std::to_string(running_states) + ") exceeds "
             "the maximum allowed number of " + std::to_string(active_states_limit) + " and will "
             "be reduced to limit memory usage.", "ConformerControls", "validateStateCounts");
      running_states = active_states_limit;
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (maximum_trial_limit > active_states_limit) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The number of energy minimizations for a single system (" +
            std::to_string(maximum_trial_limit) + ") cannot exceed the maximum allowed number "
            "of " + std::to_string(active_states_limit) + ".", "ConformerControls",
            "validateStateCounts");
    case ExceptionResponse::WARN:
      rtWarn("The number of energy minimizations for a single system (" +
             std::to_string(maximum_trial_limit) + ") cannot exceed the maximum allowed number "
             "of " + std::to_string(active_states_limit) + " and will be reduced.",
             "ConformerControls", "validateStateCounts");
      maximum_trial_limit = active_states_limit;
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (final_states > maximum_trial_limit) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The number of final states (" + std::to_string(final_states) + ") cannot exceed the "
            "number of energy minimizations evaluated per system (" +
            std::to_string(maximum_trial_limit) + ").", "ConformerControls",
            "validateStateCounts");
    case ExceptionResponse::WARN:
      rtWarn("The number of final states (" + std::to_string(final_states) + ") cannot exceed the "
             "number of energy minimizations evaluated per system (" +
             std::to_string(maximum_trial_limit) + ") and will be reduced.", "ConformerControls",
             "validateStateCounts");
      final_states = maximum_trial_limit;
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (rmsd_tolerance < 0.0) {
    switch (policy) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      rtWarn("A negative positional root mean squared deviation criterion to distinguish unique "
             "conformers is non-sensical.  A value of zero will be applied.", "ConformerControls",
             "validateStateCounts");
      rmsd_tolerance = 0.0;
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validateSeedingAttempts() {
  if (max_seeding_attempts < 0 || max_seeding_attempts >= seeding_attempts_limit) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The number of conformer seeding attempts (" + std::to_string(max_seeding_attempts) +
            ") is out of range (0 - 100).", "ConformerControls", "validateSeedingAttempts");
    case ExceptionResponse::WARN:
      rtWarn("The number of conformer seeding attempts (" + std::to_string(max_seeding_attempts) +
             ") is out of range (0 - 100).  The value will be clamped within the appropriate "
             "range.", "ConformerControls", "validateSeedingAttempts");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    if (max_seeding_attempts < 0) {
      max_seeding_attempts = 0;
    }
    else {
      max_seeding_attempts = seeding_attempts_limit;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validateClashCounts() {
  if (clash_pair_tolerance < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The allowed number of clashing atom pairs in any given structure is invalid (" +
            std::to_string(clash_pair_tolerance) + ").", "ConformerControls",
            "validateClashCounts");
    case ExceptionResponse::WARN:
      rtWarn("The allowed number of clashing atom pairs in any given structure is invalid (" +
             std::to_string(clash_pair_tolerance) + ") and will be reset to zero.",
             "ConformerControls", "validateClashCounts");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    clash_pair_tolerance = 0;
  }
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validateGroupingMethod() {
  if (strcmpCased(group_method, "system", CaseSensitivity::NO) == false &&
      strcmpCased(group_method, "source", CaseSensitivity::NO) == false &&
      strcmpCased(group_method, "sys", CaseSensitivity::NO) == false &&
      strcmpCased(group_method, "topology", CaseSensitivity::NO) == false &&
      strcmpCased(group_method, "label", CaseSensitivity::NO) == false) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Invalid system grouping method \"" + group_method + "\".", "ConformerControls",
            "validateGroupingMethod");
    case ExceptionResponse::WARN:
      rtWarn("Invalid system grouping method \"" + group_method + "\".  System grouping will be "
             "set to \"" + std::string(default_conf_output_grouping) + "\".", "ConformerControls",
             "validateGroupingMethod");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    group_method = std::string(default_conf_output_grouping);
  }
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validateSamplingIntensity() {
  if (strcmpCased(sample_effort, "minimal", CaseSensitivity::NO) == false &&
      strcmpCased(sample_effort, "light", CaseSensitivity::NO) == false &&
      strcmpCased(sample_effort, "heavy", CaseSensitivity::NO) == false &&
      strcmpCased(sample_effort, "exhaustive", CaseSensitivity::NO) == false) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Invalid sampling intensity setting \"" + sample_effort + "\".", "ConformerControls",
            "validateSamplingIntensity");
    case ExceptionResponse::WARN:
      rtWarn("Invalid sampling intensity setting \"" + sample_effort + "\".  Sampling effort will "
             "be set to \"" + std::string(default_conf_sampling_effort) + "\".",
             "ConformerControls", "validateSamplingIntensity");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    sample_effort = std::string(default_conf_sampling_effort);
  }
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validateTorsionSnapping(double *snap_setting,
                                                const std::string &desc) const {
  if (*snap_setting < 0.0 || *snap_setting > symbols::twopi) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A snapping threshold of " +
            realToString(*snap_setting, 9, 4, NumberFormat::STANDARD_REAL) + " is invalid for " +
            desc + " torsion angles.", "ConformerControls", "validateTorsionSnapping");
    case ExceptionResponse::WARN:
      rtWarn("A snapping threshold of " +
             realToString(*snap_setting, 9, 4, NumberFormat::STANDARD_REAL) + " is invalid for " +
             desc + " torsion angles and will be clamped in the range [ 0.0, 360.0 ) degrees.",
             "ConformerControls", "validateTorsionSnapping");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    if (*snap_setting < 0.0) {
      *snap_setting = 0.0;
    }
    else {
      *snap_setting = symbols::twopi;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::validateTorsionAdjustmentProtocol() {
  if (strcmpCased(adjustment_method, "adjust", CaseSensitivity::NO) == false &&
      strcmpCased(adjustment_method, "restrict", CaseSensitivity::NO) == false &&
      strcmpCased(adjustment_method, "cluster", CaseSensitivity::NO) == false &&
      strcmpCased(adjustment_method, "none", CaseSensitivity::NO) == false) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Invalid torsion angle adjustment protocol \"" + adjustment_method + "\".",
            "ConformerControls", "validateTorsionAdjustmentProtocol");
    case ExceptionResponse::WARN:
      rtWarn("Invalid torsion angle adjustment protocol \"" + adjustment_method + "\".  The "
             "protocol will be set to \"" + std::string(default_conf_adjustment_method) + "\".",
             "ConformerControls", "validateTorsionAdjustmentProtocol");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    adjustment_method = std::string(default_conf_adjustment_method);
  }
}

//-------------------------------------------------------------------------------------------------
void ConformerControls::processSamplingValues(int *sample_count, const std::string &count_keyword,
                                              std::vector<double> *sample_values,
                                              const std::string &value_keyword,
                                              const double value_stride, const double value_notch,
                                              const NamelistEmulator &t_nml) {
  double* val_ptr;
  
  // If the defaults for the rotation sample values are present and the user has specified a number
  // of rotation samples other than the default of three, there is an incongruity in the array of
  // rotatable bond angle values and the number of values to sample.  Fill that gap by ignoring the
  // default values and selecting a range of values commensurate with the user's specifications.
  // If, instead, the user has specified specific values to sample but not the number of them, make
  // the number equal to the length of the specified array.  If the user has specified neither, the
  // defaults line up.  If the user has specified both, and the number is not the length of the
  // array, raise an exception or warning depending on the policy.
  if (t_nml.getKeywordStatus(count_keyword) == InputStatus::DEFAULT) {
    *sample_count = t_nml.getKeywordEntries(value_keyword);
    sample_values->resize(*sample_count);
    val_ptr = sample_values->data();
    for (int i = 0; i < *sample_count; i++) {
      val_ptr[i] = t_nml.getRealValue(value_keyword, i);
    }
  }
  else {
    *sample_count = t_nml.getIntValue(count_keyword);
    if (t_nml.getKeywordStatus(value_keyword) == InputStatus::DEFAULT) {
      sample_values->resize(*sample_count);
      double curr_val = 180.0;
      double cval_inc = 0.0;
      const int stride_count = round(360.0 / value_stride);
      val_ptr = sample_values->data();
      for (int i = 0; i < *sample_count; i++) {
        val_ptr[i] = curr_val;
        curr_val += value_stride;
        if (i % stride_count == 0) {
          curr_val += cval_inc;
          cval_inc = -(cval_inc + value_notch);
        }
      }
    }
    else {
      if (t_nml.getKeywordEntries(value_keyword) != *sample_count) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("The number of samples specified (" + std::to_string(*sample_count) + ") does not "
                "match the number of specific sampling values provided (" +
                std::to_string(t_nml.getKeywordEntries(value_keyword)) + ") for \"" +
                count_keyword + "\" / \"" + value_keyword + "\".", "ConformerControls");
        case ExceptionResponse::WARN:
          rtWarn("The number of samples specified (" + std::to_string(*sample_count) + ") does "
                 "not match the number of specific sampling values provided (" +
                 std::to_string(t_nml.getKeywordEntries(value_keyword)) + ") for \"" +
                 count_keyword + "\" / \"" + value_keyword + "\".  Values from the list will be "
                 "taken.", "ConformerControls");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }
      *sample_count = t_nml.getKeywordEntries(value_keyword);
      sample_values->resize(*sample_count);
      val_ptr = sample_values->data();
      for (int i = 0; i < *sample_count; i++) {
        val_ptr[i] = t_nml.getRealValue(value_keyword, i);
      }
    }
  }

  // Convert all rotation sampling values to radians for compatibility with internal units.
  // Re-image the values for clarity.  By this point, the length of the values array and the
  // recorded count should be the same.
  val_ptr = sample_values->data();
  for (size_t i = 0; i < *sample_count; i++) {
    val_ptr[i] *= symbols::pi / 180.0;
    val_ptr[i] = imageValue(val_ptr[i], symbols::twopi, ImagingMethod::MINIMUM_IMAGE);    
  }
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator conformerInput(const TextFile &tf, int *start_line, bool *found,
                                const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("conformer", CaseSensitivity::AUTOMATIC, policy,
                         "Collects instructions for conformer sampling in STORMM.");
  const std::string core_help("Define a subset of atoms that will be prevented from moving, or "
                              "constrained to move in proximity to their initial positions.");
  const std::vector<std::string> core_keys_help = {
    "A data item in a system's SD file to seek out, which will contain one of the following: a "
    "list of atom indices, a list of atom names, or an atom mask specific to that file.  If "
    "the coordinates are not present in an SD file or the SD file does not contain such a data "
    "item, no core atoms will be defined.",
    "An atom mask defining the core atoms.  If given, this will supercede any data items in SD "
    "files and apply to all systems in the calculation.",
    "The restraint penalty, in kcal/mol-Angstrom^2, defining a harmonic repulsive potential that "
    "pushes core atoms away from their initial coordinates.",
    "Alias for 'rk2' in this context.  If both are supplied, the value of 'rk2' will take "
    "precedence.",
    "The harmonic restraint stiffness, in units of kcal/mol-Angstrom^2, preventing atoms from "
    "wandering away from their initial positions.",
    "Alias for 'rk3'.  If both are defined, the value of 'rk3' will take precedence.",
    "Alias for 'rk2' and 'rk3', in units of kcal/mol-Angstrom^2.  This value will apply to both "
    "stiffness constants, but specifying either 'rk2', 'repulsion', "
    "'rk3', or 'attraction' will override the effect.",
    "The distance, in Angstroms, at which to stop applying a repulsive potential pushing atoms "
    "away from their initial positions.",
    "Alias for 'r2', in units of Angstroms as all core restraints are distance-based.  If 'r2' is "
    "supplied, that value will take precedence.",
    "The distance, in Angstroms, at which to begin applying an attractive potential that keeps "
    "atoms from wandering away from their initial positions.",
    "Alias for 'r3', in units of Angstroms.  If 'r3' is supplied, that value will take "
    "precedence." };
  t_nml.addKeyword(NamelistElement("core_mask",
                                   { "data_item", "atoms", "rk2", "repulsion", "rk3", "attraction",
                                     "stiffness", "r2", "demand", "r3", "grace" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::REAL, NamelistType::REAL, NamelistType::REAL,
                                     NamelistType::REAL, NamelistType::REAL, NamelistType::REAL,
                                     NamelistType::REAL, NamelistType::REAL, NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), "16.0", std::string(""),
                                     "0.0", std::string(""), "0.0", "0.0" },
                                   DefaultIsObligatory::NO,
                                   InputRepeats::NO, core_help, core_keys_help,
                                   { KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL }));
  t_nml.addKeyword(NamelistElement("anchor_conf", NamelistType::STRING, std::string("")));
  t_nml.addKeyword(NamelistElement("sample_chirality", NamelistType::STRING,
                                   std::string(default_conf_chirality)));
  t_nml.addKeyword(NamelistElement("sample_cis_trans", NamelistType::STRING,
                                   std::string(default_conf_cis_trans)));
  t_nml.addKeyword(NamelistElement("prevent_hbonds", NamelistType::STRING,
                                   std::string(default_conf_stop_hbonds)));
  t_nml.addKeyword(NamelistElement("running_states", NamelistType::INTEGER,
                                   std::to_string(default_conf_running_states)));
  t_nml.addKeyword(NamelistElement("final_states", NamelistType::INTEGER,
                                   std::to_string(default_conf_final_states)));
  t_nml.addKeyword(NamelistElement("rotation_sample_count", NamelistType::INTEGER,
                                   std::to_string(default_conf_rotation_samples)));
  t_nml.addKeyword(NamelistElement("max_rotatable_bonds", NamelistType::INTEGER,
                                   std::to_string(default_conf_max_rotatable_bonds)));
  t_nml.addKeyword(NamelistElement("cis_trans_sample_count", NamelistType::INTEGER,
                                   std::to_string(default_conf_cis_trans_samples)));
  t_nml.addKeyword(NamelistElement("rotation_sample", NamelistType::REAL,
                                   std::string(default_conf_rotation_set_zero),
                                   DefaultIsObligatory::NO, InputRepeats::YES));
  t_nml.addKeyword(NamelistElement("cis_trans_sample", NamelistType::REAL,
                                   std::string(default_conf_cis_trans_set_zero),
                                   DefaultIsObligatory::NO, InputRepeats::YES));
  t_nml.addKeyword(NamelistElement("max_seeding_attempts", NamelistType::INTEGER,
                                   std::to_string(default_conf_max_seeding_attempts)));
  t_nml.addKeyword(NamelistElement("clash_pair_tol", NamelistType::INTEGER,
                                   std::to_string(default_conf_clash_pairs)));
  t_nml.addKeyword(NamelistElement("trial_limit", NamelistType::INTEGER,
                                   std::to_string(default_conf_max_system_trials)));
  t_nml.addKeyword(NamelistElement("local_trial_limit", NamelistType::INTEGER,
                                   std::to_string(default_conf_sample_trials)));
  t_nml.addKeyword(NamelistElement("rmsd_tol", NamelistType::REAL,
                                   std::to_string(default_conf_rmsd_tolerance)));
  t_nml.addKeyword(NamelistElement("rotamer_adjustment", NamelistType::STRING,
                                   std::string(default_conf_adjustment_method)));
  t_nml.addKeyword(NamelistElement("grouping", NamelistType::STRING,
                                   std::string(default_conf_output_grouping)));
  t_nml.addKeyword(NamelistElement("effort", NamelistType::STRING,
                                   std::string(default_conf_sampling_effort)));
  t_nml.addHelp("core_mask", "Atom mask for common core atoms.  These atoms will be held in a "
                "rigid configuration during energy minimization and other sampling operations.");
  t_nml.addHelp("anchor_conf", "An exemplary ligand structure used in aligning the common core "
                "atoms.");
  t_nml.addHelp("sample_chirality", "Sample chiral states of identifiable chiral centers.  "
                "Specify 'yes' / 'true' to sample or 'no' / 'false' to decline.");
  t_nml.addHelp("sample_cis_trans", "Sample cis and trans states of double bonds.  Specify "
                "'yes' / 'true' to sample or 'no' / 'false' to decline.");
  t_nml.addHelp("prevent_hbonds", "A quick way to have STORMM prevent hydrogen bonding between "
                "donors and acceptor atoms that it can identify in the molecule(s).  This will "
                "establish a restraint ensemble for each case with default parameters to prevent "
                "donor and acceptor pairs from coming too close.");
  t_nml.addHelp("running_states", "Number of energy-minimizations to carry out at one time in "
                "order to generate a smaller set of final states.");
  t_nml.addHelp("final_states", "Number of final, energy-minimized states to accept as unique "
                "conformers.");
  t_nml.addHelp("rotation_sample_count", "Number of samples to apply to each rotatable bond.  "
                "Locations for the samples will be chosen based on the detected minima along each "
                "bond's rotation profile.");
  t_nml.addHelp("max_rotatable_bonds", "The maximum number of rotatable bonds to explicitly "
                "sample.  This quickly runs into a combinatorial problem, but there is a "
                "guardrail in the max_system_trials keyword.");
  t_nml.addHelp("cis_trans_sample_count", "Number of samples to apply to each cis-trans isomeric "
                "bond.  Locations for the samples will be chosen based on the detected minima "
                "along each bond's rotation profile.");
  t_nml.addHelp("rotation_sample", "A specific value of rotatable bonds to sample, given in units "
                "of degrees.  Repeated entries are accepted.  The default values sample basic "
                "rotamers expected for a bond between sp3 carbon atoms.");
  t_nml.addHelp("cis_trans_sample", "A specific value of cis-trans isomeric bonds to sample, "
                "given in units of degrees.  Repeated entries are accepted.  The default values "
                "sample basic cis- and trans- configurations for a bond between two sp2 carbon "
                "atoms.");
  t_nml.addHelp("max_seeding_attempts", "If a conformer's initial configuration contains a clash "
                "between atoms (their van-der Waals radii are violated), randomization of the "
                "configuration will occur for this number of attempts.  If, after exhausting this "
                "provision, a stable configuration still cannot be found, the input configuration "
                "will be accepted as the initial coordinates for subsequent energy "
                "minimizations.");
  t_nml.addHelp("clash_pair_tol", "In order to declare a seeded conformer 'clashing', the "
                "internal conflicts must be great enough that guided minimization with a "
                "softcore potential would not be able to resolve the conflicts.  This is the "
                "number of clashing atom pairs that will be assumed resolvable by the softcore "
                "energy minimization.");
  t_nml.addHelp("trial_limit", "The maximum number of trials that will be made for each "
                "system.  Explicit sampling of chirality, cis-trans isomers, and then rotatable "
                "bonds will proceed in that priority, but the maximum number of sampled states "
                "will be cut off at this value.");
  t_nml.addHelp("local_trial_limit", "The maximum number of trials, per system, to employ when "
                "doing local sampling of each system to ascertain its preferred torsion angles "
                "and (possibly) accessible chiral states.");
  t_nml.addHelp("rmsd_tol", "Positional, mass-weighted root-mean squared deviation between "
                "conformers required to declare uniqueness.  Units of Angstroms.");
  t_nml.addHelp("rotamer_adjustment", "Method for adjusting rotamer settings in light of known "
                "energy-minimized structures, or other examples of each ligand system.");
  t_nml.addHelp("grouping", "An indication of how to group systems when selecting the best "
                "conformers for output.  Acceptable values include \"system\", \"source\", or "
                "\"sys\" (produce outputs for each system defined by a '-sys' keyword in the "
                "&files namelist), \"topology\" (group all systems sharing the same topology "
                "file), or \"label\" ( group all systems marked with the same label group, as "
                "defined in the &files namelist).");
  t_nml.addHelp("effort", "Specifies the approximate level of effort, measured in terms of the "
                "number of independent minimization attempts starting from unique configurations, "
                "that will be applied to find optimal structures in each molecule.  Options "
                "include MINIMAL, LIGHT, HEAVY, and EXHAUSTIVE, with higher and higher orders "
                "of the combinatorial space being sampled with each setting.");

  // Add more default values to the rotation_sample and cis_trans_sample keywords.
  t_nml.addDefaultValue("rotation_sample", std::string(default_conf_rotation_set_one));
  t_nml.addDefaultValue("rotation_sample", std::string(default_conf_rotation_set_two));
  t_nml.addDefaultValue("cis_trans_sample", std::string(default_conf_cis_trans_set_one));
  
  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  return t_nml;
}

} // namespace namelist
} // namespace stormm
