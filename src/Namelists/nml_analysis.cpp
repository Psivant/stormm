#include "copyright.h"
#include "Constants/symbol_values.h"
#include "Parsing/parse.h"
#include "namelist_common.h"
#include "namelist_element.h"
#include "nml_analysis.h"

namespace stormm {
namespace namelist {

using parse::NumberFormat;
using parse::realToString;
using parse::TextOrigin;
using symbols::pi;

//-------------------------------------------------------------------------------------------------
AnalysisControls::AnalysisControls(const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    policy{policy_in},
    identifier{std::string("")},
    multiple_analyses_present{false},
    general_statistical_blocks{default_analysis_stat_blocks},
    general_sample_frequency{default_analysis_sample_interval},
    general_initiation_step{default_analysis_init_step},
    hbond_mask_count{0},
    hbond_statistical_blocks{default_analysis_stat_blocks},
    hbond_sample_frequency{default_analysis_sample_interval},
    hbond_initiation_step{default_analysis_init_step},
    hbond_mask_i{}, hbond_mask_ii{}, hbond_mask_labels{},
    hbond_variable_base{default_hbond_variable_base},
    hbond_variable_composite{std::string(default_hbond_variable_base) + std::string("_") +
                             std::string(default_analysis_identifier) + std::to_string(0)},
    hbond_max_separation{default_hbond_da_max_separation},
    hbond_min_angle{default_hbond_dha_min_angle},
    hbond_candidacy_proximity{default_hbond_initial_proximity},
    hbond_occupancy_threshold{default_hbond_occ_threshold},
    hbond_activity_threshold{default_hbond_act_threshold},
    nml_transcript{"analysis"}
{
  // Load in a blank namelist so that certain keywords will be present, as if this were the means
  // by which the data was loaded.
  std::string tfs("&analysis\n&end\n");
  TextFile tf(tfs, TextOrigin::RAM);
  int start_line = 0;
  bool found;
  nml_transcript = analysisInput(tf, &start_line, &found, ExceptionResponse::SILENT);
}

//-------------------------------------------------------------------------------------------------
AnalysisControls::AnalysisControls(const TextFile &tf, int *start_line, bool *found_nml,
                                   const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    AnalysisControls(policy_in)
{
  NamelistEmulator t_nml = analysisInput(tf, start_line, found_nml, policy, wrap);
  nml_transcript = t_nml;

  // Interpret keywords for general activities
  t_nml.assignVariable(&general_statistical_blocks, "gen_stat_blocks");
  t_nml.assignVariable(&general_sample_frequency, "gen_ntpr");
  t_nml.assignVariable(&general_initiation_step, "gen_stepi");
  if (t_nml.getKeywordStatus("tag") != InputStatus::MISSING) {
    setIdentifier(t_nml.getStringValue("tag"));
  }

  // Interpret keywords for hydrogen bond analysis
  hbond_mask_count = t_nml.getKeywordEntries("hbond");
  hbond_mask_i.reserve(hbond_mask_count);
  hbond_mask_ii.reserve(hbond_mask_count);
  for (int i = 0; i < hbond_mask_count; i++) {
    hbond_mask_i.push_back(t_nml.getStringValue("hbond", "mask1", i));
    if (t_nml.getKeywordStatus("hbond", "mask2", i) == InputStatus::USER_SPECIFIED) {
      hbond_mask_ii.push_back(t_nml.getStringValue("hbond", "mask2", i));
    }
    else {
      hbond_mask_ii.push_back(hbond_mask_i.back());
    }
    hbond_mask_labels.push_back(t_nml.getStringValue("hbond", "-label", i));
  }
  t_nml.assignVariable(&hbond_variable_base, "hb_var");
  t_nml.assignVariable(&hbond_max_separation, "hb_range");
  t_nml.assignVariable(&hbond_min_angle, "hb_angle");
  t_nml.assignVariable(&hbond_candidacy_proximity, "hb_init_proximity");
  t_nml.assignVariable(&hbond_statistical_blocks, "hb_stat_blocks");
  t_nml.assignVariable(&hbond_sample_frequency, "hb_ntpr");
  t_nml.assignVariable(&hbond_initiation_step, "hb_stepi");
  t_nml.assignVariable(&hbond_occupancy_threshold, "hb_min_occ");
  t_nml.assignVariable(&hbond_activity_threshold, "hb_min_act");
  
  // Convert some parameters from the input units, such as degrees, into the internal units (e.g.
  // radians).
  hbond_min_angle *= pi / 180.0;

  // Apply general analysis directives
  applyGeneralStatBlocks();
  applyGeneralSampleFrequency();
  applyGeneralInitiationStep();
  
  // Validate inputs
  validateStatBlocks(&general_statistical_blocks, "general analyses");
  validateSamplingFrequency(&general_sample_frequency, "general analyses");
  validateInitiationStep(&general_initiation_step, "general analyses");
  validateHBondMaxSeparation();
  validateHBondMinAngle();
  validateHBondCandidacyProximity();
  validateStatBlocks(&hbond_statistical_blocks, "hydrogen bond analysis");
  validateSamplingFrequency(&hbond_sample_frequency, "hydrogen bond analysis");
  validateInitiationStep(&hbond_initiation_step, "hydrogen bond analysis");
}

//-------------------------------------------------------------------------------------------------
const std::string& AnalysisControls::getIdentifier() const {
  return identifier;
}

//-------------------------------------------------------------------------------------------------
bool AnalysisControls::multipleAnalysesPresent() const {
  return multiple_analyses_present;
}

//-------------------------------------------------------------------------------------------------
int AnalysisControls::getGeneralStatBlocks() const {
  return general_statistical_blocks;
}
  
//-------------------------------------------------------------------------------------------------
int AnalysisControls::getGeneralSamplingFrequency() const {
  return general_sample_frequency;
}
  
//-------------------------------------------------------------------------------------------------
int AnalysisControls::getGeneralInitiationStep() const {
  return general_initiation_step;
}

//-------------------------------------------------------------------------------------------------
int AnalysisControls::getHBondMaskCount() const {
  return hbond_mask_count;
}

//-------------------------------------------------------------------------------------------------
double AnalysisControls::getHBondMaxSeparation() const {
  return hbond_max_separation;
}

//-------------------------------------------------------------------------------------------------
double AnalysisControls::getHBondMinAngle() const {
  return hbond_min_angle;
}

//-------------------------------------------------------------------------------------------------
double AnalysisControls::getHBondCandidacyProximity() const {
  return hbond_candidacy_proximity;
}

//-------------------------------------------------------------------------------------------------
int AnalysisControls::getHBondStatBlocks() const {
  return hbond_statistical_blocks;
}

//-------------------------------------------------------------------------------------------------
int AnalysisControls::getHBondSamplingFrequency() const {
  return hbond_sample_frequency;
}

//-------------------------------------------------------------------------------------------------
int AnalysisControls::getHBondInitiationStep() const {
  return hbond_initiation_step;
}

//-------------------------------------------------------------------------------------------------
double AnalysisControls::getHBondOccupancyThreshold() const {
  return hbond_occupancy_threshold;
}

//-------------------------------------------------------------------------------------------------
double AnalysisControls::getHBondActivityThreshold() const {
  return hbond_activity_threshold;
}

//-------------------------------------------------------------------------------------------------
const std::string& AnalysisControls::getHBondMask(const int mask_index, const int array_id) const {
  if (mask_index < 0 || mask_index >= hbond_mask_count) {
    rtErr("Mask index " + std::to_string(mask_index) + " is invalid for a collection of " +
          std::to_string(hbond_mask_count) + " mask pairs.", "AnalysisControls", "getHBondMask");
  }
  if (array_id == 0) {
    return hbond_mask_i[mask_index];
  }
  else if (array_id == 1) {
    return hbond_mask_ii[mask_index];
  }
  else {
    rtErr("Array index " + std::to_string(array_id) + " is invalid for a hydrogen bonding matrix.",
          "AnalysisControls", "getHBondMask");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const std::string& AnalysisControls::getHBondMaskLabel(const int label_index) const {
  if (label_index < 0 || label_index >= hbond_mask_count) {
    rtErr("Label index " + std::to_string(label_index) + " is invalid for a collection of " +
          std::to_string(hbond_mask_count) + " mask pairs.", "AnalysisControls",
          "getHBondMaskLabel");
  }
  return hbond_mask_labels[label_index];
}

//-------------------------------------------------------------------------------------------------
const std::string& AnalysisControls::getHBondVariableBase() const {
  if (multiple_analyses_present) {
    return hbond_variable_composite;
  }
  else {
    return hbond_variable_base;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const NamelistEmulator& AnalysisControls::getTranscript() const {
  return nml_transcript;
}
  
//-------------------------------------------------------------------------------------------------
void AnalysisControls::setIdentifier(const std::string &identifier_in) {
  identifier = identifier_in;

  // Set composite variable base names
  hbond_variable_composite = hbond_variable_base + "_" + identifier;
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::multipleAnalysesFound(const bool found) {
  multiple_analyses_present = found;

  // Set composite variable base names with the known identifiers
  hbond_variable_composite = hbond_variable_base + "_" + identifier;
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::setGeneralStatBlocks(const int blocks_in) {
  general_statistical_blocks = blocks_in;
  validateStatBlocks(&general_statistical_blocks, "general analyses");
  applyGeneralStatBlocks();
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::setGeneralSamplingFrequency(const int frequency_in) {
  general_sample_frequency = frequency_in;
  validateSamplingFrequency(&general_sample_frequency, "general analyses");
  applyGeneralSampleFrequency();
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::setGeneralInitiationStep(const int step_in) {
  general_initiation_step = step_in;
  validateInitiationStep(&general_initiation_step, "general analyses");
  applyGeneralInitiationStep();
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::setHBondMatrix(const std::string &mask_i, const std::string &mask_ii,
                                      const std::string &label) {
  hbond_mask_i.push_back(mask_i);
  hbond_mask_ii.push_back(mask_ii);
  hbond_mask_labels.push_back(label);
}
  
//-------------------------------------------------------------------------------------------------
void AnalysisControls::setHBondMaxSeparation(const double separation_in) {
  hbond_max_separation = separation_in;
}
  
//-------------------------------------------------------------------------------------------------
void AnalysisControls::setHBondMinAngle(const double angle_in) {
  hbond_min_angle = angle_in * pi / 180.0;
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::setHBondInitialProximity(const double proximity_in) {
  hbond_candidacy_proximity = proximity_in;
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::setHBondStatBlocks(const int blocks_in) {
  hbond_statistical_blocks = blocks_in;
  validateStatBlocks(&hbond_statistical_blocks, "hydrogen bond analysis");
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::setHBondSamplingFrequency(const int frequency_in) {
  hbond_sample_frequency = frequency_in;
  validateSamplingFrequency(&hbond_sample_frequency, "hydrogen bond analysis");
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::setHBondInitiationStep(const int step_in) {
  hbond_initiation_step = step_in;
  validateInitiationStep(&hbond_initiation_step, "hydrogen bond analysis");
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::setHBondOccupancyThreshold(const double occupancy_in) {
  hbond_occupancy_threshold = occupancy_in;
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::setHBondActivityThreshold(const double activity_in) {
  hbond_activity_threshold = activity_in;
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::setHBondVariableBase(const std::string &varname_in) {
  hbond_variable_base = varname_in;
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::validateStatBlocks(int *setting, const char* desc) {

  // Special case: setting the block count to zero may be interpreted as foregoing statistical
  // block averaging.
  if (*setting == 0) {
    *setting = 1;
  }
  if (*setting < 1 || *setting == 2) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A block averaging strategy with " + std::to_string(*setting) + " blocks is invalid "
            "for " + std::string(desc) + ".", "AnalysisControls", "validateStatBlocks");
    case ExceptionResponse::WARN:
      rtWarn("A block averaging strategy with " + std::to_string(*setting) + " blocks is invalid "
             "for " + std::string(desc) + ".  Block averaging will not be implemented in this "
             "case.", "AnalysisControls", "validateStatBlocks");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    *setting = 1;
  }
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::validateSamplingFrequency(int *setting, const char* desc) {
  if (*setting < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A sampling frequency of " + std::to_string(*setting) + " steps is invalid for " +
            std::string(desc) + ".", "AnalysisControls", "validateSamplingFrequency");
    case ExceptionResponse::WARN:
      rtWarn("A sampling frequency of " + std::to_string(*setting) + " steps is invalid for " +
             std::string(desc) + ".  The frequency will be set to zero and no sampling will "
             "occur.", "AnalysisControls", "validateSamplingFrequency");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    *setting = 0;
  }
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::validateInitiationStep(int *setting, const char* desc) {
  if (*setting < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("An initiation step number of " + std::to_string(*setting) + " is invalid for " +
            std::string(desc) + ".", "AnalysisControls", "validateInitiationStep");
    case ExceptionResponse::WARN:
      rtWarn("A initiation step number of " + std::to_string(*setting) + " is invalid for " +
             std::string(desc) + ".  The value will be reset to zero and analysis will occur "
             "throughout the simulation.", "AnalysisControls", "validateInitiationStep");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    *setting = 0;
  }
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::validateHBondMaxSeparation() {
  if (hbond_max_separation < 1.0 || hbond_max_separation >= 20.0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A maximum separation of " +
            realToString(hbond_max_separation, 9, 4, NumberFormat::STANDARD_REAL) +
            " is unreasonable for detecting hydrogen bonds.", "AnalysisControls",
            "valdiateHBondMaxSeparation");
    case ExceptionResponse::WARN:
      rtWarn("A maximum separation of " +
             realToString(hbond_max_separation, 9, 4, NumberFormat::STANDARD_REAL) +
             " is unreasonable for detecting hydrogen bonds.  The default value of " +
             realToString(hbond_max_separation, 9, 4, NumberFormat::STANDARD_REAL) + " will be "
             "applied.", "AnalysisControls", "validateHBondMaxSeparation");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    hbond_max_separation = default_hbond_da_max_separation;
  }
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::validateHBondMinAngle() {
  if (hbond_min_angle >= 2.0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A minimum angle of " +
            realToString(hbond_min_angle * 180.0 / pi, 9, 4, NumberFormat::STANDARD_REAL) +
            " is unattainable for three points.", "AnalysisControls", "valdiateHBondMinAngle");
    case ExceptionResponse::WARN:
      rtWarn("A minimum angle of " +
             realToString(hbond_min_angle * 180.0 / pi, 9, 4, NumberFormat::STANDARD_REAL) +
             " is unattainable for three points.  The default of " +
             realToString(default_hbond_dha_min_angle * 180.0 / pi, 9, 4,
                          NumberFormat::STANDARD_REAL) + " will be applied.", "AnalysisControls",
             "valdiateHBondMinAngle");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    hbond_min_angle = default_hbond_dha_min_angle * pi / 180.0;
  }
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::validateHBondCandidacyProximity() {
  if (hbond_candidacy_proximity < 1.0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("An initial separation of no more than " +
            realToString(hbond_candidacy_proximity, 9, 4, NumberFormat::STANDARD_REAL) +
            " is unreasonable and would preclude too many potential hydrogen bonds from being "
            "searched in a simulation.", "AnalysisControls", "validateHBondCandidacyProximity");
    case ExceptionResponse::WARN:
      rtErr("An initial separation of no more than " +
            realToString(hbond_candidacy_proximity, 9, 4, NumberFormat::STANDARD_REAL) +
            " is unreasonable and would preclude too many potential hydrogen bonds from being "
            "searched in a simulation.  The default of " +
            realToString(default_hbond_initial_proximity, 9, 4, NumberFormat::STANDARD_REAL) +
            " will be applied.", "AnalysisControls", "validateHBondCandidacyProximity");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    hbond_candidacy_proximity = default_hbond_initial_proximity;
  }
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::applyGeneralStatBlocks() {
  if (nml_transcript.getKeywordStatus("gen_stat_blocks") == InputStatus::USER_SPECIFIED) {
    if (nml_transcript.getKeywordStatus("hb_stat_blocks") == InputStatus::DEFAULT) {
      hbond_statistical_blocks = general_statistical_blocks;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::applyGeneralSampleFrequency() {
  if (nml_transcript.getKeywordStatus("gen_ntpr") == InputStatus::USER_SPECIFIED) {
    if (nml_transcript.getKeywordStatus("hb_ntpr") == InputStatus::DEFAULT) {
      hbond_sample_frequency = general_sample_frequency;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AnalysisControls::applyGeneralInitiationStep() {
  if (nml_transcript.getKeywordStatus("gen_stepi") == InputStatus::USER_SPECIFIED) {
    if (nml_transcript.getKeywordStatus("hb_stepi") == InputStatus::DEFAULT) {
      hbond_initiation_step = general_initiation_step;
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
NamelistEmulator analysisInput(const TextFile &tf, int *start_line, bool *found,
                               const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("analysis", CaseSensitivity::AUTOMATIC, policy, "Wraps directives needed "
                         "to analyze molecular dynamics simulations in-flight.");

  // Add namelist-specific keywords for general analyses
  t_nml.addKeyword("tag", NamelistType::STRING, std::string(""));
  t_nml.addKeyword("gen_stat_blocks", NamelistType::INTEGER,
                   std::to_string(default_analysis_stat_blocks));
  t_nml.addKeyword("gen_ntpr", NamelistType::INTEGER,
                   std::to_string(default_analysis_sample_interval));
  t_nml.addKeyword("gen_stepi", NamelistType::INTEGER, std::to_string(default_analysis_init_step));
  
  // Add namelist-specific keywords for hydrogen bond analysis
  t_nml.addKeyword("hb_range", NamelistType::REAL,
                   realToString(default_hbond_da_max_separation, 9, 4,
                                NumberFormat::STANDARD_REAL));
  t_nml.addKeyword("hb_angle", NamelistType::REAL,
                   realToString(default_hbond_dha_min_angle, 9, 4,
                                NumberFormat::STANDARD_REAL));
  t_nml.addKeyword("hb_init_proximity", NamelistType::REAL,
                   realToString(default_hbond_initial_proximity, 9, 4,
                                NumberFormat::STANDARD_REAL));
  t_nml.addKeyword("hb_stat_blocks", NamelistType::INTEGER,
                   std::to_string(default_analysis_stat_blocks));
  t_nml.addKeyword("hb_ntpr", NamelistType::INTEGER,
                   std::to_string(default_analysis_sample_interval));
  t_nml.addKeyword("hb_stepi", NamelistType::INTEGER, std::to_string(default_analysis_init_step));
  t_nml.addKeyword("hb_var", NamelistType::STRING, std::string(default_hbond_variable_base));
  const std::string hbond_help("Call for hydrogen bonding analysis between donors and acceptors "
                               "found among up to two distinct groups of atoms.  If only one "
                               "group of atoms is specified, hydrogen bonding within that "
                               "selection will be assessed.  If two groups are specified, "
                               "hydrogen bonds formed between donors or acceptors in the "
                               "first group and acceptors in the second group will be assessed.");
  const std::vector<std::string> hbond_keys_help = {
    "Atom selection for the first (or sole) group, given in ambmask format.",
    "Atom selection for the second group, given in ambmask format.",
    "System label corresponding to one or more systems listed in the &files namelist.  The "
    "default value of 'all' will apply the atom selection criteria and attempt to find groups of "
    "donors and acceptors in every system of the problem." };
  t_nml.addKeyword("hbond", { "mask1", "mask2", "-label" },
                   { NamelistType::STRING, NamelistType::STRING, NamelistType::STRING },
                   { std::string(""), std::string(""), std::string("all") },
                   DefaultIsObligatory::NO, InputRepeats::YES, hbond_help, hbond_keys_help,
                   { KeyRequirement::REQUIRED, KeyRequirement::OPTIONAL,
                     KeyRequirement::OPTIONAL });
  t_nml.addKeyword("hb_min_occ", NamelistType::REAL, std::to_string(default_hbond_occ_threshold));
  t_nml.addKeyword("hb_min_act", NamelistType::REAL, std::to_string(default_hbond_act_threshold));

  // Add general analysis keyword helpmessage
  t_nml.addHelp("tag", "A general identification label to place on all processes within this "
                "&analysis namelist control block.  Multiple &analysis blocks may be included "
                "in a single input file, with unique tags being one way to distinguish their "
                "respective results.  If only one &analysis block is present, the tag can be "
                "left unset, but if multiple &analysis blocks are present with their tags unset "
                "then the program will scan through them all and try to assign a unique tag to "
                "each.");
  t_nml.addHelp("gen_stat_blocks", "The number of statistical blocks into which to divide "
                "vairous analyses.  If specified, this general setting can be superceded by "
                "specific block averaging directives in particular analyses.  The default "
                "behavior is to forego block averaging.");
  t_nml.addHelp("gen_ntpr", "General frequency, expressed in simulation time steps, at which to "
                "perform analysis.  If specified, this general setting can be superceded by "
                "specific sampling frequencies in particular analyses.");
  t_nml.addHelp("gen_stepi", "The simulation time step at which to initiate general analyses.  If "
                "specified, this general setting can be superceded by specific sampling "
                "frequencies in particular analyses.");

  // Add namelist help messages for hydrogen bonding analysis
  t_nml.addHelp("hb_range", "The maximum separation, in units of Angstroms, at which a hydrogen "
                "bond donor and acceptor could be considered to be participating in a hydrogen "
                "bond.  This is not the only criteria for confirming the presence of a hydrogen "
                "bond: necessary, but not sufficient.");
  t_nml.addHelp("hb_angle", "The minimum angle, given in degrees, at which a hydrogen bond donor, "
                "its attached polar hydrogen (proton), and a hydrogen bond proton acceptor could "
                "form to be considered to form a hydrogen bond.  This criterion, in conjunction "
                "with the maximum viable distance, establishes the existence of a hydrogen bond "
                "between atoms in the groups specified in an hbond keyword entry.");
  t_nml.addHelp("hb_init_proximity", "Proximity between hydrogen bond donor and acceptor atoms in "
                "the initial configuration needed to activate the search for a hydrogen bond "
                "between the pair throughout the simulation.  Units of Angstroms.");
  t_nml.addHelp("hb_stat_blocks", "The number of statistical blocks into which ti divid hydrogen "
                "bonding analysis.  The default behavior is to forego block averaging.");
  t_nml.addHelp("hb_ntpr", "Sampling frequency, expressed in simulation time steps, at which to "
                "perform hydrogen bond analysis.");
  t_nml.addHelp("hb_stepi", "The simulation time step at which to initiate hydrogen bond "
                "analysis.");
  t_nml.addHelp("hb_var", "The base name with which to output results from hydrogen bonding "
                "analysis.  Multiple matrix-style variables will be created in the output file "
                "with this name, plus unique suffixes.");
  t_nml.addHelp("hb_min_occ", "Occupancy required to trigger reporting of a specific hydrogen "
                "bonding donor :: acceptor arrangement");
  t_nml.addHelp("hb_min_act", "Total average presence of hydrogen bonding required to trigger "
                "the report of a system as displaying significant hydrogen bonding between the "
                "first and second masks");
  
  // Search the input file.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  return t_nml;
}
  
} // namespace namelist
} // namespace stormm
