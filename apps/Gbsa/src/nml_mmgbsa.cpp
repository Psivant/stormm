#include "../../../src/copyright.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Reporting/error_format.h"
#include "nml_mmgbsa.h"

namespace mmgbsa {

using stormm::constants::CaseSensitivity;
using stormm::errors::rtErr;
using stormm::errors::rtWarn;
using stormm::namelist::DefaultIsObligatory;
using stormm::namelist::InputRepeats;
using stormm::namelist::InputStatus;
using stormm::namelist::KeyRequirement;
using stormm::namelist::NamelistElement;
using stormm::namelist::NamelistType;
using stormm::parse::NumberFormat;
using stormm::parse::realToString;
  
//-------------------------------------------------------------------------------------------------
MMGBSAControls::MMGBSAControls(const ExceptionResponse policy_in) :
    policy{policy_in}, ligand_max_poses{0}, best_poses{0},
    averaging{WeightingScheme::FLAT},
    temperature{default_mmgbsa_temperature},
    carveout_cutoff{default_carveout_cutoff},
    carveout_halo{default_carveout_halo},
    proximity_carveout{false},
    receptor_from_complex{false},
    report_depth{EnergyReport::COMPLETE},
    decompose_energy{false},
    print_structures{false},
    structure_format{PrintedPoseFormat::PDB},
    structure_base_name{std::string(default_structure_base_name)},
    carveout_extra_mask{default_carveout_mask},
    carveout_rotating_group_limit{default_rotating_group_limit},
    carveout_whole_residues{false},
    carveout_ligand_reference{std::string("all")},
    nml_transcript{"mmgbsa"}
{}
  
//-------------------------------------------------------------------------------------------------
MMGBSAControls::MMGBSAControls(const TextFile &tf, int *start_line, bool *found_nml,
                               const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    MMGBSAControls(policy_in)
{
  const NamelistEmulator t_nml = mmgbsaInput(tf, start_line, found_nml, policy, wrap);
  nml_transcript = t_nml;
  ligand_max_poses = t_nml.getIntValue("ligand_max_poses");
  best_poses = t_nml.getIntValue("best_poses");
  averaging = translateWeightingScheme(t_nml.getStringValue("weights"));
  report_depth = translateEnergyReport(t_nml.getStringValue("depth"));
  temperature = t_nml.getRealValue("temperature");
  validateDistanceCutoff(&carveout_cutoff, t_nml.getRealValue("carveout", "cutoff"),
                         "MMGBSAControls");
  validateDistanceCutoff(&carveout_halo, t_nml.getRealValue("carveout", "halo"),
                         "MMGBSAControls");
  carveout_extra_mask = t_nml.getStringValue("carveout", "mask");
  if (t_nml.getKeywordStatus("carveout", "cutoff") == InputStatus::DEFAULT &&
      t_nml.getKeywordStatus("carveout", "mask") == InputStatus::USER_SPECIFIED) {
    proximity_carveout = false;

    // Checks on input sanity: if a rotating group size is specified, that is an input error as
    // no proximity-based atom selection is in place.
    if (t_nml.getIntValue("carveout", "add_rotg") > 0) {
      const std::string err_msg("A nontrivial size of rotating groups has been specified (" +
                                std::to_string(t_nml.getIntValue("carveout", "add_rotg")) +
                                ").  This is intended to include all particles of rotatable "
                                "groups smaller than the stated limit in the receptor carveout "
                                "(binding site) selection if one or more of their atoms are "
                                "affected by the distance-based cutoff criterion.  However, no "
                                "distance criterion was specified, and with an atom mask in play "
                                "('" + t_nml.getStringValue("carveout", "mask") + "'), the atom "
                                "mask will be the sole determinant of the carveout selection.");
      switch (policy_in) {
      case ExceptionResponse::DIE:
        rtErr(err_msg, "MMGBSAControls");
      case ExceptionResponse::WARN:
        rtWarn(err_msg, "MMGBSAControls");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
    if (t_nml.getBoolValue("whole_residues")) {
      const std::string err_msg("Specifying that whole residues be included (if any of their "
                                "atoms are affected by other selection criteria) is only valid, "
                                "and only applies to, atoms selected by a distance-based cutoff "
                                "criterion.  No such distance has been specified.");
      switch (policy_in) {
      case ExceptionResponse::DIE:
        rtErr(err_msg, "MMGBSAControls");
      case ExceptionResponse::WARN:
        rtWarn(err_msg, "MMGBSAControls");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
  }
  else {
    proximity_carveout = true;
  }
  receptor_from_complex = t_nml.getBoolValue("rcpt_from_cmpx");
  carveout_rotating_group_limit = t_nml.getIntValue("carveout", "add_rotg");
  carveout_whole_residues = t_nml.getBoolValue("whole_residues");
  decompose_energy = t_nml.getBoolValue("nrg_decomp");
  print_structures = t_nml.getBoolValue("final_structures");
  structure_format = translatePrintedPoseFormat(t_nml.getStringValue("structure_fmt"));
  structure_base_name = t_nml.getStringValue("structure_base_name");
}

//-------------------------------------------------------------------------------------------------
int MMGBSAControls::getPoseReportCount() const {
  return ligand_max_poses;
}

//-------------------------------------------------------------------------------------------------
WeightingScheme MMGBSAControls::getWeightingScheme() const {
  return averaging;
}

//-------------------------------------------------------------------------------------------------
double MMGBSAControls::getTemperature() const {
  return temperature;
}
  
//-------------------------------------------------------------------------------------------------
double MMGBSAControls::getCarveoutCutoff() const {
  return carveout_cutoff;
}
  
//-------------------------------------------------------------------------------------------------
double MMGBSAControls::getCarveoutHalo() const {
  return carveout_halo;
}
  
//-------------------------------------------------------------------------------------------------
bool MMGBSAControls::proximityCarveout() const {
  return proximity_carveout;
}

//-------------------------------------------------------------------------------------------------
EnergyReport MMGBSAControls::getEnergyReportDepth() const {
  return report_depth;
}

//-------------------------------------------------------------------------------------------------
bool MMGBSAControls::printEnergyDecomposition() const {
  return decompose_energy;
}

//-------------------------------------------------------------------------------------------------
bool MMGBSAControls::printFinalStructures() const {
  return print_structures;
}

//-------------------------------------------------------------------------------------------------
PrintedPoseFormat MMGBSAControls::getPrintedStructureFormat() const {
  return structure_format;
}

//-------------------------------------------------------------------------------------------------
const std::string& MMGBSAControls::getOutputStructureBase() const {
  return structure_base_name;
}
  
//-------------------------------------------------------------------------------------------------
const std::string& MMGBSAControls::getCarveoutExtraMask() const {
  return carveout_extra_mask;
}

//-------------------------------------------------------------------------------------------------
int MMGBSAControls::getRotatingGroupLimit() const {
  return carveout_rotating_group_limit;
}

//-------------------------------------------------------------------------------------------------
bool MMGBSAControls::completeResidueCarveout() const {
  return carveout_whole_residues;
}

//-------------------------------------------------------------------------------------------------
bool MMGBSAControls::takeReceptorFromComplex() const {
  return receptor_from_complex;
}

//-------------------------------------------------------------------------------------------------
const std::string& MMGBSAControls::getCarveoutLigandReference() const {
  return carveout_ligand_reference;
}
  
//-------------------------------------------------------------------------------------------------
void MMGBSAControls::setPoseReportCount(const int ligand_max_poses_in) {
  ligand_max_poses = ligand_max_poses_in;
}

//-------------------------------------------------------------------------------------------------
void MMGBSAControls::setWeightingScheme(const std::string &averaging_in) {
  averaging = translateWeightingScheme(averaging_in);
}

//-------------------------------------------------------------------------------------------------
void MMGBSAControls::setWeightingScheme(const WeightingScheme averaging_in) {
  averaging = averaging_in;
}

//-------------------------------------------------------------------------------------------------
void MMGBSAControls::setTemperature(const double temperature_in) {
  temperature = temperature_in;
  validateTemperature();
}

//-------------------------------------------------------------------------------------------------
void MMGBSAControls::setEnergyReportDepth(const std::string &report_depth_in) {
  report_depth = translateEnergyReport(report_depth_in);
}

//-------------------------------------------------------------------------------------------------
void MMGBSAControls::setEnergyReportDepth(const EnergyReport report_depth_in) {
  report_depth = report_depth_in;
}

//-------------------------------------------------------------------------------------------------
void MMGBSAControls::setEnergyDecompositionReporting(const bool decompose_energy_in) {
  decompose_energy = decompose_energy_in;
}

//-------------------------------------------------------------------------------------------------
void MMGBSAControls::setStructurePrinting(const bool print_structures_in) {
  print_structures = print_structures_in;
}

//-------------------------------------------------------------------------------------------------
void MMGBSAControls::setFinalStructureFormat(const std::string &structure_format_in) {
  structure_format = translatePrintedPoseFormat(structure_format_in);
}

//-------------------------------------------------------------------------------------------------
void MMGBSAControls::setFinalStructureFormat(const PrintedPoseFormat structure_format_in) {
  structure_format = structure_format_in;
}

//-------------------------------------------------------------------------------------------------
void MMGBSAControls::setCarveoutCutoff(const int carveout_cutoff_in) {
  validateDistanceCutoff(&carveout_cutoff, carveout_cutoff_in, "setCarveoutCutoff");
}

//-------------------------------------------------------------------------------------------------
void MMGBSAControls::validateTemperature() {
  if (temperature < 0.0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Negative temperatures are forbidden.", "MMGBSAControls", "validateTemperature");
    case ExceptionResponse::WARN:
      rtWarn("Negative temperatures are forbidden.  The default temperature of " +
             realToString(default_mmgbsa_temperature) + " will be restored.", "MMGBSAControls",
             "validateTemperature");
      temperature = default_mmgbsa_temperature;
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MMGBSAControls::validateDistanceCutoff(double *cutoff, const double cutoff_in,
                                            const char* caller) {
  if (cutoff_in < 0.0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A cutoff distance of " + realToString(cutoff_in, 7, 4, NumberFormat::STANDARD_REAL) +
            " is unrealistic.", "MMGBSAControls", caller);
    case ExceptionResponse::WARN:
      rtWarn("A cutoff distance of " + realToString(cutoff_in, 7, 4, NumberFormat::STANDARD_REAL) +
             " is unrealistic.  The current value of " +
             realToString(*cutoff, 7, 4, NumberFormat::STANDARD_REAL) + " will be maintained.",
             "MMGBSAControls", caller);
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  else {
    *cutoff = cutoff_in;
  }
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator mmgbsaInput(const TextFile &tf, int *start_line, bool *found,
                             const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("mmgbsa", CaseSensitivity::AUTOMATIC, policy, "Collects instructions for "
                         "MM-GBSA resampling in STORMM.");
  const std::string carve_help("A binding site carveout helps define the problem and focus the "
                               "molecular mechanics calculations on the region(s) of the "
                               "receptor that matter.  The carveout may be subdivided into two "
                               "functional parts: the 'active' atoms, which are mobile, and the "
                               "'halo' atoms, which are explicity represented but immobile.");
  const std::vector<std::string> carve_keys_help = {
    "The absolute distance from any of the relevant ligands' atoms by which an atom of the "
    "receptor may be deemed part of the active binding site",
    "A mask of receptor atoms which, irrespective of their distance from ligand atoms, will be "
    "included in the active binding site",
    "The halo region of the binding site, the non-active region where atoms are explicitly "
    "represented but held frozen, may be defined by a number of bonded connections from any "
    "active atoms.  If specified, this number will be used in conjunction with any distance "
    "criterion.",
    "The halo region of the binding site may be defined by an additional layer with this "
    "thickness.  Atoms of the receptor at least the cutoff from any relevant ligands' atoms but "
    "less than this distance will be included.",
    "Specify the file name (or base name) of one of the ligands to make that the reference for "
    "selecting atoms of the receptor.  By default, all ligands will take part in the selection.",
    "The maximum number of additional atoms in a rotatable group that will also be included in "
    "the carveout if at least one atom of the rotatable group is included in the distance-based "
    "definition of the carveout" };
  t_nml.addKeyword("ligand_max_poses", NamelistType::INTEGER,
                   std::to_string(default_ligand_max_poses));
  t_nml.addKeyword("best_poses", NamelistType::INTEGER, std::to_string(default_best_poses));
  t_nml.addKeyword("weights", NamelistType::STRING, std::string(default_weighting_scheme));
  t_nml.addKeyword("carveout", { "cutoff", "halo", "mask", "reference", "add_rotg" },
                   { NamelistType::REAL, NamelistType::REAL, NamelistType::STRING,
                     NamelistType::STRING, NamelistType::INTEGER },
                   { std::to_string(default_carveout_cutoff),
                     std::to_string(default_carveout_halo), std::string(default_carveout_mask),
                     std::string("all"), std::to_string(default_rotating_group_limit) },
                   DefaultIsObligatory::NO, InputRepeats::NO, carve_help, carve_keys_help,
                   { KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL });
  t_nml.addKeyword("whole_residues", NamelistType::BOOLEAN);
  t_nml.addKeyword("rcpt_from_cmpx", NamelistType::BOOLEAN);
  t_nml.addKeyword("temperature", NamelistType::REAL, std::to_string(default_mmgbsa_temperature));
  t_nml.addKeyword("depth", NamelistType::STRING, std::string(default_energy_report_depth));
  t_nml.addKeyword("nrg_decomp", NamelistType::BOOLEAN);
  t_nml.addKeyword("final_structures", NamelistType::BOOLEAN);
  t_nml.addKeyword("structure_fmt", NamelistType::STRING, std::string(default_structure_format));
  t_nml.addKeyword("structure_base_name", NamelistType::STRING,
                   std::string(default_structure_base_name));
  t_nml.addHelp("ligand_max_poses", "Enter the maximum number of poses that should be considered "
                "when making an estimate of the overall binding energy of any ligand.  If this "
                "number is zero or exceeds the number of poses provided for an individual ligand, "
                "the entirety of the poses will be made available for the estimate.  If this "
                "number is less than the number of provided poses, the best poses obtained after "
                "energy minimization will be incorporated into the average.");
  t_nml.addHelp("best_poses", "The number of best poses for which to report energies and, if also "
                "requested, final structures.  This input is required if the \"depth\" keyword "
                "carries a specification of BEST.");
  t_nml.addHelp("weights", "Indicate how multiple poses should be weighted when taking their "
                "average.  The averaging applies to the best scoring poses, as specified by the "
                "\"ligand_max_poses\" keyword.  Options include FLAT / EVEN / EQUAL / NONE, "
                "meaning that all poses of a compound get equal weight when taking an average for "
                "the final predicted binding energy, or BOLTZMANN, indicating that a Boltzmann "
                "weight should be assigned tp each pose result when taking the average.");
  t_nml.addHelp("whole_residues", "Indicate that complete residues should be included in the "
                "carveout if any of their atoms are selected by one of the other rules governing "
                "the active (mobile) atom subset.");
  t_nml.addHelp("rcpt_from_cmpx", "Indicate that the energy-minimized structure of the receptor "
                "should be taken from the energy-minimized structure of the complex, not from an "
                "independent energy minimization of the unbound receptor.  This may eliminate "
                "many extraneous degrees of freedom in the receptor, leading to a better "
                "converged calculation.");
  t_nml.addHelp("temperature", "The temperature, in Kelvin, at which to compute "
                "Boltzmann-weighted averages and, if applicable, dissociation constants.");
  t_nml.addHelp("depth", "Indicate the detail in which energies from individual poses should be "
                "reported.  Options include COMPLETE (the energies of all poses will be reported, "
                "in addition to averages if there are any ligands with more than one pose), "
                "AVERAGES (only the averages will be reported, calculated according to whatever "
                "scheme indicated by the keyword \"weights\"), or BEST (in addition to a list of "
                "all ligands' averages, the best binding results up to the limit given by the "
                "keyword \"best_poses\" will be listed).");
  t_nml.addHelp("nrg_decomp", "Trigger a printout of molecular mechanics energy components in the "
                "final analysis.  Specific components can be named in the accompanying &energy "
                "namelist, and if energy components are specified there then this flag will also "
                "be activated.");
  t_nml.addHelp("final_structures", "Request that the final structures of poses for which "
                "energies are reported also be printed.  These structures will be printed with "
                "the ligand and binding site carveout atoms in the format requested through the "
                "\"structure_fmt\" keyword.");
  t_nml.addHelp("structure_fmt", "The format in which to print final structures of the bound "
                "ligands.  Options include PDB and AMBER (topology and input coordinates will be "
                "printed).");
  t_nml.addHelp("structure_base_name", "Base name of output structures.  The names of each "
                "printed pose will have the form "
                "\"[base][ligand name string, (optional)][index].[ext]\", where the ligand name "
                "string is the unique portion taken from its topology file name and the extension "
                "is detemrined by the chosen output format (structure_fmt)");
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  return t_nml;
}

} // namespace mmgbsa
