#include <cstdlib>
#include <string>
#include <vector>
#include "copyright.h"
#include "FileManagement/file_listing.h"
#include "Namelists/namelist_emulator.h"
#include "Namelists/namelist_element.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Reporting/error_format.h"
#include "Restraints/restraint_enumerators.h"
#include "Trajectory/trajectory_enumerators.h"
#include "user_settings.h"

namespace stormm {
namespace namelist {

using constants::CaseSensitivity;
using constants::translateExceptionResponse;
using diskutil::DrivePathType;
using diskutil::getDrivePathType;
using errors::rtErr;
using errors::rtWarn;
using parse::NumberFormat;
using parse::TextOrigin;
using parse::verifyNumberFormat;
using parse::WrapTextSearch;
using restraints::RestraintEnsemble;
using restraints::translateRestraintEnsemble;
using trajectory::getEnumerationName;
using trajectory::translateCoordinateFileKind;
  
//-------------------------------------------------------------------------------------------------
UserSettings::UserSettings(const CommandLineParser &clip,
                           const std::vector<std::string> &sys_reqs) :
    policy{ExceptionResponse::DIE}, print_policy{default_file_writing_directive},
    has_files_nml{false}, has_debug_nml{false}, has_minimize_nml{false}, has_solvent_nml{false},
    has_random_nml{false}, has_precision_nml{false}, has_conformer_nml{false},
    has_receptor_nml{false}, has_pppm_nml{false}, has_dynamics_nml{false}, has_remd_nml{false},
    has_analysis_nml{false}, has_ffmorph_nml{false}, has_emulator_nml{false},
    has_report_nml{false}, restraint_nml_count{0}, analysis_nml_count{0},
    input_file{std::string(default_conformer_input_file)},
    file_io_input{}, debug_io_input{}, line_min_input{}, solvent_input{}, prng_input{},
    conf_input{}, receptor_input{}, pppm_input{}, dyna_input{}, remd_input{},
    ffmod_input{}, emul_input{}, diagnostic_input{}, rstr_inputs{}, analysis_inputs{}
{
  // Local variables to store command line arguments
  int cval_igseed = 0;
  std::string cval_report_file, cval_traj_file_name, cval_input_transcript_file;
  std::vector<std::string> cval_topology_file_names;
  std::vector<std::string> cval_coordinate_file_names;
  
  // Detect command line arguments, and note that their presence overrides similar directives
  // in the input deck.
  bool cli_inpfile          = false;
  bool cli_igseed           = false;
  bool cli_report           = false;
  bool cli_input_transcript = false;
  bool cli_trajname         = false;
  CoordinateFileKind c_kind = default_filecon_inpcrd_type;
  CoordinateFileKind x_kind = default_filecon_outcrd_type;
  CoordinateFileKind r_kind = default_filecon_chkcrd_type;
  const NamelistEmulator *t_nml = clip.getNamelistPointer();
  const InputStatus user_spec = InputStatus::USER_SPECIFIED;
  if (t_nml->hasKeyword("-i") && t_nml->getKeywordStatus("-i") == user_spec) {
    input_file = t_nml->getStringValue("-i");
    cli_inpfile = true;
  }
  if (t_nml->hasKeyword("-ig_seed") &&
      t_nml->getKeywordStatus("-ig_seed") == user_spec) {
    cval_igseed = t_nml->getIntValue("-ig_seed");
    cli_igseed = true;
  }
  if (t_nml->hasKeyword("-p") && t_nml->getKeywordStatus("-p") == user_spec) {
    cval_topology_file_names = t_nml->getAllStringValues("-p");
  }
  if (t_nml->hasKeyword("-c") && t_nml->getKeywordStatus("-c") == user_spec) {
    cval_coordinate_file_names = t_nml->getAllStringValues("-c");
  }
  if (t_nml->hasKeyword("-o") && t_nml->getKeywordStatus("-o") == user_spec) {
    cval_report_file = t_nml->getStringValue("-o");
    cli_report = true;
  }
  if (t_nml->hasKeyword("-t") && t_nml->getKeywordStatus("-t") == user_spec) {
    cval_input_transcript_file = t_nml->getStringValue("-t");
    cli_input_transcript = true;
  }
  if (t_nml->hasKeyword("-x") && t_nml->getKeywordStatus("-x") == user_spec) {
    cval_traj_file_name = t_nml->getStringValue("-x");
    cli_trajname = true;
  }
  if (t_nml->hasKeyword("-c_kind") && t_nml->getKeywordStatus("-c_kind") == user_spec) {
    c_kind = translateCoordinateFileKind(t_nml->getStringValue("-c_kind"));
  }
  if (t_nml->hasKeyword("-x_kind") && t_nml->getKeywordStatus("-x_kind") == user_spec) {
    x_kind = translateCoordinateFileKind(t_nml->getStringValue("-x_kind"));
  }
  if (t_nml->hasKeyword("-r_kind") && t_nml->getKeywordStatus("-r_kind") == user_spec) {
    r_kind = translateCoordinateFileKind(t_nml->getStringValue("-r_kind"));
  }
  if (t_nml->hasKeyword("-except") && t_nml->getKeywordStatus("-except") == user_spec) {
    policy = translateExceptionResponse(t_nml->getStringValue("-except"));
  }
  if (t_nml->hasKeyword("-O") && t_nml->getBoolValue("-O")) {
    print_policy = PrintSituation::OVERWRITE;
  }
  
  // Process the input file.  Take only the first instance of each namelist, as found by searching
  // from the beginning.
  if (getDrivePathType(input_file) != DrivePathType::FILE) {
    const std::string descriptor = (cli_inpfile) ? std::string("user specified") :
                                                   std::string("default");
    rtErr("The " + descriptor + " input file " + input_file + " was not found or could not be "
          "read.", "UserSettings");
  }
  TextFile inp_tf(input_file, TextOrigin::DISK, "Input deck for STORMM executable",
                  "UserSettings");

  std::vector<std::string> alternatives = {
    "coordinate_input_format",      getEnumerationName(c_kind),
    "coordinate_output_format",     getEnumerationName(x_kind),
    "coordinate_checkpoint_format", getEnumerationName(r_kind)
  };
  int start_line = 0;
  file_io_input = FilesControls(inp_tf, &start_line, &has_files_nml, policy, WrapTextSearch::NO,
                                alternatives, sys_reqs);
  start_line = 0;
  debug_io_input = DebugControls(inp_tf, &start_line, &has_debug_nml, policy);
  start_line = 0;
  line_min_input = MinimizeControls(inp_tf, &start_line, &has_minimize_nml, policy);
  start_line = 0;
  solvent_input = SolventControls(inp_tf, &start_line, &has_solvent_nml, policy);
  start_line = 0;
  prng_input = RandomControls(inp_tf, &start_line, &has_random_nml, policy);
  start_line = 0;
  prec_input = PrecisionControls(inp_tf, &start_line, &has_precision_nml, policy);
  start_line = 0;
  conf_input = ConformerControls(inp_tf, &start_line, &has_conformer_nml, policy);
  start_line = 0;
  receptor_input = ReceptorControls(inp_tf, &start_line, &has_receptor_nml, policy);
  start_line = 0;
  pppm_input = PPPMControls(inp_tf, &start_line, &has_pppm_nml, policy);
  start_line = 0;
  dyna_input = DynamicsControls(inp_tf, &start_line, &has_dynamics_nml, policy);
  start_line = 0;
  remd_input = RemdControls(inp_tf, &start_line, &has_remd_nml, policy);
  start_line = 0;
  ffmod_input = FFMorphControls(inp_tf, &start_line, &has_ffmorph_nml, policy);
  start_line = 0;
  emul_input = EmulatorControls(inp_tf, &start_line, &has_emulator_nml, policy);
  start_line = 0;
  diagnostic_input = ReportControls(inp_tf, &start_line, &has_report_nml, policy);
  start_line = 0;
  while (start_line < inp_tf.getLineCount()) {
    bool restraint_nml_found = false;
    RestraintControls tmp_rstr_input = RestraintControls(inp_tf, &start_line, &restraint_nml_found,
                                                         policy);
    if (restraint_nml_found) {
      restraint_nml_count += 1;
      rstr_inputs.push_back(tmp_rstr_input);
    }
  }
  start_line = 0;
  while (start_line < inp_tf.getLineCount()) {
    bool analysis_nml_found = false;
    AnalysisControls tmp_anls_input = AnalysisControls(inp_tf, &start_line, &analysis_nml_found,
                                                       policy);
    if (analysis_nml_found) {
      analysis_nml_count += 1;
      analysis_inputs.push_back(tmp_anls_input);
    }
  }


  // If multiple &analysis control blocks are present, make a note of this within the object
  // storing user input from each such block.  The AnalysisControls objects are then able to
  // modify the analyses they help create so that output reflects the entirety of user input.
  if (analysis_nml_count > 1) {
    for (int i = 0; i < analysis_nml_count; i++) {
      analysis_inputs[i].multipleAnalysesFound(true);
    }
    std::string identifier_root(default_analysis_identifier);
    int next_unique_idno = 0;
    for (int i = 0; i < analysis_nml_count; i++) {
      const NamelistEmulator& i_nml = analysis_inputs[i].getTranscript();
      if (i_nml.getKeywordStatus("tag") == InputStatus::MISSING) {
        bool tagged = false;
        do {
          const std::string next_unique_identifier = identifier_root +
                                                     std::to_string(next_unique_idno);
          bool matched = false;
          int j = 0;
          while (j < analysis_nml_count && matched == false) {
            const NamelistEmulator& j_nml = analysis_inputs[j].getTranscript();
            if (j_nml.getKeywordStatus("tag") != InputStatus::MISSING &&
                analysis_inputs[j].getIdentifier() == next_unique_identifier) {
              matched = true;
            }
            j++;
          }
          if (matched == false) {
            analysis_inputs[i].setIdentifier(next_unique_identifier);
            tagged = true;
          }
          next_unique_idno++;
        } while (tagged == false);
      }
    }
  }
  
  // Superimpose, or contribute, command line directives
  if (cli_igseed) {
    prng_input.setRandomSeed(cval_igseed);    
  }
  if (cli_report) {
    file_io_input.setReportFileName(cval_report_file);
  }
  if (cli_trajname) {
    file_io_input.setGeneralTrajectoryFileName(cval_traj_file_name);
  }
  if (cli_input_transcript) {
    file_io_input.setInputTranscriptFileName(cval_input_transcript_file);
  }
  if (cval_topology_file_names.size() > 0LLU) {
    for (size_t i = 0; i < cval_topology_file_names.size(); i++) {
      file_io_input.addFreeTopologyName(cval_topology_file_names[i]);
    }
  }
  if (cval_coordinate_file_names.size() > 0LLU) {
    for (size_t i = 0; i < cval_coordinate_file_names.size(); i++) {
      file_io_input.addFreeCoordinateName(cval_coordinate_file_names[i]);
    }
  }

  // Impose checks on user input.  If positional restraints are imposed, the system should not be
  // recentered, as this would generate new forces of its own.
  if (rstr_inputs.size() > 0 && has_dynamics_nml) {
    int nposn_restraints = 0;
    for (size_t i = 0; i < rstr_inputs.size(); i++) {
      if (rstr_inputs[i].getOrder() == 0) {
        const NamelistEmulator &t_nml = rstr_inputs[i].getTranscript();
        try {
          switch (translateRestraintEnsemble(t_nml.getStringValue("ensemble"))) {
          case RestraintEnsemble::SPECIFIC_ATOMS:
          case RestraintEnsemble::PRESERVE_POSITIONS:
            nposn_restraints++;
            break;
          case RestraintEnsemble::PREVENT_HBONDS:
          case RestraintEnsemble::PRESERVE_HEAVY_DIHEDRALS:
          case RestraintEnsemble::PRESERVE_DISTANCES:
            break;
          }
        }
        catch (std::runtime_error) {
          switch (policy) {
          case ExceptionResponse::DIE:
            rtErr("An unrecognized restraint ensemble " + t_nml.getStringValue("ensemble") +
                  " was found in a &restraint namelist control block during late-stage input "
                  "checks.", "UserSettings");
          case ExceptionResponse::WARN:
            rtWarn("An unrecognized restraint ensemble " + t_nml.getStringValue("ensemble") +
                   " was found in a &restraint namelist control block during late-stage input "
                   "checks.", "UserSettings");
            break;
          case ExceptionResponse::SILENT:
            break;
          }
        }
      }
      else if (rstr_inputs[i].getOrder() == 1) {
        nposn_restraints++;
      }
    }
    if (nposn_restraints > 0) {
      if (dyna_input.getCenterOfMassMotionPurgeFrequency() > 0) {
        switch (dyna_input.getTranscript().getKeywordStatus("nscm")) {
        case InputStatus::USER_SPECIFIED:
          switch (policy) {
          case ExceptionResponse::DIE:
            rtErr("Positional restraints were specified for a simulation with center of mass "
                  "zeroing.  This will introduce unnatural forces and conflicts every time the "
                  "center of mass is returned to the origin.  Rely on the positional restraints "
                  "to keep the system roughly in the space that is expected.", "UserSettings");
          case ExceptionResponse::WARN:
            rtWarn("Positional restraints were specified for a simulation with center of mass "
                   "zeroing.  This will introduce unnatural forces and conflicts every time the "
                   "center of mass is returned to the origin.  Center of mass repositioning and "
                   "momentum removal will be ignored.", "UserSettings");
            break;
          case ExceptionResponse::SILENT:
            break;
          }
        case InputStatus::DEFAULT:
          switch (policy) {
          case ExceptionResponse::DIE:
          case ExceptionResponse::WARN:
            rtWarn("Center of mass drift and momentum removal will be disabled for a simulation "
                   "making use of positional restraints.", "UserSettings");
          case ExceptionResponse::SILENT:
            break;
          }
        case InputStatus::MISSING:
          break;
        }
        dyna_input.setCenterOfMassMotionPurgeFrequency(0);
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
ExceptionResponse UserSettings::getExceptionBehavior() const {
  return policy;
}
  
//-------------------------------------------------------------------------------------------------
const std::string& UserSettings::getInputFileName() const {
  return input_file;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getFilesPresence() const {
  return has_files_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getDebugPresence() const {
  return has_debug_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getMinimizePresence() const {
  return has_minimize_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getSolventPresence() const {
  return has_solvent_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getRandomPresence() const {
  return has_random_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getPrecisionPresence() const {
  return has_precision_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getConformerPresence() const {
  return has_conformer_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getPPPMPresence() const {
  return has_pppm_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getDynamicsPresence() const {
  return has_dynamics_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getRemdPresence() const {
  return has_remd_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getAnalysisPresence() const {
  return has_analysis_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getFFMorphPresence() const {
  return has_ffmorph_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getEmulatorPresence() const {
  return has_emulator_nml;
}

//-------------------------------------------------------------------------------------------------
bool UserSettings::getReportPresence() const {
  return has_report_nml;
}

//-------------------------------------------------------------------------------------------------
const FilesControls& UserSettings::getFilesNamelistInfo() const {
  return file_io_input;
}

//-------------------------------------------------------------------------------------------------
const DebugControls& UserSettings::getDebugNamelistInfo() const {
  return debug_io_input;
}

//-------------------------------------------------------------------------------------------------
const MinimizeControls& UserSettings::getMinimizeNamelistInfo() const {
  return line_min_input;
}

//-------------------------------------------------------------------------------------------------
const SolventControls& UserSettings::getSolventNamelistInfo() const {
  return solvent_input;
}

//-------------------------------------------------------------------------------------------------
const RandomControls& UserSettings::getRandomNamelistInfo() const {
  return prng_input;
}

//-------------------------------------------------------------------------------------------------
const PrecisionControls& UserSettings::getPrecisionNamelistInfo() const {
  return prec_input;
}

//-------------------------------------------------------------------------------------------------
const ConformerControls& UserSettings::getConformerNamelistInfo() const {
  return conf_input;
}

//-------------------------------------------------------------------------------------------------
const ReceptorControls& UserSettings::getReceptorNamelistInfo() const {
  return receptor_input;
}

//-------------------------------------------------------------------------------------------------
const PPPMControls& UserSettings::getPPPMNamelistInfo() const {
  return pppm_input;
}

//-------------------------------------------------------------------------------------------------
const DynamicsControls& UserSettings::getDynamicsNamelistInfo() const {
  return dyna_input;
}

//-------------------------------------------------------------------------------------------------
const RemdControls& UserSettings::getRemdNamelistInfo() const {
  return remd_input;
}

//-------------------------------------------------------------------------------------------------
const FFMorphControls& UserSettings::getFFMorphNamelistInfo() const {
  return ffmod_input;
}

//-------------------------------------------------------------------------------------------------
const EmulatorControls& UserSettings::getEmulatorNamelistInfo() const {
  return emul_input;
}

//-------------------------------------------------------------------------------------------------
const ReportControls& UserSettings::getReportNamelistInfo() const {
  return diagnostic_input;
}

//-------------------------------------------------------------------------------------------------
int UserSettings::getRestraintNamelistCount() const {
  return restraint_nml_count;
}
  
//-------------------------------------------------------------------------------------------------
const std::vector<RestraintControls>& UserSettings::getRestraintNamelistInfo() const {
  return rstr_inputs;
}

//-------------------------------------------------------------------------------------------------
const RestraintControls& UserSettings::getRestraintNamelistInfo(const int index) const {
  if (index < 0 || index >= restraint_nml_count) {
    rtErr("The input contained " + std::to_string(restraint_nml_count) + " &restraint namelists.  "
          "Index " + std::to_string(index) + " is invalid.", "UserSettings",
          "getRestraintNamelistInfo");
  }
  return rstr_inputs[index];
}

//-------------------------------------------------------------------------------------------------
int UserSettings::getAnalysisNamelistCount() const {
  return analysis_nml_count;
}
  
//-------------------------------------------------------------------------------------------------
const std::vector<AnalysisControls>& UserSettings::getAnalysisNamelistInfo() const {
  return analysis_inputs;
}

//-------------------------------------------------------------------------------------------------
const AnalysisControls& UserSettings::getAnalysisNamelistInfo(const int index) const {
  if (index < 0 || index >= analysis_nml_count) {
    rtErr("The input contained " + std::to_string(analysis_nml_count) + " &analysis namelists.  "
          "Index " + std::to_string(index) + " is invalid.", "UserSettings",
          "getAnalysisNamelistInfo");
  }
  return analysis_inputs[index];
}

//-------------------------------------------------------------------------------------------------
PrintSituation UserSettings::getPrintingPolicy() const {
  return print_policy;
}

} // namespace namelist
} // namespace stormm
