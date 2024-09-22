#include <cstdio>
#include <cstdlib>
#include "copyright.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/textfile.h"
#include "Potential/pme_util.h"
#include "command_line_parser.h"
#include "input.h"

namespace stormm {
namespace namelist {

using energy::default_pme_cutoff;
using parse::findStringInVector;
using parse::TextOrigin;
using parse::TextFile;
using parse::vectorOfStrings;

//-------------------------------------------------------------------------------------------------
CommandLineParser::CommandLineParser(const std::string &program_name_in,
                                     const std::string &program_description,
                                     const std::vector<std::string> &noted_imports_in,
                                     const ExceptionResponse policy_in) :
    policy{policy_in},
    arg_count{0},
    program_name{program_name_in},
    executable{},
    cli_nml{program_name, CaseSensitivity::YES, policy_in, program_description, true},
    help_on_no_args{true}, exit_on_help{true}, lead_parser{true},
    command_line_text{},
    coordinations{},
    excluded_keys{},
    noted_imports{noted_imports_in}
{
  cli_nml.addKeyword("--help", NamelistType::BOOLEAN);
  cli_nml.addHelp("--help", "List command line arguments with descriptions.");
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::activateHelpOnNoArgs() {
  help_on_no_args = true;
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::suppressHelpOnNoArgs() {
  help_on_no_args = false;
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::activateExitOnHelp() {
  exit_on_help = true;
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::preventExitOnHelp() {
  exit_on_help = false;
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::addStandardAmberInputs(const std::vector<std::string> &cli_keys) {
  const size_t n_keys = cli_keys.size();
  for (size_t i = 0; i < n_keys; i++) {
    if (cli_keys[i] == "-i") {
      cli_nml.addKeyword("-i", NamelistType::STRING, std::string(default_stormm_input_deck));
      cli_nml.addHelp("-i", "The primary input file, equivalent to Amber's mdin.");
    }
    else if (cli_keys[i] == "-O") {
      cli_nml.addKeyword("-O", NamelistType::STRING, std::string(""));
      cli_nml.addHelp("-O", "Flag to activate overwriting of existing output files.");
    }
    else if (cli_keys[i] == "-p") {
      cli_nml.addKeyword("-p", NamelistType::STRING, std::string(default_stormm_input_topology));
      cli_nml.addHelp("-p", "A primary topology file to use in generating a synthesis of "
                      "systems.  This keyword may be repeated to include multiple systems.");
    }
    else if (cli_keys[i] == "-c") {
      cli_nml.addKeyword("-c", NamelistType::STRING,
                         std::string(default_stormm_input_coordinates));
      cli_nml.addHelp("-c", "A primary input coordinates file to use in generating a synthesis "
                      "of systems.  This keyword may be repeated.  Similar file names and, after "
                      "that, similar numbers of atoms with reasonable bond and angle energies "
                      "will be used to correlate one set of input coordinates with each topology "
                      "if systems are defined on the command line by this route.  For explicit "
                      "control of which topology will describe which coordinate set, use a "
                      "&files namelist in the input deck (-i).");
    }
    else if (cli_keys[i] == "-o") {
      cli_nml.addKeyword("-o", NamelistType::STRING,
                         std::string(default_stormm_report_file));
      cli_nml.addHelp("-o", "The main report file to write at the conclusion of calculations.  "
                      "Similar in nature to Amber's mdout.");
    }
    else if (cli_keys[i] == "-x") {
      cli_nml.addKeyword("-x", NamelistType::STRING,
                         std::string(default_stormm_output_trajectory));
      cli_nml.addHelp("-x", "The base name for output trajectories.  For explicit control of the "
                      "trajectory files that will be linked to each system in the synthesis, "
                      "use the &files namelist.");
    }
    else if (cli_keys[i] == "-r") {
      cli_nml.addKeyword("-r", NamelistType::STRING,
                         std::string(default_stormm_output_checkpoint));
      cli_nml.addHelp("-r", "The base name for output checkpoint files.  For explicit control "
                      "of the checkpoint files that will be linked to each system in the "
                      "synthesis, use the &files namelist.");
    }
    else if (cli_keys[i] == "-ig_seed") {
      cli_nml.addKeyword("-ig_seed", NamelistType::INTEGER, std::to_string(629295034));
      cli_nml.addHelp("-ig_seed", "Seed for the random number gnenrator which will create "
                      "perturbations in each replica's coordinates.");
    }
    else {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Invalid AMBER command line input keyword \"" + cli_keys[i] + "\".",
              "CommandLineParser", "addStandardAmberInput");
      case ExceptionResponse::WARN:
        rtWarn("Invalid AMBER command line input keyword \"" + cli_keys[i] + "\".  This input "
               "will be ignored.", "CommandLineParser", "addStandardAmberInput");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::addStandardAmberInputs(const char* key_a, const char* key_b,
                                               const char* key_c, const char* key_d) {
  addStandardAmberInputs(vectorOfStrings(key_a, key_b, key_c, key_d));
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::addStandardBenchmarkingInputs(const std::vector<std::string> &cli_keys) {
  const size_t n_keys = cli_keys.size();
  for (size_t i = 0; i < n_keys; i++) {
    if (cli_keys[i] == "-iter") {
      cli_nml.addKeyword("-iter", NamelistType::INTEGER, std::to_string(100));
      cli_nml.addHelp("-iter", "The number of iterations within each trial segment.  Choose a "
                      "value that will make the time cost of all iterations much greater than the "
                      "cost of synchronization between the CPU host and GPU device.");
    }
    else if (cli_keys[i] == "-replicas") {
      cli_nml.addKeyword("-replicas", NamelistType::INTEGER, std::to_string(1));
      cli_nml.addHelp("-replicas", "The number of replicas of a given system that will be used to "
                      "create the testing setup.");
    }
    else if (cli_keys[i] == "-trials") {
      cli_nml.addKeyword("-trials", NamelistType::INTEGER, std::to_string(4));
      cli_nml.addHelp("-trials", "The number of separate trials to conduct.  Timings will be "
                      "reported as an average of all trials, each of which may involve multiple "
                      "iterations of the calculation.");
    }
    else if (cli_keys[i] == "-cutoff") {
      cli_nml.addKeyword("-cutoff", NamelistType::REAL, std::to_string(default_pme_cutoff));
      cli_nml.addHelp("-cutoff", "The cutoff, in units of Angstroms, that will be applied to "
                      "particle-particle interactions.  This can specify the cutoff for both "
                      "electrostatic and van-der Waals interactions.");
    }
    else if (cli_keys[i] == "-elec_cutoff") {
      cli_nml.addKeyword("-elec_cutoff", NamelistType::REAL, std::string(""));
      cli_nml.addHelp("-elec_cutoff", "The cutoff, in units of Angstroms, that will be applied to "
                      "electrostatic particle-particle interactions.  This will override a value "
                      "specified with a more general keyword (e.g. -cutoff).");
    }
    else if (cli_keys[i] == "-vdw_cutoff") {
      cli_nml.addKeyword("-vdw_cutoff", NamelistType::REAL, std::string(""));
      cli_nml.addHelp("-vdw_cutoff", "The cutoff, in units of Angstroms, that will be applied to "
                      "van-der Waals particle-particle interactions.  This will override a value "
                      "specified with a more general keyword (e.g. -cutoff).");
    }
    else if (cli_keys[i] == "-pad") {
      cli_nml.addKeyword("-pad", NamelistType::REAL, std::to_string(0.05));
      cli_nml.addHelp("-pad", "The minimum padding, in units of Angstroms, added to the width of "
                      "each neighbor list cell when subdividing the simulation unit cell.");
    }
  }
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::addStandardBenchmarkingInputs(const char* key_a, const char* key_b,
                                                      const char* key_c, const char* key_d) {
  addStandardBenchmarkingInputs(vectorOfStrings(key_a, key_b, key_c, key_d));
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator* CommandLineParser::getNamelistPointer() {
  return &cli_nml;
}

//-------------------------------------------------------------------------------------------------
const NamelistEmulator* CommandLineParser::getNamelistPointer() const {
  return &cli_nml;
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::parseUserInput(const int argc, const char* argv[]) {
  if (argc == 0) {

    // There should always be at least one command line argument (the program name), but take an
    // empty list as an indication that the help message should be printed.
    cli_nml.printHelp();
    return;
  }
  executable = std::string(argv[0]);
  command_line_text = vectorOfStrings(&argv[1], argc - 1);
  std::string tmp_cli("&" + program_name + " ");
  const int nexcl = excluded_keys.size();
  for (int i = 0; i < argc - 1; i++) {
    
    // Check for excluded keywords and their types
    if (excluded_keys.find(command_line_text[i]) != excluded_keys.end()) {
      switch (excluded_keys.at(command_line_text[i])) {
      case NamelistType::BOOLEAN:
        break;
      case NamelistType::INTEGER:
      case NamelistType::REAL:
      case NamelistType::STRING:

        // Skip the next keyword as well, if it is also not in the object's namelist.  If it is,
        // then there is an error either because the user is taking a reserved word as an input
        // value, or because the user should have included a value after the keyword but forgot.
        if (i <  argc - 2 && cli_nml.hasKeyword(command_line_text[i]) == false) {
          i++;
        }
        break;
      case NamelistType::STRUCT:
        if (i < argc - 2 && command_line_text[i + 1].size() == 1 &&
            (command_line_text[i + 1][0] == '[' || command_line_text[i + 1][0] == '(' ||
             command_line_text[i + 1][0] == '{')) {
          i++;
          while (i < argc - 1 && (command_line_text[i].size() > 1 ||
                                  (command_line_text[i][0] != ']' &&
                                   command_line_text[i][0] != ')' &&
                                   command_line_text[i][0] != '}'))) {
            i++;
          }
        }
        break;
      }
    }
    else {
      tmp_cli += command_line_text[i] + " ";
    }
  }
  tmp_cli += "&end";
  const TextFile tf(tmp_cli, TextOrigin::RAM);
  bool found = false;
  readNamelist(tf, &cli_nml, 0, WrapTextSearch::NO, tf.getLineCount(), &found);

  // Print help messages if the developer indicates and there are no arguments, or if the user
  // has specified "--help" as one of the command line keywords.
  if ((help_on_no_args && argc == 1) || cli_nml.getBoolValue("--help")) {
    if (lead_parser) {
      if (coordinations.size() == 0) {
        cli_nml.printHelp();
      }
      else {
        NamelistEmulator tmp_nml = cli_nml;
        const int n_imports = noted_imports.size();
        for (size_t i = 0; i < coordinations.size(); i++) {
          const NamelistEmulator *oth_nml = coordinations[i]->getNamelistPointer();
          for (int j = 0; j < oth_nml->getKeywordCount(); j++) {
            const std::string& oth_key = oth_nml->getKeyword(j);
            if (findStringInVector(noted_imports, oth_key) < n_imports &&
                tmp_nml.hasKeyword(oth_key) == false) {
              tmp_nml.addKeyword(oth_nml, oth_key);
            }
          }
        }
        tmp_nml.printHelp();
      }
      if (exit_on_help) {
        exit(0);
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::coordinateWithPartner(CommandLineParser *other) {
  lead_parser = false;
  
  // If this object is already in the list of coordinate objects, skip.
  const int n_coordinated = coordinations.size();
  bool already_coordinated = false;
  for (int i = 0; i < n_coordinated; i++) {
    already_coordinated = (already_coordinated || (coordinations[i] == other));
  }
  if (already_coordinated) {
    return;
  }
  
  // Check all keywords of the partner, taking note of their types.
  const int partner_kw_count = other->cli_nml.getKeywordCount();
  int user_spec_keys = 0;
  for (int i = 0; i < partner_kw_count; i++) {
    const std::string& partner_key = other->cli_nml.getKeyword(i);
    const NamelistType partner_key_kind = other->cli_nml.getKeywordKind(partner_key);
    if (cli_nml.hasKeyword(partner_key, partner_key_kind) == false) {

      // Add to the object's own excluded keys list
      if (excluded_keys.find(partner_key) == excluded_keys.end()) {
        excluded_keys[partner_key] = partner_key_kind;
      }
    }
    
    // Loop over other coordinated objects and add the new partner's keys to their excluded
    // key lists, as appropriate.
    for (int j = 0; j < n_coordinated; j++) {
      if (coordinations[j]->cli_nml.hasKeyword(partner_key) == false) {
        if (coordinations[j]->excluded_keys.find(partner_key) ==
            coordinations[j]->excluded_keys.end()) {
          coordinations[j]->excluded_keys[partner_key] = partner_key_kind;
        }
      }
    }
  }

  // Check all keywords of this object, adding them to the new coordinated partner's excluded
  // keys if they are not also keywords in the partner CommandLineParser.
  const int self_kw_count = cli_nml.getKeywordCount();
  for (int i = 0; i < self_kw_count; i++) {
    const std::string& self_key = cli_nml.getKeyword(i);
    const NamelistType self_key_kind = cli_nml.getKeywordKind(self_key);
    if (other->cli_nml.hasKeyword(self_key, self_key_kind) == false) {
      if (other->excluded_keys.find(self_key) == other->excluded_keys.end()) {
        other->excluded_keys[self_key] = self_key_kind;
      }
    }
  }

  // Check all keywords of other coordinated CommandLineParsers, adding them to the new partner's
  // excluded keys as was just done for this object.
  for (int i = 0; i < n_coordinated; i++) {
    const int icoord_key_count = coordinations[i]->cli_nml.getKeywordCount();
    for (int j = 0; j < icoord_key_count; j++) {
      const std::string& ij_key = coordinations[i]->cli_nml.getKeyword(j);
      const NamelistType ij_key_kind = coordinations[i]->cli_nml.getKeywordKind(ij_key);
      if (other->cli_nml.hasKeyword(ij_key, ij_key_kind) == false) {
        if (other->excluded_keys.find(ij_key) == other->excluded_keys.end()) {
          other->excluded_keys[ij_key] = ij_key_kind;
        }
      }
    }
  }

  // Add the new partner to the list of coordinated objects
  for (int i = 0; i < n_coordinated; i++) {
    coordinations[i]->coordinations.push_back(other);
    other->coordinations.push_back(coordinations[i]);
  }
  coordinations.push_back(other);
  other->coordinations.push_back(this);
}

} // namespace namelist
} // namespace stormm
