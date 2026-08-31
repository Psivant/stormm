#include <cstdio>
#include <cstdlib>
#include "copyright.h"
#include "Namelists/nml_conformer.h"
#include "Namelists/nml_dynamics.h"
#include "Namelists/nml_ffmorph.h"
#include "Namelists/nml_files.h"
#include "Namelists/nml_mesh.h"
#include "Namelists/nml_minimize.h"
#include "Namelists/nml_nice.h"
#include "Namelists/nml_pppm.h"
#include "Namelists/nml_precision.h"
#include "Namelists/nml_random.h"
#include "Namelists/nml_report.h"
#include "Namelists/nml_restraint.h"
#include "Namelists/nml_solvent.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/textfile.h"
#include "Potential/pme_util.h"
#include "Reporting/reporting_enumerators.h"
#include "Reporting/report_table.h"
#include "Reporting/summary_file.h"
#include "command_line_parser.h"
#include "input.h"
#include "namelist_inventory.h"
#include "nml_dynamics.h"
#include "nml_files.h"

namespace stormm {
namespace namelist {

using energy::default_pme_cutoff;
using namelist::conformerInput;
using namelist::dynamicsInput;
using namelist::filesInput;
using namelist::ffmorphInput;
using namelist::meshInput;
using namelist::minimizeInput;
using namelist::niceInput;
using namelist::pppmInput;
using namelist::precisionInput;
using namelist::randomInput;
using namelist::reportInput;
using namelist::restraintInput;
using namelist::solventInput;
using parse::findStringInVector;
using parse::strcmpCased;
using parse::JustifyText;
using parse::TextOrigin;
using parse::TextFile;
using parse::vectorOfStrings;
using review::findFormatWidth;
using review::printProtectedText;
using review::protectText;
using review::ReportTable;
using review::OutputSyntax;
  
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
    control_blocks{},
    coordinations{},
    excluded_keys{},
    noted_imports{noted_imports_in},
    custom_namelists{}
{
  const std::string help_msg("List command line arguments with descriptions.");
  cli_nml.addKeyword("--help", NamelistType::BOOLEAN);
  cli_nml.addHelp("--help", help_msg);
  cli_nml.addKeyword("-help", NamelistType::BOOLEAN);
  cli_nml.addHelp("-help", help_msg);
}

//-------------------------------------------------------------------------------------------------
CommandLineParser::CommandLineParser(const std::string &program_name_in,
                                     const std::string &program_description,
                                     const ExceptionResponse policy_in) :
    CommandLineParser(program_name_in, program_description, {}, policy_in)
{}

//-------------------------------------------------------------------------------------------------
std::string CommandLineParser::getInputAsString(const int width) const {
  std::string result = executable;
  const int n_args = command_line_text.size();
  for (int i = 0; i < n_args; i++) {
    result += " ";
    result += command_line_text[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
bool CommandLineParser::doesProgramExitOnHelp() const {
  return exit_on_help;
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
      cli_nml.addKeyword("-O", NamelistType::BOOLEAN);
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
                      "that, similar numbers of atoms with reasonable bond and angle energies, "
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
    else if (cli_keys[i] == "-temp0") {
      cli_nml.addKeyword("-temp0", NamelistType::REAL,
                         std::to_string(default_simulation_temperature));
      cli_nml.addHelp("-temp0", "The temperature at which to stage systems or calculations");
    }
    else {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Invalid AMBER command line input keyword \"" + cli_keys[i] + "\".",
              "CommandLineParser", "addStandardAmberInputs");
      case ExceptionResponse::WARN:
        rtWarn("Invalid AMBER command line input keyword \"" + cli_keys[i] + "\".  This input "
               "will be ignored.", "CommandLineParser", "addStandardAmberInputs");
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
void CommandLineParser::addStandardAmberInputs() {
  addStandardAmberInputs({ "-i", "-O", "-p", "-c", "-o", "-x", "-r", "-ig_seed", "-temp0" });
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
    else if (cli_keys[i] == "-precision") {
      cli_nml.addKeyword("-precision", NamelistType::STRING, std::string("single"));
      cli_nml.addHelp("-precision", "The precision model in whcih to perform calculations.  Note "
                      "that various coordinates (positions, velocities, forces) will likely "
                      "remain in a fixed-precision representation.");
    }
    else if (cli_keys[i] == "-skip_cpu_check") {
      cli_nml.addKeyword("-skip_cpu_check", NamelistType::BOOLEAN);
      cli_nml.addHelp("-skip_cpu_check", "Skip a CPU_based check on the forces computed by the "
                      "GPU.");
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
    else if (cli_keys[i] == "-eval_nrg") {
      cli_nml.addKeyword("-eval_nrg", NamelistType::BOOLEAN);
      cli_nml.addHelp("-eval_nrg", "Request that the relevant system energy components be "
                      "evaluated");
    }
    else if (cli_keys[i] == "-eval_frc") {
      cli_nml.addKeyword("eval_frc", NamelistType::BOOLEAN);
      cli_nml.addHelp("-eval_frc", "Request that forces on particles due to relevant interactions "
                      "be evaluated");
    }
    else if (cli_keys[i] == "-omit_frc") {
      cli_nml.addKeyword("-omit_frc", NamelistType::BOOLEAN);
      cli_nml.addHelp("-omit_frc", "Request that the forces on particles be omitted from the GPU "
                      "calculation.  This will also omit any CPU check on forces, and compel an "
                      "evaluation of the energy so as to have at least one quantity for the GPU "
                      "to compute.");
    }
    else if (cli_keys[i] == "-fp_bits") {
      cli_nml.addKeyword("-fp_bits", NamelistType::INTEGER, std::to_string(24));
      cli_nml.addHelp("-fp_bits", "The number of fixed precision bits after the point (values in "
                      "kcal/mol-A) with which to accumulate forces on all particles.");
    }
    else if (cli_keys[i] == "-blur") {
      cli_nml.addKeyword("-blur", NamelistType::REAL, std::to_string(0.04));
      cli_nml.addHelp("-blur", "The Gaussian width by which to blur particle positions in each "
                      "replica of the system and each iteration of the force or energy "
                      "calculation, in units of Angstroms.");
    }
    else {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Invalid command line input keyword \"" + cli_keys[i] + "\".",
              "CommandLineParser", "addStandardBenchmarkingInputs");
      case ExceptionResponse::WARN:
        rtWarn("Invalid command line input keyword \"" + cli_keys[i] + "\".  This input "
               "will be ignored.", "CommandLineParser", "addStandardBenchmarkingInputs");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::addStandardBenchmarkingInputs(const char* key_a, const char* key_b,
                                                      const char* key_c, const char* key_d) {
  addStandardBenchmarkingInputs(vectorOfStrings(key_a, key_b, key_c, key_d));
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::addStandardBenchmarkingInputs() {
  addStandardBenchmarkingInputs({ "-iter", "-trials", "-replicas", "-precision", "-skip_cpu_check",
                                  "-cutoff", "-elec_cutoff", "-vdw_cutoff", "-pad", "-eval_nrg",
                                  "-eval_frc", "-omit_frc", "-fp_bits", "-blur" });
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::addStandardApplicationInputs(const std::vector<std::string> &cli_keys) {

  // Shunt some imputs to the Amber loader
  const size_t nkeys = cli_keys.size();
  std::vector<std::string> amb_keys, nonamb_keys;
  for (size_t i = 0; i < nkeys; i++) {
    if (cli_keys[i] == "-i" || cli_keys[i] == "-ig_seed" || cli_keys[i] == "-p" ||
        cli_keys[i] == "-c" || cli_keys[i] == "-o" || cli_keys[i] == "-O") {
      amb_keys.push_back(cli_keys[i]);
    }
    else {
      nonamb_keys.push_back(cli_keys[i]);
    }
  }
  if (amb_keys.size() > 0) {
    addStandardAmberInputs(amb_keys);
  }

  // Parse the remaining inputs
  const size_t nrem = nonamb_keys.size();
  for (size_t i = 0; i < nrem; i++) {
    if (nonamb_keys[i] == "-t") {
      cli_nml.addKeyword("-t", NamelistType::STRING, std::string(""));
      cli_nml.addHelp("-t", "Name of the input transcript file, holding a detailed record of all "
                      "inputs given by the user as well as inputs which were possible with the "
                      "available namelist blocks");
    }
    else if (nonamb_keys[i] == "-a") {
      cli_nml.addKeyword("-a", NamelistType::STRING, std::string(""));
      cli_nml.addHelp("-a", "Name of the analysis output file, containing the results of inline "
                      "analyses performed by the program during the simulation / calculation.");
    }
    else if (nonamb_keys[i] == "-c_kind") {
      cli_nml.addKeyword("-c_kind", NamelistType::STRING,
                         getEnumerationName(default_filecon_inpcrd_type));
      cli_nml.addHelp("-c_kind", "Expected format of the input coordinates file");
    }
    else if (nonamb_keys[i] == "-x_kind") {
      cli_nml.addKeyword("-x_kind", NamelistType::STRING,
                         getEnumerationName(default_filecon_outcrd_type));
      cli_nml.addHelp("-x_kind", "Requested format of the output coordinates file");
    }
    else if (nonamb_keys[i] == "-r_kind") {
      cli_nml.addKeyword("-r_kind", NamelistType::STRING,
                         getEnumerationName(default_filecon_outcrd_type));
      cli_nml.addHelp("-r_kind", "Requested format of the checkpoint file");
    }
    else if (nonamb_keys[i] == "-except") {
      cli_nml.addKeyword("-except", NamelistType::STRING,
                         getEnumerationName(ExceptionResponse::DIE));
      cli_nml.addHelp("-except", "Action to take if invalid command-line input is encountered "
                      "(this may be passed down to other control blocks at the discretion of the "
                      "program's developer)");
    }
    else {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Unrecognized command line input keyword \"" + nonamb_keys[i] + "\".",
              "CommandLineParser", "addStandardApplicationInputs");
      case ExceptionResponse::WARN:
        rtWarn("Unrecognized command line input keyword \"" + nonamb_keys[i] + "\".  This input "
               "will be ignored.", "CommandLineParser", "addStandardApplicationInputs");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::addStandardApplicationInputs(const char* key_a, const char* key_b,
                                                     const char* key_c, const char* key_d) {
  addStandardApplicationInputs(vectorOfStrings(key_a, key_b, key_c, key_d));
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::addStandardApplicationInputs() {
  addStandardAmberInputs({ "-i", "-ig_seed", "-p", "-c", "-o", "-O", "-x" });
  addStandardApplicationInputs({ "-t", "-a", "-c_kind", "-x_kind", "-r_kind", "-except" });
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
  if ((help_on_no_args && argc == 1) || cli_nml.getBoolValue("--help") ||
      cli_nml.getBoolValue("-help")) {
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
      
      // Print the list of associated namelists, control blocks that can go in the program's
      // input file.
      if (control_blocks.size() > 0) {
        const int console_width = findFormatWidth(&std::cout);
        std::vector<std::string> all_blk = control_blocks;
        const TextFile mock_tf(std::string("Some text"), TextOrigin::RAM);
        int start_line = 0;
        bool found;
        for (size_t i = 0; i < control_blocks.size(); i++) {
          for (size_t j = 0; j < namelist_inventory.size(); j++) {
            if (strcmpCased(control_blocks[i], namelist_inventory[j].getTitle(),
                            CaseSensitivity::YES)) {
              const NamelistEmulator t_nml = namelist_inventory[j].invoke(mock_tf, &start_line,
                                                                          &found);
              all_blk.push_back(t_nml.getHelp());
            }
          }
          for (size_t j = 0; j < custom_namelists.size(); j++) {
            if (strcmpCased(control_blocks[i], custom_namelists[j].getTitle(),
                            CaseSensitivity::YES)) {
              const NamelistEmulator t_nml = custom_namelists[j].invoke(mock_tf, &start_line,
                                                                        &found);
              all_blk.push_back(t_nml.getHelp());
            }
          }
        }
        ReportTable tab_of_nml(all_blk, { "Namelist", "Description" }, std::string(""),
                               console_width, { JustifyText::LEFT, JustifyText::LEFT }, true);
        const std::string nml_msg = protectText("\nApplicable namelist control blocks (re-run "
                                                "with one of these titles as the command-line "
                                                "argument, in quotes if the leading '&' is "
                                                "included, for a full description of all keywords "
                                                "in the namelist):\n" +
                                                tab_of_nml.printTable(OutputSyntax::STANDALONE),
                                                ' ', console_width);
        std::cout << nml_msg << std::endl;
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

//-------------------------------------------------------------------------------------------------
void CommandLineParser::addControlBlocks(const std::vector<std::string> &list) {
  control_blocks.insert(control_blocks.end(), list.begin(), list.end());

  // Insert an ampersand & at the front of each namelist title
  const size_t nblk = control_blocks.size();
  for (size_t i = 0; i < nblk; i++) {
    if (control_blocks[i][0] != '&') {
      control_blocks[i] = "&" + control_blocks[i];
    }
  }
}

//-------------------------------------------------------------------------------------------------
void CommandLineParser::addCustomNamelists(const std::vector<NamelistToken> &custom_namelists_in) {
  custom_namelists.insert(custom_namelists.end(), custom_namelists_in.begin(),
                          custom_namelists_in.end());
}

} // namespace namelist
} // namespace stormm
