// -*-c++-*-
#ifndef STORMM_COMMAND_LINE_PARSER_H
#define STORMM_COMMAND_LINE_PARSER_H

#include <map>
#include <string>
#include <vector>
#include "copyright.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

/// \brief Default values for STORMM command line inputs
/// \{
const char default_stormm_input_deck[] = "stormm.in";
const char default_stormm_input_topology[] = "prmtop";
const char default_stormm_input_coordinates[] = "inpcrd";
const char default_stormm_report_file[] = "stormm_report.m";
const char default_stormm_output_trajectory[] = "mdcrd";
const char default_stormm_output_checkpoint[] = "mdrst";  
/// \}

/// \brief A class for collecting command line information.  The class object will function
///        somewhat like a namelist in that it gets loaded with keywords and can return the help
///        messages for its contents upon receiving such a directive, but unlike namelists this
///        one object is designed to be reconfigured for specific programs (even test programs,
///        when it appears nested within the TestEnvironment class).
class CommandLineParser {
public:

  /// \brief The constructor is passed the list of command-line arguments.  The object will take
  ///         the first string in the list as the name of the calling program itself.
  ///
  /// \param program_description  A help message tailored for the program as a whole
  CommandLineParser(const std::string &program_name_in, const std::string &program_description,
                    const std::vector<std::string> &noted_imports_in = {},
                    ExceptionResponse policy_in = ExceptionResponse::WARN);

  /// \brief With no pointers to repair and all members coming from or built using Standard
  ///        Template Library objects, the default copy and move constructors, as well as copy and
  ///        move assignment operators, are valid.
  ///
  /// \param original  THe original object to copy or move
  /// \param other     Another object placed on the right hand side of the assignment statement
  /// \{
  CommandLineParser(const CommandLineParser &original) = default;
  CommandLineParser(CommandLineParser &&original) = default;
  CommandLineParser& operator=(const CommandLineParser &original) = default;
  CommandLineParser& operator=(CommandLineParser &&original) = default;
  /// \}

  /// \brief Activate the printout of help messages if no command line arguments are provided.
  void activateHelpOnNoArgs();

  /// \brief Suppress the printout of help messages if no command line arguments are provided.
  void suppressHelpOnNoArgs();

  /// \brief Order the program to exit after displaying general help documentation, with a return
  ///        value of zero.
  void activateExitOnHelp();

  /// \brief Prevent the program from exiting after displaying general help documentation, with a
  ///        return value of zero.
  void preventExitOnHelp();

  /// \brief Impart the typical comand-line keywords inspired by Amber.
  ///
  /// Overloaded:
  ///   - Provide a vector of many Amber input keys
  ///   - Provide a single key or a handful of keys (for developers who forget their { } braces)
  ///
  /// \param cli_keys  The list of keywords, each of which must match a recognized Amber command
  ///                  line input keyword from one of the major engines
  /// \param key_a     The first of up to four specific keys
  /// \param key_b     The second of up to four specific keys
  /// \param key_c     The third of up to four specific keys
  /// \param key_d     The fourth of up to four specific keys
  /// \{
  void addStandardAmberInputs(const std::vector<std::string> &cli_keys);
  void addStandardAmberInputs(const char* key_a, const char* key_b = nullptr,
                              const char* key_c = nullptr, const char* key_d = nullptr);
  /// \}

  /// \brief Inpart some common benchmarking keywords to the command line interface.  Overloading
  ///        and descriptions of input parameters follow from addStandardAmberInputs(), above.
  /// \{
  void addStandardBenchmarkingInputs(const std::vector<std::string> &cli_keys);
  void addStandardBenchmarkingInputs(const char* key_a, const char* key_b = nullptr,
                                     const char* key_c = nullptr, const char* key_d = nullptr);  
  /// \}
  
  /// \brief Return a pointer to the internal namelist.
  ///
  /// Overloaded:
  ///   - Return a const pointer to the namelist emulator for a const object.
  ///   - Return a mutable pointer to the namelist emulator for a non-const object.
  /// \{
  const NamelistEmulator* getNamelistPointer() const;
  NamelistEmulator* getNamelistPointer();
  /// \}

  /// \brief Supply the internal namelist with command-line input as if parsing from an input
  ///        file.
  ///
  /// \param argc  The number of command-line arguments
  /// \param argv  List of character strings obtained from the command line
  void parseUserInput(int argc, const char* argv[]);

  /// \brief The nature of the command line parser as a unique sort of namelist which can be
  ///        assembled and recombined from parts also suggests a utility in making multiple
  ///        CommandLineParser objects in different parts of a program, which might all look at
  ///        the same command line inputs but take in different information.  In order for this to
  ///        work, each object must be able to know that some other object will handle certain
  ///        inputs, which should then not trigger warnings or errors.
  ///
  /// \param other  Another object with which the CommandLineParser should coordinate when parsing
  ///               command line inputs
  void coordinateWithPartner(CommandLineParser *other);
  
private:

  ExceptionResponse policy;  ///< The course of action to take in the event of bad user input, or
                             ///<   even bad developer programming.  This will filter down to the
                             ///<   underlying namelist object.
  int arg_count;             ///< Number of command line arguments (less the program name itself)
  std::string program_name;  ///< The name of the program itself, as provided by the developer
  std::string executable;    ///< The path to the program, as found in the command line call
  NamelistEmulator cli_nml;  ///< Namelist created to process all command-line data
  bool help_on_no_args;      ///< Indicate whether help messages should be produced if no command
                             ///<   line input (other that the program name) is provided
  bool exit_on_help;         ///< Indicate that the program should exit after displaying a help
                             ///<   message, with a return value of zero
  bool lead_parser;          ///< Indicate that this is the leader of many coordinated parsers.
                             ///<   This flag will be set to TRUE on construction and FALSE if any 
                             ///<   other CommandLineParser object calls its
                             ///<   coordinateWithPartner() member function against this one.
  
  /// A list of all command line arguments, kept for documenting user activity in report files
  std::vector<std::string> command_line_text;

  /// List of all other class objects to which this object is coordinated
  std::vector<CommandLineParser*> coordinations;

  /// List of all excluded keywords
  std::map<std::string, NamelistType> excluded_keys;

  /// List of all notable keywords from coordinating CommandLineParser objects
  std::vector<std::string> noted_imports;
};
  
} // namespace namelist
} // namespace stormm

#endif
