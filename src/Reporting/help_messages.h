// -*-c++-*-
#ifndef STORMM_HELP_MESSAGES_H
#define STORMM_HELP_MESSAGES_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Namelists/namelist_emulator.h"
#include "reporting_enumerators.h"

namespace stormm {
namespace display {

using namelist::NamelistEmulator;
  
/// \brief Detect a command line help signal.  This is then relayed to the calling program so that
///        its internal help commands can be displayed.
///
/// \param argc        Number of command-line arguments (passed in from main)
/// \param argv        List of command-line arguments (passed in from main)
/// \param kind        Type of signal that is authorized to lead to a help message
/// \param help_words  Keywords that can be used to signal the call for help.  Custom keywords will
///                    be assumed to be case-sensitive.  If no list is provided, the words "help",
///                    "-help", and "--help" will trigger a message.
bool detectHelpSignal(const int argc, const char* argv[],
                      HelpSignalKind kind = HelpSignalKind::KEYWORD,
                      const std::vector<std::string> &help_words = {});

/// \brief Display the documentation for one of the namelists in the STORMM libraries.  This can
///        be called with a simple keyword corresponding to the name of the namelist, the name
///        being case-sensitive.
///
/// Overloaded:
///   - Provide the name of a namelist (must correspond to one of the namelists in STORMM's main
///     libraries)
///   - Provide command line inputs and the name of a namelist
///   - Provide command line inputs and a list of namelist names
///
/// \param module_name  Name of the namelist or namelists (these must be preceded by the '&' sign,
///                     to allow as much leeway as possible in other command-line arguments)
/// \param t_nml        An actual example of the namelist, provided to support functionality in
/// \{
bool displayNamelistHelp(const int argc, const char* argv[],
                         const std::vector<std::string> &module_name);

bool displayNamelistHelp(const int argc, const char* argv[], const std::string &module_name);

bool displayNamelistHelp(const std::string &module_name);
/// \}

} // namespace display
} // namespace stormm

#endif
