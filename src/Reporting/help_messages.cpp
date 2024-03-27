#include <cstring>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Namelists/nml_conformer.h"
#include "Namelists/nml_dynamics.h"
#include "Namelists/nml_ffmorph.h"
#include "Namelists/nml_files.h"
#include "Namelists/nml_minimize.h"
#include "Namelists/nml_precision.h"
#include "Namelists/nml_random.h"
#include "Namelists/nml_report.h"
#include "Namelists/nml_restraint.h"
#include "Namelists/nml_solvent.h"
#include "Parsing/parse.h"
#include "Parsing/textfile.h"
#include "help_messages.h"

namespace stormm {
namespace display {

using constants::CaseSensitivity;
using namelist::conformerInput;
using namelist::dynamicsInput;
using namelist::ffmorphInput;
using namelist::filesInput;
using namelist::minimizeInput;
using namelist::precisionInput;
using namelist::randomInput;
using namelist::reportInput;
using namelist::restraintInput;
using namelist::solventInput;
using namelist::KeyRequirement;
using parse::strcmpCased;
using parse::TextFile;
using parse::TextOrigin;

//-------------------------------------------------------------------------------------------------
bool detectHelpSignal(const int argc, const char* argv[], const HelpSignalKind kind,
                      const std::vector<std::string> &help_words) {
  const int nh_words = help_words.size();

  // If the only way to get help is by supplying no arguments, return TRUE or FALSE based on
  // whether there are no command line arguments.  If the help is triggered by a keyword, but not
  // keyword only, assume that running the program with no arguments is also intended to trigger
  // the help messages.
  switch (kind) {
  case HelpSignalKind::NO_ARGS_ONLY:
    return (argc == 1);
  case HelpSignalKind::NO_ARGS:    
  case HelpSignalKind::KEYWORD:
    if (argc == 1) {
      return true;
    }
    break;
  case HelpSignalKind::KEYWORD_ONLY:
    break;
  }

  // If there are keywords to parse, cases where the program is supposed to print help messages
  // only in the case that no command line arguments are supplied would already have receive a
  // signal of FALSE.  However, if the program is set to provide help messages with a particular
  // command line argument, parse each argument here.
  switch (kind) {
  case HelpSignalKind::NO_ARGS_ONLY:
    break;
  case HelpSignalKind::NO_ARGS:
  case HelpSignalKind::KEYWORD:
  case HelpSignalKind::KEYWORD_ONLY:  
    for (int i = 0; i < argc; i++) {
      if (nh_words == 0) {
        if (strcmpCased(argv[i], "help", CaseSensitivity::NO) ||
            strcmpCased(argv[i], "-help", CaseSensitivity::NO) ||
            strcmpCased(argv[i], "--help", CaseSensitivity::NO)) {
          return true;
        }
      }
      else {
        for (int j = 0; j < nh_words; j++) {
          if (strcmpCased(argv[i], help_words[j], CaseSensitivity::YES)) {
            return true;
          }
        }
      }
    }
    break;
  }
  return false;
}

//-------------------------------------------------------------------------------------------------
bool displayNamelistHelp(const int argc, const char* argv[],
                         const std::vector<std::string> &module_name) {
  const int n_names = module_name.size();
  bool found = false;
  for (int i = 0; i < n_names; i++) {
    for (int j = 0; j < argc; j++) {
      if (strcmp(argv[j], module_name[i].c_str()) == 0) {
        found = (found || displayNamelistHelp(module_name[i]));
      }
    }
  }
  return found;
}

//-------------------------------------------------------------------------------------------------
bool displayNamelistHelp(const int argc, const char* argv[], const std::string &module_name) {
  return displayNamelistHelp(argc, argv, std::vector<std::string>(1, module_name));
}

//-------------------------------------------------------------------------------------------------
bool displayNamelistHelp(const std::string &module_name) {
  const TextFile tf(std::string(""), TextOrigin::RAM);
  int start_line = 0;
  if (strcmpCased(module_name, "&files", CaseSensitivity::YES)) {
    const std::vector<KeyRequirement> sys_keyword_reqs(1, KeyRequirement::REQUIRED);
    const NamelistEmulator t_nml = filesInput(tf, &start_line, nullptr, sys_keyword_reqs);
    t_nml.printHelp();
    return true;
  }
  else if (strcmpCased(module_name, "&restraint", CaseSensitivity::YES)) {
    const NamelistEmulator t_nml = restraintInput(tf, &start_line, nullptr);
    t_nml.printHelp();
    return true;
  }
  else if (strcmpCased(module_name, "&solvent", CaseSensitivity::YES)) {
    const NamelistEmulator t_nml = solventInput(tf, &start_line, nullptr);
    t_nml.printHelp();
    return true;
  }
  else if (strcmpCased(module_name, "&minimize", CaseSensitivity::YES)) {
    const NamelistEmulator t_nml = minimizeInput(tf, &start_line, nullptr);
    t_nml.printHelp();
    return true;
  }
  else if (strcmpCased(module_name, "&conformer", CaseSensitivity::YES)) {
    const NamelistEmulator t_nml = conformerInput(tf, &start_line, nullptr);
    t_nml.printHelp();
    return true;
  }
  else if (strcmpCased(module_name, "&dynamics", CaseSensitivity::YES)) {
    const NamelistEmulator t_nml = dynamicsInput(tf, &start_line, nullptr);
    t_nml.printHelp();
    return true;
  }
  else if (strcmpCased(module_name, "&report", CaseSensitivity::YES)) {
    const NamelistEmulator t_nml = reportInput(tf, &start_line, nullptr);
    t_nml.printHelp();
    return true;
  }
  else if (strcmpCased(module_name, "&random", CaseSensitivity::YES)) {
    const NamelistEmulator t_nml = randomInput(tf, &start_line, nullptr);
    t_nml.printHelp();
    return true;
  }
  else if (strcmpCased(module_name, "&precision", CaseSensitivity::YES)) {
    const NamelistEmulator t_nml = precisionInput(tf, &start_line, nullptr);
    t_nml.printHelp();
    return true;
  }
  else if (strcmpCased(module_name, "&ffmorph", CaseSensitivity::YES)) {
    const NamelistEmulator t_nml = ffmorphInput(tf, &start_line, nullptr);
    t_nml.printHelp();
    return true;
  }
  else {
    rtErr("No namelist " + module_name + " is known in the STORMM libraries.",
          "displayNamelistHelp");
  }
  __builtin_unreachable();
}

} // namespace display
} // namespace stormm
