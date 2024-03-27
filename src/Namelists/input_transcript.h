// -*-c++-*-
#ifndef STORMM_INPUT_TRANSCRIPT_H
#define STORMM_INPUT_TRANSCRIPT_H

#include <limits.h>
#include <string>
#include <vector>
#include "copyright.h"
#include "user_settings.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

/// \brief Form a string for printing the command line based on the series of CLI arguments.
///
/// \param ui     Primary compendium of user input, containing a list of the command line arguments
/// \param width  Width of horizontal rule to print ahead of the command line.  If set to less than
///               or equal to zero, no horizontal rule will be printed.
std::string commandLineAsString(const UserSettings &ui, int width = 0);
  
/// \brief Transcribe input values from all applicable namelists to a file that will aid users in
///        understanding how their specific inputs may have interacted with pre-existing defaults.
///
/// Overloaded:
///   - Write to a string or file
///   - Assume that the current terminal width is the desired file width
///   - Indicate a preferred file width
///
/// \param ui               The primary compendium of user input, containing most namelists
///                         commonly used by various STORMM applications
/// \param width            The intended width of the output file
/// \param max_repetitions  The maximum number of repetitions that will be reported for any one
///                         keyword entry
/// \param extra_nml        Additional namelists, such as those defined for a specific application
/// \{
std::string prepareInputTranscript(const UserSettings &ui, int width,
                                   int max_repetitions = INT_MAX,
                                   const std::vector<NamelistEmulator> &extra_nml = {});

void writeInputTranscript(const UserSettings &ui, int width, int max_repetitions = INT_MAX,
                          const std::vector<NamelistEmulator> &extra_nml = {});

void writeInputTranscript(const UserSettings &ui, int max_repetitions = INT_MAX,
                          const std::vector<NamelistEmulator> &extra_nml = {});
/// \}

} // namespace namelist
} // namespace stormm

#endif
