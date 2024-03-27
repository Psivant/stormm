// -*-c++-*-
#ifndef STORMM_INPUT_H
#define STORMM_INPUT_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/textfile.h"
#include "Parsing/textguard.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using parse::TextFile;
using parse::TextGuard;
using parse::WrapTextSearch;

/// \brief Pad special characters including [], {}, and () enclosures or separators '=', ',', and
///        ';' (so long as they are not protected by quotations or comments).  The data comes in as
///        a C-style character string (for ease of calling the function based on TextFile data) and
///        the result is returned as a C++ string.
/// 
/// \param text    Raw character string of the text to parse
/// \param n_char  Number of characters in the string
std::string padNamelistTuples(const char* text, const int n_char);

/// \brief Prepare a list of std::strings representing a sequence of delimited text strings
///        found in a TextFile object.  Contiguous blocks of quoted or commented text represent a
///        single such text string, and the quotation (" (...) ", ' (...) ') or comment character
///        options (#, !, %, //, /* (...) */ ) for STORMM namelists are hard-coded herein.  In this
///        notation (...) is not a reserved character sequence but indicates that the quote or
///        comment feature can span multiple lines.  Other special characters '[', '{', '(', ')',
///        '}', and ']', representing the limits of argument tuples, will be padded with
///        whitespace for this parsing procedure unless they are found in quotes or comments in
///        the original TextFile.  Finally, '=' and ',' characters, so long as they do not appear
///        in quotes or comment blocks, will be removed so as not to appear in the output.  The
///        result can then be interpreted as a series of "keyword -> value" or "keyword -> list of
///        values in [], {}, or ()" pairs.
///
/// \param tf          Text of some input deck committed to RAM
/// \param nml         The namelist to search for and then fill
/// \param start_line  The first line of the namelist (output from detectNamelist() above).
///                    Defaults to zero to begin at the head of the file and find the first
///                    namelist matching the appropriate title.
/// \param wrap        Indicator of whether to start searching the beginning of the file for a
///                    namelist.  Defaults to NO, as it may be the case that other such namelists
///                    of the same name have already been found and should not be found again.
/// \param end_line    The line at which to terminate the search for a namelist.  Defaults to -1,
///                    indicating the end of the file.  Developers will probably not often use this
///                    parameter as its primary purpose is to stop a namelist search at the
///                    original start line after wrapping back around to the beginning of the file.
std::vector<std::string> pullNamelist(const TextFile &tf, const NamelistEmulator &nml,
                                      int start_line, WrapTextSearch wrap, int end_line);

/// \brief Load a namelist with user input obtained from a TextFile object.  Returns the next line
///        after the end of the first example of the namelist.
///
/// \param tf          Text of some input deck committed to RAM
/// \param nml         The namelist to search for and then fill
/// \param start_line  The first line of the namelist (output from detectNamelist() above).
///                    Defaults to zero to begin at the head of the file and find the first
///                    namelist matching the appropriate title.
/// \param wrap        Indicator of whether to start searching the beginning of the file for a
///                    namelist.  Defaults to NO, as it may be the case that other such namelists
///                    of the same name have already been found and should not be found again.
/// \param end_line    The line at which to terminate the search for a namelist.  Defaults to -1,
///                    indicating the end of the file.  Developers will probably not often use this
///                    parameter as its primary purpose is to stop a namelist search at the
///                    original start line after wrapping back around to the beginning of the file.
/// \param found       A way to pass back the information that a namelist was found in the input
///                    (much simpler and safer than trying to infer such a fact from changes in
///                    start_line).  Leaving the variable set to nullptr disables this feature.
int readNamelist(const TextFile &tf, NamelistEmulator *nml, int start_line = 0,
                 WrapTextSearch wrap = WrapTextSearch::NO, int end_line = -1,
                 bool *found = nullptr);
  
} // namespace namelist
} // namespace stormm

#endif
