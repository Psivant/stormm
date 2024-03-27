// -*-c++-*-
#ifndef STORMM_DISSECT_TEXTFILE_H
#define STORMM_DISSECT_TEXTFILE_H

#include <algorithm>
#include <regex>
#include <string>
#include "copyright.h"
#include "Constants/behavior.h"
#include "DataTypes/stormm_vector_types.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/textfile.h"
#include "Parsing/textguard.h"
#include "Reporting/reporting_enumerators.h"

namespace stormm {
namespace testing {

using constants::CaseSensitivity;
#ifndef STORMM_USE_HPC
using data_types::double2;
#endif
using parse::TextFile;
using parse::TextGuard;
using parse::TextOrigin;
using review::TextEnds;

/// \brief Count the average and standard deviation of the number of words separated by white space
///        per line in a TextFile object.
///
/// \param tf  The text of interest
double2 meanWordsPerLine(const TextFile &tf);

/// \brief Count the number of instances of a word or sequence of words in a text file.  This
///        function can be set to ignore specific delimiters.
///
/// Overloaded:
///   - Count the number of instances of a single word
///   - Count the number of instances of a sequence of words
///
/// \param tf          The text of interest
/// \param sequence    The word or sequence of words
/// \param separators  Separators to identify when searching for wrods.  Separators are not ignored
///                    in the search, but instead given special consideration and set apart with
///                    white space to ensure that they are noticeable and do not count against
///                    words on either side.
/// \{
int countInstances(const TextFile &tf, const std::string &sequence,
                   const std::vector<char> &separators = {},
                   const std::vector<char> &delimiters = {},
                   const std::vector<TextGuard> &quote_marks = {},
                   const std::vector<TextGuard> &comment_marks = {});

int countInstances(const TextFile &tf, const std::vector<std::string> &sequence,
                   const std::vector<char> &separators = {},
                   const std::vector<char> &delimiters = {},
                   const std::vector<TextGuard> &quote_marks = {},
                   const std::vector<TextGuard> &comment_marks = {});
/// \}

/// \brief Operate on a TextFile in the program memory as if using the standard Unix / Linux grep
///        command.  Return a new TextFile containing the line-by-line results.
///
/// Overloaded:
///   - Provide a Standard Tempalte Library regex object as the basis of the search
///   - Provide a Standard Template Library string object as the basis of the search
///   - Provide a C-style character array as the basis of the search
///
/// \param tf         The original text from within program memory
/// \param criterion  Formulation of the regular expression to search for
/// \param length     Trusted length of criterion, if provided as a C-style array
/// \param before     In addition to each line on which the expression is matched, return this many
///                   lines appearing before the matching line
/// \param after      In addition to each line on which the expression is matched, return this many
///                   lines appearing after the matching line
/// \param cased      Specify whether the search will be case sensitive
/// \{
TextFile grep(const TextFile &tf, const std::regex &criterion, int before = 0, int after = 0,
              CaseSensitivity cased = CaseSensitivity::YES,
              const std::vector<TextGuard> &quote_marks = {},
              const std::vector<TextGuard> &comment_marks = {});

TextFile grep(const TextFile &tf, const std::string &criterion, int before = 0, int after = 0,
              CaseSensitivity cased = CaseSensitivity::YES,
              const std::vector<TextGuard> &quote_marks = {},
              const std::vector<TextGuard> &comment_marks = {});

TextFile grep(const TextFile &tf, const char* criterion, const int length, int before = 0,
              int after = 0, CaseSensitivity cased = CaseSensitivity::YES,
              const std::vector<TextGuard> &quote_marks = {},
              const std::vector<TextGuard> &comment_marks = {});

} // namespace testing
} // namespace stormm

#endif
