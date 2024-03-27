// -*-c++-*-
#ifndef STORMM_PARSE_H
#define STORMM_PARSE_H

#include <vector>
#include <string>
#include "copyright.h"
#include "Constants/behavior.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "parsing_enumerators.h"
#include "polynumeric.h"
#include "textfile.h"
#include "textguard.h"

namespace stormm {
namespace parse {

using constants::ExceptionResponse;
using constants::CaseSensitivity;
using data_types::isSignedIntegralScalarType;
using data_types::isUnsignedIntegralScalarType;
using data_types::isFloatingPointScalarType;
  
/// \brief Constant to define the default formatting for integer and real number representations.
///        If these values are seen then the program knows to do free formatting.
constexpr int free_number_format = -32788;

/// \brief Overload the + operator to extend a string with a char4 object.  The x, y, z, and w
///        components of the char4 will continue to extend the string until one of them is the
///        terminating null character.
///
/// \param lhs  The string to extend
/// \param rhs  The char tuple to add
std::string operator+(const std::string &lhs, const char4 rhs);

/// \brief Convert a character tuple to a string.  This is slightly faster than creating an empty
///        string and then using the += operator.  As above, the x, y, z, and w components of the
///        character tuple will continue to extend the string unless a null character is reached.
///
/// Overloaded:
///   - Accept a char2 tuple
///   - Accept a char3 tuple
///   - Accept a char4 tuple
///
/// \param value  The character tuple to convert
/// \{
std::string char2ToString(const char2 value);
std::string char3ToString(const char3 value);
std::string char4ToString(const char4 value);
/// \}

/// \brief Convert the first four characters of a string to char4, ignoring subsequent characters
///        and padding with whitespace if the string is shorter than four characters.
///
/// \param value  The string to convert
char4 stringToChar4(const std::string &value);
  
/// \brief Convert a positive integer to an alphabetic string, essentially a base-26 number.
///
/// \param input  The integer to convert
std::string alphabetNumber(ullint input);

/// \brief Convert a series of unsigned integers to a hexadecimal code indicative of an RGB color.
///
/// \param color_in  The RGB code to convert, with the 'R', 'G', and 'B' colors in the "x", "y",
///                  and "z" members, respectively.  The "w" member is ignored.
std::string rgbHexCode(const uchar4 color_in);
  
/// \brief Determine whether a character string can qualify as an integer, real, or char4
///
/// \param a           The character string
/// \param cform       The expected format of the number
/// \param read_begin  Point in the string at which to being reading the number (defaults to the
///                    start of the string)
/// \param len         Expected length of the number (default 0, continues until the end of the
///                    string)
bool verifyNumberFormat(const char* a, NumberFormat cform, int read_begin = 0, int len = 0);

/// \brief Some text or input files may contain missing or blank entries.  Check that the space
///        to be read for input contains data of the appropriate type.  This function feeds into
///        verifyNumberFormat() above, but may apply additional checks on the length of the
///        available data and return false if there is nothing to examine.
///
/// Overloaded:
///   - Inspect a C-style array of characters
///   - Inspect a Standard Template Library string object
///   - Inspect a TextFile object, or its abstract, with a line number, starting position, and
///     length
///
/// \param line       The line data as a std::string object, a C-style array of characters, or the
///                   number of the line in a TextFile object
/// \param tf         Text file object committed to RAM with the data of interest
/// \param start_pos  Position of the first character on the line to read
/// \param length     Length of the character string to analyze
/// \param fmt        Format of the text to search for
/// \{
bool verifyContents(const char* line, int start_pos, int length,
                    NumberFormat fmt = NumberFormat::CHAR4);
bool verifyContents(const std::string &line, NumberFormat fmt);
bool verifyContents(const std::string &line, int start_pos, int length,
                    NumberFormat fmt = NumberFormat::CHAR4);
bool verifyContents(const TextFile &tf, int line, int start_pos, int length,
                    NumberFormat fmt = NumberFormat::CHAR4);
bool verifyContents(const TextFileReader &tfr, int line, int start_pos, int length,
                    NumberFormat fmt = NumberFormat::CHAR4);
/// \}

/// \brief Convert one or more characters to uppercase, if they are letters between 'a' and 'z'.
///        Otherwise, return the original character(s).
///
/// Overloaded:
///   - Convert a single character
///   - Convert a C++ string
///   - Convert a C-style string (given a limit), in place
///   - Convert a constant C-style string, returning a C++ string
///
/// \param tc      The character in question
/// \param tcs     The character string in question
/// \param n_char  The number of characters to convert in the C-style string (trusted that the
///                string has this many characters)
/// \{
char uppercase(const char tc);
std::string uppercase(const std::string &ts);
void uppercase(char* tcs, size_t n_char = 0);
std::string uppercase(const char* tcs);
/// \}

/// \brief Convert one or more characters to lowercase, if they are letters between 'a' and 'z'.
///        Otherwise, return the original character(s).
///
/// Overloaded:
///   - Convert a single character
///   - Convert a C++ string
///   - Convert a C-style string (given a limit), in place
///   - Convert a constant C-style string, returning a C++ string
///
/// \param tc   The character in question
/// \param tcs     The character string in question
/// \param n_char  The number of characters to convert in the C-style string (trusted that the
///                string has this many characters)
/// \{
char lowercase(const char tc);
std::string lowercase(const std::string &ts);
void lowercase(char* tcs, size_t n_char = 0);
std::string lowercase(const char* tcs);
/// \}

/// \brief Match two strings based on a particular case sensitivity setting.
///
/// Overloaded:
///   - Take two C-strings (the fundamental case)
///   - Take a C-string and a C++ string (in either order)
///   - Take two C++ strings
///
/// \param sa      The first string
/// \param sb      The second string
/// \param csen    Case sensitivity setting
/// \{
bool strcmpCased(const char* sa, const char* sb,
                 const CaseSensitivity csen = CaseSensitivity::YES);
  
bool strcmpCased(const std::string &sa, const char* sb,
                 const CaseSensitivity csen = CaseSensitivity::YES);

bool strcmpCased(const char* sa, const std::string &sb,
                 const CaseSensitivity csen = CaseSensitivity::YES);

bool strcmpCased(const std::string &sa, const std::string &sb,
                 const CaseSensitivity csen = CaseSensitivity::YES);
/// \}

/// \brief Match two strings, up to a given length, based on a particular case sensitivity setting.
///
/// Overloaded:
///   - Take two C-strings (the fundamental case)
///   - Take a C-string and a C++ string (in either order)
///   - Take two C++ strings
///
/// \param sa      The first string
/// \param sb      The second string
/// \param csen    Case sensitivity setting
/// \param length  The number of characters to compare (if 0 and a C++ string is provided, the
///                length of the string will be taken)
/// \{
bool strncmpCased(const char* sa, const char* sb, int length,
                  CaseSensitivity csen = CaseSensitivity::YES);

bool strncmpCased(const std::string &sa, const char* sb,
                  CaseSensitivity csen = CaseSensitivity::YES, int length = -1);

bool strncmpCased(const char* sa, const std::string &sb,
                  CaseSensitivity csen = CaseSensitivity::YES, int length = -1);

bool strncmpCased(const std::string &sa, const std::string &sb,
                  CaseSensitivity csen = CaseSensitivity::YES, int length = -1);
/// \}

/// \brief Match two strings, the second of which may contain wildcards.
///
/// \param target     The literal string, containing no wildcards, which must be matched
/// \param query      The match string, which may contain wildcards
/// \param wildcards  Array of indicators as to whether the indices in query can be treated as
///                   wildcards, and if so of what type
bool strcmpWildCard(const std::string &target, const std::string &query,
                    const std::vector<WildCardKind> &wildcards);

/// \brief Determine the appropriate number of decimal places in which to display a real number
///        without superfluous zeros.
///
/// \param value  The number to display
/// \param limit  The maximum number of decimal places to report (optional)
int realDecimalPlaces(double value, int limit = 10);

/// \brief Add leading white space to a string and return the new string as the result.
///
/// Overloaded:
///   - Return the resulting string as a new string
///   - Modify an existing string
///
/// \param input          The string to pre-pend
/// \param target_length  The length that the string should have with added white space.  If the
///                       string is already this length or longer, no white space will be added.
/// \{
std::string addLeadingWhiteSpace(const std::string &input, const size_t target_length);
void addLeadingWhiteSpace(std::string *input, const size_t target_length);
/// \}

/// \brief Add tailing white space to a string and return the new string as the result.
///
/// Overloaded:
///   - Return the resulting string as a new string
///   - Modify an existing string
///
/// \param input          The string to pre-pend
/// \param target_length  The length that the string should have with added white space.  If the
///                       string is already this length or longer, no white space will be added.
/// \{
std::string addTailingWhiteSpace(const std::string &input, const size_t target_length);
void addTailingWhiteSpace(std::string *input, const size_t target_length);
/// \}

/// \brief Remove leading white space from a string and return a new string as the result.
///
/// Overloaded:
///   - Return the resulting string as a new string
///   - Modify an existing string
///
/// \param input  The string to prune
/// \{
std::string removeLeadingWhiteSpace(const std::string &input);
void removeLeadingWhiteSpace(std::string *input);
/// \}

/// \brief Remove tailing white space from a string and return a new string as the result.
///
/// Overloaded:
///   - Return the resulting string as a new string
///   - Modify an existing string
///
/// \param input  The string to prune
/// \{
std::string removeTailingWhiteSpace(const std::string &input);
void removeTailingWhiteSpace(std::string *input);
/// \}

/// \brief Add leading white space to a series of strings.  White space will be pre-pended to the
///        front or back, insofar as it does not add excessive blank text to the output.  The
///        consensus width of all default parameter strings is returned.
///
/// Overloaded:
///   - Operate on a range of strings within the array
///   - Operate on the entrie array
///
/// \param ds           List of strings, such as keyword parameters
/// \param lower_limit  The lower limit of strings to operate upon
/// \param upper_limit  The upper limit of strings to operate upon
/// \param max_length   The maximum length to raise all strings to without further consideration.
///                     If some strings are above this length, at least long_count members of the
///                     list must aproach the upper end of the length spectrum before all strings
///                     will be lengthened.
/// \param long_count   The number of strings that must come within near_limit of the longest
///                     string in order for all strings to be lengthed to the maximum length found.
/// \param near_limit   The number of characters that at least long_count - 1 strings, other than
///                     the longest string found, must come within the length of the longest string
///                     in order for all strings to be lengthend to max_length.  
/// \{
size_t justifyStrings(std::vector<std::string> *ds, int lower_limit, int upper_limit,
                      JustifyText fside = JustifyText::LEFT, size_t max_length = 10,
                      int long_count = 2, int near_limit = 4);

size_t justifyStrings(std::vector<std::string> *ds, JustifyText fside = JustifyText::LEFT,
                      size_t max_length = 10, int long_count = 2, int near_limit = 4);
/// \}

/// \brief Compute the maximum length of an individual line, if a string were to be printed to the
///        screen or a file with its carriage returns intact.
///
/// \param input  The string to examine
size_t maximumLineLength(const std::string &input);
  
/// \brief Convert a real number to a formatted string.
///
/// Overloaded:
///   - Provide both format specifiers.  Printing method defaults to %(a).(b)lf.
///   - Provide just one format specifier (must also provide the printing method)
///
/// \param value     The number to convert
/// \param format_a  First of up to two format specifiers.  If this is given and the second format
///                  specifier is blank, it refers to the number of decimal places to print for
///                  the number.  If the second format specifier is also present, this becomes the
///                  total width of the general format number.
/// \param format_b  If given, this will convert the first format specifier into the total width
///                  and then specify the number of digits after the decimal.
/// \param style     Whether to print leading zeros when there is more format space than digits
/// \{
std::string realToString(double value, int format_a = free_number_format,
                         int format_b = free_number_format,
                         NumberFormat method = NumberFormat::STANDARD_REAL,
                         NumberPrintStyle = NumberPrintStyle::STANDARD);

std::string realToString(const double value, const int format_a, const NumberFormat method,
                         const NumberPrintStyle style = NumberPrintStyle::STANDARD);
/// \}

/// \brief Return the shortest string able to represent a real number to a relative precision,
///        preferring standard notation over scientific in the case that both representations
///        would be the same length.
///
/// \param value  The value to convert
/// \param rel    Precision to which the value should be represented
std::string minimalRealFormat(const double value, const double rel);

/// \brief Determine the format width needed to display a series of numbers with consistent
///        alignment.
///
/// Overloaded:
///   - Operate on a C-style array with trusted length
///   - Operate on a Standard Template Library vector object
///
/// \param numbers         The series of numbers to be printed in ASCII, general numeric format
///                        (this is STANDARD_REAL in the NumberFormat enumerations)
/// \param count           Trusted length of numbers, if providing a C-style array
/// \param decimal_places  The number of digits after the decimal
/// \param advance         Advance this many numbers with each step through the sequence
/// \{
template <typename T>
int findAlignmentWidth(const T* numbers, size_t count, int decimal_places, size_t advance = 1);

template <typename T>
int findAlignmentWidth(const std::vector<T> &numbers, int decimal_places, size_t advance = 1);
/// \}
  
/// Convert an integer to a formatted string.
///
/// \param value  The number to convert
/// \param width  The format specifier, showing how many characters to give the integer
/// \param style  Whether to print leading zeros when there is more format space than digits
std::string intToString(llint value, int width = free_number_format,
                        NumberPrintStyle = NumberPrintStyle::STANDARD);

/// \brief Overload the + operator to make it easier to append sequences of TextGuard objects in
///        std::vectors.
const std::vector<TextGuard> operator+(const std::vector<TextGuard> &lhs,
                                       const std::vector<TextGuard> &rhs);

/// \brief Decide whether the sequence of characters in a TextFile beginning at a given position
///        completes a guarded section of the text, given the guard sequence that started the
///        section.
///
/// Overloaded:
///   - Detect guards within a text file
///   - Detect guards within a character array of trusted length
///
/// \param tfr        Text file reading abstract
/// \param line_idx   Index at which to seek the putative termination sequence
/// \param pos_idx    Index at which to seek the putative termination sequence
/// \param
/// \param guard_seq  Text guard containing the termination sequence
/// \{
bool detectGuard(const TextFileReader &tfr, const int line_idx, const int pos_idx,
                 const std::string &guard_seq);

bool detectGuard(const char* line, const int pos_idx, const std::string &guard_seq);

bool detectGuard(const char* line, size_t pos_idx, size_t line_limit,
                 const std::string &guard_seq);
/// \}

/// \brief Test all left-hand guard sequences to determine whether one of them matches the present
///        text sequence.
///
/// \param tfr       Text file reading abstract
/// \param line_idx  Index at which to seek the putative termination sequence
/// \param pos_idx   Index at which to seek the putative termination sequence
/// \param markers   List of all text guards
/// \{
int applyGuard(const char* line, const size_t pos_idx, const size_t line_limit,
               const std::vector<TextGuard> &markers);

int applyGuard(const TextFileReader &tfr, const int line_idx, const int pos_idx,
               const std::vector<TextGuard> &markers);
/// \}

/// \brief Apply escape sequences to the text.  Escape sequences bind before any other guards such
///        as comments or quotations.
///
/// Overloaded:
///   - Mark all escaped characters in a C-style array of characters
///   - Mark all escaped characters in a Standard Template Library string
///   - Mark all escaped characters (not escape sequences, but escaped characters) in a TextFile
///
/// \param tfr      Text file reading abstract
/// \param textstr  The string (or, an array of characters) to anaylze
/// \param n_char   Trusted length of textstr, if provided as a C-style array
/// \param escapes  Vector of escape sequences
/// \{
std::vector<bool> markEscapedCharacters(const char* textstr, size_t n_char,
                                        const std::vector<TextGuard> &escapes);

std::vector<bool> markEscapedCharacters(const std::string &textstr,
                                        const std::vector<TextGuard> &escapes);

std::vector<bool> markEscapedCharacters(const TextFileReader &tfr,
                                        const std::vector<TextGuard> &escapes);
/// \}

/// \brief Scan a text file and mark off text that is in quotations or commented out.
///
/// Overloaded:
///   - Provide a C-style array of characters (carriage returns will count as line breaks)
///   - Provide a Standard Template Library string (carriage returns in the string will count as
///     line breaks)
///   - Provide a TextFile object (line breaks marked in this object will act as barrier between
///     interpretations of some comment blocks)
///
/// \param tfr           Read-only pointers to the text file of interest once it is read into RAM
/// \param markers       List of markers denoting guards for the text of interest
/// \param alternatives  List of markers denoting guards for other text, binding at the same level
///                      as markers
/// \{
std::vector<bool> markGuardedText(const char* textstr, size_t n_char, const size_t* line_limits,
                                  int line_count, const std::vector<TextGuard> &markers,
                                  const std::vector<TextGuard> &alternatives,
                                  const std::vector<TextGuard> &escapes = { TextGuard("\\") });

std::vector<bool> markGuardedText(const char* textstr, const size_t n_char,
                                  const std::vector<TextGuard> &markers,
                                  const std::vector<TextGuard> &alternatives,
                                  const std::vector<TextGuard> &escapes = { TextGuard("\\") });

std::vector<bool> markGuardedText(const std::string &str,
                                  const std::vector<TextGuard> &markers,
                                  const std::vector<TextGuard> &alternatives,
                                  const std::vector<TextGuard> &escapes = { TextGuard("\\") });

std::vector<bool> markGuardedText(const TextFileReader &tfr,
                                  const std::vector<TextGuard> &markers,
                                  const std::vector<TextGuard> &alternatives,
                                  const std::vector<TextGuard> &escapes = { TextGuard("\\") });
/// \}

/// \brief Bearing some similarities to markGuardedText, this function operates on strings only,
///        marking each guarded segment as one level "deeper" than the one enclosing it.  The scope
///        depth starts at zero, incrementing with each opening guard encountered and decrementing
///        when it reaches a relevant closing guard for the most recent opening guard.
///
/// Overloaded:
///   - Determine scopes in a C-style string, with the option to specify length
///   - Determine scopes in a C++ string
///
/// \param input_text    The text to parse for scope depth
/// \param scope_guards  Scope definitions, all must have closing guards.  Line limits are
///                      irrelevant as the sbstrate is merely a string.  The input text can include
///                      Quoted multi-line strings from input files, having been 
/// \{
std::vector<int> resolveScopes(const char* input_text, int length,
                               const std::vector<TextGuard> scope_guards,
                               ExceptionResponse policy = ExceptionResponse::DIE);

std::vector<int> resolveScopes(const std::string &input_text,
                               const std::vector<TextGuard> scope_guards,
                               ExceptionResponse policy = ExceptionResponse::DIE);
/// \}

/// \brief A quick diagnostic to assess whether a string contains characters that might separate or
///        otherwise annotate text in special ways.  Delimiters masked by quotations or comments
///        will not be counted.
///
/// Overloaded:
///   - Provide the text as a Standard Template Library string and the list of delimiters as a
///     Standard Template Library vector of chars, with a trusted starting position and total
///     length of text to analyze 
///   - Provide the text and list of delimiters as C-style char arrays with trusted lengths
///
/// \param text          The text to analyze
/// \param start_pos     Starting position of the text to analyze
/// \param length        Trusted number of characters to search, if text is provided as a C-style
///                      character array
/// \param delimiters    A list of delimiters (one character each)
/// \param ndelim        Trusted length of delimiters, if provided as a C-style character array
/// \param quote_mask    Indicate whether each character is part of a quotation
/// \param comment_mask  Indicate whether each character is part of a comment
/// \{
int countDelimiters(const char* text, size_t start_pos, size_t length, const char* delimiters,
                    int ndelim, const std::vector<bool> &quote_mask = {},
                    const std::vector<bool> &comment_mask = {});

int countDelimiters(const std::string &text, const std::vector<char> &delimiters,
                    const std::vector<bool> &quote_mask = {},
                    const std::vector<bool> &comment_mask = {});
/// \}

/// \brief Count the number of words delimited by white space in an array of characters.
///
/// Overloaded:
///   - Accept a C-style array of char data
///   - Accept a Standard Template Library string object
///   - Accept a TextFile object with a line, starting position, and range
///
/// \param line       The line data as a std::string object, a C-style array of characters, or the
///                   number of the line in a TextFile object
/// \param tf         Text file object committed to RAM with the data of interest
/// \param start_pos  Position of the first character on the line to read
/// \param length     Length of the character string to analyze
/// \{
int countWords(const char* line, const int start_pos, const int length);
int countWords(const std::string &line, const int start_pos, const int length);
int countWords(const TextFile &tf, const int line, const int start_pos, const int length);
/// \}

/// \brief Add space around key characters in a segment of text and convert selected characters to
///        white space.  The result is computed out-of-place and returned.
///
/// Overloaded:
///   - Accept a C-style array of char data with a trusted starting position and a length of
///     characters to process
///   - Accept a Standard Template Library string and operate on its entirety 
///   - Return the result as a new Standard Template Library string
///   - Return the result in a pre-allocated string or character array
///
/// \param text
/// \param start_pos
/// \param length
/// \param separators
/// \param n_separators
/// \param delimiters
/// \param n_delimiters
/// \param output
/// \param quote_mask
/// \param comment_mask
/// \param output_quoted
/// \param output_commented
/// \{
void digestText(const char* text, size_t start_pos, size_t length, const char* separators,
                int n_separators, const char* delimiters, int n_delimiters, char* output,
                const std::vector<bool> &quote_mask = {},
                const std::vector<bool> &comment_mask = {},
                std::vector<bool> *output_quoted = nullptr,
                std::vector<bool> *output_commented = nullptr);

std::string digestText(const char* text, size_t start_pos, size_t length, const char* separators,
                       int n_separators, const char* delimiters, int n_delimiters,
                       const std::vector<bool> &quote_mask = {},
                       const std::vector<bool> &comment_mask = {},
                       std::vector<bool> *output_quoted = nullptr,
                       std::vector<bool> *output_commented = nullptr);

void digestText(const std::string &text, const std::vector<char> &separators,
                const std::vector<char> &delimiters, std::string *output,
                const std::vector<bool> &quote_mask = {},
                const std::vector<bool> &comment_mask = {},
                std::vector<bool> *output_quoted = nullptr,
                std::vector<bool> *output_commented = nullptr);

std::string digestText(const std::string &text, const std::vector<char> &separators,
                       const std::vector<char> &delimiters,
                       const std::vector<bool> &quote_mask = {},
                       const std::vector<bool> &comment_mask = {},
                       std::vector<bool> *output_quoted = nullptr,
                       std::vector<bool> *output_commented = nullptr);
/// \}
  
/// \brief Separate a sequece of ascii characters into words, using whitespace as the delimiter.
///
/// Overloaded:
///   - Operate on a C-string
///   - Operate on a C++ string
///   - Operate on a C++ string with custom delimiters and no quotations or comments
///   - Operate on a TextFile object
///
/// \param text            The text to parse
/// \param n_char          Number of characters in the text
/// \param comment_mask    A mask of booleans indicating whether text is protected by comment
///                        symbols (optional).  Commented text will be ignored and count as a
///                        delimiter when parsing text into strings.
/// \param quotation_mask  A mask of booleans indicating whether text is protected by quotation
///                        marks (optional).  Quoter text will count as a single string in the
///                        final output.
/// \param quote_marks     A list of tex tguards that represent quotes.  All are given as SINGLE-
///                        line guards as the string is interpreted as only a single line (even if
///                        it contains carriage returns--the multi-line detail pertains to TextFile
///                        objects, where lines are separated by an auxiliary array of limits and
///                        there are no carrieage returns).  The quote marks are needed in order to
///                        shave quote characters off of any substrings that come from stretches of
///                        quoted text.
/// \param delimiters      List of extra delimiters (common choices would include "=" and ",")
/// \{
std::vector<std::string>
separateText(const char* text, int n_char, const std::vector<bool> &comment_mask = {},
             const std::vector<bool> &quotation_mask = {},
             const std::vector<TextGuard> &quote_marks = { TextGuard("\"", "\""),
                                                           TextGuard("'", "'") },
             const std::vector<std::string> &delimiters = {},
             const std::vector<TextGuard> &escapes = { TextGuard("\\") });

std::vector<std::string>
separateText(const std::string &text, const std::vector<bool> &comment_mask = {},
             const std::vector<bool> &quotation_mask = {},
             const std::vector<TextGuard> &quote_marks = { TextGuard("\"", "\""),
                                                           TextGuard("'", "'") },
             const std::vector<std::string> &delimiters = {},
             const std::vector<TextGuard> &escapes = { TextGuard("\\") });

std::vector<std::string>
separateText(const std::string &text, const std::vector<std::string> &delimiters,
             const std::vector<TextGuard> &escapes = { TextGuard("\\") });

std::vector<std::string>
separateText(const TextFile &text, const std::vector<bool> &comment_mask = {},
             const std::vector<bool> &quotation_mask = {},
             const std::vector<TextGuard> &quote_marks = { TextGuard("\"", "\""),
                                                           TextGuard("'", "'") },
             const std::vector<std::string> &delimiters = {},
             const std::vector<TextGuard> &escapes = { TextGuard("\\") });
/// \}

/// \brief Find one string within a vector of other strings.  Return the index, or the length of
///        the input vector (akind to std::find) if the index is not found.
///
/// \param vec  The vector of strings to search
/// \param str  The string to search for
int findStringInVector(const std::vector<std::string> &vec, const std::string &query);

/// \brief Convert a vector of strings into a vector of integers, checking each for validity and
///        reporting the result or substituting a zero if it looks suspicious.
///
/// \param sv      Vector of strings to convert
/// \param policy  The course of action if a malformed number is found
std::vector<int> vectorStrtol(const std::vector<std::string> &sv,
                              ExceptionResponse policy = ExceptionResponse::DIE);

/// \brief Convert a vector of strings into a vector of double-precision reals, checking each for
///        validity and reporting the result or substituting a zero if it looks suspicious.
///
/// \param sv      Vector of strings to convert
/// \param policy  The course of action if a malformed number is found
std::vector<double> vectorStrtod(const std::vector<std::string> &sv,
                                 ExceptionResponse policy = ExceptionResponse::DIE);

} // namespace parse
} // namespace stormm

#include "parse.tpp"

#endif
