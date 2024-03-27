// -*-c++-*-
#ifndef STORMM_TEXTFILE_H
#define STORMM_TEXTFILE_H

#include <string>
#include <vector>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "FileManagement/file_enumerators.h"
#include "Reporting/reporting_enumerators.h"
#include "parsing_enumerators.h"

namespace stormm {
namespace parse {

using diskutil::PrintSituation;
using review::TextEnds;

/// \brief Abstract for the TextFile object, providing read-only access
struct TextFileReader {

  /// \brief The constructor takes length constants and constant pointers.
  TextFileReader(int line_count_in, const int* line_lengths_in, const size_t* line_limits_in,
                 const char* text_in, const std::string file_name_in);

  /// \brief Take the default copy and move constructors.  The assignment operators will get
  ///        implicitly deleted as this is just a collection of constants.
  /// \{
  TextFileReader(const TextFileReader &original) = default;
  TextFileReader(TextFileReader &&original) = default;
  /// \}

  const int line_count;
  const int* line_lengths;
  const size_t* line_limits;
  const char* text;
  const std::string file_name;
};
  
/// \brief Structure for translating a text file into a compact, rapidly parsable vector of
///        characters in CPU RAM.  The struct contains two nested struct definitions, for a Reader
///        and a Writer.  Private objects of these structs can then be accessed with eponymous
///        data() getter functions, returning the appropriate kind of access depending on the
///        const-ness of the TextFile object itself.
class TextFile {
public:

  /// \brief Constructor for taking an ascii file or a very long, formatted string and transforming
  ///        it into a std::vector of characters with line limits recorded.
  ///
  /// Overloaded:
  ///   - Construct an empty object after taking no arguments
  ///   - Construct a complete object from a named file (throws an exception if the file does not
  ///     exist)
  ///
  /// \param file_name   Name of the input file
  /// \param source      Origin of the text--disk or RAM
  /// \param content     Content for the TextFile, and perhaps later an ASCII text file, to hold
  /// \param caller      (Optional) name of the calling function
  /// \{
  TextFile();

  TextFile(const std::string &file_name, TextOrigin source = TextOrigin::DISK,
           const std::string &content = std::string(""),
           const std::string &caller = std::string(""));

  TextFile(const char* content, const size_t length, const std::string &caller = std::string(""));
  /// \}

  /// \brief The default copy and move constructors as well as assignment operators are acceptable.
  /// \{
  TextFile(const TextFile &original) = default;
  TextFile(TextFile &&original) = default;
  TextFile& operator=(const TextFile &other) = default;
  TextFile& operator=(TextFile &&other) = default;
  /// \}
  
  /// \brief Default destructor
  ~TextFile() = default;

  /// \brief Get the name of the original file.
  std::string getFileName() const;

  /// \brief Get the line count of a text file after converting it to a character vector in memory.
  int getLineCount() const;

  /// \brief Get one of the line limits of a text file converted to a character vector in memory.
  ///
  /// \param index   The line index.  The array of limits contains one more indices than the text
  ///                file has lines, to allow the end of the last line to be determined.
  int getLineLimits(int index) const;

  /// \brief Get the length of a line.  A bounds check will be applied to the line index number.
  ///
  /// \param index   The line index.  The array of limits contains one more indices than the text
  ///                file has lines, to allow the end of the last line to be determined.
  int getLineLength(int index) const;

  /// \brief Get the length of the longest line in the file.
  int getLongestLineLength() const;

  /// \brief Get the number of text characters in the object's buffer (this will not count implicit
  ///        carriage returns between lines).
  size_t getTextSize() const;
  
  /// \brief Get one character of a text file after converting it to a character vector in memory.
  ///
  /// \param index   The character index, as ascertained by line limits and some offset
  char getChar(int index) const;

  /// \brief Get a line from the file as a string.
  ///
  /// \param line_index  The line of interest, indexing starting from zero
  std::string getLineAsString(int line_index) const;
  
  /// \brief Get a char pointer to a specific index in the object.
  ///
  /// \param index   The character index, as ascertained by line limits and some offset
  const char* getTextPointer(int index) const;

  /// \brief Get a char pointer to a specific line in the object.
  ///
  /// \param line_index  The line of interest, indexing starting from zero
  const char* getLinePointer(int line_index) const;
  
  /// \brief Get an abstract of a text file's CPU-RAM representation, for ease of use.
  const TextFileReader data() const;

  /// \brief Extract a string based on a line number, starting position, and length.
  ///
  /// \param line_number    The line number where the text shall be found
  /// \param start_pos      Starting position on the line to begin reading (defaults to the start
  ///                       of the line)
  /// \param string_length  Length of the string to extract.  A value of -1 indicates that reading
  ///                       shall continue until the end of the line.
  std::string extractString(int line_number, int start_pos = 0, int string_length = -1) const;

  /// \brief Extract a char4 tuple based on a line number, starting position, and length.
  ///
  /// \param line_number    The line number where the text shall be found
  /// \param start_pos      Starting position on the line to begin reading (defaults to the start
  ///                       of the line)
  /// \param string_length  Length of the string to extract.  The default value of 4 will try to
  ///                       read four characters, but shorter values will leave blank space at the
  ///                       end of the tuple.
  char4 extractChar4(int line_number, int start_pos, int string_length) const;

  /// \brief Convert all content to a string, with the option of making line endings carriage
  ///        returns or simple spaces (fused line endings, which could create combined words out
  ///        of the last word on one line and the first word on another, are not accepted).
  ///
  /// \param line_endings  Specify whether to insert newlines or white space between lines
  std::string toString(TextEnds line_endings = TextEnds::NEWLINE) const;

  /// \brief Write the contents of the object to disk.
  ///
  /// Overloaded:
  ///   - Supply an open file pointer
  ///   - Supply a new file name, distinct from the original name stored in the object, and a file
  ///     status expectation
  ///   - Supply just a file status expectation and (overwrite, or perhaps append in some clever
  ///     implementation) the original file
  ///
  /// \param foutp         Open output stream
  /// \param new_filename  Name of a distinct file to write
  /// \param expectation   Expected state of the output file
  /// \{
  void write(std::ofstream *foutp) const;
  void write(const std::string &new_filename,
             PrintSituation expectation = PrintSituation::OPEN_NEW) const;
  void write(PrintSituation expectation = PrintSituation::OPEN_NEW) const;
  /// \}
  
private:

  /// Name of the file that was read
  std::string orig_file;

  /// The number of lines detected in the file
  int line_count;

  /// Lengths for each line in the concatenated character array
  std::vector<int> line_lengths;
  
  /// Limits for each line's text in the concatenated character array
  std::vector<size_t> line_limits;
  
  /// The text, sans carriage returns (see line limits to determine their locations)
  std::vector<char> text;

  /// \brief Validate a requested line index.  Raise an error, along with the calling function name
  ///        for tracing purposes, if the index is not valid.
  ///
  /// \param index   The line of interest (indexing starts at zero)
  /// \param caller  Name of the calling function
  void validateLineIndex(int index, const char* caller) const;
  
  /// \brief Set the name of the corresponding file based on the characteristics of several
  ///        constructor input variables.
  ///
  /// \param file_name   Name of the input file (if the source is DISK--otherwise if the source
  ///                    is RAM and the content is blank, then this is assumed to be content and
  ///                    a blank corresponding file name is returned)
  /// \param source      Origin of the text--disk or RAM
  /// \param content     Content for the TextFile, and perhaps later an ASCII text file, to hold
  std::string setFileName(const std::string &file_name, TextOrigin source,
                          const std::string &content);

  /// \brief Break a large, formatted string into separate lines based on carriage returns.
  ///
  /// Overloaded:
  ///   - Accept a C-style character array with trusted length
  ///   - Accept a Standard Template Library string
  ///
  /// \param text_in  The text to parse.  See setFileName and above in the TextFile documentation
  ///                 to understand how this could be either of two constructor inputs.
  /// \param n_char   Trusted number of characters in the input text (if provided as a C-style
  ///                 array)
  /// \{
  void linesFromString(const char* text_in, const size_t n_char);
  void linesFromString(const std::string &text_in);
  /// \}
  
  /// \brief Check that a certain number of characters can be extracted from one line at the
  ///        specified positions.
  ///
  /// \param line_number    Number of the line to access
  /// \param start_pos      Starting position on the line to begin reading (defaults to the start
  ///                       of the line)
  /// \param string_length  Length of the character string to read (default -1 implies reading to
  ///                       the end of the line)
  int checkAvailableLength(int line_number, int start_pos = 0, int string_length = -1) const;
};

} // namespace parse
} // namespace stormm

#endif
