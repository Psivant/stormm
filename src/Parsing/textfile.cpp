#include <iostream>
#include <fstream>
#include "copyright.h"
#include "textfile.h"
#include "FileManagement/file_util.h"
#include "Reporting/error_format.h"

namespace stormm {
namespace parse {

using diskutil::openOutputFile;

//-------------------------------------------------------------------------------------------------
TextFileReader::TextFileReader(const int line_count_in, const int* line_lengths_in,
                               const size_t* line_limits_in, const char* text_in,
                               const std::string file_name_in) :
  line_count{line_count_in}, line_lengths{line_lengths_in}, line_limits{line_limits_in},
  text{text_in}, file_name{file_name_in}
{}

//-------------------------------------------------------------------------------------------------
TextFile::TextFile() :
    orig_file{std::string("")},
    line_count{0},
    line_lengths{},
    line_limits{std::vector<size_t>(1, 0)},
    text{std::vector<char>()}
{}

//-------------------------------------------------------------------------------------------------
TextFile::TextFile(const std::string &file_name, const TextOrigin source,
                   const std::string &content, const std::string &caller) :
  orig_file{setFileName(file_name, source, content)},
  line_count{0},
  line_lengths{},
  line_limits{},
  text{}
{
  switch (source) {
  case TextOrigin::DISK:
    {
      // Start reading from a file.  The content variable will be ignored.
      std::ifstream finp;
      finp.open(file_name.c_str());
      if (finp.is_open() == false) {
        if (caller.size() == 0) {
          rtErr(file_name + " was not found.", "TextFile");
        }
        else {
          rtErr(file_name + " was not found when called from " + caller + ".", "TextFile");
        }
      }
      std::string line;
      int total_chars = 0;
      int line_counter = 0;
      text.resize(0);
      line_limits.resize(0);
      line_limits.push_back(total_chars);
      while (std::getline(finp, line)) {
        const int line_length = line.size();
        for (int i = 0; i < line_length; i++) {
          text.push_back(line[i]);
        }
        total_chars += line_length;
        line_limits.push_back(total_chars);
        line_counter++;
      }
      line_count = line_counter;
      line_lengths.resize(line_count);
      for (int i = 0; i < line_count; i++) {
        line_lengths[i] = line_limits[i + 1] - line_limits[i];
      }

      // Close input file
      finp.close();
    }
    break;
  case TextOrigin::RAM:
    {
      // If ther content variable is blank, this implies that the file_name variable is the
      // actual content.  Either way, scan the content for line breaks (carriage returns) and
      // assign a file name as appropriate.
      if (content.size() > 0) {
        linesFromString(content);
      }
      else {
        linesFromString(file_name);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
TextFile::TextFile(const char* content, const size_t length, const std::string &caller) :
  orig_file{"None"},
  line_count{0},
  line_lengths{},
  line_limits{},
  text{}
{
  linesFromString(content, length);
}

//-------------------------------------------------------------------------------------------------
std::string TextFile::getFileName() const {
  if (orig_file.size() > 0LLU) {
    return orig_file;
  }
  else {
    
    // If there is no file name, assume that the text was assembled in RAM.  A file name will need
    // to be assigned for writing purposes, but to the rest of the program the file will report
    // this as its origin.
    return std::string("[RAM-assembled text]");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int TextFile::getLineCount() const {
  return line_count;
}

//-------------------------------------------------------------------------------------------------
int TextFile::getLineLimits(const int index) const {
  if (index < 0 || index > line_count) {
    rtErr("Text file " + orig_file + " has " + std::to_string(line_count) + " lines, " +
          std::to_string(index) + " requested.", "TextFile", "getLineLimits");
  }
  return line_limits[index];
}

//-------------------------------------------------------------------------------------------------
int TextFile::getLineLength(const int index) const {
  validateLineIndex(index, "getLineLength");
  return line_limits[index + 1] - line_limits[index];
}

//-------------------------------------------------------------------------------------------------
int TextFile::getLongestLineLength() const {
  size_t result = 0;
  for (int i = 0; i < line_count; i++) {
    result = std::max(result, line_limits[i + 1] - line_limits[i]);
  }
  return static_cast<int>(result);
}

//-------------------------------------------------------------------------------------------------
size_t TextFile::getTextSize() const {
  return line_limits[line_count];
}
  
//-------------------------------------------------------------------------------------------------
char TextFile::getChar(const int index) const {
  return text[index];
}

//-------------------------------------------------------------------------------------------------
std::string TextFile::getLineAsString(const int line_index) const {
  validateLineIndex(line_index, "getLineAsString");
  std::string result;
  const int llim = line_limits[line_index];
  const int hlim = line_limits[line_index + 1];
  result.reserve(hlim - llim);
  for (int i = llim; i < hlim; i++) {
    result += text[i];
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
const char* TextFile::getTextPointer(const int index) const {
  return &text[index];
}

//-------------------------------------------------------------------------------------------------
const char* TextFile::getLinePointer(const int line_index) const {
  validateLineIndex(line_index, "getLinePointer");
  return &text[line_limits[line_index]];
}

//-------------------------------------------------------------------------------------------------
const TextFileReader TextFile::data() const {
  return TextFileReader(line_count, line_lengths.data(), line_limits.data(), text.data(),
                        orig_file);
}

//-------------------------------------------------------------------------------------------------
std::string TextFile::extractString(const int line_number, const int start_pos,
                                    const int string_length) const {
  const int actual_length = checkAvailableLength(line_number, start_pos, string_length);
  std::string buffer;
  buffer.resize(actual_length);
  const char* text_ptr = &text[line_limits[line_number] + start_pos];
  for (int i = 0; i < actual_length; i++) {
    buffer[i] = text_ptr[i];
  }
  return buffer;
}

//-------------------------------------------------------------------------------------------------
char4 TextFile::extractChar4(const int line_number, const int start_pos,
                             const int string_length) const {
  const int actual_length = checkAvailableLength(line_number, start_pos, string_length);
  const char* text_ptr = &text[line_limits[line_number] + start_pos];
  return { (actual_length > 0) ? text_ptr[0] : ' ', (actual_length > 1) ? text_ptr[1] : ' ',
           (actual_length > 2) ? text_ptr[2] : ' ', (actual_length > 3) ? text_ptr[3] : ' ' };
}

//-------------------------------------------------------------------------------------------------
void TextFile::validateLineIndex(const int index, const char* caller) const {
  if (index < 0 || index >= line_count) {
    rtErr("Text file " + orig_file + " has " + std::to_string(line_count) + " lines, " +
          std::to_string(index) + " requested.", "TextFile", caller);
  }
}

//-------------------------------------------------------------------------------------------------
std::string TextFile::setFileName(const std::string &file_name, const TextOrigin source,
                                  const std::string &content) {
  switch (source) {
  case TextOrigin::DISK:
    return file_name;
  case TextOrigin::RAM:
    return (content.size() > 0) ? file_name : std::string("");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void TextFile::linesFromString(const char* text_in, const size_t n_char) {
  int n_br = 0;
  for (size_t i = 0; i < n_char; i++) {
    n_br += (text_in[i] == '\n');
  }
  line_count = n_br + (text_in[n_char - 1] != '\n');
  line_lengths.resize(line_count);
  line_limits.resize(line_count + 1);
  text.resize(n_char - n_br + 1);
  size_t n_tx = 0;
  n_br = 0;
  for (size_t i = 0; i < n_char; i++) {
    text[n_tx] = text_in[i];
    if (text_in[i] == '\n') {
      n_br++;
      line_limits[n_br] = n_tx;
    }
    else {
      n_tx++;
    }
  }

  // Determine line lengths.  Set the upper line limit index to the total number of characters
  // in case the string does not end in a carriage return.
  line_limits[line_count] = n_tx;
  for (int i = 0; i < line_count; i++) {
    line_lengths[i] = line_limits[i + 1] - line_limits[i];
  }

  // Ensure that the final line limit is catalogged
  line_limits[line_count] = n_tx;
}

//-------------------------------------------------------------------------------------------------
void TextFile::linesFromString(const std::string &text_in) {
  linesFromString(text_in.data(), text_in.size());
}

//-------------------------------------------------------------------------------------------------
int TextFile::checkAvailableLength(const int line_number, const int start_pos,
                                   const int string_length) const {
  validateLineIndex(line_number, "checkAvailableLength");
  const int available_chars = line_limits[line_number + 1] - line_limits[line_number];
  if (start_pos > available_chars || start_pos + string_length > available_chars) {
    rtErr("Line " + std::to_string(line_number) + " of text file " + orig_file + " has " +
          std::to_string(available_chars) + " characters (requested starting position " +
          std::to_string(start_pos) + ", length " + std::to_string(string_length) + ").",
          "TextFile", "checkAvailableLength");
  }
  return ((string_length < 0) ? available_chars - start_pos : string_length);
}

//-------------------------------------------------------------------------------------------------
std::string TextFile::toString(const TextEnds line_endings) const {
  std::string buffer;
  switch (line_endings) {
  case TextEnds::AS_IS:
    buffer.resize(line_limits[line_count] +
                  static_cast<size_t>((line_count > 0) ? line_count - 1 : 0));
    break;
  case TextEnds::NEWLINE:
    buffer.resize(line_limits[line_count] + static_cast<size_t>(line_count));
    break;
  }
  size_t counter = 0;
  for (int i = 0; i < line_count; i++) {
    for (size_t j = line_limits[i]; j < line_limits[i + 1]; j++) {
      buffer[counter] = text[j];
      counter++;
    }
    switch (line_endings) {
    case TextEnds::AS_IS:
      if (i < line_count - 1) {
        buffer[counter] = ' ';
      }
      break;
    case TextEnds::NEWLINE:
      buffer[counter] = '\n';
      break;
    }
    counter++;
  }
  buffer.resize(counter);
  return buffer;
}

//-------------------------------------------------------------------------------------------------
void TextFile::write(std::ofstream *foutp) const {
  for (int i = 0; i < line_count; i++) {
    std::string buffer(line_lengths[i] + 1, '\n');
    const size_t hlim = line_limits[i + 1];
    size_t k = 0;
    for (size_t j = line_limits[i]; j < hlim; j++) {
      buffer[k] = text[j];
      k++;
    }
    foutp->write(buffer.c_str(), line_lengths[i] + 1);
  }
}

//-------------------------------------------------------------------------------------------------
void TextFile::write(const std::string &new_filename, const PrintSituation expectation) const {
  if (new_filename.size() == 0LLU && orig_file.size() == 0LLU) {
    rtErr("No file name is available.", "TextFile", "write");
  }
  const std::string& fname = (new_filename.size() > 0LLU) ? new_filename : orig_file;
  std::ofstream foutp = openOutputFile(fname, expectation, "Open an ASCII file for writing a "
                                       "TextFile object to disk.");
  write(&foutp);
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
void TextFile::write(const PrintSituation expectation) const {
  write(orig_file, expectation);
}

} // namespace parse
} // namespace stormm
