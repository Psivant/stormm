#include <cstring>
#include <sys/ioctl.h>
#include <unistd.h>
#include "copyright.h"
#include "FileManagement/file_util.h"
#include "display.h"
#include "error_format.h"
#include "summary_file.h"

namespace stormm {
namespace review {

using diskutil::openOutputFile;
using display::horizontalRule;
using display::terminalHorizontalRule;
using errors::RTMessageKind;
using errors::terminalFormat;

//-------------------------------------------------------------------------------------------------
int findFormatWidth(const std::ostream *foutp) {
  int recommended_width = default_output_file_width;
  if (foutp->rdbuf() == std::cout.rdbuf()) {
    struct winsize console_dims;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &console_dims);
    recommended_width = (console_dims.ws_col > 2) ? console_dims.ws_col - 1 :
                                                    default_output_file_width;
  }
  return recommended_width;  
}

//-------------------------------------------------------------------------------------------------
std::string stormmSplash(int width) {
  std::string result = horizontalRule("+", "+", width);
  result += "Copyright 2024, PsiVant Therapeutics\n\n";
  result += "Permission is hereby granted, free of charge, to any person obtaining a copy of "
            "this software and associated documentation files (the \"Software\"), to deal in the "
            "Software without restriction, including without limitation the rights to use, copy, "
            "modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, "
            "and to permit persons to whom the Software is furnished to do so, subject to the "
            "following conditions:\n\n"
            "The above copyright notice and this permission notice shall be included in all "
            "copies or substantial portions of the Software.\n\n"
            "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR "
            "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR "
            "A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT "
            "HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF "
            "CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE "
            "OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.\n";
  result += horizontalRule("+", "+", width);
  return result;
}

//-------------------------------------------------------------------------------------------------
void stormmSplash(std::ofstream *foutp, const int width) {
  std::string buffer = stormmSplash(width);
  buffer = indentText(buffer, 0, width);
  foutp->write(buffer.c_str(), buffer.size());
}

//-------------------------------------------------------------------------------------------------
void stormmSplash(std::ostream *foutp, const int width) {
  const int actual_width = (width <= 0) ? findFormatWidth(foutp) : width;
  std::string buffer = stormmSplash(actual_width);
  buffer = indentText(buffer, 0, width);
  foutp->write(buffer.c_str(), buffer.size());
}

//-------------------------------------------------------------------------------------------------
std::string stormmWatermark(int width) {
  std::string result = horizontalRule("+", "+", width, '-');
  result += "This watermark signifies complete writing of the output, as intended, and should be "
            "taken as an indicator of a successful process.\n";
  result += horizontalRule("+", "+", width, '-');
  return result;
}

//-------------------------------------------------------------------------------------------------
void stormmWatermark(std::ofstream *foutp, const int width) {
  std::string buffer = stormmWatermark(width);
  buffer = indentText(buffer, 0, width);
  foutp->write(buffer.data(), buffer.size());  
}
  
//-------------------------------------------------------------------------------------------------
void stormmWatermark(std::ostream *foutp, const int width) {
  const int actual_width = (width <= 0) ? findFormatWidth(foutp) : width;
  std::string buffer = stormmWatermark(actual_width);
  buffer = indentText(buffer, 0, actual_width);
  foutp->write(buffer.data(), buffer.size());  
}
  
//-------------------------------------------------------------------------------------------------
void summaryHeader(const std::string &description, std::ofstream *foutp, const int file_width) {

  // Lead with the STORMM authorship and disclaimer
  stormmSplash(foutp);
  foutp->write("\n", 1);

  // Write the description, whether of a program itself or the program feature in use
  const int recommended_width = (file_width == -1) ? findFormatWidth(foutp) : file_width;
  terminalHorizontalRule("+", "+", recommended_width, '-', foutp);
  const std::string fmt_desc = terminalFormat(description, nullptr, nullptr, 0, 0, 0,
                                              recommended_width, RTMessageKind::ERROR) + "\n";
  foutp->write(fmt_desc.c_str(), fmt_desc.size());
  terminalHorizontalRule("+", "+", recommended_width, '-', foutp);
}

//-------------------------------------------------------------------------------------------------
char commentSymbol(OutputSyntax format) {
  switch (format) {
  case OutputSyntax::MATPLOTLIB:
    return '#';
  case OutputSyntax::MATRIX_PKG:
    return '%';
  case OutputSyntax::STANDALONE:
    return '|';
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
char commentSymbol(GridFileSyntax format) {
  switch (format) {
  case GridFileSyntax::MATPLOTLIB:
    return commentSymbol(OutputSyntax::MATPLOTLIB);
  case GridFileSyntax::MATRIX_PKG:
    return commentSymbol(OutputSyntax::MATRIX_PKG);
  case GridFileSyntax::OPEN_DX:
    return '#';
  case GridFileSyntax::CUBEGEN:
    return '!';
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string indentText(const std::string &text, const std::string &marker, const int indent,
                       const int width, const bool mark_on_all_lines,
                       const TextEnds block_ending) {
  const int actual_width = (width <= 0) ? default_output_file_width : width;
  size_t adv = 0;
  const size_t maxchar = text.size();
  std::string carriage_return("\n");
  carriage_return.resize(1 + indent);
  std::string result;
  result.reserve(maxchar + ((maxchar / actual_width) *
                            (marker.size() + indent + (marker.size() > 0))));
  result.resize(indent);
  for (int i = 0; i < indent; i++) {
    result[i] = ' ';
    carriage_return[i + 1] = ' ';
  }
  int total_indent;
  if (marker.size() > 0) {
    total_indent = indent + marker.size() + 1;
    result += marker + " ";
    if (mark_on_all_lines) {
      carriage_return += marker + " ";
    }
    else {
      carriage_return += std::string(marker.size() + 1, ' ');
    }
  }
  else {
    total_indent = indent;
  }
  int lcon = total_indent;
  bool first_word = true;
  while (adv < maxchar) {
    
    // Advance to the next word.
    while (adv < maxchar && text[adv] == ' ') {
      result += " ";
      lcon++;
      adv++;
    }
    bool bare_line = first_word;
    first_word = false;
    
    // Transfer carriage returns explicitly.
    while (text[adv] == '\n') {

      // Remove trailing white space
      const size_t res_length = result.size();
      size_t last_nonspace = res_length - 1;
      while (last_nonspace < res_length && result[last_nonspace] == ' ') {
        last_nonspace--;
      }
      result.resize(last_nonspace + 1);

      // Add the newline and indent further as appropriate.
      result += carriage_return;
      lcon = total_indent;
      adv++;
      bare_line = true;
      while (text[adv] == ' ') {
        result += " ";
        lcon++;
        adv++;
      }
    }

    // Count the number of characters in the next word
    size_t padv = adv;
    while (padv < maxchar && text[padv] != ' ' && text[padv] != '\n') {
      padv++;
    }

    // Can the next word fit in the space left on the current line of output?  If the output line
    // is just starting, the word must be contributed.
    const int wordlen = padv - adv;
    if (bare_line) {
      result += text.substr(adv, wordlen);
      lcon += wordlen;
    }
    else if (lcon + wordlen <= actual_width) {
      result += text.substr(adv, wordlen);
      lcon += wordlen;
    }
    else {

      // Again, remove trailing white space
      const size_t res_length = result.size();
      size_t last_nonspace = res_length - 1;
      while (last_nonspace < res_length && result[last_nonspace] == ' ') {
        last_nonspace--;
      }
      result.resize(last_nonspace + 1);

      // Make a carriage return and add the next word at the start of the following line.
      result += carriage_return + text.substr(adv, wordlen);
      lcon = total_indent + wordlen;
    }
    adv = padv;
  }

  // Do not print a hanging, commented or indented line.
  if (total_indent > 0) {
    while (static_cast<int>(result.size()) - total_indent - 1 >= 0 &&
           result[result.size() - total_indent - 1] == '\n') {
      result.resize(result.size() - total_indent - 1);
    }
  }
  else {
    while (static_cast<int>(result.size()) - 2 >= 0 && result[result.size() - 1] == '\n' &&
           result[result.size() - 2] == '\n') {
      result.resize(result.size() - 1);
    }    
  }
    
  // Do print a carriage return if appropriate to terminate a nontrivial line.
  switch (block_ending) {
  case TextEnds::AS_IS:
    break;
  case TextEnds::NEWLINE:
    const size_t rsz = result.size();
    size_t rbk = rsz - 1;
    while (rbk < rsz && result[rbk] == ' ') {
      rbk--;
    }
    if (marker.size() > 0) {
      if (rbk < rsz && rbk >= marker.size() - 1 &&
          strncmp(&result[rbk - marker.size() + 1], marker.data(), marker.size()) != 0) {
        result += "\n";
      }
    }
    else {
      if (result[result.size() - 1] != '\n') {
        result.append("\n");
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string indentText(const std::string &text, const int indent, const int width,
                       const TextEnds block_ending) {
  return indentText(text, std::string(""), indent, width, false, block_ending);
}

//-------------------------------------------------------------------------------------------------
std::string protectText(const std::string &text, const char protect, const int width,
                        const TextEnds block_ending) {
  return indentText(text, std::string(1, protect), 0, width, true, block_ending);
}

//-------------------------------------------------------------------------------------------------
void printProtectedText(const std::string &text, const char protect, std::ofstream *foutp,
                        const int file_width, const TextEnds block_ending) {
  const std::string product = protectText(text, protect, file_width, block_ending);
  foutp->write(product.data(), product.size());
}

//-------------------------------------------------------------------------------------------------
void printProtectedText(const std::string &text, const char protect,
                        const std::string &output_file, const int file_width,
                        const PrintSituation expectation, const TextEnds block_ending) {
  std::ofstream foutp = openOutputFile(output_file, expectation, "open file to write a block of "
                                       "protected text");
  printProtectedText(text, protect, &foutp, file_width, block_ending);
  foutp.close();
}

} // namespace review
} // namespace stormm
