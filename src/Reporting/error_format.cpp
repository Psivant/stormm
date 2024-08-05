#include <sys/ioctl.h>
#include <unistd.h>
#include <cstring>
#include <vector>
#include <stdexcept>
#include "copyright.h"
#include "error_format.h"
#include "FileManagement/file_listing.h"

namespace stormm {
namespace errors {

//-------------------------------------------------------------------------------------------------
std::string terminalFormat(const std::string &message, const char* class_caller,
                           const char* method_caller, const int implicit_indent,
                           const int first_indent, const int subsq_indent, const int width_in,
			   const RTMessageKind style) {

  using diskutil::OperatingSystem;
  using diskutil::detected_os;

  // Form the starting string
  std::string msg;
  if (class_caller != nullptr &&
      ((strlen(class_caller) >= 10 && strncmp(class_caller, " Warning: ", 10) == 0) ||
       (strlen(class_caller) >= 8 && strncmp(class_caller, " Alert: ", 8) == 0))) {
    msg = std::string(class_caller);
    if (method_caller != nullptr) {
      msg += " :: (" + std::string(method_caller) + ") ";
    }
    else {
      msg += " :: ";
    }
  }
  else if (class_caller != nullptr && strlen(class_caller) > 0 && method_caller != nullptr &&
           strlen(method_caller) > 0) {
    msg = std::string(class_caller) + " :: (" + std::string(method_caller) + ") ";
  }
  else if (class_caller != nullptr && strlen(class_caller) > 0) {
    msg = std::string(class_caller) + " :: ";
  }
  else if (method_caller != nullptr && strlen(method_caller) > 0) {
    msg = std::string(method_caller) + " :: ";
  }
  msg += message;

  // Obtain the console size
  struct winsize console_dims;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &console_dims);

  // Parse the message word by word
  const int msize = msg.size();
  std::vector<std::string> words;
  std::vector<int> white_space;
  int i = 0;
  if (msize > 0 && msg[0] == ' ') {
    words.push_back(std::string(""));
  }
  while (i < msize) {

    // Extract a new word, exclusive of white space and carriage returns
    std::string new_word;
    int nletter = 0;
    while (i < msize && msg[i] != ' ' && msg[i] != '\n') {
      new_word += msg[i];
      nletter++;
      i++;
    }
    if (nletter > 0) {
      words.push_back(new_word);
    }

    // Test whether a carriage return was detected
    if (i < msize && msg[i] == '\n') {
      if (nletter > 0) {
        white_space.push_back(0);
      }
      words.push_back("\n");
      i++;
    }

    // Look for a stretch of white space
    int nspc = 0;
    while (i < msize && msg[i] == ' ') {
      nspc++;
      i++;
    }
    white_space.push_back(nspc);
  }

  // The std::runtime error printout in Linux begins with "  what():  ", which is eleven characters
  // of indentation.  On a Macintosh, it is such a long line as to make further printing not
  // worthwhile (just make a new line).  Take all of that into account when crafting a message to
  // keep text fom wrapping and arranged in a clean block on any given OS.
  const int word_count = words.size();
  const int width = (width_in > 0) ? width_in : (console_dims.ws_col > 2) ?
                                                console_dims.ws_col - 1 : 80;
  std::string parsed_msg;
  switch (style) {
  case RTMessageKind::ERROR:
    switch (detected_os) {
    case OperatingSystem::LINUX:
      for (int i = 0; i < first_indent; i++) {
        parsed_msg += " ";
      }
    case OperatingSystem::UNIX:
    case OperatingSystem::WINDOWS:
      break;
    case OperatingSystem::MAC_OS:
      for (int i = 0; i < implicit_indent + first_indent; i++) {
        parsed_msg += " ";
      }
      parsed_msg += "\n";
      break;
    }
    break;
  case RTMessageKind::TABULAR:
    break;
  }
  std::string subsq_spacer;
  for (int i = 0; i < subsq_indent; i++) {
    subsq_spacer += " ";
  }
  int nchar = implicit_indent + first_indent;
  int line_count = 0;
  for (int i = 0; i < word_count; i++) {
    const int word_length = words[i].size();

    // If this is the baseline indentation, always add the next word.
    // If the next word will fit, add it to the current line.
    if ((line_count == 0 && nchar == implicit_indent + first_indent) ||
        (line_count > 0 && nchar == subsq_indent) || nchar + word_length < width) {

      // When adding the word, check for carriage returns
      parsed_msg += words[i];
      if (words[i].size() == 1 && words[i][0] == '\n') {
        parsed_msg += subsq_spacer;
        nchar = subsq_indent;
      }
      else {
        nchar += word_length;
      }
    }

    // Test the spacing and the next word.  Add spacing if the next word will fit on the line
    // after it.  Otherwise start a new line and assume that takes care of all such spacing.
    if (i + 1 < word_count) {
      if (nchar + white_space[i] + static_cast<int>(words[i + 1].size()) < width) {
        for (int j = 0; j < white_space[i]; j++) {
          parsed_msg += " ";
        }
        nchar += white_space[i];
      }
      else {
        parsed_msg += "\n" + subsq_spacer;
        nchar = subsq_indent;
      }
    }
  }

  // Pad the result with white space if there was a preset width
  if (width_in > 0) {
    for (int i = nchar; i < width; i++) {
      parsed_msg += " ";
    }
  }

  return parsed_msg;
}

//-------------------------------------------------------------------------------------------------
void rtErr(const std::string &message, const char* class_caller, const char* method_caller) {
  const std::string parsed_msg = terminalFormat(message, class_caller, method_caller, 11, 0, 11);
  throw std::runtime_error(parsed_msg);
}

//-------------------------------------------------------------------------------------------------
void rtWarn(const std::string &message, const char* class_caller, const char* method_caller) {
  std::string ccall(" Warning: ");
  if (class_caller != nullptr) {
    ccall += std::string(class_caller);
  }
  const std::string parsed_msg = terminalFormat(message, ccall.c_str(), method_caller, 0, 0, 10);
  printf("%s\n", parsed_msg.c_str());
}

//-------------------------------------------------------------------------------------------------
void rtAlert(const std::string &message, const char* class_caller, const char* method_caller) {
  std::string ccall(" Alert: ");
  if (class_caller != nullptr) {
    ccall += std::string(class_caller);
  }
  const std::string parsed_msg = terminalFormat(message, ccall.c_str(), method_caller, 0, 0, 8);
  printf("%s\n", parsed_msg.c_str());
}

//-------------------------------------------------------------------------------------------------
std::string listSeparator(const int current_item, const int item_count) {
  if (item_count > 2) {
    if (current_item < item_count - 2) {
      return ", ";
    }
    else {
      return ", and ";
    }
  }
  else if (current_item == 2) {
    return ", ";
  }
  else {
    return "";
  }
  __builtin_unreachable();
}
  
} // namespace errors
} // namespace stormm
