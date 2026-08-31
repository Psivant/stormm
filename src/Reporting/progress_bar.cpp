#include <unistd.h>
#include <sys/ioctl.h>
#include "copyright.h"
#include "error_format.h"
#include "summary_file.h"
#include "progress_bar.h"

namespace stormm {
namespace reporting {

using review::default_output_file_width;
  
//-------------------------------------------------------------------------------------------------
ProgressBar::ProgressBar(const std::string &title_in, const std::string &open_bracket_in,
                         const std::string &close_bracket_in, const char finished_mark_in,
                         const char todo_mark_in, const int cycle_count_in,
                         const ProgBarStyle style_in, std::ostream &output_in) :
    progress{0}, cycle_count{cycle_count_in}, last_percent{0},
    terminal_width{default_output_file_width},
    title_width{static_cast<int>(title_in.size())}, bar_width{0},
    style{style_in},
    update_called{false}, title_within_bar{false},
    title{title_in},
    finished_mark{finished_mark_in},
    todo_mark{todo_mark_in}, 
    open_bracket{open_bracket_in},
    close_bracket{close_bracket_in},
    output{&output_in},
    bar_contents{}
{
  updateTerminalWidth();
}

//-------------------------------------------------------------------------------------------------
ProgressBar::ProgressBar(const std::string &title_in, const int cycle_count_in,
                         const ProgBarStyle style_in, std::ostream &output_in) :
    ProgressBar(title_in, std::string("["), std::string("]"), '#', ' ', cycle_count_in,
                style_in, output_in)
{}

//-------------------------------------------------------------------------------------------------
ProgressBar::ProgressBar(const int cycle_count_in, const ProgBarStyle style_in,
                         std::ostream &output_in) :
    ProgressBar(std::string(""), std::string("["), std::string("]"), '#', ' ', cycle_count_in,
                style_in, output_in)
{}

//-------------------------------------------------------------------------------------------------
const std::string& ProgressBar::getState() const {
  return bar_contents;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::displayTitle() {
  updateTerminalWidth();
  if (title_within_bar == false) {
    printf("%s\n", title.c_str());
  }
}
  
//-------------------------------------------------------------------------------------------------
void ProgressBar::reset() {
  progress = 0;
  update_called = false;
  last_percent = 0;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setCycleCount(const int cc_in) {
  if (cc_in <= 0) {
    rtErr("The number of iterations must be positive.", "ProgressBar", "setIterations");
  }
  cycle_count = cc_in;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setDoneChar(const char sym) {
  finished_mark = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setTodoChar(const char sym) {
  todo_mark = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setOpeningBracket(const std::string &sym) {
  open_bracket = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setClosingBracket(const std::string &sym) {
  close_bracket = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setStyle(const ProgBarStyle style_in) {
  style = style_in;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setTitle(const std::string &title_in) {
  title = title_in;
  title_width = title.size();
  updateTerminalWidth();
}
  
//-------------------------------------------------------------------------------------------------
void ProgressBar::setOutputStream(std::ostream &stream) {
  output = &stream;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setTerminalWidth(const int width) {
  if (width > 0) {
    terminal_width = width;
    allocateBarContents();
  } else {
    rtErr("Terminal width must be positive.", "ProgressBar", "setTerminalWidth");
  }
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::allocateBarContents() {

  // Resize the contents of bar_contents to fit terminal width
  bar_contents.clear();
  bar_contents.resize(terminal_width);

  // Initialize with the structure: [###     ] 50%
  size_t pos = 0;

  // Opening bracket
  for (char c : open_bracket) {
    bar_contents[pos++] = c;
  }

  // Title, if it is to be placed on the same line as the progress
  if (title_within_bar) {
    bar_contents[pos++] = ' ';
    for (char c : title) {
      bar_contents[pos++] = c;
    }
    bar_contents[pos++] = ' ';
    bar_contents[pos++] = '|';
    bar_contents[pos++] = ' ';
  }
  
  // Fill with initial todo marks (space by default)
  for (int i = 0; i < bar_width; ++i) {
    bar_contents[pos++] = todo_mark;
  }
  
  // Closing bracket
  for (char c : close_bracket) {
    bar_contents[pos++] = c;
  }

  // Add percentage placeholder with leading spaces.  Reserve 4 characters for the percentage
  // display.
  bar_contents[pos++] = ' ';
  for (int i = 0; i < 4; i++) {
    bar_contents[pos++] = ' ';
  }
}
  
//-------------------------------------------------------------------------------------------------
void ProgressBar::update() {
  if (cycle_count == 0) {
    rtErr("The number of cycles has not been set.", "ProgressBar", "update");
  }

  if (!update_called) {
    update_called = true;
  }

  // Calculate the percentage of progress
  int percent = (cycle_count > 1) ? progress * 100 / (cycle_count - 1) : 100;
  percent = std::min(percent, 100);
  
  // Only update if the percentage has changed
  if (percent != last_percent) {
    size_t num_done_chars = percent * bar_width / 100;
    size_t prev_done_chars = last_percent * bar_width / 100;

    // Update only the difference between the old and new percentage.
    // Update only if the bar is to be shown.
    size_t fill_start;
    if (title_within_bar) {
      fill_start = open_bracket.length() + 4 + title_width;
    }
    else {
      fill_start = open_bracket.length();
    }
    switch (style) {
    case ProgBarStyle::FULL:
      for (size_t i = prev_done_chars; i < num_done_chars; i++) {
        bar_contents[fill_start + i] = finished_mark;
      }
      for (size_t i = num_done_chars; i < bar_width; i++) {
        bar_contents[fill_start + i] = todo_mark;
      }
      break;
    case ProgBarStyle::PERCENT:
    case ProgBarStyle::NONE:
      break;
    }
    switch (style) {
    case ProgBarStyle::FULL:
    case ProgBarStyle::PERCENT:
      {            
        // Update the percentage display with right justification
        const std::string percent_string = std::to_string(percent) + "%";
        const size_t bar_total_len = bar_contents.size();
        size_t percent_pos = bar_total_len - percent_string.length();

        // First clear the percentage area
        for (size_t i = bar_total_len - 4; i < bar_total_len; i++) {
          bar_contents[i] = ' ';
        }

        // Then write the new percentage right-justified
        for (char c : percent_string) {
          bar_contents[percent_pos++] = c;
        }
        *output << "\r" << bar_contents << std::flush;
      }
      break;
    case ProgBarStyle::NONE:
      break;
    }
    last_percent = percent;
  }

  // Increment progress
  progress++;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::finalizeTerminalOutput() const {
  printf("\n");
}
  
//-------------------------------------------------------------------------------------------------
void ProgressBar::updateTerminalWidth() {
  if (isatty(fileno(stdout))) {
    struct winsize ws;
    if (ioctl(fileno(stdout), TIOCGWINSZ, &ws) == 0) {
      terminal_width = ws.ws_col;
    }
    else {
      terminal_width = default_output_file_width;
    }
  }
  else {
    terminal_width = default_output_file_width;
  }

  // Check the title width and update the bar width
  switch (style) {
  case ProgBarStyle::FULL:
    title_within_bar = (terminal_width - title_width >= 44);
    break;
  case ProgBarStyle::PERCENT:
    title_within_bar = (terminal_width - title_width >= 11);
    break;
  case ProgBarStyle::NONE:
    title_within_bar = true;
    break;
  }

  // Refill completion symbols based on last_percent
  if (title_within_bar) {
    bar_width = terminal_width - (title_width +
                                  open_bracket.length() + close_bracket.length() + 11);
  }
  else {
    bar_width = terminal_width - (open_bracket.length() + close_bracket.length() + 7);
  }
  
  // Re-allocate the bar_contents to match new terminal width
  allocateBarContents();

  int num_done_chars = last_percent * bar_width / 100;
  if (title_within_bar) {
    for (int i = 0; i < num_done_chars; i++) {
      bar_contents[open_bracket.length() + title_width + 4 + i] = finished_mark;
    }
  }
  else {
    for (int i = 0; i < num_done_chars; i++) {
      bar_contents[open_bracket.length() + i] = finished_mark;
    }
  }

  // Update the percentage text with right justification
  std::string percentString = std::to_string(last_percent) + "%";

  // First clear the percentage area
  for (size_t i = bar_contents.size() - 4; i < bar_contents.size(); ++i) {
    bar_contents[i] = ' ';
  }
  
  // Then write the new percentage right-justified
  size_t percentPos = bar_contents.size() - percentString.length();
  for (char c : percentString) {
    bar_contents[percentPos++] = c;
  }
}

} // namespace reporting
} // namespace stormm
