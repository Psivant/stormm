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
ProgressBar::ProgressBar(int n, bool showbar, std::ostream& out) :
    progress(0), nCycles(n), lastPercent(0), doShowBar(showbar), 
    updateIsCalled(false), doneChar("#"), todoChar(" "), 
    openingBracketChar("["), closingBracketChar("]"), output(&out)
{
  updateTerminalWidth();
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::initialize(int n, bool showbar, std::ostream& out) {
  progress = 0;
  nCycles = n;
  lastPercent = 0;
  doShowBar = showbar;
  updateIsCalled = false;
  doneChar = "#";
  todoChar = " ";
  openingBracketChar = "[";
  closingBracketChar = "]";
  output = &out;
  updateTerminalWidth();
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::reset() {
  progress = 0;
  updateIsCalled = false;
  lastPercent = 0;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setIterations(int iter) {
  if (iter <= 0) {
    rtErr("The number of iterations must be positive.", "ProgressBar", "setIterations");
  }
  nCycles = iter;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setDoneChar(const std::string& sym) {
  doneChar = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setTodoChar(const std::string& sym) {
  todoChar = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setOpeningBracketChar(const std::string& sym) {
  openingBracketChar = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setClosingBracketChar(const std::string& sym) {
  closingBracketChar = sym;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::showBar(bool flag) {
  doShowBar = flag;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::setOutputStream(std::ostream& stream) {
  output = &stream;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::update() {
  if (nCycles == 0) {
    rtErr("The number of cycles has not been set.", "ProgressBar", "update");
  }
  if (updateIsCalled == false) {
    updateIsCalled = true;
  }
  
  // Subtract characters for brackets and the percentage filled.
  const int percent = progress * 100 / (nCycles - 1);
  int barWidth = termWidth - (openingBracketChar.length() + closingBracketChar.length() + 7);

  // Ensure barWidth is at least 10 to accommodate the minimum length of the bar.
  barWidth = std::max(barWidth, 10);

  // Redraw the entire progress bar for each update.
  if (doShowBar) {

    // Move to the beginning of the line.
    *output << "\r" << openingBracketChar;

    // Draw the progress part of the bar.
    int numDoneChars = percent * barWidth / 100;
    for (int i = 0; i < numDoneChars; ++i) {
      *output << doneChar;
    }
    
    // Draw the remaining part of the bar.
    for (int i = numDoneChars; i < barWidth; ++i) {
      *output << todoChar;
    }
    
    // Close the bar and display the percentage.
    *output << closingBracketChar << ' ' << percent << '%';
  }
  else {

    // Just print the percentage if the bar is not being shown.
    if (percent < 10) {
      *output << "\r  " << percent << '%';
    }
    else if (percent < 100) {
      *output << "\r " << percent << '%';
    }
    else {
      *output << "\r" << percent << '%';
    }
  }

  // Flush the output to update the display and increment the internal counters.
  *output << std::flush;
  lastPercent = percent;
  ++progress;
}

//-------------------------------------------------------------------------------------------------
void ProgressBar::updateTerminalWidth() {
  if (isatty(fileno(stdout))) {
    struct winsize ws;
    if (ioctl(fileno(stdout), TIOCGWINSZ, &ws) == 0) {
      termWidth = ws.ws_col;
    }
    else {
      termWidth = default_output_file_width;
    }
  }
  else {
    termWidth = default_output_file_width;
  }
}

} // namespace reporting
} // namespace stormm
