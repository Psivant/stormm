// -*-c++-*-
#ifndef STORMM_PROGRESSBAR_H
#define STORMM_PROGRESSBAR_H

//-------------------------------------------------------------------------------------------------
// Copyright (c) 2017 Luigi Pertoldi
//
// Modified by Psivant Therapeutics, 2024
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software
// and associated documentation files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or
// substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
// BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
// DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
// Progress bar implementation, originally created by Luigi Pertoldi
// Modified by Kush Srivastava, Psivant Therapeutics to suit the STORMM framework
// 
// Original: https://github.com/gipert/progressbar/tree/master
//-------------------------------------------------------------------------------------------------

#include <iostream>
#include <ostream>
#include <string>
#include <stdexcept>

namespace stormm {
namespace reporting {

class ProgressBar {
public:
  /// \brief The constructor for a ProgressBar.
  ///
  /// \param n        Number of iterations the ProgressBar is divided into.
  /// \param showbar  Indicate whether to show a visual progress bar or just a percentage
  /// \param out      Indicate the C standard output stream to use.  The default is cout, but the
  ///                 bar can use cerr or C++ methods such as printf.
  ProgressBar(int n = 0, bool showbar = true, std::ostream &out = std::cout);
  
  /// \brief  Initializes a new Progress Bar from scratch, taking in user-specified information.
  ///         To be used when we want to do a hard-refresh on the existing ProgressBar object.
  ///
  /// \param  n       Number of iterations the ProgressBar is divided into.
  /// \param  showbar Determines if we want to show a visual progress bar or just a percentage
  /// \param  out     Determines which C standard output stream to use. Default is cout, can
  ///                 use cerr or C++ methods such as printf.
  void initialize(int n = 0, bool showbar = true, std::ostream &out = std::cout);

  /// \brief  Resets the current instance of the ProgressBar object (sets percentage to 0).
  ///         To be used before a new loop, with every other setting remaining the same.
  void reset();
  
  /// \brief  Set a new number of iterations for an existing ProgressBar object.
  ///
  /// \param iter The new number of iterations for the ProgressBar object
  void setIterations(int iter);

  /// \brief  Set a new "done" char for a ProgressBar object.  This is the char that appears when
  ///         a percentage is done in the ProgressBar.
  ///
  /// \param sym  The string to use for a Done char (default is '#')
  void setDoneChar(const std::string &sym);

  /// \brief  Set a new "todo" char for a ProgressBar object.  This is the char that populates the
  ///         remaining of the ProgressBar object.
  ///
  /// \param sym  The string to use for a Todo char (default is ' ')
  void setTodoChar(const std::string &sym);

  /// \brief  Set a new opening bracket for a ProgressBar object.  This is the char during the
  ///         start of a ProgressBar object.
  ///
  /// \param sym  The string to use for an opening char of the bar (default is '[')
  void setOpeningBracketChar(const std::string &sym);
  
  /// \brief  Set a new closing bracket for a ProgressBar object.  This is the char to use for a
  ///         closing char of the bar (default is ']')
  void setClosingBracketChar(const std::string &sym);

  /// \brief  Toggles the visibility of the ProgressBar. Default is true.
  ///
  /// \param flag  Boolean value to show/hide the ProgressBar. If hidden, the percentage value is
  ///              still shown.
  void showBar(bool flag = true);
  
  /// \brief  Function to set the output stream of the current ProgressBar object.  Default is
  ///         cout, can be changed to any std::ostream &object.
  ///
  /// \param stream   The standard output stream to use for rendering the progress bar.
  void setOutputStream(std::ostream &stream);
  
  /// \brief  Function to update the ProgressBar, incrementing the number of iterations by 1,
  ///         calculating the appropriate percentage, and rendering the ProgressBar in the
  ///         terminal.
  void update();

private:

  int progress;                   ///< The amount of iterations that have been completed by 
                                  ///  the ProgressBar
  int nCycles;                    ///< The number of iterations the ProgressBar has to go 
                                  ///  through in total
  int lastPercent;                ///< The percentage value that the ProgressBar calculated so far
  bool doShowBar;                 ///< Boolean value determining if we want to show the 
                                  ///  visual ProgressBar
  bool updateIsCalled;            ///< Boolean flag to determine if the bar has been rendered 
                                  ///  at least once
  std::string doneChar;           ///< String that appears when a percentage is done
  std::string todoChar;           ///< String that appears for the rest of the ProgressBar
  std::string openingBracketChar; ///< String at the start of a ProgressBar
  std::string closingBracketChar; ///< String at the end of a ProgressBar

  std::ostream* output;           ///< C std output stream to render the bar through.

  int termWidth;                  ///< Stores the width of the terminal screen
  
  /// \brief  Function to dynamically infer terminal width to draw an appropriate ProgressBar.
  ///         This function is called implicitly and does not require user input.
  void updateTerminalWidth();
};

} // namespace reporting
} // namespace stormm

#endif // STORMM_PROGRESSBAR_H
