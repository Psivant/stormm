// -*-c++-*-
#ifndef STORMM_PROGRESSBAR_H
#define STORMM_PROGRESSBAR_H

#include <iostream>
#include <ostream>
#include <string>
#include <stdexcept>
#include "copyright.h"
#include "reporting_enumerators.h"

namespace stormm {
namespace reporting {

using display::ProgBarStyle;
  
class ProgressBar {
public:

  /// \brief The constructors for a ProgressBar may include the title and all visual elements.
  ///        All aspects of the object may be modified after construction.
  /// \{
  ProgressBar(const std::string &title_in, const std::string &open_bracket_in,
              const std::string &close_bracket_in, char finished_mark_in = '#',
              char todo_mark_in = ' ', int cycle_count_in = 0,
              ProgBarStyle style_in = ProgBarStyle::FULL, std::ostream &output_in = std::cout);

  ProgressBar(const std::string &title_in, int cycle_count_in = 0,
              ProgBarStyle style_in = ProgBarStyle::FULL, std::ostream &output_in = std::cout);

  ProgressBar(int cycle_count_in = 0, ProgBarStyle style_in = ProgBarStyle::FULL,
              std::ostream &output_in = std::cout);

  /// \brief The default copy and move constructors, as well as copy and move assignment operators,
  ///        are valid.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another ProgressBar object found on the right hand side of the assignment
  ///                  statement
  /// \{
  ProgressBar(const ProgressBar &original) = default;
  ProgressBar(ProgressBar &&original) = default;
  ProgressBar& operator=(const ProgressBar &original) = default;
  ProgressBar& operator=(ProgressBar &&original) = default;
  /// \}
  
  /// \brief Get the current state of the progress bar, for inspection.
  const std::string& getState() const;

  /// \brief Present an opportunity to display the progress bar's title, if it will not fit
  ///        within the bar itself (on the left hand side, next to the region that is filling
  ///        over time).  This will not print anything to the terminal unless the title is too
  ///        wide to display within the terminal next to the progress indicator.  However, calling
  ///        this routine will always update the object's sense of how wide the terminal is, and
  ///        any consequences thereof.
  void displayTitle();
  
  /// \brief Initializes a new Progress Bar from scratch, taking in user-specified information.
  ///        To be used when we want to do a hard-refresh on the existing ProgressBar object.
  ///
  /// \param  out     Determines which C standard output stream to use. Default is cout, can
  ///                 use cerr or C++ methods such as printf.
  void initialize(int n = 0, bool show_bar_in = true, std::ostream &out = std::cout);

  /// \brief Resets the current instance of the ProgressBar object (sets percentage to 0).
  ///        To be used before a new loop, with every other setting remaining the same.
  void reset();
  
  /// \brief Set a new number of iterations for an existing ProgressBar object.
  ///
  /// \param cc_in  The new number of iterations for the ProgressBar object
  void setCycleCount(int iter);

  /// \brief Set a new "done" char for a ProgressBar object.  This is the char that appears when
  ///        a percentage is done in the ProgressBar.
  ///
  /// \param sym  The symbol to use for completed work (default '#')
  void setDoneChar(char sym);

  /// \brief Set a new "todo" char for a ProgressBar object.  This is the char that populates the
  ///        remaining of the ProgressBar object.
  ///
  /// \param sym  The symbol to use for a work yet undone (default ' ')
  void setTodoChar(char sym);

  /// \brief Set a new opening bracket for a ProgressBar object.  This is the char during the
  ///        start of a ProgressBar object.
  ///
  /// \param sym  The string to use for an opening char of the bar (default is '[')
  void setOpeningBracket(const std::string &sym);
  
  /// \brief Set a new closing bracket for a ProgressBar object.  This is the char to use for a
  ///        closing char of the bar (default is ']')
  void setClosingBracket(const std::string &sym);

  /// \brief Change the display style of the ProgressBar.
  ///
  /// \param style  Indicate the detail in which progress is displayed
  void setStyle(ProgBarStyle style);

  /// \brief Set the title of a progress bar.  If there is room, the title will be displayed
  ///        ahead of the filling bar, separated by a '|' character.
  void setTitle(const std::string &title_in);
  
  /// \brief Function to set the output stream of the current ProgressBar object.  Default is
  ///        cout, can be changed to any std::ostream &object.
  ///
  /// \param stream   The standard output stream to use for rendering the progress bar.
  void setOutputStream(std::ostream &stream);

  /// \brief Allows the user to set a custom terminal width for the ProgressBar.
  ///        If this function is called, it will override the dynamically inferred width.
  ///
  /// \param width  The new width of the terminal in characters
  void setTerminalWidth(int width);

  /// \brief Function to allocate the bar_contents string for the ProgressBar object.
  ///        This is called when the bar is first initialized or when 
  ///        the terminal width is changed.
  void allocateBarContents();

  /// \brief Function to update the ProgressBar, incrementing the number of iterations by 1,
  ///        calculating the appropriate percentage, and rendering the ProgressBar in the
  ///        terminal.
  void update();

  /// \brief Finalize the last image of the progress bar and do not allow it to be overwritten in
  ///        the terminal by printing a carriage return.
  void finalizeTerminalOutput() const;
  
private:

  int progress;               ///< The amount of iterations that have been completed by 
                              ///<   the ProgressBar
  int cycle_count;            ///< The number of iterations the ProgressBar has to go 
                              ///<   through in total
  int last_percent;           ///< The percentage value that the ProgressBar calculated so far
  int terminal_width;         ///< The detected width of the terminal screen
  int title_width;            ///< The width of the title
  int bar_width;              ///< The width of the progress bar's fillable region
  ProgBarStyle style;         ///< Enumerated value determining the visual display.  Options
                              ///<   include to show a filling bar plus percentage, a
                              ///<   percentage only, or to silence the output altogether.
  bool update_called;         ///< Boolean flag to determine if the bar has been rendered 
                              ///<   at least once
  bool title_within_bar;      ///< Flag to indicate that the title fits within the terminal window
                              ///<   next to the displayed progress bar
  std::string title;          ///< Title of the progress bar.  This will be displayed on the
                              ///<   same line (part of the bar itself) if room is available,
                              ///<   or on a preceding line if there is non sufficient room
                              ///<   to print a meaningful progress bar beside the title.
  char finished_mark;         ///< Define the repeating character that fills the ProgressBar
  char todo_mark;             ///< Define the repeating character signifying unfinished work
  std::string open_bracket;   ///< String at the start of a ProgressBar
  std::string close_bracket;  ///< String at the end of a ProgressBar
  std::ostream* output;       ///< C std output stream throuth which to render the bar
  std::string bar_contents;   ///< The current state of the progress bar.  Updates that edit
                              ///<   the state of the progress bar will alter a minimal number
                              ///<   of characters in this string, and any display of the
                              ///<   progress bar will dump the entire contents of this string
                              ///<   to the chosen output mode.
  
  /// \brief  Function to dynamically infer terminal width to draw an appropriate ProgressBar.
  ///         This function is called by the constructor and the bar update function.  It does not
  ///         trigger based on user input.
  void updateTerminalWidth();
};

} // namespace reporting
} // namespace stormm

#endif // STORMM_PROGRESSBAR_H

