// -*-c++-*-
#ifndef STORMM_SUMMARY_FILE_H
#define STORMM_SUMMARY_FILE_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "copyright.h"
#include "FileManagement/file_enumerators.h"
#include "Parsing/textfile.h"
#include "reporting_enumerators.h"

namespace stormm {
namespace review {

using diskutil::PrintSituation;
  
/// \brief Default width for STORMM output files intended for human parsing.
constexpr int default_output_file_width = 80;
  
/// \brief Find the ideal format width for an output stream.  This decides based on whether the
///        output is going to the terminal (std::cout) or a file.  In the former case, STORMM
///        will try to print to the terminal's full width.  In the latter case, STORMM will use a
///        preset width given by the default_output_file_width constant.
int findFormatWidth(const std::ostream *foutp);
  
/// \brief Print the "splash screen" to stdout or to a file.  The splash screen provides the
///        authorship, origin of the software, and customary disclaimer.
///
/// Overloaded:
///   - Produce a string containing the splash screen contents, formatted for the width indicated
///   - Print a block of text to a file indicated by a file pointer
///
/// \param width  The width of screen or file to contain the information and disclaimer
/// \param foutp  Writeable file stream to accept the output
/// \{
std::string stormmSplash(int width);
void stormmSplash(std::ofstream *foutp, int width = default_output_file_width);
void stormmSplash(std::ostream *foutp = &std::cout, int width = default_output_file_width);
/// \}

/// \brief Print the "watermark" signifying the completion of writing of a STORMM output file.
///
/// Overloaded:
///   - Produce a string containing the splash screen contents, formatted for the width indicated
///   - Print a block of text to a file indicated by a file pointer
///
/// \param width  The The width of screen or file to contain the information and disclaimer
/// \param foutp  Writeable file stream to accept the output
/// \{
std::string stormmWatermark(int width);
void stormmWatermark(std::ofstream *foutp, int width = default_output_file_width);
void stormmWatermark(std::ostream *foutp = &std::cout, int width = default_output_file_width);
/// \}

/// \brief Summarize the work done in a particular STORMM run.
///
/// \param description  Content of the header to write, a description of the program or activity
/// \param foutp        Writeable file stream for the output file
/// \param file_width   Width of the file to write (default -1 will choose the default file width)
void summaryHeader(const std::string &description, std::ofstream *foutp, int file_width = -1);

/// \brief Obtain the proper comment symbol for a given output syntax.  The comment symbol will be
///        used to shield commentary from interpretation, given the assumption that a particular
///        program or sort of program will be used to run the output file as a script and plot the
///        data.
///
/// Overloaded:
///   - Produce comment symbols for report file syntax
///   - Produce comment symbols for grid files
///
/// \param syntax  Indicate the type of program which will interpret the output, or the file type
/// \{
char commentSymbol(OutputSyntax format);
char commentSymbol(GridFileSyntax format);
/// \}
  
/// \brief Indent text, possibly behind a special marker symbol.  The protectText() function below
///        calls this with zero indentation to produce text protect behind its chosen marker
///        symbol.
///
/// Overloaded:
///   - Indent with a special marker symbol included
///   - Indent with no marker symbol
///
/// \param text               The text to indent and format
/// \param marker             A marker to include ahead of any text.  If the marker has nonzero
///                           length, an additional space will be placed between the marker and any
///                           text, and the marker will extend the overall indentation by its own
///                           length.
/// \param indent             The number of white space characters to print prior to the marker
///                           on the first line and on all subsequent lines.  
/// \param width              The total width of formatted text, including the indentation
/// \param mark_on_all_lines  Indicate whether the marker should be included as part of the
///                           indentation on all lines.  If set to false, indentation on all
///                           subsequent lines will be increased by the length of the marker plus
///                           an additional space if the marker has nonzero length.
/// \{
std::string indentText(const std::string &text, const std::string &marker, int indent, int width,
                       bool mark_on_all_lines, TextEnds block_ending = TextEnds::AS_IS);

std::string indentText(const std::string &text, int indent, int width,
                       TextEnds block_ending = TextEnds::AS_IS);
/// \}

/// \brief Prepare a text string for output in a fixed format width with a protective character at
///        the start of each line.
///
/// \param text     The text to print
/// \param protect  Comment character used to protect the text
/// \param width    Total width of the file to produce (-1 or any other value <= 0 will choose the
///                 default file width)
std::string protectText(const std::string &text, char protect, int width = -1,
                        TextEnds block_ending = TextEnds::NEWLINE);

/// \brief Print a block of text, within a format limit, protected by a comment symbol.
///
/// Overloaded:
///   - Print a block of text to a file indicated by a file pointer
///   - Print a block of text to a named file, opened according to the provided protocol
///
/// \param text         The text to print
/// \param protect      Comment character used to protect the text
/// \param foutp        The file pointer to direct text into
/// \param ouptut_file  The name of the file to open
/// \param file_width   Total width of the file to produce (-1 or any other value <= 0 will choose
///                     the default file width)
/// \param expectation  State in which the file is to be written (default is to append)
/// \{
void printProtectedText(const std::string &text, char protect, std::ofstream *foutp,
                        int file_width = -1, TextEnds block_ending = TextEnds::NEWLINE);
  
void printProtectedText(const std::string &text, char protect, const std::string &output_file,
                        int file_width = -1, PrintSituation expectation = PrintSituation::APPEND,
                        TextEnds block_ending = TextEnds::NEWLINE);
/// \}

} // namespace review
} // namespace stormm


#endif
