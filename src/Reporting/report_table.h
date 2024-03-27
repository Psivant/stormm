// -*-c++-*-
#ifndef STORMM_REPORT_TABLE_H
#define STORMM_REPORT_TABLE_H

#include <string>
#include <vector>
#include "copyright.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "FileManagement/file_enumerators.h"
#include "Parsing/parsing_enumerators.h"
#include "reporting_enumerators.h"
#include "summary_file.h"

namespace stormm {
namespace review {

using diskutil::PrintSituation;
using parse::JustifyText;
  
/// \brief Tables in output file sections work differently than tables dumped to the terminal
///        window.  These tables will take in an array of data with some dimensions and then find
///        the best way to organize it, storing all columns and rows as a series of strings.
class ReportTable {
public:

  /// \brief The constructor takes multiple input formats but converts all data to strings for
  ///        internal storage, negating the need for templating.
  /// \{
  ReportTable(const std::vector<std::string> &data_in,
              const std::vector<std::string> &column_headings,
              const std::string &variable_name_in = std::string(""),
              int format_width_in = default_output_file_width,
              const std::vector<JustifyText> &alignments = {});

  ReportTable(const std::vector<double> &data_in, const std::vector<std::string> &column_headings,
              double precision = 1.0e-6, const std::string &variable_name_in = std::string(""),
              int format_width_in = default_output_file_width,
              const std::vector<JustifyText> &alignments = {});

  ReportTable(const std::vector<double> &data_in, const std::vector<std::string> &column_headings,
              const std::vector<double> &precision,
              const std::string &variable_name_in = std::string(""),
              int format_width_in = default_output_file_width,
              const std::vector<JustifyText> &alignments = {});

  ReportTable(const std::vector<double> &data_in, const std::vector<std::string> &column_headings,
              const std::vector<int> &decimal_places,
              const std::string &variable_name_in = std::string(""),
              int format_width_in = default_output_file_width,
              const std::vector<JustifyText> &alignments = {});

  ReportTable(const std::vector<int> &data_in, const std::vector<std::string> &column_headings,
              const std::string &variable_name_in = std::string(""),
              int format_width_in = default_output_file_width,
              const std::vector<JustifyText> &alignments = {});

  ReportTable(const std::vector<llint> &data_in, const std::vector<std::string> &column_headings,
              const std::string &variable_name_in = std::string(""),
              int format_width_in = default_output_file_width,
              const std::vector<JustifyText> &alignments = {});
  
  ReportTable(const std::vector<char4> &data_in, const std::vector<std::string> &column_headings,
              const std::string &variable_name_in = std::string(""),
              int format_width_in = default_output_file_width,
              const std::vector<JustifyText> &alignments = {});
  /// \}

  /// \brief With no const members or pointers to repair and all Standard Template Library
  ///        components, the ReportTable can be copied or moved using default constructors and
  ///        assignement operators.
  /// \{
  ReportTable(const ReportTable &original) = default;
  ReportTable(ReportTable &&original) = default;
  ReportTable& operator=(const ReportTable &original) = default;
  ReportTable& operator=(ReportTable &&original) = default;
  /// \}

  /// \brief Get the number of columns.
  int getColumnCount() const;

  /// \brief Get the number of header rows.
  int getHeaderRowCount() const;

  /// \brief Get the number of data rows.
  int getDataRowCount() const;

  /// \brief Get the width of a particular column.
  ///
  /// \param column_index  Index of the column of interest, starting from zero
  int getColumnWidth(int column_index) const;

  /// \brief Get the internal format width setting.
  int getTableWidth() const;

  /// \brief Get the type of content in the table.  All content is stored as strings, but depending
  ///        on whether it started as real numbers, integers, or string data there may be
  ///        significance.
  TableContentKind getContentKind() const;
  
  /// \brief Get the rendered data for a particular column and row.
  ///
  /// \param row_index     Index of the row of interest, starting from zero
  /// \param column_index  Index of the column of interest, starting from zero
  const std::string& getValue(size_t row_index, size_t column_index) const;

  /// \brief Print the table based on some external directives.
  ///
  /// Overloaded:
  ///   - Print to a formatted string containing appropriate carriage returns
  ///   - Print to a named file based on an expectation about the state in which that file will
  ///     be found
  ///   - Print to an open file using a file pointer
  ///
  /// \param variable_name   Name of the variable to assign at the beginning of table data rows.
  ///                        The entire contents of the table will fall under this variable name.
  /// \param style           One of several options for formatting the non-protected portion of
  ///                        the table to be amenable to certain plotting programs
  /// \param width           Width of the output file (ASCII characters per line) to be respected
  ///                        when printing the table
  /// \param data_row_start  The first data row to print
  /// \param data_row_end    The upper limit of data rows to print (this row index is the first
  ///                        that will not be printed).  The default of -1 implies printing all
  ///                        rows.
  /// \param foutp           File stream pointer for the output
  /// \param file_name       Name of the file to open and write the table into
  /// \param expectation     The state in which the output file is to be found to permit writing
  /// \{
  std::string printTable(OutputSyntax style, int width = -1, int data_row_start = 0,
                         int data_row_end = -1) const;

  void printTable(std::ofstream *foutp, OutputSyntax style, int width = -1, int data_row_start = 0,
                  int data_row_end = -1) const;

  void printTable(const std::string &file_name, PrintSituation expectation,
                  OutputSyntax style, int width = -1, int data_row_start = 0,
                  int data_row_end = -1) const;
  /// \}

  /// \brief Set the variable name associated with this table
  ///
  /// \param variable_name_in  The variable name to assign
  void setVariableName(const std::string &variable_name_in);
  
private:
  int column_count;                          ///< The number of columns defined for the table
  int data_row_count;                        ///< The number of data rows defined for the table
  int header_row_count;                      ///< Number of rows defined for the table column
                                             ///<   headings.  Multiple header rows may be alotted
                                             ///<   if the headers contain multiple words and the
                                             ///<   format can benefit from compacting.
  int format_width;                          ///< The width with which the table is to be printed.
                                             ///<   A new, positive parameter fed to one of the
                                             ///<   printer functions will override this setting.
  TableContentKind data_kind;                ///< The type of data reported in the table data
  std::vector<std::string> rendered_data;    ///< Tabulated data, converted to strings with the
                                             ///<   appropriate number of significant figures
  std::vector<std::string> column_headings;  ///< Column headings, presenting verbatim the
                                             ///<   information taken in by the constructor but
                                             ///<   allowing for multiple rows if the headings
                                             ///<   contain multiple words and the data format
                                             ///<   could be tightened by squeezing the header.
  std::vector<int> column_widths;            ///< Widths of each column needed to express the data
                                             ///<   (all numerical data will be right-justified,
                                             ///<   while column headings will be centered)
  std::string variable_name;                 ///< Name of the variable associated with this table.
                                             ///<   Not relevant to tables printed in "STANDALONE"
                                             ///<   style.  Matrix packages and plotting programs
                                             ///<   will read the table as a matrix associated
                                             ///<   with this variable name.
  
  /// \brief Format the column headings, taking cues from the data widths and overall table width.
  ///        This function will translate the input column headings into the internally stored
  ///        column_headings array, which will be header_row_count * column_count in length.
  ///
  /// \param column_headings_in  The input column headings, each a single string
  /// \param data_widths         Widths of each data column (all members of each data column will
  ///                            be strings of the same length after the initial rendering process)
  void formatColumnHeadings(const std::vector<std::string> &column_headings_in,
                            const std::vector<size_t> &data_widths);

  /// \brief Determine the appropriate column widths for the table, based on the longest element in
  ///        each column's data strings and the longest word in each column's header strings.
  ///
  /// \param data_widths  Widths of each data column (all members of each data column will be
  ///                     strings of the same length after the initial rendering process)
  /// \param alignments   The manner in which to align each column
  void findColumnWidths(const std::vector<size_t> &data_widths,
                        const std::vector<JustifyText> &alignments = {});
  
  /// \brief Format the data entries.  Various overloads cover the types of data that might be
  ///        encountered, many with different treatments, such that templating is neither necessary
  ///        nor effective.
  ///
  /// Overloaded:
  ///   - Take real-valued input
  ///   - Take signed integer input
  ///   - Take char4 tuple input
  ///
  /// \param data_in         The input data, in various formats
  /// \param ncol_in         The number of columns into which the data is divided
  /// \param precision       Relative precision used to represent real-valued data (if this is
  ///                        given as a vector, its length will be assumed to be the number of
  ///                        columns)
  /// \param decimal_places  Number of digits after the decimal used to represent real-valued data
  /// \{
  std::vector<std::string> dataAsStrings(const std::vector<double> &data_in, int ncol,
                                         double precision);

  std::vector<std::string> dataAsStrings(const std::vector<double> &data_in,
                                         const std::vector<double> &precision);

  std::vector<std::string> dataAsStrings(const std::vector<double> &data_in,
                                         const std::vector<int> &decimal_places);

  std::vector<std::string> dataAsStrings(const std::vector<int> &data_in, int ncol);
  
  std::vector<std::string> dataAsStrings(const std::vector<llint> &data_in, int ncol);

  std::vector<std::string> dataAsStrings(const std::vector<char4> &data_in, int ncol);
  /// \}
};

} // namespace review
} // namespace stormm
  
#endif
