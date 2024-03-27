#include "copyright.h"
#include "FileManagement/file_util.h"
#include "Math/summation.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Parsing/textfile.h"
#include "error_format.h"
#include "report_table.h"
#include "summary_file.h"

namespace stormm {
namespace review {

using diskutil::openOutputFile;
using parse::addTailingWhiteSpace;
using parse::char4ToString;
using parse::minimalRealFormat;
using parse::NumberFormat;
using parse::realToString;
using parse::TextFile;
using parse::TextOrigin;
using stmath::sum;

//-------------------------------------------------------------------------------------------------
ReportTable::ReportTable(const std::vector<std::string> &data_in,
                         const std::vector<std::string> &column_headings_in,
                         const std::string &variable_name_in, const int format_width_in,
                         const std::vector<JustifyText> &alignments) :
    column_count{static_cast<int>(column_headings_in.size())},
    data_row_count{static_cast<int>(data_in.size() / column_headings_in.size())},
    header_row_count{1},
    format_width{format_width_in},
    data_kind{TableContentKind::STRING},
    rendered_data{data_in},
    column_headings{},
    column_widths{},
    variable_name{variable_name_in}
{
  // Check that the data size factorizes into the row and column counts
  if (static_cast<size_t>(column_count) * static_cast<size_t>(data_row_count) != data_in.size()) {
    rtErr("A data set of " + std::to_string(data_in.size()) + " points arranges into " +
          std::to_string(column_count) + " columns with a remainder of " +
          std::to_string(data_in.size() - (static_cast<size_t>(column_count) *
                                           static_cast<size_t>(data_row_count))) + ".",
          "ReportTable");
  }
  std::vector<size_t> data_widths(column_count);
  for (int i = 0; i < column_count; i++) {
    for (int j = 0; j < data_row_count; j++) {
      data_widths[i] = std::max(data_widths[i], rendered_data[(i * data_row_count) + j].size());
    }
  }
  
  // Format the column headings, taking cues from the data widths
  formatColumnHeadings(column_headings_in, data_widths);
    
  // Fill in the column widths based on the formatting determined for the headings
  findColumnWidths(data_widths, alignments);
}

//-------------------------------------------------------------------------------------------------
ReportTable::ReportTable(const std::vector<double> &data_in,
                         const std::vector<std::string> &column_headings_in,
                         const double precision, const std::string &variable_name_in,
                         const int format_width_in, const std::vector<JustifyText> &alignments) :
    ReportTable(dataAsStrings(data_in, column_headings_in.size(), precision), column_headings_in,
                variable_name_in, format_width_in, alignments)
{
  data_kind = TableContentKind::REAL;
}

//-------------------------------------------------------------------------------------------------
ReportTable::ReportTable(const std::vector<double> &data_in,
                         const std::vector<std::string> &column_headings_in,
                         const std::vector<double> &precision,
                         const std::string &variable_name_in, const int format_width_in,
                         const std::vector<JustifyText> &alignments) :
    ReportTable(dataAsStrings(data_in, precision), column_headings_in, variable_name_in,
                format_width_in, alignments)
{
  data_kind = TableContentKind::REAL;
}

//-------------------------------------------------------------------------------------------------
ReportTable::ReportTable(const std::vector<double> &data_in,
                         const std::vector<std::string> &column_headings_in,
                         const std::vector<int> &decimal_places,
                         const std::string &variable_name_in, const int format_width_in,
                         const std::vector<JustifyText> &alignments) :
    ReportTable(dataAsStrings(data_in, decimal_places), column_headings_in, variable_name_in,
                format_width_in, alignments)
{
  data_kind = TableContentKind::REAL;
}

//-------------------------------------------------------------------------------------------------
ReportTable::ReportTable(const std::vector<int> &data_in,
                         const std::vector<std::string> &column_headings_in,
                         const std::string &variable_name_in, const int format_width_in,
                         const std::vector<JustifyText> &alignments) :
    ReportTable(dataAsStrings(data_in, column_headings_in.size()), column_headings_in,
                variable_name_in, format_width_in, alignments)
{
  data_kind = TableContentKind::INTEGER;
}

//-------------------------------------------------------------------------------------------------
ReportTable::ReportTable(const std::vector<llint> &data_in,
                         const std::vector<std::string> &column_headings_in,
                         const std::string &variable_name_in, const int format_width_in,
                         const std::vector<JustifyText> &alignments) :
    ReportTable(dataAsStrings(data_in, column_headings_in.size()), column_headings_in,
                variable_name_in, format_width_in, alignments)
{
  data_kind = TableContentKind::INTEGER;
}

//-------------------------------------------------------------------------------------------------
ReportTable::ReportTable(const std::vector<char4> &data_in,
                         const std::vector<std::string> &column_headings_in,
                         const std::string &variable_name_in, const int format_width_in,
                         const std::vector<JustifyText> &alignments) :
    ReportTable(dataAsStrings(data_in, column_headings_in.size()), column_headings_in,
                variable_name_in, format_width_in, alignments)
{}

//-------------------------------------------------------------------------------------------------
int ReportTable::getColumnCount() const {
  return column_count;
}

//-------------------------------------------------------------------------------------------------
int ReportTable::getHeaderRowCount() const {
  return header_row_count;
}

//-------------------------------------------------------------------------------------------------
int ReportTable::getDataRowCount() const {
  return data_row_count;
}

//-------------------------------------------------------------------------------------------------
int ReportTable::getColumnWidth(const int column_index) const {
  if (column_index >= column_count) {
    rtErr("Column index " + std::to_string(column_index) + " is invalid for a table with " +
          std::to_string(column_count) + " columns.", "ReportTable", "getColumnWidth");
  }
  return column_widths[column_index];
}

//-------------------------------------------------------------------------------------------------
int ReportTable::getTableWidth() const {
  return format_width;
}

//-------------------------------------------------------------------------------------------------
TableContentKind ReportTable::getContentKind() const {
  return data_kind;
}

//-------------------------------------------------------------------------------------------------
const std::string& ReportTable::getValue(const size_t row_index, const size_t column_index) const {
  if (row_index >= data_row_count || column_index >= column_count) {
    rtErr("Table index (" + std::to_string(row_index) + ", " + std::to_string(column_index) +
          ") is invalid for a table with " + std::to_string(data_row_count) + " data rows and " +
          std::to_string(column_count) + " columns.", "ReportTable", "getValue");
  }
  return rendered_data[row_index * column_index];
}

//-------------------------------------------------------------------------------------------------
std::string ReportTable::printTable(const OutputSyntax style, const int width,
                                    const int data_row_start, const int data_row_end) const {

  // Determine the width with which to print the table.
  const int actual_width = (width < 0) ? format_width : width;
  
  // Determine the appropriate indentation for data rows based on the presence and length of a
  // protective marker.
  const char protect = commentSymbol(style);
  std::string comment_block(2, ' '), mpl_comment_block(4, ' ');
  const std::string space_block(2, ' '), carriage_return(1, '\n');
  const std::string mpl_bracket_block("  [ "), mpl_comma_block(", "), mpl_carriage_return(" ],\n");
  comment_block[0] = protect;
  mpl_comment_block[0] = protect;

  // Detect and validate the limits of data rows to print
  if (data_row_start < 0 || data_row_start >= data_row_count) {
    rtErr("Data row index " + std::to_string(data_row_start) + " is invalid for a table with " +
          std::to_string(data_row_count) + " rows.", "ReportTable", "printTable");
  }
  int actual_row_end;
  if (data_row_end <= 0) {
    actual_row_end = data_row_count;
  }
  else if (data_row_end >= data_row_count) {
    rtErr("Data row index " + std::to_string(data_row_end) + " is invalid for a table with " +
          std::to_string(data_row_count) + " rows.", "ReportTable", "printTable");
  }
  else {
    actual_row_end = data_row_end;
  }
  
  // Loop over all columns until the entire table has been printed.  Print in multiple stages if
  // not all lines can fit side-by-side in the output format.
  std::string result, hstack_line, hstack_declaration;
  const size_t est_line = sum<int>(column_widths) + (2 * column_widths.size()) + 1;
  size_t est_chars;
  switch (style) {
  case OutputSyntax::MATPLOTLIB:
    est_chars= static_cast<size_t>(header_row_count + 4 +
                                   (actual_row_end - data_row_start)) * (est_line + 5);
    break;
  case OutputSyntax::MATRIX_PKG:
  case OutputSyntax::STANDALONE:
    est_chars= static_cast<size_t>(header_row_count + 3 +
                                   (actual_row_end - data_row_start)) * est_line;
    break;
  }
  result.reserve(est_chars);
  int col_idx = 0;
  while (col_idx < column_count) {
    const int set_start = col_idx;
    int char_used = 0;
    bool first_in_set = true;
    int eol_chars;
    switch (style) {
    case OutputSyntax::MATPLOTLIB:
      eol_chars = 3;
      break;
    case OutputSyntax::MATRIX_PKG:
    case OutputSyntax::STANDALONE:
      eol_chars = 0;
      break;
    }
    while (col_idx < column_count &&
           (first_in_set || char_used + column_widths[col_idx] + 2 + eol_chars < actual_width)) {
      char_used += (first_in_set) ? column_widths[col_idx] + 4 : column_widths[col_idx] + 2;
      first_in_set = false;
      col_idx++;
    }
    const int set_end = col_idx;

    // Print the header elements.  The STANDALONE syntax places the variable name at the head of
    // the table, as it will not be displayed elsewhere.
    if (set_start == 0) {
      switch (style) {
      case OutputSyntax::MATPLOTLIB:
      case OutputSyntax::MATRIX_PKG:
        break;
      case OutputSyntax::STANDALONE:
        result += comment_block + variable_name + "\n";
        break;
      }
    }
    for (int i = 0; i < header_row_count; i++) {
      for (int j = set_start; j < set_end; j++) {
        switch (style) {
        case OutputSyntax::MATPLOTLIB:
          result += (j == set_start) ? mpl_comment_block : space_block;
          result += column_headings[(j * header_row_count) + i];
          break;
        case OutputSyntax::MATRIX_PKG:
        case OutputSyntax::STANDALONE:
          result += (j == set_start) ? comment_block : space_block;
          result += column_headings[(j * header_row_count) + i];
          break;
        }
      }
      result += "\n";
    }
    for (int i = set_start; i < set_end; i++) {
      switch (style) {
      case OutputSyntax::MATPLOTLIB:
        result += (i == set_start) ? mpl_comment_block : space_block;
        break;
      case OutputSyntax::MATRIX_PKG:
      case OutputSyntax::STANDALONE:
        result += (i == set_start) ? comment_block : space_block;
        break;
      }
      result += std::string(column_widths[i], '-');
    }
    result += "\n";

    // String-based table data is printed as part of the same comment block.  Numerical data
    // adjusts to the chosen output format.
    switch (data_kind) {
    case TableContentKind::INTEGER:
    case TableContentKind::REAL:
      
      // Print the variable name
      switch (style) {
      case OutputSyntax::MATPLOTLIB:
        if (set_start == 0 && set_end == column_count) {
          result += variable_name + " = np.array([\n";
        }
        else {
          result += variable_name + "_" + std::to_string(set_start + 1) + "_" +
                    std::to_string(set_end) + " = np.array([\n";
        }
        break;
      case OutputSyntax::MATRIX_PKG:
        if (set_start == 0 && set_end == column_count) {
          result += variable_name + " = [\n";
        }
        else {

          // Pre-allocate to optimize script performance in the report file
          if (set_start == 0) {
            result += variable_name + " = zeros(" +
                      std::to_string(actual_row_end - data_row_start) + ", " +
                      std::to_string(column_count) + ");\n";
          }
          result += variable_name + "(:, " + std::to_string(set_start + 1) + ":" +
                    std::to_string(set_end) + ") = [\n";
        }
        break;
      case OutputSyntax::STANDALONE:
        break;
      }

      // Print the data
      switch (style) {
      case OutputSyntax::MATPLOTLIB:
        for (int i = data_row_start; i < actual_row_end; i++) {
          for (int j = set_start; j < set_end; j++) {
            result += (j == set_start) ? mpl_bracket_block : mpl_comma_block;
            result += rendered_data[(j * data_row_count) + i];
          }
          result += mpl_carriage_return;
        }
        break;
      case OutputSyntax::MATRIX_PKG:
      case OutputSyntax::STANDALONE:
        for (int i = data_row_start; i < actual_row_end; i++) {
          for (int j = set_start; j < set_end; j++) {
            result += space_block;
            result += rendered_data[(j * data_row_count) + i];
          }
          result += carriage_return;
        }
        break;
      }

      // Print the closure
      switch (style) {
      case OutputSyntax::MATPLOTLIB:
        result += "]) " + comment_block + variable_name + "\n";
        if (set_start != 0 || set_end != column_count) {
          if (hstack_line.size() == 0) {
            hstack_declaration = variable_name + " = np.hstack((";
          }
          else {
            hstack_line += ", ";
          }
          hstack_line += variable_name + "_" + std::to_string(set_start + 1) + "_" +
                         std::to_string(set_end);
        }
        break;
      case OutputSyntax::MATRIX_PKG:
        result += "]; " + comment_block + variable_name + "\n";
        break;
      case OutputSyntax::STANDALONE:
        break;
      }
      if (set_end < column_count) {
        result += comment_block + carriage_return;
      }
      break;
    case TableContentKind::STRING:

      // Print the data behind comment characters
      for (int i = data_row_start; i < actual_row_end; i++) {
        for (int j = set_start; j < set_end; j++) {
          result += (j == set_start) ? comment_block : space_block;
          result += rendered_data[(j * data_row_count) + i];
        }
        result += carriage_return;
      }
      break;
    }
  }

  // Some format options require horizontal stacking of the arrays for different clusters of
  // matrix columns.
  switch (style) {
  case OutputSyntax::MATPLOTLIB:
    if (hstack_line.size() > 0) {
      hstack_line += " ))";
      result += indentText(hstack_line, hstack_declaration, 0, actual_width, false,
                           TextEnds::NEWLINE);
    }
    break;
  case OutputSyntax::MATRIX_PKG:
  case OutputSyntax::STANDALONE:
    break;
  }
  
  return result;
}

//-------------------------------------------------------------------------------------------------
void ReportTable::printTable(std::ofstream *foutp, const OutputSyntax style, const int width,
                             const int data_row_start, const int data_row_end) const {
  const std::string result = printTable(style, width, data_row_start, data_row_end);
  foutp->write(result.data(), result.size());
}

//-------------------------------------------------------------------------------------------------
void ReportTable::printTable(const std::string &file_name, const PrintSituation expectation,
                             const OutputSyntax style, const int width, const int data_row_start,
                             const int data_row_end) const {
  std::ofstream foutp = openOutputFile(file_name, expectation, "print a table to a file");
  printTable(&foutp, style, width, data_row_start, data_row_end);
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
void ReportTable::setVariableName(const std::string &variable_name_in) {
  variable_name = variable_name_in;
}

//-------------------------------------------------------------------------------------------------
void ReportTable::formatColumnHeadings(const std::vector<std::string> &column_headings_in,
                                       const std::vector<size_t> &data_widths) {
  std::vector<TextFile> heading_tf;
  heading_tf.reserve(column_count);
  int max_ch_lines = 0;
  for (int i = 0; i < column_count; i++) {
    const std::string ihead = indentText(column_headings_in[i], 0, data_widths[i]);
    heading_tf.emplace_back(ihead, TextOrigin::RAM);
    max_ch_lines = std::max(heading_tf.back().getLineCount(), max_ch_lines);
  }
  column_headings.resize(max_ch_lines * column_headings_in.size());
  for (int i = 0; i < column_count; i++) {
    const int ni_ln = heading_tf[i].getLineCount();
    for (int j = 0; j < ni_ln; j++) {
      column_headings[(i * max_ch_lines) + j] = heading_tf[i].getLineAsString(j);
    }
    for (int j = ni_ln; j < max_ch_lines; j++) {
      column_headings[(i * max_ch_lines) + j] = std::string("");
    }
  }
  header_row_count = max_ch_lines;
}

//-------------------------------------------------------------------------------------------------
void ReportTable::findColumnWidths(const std::vector<size_t> &data_widths,
                                   const std::vector<JustifyText> &alignments) {
  column_widths.resize(column_count);
  for (int i = 0; i < column_count; i++) {
    size_t max_header_width = 0;
    for (int j = 0; j < header_row_count; j++) {
      max_header_width = std::max(max_header_width,
                                  column_headings[(i * header_row_count) + j].size());
    }
    column_widths[i] = std::max(data_widths[i], max_header_width);
    for (int j = 0; j < header_row_count; j++) {
      addTailingWhiteSpace(&column_headings[(i * header_row_count) + j], column_widths[i]);
    }
    for (int j = 0; j < data_row_count; j++) {
      addTailingWhiteSpace(&rendered_data[(i * data_row_count) + j], column_widths[i]);
    }
    justifyStrings(&column_headings, i * header_row_count, (i + 1) * header_row_count,
                   JustifyText::CENTER, column_widths[i]);
    const JustifyText ialign = (alignments.size() > i) ? alignments[i] : JustifyText::RIGHT;
    justifyStrings(&rendered_data, i * data_row_count, (i + 1) * data_row_count, ialign,
                   column_widths[i]);
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> ReportTable::dataAsStrings(const std::vector<double> &data_in,
                                                    const int ncol, const double precision) {
  const size_t ndata = data_in.size();
  std::vector<std::string> result(ndata);
  for (size_t i = 0; i < ndata; i++) {
    result[i] = minimalRealFormat(data_in[i], precision);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> ReportTable::dataAsStrings(const std::vector<double> &data_in,
                                                    const std::vector<double> &precision) {
  const size_t ndata = data_in.size();
  std::vector<std::string> result(ndata);
  const int ncol = precision.size();
  for (size_t i = 0; i < ndata; i++) {
    const int col_idx = i % ncol;
    result[i] = minimalRealFormat(data_in[i], precision[col_idx]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> ReportTable::dataAsStrings(const std::vector<double> &data_in,
                                                    const std::vector<int> &decimal_places) {
  const size_t ndata = data_in.size();
  std::vector<std::string> result(ndata);
  const size_t ncol = decimal_places.size();
  const size_t nrow = ndata / ncol;
  for (size_t i = 0; i < ndata; i++) {
    const int col_idx = i / nrow;
    result[i] = realToString(data_in[i], decimal_places[col_idx] + 2, decimal_places[col_idx],
                             NumberFormat::STANDARD_REAL);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> ReportTable::dataAsStrings(const std::vector<int> &data_in,
                                                    const int ncol) {
  const size_t ndata = data_in.size();
  std::vector<std::string> result(ndata);
  for (size_t i = 0; i < ndata; i++) {
    result[i] = std::to_string(data_in[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> ReportTable::dataAsStrings(const std::vector<llint> &data_in,
                                                    const int ncol) {
  const size_t ndata = data_in.size();
  std::vector<std::string> result(ndata);
  for (size_t i = 0; i < ndata; i++) {
    result[i] = std::to_string(data_in[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> ReportTable::dataAsStrings(const std::vector<char4> &data_in,
                                                    const int ncol) {
  const size_t ndata = data_in.size();
  std::vector<std::string> result(ndata);
  for (size_t i = 0; i < ndata; i++) {
    result[i] = char4ToString(data_in[i]);
  }

  return result;
}

} // namespace review
} // namespace stormm
