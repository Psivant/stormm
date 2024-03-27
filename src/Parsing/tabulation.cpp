#include <algorithm>
#include <cmath>
#include "copyright.h"
#include "Math/summation.h"
#include "Reporting/error_format.h"
#include "parse.h"
#include "tabulation.h"

namespace stormm {
namespace parse {

using stmath::sum;

//-------------------------------------------------------------------------------------------------
void tableHorizontalRule(const std::vector<int> &column_widths, const BorderFormat borders,
                         const int indentation, const int column_spacing,
                         const bool have_categories) {
  const int n_cols = column_widths.size();
  int table_width = sum<int>(column_widths) + ((n_cols - 2) * column_spacing);
  switch (borders) {
  case BorderFormat::NONE:
    table_width += column_spacing;
    break;
  case BorderFormat::LIGHT:
    table_width += 3;
    break;
  case BorderFormat::FULL:
    table_width += 7;
    break;
  }
  std::string ruler(table_width, '-');
  switch (borders) {
  case BorderFormat::NONE:
  case BorderFormat::LIGHT:
    {
      int running_sum = 0;
      if (have_categories) {
	if (borders == BorderFormat::NONE) {
          for (int i = 0; i < column_spacing; i++) {
	    ruler[column_widths[0] + i] = ' ';
	  }
	  running_sum += column_widths[0] + column_spacing;
	}
	else {
          ruler[column_widths[0] + 1] = '+';
          ruler[column_widths[0] + 2] = ' ';
          running_sum += column_widths[0] + 3;
        }
      }
      for (int i = have_categories; i < n_cols - 1; i++) {
	running_sum += column_widths[i];
	for (int j = 0; j < column_spacing; j++) {
	  ruler[running_sum + j] = ' ';
	}
	running_sum += column_spacing;
      }
    }
    break;
  case BorderFormat::FULL:
    ruler[0] = '+';
    ruler[table_width - 1] = '+';
    ruler[column_widths[0] + 3] = '+';
    break;
  }
  for (int i = 0; i < indentation; i++) {
    printf(" ");
  }
  printf("%s\n", ruler.c_str());
}

//-------------------------------------------------------------------------------------------------
void tableHeader(const std::vector<std::vector<std::string>> &heading_words,
                 const std::vector<int> &column_widths, const BorderFormat borders,
                 const int indentation, const int column_spacing, const bool have_categories) {
  const int n_cols = static_cast<int>(column_widths.size());
  std::vector<int> words_used(n_cols, 0);
  bool print_on = true;
  while (print_on) {

    // Print as many words in each column as will fit
    for (int i = 0; i < indentation; i++) {
      printf(" ");
    }
    if (borders == BorderFormat::FULL) {
      printf("| ");
    }
    for (int i = 0; i < n_cols; i++) {
      std::string col_string = "";
      int chars_used = 0;
      int wu_i = words_used[i];
      while (wu_i < static_cast<int>(heading_words[i].size()) &&
             chars_used + (chars_used > 0) + static_cast<int>(heading_words[i][wu_i].size()) <=
             column_widths[i]) {
	if (chars_used > 0) {
          chars_used += 1 + heading_words[i][wu_i].size();
          col_string += " " + heading_words[i][wu_i];
	}
	else {
          chars_used += heading_words[i][wu_i].size();
          col_string += heading_words[i][wu_i];
	}
        wu_i++;
      }
      words_used[i] = wu_i;
      const int ccount = (column_widths[i] - col_string.size()) / 2;
      for (int j = 0; j < ccount; j++) {
	printf(" ");
      }
      printf("%s", col_string.c_str());
      for (int j = ccount + col_string.size(); j < column_widths[i]; j++) {
        printf(" ");
      }
      if (i < have_categories) {
	switch (borders) {
	case BorderFormat::NONE:
	  for (int j = 0; j < column_spacing; j++) {
	    printf(" ");
	  }
	  break;
	case BorderFormat::LIGHT:
	case BorderFormat::FULL:
          printf(" | ");
          break;
	}
      }
      else if (i < n_cols - 1) {
        for (int j = 0; j < column_spacing; j++) {
          printf(" ");
	}
      }
    }
    switch (borders) {
    case BorderFormat::FULL:
      printf(" |\n");
      break;
    case BorderFormat::LIGHT:
    case BorderFormat::NONE:
      printf("\n");
      break;
    }
    
    // Check whether all columns' complete header text been printed
    print_on = false;
    for (int i = 0; i < n_cols; i++) {
      print_on = (words_used[i] < static_cast<int>(heading_words[i].size()) || print_on);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void printTable(const std::vector<std::string> &headings,
                const std::vector<std::string> &categories,
                const std::vector<std::vector<PolyNumeric>> &column_data,
                const std::vector<NumberFormat> &column_formats, const std::string &column_sort,
                const BorderFormat borders, const int indentation, const int column_spacing,
                const int max_decimal_places_in, const double precision,
                const SortDirection sorting_direction) {

  // Checks on input
  const bool have_categories = (categories.size() > 0);
  const int n_rows = (have_categories) ? categories.size() : column_data[0].size();
  const int n_cols = column_data.size() + have_categories;
  if (headings.size() != n_cols) {
    std::string category_names_active = (have_categories) ?
                                        "  A column of category names is also present." : "";
    rtErr("Tabulation cannot proceed with " + std::to_string(headings.size()) + " headings and " +
	  std::to_string(column_data.size()) + " data fields." + category_names_active,
          "printTable");
  }
  if (have_categories) {
    if (static_cast<int>(column_formats.size()) != n_cols - 1) {
      rtErr("Tabulation requires " + std::to_string(n_cols - 1) + " formats, but only " +
	    std::to_string(column_formats.size()) + " were supplied.", "printTable");
    }
  }
  for (int i = 0; i < n_cols - have_categories; i++) {
    if (static_cast<int>(column_data[i].size()) != n_rows) {
      rtErr("Column " + std::to_string(i + 1) + " contains data for " +
            std::to_string(column_data[i].size()) + "rows, but the table expects " +
            std::to_string(n_rows), "printTable");
    }
  }

  // Compute the numerical format width necessary for each column
  const int max_decimal_places = std::min(static_cast<int>(ceil(log10(precision))),
                                          max_decimal_places_in);
  std::vector<int> column_widths(n_cols, 0);
  std::vector<int> column_decimals(n_cols, 0);
  if (have_categories) {
    for (int i = 0; i < n_rows; i++) {
      int k = 0;
      for (int j = 0; j < categories[i].size(); j++) {
	if (categories[i][j] != '\n') {
	  k++;
	}
	else {
	  column_widths[0] = std::max(k, column_widths[0]);
	  k = 0;
	}
      }
      if (k > 0) {
        column_widths[0] = std::max(k, column_widths[0]);
      }
    }
  }
  for (int i = have_categories; i < n_cols; i++) {
    const int data_idx = i - have_categories;
    int format_spacing;
    switch (column_formats[data_idx]) {
    case NumberFormat::SCIENTIFIC:
      format_spacing = 6 + max_decimal_places;
      column_decimals[data_idx] = max_decimal_places;
      break;
    case NumberFormat::STANDARD_REAL:
      format_spacing = 1;
      break;
    case NumberFormat::INTEGER:
    case NumberFormat::LONG_LONG_INTEGER:
    case NumberFormat::UNSIGNED_INTEGER:
    case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
    case NumberFormat::CHAR4:
      format_spacing = 0;
      break;
    }
    for (int j = 0; j < n_rows; j++) {
      double ij_number;
      switch(column_formats[data_idx]) {
      case NumberFormat::STANDARD_REAL:
      case NumberFormat::SCIENTIFIC:
        ij_number = column_data[data_idx][j].d;
        break;
      case NumberFormat::INTEGER:
        ij_number = column_data[data_idx][j].i;
        break;
      case NumberFormat::LONG_LONG_INTEGER:
        ij_number = column_data[data_idx][j].lli;
        break;
      case NumberFormat::UNSIGNED_INTEGER:
        ij_number = column_data[data_idx][j].ui;
        break;
      case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
        ij_number = column_data[data_idx][j].ulli;
        break;
      case NumberFormat::CHAR4:
	break;
      }
      int ij_width = format_spacing;

      // A minus sign will add to the width for any format
      ij_width += (ij_number < 0.0);
      ij_number = fabs(ij_number);

      // Extend the format width as necessary to accommodate digits after the decimal
      switch(column_formats[data_idx]) {
      case NumberFormat::STANDARD_REAL:
	{
          const double ij_part = ij_number / precision;          
	  int ij_dec;
          if (ij_number >= 1.0001) {
            ij_width += static_cast<int>(ceil(log10(ij_number)));
          }
          else {
            ij_width += 1;
          }
	  if (ij_part >= 1.0) {
	    ij_dec = 1;
	  }
	  else {
            ij_dec = std::min(static_cast<int>(ceil(-log10(ij_part))), max_decimal_places);
	  }
	  ij_width += ij_dec;
          column_decimals[i] = std::max(column_decimals[i], ij_dec);
	}
	break;
      case NumberFormat::SCIENTIFIC:
	break;
      case NumberFormat::INTEGER:
      case NumberFormat::LONG_LONG_INTEGER:
      case NumberFormat::UNSIGNED_INTEGER:
      case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
	if (ij_number >= 1.0001) {
	  ij_width += static_cast<int>(ceil(log10(ij_number)));
	}
	else {
	  ij_width += 1;
	}
	break;
      case NumberFormat::CHAR4:
	ij_width = 4;
      }
      column_widths[i] = std::max(column_widths[i], ij_width);
    }
  }

  // Expand the numerical widths if required by long words in the headings
  std::vector<std::vector<std::string>> heading_words;
  for (int i = 0; i < n_cols; i++) {
    heading_words.push_back(separateText(headings[i]));
    const int word_count = heading_words[i].size();
    int longest_word_length = 0;
    for (int j = 0; j < word_count; j++) {
      const int extend = (j < word_count - 1 && heading_words[i][j + 1].size() == 1) ? 2 : 0;
      longest_word_length = std::max(static_cast<int>(heading_words[i][j].size() + extend),
                                     longest_word_length);
      if (extend > 0) {
	j++;
      }
    }
    column_widths[i] = std::max(column_widths[i], longest_word_length);
  }

  // Print the table header
  switch (borders) {
  case BorderFormat::NONE:
  case BorderFormat::LIGHT:
    break;
  case BorderFormat::FULL:
    tableHorizontalRule(column_widths, borders, indentation, column_spacing, have_categories);
    break;
  }
  tableHeader(heading_words, column_widths, borders, indentation, column_spacing, have_categories);
  tableHorizontalRule(column_widths, borders, indentation, column_spacing, have_categories);

  // Print the table data
  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < indentation; j++) {
      printf(" ");
    }
    if (borders == BorderFormat::FULL) {
      printf("| ");
    }
    if (have_categories) {
      printf("%-*s", column_widths[0], categories[i].c_str());
    }
    for (int j = have_categories; j < n_cols; j++) {
      const int data_idx = j - have_categories;
      if (have_categories && j == 1) {
	switch (borders) {
	case BorderFormat::NONE:
	  for (int k = 0; k < column_spacing; k++) {
	    printf(" ");
	  }
	  break;
	case BorderFormat::LIGHT:
	case BorderFormat::FULL:
          printf(" | ");
          break;
	}       
      }
      else {
	for (int k = 0; k < column_spacing; k++) {
          printf(" ");
	}
      }
      switch(column_formats[data_idx]) {
      case NumberFormat::STANDARD_REAL:
	printf("%*.*lf", column_widths[j], column_decimals[j], column_data[data_idx][i].d);
	break;
      case NumberFormat::SCIENTIFIC:
	printf("%*.*e", column_widths[j], column_decimals[j], column_data[data_idx][i].d);
	break;
      case NumberFormat::INTEGER:
	printf("%*d", column_widths[j], column_data[data_idx][i].i);
	break;
      case NumberFormat::LONG_LONG_INTEGER:
	printf("%*lld", column_widths[j], column_data[data_idx][i].lli);
	break;
      case NumberFormat::UNSIGNED_INTEGER:
	printf("%*u", column_widths[j], column_data[data_idx][i].ui);
	break;
      case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
	printf("%*llu", column_widths[j], column_data[data_idx][i].ulli);
	break;
      case NumberFormat::CHAR4:
	for (int k = 0; k < column_widths[j] - 4; k++) {
	  printf(" ");
	}
	printf("%c%c%c%c", column_data[data_idx][i].c4.x, column_data[data_idx][i].c4.y,
	       column_data[data_idx][i].c4.z, column_data[data_idx][i].c4.w);
	break;
      }
    }
    switch (borders) {
    case BorderFormat::FULL:
      printf(" |\n");
      break;
    case BorderFormat::LIGHT:
    case BorderFormat::NONE:
      printf("\n");
      break;
    }
  }

  // One final horizontal rule if necessary
  switch (borders) {
  case BorderFormat::NONE:
  case BorderFormat::LIGHT:
    break;
  case BorderFormat::FULL:
    tableHorizontalRule(column_widths, borders, indentation, column_spacing, have_categories);
    break;
  }
}

} // namespace parse
} // namespace stormm
