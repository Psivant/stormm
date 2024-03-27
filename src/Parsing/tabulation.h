// -*-c++-*-
#ifndef STORMM_TABULATION_H
#define STORMM_TABULATION_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Math/sorting_enumerators.h"
#include "parsing_enumerators.h"
#include "polynumeric.h"

namespace stormm {
namespace parse {

using stmath::SortDirection;

/// \brief Print a horizontal rule for this table.
///
/// \param borders                Style of borders to print around the table contents
/// \param indentation            Indentation for the table as a whole
/// \param column_spacing         Number of spaces to place between data columns (a special
///                               delimiter is placed between the category and data columns if
///                               this
/// \param have_categories        Flag to indicate that a categories column exists on the leftmost
///                               side of the table  
void tableHorizontalRule(const std::vector<int> &column_widths, BorderFormat borders,
                         int indentation, int column_spacing, bool have_categories);

/// \brief Print the headings for this table.
///
/// \param borders                Style of borders to print around the table contents
/// \param indentation            Indentation for the table as a whole
/// \param column_spacing         Number of spaces to place between data columns (a special
///                               delimiter is placed between the category and data columns if
///                               this
/// \param have_categories        Flag to indicate that a categories column exists on the leftmost
///                               side of the table
void tableHeader(const std::vector<std::vector<std::string>> &heading_words,
                 const std::vector<int> &column_widths, BorderFormat borders, int indentation,
                 int column_spacing, bool have_categories);

/// \brief Print a formatted table with a flexible assortment of data.  Each column of the data
///        must have as many elements as there are categories.  Sorting by one column is feasible.
///
/// \param headings               Names of each column (including the categoies column itself)
/// \param categories             Names of each row
/// \param column_data            Two-dimensional array of PolyNumeric multi-format numbers
/// \param column_formats         Formats for all number vectors in column_data.  This indiates
///                               how to read the PolyNumeric unions.
/// \param column_sort            Name of the column to sort (must match a column header)
/// \param borders                Style of borders to print around the table contents
/// \param indentation            Indentation for the table as a whole
/// \param column_spacing         Number of spaces to place between data columns (a special
///                               delimiter is placed between the category and data columns if
///                               this
/// \param max_decimal_places_in  The maximum number of digits after the decimal.  This may be
///                               reduced given the nature of the data.
/// \param precision              Results for real numbers will be reported to a precision of at
///                               least one part in this many.
/// \param sorting_direction      Order in which to arrange data based on the sorting column (see
///                               column_sort above).  Default is "AUTOMATIC", which orders
///                               categories in alphabetical order on successive lines (if the
///                               first column is selected), or numerical columns in decreasing
///                               order on successive lines.
void printTable(const std::vector<std::string> &headings,
                const std::vector<std::string> &categories,
		const std::vector<std::vector<PolyNumeric>> &column_data,
                const std::vector<NumberFormat> &column_formats,
                const std::string &column_sort = std::string(""),
                BorderFormat borders = BorderFormat::LIGHT, int indentation = 1,
                int column_spacing = 2, int max_decimal_places_in = 4, double precision = 1.0e6,
                SortDirection sorting_direction = SortDirection::AUTOMATIC);

} // namspace parse
} // namespace stormm

#endif
