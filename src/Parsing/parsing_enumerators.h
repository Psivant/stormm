// -*-c++-*-
#ifndef STORMM_PARSING_ENUMERATORS_H
#define STORMM_PARSING_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace parse {

/// \brief An enumerator for all types of numbers that might be encountered in an ascii text file.
///        This goes hand-in-hand with the PolyNumeric union below.
enum class NumberFormat {
  SCIENTIFIC,                  ///< Exponential, base ten notation
  STANDARD_REAL,               ///< General accounting notation
  INTEGER,                     ///< Print signed integers
  LONG_LONG_INTEGER,           ///< Print very large signed integers
  UNSIGNED_INTEGER,            ///< Print unsigned integers
  UNSIGNED_LONG_LONG_INTEGER,  ///< Print very large unsigned integers
  CHAR4                        ///< Print four-tuples of ASCII test characters
};
  
/// \brief Specify the way in which to print a standard real number or integer (with or without
///        leading zeros)
enum class NumberPrintStyle {
  STANDARD,      ///< Print without leading zeros
  LEADING_ZEROS  ///< Print with leading zeros
};

/// \brief Enumerate the types of wildcard characters.  The exact nature of the characters,
///        i.e. "*" or ".", may be defined elsewhere.
enum class WildCardKind {
  NONE,            ///< Not a wildcard: direct match is the only way to success
  FREE_CHARACTER,  ///< Matches any one character, but must match a character
  FREE_STRETCH     ///< Matches any number of consecutive characters, including no characters
};

/// \brief Specify the manner in which text is formatter to fit within the length of string
///        containing it.
enum class JustifyText {
  LEFT,    ///< The left side
  CENTER,  ///< Produces centered text
  RIGHT    ///< The right side
};

/// \brief List the possible types of header lines for a formatted table.  This class and the
///        functions below output strings that can later be fed into any file or output stream.
enum class TableHeadingLine {
  HEADER,          ///< A line containing the column headers (multiline headers are presently
                   ///<   formatted on separate lines by separate calls to buildTableHeader)
  HORIZONTAL_RULE  ///< A horiizontale rule to underline all table categories
};

/// \brief List the different border styles for decorating ASCII-text tables.
enum class BorderFormat {
  NONE,   ///< No borders or horizontal rules will accompany the table
  LIGHT,  ///< A horizontal rule will be printed to underline each of the category headings and a
          ///<   vertical rule will separate the leftmost column from those to the right
  FULL    ///< The table will be completely surrounded by borders, all columns will be separated
          ///<   by vertical rules, and a horizontal rule will separate the header from the table
          ///<   body.
};

/// \brief Many searches in TextFile objects will begin at a particular line.  If the query is not
///        found, it may be necessary to wrap the search back to the beginning and continue until
///        the original starting line.  This will indicate whether to do that.
enum class WrapTextSearch {
  NO, YES
};

/// \brief Differentiate between text data originating on a disk and in RAM.
enum class TextOrigin {
  DISK, RAM
};

/// \brief Produce a human-readable string corresponding to each enumeration.  Various overloads
///        of this function in this and other libraries and namespaces handle different
///        enumerators.
///
/// \param input  The enumeration of interest
/// \{
std::string getEnumerationName(NumberFormat input);
std::string getEnumerationName(NumberPrintStyle input);
std::string getEnumerationName(WildCardKind input);
std::string getEnumerationName(JustifyText input);
std::string getEnumerationName(TableHeadingLine input);
std::string getEnumerationName(BorderFormat input);
std::string getEnumerationName(WrapTextSearch input);
std::string getEnumerationName(TextOrigin input);
/// \}

} // namespace parse
} // namespace stormm

#endif
