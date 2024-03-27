// -*-c++-*-
#ifndef STORMM_ASCII_NUMBERS_H
#define STORMM_ASCII_NUMBERS_H

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "copyright.h"
#include "polynumeric.h"
#include "parse.h"

namespace stormm {
namespace parse {

/// \brief Print a series of formatted numbers for a fixed-column output.  This function will
///        be faster than a typical fprintf() call for general format real numbers, not least
///        because it vectorizes the printing.  Protections are provided for extreme cases that
///        would break a typical 32-bit integer format.
///
/// Overloaded:
///   - Print real numbers (double input, two integer format specifiers)
///   - Print integers (integer input, one format specifier)
///
/// \param foutp            Output file stream
/// \param values           The values to print
/// \param values_per_line  The number of values to print to any given line
/// \param width            Total width of an individual number's format
/// \param decimal          The number of decimal places to print
/// \param format           Format of numbers to print
/// \param caller           Name of the calling function
/// \param task             Description of the task being performed, for error reporting
/// \{
void printNumberSeries(std::ofstream *foutp, const std::vector<PolyNumeric> &values,
                       const int values_per_line, const int width, const int decimal,
                       const NumberFormat format, const std::string &caller = std::string(""),
                       const std::string &task = std::string(""));

void printNumberSeries(std::ofstream *foutp, const std::vector<PolyNumeric> &values,
                       const int values_per_line, const int width, const NumberFormat format,
                       const std::string &caller = std::string(""),
                       const std::string &task = std::string(""));
/// \}

/// \brief Read a series of numbers from a TextFile object.  The numbers are expected to be found
///        in a fixed-column format described by the input arguments.
///
/// \param tf               Text file, committed to RAM, containing the data to read
/// \param n_values         The total number of values to read (further values will be ignored)
/// \param values_per_line  The number of values to print to any given line
/// \param width            Total width of an individual number's format
/// \param decimal          The number of decimal places to print
/// \param format           Format of numbers to read
/// \param caller           Name of the calling function
/// \param task             Description of the task being performed, for error reporting
std::vector<PolyNumeric> readNumberSeries(const TextFile &tf, const int start_line,
                                          const int n_values, const int values_per_line,
                                          const int width, const int decimal,
                                          const NumberFormat format,
                                          const std::string &caller = std::string(""),
                                          const std::string &task = std::string(""));

/// \brief Extract an integer number from a character string.  This is not intended for the
///        highest-performance applications and should be used in contexts where only a few such
///        integers must be read in specialized formats.
///
/// Overloaded:
///   - Accept a C-style character string
///   - Accept a C++ string object
///   - Accept a TextFile object
///
/// \param input_text     Text string containing
/// \param line_index     The line of a text file on which to find the number (only used when
///                       parsing a TextFile object)
/// \param start_index    The starting index for reading the formatted number
/// \param number_length  Number of characters to read and take as the numerical value
/// \{
llint readIntegerValue(const char* number_text, int start_index, int number_length);

llint readIntegerValue(const std::string &number_text, int start_index, int number_length);

llint readIntegerValue(const TextFile &tf, int line_idex, int start_index, int number_length);
/// \}

/// \brief Read a real-valued number from a character string.  This is not intended for the
///        highest-performance applications and should be used in contexts where only a few such
///        integers must be read in specialized formats.  Overloading and parameters for this
///        function follow from readIntegerValue() above.
/// \{
double readRealValue(const char* number_text, int start_index, int number_length);

double readRealValue(const std::string &number_text, int start_index, int number_length);

double readRealValue(const TextFile &tf, int line_idex, int start_index, int number_length);
/// \}

} // namespace parse
} // namespace stormm

#endif
