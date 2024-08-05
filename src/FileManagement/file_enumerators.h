// -*-c++-*-
#ifndef STORMM_FILE_ENUMERATORS_H
#define STORMM_FILE_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace diskutil {

/// \brief Enumerate the situations that might be encountered when writing a trajectory file
enum class PrintSituation {
  OPEN_NEW,   ///< Expect to open a new file, and only open a new file if no file by the given name
              ///<   already exists
  APPEND,     ///< Append an existing file, or open a new file if no file by the given name exists
  OVERWRITE,  ///< Always open a new file, overwriting any existing file by the same name
  UNKNOWN     ///< Option to let the program decide what printing behavior to use, in lieu of user
              ///<   or developer input
};

/// \brief Format of the file (just a fancy way of coding "binary" or "ascii")
enum class DataFormat {
  ASCII,   ///< Ascii, formatted and human-readable output
  BINARY   ///< Binary, direct-to-disk copy of numbers
};

/// \brief Various means of compressing trajectory or restart files.
enum class TrajectoryCompression {
  NONE,   ///< No compression is present
  POSIT   ///< A format involving some loss of information but which stores coordinates in an
          ///<   enhanced, trajectory-optimized format which maximizes the number of bits devoted
          ///<   to the mantissa.  In this format, arrays of unsigned integers must be interpreted
          ///<   according to a "key" byte with the number of mantissa bits in its low four bits
          ///<   and the shift of the exponent in the high four bits.  Further details of the
          ///<   implementation can be found in binary_encoding.tpp.  The POSIT format is a regular
          ///<   width format and will result in all trajectory snapshots having the same size.
};

/// \brief Convert the enumerations above into human-readable strings.  Various overloads serve
///        each enumerator.
/// \{
std::string getEnumerationName(PrintSituation input);
std::string getEnumerationName(DataFormat input);
/// \}

} // namespace diskutil
} // namespace stormm

#endif
