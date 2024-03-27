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

/// \brief Convert the enumerations above into human-readable strings.  Various overloads server
///        each enumerator.
/// \{
std::string getEnumerationName(PrintSituation input);
std::string getEnumerationName(DataFormat input);
/// \}

} // namespace diskutil
} // namespace stormm

#endif
