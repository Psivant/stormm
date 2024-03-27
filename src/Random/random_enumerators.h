// -*-c++-*-
#ifndef STORMM_RANDOM_ENUMERATORS_H
#define STORMM_RANDOM_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace random {

/// \brief Enumerate the types of random numbers that can be generated.
enum class RandomNumberKind {
  UNIFORM,   ///< Uniform random number distribution
  GAUSSIAN   ///< Normal distribution of random numbers
};

/// \brief List the random number generator types that can power a RandomNumberMill object.
enum class RandomAlgorithm {
  XOROSHIRO_128P,  ///< Xoroshiro128+ generator (see random.h, fails BigCrush and not advised for
                   ///<   powering mills with > 1024 generator streams)
  XOSHIRO_256PP    ///< Xoshiro256++ generator (see random.h--high quality generator)
};

/// \brief Define the order in which various random number generators will fill a matrix.
enum class RngFillMode {
  COLUMNS,  ///< Fill the column-major matrix one column at a time
  ROWS      ///< Fill the column-major matrix in row-major order, one row at a time
};


/// \brief Produce human-readable strings for each enumeration.  Overloads of this function here
///        and in other libraries server different enum classes.
///
/// \param input  The enumerated value to translate
/// \{
std::string getEnumerationName(RandomNumberKind input);
std::string getEnumerationName(RandomAlgorithm input);
std::string getEnumerationName(RngFillMode input);
/// \}

} // namespace random
} // namespace stormm

#endif
