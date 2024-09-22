// -*-c++-*-
#ifndef STORMM_STATISTICAL_ENUMERATORS_H
#define STORMM_STATISTICAL_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace stmath {

/// \brief Enumerate various statistical methods for computing the variance of a sequence of
///        n numbers (i.e. as presented in a std::vector).
enum class VarianceMethod {
  VARIANCE,                    ///< Basic variance, sum_i((x_i - <x>)^2)
  STANDARD_DEVIATION,          ///< Normalized root variance, sqrt(sum_i((x_i - <x>)^2) / (n - 1))
  ROOT_MEAN_SQUARED_DEVIATION, ///< Unnormalized root variance, sqrt(sum_i((x_i - <x>)^2) / n)
  COEFFICIENT_OF_VARIATION,    ///< Standard deviation divided by mean absolute value <|x|>
  NORMALIZED_RMSD              ///< Root mean squared deviation divided by mean absolute
                               ///<   value <|x|>
};

/// \brief Enumerate possible arrangements of data that might be encountered.
enum class DataOrder {
  ASCENDING,   ///< The data is found in ascending order, element i always <= element i + 1
  DESCENDING,  ///< The data is found in descending order, element i always >= element i + 1
  NONE         ///< The data is found in no discernible order
};

/// \brief Enumerate resampling strategies
enum class ResamplingMethod {
  JACKKNIFE,      ///< Random resampling, without replacement
  BOOTSTRAP       ///< Random resampling, without replacement
};

/// \brief Produce human-readable strings corresponding to each of the enumerations above.  Various
///        overloads in this and other libraries serve each enum class.
///
/// \param input  The enumeration to translate
/// \{
std::string getEnumerationName(VarianceMethod input);
std::string getEnumerationName(DataOrder input);
std::string getEnumerationName(ResamplingMethod input);
/// \}

} // namespace stmath
} // namespace stormm

#endif
