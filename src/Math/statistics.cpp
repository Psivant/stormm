#include <cmath>
#include "copyright.h"
#include "statistics.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
double running_variance(const double sum_of_squares, const double sum_of_values,
                        const int sample_count, const VarianceMethod method) {
  
  // Catch bad inputs
  if (sample_count == 1) {
    switch (method) {
    case VarianceMethod::VARIANCE:
      return sum_of_squares;
    case VarianceMethod::ROOT_MEAN_SQUARED_DEVIATION:
    case VarianceMethod::NORMALIZED_RMSD:
      return sqrt(sum_of_squares);
    case VarianceMethod::STANDARD_DEVIATION:
    case VarianceMethod::COEFFICIENT_OF_VARIATION:
      return 0.0;
    }
  }

  const double ds_count = static_cast<double>(sample_count);
  const double variance = ds_count * sum_of_squares - sum_of_values * sum_of_values;
  switch (method) {
  case VarianceMethod::VARIANCE:
    return variance;
  case VarianceMethod::ROOT_MEAN_SQUARED_DEVIATION:
    return sqrt(variance) / ds_count;
  case VarianceMethod::NORMALIZED_RMSD:
    return sqrt(variance) * sum_of_values / (ds_count * ds_count);
  case VarianceMethod::STANDARD_DEVIATION:
    return sqrt(variance / (ds_count * (ds_count - 1.0)));
  case VarianceMethod::COEFFICIENT_OF_VARIATION:
    return sqrt(variance / (ds_count * (ds_count - 1.0))) * sum_of_values / ds_count;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double running_stdev(const double sum_of_squares, const double sum_of_values,
                     const int sample_count) {
  return running_variance(sum_of_squares, sum_of_values, sample_count,
                          VarianceMethod::STANDARD_DEVIATION);
}

//-------------------------------------------------------------------------------------------------
double running_rmsd(const double sum_of_squares, const double sum_of_values,
                    const int sample_count) {
  return running_variance(sum_of_squares, sum_of_values, sample_count,
                          VarianceMethod::ROOT_MEAN_SQUARED_DEVIATION);
}

//-------------------------------------------------------------------------------------------------
double running_coefficient_of_variation(const double sum_of_squares, const double sum_of_values,
                                        const int sample_count) {
  return running_variance(sum_of_squares, sum_of_values, sample_count,
                          VarianceMethod::COEFFICIENT_OF_VARIATION);
}

//-------------------------------------------------------------------------------------------------
double running_normalized_rmsd(const double sum_of_squares, const double sum_of_values,
                               const int sample_count) {
  return running_variance(sum_of_squares, sum_of_values, sample_count,
                          VarianceMethod::NORMALIZED_RMSD);
}

} // namepsace math
} // namespace stormm
