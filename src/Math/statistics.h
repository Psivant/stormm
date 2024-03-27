// -*-c++-*-
#ifndef STORMM_STATISTICS_H
#define STORMM_STATISTICS_H

#include "copyright.h"
#include "statistical_enumerators.h"

namespace stormm {
namespace stmath {

/// \brief Compute the variance of a set of data as a running quantity based on just the power
///        series of the one- and two-moments of the data.  Different methods fed to this function
///        can also produce standard deviation, root mean squared deviation, and either quantity
///        normalized by the mean of the data.
///
/// \param sum_of_squares  The sum of the data's second moment
/// \param sum_of_values   The sum of the data's first moment
/// \param sample_count    The number of samples in the data
/// \param method          Statistical analysis method
double running_variance(double sum_of_squares, double sum_of_values, int sample_count,
                        VarianceMethod method = VarianceMethod::VARIANCE);

/// \brief Compute the standard deviation of a set of data based on its first two moments.  Returns
///        zero if there are less than two data points.
///
/// \param sum_of_squares  The sum of the data's second moment
/// \param sum_of_values   The sum of the data's first moment
/// \param sample_count    The number of samples in the data
double running_stdev(double sum_of_squares, double sum_of_values, int sample_count);

/// \brief Compute the root mean squared deviation of a set of data based on its first two moments.
///
/// \param sum_of_squares  The sum of the data's second moment
/// \param sum_of_values   The sum of the data's first moment
/// \param sample_count    The number of samples in the data
double running_rmsd(double sum_of_squares, double sum_of_values, int sample_count);

/// \brief Compute the mean-normalized standard deviation of a set of data based on its first two
///        moments.
///
/// \param sum_of_squares  The sum of the data's second moment
/// \param sum_of_values   The sum of the data's first moment
/// \param sample_count    The number of samples in the data
double running_coefficient_of_variation(double sum_of_squares, double sum_of_values,
                                        int sample_count);

/// \brief Compute the normalized root mean squared deviation of a set of data based on its first
///        two moments.
///
/// \param sum_of_squares  The sum of the data's second moment
/// \param sum_of_values   The sum of the data's first moment
/// \param sample_count    The number of samples in the data
double running_normalized_rmsd(double sum_of_squares, double sum_of_values, int sample_count);

} // namespace stmath
} // namespace stormm

#endif
