// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace testing {

//-------------------------------------------------------------------------------------------------
template <typename T>
CheckResult check(const std::vector<T> &lhs, const RelationalOperator relationship,
                  const Approx &rhs, const std::string &error_message,
                  const TestPriority urgency) {

  using parse::free_number_format;

  // Abort immediately if this test is not runnable
  if (urgency == TestPriority::ABORT) {
    return check(false, error_message, urgency);
  }
  
  // Inequality vector comparisons are only valid for absolute relationships
  switch (relationship) {
  case RelationalOperator::EQUAL:
  case RelationalOperator::EQ:
  case RelationalOperator::NOT_EQUAL:
  case RelationalOperator::NE:
    break;
  case RelationalOperator::GREATER_THAN:
  case RelationalOperator::GT:
  case RelationalOperator::LESS_THAN:
  case RelationalOperator::LT:
  case RelationalOperator::GREATER_THAN_OR_EQUAL:
  case RelationalOperator::GE:
  case RelationalOperator::LESS_THAN_OR_EQUAL:
  case RelationalOperator::LE:
    switch (rhs.getStyle()) {
    case ComparisonType::ABSOLUTE:
      break;
    case ComparisonType::RELATIVE:
    case ComparisonType::MEAN_UNSIGNED_ERROR:
    case ComparisonType::RELATIVE_RMS_ERROR:
      rtErr("Vectorized comparisons are invalid for non-absolute " +
            getRelationalOperatorString(relationship) + "relationships.", "check");
    }
    break;
  }

  // Enumerate the possible relations between the two values
  switch (relationship) {
  case RelationalOperator::EQUAL:
  case RelationalOperator::EQ:
    if (lhs == rhs) {
      gbl_test_results.logResult(CheckResult::SUCCESS);
      return CheckResult::SUCCESS;
    }
    break;
  case RelationalOperator::NOT_EQUAL:
  case RelationalOperator::NE:
    if (lhs != rhs) {
      gbl_test_results.logResult(CheckResult::SUCCESS);
      return CheckResult::SUCCESS;
    }
    break;
  case RelationalOperator::GREATER_THAN:
  case RelationalOperator::GT:
    if (lhs > rhs) {
      gbl_test_results.logResult(CheckResult::SUCCESS);
      return CheckResult::SUCCESS;
    }
    break;
  case RelationalOperator::LESS_THAN:
  case RelationalOperator::LT:
    if (lhs < rhs) {
      gbl_test_results.logResult(CheckResult::SUCCESS);
      return CheckResult::SUCCESS;
    }
    break;
  case RelationalOperator::GREATER_THAN_OR_EQUAL:
  case RelationalOperator::GE:
    if (lhs >= rhs) {
      gbl_test_results.logResult(CheckResult::SUCCESS);
      return CheckResult::SUCCESS;
    }
    break;
  case RelationalOperator::LESS_THAN_OR_EQUAL:
  case RelationalOperator::LE:
    if (lhs <= rhs) {
      gbl_test_results.logResult(CheckResult::SUCCESS);
      return CheckResult::SUCCESS;
    }
    break;
  }
  
  // If this point is reached, the check must have failed.  Investigate why and report based on
  // the type of input vector and approximate comparison being performed.  The right hand side
  // vector WILL be of a scalar type and the same length as the left hand side vector, as the
  // above comparisons would have led to compiler errors were it not.  Different results may be
  // desirable if the vector is integral, unsigned integral, or real.
  std::string error_edit("  ");
  const std::string vofsz = "comparison of vectors of size " + std::to_string(lhs.size());
  std::string reasoning;
  switch (relationship) {
  case RelationalOperator::EQUAL:
  case RelationalOperator::EQ:
    switch (rhs.getStyle()) {
    case ComparisonType::ABSOLUTE:
      reasoning = " fails with a maximum deviation of ";
      break;
    case ComparisonType::RELATIVE:
      reasoning = " fails with a maximum relative difference of ";
      break;
    case ComparisonType::MEAN_UNSIGNED_ERROR:
      reasoning = " fails due to a mean unsigned differences of ";
      break;
    case ComparisonType::RELATIVE_RMS_ERROR:
      reasoning = " fails due to a relative root mean squared error of ";
      break;
    }
    break;
  case RelationalOperator::NOT_EQUAL:
  case RelationalOperator::NE:
    switch (rhs.getStyle()) {
    case ComparisonType::ABSOLUTE:
      reasoning = " was expected to fail but carries a maximum deviation of ";
      break;
    case ComparisonType::RELATIVE:
      reasoning = " was expected to fail but carries a maximum relative difference of ";
      break;
    case ComparisonType::MEAN_UNSIGNED_ERROR:
      reasoning = " was expected to fail but carries a mean unsigned differences of ";
      break;
    case ComparisonType::RELATIVE_RMS_ERROR:
      reasoning = " was expected to fail but carries a relative root mean squared error of ";
      break;
    }
    break;
  case RelationalOperator::GREATER_THAN:
  case RelationalOperator::GT:
  case RelationalOperator::LESS_THAN:
  case RelationalOperator::LT:
  case RelationalOperator::GREATER_THAN_OR_EQUAL:
  case RelationalOperator::GE:
  case RelationalOperator::LESS_THAN_OR_EQUAL:
  case RelationalOperator::LE:

    // Only the absolute comparison case will have made it through traps
    reasoning = " fails to satisfy the inequality with a maximum crossover of ";
    break;
  }

  // Check that the vector comparison is valid.  Trap this case to protect operations later on
  // in this function that might require their inputs to be the same length.
  if (verifyVectorApproxCompatibility(lhs, rhs) == false) {
    error_edit += "The vectors are of different lengths (" + std::to_string(lhs.size()) + " and " +
                  std::to_string(rhs.size()) + ").";
    return check(false, error_message + error_edit, urgency);
  }

  // Pre-compute the relevant deviations
  double max_deviation, mue, rel_deviation, rel_rmse;
  const std::vector<double> dlhs(lhs.begin(), lhs.end());
  const size_t npts = lhs.size();
  switch (relationship) {
  case RelationalOperator::EQUAL:
  case RelationalOperator::EQ:
  case RelationalOperator::NOT_EQUAL:
  case RelationalOperator::NE:
    switch (rhs.getStyle()) {
    case ComparisonType::ABSOLUTE:
    case ComparisonType::MEAN_UNSIGNED_ERROR:
      max_deviation = fabs(maxAbsoluteDifference<double>(dlhs, rhs.getValues()));
      mue = meanUnsignedError<double>(dlhs, rhs.getValues());
      break;
    case ComparisonType::RELATIVE:
    case ComparisonType::RELATIVE_RMS_ERROR:
      rel_deviation = fabs(maxRelativeDifference<double>(dlhs, rhs.getValues()));
      rel_rmse = relativeRmsError<double>(dlhs, rhs.getValues()) * 100.0;
      break;
    }
    break;
  case RelationalOperator::GREATER_THAN:
  case RelationalOperator::GT:
    {
      std::vector<double> drhs = rhs.getValues();
      const double rhs_margin = rhs.getMargin();
      for (size_t i = 0; i < npts; i++) {
        drhs[i] += rhs_margin;
      }
      max_deviation = 0.0;
      mue = 0.0;
      for (size_t i = 0; i < npts; i++) {
        if (dlhs[i] <= drhs[i]) {
          max_deviation = std::max(max_deviation, drhs[i] - dlhs[i]);
          mue += drhs[i] - dlhs[i];
        }
      }
      mue /= static_cast<double>(npts);
    }
    break;
  case RelationalOperator::LESS_THAN:
  case RelationalOperator::LT:
    {
      std::vector<double> drhs = rhs.getValues();
      const double rhs_margin = rhs.getMargin();
      for (size_t i = 0; i < npts; i++) {
        drhs[i] -= rhs_margin;
      }
      max_deviation = 0.0;
      mue = 0.0;
      for (size_t i = 0; i < npts; i++) {
        if (dlhs[i] >= drhs[i]) {
          max_deviation = std::max(max_deviation, dlhs[i] - drhs[i]);
          mue += dlhs[i] - drhs[i];
        }
      }
      mue /= static_cast<double>(npts);
    }
    break;
  case RelationalOperator::GREATER_THAN_OR_EQUAL:
  case RelationalOperator::GE:
    {
      std::vector<double> drhs = rhs.getValues();
      const double rhs_margin = rhs.getMargin();
      for (size_t i = 0; i < npts; i++) {
        drhs[i] -= rhs_margin;
      }
      max_deviation = 0.0;
      mue = 0.0;
      for (size_t i = 0; i < npts; i++) {
        if (dlhs[i] < drhs[i]) {
          max_deviation = std::max(max_deviation, drhs[i] - dlhs[i]);
          mue += drhs[i] - dlhs[i];
        }
      }
      mue /= static_cast<double>(npts);
    }
    break;
  case RelationalOperator::LESS_THAN_OR_EQUAL:
  case RelationalOperator::LE:
    {
      std::vector<double> drhs = rhs.getValues();
      const double rhs_margin = rhs.getMargin();
      for (size_t i = 0; i < npts; i++) {
        drhs[i] += rhs_margin;
      }
      max_deviation = 0.0;
      mue = 0.0;
      for (size_t i = 0; i < npts; i++) {
        if (dlhs[i] > drhs[i]) {
          max_deviation = std::max(max_deviation, dlhs[i] - drhs[i]);
          mue += dlhs[i] - drhs[i];
        }
      }
      mue /= static_cast<double>(npts);
    }
    break;
  }

  // Adjust settings based on the type of input
  int digits;
  NumberFormat nfmt;
  if (isSignedIntegralScalarType<T>()) {
    digits = 0;
    nfmt = NumberFormat::STANDARD_REAL;
  }
  else if (isUnsignedIntegralScalarType<T>()) {
    digits = 0;
    nfmt = NumberFormat::STANDARD_REAL;
  }
  else if (isFloatingPointScalarType<T>()) {
    digits = 4;
    nfmt = NumberFormat::SCIENTIFIC;
  }
  else {
    rtErr("Unknown input vector data type " + std::string(typeid(T).name()) + ".", "check");
  }

  // Print the appropriate message by recursively calling the simplest form of check()
  switch (rhs.getStyle()) {
  case ComparisonType::ABSOLUTE:
    switch (relationship) {
    case RelationalOperator::EQUAL:
    case RelationalOperator::EQ:
    case RelationalOperator::NOT_EQUAL:
    case RelationalOperator::NE:
      error_edit += "Absolute " + vofsz + reasoning +
                    realToString(max_deviation, free_number_format, digits, nfmt) +
                    " and mean unsigned error of " +
                    realToString(mue, 11, 4, NumberFormat::SCIENTIFIC) + ".  ";
      break;
    case RelationalOperator::GREATER_THAN:
    case RelationalOperator::GT:
    case RelationalOperator::LESS_THAN:
    case RelationalOperator::LT:
    case RelationalOperator::GREATER_THAN_OR_EQUAL:
    case RelationalOperator::GE:
    case RelationalOperator::LESS_THAN_OR_EQUAL:
    case RelationalOperator::LE:
      error_edit += "Absolute " + vofsz + reasoning +
                    realToString(max_deviation, free_number_format, digits, nfmt) +
                    " and mean crossover of " +
                    realToString(mue, 11, 4, NumberFormat::SCIENTIFIC) + ".  ";
      break;
    }
    break;
  case ComparisonType::RELATIVE:
    error_edit += "Relative " + vofsz + reasoning +
                  realToString(rel_deviation, free_number_format, digits, nfmt) +
                  " and relative root mean squared error of " +
                  realToString(rel_rmse, 11, 4, NumberFormat::SCIENTIFIC) + "%.  ";
    break;
  case ComparisonType::MEAN_UNSIGNED_ERROR:
    error_edit += "Absolute " + vofsz + reasoning +
                  realToString(mue, 11, 4, NumberFormat::SCIENTIFIC) +
                  " (maximum deviation " +
                  realToString(max_deviation, free_number_format, digits, nfmt) + ").  ";
    break;
  case ComparisonType::RELATIVE_RMS_ERROR:
    error_edit += "Relative " + vofsz + reasoning +
                  realToString(rel_rmse, 11, 4, NumberFormat::SCIENTIFIC) +
                  "% (maximum relative difference " +
                  realToString(rel_deviation, free_number_format, digits, nfmt) + ").  ";
    break;
  }

  // The alignment report must deal with the vectors in their double-precision representation, as
  // the Approx object has already converted the right hand side to that format.
  error_edit += vectorAlignmentReport(polyNumericVector(dlhs), polyNumericVector(rhs.getValues()),
                                      NumberFormat::STANDARD_REAL, rhs.getMargin());
  return check(false, error_message + error_edit, urgency);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CheckResult check(const Approx &lhs, const RelationalOperator relationship,
                  const std::vector<T> &rhs, const std::string &error_message,
                  const TestPriority urgency) {
  switch (relationship) {
  case RelationalOperator::EQUAL:
  case RelationalOperator::EQ:
  case RelationalOperator::NOT_EQUAL:
  case RelationalOperator::NE:
    return check(rhs, relationship, lhs, error_message, urgency);
  case RelationalOperator::GREATER_THAN:
  case RelationalOperator::GT:
    return check(rhs, RelationalOperator::LESS_THAN, lhs, error_message, urgency);
  case RelationalOperator::LESS_THAN:
  case RelationalOperator::LT:
    return check(rhs, RelationalOperator::GREATER_THAN, lhs, error_message, urgency);
  case RelationalOperator::GREATER_THAN_OR_EQUAL:
  case RelationalOperator::GE:
    return check(rhs, RelationalOperator::LESS_THAN_OR_EQUAL, lhs, error_message, urgency);
  case RelationalOperator::LESS_THAN_OR_EQUAL:
  case RelationalOperator::LE:
    return check(rhs, RelationalOperator::GREATER_THAN_OR_EQUAL, lhs, error_message, urgency);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T1, typename T2>
CheckResult check(const std::vector<T1> &lhs, const RelationalOperator relationship,
                  const std::vector<T2> &rhs, const std::string &error_message,
                  const TestPriority urgency) {
  return check(lhs, relationship, Approx(rhs), error_message, urgency);
}

} // namespace testing
} // namespace stormm
