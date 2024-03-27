// -*-c++-*-
#ifndef STORMM_UNIT_TEST_ENUMERATORS_H
#define STORMM_UNIT_TEST_ENUMERATORS_H

#include "copyright.h"

namespace stormm {
namespace testing {

/// \brief Enumerate the types of comparisons to make between two numbers
enum class ComparisonType {
  ABSOLUTE,            ///< Maximum tolerated absolute deviation from some target value for a
                       ///<   scalar, or absolute deviation of the largest outlier for a
                       ///<   vector of numbers
  RELATIVE,            ///< Maximum tolerated relative deviation from some target value for a
                       ///<   scalar, or relative deviation of the largest outlier for a
                       ///<   vector of numbers
  MEAN_UNSIGNED_ERROR, ///< Maximum tolerated mean unsigned error for a vector of numbers
                       ///<   (works as ABSOLUTE if attempted with a scalar)
  RELATIVE_RMS_ERROR   ///< Akin to the coefficient of variation (relative standard deviation) in
                       ///<   the limit of large sequences, this computes the root mean squared
                       ///<   error as a function of the unsigned mean value of the reference data
};

/// \brief Enumerate the possible results of the function Check
enum class CheckResult {
  SUCCESS, ///< Test passes
  SKIPPED, ///< Test was skipped
  IGNORED, ///< Test failure was ignored
  FAILURE  ///< Test fails
};

/// \brief Enumerate priorities by which a test shall run, which may modulate the results
enum class TestPriority {
  CRITICAL,     ///< The test is assumed to be feasible and must be run, and will return SUCCESS or
                ///<   FAILURE when it does run
  NON_CRITICAL, ///< The test is not essential and will return SUCCESS if it passes but IGNORED
                ///<   rather than FAILURE if it fails
  ABORT,        ///< The test should be aborted and will return SKIPPED when run
};

/// \brief Enmerate the relational and comparison operators
enum class RelationalOperator {
  EQUAL, EQ, NOT_EQUAL, NE,
  GREATER_THAN, GT, LESS_THAN, LT,
  GREATER_THAN_OR_EQUAL, GE, LESS_THAN_OR_EQUAL, LE
};

/// \brief Enumerate different levels of verbosity
enum class TestVerbosity {
  FULL,         ///< Print all available messages for the user
  COMPACT,      ///< Print a subset of information for the user
  FAILURE_ONLY  ///< Only alert the user of failures
};

/// \brief Different snapshotting operations that can be performed for checking a series of
///        numbers (in a std::vector) against reference (data in a file)
enum class SnapshotOperation {
  COMPARE,  ///< See if the data in memory matches the file on record
  SNAPSHOT  ///< Write data in memory to a file on record
};

/// \brief Describe whether the temporary directory is required for testing.
enum class TmpdirStatus {
  NOT_REQUIRED,
  REQUIRED
};

/// \brief Get a string corresponding to a particular test priority enumeration.  Various overloads
///        of this function in this and other libraries and namespaces serve different enumerators.
///
/// \param input  The enumeration to translate
/// \{
std::string getEnumerationName(ComparisonType input);
std::string getEnumerationName(CheckResult input);
std::string getEnumerationName(TestPriority input);
std::string getEnumerationName(RelationalOperator input);
std::string getEnumerationName(TestVerbosity input);
std::string getEnumerationName(SnapshotOperation input);
std::string getEnumerationName(TmpdirStatus input);
/// \}

} // namespace testing
} // namespace stormm

#endif
