// -*-c++-*-
#ifndef STORMM_UNIT_TEST_H
#define STORMM_UNIT_TEST_H

#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>
#include "copyright.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Parsing/polynumeric.h"
#include "Parsing/textfile.h"
#include "Reporting/error_format.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "approx.h"
#include "checklist.h"
#include "file_snapshot.h"
#include "test_environment.h"
#include "vector_report.h"

/// \brief Check that some action throws an error
///
/// \param test_code      Any amount of code that the developer wishes to test.  This will get a
///                       semicolon (;) appended to it and wrapped with a try { } scope.
/// \param error_message  The error message to print if the test code passes without throwing an
///                       exception. (This must be a string.)
#define CHECK_THROWS(test_code, error_message) {\
  bool exception_thrown = false;\
  try {\
    test_code;\
  }\
  catch (std::runtime_error) {\
    exception_thrown = true;\
  }\
  check(exception_thrown, error_message);\
}

/// \brief Check that some action throws an error
///
/// \param test_code      Any amount of code that the developer wishes to test.  This will get a
///                       semicolon (;) appended to it and wrapped with a try { } scope.
/// \param error_message  The error message to print if the test code passes without throwing an
///                       exception. (This must be a string.)
/// \param priority       The test priority.  This will pass on instructions to ignore a failed
///                       test, or skip a test altogether if it is not going to be feasible.
///                       (This must be a TestPriority enumeration type value.)
#define CHECK_THROWS_SOFT(test_code, error_message, priority) {\
  if (priority == stormm::testing::TestPriority::ABORT) {\
    check(true, error_message, priority);\
  }\
  else {\
    bool exception_thrown = false;\
    try {\
      test_code;\
    }\
    catch (std::runtime_error) {\
      exception_thrown = true;\
    }\
    check(exception_thrown, error_message, priority);\
  }\
}

namespace stormm {
namespace testing {

using errors::rtErr;
using errors::terminalFormat;
using data_types::isScalarType;
using data_types::isSignedIntegralScalarType;
using data_types::isUnsignedIntegralScalarType;
using data_types::isFloatingPointScalarType;
using stmath::sum;
using stmath::maxAbsoluteDifference;
using stmath::maxRelativeDifference;
using stmath::meanUnsignedError;
using stmath::relativeRmsError;
using parse::NumberFormat;
using parse::polyNumericVector;
using parse::TextFile;

/// \brief Free function for changing the section of the global test results object.
///
/// Overloaded:
///   - Jump to a section by name
///   - Jump to a section by index
///
/// \param section_name  The name of the section to add or jump to
/// \{
void section(const std::string &section_name);
void section(const int section_index);
/// \}

/// \brief Get the index of a section from the global test results, based on its name.
///
/// \param name  The name of the section of interest
int unitTestSectionIndex(const std::string &name);

/// \brief Get the name of a section from the global test results, based on its index.
///
/// \param index  The index of the section of interest
std::string unitTestSectionName(int index);

/// \brief Check that some statement is true, and report a warning or error message if the
///        statement is not.  As with the Approx object, the naming of this function is designed
///        to be reminiscent of Catch2.
///
/// Overloaded:
///   - Basic form including a boolean expression and an error message to print in the event that
///     the expression evluates to false
///   - General form including left hand side, comparison (relational) operator, and right hand
///     side with a message to be enhanced by the evaluations of each side
///   - Approximate comparison with an input value (the input value is expected on the left hand
///     side, with the approximate quality to the right of the relational operator, but if these
///     positions are reversed then the overloads cascade into the real-valued left vs. Approx
///     right scenario, with reversal of the relational operator as needed).
///   - Direct comparison between two strings or char4 tuples
///   - Direct comparison between vectors of char4 tuples
///
/// \param statement      The simplest, pre-evaluated boolean statement
/// \param error_message  The error message to print if the statement is false
/// \param lhs            Left hand side of some relation
/// \param relationship   Relational operator (in most cases ==, "is equal to")
/// \param rhs            Right hand side of some relation
/// \{
CheckResult check(const bool statement, const std::string &error_message = "",
                  TestPriority urgency = TestPriority::CRITICAL);

CheckResult check(double lhs, RelationalOperator relationship, const Approx &rhs,
                  const std::string &error_message, TestPriority urgency = TestPriority::CRITICAL);

CheckResult check(const Approx &lhs, RelationalOperator relationship, const double rhs,
                  const std::string &error_message, TestPriority urgency = TestPriority::CRITICAL);

CheckResult check(const double lhs, RelationalOperator relationship, const double rhs,
                  const std::string &error_message, TestPriority urgency = TestPriority::CRITICAL);

CheckResult check(const std::string &lhs, RelationalOperator relationship,
                  const std::string &rhs, const std::string &error_message,
                  TestPriority urgency = TestPriority::CRITICAL);

CheckResult check(const char4 &lhs, RelationalOperator relationship, const char4 &rhs,
                  const std::string &error_message, TestPriority urgency = TestPriority::CRITICAL);

CheckResult check(const std::vector<char4> &lhs, RelationalOperator relationship,
                  const std::vector<char4> &rhs, const std::string &error_message,
                  TestPriority urgency = TestPriority::CRITICAL);
/// \}

/// \brief Comparison of two vectors when the approximate object is listed second.
///
/// \param lhs            The vector of input values
/// \param relationship   The relation between the values, i.e. == or !=
/// \param rhs            Vector of reference values with a tolerance and comparison style
/// \param error_message  Error message to print
template <typename T>
CheckResult check(const std::vector<T> &lhs, RelationalOperator relationship,
                  const Approx &rhs, const std::string &error_message,
                  TestPriority urgency = TestPriority::CRITICAL);

/// \brief Comparison of two vectors when the approximate object is listed first.  As with other
///        overloads, this calls the previous form of the same comparison.
///
/// \param lhs            Vector of reference values with a tolerance and comparison style
/// \param relationship   The relation between the values, i.e. == or !=
/// \param rhs            The vector of input values
/// \param error_message  Error message to print
template <typename T>
CheckResult check(const Approx &lhs, const RelationalOperator relationship,
                  const std::vector<T> &rhs, const std::string &error_message,
                  TestPriority urgency = TestPriority::CRITICAL);

/// \brief Comparison of two vectors of any scalar types.  This feeds into the canonical form of
///        the check, making an Approx object from the second vector and then proceeding with the
///        comparison.  All scalar data types become double precision real.
template <typename T1, typename T2>
CheckResult check(const std::vector<T1> &lhs, const RelationalOperator relationship,
                  const std::vector<T2> &rhs, const std::string &error_message,
                  TestPriority urgency = TestPriority::CRITICAL);

/// \brief Compare a series of numbers to a file, or record them in a file of the given name.
///        The function returns SUCCESS if the file was written or if the data matched what was
///        recorded in the file.
///
/// Overloaded:
///   - Takes a tolerance, output precision, data format, and activity (commit data, or read data
///     as reference) in addition to the customary file name, content vector, label, and
///     error message
///   - Takes a format and activity and other optional factors (useful for cases in which the
///     primary use of snapshotting is to read and compare)
///
/// \param filename              Name of the file containing reference data
/// \param content               Data in RAM
/// \param label                 Label on the data to snapshot (serves as an identifier for human
///                              readers and a check that the intended data is being read)
/// \param error_message         Message to display if the snapshots fail to match (this functions
///                              in the same way as an error message from one of the check()
///                              functions, but displays "Snapshot MISMATCH:" on the left rather
///                              than "Check FAILED:".
/// \param comparison_tolerance  Tolerance by which to make comparisons between the curren data
///                              and reference data.  This can be a delta between observed and
///                              expected values or, in the case of text comparisons, a number of
///                              characters that may differ and still return a passing result.
/// \param activity              Whether to take or compare against the snapshot
/// \param output_precision      Absolute precision of the output numbers
/// \param data_format           Format of the data to write
/// \param expectation           Expected condition of the output file (i.e. "does not exist")
/// \param urgency               Indicate whether the test is possible and whether a failure should
///                              be considered blocking
/// \{
CheckResult snapshot(const std::string &filename, const std::vector<PolyNumeric> &content,
                     const std::string &label = std::string(""),
                     double comparison_tolerance = 1.0e-4,
                     const std::string &error_message = std::string(""),
                     SnapshotOperation activity = SnapshotOperation::COMPARE,
                     double output_precision = 1.0e-8,
                     NumberFormat data_format = NumberFormat::STANDARD_REAL,
                     PrintSituation expectation = PrintSituation::OVERWRITE,
                     TestPriority urgency = TestPriority::CRITICAL);

CheckResult snapshot(const std::string &filename, const std::vector<PolyNumeric> &content,
                     const std::string &label = std::string(""),
                     NumberFormat data_format = NumberFormat::STANDARD_REAL,
                     const std::string &error_message = std::string(""),
                     SnapshotOperation activity = SnapshotOperation::COMPARE,
                     double comparison_tolerance = 1.0e-4, double output_precision = 1.0e-8,
                     PrintSituation expectation = PrintSituation::OVERWRITE,
                     TestPriority urgency = TestPriority::CRITICAL);

CheckResult snapshot(const std::string &filename, const TextFile &content,
                     const std::string &label = std::string(""),
                     const std::string &error_message = std::string(""),
                     const SnapshotOperation activity = SnapshotOperation::COMPARE,
                     const PrintSituation expectation = PrintSituation::OVERWRITE,
                     const TestPriority urgency = TestPriority::CRITICAL,
                     int comparison_tolerance = 0);

CheckResult snapshot(const std::string &filename, const std::string &content,
                     const std::string &label = std::string(""),
                     const std::string &error_message = std::string(""),
                     const SnapshotOperation activity = SnapshotOperation::COMPARE,
                     const PrintSituation expectation = PrintSituation::OVERWRITE,
                     const TestPriority urgency = TestPriority::CRITICAL,
                     int comparison_tolerance = 0);
/// \}

/// \brief Print the summary for the global test result struct.
///
/// \param tv  The verboseness with which to report results
void printTestSummary(TestVerbosity tv);

/// \brief Count the overall number of failures in the global test result struct.
int countGlobalTestFailures();
  
/// \brief Get a string representing the type of numerical relationship, i.e. "==" or ">".
///
/// \param relationship  The operator in question
std::string getRelationalOperatorString(RelationalOperator ro);

} // namespace testing
} // namespace stormm

/// \brief ***Global*** instance of the checklist, analogous to the Ledger gbl_mem_balance_sheet
///        for tracking Hybrid array allocations (see src/Accelerator/hybrid.h)
extern stormm::testing::CheckList gbl_test_results;

#include "unit_test.tpp"

#endif
