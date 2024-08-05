#include <algorithm>
#include <cstdio>
#include <cstring>
#include <typeinfo>
#include <typeindex>
#include <unistd.h>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "FileManagement/directory_util.h"
#include "FileManagement/file_listing.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/tabulation.h"
#include "unit_test.h"

//-------------------------------------------------------------------------------------------------
stormm::testing::CheckList gbl_test_results;

namespace stormm {
namespace testing {

using constants::ExceptionResponse;
using diskutil::DrivePathType;
using diskutil::findStormmPath;
using diskutil::getDrivePathType;
using diskutil::makePathAbsolute;
using diskutil::stormmMkdir;
using diskutil::stormmBatchRmdir;
using diskutil::stormmRmdir;
using diskutil::openOutputFile;
using diskutil::osSeparator;
using diskutil::removeFile;
using diskutil::separatePath;
using errors::terminalFormat;
using parse::BorderFormat;
using parse::char4ToString;
using parse::findStringInVector;
using parse::intToString;
using parse::NumberFormat;
using data_types::operator==;
using parse::polyNumericVector;
using parse::printTable;
using parse::realDecimalPlaces;
using parse::realToString;
using parse::TextFileReader;
using parse::TextOrigin;

//-------------------------------------------------------------------------------------------------
void section(const std::string &section_name) {
  gbl_test_results.changeSection(section_name);
}

//-------------------------------------------------------------------------------------------------
void section(const int section_index) {
  gbl_test_results.changeSection(section_index);
}

//-------------------------------------------------------------------------------------------------
int unitTestSectionIndex(const std::string &name) {
  return gbl_test_results.getSectionIndex(name);
}

//-------------------------------------------------------------------------------------------------
std::string unitTestSectionName(const int index) {
  return gbl_test_results.getSectionName(index);
}

//-------------------------------------------------------------------------------------------------
CheckResult check(const bool statement, const std::string &error_message,
                  const TestPriority urgency) {

  // Skip the test if some prior condition negates running it or would invalidate the result.
  if (urgency == TestPriority::ABORT) {
    gbl_test_results.logResult(CheckResult::SKIPPED);
    return CheckResult::SKIPPED;
  }

  // Evaluate the result
  if (statement) {

    // Update the global test registry
    gbl_test_results.logResult(CheckResult::SUCCESS);
    return CheckResult::SUCCESS;
  }
  else {
    printf("Check FAILED: ");
    std::string criticality("");
    switch (urgency) {
    case TestPriority::CRITICAL:
      break;
    case TestPriority::ABORT:
    case TestPriority::NON_CRITICAL:
      criticality = "[Non-critical] ";
      break;
    }
    if (error_message.size() > 0) {
      const std::string parsed_msg = terminalFormat(criticality + error_message, "", "", 14, 0,
                                                    14);
      printf("%s\n", parsed_msg.c_str());
    }
    else {
      printf("%sNo description provided.\n", criticality.c_str());
    }
    if (urgency == TestPriority::NON_CRITICAL) {
      gbl_test_results.logResult(CheckResult::IGNORED);
      return CheckResult::IGNORED;
    }
    else if (urgency == TestPriority::CRITICAL) {
      gbl_test_results.logResult(CheckResult::FAILURE);
      return CheckResult::FAILURE;
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
CheckResult check(const double lhs, const RelationalOperator relationship, const Approx &rhs,
                  const std::string &error_message, const TestPriority urgency) {

  // Skip the test if some prior condition negates running it or would invalidate the result.
  if (urgency == TestPriority::ABORT) {
    gbl_test_results.logResult(CheckResult::SKIPPED);
    return CheckResult::SKIPPED;
  }

  // Check that the comparison is valid
  switch (relationship) {
  case RelationalOperator::EQUAL:
  case RelationalOperator::EQ:
  case RelationalOperator::NOT_EQUAL:
  case RelationalOperator::NE:

    // All comparisons are valid
    break;
  case RelationalOperator::GREATER_THAN:
  case RelationalOperator::GT:
  case RelationalOperator::LESS_THAN:
  case RelationalOperator::LT:
  case RelationalOperator::GREATER_THAN_OR_EQUAL:
  case RelationalOperator::GE:
  case RelationalOperator::LESS_THAN_OR_EQUAL:
  case RelationalOperator::LE:

    // Only absolute comparisons are valid (MEAN_UNSIGNED_ERROR works as absolute for scalars)
    switch (rhs.getStyle()) {
    case ComparisonType::ABSOLUTE:
    case ComparisonType::MEAN_UNSIGNED_ERROR:
      break;
    case ComparisonType::RELATIVE:
    case ComparisonType::RELATIVE_RMS_ERROR:
      rtErr("Relative comparisons are not permitted in combination with " +
            getRelationalOperatorString(relationship) + " relationships.", "check");
    }
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
  // the decimal precision of the target value.
  int idec = realDecimalPlaces(rhs.getTol());
  int lhs_dec = realDecimalPlaces(lhs);
  int rhs_dec = realDecimalPlaces(rhs.getValue());
  bool both_are_integers = (lhs_dec == 0 && rhs_dec == 0);
  idec = (idec > 9 && both_are_integers) ? 0 : idec;

  // Report the error.  Add an alert if the test is considered non-critical.
  std::string error_edit("  ");
  switch (rhs.getStyle()) {
  case ComparisonType::ABSOLUTE:
  case ComparisonType::MEAN_UNSIGNED_ERROR:
    switch (relationship) {
    case RelationalOperator::EQUAL:
    case RelationalOperator::EQ:
      if (idec > 0) {
        error_edit += "Absolute comparison of " + realToString(lhs, idec) + " to " +
                      realToString(rhs.getValue(), idec) + " +/- " +
                      realToString(rhs.getTol(), idec) + " fails by " +
                      realToString(fabs(lhs - rhs.getValue()), idec) + ".";
      }
      else {
        error_edit += "Absolute comparison of " + realToString(lhs, idec) + " to " +
                      realToString(rhs.getValue(), idec) + " fails.";
      }
      break;
    case RelationalOperator::NOT_EQUAL:
    case RelationalOperator::NE:
      if (idec > 0) {
        error_edit += "Absolute comparison of " + realToString(lhs, idec) + " to " +
                      realToString(rhs.getValue(), idec) + " was expected to fail but falls "
                      "within the tolerance of " + realToString(rhs.getTol(), idec) + ".";
      }
      else {
        error_edit += "Absolute comparison of " + realToString(lhs, idec) + " to " +
                      realToString(rhs.getValue(), idec) + " was expected to fail.";
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
      error_edit += "Absolute comparison " + realToString(lhs, 11, 4, NumberFormat::SCIENTIFIC) +
                    " " + getRelationalOperatorString(relationship) + " " +
                    realToString(rhs.getValue(), 11, 4, NumberFormat::SCIENTIFIC) +
                    " failed.  The inequality could not be supported beyond the marginal value "
                    "associated with the reference.";
      break;
    }
    break;
  case ComparisonType::RELATIVE:
  case ComparisonType::RELATIVE_RMS_ERROR:
    switch (relationship) {
    case RelationalOperator::EQUAL:
    case RelationalOperator::EQ:
      if (idec > 0) {
        error_edit += "Relative comparison of " + realToString(lhs, idec) + " to " +
                      realToString(rhs.getValue(), idec) + " +/- " +
                      realToString(rhs.getTol(), idec) + " fails by " +
                      realToString(100.0 * fabs((lhs - rhs.getValue()) / (rhs.getTol())), idec) +
                      "%.";
      }
      else {
        error_edit += "Relative comparison of " + realToString(lhs, idec) + " to " +
                      realToString(rhs.getValue(), idec) + " +/- " +
                      realToString(rhs.getTol(), idec) + " fails.";
      }
      break;
    case RelationalOperator::NOT_EQUAL:
    case RelationalOperator::NE:
      if (idec > 0) {
        error_edit += "Relative comparison of " + realToString(lhs, idec) + " to " +
                      realToString(rhs.getValue(), idec) + " was expected to fail but falls "
                      "within the tolerance of " + realToString(rhs.getTol(), idec) + ".";
      }
      else {
        error_edit += "Relative comparison of " + realToString(lhs, idec) + " to " +
                      realToString(rhs.getValue(), idec) + " was expected to fail.";
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

      // These relationships are invalid for non-absolute comparisons, and trapped earlier.
      break;
    }
    break;
  }
  return check(false, error_message + error_edit, urgency);
}

//-------------------------------------------------------------------------------------------------
CheckResult check(const Approx &lhs, const RelationalOperator relationship, const double rhs,
                  const std::string &error_message, const TestPriority urgency) {
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
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
CheckResult check(const double lhs, const RelationalOperator relationship, const double rhs,
                  const std::string &error_message, const TestPriority urgency) {
  if (fabs(rhs) > 1.0) {
    return check(lhs, relationship, Approx(rhs).margin(fabs(rhs) * 1.0e-6), error_message,
                 urgency);
  }
  else if (fabs(rhs >= 0.01)) {
    return check(lhs, relationship, Approx(rhs).margin(constants::small), error_message, urgency);
  }
  else if (fabs(rhs >= 0.0001)) {
    return check(lhs, relationship, Approx(rhs).margin(constants::tiny), error_message, urgency);
  }
  else {
    return check(lhs, relationship, Approx(rhs).margin(constants::verytiny), error_message,
                 urgency);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
CheckResult check(const std::string &lhs, const RelationalOperator relationship,
                  const std::string &rhs, const std::string &error_message,
                  const TestPriority urgency) {

  // Skip the test if some prior condition negates running it or would invalidate the result.
  if (urgency == TestPriority::ABORT) {
    gbl_test_results.logResult(CheckResult::SKIPPED);
    return CheckResult::SKIPPED;
  }

  // String comparisons are only valid for EQUAL and NOT_EQUAL relationships
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
    rtErr("String comparisons are invalid for " + getRelationalOperatorString(relationship) +
          "relationships.", "check");
  }

  // Enumerate the possible relations between the two values
  std::string error_edit("  ");
  switch (relationship) {
  case RelationalOperator::EQUAL:
  case RelationalOperator::EQ:
    if (lhs == rhs) {
      gbl_test_results.logResult(CheckResult::SUCCESS);
      return CheckResult::SUCCESS;
    }
    error_edit += "Comparison of \"" + lhs + "\" to \"" + rhs + "\" fails.";
    break;
  case RelationalOperator::NOT_EQUAL:
  case RelationalOperator::NE:
    if (lhs != rhs) {
      gbl_test_results.logResult(CheckResult::SUCCESS);
      return CheckResult::SUCCESS;
    }
    error_edit += "Comparison of " + lhs + " to " + rhs + " was expected to fail.";
    break;
  case RelationalOperator::GREATER_THAN:
  case RelationalOperator::GT:
  case RelationalOperator::LESS_THAN:
  case RelationalOperator::LT:
  case RelationalOperator::GREATER_THAN_OR_EQUAL:
  case RelationalOperator::GE:
  case RelationalOperator::LESS_THAN_OR_EQUAL:
  case RelationalOperator::LE:

    // These cases can no longer be reached, after traps earlier in this function
    break;
  }
  return check(false, error_message + error_edit, urgency);
}

//-------------------------------------------------------------------------------------------------
CheckResult check(const char4 &lhs, const RelationalOperator relationship, const char4 &rhs,
                  const std::string &error_message, const TestPriority urgency) {
  return check(char4ToString(lhs), relationship, char4ToString(rhs), error_message, urgency);
}

//-------------------------------------------------------------------------------------------------
CheckResult check(const std::vector<char4> &lhs, const RelationalOperator relationship,
                  const std::vector<char4> &rhs, const std::string &error_message,
                  const TestPriority urgency) {

  // Skip the test if some prior condition negates running it or would invalidate the result.
  if (urgency == TestPriority::ABORT) {
    gbl_test_results.logResult(CheckResult::SKIPPED);
    return CheckResult::SKIPPED;
  }

  // Check each entry of the two vectors
  std::string error_edit("  ");
  const size_t nlhs = lhs.size();
  if (nlhs != rhs.size()) {
    error_edit += "Comparison of char4 vectors of different lengths (" + std::to_string(nlhs) +
                  " and " + std::to_string(rhs.size()) + ") is nonsensical.";
    return check(false, error_message + error_edit, urgency);
  }
  bool same = true;
  size_t nfailures = 0;
  for (size_t i = 0; i < nlhs; i++) {
    const bool match = (lhs[i] == rhs[i]);
    same = (same && match);
    nfailures += match;
  }
  
  // String comparisons are only valid for EQUAL and NOT_EQUAL relationships
  switch (relationship) {
  case RelationalOperator::EQUAL:
  case RelationalOperator::EQ:
    error_edit += "Comparison of char4 vectors of length " + std::to_string(nlhs) + " fails at " +
                  std::to_string(nfailures) + " points.";
    return check(same, error_message + error_edit, urgency);
  case RelationalOperator::NOT_EQUAL:
  case RelationalOperator::NE:
    error_edit += "Comparison of char4 vectors of length " + std::to_string(nlhs) +
                  " was expected to fail.";
    return check(same == false, error_message + error_edit, urgency);
  case RelationalOperator::GREATER_THAN:
  case RelationalOperator::GT:
  case RelationalOperator::LESS_THAN:
  case RelationalOperator::LT:
  case RelationalOperator::GREATER_THAN_OR_EQUAL:
  case RelationalOperator::GE:
  case RelationalOperator::LESS_THAN_OR_EQUAL:
  case RelationalOperator::LE:
    rtErr("Vector comparisons of char4 are invalid for " +
          getRelationalOperatorString(relationship) + "relationships.", "check");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
CheckResult snapshot(const std::string &filename, const std::vector<PolyNumeric> &content,
                     const std::string &label, const double comparison_tolerance,
                     const std::string &error_message, const SnapshotOperation activity,
                     const double output_precision, const NumberFormat data_format,
                     const PrintSituation expectation, const TestPriority urgency) {

  // Abort if some prior condition has made this test impossible to run
  if (urgency == TestPriority::ABORT) {

    // Check first that the point is not to write the snapshot
    switch (activity) {
    case SnapshotOperation::COMPARE:
      gbl_test_results.logResult(CheckResult::SKIPPED);
      return CheckResult::SKIPPED;
    case SnapshotOperation::SNAPSHOT:
      writeSnapshot(filename, content, label, output_precision, data_format, expectation);
      gbl_test_results.logResult(CheckResult::SUCCESS);
      return CheckResult::SUCCESS;
    }
  }

  // To write the file or to read the file... that is the question.
  switch (activity) {
  case SnapshotOperation::COMPARE:
    {
      const std::vector<PolyNumeric> ref_content = readSnapshot(filename, label);
      bool mismatch = false;
      if (ref_content.size() != content.size()) {
        mismatch = true;
      }
      else {
        switch (data_format) {
        case NumberFormat::SCIENTIFIC:
        case NumberFormat::STANDARD_REAL:
          mismatch = (doubleFromPolyNumeric(content) !=
                      Approx(doubleFromPolyNumeric(ref_content)).margin(comparison_tolerance));
          break;
        case NumberFormat::INTEGER:
          mismatch = (intFromPolyNumeric(content) !=
                      Approx(intFromPolyNumeric(ref_content)).margin(comparison_tolerance));
          break;
        case NumberFormat::LONG_LONG_INTEGER:
          mismatch = (llintFromPolyNumeric(content) !=
                      Approx(llintFromPolyNumeric(ref_content)).margin(comparison_tolerance));
          break;
        case NumberFormat::UNSIGNED_INTEGER:
          mismatch = (uintFromPolyNumeric(content) !=
                      Approx(uintFromPolyNumeric(ref_content)).margin(comparison_tolerance));
          break;
        case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
          mismatch = (ullintFromPolyNumeric(content) !=
                      Approx(ullintFromPolyNumeric(ref_content)).margin(comparison_tolerance));
          break;
        case NumberFormat::CHAR4:
          break;
        }
      }
      if (mismatch) {
        printf("Snapshot MISMATCH: ");
        std::string error_edit("  ");
        error_edit += vectorAlignmentReport(ref_content, content, data_format,
                                            comparison_tolerance);
        std::string criticality;
        switch (urgency) {
        case TestPriority::CRITICAL:
          break;
        case TestPriority::ABORT:
        case TestPriority::NON_CRITICAL:
          criticality = "[Non-critical] ";
          break;
        }
        const std::string parsed_msg = terminalFormat(criticality + error_message + error_edit,
                                                      "", "", 19, 0, 19);
        printf("%s\n", parsed_msg.c_str());
        if (urgency == TestPriority::NON_CRITICAL) {
          gbl_test_results.logResult(CheckResult::IGNORED);
          return CheckResult::IGNORED;
        }
        else if (urgency == TestPriority::CRITICAL) {
          gbl_test_results.logResult(CheckResult::FAILURE);
          return CheckResult::FAILURE;
        }
      }
      else {
        gbl_test_results.logResult(CheckResult::SUCCESS);
        return CheckResult::SUCCESS;
      }
    }
    break;
  case SnapshotOperation::SNAPSHOT:
    writeSnapshot(filename, content, label, output_precision, data_format, expectation);
    gbl_test_results.logResult(CheckResult::SUCCESS);
    return CheckResult::SUCCESS;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
CheckResult snapshot(const std::string &filename, const std::vector<PolyNumeric> &content,
                     const std::string &label, const NumberFormat data_format,
                     const std::string &error_message, const SnapshotOperation activity,
                     const double comparison_tolerance, const double output_precision,
                     const PrintSituation expectation, const TestPriority urgency) {
  return snapshot(filename, content, label, comparison_tolerance, error_message, activity,
                  output_precision, data_format, expectation, urgency);
}

//-------------------------------------------------------------------------------------------------
CheckResult snapshot(const std::string &filename, const TextFile &content,
                     const std::string &label, const std::string &error_message,
                     const SnapshotOperation activity, const PrintSituation expectation,
                     const TestPriority urgency, const int comparison_tolerance) {
  
  // Abort if some prior condition has made this test impossible to run
  if (urgency == TestPriority::ABORT) {

    // Check first that the point is not to write the snapshot
    switch (activity) {
    case SnapshotOperation::COMPARE:
      gbl_test_results.logResult(CheckResult::SKIPPED);
      return CheckResult::SKIPPED;
    case SnapshotOperation::SNAPSHOT:
      writeTextSnapshot(filename, content, label, expectation);
      gbl_test_results.logResult(CheckResult::SUCCESS);
      return CheckResult::SUCCESS;
    }
  }

  // For a non-abortive test, perform the requested snapshotting activity.
  switch (activity) {
  case SnapshotOperation::COMPARE:
    {
      const TextFile ref_content = readTextSnapshot(filename, label);
      const TextFileReader rtrial = content.data();
      const TextFileReader rref = ref_content.data();
      bool mismatch = (rtrial.line_count != rref.line_count);
      const bool line_count_difference = mismatch;
      int i = 0;
      std::vector<int> line_length_differences;
      while (i < rtrial.line_count) {
        if (rtrial.line_limits[i + 1] - rtrial.line_limits[i] !=
            rref.line_limits[i + 1] - rref.line_limits[i]) {
          line_length_differences.push_back(i + 1);
        }
        i++;
      }
      mismatch = (mismatch || line_length_differences.size() > 0);
      size_t altered_chars = 0;
      std::vector<int> line_content_differences;
      i = 0;
      while (i < rtrial.line_count) {
        bool line_noted = false;
        for (size_t ic = rtrial.line_limits[i]; ic < rtrial.line_limits[i + 1]; ic++) {
          if (rtrial.text[ic] != rref.text[ic]) {
            altered_chars++;
            if (line_noted == false) {
              line_content_differences.push_back(i + 1);
              line_noted = true;
            }
          }
        }
        i++;
      }
      mismatch = (mismatch ||
                  (line_content_differences.size() > 0 && altered_chars > comparison_tolerance));
      if (mismatch) {
        printf("Snapshot MISMATCH: ");
        std::string error_edit("  ");
        if (line_count_difference) {
          error_edit += "The text samples contain different line counts (reference " +
                        std::to_string(rref.line_count) + ", comparison " +
                        std::to_string(rtrial.line_count) + ").";
        }
        else if (line_length_differences.size() > 0) {
          const int ndiff = line_length_differences.size();
          error_edit += "In all, " + std::to_string(ndiff) + " out of " +
                        std::to_string(rtrial.line_count) + " lines differ in length between "
                        "the two samples (entry \"" + label + "\"): [ ";
          int nrep = 0;
          bool elipsis_printed = false;
          for (int j = 0; j < ndiff; j++) {
            if (nrep < 7 || j == ndiff - 1) {
              error_edit += std::to_string(line_length_differences[j]);
              if (j != ndiff - 1) {
                error_edit += " ";
              }
            }
            else if (elipsis_printed == false) {
              elipsis_printed = true;
              error_edit += "... ";
            }
          }
          error_edit += " ].";
        }
        else if (line_content_differences.size() > 0) {
          const int ndiff = line_content_differences.size();
          error_edit += "In all, " + std::to_string(line_length_differences.size()) + " out of " +
                        std::to_string(rtrial.line_count) + " lines differ between the two "
                        "samples (entry \"" + label + "\"), with " +
                        std::to_string(altered_chars) + " total differences between the "
                        "characters: [ ";
          int nrep = 0;
          bool elipsis_printed = false;
          for (int j = 0; j < ndiff; j++) {
            if (nrep < 7 || j == ndiff - 1) {
              error_edit += std::to_string(line_content_differences[j]);
              if (j != ndiff - 1) {
                error_edit += " ";
              }
            }
            else if (elipsis_printed == false) {
              elipsis_printed = true;
              error_edit += "... ";
            }
          }
          error_edit += " ].";
        }
        std::string criticality("");
        switch (urgency) {
        case TestPriority::CRITICAL:
          break;
        case TestPriority::ABORT:
        case TestPriority::NON_CRITICAL:
          criticality = "[Non-critical] ";
          break;
        }
        const std::string parsed_msg = terminalFormat(criticality + error_message + error_edit,
                                                      "", "", 19, 0, 19);
        printf("%s\n", parsed_msg.c_str());
        if (urgency == TestPriority::NON_CRITICAL) {
          gbl_test_results.logResult(CheckResult::IGNORED);
          return CheckResult::IGNORED;
        }
        else if (urgency == TestPriority::CRITICAL) {
          gbl_test_results.logResult(CheckResult::FAILURE);
          return CheckResult::FAILURE;
        }
      }
      else {
        gbl_test_results.logResult(CheckResult::SUCCESS);
        return CheckResult::SUCCESS;
      }
    }
    break;
  case SnapshotOperation::SNAPSHOT:
    writeTextSnapshot(filename, content, label, expectation);
    gbl_test_results.logResult(CheckResult::SUCCESS);
    return CheckResult::SUCCESS;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
CheckResult snapshot(const std::string &filename, const std::string &content,
                     const std::string &label, const std::string &error_message,
                     const SnapshotOperation activity, const PrintSituation expectation,
                     const TestPriority urgency, const int comparison_tolerance) {
  return snapshot(filename, TextFile(content, TextOrigin::RAM), label, error_message, activity,
                  expectation, urgency, comparison_tolerance);
}

//-------------------------------------------------------------------------------------------------
void printTestSummary(const TestVerbosity tv) {
  gbl_test_results.printSummary(tv);
}

//-------------------------------------------------------------------------------------------------
int countGlobalTestFailures() {
  return gbl_test_results.getOverallFailureCount() + gbl_test_results.getOverallSkipCount();
}

//-------------------------------------------------------------------------------------------------
std::string getRelationalOperatorString(const RelationalOperator ro) {
  switch (ro) {
  case RelationalOperator::EQUAL:
  case RelationalOperator::EQ:
    return "==";
  case RelationalOperator::NOT_EQUAL:
  case RelationalOperator::NE:
    return "!=";
  case RelationalOperator::GREATER_THAN:
  case RelationalOperator::GT:
    return ">";
  case RelationalOperator::LESS_THAN:
  case RelationalOperator::LT:
    return "<";
  case RelationalOperator::GREATER_THAN_OR_EQUAL:
  case RelationalOperator::GE:
    return ">=";
  case RelationalOperator::LESS_THAN_OR_EQUAL:
  case RelationalOperator::LE:
    return "<=";
  }
  __builtin_unreachable();
}

} // namespace testing
} // namespace stormm

