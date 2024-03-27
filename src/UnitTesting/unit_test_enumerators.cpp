#include "copyright.h"
#include <string>
#include "unit_test_enumerators.h"

namespace stormm {
namespace testing {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ComparisonType input) {
  switch (input) {
  case ComparisonType::ABSOLUTE:
    return std::string("ABSOLUTE");
  case ComparisonType::RELATIVE:
    return std::string("RELATIVE");
  case ComparisonType::MEAN_UNSIGNED_ERROR:
    return std::string("MEAN_UNSIGNED_ERROR");
  case ComparisonType::RELATIVE_RMS_ERROR:
    return std::string("RELATIVE_RMS_ERROR");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const CheckResult input) {
  switch (input) {
  case CheckResult::SUCCESS:
    return std::string("SUCCESS");
  case CheckResult::SKIPPED:
    return std::string("SKIPPED");
  case CheckResult::IGNORED:
    return std::string("IGNORED");
  case CheckResult::FAILURE:
    return std::string("FAILURE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TestPriority input) {
  switch (input) {
  case TestPriority::CRITICAL:
    return std::string("CRITICAL");
  case TestPriority::NON_CRITICAL:
    return std::string("NON_CRITICAL");
  case TestPriority::ABORT:
    return std::string("ABORT");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RelationalOperator input) {
  switch (input) {
  case RelationalOperator::EQUAL:
  case RelationalOperator::EQ:
    return std::string("EQUAL");
  case RelationalOperator::NOT_EQUAL:
  case RelationalOperator::NE:
    return std::string("NOT_EQUAL");
  case RelationalOperator::GREATER_THAN:
  case RelationalOperator::GT:
    return std::string("GREATER_THAN");
  case RelationalOperator::LESS_THAN:
  case RelationalOperator::LT:
    return std::string("LESS_THAN");
  case RelationalOperator::GREATER_THAN_OR_EQUAL:
  case RelationalOperator::GE:
    return std::string("GREATER_THAN_OR_EQUAL");
  case RelationalOperator::LESS_THAN_OR_EQUAL:
  case RelationalOperator::LE:
    return std::string("LESS_THAN_OR_EQUAL");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TestVerbosity input) {
  switch (input) {
  case TestVerbosity::FULL:
    return std::string("FULL");
  case TestVerbosity::COMPACT:
    return std::string("COMPACT");
  case TestVerbosity::FAILURE_ONLY:
    return std::string("FAILURE_ONLY");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const SnapshotOperation input) {
  switch (input) {
  case SnapshotOperation::COMPARE:
    return std::string("COMPARE");
  case SnapshotOperation::SNAPSHOT:
    return std::string("SNAPSHOT");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TmpdirStatus input) {
  switch (input) {
  case TmpdirStatus::NOT_REQUIRED:
    return std::string("NOT_REQUIRED");
  case TmpdirStatus::REQUIRED:
    return std::string("REQUIRED");
  }
  __builtin_unreachable();
}

} // namespace testing
} // namespace stormm
