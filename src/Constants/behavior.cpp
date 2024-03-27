#include "copyright.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "behavior.h"

namespace stormm {
namespace constants {

using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
ExceptionResponse translateExceptionResponse(const std::string &policy) {
  if (strcmpCased(policy, std::string("die")) || strcmpCased(policy, std::string("abort"))) {
    return ExceptionResponse::DIE;
  }
  else if (strcmpCased(policy, std::string("warn")) || strcmpCased(policy, std::string("alert")) ||
           strcmpCased(policy, std::string("advise"))) {
    return ExceptionResponse::WARN;
  }
  else if (strcmpCased(policy, std::string("none")) ||
           strcmpCased(policy, std::string("silent"))) {
    return ExceptionResponse::SILENT;
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ExceptionResponse policy) {
  switch (policy) {
  case ExceptionResponse::DIE:
    return std::string("DIE");
  case ExceptionResponse::WARN:
    return std::string("WARN");
  case ExceptionResponse::SILENT:
    return std::string("SILENT");
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
PrecisionModel translatePrecisionModel(const std::string &choice) {
  if (strcmpCased(choice, std::string("single"), CaseSensitivity::NO) ||
      strcmpCased(choice, std::string("float32_t"), CaseSensitivity::NO)) {
    return PrecisionModel::SINGLE;
  }
  else if (strcmpCased(choice, std::string("double"), CaseSensitivity::NO) ||
           strcmpCased(choice, std::string("float64_t"), CaseSensitivity::NO)) {
    return PrecisionModel::DOUBLE;
  }
  else {
    rtErr("Invalid request for precision level " + choice + ".", "translatePrecisionModel");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const PrecisionModel plevel) {
  switch (plevel) {
  case PrecisionModel::SINGLE:
    return std::string("SINGLE");
  case PrecisionModel::DOUBLE:
    return std::string("DOUBLE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(CartesianDimension axis) {
  switch (axis) {
  case CartesianDimension::X:
    return std::string("X");
  case CartesianDimension::Y:
    return std::string("Y");
  case CartesianDimension::Z:
    return std::string("Z");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(UnitCellAxis axis) {
  switch (axis) {
  case UnitCellAxis::A:
    return std::string("a");
  case UnitCellAxis::B:
    return std::string("b");
  case UnitCellAxis::C:
    return std::string("c");
  }
  __builtin_unreachable();
}

} // namespace constants
} // namespace stormm
