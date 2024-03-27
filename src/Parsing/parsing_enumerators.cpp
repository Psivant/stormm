#include "copyright.h"
#include "parsing_enumerators.h"

namespace stormm {
namespace parse {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const NumberFormat input) {
  switch (input) {
  case NumberFormat::SCIENTIFIC:
    return std::string("SCIENTIFIC");
  case NumberFormat::STANDARD_REAL:
    return std::string("STANDARD_REAL");
  case NumberFormat::INTEGER:
    return std::string("INTEGER");
  case NumberFormat::LONG_LONG_INTEGER:
    return std::string("LONG_LONG_INTEGER");
  case NumberFormat::UNSIGNED_INTEGER:
    return std::string("UNSIGNED_INTEGER");
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
    return std::string("UNSIGNED_LONG_LONG_INTEGER");
  case NumberFormat::CHAR4:
    return std::string("CHAR4");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const NumberPrintStyle input) {
  switch (input) {
  case NumberPrintStyle::STANDARD:
    return std::string("STANDARD");
  case NumberPrintStyle::LEADING_ZEROS:
    return std::string("LEADING_ZEROS");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const WildCardKind input) {
  switch (input) {
  case WildCardKind::NONE:
    return std::string("NONE");
  case WildCardKind::FREE_CHARACTER:
    return std::string("FREE_CHARACTER");
  case WildCardKind::FREE_STRETCH:
    return std::string("FREE_STRETCH");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const JustifyText input) {
  switch (input) {
  case JustifyText::LEFT:
    return std::string("LEFT");
  case JustifyText::CENTER:
    return std::string("RIGHT");
  case JustifyText::RIGHT:
    return std::string("RIGHT");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TableHeadingLine input) {
  switch (input) {
  case TableHeadingLine::HEADER:
    return std::string("HEADER");
  case TableHeadingLine::HORIZONTAL_RULE:
    return std::string("HORIZONTAL_RULE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const BorderFormat input) {
  switch (input) {
  case BorderFormat::NONE:
    return std::string("NONE");
  case BorderFormat::LIGHT:
    return std::string("LIGHT");
  case BorderFormat::FULL:
    return std::string("FULL");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const WrapTextSearch input) {
  switch (input) {
  case WrapTextSearch::NO:
    return std::string("NO");
  case WrapTextSearch::YES:
    return std::string("YES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TextOrigin input) {
  switch (input) {
  case TextOrigin::DISK:
    return std::string("DISK");
  case TextOrigin::RAM:
    return std::string("RAM");
  }
  __builtin_unreachable();
}

} // namespace parse
} // namespace stormm
