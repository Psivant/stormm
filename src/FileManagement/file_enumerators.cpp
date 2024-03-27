#include "copyright.h"
#include "file_enumerators.h"

namespace stormm {
namespace diskutil {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const PrintSituation input) {
  switch (input) {
  case PrintSituation::OPEN_NEW:
    return std::string("OPEN_NEW");
  case PrintSituation::APPEND:
    return std::string("APPEND");
  case PrintSituation::OVERWRITE:
    return std::string("OVERWRITE");
  case PrintSituation::UNKNOWN:
    return std::string("UNKNOWN");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const DataFormat input) {
  switch (input) {
  case DataFormat::ASCII:
    return std::string("ASCII");
  case DataFormat::BINARY:
    return std::string("BINARY");
  }
  __builtin_unreachable();
}

} // namespace diskutil
} // namespace stormm
