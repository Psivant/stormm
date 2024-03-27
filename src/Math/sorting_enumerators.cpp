#include "copyright.h"
#include "sorting_enumerators.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(SortDirection input) {
  switch (input) {
  case SortDirection::ASCENDING:
    return std::string("ASCENDING");
  case SortDirection::DESCENDING:
    return std::string("DESCENDING");
  case SortDirection::AUTOMATIC:
    return std::string("AUTOMATIC");
  }
  __builtin_unreachable();
};

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(UniqueValueHandling input) {
  switch (input) {
  case UniqueValueHandling::UNIQUE_VALUES_ONLY:
    return std::string("UNIQUE_VALUES_ONLY");    
  case UniqueValueHandling::CONFIRM_ALL_COPIES:
    return std::string("CONFIRM_ALL_COPIES");
  }
  __builtin_unreachable();
}

} // namespace stmath
} // namespace stormm
