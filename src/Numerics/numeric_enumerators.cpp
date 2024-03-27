#include "numeric_enumerators.h"

namespace stormm {
namespace numerics {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const AccumulationMethod input) {
  switch (input) {
  case AccumulationMethod::SPLIT:
    return std::string("SPLIT");
  case AccumulationMethod::WHOLE:
    return std::string("WHOLE");
  case AccumulationMethod::AUTOMATIC:
    return std::string("AUTOMATIC");
  }
  __builtin_unreachable();
}
  
} // namespace numerics
} // namespace stormm
