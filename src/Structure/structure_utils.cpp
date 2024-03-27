#include <string>
#include "copyright.h"
#include "Reporting/error_format.h"
#include "structure_utils.h"

namespace stormm {
namespace structure {
  
//-------------------------------------------------------------------------------------------------
void coordinateBoundsCheck(const int lower_limit, const int upper_limit, const int natom,
                           const char* caller) {
  if (lower_limit < 0 || upper_limit < 0 || lower_limit >= natom || upper_limit > natom) {
    rtErr("A coordinate set with " + std::to_string(natom) + " particles cannot rotate "
          "indices " + std::to_string(lower_limit) + " to " + std::to_string(upper_limit) + ".",
          caller);
  }
}

} // namespace structure
} // namespace stormm
