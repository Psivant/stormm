#include <cmath>
#include "copyright.h"
#include "debug_nl_util.h"

namespace stormm {
namespace debug {

//-------------------------------------------------------------------------------------------------
bool cellAssignmentFailure(const int test_c, const int base_c, const int max_c,
                           const double test_fx, const double base_fx) {
  if (base_c != test_c) {

    // Check against a false positive
    if ((base_c - test_c == 1 && fabs(base_fx) < 0.001 && fabs(test_fx - 1.0) < 0.001) ||
        (test_c - base_c == 1 && fabs(base_fx - 1.0) < 0.001 && fabs(test_fx) < 0.001) ||
        (base_c == max_c - 1 && test_c == 0 &&
         fabs(base_fx - 1.0) < 0.001 && fabs(test_fx) < 0.001) ||
        (test_c == max_c - 1 && base_c == 0 &&
         fabs(base_fx) < 0.001 && fabs(test_fx - 1.0) < 0.001)) {
      return false;
    }
    else {
      return true;
    }
  }
  else {
    return false;
  }
  __builtin_unreachable();
}

} // namespace debug
} // namespace stormm
