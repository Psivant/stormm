// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace debug {

/// \brief Detect a failed assignment in the neighbor list along one unit cell axis.
///
/// \param test_c   Cell index of the assignment in the rebuilt neighbor list
/// \param base_c   Cell index of the assignment in the original simulation neighbor list
/// \param max_c    Number of neighbor list cells for the system of interest along the same axis
/// \param test_fx  Fractional coordinate of the particle in the rebuilt neighbor list
/// \param base_fx  Fractional coordinate of the particle in the original simulation neighbor list
bool cellAssignmentFailure(int test_c, int base_c, int max_c, double test_fx, double base_fx);

} // namespace debug
} // namespace stormm
