#include "soft_core_potentials.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
double evaluateQuarticFirstDerivative(const double4 coeffs, const double r) {
  return (((((4.0 * coeffs.x * r) + (3.0 * coeffs.y)) * r) + (2.0 * coeffs.z)) * r) + coeffs.w;
}
  
//-------------------------------------------------------------------------------------------------
void adjustIntervalTarget(std::vector<double> *interval_targets, const int index,
                          const double f_rswitch, const double target_zero,
                          const double move_increment) {

  // Examine whether the intervals are increasing or decreasing as the inter-particle distance
  // goes from zero to the switching point.  Most softcore potentials will havve decreasing values,
  // as the point is to attain their highest potentials at r = 0.  Otherwise, the potential would
  // have a well where particles collide.  However, for generality, the code can handle such cases.
  const int npts = interval_targets->size();
  const size_t last_point = npts - 1;
  double* idata = interval_targets->data();
  bool rising = (f_rswitch > target_zero);
  const double min_move_incr = 0.1 + fabs(0.015625 * (f_rswitch - target_zero));
  double update;
  bool still_in_bounds;
  if (rising) {
    
    // If the point in question is already further from the potential at the switching point than
    // the intended target at r = 0, keep pushing it in that direction.  Otherwise, adjust it 25%
    // in the direction of the target at r = 0.  In either case, include a small increment that
    // ensures movement is never miniscule.
    if (idata[index] < target_zero) {
      still_in_bounds = false;
      update = -move_increment;
    }
    else {
      still_in_bounds = true;
      update = -move_increment;
    }
  }
  else {

    // The mirror of the other branch.  If the point in question is already greater than the
    // target potential as the inter-particle distance reaches zero, then keep pushing it in the
    // same direction.
    if (idata[index] > target_zero) {
      still_in_bounds = false;
      update = move_increment;
    }
    else {
      still_in_bounds = true;
      update = move_increment;
    }
  }
  if (still_in_bounds && target_zero - idata[index] > 0.01) {

    // Adjust spline segments upstream in proportion to the update on the current spline
    const double adj_proportion = update / (target_zero - idata[index]);
    for (int i = index; i >= 0; i--) {
      idata[i] += adj_proportion * (target_zero - idata[i]);
    }
  }
  else {

    // Add the same update to all spline segments upstream
    for (int i = index; i >= 0; i--) {
      idata[i] += update;
    }
  }
}

} // namespace energy
} // namespace stormm
