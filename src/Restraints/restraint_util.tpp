// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace restraints {

//-------------------------------------------------------------------------------------------------
template <typename T>
Vec3<T> restraintDelta(const Vec2<T> init_k, const Vec2<T> final_k, const Vec4<T> init_r,
                       const Vec4<T> final_r, const Vec2<T> mixwt, const T dr) {
  const T r1 = (mixwt.x * init_r.x) + (mixwt.y * final_r.x);
  const T r2 = (mixwt.x * init_r.y) + (mixwt.y * final_r.y);
  const T r3 = (mixwt.x * init_r.z) + (mixwt.y * final_r.z);
  const T r4 = (mixwt.x * init_r.w) + (mixwt.y * final_r.w);
  const T k2 = (mixwt.x * init_k.x) + (mixwt.y * final_k.x);
  const T k3 = (mixwt.x * init_k.y) + (mixwt.y * final_k.y);
  T dl, du, keq;
  if (dr < r1) {
    dl = r1 - r2;
    du = k2 * ((dl * dl) + (2.0 * dl * (dr - r1)));
    keq = k2;
  }
  else if (dr < r2) {
    dl = dr - r2;
    du = k2 * dl * dl;
    keq = k2;
  }
  else if (dr < r3) {
    dl = 0.0;
    du = 0.0;
    keq = 0.0;
  }
  else if (dr < r4) {
    dl = dr - r3;
    du = k3 * dl * dl;
    keq = k3;
  }
  else {
    dl = r4 - r3;
    du = k3 * ((dl * dl) + (2.0 * dl * (dr - r4)));
    keq = k3;
  }
  return Vec3<T>(keq, dl, du);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
Vec2<T> computeRestraintMixture(const int step_number, const int init_step, const int final_step) {
  if (step_number < init_step) {

    // If the restraint has not yet engaged, neither its initial or final values have any weight
    return Vec2<T>(0.0, 0.0);
  }
  else if (init_step == final_step) {

    // The step count is far enough along that the restraint has been engaged, and it is constant.
    // Only the initial value matters.
    return Vec2<T>(1.0, 0.0);
  }
  else if (step_number < final_step) {
    const T wslide = static_cast<T>(step_number - init_step) /
                     static_cast<T>(final_step - init_step);

    // The difference between the initial and final steps is nonzero.  The mixture is a linear
    // combination of the two end points.
    return Vec2<T>(1.0 - wslide, wslide);
  }
  else {

    // The step number has advanced beyond the point at which the restraint is mature.
    return Vec2<T>(0.0, 1.0);
  }
  __builtin_unreachable();
}

} // namespace restraints
} // namespace stormm
