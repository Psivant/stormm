// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename T> T roundUp(T jagged, T increment) {
  return ((jagged + increment - static_cast<T>(1)) / increment) * increment;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc angleVerification(const Tcalc costheta, const Tcalc* crabbc, const Tcalc* crbccd,
                        const Tcalc* bc, const Tcalc* scr) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  if (tcalc_ct == double_type_index) {
    if (std::abs(costheta) >= near_to_one_lf) {

      // The Tcalcing-point representation of costheta is numerically ill-conditioned.  Compute
      // the distance from atom I to the plane of atoms J, K, and L to get the angle by the
      // arcsin of an extremely acute angle.
      const Tcalc mg_crabbc = 1.0 / sqrt(crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] +
                                         crabbc[2]*crabbc[2]);
      const Tcalc mg_crbccd = 1.0 / sqrt(crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] +
                                         crbccd[2]*crbccd[2]);
      const Tcalc nx_abbc = crabbc[0] * mg_crabbc;
      const Tcalc ny_abbc = crabbc[1] * mg_crabbc;
      const Tcalc nz_abbc = crabbc[2] * mg_crabbc;
      const Tcalc nx_bccd = crbccd[0] * mg_crbccd;
      const Tcalc ny_bccd = crbccd[1] * mg_crbccd;
      const Tcalc nz_bccd = crbccd[2] * mg_crbccd;
      Tcalc rdx = nx_bccd - nx_abbc;
      Tcalc rdy = ny_bccd - ny_abbc;
      Tcalc rdz = nz_bccd - nz_abbc;
      Tcalc rs = sqrt((rdx * rdx) + (rdy * rdy) + (rdz * rdz));
      if (fabs(rs) > 1.0) {
        rdx = nx_bccd + nx_abbc;
        rdy = ny_bccd + ny_abbc;
        rdz = nz_bccd + nz_abbc;
        rs = pi - sqrt((rdx * rdx) + (rdy * rdy) + (rdz * rdz));
      }
      return (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0f) ? rs : -rs;
    }
    else {
      return (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0) ?  acos(costheta) :
                                                                  -acos(costheta);
    }
  }
  else {
    if (std::abs(costheta) >= near_to_one_f) {

      // The Tcalcing-point representation of costheta is numerically ill-conditioned.  Compute
      // the distance from atom I to the plane of atoms J, K, and L to get the angle by the
      // arcsin of an extremely acute angle.
      const Tcalc mg_crabbc = 1.0f / sqrtf(crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] +
                                           crabbc[2]*crabbc[2]);
      const Tcalc mg_crbccd = 1.0f / sqrtf(crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] +
                                           crbccd[2]*crbccd[2]);
      const Tcalc nx_abbc = crabbc[0] * mg_crabbc;
      const Tcalc ny_abbc = crabbc[1] * mg_crabbc;
      const Tcalc nz_abbc = crabbc[2] * mg_crabbc;
      const Tcalc nx_bccd = crbccd[0] * mg_crbccd;
      const Tcalc ny_bccd = crbccd[1] * mg_crbccd;
      const Tcalc nz_bccd = crbccd[2] * mg_crbccd;
      Tcalc rdx = nx_bccd - nx_abbc;
      Tcalc rdy = ny_bccd - ny_abbc;
      Tcalc rdz = nz_bccd - nz_abbc;
      Tcalc rs = sqrtf((rdx * rdx) + (rdy * rdy) + (rdz * rdz));
      if (std::abs(rs) > 1.0f) {
        rdx = nx_bccd + nx_abbc;
        rdy = ny_bccd + ny_abbc;
        rdz = nz_bccd + nz_abbc;
        rs = pi_f - sqrtf((rdx * rdx) + (rdy * rdy) + (rdz * rdz));
      }
      return (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0f) ? rs : -rs;
    }
    else {
      return (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0f) ?  acosf(costheta) :
                                                                   -acosf(costheta);
    }
  }
  __builtin_unreachable();
}

} // namespace stmath
} // namespace stormm
