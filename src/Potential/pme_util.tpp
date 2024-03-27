// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc elecPMEDirectSpace(const double ew_coeff, const Tcalc kcoul, const Tcalc r,
                         const Tcalc r2, const int order) {
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const Tcalc ew_t = ew_coeff;
  const Tcalc value_one = 1.0;
  const Tcalc invr = value_one / r;
  const Tcalc u_quant = (tcalc_is_double) ? erfc(ew_t * r) * invr : erfcf(ew_t * r) * invr;
  if (order == 0) {
    return kcoul * u_quant;
  }
  const Tcalc bfac = 2.0 * ew_coeff / sqrt(symbols::pi);
  const Tcalc exp_quant = (tcalc_is_double) ? bfac * exp(-ew_t * ew_t * r2) :
                                              bfac * expf(-ew_t * ew_t * r2);
  Tcalc du, d2u;
  if (order == 1) {
    return -kcoul * (exp_quant + u_quant) * invr;
  }
  else if (order == 2) {
    const Tcalc value_two = 2.0;
    const Tcalc value_three = 3.0;
    const Tcalc ew_t2 = ew_coeff * ew_coeff;
    const Tcalc invr2 = value_one / r2;

    // (d/dr) [ -K (B exp(-a^2 r^2) + (erfc(a r) / r)) / r ] =
    // -K (d/dr) [ (B exp(-a^2 r^2) + (erfc(a r) / r)) / r ] =
    //
    // -K ( (d/dr) [ B exp(-a^2 r^2) / r ] +
    //      (d/dr) [ (erfc(a r) / r) / r ]    ) =
    //
    // -K ( [ -2 a^2 B r exp(-a^2 r^2) / r - B exp(-a^2 r^2) / r^2 ] +
    //      (d/dr) [ (erfc(a r) / r) / r ]    ) =
    //
    // -K ( [ -2 a^2 B exp(-a^2 r^2) - B exp(-a^2 r^2) / r^2 ] +
    //      (d/dr) [ (erfc(a r) / r) / r ]    ) =
    //
    // -K ( [ -(2 a^2 + 1 / r^2) B exp(-a^2 r^2) ] +
    //      (d/dr) [ (erfc(a r) / r) / r ]    ) =
    //
    // -K ( [ -(2 a^2 + 1 / r^2) B exp(-a^2 r^2) ] +
    //      [ -((B exp(-a^2 r^2) + (erfc(a r) / r)) / r^2) - ((erfc(a r) / r) / r^2) ] =
    //
    // -K ( [ -(2 a^2 + 1 / r^2) B exp(-a^2 r^2) ] +
    //      [ -((B exp(-a^2 r^2) / r^2) + (erfc(a r) / r^3)) - (erfc(a r) / r^3) ] =
    //
    // -K ( [ -(2 a^2 + 1 / r^2) B exp(-a^2 r^2) ] +
    //      [ -B exp(-a^2 r^2) / r^2 - (erfc(a r) / r^3) - (erfc(a r) / r^3) ] =
    //
    // -K ( [ -(2 a^2 + 1 / r^2) B exp(-a^2 r^2) ] +
    //      [ -B exp(-a^2 r^2) / r^2 - 2 (erfc(a r) / r^3) ] =
    //
    // -K ( [ -2 (a^2 + 1/r^2) B exp(-a^2 r^2) - 2 (erfc(a r) / r) / r^2 ] =
    //
    // 2K ( [ (a^2 + 1/r^2) B exp(-a^2 r^2) + (erfc(a r) / r) / r^2 ] =
    return kcoul * value_two * (((ew_t2 + invr2) * exp_quant) + (u_quant * invr2));
  }
  else if (order == 3) {
    const Tcalc value_two = 2.0;
    const Tcalc value_three = 3.0;
    const Tcalc value_four = 4.0;
    const Tcalc ew_t2 = ew_coeff * ew_coeff;
    const Tcalc ew_t4 = pow(ew_coeff, 4.0);
    const Tcalc invr2 = value_one / r2;
    const Tcalc invr3 = invr * invr2;

    // (d/dr) [ (a^2 + 1/r^2) B exp(-a^2 r^2) ] =
    // (-2 / r^3) B exp(-a^2 r^2) - (2a^2 r) (a^2 + 1/r^2) B exp(-a^2 r^2) =
    // -2B exp(-a^2 r^2) (1/r^3 + (a^4)r + (a^2)/r)
    //
    // (d/dr) [ (erfc(a r) / r) / r^2 ] =
    // (-(B exp(a^2 r^2) + (erfc(a r) / r)) / r^3) - ((2 erfc(a r) / r) / r^3) =
    // 
    // => (d/dr) [ Electrostatic PME second derivative ] =
    // 2K ( (-2B exp(-a^2 r^2) (1/r^3 + (a^4)r + (a^2)/r)) -
    //      ((B exp(a^2 r^2) + 3 (erfc(a r) / r)) / r^3)
    return value_two * kcoul * ((-value_two * exp_quant * (invr3 + (ew_t4 * r) + (ew_t2 * invr))) -
                                (invr3 * (exp_quant + (value_three * u_quant))));
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename T3>
Tcalc elecPMEDirectSpace(const T3 pa, const T3 pb, const double ew_coeff, const Tcalc kcoul,
                         const FunctionLevel order) {
  const size_t ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (ct == double_type_index);
  const Tcalc dx = pb.x - pa.x;
  const Tcalc dy = pb.y - pa.y;
  const Tcalc dz = pb.z - pa.z;
  const Tcalc r2 = (dx * dx) + (dy * dy) + (dz * dz);
  const Tcalc r = (tcalc_is_double) ? sqrt(r2) : sqrtf(r2);

  // Evaluate the switch once to obtain the necessary derivatives.
  Tcalc du, d2u, d3u;
  switch (order) {
  case FunctionLevel::VALUE:
    return elecPMEDirectSpace(ew_coeff, kcoul, r, r2, 0);
  case FunctionLevel::DX:
  case FunctionLevel::DY:
  case FunctionLevel::DZ:
    du = elecPMEDirectSpace(ew_coeff, kcoul, r, r2, 1);
    break;
  case FunctionLevel::DXX:
  case FunctionLevel::DXY:
  case FunctionLevel::DXZ:
  case FunctionLevel::DYY:
  case FunctionLevel::DYZ:
  case FunctionLevel::DZZ:
    d2u = elecPMEDirectSpace(ew_coeff, kcoul, r, r2, 2);
    break;
  case FunctionLevel::DXXX:
  case FunctionLevel::DXXY:
  case FunctionLevel::DXXZ:
  case FunctionLevel::DXYY:
  case FunctionLevel::DXYZ:
  case FunctionLevel::DXZZ:
  case FunctionLevel::DYYY:
  case FunctionLevel::DYYZ:
  case FunctionLevel::DYZZ:
  case FunctionLevel::DZZZ:
    d3u = elecPMEDirectSpace(ew_coeff, kcoul, r, r2, 3);
    break;
  }

  // A second switch will obtain the partial derivatives using the radially symmetric function
  // derivative calculators.
  return radialPartialDerivative(du, d2u, d3u, dx, dy, dz, r, r2, order);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename T3>
std::vector<Tcalc> elecPMEDirectSpace(const T3 pa, const T3 pb, const double ew_coeff,
                                      const Tcalc kcoul,
                                      const std::vector<FunctionLevel> &orders) {
  const size_t ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (ct == double_type_index);
  const Tcalc dx = pb.x - pa.x;
  const Tcalc dy = pb.y - pa.y;
  const Tcalc dz = pb.z - pa.z;
  const Tcalc r2 = (dx * dx) + (dy * dy) + (dz * dz);
  const Tcalc r = (tcalc_is_double) ? sqrt(r2) : sqrtf(r2);

  // Evaluate the switch for all orders to obtain the necessary derivatives.
  Tcalc u, du, d2u, d3u;
  bool u_found = false;
  bool du_found = false;
  bool d2u_found = false;
  bool d3u_found = false;
  const size_t n_order = orders.size();
  for (size_t i = 0; i < n_order; i++) {
    switch (orders[i]) {
    case FunctionLevel::VALUE:
      if (u_found == false) {
        u = elecPMEDirectSpace(ew_coeff, kcoul, r, r2, 0);
        u_found = true;
      }
      break;
    case FunctionLevel::DX:
    case FunctionLevel::DY:
    case FunctionLevel::DZ:
      if (du_found == false) {
        du = elecPMEDirectSpace(ew_coeff, kcoul, r, r2, 1);
        du_found = true;
      }
      break;
    case FunctionLevel::DXX:
    case FunctionLevel::DXY:
    case FunctionLevel::DXZ:
    case FunctionLevel::DYY:
    case FunctionLevel::DYZ:
    case FunctionLevel::DZZ:
      if (d2u_found == false) {
        d2u = elecPMEDirectSpace(ew_coeff, kcoul, r, r2, 2);
        d2u_found = true;
      }
      break;
    case FunctionLevel::DXXX:
    case FunctionLevel::DXXY:
    case FunctionLevel::DXXZ:
    case FunctionLevel::DXYY:
    case FunctionLevel::DXYZ:
    case FunctionLevel::DXZZ:
    case FunctionLevel::DYYY:
    case FunctionLevel::DYYZ:
    case FunctionLevel::DYZZ:
    case FunctionLevel::DZZZ:
      if (d3u_found == false) {
        d3u = elecPMEDirectSpace(ew_coeff, kcoul, r, r2, 3);
        d3u_found = true;
      }
      break;
    }
  }

  // Allocate the result
  std::vector<Tcalc> result(n_order);
  for (size_t i = 0; i < n_order; i++) {
    result[i] = radialPartialDerivative(du, d2u, d3u, dx, dy, dz, r, r2, orders[i]);
  }
  return result;
}

} // namespace energy
} // namespace stormm
