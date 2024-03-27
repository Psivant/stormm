// -*-c++-*-
#include "copyright.h"
#include "radial_derivatives.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename T> T radialFirstDerivative(const T dfunc, const T disp, const T r) {

  // The first derivatives follow from a simple application of the chain rule.          
  // (d/dx) [u(r)] = (du/dr)(dr/dx) = (du/dr)(x / r)
  if (r < constants::tiny) {
    return dfunc * disp / constants::tiny;
  }
  else {
    return dfunc * disp / r;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T radialSecondDerivative(const T dfunc, const T ddfunc, const T disp,
                                               const T r) {

  // The matrix of second partial derivatives is more complicated, involving applications of the
  // product and quotient rules of derivatives.
  //
  // (d/dx) [ (du/dr)(x / r) ] = (d/dx) [ (du/dr) ] (x / r) + (du/dr) [ (d/dx) [ x / r ] ]
  //                           = (d2u / dr2)(dr / dx) (x / r) + (du/dr) [ (d/dx) [ x / r ] ]
  //                           = (d2u / dr2)(x / r)^2 + (du/dr) [ (d/dx) [ x / r ] ]
  //                           = (d2u/dr2)(x / r)^2 + (du/dr) [ (r - (x/r)(x)) / r^2 ]
  //                           = (d2u/dr2)(x / r)^2 + (du/dr) [ (r - (x^2 / r)) / r^2 ]
  //                           = (d2u/dr2)(x / r)^2 + (du/dr) [ (1 - (x / r)^2) / r ]
  if (r < constants::small) {
    const T dpr_sq = (disp / constants::small) * (disp / constants::small);
    return (ddfunc * dpr_sq) + (dfunc * (static_cast<T>(1.0) - dpr_sq) / constants::small);    
  }
  else {
    const T dpr_sq = (disp / r) * (disp / r);
    return (ddfunc * dpr_sq) + (dfunc * (static_cast<T>(1.0) - dpr_sq) / r);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T radialSecondDerivative(const T dfunc, const T ddfunc, const T disp_x,
                                               const T disp_y, const T r, const T r2) {

  // Mixed partial derivatives are generally easier to compute than partial derivatives with
  // multiple differentiations along the same direction.
  //
  // (d/dy) [ (d/dx) [ u(r) ] ] = (d/dy) [ (du/dr)(x / r) ]
  //                            = (d/dy) [ (du/dr) ] (x / r) +
  ///                             (d/dy) [ (x / r) ] (du/dr)
  //                            = (d2u / dr2) (dr / dy) (x / r) +
  //                              (d / dr) (x / r) (dr / dy) (du / dr)
  //                            = (d2u / dr2) (y / r) (x / r) +
  //                              (-x / r^2) (y / r) (du / dr)
  //                            = (d2u / dr2)(x * y / r^2) - (du / dr)(x * y / r^3)
  //                            = (x * y) [ (d2u / dr2)(1 / r^2) - (du / dr)(1 / r^3) ]
  if (r < constants::small) {
    return disp_x * disp_y * (ddfunc - (dfunc / constants::small)) / (constants::small *
                                                                      constants::small);
  }
  else {
    const T use_r2 = (r2 >= constants::small) ? r2 : r * r;
    return disp_x * disp_y * (ddfunc - (dfunc / r)) / use_r2;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T radialThirdDerivative(const T dfunc, const T ddfunc, const T dddfunc,
                                              const T disp, const T r, const T r2) {

  // The tensor of third partial derivatives is even more complicated, but there are only three
  // base cases to solve: (d3u / dx3), (d3u / dx2y), and (d3u / dxyz).  All other mixed partial
  // derivatives can be found by substituting different dimensions.
  //
  // (d/dx) [ (d2u/dr2)(x / r)^2 + (du/dr) [ (1 - (x / r)^2) / r ] ]
  //   = (d3u / dr3)(x / r)^3 +
  //     (d2u / dr2)[ (d/dx) [ x^2 / r^2 ] ] +
  //     ((1 - (x / r)^2) / r) (d/dx) [ (du / dr) ] +
  //     (du/dr) (d/dx) [ (1 - (x / r)^2) / r]
  //
  //   = (d3u / dr3)(x / r)^3 +
  //     (d2u / dr2) ((2x)(r^2) - (2r)(x / r)(x^2)) / (r^4) +
  //     (d2u / dr2) (x / r) [ (1 - (x / r)^2) / r ] +
  //     (du/dr) ( (d/dx) [ (1 / r) ] - (d/dx) [ x^2 / r^3 ] )
  //
  //   = (d3u / dr3)(x / r)^3 +
  //     (d2u / dr2) ((2 x r^2 - 2 x^3) / (r^4)) +
  //     (d2u / dr2) (x) ( (1 - (x / r)^2) / r^2 ) +
  //     (du/dr) ( (d/dx) [ (1 / r) ] - (d/dx) [ x^2 / r^3 ] )
  //
  //   = (d3u / dr3)(x / r)^3 +
  //     (d2u / dr2) (x) ((2 r^2 - 2 x^2) / (r^4)) +
  //     (d2u / dr2) (x) ( (r^2 / r^4) - (x^2 / r^4) ) +
  //     (du/dr) ( (d/dx) [ (1 / r) ] - (d/dx) [ x^2 / r^3 ] )
  //
  //   = (d3u / dr3)(x / r)^3 +
  //     (d2u / dr2) (x) ((3 r^2 - 3 x^2) / (r^4)) +
  //     (du/dr) ( (d/dx) [ (1 / r) ] - (d/dx) [ x^2 / r^3 ] )
  //
  //   = (d3u / dr3)(x / r)^3 +
  //     (d2u / dr2) (x) ((3 r^2 - 3 x^2) / (r^4)) +
  //     (du/dr) ((3 x^3 - 3 x r^2) / r^5)
  //
  //   = (d3u / dr3)(x / r)^3 +
  //     (d2u / dr2) (x) ((3 r^2 - 3 x^2) / (r^4)) +
  //     (du/dr) (x) ((3 x^2 - 3 r^2) / r^5)
  const T value_one   = 1.0;
  const T value_three = 3.0;
  const T invr = (r < constants::small) ? value_one / constants::small : value_one / r;
  const T invr2 = value_one / r2;
  const T disp2 = disp * disp;
  return disp * ((dddfunc * disp2) +
                 (ddfunc * ((value_three * r2) - (value_three * disp2)) * invr) +
                 (dfunc * ((value_three * disp2) - (value_three * r2)) *
                  invr2)) * invr2 * invr;
}

//-------------------------------------------------------------------------------------------------
template <typename T> T radialThirdDerivative(const T dfunc, const T ddfunc, const T dddfunc,
                                              const T disp_x, const T disp_y, const T r,
                                              const T r2) {

  // The two-and-one mixed partial is perhaps the most difficult to compute.
  // 
  // (d/dy) [ (d2u/dr2)(x / r)^2 + (du/dr) [ (1 - (x / r)^2) / r ] ]
  //   = (x^2 / r^2) (d/dy) [ (d2u / dr2) ] +
  //     (d2u / dr2) (d/dy) [ (x^2 / r^2) ] +
  //     ( (1 - (x / r)^2) / r ) (d/dy) [ (du / dr) ] +
  //     (du / dr) (d/dy) [ (1 - (x / r)^2) / r ]
  //
  //   = (x^2 / r^2) (y / r) (d3u / dr3) +
  //     (d2u / dr2) ( (-2y) (x^2) / r^4 ) +
  //     ( (1 - (x / r)^2) / r ) (y / r) (d2u / dr2) +
  //     (du / dr) ( (y) (3x^2 - r^2) / r^5 )
  //
  //   = (d3u / dr3) (x * x * y / r^3) +
  //     (d2u / dr2) ( (1 - (x / r)^2) / r ) (y / r) - (2y) (x^2) / r^4 ) +
  //     (du / dr) ( (y) (3x^2 - r^2) / r^5 )
  //
  //   = (d3u / dr3) (x * x * y / r^3) +
  //     (d2u / dr2) ( (y) ((1 / r)  - (x * x / r^3)) / r ) - (2y) (x^2) / r^4 ) +
  //     (du / dr) ( (y) (3x^2 - r^2) / r^5 )
  //
  //   = (d3u / dr3) (x * x * y / r^3) +
  //     (d2u / dr2) ( (y) ((1 / r)  - (x * x / r^3)) / r ) - (2y) (x^2) / r^4 ) +
  //     (du / dr) ( (y) (3x^2 - r^2) / r^5 )
  //
  //   = (d3u / dr3) (x * x * y / r^3) +
  //     (d2u / dr2) ( (y) ((1 / r^2)  - (x * x / r^4)) - (2y) (x^2) / r^4 ) +
  //     (du / dr) ( (y) (3x^2 - r^2) / r^5 )
  //
  //   = (d3u / dr3) (x * x * y / r^3) +
  //     (d2u / dr2) (y) ( (r^2 / r^4)  - (x * x / r^4) - (2 x * x) / r^4 ) +
  //     (du / dr) ( (y) (3x^2 - r^2) / r^5 )
  //
  //   = (d3u / dr3) (x * x * y / r^3) +
  //     (d2u / dr2) (y) ( (r^2 - (3 * x * x) ) / r^4) +
  //     (du / dr) ( (y) (3x^2 - r^2) / r^5 )
  const T value_one   = 1.0;
  const T value_three = 3.0;
  const T invr = (r < constants::small) ? value_one / constants::small : value_one / r;
  const T invr2 = value_one / r2;
  const T disp_x2 = disp_x * disp_x;
  return disp_y * ((dddfunc * disp_x2) + (ddfunc * (r2 - (value_three * disp_x2)) * invr) + 
                   (dfunc * ((value_three * disp_x2) - r2) * invr2)) * invr2 * invr;
}

//-------------------------------------------------------------------------------------------------
template <typename T> T radialThirdDerivative(const T dfunc, const T ddfunc, const T dddfunc,
                                              const T disp_x, const T disp_y, const T disp_z,
                                              const T r, const T r2) {

  // The full mixed partial, where differentiation occurs only once in each direction, is the
  // simplest to compute.
  //
  // (d/dx) [ (y * z) [ (d2u / dr2)(1 / r^2) - (du / dr)(1 / r^3) ] ]
  //   = (y * z) (d/dx) [ (d2u / dr2)(1 / r^2) - (du / dr)(1 / r^3) ]
  //   = (x * y * z) [ (1 / r^3) (d3u / dr3) - (d2u / dr2) (3 / r^4) + (du / dr) (3 / r^5) ]
  const T value_one   = 1.0;
  const T value_three = 3.0;
  const T invr = (r < constants::small) ? value_one / constants::small : value_one / r;
  const T invr2 = value_one / r2;
  return (disp_x * disp_y * disp_z * invr2 * invr) *
         (dddfunc - (value_three * ddfunc * invr) + (value_three * dfunc * invr2));
}

//-------------------------------------------------------------------------------------------------
template <typename T> T radialPartialDerivative(const T du, const T d2u, const T d3u, const T dx,
                                                const T dy, const T dz, const T r, const T r2,
                                                const FunctionLevel order) {
  switch (order) {
  case FunctionLevel::VALUE:
    rtErr("The value of the function is not a valid derivative.", "radialPartialDerivative");
  case FunctionLevel::DX:
    return radialFirstDerivative(du, dx, r);
  case FunctionLevel::DY:
    return radialFirstDerivative(du, dy, r);
  case FunctionLevel::DZ:
    return radialFirstDerivative(du, dz, r);
  case FunctionLevel::DXX:
    return radialSecondDerivative(du, d2u, dx, r);
  case FunctionLevel::DXY:
    return radialSecondDerivative(du, d2u, dx, dy, r, r2);
  case FunctionLevel::DXZ:
    return radialSecondDerivative(du, d2u, dx, dz, r, r2);
  case FunctionLevel::DYY:
    return radialSecondDerivative(du, d2u, dy, r);
  case FunctionLevel::DYZ:
    return radialSecondDerivative(du, d2u, dy, dz, r, r2);
  case FunctionLevel::DZZ:
    return radialSecondDerivative(du, d2u, dz, r);
  case FunctionLevel::DXXX:
    return radialThirdDerivative(du, d2u, d3u, dx, r, r2);
  case FunctionLevel::DXXY:
    return radialThirdDerivative(du, d2u, d3u, dx, dy, r, r2);
  case FunctionLevel::DXXZ:
    return radialThirdDerivative(du, d2u, d3u, dx, dz, r, r2);
  case FunctionLevel::DXYY:
    return radialThirdDerivative(du, d2u, d3u, dy, dx, r, r2);
  case FunctionLevel::DXYZ:
    return radialThirdDerivative(du, d2u, d3u, dx, dy, dz, r, r2);
  case FunctionLevel::DXZZ:
    return radialThirdDerivative(du, d2u, d3u, dz, dx, r, r2);
  case FunctionLevel::DYYY:
    return radialThirdDerivative(du, d2u, d3u, dy, r, r2);
  case FunctionLevel::DYYZ:
    return radialThirdDerivative(du, d2u, d3u, dy, dz, r, r2);
  case FunctionLevel::DYZZ:
    return radialThirdDerivative(du, d2u, d3u, dz, dy, r, r2);
  case FunctionLevel::DZZZ:
    return radialThirdDerivative(du, d2u, d3u, dz, r, r2);
  }
  __builtin_unreachable();
}

} // namespace stmath
} // namespace stormm
