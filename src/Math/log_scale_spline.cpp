#include "copyright.h"
#include "log_scale_spline.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
ullint doublePrecisionSplineDetailMask(const int mantissa_bits) {
  ullint result = 0LLU;
  for (int i = 0; i < 52 - mantissa_bits; i++) {
    result |= (0x1LLU << i);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
uint singlePrecisionSplineDetailMask(const int mantissa_bits) {
  uint result = 0U;
  for (int i = 0; i < 23 - mantissa_bits; i++) {
    result |= (0x1U << i);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double doublePrecisionSplineArgument(const double arg, const int mantissa_bits,
                                     const ullint detail_mask) {
  Ecumenical8 xfrm = { .d = arg };
  xfrm.ulli = (((xfrm.ulli & detail_mask) << mantissa_bits) | 0x3ff0000000000000ULL);
  return xfrm.d;
}

//-------------------------------------------------------------------------------------------------
float singlePrecisionSplineArgument(const float arg, const int mantissa_bits,
                                    const uint detail_mask) {
  Ecumenical4 xfrm = { .f = arg };
  xfrm.ui = (((xfrm.ui & detail_mask) << mantissa_bits) | 0x3f800000U);
  return xfrm.f;
}

} // namespace stmath
} // namespace stormm
