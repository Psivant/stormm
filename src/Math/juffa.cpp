#include <cmath>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "juffa.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
JfErfc::JfErfc(const float alpha_in) :
  alpha{alpha_in}
{
  c0 = -1.6488499458192755E-006 * pow(alpha, 10.0);
  c1 =  2.9524665006554534E-005 * pow(alpha,  9.0);
  c2 = -2.3341951153749626E-004 * pow(alpha,  8.0);
  c3 =  1.0424943374047289E-003 * pow(alpha,  7.0);
  c4 = -2.5501426008983853E-003 * pow(alpha,  6.0);
  c5 =  3.1979939710877236E-004 * pow(alpha,  5.0);
  c6 =  2.7605379075746249E-002 * pow(alpha,  4.0);
  c7 = -1.4827402067461906E-001 * pow(alpha,  3.0);
  c8 = -9.1844764013203406E-001 * alpha * alpha;
  c9 = -1.6279070384382459E+000 * alpha;
}

//-------------------------------------------------------------------------------------------------
float JfErfc::erfcJf(const float x) {

  // Approximate log(erfc(x)) with rel. error < 7e-9
  float t;
  t = c0;
  t = (t * x) + c1;
  t = (t * x) + c2;
  t = (t * x) + c3;
  t = (t * x) + c4;
  t = (t * x) + c5;
  t = (t * x) + c6;
  t = (t * x) + c7;
  t = (t * x) + c8;
  t = (t * x) + c9;
  t = t * x;
  return exp2f(t);
}

//-------------------------------------------------------------------------------------------------
float expJf(const float x) {

  // Compute exp(x) = 2**i * exp(f); i = rintf (x / log(2))
  float j = (1.442695f * x) + 12582912.0f - 12582912.0f; // 0x1.715476p0, 0x1.8p23
  float f = (j * -0.693145752f) + x;                     // -0x1.62e400p-1  // log_2_hi 
  f = (j * -1.42860677e-6f) + f;                         // -0x1.7f7d1cp-20 // log_2_lo 
  const int i = static_cast<int>(j);

  // Approximate r = exp(f) on interval [-log(2)/2, +log(2)/2]
  float r = 0.00137805939f;      // 0x1.694000p-10
  r = (r * f) + 0.00837312452f;  // 0x1.125edcp-7
  r = (r * f) + 0.0416695364f;   // 0x1.555b5ap-5
  r = (r * f) + 0.166664720f;    // 0x1.555450p-3
  r = (r * f) + 0.499999851f;    // 0x1.fffff6p-2
  r = (r * f) + 1.00000000f;     // 0x1.000000p+0
  r = (r * f) + 1.00000000f;     // 0x1.000000p+0

  // exp(a) = 2**i * r
  const int ix = (i > 0) ? 0 : 0x83000000;
  Ecumenical4 workspc;
  workspc.i = 0x7f000000 + ix;
  const float s = workspc.f;
  workspc.i = (i << 23) - ix;
  const float t = workspc.f;
  r = r * s;
  r = r * t;

  // Handle special cases: severe overflow / underflow, then return
  if (fabsf(x) >= 104.0f) {
    r = s * s;
  }
  return r;
}

} // namespace stmath
} // namespacee stormm
