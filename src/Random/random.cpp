#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include "copyright.h"
#include "Constants/symbol_values.h"
#include "Parsing/polynumeric.h"
#include "Reporting/error_format.h"
#include "random.h"

namespace stormm {
namespace random {

using parse::PolyNumeric;

//-------------------------------------------------------------------------------------------------
Ran2Generator::Ran2Generator(const int igseed) :
  iunstable_a{-1},
  iunstable_b{123456789},
  state_sample{0},
  state_vector{0}
{
  // Initialize by calling the uniform generator
  uniformRandomNumber();

  // Compute additional random numbers based on a pseudo-random seed.  In order to make this
  // generator diverge from others based on different seeds, we must first run through about
  // sixteen cycles of the uniform generator in order to run through correlated results left
  // over after the common initialization.  The entire state vector needs to change in order
  // to produce a decorrelated random number.
  iunstable_a = igseed;
  for (int i = 0; i < 15; i++) {
    uniformRandomNumber();
  }
}

//-------------------------------------------------------------------------------------------------
double Ran2Generator::uniformRandomNumber() {

  // Local variable to avoid dereferencing the self pointer
  int lcl_iunstbl_a  = iunstable_a;
  int lcl_iunstbl_b  = iunstable_b;
  int lcl_state_sample = state_sample;

  // Populate the state vector
  int unstbl_quotient;
  if (lcl_iunstbl_a <= 0) {
    lcl_iunstbl_a = (-lcl_iunstbl_a < 1) ? 1 : -lcl_iunstbl_a;
    lcl_iunstbl_b = lcl_iunstbl_a;
    for (int j = ntab + 7; j >= 0 ; j--) {
      unstbl_quotient = lcl_iunstbl_a / iq1;
      lcl_iunstbl_a = ia1 * (lcl_iunstbl_a - unstbl_quotient * iq1) - (unstbl_quotient * ir1);
      lcl_iunstbl_a += (lcl_iunstbl_a < 0) * im1;
      state_vector[j] = lcl_iunstbl_a;
    }
    lcl_state_sample = state_vector[0];
  }
  unstbl_quotient = lcl_iunstbl_a / iq1;
  lcl_iunstbl_a = ia1 * (lcl_iunstbl_a - unstbl_quotient * iq1) - unstbl_quotient*ir1;
  lcl_iunstbl_a += (lcl_iunstbl_a < 0) * im1;
  unstbl_quotient = lcl_iunstbl_b / iq2;
  lcl_iunstbl_b = ia2 * (lcl_iunstbl_b - unstbl_quotient * iq2) - unstbl_quotient * ir2;
  lcl_iunstbl_b += (lcl_iunstbl_b < 0) * im2;
  int swap_index = lcl_state_sample / ndiv;
  lcl_state_sample = state_vector[swap_index] - lcl_iunstbl_b;
  state_vector[swap_index] = lcl_iunstbl_a;

  // Generate the alternative number so as to return the lesser of this or rand2_max
  lcl_state_sample += (lcl_state_sample < 1) * im1_minus_one;
  double temp = am * lcl_state_sample;

  // Copy lcl values back to modifiable struct member variable
  iunstable_a = lcl_iunstbl_a;
  iunstable_b = lcl_iunstbl_b;
  state_sample = lcl_state_sample;

  return std::min(ran2_max, temp);
}

//-------------------------------------------------------------------------------------------------
float Ran2Generator::spUniformRandomNumber() {
  return uniformRandomNumber();
}
  
//-------------------------------------------------------------------------------------------------
double Ran2Generator::gaussianRandomNumber() {
  const double x1 = sqrt(-2.0 * log(1.0 - uniformRandomNumber()));
  const double x2 = sin(symbols::twopi * uniformRandomNumber());
  return x1 * x2;
}

//-------------------------------------------------------------------------------------------------
float Ran2Generator::spGaussianRandomNumber() {
  const float x1 = sqrtf(-2.0f * static_cast<float>(log(1.0 - uniformRandomNumber())));
  const float x2 = sinf(symbols::twopi_f * uniformRandomNumber());
  return x1 * x2;
}

//-------------------------------------------------------------------------------------------------
Xoroshiro128pGenerator::Xoroshiro128pGenerator(const int igseed, const int niter) :
    state{seed128(igseed)}
{
  // Trap cases where ullint is not 64 bit (this should never be an issue)
  if (sizeof(ullint) != 8) {
    rtErr("The size of an unsigned long long integer (ullint) must be 64 bits in order to "
          "instantiate the xoroshiro128+ random number generator.", "Xoroshiro128pGenerator");
  }

  // Trap cases where double is not 64 bit or float is not 32 bit (this should never be an issue)
  if (sizeof(double) != 8) {
    rtErr("The size of a double-precision real number (double) must by 64 bits in order to "
          "instantiate the xoroshiro128+ random number generator.", "Xoroshiro128pGenerator");
  }
  if (sizeof(float) != 4) {
    rtErr("The size of a single-precision real number (float) must by 32 bits in order to "
          "instantiate the xoroshiro128+ random number generator.", "Xoroshiro128pGenerator");
  }

  // Run some iterations to get the bit occupancy to roughly half zeros and half ones
  for (int i = 0; i < niter; i++) {
    next();
  }
}

//-------------------------------------------------------------------------------------------------
Xoroshiro128pGenerator::Xoroshiro128pGenerator(const ullint2 state_in) :
  state{state_in}
{}
    
//-------------------------------------------------------------------------------------------------
double Xoroshiro128pGenerator::uniformRandomNumber() {
  const ullint rndbits = next();
  return static_cast<double>(rndbits >> 11) * rng_unit_bin_offset;
}

//-------------------------------------------------------------------------------------------------
float Xoroshiro128pGenerator::spUniformRandomNumber() {
  const ullint rndbits = next();
  return static_cast<float>(rndbits >> 40) * rng_unit_bin_offset_f;
}

//-------------------------------------------------------------------------------------------------
double Xoroshiro128pGenerator::gaussianRandomNumber() {
  const double x1 = sqrt(-2.0 * log(1.0 - uniformRandomNumber()));
  const double x2 = sin(symbols::twopi * uniformRandomNumber());
  return x1 * x2;
}

//-------------------------------------------------------------------------------------------------
float Xoroshiro128pGenerator::spGaussianRandomNumber() {
  const float x1 = sqrtf(-2.0f * static_cast<float>(log(1.0 - uniformRandomNumber())));
  const float x2 = sinf(symbols::twopi_f * spUniformRandomNumber());
  return x1 * x2;
}

//-------------------------------------------------------------------------------------------------
void Xoroshiro128pGenerator::jump() {
  const ullint2 stride = { xrs128p_jump_i, xrs128p_jump_ii };
  fastForward(stride);
}

//-------------------------------------------------------------------------------------------------
void Xoroshiro128pGenerator::longJump() {
  const ullint2 stride = { xrs128p_longjump_i, xrs128p_longjump_ii };
  fastForward(stride);
}

//-------------------------------------------------------------------------------------------------
ullint2 Xoroshiro128pGenerator::revealState() const {
  return state;
}

//-------------------------------------------------------------------------------------------------
ullint Xoroshiro128pGenerator::revealBitString() const {
  return state.x + state.y;
}

//-------------------------------------------------------------------------------------------------
void Xoroshiro128pGenerator::setState(const ullint2 state_in) {
  state = state_in;
}

//-------------------------------------------------------------------------------------------------
ullint2 Xoroshiro128pGenerator::seed128(const int igseed) {

  // It is within the C and C++ standards to set an unsigned integer equal to a signed integer,
  // even a negative value.  It changes the interpretation, no more.
  const ullint uiseed = igseed;
  const int nbits = sizeof(uint) * 8;
  const int nfill = sizeof(ullint) * 8;
  const int nspc  = nfill / nbits;
  int offset = 0;
  ullint2 result = { 0LLU, 0LLU };
  for (int i = 0; i < nbits; i++) {
    result.x |= (((uiseed >> i) & 0x1LLU) << ((nspc * i) + offset));
    offset++;
    offset *= (offset != nspc);
    result.y |= (((uiseed >> i) & 0x1LLU) << ((nspc * i) + offset));
    offset++;
    offset *= (offset != nspc);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
ullint Xoroshiro128pGenerator::next() {
  const ullint s0 = state.x;
  ullint       s1 = state.y;
  const ullint result = s0 + s1;
  s1 ^= s0;
  state.x = (((s0 << 24) | (s0 >> (64 - 24))) ^ s1 ^ (s1 << 16));
  state.y =  ((s1 << 37) | (s1 >> (64 - 37)));
  return result;
}

//-------------------------------------------------------------------------------------------------
void Xoroshiro128pGenerator::fastForward(const ullint2 stride) {
  ullint s0 = 0LLU;
  ullint s1 = 0LLU;
  for (int b = 0; b < 64; b++) {
    if (stride.x & (0x1LLU << b)) {
      s0 ^= state.x;
      s1 ^= state.y;
    }
    next();
  }
  for (int b = 0; b < 64; b++) {
    if (stride.y & (0x1LLU << b)) {
      s0 ^= state.x;
      s1 ^= state.y;
    }
    next();
  }
  state.x = s0;
  state.y = s1;
}

//-------------------------------------------------------------------------------------------------
Xoshiro256ppGenerator::Xoshiro256ppGenerator(const int igseed, const int niter) :
    state{seed256(igseed)}
{
  // Trap cases where ullint is not 64 bit (this should never be an issue)
  if (sizeof(ullint) != 8) {
    rtErr("The size of an unsigned long long integer (ullint) must be 64 bits in order to "
          "instantiate the xoshiro256++ random number generator.", "Xoshiro256ppGenerator");
  }

  // Trap cases where double is not 64 bit or float is not 32 bit (this should never be an issue)
  if (sizeof(double) != 8) {
    rtErr("The size of a double-precision real number (double) must by 64 bits in order to "
          "instantiate the xoroshiro256++ random number generator.", "Xoshiro256ppGenerator");
  }
  if (sizeof(float) != 4) {
    rtErr("The size of a single-precision real number (float) must by 32 bits in order to "
          "instantiate the xoroshiro256++ random number generator.", "Xoshiro256ppGenerator");
  }

  // Run some iterations to get the bit occupancy to roughly half zeros and half ones
  for (int i = 0; i < niter; i++) {
    next();
  }
}

//-------------------------------------------------------------------------------------------------
Xoshiro256ppGenerator::Xoshiro256ppGenerator(const ullint4 state_in) :
    state{state_in}
{}

//-------------------------------------------------------------------------------------------------
double Xoshiro256ppGenerator::uniformRandomNumber() {
  const ullint rndbits = next();
  return static_cast<double>(rndbits >> 11) * rng_unit_bin_offset;
}

//-------------------------------------------------------------------------------------------------
float Xoshiro256ppGenerator::spUniformRandomNumber() {
  const ullint rndbits = next();
  return static_cast<float>(rndbits >> 40) * rng_unit_bin_offset_f;
}

//-------------------------------------------------------------------------------------------------
double Xoshiro256ppGenerator::gaussianRandomNumber() {
  const double x1 = sqrt(-2.0 * log(1.0 - uniformRandomNumber()));
  const double x2 = sin(symbols::twopi * uniformRandomNumber());
  return x1 * x2;
}

//-------------------------------------------------------------------------------------------------
float Xoshiro256ppGenerator::spGaussianRandomNumber() {
  const float x1 = sqrtf(-2.0f * static_cast<float>(log(1.0 - spUniformRandomNumber())));
  const float x2 = sinf(symbols::twopi_f * spUniformRandomNumber());
  return x1 * x2;
}

//-------------------------------------------------------------------------------------------------
void Xoshiro256ppGenerator::jump() {
  const ullint4 stride = { xrs256pp_jump_i,   xrs256pp_jump_ii,
                           xrs256pp_jump_iii, xrs256pp_jump_iv };
  fastForward(stride);
}

//-------------------------------------------------------------------------------------------------
void Xoshiro256ppGenerator::longJump() {
  const ullint4 stride = { xrs256pp_longjump_i,   xrs256pp_longjump_ii,
                           xrs256pp_longjump_iii, xrs256pp_longjump_iv };
  fastForward(stride);
}

//-------------------------------------------------------------------------------------------------
ullint4 Xoshiro256ppGenerator::revealState() const {
  return state;
}

//-------------------------------------------------------------------------------------------------
ullint Xoshiro256ppGenerator::revealBitString() const {
  const ullint sxsw = state.x + state.w;
  return state.x + ((sxsw << 23) | (sxsw >> (64 - 23)));
}

//-------------------------------------------------------------------------------------------------
void Xoshiro256ppGenerator::setState(const ullint4 state_in) {
  state = state_in;
}

//-------------------------------------------------------------------------------------------------
ullint4 Xoshiro256ppGenerator::seed256(const int igseed) {

  // It is within the C and C++ standards to set an unsigned integer equal to a signed integer,
  // even a negative value.  It changes the interpretation, no more.
  const ullint uiseed = igseed;
  const ullint ujseed = (igseed ^ -1);
  const int nbits = sizeof(uint) * 8;
  const int nfill = sizeof(ullint) * 8;
  const int nspc  = nfill / nbits;
  ullint4 result = { 0LLU, 0LLU, 0LLU, 0LLU };
  for (int i = 0; i < nbits; i++) {
    result.x |= (((uiseed >> i) & 0x1LLU) << ((nspc * i)));
    result.y |= (((uiseed >> i) & 0x1LLU) << ((nspc * i) + 1));
    result.z |= (((ujseed >> i) & 0x1LLU) << ((nspc * i)));
    result.w |= (((ujseed >> i) & 0x1LLU) << ((nspc * i) + 1));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
ullint Xoshiro256ppGenerator::next() {
  const ullint sxsw = state.x + state.w;
  const ullint result = state.x + ((sxsw << 23) | (sxsw >> (64 - 23)));
  const ullint t = (state.y << 17);
  state.z ^= state.x;
  state.w ^= state.y;
  state.y ^= state.z;
  state.x ^= state.w;
  state.z ^= t;
  state.w = ((state.w << 45) | (state.w >> (64 - 45)));
  return result;
}

//-------------------------------------------------------------------------------------------------
void Xoshiro256ppGenerator::fastForward(ullint4 stride) {
  ullint s0 = 0LLU;
  ullint s1 = 0LLU;
  ullint s2 = 0LLU;
  ullint s3 = 0LLU;
  for (int b = 0; b < 64; b++) {
    if (stride.x & (0x1LLU << b)) {
      s0 ^= state.x;
      s1 ^= state.y;
      s2 ^= state.z;
      s3 ^= state.w;
    }
    next();
  }
  for (int b = 0; b < 64; b++) {
    if (stride.y & (0x1LLU << b)) {
      s0 ^= state.x;
      s1 ^= state.y;
      s2 ^= state.z;
      s3 ^= state.w;
    }
    next();
  }
  for (int b = 0; b < 64; b++) {
    if (stride.z & (0x1LLU << b)) {
      s0 ^= state.x;
      s1 ^= state.y;
      s2 ^= state.z;
      s3 ^= state.w;
    }
    next();
  }
  for (int b = 0; b < 64; b++) {
    if (stride.w & (0x1LLU << b)) {
      s0 ^= state.x;
      s1 ^= state.y;
      s2 ^= state.z;
      s3 ^= state.w;
    }
    next();
  }
  state.x = s0;
  state.y = s1;
  state.z = s2;
  state.w = s3;
}

//-------------------------------------------------------------------------------------------------
void initXoroshiro128pArray(ullint2* state_vector, const int n_generators, const int igseed,
                            const int scrub_cycles) {
  Xoroshiro128pGenerator prng(igseed, scrub_cycles);
  const int n_seeds = std::min(max_xo_long_jumps, n_generators);
  for (int i = 0; i < n_seeds; i++) {
    state_vector[i] = prng.revealState();
    for (int j = n_seeds + i; j < n_generators; j += n_seeds) {
      prng.jump();
      state_vector[j] = prng.revealState();
    }
    prng.setState(state_vector[i]);
    prng.longJump();
  }
}

//-------------------------------------------------------------------------------------------------
void initXoroshiro128pArray(std::vector<ullint2> *state_vector, const int igseed,
                            const int scrub_cycles) {
  initXoroshiro128pArray(state_vector->data(), state_vector->size(), igseed, scrub_cycles);
}

//-------------------------------------------------------------------------------------------------
void initXoroshiro128pArray(Hybrid<ullint2> *state_vector, const int igseed,
                            const int scrub_cycles) {
  initXoroshiro128pArray(state_vector->data(), state_vector->size(), igseed, scrub_cycles);
}

//-------------------------------------------------------------------------------------------------
void initXoshiro256ppArray(ullint2* state_xy, ullint2* state_zw, const int n_generators,
                           const int igseed, const int scrub_cycles) {
  Xoshiro256ppGenerator prng(igseed, scrub_cycles);
  const int n_seeds = std::min(max_xo_long_jumps, n_generators);
  for (int i = 0; i < n_seeds; i++) {
    const ullint4 seeded_state = prng.revealState();
    state_xy[i] = { seeded_state.x, seeded_state.y };
    state_zw[i] = { seeded_state.z, seeded_state.w };
    for (int j = n_seeds + i; j < n_generators; j += n_seeds) {
      prng.jump();
      const ullint4 jumped_state = prng.revealState();
      state_xy[j] = { jumped_state.x, jumped_state.y };
      state_zw[j] = { jumped_state.z, jumped_state.w };
    }
    prng.setState(seeded_state);
    prng.longJump();
  }
}

//-------------------------------------------------------------------------------------------------
void initXoshiro256ppArray(std::vector<ullint2> *state_xy, std::vector<ullint2> *state_zw,
                           const int igseed, const int scrub_cycles) {
  initXoshiro256ppArray(state_xy->data(), state_zw->data(), state_xy->size(), igseed,
                        scrub_cycles);
}

//-------------------------------------------------------------------------------------------------
void initXoshiro256ppArray(Hybrid<ullint2> *state_xy, Hybrid<ullint2> *state_zw, const int igseed,
                           const int scrub_cycles) {
  initXoshiro256ppArray(state_xy->data(), state_zw->data(), state_xy->size(), igseed,
                        scrub_cycles);
}

} // namespace random
} // namespace stormm
