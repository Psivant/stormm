// -*-c++-*-
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/bspline.h"
#include "../../src/Math/summation.h"
#include "../../src/Math/hpc_summation.cuh"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Random/random.h"
#include "../../src/Random/hpc_random.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::card::HpcConfig;
using stormm::card::Hybrid;
using stormm::card::HybridFormat;
using stormm::card::HybridTargetLevel;
using stormm::card::GpuDetails;
using stormm::data_types::int95_t;
using stormm::data_types::isFloatingPointScalarType;
using stormm::data_types::llint;
using stormm::data_types::ullint;
using stormm::data_types::ullint2;
using stormm::data_types::ullint4;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::errors::rtWarn;
using stormm::random::max_xo_long_jumps;
using stormm::random::initXoroshiro128pArray;
using stormm::random::initXoshiro256ppArray;
using stormm::random::Ran2Generator;
using stormm::random::Xoroshiro128pGenerator;
using stormm::random::Xoshiro256ppGenerator;
using stormm::random::default_xoroshiro128p_scrub;
using stormm::random::default_xoshiro256pp_scrub;
using namespace stormm::constants;
using namespace stormm::stmath;
using namespace stormm::hpc_math;
using namespace stormm::testing;
using namespace stormm::numerics;

#include "../../src/Numerics/accumulation.cui"
#include "../../src/Math/bspline.cui"

//-------------------------------------------------------------------------------------------------
// Enumerator for specifying the operator to be used during testing Split Fixed Precision Math
//
//   Operator:        The operator to be used (summation or subtraction)
//-------------------------------------------------------------------------------------------------
enum class Operator {
  SUM,        ///< Split Fixed Integer Summation
  SUBTRACT    ///< Split Fixed Integer Subtraction
};

//-------------------------------------------------------------------------------------------------
// Enumerator for specifying the types of pairs to be used during Split Fixed Math operations
//
//   FusionStyle: Specifies the type of Split Fixed Precision numbers each of the operands are
//-------------------------------------------------------------------------------------------------
enum class FusionStyle {
  FUSED_DECPL,  ///< The first number is fused, and the second number is a pair of free integers.
  FUSED_FLOAT,  ///< The first number is fused, and the second number is fused in float type.
  DECPL_DECPL,  ///< Both numbers are pairs of free integers.
  DECPL_FLOAT,  ///< The first number is a pair of free integers, and the second number is float.
  FUSED_FUSED   ///< Both the numbers are fused.
};

//-------------------------------------------------------------------------------------------------
// Helper functions to return strings versions of enumerators above.
//-------------------------------------------------------------------------------------------------
std::string getOperatorName(const Operator& op) {
  switch (op) {
    case Operator::SUBTRACT: return "SUBTRACT";
    case Operator::SUM: return "SUM";
  }
  return "";
}

std::string getFusionStyleName(const FusionStyle& fusion) {
  switch (fusion) {
    case FusionStyle::FUSED_DECPL: return "FUSED_DECPL";
    case FusionStyle::FUSED_FLOAT: return "FUSED_FLOAT";
    case FusionStyle::DECPL_DECPL: return "DECPL_DECPL";
    case FusionStyle::DECPL_FLOAT: return "DECPL_FLOAT";
    case FusionStyle::FUSED_FUSED: return "FUSED_FUSED";
  }
  return "";
}

//-------------------------------------------------------------------------------------------------
// Load a vector of random numbers using a CPU-based Ran2 generator.
//
// Arguments:
//   va:    The Hybrid object to load (modified and returned)
//   prng:  Ran2 generator object
//-------------------------------------------------------------------------------------------------
void loadRan2ByCPU(Hybrid<double> *va, Ran2Generator *prng) {
  double* vdata = va->data();
  const int nval = va->size();
  for (int i = 0; i < nval; i++) {
    vdata[i] = (prng->uniformRandomNumber() - 0.5) * sqrt(static_cast<double>(i + 1));
  }
}

//-------------------------------------------------------------------------------------------------
// Load a vector of random numbers using a CPU-based xoroshiro128+ generator.
//
// Arguments:
//   va:    The Hybrid object to load (modified and returned)
//   prng:  Xoroshiro128+ generator object
//-------------------------------------------------------------------------------------------------
void loadXoroshiro128pByCPU(Hybrid<double> *va, Xoroshiro128pGenerator *prng) {
  double* vdata = va->data();
  const int nval = va->size();
  for (int i = 0; i < nval; i++) {
    vdata[i] = (prng->uniformRandomNumber() - 0.5) * sqrt(static_cast<double>(i + 1));
  }
}

//-------------------------------------------------------------------------------------------------
// Load a vector of random numbers using a CPU-based xoshiro256++ generator.
//
// Arguments:
//   va:    The Hybrid object to load (modified and returned)
//   prng:  Xoroshiro128+ generator object
//-------------------------------------------------------------------------------------------------
void loadXoshiro256ppByCPU(Hybrid<double> *va, Xoshiro256ppGenerator *prng) {
  double* vdata = va->data();
  const int nval = va->size();
  for (int i = 0; i < nval; i++) {
    vdata[i] = (prng->uniformRandomNumber() - 0.5) * sqrt(static_cast<double>(i + 1));
  }
}

//-------------------------------------------------------------------------------------------------
// Evaluate random numbers from a plethora of xoroshiro128+ generator states.
//
// Arguments:
//   state_vector:  The vector of all states
//   n_generators:  Total number of states
//   samples:       Vector of states and iteration counts to watch out for
//   random_ouput:  Random numbers produced by the requested states and iterations
//   n_samples:     The number of samples (length of samples and random_output)
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kEvalXoroshiro128p(ullint2* state_vector, const int n_generators, const int2* samples,
                   double* random_output, const int n_samples) {
  __shared__ volatile int smax_iter;
  
  // Find the maximum iteration and cache the samples
  if (threadIdx.x == 0) {
    smax_iter = 0;
  }
  int max_iter = 0;
  for (int pos = threadIdx.x; pos < n_samples; pos += blockDim.x) {
    const int2 tmp_sample = samples[pos];
    if (tmp_sample.y > max_iter) {
      max_iter = tmp_sample.y;
    }
  }
#ifdef STORMM_USE_HIP
  const int init_shfl_stride = 32;
#else
  const int init_shfl_stride = 16;
#endif
  for (int i = init_shfl_stride; i > 0; i >>= 1) {
    const int max_comp = SHFL_DOWN(max_iter, i);
    if (max_comp > max_iter) {
      max_iter = max_comp;
    }
  }
  __syncthreads();
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (lane_idx == 0 && max_iter > 0) {
    atomicMax((int*)&smax_iter, max_iter);
  }
  __syncthreads();
  const int ginc = blockDim.x * gridDim.x;
  for (int pos = threadIdx.x + (blockIdx.x * blockDim.x); pos < n_generators; pos += ginc) {
    ullint2 tmp_state = state_vector[pos];
    for (int i = 0; i <= smax_iter; i++) {

      // Get a uniform random number in the range [0, 1).  Casting the unsigned long long int to
      // a (signed) long long int works as anything great than 2^3 just becomes negative.  The
      // bit string is unchanged.
      const ullint s0 = tmp_state.x;
      ullint       s1 = tmp_state.y;
      const ullint rndbits = s0 + s1;
      const llint work = (((rndbits >> 12) & 0xfffffffffffff) | 0x3ff0000000000000);
      const double rn_out = __longlong_as_double(work) - 1.0;

      // Horribly inefficient loop and access pattern, but this is just
      // a test program.  All of this is L1-cached, too.
      for (int j = 0; j < n_samples; j++) {
        const int2 tmp_sample = samples[j];
        if (pos == tmp_sample.x && i == tmp_sample.y) {
          random_output[j] = rn_out;
        }
      }

      // Push the state forward
      s1 ^= s0;
      tmp_state.x = (((s0 << 24) | (s0 >> (64 - 24))) ^ s1 ^ (s1 << 16));
      tmp_state.y =  ((s1 << 37) | (s1 >> (64 - 37)));      
    }

    // In a simulation, the state vector would be updated after computations
    // for this thread or atom are complete.
    state_vector[pos] = tmp_state;    
  }
}

//-------------------------------------------------------------------------------------------------
// Evaluate random numbers from a plethora of xoshiro256++ generator states.
//
// Arguments:
//   state_vector:  The vector of all states
//   n_generators:  Total number of states
//   samples:       Vector of states and iteration counts to watch out for
//   random_ouput:  Random numbers produced by the requested states and iterations
//   n_samples:     The number of samples (length of samples and random_output)
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kEvalXoshiro256pp(ullint2* state_xy, ullint2* state_zw, const int n_generators,
                  const int2* samples, double* random_output, const int n_samples) {
  __shared__ volatile int smax_iter;
  
  // Find the maximum iteration and cache the samples
  if (threadIdx.x == 0) {
    smax_iter = 0;
  }
  int max_iter = 0;
  for (int pos = threadIdx.x; pos < n_samples; pos += blockDim.x) {
    const int2 tmp_sample = samples[pos];
    if (tmp_sample.y > max_iter) {
      max_iter = tmp_sample.y;
    }
  }
#ifdef STORMM_USE_HIP
  const int init_shfl_stride = 32;
#else
  const int init_shfl_stride = 16;
#endif
  for (int i = init_shfl_stride; i > 0; i >>= 1) {
    const int max_comp = SHFL_DOWN(max_iter, i);
    if (max_comp > max_iter) {
      max_iter = max_comp;
    }
  }
  __syncthreads();
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (lane_idx == 0 && max_iter > 0) {
    atomicMax((int*)&smax_iter, max_iter);
  }
  __syncthreads();
  const int ginc = blockDim.x * gridDim.x;
  for (int pos = threadIdx.x + (blockIdx.x * blockDim.x); pos < n_generators; pos += ginc) {
    const ullint2 txy_state = state_xy[pos];
    const ullint2 tzw_state = state_zw[pos];
    ullint4 tmp_state = { txy_state.x, txy_state.y, tzw_state.x, tzw_state.y };
    for (int i = 0; i <= smax_iter; i++) {

      // Get a uniform random number in the range [0, 1).  Casting the unsigned long long int to
      // a (signed) long long int works as anything great than 2^3 just becomes negative.  The
      // bit string is unchanged.
      const ullint sxsw = tmp_state.x + tmp_state.w;
      const ullint rndbits = tmp_state.x + ((sxsw << 23) | (sxsw >> (64 - 23)));
      const llint work = (((rndbits >> 12) & 0xfffffffffffff) | 0x3ff0000000000000);
      const double rn_out = __longlong_as_double(work) - 1.0;

      // Horribly inefficient loop and access pattern, but this is just
      // a test program.  All of this is L1-cached, too.
      for (int j = 0; j < n_samples; j++) {
        const int2 tmp_sample = samples[j];
        if (pos == tmp_sample.x && i == tmp_sample.y) {
          random_output[j] = rn_out;
        }
      }

      // Push the state forward
      const ullint t = (tmp_state.y << 17);
      tmp_state.z ^= tmp_state.x;
      tmp_state.w ^= tmp_state.y;
      tmp_state.y ^= tmp_state.z;
      tmp_state.x ^= tmp_state.w;
      tmp_state.z ^= t;
      tmp_state.w = ((tmp_state.w << 45) | (tmp_state.w >> (64 - 45)));
    }

    // In a simulation, the state vector would be updated after computations
    // for this thread or atom are complete.
    state_xy[pos] = { tmp_state.x, tmp_state.y };
    state_zw[pos] = { tmp_state.z, tmp_state.w };
  }
}

//-------------------------------------------------------------------------------------------------
// Reproduce the random number generated by one of a series of xoroshiro128+ state vectors.
//
// Arguments:
//   rng_states:     List of xoroshiro128+ long-jump generator states (the dimension of this array
//                   implies the number of long-jumps taken when seeding the various states)
//   generator_idx:  Index of the generator state to query.  In a simulation, this might correspond
//                   to the atom index, or perhaps to the thread index within a particular launch
//                   grid.
//   iteration:      Produce the pseudo-random number for this point in the sequence of the
//                   particular atom or thread.
//-------------------------------------------------------------------------------------------------
double pinpointXoroshiro128p(std::vector<ullint2> &rng_states, const int generator_idx,
                             const uint iteration) {
  
  // Determine the initial generator state, possibly after having taken some short jumps.
  // The rng_states vector was created by taking up to max_xo_long_jumps long jumps of a single
  // state.  A given generator state is then created by tiling this series of initial, long jump
  // states with 1, 2, ..., n additional jumps, up to the total number of generators needed.  The
  // total number of generators is irrelevant.  This calculation requires the index for just one.
  const int n_seeds       = rng_states.size();
  const int n_short_jumps = generator_idx / n_seeds;
  const int seed_idx      = generator_idx - (n_short_jumps * n_seeds);

  // Recover the generator initial state of interest.  In a simulation, an initial generator X can
  // be jumped forward by P long jumps and N short jumps to arrive at a subsidiary generator X'.
  // Advancing X' by K iterations will produce the same output random number as advancing X by K
  // iterations, then jumping forward by P long jumps and N short jumps.
  Xoroshiro128pGenerator tgen(rng_states[seed_idx]);
  for (int i = 0; i < n_short_jumps; i++) {
    tgen.jump();
  }
  
  // Get the double-precision random number resulting from the requested iteration
  double result;
  for (uint i = 0; i <= iteration; i++) {
    result = tgen.uniformRandomNumber();
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Reproduce the random number generated by one of a series of xoshiro256++ state vectors.
//
// Arguments:
//   rng_xy_states:  List of first halves of xoshiro256++ long-jump generator states (the dimension
//                   of this array implies the number of long-jumps taken when seeding the various
//                   states)
//   rng_zw_states:  List of second halves of xoshiro256++ long-jump generator states
//   generator_idx:  Index of the generator state to query.  In a simulation, this might correspond
//                   to the atom index, or perhaps to the thread index within a particular launch
//                   grid.
//   iteration:      Produce the pseudo-random number for this point in the sequence of the
//                   particular atom or thread.
//-------------------------------------------------------------------------------------------------
double pinpointXoshiro256pp(std::vector<ullint2> &rng_xy_states,
                            std::vector<ullint2> &rng_zw_states, const int generator_idx,
                            const uint iteration) {
  
  // Determine the initial generator state, possibly after having taken some short jumps.
  // The rng_states vector was created by taking up to max_xo_long_jumps long jumps of a single
  // state.  A given generator state is then created by tiling this series of initial, long jump
  // states with 1, 2, ..., n additional jumps, up to the total number of generators needed.  The
  // total number of generators is irrelevant.  This calculation requires the index for just one.
  const int n_seeds       = rng_xy_states.size();
  const int n_short_jumps = generator_idx / n_seeds;
  const int seed_idx      = generator_idx - (n_short_jumps * n_seeds);

  // Recover the generator initial state of interest.  In a simulation, an initial generator X can
  // be jumped forward by P long jumps and N short jumps to arrive at a subsidiary generator X'.
  // Advancing X' by K iterations will produce the same output random number as advancing X by K
  // iterations, then jumping forward by P long jumps and N short jumps.
  const ullint4 seeded_state = { rng_xy_states[seed_idx].x, rng_xy_states[seed_idx].y,
                                 rng_zw_states[seed_idx].x, rng_zw_states[seed_idx].y };
  Xoshiro256ppGenerator tgen(seeded_state);
  for (int i = 0; i < n_short_jumps; i++) {
    tgen.jump();
  }
  
  // Get the double-precision random number resulting from the requested iteration
  double result;
  for (uint i = 0; i <= iteration; i++) {
    result = tgen.uniformRandomNumber();
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// This kernel will fill in some arrays with B-spline coefficients based on values of the delta.
// Derivatives can be computed by providing a negative value of the interpolation order.
//
// Arguments:
//   dx:              The array of deltas for computing B-splines
//   n:               The trusted length of dx
//   order:           The order of B-spline coefficients to compute.  Specifying negative values
//                    will have derivatives computed.
//   coefficients:    Array of coefficients, filled and returned, ordered in stretches of values
//                    for dx(0), dx(1), ..., dx(n)
//   derivatives:     Array of derivatives, filled and returned, ordered in stretches of values
//                    for dx(0), dx(1), ..., dx(n)
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kCalculateBSplines(const double* dx, const int n, const int order, double* coefficients,
                   double* derivatives) {
  double bspln_knots[6], bspln_dervs[6];
  int pos = threadIdx.x + (blockIdx.x * blockDim.x);
  while (pos < n) {
    if (order == 4) {
      devcBSpline4(dx[pos], bspln_knots);
    }
    else if (order == 5) {
      devcBSpline5(dx[pos], bspln_knots);
    }
    else if (order == 6) {
      devcBSpline6(dx[pos], bspln_knots);
    }
    else if (order == -4) {
      devcBSpline4(dx[pos], bspln_knots, bspln_dervs);
    }
    else if (order == -5) {
      devcBSpline5(dx[pos], bspln_knots, bspln_dervs);
    }
    else if (order == -6) {
      devcBSpline6(dx[pos], bspln_knots, bspln_dervs);
    }

    // Transfer the results to global memory
    const int abs_order = abs(order);
    for (int i = 0; i < abs_order; i++) {
      coefficients[(pos * abs_order) + i] = bspln_knots[i];
    }
    if (order < 0) {
      for (int i = 0; i < abs_order; i++) {
        derivatives[(pos * abs_order) + i] = bspln_dervs[i];
      }
    }
    pos += (blockDim.x * gridDim.x);
  }
}

//-------------------------------------------------------------------------------------------------
// Test various B-spline device functions.
//-------------------------------------------------------------------------------------------------
void testBSplineDeviceFuncs() {
  const int npts = 2048;
  Hybrid<double> dx(npts), coefficients(6 * npts), derivatives(6 * npts);
  for (int i = 0; i < npts; i++) {
    dx.putHost((static_cast<double>(i) + 0.28) / static_cast<double>(npts), i);
  }
  dx.upload();
  const HybridTargetLevel devc_layer = HybridTargetLevel::DEVICE;
  for (int ordr = 4; ordr < 7; ordr++) {

    // Compute B-spline knots only
    kCalculateBSplines<<<1, 1024>>>(dx.data(devc_layer), npts, ordr, coefficients.data(devc_layer),
                                    derivatives.data(devc_layer));
    coefficients.download();
    
    // Check the result against the CPU function
    std::vector<double> host_knots(ordr * npts);
    double* hkn_ptr = host_knots.data();
    for (int i = 0; i < npts; i++) {
      bSpline(dx.readHost(i), ordr, &hkn_ptr[i * ordr]);
    }
    check(coefficients.readHost(0, ordr * npts), RelationalOperator::EQUAL,
          Approx(host_knots).margin(1.0e-8), "B-spline knots of order " + std::to_string(ordr) +
          " were not computed correctly by the GPU device function.");
    
    // Compute B-spline knots and derivatives
    kCalculateBSplines<<<1, 1024>>>(dx.data(devc_layer), npts, -ordr,
                                    coefficients.data(devc_layer), derivatives.data(devc_layer));
    coefficients.download();
    derivatives.download();
    
    // Check the result against the CPU function
    std::vector<double> host_dervs(ordr * npts);
    double* hdv_ptr = host_dervs.data();
    for (int i = 0; i < npts; i++) {
      bSpline(dx.readHost(i), ordr, &hkn_ptr[i * ordr], &hdv_ptr[i * ordr]);
    }
    check(coefficients.readHost(0, ordr * npts), RelationalOperator::EQUAL,
          Approx(host_knots).margin(1.0e-8), "B-spline knots of order " + std::to_string(ordr) +
          " were not computed correctly by the GPU device function when derivatives are "
          "requested.");
    check(derivatives.readHost(0, ordr * npts), RelationalOperator::EQUAL,
          Approx(host_dervs).margin(1.0e-8), "B-spline derivatives of order " +
          std::to_string(ordr) + " were not computed correctly by the GPU device function when "
          "derivatives are requested.");
  }
}

//-------------------------------------------------------------------------------------------------
// CUDA kernels to perform split fixed precision addition & subtraction on different
// variable types. 
//
// Arguments:
//   arr_a:         Split Fixed Precision Array
//   arr_a_ovrf:    Split Fixed Precision Array Overflow for arr_a
//   arr_b:         Split Fixed Precision Array
//   arr_b_ovrf:    Split Fixed Precision Array Overflow for arr_b
//   results:       Array carrying results of an operation on arr_a and arr_b
//   results_ovrf:  Split Fixed Precision Array Overflow for results
//   cur_operator:  Enumerator detailing the current Math operator to apply
//   nums:          The number of total elements in each array
//   
//-------------------------------------------------------------------------------------------------
__global__ void testSplitFixedPrecision(int *arr_a, int *arr_a_ovrf, int *arr_b, int *arr_b_ovrf,
                                        int *results, int *results_ovrf,
                                        const Operator cur_operator, const FusionStyle method,
                                        const int nums) {
  int tid = threadIdx.x + (blockIdx.x * blockDim.x);
  while (tid < nums) {
    int2 cur_result;
    switch (cur_operator) {
    case Operator::SUBTRACT: 
      switch(method) {
      case FusionStyle::FUSED_FUSED: 
        {
          const int2 lh_number = { arr_a[tid], arr_a_ovrf[tid] };
          const int2 rh_number = { arr_b[tid], arr_b_ovrf[tid] };
          cur_result = splitFPSubtract(lh_number, rh_number);
        } 
        break;
      case FusionStyle::FUSED_DECPL:
        {
          const int2 lh_number = { arr_a[tid], arr_a_ovrf[tid] };
          cur_result = splitFPSubtract(lh_number, arr_b[tid], arr_b_ovrf[tid]);
        }
        break;
      case FusionStyle::DECPL_DECPL:
        cur_result = int63Subtract(arr_a[tid], arr_a_ovrf[tid], arr_b[tid], arr_b_ovrf[tid]);
        break;
      }
      break;

    case Operator::SUM:
      switch(method) {
      case FusionStyle::FUSED_FUSED: 
        {
          const int2 lh_number = { arr_a[tid], arr_a_ovrf[tid] };
          const int2 rh_number = { arr_b[tid], arr_b_ovrf[tid] };
          cur_result = splitFPSum(lh_number, rh_number);
        } 
        break;
      case FusionStyle::FUSED_DECPL:
        {
          const int2 lh_number = { arr_a[tid], arr_a_ovrf[tid] };
          cur_result = splitFPSum(lh_number, arr_b[tid], arr_b_ovrf[tid]);
        }
        break;
      case FusionStyle::FUSED_FLOAT:
        {
          const int2 lh_number = { arr_a[tid], arr_a_ovrf[tid] };
          const float rh_number = int63ToFloat(arr_b[tid], arr_b_ovrf[tid]);
          cur_result = splitFPSum(lh_number, rh_number);
        }
        break;
      case FusionStyle::DECPL_DECPL:
        cur_result = int63Sum(arr_a[tid], arr_a_ovrf[tid], arr_b[tid], arr_b_ovrf[tid]);
        break;
      case FusionStyle::DECPL_FLOAT:
        {
          const float rh_number = int63ToFloat(arr_b[tid], arr_b_ovrf[tid]);
          cur_result = int63Sum(arr_a[tid], arr_a_ovrf[tid], rh_number);
        }
        break;
      }
      break;
    }
    results[tid] = cur_result.x;
    results_ovrf[tid] = cur_result.y;
    tid += (blockDim.x * gridDim.x);
  }
}

__global__ void testSplitFixedPrecision(llint *arr_a, int *arr_a_ovrf,
                                        llint *arr_b, int *arr_b_ovrf, 
                                        llint *results, int *results_ovrf,
                                        const Operator cur_operator, const FusionStyle method,
                                        const int nums) {
  int tid = threadIdx.x + (blockIdx.x * blockDim.x);
  while (tid < nums) {
    int95_t cur_result;
    switch (cur_operator) {
    case Operator::SUBTRACT: 
      switch(method) {
      case FusionStyle::FUSED_FUSED: 
        {
          const int95_t lh_number = { arr_a[tid], arr_a_ovrf[tid] };
          const int95_t rh_number = { arr_b[tid], arr_b_ovrf[tid] };
          cur_result = splitFPSubtract(lh_number, rh_number);
        } 
        break;
      case FusionStyle::FUSED_DECPL:
        {
          const int95_t lh_number = { arr_a[tid], arr_a_ovrf[tid] };
          cur_result = splitFPSubtract(lh_number, arr_b[tid], arr_b_ovrf[tid]);
        }
        break;
      case FusionStyle::DECPL_DECPL:
        cur_result = int95Subtract(arr_a[tid], arr_a_ovrf[tid], arr_b[tid], arr_b_ovrf[tid]);
        break;
      }
      break;

    case Operator::SUM:
      switch(method) {
      case FusionStyle::FUSED_FUSED: 
        {
          const int95_t lh_number = { arr_a[tid], arr_a_ovrf[tid] };
          const int95_t rh_number = { arr_b[tid], arr_b_ovrf[tid] };
          cur_result = splitFPSum(lh_number, rh_number);
        } 
        break;
      case FusionStyle::FUSED_DECPL:
        {
          const int95_t lh_number = { arr_a[tid], arr_a_ovrf[tid] };
          cur_result = splitFPSum(lh_number, arr_b[tid], arr_b_ovrf[tid]);
        }
        break;
      case FusionStyle::FUSED_FLOAT:
        {
          const int95_t lh_number = { arr_a[tid], arr_a_ovrf[tid] };
          const double rh_number = int95ToDouble(arr_b[tid], arr_b_ovrf[tid]);
          cur_result = splitFPSum(lh_number, rh_number);
        }
        break;
      case FusionStyle::DECPL_DECPL:
        cur_result = int95Sum(arr_a[tid], arr_a_ovrf[tid], arr_b[tid], arr_b_ovrf[tid]);
        break;
      case FusionStyle::DECPL_FLOAT:
        {
          const double rh_number = int95ToDouble(arr_b[tid], arr_b_ovrf[tid]);
          cur_result = int95Sum(arr_a[tid], arr_a_ovrf[tid], rh_number);
        }
        break;
      }
      break;
    }
    results[tid] = cur_result.x;
    results_ovrf[tid] = cur_result.y;
    tid += (blockDim.x * gridDim.x);
  }
}

//-------------------------------------------------------------------------------------------------
// Perform multiplication on split fixed-precision numbers.  The results will bedeposited in a pair
// of holding arrays which can then be downloaded for inspection.
//
// Overloaded:
//   - Operate on int63_t (int2) factors
//   - Operate on int95_t factors
//
// Arguments:
//   primary:      Array of primary accumulators
//   overflow:     Array of overflow accumulators
//   multipliers:  Array of multiplier values
//   result_prim:  Array of primary accumulators for the results
//   result_ovrf:  Array of overflow accumulators for the results
//   ntest:        The number of tests to vectorize
//   fuse_inputs:  Flag to have input values fused into a tuple or remain separate (to engage
//                 distinct inline __device__ functions) 
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(small_block_size, 4)
testSFPMultiplication(const int* primary, const int* overflow, const int* multipliers,
                      int* results_prim, int* results_ovrf, const int ntest,
                      const bool fuse_inputs) {
  for (int i = threadIdx.x; i < ntest; i += blockDim.x * gridDim.x) {
    int2 r;
    if (fuse_inputs) {
      const int2 base = { primary[i], overflow[i] };
      r = splitFPMult(base, multipliers[i]);
    }
    else {
      r = int63Mult(primary[i], overflow[i], multipliers[i]);
    }
    results_prim[i] = r.x;
    results_ovrf[i] = r.y;
  }
}

__global__ void __launch_bounds__(small_block_size, 4)
testSFPMultiplication(const llint* primary, const int* overflow, const int* multipliers,
                      llint* results_prim, int* results_ovrf, const int ntest,
                      const bool fuse_inputs) {
  for (int i = threadIdx.x; i < ntest; i += blockDim.x * gridDim.x) {
    int95_t r;
    if (fuse_inputs) {
      const int95_t base = { primary[i], overflow[i] };
      r = splitFPMult(base, multipliers[i]);
    }
    else {
      r = int95Mult(primary[i], overflow[i], multipliers[i]);
    }
    results_prim[i] = r.x;
    results_ovrf[i] = r.y;
  }
}

//-------------------------------------------------------------------------------------------------
// Test split fixed-precision addition and subtraction on HPC resources in extreme value cases.
//
// Arguments:
//   nsmp:  The number of stremaing multiprocessors detected on the GPU
//   xrs:   Source of random numbers for additional tests
//-------------------------------------------------------------------------------------------------
void testSplitFixedPrecision(const int nsmp, Xoroshiro128pGenerator *xrs) {
  const int nums = 1024;
  const HybridTargetLevel devc_layer = HybridTargetLevel::DEVICE;
  Hybrid<int> arr_a(nums), arr_a_ovrf(nums), arr_b(nums), arr_b_ovrf(nums);
  Hybrid<int> results(nums), results_ovrf(nums);
  Hybrid<llint> llarr_a(nums), llarr_b(nums), llresults(nums);
  int* arr_a_pointer = arr_a.data();
  int* arr_a_ovrf_pointer = arr_a_ovrf.data();
  int* arr_b_pointer = arr_b.data();
  int* arr_b_ovrf_pointer = arr_b_ovrf.data();
  llint* llarr_a_pointer = llarr_a.data();
  llint* llarr_b_pointer = llarr_b.data();
  std::vector<double> rdc_a(nums), rdc_b(nums);
  std::vector<double> llrdc_a(nums), llrdc_b(nums);
  std::vector<double> target_ab_subtract(nums), target_ab_sum(nums);
  std::vector<double> lltarget_ab_subtract(nums), lltarget_ab_sum(nums);

  // Stage testing data and targets.
  const double scale_down = pow(2.0, -31.0);
  const double scale_down_ll = pow(2, -63);
  int ijkm = 0;
  for (int i = -2; i <= 2; ++i) {
    for (int j = -2; j <= 2; ++j) {
      for (int k = -2; k <= 2; ++k) {
        for (int m = -2; m <= 2; ++m) {
          if (i < 0) {
           arr_a_pointer[ijkm] = INT_MAX - (-1 - i);
           arr_b_pointer[ijkm] = INT_MAX - (-1 - k);
           llarr_a_pointer[ijkm] = LLONG_MAX - (-1 - i);
           llarr_b_pointer[ijkm] = LLONG_MAX - (-1 - k);
          }
          else {
           arr_a_pointer[ijkm] = INT_MIN + i;
           arr_b_pointer[ijkm] = INT_MIN + k;
           llarr_a_pointer[ijkm] = LLONG_MIN + i;
           llarr_b_pointer[ijkm] = LLONG_MIN + k;
          }
          arr_a_ovrf_pointer[ijkm] = j;
          arr_b_ovrf_pointer[ijkm] = m;
          rdc_a[ijkm] = hostInt63ToDouble(arr_a_pointer[ijkm], arr_a_ovrf_pointer[ijkm]);
          rdc_b[ijkm] = hostInt63ToDouble(arr_b_pointer[ijkm], arr_b_ovrf_pointer[ijkm]);
          llrdc_a[ijkm] = hostInt95ToDouble(llarr_a_pointer[ijkm], arr_a_ovrf_pointer[ijkm]);
          llrdc_b[ijkm] = hostInt95ToDouble(llarr_b_pointer[ijkm], arr_b_ovrf_pointer[ijkm]);
          target_ab_subtract[ijkm] = (rdc_a[ijkm] - rdc_b[ijkm]) * scale_down;
          target_ab_sum[ijkm] = (rdc_a[ijkm] + rdc_b[ijkm]) * scale_down;
          lltarget_ab_subtract[ijkm] = (llrdc_a[ijkm] - llrdc_b[ijkm]) * scale_down_ll;
          lltarget_ab_sum[ijkm] = (llrdc_a[ijkm] + llrdc_b[ijkm]) * scale_down_ll;
          ijkm++;
        }
      }
    }
  } 
  arr_a.upload();
  arr_a_ovrf.upload();
  arr_b.upload();
  arr_b_ovrf.upload();
  llarr_a.upload();
  llarr_b.upload();
  std::vector<Operator> operators = { Operator::SUBTRACT, Operator::SUM };
  std::vector<FusionStyle> fusions = { FusionStyle::FUSED_FUSED, FusionStyle::FUSED_DECPL,
                                       FusionStyle::FUSED_FLOAT, FusionStyle::DECPL_DECPL,
                                       FusionStyle::DECPL_FLOAT };
  for (Operator& cur_operator : operators) {
    for (FusionStyle& cur_fusion : fusions) {
      if (cur_operator == Operator::SUBTRACT &&
          (cur_fusion == FusionStyle::FUSED_FLOAT || cur_fusion == FusionStyle::DECPL_FLOAT)) {
        continue;
      }
      testSplitFixedPrecision<<<4 * nsmp, small_block_size>>>(arr_a.data(devc_layer),
                                                              arr_a_ovrf.data(devc_layer), 
                                                              arr_b.data(devc_layer),
                                                              arr_b_ovrf.data(devc_layer), 
                                                              results.data(devc_layer),
                                                              results_ovrf.data(devc_layer),
                                                              cur_operator, cur_fusion, nums);
      results.download();
      results_ovrf.download();
      std::vector<double> rdc_gpu(nums);
      for (int i = 0; i < nums; i++) {
        rdc_gpu[i] = hostInt63ToDouble(results.readHost(i), results_ovrf.readHost(i));
        rdc_gpu[i] *= scale_down;
      }
      switch(cur_operator) {
      case Operator::SUBTRACT:
        check(rdc_gpu, RelationalOperator::EQUAL, target_ab_subtract, 
              "Operation involving int63 fails. Operator: " + 
              getOperatorName(cur_operator) + " | FusionStyle: " +
              getFusionStyleName(cur_fusion) + "\n");
        break;
      case Operator::SUM:
        check(rdc_gpu, RelationalOperator::EQUAL, target_ab_sum, 
              "Operation involving int63 fails. Operator: " + 
              getOperatorName(cur_operator) + " | FusionStyle: " +
              getFusionStyleName(cur_fusion) + "\n");

        break;
      }
      testSplitFixedPrecision<<<4 * nsmp, small_block_size>>>(llarr_a.data(devc_layer),
                                                              arr_a_ovrf.data(devc_layer), 
                                                              llarr_b.data(devc_layer),
                                                              arr_b_ovrf.data(devc_layer), 
                                                              llresults.data(devc_layer),
                                                              results_ovrf.data(devc_layer),
                                                              cur_operator, cur_fusion, nums);
      llresults.download();
      results_ovrf.download();
      for (int i = 0; i < nums; i++) {
        rdc_gpu[i] = hostInt95ToDouble(llresults.readHost(i), results_ovrf.readHost(i));
        rdc_gpu[i] *= scale_down_ll;
      }
      switch (cur_operator) {
      case Operator::SUBTRACT:
        check(rdc_gpu, RelationalOperator::EQUAL, lltarget_ab_subtract, 
              "Operation involving int95_t fails. Operator: " + 
              getOperatorName(cur_operator) + " | FusionStyle: " +
              getFusionStyleName(cur_fusion) + "\n");
        break;
      case Operator::SUM:
        check(rdc_gpu, RelationalOperator::EQUAL, lltarget_ab_sum, 
              "Operation involving int95_t fails. Operator: " + 
              getOperatorName(cur_operator) + " | FusionStyle: " +
              getFusionStyleName(cur_fusion) + "\n");
        break;
      }
    }
  }

  // Test multiplication operations
  const int n_mult_test = 216;
  Hybrid<int> primary(n_mult_test), overflow(n_mult_test), resmul_prim(n_mult_test);
  Hybrid<llint> primary_ll(n_mult_test), resmul_ll_prim(n_mult_test);
  Hybrid<int> multipliers(n_mult_test), resmul_ovrf(n_mult_test);
  int* primary_ptr = primary.data();
  llint* primary_ll_ptr = primary_ll.data();
  int* overflow_ptr = overflow.data();
  int* multipliers_ptr = multipliers.data();
  for (int i = 0; i < 5; i++) {
    const int ext_val = (i < 2) ? INT_MAX - i : INT_MIN + (i - 2);
    const llint ext_ll_val = (i < 2) ? LLONG_MAX - static_cast<llint>(i) :
                                       LLONG_MIN + static_cast<llint>(i - 2);
    for (int j = 0; j < 5; j++) {
      const int ovrf_val = (j < 2) ? INT_MAX - j : INT_MIN + (j - 2);
      for (int k = 0; k < 5; k++) {
        const size_t test_idx = (((k * 6) + j) * 6) + i;
        const int mult_val = (j < 2) ? INT_MAX - j : INT_MIN + (j - 2);
        primary_ptr[test_idx] = ext_val;
        primary_ll_ptr[test_idx] = ext_ll_val;
        overflow_ptr[test_idx] = ovrf_val;
        multipliers_ptr[test_idx] = mult_val;
      }
    }
  }
  for (int i = 0; i < 6; i++) {

    // If the random numbers take these values far to the extremes, there may be truncation, but
    // the test only needs values scattered across the number line. 
    const int prim_val = (xrs->uniformRandomNumber() - 0.5) * static_cast<double>(UINT_MAX);
    const llint prim_ll_val = (xrs->uniformRandomNumber() - 0.5) * static_cast<double>(ULLONG_MAX);
    for (int j = 0; j < 6; j++) {
      const int ovrf_val = (xrs->uniformRandomNumber() - 0.5) * static_cast<double>(256.0);
      for (int k = 0; k < 6; k++) {
        if (i < 5 && j < 5 && k < 5) {
          continue;
        }
        const int mult_val = (xrs->uniformRandomNumber() - 0.5) * 256.0;
        const size_t test_idx =	(((k * 6) + j) * 6) + i;
        primary_ptr[test_idx] = prim_val;
        primary_ll_ptr[test_idx] = prim_ll_val;
        overflow_ptr[test_idx] = ovrf_val;
        multipliers_ptr[test_idx] = mult_val;
      }
    }
  }
  std::vector<double> chk_multres(n_mult_test), chk_ll_multres(n_mult_test);
  std::vector<double> multres(n_mult_test), ll_multres(n_mult_test);
  for (int i = 0; i < n_mult_test; i++) {
    chk_multres[i] = hostInt63ToDouble(hostInt63Mult(primary_ptr[i], overflow_ptr[i],
                                                     multipliers_ptr[i])) * scale_down;
    chk_ll_multres[i] = hostInt95ToDouble(hostInt95Mult(primary_ll_ptr[i], overflow_ptr[i],
                                                        multipliers_ptr[i])) * scale_down_ll;
  }
  primary.upload();
  primary_ll.upload();
  overflow.upload();
  multipliers.upload();
  const std::vector<bool> fuse_or_not = { true, false };
  const int* res_ptr = resmul_prim.data();
  const llint* res_ll_ptr = resmul_ll_prim.data();
  const int* res_ovrf_ptr = resmul_ovrf.data();
  for (size_t i = 0; i < 2; i++) {
    testSFPMultiplication<<<4 * nsmp, small_block_size>>>(primary.data(devc_layer),
                                                          overflow.data(devc_layer),
                                                          multipliers.data(devc_layer),
                                                          resmul_prim.data(devc_layer),
                                                          resmul_ovrf.data(devc_layer),
                                                          n_mult_test, fuse_or_not[i]);
    resmul_prim.download();
    resmul_ovrf.download();
    for (int j = 0; j < n_mult_test; j++) {
      multres[j] = hostInt63ToDouble(res_ptr[j], res_ovrf_ptr[j]) * scale_down;
    }
    testSFPMultiplication<<<4 * nsmp, small_block_size>>>(primary_ll.data(devc_layer),
                                                          overflow.data(devc_layer),
                                                          multipliers.data(devc_layer),
                                                          resmul_ll_prim.data(devc_layer),
                                                          resmul_ovrf.data(devc_layer),
                                                          n_mult_test, fuse_or_not[i]);
    resmul_ll_prim.download();
    resmul_ovrf.download();
    for (int j = 0; j < n_mult_test; j++) {
      multres[j] = hostInt95ToDouble(res_ll_ptr[j], res_ovrf_ptr[j]) * scale_down_ll;
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Kernel to encapsulate the block-wide prefix sum.  The scratch arrays are held in __shared__
// memory while the main prefix sum is computed in __global__ memory.
//-------------------------------------------------------------------------------------------------
template <typename T>
__global__ void __launch_bounds__(large_block_size, 1) kDoBlockPrefix(T* p, const int ntest) {
  __shared__ T s1[warp_size_int * warp_size_int], s2[warp_size_int];
  if (blockIdx.x == 0) {
    blockExclusivePrefixSum(p, s1, s2, ntest);
  }
}

//-------------------------------------------------------------------------------------------------
// Test the block-wide prefix sum.
//
// Arguments:
//   ntest:  The quantity of numbers in the prefix sum series.
//   xrs:    Source of random numbers for prefix sums
//-------------------------------------------------------------------------------------------------
template <typename T>
void testBlockPrefixSum(const int ntest, Xoroshiro128pGenerator *xrs) {
  Hybrid<T> p(ntest + 1, "prefix_array");
  std::vector<T> pv(ntest + 1, static_cast<T>(0));
  if (isFloatingPointScalarType<T>()) {
    for (int i = 0; i < ntest; i++) {
      pv[i] = xrs->uniformRandomNumber() * 50.0;
    }
  }
  else {
    for (int i = 0; i < ntest; i++) {
      pv[i] = round(xrs->uniformRandomNumber() * 50.0);
    }
  }
  p.putHost(pv, 0, ntest + 1);
  p.upload();
  kDoBlockPrefix<T><<<1, large_block_size>>>(p.data(HybridTargetLevel::DEVICE), ntest + 1);
  prefixSumInPlace<T>(&pv, PrefixSumType::EXCLUSIVE);
  const std::vector<T> pd = p.readDevice();
  check(pd, RelationalOperator::EQUAL, Approx(pd).margin(1.0e-4), "A prefix sum of " +
        getStormmScalarTypeName<T>() + " computed by the parallel GPU method did not match that "
        "computed by the serial CPU method.");
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  HpcConfig gpu_config(ExceptionResponse::WARN);
  std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  
  // Section 1
  section("Vector processing capabilities");
  
  // Section 2
  section("GPU-based Xoroshiro128+ PRNG");

  // Section 3
  section("Included device functions");
  
  // Perform a summation over a double-precision real vector using the GPU
  section(1);
  const int n_tiny  = 16;
  const int n_small = 128;
  const int n_chunk = 517;
  const int n_block = 1024;
  const int n_large = 23552;
  const int n_giant = 2500000;
  Hybrid<double> tiny_set(n_tiny, "tiny_vector");
  Hybrid<double> small_set(n_small, "small_vector");
  Hybrid<double> chunk_set(n_chunk, "chunk_vector");
  Hybrid<double> block_set(n_block, "block_vector");
  Hybrid<double> large_set(n_large, "large_vector");
  Hybrid<double> giant_set(n_giant, "giant_vector");
  Hybrid<double> tb_buffer(gpu.getSMPCount(), "sum_accumulators", HybridFormat::HOST_MOUNTED);
  Ran2Generator prng(oe.getRandomSeed());
  loadRan2ByCPU(&tiny_set, &prng);
  loadRan2ByCPU(&small_set, &prng);
  loadRan2ByCPU(&chunk_set, &prng);
  loadRan2ByCPU(&block_set, &prng);
  tiny_set.upload();
  loadRan2ByCPU(&tiny_set, &prng);
  const double gpu_tiny_sum  = sum(tiny_set, &tb_buffer, gpu);
  const double cpu_tiny_sum  = sum(tiny_set, &tb_buffer, gpu, HybridTargetLevel::HOST);
  check(gpu_tiny_sum, RelationalOperator::EQUAL, Approx(-0.4835516929).margin(1.0e-8),
        "The tiniest vector did not sum correctly on the GPU.\n");
  check(cpu_tiny_sum, RelationalOperator::EQUAL, Approx(-1.1269689474).margin(1.0e-8),
        "The tiniest vector did not sum correctly on the CPU.\n");
  small_set.upload();
  chunk_set.upload();
  block_set.upload();
  check(sum(small_set, &tb_buffer, gpu), RelationalOperator::EQUAL, sum<double>(small_set),
        "The small vector did not sum correctly on the GPU.\n");
  check(sum(chunk_set, &tb_buffer, gpu), RelationalOperator::EQUAL, sum<double>(chunk_set),
        "The medium, odd-sized vector did not sum correctly on the GPU.\n");
  check(sum(block_set, &tb_buffer, gpu), RelationalOperator::EQUAL, sum<double>(block_set),
        "The full block-sized vector did not sum correctly on the GPU.\n");
  Xoroshiro128pGenerator fast_prng(78172);
  loadXoroshiro128pByCPU(&large_set, &fast_prng);
  loadXoroshiro128pByCPU(&giant_set, &fast_prng);
  large_set.upload();
  giant_set.upload();
  check(sum(large_set, &tb_buffer, gpu), RelationalOperator::EQUAL, sum<double>(large_set),
        "The large vector did not sum correctly on the GPU.\n");
  check(sum(giant_set, &tb_buffer, gpu), RelationalOperator::EQUAL, sum<double>(giant_set),
        "The giant vector did not sum correctly on the GPU.\n");
  
  // Test the GPU-base random number seeding and synchronized CPU/GPU generation
  const int n_cellulose_atoms = 408609;
  Hybrid<ullint2> rng128p_states(n_cellulose_atoms, "xoroshiro128p_state");
  Hybrid<ullint2> rng256pp_xy_states(n_cellulose_atoms, "xoroshiro256pp_sxy");
  Hybrid<ullint2> rng256pp_zw_states(n_cellulose_atoms, "xoroshiro256pp_szw");
  initXoroshiro128pArray(&rng128p_states, 8773925, default_xoroshiro128p_scrub, gpu);
  initXoshiro256ppArray(&rng256pp_xy_states, &rng256pp_zw_states, 4091832,
                        default_xoshiro256pp_scrub, gpu);
  std::vector<ullint2>cpu_rng128p_states(n_cellulose_atoms);
  std::vector<ullint2> cpu_rng256pp_xy_states(n_cellulose_atoms);
  std::vector<ullint2> cpu_rng256pp_zw_states(n_cellulose_atoms);
  initXoroshiro128pArray(&cpu_rng128p_states, 8773925, default_xoroshiro128p_scrub);
  initXoshiro256ppArray(&cpu_rng256pp_xy_states, &cpu_rng256pp_zw_states, 4091832,
                        default_xoshiro256pp_scrub);
  int xrs128p_deviations = 0;
  int xrs256pp_deviations = 0;
  const ullint2* rng128p_st_ptr  = rng128p_states.data();
  const ullint2* rng256pp_st_xy_ptr = rng256pp_xy_states.data();
  const ullint2* rng256pp_st_zw_ptr = rng256pp_zw_states.data();
  for (int i = 0; i < n_cellulose_atoms; i++) {
    xrs128p_deviations  += (cpu_rng128p_states[i].x != rng128p_st_ptr[i].x ||
                            cpu_rng128p_states[i].y != rng128p_st_ptr[i].y);
    xrs256pp_deviations += (cpu_rng256pp_xy_states[i].x != rng256pp_st_xy_ptr[i].x ||
                            cpu_rng256pp_xy_states[i].y != rng256pp_st_xy_ptr[i].y ||
                            cpu_rng256pp_zw_states[i].x != rng256pp_st_zw_ptr[i].x ||
                            cpu_rng256pp_zw_states[i].y != rng256pp_st_zw_ptr[i].y);
  }
  check(xrs128p_deviations, RelationalOperator::EQUAL, 0, "Deviations were found between "
        "CPU-initialized and GPU-initialized Xoroshiro128+ random number generator arrays.");
  check(xrs256pp_deviations, RelationalOperator::EQUAL, 0, "Deviations were found between "
        "CPU-initialized and GPU-initialized Xoshiro256++ random number generator arrays.");
  rng128p_states.download();
  rng256pp_xy_states.download();
  rng256pp_zw_states.download();
  Xoroshiro128pGenerator xrs128p_check(8773925);
  Xoshiro256ppGenerator xrs256pp_check(4091832);
  const int n_seeds_made = std::min(max_xo_long_jumps, n_cellulose_atoms);
  std::vector<ullint2> cpu_128p_seeds(n_seeds_made);
  std::vector<ullint2> cpu_256pp_xy_seeds(n_seeds_made);
  std::vector<ullint2> cpu_256pp_zw_seeds(n_seeds_made);
  for (int i = 0; i < n_seeds_made; i++) {
    cpu_128p_seeds[i] = xrs128p_check.revealState();
    const ullint4 cpu_found_state = xrs256pp_check.revealState();
    cpu_256pp_xy_seeds[i] = { cpu_found_state.x, cpu_found_state.y };
    cpu_256pp_zw_seeds[i] = { cpu_found_state.z, cpu_found_state.w };
    xrs128p_check.longJump();
    xrs256pp_check.longJump();
  }

  // Create a smattering of generator indices (within the bounds of the rng128p_states above) and
  // some (low) iteration counts at which to test each of them.  Predict the results on the CPU,
  // then compute them on the GPU.
  const int n_samples = 16;
  Hybrid<int2> samples(n_samples, "assorted_points");
  int2* samp_ptr = samples.data();
  for (int i = 0; i < n_samples; i++) {
    samp_ptr[i] = { 918 * i, (2 * i) + 5 };
  }
  Hybrid<double> random_output(n_samples, "random_pluckings");
  double* ro_ptr = random_output.data();
  for (int i = 0; i < n_samples; i++) {
    ro_ptr[i] = pinpointXoroshiro128p(cpu_128p_seeds, samp_ptr[i].x, samp_ptr[i].y);
  }
  samples.upload();
  const int nsmp = gpu.getSMPCount();
  const int nthr = gpu.getMaxThreadsPerBlock();
  kEvalXoroshiro128p<<<nsmp, nthr>>>(rng128p_states.data(HybridTargetLevel::DEVICE),
                                     n_cellulose_atoms, samples.data(HybridTargetLevel::DEVICE),
                                     random_output.data(HybridTargetLevel::DEVICE), n_samples);
  check(random_output.readHost(), RelationalOperator::EQUAL, random_output.readDevice(), "Random "
        "numbers from an array of Xoroshiro128+ generators computed on the CPU and GPU do not "
        "agree.");
  kEvalXoshiro256pp<<<nsmp, nthr>>>(rng256pp_xy_states.data(HybridTargetLevel::DEVICE),
                                    rng256pp_zw_states.data(HybridTargetLevel::DEVICE),
                                    n_cellulose_atoms, samples.data(HybridTargetLevel::DEVICE),
                                    random_output.data(HybridTargetLevel::DEVICE), n_samples);
  for (int i = 0; i < n_samples; i++) {
    ro_ptr[i] = pinpointXoshiro256pp(cpu_256pp_xy_seeds, cpu_256pp_zw_seeds, samp_ptr[i].x,
                                     samp_ptr[i].y);
  }
  check(random_output.readHost(), RelationalOperator::EQUAL, random_output.readDevice(), "Random "
        "numbers from an array of Xoshiro256++ generators computed on the CPU and GPU do not "
        "agree.");

  // Test other device operations
  section(3);
  testBSplineDeviceFuncs();
  testSplitFixedPrecision(nsmp, &xrs128p_check);
  const std::vector<int> prfx_lengths = { 29, 30, 31, 32, 33, 34, 62, 63, 64, 65, 66, 1022, 1023,
                                          1024, 1025, 1026, 2046, 2047, 2048, 2049, 2050 };
  for (size_t i = 0; i < prfx_lengths.size(); i++) {
    testBlockPrefixSum<int>(prfx_lengths[i], &fast_prng);
    testBlockPrefixSum<double>(prfx_lengths[i], &fast_prng);
  }
  
  // Print results
  printTestSummary(oe.getVerbosity());
  return countGlobalTestFailures();
}
