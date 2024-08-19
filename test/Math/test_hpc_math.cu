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
using stormm::data_types::llint;
using stormm::data_types::int95_t;
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
// Enumerator for testing Split Fixed Precision Math on Device
//
//   Operator:        The operator to be used (summation or subtraction)
//-------------------------------------------------------------------------------------------------
enum class Operator {
  SUM,          ///< Split Fixed Integer Summation
  SUBTRACT      ///< Split Fixed Integer Subtraction
};

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
                                        const Operator cur_operator, int nums) {
  int tid = threadIdx.x + (blockIdx.x * blockDim.x);
  while (tid < nums) {
    switch (cur_operator) {
    case Operator::SUBTRACT: 
      {
        int2 cur_result = int63Subtract(arr_a[tid], arr_a_ovrf[tid], arr_b[tid], arr_b_ovrf[tid]);
        results[tid] = cur_result.x;
        results_ovrf[tid] = cur_result.y;
      }
      break;
    case Operator::SUM:
      {
        int2 cur_result = int63Sum(arr_a[tid], arr_a_ovrf[tid], arr_b[tid], arr_b_ovrf[tid]);
        results[tid] = cur_result.x;
        results_ovrf[tid] = cur_result.y;
      }
      break;
    }
    tid += (blockDim.x * gridDim.x);
  }
}

__global__ void testSplitFixedPrecision(llint *arr_a, int *arr_a_ovrf,
                                        llint *arr_b, int *arr_b_ovrf, 
                                        llint *results, int *results_ovrf,
                                        Operator cur_operator, int nums) {
  int tid = threadIdx.x + (blockIdx.x * blockDim.x);
  while (tid < nums) {
    switch (cur_operator) {
    case Operator::SUBTRACT:
      {
        int95_t cur_result = int95Subtract( arr_a[tid], arr_a_ovrf[tid],
                                            arr_b[tid], arr_b_ovrf[tid]);
        results[tid] = cur_result.x;
        results_ovrf[tid] = cur_result.y;
      }
      break;
    case Operator::SUM:
      {
        int95_t cur_result = int95Sum( arr_a[tid], arr_a_ovrf[tid],
                                       arr_b[tid], arr_b_ovrf[tid]);
        results[tid] = cur_result.x;
        results_ovrf[tid] = cur_result.y;
      }
      break;
    }
    tid += (blockDim.x * gridDim.x);
  }
}

template<typename T, typename T_x, typename T_y>
__global__ void testSplitFixedPrecision(T *arr_a, T *arr_b, 
                                        T_x *results, T_y *results_ovrf, Operator cur_operator,
                                        int nums) {
  int tid = threadIdx.x + (blockIdx.x * blockDim.x);
  while (tid < nums) {
    switch (cur_operator) {
    case Operator::SUBTRACT:
      {
        T cur_result = splitFPSubtract(arr_a[tid], arr_b[tid]);
        results[tid] = cur_result.x;
        results_ovrf[tid] = cur_result.y;
      }  
      break;
    case Operator::SUM:
      {
        T cur_result = splitFPSum(arr_a[tid], arr_b[tid]);
        results[tid] = cur_result.x;
        results_ovrf[tid] = cur_result.y;
      }
      break;
    }
    tid += (blockDim.x * gridDim.x);
  }
}

template<typename T, typename T_x, typename T_y>
__global__ void testSplitFixedPrecision(T *arr_a, T_x *arr_b, T_y *arr_b_ovrf, 
                                        T_x *results, T_y *results_ovrf, Operator cur_operator,
                                        int nums) {
  int tid = threadIdx.x + (blockIdx.x * blockDim.x);
  while (tid < nums) {
    switch (cur_operator) {
    case Operator::SUBTRACT:
      {
        T cur_result = splitFPSubtract(arr_a[tid], arr_b[tid], arr_b_ovrf[tid]);
        results[tid] = cur_result.x;
        results_ovrf[tid] = cur_result.y;
      }
      break;
    case Operator::SUM:
      {
        T cur_result = splitFPSum(arr_a[tid], arr_b[tid], arr_b_ovrf[tid]);
        results[tid] = cur_result.x;
        results_ovrf[tid] = cur_result.y;
      }
      break;
    }
    tid += (blockDim.x + gridDim.x);
  }
}

//-------------------------------------------------------------------------------------------------
// Test split fixed-precision addition and subtraction on HPC resources in extreme value cases.
//
// Arguments:
//   nsmp:  The number of stremaing multiprocessors detected on the GPU
//-------------------------------------------------------------------------------------------------
void testSplitFixedPrecision(const int nsmp) {
  const int nums = 1024;
  const HybridTargetLevel devc_layer = HybridTargetLevel::DEVICE;
  Hybrid<int> arr_a(nums), arr_a_ovrf(nums), arr_b(nums), arr_b_ovrf(nums);
  Hybrid<int2> arr_a_tuple, arr_b_tuple;
  Hybrid<int> results(nums), results_ovrf(nums);

  // Test 1: Basic tests on Float -> Int63
  for(int i = 0; i < nums; ++i) {
    hostFloatToInt63(INT_MIN, &arr_a.data()[i], &arr_a_ovrf.data()[i]);
    hostFloatToInt63(0, &arr_b.data()[i], &arr_b_ovrf.data()[i]);
  } 
  arr_a.upload();
  arr_a_ovrf.upload();
  arr_b.upload();
  arr_b_ovrf.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a.data(devc_layer),
                                                      arr_a_ovrf.data(devc_layer), 
                                                      arr_b.data(devc_layer),
                                                      arr_b_ovrf.data(devc_layer), 
                                                      results.data(devc_layer),
                                                      results_ovrf.data(devc_layer),
                                                      Operator::SUBTRACT, nums);

  // Download results from device
  results.download();
  results_ovrf.download();
  std::vector<float> fl_result, fl_expected;
  for(int i = 0; i < nums; i++){
    fl_result.push_back(hostInt63ToFloat(results.readHost(i), results_ovrf.readHost(i))); 
    fl_expected.push_back(static_cast<float>(INT_MIN));
  }

  check(fl_result, RelationalOperator::EQUAL, fl_expected, "Float base case (subtraction with 0) "
        "was not calculated properly");

  // Test 2: Tests on Float 63 Underflows and Overflows
  // Test 2.1: Testing Float Underflow
  for(int i = 0; i < nums; ++i) {
    arr_a_tuple.pushBack(hostFloatToInt63(INT_MIN));
    hostFloatToInt63(1, &arr_b.data()[i], &arr_b_ovrf.data()[i]);
  }
  arr_a_tuple.upload();
  arr_b.upload();
  arr_b_ovrf.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a_tuple.data(devc_layer),
                                                      arr_b.data(devc_layer),
                                                      arr_b_ovrf.data(devc_layer),
                                                      results.data(devc_layer),
                                                      results_ovrf.data(devc_layer),
                                                      Operator::SUBTRACT, nums);

  // Download results from device
  results.download();
  results_ovrf.download();
  fl_result.clear();
  fl_expected.clear();
  for(int i = 0; i < nums; i++){
    fl_result.push_back(hostInt63ToFloat(results.readHost(i), results_ovrf.readHost(i)));
    fl_expected.push_back(static_cast<float>(INT_MIN) - static_cast<float>(1));
  }

  check(fl_result, RelationalOperator::EQUAL, fl_expected, "Error in calculating Float -> Int64 "
        "Underflows.");

  // Test 2.1.1: Testing yet another underflow, with the ovrf populated
  for(int i = 0; i < nums; ++i) {
    arr_b_tuple.pushBack(hostFloatToInt63(1));
  }
  arr_a_tuple.upload();
  arr_b_tuple.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a_tuple.data(devc_layer),
                                                      arr_b_tuple.data(devc_layer),
                                                      results.data(devc_layer),
                                                      results_ovrf.data(devc_layer),
                                                      Operator::SUBTRACT, nums);
  
  // Download results from device
  results.download();
  results_ovrf.download();
  fl_result.clear();
  fl_expected.clear();
  for(int i = 0; i < nums; i++){
    fl_result.push_back(hostInt63ToFloat(results.readHost(i), results_ovrf.readHost(i)));
    fl_expected.push_back(static_cast<float>(INT_MIN) - static_cast<float>(1));
  }
  check(fl_result, RelationalOperator::EQUAL, fl_expected, "Error in calculating Float -> Int64 "
        "Underflows (second overload).");
  
  // Test 2.2: Testing Float Overflow 
  for(int i = 0; i < nums; ++i) {
    hostFloatToInt63(INT_MAX, &arr_a.data()[i], &arr_a_ovrf.data()[i]);
    hostFloatToInt63(1, &arr_b.data()[i], &arr_b_ovrf.data()[i]);
  }
  arr_a.upload();
  arr_a_ovrf.upload();
  arr_b.upload();
  arr_b_ovrf.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a.data(devc_layer),
                                                      arr_a_ovrf.data(devc_layer), 
                                                      arr_b.data(devc_layer),
                                                      arr_b_ovrf.data(devc_layer),
                                                      results.data(devc_layer),
                                                      results_ovrf.data(devc_layer), 
                                                      Operator::SUM, nums);

  // Download results from device
  results.download();
  results_ovrf.download();
  fl_result.clear();
  fl_expected.clear();
  for(int i = 0; i < nums; i++){
    fl_result.push_back(hostInt63ToFloat(results.readHost(i), results_ovrf.readHost(i)));
    fl_expected.push_back(static_cast<float>(INT_MAX) + static_cast<float>(1));
  }

  check(fl_result, RelationalOperator::EQUAL, fl_expected, "Error in calculating Float -> Int64 "
        "Overflows.");

  // Test 2.2.1: Testing yet another overflow, with the ovrf populated
  for(int i = 0; i < nums; ++i) {
    arr_a_tuple.data()[i] = (hostFloatToInt63(INT_MAX));
  }
  arr_a_tuple.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a_tuple.data(devc_layer), 
                                                      arr_b_tuple.data(devc_layer),
                                                      results.data(devc_layer),
                                                      results_ovrf.data(devc_layer),
                                                      Operator::SUM, nums);
  
  // Download results from device
  results.download();
  results_ovrf.download();
  fl_result.clear();
  fl_expected.clear();
  for(int i = 0; i < nums; i++){
    fl_result.push_back(hostInt63ToFloat(results.readHost(i), results_ovrf.readHost(i)));
    fl_expected.push_back(static_cast<float>(INT_MAX) + static_cast<float>(1));
  }
  check(fl_result, RelationalOperator::EQUAL, fl_expected, "Error in calculating Float -> Int64 "
        "Overflows (second overload).");

  // Test 3: Basic tests on Llint -> Int63
  for(int i = 0; i < nums; ++i) {
    hostLongLongToInt63(INT_MIN, &arr_a.data()[i], &arr_a_ovrf.data()[i]);
    hostLongLongToInt63(0, &arr_b.data()[i], &arr_b_ovrf.data()[i]);
  } 
  arr_a.upload();
  arr_a_ovrf.upload();
  arr_b.upload();
  arr_b_ovrf.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a.data(devc_layer),
                                                      arr_a_ovrf.data(devc_layer), 
                                                      arr_b.data(devc_layer),
                                                      arr_b_ovrf.data(devc_layer),
                                                      results.data(devc_layer),
                                                      results_ovrf.data(devc_layer),
                                                      Operator::SUBTRACT, nums);

  // Download results from device
  results.download();
  results_ovrf.download();
  std::vector<llint> llint_result, llint_expected;
  for(int i = 0; i < nums; i++){
    llint_result.push_back(hostInt63ToLongLong(results.readHost(i), results_ovrf.readHost(i)));
    llint_expected.push_back(static_cast<llint>(INT_MIN));
  }

  check(llint_result, RelationalOperator::EQUAL, llint_expected, "Llint Base Case "
        "(Subtraction with 0) was not calculated properly.");

  // Test 4: Tests on Llint 63 Underflows and Overflows
  // Test 4.1: Testing Llint Underflow
  for(int i = 0; i < nums; ++i) {
    hostLongLongToInt63(INT_MIN, &arr_a.data()[i], &arr_a_ovrf.data()[i]);
    hostLongLongToInt63(1, &arr_b.data()[i], &arr_b_ovrf.data()[i]);
  }
  arr_a.upload();
  arr_a_ovrf.upload();
  arr_b.upload();
  arr_b_ovrf.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a.data(devc_layer),
                                                      arr_a_ovrf.data(devc_layer), 
                                                      arr_b.data(devc_layer),
                                                      arr_b_ovrf.data(devc_layer),
                                                      results.data(devc_layer),
                                                      results_ovrf.data(devc_layer),
                                                      Operator::SUBTRACT, nums);

  // Download results from device
  results.download();
  results_ovrf.download();
  llint_result.clear();
  llint_expected.clear();
  for(int i = 0; i < nums; i++){
    llint_result.push_back(hostInt63ToLongLong(results.readHost(i), results_ovrf.readHost(i)));
    llint_expected.push_back(static_cast<llint>(INT_MIN) - static_cast<llint>(1));
  }

  check(llint_result, RelationalOperator::EQUAL, llint_expected, "Error in calculating "
        "Llint -> Int64 Underflows.");
  
  // Test 4.1.1: Testing yet another underflow, with the ovrf populated
  for(int i = 0; i < nums; ++i) {
    arr_a_tuple.data()[i] = (hostLongLongToInt63(INT_MIN));
    arr_b_tuple.data()[i] = (hostLongLongToInt63(1));
  }
  arr_a_tuple.upload();
  arr_b_tuple.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a_tuple.data(devc_layer), 
                                                      arr_b.data(devc_layer),
                                                      arr_b_ovrf.data(devc_layer),
                                                      results.data(devc_layer),
                                                      results_ovrf.data(devc_layer),
                                                      Operator::SUBTRACT, nums);
  
  // Download results from device
  results.download();
  results_ovrf.download();
  llint_result.clear();
  llint_expected.clear();
  for(int i = 0; i < nums; i++){
    llint_result.push_back(hostInt63ToLongLong(results.readHost(i), results_ovrf.readHost(i)));
    llint_expected.push_back(static_cast<llint>(INT_MIN) - static_cast<llint>(1));
  }

  check(llint_result, RelationalOperator::EQUAL, llint_expected, "Error in calculating "
        "Llint -> Int64 Underflows. (second overload)");

  // Test 4.2: Testing Llint Overflow 
  for(int i = 0; i < nums; ++i) {
    hostLongLongToInt63(INT_MAX, &arr_a.data()[i], &arr_a_ovrf.data()[i]);
    hostLongLongToInt63(1, &arr_b.data()[i], &arr_b_ovrf.data()[i]);
  }
  arr_a.upload();
  arr_a_ovrf.upload();
  arr_b.upload();
  arr_b_ovrf.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a.data(devc_layer),
                                                      arr_a_ovrf.data(devc_layer), 
                                                      arr_b.data(devc_layer),
                                                      arr_b_ovrf.data(devc_layer),
                                                      results.data(devc_layer),
                                                      results_ovrf.data(devc_layer),
                                                      Operator::SUM, nums);

  // Download results from device
  results.download();
  results_ovrf.download();
  llint_result.clear();
  llint_expected.clear();
  for(int i = 0; i < nums; i++){
    llint_result.push_back(hostInt63ToLongLong(results.readHost(i), results_ovrf.readHost(i)));
    llint_expected.push_back(static_cast<llint>(INT_MAX) + static_cast<llint>(1));
  }

  check(llint_result, RelationalOperator::EQUAL, llint_expected, "Error in calculating "
        "Llint -> Int64 Overflows.");
  
  // Test 4.2.1: Testing yet another overflow, with the ovrf populated
  for(int i = 0; i < nums; ++i) {
    arr_a_tuple.data()[i] = (hostLongLongToInt63(INT_MAX));
    arr_b_tuple.data()[i] = (hostLongLongToInt63(1));
  }
  arr_a_tuple.upload();
  arr_b_tuple.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a_tuple.data(devc_layer), 
                                                      arr_b_tuple.data(devc_layer),
                                                      results.data(devc_layer),
                                                      results_ovrf.data(devc_layer),
                                                      Operator::SUM, nums);
  // Download results from device
  results.download();
  results_ovrf.download();
  llint_result.clear();
  llint_expected.clear();
  for(int i = 0; i < nums; i++){
    llint_result.push_back(hostInt63ToLongLong(results.readHost(i), results_ovrf.readHost(i)));
    llint_expected.push_back(static_cast<llint>(INT_MAX) + static_cast<llint>(1));
  }

  check(llint_result, RelationalOperator::EQUAL, llint_expected, "Error in calculating "
        "Llint -> Int64 Overflows. (second overload)");

  // Test 5: Basic tests on Double -> Int63
  for(int i = 0; i < nums; ++i) {
    hostDoubleToInt63(INT_MIN, &arr_a.data()[i], &arr_a_ovrf.data()[i]);
    hostDoubleToInt63(0, &arr_b.data()[i], &arr_b_ovrf.data()[i]);
  } 
  arr_a.upload();
  arr_a_ovrf.upload();
  arr_b.upload();
  arr_b_ovrf.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a.data(devc_layer),
                                                      arr_a_ovrf.data(devc_layer), 
                                                      arr_b.data(devc_layer),
                                                      arr_b_ovrf.data(devc_layer),
                                                      results.data(devc_layer),
                                                      results_ovrf.data(devc_layer),
                                                      Operator::SUBTRACT, nums);

  // Download results from device
  results.download();
  results_ovrf.download();
  std::vector<double> dbl_result, dbl_expected;
  for(int i = 0; i < nums; i++){
    dbl_result.push_back(hostInt63ToDouble(results.readHost(i), results_ovrf.readHost(i)));
    dbl_expected.push_back(static_cast<double>(INT_MIN)); 
  }

  check(dbl_result, RelationalOperator::EQUAL, dbl_expected, "Double Base Case was not "
        "calculated properly.");

  // Test 6: Tests on Double 63 Underflows and Overflows
  // Test 6.1: Testing Double Underflow
  for(int i = 0; i < nums; ++i) {
    hostDoubleToInt63(INT_MIN, &arr_a.data()[i], &arr_a_ovrf.data()[i]);
    hostDoubleToInt63(1, &arr_b.data()[i], &arr_b_ovrf.data()[i]);
  }
  arr_a.upload();
  arr_a_ovrf.upload();
  arr_b.upload();
  arr_b_ovrf.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a.data(devc_layer),
                                                      arr_a_ovrf.data(devc_layer), 
                                                      arr_b.data(devc_layer),
                                                      arr_b_ovrf.data(devc_layer),
                                                      results.data(devc_layer),
                                                      results_ovrf.data(devc_layer),
                                                      Operator::SUBTRACT, nums);

  // Download results from device
  results.download();
  results_ovrf.download();
  dbl_result.clear();
  dbl_expected.clear();
  for(int i = 0; i < nums; i++){
    dbl_result.push_back(hostInt63ToDouble(results.readHost(i), results_ovrf.readHost(i))), 
    dbl_expected.push_back(static_cast<double>(INT_MIN) - static_cast<double>(1));
  }

  check(dbl_result, RelationalOperator::EQUAL, dbl_expected, "Error in calculating "
        "Double -> Int64 Underflows.");
  
  // Test 6.1.1: Testing yet another overflow, with the ovrf populated
  for(int i = 0; i < nums; ++i) {
    arr_a_tuple.data()[i] = (hostDoubleToInt63(INT_MIN));
    arr_b_tuple.data()[i] = (hostLongLongToInt63(1));
  }
  arr_a_tuple.upload();
  arr_b_tuple.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a_tuple.data(devc_layer), 
                                                      arr_b_tuple.data(devc_layer),
                                                      results.data(devc_layer),
                                                      results_ovrf.data(devc_layer),
                                                      Operator::SUBTRACT, nums);

  // Download results from device
  results.download();
  results_ovrf.download();
  dbl_result.clear();
  dbl_expected.clear();
  for(int i = 0; i < nums; i++){
    dbl_result.push_back(hostInt63ToDouble(results.readHost(i), results_ovrf.readHost(i))), 
    dbl_expected.push_back(static_cast<double>(INT_MIN) - static_cast<double>(1));
  }

  check(dbl_result, RelationalOperator::EQUAL, dbl_expected, "Error in calculating "
        "Double -> Int64 Underflows. (second overload)");

  // Test 6.2: Testing Double Overflow 
  for(int i = 0; i < nums; ++i) {
    hostDoubleToInt63(INT_MAX, &arr_a.data()[i], &arr_a_ovrf.data()[i]);
    hostDoubleToInt63(1, &arr_b.data()[i], &arr_b_ovrf.data()[i]);
  }
  arr_a.upload();
  arr_a_ovrf.upload();
  arr_b.upload();
  arr_b_ovrf.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a.data(devc_layer),
                                                      arr_a_ovrf.data(devc_layer), 
                                                      arr_b.data(devc_layer),
                                                      arr_b_ovrf.data(devc_layer),
                                                      results.data(devc_layer),
                                                      results_ovrf.data(devc_layer),
                                                      Operator::SUM, nums);

  // Download results from device
  results.download();
  results_ovrf.download();
  dbl_result.clear();
  dbl_expected.clear();
  for(int i = 0; i < nums; i++){
    dbl_result.push_back(hostInt63ToDouble(results.readHost(i), results_ovrf.readHost(i)));
    dbl_expected.push_back(static_cast<double>(INT_MAX) + static_cast<double>(1));
  }

  check(dbl_result, RelationalOperator::EQUAL, dbl_expected, "Error in calculating "
        "Double -> Int64 Overflow.");
  
  // Test 6.2.1: Testing yet another overflow, with the ovrf populated
  for(int i = 0; i < nums; ++i) {
    arr_a_tuple.data()[i] = (hostDoubleToInt63(INT_MAX));
  }
  arr_a_tuple.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a_tuple.data(devc_layer), 
                                                      arr_b.data(devc_layer),
                                                      arr_b_ovrf.data(devc_layer),
                                                      results.data(devc_layer),
                                                      results_ovrf.data(devc_layer),
                                                      Operator::SUM, nums);
  // Download results from device
  results.download();
  results_ovrf.download();
  dbl_result.clear();
  dbl_expected.clear();
  for(int i = 0; i < nums; i++){
    dbl_result.push_back(hostInt63ToDouble(results.readHost(i), results_ovrf.readHost(i))), 
    dbl_expected.push_back(static_cast<double>(INT_MAX) + static_cast<double>(1));
  }

  check(dbl_result, RelationalOperator::EQUAL, dbl_expected, "Error in calculating "
        "Double -> Int64 Overflows. (second overload)");

  // // Test 7: Basic tests on Double -> Int95
  Hybrid<llint> arr_a_95(nums), arr_b_95(nums), results_95(nums);
  Hybrid<int> arr_a_95_ovrf(nums), arr_b_95_ovrf(nums), results_95_ovrf(nums);
  llint * arr_a_95_pointer = arr_a_95.data();
  llint * arr_b_95_pointer = arr_b_95.data();
  for(int i = 0; i < nums; ++i) {
    hostDoubleToInt95(LLONG_MIN, &arr_a_95.data()[i], &arr_a_95_ovrf.data()[i]);
    hostDoubleToInt95(0, &arr_b_95.data()[i], &arr_b_95_ovrf.data()[i]);
  } 
  arr_a_95.upload();
  arr_a_95_ovrf.upload();
  arr_b_95.upload();
  arr_b_95_ovrf.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a_95.data(devc_layer),
                                                      arr_a_95_ovrf.data(devc_layer), 
                                                      arr_b_95.data(devc_layer),
                                                      arr_b_95_ovrf.data(devc_layer),
                                                      results_95.data(devc_layer),
                                                      results_95_ovrf.data(devc_layer),
                                                      Operator::SUBTRACT, nums);

  // Download results from device
  results_95.download();
  results_95_ovrf.download();
  dbl_result.clear();
  dbl_expected.clear();
  for(int i = 0; i < nums; i++){
    dbl_result.push_back(hostInt95ToDouble(results_95.readHost(i), results_95_ovrf.readHost(i)));
    dbl_expected.push_back(static_cast<double>(LLONG_MIN));
  }

  check(dbl_result, RelationalOperator::EQUAL, dbl_expected, "Double -> Int95 Base "
        "Case fails.");

  // Test 8: Tests on Double 95 Underflows and Overflows
  // Test 8.1: Testing Double Underflow
  for(int i = 0; i < nums; ++i) {
    hostDoubleToInt95(LLONG_MIN, &arr_a_95_pointer[i], &arr_a_ovrf.data()[i]);
    hostDoubleToInt95(1, &arr_b_95_pointer[i], &arr_b_ovrf.data()[i]);
  }
  arr_a_95.upload();
  arr_a_ovrf.upload();
  arr_b_95.upload();
  arr_b_ovrf.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a_95.data(devc_layer),
                                                      arr_a_ovrf.data(devc_layer), 
                                                      arr_b_95.data(devc_layer),
                                                      arr_b_ovrf.data(devc_layer),
                                                      results_95.data(devc_layer),
                                                      results_95_ovrf.data(devc_layer),
                                                      Operator::SUBTRACT, nums);

  // Download results from device
  results_95.download();
  results_95_ovrf.download();
  dbl_result.clear();
  dbl_expected.clear();
  for(int i = 0; i < nums; i++){
    dbl_result.push_back(hostInt95ToDouble(results_95.readHost(i), results_95_ovrf.readHost(i)));
    dbl_expected.push_back(static_cast<double>(LLONG_MIN) - static_cast<double>(1));
  }

  check(dbl_result, RelationalOperator::EQUAL, dbl_expected, "Error in calculating "
        "Double -> Int95 Underflow.");

  // Test 8.2: Testing Double Overflow 
  for(int i = 0; i < nums; ++i) {
    hostDoubleToInt95(LLONG_MAX, &arr_a_95_pointer[i], &arr_a_ovrf.data()[i]);
    hostDoubleToInt95(1, &arr_b_95_pointer[i], &arr_b_ovrf.data()[i]);
  }
  arr_a_95.upload();
  arr_a_ovrf.upload();
  arr_b_95.upload();
  arr_b_ovrf.upload();
  testSplitFixedPrecision<<<nsmp, small_block_size>>>(arr_a_95.data(devc_layer),
                                                      arr_a_ovrf.data(devc_layer), 
                                                      arr_b_95.data(devc_layer),
                                                      arr_b_ovrf.data(devc_layer),
                                                      results_95.data(devc_layer),
                                                      results_95_ovrf.data(devc_layer),
                                                      Operator::SUM, nums);

  // Download results from device
  results_95.download();
  results_95_ovrf.download();
  dbl_result.clear();
  dbl_expected.clear();
  for(int i = 0; i < nums; i++){
    dbl_result.push_back(hostInt95ToDouble(results_95.readHost(i), results_95_ovrf.readHost(i)));
    dbl_expected.push_back(static_cast<double>(LLONG_MAX) + static_cast<double>(1));
  }
  check(dbl_result, RelationalOperator::EQUAL, dbl_expected, "Error in calculating "
        "Double -> Int95 Overflow.");
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
  testSplitFixedPrecision(nsmp);
  
  // Print results
  printTestSummary(oe.getVerbosity());
  
  return 0;
}
