// -*-c++-*-
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "hpc_random.h"
#include "hpc_random.cuh"
#include "random.h"

namespace stormm {
namespace random {

using card::HybridFormat;
using card::HybridTargetLevel;
  
#include "xor_shift_rng.cui"

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kInitXoroshiro128pArray(ullint2* state_vector, const int n_seeds, const int n_generators) {
  
  // The only threads that seed additional generators are those that can take a seed
  int seed_pos = threadIdx.x + (blockIdx.x * blockDim.x);
  if (seed_pos < n_seeds) {

    // Read the seed from the array
    ullint2 my_state = state_vector[seed_pos];
    state_vector[seed_pos] = my_state;
    for (int gen_pos = seed_pos + n_seeds; gen_pos < n_generators; gen_pos += n_seeds) {
      my_state = xoroshiro128p_jump(my_state);
      state_vector[gen_pos] = my_state;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void initXoroshiro128pArray(Hybrid<ullint2> *state_vector, const int igseed,
                                   const int scrub_cycles, const GpuDetails &gpu) {

  // Sanity check on the type of hybrid object
  ullint2* svdata;
  std::vector<ullint2> staging_space;
  const int n_generators = state_vector->size();
  const int n_seeds = std::min(max_xo_long_jumps, n_generators);
  const HybridFormat svfmt = state_vector->getFormat();
  switch (svfmt) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
    svdata = state_vector->data();
    break;
  case HybridFormat::HOST_ONLY:

    // If there is no allocated device data and the GPU cannot communicate directly into host
    // memory, this function was called in error.
    rtErr("Pageable host memory cannot be accessed by GPU kernels.  Use the " +
          getEnumerationName(HybridFormat::HOST_MOUNTED) + " memory format to allocate host-bound "
          "memory for the GPU to fill, or any of the formats with GPU-resident data.",
          "initXoroshiro128pArray");
  case HybridFormat::DEVICE_ONLY:

    // If there is no allocated host data, create a vector for staging the work
    staging_space.resize(n_seeds);
    svdata = staging_space.data();
    break;
  }
  Xoroshiro128pGenerator prng(igseed, scrub_cycles);
  svdata[0] = prng.revealState();
  for (int i = 1; i < n_seeds; i++) {
    prng.longJump();
    svdata[i] = prng.revealState();
  }
  switch (svfmt) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    state_vector->upload();
    break;
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
    break;
  case HybridFormat::DEVICE_ONLY:
    state_vector->putDevice(staging_space);
    break;
  }
  const int nsmp = gpu.getSMPCount();
  const int nthr = gpu.getMaxThreadsPerBlock();
  kInitXoroshiro128pArray<<<nsmp, nthr>>>(state_vector->data(HybridTargetLevel::DEVICE), n_seeds,
                                          n_generators);
  switch (svfmt) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    state_vector->download();
    break;
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::DEVICE_ONLY:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kInitXoshiro256ppArray(ullint2* state_xy, ullint2* state_zw, const int n_seeds,
                       const int n_generators) {
  
  // The only threads that seed additional generators are those that can take a seed
  int seed_pos = threadIdx.x + (blockIdx.x * blockDim.x);
  if (seed_pos < n_seeds) {
    
    // Read the seed from the array
    const ullint2 xy_state = state_xy[seed_pos];
    const ullint2 zw_state = state_zw[seed_pos];
    ullint4 my_state = { xy_state.x, xy_state.y, zw_state.x, zw_state.y };
    for (int gen_pos = seed_pos + n_seeds; gen_pos < n_generators; gen_pos += n_seeds) {
      my_state = xoshiro256pp_jump(my_state);
      state_xy[gen_pos] = { my_state.x, my_state.y };
      state_zw[gen_pos] = { my_state.z, my_state.w };
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void initXoshiro256ppArray(Hybrid<ullint2> *state_xy, Hybrid<ullint2> *state_zw,
                                  const int igseed, const int scrub_cycles,
                                  const GpuDetails &gpu) {

  // Sanity check on the type of hybrid object
  ullint2 *sv_xy_ptr, *sv_zw_ptr;
  std::vector<ullint2> staging_xy_space, staging_zw_space;
  const int n_generators = state_xy->size();
  const int n_seeds = std::min(max_xo_long_jumps, n_generators);
  const HybridFormat svfmt = state_xy->getFormat();
  switch (svfmt) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
    sv_xy_ptr = state_xy->data();
    sv_zw_ptr = state_zw->data();
    break;
  case HybridFormat::HOST_ONLY:
    rtErr("Pageable host memory cannot be accessed by GPU kernels.  Use the " +
          getEnumerationName(HybridFormat::HOST_MOUNTED) + " memory format to allocate host-bound "
          "memory for the GPU to fill, or any of the formats with GPU-resident data.",
          "initXoshiro256ppArray");
  case HybridFormat::DEVICE_ONLY:

    // If there is no allocated host data, create a vector for staging the work
    staging_xy_space.resize(n_seeds);
    staging_zw_space.resize(n_seeds);
    sv_xy_ptr = staging_xy_space.data();
    sv_zw_ptr = staging_zw_space.data();
    break;
  }
  Xoshiro256ppGenerator prng(igseed, scrub_cycles);
  const ullint4 seeded_state = prng.revealState();
  sv_xy_ptr[0] = { seeded_state.x, seeded_state.y };
  sv_zw_ptr[0] = { seeded_state.z, seeded_state.w };
  for (int i = 1; i < n_seeds; i++) {
    prng.longJump();
    const ullint4 jumped_state = prng.revealState();
    sv_xy_ptr[i] = { jumped_state.x, jumped_state.y };
    sv_zw_ptr[i] = { jumped_state.z, jumped_state.w };
  }
  switch (svfmt) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    state_xy->upload();
    state_zw->upload();
    break;
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
    break;
  case HybridFormat::DEVICE_ONLY:
    state_xy->putDevice(staging_xy_space);
    state_zw->putDevice(staging_zw_space);
    break;
  }
  const int nsmp = gpu.getSMPCount();
  const int nthr = gpu.getMaxThreadsPerBlock();
  kInitXoshiro256ppArray<<<nsmp, nthr>>>(state_xy->data(HybridTargetLevel::DEVICE),
                                         state_zw->data(HybridTargetLevel::DEVICE), n_seeds,
                                         n_generators);
  switch (svfmt) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    state_xy->download();
    state_zw->download();
    break;
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::DEVICE_ONLY:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kFillRandomCache(ullint2* state_xy, ullint2* state_zw, double* cache, const size_t length,
                 const size_t depth, const RandomAlgorithm method, const RandomNumberKind product,
                 const size_t index_start, const size_t index_end) {
  size_t pos = index_start + threadIdx.x + (blockIdx.x * blockDim.x);
  while (pos < index_end) {
    switch (method) {
    case RandomAlgorithm::XOROSHIRO_128P:
      {
        ullint2 my_state = state_xy[pos];
        for (size_t i = 0; i < depth; i++) {
          switch (product) {
          case RandomNumberKind::GAUSSIAN:
            cache[(i * length) + pos] = xoroshiro128p_normal(&my_state);
            break;
          case RandomNumberKind::UNIFORM:
            cache[(i * length) + pos] = xoroshiro128p_uniform(&my_state);
            break;
          }
        }
        state_xy[pos] = my_state;
      }
      break;
    case RandomAlgorithm::XOSHIRO_256PP:
      {
        const ullint2 xy_read = state_xy[pos];
        const ullint2 zw_read = state_zw[pos];
        ullint4 my_state = { xy_read.x, xy_read.y, zw_read.x, zw_read.y };
        for (size_t i = 0; i < depth; i++) {
          switch (product) {
          case RandomNumberKind::GAUSSIAN:
            cache[(i * length) + pos] = xoshiro256pp_normal(&my_state);
            break;
          case RandomNumberKind::UNIFORM:
            cache[(i * length) + pos] = xoshiro256pp_uniform(&my_state);
            break;
          }
        }
        const ullint2 xy_write = { my_state.x, my_state.y };
        const ullint2 zw_write = { my_state.z, my_state.w };
        state_xy[pos] = xy_write;
        state_zw[pos] = zw_write;
      }
      break;
    }
    pos += (blockDim.x * gridDim.x);
  }
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kFillRandomCache(ullint2* state_xy, ullint2* state_zw, float* cache, const size_t length,
                 const size_t depth, const RandomAlgorithm method, const RandomNumberKind product,
                 const size_t index_start, const size_t index_end) {
  size_t pos = index_start + threadIdx.x + (blockIdx.x * blockDim.x);
  while (pos < index_end) {
    switch (method) {
    case RandomAlgorithm::XOROSHIRO_128P:
      {
        ullint2 my_state = state_xy[pos];
        for (size_t i = 0; i < depth; i++) {
          switch (product) {
          case RandomNumberKind::GAUSSIAN:
            cache[(i * length) + pos] = xoroshiro128p_normalf(&my_state);
            break;
          case RandomNumberKind::UNIFORM:
            cache[(i * length) + pos] = xoroshiro128p_uniformf(&my_state);
            break;
          }
        }
        state_xy[pos] = my_state;
      }
      break;
    case RandomAlgorithm::XOSHIRO_256PP:
      {
        const ullint2 xy_read = state_xy[pos];
        const ullint2 zw_read = state_zw[pos];
        ullint4 my_state = { xy_read.x, xy_read.y, zw_read.x, zw_read.y };
        for (size_t i = 0; i < depth; i++) {
          switch (product) {
          case RandomNumberKind::GAUSSIAN:
            cache[(i * length) + pos] = xoshiro256pp_normalf(&my_state);
            break;
          case RandomNumberKind::UNIFORM:
            cache[(i * length) + pos] = xoshiro256pp_uniformf(&my_state);
            break;
          }
        }
        const ullint2 xy_write = { my_state.x, my_state.y };
        const ullint2 zw_write = { my_state.z, my_state.w };
        state_xy[pos] = xy_write;
        state_zw[pos] = zw_write;
      }
      break;
    }
    pos += (blockDim.x * gridDim.x);
  }
}

//-------------------------------------------------------------------------------------------------
extern void fillRandomCache(ullint2* state_xy, ullint2* state_zw, double* cache,
                            const size_t length, const size_t depth, const RandomAlgorithm method,
                            const RandomNumberKind product, const size_t index_start,
                            const size_t index_end, const GpuDetails &gpu) {
  kFillRandomCache<<<gpu.getSMPCount(), large_block_size>>>(state_xy, state_zw, cache, length,
                                                            depth, method, product, index_start,
                                                            index_end);
}

//-------------------------------------------------------------------------------------------------
extern void fillRandomCache(ullint2* state_xy, ullint2* state_zw, float* cache,
                            const size_t length, const size_t depth, const RandomAlgorithm method,
                            const RandomNumberKind product, const size_t index_start,
                            const size_t index_end, const GpuDetails &gpu) {
  kFillRandomCache<<<gpu.getSMPCount(), large_block_size>>>(state_xy, state_zw, cache, length,
                                                            depth, method, product, index_start,
                                                            index_end);
}

//-------------------------------------------------------------------------------------------------
extern void fillRandomCache(Hybrid<ullint2> *state_xy, Hybrid<ullint2> *state_zw,
                            Hybrid<double> *cache, const size_t length, const size_t depth,
                            const RandomAlgorithm method, const RandomNumberKind product,
                            const size_t index_start, const size_t index_end,
                            const GpuDetails &gpu) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  kFillRandomCache<<<gpu.getSMPCount(),
                     large_block_size>>>(state_xy->data(tier), state_zw->data(tier),
                                         cache->data(tier), length, depth, method, product,
                                         index_start, index_end);
}

//-------------------------------------------------------------------------------------------------
extern void fillRandomCache(Hybrid<ullint2> *state_xy, Hybrid<ullint2> *state_zw,
                            Hybrid<float> *cache, const size_t length, const size_t depth,
                            const RandomAlgorithm method, const RandomNumberKind product,
                            const size_t index_start, const size_t index_end,
                            const GpuDetails &gpu) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  kFillRandomCache<<<gpu.getSMPCount(),
                     large_block_size>>>(state_xy->data(tier), state_zw->data(tier),
                                         cache->data(tier), length, depth, method, product,
                                         index_start, index_end);
}

} // namespace random
} // namespace stormm
