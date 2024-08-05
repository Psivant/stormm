// -*-c++-*-
#ifndef STORMM_SUMMATION_HPC_H
#define STORMM_SUMMATION_HPC_H

#include <cuda.h>
#include <cuda_runtime.h>
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/ptx_macros.h"
#include "Math/summation.h"
#include "Reporting/error_format.h"

namespace stormm {
namespace hpc_math {

using card::Hybrid;
using card::HybridFormat;
using card::HybridTargetLevel;
using card::GpuDetails;
  
/// \brief Kernel for summing a vector of scalar elements (i.e. double, float, short unsigned int)
///
/// \param vdata   Pointer to the data array allocated on the GPU
/// \param length  Length of the data array (trusted)
/// \param result  Pointer to an array of block accumulators (overwritten by this kernel)
template <typename TSum, typename TBase>
__global__ void __launch_bounds__(large_block_size, 1)
kSumVector(const TBase* vdata, const size_t length, TSum* result) {
  __shared__ TSum warp_sums[(large_block_size >> warp_bits)];
  TSum thread_sum = (TSum)(0);
  const size_t ginc = gridDim.x * blockDim.x;
  for (size_t pos = (blockIdx.x * blockDim.x) + threadIdx.x; pos < length; pos += ginc) {
    thread_sum += static_cast<TSum>(vdata[pos]);
  }
#if STORMM_USE_HIP
  thread_sum += SHFL_DOWN(thread_sum, 32);
#endif
  thread_sum += SHFL_DOWN(thread_sum, 16);
  thread_sum += SHFL_DOWN(thread_sum,  8);
  thread_sum += SHFL_DOWN(thread_sum,  4);
  thread_sum += SHFL_DOWN(thread_sum,  2);
  thread_sum += SHFL_DOWN(thread_sum,  1);
  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (lane_idx == 0) {
    warp_sums[warp_idx] = thread_sum;
  }
  __syncthreads();
  if (warp_idx == 0) {

    // With large_block_size set as it is and warp sizes being what they are, HIP can skip
    // one of these shuffle operations.  HIP has an extra shuffle operation above, CUDA has
    // an extra one down here.  A check on the warp index is still needed to avoid overrunning
    // the warp_sums array.
    thread_sum = (lane_idx < (large_block_size >> warp_bits)) ? warp_sums[lane_idx] : (TSum)(0);
#ifdef STORMM_USE_CUDA
    thread_sum += SHFL_DOWN(thread_sum, 16);
#endif
    thread_sum += SHFL_DOWN(thread_sum,  8);
    thread_sum += SHFL_DOWN(thread_sum,  4);
    thread_sum += SHFL_DOWN(thread_sum,  2);
    thread_sum += SHFL_DOWN(thread_sum,  1);
    if (lane_idx == 0) {
      result[blockIdx.x] = thread_sum;
    }
  }
}

/// \brief Launch the appropriately templated vector summation kernel (scalar data elements)
///
/// \param hb      The data vector to reduce
/// \param buffer  Temporary, device-mappable array of host memory (HOST_MOUNTED) to hold the
///                results coming out of each thread block
/// \param gpu     Details of the GPU
template <typename TSum, typename TBase>
TSum launchSumVector(const Hybrid<TBase> &hb, Hybrid<TSum> *buffer, const GpuDetails &gpu) {

  // Check the provided result vector: does it have the correct host-mapped memory format?
  if (buffer->getFormat() != HybridFormat::HOST_MOUNTED) {
    rtErr("The results of an HPC summation must be buffered in a GPU-mappable Hybrid object.  "
          "Choose HOST_MOUNTED as the format for " + std::string(buffer->getLabel().name) + ".",
          "launchSumVector");
  }

  // Use the device to compute a series of sums for each block.  The results will be sent back to
  // the host-mapped memory via the PCIE bus.
  const int nsmp = gpu.getSMPCount();
  const int nthr = gpu.getMaxThreadsPerBlock();
  kSumVector<TSum><<<nsmp, nthr>>>(hb.data(HybridTargetLevel::DEVICE), hb.size(),
                                   buffer->getDeviceValidHostPointer());
  cudaDeviceSynchronize();
  
  // Sum the results on the host rather than launching a new kernel (which would have much more
  // latency).
  TSum result = static_cast<TSum>(0);
  TSum* buffdata = buffer->data();
  for (int i = 0; i < nsmp; i++) {
    result += buffdata[i];
  }
  return result;
}

/// \brief Templated function for taking a sum over scalar data elements on the GPU.  This will
///        call the equivalent C++ function if the host data, rather than that which is found on
///        the device, is to be summed.
///
/// \param  hb      The data to sum
/// \param  buffer  Array of accumulators
/// \param  gpu     Details of the GPU to utilize (which must also be the one where any device
///                 data resides)
/// \param  tier    Level at which to perform the summation
template <typename TSum, typename TBase>
TSum sum(const Hybrid<TBase> &hb, Hybrid<TSum> *buffer, const GpuDetails &gpu,
         const HybridTargetLevel tier = HybridTargetLevel::DEVICE) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return stmath::sum<TSum>(hb);
  case HybridTargetLevel::DEVICE:

    // Test the provided GPU: is it real?
    if (gpu.getGpuSupported()) {
      if (buffer == nullptr || buffer->size() < gpu.getSMPCount()) {
        Hybrid<TSum> tmp_buffer(gpu.getSMPCount(), "gpu_block_sums");
        return launchSumVector<TSum>(hb, &tmp_buffer, gpu);
      }
      else {
        return launchSumVector<TSum>(hb, buffer, gpu);
      }
    }
    else {
      rtErr("Unsupported GPU.\n");
    }
  }
  __builtin_unreachable();
}

/// \brief Kernel for summing a vector of two-tuples (i.e. double2, float2, ushort2)
///
/// \param vdata   Pointer to the data array allocated on the GPU
/// \param length  Length of the data array (trusted)
/// \param result  Pointer to an array of block accumulators (overwritten by this kernel)
template <typename TSum, typename TBase>
__global__ void __launch_bounds__(large_block_size, 1)
kSumVectorTuple2(const TBase* vdata, const size_t length, TSum* result) {
  __shared__ TSum warp_sums[(large_block_size >> warp_bits)];
  TSum thread_sum = { (TSum)(0), (TSum)(0) };
  const size_t ginc = gridDim.x * blockDim.x;
  for (size_t pos = (blockIdx.x * blockDim.x) + threadIdx.x; pos < length; pos += ginc) {
    TBase vdpos = vdata[pos];
    thread_sum.x += vdpos.x;
    thread_sum.y += vdpos.y;
  }
#if STORMM_USE_HIP
  thread_sum.x += SHFL_DOWN(thread_sum.x, 32);
  thread_sum.y += SHFL_DOWN(thread_sum.y, 32);
#endif
  thread_sum.x += SHFL_DOWN(thread_sum.x, 16);
  thread_sum.y += SHFL_DOWN(thread_sum.y, 16);
  thread_sum.x += SHFL_DOWN(thread_sum.x,  8);
  thread_sum.y += SHFL_DOWN(thread_sum.y,  8);
  thread_sum.x += SHFL_DOWN(thread_sum.x,  4);
  thread_sum.y += SHFL_DOWN(thread_sum.y,  4);
  thread_sum.x += SHFL_DOWN(thread_sum.x,  2);
  thread_sum.y += SHFL_DOWN(thread_sum.y,  2);
  thread_sum.x += SHFL_DOWN(thread_sum.x,  1);
  thread_sum.y += SHFL_DOWN(thread_sum.y,  1);
  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (lane_idx == 0) {
    warp_sums[warp_idx].x = thread_sum.x;
    warp_sums[warp_idx].y = thread_sum.y;
  }
  __syncthreads();
  if (warp_idx == 0) {

    // With large_block_size set as it is and warp sizes being what they are, HIP can skip
    // one of these shuffle operations.  HIP has an extra shuffle operation above, CUDA has
    // an extra one down here.  A check on the warp index is still needed to avoid overrunning
    // the warp_sums array.
    if (lane_idx < (large_block_size >> warp_bits)) {
      thread_sum.x = warp_sums[lane_idx].x;
      thread_sum.y = warp_sums[lane_idx].y;
    }
    else {
      thread_sum.x = (TSum)(0);
      thread_sum.y = (TSum)(0);
    }
#ifdef STORMM_USE_CUDA
    thread_sum.x += SHFL_DOWN(thread_sum.x, 16);
    thread_sum.y += SHFL_DOWN(thread_sum.y, 16);
#endif
    thread_sum.x += SHFL_DOWN(thread_sum.x,  8);
    thread_sum.y += SHFL_DOWN(thread_sum.y,  8);
    thread_sum.x += SHFL_DOWN(thread_sum.x,  4);
    thread_sum.y += SHFL_DOWN(thread_sum.y,  4);
    thread_sum.x += SHFL_DOWN(thread_sum.x,  2);
    thread_sum.y += SHFL_DOWN(thread_sum.y,  2);
    thread_sum.x += SHFL_DOWN(thread_sum.x,  1);
    thread_sum.y += SHFL_DOWN(thread_sum.y,  1);
    if (lane_idx == 0) {
      TSum tmp_result = { thread_sum.x, thread_sum.y };
      result[blockIdx.x] = tmp_result;
    }
  }
}

/// \brief Launch the appropriately templated vector summation kernel (data elements are tuples of
///        two scalars, i.e. double2, ushort2, int2)
///
/// \param hb      The data vector to reduce
/// \param buffer  Temporary, device-mappable array of host memory (HOST_MOUNTED) to hold the
///                results coming out of each thread block
/// \param gpu     Details of the GPU
template <typename TSum, typename TBase>
TSum launchSumVectorTuple2(const Hybrid<TBase> &hb, Hybrid<TSum> *buffer, const GpuDetails &gpu) {

  // Check the provided result vector: does it have the correct host-mapped memory format?
  if (buffer->getFormat() != HybridFormat::HOST_MOUNTED) {
    rtErr("The results of an HPC summation must be buffered in a GPU-mappable Hybrid object.  "
          "Choose HOST_MOUNTED as the format for " + std::string(buffer->getLabel().name) + ".",
          "launchSumVectorTuple2");
  }

  // Use the device to compute a series of sums for each block.  The results will be sent back to
  // the host-mapped memory via the PCIE bus.
  const int nsmp = gpu.getSMPCount();
  const int nthr = gpu.getMaxThreadsPerBlock();
  kSumVectorTuple2<TSum><<<nsmp, nthr>>>(hb.data(HybridTargetLevel::DEVICE), hb.size(),
                                         buffer->data(HybridTargetLevel::DEVICE));
  cudaDeviceSynchronize();
  
  // Sum the results on the host rather than launching a new kernel (which would have much more
  // latency).
  TSum result = { static_cast<TSum>(0), static_cast<TSum>(0) };
  TSum* buffdata = buffer->data();
  for (int i = 0; i < nsmp; i++) {
    result.x += buffdata[i].x;
    result.y += buffdata[i].y;
  }
  return result;
}

/// \brief Templated function for taking a sum over scalar data elements on the GPU.  This will
///        call the equivalent C++ function if the host data, rather than that which is found on
///        the device, is to be summed.
///
/// \param  hb      The data to sum
/// \param  buffer  Array of accumulators
/// \param  gpu     Details of the GPU to utilize (which must also be the one where any device
///                 data resides)
/// \param  tier    Level at which to perform the summation
template <typename TSum, typename TBase>
TSum sumTuple2(const Hybrid<TBase> &hb, Hybrid<TSum> *buffer, const GpuDetails &gpu,
               const HybridTargetLevel tier = HybridTargetLevel::DEVICE) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return stmath::sumTuple2<TSum>(hb);
  case HybridTargetLevel::DEVICE:

    // Test the provided GPU: is it real?
    if (gpu.getGpuSupported()) {
      if (buffer == nullptr || buffer->size() < gpu.getSMPCount()) {
        Hybrid<TSum> tmp_buffer(gpu.getSMPCount(), "gpu_block_sums");
        return launchSumVectorTuple2<TSum>(hb, &tmp_buffer, gpu);
      }
      else {
        return launchSumVectorTuple2<TSum>(hb, buffer, gpu);
      }
    }
    else {
      rtErr("Unsupported GPU.\n");
    }
  }
  __builtin_unreachable();
}
  
/// \brief Kernel for summing a vector of three-tuples (i.e. double3, float3, ushort3)
///
/// \param vdata   Pointer to the data array allocated on the GPU
/// \param length  Length of the data array (trusted)
/// \param result  Pointer to an array of block accumulators (overwritten by this kernel)
template <typename TSum, typename TBase>
__global__ void __launch_bounds__(large_block_size, 1)
kSumVectorTuple3(const TBase* vdata, const size_t length, TSum* result) {
  __shared__ TSum warp_sums[(large_block_size >> warp_bits)];
  TSum thread_sum = { (TSum)(0), (TSum)(0), (TSum)(0) };
  const size_t ginc = gridDim.x * blockDim.x;
  for (size_t pos = (blockIdx.x * blockDim.x) + threadIdx.x; pos < length; pos += ginc) {
    TBase vdpos = vdata[pos];
    thread_sum.x += vdpos.x;
    thread_sum.y += vdpos.y;
    thread_sum.z += vdpos.z;
  }
#if STORMM_USE_HIP
  thread_sum.x += SHFL_DOWN(thread_sum.x, 32);
  thread_sum.y += SHFL_DOWN(thread_sum.y, 32);
  thread_sum.z += SHFL_DOWN(thread_sum.z, 32);
#endif
  thread_sum.x += SHFL_DOWN(thread_sum.x, 16);
  thread_sum.y += SHFL_DOWN(thread_sum.y, 16);
  thread_sum.z += SHFL_DOWN(thread_sum.z, 16);
  thread_sum.x += SHFL_DOWN(thread_sum.x,  8);
  thread_sum.y += SHFL_DOWN(thread_sum.y,  8);
  thread_sum.z += SHFL_DOWN(thread_sum.z,  8);
  thread_sum.x += SHFL_DOWN(thread_sum.x,  4);
  thread_sum.y += SHFL_DOWN(thread_sum.y,  4);
  thread_sum.z += SHFL_DOWN(thread_sum.z,  4);
  thread_sum.x += SHFL_DOWN(thread_sum.x,  2);
  thread_sum.y += SHFL_DOWN(thread_sum.y,  2);
  thread_sum.z += SHFL_DOWN(thread_sum.z,  2);
  thread_sum.x += SHFL_DOWN(thread_sum.x,  1);
  thread_sum.y += SHFL_DOWN(thread_sum.y,  1);
  thread_sum.z += SHFL_DOWN(thread_sum.z,  1);
  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (lane_idx == 0) {
    warp_sums[warp_idx].x = thread_sum.x;
    warp_sums[warp_idx].y = thread_sum.y;
    warp_sums[warp_idx].z = thread_sum.z;
  }
  __syncthreads();
  if (warp_idx == 0) {

    // With large_block_size set as it is and warp sizes being what they are, HIP can skip
    // one of these shuffle operations.  HIP has an extra shuffle operation above, CUDA has
    // an extra one down here.  A check on the warp index is still needed to avoid overrunning
    // the warp_sums array.
    if (lane_idx < (large_block_size >> warp_bits)) {
      thread_sum.x = warp_sums[lane_idx].x;
      thread_sum.y = warp_sums[lane_idx].y;
      thread_sum.z = warp_sums[lane_idx].z;
    }
    else {
      thread_sum.x = (TSum)(0);
      thread_sum.y = (TSum)(0);
      thread_sum.z = (TSum)(0);
    }
#ifdef STORMM_USE_CUDA
    thread_sum.x += SHFL_DOWN(thread_sum.x, 16);
    thread_sum.y += SHFL_DOWN(thread_sum.y, 16);
    thread_sum.z += SHFL_DOWN(thread_sum.z, 16);
#endif
    thread_sum.x += SHFL_DOWN(thread_sum.x,  8);
    thread_sum.y += SHFL_DOWN(thread_sum.y,  8);
    thread_sum.z += SHFL_DOWN(thread_sum.z,  8);
    thread_sum.x += SHFL_DOWN(thread_sum.x,  4);
    thread_sum.y += SHFL_DOWN(thread_sum.y,  4);
    thread_sum.z += SHFL_DOWN(thread_sum.z,  4);
    thread_sum.x += SHFL_DOWN(thread_sum.x,  2);
    thread_sum.y += SHFL_DOWN(thread_sum.y,  2);
    thread_sum.z += SHFL_DOWN(thread_sum.z,  2);
    thread_sum.x += SHFL_DOWN(thread_sum.x,  1);
    thread_sum.y += SHFL_DOWN(thread_sum.y,  1);
    thread_sum.z += SHFL_DOWN(thread_sum.z,  1);
    if (lane_idx == 0) {
      TSum tmp_result = { thread_sum.x, thread_sum.y, thread_sum.z };
      result[blockIdx.x] = tmp_result;
    }
  }
}

/// \brief Launch the appropriately templated vector summation kernel (data elements are tuples of
///        three scalars, i.e. double3, ushort3, int3)
///
/// \param hb      The data vector to reduce
/// \param buffer  Temporary, device-mappable array of host memory (HOST_MOUNTED) to hold the
///                results coming out of each thread block
/// \param gpu     Details of the GPU
template <typename TSum, typename TBase>
TSum launchSumVectorTuple3(const Hybrid<TBase> &hb, Hybrid<TSum> *buffer, const GpuDetails &gpu) {

  // Check the provided result vector: does it have the correct host-mapped memory format?
  if (buffer->getFormat() != HybridFormat::HOST_MOUNTED) {
    rtErr("The results of an HPC summation must be buffered in a GPU-mappable Hybrid object.  "
          "Choose HOST_MOUNTED as the format for " + std::string(buffer->getLabel().name) + ".",
          "launchSumVectorTuple3");
  }

  // Use the device to compute a series of sums for each block.  The results will be sent back to
  // the host-mapped memory via the PCIE bus.
  const int nsmp = gpu.getSMPCount();
  const int nthr = gpu.getMaxThreadsPerBlock();
  kSumVectorTuple3<TSum><<<nsmp, nthr>>>(hb.data(HybridTargetLevel::DEVICE), hb.size(),
                                         buffer->data(HybridTargetLevel::DEVICE));
  cudaDeviceSynchronize();
  
  // Sum the results on the host rather than launching a new kernel (which would have much more
  // latency).
  TSum result = { static_cast<TSum>(0), static_cast<TSum>(0), static_cast<TSum>(0) };
  TSum* buffdata = buffer->data();
  for (int i = 0; i < nsmp; i++) {
    result.x += buffdata[i].x;
    result.y += buffdata[i].y;
    result.z += buffdata[i].z;
  }
  return result;
}

/// \brief Templated function for taking a sum over scalar data elements on the GPU.  This will
///        call the equivalent C++ function if the host data, rather than that which is found on
///        the device, is to be summed.
///
/// \param  hb      The data to sum
/// \param  buffer  Array of accumulators
/// \param  gpu     Details of the GPU to utilize (which must also be the one where any device
///                 data resides)
/// \param  tier    Level at which to perform the summation
template <typename TSum, typename TBase>
TSum sumTuple3(const Hybrid<TBase> &hb, Hybrid<TSum> *buffer, const GpuDetails &gpu,
               const HybridTargetLevel tier = HybridTargetLevel::DEVICE) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return stmath::sumTuple3<TSum>(hb);
  case HybridTargetLevel::DEVICE:

    // Test the provided GPU: is it real?
    if (gpu.getGpuSupported()) {
      if (buffer == nullptr || buffer->size() < gpu.getSMPCount()) {
        Hybrid<TSum> tmp_buffer(gpu.getSMPCount(), "gpu_block_sums");
        return launchSumVectorTuple3<TSum>(hb, &tmp_buffer, gpu);
      }
      else {
        return launchSumVectorTuple3<TSum>(hb, buffer, gpu);
      }
    }
    else {
      rtErr("Unsupported GPU.\n");
    }
  }
  __builtin_unreachable();
}

/// \brief Kernel for summing a vector of four-tuples (i.e. double4, float4, ushort4)
///
/// \param vdata   Pointer to the data array allocated on the GPU
/// \param length  Length of the data array (trusted)
/// \param result  Pointer to an array of block accumulators (overwritten by this kernel)
template <typename TSum, typename TBase>
__global__ void __launch_bounds__(large_block_size, 1)
kSumVectorTuple4(const TBase* vdata, const size_t length, TSum* result) {
  __shared__ TSum warp_sums[(large_block_size >> warp_bits)];
  TSum thread_sum = { (TSum)(0), (TSum)(0), (TSum)(0) };
  const size_t ginc = gridDim.x * blockDim.x;
  for (size_t pos = (blockIdx.x * blockDim.x) + threadIdx.x; pos < length; pos += ginc) {
    TBase vdpos = vdata[pos];
    thread_sum.x += vdpos.x;
    thread_sum.y += vdpos.y;
    thread_sum.z += vdpos.z;
    thread_sum.w += vdpos.w;
  }
#if STORMM_USE_HIP
  thread_sum.x += SHFL_DOWN(thread_sum.x, 32);
  thread_sum.y += SHFL_DOWN(thread_sum.y, 32);
  thread_sum.z += SHFL_DOWN(thread_sum.z, 32);
  thread_sum.w += SHFL_DOWN(thread_sum.w, 32);
#endif
  thread_sum.x += SHFL_DOWN(thread_sum.x, 16);
  thread_sum.y += SHFL_DOWN(thread_sum.y, 16);
  thread_sum.z += SHFL_DOWN(thread_sum.z, 16);
  thread_sum.w += SHFL_DOWN(thread_sum.w, 16);
  thread_sum.x += SHFL_DOWN(thread_sum.x,  8);
  thread_sum.y += SHFL_DOWN(thread_sum.y,  8);
  thread_sum.z += SHFL_DOWN(thread_sum.z,  8);
  thread_sum.w += SHFL_DOWN(thread_sum.w,  8);
  thread_sum.x += SHFL_DOWN(thread_sum.x,  4);
  thread_sum.y += SHFL_DOWN(thread_sum.y,  4);
  thread_sum.z += SHFL_DOWN(thread_sum.z,  4);
  thread_sum.w += SHFL_DOWN(thread_sum.w,  4);
  thread_sum.x += SHFL_DOWN(thread_sum.x,  2);
  thread_sum.y += SHFL_DOWN(thread_sum.y,  2);
  thread_sum.z += SHFL_DOWN(thread_sum.z,  2);
  thread_sum.w += SHFL_DOWN(thread_sum.w,  2);
  thread_sum.x += SHFL_DOWN(thread_sum.x,  1);
  thread_sum.y += SHFL_DOWN(thread_sum.y,  1);
  thread_sum.z += SHFL_DOWN(thread_sum.z,  1);
  thread_sum.w += SHFL_DOWN(thread_sum.w,  1);
  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (lane_idx == 0) {
    warp_sums[warp_idx].x = thread_sum.x;
    warp_sums[warp_idx].y = thread_sum.y;
    warp_sums[warp_idx].z = thread_sum.z;
    warp_sums[warp_idx].w = thread_sum.w;
  }
  __syncthreads();
  if (warp_idx == 0) {

    // With large_block_size set as it is and warp sizes being what they are, HIP can skip
    // one of these shuffle operations.  HIP has an extra shuffle operation above, CUDA has
    // an extra one down here.  A check on the warp index is still needed to avoid overrunning
    // the warp_sums array.
    if (lane_idx < (large_block_size >> warp_bits)) {
      thread_sum.x = warp_sums[lane_idx].x;
      thread_sum.y = warp_sums[lane_idx].y;
      thread_sum.z = warp_sums[lane_idx].z;
      thread_sum.w = warp_sums[lane_idx].w;
    }
    else {
      thread_sum.x = (TSum)(0);
      thread_sum.y = (TSum)(0);
      thread_sum.z = (TSum)(0);
      thread_sum.w = (TSum)(0);
    }
#ifdef STORMM_USE_CUDA
    thread_sum.x += SHFL_DOWN(thread_sum.x, 16);
    thread_sum.y += SHFL_DOWN(thread_sum.y, 16);
    thread_sum.z += SHFL_DOWN(thread_sum.z, 16);
    thread_sum.w += SHFL_DOWN(thread_sum.w, 16);
#endif
    thread_sum.x += SHFL_DOWN(thread_sum.x,  8);
    thread_sum.y += SHFL_DOWN(thread_sum.y,  8);
    thread_sum.z += SHFL_DOWN(thread_sum.z,  8);
    thread_sum.w += SHFL_DOWN(thread_sum.w,  8);
    thread_sum.x += SHFL_DOWN(thread_sum.x,  4);
    thread_sum.y += SHFL_DOWN(thread_sum.y,  4);
    thread_sum.z += SHFL_DOWN(thread_sum.z,  4);
    thread_sum.w += SHFL_DOWN(thread_sum.w,  4);
    thread_sum.x += SHFL_DOWN(thread_sum.x,  2);
    thread_sum.y += SHFL_DOWN(thread_sum.y,  2);
    thread_sum.z += SHFL_DOWN(thread_sum.z,  2);
    thread_sum.w += SHFL_DOWN(thread_sum.w,  2);
    thread_sum.x += SHFL_DOWN(thread_sum.x,  1);
    thread_sum.y += SHFL_DOWN(thread_sum.y,  1);
    thread_sum.z += SHFL_DOWN(thread_sum.z,  1);
    thread_sum.w += SHFL_DOWN(thread_sum.w,  1);
    if (lane_idx == 0) {
      TSum tmp_result = { thread_sum.x, thread_sum.y, thread_sum.z, thread_sum.w };
      result[blockIdx.x] = tmp_result;
    }
  }
}

/// \brief Launch the appropriately templated vector summation kernel (data elements are tuples of
///        four scalars, i.e. double4, ushort4, int4)
///
/// \param hb      The data vector to reduce
/// \param buffer  Temporary, device-mappable array of host memory (HOST_MOUNTED) to hold the
///                results coming out of each thread block
/// \param gpu     Details of the GPU
template <typename TSum, typename TBase>
TSum launchSumVectorTuple4(const Hybrid<TBase> &hb, Hybrid<TSum> *buffer, const GpuDetails &gpu) {

  // Check the provided result vector: does it have the correct host-mapped memory format?
  if (buffer->getFormat() != HybridFormat::HOST_MOUNTED) {
    rtErr("The results of an HPC summation must be buffered in a GPU-mappable Hybrid object.  "
          "Choose HOST_MOUNTED as the format for " + std::string(buffer->getLabel().name) + ".",
          "launchSumVectorTuple4");
  }

  // Use the device to compute a series of sums for each block.  The results will be sent back to
  // the host-mapped memory via the PCIE bus.
  const int nsmp = gpu.getSMPCount();
  const int nthr = gpu.getMaxThreadsPerBlock();
  kSumVectorTuple4<TSum><<<nsmp, nthr>>>(hb.data(HybridTargetLevel::DEVICE), hb.size(),
                                         buffer->data(HybridTargetLevel::DEVICE));
  cudaDeviceSynchronize();
  
  // Sum the results on the host rather than launching a new kernel (which would have much more
  // latency).
  TSum result = { static_cast<TSum>(0), static_cast<TSum>(0), static_cast<TSum>(0),
                  static_cast<TSum>(0) };
  TSum* buffdata = buffer->data();
  for (int i = 0; i < nsmp; i++) {
    result.x += buffdata[i].x;
    result.y += buffdata[i].y;
    result.z += buffdata[i].z;
    result.w += buffdata[i].w;
  }
  return result;
}

/// \brief Templated function for taking a sum over scalar data elements on the GPU.  This will
///        call the equivalent C++ function if the host data, rather than that which is found on
///        the device, is to be summed.
///
/// \param  hb      The data to sum
/// \param  buffer  Array of accumulators
/// \param  gpu     Details of the GPU to utilize (which must also be the one where any device
///                 data resides)
/// \param  tier    Level at which to perform the summation
template <typename TSum, typename TBase>
TSum sumTuple4(const Hybrid<TBase> &hb, Hybrid<TSum> *buffer, const GpuDetails &gpu,
               const HybridTargetLevel tier = HybridTargetLevel::DEVICE) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return stmath::sumTuple4<TSum>(hb);
  case HybridTargetLevel::DEVICE:

    // Test the provided GPU: is it real?
    if (gpu.getGpuSupported()) {
      if (buffer == nullptr || buffer->size() < gpu.getSMPCount()) {
        Hybrid<TSum> tmp_buffer(gpu.getSMPCount(), "gpu_block_sums");
        return launchSumVectorTuple4<TSum>(hb, &tmp_buffer, gpu);
      }
      else {
        return launchSumVectorTuple4<TSum>(hb, buffer, gpu);
      }
    }
    else {
      rtErr("Unsupported GPU.\n");
    }
  }
  __builtin_unreachable();
}
  
} // namespace hpc_math
} // namespace stormm

#endif
