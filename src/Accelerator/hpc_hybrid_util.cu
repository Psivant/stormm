// -*-c++-*-
#include "copyright.h"
#include "gpu_details.h"
#include "hpc_hybrid_util.cuh"
#include "hpc_hybrid_util.h"
#include "hybrid.h"

namespace stormm {
namespace card {

//-------------------------------------------------------------------------------------------------
__global__ __launch_bounds__(medium_block_size, 2)
void kDeepCopy(llint* dest_host, int* dest_host_ovrf, llint* dest_devc, int* dest_devc_ovrf,
               const llint* orig_host, const int* orig_host_ovrf, const llint* orig_devc,
               const int* orig_devc_ovrf, const size_t dest_offset, const size_t orig_offset,
               const size_t length, const int dest_bits, const int orig_bits, const bool do_hdc,
               const bool do_dhc, const bool do_ddc) {
  const size_t stride = blockDim.x * gridDim.x;
  for (size_t tid = threadIdx.x + (blockIdx.x * blockDim.x); tid < length; tid += stride) {
    if (do_hdc) {
      const int95_t val = changeFPBits(orig_host[orig_offset + tid],
                                       orig_host_ovrf[orig_offset + tid], orig_bits, dest_bits);
      dest_devc[dest_offset + tid]      = val.x;
      dest_devc_ovrf[dest_offset + tid] = val.y;
    }
    if (do_dhc) {
      const int95_t val = changeFPBits(orig_devc[orig_offset + tid],
                                       orig_devc_ovrf[orig_offset + tid], orig_bits, dest_bits);
      dest_host[dest_offset + tid]      = val.x;
      dest_host_ovrf[dest_offset + tid] = val.y;
    }
    if (do_ddc) {
      const int95_t val = changeFPBits(orig_devc[orig_offset + tid],
                                       orig_devc_ovrf[orig_offset + tid], orig_bits, dest_bits);
      dest_devc[dest_offset + tid]      = val.x;
      dest_devc_ovrf[dest_offset + tid] = val.y;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void launchDeepCopy(void* vdest_host, void* vdest_devc, const void* vorig_host,
                    const void* vorig_devc, const size_t dest_offset, const size_t orig_offset,
                    const size_t length, const size_t ct, const bool do_hdc, const bool do_dhc,
                    const bool do_ddc, const GpuDetails &gpu, const int dest_bits,
                    const int orig_bits) {
  if (ct == float_type_index) {
    unrollLaunchDeepCopy<float>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == double_type_index) {
    unrollLaunchDeepCopy<double>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == bool_type_index) {
    unrollLaunchDeepCopy<bool>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                               orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == char_type_index) {
    unrollLaunchDeepCopy<char>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                               orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == short_type_index) {
    unrollLaunchDeepCopy<short>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                orig_offset, length, do_hdc, do_dhc, do_ddc, gpu, dest_bits,
                                orig_bits);
  }
  else if (ct == int_type_index) {
    unrollLaunchDeepCopy<int>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                              orig_offset, length, do_hdc, do_dhc, do_ddc, gpu, dest_bits,
                              orig_bits);
  }
  else if (ct == llint_type_index) {
    unrollLaunchDeepCopy<llint>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                orig_offset, length, do_hdc, do_dhc, do_ddc, gpu, dest_bits,
                                orig_bits);
  }
  else if (ct == uchar_type_index) {
    unrollLaunchDeepCopy<uchar>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == ushort_type_index) {
    unrollLaunchDeepCopy<ushort>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == uint_type_index) {
    unrollLaunchDeepCopy<uint>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                               orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == size_t_type_index) {
    unrollLaunchDeepCopy<size_t>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == ullint_type_index) {
    unrollLaunchDeepCopy<ullint>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == float2_type_index) {
    unrollLaunchDeepCopy<float2>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == float3_type_index) {
    unrollLaunchDeepCopy<float3>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == float4_type_index) {
    unrollLaunchDeepCopy<float4>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == double2_type_index) {
    unrollLaunchDeepCopy<double2>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                  orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == double3_type_index) {
    unrollLaunchDeepCopy<double3>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                  orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == double4_type_index) {
    unrollLaunchDeepCopy<double4>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                  orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == char2_type_index) {
    unrollLaunchDeepCopy<char2>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == char3_type_index) {
    unrollLaunchDeepCopy<char3>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == char4_type_index) {
    unrollLaunchDeepCopy<char4>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == short2_type_index) {
    unrollLaunchDeepCopy<short2>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == short3_type_index) {
    unrollLaunchDeepCopy<short3>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == short4_type_index) {
    unrollLaunchDeepCopy<short4>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == int2_type_index) {
    unrollLaunchDeepCopy<int2>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                               orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == int3_type_index) {
    unrollLaunchDeepCopy<int3>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                               orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == int4_type_index) {
    unrollLaunchDeepCopy<int4>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                               orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == longlong2_type_index) {
    unrollLaunchDeepCopy<llint2>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == longlong3_type_index) {
    unrollLaunchDeepCopy<llint3>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == longlong4_type_index) {
    unrollLaunchDeepCopy<llint4>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == uchar2_type_index) {
    unrollLaunchDeepCopy<uchar2>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == uchar3_type_index) {
    unrollLaunchDeepCopy<uchar3>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == uchar4_type_index) {
    unrollLaunchDeepCopy<uchar4>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == ushort2_type_index) {
    unrollLaunchDeepCopy<ushort2>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                  orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == ushort3_type_index) {
    unrollLaunchDeepCopy<ushort3>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                  orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == ushort4_type_index) {
    unrollLaunchDeepCopy<ushort4>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                  orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == uint2_type_index) {
    unrollLaunchDeepCopy<uint2>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == uint3_type_index) {
    unrollLaunchDeepCopy<uint3>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == uint4_type_index) {
    unrollLaunchDeepCopy<uint4>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == ulonglong2_type_index) {
    unrollLaunchDeepCopy<ullint2>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                  orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == ulonglong3_type_index) {
    unrollLaunchDeepCopy<ullint3>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                  orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
  else if (ct == ulonglong4_type_index) {
    unrollLaunchDeepCopy<ullint4>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                  orig_offset, length, do_hdc, do_dhc, do_ddc, gpu);
  }
}

//-------------------------------------------------------------------------------------------------
void launchDeepCopy(void* vdest_host, void* vdest_devc, const void* vorig_host,
                    const void* vorig_devc, const size_t dest_offset, const size_t orig_offset,
                    const size_t length, const size_t dest_ct, const size_t orig_ct,
                    const bool do_hdc, const bool do_dhc, const bool do_ddc, const GpuDetails &gpu,
                    const int dest_bits, const int orig_bits) {

  // Unroll along the original array type so that the templated unrollLaunchDeepCopy() (overloaded
  // variant including the type code for the destination arrays) can
  if (orig_ct == float_type_index) {
    unrollLaunchDeepCopy<float>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                orig_offset, length, dest_ct, do_hdc, do_dhc, do_ddc, gpu,
                                dest_bits, orig_bits);
  }
  else if (orig_ct == double_type_index) {
    unrollLaunchDeepCopy<double>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, dest_ct, do_hdc, do_dhc, do_ddc, gpu,
                                 dest_bits, orig_bits);
  }
  else if (orig_ct == char_type_index) {
    unrollLaunchDeepCopy<char>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                               orig_offset, length, dest_ct, do_hdc, do_dhc, do_ddc, gpu,
                               dest_bits, orig_bits);
  }
  else if (orig_ct == short_type_index) {
    unrollLaunchDeepCopy<short>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                orig_offset, length, dest_ct, do_hdc, do_dhc, do_ddc, gpu,
                                dest_bits, orig_bits);
  }
  else if (orig_ct == int_type_index) {
    unrollLaunchDeepCopy<int>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                              orig_offset, length, dest_ct, do_hdc, do_dhc, do_ddc, gpu, dest_bits,
                              orig_bits);
  }
  else if (orig_ct == llint_type_index) {
    unrollLaunchDeepCopy<llint>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                orig_offset, length, dest_ct, do_hdc, do_dhc, do_ddc, gpu,
                                dest_bits, orig_bits);
  }
  else if (orig_ct == uchar_type_index) {
    unrollLaunchDeepCopy<uchar>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                orig_offset, length, dest_ct, do_hdc, do_dhc, do_ddc, gpu,
                                dest_bits, orig_bits);
  }
  else if (orig_ct == ushort_type_index) {
    unrollLaunchDeepCopy<ushort>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, dest_ct, do_hdc, do_dhc, do_ddc, gpu,
                                 dest_bits, orig_bits);
  }
  else if (orig_ct == uint_type_index) {
    unrollLaunchDeepCopy<uint>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                               orig_offset, length, dest_ct, do_hdc, do_dhc, do_ddc, gpu,
                               dest_bits, orig_bits);
  }
  else if (orig_ct == ullint_type_index) {
    unrollLaunchDeepCopy<ullint>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, dest_ct, do_hdc, do_dhc, do_ddc, gpu,
                                 dest_bits, orig_bits);
  }
  else if (orig_ct == bool_type_index) {
    unrollLaunchDeepCopy<bool>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                               orig_offset, length, dest_ct, do_hdc, do_dhc, do_ddc, gpu,
                               dest_bits, orig_bits);
  }
  else if (orig_ct == size_t_type_index) {
    unrollLaunchDeepCopy<size_t>(vdest_host, vdest_devc, vorig_host, vorig_devc, dest_offset,
                                 orig_offset, length, dest_ct, do_hdc, do_dhc, do_ddc, gpu,
                                 dest_bits, orig_bits);
  }
  else {
    rtErr("Only scalar data types are sanctioned for interconversion.", "launchDeepCopy");
  }
}

} // namespace card
} // namespace stormm
