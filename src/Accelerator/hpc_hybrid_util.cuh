// -*-c++-*-
#ifndef STORMM_HPC_HYBRID_CUH
#define STORMM_HPC_HYBRID_CUH

#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/stormm_vector_types.h"
#include "Numerics/split_fixed_precision.h"

namespace stormm {
namespace card {
  
#include "Numerics/accumulation.cui"

/// \brief Launch a kernel to perform elementwise copying between two Hybrid arrays, accomplishing
///        everything but host-to-host copying if possible.  In situations where two Hybrid arrays
///        represent a single set of values in split fixed-precision notation, the kernel will
///        perform the necessary conversion.  The conversions performed by the second, third, and
///        fourth templated overload of this function may not be bitwise-conservative in particular
///        cases of integer-to-integer conversions, for which the overloads supplied with two
///        fixed-precision bit counts should be called.
///
/// Overloaded:
///   - Directly assign one set of values to another set of values
///   - Assign one set of values as a real-valued or fixed-precision interpretation of the other
///   - Provide a single, real-valued conversion factor for multiplying original values when
///     transcribing them to the destination array(s)
///   - Provide a pair of integer bit counts for fixed-precision models when transcribing
///     information between integer types
///
/// \param dest_host    Pointer to the host-side array of the destination Hybrid object
/// \param dest_devc    Pointer to the device-side array of the destination Hybrid object
/// \param orig_host    Pointer to the host-side array of the original Hybrid object
/// \param orig_devc    Pointer to the device-side array of the original Hybrid object
/// \param dest_offset  Offset for writing indices to the destination array(s)
/// \param orig_offset  Offset for reading indices to the original array(s)
/// \param length       The total nuber of array elements, whether on the device or on the host, to
///                     transfer
/// \param conversion   The real-valued conversion factor for multiplying original values before
///                     placing them in the destination array
/// \param dest_bits    Integer-valued number of bits after the point in the destination array(s)
/// \param orig_bits    Integer-valued number of bits after the point in the original array(s)
/// \param do_hdc       Request host-to-device copying by setting this to TRUE
/// \param do_dhc       Request device-to-host copying by setting this to TRUE
/// \param do_ddc       Request device-to-device copying by setting this to TRUE
/// \{
template <typename Tdest, typename Torig> __global__ __launch_bounds__(medium_block_size, 2)
void kDeepCopy(Tdest* dest_host, Tdest* dest_devc, const Torig* orig_host, const Torig* orig_devc,
               const size_t dest_offset, const size_t orig_offset, const size_t length,
               const bool do_hdc, const bool do_dhc, const bool do_ddc) {
  const size_t stride = blockDim.x * gridDim.x;
  for (size_t tid = threadIdx.x + (blockIdx.x * blockDim.x); tid < length; tid += stride) {
    if (do_hdc) {
      dest_devc[dest_offset + tid] = orig_host[orig_offset + tid];
    }
    if (do_dhc) {
      dest_host[dest_offset + tid] = orig_devc[orig_offset + tid];
    }
    if (do_ddc) {
      dest_devc[dest_offset + tid] = orig_devc[orig_offset + tid];
    }
  }
}

template <typename Tdest, typename Torig> __global__ __launch_bounds__(medium_block_size, 2)
void kDeepCopy(Tdest* dest_host, Tdest* dest_devc, const Torig* orig_host, const Torig* orig_devc,
               const size_t dest_offset, const size_t orig_offset, const size_t length,
               const double conversion, const bool do_hdc, const bool do_dhc, const bool do_ddc) {
  const size_t stride = blockDim.x * gridDim.x;
  if (conversion > 1.0) {
    for (size_t tid = threadIdx.x + (blockIdx.x * blockDim.x); tid < length; tid += stride) {
      if (do_hdc) {
        dest_devc[dest_offset + tid] = __double2ll_rn((double)(orig_host[orig_offset + tid]) *
                                                    conversion);
      }
      if (do_dhc) {
        dest_host[dest_offset + tid] = __double2ll_rn((double)(orig_devc[orig_offset + tid]) *
                                                      conversion);
      }
      if (do_ddc) {
        dest_devc[dest_offset + tid] = __double2ll_rn((double)(orig_devc[orig_offset + tid]) *
                                                      conversion);
      }
    }
  }
  else {
    for (size_t tid = threadIdx.x + (blockIdx.x * blockDim.x); tid < length; tid += stride) {
      if (do_hdc) {
        dest_devc[dest_offset + tid] = (double)(orig_host[orig_offset + tid]) * conversion;
      }
      if (do_dhc) {
        dest_host[dest_offset + tid] = (double)(orig_devc[orig_offset + tid]) * conversion;
      }
      if (do_ddc) {
        dest_devc[dest_offset + tid] = (double)(orig_devc[orig_offset + tid]) * conversion;
      }
    }
  }
}

template <typename Tdest, typename Torig> __global__ __launch_bounds__(medium_block_size, 2)
void kDeepCopy(Tdest* dest_host, Tdest* dest_devc, const Torig* orig_host,
               const int* orig_host_ovrf, const Torig* orig_devc, const int* orig_devc_ovrf,
               const size_t dest_offset, const size_t orig_offset, const size_t length,
               const double conversion, const bool do_hdc, const bool do_dhc, const bool do_ddc) {
  const size_t stride = blockDim.x * gridDim.x;
  for (size_t tid = threadIdx.x + (blockIdx.x * blockDim.x); tid < length; tid += stride) {
    if (do_hdc) {
      dest_devc[dest_offset + tid] = int95ToDouble(orig_host[orig_offset + tid],
                                                   orig_host_ovrf[orig_offset + tid]) * conversion;
    }
    if (do_dhc) {
      dest_host[dest_offset + tid] = int95ToDouble(orig_devc[orig_offset + tid],
                                                   orig_devc_ovrf[orig_offset + tid]) * conversion;
    }
    if (do_ddc) {
      dest_devc[dest_offset + tid] = int95ToDouble(orig_devc[orig_offset + tid],
                                                   orig_devc_ovrf[orig_offset + tid]) * conversion;
    }
  }
}

template <typename Tdest, typename Torig> __global__ __launch_bounds__(medium_block_size, 2)
void kDeepCopy(Tdest* dest_host, int* dest_host_ovrf, Tdest* dest_devc, int* dest_devc_ovrf,
               const Torig* orig_host, const Torig* orig_devc, const size_t dest_offset,
               const size_t orig_offset, const size_t length, const double conversion,
               const bool do_hdc, const bool do_dhc, const bool do_ddc) {
  const size_t stride = blockDim.x * gridDim.x;
  for (size_t tid = threadIdx.x + (blockIdx.x * blockDim.x); tid < length; tid += stride) {
    if (do_hdc) {
      const int95_t val = doubleToInt95((double)(orig_host[orig_offset + tid]) * conversion);
      dest_devc[dest_offset + tid]      = val.x;
      dest_devc_ovrf[dest_offset + tid] = val.y;
    }
    if (do_dhc) {
      const int95_t val = doubleToInt95((double)(orig_devc[orig_offset + tid]) * conversion);
      dest_host[dest_offset + tid]      = val.x;
      dest_host_ovrf[dest_offset + tid] = val.y;
    }
    if (do_ddc) {
      const int95_t val = doubleToInt95((double)(orig_devc[orig_offset + tid]) * conversion);
      dest_devc[dest_offset + tid]      = val.x;
      dest_devc_ovrf[dest_offset + tid] = val.y;
    }
  }
}

__global__ __launch_bounds__(medium_block_size, 2)
void kDeepCopy(llint* dest_host, int* dest_host_ovrf, llint* dest_devc, int* dest_devc_ovrf,
               const llint* orig_host, const int* orig_host_ovrf, const llint* orig_devc,
               const int* orig_devc_int, const size_t dest_offset, const size_t orig_offset,
               const size_t length, const int dest_bits, const int orig_bits, const bool do_hdc,
               const bool do_dhc, const bool do_ddc);

template <typename Tdest>  __global__ __launch_bounds__(medium_block_size, 2)
void kDeepCopy(Tdest* dest_host, Tdest* dest_devc, const llint* orig_host,
               const int* orig_host_ovrf, const llint* orig_devc, const int* orig_devc_ovrf,
               const size_t dest_offset, const size_t orig_offset, const size_t length,
               const int dest_bits, const int orig_bits, const bool do_hdc, const bool do_dhc,
               const bool do_ddc) {
  const size_t stride = blockDim.x * gridDim.x;
  for (size_t tid = threadIdx.x + (blockIdx.x * blockDim.x); tid < length; tid += stride) {
    if (do_hdc) {
      const int95_t val = changeFPBits(orig_host[orig_offset + tid],
                                       orig_host_ovrf[orig_offset + tid], orig_bits, dest_bits);
      dest_devc[dest_offset + tid]      = val.x;
    }
    if (do_dhc) {
      const int95_t val = changeFPBits(orig_devc[orig_offset + tid],
                                       orig_devc_ovrf[orig_offset + tid], orig_bits, dest_bits);
      dest_host[dest_offset + tid]      = val.x;
    }
    if (do_ddc) {
      const int95_t val = changeFPBits(orig_devc[orig_offset + tid],
                                       orig_devc_ovrf[orig_offset + tid], orig_bits, dest_bits);
      dest_devc[dest_offset + tid]      = val.x;
    }
  }
}

template <typename Torig>  __global__ __launch_bounds__(medium_block_size, 2)
void kDeepCopy(llint* dest_host, int* dest_host_ovrf, llint* dest_devc, int* dest_devc_ovrf,
               const Torig* orig_host, const Torig* orig_devc, const size_t dest_offset,
               const size_t orig_offset, const size_t length, const int dest_bits,
               const int orig_bits, const bool do_hdc, const bool do_dhc, const bool do_ddc) {
  const size_t stride = blockDim.x * gridDim.x;
  for (size_t tid = threadIdx.x + (blockIdx.x * blockDim.x); tid < length; tid += stride) {
    if (do_hdc) {
      const int95_t val = changeFPBits(orig_host[orig_offset + tid], 0, orig_bits, dest_bits);
      dest_devc[dest_offset + tid]      = val.x;
      dest_devc_ovrf[dest_offset + tid] = val.y;
    }
    if (do_dhc) {
      const int95_t val = changeFPBits(orig_devc[orig_offset + tid], 0, orig_bits, dest_bits);
      dest_host[dest_offset + tid]      = val.x;
      dest_host_ovrf[dest_offset + tid] = val.y;
    }
    if (do_ddc) {
      const int95_t val = changeFPBits(orig_devc[orig_offset + tid], 0, orig_bits, dest_bits);
      dest_devc[dest_offset + tid]      = val.x;
      dest_devc_ovrf[dest_offset + tid] = val.y;
    }
  }
}

template <typename Tdest, typename Torig>  __global__ __launch_bounds__(medium_block_size, 2)
void kDeepCopy(Tdest* dest_host, Tdest* dest_devc,  const Torig* orig_host, const Torig* orig_devc,
               const size_t dest_offset, const size_t orig_offset, const size_t length,
               const int dest_bits, const int orig_bits, const bool do_hdc, const bool do_dhc,
               const bool do_ddc) {
  const size_t stride = blockDim.x * gridDim.x;
  const llint iconv = pow(2.0, abs(orig_bits - dest_bits));
  for (size_t tid = threadIdx.x + (blockIdx.x * blockDim.x); tid < length; tid += stride) {
    if (do_hdc) {
      const llint val = orig_host[orig_offset + tid];
      if (orig_bits > dest_bits) {
        dest_devc[dest_offset + tid] = val / iconv;
      }
      else if (orig_bits < dest_bits) {
        dest_devc[dest_offset + tid] = val * iconv;
      }
      else {
        dest_devc[dest_offset + tid] = val;
      }
    }
    if (do_dhc) {
      const llint val = orig_devc[orig_offset + tid];
      if (orig_bits > dest_bits) {
        dest_host[dest_offset + tid] = val / iconv;
      }
      else if (orig_bits < dest_bits) {
        dest_host[dest_offset + tid] = val * iconv;
      }
      else {
        dest_host[dest_offset + tid] = val;
      }
    }
    if (do_ddc) {
      const llint val = orig_devc[orig_offset + tid];
      if (orig_bits > dest_bits) {
        dest_devc[dest_offset + tid] = val / iconv;
      }
      else if (orig_bits < dest_bits) {
        dest_devc[dest_offset + tid] = val * iconv;
      }
      else {
        dest_devc[dest_offset + tid] = val;
      }
    }
  }
}
/// \}

/// \brief Unroll the launchDeepCopy kernel launch by establishing one of the two templated types.
///        Descriptions of input parameters follow from kDeepCopy(), above, with analogous names:
///
/// \param vdest_host  Pointer to the void-casted host-side array of the destination Hybrid
///                    object's host-side data
/// \param vdest_devc  Pointer to the void-casted host-side array of the destination Hybrid
///                    object's device-side data
/// \{
template <typename T>
void unrollLaunchDeepCopy(void* vdest_host, void* vdest_devc, const void* vorig_host,
                          const void* vorig_devc, const size_t dest_offset,
                          const size_t orig_offset, const size_t length, const bool do_hdc,
                          const bool do_dhc, const bool do_ddc, const GpuDetails &gpu) {
  T* dest_host = reinterpret_cast<T*>(vdest_host);
  T* dest_devc = reinterpret_cast<T*>(vdest_devc);
  const T* orig_host = reinterpret_cast<const T*>(vorig_host);
  const T* orig_devc = reinterpret_cast<const T*>(vorig_devc);
  const int nblk = 2 * gpu.getSMPCount();
  const int nthr = medium_block_size;
  kDeepCopy<T, T><<<nblk, nthr>>>(dest_host, dest_devc, orig_host, orig_devc, dest_offset,
                                  orig_offset, length, do_hdc, do_dhc, do_ddc);
}

template <typename T>
void unrollLaunchDeepCopy(void* vdest_host, void* vdest_devc, const void* vorig_host,
                          const void* vorig_devc, const size_t dest_offset,
                          const size_t orig_offset, const size_t length, const bool do_hdc,
                          const bool do_dhc, const bool do_ddc, const GpuDetails &gpu,
                          const int dest_bits, const int orig_bits) {
  T* dest_host = reinterpret_cast<T*>(vdest_host);
  T* dest_devc = reinterpret_cast<T*>(vdest_devc);
  const T* orig_host = reinterpret_cast<const T*>(vorig_host);
  const T* orig_devc = reinterpret_cast<const T*>(vorig_devc);
  const int nblk = 2 * gpu.getSMPCount();
  const int nthr = medium_block_size;
  kDeepCopy<T, T><<<nblk, nthr>>>(dest_host, dest_devc, orig_host, orig_devc, dest_offset,
                                  orig_offset, length, dest_bits, orig_bits, do_hdc, do_dhc,
                                  do_ddc);
}

template <typename Torig>
void unrollLaunchDeepCopy(void* vdest_host, void* vdest_devc, const void* vorig_host,
                          const void* vorig_devc, const size_t dest_offset,
                          const size_t orig_offset, const size_t length, const size_t ct_dest,
                          const bool do_hdc, const bool do_dhc, const bool do_ddc,
                          const GpuDetails &gpu, const int dest_bits, const int orig_bits) {
  const Torig* orig_host = reinterpret_cast<const Torig*>(vorig_host);
  const Torig* orig_devc = reinterpret_cast<const Torig*>(vorig_devc);
  const int nblk = 2 * gpu.getSMPCount();
  const int nthr = medium_block_size;
  if (orig_bits == 0) {
    if (dest_bits == 0) {
      if (ct_dest == float_type_index) {
        float* dest_host = reinterpret_cast<float*>(vdest_host);
        float* dest_devc = reinterpret_cast<float*>(vdest_devc);
        kDeepCopy<float, Torig><<<nblk, nthr>>>(dest_host, dest_devc, orig_host, orig_devc,
                                                dest_offset, orig_offset, length, do_hdc, do_dhc,
                                                do_ddc);
      }
      else if (ct_dest == double_type_index) {
        double* dest_host = reinterpret_cast<double*>(vdest_host);
        double* dest_devc = reinterpret_cast<double*>(vdest_devc);
        kDeepCopy<double, Torig><<<nblk, nthr>>>(dest_host, dest_devc, orig_host, orig_devc,
                                                 dest_offset, orig_offset, length, do_hdc, do_dhc,
                                                 do_ddc);
      }
      else {
        rtErr("Real-valued numbers are only sanctioned for conversion to other real-valued "
              "numbers unless a nontrivial scaling factor is applied.", "unrollDeepCopy");
      }
    }
    else {
      const double factor = pow(2.0, dest_bits);
      if (ct_dest == short_type_index) {
        short* dest_host = reinterpret_cast<short*>(vdest_host);
        short* dest_devc = reinterpret_cast<short*>(vdest_devc);
        kDeepCopy<short, Torig><<<nblk, nthr>>>(dest_host, dest_devc, orig_host, orig_devc,
                                                dest_offset, orig_offset, length, factor, do_hdc,
                                                do_dhc, do_ddc);
      }
      else if (ct_dest == int_type_index) {
        int* dest_host = reinterpret_cast<int*>(vdest_host);
        int* dest_devc = reinterpret_cast<int*>(vdest_devc);
        kDeepCopy<int, Torig><<<nblk, nthr>>>(dest_host, dest_devc, orig_host, orig_devc,
                                              dest_offset, orig_offset, length, factor, do_hdc,
                                              do_dhc, do_ddc);
      }
      else if (ct_dest == llint_type_index) {
        llint* dest_host = reinterpret_cast<llint*>(vdest_host);
        llint* dest_devc = reinterpret_cast<llint*>(vdest_devc);
        kDeepCopy<llint, Torig><<<nblk, nthr>>>(dest_host, dest_devc, orig_host, orig_devc,
                                                dest_offset, orig_offset, length, factor, do_hdc,
                                                do_dhc, do_ddc);
      }
      else {
        rtErr("Only signed integral scalar types are sanctioned for real-valued to "
              "fixed-precision conversions.", "unrollLaunchDeepCopy");
      }
    }
  }
  else {
    if (dest_bits == 0) {
      const double factor = pow(2.0, -orig_bits);
      if (ct_dest == float_type_index) {
        float* dest_host = reinterpret_cast<float*>(vdest_host);
        float* dest_devc = reinterpret_cast<float*>(vdest_devc);
        kDeepCopy<float, Torig><<<nblk, nthr>>>(dest_host, dest_devc, orig_host, orig_devc,
                                                dest_offset, orig_offset, length, factor, do_hdc,
                                                do_dhc, do_ddc);
      }
      else if (ct_dest == double_type_index) {
        double* dest_host = reinterpret_cast<double*>(vdest_host);
        double* dest_devc = reinterpret_cast<double*>(vdest_devc);
        kDeepCopy<double, Torig><<<nblk, nthr>>>(dest_host, dest_devc, orig_host, orig_devc,
                                                 dest_offset, orig_offset, length, factor, do_hdc,
                                                 do_dhc, do_ddc);
      }
      else {
        rtErr("Signed integral scalar types are sanctioned for fixed-precision to real-valued "
              "conversions only with real-valued data types.", "unrollLaunchDeepCopy");
      }
    }
    else {
      if (ct_dest == short_type_index) {
        short* dest_host = reinterpret_cast<short*>(vdest_host);
        short* dest_devc = reinterpret_cast<short*>(vdest_devc);
        kDeepCopy<short, Torig><<<nblk, nthr>>>(dest_host, dest_devc, orig_host, orig_devc,
                                                dest_offset, orig_offset, length, dest_bits,
                                                orig_bits, do_hdc, do_dhc, do_ddc);
      }
      else if (ct_dest == int_type_index) {
        int* dest_host = reinterpret_cast<int*>(vdest_host);
        int* dest_devc = reinterpret_cast<int*>(vdest_devc);
        kDeepCopy<int, Torig><<<nblk, nthr>>>(dest_host, dest_devc, orig_host, orig_devc,
                                              dest_offset, orig_offset, length, dest_bits,
                                              orig_bits, do_hdc, do_dhc, do_ddc);
      }
      else if (ct_dest == llint_type_index) {
        llint* dest_host = reinterpret_cast<llint*>(vdest_host);
        llint* dest_devc = reinterpret_cast<llint*>(vdest_devc);
        kDeepCopy<llint, Torig><<<nblk, nthr>>>(dest_host, dest_devc, orig_host, orig_devc,
                                                dest_offset, orig_offset, length, dest_bits,
                                                orig_bits, do_hdc, do_dhc, do_ddc);
      }
      else {
        rtErr("Only signed integral scalar types are sanctioned for fixed-precision scaling "
              "factor conversions.", "unrollLaunchDeepCopy");
      }
    }
  }
}
/// \}

} // namespace card
} // namespace stormm

#endif
