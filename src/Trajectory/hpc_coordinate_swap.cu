// -*-c++-*-
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "hpc_coordinate_swap.h"
#include "hpc_coordinate_swap.cuh"

namespace stormm {
namespace trajectory {

using constants::PrecisionModel;
using numerics::globalpos_scale_nonoverflow_bits;
using numerics::velocity_scale_nonoverflow_bits;
using numerics::force_scale_nonoverflow_bits;
  
/// \brief Swap a set of Cartesian X, Y, or Z coordinates.  Following the related, overloaded
///        __device__ function copyCoordinateSet() in hpc_coordinate_copy.cuh, this function
///        increments a counter in performing its work, to allow warps of any given block to
///        make simultaneous progress over each coordinate axis.
///
/// \param first             The first array containing coordinates to swap
/// \param first_ovrf        Overflow bits for the first array
/// \param first_start_idx   Starting index at which to begin swapping from the first array
/// \param second            The second array containing coordinates to swap
/// \param second_ovrf       Overflow bits for the first array
/// \param second_start_idx  Starting index at which to begin swapping from the second array
/// \param count             The number of particle coordinates to copy
/// \param pos               Internal counter passed in from the calling function
/// \param iter              The number of passes through this array, each incrementing the counter
///                          behind pos
/// \param stride            The padded form of count (pre-computed and passed in for convenience
///                          and probably reduced register pressure)
/// \param advance           The amount to increment an internal variant of the pos counter (this
///                          can be the width of a warp, the width of a block, or the width of the
///                          entire kernel launch grid depending on the context in which this copy
///                          routine is called)
__device__ __forceinline__
size_t swapCoordinateSet(llint* first, int* first_ovrf, const size_t first_start_idx,
                         llint* second, int* second_ovrf, const size_t second_start_idx,
                         const int count, const size_t pos, const size_t iter, const size_t stride,
                         const size_t advance) {
  const size_t pos_limit = (iter + 1) * stride;
  size_t ipos = pos;
  while (ipos < pos_limit) {
    const int rel_pos = ipos - (iter * stride);
    if (rel_pos < count) {
      const size_t a_idx = first_start_idx + rel_pos;
      const size_t b_idx = second_start_idx + rel_pos;
      const llint tmp_ll = __ldcv(&first[a_idx]);
      __stwt(&first[a_idx], __ldcv(&second[b_idx]));
      __stwt(&second[b_idx], tmp_ll);
      const int tmp = __ldcv(&first_ovrf[a_idx]);
      __stwt(&first_ovrf[a_idx], __ldcv(&second_ovrf[b_idx]));
      __stwt(&second_ovrf[b_idx], tmp);
    }
    ipos += advance;
  }
  return ipos;
}

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__
void kSwapBoxes(double* first_umat, double* first_invu, double* first_bdim, double* second_umat,
                double* second_invu, double* second_bdim, const size_t frm_first,
                const size_t frm_second) {
  const int warp_idx = threadIdx.x >> warp_bits;
  if (warp_idx == 0) {
    const size_t first_xfrm_offset  = devcRoundUp(9, warp_size_zu) * frm_first;
    const size_t second_xfrm_offset = devcRoundUp(9, warp_size_zu) * frm_second;
    size_t pos = threadIdx.x;
    while (pos < 9) {
      const size_t a_idx = first_xfrm_offset + pos;
      const size_t b_idx = second_xfrm_offset + pos;
      const double tmp = __ldcv(&first_umat[a_idx]);
      __stwt(&first_umat[a_idx], __ldcv(&second_umat[b_idx]));
      __stwt(&second_umat[b_idx], tmp);
      pos += warp_size_zu;
    }
  }
  else if (warp_idx == 1) {
    const size_t first_xfrm_offset  = devcRoundUp(9, warp_size_zu) * frm_first;
    const size_t second_xfrm_offset = devcRoundUp(9, warp_size_zu) * frm_second;
    size_t pos = (threadIdx.x & warp_bits_mask_int);
    while (pos < 9) {
      const size_t a_idx = first_xfrm_offset + pos;
      const size_t b_idx = second_xfrm_offset + pos;
      const double tmp = __ldcv(&first_invu[a_idx]);
      __stwt(&first_invu[a_idx], __ldcv(&second_invu[b_idx]));
      __stwt(&second_invu[b_idx], tmp);
      pos += warp_size_zu;
    }
  }
  else if (warp_idx == 2) {
    const size_t first_xfrm_offset  = devcRoundUp(6, warp_size_zu) * frm_first;
    const size_t second_xfrm_offset = devcRoundUp(6, warp_size_zu) * frm_second;
    size_t pos = (threadIdx.x & warp_bits_mask_int);
    while (pos < 6) {
      const size_t a_idx = first_xfrm_offset + pos;
      const size_t b_idx = second_xfrm_offset + pos;
      const double tmp = __ldcv(&first_bdim[a_idx]);
      __stwt(&first_bdim[a_idx], __ldcv(&second_bdim[b_idx]));
      __stwt(&second_bdim[b_idx], tmp);
      pos += warp_size_zu;
    }
  }
}
  
/// \brief Swap coordinates in between two coordinate synthesis objects.
///
/// Overloaded:
///   - Swap between PhaseSpaceSynthesis objects or Condensates
///   - Swap one system or multiple systems
///
/// \param first       The first coordinate synthesis object
/// \param second      The second coordinate synthesis object
/// \param frm_first   The system index within the first object to swap
/// \param frm_first   The system index within the second object to swap
/// \param frm_pairs   Array of system index pairs to swap.  Any given launchSwapCoordinates
///                    function will either have this array with its length pair_count, or both of
///                    frm_first and frm_second.
/// \param pair_count  The trusted length of frm_pairs  
/// \{
__global__ void __launch_bounds__(small_block_size, 4)
kSwapCoordinates(PsSynthesisWriter first, const size_t frm_first, PsSynthesisWriter second,
                 const size_t frm_second) {
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);

  // The atom count, as well as fixed-precision bit settings, have been verified to correspond
  // between the two objects.
  const int natom = first.atom_counts[frm_first];
  const size_t stride = devcRoundUp(natom, warp_size_int);
  const size_t advance = (blockDim.x * gridDim.x);
  const size_t first_offset = first.atom_starts[frm_first];
  const size_t second_offset = second.atom_starts[frm_second];

  // If the precision of positions, velocities, or forces is such that the overflow would not be
  // needed, save the communication bandwidth by calling the templated version of the swap
  // function.  Otherwise, use the full-precision int95_t swap function.
  if (first.gpos_bits <= globalpos_scale_nonoverflow_bits) {
    pos = swapCoordinateSet(first.xcrd, first_offset, second.xcrd, second_offset, natom, pos, 0,
                            stride, advance);
    pos = swapCoordinateSet(first.ycrd, first_offset, second.ycrd, second_offset, natom, pos, 1,
                            stride, advance);
    pos = swapCoordinateSet(first.zcrd, first_offset, second.zcrd, second_offset, natom, pos, 2,
                            stride, advance);
  }
  else {
    pos = swapCoordinateSet(first.xcrd, first.xcrd_ovrf, first_offset, second.xcrd,
                            second.xcrd_ovrf, second_offset, natom, pos, 0, stride, advance);
    pos = swapCoordinateSet(first.ycrd, first.ycrd_ovrf, first_offset, second.ycrd,
                            second.ycrd_ovrf, second_offset, natom, pos, 1, stride, advance);
    pos = swapCoordinateSet(first.zcrd, first.zcrd_ovrf, first_offset, second.zcrd,
                            second.zcrd_ovrf, second_offset, natom, pos, 2, stride, advance);
  }
  if (first.vel_bits <= velocity_scale_nonoverflow_bits) {
    pos = swapCoordinateSet(first.xvel, first_offset, second.xvel, second_offset, natom, pos, 3,
                            stride, advance);
    pos = swapCoordinateSet(first.yvel, first_offset, second.yvel, second_offset, natom, pos, 4,
                            stride, advance);
    pos = swapCoordinateSet(first.zvel, first_offset, second.zvel, second_offset, natom, pos, 5,
                            stride, advance);
  }
  else {
    pos = swapCoordinateSet(first.xvel, first.xvel_ovrf, first_offset, second.xvel,
                            second.xvel_ovrf, second_offset, natom, pos, 3, stride, advance);
    pos = swapCoordinateSet(first.yvel, first.yvel_ovrf, first_offset, second.yvel,
                            second.yvel_ovrf, second_offset, natom, pos, 4, stride, advance);
    pos = swapCoordinateSet(first.zvel, first.zvel_ovrf, first_offset, second.zvel,
                            second.zvel_ovrf, second_offset, natom, pos, 5, stride, advance);
  }
  if (first.frc_bits <= force_scale_nonoverflow_bits) {
    pos = swapCoordinateSet(first.xfrc, first_offset, second.xfrc, second_offset, natom, pos, 6,
                            stride, advance);
    pos = swapCoordinateSet(first.yfrc, first_offset, second.yfrc, second_offset, natom, pos, 7,
                            stride, advance);
    pos = swapCoordinateSet(first.zfrc, first_offset, second.zfrc, second_offset, natom, pos, 8,
                            stride, advance);
  }
  else {
    pos = swapCoordinateSet(first.xfrc, first.xfrc_ovrf, first_offset, second.xfrc,
                            second.xfrc_ovrf, second_offset, natom, pos, 6, stride, advance);
    pos = swapCoordinateSet(first.yfrc, first.yfrc_ovrf, first_offset, second.yfrc,
                            second.yfrc_ovrf, second_offset, natom, pos, 7, stride, advance);
    pos = swapCoordinateSet(first.zfrc, first.zfrc_ovrf, first_offset, second.zfrc,
                            second.zfrc_ovrf, second_offset, natom, pos, 8, stride, advance);
  }
  if (blockIdx.x == 0) {
    kSwapBoxes(first.umat, first.invu, first.boxdims, second.umat, second.invu, second.boxdims,
               frm_first, frm_second);
  }
}

__global__ void __launch_bounds__(small_block_size, 4)
kSwapCoordinates(PsSynthesisWriter first, PsSynthesisWriter second, const int2* frm_pairs,
                 const size_t pair_count) {
  const int blocks_per_frm = ((3 * gridDim.x) + pair_count - 1) / pair_count;
  int work_idx = blockIdx.x;
  while (work_idx < pair_count * blocks_per_frm) {
    const int pair_idx = work_idx / blocks_per_frm;
    const int2 swap_x_to_y = __ldca(&frm_pairs[pair_idx]);

    // Skip swap requests for invalid systems.  The requests will be checked on the CPU, where an
    // error will be raised.  The precision models of each synthesis will have been checked prior
    // to launching this kernel, but the list of all system indices and their respective atom
    // counts will be checked on the CPU after launching the kernel, to avoid making the GPU wait
    // on what is a potentially longer process and instead mask it with the kernel processing
    // time.
    if (swap_x_to_y.x < 0 || swap_x_to_y.x >= first.system_count ||
        swap_x_to_y.y < 0 || swap_x_to_y.y >= second.system_count ||
        __ldca(&first.atom_counts[swap_x_to_y.x]) != __ldca(&second.atom_counts[swap_x_to_y.y])) {
      work_idx += gridDim.x;
      continue;
    }

    // As in the single-system case above, avoid swapping overflow bits if they will not be needed.
    const int natom = first.atom_counts[swap_x_to_y.x];
    const size_t padded_natom = devcRoundUp(natom, warp_size_int);
    const size_t first_offset = first.atom_starts[swap_x_to_y.x];
    const size_t second_offset = second.atom_starts[swap_x_to_y.y];
    const int wu_in_pair = (work_idx % blocks_per_frm);
    size_t pos = threadIdx.x + (wu_in_pair * blockDim.x);
    if (first.gpos_bits <= globalpos_scale_nonoverflow_bits) {
      pos = swapCoordinateSet(first.xcrd, first_offset, second.xcrd, second_offset, natom, pos, 0,
                              padded_natom, blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.ycrd, first_offset, second.ycrd, second_offset, natom, pos, 1,
                              padded_natom, blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.zcrd, first_offset, second.zcrd, second_offset, natom, pos, 2,
                              padded_natom, blockDim.x * blocks_per_frm);
    }
    else {
      pos = swapCoordinateSet(first.xcrd, first.xcrd_ovrf, first_offset, second.xcrd,
                              second.xcrd_ovrf, second_offset, natom, pos, 0, padded_natom,
                              blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.ycrd, first.ycrd_ovrf, first_offset, second.ycrd,
                              second.ycrd_ovrf, second_offset, natom, pos, 1, padded_natom,
                              blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.zcrd, first.zcrd_ovrf, first_offset, second.zcrd,
                              second.zcrd_ovrf, second_offset, natom, pos, 2, padded_natom,
                              blockDim.x * blocks_per_frm);
    }
    if (first.vel_bits <= velocity_scale_nonoverflow_bits) {
      pos = swapCoordinateSet(first.xvel, first_offset, second.xvel, second_offset, natom, pos, 3,
                              padded_natom, blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.yvel, first_offset, second.yvel, second_offset, natom, pos, 4,
                              padded_natom, blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.zvel, first_offset, second.zvel, second_offset, natom, pos, 5,
                              padded_natom, blockDim.x * blocks_per_frm);
    }
    else {
      pos = swapCoordinateSet(first.xvel, first.xvel_ovrf, first_offset, second.xvel,
                              second.xvel_ovrf, second_offset, natom, pos, 3, padded_natom,
                              blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.yvel, first.yvel_ovrf, first_offset, second.yvel,
                              second.yvel_ovrf, second_offset, natom, pos, 4, padded_natom,
                              blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.zvel, first.zvel_ovrf, first_offset, second.zvel,
                              second.zvel_ovrf, second_offset, natom, pos, 5, padded_natom,
                              blockDim.x * blocks_per_frm);
    }
    if (first.frc_bits <= force_scale_nonoverflow_bits) {
      pos = swapCoordinateSet(first.xfrc, first_offset, second.xfrc, second_offset, natom, pos, 6,
                              padded_natom, blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.yfrc, first_offset, second.yfrc, second_offset, natom, pos, 7,
                              padded_natom, blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.zfrc, first_offset, second.zfrc, second_offset, natom, pos, 8,
                              padded_natom, blockDim.x * blocks_per_frm);
    }
    else {
      pos = swapCoordinateSet(first.xfrc, first.xfrc_ovrf, first_offset, second.xfrc,
                              second.xfrc_ovrf, second_offset, natom, pos, 6, padded_natom,
                              blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.yfrc, first.yfrc_ovrf, first_offset, second.yfrc,
                              second.yfrc_ovrf, second_offset, natom, pos, 7, padded_natom,
                              blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.zfrc, first.zfrc_ovrf, first_offset, second.zfrc,
                              second.zfrc_ovrf, second_offset, natom, pos, 8, padded_natom,
                              blockDim.x * blocks_per_frm);
    }
    if (wu_in_pair == 0) {
      kSwapBoxes(first.umat, first.invu, first.boxdims, second.umat, second.invu, second.boxdims,
                 swap_x_to_y.x, swap_x_to_y.y);
    }
    
    // Advance to the next work unit, perhaps a new system to swap
    work_idx += gridDim.x;
  }
}

__global__ void __launch_bounds__(small_block_size, 4)
kSwapCoordinates(CondensateWriter first, const size_t frm_first, CondensateWriter second,
                 const size_t frm_second) {

  // The atom count, as well as fixed-precision bit settings, have been verified to correspond
  // between the two objects.
  const int natom = first.atom_counts[frm_first];
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t stride = devcRoundUp(natom, warp_size_int);
  const size_t advance = (blockDim.x * gridDim.x);
  const size_t first_offset = first.atom_starts[frm_first];
  const size_t second_offset = second.atom_starts[frm_second];

  // If the precision of positions, velocities, or forces is such that the overflow would not be
  // needed, save the communication bandwidth by calling the templated version of the swap
  // function.  Otherwise, use the full-precision int95_t swap function.
  switch (first.mode) {
  case PrecisionModel::DOUBLE:
    pos = swapCoordinateSet(first.xcrd, first_offset, second.xcrd, second_offset, natom, pos, 0,
                            stride, advance);
    pos = swapCoordinateSet(first.ycrd, first_offset, second.ycrd, second_offset, natom, pos, 1,
                            stride, advance);
    pos = swapCoordinateSet(first.zcrd, first_offset, second.zcrd, second_offset, natom, pos, 2,
                            stride, advance);
    break;
  case PrecisionModel::SINGLE:
    pos = swapCoordinateSet(first.xcrd_sp, first_offset, second.xcrd_sp, second_offset, natom,
                            pos, 0, stride, advance);
    pos = swapCoordinateSet(first.ycrd_sp, first_offset, second.ycrd_sp, second_offset, natom,
                            pos, 1, stride, advance);
    pos = swapCoordinateSet(first.zcrd_sp, first_offset, second.zcrd_sp, second_offset, natom,
                            pos, 2, stride, advance);
    break;
  }
  if (blockIdx.x == 0) {
    kSwapBoxes(first.umat, first.invu, first.boxdims, second.umat, second.invu, second.boxdims,
               frm_first, frm_second);
  }
}

__global__ void __launch_bounds__(small_block_size, 4)
kSwapCoordinates(CondensateWriter first, CondensateWriter second, const int2* frm_pairs,
                 const size_t pair_count) {
  const int blocks_per_frm = (gridDim.x + pair_count - 1) / pair_count;
  int work_idx = blockIdx.x;
  while (work_idx < pair_count * blocks_per_frm) {
    const int pair_idx = work_idx / blocks_per_frm;
    const int2 swap_x_to_y = __ldca(&frm_pairs[pair_idx]);

    // Skip swap requests for invalid systems.  The requests will be checked on the CPU, where an
    // error will be raised.  The precision modes of each Condensate object will have been checked
    // prior to launching this kernel.
    if (swap_x_to_y.x < 0 || swap_x_to_y.x >= first.system_count ||
        swap_x_to_y.y < 0 || swap_x_to_y.y >= second.system_count ||
        __ldca(&first.atom_counts[swap_x_to_y.x]) != __ldca(&second.atom_counts[swap_x_to_y.y])) {
      work_idx += gridDim.x;
      continue;
    }

    // Compute the atom counts and offsets, then loop over X, Y, and Z coordinates
    const int natom = first.atom_counts[swap_x_to_y.x];
    const size_t padded_natom = devcRoundUp(natom, warp_size_int);
    const size_t first_offset = first.atom_starts[swap_x_to_y.x];
    const size_t second_offset = second.atom_starts[swap_x_to_y.y];
    const int wu_in_pair = (work_idx % blocks_per_frm);
    size_t pos = threadIdx.x + (wu_in_pair * blockDim.x);
    switch (first.mode) {
    case PrecisionModel::DOUBLE:
      pos = swapCoordinateSet(first.xcrd, first_offset, second.xcrd, second_offset, natom, pos, 0,
                              padded_natom, blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.ycrd, first_offset, second.ycrd, second_offset, natom, pos, 1,
                              padded_natom, blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.zcrd, first_offset, second.zcrd, second_offset, natom, pos, 2,
                              padded_natom, blockDim.x * blocks_per_frm);
      break;
    case PrecisionModel::SINGLE:
      pos = swapCoordinateSet(first.xcrd_sp, first_offset, second.xcrd_sp, second_offset, natom,
                              pos, 0, padded_natom, blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.ycrd_sp, first_offset, second.ycrd_sp, second_offset, natom,
                              pos, 1, padded_natom, blockDim.x * blocks_per_frm);
      pos = swapCoordinateSet(first.zcrd_sp, first_offset, second.zcrd_sp, second_offset, natom,
                              pos, 2, padded_natom, blockDim.x * blocks_per_frm);
      break;
    }
    if (wu_in_pair == 0) {
      kSwapBoxes(first.umat, first.invu, first.boxdims, second.umat, second.invu, second.boxdims,
                 swap_x_to_y.x, swap_x_to_y.y);
    }

    // Advance to the next work unit in the swap
    work_idx += gridDim.x;
  }
}
/// \}

//-------------------------------------------------------------------------------------------------
void launchSwapCoordinates(CoordinateSeriesWriter<void> *v_first, const size_t frm_first,
                           CoordinateSeriesWriter<void> *v_second, const size_t frm_second,
                           const size_t ct, const GpuDetails &gpu) {
  if (ct == double_type_index) {
    unrollSwapCoordLaunch<double>(v_first, frm_first, v_second, frm_second, gpu);
  }
  else if (ct == float_type_index) {
    unrollSwapCoordLaunch<float>(v_first, frm_first, v_second, frm_second, gpu);
  }
  else if (ct == llint_type_index) {
    unrollSwapCoordLaunch<llint>(v_first, frm_first, v_second, frm_second, gpu);
  }
  else if (ct == int_type_index) {
    unrollSwapCoordLaunch<int>(v_first, frm_first, v_second, frm_second, gpu);
  }
  else if (ct == short_type_index) {
    unrollSwapCoordLaunch<short int>(v_first, frm_first, v_second, frm_second, gpu);
  }
}

//-------------------------------------------------------------------------------------------------
void launchSwapCoordinates(CoordinateSeriesWriter<void> *v_first,
                           CoordinateSeriesWriter<void> *v_second, const int2* frm_pairs,
                           const int pair_count, const size_t ct, const GpuDetails &gpu) {
  if (ct == double_type_index) {
    unrollSwapCoordLaunch<double>(v_first, v_second, frm_pairs, pair_count, gpu);
  }
  else if (ct == float_type_index) {
    unrollSwapCoordLaunch<float>(v_first, v_second, frm_pairs, pair_count, gpu);
  }
  else if (ct == llint_type_index) {
    unrollSwapCoordLaunch<llint>(v_first, v_second, frm_pairs, pair_count, gpu);
  }
  else if (ct == int_type_index) {
    unrollSwapCoordLaunch<int>(v_first, v_second, frm_pairs, pair_count, gpu);
  }
  else if (ct == short_type_index) {
    unrollSwapCoordLaunch<short int>(v_first, v_second, frm_pairs, pair_count, gpu);
  }
}

//-------------------------------------------------------------------------------------------------
void launchSwapCoordinates(PsSynthesisWriter *first, const size_t frm_first,
                           PsSynthesisWriter *second, const size_t frm_second,
                           const GpuDetails &gpu) {  
  kSwapCoordinates<<<4 * gpu.getSMPCount(), small_block_size>>>(*first, frm_first, *second,
                                                                frm_second);
}

//-------------------------------------------------------------------------------------------------
void launchSwapCoordinates(PsSynthesisWriter *first, PsSynthesisWriter *second,
                           const int2* frm_pairs, const int pair_count, const GpuDetails &gpu) {
  kSwapCoordinates<<<4 * gpu.getSMPCount(), small_block_size>>>(*first, *second, frm_pairs,
                                                                pair_count);
}

//-------------------------------------------------------------------------------------------------
void launchSwapCoordinates(CondensateWriter *first, const size_t frm_first,
                           CondensateWriter *second, const size_t frm_second,
                           const GpuDetails &gpu) {
  kSwapCoordinates<<<4 *gpu.getSMPCount(), small_block_size>>>(*first, frm_first, *second,
                                                               frm_second);
}

//-------------------------------------------------------------------------------------------------
void launchSwapCoordinates(CondensateWriter *first, CondensateWriter *second,
                           const int2* frm_pairs, const int pair_count, const GpuDetails &gpu) {
  kSwapCoordinates<<<4 * gpu.getSMPCount(), small_block_size>>>(*first, *second, frm_pairs,
                                                                pair_count);
}                           

} // namespace trajectory
} // namespace stormm
