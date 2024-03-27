// -*-c++-*-
#ifndef STORMM_HPC_COORDINATE_SWAP_CUH
#define STORMM_HPC_COORDINATE_SWAP_CUH

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "coordinate_series.h"

namespace stormm {
namespace trajectory {

using card::GpuDetails;

#include "Math/rounding.cui"

/// \brief Swap a set of Cartesian X, Y, or Z coordinates.  Following the related, overloaded
///        __device__ function copyCoordinateSet() in hpc_coordinate_copy.cuh, these functions
///        increment a counter as they perform their work, to allow warps of any given block to
///        make simultaneous progress over each coordinate axis.  This is a templated overload of
///        the variant swapping int95_t data.
///
/// \param first             The first array containing coordinates to swap
/// \param first_start_idx   Starting index at which to begin swapping from the first array
/// \param second            The second array containing coordinates to swap
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
template <typename T> __device__ __forceinline__
size_t swapCoordinateSet(T* first, const size_t first_start_idx, T* second,
                         const size_t second_start_idx, const int count, const size_t pos,
                         const size_t iter, const size_t stride, const size_t advance) {
  const size_t pos_limit = (iter + 1) * stride;
  size_t ipos = pos;
  while (ipos < pos_limit) {
    const int rel_pos = ipos - (iter * stride);
    if (rel_pos < count) {
      const size_t a_idx = first_start_idx + rel_pos;
      const size_t b_idx = second_start_idx + rel_pos;
      const T tmp = __ldcv(&first[a_idx]);
      __stwt(&first[a_idx], __ldcv(&second[b_idx]));
      __stwt(&second[b_idx], tmp);
    }
    ipos += advance;
  }
  return ipos;
}

/// \brief Swap box information (transformation matrices into and out of box space as well as box
///        dimensions).  This function assumes that each thread block contains at least three
///        warps.
///
/// \param first_umat   Box space transformation matrices for the first collection of systems,
///                     taking Cartesian coordinates into unit cell fraction space
/// \param first_invu   Inverse transformation matrices for the first collections of systems,
///                     taking unit cell fractional coordinates into Cartesian space
/// \param first_bdim   Sets of box dimensions in the first collection of systems
/// \param second_umat  Box space transformation matrices for the first collection of systems
/// \param second_invu  Inverse transformation matrices for the second collection of systems
/// \param second_bdim  Sets of box dimensions in the second collection of systems
/// \param frm_first    Frame or system index to swap from the first object
/// \param frm_second   Frame or system index to swap from the second object
__device__ __forceinline__
void kSwapBoxes(double* first_umat, double* first_invu, double* first_bdim, double* second_umat,
                double* second_invu, double* second_bdim, const size_t frm_first,
                const size_t frm_second);
  
/// \brief Swap the contents of frames between two coordinate series.
///
/// Overloaded:
///   - Swap one frame between numbered indices of each series
///   - Swap multiple frames between each series
///
/// \param first       Abstract of the first coordinate series
/// \param second      Abstract of the second coordinate series
/// \param frm_first   Frame number to swap from the first object
/// \param frm_second  Frame number to swap from the second object
/// \param frm_pairs   List of pairs of frames to swap, indexing into the first object in each
///                    tuple's "x" member and indexing into the second object in each tuple's "y"
///                    member
/// \param pair_count  The trusted length of frm_pairs
/// \{
template <typename T> __global__ void __launch_bounds__(small_block_size, 4)
kSwapCoordinates(CoordinateSeriesWriter<T> first, const size_t frm_first,
                 CoordinateSeriesWriter<T> second, const size_t frm_second) {
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t advance = (blockDim.x * gridDim.x);
  const size_t padded_natom = devcRoundUp(first.natom, warp_size_int);
  const size_t first_offset = padded_natom * frm_first;
  const size_t second_offset = padded_natom * frm_second;
  pos = swapCoordinateSet(first.xcrd, first_offset, second.xcrd, second_offset, first.natom, pos,
                          0, padded_natom, advance);
  pos = swapCoordinateSet(first.ycrd, first_offset, second.ycrd, second_offset, first.natom, pos,
                          1, padded_natom, advance);
  pos = swapCoordinateSet(first.zcrd, first_offset, second.zcrd, second_offset, first.natom, pos,
                          2, padded_natom, advance);
  if (blockIdx.x == 0) {
    kSwapBoxes(first.umat, first.invu, first.boxdim, second.umat, second.invu, second.boxdim,
               frm_first, frm_second);
  }
}

template <typename T> __global__ void __launch_bounds__(small_block_size, 4)
kSwapCoordinates(CoordinateSeriesWriter<T> first, CoordinateSeriesWriter<T> second,
                 const int2* frm_pairs, int pair_count) {
  const int blocks_per_frm = (gridDim.x + pair_count - 1) / pair_count;
  const size_t padded_natom = devcRoundUp(first.natom, warp_size_int);
  int work_idx = blockIdx.x;
  while (work_idx < pair_count * blocks_per_frm) {
    const int pair_idx = work_idx / blocks_per_frm;
    const int2 swap_x_to_y = __ldca(&frm_pairs[pair_idx]);

    // Skip swap requests for invalid systems.  The requests will be checked on the CPU, where an
    // error will be raised.  The precision modes of each Condensate object will have been checked
    // prior to launching this kernel.
    if (swap_x_to_y.x < 0 || swap_x_to_y.x >= first.nframe ||
        swap_x_to_y.y < 0 || swap_x_to_y.y >= second.nframe) {
      work_idx += gridDim.x;
      continue;
    }

    // Compute the atom counts and offsets, then loop over X, Y, and Z coordinates
    const size_t first_offset = padded_natom * (size_t)(swap_x_to_y.x);
    const size_t second_offset = padded_natom * (size_t)(swap_x_to_y.y);
    const int wu_in_pair = (work_idx % blocks_per_frm);
    size_t pos = threadIdx.x + (wu_in_pair * blockDim.x);
    pos = swapCoordinateSet(first.xcrd, first_offset, second.xcrd, second_offset, first.natom, pos,
                            0, padded_natom, blockDim.x * blocks_per_frm);
    pos = swapCoordinateSet(first.ycrd, first_offset, second.ycrd, second_offset, first.natom, pos,
                            1, padded_natom, blockDim.x * blocks_per_frm);
    pos = swapCoordinateSet(first.zcrd, first_offset, second.zcrd, second_offset, first.natom, pos,
                            2, padded_natom, blockDim.x * blocks_per_frm);
    if (wu_in_pair == 0) {
      kSwapBoxes(first.umat, first.invu, first.boxdim, second.umat, second.invu, second.boxdim,
                 swap_x_to_y.x, swap_x_to_y.y);
    }

    // Advance to the next work unit in the swap
    work_idx += gridDim.x;
  }
}
/// \}

/// \brief Unroll the call to launch a coordinate swap within a coordinate series after detecting
///        and unrolling the template.  Overloading and descriptions of input parameters follow
///        from the corresponding kernel kSwapCoordinates(), above, with the addition of:
///
/// \param gpu         Details of the available GPU
/// \{
template <typename T>
void unrollSwapCoordLaunch(CoordinateSeriesWriter<void> *v_first, const size_t frm_first,
                           CoordinateSeriesWriter<void> *v_second, const size_t frm_second,
                           const GpuDetails &gpu) {
  CoordinateSeriesWriter<T> first  = restoreType<T>(v_first);
  CoordinateSeriesWriter<T> second = restoreType<T>(v_second);
  kSwapCoordinates<T><<<4 * gpu.getSMPCount(), small_block_size>>>(first, frm_first, second,
                                                                   frm_second);
}

template <typename T>
void unrollSwapCoordLaunch(CoordinateSeriesWriter<void> *v_first,
                           CoordinateSeriesWriter<void> *v_second, const int2* frm_pairs,
                           const int pair_count, const GpuDetails &gpu) {
  CoordinateSeriesWriter<T> first  = restoreType<T>(v_first);
  CoordinateSeriesWriter<T> second = restoreType<T>(v_second);
  kSwapCoordinates<T><<<4 * gpu.getSMPCount(), small_block_size>>>(first, second, frm_pairs,
                                                                   pair_count);
}
/// \}

} // namespace trajectory
} // namespace stormm

#endif
