// -*-c++-*-
#ifndef STORMM_HPC_PME_POTENTIAL_CUH
#define STORMM_HPC_PME_POTENTIAL_CUH

#include "copyright.h"

// Identify GPUs with large chip caches, hence more __shared__ memory available
#ifdef STORMM_USE_CUDA
#  if __CUDA_ARCH__ >= 750
#    define LARGE_CHIP_CACHE
#  endif
#endif

namespace stormm {
namespace energy {

#include "Numerics/accumulation.cui"

/// \brief An inline __device__ function to encapsulate force reductions for receiving atoms
///        across the warp, whether based on thread-specific accumulators in __shared__ memory or
///        in registers.  Atomic operations will store the accumulated values in the relevant
///        CellGrid neighbor list arrays.
///
/// \param batch_size        Number of atoms currently handled by the warp (as many as one atom
///                          per thread)
/// \param lane_idx          Index of the thread lane within the warp
/// \param base_plan_index   The plan index for sending atoms
/// \param lane_assignments  Array from which to draw lane shuffle reading assignments prior to
///                          reducing and storing the accumulated forces (if set to nullptr, no
///                          shuffle operations to position atoms will occur)
/// \param acc_fx            Accumulated force acting on the thread's atom in the X direction
/// \param acc_fy            Accumulated force acting on the thread's atom in the Y direction
/// \param acc_fz            Accumulated force acting on the thread's atom in the Z direction
/// \param frc_scale         Scaling factor to be applied to all forces
/// \param img_idx           Index of the atom in the relevant neighbor list
/// \param gbl_xfrc          Neighbor list accumulator array for Cartesian X forces
/// \param gbl_yfrc          Neighbor list accumulator array for Cartesian Y forces
/// \param gbl_zfrc          Neighbor list accumulator array for Cartesian Z forces
/// \param gbl_xfrc_ovrf     Overflow bits for forces acting along the Cartesian X direction
/// \param gbl_yfrc_ovrf     Overflow bits for forces acting along the Cartesian Y direction
/// \param gbl_zfrc_ovrf     Overflow bits for forces acting along the Cartesian Z direction
template <typename Tcoord, typename Tacc, typename  Tcalc, typename Tcoord4>
__device__ __forceinline__
void storeRecvForces(const int batch_size, const int lane_idx, const int base_plan_index,
                     const int* lane_assignments, const Tcalc acc_fx, const Tcalc acc_fy,
                     const Tcalc acc_fz, const uint img_idx,
                     CellGridWriter<Tcoord, Tacc, Tcalc, Tcoord4> cgw) {
  int reorg_idx, reduction_steps, lane_max;
  if (batch_size > half_warp_size_int) {
    lane_max = warp_size_int;
    reorg_idx = base_plan_index;
    reduction_steps = 0;
  }
  else if (batch_size > quarter_warp_size_int) {
    lane_max = half_warp_size_int;
    reorg_idx = base_plan_index + 1;
    reduction_steps = 1;
  }
  else {
    lane_max = quarter_warp_size_int;
    reorg_idx = base_plan_index + 2;
    reduction_steps = 2;
  }
  const int reduce_lane = lane_assignments[(reorg_idx * warp_size_int) + lane_idx];
  Tcalc my_acc_fx = SHFL(acc_fx, reduce_lane);
  Tcalc my_acc_fy = SHFL(acc_fy, reduce_lane);
  Tcalc my_acc_fz = SHFL(acc_fz, reduce_lane);
  if (reduction_steps == 1) {
    my_acc_fx += SHFL(my_acc_fx, lane_idx + half_warp_size_int);
    my_acc_fy += SHFL(my_acc_fy, lane_idx + half_warp_size_int);
    my_acc_fz += SHFL(my_acc_fz, lane_idx + half_warp_size_int);
  }
  else if (reduction_steps == 2) {
    my_acc_fx += SHFL(my_acc_fx, lane_idx + half_warp_size_int);
    my_acc_fy += SHFL(my_acc_fy, lane_idx + half_warp_size_int);
    my_acc_fz += SHFL(my_acc_fz, lane_idx + half_warp_size_int);
    my_acc_fx += SHFL(my_acc_fx, lane_idx + quarter_warp_size_int);
    my_acc_fy += SHFL(my_acc_fy, lane_idx + quarter_warp_size_int);
    my_acc_fz += SHFL(my_acc_fz, lane_idx + quarter_warp_size_int);
  }
  if (lane_idx < lane_max && img_idx != 0xffffffff) {
    atomicSplit(my_acc_fx * cgw.frc_scale, img_idx, cgw.xfrc, cgw.xfrc_ovrf);
    atomicSplit(my_acc_fy * cgw.frc_scale, img_idx, cgw.yfrc, cgw.yfrc_ovrf);
    atomicSplit(my_acc_fz * cgw.frc_scale, img_idx, cgw.zfrc, cgw.zfrc_ovrf);
  }
}

/// \brief Another inline function to store sending atom forces in either of the split
///        fixed-precision formats.  Descriptions of input parameters follow from
///        storeRecvForces(), above.
///
/// Overloaded:
///   - Accept int63_t (int, int) inputs
///   - Accept int95_t (long long int, int) inputs
/// \{
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
__device__ __forceinline__
void storeSendForces(const int lane_idx, const int base_plan_index, const llint acc_fx,
                     const llint acc_fy, const llint acc_fz, const int acc_fx_ovrf,
                     const int acc_fy_ovrf, const int acc_fz_ovrf, const uint img_idx,
                     CellGridWriter<Tcoord, Tacc, Tcalc, Tcoord4> cgw) {
  int95_t my_acc_fx = { acc_fx, acc_fx_ovrf };
  int95_t my_acc_fy = { acc_fy, acc_fy_ovrf };
  int95_t my_acc_fz = { acc_fz, acc_fz_ovrf };
  int lane_max = warp_size_int;
  if (base_plan_index >= 3) {
    const int shfl_lane = lane_idx + half_warp_size_int;
    my_acc_fx = splitFPSum(my_acc_fx, SHFL(my_acc_fx.x, shfl_lane), SHFL(my_acc_fx.y, shfl_lane));
    my_acc_fy = splitFPSum(my_acc_fy, SHFL(my_acc_fy.x, shfl_lane), SHFL(my_acc_fy.y, shfl_lane));
    my_acc_fz = splitFPSum(my_acc_fz, SHFL(my_acc_fz.x, shfl_lane), SHFL(my_acc_fz.y, shfl_lane));
    lane_max >>= 1;
  }
  if (base_plan_index == 6) {
    const int shfl_lane = lane_idx + quarter_warp_size_int;
    my_acc_fx = splitFPSum(my_acc_fx, SHFL(my_acc_fx.x, shfl_lane), SHFL(my_acc_fx.y, shfl_lane));
    my_acc_fy = splitFPSum(my_acc_fy, SHFL(my_acc_fy.x, shfl_lane), SHFL(my_acc_fy.y, shfl_lane));
    my_acc_fz = splitFPSum(my_acc_fz, SHFL(my_acc_fz.x, shfl_lane), SHFL(my_acc_fz.y, shfl_lane));
    lane_max >>= 1;
  }
  if (lane_idx < lane_max && img_idx != 0xffffffff) {
    atomicSplit(my_acc_fx, img_idx, cgw.xfrc, cgw.xfrc_ovrf);
    atomicSplit(my_acc_fy, img_idx, cgw.yfrc, cgw.yfrc_ovrf);
    atomicSplit(my_acc_fz, img_idx, cgw.zfrc, cgw.zfrc_ovrf);
  }
}

template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
__device__ __forceinline__
void storeSendForces(const int lane_idx, const int base_plan_index, const int acc_fx,
                     const int acc_fy, const int acc_fz, const int acc_fx_ovrf,
                     const int acc_fy_ovrf, const int acc_fz_ovrf, const uint img_idx,
                     CellGridWriter<Tcoord, Tacc, Tcalc, Tcoord4> cgw) {
  int2 my_acc_fx = { acc_fx, acc_fx_ovrf };
  int2 my_acc_fy = { acc_fy, acc_fy_ovrf };
  int2 my_acc_fz = { acc_fz, acc_fz_ovrf };
  int lane_max = warp_size_int;
  if (base_plan_index >= 3) {
    const int shfl_lane = lane_idx + half_warp_size_int;
    my_acc_fx = splitFPSum(my_acc_fx, SHFL(my_acc_fx.x, shfl_lane), SHFL(my_acc_fx.y, shfl_lane));
    my_acc_fy = splitFPSum(my_acc_fy, SHFL(my_acc_fy.x, shfl_lane), SHFL(my_acc_fy.y, shfl_lane));
    my_acc_fz = splitFPSum(my_acc_fz, SHFL(my_acc_fz.x, shfl_lane), SHFL(my_acc_fz.y, shfl_lane));
    lane_max >>= 1;
  }
  if (base_plan_index == 6) {
    const int shfl_lane = lane_idx + quarter_warp_size_int;
    my_acc_fx = splitFPSum(my_acc_fx, SHFL(my_acc_fx.x, shfl_lane), SHFL(my_acc_fx.y, shfl_lane));
    my_acc_fy = splitFPSum(my_acc_fy, SHFL(my_acc_fy.x, shfl_lane), SHFL(my_acc_fy.y, shfl_lane));
    my_acc_fz = splitFPSum(my_acc_fz, SHFL(my_acc_fz.x, shfl_lane), SHFL(my_acc_fz.y, shfl_lane));
    lane_max >>= 1;
  }
  if (lane_idx < lane_max && img_idx != 0xffffffff) {
    atomicSplit(my_acc_fx, img_idx, cgw.xfrc, cgw.xfrc_ovrf);
    atomicSplit(my_acc_fy, img_idx, cgw.yfrc, cgw.yfrc_ovrf);
    atomicSplit(my_acc_fz, img_idx, cgw.zfrc, cgw.zfrc_ovrf);
  }
}
/// \}

/// \brief Store the accumulated sending atom forces in temporary arrays.  This invokes a
///        conversion between floating point and split fixed-preicison integer representations of
///        the forces.  Overloading follows from storeSendForces(), above.  If there is a large
///        chip cache (detected by a CUDA architecture >= 750) then all storage is assumed to be in
///        __shared__ memory.  Otherwise, the overflow arrays are assumed to point to
///        block-exclusive regions of __global__ memory and a separate index is calculated.  This
///        function serves one dimension at a time to conserve register pressure.  The extra work
///        to calculate a new index is only undertaken if a card without sufficient L1 cache on
///        each streaming multiprocessor (SM) is in use.  Newer cards tend to have more L1 cache.
///
/// \param acc_frc       Accumulated, real-valued Cartesian X, Y, or Z direction force on each
///                      thread's particle
/// \param frc_scale     Scaling factor to apply to input real-valued forces
/// \param sh_frc        Array in __shared__ memory, holding primary accumulators of Cartesian X,
///                      Y, or Z direction forces
/// \param var_frc_ovrf  Variable nature (__shared__ or __global__) array for holding overflow
///                      bits of Cartesian X, Y, or Z force accumulations
/// \{
template <typename Tcalc> __device__ __forceinline__
void cacheSendForces(const Tcalc acc_frc, const Tcalc frc_scale, llint* sh_frc,
                     int* var_frc_ovrf) {
#ifdef LARGE_CHIP_CACHE
  const size_t var_idx = threadIdx.x;
#else
  const size_t var_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
#endif
  const size_t thr_idx_zu = threadIdx.x;
  const int95_t ifrc = int95Sum(sh_frc[thr_idx_zu], var_frc_ovrf[var_idx], acc_frc * frc_scale);
  sh_frc[thr_idx_zu] = ifrc.x;
  var_frc_ovrf[var_idx] = ifrc.y;
}

template <typename Tcalc> __device__ __forceinline__
void cacheSendForces(const Tcalc acc_frc, const Tcalc frc_scale, int* sh_frc, int* var_frc_ovrf) {
#ifdef LARGE_CHIP_CACHE
  const size_t var_idx = threadIdx.x;
#else
  const size_t var_idx = (blockIdx.x * blockDim.x) + threadIdx.x;
#endif
  const size_t thr_idx_zu = threadIdx.x;
  const int2 ifrc = int63Sum(sh_frc[thr_idx_zu], var_frc_ovrf[var_idx], acc_frc * frc_scale);
  sh_frc[thr_idx_zu] = ifrc.x;
  var_frc_ovrf[var_idx] = ifrc.y;
}
/// \}

} // namespace energy
} // namespace stormm

// Clear hardware-dependent definitions
#ifdef STORMM_USE_CUDA
#  if __CUDA_ARCH__ >= 750
#    undef LARGE_CHIP_CACHE
#  endif
#endif

#endif
