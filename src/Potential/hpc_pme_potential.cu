// -*-c++-*-
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "MolecularMechanics/mm_controls.h"
#include "Synthesis/synthesis_abstracts.h"
#include "pmigrid.h"
#include "hpc_pme_potential.h"
#include "hpc_pme_potential.cuh"

namespace stormm {
namespace energy {

using mm::MMControlKit;
using mm::MolecularMechanicsControls;
using synthesis::SyNonbondedKit;
  
#include "Accelerator/syncwarp.cui"
#include "Math/rounding.cui"
#include "Structure/local_arrangement.cui"

/// \brief A device function for evaluating a local exclusion mask based on a topology index delta.
///        This is equivalent to evaluateLocalMask() in the LocalExclusionMask library (see
///        local_exclusionmask.h).
///
__device__ __forceinline__ bool devcEvaluateLocalMask(int atom_i, int atom_j, ullint prof,
                                                      const uint2* secondary_ptr) {
  const int del_ij = atom_j - atom_i;
  switch (prof & lmask_mode_bitmask) {
  case lmask_mode_a:
    {
      return (abs(del_ij) <= lmask_long_local_span &&
              ((prof >> (lmask_long_local_span + del_ij)) & 0x1));
    }
    break;
  case lmask_mode_b:
    {
      const int abs_dij = abs(del_ij);
      if (abs_dij > lmask_b_max_reach) {
        return false;
      }
      else if (abs_dij <= lmask_short_local_span) {
        return ((prof >> (lmask_short_local_span + del_ij)) & 0x1);
      }
      else if (del_ij > 0) {
        const int upper_shft = ((prof & lmask_b_upper_shft) >> lmask_b_upper_shft_pos);
        const int upper_mask_start = lmask_short_local_span + upper_shft;
        const int rel_ij = del_ij - upper_mask_start - 1;
        return (rel_ij >= 0 && rel_ij < lmask_short_extra_span &&
                ((prof >> (rel_ij + lmask_b_upper_mask_pos)) & 0x1));
      }
      else {

        // The only remaining case is that del_ij < 0
        const int lower_shft = ((prof & lmask_b_lower_shft) >> lmask_b_lower_shft_pos);
        const int lower_mask_start = -lmask_short_local_span - lower_shft - lmask_short_extra_span;
        const int rel_ij = del_ij - lower_mask_start;
        return (rel_ij >= 0 && rel_ij < lmask_short_extra_span &&
                ((prof >> (rel_ij + lmask_b_lower_mask_pos)) & 0x1));
      }
    }  
    break;
  case lmask_mode_c:
    {
      if (abs(del_ij) <= lmask_short_local_span) {
        return ((prof >> (lmask_short_local_span + del_ij)) & 0x1);
      }

      // Forming the unsigned long long int on the r.h.s. and then converting it to a signed
      // short int will translate the bit string appropriately.
      const int alt_mask_shft = static_cast<short int>((prof & lmask_c_shft) >> lmask_c_shft_pos);

      // Run the shift in terms of the index atom
      const int rel_ij = del_ij - alt_mask_shft;
      return (rel_ij >= 0 && rel_ij < lmask_long_extra_span &&
              ((prof >> (rel_ij + lmask_c_alt_mask_pos)) & 0x1));
    }
    break;
  case lmask_mode_d:
    {
      if (abs(del_ij) <= lmask_short_local_span) {
        return ((prof >> (lmask_short_local_span + del_ij)) & 0x1);
      }

      // This is the best possible path.  Obtain the number of masks and loop over all of them.
      const size_t nmasks = ((prof & lmask_d_array_cnt) >> lmask_d_array_cnt_pos);
      const size_t start_idx = ((prof & lmask_d_array_idx) >> lmask_d_array_idx_pos);
      for (size_t i = 0; i < nmasks; i++) {
        const uint2 tmask = secondary_ptr[start_idx + i];
        const int tmask_x = tmask.x;
        if (del_ij >= tmask_x && del_ij < tmask_x + 32 &&
            ((tmask.y >> (del_ij - tmask_x)) & 0x1)) {
          return true;
        }
      }  
      return false;
    }
    break;
  case lmask_mode_e:
    {
      // This is the best possible path and there is no local exclusion arrangement to test.  Loop
      // over all the masks.
      const size_t nmasks = ((prof & lmask_e_array_cnt) >> lmask_e_array_cnt_pos);
      const size_t start_idx = (prof & lmask_e_array_idx);
      for (size_t i = 0; i < nmasks; i++) {
        const uint2 tmask = secondary_ptr[start_idx + i];
        const int tmask_x = tmask.x;
        if (del_ij >= tmask_x && del_ij < tmask_x + 32 &&
            ((tmask.y >> (del_ij - tmask_x)) & 0x1)) {
          return true;
        }
      }
      return false;
    }
    break;
  case lmask_mode_f:
    break;
  default:
    break;
  }
  __builtin_unreachable();
}

// Detect the chip cache size.  Turing chips count as having a "large chip cache" not because
// they can allocate up to 100 kB of __shared__ memory but because they can only have 1024 threads
// activ eon a single block, which would allocate somewhat less than Maxwell, Pascal, Volta,
// Ampere, or Lovelace / Hopper cards would require.
#ifdef STORMM_USE_CUDA
#  if __CUDA_ARCH__ >= 750
#    define LARGE_CHIP_CACHE
#  endif
#endif

#define PMENB_BLOCK_MULTIPLICITY 2

// Single-precision tile evaluation
#define TCALC float
#  define TCALC2 float2
#  define TCALC3 float3
#  define TCALC4 float4
#  define TCALC_IS_SINGLE
#  define TCOORD_IS_REAL

// Begin with a float32_t coordinate representation, natural for the float32_t arithmetic mode.
#  define TCOORD  float
#  define TCOORD4 float4
#  define TACC    int

// Other definitions associated with 32-bit floating-point arithmetic
#  define SQRT_FUNC sqrtf
#  define LLCONV_FUNC __float2ll_rn
  
// Compile the kernels with or without energy and force computations, dual neighbor lists,
// provisions for small box sizes, and clash forgiveness.
#  define COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define DUAL_GRIDS
#        define TINY_BOX
#          define CLASH_FORGIVENESS
#            define PMENB_WARPS_PER_BLOCK 16
#            define KERNEL_NAME kffTowerPlateFEDualTinyNonClash
#              include "tower_plate_pairs.cui"
#            undef KERNEL_NAME
#            undef PMENB_WARPS_PER_BLOCK
#          undef CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 16
#          define KERNEL_NAME kffTowerPlateFEDualTiny
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 16
#          define KERNEL_NAME kffTowerPlateFEDualNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kffTowerPlateFEDual
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 16
#          define KERNEL_NAME kffTowerPlateFETinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kffTowerPlateFETiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kffTowerPlateFENonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 16
#      define KERNEL_NAME kffTowerPlateFE
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 16
#          define KERNEL_NAME kffTowerPlateFXDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kffTowerPlateFXDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kffTowerPlateFXDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 18
#      define KERNEL_NAME kffTowerPlateFXDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kffTowerPlateFXTinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 16
#      define KERNEL_NAME kffTowerPlateFXTiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef TINY_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 16
#      define KERNEL_NAME kffTowerPlateFXNonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 18
#    define KERNEL_NAME kffTowerPlateFX
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 16
#          define KERNEL_NAME kffTowerPlateXEDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kffTowerPlateXEDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kffTowerPlateXEDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 18
#      define KERNEL_NAME kffTowerPlateXEDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kffTowerPlateXETinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 16
#      define KERNEL_NAME kffTowerPlateXETiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef TINY_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 16
#      define KERNEL_NAME kffTowerPlateXENonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 16
#    define KERNEL_NAME kffTowerPlateXE
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#  undef COMPUTE_ENERGY

#  undef TCOORD
#  undef TCOORD4
#  undef TACC

// Compile additional kernels for the float64_t coordinate representation.
#  define TCOORD double
#  define TCOORD4 double4
#  define TACC   llint
#  define TCOORD_IS_LONG

// Compile the kernels with or without energy and force computations, dual neighbor lists,
// provisions for small box sizes, and clash forgiveness.
#  define COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define DUAL_GRIDS
#        define TINY_BOX
#          define CLASH_FORGIVENESS
#            define PMENB_WARPS_PER_BLOCK 16
#            define KERNEL_NAME kdfTowerPlateFEDualTinyNonClash
#              include "tower_plate_pairs.cui"
#            undef KERNEL_NAME
#            undef PMENB_WARPS_PER_BLOCK
#          undef CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 16
#          define KERNEL_NAME kdfTowerPlateFEDualTiny
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 16
#          define KERNEL_NAME kdfTowerPlateFEDualNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kdfTowerPlateFEDual
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 16
#          define KERNEL_NAME kdfTowerPlateFETinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kdfTowerPlateFETiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kdfTowerPlateFENonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 16
#      define KERNEL_NAME kdfTowerPlateFE
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 16
#          define KERNEL_NAME kdfTowerPlateFXDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kdfTowerPlateFXDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kdfTowerPlateFXDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 16
#      define KERNEL_NAME kdfTowerPlateFXDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kdfTowerPlateFXTinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 16
#      define KERNEL_NAME kdfTowerPlateFXTiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef TINY_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 16
#      define KERNEL_NAME kdfTowerPlateFXNonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 16
#    define KERNEL_NAME kdfTowerPlateFX
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 16
#          define KERNEL_NAME kdfTowerPlateXEDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kdfTowerPlateXEDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kdfTowerPlateXEDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 16
#      define KERNEL_NAME kdfTowerPlateXEDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 16
#        define KERNEL_NAME kdfTowerPlateXETinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 16
#      define KERNEL_NAME kdfTowerPlateXETiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef TINY_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 16
#      define KERNEL_NAME kdfTowerPlateXENonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 16
#    define KERNEL_NAME kdfTowerPlateXE
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#  undef COMPUTE_ENERGY

#  undef LLCONV_FUNC
#  undef SQRT_FUNC
  
#  undef TCOORD
#  undef TCOORD4
#  undef TACC
#  undef TCOORD_IS_LONG
  
#  undef TCOORD_IS_REAL
#  undef TCALC_IS_SINGLE
#  undef TCALC2
#  undef TCALC3
#  undef TCALC4
#undef TCALC
  
// Double-precision tile evaluation
#define TCALC double
#  define TCALC2 double2
#  define TCALC3 double3
#  define TCALC4 double4
#  define TCOORD_IS_REAL

// Begin with a float64_t coordinate representation, appropriate for the float64_t arithmetic mode.
#  define TCOORD double
#  define TCOORD4 double4
#  define TACC   llint
#  define TCOORD_IS_LONG

// Other definitions associated with 64-bit floating-point arithmetic
#  define SQRT_FUNC sqrt
#  define LLCONV_FUNC __double2ll_rn

// Compile the kernels with or without energy and force computations, dual neighbor lists,
// provisions for small box sizes, and clash forgiveness.
#  define COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define DUAL_GRIDS
#        define TINY_BOX
#          define CLASH_FORGIVENESS
#            define PMENB_WARPS_PER_BLOCK 12
#            define KERNEL_NAME kddTowerPlateFEDualTinyNonClash
#              include "tower_plate_pairs.cui"
#            undef KERNEL_NAME
#            undef PMENB_WARPS_PER_BLOCK
#          undef CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 12
#          define KERNEL_NAME kddTowerPlateFEDualTiny
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 12
#          define KERNEL_NAME kddTowerPlateFEDualNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kddTowerPlateFEDual
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 12
#          define KERNEL_NAME kddTowerPlateFETinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kddTowerPlateFETiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kddTowerPlateFENonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 12
#      define KERNEL_NAME kddTowerPlateFE
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 12
#          define KERNEL_NAME kddTowerPlateFXDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kddTowerPlateFXDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kddTowerPlateFXDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 12
#      define KERNEL_NAME kddTowerPlateFXDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kddTowerPlateFXTinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 12
#      define KERNEL_NAME kddTowerPlateFXTiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef TINY_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 12
#      define KERNEL_NAME kddTowerPlateFXNonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 12
#    define KERNEL_NAME kddTowerPlateFX
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 12
#          define KERNEL_NAME kddTowerPlateXEDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kddTowerPlateXEDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kddTowerPlateXEDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 12
#      define KERNEL_NAME kddTowerPlateXEDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kddTowerPlateXETinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 12
#      define KERNEL_NAME kddTowerPlateXETiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef TINY_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 12
#      define KERNEL_NAME kddTowerPlateXENonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 12
#    define KERNEL_NAME kddTowerPlateXE
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#  undef COMPUTE_ENERGY

#  undef TCOORD
#  undef TCOORD4
#  undef TACC
#  undef TCOORD_IS_LONG

// Continue with float32_t representations of the coordinates, to allow for high-precision
// calculations to be applied to otherwise low-precision coordinate representations.
#  define TCOORD float
#  define TCOORD4 float4
#  define TACC   int

// Compile the kernels with or without energy and force computations, dual neighbor lists,
// provisions for small box sizes, and clash forgiveness.
#  define COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define DUAL_GRIDS
#        define TINY_BOX
#          define CLASH_FORGIVENESS
#            define PMENB_WARPS_PER_BLOCK 12
#            define KERNEL_NAME kfdTowerPlateFEDualTinyNonClash
#              include "tower_plate_pairs.cui"
#            undef KERNEL_NAME
#            undef PMENB_WARPS_PER_BLOCK
#          undef CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 12
#          define KERNEL_NAME kfdTowerPlateFEDualTiny
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 12
#          define KERNEL_NAME kfdTowerPlateFEDualNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kfdTowerPlateFEDual
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 12
#          define KERNEL_NAME kfdTowerPlateFETinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kfdTowerPlateFETiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kfdTowerPlateFENonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 12
#      define KERNEL_NAME kfdTowerPlateFE
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 12
#          define KERNEL_NAME kfdTowerPlateFXDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kfdTowerPlateFXDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kfdTowerPlateFXDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 12
#      define KERNEL_NAME kfdTowerPlateFXDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kfdTowerPlateFXTinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 12
#      define KERNEL_NAME kfdTowerPlateFXTiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef TINY_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 12
#      define KERNEL_NAME kfdTowerPlateFXNonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 12
#    define KERNEL_NAME kfdTowerPlateFX
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 12
#          define KERNEL_NAME kfdTowerPlateXEDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kfdTowerPlateXEDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kfdTowerPlateXEDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 12
#      define KERNEL_NAME kfdTowerPlateXEDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 12
#        define KERNEL_NAME kfdTowerPlateXETinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 12
#      define KERNEL_NAME kfdTowerPlateXETiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef TINY_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 12
#      define KERNEL_NAME kfdTowerPlateXENonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 12
#    define KERNEL_NAME kfdTowerPlateXE
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#  undef COMPUTE_ENERGY

#  undef LLCONV_FUNC
#  undef SQRT_FUNC

#  undef TCOORD
#  undef TCOORD4
#  undef TACC
  
#  undef TCOORD_IS_REAL
#  undef TCALC2
#  undef TCALC3
#  undef TCALC4
#undef TCALC

#undef PMENB_BLOCK_MULTIPLICITY

// Clear hardware-dependent definitions
#ifdef STORMM_USE_CUDA
#  if __CUDA_ARCH__ >= 750
#    undef LARGE_CHIP_CACHE
#  endif
#endif

#ifdef STORMM_USE_CUDA
//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes queryPMEPairsKernelRequirements(const PrecisionModel coord_prec,
                                                          const PrecisionModel calc_prec,
                                                          const NeighborListKind neighbor_list,
                                                          const EvaluateForce eval_frc,
                                                          const EvaluateEnergy eval_nrg,
                                                          const TinyBoxPresence has_tiny_box,
                                                          const ClashResponse clash_handling) {

  // As with other kernel querying functions, the kernel manager calling this function will have
  // specifications of the GPU in use.  It is the overall thread occupancy and multiplicity of each
  // kernel that this function must return.
  cudaFuncAttributes result;
  cudaError_t cfa;
  switch (clash_handling) {
  case ClashResponse::NONE:
    switch (coord_prec) {
    case PrecisionModel::DOUBLE:
      switch (calc_prec) {
      case PrecisionModel::DOUBLE:
        switch (eval_frc) {
        case EvaluateForce::YES:
          switch (eval_nrg) {
          case EvaluateEnergy::YES:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFEDualTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFEDual);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFETiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFE);
                break;
              }
              break;
            }
            break;
          case EvaluateEnergy::NO:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFXDualTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFXDual);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFXTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFX);
                break;
              }
              break;
            }
            break;
          }
          break;
        case EvaluateForce::NO:

          // If the force is not being evaluated, the energy must be required
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kddTowerPlateXEDualTiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kddTowerPlateXEDual);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kddTowerPlateXETiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kddTowerPlateXE);
              break;
            }
            break;
          }
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (eval_frc) {
        case EvaluateForce::YES:
          switch (eval_nrg) {
          case EvaluateEnergy::YES:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFEDualTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFEDual);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFETiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFE);
                break;
              }
              break;
            }
            break;
          case EvaluateEnergy::NO:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFXDualTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFXDual);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFXTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFX);
                break;
              }
              break;
            }
            break;
          }
          break;
        case EvaluateForce::NO:

          // If the force is not being evaluated, the energy must be required
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kdfTowerPlateXEDualTiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kdfTowerPlateXEDual);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kdfTowerPlateXETiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kdfTowerPlateXE);
              break;
            }
            break;
          }
          break;
        }
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      switch (calc_prec) {
      case PrecisionModel::DOUBLE:
        switch (eval_frc) {
        case EvaluateForce::YES:
          switch (eval_nrg) {
          case EvaluateEnergy::YES:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFEDualTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFEDual);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFETiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFE);
                break;
              }
              break;
            }
            break;
          case EvaluateEnergy::NO:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFXDualTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFXDual);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFXTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFX);
                break;
              }
              break;
            }
            break;
          }
          break;
        case EvaluateForce::NO:

          // If the force is not being evaluated, the energy must be required
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kfdTowerPlateXEDualTiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kfdTowerPlateXEDual);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kfdTowerPlateXETiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kfdTowerPlateXE);
              break;
            }
            break;
          }
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (eval_frc) {
        case EvaluateForce::YES:
          switch (eval_nrg) {
          case EvaluateEnergy::YES:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFEDualTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFEDual);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFETiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFE);
                break;
              }
              break;
            }
            break;
          case EvaluateEnergy::NO:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFXDualTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFXDual);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFXTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFX);
                break;
              }
              break;
            }
            break;
          }
          break;
        case EvaluateForce::NO:

          // If the force is not being evaluated, the energy must be required
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kffTowerPlateXEDualTiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kffTowerPlateXEDual);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kffTowerPlateXETiny);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kffTowerPlateXE);
              break;
            }
            break;
          }
          break;
        }
        break;
      }
      break;
    }
    break;
  case ClashResponse::FORGIVE:
    switch (coord_prec) {
    case PrecisionModel::DOUBLE:
      switch (calc_prec) {
      case PrecisionModel::DOUBLE:
        switch (eval_frc) {
        case EvaluateForce::YES:
          switch (eval_nrg) {
          case EvaluateEnergy::YES:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFEDualTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFEDualNonClash);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFETinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFENonClash);
                break;
              }
              break;
            }
            break;
          case EvaluateEnergy::NO:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFXDualTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFXDualNonClash);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFXTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kddTowerPlateFXNonClash);
                break;
              }
              break;
            }
            break;
          }
          break;
        case EvaluateForce::NO:

          // If the force is not being evaluated, the energy must be required
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kddTowerPlateXEDualTinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kddTowerPlateXEDualNonClash);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kddTowerPlateXETinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kddTowerPlateXENonClash);
              break;
            }
            break;
          }
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (eval_frc) {
        case EvaluateForce::YES:
          switch (eval_nrg) {
          case EvaluateEnergy::YES:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFEDualTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFEDualNonClash);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFETinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFENonClash);
                break;
              }
              break;
            }
            break;
          case EvaluateEnergy::NO:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFXDualTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFXDualNonClash);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFXTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kdfTowerPlateFXNonClash);
                break;
              }
              break;
            }
            break;
          }
          break;
        case EvaluateForce::NO:

          // If the force is not being evaluated, the energy must be required
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kdfTowerPlateXEDualTinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kdfTowerPlateXEDualNonClash);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kdfTowerPlateXETinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kdfTowerPlateXENonClash);
              break;
            }
            break;
          }
          break;
        }
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      switch (calc_prec) {
      case PrecisionModel::DOUBLE:
        switch (eval_frc) {
        case EvaluateForce::YES:
          switch (eval_nrg) {
          case EvaluateEnergy::YES:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFEDualTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFEDualNonClash);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFETinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFENonClash);
                break;
              }
              break;
            }
            break;
          case EvaluateEnergy::NO:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFXDualTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFXDualNonClash);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFXTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kfdTowerPlateFXNonClash);
                break;
              }
              break;
            }
            break;
          }
          break;
        case EvaluateForce::NO:

          // If the force is not being evaluated, the energy must be required
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kfdTowerPlateXEDualTinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kfdTowerPlateXEDualNonClash);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kfdTowerPlateXETinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kfdTowerPlateXENonClash);
              break;
            }
            break;
          }
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (eval_frc) {
        case EvaluateForce::YES:
          switch (eval_nrg) {
          case EvaluateEnergy::YES:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFEDualTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFEDualNonClash);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFETinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFENonClash);
                break;
              }
              break;
            }
            break;
          case EvaluateEnergy::NO:
            switch (neighbor_list) {
            case NeighborListKind::DUAL:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFXDualTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFXDualNonClash);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFXTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kffTowerPlateFXNonClash);
                break;
              }
              break;
            }
            break;
          }
          break;
        case EvaluateForce::NO:

          // If the force is not being evaluated, the energy must be required
          switch (neighbor_list) {
          case NeighborListKind::DUAL:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kffTowerPlateXEDualTinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kffTowerPlateXEDualNonClash);
              break;
            }
            break;
          case NeighborListKind::MONO:
            switch (has_tiny_box) {
            case TinyBoxPresence::YES:
              cfa = cudaFuncGetAttributes(&result, kffTowerPlateXETinyNonClash);
              break;
            case TinyBoxPresence::NO:
              cfa = cudaFuncGetAttributes(&result, kffTowerPlateXENonClash);
              break;
            }
            break;
          }
          break;
        }
        break;
      }
      break;
    }
    break;
  }

  // Check for errors
  if (cfa != cudaSuccess) {

    // Construct the appropriate error message
    std::string error_message("Error obtaining attributes from kernel k");
    switch (coord_prec) {
    case PrecisionModel::DOUBLE:
      error_message += "d";
      break;
    case PrecisionModel::SINGLE:
      error_message += "s";
      break;
    }
    switch (calc_prec) {
    case PrecisionModel::DOUBLE:
      error_message += "d";
      break;
    case PrecisionModel::SINGLE:
      error_message += "s";
      break;
    }
    error_message += "TowerPlate";
    switch (eval_frc) {
    case EvaluateForce::YES:
      error_message += "Force";
      break;
    case EvaluateForce::NO:
      break;
    }
    switch (eval_nrg) {
    case EvaluateEnergy::YES:
      error_message += "Energy";
      break;
    case EvaluateEnergy::NO:
      break;
    }
    switch (neighbor_list) {
    case NeighborListKind::DUAL:
      error_message += "Dual";
      break;
    case NeighborListKind::MONO:
      break;
    }
    switch (clash_handling) {
    case ClashResponse::FORGIVE:
      error_message += "NonClash";
      break;
    case ClashResponse::NONE:
      break;
    }

    // Report the error
    rtErr(error_message, "queryPMEPairsKernelRequirements");
  }
  return result;
}
#endif

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<double, double4> &nrg_tab,
                           CellGridWriter<double, llint, double, double4> *cgw, TilePlan *tlpn,
                           ScoreCardWriter *scw, MMControlKit<double> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // DOUBLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kddTowerPlateFENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *scw, *cgw,
                                                      *ctrl);
        break;
      case EvaluateEnergy::NO:
        kddTowerPlateFXNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *cgw, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kddTowerPlateXENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                    clash_distance, clash_ratio, *scw, *cgw,
                                                    *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kddTowerPlateFE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kddTowerPlateFX<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kddTowerPlateXE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<double, double4> &nrg_tab,
                           const PsSynthesisBorders &sysbrd,
                           CellGridWriter<double, llint, double, double4> *cgw, TilePlan *tlpn,
                           ScoreCardWriter *scw, MMControlKit<double> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // DOUBLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kddTowerPlateFETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kddTowerPlateFXTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                          clash_distance, clash_ratio, *cgw,
                                                          *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kddTowerPlateXETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                        clash_distance, clash_ratio, *scw, *cgw,
                                                        *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kddTowerPlateFETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw,
                                                  *cgw, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kddTowerPlateFXTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *cgw,
                                                  *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kddTowerPlateXETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw, *cgw,
                                                *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<double, double4> &nrg_tab,
                           CellGridWriter<float, int, float, float4> *cgw, TilePlan *tlpn,
                           ScoreCardWriter *scw, MMControlKit<double> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // DOUBLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kfdTowerPlateFENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *scw, *cgw,
                                                      *ctrl);
        break;
      case EvaluateEnergy::NO:
        kfdTowerPlateFXNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *cgw, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kfdTowerPlateXENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                    clash_distance, clash_ratio, *scw, *cgw,
                                                    *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kfdTowerPlateFE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kfdTowerPlateFX<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kfdTowerPlateXE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<double, double4> &nrg_tab,
                           const PsSynthesisBorders &sysbrd,
                           CellGridWriter<float, int, float, float4> *cgw, TilePlan *tlpn,
                           ScoreCardWriter *scw, MMControlKit<double> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // DOUBLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kfdTowerPlateFETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kfdTowerPlateFXTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                          clash_distance, clash_ratio, *cgw,
                                                          *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kfdTowerPlateXETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                        clash_distance, clash_ratio, *scw, *cgw,
                                                        *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kfdTowerPlateFETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw,
                                                  *cgw, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kfdTowerPlateFXTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *cgw,
                                                  *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kfdTowerPlateXETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw, *cgw,
                                                *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<double, double4> &nrg_tab,
                           CellGridWriter<double, llint, double, double4> *cgw_qq,
                           CellGridWriter<double, llint, double, double4> *cgw_lj,
                           TilePlan *tlpn, ScoreCardWriter *scw, MMControlKit<double> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // DOUBLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kddTowerPlateFEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw_qq, *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kddTowerPlateFXDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *cgw_qq,
                                                          *cgw_lj, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kddTowerPlateXEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *scw,
                                                        *cgw_qq, *cgw_lj, *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kddTowerPlateFEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                  *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kddTowerPlateFXDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq, *cgw_lj,
                                                  *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kddTowerPlateXEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                *cgw_lj, *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<double, double4> &nrg_tab,
                           const PsSynthesisBorders &sysbrd,
                           CellGridWriter<double, llint, double, double4> *cgw_qq,
                           CellGridWriter<double, llint, double, double4> *cgw_lj,
                           TilePlan *tlpn, ScoreCardWriter *scw, MMControlKit<double> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // DOUBLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kddTowerPlateFEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              sysbrd, clash_distance, clash_ratio,
                                                              *scw, *cgw_qq, *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kddTowerPlateFXDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              sysbrd, clash_distance, clash_ratio,
                                                              *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kddTowerPlateXEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            sysbrd, clash_distance, clash_ratio,
                                                            *scw, *cgw_qq, *cgw_lj, *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kddTowerPlateFEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kddTowerPlateFXDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kddTowerPlateXEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw,
                                                    *cgw_qq, *cgw_lj, *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<double, double4> &nrg_tab,
                           CellGridWriter<float, int, float, float4> *cgw_qq,
                           CellGridWriter<float, int, float, float4> *cgw_lj, TilePlan *tlpn,
                           ScoreCardWriter *scw, MMControlKit<double> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // DOUBLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kfdTowerPlateFEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw_qq, *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kfdTowerPlateFXDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *cgw_qq,
                                                          *cgw_lj, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kfdTowerPlateXEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *scw,
                                                        *cgw_qq, *cgw_lj, *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kfdTowerPlateFEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                  *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kfdTowerPlateFXDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                  *cgw_lj, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kfdTowerPlateXEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                *cgw_lj, *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<double, double4> &nrg_tab,
                           const PsSynthesisBorders &sysbrd,
                           CellGridWriter<float, int, float, float4> *cgw_qq,
                           CellGridWriter<float, int, float, float4> *cgw_lj, TilePlan *tlpn,
                           ScoreCardWriter *scw, MMControlKit<double> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // DOUBLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kfdTowerPlateFEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              sysbrd, clash_distance, clash_ratio,
                                                              *scw, *cgw_qq, *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kfdTowerPlateFXDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              sysbrd, clash_distance, clash_ratio,
                                                              *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kfdTowerPlateXEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw_qq, *cgw_lj, *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kfdTowerPlateFEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kfdTowerPlateFXDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kfdTowerPlateXEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw,
                                                    *cgw_qq, *cgw_lj, *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<float, float4> &nrg_tab,
                           CellGridWriter<float, int, float, float4> *cgw, TilePlan *tlpn,
                           ScoreCardWriter *scw, MMControlKit<float> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // SINGLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kffTowerPlateFENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *scw, *cgw,
                                                      *ctrl);
        break;
      case EvaluateEnergy::NO:
        kffTowerPlateFXNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *cgw, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kffTowerPlateXENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                    clash_distance, clash_ratio, *scw, *cgw,
                                                    *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kffTowerPlateFE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kffTowerPlateFX<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kffTowerPlateXE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<float, float4> &nrg_tab, const PsSynthesisBorders &sysbrd,
                           CellGridWriter<float, int, float, float4> *cgw, TilePlan *tlpn,
                           ScoreCardWriter *scw, MMControlKit<float> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // SINGLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kffTowerPlateFETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kffTowerPlateFXTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                          clash_distance, clash_ratio, *cgw,
                                                          *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kffTowerPlateXETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                        clash_distance, clash_ratio, *scw, *cgw,
                                                        *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kffTowerPlateFETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw,
                                                  *cgw, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kffTowerPlateFXTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *cgw,
                                                  *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kffTowerPlateXETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw, *cgw,
                                                *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<float, float4> &nrg_tab,
                           CellGridWriter<double, llint, double, double4> *cgw, TilePlan *tlpn,
                           ScoreCardWriter *scw, MMControlKit<float> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // SINGLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kdfTowerPlateFENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *scw, *cgw,
                                                      *ctrl);
        break;
      case EvaluateEnergy::NO:
        kdfTowerPlateFXNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *cgw, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kdfTowerPlateXENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                    clash_distance, clash_ratio, *scw, *cgw,
                                                    *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kdfTowerPlateFE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kdfTowerPlateFX<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kdfTowerPlateXE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<float, float4> &nrg_tab, const PsSynthesisBorders &sysbrd,
                           CellGridWriter<double, llint, double, double4> *cgw, TilePlan *tlpn,
                           ScoreCardWriter *scw, MMControlKit<float> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // SINGLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kdfTowerPlateFETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kdfTowerPlateFXTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                          clash_distance, clash_ratio, *cgw,
                                                          *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kdfTowerPlateXETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                        clash_distance, clash_ratio, *scw, *cgw,
                                                        *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kdfTowerPlateFETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw,
                                                  *cgw, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kdfTowerPlateFXTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *cgw,
                                                  *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kdfTowerPlateXETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw, *cgw,
                                                *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<float, float4> &nrg_tab,
                           CellGridWriter<float, int, float, float4> *cgw_qq,
                           CellGridWriter<float, int, float, float4> *cgw_lj,
                           TilePlan *tlpn, ScoreCardWriter *scw, MMControlKit<float> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // SINGLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kffTowerPlateFEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw_qq, *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kffTowerPlateFXDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *cgw_qq,
                                                          *cgw_lj, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kffTowerPlateXEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *scw, *cgw_qq,
                                                        *cgw_lj, *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kffTowerPlateFEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                  *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kffTowerPlateFXDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                  *cgw_lj, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kffTowerPlateXEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                *cgw_lj, *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<float, float4> &nrg_tab, const PsSynthesisBorders &sysbrd,
                           CellGridWriter<float, int, float, float4> *cgw_qq,
                           CellGridWriter<float, int, float, float4> *cgw_lj,
                           TilePlan *tlpn, ScoreCardWriter *scw, MMControlKit<float> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // SINGLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kffTowerPlateFEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              sysbrd, clash_distance, clash_ratio,
                                                              *scw, *cgw_qq, *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kffTowerPlateFXDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              sysbrd, clash_distance, clash_ratio,
                                                              *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kffTowerPlateXEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            sysbrd, clash_distance, clash_ratio,
                                                            *scw, *cgw_qq, *cgw_lj, *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kffTowerPlateFEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kffTowerPlateFXDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kffTowerPlateXEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw,
                                                    *cgw_qq, *cgw_lj, *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<float, float4> &nrg_tab,
                           CellGridWriter<double, llint, double, double4> *cgw_qq,
                           CellGridWriter<double, llint, double, double4> *cgw_lj,
                           TilePlan *tlpn, ScoreCardWriter *scw, MMControlKit<float> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // SINGLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kdfTowerPlateFEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw_qq, *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kdfTowerPlateFXDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *cgw_qq,
                                                          *cgw_lj, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kdfTowerPlateXEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *scw,
                                                        *cgw_qq, *cgw_lj, *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kdfTowerPlateFEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                  *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kdfTowerPlateFXDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                  *cgw_lj, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kdfTowerPlateXEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                *cgw_lj, *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<float, float4> &nrg_tab, const PsSynthesisBorders &sysbrd,
                           CellGridWriter<double, llint, double, double4> *cgw_qq,
                           CellGridWriter<double, llint, double, double4> *cgw_lj,
                           TilePlan *tlpn, ScoreCardWriter *scw, MMControlKit<float> *ctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const int2 bt_tp, const double clash_distance,
                           const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // SINGLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kdfTowerPlateFEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              sysbrd, clash_distance, clash_ratio,
                                                              *scw, *cgw_qq, *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kdfTowerPlateFXDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              sysbrd, clash_distance, clash_ratio,
                                                              *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kdfTowerPlateXEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            sysbrd, clash_distance, clash_ratio,
                                                            *scw, *cgw_qq, *cgw_lj, *ctrl);
      break;
    }
  }
  else {
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        kdfTowerPlateFEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        break;
      case EvaluateEnergy::NO:
        kdfTowerPlateFXDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case EvaluateForce::NO:
      kdfTowerPlateXEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, sysbrd, *scw,
                                                    *cgw_qq, *cgw_lj, *ctrl);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const PrecisionModel prec, const LocalExclusionMask &lem,
                           const PPITable &pairs_tbl,
                           CellGrid<double, llint, double, double4> *cg, TileManager *tlmn,
                           ScoreCard *sc, MolecularMechanicsControls *mmctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const CoreKlManager &launcher, const double clash_distance,
                           const double clash_ratio) {

  // Obtain the correct kernel launch parameters
  const ClashResponse mitigate_clash = (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) ?
                                       ClashResponse::FORGIVE : ClashResponse::NONE;
  const AtomGraphSynthesis *poly_ag = cg->getTopologySynthesisPointer();
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  const LocalExclusionMaskReader lemr = lem.data(devc);
  TilePlan tlpn = tlmn->data(devc);
  ScoreCardWriter scw = sc->data(devc);
  CellGridWriter<double, llint, double, double4> cgw = cg->data(devc);
  const TinyBoxPresence has_tiny_box = cg->getTinyBoxPresence();
  const int2 bt_tp = launcher.getPMEPairsKernelDims(PrecisionModel::DOUBLE, prec, 
                                                    NeighborListKind::MONO, has_tiny_box,
                                                    eval_frc, eval_nrg, mitigate_clash);
  const PsSynthesisBorders sysbrd = cg->getUnitCellTransforms(devc);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag->getDoublePrecisionNonbondedKit(devc);
      const PPIKit<double, double4> nrg_tab = pairs_tbl.dpData();
      MMControlKit<double> ctrl = mmctrl->dpData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw, &tlpn, &scw, &ctrl, eval_frc, eval_nrg,
                       bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag->getSinglePrecisionNonbondedKit(devc);
      const PPIKit<float, float4> nrg_tab = pairs_tbl.spData();
      MMControlKit<float> ctrl = mmctrl->spData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw, &tlpn, &scw, &ctrl, eval_frc, eval_nrg,
                       bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const PrecisionModel prec, const LocalExclusionMask &lem,
                           const PPITable &pairs_tbl,
                           CellGrid<double, llint, double, double4> *cg_qq,
                           CellGrid<double, llint, double, double4> *cg_lj, TileManager *tlmn,
                           ScoreCard *sc, MolecularMechanicsControls *mmctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const CoreKlManager &launcher, const double clash_distance,
                           const double clash_ratio) {

  // Obtain the correct kernel launch parameters
  const ClashResponse mitigate_clash = (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) ?
                                       ClashResponse::FORGIVE : ClashResponse::NONE;
  const AtomGraphSynthesis *poly_ag = cg_qq->getTopologySynthesisPointer();
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  const LocalExclusionMaskReader lemr = lem.data(devc);
  TilePlan tlpn = tlmn->data(devc);
  ScoreCardWriter scw = sc->data(devc);
  CellGridWriter<double, llint, double, double4> cgw_qq = cg_qq->data(devc);
  CellGridWriter<double, llint, double, double4> cgw_lj = cg_lj->data(devc);
  const TinyBoxPresence has_tiny_box = (cg_qq->getTinyBoxPresence() == TinyBoxPresence::YES ||
                                        cg_lj->getTinyBoxPresence() == TinyBoxPresence::YES) ?
                                       TinyBoxPresence::YES : TinyBoxPresence::NO;
  const int2 bt_tp = launcher.getPMEPairsKernelDims(PrecisionModel::DOUBLE, prec, 
                                                    NeighborListKind::DUAL, has_tiny_box,
                                                    eval_frc, eval_nrg, mitigate_clash);
  const PsSynthesisBorders sysbrd = cg_qq->getUnitCellTransforms(devc);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag->getDoublePrecisionNonbondedKit(devc);
      const PPIKit<double, double4> nrg_tab = pairs_tbl.dpData();
      MMControlKit<double> ctrl = mmctrl->dpData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl,
                       eval_frc, eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag->getSinglePrecisionNonbondedKit(devc);
      const PPIKit<float, float4> nrg_tab = pairs_tbl.spData();
      MMControlKit<float> ctrl = mmctrl->spData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl,
                       eval_frc, eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const PrecisionModel prec, const LocalExclusionMask &lem,
                           const PPITable &pairs_tbl, CellGrid<float, int, float, float4> *cg,
                           TileManager *tlmn, ScoreCard *sc, MolecularMechanicsControls *mmctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const CoreKlManager &launcher, const double clash_distance,
                           const double clash_ratio) {

  // Obtain the correct kernel launch parameters
  const ClashResponse mitigate_clash = (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) ?
                                       ClashResponse::FORGIVE : ClashResponse::NONE;
  const AtomGraphSynthesis *poly_ag = cg->getTopologySynthesisPointer();
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  const LocalExclusionMaskReader lemr = lem.data(devc);
  TilePlan tlpn = tlmn->data(devc);
  ScoreCardWriter scw = sc->data(devc);
  CellGridWriter<float, int, float, float4> cgw = cg->data(devc);
  const TinyBoxPresence has_tiny_box = cg->getTinyBoxPresence();
  const int2 bt_tp = launcher.getPMEPairsKernelDims(PrecisionModel::SINGLE, prec,
                                                    NeighborListKind::MONO,
                                                    cg->getTinyBoxPresence(), eval_frc, eval_nrg,
                                                    mitigate_clash);
  const PsSynthesisBorders sysbrd = cg->getUnitCellTransforms(devc);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag->getDoublePrecisionNonbondedKit(devc);
      const PPIKit<double, double4> nrg_tab = pairs_tbl.dpData();
      MMControlKit<double> ctrl = mmctrl->dpData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw, &tlpn, &scw, &ctrl, eval_frc, eval_nrg,
                       bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag->getSinglePrecisionNonbondedKit(devc);
      const PPIKit<float, float4> nrg_tab = pairs_tbl.spData();
      MMControlKit<float> ctrl = mmctrl->spData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw, &tlpn, &scw, &ctrl, eval_frc, eval_nrg,
                       bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const PrecisionModel prec, const LocalExclusionMask &lem,
                           const PPITable &pairs_tbl, CellGrid<float, int, float, float4> *cg_qq,
                           CellGrid<float, int, float, float4> *cg_lj, TileManager *tlmn,
                           ScoreCard *sc, MolecularMechanicsControls *mmctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const CoreKlManager &launcher, const double clash_distance,
                           const double clash_ratio) {

  // Obtain the correct kernel launch parameters
  const ClashResponse mitigate_clash = (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) ?
                                       ClashResponse::FORGIVE : ClashResponse::NONE;
  const AtomGraphSynthesis *poly_ag = cg_qq->getTopologySynthesisPointer();
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  const LocalExclusionMaskReader lemr = lem.data(devc);
  TilePlan tlpn = tlmn->data(devc);
  ScoreCardWriter scw = sc->data(devc);
  CellGridWriter<float, int, float, float4> cgw_qq = cg_qq->data(devc);
  CellGridWriter<float, int, float, float4> cgw_lj = cg_lj->data(devc);
  const TinyBoxPresence has_tiny_box = (cg_qq->getTinyBoxPresence() == TinyBoxPresence::YES ||
                                        cg_lj->getTinyBoxPresence() == TinyBoxPresence::YES) ?
                                       TinyBoxPresence::YES : TinyBoxPresence::NO;
  const int2 bt_tp = launcher.getPMEPairsKernelDims(PrecisionModel::SINGLE, prec,
                                                    NeighborListKind::DUAL, has_tiny_box,
                                                    eval_frc, eval_nrg, mitigate_clash);
  const PsSynthesisBorders sysbrd = cg_qq->getUnitCellTransforms(devc);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag->getDoublePrecisionNonbondedKit(devc);
      const PPIKit<double, double4> nrg_tab = pairs_tbl.dpData();
      MMControlKit<double> ctrl = mmctrl->dpData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl,
                       eval_frc, eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag->getSinglePrecisionNonbondedKit(devc);
      const PPIKit<float, float4> nrg_tab = pairs_tbl.spData();
      MMControlKit<float> ctrl = mmctrl->spData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl,
                       eval_frc, eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  }
}

} // namespace energy
} // namespace stormm
