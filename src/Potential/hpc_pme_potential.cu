// -*-c++-*-
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Constants/hpc_bounds.h"
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

#define PMENB_BLOCK_MULTIPLICITY 4

// Single-precision tile evaluation
#define TCALC float
#  define TCALC2 float2
#  define TCALC4 float4
#  define TCALC_IS_SINGLE
#  define TCOORD_IS_REAL

// Begin with a float32_t coordinate representation, natural for the float32_t arithmetic mode.
#  define TCOORD  float
#  define TCOORD4 float4
#  define TACC    int

// Compile the kernels with or without energy and force computations, dual neighbor lists,
// provisions for small box sizes, and clash forgiveness.
#  define COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define DUAL_GRIDS
#        define TINY_BOX
#          define CLASH_FORGIVENESS
#            define PMENB_WARPS_PER_BLOCK 8
#            define KERNEL_NAME kffTowerPlateFEDualTinyNonClash
#              include "tower_plate_pairs.cui"
#            undef KERNEL_NAME
#            undef PMENB_WARPS_PER_BLOCK
#            define PMENB_WARPS_PER_BLOCK 8
#            define KERNEL_NAME kffTowerTowerFEDualTinyNonClash
#              include "tower_tower_pairs.cui"
#            undef KERNEL_NAME
#            undef PMENB_WARPS_PER_BLOCK
#          undef CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kffTowerPlateFEDualTiny
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kffTowerTowerFEDualTiny
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef SMALL_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kffTowerPlateFEDualNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kffTowerTowerFEDualNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerPlateFEDual
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerTowerFEDual
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kffTowerPlateFETinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kffTowerTowerFETinyNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerPlateFETiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerTowerFETiny
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerPlateFENonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerTowerFENonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kffTowerPlateFE
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kffTowerTowerFE
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kffTowerPlateFXDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kffTowerTowerFXDualTinyNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerPlateFXDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerTowerFXDualTiny
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerPlateFXDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerTowerFXDualNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 10
#      define KERNEL_NAME kffTowerPlateFXDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 9
#      define KERNEL_NAME kffTowerTowerFXDual
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerPlateFXTinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerTowerFXTinyNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kffTowerPlateFXTiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kffTowerTowerFXTiny
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef SMALL_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kffTowerPlateFXNonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kffTowerTowerFXNonClash
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 10
#    define KERNEL_NAME kffTowerPlateFX
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#    define PMENB_WARPS_PER_BLOCK 9
#    define KERNEL_NAME kffTowerTowerFX
#      include "tower_tower_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kffTowerPlateXEDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kffTowerTowerXEDualTinyNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerPlateXEDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerTowerXEDualTiny
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerPlateXEDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerTowerXEDualNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 10
#      define KERNEL_NAME kffTowerPlateXEDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 9
#      define KERNEL_NAME kffTowerTowerXEDual
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerPlateXETinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kffTowerTowerXETinyNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kffTowerPlateXETiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kffTowerTowerXETiny
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef SMALL_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kffTowerPlateXENonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kffTowerTowerXENonClash
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 8
#    define KERNEL_NAME kffTowerPlateXE
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#    define PMENB_WARPS_PER_BLOCK 8
#    define KERNEL_NAME kffTowerTowerXE
#      include "tower_tower_pairs.cui"
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
#            define PMENB_WARPS_PER_BLOCK 8
#            define KERNEL_NAME kdfTowerPlateFEDualTinyNonClash
#              include "tower_plate_pairs.cui"
#            undef KERNEL_NAME
#            undef PMENB_WARPS_PER_BLOCK
#            define PMENB_WARPS_PER_BLOCK 8
#            define KERNEL_NAME kdfTowerTowerFEDualTinyNonClash
#              include "tower_tower_pairs.cui"
#            undef KERNEL_NAME
#            undef PMENB_WARPS_PER_BLOCK
#          undef CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kdfTowerPlateFEDualTiny
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kdfTowerTowerFEDualTiny
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef SMALL_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kdfTowerPlateFEDualNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kdfTowerTowerFEDualNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerPlateFEDual
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerTowerFEDual
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kdfTowerPlateFETinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kdfTowerTowerFETinyNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerPlateFETiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerTowerFETiny
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerPlateFENonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerTowerFENonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kdfTowerPlateFE
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kdfTowerTowerFE
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kdfTowerPlateFXDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kdfTowerTowerFXDualTinyNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerPlateFXDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerTowerFXDualTiny
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerPlateFXDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerTowerFXDualNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kdfTowerPlateFXDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kdfTowerTowerFXDual
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerPlateFXTinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerTowerFXTinyNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kdfTowerPlateFXTiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kdfTowerTowerFXTiny
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef SMALL_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kdfTowerPlateFXNonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kdfTowerTowerFXNonClash
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 8
#    define KERNEL_NAME kdfTowerPlateFX
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#    define PMENB_WARPS_PER_BLOCK 8
#    define KERNEL_NAME kdfTowerTowerFX
#      include "tower_tower_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kdfTowerPlateXEDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 8
#          define KERNEL_NAME kdfTowerTowerXEDualTinyNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerPlateXEDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerTowerXEDualTiny
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerPlateXEDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerTowerXEDualNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kdfTowerPlateXEDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kdfTowerTowerXEDual
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerPlateXETinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 8
#        define KERNEL_NAME kdfTowerTowerXETinyNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kdfTowerPlateXETiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kdfTowerTowerXETiny
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef SMALL_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kdfTowerPlateXENonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 8
#      define KERNEL_NAME kdfTowerTowerXENonClash
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 8
#    define KERNEL_NAME kdfTowerPlateXE
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#    define PMENB_WARPS_PER_BLOCK 8
#    define KERNEL_NAME kdfTowerTowerXE
#      include "tower_tower_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#  undef COMPUTE_ENERGY

#  undef TCOORD
#  undef TCOORD4
#  undef TACC
#  undef TCOORD_IS_LONG
  
#  undef TCOORD_IS_REAL
#  undef TCALC_IS_SINGLE
#  undef TCALC2
#  undef TCALC4
#undef TCALC
  
// Double-precision tile evaluation
#define TCALC double
#  define TCALC2 double2
#  define TCALC4 double4
#  define TCOORD_IS_REAL

// Begin with a float64_t coordinate representation, appropriate for the float64_t arithmetic mode.
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
#            define PMENB_WARPS_PER_BLOCK 6
#            define KERNEL_NAME kddTowerPlateFEDualTinyNonClash
#              include "tower_plate_pairs.cui"
#            undef KERNEL_NAME
#            undef PMENB_WARPS_PER_BLOCK
#            define PMENB_WARPS_PER_BLOCK 6
#            define KERNEL_NAME kddTowerTowerFEDualTinyNonClash
#              include "tower_tower_pairs.cui"
#            undef KERNEL_NAME
#            undef PMENB_WARPS_PER_BLOCK
#          undef CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kddTowerPlateFEDualTiny
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kddTowerTowerFEDualTiny
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef SMALL_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kddTowerPlateFEDualNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kddTowerTowerFEDualNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerPlateFEDual
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerTowerFEDual
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kddTowerPlateFETinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kddTowerTowerFETinyNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerPlateFETiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerTowerFETiny
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerPlateFENonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerTowerFENonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kddTowerPlateFE
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kddTowerTowerFE
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kddTowerPlateFXDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kddTowerTowerFXDualTinyNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerPlateFXDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerTowerFXDualTiny
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerPlateFXDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerTowerFXDualNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kddTowerPlateFXDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kddTowerTowerFXDual
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerPlateFXTinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerTowerFXTinyNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kddTowerPlateFXTiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kddTowerTowerFXTiny
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef SMALL_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kddTowerPlateFXNonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kddTowerTowerFXNonClash
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 6
#    define KERNEL_NAME kddTowerPlateFX
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#    define PMENB_WARPS_PER_BLOCK 6
#    define KERNEL_NAME kddTowerTowerFX
#      include "tower_tower_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kddTowerPlateXEDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kddTowerTowerXEDualTinyNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerPlateXEDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerTowerXEDualTiny
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerPlateXEDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerTowerXEDualNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kddTowerPlateXEDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kddTowerTowerXEDual
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerPlateXETinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kddTowerTowerXETinyNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kddTowerPlateXETiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kddTowerTowerXETiny
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef SMALL_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kddTowerPlateXENonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kddTowerTowerXENonClash
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 6
#    define KERNEL_NAME kddTowerPlateXE
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#    define PMENB_WARPS_PER_BLOCK 6
#    define KERNEL_NAME kddTowerTowerXE
#      include "tower_tower_pairs.cui"
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
#            define PMENB_WARPS_PER_BLOCK 6
#            define KERNEL_NAME kfdTowerPlateFEDualTinyNonClash
#              include "tower_plate_pairs.cui"
#            undef KERNEL_NAME
#            undef PMENB_WARPS_PER_BLOCK
#            define PMENB_WARPS_PER_BLOCK 6
#            define KERNEL_NAME kfdTowerTowerFEDualTinyNonClash
#              include "tower_tower_pairs.cui"
#            undef KERNEL_NAME
#            undef PMENB_WARPS_PER_BLOCK
#          undef CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kfdTowerPlateFEDualTiny
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kfdTowerTowerFEDualTiny
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef SMALL_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kfdTowerPlateFEDualNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kfdTowerTowerFEDualNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerPlateFEDual
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerTowerFEDual
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kfdTowerPlateFETinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kfdTowerTowerFETinyNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerPlateFETiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerTowerFETiny
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerPlateFENonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerTowerFENonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kfdTowerPlateFE
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kfdTowerTowerFE
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kfdTowerPlateFXDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kfdTowerTowerFXDualTinyNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerPlateFXDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerTowerFXDualTiny
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerPlateFXDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerTowerFXDualNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kfdTowerPlateFXDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kfdTowerTowerFXDual
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerPlateFXTinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerTowerFXTinyNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kfdTowerPlateFXTiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kfdTowerTowerFXTiny
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef SMALL_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kfdTowerPlateFXNonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kfdTowerTowerFXNonClash
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 6
#    define KERNEL_NAME kfdTowerPlateFX
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#    define PMENB_WARPS_PER_BLOCK 6
#    define KERNEL_NAME kfdTowerTowerFX
#      include "tower_tower_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kfdTowerPlateXEDualTinyNonClash
#            include "tower_plate_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#          define PMENB_WARPS_PER_BLOCK 6
#          define KERNEL_NAME kfdTowerTowerXEDualTinyNonClash
#            include "tower_tower_pairs.cui"
#          undef KERNEL_NAME
#          undef PMENB_WARPS_PER_BLOCK
#        undef CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerPlateXEDualTiny
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerTowerXEDualTiny
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef SMALL_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerPlateXEDualNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerTowerXEDualNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kfdTowerPlateXEDual
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kfdTowerTowerXEDual
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef DUAL_GRIDS
#    define TINY_BOX
#      define CLASH_FORGIVENESS
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerPlateXETinyNonClash
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#        define PMENB_WARPS_PER_BLOCK 6
#        define KERNEL_NAME kfdTowerTowerXETinyNonClash
#          include "tower_tower_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kfdTowerPlateXETiny
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kfdTowerTowerXETiny
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef SMALL_BOX
#    define CLASH_FORGIVENESS
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kfdTowerPlateXENonClash
#        include "tower_plate_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#      define PMENB_WARPS_PER_BLOCK 6
#      define KERNEL_NAME kfdTowerTowerXENonClash
#        include "tower_tower_pairs.cui"
#      undef KERNEL_NAME
#      undef PMENB_WARPS_PER_BLOCK
#    undef CLASH_FORGIVENESS
#    define PMENB_WARPS_PER_BLOCK 6
#    define KERNEL_NAME kfdTowerPlateXE
#      include "tower_plate_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#    define PMENB_WARPS_PER_BLOCK 6
#    define KERNEL_NAME kfdTowerTowerXE
#      include "tower_tower_pairs.cui"
#    undef KERNEL_NAME
#    undef PMENB_WARPS_PER_BLOCK
#  undef COMPUTE_ENERGY

#  undef TCOORD
#  undef TCOORD4
#  undef TACC
  
#  undef TCOORD_IS_REAL
#  undef TCALC2
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
                                                          const ClashResponse clash_handling,
                                                          const PairStance span) {

  // As with other kernel querying functions, the kernel manager calling this function will have
  // specifications of the GPU in use.  It is the overall thread occupancy and multiplicity of each
  // kernel that this function must return.
  cudaFuncAttributes result;
  cudaError_t cfa;
  switch (span) {
  case PairStance::TOWER_PLATE:
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
    break;
  case PairStance::TOWER_TOWER:
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
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFEDualTiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFEDual);
                  break;
                }
                break;
              case NeighborListKind::MONO:
               switch (has_tiny_box) {
                 case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFETiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFE);
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
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFXDualTiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFXDual);
                  break;
                }
                break;
              case NeighborListKind::MONO:
                switch (has_tiny_box) {
                case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFXTiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFX);
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
                cfa = cudaFuncGetAttributes(&result, kddTowerTowerXEDualTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kddTowerTowerXEDual);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kddTowerTowerXETiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kddTowerTowerXE);
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
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFEDualTiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFEDual);
                  break;
                }
                break;
              case NeighborListKind::MONO:
                switch (has_tiny_box) {
                case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFETiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFE);
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
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFXDualTiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFXDual);
                  break;
                }
                break;
              case NeighborListKind::MONO:
                switch (has_tiny_box) {
                case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFXTiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFX);
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
                cfa = cudaFuncGetAttributes(&result, kdfTowerTowerXEDualTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kdfTowerTowerXEDual);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kdfTowerTowerXETiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kdfTowerTowerXE);
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
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFEDualTiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFEDual);
                  break;
                }
                break;
              case NeighborListKind::MONO:
                switch (has_tiny_box) {
                case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFETiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFE);
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
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFXDualTiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFXDual);
                  break;
                }
                break;
              case NeighborListKind::MONO:
                switch (has_tiny_box) {
                case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFXTiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFX);
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
              cfa = cudaFuncGetAttributes(&result, kfdTowerTowerXEDualTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kfdTowerTowerXEDual);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kfdTowerTowerXETiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kfdTowerTowerXE);
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
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFEDualTiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFEDual);
                  break;
                }
                break;
              case NeighborListKind::MONO:
                switch (has_tiny_box) {
                case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFETiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFE);
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
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFXDualTiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFXDual);
                  break;
                }
                break;
              case NeighborListKind::MONO:
                switch (has_tiny_box) {
                case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFXTiny);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFX);
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
                cfa = cudaFuncGetAttributes(&result, kffTowerTowerXEDualTiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kffTowerTowerXEDual);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kffTowerTowerXETiny);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kffTowerTowerXE);
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
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFEDualTinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFEDualNonClash);
                  break;
                }
                break;
              case NeighborListKind::MONO:
                switch (has_tiny_box) {
                case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFETinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFENonClash);
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
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFXDualTinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFXDualNonClash);
                  break;
                }
                break;
              case NeighborListKind::MONO:
                switch (has_tiny_box) {
                case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFXTinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kddTowerTowerFXNonClash);
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
                cfa = cudaFuncGetAttributes(&result, kddTowerTowerXEDualTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kddTowerTowerXEDualNonClash);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kddTowerTowerXETinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kddTowerTowerXENonClash);
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
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFEDualTinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFEDualNonClash);
                  break;
                }
                break;
              case NeighborListKind::MONO:
                switch (has_tiny_box) {
                case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFETinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFENonClash);
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
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFXDualTinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFXDualNonClash);
                  break;
                }
                break;
              case NeighborListKind::MONO:
                switch (has_tiny_box) {
                case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFXTinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kdfTowerTowerFXNonClash);
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
                cfa = cudaFuncGetAttributes(&result, kdfTowerTowerXEDualTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kdfTowerTowerXEDualNonClash);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kdfTowerTowerXETinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kdfTowerTowerXENonClash);
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
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFEDualTinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFEDualNonClash);
                  break;
                }
                break;
              case NeighborListKind::MONO:
                switch (has_tiny_box) {
                case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFETinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFENonClash);
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
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFXDualTinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFXDualNonClash);
                  break;
                }
                break;
              case NeighborListKind::MONO:
                switch (has_tiny_box) {
                case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFXTinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kfdTowerTowerFXNonClash);
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
                cfa = cudaFuncGetAttributes(&result, kfdTowerTowerXEDualTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kfdTowerTowerXEDualNonClash);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kfdTowerTowerXETinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kfdTowerTowerXENonClash);
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
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFEDualTinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFEDualNonClash);
                  break;
                }
                break;
              case NeighborListKind::MONO:
                switch (has_tiny_box) {
                case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFETinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFENonClash);
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
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFXDualTinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFXDualNonClash);
                  break;
                }
                break;
              case NeighborListKind::MONO:
                switch (has_tiny_box) {
                case TinyBoxPresence::YES:
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFXTinyNonClash);
                  break;
                case TinyBoxPresence::NO:
                  cfa = cudaFuncGetAttributes(&result, kffTowerTowerFXNonClash);
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
                cfa = cudaFuncGetAttributes(&result, kffTowerTowerXEDualTinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kffTowerTowerXEDualNonClash);
                break;
              }
              break;
            case NeighborListKind::MONO:
              switch (has_tiny_box) {
              case TinyBoxPresence::YES:
                cfa = cudaFuncGetAttributes(&result, kffTowerTowerXETinyNonClash);
                break;
              case TinyBoxPresence::NO:
                cfa = cudaFuncGetAttributes(&result, kffTowerTowerXENonClash);
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
    std::string stance_name;
    switch (span) {
    case PairStance::TOWER_PLATE:
      stance_name += "TowerPlate";
      break;
    case PairStance::TOWER_TOWER:
      stance_name = "TowerTower";
      break;
    }
    error_message += stance_name;
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
                           const TinyBoxPresence has_tiny_box, const int2 bt_tp, const int2 bt_tt,
                           const double clash_distance, const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // DOUBLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kddTowerPlateFETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw, *ctrl);
          kddTowerTowerFETinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kddTowerPlateFXTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw,
                                                            *ctrl);
          kddTowerTowerFXTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw,
                                                            *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kddTowerPlateXETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw, *cgw,
                                                          *ctrl);
        kddTowerTowerXETinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw, *cgw,
                                                          *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kddTowerPlateFENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *scw, *cgw,
                                                        *ctrl);
          kddTowerTowerFENonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *scw, *cgw,
                                                        *ctrl);
          break;
        case EvaluateEnergy::NO:
          kddTowerPlateFXNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *cgw, *ctrl);
          kddTowerTowerFXNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kddTowerPlateXENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *scw, *cgw,
                                                      *ctrl);
        kddTowerTowerXENonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *scw, *cgw,
                                                      *ctrl);
        break;
      }
      break;
    }
  }
  else {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kddTowerPlateFETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                    *ctrl);
          kddTowerTowerFETiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                    *ctrl);
          break;
        case EvaluateEnergy::NO:
          kddTowerPlateFXTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          kddTowerTowerFXTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kddTowerPlateXETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                  *ctrl);
        kddTowerTowerXETiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                  *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kddTowerPlateFE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
          kddTowerTowerFE<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kddTowerPlateFX<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          kddTowerTowerFX<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kddTowerPlateXE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
        kddTowerTowerXE<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
        break;
      }
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
                           const TinyBoxPresence has_tiny_box, const int2 bt_tp, const int2 bt_tt,
                           const double clash_distance, const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // DOUBLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfdTowerPlateFETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw, *ctrl);
          kfdTowerTowerFETinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfdTowerPlateFXTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw,
                                                            *ctrl);
          kfdTowerTowerFXTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw,
                                                            *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfdTowerPlateXETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw, *cgw,
                                                          *ctrl);
        kfdTowerTowerXETinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw, *cgw,
                                                          *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfdTowerPlateFENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *scw, *cgw,
                                                        *ctrl);
          kfdTowerTowerFENonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *scw, *cgw,
                                                        *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfdTowerPlateFXNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *cgw, *ctrl);
          kfdTowerTowerFXNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfdTowerPlateXENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *scw, *cgw,
                                                      *ctrl);
        kfdTowerTowerXENonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *scw, *cgw,
                                                      *ctrl);
        break;
      }
      break;
    }
  }
  else {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfdTowerPlateFETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                    *ctrl);
          kfdTowerTowerFETiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                    *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfdTowerPlateFXTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          kfdTowerTowerFXTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfdTowerPlateXETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                  *ctrl);
        kfdTowerTowerXETiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                  *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfdTowerPlateFE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
          kfdTowerTowerFE<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfdTowerPlateFX<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          kfdTowerTowerFX<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfdTowerPlateXE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
        kfdTowerTowerXE<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
        break;
      }
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
                           const TinyBoxPresence has_tiny_box, const int2 bt_tp, const int2 bt_tt,
                           const double clash_distance, const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // DOUBLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kddTowerPlateFEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio, *scw,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          kddTowerTowerFEDualTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio, *scw,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kddTowerPlateFXDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          kddTowerTowerFXDualTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kddTowerPlateXEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              clash_distance, clash_ratio, *scw,
                                                              *cgw_qq, *cgw_lj, *ctrl);
        kddTowerTowerXEDualTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              clash_distance, clash_ratio, *scw,
                                                              *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kddTowerPlateFEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw_qq, *cgw_lj, *ctrl);
          kddTowerTowerFEDualNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kddTowerPlateFXDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw_qq,
                                                            *cgw_lj, *ctrl);
          kddTowerTowerFXDualNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw_qq,
                                                            *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kddTowerPlateXEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw_qq, *cgw_lj, *ctrl);
        kddTowerTowerXEDualNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    }
  }
  else {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kddTowerPlateFEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                        *cgw_qq, *cgw_lj, *ctrl);
          kddTowerTowerFEDualTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                        *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kddTowerPlateFXDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                        *cgw_lj, *ctrl);
          kddTowerTowerFXDualTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                        *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kddTowerPlateXEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        kddTowerTowerXEDualTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kddTowerPlateFEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          kddTowerTowerFEDual<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kddTowerPlateFXDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          kddTowerTowerFXDual<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kddTowerPlateXEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                  *cgw_lj, *ctrl);
        kddTowerTowerXEDual<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                  *cgw_lj, *ctrl);
        break;
      }
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
                           const TinyBoxPresence has_tiny_box, const int2 bt_tp, const int2 bt_tt,
                           const double clash_distance, const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // DOUBLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfdTowerPlateFEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio, *scw,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          kfdTowerTowerFEDualTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio, *scw,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfdTowerPlateFXDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          kfdTowerTowerFXDualTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfdTowerPlateXEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              clash_distance, clash_ratio, *scw,
                                                              *cgw_qq, *cgw_lj, *ctrl);
        kfdTowerTowerXEDualTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              clash_distance, clash_ratio, *scw,
                                                              *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfdTowerPlateFEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw_qq, *cgw_lj, *ctrl);
          kfdTowerTowerFEDualNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfdTowerPlateFXDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw_qq,
                                                            *cgw_lj, *ctrl);
          kfdTowerTowerFXDualNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw_qq,
                                                            *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfdTowerPlateXEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw_qq, *cgw_lj, *ctrl);
        kfdTowerTowerXEDualNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    }
  }
  else {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfdTowerPlateFEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                        *cgw_qq, *cgw_lj, *ctrl);
          kfdTowerTowerFEDualTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                        *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfdTowerPlateFXDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                        *cgw_lj, *ctrl);
          kfdTowerTowerFXDualTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                        *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfdTowerPlateXEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        kfdTowerTowerXEDualTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kfdTowerPlateFEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          kfdTowerTowerFEDual<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kfdTowerPlateFXDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          kfdTowerTowerFXDual<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kfdTowerPlateXEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                  *cgw_lj, *ctrl);
        kfdTowerTowerXEDual<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                  *cgw_lj, *ctrl);
        break;
      }
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
                           const TinyBoxPresence has_tiny_box, const int2 bt_tp, const int2 bt_tt,
                           const double clash_distance, const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // SINGLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kffTowerPlateFETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw, *ctrl);
          kffTowerTowerFETinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kffTowerPlateFXTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw,
                                                            *ctrl);
          kffTowerTowerFXTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw,
                                                            *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kffTowerPlateXETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw, *cgw,
                                                          *ctrl);
        kffTowerTowerXETinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw, *cgw,
                                                          *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kffTowerPlateFENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *scw, *cgw,
                                                        *ctrl);
          kffTowerTowerFENonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *scw, *cgw,
                                                        *ctrl);
          break;
        case EvaluateEnergy::NO:
          kffTowerPlateFXNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *cgw, *ctrl);
          kffTowerTowerFXNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kffTowerPlateXENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *scw, *cgw,
                                                      *ctrl);
        kffTowerTowerXENonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *scw, *cgw,
                                                      *ctrl);
        break;
      }
      break;
    }
  }
  else {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kffTowerPlateFETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                    *ctrl);
          kffTowerTowerFETiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                    *ctrl);
          break;
        case EvaluateEnergy::NO:
          kffTowerPlateFXTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          kffTowerTowerFXTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kffTowerPlateXETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                  *ctrl);
        kffTowerTowerXETiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                  *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kffTowerPlateFE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
          kffTowerTowerFE<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kffTowerPlateFX<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          kffTowerTowerFX<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kffTowerPlateXE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
        kffTowerTowerXE<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
        break;
      }
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
                           const TinyBoxPresence has_tiny_box, const int2 bt_tp, const int2 bt_tt,
                           const double clash_distance, const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // SINGLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdfTowerPlateFETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw, *ctrl);
          kdfTowerTowerFETinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdfTowerPlateFXTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw,
                                                            *ctrl);
          kdfTowerTowerFXTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw,
                                                            *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdfTowerPlateXETinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw, *cgw,
                                                          *ctrl);
        kdfTowerTowerXETinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw, *cgw,
                                                          *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdfTowerPlateFENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *scw, *cgw,
                                                        *ctrl);
          kdfTowerTowerFENonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *scw, *cgw,
                                                        *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdfTowerPlateFXNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *cgw, *ctrl);
          kdfTowerTowerFXNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                        clash_distance, clash_ratio, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdfTowerPlateXENonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *scw, *cgw,
                                                      *ctrl);
        kdfTowerTowerXENonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                      clash_distance, clash_ratio, *scw, *cgw,
                                                      *ctrl);
        break;
      }
      break;
    }
  }
  else {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdfTowerPlateFETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                    *ctrl);
          kdfTowerTowerFETiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                    *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdfTowerPlateFXTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          kdfTowerTowerFXTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdfTowerPlateXETiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                  *ctrl);
        kdfTowerTowerXETiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw,
                                                  *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdfTowerPlateFE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
          kdfTowerTowerFE<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdfTowerPlateFX<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          kdfTowerTowerFX<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdfTowerPlateXE<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
        kdfTowerTowerXE<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw, *ctrl);
        break;
      }
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
                           const TinyBoxPresence has_tiny_box, const int2 bt_tp, const int2 bt_tt,
                           const double clash_distance, const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // SINGLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kffTowerPlateFEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio, *scw,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          kffTowerTowerFEDualTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio, *scw,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kffTowerPlateFXDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          kffTowerTowerFXDualTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kffTowerPlateXEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              clash_distance, clash_ratio, *scw,
                                                              *cgw_qq, *cgw_lj, *ctrl);
        kffTowerTowerXEDualTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              clash_distance, clash_ratio, *scw,
                                                              *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kffTowerPlateFEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw_qq, *cgw_lj, *ctrl);
          kffTowerTowerFEDualNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kffTowerPlateFXDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw_qq,
                                                            *cgw_lj, *ctrl);
          kffTowerTowerFXDualNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw_qq,
                                                            *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kffTowerPlateXEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw_qq, *cgw_lj, *ctrl);
        kffTowerTowerXEDualNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    }
  }
  else {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kffTowerPlateFEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                        *cgw_qq, *cgw_lj, *ctrl);
          kffTowerTowerFEDualTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                        *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kffTowerPlateFXDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                        *cgw_lj, *ctrl);
          kffTowerTowerFXDualTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                        *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kffTowerPlateXEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        kffTowerTowerXEDualTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kffTowerPlateFEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          kffTowerTowerFEDual<<<1, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kffTowerPlateFXDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          kffTowerTowerFXDual<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kffTowerPlateXEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                  *cgw_lj, *ctrl);
        kffTowerTowerXEDual<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                  *cgw_lj, *ctrl);
        break;
      }
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
                           const TinyBoxPresence has_tiny_box, const int2 bt_tp, const int2 bt_tt,
                           const double clash_distance, const double clash_ratio) {

  // Clash dampening is detected by the values of the parameters rather than an explicit
  // enumeration.  All launches from this overloaded variant of launchPMEPairs() will invoke
  // SINGLE precision calculations, coordinates, and parameter sets.  The single neighbor list grid
  // also implies a branch of the kernels taking unified neighbor lists.  The neighbor list itself
  // will be queried for a condition to see whether there is a "tiny" simulation box with 4 cells
  // along any one dimension.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdfTowerPlateFEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio, *scw,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          kdfTowerTowerFEDualTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio, *scw,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdfTowerPlateFXDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          kdfTowerTowerFXDualTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                                clash_distance, clash_ratio,
                                                                *cgw_qq, *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdfTowerPlateXEDualTinyNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              clash_distance, clash_ratio, *scw,
                                                              *cgw_qq, *cgw_lj, *ctrl);
        kdfTowerTowerXEDualTinyNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                              clash_distance, clash_ratio, *scw,
                                                              *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdfTowerPlateFEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw_qq, *cgw_lj, *ctrl);
          kdfTowerTowerFEDualNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *scw,
                                                            *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdfTowerPlateFXDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw_qq,
                                                            *cgw_lj, *ctrl);
          kdfTowerTowerFXDualNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                            clash_distance, clash_ratio, *cgw_qq,
                                                            *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdfTowerPlateXEDualNonClash<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw_qq, *cgw_lj, *ctrl);
        kdfTowerTowerXEDualNonClash<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab,
                                                          clash_distance, clash_ratio, *scw,
                                                          *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    }
  }
  else {
    switch (has_tiny_box) {
    case TinyBoxPresence::YES:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdfTowerPlateFEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                        *cgw_qq, *cgw_lj, *ctrl);
          kdfTowerTowerFEDualTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                        *cgw_qq, *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdfTowerPlateFXDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                        *cgw_lj, *ctrl);
          kdfTowerTowerFXDualTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                        *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdfTowerPlateXEDualTiny<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        kdfTowerTowerXEDualTiny<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw,
                                                      *cgw_qq, *cgw_lj, *ctrl);
        break;
      }
      break;
    case TinyBoxPresence::NO:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          kdfTowerPlateFEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          kdfTowerTowerFEDual<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          break;
        case EvaluateEnergy::NO:
          kdfTowerPlateFXDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          kdfTowerTowerFXDual<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *cgw_qq,
                                                    *cgw_lj, *ctrl);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdfTowerPlateXEDual<<<bt_tp.x, bt_tp.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                  *cgw_lj, *ctrl);
        kdfTowerTowerXEDual<<<bt_tt.x, bt_tt.y>>>(poly_nbk, lemr, *tlpn, nrg_tab, *scw, *cgw_qq,
                                                  *cgw_lj, *ctrl);
        break;
      }
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
                                                    eval_frc, eval_nrg, mitigate_clash,
                                                    PairStance::TOWER_PLATE);
  const int2 bt_tt = launcher.getPMEPairsKernelDims(PrecisionModel::DOUBLE, prec, 
                                                    NeighborListKind::MONO, has_tiny_box,
                                                    eval_frc, eval_nrg, mitigate_clash,
                                                    PairStance::TOWER_TOWER);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag->getDoublePrecisionNonbondedKit(devc);
      const PPIKit<double, double4> nrg_tab = pairs_tbl.dpData();
      MMControlKit<double> ctrl = mmctrl->dpData(devc);
      launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw, &tlpn, &scw, &ctrl, eval_frc, eval_nrg,
                     has_tiny_box, bt_tp, bt_tt, clash_distance, clash_ratio);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag->getSinglePrecisionNonbondedKit(devc);
      const PPIKit<float, float4> nrg_tab = pairs_tbl.spData();
      MMControlKit<float> ctrl = mmctrl->spData(devc);
      launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw, &tlpn, &scw, &ctrl, eval_frc, eval_nrg,
                     has_tiny_box, bt_tp, bt_tt, clash_distance, clash_ratio);
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
                                                    eval_frc, eval_nrg, mitigate_clash,
                                                    PairStance::TOWER_PLATE);
  const int2 bt_tt = launcher.getPMEPairsKernelDims(PrecisionModel::DOUBLE, prec,
                                                    NeighborListKind::DUAL, has_tiny_box,
                                                    eval_frc, eval_nrg, mitigate_clash,
                                                    PairStance::TOWER_TOWER);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag->getDoublePrecisionNonbondedKit(devc);
      const PPIKit<double, double4> nrg_tab = pairs_tbl.dpData();
      MMControlKit<double> ctrl = mmctrl->dpData(devc);
      launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl, eval_frc,
                     eval_nrg, has_tiny_box, bt_tp, bt_tt, clash_distance, clash_ratio);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag->getSinglePrecisionNonbondedKit(devc);
      const PPIKit<float, float4> nrg_tab = pairs_tbl.spData();
      MMControlKit<float> ctrl = mmctrl->spData(devc);
      launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl, eval_frc,
                     eval_nrg, has_tiny_box, bt_tp, bt_tt, clash_distance, clash_ratio);
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
                                                    mitigate_clash, PairStance::TOWER_PLATE);
  const int2 bt_tt = launcher.getPMEPairsKernelDims(PrecisionModel::SINGLE, prec,
                                                    NeighborListKind::MONO,
                                                    cg->getTinyBoxPresence(), eval_frc, eval_nrg,
                                                    mitigate_clash, PairStance::TOWER_TOWER);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag->getDoublePrecisionNonbondedKit(devc);
      const PPIKit<double, double4> nrg_tab = pairs_tbl.dpData();
      MMControlKit<double> ctrl = mmctrl->dpData(devc);
      launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw, &tlpn, &scw, &ctrl, eval_frc, eval_nrg,
                     has_tiny_box, bt_tp, bt_tt, clash_distance, clash_ratio);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag->getSinglePrecisionNonbondedKit(devc);
      const PPIKit<float, float4> nrg_tab = pairs_tbl.spData();
      MMControlKit<float> ctrl = mmctrl->spData(devc);
      launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw, &tlpn, &scw, &ctrl, eval_frc, eval_nrg,
                     has_tiny_box, bt_tp, bt_tt, clash_distance, clash_ratio);
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
                                                    eval_frc, eval_nrg, mitigate_clash,
                                                    PairStance::TOWER_PLATE);
  const int2 bt_tt = launcher.getPMEPairsKernelDims(PrecisionModel::SINGLE, prec,
                                                    NeighborListKind::DUAL, has_tiny_box,
                                                    eval_frc, eval_nrg, mitigate_clash,
                                                    PairStance::TOWER_TOWER);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag->getDoublePrecisionNonbondedKit(devc);
      const PPIKit<double, double4> nrg_tab = pairs_tbl.dpData();
      MMControlKit<double> ctrl = mmctrl->dpData(devc);
      launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl, eval_frc,
                     eval_nrg, has_tiny_box, bt_tp, bt_tt, clash_distance, clash_ratio);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag->getSinglePrecisionNonbondedKit(devc);
      const PPIKit<float, float4> nrg_tab = pairs_tbl.spData();
      MMControlKit<float> ctrl = mmctrl->spData(devc);
      launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl, eval_frc,
                     eval_nrg, has_tiny_box, bt_tp, bt_tt, clash_distance, clash_ratio);
    }
    break;
  }
}

} // namespace energy
} // namespace stormm
