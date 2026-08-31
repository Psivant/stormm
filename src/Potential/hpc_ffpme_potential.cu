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
#include "evaluate_localmask.cui"

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

#  undef LLCONV_FUNC
#  undef SQRT_FUNC
  
#  undef TCOORD
#  undef TCOORD4
#  undef TACC

#  undef TCOORD_IS_REAL
#  undef TCALC_IS_SINGLE
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
extern cudaFuncAttributes queryFFPMEPairsKernelRequirements(const NeighborListKind neighbor_layout,
                                                            const EvaluateForce eval_frc,
                                                            const EvaluateEnergy eval_nrg,
                                                            const TinyBoxPresence has_tiny_box,
                                                            const ClashResponse clash_handling) {

  // As with other kernel querying functions, the kernel manager calling this function will have
  // specifications of the GPU in use.  It is the overall thread occupancy and multiplicity of each
  // kernel that this function must return.
  cudaFuncAttributes result;
  cudaError_t cfa = cudaErrorInvalidValue;
  switch (clash_handling) {
  case ClashResponse::NONE:
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (neighbor_layout) {
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
        switch (neighbor_layout) {
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
      switch (neighbor_layout) {
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
  case ClashResponse::FORGIVE:
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (neighbor_layout) {
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
        switch (neighbor_layout) {
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
      switch (neighbor_layout) {
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

  // Check for errors
  if (cfa != cudaSuccess) {

    // Construct the appropriate error message
    std::string error_message("Error obtaining attributes from kernel kff");
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
    switch (neighbor_layout) {
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
    rtErr(error_message, "queryFFPMEPairsKernelRequirements");
  }
  return result;
}
#endif

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

} // namespace energy
} // namespace stormm
