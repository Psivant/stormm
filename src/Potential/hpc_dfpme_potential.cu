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

// Other definitions associated with 64-bit floating-point arithmetic
#  define SQRT_FUNC sqrtf
#  define LLCONV_FUNC __float2ll_rn

// Compile additional kernels for the float64_t coordinate representation.
#  define TCOORD double
#  define TCOORD4 double4_16a
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
  
#  undef TCOORD
#  undef TCOORD4
#  undef TACC
#  undef TCOORD_IS_LONG

#  undef LLCONV_FUNC
#  undef SQRT_FUNC
  
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
extern cudaFuncAttributes queryDFPMEPairsKernelRequirements(const NeighborListKind neighbor_layout,
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
        switch (neighbor_layout) {
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
      switch (neighbor_layout) {
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
  case ClashResponse::FORGIVE:
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (neighbor_layout) {
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
        switch (neighbor_layout) {
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
      switch (neighbor_layout) {
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

  // Check for errors
  if (cfa != cudaSuccess) {

    // Construct the appropriate error message
    std::string error_message("Error obtaining attributes from kernel kdf");
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
    rtErr(error_message, "queryDFPMEPairsKernelRequirements");
  }
  return result;
}
#endif

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<float, float4> &nrg_tab,
                           CellGridWriter<double, llint, double, double4_16a> *cgw, TilePlan *tlpn,
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
                           CellGridWriter<double, llint, double, double4_16a> *cgw, TilePlan *tlpn,
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
                           CellGridWriter<double, llint, double, double4_16a> *cgw_qq,
                           CellGridWriter<double, llint, double, double4_16a> *cgw_lj,
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
                           CellGridWriter<double, llint, double, double4_16a> *cgw_qq,
                           CellGridWriter<double, llint, double, double4_16a> *cgw_lj,
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

} // namespace energy
} // namespace stormm
