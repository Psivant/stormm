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

// Double-precision tile evaluation
#define TCALC double
#  define TCALC2 double2
#  define TCALC3 double3
#  define TCALC4 double4_16a
#  define TCOORD_IS_REAL

// Other definitions associated with 64-bit floating-point arithmetic
#  define SQRT_FUNC sqrt
#  define LLCONV_FUNC __double2ll_rn

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
#        define PMENB_WARPS_PER_BLOCK 10
#        define KERNEL_NAME kfdTowerPlateFEDual
#          include "tower_plate_pairs.cui"
#        undef KERNEL_NAME
#        undef PMENB_WARPS_PER_BLOCK
#      undef DUAL_GRIDS
#      define TINY_BOX
#        define CLASH_FORGIVENESS
#          define PMENB_WARPS_PER_BLOCK 10
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

#  undef TCOORD
#  undef TCOORD4
#  undef TACC

#  undef LLCONV_FUNC
#  undef SQRT_FUNC
  
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
extern cudaFuncAttributes queryFDPMEPairsKernelRequirements(const NeighborListKind neighbor_layout,
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
        switch (neighbor_layout) {
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
      switch (neighbor_layout) {
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
  case ClashResponse::FORGIVE:
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (neighbor_layout) {
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
        switch (neighbor_layout) {
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
      switch (neighbor_layout) {
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
  }

  // Check for errors
  if (cfa != cudaSuccess) {

    // Construct the appropriate error message
    std::string error_message("Error obtaining attributes from kernel kfd");
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
    rtErr(error_message, "queryFDPMEPairsKernelRequirements");
  }
  return result;
}
#endif

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<double, double4_16a> &nrg_tab,
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
                           const PPIKit<double, double4_16a> &nrg_tab,
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
                           const PPIKit<double, double4_16a> &nrg_tab,
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
                           const PPIKit<double, double4_16a> &nrg_tab,
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

} // namespace energy
} // namespace stormm
