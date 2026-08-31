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

// Begin with a float64_t coordinate representation, appropriate for the float64_t arithmetic mode.
#  define TCOORD double
#  define TCOORD4 double4_16a
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

#  undef LLCONV_FUNC
#  undef SQRT_FUNC

#  undef TCOORD
#  undef TCOORD4
#  undef TACC
#  undef TCOORD_IS_LONG
  
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
extern cudaFuncAttributes queryDDPMEPairsKernelRequirements(const NeighborListKind neighbor_layout,
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
        switch (neighbor_layout) {
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
      switch (neighbor_layout) {
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
  case ClashResponse::FORGIVE:
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (neighbor_layout) {
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
        switch (neighbor_layout) {
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
      switch (neighbor_layout) {
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
  }

  // Check for errors
  if (cfa != cudaSuccess) {

    // Construct the appropriate error message
    std::string error_message("Error obtaining attributes from kernel kdd");
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
    rtErr(error_message, "queryDDPMEPairsKernelRequirements");
  }

  return result;
}
#endif

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                           const LocalExclusionMaskReader &lemr,
                           const PPIKit<double, double4_16a> &nrg_tab,
                           CellGridWriter<double, llint, double, double4_16a> *cgw, TilePlan *tlpn,
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
                           const PPIKit<double, double4_16a> &nrg_tab,
                           const PsSynthesisBorders &sysbrd,
                           CellGridWriter<double, llint, double, double4_16a> *cgw, TilePlan *tlpn,
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
                           const PPIKit<double, double4_16a> &nrg_tab,
                           CellGridWriter<double, llint, double, double4_16a> *cgw_qq,
                           CellGridWriter<double, llint, double, double4_16a> *cgw_lj,
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
                           const PPIKit<double, double4_16a> &nrg_tab,
                           const PsSynthesisBorders &sysbrd,
                           CellGridWriter<double, llint, double, double4_16a> *cgw_qq,
                           CellGridWriter<double, llint, double, double4_16a> *cgw_lj,
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

} // namespace energy
} // namespace stormm
