// -*-c++-*-
#include "copyright.h"
#include "Constants/fixed_precision.h"
#include "Math/hpc_summation.cuh"
#include "hpc_cellgrid.cuh"

namespace stormm {
namespace energy {

using card::GpuDetails;
using numerics::globalpos_scale_nonoverflow_bits;
using hpc_math::blockExclusivePrefixSum;
using synthesis::PsSynthesisReader;

#include "Accelerator/syncwarp.cui"
#include "Math/rounding.cui"

/// \brief A device function to handle transfer of atoms between one image and the alternate image
///        given a mutable CellGrid abstract.
///
/// \{
__device__ __forceinline__
void atomToNextImage(CellGridWriter<double, llint, double, double4> cgw, const size_t next_pos,
                     const size_t orig_pos) {
  const int topl_idx = __ldg(&cgw.nonimg_atom_idx[orig_pos]);
  __stwt(&cgw.nonimg_atom_idx_alt[next_pos], topl_idx);
  cgw.image_alt[next_pos] = cgw.image[orig_pos];
  __stwt(&cgw.img_atom_idx_alt[topl_idx], next_pos);
  __stwt(&cgw.img_atom_chn_cell_alt[next_pos], __ldcv(&cgw.img_atom_chn_cell[orig_pos]));
  __stwt(&cgw.xfrc[next_pos], 0LL);
  __stwt(&cgw.yfrc[next_pos], 0LL);
  __stwt(&cgw.zfrc[next_pos], 0LL);
  __stwt(&cgw.xfrc_ovrf[next_pos], 0);
  __stwt(&cgw.yfrc_ovrf[next_pos], 0);
  __stwt(&cgw.zfrc_ovrf[next_pos], 0);
}

__device__ __forceinline__
void atomToNextImage(CellGridWriter<float, int, float, float4> cgw, const size_t next_pos,
                     const size_t orig_pos) {
  const int topl_idx = __ldg(&cgw.nonimg_atom_idx[orig_pos]);
  __stwt(&cgw.nonimg_atom_idx_alt[next_pos], topl_idx);
  __stwt(&cgw.image_alt[next_pos], __ldcv(&cgw.image[orig_pos]));
  __stwt(&cgw.img_atom_idx_alt[topl_idx], next_pos);
  __stwt(&cgw.img_atom_chn_cell_alt[next_pos], __ldcv(&cgw.img_atom_chn_cell[orig_pos]));
  __stwt(&cgw.xfrc[next_pos], 0);
  __stwt(&cgw.yfrc[next_pos], 0);
  __stwt(&cgw.zfrc[next_pos], 0);
  __stwt(&cgw.xfrc_ovrf[next_pos], 0);
  __stwt(&cgw.yfrc_ovrf[next_pos], 0);
  __stwt(&cgw.zfrc_ovrf[next_pos], 0);
}
/// \}  

// Enumerate migration kernels involving twin electrostatic and van-der Waals neighbor lists
#define DUAL_GRIDS
#  define TCOORD double
#  define TACC llint
#  define TCOORD4 double4
#  define TCOORD_IS_LONG
#    define FINE_COORDINATES
#      define KERNEL_NAME kdMigrationOneDualFine
#        include "migration_i.cui"
#      undef KERNEL_NAME
#    undef FINE_COORDINATES
#    define KERNEL_NAME kdMigrationOneDual
#      include "migration_i.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kdMigrationTwoDual
#      include "migration_ii.cui"
#    undef KERNEL_NAME
#  undef TCOORD_IS_LONG
#  undef TCOORD4
#  undef TACC
#  undef TCOORD
#  define TCOORD float
#  define TACC int
#  define TCOORD4 float4
#    define FINE_COORDINATES
#      define KERNEL_NAME kfMigrationOneDualFine
#        include "migration_i.cui"
#      undef KERNEL_NAME
#    undef FINE_COORDINATES
#    define KERNEL_NAME kfMigrationOneDual
#      include "migration_i.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfMigrationTwoDual
#      include "migration_ii.cui"
#    undef KERNEL_NAME
#  undef TCOORD4
#  undef TACC
#  undef TCOORD
#undef DUAL_GRIDS

// Enumerate migration kernels involving one unified neighbor list
#define TCOORD double
#define TACC llint
#define TCOORD4 double4
#define TCOORD_IS_LONG
#  define FINE_COORDINATES
#    define KERNEL_NAME kdMigrationOneFine
#      include "migration_i.cui"
#    undef KERNEL_NAME
#  undef FINE_COORDINATES
#  define KERNEL_NAME kdMigrationOne
#    include "migration_i.cui"
#  undef KERNEL_NAME
#  define KERNEL_NAME kdMigrationTwo
#    include "migration_ii.cui"
#  undef KERNEL_NAME
#undef TCOORD_IS_LONG
#undef TCOORD4
#undef TACC
#undef TCOORD
#define TCOORD float
#define TACC int
#define TCOORD4 float4
#  define FINE_COORDINATES
#    define KERNEL_NAME kfMigrationOneFine
#      include "migration_i.cui"
#    undef KERNEL_NAME
#  undef FINE_COORDINATES
#  define KERNEL_NAME kfMigrationOne
#    include "migration_i.cui"
#  undef KERNEL_NAME
#  define KERNEL_NAME kfMigrationTwo
#    include "migration_ii.cui"
#  undef KERNEL_NAME
#undef TCOORD4
#undef TACC
#undef TCOORD

//-------------------------------------------------------------------------------------------------
cudaFuncAttributes queryMigrationKernelRequirements(const PrecisionModel coord_prec,
                                                    const NeighborListKind neighbor_list,
                                                    const int stage, const int gpos_bits) {
  const bool fine_coords = (gpos_bits <= 0 || gpos_bits > globalpos_scale_nonoverflow_bits);
  cudaFuncAttributes result;
  cudaError_t cfa;
  bool stage_error = false;
  switch (coord_prec) {
  case PrecisionModel::DOUBLE:
    switch (neighbor_list) {
    case NeighborListKind::DUAL:
      if (fine_coords) {
        switch (stage) {
        case 1:
          cfa = cudaFuncGetAttributes(&result, kdMigrationOneDualFine);
          break;
        case 2:
        default:
          stage_error = true;
          break;
        }
      }
      else {
        switch (stage) {
        case 1:
          cfa = cudaFuncGetAttributes(&result, kdMigrationOneDual);
          break;
        case 2:
          cfa = cudaFuncGetAttributes(&result, kdMigrationTwoDual);
          break;
        default:
          stage_error = true;
          break;
        }
      }
      break;
    case NeighborListKind::MONO:
      if (fine_coords) {
        switch (stage) {
        case 1:
          cfa = cudaFuncGetAttributes(&result, kdMigrationOneFine);
          break;
        case 2:
        default:
          stage_error = true;
          break;
        }
      }
      else {
        switch (stage) {
        case 1:
          cfa = cudaFuncGetAttributes(&result, kdMigrationOne);
          break;
        case 2:
          cfa = cudaFuncGetAttributes(&result, kdMigrationTwo);
          break;
        default:
          stage_error = true;
          break;
        }
      }
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (neighbor_list) {
    case NeighborListKind::DUAL:
      if (fine_coords) {
        switch (stage) {
        case 1:
          cfa = cudaFuncGetAttributes(&result, kfMigrationOneDualFine);
          break;
        case 2:
        default:
          stage_error = true;
          break;
        }
      }
      else {
        switch (stage) {
        case 1:
          cfa = cudaFuncGetAttributes(&result, kfMigrationOneDual);
          break;
        case 2:
          cfa = cudaFuncGetAttributes(&result, kfMigrationTwoDual);
          break;
        default:
          stage_error = true;
          break;
        }
      }
      break;
    case NeighborListKind::MONO:
      if (fine_coords) {
        switch (stage) {
        case 1:
          cfa = cudaFuncGetAttributes(&result, kfMigrationOneFine);
          break;
        case 2:
        default:
          stage_error = true;
          break;
        }
      }
      else {
        switch (stage) {
        case 1:
          cfa = cudaFuncGetAttributes(&result, kfMigrationOne);
          break;
        case 2:
          cfa = cudaFuncGetAttributes(&result, kfMigrationTwo);
          break;
        default:
          stage_error = true;
          break;
        }
      }
      break;
    }
    break;
  }

  // Check for errors
  if (stage_error) {
    rtErr("An invalid migration stage " + std::to_string(stage) + " was set for a " +
          getEnumerationName(coord_prec) + " coordinate representation and a " +
          getEnumerationName(neighbor_list) + " neighbor list configuration.  The "
          "fixed-precision coordinate representation entails " + std::to_string(gpos_bits) +
          " bits after the point.", "queryMigrationKernelRequirements");
  }
  if (cfa != cudaSuccess) {

    // Construct the appropriate error message
    std::string error_message("Error obtaining attributes for kernel k");
    switch (coord_prec) {
    case PrecisionModel::DOUBLE:
      error_message += "d";
      break;
    case PrecisionModel::SINGLE:
      error_message += "f";
      break;
    }
    error_message += "Migration";
    switch (stage) {
    case 1:
      error_message += "One";
      break;
    case 2:
      error_message += "Two";
      break;
    default:
      error_message += "_" + std::to_string(stage) + "_";
    }
    switch (neighbor_list) {
    case NeighborListKind::DUAL:
      error_message += "Dual";
      break;
    case NeighborListKind::MONO:
      break;
    }
    if (fine_coords) {
      error_message += "Fine";
    }

    // Report the error
    rtErr(error_message, "queryMigrationKernelRequirements");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void launchMigration(CellGridWriter<double, llint, double, double4> *cgw,
                     const CellOriginsReader &corg, const PsSynthesisReader &poly_psr,
                     const int2 bt_i, const int2 bt_ii) {
  if (poly_psr.gpos_bits > globalpos_scale_nonoverflow_bits) {
    kdMigrationOneFine<<<bt_i.x, bt_i.y>>>(poly_psr, *cgw, corg);
  }
  else {
    kdMigrationOne<<<bt_i.x, bt_i.y>>>(poly_psr, *cgw, corg);
  }
  kdMigrationTwo<<<bt_ii.x, bt_ii.y>>>(*cgw);
}

//-------------------------------------------------------------------------------------------------
void launchMigration(CellGridWriter<double, llint, double, double4> *cgw_qq,
                     CellGridWriter<double, llint, double, double4> *cgw_lj,
                     const CellOriginsReader &corg_qq, const CellOriginsReader &corg_lj,
                     const PsSynthesisReader &poly_psr, const int2 bt_i, const int2 bt_ii) {
  if (poly_psr.gpos_bits > globalpos_scale_nonoverflow_bits) {
    kdMigrationOneDualFine<<<bt_i.x, bt_i.y>>>(poly_psr, *cgw_qq, *cgw_lj, corg_qq, corg_lj);
  }
  else {
    kdMigrationOneDual<<<bt_i.x, bt_i.y>>>(poly_psr, *cgw_qq, *cgw_lj, corg_qq, corg_lj);
  }
  kdMigrationTwoDual<<<bt_ii.x, bt_ii.y>>>(*cgw_qq, *cgw_lj);
}

//-------------------------------------------------------------------------------------------------
void launchMigration(CellGridWriter<float, int, float, float4> *cgw,
                     const CellOriginsReader &corg, const PsSynthesisReader &poly_psr,
                     const int2 bt_i, const int2 bt_ii) {
  if (poly_psr.gpos_bits > globalpos_scale_nonoverflow_bits) {
    kfMigrationOneFine<<<bt_i.x, bt_i.y>>>(poly_psr, *cgw, corg);
  }
  else {
    kfMigrationOne<<<bt_i.x, bt_i.y>>>(poly_psr, *cgw, corg);
  }
  kfMigrationTwo<<<bt_ii.x, bt_ii.y>>>(*cgw);
}

//-------------------------------------------------------------------------------------------------
void launchMigration(CellGridWriter<float, int, float, float4> *cgw_qq,
                     CellGridWriter<float, int, float, float4> *cgw_lj,
                     const CellOriginsReader &corg_qq, const CellOriginsReader &corg_lj,
                     const PsSynthesisReader &poly_psr, const int2 bt_i, const int2 bt_ii) {
  if (poly_psr.gpos_bits > globalpos_scale_nonoverflow_bits) {
    kfMigrationOneDualFine<<<bt_i.x, bt_i.y>>>(poly_psr, *cgw_qq, *cgw_lj, corg_qq, corg_lj);
  }
  else {
    kfMigrationOneDual<<<bt_i.x, bt_i.y>>>(poly_psr, *cgw_qq, *cgw_lj, corg_qq, corg_lj);
  }
  kfMigrationTwoDual<<<bt_ii.x, bt_ii.y>>>(*cgw_qq, *cgw_lj);
}

//-------------------------------------------------------------------------------------------------
void launchMigration(CellGrid<double, llint, double, double4> *cg,
                     const PhaseSpaceSynthesis &poly_ps, const CoreKlManager &launcher) {
  CellGridWriter<double, llint, double, double4> cgw = cg->data(HybridTargetLevel::DEVICE);
  const CoordinateCycle next_ori = getNextCyclePosition(cg->getCyclePosition());
  const CellOriginsReader corg = cg->getRulers(next_ori, HybridTargetLevel::DEVICE);
  const PsSynthesisReader poly_psr = poly_ps.data(HybridTargetLevel::DEVICE);
  const int2 bt_i = launcher.getMigrationKernelDims(PrecisionModel::DOUBLE, NeighborListKind::MONO,
                                                    1, poly_psr.gpos_bits, cgw.total_chain_count);
  const int2 bt_ii = launcher.getMigrationKernelDims(PrecisionModel::DOUBLE,
                                                     NeighborListKind::MONO, 2, poly_psr.gpos_bits,
                                                     cgw.total_chain_count);
  launchMigration(&cgw, corg, poly_psr, bt_i, bt_ii);
}

//-------------------------------------------------------------------------------------------------
void launchMigration(CellGrid<double, llint, double, double4> *cg_qq,
                     CellGrid<double, llint, double, double4> *cg_lj,
                     const PhaseSpaceSynthesis &poly_ps, const CoreKlManager &launcher) {
  CellGridWriter<double, llint, double, double4> cgw_qq = cg_qq->data(HybridTargetLevel::DEVICE);
  CellGridWriter<double, llint, double, double4> cgw_lj = cg_lj->data(HybridTargetLevel::DEVICE);
  const CoordinateCycle next_ori = getNextCyclePosition(cg_qq->getCyclePosition());
  const CellOriginsReader corg_qq = cg_qq->getRulers(next_ori, HybridTargetLevel::DEVICE);
  const CellOriginsReader corg_lj = cg_lj->getRulers(next_ori, HybridTargetLevel::DEVICE);
  const PsSynthesisReader poly_psr = poly_ps.data(HybridTargetLevel::DEVICE);
  const int chain_count = cgw_qq.total_chain_count + cgw_lj.total_chain_count;
  const int2 bt_i = launcher.getMigrationKernelDims(PrecisionModel::DOUBLE, NeighborListKind::DUAL,
                                                    1, poly_psr.gpos_bits, chain_count);
  const int2 bt_ii = launcher.getMigrationKernelDims(PrecisionModel::DOUBLE,
                                                     NeighborListKind::DUAL, 2, poly_psr.gpos_bits,
                                                     chain_count);
  launchMigration(&cgw_qq, &cgw_lj, corg_qq, corg_lj, poly_psr, bt_i, bt_ii);
}

//-------------------------------------------------------------------------------------------------
void launchMigration(CellGrid<float, int, float, float4> *cg,
                     const PhaseSpaceSynthesis &poly_ps, const CoreKlManager &launcher) {
  CellGridWriter<float, int, float, float4> cgw = cg->data(HybridTargetLevel::DEVICE);
  const CoordinateCycle next_ori = getNextCyclePosition(cg->getCyclePosition());
  const CellOriginsReader corg = cg->getRulers(next_ori, HybridTargetLevel::DEVICE);
  const PsSynthesisReader poly_psr = poly_ps.data(HybridTargetLevel::DEVICE);
  const int2 bt_i = launcher.getMigrationKernelDims(PrecisionModel::SINGLE, NeighborListKind::MONO,
                                                    1, poly_psr.gpos_bits, cgw.total_chain_count);
  const int2 bt_ii = launcher.getMigrationKernelDims(PrecisionModel::SINGLE,
                                                     NeighborListKind::MONO, 2, poly_psr.gpos_bits,
                                                     cgw.total_chain_count);
  launchMigration(&cgw, corg, poly_psr, bt_i, bt_ii);
}

//-------------------------------------------------------------------------------------------------
void launchMigration(CellGrid<float, int, float, float4> *cg_qq,
                     CellGrid<float, int, float, float4> *cg_lj,
                     const PhaseSpaceSynthesis &poly_ps, const CoreKlManager &launcher) {
  CellGridWriter<float, int, float, float4> cgw_qq = cg_qq->data(HybridTargetLevel::DEVICE);
  CellGridWriter<float, int, float, float4> cgw_lj = cg_lj->data(HybridTargetLevel::DEVICE);
  const CoordinateCycle next_ori = getNextCyclePosition(cg_qq->getCyclePosition());
  const CellOriginsReader corg_qq = cg_qq->getRulers(next_ori, HybridTargetLevel::DEVICE);
  const CellOriginsReader corg_lj = cg_lj->getRulers(next_ori, HybridTargetLevel::DEVICE);
  const PsSynthesisReader poly_psr = poly_ps.data(HybridTargetLevel::DEVICE);
  const int chain_count = cgw_qq.total_chain_count + cgw_lj.total_chain_count;
  const int2 bt_i = launcher.getMigrationKernelDims(PrecisionModel::SINGLE, NeighborListKind::DUAL,
                                                    1, poly_psr.gpos_bits, chain_count);
  const int2 bt_ii = launcher.getMigrationKernelDims(PrecisionModel::SINGLE,
                                                     NeighborListKind::DUAL, 2, poly_psr.gpos_bits,
                                                     chain_count);
  launchMigration(&cgw_qq, &cgw_lj, corg_qq, corg_lj, poly_psr, bt_i, bt_ii);
}

//-------------------------------------------------------------------------------------------------
void launchCellGridAction(CellGridWriter<void, void, void, void> *cgw, const size_t tc_mat,
                          const size_t tc_acc, const GpuDetails &gpu,
                          const CellGridAction process) {

  // Restoring the type of the CellGrid object carries many possibilities, although only the type
  // of the accumulators is essential to their initialization.
  if (tc_mat == double_type_index) {
    unrollLaunchCellGridAction<double, double, double4>(cgw, tc_acc, gpu, process);
  }
  else if (tc_mat == float_type_index) {
    unrollLaunchCellGridAction<float, float, float4>(cgw, tc_acc, gpu, process);
  }
  else if (tc_mat == int_type_index) {
    unrollLaunchCellGridAction<int, float, int4>(cgw, tc_acc, gpu, process);
  }
  else if (tc_mat == llint_type_index) {
    unrollLaunchCellGridAction<llint, double, llint4>(cgw, tc_acc, gpu, process);
  }
  else {
    rtErr("The CellGrid object must take either 32-bit or 64-bit signed integers as the primary "
          "accumulator.", "launchCellGridAction");
  }
}

//-------------------------------------------------------------------------------------------------
void launchCellGridAction(CellGridWriter<void, void, void, void> *cgw, const size_t tc_mat,
                          const size_t tc_acc, const PsSynthesisReader &poly_psr,
                          const GpuDetails &gpu, const CellGridAction process) {

  // Restoring the type of the CellGrid object carries many possibilities, although only the type
  // of the accumulators is essential to their initialization.
  if (tc_mat == double_type_index) {
    unrollLaunchCellGridAction<double, double, double4>(cgw, tc_acc, poly_psr, gpu, process);
  }
  else if (tc_mat == float_type_index) {
    unrollLaunchCellGridAction<float, float, float4>(cgw, tc_acc, poly_psr, gpu, process);
  }
  else if (tc_mat == int_type_index) {
    unrollLaunchCellGridAction<int, float, int4>(cgw, tc_acc, poly_psr, gpu, process);
  }
  else if (tc_mat == llint_type_index) {
    unrollLaunchCellGridAction<llint, double, llint4>(cgw, tc_acc, poly_psr, gpu, process);
  }
  else {
    rtErr("The CellGrid object must take either 32-bit or 64-bit signed integers as the primary "
          "accumulator.", "launchCellGridAction");
  }
}

//-------------------------------------------------------------------------------------------------
void launchCellGridAction(const CellGridReader<void, void, void, void> &cgr, const size_t tc_mat,
                          const size_t tc_acc, PsSynthesisWriter *poly_psw, const GpuDetails &gpu,
                          const CellGridAction process) {

  // Restoring the type of the CellGrid object carries many possibilities, although only the type
  // of the accumulators is essential to their initialization.
  if (tc_mat == double_type_index) {
    unrollLaunchCellGridAction<double, double, double4>(cgr, tc_acc, poly_psw, gpu, process);
  }
  else if (tc_mat == float_type_index) {
    unrollLaunchCellGridAction<float, float, float4>(cgr, tc_acc, poly_psw, gpu, process);
  }
  else if (tc_mat == int_type_index) {
    unrollLaunchCellGridAction<int, float, int4>(cgr, tc_acc, poly_psw, gpu, process);
  }
  else if (tc_mat == llint_type_index) {
    unrollLaunchCellGridAction<llint, double, llint4>(cgr, tc_acc, poly_psw, gpu, process);
  }
  else {
    rtErr("The CellGrid object must take either 32-bit or 64-bit signed integers as the primary "
          "accumulator.", "launchCellGridAction");
  }
}

} // namespace energy
} // namespace stormm
