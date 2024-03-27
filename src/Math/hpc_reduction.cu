// -*-c++-*-
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "DataTypes/stormm_vector_types.h"
#include "Numerics/split_fixed_precision.h"
#include "hpc_reduction.h"

namespace stormm {
namespace stmath {

using constants::twice_warp_bits_mask_int;

#include "Numerics/accumulation.cui"

//-------------------------------------------------------------------------------------------------
// Perform an accumulation over the conjugate gradient work units' relevant atoms to obtain the
// sum of squared gradients (gg) and the evolution of the gradient (dgg).
//
// Overloaded:
//   - Operate on standard fixed-precision data in long long int format
//   - Operate on extended fixed-precision data in { long long int, int } format
//
// Arguments:
//   gg_collector:     __shared__ memory resource for storing gg accumulation per warp
//   dgg_collector:    __shared__ memory resource for storing dgg accumulation per warp
//   start_pos:        First atom affected by this work unit in the following global arrays
//   end_pos:          Upper limit of atoms affected by this work unit
//   xfrc:             Current forces on particles in the Cartesian X direction
//   yfrc:             Current forces on particles in the Cartesian Y direction
//   zfrc:             Current forces on particles in the Cartesian Z direction
//   [x,y,z]frc_ovrf:  Overflow in X, Y, or Z forces
//   xprv:             Prior iteration forces on particles in the Cartesian X direction
//   yprv:             Prior iteration forces on particles in the Cartesian Y direction
//   zprv:             Prior iteration forces on particles in the Cartesian Z direction
//   [x,y,z]prv_ovrf:  Overflow in prior iteration X, Y, or Z forces
//-------------------------------------------------------------------------------------------------
__device__ __forceinline__
void conjGradCoreGather(double* gg_collector, double* dgg_collector, const int start_pos,
                        const int end_pos, const double inv_frc_scale, const llint* xfrc,
                        const llint* yfrc, const llint* zfrc, const llint* xprv, const llint* yprv,
                        const llint* zprv) {
  double gg  = 0.0;
  double dgg = 0.0;
  const int natom = end_pos - start_pos;
  const int padded_natom = (((natom + warp_size_int - 1) >> warp_bits) << warp_bits);
  int pos = threadIdx.x;
  while (pos < padded_natom) {
    if (pos < natom) {
      const llint llpx = xprv[start_pos + pos];
      const llint llfx = xfrc[start_pos + pos];
      const llint lldx = llfx - llpx;
      const double dpx = (double)(llpx);
      const double dfx = (double)(llfx);
      const double ddx = (double)(lldx);
      gg += (dpx * dpx);
      dgg += (ddx * dfx);
    }
    pos += blockDim.x;
  }
  while (pos < 2 * padded_natom) {
    const int rel_pos = pos - padded_natom;
    if (rel_pos < natom) {
      const llint llpy = yprv[start_pos + rel_pos];
      const llint llfy = yfrc[start_pos + rel_pos];
      const llint lldy = llfy - llpy;
      const double dpy = (double)(llpy);
      const double dfy = (double)(llfy);
      const double ddy = (double)(lldy);
      gg += (dpy * dpy);
      dgg += (ddy * dfy);
    }
    pos += blockDim.x;
  }
  while (pos < 3 * padded_natom) {
    const int rel_pos = pos - (2 * padded_natom);
    if (rel_pos < natom) {
      const llint llpz = zprv[start_pos + rel_pos];
      const llint llfz = zfrc[start_pos + rel_pos];
      const llint lldz = llfz - llpz;
      const double dpz = (double)(llpz);
      const double dfz = (double)(llfz);
      const double ddz = (double)(lldz);
      gg += (dpz * dpz);
      dgg += (ddz * dfz);
    }
    pos += blockDim.x;
  }
  WARP_REDUCE_DOWN(gg);
  WARP_REDUCE_DOWN(dgg);
  if ((threadIdx.x & warp_bits_mask_int) == 0) {
    const size_t warp_idx = (threadIdx.x >> warp_bits);
    gg_collector[warp_idx]  =  gg;
    dgg_collector[warp_idx] = dgg;
  }  
}

__device__ __forceinline__
void conjGradCoreGather(double* gg_collector, double* dgg_collector, const int start_pos,
                        const int end_pos, const double inv_frc_scale, const llint* xfrc,
                        const int* xfrc_ovrf, const llint* yfrc, const int* yfrc_ovrf,
                        const llint* zfrc, const int* zfrc_ovrf, const llint* xprv,
                        const int* xprv_ovrf, const llint* yprv, const int* yprv_ovrf,
                        const llint* zprv, const int* zprv_ovrf) {
  double gg  = 0.0;
  double dgg = 0.0;
  for (int tpos = start_pos + threadIdx.x; tpos < end_pos; tpos += blockDim.x) {
    const double dpx = (((double)(xprv_ovrf[tpos]) * max_llint_accumulation) +
                        (double)(xprv[tpos])) * inv_frc_scale;
    const double dpy = (((double)(yprv_ovrf[tpos]) * max_llint_accumulation) +
                        (double)(yprv[tpos])) * inv_frc_scale;
    const double dpz = (((double)(zprv_ovrf[tpos]) * max_llint_accumulation) +
                        (double)(zprv[tpos])) * inv_frc_scale;
    gg += (dpx * dpx) + (dpy * dpy) + (dpz * dpz);
    const double dfx = (((double)(xfrc_ovrf[tpos]) * max_llint_accumulation) +
                        (double)(xfrc[tpos])) * inv_frc_scale;
    const double dfy = (((double)(yfrc_ovrf[tpos]) * max_llint_accumulation) +
                        (double)(yfrc[tpos])) * inv_frc_scale;
    const double dfz = (((double)(zfrc_ovrf[tpos]) * max_llint_accumulation) +
                        (double)(zfrc[tpos])) * inv_frc_scale;
    const double ddx = dfx - dpx;
    const double ddy = dfy - dpy;
    const double ddz = dfz - dpz;
    dgg += (ddx * dfx) + (ddy * dfy) + (ddz * dfz);
  }
  WARP_REDUCE_DOWN(gg);
  WARP_REDUCE_DOWN(dgg);
  if ((threadIdx.x & warp_bits_mask_int) == 0) {
    const size_t warp_idx = (threadIdx.x >> warp_bits);
    gg_collector[warp_idx]  =  gg;
    dgg_collector[warp_idx] = dgg;
  }
}

//-------------------------------------------------------------------------------------------------
// Distribute the results of the conjugate gradient calculation to normalize the forces and update
// prior holdings for the next cycle.
//
// Arguments:
//   gam:              Contribution factor applied to prior conjugate gradient results
//   start_pos:        First atom affected by this work unit in the following global arrays
//   end_pos:          Upper limit of atoms affected by this work unit
//   xfrc:             Current forces on particles in the Cartesian X direction
//   yfrc:             Current forces on particles in the Cartesian Y direction
//   zfrc:             Current forces on particles in the Cartesian Z direction
//   [x,y,z]frc_ovrf:  Overflow in X, Y, or Z forces
//   xprv:             Prior iteration forces on particles in the Cartesian X direction
//   yprv:             Prior iteration forces on particles in the Cartesian Y direction
//   zprv:             Prior iteration forces on particles in the Cartesian Z direction
//   [x,y,z]prv_ovrf:  Overflow in prior iteration X, Y, or Z forces
//-------------------------------------------------------------------------------------------------
__device__ __forceinline__
double conjGradScatter(const double gam, double* msum_collector, const int atom_start_pos,
                       const int atom_end_pos, llint* xfrc, llint* yfrc, llint* zfrc, llint* xprv,
                       llint* yprv, llint* zprv, llint* x_cg_temp, llint* y_cg_temp,
                       llint* z_cg_temp) {
  double msum = 0.0;
  const int natom = atom_end_pos - atom_start_pos;
  const int padded_natom = (((natom + warp_size_int - 1) >> warp_bits) << warp_bits);
  int pos = threadIdx.x;
  while (pos < padded_natom) {
    if (pos < natom) {
      const size_t gbl_pos = atom_start_pos + pos;
      const llint ifx = xfrc[gbl_pos];
      xprv[gbl_pos] = ifx;
      if (fabs(gam) > constants::verytiny) {
        const double dfcgx = (double)(ifx) + (gam * (double)(x_cg_temp[gbl_pos]));
        msum += (dfcgx * dfcgx);
        const llint cg_x = __double2ll_rn(dfcgx);
        x_cg_temp[gbl_pos] = cg_x;
        xfrc[gbl_pos] = cg_x;
      }
      else {
        const double dfcgx = (double)(ifx);
        msum += (dfcgx * dfcgx);
        x_cg_temp[gbl_pos] = ifx;
      }
    }
    pos += blockDim.x;
  }
  while (pos < 2 * padded_natom) {
    const int rel_pos = pos - padded_natom;
    if (rel_pos < natom) {
      const size_t gbl_pos = atom_start_pos + rel_pos;
      const llint ify = yfrc[gbl_pos];
      yprv[gbl_pos] = ify;
      if (fabs(gam) > constants::verytiny) {
        const double dfcgy = (double)(ify) + (gam * (double)(y_cg_temp[gbl_pos]));
        msum += (dfcgy * dfcgy);
        const llint cg_y = __double2ll_rn(dfcgy);
        y_cg_temp[gbl_pos] = cg_y;
        yfrc[gbl_pos] = cg_y;
      }
      else {
        const double dfcgy = (double)(ify);
        msum += (dfcgy * dfcgy);
        y_cg_temp[gbl_pos] = ify;
      }
    }
    pos += blockDim.x;
  }
  while (pos < 3 * padded_natom) {
    const int rel_pos = pos - (2 * padded_natom);
    if (rel_pos < natom) {
      const size_t gbl_pos = atom_start_pos + rel_pos;
      const llint ifz = zfrc[gbl_pos];
      zprv[gbl_pos] = ifz;
      if (fabs(gam) > constants::verytiny) {
        const double dfcgz = (double)(ifz) + (gam * (double)(z_cg_temp[gbl_pos]));
        msum += (dfcgz * dfcgz);
        const llint cg_z = __double2ll_rn(dfcgz);
        z_cg_temp[gbl_pos] = cg_z;
        zfrc[gbl_pos] = cg_z;
      }
      else {
        const double dfcgz = (double)(ifz);
        msum += (dfcgz * dfcgz);
        z_cg_temp[gbl_pos] = ifz;
      }
    }
    pos += blockDim.x;
  }
  WARP_REDUCE_DOWN(msum);
  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (lane_idx == 0) {
    msum_collector[warp_idx] = msum;
  }
  __syncthreads();

  // All warps will copy the final summation over msum so that the correct accumulated value may
  // be returned on all threads.
  msum = (lane_idx < (blockDim.x >> warp_bits)) ? msum_collector[lane_idx] : 0.0;
  if (blockDim.x == 4 * warp_size_int) {
    msum += SHFL_DOWN(msum, 2);
    msum += SHFL_DOWN(msum, 1);
  }
  else {
    WARP_REDUCE_DOWN(msum);
  }
  msum = SHFL(msum, 0);
  return msum;
}

__device__ __forceinline__
double conjGradScatter(const double gam, double* msum_collector, const int atom_start_pos,
                       const int atom_end_pos, llint* xfrc, int* xfrc_ovrf, llint* yfrc,
                       int* yfrc_ovrf, llint* zfrc, int* zfrc_ovrf, llint* xprv, int* xprv_ovrf,
                       llint* yprv, int* yprv_ovrf, llint* zprv, int* zprv_ovrf, llint* x_cg_temp,
                       int* x_cg_temp_ovrf, llint* y_cg_temp, int* y_cg_temp_ovrf,
                       llint* z_cg_temp, int* z_cg_temp_ovrf) {
  double msum = 0.0;
  for (int tpos = atom_start_pos + threadIdx.x; tpos < atom_end_pos; tpos += blockDim.x) {
    const llint ifx = xfrc[tpos];
    const llint ify = yfrc[tpos];
    const llint ifz = zfrc[tpos];
    xprv[tpos] = ifx;
    yprv[tpos] = ify;
    zprv[tpos] = ifz;
    const int ifx_ovrf = xfrc_ovrf[tpos];
    const int ify_ovrf = yfrc_ovrf[tpos];
    const int ifz_ovrf = zfrc_ovrf[tpos];
    xprv_ovrf[tpos] = ifx_ovrf;
    yprv_ovrf[tpos] = ify_ovrf;
    zprv_ovrf[tpos] = ifz_ovrf;
    const double  fx_part = ((double)(ifx_ovrf) * max_llint_accumulation) + (double)(ifx);
    const double cgx_part = ((double)(x_cg_temp_ovrf[tpos]) * max_llint_accumulation) +
                            (double)(x_cg_temp[tpos]);
    const double  fy_part = ((double)(ify_ovrf) * max_llint_accumulation) + (double)(ify);
    const double cgy_part = ((double)(y_cg_temp_ovrf[tpos]) * max_llint_accumulation) +
                            (double)(y_cg_temp[tpos]);
    const double  fz_part = ((double)(ifz_ovrf) * max_llint_accumulation) + (double)(ifz);
    const double cgz_part = ((double)(z_cg_temp_ovrf[tpos]) * max_llint_accumulation) +
                            (double)(z_cg_temp[tpos]);
    if (fabs(gam) > constants::verytiny) {
      const double dfcgx = fx_part + (gam * cgx_part);
      const double dfcgy = fy_part + (gam * cgy_part);
      const double dfcgz = fz_part + (gam * cgz_part);
      msum += (dfcgx * dfcgx) + (dfcgy * dfcgy) + (dfcgz * dfcgz);
      const int95_t cg_x = doubleToInt95(dfcgx);
      const int95_t cg_y = doubleToInt95(dfcgy);
      const int95_t cg_z = doubleToInt95(dfcgz);
      x_cg_temp[tpos] = cg_x.x;
      y_cg_temp[tpos] = cg_y.x;
      z_cg_temp[tpos] = cg_z.x;
      x_cg_temp_ovrf[tpos] = cg_x.y;
      y_cg_temp_ovrf[tpos] = cg_y.y;
      z_cg_temp_ovrf[tpos] = cg_z.y;
      xfrc[tpos] = cg_x.x;
      yfrc[tpos] = cg_y.x;
      zfrc[tpos] = cg_z.x;
      xfrc_ovrf[tpos] = cg_x.y;
      yfrc_ovrf[tpos] = cg_y.y;
      zfrc_ovrf[tpos] = cg_z.y;
    }
    else {
      msum += ((fx_part * fx_part) + (fy_part * fy_part) + (fz_part * fz_part));
      x_cg_temp[tpos] = ifx;
      y_cg_temp[tpos] = ify;
      z_cg_temp[tpos] = ifz;
      x_cg_temp_ovrf[tpos] = ifx_ovrf;
      y_cg_temp_ovrf[tpos] = ify_ovrf;
      z_cg_temp_ovrf[tpos] = ifz_ovrf;
    }
  }
  WARP_REDUCE_DOWN(msum);
  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (lane_idx == 0) {
    msum_collector[warp_idx] = msum;
  }
  __syncthreads();

  // All warps will copy the final summation over msum so that the correct accumulated value may
  // be returned on all threads.
  msum = (lane_idx < (blockDim.x >> warp_bits)) ? msum_collector[lane_idx] : 0.0;
  WARP_REDUCE_DOWN(msum);
  msum = SHFL(msum, 0);
  return msum;
}

// Double-precision floating point conjugate gradient definitions
#define TCALC            double
#define TCALC_IS_DOUBLE
#define KGATHER_NAME     kdgtConjGrad
#define KSCATTER_NAME    kdscConjGrad
#define KALLREDUCE_NAME  kdrdConjGrad
#define KRESCALE_NAME    kdrsConjGrad
#include "conjugate_gradient.cui"
#undef KGATHER_NAME
#undef KSCATTER_NAME
#undef KALLREDUCE_NAME
#undef KRESCALE_NAME
#undef TCALC
#undef TCALC_IS_DOUBLE
  
// Single-precision floating point conjugate gradient definitions.  The single-precision form
// still does its accumulation in double-precision, but does not store its data in the extended
// fixed-precision format which the double-precision forms of the kernels read and write.
#define TCALC            float
#define KGATHER_NAME     kfgtConjGrad
#define KSCATTER_NAME    kfscConjGrad
#define KALLREDUCE_NAME  kfrdConjGrad
#define KRESCALE_NAME    kfrsConjGrad
#include "conjugate_gradient.cui"
#undef KGATHER_NAME
#undef KSCATTER_NAME
#undef KALLREDUCE_NAME
#undef KRESCALE_NAME
#undef TCALC

//-------------------------------------------------------------------------------------------------
extern void reductionKernelSetup() {
  const cudaSharedMemConfig sms_eight = cudaSharedMemBankSizeEightByte;
  if (cudaFuncSetSharedMemConfig(kdgtConjGrad, sms_eight) != cudaSuccess) {
    rtErr("Error setting kdgtConjGrad __shared__ memory bank size to eight bytes.",
          "reductionKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kdscConjGrad, sms_eight) != cudaSuccess) {
    rtErr("Error setting kdscConjGrad __shared__ memory bank size to eight bytes.",
          "reductionKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kdrsConjGrad, sms_eight) != cudaSuccess) {
    rtErr("Error setting kdrsConjGrad __shared__ memory bank size to eight bytes.",
          "reductionKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kdrdConjGrad, sms_eight) != cudaSuccess) {
    rtErr("Error setting kdrdConjGrad __shared__ memory bank size to eight bytes.",
          "reductionKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kfgtConjGrad, sms_eight) != cudaSuccess) {
    rtErr("Error setting kfgtConjGrad __shared__ memory bank size to eight bytes.",
          "reductionKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kfscConjGrad, sms_eight) != cudaSuccess) {
    rtErr("Error setting kfscConjGrad __shared__ memory bank size to eight bytes.",
          "reductionKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kfrsConjGrad, sms_eight) != cudaSuccess) {
    rtErr("Error setting kfrsConjGrad __shared__ memory bank size to eight bytes.",
          "reductionKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kfrdConjGrad, sms_eight) != cudaSuccess) {
    rtErr("Error setting kfrdConjGrad __shared__ memory bank size to eight bytes.",
          "reductionKernelSetup");
  }
}

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes queryReductionKernelRequirements(const PrecisionModel prec,
                                                           const ReductionGoal purpose,
                                                           const ReductionStage process) {
  cudaFuncAttributes result;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    switch (purpose) {
    case ReductionGoal::NORMALIZE:
    case ReductionGoal::CENTER_ON_ZERO:
      break;
    case ReductionGoal::CONJUGATE_GRADIENT:
      switch (process) {
      case ReductionStage::GATHER:
        if (cudaFuncGetAttributes(&result, kdgtConjGrad) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel kdgtConjGrad.",
                "queryReductionKernelRequirements");
        }
        break;
      case ReductionStage::SCATTER:
        if (cudaFuncGetAttributes(&result, kdscConjGrad) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel kdscConjGrad.",
                "queryReductionKernelRequirements");
        }
        break;
      case ReductionStage::RESCALE:
        if (cudaFuncGetAttributes(&result, kdrsConjGrad) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel kdrsConjGrad.",
                "queryReductionKernelRequirements");
        }
        break;
      case ReductionStage::ALL_REDUCE:
        if (cudaFuncGetAttributes(&result, kdrdConjGrad) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel kdrdConjGrad.",
                "queryReductionKernelRequirements");
        }
        break;
      }
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (purpose) {
    case ReductionGoal::NORMALIZE:
    case ReductionGoal::CENTER_ON_ZERO:
      break;
    case ReductionGoal::CONJUGATE_GRADIENT:
      switch (process) {
      case ReductionStage::GATHER:
        if (cudaFuncGetAttributes(&result, kfgtConjGrad) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel kfgtConjGrad.",
                "queryReductionKernelRequirements");
        }
        break;
      case ReductionStage::SCATTER:
        if (cudaFuncGetAttributes(&result, kfscConjGrad) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel kfscConjGrad.",
                "queryReductionKernelRequirements");
        }
        break;
      case ReductionStage::RESCALE:
        if (cudaFuncGetAttributes(&result, kfrsConjGrad) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel kfrsConjGrad.",
                "queryReductionKernelRequirements");
        }
        break;
      case ReductionStage::ALL_REDUCE:
        if (cudaFuncGetAttributes(&result, kfrdConjGrad) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel kfrdConjGrad.",
                "queryReductionKernelRequirements");
        }
        break;
      }
      break;
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
extern void launchConjugateGradient(const ReductionKit &redk, ConjGradSubstrate *cgsbs,
                                    MMControlKit<double> *ctrl, const int2 bt) {

  // All conjugate gradient kernels take the same launch parameters.
  switch (redk.rps) {
  case RdwuPerSystem::ONE:
    kdrdConjGrad<<<bt.x, bt.y>>>(redk, *cgsbs, *ctrl);
    break;
  case RdwuPerSystem::MULTIPLE:
    kdgtConjGrad<<<bt.x, bt.y>>>(redk, *cgsbs, *ctrl);
    kdscConjGrad<<<bt.x, bt.y>>>(redk, *cgsbs, *ctrl);
    kdrsConjGrad<<<bt.x, bt.y>>>(redk, *cgsbs, *ctrl);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchConjugateGradient(const ReductionKit &redk, ConjGradSubstrate *cgsbs,
                                    MMControlKit<float> *ctrl, const int2 bt) {

  // All conjugate gradient kernels take the same launch parameters.
  switch (redk.rps) {
  case RdwuPerSystem::ONE:
    kfrdConjGrad<<<bt.x, bt.y>>>(redk, *cgsbs, *ctrl);
    break;
  case RdwuPerSystem::MULTIPLE:
    kfgtConjGrad<<<bt.x, bt.y>>>(redk, *cgsbs, *ctrl);
    kfscConjGrad<<<bt.x, bt.y>>>(redk, *cgsbs, *ctrl);
    kfrsConjGrad<<<bt.x, bt.y>>>(redk, *cgsbs, *ctrl);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchConjugateGradient(const PrecisionModel prec, const AtomGraphSynthesis poly_ag,
                                    PhaseSpaceSynthesis *poly_ps, ReductionBridge *rbg,
                                    MolecularMechanicsControls *mmctrl,
                                    const CoreKlManager &launcher) {
  ReductionKit redk(poly_ag, HybridTargetLevel::DEVICE);
  ConjGradSubstrate cgsbs(poly_ps, rbg, HybridTargetLevel::DEVICE);
  const int2 bt = launcher.getReductionKernelDims(prec, ReductionGoal::CONJUGATE_GRADIENT,
                                                  ReductionStage::ALL_REDUCE);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      MMControlKit<double> ctrl = mmctrl->dpData();
      launchConjugateGradient(redk, &cgsbs, &ctrl, bt);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      MMControlKit<float> ctrl = mmctrl->spData();
      launchConjugateGradient(redk, &cgsbs, &ctrl, bt);
    }
    break;
  }
}

} // namespace synthesis
} // namespace stormm
