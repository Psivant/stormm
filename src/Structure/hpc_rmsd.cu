// -*-c++-*-
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Constants/hpc_bounds.h"
#include "Math/matrix_ops.h"
#include "Numerics/split_fixed_precision.h"
#include "Reporting/error_format.h"
#include "hpc_rmsd.h"
#include "rmsd_plan.h"

namespace stormm {
namespace structure {

using analysis::CompGuideKit;
using constants::PrecisionModel;
using constants::warp_bits;
using constants::warp_bits_mask_int;
using constants::warp_size_int;
using numerics::max_llint_accumulation;
using stmath::maximum_ql_iterations;

//-------------------------------------------------------------------------------------------------
#define TCALC double
#define TCALC4 double4
#define TCALC_IS_DOUBLE
#define RMSD_REF_KERNEL_NAME kdComputeRMSDToReference
#define RMSD_MAT_KERNEL_NAME kdComputeRMSDMatrix
#define SQRT_FUNC sqrt
#define FABS_FUNC fabs
#include "rmsd_calculation.cui"
#undef TCALC
#undef TCALC4
#undef RMSD_REF_KERNEL_NAME
#undef RMSD_MAT_KERNEL_NAME
#undef SQRT_FUNC
#undef FABS_FUNC

#define TCALC float
#define TCALC4 float4
#define RMSD_REF_KERNEL_NAME kfComputeRMSDToReference
#define RMSD_MAT_KERNEL_NAME kfComputeRMSDMatrix
#define SQRT_FUNC sqrtf
#define FABS_FUNC fabsf
#include "rmsd_calculation.cui"
#undef TCALC
#undef TCALC4
#undef RMSD_REF_KERNEL_NAME
#undef RMSD_MAT_KERNEL_NAME
#undef SQRT_FUNC
#undef FABS_FUNC

//-------------------------------------------------------------------------------------------------
cudaFuncAttributes queryRMSDKernelRequirements(const PrecisionModel prec, const RMSDTask order) {
  cudaFuncAttributes result;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    switch (order) {
    case RMSDTask::REFERENCE:
      if (cudaFuncGetAttributes(&result, kdComputeRMSDToReference) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel kdComputeRMSDToReference.",
              "queryRMSDKernelRequirements");
      }
      break;
    case RMSDTask::MATRIX:
      if (cudaFuncGetAttributes(&result, kdComputeRMSDMatrix) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel kdComputeRMSDMatrix.",
              "queryRMSDKernelRequirements");
      }
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (order) {
    case RMSDTask::REFERENCE:
      if (cudaFuncGetAttributes(&result, kfComputeRMSDToReference) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel kfComputeRMSDToReference.",
              "queryRMSDKernelRequirements");
      }
      break;
    case RMSDTask::MATRIX:
      if (cudaFuncGetAttributes(&result, kfComputeRMSDMatrix) != cudaSuccess) {
        rtErr("Error obtaining attributes for kernel kfComputeRMSDMatrix.",
              "queryRMSDKernelRequirements");
      }
      break;
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &poly_ps,
          const Hybrid<int> &reference_frames, Hybrid<double> *result,
          const CoreKlManager &launcher) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  const int2 lp = launcher.getRMSDKernelDims(PrecisionModel::DOUBLE, RMSDTask::REFERENCE);
  kdComputeRMSDToReference<<<lp.x, lp.y>>>(cg.data(tier), rplan.dpData(tier), poly_ps.data(tier),
                                           reference_frames.data(tier), result->data(tier));
}

//-------------------------------------------------------------------------------------------------
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &poly_ps,
          const Hybrid<int> &reference_frames, Hybrid<float> *result,
          const CoreKlManager &launcher) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  const int2 lp = launcher.getRMSDKernelDims(PrecisionModel::SINGLE, RMSDTask::REFERENCE);
  kfComputeRMSDToReference<<<lp.x, lp.y>>>(cg.data(tier), rplan.spData(tier), poly_ps.data(tier),
                                           reference_frames.data(tier), result->data(tier));
}

//-------------------------------------------------------------------------------------------------
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &poly_ps,
          Hybrid<double> *result, const CoreKlManager &launcher) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  const int2 lp = launcher.getRMSDKernelDims(PrecisionModel::DOUBLE, RMSDTask::MATRIX);
  kdComputeRMSDMatrix<<<lp.x, lp.y>>>(cg.data(tier), rplan.dpData(tier), poly_ps.data(tier),
                                      result->data(tier));
}

//-------------------------------------------------------------------------------------------------
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &poly_ps,
          Hybrid<float> *result, const CoreKlManager &launcher) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  const int2 lp = launcher.getRMSDKernelDims(PrecisionModel::SINGLE, RMSDTask::MATRIX);
  kfComputeRMSDMatrix<<<lp.x, lp.y>>>(cg.data(tier), rplan.spData(tier), poly_ps.data(tier),
                                      result->data(tier));
}

} // namespace structure
} // namespace stormm
