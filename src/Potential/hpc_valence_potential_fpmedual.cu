// -*-c++-*-
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Accelerator/gpu_details.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/rounding.h"
#include "Numerics/numeric_enumerators.h"
#include "Numerics/split_fixed_precision.h"
#include "Potential/cellgrid.h"
#include "Potential/energy_enumerators.h"
#include "Synthesis/valence_workunit.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/trajectory_enumerators.h"
#include "hpc_valence_potential.h"

namespace stormm {
namespace energy {

using constants::large_block_size;
using constants::medium_block_size;
using constants::small_block_size;
using constants::twice_warp_bits_mask_int;
using constants::twice_warp_size_int;
using constants::warp_size_int;
using constants::warp_bits;
using constants::warp_bits_mask_int;
using numerics::chooseAccumulationMethod;
using numerics::getEnumerationName;
using stmath::roundUp;
using symbols::asymptotic_to_one_f;
using symbols::asymptotic_to_one_lf;
using symbols::boltzmann_constant_f;
using symbols::gafs_to_kcal_f;
using symbols::inverse_one_minus_asymptote_f;
using symbols::inverse_one_minus_asymptote_lf;
using symbols::inverse_twopi_f;
using symbols::kcal_to_gafs_f;
using symbols::near_to_one_f;
using symbols::near_to_one_lf;
using symbols::pi;
using symbols::pi_f;
using symbols::twopi;
using symbols::twopi_f;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::half_valence_work_unit_atoms;
using synthesis::quarter_valence_work_unit_atoms;
using synthesis::eighth_valence_work_unit_atoms;
using synthesis::VwuAbstractMap;
using synthesis::vwu_abstract_length;
using trajectory::ThermostatKind;
using trajectory::ThermostatPartition;
using topology::TorsionKind;
using topology::VirtualSiteKind;

#include "Accelerator/syncwarp.cui"
#include "Math/rounding.cui"
#include "Math/vector_formulas.cui"
#include "Numerics/accumulation.cui"
#include "Trajectory/thermostat_utilities.cui"
#include "valence_util.cui"

// Single-precision floating point definitions
#define TCALC float
#  define TCALC2 float2
#  define TCALC3 float3
#  define TCALC4 float4
#  define LLCONV_FUNC __float2ll_rn
#  define SPLITCONV_FUNC floatToInt63
#  define SPLIT_TYPE int2
#  define SQRT_FUNC sqrtf
#  define CBRT_FUNC cbrtf
#  define ACOS_FUNC acosf
#  define COS_FUNC  cosf
#  define SIN_FUNC  sinf
#  define ABS_FUNC  fabsf
#  define MIX_FUNC  computeRestraintMixtureF
#  define TCALC_IS_SINGLE
  
// Make additional kernels for PME-based simulations, where the atom updates will also draw upon a
// neighbor list cell grid, or even a pair of such grids.
#  define PME_COMPATIBLE
#  define DUAL_GRIDS
#  define VALENCE_BLOCK_MULTIPLICITY 2
#  define UPDATE_ATOMS
#  define COMPUTE_FORCE
#  define TCOORD double
#  define TACC llint
#  define TCOORD4 double4_16a
#  define TCOORD_IS_LONG
#    define SPLIT_FORCE_ACCUMULATION
#      define VALENCE_KERNEL_THREAD_COUNT 384
#        define KERNEL_NAME kfsdPmeDualValenceAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT 320
#          define KERNEL_NAME kfsdPmeDualValenceEnergyAtomUpdate
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef COMPUTE_ENERGY
#    undef SPLIT_FORCE_ACCUMULATION
#    define VALENCE_KERNEL_THREAD_COUNT 384
#      define KERNEL_NAME kfdPmeDualValenceAtomUpdate
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define COMPUTE_ENERGY
#      define VALENCE_KERNEL_THREAD_COUNT 320
#        define KERNEL_NAME kfdPmeDualValenceEnergyAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#    undef COMPUTE_ENERGY
  
// Make new kernels with a clash forgiveness check.
#    define CLASH_FORGIVENESS
#      define SPLIT_FORCE_ACCUMULATION
#        define VALENCE_KERNEL_THREAD_COUNT 384
#          define KERNEL_NAME kfsdPmeDualValenceAtomUpdateNonClash
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define COMPUTE_ENERGY
#          define VALENCE_KERNEL_THREAD_COUNT 320
#            define KERNEL_NAME kfsdPmeDualValenceEnergyAtomUpdateNonClash
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#          undef VALENCE_KERNEL_THREAD_COUNT
#        undef COMPUTE_ENERGY
#      undef SPLIT_FORCE_ACCUMULATION
#      define VALENCE_KERNEL_THREAD_COUNT 384
#        define KERNEL_NAME kfdPmeDualValenceAtomUpdateNonClash
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT 320
#          define KERNEL_NAME kfdPmeDualValenceEnergyAtomUpdateNonClash
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef COMPUTE_ENERGY
#    undef CLASH_FORGIVENESS
#  undef TCOORD_IS_LONG
#  undef TCOORD4
#  undef TACC
#  undef TCOORD

#  define TCOORD float
#  define TACC int
#  define TCOORD4 float4
#    define SPLIT_FORCE_ACCUMULATION
#      define VALENCE_KERNEL_THREAD_COUNT 384
#        define KERNEL_NAME kfsfPmeDualValenceAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT 320
#          define KERNEL_NAME kfsfPmeDualValenceEnergyAtomUpdate
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef COMPUTE_ENERGY
#    undef SPLIT_FORCE_ACCUMULATION
#    define VALENCE_KERNEL_THREAD_COUNT 384
#      define KERNEL_NAME kffPmeDualValenceAtomUpdate
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define COMPUTE_ENERGY
#      define VALENCE_KERNEL_THREAD_COUNT 320
#        define KERNEL_NAME kffPmeDualValenceEnergyAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#    undef COMPUTE_ENERGY
  
// Make new kernels with a clash forgiveness check.
#    define CLASH_FORGIVENESS
#      define SPLIT_FORCE_ACCUMULATION
#        define VALENCE_KERNEL_THREAD_COUNT 384
#          define KERNEL_NAME kfsfPmeDualValenceAtomUpdateNonClash
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define COMPUTE_ENERGY
#          define VALENCE_KERNEL_THREAD_COUNT 320
#            define KERNEL_NAME kfsfPmeDualValenceEnergyAtomUpdateNonClash
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#          undef VALENCE_KERNEL_THREAD_COUNT
#        undef COMPUTE_ENERGY
#      undef SPLIT_FORCE_ACCUMULATION
#      define VALENCE_KERNEL_THREAD_COUNT 384
#        define KERNEL_NAME kffPmeDualValenceAtomUpdateNonClash
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT 320
#          define KERNEL_NAME kffPmeDualValenceEnergyAtomUpdateNonClash
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef COMPUTE_ENERGY
#    undef CLASH_FORGIVENESS
#  undef TCOORD4
#  undef TACC
#  undef TCOORD

#  undef COMPUTE_FORCE
#  undef UPDATE_ATOMS
#  undef VALENCE_BLOCK_MULTIPLICITY
#  undef DUAL_GRIDS
#  undef PME_COMPATIBLE

// Clear single-precision floating point definitions
#  undef TCALC2
#  undef TCALC3
#  undef TCALC4
#  undef LLCONV_FUNC
#  undef SPLITCONV_FUNC
#  undef SPLIT_TYPE
#  undef SQRT_FUNC
#  undef CBRT_FUNC
#  undef ACOS_FUNC
#  undef COS_FUNC
#  undef SIN_FUNC
#  undef ABS_FUNC
#  undef MIX_FUNC
#  undef TCALC_IS_SINGLE
#undef TCALC

//-------------------------------------------------------------------------------------------------
#ifdef STORMM_USE_CUDA
extern cudaFuncAttributes
queryValenceKernelRequirementsFPMEDual(const EvaluateEnergy eval_nrg,
                                       const AccumulationMethod acc_meth,
                                       const ClashResponse collision_handling,
                                       const PrecisionModel neighbor_prec) {
  
  // The kernel manager will have information about the GPU to use--look at the work units from
  // the perspective of overall occupancy on the GPU.
  cudaFuncAttributes result;
  cudaError_t cfa = cudaErrorInvalidValue;
  switch (collision_handling) {
  case ClashResponse::NONE:
    switch (eval_nrg) {
    case EvaluateEnergy::YES:
      switch (acc_meth) {
      case AccumulationMethod::SPLIT:
        switch (neighbor_prec) {
        case PrecisionModel::DOUBLE:
          cfa = cudaFuncGetAttributes(&result, kfsdPmeDualValenceEnergyAtomUpdate);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kfsfPmeDualValenceEnergyAtomUpdate);
          break;
        }
        break;
      case AccumulationMethod::WHOLE:
        switch (neighbor_prec) {
        case PrecisionModel::DOUBLE:
          cfa = cudaFuncGetAttributes(&result, kfdPmeDualValenceEnergyAtomUpdate);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kffPmeDualValenceEnergyAtomUpdate);
          break;
        }
        break;
      }
      break;
    case EvaluateEnergy::NO:
      switch (acc_meth) {
      case AccumulationMethod::SPLIT:
        switch (neighbor_prec) {
        case PrecisionModel::DOUBLE:
          cfa = cudaFuncGetAttributes(&result, kfsdPmeDualValenceAtomUpdate);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kfsfPmeDualValenceAtomUpdate);
          break;
        }
        break;
      case AccumulationMethod::WHOLE:
        switch (neighbor_prec) {
        case PrecisionModel::DOUBLE:
          cfa = cudaFuncGetAttributes(&result, kfdPmeDualValenceAtomUpdate);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kfdPmeDualValenceAtomUpdate);
          break;
        }
        break;
      }
      break;
    }
    break;
  case ClashResponse::FORGIVE:
    switch (eval_nrg) {
    case EvaluateEnergy::YES:
      switch (acc_meth) {
      case AccumulationMethod::SPLIT:
        switch (neighbor_prec) {
        case PrecisionModel::DOUBLE:
          cfa = cudaFuncGetAttributes(&result, kfsdPmeDualValenceEnergyAtomUpdateNonClash);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kfsfPmeDualValenceEnergyAtomUpdateNonClash);
          break;
        }
        break;
      case AccumulationMethod::WHOLE:
        switch (neighbor_prec) {
        case PrecisionModel::DOUBLE:
          cfa = cudaFuncGetAttributes(&result, kfdPmeDualValenceEnergyAtomUpdateNonClash);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kffPmeDualValenceEnergyAtomUpdateNonClash);
          break;
        }
        break;
      }
      break;
    case EvaluateEnergy::NO:
      switch (acc_meth) {
      case AccumulationMethod::SPLIT:
        switch (neighbor_prec) {
        case PrecisionModel::DOUBLE:
          cfa = cudaFuncGetAttributes(&result, kfsdPmeDualValenceAtomUpdateNonClash);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kfsfPmeDualValenceAtomUpdateNonClash);
          break;
        }
        break;
      case AccumulationMethod::WHOLE:
        switch (neighbor_prec) {
        case PrecisionModel::DOUBLE:
          cfa = cudaFuncGetAttributes(&result, kfdPmeDualValenceAtomUpdateNonClash);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kffPmeDualValenceAtomUpdateNonClash);
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
    std::string error_message("Error obtaining attributes for kernel kf");
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
      error_message += "s";
      break;
    case AccumulationMethod::WHOLE:
      break;
    case AccumulationMethod::AUTOMATIC:
      rtErr("Kernels do not accept " + getEnumerationName(acc_meth) + " accumulation.",
            "queryValenceKernelRequirementsFPMEDual");
      break;
    }
    switch (neighbor_prec) {
    case PrecisionModel::DOUBLE:
      error_message += "d";
      break;
    case PrecisionModel::SINGLE:
      error_message += "f";
      break;
    }
    error_message += "PmeValence";
    switch (eval_nrg) {
    case EvaluateEnergy::YES:
      error_message += "EnergyAtomUpdate";
      break;
    case EvaluateEnergy::NO:
      error_message += "AtomUpdate";
      break;
    }
    error_message += ".";

    // Report the error
    rtErr(error_message, "queryValenceKernelRequirementsFPMEDual");
  }
  
  return result;
}
#endif

//-------------------------------------------------------------------------------------------------
extern void launchValence(const SyValenceKit<float> &poly_vk,
                          const SyRestraintKit<float, float2, float4> &poly_rk,
                          const CellGridReader<double, llint, double, double4_16a> &cgr_qq,
                          const CellGridReader<double, llint, double, double4_16a> &cgr_lj,
                          MMControlKit<float> *ctrl, PsSynthesisWriter *poly_psw,
                          const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                          ThermostatWriter<float> *tstw, ScoreCardWriter *scw,
                          CacheResourceKit<float> *gmem_r, const EvaluateForce eval_force,
                          const EvaluateEnergy eval_energy, const VwuGoal purpose,
                          const AccumulationMethod force_sum, const int2 bt,
                          const float clash_distance, const float clash_ratio) {
  switch (purpose) {
  case VwuGoal::ACCUMULATE:
    launchValence(poly_vk, poly_rk, ctrl, poly_psw, poly_auk, tstw, scw, gmem_r, eval_force,
                  eval_energy, purpose, force_sum, bt, clash_distance, clash_ratio);
    break;
  case VwuGoal::MOVE_PARTICLES:
    if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        switch (force_sum) {
        case AccumulationMethod::SPLIT:
          kfsdPmeDualValenceEnergyAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                     *poly_psw, cgr_qq, cgr_lj,
                                                                     clash_distance, clash_ratio,
                                                                     poly_auk, *tstw, *scw,
                                                                     *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kfdPmeDualValenceEnergyAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                    *poly_psw, cgr_qq, cgr_lj,
                                                                    clash_distance, clash_ratio,
                                                                    poly_auk, *tstw, *scw,
                                                                    *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (force_sum) {
        case AccumulationMethod::SPLIT:
          kfsdPmeDualValenceAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                               cgr_qq, cgr_lj, clash_distance,
                                                               clash_ratio, poly_auk, *tstw,
                                                               *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kfdPmeDualValenceAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                              cgr_qq, cgr_lj, clash_distance,
                                                              clash_ratio, poly_auk, *tstw,
                                                              *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      }
    }
    else {
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        switch (force_sum) {
        case AccumulationMethod::SPLIT:
          kfsdPmeDualValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                             cgr_qq, cgr_lj, poly_auk, *tstw, *scw,
                                                             *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kfdPmeDualValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                            cgr_qq, cgr_lj, poly_auk, *tstw, *scw,
                                                            *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (force_sum) {
        case AccumulationMethod::SPLIT:
          kfsdPmeDualValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr_qq,
                                                       cgr_lj, poly_auk, *tstw, *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kfdPmeDualValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr_qq,
                                                      cgr_lj, poly_auk, *tstw, *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchValence(const SyValenceKit<float> &poly_vk,
                          const SyRestraintKit<float, float2, float4> &poly_rk,
                          const CellGridReader<float, int, float, float4> &cgr_qq,
                          const CellGridReader<float, int, float, float4> &cgr_lj,
                          MMControlKit<float> *ctrl, PsSynthesisWriter *poly_psw,
                          const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                          ThermostatWriter<float> *tstw, ScoreCardWriter *scw,
                          CacheResourceKit<float> *gmem_r, const EvaluateForce eval_force,
                          const EvaluateEnergy eval_energy, const VwuGoal purpose,
                          const AccumulationMethod force_sum, const int2 bt,
                          const float clash_distance, const float clash_ratio) {
  switch (purpose) {
  case VwuGoal::ACCUMULATE:
    launchValence(poly_vk, poly_rk, ctrl, poly_psw, poly_auk, tstw, scw, gmem_r, eval_force,
                  eval_energy, purpose, force_sum, bt, clash_distance, clash_ratio);
    break;
  case VwuGoal::MOVE_PARTICLES:
    if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        switch (force_sum) {
        case AccumulationMethod::SPLIT:
          kfsfPmeDualValenceEnergyAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                     *poly_psw, cgr_qq, cgr_lj,
                                                                     clash_distance, clash_ratio,
                                                                     poly_auk, *tstw, *scw,
                                                                     *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kffPmeDualValenceEnergyAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                    *poly_psw, cgr_qq, cgr_lj,
                                                                    clash_distance, clash_ratio,
                                                                    poly_auk, *tstw, *scw,
                                                                    *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (force_sum) {
        case AccumulationMethod::SPLIT:
          kfsfPmeDualValenceAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                               cgr_qq, cgr_lj, clash_distance,
                                                               clash_ratio, poly_auk, *tstw,
                                                               *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kffPmeDualValenceAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                              cgr_qq, cgr_lj, clash_distance,
                                                              clash_ratio, poly_auk, *tstw,
                                                              *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      }
    }
    else {
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        switch (force_sum) {
        case AccumulationMethod::SPLIT:
          kfsfPmeDualValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                             cgr_qq, cgr_lj, poly_auk, *tstw, *scw,
                                                             *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kffPmeDualValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                            cgr_qq, cgr_lj, poly_auk, *tstw, *scw,
                                                            *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (force_sum) {
        case AccumulationMethod::SPLIT:
          kfsfPmeDualValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr_qq,
                                                       cgr_lj, poly_auk, *tstw, *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kffPmeDualValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr_qq,
                                                      cgr_lj, poly_auk, *tstw, *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      }
    }
    break;
  }
}
  
} // namespace energy
} // namespace stormm
