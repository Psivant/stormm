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
#  define VALENCE_BLOCK_MULTIPLICITY 2
#  define UPDATE_ATOMS
#  define COMPUTE_FORCE
#  define TCOORD double
#  define TACC llint
#  define TCOORD4 double4_16a
#  define TCOORD_IS_LONG
#    define SPLIT_FORCE_ACCUMULATION
#      define VALENCE_KERNEL_THREAD_COUNT 384
#        define KERNEL_NAME kfsdPmeValenceAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT 320
#          define KERNEL_NAME kfsdPmeValenceEnergyAtomUpdate
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef COMPUTE_ENERGY
#    undef SPLIT_FORCE_ACCUMULATION
#    define VALENCE_KERNEL_THREAD_COUNT 384
#      define KERNEL_NAME kfdPmeValenceAtomUpdate
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define COMPUTE_ENERGY
#      define VALENCE_KERNEL_THREAD_COUNT 320
#        define KERNEL_NAME kfdPmeValenceEnergyAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#    undef COMPUTE_ENERGY
  
// Make new kernels with a clash forgiveness check.
#    define CLASH_FORGIVENESS
#      define SPLIT_FORCE_ACCUMULATION
#        define VALENCE_KERNEL_THREAD_COUNT 384
#          define KERNEL_NAME kfsdPmeValenceAtomUpdateNonClash
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define COMPUTE_ENERGY
#          define VALENCE_KERNEL_THREAD_COUNT 320
#            define KERNEL_NAME kfsdPmeValenceEnergyAtomUpdateNonClash
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#          undef VALENCE_KERNEL_THREAD_COUNT
#        undef COMPUTE_ENERGY
#      undef SPLIT_FORCE_ACCUMULATION
#      define VALENCE_KERNEL_THREAD_COUNT 384
#        define KERNEL_NAME kfdPmeValenceAtomUpdateNonClash
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT 320
#          define KERNEL_NAME kfdPmeValenceEnergyAtomUpdateNonClash
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
#        define KERNEL_NAME kfsfPmeValenceAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT 320
#          define KERNEL_NAME kfsfPmeValenceEnergyAtomUpdate
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef COMPUTE_ENERGY
#    undef SPLIT_FORCE_ACCUMULATION
#    define VALENCE_KERNEL_THREAD_COUNT 384
#      define KERNEL_NAME kffPmeValenceAtomUpdate
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define COMPUTE_ENERGY
#      define VALENCE_KERNEL_THREAD_COUNT 320
#        define KERNEL_NAME kffPmeValenceEnergyAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#    undef COMPUTE_ENERGY
  
// Make new kernels with a clash forgiveness check.
#    define CLASH_FORGIVENESS
#      define SPLIT_FORCE_ACCUMULATION
#        define VALENCE_KERNEL_THREAD_COUNT 384
#          define KERNEL_NAME kfsfPmeValenceAtomUpdateNonClash
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define COMPUTE_ENERGY
#          define VALENCE_KERNEL_THREAD_COUNT 320
#            define KERNEL_NAME kfsfPmeValenceEnergyAtomUpdateNonClash
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#          undef VALENCE_KERNEL_THREAD_COUNT
#        undef COMPUTE_ENERGY
#      undef SPLIT_FORCE_ACCUMULATION
#      define VALENCE_KERNEL_THREAD_COUNT 384
#        define KERNEL_NAME kffPmeValenceAtomUpdateNonClash
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT 320
#          define KERNEL_NAME kffPmeValenceEnergyAtomUpdateNonClash
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
queryValenceKernelRequirementsFPME(const EvaluateEnergy eval_nrg,
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
          cfa = cudaFuncGetAttributes(&result, kfsdPmeValenceEnergyAtomUpdate);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kfsfPmeValenceEnergyAtomUpdate);
          break;
        }
        break;
      case AccumulationMethod::WHOLE:
        switch (neighbor_prec) {
        case PrecisionModel::DOUBLE:
          cfa = cudaFuncGetAttributes(&result, kfdPmeValenceEnergyAtomUpdate);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kffPmeValenceEnergyAtomUpdate);
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
          cfa = cudaFuncGetAttributes(&result, kfsdPmeValenceAtomUpdate);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kfsfPmeValenceAtomUpdate);
          break;
        }
        break;
      case AccumulationMethod::WHOLE:
        switch (neighbor_prec) {
        case PrecisionModel::DOUBLE:
          cfa = cudaFuncGetAttributes(&result, kfdPmeValenceAtomUpdate);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kfdPmeValenceAtomUpdate);
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
          cfa = cudaFuncGetAttributes(&result, kfsdPmeValenceEnergyAtomUpdateNonClash);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kfsfPmeValenceEnergyAtomUpdateNonClash);
          break;
        }
        break;
      case AccumulationMethod::WHOLE:
        switch (neighbor_prec) {
        case PrecisionModel::DOUBLE:
          cfa = cudaFuncGetAttributes(&result, kfdPmeValenceEnergyAtomUpdateNonClash);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kffPmeValenceEnergyAtomUpdateNonClash);
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
          cfa = cudaFuncGetAttributes(&result, kfsdPmeValenceAtomUpdateNonClash);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kfsfPmeValenceAtomUpdateNonClash);
          break;
        }
        break;
      case AccumulationMethod::WHOLE:
        switch (neighbor_prec) {
        case PrecisionModel::DOUBLE:
          cfa = cudaFuncGetAttributes(&result, kfdPmeValenceAtomUpdateNonClash);
          break;
        case PrecisionModel::SINGLE:
          cfa = cudaFuncGetAttributes(&result, kffPmeValenceAtomUpdateNonClash);
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
            "queryValenceKernelRequirementsFPME");
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
    rtErr(error_message, "queryValenceKernelRequirementsFPME");
  }
  
  return result;
}
#endif

//-------------------------------------------------------------------------------------------------
extern void launchValence(const SyValenceKit<float> &poly_vk,
                          const SyRestraintKit<float, float2, float4> &poly_rk,
                          const CellGridReader<double, llint, double, double4_16a> &cgr,
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
          kfsdPmeValenceEnergyAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                 *poly_psw, cgr, clash_distance,
                                                                 clash_ratio, poly_auk, *tstw,
                                                                 *scw, *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kfdPmeValenceEnergyAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                *poly_psw, cgr, clash_distance,
                                                                clash_ratio, poly_auk, *tstw,
                                                                *scw, *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (force_sum) {
        case AccumulationMethod::SPLIT:
          kfsdPmeValenceAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr,
                                                           clash_distance, clash_ratio, poly_auk,
                                                           *tstw, *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kfdPmeValenceAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr,
                                                          clash_distance, clash_ratio, poly_auk,
                                                          *tstw, *gmem_r);
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
          kfsdPmeValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr,
                                                         poly_auk, *tstw, *scw, *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kfdPmeValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr,
                                                        poly_auk, *tstw, *scw, *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (force_sum) {
        case AccumulationMethod::SPLIT:
          kfsdPmeValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr,
                                                   poly_auk, *tstw, *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kfdPmeValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr,
                                                  poly_auk, *tstw, *gmem_r);
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
                          const CellGridReader<float, int, float, float4> &cgr,
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
          kfsfPmeValenceEnergyAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                 *poly_psw, cgr, clash_distance,
                                                                 clash_ratio, poly_auk, *tstw,
                                                                 *scw, *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kffPmeValenceEnergyAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                *poly_psw, cgr, clash_distance,
                                                                clash_ratio, poly_auk, *tstw,
                                                                *scw, *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (force_sum) {
        case AccumulationMethod::SPLIT:
          kfsfPmeValenceAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr,
                                                           clash_distance, clash_ratio, poly_auk,
                                                           *tstw, *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kffPmeValenceAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr,
                                                          clash_distance, clash_ratio, poly_auk,
                                                          *tstw, *gmem_r);
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
          kfsfPmeValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr,
                                                         poly_auk, *tstw, *scw, *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kffPmeValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr,
                                                        poly_auk, *tstw, *scw, *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (force_sum) {
        case AccumulationMethod::SPLIT:
          kfsfPmeValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr,
                                                   poly_auk, *tstw, *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          kffPmeValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr,
                                                  poly_auk, *tstw, *gmem_r);
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
