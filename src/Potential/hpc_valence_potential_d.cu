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

// Double-precision floating point definitions
#define TCALC double
#  define VALENCE_BLOCK_MULTIPLICITY  2
#  define TCALC2 double2
#  define TCALC3 double3
#  define TCALC4 double4_16a
#  define LLCONV_FUNC __double2ll_rn
#  define SPLITCONV_FUNC doubleToInt95
#  define SPLIT_TYPE int95_t
#  define SQRT_FUNC sqrt
#  define CBRT_FUNC cbrt
#  define ACOS_FUNC acos
#  define COS_FUNC  cos
#  define SIN_FUNC  sin
#  define ABS_FUNC  fabs
#  define MIX_FUNC  computeRestraintMixtureD
#  define SPLIT_FORCE_ACCUMULATION

// Compile the standard kernels with all combinations of energy, and force accumulation methods.
#  define COMPUTE_FORCE
#    define VALENCE_KERNEL_THREAD_COUNT 256
#      define KERNEL_NAME kdsValenceForceAccumulation
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define VALENCE_KERNEL_THREAD_COUNT 192
#      define UPDATE_ATOMS
#        define KERNEL_NAME kdsValenceAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef UPDATE_ATOMS
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define COMPUTE_ENERGY
#      define VALENCE_KERNEL_THREAD_COUNT 256
#        define KERNEL_NAME kdsValenceForceEnergyAccumulation
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define UPDATE_ATOMS
#        define VALENCE_KERNEL_THREAD_COUNT 192
#          define KERNEL_NAME kdsValenceEnergyAtomUpdate
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef UPDATE_ATOMS
#    undef  COMPUTE_ENERGY
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define VALENCE_KERNEL_THREAD_COUNT 256
#      define KERNEL_NAME kdsValenceEnergyAccumulation
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef VALENCE_KERNEL_THREAD_COUNT
#  undef COMPUTE_ENERGY

// Make new kernels with a clash forgiveness check.
#  define CLASH_FORGIVENESS
#    define COMPUTE_FORCE
#      define VALENCE_KERNEL_THREAD_COUNT 256
#        define KERNEL_NAME kdsValenceForceAccumulationNonClash
#          include "valence_potential.cui"
#        undef KERNEL_NAME  
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 192
#        define UPDATE_ATOMS
#          define KERNEL_NAME kdsValenceAtomUpdateNonClash
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef UPDATE_ATOMS
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT 256
#          define KERNEL_NAME kdsValenceForceEnergyAccumulationNonClash
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define UPDATE_ATOMS
#          define VALENCE_KERNEL_THREAD_COUNT 192
#            define KERNEL_NAME kdsValenceEnergyAtomUpdateNonClash
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#          undef VALENCE_KERNEL_THREAD_COUNT
#        undef UPDATE_ATOMS
#      undef  COMPUTE_ENERGY
#    undef COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define VALENCE_KERNEL_THREAD_COUNT 256
#        define KERNEL_NAME kdsValenceEnergyAccumulationNonClash
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#    undef COMPUTE_ENERGY
#  undef CLASH_FORGIVENESS

// Clear double-precision floating point definitions
#  undef VALENCE_BLOCK_MULTIPLICITY
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
#  undef SPLIT_FORCE_ACCUMULATION
#undef TCALC

//-------------------------------------------------------------------------------------------------
#ifdef STORMM_USE_CUDA
extern cudaFuncAttributes queryValenceKernelRequirementsD(const EvaluateForce eval_frc,
                                                          const EvaluateEnergy eval_nrg,
                                                          const VwuGoal purpose,
                                                          const ClashResponse collision_handling) {
  
  // The kernel manager will have information about the GPU to use--look at the work units from
  // the perspective of overall occupancy on the GPU.
  cudaFuncAttributes result;
  cudaError_t cfa = cudaErrorInvalidValue;
  switch (collision_handling) {
  case ClashResponse::NONE:
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (purpose) {
        case VwuGoal::ACCUMULATE:
          cfa = cudaFuncGetAttributes(&result, kdsValenceForceEnergyAccumulation);
          break;
        case VwuGoal::MOVE_PARTICLES:
          cfa = cudaFuncGetAttributes(&result, kdsValenceEnergyAtomUpdate);
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (purpose) {
        case VwuGoal::ACCUMULATE:
          cfa = cudaFuncGetAttributes(&result, kdsValenceForceAccumulation);
          break;
        case VwuGoal::MOVE_PARTICLES:
          cfa = cudaFuncGetAttributes(&result, kdsValenceAtomUpdate);
          break;
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      cfa = cudaFuncGetAttributes(&result, kdsValenceEnergyAccumulation);
      break;
    }
    break;
  case ClashResponse::FORGIVE:
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (purpose) {
        case VwuGoal::ACCUMULATE:
          cfa = cudaFuncGetAttributes(&result, kdsValenceForceEnergyAccumulationNonClash);
          break;
        case VwuGoal::MOVE_PARTICLES:
          cfa = cudaFuncGetAttributes(&result, kdsValenceEnergyAtomUpdateNonClash);
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (purpose) {
        case VwuGoal::ACCUMULATE:
          cfa = cudaFuncGetAttributes(&result, kdsValenceForceAccumulationNonClash);
          break;
        case VwuGoal::MOVE_PARTICLES:
          cfa = cudaFuncGetAttributes(&result, kdsValenceAtomUpdateNonClash);
          break;
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      cfa = cudaFuncGetAttributes(&result, kdsValenceEnergyAccumulationNonClash);
      break;
    }
    break;
  }
  
  // Check for errors
  if (cfa != cudaSuccess) {

    // Construct the appropriate error message
    std::string error_message("Error obtaining attributes for kernel kds");
    error_message += "Valence";
    switch (eval_frc) {
    case EvaluateForce::YES:
      switch (eval_nrg) {
      case EvaluateEnergy::YES:
        switch (purpose) {
        case VwuGoal::ACCUMULATE:
          error_message += "ForceEnergyAccumulation";
          break;
        case VwuGoal::MOVE_PARTICLES:
          error_message += "EnergyAtomUpdate";
          break;
        }
        break;
      case EvaluateEnergy::NO:
        switch (purpose) {
        case VwuGoal::ACCUMULATE:
          error_message += "ForceAccumulation";
          break;
        case VwuGoal::MOVE_PARTICLES:
          error_message += "AtomUpdate";
          break;
        }
        break;
      }
      break;
    case EvaluateForce::NO:
      error_message += "EnergyAccumulation";
      break;
    }
    error_message += ".";

    // Report the error
    rtErr(error_message, "queryValenceKernelRequirementsD");
  }
  
  return result;
}
#endif

//-------------------------------------------------------------------------------------------------
extern void launchValence(const SyValenceKit<double> &poly_vk,
                          const SyRestraintKit<double, double2, double4_16a> &poly_rk,
                          MMControlKit<double> *ctrl, PsSynthesisWriter *poly_psw,
                          const SyAtomUpdateKit<double, double2, double4_16a> &poly_auk,
                          ThermostatWriter<double> *tstw, ScoreCardWriter *scw,
                          CacheResourceKit<double> *gmem_r, const EvaluateForce eval_force,
                          const EvaluateEnergy eval_energy, const VwuGoal purpose, const int2 bt,
                          const double clash_distance, const double clash_ratio) {
  
  // Rather than a switch over cases of the ClashResponse enumerator, just use the nonzero values
  // of either parameter to indicate that clash damping has been requested.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (purpose) {
    case VwuGoal::ACCUMULATE:

      // When the goal is to accumulate energies, forces, or both, the force accumulation method
      // is set to use int64 data.  A 95-bit method that splits the accumulation with overflow into
      // a secondary 32-bit int may be added, and likewise become the sole option for
      // double-precision computations.
      switch (eval_force) {
      case EvaluateForce::YES:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          kdsValenceForceEnergyAccumulationNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                    *poly_psw, clash_distance,
                                                                    clash_ratio, *scw, *gmem_r);
          break;
        case EvaluateEnergy::NO:
          kdsValenceForceAccumulationNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                              clash_distance, clash_ratio,
                                                              *gmem_r);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdsValenceEnergyAccumulationNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                             clash_distance, clash_ratio, *scw,
                                                             *gmem_r);
        break;
      }
     break;
    case VwuGoal::MOVE_PARTICLES:

      // When the goal is to move particles, evaluating the force is obligatory, but the manner in
      // which forces are accumulated is still important.  Whether to accumulate energies while
      // evaluating forces and moving the particles remains a consideration in choosing the proper
      // kernel.
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        kdsValenceEnergyAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           clash_distance, clash_ratio, poly_auk,
                                                           *tstw, *scw, *gmem_r);
        break;
      case EvaluateEnergy::NO:
        kdsValenceAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                     clash_distance, clash_ratio, poly_auk, *tstw,
                                                     *gmem_r);
        break;
      }
      break;
    }
  }
  else {
    switch (purpose) {
    case VwuGoal::ACCUMULATE:
      switch (eval_force) {
      case EvaluateForce::YES:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          kdsValenceForceEnergyAccumulation<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                            *scw, *gmem_r);
          break;
        case EvaluateEnergy::NO:
          kdsValenceForceAccumulation<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *gmem_r);
          break;
        }
        break;
      case EvaluateForce::NO:
        kdsValenceEnergyAccumulation<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                     *gmem_r);
        break;
      }
     break;
    case VwuGoal::MOVE_PARTICLES:
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        kdsValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, poly_auk,
                                                   *tstw, *scw, *gmem_r);
        break;
      case EvaluateEnergy::NO:
        kdsValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, poly_auk, *tstw,
                                             *gmem_r);
        break;
      }
      break;
    }
  }
}

} // namespace energy
} // namespace stormm
