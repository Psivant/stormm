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

// Define PME-compatible variants of kernels invoking double-precision arithmetic, beginning with
// single neighbor list cell grids.  TCALC is inherent to the valence kernel, but the coordinate
// type TCOORD and associated TACC may be defined independently.  Begin with double-precision
// coordinates and 95-bit force accumulation in the neighbor list cell grids.
#  define PME_COMPATIBLE
#  define VALENCE_BLOCK_MULTIPLICITY  2
#  define UPDATE_ATOMS
#  define COMPUTE_FORCE
#  define DUAL_GRIDS
#  define TCOORD double
#  define TACC llint
#  define TCOORD4 double4_16a
#  define TCOORD_IS_LONG
#    define VALENCE_KERNEL_THREAD_COUNT 192
#      define KERNEL_NAME kdsdPmeDualValenceAtomUpdate
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define VALENCE_KERNEL_THREAD_COUNT 192
#      define COMPUTE_ENERGY
#        define KERNEL_NAME kdsdPmeDualValenceEnergyAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef COMPUTE_ENERGY
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define CLASH_FORGIVENESS
#      define VALENCE_KERNEL_THREAD_COUNT 192
#        define KERNEL_NAME kdsdPmeDualValenceAtomUpdateNonClash
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT 192
#          define KERNEL_NAME kdsdPmeDualValenceEnergyAtomUpdateNonClash
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef  COMPUTE_ENERGY
#    undef CLASH_FORGIVENESS
#  undef TCOORD_IS_LONG
#  undef TCOORD4
#  undef TACC
#  undef TCOORD
  
// Define additional PME-compatible kernels for use with float coordinates and int63 accumulation
// in the neighbor list's non-bonded forces.
#  define TCOORD float
#  define TACC int
#  define TCOORD4 float4
#    define VALENCE_KERNEL_THREAD_COUNT 192
#      define KERNEL_NAME kdsfPmeDualValenceAtomUpdate
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define VALENCE_KERNEL_THREAD_COUNT 192
#      define COMPUTE_ENERGY
#        define KERNEL_NAME kdsfPmeDualValenceEnergyAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef  COMPUTE_ENERGY
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define CLASH_FORGIVENESS
#      define VALENCE_KERNEL_THREAD_COUNT 192
#        define KERNEL_NAME kdsfPmeDualValenceAtomUpdateNonClash
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT 192
#          define KERNEL_NAME kdsfPmeDualValenceEnergyAtomUpdateNonClash
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef  COMPUTE_ENERGY
#    undef CLASH_FORGIVENESS
#  undef TCOORD4
#  undef TACC
#  undef TCOORD

// Clear double-precision floating point definitions
#  undef VALENCE_BLOCK_MULTIPLICITY
#  undef UPDATE_ATOMS
#  undef COMPUTE_FORCE
#  undef DUAL_GRIDS
#  undef PME_COMPATIBLE
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
extern cudaFuncAttributes
queryValenceKernelRequirementsDPMEDual(const EvaluateEnergy eval_nrg,
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
      switch (neighbor_prec) {
      case PrecisionModel::DOUBLE:
        cfa = cudaFuncGetAttributes(&result, kdsdPmeDualValenceEnergyAtomUpdate);
        break;
      case PrecisionModel::SINGLE:
        cfa = cudaFuncGetAttributes(&result, kdsfPmeDualValenceEnergyAtomUpdate);
        break;
      }
      break;
    case EvaluateEnergy::NO:
      switch (neighbor_prec) {
      case PrecisionModel::DOUBLE:
        cfa = cudaFuncGetAttributes(&result, kdsdPmeDualValenceAtomUpdate);
        break;
        break;
      case PrecisionModel::SINGLE:
        cfa = cudaFuncGetAttributes(&result, kdsfPmeDualValenceAtomUpdate);
        break;
      }
      break;
    }
    break;
  case ClashResponse::FORGIVE:
    switch (eval_nrg) {
    case EvaluateEnergy::YES:
      switch (neighbor_prec) {
      case PrecisionModel::DOUBLE:
        cfa = cudaFuncGetAttributes(&result, kdsdPmeDualValenceEnergyAtomUpdateNonClash);
        break;
      case PrecisionModel::SINGLE:
        cfa = cudaFuncGetAttributes(&result, kdsfPmeDualValenceEnergyAtomUpdateNonClash);
        break;
      }
      break;
    case EvaluateEnergy::NO:
      switch (neighbor_prec) {
      case PrecisionModel::DOUBLE:
        cfa = cudaFuncGetAttributes(&result, kdsdPmeDualValenceAtomUpdateNonClash);
        break;
      case PrecisionModel::SINGLE:
        cfa = cudaFuncGetAttributes(&result, kdsfPmeDualValenceAtomUpdateNonClash);
        break;
      }
      break;
    }
    break;
  }
  
  // Check for errors
  if (cfa != cudaSuccess) {

    // Construct the appropriate error message
    std::string error_message("Error obtaining attributes for kernel kds");
    switch (neighbor_prec) {
    case PrecisionModel::DOUBLE:
      error_message += "d";
      break;
    case PrecisionModel::SINGLE:
      error_message += "f";
      break;
    }
    error_message += "PmeDualValence";
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
    rtErr(error_message, "queryValenceKernelRequirementsDPMEDual");
  }
  
  return result;
}
#endif

//-------------------------------------------------------------------------------------------------
extern void launchValence(const SyValenceKit<double> &poly_vk,
                          const SyRestraintKit<double, double2, double4_16a> &poly_rk,
                          const CellGridReader<double, llint, double, double4_16a> &cgr_qq,
                          const CellGridReader<double, llint, double, double4_16a> &cgr_lj,
                          MMControlKit<double> *ctrl, PsSynthesisWriter *poly_psw,
                          const SyAtomUpdateKit<double, double2, double4_16a> &poly_auk,
                          ThermostatWriter<double> *tstw, ScoreCardWriter *scw,
                          CacheResourceKit<double> *gmem_r, const EvaluateForce eval_force,
                          const EvaluateEnergy eval_energy, const VwuGoal purpose, const int2 bt,
                          const double clash_distance, const double clash_ratio) {
  switch (purpose) {
  case VwuGoal::ACCUMULATE:
    launchValence(poly_vk, poly_rk, ctrl, poly_psw, poly_auk, tstw, scw, gmem_r, eval_force,
                  eval_energy, purpose, bt, clash_distance, clash_ratio);
    break;
  case VwuGoal::MOVE_PARTICLES:
    if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        kdsdPmeDualValenceEnergyAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                   *poly_psw, cgr_qq, cgr_lj,
                                                                   clash_distance, clash_ratio,
                                                                   poly_auk, *tstw, *scw, *gmem_r);
        break;
      case EvaluateEnergy::NO:
        kdsdPmeDualValenceAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                             cgr_qq, cgr_lj, clash_distance,
                                                             clash_ratio, poly_auk, *tstw,
                                                             *gmem_r);
        break;
      }
    }
    else {
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        kdsdPmeDualValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           cgr_qq, cgr_lj, poly_auk, *tstw, *scw,
                                                           *gmem_r);
        break;
      case EvaluateEnergy::NO:
        kdsdPmeDualValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr_qq,
                                                     cgr_lj, poly_auk, *tstw, *gmem_r);
        break;
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchValence(const SyValenceKit<double> &poly_vk,
                          const SyRestraintKit<double, double2, double4_16a> &poly_rk,
                          const CellGridReader<float, int, float, float4> &cgr_qq,
                          const CellGridReader<float, int, float, float4> &cgr_lj,
                          MMControlKit<double> *ctrl, PsSynthesisWriter *poly_psw,
                          const SyAtomUpdateKit<double, double2, double4_16a> &poly_auk,
                          ThermostatWriter<double> *tstw, ScoreCardWriter *scw,
                          CacheResourceKit<double> *gmem_r, const EvaluateForce eval_force,
                          const EvaluateEnergy eval_energy, const VwuGoal purpose, const int2 bt,
                          const double clash_distance, const double clash_ratio) {
  switch (purpose) {
  case VwuGoal::ACCUMULATE:
    launchValence(poly_vk, poly_rk, ctrl, poly_psw, poly_auk, tstw, scw, gmem_r, eval_force,
                  eval_energy, purpose, bt, clash_distance, clash_ratio);
    break;
  case VwuGoal::MOVE_PARTICLES:
    if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        kdsfPmeDualValenceEnergyAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                   *poly_psw, cgr_qq, cgr_lj,
                                                                   clash_distance, clash_ratio,
                                                                   poly_auk, *tstw, *scw, *gmem_r);
        break;
      case EvaluateEnergy::NO:
        kdsfPmeDualValenceAtomUpdateNonClash<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                             cgr_qq, cgr_lj, clash_distance,
                                                             clash_ratio, poly_auk, *tstw,
                                                             *gmem_r);
        break;
      }
    }
    else {
      switch (eval_energy) {
      case EvaluateEnergy::YES:
        kdsfPmeDualValenceEnergyAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           cgr_qq, cgr_lj, poly_auk, *tstw, *scw,
                                                           *gmem_r);
        break;
      case EvaluateEnergy::NO:
        kdsfPmeDualValenceAtomUpdate<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, cgr_qq,
                                                     cgr_lj, poly_auk, *tstw, *gmem_r);
        break;
      }
    }
    break;
  }
}

} // namespace energy
} // namespace stormm
