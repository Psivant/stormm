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
#include "Potential/energy_enumerators.h"
#include "Synthesis/valence_workunit.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/trajectory_enumerators.h"
#include "hpc_valence_potential.h"

namespace stormm {
namespace energy {

using card::GpuDetails;
using card::CoreKlManager;
using constants::PrecisionModel;
using constants::large_block_size;
using constants::medium_block_size;
using constants::small_block_size;
using constants::twice_warp_bits_mask_int;
using constants::twice_warp_size_int;
using constants::warp_size_int;
using constants::warp_bits;
using constants::warp_bits_mask_int;
using numerics::chooseAccumulationMethod;
using numerics::AccumulationMethod;
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
using synthesis::VwuGoal;
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

//-------------------------------------------------------------------------------------------------
// Compute an angle based on the value of its cosine, with the understanding that a fallback
// method is appropriate when the angle is too actue for acos to be numerically well-conditioned.
//
// Overloaded:
//   - Single-precision version
//   - Double-precision version
//
// Arguments:
//   costheta:   Cosine value of the angle of interest
//   crabbc:     The first of two vectors decsribing the displacements that determine the angle
//   crbccd:     The second of two vectors decsribing the displacements that determine the angle
//   bc:         Vector defining the directionality of the angle
//   scr:        Second vector defining the directionality of the angle
//-------------------------------------------------------------------------------------------------
__device__ __forceinline__ float devcAngleVerification(const float costheta, const float3 crabbc,
                                                       const float3 crbccd, const float3 bc,
                                                       const float3 scr) {
  if (fabsf(costheta) >= near_to_one_f) {

    // The floating-point representation of costheta is numerically ill-conditioned.  Compute the
    // distance from atom I to the plane of atoms J, K, and L to get the angle by the arcsin of an
    // extremely acute angle.
    const float mg_crabbc = 1.0f / sqrtf((crabbc.x * crabbc.x) + (crabbc.y * crabbc.y) +
                                         (crabbc.z * crabbc.z));
    const float mg_crbccd = 1.0f / sqrtf((crbccd.x * crbccd.x) + (crbccd.y * crbccd.y) +
                                         (crbccd.z * crbccd.z));
    const float nx_abbc = crabbc.x * mg_crabbc;
    const float ny_abbc = crabbc.y * mg_crabbc;
    const float nz_abbc = crabbc.z * mg_crabbc;
    const float nx_bccd = crbccd.x * mg_crbccd;
    const float ny_bccd = crbccd.y * mg_crbccd;
    const float nz_bccd = crbccd.z * mg_crbccd;
    float rdx = nx_bccd - nx_abbc;
    float rdy = ny_bccd - ny_abbc;
    float rdz = nz_bccd - nz_abbc;
    float rs = sqrtf((rdx * rdx) + (rdy * rdy) + (rdz * rdz));
    if (fabsf(rs) > 1.0f) {
      rdx = nx_bccd + nx_abbc;
      rdy = ny_bccd + ny_abbc;
      rdz = nz_bccd + nz_abbc;
      rs = pi_f - sqrtf((rdx * rdx) + (rdy * rdy) + (rdz * rdz));
    }
    return ((scr.x * bc.x) + (scr.y * bc.y) + (scr.z * bc.z) > 0.0f) ? rs : -rs;
  }
  else {
    return ((scr.x * bc.x) + (scr.y * bc.y) + (scr.z * bc.z) > 0.0f) ?
            acosf(costheta) : -acosf(costheta);
  }
}

__device__ __forceinline__ double devcAngleVerification(const double costheta,
                                                        const double3 crabbc, const double3 crbccd,
                                                        const double3 bc, const double3 scr) {
  if (fabs(costheta) >= near_to_one_lf) {

    // The double-precision arccosine function is also vulnerable to numerical instability near
    // zero, so planar dihedral angles can still generate divergent forces on the order of 3.0e-7
    // kcal/mol-A.  Correct this with a similar strategy to the single-precision case.
    const double mg_crabbc = 1.0 / sqrt((crabbc.x * crabbc.x) + (crabbc.y * crabbc.y) +
                                        (crabbc.z * crabbc.z));
    const double mg_crbccd = 1.0 / sqrt((crbccd.x * crbccd.x) + (crbccd.y * crbccd.y) +
                                        (crbccd.z * crbccd.z));
    const double nx_abbc = crabbc.x * mg_crabbc;
    const double ny_abbc = crabbc.y * mg_crabbc;
    const double nz_abbc = crabbc.z * mg_crabbc;
    const double nx_bccd = crbccd.x * mg_crbccd;
    const double ny_bccd = crbccd.y * mg_crbccd;
    const double nz_bccd = crbccd.z * mg_crbccd;
    double rdx = nx_bccd - nx_abbc;
    double rdy = ny_bccd - ny_abbc;
    double rdz = nz_bccd - nz_abbc;
    double rs = sqrt((rdx * rdx) + (rdy * rdy) + (rdz * rdz));
    if (fabs(rs) > 1.0) {
      rdx = nx_bccd + nx_abbc;
      rdy = ny_bccd + ny_abbc;
      rdz = nz_bccd + nz_abbc;
      rs = pi - sqrt((rdx * rdx) + (rdy * rdy) + (rdz * rdz));
    }
    return ((scr.x * bc.x) + (scr.y * bc.y) + (scr.z * bc.z) > 0.0) ? rs : -rs;
  }
  else {
    return ((scr.x * bc.x) + (scr.y * bc.y) + (scr.z * bc.z) > 0.0) ?
            acos(costheta) : -acos(costheta);
  }
}

//-------------------------------------------------------------------------------------------------
// Compute critical elements of the restraining potential: its difference from the target value
// that determines some harmonic stiffness penalty, the harmonic penalty stiffness, and the energy
// contribution.
//
// Overloaded:
//   - Single-precision version
//   - Double-precision version
//
// Arguments:
//   init_k   Initial stiffness parameters
//   final_k  Final stiffness parameters
//   init_r   Initial displacement parameters
//   final_r  Final displacement parameters
//   mixwt    Pre-calculated mixing factor for combining initial and final parameters
//   dr       The measured value of the restraint coordinate among its participating atoms
//-------------------------------------------------------------------------------------------------
__device__ __forceinline__
double3 restraintDelta(const double2 init_k, const double2 final_k, const double4 init_r,
                      const double4 final_r, const double2 mixwt, const double dr) {
  const double r1 = (mixwt.x * init_r.x) + (mixwt.y * final_r.x);
  const double r2 = (mixwt.x * init_r.y) + (mixwt.y * final_r.y);
  const double r3 = (mixwt.x * init_r.z) + (mixwt.y * final_r.z);
  const double r4 = (mixwt.x * init_r.w) + (mixwt.y * final_r.w);
  const double k2 = (mixwt.x * init_k.x) + (mixwt.y * final_k.x);
  const double k3 = (mixwt.x * init_k.y) + (mixwt.y * final_k.y);
  double dl, du, keq;
  if (dr < r1) {
    dl = r1 - r2;
    du = k2 * ((dl * dl) + (2.0 * dl * (dr - r1)));
    keq = k2;
  }
  else if (dr < r2) {
    dl = dr - r2;
    du = k2 * dl * dl;
    keq = k2;
  }
  else if (dr < r3) {
    dl = 0.0;
    du = 0.0;
    keq = 0.0;
  }
  else if (dr < r4) {
    dl = dr - r3;
    du = k3 * dl * dl;
    keq = k3;
  }
  else {
    dl = r4 - r3;
    du = k3 * ((dl * dl) + (2.0 * dl * (dr - r4)));
    keq = k3;
  }
  return { keq, dl, du };
}

__device__ __forceinline__
float3 restraintDelta(const float2 init_k, const float2 final_k, const float4 init_r,
                      const float4 final_r, const float2 mixwt, const float dr) {
  const float r1 = (mixwt.x * init_r.x) + (mixwt.y * final_r.x);
  const float r2 = (mixwt.x * init_r.y) + (mixwt.y * final_r.y);
  const float r3 = (mixwt.x * init_r.z) + (mixwt.y * final_r.z);
  const float r4 = (mixwt.x * init_r.w) + (mixwt.y * final_r.w);
  const float k2 = (mixwt.x * init_k.x) + (mixwt.y * final_k.x);
  const float k3 = (mixwt.x * init_k.y) + (mixwt.y * final_k.y);
  float dl, du, keq;
  if (dr < r1) {
    dl = r1 - r2;
    du = k2 * ((dl * dl) + (2.0 * dl * (dr - r1)));
    keq = k2;
  }
  else if (dr < r2) {
    dl = dr - r2;
    du = k2 * dl * dl;
    keq = k2;
  }
  else if (dr < r3) {
    dl = 0.0;
    du = 0.0;
    keq = 0.0;
  }
  else if (dr < r4) {
    dl = dr - r3;
    du = k3 * dl * dl;
    keq = k3;
  }
  else {
    dl = r4 - r3;
    du = k3 * ((dl * dl) + (2.0 * dl * (dr - r4)));
    keq = k3;
  }
  return { keq, dl, du };
}

//-------------------------------------------------------------------------------------------------
// Compute the mixture of end-point values that will determine the actual strength and displacement
// settings of a flat-bottom bimodal harmonic restraint.  The flag about a RestraintApparatus
// having time-dependent restraints is mostly for convenience, a way to tell whether there is any
// time-dependent restraint in the collection at all.  Initial and final settings of the steps for
// each restraint encode whether there is actual time dependence in the result.
//
// Overloaded:
//   - Single-precision version
//   - Double-precision version
//
// Arguments:
//   step_number  The current step number of the simulation (may include energy minimization step
//                counts)
//   init_step    The initial step at which the restraint engages
//   final_step   The final step at which the restraint becomes mature
//-------------------------------------------------------------------------------------------------
__device__ __forceinline__
double2 computeRestraintMixtureD(const int step_number, const int init_step,
                                 const int final_step) {
  if (step_number < init_step) {

    // If the restraint has not yet engaged, neither its initial or final values have any weight
    return { (double)(0.0), (double)(0.0) };
  }
  else if (init_step == final_step) {

    // The step count is far enough along that the restraint has been engaged, and it is constant.
    // Only the initial value matters.
    return { (double)(1.0), (double)(0.0) };
  }
  else if (step_number < final_step) {
    const double wslide = (double)(step_number - init_step) / (double)(final_step - init_step);

    // The difference between the initial and final steps is nonzero.  The mixture is a linear
    // combination of the two end points.
    return { (double)(1.0) - wslide, wslide };
  }

  // The step number has advanced beyond the point at which the restraint is mature.
  return { (double)(0.0), (double)(1.0) };
}

__device__ __forceinline__
float2 computeRestraintMixtureF(const int step_number, const int init_step, const int final_step) {
  if (step_number < init_step) {

    // If the restraint has not yet engaged, neither its initial or final values have any weight
    return { (float)(0.0), (float)(0.0) };
  }
  else if (init_step == final_step) {

    // The step count is far enough along that the restraint has been engaged, and it is constant.
    // Only the initial value matters.
    return { (float)(1.0), (float)(0.0) };
  }
  else if (step_number < final_step) {
    const float wslide = (float)(step_number - init_step) / (float)(final_step - init_step);

    // The difference between the initial and final steps is nonzero.  The mixture is a linear
    // combination of the two end points.
    return { (float)(1.0) - wslide, wslide };
  }

  // The step number has advanced beyond the point at which the restraint is mature.
  return { (float)(0.0), (float)(1.0) };
}

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

// Compile the standard kernels with all combinations of energy, and force accumulation methods.
#  define COMPUTE_FORCE
#    define SPLIT_FORCE_ACCUMULATION
#      define VALENCE_KERNEL_THREAD_COUNT 512
#        define VALENCE_BLOCK_MULTIPLICITY 2
#        define KERNEL_NAME kfsValenceForceAccumulationXL
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 256
#        define VALENCE_BLOCK_MULTIPLICITY 4
#        define KERNEL_NAME kfsValenceForceAccumulationLG
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 128
#        define VALENCE_BLOCK_MULTIPLICITY 8
#        define KERNEL_NAME kfsValenceForceAccumulationMD
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 64
#        define VALENCE_BLOCK_MULTIPLICITY 16
#        define KERNEL_NAME kfsValenceForceAccumulationSM
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define UPDATE_ATOMS
#        define VALENCE_KERNEL_THREAD_COUNT 512
#          define VALENCE_BLOCK_MULTIPLICITY 2
#          define KERNEL_NAME kfsValenceAtomUpdateXL
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 256
#          define VALENCE_BLOCK_MULTIPLICITY 4
#          define KERNEL_NAME kfsValenceAtomUpdateLG
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 128
#          define VALENCE_BLOCK_MULTIPLICITY 8
#          define KERNEL_NAME kfsValenceAtomUpdateMD
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 64
#          define VALENCE_BLOCK_MULTIPLICITY 16
#          define KERNEL_NAME kfsValenceAtomUpdateSM
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef UPDATE_ATOMS
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT 448
#          define VALENCE_BLOCK_MULTIPLICITY 2
#          define KERNEL_NAME kfsValenceForceEnergyAccumulationXL
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 224
#          define VALENCE_BLOCK_MULTIPLICITY 4
#          define KERNEL_NAME kfsValenceForceEnergyAccumulationLG
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 128
#          define VALENCE_BLOCK_MULTIPLICITY 7
#          define KERNEL_NAME kfsValenceForceEnergyAccumulationMD
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 64
#          define VALENCE_BLOCK_MULTIPLICITY 14
#          define KERNEL_NAME kfsValenceForceEnergyAccumulationSM
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define UPDATE_ATOMS
#          define VALENCE_KERNEL_THREAD_COUNT 384
#            define VALENCE_BLOCK_MULTIPLICITY 2
#            define KERNEL_NAME kfsValenceEnergyAtomUpdateXL
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#          define VALENCE_KERNEL_THREAD_COUNT 192
#            define VALENCE_BLOCK_MULTIPLICITY 4
#            define KERNEL_NAME kfsValenceEnergyAtomUpdateLG
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#          define VALENCE_KERNEL_THREAD_COUNT 96
#            define VALENCE_BLOCK_MULTIPLICITY 8
#            define KERNEL_NAME kfsValenceEnergyAtomUpdateMD
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#          define VALENCE_KERNEL_THREAD_COUNT 64
#            define VALENCE_BLOCK_MULTIPLICITY 12
#            define KERNEL_NAME kfsValenceEnergyAtomUpdateSM
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#        undef UPDATE_ATOMS
#      undef COMPUTE_ENERGY
#      undef VALENCE_KERNEL_THREAD_COUNT
#    undef SPLIT_FORCE_ACCUMULATION
#    define VALENCE_KERNEL_THREAD_COUNT 512
#      define VALENCE_BLOCK_MULTIPLICITY 2
#      define KERNEL_NAME kfValenceForceAccumulationXL
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#      undef VALENCE_BLOCK_MULTIPLICITY
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define VALENCE_KERNEL_THREAD_COUNT 256
#      define VALENCE_BLOCK_MULTIPLICITY 4
#      define KERNEL_NAME kfValenceForceAccumulationLG
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#      undef VALENCE_BLOCK_MULTIPLICITY
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define VALENCE_KERNEL_THREAD_COUNT 128
#      define VALENCE_BLOCK_MULTIPLICITY 8
#      define KERNEL_NAME kfValenceForceAccumulationMD
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#      undef VALENCE_BLOCK_MULTIPLICITY
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define VALENCE_KERNEL_THREAD_COUNT 64
#      define VALENCE_BLOCK_MULTIPLICITY 16
#      define KERNEL_NAME kfValenceForceAccumulationSM
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#      undef VALENCE_BLOCK_MULTIPLICITY
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define UPDATE_ATOMS
#      define VALENCE_KERNEL_THREAD_COUNT 448
#        define VALENCE_BLOCK_MULTIPLICITY 2
#        define KERNEL_NAME kfValenceAtomUpdateXL
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 224
#        define VALENCE_BLOCK_MULTIPLICITY 4
#        define KERNEL_NAME kfValenceAtomUpdateLG
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 128
#        define VALENCE_BLOCK_MULTIPLICITY 7
#        define KERNEL_NAME kfValenceAtomUpdateMD
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 64
#        define VALENCE_BLOCK_MULTIPLICITY 14
#        define KERNEL_NAME kfValenceAtomUpdateSM
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#    undef UPDATE_ATOMS
#    define COMPUTE_ENERGY
#      define VALENCE_KERNEL_THREAD_COUNT 448
#        define VALENCE_BLOCK_MULTIPLICITY 2
#        define KERNEL_NAME kfValenceForceEnergyAccumulationXL
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 224
#        define VALENCE_BLOCK_MULTIPLICITY 4
#        define KERNEL_NAME kfValenceForceEnergyAccumulationLG
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 128
#        define VALENCE_BLOCK_MULTIPLICITY 7
#        define KERNEL_NAME kfValenceForceEnergyAccumulationMD
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 64
#        define VALENCE_BLOCK_MULTIPLICITY 14
#        define KERNEL_NAME kfValenceForceEnergyAccumulationSM
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define UPDATE_ATOMS
#        define VALENCE_KERNEL_THREAD_COUNT 384
#          define VALENCE_BLOCK_MULTIPLICITY 2
#          define KERNEL_NAME kfValenceEnergyAtomUpdateXL
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 192
#          define VALENCE_BLOCK_MULTIPLICITY 4
#          define KERNEL_NAME kfValenceEnergyAtomUpdateLG
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 96
#          define VALENCE_BLOCK_MULTIPLICITY 8
#          define KERNEL_NAME kfValenceEnergyAtomUpdateMD
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 64
#          define VALENCE_BLOCK_MULTIPLICITY 12
#          define KERNEL_NAME kfValenceEnergyAtomUpdateSM
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef UPDATE_ATOMS
#    undef COMPUTE_ENERGY
#    undef VALENCE_KERNEL_THREAD_COUNT
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define VALENCE_KERNEL_THREAD_COUNT 512
#      define VALENCE_BLOCK_MULTIPLICITY 2
#      define KERNEL_NAME kfValenceEnergyAccumulationXL
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#      undef VALENCE_BLOCK_MULTIPLICITY
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define VALENCE_KERNEL_THREAD_COUNT 256
#      define VALENCE_BLOCK_MULTIPLICITY 4
#      define KERNEL_NAME kfValenceEnergyAccumulationLG
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#      undef VALENCE_BLOCK_MULTIPLICITY
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define VALENCE_KERNEL_THREAD_COUNT 128
#      define VALENCE_BLOCK_MULTIPLICITY 8
#      define KERNEL_NAME kfValenceEnergyAccumulationMD
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#      undef VALENCE_BLOCK_MULTIPLICITY
#    undef VALENCE_KERNEL_THREAD_COUNT
#    define VALENCE_KERNEL_THREAD_COUNT 64
#      define VALENCE_BLOCK_MULTIPLICITY 16
#      define KERNEL_NAME kfValenceEnergyAccumulationSM
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#      undef VALENCE_BLOCK_MULTIPLICITY
#    undef VALENCE_KERNEL_THREAD_COUNT
#  undef COMPUTE_ENERGY

// Make new kernels with a clash forgiveness check.
#  define CLASH_FORGIVENESS
#    define COMPUTE_FORCE
#      define SPLIT_FORCE_ACCUMULATION
#        define VALENCE_KERNEL_THREAD_COUNT 512
#          define VALENCE_BLOCK_MULTIPLICITY 2
#          define KERNEL_NAME kfsValenceForceAccumulationNonClashXL
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 256
#          define VALENCE_BLOCK_MULTIPLICITY 4
#          define KERNEL_NAME kfsValenceForceAccumulationNonClashLG
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 128
#          define VALENCE_BLOCK_MULTIPLICITY 8
#          define KERNEL_NAME kfsValenceForceAccumulationNonClashMD
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 64
#          define VALENCE_BLOCK_MULTIPLICITY 16
#          define KERNEL_NAME kfsValenceForceAccumulationNonClashSM
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define UPDATE_ATOMS
#          define VALENCE_KERNEL_THREAD_COUNT 448
#            define VALENCE_BLOCK_MULTIPLICITY 2
#            define KERNEL_NAME kfsValenceAtomUpdateNonClashXL
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#          define VALENCE_KERNEL_THREAD_COUNT 224
#            define VALENCE_BLOCK_MULTIPLICITY 4
#            define KERNEL_NAME kfsValenceAtomUpdateNonClashLG
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#          define VALENCE_KERNEL_THREAD_COUNT 128
#            define VALENCE_BLOCK_MULTIPLICITY 7
#            define KERNEL_NAME kfsValenceAtomUpdateNonClashMD
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#          define VALENCE_KERNEL_THREAD_COUNT 64
#            define VALENCE_BLOCK_MULTIPLICITY 14
#            define KERNEL_NAME kfsValenceAtomUpdateNonClashSM
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#        undef UPDATE_ATOMS
#        define COMPUTE_ENERGY
#          define VALENCE_KERNEL_THREAD_COUNT 384
#            define VALENCE_BLOCK_MULTIPLICITY 2
#            define KERNEL_NAME kfsValenceForceEnergyAccumulationNonClashXL
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#          define VALENCE_KERNEL_THREAD_COUNT 192
#            define VALENCE_BLOCK_MULTIPLICITY 4
#            define KERNEL_NAME kfsValenceForceEnergyAccumulationNonClashLG
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#          define VALENCE_KERNEL_THREAD_COUNT 96
#            define VALENCE_BLOCK_MULTIPLICITY 8
#            define KERNEL_NAME kfsValenceForceEnergyAccumulationNonClashMD
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#          define VALENCE_KERNEL_THREAD_COUNT 64
#            define VALENCE_BLOCK_MULTIPLICITY 12
#            define KERNEL_NAME kfsValenceForceEnergyAccumulationNonClashSM
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#          define UPDATE_ATOMS
#            define VALENCE_KERNEL_THREAD_COUNT 320
#              define VALENCE_BLOCK_MULTIPLICITY 2
#              define KERNEL_NAME kfsValenceEnergyAtomUpdateNonClashXL
#                include "valence_potential.cui"
#              undef KERNEL_NAME
#              undef VALENCE_BLOCK_MULTIPLICITY
#            undef VALENCE_KERNEL_THREAD_COUNT
#            define VALENCE_KERNEL_THREAD_COUNT 160
#              define VALENCE_BLOCK_MULTIPLICITY 4
#              define KERNEL_NAME kfsValenceEnergyAtomUpdateNonClashLG
#                include "valence_potential.cui"
#              undef KERNEL_NAME
#              undef VALENCE_BLOCK_MULTIPLICITY
#            undef VALENCE_KERNEL_THREAD_COUNT
#            define VALENCE_KERNEL_THREAD_COUNT 128
#              define VALENCE_BLOCK_MULTIPLICITY 5
#              define KERNEL_NAME kfsValenceEnergyAtomUpdateNonClashMD
#                include "valence_potential.cui"
#              undef KERNEL_NAME
#              undef VALENCE_BLOCK_MULTIPLICITY
#            undef VALENCE_KERNEL_THREAD_COUNT
#            define VALENCE_KERNEL_THREAD_COUNT 64
#              define VALENCE_BLOCK_MULTIPLICITY 10
#              define KERNEL_NAME kfsValenceEnergyAtomUpdateNonClashSM
#                include "valence_potential.cui"
#              undef KERNEL_NAME
#              undef VALENCE_BLOCK_MULTIPLICITY
#            undef VALENCE_KERNEL_THREAD_COUNT
#          undef UPDATE_ATOMS
#        undef COMPUTE_ENERGY
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef SPLIT_FORCE_ACCUMULATION
#      define VALENCE_KERNEL_THREAD_COUNT 512
#        define VALENCE_BLOCK_MULTIPLICITY 2
#        define KERNEL_NAME kfValenceForceAccumulationNonClashXL
#          include "valence_potential.cui"
#        undef KERNEL_NAME  
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 256
#        define VALENCE_BLOCK_MULTIPLICITY 4
#        define KERNEL_NAME kfValenceForceAccumulationNonClashLG
#          include "valence_potential.cui"
#        undef KERNEL_NAME  
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 128
#        define VALENCE_BLOCK_MULTIPLICITY 8
#        define KERNEL_NAME kfValenceForceAccumulationNonClashMD
#          include "valence_potential.cui"
#        undef KERNEL_NAME  
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 64
#        define VALENCE_BLOCK_MULTIPLICITY 16
#        define KERNEL_NAME kfValenceForceAccumulationNonClashSM
#          include "valence_potential.cui"
#        undef KERNEL_NAME  
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define UPDATE_ATOMS
#        define VALENCE_KERNEL_THREAD_COUNT 448
#          define VALENCE_BLOCK_MULTIPLICITY 2
#          define KERNEL_NAME kfValenceAtomUpdateNonClashXL
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 224
#          define VALENCE_BLOCK_MULTIPLICITY 4
#          define KERNEL_NAME kfValenceAtomUpdateNonClashLG
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 128
#          define VALENCE_BLOCK_MULTIPLICITY 7
#          define KERNEL_NAME kfValenceAtomUpdateNonClashMD
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 64
#          define VALENCE_BLOCK_MULTIPLICITY 14
#          define KERNEL_NAME kfValenceAtomUpdateNonClashSM
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#      undef UPDATE_ATOMS
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT 384
#          define VALENCE_BLOCK_MULTIPLICITY 2
#          define KERNEL_NAME kfValenceForceEnergyAccumulationNonClashXL
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 192
#          define VALENCE_BLOCK_MULTIPLICITY 4
#          define KERNEL_NAME kfValenceForceEnergyAccumulationNonClashLG
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 96
#          define VALENCE_BLOCK_MULTIPLICITY 8
#          define KERNEL_NAME kfValenceForceEnergyAccumulationNonClashMD
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define VALENCE_KERNEL_THREAD_COUNT 64
#          define VALENCE_BLOCK_MULTIPLICITY 12
#          define KERNEL_NAME kfValenceForceEnergyAccumulationNonClashSM
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_BLOCK_MULTIPLICITY
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define UPDATE_ATOMS
#          define VALENCE_KERNEL_THREAD_COUNT 320
#            define VALENCE_BLOCK_MULTIPLICITY 2
#            define KERNEL_NAME kfValenceEnergyAtomUpdateNonClashXL
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#          define VALENCE_KERNEL_THREAD_COUNT 160
#            define VALENCE_BLOCK_MULTIPLICITY 4
#            define KERNEL_NAME kfValenceEnergyAtomUpdateNonClashLG
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#          define VALENCE_KERNEL_THREAD_COUNT 128
#            define VALENCE_BLOCK_MULTIPLICITY 5
#            define KERNEL_NAME kfValenceEnergyAtomUpdateNonClashMD
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#          define VALENCE_KERNEL_THREAD_COUNT 64
#            define VALENCE_BLOCK_MULTIPLICITY 10
#            define KERNEL_NAME kfValenceEnergyAtomUpdateNonClashSM
#              include "valence_potential.cui"
#            undef KERNEL_NAME
#            undef VALENCE_BLOCK_MULTIPLICITY
#          undef VALENCE_KERNEL_THREAD_COUNT
#        undef UPDATE_ATOMS
#      undef COMPUTE_ENERGY
#      undef VALENCE_KERNEL_THREAD_COUNT
#    undef COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define VALENCE_KERNEL_THREAD_COUNT 512
#        define VALENCE_BLOCK_MULTIPLICITY 2
#        define KERNEL_NAME kfValenceEnergyAccumulationNonClashXL
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 256
#        define VALENCE_BLOCK_MULTIPLICITY 4
#        define KERNEL_NAME kfValenceEnergyAccumulationNonClashLG
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 128
#        define VALENCE_BLOCK_MULTIPLICITY 8
#        define KERNEL_NAME kfValenceEnergyAccumulationNonClashMD
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define VALENCE_KERNEL_THREAD_COUNT 64
#        define VALENCE_BLOCK_MULTIPLICITY 16
#        define KERNEL_NAME kfValenceEnergyAccumulationNonClashSM
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_BLOCK_MULTIPLICITY
#      undef VALENCE_KERNEL_THREAD_COUNT
#    undef COMPUTE_ENERGY
#  undef CLASH_FORGIVENESS

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

// Double-precision floating point definitions
#define TCALC double
#  define VALENCE_BLOCK_MULTIPLICITY  2
#  define TCALC2 double2
#  define TCALC3 double3
#  define TCALC4 double4
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
#  define VALENCE_KERNEL_THREAD_COUNT 256
#  define COMPUTE_FORCE
#    define KERNEL_NAME kdsValenceForceAccumulation
#      include "valence_potential.cui"
#    undef KERNEL_NAME  
#    define UPDATE_ATOMS
#      define KERNEL_NAME kdsValenceAtomUpdate
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#    undef UPDATE_ATOMS
#    define COMPUTE_ENERGY
#      define KERNEL_NAME kdsValenceForceEnergyAccumulation
#        include "valence_potential.cui"
#      undef KERNEL_NAME
#      define UPDATE_ATOMS
#        define KERNEL_NAME kdsValenceEnergyAtomUpdate
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef UPDATE_ATOMS
#    undef  COMPUTE_ENERGY
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define KERNEL_NAME kdsValenceEnergyAccumulation
#      include "valence_potential.cui"
#    undef KERNEL_NAME
#  undef COMPUTE_ENERGY
#  undef VALENCE_KERNEL_THREAD_COUNT

// Make new kernels with a clash forgiveness check.
#  define CLASH_FORGIVENESS
#    define COMPUTE_FORCE
#      define VALENCE_KERNEL_THREAD_COUNT 256
#      define KERNEL_NAME kdsValenceForceAccumulationNonClash
#        include "valence_potential.cui"
#      undef KERNEL_NAME  
#      define UPDATE_ATOMS
#        define KERNEL_NAME kdsValenceAtomUpdateNonClash
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#      undef UPDATE_ATOMS
#      undef VALENCE_KERNEL_THREAD_COUNT
#      define COMPUTE_ENERGY
#        define VALENCE_KERNEL_THREAD_COUNT small_block_size
#        define KERNEL_NAME kdsValenceForceEnergyAccumulationNonClash
#          include "valence_potential.cui"
#        undef KERNEL_NAME
#        undef VALENCE_KERNEL_THREAD_COUNT
#        define UPDATE_ATOMS
#          define VALENCE_KERNEL_THREAD_COUNT 192
#          define KERNEL_NAME kdsValenceEnergyAtomUpdateNonClash
#            include "valence_potential.cui"
#          undef KERNEL_NAME
#          undef VALENCE_KERNEL_THREAD_COUNT
#        undef UPDATE_ATOMS
#      undef  COMPUTE_ENERGY
#    undef COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define VALENCE_KERNEL_THREAD_COUNT small_block_size
#      define KERNEL_NAME kdsValenceEnergyAccumulationNonClash
#        include "valence_potential.cui"
#      undef KERNEL_NAME
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
extern cudaFuncAttributes queryValenceKernelRequirements(const PrecisionModel prec,
                                                         const EvaluateForce eval_frc,
                                                         const EvaluateEnergy eval_nrg,
                                                         const AccumulationMethod acc_meth,
                                                         const VwuGoal purpose,
                                                         const ClashResponse collision_handling,
                                                         const ValenceKernelSize kwidth) {

  // The kernel manager will have information about the GPU to use--look at the work units from
  // the perspective of overall occupancy on the GPU.
  cudaFuncAttributes result;
  cudaError_t cfa;
  switch (collision_handling) {
  case ClashResponse::NONE:
    switch (prec) {
    case PrecisionModel::DOUBLE:
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
      case EvaluateForce::NO:
        cfa = cudaFuncGetAttributes(&result, kdsValenceEnergyAccumulation);
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          switch (acc_meth) {
          case AccumulationMethod::SPLIT:
            switch (purpose) {
            case VwuGoal::ACCUMULATE:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceEnergyAccumulationXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceEnergyAccumulationLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceEnergyAccumulationMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceEnergyAccumulationSM);
                break;
              }
              break;
            case VwuGoal::MOVE_PARTICLES:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfsValenceEnergyAtomUpdateXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfsValenceEnergyAtomUpdateLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfsValenceEnergyAtomUpdateMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfsValenceEnergyAtomUpdateSM);
                break;
              }
              break;
            }
            break;
          case AccumulationMethod::WHOLE:
            switch (purpose) {
            case VwuGoal::ACCUMULATE:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceEnergyAccumulationXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceEnergyAccumulationLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceEnergyAccumulationMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceEnergyAccumulationSM);
                break;
              }
              break;
            case VwuGoal::MOVE_PARTICLES:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAtomUpdateXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAtomUpdateLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAtomUpdateMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAtomUpdateSM);
                break;
              }
              break;
            }
            break;
          }
          break;
        case EvaluateEnergy::NO:
          switch (acc_meth) {
          case AccumulationMethod::SPLIT:
            switch (purpose) {
            case VwuGoal::ACCUMULATE:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceAccumulationXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceAccumulationLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceAccumulationMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceAccumulationSM);
                break;
              }
              break;
            case VwuGoal::MOVE_PARTICLES:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfsValenceAtomUpdateXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfsValenceAtomUpdateLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfsValenceAtomUpdateMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfsValenceAtomUpdateSM);
                break;
              }
              break;
            }
            break;
          case AccumulationMethod::WHOLE:
            switch (purpose) {
            case VwuGoal::ACCUMULATE:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceAccumulationXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceAccumulationLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceAccumulationMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceAccumulationSM);
                break;
              }
              break;
            case VwuGoal::MOVE_PARTICLES:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfValenceAtomUpdateXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfValenceAtomUpdateLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfValenceAtomUpdateMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfValenceAtomUpdateSM);
                break;
              }
              break;
            }
            break;
          }
          break;
        }
        break;
      case EvaluateForce::NO:
        switch (kwidth) {
        case ValenceKernelSize::XL:
          cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAccumulationXL);
          break;
        case ValenceKernelSize::LG:
          cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAccumulationLG);
          break;
        case ValenceKernelSize::MD:
          cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAccumulationMD);
          break;
        case ValenceKernelSize::SM:
          cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAccumulationSM);
          break;
        }
        break;
      }
      break;
    }
    break;
  case ClashResponse::FORGIVE:
    switch (prec) {
    case PrecisionModel::DOUBLE:
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
      case EvaluateForce::NO:
        cfa = cudaFuncGetAttributes(&result, kdsValenceEnergyAccumulationNonClash);
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      switch (eval_frc) {
      case EvaluateForce::YES:
        switch (eval_nrg) {
        case EvaluateEnergy::YES:
          switch (acc_meth) {
          case AccumulationMethod::SPLIT:
            switch (purpose) {
            case VwuGoal::ACCUMULATE:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceEnergyAccumulationNonClashXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceEnergyAccumulationNonClashLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceEnergyAccumulationNonClashMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceEnergyAccumulationNonClashSM);
                break;
              }
              break;
            case VwuGoal::MOVE_PARTICLES:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfsValenceEnergyAtomUpdateNonClashXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfsValenceEnergyAtomUpdateNonClashLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfsValenceEnergyAtomUpdateNonClashMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfsValenceEnergyAtomUpdateNonClashSM);
                break;
              }
              break;
            }
            break;
          case AccumulationMethod::WHOLE:
            switch (purpose) {
            case VwuGoal::ACCUMULATE:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceEnergyAccumulationNonClashXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceEnergyAccumulationNonClashLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceEnergyAccumulationNonClashMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceEnergyAccumulationNonClashSM);
                break;
              }
              break;
            case VwuGoal::MOVE_PARTICLES:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAtomUpdateNonClashXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAtomUpdateNonClashLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAtomUpdateNonClashMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAtomUpdateNonClashSM);
                break;
              }
              break;
            }
            break;
          }
          break;
        case EvaluateEnergy::NO:
          switch (acc_meth) {
          case AccumulationMethod::SPLIT:
            switch (purpose) {
            case VwuGoal::ACCUMULATE:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceAccumulationNonClashXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceAccumulationNonClashLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceAccumulationNonClashMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfsValenceForceAccumulationNonClashSM);
                break;
              }
              break;
            case VwuGoal::MOVE_PARTICLES:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfsValenceAtomUpdateNonClashXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfsValenceAtomUpdateNonClashLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfsValenceAtomUpdateNonClashMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfsValenceAtomUpdateNonClashSM);
                break;
              }
              break;
            }
            break;
          case AccumulationMethod::WHOLE:
            switch (purpose) {
            case VwuGoal::ACCUMULATE:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceAccumulationNonClashXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceAccumulationNonClashLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceAccumulationNonClashMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfValenceForceAccumulationNonClashSM);
                break;
              }
              break;
            case VwuGoal::MOVE_PARTICLES:
              switch (kwidth) {
              case ValenceKernelSize::XL:
                cfa = cudaFuncGetAttributes(&result, kfValenceAtomUpdateNonClashXL);
                break;
              case ValenceKernelSize::LG:
                cfa = cudaFuncGetAttributes(&result, kfValenceAtomUpdateNonClashLG);
                break;
              case ValenceKernelSize::MD:
                cfa = cudaFuncGetAttributes(&result, kfValenceAtomUpdateNonClashMD);
                break;
              case ValenceKernelSize::SM:
                cfa = cudaFuncGetAttributes(&result, kfValenceAtomUpdateNonClashSM);
                break;
              }
              break;
            }
            break;
          }
          break;
        }
        break;
      case EvaluateForce::NO:
        switch (kwidth) {
        case ValenceKernelSize::XL:
          cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAccumulationNonClashXL);
          break;
        case ValenceKernelSize::LG:
          cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAccumulationNonClashLG);
          break;
        case ValenceKernelSize::MD:
          cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAccumulationNonClashMD);
          break;
        case ValenceKernelSize::SM:
          cfa = cudaFuncGetAttributes(&result, kfValenceEnergyAccumulationNonClashSM);
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
    std::string error_message("Error obtaining attributes for kernel k");
    switch (prec) {
    case PrecisionModel::DOUBLE:
      error_message += "ds";
      break;
    case PrecisionModel::SINGLE:
      error_message += "f";
      switch (acc_meth) {
      case AccumulationMethod::SPLIT:
        error_message += "s";
        break;
      case AccumulationMethod::WHOLE:
        break;
      case AccumulationMethod::AUTOMATIC:
        rtErr("Kernels do not accept " + getEnumerationName(acc_meth) + " accumulation.",
              "queryValenceKernelRequirements");
        break;
      }
      break;
    }
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
    switch (kwidth) {
    case ValenceKernelSize::XL:
      error_message += "XL";
      break;
    case ValenceKernelSize::LG:
      error_message += "LG";
      break;
    case ValenceKernelSize::MD:
      error_message += "MD";
      break;
    case ValenceKernelSize::SM:
      error_message += "SM";
      break;
    }
    error_message += ".";

    // Report the error
    rtErr(error_message, "queryValenceKernelRequirements");
  }
  return result;
}
#endif

//-------------------------------------------------------------------------------------------------
extern void launchValence(const SyValenceKit<double> &poly_vk,
                          const SyRestraintKit<double, double2, double4> &poly_rk,
                          MMControlKit<double> *ctrl, PsSynthesisWriter *poly_psw,
                          const SyAtomUpdateKit<double, double2, double4> &poly_auk,
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

//-------------------------------------------------------------------------------------------------
extern void launchValence(const SyValenceKit<float> &poly_vk,
                          const SyRestraintKit<float, float2, float4> &poly_rk,
                          MMControlKit<float> *ctrl, PsSynthesisWriter *poly_psw,
                          const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                          ThermostatWriter<float> *tstw, ScoreCardWriter *scw,
                          CacheResourceKit<float> *gmem_r, const EvaluateForce eval_force,
                          const EvaluateEnergy eval_energy, const VwuGoal purpose,
                          const AccumulationMethod force_sum, const int2 bt,
                          const float clash_distance, const float clash_ratio) {
  AccumulationMethod refined_force_sum;
  switch (force_sum) {
  case AccumulationMethod::SPLIT:
  case AccumulationMethod::WHOLE:
    refined_force_sum = force_sum;
    break;
  case AccumulationMethod::AUTOMATIC:
    refined_force_sum = chooseAccumulationMethod(poly_psw->frc_bits);
    break;
  }
  
  // Rather than a switch over cases of the ClashResponse enumerator, just use the nonzero values
  // of either parameter to indicate that clash damping has been requested.
  if (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (purpose) {
    case VwuGoal::ACCUMULATE:
    
      // When the goal is to accumulate energies, forces, or both, the force accumulation method
      // becomes a critical detail when choosing the kernel.
      switch (eval_force) {
      case EvaluateForce::YES:
        switch (refined_force_sum) {
        case AccumulationMethod::SPLIT:
          switch (eval_energy) {
          case EvaluateEnergy::YES:

            // Use the launch grid to determine what size the kernel is.  This requires that all
            // "XL" kernels have > 256 threads per block, "LG" kernels have > 128 threads per
            // block, and "MD" kernels have > 64 threads per block.  All "SM" kernels launch with
            // 64 threads per block.
            if (bt.y > 256) {
              kfsValenceForceEnergyAccumulationNonClashXL<<<bt.x,
                                                            bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                    *poly_psw, clash_distance,
                                                                    clash_ratio, *scw, *gmem_r);
            }
            else if (bt.y > 128) {
              kfsValenceForceEnergyAccumulationNonClashLG<<<bt.x,
                                                            bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                    *poly_psw, clash_distance,
                                                                    clash_ratio, *scw, *gmem_r);
            }
            else if (bt.y > 64) {
              kfsValenceForceEnergyAccumulationNonClashMD<<<bt.x,
                                                            bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                    *poly_psw, clash_distance,
                                                                    clash_ratio, *scw, *gmem_r);
            }
            else {
              kfsValenceForceEnergyAccumulationNonClashSM<<<bt.x,
                                                            bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                    *poly_psw, clash_distance,
                                                                    clash_ratio, *scw, *gmem_r);
            }
            break;
          case EvaluateEnergy::NO:
            if (bt.y > 256) {
              kfsValenceForceAccumulationNonClashXL<<<bt.x,
                                                      bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                              clash_distance, clash_ratio,
                                                              *gmem_r);
            }
            else if (bt.y > 128) {
              kfsValenceForceAccumulationNonClashLG<<<bt.x,
                                                      bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                              clash_distance, clash_ratio,
                                                              *gmem_r);
            }
            else if (bt.y > 64) {
              kfsValenceForceAccumulationNonClashMD<<<bt.x,
                                                      bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                              clash_distance, clash_ratio,
                                                              *gmem_r);
            }
            else {
              kfsValenceForceAccumulationNonClashSM<<<bt.x,
                                                      bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                              clash_distance, clash_ratio,
                                                              *gmem_r);
            }
            break;
          }
          break;
        case AccumulationMethod::WHOLE:
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            if (bt.y > 256) {
              kfValenceForceEnergyAccumulationNonClashXL<<<bt.x,
                                                           bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                   *poly_psw, clash_distance,
                                                                   clash_ratio, *scw, *gmem_r);
            }
            else if (bt.y > 128) {
              kfValenceForceEnergyAccumulationNonClashLG<<<bt.x,
                                                           bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                   *poly_psw, clash_distance,
                                                                   clash_ratio, *scw, *gmem_r);
            }
            else if (bt.y > 64) {
              kfValenceForceEnergyAccumulationNonClashMD<<<bt.x,
                                                           bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                   *poly_psw, clash_distance,
                                                                   clash_ratio, *scw, *gmem_r);
            }
            else {
              kfValenceForceEnergyAccumulationNonClashSM<<<bt.x,
                                                           bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                   *poly_psw, clash_distance,
                                                                   clash_ratio, *scw, *gmem_r);
            }
            break;
          case EvaluateEnergy::NO:
            if (bt.y > 256) {
              kfValenceForceAccumulationNonClashXL<<<bt.x,
                                                     bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                             clash_distance, clash_ratio, *gmem_r);
            }
            else if (bt.y > 128) {
              kfValenceForceAccumulationNonClashLG<<<bt.x,
                                                     bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                             clash_distance, clash_ratio, *gmem_r);
            }
            else if (bt.y > 64) {
              kfValenceForceAccumulationNonClashMD<<<bt.x,
                                                     bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                             clash_distance, clash_ratio, *gmem_r);
            }
            else {
              kfValenceForceAccumulationNonClashSM<<<bt.x,
                                                     bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                             clash_distance, clash_ratio, *gmem_r);
            }
            break;
          }
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateForce::NO:
        if (bt.y > 256) {
          kfValenceEnergyAccumulationNonClashXL<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                                clash_distance, clash_ratio, *scw,
                                                                *gmem_r);
        }
        else if (bt.y > 128) {
          kfValenceEnergyAccumulationNonClashLG<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                                clash_distance, clash_ratio, *scw,
                                                                *gmem_r);
        }
        else if (bt.y > 64) {
          kfValenceEnergyAccumulationNonClashMD<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                                clash_distance, clash_ratio, *scw,
                                                                *gmem_r);
        }
        else {
          kfValenceEnergyAccumulationNonClashSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                                clash_distance, clash_ratio, *scw,
                                                                *gmem_r);
        }
        break;
      }
      break;
    case VwuGoal::MOVE_PARTICLES:
    
      // When the goal is to move particles, evaluating the force is obligatory, but the manner in
      // which forces are accumulated is still important.  Whether to accumulate energies while
      // evaluating forces and moving the particles remains a consideration in choosing the proper
      // kernel.
      switch (refined_force_sum) {
      case AccumulationMethod::SPLIT:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          if (bt.y > 256) {
            kfsValenceEnergyAtomUpdateNonClashXL<<<bt.x,
                                                   bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           clash_distance, clash_ratio, poly_auk,
                                                           *tstw, *scw, *gmem_r);
          }
          else if (bt.y > 128) {
            kfsValenceEnergyAtomUpdateNonClashLG<<<bt.x,
                                                   bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           clash_distance, clash_ratio, poly_auk,
                                                           *tstw, *scw, *gmem_r);
          }
          else if (bt.y > 64) {
            kfsValenceEnergyAtomUpdateNonClashMD<<<bt.x,
                                                   bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           clash_distance, clash_ratio, poly_auk,
                                                           *tstw, *scw, *gmem_r);
          }
          else {
            kfsValenceEnergyAtomUpdateNonClashSM<<<bt.x,
                                                   bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           clash_distance, clash_ratio, poly_auk,
                                                           *tstw, *scw, *gmem_r);
          }
          break;
        case EvaluateEnergy::NO:
          if (bt.y > 256) {
            kfsValenceAtomUpdateNonClashXL<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           clash_distance, clash_ratio, poly_auk,
                                                           *tstw, *gmem_r);
          }
          else if (bt.y > 128) {
            kfsValenceAtomUpdateNonClashLG<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           clash_distance, clash_ratio, poly_auk,
                                                           *tstw, *gmem_r);
          }
          else if (bt.y > 64) {
            kfsValenceAtomUpdateNonClashMD<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           clash_distance, clash_ratio, poly_auk,
                                                           *tstw, *gmem_r);
          }
          else {
            kfsValenceAtomUpdateNonClashSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           clash_distance, clash_ratio, poly_auk,
                                                           *tstw, *gmem_r);
          }
          break;
        }
        break;
      case AccumulationMethod::WHOLE:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          if (bt.y > 256) {
            kfValenceEnergyAtomUpdateNonClashXL<<<bt.x,
                                                  bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                          clash_distance, clash_ratio, poly_auk,
                                                          *tstw, *scw, *gmem_r);
          }
          else if (bt.y > 128) {
            kfValenceEnergyAtomUpdateNonClashLG<<<bt.x,
                                                  bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                          clash_distance, clash_ratio, poly_auk,
                                                          *tstw, *scw, *gmem_r);
          }
          else if (bt.y > 64) {
            kfValenceEnergyAtomUpdateNonClashMD<<<bt.x,
                                                  bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                          clash_distance, clash_ratio, poly_auk,
                                                          *tstw, *scw, *gmem_r);
          }
          else {
            kfValenceEnergyAtomUpdateNonClashSM<<<bt.x,
                                                  bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                          clash_distance, clash_ratio, poly_auk,
                                                          *tstw, *scw, *gmem_r);
          }
          break;
        case EvaluateEnergy::NO:
          if (bt.y > 256) {
            kfValenceAtomUpdateNonClashXL<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                          clash_distance, clash_ratio, poly_auk,
                                                          *tstw,*gmem_r);
          }
          else if (bt.y > 128) {
            kfValenceAtomUpdateNonClashLG<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                          clash_distance, clash_ratio, poly_auk,
                                                          *tstw,*gmem_r);
          }
          else if (bt.y > 64) {
            kfValenceAtomUpdateNonClashMD<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                          clash_distance, clash_ratio, poly_auk,
                                                          *tstw,*gmem_r);
          }
          else {
            kfValenceAtomUpdateNonClashSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                          clash_distance, clash_ratio, poly_auk,
                                                          *tstw,*gmem_r);
          }
          break;
        }
        break;
      case AccumulationMethod::AUTOMATIC:
        break;
      }
      break;
    }
  }
  else {
    switch (purpose) {
    case VwuGoal::ACCUMULATE:

      // See above for the rationale on whether forces or energies are evaluated in each context.
      switch (eval_force) {
      case EvaluateForce::YES:
        switch (refined_force_sum) {
        case AccumulationMethod::SPLIT:
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            if (bt.y > 256) {
              kfsValenceForceEnergyAccumulationXL<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                  *poly_psw, *scw, *gmem_r);
            }
            else if (bt.y > 128) {
              kfsValenceForceEnergyAccumulationLG<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                  *poly_psw, *scw, *gmem_r);
            }
            else if (bt.y > 64) {
              kfsValenceForceEnergyAccumulationMD<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                  *poly_psw, *scw, *gmem_r);
            }
            else {
              kfsValenceForceEnergyAccumulationSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                  *poly_psw, *scw, *gmem_r);
            }
            break;
          case EvaluateEnergy::NO:
            if (bt.y > 256) {
              kfsValenceForceAccumulationXL<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                            *gmem_r);
            }
            else if (bt.y > 128) {
              kfsValenceForceAccumulationLG<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                            *gmem_r);
            }
            else if (bt.y > 64) {
              kfsValenceForceAccumulationMD<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                            *gmem_r);
            }
            else {
              kfsValenceForceAccumulationSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                            *gmem_r);
            }
            break;
          }
          break;
        case AccumulationMethod::WHOLE:
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            if (bt.y > 256) {
              kfValenceForceEnergyAccumulationXL<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                 *poly_psw, *scw, *gmem_r);
            }
            else if (bt.y > 128) {
              kfValenceForceEnergyAccumulationLG<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                 *poly_psw, *scw, *gmem_r);
            }
            else if (bt.y > 64) {
              kfValenceForceEnergyAccumulationMD<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                 *poly_psw, *scw, *gmem_r);
            }
            else {
              kfValenceForceEnergyAccumulationSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl,
                                                                 *poly_psw, *scw, *gmem_r);
            }
            break;
          case EvaluateEnergy::NO:
            if (bt.y > 256) {
              kfValenceForceAccumulationXL<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           *gmem_r);
            }
            else if (bt.y > 128) {
              kfValenceForceAccumulationLG<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           *gmem_r);
            }
            else if (bt.y > 64) {
              kfValenceForceAccumulationMD<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           *gmem_r);
            }
            else {
              kfValenceForceAccumulationSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                           *gmem_r);
            }
            break;
          }
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case EvaluateForce::NO:
        if (bt.y > 256) {
          kfValenceEnergyAccumulationXL<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                        *gmem_r);
        }
        else if (bt.y > 128) {
          kfValenceEnergyAccumulationLG<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                        *gmem_r);
        }
        else if (bt.y > 64) {
          kfValenceEnergyAccumulationMD<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                        *gmem_r);
        }
        else {
          kfValenceEnergyAccumulationSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, *scw,
                                                        *gmem_r);
        }
        break;
      }
      break;
    case VwuGoal::MOVE_PARTICLES:
    
      // See above for the rationale on the choice of each kernel.
      switch (refined_force_sum) {
      case AccumulationMethod::SPLIT:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          if (bt.y > 256) {
            kfsValenceEnergyAtomUpdateXL<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                         poly_auk, *tstw, *scw, *gmem_r);
          }
          else if (bt.y > 128) {
            kfsValenceEnergyAtomUpdateLG<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                         poly_auk, *tstw, *scw, *gmem_r);
          }
          else if (bt.y > 64) {
            kfsValenceEnergyAtomUpdateMD<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                         poly_auk, *tstw, *scw, *gmem_r);
          }
          else {
            kfsValenceEnergyAtomUpdateSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                         poly_auk, *tstw, *scw, *gmem_r);
          }
          break;
        case EvaluateEnergy::NO:
          if (bt.y > 256) {
            kfsValenceAtomUpdateXL<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, poly_auk,
                                                   *tstw, *gmem_r);
          }
          else if (bt.y > 128) {
            kfsValenceAtomUpdateLG<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, poly_auk,
                                                   *tstw, *gmem_r);
          }
          else if (bt.y > 64) {
            kfsValenceAtomUpdateMD<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, poly_auk,
                                                   *tstw, *gmem_r);
          }
          else {
            kfsValenceAtomUpdateSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, poly_auk,
                                                   *tstw, *gmem_r);
          }
          break;
        }
        break;
      case AccumulationMethod::WHOLE:
        switch (eval_energy) {
        case EvaluateEnergy::YES:
          if (bt.y > 256) {
            kfValenceEnergyAtomUpdateXL<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                        poly_auk, *tstw, *scw, *gmem_r);
          }
          else if (bt.y > 128) {
            kfValenceEnergyAtomUpdateLG<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                        poly_auk, *tstw, *scw, *gmem_r);
          }
          else if (bt.y > 64) {
            kfValenceEnergyAtomUpdateMD<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                        poly_auk, *tstw, *scw, *gmem_r);
          }
          else {
            kfValenceEnergyAtomUpdateSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw,
                                                        poly_auk, *tstw, *scw, *gmem_r);
          }
          break;
        case EvaluateEnergy::NO:
          if (bt.y > 256) {
            kfValenceAtomUpdateXL<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, poly_auk,
                                                  *tstw, *gmem_r);
          }
          else if (bt.y > 128) {
            kfValenceAtomUpdateLG<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, poly_auk,
                                                  *tstw, *gmem_r);
          }
          else if (bt.y > 64) {
            kfValenceAtomUpdateMD<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, poly_auk,
                                                  *tstw, *gmem_r);
          }
          else {
            kfValenceAtomUpdateSM<<<bt.x, bt.y>>>(poly_vk, poly_rk, *ctrl, *poly_psw, poly_auk,
                                                  *tstw, *gmem_r);
          }
          break;
        }
        break;
      case AccumulationMethod::AUTOMATIC:
        break;
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchValence(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                          MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                          Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                          const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                          const VwuGoal purpose, const AccumulationMethod force_sum,
                          const CoreKlManager &launcher, const double clash_distance,
                          const double clash_ratio) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  const int2 bt = launcher.getValenceKernelDims(prec, eval_force, eval_energy,
                                                AccumulationMethod::SPLIT, purpose,
                                                ClashResponse::NONE);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(tier);
      const SyRestraintKit<double, double2, double4> poly_rk =
        poly_ag.getDoublePrecisionRestraintKit(tier);
      const SyAtomUpdateKit<double, double2, double4> poly_auk =
        poly_ag.getDoublePrecisionAtomUpdateKit(tier);
      MMControlKit<double> ctrl = mmctrl->dpData(tier);
      ThermostatWriter tstw = heat_bath->dpData(tier);
      CacheResourceKit<double> gmem_r = tb_space->dpData(tier);
      launchValence(poly_vk, poly_rk, &ctrl, &poly_psw, poly_auk, &tstw, &scw, &gmem_r, eval_force,
                    eval_energy, purpose, bt, clash_distance, clash_ratio);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(tier);
      const SyRestraintKit<float, float2, float4> poly_rk =
        poly_ag.getSinglePrecisionRestraintKit(tier);
      const SyAtomUpdateKit<float, float2, float4> poly_auk =
        poly_ag.getSinglePrecisionAtomUpdateKit(tier);
      MMControlKit<float> ctrl = mmctrl->spData(tier);
      ThermostatWriter tstw = heat_bath->spData(tier);
      CacheResourceKit<float> gmem_r = tb_space->spData(tier);
      launchValence(poly_vk, poly_rk, &ctrl, &poly_psw, poly_auk, &tstw, &scw, &gmem_r, eval_force,
                    eval_energy, purpose, force_sum, bt, clash_distance, clash_ratio);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchValence(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                          MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                          Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                          const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                          const VwuGoal purpose, const CoreKlManager &launcher,
                          const double clash_distance, const double clash_ratio) {
  if (prec == PrecisionModel::DOUBLE || poly_ps->getForceAccumulationBits() <= 24) {
    launchValence(prec, poly_ag, mmctrl, poly_ps, heat_bath, sc, tb_space, eval_force, eval_energy,
                  purpose, AccumulationMethod::SPLIT, launcher, clash_distance, clash_ratio);
  }
  else {
    launchValence(prec, poly_ag, mmctrl, poly_ps, heat_bath, sc, tb_space, eval_force, eval_energy,
                  purpose, AccumulationMethod::WHOLE, launcher, clash_distance, clash_ratio);
  }
}

} // namespace energy
} // namespace stormm
