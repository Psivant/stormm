// -*-c++-*-

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/ptx_macros.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Synthesis/valence_workunit.h"
#include "hpc_integration.h"

namespace stormm {
namespace trajectory {

using card::HybridTargetLevel;
using energy::ValenceKernelSize;
using numerics::AccumulationMethod;
using numerics::getEnumerationName;
using synthesis::AtomGraphSynthesis;
using synthesis::eighth_valence_work_unit_atoms;
using synthesis::half_valence_work_unit_atoms;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::quarter_valence_work_unit_atoms;
using synthesis::vwu_abstract_length;
using synthesis::VwuAbstractMap;

#include "Accelerator/syncwarp.cui"
#include "Math/rounding.cui"
#include "Numerics/accumulation.cui"

#define CONSTRAINT_STANDALONE

// Double-precision floating point arithmetic
#define TCALC double
#  define SPLIT_FORCE_ACCUMULATION
#  define TCALC2 double2
#  define TCALC4 double4
#  define SQRT_FUNC sqrt
#  define ABS_FUNC fabs
#    define INTEG_KERNEL_THREAD_COUNT 512
#    define INTEG_BLOCK_MULTIPLICITY 2
#    define KERNEL_NAME kdsIntegrationVelCnst
#      include "velocity_constraints.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kdsIntegrationGeomCnst
#      include "geometry_constraints.cui"
#    undef KERNEL_NAME
#    undef INTEG_BLOCK_MULTIPLICITY
#    undef INTEG_KERNEL_THREAD_COUNT
#  undef SPLIT_FORCE_ACCUMULATION
#  undef TCALC2
#  undef TCALC4
#  undef SQRT_FUNC
#  undef ABS_FUNC
#undef TCALC

// Single-precision floating point arithmetic, split accumulation
#define TCALC float
#  define TCALC_IS_SINGLE
#  define TCALC2 float2
#  define TCALC4 float4
#  define SQRT_FUNC sqrtf
#  define ABS_FUNC fabsf
#  define LLCONV_FUNC __float2ll_rn
#  define SPLIT_FORCE_ACCUMULATION
#    define INTEG_KERNEL_THREAD_COUNT 512
#    define INTEG_BLOCK_MULTIPLICITY 2
#    define KERNEL_NAME kfsIntegrationVelCnstXL
#      include "velocity_constraints.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfsIntegrationGeomCnstXL
#      include "geometry_constraints.cui"
#    undef KERNEL_NAME
#    undef INTEG_BLOCK_MULTIPLICITY
#    undef INTEG_KERNEL_THREAD_COUNT
#    define INTEG_KERNEL_THREAD_COUNT 256
#    define INTEG_BLOCK_MULTIPLICITY 4
#    define KERNEL_NAME kfsIntegrationVelCnstLG
#      include "velocity_constraints.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfsIntegrationGeomCnstLG
#      include "geometry_constraints.cui"
#    undef KERNEL_NAME
#    undef INTEG_BLOCK_MULTIPLICITY
#    undef INTEG_KERNEL_THREAD_COUNT
#    define INTEG_KERNEL_THREAD_COUNT 128
#    define INTEG_BLOCK_MULTIPLICITY 8
#    define KERNEL_NAME kfsIntegrationVelCnstMD
#      include "velocity_constraints.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfsIntegrationGeomCnstMD
#      include "geometry_constraints.cui"
#    undef KERNEL_NAME
#    undef INTEG_BLOCK_MULTIPLICITY
#    undef INTEG_KERNEL_THREAD_COUNT
#    define INTEG_KERNEL_THREAD_COUNT 64
#    define INTEG_BLOCK_MULTIPLICITY 16
#    define KERNEL_NAME kfsIntegrationVelCnstSM
#      include "velocity_constraints.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfsIntegrationGeomCnstSM
#      include "geometry_constraints.cui"
#    undef KERNEL_NAME
#    undef INTEG_BLOCK_MULTIPLICITY
#    undef INTEG_KERNEL_THREAD_COUNT
#  undef SPLIT_FORCE_ACCUMULATION
#  define INTEG_KERNEL_THREAD_COUNT 512
#  define INTEG_BLOCK_MULTIPLICITY 2
#  define KERNEL_NAME kfIntegrationVelCnstXL
#    include "velocity_constraints.cui"
#  undef KERNEL_NAME
#  define KERNEL_NAME kfIntegrationGeomCnstXL
#    include "geometry_constraints.cui"
#  undef KERNEL_NAME
#  undef INTEG_BLOCK_MULTIPLICITY
#  undef INTEG_KERNEL_THREAD_COUNT
#  define INTEG_KERNEL_THREAD_COUNT 256
#  define INTEG_BLOCK_MULTIPLICITY 4
#  define KERNEL_NAME kfIntegrationVelCnstLG
#    include "velocity_constraints.cui"
#  undef KERNEL_NAME
#  define KERNEL_NAME kfIntegrationGeomCnstLG
#    include "geometry_constraints.cui"
#  undef KERNEL_NAME
#  undef INTEG_BLOCK_MULTIPLICITY
#  undef INTEG_KERNEL_THREAD_COUNT
#  define INTEG_KERNEL_THREAD_COUNT 128
#  define INTEG_BLOCK_MULTIPLICITY 8
#  define KERNEL_NAME kfIntegrationVelCnstMD
#    include "velocity_constraints.cui"
#  undef KERNEL_NAME
#  define KERNEL_NAME kfIntegrationGeomCnstMD
#    include "geometry_constraints.cui"
#  undef KERNEL_NAME
#  undef INTEG_BLOCK_MULTIPLICITY
#  undef INTEG_KERNEL_THREAD_COUNT
#  define INTEG_KERNEL_THREAD_COUNT 64
#  define INTEG_BLOCK_MULTIPLICITY 16
#  define KERNEL_NAME kfIntegrationVelCnstSM
#    include "velocity_constraints.cui"
#  undef KERNEL_NAME
#  define KERNEL_NAME kfIntegrationGeomCnstSM
#    include "geometry_constraints.cui"
#  undef KERNEL_NAME
#  undef INTEG_BLOCK_MULTIPLICITY
#  undef INTEG_KERNEL_THREAD_COUNT
#  undef TCALC_IS_SINGLE
#  undef TCALC2
#  undef TCALC4
#  undef SQRT_FUNC
#  undef ABS_FUNC
#  undef LLCONV_FUNC
#undef TCALC

#undef CONSTRAINT_STANDALONE
  
//-------------------------------------------------------------------------------------------------
cudaFuncAttributes queryIntegrationKernelRequirements(const PrecisionModel prec,
                                                      const AccumulationMethod acc_meth,
                                                      const ValenceKernelSize kwidth,
                                                      const IntegrationStage process) {
  cudaFuncAttributes result;
  cudaError_t cfa;

  switch (prec) {
  case PrecisionModel::DOUBLE:
    switch (process) {
    case IntegrationStage::VELOCITY_ADVANCE:
      break;
    case IntegrationStage::VELOCITY_CONSTRAINT:
      cfa = cudaFuncGetAttributes(&result, kdsIntegrationVelCnst);
      break;
    case IntegrationStage::POSITION_ADVANCE:
      break;
    case IntegrationStage::GEOMETRY_CONSTRAINT:
      cfa = cudaFuncGetAttributes(&result, kdsIntegrationGeomCnst);
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
      switch (kwidth) {
      case ValenceKernelSize::XL:
        switch (process) {
        case IntegrationStage::VELOCITY_ADVANCE:
          break;
        case IntegrationStage::VELOCITY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfsIntegrationVelCnstXL);
          break;
        case IntegrationStage::POSITION_ADVANCE:
          break;
        case IntegrationStage::GEOMETRY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfsIntegrationGeomCnstXL);
          break;
        }
        break;
      case ValenceKernelSize::LG:
        switch (process) {
        case IntegrationStage::VELOCITY_ADVANCE:
          break;
        case IntegrationStage::VELOCITY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfsIntegrationVelCnstLG);
          break;
        case IntegrationStage::POSITION_ADVANCE:
          break;
        case IntegrationStage::GEOMETRY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfsIntegrationGeomCnstLG);
          break;
        }
        break;
      case ValenceKernelSize::MD:
        switch (process) {
        case IntegrationStage::VELOCITY_ADVANCE:
          break;
        case IntegrationStage::VELOCITY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfsIntegrationVelCnstMD);
          break;
        case IntegrationStage::POSITION_ADVANCE:
          break;
        case IntegrationStage::GEOMETRY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfsIntegrationGeomCnstMD);
          break;
        }
        break;
      case ValenceKernelSize::SM:
        switch (process) {
        case IntegrationStage::VELOCITY_ADVANCE:
          break;
        case IntegrationStage::VELOCITY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfsIntegrationVelCnstSM);
          break;
        case IntegrationStage::POSITION_ADVANCE:
          break;
        case IntegrationStage::GEOMETRY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfsIntegrationGeomCnstSM);
          break;
        }
        break;
      }
      break;
    case AccumulationMethod::WHOLE:
      switch (kwidth) {
      case ValenceKernelSize::XL:
        switch (process) {
        case IntegrationStage::VELOCITY_ADVANCE:
          break;
        case IntegrationStage::VELOCITY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfIntegrationVelCnstXL);
          break;
        case IntegrationStage::POSITION_ADVANCE:
          break;
        case IntegrationStage::GEOMETRY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfIntegrationGeomCnstXL);
          break;
        }
        break;
      case ValenceKernelSize::LG:
        switch (process) {
        case IntegrationStage::VELOCITY_ADVANCE:
          break;
        case IntegrationStage::VELOCITY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfIntegrationVelCnstLG);
          break;
        case IntegrationStage::POSITION_ADVANCE:
          break;
        case IntegrationStage::GEOMETRY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfIntegrationGeomCnstLG);
          break;
        }
        break;
      case ValenceKernelSize::MD:
        switch (process) {
        case IntegrationStage::VELOCITY_ADVANCE:
          break;
        case IntegrationStage::VELOCITY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfIntegrationVelCnstMD);
          break;
        case IntegrationStage::POSITION_ADVANCE:
          break;
        case IntegrationStage::GEOMETRY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfIntegrationGeomCnstMD);
          break;
        }
        break;
      case ValenceKernelSize::SM:
        switch (process) {
        case IntegrationStage::VELOCITY_ADVANCE:
          break;
        case IntegrationStage::VELOCITY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfIntegrationVelCnstSM);
          break;
        case IntegrationStage::POSITION_ADVANCE:
          break;
        case IntegrationStage::GEOMETRY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kfIntegrationGeomCnstSM);
          break;
        }
        break;
      }
      break;
    case AccumulationMethod::AUTOMATIC:
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
      case AccumulationMethod::AUTOMATIC:
        rtErr("No kernel is available for " + getEnumerationName(acc_meth) +
              " force accumulation.", "queryIntegrationKernelRequirements");
        break;
      }
      break;
    }
    error_message += "Integration";
    switch (process) {
    case IntegrationStage::VELOCITY_ADVANCE:
      error_message += "VelAdv";
      break;
    case IntegrationStage::VELOCITY_CONSTRAINT:
      error_message += "VelCnst";
      break;
    case IntegrationStage::POSITION_ADVANCE:
      error_message += "PosAdv";
      break;
    case IntegrationStage::GEOMETRY_CONSTRAINT:
      error_message += "GeomCnst";
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

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *tb_resw,
                              MMControlKit<double> *ctrl, const SyValenceKit<double> &poly_vk,
                              const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                              const ThermostatWriter<double> &tstw, const int2 lp,
                              const IntegrationStage process) {
  switch (process) {
  case IntegrationStage::VELOCITY_ADVANCE:
    break;
  case IntegrationStage::VELOCITY_CONSTRAINT:
    kdsIntegrationVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
    break;
  case IntegrationStage::POSITION_ADVANCE:
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    kdsIntegrationGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *tb_resw,
                              MMControlKit<float> *ctrl, const SyValenceKit<float> &poly_vk,
                              const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                              const ThermostatWriter<float> &tstw, const int2 lp,
                              const AccumulationMethod acc_meth, const ValenceKernelSize kwidth,
                              const IntegrationStage process) {
  switch (process) {
  case IntegrationStage::VELOCITY_ADVANCE:
    break;
  case IntegrationStage::VELOCITY_CONSTRAINT:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      switch (kwidth) {
      case ValenceKernelSize::XL:
        kfsIntegrationVelCnstXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                                *tb_resw);
        break;
      case ValenceKernelSize::LG:
        kfsIntegrationVelCnstLG<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                                *tb_resw);
        break;
      case ValenceKernelSize::MD:
        kfsIntegrationVelCnstMD<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                                *tb_resw);
        break;
      case ValenceKernelSize::SM:
        kfsIntegrationVelCnstSM<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                                *tb_resw);
        break;
      }
      break;
    case AccumulationMethod::WHOLE:
      switch (kwidth) {
      case ValenceKernelSize::XL:
        kfIntegrationVelCnstXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                               *tb_resw);
        break;
      case ValenceKernelSize::LG:
        kfIntegrationVelCnstLG<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                               *tb_resw);
        break;
      case ValenceKernelSize::MD:
        kfIntegrationVelCnstMD<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                               *tb_resw);
        break;
      case ValenceKernelSize::SM:
        kfIntegrationVelCnstSM<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                               *tb_resw);
        break;
      }
      break;
    }
    break;
  case IntegrationStage::POSITION_ADVANCE:
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      switch (kwidth) {
      case ValenceKernelSize::XL:
        kfsIntegrationGeomCnstXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                                 *tb_resw);
        break;
      case ValenceKernelSize::LG:
        kfsIntegrationGeomCnstLG<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                                 *tb_resw);
        break;
      case ValenceKernelSize::MD:
        kfsIntegrationGeomCnstMD<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                                 *tb_resw);
        break;
      case ValenceKernelSize::SM:
        kfsIntegrationGeomCnstSM<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                                 *tb_resw);
        break;
      }
      break;
    case AccumulationMethod::WHOLE:
      switch (kwidth) {
      case ValenceKernelSize::XL:
        kfIntegrationGeomCnstXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                                *tb_resw);
        break;
      case ValenceKernelSize::LG:
        kfIntegrationGeomCnstLG<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                                *tb_resw);
        break;
      case ValenceKernelSize::MD:
        kfIntegrationGeomCnstMD<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                                *tb_resw);
        break;
      case ValenceKernelSize::SM:
        kfIntegrationGeomCnstSM<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw,
                                                *tb_resw);
        break;
      }
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PhaseSpaceSynthesis *poly_ps, CacheResource *tb_res, Thermostat *tst,
                              MolecularMechanicsControls *mmctrl,
                              const AtomGraphSynthesis &poly_ag, const CoreKlManager &launcher,
                              const PrecisionModel prec, const AccumulationMethod acc_meth,
                              const IntegrationStage process) {
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(devc_tier);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      CacheResourceKit<double> tb_resw = tb_res->dpData(devc_tier);
      MMControlKit<double> ctrl = mmctrl->dpData(devc_tier);
      SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(devc_tier);
      SyAtomUpdateKit<double,
                      double2,
                      double4> poly_auk = poly_ag.getDoublePrecisionAtomUpdateKit(devc_tier);
      ThermostatWriter<double> tstw = tst->dpData(devc_tier);
      const int2 lp = launcher.getIntegrationKernelDims(prec, acc_meth, process);
      launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, poly_vk, poly_auk, tstw, lp, process);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      CacheResourceKit<float> tb_resw = tb_res->spData(devc_tier);
      MMControlKit<float> ctrl = mmctrl->spData(devc_tier);
      SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(devc_tier);
      SyAtomUpdateKit<float,
                      float2,
                      float4> poly_auk = poly_ag.getSinglePrecisionAtomUpdateKit(devc_tier);
      ThermostatWriter<float> tstw = tst->spData(devc_tier);
      const int2 lp = launcher.getIntegrationKernelDims(prec, acc_meth, process);
      const ValenceKernelSize kwidth = poly_ag.getValenceThreadBlockSize();
      launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, poly_vk, poly_auk, tstw, lp, acc_meth,
                               kwidth, process);
    }
    break;
  }
}

} // namespace trajectory
} // namespace stormm
