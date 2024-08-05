// -*-c++-*-
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/stormm_vector_types.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Synthesis/valence_workunit.h"
#include "Topology/atomgraph_enumerators.h"
#include "hpc_virtual_site_handling.h"

namespace stormm {
namespace structure {

using card::HybridTargetLevel;
using synthesis::eighth_valence_work_unit_atoms;
using synthesis::half_valence_work_unit_atoms;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::quarter_valence_work_unit_atoms;
using synthesis::vwu_abstract_length;
using synthesis::VwuAbstractMap;
using topology::VirtualSiteKind;
using trajectory::CoordinateCycle;
  
#include "Math/rounding.cui"
#include "Math/vector_formulas.cui"
#include "Numerics/accumulation.cui"
  
#define VSITE_STANDALONE
#define TCALC  float
#  define TCALC2 float2
#  define TCALC3 float3
#  define TCALC4 float4
#  define TCALC_IS_SINGLE
#  define COS_FUNC cosf
#  define SIN_FUNC sinf
#  define SQRT_FUNC sqrtf
#  define LLCONV_FUNC __float2ll_rn
#  define VSITE_STANDALONE_BLOCK_MULTIPLIER 2
#    define VSITE_STANDALONE_THREAD_COUNT 512
#      define KERNEL_NAME kfPlaceVirtualSitesXL
#        include "virtual_site_placement.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfTransmitVSiteForcesXL
#        include "virtual_site_transmission.cui"
#      undef KERNEL_NAME
#    undef VSITE_STANDALONE_THREAD_COUNT
#  undef VSITE_STANDALONE_BLOCK_MULTIPLIER
#  define VSITE_STANDALONE_BLOCK_MULTIPLIER 4
#    define VSITE_STANDALONE_THREAD_COUNT 256
#      define KERNEL_NAME kfPlaceVirtualSitesLG
#        include "virtual_site_placement.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfTransmitVSiteForcesLG
#        include "virtual_site_transmission.cui"
#      undef KERNEL_NAME
#    undef VSITE_STANDALONE_THREAD_COUNT
#  undef VSITE_STANDALONE_BLOCK_MULTIPLIER
#  define VSITE_STANDALONE_BLOCK_MULTIPLIER 8
#    define VSITE_STANDALONE_THREAD_COUNT 128
#      define KERNEL_NAME kfPlaceVirtualSitesMD
#        include "virtual_site_placement.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfTransmitVSiteForcesMD
#        include "virtual_site_transmission.cui"
#      undef KERNEL_NAME
#    undef VSITE_STANDALONE_THREAD_COUNT
#  undef VSITE_STANDALONE_BLOCK_MULTIPLIER
#  define VSITE_STANDALONE_BLOCK_MULTIPLIER 16
#    define VSITE_STANDALONE_THREAD_COUNT 64
#      define KERNEL_NAME kfPlaceVirtualSitesSM
#        include "virtual_site_placement.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfTransmitVSiteForcesSM
#        include "virtual_site_transmission.cui"
#      undef KERNEL_NAME
#    undef VSITE_STANDALONE_THREAD_COUNT
#  undef VSITE_STANDALONE_BLOCK_MULTIPLIER
#  undef TCALC2
#  undef TCALC3
#  undef TCALC4
#  undef TCALC_IS_SINGLE
#  undef COS_FUNC
#  undef SIN_FUNC
#  undef SQRT_FUNC
#  undef LLCONV_FUNC
#undef TCALC

#define TCALC  double
#  define TCALC2 double2
#  define TCALC3 double3
#  define TCALC4 double4
#  define COS_FUNC cos
#  define SIN_FUNC sin
#  define SQRT_FUNC sqrt
#  define SPLIT_FORCE_ACCUMULATION
#  define VSITE_STANDALONE_BLOCK_MULTIPLIER 2
#  define VSITE_STANDALONE_THREAD_COUNT 256
#    define KERNEL_NAME kdPlaceVirtualSites
#      include "virtual_site_placement.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kdTransmitVSiteForces
#      include "virtual_site_transmission.cui"
#    undef KERNEL_NAME
#  undef VSITE_STANDALONE_THREAD_COUNT
#  undef VSITE_STANDALONE_BLOCK_MULTIPLIER
#  undef SPLIT_FORCE_ACCUMULATION
#  undef TCALC2
#  undef TCALC3
#  undef TCALC4
#  undef COS_FUNC
#  undef SIN_FUNC
#  undef SQRT_FUNC
#undef TCALC
#undef VSITE_STANDALONE

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes queryVirtualSiteKernelRequirements(const PrecisionModel prec,
                                                             const VirtualSiteActivity purpose,
                                                             const ValenceKernelSize kwidth) {
  cudaFuncAttributes result;
  cudaError_t cfa;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    switch (purpose) {
    case VirtualSiteActivity::PLACEMENT:
      cfa = cudaFuncGetAttributes(&result, kdPlaceVirtualSites);
      break;
    case VirtualSiteActivity::TRANSMIT_FORCES:
      cfa = cudaFuncGetAttributes(&result, kdTransmitVSiteForces);
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (purpose) {
    case VirtualSiteActivity::PLACEMENT:
      switch (kwidth) {
      case ValenceKernelSize::XL:
        cfa = cudaFuncGetAttributes(&result, kfPlaceVirtualSitesXL);
        break;
      case ValenceKernelSize::LG:
        cfa = cudaFuncGetAttributes(&result, kfPlaceVirtualSitesLG);
        break;
      case ValenceKernelSize::MD:
        cfa = cudaFuncGetAttributes(&result, kfPlaceVirtualSitesMD);
        break;
      case ValenceKernelSize::SM:
        cfa = cudaFuncGetAttributes(&result, kfPlaceVirtualSitesSM);
        break;
      }
      break;
    case VirtualSiteActivity::TRANSMIT_FORCES:
      switch (kwidth) {
      case ValenceKernelSize::XL:
        cfa = cudaFuncGetAttributes(&result, kfTransmitVSiteForcesXL);
        break;
      case ValenceKernelSize::LG:
        cfa = cudaFuncGetAttributes(&result, kfTransmitVSiteForcesLG);
        break;
      case ValenceKernelSize::MD:
        cfa = cudaFuncGetAttributes(&result, kfTransmitVSiteForcesMD);
        break;
      case ValenceKernelSize::SM:
        cfa = cudaFuncGetAttributes(&result, kfTransmitVSiteForcesSM);
        break;
      }
      break;
    }
    break;
  }
  if (cfa != cudaSuccess) {
    std::string error_message("Error obtaining attributes for kernel k");
    switch (prec) {
    case PrecisionModel::DOUBLE:
      error_message += "d";
      break;
    case PrecisionModel::SINGLE:
      error_message += "f";
      break;
    }
    switch (purpose) {
    case VirtualSiteActivity::PLACEMENT:
      error_message += "PlaceVirtualSite";
      break;
    case VirtualSiteActivity::TRANSMIT_FORCES:
      error_message += "TransmitVSiteForces";
      break;
    }
    switch (prec) {
    case PrecisionModel::DOUBLE:
      break;
    case PrecisionModel::SINGLE:
      error_message += getEnumerationName(kwidth);
      break;
    }
    error_message += ".";
    rtErr(error_message, "queryVirtualSiteKernelRequirements");
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
void launchVirtualSitePlacement(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *gmem_r,
                                const SyValenceKit<double> &poly_vk,
                                const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                                const int2 bt) {
  kdPlaceVirtualSites<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
}

//-------------------------------------------------------------------------------------------------
void launchVirtualSitePlacement(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *gmem_r,
                                const SyValenceKit<float> &poly_vk,
                                const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                                const int2 bt) {

  // Like the valence kernels themselves, the launch grid can be used to determine the appropriate
  // kernel to call.
  if (bt.y > 256) {
    kfPlaceVirtualSitesXL<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
  }
  else if (bt.y > 128) {
    kfPlaceVirtualSitesLG<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
  }
  else if (bt.y > 64) {
    kfPlaceVirtualSitesMD<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
  }
  else {
    kfPlaceVirtualSitesSM<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
  }
}

//-------------------------------------------------------------------------------------------------
void launchTransmitVSiteForces(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *gmem_r,
                               const SyValenceKit<double> &poly_vk,
                               const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                               const int2 bt) {
  kdTransmitVSiteForces<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
}

//-------------------------------------------------------------------------------------------------
void launchTransmitVSiteForces(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *gmem_r,
                               const SyValenceKit<float> &poly_vk,
                               const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                               const int2 bt) {
  if (bt.y > 256) {
    kfTransmitVSiteForcesXL<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
  }
  else if (bt.y > 128) {
    kfTransmitVSiteForcesLG<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
  }
  else if (bt.y > 64) {
    kfTransmitVSiteForcesMD<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
  }
  else {
    kfTransmitVSiteForcesSM<<<bt.x, bt.y>>>(*poly_psw, poly_vk, poly_auk, *gmem_r);
  }
}

//-------------------------------------------------------------------------------------------------
void launchVirtualSiteHandling(const PrecisionModel prec, const VirtualSiteActivity purpose,
                               PhaseSpaceSynthesis *poly_ps, CacheResource *tb_space,
                               const AtomGraphSynthesis &poly_ag, const CoreKlManager &launcher) {

  // Bail out if there are no virtual sites in the synthesis
  if (poly_ag.getVirtualSiteCount() == 0) {
    return;
  }
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(devc_tier);
  const int2 lp = launcher.getVirtualSiteKernelDims(prec, VirtualSiteActivity::PLACEMENT);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      CacheResourceKit<double> gmem_r = tb_space->dpData(devc_tier);
      const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(devc_tier);
      const SyAtomUpdateKit<double, double2, double4> poly_auk =
        poly_ag.getDoublePrecisionAtomUpdateKit(devc_tier);
      switch (purpose) {
      case VirtualSiteActivity::PLACEMENT:
        launchVirtualSitePlacement(&poly_psw, &gmem_r, poly_vk, poly_auk, lp);
        break;
      case VirtualSiteActivity::TRANSMIT_FORCES:
        launchTransmitVSiteForces(&poly_psw, &gmem_r, poly_vk, poly_auk, lp);
        break;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      CacheResourceKit<float> gmem_r = tb_space->spData(devc_tier);
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(devc_tier);
      const SyAtomUpdateKit<float, float2, float4> poly_auk =
        poly_ag.getSinglePrecisionAtomUpdateKit(devc_tier);
      switch (purpose) {
      case VirtualSiteActivity::PLACEMENT:
        launchVirtualSitePlacement(&poly_psw, &gmem_r, poly_vk, poly_auk, lp);
        break;
      case VirtualSiteActivity::TRANSMIT_FORCES:
        launchTransmitVSiteForces(&poly_psw, &gmem_r, poly_vk, poly_auk, lp);
        break;
      }
    }
    break;
  }
}

} // namespace structure
} // namespace stormm
