#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/vector_ops.h"
#include "Trajectory/coordinate_copy.h"
#include "cull_synthesis.h"

namespace stormm {
namespace synthesis {

using stmath::minValue;
using stmath::maxValue;
using trajectory::coordCopy;
  
//-------------------------------------------------------------------------------------------------
void validateCulling(const int system_count, const std::vector<int> &survivors) {
  const int min_idx = minValue(survivors);
  const int max_idx = maxValue(survivors);
  if (min_idx < 0 || max_idx >= system_count) {
    const int offender = (min_idx < 0) ? min_idx : max_idx;
    rtErr("A system index of " + std::to_string(offender) + " is invalid for a synthesis of " +
          std::to_string(system_count) + " systems.", "cullSynthesis");
  }
}
  
//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis cullSynthesis(const AtomGraphSynthesis *original,
                                 const std::vector<int> &survivors, const GpuDetails &gpu,
                                 const ExceptionResponse policy) {
  validateCulling(original->getSystemCount(), survivors);
  const size_t nsurv = survivors.size();
  std::vector<AtomGraph*> topl_list;
  std::vector<RestraintApparatus*> rstr_list;
  topl_list.reserve(nsurv);
  rstr_list.reserve(nsurv);
  for (size_t i = 0; i < nsurv; i++) {
    topl_list.push_back(const_cast<AtomGraph*>(original->getSystemTopologyPointer(survivors[i])));
    rstr_list.push_back(original->getSystemRestraintPointer(survivors[i]));
  }
  return AtomGraphSynthesis(topl_list, rstr_list, policy, gpu);
}

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis cullSynthesis(const AtomGraphSynthesis &original,
                                 const std::vector<int> &survivors, const GpuDetails &gpu,
                                 const ExceptionResponse policy) {
  return cullSynthesis(original.getSelfPointer(), survivors, gpu, policy);
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis cullSynthesis(const PhaseSpaceSynthesis *original,
                                  const std::vector<int> &survivors, const GpuDetails &gpu) {
  validateCulling(original->getSystemCount(), survivors);
  const int nsurv = survivors.size();
  std::vector<AtomGraph*> topl_list;
  std::vector<CoordinateFrame> coord_list;
  topl_list.reserve(nsurv);
  coord_list.reserve(nsurv);
  for (int i = 0; i < nsurv; i++) {
    topl_list.push_back(const_cast<AtomGraph*>(original->getSystemTopologyPointer(survivors[i])));

    // The coordinate list will be created in host-bound memory to conserve GPU resources, if
    // applicable.  The result will come in the format of the original synthesis object.  In
    // essence, the CoordinateFrame array serves only to size the result.
    coord_list.push_back(original->exportCoordinates(survivors[i], HybridFormat::HOST_ONLY));
  }
  std::vector<int2> system_pairs(nsurv);
  for (int i = 0; i < nsurv; i++) {
    system_pairs[i] = { survivors[i], i };
  }
  PhaseSpaceSynthesis result(coord_list, topl_list, original->getGlobalPositionBits(),
                             original->getLocalPositionBits(), original->getVelocityBits(),
                             original->getForceAccumulationBits(), original->getFormat(), gpu);
  std::vector<HybridTargetLevel> tier_list;  
  switch (original->getFormat()) {
  case HybridFormat::HOST_ONLY:
    tier_list.push_back(HybridTargetLevel::HOST);
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::UNIFIED:
    tier_list.push_back(HybridTargetLevel::HOST);
    break;
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    tier_list.push_back(HybridTargetLevel::HOST);
    tier_list.push_back(HybridTargetLevel::DEVICE);
    break;
  case HybridFormat::DEVICE_ONLY:
    tier_list.push_back(HybridTargetLevel::DEVICE);
    break;
#endif
  }

  // Copy coordinates from the original object to the result, preserving all information at a
  // bitwise level.  The conversion to double-precision in creating the list of CoordinateFrame
  // objects has no bearing on the information thus conveyed.
  for (size_t i = 0; i < tier_list.size(); i++) {
    coordCopy(&result, *original, system_pairs, tier_list[i], tier_list[i], gpu);
  }

  // Fast-forward the result to match the original object's stage in the coordinate cycle.
  result.updateCyclePosition(original->getCyclePosition());
  return result;
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis cullSynthesis(const PhaseSpaceSynthesis &original,
                                  const std::vector<int> &survivors, const GpuDetails &gpu) {
  return cullSynthesis(original.getSelfPointer(), survivors, gpu);
}

//-------------------------------------------------------------------------------------------------
SynthesisCacheMap remapSynthesis(const SynthesisCacheMap &original,
                                 const std::vector<int> &survivors) {

  // Determine which system cache index each member of the survivor population corresponds to.
  const size_t nsys = survivors.size();
  std::vector<int> new_correspondence(nsys);
  for (size_t i = 0; i < nsys; i++) {
    new_correspondence[i] = original.getSystemCacheIndex(survivors[i]);
  }
  return SynthesisCacheMap(new_correspondence, original.getCachePointer(),
                           original.getTopologySynthesisPointer(),
                           original.getCoordinateSynthesisPointer());
}
  
} // namespace synthesis
} // namespace stormm
