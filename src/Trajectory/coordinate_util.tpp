// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void checkCopyValidity(const Tdest *destination, const Torig &origin,
                       const HybridTargetLevel destination_tier,
                       const HybridTargetLevel origin_tier) {
  const std::string dest_obj = nameCoordinateType(std::type_index(typeid(Tdest)).hash_code());
  const std::string orig_obj = nameCoordinateType(std::type_index(typeid(Torig)).hash_code());
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    confirmCpuMemory(destination->getFormat(), "The destination object (" + dest_obj +
                     ") does not have memory allocated on the CPU host to receive coordinates.",
                     "coordCopy");
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    confirmGpuMemory(destination->getFormat(), "The destination object (" + dest_obj +
                     ") does not have memory allocated on the GPU device to receive coordinates.",
                     "coordCopy");
    break;
#endif
  }
  switch (origin_tier) {
  case HybridTargetLevel::HOST:
    confirmCpuMemory(origin.getFormat(), "The origin object (" + orig_obj + ") does not have "
                     "memory allocated on the CPU host to provide coordinates.", "coordCopy");
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    confirmGpuMemory(origin.getFormat(), "The origin object (" + orig_obj + ") does not have "
                     "memory allocated on the GPU device to provide coordinates.", "coordCopy");
    break;
#endif
  }
}

} // namespace trajectory
} // namespace stormm
