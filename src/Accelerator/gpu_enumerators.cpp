#include "copyright.h"
#include "gpu_enumerators.h"

namespace stormm {
namespace card {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const HybridKind input) {
  switch (input) {
  case HybridKind::ARRAY:
    return std::string("ARRAY");
  case HybridKind::POINTER:
    return std::string("POINTER");
  }   
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const HybridFormat input) {
  switch (input) {
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
    return std::string("EXPEDITED");
  case HybridFormat::DECOUPLED:
    return std::string("DECOUPLED");
  case HybridFormat::UNIFIED:
    return std::string("UNIFIED");
  case HybridFormat::HOST_ONLY:
    return std::string("HOST_ONLY");
  case HybridFormat::DEVICE_ONLY:
    return std::string("DEVICE_ONLY");
  case HybridFormat::HOST_MOUNTED:
    return std::string("HOST_MOUNTED");
#else
  case HybridFormat::HOST_ONLY:
    return std::string("HOST_ONLY");
#endif
  }   
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const HybridTargetLevel input) {
  switch (input) {
  case HybridTargetLevel::HOST:
    return std::string("HOST");
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return std::string("DEVICE");
#endif
  }   
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const HpcKernelSync input) {
  switch (input) {
  case HpcKernelSync::BEFORE:
    return std::string("BEFORE");
  case HpcKernelSync::AFTER:
    return std::string("AFTER");
  case HpcKernelSync::BEFORE_AND_AFTER:
    return std::string("BEFORE_AND_AFTER");
  case HpcKernelSync::NO_SYNC:
    return std::string("NO_SYNC");
  case HpcKernelSync::MEMORY_AUTO:
    return std::string("MEMORY_AUTO");
  }
  __builtin_unreachable();
}

} // namespace card
} // namespace stormm
