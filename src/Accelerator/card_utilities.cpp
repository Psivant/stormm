#include "copyright.h"
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
#    include <cuda_runtime.h>
#  endif
#endif
#include "Reporting/error_format.h"
#include "card_utilities.h"

namespace stormm {
namespace card {

//-------------------------------------------------------------------------------------------------
void launchPreparation(const HpcKernelSync sync) {
#ifdef STORMM_USE_HPC
  switch (sync) {
  case HpcKernelSync::BEFORE:
  case HpcKernelSync::BEFORE_AND_AFTER:
#  ifdef STORMM_USE_CUDA
    cudaDeviceSynchronize();
#endif
    break;
  case HpcKernelSync::AFTER:
  case HpcKernelSync::NO_SYNC:
    break;
  case HpcKernelSync::MEMORY_AUTO:
    rtErr("The " + getEnumerationName(sync) + " option for synchronization requires specification "
          "of the intended memory destination and origin.  Use an alternate, overloaded form of "
          "this function to invoke such behavior.", "launchPreparation");
  }
#endif
}

//-------------------------------------------------------------------------------------------------
void launchPreparation(const HpcKernelSync sync, const HybridTargetLevel memory_dest,
                       const HybridTargetLevel memory_orig) {
#ifdef STORMM_USE_HPC
  switch (sync) {
  case HpcKernelSync::BEFORE:
  case HpcKernelSync::BEFORE_AND_AFTER:
#  ifdef STORMM_USE_CUDA
    cudaDeviceSynchronize();
#endif
    break;
  case HpcKernelSync::AFTER:
  case HpcKernelSync::NO_SYNC:
    break;
  case HpcKernelSync::MEMORY_AUTO:
    switch (memory_dest) {
    case HybridTargetLevel::HOST:
      break;
    case HybridTargetLevel::DEVICE:
      switch (memory_orig) {
      case HybridTargetLevel::HOST:
#  ifdef STORMM_USE_CUDA
        cudaDeviceSynchronize();
#  endif
        break;
      case HybridTargetLevel::DEVICE:
        break;
      }
      break;
    }
  }
#endif
}

//-------------------------------------------------------------------------------------------------
void launchResolution(const HpcKernelSync sync) {
#ifdef STORMM_USE_HPC
  switch (sync) {
  case HpcKernelSync::AFTER:
  case HpcKernelSync::BEFORE_AND_AFTER:
#  ifdef STORMM_USE_CUDA
    cudaDeviceSynchronize();
#endif
    break;
  case HpcKernelSync::BEFORE:
  case HpcKernelSync::NO_SYNC:
    break;
  case HpcKernelSync::MEMORY_AUTO:
    rtErr("The " + getEnumerationName(sync) + " option for synchronization requires specification "
          "of the intended memory destination and origin.  Use an alternate, overloaded form of "
          "this function to invoke such behavior.", "launchPreparation");
  }
#endif
}

//-------------------------------------------------------------------------------------------------
void launchResolution(const HpcKernelSync sync, const HybridTargetLevel memory_dest,
                      const HybridTargetLevel memory_orig) {
#ifdef STORMM_USE_HPC
  switch (sync) {
  case HpcKernelSync::AFTER:
  case HpcKernelSync::BEFORE_AND_AFTER:
#  ifdef STORMM_USE_CUDA
    cudaDeviceSynchronize();
#endif
    break;
  case HpcKernelSync::BEFORE:
  case HpcKernelSync::NO_SYNC:
    break;
  case HpcKernelSync::MEMORY_AUTO:
    switch (memory_dest) {
    case HybridTargetLevel::HOST:
      switch (memory_orig) {
      case HybridTargetLevel::HOST:
        break;
      case HybridTargetLevel::DEVICE:
#  ifdef STORMM_USE_CUDA
        cudaDeviceSynchronize();
#  endif
        break;
      }
      break;
    case HybridTargetLevel::DEVICE:
      break;
    }
    break;
  }
#endif
}

} // namespace card
} // namespace stormm
