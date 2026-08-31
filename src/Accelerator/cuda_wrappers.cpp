#include <string>
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CDUA
#    include <cuda_runtime.h>
#    include <cufft.h>
#  endif
#endif
#include "copyright.h"
#include "cuda_wrappers.h"

#ifdef STORMM_USE_HPC

namespace stormm {
namespace card {

#  ifdef STORMM_USE_CUDA
  
//-------------------------------------------------------------------------------------------------
std::string getHpcErrorString(const cufftResult event, const HpcErrorVerbosity style) {
  switch (style) {
  case HpcErrorVerbosity::REASONING:
    switch (event) {
    case CUFFT_SUCCESS:
      return std::string("A cuFFT success");
    case CUFFT_INVALID_PLAN:
      return std::string("The plan parameter is not a valid handle");
    case CUFFT_ALLOC_FAILED:
      return std::string("The allocation of GPU or CPU memory for the plan failed");
    case CUFFT_INVALID_TYPE:
      return std::string("The data was of an invalid type for cuFFT");
    case CUFFT_INVALID_VALUE:
      return std::string("One or more invalid parameters were passed to the API");
    case CUFFT_INTERNAL_ERROR:
      return std::string("An internal driver error was detected");
    case CUFFT_EXEC_FAILED:
      return std::string("cuFFT failed to execute the transform on the GPU");
    case CUFFT_SETUP_FAILED:
      return std::string("The cuFFT library failed to initialize");
    case CUFFT_INVALID_SIZE:
      return std::string("One or more of the parameters is not a supported size");
    case CUFFT_UNALIGNED_DATA:
      return std::string("The data was not aligned for cuFFT");
    case CUFFT_INCOMPLETE_PARAMETER_LIST:
      return std::string("Missing parameters in call");
    case CUFFT_INVALID_DEVICE:
      return std::string("An invalid GPU index was specified in a descriptor or Execution of a "
                         "plan was on different GPU than plan creation");
    case CUFFT_PARSE_ERROR:
      return std::string("Internal plan database error");
    case CUFFT_NO_WORKSPACE:
      return std::string("No workspace has been provided prior to plan execution");
    case CUFFT_NOT_IMPLEMENTED:
      return std::string("Function does not implement functionality for parameters given");
    case CUFFT_LICENSE_ERROR:
      return std::string("Used in previous versions");
    case CUFFT_NOT_SUPPORTED:
      return std::string("Operation is not supported for parameters given");
    }
    break;
  case HpcErrorVerbosity::ENUMERATION:
    switch (event) {
    case CUFFT_SUCCESS:
      return std::string("CUFFT_SUCCESS");
    case CUFFT_INVALID_PLAN:
      return std::string("CUFFT_INVALID_PLAN");
    case CUFFT_ALLOC_FAILED:
      return std::string("CUFFT_ALLOC_FAILED");
    case CUFFT_INVALID_TYPE:
      return std::string("CUFFT_INVALID_TYPE");
    case CUFFT_INVALID_VALUE:
      return std::string("CUFFT_INVALID_VALUE");
    case CUFFT_INTERNAL_ERROR:
      return std::string("CUFFT_INTERNAL_ERROR");
    case CUFFT_EXEC_FAILED:
      return std::string("CUFFT_EXEC_FAILED");
    case CUFFT_SETUP_FAILED:
      return std::string("CUFFT_SETUP_FAILED");
    case CUFFT_INVALID_SIZE:
      return std::string("CUFFT_INVALID_SIZE");
    case CUFFT_UNALIGNED_DATA:
      return std::string("CUFFT_UNALIGNED_DATA");
    case CUFFT_INCOMPLETE_PARAMETER_LIST:
      return std::string("CUFFT_INCOMPLETE_PARAMETER_LIST");
    case CUFFT_INVALID_DEVICE:
      return std::string("CUFFT_INVALID_DEVICE");
    case CUFFT_PARSE_ERROR:
      return std::string("CUFFT_PARSE_ERROR");
    case CUFFT_NO_WORKSPACE:
      return std::string("CUFFT_NO_WORKSPACE");
    case CUFFT_NOT_IMPLEMENTED:
      return std::string("CUFFT_NOT_IMPLEMENTED");
    case CUFFT_LICENSE_ERROR:
      return std::string("CUFFT_LICENSE_ERROR");
    case CUFFT_NOT_SUPPORTED:
      return std::string("CUFFT_NOT_SUPPORTED");
    }
    break;
  }
  __builtin_unreachable();
}

#endif // STORMM_USE_CUDA

} // namespace card
} // namespace stormm

#endif
