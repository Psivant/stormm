#include "copyright.h"
#include "hybrid_util.h"

namespace stormm {
namespace card {

//-------------------------------------------------------------------------------------------------
void confirmCpuMemory(const HybridFormat tfmt, const std::string &message,
                      const char* class_caller, const char* method_caller) {
  switch (tfmt) {
  case HybridFormat::HOST_ONLY:
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
    break;
  case HybridFormat::DEVICE_ONLY:
    rtErr(message, class_caller, method_caller);
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void confirmGpuMemory(const HybridFormat tfmt, const std::string &message,
                      const char* class_caller, const char* method_caller) {
  switch (tfmt) {
  case HybridFormat::HOST_ONLY:
    rtErr(message, class_caller, method_caller);
#ifdef STORMM_USE_HPC
  case HybridFormat::HOST_MOUNTED:
    rtErr(message, class_caller, method_caller);
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::DEVICE_ONLY:
    break;
#endif
  }
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void confirmHostVisibleToGpu(const HybridFormat tfmt, const std::string &message,
                             const char* class_caller, const char* method_caller) {
  switch (tfmt) {
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::EXPEDITED:
  case HybridFormat::UNIFIED:
    break;
  case HybridFormat::HOST_ONLY:
  case HybridFormat::DEVICE_ONLY:
  case HybridFormat::DECOUPLED:
    rtErr(message, class_caller, method_caller);
  }
}

//-------------------------------------------------------------------------------------------------
void markCopyInstructions(const HybridFormat dest_format, const HybridFormat orig_format,
                          bool *do_hhc, bool *do_hdc, bool *do_dhc, bool *do_ddc,
                          bool *use_kernel) {
  *do_hhc = false;
  *do_hdc = false;
  *do_dhc = false;
  *do_ddc = false;
  switch (orig_format) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
    switch (dest_format) {
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
      *do_hhc = true;
      *do_ddc = true;
      break;
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_ONLY:
    case HybridFormat::HOST_MOUNTED:

      // Take the host data as the basis of the copy if the original has data at both levels and
      // the destination entails unified memory.
      *do_hhc = true;
      break;
    case HybridFormat::DEVICE_ONLY:
      *do_ddc = true;
      break;
    }
    break;
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
    switch (dest_format) {
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
      *do_hhc = true;
      *do_hdc = true;
      break;
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_ONLY:
    case HybridFormat::HOST_MOUNTED:
      *do_hhc = true;
      break;
    case HybridFormat::DEVICE_ONLY:
      *do_hdc = true;
      break;
    }
    break;
  case HybridFormat::DEVICE_ONLY:
    switch (dest_format) {
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
      *do_dhc = true;
      *do_ddc = true;
      break;
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_ONLY:
    case HybridFormat::HOST_MOUNTED:
      *do_dhc = true;
      break;
    case HybridFormat::DEVICE_ONLY:
      *do_ddc = true;
      break;
    }
    break;
  }
  *use_kernel = (*do_hdc || *do_dhc || *do_ddc);
  switch (orig_format) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::DEVICE_ONLY:
    break;
  case HybridFormat::HOST_ONLY:
    *use_kernel = false;
    break;
  }
  switch (dest_format) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::DEVICE_ONLY:
    break;
  case HybridFormat::HOST_ONLY:
    *use_kernel = false;
    break;
  }
}
#endif

//-------------------------------------------------------------------------------------------------
void checkFormatCompatibility(const HybridTargetLevel request, const HybridFormat native,
                              const char* class_caller, const char* method_caller) {
  switch (request) {
  case HybridTargetLevel::HOST:
    confirmCpuMemory(native, "A request for " + getEnumerationName(request) + " data is "
                     "incompatible with memory in format " + getEnumerationName(native) + ".",
                     class_caller, method_caller);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    confirmGpuMemory(native, "A request for " + getEnumerationName(request) + " data is "
                     "incompatible with memory in format " + getEnumerationName(native) + ".",
                     class_caller, method_caller);
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
size_t screenHybridTransfer(const size_t destination_size, const size_t original_size,
                            const size_t length, const size_t destination_offset,
                            const size_t original_offset, const char* caller) {
  if (destination_offset >= destination_size) {
    rtErr("An offset of " + std::to_string(destination_offset) + " is inaccessible to a target "
          "array of " + std::to_string(destination_size) + " elements.", caller);
  }
  if (original_offset >= original_size) {
    rtErr("An offset of " + std::to_string(original_offset) + " is inaccessible to a base "
          "array of " + std::to_string(original_size) + " elements.", caller);
  }
  const size_t result = (length == 0) ? original_size - original_offset : length;
  if (destination_size < destination_offset + result) {
    rtErr("The target array of " + std::to_string(destination_size) + " cannot accommodate " +
          std::to_string(result) + " elements copied with a starting index of " +
          std::to_string(destination_offset) + ".", caller);
  }
  return result;
}

} // namespace card
} // namespace stormm
