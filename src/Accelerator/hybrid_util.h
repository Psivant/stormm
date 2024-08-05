// -*-c++-*-
#ifndef STORMM_HYBRID_UTIL_H
#define STORMM_HYBRID_UTIL_H

#include "copyright.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "gpu_details.h"
#include "hybrid.h"
#ifdef STORMM_USE_HPC
#  include "hpc_hybrid_util.h"
#endif

namespace stormm {
namespace card {

using data_types::isHpcVectorType;
using data_types::isFloatingPointScalarType;
using data_types::isSignedIntegralScalarType;
using data_types::isUnsignedIntegralScalarType;
  
/// \brief A free function to test whether a given Hybrid format contains a valid CPU host memory
///        component and report what could be a more descriptive error message.  An error would
///        also be returned by a given Hybrid array, if accessed at an adress on the GPU which is
///        unavailable.  However, this would only return the char string associated with the
///        Hybrid itself and might not offer a sufficient clue as to what is wrong.
///
/// \param tfmt           The relevant format to test
/// \param message        The basic error message to display in the event of a failure
/// \param class_caller   The class of the object calling the function
/// \param method_caller  The method of the class object calling the function
void confirmCpuMemory(HybridFormat tfmt, const std::string &message,
                      const char* class_caller = nullptr,
                      const char* method_caller = nullptr);

/// \brief A free function to test whether a given Hybrid format contains a valid HPC component
///        and report what could be a more descriptive error message.  Descriptions of input
///        parameters follow from confirmCpuMemory(), above.
void confirmGpuMemory(HybridFormat tfmt, const std::string &message,
                      const char* class_caller = nullptr,
                      const char* method_caller = nullptr);

#ifdef STORMM_USE_HPC
/// \brief A free function to test whether a block of host memory might be visible to the GPU.
///        Descriptions of input parameters follow from confirmCpuMemory(), above.
void confirmHostVisibleToGpu(HybridFormat tfmt, const std::string &message,
                             const char* class_caller, const char* method_caller);

/// \brief Assess the copying between two Hybrid arrays as required to fulfill the work in
///        deepCopy() or deepRecast().  All boolean pointers represent single values which are
///        modified and returned.
///
/// \param dest_format  Memory layout of the destination Hybrid array
/// \param orig_format  Memory layout of the original Hybrid array
/// \param do_hhc       Return as true if host-to-host copying is needed
/// \param do_hdc       Return as true if host-to-device copying is needed
/// \param do_dhc       Return as true if device-to-host copying is needed
/// \param do_ddc       Return as true if device-to-device copying is needed
/// \param use_kernel   Return as true if all memory is visible to the GPU, permitting the use of
///                     a single kernel rather than multiple calls to cudaMemcpy() to perform the
///                     memory replication
void markCopyInstructions(HybridFormat dest_format, HybridFormat orig_format, bool *do_hhc,
                          bool *do_hdc, bool *do_dhc, bool *do_ddc, bool *use_kernel);
#endif

/// \brief A free function to test whether a particular request for memory at one level is valid.
///        This encapsulates the error message and can relay the calling object's class and member
///        function.  Descriptions of input parameters follow from confirmCpuMemory(), in addition
///        to:
///
/// \param request  The memory level at which data is requested
/// \param native   The memory layout with which data is present in the object
void checkFormatCompatibility(HybridTargetLevel request, HybridFormat native,
                              const char* class_caller, const char* method_caller);
  
/// \brief Check that the counds and offset of the arrays are valid for performing the requested
///        transfer of data between two Hybrid objects.
///
/// \param destination_size    Available size of the destination array
/// \param original_size       Size of the original array
/// \param length              The number of array indices to copy from the original array to the
///                            destination
/// \param destination_offset  Begin copying memory to this index in the destination array
/// \param original_offset     Begin copying memory from this index in the original array
size_t screenHybridTransfer(size_t destination_size, size_t original_size, size_t length,
                            size_t destination_offset, size_t original_offset, const char* caller);

/// \brief Perform a deep copy of the information between Hybrid objects which have the same
///        templated data type but may be of different memory layouts.
///
/// Overloaded:
///   - Provide the GPU information to initiate what is otherwise, by default, a complete copy of
///     the original array to the new array.
///   - Provide information about the copy bounds up to and including the GPU specifications.
///
/// \param destination         The destination array into which data shall be copied
/// \param original            The original array from which data will be taken
/// \param length              The number of array indices to copy from the original array to the
///                            destination
/// \param destination_offset  Begin copying memory to this index in the destination array
/// \param original_offset     Begin copying memory from this index in the original array
/// \{
template <typename T>
void deepCopy(Hybrid<T> *destination, const Hybrid<T> &original, size_t length = 0,
              size_t destination_offset = 0, size_t original_offset = 0,
              const GpuDetails &gpu = null_gpu);

template <typename T>
void deepCopy(Hybrid<T> *destination, const Hybrid<T> &original, const GpuDetails &gpu,
              size_t length = 0, size_t destination_offset = 0, size_t original_offset = 0);
/// \}

/// \brief Perform a deep copy of the information in one or a pair of Hybrid objects which have
///        different templated memory data types and may also be of different memory layouts.
///        This routine uses the same priority system for adapting memory formats as deepCopy(),
///        above.
///        
/// Overloaded:
///   - Accept a single Hybrid array as an input
///   - Accept a pair of Hybrid arrays, assumed to be the pimary and overflow accumulators of
///     split fixed-precision data
///   - Perform the same transfers but with triplets of inputs and outputs
///
/// Descriptions of input parameters follow from deepCopy(), above, in addition to:
///
/// \param overflow  Overflow accumulators for what is assumed to be split fixed-precision data
/// \{
template <typename Tdest, typename Torig>
void deepRecast(Hybrid<Tdest> *destination, const Hybrid<Torig> &original, double scale,
                size_t length = 0, size_t destination_offset = 0, size_t original_offset = 0,
                const GpuDetails &gpu = null_gpu);

template <typename Tdest, typename Torig>
void deepRecast(Hybrid<Tdest> *destination, const Hybrid<Torig> &original,
                const Hybrid<int> &overflow, double scale, size_t length = 0,
                size_t destination_offset = 0, size_t original_offset = 0,
                const GpuDetails &gpu = null_gpu);
/// \}

} // namespace card
} // namespace stormm

#include "hybrid_util.tpp"

#endif
