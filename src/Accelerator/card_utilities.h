// -*-c++-*-
#ifndef STORMM_CARD_UTILITIES_H
#define STORMM_CARD_UTILITIES_H

#include <vector>
#include "copyright.h"
#include "gpu_enumerators.h"

namespace stormm {
namespace card {

/// \brief Synchronize the CPU host and GPU device at the start of a function that will launch a
///        kernel.
///
/// Overloaded:
///   - Provide a synchronization order which will be obeyed with no other considerations (the
///     MEMORY_AUTO enumeration will error out in this case)
///   - Provide indicators of the intended memory origin and destination to put the synchronization
///     order in context
///
/// \param sync         The synchronization order
/// \param memory_dest  Destination tier for memory that will be handled by the kernel to be
///                     launched
/// \param memory_orig  The origin tier of memory that will be handled by the kernel to be launched
/// \{
void launchPreparation(HpcKernelSync sync);

void launchPreparation(HpcKernelSync sync, HybridTargetLevel memory_dest,
                       HybridTargetLevel memory_orig);
/// \}

/// \brief Synchronize the CPU host and GPU device at the end of a function that has launched a
///        kernel.  Overloads and descriptions of parameters follow from launchPreparation() above.
/// \{
void launchResolution(HpcKernelSync sync);

void launchResolution(HpcKernelSync sync, HybridTargetLevel memory_dest,
                      HybridTargetLevel memory_orig);
/// \}

/// \brief Wrap the native memory copy function for the GPU device at hand (e.g. cudaMemcpy() for
///        an NVIDIA GPU, hipMemcpy for an AMD GPU), using STORMM's native enumerations.  Other
///        overloads of deepCopy() can be found for Hybrid objects in the hybrid_util.h library,
///        whereas this one works with any pointers and assumes no qualifications such as
///        restricted paging in the underlying memory.
///
/// Overloaded:
///   - Provide a mutable pointer to the destination, with element size and element count, to copy
///     into a trusted, pre-allocated space
///   - Copy a block of memory byte-for-byte into the contents of a Standard Template Library
///     vector, which is then returned.  The size of the element is implied by the template type
///     in this case, and the destination is implied to be on the CPU host.  Only the number of
///     elements and the location of the original memory are needed.
///
/// \param destination       Pointer to the memory which will be populated
/// \param origin            Pointer to the memory which will be used to populate the destination
/// \param element_size      Size of the data type being copied, in bytes
/// \param element_count     Number of elements in the array to be copied.  The total transfer 
///                          will be element_size * element_count bytes.
/// \param destination_tier  Indicate whether to copy data to the CPU host or the GPU device
/// \param origin_tier       Indicate whether to copy data from the CPU host or the GPU device
/// \param desc              Optional description of the purpose of the transfer
/// \{
void deepCopy(void* destination, const void* origin, size_t element_size, size_t element_count,
              HybridTargetLevel destination_tier, HybridTargetLevel origin_tier,
              const char* desc = nullptr);

template <typename T>
std::vector<T> deepCopy(const void* origin, size_t element_count, HybridTargetLevel origin_tier,
                        const char* desc = nullptr);
/// \}

} // namespace card
} // namespace stormm

// Include the launch guards in any other STORMM libraries
namespace stormm {
  using card::launchPreparation;
  using card::launchResolution;
  using card::deepCopy;
} // namespace stormm

#include "card_utilities.tpp"

#endif
