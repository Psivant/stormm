// -*-c++-*-
#ifndef STORMM_CARD_UTILITIES_H
#define STORMM_CARD_UTILITIES_H

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

} // namespace card
} // namespace stormm

// Include the launch guards in any other STORMM libraries
namespace stormm {
  using card::launchPreparation;
  using card::launchResolution;
} // namespace stormm

#endif
