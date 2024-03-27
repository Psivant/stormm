// -*-c++-*-
#ifndef STORMM_GPU_ENNUMERATORS_H
#define STORMM_GPU_ENNUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace card {

/// \brief Indication of whether a Hybrid object is an actual array or just contains pointers to
///        another Hybrid object that actually is an array.
enum class HybridKind {
  ARRAY,   ///< The Hybrid object has its own memory tied to its pointers.
  POINTER  ///< The Hybrid object does not own any memory and instead targets other ARRAY-kind
           ///<   Hybrid objects.
};

/// \brief Available formats for Hybrid data arrays.
enum class HybridFormat {
#ifdef STORMM_USE_HPC
  EXPEDITED,   ///< The default, best general-purpose format for side-by-side data on the host and
               ///<   device: pinned (page-locked) memory on the host and a separate array on the
               ///<   GPU.  40% higher upload and download bandwidth than DECOUPLED memory with the
               ///<   only constraint being that host memory is page-locked, unswappable to disk.
  DECOUPLED,   ///< Pageable memory on the host and a separate array on the device.
  UNIFIED,     ///< Managed memory, one array that the developer sees which is automatically
               ///<   updated based on changes at either the host or device levels.
  HOST_ONLY,   ///< Page-locked memory allocated by cudaHostAlloc(), accessible to the GPU by a
               ///<   device pointer and transfer through the PCIE bus
  DEVICE_ONLY  ///< A cudaMalloc() allocated array on the device with no equivalent on the host.
               ///<   accessible to the host only through cudaMemcpy or one of the wrappers like
               ///<   readHost() and putHost().
#else
  HOST_ONLY    ///< A new[] allocated array when CUDA or HIP is not used in the compilation
#endif
};

/// \brief Two levels to which the pointers might go, depending on whether an HPC language is
///        compiled
enum class HybridTargetLevel {
#ifdef STORMM_USE_HPC
  HOST, DEVICE
#else
  HOST
#endif
};

/// \brief List the different protocols for synchronizing the device and the host when launching
///        HPC kernels.  The CPU will launch a kernel and then continue working through its own
///        code, and can even launch a large number of kernels back-to-back through iterations of
///        a loop.  If Kernels are all launched in the same stream and do not require CPU
///        communication, there is generally no need to synchronize.  However, if the Kernels
///        depend on one another and launch in separate streams, it can be essential to have the
///        host wait until certain kernels are finished before launching a new kernel.  More often,
///        a kernel will read or write data on host resources which are critical to the CPU
///        process, and this requires synchronization as well.
enum class HpcKernelSync {
  BEFORE,            ///< Wait until other kernels are finished before launching a new kernel
  AFTER,             ///< Launch the kernel and then have the CPU wait until it finishes
  BEFORE_AND_AFTER,  ///< Synchronize before and after launching the kernel
  NO_SYNC,           ///< Synchronization is declared to be unnecessary
  MEMORY_AUTO        ///< Triggers common behavior relating to memory tranfers between the host and
                     ///<   device.  When combined with instructions about the origin and
                     ///<   destination of the memory, this will trigger synchronization at the
                     ///<   start of a kernel launcher for transferring memory from the host to the
                     ///<   device, at the end of a kernel launcher for transferring memory from
                     ///<   the device to the host, and no synchronization in the case that the
                     ///<   memory transfers from host to host or device to device.
};

/// \brief Produce human-readable strings corresponding to each enumerated value.  Overloads of
///        this function here and in other libraries are provided for each enumerator.
///
/// \param input  The enumerated value of interest
/// \{
std::string getEnumerationName(HybridKind input);
std::string getEnumerationName(HybridFormat input);
std::string getEnumerationName(HybridTargetLevel input);
std::string getEnumerationName(HpcKernelSync input);
/// \}
  
} // namespace card
} // namespace stormm

#endif
