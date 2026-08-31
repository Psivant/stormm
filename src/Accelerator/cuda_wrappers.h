// -*-c++-*-
#ifndef STORMM_KFUNC_POINTERS_H
#define STORMM_KFUNC_POINTERS_H

#include <string>
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
#    include <cuda_runtime.h>
#    include <cufft.h>
#  endif
#endif
#include "copyright.h"
#include "gpu_enumerators.h"

#ifdef STORMM_USE_HPC

namespace stormm {
namespace card {

#  ifdef STORMM_USE_CUDA
/// \brief A templated function to accept any function pointer and reinterpret it as a constant
///        char pointer for submission to a standard cudaFuncGetAttributes() call.  This will
///        return a CUDA error type that should read "cudaSuccess."
///
/// \param attrib
/// \param 
template <typename T> cudaError_t wrapCudaFuncGetAttributes(cudaFuncAttributes *attrib, T ptr);

/// \brief Produce detailed, comprehensible error messages based on an error code returned by one
///        of the CUDA FFT functions.
///
/// \param event  The cuFFT-related error code
std::string getHpcErrorString(cufftResult event,
                              HpcErrorVerbosity style = HpcErrorVerbosity::REASONING);
#  endif
  
} // namespace card
} // namespace stormm

#  include "cuda_wrappers.tpp"

#endif // STORMM_USE_HPC

#endif
