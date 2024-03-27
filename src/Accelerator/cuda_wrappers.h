// -*-c++-*-
#ifndef STORMM_KFUNC_POINTERS_H
#define STORMM_KFUNC_POINTERS_H

#include "copyright.h"

#  ifdef STORMM_USE_HPC
#  include <cuda_runtime.h>

namespace stormm {
namespace card {

#    ifdef STORMM_USE_CUDA
/// \brief A templated function to accept any function pointer and reinterpret it as a constant
///        char pointer for submission to a standard cudaFuncGetAttributes() call.  This will
///        return a CUDA error type that should read "cudaSuccess."
///
/// \param attrib
/// \param 
template <typename T> cudaError_t wrapCudaFuncGetAttributes(cudaFuncAttributes *attrib, T ptr);
#    endif

} // namespace card
} // namespace stormm

#  include "cuda_wrappers.tpp"

#  endif
#endif
