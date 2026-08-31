// -*-c++-*-
#ifndef STORMM_HPC_GATHER_FORCES_H
#define STORMM_HPC_GATHER_FORCES_H

#ifdef STORMM_USE_CUDA
#  include <cuda_runtime.h>
#endif
#include "copyright.h"
#include "Constants/behavior.h"

namespace stormm {
namespace energy {

/// \brief Get the kernel attributes for one of the general-purpose force gathering kernels, to
///        fill out tables in the core kernel manager.
///
/// \param calc_prec  Indicate whether to carry out calculations in single- or double-precision
/// \param acc_prec   Indicate whether to carry out force accumulation in 63- or 95-bit precision
/// \param order      The order of particle-mesh interpolation
cudaFuncAttributes queryGeneralForceGatheringKernelRequirements(PrecisionModel prec,
                                                                size_t cg_tmat, int order);
  
} // namespace energy
} // namespace stormm

#endif
