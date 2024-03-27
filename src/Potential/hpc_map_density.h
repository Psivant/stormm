// -*-c++-*-
#ifndef STORMM_HPC_MAP_DENSITY_H
#define STORMM_HPC_MAP_DENSITY_H

#ifdef STORMM_USE_CUDA
#  include <cuda_runtime.h>
#endif
#include "copyright.h"
#include "Constants/behavior.h"

namespace stormm {
namespace energy {

using constants::PrecisionModel;

/// \brief Get the kernel attributes for one of the __shared__ accumulation density mapping
///        kernels, to fill out tables in the core kernel manager.
///
/// \param calc_prec  Indicate whether to carry out calculations in single- or double-precision
/// \param acc_prec   Indicate whether to carry out accumulation in 63- or 95-bit precision
cudaFuncAttributes queryShrAccQMapKernelRequirements(PrecisionModel calc_prec,
                                                     PrecisionModel acc_prec, bool overflow_needed,
                                                     size_t cg_tmat, int order);

/// \brief Get the kernel attributes for one of the general-purpose density mapping kernels, to
///        fill out tables in the core kernel manager.  Descriptions of input parameters follow
///        from queryRegAccQMapKernelRequirements(), above.
cudaFuncAttributes queryGeneralQMapKernelRequirements(PrecisionModel prec, size_t cg_tmat,
                                                      int order);

} // namespace energy
} // namespace stormm

#endif
