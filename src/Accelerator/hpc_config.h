// -*-c++-*-
#ifndef STORMM_HPC_CONFIG_H
#define STORMM_HPC_CONFIG_H

#ifdef STORMM_USE_CUDA
#include <cusolverDn.h>
#include "cublas_v2.h"
#endif
#include "copyright.h"
#include "Constants/behavior.h"
#include "gpu_details.h"

namespace stormm {
namespace card {

using constants::ExceptionResponse;

/// \brief Take an image of the available system hardware.  One or more GPUs can be assigned to a
///        thread based on this information.
class HpcConfig {
public:
  
  /// \brief Constructor for an HpcConfig object.  One such object should be present in any given
  ///        STORMM executable.
  HpcConfig(ExceptionResponse policy = ExceptionResponse::DIE);

  /// \brief Destructor encapsulates HPC shutdown protocols
  ~HpcConfig();

  /// \brief Return the total number of GPUs in the server or workstation, whether they are
  ///        supported by STORMM or not, whether they are available or not
  int getOverallGpuCount() const;

  /// \brief Return the count of available and supported GPUs in the server or workstation
  int getAvailableGpuCount() const;

  /// \brief Return the count of supported and supported GPUs in the server or workstation.  The
  ///        available GPUs are a subset of the supported GPUs.
  int getSupportedGpuCount() const;

  /// \brief Return information on a particular GPU in the server or workstation
  ///
  /// \param gpu_index  Index of the GPU of interest
  GpuDetails getGpuInfo(int gpu_index) const;

  /// \brief Return the indices and specs of one or more GPU devices
  ///
  /// \param requested_count  The number of available GPUs sought for this program's runtime
  std::vector<int> getGpuDevice(int requested_count) const;

#ifdef STORMM_USE_HPC
  /// \brief Return the cuBLAS handle, stored for the lifetime of this HpcConfig object
  cublasHandle_t getCuBlasHandle() const;

  /// \brief Return the cuSolver handle, stored for the lifetime of this HpcConfig object
  cusolverDnHandle_t getCuSolverHandle() const;
#endif
  
private:
  int overall_gpu_count;              ///< The physical number of GPUs detected in the server
  int available_gpu_count;            ///< The number of available GPUs
  int supported_gpu_count;            ///< The number of supported GPUs  
  std::vector<GpuDetails> gpu_list;   ///< Details an availability of each GPU in the system
#ifdef STORMM_USE_HPC
  cublasHandle_t cublas_handle;       ///< The cuBLAS handle, initialized during construction
  cusolverDnHandle_t cusolver_handle; ///< The cuSolver handle, initialized during construction
#endif
};

/// \brief Return a list of all viable GPUs which will support STORMM's code and are not occupied
///        with some other activity.
///
/// \param policy  The behavior to take if no viable GPUs are found
const std::vector<GpuDetails> queryGpuStats(ExceptionResponse policy = ExceptionResponse::DIE);

} // namespace card
} // namespace stormm

#endif
