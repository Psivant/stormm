// -*-c++-*-
#include <cuda_runtime.h>
#include "Accelerator/gpu_details.h"
#include "DataTypes/stormm_vector_types.h"
#include "copyright.h"

namespace stormm {
namespace stmath {

using card::GpuDetails;
  
/// \brief Perform a normalization on a data array, dividng each of its values by the total number
///        of points.  Each overloaded variant of this functions launches one of the corresponding
///        kernels, all of which are noted to be inefficient outside of prototyping.
///
/// Overloaded:
///   - Normalize double-precision real data
///   - Normalize single-precision real data
///   - Normalize double-precision complex data
///   - Normalize single-precision complex data
///
/// \param data          The data to normalize
/// \param n_points      The number of points in the data array.  For multi-dimensional arrays,
///                      this is the product of the lengths along all dimensions.
/// \param n_batch       The number of problems in the batch
/// \param batch_stride  The stride between signals in the FFT batch
/// \{
void normalizeFFT(double* data, size_t n_points, size_t n_batch, size_t batch_stride,
                  const GpuDetails &gpu);
void normalizeFFT(float* data, size_t n_points, size_t n_batch, size_t batch_stride,
                  const GpuDetails &gpu);
void normalizeFFT(double2* data, size_t n_points, size_t n_batch, size_t batch_stride,
                  const GpuDetails &gpu);
void normalizeFFT(float2* data, size_t n_points, size_t n_batch, size_t batch_stride,
                  const GpuDetails &gpu);
/// \}

} // namespace stmath
} // namespace stormm
