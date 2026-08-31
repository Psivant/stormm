// -*-c++-*-
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "hpc_fft_stage.h"

namespace stormm {
namespace stmath {

/// \brief Perform a normalization on a data array.  These kernels pull data from main memory and
///        perform a single multiplication operation on it before putting the data right back.
///        Hence, they are inefficient and should not be used beyond prototyping applications.
///
/// Overloaded:
///   - Normalize double-precision real data
///   - Normalize single-precision real data
///   - Normalize double-precision complex data
///   - Normalize single-precision complex data
///
/// \param data      The data to normalize
/// \param n_points  The number of points in the data array.  For multi-dimensional arrays, this
///                  is the product of the lengths along all dimensions.
/// \{
__global__ void __launch_bounds__(large_block_size, 1)
kNormalizeFFT(double* data, const size_t n_points, const size_t n_batch,
              const size_t batch_stride) {
  const double dn = 1.0 / (double)(n_points);
  const size_t grid_stride = gridDim.x * blockDim.x;
  const size_t total_work = n_points * n_batch;
  for (size_t i = (blockIdx.x * blockDim.x) + threadIdx.x; i < total_work; i += grid_stride) {
    const size_t problem_id = i / n_points;
    const size_t sub_idx = i - (problem_id * n_points);
    data[(problem_id * batch_stride) + sub_idx] *= dn;
  }
}

__global__ void __launch_bounds__(large_block_size, 1)
  kNormalizeFFT(float* data, const size_t n_points, const size_t n_batch,
                const size_t batch_stride) {
  const float dn = 1.0 / (float)(n_points);
  const size_t grid_stride = gridDim.x * blockDim.x;
  const size_t total_work = n_points * n_batch;
  for (size_t i = (blockIdx.x * blockDim.x) + threadIdx.x; i < total_work; i += grid_stride) {
    const size_t problem_id = i / n_points;
    const size_t sub_idx = i - (problem_id * n_points);
    data[(problem_id * batch_stride) + sub_idx] *= dn;
  }
}

__global__ void __launch_bounds__(large_block_size, 1)
kNormalizeFFT(double2* data, const size_t n_points, const size_t n_batch,
              const size_t batch_stride) {
  const double dn = 1.0 / (double)(n_points);
  const size_t grid_stride = gridDim.x * blockDim.x;
  const size_t total_work = n_points * n_batch;
  for (size_t i = (blockIdx.x * blockDim.x) + threadIdx.x; i < total_work; i += grid_stride) {
    const size_t problem_id = i / n_points;
    const size_t sub_idx = i - (problem_id * n_points);
    double2 tval = data[(problem_id * batch_stride) + sub_idx];
    tval.x *= dn;
    tval.y *= dn;
    data[(problem_id * batch_stride) + sub_idx] = tval;
  }
}

__global__ void __launch_bounds__(large_block_size, 1)
kNormalizeFFT(float2* data, const size_t n_points, const size_t n_batch,
              const size_t batch_stride) {
  const float dn = 1.0 / (float)(n_points);
  const size_t grid_stride = gridDim.x * blockDim.x;
  const size_t total_work = n_points * n_batch;
  for (size_t i = (blockIdx.x * blockDim.x) + threadIdx.x; i < total_work; i += grid_stride) {
    const size_t problem_id = i / n_points;
    const size_t sub_idx = i - (problem_id * n_points);
    float2 tval = data[(problem_id * batch_stride) + sub_idx];
    tval.x *= dn;
    tval.y *= dn;
    data[(problem_id * batch_stride) + sub_idx] = tval;
  }
}
/// \}

//-------------------------------------------------------------------------------------------------
extern void normalizeFFT(double* data, const size_t n_points, const size_t n_batch,
                         const size_t batch_stride, const GpuDetails &gpu) {
  kNormalizeFFT<<<gpu.getSMPCount(), large_block_size>>>(data, n_points, n_batch, batch_stride);
}

//-------------------------------------------------------------------------------------------------
extern void normalizeFFT(float* data, const size_t n_points, const size_t n_batch,
                         const size_t batch_stride, const GpuDetails &gpu) {
  kNormalizeFFT<<<gpu.getSMPCount(), large_block_size>>>(data, n_points, n_batch, batch_stride);
}

//-------------------------------------------------------------------------------------------------
extern void normalizeFFT(double2* data, const size_t n_points, const size_t n_batch,
                         const size_t batch_stride, const GpuDetails &gpu) {
  kNormalizeFFT<<<gpu.getSMPCount(), large_block_size>>>(data, n_points, n_batch, batch_stride);
}

//-------------------------------------------------------------------------------------------------
extern void normalizeFFT(float2* data, const size_t n_points, const size_t n_batch,
                         const size_t batch_stride, const GpuDetails &gpu) {
  kNormalizeFFT<<<gpu.getSMPCount(), large_block_size>>>(data, n_points, n_batch, batch_stride);
}

} // namespace stmath
} // namespace stormm
