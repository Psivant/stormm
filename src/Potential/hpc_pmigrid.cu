// -*-c++-*-
#include "copyright.h"
#include "pmigrid.h"

namespace stormm {
namespace energy {

#include "Numerics/accumulation.cui"

/// \brief Initialize a series of particle-mesh interaction grids by setting all density values to
///        zero.  This is only provided for fixed-precision accumulation on the GPU.  If the grids
///        are stored as real-valued numbers there is no guard against non-associative addition in
///        a parallel computing environment.  It is still possible to work with real-valued grids
///        on the GPU--it is only that one thread (or a fixed series of threads performing a
///        consistent reduction) must write the exact, complete value into the grids.  Such a
///        process would not require initialization of the grid values themselves as the
///        accumulation is done elsewhere.
/// 
/// \param pm_acc  Abstract of the particle-mesh interaction grid object, containing sizing arrays,
///                precision mode, and pointers to the array data
__global__ void __launch_bounds__(large_block_size, 1)
kInitialize(PMIGridAccumulator pm_acc) {
  const uint4 fdims = pm_acc.dims[pm_acc.nsys - 1];
  size_t final_index;
  switch (pm_acc.fftm) {
  case FFTMode::IN_PLACE:
    final_index = fdims.w + (2 * ((fdims.x / 2) + 1) * fdims.y * fdims.z);
    break;
  case FFTMode::OUT_OF_PLACE:
    final_index = fdims.w + (fdims.x * fdims.y * fdims.z);
    break;
  }
  const size_t stride = blockDim.x * gridDim.x;
  switch (pm_acc.mode) {
  case PrecisionModel::DOUBLE:
    for (size_t pos = threadIdx.x + (blockDim.x * blockIdx.x); pos < final_index; pos += stride) {
      pm_acc.lldata[pos] = 0LL;
      if (pm_acc.use_overflow) {
        pm_acc.overflow[pos] = 0;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    for (size_t pos = threadIdx.x + (blockDim.x * blockIdx.x); pos < final_index; pos += stride) {
      pm_acc.idata[pos] = 0;
      if (pm_acc.use_overflow) {
        pm_acc.overflow[pos] = 0;
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMIGridInitialization(PMIGridAccumulator *pm_acc, const GpuDetails &gpu) {
  kInitialize<<<gpu.getSMPCount(), large_block_size>>>(*pm_acc);
}

//-------------------------------------------------------------------------------------------------
extern void launchPMIGridInitialization(PMIGridAccumulator *pm_acc, const int block_count) {
  kInitialize<<<block_count, large_block_size>>>(*pm_acc);
}

/// \brief Convert fixed-precision integer data on a series of particle-mesh interaction grids to
///        real-valued format.
///
/// \param pm_acc  Fixed-precision accumulation abstract for the particle-mesh interaction grids
/// \param pm_wrt  Real-valued writeable abstract for the particle-mesh interaction grids
__global__ void __launch_bounds__(large_block_size, 1)
kConvertToReal(PMIGridWriter pm_wrt, const PMIGridAccumulator pm_acc) {
  const uint4 fdims = pm_acc.dims[pm_acc.nsys - 1];
  size_t final_index;
  switch (pm_acc.fftm) {
  case FFTMode::IN_PLACE:
    final_index = fdims.w + (2 * ((fdims.x / 2) + 1) * fdims.y * fdims.z);
    break;
  case FFTMode::OUT_OF_PLACE:
    final_index = fdims.w + (fdims.x * fdims.y * fdims.z);
    break;
  }
  const size_t stride = blockDim.x * gridDim.x;
  const size_t first_index = threadIdx.x + (blockDim.x * blockIdx.x);
  switch (pm_acc.mode) {
  case PrecisionModel::DOUBLE:
    {
      const double conv_scale = 1.0 / pm_acc.fp_scale;
      if (pm_acc.use_overflow) {
        for (size_t pos = first_index; pos < final_index; pos += stride) {
          pm_wrt.ddata[pos] = int95ToDouble(__ldca(&pm_acc.lldata[pos]),
                                            __ldcv(&pm_acc.overflow[pos])) * conv_scale;
        }
      }
      else {
        for (size_t pos = first_index; pos < final_index; pos += stride) {
          pm_wrt.ddata[pos] = (double)(pm_acc.lldata[pos]) * conv_scale;
        }
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const float conv_scale = (float)(1.0) / pm_acc.fp_scale;
      if (pm_acc.use_overflow) {
        for (size_t pos = first_index; pos < final_index; pos += stride) {
          pm_wrt.fdata[pos] = int63ToFloat(__ldca(&pm_acc.idata[pos]),
                                           __ldcv(&pm_acc.overflow[pos])) * conv_scale;
        }
      }
      else {
        for (size_t pos = first_index; pos < final_index; pos += stride) {
          pm_wrt.fdata[pos] = (float)(pm_acc.idata[pos]) * conv_scale;
        }
      }
    }
    break;
  }
}
  
//-------------------------------------------------------------------------------------------------
extern void launchPMIGridRealConversion(PMIGridWriter *pm_wrt, const PMIGridAccumulator &pm_acc,
                                        const GpuDetails &gpu) {
  kConvertToReal<<<gpu.getSMPCount(), large_block_size>>>(*pm_wrt, pm_acc);
}

//-------------------------------------------------------------------------------------------------
extern void launchPMIGridRealConversion(PMIGridWriter *pm_wrt, const PMIGridAccumulator &pm_acc,
                                        const int block_count) {
  kConvertToReal<<<block_count, large_block_size>>>(*pm_wrt, pm_acc);
}

} // namespace energy
} // namespace stormm
