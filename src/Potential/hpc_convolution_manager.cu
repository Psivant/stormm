// -*-c++-*-
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "Constants/symbol_values.h"
#include "convolution_manager.h"

namespace stormm {
namespace energy {

using symbols::pi;
using symbols::pi_f;

#include "Accelerator/ptx_macros.h"
#include "Math/rounding.cui"
  
#define TCALC double
#define TCALC2 double2
#define EXP_FUNC exp
#  define KERNEL_NAME kdPMEGreensFunction
#    include "pme_greens_function.cui"
#  undef KERNEL_NAME
#  define COMPUTE_ENERGY
#    define KERNEL_NAME kdePMEGreensFunction
#      include "pme_greens_function.cui"
#    undef KERNEL_NAME
#  undef COMPUTE_ENERGY
#undef EXP_FUNC
#undef TCALC2
#undef TCALC

#define TCALC float
#define TCALC2 float2
#define TCALC_IS_SINGLE
#define EXP_FUNC expf
#  define KERNEL_NAME kfPMEGreensFunction
#    include "pme_greens_function.cui"
#  undef KERNEL_NAME
#  define COMPUTE_ENERGY
#    define KERNEL_NAME kfePMEGreensFunction
#      include "pme_greens_function.cui"
#    undef KERNEL_NAME
#  undef COMPUTE_ENERGY
#undef EXP_FUNC
#undef TCALC_IS_SINGLE
#undef TCALC2
#undef TCALC
//-------------------------------------------------------------------------------------------------
void pmeGreensFunction(ConvolutionWriter<double, double2> *cvolw, const PsSynthesisBorders &pssb,
                       const PMIGridReader &pmigr, const GpuDetails &gpu, ScoreCardWriter *scw) {
  if (scw == nullptr) {
    kdPMEGreensFunction<<<4 * gpu.getSMPCount(), small_block_size>>>(*cvolw, pssb, pmigr);
  }
  else {
    kdePMEGreensFunction<<<4 * gpu.getSMPCount(), small_block_size>>>(*cvolw, pssb, pmigr, *scw);
  }
}

//-------------------------------------------------------------------------------------------------
void pmeGreensFunction(ConvolutionWriter<float, float2> *cvolw, const PsSynthesisBorders &pssb,
                       const PMIGridReader &pmigr, const GpuDetails &gpu, ScoreCardWriter *scw) {
  if (scw == nullptr) {
    kfPMEGreensFunction<<<4 * gpu.getSMPCount(), small_block_size>>>(*cvolw, pssb, pmigr);
  }
  else {
    kfePMEGreensFunction<<<4 * gpu.getSMPCount(), small_block_size>>>(*cvolw, pssb, pmigr, *scw);
  }
}

//-------------------------------------------------------------------------------------------------
void applyConvolution(ConvolutionWriter<double, double2> *cvolw,
                      const PsSynthesisBorders &pssb, const PMIGridReader &pmigr,
                      const GpuDetails &gpu, ScoreCardWriter *scw) {
  const size_t n_fft = cvolw->fft_ops->size();
  if (n_fft == 0) {
    rtErr("No FFT groups were found in the convolution(s) for " +
          std::to_string(cvolw->system_count) + " systems.\n", "applyConvolution");
  }
  for (size_t i = 0; i < n_fft; i++) {
    cvolw->fft_ops->at(i).forwardFFT();
  }
  if (scw == nullptr) {
    kdPMEGreensFunction<<<4 * gpu.getSMPCount(), small_block_size>>>(*cvolw, pssb, pmigr);
  }
  else {
    kdePMEGreensFunction<<<4 * gpu.getSMPCount(), small_block_size>>>(*cvolw, pssb, pmigr, *scw);
  }
  for (size_t i = 0; i < n_fft; i++) {
    cvolw->fft_ops->at(i).backwardFFT();
  }
}

//-------------------------------------------------------------------------------------------------
void applyConvolution(ConvolutionWriter<float, float2> *cvolw,
                      const PsSynthesisBorders &pssb, const PMIGridReader &pmigr,
                      const GpuDetails &gpu, ScoreCardWriter *scw) {
  const size_t n_fft = cvolw->fft_ops->size();
  if (n_fft == 0) {
    rtErr("No FFT groups were found in the convolution(s) for " +
          std::to_string(cvolw->system_count) + " systems.\n", "applyConvolution");
  }
  for (size_t i = 0; i < n_fft; i++) {
    cvolw->fft_ops->at(i).forwardFFT();
  }
  if (scw == nullptr) {
    kfPMEGreensFunction<<<4 * gpu.getSMPCount(), small_block_size>>>(*cvolw, pssb, pmigr);
  }
  else {
    kfePMEGreensFunction<<<4 * gpu.getSMPCount(), small_block_size>>>(*cvolw, pssb, pmigr, *scw);
  }
  for (size_t i = 0; i < n_fft; i++) {
    cvolw->fft_ops->at(i).backwardFFT();
  }
}

//-------------------------------------------------------------------------------------------------
void applyConvolution(ConvolutionManager *cvol, const PhaseSpaceSynthesis *poly_ps,
                      const PMIGrid *pmig, const GpuDetails &gpu, ScoreCard *sc) {
  const HybridTargetLevel op_tier = cvol->getOperatingTier();
  const PsSynthesisBorders pssb = poly_ps->borders(op_tier);
  const PMIGridReader pmigr = pmig->data(op_tier);
  switch (op_tier) {
  case HybridTargetLevel::HOST:

    // A GPU-enabled program will be able to execute convolutions using CPU resources if the
    // ConvolutionManager is configured in that way, but a vector of energies output by the
    // function designed for operations on the CPU will not be available if called via this route.
    applyConvolution(cvol, poly_ps, pmig, sc);
    break;
  case HybridTargetLevel::DEVICE:

    // The kernel is launched according to the number of streaming multiprocessors on the GPU
    cvol->forwardFFT();
    switch (pmig->getMode()) {
    case PrecisionModel::DOUBLE:
      {
        ConvolutionWriter<double, double2> cvolw = cvol->dpData();
        if (sc == nullptr) {
          kdPMEGreensFunction<<<4 * gpu.getSMPCount(), small_block_size>>>(cvolw, pssb, pmigr);
        }
        else {
          ScoreCardWriter scw = sc->data(HybridTargetLevel::DEVICE);
          kdePMEGreensFunction<<<4 * gpu.getSMPCount(), small_block_size>>>(cvolw, pssb, pmigr,
                                                                            scw);
        }
      }
      break;
    case PrecisionModel::SINGLE:
      {
        ConvolutionWriter<float, float2> cvolw = cvol->spData();
        if (sc == nullptr) {
          kfPMEGreensFunction<<<4 * gpu.getSMPCount(), small_block_size>>>(cvolw, pssb, pmigr);
        }
        else {
          ScoreCardWriter scw = sc->data(HybridTargetLevel::DEVICE);
          kfePMEGreensFunction<<<4 * gpu.getSMPCount(), small_block_size>>>(cvolw, pssb, pmigr,
                                                                            scw);
        }
      }
      break;
    }
    cvol->backwardFFT();
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void applyConvolution(ConvolutionManager *cvol, const PhaseSpaceSynthesis &poly_ps,
                      const PMIGrid &pmig, const GpuDetails &gpu, ScoreCard *sc) {
  applyConvolution(cvol, poly_ps.getSelfPointer(), pmig.getSelfPointer(), gpu, sc);
}

} // namespace energy
} // namespace stormm
