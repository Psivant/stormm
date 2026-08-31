#include "copyright.h"
#include "Accelerator/cuda_wrappers.h"
#include "Math/series_ops.h"
#include "Reporting/error_format.h"
#include "fft_stage.h"
#ifdef STORMM_USE_HPC
#  include "hpc_fft_stage.h"
#endif

namespace stormm {
namespace stmath {
  
#ifdef STORMM_USE_HPC
using card::getHpcErrorString;
#endif
using stmath::incrementingSeries;
using namespace pocketfft;
  
//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double* d_signal_in, double2* z_frequency_in, const HybridTargetLevel tier_in,
                   const Normalization normalize_in, const FFTMode mode_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    mode{mode_in}, dimensionality{0}, tier{tier_in}, nx{nx_in}, ny{ny_in}, nz{nz_in}, nw{nw_in},
    problem_size{1}, batch_count{batch_count_in}, prec{PrecisionModel::DOUBLE},
    normalize{normalize_in},
    signal_batch_stride{signal_batch_stride_in},
    frequency_batch_stride{frequency_batch_stride_in},
    external_frequency_data{z_frequency_in != nullptr},
    signal_shape{}, frequency_shape{}, host_signal_shape{},
    host_signal_stride{}, host_frequency_stride{}, axes{},
    c_frequency_data{HybridKind::ARRAY, "cfrequency_fft_data", setFrequencyFormat(tier_in)},
    z_frequency_data{HybridKind::ARRAY, "zfrequency_fft_data", setFrequencyFormat(tier_in)},
    d_signal{d_signal_in}, f_signal{nullptr}, z_signal{nullptr}, c_signal{nullptr},
#ifdef STORMM_USE_HPC
    z_frequency{nullptr}, c_frequency{nullptr}, forward_plan{}, backward_plan{}
#else
    z_frequency{nullptr}, c_frequency{nullptr}
#endif
{
  setDimensions(z_frequency_in, nullptr);
  setFrequencyDataSpace(z_frequency_in);
#ifdef STORMM_USE_HPC
  makePlans();
#endif
}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double2* z_signal_in, double2* z_frequency_in,
                   const HybridTargetLevel tier_in, const Normalization normalize_in,
                   const FFTMode mode_in, const size_t nx_in, const size_t ny_in,
                   const size_t nz_in, const size_t nw_in, const int batch_count_in,
                   const size_t signal_batch_stride_in, const size_t frequency_batch_stride_in) :
    mode{mode_in}, dimensionality{0}, tier{tier_in}, nx{nx_in}, ny{ny_in}, nz{nz_in}, nw{nw_in},
    problem_size{1}, batch_count{batch_count_in}, prec{PrecisionModel::DOUBLE},
    normalize{normalize_in},
    signal_batch_stride{signal_batch_stride_in},
    frequency_batch_stride{frequency_batch_stride_in},
    external_frequency_data{z_frequency_in != nullptr},
    signal_shape{}, frequency_shape{}, host_signal_shape{},
    host_signal_stride{}, host_frequency_stride{}, axes{},
    c_frequency_data{HybridKind::ARRAY, "cfrequency_fft_data", setFrequencyFormat(tier_in)},
    z_frequency_data{HybridKind::ARRAY, "zfrequency_fft_data", setFrequencyFormat(tier_in)},
    d_signal{nullptr}, f_signal{nullptr}, z_signal{z_signal_in}, c_signal{nullptr},
#ifdef STORMM_USE_HPC
    z_frequency{nullptr}, c_frequency{nullptr}, forward_plan{}, backward_plan{}
#else
    z_frequency{nullptr}, c_frequency{nullptr}
#endif
{
  setDimensions(z_frequency_in, nullptr);
  setFrequencyDataSpace(z_frequency_in);
#ifdef STORMM_USE_HPC
  makePlans();
#endif
}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double* d_signal_in, const HybridTargetLevel tier_in,
                   const Normalization normalize_in, const FFTMode mode_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(d_signal_in, nullptr, tier_in, mode_in, nx_in, ny_in, nz_in, nw_in, batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double2* z_signal_in, const HybridTargetLevel tier_in,
                   const Normalization normalize_in, const FFTMode mode_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(z_signal_in, nullptr, tier_in, mode_in, nx_in, ny_in, nz_in, nw_in, batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float* f_signal_in, float2* c_frequency_in, const HybridTargetLevel tier_in,
                   const Normalization normalize_in, const FFTMode mode_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    mode{mode_in}, dimensionality{0}, tier{tier_in}, nx{nx_in}, ny{ny_in}, nz{nz_in}, nw{nw_in},
    problem_size{1}, batch_count{batch_count_in}, prec{PrecisionModel::SINGLE},
    normalize{normalize_in},
    signal_batch_stride{signal_batch_stride_in},
    frequency_batch_stride{frequency_batch_stride_in},
    external_frequency_data{c_frequency_in != nullptr},
    signal_shape{}, frequency_shape{}, host_signal_shape{},
    host_signal_stride{}, host_frequency_stride{}, axes{},
    c_frequency_data{HybridKind::ARRAY, "cfrequency_fft_data", setFrequencyFormat(tier_in)},
    z_frequency_data{HybridKind::ARRAY, "zfrequency_fft_data", setFrequencyFormat(tier_in)},
    d_signal{nullptr}, f_signal{f_signal_in}, z_signal{nullptr}, c_signal{nullptr},
#ifdef STORMM_USE_HPC
    z_frequency{nullptr}, c_frequency{nullptr}, forward_plan{}, backward_plan{}
#else
    z_frequency{nullptr}, c_frequency{nullptr}
#endif
{
  setDimensions(nullptr, c_frequency_in);
  setFrequencyDataSpace(c_frequency_in);
#ifdef STORMM_USE_HPC
  makePlans();
#endif
}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float2* c_signal_in, float2* c_frequency_in, const HybridTargetLevel tier_in,
                   const Normalization normalize_in, const FFTMode mode_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    mode{mode_in}, dimensionality{0}, tier{tier_in}, nx{nx_in}, ny{ny_in}, nz{nz_in}, nw{nw_in},
    problem_size{1}, batch_count{batch_count_in}, prec{PrecisionModel::SINGLE},
    normalize{normalize_in},
    signal_batch_stride{signal_batch_stride_in},
    frequency_batch_stride{frequency_batch_stride_in},
    external_frequency_data{c_frequency_in != nullptr},
    signal_shape{}, frequency_shape{}, host_signal_shape{},
    host_signal_stride{}, host_frequency_stride{}, axes{},
    c_frequency_data{HybridKind::ARRAY, "cfrequency_fft_data", setFrequencyFormat(tier_in)},
    z_frequency_data{HybridKind::ARRAY, "zfrequency_fft_data", setFrequencyFormat(tier_in)},
    d_signal{nullptr}, f_signal{nullptr}, z_signal{nullptr}, c_signal{c_signal_in},
#ifdef STORMM_USE_HPC
    z_frequency{nullptr}, c_frequency{nullptr}, forward_plan{}, backward_plan{}
#else
    z_frequency{nullptr}, c_frequency{nullptr}
#endif
{
  setDimensions(nullptr, c_frequency_in);
  setFrequencyDataSpace(c_frequency_in);
#ifdef STORMM_USE_HPC
  makePlans();
#endif
}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float* f_signal_in, const HybridTargetLevel tier_in,
                   const Normalization normalize_in, const FFTMode mode_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(f_signal_in, nullptr, tier_in, mode_in, nx_in, ny_in, nz_in, nw_in, batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float2* c_signal_in, const HybridTargetLevel tier_in,
                   const Normalization normalize_in, const FFTMode mode_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(c_signal_in, nullptr, tier_in, mode_in, nx_in, ny_in, nz_in, nw_in, batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double* d_signal_in, double2* z_frequency_in, const HybridTargetLevel tier_in,
                   const FFTMode mode_in, const size_t nx_in, const size_t ny_in,
                   const size_t nz_in, const size_t nw_in, const int batch_count_in,
                   const size_t signal_batch_stride_in, const size_t frequency_batch_stride_in) :
    FFTStage(d_signal_in, z_frequency_in, tier_in, Normalization::NO, mode_in, nx_in, ny_in, nz_in,
             nw_in, batch_count_in, signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double2* z_signal_in, double2* z_frequency_in, const HybridTargetLevel tier_in,
                   const FFTMode mode_in, const size_t nx_in, const size_t ny_in,
                   const size_t nz_in, const size_t nw_in, const int batch_count_in,
                   const size_t signal_batch_stride_in, const size_t frequency_batch_stride_in) :
    FFTStage(z_signal_in, z_frequency_in, tier_in, Normalization::NO, mode_in, nx_in, ny_in, nz_in,
             nw_in, batch_count_in, signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double* d_signal_in, const HybridTargetLevel tier_in, const FFTMode mode_in,
                   const size_t nx_in, const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(d_signal_in, nullptr, tier_in, Normalization::NO, mode_in, nx_in, ny_in, nz_in,
             nw_in, batch_count_in, signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double2* z_signal_in, const HybridTargetLevel tier_in, const FFTMode mode_in,
                   const size_t nx_in, const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(z_signal_in, nullptr, tier_in, Normalization::NO, mode_in, nx_in, ny_in, nz_in,
             nw_in, batch_count_in, signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float* f_signal_in, float2* c_frequency_in, const HybridTargetLevel tier_in,
                   const FFTMode mode_in, const size_t nx_in, const size_t ny_in,
                   const size_t nz_in, const size_t nw_in, const int batch_count_in,
                   const size_t signal_batch_stride_in, const size_t frequency_batch_stride_in) :
    FFTStage(f_signal_in, c_frequency_in, tier_in, Normalization::NO, mode_in, nx_in, ny_in, nz_in,
             nw_in, batch_count_in, signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float2* c_signal_in, float2* c_frequency_in, const HybridTargetLevel tier_in,
                   const FFTMode mode_in, const size_t nx_in, const size_t ny_in,
                   const size_t nz_in, const size_t nw_in, const int batch_count_in,
                   const size_t signal_batch_stride_in, const size_t frequency_batch_stride_in) :
    FFTStage(c_signal_in, c_frequency_in, tier_in, Normalization::NO, mode_in, nx_in, ny_in, nz_in,
             nw_in, batch_count_in, signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float* f_signal_in, const HybridTargetLevel tier_in, const FFTMode mode_in,
                   const size_t nx_in, const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(f_signal_in, nullptr, tier_in, Normalization::NO, mode_in, nx_in, ny_in, nz_in,
             nw_in, batch_count_in, signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float2* c_signal_in, const HybridTargetLevel tier_in, const FFTMode mode_in,
                   const size_t nx_in, const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(c_signal_in, nullptr, tier_in, Normalization::NO, mode_in, nx_in, ny_in, nz_in,
             nw_in, batch_count_in, signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double* d_signal_in, std::complex<double>* z_frequency_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(d_signal_in, reinterpret_cast<double2*>(z_frequency_in), HybridTargetLevel::HOST,
             Normalization::NO, FFTMode::OUT_OF_PLACE, nx_in, ny_in, nz_in, nw_in, batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(std::complex<double>* z_signal_in, std::complex<double>* z_frequency_in,
                   const size_t nx_in, const size_t ny_in, const size_t nz_in,
                   const size_t nw_in, const int batch_count_in,
                   const size_t signal_batch_stride_in, const size_t frequency_batch_stride_in) :
    FFTStage(reinterpret_cast<double2*>(z_signal_in), reinterpret_cast<double2*>(z_frequency_in),
             HybridTargetLevel::HOST, Normalization::NO, FFTMode::OUT_OF_PLACE, nx_in, ny_in,
             nz_in, nw_in, signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(std::complex<double>* z_signal_in, const size_t nx_in, const size_t ny_in,
                   const size_t nz_in, const size_t nw_in, const int batch_count_in,
                   const size_t signal_batch_stride_in, const size_t frequency_batch_stride_in) :
    FFTStage(reinterpret_cast<double2*>(z_signal_in), nullptr, HybridTargetLevel::HOST,
             Normalization::NO, FFTMode::OUT_OF_PLACE, nx_in, ny_in, nz_in, nw_in, batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float* f_signal_in, std::complex<float>* c_frequency_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(f_signal_in, reinterpret_cast<float2*>(c_frequency_in), HybridTargetLevel::HOST,
             Normalization::NO, FFTMode::OUT_OF_PLACE, nx_in, ny_in, nz_in, nw_in, batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(std::complex<float>* c_signal_in, std::complex<float>* c_frequency_in,
                   const size_t nx_in, const size_t ny_in, const size_t nz_in,
                   const size_t nw_in, const int batch_count_in,
                   const size_t signal_batch_stride_in, const size_t frequency_batch_stride_in) :
    FFTStage(reinterpret_cast<float2*>(c_signal_in), reinterpret_cast<float2*>(c_frequency_in),
             HybridTargetLevel::HOST, Normalization::NO, FFTMode::OUT_OF_PLACE, nx_in, ny_in,
             nz_in, nw_in,  batch_count_in, signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(std::complex<float>* c_signal_in, const size_t nx_in, const size_t ny_in,
                   const size_t nz_in, const size_t nw_in, const int batch_count_in,
                   const size_t signal_batch_stride_in, const size_t frequency_batch_stride_in) :
    FFTStage(reinterpret_cast<float2*>(c_signal_in), nullptr, HybridTargetLevel::HOST,
             Normalization::NO, FFTMode::OUT_OF_PLACE, nx_in, ny_in, nz_in, nw_in,  batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double* d_signal_in, double2* z_frequency_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(d_signal_in, z_frequency_in, HybridTargetLevel::HOST, Normalization::NO,
             FFTMode::OUT_OF_PLACE, nx_in, ny_in, nz_in, nw_in, batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double2* z_signal_in, double2* z_frequency_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in):
    FFTStage(z_signal_in, z_frequency_in, HybridTargetLevel::HOST, Normalization::NO,
             FFTMode::OUT_OF_PLACE, nx_in, ny_in, nz_in, nw_in, batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double* d_signal_in, const size_t nx_in, const size_t ny_in, const size_t nz_in,
                   const size_t nw_in, const int batch_count_in,
                   const size_t signal_batch_stride_in, const size_t frequency_batch_stride_in):
    FFTStage(d_signal_in, nullptr, HybridTargetLevel::HOST, Normalization::NO,
             FFTMode::OUT_OF_PLACE, nx_in, ny_in, nz_in, nw_in, batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double2* z_signal_in, const size_t nx_in, const size_t ny_in,
                   const size_t nz_in, const size_t nw_in, const int batch_count_in,
                   const size_t signal_batch_stride_in, const size_t frequency_batch_stride_in):
    FFTStage(z_signal_in, nullptr, HybridTargetLevel::HOST, Normalization::NO,
             FFTMode::OUT_OF_PLACE, nx_in, ny_in, nz_in, nw_in, batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float* f_signal_in, float2* c_frequency_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in):
    FFTStage(f_signal_in, c_frequency_in, HybridTargetLevel::HOST, Normalization::NO,
             FFTMode::OUT_OF_PLACE, nx_in, ny_in, nz_in, nw_in, batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float2* c_signal_in, float2* c_frequency_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in):
    FFTStage(c_signal_in, c_frequency_in, HybridTargetLevel::HOST, Normalization::NO,
             FFTMode::OUT_OF_PLACE, nx_in, ny_in, nz_in, nw_in, batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float* f_signal_in, const size_t nx_in, const size_t ny_in, const size_t nz_in,
                   const size_t nw_in, const int batch_count_in,
                   const size_t signal_batch_stride_in, const size_t frequency_batch_stride_in):
    FFTStage(f_signal_in, nullptr, HybridTargetLevel::HOST, Normalization::NO,
             FFTMode::OUT_OF_PLACE, nx_in, ny_in, nz_in, nw_in, batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float2* c_signal_in, const size_t nx_in, const size_t ny_in, const size_t nz_in,
                   const size_t nw_in, const int batch_count_in,
                   const size_t signal_batch_stride_in, const size_t frequency_batch_stride_in):
    FFTStage(c_signal_in, nullptr, HybridTargetLevel::HOST, Normalization::NO,
             FFTMode::OUT_OF_PLACE, nx_in, ny_in, nz_in, nw_in, batch_count_in,
             signal_batch_stride_in, frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double* d_signal_in, double2* z_frequency_in, const HybridTargetLevel tier_in,
                   const size_t nx_in, const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(d_signal_in, z_frequency_in, tier_in, Normalization::NO, FFTMode::OUT_OF_PLACE,
             nx_in, ny_in, nz_in, nw_in, batch_count_in, signal_batch_stride_in,
             frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double2* z_signal_in, double2* z_frequency_in, const HybridTargetLevel tier_in,
                   const size_t nx_in, const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(z_signal_in, z_frequency_in, tier_in, Normalization::NO, FFTMode::OUT_OF_PLACE,
             nx_in, ny_in, nz_in, nw_in, batch_count_in, signal_batch_stride_in,
             frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double* d_signal_in, const HybridTargetLevel tier_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(d_signal_in, nullptr, tier_in, Normalization::NO, FFTMode::OUT_OF_PLACE,
             nx_in, ny_in, nz_in, nw_in, batch_count_in, signal_batch_stride_in,
             frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(double2* z_signal_in, const HybridTargetLevel tier_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(z_signal_in, nullptr, tier_in, Normalization::NO, FFTMode::OUT_OF_PLACE,
             nx_in, ny_in, nz_in, nw_in, batch_count_in, signal_batch_stride_in,
             frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float* f_signal_in, float2* c_frequency_in, const HybridTargetLevel tier_in,
                   const size_t nx_in, const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(f_signal_in, c_frequency_in, tier_in, Normalization::NO, FFTMode::OUT_OF_PLACE, nx_in,
             ny_in, nz_in, nw_in, batch_count_in, signal_batch_stride_in,
             frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float2* c_signal_in, float2* c_frequency_in, const HybridTargetLevel tier_in,
                   const size_t nx_in, const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(c_signal_in, c_frequency_in, tier_in, Normalization::NO, FFTMode::OUT_OF_PLACE,
             nx_in, ny_in, nz_in, nw_in, batch_count_in, signal_batch_stride_in,
             frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float* f_signal_in, const HybridTargetLevel tier_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(f_signal_in, nullptr, tier_in, Normalization::NO, FFTMode::OUT_OF_PLACE,
             nx_in, ny_in, nz_in, nw_in, batch_count_in, signal_batch_stride_in,
             frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::FFTStage(float2* c_signal_in, const HybridTargetLevel tier_in, const size_t nx_in,
                   const size_t ny_in, const size_t nz_in, const size_t nw_in,
                   const int batch_count_in, const size_t signal_batch_stride_in,
                   const size_t frequency_batch_stride_in) :
    FFTStage(c_signal_in, nullptr, tier_in, Normalization::NO, FFTMode::OUT_OF_PLACE,
             nx_in, ny_in, nz_in, nw_in, batch_count_in, signal_batch_stride_in,
             frequency_batch_stride_in)
{}

//-------------------------------------------------------------------------------------------------
FFTStage::~FFTStage() {
  switch (tier) {
  case HybridTargetLevel::HOST:
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
#  ifdef STORMM_USE_CUDA
    cufftDestroy(forward_plan);
    cufftDestroy(backward_plan);
#  endif
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
FFTMode FFTStage::getMode() const {
  return mode;
}

//-------------------------------------------------------------------------------------------------
size_t FFTStage::getProblemSize(const int dimension_id) const {
  validateDimension(dimension_id, "getProblemLength");
  return signal_shape[dimension_id];
}

//-------------------------------------------------------------------------------------------------
size_t FFTStage::getProblemSize() const {
  return problem_size;
}

//-------------------------------------------------------------------------------------------------
int FFTStage::getBatchCount() const {
  return batch_count;
}

//-------------------------------------------------------------------------------------------------
bool FFTStage::hasExternalFrequencyData() const {
  return external_frequency_data;
}

//-------------------------------------------------------------------------------------------------
size_t FFTStage::getFrequencyCount(const int dimension_id) const {
  validateDimension(dimension_id, "getFrequencyCount");
  return frequency_shape[dimension_id];
}

//-------------------------------------------------------------------------------------------------
PrecisionModel FFTStage::getPrecision() const {
  return prec;
}

//-------------------------------------------------------------------------------------------------
HybridTargetLevel FFTStage::getTier() const {
  return tier;
}

//-------------------------------------------------------------------------------------------------
bool FFTStage::normalizedRoundTrip() const {
  switch (normalize) {
  case Normalization::YES:
    return true;
  case Normalization::NO:
    return false;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
size_t FFTStage::getSignalBatchStride() const {
  return signal_batch_stride;
}

//-------------------------------------------------------------------------------------------------
size_t FFTStage::getFrequencyBatchStride() const {
  return frequency_batch_stride;
}

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
//-------------------------------------------------------------------------------------------------
cufftHandle FFTStage::getForwardPlan() const {
  return forward_plan;
}

//-------------------------------------------------------------------------------------------------
cufftHandle FFTStage::getBackwardPlan() const {
  return backward_plan;
}
#  endif
#endif

//-------------------------------------------------------------------------------------------------
double2 FFTStage::getFrequencyDataFromProblem(const size_t problem_id, const size_t x_pos) const {
  if (dimensionality != 1) {
    rtErr("Exactly " + std::to_string(dimensionality) + " indices must be provided for a problem "
          "of this dimension.",	"FFTStage", "getFrequencyData");
  }

  // Apply periodic considerations.  If the FFT is real-to-complex, the following filter will
  // detect a request for information that might be outside the explicit data array.  A
  // complex-to-complex transform will have frequencies for each piece of real data.
  size_t x_act = x_pos - (signal_shape[0] * (x_pos / signal_shape[0]));
  double2 result;
  bool invoke_symmetry;
  if (x_act >= frequency_shape[0]) {
    x_act = signal_shape[0] - x_act;
    invoke_symmetry = true;
  }
  else {
    invoke_symmetry = false;
  }
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      result = z_frequency[x_act];
      break;
    case PrecisionModel::SINGLE:
      {
        const float2 tmp_rslt = c_frequency[x_act + (problem_id * frequency_batch_stride)];
        result = { static_cast<double>(tmp_rslt.x), static_cast<double>(tmp_rslt.y) };
      }
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    result = downloadComplexNumber(x_act + (problem_id * frequency_batch_stride),
                                   "getFrequencyDataFromProblem");
    break;
#endif
  }
  if (invoke_symmetry) {
    result.y = -result.y;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double2 FFTStage::getFrequencyDataFromProblem(const size_t problem_id, const size_t x_pos,
                                              const size_t y_pos) const {
  if (dimensionality != 2) {
    rtErr("Exactly " + std::to_string(dimensionality) + " indices must be provided for a problem "
          "of this dimension.",	"FFTStage", "getFrequencyData");
  }

  // Apply periodic considerations.  If the FFT is real-to-complex, the following filter will
  // detect a request for information that might be outside the explicit data array.  A
  // complex-to-complex transform will have frequencies for each piece of real data.
  size_t x_act = x_pos - (signal_shape[0] * (x_pos / signal_shape[0]));
  size_t y_act = y_pos - (signal_shape[1] * (y_pos / signal_shape[1]));
  bool invoke_symmetry;
  if (x_act >= frequency_shape[1]) {
    x_act = signal_shape[0] - x_act;
    if (y_act > 0) {
      y_act = signal_shape[1] - y_act;
    }
    invoke_symmetry = true;
  }
  else {
    invoke_symmetry = false;
  }
  const size_t xy_act = (y_act * frequency_shape[1]) + x_act +
                        (problem_id * frequency_batch_stride);
  double2 result;
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      result = z_frequency[xy_act];
      break;
    case PrecisionModel::SINGLE:
      {
        const float2 tmp_rslt = c_frequency[xy_act];
        result = { static_cast<double>(tmp_rslt.x), static_cast<double>(tmp_rslt.y) };
      }
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
#  ifdef STORMM_USE_CUDA
    result = downloadComplexNumber(xy_act, "getFrequencyData");
#  endif
    break;
#endif
  }
  if (invoke_symmetry) {
    result.y = -result.y;
  }
  return result;  
}

//-------------------------------------------------------------------------------------------------
double2 FFTStage::getFrequencyDataFromProblem(const size_t problem_id, const size_t x_pos,
                                              const size_t y_pos, const size_t z_pos) const {
  if (dimensionality != 3) {
    rtErr("Exactly " + std::to_string(dimensionality) + " indices must be provided for a problem "
          "of this dimension.",	"FFTStage", "getFrequencyData");
  }

  // Apply periodic considerations.  If the FFT is real-to-complex, the following filter will
  // detect a request for information that might be outside the explicit data array.  A
  // complex-to-complex transform will have frequencies for each piece of real data.
  size_t x_act = x_pos - (signal_shape[0] * (x_pos / signal_shape[0]));
  size_t y_act = y_pos - (signal_shape[1] * (y_pos / signal_shape[1]));
  size_t z_act = z_pos - (signal_shape[2] * (z_pos / signal_shape[2]));
  size_t xyz_act;
  double2 result;
  bool invoke_symmetry;
  if (x_act >= frequency_shape[2]) {
    x_act = signal_shape[0] - x_act;
    if (y_act > 0) {
      y_act = signal_shape[1] - y_act;
    }
    if (z_act > 0) {
      z_act = signal_shape[2] - z_act;
    }
    invoke_symmetry = true;
  }
  else {
    invoke_symmetry = false;
  }
  xyz_act = (((z_act * frequency_shape[1]) + y_act) * frequency_shape[2]) + x_act +
            (problem_id * frequency_batch_stride);
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      result = z_frequency[xyz_act];
      break;
    case PrecisionModel::SINGLE:
      {
        const float2 tmp_rslt = c_frequency[xyz_act];
        result = { static_cast<double>(tmp_rslt.x), static_cast<double>(tmp_rslt.y) };
      }
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
#  ifdef STORMM_USE_CUDA
    result = downloadComplexNumber(xyz_act, "getFrequencyData");
#  endif
    break;
#endif
  }
  if (invoke_symmetry) {
    result.y = -result.y;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double2 FFTStage::getFrequencyDataFromProblem(const size_t problem_id, const size_t x_pos,
                                              const size_t y_pos, const size_t z_pos,
                                              const size_t w_pos) const {
  if (dimensionality != 4) {
    rtErr("Exactly " + std::to_string(dimensionality) + " indices must be provided for a problem "
          "of this dimension.",	"FFTStage", "getFrequencyData");
  }

  // Apply periodic considerations.  If the FFT is real-to-complex, the following filter will
  // detect a request for information that might be outside the explicit data array.  A
  // complex-to-complex transform will have frequencies for each piece of real data.
  size_t x_act = x_pos - (signal_shape[0] * (x_pos / signal_shape[0]));
  size_t y_act = y_pos - (signal_shape[1] * (y_pos / signal_shape[1]));
  size_t z_act = z_pos - (signal_shape[2] * (z_pos / signal_shape[2]));
  size_t w_act = w_pos - (signal_shape[3] * (z_pos / signal_shape[3]));
  size_t xyzw_act;
  double2 result;
  bool invoke_symmetry;
  if (x_act >= frequency_shape[3]) {
    x_act = signal_shape[0] - x_act;
    if (y_act > 0) {
      y_act = signal_shape[1] - y_act;
    }
    if (z_act > 0) {
      z_act = signal_shape[2] - z_act;
    }
    if (w_act > 0) {
      w_act = signal_shape[3] - w_act;
    }
    invoke_symmetry = true;
  }
  else {
    invoke_symmetry = false;
  }
  xyzw_act = (((((w_act * frequency_shape[1]) + z_act) * frequency_shape[2]) + y_act) *
              frequency_shape[3]) + x_act + (problem_id * frequency_batch_stride);
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      result = z_frequency[xyzw_act];
      break;
    case PrecisionModel::SINGLE:
      {
        const float2 tmp_rslt = c_frequency[xyzw_act];
        result = { static_cast<double>(tmp_rslt.x), static_cast<double>(tmp_rslt.y) };
      }
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
#  ifdef STORMM_USE_CUDA
    result = downloadComplexNumber(xyzw_act, "getFrequencyData");
    break;
#  endif
#endif
  }
  if (invoke_symmetry) {
    result.y = -result.y;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double2 FFTStage::getFrequencyData(const size_t x_pos) const {
  return getFrequencyDataFromProblem(0, x_pos);
}

//-------------------------------------------------------------------------------------------------
double2 FFTStage::getFrequencyData(const size_t x_pos, const size_t y_pos) const {
  return getFrequencyDataFromProblem(0, x_pos, y_pos);
}

//-------------------------------------------------------------------------------------------------
double2 FFTStage::getFrequencyData(const size_t x_pos, const size_t y_pos,
                                   const size_t z_pos) const {
  return getFrequencyDataFromProblem(0, x_pos, y_pos, z_pos);
}

//-------------------------------------------------------------------------------------------------
double2 FFTStage::getFrequencyData(const size_t x_pos, const size_t y_pos,
                                   const size_t z_pos, const size_t w_pos) const {
  return getFrequencyDataFromProblem(0, x_pos, y_pos, z_pos, w_pos);
}

//-------------------------------------------------------------------------------------------------
void FFTStage::forwardFFT() {
  switch (tier) {
  case HybridTargetLevel::HOST:
    for (int i = 0; i < batch_count; i++) {
      const size_t zi = i;
      const size_t signal_offset = zi * signal_batch_stride;
      const size_t transform_offset = zi * frequency_batch_stride;
      switch (prec) {
      case PrecisionModel::DOUBLE:
        {
          std::complex<double>* tmp_zfreq = reinterpret_cast<std::complex<double>*>(z_frequency);
          if (d_signal != nullptr) {
            r2c<double>(host_signal_shape, host_signal_stride, host_frequency_stride, axes,
                        FORWARD, &d_signal[signal_offset], &tmp_zfreq[transform_offset], 1.0);
          }
          else {
            std::complex<double>* tmp_zsig = reinterpret_cast<std::complex<double>*>(z_signal);
            c2c<double>(host_signal_shape, host_signal_stride, host_frequency_stride, axes,
                        FORWARD, &tmp_zsig[signal_offset], &tmp_zfreq[transform_offset], 1.0);
          }
        }
        break;
      case PrecisionModel::SINGLE:
        {
          std::complex<float>* tmp_cfreq = reinterpret_cast<std::complex<float>*>(c_frequency);
          if (f_signal != nullptr) {
            r2c<float>(host_signal_shape, host_signal_stride, host_frequency_stride, axes, FORWARD,
                       &f_signal[signal_offset], &tmp_cfreq[transform_offset], 1.0f);
          }
          else {
            std::complex<float>* tmp_csig = reinterpret_cast<std::complex<float>*>(c_signal);
            c2c<float>(host_signal_shape, host_signal_stride, host_frequency_stride, axes, FORWARD,
                       &tmp_csig[signal_offset], &tmp_cfreq[transform_offset], 1.0f);
          }
        }
        break;
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (prec) {
    case PrecisionModel::DOUBLE:
#  ifdef STORMM_USE_CUDA
      if (d_signal != nullptr) {
        const cufftResult event = cufftExecD2Z(forward_plan, d_signal,
                                               reinterpret_cast<cufftDoubleComplex*>(z_frequency));
        if (event != CUFFT_SUCCESS) {
          rtErr("Error in cufftExecD2Z: " + getHpcErrorString(event) + ".", "FFTStage",
                "forwardFFT");
        }
      }
      else {
        const cufftResult event = cufftExecZ2Z(forward_plan,
                                               reinterpret_cast<cufftDoubleComplex*>(z_signal),
                                               reinterpret_cast<cufftDoubleComplex*>(z_frequency),
                                               CUFFT_FORWARD);
        if (event != CUFFT_SUCCESS) {
          rtErr("Error in cufftExecD2Z: " + getHpcErrorString(event) + ".", "FFTStage",
                "forwardFFT");
        }
      }
#  endif
      break;
    case PrecisionModel::SINGLE:
#  ifdef STORMM_USE_CUDA
      if (f_signal != nullptr) {
        const cufftResult event = cufftExecR2C(forward_plan, f_signal,
                                               reinterpret_cast<cufftComplex*>(c_frequency));
        if (event != CUFFT_SUCCESS) {
          rtErr("Error in cufftExecR2C: " + getHpcErrorString(event) + ".", "FFTStage",
                "forwardFFT");
        }
      }
      else {
        const cufftResult event = cufftExecC2C(forward_plan,
                                               reinterpret_cast<cufftComplex*>(c_signal),
                                               reinterpret_cast<cufftComplex*>(c_frequency),
                                               CUFFT_FORWARD);
        if (event != CUFFT_SUCCESS) {
          rtErr("Error in cufftExecC2C: " + getHpcErrorString(event) + ".", "FFTStage",
                "forwardFFT");
        }
      }
#  endif
      break;
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void FFTStage::backwardFFT(const GpuDetails &gpu) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    for (int i = 0; i < batch_count; i++) {
      const size_t zi = i;
      const size_t signal_offset = zi * signal_batch_stride;
      const size_t transform_offset = zi * frequency_batch_stride;
      double norm_factor;
      switch (normalize) {
      case Normalization::YES:
        norm_factor = 1.0 / static_cast<double>(problem_size);
        break;
      case Normalization::NO:
        norm_factor = 1.0;
        break;
      }
      switch (prec) {
      case PrecisionModel::DOUBLE:
        {
          std::complex<double>* tmp_zfreq = reinterpret_cast<std::complex<double>*>(z_frequency);
          if (d_signal != nullptr) {
            c2r<double>(host_signal_shape, host_frequency_stride, host_signal_stride, axes,
                        BACKWARD, &tmp_zfreq[transform_offset], &d_signal[signal_offset],
                        norm_factor);
          }
          else {
            std::complex<double>* tmp_zsig = reinterpret_cast<std::complex<double>*>(z_signal);
            c2c<double>(host_signal_shape, host_frequency_stride, host_signal_stride, axes,
                        BACKWARD, &tmp_zfreq[transform_offset], &tmp_zsig[signal_offset],
                        norm_factor);
          }
        }
        break;
      case PrecisionModel::SINGLE:
        {
          std::complex<float>* tmp_cfreq = reinterpret_cast<std::complex<float>*>(c_frequency);
          if (f_signal != nullptr) {
            c2r<float>(host_signal_shape, host_frequency_stride, host_signal_stride, axes,
                       BACKWARD, &tmp_cfreq[transform_offset], &f_signal[signal_offset],
                       static_cast<float>(norm_factor));
          }
          else {
            std::complex<float>* tmp_csig = reinterpret_cast<std::complex<float>*>(c_signal);
            c2c<float>(host_signal_shape, host_frequency_stride, host_signal_stride, axes,
                       BACKWARD, &tmp_cfreq[transform_offset], &tmp_csig[signal_offset],
                       static_cast<float>(norm_factor));
          }
        }
        break;
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (prec) {
    case PrecisionModel::DOUBLE:
#  ifdef STORMM_USE_CUDA
      if (d_signal != nullptr) {
        cufftExecZ2D(backward_plan, reinterpret_cast<cufftDoubleComplex*>(z_frequency), d_signal);
        switch (normalize) {
        case Normalization::YES:
          normalizeFFT(d_signal, problem_size, batch_count, signal_batch_stride, gpu);
          break;
        case Normalization::NO:
          break;
        }
      }
      else {
        cufftExecZ2Z(backward_plan, reinterpret_cast<cufftDoubleComplex*>(z_frequency),
                     reinterpret_cast<cufftDoubleComplex*>(z_signal), CUFFT_INVERSE);
        switch (normalize) {
        case Normalization::YES:
          normalizeFFT(z_signal, problem_size, batch_count, signal_batch_stride, gpu);
          break;
        case Normalization::NO:
          break;
        }
      }
#  endif
      break;
    case PrecisionModel::SINGLE:
#  ifdef STORMM_USE_CUDA
      if (f_signal != nullptr) {
        cufftExecC2R(backward_plan, reinterpret_cast<cufftComplex*>(c_frequency), f_signal);
        switch (normalize) {
        case Normalization::YES:
          normalizeFFT(f_signal, problem_size, batch_count, signal_batch_stride, gpu);
          break;
        case Normalization::NO:
          break;
        }
      }
      else {
        cufftComplex* tmp_csig = reinterpret_cast<cufftComplex*>(z_signal);
        cufftExecC2C(backward_plan, reinterpret_cast<cufftComplex*>(c_frequency),
                     reinterpret_cast<cufftComplex*>(c_signal), CUFFT_INVERSE);
        switch (normalize) {
        case Normalization::YES:
          normalizeFFT(c_signal, problem_size, batch_count, signal_batch_stride, gpu);
          break;
        case Normalization::NO:
          break;
        }
      }
#  endif
      break;
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void FFTStage::validateProblemIndex(const size_t problem_id, const char* caller) const {
  if (problem_id > 0 && problem_id >= batch_count) {
    rtErr("Problem index " + std::to_string(problem_id) + " is invalid for a set of " +
          std::to_string(batch_count) + " problems.", "FFTStage", caller);
  }
}

//-------------------------------------------------------------------------------------------------
HybridFormat FFTStage::setFrequencyFormat(const HybridTargetLevel tier_in) {
  switch (tier_in) {
  case HybridTargetLevel::HOST:
    return HybridFormat::HOST_ONLY;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return HybridFormat::DEVICE_ONLY;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void FFTStage::validateDimension(const int dimension_id, const char* caller) const {
  if (dimension_id < 0 || dimension_id >= dimensionality) {
    rtErr("Dimension index " + std::to_string(dimension_id) + " is invalid for a problem with " +
          std::to_string(dimensionality) + " dimensions.", "FFTStage", caller);
  }
}

//-------------------------------------------------------------------------------------------------
void FFTStage::checkSignalFormat(const size_t ct) const {

  // Check the base type for complex-valued data.
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  if ((ct == stdz_type_index || ct == double2_type_index || ct == cufftz_type_index) &&
      z_signal == nullptr) {
#  endif
#else
  if ((ct == stdz_type_index || ct == double2_type_index) && z_signal == nullptr) {
#endif
    rtErr("A request for complex, double-precision data was submitted to a problem not based on "
          "the same.", "FFTStage", "checkSignalFormat");
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  }
#  endif
#else
  }
#endif
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  else if ((ct == stdc_type_index || ct == float2_type_index  || ct == cufftc_type_index) &&
           c_signal == nullptr) {
#  endif
#else
  else if ((ct == stdc_type_index || ct == float2_type_index) && c_signal == nullptr) {
#endif
    rtErr("A request for complex, single-precision data was submitted to a problem not based on "
          "the same.", "FFTStage", "checkSignalFormat");
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  }
#  endif
#else
  }
#endif
  else if (ct == double_type_index && d_signal == nullptr) {
    rtErr("A request for real, double-precision data was submitted to a problem not based on the "
          "same.", "FFTStage", "checkSignalFormat");
  }
  else if (ct == float_type_index && f_signal == nullptr) {
    rtErr("A request for real, single-precision data was submitted to a problem not based on the "
          "same.", "FFTStage", "checkSignalFormat");
  }
}

//-------------------------------------------------------------------------------------------------
void FFTStage::setDimensions(const double2* z_frequency_in, const float2* c_frequency_in) {

  // Check that the dimensions are listed in sensible order.
  if (nx <= 0 || (nz > 0 && ny <= 0) || (nw > 0 && (ny <= 0 || nz <= 0))) {
    std::string dim_str = std::string("[ ") + std::to_string(nx);
    int max_dim = 1;
    if (ny >= 0) max_dim = 2;
    if (nz >= 0) max_dim = 3;
    if (nw >= 0) max_dim = 4;
    if (max_dim >= 2) {
      dim_str += std::string(" ") + std::to_string(ny);
    }
    if (max_dim >= 3) {
      dim_str += std::string(" ") + std::to_string(nz);
    }
    if (max_dim == 4) {
      dim_str += std::string(" ") + std::to_string(nw);
    }
    dim_str += std::string(" ]");
    rtErr("Dimensions " + dim_str + " are invalid.", "FFTStage", "setDimensions");
  }
  dimensionality = (nx > 0) + (ny > 0) + (nz > 0) + (nw > 0);

  // For real-to-complex transforms (or complex-to real inverse transforms), which are expected to
  // be performed if the class object is initialized based on real data, out-of-place transforms
  // are all that is permitted.
  switch (mode) {
  case FFTMode::IN_PLACE:
    switch (tier) {
    case HybridTargetLevel::HOST:
      rtErr("In-place transforms are not permitted in CPU-based operations.", "FFTStage",
            "setDimensions");
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
  case FFTMode::OUT_OF_PLACE:
    break;
  }
  
  // Set the shape and strides for C++ computations.
  axes = incrementingSeries<size_t>(0, dimensionality);
  signal_shape.resize(dimensionality);
  frequency_shape.resize(dimensionality);
  host_signal_shape.resize(dimensionality);
  signal_shape[0] = nx;
  if (dimensionality > 1) signal_shape[1] = ny;
  if (dimensionality > 2) signal_shape[2] = nz;
  if (dimensionality > 3) signal_shape[3] = nw;
  for (int i = 0; i < dimensionality; i++) {
    problem_size *= signal_shape[i];
  }
  size_t tmp_orig, tmp_freq;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    tmp_orig = (d_signal != nullptr) ? sizeof(double) : sizeof(double2);
    tmp_freq = sizeof(double2);
    break;
  case PrecisionModel::SINGLE:
    tmp_orig = (f_signal != nullptr) ? sizeof(float) : sizeof(float2);
    tmp_freq = sizeof(float2);
    break;
  }
  host_signal_stride.resize(dimensionality);
  host_frequency_stride.resize(dimensionality);
  for (int i = 0; i < dimensionality; i++) {
    host_signal_shape[i] = signal_shape[dimensionality - 1 - i];
    frequency_shape[i] = (i == dimensionality - 1) ? (host_signal_shape[i] / 2) + 1 :
                                                     host_signal_shape[i];
  }
  for (int i = dimensionality - 1; i >= 0; i--) {
    host_signal_stride[i] = tmp_orig;
    tmp_orig *= host_signal_shape[i];
    host_frequency_stride[i] = tmp_freq;
    tmp_freq *= frequency_shape[i];
  }
  
  // Check that the strides are reasonable, if there is a non-unitary batch count.
  if (batch_count > 1) {
    size_t shape_vol = 1;
    size_t fshape_vol = 1;
    for (int i = 0; i < dimensionality; i++) {
      shape_vol *= signal_shape[i];
      fshape_vol *= frequency_shape[i];
    }

    // If the frequency data space is not yet set, set the frequency data stride based on the
    // total frequency volume.
    if (c_frequency_in == nullptr && z_frequency_in == nullptr) {
      frequency_batch_stride = fshape_vol;
    }
    if (signal_batch_stride < shape_vol || frequency_batch_stride < fshape_vol) {
      std::string dim_str = std::string("[ ") + std::to_string(nx);
      int max_dim = 1;
      if (dimensionality >= 2) {
        dim_str += std::string(" ") + std::to_string(ny);
      }
      if (dimensionality >= 3) {
        dim_str += std::string(" ") + std::to_string(nz);
      }
      if (dimensionality == 4) {
        dim_str += std::string(" ") + std::to_string(nw);
      }
      dim_str += std::string(" ]");
      if (signal_batch_stride < shape_vol) {
        if (signal_batch_stride == 0) {
          rtErr("A batch stride must be provided for a non-unitary batch count (" +
                std::to_string(batch_count) + ").", "FFTStage", "setDimensions");
        }
        else {
          rtErr("A signal data batch stride of " + std::to_string(signal_batch_stride) +
                " is insufficient for problems of size " + dim_str + ".", "FFTStage",
                "setDimensions");
        }
      }
      if (frequency_batch_stride < fshape_vol) {
        if (frequency_batch_stride == 0) {
          rtErr("A frequency space batch stride must be provided for a non-unitary batch count (" +
                std::to_string(batch_count) + ").", "FFTStage", "setDimensions");
        }
        else {
          rtErr("A frequency data batch stride of " + std::to_string(frequency_batch_stride) +
                " is insufficient for problems of size " + dim_str + ".", "FFTStage",
                "setDimensions");
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void FFTStage::setFrequencyDataSpace(double2* frequency_in) {
  if (frequency_in == nullptr) {
    size_t freq_size = 1;
    for (int i = 0; i < dimensionality; i++) {
      freq_size *= frequency_shape[i];
    }
    if (batch_count > 1) {
      freq_size *= frequency_batch_stride;
    }
    z_frequency_data.resize(freq_size);
    z_frequency = z_frequency_data.data(tier);
  }
  else {
    z_frequency = frequency_in;
  }
}
  
//-------------------------------------------------------------------------------------------------
void FFTStage::setFrequencyDataSpace(float2* frequency_in) {
  if (frequency_in == nullptr) {
    size_t freq_size = 1;
    for (int i = 0; i < dimensionality; i++) {
      freq_size *= frequency_shape[i];
    }
    if (batch_count > 1) {
      freq_size *= frequency_batch_stride;
    }
    c_frequency_data.resize(freq_size);
    c_frequency = c_frequency_data.data(tier);
  }
  else {
    c_frequency = frequency_in;
  }
}

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
//-------------------------------------------------------------------------------------------------
void FFTStage::makePlans() {

  // Set the operation types based on the nature of the signal data (which implies the nature of
  // the frequency data).
  cufftType fwd_op, bkwd_op;
  if (d_signal != nullptr) {
    fwd_op  = CUFFT_D2Z;
    bkwd_op = CUFFT_Z2D;
  }
  else if (f_signal != nullptr) {
    fwd_op  = CUFFT_R2C;
    bkwd_op = CUFFT_C2R;
  }
  else if (z_signal != nullptr) {
    fwd_op  = CUFFT_Z2Z;
    bkwd_op  = CUFFT_Z2Z;
  }
  else if (c_signal != nullptr) {
    fwd_op  = CUFFT_C2C;
    bkwd_op  = CUFFT_C2C;
  }

  // Make plans for single transformations or multiple transforms, in one or more dimensions.
  switch (tier) {
  case HybridTargetLevel::HOST:
    break;
  case HybridTargetLevel::DEVICE:
    {
      cufftResult event;
      bool forward_problem = false;
      bool backward_problem = false;
      if (dimensionality == 1 && batch_count == 1) {
        event = cufftPlan1d(&forward_plan, nx, fwd_op, 1);
        if (event != CUFFT_SUCCESS) {
          forward_problem = true;
        }
        event = cufftPlan1d(&backward_plan, nx, bkwd_op, 1);
        if (event != CUFFT_SUCCESS) {
          backward_problem = true;
        }
      }
      else {
        std::vector<int> cufft_shape(dimensionality);

        // NVIDIA's cuFFT expects multi-dimensional arrays to be laid out in "C" order (row-major).
        // Feeding the dimensions to the program in reverse order will resolve the conflict with
        // Fortran (column-major) order.
        for (int i = 0; i < dimensionality; i++) {
          cufft_shape[i] = signal_shape[dimensionality - 1 - i];
        }
        event = cufftPlanMany(&forward_plan, dimensionality, cufft_shape.data(), nullptr, 1,
                              signal_batch_stride, nullptr, 1, frequency_batch_stride, fwd_op,
                              batch_count);
        if (event != CUFFT_SUCCESS) {
          forward_problem = true;
        }
        event = cufftPlanMany(&backward_plan, dimensionality, cufft_shape.data(), nullptr, 1,
                              frequency_batch_stride, nullptr, 1, signal_batch_stride, bkwd_op,
                              batch_count);
        if (event != CUFFT_SUCCESS) {
          backward_problem = true;
        }
      }
      if (forward_problem || backward_problem) {
        std::string len_string = (dimensionality > 1) ? std::string("s { ") : std::string(" { ");
        for (int i = 0; i < dimensionality; i++) {
          len_string += std::to_string(signal_shape[i]);
          if (i < dimensionality - 1) {
            len_string += std::string(", ");
          }
        }
        len_string += std::string(" },");
        std::string direction;
        if (forward_problem && backward_problem) {
          direction = std::string("both forward and backward plans");
        }
        else if (forward_problem) {
          direction = std::string("the forward plan");
        }
        else if (backward_problem) {
          direction = std::string("the backward plan");
        }
        rtErr("Creation of " + direction + " was unsuccessful for a problem of dimension " +
              std::to_string(dimensionality) + ", length" + len_string + ".  Reason: " +
              getHpcErrorString(event) + ".", "FFTStage", "makePlans");
      }
    }
    break;
  }
}
#  endif

//-------------------------------------------------------------------------------------------------
double2 FFTStage::downloadComplexNumber(const size_t pos, const char* caller) const {
  double result[2];
  switch (prec) {
  case PrecisionModel::DOUBLE:
#  ifdef STORMM_USE_CUDA
    if (cudaMemcpy(&result, z_frequency + pos, 2 * sizeof(double), cudaMemcpyDeviceToHost) !=
        cudaSuccess) {
      rtErr("Failure in cudaMemcpy for " + getEnumerationName(prec) + "-precision complex "
            "element " + std::to_string(pos) + ".", "FFTStage", caller);
    }
#  endif
    break;
  case PrecisionModel::SINGLE:
    {
      float result_f[2];
#  ifdef STORMM_USE_CUDA
      if (cudaMemcpy(&result_f, c_frequency + pos, 2 * sizeof(float), cudaMemcpyDeviceToHost) !=
          cudaSuccess) {
        rtErr("Failure in cudaMemcpy for " + getEnumerationName(prec) + "-precision complex "
              "element " + std::to_string(pos) + ".", "FFTStage", "getFrequencyData");
      }
#  endif
      result[0] = result_f[0];
      result[1] = result_f[1];
    }
    break;
  }
  return { result[0], result[1] };
}
#endif

} // namespace stmath
} // namespace stormm
