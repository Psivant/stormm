// -*-c++-*-
#ifndef STORMM_FFT_STAGE_H
#define STORMM_FFT_STAGE_H

#include <complex>
#include "copyright.h"
#include "pocketfft_hdronly.h"
#include "Accelerator/gpu_enumerators.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "DataTypes/stormm_vector_types.h"
#include "math_enumerators.h"

namespace stormm {
namespace stmath {

using card::GpuDetails;
using card::Hybrid;
using card::HybridFormat;
using card::HybridKind;
using card::HybridTargetLevel;
using constants::PrecisionModel;

/// \brief Prepare the plans, (optional) out-of-place memory needed for CPU-based FFTs, and
///        encapsulate the processes for performing Fast Fourier Transforms.
class FFTStage {
public:

  /// \brief The constructor will prepare stages for FFTs on the CPU, GPU, or both.  Dimensions of
  ///        the problem are required.  Up to four dimensions will be considered.  Overloads are
  ///        provided to automatically set the object to work in single- or double-precision mode,
  ///        performing real-to-complex or complex-to-complex transformations, based on the format
  ///        of the original (signal) data.
  /// \{
  FFTStage(double* d_signal_in, double2* z_frequency_in, HybridTargetLevel tier_in,
           Normalization normalize_in, FFTMode mode_in, size_t nx_in, size_t ny_in = 0,
           size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(double2* z_signal_in, double2* z_frequency_in, HybridTargetLevel tier_in,
           Normalization normalize_in, FFTMode mode_in, size_t nx_in, size_t ny_in = 0,
           size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(double* d_signal_in, HybridTargetLevel tier_in, Normalization normalize_in,
           FFTMode mode_in, size_t nx_in, size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0,
           int batch_count_in = 1, size_t signal_batch_stride_in = 0,
           size_t frequency_batch_stride_in = 0);

  FFTStage(double2* z_signal_in, HybridTargetLevel tier_in, Normalization normalize_in,
           FFTMode mode_in, size_t nx_in, size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0,
           int batch_count_in = 1, size_t signal_batch_stride_in = 0,
           size_t frequency_batch_stride_in = 0);

  FFTStage(float* f_signal_in, float2* c_frequency_in, HybridTargetLevel tier_in,
           Normalization normalize_in, FFTMode mode_in, size_t nx_in, size_t ny_in = 0,
           size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(float2* c_signal_in, float2* c_frequency_in, HybridTargetLevel tier_in,
           Normalization normalize_in, FFTMode mode_in, size_t nx_in, size_t ny_in = 0,
           size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(float* f_signal_in, HybridTargetLevel tier_in, Normalization normalize_in,
           FFTMode mode_in, size_t nx_in, size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0,
           int batch_count_in = 1, size_t signal_batch_stride_in = 0,
           size_t frequency_batch_stride_in = 0);

  FFTStage(float2* c_signal_in, HybridTargetLevel tier_in, Normalization normalize_in,
           FFTMode mode_in, size_t nx_in, size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0,
           int batch_count_in = 1, size_t signal_batch_stride_in = 0,
           size_t frequency_batch_stride_in = 0);
  
  FFTStage(double* d_signal_in, double2* z_frequency_in, HybridTargetLevel tier_in,
           FFTMode mode_in, size_t nx_in, size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0,
           int batch_count_in = 1, size_t signal_batch_stride_in = 0,
           size_t frequency_batch_stride_in = 0);

  FFTStage(double2* z_signal_in, double2* z_frequency_in, HybridTargetLevel tier_in,
           FFTMode mode_in, size_t nx_in, size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0,
           int batch_count_in = 1, size_t signal_batch_stride_in = 0,
           size_t frequency_batch_stride_in = 0);

  FFTStage(double* d_signal_in, HybridTargetLevel tier_in, FFTMode mode_in, size_t nx_in,
           size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(double2* z_signal_in, HybridTargetLevel tier_in, FFTMode mode_in, size_t nx_in,
           size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(float* f_signal_in, float2* c_frequency_in, HybridTargetLevel tier_in, FFTMode mode_in,
           size_t nx_in, size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0,
           int batch_count_in = 1, size_t signal_batch_stride_in = 0,
           size_t frequency_batch_stride_in = 0);

  FFTStage(float2* c_signal_in, float2* c_frequency_in, HybridTargetLevel tier_in, FFTMode mode_in,
           size_t nx_in, size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0,
           int batch_count_in = 1, size_t signal_batch_stride_in = 0,
           size_t frequency_batch_stride_in = 0);

  FFTStage(float* f_signal_in, HybridTargetLevel tier_in, FFTMode mode_in, size_t nx_in,
           size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(float2* c_signal_in, HybridTargetLevel tier_in, FFTMode mode_in, size_t nx_in,
           size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);
  
  FFTStage(double* d_signal_in, std::complex<double>* z_frequency_in, size_t nx_in,
           size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(std::complex<double>* z_signal_in, std::complex<double>* z_frequency_in, size_t nx_in,
           size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(std::complex<double>* z_signal_in, size_t nx_in, size_t ny_in = 0, size_t nz_in = 0,
           size_t nw_in = 0, int batch_count_in = 1, size_t signal_batch_stride_in = 0,
           size_t frequency_batch_stride_in = 0);

  FFTStage(float* f_signal_in, std::complex<float>* c_frequency_in, size_t nx_in,
           size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(std::complex<float>* c_signal_in, std::complex<float>* c_frequency_in, size_t nx_in,
           size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(std::complex<float>* c_signal_in, size_t nx_in, size_t ny_in = 0, size_t nz_in = 0,
           size_t nw_in = 0, int batch_count_in = 1, size_t signal_batch_stride_in = 0,
           size_t frequency_batch_stride_in = 0);
  
  FFTStage(double* d_signal_in, double2* z_frequency_in, size_t nx_in, size_t ny_in = 0,
           size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(double2* z_signal_in, double2* z_frequency_in, size_t nx_in, size_t ny_in = 0,
           size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(double* d_signal_in, size_t nx_in, size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0,
           int batch_count_in = 1, size_t signal_batch_stride_in = 0,
           size_t frequency_batch_stride_in = 0);

  FFTStage(double2* z_signal_in, size_t nx_in, size_t ny_in = 0, size_t nz_in = 0,
           size_t nw_in = 0, int batch_count_in = 1, size_t signal_batch_stride_in = 0,
           size_t frequency_batch_stride_in = 0);

  FFTStage(float* f_signal_in, float2* c_frequency_in, size_t nx_in, size_t ny_in = 0,
           size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(float2* c_signal_in, float2* c_frequency_in, size_t nx_in, size_t ny_in = 0,
           size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(float* f_signal_in, size_t nx_in, size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0,
           int batch_count_in = 1, size_t signal_batch_stride_in = 0,
           size_t frequency_batch_stride_in = 0);

  FFTStage(float2* c_signal_in, size_t nx_in, size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0,
           int batch_count_in = 1, size_t signal_batch_stride_in = 0,
           size_t frequency_batch_stride_in = 0);
  
  FFTStage(double* d_signal_in, double2* z_frequency_in, HybridTargetLevel tier_in, size_t nx_in,
           size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(double2* z_signal_in, double2* z_frequency_in, HybridTargetLevel tier_in,
           size_t nx_in, size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0,
           int batch_count_in = 1, size_t signal_batch_stride_in = 0,
           size_t frequency_batch_stride_in = 0);

  FFTStage(double* d_signal_in, HybridTargetLevel tier_in, size_t nx_in, size_t ny_in = 0,
           size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(double2* z_signal_in, HybridTargetLevel tier_in, size_t nx_in, size_t ny_in = 0,
           size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(float* f_signal_in, float2* c_frequency_in, HybridTargetLevel tier_in, size_t nx_in,
           size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(float2* c_signal_in, float2* c_frequency_in, HybridTargetLevel tier_in, size_t nx_in,
           size_t ny_in = 0, size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(float* f_signal_in, HybridTargetLevel tier_in, size_t nx_in, size_t ny_in = 0,
           size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);

  FFTStage(float2* c_signal_in, HybridTargetLevel tier_in, size_t nx_in, size_t ny_in = 0,
           size_t nz_in = 0, size_t nw_in = 0, int batch_count_in = 1,
           size_t signal_batch_stride_in = 0, size_t frequency_batch_stride_in = 0);
  /// \}

  /// \brief With no const elements, the default copy and move constructors, as well as copy and
  ///        move assignment operators, remain valid.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right hand side of an assignment statement
  /// \{
  FFTStage(const FFTStage &original) = default;
  FFTStage(FFTStage &&original) = default;
  FFTStage& operator=(const FFTStage &original) = default;
  FFTStage& operator=(FFTStage &&original) = default;
  /// \}
  
  /// \brief The destructor must be written out in order to handle destruction of the FFT plans in
  ///        specific cases.
  ~FFTStage();

  /// \brief Get the mode in which FFTs are to run.
  FFTMode getMode() const;
  
  /// \brief Get the length of the problem in one dimension.
  ///
  /// Overloaded:
  ///   - Get the size along a particular axis
  ///   - Get the total volume of the problem, the product of sizes along all axes
  ///
  /// \param dimension_id  The dimension of interest.  This will be checked for validity.
  /// \{
  size_t getProblemSize(int dimension_id) const;
  size_t getProblemSize() const;
  /// \}

  /// \brief Get the number of problems in the managed batch.
  int getBatchCount() const;

  /// \brief Query whether the object references external data to store the frequency output of
  ///        its transforms.
  bool hasExternalFrequencyData() const;
  
  /// \brief Get the length of the problem in one dimension.
  ///
  /// \param dimension_id  The dimension of interest.  This will be checked for validity.
  size_t getFrequencyCount(int dimension_id) const;

  /// \brief Get the mode in which the object is operating.
  PrecisionModel getPrecision() const;

  /// \brief Get an indication of whether FFTs are performed on the CPU host or the GPU device.
  HybridTargetLevel getTier() const;

  /// \brief Get an indication of whether the data is (re)normalized after a backward / inverse
  ///        FFT.
  bool normalizedRoundTrip() const;

  /// \brief Get the signal batch stride, the distance between separate signals in the data.  This
  ///        class expects that individual, distinct signals are contiguous in memory although
  ///        there could be padding between successive signals.  The stride will be stated in terms
  ///        of the number of elements of the signal data type.
  size_t getSignalBatchStride() const;

  /// \brief Get the frequency batch stride, the distance between separate frequency projections of
  ///        each signal.  This stride will be set automatically if new data is allocated by the
  ///        object itself to hold the complex result.  The stride will be given in terms of the
  ///        number of elements of the complex frequency data type.
  size_t getFrequencyBatchStride() const;

  /// \brief Get a value of the original data from a specific problem within a batch of many
  ///        similar transforms.  This is the most general case of getSignalData(), below.
  ///
  /// Overloaded:
  ///   - Produce signal data for problems of one, two, three, or four dimensions
  ///
  /// \param problem_id  Index of the problem of interest
  /// \param x_pos       Index along the array, perhaps its first dimension, at which to extract
  ///                    data
  /// \param y_pos       Index along the second dimension
  /// \param z_pos       Index along the third dimension
  /// \param w_pos       Index along the fourth dimension
  /// \{
  template <typename T> T getSignalDataFromProblem(size_t problem_id, size_t x_pos) const;
  template <typename T> T getSignalDataFromProblem(size_t problem_id, size_t x_pos,
                                                   size_t y_pos) const;
  template <typename T> T getSignalDataFromProblem(size_t problem_id, size_t x_pos, size_t y_pos,
                                                   size_t z_pos) const;
  template <typename T> T getSignalDataFromProblem(size_t problem_id, size_t x_pos, size_t y_pos,
                                                   size_t z_pos, size_t w_pos) const;
  /// \}
  
  /// \brief Get a pointer to or value of the original data, from the first problem of the batch
  ///        if the FFT involves multiple problems.  When getting a pointer to the full data, the  
  ///        developer may specify multiple formats for the return type, but each will be checked
  ///        against the object's actual contents.  An exception will be raised in the event of an
  ///        invalid return type.  This function applies only to the first 
  ///
  /// Overloaded:
  ///   - Produce a const-qualified pointer from a const object.
  ///   - Produce a pointer to mutable data from a non-const object.
  ///   - Produce the value of the complex data at a specific element.  The number of provided
  ///     indices and positions along each of them will be checked for validity.  If the transform
  ///     is real-to-complex and the requested index is beyond the stored data, the return value
  ///     will be calculated by symmetry.
  ///
  /// \param x_pos  Index along the array, perhaps its first dimension, at which to extract data.
  /// \param y_pos  Index along the second dimension
  /// \param z_pos  Index along the third dimension
  /// \param w_pos  Index along the fourth dimension
  /// \{
  template <typename T> const T* getSignalData() const;
  template <typename T> T* getSignalData();
  template <typename T> T getSignalData(size_t x_pos) const;
  template <typename T> T getSignalData(size_t x_pos, size_t y_pos) const;
  template <typename T> T getSignalData(size_t x_pos, size_t y_pos, size_t z_pos) const;
  template <typename T> T getSignalData(size_t x_pos, size_t y_pos, size_t z_pos,
                                        size_t w_pos) const;
  /// \}

#ifdef STORMM_USE_HPC
  /// \brief Get the handle to the plan for performing forward transforms.
#  ifdef STORMM_USE_CUDA
  cufftHandle getForwardPlan() const;
#  endif
  /// \brief Get the handle to the plan for performing backward transforms.
#  ifdef STORMM_USE_CUDA
  cufftHandle getBackwardPlan() const;
#  endif
#endif
  /// \brief Get a pointer to or value of the frequency-space data for a specific problem out of a
  ///        batch.  Overloading and descriptions of input parameters follow from
  ///        getSignalDataFromProblem(), above.
  /// \{
  double2 getFrequencyDataFromProblem(size_t problem_id, size_t x_pos) const;
  double2 getFrequencyDataFromProblem(size_t problem_id, size_t x_pos, size_t y_pos) const;
  double2 getFrequencyDataFromProblem(size_t problem_id, size_t x_pos, size_t y_pos,
                                      size_t z_pos) const;
  double2 getFrequencyDataFromProblem(size_t problem_id, size_t x_pos, size_t y_pos, size_t z_pos,
                                      size_t w_pos) const;
  /// \}

  /// \brief Get a pointer to or value of the frequency-space data.  Overloading and descriptions
  ///        of input parameters follow from getSignalData(), above.
  /// \{
  template <typename T> const T* getFrequencyData() const;
  template <typename T> T* getFrequencyData();
  double2 getFrequencyData(size_t x_pos) const;
  double2 getFrequencyData(size_t x_pos, size_t y_pos) const;
  double2 getFrequencyData(size_t x_pos, size_t y_pos, size_t z_pos) const;
  double2 getFrequencyData(size_t x_pos, size_t y_pos, size_t z_pos, size_t w_pos) const;
  /// \}
  
  /// \brief Carry out a forward FFT, transforming the signal data into frequencies.
  void forwardFFT();

  /// \brief Carry out a backward FFT, transforming the frequency data back into the same space as
  ///        the signal data.
  ///
  /// \param gpu  Details of the GPU that should be used to carry out a post-hoc normalization of
  ///             data on the DEVICE memory.  If the FFTStage object is not tasked with
  ///             normalization, this parameter is ignored.  No GPU details need be passed to the
  ///             the object constructor, as cuFFT plans are created irrespective of such details.
  void backwardFFT(const GpuDetails &gpu = null_gpu);

private:

  FFTMode mode;             ///< Critical parameter governing whether the FFTs will be done "in
                            ///<   place" or "out of place." With PocketFFT on the CPU, "out of
                            ///<   place" is the only viable option.
  int dimensionality;       ///< Number of dimensions in the problem, obtained by inspecting each
                            ///<   of the lengths
  HybridTargetLevel tier;   ///< Indicator of where the object is prepared to carry out FFTs: on
                            ///<   the CPU host, on the GPU device, or both
  size_t nx;                ///< Number of array elements in the first (of up to four) dimensions
  size_t ny;                ///< Number of array elements in the second (of up to four) dimensions.
                            ///<   A value of zero means that there is but one dimension.
  size_t nz;                ///< Number of array elements in the third (of up to four) dimensions.
                            ///<   A value of zero means that there are at most two dimensions.
  size_t nw;                ///< Number of array elements in the fourth (of up to four) dimensions.
  size_t problem_size;      ///< Total number of points in the grid for one problem, convenient for
                            ///<   normalizations
  int batch_count;          ///< The number of simultaneous FFTs to perform in this set, each
                            ///<   problem arranged with a consistent stride in the same array
  PrecisionModel prec;      ///< The precision of calculations and an indication of which array
                            ///<   pointers held by the object will point to meaningful data.  This
                            ///<   is stored for ease of access.
  Normalization normalize;  ///< Task the manager with normalizing the results of the backward FFT.
                            ///<   While useful in calculations of a convolution, the normalization
                            ///<   performed in this way is not efficient for GPU operations as it
                            ///<   triggers an extra memory "round trip" to cap the signal >>
                            ///<   frequency decomposition >> signal round trip and is easily
                            ///<   incorporated by hand into other operations.  The normalization
                            ///<   option is intended for prototyping GPU applications.  The effect
                            ///<   is equivalent but not as pronounced with CPU work.  Default
                            ///<   "NO."

  /// The strides between data sets for each problem in a batch.  Default zero.  If the batch count
  /// is greater than 1, a non-zero stride must also be provided in the call to the constructor.
  /// \{
  size_t signal_batch_stride;
  size_t frequency_batch_stride;
  /// \}

  /// Indicator of whether external data arrays were provided with pre-allocated space to hold
  /// frequency-space data, the result of a forward transform operation on the signal data (which
  /// must be pre-allocated).
  bool external_frequency_data;

  /// The shape of the problem, incorporating all lengths provided to the constructor as the
  /// "signal" shape and then computing the appropriate shape of the problem in frequency space
  /// \{
  std::vector<size_t> signal_shape;
  std::vector<size_t> frequency_shape;
  /// \}

  /// The shape of the problem, as adapted for evaluation the host with PocketFFT.  The critical
  /// difference is that the sequence of lengths is reversed to evaluate the Fortran-ordered
  /// arrays.
  /// \{
  std::vector<size_t> host_signal_shape;
  /// \}
  
  /// The strides to be taken if PocketFFT is used on the C++ layer
  /// \{
  std::vector<ptrdiff_t> host_signal_stride;
  std::vector<ptrdiff_t> host_frequency_stride;
  /// \}

  /// The axis indexing for PocketFFT
  std::vector<size_t> axes;
  
  /// Dedicated data arrays for transformations into and out of frequency space, allocated only if
  /// external sources are unavailable.
  /// \{
  Hybrid<double2> z_frequency_data;
  Hybrid<float2> c_frequency_data;
  /// \}

  /// Pointers external arrays containing the signal data.  The first letter of each member
  /// variable follows a convention in many FFT packages, f = float (real), d = double (real),
  /// c = float2 (complex), z = double2 (complex).
  /// \{
  double* d_signal;
  float* f_signal;
  double2* z_signal;
  float2* c_signal;
  /// \}

  /// Pointers to external arrays containing frequency-space data.  These arrays are optional.
  /// \{
  double2* z_frequency;
  float2* c_frequency;
  /// \}

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  cufftHandle forward_plan;   ///< Plan to carry out the transformation into frequency space
  cufftHandle backward_plan;  ///< Plan to carry out the transformation back to original space
#  endif
#endif

  /// \brief Validate a requested problem index against the actual number of problems in the FFT
  ///        batch.
  ///
  /// \param problem_id  Index of the problem of interest
  /// \param caller      Name of the calling function, for backtracing purposes
  void validateProblemIndex(const size_t problem_id, const char* caller) const;
  
  /// \brief Set the frequency format.  This is used during initialization.
  ///
  /// \param tier_in  Specify whether transforms are to be performed on the CPU host or GPU device
  HybridFormat setFrequencyFormat(HybridTargetLevel tier_in);

  /// \brief Validate a selection of one of the problem's dimensions.
  ///
  /// \param  dimension_id  The dimension of interest
  /// \param  caller        Name of the calling function (for error tracing purposes)
  void validateDimension(int dimension_id, const char* caller) const;

  /// \brief Check the data format of the problem against the format of a request for its signal
  ///        data.  This will guard against requests for real data from problems with a complex
  ///        signal and vice-versa.
  ///
  /// \param ct  Codified type format of the return value.  The problem's signal data is to be
  ///            rendered in this format (if reasonable).
  void checkSignalFormat(size_t ct) const;
  
  /// \brief Handle the size settings for the object, regardless of the precision model.
  ///
  /// \param z_frequency_in  Input from the constructor, indicating the double-precision frequency
  ///                        space data pointer.  If null, no double-precision data is available
  ///                        outside of what the object may allocate for itself.  This is provided
  ///                        so that the frequency stride may be set automatically, if needed.
  /// \param c_frequency_in  Input from the constructor, indicating the single-precision frequency
  ///                        space data pointer
  void setDimensions(const double2* z_frequency_in, const float2* c_frequency_in);

  /// \brief Set the frequency space data pointer, after allocating private data for the object,
  ///        if needed.  The frequency_shape array member variable must be established prior to
  ///        calling this function.
  ///
  /// Overloaded:
  ///   - Allocate, if needed, and set the pointer for double-precision complex data
  ///   - Allocate, if needed, and set the pointer for single-precision complex data
  ///
  /// \param frequency_in   A pointer to pre-allocated memory for the complex data filled by the
  ///                       frequency space representation of the signal data.  If this is
  ///                       submitted as nullptr it is assumed that no allocation was available.
  /// \{
  void setFrequencyDataSpace(double2* frequency_in);
  void setFrequencyDataSpace(float2* frequency_in);
  /// \}

  /// \brief Lay out the plans for FFTs in the HPC layer.  This will wrap the APIs for FFTs in
  ///        high-performance computing packages.
  void makePlans();

  /// \brief Download a specifc complex number from the frequency data of the object's managed
  ///        problem.
  ///
  /// \param pos     Index of the element to download
  /// \param caller  Name of the calling function
  double2 downloadComplexNumber(const size_t pos, const char* caller) const;
};

} // namespace stmath
} // namespace stormm

#include "fft_stage.tpp"

#endif
