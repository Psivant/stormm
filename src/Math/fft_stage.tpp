// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename T> const T* FFTStage::getSignalData() const {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  checkSignalFormat(ct);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    if (d_signal != nullptr) {
      return reinterpret_cast<const T*>(d_signal);
    }
    else {
      return reinterpret_cast<const T*>(z_signal);
    }
    break;
  case PrecisionModel::SINGLE:
    if (f_signal != nullptr) {
      return reinterpret_cast<const T*>(f_signal);
    }
    else {
      return reinterpret_cast<const T*>(c_signal);
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T* FFTStage::getSignalData() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  checkSignalFormat(ct);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    if (d_signal != nullptr) {
      return reinterpret_cast<T*>(d_signal);
    }
    else {
      return reinterpret_cast<T*>(z_signal);
    }
    break;
  case PrecisionModel::SINGLE:
    if (f_signal != nullptr) {
      return reinterpret_cast<T*>(f_signal);
    }
    else {
      return reinterpret_cast<T*>(c_signal);
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T FFTStage::getSignalDataFromProblem(const size_t problem_id,
                                                           const size_t x_pos) const {
  validateProblemIndex(problem_id, "getSignalDataFromProblem");
  if (dimensionality != 1) {
    rtErr("Exactly " + std::to_string(dimensionality) + " indices must be provided for a problem "
          "of this dimension.", "FFTStage", "checkSignalFormat");
  }
  const T* tmp_ptr = getSignalData<T>();
  const size_t x_act = x_pos + (problem_id * signal_batch_stride);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return tmp_ptr[x_act];
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      std::vector<T> tmp_val(1);
      if (cudaMemcpy(tmp_val.data(), &tmp_ptr[x_act], sizeof(T), cudaMemcpyDeviceToHost) !=
          cudaSuccess) {
        rtErr("Failure in cudaMemcpy for position " + std::to_string(x_act) + ".  This pertains "
              "to problem " + std::to_string(problem_id) + " out of a batch of " +
              std::to_string(batch_count) + " problems, position " + std::to_string(x_pos) + ".",
              "FFTStage", "getSignalDataFromProblem");
      }
      return tmp_val[0];
    }
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T FFTStage::getSignalDataFromProblem(const size_t problem_id,
                                                           const size_t x_pos,
                                                           const size_t y_pos) const {
  validateProblemIndex(problem_id, "getSignalDataFromProblem");
  if (dimensionality != 2) {
    rtErr("Exactly " + std::to_string(dimensionality) + " indices must be provided for a problem "
          "of this dimension.", "FFTStage", "checkSignalFormat");
  }
  const T* tmp_ptr = getSignalData<T>();
  const size_t xy_act = (y_pos * signal_shape[0]) + x_pos + (problem_id * signal_batch_stride);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return tmp_ptr[xy_act];
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      std::vector<T> tmp_val(1);
      if (cudaMemcpy(tmp_val.data(), &tmp_ptr[xy_act], sizeof(T), cudaMemcpyDeviceToHost) !=
          cudaSuccess) {
        rtErr("Failure in cudaMemcpy for position " + std::to_string(xy_act) + ".  This pertains "
              "to problem " + std::to_string(problem_id) + " out of a batch of " +
              std::to_string(batch_count) + " problems, position (" + std::to_string(x_pos) +
              ", " + std::to_string(y_pos) + ").", "FFTStage", "getSignalDataFromProblem");
      }
      return tmp_val[0];
    }
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T FFTStage::getSignalDataFromProblem(const size_t problem_id,
                                                           const size_t x_pos, const size_t y_pos,
                                                           const size_t z_pos) const {
  validateProblemIndex(problem_id, "getSignalDataFromProblem");
  if (dimensionality != 3) {
    rtErr("Exactly " + std::to_string(dimensionality) + " indices must be provided for a problem "
          "of this dimension.", "FFTStage", "checkSignalFormat");
  }
  const T* tmp_ptr = getSignalData<T>();
  const size_t xyz_act = (((z_pos * signal_shape[1]) + y_pos) * signal_shape[0]) + x_pos +
                         (problem_id * signal_batch_stride);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return tmp_ptr[xyz_act];
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      std::vector<T> tmp_val(1);
      if (cudaMemcpy(tmp_val.data(), &tmp_ptr[xyz_act], sizeof(T), cudaMemcpyDeviceToHost) !=
          cudaSuccess) {
        rtErr("Failure in cudaMemcpy for position " + std::to_string(xyz_act) + ".  This pertains "
              "to problem " + std::to_string(problem_id) + " out of a batch of " +
              std::to_string(batch_count) + " problems, position (" + std::to_string(x_pos) +
              ", " + std::to_string(y_pos) + ", " + std::to_string(z_pos) + ").", "FFTStage",
              "getSignalDataFromProblem");
      }
      return tmp_val[0];
    }
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T FFTStage::getSignalDataFromProblem(const size_t problem_id,
                                                           const size_t x_pos, const size_t y_pos,
                                                      	   const size_t z_pos,
                                                           const size_t w_pos) const {
  validateProblemIndex(problem_id, "getSignalDataFromProblem");
  if (dimensionality != 4) {
    rtErr("Exactly " + std::to_string(dimensionality) + " indices must be provided for a problem "
          "of this dimension.", "FFTStage", "checkSignalFormat");
  }
  const T* tmp_ptr = getSignalData<T>();
  return tmp_ptr[(((((w_pos * signal_shape[2]) + z_pos) * signal_shape[1]) + y_pos) *
                  signal_shape[0]) + x_pos + (problem_id * signal_batch_stride)];
}

//-------------------------------------------------------------------------------------------------
template <typename T> T FFTStage::getSignalData(const size_t x_pos) const {
  return getSignalDataFromProblem<T>(0, x_pos);
}

//-------------------------------------------------------------------------------------------------
template <typename T> T FFTStage::getSignalData(const size_t x_pos, const size_t y_pos) const {
  return getSignalDataFromProblem<T>(0, x_pos, y_pos);
}

//-------------------------------------------------------------------------------------------------
template <typename T> T FFTStage::getSignalData(const size_t x_pos, const size_t y_pos,
                                                const size_t z_pos) const {
  return getSignalDataFromProblem<T>(0, x_pos, y_pos, z_pos);
}

//-------------------------------------------------------------------------------------------------
template <typename T> T FFTStage::getSignalData(const size_t x_pos, const size_t y_pos,
                                                const size_t z_pos, const size_t w_pos) const {
  return getSignalDataFromProblem<T>(0, x_pos, y_pos, z_pos, w_pos);
}

//-------------------------------------------------------------------------------------------------
template <typename T> const T* FFTStage::getFrequencyData() const {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  switch (prec) {
  case PrecisionModel::DOUBLE:
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
    if (ct == stdz_type_index || ct == double2_type_index || ct == cufftz_type_index) {
#  endif
#else
    if (ct == stdz_type_index || ct == double2_type_index) {
#endif
      return reinterpret_cast<T*>(z_frequency);
#ifdef STORMM_USE_HPC
    }
#else
    }
#endif
    else {
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
      if (ct == stdc_type_index || ct == float2_type_index || ct == cufftc_type_index) {
#  endif
#else
      if (ct == stdc_type_index || ct == float2_type_index) {
#endif
        rtErr("A request for a single-precision complex data pointer was sent to an object "
              "holding double-precision complex data.", "FFTStage", "getFrequencyData");
#ifdef STORMM_USE_HPC
      }
#else
      }
#endif
      else {
        rtErr("A request for an invalid complex data type was received.", "FFTStage",
              "getFrequencyData");
      }
    }
    break;
  case PrecisionModel::SINGLE:
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
    if (ct == stdc_type_index || ct == float2_type_index || ct == cufftc_type_index) {
#  endif
#else
    if (ct == stdc_type_index || ct == float2_type_index) {
#endif
      return reinterpret_cast<T*>(c_frequency);
#ifdef STORMM_USE_HPC
    }
#else
    }
#endif
    else {
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
      if (ct == stdz_type_index || ct == double2_type_index || ct == cufftz_type_index) {
#  endif
#else
      if (ct == stdz_type_index || ct == double2_type_index) {
#endif
        rtErr("A request for a double-precision complex data pointer was sent to an object "
              "holding single-precision complex data.", "FFTStage", "getFrequencyData");
#ifdef STORMM_USE_HPC
      }
#else
      }
#endif
      else {
        rtErr("A request for an invalid complex data type was received.", "FFTStage",
              "getFrequencyData");
      }
    }
    break;
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> T* FFTStage::getFrequencyData() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  switch (prec) {
  case PrecisionModel::DOUBLE:
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
    if (ct == stdz_type_index || ct == double2_type_index || ct == cufftz_type_index) {
#  endif
#else
    if (ct == stdc_type_index || ct == double2_type_index) {
#endif
      return reinterpret_cast<T*>(z_frequency);
#ifdef STORMM_USE_HPC
    }
#else
    }
#endif
    else {
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
      if (ct == stdc_type_index || ct == float2_type_index || ct == cufftc_type_index) {
#  endif
#else
      if (ct == stdc_type_index || ct == float2_type_index) {
#endif
        rtErr("A request for a single-precision complex data pointer was sent to an object "
              "holding double-precision complex data.", "FFTStage", "getFrequencyData");
#ifdef STORMM_USE_HPC
      }
#else
      }
#endif
      else {
        rtErr("A request for an invalid complex data type was received.", "FFTStage",
              "getFrequencyData");
      }
    }
    break;
  case PrecisionModel::SINGLE:
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
    if (ct == stdc_type_index || ct == float2_type_index || ct == cufftc_type_index) {
#  endif
#else
    if (ct == stdc_type_index || ct == float2_type_index) {
#endif
      return reinterpret_cast<T*>(c_frequency);
#ifdef STORMM_USE_HPC
    }
#else
    }
#endif
    else {
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
      if (ct == stdz_type_index || ct == double2_type_index || ct == cufftz_type_index) {
#  endif
#else
      if (ct == stdz_type_index || ct == double2_type_index) {
#endif
        rtErr("A request for a double-precision complex data pointer was sent to an object "
              "holding single-precision complex data.", "FFTStage", "getFrequencyData");
#ifdef STORMM_USE_HPC
      }
#else
      }
#endif
      else {
        rtErr("A request for an invalid complex data type was received.", "FFTStage",
              "getFrequencyData");
      }
    }
    break;
  }
  __builtin_unreachable();
}

} // namespace stmath
} // namespace stormm
