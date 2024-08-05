// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace card {

//-------------------------------------------------------------------------------------------------
template <typename T>
void deepCopy(Hybrid<T> *destination, const Hybrid<T> &original, const size_t length,
              const size_t destination_offset, const size_t original_offset,
              const GpuDetails &gpu) {
  const size_t actual_length = screenHybridTransfer(destination->size(), original.size(), length,
                                                    destination_offset, original_offset,
                                                    "deepCopy");
  if (actual_length == 0) {
    return;
  }
#ifdef STORMM_USE_HPC
  T* dest_host = destination->data(HybridTargetLevel::HOST);
  T* dest_devc = destination->data(HybridTargetLevel::DEVICE);
  const T* orig_host = original.data(HybridTargetLevel::HOST);
  const T* orig_devc = original.data(HybridTargetLevel::DEVICE);
  bool do_hhc, do_hdc, do_dhc, do_ddc, use_kernel;
  markCopyInstructions(destination->getFormat(), original.getFormat(), &do_hhc, &do_hdc, &do_dhc,
                       &do_ddc, &use_kernel);
  if (use_kernel && false) {
    const void* vorig_host = (do_hdc) ? reinterpret_cast<const void*>(orig_host) : nullptr;
    const void* vorig_devc = (do_dhc || do_ddc) ? reinterpret_cast<const void*>(orig_devc) :
                                                  nullptr;
    void* vdest_host = (do_dhc) ? reinterpret_cast<void*>(dest_host) : nullptr;
    void* vdest_devc = (do_hdc || do_ddc) ? reinterpret_cast<void*>(dest_devc) : nullptr;
    const size_t ct_data = std::type_index(typeid(T)).hash_code();
    launchDeepCopy(vdest_host, vdest_devc, vorig_host, vorig_devc, destination_offset,
                   original_offset, actual_length, ct_data, do_hdc, do_dhc, do_ddc, gpu, 0, 0);
  }
  else {
    if (do_dhc && cudaMemcpy(&dest_host[destination_offset], &orig_devc[original_offset],
                             actual_length * sizeof(T), cudaMemcpyDeviceToHost) != cudaSuccess) {
      rtErr("Error in cudaMemcpy for device-to-host transfer of " +
            std::to_string(actual_length) + " elements of " + std::to_string(sizeof(T)) +
            " bytes each.  Hybrid object " + std::string(destination->getLabel().name) +
            " (format " + getEnumerationName(destination->getFormat()) + ") is receiving from " +
            std::string(original.getLabel().name) + " (format " +
            getEnumerationName(original.getFormat()) + ").", "deepCopy");
    }
    if (do_ddc && cudaMemcpy(&dest_devc[destination_offset], &orig_devc[original_offset],
                             actual_length * sizeof(T), cudaMemcpyDeviceToDevice) != cudaSuccess) {
      rtErr("Error in cudaMemcpy for device-to-device transfer of " +
            std::to_string(actual_length) + " elements of " + std::to_string(sizeof(T)) +
            " bytes each.  Hybrid object " + std::string(destination->getLabel().name) +
            " (format " + getEnumerationName(destination->getFormat()) + ") is receiving from " +
            std::string(original.getLabel().name) + " (format " +
            getEnumerationName(original.getFormat()) + ").", "deepCopy");
    }
    if (do_hdc && cudaMemcpy(&dest_devc[destination_offset], &orig_host[original_offset],
                             actual_length * sizeof(T), cudaMemcpyHostToDevice) != cudaSuccess) {
      rtErr("Error in cudaMemcpy for host-to-device transfer of " +
            std::to_string(actual_length) + " elements of " + std::to_string(sizeof(T)) +
            " bytes each.  Hybrid object " + std::string(destination->getLabel().name) +
            " (format " + getEnumerationName(destination->getFormat()) + ") is receiving from " +
            std::string(original.getLabel().name) + " (format " +
            getEnumerationName(original.getFormat()) + ").", "deepCopy");
    }
  }
  if (do_hhc) {
    memcpy(&dest_host[destination_offset], &orig_host[original_offset], actual_length * sizeof(T));
  }
  if (use_kernel) {
#ifdef STORMM_USE_CUDA
    cudaDeviceSynchronize();
#endif
  }
#else
  T* dest_host = destination->data();
  const T* orig_host = original.data();
  memcpy(&dest_host[destination_offset], &orig_host[original_offset], actual_length * sizeof(T));
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void deepCopy(Hybrid<T> *destination, const Hybrid<T> &original, const GpuDetails &gpu,
              const size_t length, const size_t destination_offset, const size_t original_offset) {
  deepCopy(destination, original, length, destination_offset, original_offset, gpu);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void deepRecast(Hybrid<Tdest> *destination, const Hybrid<Torig> &original, const size_t length,
                const size_t destination_offset, const size_t original_offset, const int dest_bits,
                const int orig_bits, const GpuDetails &gpu) {
  const size_t actual_length = screenHybridTransfer(destination->size(), original.size(), length,
                                                    destination_offset, original_offset,
                                                    "deepRecast");
  if (isHpcVectorType<Tdest>() || isHpcVectorType<Torig>()) {
    rtErr("Conversions involving HPC vector tuples are not allowed (request to recast " +
          std::string(original.getLabel().name) + " to " +
          std::string(destination->getLabel().name) + ").", "deepRecast");
  }
  if (actual_length == 0) {
    return;
  }
#ifdef STORMM_USE_HPC
  Tdest* dest_host = destination->data(HybridTargetLevel::HOST);
  Tdest* dest_devc = destination->data(HybridTargetLevel::DEVICE);
  const Torig* orig_host = original.data(HybridTargetLevel::HOST);
  const Torig* orig_devc = original.data(HybridTargetLevel::DEVICE);
  bool do_hhc, do_hdc, do_dhc, do_ddc, use_kernel;
  markCopyInstructions(destination->getFormat(), original.getFormat(), &do_hhc, &do_hdc, &do_dhc,
                       &do_ddc, &use_kernel);
  if (use_kernel) {
    const void* vorig_host = (do_hdc) ? reinterpret_cast<const void*>(orig_host) : nullptr;
    const void* vorig_devc = (do_dhc || do_ddc) ? reinterpret_cast<const void*>(orig_devc) :
                                                  nullptr;
    void* vdest_host = (do_dhc) ? reinterpret_cast<void*>(dest_host) : nullptr;
    void* vdest_devc = (do_hdc || do_ddc) ? reinterpret_cast<void*>(dest_devc) : nullptr;
    const size_t dest_ct = std::type_index(typeid(Tdest)).hash_code();
    const size_t orig_ct = std::type_index(typeid(Torig)).hash_code();
    launchDeepCopy(vdest_host, vdest_devc, vorig_host, vorig_devc, destination_offset,
                   original_offset, actual_length, dest_ct, orig_ct, do_hdc, do_dhc, do_ddc, gpu,
                   dest_bits, orig_bits);
  }
  else {
  }
  if (do_hhc) {
    Tdest* dest_ptr = destination->data(HybridTargetLevel::HOST);
    const Torig* orig_ptr = original.data(HybridTargetLevel::HOST);
    bool no_conversion = false;
    if (isSignedIntegralScalarType<Torig>()) {
      if (isSignedIntegralScalarType<Tdest>()) {
        if (orig_bits > dest_bits) {
          const llint cfac = pow(2.0, orig_bits - dest_bits);
          for (size_t i = 0; i < length; i++) {
            const llint tmp = orig_ptr[original_offset + i];
            dest_ptr[destination_offset + i] = tmp / cfac;
          }
        }
        else if (orig_bits < dest_bits) {
          const llint cfac = pow(2.0, dest_bits - orig_bits);
          for (size_t i = 0; i < length; i++) {
            const llint tmp = orig_ptr[original_offset + i];
            dest_ptr[destination_offset + i] = tmp * cfac;
          }
        }
        else {
          no_conversion = true;
        }
      }
      else if (isUnsignedIntegralScalarType<Tdest>()) {
        if (orig_bits != dest_bits) {
          rtErr("Conversion between signed and unsigned fixed-precision models of different bit "
                "counts is forbidden (request to recast " + std::string(original.getLabel().name) +
                " to " + std::string(destination->getLabel().name) + " with " +
                std::to_string(orig_bits) + " and " + std::to_string(orig_bits) + " bits after "
                "the point, respectively.", "deepRecast");
        }
        no_conversion = true;
      }
      else if (isFloatingPointScalarType<Tdest>()) {
        const Tdest cfac = pow(2.0, -dest_bits);
        for (size_t i = 0; i < length; i++) {
          dest_ptr[destination_offset + i] = static_cast<Tdest>(orig_ptr[original_offset + i]) *
                                             cfac;
        }
      }
    }
    else if (isUnsignedIntegralScalarType<Torig>()) {
      if (isSignedIntegralScalarType<Tdest>()) {
        if (orig_bits != dest_bits) {
          rtErr("Conversion between signed and unsigned fixed-precision models of different bit "
                "counts is forbidden (request to recast " + std::string(original.getLabel().name) +
                " to " + std::string(destination->getLabel().name) + " with " +
                std::to_string(orig_bits) + " and " + std::to_string(orig_bits) + " bits after "
                "the point, respectively.", "deepRecast");
        }
        no_conversion = true;
      }
      else if (isUnsignedIntegralScalarType<Tdest>()) {
        if (orig_bits > dest_bits) {
          const ullint cfac = pow(2.0, orig_bits - dest_bits);
          for (size_t i = 0; i < length; i++) {
            const ullint tmp = orig_ptr[original_offset + i];
            dest_ptr[destination_offset + i] = tmp / cfac;
          }
        }
        else if (orig_bits < dest_bits) {
          const ullint cfac = pow(2.0, dest_bits - orig_bits);
          for (size_t i = 0; i < length; i++) {
            const ullint tmp = orig_ptr[original_offset + i];
            dest_ptr[destination_offset + i] = tmp * cfac;
          }
        }
        else {
          no_conversion = true;
        }
      }
      else if (isFloatingPointScalarType<Tdest>()) {
        const Tdest cfac = pow(2.0, -orig_bits);
        for (size_t i = 0; i < length; i++) {
          dest_ptr[destination_offset + i] = static_cast<Tdest>(orig_ptr[original_offset + i]) *
                                             cfac;
        }
      }
    }
    else if (isFloatingPointScalarType<Torig>()) {
      if (isSignedIntegralScalarType<Tdest>()) {
        const Torig cfac = pow(2.0, dest_bits);
        for (size_t i = 0; i < length; i++) {
          dest_ptr[destination_offset + i] = llround(orig_ptr[original_offset + i] * cfac); 
        }
      }
      else if (isUnsignedIntegralScalarType<Tdest>()) {
        rtErr("Converting floating point types to unsigned integral fixed-precision models is "
              "forbidden (request to recast " + std::string(original.getLabel().name) +
              " to " + std::string(destination->getLabel().name) + ").", "deepRecast");
      }
      else if (isFloatingPointScalarType<Tdest>()) {
        no_conversion = true;
      }
    }
    if (no_conversion) {
      for (size_t i = 0; i < length; i++) {
        dest_ptr[destination_offset + i] = orig_ptr[original_offset + i];
      }
    }
  }
  if (use_kernel) {
#ifdef STORMM_USE_CUDA
    cudaDeviceSynchronize();
#endif
  }
#else
#endif
}

} // namespace card
} // namespace stormm
