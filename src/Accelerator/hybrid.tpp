// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace card {

//-------------------------------------------------------------------------------------------------
template <typename T> Hybrid<T>::Hybrid(const size_t length_in, const char* tag_in,
                                        const HybridFormat format_in, const HybridKind kind_in) :
    kind{kind_in},
    format{format_in},
    length{length_in},
    element_size{sizeof(T)},
    growth_increment{stmath::getSmallestLot(sizeof(T), hybrid_byte_increment)},
    max_capacity{(length == 0 ? growth_increment : roundUp<size_t>(length, growth_increment))},
    pointer_index{0},
    host_data{nullptr},
#ifdef STORMM_USE_HPC
    devc_data{nullptr},
#endif
    label{assignLabel(tag_in)},
    allocations{0},
    target_serial_number{-1}
{
  // This primary version of the constructor is called by all subsequent versions.  Include
  // checks on the data type and logging for this object in the central ledger here, to have
  // them included in all other Hybrid constructors calling this one in their initializer lists.
  enforceDataTypeLimits();
  gbl_mem_balance_sheet.setEntry(label.serial_number, kind, label, length, element_size,
                                 allocations);
  switch (kind) {
  case HybridKind::ARRAY:
    allocate();
    break;
  case HybridKind::POINTER:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> Hybrid<T>::Hybrid(const std::vector<T> &S_in, const char* tag_in,
                                        const HybridFormat format_in) :
    Hybrid(S_in.size(), tag_in, format_in)
{
#ifdef STORMM_USE_HPC
  switch (format) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    memcpy(host_data, S_in.data(), S_in.size() * sizeof(T));
    upload();
    break;
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::UNIFIED:
    memcpy(host_data, S_in.data(), S_in.size() * sizeof(T));
    break;
  case HybridFormat::DEVICE_ONLY:
    if (cudaMemcpy(devc_data, S_in.data(), S_in.size() * sizeof(T), cudaMemcpyHostToDevice) !=
        cudaSuccess) {
      rtErr("cudaMemcpy failed in std::vector constructor for object " + std::string(label.name) +
            ".", "Hybrid");
    }
    break;
  }
#else
  memcpy(host_data, S_in.data(), S_in.size() * sizeof(T));
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename T> Hybrid<T>::Hybrid(const HybridKind kind_in, const char* tag_in,
                                        const HybridFormat format_in, const size_t length_in) :
    Hybrid(length_in, tag_in, format_in, kind_in)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
Hybrid<T>::~Hybrid() {
  switch (kind) {
  case HybridKind::ARRAY:
    deallocate();
    break;
  case HybridKind::POINTER:
    break;
  }
  gbl_mem_balance_sheet.unsetEntry(label.serial_number);
}

//-------------------------------------------------------------------------------------------------
template <typename T> Hybrid<T>::Hybrid(const Hybrid<T> &original) :
    kind{original.kind},
    format{original.format},
    length{original.length},
    element_size{original.element_size},
    growth_increment{original.growth_increment},
    max_capacity{original.max_capacity},
    pointer_index{original.pointer_index},
    host_data{nullptr},
#ifdef STORMM_USE_HPC
    devc_data{nullptr},
#endif
    label{assignLabel(original.label.name)},
    allocations{0},
    target_serial_number{-1}
{
  // The copy of any Hybrid object (POINTER- or ARRAY-kind) gets its own slot in the global
  // memory Ledger.  From this point forward, it is its own object.
  gbl_mem_balance_sheet.setEntry(label.serial_number, kind, label, length, element_size,
                                 allocations);
  
  // No enforcement of data type limits is needed, as the type of the original has been verified
  switch (kind) {
  case HybridKind::ARRAY:

    // Allocate and copy memory from the original
    allocate();
    switch(format) {
#ifdef STORMM_USE_HPC
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
      memcpy(host_data, original.host_data, length * sizeof(T));
      cudaMemcpy(devc_data, original.devc_data, length * sizeof(T),
                 cudaMemcpyDeviceToDevice);
      break;
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_MOUNTED:
      memcpy(host_data, original.host_data, length * sizeof(T));
      break;
    case HybridFormat::DEVICE_ONLY:
      cudaMemcpy(devc_data, original.devc_data, length * sizeof(T),
                 cudaMemcpyDeviceToDevice);
      break;
#endif
    case HybridFormat::HOST_ONLY:
      memcpy(host_data, original.host_data, length * sizeof(T));
      break;
    }
    break;
  case HybridKind::POINTER:

    // Replicate the remaining effects of setPointer(), where valid to do so.  Unlike the
    // actual setPointer() member function, this is copying the POINTER-kind object's settings.
    // A POINTER-kind object may thus be copied, although still not set to another POINTER-kind
    // object.
    target_serial_number = original.getTargetSerialNumber();
    allocations = original.getAllocations();
    switch (format) {
#ifdef STORMM_USE_HPC
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
      host_data = const_cast<T*>(original.host_data);
      devc_data = const_cast<T*>(original.devc_data);
      break;
    case HybridFormat::HOST_ONLY:
    case HybridFormat::HOST_MOUNTED:
      devc_data = nullptr;
      host_data = const_cast<T*>(original.host_data);
      break;
    case HybridFormat::DEVICE_ONLY:
      host_data = nullptr;
      devc_data = const_cast<T*>(original.devc_data);
      break;
#else
    case HybridFormat::HOST_ONLY:
      host_data = const_cast<T*>(original.host_data);
      break;
#endif
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> Hybrid<T>::Hybrid(Hybrid<T> &&original) :
    kind{original.kind},
    format{original.format},
    length{original.length},
    element_size{original.element_size},
    growth_increment{original.growth_increment},
    max_capacity{original.max_capacity},
    pointer_index{original.pointer_index},
    host_data{original.host_data},
#ifdef STORMM_USE_HPC
    devc_data{original.devc_data},
#endif
    label{assignLabel(original.label.name)},
    allocations{0},
    target_serial_number{-1}
{
  // The copy of any Hybrid object (POINTER- or ARRAY-kind) gets its own slot in the global
  // memory Ledger.  While it would be nice to transfer the index here and now, it is risky
  // to unset the original object's place in the Ledger first.
  gbl_mem_balance_sheet.setEntry(label.serial_number, kind, label, length, element_size,
                                 allocations);

  // No enforcement of data type limits is needed, as the type of the original has been verified.
  // For the move constructor, the imperative is not to copy the original object's data or
  // behavior (that has been accomplished in the initializer list), but rather to complete the
  // hand-off by resetting the original object's information to something that can be readily
  // destroyed without consequence.
  switch (kind) {
  case HybridKind::ARRAY:

    // Prepare the original for clean destruction
    switch(original.format) {
#ifdef STORMM_USE_HPC
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
      original.host_data = nullptr;
      original.devc_data = nullptr;
      break;
    case HybridFormat::HOST_MOUNTED:
      original.host_data = nullptr;
      break;
    case HybridFormat::DEVICE_ONLY:
      original.devc_data = nullptr;
      break;
#endif
    case HybridFormat::HOST_ONLY:
      original.host_data = nullptr;
      break;
    }
    original.length = 0;
    original.allocate();
    break;
  case HybridKind::POINTER:

    // Not knowing whether the object was POINTER- or ARRAY-kind, the target serial number was
    // set to an invalid value of -1.  Given that this is a POINTER-kind array, set the target
    // serial number now.  Executing this logic up in the initializer list would have been dicey.
    target_serial_number = original.target_serial_number;
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> Hybrid<T>& Hybrid<T>::operator=(const Hybrid<T> &other) {

  // Guard against self-assignment
  if (label.serial_number == other.label.serial_number) {
    return *this;
  }

  // Empty the existing Hybrid object
  if (kind == HybridKind::ARRAY && allocations > 0) {
#ifdef STORMM_USE_HPC
    if (host_data != nullptr || devc_data != nullptr) {
      deallocate();
    }
#else
    if (host_data != nullptr) {
      deallocate();
    }
#endif
  }
  gbl_mem_balance_sheet.unsetEntry(label.serial_number);

  // Retrace the inline initialization to rebuild the Hybrid object in the image of the other.
  format = other.format;
  length = other.length;
  element_size = other.element_size;
  growth_increment = other.growth_increment;
  max_capacity = other.max_capacity;
  pointer_index = other.pointer_index;
  label = assignLabel(other.label.name);
  allocations = 0;
  target_serial_number = -1;

  // Replace this object in the global memory ledger
  gbl_mem_balance_sheet.setEntry(label.serial_number, kind, label, length, element_size,
                                 allocations);

  // No enforcement of data type limits is needed, as the type of the other has been verified
  switch (kind) {
  case HybridKind::ARRAY:

    // Allocate and copy memory from the other
    allocate();
    switch(format) {
#ifdef STORMM_USE_HPC
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
      memcpy(host_data, other.host_data, length * sizeof(T));
      cudaMemcpy(devc_data, other.devc_data, length * sizeof(T),
                 cudaMemcpyDeviceToDevice);
      break;
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_MOUNTED:
      memcpy(host_data, other.host_data, length * sizeof(T));
      break;
    case HybridFormat::DEVICE_ONLY:
      cudaMemcpy(devc_data, other.devc_data, length * sizeof(T),
                 cudaMemcpyDeviceToDevice);
      break;
#endif
    case HybridFormat::HOST_ONLY:
      memcpy(host_data, other.host_data, length * sizeof(T));
      break;
    }
    break;
  case HybridKind::POINTER:

    // Replicate the remaining effects of setPointer(), where valid to do so.  Unlike the
    // actual setPointer() member function, this is copying the POINTER-kind object's settings.
    // A POINTER-kind object may thus be copied, although still not set to another POINTER-kind
    // object.
    target_serial_number = other.getTargetSerialNumber();
    allocations = other.getAllocations();
    switch (format) {
#ifdef STORMM_USE_HPC
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
      host_data = const_cast<T*>(other.host_data);
      devc_data = const_cast<T*>(other.devc_data);
      break;
    case HybridFormat::HOST_ONLY:
    case HybridFormat::HOST_MOUNTED:
      host_data = const_cast<T*>(other.host_data);
      devc_data = nullptr;
      break;
    case HybridFormat::DEVICE_ONLY:
      host_data = nullptr;
      devc_data = const_cast<T*>(other.devc_data);
      break;
#else
    case HybridFormat::HOST_ONLY:
      host_data = const_cast<T*>(other.host_data);
      break;
#endif
    }
    break;
  }
  return *this;
}

//-------------------------------------------------------------------------------------------------
template <typename T> Hybrid<T>& Hybrid<T>::operator=(Hybrid<T> &&other) {
  
  // Guard against self-assignment
  if (label.serial_number == other.label.serial_number) {
    return *this;
  }

  // Empty the existing Hybrid object
  if (kind == HybridKind::ARRAY && allocations > 0) {
#ifdef STORMM_USE_HPC
    if (host_data != nullptr || devc_data != nullptr) {
      deallocate();
    }
#else
    if (host_data != nullptr) {
      deallocate();
    }
#endif
  }
  gbl_mem_balance_sheet.unsetEntry(label.serial_number);

  // Retrace the inline initialization to rebuild the Hybrid object in the image of the other.
  format = other.format;
  length = other.length;
  element_size = other.element_size;
  growth_increment = other.growth_increment;
  max_capacity = other.max_capacity;
  pointer_index = other.pointer_index;
  host_data = other.host_data;
#ifdef STORMM_USE_HPC
  devc_data = other.devc_data;
#endif
  label = assignLabel(other.label.name);
  allocations = other.allocations;
  
  // Replace this object in the global memory ledger
  gbl_mem_balance_sheet.setEntry(label.serial_number, kind, label, length, element_size,
                                 allocations);

  switch (kind) {
  case HybridKind::ARRAY:

    // Prepare the original for clean destruction
    switch(other.format) {
#ifdef STORMM_USE_HPC
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
      other.host_data = nullptr;
      other.devc_data = nullptr;
      break;
    case HybridFormat::HOST_MOUNTED:
      other.host_data = nullptr;
      break;
    case HybridFormat::DEVICE_ONLY:
      other.devc_data = nullptr;
      break;
#endif
    case HybridFormat::HOST_ONLY:
      other.host_data = nullptr;
      break;
    }
    other.length = 0;
    other.allocate();
    break;
  case HybridKind::POINTER:

    // Not knowing whether the object was POINTER- or ARRAY-kind, the target serial number was
    // set to an invalid value of -1.  Given that this is a POINTER-kind array, set the target
    // serial number now.  Executing this logic up in the initializer list would have been dicey.
    target_serial_number = other.target_serial_number;
    break;
  }
  return *this;
}

//-------------------------------------------------------------------------------------------------
template <typename T> HybridKind Hybrid<T>::getKind() const {
  return kind;
}

//-------------------------------------------------------------------------------------------------
template <typename T> HybridFormat Hybrid<T>::getFormat() const {
  return format;
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t Hybrid<T>::size() const {
  return length;
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t Hybrid<T>::capacity() const {
  return max_capacity;
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t Hybrid<T>::getPointerIndex() const {
  return pointer_index;
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t Hybrid<T>::getElementSize() const {
  return element_size;
}

//-------------------------------------------------------------------------------------------------
template <typename T> HybridLabel Hybrid<T>::getLabel() const {
  return label;
}

//-------------------------------------------------------------------------------------------------
template <typename T> int Hybrid<T>::getSerialNumber() const {
  return label.serial_number;
}

//-------------------------------------------------------------------------------------------------
template <typename T> int Hybrid<T>::getAllocations() const {
  return allocations;
}

//-------------------------------------------------------------------------------------------------
template <typename T> int Hybrid<T>::getTargetSerialNumber() const {
  return target_serial_number;
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool Hybrid<T>::verifyTarget() const {

  // First test: does the memory ledger contain an active
  // array in the slot that this pointer expects?
  if (gbl_mem_balance_sheet.testActive(target_serial_number) == false) {
    return false;
  }

  // Second test: does the number of allocations recorded for the target match the number
  // first seen by this pointer when it was set?
  return (gbl_mem_balance_sheet.getAllocations(target_serial_number) == allocations);
}

//-------------------------------------------------------------------------------------------------
template <typename T> const T* Hybrid<T>::data(const HybridTargetLevel tier) const {
#ifdef STORMM_USE_HPC
  switch (tier) {
  case HybridTargetLevel::HOST:
    return host_data;
  case HybridTargetLevel::DEVICE:
    return devc_data;
  }
  __builtin_unreachable();
#else
  return host_data;
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename T> T* Hybrid<T>::data(const HybridTargetLevel tier) {
#ifdef STORMM_USE_HPC
  switch (tier) {
  case HybridTargetLevel::HOST:
    return host_data;
  case HybridTargetLevel::DEVICE:
    return devc_data;
  }
  __builtin_unreachable();
#else
  return host_data;
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename T> T Hybrid<T>::readHost(const size_t index) const {

  // Check that there is data on the host
  if (host_data == nullptr) {
    rtErr("No host data exists in Hybrid object " + std::string(label.name) + " (format " +
          getEnumerationName(format) + ", kind " + getEnumerationName(kind) + ").", "Hybrid",
          "readHost");
  }
  if (index >= length) {
    rtErr("Index " + std::to_string(index) + " is unreadable for array of size " +
          std::to_string(length) + " in Hybrid object " + std::string(label.name) + ".", "Hybrid",
          "readHost");
  }
  return host_data[index];
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> Hybrid<T>::readHost(const size_t offset, const size_t count) const {

  // Check that there is data on the host
  if (host_data == nullptr) {
    rtErr("No host data exists in Hybrid object " + std::string(label.name) + " (format " +
          getEnumerationName(format) + ", kind " + getEnumerationName(kind) + ").", "Hybrid",
          "readHost");
  }
  if (offset > length) {
    rtErr("An offset of " + std::to_string(offset) + " is not accessible in an array of " +
          std::to_string(length) + " elements in Hybrid object " + std::string(label.name) + ".",
          "Hybrid", "readHost");
  }
  if (offset + count > length) {
    rtErr("An offset of " + std::to_string(offset) + " cannot provide access to " +
          std::to_string(count) + " elements in Hybrid object " + std::string(label.name) +
          " (maximum capacity " + std::to_string(max_capacity) + ", length " +
          std::to_string(length) + ").", "Hybrid", "readHost");
  }
  std::vector<T> v;
  v.resize(count);
  memcpy(v.data(), &host_data[offset], count * sizeof(T));
  return v;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> Hybrid<T>::readHost() const {
  return readHost(0, length);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Hybrid<T>::readHost(T* v, const size_t offset, const size_t count) const {

  // Check that there is data on the host
  if (host_data == nullptr) {
    rtErr("No host data exists in Hybrid object " + std::string(label.name) + " (format " +
          getEnumerationName(format) + ", kind " + getEnumerationName(kind) + ").", "Hybrid",
          "readHost");
  }
  if (offset > length) {
    rtErr("An offset of " + std::to_string(offset) + " is not accessible in an array of " +
          std::to_string(length) + " elements in Hybrid object " + std::string(label.name) + ".",
          "Hybrid", "readHost");
  }
  if (offset + count > length) {
    rtErr("An offset of " + std::to_string(offset) + " cannot provide access to " +
          std::to_string(count) + " elements in Hybrid object " + std::string(label.name) +
          " (maximum capacity " + std::to_string(max_capacity) + ", length " +
          std::to_string(length) + ").", "Hybrid", "readHost");
  }
  memcpy(v, &host_data[offset], count * sizeof(T));
}

//-------------------------------------------------------------------------------------------------
template <typename T> void Hybrid<T>::putHost(const T value, const size_t index) {

  // Check that there is data on the host
  if (host_data == nullptr) {
    rtErr("No host data exists in Hybrid object " + std::string(label.name) + " (format " +
          getEnumerationName(format) + ", kind " + getEnumerationName(kind) + ").", "Hybrid",
          "putHost");
  }
  if (index >= length) {
    rtErr("Index " + std::to_string(index) + " is unwritable for array of size " +
          std::to_string(length) + " in Hybrid object " + std::string(label.name) + ".", "Hybrid",
          "putHost");
  }
  memcpy(&host_data[index], &value, sizeof(T));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Hybrid<T>::putHost(const std::vector<T> &values, const size_t offset, const size_t count) {

  // Check that there is data on the host
  if (host_data == nullptr) {
    rtErr("No host data exists in Hybrid object " + std::string(label.name) + " (format " +
          getEnumerationName(format) + ", kind " + getEnumerationName(kind) + ").", "Hybrid",
          "putHost");
  }
  if (count > values.size()) {
    rtErr("There is only enough data provided to import " + std::to_string(values.size()) +
          " data elements out of " + std::to_string(count) + " requested into Hybrid object " +
          std::string(label.name), "Hybrid", "putHost");
  }
  if (offset + count > length) {
    rtErr("There is only enough space available to import " + std::to_string(length - offset) +
          " data elements out of " + std::to_string(count) + " requested into Hybrid object " +
          std::string(label.name), "Hybrid", "putHost");
  }
  memcpy(&host_data[offset], values.data(), count * sizeof(T));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Hybrid<T>::putHost(const std::vector<T> &values) {
  putHost(values, 0, values.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T>
size_t Hybrid<T>::putHost(Hybrid<T> *target, const std::vector<T> &values, const size_t offset,
                          const size_t padding, const size_t count) {
  const size_t actual_count = (count == 0) ? values.size() : count;
  const size_t padded_count = (padding > 0) ? roundUp<size_t>(actual_count, padding) :
                                              actual_count;
  if (offset > target->length || offset + padded_count > target->length) {
    rtErr("The Hybrid object " + std::string(target->label.name) + " (targeted by " +
          std::string(label.name) + " at index " + std::to_string(offset) + ") holds only " +
          std::to_string(target->length) + " elements, not enough to put " +
          std::to_string(actual_count) + " with " + std::to_string(padded_count - actual_count) +
          " zero padding.", "Hybrid", "putHost");
  }
  setPointer(target, offset, actual_count);

  // Check that there is data on the host
  if (host_data == nullptr) {
    rtErr("No host data exists in Hybrid object " + std::string(label.name), "Hybrid",
          "putHost");
  }
  memcpy(host_data, values.data(), actual_count * sizeof(T));
  memset(&host_data[actual_count], 0, (padded_count - actual_count) * sizeof(T));
  return offset + padded_count;
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
template <typename T> void Hybrid<T>::upload(const size_t start_pos, const size_t length_copy) {
  const size_t n_copy = (length_copy == 0) ? length - start_pos : length_copy;
  if (start_pos > length) {
    rtErr("Cannot upload elements of " + std::string(label.name) + " starting at index " +
          std::to_string(start_pos) + " (the array has " + std::to_string(length) +
          " elements in all).", "Hybrid", "upload");
  }
  if (start_pos + n_copy > length) {
    rtErr("Cannot upload " + std::to_string(n_copy) + " elements of " +
          std::string(label.name) + " starting at index " + std::to_string(start_pos) +
          " (the array had " + std::to_string(length) + " elements in all).", "Hybrid",
          "upload");
  }
  switch (format) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    if (cudaMemcpy(devc_data + start_pos, host_data + start_pos, n_copy * sizeof(T),
                   cudaMemcpyHostToDevice) != cudaSuccess) {
      rtErr("Failure in cudaMemcpy for " + std::string(label.name), "Hybrid", "upload");
    }
    break;
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::DEVICE_ONLY:
  case HybridFormat::HOST_MOUNTED:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void Hybrid<T>::download(const size_t start_pos, const size_t length_copy) {
  const size_t n_copy = (length_copy == 0) ? length - start_pos : length_copy;
  if (start_pos > length) {
    std::string h_name(label.name);
    rtErr("Cannot download elements of " + std::string(label.name) + " starting at index " +
          std::to_string(start_pos) + " (the array has " + std::to_string(length) +
          " elements).", "Hybrid", "download");
  }
  if (start_pos + n_copy > length) {
    rtErr("Cannot download " + std::to_string(n_copy) + " elements of " +
          std::string(label.name) + " starting at index " + std::to_string(start_pos) + " (the "
          "array had " + std::to_string(length) + " elements in all).", "Hybrid", "download");
  }
  switch (format) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    if (cudaMemcpy(host_data + start_pos, devc_data + start_pos, n_copy * sizeof(T),
                   cudaMemcpyDeviceToHost) != cudaSuccess) {
      rtErr("Failure in cudaMemcpy for " + std::string(label.name), "Hybrid", "download");
    }
    break;
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::DEVICE_ONLY:
  case HybridFormat::HOST_MOUNTED:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> T Hybrid<T>::readDevice(const size_t index) const {

  // Check that there is data on the device
  if (devc_data == nullptr) {
    rtErr("No device data exists in Hybrid object " + std::string(label.name) + ".", "Hybrid",
          "readDevice");
  }
  if (index > length) {
    rtErr("Hybrid object " + std::string(label.name) + " does not have " + std::to_string(index) +
          " elements.  Its length is " + std::to_string(length) + ".", "Hybrid", "readDevice");
  }
  T rval;
  if (cudaMemcpy(&rval, &devc_data[index], sizeof(T), cudaMemcpyDeviceToHost) != cudaSuccess) {
    rtErr("Error in cudaMemcpy (downloading device data in " + std::string(label.name) + ").",
          "Hybrid", "readDevice");
  }
  return rval;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> Hybrid<T>::readDevice(const size_t offset, const size_t count) const {

  // Check that there is data on the device
  if (devc_data == nullptr) {
    rtErr("No device data exists in Hybrid object " + std::string(label.name) + ".", "Hybrid",
          "readDevice");
  }
  if (offset + count > length) {
    rtErr("Hybrid object " + std::string(label.name) + " does not have " + std::to_string(count) +
          " elements after index " + std::to_string(offset) + ".  Its length is " +
          std::to_string(length) + ".", "Hybrid", "readDevice");
  }
  std::vector<T> v;
  v.resize(count);
  if (cudaMemcpy(v.data(), &devc_data[offset], count * sizeof(T), cudaMemcpyDeviceToHost) !=
      cudaSuccess) {
    rtErr("Error in cudaMemcpy (downloading device data in " + std::string(label.name) + ").",
          "Hybrid", "readDevice");
  }
  return v;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> Hybrid<T>::readDevice() const {
  return readDevice(0, length);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Hybrid<T>::readDevice(T* v, const size_t offset, const size_t count) const {

  // Check that there is data on the device
  if (devc_data == nullptr) {
    rtErr("No device data exists in Hybrid object " + std::string(label.name) + ".", "Hybrid",
          "readDevice");
  }
  if (offset + count > length) {
    rtErr("Hybrid object " + std::string(label.name) + " does not have " + std::to_string(count) +
          " elements after index " + std::to_string(offset) + ".  Its length is " +
          std::to_string(length) + ".", "Hybrid", "readDevice");
  }
  if (cudaMemcpy(v, &devc_data[offset], count * sizeof(T), cudaMemcpyDeviceToHost) !=
      cudaSuccess) {
    rtErr("Error in cudaMemcpy (downloading device data in " + std::string(label.name) + ").",
          "Hybrid", "readDevice");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void Hybrid<T>::putDevice(const T value, const size_t index) {

  // Check that there is data on the host
  if (devc_data == nullptr) {
    rtErr("No device data exists in Hybrid object " + std::string(label.name), "Hybrid",
          "putDevice");
  }
  if (index >= length) {
    rtErr("Index " + std::to_string(index) + " is unwritable for array of size " +
          std::to_string(length) + ".", "Hybrid", "putDevice");
  }
  cudaMemcpy(&devc_data[index], &value, sizeof(T), cudaMemcpyHostToDevice);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Hybrid<T>::putDevice(const std::vector<T> &values, const size_t offset, const size_t count) {

  // Check that there is data on the device
  if (devc_data == nullptr) {
    rtErr("No device data exists in Hybrid object " + std::string(label.name), "Hybrid",
          "putDevice");
  }
  if (count > values.size()) {
    rtErr("There is only enough data provided to import " + std::to_string(values.size()) +
          " data elements out of " + std::to_string(count) + " requested into Hybrid object " +
          std::string(label.name), "Hybrid", "putDevice");
  }
  if (offset + count > length) {
    rtErr("There is only enough space available to import " + std::to_string(length - offset) +
          " data elements out of " + std::to_string(count) + " requested into Hybrid object " +
          std::string(label.name), "Hybrid", "putDevice");
  }
  cudaMemcpy(&devc_data[offset], values.data(), count * sizeof(T), cudaMemcpyHostToDevice);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Hybrid<T>::putDevice(const std::vector<T> &values) {
  putDevice(values, 0, values.size());
}
#endif

//-------------------------------------------------------------------------------------------------
template <typename T> void Hybrid<T>::shrinkToFit() {

  // Check that this is an actual array
  switch (kind) {
  case HybridKind::POINTER:
    max_capacity = length;
    break;
  case HybridKind::ARRAY:
    reallocate(length, length);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void Hybrid<T>::pushBack(const T element) {
  if (length < max_capacity) {
    host_data[length] = element;
    length += 1;
  }
  else {
    switch (kind) {
    case HybridKind::POINTER:
      rtErr("pointer object cannot assign element " + std::to_string(length) +
            " beyond its capacity.", "Hybrid", "pushBack");
    case HybridKind::ARRAY:
      reallocate(length + 1, 0);
      break;
    }
    host_data[length - 1] = element;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void Hybrid<T>::pushBack(const T* elements, const size_t element_count) {
  if (length + element_count < max_capacity) {
    for (size_t i = 0; i < element_count; i++) {
      host_data[length + i] = elements[i];
    }
    length += element_count;
  }
  else {
    const size_t old_length = length;
    switch (kind) {
    case HybridKind::POINTER:
      rtErr("A " + getEnumerationName(HybridKind::POINTER) + "-kind object cannot assign "
            "element " + std::to_string(length) + " beyond its capacity.", "Hybrid", "pushBack");
    case HybridKind::ARRAY:
      reallocate(length + element_count, 0);
      break;
    }
    for (size_t i = 0; i < element_count; i++) {
      host_data[old_length + i] = elements[i];
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void Hybrid<T>::pushBack(const std::vector<T> &elements) {
  pushBack(elements.data(), elements.size());
}

//-------------------------------------------------------------------------------------------------
template <typename T> void Hybrid<T>::pushBack(const Hybrid<T> &elements) {
  if (label.serial_number == elements.label.serial_number) {

    // Take care when pushing a Hybrid object to the back of itself.  This requires a special
    // procedure but a full out-of-place copy is not necessary.
    if (length + length < max_capacity) {
      pushBack(host_data, length);
    }
    else {
      const size_t old_length = length;
      switch (kind) {
      case HybridKind::POINTER:
        rtErr("A " + getEnumerationName(HybridKind::POINTER) + "-kind object cannot assign "
              "element " + std::to_string(length) + " beyond its capacity.", "Hybrid", "pushBack");
      case HybridKind::ARRAY:
        reallocate(length + length, 0);
        break;
      }
      for (size_t i = 0; i < old_length; i++) {
        host_data[old_length + i] = host_data[i];
      }
    }
  }
  else {
    pushBack(elements.data(), elements.size());
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void Hybrid<T>::resize(const size_t new_length) {
  if (new_length <= max_capacity) {
    if (new_length > length) {
      const std::string errmsg = std::to_string(length) + " to " + std::to_string(new_length) +
                                 " x " + std::to_string(sizeof(T)) + " bytes in Hybrid object " +
                                 std::string(label.name) + ".";
      const size_t delta_len = new_length - length;
      switch (format) {
#ifdef STORMM_USE_HPC
      case HybridFormat::EXPEDITED:
        if (cudaMemset((void *)(&devc_data[length]), 0, delta_len * sizeof(T)) != cudaSuccess) {
          rtErr("cudaMemset unsuccessful for " + errmsg, "Hybrid", "allocate");
        }
        memset(&host_data[length], 0, delta_len * sizeof(T));
        break;
      case HybridFormat::DECOUPLED:
        if (cudaMemset((void *)(&devc_data[length]), 0, delta_len * sizeof(T)) != cudaSuccess) {
          rtErr("cudaMemset unsuccessful for " + errmsg, "Hybrid", "allocate");
        }
        memset(&host_data[length], 0, delta_len * sizeof(T));
        break;
      case HybridFormat::UNIFIED:
        memset(&host_data[length], 0, delta_len * sizeof(T));
        break;
      case HybridFormat::HOST_ONLY:
        memset(&host_data[length], 0, delta_len * sizeof(T));
        break;
      case HybridFormat::DEVICE_ONLY:
        if (cudaMemset((void *)(&devc_data[length]), 0, delta_len * sizeof(T)) != cudaSuccess) {
          rtErr("cudaMemset unsuccessful for " + errmsg, "Hybrid", "allocate");
        }
        break;
      case HybridFormat::HOST_MOUNTED:
        memset(&host_data[length], 0, delta_len * sizeof(T));
        break;
#else
      case HybridFormat::HOST_ONLY:
        memset(&host_data[length], 0, delta_len * sizeof(T));
        break;
#endif
          
      }
    }
    length = new_length;
  }
  else {
    switch (kind) {
    case HybridKind::POINTER:
      rtErr("POINTER-kind Hybrid object cannot resize to " + std::to_string(new_length) +
            " elements beyond its capacity of " + std::to_string(max_capacity) + ".  Reallocate "
            "the target array and refresh the pointer.", "Hybrid", "resize");
    case HybridKind::ARRAY:
      reallocate(new_length, new_length);
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void Hybrid<T>::resize(const size_t new_length, const T value) {
  const size_t original_length = length;
  resize(new_length);
  for (size_t i = original_length; i < new_length; i++) {
    host_data[i] = value;
  }
#ifdef STORMM_USE_HPC
  if (new_length > original_length) {
    upload(original_length, new_length - original_length);
  }
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Hybrid<T>::setPointer(Hybrid<T> *target, const size_t position, const llint new_length) {

  // Check that this Hybrid object is a pointer, and that the target is an array
  const HybridLabel tlbl = target->label;
  if (kind != HybridKind::POINTER) {
    rtErr("Hybrid object " + std::string(label.name) + " is not a pointer (targeted by " +
          std::string(label.name) + ").", "Hybrid", "setPointer");
  }
  if (target->kind != HybridKind::ARRAY) {
    rtErr("Hybrid object " + std::string(tlbl.name) + " is not an array (targeted by " +
          std::string(label.name) + ").", "Hybrid", "setPointer");
  }

  // Pointers can only be set to a Hybrid of similar memory format.
  if (format != target->format) {
    rtErr("Hybrid object " + std::string(tlbl.name) + " has a different format (" +
          getEnumerationName(target->format) + ") than " + std::string(label.name) + " (" +
          getEnumerationName(format) + ").", "Hybrid", "setPointer");
  }

  // Confirm that the location is within the bounds of the target
  if (position > target->max_capacity || (position == target->max_capacity && new_length > 0)) {
    rtErr("object " + std::string(label.name) + " cannot point to element " +
          std::to_string(position) + " of target " + std::string(tlbl.name) + " (length " +
          std::to_string(target->capacity()) + ").", "Hybrid", "setPointer");
  }

  // Set the target, length, and capacity of the pointer
  pointer_index = position;
  if (new_length >= 0) {
    if (new_length > target->max_capacity - position) {
      rtErr("object " + std::string(tlbl.name) + " of length " +
            std::to_string(target->max_capacity) + " cannot accept a pointer from " +
            std::string(label.name) + " at position " + std::to_string(position) +
            " of length " + std::to_string(new_length), "Hybrid", "setPointer");
    }
    length = new_length;
    max_capacity = new_length;
  }
  else {
    length = target->length - position;
    max_capacity = target->max_capacity - position;
  }

  // Set the format of the pointer to that of the target
  format = target->format;

  // Record the target's serial number, and take the target's number of allocations as this
  // object's number of allocations.  If the validity of this pointer needs to be checked in
  // the future, a match between the allocations of this object and its target is a strong
  // indicator that the target data of these pointers is current.
  target_serial_number = target->label.serial_number;
  allocations = target->allocations;

  // Set the host and device data pointers as appropriate
  switch (format) {
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
    host_data = &target->data(HybridTargetLevel::HOST)[position];
    devc_data = &target->data(HybridTargetLevel::DEVICE)[position];
    break;
  case HybridFormat::HOST_ONLY:
    host_data = &target->data(HybridTargetLevel::HOST)[position];
    devc_data = nullptr;
    break;
  case HybridFormat::DEVICE_ONLY:
    host_data = nullptr;
    devc_data = &target->data(HybridTargetLevel::DEVICE)[position];
    break;
#else
  case HybridFormat::HOST_ONLY:
    host_data = &target->data(HybridTargetLevel::HOST)[position];
    break;
#endif
  }

  // Record changes in the ledger
  gbl_mem_balance_sheet.setEntry(label.serial_number, kind, label, length, element_size,
                                 allocations);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Hybrid<T>::swapTarget(Hybrid<T> *new_target, const ExceptionResponse policy) {
  if (kind != HybridKind::POINTER) {
    rtErr("The non-POINTER kind Hybrid object " + std::string(label.name) + " cannot be set to "
          "target any other Hybrid object.", "Hybrid", "swapTarget");
  }
  if (new_target->kind != HybridKind::ARRAY) {
    rtErr("The Hybrid object " + std::string(new_target->label.name) + " (targeted by " +
          std::string(label.name) + ") must be ARRAY-kind.", "Hybrid", "swapTarget");
  }
  if (new_target->format != format) {
    rtErr("Hybrid object " + std::string(new_target->label.name) + " has a different format (" +
          getEnumerationName(new_target->format) + ") than " + std::string(label.name) + " (" +
          getEnumerationName(format) + ").", "Hybrid", "swapTarget");
  }
  if (target_serial_number < 0) {

    // Is there is no current target, then no swapping can occur.  This may or may not lead to
    // errors later in the program.  A range of responses is offered, with a SILENT policy of
    // returning immediately being the default.
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Hybrid object " + std::string(label.name) + " does not currently point to an array.",
            "Hybrid", "swapTarget");
    case ExceptionResponse::WARN:
      rtWarn("Hybrid object " + std::string(label.name) + " does not currently point to an "
             "array.  No target swap will occur.", "Hybrid", "swapTarget");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return;
  }
  setPointer(new_target, pointer_index, length);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const Hybrid<T> Hybrid<T>::getPointer(const size_t position, const size_t new_length,
                                      const char* tag_in) const {
  if (kind != HybridKind::ARRAY) {
    rtErr("Error creating a POINTER-kind Hybrid object to " + std::string(label.name) + ".  Only "
          "an ARRAY-kind Hybrid object may produce a pointer to some point in its data.", "Hybrid",
          "getPointer");
  }

  // Check that the target position is within bounds
  if (position >= max_capacity) {
    rtErr("Position " + std::to_string(position) + " is invalid for Hybrid object " +
          std::string(label.name) + " (capacity " + std::to_string(max_capacity) + ").", "Hybrid",
          "getPointer");
  }

  // Check the length
  size_t actual_length;
  if (new_length == 0) {
    actual_length = (length >= position) ? max_capacity - position : length - position;
  }
  else {
    if (position + length > max_capacity) {
      rtErr("A POINTER-king Hybrid of length " + std::to_string(new_length) + " cannot be created "
            "from position " + std::to_string(position) + " of Hybrid object " +
            std::string(label.name) + " (maximum capacity " + std::to_string(max_capacity) + ").",
            "Hybrid", "getPointer");
    }
    actual_length = length;
  }
  
  // Create and return the pointer object.  This must mimic the effects of setPointer, but with
  // temporary const-casting (protected within the context of this function) to set the array
  // pointers.  The const casting will be safe because of the return type's guaranteed const-ness.
  Hybrid<T> result(HybridKind::POINTER, tag_in, format);

  // Piece together the new Hybrid array, as a POINTER-kind object into this object
  if (new_length > 0) {
    result.max_capacity = new_length;
  }
  else {
    result.max_capacity = max_capacity - position;
  }
  result.pointer_index = position;
  result.allocations = allocations;
  result.target_serial_number = label.serial_number;

  // Set the host and device data pointers as appropriate
  switch (format) {
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
    host_data = const_cast<T*>(&host_data[position]);
    devc_data = const_cast<T*>(&devc_data[position]);
    break;
  case HybridFormat::HOST_ONLY:
    host_data = const_cast<T*>(&host_data[position]);
    devc_data = nullptr;
    break;
  case HybridFormat::DEVICE_ONLY:
    host_data = nullptr;
    devc_data = const_cast<T*>(&devc_data[position]);
    break;
#else
  case HybridFormat::HOST_ONLY:
    host_data = const_cast<T*>(&host_data[position]);
    break;
#endif
  }

  // Record changes in the ledger
  gbl_mem_balance_sheet.setEntry(result.serial_number, result.kind, result.label, result.length,
                                 result.element_size, result.allocations);
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
Hybrid<T> Hybrid<T>::getPointer(const size_t position, const size_t new_length,
                                const char* tag_in) {
  Hybrid<T> result(HybridKind::POINTER, tag_in, format);
  result.setPointer(this, position, new_length);
  return result;
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
template <typename T> const T* Hybrid<T>::getDeviceValidHostPointer() const {
  T* proto;
  switch (format) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::UNIFIED:
#  ifdef STORMM_USE_CUDA
    if (cudaHostGetDevicePointer((void **)&proto, (void *)host_data, 0) != cudaSuccess) {
      rtErr("Error in retrieving device-accessible pointer to host data in Hybrid object " +
            std::string(label.name) + " of type " +
            std::string(std::type_index(typeid(T)).name()) + ".", "Hybrid",
            "getDeviceValidHostPointer");
    }
#  endif
    break;
  case HybridFormat::DECOUPLED:
    rtErr("In order to retrieve a pointer to host memory valid on the GPU device, the host "
          "memory must be allocated in page-locked format with cudaHostAlloc().  " +
          getEnumerationName(format) + " is not compatible.", "Hybrid",
          "getDeviceValidHostPointer");
  case HybridFormat::DEVICE_ONLY:
    rtErr("No host memory is available in " + getEnumerationName(format) + " allocation format.",
          "Hybrid", "getDeviceValidHostPointer");
  }
  return proto;
}

//-------------------------------------------------------------------------------------------------
template <typename T> T* Hybrid<T>::getDeviceValidHostPointer() {
  T* proto;
  switch (format) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::UNIFIED:
#  ifdef STORMM_USE_CUDA
    if (cudaHostGetDevicePointer((void **)&proto, (void *)host_data, 0) != cudaSuccess) {
      rtErr("Error in retrieving device-accessible pointer to host data in Hybrid object " +
            std::string(label.name) + " of type " +
            std::string(std::type_index(typeid(T)).name()) + ".", "Hybrid",
            "getDeviceValidHostPointer");
    }
#  endif
    break;
  case HybridFormat::DECOUPLED:
    rtErr("In order to retrieve a pointer to host memory valid on the GPU device, the host "
          "memory must be allocated in page-locked format with cudaHostAlloc().  " +
          getEnumerationName(format) + " is not compatible.", "Hybrid",
          "getDeviceValidHostPointer");
  case HybridFormat::HOST_ONLY:
    rtErr("Pageable host memory in a " + getEnumerationName(format) + " allocation is "
          "inaccessible to the GPU.", "Hybrid", "getDeviceValidHostPointer");
  case HybridFormat::DEVICE_ONLY:
    rtErr("No host memory is available in " + getEnumerationName(format) + " allocation format.",
          "Hybrid", "getDeviceValidHostPointer");
  }
  return proto;
}
#endif

//-------------------------------------------------------------------------------------------------
template <typename T> void Hybrid<T>::enforceDataTypeLimits() {

  using data_types::isScalarType;
  using data_types::isHpcVectorType;

  if (isScalarType<T>() || isHpcVectorType<T>()) {
    return;
  }
  else {
    rtErr("data type " + std::string(std::type_index(typeid(T)).name()) + " is not a valid "
          "Hybrid data type.  Restrict allocations of this template object to POD types, i.e. "
          "double, and GPU-purposed vector types, i.e. int4.", "Hybrid",
          "enforceDataTypeLimits");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> void Hybrid<T>::allocate() {

  // Make sure that this is not a pointer
  assert (kind == HybridKind::ARRAY);

  // No memory should currently be allocated
  assert (host_data == nullptr);
#ifdef STORMM_USE_HPC
  assert (devc_data == nullptr);
#endif
  // Allocate memory
  std::string errmsg = std::to_string(length) + " x " + std::to_string(sizeof(T)) +
                       " bytes in Hybrid object " + std::string(label.name) + ".";
  switch (format) {
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
    if (cudaMalloc((void **)&devc_data, max_capacity * sizeof(T)) != cudaSuccess) {
      rtErr("cudaMalloc unsuccessful for " + errmsg, "Hybrid", "allocate");
    }
    if (cudaMemset((void *)devc_data, 0, length * sizeof(T)) != cudaSuccess) {
      rtErr("cudaMemset unsuccessful for " + errmsg, "Hybrid", "allocate");
    }
    if (cudaHostAlloc((void **)&host_data, max_capacity * sizeof(T), cudaHostAllocMapped) !=
        cudaSuccess) {
      rtErr("cudaHostAlloc unsuccessful for " + errmsg, "Hybrid", "allocate");
    }
    memset(host_data, 0, length * sizeof(T));
    break;
  case HybridFormat::DECOUPLED:
    if (cudaMalloc((void **)&devc_data, max_capacity * sizeof(T)) != cudaSuccess) {
      rtErr("cudaMalloc unsuccessful for " + errmsg, "Hybrid", "allocate");
    }
    if (cudaMemset((void *)devc_data, 0, length * sizeof(T)) != cudaSuccess) {
      rtErr("cudaMemset unsuccessful for " + errmsg, "Hybrid", "allocate");
    }
    host_data = new T[max_capacity];
    memset(host_data, 0, length * sizeof(T));
    break;
  case HybridFormat::UNIFIED:
    if (cudaMallocManaged(&host_data, max_capacity * sizeof(T)) != cudaSuccess) {
      rtErr("cudaMallocManaged unsuccessful for " + errmsg, "Hybrid", "allocate");
    }
    memset(host_data, 0, length * sizeof(T));
    devc_data = host_data;
    break;
  case HybridFormat::HOST_ONLY:
    host_data = new T[max_capacity];
    memset(host_data, 0, length * sizeof(T));
    break;
  case HybridFormat::DEVICE_ONLY:
    if (cudaMalloc((void **)&devc_data, max_capacity * sizeof(T)) != cudaSuccess) {
      rtErr("cudaMalloc unsuccessful for " + errmsg, "Hybrid", "allocate");
    }
    if (cudaMemset((void *)devc_data, 0, length * sizeof(T)) != cudaSuccess) {
      rtErr("cudaMemset unsuccessful for " + errmsg, "Hybrid", "allocate");
    }
    break;
  case HybridFormat::HOST_MOUNTED:
    if (cudaHostAlloc((void **)&host_data, max_capacity * sizeof(T), cudaHostAllocMapped) !=
        cudaSuccess) {
      rtErr("cudaHostAlloc unsuccessful for " + errmsg, "Hybrid", "allocate");
    }
    memset(host_data, 0, length * sizeof(T));
    break;
#else
  case HybridFormat::HOST_ONLY:
    host_data = new T[max_capacity];
    memset(host_data, 0, length * sizeof(T));
    break;
#endif
  }
  
  // Increment this object's allocations and record the result in the ledger
  allocations += 1;
  gbl_mem_balance_sheet.logMemory(max_capacity, element_size, format, label, allocating);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void Hybrid<T>::deallocate() {

  // Make sure that this is not a pointer
  assert (kind == HybridKind::ARRAY);

  // Memory should currently be allocated
  switch (format) {
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
    assert (devc_data != nullptr);
    assert (host_data != nullptr);
    break;
  case HybridFormat::HOST_ONLY:
    assert (host_data != nullptr);
    break;
  case HybridFormat::DEVICE_ONLY:
    assert (devc_data != nullptr);
    break;
  case HybridFormat::HOST_MOUNTED:
    assert (host_data != nullptr);
    break;
#else
  case HybridFormat::HOST_ONLY:
    assert (host_data != nullptr);
    break;
#endif
  }

  // If no memory is currently allocated, do nothing
  switch (format) {
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
    if (cudaFreeHost(host_data) != cudaSuccess) {
      rtErr("cudaFreeHost failed in object " + std::string(label.name) + ".", "Hybrid",
            "deallocate");
    }
    if (cudaFree(devc_data) != cudaSuccess) {
      rtErr("cudaFree failed in object " + std::string(label.name) + ".", "Hybrid",
            "deallocate");
    }
    break;
  case HybridFormat::DECOUPLED:
    delete[] host_data;
    if (cudaFree(devc_data) != cudaSuccess) {
      rtErr("cudaFree failed in object " + std::string(label.name) + ".", "Hybrid",
            "deallocate");
    }
    break;
  case HybridFormat::UNIFIED:
    if (cudaFree(host_data) != cudaSuccess) {
      rtErr("cudaFree failed in object " + std::string(label.name) + ".", "Hybrid",
            "deallocate");
    }
    devc_data = nullptr;
    break;
  case HybridFormat::HOST_ONLY:
    delete[] host_data;
    break;
  case HybridFormat::DEVICE_ONLY:
    if (cudaFree(devc_data) != cudaSuccess) {
      rtErr("cudaFree failed in object " + std::string(label.name) + ".", "Hybrid",
            "deallocate");
    }
    break;
  case HybridFormat::HOST_MOUNTED:
    if (cudaFreeHost(host_data) != cudaSuccess) {
      rtErr("cudaFreeHost failed in object " + std::string(label.name) + ".", "Hybrid",
            "deallocate");
    }
    break;
#else
  case HybridFormat::HOST_ONLY:
    delete[] host_data;
    break;
#endif
  }

  // Set data pointers to NULL
  host_data = nullptr;
#ifdef STORMM_USE_HPC
  devc_data = nullptr;
#endif
  // Record the result in the ledger, then zero the capacity and current size.
  gbl_mem_balance_sheet.logMemory(max_capacity, element_size, format, label, deallocating);
  length = 0;
  max_capacity = 0;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Hybrid<T>::reallocate(const size_t new_length, const size_t new_capacity) {

  // Make sure that this is not a pointer
  assert (kind == HybridKind::ARRAY);

  // Determine the overlap of old and new arrays
  const size_t overlap = std::min(length, new_length);

  // Allocate buffers for existing memory according to the appropriate format
  T* host_data_buffer;
#ifdef STORMM_USE_HPC
  T* devc_data_buffer;
  switch (format) {
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    if (overlap != 0) {
      host_data_buffer = new T[overlap];
      memcpy(host_data_buffer, host_data, overlap * sizeof(T));
    }
    break;
  case HybridFormat::DEVICE_ONLY:
    break;
  }
  switch (format) {
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
    break;
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::DEVICE_ONLY:
    if (overlap != 0) {
      if (cudaMalloc((void **)&devc_data_buffer, overlap * sizeof(T)) != cudaSuccess) {
        rtErr("cudaMalloc failed to allocate a device-side buffer array of size " +
              std::to_string(overlap) + " for transferring overlapping data in object " +
              std::string(label.name) + ".", "Hybrid", "reallocate");
      }
      if (cudaMemcpy(devc_data_buffer, devc_data, overlap * sizeof(T), cudaMemcpyDeviceToDevice) !=
          cudaSuccess) {
        rtErr("cudaMemcpy failed to transfer overlap for arrays of " + std::to_string(length) +
              " and " + std::to_string(new_length) + " with elements of " +
              std::to_string(sizeof(T)) + " bytes in object " + std::string(label.name) + ".",
              "Hybrid", "reallocate");
      }
    }
    break;
  }
  const size_t old_length = length;
#else
  if (overlap != 0) {
    host_data_buffer = new T[overlap];
    memcpy(host_data_buffer, host_data, overlap * sizeof(T));
  }
#endif
  // Free existing memory
  deallocate();

  // Determine the new maximum array length
  size_t local_cap;
  if (new_capacity == 0) {
    local_cap = getPaddedMemorySize(new_length, growth_increment, sizeof(T));
  }
  else {
    local_cap = (new_capacity < new_length) ? roundUp<size_t>(new_length, growth_increment) :
                                              roundUp<size_t>(new_capacity, growth_increment);
  }

  // Allocate for new array lengths
  length = new_length;
  max_capacity = local_cap;
  allocate();

  // Transfer buffered data as appropriate, then free the buffered data
#ifdef STORMM_USE_HPC
  switch (format) {
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY: 
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    if (overlap != 0) {
      memcpy(host_data, host_data_buffer, overlap * sizeof(T));
      delete[] host_data_buffer;
    }
    break;
  case HybridFormat::DEVICE_ONLY:
    break;
  }
  switch (format) {
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY: 
  case HybridFormat::HOST_MOUNTED:
    break;  
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::DEVICE_ONLY:
    if (overlap != 0) {
      if (cudaMemcpy(devc_data, devc_data_buffer, overlap * sizeof(T), cudaMemcpyDeviceToDevice) !=
        cudaSuccess) {
        rtErr("cudaMemcpy failed to merge arrays of " + std::to_string(old_length) + " and " +
              std::to_string(length) + " with elements of " + std::to_string(sizeof(T)) +
              " bytes in object " + std::string(label.name) + ".", "Hybrid", "reallocate");
      }
      if (cudaFree(devc_data_buffer) != cudaSuccess) {
        rtErr("cudaFree of overlap transfer buffer failed in object " + std::string(label.name) +
              ".", "Hybrid", "reallocate");
      }
    }
    break;
  }
#else
  if (overlap != 0) {
    memcpy(host_data, host_data_buffer, overlap * sizeof(T));
    delete[] host_data_buffer;
  }
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename T> HybridLabel Hybrid<T>::assignLabel(const char* tag) {
  
  // If there is no tag provided, make one up
  HybridLabel hlbl;
  if (tag == nullptr) {
    const ulint obj_count = gbl_mem_balance_sheet.getActiveEntryCount();
    snprintf(hlbl.name, 23, "array%lu", obj_count);
  }
  else {
    if (strlen(tag) > 22) {
      std::string msg(tag);
      rtErr("Hybrid struct name " + msg + " is too long.", "assignLabel");
    }
    strcpy(hlbl.name, tag);
  }
  switch (format) {
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
    hlbl.format = expedited_code;
    break;
  case HybridFormat::DECOUPLED:
    hlbl.format = decoupled_code;
    break;
  case HybridFormat::UNIFIED:
    hlbl.format = unified_code;
    break;
  case HybridFormat::HOST_ONLY:
    hlbl.format = host_only_code;
    break;
  case HybridFormat::DEVICE_ONLY:
    hlbl.format = devc_only_code;
    break;
  case HybridFormat::HOST_MOUNTED:
    hlbl.format = host_mounted_code;
    break;
#else
  case HybridFormat::HOST_ONLY:
    hlbl.format = host_mounted_code;
    break;
#endif
  }

  // Scan available tags to decide on a serial number.
  hlbl.serial_number = gbl_mem_balance_sheet.getSerialNumber();

  return hlbl;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const Hybrid<T> setPointer(const Hybrid<T> *target, const size_t offset, const size_t length,
                           const char* output_name) {
  Hybrid<T> result(HybridKind::POINTER,
                   (output_name == nullptr) ? target->getLabel().name : output_name);
  result.setPointer(const_cast<Hybrid<T>*>(target), offset, length);
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
Hybrid<T> setPointer(Hybrid<T> *target, const size_t offset, const size_t length,
                     const char* output_name) {
  Hybrid<T> result(HybridKind::POINTER,
                   (output_name == nullptr) ? target->getLabel().name : output_name);
  result.setPointer(target, offset, length);
  return result;
}

} // namespace card
} // namespace stormm
