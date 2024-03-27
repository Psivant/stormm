// -*-c++-*-
#ifndef STORMM_HYBRID_H
#define STORMM_HYBRID_H

#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <cassert>
#include <typeinfo>
#include <typeindex>
#include <vector>
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
#include <cuda_runtime.h>
#  endif
#endif
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/rounding.h"
#include "Parsing/polynumeric.h"
#include "Reporting/error_format.h"
#include "gpu_enumerators.h"

namespace stormm {
namespace card {

using constants::ExceptionResponse;
using constants::mega;
using errors::rtErr;
using stmath::roundUp;
using stmath::getPaddedMemorySize;
using parse::PolyNumeric;

/// \brief Symbols for each of the Hybrid array formats
/// \{
constexpr char expedited_code = 'X';
constexpr char decoupled_code = 'L';
constexpr char unified_code   = 'U';
constexpr char devc_only_code = 'D';
constexpr char host_only_code = 'H';
/// \}

/// Symbols for laying out or freeing memory
constexpr llint allocating = 1;
constexpr llint deallocating = -1;
constexpr ulint hybrid_byte_increment = 128;
#ifdef STORMM_USE_HPC
constexpr HybridFormat default_hpc_format = HybridFormat::EXPEDITED;
#else
constexpr HybridFormat default_hpc_format = HybridFormat::HOST_ONLY;
#endif

/// \brief Byte lengths of specific data types, or components thereof
constexpr ulint tag_name_length = 23;

/// \brief A Hybrid object's immutable identifying information, stored as a const member variable
///        within the Hybrid object itself.
struct HybridLabel {

  /// A descriptive name for the Hybrid object
  char name[tag_name_length];

  /// Shorthand for the Hybrid object's format
  char format;

  /// Serial number of the Hybrid object in the global ledger, never to be changed so long as the
  /// Hybrid object exists
  int serial_number;
};

/// \brief Entries of the memory ledger recording not just a copy of the immutable type data but
///        also allocation sizes and allocation states.
struct LedgerEntry {
  HybridKind kind;
  HybridLabel label;
  size_t length;
  size_t element_size;
  int allocations;
  bool active;
};

/// \brief A struct for storing records of all active Hybrid memory arrays.  This will also store
///        tables of identifiers that indicate the most recent state of a group of Hybrid arrays,
///        which can be accessed by structs of pointers to make sure that they reflect the most
///        recent state of all their underlying arrays.
class Ledger {
public:

  /// \brief Constructor initializes the number of entries and allocations
  Ledger();

  /// \brief Default destructor
  ~Ledger() = default;

  /// \brief Retrieve the number of entries for Hybrid memory objects currently stored.
  int getActiveEntryCount() const;

  /// \brief Produce a number corresponding to some unused entry in the ledger.  This will be
  ///        either the lowest positive integer not yet assigned, or if some integers once assigned
  ///        to destroyed Hybrid objects have been freed up, the first such integer from that list.
  int getSerialNumber();

  /// \brief Retrieve the total allocation of decoupled memory, in bytes
  llint getTotalExpedited() const;

  /// \brief Retrieve the total allocation of decoupled memory, in bytes
  llint getTotalDecoupled() const;

  /// \brief Retrieve the total allocation of unified memory, in bytes
  llint getTotalUnified() const;

  /// \brief Retrieve the total allocation of host-only memory, in bytes
  llint getTotalHostOnly() const;

  /// \brief Retrieve the total allocation of device-only memory, in bytes
  llint getTotalDevcOnly() const;

  /// \brief Begin a catalog for the place of a Hybrid object among many others.  This is to be
  ///        called once for each object, in its constructor.
  ///
  /// \param source  The Hybrid object to be catalogged
  void setEntry(int serno, HybridKind kind,  const HybridLabel hlbl, size_t length,
                size_t element_size, int allocations);

  /// \brief Remove a Hybrid object from the catalog.  This is to be called by the Hybrid object's
  ///        destructor.
  ///
  /// \param source  The Hybrid object to be catalogged
  void unsetEntry(int serno);

  /// \brief Get information about a single entry of a ledger
  ///
  /// Overloaded:
  ///   - Produce the entry from a specific index
  ///   - Produce all entries corresponding to a particular label name
  /// \{
  LedgerEntry getEntry(int index) const;
  std::vector<LedgerEntry> getEntry(const std::string &name) const;
  /// \}

  /// \brief Log each tagged Hybrid array in the ledger
  void logMemory(size_t capacity, size_t element_size, HybridFormat fmt, const HybridLabel hlbl,
                 llint multiplier);

  /// \brief Query an index in the ledger to determine whether it catalogs an active Hybrid object.
  ///
  /// \param serno   Serial number of the Hybrid object in question
  bool testActive(int serno);

  /// \brief Query the number of times the Hybrid object catalogged in a particular entry has had
  ///        its data arrays allocated.
  int getAllocations(int serno);

  /// \brief Print the current memory usage of all Hybrid objects.  The function will print totals
  ///        followed by specific usages.
  ///
  /// \param n_display          Display at least this many usages.  Show all Hybrid objects if this
  ///                           value is less than zero or greater than the number of unique
  ///                           arrays.
  /// \param display_threshold  In addition, print any Hybrid objects allocated to more than this
  ///                           value on either the host or device
  void printMemoryProfile(int n_display = 8, llint display_threshold = mega);

private:

  /// Total memory allocated in each of four types
  /// \{
  llint total_expedited;
  llint total_decoupled;
  llint total_unified;
  llint total_host_only;
  llint total_devc_only;
  /// \}

  /// A list of Hybrid array entries, tracking memory used and allocation states
  std::vector<LedgerEntry> entries;

  /// A list of available Hybrid array serial numbers
  std::vector<int> free_slots;
};

/// \brief An evolution of GpuBuffer in pmemd.cuda, the Composite array has elements that are
///        accessible from either the GPU or the CPU.  In unified mode, the two data are one
///        insofar as the programmer sees it, and the page migration engine engaged by
///        cudaMallocManaged handles the two physical memory spaces at a very low level.  In other
///        modes, composite memory must be uploaded or downloaded explicitly in order to maintain
///        synchrony between host and device memory spaces.  This struct mimics a lot of the most
///        noteworthy behavior of std::vector, but with the potential to behave as a pointer as
///        well.  In this manner, it is intended to confer the convenience of C++ with the freedom
///        of classic C programming.
template <typename T> class Hybrid {
public:

  /// \brief Constructors vary based on whether CUDA is part of the compilation.
  ///
  /// Overloaded:
  ///   - Simple constructor takes an optional size and an optional format
  ///   - Secondary constructor takes a std::vector of the intended type and allocates based on its
  ///     size, again with an optional format
  ///
  /// \param length_in   The input length of the array
  /// \param tag_in      A human-readable tag by which to call this array
  /// \param format_in   Format of the resulting hybrid data structure
  /// \param kind_in     The kind of Hybrid object this will be
  /// \{
  explicit Hybrid(size_t length_in = 0, const char* tag_in = nullptr,
                  const HybridFormat format_in = default_hpc_format,
                  const HybridKind kind_in = HybridKind::ARRAY);

  explicit Hybrid(const std::vector<T> &S_in, const char* tag_in = nullptr,
                  const HybridFormat format_in = default_hpc_format);

  explicit Hybrid(const HybridKind kind_in, const char* tag_in = nullptr,
                  const HybridFormat format_in = default_hpc_format, size_t length_in = 0);
  /// \}

  /// Destructor frees data with no re-allocation
  ~Hybrid();

  /// \brief Copy constructor handles the reassignment of the underlying raw pointers
  ///
  /// \param original  The Hybrid object to copy
  Hybrid(const Hybrid &original);

  /// \brief The move constructor handles migration of a Hybrid object.
  ///
  /// \param original  The Hybrid object to move
  Hybrid(Hybrid &&original);  
  
  /// \brief Copy assignment constructor handles the reassignment of the underlying raw pointers
  ///
  /// \param other     The Hybrid object to copy (a different name for a better semantic fit in
  ///                  the context of the = sign)
  Hybrid& operator=(const Hybrid &other);

  /// \brief The move assignment operator must likewise handle transfer of the underlying data.
  ///
  /// \param other     The Hybrid object to move (a different name for a better semantic fit in
  ///                  the context of the = sign)
  Hybrid& operator=(Hybrid &&other);
  
  /// \brief Get the object kind (pointer or array)
  HybridKind getKind() const;

  /// \brief Get the memory format
  HybridFormat getFormat() const;

  /// \brief Return the number of elements in the Hybrid object's data array(s)
  size_t size() const;

  /// \brief Return the number of elements in the Hybrid object's data array(s)
  size_t capacity() const;

  /// \brief Return the pointer start position, targeting another ARRAY-kind Hybrid object
  ///        (relevant for POINTER-kind Hybrid objects only)
  size_t getPointerIndex() const;

  /// \brief Return the size of an individual element in the data array(s)
  size_t getElementSize() const;

  /// \brief Return the object's identifying label
  HybridLabel getLabel() const;

  /// \brief Return the object's identifying serial number
  int getSerialNumber() const;

  /// \brief Return the number of times this object has been allocated
  int getAllocations() const;
  
  /// \brief Produce the serial number of the target (only valid for POINTER-kind objects).  This
  ///        function is needed for the copy constructor in the POINTER-kind case.
  int getTargetSerialNumber() const;

  /// \brief Verify that the target of a HybridKind::POINTER object remains in the state it was
  ///        in when the pointer was first set.
  bool verifyTarget() const;

  /// \brief Get a pointer directly to the GpuArray's data on either the host or the device
  ///
  /// Overloaded:
  ///   - Get a const pointer if the object is const
  ///   - Get a pointer that can modify the underlying data if the object is non-const
  /// 
  /// \param tier  The level at which to set the pointer
  /// \{
  const T* data(const HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  T* data(const HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Get data from a given index of the host_data array.  This substitutes for direct array
  ///        access via the [ ] operator, as the host_data and devc_data pointers are private.
  ///
  /// Overloaded:
  ///   - Read a single element
  ///   - Read a subset of the elements, starting at an arbitrary point in the host_data array
  ///   - Read the entire host_data array, up to its length, into a new (out-of-place) array
  ///   - Read data into a pre-allocated, trusted array
  ///
  /// \param index   Index of the host_data array to read
  /// \param offset  Starting index of the host_data array to read
  /// \param count   The number of elements to access and return
  /// \{
  T readHost(size_t index) const;
  std::vector<T> readHost(size_t offset, size_t count) const;
  std::vector<T> readHost() const;
  void readHost(T* v, size_t offset, size_t count) const;
  /// \}

  /// \brief Put data into the host_data array.  This substitutes for direct array access via the
  ///        [ ] operator, as the host_data and devc_data pointers are private.
  ///
  /// Overloaded:
  ///   - Put a single element into the host_data at a specific index
  ///   - Put an array of elements into the host_data, starting at a specific offset index
  ///   - Target this (POINTER-kind) Hybrid object to another (ARRAY-kind) Hybrid object at a
  ///     specific offset, then load this object's host_data with an array of elements and pad
  ///     additional space in the target with zeros.  Return the next unpadded index in the target
  ///     for convenience.
  ///
  /// \param value   The data to write into host_data
  /// \param values  The data to write into host_data
  /// \param target  Another Hybrid object to target
  /// \param index   The index of host_data to write
  /// \param offset  The index at which to begin writing values into host_data 
  /// \param count   If the data comes as an array, write only the first count values into the
  ///                Hybrid object's host_data.  The default of zero implies that the entirety of
  ///                values shall be written.
  /// \{
  void putHost(const T value, size_t index);
  void putHost(const std::vector<T> &values, size_t offset, size_t count);
  void putHost(const std::vector<T> &values);
  size_t putHost(Hybrid<T> *target, const std::vector<T> &values, size_t offset = 0,
                 size_t padding = 0, size_t count = 0);
  /// \}
#ifdef STORMM_USE_HPC
  /// \brief Copy data from the host to the device, if necessary and if possible
  void upload(size_t start_pos = 0, size_t length_copy = 0);

  /// \brief Copy data from the device to the host, if necessary and if possible
  void download(size_t start_pos = 0, size_t length_copy = 0);

  /// \brief Copy the device-side data to be exported as a std::vector.  This operation will not
  ///        invoke the download() member function or modify the host-side data.
  ///
  /// Overloaded:
  ///   - Read and return data at a specific index
  ///   - Read and return an array of data, starting at a specific offset with a specific length
  ///   - Read and return the entire array of data from the device as a Standard Template Library
  ///     vector (out-of-place)
  ///   - Read data into a pre-allocated, trusted array
  ///
  /// \param index   Index of the device data to read
  /// \param offset  Index of the device data at which to begin reading
  /// \param count   The number of elements to read
  /// \param v       An array into which the results will be placed (this must be allocated with
  ///                enough space to hold the results)
  /// \{
  T readDevice(size_t index) const;
  std::vector<T> readDevice(size_t offset, size_t count) const;
  std::vector<T> readDevice() const;
  void readDevice(T* v, size_t offset, size_t count) const;
  /// \}

  /// \brief Put data into the devc_data array.  This allows direct host-to-device transfer of
  ///        data, bypassing the host_data array.
  ///
  /// Overloaded:
  ///   - Post data to a specific index
  ///   - Post part of a std::vector of data, starting at a specific offset of the Hybrid object's
  ///     device data with a specific length
  ///   - Post an entire std::vector of data to the start of the Hybrid object's device data
  ///
  /// \param value   The data to write into host_data
  /// \param values  The data to write into host_data
  /// \param index   Index of the device data to read
  /// \param offset  Index of the device data at which to begin reading
  /// \param count   The number of elements to read
  /// \{
  void putDevice(const T value, size_t index);
  void putDevice(const std::vector<T> &values, size_t offset, size_t count);
  void putDevice(const std::vector<T> &values);
  /// \}
#endif
  /// \brief Trim a Hybrid array to the exact size of its data
  void shrinkToFit();

  /// \brief Mimic the C++ std::vector push_back functionality.  A bounds check is followed by
  ///        extension of the data arrays in the current format if necessary.  The new element is
  ///        placed at the end of the host_data array.  The Hybrid object must be ARRAY-kind or a
  ///        POINTER-kind with sufficient maximum capacity (no reallocation needed).  Additions
  ///        will go to the CPU host memory.
  ///
  /// Overloaded:
  ///   - Add a single element
  ///   - Add multiple elements from a C-style array, Standard Template Library vector, or another
  ///     Hybrid object.
  ///
  /// \param element        The new item to add
  /// \param elements       An array of new elements to add
  /// \param element_count  The trusted length of elements, if adding from a C-style array
  /// \{
  void pushBack(const T element);
  void pushBack(const T* elements, const size_t element_count);
  void pushBack(const std::vector<T> &elements);
  void pushBack(const Hybrid<T> &elements);
  /// \}

  /// \brief A std::vector-like resize() feature.
  ///
  /// Overloaded:
  ///   - Resize and fill the new space with zero
  ///   - Resize and fill the new space with a constant value
  ///
  /// \param new_length  The new length.  If less than the current length, all that will change is
  ///                    the counted number of elements.  The original data will persist, as will
  ///                    the array(s) storing it.
  /// \param new_value   The value with which to populate new space.
  /// \{
  void resize(size_t new_length);
  void resize(size_t new_length, T value);
  /// \}

  /// \brief Set a Hybrid pointer struct to a segment of a Hybrid array struct.
  ///
  /// \param target      Pointer to the Hybrid object to point into.  This is a pointer because
  ///                    otherwise it would have to be a non-const reference.
  /// \param position    The location in target where this Hybrid pointer shall direct its own
  ///                    data array(s)
  /// \param new_length  The new, maximum length of data that this pointer shall occupy.  If
  ///                    left unset or set to a negative value, it will default to the remaining
  ///                    elements of the Hybrid array struct as it is currently found
  ///                    (target.size() - position).  The pointer's maximum capacity will likewise
  ///                    be set to new_length, or if new_length is not set it will default to
  ///                    (target.capacity() - position).
  void setPointer(Hybrid<T> *target, size_t position = 0, llint new_length = -1LL);

  /// \brief Swap the target of a POINTER-kind Hybrid object, transferring its starting index and
  ///        current length to apply instead to the new target.  Bounds checks will still be
  ///        applied.  Not valid for ARRAY-kind Hybrid objects or for POINTER-kind objects that
  ///        have no current target.
  ///
  /// \param target  Pointer to the new Hybrid object to point into.  This is a pointer because
  ///                otherwise it would have to be a non-const reference.
  /// \param policy  Indicate different actions to take if the Hybrid object does not yet point to
  ///                any target.
  void swapTarget(Hybrid<T> *new_target, ExceptionResponse policy = ExceptionResponse::SILENT);
  
  /// \brief Get a pointer to a Hybrid object at a specific location.  This is the way to get
  ///        const Hybrid POINTER-kind objects to const ARRAY-kind Hybrid objects.  This function
  ///        is only available to ARRAY-kind Hybrid objects.
  ///
  /// Overloaded:
  ///   - Get a const Hybrid POINTER-kind object to a const Hybrid array
  ///   - Get a non-const Hybrid POINTER-kind object to a non-const array
  ///
  /// \param position    The position in the array to which the POINTER-kind Hybrid object targets
  /// \param new_length  The length that the POINTER-kind Hybrid object will have (default 0
  ///                    implies it will be the length of the current array minus the POINTER-kind
  ///                    object's starting position)
  /// \{
  const Hybrid<T> getPointer(size_t position = 0, size_t new_length = 0,
                             const char* tag_in = nullptr) const;
  Hybrid<T> getPointer(size_t position = 0, size_t new_length = 0, const char* tag_in = nullptr);
  /// \}
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  /// \brief Get a device-readable and -writeable pointer to the host mapped data of a Hybrid
  ///        object.  For many GPUs, it is possible to use the host pointer for registered, pinned
  ///        memory, but the safest approach is to always get this pointer.
  ///
  /// Overloaded:
  ///   - Get a const Hybrid POINTER-kind object to a const Hybrid array
  ///   - Get a non-const Hybrid POINTER-kind object to a non-const array
  ///
  /// \{
  const T* getDeviceValidHostPointer() const;
  T* getDeviceValidHostPointer();
  /// \}
#  endif
#endif
private:

  /// Indicate whether this struct is a pointer to aspects of another Hybrid struct.  If so, it
  /// has no memory arrays associated with it, and cannot allocate, deallocate, or reallocate its
  /// memory.
  HybridKind kind;

  /// Indication of the design of this Hybrid array.  Each format has its uses: "_PINNED" memory is
  /// most convenient for transfers between host and device (in fact, it is the only format in
  /// which this process can happen directly--otherwise separate "pinned" memory must be allocated
  /// by cudaMemcpy), "decoupled" memory offers device memory in the fastest format for the GPU
  /// device to do computations on, "unified" memory is strong in both respects although not quite
  /// optimal in either, "devc_only" data will save host memory space, and "host_only" will serve
  /// as a placeholder when CUDA is not compiled.
  HybridFormat format;

  /// Length of the allocated arrays in either location
  size_t length;

  /// Size of an individual element
  size_t element_size;

  /// Smallest number of template type elements divisible by hybrid_byte_increment
  size_t growth_increment;

  /// Maximum number of storable elements in the array
  size_t max_capacity;

  /// Position to which the start of this object's data arrays point in another Hybrid object
  /// (relevant only for POINTER-kind objects)
  size_t pointer_index;

  /// Data array on the host.  This raw pointer is amenable to cudaHostAlloc() when pinned memory
  /// is desired.
  T* host_data;

  /// Data allocated on the device by cudaMalloc() or (in pinned fashion) with cudaHostAlloc()
#ifdef STORMM_USE_HPC
  T* devc_data;
#endif

  /// The immutable identifying information for this Hybrid object
  HybridLabel label;

  /// Number of times in which this object has had memory allocated
  int allocations;

  /// Serial number of the Hybrid array to which this object points (only valid for POINTER kind
  /// Hybrid objects)
  int target_serial_number;

  /// \brief Check the type of the object to be allocated.  Only POD types and HPC-tailored vector
  ///        types will be permitted.  It does not make sense to have a Hybrid<std::vector> as such
  ///        STL containers have no counterparts on the GPU.  Likewise, it does not make sense to
  ///        have struct MyStruct { float fx; int iy; double dz; }; Hybrid<MyStruct>.  Rather, the
  ///        GPU prefers MyStruct { Hybrid<float> fx; Hybrid<int> iy; Hybrid<double> dz }; for
  ///        optimal memory coalescence.  Trap use cases that violate these norms.
  void enforceDataTypeLimits();

  /// \brief Lay out the memory at each level of this object.
  void allocate();

  /// \brief Free memory associated with this object.
  void deallocate();

  /// \brief Re-allocate memory associated with this object.
  void reallocate(size_t new_length, size_t new_capacity = 0);

  /// \brief Function for creating a tag for a Hybrid struct.  If not provided in the constructor,
  ///        a tag will be created based on the number of previously created Hybrid objects.
  ///
  /// \param tag  The name for this Hybrid object
  HybridLabel assignLabel(const char* tag);
#if 0
  /// \brief Making the reinterpretHybridData function a friend of the Hybrid class does not
  ///        represent a risk of significant code bloat, any more than making an extra member
  ///        function for accomplishing the same thing.  Creating the free function provides a
  ///        clearer API, as to type {Hybrid Name}.template <template type>reinterpretCast() or
  ///        the like is an uncommon expression, particularly for novice C++ programmers.
  /// \{
  template <typename Treturn>
  friend const Hybrid<Treturn> reinterpretCast(const Hybrid *target, size_t offset, size_t length,
                                               const char* output_name);

  template <typename Treturn>
  friend Hybrid<Treturn> reinterpretCast(Hybrid *target, size_t offset, size_t length,
                                         const char* output_name);
  /// \}
#endif
};

/// \brief A free function can create a const POINTER-kind Hybrid set to target another const
///        Hybrid, even though a const constructor is not possible and a non-const Hybrid cannot
///        be set to target anything.  The const casting is handled internally.
///
/// Overloaded:
///   - Return a const POINTER-kind Hybrid targeting an existing const Hybrid of the same type.
///   - Return a non-const POINTER-kind Hybrid targeting a non-const Hybrid of the same type.
///     This can also be accomplished by using an existing Hybrid's setPointer() member function,
///     given compatible template types.
///
/// \param target       The existing Hybrid object to target
/// \param offset       The offset in the target array (this will be checked for validity by the
///                     Hybrid class's setPointer() member function)
/// \param length       The length of the resulting pointer array (also checked for validity)
/// \param output_name  New label to apply to the output Hybrid object (if unspecified, the name
///                     of the input Hybrid object will be taken instead)
/// \{
template <typename T>
const Hybrid<T> setPointer(const Hybrid<T> *target, size_t offset, size_t length,
                           const char* output_name = nullptr);

template <typename T>
Hybrid<T> setPointer(Hybrid<T> *target, size_t offset, size_t length,
                     const char* output_name = nullptr);
/// \}

#if 0
/// \brief Re-interpret the data in one ARRAY-kind Hybrid object as another data type using a
///        POINTER-kind Hybrid.  The data types must be of the same sizes.  Overloading and
///        descriptions of parameters follow from setPointer() above.
/// \{
template <typename Treturn, typename Ttarget>
const Hybrid<Treturn> reinterpretCast(const Hybrid<Ttarget> *target, size_t offset,
                                      size_t length, const char* output_name = nullptr);

template <typename Treturn, typename Ttarget>
Hybrid<Treturn> reinterpretCast(Hybrid<Ttarget> *target, size_t offset, size_t length,
                                const char* output_name = nullptr);
/// \}
#endif
} // namespace card
} // namespace stormm

/// \brief ***Global*** memory Ledger instance that Hybrid allocate / deallocate functions will
///        reference.  While this should be the only instance a programmer needs, it is feasible
///        to create additional Ledger objects for some other purpose. 
extern stormm::card::Ledger gbl_mem_balance_sheet;

#include "hybrid.tpp"

#endif
