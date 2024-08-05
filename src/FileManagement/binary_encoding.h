// -*-c++-*-
#ifndef STORMM_BINARY_ENCODING_H
#define STORMM_BINARY_ENCODING_H

#include <cstdint>
#include <string>
#include <vector>
#include "copyright.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"

namespace stormm {
namespace diskutil {

using card::Hybrid;
using card::HybridTargetLevel;
using synthesis::PhaseSpaceSynthesis;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateSeries;
using trajectory::PhaseSpace;

/// \brief Default descriptors of various coordinate objects for binary file encoding.
/// \{
const char default_binary_file_descriptor[] = "No description provided";
const char default_cf_descriptor[] = "CoordinateFrame with positions";
const char default_ps_descriptor[] = "PhaseSpace with positions and velocities";
const char default_cs_descriptor[] = "CoordinateSeries with positions";
const char default_polyps_descriptor[] = "PhaseSpaceSynthesis with positions";
/// \}
  
/// \brief
class BinaryFileContent {
public:

  /// \brief The constructor requires the type of data, its length, some way to make a pointer to
  ///        the data itself, and the description.  The lower and upper bounds of the bytes for the
  ///        actual content will be calculated and stored in the class object once the list of
  ///        binary file content is compiled.
  ///
  /// Overloaded:
  ///   - Provide a scalar item for content, the type of which will be detected to translate into
  ///     the type_code member variable
  ///   - Provide a void-casted pointer to an array of data for content, the type of which must be
  ///     indicated by a prior-generated hash code
  /// \{
  template <typename T>
  BinaryFileContent(T content_in, bool repeating_in, const std::string description_in,
                    HybridTargetLevel tier_in = HybridTargetLevel::HOST);

  BinaryFileContent(size_t ct_vptr_in, size_t length_in, const void* vptr_in, bool repeating_in,
                    const std::string &description_in,
                    HybridTargetLevel tier_in = HybridTargetLevel::HOST);
  /// \}

  /// \brief Set the lower bound of memory for this content within a binary file.
  ///
  /// \param lower_bound_in  The lower bound of bytes for this data entry in the binary file
  void setLowerBound(llint lower_bound_in);

  /// \brief Set the upper bound of memory for this content within a binary file.
  ///
  /// \param upper_bound_in  The upper bound of bytes for this data entry in the binary file
  void setUpperBound(llint upper_bound_in);
  
private:
  size_t ct_vptr;           ///< The output of the std::type_index(typeid(T)).hash_code(), where
                            ///<   T is the type of the original data.  For each such type index
                            ///<   hash code there is a unique type_code that can be generated.
  uint16_t type_code;       ///< The coded type of the data to be stored in the binary file.  This
                            ///<   is a bit-packed string, encoding in the lowest two bits
                            ///<   whether the data type is integral (00), floating-point
                            ///<   (01), or a character (10).  The third bit indicates whether the
                            ///<   type is signed (1) or unsigned (0) (this will always be 1 for
                            ///<   floating point types).  Bits 4 through 7 encode the number of
                            ///<   elements in the tuple: a non-tuple such as a standard float gets
                            ///<   the code 0001, a tuple of two values would get the code 0010.
                            ///<   The eighth bit indicates whether the content pertains to a
                            ///<   scalar value (0) or an array(1).  The ninth bit indicates
                            ///<   whether the content can be entered repeatedly.  The high seven
                            ///<   bits provide the total size of the data's primary element, in
                            ///<   bytes: a uint_16 would get the code 0000010, a tuple of four
                            ///<   doubles the code 0001000.  Much larger data types are possible,
                            ///<   in theory, as are tuples beyond four values, although no such
                            ///<   types are currently supported by the rest of the code base.
  size_t length;            ///< The total length of the associated data array
  std::string description;  ///< A description of the array
  void* vptr;               ///< Void-casted pointer to the data itself, set to nullptr in the case
                            ///<   that the content is a single value.
  HybridTargetLevel tier;   ///< Indicate whether the data to be written to the binary file is
                            ///<   present in the CPU host or GPU device RAM.
  llint lower_bound;        ///< The first byte in the binary file containing the actual data
                            ///<   laid out in this contents entry
  llint upper_bound;        ///< The upper bound of indexed bytes in the binary file containing the
                            ///<   actual data laid out in this contents entry
};

/// \brief The table of contents keeps a list of arrays and specific values as well as text-based
///        strings describing them.  The strings can come from Hybrid objects' labels, for example,
///        but can also accept arbitrary character strings.
class BinaryFileKey {
public:

  /// \brief The constructor can take one or more arrays as well as specific objects, for
  ///        convenience.
  ///
  /// Overloaded:
  ///   - Prepare a blank table of contents with just a descriptor for the overall file and a
  ///     notion of the origin of the information
  ///   - Prepare a table of contents for a single array.  More can be appended later, modifying
  ///     the table of contents, but the binary file must be written in order after the table of
  ///     contents has been decided.
  ///   - Provide a container (Hybrid objects or Standard Template Library vectors) or a C-style
  ///     arrays with a trusted length
  ///   - Prepare a table of contents for specific objects.
  ///
  /// \param content   An array or STORMM object with data to write to the binary file
  /// \param length    Trusted length of content, if provided as a C-style array
  /// \param desc      The description of the array initializing the table of contents
  /// \{
  BinaryFileKey(const std::string &description_in = std::string(""));
  
  template <typename T> BinaryFileKey(const T* content, size_t length,
                                      const std::string &description_in = std::string(""),
                                      HybridTargetLevel tier_in = HybridTargetLevel::HOST);

  template <typename T> BinaryFileKey(const std::vector<T> *content,
                                      const std::string &description_in = std::string(""));

  template <typename T> BinaryFileKey(const std::vector<T> &content,
                                      const std::string &description_in = std::string(""));

  template <typename T> BinaryFileKey(const Hybrid<T> *content,
                                      const std::string &description_in = std::string(""),
                                      HybridTargetLevel tier_in = HybridTargetLevel::HOST);

  template <typename T> BinaryFileKey(const Hybrid<T> &content,
                                      const std::string &description_in = std::string(""),
                                      HybridTargetLevel tier_in = HybridTargetLevel::HOST);

  BinaryFileKey(const CoordinateFrame *cf,
                const std::string &description_in = std::string(default_cf_descriptor),
                HybridTargetLevel tier_in = HybridTargetLevel::HOST);

  BinaryFileKey(const CoordinateFrame &cf,
                const std::string &description_in = std::string(default_cf_descriptor),
                HybridTargetLevel tier_in = HybridTargetLevel::HOST);

  BinaryFileKey(const PhaseSpace *ps,
                const std::string &description_in = std::string(default_ps_descriptor),
                HybridTargetLevel tier_in = HybridTargetLevel::HOST);

  BinaryFileKey(const PhaseSpace &ps,
                const std::string &description_in = std::string(default_ps_descriptor),
                HybridTargetLevel tier_in = HybridTargetLevel::HOST);

  template <typename T>
  BinaryFileKey(const CoordinateSeries<T> *content,
                const std::string &description_in = std::string(default_cs_descriptor),
                HybridTargetLevel tier_in = HybridTargetLevel::HOST);

  template <typename T>
  BinaryFileKey(const CoordinateSeries<T> &content,
                const std::string &description_in = std::string(default_cs_descriptor),
                HybridTargetLevel tier_in = HybridTargetLevel::HOST);

  BinaryFileKey(const PhaseSpaceSynthesis *content,
                const std::string &description_in = std::string(default_polyps_descriptor),
                HybridTargetLevel tier_in = HybridTargetLevel::HOST);

  BinaryFileKey(const PhaseSpaceSynthesis &content,
                const std::string &description_in = std::string(default_polyps_descriptor),
                HybridTargetLevel tier_in = HybridTargetLevel::HOST);
  /// \}

  /// \brief Add a new array to the current table of contents.  The array will be appended by
  ///        default, but specifying a position for the array among the list is possible.
  ///        Overloading and descriptions of input parameters follow from the constructor, above,
  ///        in addition to:
  ///
  /// \param vptr  The void-casted pointer of whatever data array in its raw form.  The pointer may
  ///              address memory on the CPU host or GPU device.
  /// \param 
  /// \{
  void addArray(const void *vptr, size_t length, size_t ct_vptr, bool repeating_in,
                const std::string &description_in = std::string(""),
                HybridTargetLevel tier_in = HybridTargetLevel::HOST);

  template <typename T>
  void addArray(const T* content, size_t length, bool repeating_in,
                const std::string &description_in,
                HybridTargetLevel tier_in = HybridTargetLevel::HOST);
  
  template <typename T>
  void addArray(const std::vector<T> *content, bool repeating_in,
                const std::string &description_in = std::string(""));

  template <typename T>
  void addArray(const std::vector<T> &content, bool repeating_in,
                const std::string &description_in = std::string(""));

  template <typename T>
  void addArray(const Hybrid<T> *content, bool repeating_in,
                const std::string &description_in = std::string(""),
                HybridTargetLevel tier_in = HybridTargetLevel::HOST);

  template <typename T>
  void addArray(const Hybrid<T> &content, bool repeating_in,
                const std::string &description_in = std::string(""),
                HybridTargetLevel tier_in = HybridTargetLevel::HOST);
  /// \}

  /// \brief Add a critical scalar constant to the table of contents.  Every added array gets its
  ///        own sizing constants (size of the element, number of elements), and the number of
  ///        arrays also gets automatic inclusion.  Constants like the time step size, simulation
  ///        time constant, or fixed-precision bit counts are also critical to interpreting a
  ///        binary file.  Only data scalar types may be included in the table of contents, as they
  ///        will be decomposed into a program-specific bit string describing their attributes.
  ///
  /// \param parm  The parameter to include in the table of contents
  /// \param desc  The description of the parameter
  template <typename T> void addScalar(T parm,
                                       const std::string &description_in = std::string(""));
  /// \}

  /// \brief Add a two-component tuple to the table of contents.  Descriptions of input parameters
  ///        follow from addScalar(), above.
  template <typename T2> void addTuple2(T2 parm,
                                        const std::string &description_in = std::string(""));

  /// \brief Add a three-component tuple to the table of contents.  Descriptions of input
  ///        parameters follow from addScalar(), above.
  template <typename T3> void addTuple3(T3 parm,
                                        const std::string &description_in = std::string(""));

  /// \brief Add a four-component tuple to the table of contents.  Descriptions of input
  ///        parameters follow from addScalar(), above.
  template <typename T4> void addTuple4(T4 parm,
                                        const std::string &description_in = std::string(""));

private:
  std::string description;                         ///< Human-readable string decribing the overall
                                                   ///<   contents of the file.  Various arrays and
                                                   ///<   critical values will also have their own
                                                   ///<   descriptions.
  std::vector<char> tiers;                         ///< Character encoding the origin of the data,
                                                   ///<   whether on the CPU host or GPU device.
                                                   ///<   This is an explicit codification of the
                                                   ///<   HybridTargetLevel enumeration for binary
                                                   ///<   files to ensure backward compatibility
                                                   ///<   should another "tier" ever be added.
                                                   ///<   Other enumerators likewise must be
                                                   ///<   translated into explicit values to ensure
                                                   ///<   backwards compatibility, akin to the way
                                                   ///<   that type indices must be codified for
                                                   ///<   cross-platform compatibility.
  std::vector<BinaryFileContent> singleton_items;  ///< A list of scalars held within the table of
                                                   ///<   contents and considered part of its
                                                   ///<   overall byte footprint in the binary file
  std::vector<BinaryFileContent> array_items;      ///< A list of array items to be held within the
                                                   ///<   binary file
  size_t stride;                                   ///< If the table of contents includes repeating
                                                   ///<   data (one or more items could be written
                                                   ///<   multiple times into the file, in a
                                                   ///<   repeating sequence), this number will be
                                                   ///<   nonzero and indicate the length of the
                                                   ///<   repeating components.
  size_t toc_bytes;                                ///< The number of bytes occupied by the table
                                                   ///<   of contents itself, calculated and stored
                                                   ///<   for convenience.
  size_t nonrepeating_bytes;                       ///< The length of non-repeating data, including
                                                   ///<   the table of contents.
  llint total_expected_bytes;                      ///< The total number of bytes expected to be
                                                   ///<   written to the file, acts in conjunction
                                                   ///<   with the final checksum.  If the total
                                                   ///<   number of bytes is not known in advance,
                                                   ///<   the total expected bytes will be recorded
                                                   ///<   as a negative 64-bit integer.
  llint total_written_bytes;                       ///< The total number of bytes committed to the
                                                   ///<   file.  This will be written to the end
                                                   ///<   of the file, one of three 64-bit integers
                                                   ///<   making up the final checksum.
  llint toc_checksum;                              ///< The sum of all elements of the table of
                                                   ///<   contents.  This will be written at the
                                                   ///<   end of all data, after a watermark.
  llint content_checksum;                          ///< The sum of all 16-bit elements of the
                                                   ///<   contents.  If any of the content contains
                                                   ///<   an odd number of bytes, the sum of the
                                                   ///<   final byte will be included in the
                                                   ///<   checksum.
  TrajectoryCompression packed_format;             ///< Method of compressing the file if it is
                                                   ///<   some type of trajectory
};

/// \brief Translate a HybridTargetLevel enumeration into a character for writing to the binary
///        file.
char codifyHybridTargetLevel(HybridTargetLevel tier);
  
/// \brief Translate a type index produced by the C++ Standard Library into a 16-bit code for
///        storage in a BinaryFileContent class object and, later, a binary file's table of
///        contents.
///
/// \param ct             The type index found by C++ Standard Library features
/// \param is_array       Indicate whether the data entry is a singleton value or an array
/// \param is_repeatable  Indicate whether the data entry can be entered repeatedly in the binary
///                       output.  If a binary file has repeating entries A, B, and E with
///                       non-repeating entries C, D, and F, then all non-repeating entries will be
///                       written at the head of the file (C, D, and F), followed by A1, B1, C1,
///                       A2, B2, C2, ..., An, Bn, Cn.
uint16_t codifyTypeIndex(size_t ct, bool is_array, bool is_repeatable);
  
} // namespace diskutil
} // namespace stormm

#include "binary_encoding.tpp"

#endif
