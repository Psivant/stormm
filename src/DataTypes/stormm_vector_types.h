// -*-c++-*-
#ifndef STORMM_VECTOR_TYPES_H
#define STORMM_VECTOR_TYPES_H

#include <complex>
#include <cstddef>
#include <typeinfo>
#include <typeindex>
#include <vector>
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
#    include <vector_types.h>
#    include <cufft.h>
#  endif
#endif
#include "copyright.h"
#include "Reporting/error_format.h"
#include "common_types.h"

namespace stormm {
namespace data_types {

#ifndef STORMM_USE_HPC
// If no HPC language is defined, then some of the CUDA POD vector types will need to be
// delineated in order to run the program in CPU code.
struct int2 {
  int x;
  int y;
};

struct int3 {
  int x;
  int y;
  int z;
};

struct int4 {
  int x;
  int y;
  int z;
  int w;
};

struct uint2 {
  unsigned int x;
  unsigned int y;
};

struct uint3 {
  unsigned int x;
  unsigned int y;
  unsigned int z;
};

struct uint4 {
  unsigned int x;
  unsigned int y;
  unsigned int z;
  unsigned int w;
};

struct float2 {
  float x;
  float y;
};

struct float3 {
  float x;
  float y;
  float z;
};

struct float4 {
  float x;
  float y;
  float z;
  float w;
};

struct longlong2 {
  long long int x;
  long long int y;
};

struct longlong3 {
  long long int x;
  long long int y;
  long long int z;
};

struct ulonglong2 {
  unsigned long long int x;
  unsigned long long int y;
};

struct ulonglong3 {
  unsigned long long int x;
  unsigned long long int y;
  unsigned long long int z;
};

struct double2 {
  double x;
  double y;
};

struct double3 {
  double x;
  double y;
  double z;
};

struct char2 {
  char x;
  char y;
};

struct char3 {
  char x;
  char y;
  char z;
};

struct char4 {
  char x;
  char y;
  char z;
  char w;
};

struct uchar2 {
  unsigned char x;
  unsigned char y;
};

struct uchar3 {
  unsigned char x;
  unsigned char y;
  unsigned char z;
};

struct uchar4 {
  unsigned char x;
  unsigned char y;
  unsigned char z;
  unsigned char w;
};

struct short2 {
  short int x;
  short int y;
};

struct short3 {
  short int x;
  short int y;
  short int z;
};

struct short4 {
  short int x;
  short int y;
  short int z;
  short int w;
};

struct ushort2 {
  unsigned short int x;
  unsigned short int y;
};

struct ushort3 {
  unsigned short int x;
  unsigned short int y;
  unsigned short int z;
};

struct ushort4 {
  unsigned short int x;
  unsigned short int y;
  unsigned short int z;
  unsigned short int w;
};

// Various 32-byte tuples require explicit alignment specifications in CUDA 13.0 and beyond.  To
// ensure compatibility between host code and GPU-compiled code, cognates of these types are
// defined whenever an HPC compiler is not is use / the HPC language's vector_types.h is not
// available.
#  if defined(__linux__) || defined(__unix) || defined(__unix__)
typedef struct double4_16a_base {
  double x;
  double y;
  double z;
  double w;
} __attribute__((aligned(16))) double4_16a;

typedef struct longlong4_16a_base {
  long long int x;
  long long int y;
  long long int z;
  long long int w;
} __attribute__((aligned(16))) longlong4_16a;

typedef struct ulonglong4_16a_base {
  unsigned long long int x;
  unsigned long long int y;
  unsigned long long int z;
  unsigned long long int w;
} __attribute__((aligned(16))) ulonglong4_16a;
#  elif defined(__APPLE__) || defined(__MACH__)
typedef struct double4_16a_base {
  double x;
  double y;
  double z;
  double w;
} __attribute__((aligned(16))) double4_16a;

typedef struct longlong4_16a_base {
  long long int x;
  long long int y;
  long long int z;
  long long int w;
} __attribute__((aligned(16))) longlong4_16a;

typedef struct ulonglong4_16a_base {
  unsigned long long int x;
  unsigned long long int y;
  unsigned long long int z;
  unsigned long long int w;
} __attribute__((aligned(16))) ulonglong4_16a;
#  elif defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
typedef __declspec(align(16)) struct double4_16a_base {
  double x;
  double y;
  double z;
  double w;
} double4_16a;

typedef __declspec(align(16)) struct longlong4_16a_base {
  long long int x;
  long long int y;
  long long int z;
  long long int w;
} longlong4_16a;

typedef __declspec(align(16)) struct ulonglong4_16a_base {
  unsigned long long int x;
  unsigned long long int y;
  unsigned long long int z;
  unsigned long long int w;
} ulonglong4_16a;
#  endif // Various operating system cases
#else // STORMM_USE_HPC

// Certain vector types which have been made the standard in CUDA 13.0 and beyond have had their
// nomenclature hard-wired into the code base (e.g. double4_16a, rather than double4).  In other
// HPC languages and even prior versions of CUDA, these aligned types may not be defined, but
// their primitive names are more likely to be available.  Revert to those names in such cases.
#  if (CUDART_VERSION < 13000)
typedef double4 double4_16a;
typedef longlong4 longlong4_16a;
typedef ulonglong4 ulonglong4_16a;
#  else 
#    ifndef __CUDACC__
#      if defined(__linux__) || defined(__unix) || defined(__unix__)
typedef struct double4_16a_base {
  double x;
  double y;
  double z;
  double w;
} __attribute__((aligned(16))) double4_16a;

typedef struct longlong4_16a_base {
  long long int x;
  long long int y;
  long long int z;
  long long int w;
} __attribute__((aligned(16))) longlong4_16a;

typedef struct ulonglong4_16a_base {
  unsigned long long int x;
  unsigned long long int y;
  unsigned long long int z;
  unsigned long long int w;
} __attribute__((aligned(16))) ulonglong4_16a;
#      elif defined(__APPLE__) || defined(__MACH__)
typedef struct double4_16a_base {
  double x;
  double y;
  double z;
  double w;
} __attribute__((aligned(16))) double4_16a;

typedef struct longlong4_16a_base {
  long long int x;
  long long int y;
  long long int z;
  long long int w;
} __attribute__((aligned(16))) longlong4_16a;

typedef struct ulonglong4_16a_base {
  unsigned long long int x;
  unsigned long long int y;
  unsigned long long int z;
  unsigned long long int w;
} __attribute__((aligned(16))) ulonglong4_16a;
#      elif defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
typedef __declspec(align(16)) struct double4_16a_base {
  double x;
  double y;
  double z;
  double w;
} double4_16a;

typedef __declspec(align(16)) struct longlong4_16a_base {
  long long int x;
  long long int y;
  long long int z;
  long long int w;
} longlong4_16a;

typedef __declspec(align(16)) struct ulonglong4_16a_base {
  unsigned long long int x;
  unsigned long long int y;
  unsigned long long int z;
  unsigned long long int w;
} ulonglong4_16a;
#      endif // Various operating system cases
#    endif // CUDA compiler is not in use
#  endif // CUDART_VERSION < 13000
#endif // STORMM_USE_HPC

// Type definitions for consistent nomenclature
typedef longlong2 llint2;
typedef longlong3 llint3;
typedef longlong4_16a llint4;
typedef ulonglong2 ullint2;
typedef ulonglong3 ullint3;
typedef ulonglong4_16a ullint4;

// CUDA- and HIP-defined vector types corresponding to each of the above POD types
static const size_t int2_type_index = std::type_index(typeid(int2)).hash_code();
static const size_t int3_type_index = std::type_index(typeid(int3)).hash_code();
static const size_t int4_type_index = std::type_index(typeid(int4)).hash_code();
static const size_t double2_type_index = std::type_index(typeid(double2)).hash_code();
static const size_t double3_type_index = std::type_index(typeid(double3)).hash_code();
static const size_t double4_type_index = std::type_index(typeid(double4_16a)).hash_code();
static const size_t float2_type_index = std::type_index(typeid(float2)).hash_code();
static const size_t float3_type_index = std::type_index(typeid(float3)).hash_code();
static const size_t float4_type_index = std::type_index(typeid(float4)).hash_code();
static const size_t char2_type_index = std::type_index(typeid(char2)).hash_code();
static const size_t char3_type_index = std::type_index(typeid(char3)).hash_code();
static const size_t char4_type_index = std::type_index(typeid(char4)).hash_code();
static const size_t uchar2_type_index = std::type_index(typeid(uchar2)).hash_code();
static const size_t uchar3_type_index = std::type_index(typeid(uchar3)).hash_code();
static const size_t uchar4_type_index = std::type_index(typeid(uchar4)).hash_code();
static const size_t uint2_type_index = std::type_index(typeid(uint2)).hash_code();
static const size_t uint3_type_index = std::type_index(typeid(uint3)).hash_code();
static const size_t uint4_type_index = std::type_index(typeid(uint4)).hash_code();
static const size_t longlong2_type_index = std::type_index(typeid(longlong2)).hash_code();
static const size_t longlong3_type_index = std::type_index(typeid(longlong3)).hash_code();
static const size_t longlong4_type_index = std::type_index(typeid(longlong4_16a)).hash_code();
static const size_t ulonglong2_type_index = std::type_index(typeid(ulonglong2)).hash_code();
static const size_t ulonglong3_type_index = std::type_index(typeid(ulonglong3)).hash_code();
static const size_t ulonglong4_type_index = std::type_index(typeid(ulonglong4_16a)).hash_code();
static const size_t short2_type_index = std::type_index(typeid(short2)).hash_code();
static const size_t short3_type_index = std::type_index(typeid(short3)).hash_code();
static const size_t short4_type_index = std::type_index(typeid(short4)).hash_code();
static const size_t ushort2_type_index = std::type_index(typeid(ushort2)).hash_code();
static const size_t ushort3_type_index = std::type_index(typeid(ushort3)).hash_code();
static const size_t ushort4_type_index = std::type_index(typeid(ushort4)).hash_code();
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
static const size_t cufftz_type_index = std::type_index(typeid(cufftDoubleComplex)).hash_code();
static const size_t cufftc_type_index = std::type_index(typeid(cufftComplex)).hash_code();
#  endif
#endif
static const size_t stdz_type_index = std::type_index(typeid(std::complex<double>)).hash_code();
static const size_t stdc_type_index = std::type_index(typeid(std::complex<float>)).hash_code();

/// \brief Test whether some data type is a recognized HPC (CUDA or HIP) vector type.
template <typename T> bool isHpcVectorType();

/// \brief Test whether some data types is a signed integral HPC vector type
template <typename T> bool isSignedIntegralHpcVectorType();

/// \brief Test whether some data types is an unsigned integral HPC vector type
template <typename T> bool isUnsignedIntegralHpcVectorType();

/// \brief Test whether some data types is a real number HPC vector type
template <typename T> bool isFloatingPointHpcVectorType();

/// \brief Get the size of an HPC vector, in terms of the number of member variables
template <typename T> int getHpcVectorTypeSize();

/// \brief Produce a platform-independent name by which to identify one of the scalar data types.
///
/// Overloaded:
///   - Provide the type as a template parameter
///   - Provide the type as an indexed code
///
/// \param index_code  The codified integer assigned to the data type at runtime
/// \{
template <typename T> std::string getHpcVectorTypeName();
std::string getHpcVectorTypeName(size_t index_code);
/// \}

/// \brief A mixed tuple for 95-bit integer accumulation.  This is the proper way to take
///        double-precision floating point arithmetic into fixed precision format with minimal
///        loss.
struct int95_t {
  long long int x;  ///< The primary accumulator.  If the fixed precision takes 55 bits after
                    ///<   the decimal, this accumulator is enough to capture all the
                    ///<   information in a 64-bit floating point number in the range +/-256.
  int y;            ///< Overflow accumulator, used when the primary value would be exceeded
};

// Also define the fixed-precision type for 95-bit storage of double-precision numbers.
static const size_t int95t_type_index = std::type_index(typeid(int95_t)).hash_code();

/// \brief A template-compatible class for implementing vectorized two-tuples in C++ code, based
///        on a single scalar type.
template <typename T> struct Vec2 {

  /// \brief Constructors will accept two scalars or one of the defined vector tuples.
  /// \{
  explicit Vec2(T x_in, T y_in);
  template <typename Tinput> Vec2(const Tinput v_in);
  /// \}

  // Public member variables
  T x;  ///< The first member, equivalent to the x member of a float2 or longlong2 tuple
  T y;  ///< The second member, equivalent to the y member of a float2 or longlong2 tuple
};

/// \brief A template-compatible class for implementing vectorized three-tuples in C++ code, based
///        on a single scalar type.
template <typename T> struct Vec3 {

  /// \brief Constructors will accept three scalars or one of the defined vector tuples.
  /// \{
  explicit Vec3(T x_in, T y_in, T z_in);
  template <typename Tinput> Vec3(const Tinput v_in);
  /// \}
  
  T x;  ///< The first member, equivalent to the x member of a double3 or uint3 tuple
  T y;  ///< The second member, equivalent to the y member of a double3 or uint3 tuple
  T z;  ///< The third member, equivalent to the z member of a double3 or uint3 tuple
};

/// \brief A template-compatible class for implementing vectorized four-tuples in C++ code, based
///        on a single scalar type.
template <typename T> struct Vec4 {

  /// \brief Constructors will accept four scalars or one of the defined vector tuples.
  /// \{
  explicit Vec4(T x_in, T y_in, T z_in, T w_in);
  template <typename Tinput> Vec4(const Tinput v_in);
  /// \}
  
  T x;  ///< The first member, equivalent to the x member of a float4 or int4 tuple
  T y;  ///< The second member, equivalent to the y member of a float4 or int4 tuple
  T z;  ///< The third member, equivalent to the z member of a float4 or int4 tuple
  T w;  ///< The fourth member, equivalent to the w member of a float4 or int4 tuple
};
  
/// \brief Provide options for selected vector type conversions, and vector forms thereof:
///        double(2,3,4) <--> float(2,3,4), double(2,3,4) <-- int(2,3,4),
///        double(2,3,4) <-- uint(2,3,4), double(2,3,4) <-- llint(2,3,4),
///        double(2,3,4) <-- ullint(2,3,4), float(2,3,4) <-- int(2,3,4),
///        float(2,3,4) <-- uint(2,3,4), float(2,3,4) <-- llint(2,3,4),
///        float(2,3,4) <-- ullint(2,3,4), double(2,3,4) <-- short(2,3,4),
///        double(2,3,4) <-- ushort(2,3,4), float(2,3,4) <-- short(2,3,4),
///        float(2,3,4) <-- ushort(2,3,4)
///
/// \param rhs  Right-hand argument, basis for the assignment
/// \{
template <typename T> double2 vtConv2(const T rhs);
template <typename T> double3 vtConv3(const T rhs);
template <typename T> double4_16a vtConv4(const T rhs);
template <typename T> float2 vtConv2f(const T rhs);
template <typename T> float3 vtConv3f(const T rhs);
template <typename T> float4 vtConv4f(const T rhs);
template <typename T> std::vector<double2> vtConv2(const std::vector<T> &rhs);
template <typename T> std::vector<double3> vtConv3(const std::vector<T> &rhs);
template <typename T> std::vector<double4_16a> vtConv4(const std::vector<T> &rhs);
template <typename T> std::vector<float2> vtConv2f(const std::vector<T> &rhs);
template <typename T> std::vector<float3> vtConv3f(const std::vector<T> &rhs);
template <typename T> std::vector<float4> vtConv4f(const std::vector<T> &rhs);
/// \}

/// \brief Provide a means to convert various 32-bit tuple types to unsigned 32-bit integers.
///
/// \param ctuple  The character tuple to convert
/// \{
uint short2ToUint(const short2 ituple);
uint ushort2ToUint(const ushort2 ituple);
uint char4ToUint(const char4 ctuple);
uint uchar4ToUint(const uchar4 ctuple);
/// \}

/// \brief Provide a means to convert various 32-bit tuple types to unsigned 32-bit integers.
///
/// \param ctuple  The character tuple to convert
/// \{
short2 uintToShort2(const uint val);
ushort2 uintToUshort2(const uint val);
char4 uintToChar4(const uint val);
uchar4 uintToUchar4(const uint val);
/// \}

/// \brief Operators for particular vector types.
///
/// \param lhs  Left-hand argument, basis for comparison
/// \param rhs  Right-hand argument, basis for comparison
/// \{
bool operator==(const short2 lhs, const short2 rhs);
bool operator!=(const short2 lhs, const short2 rhs);
bool operator==(const short3 lhs, const short3 rhs);
bool operator!=(const short3 lhs, const short3 rhs);
bool operator==(const short4 lhs, const short4 rhs);
bool operator!=(const short4 lhs, const short4 rhs);

bool operator==(const ushort2 lhs, const ushort2 rhs);
bool operator!=(const ushort2 lhs, const ushort2 rhs);
bool operator==(const ushort3 lhs, const ushort3 rhs);
bool operator!=(const ushort3 lhs, const ushort3 rhs);
bool operator==(const ushort4 lhs, const ushort4 rhs);
bool operator!=(const ushort4 lhs, const ushort4 rhs);

bool operator==(const int2 lhs, const int2 rhs);
bool operator!=(const int2 lhs, const int2 rhs);
bool operator==(const int3 lhs, const int3 rhs);
bool operator!=(const int3 lhs, const int3 rhs);
bool operator==(const int4 lhs, const int4 rhs);
bool operator!=(const int4 lhs, const int4 rhs);

bool operator==(const uint2 lhs, const uint2 rhs);
bool operator!=(const uint2 lhs, const uint2 rhs);
bool operator==(const uint3 lhs, const uint3 rhs);
bool operator!=(const uint3 lhs, const uint3 rhs);
bool operator==(const uint4 lhs, const uint4 rhs);
bool operator!=(const uint4 lhs, const uint4 rhs);

bool operator==(const llint2 lhs, const llint2 rhs);
bool operator!=(const llint2 lhs, const llint2 rhs);
bool operator==(const llint3 lhs, const llint3 rhs);
bool operator!=(const llint3 lhs, const llint3 rhs);
bool operator==(const llint4 lhs, const llint4 rhs);
bool operator!=(const llint4 lhs, const llint4 rhs);

bool operator==(const ullint2 lhs, const ullint2 rhs);
bool operator!=(const ullint2 lhs, const ullint2 rhs);
bool operator==(const ullint3 lhs, const ullint3 rhs);
bool operator!=(const ullint3 lhs, const ullint3 rhs);
bool operator==(const ullint4 lhs, const ullint4 rhs);
bool operator!=(const ullint4 lhs, const ullint4 rhs);

bool operator<(const char2 lhs, const char2 rhs);
bool operator<=(const char2 lhs, const char2 rhs);
bool operator>(const char2 lhs, const char2 rhs);
bool operator>=(const char2 lhs, const char2 rhs);
bool operator==(const char2 lhs, const char2 rhs);
bool operator!=(const char2 lhs, const char2 rhs);

bool operator<(const char3 lhs, const char3 rhs);
bool operator<=(const char3 lhs, const char3 rhs);
bool operator>(const char3 lhs, const char3 rhs);
bool operator>=(const char3 lhs, const char3 rhs);
bool operator==(const char3 lhs, const char3 rhs);
bool operator!=(const char3 lhs, const char3 rhs);

bool operator<(const char4 lhs, const char4 rhs);
bool operator<=(const char4 lhs, const char4 rhs);
bool operator>(const char4 lhs, const char4 rhs);
bool operator>=(const char4 lhs, const char4 rhs);
bool operator==(const char4 lhs, const char4 rhs);
bool operator!=(const char4 lhs, const char4 rhs);
/// \}

/// \brief Mimic critical CUDA typecasting intrinsics with simple unions
/// {
union Ecumenical2 {
  short si;
  ushort usi;
  char2 c2;
  uchar2 uc2;
};

union Ecumenical4 {
  int i;
  uint ui;
  float f;
  char4 c4;
  uchar4 uc4;
};

union Ecumenical8 {
  llint lli;
  ullint ulli;
  double d;
  int2 i2;
  uint2 ui2;
  float2 f2;
};
/// \}

} // namespace data_types
} // namespace stormm

#include "stormm_vector_types.tpp"

namespace stormm {

// Make the HPC vectorized tuple types available through STORMM namespaces, if they were defined
// without an overarching HPC library.
#ifndef STORMM_USE_HPC
using data_types::int2;
using data_types::int3;
using data_types::int4;
using data_types::double2;
using data_types::double3;
using data_types::double4_16a;
using data_types::float2;
using data_types::float3;
using data_types::float4;
using data_types::char2;
using data_types::char3;
using data_types::char4;
using data_types::uchar2;
using data_types::uchar3;
using data_types::uchar4;
using data_types::uint2;
using data_types::uint3;
using data_types::uint4;
using data_types::longlong2;
using data_types::longlong3;
using data_types::longlong4_16a;
using data_types::ulonglong2;
using data_types::ulonglong3;
using data_types::ulonglong4_16a;
using data_types::short2;
using data_types::short3;
using data_types::short4;
using data_types::ushort2;
using data_types::ushort3;
using data_types::ushort4;
#else // STORMM_USE_HPC
#  if (!defined(__CUDACC__)) || (CUDART_VERSION < 13000)
using data_types::double4_16a;
using data_types::longlong4_16a;
using data_types::ulonglong4_16a;
#  endif
#endif // STORMM_USE_HPC

// Alternate type definitions need to be carried through regardless of HPC compilation
using data_types::llint2;
using data_types::llint3;
using data_types::llint4;
using data_types::ullint2;
using data_types::ullint3;
using data_types::ullint4;

// Also make the int95_t data type generally available by including this header file.
using data_types::int95_t;
  
// Make type index identifiers available in all STORMM namespaces
using data_types::int2_type_index;
using data_types::int3_type_index;
using data_types::int4_type_index;
using data_types::double2_type_index;
using data_types::double3_type_index;
using data_types::double4_type_index;
using data_types::float2_type_index;
using data_types::float3_type_index;
using data_types::float4_type_index;
using data_types::char2_type_index;
using data_types::char3_type_index;
using data_types::char4_type_index;
using data_types::uchar2_type_index;
using data_types::uchar3_type_index;
using data_types::uchar4_type_index;
using data_types::uint2_type_index;
using data_types::uint3_type_index;
using data_types::uint4_type_index;
using data_types::longlong2_type_index;
using data_types::longlong3_type_index;
using data_types::longlong4_type_index;
using data_types::ulonglong2_type_index;
using data_types::ulonglong3_type_index;
using data_types::ulonglong4_type_index;
using data_types::short2_type_index;
using data_types::short3_type_index;
using data_types::short4_type_index;
using data_types::ushort2_type_index;
using data_types::ushort3_type_index;
using data_types::ushort4_type_index;
using data_types::int95t_type_index;
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
using data_types::cufftz_type_index;
using data_types::cufftc_type_index;
#  endif
#endif
using data_types::stdz_type_index;
using data_types::stdc_type_index;

// Template-compatible vector types for C++ code
using data_types::Vec2;
using data_types::Vec3;
using data_types::Vec4;

// Templated conversion functions for getting anything of the right tuple complexity into a float
// or double tuple
using data_types::vtConv2;
using data_types::vtConv3;
using data_types::vtConv4;
using data_types::vtConv2f;
using data_types::vtConv3f;
using data_types::vtConv4f;
using data_types::short2ToUint;
using data_types::ushort2ToUint;
using data_types::char4ToUint;
using data_types::uchar4ToUint;
using data_types::uintToShort2;
using data_types::uintToUshort2;
using data_types::uintToChar4;
using data_types::uintToUchar4;

/// Include the ecumenical types for byte-wise conversions (in absence of various HPC intrinsics)
using data_types::Ecumenical2;
using data_types::Ecumenical4;
using data_types::Ecumenical8;
  
// Overloaded assignment and comparison operators
using data_types::operator<;
using data_types::operator<=;
using data_types::operator>;
using data_types::operator>=;
using data_types::operator==;
using data_types::operator!=;

} // namespace stormm

#endif
