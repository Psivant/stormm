// -*-c++-*-
#ifndef STORMM_COMMON_TYPES_H
#define STORMM_COMMON_TYPES_H

#include <string>
#include <typeinfo>
#include <typeindex>
#include <sys/types.h>
#include "copyright.h"
#include "Reporting/error_format.h"

namespace stormm {
namespace data_types {

/// \brief Integral type casts
/// \{
typedef unsigned int uint;
typedef unsigned long int ulint;
typedef unsigned long int ulong;
typedef long long int llint;
typedef unsigned long long int ullint;
typedef unsigned short int ushort;
typedef unsigned char uchar;
/// \}

/// \brief Type indices for common POD types, from which Hybrid objects might be composed
/// \{
static const size_t int_type_index = std::type_index(typeid(int)).hash_code();
static const size_t longdouble_type_index = std::type_index(typeid(long double)).hash_code();
static const size_t double_type_index = std::type_index(typeid(double)).hash_code();
static const size_t float_type_index = std::type_index(typeid(float)).hash_code();
static const size_t char_type_index = std::type_index(typeid(char)).hash_code();
static const size_t uchar_type_index = std::type_index(typeid(uchar)).hash_code();
static const size_t uint_type_index = std::type_index(typeid(uint)).hash_code();
static const size_t ulint_type_index = std::type_index(typeid(ulint)).hash_code();
static const size_t llint_type_index = std::type_index(typeid(llint)).hash_code();
static const size_t ullint_type_index = std::type_index(typeid(ullint)).hash_code();
static const size_t short_type_index = std::type_index(typeid(short int)).hash_code();
static const size_t ushort_type_index = std::type_index(typeid(short unsigned int)).hash_code();
static const size_t bool_type_index = std::type_index(typeid(bool)).hash_code();
static const size_t size_t_type_index = std::type_index(typeid(size_t)).hash_code();
/// \}
  
/// \brief Test whether some data type is a recognized scalar.
template <typename T> bool isScalarType();

/// \brief Test whether some data type is a recognized (signed) integral scalar.
template <typename T> bool isSignedIntegralScalarType();

/// \brief Test whether some data type is a recognized (unsigned) integral scalar.
template <typename T> bool isUnsignedIntegralScalarType();

/// \brief Test whether some data type is a recongized floating point number.
template <typename T> bool isFloatingPointScalarType();

/// \brief Produce a platform-independent name by which to identify one of the scalar data types.
template <typename T> std::string getStormmScalarTypeName();

} // namespace data_types
} // namespace stormm

#include "common_types.tpp"

namespace stormm {
using data_types::ulint;
using data_types::ulong;
using data_types::llint;
using data_types::ullint;
using data_types::ushort;
using data_types::uchar;
using data_types::int_type_index;
using data_types::longdouble_type_index;
using data_types::double_type_index;
using data_types::float_type_index;
using data_types::char_type_index;
using data_types::ushort_type_index;
using data_types::uchar_type_index;
using data_types::ulint_type_index;
using data_types::llint_type_index;
using data_types::ullint_type_index;
using data_types::short_type_index;
using data_types::ushort_type_index;
using data_types::bool_type_index;
using data_types::size_t_type_index;
} // namespace stormm

#endif
