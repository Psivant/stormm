// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace data_types {

//-------------------------------------------------------------------------------------------------
template <typename T> bool isScalarType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  return (ct == int_type_index || ct == uint_type_index || ct == double_type_index ||
          ct == float_type_index || ct == longdouble_type_index || ct == char_type_index ||
          ct == uchar_type_index || ct == llint_type_index || ct == ullint_type_index ||
          ct == short_type_index || ct == ushort_type_index || ct == ulint_type_index ||
          ct == bool_type_index || ct == size_t_type_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool isSignedIntegralScalarType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  return (ct == int_type_index || ct == char_type_index || ct == llint_type_index ||
          ct == short_type_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool isUnsignedIntegralScalarType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  return (ct == uint_type_index || ct == ulint_type_index || ct == uchar_type_index ||
          ct == ullint_type_index || ct == ushort_type_index || ct == bool_type_index ||
          ct == size_t_type_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool isFloatingPointScalarType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  return (ct == double_type_index || ct == float_type_index || ct == longdouble_type_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::string getStormmScalarTypeName() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == int_type_index) return "int";
  else if (ct == uint_type_index) return "unsigned_int";
  else if (ct == ulint_type_index) return "unsigned_long_int";
  else if (ct == longdouble_type_index) return "long double";
  else if (ct == double_type_index) return "double";
  else if (ct == float_type_index) return "float";
  else if (ct == char_type_index) return "char";
  else if (ct == uchar_type_index) return "unsigned_char";
  else if (ct == llint_type_index) return "long_long_int";
  else if (ct == ullint_type_index) return "unsigned_long_long_int";
  else if (ct == short_type_index) return "short_int";
  else if (ct == ushort_type_index) return "unsigned_short_int";
  else if (ct == bool_type_index) return "bool";
  else if (ct == size_t_type_index) return "size_t";
  else {
    rtErr("Data type " + std::string(std::type_index(typeid(T)).name()) + " is not a recognized "
          "scalar type.", "getStormmScalarTypeName");
  }
  __builtin_unreachable();
}

} // namespace data_types
} // namespace stormm
