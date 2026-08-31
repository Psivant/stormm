#include "copyright.h"
#include "common_types.h"

namespace stormm {
namespace data_types {

//-------------------------------------------------------------------------------------------------
std::string getStormmScalarTypeName(const size_t index_code) {
  if (index_code == int_type_index) {
    return getStormmScalarTypeName<int>();
  }
  else if (index_code == double_type_index) {
    return getStormmScalarTypeName<double>();
  }
  else if (index_code == longdouble_type_index) {
    return getStormmScalarTypeName<long double>();
  }
  else if (index_code == float_type_index) {
    return getStormmScalarTypeName<float>();
  }
  else if (index_code == char_type_index) {
    return getStormmScalarTypeName<char>();
  }
  else if (index_code == uchar_type_index) {
    return getStormmScalarTypeName<unsigned char>();
  }
  else if (index_code == uint_type_index) {
    return getStormmScalarTypeName<unsigned int>();
  }
  else if (index_code == ulint_type_index) {
    return getStormmScalarTypeName<unsigned long int>();
  }
  else if (index_code == llint_type_index) {
    return getStormmScalarTypeName<long long int>();
  }
  else if (index_code == ullint_type_index) {
    return getStormmScalarTypeName<unsigned long long int>();
  }
  else if (index_code == short_type_index) {
    return getStormmScalarTypeName<short>();
  }
  else if (index_code == ushort_type_index) {
    return getStormmScalarTypeName<unsigned short>();
  }
  else if (index_code == bool_type_index) {
    return getStormmScalarTypeName<bool>();
  }
  else if (index_code == size_t_type_index) {
    return getStormmScalarTypeName<size_t>();
  }
  else {
    return std::string("unknown");
  }
  __builtin_unreachable();
}

} // namespace data_types
} // namespace stormm
