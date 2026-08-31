// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace data_types {

//-------------------------------------------------------------------------------------------------
template <typename T> bool isScalarType() {
  return (std::is_same_v<T, int>                    || std::is_same_v<T, unsigned int> ||
          std::is_same_v<T, double>                 || std::is_same_v<T, float>  ||
          std::is_same_v<T, long double>            || std::is_same_v<T, char> ||
          std::is_same_v<T, unsigned char>          || std::is_same_v<T, long long int> ||
          std::is_same_v<T, unsigned long long int> || std::is_same_v<T, short int> ||
          std::is_same_v<T, unsigned short int>     || std::is_same_v<T, unsigned long int> ||
          std::is_same_v<T, bool>                   || std::is_same_v<T, size_t>);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool isSignedIntegralScalarType() {
  return (std::is_same_v<T, int>           || std::is_same_v<T, char> ||
          std::is_same_v<T, long long int> || std::is_same_v<T, short int>);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool isUnsignedIntegralScalarType() {
  return (std::is_same_v<T, unsigned int>       || std::is_same_v<T, unsigned long int> ||
          std::is_same_v<T, unsigned char>      || std::is_same_v<T, unsigned long long int> ||
          std::is_same_v<T, unsigned short int> || std::is_same_v<T, bool> ||
          std::is_same_v<T, size_t>);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool isFloatingPointScalarType() {
  return (std::is_same_v<T, double> || std::is_same_v<T, float> || std::is_same_v<T, long double>);
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::string getStormmScalarTypeName() {
  if (std::is_same_v<T, int>) return "int";
  else if (std::is_same_v<T, unsigned int>) return "unsigned_int";
  else if (std::is_same_v<T, unsigned long int>) return "unsigned_long_int";
  else if (std::is_same_v<T, long double>) return "long double";
  else if (std::is_same_v<T, double>) return "double";
  else if (std::is_same_v<T, float>) return "float";
  else if (std::is_same_v<T, char>) return "char";
  else if (std::is_same_v<T, unsigned char>) return "unsigned_char";
  else if (std::is_same_v<T, long long int>) return "long_long_int";
  else if (std::is_same_v<T, unsigned long long int>) return "unsigned_long_long_int";
  else if (std::is_same_v<T, short int>) return "short_int";
  else if (std::is_same_v<T, unsigned short int>) return "unsigned_short_int";
  else if (std::is_same_v<T, bool>) return "bool";
  else if (std::is_same_v<T, size_t>) return "size_t";
  else {
    rtErr("Data type " + std::string(std::type_index(typeid(T)).name()) + " is not a recognized "
          "scalar type.", "getStormmScalarTypeName");
  }
  __builtin_unreachable();
}

} // namespace data_types
} // namespace stormm
