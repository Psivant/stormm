// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace parse {

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<PolyNumeric> polyNumericVector(const std::vector<T> &values) {

  // Pre-allocate the result
  std::vector<PolyNumeric> result;
  const size_t vsize = values.size();
  result.resize(vsize);
  
  // Test the type of input, and prevent this from operating on a forbidden type
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (isHpcVectorType<T>()) {
    if (ct == char4_type_index) {
      rtErr("A vector of char4 data was fed into the template function, whereas the compiler "
            "should have found the special overloaded case.", "polyNumericVector");
    }
    else {
      rtErr("Data type " + getStormmHpcVectorTypeName<T>() + " cannot be converted into "
            "PolyNumeric.  Use a scalar type or char4.", "polyNumericVector");
    }
  }
  else if (isScalarType<T>()) {
    if (isSignedIntegralScalarType<T>()) {
      if (ct == llint_type_index) {
        for (size_t i = 0; i < vsize; i++) {
          result[i].lli = values[i];
        }
      }
      else {
        for (size_t i = 0; i < vsize; i++) {
          result[i].i = values[i];
        }
      }
    }
    else if (isUnsignedIntegralScalarType<T>()) {
      if (ct == ullint_type_index) {
        for (size_t i = 0; i < vsize; i++) {
          result[i].ulli = values[i];
        }
      }
      else {
        for (size_t i = 0; i < vsize; i++) {
          result[i].ui = values[i];
        }
      }
    }
    else if (isFloatingPointScalarType<T>()) {
      for (size_t i = 0; i < vsize; i++) {
        result[i].d = values[i];
      }
    }
  }
  else {
    rtErr("Data type " + std::string(std::type_index(typeid(T)).name()) + " cannot be converted "
          "into PolyNumeric.  Use a scalar type or char4.", "polyNumericVector");
  }

  return result;
}

} // namespace parse
} // namespace stormm
