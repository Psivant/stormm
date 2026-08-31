// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace data_types {

//-------------------------------------------------------------------------------------------------
template <typename T> std::string getStormmTypeName() {
  if (isScalarType<T>()) {
    return getStormmScalarTypeName<T>();
  }
  else if (isHpcVectorType<T>()) {
    return getHpcVectorTypeName<T>();
  }
  else {
    return std::string("unknown");
  }
  __builtin_unreachable();
}

} // namespace data_types
} // namespace stormm
