// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace parse {

//-------------------------------------------------------------------------------------------------
template <typename T>
int findAlignmentWidth(const T* numbers, const size_t count, int decimal_places, size_t advance) {
  int result = 0;
  if (isFloatingPointScalarType<T>() || isSignedIntegralScalarType<T>()) {
    for (size_t i = 0; i < count; i += advance) {
      if (fabs(numbers[i]) < 1.0) {
        result = std::max(result, 1 + (numbers[i] < 0.0));
      }
      else {
        result = std::max(result, static_cast<int>(ceil(log10(fabs(numbers[i]))) +
                                                   (numbers[i] < 0.0)));
      }
    }
  }
  else if (isUnsignedIntegralScalarType<T>()) {
    for (size_t i = 0; i < count; i += advance) {
      result = std::max(result, static_cast<int>(ceil(log10(numbers[i]))));
    }
  }
  else {
    rtErr("A recognized scalar data type is required.", "findAlignmentWidth");
  }
  return result + decimal_places + (decimal_places > 0);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
int findAlignmentWidth(const std::vector<T> &numbers, int decimal_places, size_t advance) {
  return findAlignmentWidth<T>(numbers.data(), numbers.size(), decimal_places, advance);
}

} // namespace parse
} // namespace stormm
