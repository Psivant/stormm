// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace card {

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> deepCopy(const void* origin, const size_t element_count,
                        const HybridTargetLevel origin_tier, const char* desc) {
  std::vector<T> result(element_count);
  deepCopy((void*)(result.data()), origin, sizeof(T), element_count, HybridTargetLevel::HOST,
           origin_tier, desc);
  return result;
}

} // namespace card
} // namespace stormm
