#include "copyright.h"
#include "casting_ops.h"

namespace stormm {
namespace data_types {

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph*> constCastVector(const std::vector<const AtomGraph*> &va) {
  std::vector<AtomGraph*> result;
  const size_t nva = va.size();
  result.reserve(nva);
  for (size_t i = 0; i < nva; i++) {
    result.push_back(const_cast<AtomGraph*>(va[i]));
  }
  return result;
}

} // namespace data_types
} // namespace stormm
