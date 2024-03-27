// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace card {

//-------------------------------------------------------------------------------------------------
template <typename T> cudaError_t wrapCudaFuncGetAttributes(cudaFuncAttributes *attrib, T ptr) {
  return cudaFuncGetAttributes(attrib, reinterpret_cast<const char*>(ptr));
}

} // namespace card
} // namespace stormm
