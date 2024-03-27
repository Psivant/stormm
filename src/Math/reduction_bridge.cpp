#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "reduction_bridge.h"

namespace stormm {
namespace stmath {

using card::HybridKind;

//-------------------------------------------------------------------------------------------------
ReductionBridge::ReductionBridge(const size_t n_values) :
  x_buffer{HybridKind::POINTER, "bridge_xbuff"},
  y_buffer{HybridKind::POINTER, "bridge_ybuff"},
  z_buffer{HybridKind::POINTER, "bridge_zbuff"},
  storage{3LLU * roundUp(n_values, warp_size_zu)}
{
  const size_t padded_nval = roundUp(n_values, warp_size_zu);
  x_buffer.setPointer(&storage,                  0, n_values);
  y_buffer.setPointer(&storage,        padded_nval, n_values);
  z_buffer.setPointer(&storage, 2LLU * padded_nval, n_values);
}

//-------------------------------------------------------------------------------------------------
ReductionBridge::ReductionBridge(const ReductionBridge &original) :
  x_buffer{original.x_buffer},
  y_buffer{original.y_buffer},
  z_buffer{original.z_buffer},
  storage{original.storage}
{
  x_buffer.swapTarget(&storage);
  y_buffer.swapTarget(&storage);
  z_buffer.swapTarget(&storage);
}

//-------------------------------------------------------------------------------------------------
ReductionBridge& ReductionBridge::operator=(const ReductionBridge &other) {
  if (this == &other) {
    return *this;
  }
  x_buffer = other.x_buffer;
  y_buffer = other.y_buffer;
  z_buffer = other.z_buffer;
  storage = other.storage;
  x_buffer.swapTarget(&storage);
  y_buffer.swapTarget(&storage);
  z_buffer.swapTarget(&storage);
  return *this;
}

//-------------------------------------------------------------------------------------------------
size_t ReductionBridge::size() const {
  return x_buffer.size();
}
  
//-------------------------------------------------------------------------------------------------
const double* ReductionBridge::getPointer(const CartesianDimension cdim,
                                          const HybridTargetLevel tier) const {
  switch (cdim) {
  case CartesianDimension::X:
    return x_buffer.data(tier);
  case CartesianDimension::Y:
    return y_buffer.data(tier);
  case CartesianDimension::Z:
    return z_buffer.data(tier);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double* ReductionBridge::getPointer(const CartesianDimension cdim, const HybridTargetLevel tier) {
  switch (cdim) {
  case CartesianDimension::X:
    return x_buffer.data(tier);
  case CartesianDimension::Y:
    return y_buffer.data(tier);
  case CartesianDimension::Z:
    return z_buffer.data(tier);
  }
  __builtin_unreachable();
}

} // namespace stmath
} // namespace stormm
