// -*-c++-*-
#ifndef STORMM_MIXED_TYPES_H
#define STORMM_MIXED_TYPES_H

#include "copyright.h"

namespace stormm {
namespace data_types {

/// \brief A combination of integer and double data.  Like other mixed tuples in this library,
///        this will maintain the x, y, z, and w member variable naming conventions.
struct CombineIDp {
  int x;
  double y;
};

/// \brief a templated, combined type for tagging any data type with an associated integer count
template <typename T> struct ValueWithCounter {
  T value;
  int count;
};
  
} // namespace data_types
} // namespace stormm

namespace stormm {
using data_types::CombineIDp;
using data_types::ValueWithCounter;
}

#endif
