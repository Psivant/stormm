// -*-c++-*-
#ifndef STORMM_STORMM_TYPE_NAME_H
#define STORMM_STORMM_TYPE_NAME_H

#include <string>
#include "copyright.h"
#include "common_types.h"
#include "stormm_vector_types.h"

namespace stormm {
namespace data_types {

/// \brief Make an exhaustive effort to assign a human-readable string to a given data type.  This
///        function will test the templated parameter against all tracked data types in order to
///        find a match, or return "unknown" if no such match could be found.
template <typename T> std::string getStormmTypeName();

} // namespace data_types
} // namespace stormm

#endif

