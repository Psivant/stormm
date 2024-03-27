// -*-c++-*-
#ifndef STORMM_CASTING_OPS_H
#define STORMM_CASTING_OPS_H

#include <vector>
#include "copyright.h"
#include "Topology/atomgraph.h"

namespace stormm {
namespace data_types {

using topology::AtomGraph;
  
/// \brief Cast all elements of a vector to non-const.  Respect the const-ness of the vector
///        itself: this moves the const-ness to a single qualifier on the container rather than
///        individual qualifiers on all elements.
///
/// Overloaded:
///   - Operate on AtomGraph pointers
///
/// \param va  The vector of const items to consolidate as non-const inside a const container
/// \{
const std::vector<AtomGraph*> constCastVector(const std::vector<const AtomGraph*> &va);
/// \}

} // namespace data_types
} // namespace stormm

#endif
