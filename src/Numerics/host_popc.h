// -*-c++-*-
#ifndef STORMM_HOSE_POPC_H
#define STORMM_HOSE_POPC_H

#include "copyright.h"
#include "DataTypes/common_types.h"

namespace stormm {
namespace numerics {

/// \brief Count the number of active bits in a short unsigned integer.
///
/// \param x  The integer to analyze
int hostPopcs(ushort x);

/// \brief Count the number of active bits in a short unsigned integer.
///
/// \param x  The integer to analyze
int hostPopc(uint x);

/// \brief Count the number of active bits in a short unsigned integer.
///
/// \param x  The integer to analyze
int hostPopcll(ullint x);

} // namespace numerics
} // namespace stormm

#endif
