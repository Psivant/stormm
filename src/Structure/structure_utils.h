// -*-c++-*-
#ifndef STORMM_STRUCTURE_UTILS_H
#define STORMM_STRUCTURE_UTILS_H

#include "copyright.h"

namespace stormm {
namespace structure {

/// \brief Check the supplied bounds for a partial scan or update of coordinates.
///
/// \param lower_limit  The lower limit of atoms to scan
/// \param upper_limit  The (proposed) upper limit of atoms to scan
/// \param natom        The number of atoms actually available in the coordinate set of interest
/// \param caller       Name of the calling function
void coordinateBoundsCheck(const int lower_limit, const int upper_limit, const int natom,
                           const char* caller = nullptr);

} // namespace structure
} // namespace stormm

#endif
