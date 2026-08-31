// -*-c++-*-
#ifndef STORMM_DEBUGGING_ENUMERATORS_H
#define STORMM_DEBUGGING_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace debug {

/// \brief List the places in which a large force, high speed, or unusual position might be
///        observed.
enum class AnomalyContext {
  PHASESPACE_SYNTHESIS,   ///< The anomalous event was recorded in a coordinate synthesis
  CELLGRID,               ///< The anomalous event was recorded in a neighbor list object
                          ///<   comprising all types of particle
  QQ_CELLGRID,            ///< The anomalous event was recorded in a neighbor list object
                          ///<   comprising particles with nonzero partial charges
  LJ_CELLGRID,            ///< The anomalous event was recorded in a neighbor list object
                          ///<   comprising particles with van-der Waals (Lennard-Jones) properties
};

/// \brief Get a human-readable name for the enumerations detailed above.
///
/// \param input  The enumeration of interest
/// \{
std::string getEnumerationName(AnomalyContext input);
/// \}

} // namespace debug
} // namespace stormm

#endif
