// -*-c++-*-
#ifndef STORMM_NUMERIC_ENUMERATORS_H
#define STORMM_NUMERIC_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace numerics {

/// \brief Enumerate the choices for carrying out fixed-precision accumulation
enum class AccumulationMethod {
  SPLIT,     ///< Use split accumulation, stashing the low 32 bits in a locally cached int and the
             ///<   high 32 bits in a secondary accumulator probably located further away in main
             ///<   memory.  So long as most of the work happens in the low 32 bits, this reduces
             ///<   local memory demand and overall memory bandwidth by a factor of two, lowers
             ///<   GPU kernel register pressure on many architectures, and has shown 2.2 - 2.9x
             ///<   the speed of accumulating in int64.
  WHOLE,     ///< Sum fixed-precision numbers in int64 accumulators.  This is needed when the
             ///<   fixed-precision work cannot be mostly confined to the low 32 bits.
  AUTOMATIC  ///< Determine the accumulation method by looking at the number of fixed-precision
             ///<   bits after the decimal and making some assumptions about typical molecular
             ///<   mechanics forces.
};

/// \brief Get a human-readable name for the enumerations detailed above.
///
/// \param input  The enumeration of interest
/// \{
std::string getEnumerationName(AccumulationMethod input);
/// \}

} // namespace numerics
} // namespace stormm

#endif
