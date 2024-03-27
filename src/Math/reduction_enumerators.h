// -*-c++-*-
#ifndef STORMM_REDUCTION_ENUMERATORS_H
#define STORMM_REDUCTION_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace stmath {

/// \brief Enumerate stages of a reduction operation, which can be taken piecemeal or all at once.
enum class ReductionStage {
  GATHER,     ///< Gather information from many particles or other sources.
  SCATTER,    ///< Scatter the accumulated information (this can also include a secondary gather
              ///<   on block-wide accumulated results).
  RESCALE,    ///< Rescale some value, essentially a secondary scattering process that can be done
              ///<   after some accumulation occurs during the first scattering process.
  ALL_REDUCE  ///< Perform the entire reduction in one step, when no intermediate accumulation is
              ///<   needed.
};

/// \brief Indicate the action to take when scattering the fully gathered sum (which is always a
///        sum) to the writeable array positions.  If further manipulations are needed, i.e. the
///        gathered sum is to be multiplied by two and then used as a divisor, those must be done
///        explicitly on the buffered values at the end of a gathering operation, so that the
///        scatter can then be applied in a subsequent step.
enum class ReductionGoal {
  NORMALIZE,          ///< Compute the norm of the data, whether the norm of the vector of scalar
                      ///<   values present at all atoms or the norm of the collection of forces
                      ///<   acting on each atom.
  CENTER_ON_ZERO,     ///< Compute the means of the data in up to three dimensions and subtract
                      ///<   these means from all values.
  CONJUGATE_GRADIENT  ///< Perform a conjugate gradient transformation of the forces acting on all
                      ///<   particles based on current and prior information.
};

/// \brief Indicate whether parallel reductions can be performed with one work unit per system.
enum class RdwuPerSystem {
  ONE,      ///< Each system's gather and scatter operations are performed by a single work unit,
            ///<   making one-kernel all-reduce possible.
  MULTIPLE  ///< One or more systems' gather and scatter operations are shared across multiple
            ///<   work units, forcing separate gather and scatter kernels.
};

/// \brief Enumerate the features of a reduction work unit abstract, as encoded in the synthesis
///        of topologies.
enum class RdwuAbstractMap {
  ATOM_START = 0, ///< Absolute starting point of atoms in various (concatenated, in the case of a
                  ///<   synthesis of many systems) arrays
  ATOM_END,       ///< Absolute upper limit of atoms to read or write in various arrays
  RESULT_INDEX,   ///< Location in one or more holding arrays for intermediate results of gathering
  DEPN_START,     ///< Start of dependencies in one or more holding arrays for making a final
                  ///<   assembly of the results prior to scattering
  DEPN_END,       ///< Upper limit of dependencies--if this is only one greater than the start, it
                  ///<   means that the reduction can be done by a single work unit
  SYSTEM_ID       ///< Index of the system to which this reduction pertains
};

/// \brief Produce a human-readable string to describe each enumeration.  Overloads of this
///        function are found in other libraries as well.
///
/// \param input  The enumerated value to translate
/// \{
std::string getEnumerationName(ReductionStage input);
std::string getEnumerationName(ReductionGoal input);
std::string getEnumerationName(RdwuPerSystem input);
std::string getEnumerationName(RdwuAbstractMap input);
/// \}
  
} // namespace stmath
} // namespace stormm

#endif
