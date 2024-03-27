// -*-c++-*-
#ifndef STORMM_BEHAVIOR_H
#define STORMM_BEHAVIOR_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace constants {

/// \brief Enumerate the ways in which the code can respond to an exception under a variety of
///        circumstances: die, warn (and continue), or stay silent and continue
enum class ExceptionResponse {
  DIE, WARN, SILENT
};

/// \brief Enumerate the degree of modification to make to user input that, while it may not meet
///        a strict standard, is nevertheless tolerable and salvageable.
enum class ModificationPolicy {
  MODIFY, DO_NOT_MODIFY
};
  
/// \brief Express choices for the case-sensitivity of various inputs
enum class CaseSensitivity {
  YES,      ///< Case matters
  NO,       ///< Evaluate without case sensitivity--convert everything to uppercase before parsing
  AUTOMATIC ///< Do not impose one style or another--defer to local, default behavior for
            ///<   individual inputs
};

/// \brief Enumerate precision models
enum class PrecisionModel {
  SINGLE,  ///< Evaluates most interactions in 32-bit integer and floating point arithmetic.
           ///< References single-precision data arrays if available.
  DOUBLE,  ///< Evaluates most interactions in 64-bit integer and floating point arithmetic.
           ///< References double-precision data arrays if available.
};

/// \brief Enumerate the Cartesian dimensions
enum class CartesianDimension {
  X = 0, Y, Z
};

/// \brief Enumerate the unit cell axes
enum class UnitCellAxis {
  A = 0, B, C
};

/// \brief Convert a string object, possibly supplied by the user, into an ExceptionResponse
///        enumeration to indicate the program's intended behavior.
///
/// \param policy  Human-readable indication of the exception response behavior
ExceptionResponse translateExceptionResponse(const std::string &policy);
  
/// \brief Translate the exception response behavior into a human-readable string.
///
/// \param choice  The named exception response behavior
std::string getEnumerationName(ExceptionResponse policy);
  
/// \brief Translate a string into a known precision level enumeration.
///
/// \param choice  The named precision model (will be checked for validity)
PrecisionModel translatePrecisionModel(const std::string &choice);

/// \brief Get a descriptive string corresponding to each enumerated compute precision model.
///
/// \param pmodel  The precision model to name
std::string getEnumerationName(PrecisionModel pmodel);

/// \brief Translate the name of a Cartesian axis.
///
/// \param axis  The axis of interest
std::string getEnumerationName(CartesianDimension axis);

/// \brief Translate the name of a unit cell axis.
///
/// \param axis  The axis of interest
std::string getEnumerationName(UnitCellAxis axis);

} // namespace constants
} // namespace stormm

#endif
