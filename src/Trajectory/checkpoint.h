// -*-c++-*-
#ifndef STORMM_CHECKPOINT_H
#define STORMM_CHECKPOINT_H

#include <string>
#include "copyright.h"
#include "Synthesis/phasespace_synthesis.h"
#include "phasespace.h"

namespace stormm {
namespace trajectory {

using synthesis::PhaseSpaceSynthesis;

/// \brief Write a checkpoint file for a set of coordinates and velocities.  The fixed precision
///        model for coordinates as well as velocities will be encoded in the file header.  Systems
///        that use double-precision (float64_t) coordinates will have their information encoded
///        as int64_t data, with the appropriate flag checked to indicate that the conversion
///        should be reversed when reading the file.  STORMM  uses a binary checkpoint file format
///        that seeks to conserve file size in various ways.  Each checkpoint file applies to one
///        and only one system.  A collection of systems can therefore be split up and run in
///        different ways from segment to segment.
///
/// Overloaded:
///   - Provide the coordinates and leading atom random number generator as const pointers
///   - Provide the coordinates and leading atom random number generator by const reference
///   - Provide a single system or multiple systems in a synthesis
///
/// \param ps            System coordinates
/// \param poly_ps       Synthesis of system coordinates
/// \param system_index  Index of the system from poly_ps to take
/// \param file_name     Name of the file to write
/// \param xrs           The random number generator associated with the first atom of the system.
///                      Its state will be read and recorded in the checkpoint file.  Upon restart,
///                      random number generators for other atoms in the system will be inferred
///                      by the generator's innate long jump function.
/// \{
template <typename Trng>
void writeCheckpointFile(const PhaseSpace *ps, const std::string &file_name, const Trng *xrs);

template <typename Trng>
void writeCheckpointFile(const PhaseSpace &ps, const std::string &file_name, const Trng &xrs);

template <typename Trng>
void writeCheckpointFile(const PhaseSpaceSynthesis *poly_ps, int system_index,
                         const std::string &file_name, const Trng *xrs);

template <typename Trng>
void writeCheckpointFile(const PhaseSpaceSynthesis &poly_ps, int system_index,
                         const std::string &file_name, const Trng &xrs);
/// \}

} // namespace trajectory
} // namespace stormm

#include "checkpoint.tpp"

#endif
