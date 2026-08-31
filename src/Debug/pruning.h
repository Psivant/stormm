// -*-c++-*-
#ifndef STORMM_PRUNING_H
#define STORMM_PRUNING_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Namelists/nml_debug.h"
#include "Synthesis/phasespace_synthesis.h"

namespace stormm {
namespace debug {

using namelist::DebugControls;
using synthesis::PhaseSpaceSynthesis;

/// \brief Loop over all structures using CPU resources, evaluating all bond and bond angle terms
///        one at a time to check for high strain.  Report the system and atom indices of any
///        suspicous terms.  Return a list of systems which had no excessive strain in their bonds
///        and bond angles.  This function works with CPU host resources and assumes that the
///        relevant coordinates and structures ar available there.
///
/// Overloaded:
///   - Provide the original synthesis objects by const pointer
///   - Provide the original synthesis objects by const reference
///
/// \param poly_ps           The original coordinate synthesis
/// \param dbgcon            Derived from information in a &debug namelist.  This contains settings
///                          for the maximum bond and angle strain.
/// \{
std::vector<int> checkBondAndAngleSanity(const PhaseSpaceSynthesis *poly_ps,
                                         const DebugControls &dbgcon);

std::vector<int> checkBondAndAngleSanity(const PhaseSpaceSynthesis &poly_ps,
                                         const DebugControls &dbgcon);
/// \}
  
} // namespace debug
} // namespace stormm

#endif
