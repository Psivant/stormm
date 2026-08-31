// -*-c++-*-
#ifndef STORMM_PRESENT_DEBUG_H
#define STORMM_PRESENT_DEBUG_H

#include "copyright.h"
#include "MolecularMechanics/dynamics_intervention.h"
#include "Namelists/user_settings.h"

namespace stormm {
namespace review {

using mm::DynamicsIntervention;
using mm::dyna_tk;
using namelist::UserSettings;

/// \brief Collect the results of any and all analyses stored in a DynamicsIntervention class
///        object and print them to a file.
void createDebugReport(const UserSettings &ui, const DynamicsIntervention &mm_intv = dyna_tk);
 
} // namespace review
} // namespace stormm

#endif
