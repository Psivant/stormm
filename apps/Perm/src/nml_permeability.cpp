#include "../../../src/copyright.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Reporting/error_format.h"
#include "nml_permeability.h"

namespace permeability {

using stormm::constants::CaseSensitivity;
using stormm::errors::rtErr;
using stormm::errors::rtWarn;
using stormm::namelist::DefaultIsObligatory;
using stormm::namelist::InputRepeats;
using stormm::namelist::InputStatus;
using stormm::namelist::KeyRequirement;
using stormm::namelist::NamelistElement;
using stormm::namelist::NamelistType;
using stormm::parse::NumberFormat;
using stormm::parse::realToString;

//-------------------------------------------------------------------------------------------------
PermeabilityControls::PermeabilityControls(const ExceptionResponse policy_in) :
    policy{policy_in}, train_model{false}, simulation_count{default_simulation_count},
    discard_time{default_discard_time}, high_dielectric{default_high_dielectric},
    low_dielectric{default_low_dielectric},
{}

//-------------------------------------------------------------------------------------------------
bool PermeabilityControls::trainModel() const {
  return train_model;
}

//-------------------------------------------------------------------------------------------------
int PermeabilityControls::getSimulationCount() const {
  return simulation_count;
}

//-------------------------------------------------------------------------------------------------
double PermeabilityControls::getDiscardTime() const {
  return discard_time;
}

//-------------------------------------------------------------------------------------------------
double PermeabilityControls::getHighDielectric() const {
  return high_dielectric;
}

//-------------------------------------------------------------------------------------------------
double PermeabilityControls::getLowDielectric() const {
  return low_dielectric;
}

} // namespace permeability
