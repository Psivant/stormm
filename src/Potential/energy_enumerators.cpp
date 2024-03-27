#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "energy_enumerators.h"

namespace stormm {
namespace energy {

using constants::CaseSensitivity;
using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const EvaluateForce input) {
  switch (input) {
  case EvaluateForce::NO:
    return std::string("NO");
  case EvaluateForce::YES:
    return std::string("YES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const EvaluateEnergy input) {
  switch (input) {
  case EvaluateEnergy::NO:
    return std::string("NO");
  case EvaluateEnergy::YES:
    return std::string("YES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const EvaluateVirial input) {
  switch (input) {
  case EvaluateVirial::NO:
    return std::string("NO");
  case EvaluateVirial::YES:
    return std::string("YES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const DihedralStyle input) {
  switch (input) {
  case DihedralStyle::COSINE:
    return std::string("COSINE");
  case DihedralStyle::HARMONIC:
    return std::string("HARMONIC");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const StateVariable input) {
  switch (input) {
  case StateVariable::BOND:
    return std::string("BOND");
  case StateVariable::ANGLE:
    return std::string("ANGLE");
  case StateVariable::PROPER_DIHEDRAL:
    return std::string("PROPER_DIHEDRAL");
  case StateVariable::IMPROPER_DIHEDRAL:
    return std::string("IMPROPER_DIHEDRAL");
  case StateVariable::UREY_BRADLEY:
    return std::string("UREY_BRADLEY");
  case StateVariable::CHARMM_IMPROPER:
    return std::string("CHARMM_IMPROPER");
  case StateVariable::CMAP:
    return std::string("CMAP");
  case StateVariable::VDW:
    return std::string("VDW");
  case StateVariable::VDW_ONE_FOUR:
    return std::string("VDW_ONE_FOUR");
  case StateVariable::ELECTROSTATIC:
    return std::string("ELECTROSTATIC");
  case StateVariable::ELEC_ONE_FOUR:
    return std::string("ELEC_ONE_FOUR");
  case StateVariable::GENERALIZED_BORN:
    return std::string("GENERALIZED_BORN");
  case StateVariable::RESTRAINT:
    return std::string("RESTRAINT");
  case StateVariable::KINETIC:
    return std::string("KINETIC");
  case StateVariable::PRESSURE:
    return std::string("PRESSURE");
  case StateVariable::VIRIAL_11:
    return std::string("VIRIAL_11");
  case StateVariable::VIRIAL_12:
    return std::string("VIRIAL_12");
  case StateVariable::VIRIAL_22:
    return std::string("VIRIAL_22");
  case StateVariable::VIRIAL_13:
    return std::string("VIRIAL_13");
  case StateVariable::VIRIAL_23:
    return std::string("VIRIAL_23");
  case StateVariable::VIRIAL_33:
    return std::string("VIRIAL_33");
  case StateVariable::VOLUME:
    return std::string("VOLUME");
  case StateVariable::TEMPERATURE_ALL:
    return std::string("TEMPERATURE_ALL");
  case StateVariable::TEMPERATURE_PROTEIN:
    return std::string("TEMPERATURE_PROTEIN");
  case StateVariable::TEMPERATURE_LIGAND:
    return std::string("TEMPERATURE_LIGAND");
  case StateVariable::TEMPERATURE_SOLVENT:
    return std::string("TEMPERATURE_SOLVENT");
  case StateVariable::DU_DLAMBDA:
    return std::string("DU_DLAMBDA");
  case StateVariable::POTENTIAL_ENERGY:
    return std::string("POTENTIAL_ENERGY");
  case StateVariable::TOTAL_ENERGY:
    return std::string("TOTAL_ENERGY");
  case StateVariable::ALL_STATES:
    return std::string("ALL_STATES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const NonbondedPotential input) {
  switch (input) {
  case NonbondedPotential::ELECTROSTATIC:
    return std::string("ELECTROSTATIC");
  case NonbondedPotential::VAN_DER_WAALS:
    return std::string("VAN_DER_WAALS");
  case NonbondedPotential::CLASH:
    return std::string("CLASH");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const NonbondedTheme input) {
  switch (input) {
  case NonbondedTheme::ELECTROSTATIC:
    return std::string("ELECTROSTATIC");
  case NonbondedTheme::VAN_DER_WAALS:
    return std::string("VAN_DER_WAALS");
  case NonbondedTheme::ALL:
    return std::string("ALL");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const DecomposablePotential input) {
  switch (input) {
  case DecomposablePotential::ELECTROSTATIC:
    return std::string("ELECTROSTATIC");
  case DecomposablePotential::DISPERSION:
    return std::string("DISPERSION");
  case DecomposablePotential::ELEC_PME_DIRECT:
    return std::string("ELEC_PME_DIRECT");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const VdwCombiningRule input) {
  switch (input) {
  case VdwCombiningRule::LORENTZ_BERTHELOT:
    return std::string("LORENTZ_BERTHELOT");
  case VdwCombiningRule::GEOMETRIC:
    return std::string("GEOMETRIC");
  case VdwCombiningRule::NBFIX:
    return std::string("NBFIX");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ClashResponse input) {
  switch (input) {
  case ClashResponse::NONE:
    return std::string("NONE");
  case ClashResponse::FORGIVE:
    return std::string("FORGIVE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const EnergySample input) {
  switch (input) {
  case EnergySample::TIME_SERIES:
    return std::string("TIME_SERIES");
  case EnergySample::FINAL:
    return std::string("FINAL");
  case EnergySample::TIME_AVERAGE:
    return std::string("TIME_AVERAGE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const SplineScaffold input) {
  switch (input) {
  case SplineScaffold::LINEAR:
    return std::string("LINEAR");
  case SplineScaffold::SIGMOIDAL:
    return std::string("SIGMOIDAL");
  case SplineScaffold::QUADRATIC_HILL:
    return std::string("QUADRATIC_HILL");
  case SplineScaffold::QUADRATIC_VALLEY:
    return std::string("QUADRATIC_VALLEY");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const CellGridAction input) {
  switch (input) {
  case CellGridAction::INIT_FORCES:
    return std::string("INIT_FORCES");
  case CellGridAction::XFER_FORCES:
    return std::string("XFER_FORCES");
  case CellGridAction::UPDATE_IMG_COORD:
    return std::string("UPDATE_IMG_COORD");
  case CellGridAction::UPDATE_IMG_CELLS:
    return std::string("UPDATE_IMG_CELLS");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const QMapMethod input) {
  switch (input) {
  case QMapMethod::ACC_SHARED:
    return std::string("ACC_SHARED");
  case QMapMethod::GENERAL_PURPOSE:
    return std::string("GENERAL_PURPOSE");
  case QMapMethod::AUTOMATIC:
    return std::string("AUTOMATIC");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const PMIStrategy input) {
  switch (input) {
  case PMIStrategy::RECOMMENDED:
    return std::string("RECOMMENDED");
  case PMIStrategy::TIGHT:
    return std::string("TIGHT");
  case PMIStrategy::VERY_TIGHT:
    return std::string("VERY_TIGHT");
  case PMIStrategy::RECOMMENDED_PM_HEAVY:
    return std::string("RECOMMENDED_PM_HEAVY");
  case PMIStrategy::RECOMMENDED_PP_HEAVY:
    return std::string("RECOMMENDED_PP_HEAVY");
  case PMIStrategy::TIGHT_PM_HEAVY:
    return std::string("TIGHT_PM_HEAVY");
  case PMIStrategy::TIGHT_PP_HEAVY:
    return std::string("TIGHT_PP_HEAVY");
  case PMIStrategy::NO_AUTOMATION:
    return std::string("NO_AUTOMATION");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ValenceKernelSize input) {
  switch (input) {
  case ValenceKernelSize::XL:
    return std::string("XL");
  case ValenceKernelSize::LG:
    return std::string("LG");
  case ValenceKernelSize::MD:
    return std::string("MD");
  case ValenceKernelSize::SM:
    return std::string("SM");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
NonbondedPotential translateNonbondedPotential(const std::string &input) {
  if (strcmpCased(input, std::string("electrostatic"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("elec"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("charge-charge"), CaseSensitivity::NO)) {
    return NonbondedPotential::ELECTROSTATIC;
  }
  else if (strcmpCased(input, std::string("vdw"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("lennard-jones"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("van-der waals"), CaseSensitivity::NO)) {
    return NonbondedPotential::VAN_DER_WAALS;
  }
  else if (strcmpCased(input, std::string("clash"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("occlusion"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("exclusion"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("bump"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("steric"), CaseSensitivity::NO)) {
    return NonbondedPotential::CLASH;
  }
  else {
    rtErr("The input \"" + input + "\" does not have a valid enumeration.",
          "translateNonbondedPotential");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
NonbondedTheme translateNonbondedTheme(const std::string &input) {
  if (strcmpCased(input, std::string("electrostatic"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("elec"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("charge-charge"), CaseSensitivity::NO)) {
    return NonbondedTheme::ELECTROSTATIC;
  }
  else if (strcmpCased(input, std::string("vdw"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("lennard-jones"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("van-der waals"), CaseSensitivity::NO)) {
    return NonbondedTheme::VAN_DER_WAALS;
  }
  else if (strcmpCased(input, std::string("all"), CaseSensitivity::NO)) {
    return NonbondedTheme::ALL;
  }
  else {
    rtErr("The input \"" + input + "\" does not have a valid enumeration.",
          "translateNonbondedTheme");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
EnergySample translateEnergySample(const std::string &input) {
  if (strcmpCased(input, std::string("time_series"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("series"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("all"), CaseSensitivity::NO)) {
    return EnergySample::TIME_SERIES;
  }
  else if (strcmpCased(input, std::string("mean"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("average"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("timeaverage"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("time_average"), CaseSensitivity::NO)) {
    return EnergySample::TIME_AVERAGE;
  }
  else if (strcmpCased(input, std::string("final"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("last"), CaseSensitivity::NO)) {
    return EnergySample::FINAL;
  }
  else {
    rtErr("The input \"" + input + "\" does not have a valid enumeration.",
          "translateEnergySample");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
PMIStrategy translatePMIStrategy(const std::string &input) {
  if (strcmpCased(input, std::string("recommended"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("standard"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("default"), CaseSensitivity::NO)) {
    return PMIStrategy::RECOMMENDED;
  }
  else if (strcmpCased(input, std::string("tight"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("high"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("fine"), CaseSensitivity::NO)) {
    return PMIStrategy::TIGHT;
  }
  else if (strcmpCased(input, std::string("verytight"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("veryhigh"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("veryfine"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("very_tight"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("very_high"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("very_fine"), CaseSensitivity::NO)) {
    return PMIStrategy::VERY_TIGHT;
  }
  else if (strcmpCased(input, std::string("gridheavy"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("pmheavy"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("grid_heavy"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("pm_heavy"), CaseSensitivity::NO)) {
    return PMIStrategy::RECOMMENDED_PM_HEAVY;
  }
  else if (strcmpCased(input, std::string("pairheavy"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("ppheavy"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("pair_heavy"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("pp_heavy"), CaseSensitivity::NO)) {
    return PMIStrategy::RECOMMENDED_PP_HEAVY;
  }
  else if (strcmpCased(input, std::string("tight_gridheavy"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("tight_pmheavy"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("tight_grid_heavy"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("tight_pm_heavy"), CaseSensitivity::NO)) {
    return PMIStrategy::TIGHT_PM_HEAVY;
  }
  else if (strcmpCased(input, std::string("tight_pairheavy"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("tight_ppheavy"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("tight_pair_heavy"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("tight_pp_heavy"), CaseSensitivity::NO)) {
    return PMIStrategy::TIGHT_PP_HEAVY;
  }
  else if (strcmpCased(input, std::string("no_auto"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("none"), CaseSensitivity::NO)) {
    return PMIStrategy::NO_AUTOMATION;
  }
  else {
    rtErr("The input \"" + input + "\" does not have a valid enumeration.",
          "translatePMIStrategy");
  }
  __builtin_unreachable();
}

} // namespace energy
} // namespace stormm
