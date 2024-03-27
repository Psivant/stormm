#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "restraint_enumerators.h"

namespace stormm {
namespace restraints {

using constants::CaseSensitivity;
using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RestraintEnsemble input) {
  switch (input) {
  case RestraintEnsemble::SPECIFIC_ATOMS:
    return std::string("SPECIFIC_ATOMS");
  case RestraintEnsemble::PREVENT_HBONDS:
    return std::string("PREVENT_HBONDS");
  case RestraintEnsemble::PRESERVE_HEAVY_DIHEDRALS:
    return std::string("PRESERVE_HEAVY_DIHEDRALS");
  case RestraintEnsemble::PRESERVE_POSITIONS:
    return std::string("PRESERVE_POSITIONS");
  case RestraintEnsemble::PRESERVE_DISTANCES:
    return std::string("PRESERVE_DISTANCES");
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
RestraintEnsemble translateRestraintEnsemble(const std::string &rst_group) {
  if (strcmpCased(rst_group, "hydrogen_bonds", CaseSensitivity::NO) ||
      strcmpCased(rst_group, "hydrogenbonds", CaseSensitivity::NO) ||
      strcmpCased(rst_group, "hbonds", CaseSensitivity::NO) ||
      strcmpCased(rst_group, "h-bonds", CaseSensitivity::NO) ||
      strcmpCased(rst_group, "prevent_hbonds", CaseSensitivity::NO) ||
      strcmpCased(rst_group, "prevent_hydrogen_bonds", CaseSensitivity::NO) ||
      strcmpCased(rst_group, "preventhydrogenbonds", CaseSensitivity::NO)) {
    return RestraintEnsemble::PREVENT_HBONDS;
  }
  else if (strcmpCased(rst_group, "heavy_dihedrals", CaseSensitivity::NO) ||
           strcmpCased(rst_group, "heavydihedrals", CaseSensitivity::NO) ||
           strcmpCased(rst_group, "heavy_dihe", CaseSensitivity::NO) ||
           strcmpCased(rst_group, "heavydihe", CaseSensitivity::NO) ||
           strcmpCased(rst_group, "preserve_heavy_dihedrals", CaseSensitivity::NO) ||
           strcmpCased(rst_group, "preserveheavydihedrals", CaseSensitivity::NO) ||
           strcmpCased(rst_group, "hold_heavy_dihedrals", CaseSensitivity::NO) ||
           strcmpCased(rst_group, "holdheavydihedrals", CaseSensitivity::NO)) {
    return RestraintEnsemble::PRESERVE_HEAVY_DIHEDRALS;
  }
  else if (strcmpCased(rst_group, "heavy_distances", CaseSensitivity::NO) ||
           strcmpCased(rst_group, "heavydistances", CaseSensitivity::NO) ||
           strcmpCased(rst_group, "preserve_heavy_distances", CaseSensitivity::NO) ||
           strcmpCased(rst_group, "preserveheavydistances", CaseSensitivity::NO) ||
           strcmpCased(rst_group, "hold_heavy_distances", CaseSensitivity::NO) ||
           strcmpCased(rst_group, "holdheavydistances", CaseSensitivity::NO)) {
    return RestraintEnsemble::PRESERVE_DISTANCES;
  }
  else if (strcmpCased(rst_group, "positions")) {
    return RestraintEnsemble::PRESERVE_POSITIONS;
  }
  else {
    rtErr("Unknown restraint ensemble code \"" + rst_group + "\".", "translateRestraintEnsemble");
  }
  __builtin_unreachable();
}

} // namespace restraints
} // namespace stormm
