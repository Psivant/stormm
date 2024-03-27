#include <string>
#include "copyright.h"
#include "forcefield_enumerators.h"

namespace stormm {
namespace modeling {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ParameterKind kind) {
  switch (kind) {
  case ParameterKind::BOND:
    return std::string("harmonic bond");
  case ParameterKind::ANGLE:
    return std::string("harmonic angle");
  case ParameterKind::DIHEDRAL:
    return std::string("cosine-based dihedral");
  case ParameterKind::UREY_BRADLEY:
    return std::string("Urey-Bradley stretching angle");
  case ParameterKind::CHARMM_IMPROPER:
    return std::string("CHARMM improper dihedral");
  case ParameterKind::CMAP:
    return std::string("CMAP 2D bicubic spline surface");
  case ParameterKind::ATTN_14_SCALE:
    return std::string("attenuated 1:4 interaction");
  case ParameterKind::CHARGE:
    return std::string("atomic partial charge");
  case ParameterKind::LENNARD_JONES:
    return std::string("Lennard-Jones potential");
  case ParameterKind::BUCKINGHAM:
    return std::string("Buckingham potential");
  case ParameterKind::VIRTUAL_SITE_FRAME:
    return std::string("virtual site frame");
  case ParameterKind::NONE:
    return std::string("no parameter");
  }
  __builtin_unreachable();
}
  
} // namespace modeling
} // namespace stormm
