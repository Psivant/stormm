#include "copyright.h"
#include "chemistry_enumerators.h"

namespace stormm {
namespace chemistry {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(ChiralOrientation input) {
  switch (input) {
  case ChiralOrientation::RECTUS:
    return std::string("RECTUS");
  case ChiralOrientation::SINISTER:
    return std::string("SINISTER");
  case ChiralOrientation::NONE:
    return std::string("NONE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(ChiralInversionProtocol input) {
  switch (input) {
  case ChiralInversionProtocol::ROTATE:
    return std::string("ROTATE");
  case ChiralInversionProtocol::REFLECT:
    return std::string("REFLECT");
  case ChiralInversionProtocol::DO_NOT_INVERT:
    return std::string("DO_NOT_INVERT");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(MapRotatableGroups input) {
  switch (input) {
  case MapRotatableGroups::YES:
    return std::string("YES");
  case MapRotatableGroups::NO:
    return std::string("NO");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(RotatorCriterion input) {
  switch (input) {
  case RotatorCriterion::COM_PROXIMITY:
    return std::string("COM_PROXIMITY");
  case RotatorCriterion::GROUP_SIZE:
    return std::string("GROUP_SIZE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(ConformationEdit input) {
  switch (input) {
  case ConformationEdit::BOND_ROTATION:
    return std::string("BOND_ROTATION");
  case ConformationEdit::CIS_TRANS_FLIP:
    return std::string("CIS_TRANS_FLIP");
  case ConformationEdit::CHIRAL_INVERSION:
    return std::string("CHIRAL_INVERSION");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(EquivalenceSwap input) {
  switch (input) {
  case EquivalenceSwap::FREE_FOR_ALL:
    return std::string("FREE_FOR_ALL");
  case EquivalenceSwap::ROTARY:
    return std::string("ROTARY");
  }
  __builtin_unreachable();
}

} // namespace chemistry
} // namespace stormm
