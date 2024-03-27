#include <string>
#include "copyright.h"
#include "Parsing/parse.h"
#include "molecule_format_enumerators.h"

namespace stormm {
namespace structure {

using data_types::operator==;
  
//-------------------------------------------------------------------------------------------------
MdlMolPropertyKind translateMdlMolPropertyKind(char4 input) {
  if      (input == char4({ ' ', ' ', ' ', 'A' })) return MdlMolPropertyKind::ATOM_ALIAS;
  else if (input == char4({ ' ', ' ', ' ', 'V' })) return MdlMolPropertyKind::ATOM_VALUE;
  else if (input == char4({ ' ', ' ', ' ', 'G' })) return MdlMolPropertyKind::GROUP_ABBREVIATION;
  else if (input == char4({ 'C', 'H', 'G', 'M' })) return MdlMolPropertyKind::CHARGE;
  else if (input == char4({ 'R', 'A', 'D', 'M' })) return MdlMolPropertyKind::RADICAL;
  else if (input == char4({ 'I', 'S', 'O', 'M' })) return MdlMolPropertyKind::ISOTOPE;
  else if (input == char4({ 'R', 'B', 'C', 'M' })) return MdlMolPropertyKind::RING_BOND_COUNT;
  else if (input == char4({ 'S', 'U', 'B', 'M' })) return MdlMolPropertyKind::SUBSTITUTION_COUNT;
  else if (input == char4({ 'U', 'N', 'S', 'M' })) return MdlMolPropertyKind::UNSATURATED_COUNT;
  else if (input == char4({ 'L', 'I', 'N', 'M' })) return MdlMolPropertyKind::LINK_ATOM;
  else if (input == char4({ 'A', 'L', 'S', 'M' })) return MdlMolPropertyKind::ATOM_LIST;
  else if (input == char4({ 'A', 'P', 'O', 'M' })) return MdlMolPropertyKind::ATTACHMENT_POINT;
  else if (input == char4({ 'A', 'A', 'L', 'M' })) return MdlMolPropertyKind::ATTACHMENT_ORDER;
  else if (input == char4({ 'R', 'G', 'P', 'M' })) {
    return MdlMolPropertyKind::RGROUP_LABEL_LOCATION;
  }
  else if (input == char4({ 'L', 'O', 'G', 'M' })) return MdlMolPropertyKind::RGROUP_LOGIC;
  else if (input == char4({ 'S', 'T', 'Y', 'M' })) return MdlMolPropertyKind::SGROUP_TYPE;
  else if (input == char4({ 'S', 'S', 'T', 'M' })) return MdlMolPropertyKind::SGROUP_SUBTYPE;
  else if (input == char4({ 'S', 'L', 'B', 'M' })) return MdlMolPropertyKind::SGROUP_LABELS;
  else if (input == char4({ 'S', 'C', 'N', 'M' })) return MdlMolPropertyKind::SGROUP_CONNECTIVITY;
  else if (input == char4({ 'S', 'D', 'S', 'M' })) return MdlMolPropertyKind::SGROUP_EXPANSION;
  else if (input == char4({ 'S', 'A', 'L', 'M' })) return MdlMolPropertyKind::SGROUP_ATOM_LIST;
  else if (input == char4({ 'S', 'B', 'L', 'M' })) return MdlMolPropertyKind::SGROUP_BOND_LIST;
  else if (input == char4({ 'S', 'P', 'A', 'M' })) return MdlMolPropertyKind::MG_PARENT_ATOM_LIST;
  else if (input == char4({ 'S', 'M', 'T', 'M' })) return MdlMolPropertyKind::SGROUP_SUBSCRIPT;
  else if (input == char4({ 'C', 'R', 'S', 'M' })) {
    return MdlMolPropertyKind::SGROUP_CORRESPONDENCE;
  }
  else if (input == char4({ 'S', 'D', 'I', 'M' })) return MdlMolPropertyKind::SGROUP_DISPLAY_INFO;
  else if (input == char4({ 'S', 'B', 'V', 'M' })) return MdlMolPropertyKind::SGROUP_BOND_VECTOR;
  else if (input == char4({ 'S', 'D', 'T', 'M' })) return MdlMolPropertyKind::SGROUP_FIELD;
  else if (input == char4({ 'S', 'D', 'D', 'M' })) return MdlMolPropertyKind::SGROUP_DISPLAY;
  else if (input == char4({ 'S', 'C', 'D', 'M' }) || input == char4({ 'S', 'E', 'D', 'M' })) {
    return MdlMolPropertyKind::SGROUP_DATA;
  }
  else if (input == char4({ 'S', 'P', 'L', 'M' })) return MdlMolPropertyKind::SGROUP_HIERARCHY;
  else if (input == char4({ 'S', 'N', 'C', 'M' })) return MdlMolPropertyKind::SGROUP_COMP_NUMBER;
  else if (input == char4({ '$', '3', 'D', 'M' })) return MdlMolPropertyKind::SPATIAL_FEATURE;
  else if (input == char4({ 'P', 'X', 'A', 'M' })) return MdlMolPropertyKind::PHANTOM_ATOM;
  else if (input == char4({ 'S', 'A', 'P', 'M' })) return MdlMolPropertyKind::SGROUP_ATTACH_POINT;
  else if (input == char4({ 'S', 'C', 'L', 'M' })) return MdlMolPropertyKind::SGROUP_CLASS;
  else if (input == char4({ 'R', 'E', 'G', 'M' })) return MdlMolPropertyKind::LARGE_REGNO;
  else if (input == char4({ 'S', 'B', 'T', 'M' })) return MdlMolPropertyKind::SGROUP_BRACKET_STYLE;
  else if (input == char4({ 'S', 'K', 'P', 'S' })) return MdlMolPropertyKind::SKIP;
  else if (input == char4({ ' ', ' ', ' ', ' ' })) return MdlMolPropertyKind::NONE;
  else {
    const std::string str_code = std::to_string(input.w) + "  " + std::to_string(input.x) +
                                 std::to_string(input.y) + std::to_string(input.z);
    rtErr("The MDL MOL property code " + str_code + " is unrecognized.",
          "translateMdlMolPropertyKind");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MdlMolVersion input) {
  switch (input) {
  case MdlMolVersion::V2000:
    return std::string("V2000");
  case MdlMolVersion::V3000:
    return std::string("V3000");
  case MdlMolVersion::UNKNOWN:
    return std::string("UNKNOWN");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MolObjIndexKind input) {
  switch (input) {
  case MolObjIndexKind::ATOM:
    return std::string("ATOM");
  case MolObjIndexKind::BOND:
    return std::string("BOND");
  case MolObjIndexKind::OTHER:
    return std::string("OTHER");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RadicalState input) {
  switch (input) {
  case RadicalState::NONE:
    return std::string("NONE");
  case RadicalState::SINGLET:
    return std::string("SINGLET");
  case RadicalState::DOUBLET:
    return std::string("DOUBLET");
  case RadicalState::TRIPLET:
    return std::string("TRIPLET");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MolObjChirality input) {
  switch (input) {
  case MolObjChirality::ACHIRAL:
    return std::string("ACHIRAL");
  case MolObjChirality::CHIRAL:
    return std::string("CHIRAL");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const HydrogenAssignment input) {
  switch (input) {
  case HydrogenAssignment::VALENCE_SHELL:
    return std::string("VALENCE_SHELL");
  case HydrogenAssignment::DO_NOT_HYDROGENATE:
    return std::string("DO_NOT_HYDROGENATE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MdlMolBondOrder input) {
  switch (input) {
  case MdlMolBondOrder::SINGLE:
    return std::string("SINGLE");
  case MdlMolBondOrder::DOUBLE:
    return std::string("DOUBLE");
  case MdlMolBondOrder::TRIPLE:
    return std::string("TRIPLE");
  case MdlMolBondOrder::AROMATIC:
    return std::string("AROMATIC");
  case MdlMolBondOrder::SINGLE_OR_DOUBLE:
    return std::string("SINGLE_OR_DOUBLE");
  case MdlMolBondOrder::SINGLE_OR_AROMATIC:
    return std::string("SINGLE_OR_AROMATIC");
  case MdlMolBondOrder::DOUBLE_OR_AROMATIC:
    return std::string("DOUBLE_OR_AROMATIC");
  case MdlMolBondOrder::ANY:
    return std::string("ANY");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MdlMolBondStereo input) {
  switch (input) {
  case MdlMolBondStereo::NOT_STEREO:
    return std::string("NOT_STEREO");
  case MdlMolBondStereo::UP:
    return std::string("UP");
  case MdlMolBondStereo::CIS_OR_TRANS:
    return std::string("CIS_OR_TRANS");
  case MdlMolBondStereo::EITHER:
    return std::string("EITHER");
  case MdlMolBondStereo::DOWN:
    return std::string("DOWN");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MolObjAtomStereo input) {
  switch (input) {
  case MolObjAtomStereo::NOT_STEREO:
    return std::string("NOT_STEREO");
  case MolObjAtomStereo::ODD:
    return std::string("ODD");
  case MolObjAtomStereo::EVEN:
    return std::string("EVEN");
  case MolObjAtomStereo::UNMARKED:
    return std::string("UNMARKED");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MolObjRingState input) {
  switch (input) {
  case MolObjRingState::EITHER:
    return std::string("EITHER");
  case MolObjRingState::RING:
    return std::string("RING");
  case MolObjRingState::CHAIN:
    return std::string("CHAIN");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MolObjReactionCenter input) {
  switch (input) {
  case MolObjReactionCenter::NON_CENTER:
    return std::string("NON_CENTER");
  case MolObjReactionCenter::UNMARKED:
    return std::string("UNMARKED");
  case MolObjReactionCenter::CENTER:
    return std::string("CENTER");
  case MolObjReactionCenter::UNREACTIVE:
    return std::string("UNREACTIVE");
  case MolObjReactionCenter::BOND_MADE_OR_BROKEN:
    return std::string("BOND_MADE_OR_BROKEN");
  case MolObjReactionCenter::CENTER_WITH_FORMATION:
    return std::string("CENTER_WITH_FORMATION");
  case MolObjReactionCenter::BOND_ORDER_CHANGE:
    return std::string("BOND_ORDER_CHANGE");
  case MolObjReactionCenter::CENTER_WITH_ORDER_CHANGE:
    return std::string("CENTER_WITH_ORDER_CHANGE");
  case MolObjReactionCenter::BOND_FORMATION_AND_ORDER_CHANGE:
    return std::string("BOND_FORMATION_AND_ORDER_CHANGE");
  case MolObjReactionCenter::CENTER_WITH_FORMATION_AND_ORDER_CHANGE:
    return std::string("CENTER_WITH_FORMATION_AND_ORDER_CHANGE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const StereoRetention input) {
  switch (input) {
  case StereoRetention::NOT_APPLIED:
    return std::string("NOT_APPLIED");
  case StereoRetention::INVERTED:
    return std::string("INVERTED");
  case StereoRetention::RETAINED:
    return std::string("RETAINED");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MolObjPropField input) {
  switch (input) {
  case MolObjPropField::INTEGER:
    return std::string("INTEGER");
  case MolObjPropField::REAL:
    return std::string("REAL");
  case MolObjPropField::CHAR4:
    return std::string("CHAR4");
  case MolObjPropField::STRING:
    return std::string("STRING");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const DataRequestKind input) {
  switch (input) {
  case DataRequestKind::STATE_VARIABLE:
    return std::string("STATE_VARIABLE");
  case DataRequestKind::ATOM_INFLUENCES:
    return std::string("ATOM_INLFUENCES");
  case DataRequestKind::TOPOLOGY_PARAMETER:
    return std::string("TOPOLOGY_PARAMETER");
  case DataRequestKind::STRING:
    return std::string("STRING");
  case DataRequestKind::ALL_KINDS:
    return std::string("ALL_KINDS");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MdlMolDataItemKind input) {
  switch (input) {
  case MdlMolDataItemKind::STATE_VARIABLE:
    return std::string("STATE_VARIABLE");
  case MdlMolDataItemKind::ATOM_INFLUENCES:
    return std::string("ATOM_INFLUENCES");
  case MdlMolDataItemKind::TOPOLOGY_PARAMETER:
    return std::string("TOPOLOGY_PARAMETER");
  case MdlMolDataItemKind::STRING:
    return std::string("STRING");
  case MdlMolDataItemKind::NATIVE:
    return std::string("NATIVE");
  case MdlMolDataItemKind::NONE:
    return std::string("NONE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MdlMolPropertyKind input) {
  switch (input) {
  case MdlMolPropertyKind::ATOM_ALIAS:
    return std::string("ATOM_ALIAS");
  case MdlMolPropertyKind::ATOM_VALUE:
    return std::string("ATOM_VALUE");
  case MdlMolPropertyKind::GROUP_ABBREVIATION:
    return std::string("GROUP_ABBREVIATION");
  case MdlMolPropertyKind::CHARGE:
    return std::string("CHARGE");
  case MdlMolPropertyKind::RADICAL:
    return std::string("RADICAL");
  case MdlMolPropertyKind::ISOTOPE:
    return std::string("ISOTOPE");
  case MdlMolPropertyKind::RING_BOND_COUNT:
    return std::string("RING_BOND_COUNT");
  case MdlMolPropertyKind::SUBSTITUTION_COUNT:
    return std::string("SUBSTITUTION_COUNT");
  case MdlMolPropertyKind::UNSATURATED_COUNT:
    return std::string("UNSATURATED_COUNT");
  case MdlMolPropertyKind::LINK_ATOM:
    return std::string("LINK_ATOM");
  case MdlMolPropertyKind::ATOM_LIST:
    return std::string("ATOM_LIST");
  case MdlMolPropertyKind::ATTACHMENT_POINT:
    return std::string("ATTACHMENT_POINT");
  case MdlMolPropertyKind::ATTACHMENT_ORDER:
    return std::string("ATTACHMENT_ORDER");
  case MdlMolPropertyKind::RGROUP_LABEL_LOCATION:
    return std::string("RGROUP_LABEL_LOCATION");
  case MdlMolPropertyKind::RGROUP_LOGIC:
    return std::string("RGROUP_LOGIC");
  case MdlMolPropertyKind::SGROUP_TYPE:
    return std::string("SGROUP_TYPE");
  case MdlMolPropertyKind::SGROUP_SUBTYPE:
    return std::string("SGROUP_SUBTYPE");
  case MdlMolPropertyKind::SGROUP_LABELS:
    return std::string("SGROUP_LABELS");
  case MdlMolPropertyKind::SGROUP_CONNECTIVITY:
    return std::string("SGROUP_CONNECTIVITY");
  case MdlMolPropertyKind::SGROUP_EXPANSION:
    return std::string("SGROUP_EXPANSION");
  case MdlMolPropertyKind::SGROUP_ATOM_LIST:
    return std::string("SGROUP_ATOM_LIST");
  case MdlMolPropertyKind::SGROUP_BOND_LIST:
    return std::string("SGROUP_BOND_LIST");
  case MdlMolPropertyKind::MG_PARENT_ATOM_LIST:
    return std::string("MG_PARENT_ATOM_LIST");
  case MdlMolPropertyKind::SGROUP_SUBSCRIPT:
    return std::string("SGROUP_SUBSCRIPT");
  case MdlMolPropertyKind::SGROUP_CORRESPONDENCE:
    return std::string("SGROUP_CORRESPONDENCE");
  case MdlMolPropertyKind::SGROUP_DISPLAY_INFO:
    return std::string("SGROUP_DISPLAY_INFO");
  case MdlMolPropertyKind::SGROUP_BOND_VECTOR:
    return std::string("SGROUP_BOND_VECTOR");
  case MdlMolPropertyKind::SGROUP_FIELD:
    return std::string("SGROUP_FIELD");
  case MdlMolPropertyKind::SGROUP_DISPLAY:
    return std::string("SGROUP_DISPLAY");
  case MdlMolPropertyKind::SGROUP_DATA:
    return std::string("SGROUP_DATA");
  case MdlMolPropertyKind::SGROUP_HIERARCHY:
    return std::string("SGROUP_HIERARCHY");
  case MdlMolPropertyKind::SGROUP_COMP_NUMBER:
    return std::string("SGROUP_COMP_NUMBER");
  case MdlMolPropertyKind::SPATIAL_FEATURE:
    return std::string("SPATIAL_FEATURE");
  case MdlMolPropertyKind::PHANTOM_ATOM:
    return std::string("PHANTOM_ATOM");
  case MdlMolPropertyKind::SGROUP_ATTACH_POINT:
    return std::string("SGROUP_ATTACH_POINT");
  case MdlMolPropertyKind::SGROUP_CLASS:
    return std::string("SGROUP_CLASS");
  case MdlMolPropertyKind::LARGE_REGNO:
    return std::string("LARGE_REGNO");
  case MdlMolPropertyKind::SGROUP_BRACKET_STYLE:
    return std::string("SGROUP_BRACKET_STYLE");
  case MdlMolPropertyKind::SKIP:
    return std::string("SKIP");
  case MdlMolPropertyKind::NONE:
    return std::string("NONE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MolObjFeatureKind input) {
  switch (input) {
  case MolObjFeatureKind::PT_FD_2:
    return std::string("PT_FD_2");
  case MolObjFeatureKind::PT_FLEX_2:
    return std::string("PT_FLEX_2");
  case MolObjFeatureKind::PT_PN_NORMAL:
    return std::string("PT_N_NORMAL");
  case MolObjFeatureKind::LN_BEST_FIT:
    return std::string("LN_BEST_FIT");
  case MolObjFeatureKind::PN_BEST_FIT:
    return std::string("PN_BEST_FIT");
  case MolObjFeatureKind::PN_PT_LN:
    return std::string("PN_PT_LN");
  case MolObjFeatureKind::CENTROID:
    return std::string("CENTROID");
  case MolObjFeatureKind::LN_PN_NORMAL:
    return std::string("LN_PN_NORMAL");
  case MolObjFeatureKind::DISTANCE_PT_PT:
    return std::string("DISTANCE_PT_PT");
  case MolObjFeatureKind::DISTANCE_PT_LN:
    return std::string("DISTANCE_PT_LN");
  case MolObjFeatureKind::DISTANCE_PT_PN:
    return std::string("DISTANCE_PT_PN");
  case MolObjFeatureKind::ANGLE_PT_PT_PT:
    return std::string("ANGLE_PT_PT_PT");
  case MolObjFeatureKind::ANGLE_LN_LN:
    return std::string("ANGLE_LN_LN");
  case MolObjFeatureKind::ANGLE_PN_PN:
    return std::string("ANGLE_PN_PN");
  case MolObjFeatureKind::SPHERE_PT:
    return std::string("SPHERE_PT");
  case MolObjFeatureKind::FIXED_ATOMS:
    return std::string("FIXED_ATOMS");
  case MolObjFeatureKind::ATOM_CONSTRAINT:
    return std::string("ATOM_CONSTRAINT");
  case MolObjFeatureKind::PAIR_CONSTRAINT:
    return std::string("PAIR_CONSTRAINT");
  }
  __builtin_unreachable();
}

} // namespace structure
} // namespace stormm
