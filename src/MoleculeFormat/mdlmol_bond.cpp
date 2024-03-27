#include "copyright.h"
#include "Parsing/ascii_numbers.h"
#include "Parsing/parse.h"
#include "FileManagement/file_listing.h"
#include "mdlmol_bond.h"

namespace stormm {
namespace structure {

using diskutil::getBaseName;
using parse::readIntegerValue;
using parse::verifyContents;
  
//-------------------------------------------------------------------------------------------------
MdlMolBond::MdlMolBond() :
    i_atom{-1}, j_atom{-1}, order{MdlMolBondOrder::SINGLE}, stereo{MdlMolBondStereo::NOT_STEREO},
    ring_state{MolObjRingState::EITHER}, reactivity{MolObjReactionCenter::NON_CENTER}
{}

//-------------------------------------------------------------------------------------------------
MdlMolBond::MdlMolBond(const int i_atom_in, const int j_atom_in) :
    i_atom{i_atom_in}, j_atom{j_atom_in}, order{MdlMolBondOrder::SINGLE},
    stereo{MdlMolBondStereo::NOT_STEREO}, ring_state{MolObjRingState::EITHER},
    reactivity{MolObjReactionCenter::NON_CENTER}
{}

//-------------------------------------------------------------------------------------------------
MdlMolBond::MdlMolBond(const int i_atom_in, const int j_atom_in, const MdlMolBondOrder order_in,
                       const MdlMolBondStereo stereo_in, const MolObjRingState ring_state_in,
                       const MolObjReactionCenter reactivity_in) :
    i_atom{i_atom_in}, j_atom{j_atom_in}, order{order_in}, stereo{stereo_in},
    ring_state{ring_state_in}, reactivity{reactivity_in}
{}

//-------------------------------------------------------------------------------------------------
MdlMolBond::MdlMolBond(const TextFile &tf, const int line_number, const std::string &title) :
    MdlMolBond()
{
  const char* bond_line_ptr = tf.getLinePointer(line_number);
  if (tf.getLineLength(line_number) < 6) {
    rtErr("Line " + std::to_string(line_number) + " of " + getBaseName(tf.getFileName()) +
          " cannot contain MDL MOL bond information due to its length being only " +
          std::to_string(tf.getLineLength(line_number)));
  }

  // Subtract 1 from the atom indexing given in the MDL MOL format, which starts at 1.
  i_atom = readIntegerValue(bond_line_ptr, 0, 3) - 1;
  j_atom = readIntegerValue(bond_line_ptr, 3, 3) - 1;
  order = (verifyContents(tf, line_number, 6, 3)) ?
          interpretBondOrder(readIntegerValue(bond_line_ptr, 6, 3), title) :
          default_mdl_bond_order;
  stereo = (verifyContents(tf, line_number, 9, 3)) ?
           interpretBondStereochemistry(readIntegerValue(bond_line_ptr, 9, 3), title) :
           default_mdl_bond_stereochemistry;
  ring_state = (verifyContents(tf, line_number, 12, 3)) ?
               interpretRingState(readIntegerValue(bond_line_ptr, 12, 3), title) :
               default_mdl_ring_status;
  reactivity = (verifyContents(tf, line_number, 15, 3)) ?
               interpretBondReactivePotential(readIntegerValue(bond_line_ptr, 15, 3), title) :
               default_mdl_bond_reactivity;
}

//-------------------------------------------------------------------------------------------------
MdlMolBondOrder MdlMolBond::interpretBondOrder(const int code_in, const std::string &title) {
  switch (code_in) {
  case 1:
    return MdlMolBondOrder::SINGLE;
  case 2:
    return MdlMolBondOrder::DOUBLE;
  case 3:
    return MdlMolBondOrder::TRIPLE;
  case 4:
    return MdlMolBondOrder::AROMATIC;
  case 5:
    return MdlMolBondOrder::SINGLE_OR_DOUBLE;
  case 6:
    return MdlMolBondOrder::SINGLE_OR_AROMATIC;
  case 7:
    return MdlMolBondOrder::DOUBLE_OR_AROMATIC;
  case 8:
    return MdlMolBondOrder::ANY;
  default:
    rtErr("An invalid bond order code " + std::to_string(code_in) + " was specified in an MDL MOL "
          "entry titled \"" + title + "\".", "MdlMolObj", "interpretBondOrder");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
MdlMolBondStereo MdlMolBond::interpretBondStereochemistry(const int code_in,
                                                          const std::string &title) {
  switch (code_in) {
  case 0:
    return MdlMolBondStereo::NOT_STEREO;
  case 1:
    return MdlMolBondStereo::UP;
  case 3:
    return MdlMolBondStereo::CIS_OR_TRANS;
  case 4:
    return MdlMolBondStereo::EITHER;
  case 6:
    return MdlMolBondStereo::DOWN;
  default:
    rtErr("An invalid bond stereochemistry code " + std::to_string(code_in) + " was specified in "
          "an MDL MOL entry titled \"" + title + "\".", "MdlMolObj",
          "interpretBondStereoChemistry");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
MolObjRingState MdlMolBond::interpretRingState(const int code_in, const std::string &title) {
  switch (code_in) {
  case 0:
    return MolObjRingState::EITHER;
  case 1:
    return MolObjRingState::RING;
  case 2:
    return MolObjRingState::CHAIN;
  default:
    rtErr("An invalid bond ring status code " + std::to_string(code_in) + " was specified in an "
          "MDL MOL entry titled \"" + title + "\".", "MdlMolObj", "interpretRingState");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
MolObjReactionCenter MdlMolBond::interpretBondReactivePotential(const int code_in,
                                                                const std::string &title) {
  switch (code_in) {
  case -1:
    return MolObjReactionCenter::NON_CENTER;
  case 0:
    return MolObjReactionCenter::UNMARKED;
  case 1:
    return MolObjReactionCenter::CENTER;
  case 2:
    return MolObjReactionCenter::UNREACTIVE;
  case 4:
    return MolObjReactionCenter::BOND_MADE_OR_BROKEN;
  case 5:
    return MolObjReactionCenter::CENTER_WITH_FORMATION;
  case 8:
    return MolObjReactionCenter::BOND_ORDER_CHANGE;
  case 9:
    return MolObjReactionCenter::CENTER_WITH_ORDER_CHANGE;
  case 12:
    return MolObjReactionCenter::BOND_FORMATION_AND_ORDER_CHANGE;
  case 13:
    return MolObjReactionCenter::CENTER_WITH_FORMATION_AND_ORDER_CHANGE;
  default:
    rtErr("An invalid bond reactive potential code " + std::to_string(code_in) + " was specified "
          "in an MDL MOL entry titled \"" + title + "\".", "MdlMolObj",
          "interpretBondReactivePotential");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int MdlMolBond::getFirstAtom() const {
  return i_atom;
}

//-------------------------------------------------------------------------------------------------
int MdlMolBond::getSecondAtom() const {
  return j_atom;
}

//-------------------------------------------------------------------------------------------------
MdlMolBondOrder MdlMolBond::getOrder() const {
  return order;
}

//-------------------------------------------------------------------------------------------------
MdlMolBondStereo MdlMolBond::getStereochemistry() const {
  return stereo;
}

//-------------------------------------------------------------------------------------------------
MolObjRingState MdlMolBond::getRingStatus() const {
  return ring_state;
}

//-------------------------------------------------------------------------------------------------
MolObjReactionCenter MdlMolBond::getReactivity() const {
  return reactivity;
}

//-------------------------------------------------------------------------------------------------
void MdlMolBond::setFirstAtom(const int index_in) {
  i_atom = index_in;
}

//-------------------------------------------------------------------------------------------------
void MdlMolBond::setSecondAtom(const int index_in) {
  j_atom = index_in;
}

//-------------------------------------------------------------------------------------------------
void MdlMolBond::setOrder(const MdlMolBondOrder order_in) {
  order = order_in;
}

//-------------------------------------------------------------------------------------------------
void MdlMolBond::setStereochemistry(const MdlMolBondStereo stereo_in) {
  stereo = stereo_in;
}

//-------------------------------------------------------------------------------------------------
void MdlMolBond::setRingStatus(const MolObjRingState status_in) {
  ring_state = status_in;
}

//-------------------------------------------------------------------------------------------------
void MdlMolBond::setReactivity(const MolObjReactionCenter potential_in) {
  reactivity = potential_in;
}

} // namespace structure
} // namespace stormm
