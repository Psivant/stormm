// -*-c++-*-
#ifndef STORMM_MOLOBJ_BOND_H
#define STORMM_MOLOBJ_BOND_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Parsing/textfile.h"
#include "molecule_format_enumerators.h"

namespace stormm {
namespace structure {

using parse::TextFile;
  
/// \brief Default settings for the MDL MOL bond initializations
/// \{
constexpr MdlMolBondOrder default_mdl_bond_order = MdlMolBondOrder::SINGLE;
constexpr MdlMolBondStereo default_mdl_bond_stereochemistry = MdlMolBondStereo::NOT_STEREO;
constexpr MolObjRingState default_mdl_ring_status = MolObjRingState::EITHER;
constexpr MolObjReactionCenter default_mdl_bond_reactivity = MolObjReactionCenter::NON_CENTER;
/// \}
  
/// \brief A bond, as presented in the MDL molecule file format.  This unguarded struct will be
///        returned to the developer from a private array inside of the MdlMolObj object, so
///        further protection would be a hindrance.
class MdlMolBond {
public:

  /// \brief The constructor can take all member variables, or just the atoms so that more
  ///        information can be filled in later.
  ///
  /// Overloaded:
  ///   - Construct a blank object with -1 atom indices to indicate its invalid nature
  ///   - Construct an object with just the atom indices
  ///   - Construct a complete object with all details pre-loaded
  ///   - Construct the object based on text input of trusted bounds
  ///
  /// \param tf           Text of the original .sdf or .mol file, read into RAM
  /// \param line_number  Number of the line on which to read the data
  /// \param title        The title of the structure, if known, for error tracing purposes
  /// \{
  MdlMolBond();
  MdlMolBond(int i_atom_in, int j_atom_in);
  MdlMolBond(int i_atom_in, int j_atom_in, MdlMolBondOrder order_in, MdlMolBondStereo stereo_in,
             MolObjRingState ring_state_in, MolObjReactionCenter reactivity_in);
  MdlMolBond(const TextFile &tf, int line_number, const std::string &title = std::string(""));
  /// \}

  /// \brief The default copy and move constructors as well as assignment operators are adequate.
  /// \{
  MdlMolBond(const MdlMolBond &original) = default;
  MdlMolBond(MdlMolBond &&original) = default;
  MdlMolBond& operator=(const MdlMolBond &other) = default;
  MdlMolBond& operator=(MdlMolBond &&other) = default;
  /// \}

  /// \brief Get the first atom in the bond.  Having separate functions for each atom is a more
  ///        intuitive way to offer the getter functions, whereas in a BoundedRestraint object the
  ///        atom index getter takes an argument for the first, second, third, or fourth atom.
  int getFirstAtom() const;

  /// \brief Get the second atom in the bond.
  int getSecondAtom() const;

  /// \brief Get the order of the bond.
  MdlMolBondOrder getOrder() const;

  /// \brief Get the order of the bond.
  MdlMolBondStereo getStereochemistry() const;

  /// \brief Get the ring status--is the bond known to be part of a ring?
  MolObjRingState getRingStatus() const;

  /// \brief Get the reactive potential of the bond.
  MolObjReactionCenter getReactivity() const;

  /// \brief Set the index of the first atom in the bond.
  void setFirstAtom(int index_in);

  /// \brief Set the index of the second atom in the bond.
  void setSecondAtom(int index_in);

  /// \brief Set the order of the bond, perhaps after computations with an associated
  ///        ChemicalFeatures object.
  void setOrder(MdlMolBondOrder order_in);

  /// \brief Set the stereochemical details of the bond.
  void setStereochemistry(MdlMolBondStereo stereo_in);

  /// \brief Mark the status of the bond with respect to any ring features.
  void setRingStatus(MolObjRingState status_in);

  /// \brief Mark the reactive potential of the bond.
  void setReactivity(MolObjReactionCenter potential_in);

private:
  int i_atom;                       ///< The first atom in the bond
  int j_atom;                       ///< The second atom in the bond
  MdlMolBondOrder order;            ///< The bond order (single, double, aromatic, etc.)
  MdlMolBondStereo stereo;          ///< Indicator of the bond stereochemistry
  MolObjRingState ring_state;       ///< Indicator of whether the atom is part of a ring
  MolObjReactionCenter reactivity;  ///< Indicator of a bond as a center of reactivity

  /// \brief Interpret a code for the order of a bond in the structure.
  ///
  /// \param code_in  A numeric code to be translated into the bond order
  /// \param title    Title of the parent structure, for error tracing purposes
  MdlMolBondOrder interpretBondOrder(int code_in, const std::string &title);

  /// \brief Interpret a code for the stereochemistry of a bond in the structure.
  ///
  /// \param code_in  A numeric code to be translated into the bond stereochemistry
  /// \param title    Title of the parent structure, for error tracing purposes
  MdlMolBondStereo interpretBondStereochemistry(int code_in, const std::string &title);

  /// \brief Interpret a code for the ring status of a bond in the structure.
  ///
  /// \param code_in  A numeric code to be translated into the ring status
  /// \param title    Title of the parent structure, for error tracing purposes
  MolObjRingState interpretRingState(int code_in, const std::string &title);

  /// \brief Interpret a code for the reactive potential of a bond in the structure.
  ///
  /// \param code_in  A numeric code to be translated into the reactive potential
  /// \param title    Title of the parent structure, for error tracing purposes
  MolObjReactionCenter interpretBondReactivePotential(int code_in, const std::string &title);
};

/// \brief Overload the + operator to concatenate vectors of MDL and SDF bonds.
std::vector<MdlMolBond> operator+(const std::vector<MdlMolBond> &lhs,
                                  const std::vector<MdlMolBond> &rhs);

} // namespace structure
} // namespace stormm

#endif
