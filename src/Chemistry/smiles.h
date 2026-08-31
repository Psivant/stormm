// -*-c++-*-
#ifndef STORMM_SMILES_H
#define STORMM_SMILES_H

#include <string>
#include <vector>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "chemistry_enumerators.h"

namespace stormm {
namespace chemistry {

/// \brief Encapsulate the methods for determining a single atom, or a heavy atom with one or more
///        attached hydrogen atoms.
class SmilesAtom {
public:

  /// \brief The constructor takes a string of characters and will not subdivide it further.
  SmilesAtom(const std::string &content_in, const std::vector<int> &major_bond_orders_in);

  /// \brief Get the atomic number.
  int getAtomicNumber() const;

  /// \brief Get the number of attached hydrogens.
  int getAttachedHydrogenCount() const;

  /// \brief Get the formal charge of the (major) atom.
  int getFormalCharge() const;

  /// \brief Get the isotopic number of the (major) atom.
  int getIsotopicNumber() const;
  
  /// \brief Get the chirality of the atom.
  ChiralOrientation getChirality() const;

  /// \brief Indicate whether the atom is aromatic.
  bool isAromatic() const;

  /// \brief Adjust the major bound count
  ///
  /// \param major_bond_orders_in  The new list of major bond orders
  void setMajorBondOrders(const std::vector<int> &major_bond_orders_in);

  /// \brief Add a bond to another major atom
  ///
  /// \param order  The order of the new bond to add
  void addMajorBond(int order);
  
private:
  int atomic_number;             ///< The atomic number of the atom, indicating its element
  int attached_hydrogens;        ///< The number of attached hydrogens
  bool inferred_hydrogen_count;  ///< Indication of whether the number of attached hydrogens has
                                 ///<   been inferred, and thus should be adjusted in response to
                                 ///<   any updates in the bond count
  int formal_charge;             ///< Formal charge of the (major) atom
  int isotopic_number;           ///< The isotopic number of the (major) atom
  bool is_aromatic;              ///< Indicate whether the atom is part of an aromatic group
  ChiralOrientation chirality;   ///< Chirality of the atom center (this can be more complex than
                                 ///<   mere R- and S- organic carbon centers)
  int chiral_order;              ///< The order of the chirality, indicating the number of branches
                                 ///<   involved in the center.  Four indicates tetrahedral,
                                 ///<   negative four square planar, five trigonal bipyramidal,
                                 ///<   and six octahedral (square bipyramidal).
  int chiral_id;                 ///< Identifying number associated with the chirality type
  std::string content;           ///< The textual basis for the atom, recorded for error tracing
                                 ///<   purposes

  /// A list of the orders for bonds connecting the atom to its surroundings
  std::vector<int> major_bond_orders;

  /// \brief Set the number of attached hydrogens, in essence the number of valence electrons for
  ///        the atom.  This is done automatically during construction and updates to the 
  ///        the class object that might have bearing on
  ///        the number of valence electrons.
  void updateAttachedHydrogenCount();

  /// \brief Trim brackets from the content so as to make a clean presentation to put errors or
  ///        other messages.
  std::string bracketlessContent() const;
};
  
/// \brief Collect the methods and annotations needed to parse a SMILES string.
class SmilesString {
public:

  /// \brief The constructor takes a Standard Template Library string, or something that could be
  ///        used to create one.
  /// \{
  SmilesString(const std::string &basis_in);
  SmilesString(const char* basis_in);

  /// \brief With Standard Template Library classes to handle all components, the default copy and
  ///        move constructors, as well as copy and move assignment operators, are valid.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right and side of the assignment statement
  /// \{
  SmilesString(const SmilesString &original) = default;
  SmilesString(SmilesString &&original) = default;
  SmilesString& operator=(const SmilesString &original) = default;
  SmilesString& operator=(SmilesString &&original) = default;
  /// \}

  /// \brief Get the total number of atoms in the system.  The original SMILES string may only
  ///        present heavy atom content, but this count will include all hydrogens.
  int getAtomCount() const;

  /// \brief Get the total number of involving one atom or the total number of bonds in the system,
  ///        both between atoms explicitly stated in the original basis as well as bonds between
  ///        those atoms and any implied hydrogens.
  ///
  /// \param index  Index of an atom of interest
  /// \{
  int getBondCount() const;
  int getBondCount(int index) const;
  /// \}

  /// \brief Get the Z number of one of the atoms in the system.
  ///
  /// \param index  Index of the atom of interest
  int getAtomicNumber(int index) const;

  /// \brief Get the Z numbers of all atoms in the system.
  const std::vector<int>& getAtomicNumbers() const;

  /// \brief Get the isotopic number of one of the atoms in the system.
  ///
  /// \param index  Index of the atom of interest
  int getIsotopicNumber(int index) const;
  
  /// \brief Get the order of a bond between two atoms.  The request will be checked for validity
  ///        (that there is a bond between the two atoms) and return the value stated in the
  ///        original SMILES string.  A return value of -1 indicates that the bond is part of an
  ///        aromatic ring system.
  ///
  /// \param index_a  The index of the first atom in the bond
  /// \param index_b  The index of the second atom in the bond
  int getBondOrder(int index_a, int index_b) const;

  /// \brief Get the number of atoms to which one atom of interest is bonded.  This count will
  ///        include added hydrogens.
  ///
  /// \param index  Index of the atom of interest
  int getBondedPartnerCount(int index) const;

  /// \brief For a specified atom, get the indices of bonded atom partners within the system.
  ///
  /// \param index  Index of the atom of interest
  std::vector<int> getBondedPartnerIndices(const int index) const;

  /// \brief Present the original SMILES string.
  const std::string& getBasis() const;
  
private:
  std::string basis;                 ///< The content of the original SMILES string
  int atom_count;                    ///< The number of atoms, both in the original string and
                                     ///<   perhaps inferred (unmentioned hydrogens)
  int bond_count;                    ///< The number of bonds between all atoms, both mentioned in
                                     ///<   the SMILES string and inferred between those atoms and
                                     ///<   added hydrogens
  int net_charge;                    ///< The net formal charge of the molecule, a sum of all
                                     ///<   atoms' formal charges
  std::vector<int> atomic_numbers;   ///< Element identifiers for all atoms in order
  std::vector<int> isotopic_numbers; ///< Isotope numbers for each atom
  std::vector<int> formal_charges;   ///< Formal charges of all atoms, based on information given
                                     ///<   in the SMILES string
  std::vector<int2> bonds;           ///< Bonds between all connected atoms
  std::vector<int> bond_orders;      ///< Orders of all bonds presented in the SMILES string
                                     ///<   itself.  These  may or may not match the bond orders
                                     ///<   decided by the chemical perception in the
                                     ///<   ChemicalFeatures class.
  std::vector<bool> aromatic_bonds;  ///< Flags for every bond in the system to indicate its
                                     ///<   aromaticity
  std::vector<int> bonded_partners;  ///< The other atoms to which each atom makes a bond.  This
                                     ///<   list is reflexive, i.e. each bond is mentioned twice.
  std::vector<int> partner_bounds;   ///< The bound list for bonded_partners, above.  The partners
                                     ///<   for the kth atom are given in indices partner_bounds[k]
                                     ///<   to partner_bounds[k + 1] of bonded_partners.

  /// Indices into the array member variable bonds based on the ordering of bond partners found in
  /// bonded_partners.  This is a reverse lookup table.  To find the order of the bond between
  /// atoms k and m, search bonded_partners in elements partner_bounds[k] to partner_bounds[k + 1],
  /// find the index of bonded_partners at which the partner is m (let this be m_at_k), then get
  /// the bond index from bond_indices[m_at_k] and look up bond_orders[bond_indices[m_at_k]].
  std::vector<int> partner_bond_indices;
  
  /// \brief Validate an atom index for the system.
  ///
  /// \param index  The atom of interest
  void validateAtomIndex(int index) const;

  /// \brief Parse an atom from within a pair of square brackets.  Return the enclosed string
  ///        content.
  ///
  /// \param content    Text string providing the atom's text as well as larger context, such as
  ///                   the entire branch or even the entire molecule
  /// \param parse_pos  The position
  std::string parseBracketedAtom(const std::string &content, int *parse_pos);

  /// \brief Parse a branch from within a pair of parentheses.
  ///
  ///
  void parseBranch(const std::string &content, std::vector<std::string> *atom_content,
                   std::vector<int3> *bond_content, std::vector<int3> *ring_closures,
                   int branch_initiating_atom, int anchor_bond_order);

  /// \brief Detect whether a bond between major atoms is part of an aromatic system, and flag any
  ///        other bonds in the aromatic system as well.
  ///
  /// \param iatom                      The first atom in the bond
  /// \param jatom                      The second atom in the bond
  /// \param ij_setting                 Index of the partner (out of the entire, concatenated list
  ///                                   of bond partners) where the second atom is connected to the
  ///                                   first
  /// \param major_atoms                A list of major atoms in the chemical system
  /// \param major_bonds                A list of major bonds in the chemical system.  This is
  ///                                   needed, in fact, because some "major bonds" may be of order
  ///                                   zero to indicate that there is no bond between the
  ///                                   consecutive atoms.  The atoms at either end of a bond are
  ///                                   given in the "x" and "y" members of each tuple, while the
  ///                                   bond order is given in the "z" member.
  /// \param major_bond_indices         
  /// \param major_bond_partners
  /// \param major_bond_partner_bounds
  /// \param major_bond_aromaticity     A list of flags for each bond, TRUE if the bond is part of
  ///                                   an aromatic system and FALSE otherwise.  Modified and
  ///                                   returned.
  void detectAromaticBonds(int iatom, int jatom, int ij_setting,
                           const std::vector<SmilesAtom> &major_atoms,
                           const std::vector<int3> &major_bonds,
                           const std::vector<int> &major_bond_indices,
                           const std::vector<int> &major_bond_partners,
                           const std::vector<int> &major_bond_partner_bounds,
                           std::vector<bool> *major_bond_aromaticity);
};

} // namespace chemistry
} // namespace stormm

#endif
