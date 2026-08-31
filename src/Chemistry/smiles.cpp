#include "copyright.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "periodic_table.h"
#include "smiles.h"

namespace stormm {
namespace chemistry {

using stmath::maxValue;
using stmath::prefixSumInPlace;
using stmath::PrefixSumType;
using stmath::sum;
  
//-------------------------------------------------------------------------------------------------
SmilesAtom::SmilesAtom(const std::string &content_in,
                       const std::vector<int> &major_bond_orders_in) :
    atomic_number{0}, attached_hydrogens{0}, formal_charge{0}, isotopic_number{0},
    is_aromatic{false}, chirality{ChiralOrientation::NONE}, chiral_order{0}, chiral_id{0},
    content{content_in}, major_bond_orders{major_bond_orders_in}
{
  if (content.size() == 0) {
    rtErr("An atom cannot be created based on no content.", "SmilesAtom");
  }
  const int zero_int = '0';
  if ((content[0] == '[' && content.back() != ']') ||
      (content[0] != '[' && content.back() == ']')) {
    rtErr("Incomplete bracketing detected around atom \"" + content + "\".", "SmilesAtom");
  }
  const int max_pos = content.size() - (content.back() == ']');
  const int min_pos = (content[0] == '[');
  int parse_pos = min_pos;
  while (parse_pos < max_pos && content[parse_pos] >= '0' && content[parse_pos] <= '9') {
    parse_pos++;
  }
  if (parse_pos > 0) {
    int rpos = parse_pos - 1;
    int mult = 1;
    while (rpos >= min_pos) {
      isotopic_number += (static_cast<int>(content[rpos]) - zero_int) * mult;
      mult *= 10;
      rpos--;
    }
  }
  
  // Determine the one or two letters of the major atom
  char2 major_atom;
  int major_atom_width;
  if ((content[parse_pos] < 'a' || content[parse_pos] > 'z') &&
      (content[parse_pos] < 'A' || content[parse_pos] > 'Z')) {
    rtErr("No element was found in atom specification [" + bracketlessContent() + "].",
          "SmilesAtom");
  }
  if (parse_pos + 1 < max_pos &&
      ((content[parse_pos + 1] >= 'a' && content[parse_pos + 1] <= 'z') ||
       (content[parse_pos + 1] >= 'A' && content[parse_pos + 1] <= 'Z'))) {
    if (content[parse_pos + 1] == 'H' || content[parse_pos + 1] == 'h') {
      if (content[parse_pos] == 'R' || content[parse_pos] == 'r' ||
          content[parse_pos] == 'T' || content[parse_pos] == 't') {
        major_atom = { content[parse_pos], content[parse_pos + 1] };
        major_atom_width = 2;
      }
      else {
        major_atom = { content[parse_pos], ' ' };
        major_atom_width = 1;
      }
    }
    else {
      major_atom = { content[parse_pos], content[parse_pos + 1] };
      major_atom_width = 2;
    }
  }
  else {
    major_atom = { content[parse_pos], ' ' };
    major_atom_width = 1;
  }
  parse_pos += major_atom_width;

  // Look for chiral markers
  bool decide_chirality = false;
  if (parse_pos < max_pos && content[parse_pos] == '@') {
    parse_pos++;

    // Check for shorthand (multiple consecutive '@' symbols)
    int n_at_symbols = 1;
    while (parse_pos < max_pos && content[parse_pos] == '@') {
      n_at_symbols++;
      parse_pos++;
    }
    if (n_at_symbols == 1) {
      
      // There may be more elaborate chiralities defined in the file.  That can be determined with
      // particular letters following the '@' symbols.
      if (parse_pos < max_pos - 3) {
        if (content[parse_pos    ] == 'T' && content[parse_pos + 1] == 'H' &&
            (content[parse_pos + 2] == '1' || content[parse_pos + 2] == '2')) {
          chiral_order = 4;
          chiral_id = static_cast<int>(content[parse_pos + 2]) - zero_int;
          if (content[parse_pos + 2] == '1') {
            chirality = ChiralOrientation::RECTUS;
          }
          else {
            chirality = ChiralOrientation::SINISTER;
          }
          parse_pos += 3;
        }
        else if ((content[parse_pos] == 'T' && content[parse_pos + 1] == 'B') ||
                 (content[parse_pos] == 'S' && content[parse_pos + 1] == 'P') ||
                 (content[parse_pos] == 'O' && content[parse_pos + 1] == 'H')) {
          int chnum_pos = parse_pos + 2;
          while (chnum_pos < max_pos && content[chnum_pos] >= '0' && content[chnum_pos] <= '9') {
            chnum_pos++;
          }
          chnum_pos -= parse_pos + 2;
          if (chnum_pos == 0) {
            rtErr("No numerical input was detected after the long-form chirality specification in "
                  "atom [" + bracketlessContent() + "].", "SmilesAtom");
          }
          else if (chnum_pos == 1) {
            chiral_id = static_cast<int>(content[parse_pos + 2]) - zero_int;
          }
          else if (chnum_pos == 2) {
            chiral_id = (10 * (static_cast<int>(content[parse_pos + 3]) - zero_int)) +
                        static_cast<int>(content[parse_pos + 2]) - zero_int;
          }
          else {
            rtErr("Invalid numerical input was detected after the long-form chirality "
                  "specification in atom [" + bracketlessContent() + "].", "SmilesAtom");
          }
          parse_pos += 2 + chnum_pos;
        }
      }
      else {

        // There is one '@' symbol and it is giving a shorthand chiral designation.  Use the bond
        // count to determine the chiral type.
        chiral_id = 1;
        decide_chirality = true;
      }
    }
    else {
      chiral_id = n_at_symbols;
      decide_chirality = true;

      if (major_bond_orders.size() == 4) {
        if (chiral_id == 2) {
          chirality = ChiralOrientation::SINISTER;
        }
        else {
          rtErr("The only option for a third or higher chiral index (in atom [" +
                bracketlessContent() + "]) is to specify @SP... (explicit chiral designation).  "
                "Otherwise, tetrahedral chirality can only have two possible orientations (@ or "
                "@@ in the shorthand).", "SmilesAtom");
        }
      }
    }
  }
  
  // The next part of the atom specification may indicate attached hydrogens.
  if (parse_pos < max_pos && (content[parse_pos] == 'H' || content[parse_pos + 1] == 'h')) {
    parse_pos++;
    if (parse_pos < max_pos && content[parse_pos] >= '1' && content[parse_pos] <= '9') {
      attached_hydrogens = static_cast<int>(content[parse_pos]) - zero_int;
      parse_pos++;
    }
    else {
      attached_hydrogens = 1;
    }
  }
  
  // Check for a formal charge
  if (parse_pos < max_pos && (content[parse_pos] == '+' || content[parse_pos] == '-')) {
    const int q_mult = (content[parse_pos] == '+') ? 1 : -1;
    parse_pos++;
    if (parse_pos < max_pos) {
      if (content[parse_pos] >= '0' && content[parse_pos] <= '9') {
        formal_charge = q_mult * (static_cast<int>(content[parse_pos]) - zero_int);
        parse_pos++;
      }
      else {
        int q_count = 1;
        while (parse_pos < max_pos &&
               (content[parse_pos] == '+' || content[parse_pos] == '-')) {
          if ((q_mult ==  1 && content[parse_pos] == '-') ||
              (q_mult == -1 && content[parse_pos] == '+')) {
            rtErr("Error in formal charge assignment for atom [" + bracketlessContent() + "].",
                  "SmilesAtom");
          }
          q_count++;
          parse_pos++;
        }
        formal_charge = q_mult * q_count;
      }
    }
    else {
      formal_charge = q_mult;
    }
  }

  // Assign an atomic number based on the atom found.  Search the most common cases before
  // performing a complete loop over all known elements.  Also list out aromatic cases, as
  // lowercase letters denote aromaticity in SMILES strings but the list of elements uses the
  // typical IUPAC convention to give the first letter as a capital and any subsequent letters of
  // the element symbol in lowercase.
  if (major_atom.x == 'C' && major_atom.y == ' ') {
    atomic_number = 6;
  }
  else if (major_atom.x == 'N' && major_atom.y == ' ') {
    atomic_number = 7;
  }
  else if (major_atom.x == 'O' && major_atom.y == ' ') {
    atomic_number = 8;
  }
  else if (major_atom.x == 'H' && major_atom.y == ' ') {
    atomic_number = 1;
  }
  else if (major_atom.x == 'S' && major_atom.y == ' ') {
    atomic_number = 16;
  }
  else if (major_atom.x == 'P' && major_atom.y == ' ') {
    atomic_number = 15;
  }
  else if (major_atom.x == 'F' && major_atom.y == ' ') {
    atomic_number = 9;
  }
  else if (major_atom.x == 'C' && major_atom.y == 'l') {
    atomic_number = 17;
  }
  else if (major_atom.x == 'B' && major_atom.y == 'r') {
    atomic_number = 35;
  }
  else if (major_atom.x == 'I' && major_atom.y == ' ') {
    atomic_number = 53;
  }
  else if (major_atom.x == 'c' && major_atom.y == ' ') {
    atomic_number = 6;
    is_aromatic = true;
  }
  else if (major_atom.x == 'n' && major_atom.y == ' ') {
    atomic_number = 7;
    is_aromatic = true;
  }
  else if (major_atom.x == 'o' && major_atom.y == ' ') {
    atomic_number = 8;
    is_aromatic = true;
  }
  else if (major_atom.x == 's' && major_atom.y == ' ') {
    atomic_number = 16;
    is_aromatic = true;
  }
  else if (major_atom.x == 'b' && major_atom.y == ' ') {
    atomic_number = 5;
    is_aromatic = true;
  }
  else if (major_atom.x == 's' && major_atom.y == 'i') {
    atomic_number = 14;
    is_aromatic = true;
  }
  else if (major_atom.x == 'p' && major_atom.y == ' ') {
    atomic_number = 15;
    is_aromatic = true;
  }
  else if (major_atom.x == 'a' && major_atom.y == 's') {
    atomic_number = 33;
    is_aromatic = true;
  }
  else if (major_atom.x == 's' && major_atom.y == 'e') {
    atomic_number = 34;
    is_aromatic = true;
  }
  else {
    for (int i = 0; i < element_maximum_count + 1; i++) {
      if (major_atom == elemental_symbols[i]) {
        atomic_number = i;
        break;
      }
    }
  }

  // Set the isotope, if it has not been specified.
  if (isotopic_number == 0) {
    isotopic_number = round(elemental_masses[atomic_number]);
  }

  // Adjust the number of attached hydrogens
  if (attached_hydrogens == 0) {
    inferred_hydrogen_count = true;
    updateAttachedHydrogenCount();
  }
  else {
    inferred_hydrogen_count = false;
  }

  // Check that there are sufficient numbers of unique atoms bonded to support the chiral
  // designation.  While this is not a rigorous check (multiple groups based on other major
  // atoms bonded to this one may be identical), it will serve to filter out many bad inputs.
  if (decide_chirality) {
    const int total_bonds = attached_hydrogens + major_bond_orders.size();
    const int unique_bonds = major_bond_orders.size() + (attached_hydrogens > 0);
    if (unique_bonds < 4) {
      rtErr("There are not enough unique bonds to atom [" + bracketlessContent() + "] to create "
            "a chiral center.", "SmilesAtom");
    }
    if (chiral_id == 1) {
      if (total_bonds == 4) {
        chirality = ChiralOrientation::RECTUS;
      }
      else if (total_bonds >= 5) {
        chirality = ChiralOrientation::NONE;
      }
    }
    else if (chiral_id == 2) {
      if (total_bonds == 4) {
        chirality = ChiralOrientation::SINISTER;
      }
      else {
        chirality = ChiralOrientation::NONE;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
int SmilesAtom::getAtomicNumber() const {
  return atomic_number;
}
  
//-------------------------------------------------------------------------------------------------
int SmilesAtom::getAttachedHydrogenCount() const {
  return attached_hydrogens;
}

//-------------------------------------------------------------------------------------------------
int SmilesAtom::getFormalCharge() const {
  return formal_charge;
}

//-------------------------------------------------------------------------------------------------
int SmilesAtom::getIsotopicNumber() const {
  return isotopic_number;
}

//-------------------------------------------------------------------------------------------------
ChiralOrientation SmilesAtom::getChirality() const {
  return chirality;
}

//-------------------------------------------------------------------------------------------------
bool SmilesAtom::isAromatic() const {
  return is_aromatic;
}

//-------------------------------------------------------------------------------------------------
void SmilesAtom::setMajorBondOrders(const std::vector<int> &major_bond_orders_in) {
  major_bond_orders = major_bond_orders_in;
  updateAttachedHydrogenCount();
}

//-------------------------------------------------------------------------------------------------
void SmilesAtom::addMajorBond(const int order) {
  major_bond_orders.push_back(order);
}

//-------------------------------------------------------------------------------------------------
void SmilesAtom::updateAttachedHydrogenCount() {
  bool apply_bond_info = false;
  if (inferred_hydrogen_count) {
    switch (atomic_number) {
    case 6:
    case 7:
    case 8:
    case 9:
      attached_hydrogens = (10 - atomic_number) + formal_charge;
      apply_bond_info = true;
      break;
    case 14:
    case 15:
    case 16:
      {
        int bord_sum = 0;
        const size_t nbonds = major_bond_orders.size();
        for (size_t i = 0; i < nbonds; i++) {
          if (major_bond_orders[i] < 0) {
            is_aromatic = true;
            bord_sum += 1;
          }
        }
        if (bord_sum == 0) {
          attached_hydrogens = 2;
        }
        else if (bord_sum == 1) {
          attached_hydrogens = 1;
        }
        else {
          attached_hydrogens = 0;
        }
      }
      break;
    case 17:
      attached_hydrogens = 1 + formal_charge;
      apply_bond_info = true;
      break;
    case 34:
    case 35:
      attached_hydrogens = (18 - atomic_number) + formal_charge;
      apply_bond_info = true;
      break;
    default:
      attached_hydrogens = 0;
      break;
    }

    // Apply information about bond orders to adjust the number of hydrogens.
    if (apply_bond_info) {
      const size_t nbonds = major_bond_orders.size();
      for (size_t i = 0; i < nbonds; i++) {
        if (major_bond_orders[i] < 0) {

          // This would imply an error if the atom has an aromatic bond connecting it to something
          // else, yet the atom itself has not been marked as aromatic.  Take the order of the
          // bond to be one, and set the aromaticity regardless.
          is_aromatic = true;
          attached_hydrogens -= 1;
        }
        else {
          attached_hydrogens -= major_bond_orders[i];
        }
      }
    }

    // Subtract one attached hydrogen if all bond orders are explicitly given as one (or, -1, for
    // an aromatic bond "order") and the atom is aromatic.  If one of the bonds has already been
    // given order 2 (a Kekule representation), the number of attached hydrogens need not change.
    if (maxValue(major_bond_orders) == 1 && is_aromatic) {
      attached_hydrogens--;
    }
    
    // Check for unusual hydrogen counts
    if (attached_hydrogens < 0) {
      rtErr("An invalid number of attached hydrogens was obtained (" +
            std::to_string(attached_hydrogens) + ") for atom [" + bracketlessContent() + "].",
            "SmilesAtom");
    }
  }
}

//-------------------------------------------------------------------------------------------------
std::string SmilesAtom::bracketlessContent() const {
  std::string result;
  if (content.size() == 0) {
    return result;
  }
  if (content[0] == '[') {
    return content.substr(1, content.size() - 2);
  }
  else {
    return content;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
SmilesString::SmilesString(const std::string &basis_in) :
    basis{basis_in}, atom_count{0}, bond_count{0}, net_charge{0}, atomic_numbers{},
    isotopic_numbers{}, bonds{}, bond_orders{}, aromatic_bonds{}, bonded_partners{},
    partner_bounds{}, partner_bond_indices{}
{
  // Some baseline checks for sanity
  if (basis[0] == '(') {
    rtErr("The first character of a SMILES string cannot initiate a branch, as there is no atom "
          "to anchor the branch.", "SmilesString");
  }

  // Prepare a temporary array to hold bonds found during the search.  This int3 array holds the
  // index of the first atom in its "x" member, the index of the second atom in its "y" member, and
  // the order of the bond in its "z" member.
  std::vector<int3> major_bond_content;

  // Prepare an array of strings to hold the contents of individual atoms.  This will be assembled
  // in tandem with the temporary bond list, so as to have each (major) atom ready for
  // interpretation along with a definite list of the bonds connecting it to everything else.
  std::vector<std::string> major_atom_content;

  // Ring closures add another layer of detail, but provide a critical topological description
  // comprising a great deal of chemistry.  This array will track the atom indices closing a ring
  // in the "x" member and the ID number of the ring being closed in the "y" member (the numbering
  // begins at 1).  Connect the matching "x" members once the atom list is complete to add new
  // bonds and close all rings.
  std::vector<int3> ring_closures;
  
  // Iterate over the string, parsing branches.  Begin with the highest level and work inwards
  // through recursive function calls.
  this->parseBranch(basis, &major_atom_content, &major_bond_content, &ring_closures, -1, 0);
  
  // Apply ring closures.  Because ring closure ID numbers may be re-used, the convention is to
  // take the first such ID number that matches an "open" ring, after which the ID number can be
  // used to close another ring.
  const size_t n_closure = ring_closures.size();
  std::vector<bool> closure_coverage(ring_closures.size(), false);
  for (size_t i = 0; i < n_closure; i++) {
    if (closure_coverage[i]) {
      continue;
    }
    size_t j = i + 1;
    bool found = false;
    while (j < n_closure && (! found)) {
      if (ring_closures[i].y == ring_closures[j].y) {
        const int rcbo = std::max(ring_closures[i].z, ring_closures[j].z);
        major_bond_content.push_back({ ring_closures[i].x, ring_closures[j].x, rcbo });
        closure_coverage[j] = true;
        found = true;
      }
      j++;
    }
    if (j == n_closure && (! found)) {
      rtErr("No closure was found for ring ID " + std::to_string(ring_closures[i].y) + " in "
            "SMILES string " + basis, "SmilesString");
    }
  }

  // Create SmilesAtom class objects for each major atom identified in the original string.  Begin
  // by creating reverse lookup tables for every atom's connections, which will be reflexive
  // (double-counting each connection i => j and j => i).  For its nonstandard nature, this is not
  // done with indexingArray() from the series_ops library.
  const size_t n_major_atoms = major_atom_content.size();
  const size_t n_major_bonds = major_bond_content.size();
  std::vector<int> major_bond_partner_bounds(n_major_atoms + 1, 0);
  std::vector<int> major_bond_partners(2 * n_major_bonds);
  std::vector<int> major_bond_indices(2 * n_major_bonds);
  for (size_t i = 0; i < n_major_bonds; i++) {
    major_bond_partner_bounds[major_bond_content[i].x] += 1;
    major_bond_partner_bounds[major_bond_content[i].y] += 1;
  }
  prefixSumInPlace(&major_bond_partner_bounds, PrefixSumType::EXCLUSIVE);
  std::vector<int> major_bond_partner_counters = major_bond_partner_bounds;
  for (size_t pos = 0; pos < n_major_bonds; pos++) {
    int iatom = major_bond_content[pos].x;
    int jatom = major_bond_content[pos].y;
    const int iatom_pos = major_bond_partner_counters[iatom];
    const int jatom_pos = major_bond_partner_counters[jatom];
    major_bond_partners[iatom_pos] = jatom;
    major_bond_partners[jatom_pos] = iatom;
    major_bond_indices[iatom_pos] = pos;
    major_bond_indices[jatom_pos] = pos;
    major_bond_partner_counters[iatom] = iatom_pos + 1;
    major_bond_partner_counters[jatom] = jatom_pos + 1;
  }
  std::vector<SmilesAtom> major_atoms;
  major_atoms.reserve(n_major_atoms);
  for (size_t i = 0; i < n_major_atoms; i++) {
    const int j_llim = major_bond_partner_bounds[i];
    const int j_hlim = major_bond_partner_bounds[i + 1];
    std::vector<int> ibo(j_hlim - j_llim);
    for (int j = j_llim; j < j_hlim; j++) {
      ibo[j - j_llim] = major_bond_content[major_bond_indices[j]].z;
    }
    major_atoms.emplace_back(major_atom_content[i], ibo);
  }
  
  // Check for repeated bond specifications as a result of ring closure.  Recalculate the hydrogen
  // content if rings have changed the bonding patterns based on the first pass through the SMILES
  // string.
  if (n_closure > 0) {
    for (size_t i = 0; i < n_major_atoms; i++) {
      const int j_llim = major_bond_partner_bounds[i];
      const int j_hlim = major_bond_partner_bounds[i + 1];
      std::vector<int> updated_orders;
      updated_orders.reserve(j_hlim - j_llim);
      for (int j = j_llim; j < j_hlim; j++) {
        for (int k = j + 1; k < j_hlim; k++) {
          if (major_bond_partners[j] == major_bond_partners[k]) {
            rtErr("Ring closure has created a duplicate bond between major atoms " +
                  std::to_string(j) + " and " + std::to_string(k) + " in SMILES string " + basis);
          }
        }
        const int j_bond_idx = major_bond_indices[j];
        updated_orders.push_back(major_bond_content[j_bond_idx].z);
      }
      major_atoms[i].setMajorBondOrders(updated_orders);
    }
  }

  // Label each major bond (each bond between major atoms) as aromatic or non-aromatic.
  std::vector<bool> major_bond_aromaticity(n_major_bonds, false);
  for (size_t i = 0; i < n_major_atoms; i++) {
    if (major_atoms[i].isAromatic()) {
      const int j_llim = major_bond_partner_bounds[i];
      const int j_hlim = major_bond_partner_bounds[i + 1];
      for (int j = j_llim; j < j_hlim; j++) {
        const int jatom = major_bond_partners[j];
        if (major_atoms[jatom].isAromatic()) {

          // Look up the bond index and see whether it has already been flagged as aromatic.
          const int ij_bond_idx = major_bond_indices[j];
          if (major_bond_aromaticity[ij_bond_idx]) {
            continue;
          }
          
          // Atoms at both ends of the bond are aromatic, but this does not guarantee that the
          // bond itself is part of an aromatic system.  Consider, for example, a biphenyl
          // compound.  Check that there is a complete path through bonded aromatic atoms not
          // involving this bond.  This may sweep up some non-aromatic systems, but those will be
          // rare.
          detectAromaticBonds(i, jatom, j, major_atoms, major_bond_content, major_bond_indices,
                              major_bond_partners, major_bond_partner_bounds,
                              &major_bond_aromaticity);
        }
      }
    }
  }
  
  // Step through the major atoms and add their real atom content (including hydrogens) to the
  // master lists.
  atom_count = major_atoms.size();
  bond_count = major_bond_content.size();
  atomic_numbers.resize(atom_count);
  isotopic_numbers.resize(atom_count);
  formal_charges.resize(atom_count);
  for (size_t i = 0; i < n_major_atoms; i++) {
    const int nhyd = major_atoms[i].getAttachedHydrogenCount();
    net_charge += major_atoms[i].getFormalCharge();
    atom_count += nhyd;
    bond_count += nhyd;
  }
  int current_atom = 0;
  atomic_numbers.resize(atom_count);
  isotopic_numbers.resize(atom_count);
  formal_charges.resize(atom_count);
  std::vector<int> full_atom_from_major_atom(atom_count);
  std::vector<int> major_atom_is_full_atom(major_atoms.size());
  std::vector<bool> full_atom_is_added_hydrogen(atom_count, false);
  for (size_t i = 0; i < n_major_atoms; i++) {
    atomic_numbers[current_atom]   = major_atoms[i].getAtomicNumber();
    isotopic_numbers[current_atom] = major_atoms[i].getIsotopicNumber();
    formal_charges[current_atom]   = major_atoms[i].getFormalCharge();
    const int nhyd = major_atoms[i].getAttachedHydrogenCount();
    major_atom_is_full_atom[i] = current_atom;
    const int current_major_atom = current_atom;
    full_atom_from_major_atom[current_atom] = i;
    current_atom++;
    for (int j = 0; j < nhyd; j++) {
      atomic_numbers[current_atom] = 1;
      isotopic_numbers[current_atom] = 1;
      formal_charges[current_atom] = 0;
      full_atom_from_major_atom[current_atom] = i;
      full_atom_is_added_hydrogen[current_atom] = true;
      current_atom++;
    }
  }

  // With the atoms in place and maps established, the bonds can now be drawn.  The number of
  // bonds has been increased to include attached hydrogens.
  int current_bond = 0;
  bonds.resize(bond_count);
  bond_orders.resize(bond_count);
  aromatic_bonds.resize(bond_count);
  for (size_t i = 0; i < n_major_atoms; i++) {
    const int iatm_map = major_atom_is_full_atom[i];

    // Add bonds to attached hydrogens.  This will intersperse bonds to added hydrogens with other
    // bonds in the system, just as added hydrogens are interspersed throughout the atom order.
    // This is the best approach for STORMM's exclusion masking system, and produces molecular
    // structures that look more or less like the ordering of atoms in an AMBER topology.
    int htest_idx = iatm_map + 1;
    while (htest_idx < atom_count && full_atom_from_major_atom[htest_idx] == i) {
      bonds[current_bond] = { iatm_map, htest_idx };
      bond_orders[current_bond] = 1;
      aromatic_bonds[current_bond] = false;
      current_bond++;
      htest_idx++;
    }
  
    // Add all bonds to other major atoms emanating from this one, if the major atom index of the
    // other atom is higher.
    for (int j = major_bond_partner_bounds[i]; j < major_bond_partner_bounds[i + 1]; j++) {
      if (major_bond_partners[j] > i) {
        const int jatm_map = major_atom_is_full_atom[major_bond_partners[j]];
        const int major_bond_idx = major_bond_indices[j];
        bonds[current_bond] = {  iatm_map, jatm_map };
        bond_orders[current_bond] = bond_orders[major_bond_idx];
        aromatic_bonds[current_bond] = major_bond_aromaticity[major_bond_idx];
        current_bond++;
      }
    }
  }

  // Produce tables, similar to major_bond_partners and its bounds array, of all bonds in the
  // system with all hydrogens included.
  bonded_partners.resize(2 * bond_count);
  partner_bounds.resize(atom_count + 1, 0);
  partner_bond_indices.resize(2 * bond_count);
  for (int i = 0; i < bond_count; i++) {
    partner_bounds[bonds[i].x] += 1;
    partner_bounds[bonds[i].y] += 1;
  }
  prefixSumInPlace(&partner_bounds, PrefixSumType::EXCLUSIVE);
  std::vector<int> pcounters = partner_bounds;
  for (int pos = 0; pos < bond_count; pos++) {
    const int iatom = bonds[pos].x;
    const int jatom = bonds[pos].y;
    const int iatom_next = pcounters[iatom];
    const int jatom_next = pcounters[jatom];
    bonded_partners[iatom_next] = jatom;
    bonded_partners[jatom_next] = iatom;
    partner_bond_indices[iatom_next] = pos;
    partner_bond_indices[jatom_next] = pos;
    pcounters[iatom] = iatom_next + 1;
    pcounters[jatom] = jatom_next + 1;
  }
}

//-------------------------------------------------------------------------------------------------
SmilesString::SmilesString(const char* basis_in) :
    SmilesString(std::string(basis_in))
{}

//-------------------------------------------------------------------------------------------------
int SmilesString::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
int SmilesString::getBondCount() const {
  return bond_count;
}

//-------------------------------------------------------------------------------------------------
int SmilesString::getBondCount(const int index) const {
  validateAtomIndex(index);
  return partner_bounds[index + 1] - partner_bounds[index];
}

//-------------------------------------------------------------------------------------------------
int SmilesString::getAtomicNumber(const int index) const {
  validateAtomIndex(index);
  return atomic_numbers[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& SmilesString::getAtomicNumbers() const {
  return atomic_numbers;
}

//-------------------------------------------------------------------------------------------------
int SmilesString::getIsotopicNumber(const int index) const {
  validateAtomIndex(index);
  return isotopic_numbers[index];
}

//-------------------------------------------------------------------------------------------------
int SmilesString::getBondOrder(const int index_a, const int index_b) const {
  validateAtomIndex(index_a);
  validateAtomIndex(index_b);
  int bond_entry_id = -1;
  for (int i = partner_bounds[index_a]; i < partner_bounds[index_a + 1]; i++) {
    if (bonded_partners[i] == index_b) {
      bond_entry_id = partner_bond_indices[i];
    }
  }
  if (bond_entry_id == -1) {
    return 0;
  }
  else {
    return (aromatic_bonds[bond_entry_id]) ? -1 : bond_orders[bond_entry_id];
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int SmilesString::getBondedPartnerCount(const int index) const {
  validateAtomIndex(index);
  return partner_bounds[index + 1] - partner_bounds[index];
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SmilesString::getBondedPartnerIndices(const int index) const {
  validateAtomIndex(index);
  const int llim = partner_bounds[index];
  const int hlim = partner_bounds[index + 1];
  std::vector<int> result(hlim - llim);
  for (int i = llim; i < hlim; i++) {
    result[i - llim] = bonded_partners[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const std::string& SmilesString::getBasis() const {
  return basis;
}
  
//-------------------------------------------------------------------------------------------------
void SmilesString::validateAtomIndex(const int index) const {
  if (index < 0 || index >= atom_count) {
    rtErr("Atom index " + std::to_string(index) + " is invalid for a collection of " +
          std::to_string(atom_count) + " atoms.", "SmilesString", "validateAtomIndex");
  }
}

//-------------------------------------------------------------------------------------------------
std::string SmilesString::parseBracketedAtom(const std::string &content, int *parse_pos) {
  std::string result;
  const int max_pos = content.size();
  int lpp = *parse_pos + 1;
  while (lpp < max_pos && content[lpp] != ']') {
    lpp++;
  }
  if (lpp == max_pos) {
    rtErr("An atom entry lacks a closing bracket in \"" + content + "\".", "SmilesString",
          "parseAtom");
  }
  else if (lpp == *parse_pos + 1) {
    rtErr("An empty atom was found in \"" + content + "\".", "SmilesString", "parseAtom");
  }
  result.reserve(lpp - *parse_pos - 1);
  lpp = *parse_pos + 1;
  while (content[lpp] != ']') {
    result.push_back(content[lpp]);
    lpp++;
  }
  *parse_pos = lpp + 1;
  return result;
}

//-------------------------------------------------------------------------------------------------
void SmilesString::parseBranch(const std::string &content,
                               std::vector<std::string> *atom_content,
                               std::vector<int3> *bond_content, std::vector<int3> *ring_closures,
                               const int branch_initiating_atom, const int anchor_bond_order) {
  const int max_pos = content.size();
  int parse_pos = 0;
  int previous_atom = branch_initiating_atom;
  int current_atom;
  int confirmed_anchor_bond_order = anchor_bond_order;
  int next_bond_order = anchor_bond_order;
  const int zero_int = '0';
  while (parse_pos < max_pos) {

    // Assume that the bond to the next atom is a single bond.  Special characters will indicate
    // if it is anything else.
    if (content[parse_pos] == '[') {
      atom_content->push_back(this->parseBracketedAtom(content, &parse_pos));

      // If a branch comes next, this new atom will be the initiator.
      current_atom = atom_content->size() - 1;
      if (previous_atom >= 0) {
        bond_content->push_back({ previous_atom, current_atom, next_bond_order });
      }
      previous_atom = current_atom;
      next_bond_order = 1;
    }
    else if (content[parse_pos] == 'C') {

      // Carbon and chlorine must be distinguished
      if (parse_pos + 1 < max_pos && content[parse_pos + 1] == 'l') {
        atom_content->push_back("Cl");
        parse_pos += 2;
      }
      else {
        atom_content->push_back("C");
        parse_pos++;
      }
      current_atom = atom_content->size() - 1;
      if (previous_atom >= 0) {
        bond_content->push_back({ previous_atom, current_atom, next_bond_order });
      }
      previous_atom = current_atom;
      next_bond_order = 1;
    }
    else if (content[parse_pos] == 'B') {

      // Boron and bromine must be distinguished
      if (parse_pos + 1 < max_pos && content[parse_pos + 1] == 'r') {
        atom_content->push_back("Br");
        parse_pos += 2;
      }
      else {
        atom_content->push_back("B");
        parse_pos++;
      }
      current_atom = atom_content->size() - 1;
      if (previous_atom >= 0) {
        bond_content->push_back({ previous_atom, current_atom, next_bond_order });
      }
      previous_atom = current_atom;
      next_bond_order = 1;
    }
    else if (content[parse_pos] == 'N' || content[parse_pos] == 'O' || content[parse_pos] == 'P' ||
             content[parse_pos] == 'S' || content[parse_pos] == 'F' || content[parse_pos] == 'I' ||
             content[parse_pos] == 'b' || content[parse_pos] == 'c' || content[parse_pos] == 'n' ||
             content[parse_pos] == 'o' || content[parse_pos] == 'p' || content[parse_pos] == 's') {

      // Six more atoms in the "organic" subset can be logged based on their single-character
      // symbols in the periodic table.  In addition, the aromatic cases of various elements in the
      // organic subset must be recognized.  For the aromatic atoms, the bond orders will remain
      // set to 1, pending adjustments once the 
      atom_content->push_back(std::string(1, content[parse_pos]));
      parse_pos++;
      current_atom = atom_content->size() - 1;
      if (previous_atom >= 0) {
        bond_content->push_back({ previous_atom, current_atom, next_bond_order });
      }
      previous_atom = current_atom;
      next_bond_order = 1;
    }
    else if (content[parse_pos] == '-') {
      next_bond_order = 1;
      parse_pos++;
    }
    else if (content[parse_pos] == '=') {

      // Modify the bond to the next atom, or make a post-hoc fix to the anchoring bond if this is
      // the first character in the branch.
      next_bond_order = 2;
      parse_pos++;
    }
    else if (content[parse_pos] == '#') {
      next_bond_order = 3;
      parse_pos++;
    }
    else if (content[parse_pos] == '$') {
      next_bond_order = 4;
      parse_pos++;
    }
    else if (content[parse_pos] == ':') {

      // The next bond is said to be an aromatic bond. (Whether it is a true part of an aromatic
      // system will be reinterpreted by the ChemicalFeatures class if and when the SmilesString
      // is translated to produce such an object.)
      next_bond_order = -1;
      parse_pos++;
    }
    else if (content[parse_pos] == '.') {

      // There is no bond between the current atom and the next atom.  This will have implications
      // for the number of attached hydrogens.
      next_bond_order = 0;
      parse_pos++;
    }
    else if (content[parse_pos] >= '0' && content[parse_pos] <= '9') {
      const int number_init_pos = parse_pos;
      parse_pos++;

      // Detect a number indicative of ring closure.
      while (parse_pos < max_pos && content[parse_pos] >= '0' && content[parse_pos] <= '9') {
        parse_pos++;
      }
      const int number_final_pos = parse_pos;
      int mult = 1;
      int ring_number = 0;
      for (int i = number_final_pos - 1; i >= number_init_pos; i--) {
        ring_number += (content[i] - zero_int) * mult;
        mult *= 10;
      }
      ring_closures->push_back({ current_atom, ring_number, next_bond_order });
      next_bond_order = 1;
    }
    else if (content[parse_pos] == '(') {
      parse_pos++;
      int left_depth = 1;
      int right_depth = 0;
      const int branch_init = parse_pos;
      while (right_depth < left_depth && parse_pos < max_pos) {
        left_depth  += (content[parse_pos] == '(');
        right_depth += (content[parse_pos] == ')');
        parse_pos++;
      }
      if (right_depth < left_depth) {
        rtErr("An unterminated branch was found in SMILES string segment \"" + content + "\".",
              "SmilesString", "parseBranch");
      }
      const int branch_end = parse_pos - 1;
      this->parseBranch(content.substr(branch_init, branch_end - branch_init), atom_content,
                        bond_content, ring_closures, current_atom, next_bond_order);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void SmilesString::detectAromaticBonds(const int iatom, const int jatom, const int ij_setting,
                                       const std::vector<SmilesAtom> &major_atoms,
                                       const std::vector<int3> &major_bonds,
                                       const std::vector<int> &major_bond_indices,
                                       const std::vector<int> &major_bond_partners,
                                       const std::vector<int> &major_bond_partner_bounds,
                                       std::vector<bool> *major_bond_aromaticity) {

  // The "path" is both a list of major atom indices as well as the progess searching through each
  // major atom's bond partners for the next element.  The result of the search is to close the
  // path if the first major atom (given by the iatom input parameter) can be reached from the
  // third, fourth, or higher numbered atoms along the path without backtracking.
  std::vector<int3> path;
  path.reserve(16);
  path.push_back({ iatom, ij_setting, major_bond_partner_bounds[iatom + 1] });
  path.push_back({ jatom,
                   major_bond_partner_bounds[jatom],
                   major_bond_partner_bounds[jatom + 1] });
  bool path_open = true;
  bool path_viable = true;
  int path_length = 2;
  int path_head = 1;
  while (path_open && path_viable) {
    const int3 current_atom = path.back();
    int search_pos = current_atom.y;

    // If all options for continuing the path beyond its second atom have been exhausted, the path
    // is no longer viable.
    if (path.size() == 2 && search_pos >= current_atom.z) {
      path_viable = false;
      continue;
    }
    
    // Continue searching.
    while (search_pos < current_atom.z) {

      // If the "bond" is actually of order zero, there is no connection.  Go on.
      const int sb_idx = major_bond_indices[search_pos];
      if (major_bonds[sb_idx].z == 0) {
        search_pos++;
        continue;
      }
      
      // If the path is long enough and the bond leads back to the start, this is a closed path.
      if (path.size() > 2 && major_bond_partners[search_pos] == path[0].x) {
        path.back().y = search_pos;
        path.push_back({ iatom,
                         major_bond_partner_bounds[iatom],
                         major_bond_partner_bounds[iatom + 1] });
        path_open = false;
        break;
      }

      // If a new and unique atom can be added to the path, keep going.
      bool backtracking = false;
      const size_t path_len = path.size();
      for (size_t i = 1; i < path_len; i++) {
        backtracking = (backtracking || (major_bond_partners[search_pos] == path[i].x));
      }
      if (backtracking) {
        search_pos++;
      }
      else {
        const int next_atom = major_bond_partners[search_pos];
        path.back().y = search_pos;
        path.push_back({ next_atom,
                         major_bond_partner_bounds[next_atom],
                         major_bond_partner_bounds[next_atom + 1] });
        break;
      }
    }

    // If the search has been exhausted for the current path head, remove that atom, revert to
    // the previous atom, and continue searching.
    if (search_pos == current_atom.z) {
      path.pop_back();
      path.back().y += 1;
    }
  }

  // If a viable path was found, mark all bonds along the path as aromatic.
  if (path_viable) {
    const size_t path_len = path.size();
    for (size_t i = 0; i < path_len; i++) {
      const int ib_idx = major_bond_indices[path[i].y];
      major_bond_aromaticity->at(ib_idx) = true;
    }
  }
}
  
} // namespace chemistry
} // namespace stormm
