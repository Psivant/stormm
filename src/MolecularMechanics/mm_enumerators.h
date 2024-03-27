// -*-c++-*-
#ifndef STORMM_MM_ENUMERATORS_H
#define STORMM_MM_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace mm {

/// \brief Enumerate the dimensions of a molecular mechanics problem: atoms, bonds, angles, other
///        valence terms, numbers of non-bonded tiles (if a fixed, all-to-all system),
///        Lennard-Jones atom types, starting indices of the Lennard-Jones tables for the system,
///        and restraints.
enum class ValenceWorkUnitSpecs {
  ATOM_INDEX_START = 0,   ///< Starting position at which atom indices for this valence work unit
                          ///<   or conformer begin.  If negative, this implies that it is not
                          ///<   necessary to read a list of atom indices and then get the actual
                          ///<   atoms.  Rather, atoms in a range starting at
                          ///<   ATOM_DIRECT_READ_START and continuing for the next continguous
                          ///<   ATOM_COUNT atoms are to be read.
  ATOM_DIRECT_READ_START, ///< Starting position for reading a contiguous array of ATOM_COUNT atoms
                          ///<   that participate in the work unit.  If negative, this implies that
                          ///<   the work unit must read indices for atoms beginning at an array
                          ///<   index given by ATOM_INDEX_START and obtain the next ATOM_COUNT
                          ///<   entries in order to know the indices of actual atoms with which to
                          ///<   build the system.  ATOM_DIRECT_READ_START and ATOM_INDEX_START are
                          ///<   mutually exclusive, and although different STORMM applications
                          ///<   will look to one or the other by construction, only one of these
                          ///<   enumerations should be set.
  BOND_INSR_START,        ///< Starting position at which to begin reading bond instructions
  ANGL_INSR_START,        ///< Starting position for three-atom harmonic angle instructions
  DIHE_INSR_START,        ///< Starting position for cosine-based dihedral instructions
  ROGUE_14_INSR_START,    ///< Starting position for inferred / free 1:4 attenuated pair
                          ///<   interactions (non-bonded exclusions) which cannot be assigned to
                          ///<   any particular dihedral term
  UBRD_INSR_START,        ///< Starting position for Urey-Bradley harmonic angle instructions
  CIMP_INSR_START,        ///< Starting position for CHARMM improper instructions
  CMAP_INSR_START,        ///< Starting position for CMAP instructions
  RESTRAINT_INSR_START,   ///< Starting position for restraint instructions
  CONSTRAINT_INSR_START,  ///< Starting positions for constraint group (non-SETTLE) instructions
  SETTLE_INSR_START,      ///< Starting positions for SETTLE group instructions
  VSITE_INSR_START,       ///< Starting positions for virtual site instructions
  MOVE_MASK_START,        ///< Starting position for the movement mask.  Each bit in this mask
                          ///<   indicates whether the corresponding atom in the work unit is to
                          ///<   have its position updated once forces from all sources have been
                          ///<   accumulated.  The length of the bit mask is implied by the number
                          ///<   of atoms in the work unit.
  LJ_TABLE_START,         ///< Starting positions for Lennard Jones non-bonded tables that this
                          ///<   work unit will access.  Most valence work units will not access
                          ///<   the non-bonded Lennard-Jones parameters, they will make use of the
                          ///<   1:4 Lennard-Jones parameter tables instead.  However, small
                          ///<   molecule conformations will perform short-ranged non-bonded
                          ///<   computations and use these tables.  The Lennard-Jones A and B
                          ///<   coefficients are fused into a single table, while the C
                          ///<   coefficients, if they exist, will be packed into a separate table.
  LJ_14_TABLE_START,      ///< Starting positions for Lennard-Jones 1:4 non-bonded tables.  For
                          ///<   Amber force fields, the 1:4 Lennard-Jones tables carry the same
                          ///<   values as the standard Lennard-Jones tables (with or without
                          ///<   NB-Fix), but for CHARMM topologies the 1:4 Lennard-Jones
                          ///<   interactions are independent of their non-bonded counterparts.
  NONBONDED_TILE_START,   ///< Starting index for specifications of non-bonded interaction tiles,
                          ///<   if applicable (i.e. small molecule conformations where the valence
                          ///<   and non-bonded computations are fused into the same kernel)
  ATOM_COUNT,             ///< Number of atoms in this work unit (or small molecule conformation--a
                          ///<   small molecule conformation is essentially one work unit)
  BOND_COUNT,             ///< Number of bond instructions in this work unit
  ANGL_COUNT,             ///< Number of angle instructions in this work unit
  DIHE_COUNT,             ///< Number of cosine-based dihedral instructions in this work unit
                          ///<   (these dihedrals often imply 1:4 interactions and subsume the vast
                          ///<   majority of the attenuated pairs, but some free pairs may remain)
  ROGUE_14_COUNT,         ///< Number of free 1:4 interactions (non-bonded exclusions) which must
                          ///<   be handled on their own
  UBRD_COUNT,             ///< Number of Urey-Bradley harmonic angle instructions in this work unit
  CIMP_COUNT,             ///< Number of CHARMM improper instructions in this work unit
  CMAP_COUNT,             ///< Number of CMAP instructions in this work unit
  RESTRAINT_COUNT,        ///< Number of restraints in the system (restraints for 1, 2, 3, or 4
                          ///<   atoms are all just restraints--they will be arranged in the list
                          ///<   with padding such that all restraints of similar orders are
                          ///<   grouped together to limit warp divergence.
  CONSTRAINT_GROUP_COUNT, ///< Number of bond length constraint groups (SHAKE or RATTLE)
  SETTLE_GROUP_COUNT,     ///< Number of SETTLE groups (fast rigid waters)
  VSITE_COUNT,            ///< Number of virtual sites in the work unit
  LJ_TYPE_COUNT,          ///< Number of unique Lennard-Jones atom types, defining the size of the
                          ///<   Lennard-Jones interaction matrix for this work unit or
                          ///<   conformation
  NONBONDED_TILE_COUNT,   ///< Number of non-bonded tiles that this work unit must perform,
                          ///<   computing pair interactions for all particles to all others in
                          ///<   each tile.  This is only relevant to small molecule conformers,
                          ///<   not valence work units serving a system with periodic boundaries
                          ///<   and a namelist or a large system with Generalized Born solvent
                          ///<   (for which the non-bonded computations would comprise their own
                          ///<   kernel).
  DESCRIPTOR_COUNT        ///< Total number of descriptors in this enumerator
};

/// \brief Produce a human-readable string to describe each enumeration.  Overloads of this
///        function are found in other libraries as well.
///
/// \param input  The enumerated value to translate
/// \{
std::string getEnumerationName(ValenceWorkUnitSpecs input);
/// \}

} // namespace mm
} // namespace stormm

#endif
