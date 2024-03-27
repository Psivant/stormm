// -*-c++-*-
#ifndef STORMM_FORWARD_EXCLUSIONMASK_H
#define STORMM_FORWARD_EXCLUSIONMASK_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/stormm_vector_types.h"
#include "Topology/atomgraph.h"

namespace stormm {
namespace energy {

using card::Hybrid;
using card::HybridTargetLevel;
using topology::AtomGraph;

/// \brief The abstract for the forward exclusion mask turns out to be shorter than that for the
///        StaticExclusionMask.
struct ForwardExclusionMaskReader {

  /// \brief The constructor works as expected to translate a list of input arguments into the
  ///        object's member variables.
  ForwardExclusionMaskReader(int natom_in, int mask2_size_in, const int* primary_idx_in,
                             const uint2* mask1_in, const uint2* mask2_in);

  const int natom;         ///< The number of atoms
  const int mask2_size;    ///< Size of the secondary masks (mask2) array
  const int* primary_idx;  ///< Index of the primary mask for each atom 
  const uint2* mask1;      ///< The primary masks for each atom in the system
  const uint2* mask2;      ///< Secondary masks needed by any atoms in the system
};

/// \brief This exclusion mask relies on the reality that, in most topologies, the bulk of one
///        atom's exclusions will be found within the next 32 atoms listed in the topology.  If
///        an atom excludes other atoms beyond this range, a separate index records the position
///        at which this happens and provides either an index into a list of additional atoms or,
///        if possible, the index of one other reference index in the topology followed by a
///        second mask of 32 atoms which indicate exclusions relative to the reference index.  To
///        keep the list as compact as possible, the exclusion mask only considers "forward"
///        exclusions, meaning excluded atoms with higher index than the one for which the mask
///        is made.  There is no double-counting: if atom 1 excludes atom 5, then atom 5 does not
///        note an exclusion of atom 1.
class ForwardExclusionMask {
public:

  /// \brief The constructor requires a topology and creates a blank object if given nullptr.
  /// \{
  ForwardExclusionMask(const AtomGraph *ag_in = nullptr);
  ForwardExclusionMask(const AtomGraph &ag_in);
  /// \}

  /// \brief The default copy and move constructors as well as the copy and move assignment
  ///        operator will suffice for this object, which has no POINTER-kind Hybrid objects among
  ///        its members.
  /// \{
  ForwardExclusionMask(const ForwardExclusionMask &original) = default;
  ForwardExclusionMask(ForwardExclusionMask &&original) = default;
  ForwardExclusionMask& operator=(const ForwardExclusionMask &other) = default;
  ForwardExclusionMask& operator=(ForwardExclusionMask &&other) = default;
  /// \}

  /// \brief Get the number of atoms in the system.
  int getAtomCount() const;

  /// \brief Get the number of unique primary masks that the system requires to store all of its
  ///        non-bonded exclusions.
  int getPrimaryMaskCount() const;

  /// \brief Get the number of primary masks that involve a half-mask extension rather than a
  ///        reference to a full list of additional secondary masks.
  int getExtendedMaskCount() const;

  /// \brief Get the number of secondary masks that the system requires to store all of its
  ///        non-bonded exclusions.
  int getSecondaryMaskCount() const;

  /// \brief Get the total number of exclusions by counting all the checked bits and accounting
  ///        for excluded reference atoms.
  int getTotalExclusionsCount() const;

  /// \brief Get an atom's primary mask
  ///
  /// \param index  Topological index of the atom of interest
  uint2 getPrimaryMask(int index) const;

  /// \brief Get all secondary masks associated with an atom.  Each presents the exclusions for a
  ///        cluster of atoms following some reference atom.
  ///
  /// \param start   The index of the first secondary mask to read
  /// \param length  The number of secondary masks to read
  std::vector<uint2> getSecondaryMasks(int start, int length) const;

  /// \brief Get the collection of this object's pointers and critical constants
  ///
  /// \param tier  Level at which to get the data pointers (host CPU or GPU device)
  const ForwardExclusionMaskReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Test whether a pair of atoms is an exclusion or not
  ///
  /// \param atom_i  The first atom
  /// \param atom_j  The second atom (if atom_j < atom_i, the function will make a recursive call
  ///                with atom_i and atom_j reversed)
  bool testExclusion(int atom_i, int atom_j) const;

private:

  /// The number of atoms in the system
  int atom_count;

  /// Indices into the primary_masks array for each atom
  Hybrid<int> primary_mask_indices;

  /// The method of first resort: every atom i in the topology gets an index into this list.  The
  /// first component of the tuple (the x member) is a mask of 32 bits indicating whether atoms
  /// i + 1, i + 2, ..., i + 32 are exclusions (1) or not (0).  The tuple's second component (its y
  /// member) is a bit-packed value providing an indicator of whether there is a single cluster of
  /// additional exclusions or if there are multiple clusters of additional exclusions.  A single
  /// cluster of additional exclusions is indicated by a (0) in the y-tuple's highest bit (31), and
  /// means that bits 0 through 14 are an unsigned integer P such that the atom with topological
  /// index (i + P) is excluded by atom i, and that bits 15 through 30 indicate up to 16
  /// consecutive atoms following (i + P) which may also be exclusions of atom i (exclusions
  /// indicated by placing a (1) in the bit).  Multiple clusters of additional exclusions are
  /// indicated by placing a (1) in the y member's bit 31, in which case bits 0 through 21 of the
  /// tuple's y-member are read as an unsigned integer S and bits 22 through 30 are read as an
  /// unsigned integer T.  The program then knows to search elements [S ... S + T) of the
  /// secondary_masks array for the absolute topological indices of reference atoms and clusters of
  /// up to 32 atoms that follow them.  In this scheme, an atom may exclude up to 512 clusters of
  /// up to 33 atoms and up to 4,194,304 additional clusters of exclusions may be tracked.
  ///
  /// It is very difficult to even conceive of a system that might break this format, as it would
  /// have to be broken by the sheer size of a system with a wild arrangement of excluded atoms
  /// throughout.  Disulfide bridges in standard models of amino acids can be captured within the
  /// confines of the primary_masks array.  A protein of perhaps 1,800 residues with a circular
  /// disulfide link between its N- and C-terminal residues might exceed the 32,768 relative atom
  /// index limit encoded by bits 0 through 14 and thereby generate one entry in the
  /// secondary_masks array.  Very large assemblies of multiple proteins will not exceed the limits
  /// of the primary_masks array.  Solvent molecules of more than 33 atoms with circular bonded
  /// arrangements might each generate a single entry in the secondary_masks array, which might
  /// imply an upper particle limit somewhere in the range of 142,000,000 atoms, although it is
  /// hard to imagine a solvent molecule with such characteristics.  A situation where exclusion
  /// mapping is thoroughly scrambled and one atom bonds to eight neighbors, which each bond to
  /// eight of their own neighbors, which each bond to eight more of their own neighbors, could
  /// break the scheme.  But STORMM will just refuse to simulate that system.
  Hybrid<uint2> primary_masks;

  /// Secondary exclusion masks.  If this array is needed at all, which will be extremely rare in
  /// simulations of most biomolecules, the y-member of an atom's primary_masks tuple will provide
  /// the starting index and length of the array to access.
  Hybrid<uint2> secondary_masks;

  /// Pointer to the original topology
  AtomGraph *ag_pointer;
};

}
}

#endif
