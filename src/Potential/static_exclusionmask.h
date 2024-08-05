// -*-c++-*-
#ifndef STORMM_SIMPLE_EXCLUSIONMASK_H
#define STORMM_SIMPLE_EXCLUSIONMASK_H

#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/common_types.h"
#include "Topology/atomgraph.h"

namespace stormm {
namespace energy {

using card::Hybrid;
using card::HybridTargetLevel;
using topology::AtomGraph;

/// \brief The maxmimum number of supertiles, obtained from the limit imposed by the maximum value
///        of a 32-bit signed int.  Translates to just over 11.8 million atoms.
constexpr int max_supertile_count = 46300;

/// \brief The maximum number of unique tiles whose masks could be indexed with a 32-bit signed it
constexpr int max_unique_tiles = 67108864;

/// \brief The number of atoms in a supertile stride
constexpr int supertile_length = 256;

/// \brief The number of atoms in a tile stride (tiles are tile_length x tile_length atoms)
/// \{
constexpr int tile_length = 16;
constexpr int half_tile_length = tile_length / 2;
constexpr int three_halves_tile_length = (tile_length * 3) / 2;
constexpr int tile_length_bits_mask = tile_length - 1;
constexpr int tile_length_bits = 4;
/// \}

/// \brief The number of tile lengths that fit into one supertile length
constexpr int tile_lengths_per_supertile = supertile_length / tile_length;

/// \brief The number of square tiles that fit into a full, square supertile
constexpr int tiles_per_supertile = tile_lengths_per_supertile * tile_lengths_per_supertile;

/// \brief The abstract for a StaticExclusionMask object, read-only due to the const-ness of the
///        data() member function that returns it.
struct StaticExclusionMaskReader {

  /// Constructor works as any other abstract, taking a list of relevant constants and pointers
  StaticExclusionMaskReader(int natom_in, int supertile_stride_count_in,
                            const int* supertile_map_idx_in, const int* tile_map_idx_in,
                            const uint* mask_data_in);

  const int natom;                   ///< Number of atoms in the system
  const int supertile_stride_count;  ///< Number of supertiles in the system
  const int* supertile_map_idx;      ///< Map indices for each supertile's stretch of tile masks.
                                     ///<   The indices in this array are pre-inflated by 256 (the
                                     ///<   number of tiles in a supertile)
  const int* tile_map_idx;           ///< Map indices for each tile's stretch of atoms within a
                                     ///<   supertile.  These are pre-inflated by 32, the length
                                     ///<   of an actual tile's mask.
  const uint* mask_data;             ///< The actual exclusion mask data
};

/// \brief A simple pair list for an all-to-all calculation with exclusion masks.  The list stores
///        masks for 16 x 16 tiles of atoms indexed according to the atom order in some original
///        topology.  Each 16 x 16 tile mask will store the exclusions of 16 sending and 16
///        receiving atoms in the low 16 bits of an unsigned integer (probably 32-bit, but 16-bit
///        format is also acceptable).  The 16 x 16 tile masks are, in turn, arranged in groups of
///        16 x 16 such that each group of tiles covers the interactions of two consecutive
///        sequences of 256 atoms apiece, for a total of 65,536 unique bits in 32kB (which counts
///        the wasted high bits in each sending atom's mask as well as the reserved nature of the
///        second sixteen bits).  Many groups will have no exclusions at all, and these will all
///        index into the first group stored in a master list.  In this manner, a system of N atoms
///        can have a complete pairlist constructed with ((N + 255) / 256)^2 super-tiles, with a
///        storage requirement growing as N^2 but with a very small prefactor.  With further
///        hierarchical ordering, the pair list could even become O(N log N) in storage
///        requirements.  In the present scheme, simple pair lists can be built for very large
///        systems if the exclusions are not spread throughout: a system of 131072 atoms with
///        exclusions only between atoms +/-255 indices from one another would have space
///        requirements of at most 52MB, probably much less.  It is far better than the 2GB which
///        would be required to store every exclusion as a unique bit.
class StaticExclusionMask {
public:

  /// \brief The constructor requires a topology, which will have lists of exclusions
  ///
  /// \param ag_in  The topology to base this all-to-all pair list upon (fed in as a pointer so
  ///               that a private member variable may retain a pointer to the original topology)
  /// \{
  StaticExclusionMask(const AtomGraph *ag_in = nullptr);
  StaticExclusionMask(const AtomGraph &ag_in);
  /// \}
  
  /// \brief The default copy and move constructors as well as the copy assignment operator will
  ///        suffice for this object, which has no POINTER-kind Hybrid objects among its members.
  /// \{
  StaticExclusionMask(const StaticExclusionMask &original) = default;
  StaticExclusionMask(StaticExclusionMask &&original) = default;
  StaticExclusionMask& operator=(const StaticExclusionMask &other) = default;
  StaticExclusionMask& operator=(StaticExclusionMask &&other) = default;
  /// \}
  
  /// \brief Get the total atom count
  int getAtomCount() const;

  /// \brief Get the number of supertile strides that the system contains (how many groups of up
  ///        to 256 atoms will this system break into?)
  int getSuperTileStrideCount() const;

  /// \brief Get the number of supertile strides that the system contains (how many groups of up
  ///        to 16 atoms will this system break into?)
  int getTileStrideCount() const;

  /// \brief Get the total number of unique mask super-tiles
  int getUniqueSuperTileCount() const;

  /// \brief Get the total number of unique tiles within a given super-tile
  ///
  /// Overloaded:
  ///   - Get the total number of unique 16 x 16 atom tiles in all super-tiles (no arguments)
  ///   - Get the number of unique 16 x 16 atom tiles in a particular super-tile
  ///
  /// \param supertile_i_index  The ith group of 256 atoms is sending to collect interactions
  /// \param supertile_j_index  The jth group of 256 atoms is receiving to collect interactions
  /// \{
  int getUniqueTileCount() const;
  int getUniqueTileCount(int supertile_i_index, int supertile_j_index) const;
  /// \}

  /// \brief Get the exclusion masks for sending atoms in a particular tile
  ///
  /// \param supertile_i_index  The ith stretch of 256 atoms is interacting with the jth
  /// \param supertile_j_index  The jth stretch of 256 atoms is interacting with the ith 
  /// \param tile_i_index       The ith group of 16 atoms within the supertile is sending to
  ///                           collect interactions
  /// \param tile_j_index       The jth group of 16 atoms within the supertile is receiving to
  ///                           collect interactions
  std::vector<uint> getTileExclusions(int supertile_i_index, int supertile_j_index,
                                      int tile_i_index, int tile_j_index) const;

  /// \brief Return a pointer to the topology that built this exclusion mask
  const AtomGraph* getTopologyPointer() const;

  /// \brief Get the abstract for this exclusion mask object
  const StaticExclusionMaskReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a pointer ot the object itself (useful when passing the object by const reference)
  const StaticExclusionMask* getSelfPointer() const;
  
#ifdef STORMM_USE_HPC
  /// \brief Upload the object's information to the GPU device
  void upload();

  /// \brief Download the object's information from the GPU device
  void download();
#endif
  
  /// \brief Test whether some atom pair is an exclusion
  ///
  /// \param atom_i  The ith atom to test for a pair exclusion with atom_j
  /// \param atom_j  The jth atom to test for a pair exclusion with atom_i
  bool testExclusion(int atom_i, int atom_j) const;

private:
  int atom_count;               ///< The number of atoms in the system
  int supertile_stride_count;   ///< The number of groups of 256 atoms that can be defined for the
                                ///<   system.  (atom_count + 255) / 256.
  int tile_stride_count;        ///< The number of groups of 16 atoms that can be defined for the
                                ///<   system.  (atom_count + 15) / 16.
  int unique_supertile_count;   ///< The number of unique super-tiles.  The empty supertile is
                                ///<   unique and every unique super-tile has its own unique, empty
                                ///<   tile.
  int unique_tile_count;        ///< The number of unique 16 x 16 tiles.  Tile (i, j) is different
                                ///<   from tile (j, i) as the sending and receiving atoms are
                                ///<   reversed.  If this counter overflows, the system is likely
                                ///<   huge (at least 720k atoms, with at least one excluded pair
                                ///<   between any two groups of 16 atoms, which would be extreme).
                                ///<   Overflow in this counter will be trapped, but is unlikely to
                                ///<   occur before other aspects of the code hit practical
                                ///<   limitations on the system memory.

  // A collection of the mask groups, arranged in a concatenated fashion
  Hybrid<uint> all_masks;             ///< The actual mask data.  Each tile's exclusion map is
                                      ///<   a stretch of 32 unsigned integers in this array (16
                                      ///<   for the tile's sending atoms and 16 for the tile's
                                      ///<   receiving atoms).
  Hybrid<int> supertile_map_indices;  ///< Indicies into tile_map_indices where each supertile's
                                      ///<   list of tile map indices can be found.  These indices
                                      ///<   are pre-inflated by 256 (the number of tiles in a
                                      ///<   supertile).
  Hybrid<int> tile_map_indices;       ///< Indices into the all_masks array where each tile's
                                      ///<   exclusions are to be found.  These indices are
                                      ///<   pre-inflated by 32 (see above).
  AtomGraph *ag_pointer;              ///< Pointer to the original topology
};

} // namespace energy
} // namespace stormm

#endif
