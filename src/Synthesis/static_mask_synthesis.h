// -*-c++-*-
#ifndef STORMM_STATIC_MASK_SYNTHESIS_H
#define STORMM_STATIC_MASK_SYNTHESIS_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Potential/static_exclusionmask.h"
#include "Topology/atomgraph.h"

namespace stormm {
namespace synthesis {

using card::Hybrid;
using card::HybridTargetLevel;
using energy::StaticExclusionMask;
using topology::AtomGraph;

/// \brief The read-only abstract for a static exclusion mask compilation.  This provides access
///        in a similar format to the static exclusion mask reader for a single system.
struct SeMaskSynthesisReader {

  /// \brief The constructor takes the object's constants and data pointers on either the CPU or
  ///        HPC resources.
  SeMaskSynthesisReader(int nsys_in, const int* atom_counts_in, const int* atom_offsets_in,
                        int nsupertile_in, int ntile_in, const int* supertile_map_idx_in,
                        const int* supertile_map_bounds_in, const int* tile_map_idx_in,
                        const uint* mask_data_in);

  /// \brief Like other abstracts, this one is compatible with the default copy and move
  ///        constructors, but the copy and move assignment operators are implicitly deleted.
  /// \{
  SeMaskSynthesisReader(const SeMaskSynthesisReader &original) = default;
  SeMaskSynthesisReader(SeMaskSynthesisReader &&original) = default;
  /// \}
  
  const int nsys;                   ///< The number of systems covered by this object
  const int* atom_counts;           ///< Counts of atoms in all systems
  const int* atom_offsets;          ///< Offsets for the atoms of any given system within a
                                    ///<   corresponding coordinate or topology synthesis
  const int nsupertile;             ///< Number of unique supertiles stored by the mask
  const int ntile;                  ///< Number of unique tiles stored by the mask
  const int* supertile_map_idx;     ///< Supertile maps for all systems, stored in a bounded list
                                    ///<   delineated by supertile_map_bounds
  const int* supertile_map_bounds;  ///< Bounds array for the supertile index maps of each system
                                    ///<   stored in supertile_map_idx
  const int* tile_map_idx;          ///< Tile maps for all unique supertiles
  const uint* mask_data;            ///< Mask data for all unique tiles
};

/// \brief An exclusion mask object for a compilation of systems.  All systems are represented in
///        full detail, and a series of systems using the same topology will each get their own
///        entries in the supertile map indices array, filtering on down to tile indices and the
///        unique mask data.
class StaticExclusionMaskSynthesis {
public:

  /// \brief The constructor requires a list of pre-computed StaticExclusionMask objects, each of
  ///        which will contain pointers to the original topologies, and a list of indices stating
  ///        which masks to apply in a given order for the synthesis of systems.
  ///
  /// Overloaded:
  ///   - Create an empty object
  ///   - Accept a list of static exclusion mask object pointers
  ///   - Accept a list of static exclusion mask objects
  ///   - Accept a list of topology pointers, from which exclusion masks will be constructed
  ///   - Accept an existing synthesis of topologies
  ///
  /// \param base_masks        List of exclusion masks for all unique topologies in the synthesis
  /// \param base_topologies   List of topologies from which to construct the exclusion masks
  ///                          (providing this removes the need to generate the list of exclusion
  ///                          masks outside this function)
  /// \param topology_indices  List of indices into base_masks indicating how to compile the
  ///                          synthesis of systems
  /// \{
  StaticExclusionMaskSynthesis();
  
  StaticExclusionMaskSynthesis(const std::vector<StaticExclusionMask*> &base_masks,
                               const std::vector<int> &topology_indices);

  StaticExclusionMaskSynthesis(const std::vector<StaticExclusionMask> &base_masks,
                               const std::vector<int> &topology_indices);

  StaticExclusionMaskSynthesis(const std::vector<AtomGraph*> &base_toplogies,
                               const std::vector<int> &topology_indices);
  /// \}

  /// \brief With no const members or pointers to repair (including POINTER-kind Hybrid) objects,
  ///        this object is compatible with the default copy and move constructors as well as copy
  ///        and move assignment operators.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right hand side of the assignment statement
  /// \{
  StaticExclusionMaskSynthesis(const StaticExclusionMaskSynthesis &original) = default;
  StaticExclusionMaskSynthesis(StaticExclusionMaskSynthesis &&original) = default;
  StaticExclusionMaskSynthesis& operator=(const StaticExclusionMaskSynthesis &original) = default;
  StaticExclusionMaskSynthesis& operator=(StaticExclusionMaskSynthesis &&original) = default;  
  /// \}
  
  /// \brief Get the number of systems covered by this object
  int getSystemCount() const;
  
  /// \brief Get the number of atoms in one of the systems described by the mask.
  ///
  /// \param index  The system index for which to query the atom count
  int getAtomCount(int index = 0) const;

  /// \brief Get the starting position of atoms in one of the systems described by the mask.
  ///
  /// \param index  The system index for which to query the atom offset
  int getAtomOffset(int index = 0) const;

  /// \brief Obtain from the mask whether a combination of two atoms in a particular system
  ///        constitutes an exclusion.
  ///
  /// \param system_index  Index the system of interest
  /// \param atom_i        Index of the first atom within the system of interest
  /// \param atom_j        Index of the second atom within the system of interest
  bool testExclusion(int system_index, int atom_i, int atom_j) const;

  /// \brief Get the abstract for this static exclusion mask synthesis.
  SeMaskSynthesisReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

#ifdef STORMM_USE_HPC
  /// \brief Upload the exclusion mask synthesis to the GPU device.
  void upload();

  /// \brief Download the exclusion mask synthesis from the GPU device.
  void download();
#endif

private:
  int system_count;                   ///< The number of systems covered by this object
  Hybrid<int> atom_counts;            ///< Atom counts for all systems.  Supertile strides are not
                                      ///<   stored, as this object will be used primarily to
                                      ///<   direct the creation of instructions for non-bonded
                                      ///<   work units.  The atom count is a more meaningful
                                      ///<   quantity, from which the number of supertile strides
                                      ///<   can be readily derived.
  Hybrid<int> atom_offsets;           ///< Atomic offsets for all systems.  These will be
                                      ///<   calculated with the same warp size padding as their
                                      ///<   cognates in the PhaseSpaceSynthesis and
                                      ///<   AtomGraphSynthesis.
  int unique_supertile_count;         ///< The number of unique supertiles, essentially a sum of
                                      ///<   unique supertiles in all unique masks but omitting the
                                      ///<   zero supertile for all but the first mask.
  int unique_tile_count;              ///< Number of unique tiles compiled from all unique masks
  Hybrid<int> supertile_map_indices;  ///< Indicies into tile_map_indices where each supertile's
                                      ///<   list of tile map indices can be found.  These indices
                                      ///<   are pre-inflated by 256 (the number of tiles in a
                                      ///<   supertile).
  Hybrid<int> supertile_map_bounds;   ///< Bounds array for the supertile_map_indices array above.
                                      ///<   This is the substantive change from single-topology
                                      ///<   static exclusion masks: the bounds array provides
                                      ///<   offsets for each system, avoiding a O(N^2) explosion
                                      ///<   in the supertile indexing, even though such an
                                      ///<   explosion would be somewhat muted by the large size
                                      ///<   of supertiles.  The presence of this array also
                                      ///<   negates the need for a hypothetical array of starting
                                      ///<   indices for each system's atoms in some concatenated
                                      ///<   list, which would not match up with the concatenation
                                      ///<   in a PhaseSpaceSynthesis anyway.
  Hybrid<int> tile_map_indices;       ///< Indices into the all_masks array where each tile's
                                      ///<   exclusions are to be found.  These indices are
                                      ///<   pre-inflated by 32 (see above).
  Hybrid<uint> all_masks;             ///< The actual mask data.  Each tile's exclusion map is
                                      ///<   a stretch of 32 unsigned integers in this array (16
                                      ///<   for the tile's sending atoms and 16 for the tile's
                                      ///<   receiving atoms).  

  /// \brief Encapsulate the contents of the constructor to permit multiple pathways for creating
  ///        the object.
  ///
  /// \param base_masks        List of exclusion masks for all unique topologies in the synthesis
  /// \param topology_indices  List of indices into base_masks indicating how to compile the
  ///                          synthesis of systems
  void build(const std::vector<StaticExclusionMask*> &base_masks,
             const std::vector<int> &topology_indices);
};
  
} // namespace synthesis
} // namespace stormm

#endif
