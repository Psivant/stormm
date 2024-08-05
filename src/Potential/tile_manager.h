// -*-c++-*-
#ifndef STORMM_TILE_MANAGER_H
#define STORMM_TILE_MANAGER_H

#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/stormm_vector_types.h"

namespace stormm {
namespace energy {

using card::Hybrid;
using card::HybridTargetLevel;
  
/// \brief The abstract of the TileManager contains pointers to direct the relative indexing of
///        atom reads and the rearrangements that follow the tile evaluation needed to bring the
///        the atoms back into an arrangement suitable for reducing the forces.  The abstract is
///        writeable.
struct TilePlan {

  /// \brief The constructor follows from other abstracts and takes a list of all pointers and
  ///        critical constants.
  TilePlan(int nchoice_in, const int* read_assign_in, const int* self_assign_in,
           const int* reduce_prep_in, const int* self_prep_in, const float* scalings_in,
           const int* nt_stencil_in, int* xfrc_ovrf_in, int* yfrc_ovrf_in, int* zfrc_ovrf_in);

  /// \brief The presence of const member variables negates copy and move assignment operations,
  ///        but the default copy and move constructors are valid.
  ///
  /// \param original  The original TilePlan to copy or move
  /// \{
  TilePlan(const TilePlan &original) = default;
  TilePlan(TilePlan &&original) = default;
  /// \}

  const int nchoice;       ///< The number of choices for batch degeneracy in either sending or
                           ///<   receiving atoms.  The read_assign and reduce_prep arrays are
                           ///<   arranged such that the relative indices for receiving atoms in a
                           ///<   { 1, 1 } degeneracy (each batch is completely full, no atoms are
                           ///<   replicated) are followed by indices for a { 1, 2 } degeneracy
                           ///<   (the batch of receiving atoms is half full, each atom appears
                           ///<   twice) and then a { 1, 4 } degeneracy (each of the receiving
                           ///<   atoms appears four times in the overall batch).  If there are
                           ///<   higher levels of degeneracy, they would appear next, but the
                           ///<   degeneracy of the sending atoms then ticks forward and reading /
                           ///<   rearrangement instructions for { 2, 1 } degeneracy then appear,
                           ///<   followed by { 2, 2 } and { 2, 4 }.  Eventually, { 4, 1 } to
                           ///<   { 4, 4 } and perhaps higher degeneracies are covered.  Both
                           ///<   arrays are thus sized according to the warp size and the square
                           ///<   of nchoice.
  const int* read_assign;  ///< Reading assignments for threads, following the layout described
                           ///<   above.
  const int* self_assign;  ///< Reading assignments for threads in self-interaction tiles, arranged
                           ///<   in a manner similar to read_assign but with only M x M tiles
                           ///<   (squares, not rectangles).
  const int* reduce_prep;  ///< Reduction preparation shuffle lane assignments for threads.  There
                           ///<   are nchoice sets of assignments, each the width of a warp.
  const int* self_prep;    ///< Reduction preparation shuffle lane assignments for threads handling
                           ///<   self interaction tiles.  There are two time nchoice sets of data,
                           ///<   each the width of a warp.  The first nchoice sets handle tiles
                           ///<   where the sending and receiving atoms are truly one and the same
                           ///<   and the number of iterations in the inner loop is thereby reduced
                           ///<   by half.  The second nchoice sets handle tiles where the sending
                           ///<   and receiving atoms are not the same but nonetheless come from
                           ///<   the same central neighbor list decomposition cell.
  const float* scalings;   ///< Scaling factors to be applied to each interaction computed by the
                           ///<   threads of a warp in the first iteration of evaluating the tile
  const int* nt_stencil;   ///< Coded bitmasks instructing each thread which cell, relative to the
                           ///<   center, to access in order to complete the neutral territory
                           ///<   tower-plate arrangement.  These instructions will be combined
                           ///<   with knowledge of the system itself in order to arrive at the
                           ///<   exact cell index from which to draw atoms.
  int* xfrc_ovrf;          ///< Overflow accumulators for forces on receiving atoms in the
                           ///<   Cartesian X direction
  int* yfrc_ovrf;          ///< Overflow accumulators for receiving atom forces in the Y direction
  int* zfrc_ovrf;          ///< Overflow accumulators for receiving atom forces in the Z direction
};

/// \brief When tiles are loaded for what could be partial batches of atoms, it is critical for
///        threads of the warp to be able to quickly calculate which atoms to load and perhaps
///        replicate.  After the tile has been evaluated, a rearrangement may be necessary to
///        put atoms back in an order that the accumulated forces can be reduced and then added
///        back to global arrays.  This class will store tables of relative reading assignments
///        and the preparatory rearrangements needed for the reduction step.  The tables will be
///        read with full L1-caching protocols to occupy a few kB of L1 with frequent re-use during
///        particle-particle interaction tile evaluation.
class TileManager {
public:

  /// \brief The constructor requires specifications of the GPU and can accept a value for the
  ///        maximum degeneracy for testing purposes.
  ///
  /// \param launch_parameters  Specifications of for allocating block-specific storage arrays 
  /// \param max_deg_in         Maximum degeneracy to permit in any one atom batch fulfilling one
  ///                           side of the tile.  The default of zero or less engages the default
  ///                           degeneracies.
  TileManager(const int2 launch_parameters, int max_deg_in = 0);

  /// \brief The default copy and move constructors as well as assignment operators are all valid.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object to fulfill the right hand side of the assignment operation
  /// \{
  TileManager(const TileManager &original);
  TileManager(TileManager &&original);
  TileManager& operator=(const TileManager &original);
  TileManager& operator=(TileManager &&original);
  /// \}
  
  /// \brief Get the base-two logarithm of the maximum degeneracy in the batch atoms.
  int getMaximumBatchDegeneracy() const;

  /// \brief Get the layout of sending atoms for a given degeneracy.  This is for convenience in
  ///        testing or prototyping new methods with the tile scheme.
  ///
  /// \param sending_atom_degeneracy  Base-2 logarithm of degeneracy in the sending atom batch
  std::vector<int> getSendingAtomLayout(const int sending_atom_degeneracy) const;

  /// \brief Get the list of reading assignments (for inspection) pertaining to given degeneracies
  ///        in the sending (tile X-axis) and receiving (tile Y-axis) atoms.
  ///
  /// \param sending_atom_degeneracy  Base-2 logarithm of degeneracy in the sending atom batch
  /// \param recving_atom_degeneracy  Base-2 logarithm of degeneracy in the receiving atom batch
  std::vector<int> getReadAssignments(int sending_atom_degeneracy,
                                      int recving_atom_degeneracy) const;

  /// \brief Get the list of reading assignments (for inspection) pertaining to given degeneracies
  ///        in the tiles of the central cell, where sending and receiving atoms may be one and
  ///        the same.  Descriptions of input parameters follow from getReadAssignments(), above.
  std::vector<int> getSelfAssignments(int sending_atom_degeneracy,
                                      int recving_atom_degeneracy) const;
  
  /// \brief Get the list of reduction preparations (for inspection) pertaining to given
  ///        degeneracies in the sending (tile X-axis) and receiving (tile Y-axis) atoms.
  ///        Descriptions of input parameters follow from getReadAssignments(), above.
  std::vector<int> getReductionPreparations(int sending_atom_degeneracy,
                                            int recving_atom_degeneracy) const;

  /// \brief Get the list of reduction preparations (for inspection) as above, but for tiles that
  ///        might include self interactions.  Descriptions of input parameters follow from
  ///        getReadAssignments(), above.
  std::vector<int> getSelfPreparations(int sending_atom_degeneracy,
                                       int recving_atom_degeneracy) const;

  /// \brief Get the scaling factors to be applied to each thread's computed interaction in the
  ///        first iteration to evaluate the tile, given an atom degeneracy.  This is for developer
  ///        inspection and debugging.
  std::vector<float> getThreadScalings(int atom_degeneracy) const;
  
  /// \brief Get the map of the neutral territory decomposition, for inspection.
  std::vector<int> getNeutralTerritoryStencil() const;
  
  /// \brief Get the abstract.
  ///
  /// \param tier  Indicate whether to obtain pointers to data on the CPU host or GPU device
  TilePlan data(HybridTargetLevel tier = HybridTargetLevel::HOST);

#ifdef STORMM_USE_HPC
  /// \brief Upload tile plans to the GPU.
  void upload();

  /// \brief Download tile plans from the GPU.
  void download();
#endif

private:

  int maximum_degeneracy;           ///< Log of the maximum degeneracy in atom batches from either
                                    ///<   the sending or receiving atoms
  int block_count;                  ///< The number of thread blocks for which the scratch tables
                                    ///<   of force accumulators are prepared
  int thread_count;                 ///< The number of threads per block for which the scratch
                                    ///<   tables of force accumulators are prepared
  Hybrid<int> read_assignments;     ///< Assignments for reading atoms for the receiving batch.
                                    ///<   The sending batch assignments are simple to calculate:
                                    ///<   the relative indices simply proceed in increasing order
                                    ///<   up to the total number of atoms in the batch, then the
                                    ///<   sequence repeats until the warp is filled.
  Hybrid<int> self_assignments;     ///< Assignments for reading atoms of the receiving batch in
                                    ///<   tiles evaluated for the central cell.  These are
                                    ///<   specialized to account for the case of tiles wherein
                                    ///<   the sending and receiving atoms are the same.
  Hybrid<int> reduce_preparations;  ///< The other threads that each thread of the warp must read
                                    ///<   from (shuffle instructions) so that the degenerate batch
                                    ///<   of atoms into an order suitable for reducing accumulated
                                    ///<   forces.
  Hybrid<int> self_preparations;    ///< These assignments serve the same purpose as those in
                                    ///<   reduce_preparations, but applicable to tiles in the 
                                    ///<   central cell where sending and receiving atoms might be
                                    ///<   one and the same.  There are two halves to this array:
                                    ///<   the first half contains sets of shuffle assignments
                                    ///<   to be used when the sending and receiving atoms are the
                                    ///<   same.  The second half contains sets of atoms to be used
                                    ///<   when the sending and receiving atoms are different (this
                                    ///<   will be the same content as the trace of the
                                    ///<   reduce_preparations matrix of assignments).  Tiles of
                                    ///<   central cell are all square, and when the sending and
                                    ///<   receiving atoms are the same the number of iterations
                                    ///<   needed to evaluate the tile is cut by half, implying a
                                    ///<   different set of shuffling assignments to prepare for
                                    ///<   data dumps.
  Hybrid<float> thread_scalings;    ///< Values to be applied to interactions (force as well as
                                    ///<   energy) computed by each thread in the first iteration
                                    ///<   of tile evaluation, if requested in a self-interacting
                                    ///<   tile where sending and receiving atoms are one and the
                                    ///<   same.
  Hybrid<int> tower_plate_stencil;  ///< Stencil for relative displacements that threads of a warp
                                    ///<   should calculate in order to obtain their assigned cells
                                    ///<   when constructing indexing limits and prefix sums for
                                    ///<   the tower-and-plate neutral territory decomposition.
                                    ///<   Each element is a bit-pack integer encoding, in its low
                                    ///<   eight bits, the A-axis displacement, in its next eight
                                    ///<   bits the B-axis displacement, and in its third set of
                                    ///<   eight bits the C-axis displacement.  This will save a
                                    ///<   great deal of memory, as the work units for neighbor
                                    ///<   list cells could otherwise take significant resources.
  Hybrid<int> x_force_overflow;     ///< Overflow accumulators for Cartesian X forces, allocated
                                    ///<   for each thread of each block on the GPU
  Hybrid<int> y_force_overflow;     ///< Overflow accumulators for Cartesian Y forces
  Hybrid<int> z_force_overflow;     ///< Overflow accumulators for Cartesian Z forces
  Hybrid<int> int_data;             ///< Integer data storage array, targeted by the POINTER-kind
                                    ///<   integer Hybrid objects above

  /// \brief Allocate space for the object.
  void allocate();

  /// \brief Lay out the tower-plate neutral territory stencil.
  void planTowerPlateStencil();
  
  /// \brief Compute the order of unique receiving atom indices (relative indices for the warp to
  ///        read given the tile offset in whatever neutral-territory scheme) given the total
  ///        number of unique atoms, the number of repeats in the sending atom sequence, and the
  ///        stagger iteration.
  ///
  /// \param natom         The number of unique receiving atoms in the sequence
  /// \param sending_reps  The number of times that the batch of unique sending atoms repeats
  /// \param iter          Iteration count for laying sequences of repeating atoms across the warp.
  ///                      This can go as high as the warp size over natom.
  std::vector<int> computeStaggeredOrder(int natom, int sending_reps, int iter);
};
  
} // namespace energy
} // namespace stormm

#endif
