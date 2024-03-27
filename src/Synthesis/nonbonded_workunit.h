// -*-c++-*-
#ifndef STORMM_NONBONDED_WORKUNIT_H
#define STORMM_NONBONDED_WORKUNIT_H

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Math/summation.h"
#include "Topology/atomgraph.h"
#include "Potential/static_exclusionmask.h"
#include "static_mask_synthesis.h"
#include "synthesis_enumerators.h"

namespace stormm {
namespace synthesis {

using card::GpuDetails;
using energy::StaticExclusionMask;
using energy::tile_length;
using energy::tile_lengths_per_supertile;
using stmath::sum;
  
/// \brief The maximum number of imported atoms in a "small" nonbonded block (256 threads)
constexpr int small_block_max_atoms = 320;

/// \brief The maximum number of tile sides' worth of atoms that can be imported
constexpr int small_block_max_imports = small_block_max_atoms / tile_length;

/// \brief The number of tiles that a small block will be designed to process simultaneously
constexpr int small_block_tile_width = 8;

/// \brief Tile counts for various sizes of non-bonded work units that use the small block size
/// \{
constexpr int tiny_nbwu_tiles   =     small_block_tile_width;
constexpr int small_nbwu_tiles  = 2 * small_block_tile_width;
constexpr int medium_nbwu_tiles = 4 * small_block_tile_width;
constexpr int large_nbwu_tiles  = 8 * small_block_tile_width;
/// \}

/// \brief Tile count for non-bonded work units using the large block size
constexpr int huge_nbwu_tiles   = tile_lengths_per_supertile * tile_lengths_per_supertile;

/// \brief Length of the abstract for a non-bonded work unit based on all-to-all tile groups
constexpr int tile_groups_wu_abstract_length = 64;

/// \brief Length of the abstract for a non-bonded work unit based on all-to-all supertiles
constexpr int supertile_wu_abstract_length = 8;

/// \brief Collect a series of tiles for non-bonded computations as well as the required atom
///        imports to carry them out.  This will accomplish the task of planning the non-bonded
///        computation, given a single topology or a synthesis of topologies, to optimize the
///        thread utilization on a GPU.
class NonbondedWorkUnit {
public:

  /// \brief The constructor accepts an exclusion mask and a list of nonbonded interaction tiles
  ///        to compute.  Non-bonded interaction tiles go according to "abscissa atoms" and
  ///        "ordinate atoms," although both abscissa and ordinate indices refer to the same list
  ///        of imported atoms.  There are three cases to consider:
  ///
  ///        Static exclusion mask, tiny up to large work units:
  ///        This will be the majority of the cases with implicit-solvent systems.  In practice,
  ///        the non-bonded work unit becomes a list of 50 integers.  A complete description is
  ///        available in the documentation for the getAbstract() member function of this object.
  ///        These cases will run with a 256-thread block size.
  ///
  ///        Static exclusion mask, huge work units:
  ///        This will handle cases of implicit-solvent systems with sizes so large that millions
  ///        of smaller work units would be required to cover everything.  In practice, the
  ///        non-bonded work unit is reduced to a list of 8 integers, now representing the lower
  ///        limits of the abscissa and ordinate atoms to import, and the numbers of atoms to
  ///        import along each axis, in a supertile for which the work unit is to compute all
  ///        interactions.  There is no meaningful list of all interactions in this case, as it
  ///        might be prohibitive even to store such a thing.  Instead, the work unit will proceed
  ///        over all tiles in the supertile after computing whether it lies along the diagonal.
  ///        This will require a larger block size (512 threads minimum, up to 768 depending on
  ///        the architecture).
  ///
  ///        Forward exclusion mask:
  ///        This will handle cases of neighbor list-based non-bonded work units.  The work unit
  ///        assumes a honeycomb-packed image of all atoms in or about the primary unit cell (some
  ///        honeycomb pencils will straddle the unit cell boundary but their positions will be
  ///        known as part of the decomposition).  The work unit will consist of thirty integers:
  ///        seven starting locations of atom imports, seven bit-packed integers detailing the
  ///        lengths of each stretch of atoms (first 20 bits) and the obligatory y- and z- imaging
  ///        moves to make with such atoms (last 12 bits), seven integers detailing segments of
  ///        each stretch of imported atoms to replicate in +x (and where in the list of imported
  ///        atoms to put them), and finally seven integers detailing segments of each stretch of
  ///        imported atoms to replicate in -x (and where to put the replicas).  The final two
  ///        integers state the starting and ending indices of a list of tile instructions to
  ///        process.  The tile instructions for the neighbor list-based non-bonded work units are
  ///        much more complex than those for non-bonded work based on a static exclusion mask.
  ///
  /// \param se              Static exclusion mask, or synthesis thereof, for one or more systems
  ///                        in isolated boundary conditions
  /// \param tile_list       Paired with a static exclusion mask and non-huge tiles, the specific
  ///                        list of tiles to include for this work unit in the x- and y-members,
  ///                        plus the system index number in the z member (if more than one system
  ///                        is present in a synthesis).  Among them, the tiles must not call for
  ///                        importing more than small_block_max_atoms atoms.
  /// \param abscissa_start  The abscissa axis start of the supertile to process.  Paired with a
  ///                        static exclusion mask in extreme circumstances of very large implicit
  ///                        solvent systems.  This will call for importing 512 atoms
  ///                        (2 x supertile_length) in the most general case and will require
  ///                        larger thread blocks.  If computing for a synthesis of static
  ///                        exclusion masks, the abscissa starting point is given as a relative
  ///                        index within the local system to which the supertile belongs.
  /// \param ordinate_start  The ordinate axis start of the supertile to process.  If computing for
  ///                        a synthesis of static exclusion masks, the ordinate starting point is
  ///                        given as a relative index within the local system to which the
  ///                        supertile belongs.
  /// \{
  NonbondedWorkUnit(const StaticExclusionMask &se, const std::vector<int3> &tile_list);
  
  NonbondedWorkUnit(const StaticExclusionMask &se, int abscissa_start, int ordinate_start);

  NonbondedWorkUnit(const StaticExclusionMaskSynthesis &se,
                    const std::vector<int3> &tile_list);

  NonbondedWorkUnit(const StaticExclusionMaskSynthesis &se, int abscissa_start, int ordinate_start,
                    int system_index);
  /// \}

  /// \brief Take the default copy and move constructors as well as assignment operators for this
  ///        object, which uses only Standard Template Library member variable types.
  /// \{
  NonbondedWorkUnit(const NonbondedWorkUnit &original) = default;
  NonbondedWorkUnit(NonbondedWorkUnit &&original) = default;
  NonbondedWorkUnit& operator=(const NonbondedWorkUnit &other) = default;
  NonbondedWorkUnit& operator=(NonbondedWorkUnit &&other) = default;
  /// \}
  
  /// \brief Get the tile count of this work units.
  int getTileCount() const;

  /// \brief Get the number of tile_length atom imports needed for this work units.
  int getImportCount() const;

  /// \brief Get the mask for initializing per-atom properties on any of the imported atom groups.
  int getInitializationMask() const;

  /// \brief Get the abscissa and ordinate atom limits for a tile from within this work unit.
  ///        The abscissa limits are returned in the x and y members, the ordinate limits in the
  ///        z and w members.
  ///
  /// \param index  Index of the tile from within the work unit
  int4 getTileLimits(int index) const;

  /// \brief Get the list of tile instructions for this work unit.
  const std::vector<uint2>& getTileInstructions() const;
  
  /// \brief Get an abstract for this work unit, layed out to work within an AtomGraphSynthesis.
  ///        For TILE_GROUPS-type work units serving systems in isolated boundary conditions, the
  ///        non-bonded abstract has the following format:
  ///
  ///  Slot   Description
  /// ------  -------------------------------------------------------------------------------------
  ///      0  The total number of atom group imports for various tiles.  Each group is up to 16
  ///         atoms, and all atoms pertain to a single system in the synthesis.
  ///   1-20  Starting atom indices of each tile group in the topology / coordinate synthesis
  ///         (the extent of this segment is set by small_block_max_imports)
  ///  21-25  Counts of the numbers of atoms in each tile group, with four groups' counts packed
  ///         into each int.  In this scheme, the tile length may be extended up to 128 atoms (256
  ///         if these values are read as unsigned ints), but a tile length of 16 appears to be
  ///         optimal for NVIDIA and probably other architectures as well.
  ///  26-27  Starting and ending locations of the tile instruction list.  After reading up to 20
  ///         tile atom groups, each work unit tries to combine the groups into different
  ///         combinations of actual tiles.  The two interacting groups are given in the x-member
  ///         of a uint2 tuple, while the starting location of exclusion bit masks for the tile
  ///         is given in the y-member.  The difference between values in slots 26 and 27 indicates
  ///         the number of tiles to perform, and is designed to be a multiple of eight if possible
  ///         (anticipating blocks of four or eight warps, with one warp handling each tile).
  ///  28-47  System indices within the synthesis to which the atoms in each tile group pertain
  ///     48  Bitmask indicating whether the current work unit is the first to import a particular
  ///         group of atoms and therefore can be tasked with contributing certain per-atom effects
  ///         when accumulating forces or radii due to those atoms
  ///  49-50  Lower and upper limits of atoms for which the work unit is tasked with initializing
  ///         force or other property accumulators to have them ready for use in future iterations
  ///         of the molecular mechanics force / energy evaluation cycle.  The integer in slot
  ///         49 indicates the first atom of a contiguous list, and need to be an atom that the
  ///         work unit imported for one of its tile computations.  The integer in slot 50 is a
  ///         bit-packed value, with the low 16 bits indicating up to 65,504 (not 65,536) atoms
  ///         following the first index.  The high 16 bits of the integer in slot 50 indicate
  ///         whether to initialize force X, Y, or Z as well as psi or sum_deijda accumulators for
  ///         the next iteration of the cycle (one bit per accumulator), as well as whether to
  ///         perform random number caching for up to 15 cycles.  These two slots of the work unit
  ///         can be altered after the non-bonded work unit list is constructed.
  ///
  /// \param instruction_start  The starting point of instructions for this group of tiles, if
  ///                           the work unit will fit on a small thread block, having less than
  ///                           huge_nbwu_tiles tiles.  Otherwise no instruction start is needed.
  std::vector<int> getAbstract(int instruction_start = 0) const;

  /// \brief Set the initialization mask for atomic properties that contribute once to an
  ///        accumulated sum, i.e. contributions of baseline atomic radii to the Generalized Born
  ///        effective radii.
  ///
  /// \param mask_in  The new mask to use.  If the jth bit of this mask is set to 1, the work
  ///                 unit will perform initialization protocols on its jth set of atom imports,
  ///                 i.e. the Born radii derivatives for those atoms will be contributed to the
  ///                 force accumulators.
  void setInitializationMask(int mask_in);

  /// \brief Set the atom index at which to begin accumulator refreshing.  The actual work done
  ///        depends on the accumulator refresh code set by the next member function.
  ///
  /// \param index_in  The starting atom index relevant to this work unit
  void setRefreshAtomIndex(int index_in);
  
  /// \brief Set the accumulator refreshing instructions for this work unit, i.e. "set X force
  ///        accumulators to zero for this many atoms." The refreshing starts at the atom index
  ///        set by the preceding member function.
  ///
  /// \param code_in  The refesh code to execute
  void setRefreshWorkCode(int code_in);
  
private:
  int tile_count;                          ///< Number of tiles to be processed by this work unit
  int import_count;                        ///< Number of imported tile abscissa or ordinate atom
                                           ///<   sets to be stored locally by this work unit
  NbwuKind kind;                           ///< The type of non-bonded work unit.  For isolated
                                           ///<   boundary conditions, there is a choice between
                                           ///<   TILE_GROUPS and SUPERTILES.
  int score;                               ///< The estimated effort score of this work unit
  int init_accumulation;                   ///< Some calculations such as Generalized Born solvent
                                           ///<   effects will require that certain properties be
                                           ///<   initialized on a per-atom basis.  This mask, its
                                           ///<   jth bit being set to 1 to indicate that the
                                           ///<   property of interest should be computed for
                                           ///<   imported group j, manages the problem of
                                           ///<   initializing these properties and counting their
                                           ///<   contributions once and only once.
  int refresh_atom_start;                  ///< First atom at which to begin resetting accumulation
                                           ///<   counters in some alternate force or intermediate
                                           ///<   quantity array.  The non-bonded kernels will
                                           ///<   spread out the work of initializing accumulators
                                           ///<   for subsequent force and energy calculations.
  int refresh_code;                        ///< Instructions for this work unit on what to reset
                                           ///<   (bits 0-7), how many random number sets to
                                           ///<   compute (bits 8-15) and how many atoms (starting
                                           ///    from refresh_atom_start) to apply it to
  std::vector<int> imports;                ///< List of starting positions for atoms that must be
                                           ///<   cached in order to process all tiles in this
                                           ///<   work unit
  std::vector<int> import_system_indices;  ///< System indices of each imported tile.  All atoms
                                           ///<   in any one tile pertain to the same system, but
                                           ///<   each imported tile within the work unit may be
                                           ///<   part of a different system.
  std::vector<int> import_size_keys;       ///< Bit-packed integers with the number of atoms to
                                           ///<   import after each starting position
  std::vector<uint2> tile_instructions;    ///< Instructions for processing each tile based on an
                                           ///<   appropriate exclusion mask

  /// \brief Find the system to which a particular import belongs.  This function is only used in
  ///        the event that the non-bonded work unit pertains to synthesis of systems, otherwise
  ///        the system index is obviously zero.
  ///
  /// \param ser         Reader abstract for the exclusion mask synthesis
  /// \param atom_index  Index of the first imported atom in the tile of interest
  int getImportSystemIndex(const SeMaskSynthesisReader &ser, int atom_index);
};

/// \brief Add one tile to the growing list that will eventually define a work unit, if it is
///        possible to do so.  This encapsualtes some work at a rather high cost in terms of
///        pointers for multiple returned values.  Return true if the tile addition was successful
///        or false if not.
///
/// \param tile_list           The growing list of tiles, pre-allocated to be able to hold any
///                            additions that this function might make.  Appended and returned.
/// \param import_coverage     Array of 0's and 1's (an int, rather than boolean, array for
///                            versatility in passing to and from this function and for direct
///                            summation with other integers).  Spans all systems but does not
///                            necessarily align with the atom offsets in the corresponding phase
///                            space synthesis, topology synthesis, or mask synthesis.  Edited and
///                            returned.
/// \param system_tile_starts  Bounds array for systems within the import coverage array
/// \param import_count        The total number of imported tiles' worth of atoms.  Updated and
///                            returned.
/// \param current_tile_count  The current number of tiles, updarted and returned.
/// \param ti                  Tile index of the abscissa atoms (starting atom index divided by
///                            tile_length)
/// \param tj                  Tile index of the ordinate atoms (starting atom index divided by
///                            tile_length)
/// \param sysid               System index from which the tiles are derived
bool addTileToWorkUnitList(int3* tile_list, int* import_coverage,
                           const std::vector<int> &system_tile_starts, int *import_count,
                           int *current_tile_count, const int ti, const int tj, const int sysid);

/// \brief Set the masks directing each non-bonded work unit with exclusive rights to perform
///        certain initializations that must only be performed once per atom.
///
/// \param all_wu  A list of non-bonded work units, presumably with blank initialization masks
/// \param kind    The type of work units in the list
void setPropertyInitializationMasks(std::vector<NonbondedWorkUnit> *all_wu, NbwuKind kind);

/// \brief Obtain a code for initializing accumulators for use on subsequent cycles.
///
/// \param init_request
/// \param cache_depth   The number of random numbers for influencing atomic motions in the X, Y,
///                      and Z directions 
int nonbondedWorkUnitInitCode(const InitializationTask init_request, const int cache_depth);

/// \brief Distribute the total atom count over each non-bonded work unit to balance the effort of
///        initializing accumulators and computing random numbers for the entire system or
///        synthesis of systems.
///
/// \param result     The list of non-bonded work units, modified and returned
/// \param natom      Total number of atoms to distribute (the atoms that each work unit will
///                   initialize are not necessarily the atoms for which it computes interactions)
/// \param init_code  Only the first 16 bits of this int are relevant.  This code indicates whether
///                   the non-bonded work units are to initialize force accumulators, Generalized
///                   Born intermediate accumulators, or random values.
/// \param gpu        GPU specifications (HPC compilation only, feeding this an argument of
///                   null_gpu will not result in the program taking special considerations for the
///                   size of the workload)
void distributeInitializationRanges(std::vector<NonbondedWorkUnit> *result, int natom,
                                    int init_code, const GpuDetails &gpu = null_gpu);
  
/// \brief Create a list of work units of a consistent size so as to optimize the load balancing
///        on a given resource (GPU, on or more CPUs, etc.).
///
/// Overloaded:
///   - Take a single topology or a list of topologies
///   - Take one or more static exclusion masks corresponding to the input topologies
///   - Accept a target number of work units
///   - Accept GPU specifications from which to infer a target number of work units
///
/// \param se                   Static exclusion mask corresponding to the input topology
/// \param poly_se              Synthesis of static exclusion masks corresponding to each input
///                             topology
/// \param init_request         Indicate a pattern for the non-bonded work units to initialize
///                             accumulators for subsequent force and energy calculations.
/// \param random_cache_depth   Number of random values to store for each atom in the synthesis
/// \param target_nbwu_count    The target number of work units to produce.  The goal will be to
///                             get as close to this number without going over, or if the number
///                             must be exceeded due to absolute limits on the size of any one
///                             work unit, to produce a number of work units that comes as close
///                             to a multiple of this number as possible, without going over.
/// \param gpu                  GPU specifications (HPC compilation only)
/// \{
std::vector<NonbondedWorkUnit>
buildNonbondedWorkUnits(const StaticExclusionMaskSynthesis &poly_se,
                        InitializationTask init_request = InitializationTask::NONE,
                        int random_cache_depth = 0, const GpuDetails &gpu = null_gpu);

std::vector<NonbondedWorkUnit>
buildNonbondedWorkUnits(const StaticExclusionMask &se,
                        InitializationTask init_request = InitializationTask::NONE,
                        int random_cache_depth = 0);
/// \}

} // namespace synthesis
} // namespace stormm

#include "nonbonded_workunit.tpp"

#endif
