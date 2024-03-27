// -*-c++-*-
#ifndef STORMM_CELLGRID_H
#define STORMM_CELLGRID_H

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "Math/formulas.h"
#include "Math/math_enumerators.h"
#include "Math/matrix_ops.h"
#include "Math/rounding.h"
#include "Math/series_ops.h"
#include "Math/vector_ops.h"
#include "Numerics/split_fixed_precision.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/parse.h"
#include "Potential/energy_enumerators.h"
#include "Structure/local_arrangement.h"
#include "Structure/structure_enumerators.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_analysis.h"
#include "Topology/atomgraph_enumerators.h"
#include "Topology/topology_util.h"
#include "Trajectory/coordinate_copy.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/trajectory_enumerators.h"

namespace stormm {
namespace energy {

using card::GpuDetails;
using card::Hybrid;
using card::HybridKind;
using card::HybridTargetLevel;
using constants::CartesianDimension;
using constants::ExceptionResponse;
using constants::PrecisionModel;
using constants::UnitCellAxis;
using data_types::getStormmScalarTypeName;
using data_types::getStormmHpcVectorTypeName;
using data_types::isFloatingPointScalarType;
using numerics::force_scale_nonoverflow_bits;
using numerics::globalpos_scale_nonoverflow_bits;
using parse::NumberFormat;
using parse::realToString;
using stmath::hessianNormalWidths;
using stmath::indexingArray;
using stmath::ipowl;
using stmath::LimitApproach;
using stmath::matrixMultiply;
using stmath::maxValue;
using stmath::mean;
using stmath::nearestFactor;
using stmath::prefixSumInPlace;
using stmath::PrefixSumType;
using stmath::roundUp;
using structure::imageCoordinates;
using structure::ImagingMethod;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using synthesis::PsSynthesisWriter;
using synthesis::SyNonbondedKit;
using topology::AtomGraph;
using topology::hasVdwProperties;
using topology::inferCombiningRule;
using topology::NonbondedKit;
using topology::UnitCellType;
using trajectory::coordCopy;
using trajectory::CoordinateCycle;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameWriter;
using trajectory::TrajectoryKind;
  
/// \brief Default size constants for the cell grid
/// \{
constexpr size_t default_cellgrid_base_capacity = 64;
constexpr size_t cellgrid_tp_excl_allotment = 12;
constexpr size_t cellgrid_total_excl_allotment = 15;
constexpr int default_mesh_ticks = 4;
constexpr int maximum_spatial_decomposition_cells = 4095;
constexpr int maximum_cellgrid_cell_count = 268435456;
constexpr int maximum_xl_cell_count = 64;
constexpr double minimum_cell_width = 2.5;
constexpr double maximum_cell_width = 12.0;
constexpr int maximum_cellgrid_systems = 65536;
static const int sp_charge_index_bits = 23;
static const int dp_charge_index_bits = 32;
static const int sp_charge_index_mask = 0x7ffff;
static const llint dp_charge_index_mask = 0xffffffffLL;
/// \}

/// \brief Writeable abstract for the CellGrid object, able to receive new coordinates or
///        accumulate forces.
template <typename T, typename Tacc, typename Tcalc, typename T4> struct CellGridWriter {

  /// \brief As is other abstracts, the constructor takes a list of inputs for all member
  ///        variables.
  CellGridWriter(NonbondedTheme theme_in, int system_count_in, int total_cell_count_in,
                 int total_chain_count_in, int mesh_ticks_in, size_t cell_base_capacity_in,
                 size_t cell_excl_capacity_in, float lpos_scale_in, float lpos_inv_scale_in,
                 ullint* system_cell_grids_in, Tcalc* system_cell_umat_in, T* system_cell_invu_in,
                 T* system_pmig_invu, uint2* cell_limits_in, const uint2* cell_limits_old_in,
                 const uint* chain_limits_in, const int* system_chain_bounds_in,
                 const int* chain_system_owner_in, T4* image_in, const T4* image_old_in,
                 uint2* entry_room_in, uint2* exit_room_in, const int* transit_room_bounds_in,
                 int* entry_room_counts_in, int* exit_room_counts_in, int* nonimg_atom_idx_in,
                 uint* img_atom_idx_in, uint* exclusion_maps_in, const uint* exclusion_maps_old_in,
                 const int* nt_groups_in, Tacc* xfrc_in, Tacc* yfrc_in, Tacc* zfrc_in,
                 int* xfrc_ovrf_in, int* yfrc_ovrf_in, int* zfrc_ovrf_in, Tacc* xfrc_hw_in,
                 Tacc* yfrc_hw_in, Tacc* zfrc_hw_in, int* xfrc_hw_ovrf_in, int* yfrc_hw_ovrf_in,
                 int* zfrc_hw_ovrf_in);

  /// \brief The presence of const array sizing members implicitly deletes the copy and move
  ///        assignment operators, but the default copy and move constructors are valid.
  /// \{
  CellGridWriter(const CellGridWriter &original) = default;
  CellGridWriter(CellGridWriter &&original) = default;
  /// \}
  
  const NonbondedTheme theme;       ///< The non-bonded property for which the neighbor list is
                                    ///<   maintained (electrostatic, van-der Waals, or both)
  const int system_count;           ///< Number of systems whose spatial decompositions are stored
  const int total_cell_count;       ///< The total number of decomposition cells across all systems
  const int total_chain_count;      ///< The number of cell chains, groups of cells found in a row
                                    ///<   along the unit cell A axis of one simulation, across all
                                    ///<   systems
  const int mesh_ticks;             ///< The number of mesh grid spacings per spatial decomposition
                                    ///<   cell, pertaining to a particle-mesh interaction grid
                                    ///<   associated with the neighbor list decomposition.
  const size_t cell_base_capacity;  ///< The maximum capacity of any one cell (all cells in a row
                                    ///<   along a simulation box's A axis share their excess
                                    ///<   capacity, making a provision for one cell to take on
                                    ///<   a much larger number of atoms)
  const size_t cell_excl_capacity;  ///< Capacity allotted for exclusions associated with the
                                    ///<   cell's neutral territory interaction assignments
  const float lpos_scale;           ///< The local position scaling factor applied to coordinates
                                    ///<   in the images (expressed as a 32-bit floating point
                                    ///<   number, as this will be converted to 64-bit format with
                                    ///<   no more or less information if needed)
  const float lpos_inv_scale;       ///< The inverse scaling factor applied to coordinates in the
                                    ///<   images
  ullint* system_cell_grids;        ///< Dimensions of each cell grid and the starting point for
                                    ///<   its cells in the array.  This array has one element for
                                    ///<   each system in the underlying coordinate synthesis.
  Tcalc* system_cell_umat;          ///< Transformation matrices for taking Cartesian coordinates
                                    ///<   into the fractional space of the local coordinate system
                                    ///<   in each system's individual cells.
  T* system_cell_invu;              ///< Inverse transformation matrix for each system's individual
                                    ///<   cells.  The columns of this matrix indicate the size of
                                    ///<   the cell.
  T* system_pmig_invu;              ///< Inverse transformation matrix for each system's individual
                                    ///<   particle-mesh interaction grid elements.  The columns of
                                    ///<   this matrix indicate the size of the cell.
  uint2* cell_limits;               ///< Limits of the atoms in each cell out of the entire current
                                    ///<   image.  The "x" and member of the tuple provides
                                    ///<   an absolute index within the current image for the lower
                                    ///<   bound of atoms.  The "y" member provides, in its low
                                    ///<   16 bits, the system index (the maximum number of systems
                                    ///<   is set in maximum_cellgrid_systems).  The "y" member of
                                    ///<   the tuple also provides, in its high 16 bits, the number
                                    ///<   of atoms in the particular cell.  The system index must
                                    ///<   be used to retrieve information in system_cell_grids in
                                    ///<   order to put the cell in context and see, for example,
                                    ///<   which other cells are its neighbors.
  const uint2* cell_limits_old;     ///< Limits of atoms in each cell of the prior image.  The
                                    ///<   cell_limits member variable above describes the tuple.
  const uint* chain_limits;         ///< Immutable limits on the cell chains within each cell grid.
                                    ///<   The ith chain cannot hold more than a number of atoms
                                    ///<   given by chain_limits[i + 1] - chain_limits[i].  The
                                    ///<   atom capacity across chains in a particular system is
                                    ///<   consistent, but different systems may have different
                                    ///<   capacities.
  const int* system_chain_bounds;   ///< Bounds on the list of chains held by each system.  The
                                    ///<   process is for one or more blocks to look up these
                                    ///<   limits and then do chains from within the system before
                                    ///<   looking up new limits and moving on to the next block.
                                    ///<   With the chain bounds known, the proper limits of atoms
                                    ///<   found in chain_limits can be accessed, and from those
                                    ///<   the lists of atoms in each cell of the chain.
  const int* chain_system_owner;    ///< Indices of the systems that own each chain.  The system
                                    ///<   which owns the kth chain of cells is given by the value
                                    ///<   of chain_system_owner[k].
  T4* image;                        ///< The array of atom coordinates and property information.
                                    ///<   This is modifiable so that it can be built by one of a
                                    ///<   collection of similar kernels.
  const T4* image_old;              ///< The old array of atom coordinates, readable but taken as
                                    ///<   fixed so that the new image can be built out-of-place.
                                    ///<   See the time cycle to understand how this and image
                                    ///<   point to image or image_alt (and vice-versa) in the
                                    ///<   original CellGrid object.
  uint2* entry_room;                ///< Holding space for atoms that will enter each cell chain.
                                    ///<   The "x" member of each tuple holds the index in
                                    ///<   image_old where atom is coming from.  The "y" member
                                    ///<   indicates the destination cell within the chain.
  uint2* exit_room;                 ///< Holding space for atoms that are exiting each cell chain.
                                    ///<   The "x" member of each tuple holds the index in
                                    ///<   image_old of the atom that is leaving and shall not
                                    ///<   appear in image.  The "y" member indicates the cell
                                    ///<   within the chain which is affected by the departure.
  const int* transit_room_bounds;   ///< The fixed demarcations defining each cell chain's entry
                                    ///<   and exit areas
  int* entry_room_counts;           ///< Numbers of atoms waiting in each cell chain's entry area
  int* exit_room_counts;            ///< Numbers of atoms waiting in each cell chain's exit area
  int* nonimg_atom_idx;             ///< Indices of each atom in the image within the associated
                                    ///<   synthesis of molecular systems (a PhaseSpaceSynthesis as
                                    ///<   well as an AtomGraphSynthesis, where atoms are held in
                                    ///<   topological order).  The value of nonimg_atom_idx[k]
                                    ///<   indicates what atom out of the PhaseSpaceSynthesis and
                                    ///<   corresponding AtomGraphSyntehsis the kth atom of the
                                    ///<   image really is.
  uint* img_atom_idx;               ///< Indices where each atom of the associated synthesis of
                                    ///<   molecular systems can be found in the image.  To figure
                                    ///<   out the image array location of the kth atom of some
                                    ///<   PhaseSpaceSynthesis or AtomGraphSynthesis, look up
                                    ///<   img_atom_idx[k].
  uint* exclusion_maps;             ///< Bit masks indicating the exclusions of atom pairs in each
                                    ///<   cell based on a tower-plate Neutral Territory
                                    ///<   decomposition.
  const uint* exclusion_maps_old;   ///< Exclusion bit masks for the outdated arrangement
  const int* nt_groups;             ///< Concatenated lists of cells associated by the "tower and
                                    ///<   plate" neutral territory decomposition.  The ith group
                                    ///<   of 16 elements from this list enumerates the cells to
                                    ///<   associate with the ith cell, whose atom limits are
                                    ///<   described in the cell_limits or cell_limits_old arrays.
  Tacc* xfrc;                       ///< Accumulators for Cartesian X forces on imaged particles
  Tacc* yfrc;                       ///< Accumulators for Cartesian Y forces on imaged particles
  Tacc* zfrc;                       ///< Accumulators for Cartesian Z forces on imaged particles
  int* xfrc_ovrf;                   ///< Overflow bits for xfrc
  int* yfrc_ovrf;                   ///< Overflow bits for yfrc
  int* zfrc_ovrf;                   ///< Overflow bits for zfrc
  Tacc* xfrc_hw;                    ///< Workspace for warp-specific Cartesian X force accumulation
  Tacc* yfrc_hw;                    ///< Workspace for warp-specific Cartesian Y force accumulation
  Tacc* zfrc_hw;                    ///< Workspace for warp-specific Cartesian Z force accumulation
  int* xfrc_hw_ovrf;                ///< Overflow bits for xfrc_hw
  int* yfrc_hw_ovrf;                ///< Overflow bits for yfrc_hw
  int* zfrc_hw_ovrf;                ///< Overflow bits for zfrc_hw
};

/// \brief Read-only abstract for the CellGrid object.  This is unused in typical MD applications,
///        when it is always the case that the cell grid positions and atom content are being
///        updated or forces are being accumulated, but it may be useful when the cell grid is in
///        demand strictly as a directory of which atoms are close to one another for building or
///        detailing some other object.  The associated workspacec variables are not included in
///        this form of the abstract.
template <typename T, typename Tacc, typename Tcalc, typename T4> struct CellGridReader {

  /// \brief As is other abstracts, the constructor takes a list of inputs for all member
  ///        variables.
  /// \{
  CellGridReader(NonbondedTheme theme_in, int system_count_in, int total_cell_count_in,
                 int total_chain_count_in, int mesh_ticks_in, size_t cell_base_capacity_in,
                 size_t cell_excl_capacity_in, float lpos_scale_in, float lpos_inv_scale_in,
                 const ullint* system_cell_grids_in, const Tcalc* system_cell_umat_in,
                 const T* system_cell_invu_in, const T* system_pmig_invu_in,
                 const uint2* cell_limits_in, const uint* chain_limits_in,
                 const int* system_chain_bounds_in, const int* chain_system_owner_in,
                 const T4* image_in, const int* nonimg_atom_idx_in, const uint* img_atom_idx_in,
                 const uint* exclusion_maps_in, const int* nt_groups_in, const Tacc* xfrc_in,
                 const Tacc* yfrc_in, const Tacc* zfrc_in, const int* xfrc_ovrf_in,
                 const int* yfrc_ovrf_in, const int* zfrc_ovrf_in);

  CellGridReader(const CellGridWriter<T, Tacc, Tcalc, T4> &cgw);

  CellGridReader(const CellGridWriter<T, Tacc, Tcalc, T4> *cgw);
  /// \}

  /// \brief The const-ness of all members implicitly deletes the copy and move assignment
  ///        operators, but the default copy and move constructors are valid.
  /// \{
  CellGridReader(const CellGridReader &original) = default;
  CellGridReader(CellGridReader &&original) = default;
  /// \}

  const NonbondedTheme theme;       ///< The non-bonded property for which the neighbor list is
                                    ///<   maintained (electrostatic, van-der Waals, or both)
  const int system_count;           ///< Number of systems whose spatial decompositions are stored
  const int total_cell_count;       ///< The total number of decomposition cells across all systems
  const int total_chain_count;      ///< The number of cell chains, groups of cells found in a row
                                    ///<   along the unit cell A axis of one simulation, across all
                                    ///<   systems
  const int mesh_ticks;             ///< The number of mesh grid spacings per spatial decomposition
                                    ///<   cell, pertaining to a particle-mesh interaction grid
                                    ///<   associated with the neighbor list decomposition.
  const size_t cell_base_capacity;  ///< The maximum capacity of any one cell (all cells in a row
                                    ///<   along a simulation box's A axis share their excess
                                    ///<   capacity, making a provision for one cell to take on
                                    ///<   a much larger number of atoms)
  const size_t cell_excl_capacity;  ///< Capacity allotted for exclusions associated with the
                                    ///<   cell's neutral territory interaction assignments
  const float lpos_scale;           ///< The local position scaling factor applied to coordinates
                                    ///<   in the images (expressed as a 32-bit floating point
                                    ///<   number, as this will be converted to 64-bit format with
                                    ///<   no more or less information if needed)
  const float lpos_inv_scale;       ///< The inverse scaling factor applied to coordinates in the
                                    ///<   images
  const ullint* system_cell_grids;  ///< Dimensions of each cell grid and the starting point for
                                    ///<   its cells in the array
  const Tcalc* system_cell_umat;    ///< Transformation matrices for taking Cartesian coordinates
                                    ///<   into the fractional space of the local coordinate system
                                    ///<   in each system's individual cells.
  const T* system_cell_invu;        ///< Inverse transformation matrix for each system's individual
                                    ///<   cells.  The columns of this matrix indicate the size of
                                    ///<   the cell.
  const T* system_pmig_invu;        ///< Inverse transformation matrix for each system's individual
                                    ///<   particle-mesh interaction grid elements.  The columns of
                                    ///<   this matrix indicate the size of the cell.
  const uint2* cell_limits;         ///< Limits of the atoms in each cell out of the entire image.
                                    ///<   The "x" and member of the tuple provides an absolute
                                    ///<   index within the current image for the lower bound of
                                    ///<   atoms.  The "y" member provides, in its low 16 bits, the
                                    ///<   system index (the maximum number of systems is set in
                                    ///<   maximum_cellgrid_systems).  The "y" member of the tuple
                                    ///<   also provides, in its high 16 bits, the number of atoms
                                    ///<   in the particular cell.  The system index must be used
                                    ///<   to retrieve information in system_cell_grids in order to
                                    ///<   put the cell in context and see, for example, which
                                    ///<   other cells are its neighbors.
  const uint* chain_limits;         ///< Immutable limits on the cell chains within each cell grid.
                                    ///<   The ith chain cannot hold more than a number of atoms
                                    ///<   given by chain_limits[i + 1] - chain_limits[i].  The
                                    ///<   atom capacity across chains in a particular system is
                                    ///<   consistent, but different systems may have different
                                    ///<   capacities.
  const int* system_chain_bounds;   ///< Bounds on the list of chains held by each system.  The
                                    ///<   process is for one or more blocks to look up these
                                    ///<   limits and then do chains from within the system before
                                    ///<   looking up new limits and moving on to the next block.
                                    ///<   With the chain bounds known, the proper limits of atoms
                                    ///<   found in chain_limits can be accessed, and from those
                                    ///<   the lists of atoms in each cell of the chain.
  const int* chain_system_owner;    ///< Indices of the systems that own each chain.  The system
                                    ///<   which owns the kth chain of cells is given by the value
                                    ///<   of chain_system_owner[k].
  const T4* image;                  ///< The array of atom coordinates and property information
  const int* nonimg_atom_idx;       ///< Indices of each atom in the image within the associated
                                    ///<   synthesis of molecular systems (a PhaseSpaceSynthesis as
                                    ///<   well as an AtomGraphSynthesis, where atoms are held in
                                    ///<   topological order).  The value of nonimg_atom_idx[k]
                                    ///<   indicates what atom out of the PhaseSpaceSynthesis and
                                    ///<   corresponding AtomGraphSyntehsis the kth atom of the
                                    ///<   image really is.
  const uint* img_atom_idx;         ///< Indices where each atom of the associated synthesis of
                                    ///<   molecular systems can be found in the image.  To figure
                                    ///<   out the image array location of the kth atom of some
                                    ///<   PhaseSpaceSynthesis or AtomGraphSynthesis, look up
                                    ///<   img_atom_idx[k].
  const uint* exclusion_maps;       ///< Bit masks indicating the exclusions of atom pairs in each
                                    ///<   cell based on a tower-plate Neutral Territory
                                    ///<   decomposition.
  const int* nt_groups;             ///< Concatenated lists of cells associated by the "tower and
                                    ///<   plate" neutral territory decomposition.  The ith group
                                    ///<   of 16 elements from this list enumerates the cells to
                                    ///<   associate with the ith cell, whose atom limits are
                                    ///<   described in the cell_limits or cell_limits_old arrays.
  const Tacc* xfrc;                 ///< Accumulated Cartesian X forces on all imaged particles
  const Tacc* yfrc;                 ///< Accumulated Cartesian Y forces on all imaged particles
  const Tacc* zfrc;                 ///< Accumulated Cartesian Z forces on all imaged particles
  const int* xfrc_ovrf;             ///< Overflow bits for xfrc
  const int* yfrc_ovrf;             ///< Overflow bits for yfrc
  const int* zfrc_ovrf;             ///< Overflow bits for zfrc
};

/// \brief An object to manage the spatial decomposition of a system of particles.  The general
///        strategy is to arrange particles in a grid of cells, each at least half the direct-space
///        cutoff in width in all directions.  A work unit entails computing all interactions
///        assigned to one of the cells.  This is done using a neutral-territory decomposition
///        the lines of the "tower-plate" arrangement, which consists of seventeen cells:
///
///                  ---X--- ^                    ---X--- ^                    ------- ^
///                  ---X--- |                    ---X--- |                    ------- |
///    View along A: ---XXX- c      View along B: -XXXXX- c      View along C: ---XXX- b
///                  ---X---                      ---X---                      -XXXXX-
///                  ---X--- b->                  ---X--- a->                  -XXXXX- a->
///
/// As with other grids in STORMM, the cell grid index increments most rapidly along the A axis.
/// Furthermore, the atoms of all cells in a stack along A are arranged in contiguous memory, with
/// no spacing between them.  Each row of cells along a mesh's A axis is arranged with a degree of
/// padding so that the population of any particular row may fluctuate within some anticipated
/// bounds.  Memory access is optimized by arranging the tower along the C axis, providing
/// continuous atom flow for 3, 5, and 5 cells when loading the larger plate.  Atoms in the tower
/// will interact with the most other atoms and require one global write (atomicAdd) to commit
/// their forces back to the totals at the end of the work unit.  Therefore, it is most efficient
/// to have these atoms be the least contiguous of the whole.
///
/// In all cells, cache behavior is optimized by storing coordinates and properties as four-tuples
/// such that only eight atoms in single-precision mode (or four in double-precision mode) are
/// required to fill a cache line.
///
/// The work cycle is as follows:
///   1 .) A KERNEL COMPUTING VALENCE INTERACTIONS will accumulate all forces, including non-bonded
///        contributions originating in this kernel on a previous cycle, and move particles.
///        Following moves, this VALENCE INTERACTION KERNEL will compute new cell addresses for
///        each atom and populate the vestibules of particles entering and leaving each cell.  The
///        VALENCE INTERACTION KERNEL will also deposit the new, scaled coordinates of each cell
///        in the previous cycle's atoms array, which is no longer in use.
///   2.)  A NEIGHBOR LIST UPDATE KERNEL will fire off ahead of the non-bonded kernels.  This will
///        handle re-imaging of coordinates and updates to cells.  Taking the (now complete)
///        vestibules of entering and exiting particles, the kernel begins by performing an
///        out-of-place sort to adjust each cell in the mesh.  Each warp performs a prefix sum over
///        the adjustments in the particular cell's column in order to arrive at the limits within
///        which to transfer particles, first from the cell's contents in the old atoms array (so
///        long as they are to be retained) and then by referencing the prior indices of particles
///        that are moving in and transferring the particles from their prior cells, but will
///        likely be hidden by the memory accesses happening at the same time.  Particles exit
///        cells, causing contractions of the list of bitmasks and then of the bitmasks themselves.
///        Particles leaving the plate cell prunes the array of bitmasks while particles leaving
///        the tower prunes the bitmasks themselves--see notes in this library's template
///        implementation file.  Likewise, particles entering the tower cells will add bits to the
///        back of each mask, and particles entering cells of the plate will expand the arrays of
///        bitmasks in various places.  The exact list of particles in each cell will be in flux,
///        but with the identities of particles entering and the indices (positions within the
///        cell grid image) of particles exiting already established in the preceding kernel, the
///        non-bonded lists can be updated without risk of any race conditions.  One warp will be
///        tasked with updating each cell's atom coordinates while a second will be tasked with
///        updating its exclusion masks.  Both the atom coordinate update and the neighbor list
///        update are performed out-of-place.
///   3 .) Taking the newly updated neighbor lists, the RECIPROCAL-SPACE KERNEL will proceed to
///        assemble the the FFT pencils and use strong barriers in order to complete the
///        convolution.
///   4 .) The DIRECT-SPACE KERNEL will launch to complete the work, and perform interpolation of
///        long-ranged particle-particle interactions if that is required. 
template <typename T, typename Tacc, typename Tcalc, typename T4> class CellGrid {
public:

  /// \brief The constructor accepts an associated PhaseSpaceSynthesis object and a precision
  ///        model.  The system count and scaling factors will be set by the coordinate synthesis.
  ///        The total cell count and cell layout can be modulated by specifying a cutoff and
  ///        a padding related to how much larger each cell should be sized.
  ///
  /// \param cutoff                 The minimum distance between any two cell sides
  /// \param padding                An amount added to cutoff in the determination of the minimum
  ///                               cell size
  /// \param mesh_subdivisions_in   The number of mesh elements along each axis of one of the
  ///                               spatial decomposition cell.  This value, times the cell grid
  ///                               dimensions themselves, defines the size of a particle-mesh
  ///                               interaction grid associated with the neighbor list, i.e. the
  ///                               PME reciprocal space grid.
  /// \param theme_in               The non-bonded potential that will be computed based on the
  ///                               neighbor list represented in the cell grid
  /// \param gpu                    Details of the GPU that will utilize and manipulate the cell
  ///                               grid
  /// \param cell_base_capacity_in  The atom capacity of an individual decomposition cell.  This
  ///                               is set to be larger than the atom content is likely to go, and
  ///                               in all cases where it might be relevant multiple cells will
  ///                               pool their capacities to ensure that overflow is extremely
  ///                               unlikely.
  /// \param policy_in              Course of action in the case of bad input
  /// \{
  CellGrid(const PhaseSpaceSynthesis *poly_ps_ptr_in, const AtomGraphSynthesis *poly_ag_ptr_in,
           double cutoff, double padding, int mesh_subdivisions_in, NonbondedTheme theme_in,
           const GpuDetails &gpu = null_gpu,
           size_t cell_base_capacity_in = default_cellgrid_base_capacity,
           ExceptionResponse policy_in = ExceptionResponse::WARN);

  CellGrid(const PhaseSpaceSynthesis &poly_ps_ptr_in, const AtomGraphSynthesis &poly_ag_ptr_in,
           double cutoff, double padding, int mesh_subdivisions_in, NonbondedTheme theme_in,
           const GpuDetails &gpu = null_gpu,
           size_t cell_base_capacity_in = default_cellgrid_base_capacity,
           ExceptionResponse policy_in = ExceptionResponse::WARN);
  /// \}

  /// \brief The default copy and move constructors, as well as assignment operators, apply so long
  ///        as there are no POINTER-kind Hybrid objects to repair and no const members.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right hand side of the assignment statement
  /// \{
  CellGrid(const CellGrid &original) = default;
  CellGrid(CellGrid &&original) = default;
  CellGrid& operator=(const CellGrid &original) = default;
  CellGrid& operator=(CellGrid &&original) = default;
  /// \}

  /// \brief Get the number of systems in the object.
  int getSystemCount() const;

  /// \brief Get the total number of cells allocated in the object.
  int getTotalCellCount() const;

  /// \brief Get the number of cells allotted to any one system.
  ///
  /// Overloaded:
  ///   - Get the total number of cells in the system.
  ///   - Get the number cells along one of the system's simulation box axes.
  ///
  /// \param index  The system of interest
  /// \param axis   The axis of interest
  /// \{
  int getCellCount(int index) const;
  int getCellCount(int index, UnitCellAxis axis) const;
  int getCellCount(int index, CartesianDimension axis) const;
  /// \}

  /// \brief Get the base capacity of any particular cell in the grid.  All systems' cells share
  ///        the same base capacity.
  size_t getCellBaseCapacity() const;

  /// \brief Get the effective cutoff.
  double getEffectiveCutoff() const;

  /// \brief Get the non-bonded theme of this cell grid--does it hold atoms with only electrostatic
  ///        or van-der Waals properties?
  NonbondedTheme getTheme() const;
  
  /// \brief Get the position in the time cycle.
  CoordinateCycle getCyclePosition() const;

  /// \brief Get the number of bits in the fixed-precision format.
  int getPositionScalingBits() const;

  /// \brief Get the number of mesh elements spanning each spatial decomposition cell.  The "mesh"
  ///        refers to a particle-mesh interaction grid, e.g. the PME reciprocal space lattice.
  int getMeshSubdivisions() const;

  /// \brief Find the image location of an atom based on its system and topological atom number.
  ///
  /// Overloaded:
  ///   - Check in the most recent image
  ///   - Specify the primary or alternate image to check for the atom
  ///
  /// \param system_index  The system of interest within the associated synthesis
  /// \param atom_index    Topological atom index within the system
  /// \param orientation   Specify whether to look in the primary or
  /// \{
  uint getImageIndex(int system_index, int atom_index, CoordinateCycle orientation) const;
  uint getImageIndex(int system_index, int atom_index) const;
  /// \}

  /// \brief Find the cell grid location of an atom based on its system and topological atom
  ///        number.  The result will be returned as the index of the cell within the grid
  ///        belonging to the system of interest, with the indices along the A, B, and C axes
  ///        placed in the "x", "y", and "z" members of the tuple.  The overall index of the cell
  ///        in the entire CellGrid object (up to 2^28 cells) will be returned in the "w" member of
  ///        the tuple.  All of this information is available in bit-packed form in the
  ///        system_cell_limits array.  Overloading and descriptions of parameters follow from
  ///        getImageIndex() above.
  /// \{
  int4 getCellLocation(int system_index, int atom_index, CoordinateCycle orientation) const;
  int4 getCellLocation(int system_index, int atom_index) const;
  /// \}

  /// \brief Get the chain index of an atom based on its system and topological atom index.
  ///        The chain index is returned in the "x" member of the tuple while the offset within the
  ///        chain is returned in the "y" member.  Overloading and descriptions of parameters
  ///        follow from getImageIndex() above. 
  /// \{
  int2 getChainLocation(int system_index, int atom_index, CoordinateCycle orientation) const;
  int2 getChainLocation(int system_index, int atom_index) const;
  /// \}
  
  /// \brief Get a const pointer to the coordinate synthesis served by this object.
  const PhaseSpaceSynthesis* getCoordinateSynthesisPointer() const;

  /// \brief Get a const pointer to the topology synthesis served by this object.
  const AtomGraphSynthesis* getTopologySynthesisPointer() const;

  /// \brief Obtain the object's abstract in order to access its members in C-programming fashion,
  ///        whether on the CPU host or GPU device.
  ///
  /// Overloaded:
  ///   - Return a const reader from a const CellGrid object
  ///   - Return an editable writer from a non-const CellGrid object
  ///   - Indicate the point in the time cycle or take the device's internal setting
  ///
  /// \param orientation  The point in the time cycle at which to obtain some pointers
  /// \param tier         Indicate whether to target pointers to memory on the CPU host or GPU
  ///                     device
  /// \{
  CellGridWriter<T, Tacc, Tcalc, T4> data(CoordinateCycle orientation,
                                          HybridTargetLevel tier = HybridTargetLevel::HOST);

  CellGridWriter<T, Tacc, Tcalc, T4> data(HybridTargetLevel tier = HybridTargetLevel::HOST);

  const CellGridReader<T, Tacc, Tcalc, T4>
  data(CoordinateCycle orientation, HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  const CellGridReader<T, Tacc, Tcalc, T4>
  data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Obtain the object's abstract with templating removed.  Overloading and descriptions of
  ///        parameters follow from data() above.
  /// \{
  CellGridWriter<void, void, void, void>
  templateFreeData(CoordinateCycle orientation, HybridTargetLevel tier = HybridTargetLevel::HOST);

  CellGridWriter<void, void, void, void>
  templateFreeData(HybridTargetLevel tier = HybridTargetLevel::HOST);

  const CellGridReader<void, void, void, void>
  templateFreeData(CoordinateCycle orientation,
                   HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  const CellGridReader<void, void, void, void>
  templateFreeData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Get a pointer to the object itself, useful if the object is only available by const
  ///        reference.
  const CellGrid<T, Tacc, Tcalc, T4>* getSelfPointer() const;
  
#ifdef STORMM_USE_HPC
  /// \brief Upload all data to the device.
  void upload();

  /// \brief Download all data from the device.
  void download();
#endif

  /// \brief Initialize forces for the cell grid, on the CPU host or GPU device.  This standalone
  ///        feature provides a means for performing this activity at will, although in the most
  ///        optimized code the process it performs will be fused with some other kernel.  This
  ///        will initialize forces for the current cell layout and should only be called once the
  ///        contents of image_cell_limits or image_cell_limits_alt have been settled, whichever is
  ///        appropriate to the current time cycle.  There is only one array of forces in the 
  ///        CellGrid, as once they are accumulated in this object they are added to one of the two
  ///        time-cycle dependent arrays in the accompanying PhaseSpaceSynthesis.
  ///
  /// \param tier  Indicate whether to initialize forces on the CPU host or GPU device
  /// \param gpu   Details of the GPU that will perform the initialization of device memory
  void initializeForces(HybridTargetLevel tier = HybridTargetLevel::HOST,
                        const GpuDetails &gpu = null_gpu);

  /// \brief Contribute forces accumulated in the cell grid back to a coordinate synthesis with
  ///        its own force accumulators and forces from other sources.
  ///
  /// Overloaded:
  ///   - Contribute the forces to an arbitrary PhaseSpaceSynthesis
  ///   - Contribute the forces back to the CellGrid's internally referenced synthesis
  ///
  /// \param dest  Specify a different synthesis to receive the accumulated forces
  /// \param tier  Indicate whether to initialize forces on the CPU host or GPU device
  /// \param gpu   Details of the GPU that will perform the initialization of device memory
  /// \{
  void contributeForces(PhaseSpaceSynthesis *dest, HybridTargetLevel tier,
                        const GpuDetails &gpu) const;

  void contributeForces(HybridTargetLevel tier, const GpuDetails &gpu) const;
  /// \}

private:
  int system_count;           ///< The total number of distinct systems being simulated
  int total_cell_count;       ///< Total number of spatial decomposition cells across all systems
  int total_chain_count;      ///< The total number of columns across all systems, each of them
                              ///<   a free space in which the cells along the chain arrange their
                              ///<   atoms in fully contiguous memory
  size_t cell_base_capacity;  ///< The maximum atom load that any one cell is allocated to handle
  size_t cell_excl_capacity;  ///< The allotment of 32-bit unsigned int bit masks given to each
                              ///<   cell for tracking exclusions in its assigned neutral territory
                              ///<   decomposition
  double effective_cutoff;    ///< The combination of cutoff and padding that contributes to the
                              ///<   minimum distance measurement for each cell
  int mesh_subdivisions;      ///< The number of discretizations of the associated particle-mesh
                              ///<   grid, per spatial decomposition cell.  A common value is four.
  ExceptionResponse policy;   ///< Protocol in the event of bad input or any need to grow beyond
                              ///<   stipulated bounds
  NonbondedTheme theme;       ///< The theme of the neighbor list, to which the interactions it
                              ///<   catalogs pertains

  /// Like the PhaseSpace and PhaseSpaceSynthesis objects, this spatial decomposition incorporates
  /// the concept of a step cycle for out-of-place rebuilding.
  CoordinateCycle cycle_position;

  // The local coordinate systems for each cell are tailored to maximize the amount of information
  // that can be encoded in the fixed-precision representation with the data type implied by the
  // precision member variable (PrecisionModel::SINGLE implies 32-bit ints, PrecisionModel::DOUBLE
  // implies 64-bit long long ints).  The number of fixed precision bits is designed to accommodate
  // stacking up to two cells against the home cell along any of the simulation box axes and still
  // have any coordinates within those boxes be represented by the fixed-precision format without
  // overflow.  Furhermore, the format is designed to accommodate moderate stretching of the
  // simulation box (up to 5% in any direction).  If the simulation expands significantly, this
  // will trigger the cell grid to be rebuilt.
  int localpos_scale_bits;        ///< Number of bits after the decimal in local fixed-precision
                                  ///<   coordinate representations within each cell
  double localpos_scale;          ///< Scaling factor for taking Cartesian coordinates into the
                                  ///<   local fixed-precision model for each cell
  double localpos_inverse_scale;  ///< The inverse scaling factor for local coordinates
  
  /// The limits of cells in each system, a bounds array of sorts on sp_atoms or atoms below.
  /// The ith entry is a bit-packed string which indicates, in its low 28 bits, the point in the
  /// cell list at which cells contents begin to pertain to the ith simulation.  In the next three
  /// stretches of 12 bits, the dimensions of the cell grid along the a, b, and c box vectors are
  /// provided with a maximum dimension of 4095 cells in any given direction.
  Hybrid<ullint> system_cell_grids;

  /// Transformation matrices for taking particles into the local coordinate fractional space of
  /// the spatial decomposition cell.  These matrices are the inverse of matrices found in
  /// system_cell_invu and system_cell_invu_alt, below, and needed when computing the alignment of
  /// particles to an associated particle-mesh interaction grid.
  Hybrid<Tcalc> system_cell_umat;
  Hybrid<Tcalc> system_cell_umat_alt;
  
  /// Cell dimensions represented in the precision of the cell grid, whether real-values or
  /// fixed-precision integers.  Before each cell-to-cell comparison, particles will be shifted by
  /// adding these amounts in increments of -1, 0, or +1 depending on the relative position of the
  /// other cell relative to the home cell.
  /// \{
  Hybrid<T> system_cell_invu;
  Hybrid<T> system_cell_invu_alt;
  /// \}

  /// Dimensions of the particle-mesh interaction grid element, represented in the precision of the
  /// cell grid, whether real-values or fixed-precision integers.  These matrices can be useful for
  /// improving the precision by which particle-mesh interaction fractional coordinates are
  /// obtained, if there is a fixed-precision representation of the coordinates.
  /// \{
  Hybrid<T> system_pmig_invu;
  Hybrid<T> system_pmig_invu_alt;
  /// \}

  /// Dimensions of the element of the associated particle-mesh interaction grid, represented in
  /// the precision of the cell grid.  When representing the coordinates in fixed precision, this
  /// transformation matrix can improve the precision by which the fractional value used to compute
  /// particle-mesh mapping coefficients, i.e. B-spline coefficients, is derived.
  Hybrid<T> particle_mesh_grid_invu;
  Hybrid<T> particle_mesh_grid_invu_alt;
  
  /// The number of atoms in each cell is limited by the cell base capacity and, beyond that, the
  /// number of cells along the simulation box A axis which can collectively buffer a transient
  /// high particle count.  These arrays are each total_cell_count in length and detail, for the
  /// WHITE and BLACK images, the starting point of atoms in the image (the "x" member)
  /// followed by the total number of atoms in the cell (the high 16 bits of the "y" member) and
  /// the system index (the low 16 bits of the "y" member).
  /// \{
  Hybrid<uint2> image_cell_limits;
  Hybrid<uint2> image_cell_limits_alt;
  /// \}

  /// While the instantaneous limits of each cell will vary with time (the primary and alternate
  /// images will also have different counts), the absolute column limits are fixed when the object
  /// is allocated.  These apply to both primary and alternate images.
  Hybrid<uint> image_chain_limits;

  /// Each system can have a unique number of chains.  This array stores the offset at which each
  /// system's column counts begin, and acts as a sort of bounds array into image_column limits
  /// that provides information to supplement system_cell_grids, above.  These limits apply to both
  /// primary and alternate images.
  Hybrid<int> system_chain_bounds;

  /// It is important to be able to perform the reverse lookup of the above: to what system does
  /// the kth chain of the overall list belong?  This array, which like system_chain_bounds is
  /// constant once the object is created, will provide the system index of the kth overall chain
  /// in its kth index.
  Hybrid<int> chain_system_membership;
  
  /// The cell's atom positions and property data is encoded in one array.  Atom coordinates are
  /// stored relative to the origin of their particular spatial decomposition cell in
  /// fixed-precision format.  For each cell read, the appropriate offset will be applied to shift
  /// the coordinates +/-1 cells along the a, b, or c unit cell axes, after which they will be
  /// converted to real numbers.
  ///
  /// Each tuple contains the relative Cartesian X, Y, and Z coordinates in the "x", "y", and "z"
  /// members, plus the property in the "w" member.  Integers are used for fixed-precision
  /// representations.  Lennard-Jones type indices are represented as integers, of course, but
  /// atomic partial charges are converted to bitstrings in the appropriate precision and then to
  /// integers.
  ///
  /// Maps into the exclusions for each cell are handled by 32-bit unsigned ints, starting at some
  /// offset.  The ith integer after the offset describes exclusions involving the ith atom of the
  /// adjacent cell, and the jth bit in the uint indicates whether the interaction with the jth
  /// atom of the home cell is an exclusion.  Space is allocated to map exclusions between any cell
  /// with a normal atom load interacting with its neighbors: the cell_base_capacity variable grows
  /// in steps of 32, and the exclusions map space for the cell grows as
  ///
  ///   map space = (cell_base_capacity * 14 * cell_base_capacity / 32) uint 
  ///
  /// There are no "extra large" cells but each column is allocated to store more atoms than should
  /// be required to fit in the column as a whole at any one time.  If, by some stroke of fate, a
  /// column is left with more atoms than the allocation allows, this will trigger an error and
  /// eventually alert the user to ask for more column memory in the input.  Such situations will
  /// be extremely rare, and probably never occur in standard MD operations.
  /// \{
  Hybrid<T4> image;
  Hybrid<T4> image_alt;
  /// \}

  /// As the structure evolves over the course of a simulation, each cell will have atoms coming
  /// and going.  The key is to manage this flow, and conserve work done on neighbor lists between
  /// cells from step to step.  The following two arrays manage atoms entering and leaving each
  /// cell.  Each migrating atom will be catalogged in both arrays.  The process will be to list
  /// the indices of atoms in the other image, i.e. the WHITE image entering and leaving
  /// vestibules reference atom indices in the BLACK image, and vice-versa.  The "x" member
  /// of each tuple holds the index in question while the "y" member holds the index of the cell
  /// that is being entered or left.  Each entering and leaving array will be sized based on the
  /// expected movement of atoms and the density of atoms in each cell.
  /// \{
  Hybrid<uint2> image_chain_entering;
  Hybrid<uint2> image_chain_leaving;
  /// \}

  /// Like the chains of cells themselves, the entering and leaving vestibules of any one system
  /// are grouped together so that heavy traffic into or out of one cell does not overflow the
  /// pooled space.  The offsets of any one system's vestibules are given by the same array as its
  /// cell chains, system_chain_bounds.  The limits of each system's entering and leaving
  /// vestibules, indices in the image_cell_entering and / or image_cell_leaving arrays, are given
  /// in the array below.  Like the image_chain_limits array, this is set upon construction of the
  /// object.
  Hybrid<int> image_vestibule_limits;

  /// Like the cell chains in either image, the vestibules must keep records of their total
  /// contents. (There is only one collections of entering and leaving vestibules, serving as a
  /// transition between the primary and alternate images.) These arrays provide that tracking.
  /// \{
  Hybrid<int> image_chain_entering_counts;
  Hybrid<int> image_chain_leaving_counts;
  /// \}

  /// Atom indices in the associated PhaseSpaceSynthesis object.  The format of this array tracks
  /// the layout of atoms in the image.
  /// \{
  Hybrid<int> nonimaged_atom_indices;
  Hybrid<int> nonimaged_atom_indices_alt;
  /// \}

  /// Locations of each atom in the synthesis within the spatial decomposition.  These are
  /// absolute indices in to the image array.  The need to represent absolute indices in these
  /// arrays, while keeping memory bandwidth requirements down, guides the choice of an unsigned
  /// int here (as well as in the cell limits arrays).  Some degree of padding is needed in each
  /// cell, which combines into a buffer of the atom population for all cells along each
  /// simulation's A axis.
  /// \{
  Hybrid<uint> image_array_indices;
  Hybrid<uint> image_array_indices_alt;
  /// \}

  /// Exclusion maps for each cell's associated pair interactions
  /// \{
  Hybrid<uint> exclusion_maps;
  Hybrid<uint> exclusion_maps_alt;
  /// \}

  /// Work groups indicating the cells associated with any given cell in order to complete the
  /// neutral territory interactions by the "tower and plate" method.
  Hybrid<int> nt_work_groups;
  
  /// Charts to help in the navigation of each exclusion map for a any given cell.  There is one
  /// stride of 32 numbers (not a matter of warp size) for each home cell, which can be taken to
  /// have coordinates Ha, Hb, and Hc in the grid for its particular system.  The data in each
  /// stride describes:
  ///
  ///     Pos.
  ///   -  0- 4 : Total number of atoms in the tower cells, from bottom to top
  ///   -  5, 6 : Total number of atoms in the plate cells (Ha+1,Hb,Hc) and (Ha+2,Hb,Hc) from the
  ///             home cell located at (Ha,Hb,Hc)
  ///   -  7-11 : Total number of atoms in the plate cells (Ha-2,Hb+1,Hc) to (Ha+2,Hb+1,Hc)
  ///   - 12-16 : Total number of atoms in the plate cells (Ha-2,Hb+2,Hc) to (Ha+2,Hb+2,Hc)
  ///   -    17 : Total number of atoms in the plate
  
  
  // Force accumulation arrays track the total number of atoms in the atoms or sp_atoms list.
  // These accumulators are split, with the primary accummulators being llint for double-precision
  // and int for single-precision.  The overflow accumulators are used in both cases, even if
  // accesses to them are rare.
  Hybrid<Tacc> x_force;          ///< Main Cartesian X force accumulators in double-precision mode
  Hybrid<Tacc> y_force;          ///< Main Cartesian Y force accumulators in double-precision mode
  Hybrid<Tacc> z_force;          ///< Main Cartesian Z force accumulators in double-precision mode
  Hybrid<int> x_force_overflow;  ///< Overflow accumulators for Cartesian X forces in either mode
  Hybrid<int> y_force_overflow;  ///< Overflow accumulators for Cartesian Y forces in either mode
  Hybrid<int> z_force_overflow;  ///< Overflow accumulators for Cartesian Z forces in either mode

  // Each warp requires a workspace for accumulating forces on the home cell when running on the
  // GPU, and the object comes prepared with space to store temporary values in L1 prior to
  // dumping them back to the main accumulators above.
  Hybrid<Tacc> warp_x_work;          ///< Workspace for GPU warps in double-precision mode,
                                     ///<   accumulating Cartesian X forces on particles in the
                                     ///<   home cell
  Hybrid<Tacc> warp_y_work;          ///< Workspace for GPU warps in double-precision mode
                                     ///<   accumulating Cartesian Y forces on particles in the
                                     ///<   home cell
  Hybrid<Tacc> warp_z_work;          ///< Workspace for GPU warps in double-precision mode
                                     ///<   accumulating Cartesian Z forces on particles in the
                                     ///<   home cell
  Hybrid<int> warp_x_overflow_work;  ///< Overflow workspace for warps in either mode, handling
                                     ///<   Cartesian X forces in the home cell
  Hybrid<int> warp_y_overflow_work;  ///< Overflow workspace for warps in either mode, handling
                                     ///<   Cartesian Y forces in the home cell
  Hybrid<int> warp_z_overflow_work;  ///< Overflow workspace for warps in either mode, handling
                                     ///<   Cartesian Z forces in the home cell
  
  /// The coordinate synthesis which this object serves
  PhaseSpaceSynthesis *poly_ps_ptr;

  /// The topology synthesis describing atom proprties of the cell grid image.  While it is
  /// possible to construct the cell grid image without the topology synthesis under certain
  /// circumstances (using the pointers to individual topologies held by the coordinate synthesis),
  /// having the pre-constructed topology synthesis is necessary to ensure that it can be
  /// constructed in modes where the cell grid stores a neighbor list suitable to "all" non-bonded
  /// potential forms.
  AtomGraphSynthesis *poly_ag_ptr;

  /// \brief Enforce a system count of no more than 4096 periodic simulations.  This will always
  ///        result in an error if the limit is violated.  Enforce unit cell types with periodic
  ///        boundary conditions.
  void validateCoordinateSynthesis() const;

  /// \brief Validate the effective cutoff.  This will impose a limit on the extent of direct-space
  ///        interactions which should keep cell populations within the established bounds.
  ///
  /// \param eff_cut_in  The effective cutoff to impose
  /// \param policy_in   Course to take if the requested number is too small or too great
  double validateEffectiveCutoff(double eff_cut_in, ExceptionResponse policy_in) const;

  /// \brief Validate the index of a queried atom within a system of the associated synthesis.
  ///        This check is performed for various indexing functions peering into the cell grid
  ///        arrangement.
  ///
  /// \param system_index  The system of interest in the associated synthesis
  /// \param atom_index    The topological index of the atom within the system
  /// \param caller        Name of the calling function
  void validateAtomIndex(int system_index, int atom_index, const char* caller) const;
  
  /// \brief Compute the optimal fixed-precision model for the cell size and atomic position data
  ///        type at hand.
  ///
  /// \param invu_samples  Temporary matrix formed with the real-valued inverse transformation
  ///                      matrices for cells in each system's cell grid
  void computeFixedPrecisionModel(const std::vector<double> &invu_samples);

  /// \brief Compute the populations of each cell in the proposed layout.  This will resize the
  ///        populations array as necessary.
  ///
  /// \param cell_populations  Vector of populations for each cell in the system's grid.  Resized,
  ///                          modified, and returned.
  /// \param system_index      Index of the system from within the associated PhaseSpaceSynthesis
  /// \param na                The number of cells along the grid's A axis
  /// \param nb                The number of cells along the grid's B axis
  /// \param nc                The number of cells along the grid's C axis
  void tallyCellPopulations(std::vector<int> *cell_populations, int system_index, int na, int nb,
                            int nc);
  
  /// \brief Populate all systems in one of the images (image or image_alt).  This must be called
  ///        after boundaries on the images and cell counts have been established.  Pointers to the
  ///        appropriate image and its supplemental arrays will be set internally.
  ///
  /// \param cyc  The point in the time cycle of the image to fill
  void populateImage(CoordinateCycle cyc);

  /// \brief Prepare work groups, lists of cell indices that should be loaded in order to complete
  ///        the tower / plate configuration for a particular central (home) cell.  Each work unit
  ///        consists of 16 integers (the index of the home cells is implied by the index of the
  ///        work unit).  The relative displacements required for atoms in each cell are implied by
  ///        its place in the list:
  ///
  ///  Pos.      Description of the relative location
  ///    0    ( 0,  0, -2), the lowest cell in the tower
  ///    1    ( 0,  0, -1), just below the home cell in the tower
  ///    2    ( 0,  0, +1), just above the home cell in the tower
  ///    3    ( 0,  0, +2), just above the home cell in the tower
  ///    4    ( 1,  0,  0), just to the right of the home cell if looking down the unit cell's B
  ///         axis, in the plate
  ///    5    ( 2,  0,  0), farthest the right of the home cell in the plate
  ///  6-10   (-2, +1,  0) to (+2, +1,  0), the first row behind the tower if looking up the unit
  ///         cell's B axis, in the plate
  /// 11-15   (-2, +2,  0) to (+2, +2,  0), the second row behind the tower, in the plate
  ///
  /// The intricacy and length of serialized operations required to derive them from information in
  /// other arrays justifies making these work units, which will remain valid for the life of the
  /// object, once on the CPU.
  void prepareWorkGroups();
  
  /// \brief Construct the initial bit masks for one of the cell grid's images.
  ///
  /// \param cyc  The point in the time cycle of the image to fill
  void initializeExclusionMasks(CoordinateCycle cyc);
};

/// \brief Re-apply templated behavior to a void-casted abstract of the CellGrid object.  Various
///        overloads of this function in other libraries handle different objects.
///
/// \param rasa  The original object, a "blank slate" with templated characteristics voided
/// \{
template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridWriter<T, Tacc, Tcalc, T4> restoreType(CellGridWriter<void, void, void, void> *rasa);

template <typename T, typename Tacc, typename Tcalc, typename T4>
CellGridWriter<T, Tacc, Tcalc, T4> restoreType(CellGridWriter<void, void, void, void> &rasa);

template <typename T, typename Tacc, typename Tcalc, typename T4>
const CellGridReader<T, Tacc, Tcalc, T4>
restoreType(const CellGridReader<void, void, void, void> *rasa);

template <typename T, typename Tacc, typename Tcalc, typename T4>
const CellGridReader<T, Tacc, Tcalc, T4>
restoreType(const CellGridReader<void, void, void, void> &rasa);
/// \}

/// \brief Compute the spatial decomposition cell in which a particle belongs, and re-image the
///        particle's coordinates as appropriate.  This follows a great deal of the functionality
///        of the imageCoordinates() function in local_arrangement.h, but intercepts the
///        fractional coordinates to obtain the correct cell location before returning them to
///        real space.
///
/// Overloaded:
///   - Locate the cell for a single particle
///   - Locate the cells for a multitude of particles given in C-style arrays, Standard Template
///     Library vectors, or Hybrid objects
///
/// \param x          The particle's Cartesian X coordinate(s), modified and returned
/// \param y          The particle's Cartesian Y coordinate(s), modified and returned
/// \param z          The particle's Cartesian Z coordinate(s), modified and returned
/// \param cell_idx   The cell indices into which the particles should go (modified and returned).
///                   The index applies to a grid where the A index increments the fastest,
///                   followed by the B index, and the the C index increments the slowest (similar
///                   to other grids in STORMM).
/// \param umat       Inverse transformation matrix 
/// \param invu       Inverse transformation matrix 
/// \param cell_na    The number of cells along the simulation A axis (may be provided as double
///                   for optimization based on the manner in which it will be used)
/// \param cell_nb    The number of cells along the simulation B axis
/// \param cell_nc    The number of cells along the simulation C axis
/// \param unit_cell  Indicate the type of unit cell (ORTHORHOMBIC will exploit optimizations for
///                   rectilinear geometry)
/// \param scale      Coordinates are scaled according to this value (anything greater than 1.0
///                   will imply that there is a fixed-precision representation in effect)
/// \{
template <typename T>
int3 locateDecompositionCell(T *x, T *y, T *z, const double* umat, const double* invu,
                             double cell_na, double cell_nb, double cell_nc,
                             UnitCellType unit_cell = UnitCellType::TRICLINIC, double scale = 1.0);

template <typename T>
void locateDecompositionCell(T* x, T* y, T* z, int* cell_idx, int natom, const double* umat,
                             const double* invu, int cell_na, int cell_nb, int cell_nc,
                             UnitCellType unit_cell = UnitCellType::TRICLINIC, double scale = 1.0);

template <typename T>
void locateDecompositionCell(std::vector<T> *x, std::vector<T> *y, std::vector<T> *z,
                             std::vector<int> *cell_idx, const double* umat, const double* invu,
                             int cell_na, int cell_nb, int cell_nc,
                             UnitCellType unit_cell = UnitCellType::TRICLINIC, double scale = 1.0);

template <typename T>
void locateDecompositionCell(Hybrid<T> *x, Hybrid<T> *y, Hybrid<T> *z, Hybrid<int> *cell_idx,
                             const double* umat, const double* invu, int cell_na, int cell_nb,
                             int cell_nc, UnitCellType unit_cell = UnitCellType::TRICLINIC,
                             double scale = 1.0);
/// \}

/// \brief Compute the rate at which particles might wander into or out of a cell with the stated
///        volume.  This is intended to be an aggressive estimate, significantly more than the
///        number which could actually migrate in any one step.  The number is computed based on an
///        assumption of random diffusion in a normal distribution.
///
/// \param effective_cutoff  Estimate of the width of the cell
/// \param sigma             The Gaussian width of the particles' velocities, indicating the 
///                          amount that they migrate after one time step
double computeMigrationRate(double effective_cutoff, double sigma);

/// \brief Based on the mesh subdivision count, reduce the number of cells along any given
///        dimension to the largest product of 2, 3, 5, and 7.  These are the factors for which
///        HPC FFT libraries are optimized.  If the number of mesh subdivisions is odd, then the
///        number of cells along at least one dimension much be even.
///
/// \param cell_na       Initial estimate for the number of cells along the unit cell A axis
/// \param cell_nb       Initial estimate for the number of cells along the unit cell B axis
/// \param cell_nc       Initial estimate for the number of cells along the unit cell C axis
/// \param subdivisions  The number of subdivisions of every spatial decomposition cell, along
///                      each dimension, which will determine the particle-mesh interaction grid
///                      dimensions (cell_na * subdivisions, cell_nb * subdivisions,
///                      cell_nc * subdivisions)
int3 optimizeCellConfiguration(int cell_na, int cell_nb, int cell_nc, int subdivisions);
  
#ifdef STORMM_USE_HPC
/// \brief Launch a standalone kernel to initialize forces in the CellGrid object.
///
/// Overloaded:
///   - Provide only the CellGrid abstract for processes involving only its organization
///   - Provide a read-only abstract for the coordinate synthesis to perform selected processes
///   - Provide a writeable abstract for the coordinate synthesis to perform selected processes
///
/// \param cgw       The cell grid abstract, containing accumulator arrays and all critical bounds
/// \param cgr       The cell grid abstract, containing accumulator arrays and all critical bounds
///                  but lacking the prior image and working arrays (as the struct is read-only)
/// \param tc_mat    Data type of cell axis matrices, implying the data type of composite
///                  coordinate and atom property four-tuples
/// \param tc_acc    Data type of force accumulators
/// \param poly_psw  Writeable abstract for the coordiante synthesis (it is trusted that the bounds
///                  of systems in this object correspond to those used to construct the cell grid)
/// \param poly_psr  Read-only abstract for the coordiante synthesis (it is trusted that the bounds
///                  of systems in this object correspond to those used to construct the cell grid)
/// \param gpu       Details of the GPU that will perform the initialization
/// \param process   The action to perform (certain actions in overloads with improper inputs, e.g.
///                  lacking a writeable coordinate synthesis, will raise runtime errors)
void launchCellGridAction(CellGridWriter<void, void, void, void> *cgw, size_t tc_mat,
                          size_t tc_acc, const GpuDetails &gpu, CellGridAction process);

void launchCellGridAction(CellGridWriter<void, void, void, void> *cgw, size_t tc_mat,
                          size_t tc_acc, const PsSynthesisReader &poly_psr, const GpuDetails &gpu,
                          CellGridAction process);

void launchCellGridAction(const CellGridReader<void, void, void, void> &cgr, size_t tc_mat,
                          size_t tc_acc, PsSynthesisWriter *poly_psw, const GpuDetails &gpu,
                          CellGridAction process);
#endif
  
} // namespace energy
} // namespace stormm

#include "cellgrid.tpp"

#endif
