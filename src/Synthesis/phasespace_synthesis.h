// -*-c++-*-
#ifndef STORMM_PHASESPACE_SYNTHESIS_H
#define STORMM_PHASESPACE_SYNTHESIS_H

#include <typeinfo>
#include <typeindex>
#include <sys/types.h>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "FileManagement/file_util.h"
#include "Math/rounding.h"
#include "Numerics/split_fixed_precision.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/trajectory_enumerators.h"

namespace stormm {
namespace synthesis {

using card::GpuDetails;
using card::Hybrid;
using card::HybridTargetLevel;
using diskutil::PrintSituation;
using stmath::roundUp;
using numerics::default_globalpos_scale_bits;
using numerics::default_localpos_scale_bits;
using numerics::default_velocity_scale_bits;
using numerics::default_force_scale_bits;
using numerics::force_scale_nonoverflow_bits;
using numerics::globalpos_scale_nonoverflow_bits;
using numerics::velocity_scale_nonoverflow_bits;
using topology::AtomGraph;
using topology::UnitCellType;
using trajectory::CoordinateCycle;
using trajectory::CoordinateFileKind;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateFrameWriter;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
using trajectory::CoordinateSeriesWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceReader;
using trajectory::PhaseSpaceWriter;
using trajectory::TrajectoryKind;

/// \brief A read-only abstract for the system demarcations in the object.  This information is
///        sometimes critical, and needed on the CPU host even as the general abstract is needed on
///        the GPU device.
struct PsSynthesisBorders {

  /// \brief The constructor accepts the total number of systems as well as pointers to the number
  ///        of atoms and starting indices of each system.
  PsSynthesisBorders(int system_count_in, const int* atom_starts_in, const int* atom_counts_in);

  /// \brief Copy and move constructors--as with any object containing const members, the move
  ///        assignment operator is implicitly deleted.
  /// \{
  PsSynthesisBorders(const PsSynthesisBorders &original) = default;
  PsSynthesisBorders(PsSynthesisBorders &&other) = default;
  /// \}
  
  const int system_count;  ///< The total number of systems in the object, and the trusted length
                           ///<   of each of the arrays below
  const int* atom_starts;  ///< Starting indices for the atoms of each system
  const int* atom_counts;  ///< Atom counts for each system
};
  
/// \brief The writer for a PhaseSpaceSynthesis object, containing all of the data relevant for
///        propagating dynamics in a collection of systems.
struct PsSynthesisWriter {

  /// \brief Constructor takes a straight list of pointers and constants from the parent object
  ///        (in this case, a PhaseSpaceSynthesis), like other abstracts.
  PsSynthesisWriter(int system_count_in, int unique_topology_count_in, UnitCellType unit_cell_in,
                    const int* atom_starts_in, const int* atom_counts_in,
                    const int* common_ag_list_in, const int* common_ag_bounds_in,
                    const int* unique_ag_idx_in, const int* replica_idx_in, double gpos_scale_in,
                    double lpos_scale_in, double vel_scale_in, double frc_scale_in,
                    int gpos_bits_in, int lpos_bits_in, int vel_bits_in, int frc_bits_in,
                    llint* boxvecs_in, int* boxvec_ovrf_in, double* umat_in, double* invu_in,
                    double* boxdims_in, llint* alt_boxvecs_in, int* alt_boxvec_ovrf_in,
                    double* umat_alt_in, double* invu_alt_in, double* alt_boxdims_in,
                    llint* xcrd_in, llint* ycrd_in, llint* zcrd_in, int* xcrd_ovrf_in,
                    int* ycrd_ovrf_in, int* zcrd_ovrf_in, llint* xvel_in, llint* yvel_in,
                    llint* zvel_in, int* xvel_ovrf_in, int* yvel_ovrf_in, int* zvel_ovrf_in,
                    llint* xfrc_in, llint* yfrc_in, llint* zfrc_in, int* xfrc_ovrf_in,
                    int* yfrc_ovrf_in, int* zfrc_ovrf_in, llint* xalt_in, llint* yalt_in,
                    llint* zalt_in, int* xalt_ovrf_in, int* yalt_ovrf_in, int* zalt_ovrf_in,
                    llint* vxalt_in, llint* vyalt_in, llint* vzalt_in, int* vxalt_ovrf_in,
                    int* vyalt_ovrf_in, int* vzalt_ovrf_in, llint* fxalt_in, llint* fyalt_in,
                    llint* fzalt_in, int* fxalt_ovrf_in, int* fyalt_ovrf_in, int* fzalt_ovrf_in);

  /// \brief Copy and move constructors--as with any object containing const members, the move
  ///        assignment operator is implicitly deleted.
  /// \{
  PsSynthesisWriter(const PsSynthesisWriter &original) = default;
  PsSynthesisWriter(PsSynthesisWriter &&other) = default;
  /// \}
  
  // Even in the writer, some information may not be altered.  This includes the simulation
  // conditions, the number of systems, and the number of atoms in each system.  The layout
  // of the data must not be corrupted.
  const int system_count;               ///< The number of independent systems
  const int unique_topology_count;      ///< The number of unique topologies
  const UnitCellType unit_cell;         ///< Type of unit cells (or none) each system resides in
  const int* atom_starts;               ///< Points at which each system starts in the atom list
  const int* atom_counts;               ///< Atom counts for all systems
  const int* common_ag_list;            ///< Concatenated lists of systems in the synthesis sharing
                                        ///<   the same topology
  const int* common_ag_bounds;          ///< Bounds array for common_ag_list
  const int* unique_ag_idx;             ///< Indices of the unique topologies referenced by each
                                        ///<   set of coordinates
  const int* replica_idx;               ///< The instance number of each system in the list of
                                        ///<   all systems sharing its topology.
  
  // Scaling factors: the PhaseSpaceSynthesis permits a customizable discretization of fixed-point
  // arithmetic.
  const double gpos_scale;       ///< Global position coordinate scaling factor
  const double lpos_scale;       ///< Local position coordinate scaling factor
  const double vel_scale;        ///< Velocity coordinate scaling factor
  const double frc_scale;        ///< Scaling factor for fixed-precision force accumulation
  const double inv_gpos_scale;   ///< Inverse global coordinate scaling factor
  const double inv_lpos_scale;   ///< Inverse local coordinate scaling factor
  const double inv_vel_scale;    ///< Inverse velocity scaling factor
  const double inv_frc_scale;    ///< Inverse force scaling factor
  const float gpos_scale_f;      ///< Global position coordinate scaling factor
  const float lpos_scale_f;      ///< Local position coordinate scaling factor
  const float vel_scale_f;       ///< Velocity coordinate scaling factor
  const float frc_scale_f;       ///< Scaling factor for fixed-precision force accumulation
  const float inv_gpos_scale_f;  ///< Inverse global coordinate scaling factor
  const float inv_lpos_scale_f;  ///< Inverse local coordinate scaling factor
  const float inv_vel_scale_f;   ///< Inverse velocity scaling factor
  const float inv_frc_scale_f;   ///< Inverse force scaling factor
  const int gpos_bits;           ///< Global position coordinate bits after the decimal
  const int lpos_bits;           ///< Local position coordinate bits after the decimal
  const int vel_bits;            ///< Velocity coordinate bits after the decimal
  const int frc_bits;            ///< Force component bits after the decimal

  // Pointers to the transformations and box vectors are mutable if the systems change volume.
  llint* boxvecs;        ///< Discretized box vectors
  int* boxvec_ovrf;      ///< Overflow arrays for the discretized box vectors
  double* umat;          ///< Box (fractional) space transformation matrices, one per warp
  double* invu;          ///< Inverse transformation matrices, one per warp
  double* boxdims;       ///< Box dimensions (a, b, c, alpha, beta, gamma)
  llint* alt_boxvecs;    ///< Discretized box vectors
  int* alt_boxvec_ovrf;  ///< Overflow arrays for the discretized box vectors
  double* umat_alt;      ///< Box (fractional) space transformation matrices, one per warp
  double* invu_alt;      ///< Inverse transformation matrices, one per warp
  double* alt_boxdims;   ///< Box dimensions (a, b, c, alpha, beta, gamma)

  // Pointers to the coordinate, velocity, and force data--these are mutable for accumulating
  // forces and letting a trajectory evolve.
  llint* xcrd;      ///< Non-wrapped Cartesian X coordinates of all particles
  llint* ycrd;      ///< Non-wrapped Cartesian Y coordinates of all particles
  llint* zcrd;      ///< Non-wrapped Cartesian Z coordinates of all particles
  int* xcrd_ovrf;   ///< Non-wrapped Cartesian X coordinate overflow buffers
  int* ycrd_ovrf;   ///< Non-wrapped Cartesian Y coordinate overflow buffers
  int* zcrd_ovrf;   ///< Non-wrapped Cartesian Z coordinate overflow buffers
  llint* xvel;      ///< Cartesian X velocities
  llint* yvel;      ///< Cartesian Y velocities
  llint* zvel;      ///< Cartesian Z velocities
  int* xvel_ovrf;   ///< Cartesian X velocity overflow buffers
  int* yvel_ovrf;   ///< Cartesian Y velocity overflow buffers
  int* zvel_ovrf;   ///< Cartesian Z velocity overflow buffers
  llint* xfrc;      ///< Discretized Cartesian X forces
  llint* yfrc;      ///< Discretized Cartesian Y forces
  llint* zfrc;      ///< Discretized Cartesian Z forces
  int* xfrc_ovrf;   ///< Discretized Cartesian X force overflow buffers
  int* yfrc_ovrf;   ///< Discretized Cartesian Y force overflow buffers
  int* zfrc_ovrf;   ///< Discretized Cartesian Z force overflow buffers
  llint* xalt;      ///< Alternate Cartesian X positions of particles
  llint* yalt;      ///< Alternate Cartesian Y positions of particles
  llint* zalt;      ///< Alternate Cartesian Z positions of particles
  int* xalt_ovrf;   ///< Overflow buffers for particles' alternate Cartesian X locations
  int* yalt_ovrf;   ///< Overflow buffers for particles' alternate Cartesian Y locations
  int* zalt_ovrf;   ///< Overflow buffers for particles' alternate Cartesian Z locations
  llint* vxalt;     ///< Alternate Cartesian X velocities of particles
  llint* vyalt;     ///< Alternate Cartesian Y velocities of particles
  llint* vzalt;     ///< Alternate Cartesian Z velocities of particles
  int* vxalt_ovrf;  ///< Overflow buffers for particles' alternate Cartesian X velocities
  int* vyalt_ovrf;  ///< Overflow buffers for particles' alternate Cartesian Y velocities
  int* vzalt_ovrf;  ///< Overflow buffers for particles' alternate Cartesian Z velocities
  llint* fxalt;     ///< Alternate Cartesian X forces acting on particles
  llint* fyalt;     ///< Alternate Cartesian Y forces acting on particles
  llint* fzalt;     ///< Alternate Cartesian Z forces acting on particles
  int* fxalt_ovrf;  ///< Overflow buffers for particles' alternate Cartesian X forces
  int* fyalt_ovrf;  ///< Overflow buffers for particles' alternate Cartesian Y forces
  int* fzalt_ovrf;  ///< Overflow buffers for particles' alternate Cartesian Z force
};

/// \brief The reader for a PhaseSpaceSynthesis object, containing all of the data relevant for
///        propagating dynamics in a collection of systems.
struct PsSynthesisReader {

  /// \brief The constructor can take lists of pointers and constants from the parent object
  ///        (in this case, a PhaseSpaceSynthesis), like other abstracts.  Like the
  ///        CoordinateFrameReader, it can also convert its cognate writer by turning all of
  ///        the relevant pointers const.
  /// \{
  PsSynthesisReader(int system_count_in, int unique_topology_count_in, UnitCellType unit_cell_in,
                    const int* atom_starts_in, const int* atom_counts_in,
                    const int* common_ag_list_in, const int* common_ag_bounds_in,
                    const int* unique_ag_idx_in, const int* replica_idx_in, double gpos_scale_in,
                    double lpos_scale_in, double vel_scale_in, double frc_scale_in,
                    int gpos_bits_in, int lpos_bits_in, int vel_bits_in, int frc_bits_in,
                    const llint* boxvecs_in, const int* boxvec_ovrf_in, const double* umat_in,
                    const double* invu_in, const double* boxdims_in, const llint* alt_boxvecs_in,
                    const int* alt_boxvec_ovrf_in, const double* umat_alt_in,
                    const double* invu_alt_in, const double* alt_boxdims_in, const llint* xcrd_in,
                    const llint* ycrd_in, const llint* zcrd_in, const int* xcrd_ovrf_in,
                    const int* ycrd_ovrf_in, const int* zcrd_ovrf_in, const llint* xvel_in,
                    const llint* yvel_in, const llint* zvel_in, const int* xvel_ovrf_in,
                    const int* yvel_ovrf_in, const int* zvel_ovrf_in, const llint* xfrc_in,
                    const llint* yfrc_in, const llint* zfrc_in, const int* xfrc_ovrf_in,
                    const int* yfrc_ovrf_in, const int* zfrc_ovrf_in, const llint* xalt_in,
                    const llint* yalt_in, const llint* zalt_in, const int* xalt_ovrf_in,
                    const int* yalt_ovrf_in, const int* zalt_ovrf_in, const llint* vxalt_in,
                    const llint* vyalt_in, const llint* vzalt_in, const int* vxalt_ovrf_in,
                    const int* vyalt_ovrf_in, const int* vzalt_ovrf_in, const llint* fxalt_in,
                    const llint* fyalt_in, const llint* fzalt_in, const int* fxalt_ovrf_in,
                    const int* fyalt_ovrf_in, const int* fzalt_ovrf_in);

  PsSynthesisReader(const PsSynthesisWriter &psyw);
  /// \}

  /// \brief Copy and move constructors--as with any object containing const members, the move
  ///        assignment operator is implicitly deleted.
  /// \{
  PsSynthesisReader(const PsSynthesisReader &original) = default;
  PsSynthesisReader(PsSynthesisReader &&other) = default;
  /// \}
  
  // System sizing information
  const int system_count;               ///< The number of independent systems
  const int unique_topology_count;      ///< The number of unique topologies
  const UnitCellType unit_cell;         ///< Type of unit cells (or none) each system resides in
  const int* atom_starts;               ///< Points at which each system starts in the atom list
  const int* atom_counts;               ///< Atom counts for all systems
  const int* common_ag_list;            ///< Concatenated lists of systems in the synthesis sharing
                                        ///<   the same topology
  const int* common_ag_bounds;          ///< Bounds array for common_ag_list
  const int* unique_ag_idx;             ///< Indices of the unique topologies referenced by each
                                        ///<   set of coordinates
  const int* replica_idx;               ///< The instance number of each system in the list of
                                        ///<   all systems sharing its topology.

  // Scaling factors: the PhaseSpaceSynthesis permits a customizable discretization of fixed-point
  // arithmetic.
  const double gpos_scale;       ///< Global position coordinate scaling factor
  const double lpos_scale;       ///< Local position coordinate scaling factor
  const double vel_scale;        ///< Velocity coordinate scaling factor
  const double frc_scale;        ///< Scaling factor for fixed-precision force accumulation
  const double inv_gpos_scale;   ///< Inverse global coordinate scaling factor
  const double inv_lpos_scale;   ///< Inverse local coordinate scaling factor
  const double inv_vel_scale;    ///< Inverse velocity scaling factor
  const double inv_frc_scale;    ///< Inverse force scaling factor
  const float gpos_scale_f;      ///< Global position coordinate scaling factor
  const float lpos_scale_f;      ///< Local position coordinate scaling factor
  const float vel_scale_f;       ///< Velocity coordinate scaling factor
  const float frc_scale_f;       ///< Scaling factor for fixed-precision force accumulation
  const float inv_gpos_scale_f;  ///< Inverse global coordinate scaling factor
  const float inv_lpos_scale_f;  ///< Inverse local coordinate scaling factor
  const float inv_vel_scale_f;   ///< Inverse velocity scaling factor
  const float inv_frc_scale_f;   ///< Inverse force scaling factor
  const int gpos_bits;           ///< Global position coordinate bits after the decimal
  const int lpos_bits;           ///< Local position coordinate bits after the decimal
  const int vel_bits;            ///< Velocity coordinate bits after the decimal
  const int frc_bits;            ///< Force component bits after the decimal
  
  // Pointers to the transformations and box vectors are likewise const--once created, this
  // object is valid for a system held in constant volume.
  const llint* boxvecs;        ///< Discretized box vectors
  const int* boxvec_ovrf;      ///< Overflow arrays for the discretized box vectors
  const double* umat;          ///< Box (fractional) space transformation matrices, one per warp
  const double* invu;          ///< Inverse transformation matrices, one per warp
  const double* boxdims;       ///< Box dimensions (a, b, c, alpha, beta, gamma)
  const llint* alt_boxvecs;    ///< Discretized box vectors
  const int* alt_boxvec_ovrf;  ///< Overflow arrays for the discretized box vectors
  const double* umat_alt;      ///< Box (fractional) space transformation matrices, one per warp
  const double* invu_alt;      ///< Inverse transformation matrices, one per warp
  const double* alt_boxdims;   ///< Box dimensions (a, b, c, alpha, beta, gamma)

  // Pointers to the coordinate, velocity, and force data--these are mutable for accumulating
  // forces and letting a trajectory evolve.
  const llint* xcrd;      ///< Non-wrapped Cartesian X coordinates of all particles
  const llint* ycrd;      ///< Non-wrapped Cartesian Y coordinates of all particles
  const llint* zcrd;      ///< Non-wrapped Cartesian Z coordinates of all particles
  const int* xcrd_ovrf;   ///< Non-wrapped Cartesian X coordinate overflow buffers
  const int* ycrd_ovrf;   ///< Non-wrapped Cartesian Y coordinate overflow buffers
  const int* zcrd_ovrf;   ///< Non-wrapped Cartesian Z coordinate overflow buffers
  const llint* xvel;      ///< Cartesian X velocities
  const llint* yvel;      ///< Cartesian Y velocities
  const llint* zvel;      ///< Cartesian Z velocities
  const int* xvel_ovrf;   ///< Cartesian X velocity overflow buffers
  const int* yvel_ovrf;   ///< Cartesian Y velocity overflow buffers
  const int* zvel_ovrf;   ///< Cartesian Z velocity overflow buffers
  const llint* xfrc;      ///< Discretized Cartesian X forces
  const llint* yfrc;      ///< Discretized Cartesian Y forces
  const llint* zfrc;      ///< Discretized Cartesian Z forces
  const int* xfrc_ovrf;   ///< Discretized Cartesian X force overflow buffers
  const int* yfrc_ovrf;   ///< Discretized Cartesian Y force overflow buffers
  const int* zfrc_ovrf;   ///< Discretized Cartesian Z force overflow buffers
  const llint* xalt;      ///< Alternate Cartesian X positions of particles
  const llint* yalt;      ///< Alternate Cartesian Y positions of particles
  const llint* zalt;      ///< Alternate Cartesian Z positions of particles
  const int* xalt_ovrf;   ///< Overflow buffers for particles' alternate Cartesian X locations
  const int* yalt_ovrf;   ///< Overflow buffers for particles' alternate Cartesian Y locations
  const int* zalt_ovrf;   ///< Overflow buffers for particles' alternate Cartesian Z locations
  const llint* vxalt;     ///< Alternate Cartesian X velocities of particles
  const llint* vyalt;     ///< Alternate Cartesian Y velocities of particles
  const llint* vzalt;     ///< Alternate Cartesian Z velocities of particles
  const int* vxalt_ovrf;  ///< Overflow buffers for particles' alternate Cartesian X velocities
  const int* vyalt_ovrf;  ///< Overflow buffers for particles' alternate Cartesian Y velocities
  const int* vzalt_ovrf;  ///< Overflow buffers for particles' alternate Cartesian Z velocities
  const llint* fxalt;     ///< Alternate Cartesian X forces acting on all particles
  const llint* fyalt;     ///< Alternate Cartesian Y forces acting on all particles
  const llint* fzalt;     ///< Alternate Cartesian Z forces acting on all particles
  const int* fxalt_ovrf;  ///< Overflow buffers for particles' alternate Cartesian X forces
  const int* fyalt_ovrf;  ///< Overflow buffers for particles' alternate Cartesian Y forces
  const int* fzalt_ovrf;  ///< Overflow buffers for particles' alternate Cartesian Z forces
};

/// \brief A fixed-precision representation of coordinates, velocities, and forces to manage a set
///        of simulations.  The time steps are stored in units of femtoseconds.  Coordinates are
///        stored as long long integers to utilize the 32-bit int pipeline for difference
///        computations, bypassing fp64 computations wherever possible to get the best throughput
///        on many visualization, AI-oriented, and gaming cards.
class PhaseSpaceSynthesis {
public:

  /// \brief The constructor works from a series PhaseSpace object, importing  The LabFrame is not
  ///        a POINTER-kind object of any sort.
  ///
  /// Overloaded:
  ///   - Take arrays of PhaseSpace objects and AtomGraph pointers
  ///   - Take a SystemCache object and unpack it
  ///   - Skip the specification of thermostats, barostats, and the time step (just use defaults),
  ///     but explicitly specify at least the global position scaling bit count and possibly other
  ///     bit counts
  ///
  /// \param ps_list       Array of input coordinates, velocities, and forces objects
  /// \param ag_list       Array of pointers to input topologies
  /// \param index_key     Indices of the given topology and coordinate objects to assemble into a
  ///                      larger list of systems to be held within the resulting
  ///                      PhaseSpaceSynthesis object.
  /// \param ps_index_key  Indices of the given coordinate objects to assemble into a larger list
  ///                      of systems to be held within the resulting PhaseSpaceSynthesis object.
  /// \param ag_index_key  Indices of the given topology objects to assemble into a larger list of
  ///                      systems to be held within the resulting PhaseSpaceSynthesis object.
  /// \{
  PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list,
                      const std::vector<AtomGraph*> &ag_list,
                      int globalpos_scale_bits_in = default_globalpos_scale_bits,
                      int localpos_scale_bits_in = default_localpos_scale_bits,
                      int velocity_scale_bits_in = default_velocity_scale_bits,
                      int force_scale_bits_in = default_force_scale_bits);

  PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list,
                      const std::vector<AtomGraph*> &ag_list,
                      const std::vector<int> &index_key,
                      int globalpos_scale_bits_in = default_globalpos_scale_bits,
                      int localpos_scale_bits_in = default_localpos_scale_bits,
                      int velocity_scale_bits_in = default_velocity_scale_bits,
                      int force_scale_bits_in = default_force_scale_bits);

  PhaseSpaceSynthesis(const std::vector<PhaseSpace> &ps_list, const std::vector<int> &ps_index_key,
                      const std::vector<AtomGraph*> &ag_list, const std::vector<int> &ag_index_key,
                      int globalpos_scale_bits_in = default_globalpos_scale_bits,
                      int localpos_scale_bits_in = default_localpos_scale_bits,
                      int velocity_scale_bits_in = default_velocity_scale_bits,
                      int force_scale_bits_in = default_force_scale_bits);
  /// \}

  /// \brief Copy and move constructors work much like their counterparts in the smaller
  ///        PhaseSpace object.
  ///
  /// \param original  The original PhaseSpaceSynthesis object
  /// \{
  PhaseSpaceSynthesis(const PhaseSpaceSynthesis &original);
  PhaseSpaceSynthesis(PhaseSpaceSynthesis &&original);
  /// \}

  /// \brief Copy and move assignment operators must also be explicitly written out.
  ///
  /// \param other  The other PhaseSpaceSynthesis object
  /// \{
  PhaseSpaceSynthesis& operator=(const PhaseSpaceSynthesis &other);
  PhaseSpaceSynthesis& operator=(PhaseSpaceSynthesis &&other);
  /// \}

  /// \brief Get the number of systems in the object.
  int getSystemCount() const;
  
  /// \brief Get the number of unique topologies in the object.
  int getUniqueTopologyCount() const;
  
  /// \brief Get the unit cell type.
  UnitCellType getUnitCellType() const;

  /// \brief Get the current position in the coordinate cycle.
  CoordinateCycle getCyclePosition() const;
  
  /// \brief Get the global position scaling bit count
  int getGlobalPositionBits() const;

  /// \brief Get the local position scaling bit count
  int getLocalPositionBits() const;

  /// \brief Get the velocity scaling bit count
  int getVelocityBits() const;

  /// \brief Get the force accumulation bit count
  int getForceAccumulationBits() const;

  /// \brief Get the topology pointer for a particular system within the synthesis.
  ///
  /// Overloaded:
  ///   - Get the topology pointer for a specific system
  ///   - Get topology pointers for all systems
  ///
  /// \param system_index  Index of the system of interest within the synthesis
  /// \{
  const AtomGraph* getSystemTopologyPointer(int system_index) const;
  const std::vector<AtomGraph*>& getSystemTopologyPointer() const;
  /// \}
  
  /// \brief Get a list of unique topology pointers from this coordinate synthesis, const-casted
  ///        for accessibility to other functions.
  const std::vector<AtomGraph*>& getUniqueTopologies() const;

  /// \brief Get the starting index of atoms for one of the systems, using its index.
  ///
  /// \param system_index  Index of the system of interest
  int getAtomOffset(int system_index) const;

  /// \brief Get the number of atoms in one of the systems, using its index.
  ///
  /// \param system_index  Index of the system of interest
  int getAtomCount(int system_index) const;

  /// \brief Get the total number of atoms across all systems in the synthesis, including padding
  ///        between systems and after the end of the final system.
  int getPaddedAtomCount() const;

  /// \brief Get the index of the unique topology used by a particular system in the synthesis.
  ///
  /// \param system_index  Index of the system of interest
  int getUniqueTopologyIndex(int system_index) const;

  /// \brief Get the unique topology index of all systems in the synthesis.
  std::vector<int> getUniqueTopologyIndices() const;
  
  /// \brief Get a list of system indices from within this coordinate synthesis, providing examples
  ///        of each unique topology as presented in the order the topology pointers appear in the
  ///        output of getUniqueTopologies() above.
  std::vector<int> getUniqueTopologyExampleIndices() const;

  /// \brief Get the number of systems that make use of a particular topology.
  ///
  /// \param topology_index  Index of the unique topology of interest
  int getTopologyInstanceCount(int topology_index) const;
  
  /// \brief Get a list of system indices which all reference the same unique topology (the unique
  ///        topology index is determined by the order in which each unique topology appears in the
  ///        overall list of systems in the object).
  ///
  /// \param topology_index  Index of the unique topology of interest
  std::vector<int> getSystemIndicesByTopology(int topology_index) const;
  
  /// \brief Get a const pointer to the object itself.
  const PhaseSpaceSynthesis* getSelfPointer() const;

  /// \brief Get the reader or writer, as appropriate based on the const-ness of this object.
  ///
  /// Overloaded:
  ///   - Get a read-only abstract for a const object
  ///   - Get a writeable abstract for a mutable object
  ///   - Choose the stage in the time cycle that shall be deemed the relevant coordinates
  ///
  /// \param tier         The level (host or device) at which to get the set of pointers
  /// \param orientation  Stage in the object's time cycle to take as the current coordinates
  /// \{
  const PsSynthesisReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  const PsSynthesisReader data(CoordinateCycle orientation,
                               HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  PsSynthesisWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  PsSynthesisWriter data(CoordinateCycle orientation,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Get the read-only summary of the system sizes.  The data is never needed in as a
  ///        collection of device viewable pointers to host-side data, as this information is
  ///        constant as of the creation of the object and therefore consistent on both the CPU
  ///        host and GPU device.
  ///
  /// \param tier  The level (host or device) at which to get the set of pointers
  const PsSynthesisBorders borders(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

#ifdef STORMM_USE_HPC
  /// \brief Get a special writer which allows the device to read and write to host-mapped data.
  ///        This form of the writer can be used in kernel calls that streamline download and
  ///        upload of specific systems within the PhaseSpaceSynthesis with the launch latency
  ///        of only a single kernel call.  This offers a huge reduction in bandwidth requirements
  ///        and lower latency than multiple cudaMemcpy calls.
  ///
  /// Overloaded:
  ///   - Let the object automatically assign pointers to past, present, and future coordinates
  ///   - Specify which stage of the time cycle is to be taken as current
  ///   - Produce a read-only abstract for a const synthesis or a writeable abstract for a mutable
  ///     synthesis
  ///
  /// \param orientation  Optional argument stipulating from which stage of the cycle to take the
  ///                     current coordinates, past coordinates, and pointers to any future
  ///                     coordinates.  If not specified the object's own cycle position will be
  ///                     used.
  /// \{
  const PsSynthesisReader deviceViewToHostData(CoordinateCycle orientation) const;
  const PsSynthesisReader deviceViewToHostData() const;
  PsSynthesisWriter deviceViewToHostData(CoordinateCycle orientation);
  PsSynthesisWriter deviceViewToHostData();
  /// \}
  
  /// \brief Upload data to the device.
  ///
  /// Overloaded:
  ///   - Upload all data
  ///   - Upload data associated with a specific type of coordinates
  ///   - Upload data assoicated with a specific type of coordinates and specific systems
  ///
  /// \param kind                Choose from POSITIONS, VELOCITIES, FORCES (positional data will
  ///                            include box dimensions and is the only way to get box dimensions)
  /// \param system_index        The specific system to upload
  /// \param system_lower_bound  Upload systems in the range [ lower_bound, upper_bound )
  /// \param system_upper_bound  Upload systems in the range [ lower_bound, upper_bound )
  /// \{
  void upload();
  void upload(TrajectoryKind kind);
  void upload(TrajectoryKind kind, int system_index, const GpuDetails &gpu);
  void upload(TrajectoryKind kind, int system_lower_bound, int system_upper_bound,
              const GpuDetails &gpu);
  /// \}
  
  /// \brief Download data from the device
  ///
  /// Overloaded:
  ///   - Download all data
  ///   - Download data associated with a specific type of coordinates
  ///   - Download data assoicated with a specific type of coordinates and specific systems
  ///
  /// \param kind                Choose from POSITIONS, VELOCITIES, FORCES (positional data will
  ///                            include box dimensions and is the only way to get box dimensions)
  /// \param system_index        The specific system to download
  /// \param system_lower_bound  Download systems in the range [ lower_bound, upper_bound )
  /// \param system_upper_bound  Download systems in the range [ lower_bound, upper_bound )
  /// \{
  void download();
  void download(TrajectoryKind kind);
  void download(TrajectoryKind kind, int system_index, const GpuDetails &gpu);
  void download(TrajectoryKind kind, int system_lower_bound, int system_upper_bound,
                const GpuDetails &gpu);
  /// \}
#endif
  
  /// \brief Extract the phase space (plus forces) of a specific system within the synthesis.
  ///
  /// \param ps           Pointer to an allocated PhaseSpace object (i.e. the original) ready to
  ///                     accept data from the synthesis (which may have evolved since it was
  ///                     first loaded)
  /// \param trajkind     Type of coordinates to copy
  /// \param index        Index of the system of interest within the synthesis
  /// \param origin       The level (host or device) at which to get the data
  /// \param destination  The level (host or device) at which to get the data
  /// \param gpu          Details of the GPU in use
  /// \{
  void extractSystem(PhaseSpace *ps, int index,
                     HybridTargetLevel origin = HybridTargetLevel::HOST,
                     HybridTargetLevel destination = HybridTargetLevel::HOST,
                     const GpuDetails &gpu = null_gpu) const;
  /// \}
  
  /// \brief Export a system's coordinates, velocities, and forces to a PhaseSpace object.
  ///
  /// \param index     Index of the system of interest within the synthesis
  /// \param tier      The level (host or device) at which to get the data
  PhaseSpace exportSystem(int index, HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
  /// \brief Export a system's coordinates, velocities, or forces to a compact CoordinateFrame
  ///        object.
  ///
  /// Overloaded:
  ///   - Provide the stage of the coordinate cycle from which to obtain coordinates
  ///   - Assume that the WHITE stage of the coordinate cycle is desired
  ///
  /// \param index        Index of the system of interest within the synthesis
  /// \param trajkind     Type of coordinates to copy
  /// \param orientation  Stage of the coordinate cycle from which to obtain coordinates
  /// \param tier         The level (host or device) at which to get the data
  /// \{
  CoordinateFrame exportCoordinates(int index, CoordinateCycle orientation,
                                    TrajectoryKind trajkind = TrajectoryKind::POSITIONS,
                                    HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  CoordinateFrame exportCoordinates(int index, TrajectoryKind trajkind = TrajectoryKind::POSITIONS,
                                    HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Get the interlaced X, Y, and Z coordinates of one system in particular.  The result
  ///        is comparable to the getInterlacedCoordinates() functions held by CoordinateFrame and
  ///        PhaseSpace objects.  The result will be provided in internal units of Angstroms (for
  ///        positions), Angstroms per femtosecond (for velocities), or kcal/mol-A (for forces).
  ///
  /// Overloaded:
  ///   - Provide the stage of the coordinate cycle from which to obtain coordinates
  ///   - Assume that the WHITE stage of the coordinate cycle is desired
  ///
  /// \param index        Index of the system of interest within the synthesis
  /// \param trajkind     Type of coordinates to copy
  /// \param orientation  Stage of the coordinate cycle from which to obtain coordinates
  /// \param tier         The level (host or device) at which to get the data
  /// \{
  std::vector<double>
  getInterlacedCoordinates(int index, CoordinateCycle orientation,
                           TrajectoryKind trajkind = TrajectoryKind::POSITIONS,
                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  std::vector<double>
  getInterlacedCoordinates(int index, TrajectoryKind trajkind = TrajectoryKind::POSITIONS,
                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Modify the object's cycle position, as in the course of propagating dynamics.
  ///
  /// Overloaded:
  ///   - Increment the cycle position one step forward (this toggles between the WHITE and
  ///     BLACK states, and is the same as decrementing the cycle position unless the time
  ///     cycle takes on a third possible stage)
  ///   - Set the time cycle to a specific point
  ///
  /// \param time_point  The point in the time cycle which the object is to take as "current"
  /// \{
  void updateCyclePosition();
  void updateCyclePosition(CoordinateCycle time_point);
  /// \}
  
  /// \brief Initialize the forces within a PhaseSpaceSynthesis object.  This is the analog of the
  ///        eponymous function in the PhaseSpace object.
  ///
  /// Overloaded:
  ///   - Initialize one or more systems' forces on the host if STORMM is compiled for CPU only
  ///   - Initialize one or more systems' forces on either the host or the device if STORMM is
  ///     compiled for HPC
  ///   - Initialize forces for a specific point in the time cycle, or the current point
  ///
  /// \param orientation  Point in the time cycle at which to clear / initialize the forces
  /// \param gpu          Details of the GPU in use
  /// \param tier         The level (host or device) at which to initialize forces
  /// \param index        Index of the system of interest within the synthesis--if negative, all
  ///                     systems will have their forces initialized
  /// \{
#ifdef STORMM_USE_HPC
  void initializeForces(CoordinateCycle orientation, const GpuDetails &gpu,
                        HybridTargetLevel tier = HybridTargetLevel::HOST, int index = -1);
  void initializeForces(const GpuDetails &gpu = null_gpu,
                        HybridTargetLevel tier = HybridTargetLevel::HOST, int index = -1);
#else
  void initializeForces(CoordinateCycle orientation, int index = -1);
  void initializeForces(int index = -1);
#endif
  /// \}

  /// \brief Set the alternate positions to the current forces, and initialize velocities to zero,
  ///        as part of the first step of conjugate gradient optimization.  This primes the system
  ///        so that the alternate coordinates and velocities arrays can hold the prior forces and
  ///        temporary conjugate gradient memory.
  ///
  /// \param gpu    Details of the GPU in use
  /// \param tier   The level (host or device) at which to initialize vectors
#ifdef STORMM_USE_HPC
  void primeConjugateGradientCalculation(const GpuDetails &gpu,
                                         HybridTargetLevel tier = HybridTargetLevel::HOST);
#else
  void primeConjugateGradientCalculation();
#endif
  
  /// \brief Print a list of structures to a trajectory file.  Download will be performed when
  ///        calling this function, over the subset of relevant frames and data.
  ///
  /// \param system_indices  List of system coordinates / velocities / forces to print
  /// \param file_name       Name of the file to write, or base name of a set of files to write
  /// \param current_time    Current time progress of the group of simulations
  /// \param output_kind     The type of trajectory file to write
  /// \param expectation     The state that the output trajectory file is expected to be found in
  void printTrajectory(const std::vector<int> &system_indices, const std::string &file_name,
                       double current_time, CoordinateFileKind output_kind,
                       PrintSituation expectation) const;

  /// \brief Import a system from one of the other coordinate objects, or from a series of C-style
  ///        arrays with trusted lengths.  The imported system's size must correspond to that
  ///        expected by the atom count of the system it will replace.
  ///
  /// Overloaded:
  ///   - Provide a PhaseSpace object (all coordinates, velocities, and forces of the input object
  ///     will be transferred from the WHITE stage of the time cycle in the PhaseSpace object,
  ///     into the specified stage of the time cycle in this synthesis, and other stages of the
  ///     time cycle will be transferred accordingly).
  ///   - Provide three arrays of Cartesian X, Y, and Z coordinates, a scaling factor if the
  ///     data type is fixed-precision integral, plus indications of whether the data is for
  ///     positions, velocities, or forces, and at what stage of the time cycle the data is to
  ///     enter the PhaseSpace object.
  ///   - Provide a CoordinateFrame or CoordinateSeries object with a frame number, plus
  ///     indications of whether the object truly contains positions, velocities, or forces, and
  ///     what stage of the time cycle the data is to enter the PhaseSpaceSynthesis.
  ///   - Provide an abstract of any of the major coordinate objects, plus the other information to
  ///     target the import to the correct system, HPC tier, and place in the time cycle.
  ///
  /// \param ps                 Complete phase space and time cycle data intended to replace one
  ///                           of the systems in the synthesis
  /// \param system_index       Index of the system within this synthesis that the imported
  ///                           coordinates shall replace
  /// \param orientation        Stage of the time cycle at which the WHITE stage of an input
  ///                           PhaseSpace object is to enter the synthesis, or at which the data
  ///                           in raw arrays, a CoordinateFrame, or a CoordinateSeries object is
  ///                           to enter the synthesis
  /// \param tier               The level (host or device) at which to perform the transfer
  /// \param x_import           Input Cartesian X coordinates (these could be positions,
  ///                           velocities, or forces)
  /// \param y_import           Input Cartesian Y coordinates
  /// \param z_import           Input Cartesian Z coordinates
  /// \param box_xform_in       Transformation matrix to take coordinates into fractional space
  ///                           (for positions only--provide nullptr for velocities or forces)
  /// \param inverse_xform_in   Transformation matrix to take coordinates back to real space.  The
  ///                           units of elements in this matrix are Angstroms.
  /// \param box_dimensions_in  Dimensions of the box (redundant with the information stored in
  ///                           either of the transformation matrices, but convenient and perhaps
  ///                           best able to preserve bitwise information to pass it directly).
  ///                           The units of this array are Angstroms.
  /// \param kind               Specifies whether the Cartesian X, Y, and Z data are positions,
  ///                           velocities, or forces
  /// \param scaling_factor     Scaling factor to take the input X, Y, and Z data into internal
  ///                           units of Angstroms, Angstroms per femtosecond, or kcal/mol-A^2
  /// \param cf                 Input coordinate frame object containing X, Y, and Z data as
  ///                           double-precision objects
  /// \param cs                 Input coordinate series object with X, Y, and Z data for many
  ///                           frames, one of which will be copied over.  This object contains
  ///                           its own scaling factor.
  /// \param frame_index        Index of a CoordinateSeries object to be transferred
  /// \{
  void import(const PhaseSpaceReader &psr, int system_index, CoordinateCycle orientation,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  void import(const PhaseSpaceReader &psr, int system_index,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  void import(const PhaseSpaceWriter &psw, int system_index, CoordinateCycle orientation,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  void import(const PhaseSpaceWriter &psw, int system_index,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  void import(const PhaseSpace &ps, int system_index, CoordinateCycle orientation,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  void import(const PhaseSpace &ps, int system_index,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  void import(const CoordinateFrameReader &cfr, int system_index, CoordinateCycle orientation,
              TrajectoryKind kind = TrajectoryKind::POSITIONS,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  void import(const CoordinateFrameReader &cfr, int system_index,
              TrajectoryKind kind = TrajectoryKind::POSITIONS,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  void import(const CoordinateFrameWriter &cfw, int system_index, CoordinateCycle orientation,
              TrajectoryKind kind = TrajectoryKind::POSITIONS,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  void import(const CoordinateFrameWriter &cfw, int system_index,
              TrajectoryKind kind = TrajectoryKind::POSITIONS,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  void import(const CoordinateFrame &cf, int system_index, CoordinateCycle orientation,
              TrajectoryKind kind = TrajectoryKind::POSITIONS,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  void import(const CoordinateFrame &cf, int system_index,
              TrajectoryKind kind = TrajectoryKind::POSITIONS,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  template <typename T>
  void import(const T* x_import, const T* y_import, const T* z_import, const double* box_xform_in,
              const double* inverse_xform_in, const double* box_dimensions_in, int system_index,
              CoordinateCycle orientation, double inverse_scaling_factor = 1.0,
              TrajectoryKind kind = TrajectoryKind::POSITIONS,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  template <typename T>
  void import(const CoordinateSeriesReader<T> &csr, int frame_index, int system_index,
              CoordinateCycle orientation, TrajectoryKind kind = TrajectoryKind::POSITIONS,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  template <typename T>
  void import(const CoordinateSeriesReader<T> &csr, int frame_index, int system_index,
              TrajectoryKind kind = TrajectoryKind::POSITIONS,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  template <typename T>
  void import(const CoordinateSeriesWriter<T> &csw, int frame_index, int system_index,
              CoordinateCycle orientation, TrajectoryKind kind = TrajectoryKind::POSITIONS,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  template <typename T>
  void import(const CoordinateSeriesWriter<T> &csw, int frame_index, int system_index,
              TrajectoryKind kind = TrajectoryKind::POSITIONS,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  template <typename T>
  void import(const CoordinateSeries<T> &cs, int frame_index, int system_index,
              CoordinateCycle orientation, TrajectoryKind kind = TrajectoryKind::POSITIONS,
              HybridTargetLevel tier = HybridTargetLevel::HOST);

  template <typename T>
  void import(const CoordinateSeries<T> &cs, int frame_index, int system_index,
              TrajectoryKind kind = TrajectoryKind::POSITIONS,
              HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}
  
private:
  int system_count;               ///< The number of systems to tend at once
  int unique_topology_count;      ///< The number of unique topologies used by this synthesis
  UnitCellType unit_cell;         ///< The types of unit cells.  All unit cells must exist in
                                  ///<   isolated boundary conditions, or all cells must exist in
                                  ///<   periodic boundary conditions with a rectilinear or
                                  ///<   triclinic unit cell.
  CoordinateCycle cycle_position; ///< The stage of the time cycle (past >> present >> future) that
                                  ///<   holds the relevant, current coordinates

  // Scaling constants for fixed-precision coordinates (position, velocity) and also forces
  // in this PhaseSpaceSynthesis apply identically to all systems
  double globalpos_scale;         ///< Global position coordinate scaling factor
  double localpos_scale;          ///< Local position coordinate scaling factor
  double velocity_scale;          ///< Velocity coordinate scaling factor
  double force_scale;             ///< Scaling factor for fixed-precision force accumulation
  double inverse_globalpos_scale; ///< Inverse global coordinate scaling factor
  double inverse_localpos_scale;  ///< Inverse local coordinate scaling factor
  double inverse_velocity_scale;  ///< Inverse velocity scaling factor
  double inverse_force_scale;     ///< Inverse force scaling factor
  int globalpos_scale_bits;       ///< Global position coordinate bits after the decimal
  int localpos_scale_bits;        ///< Local position coordinate bits after the decimal
  int velocity_scale_bits;        ///< Velocity coordinate bits after the decimal
  int force_scale_bits;           ///< Force component bits after the decimal

  /// Starting positions for each system's stretch of atoms in xyz_qlj, (x,y,z)_velocities, and
  /// (x,y,z)_forces.  Atoms in each of those arrays will remain in their original orders, as
  /// presented in their respective topologies.
  Hybrid<int> atom_starts;

  /// The numbers of atoms in each system
  Hybrid<int> atom_counts;

  /// Lists of systems which make use of each unique topology.  The order of unique topologies is
  /// the same as appears in the unique_topologies array, returned by the getUniqueTopologies()
  /// function.
  Hybrid<int> shared_topology_instances;

  /// Bounds array for shared_topology_instances above
  Hybrid<int> shared_topology_instance_bounds;

  /// Indices of the unique topologies used by each system (this enables the reverse lookup of
  /// shared_topology_instances above).
  Hybrid<int> unique_topology_reference;

  /// Indices indicating the number of the replica of each system in the synthesis.  If a system
  /// is the first to make use of a particular topology, it has a value of 0 in this array.  The
  /// fourth system to make use of the same topology would have a value of 3.
  Hybrid<int> shared_topology_instance_index;
  
  // These variables are POINTER-kind Hybrid objects targeting the llint_data and int_data arrays
  Hybrid<llint> x_coordinates;          ///< Cartesian X coordinates of all particles
  Hybrid<llint> y_coordinates;          ///< Cartesian Y coordinates of all particles
  Hybrid<llint> z_coordinates;          ///< Cartesian Z coordinates of all particles
  Hybrid<int> x_coordinate_overflow;    ///< X direction coordinate overflow buffers
  Hybrid<int> y_coordinate_overflow;    ///< Y direction coordinate overflow buffers
  Hybrid<int> z_coordinate_overflow;    ///< Z direction coordinate overflow buffers
  Hybrid<llint> x_alt_coordinates;      ///< Alternate Cartesian X coordinates of all particles
  Hybrid<llint> y_alt_coordinates;      ///< Alternate Cartesian Y coordinates of all particles
  Hybrid<llint> z_alt_coordinates;      ///< Alternate Cartesian Z coordinates of all particles
  Hybrid<int> x_alt_coord_overflow;     ///< Overflow buffers for alternate Cartesian X coordinates
  Hybrid<int> y_alt_coord_overflow;     ///< Overflow buffers for alternate Cartesian Y coordinates
  Hybrid<int> z_alt_coord_overflow;     ///< Overflow buffers for alternate Cartesian Z coordinates
  Hybrid<llint> x_velocities;           ///< Cartesian X velocities of all particles
  Hybrid<llint> y_velocities;           ///< Cartesian Y velocities of all particles
  Hybrid<llint> z_velocities;           ///< Cartesian Z velocities of all particles
  Hybrid<int> x_velocity_overflow;      ///< Overflow buffers for velocities in the X direction
  Hybrid<int> y_velocity_overflow;      ///< Overflow buffers for velocities in the Y direction
  Hybrid<int> z_velocity_overflow;      ///< Overflow buffers for velocities in the Z direction
  Hybrid<llint> x_alt_velocities;       ///< Alternate Cartesian X velocities of all particles
  Hybrid<llint> y_alt_velocities;       ///< Alternate Cartesian Y velocities of all particles
  Hybrid<llint> z_alt_velocities;       ///< Alternate Cartesian Z velocities of all particles
  Hybrid<int> x_alt_velocity_overflow;  ///< Overflow buffers for alternate velocities in X
  Hybrid<int> y_alt_velocity_overflow;  ///< Overflow buffers for alternate velocities in Y
  Hybrid<int> z_alt_velocity_overflow;  ///< Overflow buffers for alternate velocities in Z
  Hybrid<llint> x_forces;               ///< Forces acting on particles in the X direction
  Hybrid<llint> y_forces;               ///< Forces acting on particles in the Y direction
  Hybrid<llint> z_forces;               ///< Forces acting on particles in the Z direction
  Hybrid<int> x_force_overflow;         ///< Force overflows acting on particles in X
  Hybrid<int> y_force_overflow;         ///< Force overflows acting on particles in Y
  Hybrid<int> z_force_overflow;         ///< Force overflows acting on particles in Z
  Hybrid<llint> x_alt_forces;           ///< Alternate accumulators for forces acting on particles
                                        ///<   in the X direction
  Hybrid<llint> y_alt_forces;           ///< Alternate forces acting in the Y direction
  Hybrid<llint> z_alt_forces;           ///< Alternate forces acting in the Z direction
  Hybrid<int> x_alt_force_overflow;     ///< Overflows for alternate forces acting in X
  Hybrid<int> y_alt_force_overflow;     ///< Overflows for alternate forces acting in Y
  Hybrid<int> z_alt_force_overflow;     ///< Overflows for alternate forces acting in Z
  Hybrid<llint> box_vectors;            ///< Scaled real space transformation matrix--moving
                                        ///<   particles between images by adding or subtracting
                                        ///<   multiples of these vectors can be expeditious and
                                        ///<   keeps coordinate representations consistent
                                        ///<   between the lab frame and the primary unit cell.
  Hybrid<int> box_vector_overflow;      ///< Overflow arrays for the discretized box vectors, to
                                        ///<   let their precision match that of the coordinate
                                        ///<   arrays (up to 72 bits after the decimal, see
                                        ///<   Constants/fixed_precision.h)
  Hybrid<llint> alt_box_vectors;        ///< Fixed-precision representations for the alternate
                                        ///<   coordinate representations' unit cells
  Hybrid<int> alt_box_vector_overflow;  ///< Overflow array for numbers in alt_box_vectors
  
  // The following are POINTER-kind Hybrid objects targeting the floating point data arrays
  Hybrid<double> box_space_transforms;    ///< Transformation matrices to take coordinates into
                                          ///<   box (fractional) space
  Hybrid<double> inverse_transforms;      ///< Transformation matrices to go back to real space
  Hybrid<double> box_dimensions;          ///< Three box lengths and the angles in the planes
                                          ///<   normal to each axis in each system (synchronized
                                          ///<   with the transformation matrices--each update of
                                          ///<   the box_dimensions triggers an update of the
                                          ///<   transformations)
  Hybrid<double> alt_box_transforms;      ///< Transformation matrices to take coordinates into
                                          ///<   box (fractional) space for the BLACK
                                          ///<   coordinate sets
  Hybrid<double> alt_inverse_transforms;  ///< Inverse transformation matrices for the BLACK
                                          ///<   coordinate sets
  Hybrid<double> alt_box_dimensions;      ///< Box lengths and the angles for the BLACK
                                          ///<   coordinate sets

  // Data arrays
  Hybrid<int> int_data;        ///< Counts of atoms and starting points for each system, plus
                               ///<   overflow data
  Hybrid<llint> llint_data;    ///< The discretized data for all of phase space and forces
  Hybrid<double> double_data;  ///< Double-precision floating point transformations--these are the
                               ///<   standard for moving coordinates into a re-imaged
                               ///<   configuration for later computations.  Once particles are
                               ///<   in a re-imaged configuration for a given box size, they can
                               ///<   be manipulated in parallel with the coordinates in the lab
                               ///<   frame.  The unit cell transformations are computed in double
                               ///<   precision, then used to construct the long long int
                               ///<   box_vectors as scaled representation, then synchronized to
                               ///<   the box_vectors to make consistent representations for
                               ///<   calculating excursions from the unit cell boundaries.

  /// Pointers to the topologies that describe each system
  std::vector<AtomGraph*> topologies;

  /// Pointers to the unique topologies used by this synthesis
  std::vector<AtomGraph*> unique_topologies;
  
  /// \brief Allocate private array data
  ///
  /// \param atom_stride  The total number of padded atoms in the system (sum over all individual
  ///                     systems with warp size padding in each of them)
  void allocate(size_t atom_stride);
  
#ifdef STORMM_USE_HPC
  /// \brief Extract a system into a pre-allocated PhaseSpace object based on information in this
  ///        PhaseSpaceSynthesis on the HPC device.
  ///
  /// \param psw     Abstract for coordinates, velocities, and forces in double-precision real
  ///                values, with an implicit indication of the level onto which the data shall be
  ///                extracted (based on whether the abstract was taken for the HOST or DEVICE)
  /// \param index   Index of the system to extract
  /// \param gpu     Details of the GPU in use
  /// \param origin  Level of the synthesis from which the data shall be extracted
  void extractSystem(PhaseSpaceWriter *psw, int index, const GpuDetails &gpu,
                     HybridTargetLevel origin = HybridTargetLevel::DEVICE) const;
#endif

  /// \brief Check whether the requested system index is valid.
  ///
  /// \param index   The numerical index of the system of interest within the synthesis
  /// \param caller  [Optional] Name of the calling function
  void validateSystemIndex(int index, const char* caller = nullptr) const;
};

/// \brief Define the type index for the PhaseSpaceSynthesis object.
static const size_t phasespace_synthesis_type_index =
  std::type_index(typeid(PhaseSpaceSynthesis)).hash_code();

} // namespace trajectory
} // namespace stormm

#include "phasespace_synthesis.tpp"

// As with common types and STORMM vector types, define the type indices for general use in the
// STORMM namespace.
namespace stormm {
using synthesis::phasespace_synthesis_type_index;
} // namespace stormm

#endif
