// -*-c++-*-
#ifndef STORMM_RMSD_PLAN_H
#define STORMM_RMSD_PLAN_H

#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/gpu_details.h"
#include "Chemistry/atom_equivalence.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_cache_map.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Synthesis/systemcache.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

using card::GpuDetails;
using card::Hybrid;
using card::HybridTargetLevel;
using chemistry::AtomEquivalence;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using synthesis::SynthesisCacheMap;
using synthesis::SystemCache;
using synthesis::SystemGrouping;
using topology::AtomGraph;
using topology::ChemicalDetailsKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateSeries;

/// \brief Assume that, if 75% of the mass (or coordinates) of a molecule are asymmetric or
///        otherwise placed in definie positions, the remainder of the molecule's symmetric atom
///        groups can be tested for best fit based on the pose that aligning the rest of the
///        molecule leaves them in.
constexpr double default_required_mass_fraction = 0.75;

/// \brief Read-only abstract for the RMSDPlan object, containing pointers to various buckets of
///        atom indices.
template <typename T> struct RMSDPlanReader {

  /// \brief The constructor takes inputs for all arguments.
  RMSDPlanReader(int plan_count_in, RMSDMethod strategy_in, double mass_fraction_in,
                 const T* masses_in, const int* atom_counts_in, const int* atom_starts_in,
                 const int* alignment_steps_in, const int* core_atoms_in,
                 const int* core_counts_in, const int* core_starts_in, const int* symm_atoms_in,
                 const int4* symm_bounds_in, const int2* symm_ranges_in);
  
  /// \brief The presence of const members will implicitly delete the copy and move assignment
  ///        operators, but the default copy and move constructors will apply.
  /// \{
  RMSDPlanReader(const RMSDPlanReader &original);
  RMSDPlanReader(RMSDPlanReader &&original);
  /// \}

  const int plan_count;        ///< The number of individual system RMSD calculation plans
  const RMSDMethod strategy;   ///< RMSD calculations are mass- or coordinate-weighted, aligned or
                               ///<   not
  const double mass_fraction;  ///< Total mass fraction needed to provide a critical alignment for
                               ///<   subsequent determination of symmetry-related atom groups
  const T* masses;             ///< Masses of particles in all systems
  const int* atom_counts;      ///< Numbers of atoms involved in the systems served by each plan
  const int* atom_starts;      ///< Starting positions of the atom series associated with each plan
                               ///<   (each plan pertains to a unique topology)
  const int* alignment_steps;  ///< Protocols for performing RMSD calculations under each plan
  const int* core_atoms;       ///< Concatenated lists of asymmetric atoms in each system
  const int* core_counts;      ///< Exact numbers of asymmetric atoms in each system (each system's
                               ///<   core atom list is padded by the warp size)
  const int* core_starts;      ///< Bounds array for core_atoms, above (use this for offsets, and
                               ///<   core counts for exact counts
  const int* symm_atoms;       ///< Concatenated lists of each system's symmetry-related atoms
  const int4* symm_bounds;     ///< Bounds for symmetry-related atom group lists in the symm_atoms
                               ///<   array.  The "x" and "y" members list the lower and upper
                               ///<   limits of all symmetry-related domains in each group, while
                               ///<   the "z" members list the domain size and the "w" member gives
                               ///<   the total number of domains.
  const int2* symm_ranges;     ///< Bounds array on the symmetric atom bounds array--indicating the
                               ///<   lower and upper limits of symmetric atom numbered groups in
                               ///<   each system in the "x" and "y" members of each tuple.
};
  
/// \brief Collect instructions for one or more systems (intend to work with any coordinate object,
///        including the PhaseSpaceSynthesis)
class RMSDPlan {
public:

  /// \brief The constructor can take any number of AtomEquivalence objects and build plans from
  ///        them.  Topology pointers will be harvested from these objects and, if necessary,
  ///        matched to topologies in a synthesis object in order to connect each system to the
  ///        appropriate instruction set.
  ///
  /// \param strategy_in   The type of RMSD calculations that this plan will guide
  /// \param rmf_in        The required mass fraction that constitutes enough of the molecule to
  ///                      perform placement of small symmetric atom groups without further
  ///                      alignment calculations
  /// \param ag_in         Topology for which to draw a positional RMSD computation plan
  /// \param cf_in         Input coordinates (used for making the necessary AtomEquivalence objects
  ///                      if they are not supplied explicitly)
  /// \param eq_in         Pre-computed breakdown of all symmetry-related atoms
  /// \param poly_ps_in    PhaseSpaceSynthesis object for which to draw plans covering all systems.
  ///                      Topology pointers will be harvested from the object.
  /// \param eq_list_in    List of symmetry-related atom breakdowns for all systems in poly_ps_in
  /// \param gpu           Details of the GPU that will be used in computing RMSD values
  /// \param low_mol_idx   Lower limit of molecules from the topology to assess for symmetry groups
  ///                      (this, and high_mol_idx, are only used if a single topology is provided
  ///                      for automatic deduction of symmetry groups)
  /// \param high_mol_idx  Upper limit of molecules from the topology to assess for symmetry
  ///                      groups.  The default of -1 indicates that only one molecule, indicated
  ///                      by the lower index, should be assessed.
  /// \{
  RMSDPlan(RMSDMethod strategy_in = RMSDMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction,
           const PhaseSpaceSynthesis *poly_ps_in = nullptr, const SystemCache *sc_in = nullptr,
           const SynthesisCacheMap *scmap_in = nullptr);

  RMSDPlan(const AtomGraph &ag_in, const CoordinateFrame &cf_in,
           RMSDMethod strategy_in = RMSDMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu,
           int low_mol_idx = 0, int high_mol_idx = -1);

  RMSDPlan(const AtomEquivalence &eq_in, RMSDMethod strategy_in = RMSDMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu);
  
  RMSDPlan(const PhaseSpaceSynthesis *poly_ps_in, RMSDMethod strategy_in = RMSDMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu,
           int low_mol_idx = 0, int high_mol_idx = -1);

  RMSDPlan(const PhaseSpaceSynthesis *poly_ps_in, const std::vector<AtomEquivalence> &eq_list_in,
           RMSDMethod strategy_in = RMSDMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu);

  RMSDPlan(const PhaseSpaceSynthesis &poly_ps_in, RMSDMethod strategy_in = RMSDMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu,
           int low_mol_idx = 0, int high_mol_idx = -1);

  RMSDPlan(const PhaseSpaceSynthesis &poly_ps_in, const std::vector<AtomEquivalence> &eq_list_in,
           RMSDMethod strategy_in = RMSDMethod::ALIGN_MASS,
           double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu);

  RMSDPlan(const PhaseSpaceSynthesis *poly_ps_in, const SystemCache *sc_in,
           const SynthesisCacheMap *scmap_in,
           const RMSDMethod strategy_in = RMSDMethod::ALIGN_MASS,
           const double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu,
           const int low_mol_idx = 0, const int high_mol_idx = -1);

  RMSDPlan(const PhaseSpaceSynthesis &poly_ps_in, const SystemCache &sc_in,
           const SynthesisCacheMap &scmap_in,
           const RMSDMethod strategy_in = RMSDMethod::ALIGN_MASS,
           const double rmf_in = default_required_mass_fraction, const GpuDetails &gpu = null_gpu,
           const int low_mol_idx = 0, const int high_mol_idx = -1);
  /// \}

  /// \brief With no POINTER-kind Hybrid or pointers to its own member variables to repair, and no
  ///        const members, the RMSDPlan can make use of default copy and move constructors as well
  ///        as copy and move assignment operators.
  ///
  /// \param original  Another object to copy or move in constructing this one
  /// \param other     The right-hand side object of an assignment
  /// \{
  RMSDPlan(const RMSDPlan &original) = default;
  RMSDPlan(RMSDPlan &&original) = default;
  RMSDPlan& operator=(const RMSDPlan &original) = default;
  RMSDPlan& operator=(RMSDPlan &&original) = default;
  /// \}
  
  /// \brief Get the number of plans kept within this object.
  int getPlanCount() const;

  /// \brief Get the plan index applicable to a specific grouping of systems (this assumes that all
  ///        systems in the group have a valid interpretation by a single plan, and if that is not
  ///        the case the result will be undefined).
  ///
  /// \param index         Index of the grouping from one of the internal lists
  /// \param organization  The manner in which systems are grouped
  int getPlanIndex(int index, SystemGrouping organization) const;
  
  /// \brief Get the general, prescribed strategy for computing RMSDs.
  RMSDMethod getGeneralStrategy() const;

  /// \brief Get the mass fraction required to test symmetry-related atom arrangements without
  ///        specific alignments.
  double getRequiredMassFraction() const;

  /// \brief Get one of the topology pointers used by the RMSD plan.
  ///
  /// \param plan_index  The system / plan of interest
  const AtomGraph* getTopologyPointer(int plan_index) const;

  /// \brief Get one of the alignment protocols.
  ///
  /// \param plan_index  The system / plan of interest
  RMSDAlignmentProtocol getAlignmentProtocol(int plan_index) const;

  /// \brief Get the read-only abstract of the system in double precision.  The masses of particles
  ///        are the only templated type.
  ///
  /// \param tier  Get pointers at the level of the CPU host or GPU device
  const RMSDPlanReader<double> dpData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the read-only abstract of the system in single precision.  The masses of particles
  ///        are the only templated type.
  ///
  /// \param tier  Get pointers at the level of the CPU host or GPU device
  const RMSDPlanReader<float> spData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the coordinate synthesis pointer.
  const PhaseSpaceSynthesis* getCoordinateSynthesisPointer() const;

  /// \brief Get the system cache pointer.
  const SystemCache* getCachePointer() const;

  /// \brief Get a pointer to the map between the cache and the coordinate synthesis.
  const SynthesisCacheMap* getSynthesisMapPointer() const;

#ifdef STORMM_USE_HPC
  /// \brief Upload the object's contents to the GPU device.
  void upload();

  /// \brief Download the object's contents from the GPU device.
  void download();
#endif

  /// \brief Attach a SystemCache and map between the cache and the synthesis to the object.
  ///
  /// Overloaded:
  ///   - Specify the cache and map by const pointer
  ///   - Specify the cache and map by const reference
  ///
  /// \param sc     The system cache, reflecting user specified input coordinates and topologies
  /// \param scmap  The map between the synthesis and the system cache
  /// \{
  void addSystemCache(const SystemCache *sc, const SynthesisCacheMap *scmap);
  void addSystemCache(const SystemCache &sc, const SynthesisCacheMap &scmap);
  /// \}
  
private:
  
  /// The number of distinct plans, one per unique topology this object is established to serve
  int plan_count;

  // Additional directives and arrays pertinent to each individual plan
  RMSDMethod general_strategy;         ///< Align systems by mass (ALIGN_MASS) or by coordinates of
                                       ///<   all particles (ALIGN_GEOM), no do not align and take
                                       ///<   the mass-weighted (NO_ALIGN_MASS) or the
                                       ///<   coordinate-averaged (NO_ALIGN_GEOM) RMSD
  double required_mass_fraction;       ///< The mass fraction (or number of atoms, if mass is not a
                                       ///<   a consideration)
  Hybrid<double> masses;               ///< Masses of each particle represented in double precision
  Hybrid<float> sp_masses;             ///< Masses of each particle represented in single precision
  Hybrid<int> atom_counts;             ///< Overall numbers of atoms in each unique topology
  Hybrid<int> atom_starts;             ///< Starting indices of each unique system in the masses
                                       ///<   and sp_masses arrays.  The masses are padded up to
                                       ///<   the warp size in each system.
  Hybrid<int> asymmetric_core_atoms;   ///< Atoms found in each system's asymmetric core, given in
                                       ///<   terms of the system's own topological indices.
  Hybrid<int> asymmetric_core_counts;  ///< Atom counts in each system's core, not part of any of
                                       ///<   its symmetry-related groups
  Hybrid<int> asymmetric_core_starts;  ///< Starting indices in asymmetric_core_atoms for each
                                       ///<   system's list of asymmetric core atoms
  Hybrid<int> symmetry_group_atoms;    ///< Atoms involved in each collection of interchangeable
                                       ///<   groups needed to compute a best-fit RMSD.  These
                                       ///<   indicies work in terms of each system's local
                                       ///<   topology numbering.  Unlike the energy calculation
                                       ///<   work units, which focus on one system at a time and
                                       ///<   contain instructions for pulling out specific atoms
                                       ///<   from the overall list, symmetry group plans will
                                       ///<   describe how to compare one snapshot of a system to
                                       ///<   another in terms of the one system's structure.
  Hybrid<int4> symmetry_group_bounds;  ///< Bounds array for symmetry_group_atoms, indicating the
                                       ///<   limits of atom indices involved in each collection of
                                       ///<   equivalent atoms.  Because each system's
                                       ///<   symmetry-related groups are defined back-to-back,
                                       ///<   without padding, each element of this array will give
                                       ///<   the starting index of symmetry_group_atoms in the "x"
                                       ///<   member and the upper limit in the "y" member.
                                       ///<   Between systems, however, the symmetry-related atom
                                       ///<   lists will be padded by the warp size.  The "z"
                                       ///<   member of each element gives the specific size of
                                       ///<   domains in each symmetry-related group (the number of
                                       ///<   symmetry-related atoms involved in each swap), while
                                       ///<   the "w" member gives the number of domains (the
                                       ///<   number of interchangeable sub-structures).
  Hybrid<int> alignment_steps;         ///< Integer cast of EquivalenceSwap enumerations for
                                       ///<   the groups of each system.  This is indexes in step
                                       ///<   with symmetry_group_bounds above and each entry can
                                       ///<   be thought of as a fifth member of the tuples in that
                                       ///<   array.
  Hybrid<int2> symmetry_group_ranges;  ///< A bounds array on symmetry_group_bounds, indicating
                                       ///<   the lower and upper limits of symmetry groups
                                       ///<   relevant to each system in its "x" and "y" members,
                                       ///<   respectively.

  /// The number of symmetry-related groups that must be placed in order to ensure that the
  /// required mass fraction of each system is satisfied before subsequent groups are placed.
  /// Because the symmetry-related groups are arranged in order of decreasing size, these numbers
  /// indicate the number of groups to take from the front of the list.  When running the
  /// combinatorial permutations of all such groups, the progress will be stored in an unsigned
  /// 32-bit integers devoting four bits to each symmetry-related group.  This will break if any
  /// of the groups have greater than 16-fold symmetry, which is unheard of in chemistry, and at
  /// most eight distinct groups can be sampled in combinatorial fashion.  At the low end, this
  /// implies at most 2^8 = 256 graph isomorphisms to test for the core alignment, and if just
  /// half of the groups have three-fold, as opposed to two-fold, symmetry, the number of
  /// isomorphism tests rises to 3^4 * 2^4 = 1296, an enormous number of computations just to lay
  /// the foundation for a single RMSD calculation.  But, it implies that a single warp running
  /// each RMSD calculation standas a good chance of working a compute-bound problem.
  Hybrid<int> alignment_prerequisites;

  /// Indicate the approach that must be taken to compute the overall RMSD.  This is an integer
  /// re-cast of the RMSDAlignmentProtocol markers for each system.
  Hybrid<int> alignment_protocols;

  /// The number of symmetry levels indicates the depth of rearrangement that any given atom could
  /// undergo while evaluating different settings for symmmetry-equivalent groups and their
  /// dependents.  The symmetry-equivalent groups to which each atom belongs at any given level of
  /// symmetry are recorded in the array symmetry_group_membership.
  Hybrid<int> symmetry_depth;

  /// The membership of atoms in each symmetry-equivalent group for a given level of evaluation.
  /// The symmetry group is given in the "x" member and the domain of the group in the "y" member.
  /// Entries of -1 in the "x" member indicate that, at that level of symmetry, the atom is not
  /// involved in any symmetry group.  There are entries of this array for every atom of the system
  /// and all possible layers of the system.  Any given atom may be a member of as many distinct
  /// symmetry groups as the system has symmetry layers, but all of those groups will be
  /// dependents, in some tree, upon one symmetry group in particular.  Each system's entries are
  /// padded by the warp size.
  Hybrid<int2> symmetry_group_membership;

  /// Bounds array for symmetry_group_membership above.  The multiplicative effect that the number
  /// of symmetry layers may have on array sizes, plus the need to store additional information
  /// about the upper bound of each atom's applicable groups, forces the type to be ullint rather
  /// than int, even though the total number of atoms in the synthesis cannot exceed INT_MAX.  Each
  /// entry stores the lower limit of the particular atom's symmetry group indices in the low 48
  /// bits.  The number of groups that the atom is included in makes up the high 16 bits.
  /// If SGMB[i] is the ith element of symmetry_group_membership_bounds, then for most atoms,
  ///
  ///   (SGMB[i] & 0xffffffffffffULL) + ((SGMB[i] >> 48) & 0xffffULL) =
  ///   (SGMB[i + 1] & 0xffffffffffffULL).
  ///
  /// In the padded spaces between systems, that relationship will not hold, but it is simply
  /// illustrative.
  Hybrid<ullint> symmetry_group_membership_bounds;
  
  /// List of all topologies described by this plan
  std::vector<AtomGraph*> ag_pointers; 

  /// A pointer to the coordinate synthesis upon which the plan is based (if this is the nullptr,
  /// then there is no synthesis for the plan and the number of unique systems will be expected to
  /// be one).
  PhaseSpaceSynthesis *poly_ps_ptr;

  /// A pointer to the systems cache upon which the synthesis might be based (if this is the
  /// nullptr, there is no cache available and meta-data on systems cannot be used in devising
  /// plans for which systems shall be compared).  This is only valid if a synthesis pointer is
  /// provided.
  SystemCache *sc_ptr;

  /// Map between the systems cache and the associated synthesis
  SynthesisCacheMap *scmap_ptr;

  /// \brief Resize the storage arrays and add groups of symmetry-related atoms to the object's
  ///        Hybrid arrays.
  ///
  /// \param eq_tables  List of all relevant AtomEquivalence objects
  void chartSymmetryGroups(const std::vector<AtomEquivalence> &eq_tables);
  
  /// \brief Compute codes for each atom which, given a series of numbers indicating the
  ///        configurations that each symmetry-related group is to take on, will guide the
  ///        creation of a molecular isomorphism suitable for a direct RMSD calculation.  Greater
  ///        description of the method is found in the implementation.
  void writeSymmetryGroupCodes();
};

} // namespace structure
} // namespace stormm

#include "rmsd_plan.tpp"

#endif
