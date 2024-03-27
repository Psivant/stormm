// -*-c++-*-
#ifndef STORMM_ATOMGRAPH_SYNTHESIS_H
#define STORMM_ATOMGRAPH_SYNTHESIS_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/gpu_details.h"
#include "Constants/behavior.h"
#include "Constants/generalized_born.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/reduction_enumerators.h"
#include "Math/reduction_workunit.h"
#include "Potential/energy_enumerators.h"
#include "Restraints/restraint_apparatus.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_enumerators.h"
#include "Topology/topology_util.h"
#include "UnitTesting/stopwatch.h"
#include "static_mask_synthesis.h"
#include "synthesis_abstracts.h"
#include "valence_workunit.h"

namespace stormm {
namespace synthesis {

using card::GpuDetails;
using card::Hybrid;
using card::HybridTargetLevel;
using constants::ExceptionResponse;
using energy::ValenceKernelSize;
using stmath::RdwuPerSystem;
using stmath::ReductionWorkUnit;
using restraints::RestraintApparatus;
using restraints::RestraintKit;
using topology::AtomGraph;
using topology::AtomicRadiusSet;
using topology::ChemicalDetailsKit;
using topology::getRealParameters;
using topology::ImplicitSolventModel;
using topology::MassForm;
using topology::SettleSetting;
using topology::ShakeSetting;
using topology::UnitCellType;
using testing::StopWatch;
using namespace generalized_born_defaults;

/// \brief A collection of one or more AtomGraph objects, with similar components arranged in
///        contiguous arrays (often padded by the GPU warp size to prevent one system from flowing
///        into another).  Work units covering all systems are laid out in additional, contiguous
///        arrays.  While individual topologies (AtomGraphs) have Hybrid POINTER-kind objects for
///        details such as charges and atom type indices and all of the POINTER-kind objects
///        targeted a single ARRAY-kind object of the correct memory type, the synthesis, which may
///        have many topologies and a great deal more overall information, stores most of its data
///        in a series of ARRAY-kind objects, one for each member variable.
class AtomGraphSynthesis {
public:

  /// \brief The constructor takes a series of topologies and NMR restraints.  The NMR restraints
  ///        point to specific topologies and thereby apply to any coordinate sets that also point
  ///        to those topologies.
  ///
  /// Overloaded:
  ///   - Prepare the energy surface with or without restraints
  ///   - Accept lists of unique topologies and restraints with auxiliary lists of indices for
  ///     replicating them, or accept explicit lists of topologies and restraints with one index
  ///     for each system.
  ///   
  /// \param topologies_in         List of input topology pointers
  /// \param restraints_in         List of input restraint collections
  /// \param topology_indices_in   Indicators of which topologies describe each system in this
  ///                              synthesis.  This allows a synthesis to describe many copies
  ///                              of each system in its work units while storing relatively small
  ///                              amounts of data.
  /// \param restraint_indices_in  Indicators of which collections of restraints supplement the
  ///                              energy surfaces of each system within this synthesis.  Different
  ///                              sets of restraints allow many systems referencing the same
  ///                              topology to evolve on different energy surfaces.
  /// \param policy                Instruction on what to do if questionable input is encountered
  /// \param vwu_atom_limit        Maximum number of atoms to include in each valence work unit,
  ///                              taken from user input or a pre-computed value based on known
  ///                              launch bounds of the kernel that will evaluate the work units
  /// \param timer_in              Timer to track the wall time of this setup
  /// \{
  AtomGraphSynthesis(const std::vector<AtomGraph*> &topologies_in,
                     const std::vector<RestraintApparatus*> &restraints_in,
                     const std::vector<int> &topology_indices_in,
                     const std::vector<int> &restraint_indices_in,
                     ExceptionResponse policy_in = ExceptionResponse::WARN,
                     const GpuDetails &gpu = null_gpu, StopWatch *timer_in = nullptr);

  AtomGraphSynthesis(const std::vector<AtomGraph*> &topologies_in,
                     const std::vector<RestraintApparatus*> &restraints_in,
                     ExceptionResponse policy_in = ExceptionResponse::WARN,
                     const GpuDetails &gpu = null_gpu,
                     StopWatch *timer_in = nullptr);
  
  AtomGraphSynthesis(const std::vector<AtomGraph*> &topologies_in,
                     const std::vector<int> &topology_indices_in,
                     ExceptionResponse policy_in = ExceptionResponse::WARN,
                     const GpuDetails &gpu = null_gpu, StopWatch *timer_in = nullptr);

  AtomGraphSynthesis(const std::vector<AtomGraph*> &topologies_in,
                     ExceptionResponse policy_in = ExceptionResponse::WARN,
                     const GpuDetails &gpu = null_gpu,
                     StopWatch *timer_in = nullptr);
  /// \}

  /// \brief The presence of POINTER-kind Hybrid objects necessitates manually written copy and
  ///        move constructors as well as copy and move assignment operators, but all are valid.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right hand side of the assignment statement
  /// \{
  AtomGraphSynthesis(const AtomGraphSynthesis &original);
  AtomGraphSynthesis(AtomGraphSynthesis &&original);
  AtomGraphSynthesis& operator=(const AtomGraphSynthesis &original);
  AtomGraphSynthesis& operator=(AtomGraphSynthesis &&original);
  /// \}
  
  /// \brief Get the number of unique topologies described by the synthesis
  int getUniqueTopologyCount() const;

  /// \brief Get a topology pointer for a specific system contained within the synthesis.
  ///
  /// Overloaded:
  ///   - Get the topology pointer for a specific system
  ///   - Get topology pointers for all systems
  ///
  /// \param system_index  Index of the system of interest
  /// \{
  const AtomGraph* getSystemTopologyPointer(int system_index) const;
  const std::vector<AtomGraph*> getSystemTopologyPointer() const;
  /// \}

  /// \brief Get a const reference to the list of all pointers for unique topologies making up
  ///        this synthesis.
  const std::vector<AtomGraph*>& getUniqueTopologies() const;

  /// \brief Get a const reference to the list of all topology indices, indicating the topologies
  ///        describing each systems contained within this synthesis.
  std::vector<int> getTopologyIndices() const;
  
  /// \brief Get a restraint apparatus pointer for a sepcific system contained within the
  ///        synthesis.
  ///
  /// \param system_index  Index of the system of interest
  RestraintApparatus* getSystemRestraintPointer(int system_index) const;
  
  /// \brief Get the number of systems that this synthesis describes
  int getSystemCount() const;

  /// \brief Get the number of atoms in one or more systems.
  ///
  /// Overloaded:
  /// - Get the total number of atoms, summed over all systems, including replicas.  This does not
  ///   include padding.
  /// - Get the total number of atoms in a specific sytem, without padding.
  ///
  /// \param system_index  index of the system of interest
  /// \{
  int getAtomCount() const;
  int getAtomCount(int system_index) const;
  /// \}

  /// \brief Get the entire padded number of atoms covering all systems.
  int getPaddedAtomCount() const;

  /// \brief Get the starting point for atoms of a specific system in the lineup of all topologies.
  ///
  /// \param system_index  index of the system of interest
  int getAtomOffset(int system_index) const;
  
  /// \brief Get the total number of virtual sites across all systems, including replicas.
  int getVirtualSiteCount() const;

  /// \brief Get the total number of bond terms across all systems, including replicas.
  int getBondTermCount() const;

  /// \brief Get the total number of bond angle terms across all systems, including replicas.
  int getAngleTermCount() const;
  
  /// \brief Get the total number of cosine-based dihedral terms across all systems and replicas.
  int getDihedralTermCount() const;
  
  /// \brief Get the total number of Urey-Bradley terms across all systems and replicas.
  int getUreyBradleyTermCount() const;
  
  /// \brief Get the total number of CHARMM improper terms across all systems and replicas.
  int getCharmmImproperTermCount() const;
  
  /// \brief Get the total number of CMAP terms across all systems and replicas.
  int getCmapTermCount() const;

  /// \brief Get the number of unique atom types (a parameter, not an extensive quantity dependent
  ///        on the number of systems).
  int getLJTypeCount() const;

  /// \brief Get the number of unique charge parameters.
  int getChargeTypeCount() const;

  /// \brief Get the number of unique harmonic bond parameter sets.
  int getBondParameterCount() const;
  
  /// \brief Get the number of unique harmonic bond angle parameter sets.
  int getAngleParameterCount() const;
  
  /// \brief Get the number of unique cosine-based dihedral parameter sets.
  int getDihedralParameterCount() const;

  /// \brief Get the number of unique Urey-Bradley harmonic angle parameter sets.
  int getUreyBradleyParameterCount() const;

  /// \brief Get the number of unique CHARMM improper parameter sets.
  int getCharmmImproperParameterCount() const;

  /// \brief Get the number of unique CMAP surfaces.
  int getCmapSurfaceCount() const;

  /// \brief Get the sizes of all individual systems as a const reference to the Hybrid array
  ///        member variable.
  const Hybrid<int>& getSystemAtomCounts() const;
  
  /// \brief Get the starting locations of all individual systems as a const reference to the
  ///        Hybrid array member variable.
  const Hybrid<int>& getSystemAtomOffsets() const;

  /// \brief Get the numbers of unconstrained degrees of freedom for all individual systems in the
  //         synthesis, as a const reference to the internal Hybrid array member variable.
  const Hybrid<int>& getDegreesOfFreedom() const;
  
  /// \brief Get the numbers of constrained degrees of freedom for all individual systems in the
  //         synthesis, as a const reference to the internal Hybrid array member variable.
  const Hybrid<int>& getConstrainedDegreesOfFreedom() const;
  
  /// \brief Get the unit cell type that will be taken for all systems (TRICLINIC subsumes
  ///        ORTHORHOMBIC in a sort of type promotion).
  UnitCellType getUnitCellType() const;

  /// \brief Get the implicit solvent model in use across all systems.
  ImplicitSolventModel getImplicitSolventModel() const;

  /// \brief Get the dielectric constant (supporting the implicit solvent model) for all systems.
  double getDielectricConstant() const;

  /// \brief Get the salt concentration (supporting the implicit solvent model) for all systems.
  double getSaltConcentration() const;

  /// \brief Get the fundamental Coulomb constant defining the electrostatics of all systems.
  double getCoulombConstant() const;

  /// \brief Get the name of the PB radii set for one or more systems.
  ///
  /// Overloaded:
  ///   - Get the PB radii set for all systems
  ///   - Get the PB radii set for a series of systems between low and high limits
  ///   - Get the PB radii set for a specific system
  /// \{
  std::vector<AtomicRadiusSet> getPBRadiiSet() const;
  std::vector<AtomicRadiusSet> getPBRadiiSet(int low_limit, int high_limit) const;
  AtomicRadiusSet getPBRadiiSet(int index) const;
  /// \}

  /// \brief Get partial charges stored within the synthesis.
  ///
  /// Overloaded:
  ///   - Get partial charges for all systems, padded by the warp size between systems
  ///   - Get partial charges for a single system
  ///   - Get partial charges for a specific range of atoms in a single system
  ///
  /// \param tier          Obtain the data from arrays on the host or the device
  /// \param system_index  Index of the specific system to query
  /// \param low_index     Lower limit of atoms in the system to query (will be validated)
  /// \param high_index    Upper limit of atoms in the system to query (will be validated)
  /// \{
  template <typename T>
  std::vector<T> getPartialCharges(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  template <typename T>
  std::vector<T> getPartialCharges(int system_index,
                                   HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  template <typename T>
  std::vector<T> getPartialCharges(int system_index, int low_index, int high_index,
                                   HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Get the masses (or inverse masses) of atoms in the synthesis.
  ///
  /// Overloaded:
  ///   - Get masses for all systems, padded by the warp size between systems
  ///   - Get masses for a single system
  ///   - Get masses for a specific range of atoms in a single system
  ///
  /// \param tier          Obtain the data from arrays on the host or the device
  /// \param system_index  Index of the specific system to query
  /// \param low_index     Lower limit of atoms in the system to query (will be validated)
  /// \param high_index    Upper limit of atoms in the system to query (will be validated)
  /// \{
  template <typename T> std::vector<T> getAtomicMasses(HybridTargetLevel tier) const;
  template <typename T> std::vector<T> getAtomicMasses(HybridTargetLevel tier,
                                                       int system_index) const;
  template <typename T> std::vector<T> getAtomicMasses(HybridTargetLevel tier, int system_index,
                                                       int low_index, int high_index) const;
  /// \}

  /// \brief Get the overall number of valence work units needed to account for interactions in
  ///        all systems.
  int getValenceWorkUnitCount() const;

  /// \brief Get the maximum size of valence work units.  This will have been either set by the
  ///        user or tailored by the automated heuristics to produce the best saturation.
  int getValenceWorkUnitSize() const;

  /// \brief Get the necessary thread block size for evaluating the valence work units.
  ValenceKernelSize getValenceThreadBlockSize() const;
  
  /// \brief Get the abstracts (condensed lists of import and instruction set limits) for the
  ///        valence work units spanning all systems in this synthesis.
  const Hybrid<int2>& getValenceWorkUnitAbstracts() const;
  
  /// \brief Get the type of non-bonded work required by systems in this synthesis.
  NbwuKind getNonbondedWorkType() const;

  /// \brief Get the number of non-bonded work units serving systems in this synthesis.
  int getNonbondedWorkUnitCount() const;

  /// \brief Get the depth of the random numbers cache that the non-bonded work units of this
  ///        topology synthesis will try to fill.  This will raise a warning and return zero if
  ///        no cache depth has yet been set.
  int getRandomCacheDepth() const;
  
  /// \brief Get the number of reduction work units spanning all systems.
  int getReductionWorkUnitCount() const;

  /// \brief Get a qualitative assessment of the number of reduction work units assigned to any
  ///        one system.
  RdwuPerSystem getRdwuPerSystem() const;

  /// \brief Get the reduction work unit abstracts (this is not required for valence or non-bonded
  ///        work units as the abstracts needed for implementing these work units are produced by
  ///        member functions of this class (see below).  The reduction work, however, is managed
  ///        by an abstract splicing together elements of this class, PhaseSpaceSynthesis, and a
  ///        small class allocating temporary storage space.
  const Hybrid<int>& getReductionWorkUnitAbstracts() const;

  /// \brief Get a minimal kit with double-precision parameter detail for computing valence
  ///        interactions for all systems based on the work units stored in this object.
  ///
  /// \param tier  Level at which to obtain pointers for the abstract
  SyValenceKit<double>
  getDoublePrecisionValenceKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a minimal kit with single-precision parameter detail for computing valence
  ///        interactions for all systems based on the work units stored in this object.
  ///
  /// \param tier  Level at which to obtain pointers for the abstract
  SyValenceKit<float>
  getSinglePrecisionValenceKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a minimal kit with double-precision parameter detail for computing valence
  ///        interactions for all systems based on the work units stored in this object.
  ///
  /// \param tier  Level at which to obtain pointers for the abstract
  SyRestraintKit<double, double2, double4>
  getDoublePrecisionRestraintKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a minimal kit with single-precision parameter detail for computing valence
  ///        interactions for all systems based on the work units stored in this object.
  ///
  /// \param tier  Level at which to obtain pointers for the abstract
  SyRestraintKit<float, float2, float4>
  getSinglePrecisionRestraintKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a minimal kit with double-precision real numbers for computing non-bonded
  ///        interactions for all systems based on work units stored in this object.
  ///
  /// \param tier  Level at which to obtain pointers for the abstract
  SyNonbondedKit<double, double2>
  getDoublePrecisionNonbondedKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
  /// \brief Get a minimal kit with single-precision real numbers for computing non-bonded
  ///        interactions for all systems based on work units stored in this object.
  ///
  /// \param tier  Level at which to obtain pointers for the abstract
  SyNonbondedKit<float, float2>
  getSinglePrecisionNonbondedKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a minimal kit with double-precision real numbers for updating atom and virtual
  ///        site positions based on work units stored in this object.
  ///
  /// \param tier  Level at which to obtain pointers for the abstract
  SyAtomUpdateKit<double, double2, double4>
  getDoublePrecisionAtomUpdateKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a minimal kit with single-precision real numbers for updating atom and virtual
  ///        site positions based on work units stored in this object.
  ///
  /// \param tier  Level at which to obtain pointers for the abstract
  SyAtomUpdateKit<float, float2, float4>
  getSinglePrecisionAtomUpdateKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a const pointer to the object itself, useful for retrieving a valid pointer when
  ///        the object is available as a const reference.
  const AtomGraphSynthesis* getSelfPointer() const;
  
#ifdef STORMM_USE_HPC
  /// \brief Upload the object
  void upload();

  /// \brief Download the object
  void download();
#endif
  
  /// \brief Construct non-bonded work units for all unique topologies (there are no restraints
  ///        for non-bonded interactions that might distinguish systems with the same topology, as
  ///        was a consideration when developing the valence work units).  Load the instructions
  ///        into the topology synthesis for availability on the GPU.
  ///
  /// \param poly_se              Synthesis of static exclusion masks for a compilation of systems
  ///                             in isolated boundary conditions.
  /// \param init_request         Indicate a pattern for the non-bonded work units to initialize
  ///                             accumulators for subsequent force and energy calculations.
  /// \param random_cache_depth   Number of random values to store for each atom in the synthesis
  ///                             (the actual number of values stored is multiplied by three, for
  ///                             Cartesian X, Y, and Z contributions)
  /// \param gpu                  Details of the GPU in use (this may change the profile of the
  ///                             workload)
  void loadNonbondedWorkUnits(const StaticExclusionMaskSynthesis &poly_se,
                              InitializationTask init_request = InitializationTask::NONE,
                              int random_cache_depth = 0, const GpuDetails &gpu = null_gpu);

  /// \brief Apply an implicit solvent model to the synthesis.  Any mode of operation other than
  ///        taking the original topologies' parameters as is (calling this member function with
  ///        no arguments) will cause changes ot the underlying topologies.
  ///
  /// Overloaded:
  ///   - Take the implicit solvent models as described in the original topologies
  ///   - Apply specific atomic radius sets to each system
  ///   - Apply a common atomic radius set ot all systems
  ///   - Accept a GB neck model, its double-precision abstract, or no such model
  ///
  /// \param igb_in  A specific, established implicit solvent model (from literature and hard-coded
  ///                into STORMM) to apply
  /// \param dielectric_in  The desired dielectric constant
  /// \param saltcon_in     The intended salt concentration (affects the GB decay parameter Kappa)
  /// \param radii_sets_in  Radii to impart to the topology (this is often coupled to the choice of
  ///                       implicit solvent model, but for purposes of experimentation or new
  ///                       model development might be flexible)
  /// \param policy         Indicator of what to do if the topology's PB radii to not meet the
  ///                       implicit solvent model requirements, or there is some other problem
  /// \{
  void setImplicitSolventModel();
  
  void setImplicitSolventModel(ImplicitSolventModel igb_in,
                               const NeckGeneralizedBornKit<double> &ngbk,
                               const std::vector<AtomicRadiusSet> &radii_sets_in,
                               double dielectric_in = 80.0, double saltcon_in = 0.0,
                               ExceptionResponse policy = ExceptionResponse::WARN);

  void setImplicitSolventModel(ImplicitSolventModel igb_in,
                               const NeckGeneralizedBornKit<double> &ngbk,
                               AtomicRadiusSet radii_set = AtomicRadiusSet::NONE,
                               double dielectric_in = 80.0, double saltcon_in = 0.0,
                               ExceptionResponse policy = ExceptionResponse::WARN);

  void setImplicitSolventModel(ImplicitSolventModel igb_in,
                               const NeckGeneralizedBornTable &ngb_tab,
                               const std::vector<AtomicRadiusSet> &radii_sets_in,
                               double dielectric_in = 80.0, double saltcon_in = 0.0,
                               ExceptionResponse policy = ExceptionResponse::WARN);

  void setImplicitSolventModel(ImplicitSolventModel igb_in,
                               const NeckGeneralizedBornTable &ngb_tab,
                               AtomicRadiusSet radii_set = AtomicRadiusSet::NONE,
                               double dielectric_in = 80.0, double saltcon_in = 0.0,
                               ExceptionResponse policy = ExceptionResponse::WARN);
  /// \}
  
private:

  // The first member variables pertain to totals across all systems: atoms, potential function
  // terms, and sizes of composite parameter arrays that span all topologies
  ExceptionResponse policy;       ///< Indicator of what to do with bad or questionable input
  int topology_count;             ///< Number of unique topologies represented in this synthesis
  int restraint_network_count;    ///< Number of unique restraint apparatuses in this synthesis
  int system_count;               ///< Number of independent coordinate sets described by this
                                  ///<   synthesis
  int total_atoms;                ///< Total number of atoms, spanning all topologies
  int total_virtual_sites;        ///< Total number of atoms, spanning all topologies
  int total_bond_terms;           ///< Total number of bonds, spanning all topologies
  int total_angl_terms;           ///< Total number of bond angles, spanning all topologies
  int total_dihe_terms;           ///< Total number of dihedrals, spanning all topologies
  int total_ubrd_terms;           ///< Total number of Urey-Bradley angles, spanning all topologies
  int total_cimp_terms;           ///< Total number of CHARMM impropers, spanning all topologies
  int total_cmap_terms;           ///< Total number of CMAPs, spanning all topologies
  int total_lj_types;             ///< Total number of unique Lennard-Jones interaction types,
                                  ///<   spanning all topologies
  int total_charge_types;         ///< Total number of unique atomic partial charges, spanning all
                                  ///<   topologies
  int total_bond_params;          ///< Total number of unique bond parameter sets
  int total_angl_params;          ///< Total number of unique bond angle parameter sets
  int total_dihe_params;          ///< Total number of unique dihedral parameter sets
  int total_ubrd_params;          ///< Total number of unique Urey-Bradley angle parameter sets
  int total_cimp_params;          ///< Total number of unique CHARMM improper parameter sets
  int total_cmap_surfaces;        ///< Total number of unique CMAP surfaces
  int total_attn14_params;        ///< Total number of unique 1:4 attenuated interaction pairs
                                  ///<   (electrostatic and van-der Waals)
  int total_vste_params;          ///< Total number of unique virtual site frames, across all
                                  ///<   systems
  int total_position_restraints;  ///< Total number of positional restraints
  int total_distance_restraints;  ///< Total number of distance restraints
  int total_angle_restraints;     ///< Total number of angle restraints
  int total_dihedral_restraints;  ///< Total number of dihedral restraints
  
  // The presence of periodic boundaries or an implicit solvent must be common to all systems, as
  // must the cutoff and other run parameters, although these are contained in other data
  // structres.  If an implicit solvent model is in effect, its parameters are taken as common to
  // all systems.  Other facets of the the system behavior, such as bond constraints, must also
  // be common to all systems.
  UnitCellType periodic_box_class;    ///< The type of unit cell.  The type can be NONE, but all
                                      ///<   systems must then have no unit cell boundaries.  If
                                      ///<   all but one system is ORTHORHOMBIC and the final
                                      ///<   system is TRICLINIC, the cell type will be TRICLINIC.
  ImplicitSolventModel gb_style;      ///< The flavor of Generalized Born (or other implicit
                                      ///<   solvent) to use, i.e. Hawkins / Cramer / Truhlar.  All
                                      ///<   systems must share the same implicit solvent model
                                      ///<   type, and the unit cell type must be NONE if any
                                      ///<   implicit solvent model is in effect.
  double dielectric_constant;         ///< Dielectric constant to take for implicit solvent
                                      ///<   calculations
  double is_kappa;                    ///< Kappa scaling parameter used in some GB models (and
                                      ///<   possibly other implicit solvents)
  double salt_concentration;          ///< Salt concentration affecting implicit solvent models
  double gb_offset;                   ///< Standard offset for taking GB radii from PB radii
  double gb_neckscale;                ///< Scaling factor inherent in "neck" GB approximations
  double gb_neckcut;                  ///< Cutoff distance for apply GB "neck" corrections
  double coulomb_constant;            ///< Coulomb's constant in units of kcal-A/mol-e^2 (Amber
                                      ///<   differs from other programs in terms of what this is,
                                      ///<   so it can be set here)
  ShakeSetting use_bond_constraints;  ///< Toggles use of bond length constraints
  SettleSetting use_settle;           ///< Toggles analytic constraints on rigid water
  int largest_constraint_group;       ///< Number of bonds involved in the largest hub-and-spoke
                                      ///<   constraint group for any system in the synthesis
  char4 water_residue_name;           ///< Name of water residue, compared to residue_names

  /// Settings for the Poisson-Boltzmann radii sets for each system, also used in GB.  While it is
  /// typical to use a prescribed radii set with a particular GB model, the choice is not limited.
  /// Therefore, each system may operate with a different set of radii, provided that they pass
  /// any particular requirements of the common GB model in use.
  std::vector<AtomicRadiusSet> pb_radii_sets;

  /// An array of pointers to the individual topologies that form the basis of this work plan.
  /// This array is topology_count in length and accessible only on the host as it is used to
  /// organize work units, not process the actual energy calculations.
  std::vector<AtomGraph*> topologies;

  /// An array of pointers to the individual restraint apparatuses that form the basis of this work
  /// plan.  This array is restraint_apparatus_count in length and, like the topologies array, is
  /// only accessible on the host.
  std::vector<RestraintApparatus*> restraint_networks;

  /// Dummy restraint apparatuses, created and stored in the object for each unique topology to
  /// stand in where external restraint apparatuses were not supplied.  This is critical to make
  /// persistent objects that will be available for the lifetime of the AtomGraphSynthesis.
  std::vector<RestraintApparatus> restraint_dummies;
  
  /// An array of indices for the source topologies guiding the motion of each system (this array
  /// is system_count in length)
  Hybrid<int> topology_indices;

  /// An array of indices for the restraint collections guiding the motion of each system (this
  /// array is system_count in length)
  Hybrid<int> restraint_indices;
  
  // The following arrays are POINTER-kind objects, each directed at a segment of a data array
  // for storing critical topology sizes.  There is one of each of these numbers for every
  // system in this work plan: each of the Hybrid arrays that follows is system_count in length.
  // All of the descriptions below follow from the eponymous member variables in the AtomGraph
  // object and should be taken to mean "... for each topology in the work plan." The pointers
  // target the ARRAY-kind object int_system_data, at intervals of the number of systems in the
  // work plan, with no padding for the warp size.  If these POINTER-kind objects list bounds, the
  // bounds themselves will reflect padding for the HPC warp size.
  Hybrid<int> atom_counts;             ///< Total number of atoms and virtual sites in each system
  Hybrid<int> residue_counts;          ///< Total number of residues, including solvent molecules
                                       ///<   in each system
  Hybrid<int> molecule_counts;         ///< Total number of molecules in each system
  Hybrid<int> largest_residue_sizes;   ///< Number of atoms in the largest residue of each system
  Hybrid<int> last_solute_residues;    ///< Last residue of the solutes, indexed according to the
                                       ///<   overall order in this synthesis rather than the
                                       ///<   original topologies
  Hybrid<int> last_solute_atoms;       ///< Last atoms of the solutes (a typical solute is a
                                       ///<   biomolecule), indexed according to the overall order
                                       ///<   in this synthesis rather than the original topologies
  Hybrid<int> first_solvent_molecules; ///< First molecule in what is deemed to be solvent in each
                                       ///<   system, indexed according to the overall order in
                                       ///<   this synthesis rather than the original topologies

  // Valence term and off-center particle quantities
  Hybrid<int> ubrd_term_counts;     ///< Number of Urey-Bradley angle stretch terms in each system
  Hybrid<int> cimp_term_counts;     ///< Number of CHARMM impropers in each system
  Hybrid<int> cmap_term_counts;     ///< Number of CMAP terms in each system
  Hybrid<int> bond_term_counts;     ///< Numbers of bonded interactions in each system
  Hybrid<int> angl_term_counts;     ///< Numbers of bond angle interactions in each system
  Hybrid<int> dihe_term_counts;     ///< Numbers of dihedral cosine terms in each system
  Hybrid<int> virtual_site_counts;  ///< Number of virtual sites / extra points out of all atoms
                                    ///<   per system

  // Restraint tallies for each system
  Hybrid<int> posn_restraint_counts;  ///< Number of positional restraints per system
  Hybrid<int> bond_restraint_counts;  ///< Number of distance restraints per system
  Hybrid<int> angl_restraint_counts;  ///< Number of three-point angle restraints per system
  Hybrid<int> dihe_restraint_counts;  ///< Number of four-point dihedral restraints per system
  
  // Information relevant to non-bonded calculations
  Hybrid<int> lj_type_counts;          ///< Number of distinct Lennard-Jones types in each system
  Hybrid<int> total_exclusion_counts;  ///< Total number of non-bonded exclusions, including 1-4

  // Information relevant to the MD propagation algorithm
  Hybrid<int> rigid_water_counts;       ///< Number of rigid water molecules subject to SETTLE
  Hybrid<int> bond_constraint_counts;   ///< Bonds with lengths constrained by SHAKE or RATTLE
  Hybrid<int> degrees_of_freedom;       ///< Total degrees of freedom in the system, in the absence
                                        ///<   of constraints (3N -6 for N particles)
  Hybrid<int> cnst_degrees_of_freedom;  ///< Total degrees of freedom, 3N - 6 - constraints in
                                        ///<   each system
  Hybrid<int> nonrigid_particle_counts; ///< A rigid water is one non-rigid particle.  A protein
                                        ///<   with N atoms and no bond length constraints is N
                                        ///<   particles.

  // The following are bounds arrays for the larger data arrays below: more segments of memory,
  // each the number of systems in length.  The contents of these arrays are, with some
  // exceptions, prefix sums over the associated counts arrays above, taking into account zero
  // padding for the stuff they are providing indices into.  A series of systems with atom counts
  // 31, 59, and 78 with an HPC warp size of 32 would generate an atom bounds array stating 0, 32,
  // and 96.  There is no capping value, and the (k+1)th bound minus the kth bounds does not
  // provide information on the exact number associated with the kth system.  Programs should read
  // the exact number of atoms or force field terms per system and the starting bound from one of
  // the arrays below.  The arrays for valence terms spanning all systems provide the global
  // indices of atoms in unified (but not re-ordered) arrays of atom descriptors and also
  // coordinates.  In this sense, these arrays are intermediates between the valence work units
  // that follow and the original AtomGraph objects used to construct this AtomGraphSynthesis.
  Hybrid<int> atom_offsets;               ///< Starting indices for each system's various atomic
                                          ///<   descriptor arrays, which will have eponymous
                                          ///<   names to their AtomGraph counterparts.
  Hybrid<int> atom_bit_offsets;           ///< Starting indices for each system's bit masks.
  Hybrid<int> residue_offsets;            ///< Starting indices for each system's residue
                                          ///<   descriptor arrays
  Hybrid<int> molecule_offsets;           ///< Starting indices for each system's molecule
                                          ///<   descriptor arrays
  Hybrid<int> ubrd_term_offsets;          ///< Starting indices for Urey-Bradley parameter index
                                          ///<   lists for each system
  Hybrid<int> cimp_term_offsets;          ///< CHARMM improper parameter index lists
  Hybrid<int> cmap_term_offsets;          ///< CMAP atom indexing lists
  Hybrid<int> bond_term_offsets;          ///< Bond stretching term parameter index lists
  Hybrid<int> angl_term_offsets;          ///< Angle bending term parameter index lists
  Hybrid<int> dihe_term_offsets;          ///< Dihedral term parameter index lists
  Hybrid<int> virtual_site_offsets;       ///< Starting indices for each system's virtual site
                                          ///<   descriptor lists
  Hybrid<int> posn_restraint_offsets;     ///< Starting indices for each system's positional
                                          ///<   restraints
  Hybrid<int> bond_restraint_offsets;     ///< Starting indices for each system's distance
                                          ///<   restraints
  Hybrid<int> angl_restraint_offsets;     ///< Starting indices for each system's three-point angle
                                          ///<   restraints
  Hybrid<int> dihe_restraint_offsets;     ///< Starting indices for each system's four-point
                                          ///<   dihedral angle restraints
  Hybrid<int> sett_group_offsets;         ///< Starting indices for SETTLE-constrained rigid waters
                                          ///<   in each system
  Hybrid<int> cnst_group_offsets;         ///< Starting indices for hub-and-spoke constraint groups
                                          ///<   in each system
  Hybrid<int> nb_exclusion_offsets;       ///< Starting indices for unified nonbonded exclusion
                                          ///<   lists
  Hybrid<int> lennard_jones_abc_offsets;  ///< Starting indices for each system's relevant
                                          ///<   Lennard-Jones tables (A + B, or C coefficients)
  
  /// Data array to hold all system size information, the target of all the above pointers
  Hybrid<int> int_system_data;          

  // Maps of each system and its distinguishable parts.  From here onwards, the AtomGraphSynthesis
  // keeps a combination of POINTER-kind and ARRAY-kind Hybrid Object member variables, whatever is
  // most convenient for grouping and storing data as it is produced.  As above, many of the
  // names match what is in the AtomGraph object.  The arrays are laid out such that each system's
  // details occur in one contiguous stretch, padded by blank elements up to a multiple of the
  // HPC warp size.  Some of these arrays are indices themselves, but will still need a bounds
  // array, i.e. the point at which to start reading residue limits for the kth system.  Bounds
  // arrays for individual systems within each of these objects are written into the collected
  // integer system data array above, accessed through one of the preceding Hybrid POINTER-kind
  // objects.

  // Atom and residue details
  Hybrid<int2> residue_limits;            ///< Atomic indices marking boundaries of individual
                                          ///<   residues: lower limit in the x index and upper
                                          ///<   limit in the y index
  Hybrid<int> atom_struc_numbers;         ///< Structure atom numbers, such as those taken from a
                                          ///<   PDB file
  Hybrid<int> residue_numbers;            ///< Residue numbers, such as those taken from a PDB file
  Hybrid<int2> molecule_limits;           ///< Atomic indices marking the boundaries of molecules
                                          ///<   in the molecule_contents array: lower limit in the
                                          ///<   x index and upper limit in the y index
  Hybrid<int> atomic_numbers;             ///< Atomic numbers for atoms in all systems
  Hybrid<int> mobile_atoms;               ///< Atom mobility masks for each system
  Hybrid<int> molecule_membership;        ///< Distinct molecules, indexed from 0 for each system,
                                          ///<   indicating molecules to which each atom belongs
  Hybrid<int> molecule_contents;          ///< Indices of atoms making up each molecule, starting
                                          ///<   from the lower bound atom index for each system,
                                          ///<   bounded by molecule_limits
  Hybrid<double> atomic_charges;          ///< Partial charges on all atoms and virtual sites
  Hybrid<double> atomic_masses;           ///< Masses of all atoms
  Hybrid<double> inverse_atomic_masses;   ///< Inverse masses for all atoms
  Hybrid<float> sp_atomic_charges;        ///< Partial charges on all atoms (single precision)
  Hybrid<float> sp_atomic_masses;         ///< Masses of all atoms (single precision)
  Hybrid<float> sp_inverse_atomic_masses; ///< Inverse masses for all atoms (single precision)
  Hybrid<char4> atom_names;               ///< Four letter names of all atoms, i.e. PDB names
  Hybrid<char4> atom_types;               ///< Four letter names of all atom types, i.e. CT
  Hybrid<char4> residue_names;            ///< Four letter names of all residues, i.e. PDB names
  Hybrid<int> chem_int_data;              ///< Chemistry-related integer data
  Hybrid<int2> chem_int2_data;            ///< Chemistry-related int2 data
  Hybrid<double> chem_double_data;        ///< Chemistry-related double-precision data
  Hybrid<float> chem_float_data;          ///< Chemistry-related single-precision data
  Hybrid<char4> chem_char4_data;          ///< Chemistry-related char4 data

  // Valence parameter maps and bounds, showing how the parameters of each unique topology map to
  // those of the consensus tables in the synthesis.
  Hybrid<int> ubrd_param_map;            ///< Urey-Bradley parameter index maps
  Hybrid<int2> ubrd_param_map_bounds;    ///< Bounds array for unique topologies in the
                                         ///<   Urey-Bradley parameter maps (dimension
                                         ///<   topology_count, a series of low and high limits
                                         ///<   for each unique topology's Urey-Bradley parameter
                                         ///<   index correspondence)
  Hybrid<int> cimp_param_map;            ///< CHARMM improper parameter index maps
  Hybrid<int2> cimp_param_map_bounds;    ///< Bounds array for CHARMM improper parameter index maps
  Hybrid<int> cmap_param_map;            ///< CMAP surface index maps
  Hybrid<int2> cmap_param_map_bounds;    ///< Bounds array for CMAP surface index maps
  Hybrid<int> bond_param_map;            ///< Harmonic bond index maps
  Hybrid<int2> bond_param_map_bounds;    ///< Bounds array for harmonic bond index maps
  Hybrid<int> angl_param_map;            ///< Harmonic angle index maps
  Hybrid<int2> angl_param_map_bounds;    ///< Bounds array for harmonic angle index maps
  Hybrid<int> dihe_param_map;            ///< Cosine-based dihedral index maps
  Hybrid<int2> dihe_param_map_bounds;    ///< Bounds array for cosine-based dihedral index maps
  Hybrid<int> attn14_param_map;          ///< Cosine-based dihedral index maps
  Hybrid<int2> attn14_param_map_bounds;  ///< Bounds array for cosine-based dihedral index maps  
  Hybrid<int> vste_param_map;            ///< Virtual site frame specifications mapping 
  Hybrid<int2> vste_param_map_bounds;    ///< Bounds array for vste_param_map
  Hybrid<int> sett_param_map;            ///< SETTLE group geometry mapping
  Hybrid<int2> sett_param_map_bounds;    ///< Buunds array for SETTLE group geometry mapping
  Hybrid<int> cnst_param_map;            ///< Constraint group parameter set mapping
  Hybrid<int2> cnst_param_map_bounds;    ///< Buunds array for constraint group parameter set
                                         ///<   mapping
  
  // CHARMM and basic force field valence term details.  Each of these objects points into one of
  // the data arrays at the bottom of the section.
  Hybrid<double> ubrd_stiffnesses;      ///< Stiffness constant of each Urey-Bradley stretch
  Hybrid<double> ubrd_equilibria;       ///< Equilibrium length of each Urey-Bradley stretch
  Hybrid<double> cimp_stiffnesses;      ///< CHARMM impropers are harmonic, too!
  Hybrid<double> cimp_phase_angles;     ///< The "equilibria" for CHARMM impropers
  Hybrid<int> cmap_surface_dimensions;  ///< Dimensions for every unique CMAP surface
  Hybrid<int> cmap_surface_bounds;      ///< Bounds array for every unique CMAP surface
  Hybrid<int> cmap_patch_bounds;        ///< Bounds array for CMAP patches
  Hybrid<double> cmap_surfaces;         ///< Concatenated, column-major format matrices for
                                        ///<   every CMAP surface term
  Hybrid<double> cmap_patches;          ///< Concatenated, column-major format matrices for each
                                        ///<   CMAP aggregated patch representation.  See the
                                        ///<   description in Topology/atomgraph.h.
  Hybrid<float> sp_ubrd_stiffnesses;    ///< Stiffness constant of each Urey-Bradley stretch,
                                        ///<   in single precision
  Hybrid<float> sp_ubrd_equilibria;     ///< Equilibrium length of each Urey-Bradley stretch,
                                        ///<   in single precision
  Hybrid<float> sp_cimp_stiffnesses;    ///< CHARMM impropers are harmonic, too!
  Hybrid<float> sp_cimp_phase_angles;   ///< The "equilibria" for CHARMM impropers
  Hybrid<float> sp_cmap_surfaces;       ///< Concatenated, column-major format matrices for every
                                        ///<   CMAP surface term, in single precision
  Hybrid<float> sp_cmap_patches;        ///< Concatenated, column-major format matrices for each
                                        ///<   CMAP aggregated patch representation in single
                                        ///<   precision (see description in Topology/atomgraph.h).

  // Basic force field valence term details, consensus tables form all topologies
  Hybrid<double> bond_stiffnesses;       ///< Stiffness of each bond stretch, kcal/mol-A^2
  Hybrid<double> bond_equilibria;        ///< Equilibrium lengths of all bonds, A
  Hybrid<double> angl_stiffnesses;       ///< Stiffness of each angle bend, kcal/mol-rad^2
  Hybrid<double> angl_equilibria;        ///< Equilibrium angle for all bending terms, radians
  Hybrid<double> dihe_amplitudes;        ///< Amplitudes of each dihedral cosine term, kcal/mol
  Hybrid<double> dihe_periodicities;     ///< Periodicity of each dihedral / torsion cosine term
  Hybrid<double> dihe_phase_angles;      ///< Phase angle of each dihedral / torsion cosine term
  Hybrid<double> attn14_elec_factors;    ///< Attenuated 1:4 non-bonded interaction electrostatic
                                         ///<   scaling factors
  Hybrid<double> attn14_vdw_factors;     ///< Attenuated 1:4 non-bonded interaction van-der Waals
                                         ///<   scaling factors
  Hybrid<float> sp_bond_stiffnesses;     ///< Stiffness of each bond stretch (single precision)
  Hybrid<float> sp_bond_equilibria;      ///< Equilibrium lengths of all bonds (single precision)
  Hybrid<float> sp_angl_stiffnesses;     ///< Angle bending stiffnesses (single precision)
  Hybrid<float> sp_angl_equilibria;      ///< Angle bending equilibria (single precision)
  Hybrid<float> sp_dihe_amplitudes;      ///< Amplitudes of torsion cosine terms (single precision)
  Hybrid<float> sp_dihe_periodicities;   ///< Periodicities of torsion terms (single precision)
  Hybrid<float> sp_dihe_phase_angles;    ///< Phase angles of torsion terms (single precision)
  Hybrid<float> sp_attn14_elec_factors;  ///< Attenuated 1:4 non-bonded interaction electrostatic
                                         ///<   scaling factors
  Hybrid<float> sp_attn14_vdw_factors;   ///< Attenuated 1:4 non-bonded interaction van-der Waals
                                         ///<   scaling factors
  Hybrid<double> valparam_double_data;   ///< Valence parameter double-pecision data
  Hybrid<float> valparam_float_data;     ///< Valence parameter single-pecision data
  Hybrid<int> valparam_int_data;         ///< Valence parameter integer data (for CMAP sizes and
                                         ///<   bounds)
  Hybrid<int2> valparam_int2_data;       ///< Valence parameter int2 data (this is for the
                                         ///<   parameter bounds arrays, i.e.
                                         ///<   bond_param_map_bounds)
  
  // Valence term indexing arrays, all of them indexing atoms in the synthesis list, updated from
  // the indexing in their original topologies and parameters in the condensed tables of the
  // synthesis.  All of the following arrays are POINTER-kind objects targeting valence_int_data.
  Hybrid<int> ubrd_i_atoms;      ///< Urey-Bradley I atoms
  Hybrid<int> ubrd_k_atoms;      ///< Urey-Bradley K atoms
  Hybrid<int> ubrd_param_idx;    ///< Urey-Bradley parameter indices
  Hybrid<int> cimp_i_atoms;      ///< CHARMM improper I atoms
  Hybrid<int> cimp_j_atoms;      ///< CHARMM improper J atoms
  Hybrid<int> cimp_k_atoms;      ///< CHARMM improper K atoms (the center atom of the quartet)
  Hybrid<int> cimp_l_atoms;      ///< CHARMM improper L atoms
  Hybrid<int> cimp_param_idx;    ///< CHARMM improper parameter indices
  Hybrid<int> cmap_i_atoms;      ///< Correction map I atoms
  Hybrid<int> cmap_j_atoms;      ///< Correction map J atoms
  Hybrid<int> cmap_k_atoms;      ///< Correction map K atoms
  Hybrid<int> cmap_l_atoms;      ///< Correction map L atoms
  Hybrid<int> cmap_m_atoms;      ///< Correction map M atoms
  Hybrid<int> cmap_param_idx;    ///< CMAP surface indices
  Hybrid<int> bond_i_atoms;      ///< Harmonic bond I atoms
  Hybrid<int> bond_j_atoms;      ///< Harmonic bond J atoms
  Hybrid<int> bond_param_idx;    ///< Bond parameter indices
  Hybrid<int> angl_i_atoms;      ///< Harmonic angle I atoms
  Hybrid<int> angl_j_atoms;      ///< Harmonic angle J atoms
  Hybrid<int> angl_k_atoms;      ///< Harmonic angle K atoms
  Hybrid<int> angl_param_idx;    ///< Bond angle parameter indices
  Hybrid<int> dihe_i_atoms;      ///< Cosine-based dihedral (proper or improper) I atoms
  Hybrid<int> dihe_j_atoms;      ///< Cosine-based dihedral (proper or improper) J atoms
  Hybrid<int> dihe_k_atoms;      ///< Cosine-based dihedral (proper or improper) K atoms
  Hybrid<int> dihe_l_atoms;      ///< Cosine-based dihedral (proper or improper) L atoms
  Hybrid<int> dihe_param_idx;    ///< Cosine-based dihedral parameter indices
  Hybrid<int> valence_int_data;  ///< Array targeted by the preceding POINTER-kind hybrid objects
                                 ///<   in this section
  
  // Non-bonded parameter indexing and van-der Waals tables
  Hybrid<int> charge_indices;                   ///< Atomic charge indices, 0 to
                                                ///<   charge_type_count - 1, pointed to target
                                                ///<   chem_int_data as a convenient and sensible
                                                ///<   way to keep the data stored
  Hybrid<int> lennard_jones_indices;            ///< Lennard-Jones indices, 0 to
                                                ///<   lj_type_count - 1, pointed to target
                                                ///<   chem_int_data
  Hybrid<double> charge_parameters;             ///< Unique charge parameters for all systems.
                                                ///<   These are pooled over all systems, and the
                                                ///<   union of all unique charges is queried
                                                ///<   for each atom according to charge_indices.
  Hybrid<double> lennard_jones_a_coeff;         ///< Lennard-Jones A coefficients, a series
                                                ///<   of tables covering all systems
  Hybrid<double> lennard_jones_b_coeff;         ///< Lennard-Jones B coefficients, a series
                                                ///<   of tables covering all systems
  Hybrid<double> lennard_jones_c_coeff;         ///< Lennard-Jones C coefficients, a series of
                                                ///<   tables covering all systems
  Hybrid<double> lennard_jones_14_a_coeff;      ///< Lennard-Jones A coefficients, a series
                                                ///<   of tables covering all systems
  Hybrid<double> lennard_jones_14_b_coeff;      ///< Lennard-Jones B coefficients, a series
                                                ///<   of tables covering all systems
  Hybrid<double> lennard_jones_14_c_coeff;      ///< Lennard-Jones C coefficients, a series of
                                                ///<   tables covering all systems
  Hybrid<double> lennard_jones_sigma;           ///< Lennard-Jones sigma values, a series of tables
                                                ///<   covering all systems
  Hybrid<double> lennard_jones_14_sigma;        ///< Lennard-Jones 1:4 sigma values, a series of
                                                ///<   tables covering all systems
  Hybrid<float> sp_charge_parameters;           ///< Unique charge parameters for all systems
                                                ///<    (single precision)
  Hybrid<float> sp_lennard_jones_a_coeff;       ///< Lennard-Jones A coefficients, a series of
                                                ///<   tables covering all systems (single
                                                ///<   precision)
  Hybrid<float> sp_lennard_jones_b_coeff;       ///< Lennard-Jones B coefficients, a series of
                                                ///<   tables covering all systems (single
                                                ///<   precision)
  Hybrid<float> sp_lennard_jones_c_coeff;       ///< Lennard-Jones C coefficients, a series of
                                                ///<   tables covering all systems (single
                                                ///<   precision)
  Hybrid<float> sp_lennard_jones_14_a_coeff;    ///< Lennard-Jones A coefficients, a series of
                                                ///<   tables covering all systems (single
                                                ///<   precision)
  Hybrid<float> sp_lennard_jones_14_b_coeff;    ///< Lennard-Jones B coefficients, a series of
                                                ///<   tables covering all systems (single
                                                ///<   precision)
  Hybrid<float> sp_lennard_jones_14_c_coeff;    ///< Lennard-Jones C coefficients, a series of
                                                ///<  tables covering all systems (single
                                                ///<   precision)
  Hybrid<float> sp_lennard_jones_sigma;         ///< Lennard-Jones sigma values, a series of tables
                                                ///<   covering all systems (single precision)
  Hybrid<float> sp_lennard_jones_14_sigma;      ///< Lennard-Jones 1:4 sigma values, a series of
                                                ///<   tables covering all systems (single
                                                ///<   precision)

  // Implicit solvent model parameters.  Like the atomic partial charges arrays, most of these will
  // target chem_double_data, chem_float_data, and (in the case of the indices) chem_int_data,
  // for convenience.  The neck tables are their own ARRAY-type objects to expedite applying a
  // new implicit solvent model to an existing synthesis.
  int neck_table_size;                   ///< Dimension of the Neck Generalized Born limit tables.
                                         ///<   All indices into the limits tables, defined in
                                         ///<   neck_gb_indices (below), must be less than this
                                         ///<   value.
  Hybrid<int> neck_gb_indices;           ///< Indicies into separation and maximum value parameter
                                         ///<   tables for Mongan's "neck" GB implementations
  Hybrid<double> atomic_pb_radii;        ///< Radii of all atoms according to the pb_radii_set
  Hybrid<double> gb_screening_factors;   ///< Generalized Born screening factors for all atoms
  Hybrid<double> gb_alpha_parameters;    ///< Generalized Born alpha parameters for each atom (set
                                         ///<   according to a particular GB scheme, not part of
                                         ///<   an Amber topology file but added later)
  Hybrid<double> gb_beta_parameters;     ///< Generalized Born beta parameters for each atom
  Hybrid<double> gb_gamma_parameters;    ///< Generalized Born gamma parameters for each atom
  Hybrid<double2> neck_limit_tables;     ///< Neck maximum separations (the maximum distance at
                                         ///<   which the neck correction applies) in the x member
                                         ///<   and the maximum value of the neck correction in the
                                         ///<   y member.  Neither table is symmetric.  Both are
                                         ///<   copied into the synthesis for the purpose of
                                         ///<   interlacing the values, which are needed in close
                                         ///<   proximity throughout the "neck" GB code.
  Hybrid<float> sp_atomic_pb_radii;      ///< P.B. Radii of all atoms (single precision)
  Hybrid<float> sp_gb_screening_factors; ///< Generalized Born screening factors (single precision)
  Hybrid<float> sp_gb_alpha_parameters;  ///< Single-precision Generalized Born alpha parameters
  Hybrid<float> sp_gb_beta_parameters;   ///< Single-precision Generalized Born beta parameters
  Hybrid<float> sp_gb_gamma_parameters;  ///< Single-precision Generalized Born gamma parameters
  Hybrid<float2> sp_neck_limit_tables;   ///< Single-precision neck GB limit tables
  
  // NMR restraint terms and details: these function exactly like other parameter sets and are
  // indexed by lists of atoms in the bond work units arrays.  They can be included in the
  // synthesis of AtomGraphs due to their nature as potential terms, whereas the original
  // topologies had to be read from files that did not contain such terms.
  Hybrid<int2> rposn_step_bounds;    ///< Initial (x member) and final (y member) step numbers for
                                     ///<   applying positional restraints
  Hybrid<int2> rbond_step_bounds;    ///< Initial (x member) and final (y member) step numbers for
                                     ///<   applying positional restraints
  Hybrid<int2> rangl_step_bounds;    ///< Initial (x member) and final (y member) step numbers for
                                     ///<   applying positional restraints
  Hybrid<int2> rdihe_step_bounds;    ///< Initial (x member) and final (y member) step numbers for
                                     ///<   applying positional restraints
  Hybrid<double2> rposn_init_k;      ///< Initial stiffnesses for time-dependent positional
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<double2> rposn_final_k;     ///< Final stiffnesses for time-dependent positional
                                     ///<   restraints (ignored for time-independent restraints)
  Hybrid<double4> rposn_init_r;      ///< Initial displacements for time-dependent positional
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<double4> rposn_final_r;     ///< Final displacments for time-dependent positional
                                     ///<   restraints (ignored for time-independent restraints)
  Hybrid<double2> rposn_init_xy;     ///< Initial X and Y Cartesian coordinates for the target
                                     ///<   location of time-dependent positional restraints, or
                                     ///<   the static values of time-independent restraints
  Hybrid<double> rposn_init_z;       ///< Initial Z Cartesian coordinates for the target
                                     ///<   location of time-dependent positional restraints, or
                                     ///<   the static values of time-independent restraints
  Hybrid<double2> rposn_final_xy;    ///< Final X and Y Cartesian coordinates for the target
                                     ///<   location of time-dependent positional restraints, or
                                     ///<   the static values of time-independent restraints
  Hybrid<double> rposn_final_z;      ///< Final Z Cartesian coordinates for the target location of
                                     ///<   time-dependent positional restraints, or the static
                                     ///<   values of time-independent restraints
  Hybrid<double2> rbond_init_k;      ///< Initial stiffnesses for time-dependent distance
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<double2> rbond_final_k;     ///< Final stiffnesses for time-dependent distance restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<double4> rbond_init_r;      ///< Initial displacements for time-dependent distance
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<double4> rbond_final_r;     ///< Final displacments for time-dependent distance restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<double2> rangl_init_k;      ///< Initial stiffnesses for time-dependent angle restraints,
                                     ///<   or the static values of time-independent restraints
  Hybrid<double2> rangl_final_k;     ///< Final stiffnesses for time-dependent angle restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<double4> rangl_init_r;      ///< Initial displacements for time-dependent angle
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<double4> rangl_final_r;     ///< Final displacments for time-dependent angle restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<double2> rdihe_init_k;      ///< Initial stiffnesses for time-dependent dihedral
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<double2> rdihe_final_k;     ///< Final stiffnesses for time-dependent dihedral restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<double4> rdihe_init_r;      ///< Initial displacements for time-dependent dihedral
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<double4> rdihe_final_r;     ///< Final displacments for time-dependent dihedral restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<float2> sp_rposn_init_k;    ///< Initial stiffnesses for time-dependent positional
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<float2> sp_rposn_final_k;   ///< Final stiffnesses for time-dependent positional
                                     ///<   restraints (ignored for time-independent restraints)
  Hybrid<float4> sp_rposn_init_r;    ///< Initial displacements for time-dependent positional
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<float4> sp_rposn_final_r;   ///< Final displacments for time-dependent positional
                                     ///<   restraints (ignored for time-independent restraints)
  Hybrid<float2> sp_rposn_init_xy;   ///< Initial X and Y Cartesian coordinates for the target
                                     ///<   location of time-dependent positional restraints, or
                                     ///<   the static values of time-independent restraints
  Hybrid<float> sp_rposn_init_z;     ///< Initial Z Cartesian coordinates for the target
                                     ///<   location of time-dependent positional restraints, or
                                     ///<   the static values of time-independent restraints
  Hybrid<float2> sp_rposn_final_xy;  ///< Final X and Y Cartesian coordinates for the target
                                     ///<   location of time-dependent positional restraints, or
                                     ///<   the static values of time-independent restraints
  Hybrid<float> sp_rposn_final_z;    ///< Final Z Cartesian coordinates for the target location of
                                     ///<   time-dependent positional restraints, or the static
                                     ///<   values of time-independent restraints
  Hybrid<float2> sp_rbond_init_k;    ///< Initial stiffnesses for time-dependent distance
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<float2> sp_rbond_final_k;   ///< Final stiffnesses for time-dependent distance restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<float4> sp_rbond_init_r;    ///< Initial displacements for time-dependent distance
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<float4> sp_rbond_final_r;   ///< Final displacments for time-dependent distance restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<float2> sp_rangl_init_k;    ///< Initial stiffnesses for time-dependent angle restraints,
                                     ///<   or the static values of time-independent restraints
  Hybrid<float2> sp_rangl_final_k;   ///< Final stiffnesses for time-dependent angle restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<float4> sp_rangl_init_r;    ///< Initial displacements for time-dependent angle
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<float4> sp_rangl_final_r;   ///< Final displacments for time-dependent angle restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<float2> sp_rdihe_init_k;    ///< Initial stiffnesses for time-dependent dihedral
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<float2> sp_rdihe_final_k;   ///< Final stiffnesses for time-dependent dihedral restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<float4> sp_rdihe_init_r;    ///< Initial displacements for time-dependent dihedral
                                     ///<   restraints, or the static values of time-independent
                                     ///<   restraints
  Hybrid<float4> sp_rdihe_final_r;   ///< Final displacments for time-dependent dihedral restraints
                                     ///<   (ignored for time-independent restraints)
  Hybrid<double> nmr_double_data;    ///< Double-precision information pertianing to (NMR)
                                     ///<   restraints of the kinds delinated above.  The Hybrid
                                     ///<   objects above are POINTER-kind, targeting arrays like
                                     ///<   this one.
  Hybrid<double2> nmr_double2_data;  ///< Double-precision double tuple information pertianing to
                                     ///<   (NMR) restraints of the kinds delinated above.  The
                                     ///<   Hybrid objects above are POINTER-kind, targeting arrays
                                     ///<   like this one.
  Hybrid<double4> nmr_double4_data;  ///< Double-precision quadruple tuple information pertianing
                                     ///<   to (NMR) restraints of the kinds delinated above.  The
                                     ///<   Hybrid objects above are POINTER-kind, targeting arrays
                                     ///<   like this one.
  Hybrid<float> nmr_float_data;      ///< Single-precision information pertianing to (NMR)
                                     ///<   restraints of the kinds delinated above.  The Hybrid
                                     ///<   objects above are POINTER-kind, targeting arrays like
                                     ///<   this one.
  Hybrid<float2> nmr_float2_data;    ///< Single-precision double tuple information pertianing to
                                     ///<   (NMR) restraints of the kinds delinated above.  The
                                     ///<   Hybrid objects above are POINTER-kind, targeting arrays
                                     ///<   like this one.
  Hybrid<float4> nmr_float4_data;    ///< Single-precision quadruple tuple information pertianing
                                     ///<   to (NMR) restraints of the kinds delinated above.  The
                                     ///<   Hybrid objects above are POINTER-kind, targeting arrays
                                     ///<   like this one.

  // Restraint indexing arrays, offering atom indices in the concatenated atom and restraint
  // parameter arrays of the synthesis.  The Hybrid objects in this section are POINTER-kind
  // objects targeting the nmr_int_data array, much like the parameters in the preceding section
  // target nmr_[type]_data.
  Hybrid<int> rposn_atoms;              ///< Atom indices involved in positional restraints
  Hybrid<int> rposn_kr_param_idx;       ///< Restraint parameters for k(2,3) and r(1,2,3,4) in
                                        ///<   each positional restraint
  Hybrid<int> rposn_xyz_param_idx;      ///< Restraint parameters for the target positions in each
                                        ///<   positional restraint
  Hybrid<int> rbond_i_atoms;            ///< Atom I indices involved in distance restraints
  Hybrid<int> rbond_j_atoms;            ///< Atom J indices involved in distance restraints
  Hybrid<int> rbond_param_idx;          ///< Restraint parameters for k(2,3) and r(1,2,3,4) in each
                                        ///<   distance restraint
  Hybrid<int> rangl_i_atoms;            ///< Atom I indices involved in angle restraints
  Hybrid<int> rangl_j_atoms;            ///< Atom J indices involved in angle restraints
  Hybrid<int> rangl_k_atoms;            ///< Atom K indices involved in angle restraints
  Hybrid<int> rangl_param_idx;          ///< Restraint parameters for k(2,3) and r(1,2,3,4) in each
                                        ///<   three-point angle restraint
  Hybrid<int> rdihe_i_atoms;            ///< Atom I indices involved in dihedral restraints
  Hybrid<int> rdihe_j_atoms;            ///< Atom J indices involved in dihedral restraints
  Hybrid<int> rdihe_k_atoms;            ///< Atom K indices involved in dihedral restraints
  Hybrid<int> rdihe_l_atoms;            ///< Atom L indices involved in dihedral restraints
  Hybrid<int> rdihe_param_idx;          ///< Restraint parameters for k(2,3) and r(1,2,3,4) in each
                                        ///<   four-point dihedral angle restraint
  Hybrid<int> rposn_kr_param_map;       ///< Positional restraint K-R parameter maps
  Hybrid<int> rposn_xyz_param_map;      ///< Positional restraint X-Y-Z parameter maps
  Hybrid<int2> rposn_param_map_bounds;  ///< Positional restraint parameter map bounds (this
                                        ///<   applies to both the K-R and X-Y-Z parameter maps as
                                        ///<   the individual topologies make every restraint a
                                        ///<   unique parameter set)
  Hybrid<int> rbond_param_map;          ///< Distance restraint parameter maps
  Hybrid<int2> rbond_param_map_bounds;  ///< Distance restraint parameter map bounds
  Hybrid<int> rangl_param_map;          ///< Three-point angle restraint parameter maps
  Hybrid<int2> rangl_param_map_bounds;  ///< Three-point angle restraint parameter map bounds
  Hybrid<int> rdihe_param_map;          ///< Four-point dihedral restraint parameter maps
  Hybrid<int2> rdihe_param_map_bounds;  ///< Four-point dihedral restraint parameter map bounds
  Hybrid<int> nmr_int_data;             ///< Array targeted by POINTER-kind objects in this section
  Hybrid<int2> nmr_int2_data;           ///< Array targeted by POINTER-kind objects in this section

  // Virtual site details: virtual sites are considered parameters, which like valence terms apply
  // to a small group of atoms and index into a table of parameters including the frame type and
  // up to three dimensional measurements.
  Hybrid<double4> virtual_site_parameters;    ///< Frame types for each virtual site parameter set
                                              ///<   (w member) plus up to three dimensions for
                                              ///<   each frame (stored in the x, y, and z members,
                                              ///<   respectively).  This array spans all unique
                                              ///<   frame types across all systems in the
                                              ///<   synthesis.  The frame type, an integer, is
                                              ///<   stored as a double-precision real, which for
                                              ///<   integers of this size is exact.
  Hybrid<float4> sp_virtual_site_parameters;  ///< Single-precision version of the virtual site
                                              ///<   frame specifications.
  Hybrid<int> virtual_site_atoms;             ///< Indices of virtual sites in the concatenated
                                              ///<   list of all atoms, taken from their original
                                              ///<   topological orders plus the system offset.
                                              ///<   There are as many entries in this array, and
                                              ///<   the frame atom indexing arrays that follow,
                                              ///<   as there are virtual sites in all systems of
                                              ///<   the synthesis.  This is at least as large,
                                              ///<   and probably much greater, than the length of
                                              ///<   the preceding frame type and dimension arrays.
  Hybrid<int> virtual_site_frame1_atoms;      ///< Frame 1 (parent) atoms for each virtual site,
                                              ///<   using a similar indexing scheme as above.
  Hybrid<int> virtual_site_frame2_atoms;      ///< Frame 2 atoms for each virtual site
  Hybrid<int> virtual_site_frame3_atoms;      ///< Frame 3 atoms for each virtual site
  Hybrid<int> virtual_site_frame4_atoms;      ///< Frame 4 atoms for each virtual site
  Hybrid<int> virtual_site_parameter_indices; ///< Parameter indices for each virtual site,
                                              ///<   indexing into the frame type and dimension
                                              ///<   arrays above.
  Hybrid<int> vsite_int_data;                 ///< Virtual site integer data, storing virtual site
                                              ///<   and frame atom indices as well as virtual
                                              ///<   site parameter sets, but not frame types

  // SETTLE constraint group parameter sets and atom indices
  Hybrid<double4> settle_group_geometry;    ///< SETTLE group geometric parameters (ra, rb, rc,
                                            ///<   and invra in the tuple's x, y, z, and w
                                            ///<   members, respectively)
  Hybrid<double2> settle_group_masses;      ///< SETTLE group inverse masses (mormt and mhrmt in
                                            ///<   the tuple's x and y members, respectively)
  Hybrid<float4> sp_settle_group_geometry;  ///< SETTLE group geometric parameters (ra, rb, rc,
                                            ///<   and invra in the tuple's x, y, z, and w
                                            ///<   members, respectively), single precision
  Hybrid<float2> sp_settle_group_masses;    ///< SETTLE group inverse masses (mormt and mhrmt in
                                            ///<   the tuple's x and y members, respectively),
                                            ///<   single precision
  Hybrid<int4> settle_group_indexing;       ///< Oxygen atom indices, first and second hydrogen
                                            ///<   atom indices, and SETTLE parameter indices in
                                            ///<   the tuple's x, y, z, and w members, respectively

  // Constraint group parameter sets and atom indices
  Hybrid<int> constraint_group_indices;       ///< Atom indices for all hub-and-spoke constraint
                                              ///<   groups, beginning with each group's central
                                              ///<   atom and then listing all peripheral atoms
  Hybrid<int2> constraint_group_bounds;       ///< Bounds for every constraint group in the
                                              ///<   constraint_group_indices array.  Due to the
                                              ///<   presence of more than one system in the
                                              ///<   synthesis, it is not possible to have the
                                              ///<   bounds simply listed one after another and
                                              ///<   also have padding between systems' group
                                              ///<   indices.  But, these bounds are not used in
                                              ///<   high-preformance contexts.
  Hybrid<int> constraint_group_param_idx;     ///< Parameter indices for all hub-and-spoke
                                              ///<   constraint groups, referencing the
                                              ///<   constraint_param_bounds array
  Hybrid<int> constraint_param_bounds;        ///< Bounds array for constraint_group_params and
                                              ///<   sp_constraint_group_params
  Hybrid<double2> constraint_group_params;    ///< Constraint parameters: the squared target bond
                                              ///<   length in the x member and inverse mass in the
                                              ///<   y member, for each constrained bond in the
                                              ///<   group
  Hybrid<float2> sp_constraint_group_params;  ///< Single precision constraint parameters: length
                                              ///<   length and inverse mass in the x and y
                                              ///<   members of the tuple, respectively
  
  // Valence work units (VWUs): work units providing instruction sets for the GPU to operate on a
  // continuous, non-rearranged list of atoms and implement all valence terms.  Each VWU pertains
  // to one and only one individual topology from the list and each topology will have at least one
  // VWU, or more depending on what is needed to cover all of its atoms.  The VWU imports a list of
  // atoms, then computes a series of bond, angle, dihedral, and other force field terms based on
  // the cached atoms.  Forces are accumulated on all atoms in __shared__ memory and then
  // contributed back to the global arrays.
  int total_valence_work_units;  ///< Total count of valence work units spanning all systems
  int valence_work_unit_size;    ///< Maximum size of the value work units

  /// The thread block size needed to launch valence work units associated with this synthesis.
  /// Larger thread blocks may be used, but this is the minimum recommended size.
  ValenceKernelSize valence_thread_block_size;
  
  /// Instruction sets for the bond work units, storing integers for the low (_L nomenclature) and
  /// high (_H nomeclature) limits of all the types of interactions in a given work unit.  Each
  /// VWU takes a stride of integers from this array.  The length of this stride is not a matter
  /// of the warp size.
  Hybrid<int2> vwu_instruction_sets;

  /// A list of atoms that the VWU shall import, with indices into the global array of atoms for
  /// all systems.  Each VWU may import up to 3/4 as many atoms as the kernel blocks have threads,
  /// and each VWU takes a stride of that many ints from this array.
  Hybrid<int> vwu_import_lists;

  /// A mask of the imported atoms that each VWU shall move (x member) and also update (y member)
  /// in terms of position and velocity.  Each bit is set to 1 if the atom must be moved or
  /// updated, 0 otherwise.  In terms of bitwise operations, for all elements of this array,
  /// the move mask completely subsumes the update mask, (x | y) = x.  This POINTER-kind object
  /// targets insr_uint2_data like other instructions sets.
  Hybrid<uint2> vwu_manipulation_masks;

  /// Instructions for bond stretching and Urey-Bradley interactions.  Each uint2 tuple contains
  /// two atom indices in the x member (bits 1-10 and 11-20) and the parameter index of the bond /
  /// Urey-Bradley term in the y member.  A flag in the 21st bit of the x member indicates whether
  /// the term is a Urey-Bradley item a harmonic bond, which direct parameter choices as well as
  /// energy accumulation.
  Hybrid<uint2> cbnd_instructions;

  /// Instructions for angle bending and three-body NMR restraint interactions.  Each uint2 tuple
  /// contains three atom indices in the x member (bits 1-10, 11-20, and 21-30) and the parameter
  /// index of the angle in the y member.
  Hybrid<uint2> angl_instructions;

  /// Instructions for dihedral and CHARMM improper dihedral interactions.  Each uint2 tuple
  /// contains three atom indices in the x member (bits 1-10, 11-20, and 21-30), plus a flag in
  /// bit 31 to indicate whether the form of the potential should be a cosine term (0) or the
  /// CHARMM harmonic improper (1) and another flag in bit 32, relevant only when energy is being
  /// computed and only if the term is not a CHARMM improper, to indicate whether to add the
  /// energy to a proper (0) or improper (1) accumulator.  A fourth atom index in the y member
  /// (bits 1-10) precedes the index for 1:4 scaling factors (bits 11-15, up to 32 unique
  /// combinations of electrostatic and Lennard-Jones scaling factors, including zero if no 1:4
  /// interaction should be computed).  The 16th bit indicates whether an additional cosine-based
  /// dihedral term overlaps with the same four atoms.  If so, its parameter index is given in an
  /// auxiliary array of uint values (cdhe_overtones), which is only accessed in cases where such
  /// overlaps occur.  The final bits (17-32) indiate the dihedral or CHARMM improper dihedral
  /// parameter index (up to 65535 of either with the 65536th slot being a "no interaction"
  /// placeholder).  These instructions can store a number of unique parameters larger than any
  /// known force field, but the limit applies to the number of unique parameters at work in a
  /// particular simulation.  In some extreme case where a model does have more than 65535 unique
  /// dihedral parameters, the parameter of the first dihedral can be set to 65535 (the placeholder
  /// 65536th parameter) and the secondary dihedral (the overtone) can then use a parameter index
  /// up to 1048576 with up to 4096 1:4 scaling factor pairs.
  Hybrid<uint2> cdhe_instructions;

  /// Overtones for various dihedrals that overlap with the same four atoms.  This allows a single
  /// dihedral computation to be extended, with 50% of the marginal read cost, little additional
  /// arithmetic given the cost of computing a single dihedral, and no additional atomic ops to
  /// L1 memory.  This POINTER-kind object targets insr_uint_data.
  Hybrid<uint> cdhe_overtones;
  
  /// Instructions for CMAP interactions.  Each uint2 tuple contains three atom indices in the x
  /// member (bits 1-10, 11-20, and 21-30, plus two more in the y member (bits 1-10 and 11-20).
  /// The CMAP parameter index is given in the final twelve bits of the y member, offering a total
  /// of 4,096 unique CMAPs that one AtomGraphSynthesis can handle.
  Hybrid<uint2> cmap_instructions;

  /// Instructions for inferred, 1:4 excluded interactions.  Each uint object holds a pair of atom
  /// indices for the I and L atoms (bits 1-10 and 11-20), plus an index into a potentially much
  /// larger list of 1:4 scaling factor pairs (up to 4096) than the 32 allowed in the standard
  /// dihedral / composite dihedral instruction.  This POINTER-kind object targets insr_uint_data.
  Hybrid<uint> infr14_instructions;
  
  /// Instructions for positional restraints.  Each uint2 tuple contains an atom index (bits 1-10),
  /// the index of the k(2,3) / r(1,2,3,4) and time dependence series (bits 11-31), and a flag to
  /// indicate whether the restraint is time-dependent (bit 32) in the x member, followed by the
  /// index of the target x, y, and z coordinates in the tuple's y member.  32 bits are needed to
  /// index the target x / y / z parameters in order to provide a sufficiently large set of unique
  /// restraints.
  Hybrid<uint2> rposn_instructions;

  /// Instructions for NMR two-body restraints.  Each uint2 tuple contains two atom indices in the
  /// x member (bits 1-10 and 11-20) followed by a flag to indicate whether the restraint is
  /// time-dependent (bit 32).  The index of the k(2,3) / r(1,2,3,4) series is stored in the
  /// tuple's y member.
  Hybrid<uint2> rbond_instructions;

  /// Instructions for NMR three-body restraints.  Each uint2 tuple contains two atom indices in
  /// the x member (bits 1-10, 11-20, and 21-30) followed by a flag to indicate whether the
  /// restraint is time-dependent (bit 32).  The index of the k(2,3) / r(1,2,3,4) series is stored
  /// in the tuple's y member.
  Hybrid<uint2> rangl_instructions;

  /// Instructions for NMR four-body restraints.  Each uint2 tuple contains three atom indices in
  /// the x member (bits 1-10, 11-20, and 21-30) followed by a flag to indicate whether the
  /// restraint is time-dependent (bit 32).  The fourth atom index is stored in the y member
  /// (bits 1-10).  The index of the k(2,3) / r(1,2,3,4) series is stored primarily in the tuple's
  /// y member (bits 11-32) for a total of 4,194,304 unique dihedral restraint k(2,3) / r(1,2,3,4)
  /// parameter sets in a single AtomGraphSynthesis.
  Hybrid<uint2> rdihe_instructions;

  /// Instructions for virtual site placement and force transfer to frame atoms.  Each uint2 tuple
  /// contains the index of the virtual site itself in the x member (bits 1-10), plus up to two
  /// frame atom indices (bits 11-20 and 21-30).  The third and fourth frame atom indices
  /// are found in the y member (if a third or fourth atom are even part of the frame).  The frame
  /// parameter index (which indicates the frame type and up to three frame dimensions) is
  /// presented in the final twelve bits (21-32) of the y member.  More can be done with the
  /// remaining bits to increase the number of available virtual site frame parameters, if needed.
  Hybrid<uint2> vste_instructions;

  /// Instructions for SETTLE group constraints.  These are relatively simple: three atom indices
  /// are encoded in bits 1-10, 11-20, and 21-30 of the x member, while the constraint parameter
  /// set (multiple SETTLE configurations are supported in one simulation) is stored in the
  /// y member.
  Hybrid<uint2> sett_instructions;

  /// Instructions for hub-and-spoke constraint groups.  These are more complex: the size of the
  /// constraint group itself is uncertain, but for practical purposes the number of participating
  /// atoms will be limited to sixteen.  Each instruction pertains to one bond of the group,
  /// providing the atom index of the central (heavy) atom in bits 1-10 of the x member and the
  /// atom index of the light atom (hydrogen) in bits 11-20 of the x member.  The lane index of the
  /// base thread for the constraint group is given in bits 21-28 and the number of constrained
  /// bonds is given in 29-32.  The parameter index indicating the target equilibrium length and
  /// inverse masses for this constraint group is given in the y member.
  Hybrid<uint2> cnst_instructions;

  // Instructions for evaluating energy terms--these are stored as extra arrays, as they will only
  // be called when evaluating energies, and offer a compact format in which each interaction in
  // the list is a bit in a bitstring, (1) indicating that the work unit should contribute the
  // interaction to its sum, (0) indicating that the work unit should omit the interaction.
  // All of these are POINTER-kind objects targeting insr_uint_data.
  Hybrid<uint> accumulate_cbnd_energy;   ///< Contribute harmonic bond or Urey-Bradley energies to
                                         ///<   the total
  Hybrid<uint> accumulate_angl_energy;   ///< Contribute harmonic angle energies to the total
  Hybrid<uint> accumulate_cdhe_energy;   ///< Contribute proper and improper cosine-based dihedral
                                         ///<   or CHARMM improper dihedral energies to the total
  Hybrid<uint> accumulate_cmap_energy;   ///< Contribute CMAP energies to the total
  Hybrid<uint> accumulate_infr14_energy; ///< Contribute CMAP energies to the total
  Hybrid<uint> accumulate_rposn_energy;  ///< Contribute positional restraint energies to the total
                                         ///<   (this information could be implicit in directives
                                         ///<   to log an update to the atom, but it is not a great
                                         ///<   bandwidth burden and is included for consistency)
  Hybrid<uint> accumulate_rbond_energy;  ///< Contribute distance restraint energies
  Hybrid<uint> accumulate_rangl_energy;  ///< Contribute three-point angle restraint energies
  Hybrid<uint> accumulate_rdihe_energy;  ///< Contribute four-point dihedral restraint energies
  
  /// Collected array of all uint instructions
  Hybrid<uint> insr_uint_data;

  /// Collected array of all uint2 instructions
  Hybrid<uint2> insr_uint2_data;

  // Non-bonded work units (NBWUs): the synthesis also stores information on how to carry out
  // non-bonded tiles (or domain decomposition workloads, depending on the nature of the boundary
  // conditions).  These are less detailed, overall, than the valence work units.
  int total_nonbonded_work_units;  ///< The total number of non-bonded work units (tile groups for
                                   ///<   simulations in isolated-boundary conditions or cell
                                   ///<   groups for simulations in periodic boundary conditions)
  NbwuKind nonbonded_work_type;    ///< Enumerator for the type of work to expect in non-bonded
                                   ///<   work units.  This will inform the nature of the kernel
                                   ///<   to call in order to compute non-bonded interactions.

  /// Abstracts for non-bonded work units.  The nature of this array varies with the setting of
  /// nonbonded_work_type above, but each kernel will know how to interpret data in this array for
  /// its own purposes.
  ///
  /// NbwuKind: TILE_GROUPS
  /// Each abstract is a series of 32 integers.  The contents include: (index 0) the total number
  /// of tiles to import (up to 20), (indices 1-20) the starting index of atoms in each tile side,
  /// (indices 21-25) the number of atoms to import in each tile, and (indices 26-27) the limits
  /// of tile instructions to carry out based on content in tile_group_instructions.
  ///
  /// NbwuKind: SUPERTILES
  /// Each abstract is a series of 4 integers: (index 0) the upper limit of atoms in the
  /// particular system (one of many within the synthesis) to which the supertile belongs,
  /// (indices 1-2) lower limits of atoms for the supertile's abscissa and ordinate axes, (index
  /// 4) the supertile map index within some associated StaticExclusionsMaskSynthesis
  /// (pre-inflated by the square of the tile lengths per supertile).  In supertile work units,
  /// all tiles pertain to the same system.
  ///
  /// NbwuKind: DOMAIN
  /// Abstracts for domain decompositions in neighbor list-based periodic systems.  Each abstract
  /// is a series of 64 integers, 50 of which are used.
  Hybrid<int> nonbonded_abstracts;

  /// Instructions for the non-bonded work.  Each element is an instruction to do one non-bonded
  /// tile, along with an index into the exclusion masks data array held within some associated
  /// mask object.  The exclusion masks are not kept as part of this topology synthesis, but
  /// abstracts for the non-bonded work will span it and the associated mask object.  
  Hybrid<uint2> nbwu_instructions;

  // Reduction work units.  Various processes in molecular mechanics require system-wide
  // accumulations of particular quantities, and a collection of any number of systems of any size
  // requires a sophisticated system to make that happen efficiently.  The reduction work units
  // are an ecumenical solution to all such reductions, but they may need to combine with different
  // substrates to do a specific process.
  int total_reduction_work_units;  ///< Number of reduction work units for the entire system
  RdwuPerSystem rdwu_per_system;   ///< Indication of the type of reduction operations that can
                                   ///<   be performed: are there just one, or multiple reduction
                                   ///<   operations to perform per system?
  
  /// Abstracts for reduction work units.  There are not a corresponding array of instructions, as
  /// each reduction work unit takes one abstract and proceeds to perform gather, scatter, or
  /// all-reduce operations on the stated range of atoms.
  Hybrid<int> reduction_abstracts;
  
  /// Timings data, for reporting purposes
  StopWatch* timer;

  /// \brief Repair pointers in the copy and copy assignment constructors.
  void rebasePointers();
  
  /// \brief Create a series of dummy restraints for all topologies, then redirect any systems
  ///        which have nullptr for the restraint network (RestraintApparatus) pointer to the
  ///        topology-specific set of dummy restraints.
  ///
  /// \param restraint_indices_in  Input restraint network indices, referencing the object member
  ///                              variable restraint_networks and passed down in the constructor
  /// \param topology_indices_in   Input topology indices, referencing the object member variable
  ///                              topologies and passed down in the constructor
  std::vector<int> createDummyRestraints(const std::vector<int> &restraint_indices_in,
                                         const std::vector<int> &topology_indices_in);
  
  /// \brief Check the central lists of topologies and system indices to ensure that the requested
  ///        synthesis is sane.  Condense the list of topologies by weeding out duplicates.
  ///        Re-align the list of systems to work with the condensed list of topologies and return
  ///        that result.
  ///
  /// \param topology_indices_in  List of topologies which describe each system involved in this
  ///                             synthesis
  std::vector<int> checkTopologyList(const std::vector<int> &topology_indices_in);

  /// \brief Check the central lists of restraint groups and system indices to ensure that the
  ///        requested synthesis always pairs restraint groups with systems referencing consistent
  ///        topologies.  Condense the list of restraint groups by weeding out duplicates.
  ///        Re-align the list of systems to work with the condensed list of restraint groups and
  ///        return that result.
  ///
  /// \param restraint_indices_in  List of restraint apparatuses (groups, collections) which
  ///                              guide each system in this synthesis
  /// \param topology_indices_in   List of topologies which describe each system involved in this
  ///                              synthesis.  Must be supplied in its original, unpruned form.
  std::vector<int> checkRestraintList(const std::vector<int> &restraint_indices_in,
                                      const std::vector<int> &topology_indices_in,
                                      const std::vector<int> &topology_index_rebase);
  
  /// \brief Check settings that must be consistent between all topologies in the synthesis.
  void checkCommonSettings();

  /// \brief Lay down arrays that will hold each system's atoms and terms: particles and
  ///        interactions, with updated indexing into the synthesis as opposed to any of the
  ///        original topologies, but as of yet without parameter indices.  This routine handles
  ///        valence terms, restraints, and virtual sites.
  ///
  /// \param topology_indices_in     Original list of system topology indices, oriented towards
  ///                                the input list of AtomGraph pointers
  /// \param topology_index_rebase   Re-arrangement of system topology indices meeting the
  ///                                condensed list of unique topologies created by
  ///                                checkTopologyList()
  /// \param new_restraint_indices   Refined list of restraint group indices, oriented towards
  ///                                the input list of RestraintApparatus pointers
  /// \param restraint_index_rebase  Re-arrangement of restraint group indices meeting the
  ///                                condensed list of unique restraints created by
  ///                                checkRestraintList()
  void buildAtomAndTermArrays(const std::vector<int> &topology_indices_in,
                              const std::vector<int> &topology_index_rebase,
                              const std::vector<int> &new_restraint_indices,
                              const std::vector<int> &restraint_index_rebase);

  /// \brief Find unique parameters for each of the valence energy components, charges, and
  ///        virtual site frames.  Then, assign terms indexing atoms in the order of the synthesis
  ///        to use parameters from consensus tables.  This routine handles valence terms, virtual
  ///        sites, and charge interactions.
  void condenseParameterTables();
  
  /// \brief Extend the Lennard-Jones tables with a specific Lennard-Jones atom type, including
  ///        pair-specific cross-terms (i.e. CHARMM NB-fix details).  All cross terms involving the
  ///        parameter with other known parameters must be accounted for.  This can be an involved
  ///        process, but if the current topology's tables completely subsume the former topology's
  ///        tables or are completely covered by the current topology's tables it is tractable.
  void extendLJMatrices();

  /// \brief Find unique restraint k(2,3) and r(1,2,3,4) series.  This encapsulates what would
  ///        otherwise be a lot of repetitious code.
  ///
  /// \param order                  Order of the restraint (1 = positional, 2 = bond, ...)
  /// \param network_table_offsets  Pre-computed table of offsets for each restraint apparatus,
  ///                               relevant to the particular order of the restraints
  /// \param synthesis_index        Indices of parameters in each restraint apparatus into the
  ///                               condensed table of k(2,3) and r(1,2,3,4) parameter sets stored
  ///                               in the AtomGraphSynthesis (must be pre-allocated to be filled
  ///                               and returned)
  /// \param filtered_step_bounds   Initial and final bounds for the set of condensed restraint
  ///                               parameters.  No pre-allocation is needed.  This will be
  ///                               assembled and returned.
  /// \param filtered_init_keq      Initial stiffness constants for the set of condensed restraint
  ///                               parameters.  No pre-allocation is needed.  This will be
  ///                               assembled and returned.
  /// \param filtered_finl_keq      Mature stiffness constants for the set of condensed restraint
  ///                               parameters.  No pre-allocation is needed.  This will be
  ///                               assembled and returned.
  /// \param filtered_init_r        Initial displacement settings for the set of condensed
  ///                               restraint parameters.  No pre-allocation is needed.  This will
  ///                               be assembled and returned.
  /// \param filtered_finl_r        Mature displacement settings for the set of condensed restraint
  ///                               parameters.  No pre-allocation is needed.  This will be
  ///                               assembled and returned.
  int mapUniqueRestraintKRSeries(int order, const std::vector<int> &network_table_offsets,
                                 std::vector<int> *synthesis_index,
                                 std::vector<int2> *filtered_step_bounds,
                                 std::vector<double2> *filtered_init_keq,
                                 std::vector<double2> *filtered_finl_keq,
                                 std::vector<double4> *filtered_init_r,
                                 std::vector<double4> *filtered_finl_r);

  /// \brief Filter and condense the restraint networks in the same manner as valence terms were
  ///        condensed (each network comes from a RestraintApparatus object, and is just another
  ///        name to avoid repeating the type name as the name of an actual variable).
  void condenseRestraintNetworks();

  /// \brief Utility function for setting the limit arrays in the abstract of each ValenceWorkUnit.
  ///        This will set the x member of the appropriate tuple (identified by slot, read as a
  ///        literal integer, plus the appropriate stride times the valence work unit index) to
  ///        the current item count, the y member to the item count plus the quantity, and then
  ///        return the item count plus quantity padded with the warp size.
  ///
  /// \param item_counter   The current enumber of items across all previous work units
  /// \param vwu_counter    The index of the work unit
  /// \param slot           Type of instruction, the place in the work unit's abstract to set
  /// \param item_quantity  The number of items in the current work unit
  int setVwuAbstractLimits(int item_counter, int vwu_counter, VwuAbstractMap slot,
                           int item_quantity);
  
  /// \brief Construct valence work units for all systems and load their instructions into the
  ///        topology synthesis for availability on the GPU.
  ///
  /// \param vwu_atom_limit  The maximum number of atoms to assign to any one valence work unit
  void loadValenceWorkUnits(int vwu_atom_limit = maximum_valence_work_unit_atoms);

  /// \brief Construct the array of reduction work units.
  ///
  /// \param gpu  Details of the GPU to employ in the calculations
  void loadReductionWorkUnits(const GpuDetails &gpu = null_gpu);

  /// \brief Import atomic parameters for implicit solvent models based on one of the underlying
  ///        topologies.
  ///
  /// \param system_index  Index of the system to import (loop over all systems to complete the
  ///                      process)
  void importImplicitSolventAtomParameters(int system_index);

  /// \brief Impart the hard-wired "Neck" Generalized Born tables to the synthesis.
  ///
  /// Overloaded:
  ///   - Provide customized tables of neck GB parameters
  ///   - Use the hard-coded default tables
  /// 
  /// \param ngb_tab  Optional, customized table of neck Generalized Born parameters
  /// \{
  void setImplicitSolventNeckParameters(const NeckGeneralizedBornTable &ngb_tab);
  void setImplicitSolventNeckParameters();
  /// \}
};

} // namespace synthesis
} // namespace stormm

#include "atomgraph_synthesis.tpp"

#endif
