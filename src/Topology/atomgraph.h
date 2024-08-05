// -*-c++-*-
#ifndef STORMM_TOPOLOGY_H
#define STORMM_TOPOLOGY_H

#include <cmath>
#include <vector>
#include <string>
#include <time.h>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/hpc_bounds.h"
#include "Constants/symbol_values.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Parsing/citation.h"
#include "atomgraph_abstracts.h"
#include "atomgraph_enumerators.h"
#include "atomgraph_refinement.h"
#include "topology_limits.h"

namespace stormm {
namespace topology {

using constants::ExceptionResponse;
using constants::PrecisionModel;
using card::Hybrid;
using card::HybridTargetLevel;
using parse::Citation;
using symbols::amber_ancient_bioq;

/// \brief Use the lowest bit of a 32-bit float representing the range [1.0, 2.0) as the default
///        input for the smoothCharges function, ensuring that all charges in the range (-2e, +2e)
///        will be expressed in increments that a floating point number can represent exactly.
///        The hexfloat would read 0x1p-23.
constexpr double charge_precision_inc = 1.1920928955078125E-7;

/// \brief Default 1:4 non-bonded screening factors are taken from various major codes.
/// \{
constexpr double amber_default_elec14_screen = 1.2;
constexpr double amber_default_vdw14_screen = 2.0;
constexpr double charmm_default_elec14_screen = 1.0;
constexpr double charmm_default_vdw14_screen = 1.0;
constexpr double glycam_default_elec14_screen = 1.0;
constexpr double glycam_default_vdw14_screen = 1.0;
/// \}

/// \brief A struct to hold information relating to an Amber topology.  This struct's member
///        functions are limited to getters for its private data.  Because the primary means of
///        constructing a topology will be complex, i.e. meticulous parsing of a file or expansion
///        of an existing topology based on some new information, the constructors will need to be
///        buried within wrapper functions that perform such procedures.  This struct contains all
///        information contained within an Amber prmtop-format file, nothng more, as described on:
///
///        http://ambermd.org/FileFormats.php
///
///        The design is intended to be both performant as well as accessible to developers.
class AtomGraph {
public:

  /// \brief The blank constructor makes a blank AtomGraph, which is used by the general-
  ///        purpose file-based constructor to delegate initialization
  AtomGraph();

  /// \brief The general-purpose constructor for file-based topology creation
  ///
  /// Overloaded:
  ///   - Take just the file name, the exception response policy, and the kind of topology.
  ///     APpropriate default values will be applied for each kind of topology.
  ///   - In addition, take specific values for Coulomb's constant, the general 1:4 screening 
  ///     parameters, and charge discretization.
  ///   - Extract a subset of the atoms from an existing AtomGraph and create a new AtomGraph
  ///     from it.  All parameters and terms from the original AtomGraph will be carried over,
  ///     provided that the new AtomGraph contains all of the atoms that each term spans.
  ///     Exclusions from the original AtomGraph will be inferred from the bonding pattern.
  ///
  /// \param file_name                 Name of the source file
  /// \param policy                    The alert level to raise if a problem is encountered
  /// \param engine_format             Format of the topology file to read
  /// \param coulomb_constant_in       Value of Coulomb's constant, in kcal-A/mol-e^2
  /// \param default_elec14_screening  The 1:4 electrostatic scaling to apply in absense of any
  ///                                  other indications
  /// \param default_vdw14_screening   The 1:4 van-der Waals scaling to apply in absense of any
  ///                                  other indications
  /// \param charge_rounding_tol       The maximum tolerance at which to initiate charge rounding
  /// \param charge_discretization     Increment with which to discretize charges
  /// \{
  AtomGraph(const std::string &file_name, ExceptionResponse policy = ExceptionResponse::WARN,
            TopologyKind engine_format = TopologyKind::AMBER);

  AtomGraph(const std::string &file_name, ExceptionResponse policy, TopologyKind engine_format,
            double coulomb_constant_in, double default_elec14_screening,
            double default_vdw14_screening, double charge_rounding_tol,
            double charge_discretization);

  AtomGraph(const AtomGraph &original, const std::vector<int> &atom_subset,
            ExceptionResponse policy = ExceptionResponse::DIE);
  /// \}

  /// \brief The default destructor is adequate
  ~AtomGraph() = default;

  /// \brief The copy constructor takes a lot of tedium upon itself to prevent more complexity in
  ///        working with AtomGraph objects downstream.
  ///
  /// \param original  The original AtomGraph.  A deep copy of all data and repair of all pointers
  ///                  in the present object will be undertaken.
  AtomGraph(const AtomGraph &original);

  /// \brief Copy constructor
  ///
  /// \param other  Another AtomGraph to form the basis for re-assigning members of this one
  AtomGraph& operator=(const AtomGraph &other);

  /// \brief The move constructor makes prodigious use of std::move for each string and Hybrid
  ///        member variable.
  ///
  /// \param original  The original AtomGraph to move into this one.  The std::move function
  ///                  handles movement of the underlying Hybrid objects.
  AtomGraph(AtomGraph &&original);

  /// \brief Move constructor
  ///
  /// \param other  Another AtomGraph to form the basis for re-assigning members of this one  
  AtomGraph& operator=(AtomGraph &&other);

  /// \brief Load the  AtomGraph's various Hybrid objects with data held in temporary CPU
  ///        std::vectors.  See the function itself for details on each argument, but the arguments
  ///        generally follow the names of member variables in the AtomGraph itself.
  void loadHybridArrays(const std::vector<int> &tmp_desc,
                        const std::vector<int> &tmp_residue_limits,
                        const std::vector<int> &tmp_atom_struc_numbers,
                        const std::vector<int> &tmp_residue_numbers,
                        const std::vector<int> &tmp_molecule_limits,
                        const std::vector<int> &tmp_atomic_numbers,
                        const std::vector<int> &tmp_molecule_membership,
                        const std::vector<int> &tmp_mobile_atoms,
                        const std::vector<int> &tmp_molecule_contents,
                        const std::vector<int> &tmp_cmap_surf_dims,
                        const std::vector<int> &tmp_cmap_surf_bounds,
                        const std::vector<int> &tmp_charge_indices,
                        const std::vector<int> &tmp_lennard_jones_indices,
                        const std::vector<int> &tmp_inferred_14_i_atoms,
                        const std::vector<int> &tmp_inferred_14_j_atoms,
                        const std::vector<int> &tmp_inferred_14_param_idx,
                        const std::vector<int> &tmp_neck_gb_indices,
                        const std::vector<int> &tmp_tree_joining_info,
                        const std::vector<int> &tmp_last_rotator_info,
                        const std::vector<double> &tmp_charges,
                        const std::vector<double> &tmp_masses,
                        const std::vector<double> &tmp_ub_stiffnesses,
                        const std::vector<double> &tmp_ub_equilibria,
                        const std::vector<double> &tmp_charmm_impr_stiffnesses,
                        const std::vector<double> &tmp_charmm_impr_phase_angles,
                        const std::vector<double> &tmp_cmap_surfaces,
                        const std::vector<double> &tmp_bond_stiffnesses,
                        const std::vector<double> &tmp_bond_equilibria,
                        const std::vector<double> &tmp_angl_stiffnesses,
                        const std::vector<double> &tmp_angl_equilibria,
                        const std::vector<double> &tmp_dihe_amplitudes,
                        const std::vector<double> &tmp_dihe_periodicities,
                        const std::vector<double> &tmp_dihe_phase_angles,
                        const std::vector<double> &tmp_charge_parameters,
                        const std::vector<double> &tmp_lj_a_values,
                        const std::vector<double> &tmp_lj_b_values,
                        const std::vector<double> &tmp_lj_c_values,
                        const std::vector<double> &tmp_lj_14_a_values,
                        const std::vector<double> &tmp_lj_14_b_values,
                        const std::vector<double> &tmp_lj_14_c_values,
                        const std::vector<double> &tmp_atomic_pb_radii,
                        const std::vector<double> &tmp_gb_screening_factors,
                        const std::vector<double> &tmp_gb_coef,
                        const std::vector<double> &tmp_solty_info,
                        const std::vector<double> &tmp_hbond_a_values,
                        const std::vector<double> &tmp_hbond_b_values,
                        const std::vector<double> &tmp_hbond_cutoffs,
                        const std::vector<char4> &tmp_atom_names,
                        const std::vector<char4> &tmp_atom_types,
                        const std::vector<char4> &tmp_residue_names,
                        const std::vector<char4> &tmp_tree_symbols,
                        const CmapAccessories &cmap_table,
                        const CondensedExclusions &cond_excl,
                        const BasicValenceTable &basic_vtable,
                        const CharmmValenceTable &charmm_vtable,
                        const AttenuationParameterSet &attn_parm,
                        const VirtualSiteTable &vsite_table,
                        const Map1234 &all_nb_excl, const ConstraintTable &cnst_table);
  
  /// \brief Build an AtomGraph form a file.  This is called by the general-purpose constructor
  ///        or also by the developer after instantiating an empty object.
  ///
  /// \param file_name  Name of the source file
  /// \param policy     Indicates the alert level to raise if a problem is encountered
  void buildFromPrmtop(const std::string &file_name,
                       const ExceptionResponse policy = ExceptionResponse::WARN,
                       double coulomb_constant_in = amber_ancient_bioq,
                       double default_elec14_screening = amber_default_elec14_screen,
                       double default_vdw14_screening = amber_default_vdw14_screen,
                       double charge_rounding_tol = 0.001,
                       double charge_discretization = charge_precision_inc);

  // Most getter functions (getThisStuff) will return the corresponding member variable.  Usually
  // a variable with a single value is an int, but it could also be an enumeration or a string.
  // For Hybrid object member variables, getThisHybridObjectStuff() will return the entirety of
  // the POINTER-kind object's host_data vector, up to its stated length.

  /// \brief Get the title of the topology (may be blank).
  std::string getTitle() const;
  
  /// \brief Get the file name where a topology originated (may be blank, indicating that the
  ///        topology was produced by some other means).
  std::string getFileName() const;

  /// \brief Get the number of atoms, using the topology's dedicated private variable rather than
  ///        a list of input dimensions that was mostly stored for developers most familiar with
  ///        Amber.
  int getAtomCount() const;

  /// \brief Get the number of residues.
  int getResidueCount() const;

  /// \brief Get the number of separate molecules in the system
  int getMoleculeCount() const;
  
  /// \brief Get the number of organic compounds in the system
  int getOrganicCompoundsCount() const;

  /// \brief Get the number of separate molecules in the system
  int getLargestResidueSize() const;

  /// \brief Get the final solute residue, indexed according to C/C++ array indexing conventions
  int getLastSoluteResidue() const;

  /// \brief Get the final solute atom, indexed according to C/C++ array indexing conventions
  int getLastSoluteAtom() const;

  /// \brief Get the first solvent molecule, indexed according to C/C++ array indexing conventions
  int getFirstSolventMolecule() const;

  /// \brief Get the first solvent atom, indexed according to C/C++ array indexing conventions
  int getFirstSolventAtom() const;

  /// \brief Get the largest molecule's size
  int getLargestMoleculeSize() const;

  /// \brief Get the total mass of all atoms in the topology (this is a computation, not a stored
  ///        value).
  double getTotalMass() const;

  /// \brief Get the number of degrees of freedom, without consideration to geometric constraints.
  int getDegreesOfFreedom() const;

  /// \brief Get the number of degrees of freedom after geometric constraints have been applied.
  int getConstrainedDegreesOfFreedom() const;
  
  /// \brief Get a descriptor from within the array of topology descriptors.  If this topology were
  ///        read from an Amber-format file, all descriptors are taken from the preamble in its
  ///        first ~15 lines.
  ///
  /// Overloaded:
  ///   - Take one of the longer, clarified enumerations from the TopologyDescriptor class
  ///   - Take one of the shorter enumerations pulled from the Amber sander program and the online
  ///     documentation (http://ambermd.org/FileFormats.php)
  ///
  /// \param choice   Index of the descriptor in the list.  Call this function using the enum
  ///                 classes TopologyDescriptors or SanderDescriptors (see above) for easy,
  ///                 annotated access.   
  /// \{
  int getDescriptor(TopologyDescriptor choice) const;
  int getDescriptor(SanderDescriptor choice) const;
  /// \}

  /// \brief Get residue limits for the beginning and ends of up to N residues
  ///
  /// Overloaded:
  ///   - Get all residue limits as a std::vector<int>
  ///   - Get the initial and final atom indices of a specific residue as an int2
  ///
  /// \param index  Index of a specific residue
  /// \{
  const Hybrid<int>& getResidueLimits() const;
  int2 getResidueLimits(int index) const;
  /// \}

  /// \brief  Get the index (topology index, not naturaWl / structure-informed residue number) of a
  ///         particular residue.
  ///
  /// Overloaded:
  ///   - Get all residue indices as a std::vector<int>
  ///   - Get the residue index of a particular atom
  ///
  /// \param atom_index  Index of the atom in question
  /// \{
  std::vector<int> getResidueIndex() const;
  int getResidueIndex(int atom_index) const;
  /// \}

  /// \brief Get the (arbitrarily defined) atom numbers for atoms in the system.  For example, a
  ///        PDB file could ahve its own atom numbering scheme.  This provides a way to transcribe
  ///        that into the AtomGraph object.
  ///
  /// Overloaded:
  ///   - Get the structural atom numbers of all atoms in the system, as a std::vector<int>
  ///   - Get structural atom numbers of a stretch of atoms in the system, as a std::vector<int>
  ///   - Get the structural atom number of a specific atom
  ///
  /// \param low_index   Index of the first atom for which to read a residue number
  /// \param high_index  Index of the last atom for which to read a residue number
  /// \param index       Index of a specific atom
  /// \{
  std::vector<int> getStructuralAtomNumber() const;
  std::vector<int> getStructuralAtomNumber(int low_index, int high_index) const;
  int getStructuralAtomNumber(int index) const;
  /// \}

  /// \brief Get the residue numbers for atoms in the system
  ///
  /// Overloaded:
  ///   - Get the residue numbers of all atoms in the system, as a std::vector<int>
  ///   - Get the residue numbers of a stretch of atoms in the system, as a std::vector<int>
  ///   - Get the residue number of a specific atom
  ///
  /// \param low_index   Index of the first atom for which to read a residue number
  /// \param high_index  Index of the last atom for which to read a residue number
  /// \param index       Index of a specific atom
  /// \{
  std::vector<int> getResidueNumber() const;
  std::vector<int> getResidueNumber(int low_index, int high_index) const;
  int getResidueNumber(int index) const;
  /// \}

  /// \brief Get the molecule limits, indices into molecule_contents, for every molecule in the
  ///        system.  This allows that molecules not be contiguous in the topology atom indexing,
  ///        but will be recoverable by molecule_contents.
  ///
  /// Overloaded:
  ///   - Get all molecule limits as a std::vector<int>
  ///   - Get the initial and final indices of a specific molecule in molecule_contents as an int2
  ///
  /// \param index  Index of a specific molecule
  /// \{
  std::vector<int> getMoleculeLimits() const;
  int2 getMoleculeLimits(int index) const;
  /// \} 

  /// \brief Get the sizes of all molecules in the system.  Rather redundant with the molecule
  ///        limits being accessible, but convenient.
  ///
  /// Overloaded:
  ///   - Get the sizes of all molecules as a std::vector<int>
  ///   - Get the size of a single molecule
  ///
  /// \param index  Index of a specific residue
  std::vector<int> getParticlesPerMolecule() const;
  int getParticlesPerMolecule(int index) const;
  /// \}

  /// \brief Get the atomic numbers for atoms in the system
  ///
  /// Overloaded:
  ///   - Get the atomic numbers of all atoms in the system, as a std::vector<int>
  ///   - Get the atomic numbers of a stretch of atoms in the system, as a std::vector<int>
  ///   - Get the atomic number of a specific atom
  ///
  /// \param low_index   Index of the first atom for which to read a Z number
  /// \param high_index  Index of the last atom for which to read a Z number
  /// \param index       Index of a specific atom
  /// \{
  std::vector<int> getAtomicNumber() const;
  std::vector<int> getAtomicNumber(int low_index, int high_index) const;
  int getAtomicNumber(int index) const;
  /// \}

  /// \brief Get the mobile atom mask for the system, as a boolean vector translated from the
  ///        bitmask stored internally.
  ///
  /// Overloaded:
  ///   - Get the mobility of all atoms in the system, as a std::vector<bool>
  ///   - Get the mobility of a stretch of atoms in the system
  ///   - Get the mobility of a particular atom
  ///
  /// \param low_index   Index of the first atom for which to read the mobility
  /// \param high_index  Index of the last atom for which to read the mobility
  /// \param index       Index of a specific atom
  /// \{
  std::vector<bool> getAtomMobility() const;
  std::vector<bool> getAtomMobility(int low_index, int high_index) const;
  bool getAtomMobility(int index) const;
  /// \}

  /// \brief Get the raw mobility bitmask over a particular range.
  ///
  /// Overloaded:
  ///   - Get the mobility of all atoms in the system, as a std::vector<uint>
  ///   - Get the mobility of a stretch of atoms in the system
  ///
  /// \param low_index   Index of the first atom for which to read the mobility
  /// \param high_index  Index of the last atom for which to read the mobility
  /// \{
  std::vector<uint> getAtomMobilityMask() const;
  std::vector<uint> getAtomMobilityMask(int low_index, int high_index) const;
  /// \}
  
  /// \brief Change an atom's mobility in the topology.
  ///
  /// Overloaded:
  ///   - Change all atoms to become mobile or not, or toggle mobility.
  ///   - Change a specific stretch of atoms to become mobile or not, or toggle mobility.
  ///   - Make a specific atom mobile or not, or toggle mobility.
  ///   - Make a scattered selection of atoms mask mobile or not, or toggle the mobility.
  ///
  /// \param movement    The desired atom mobility (OFF, ON, or TOGGLE)
  /// \param low_index   The start of a stretch of atoms
  /// \param high_index  The end of a stretch of atoms
  /// \param index       A specific atom index
  /// \{
  void modifyAtomMobility(MobilitySetting movement);
  void modifyAtomMobility(int low_index, int high_index, MobilitySetting movement);
  void modifyAtomMobility(int index, MobilitySetting movement);
  void modifyAtomMobility(const std::vector<int> &mask, MobilitySetting movement);
  /// \}

  /// \brief Get molecule membership for atoms within the system
  ///
  /// Overloaded:
  ///   - Get the molecular membership of all atoms in the system, indexed from 0
  ///   - Get the molecular membership a stretch of atoms in the system
  ///   - Get the molecule index that a particular atom is part of
  ///
  /// \param low_index   Index of the first atom for which to read the mobility
  /// \param high_index  Index of the last atom for which to read the mobility
  /// \param index       Index of a specific atom
  /// \{
  std::vector<int> getMoleculeMembership() const;
  std::vector<int> getMoleculeMembership(int low_index, int high_index) const;
  int getMoleculeMembership(int index) const;
  /// \}

  /// \brief Get the contents of a one or more molecules in the system
  ///
  /// Overloaded:
  ///   - Get the contents of all molecules
  ///   - Get the contents of one particular molecule
  ///
  /// \param index  Index of a specific molecule
  /// \{
  std::vector<int> getMoleculeContents() const;
  std::vector<int> getMoleculeContents(int index) const;
  /// \}

  /// \brief Get the atomic partial charges of atoms in the system.
  ///
  /// Overloaded:
  ///   - Get partial charges of all atoms in the system, indexed from 0
  ///   - Get partial charges of a stretch of atoms in the system
  ///   - Get the partial charge of a particular atom
  ///
  /// \param low_index   The start of a stretch of atoms
  /// \param high_index  The end of a stretch of atoms
  /// \param index       A specific atom index
  /// \{
  template <typename T> std::vector<T> getPartialCharge() const;
  template <typename T> std::vector<T> getPartialCharge(int low_index, int high_index) const;
  template <typename T> T getPartialCharge(int index) const;
  /// \}

  /// \brief Get the masses of atoms in the system.  The developer may also choose to take the
  ///        pre-computed inverse masses (this can be useful when dealing with extra points, as
  ///        inverses of zero-mass particles must also be set to zero by a conditional statement).
  ///
  /// Overloaded:
  ///   - Get masses of all atoms in the system, indexed from 0
  ///   - Get masses of a stretch of atoms in the system
  ///   - Get the mass of a particular atom
  ///
  /// \param rep         The representation of the mass to take: ORDINARY or INVERSE
  /// \param low_index   The start of a stretch of atoms
  /// \param high_index  The end of a stretch of atoms
  /// \param index       A specific atom index
  /// \{
  template <typename T> std::vector<T> getAtomicMass(MassForm rep = MassForm::ORDINARY) const;
  template <typename T> std::vector<T> getAtomicMass(int low_index, int high_index,
                                                     MassForm rep = MassForm::ORDINARY) const;
  template <typename T> T getAtomicMass(int index, MassForm rep = MassForm::ORDINARY) const;
  /// \}

  /// \brief Get the names of atoms in the system.
  ///
  /// Overloaded:
  ///   - Get names of all atoms in the system, indexed from 0
  ///   - Get names of a stretch of atoms in the system
  ///   - Get the name of a particular atom
  ///
  /// \param low_index   The start of a stretch of atoms
  /// \param high_index  The end of a stretch of atoms
  /// \param index       A specific atom index
  /// \{
  const Hybrid<char4>& getAtomName() const;
  std::vector<char4> getAtomName(int low_index, int high_index) const;
  char4 getAtomName(int index) const;
  /// \}

  /// \brief Get the types of atoms in the system.
  ///
  /// Overloaded:
  ///   - Get types of all atoms in the system, indexed from 0
  ///   - Get types of a stretch of atoms in the system
  ///   - Get the type of a particular atom
  ///
  /// \param low_index   The start of a stretch of atoms
  /// \param high_index  The end of a stretch of atoms
  /// \param index       A specific atom index
  /// \{
  const Hybrid<char4>& getAtomType() const;
  std::vector<char4> getAtomType(int low_index, int high_index) const;
  char4 getAtomType(int index) const;
  /// \}

  /// \brief Get a table with the names of all unique atom types, arranged according to their
  ///        Lennard-Jones indices as they appear in the topology (in ascending order).
  std::vector<std::vector<char4>> getAtomTypeNameTable() const;
  
  /// \brief Get the names of residues in the system.
  ///
  /// Overloaded:
  ///   - Get names of all residues in the system, indexed from 0
  ///   - Get the name of a particular residue
  ///
  /// \param index       A specific residue index
  /// \{
  const Hybrid<char4>& getResidueName() const;
  char4 getResidueName(int index) const;
  /// \}

  /// \brief Get the number of Urey-Bradley terms in the system
  int getUreyBradleyTermCount() const;

  /// \brief Get the number of CHARMM improper terms
  int getCharmmImprTermCount() const;

  /// \brief Get the number of CMAP terms (the number of CMAP parameter surfaces is different)
  int getCmapTermCount() const;

  /// \brief Return the number of unique Urey-Bradley parameter sets
  int getUreyBradleyParameterCount() const;

  /// \brief Return the number of unique CHARMM improper parameter sets
  int getCharmmImprParameterCount() const;

  /// \brief Return the number of unique CMAP surfaces
  int getCmapSurfaceCount() const;

  /// \brief Return the dimension of a particular CMAP surface (all CMAPs are assumed to be square)
  ///
  /// \param index  The CMAP surface of interest (along the list of available CMAP surfaces, not
  ///               the list of CMAP terms)
  int getCmapDimension(int index) const;

  /// \brief Collect the atom indices and actual parameters to implement one of the system's
  ///        Urey-Bradley terms.  The special-purpose, unguarded struct that holds the result also
  ///        provides the original parameter index, if of interest.  This is for computing the
  ///        potential energy and forces in the most straightforward manner possible, albeit slow.
  ///
  /// \param index  The Urey-Bradley term to query, indexed according to the list of Urey-Bradley
  ///               terms (not the general atom list)
  template <typename T> UreyBradleyTerm<T> getUreyBradleyTerm(int index) const;

  /// \brief Collect the atom indices and actual parameters to implement one of the system's
  ///        CHARMM improper terms.  The special-pupose, unguarded struct that holds the result
  ///        also provides the original parameter index, if of interest.  This is for computing the
  ///        potential energy and forces in the most straightforward manner possible, albeit slow.
  ///
  /// \param index  The CHARMM improper term to query, indexed according to the list of CHARMM
  ///               improper terms (not the general atom list)
  template <typename T> CharmmImprTerm<T> getCharmmImprTerm(int index) const; 

  /// \brief Collect the atom indices and a pointer to the surface to implement one of the system's
  ///        CMAP terms.  The special-purpose, unguarded struct that holds the result also provides
  ///        the original parameter index, if of interest.
  ///
  /// \param index  The CMAP term to query, indexed according to the list of CMAP terms (not the
  ///               general atom list)
  template <typename T> CmapTerm<T> getCmapTerm(int index) const;

  /// \brief Get the number of bond stretching terms in the system
  int getBondTermCount() const;

  /// \brief Get the number of angle bending terms in the system
  int getAngleTermCount() const;

  /// \brief Get the number of dihedral torsion terms in the system
  int getDihedralTermCount() const;

  /// \brief Return the number of unique bond stretching parameter sets
  int getBondParameterCount() const;

  /// \brief Return the number of unique bond angle bending parameter sets
  int getAngleParameterCount() const;

  /// \brief Return the number of unique dihedral parameter sets
  int getDihedralParameterCount() const;

  /// \brief Collect the atom indices and actual parameters for implement one of the system's
  ///        harmonic bond stretching terms.  
  ///
  /// \param index  The bond term to query, indexed according to the list of bond terms
  template <typename T> BondTerm<T> getBondTerm(int index) const;

  /// \brief Collect the atom indices and actual parameters for implement one of the system's
  ///        harmonic bond angle bending terms.
  ///
  /// \param index  The bond angle term to query, indexed according to the list of bond angle terms
  template <typename T> AngleTerm<T> getAngleTerm(int index) const;

  /// \brief Collect the atom indices and actual parameters for implement one of the system's
  ///        cosine-based dihedral torsion terms.
  ///
  /// \param index  The dihedral term to query, indexed according to the list of dihedral terms
  template <typename T> DihedralTerm<T> getDihedralTerm(int index) const;

  /// \brief Get the number of virtual sites in the system
  int getVirtualSiteCount() const;

  /// \brief Get the number of unique virtual site frames in the system.  For example, a system of
  ///        1024 TIP-4P waters would have one unique frame, applied to 1024 different groups of
  ///        atoms.
  int getVirtualSiteParameterSetCount() const;
  
  /// \brief Get a list of the general topological indices of one or more virtual sites.
  ///
  /// Overloaded:
  ///   - Get a list of all virtual sites, indexed by their order in the general topology
  ///   - Get virtual sites within a range, limits given by ordering in the general topology,
  ///     results given as indices into the general topology (part x) plus indices into the
  ///     array of virtual sites (part y).  For the first case, "get all virtual sites," the
  ///     index in the returned array of general topology atom indices serves as an implicit
  ///     indicator of the index into the list of virtual sites.
  ///   - Get the topological index of one virtual site out of the list of virtual sites
  ///
  /// \param low_index   The start of a stretch of atoms in the general topology
  /// \param high_index  The end of a stretch of atoms in the general topology
  /// \param index       A specific virtual site index
  /// \{
  std::vector<int> findVirtualSites() const;
  std::vector<int2> findVirtualSites(int low_index, int high_index) const;
  int findVirtualSites(int index) const;
  /// \}

  /// \brief Get the virtual site index of an atom index.  If the atom index is not a virtual site,
  ///        this function will throw an error.
  ///
  /// \param atom_index  Index of the atom in question
  /// \param policy      Action to take if the atom turns out not to be a virtual site at all
  int getVirtualSiteIndex(int atom_index, ExceptionResponse policy = ExceptionResponse::DIE) const;
  
  /// \brief Get the frame type of a particular virtual site.
  ///
  /// \param index  The virtual site of interest, indexed according to its place in the list of
  ///               virtual sites (not the general topology atom list)
  VirtualSiteKind getVirtualSiteFrameType(int index) const;

  /// \brief Get the general topology index of one of a virtual site's frame atoms
  ///
  /// \param index  The virtual site of interest, indexed according to its place in the list of
  ///               virtual sites (not the general topology atom list)
  /// \param nfrm   The frame atom to pull out (values of 1:4 for one, two, three, or four-atom
  ///               frames are acceptable)
  int getVirtualSiteFrameAtom(int index, int nfrm) const;

  /// \brief Get the dimensions (could be length in A, could be scaling factor for a distance or
  ///        cross product) of one virtual site's frame.
  ///
  /// \param index     The virtual site of interest, indexed according to its place in the list of
  ///                  virtual sites (not the general topology atom list)
  /// \param ndim      The frame dimension to pull out (values of 1-3 for frames with one, two, or
  ///                  three details are acceptable)
  template <typename T> T getVirtualSiteFrameDimension(int index, int ndim);

  /// \brief Get the number of charge types in the system
  int getChargeTypeCount() const;
  
  /// \brief Get the number of atom types in the system.
  int getLJTypeCount() const;

  /// \brief Get the number of exclusions in the entire topology.
  int getTotalExclusions() const;

  /// \brief Get the simulation cell type (isolated boundary conditions, rectilinear, triclinic)
  UnitCellType getUnitCellType() const;

  /// \brief Get the implicit solvent model type, i.e. some flavor of GB or NONE.
  ImplicitSolventModel getImplicitSolventModel() const;

  /// \brief Get the dielectric constant (part of the implicit solvent setup)
  double getDielectricConstant() const;

  /// \brief Get the salt concentration (part of the implicit solvent setup)
  double getSaltConcentration() const;

  /// \brief Get the fundamental coulomb constant that this system uses in electrostatic
  ///        calculations.
  double getCoulombConstant() const;

  /// \brief Get the general 1:4 screening factor on electrostatic interactions.
  double getElectrostatic14Screening() const;

  /// \brief Get the general 1:4 screening factor on van-der Waals interactions.
  double getVanDerWaals14Screening() const;

  /// \brief Get the PB Radii set
  std::string getPBRadiiSet() const;

  /// \brief Get the water residue name
  std::string getWaterResidueName() const;
  
  /// \brief Get the number of rigid waters in the system.  This will check through the SETTLE
  ///        groups and confirm the number that are, by their chemistry, water.
  int getRigidWaterCount() const;

  /// \brief Get the number of SETTLE groups in the system
  int getSettleGroupCount() const;

  /// \brief Get the number of bond constraints
  int getBondConstraintCount() const;

  /// \brief Get the total number of constrained groups
  int getConstraintGroupCount() const;

  /// \brief Get the atoms of a single constraint group
  ///
  /// \param index  The index of the constraint group within the topology
  std::vector<int> getConstraintGroupAtoms(int index) const;
  
  /// \brief Get the total size of the constrained group atoms list
  int getConstraintGroupTotalSize() const;
  
  /// \brief Get the number of non-rigid particles in the system
  int getNonrigidParticleCount() const;
  
  /// \brief Get the charge type index of atoms in the system.
  ///
  /// Overloaded:
  ///   - Get charge indices of all atoms in the system, indexed from 0
  ///   - Get charge indices of a stretch of atoms in the system
  ///   - Get the charge index of a particular atom
  ///
  /// \param low_index   The start of a stretch of atoms
  /// \param high_index  The end of a stretch of atoms
  /// \param index       A specific atom index
  /// \{
  std::vector<int> getChargeIndex() const;
  std::vector<int> getChargeIndex(int low_index, int high_index) const;
  int getChargeIndex(int index) const;
  /// \}

  /// \brief Get the Lennard-Jones type index of atoms in the system.
  ///
  /// Overloaded:
  ///   - Get Lennard-Jones indices of all atoms in the system, indexed from 0
  ///   - Get Lennard-Jones indices of a stretch of atoms in the system
  ///   - Get the Lennard-Jones index of a particular atom
  ///
  /// \param low_index   The start of a stretch of atoms
  /// \param high_index  The end of a stretch of atoms
  /// \param index       A specific atom index
  /// \{
  std::vector<int> getLennardJonesIndex() const;
  std::vector<int> getLennardJonesIndex(int low_index, int high_index) const;
  int getLennardJonesIndex(int index) const;
  /// \}

  /// \brief Get the exclusion list for a particular atom.  This will concatenate 1:1, 1:2, 1:3,
  ///        and 1:4 exclusions, without regard to which is which.  To get a specific type of
  ///        exclusion, use one of the getNonbondedXXExclusions() functions.  All functions of
  ///        this sort will return a complete list such that if exclusions for all atoms are
  ///        considered in sequence, all exclusions will be double-counted.
  ///
  /// \param index  The atom of interest
  std::vector<int> getAtomExclusions(int index) const;

  /// \brief Get the 1:1 exclusions for a particular atom.  This list will only be populated if
  ///        the atom is a virtual site or has virtual sites that claim it as their parent.
  ///
  /// \param index  The atom of interest
  std::vector<int> getNonbonded11Exclusions(int index) const;

  /// \brief Get all 1:2 exclusions for a particular atom.
  ///
  /// \param index  The atom of interest
  std::vector<int> getNonbonded12Exclusions(int index) const;

  /// \brief Get all 1:3 exclusions for a particular atom.
  ///
  /// \param index  The atom of interest
  std::vector<int> getNonbonded13Exclusions(int index) const;

  /// \brief Get all 1:4 exclusions for a particular atom.
  ///
  /// \param index  The atom of interest
  std::vector<int> getNonbonded14Exclusions(int index) const;

  /// \brief Get the charge parameter for a specific atom.  The quantity is derived from tables at
  ///        the specified level of precision (float or double), and equals the value obtained from
  ///        getPartialCharge<T>(index), although that function reaches into a separate array with
  ///        a great deal of redundant information.
  template <typename T> T getChargeParameter(int index) const;

  /// \brief Get the Lennard-Jones interaction sigma.  The quantity is derived from tables at the
  ///        specified level of precision (float or double), but computed in double precision
  ///        based on that information before it is returned in the template data type.
  ///
  /// Overloaded:
  ///   - Get a single atom's self-interaction sigma
  ///   - Get the sigma for a pair of atoms (may not obey combining rules, e.g. CHARMM NB-Fix)
  ///   - Get a table for the self-interaction sigma values of all types present in the topology
  ///
  /// \param index_a   Index for the first atom
  /// \param index_b   Index for the second atom
  /// \{
  template <typename T> T getLennardJonesSigma(int index_a, int index_b) const;
  template <typename T> T getLennardJonesSigma(int index_a) const;
  template <typename T> std::vector<T> getLennardJonesSigma() const;
  /// \}

  /// \brief Get the Lennard-Jones interaction epsilon value.  The quantity is derived from tables
  ///        at the specified level of precision (float or double), but computed in double
  ///        precision based on that information before it is returned in the template data type.
  ///
  /// Overloaded:
  ///   - Get a single atom's self-interaction epsilon
  ///   - Get the epsilon for a pair of atoms (may not obey combining rules, e.g. CHARMM NB-Fix)
  ///   - Get a table for the self-interaction epsilon values of all types present in the topology
  ///
  /// \param index_a   Index for the first atom
  /// \param index_b   Index for the second atom
  /// \{
  template <typename T> T getLennardJonesEpsilon(int index_a, int index_b) const;
  template <typename T> T getLennardJonesEpsilon(int index_a) const;
  template <typename T> std::vector<T> getLennardJonesEpsilon() const;
  /// \}

  /// \brief Get the Poisson-Boltzmann radius of a particular atom
  ///
  /// Overloaded:
  ///   - Get Poisson-Boltzmann radii of all atoms in the system, indexed from 0
  ///   - Get Poisson-Boltzmann radii for a stretch of atoms in the system
  ///   - Get the Poisson-Boltzmann radius of a particular atom
  ///
  /// \param low_index   The start of a stretch of atoms
  /// \param high_index  The end of a stretch of atoms
  /// \param index       A specific atom index
  /// \{
  template <typename T> std::vector<T> getAtomPBRadius() const;
  template <typename T> std::vector<T> getAtomPBRadius(int low_index, int high_index) const;
  template <typename T> T getAtomPBRadius(int index) const;
  /// \}

  /// \brief Get GB screening factors for one or more atoms.
  ///
  /// Overloaded:
  ///   - Get GB screening factors of all atoms in the system, indexed from 0
  ///   - Get GB screening factors for a stretch of atoms in the system
  ///   - Get the GB screening factor of a particular atom
  ///
  /// \param low_index   The start of a stretch of atoms
  /// \param high_index  The end of a stretch of atoms
  /// \param index       A specific atom index
  /// \{
  template <typename T> std::vector<T> getGBScreeningFactor() const;
  template <typename T> std::vector<T> getGBScreeningFactor(int low_index, int high_index) const;
  template <typename T> T getGBScreeningFactor(int index) const;
  /// \}

  /// \brief Match a long atom name and return its key in the standard list of char4 atom names
  ///
  /// \param query  The extended atom name to match
  std::vector<char4> matchOverflowAtomName(const std::string &query) const;

  /// \brief Match a long atom type and return its key in the standard list of char4 atom types
  ///
  /// \param query  The extended atom type name to match
  std::vector<char4> matchOverflowAtomType(const std::string &query) const;

  /// \brief Match a long residue name and return its key among the standard char4 residue names
  ///
  /// \param query  The extended residue name to match
  std::vector<char4> matchOverflowResidueName(const std::string &query) const;

  /// \brief Get an abstract for calculating valence terms in double precision.
  ///
  /// \param tier  Indicator of whether pointers should target the CPU or the GPU layer
  ValenceKit<double>
  getDoublePrecisionValenceKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
  /// \brief Get a single-precision abstract for calculating valence terms.
  ///
  /// \param tier  Indicator of whether pointers should target the CPU or the GPU layer
  ValenceKit<float>
  getSinglePrecisionValenceKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a double-precision abstract for calculating non-bonded energy and forces.
  ///
  /// \param tier  Indicator of whether pointers should target the CPU or the GPU layer
  NonbondedKit<double>
  getDoublePrecisionNonbondedKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a single-precision abstract for calculating non-bonded energy and forces.
  ///
  /// \param tier  Indicator of whether pointers should target the CPU or the GPU layer
  NonbondedKit<float>
  getSinglePrecisionNonbondedKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a double-precision abstract for calculating implicit solvent energy and forces.
  ///
  /// \param tier  Indicator of whether pointers should target the CPU or the GPU layer
  ImplicitSolventKit<double>
  getDoublePrecisionImplicitSolventKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a single-precision abstract for calculating implicit solvent energy and forces.
  ///
  /// \param tier  Indicator of whether pointers should target the CPU or the GPU layer
  ImplicitSolventKit<float>
  getSinglePrecisionImplicitSolventKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get an abstract of the chemical details of the system, such as atom name and type or
  ///        residue names and limits.
  ///
  /// \param tier  Indicator of whether pointers should target the CPU or the GPU layer
  ChemicalDetailsKit
  getChemicalDetailsKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a double-precision abstract for placing virtual sites and transmitting their
  ///        forces to frame atoms with mass.
  ///
  /// \param tier  Indicator of whether pointers should target the CPU or the GPU layer
  VirtualSiteKit<double>
  getDoublePrecisionVirtualSiteKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a single-precision abstract for placing virtual sites and transmitting their
  ///        forces to frame atoms with mass.
  ///
  /// \param tier  Indicator of whether pointers should target the CPU or the GPU layer
  VirtualSiteKit<float>
  getSinglePrecisionVirtualSiteKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get an abstract for managing constraints based on double-precision parameters.
  ///
  /// \param tier  Indicator of whether pointers should target the CPU or the GPU layer
  ConstraintKit<double>
  getDoublePrecisionConstraintKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get an abstract for managing constraints based on single-precision parameters.
  ///        While the overall implementation for constraints cannot be implemented entirely in
  ///        single precision, having the underlying parameters in this form is not, in principle,
  ///        a drvier of energy drift.  It will simply set the system into a slightly different
  ///        geometry, which is probably still closer to the double-precision result than many
  ///        cycles of SHAKE or RATTLE would get it anyway.
  ///
  /// \param tier  Indicator of whether pointers should target the CPU or the GPU layer
  ConstraintKit<float>
  getSinglePrecisionConstraintKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  /// \brief Get a double-precision abstract for calculating valence terms, visible on the GPU
  ///        device with data stored on the CPU host.
  ValenceKit<double> getDeviceViewToHostDPValenceKit() const;

  /// \brief Get a single-precision abstract for calculating valence terms, visible on the GPU
  ///        device with data stored on the CPU host.
  ValenceKit<float> getDeviceViewToHostSPValenceKit() const;

  /// \brief Get a double-precision abstract for calculating non-bonded terms, visible on the GPU
  ///        device with data stored on the CPU host.
  NonbondedKit<double> getDeviceViewToHostDPNonbondedKit() const;

  /// \brief Get a single-precision abstract for calculating non-bonded terms, visible on the GPU
  ///        device with data stored on the CPU host.
  NonbondedKit<float> getDeviceViewToHostSPNonbondedKit() const;

  /// \brief Get a double-precision abstract for the system's implicit solvent details, visible on
  ///        the GPU device with data stored on the CPU host.
  ImplicitSolventKit<double> getDeviceViewToHostDPImplicitSolventKit() const;

  /// \brief Get a single-precision abstract for the system's implicit solvent details, visible on
  ///        the GPU device with data stored on the CPU host.
  ImplicitSolventKit<float> getDeviceViewToHostSPImplicitSolventKit() const;

  /// \brief Get an abstract for the system's chemical details, visible on the GPU device with data
  ///        stored on the CPU host.
  ChemicalDetailsKit getDeviceViewToHostChemicalDetailsKit() const;

  /// \brief Get a double-precision abstract for placing or transmuting forces from virtual sites,
  ///        visible on the GPU device with data stored on the CPU host.
  VirtualSiteKit<double> getDeviceViewToHostDPVirtualSiteKit() const;

  /// \brief Get a single-precision abstract for placing or transmuting forces from virtual sites,
  ///        visible on the GPU device with data stored on the CPU host.
  VirtualSiteKit<float> getDeviceViewToHostSPVirtualSiteKit() const;

  /// \brief Get a double-precision abstract for implementing constraints, visible on the GPU
  ///        device with data stored on the CPU host.
  ConstraintKit<double> getDeviceViewToHostDPConstraintKit() const;

  /// \brief Get a single-precision abstract for implementing constraints, visible on the GPU
  ///        device with data stored on the CPU host.
  ConstraintKit<float> getDeviceViewToHostSPConstraintKit() const;
#  endif
#endif

  /// \brief Get a pointer to the object itself.
  ///
  /// Overloaded:
  ///   - Get a const pointer to a const AtomGraph
  ///   - Get a non-const pointer to a mutable AtomGraph
  /// \{
  const AtomGraph* getSelfPointer() const;
  AtomGraph* getSelfPointer();
  /// \}
  
#ifdef STORMM_USE_HPC
  /// \brief Push all data to the device
  void upload();

  /// \brief Pull all data back to the host
  void download();
#endif

  /// \brief Set the source file for the topology.  Primary use is to impart a destination file
  ///        for a topology to be written to disk, but can also be used to "disguise" one topology
  ///        from another even if both are from the same original file.
  void setSource(const std::string &new_source);

  /// \brief Alter the parameters for a bond parameter set.  This will alter the parameters in the
  ///        topology and all bond terms that use the parameters will then be altered.  This
  ///        function does not add or remove bond terms or parameter sets, nor does it change the
  ///        mapping of bond terms to parameter sets.
  ///
  /// Overloaded:
  ///   - Accept the atom types for each atom in the bond
  ///   - Accept the parameter set index for the bond
  ///
  /// \param new_keq   The new stiffness constant to apply
  /// \param new_leq   The new equilibrium constant to apply
  /// \param set_keq   Whether to change the stiffness constant
  /// \param set_leq   Whether to change the equilibrium
  /// \param parm_idx  Index of the bond parameter set from within the topology's table (one cannot
  ///                  depend on this index to be consistent across topologies of different systems
  ///                  making use of the same force field)
  /// \param type_a    Atom type of the first atom in the bond (this will be consistent across
  ///                  multiple topologies using the same force field)
  /// \param type_b    Atom type of the second atom in the bond
  /// \param policy    Behavior if a bond parameter set matching the index or atom types is not
  ///                  found in the topology
  /// \{
  void setBondParameters(double new_keq, double new_leq, bool set_keq, bool set_leq, int parm_idx,
                         ExceptionResponse policy = ExceptionResponse::SILENT);

  void setBondParameters(double new_keq, double new_leq, bool set_keq, bool set_leq, char4 type_a,
                         char4 type_b, ExceptionResponse policy = ExceptionResponse::SILENT);
  /// \}

  /// \brief Alter the parameters for an angle parameter set.  This will alter the parameters in
  ///        the topology and all angle terms that use the parameters will then be use the new
  ///        parameters.  This function does not add or remove angle terms or parameter sets, nor
  ///        does it change the mapping of angle terms to parameter sets.
  ///
  /// Overloaded:
  ///   - Accept the atom types for each atom in the angle
  ///   - Accept the parameter set index for the angle
  ///
  /// \param new_keq   The new stiffness constant to apply
  /// \param new_teq   The new angular equilibrium constant to apply
  /// \param set_keq   Whether to change the stiffness constant
  /// \param set_teq   Whether to change the equilibrium
  /// \param parm_idx  Index of the angle parameter set from within the topology's table (one
  ///                  cannot depend on this index to be consistent across topologies of different
  ///                  systems making use of the same force field)
  /// \param type_a    Atom type of the first atom in the angle (this will be consistent across
  ///                  multiple topologies using the same force field)
  /// \param type_b    Atom type of the second atom in the angle
  /// \param type_c    Atom type of the third atom in the angle
  /// \param policy    Behavior if an angle parameter set matching the index or atom types is not
  ///                  found in the topology
  /// \{
  void setAngleParameters(double new_keq, double new_teq, bool set_keq, bool set_teq, int parm_idx,
                          ExceptionResponse policy = ExceptionResponse::SILENT);

  void setAngleParameters(double new_keq, double new_teq, bool set_keq, bool set_teq, char4 type_a,
                          char4 type_b, char4 type_c,
                          ExceptionResponse policy = ExceptionResponse::SILENT);
  /// \}

  /// \brief Alter the parameters for a cosine-based dihedral angle parameter set.  This will alter
  ///        the parameters in the topology and all dihedral terms that use the parameters will
  ///        then be use the new parameters.  This function does not add or remove dihedral terms
  ///        or parameter sets, nor does it change the mapping of dihedral terms to parameter sets.
  ///
  /// Overloaded:
  ///   - Accept the atom types for each atom in the dihedral angle
  ///   - Accept the parameter set index for the dihedral angle
  ///
  /// \param new_amplitude    The new cosine series term amplitude constant to apply
  /// \param new_phase_angle  The new phase angle to apply
  /// \param set_amplitude    Whether to change the amplitude constant
  /// \param set_phase_angle  Whether to change the equilibrium
  /// \param parm_idx         Index of the dihedral angle parameter set from within the topology's
  ///                         table (one cannot depend on this index to be consistent across
  ///                         topologies of different systems making use of the same force field)
  /// \param type_a           Atom type of the first atom in the dihedral (this will be consistent
  ///                         across multiple topologies using the same force field)
  /// \param type_b           Atom type of the second atom in the dihedral
  /// \param type_c           Atom type of the third atom in the dihedral
  /// \param type_d           Atom type of the fourth atom in the dihedral
  /// \param periodicity      Periodicity of the dihedral parameter set (this disambiguates
  ///                         dihedral parameters that apply on top of one another to the same
  ///                         four atom types)
  /// \param policy           Behavior if a dihedral parameter set matching the index or atom types
  ///                         and periodicity is not found in the topology
  /// \{
  void setDihedralParameters(double new_amplitude, double new_phase_angle, bool set_amplitude,
                             bool set_phase_angle, int parm_idx,
                             ExceptionResponse policy = ExceptionResponse::SILENT);

  void setDihedralParameters(double new_amplitude, double new_phase_angle, bool set_amplitude,
                             bool set_phase_angle, char4 type_a, char4 type_b, char4 type_c,
                             char4 type_d, double periodicity,
                             ExceptionResponse policy = ExceptionResponse::SILENT);
  /// \}

  /// \brief Alter the parameters for a Urey-Bradley parameter set.
  ///
  /// Overloaded:
  ///   - Accept the atom types for each atom in the Urey-Bradley term
  ///   - Accept the parameter set index for the Urey-Bradley term
  ///
  /// \param new_keq   The new stiffness constant to apply
  /// \param new_leq   The new equilibrium constant to apply
  /// \param set_keq   Whether to change the stiffness constant
  /// \param set_leq   Whether to change the equilibrium
  /// \param parm_idx  Index of the bond parameter set from within the topology's table (one cannot
  ///                  depend on this index to be consistent across topologies of different systems
  ///                  making use of the same force field)
  /// \param type_a    Atom type of the first atom in the Urey-Bradley term, the I atom (this will
  ///                  be consistent across multiple topologies using the same force field)
  /// \param type_b    Atom type of the second atom in the Urey-Bradley term (the K atom)
  /// \param policy    Behavior if a Urey-Bradley parameter set matching the index or atom types is
  ///                  not found in the topology  
  /// \{
  void setUreyBradleyParameters(double new_keq, double new_leq, bool set_keq, bool set_leq,
                                int parm_idx,
                                ExceptionResponse policy = ExceptionResponse::SILENT);

  void setUreyBradleyParameters(double new_keq, double new_leq, bool set_keq, bool set_leq,
                                char4 type_a, char4 type_b, char4 type_c,
                                ExceptionResponse policy = ExceptionResponse::SILENT);
  /// \}

  /// \brief Alter the characteristics of a CHARMM improper dihedral parameter set.
  ///
  /// Overloaded:
  ///   - Accept the atom types for each atom in the CHARMM improper term
  ///   - Accept the parameter set index for the CHARMM improper term
  ///
  /// \param new_stiffness    The new stiffness constant to apply
  /// \param new_phase_angle  The new phase angle constant to apply
  /// \param set_stiffness    Whether to change the stiffness constant
  /// \param set_phase_angle  Whether to change the equilibrium
  /// \param parm_idx         Index of the CHARMM improper dihedral parameter set from within the
  ///                         topology's table (one cannot depend on this index to be consistent
  ///                         across topologies of different systems making use of the same force
  ///                         field)
  /// \param type_a           Atom type of the first atom in the improper (this will be consistent
  ///                         across multiple topologies using the same force field)
  /// \param type_b           Atom type of the second atom in the improper
  /// \param type_c           Atom type of the third atom in the improper
  /// \param type_d           Atom type of the fourth atom in the improper
  /// \param policy           Behavior if a CHARMM improper dihedral parameter set matching the
  ///                         index or atom types and periodicity is not found in the topology 
  /// \{
  void setCharmmImprParameters(double new_stiffness, double new_phase_angle, bool set_stiffness,
                               bool set_phase_angle, int parm_idx,
                               ExceptionResponse policy = ExceptionResponse::SILENT);

  void setCharmmImprParameters(double new_stiffness, double new_phase_angle, bool set_stiffness,
                               bool set_phase_angle, char4 type_a, char4 type_b, char4 type_c,
                               char4 type_d, ExceptionResponse policy = ExceptionResponse::SILENT);
  /// \}

  
  /// \brief Set the implicit solvent model.  This will leverage a lot of hard-coded constants
  ///        to fill out some allocated but otherwise blank arrays and impart one particular
  ///        Generalized Born or other implicit solvent method.
  ///
  /// \param igb_in         The implicit solvent model to impart onto the system
  /// \param dielectric_in  The desired dielectric constant
  /// \param saltcon_in     The intended salt concentration (affects the GB decay parameter Kappa)
  /// \param radii_set      Radii to impart to the topology (this is often coupled to the choice of
  ///                       implicit solvent model, but for purposes of experimentation or new
  ///                       model development might be flexible)
  /// \param policy         Indicator of what to do if the topology's PB radii to not meet the
  ///                       implicit solvent model requirements, or there is some other problem
  void setImplicitSolventModel(ImplicitSolventModel igb_in, double dielectric_in = 80.0,
                               double saltcon_in = 0.0,
                               AtomicRadiusSet radii_set = AtomicRadiusSet::NONE,
                               ExceptionResponse policy = ExceptionResponse::WARN);

  /// \brief Set the name of the water residue
  ///
  /// Overloaded:
  ///   - Set the name based on a char4 (the native representation in the topology)
  ///   - Set the name based on up to four characters from a string
  ///
  /// \param new_name  The new water residue name to use
  /// \{
  void setWaterResidueName(const char4 new_name);
  void setWaterResidueName(const std::string &new_name);
  /// \}
  
private:
  
  // Title, version, and date stamps found in this topology
  char version_stamp[16];  ///< Version stamp for the program creating the topology
  tm date;                 ///< Date on which the topology was first created
  std::string title;       ///< Title of the topology, the first named field in an Amber prmtop
  std::string source;      ///< Name of the file from which this originated

  /// A list of identifiers for any number of force fields that the topology file references
  std::vector<Citation> force_fields;

  // Counts of atoms, residues, and other parts of the system
  int atom_count;                 ///< Total number of atoms and virtual sites
  int residue_count;              ///< Total number of residues, including solvent molecules
  int molecule_count;             ///< Total number of molecules in the system
  int largest_residue_size;       ///< Number of atoms in the largest residue
  int last_solute_residue;        ///< Last residue of the solute
  int last_solute_atom;           ///< Last atom of the solute (a typical solute is a biomolecule)
  int first_solvent_molecule;     ///< First molecule in what is deemed to be solvent
  int last_atom_before_cap;       ///< The solvent CAP will be found at the end of the atom list
  int implicit_copy_count;        ///< Number of path integral molecular dynamics slices, number
                                  ///<   of elastic band beads, or other population of a system
                                  ///<   with many implicit copies
  int largest_molecule_size;      ///< Size of the largest single molecule in the system
  int unconstrained_dof;          ///< Total number of degrees of freedom in the system, without
                                  ///<   any geometric constraints
  int constrained_dof;            ///< Total number of degrees of freedom in the system, with
                                  ///<   geometry constraints enforced
  Hybrid<int> descriptors;        ///< Topology stats founds in the first ten lines of an Amber
                                  ///<   prmtop
  Hybrid<int> residue_limits;     ///< Atomic indices marking boundaries of residues or molecules
  Hybrid<int> atom_struc_numbers; ///< The number of each atom in some input structure (defaults to
                                  ///<   a series of numbers beginning 1, 2, ... for the 0th, 1st,
                                  ///<   2nd atoms...
  Hybrid<int> residue_numbers;    ///< The residue number of each atom (defaults to a series of
                                  ///<   numbers begining 1, 2, ... for residues whose indices are
                                  ///<   0, 1, ... but can be changed according to additional
                                  ///<   input)
  Hybrid<int> molecule_limits;    ///< Prefix sum of the number of atoms per molecule, akin to
                                  ///<   residue_limits and indexing into molecule_contents

  // Atom and residue details
  Hybrid<int> atomic_numbers;             ///< Atomic Z-numbers
  Hybrid<int> mobile_atoms;               ///< Belly mask indicating mobile (1) or frozen (0) atoms
  Hybrid<int> molecule_membership;        ///< Distinct molecules, indexed from 0, to which each
                                          ///<   atom belongs
  Hybrid<int> molecule_contents;          ///< Indices of atoms making up each molecule, bounded by
                                          ///<   molecule_limits
  Hybrid<double> atomic_charges;          ///< Partial charges on all atoms and virtual sites
  Hybrid<double> atomic_masses;           ///< Masses of all atoms
  Hybrid<double> inverse_atomic_masses;   ///< Inverse masses for all atoms
  Hybrid<float> sp_atomic_charges;        ///< Partial charges on all atoms (single precision)
  Hybrid<float> sp_atomic_masses;         ///< Masses of all atoms (single precision)
  Hybrid<float> sp_inverse_atomic_masses; ///< Inverse masses for all atoms (single precision)
  Hybrid<char4> atom_names;               ///< Four letter names of all atoms, i.e. PDB names
  Hybrid<char4> atom_types;               ///< Four letter names of all atom types, i.e. CT
  Hybrid<char4> residue_names;            ///< Four letter names of all residues, i.e. PDB names

  // Relevant information for CHARMM valence terms (the CHARMM force field has many terms that
  // function in the same ways as the parameters that make up all of Amber, but since Amber works
  // exclusively in these terms (harmonic angles in the theta bending dimension, cosine-based
  // dihedral terms, cosine-based impropers using the same formula) they are referred to as
  // "AMBER" force field family parameters (see ForceFieldFamily enum class above).  The "CHARMM"
  // family of parameters include Urey-Bradley angle terms, CHARMM impropers (yet another take on
  // harmonic potentials), and CHARMM CMAP terms.  The Amber ff19SB force field uses CHARMM CMAPs,
  // but such terms always belong to the "CHARMM" family of parameters.
  int urey_bradley_term_count;                ///< Total number of Urey-Bradley angle stretch terms
  int charmm_impr_term_count;                 ///< Total number of CHARMM impropers
  int cmap_term_count;                        ///< Total number of CMAP terms
  int urey_bradley_parameter_count;           ///< Number of unique Urey-Bradley parameter pairs
  int charmm_impr_parameter_count;            ///< Number of unique CHARMM improper parameter pairs
  int cmap_surface_count;                     ///< Number of unique CMAP surfaces
  int urey_bradley_pert_term_count;           ///< Number of Urey-Bradley terms affecting the
                                              ///<   perturbed region
  int charmm_impr_pert_term_count;            ///< Number of CHARMM-style improper terms affecting
                                              ///<   the perturbed region
  int cmap_pert_term_count;                   ///< Number of CMAPs affecting the perturbed region
  int urey_bradleys_in_perturbed_group;       ///< Count of UB terms fully within perturbed group
  int charmm_imprs_in_perturbed_group;        ///< Count of CHARMM impropers fully within the
                                              ///<   perturbed group
  int cmaps_in_perturbed_group;               ///< Count of CMAPs fully within the perturbed group
  Hybrid<int> urey_bradley_i_atoms;           ///< First atom in each Urey-Bradley stretching term
  Hybrid<int> urey_bradley_k_atoms;           ///< Second atom in each Urey-Bradley stretching term
  Hybrid<int> urey_bradley_parameter_indices; ///< Index of Urey-Bradley parameters for each term
  Hybrid<int> urey_bradley_assigned_atoms;    ///< Indices of other atoms for Urey-Bradley terms
                                              ///<   assigned to a specific atom controlling each
                                              ///<   of them.  The atom groups controlled by atom
                                              ///<   a_idx are [urey_bradley_assigned_bounds[a_idx]
                                              ///<   ... urey_bradley_assigned_bounds[a_idx + 1])
  Hybrid<int> urey_bradley_assigned_index;    ///< Parameter indices from Urey-Bradley terms
                                              ///<   controlled by each atom.  The parameter
                                              ///<   indices controlled by atom a_idx are given by
                                              ///< [urey_bradley_assigned_bounds[a_idx] ...
                                              ///<  urey_bradley_assigned_bounds[a_idx + 1])
  Hybrid<int> urey_bradley_assigned_terms;    ///< Indices of Urey-Bradley terms which each atom
                                              ///<   commands
  Hybrid<int> urey_bradley_assigned_bounds;   ///< Bounds arrays for atom-controlled Urey-Bradley
                                              ///<   indexing
  Hybrid<int> charmm_impr_i_atoms;            ///< First atom in each CHARMM improper term
  Hybrid<int> charmm_impr_j_atoms;            ///< Second atom in each CHARMM improper term
  Hybrid<int> charmm_impr_k_atoms;            ///< Third atom in each CHARMM improper term
  Hybrid<int> charmm_impr_l_atoms;            ///< Fourth atom in each CHARMM improper term
  Hybrid<int> charmm_impr_parameter_indices;  ///< Index of each CHARMM improper's parameter pair
  Hybrid<int> charmm_impr_assigned_atoms;     ///< Indices of other atoms for CHARMM improper terms
                                              ///<   assigned to a specific atom controlling each
                                              ///<   of them.  The atom groups controlled by atom
                                              ///<   a_idx are given by
                                              ///<   [3 * charmm_impr_assigned_bounds[a_idx] ...
                                              ///<    3 * charmm_impr_assigned_bounds[a_idx + 1])
  Hybrid<int> charmm_impr_assigned_index;     ///< Parameter indices from Urey-Bradley terms
                                              ///<   controlled by each atom.  The parameter
                                              ///<   indices controlled by atom a_idx are given by
                                              ///< [charmm_impr_assigned_bounds[a_idx] ...
                                              ///<  charmm_impr_assigned_bounds[a_idx + 1])
  Hybrid<int> charmm_impr_assigned_terms;     ///< Indices of CHARMM improper terms which each atom
                                              ///<   commands
  Hybrid<int> charmm_impr_assigned_bounds;    ///< Bounds arrays for atom-controlled CHARMM
                                              ///<   improper indexing
  Hybrid<int> cmap_i_atoms;                   ///< 1st:UNK atom of each CMAP 2-D surface term
  Hybrid<int> cmap_j_atoms;                   ///< 2nd:1st atom of each CMAP 2-D surface term
  Hybrid<int> cmap_k_atoms;                   ///< 3rd:2nd atom of each CMAP 2-D surface term
  Hybrid<int> cmap_l_atoms;                   ///< 4th:3rd atom of each CMAP 2-D surface term
  Hybrid<int> cmap_m_atoms;                   ///< UNK:4th atom of each CMAP 2-D surface term
  Hybrid<int> cmap_surface_dimensions;        ///< Sizes of each CMAP surface (this array has a
                                              ///<   size based on the number of unique CMAP terms)
  Hybrid<int> cmap_surface_bounds;            ///< Starting points for each CMAP surface (this
                                              ///<   array has a size based on the number of unique
                                              ///<   CMAP terms)
  Hybrid<int> cmap_patch_bounds;              ///< Starting points for each aggregated CMAP value /
                                              ///<   derivative patch array.  The location of patch
                                              ///<   (j, k) for the ith CMAP (of grid dimension
                                              ///<   dim_i), given j and k within the range
                                              ///<   [0, dim_i), is obtained by taking the ith
                                              ///<   offset in cmap_patch_bounds plus
                                              ///<   16 * ((k * dim_i) + j).
  Hybrid<int> cmap_surface_indices;           ///< Surface parameter index of each CMAP term
  Hybrid<int> cmap_assigned_atoms;            ///< Indices of other atoms for CMAPs assigned to a
                                              ///<   specific atom controlling each CMAP.  The
                                              ///<   atom groups from each CMAP controlled by atom
                                              ///<   a_idx are [4 * cmap_assigned_bounds[a_idx]
                                              ///<   ... 4 * cmap_assigned_bounds[a_idx + 1])
  Hybrid<int> cmap_assigned_index;            ///< Parameter indices from CMAPs controlled by each
                                              ///<   atom.  The parameter indices for CMAPs
                                              ///<   controlled by atom a_idx are given by
                                              ///< [cmap_assigned_bounds[a_idx] ...
                                              ///<  cmap_assigned_bounds[a_idx + 1])
  Hybrid<int> cmap_assigned_terms;            ///< Indices of CMAP terms which each atom commands
  Hybrid<int> cmap_assigned_bounds;           ///< Bounds arrays for atom-controlled CMAP indexing
  Hybrid<double> urey_bradley_stiffnesses;    ///< Stiffness constant of each Urey-Bradley stretch
  Hybrid<double> urey_bradley_equilibria;     ///< Equilibrium length of each Urey-Bradley stretch
  Hybrid<double> charmm_impr_stiffnesses;     ///< CHARMM impropers are harmonic, too!
  Hybrid<double> charmm_impr_phase_angles;    ///< The "equilibria" for CHARMM impropers
  Hybrid<double> cmap_surfaces;               ///< Concatenated, column-major format matrices for
                                              ///<   every CMAP surface term
  Hybrid<double> cmap_phi_derivatives;        ///< Concatenated, column-major format matrices for
                                              ///<   each CMAP first dimension first derivative
  Hybrid<double> cmap_psi_derivatives;        ///< Concatenated, column-major format matrices for
                                              ///<   each CMAP second dimension first derivative
  Hybrid<double> cmap_phi_psi_derivatives;    ///< Concatenated, column-major format matrices for
                                              ///<   each CMAP cross derivative

  /// Concatenated, column-major format matrices for each CMAP aggregated patch representation.
  /// Patch (j, k) is a collection of coefficients for the surface (S), first dimension derivative
  /// (dS/dphi), second dimension first derivative (dS/dpsi), and cross derivative (d2S/dphi-dpsi).
  /// See cmap_patch_bounds for more on the array indexing.
  ///
  ///   [         S(j  ,k  ) S(j  ,k+1)      (dS/dpsi)(j  ,k  )      (dS/dpsi)(j  ,k+1)
  ///             S(j+1,k  ) S(j+1,k+1)      (dS/dpsi)(j+1,k  )      (dS/dpsi)(j+1,k+1)
  ///     (dS/dphi)(j  ,k  ) S(j  ,k+1) (dS/dphi-dpsi)(j  ,k  ) (dS/dphi-dpsi)(j  ,k+1)
  ///     (dS/dphi)(j+1,k  ) S(j+1,k+1) (dS/dphi-dpsi)(j+1,k  ) (dS/dphi-dpsi)(j+1,k+1) ]
  Hybrid<double> cmap_patches;
  Hybrid<float> sp_urey_bradley_stiffnesses;  ///< Stiffness constant of each Urey-Bradley stretch,
                                              ///<   in single precision
  Hybrid<float> sp_urey_bradley_equilibria;   ///< Equilibrium length of each Urey-Bradley stretch,
                                              ///<   in single precision
  Hybrid<float> sp_charmm_impr_stiffnesses;   ///< CHARMM impropers are harmonic, too!
  Hybrid<float> sp_charmm_impr_phase_angles;  ///< The "equilibria" for CHARMM impropers
  Hybrid<float> sp_cmap_surfaces;             ///< Concatenated, column-major format matrices for
                                              ///<   every CMAP surface term, in single precision
  Hybrid<float> sp_cmap_phi_derivatives;      ///< Single-precision form of eponymous array above
  Hybrid<float> sp_cmap_psi_derivatives;      ///< Single-precision form of eponymous array above
  Hybrid<float> sp_cmap_phi_psi_derivatives;  ///< Single-precision form of eponymous array above
  Hybrid<float> sp_cmap_patches;              ///< Single-precision CMAP patch representation

  // Relevant information for the standard ("Amber") bonded calculation
  int bond_term_with_hydrogen;           ///< Number of bonded interactions with hydrogen
  int angl_term_with_hydrogen;           ///< Number of angle interactions with hydrogen
  int dihe_term_with_hydrogen;           ///< Number of dihedral interactions with hydrogen
  int bond_term_without_hydrogen;        ///< Number of bonded interactions without hydrogen
  int angl_term_without_hydrogen;        ///< Number of angle interactions without hydrogen
  int dihe_term_without_hydrogen;        ///< Number of dihedral interactions without hydrogen
  int bond_term_count;                   ///< Total number of bonded interactions
  int angl_term_count;                   ///< Total number of bond angle interactions
  int dihe_term_count;                   ///< Total number of dihedral cosine terms in the system
  int bond_parameter_count;              ///< The number of unique bond parameter sets
  int angl_parameter_count;              ///< The number of unique angle parameter sets
  int dihe_parameter_count;              ///< Number of unique dihedral Fourier basis functions
  int bond_perturbation_term_count;      ///< Number of bonds affecting the perturbed region
  int angl_perturbation_term_count;      ///< Number of angles affecting the perturbed region
  int dihe_perturbation_term_count;      ///< Number of dihedrals affecting the perturbed region
  int bonds_in_perturbed_group;          ///< Count of bonds fully within the perturbed group
  int angls_in_perturbed_group;          ///< Count of angles fully within the perturbed group
  int dihes_in_perturbed_group;          ///< Count of dihedrals fully within the perturbed group
  int bonded_group_count;                ///< Number of bonded atom groups
  Hybrid<double> bond_stiffnesses;       ///< Stiffness of each bond stretch, kcal/mol-A^2
  Hybrid<double> bond_equilibria;        ///< Equilibrium lengths of all bonds, A
  Hybrid<double> angl_stiffnesses;       ///< Stiffness of each angle bend, kcal/mol-rad^2
  Hybrid<double> angl_equilibria;        ///< Equilibrium angle for all bending terms, radians
  Hybrid<double> dihe_amplitudes;        ///< Amplitudes of each dihedral cosine term, kcal/mol
  Hybrid<double> dihe_periodicities;     ///< Periodicity of each dihedral / torsion cosine term
  Hybrid<double> dihe_phase_angles;      ///< Phase angle of each dihedral / torsion cosine term
  Hybrid<float> sp_bond_stiffnesses;     ///< Stiffness of each bond stretch (single precision)
  Hybrid<float> sp_bond_equilibria;      ///< Equilibrium lengths of all bonds (single precision)
  Hybrid<float> sp_angl_stiffnesses;     ///< Angle bending stiffnesses (single precision)
  Hybrid<float> sp_angl_equilibria;      ///< Angle bending equilibria (single precision)
  Hybrid<float> sp_dihe_amplitudes;      ///< Amplitudes of torsion cosine terms (single precision)
  Hybrid<float> sp_dihe_periodicities;   ///< Periodicities of torsion terms (single precision)
  Hybrid<float> sp_dihe_phase_angles;    ///< Phase angles of torsion terms (single precision)
  Hybrid<int> bond_i_atoms;              ///< The first atom of each bond stretching term
  Hybrid<int> bond_j_atoms;              ///< The second atom of each bond stretching term
  Hybrid<int> bond_parameter_indices;    ///< Parameter index of each bond stretching term
  Hybrid<int> bond_assigned_atoms;       ///< Indices of other atoms for bond stretching terms
                                         ///<   assigned to a specific atom controlling each
                                         ///<   of them.  The atom groups controlled by atom
                                         ///<   a_idx are [bond_assigned_bounds[a_idx] ...
                                         ///<   bond_assigned_bounds[a_idx + 1])
  Hybrid<int> bond_assigned_index;       ///< Parameter indices from bond terms controlled by each
                                         ///<   atom.  The parameter indices controlled by atom
                                         ///<   a_idx are given by [bond_assigned_bounds[a_idx] ...
                                         ///<   bond_assigned_bounds[a_idx + 1])
  Hybrid<int> bond_assigned_terms;       ///< Indices of bond terms which this atom commands
  Hybrid<int> bond_assigned_bounds;      ///< Bounds arrays for atom-controlled bond indexing
  Hybrid<int> angl_i_atoms;              ///< The first atom of each angle bending term
  Hybrid<int> angl_j_atoms;              ///< The second atom of each angle bending term
  Hybrid<int> angl_k_atoms;              ///< The third atom of each angle bending term
  Hybrid<int> angl_parameter_indices;    ///< Parameter index of each angle bending term
  Hybrid<int> angl_assigned_atoms;       ///< Indices of other atoms for bond angle bending terms
                                         ///<   assigned to a specific atom controlling each of
                                         ///<   them.  The atom groups controlled by atom a_idx are
                                         ///<   given by [2 * angl_assigned_bounds[a_idx] ...
                                         ///<   2 * angl_assigned_bounds[a_idx + 1])
  Hybrid<int> angl_assigned_index;       ///< Parameter indices from bond angle bending terms
                                         ///<   controlled by each atom.  The parameter indices
                                         ///<   controlled by atom a_idx are given by
                                         ///<   [angl_assigned_bounds[a_idx] ...
                                         ///<   angl_assigned_bounds[a_idx + 1])
  Hybrid<int> angl_assigned_terms;       ///< Indices of angle terms which this atom commands
  Hybrid<int> angl_assigned_bounds;      ///< Bounds arrays for atom-controlled bond angle indexing
  Hybrid<int> dihe_i_atoms;              ///< First atom of each dihedral / torsion twisting term
  Hybrid<int> dihe_j_atoms;              ///< Second atom of each dihedral / torsion twisting term
  Hybrid<int> dihe_k_atoms;              ///< Third atom of each dihedral / torsion twisting term
  Hybrid<int> dihe_l_atoms;              ///< Fourth atom of each dihedral / torsion twisting term
  Hybrid<int> dihe_parameter_indices;    ///< Parameter index of each dihedral twisting term
  Hybrid<int> dihe14_parameter_indices;  ///< Parameter index for each dihedral's 1:4 screening
                                         ///<   factor pair
  Hybrid<int> dihe_assigned_atoms;       ///< Indices of other atoms for dihedral twisting terms
                                         ///<   assigned to a specific atom controlling each of
                                         ///<   them.  The atom groups controlled by atom a_idx are
                                         ///<   given by [3 * dihe_assigned_bounds[a_idx] ...
                                         ///<   3 * dihe_assigned_bounds[a_idx + 1])
  Hybrid<int> dihe_assigned_index;       ///< Parameter indices from dihedral twisting terms
                                         ///<   controlled by each atom.  The parameter indices
                                         ///<   controlled by atom a_idx are given by
                                         ///<   [dihe_assigned_bounds[a_idx] ...
                                         ///<   dihe_assigned_bounds[a_idx + 1])
  Hybrid<int> dihe_assigned_terms;       ///< Indices of dihedral terms which this atom commands
  Hybrid<int> dihe_assigned_bounds;      ///< Bounds arrays for atom-controlled dihedral indexing
  Hybrid<char4> bond_modifiers;          ///< Enumerations for special characteristics of bonds
  Hybrid<char4> angl_modifiers;          ///< Enumerations for special characteristics of angles
  Hybrid<char4> dihe_modifiers;          ///< Enumerations for special characteristics of dihedrals
  Hybrid<char4> bond_assigned_mods;      ///< Contents of bond_modifiers, re-arranged for
                                         ///<   atom-controlled bonds and indexed as
                                         ///<   bond_assigned_index above.
  Hybrid<char4> angl_assigned_mods;      ///< Contents of angl_modifiers, re-arranged for
                                         ///<   atom-controlled bond angles and indexed as
                                         ///<   angl_assigned_index above.
  Hybrid<char4> dihe_assigned_mods;      ///< Contents of dihe_modifiers, re-arranged for
                                         ///<   atom-controlled dihedrals and indexed as
                                         ///<   dihe_assigned_index above.

  // Information relevant to virtual site placement
  int virtual_site_count;                      ///< Number of v-sites / extra points
  int virtual_site_parameter_set_count;        ///< Number of parameter sets describing unique
                                               ///<   virtual site frames (the length of
                                               ///<   virtual_site_frame_dim1 and related arrays)
  Hybrid<int> virtual_site_atoms;              ///< List of atoms which are massless virtual sites
  Hybrid<int> virtual_site_frame_types;        ///< Frame types for each virtual site
  Hybrid<int> virtual_site_frame1_atoms;       ///< Parent atoms (frame I) for each virtual site
  Hybrid<int> virtual_site_frame2_atoms;       ///< Frame atom 2 for each virtual site
  Hybrid<int> virtual_site_frame3_atoms;       ///< Frame atom 3 (optional) for each virtual site
  Hybrid<int> virtual_site_frame4_atoms;       ///< Frame atom 4 (optional) for each virtual site
  Hybrid<int> virtual_site_parameter_indices;  ///< Parameter indices for virtual sites, indicating
                                               ///<   the frame type and the values of up to three
                                               ///<   dimensional measurements.  Virtual sites
                                               ///<   operate much like valence parameters: each
                                               ///<   has a unique set of two to four frame atoms
                                               ///<   and references a set of parameters.
  Hybrid<double> virtual_site_frame_dim1;      ///< First frame dimension for each v-site
  Hybrid<double> virtual_site_frame_dim2;      ///< Second (optional) frame dimension for each
                                               ///<   v-site
  Hybrid<double> virtual_site_frame_dim3;      ///< Third (optional) frame dimension for each
                                               ///<   v-site
  Hybrid<float> sp_virtual_site_frame_dim1;    ///< First frame dimension (single precision)
  Hybrid<float> sp_virtual_site_frame_dim2;    ///< Second frame dimension (single precision)
  Hybrid<float> sp_virtual_site_frame_dim3;    ///< Third frame dimension (single precision)

  // Relevant information for the non-bonded calculation
  int charge_type_count;                 ///< Number of distinct atomic partial charges
  int lj_type_count;                     ///< Number of distinct Lennard-Jones types
  int total_exclusions;                  ///< Total number of non-bonded exclusions, including 1:4
  int attenuated_14_type_count;          ///< The number of distinct 1:4 scaling factor pairs 
  int inferred_14_attenuations;          ///< Most 1:4 exclusions are actually attenuations
                                         ///<   according to parameters found in
                                         ///<   dihe_elec_screenings and dihe_vdw_screenings.
                                         ///<   This is the number of 1:4 exclusions whose
                                         ///<   attenuated strengths were not taken directly from
                                         ///<   dihedral interactions spanning the two atoms.  In
                                         ///<   these cases the 1:4 exclusion involves a virtual
                                         ///<   site which is not, itself, subject to dihedral
                                         ///<   interactions but inherits its exclusions from a
                                         ///<   parent atom which does paticipate in dihedral
                                         ///<   interactions.
  UnitCellType periodic_box_class;       ///< Type of periodic box: 0 = nonperiodic,
                                         ///<   1 = rectilinear prism, 2 = triclinic prism
  ImplicitSolventModel gb_style;         ///< The flavor of Generalized Born (or other implicit
                                         ///<   solvent) to use, i.e. Hawkins / Cramer / Truhlar.
                                         ///<   This information does not come as part of the
                                         ///<   original topology if loaded from an Amber prmtop,
                                         ///<   but can be added by a member function later.
  double dielectric_constant;            ///< Dielectric constant to take for implicit solvent
                                         ///<   calculations
  double salt_concentration;             ///< Salt concentration affecting implicit solvent
                                         ///<   calculations
  double coulomb_constant;               ///< Coulomb's constant in units of kcal-A/mol-e^2
                                         ///<   (Amber differs from other programs in terms of
                                         ///<   what this is, so it can be set here)
  double elec14_screening_factor;        ///< General screening factor for electrostatic 1:4
                                         ///<   interactions, applied if no pair-specific factor is
                                         ///<   available
  double vdw14_screening_factor;         ///< General screening factor for 1:4 van-der Waals
                                         ///<   interactions
  std::string pb_radii_set;              ///< The Poisson-Boltzmann radii set, also used in GB
  Hybrid<int> charge_indices;            ///< Atomic charge indices, 0 to charge_type_count - 1
  Hybrid<int> lennard_jones_indices;     ///< Lennard-Jones indices, 0 to atom_type_count - 1
  Hybrid<int> atom_exclusion_bounds;     ///< Bounds of non-bonded exclusions owned by each atom
  Hybrid<int> atom_exclusion_list;       ///< The list of all non-bonded exclusions, indexed by the
                                         ///<   prefix sum of atom_exclusion_counts.  Non-reflexive
                                         ///<   (if exclusion I:J is listed, then exclusion J:I
                                         ///<   will not be).  Exclusions for nb11, nb12, nb13, and
                                         ///<   nb14 arrays are reflexive.
  Hybrid<int> nb11_exclusion_bounds;     ///< Bounds of non-bonded 1:1 exclusions to each atom
  Hybrid<int> nb11_exclusion_list;       ///< Non-bonded 1:1 exclusions (virtual site::parent atom)
  Hybrid<int> nb12_exclusion_bounds;     ///< Bounds of non-bonded 1:2 exclusions to each atom
  Hybrid<int> nb12_exclusion_list;       ///< Non-bonded 1:2 exclusions (separated by one bond)
  Hybrid<int> nb13_exclusion_bounds;     ///< Bounds of non-bonded 1:3 exclusions to each atom
  Hybrid<int> nb13_exclusion_list;       ///< Non-bonded 1:3 exclusions (separated by two bonds)
  Hybrid<int> nb14_exclusion_bounds;     ///< Bounds of non-bonded 1:4 exclusions to each atom
  Hybrid<int> nb14_exclusion_list;       ///< Non-bonded 1:4 exclusions (separated by three bonds)
  Hybrid<int> infr14_i_atoms;            ///< First atoms in 1:4 attenuated pair interactions with
                                         ///<   inferred screening factors
  Hybrid<int> infr14_l_atoms;            ///< Second atoms in 1:4 attenuated pair interactions with
                                         ///<   inferred screening factors
  Hybrid<int> infr14_parameter_indices;  ///< Parameter indices of dihedral interactions guiding
                                         ///<   1:4 exclusions with inferred attenuations (this
                                         ///<   index goes into a pair of arrays of up to 32 unique
                                         ///<   pairs of 1:4 scaling factors, the same as indexed
                                         ///<   by dihe14_parameter_indices).
  Hybrid<int> neck_gb_indices;           ///< Indicies into separation and maximum value parameter
                                         ///<   tables for Mongan's "neck" GB implementations
  Hybrid<double> charge_parameters;      ///< Atomic partial charges, condensed to a list of unique
                                         ///<   values 0 ... charge_type_count - 1 and indexed by
                                         ///<   the charge_indices array
  Hybrid<double> lj_a_values;            ///< Lennard-Jones A coefficients, U = A/r^12 - B/r^6
  Hybrid<double> lj_b_values;            ///< Lennard-Jones B coefficients, U = A/r^12 - B/r^6
  Hybrid<double> lj_c_values;            ///< Lennard-Jones C coefficients, U = A/r^12 - B/r^6
  Hybrid<double> lj_14_a_values;         ///< Lennard-Jones A coefficients for 1:4 interactions
  Hybrid<double> lj_14_b_values;         ///< Lennard-Jones B coefficients for 1:4 interactions
  Hybrid<double> lj_14_c_values;         ///< Lennard-Jones C coefficients for 1:4 interactions
  Hybrid<double> lj_sigma_values;        ///< Lennard-Jones sigma parameters, back-computed from
                                         ///<   combinations of A and B coefficients for each pair
  Hybrid<double> lj_14_sigma_values;     ///< Lennard-Jones 1:4 sigma parameters, back-computed
                                         ///<   from combinations of 1:4 A and B coefficients
  Hybrid<double> lj_type_corrections;    ///< Long-ranged homogeneity coefficients for vdW energy
  Hybrid<double> attn14_elec_factors;    ///< Array of attenuated electrostatic strengths for 1:4
                                         ///<   non-bonded interactions
  Hybrid<double> attn14_vdw_factors;     ///< Array of attenuated van-der Waals (Lennard-Jones)
                                         ///<   strengths for 1:4 non-bonded interactions
  Hybrid<double> atomic_pb_radii;        ///< Radii of all atoms according to the pb_radii_set
  Hybrid<double> gb_screening_factors;   ///< Generalized Born screening factors for all atoms
  Hybrid<double> gb_alpha_parameters;    ///< Generalized Born alpha parameters for each atom (set
                                         ///<   according to a particular GB scheme, not part of
                                         ///<   an Amber topology file but added later)
  Hybrid<double> gb_beta_parameters;     ///< Generalized Born beta parameters for each atom
  Hybrid<double> gb_gamma_parameters;    ///< Generalized Born gamma parameters for each atom
  Hybrid<float> sp_charge_parameters;    ///< Condensed, unique atomic partial charges in single
                                         ///<   precision representation
  Hybrid<float> sp_lj_a_values;          ///< Lennard-Jones A coefficients (single precision)
  Hybrid<float> sp_lj_b_values;          ///< Lennard-Jones B coefficients (single precision)
  Hybrid<float> sp_lj_c_values;          ///< Lennard-Jones C coefficients (single precision)
  Hybrid<float> sp_lj_14_a_values;       ///< Lennard-Jones 1:4 A coefficients (single precision)
  Hybrid<float> sp_lj_14_b_values;       ///< Lennard-Jones 1:4 B coefficients (single precision)
  Hybrid<float> sp_lj_14_c_values;       ///< Lennard-Jones 1:4 C coefficients (single precision)
  Hybrid<float> sp_lj_sigma_values;      ///< Lennard-Jones sigma parameters (single precision)
  Hybrid<float> sp_lj_14_sigma_values;   ///< Lennard-Jones 1:4 sigma parameters (single precision)
  Hybrid<float> sp_lj_type_corrections;  ///< Long-ranged homogeneity factors (single precision)
  Hybrid<float> sp_attn14_elec_factors;  ///< Array of attenuated electrostatic strengths for 1:4
                                         ///<   non-bonded interactions, single-precision
  Hybrid<float> sp_attn14_vdw_factors;   ///< Array of attenuated van-der Waals strengths for 1:4
                                         ///<   non-bonded interactions, single-precision
  Hybrid<float> sp_atomic_pb_radii;      ///< P.B. Radii of all atoms (single precision)
  Hybrid<float> sp_gb_screening_factors; ///< Generalized Born screening factors (single precision)
  Hybrid<float> sp_gb_alpha_parameters;  ///< Single-precision Generalized Born alpha parameters
  Hybrid<float> sp_gb_beta_parameters;   ///< Single-precision Generalized Born beta parameters
  Hybrid<float> sp_gb_gamma_parameters;  ///< Single-precision Generalized Born gamma parameters

  // MD propagation algorithm directives
  ShakeSetting use_bond_constraints;           ///< Toggles use of bond length constraints
  SettleSetting use_settle;                    ///< Toggles analytic constraints on rigid water
  PerturbationSetting use_perturbation_info;   ///< Toggles perturbations
  SolventCapSetting use_solvent_cap_option;    ///< Toggles the solvent cap option
  PolarizationSetting use_polarization;        ///< Toggles use of polarization
  char4 water_residue_name;                    ///< Name of the water residue, can be compared to
                                               ///<   residue_names
  std::string bond_constraint_mask;            ///< Atoms involved in bond length constraints
  std::string bond_constraint_omit_mask;       ///< Atoms to be excluded from bond length
                                               ///<   constraints
  int bond_constraint_count;                   ///< Bonds with lengths constrained by SHAKE or
                                               ///<   RATTLE.  This is more for informative
                                               ///<   purposes, as groups of constrained bonds
                                               ///<   which connect at the same central atom are
                                               ///<   handled as their own individual objects.
  int nonrigid_particle_count;                 ///< A rigid water is one non-rigid particle.  A
                                               ///<   protein with N atoms and no bond length
                                               ///<   constraints is N particles.
  int settle_group_count;                      ///< Number of rigid water molecules, or other
                                               ///<   molecules with the right features, subject to
                                               ///<   SETTLE
  int settle_parameter_count;                  ///< The number of unique SETTLE parameter sets in
                                               ///<   the topology (there is probably only one, but
                                               ///<   there could be more if some water is made of
                                               ///<   isotopes or there is some other triangular,
                                               ///<   rigid molecule).
  int constraint_group_count;                  ///< Number of groups of constrained bond groups
  int constraint_parameter_count;              ///< Number of unique constraint group parameter
                                               ///<   sets (each constraint group will reference
                                               ///<   one set of the appropriate size, but there
                                               ///<   could be multiple groups that use the exact
                                               ///<   same sequence of parameters)
  Hybrid<int> settle_oxygen_atoms;             ///< List of oxygen atoms involved in fast waters
  Hybrid<int> settle_hydro1_atoms;             ///< First hydrogen atoms involved in fast waters
  Hybrid<int> settle_hydro2_atoms;             ///< Second hydrogen atoms involved in fast waters
  Hybrid<int> settle_parameter_indices;        ///< Indices into a list of (perhaps more than one)
                                               ///<   set of SETTLE parameters for each
                                               ///<   SETTLE-suitable molecule in the system
  Hybrid<int> constraint_group_atoms;          ///< List of all atoms involved in constraint
                                               ///<   groups.  In each group, the central atom, to
                                               ///<   which all others bind, is listed first.  It
                                               ///<   is the first atom of any of the constrained
                                               ///<   bond.  All other atoms are the distal
                                               ///<   termini of constrained bonds.  Serial A-B-C-D
                                               ///<   constraints and constrained ring systems are
                                               ///<   not supported.
  Hybrid<int> constraint_group_bounds;         ///< Bounds array for the constrained group atoms,
                                               ///<   length constraint_group_count + 1 to provide
                                               ///<   the limits on the final group.
  Hybrid<int> constraint_parameter_indices;    ///< Parameter indices for each constraint group.
                                               ///<   This is tricky: the index indicates an
                                               ///<   element of the bounds array to read, and the
                                               ///<   size of some constraint group k is then
                                               ///<   indicated by
                                               ///<   constraint_parameter_bounds[k + 1] -
                                               ///<   constraint_parameter_bounds[k].  There are
                                               ///<   then a mass and a target length for every
                                               ///<   atom in the constraint group, given in
                                               ///<   constraint_inverse_masses and
                                               ///<   constraint_squared_lengths, respectively.
  Hybrid<int> constraint_parameter_bounds;     ///< Bounds array for constraint_inverse_masses and
                                               ///<   constraint_squared_lengths below
  Hybrid<double> settle_mormt;                 ///< Proportional mass of "oxygen" in SETTLE systems
  Hybrid<double> settle_mhrmt;                 ///< Proportional mass of "hydrogen" in SETTLE
                                               ///<   systems
  Hybrid<double> settle_ra;                    ///< Internal distance measurement of SETTLE groups
  Hybrid<double> settle_rb;                    ///< Internal distance measurement of SETTLE groups
  Hybrid<double> settle_rc;                    ///< Internal distance measurement of SETTLE groups
  Hybrid<double> settle_invra;                 ///< Internal distance measurement of SETTLE groups
  Hybrid<double> constraint_inverse_masses;    ///< Inverse masses for atoms of constrained groups
  Hybrid<double> constraint_squared_lengths;   ///< Target lengths for constrained bonds
  Hybrid<float> sp_settle_mormt;               ///< Proportional mass of "oxygen" in SETTLE systems
  Hybrid<float> sp_settle_mhrmt;               ///< Proportional mass of "hydrogen" in SETTLE
                                               ///<   systems
  Hybrid<float> sp_settle_ra;                  ///< Internal distance measurement of SETTLE groups
  Hybrid<float> sp_settle_rb;                  ///< Internal distance measurement of SETTLE groups
  Hybrid<float> sp_settle_rc;                  ///< Internal distance measurement of SETTLE groups
  Hybrid<float> sp_settle_invra;               ///< Internal distance measurement of SETTLE groups
  Hybrid<float> sp_constraint_inverse_masses;  ///< Inverse masses for atoms of constrained groups
  Hybrid<float> sp_constraint_squared_lengths; ///< Target lengths for constrained bonds

  // Atom, atom type, and residue name overflow keys
  Hybrid<char4> atom_overflow_names;    ///< Codified names of atoms which were too long
  Hybrid<char4> atom_overflow_types;    ///< Codified names of atom types which were too long
  Hybrid<char4> residue_overflow_names; ///< Codified names of residues which were too long

  // Information currently unused
  int unused_nhparm;                ///< Number of unique hydrogen bonding parameters
  int unused_nparm;                 ///< If this is an Amber prmtop, indicates whether addles
                                    ///<   created the file
  int unused_natyp;                 ///< Atom types pertaining to hydrogen bonding interactions
  int hbond_10_12_parameter_count;  ///< Number of 10-12 hybrodgen bonding parameters
  int heavy_bonds_plus_constraints; ///< Bonds without hydrogen, plus constrained bonds
  int heavy_angls_plus_constraints; ///< Angles without hydrogen, plus constrained angles
  int heavy_dihes_plus_constraints; ///< Dihedrals without hydrogen, plus constrained dihedrals
  Hybrid<int> tree_joining_info;    ///< Tree joining info, more relevant to Amber tleap than here
  Hybrid<int> last_rotator_info;    ///< Last rotatable atom depending on any given atom in a chain
  Hybrid<double> solty_info;        ///< Long ago reserved for future use in Amber
  Hybrid<double> hbond_a_values;    ///< Hydrogen bond A coefficients, similar to L-J A coeffients
  Hybrid<double> hbond_b_values;    ///< Hydrogen bond B coefficients, similar to L-J B coeffients
  Hybrid<double> hbond_cutoffs;     ///< Cutoff values for hydrogen bonding of any given type
  Hybrid<char4> tree_symbols;       ///< Tree symbols, i.e. BLA, more relevant to tleap than here

  /// Hybrid data structures store data arrays that will be accessible to accelerators.  These are
  /// contiguous arrays for each type of data: integers and doubles.  There are no unsigned
  /// integers, from the standpoint of this topology struct, but there are single-precision (fp32)
  /// representations of the topology parameters.  The char4 data, encoding atom parameters, is
  /// intended to serve for debugging purposes in kernels.
  /// \{
  Hybrid<int> int_data;
  Hybrid<double> double_data;
  Hybrid<float> float_data;
  Hybrid<char4> char4_data;
  /// \}
  
  /// \brief Function to wrap pointer re-assignments after copy constructor initialization or
  ///        a copy assignment operation.
  void rebasePointers();

  /// Allow the AtomGraphStage to access private members directly to modify and guide construction
  /// of topologies used in STORMM's production calculations.
  friend class AtomGraphStage;
};

} // namespace topology
} // namespace stormm

#include "atomgraph.tpp"

#endif
