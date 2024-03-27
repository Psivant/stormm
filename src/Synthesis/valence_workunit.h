// -*-c++-*-
#ifndef STORMM_VALENCE_WORKUNIT_H
#define STORMM_VALENCE_WORKUNIT_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/gpu_details.h"
#include "Constants/behavior.h"
#include "Potential/energy_enumerators.h"
#include "Restraints/restraint_apparatus.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "synthesis_enumerators.h"

namespace stormm {
namespace synthesis {

using card::GpuDetails;
using card::Hybrid;
using constants::PrecisionModel;
using energy::EvaluateForce;
using energy::EvaluateEnergy;
using energy::ValenceKernelSize;
using restraints::RestraintApparatus;
using restraints::RestraintKit;
using topology::AtomGraph;
using topology::ConstraintKit;
using topology::ValenceKit;
using topology::VirtualSiteKit;

/// \brief The maximum value for the atom limit in a valence work unit--any higher and the amount
///        of __shared__ memory in a block of 256 threads might need to be stretched.
constexpr int maximum_valence_work_unit_atoms = 544;
constexpr int half_valence_work_unit_atoms = maximum_valence_work_unit_atoms / 2;
constexpr int quarter_valence_work_unit_atoms = maximum_valence_work_unit_atoms / 4;
constexpr int eighth_valence_work_unit_atoms = maximum_valence_work_unit_atoms / 8;

/// \brief The minimium value for the atom limit in a valence work unit--there should be no need
///        to force the size below 72.  The size of woork units in pmemd.cuda was 128.
constexpr int minimum_valence_work_unit_atoms = eighth_valence_work_unit_atoms;

/// \brief The maximum tolerated number of recursive calls to various atom searching functions.
///        The halo region needed to determine the correct movement of any given atom cannot chain
///        indefinitely.  A value of 2 should permit all the recursive calls that might be needed,
///        but allow more to let whatever erroneous pattern become clear for debugging purposes.
constexpr int max_atom_search_stacks = 8;

/// \brief The maximum size of any constraint group.  This limit is imposed to ensure that one
///        warp can handle each constraint group, and the Intel GPU warp size is 16.
constexpr int max_constraint_group_size = 16;
  
/// \brief Object to track how different valence terms in a topology are delegated.  Valence work
///        units may evaluate a valence term without being responsible for moving both atoms, or
///        even for moving any of the atoms at all.  Each valence term is only fully delegated
///        once valence work units that evaluate it are responsible for moving all atoms that the
///        valence term contains.  In order for a work unit to move any aotm, it must evaluate all
///        valence terms that include that atom.
class ValenceDelegator {
public:

  /// \brief The object is constructed based on a single topology and oversees the construction of
  ///        an array of valence work units.
  ///
  /// \param ag_in  Pointer to the topology containing valence terms to delegate among work units
  /// \param ra_in  Pointer to the complete collection of restraints applicable to the system
  ValenceDelegator(const AtomGraph *ag_in, const RestraintApparatus *ra_in = nullptr);

  /// \brief Get the number of work units to which a particular atom is (currently) assigned.
  ///
  /// \param atom_index  The atom of interest, as indexed in the original topology
  int getAtomAssignmentCount(int atom_index) const;

  /// \brief Get the index of the first unassigned atom in the topology.
  int getFirstUnassignedAtom() const;

  /// \brief Get a list of bonds affecting a given list of atoms.  The result is returned (pruned
  ///        for unique items) as a series of term indices referencing the original topology.
  ///
  /// \param atom_indices  Topological indices of the atoms of interest
  std::vector<int> getBondAffectors(const std::vector<int> &atom_indices) const;
  
  /// \brief Get a list of harmonic angle terms affecting a given list of atoms.  The result is
  ///        returned (pruned for unique items) as a series of term indices referencing the
  ///        original topology.
  ///
  /// \param atom_indices  Topological indices of the atoms of interest
  std::vector<int> getAngleAffectors(const std::vector<int> &atom_indices) const;
 
  /// \brief Get a list of cosine-based dihedrals affecting a given list of atoms.  The result is
  ///        returned (pruned for unique items) as a series of term indices referencing the
  ///        original topology.
  ///
  /// \param atom_indices  Topological indices of the atoms of interest
  std::vector<int> getDihedralAffectors(const std::vector<int> &atom_indices) const;
  
  /// \brief Get a list of Urey-Bradley harmonic angles affecting a given list of atoms.  The
  ///        result is returned (pruned for unique items) as a series of term indices referencing
  ///        the original topology.
  ///
  /// \param atom_indices  Topological indices of the atoms of interest
  std::vector<int> getUreyBradleyAffectors(const std::vector<int> &atom_indices) const;
   
  /// \brief Get a list of CHARMM improper dihedral terms affecting a given list of atoms.  The
  ///        result is returned (pruned for unique items) as a series of term indices referencing
  ///        the original topology.
  ///
  /// \param atom_indices  Topological indices of the atoms of interest
  std::vector<int> getCharmmImproperAffectors(const std::vector<int> &atom_indices) const;

  /// \brief Get a list of CMAP terms affecting a given list of atoms.  The result is returned
  ///        (pruned for unique items) as a series of term indices referencing the original
  ///        topology.
  ///
  /// \param atom_indices  Topological indices of the atoms of interest
  std::vector<int> getCmapAffectors(const std::vector<int> &atom_indices) const;

  /// \brief Get a list of inferred 1:4 attenuated non-bonded interactions affecting a given list
  ///        of atoms.  The result is returned (pruned for unique items) as a series of term
  ///        indices referencing the original topology.
  ///
  /// \param atom_indices  Topological indices of the atoms of interest
  std::vector<int> getInferred14Affectors(const std::vector<int> &atom_indices) const;

  /// \brief Get a list of positional restraints affecting a given list of atoms.  The result is
  ///        returned (pruned for unique items) as a series of term indices referencing the
  ///        original topology.
  ///
  /// \param atom_indices  Topological indices of the atoms of interest
  std::vector<int> getPositionalRestraintAffectors(const std::vector<int> &atom_indices) const;
  
  /// \brief Get a list of distance restraints affecting a given list of atoms.  The result is
  ///        returned (pruned for unique items) as a series of term indices referencing the
  ///        original topology.
  ///
  /// \param atom_indices  Topological indices of the atoms of interest
  std::vector<int> getDistanceRestraintAffectors(const std::vector<int> &atom_indices) const;
  
  /// \brief Get a list of three-point angle restraints affecting a given list of atoms.  The
  ///        result is returned (pruned for unique items) as a series of term indices referencing
  ///        the original topology.
  ///
  /// \param atom_indices  Topological indices of the atoms of interest
  std::vector<int> getAngleRestraintAffectors(const std::vector<int> &atom_indices) const;

  /// \brief Get a list of four-point dihedral angle restraints affecting a given list of atoms.
  ///        The result is returned (pruned for unique items) as a series of term indices
  ///        referencing the original topology.
  ///
  /// \param atom_indices  Topological indices of the atoms of interest
  std::vector<int> getDihedralRestraintAffectors(const std::vector<int> &atom_indices) const;

  /// \brief Get a list of virtual site frames affecting a given list of atoms.  The result is
  ///        returned (pruned for unique items) as a series of term indices referencing the
  ///        original topology.
  ///
  /// \param atom_indices  Topological indices of the atoms of interest
  std::vector<int> getVirtualSiteAffectors(const std::vector<int> &atom_indices) const;

  /// \brief Get a list of SETTLE (fast, rigid water) constraint groups affecting a given list of
  ///        atoms.  The result is returned (pruned for unique items) as a series of term indices
  ///        referencing the original topology.
  ///
  /// \param atom_indices  Topological indices of the atoms of interest
  std::vector<int> getSettleGroupAffectors(const std::vector<int> &atom_indices) const;

  /// \brief Get a list of hub-and-spoke constraint groups affecting a given list of atoms.  The
  ///        result is returned (pruned for unique items) as a series of term indices referencing
  ///        the original topology.
  ///
  /// \param atom_indices  Topological indices of the atoms of interest
  std::vector<int> getConstraintGroupAffectors(const std::vector<int> &atom_indices) const;

  /// \brief Find all partners of a given atom such that the work unit will be able to correctly
  ///        move the atom.  This implies all other participants in any constraint group that
  ///        includes the atom, or (if the atom is a virtual site) all frame atoms and any
  ///        constraint groups which they participate in.  The returned list provides indices of
  ///        atoms as they are found in the original topology.  It is pruned to remove duplicates.
  ///
  /// \param atom_idx      Topological index of the atom of interest
  /// \param caller_stack  Cumulative stack indicaing previous recursive calls.  If this gets too
  ///                      long, it will trigger a runtime error.
  std::vector<int> findMovementPartners(int atom_idx,
                                        const std::vector<int> &caller_stack = {}) const;

  /// \brief Find all partners of a given atom such that the work unit will be able to correctly
  ///        evaluate the force need to move the atom.  This implies other atoms that participate
  ///        in any valence terms, inferred 1:4 interactions, or restraints involving the atom, as
  ///        well as other atoms which participate in a virtual site frame with which the atom of
  ///        interest is involved.  If the atom of interest is one of the frame atoms to a virtual
  ///        site, the relevant atoms also include any atoms involved in valence interactions
  ///        (i.e. inferred non-bonded 1:4 interactions) with the virtual site itself.
  ///
  /// \param atom_idx      Topological index of the atom of interest
  /// \param caller_stack  Cumulative stack indicaing previous recursive calls.  If this gets too
  ///                      long, it will trigger a runtime error.
  std::vector<int> findForcePartners(int atom_idx,
                                     const std::vector<int> &caller_stack = {}) const;

  /// \brief Accumulate a list of all atoms which have bearing on the way that a particular atom
  ///        shall move.  This routine will loop over the ValenceDelgator's bounded arrays for
  ///        valence terms, restraints, virtual sites, and various constraint groups to determine
  ///        any other atoms that must move, and any atoms which can contribute forces.
  ///
  /// \param atom_index  Topological index of the atom in question
  std::vector<int> getUpdateDependencies(const int atom_index) const;
  
  /// \brief Get the work unit currently assigned to update an atom's position and velocity.
  ///        The value of -1 signifies that no work unit has yet been assigned to handle the atom.
  ///
  /// \param atom_index  Topological index of the atom of interest
  int getUpdateWorkUnit(int atom_index) const;

  /// \brief Get the work unit currently assigned to contribute a particular bond term's energy
  ///        into the global accumulator.
  ///
  /// \param bond_index  Topological index of the bond term of interest
  int getBondAccumulatorWorkUnit(const int bond_index) const;
  
  /// \brief Get the work unit currently assigned to contribute a particular harmonic angle term's
  ///        energy into the global accumulator.
  ///
  /// \param angl_index  Topological index of the angle term of interest
  int getAngleAccumulatorWorkUnit(const int angl_index) const;

  /// \brief Get the work unit currently assigned to contribute a particular cosine-based dihedral
  ///        term's energy into the global accumulator.
  ///
  /// \param dihe_index  Topological index of the cosine-based dihedral term of interest
  int getDihedralAccumulatorWorkUnit(const int dihe_index) const;

  /// \brief Get the work unit currently assigned to contribute a particular Urey-Bradley term's
  ///        energy into the global accumulator.
  ///
  /// \param ubrd_index  Topological index of the Urey-Bradley term of interest
  int getUreyBradleyAccumulatorWorkUnit(const int ubrd_index) const;

  /// \brief Get the work unit currently assigned to contribute a particular CHARMM improper
  ///        dihedral term's energy into the global accumulator.
  ///
  /// \param cimp_index  Topological index of the CHARMM improper term of interest
  int getCharmmImproperAccumulatorWorkUnit(const int cimp_index) const;

  /// \brief Get the work unit currently assigned to contribute a particular CMAP term's energy
  ///        into the global accumulator.
  ///
  /// \param cmap_index  Topological index of the CMAP surface term of interest
  int getCmapAccumulatorWorkUnit(const int cmap_index) const;

  /// \brief Get the work unit currently assigned to contribute a particular inferred 1:4
  ///        attenuated pair interaction energy into the global accumulator.
  ///
  /// \param infr14_index  Topological index of the inferred 1:4 interaction of interest
  int getInferred14AccumulatorWorkUnit(const int infr14_index) const;

  /// \brief Get the work unit currently assigned to contribute a particular positional restraint
  ///        penalty energy into the global accumulator.
  ///
  /// \param rposn_index  Index of the positional restraint in the restraint apparatus
  int getPositionalRestraintAccumulatorWorkUnit(const int rposn_index) const;

  /// \brief Get the work unit currently assigned to contribute a particular distance restraint
  ///        penalty energy into the global accumulator.
  ///
  /// \param rbond_index  Index of the distance restraint in the restraint apparatus
  int getDistanceRestraintAccumulatorWorkUnit(const int rbond_index) const;

  /// \brief Get the work unit currently assigned to contribute a particular three-point angle
  ///        restraint penalty energy into the global accumulator.
  ///
  /// \param rangl_index  Index of the three-point angle restraint in the restraint apparatus
  int getAngleRestraintAccumulatorWorkUnit(const int rangl_index) const;

  /// \brief Get the work unit currently assigned to contribute a particular four-point dihedral
  ///        restraint penalty energy into the global accumulator.
  ///
  /// \param rdihe_index  Index of the four-point dihedral restraint in the restraint apparatus
  int getDihedralRestraintAccumulatorWorkUnit(const int rdihe_index) const;

  /// \brief Get a pointer to the topology that this delegator was built for.
  const AtomGraph* getTopologyPointer() const;
  
  /// \brief Get a pointer to the topology that this delegator was built for.
  const RestraintApparatus* getRestraintApparatusPointer() const;
    
  /// \brief Check that an atom is present in a particular work unit.
  ///
  /// \param atom_index  The topological index of the atom of interest
  /// \param vwu_index   Index of the ValenceWorkUnit in a list that this delegator is managing
  bool checkPresence(int atom_index, int vwu_index) const;

  /// \brief Mark the addition of an atom to a specific ValenceWorkUnit.
  ///
  /// \param vwu_index   Index of the ValenceWorkUnit receiving the new atom
  /// \param atom_index  Index of the atom to add, referencing the original topology
  void markAtomAddition(int vwu_index, int atom_index);

  /// \brief Mark the updates (position, velocity) of a particular atom in the topology as the
  ///        responsibility of a particular work unit in the list.  Returns TRUE if the assignment
  ///        is successful, or FALSE if the assignment has already been made to some work unit.
  ///
  /// \param atom_index  The topological index of the atom of interest
  /// \param vwu_index   Index of the ValenceWorkUnit in a list that this delegator is managing  
  bool setUpdateWorkUnit(int atom_index, int vwu_index);

  /// \brief Mark the work unit that will accumulate the potential due to a specific bond
  ///        interaction.  Returns TRUE if the assignment is successful, FALSE if the assignment
  ///        has already been made to some other work unit.
  ///
  /// \param bond_index  The topological index of the bond term of interest
  /// \param vwu_index   Index of the ValenceWorkUnit in a list that this delegator is managing
  bool setBondAccumulatorWorkUnit(int bond_index, int vwu_index);
  
  /// \brief Mark the work unit that will accumulate the potential due to a specific harmonic angle
  ///        interaction.  Returns TRUE if the assignment is successful, FALSE if the assignment
  ///        has already been made to some other work unit.
  ///
  /// \param angl_index  The topological index of the harmonic angle term of interest
  /// \param vwu_index   Index of the ValenceWorkUnit in a list that this delegator is managing
  bool setAngleAccumulatorWorkUnit(int angl_index, int vwu_index);

  /// \brief Mark the work unit that will accumulate the potential due to a specific cosine-based
  ///        dihedral interaction.  Returns TRUE if the assignment is successful, FALSE if the
  ///        assignment has already been made to some other work unit.
  ///
  /// \param dihe_index  The topological index of the dihedral term of interest
  /// \param vwu_index   Index of the ValenceWorkUnit in a list that this delegator is managing
  bool setDihedralAccumulatorWorkUnit(int dihe_index, int vwu_index);

  /// \brief Mark the work unit that will accumulate the potential due to a specific Urey-Bradley
  ///        interaction.  Returns TRUE if the assignment is successful, FALSE if the assignment
  ///        has already been made to some other work unit.
  ///
  /// \param ubrd_index  The topological index of the Urey-Bradley term of interest
  /// \param vwu_index   Index of the ValenceWorkUnit in a list that this delegator is managing
  bool setUreyBradleyAccumulatorWorkUnit(int ubrd_index, int vwu_index);

  /// \brief Mark the work unit that will accumulate the potential due to a specific CHARMM
  ///        improper dihedral interaction.  Returns TRUE if the assignment is successful, FALSE
  ///        if the assignment has already been made to some other work unit.
  ///
  /// \param cimp_index  The topological index of the CHARMM improper dihedral term of interest
  /// \param vwu_index   Index of the ValenceWorkUnit in a list that this delegator is managing
  bool setCharmmImproperAccumulatorWorkUnit(int cimp_index, int vwu_index);

  /// \brief Mark the work unit that will accumulate the potential due to a specific CMAP splined
  ///        surface interaction.  Returns TRUE if the assignment is successful, FALSE if the
  ///        assignment has already been made to some other work unit.
  ///
  /// \param cmap_index  The topological index of the CMAP term of interest
  /// \param vwu_index   Index of the ValenceWorkUnit in a list that this delegator is managing
  bool setCmapAccumulatorWorkUnit(int cmap_index, int vwu_index);

  /// \brief Mark the work unit that will accumulate the potential due to a specific inferred 1:4
  ///        non-bonded, attenuated interaction.  Returns TRUE if the assignment is successful,
  ///        FALSE if the assignment has already been made to some other work unit.
  ///
  /// \param infr14_index  The topological index of the CMAP term of interest
  /// \param vwu_index     Index of the ValenceWorkUnit in a list that this delegator is managing
  bool setInferred14AccumulatorWorkUnit(int infr14_index, int vwu_index);

  /// \brief Mark the work unit that will accumulate the potential due to a positional restraint.
  ///        Returns TRUE if the assignment is successful, FALSE if the assignment has already
  ///        been made to some other work unit.
  ///
  /// \param rposn_index  The restraint apparatus index of the restraint of interest
  /// \param vwu_index    Index of the ValenceWorkUnit in a list that this delegator is managing
  bool setPositionalRestraintAccumulatorWorkUnit(int rposn_index, int vwu_index);

  /// \brief Mark the work unit that will accumulate the potential due to a distance restraint.
  ///        Returns TRUE if the assignment is successful, FALSE if the assignment has already
  ///        been made to some other work unit.
  ///
  /// \param rbond_index  The restraint apparatus index of the restraint of interest
  /// \param vwu_index    Index of the ValenceWorkUnit in a list that this delegator is managing
  bool setDistanceRestraintAccumulatorWorkUnit(int rbond_index, int vwu_index);

  /// \brief Mark the work unit that will accumulate the potential due to a three-point angle
  ///        restraint.  Returns TRUE if the assignment is successful, FALSE if the assignment has
  ///        already been made to some other work unit.
  ///
  /// \param rangl_index  The restraint apparatus index of the restraint of interest
  /// \param vwu_index    Index of the ValenceWorkUnit in a list that this delegator is managing
  bool setAngleRestraintAccumulatorWorkUnit(int rangl_index, int vwu_index);

  /// \brief Mark the work unit that will accumulate the potential due to a four-point dihedral
  ///        restraint.  Returns TRUE if the assignment is successful, FALSE if the assignment has
  ///        already been made to some other work unit.
  ///
  /// \param rdihe_index  The restraint apparatus index of the restraint of interest
  /// \param vwu_index    Index of the ValenceWorkUnit in a list that this delegator is managing
  bool setDihedralRestraintAccumulatorWorkUnit(int rdihe_index, int vwu_index);

private:
  int atom_count;                 ///< The number of atoms in the system overall (taken from the
                                  ///<   topology)
  int first_unassigned_atom;      ///< The first atom index with no current work unit assignments
  int max_presence_allocation;    ///< The allocation size for the maximum number of work units
                                  ///<   to which any atom can be assigned.

  // The following lists are quick to construct and contain a double-, triple-, or even a
  // quintuple-counting of all information in the topology's own "atom assignment" arrays for the
  // same valence terms.  Is there a better way?
  std::vector<int> bond_affector_list;    ///< List of all harmonic bonds affecting a given atom
  std::vector<int> bond_affector_bounds;  ///< Bounds array for bond_affector_list
  std::vector<int> angl_affector_list;    ///< List of all bond angles affecting a given atom
  std::vector<int> angl_affector_bounds;  ///< Bounds array for angl_affector_list
  std::vector<int> dihe_affector_list;    ///< List of all cosine dihedrals affecting a given atom
  std::vector<int> dihe_affector_bounds;  ///< Bounds array for dihe_affector_list
  std::vector<int> ubrd_affector_list;    ///< List of Urey-Bradley terms affecting a given atom
  std::vector<int> ubrd_affector_bounds;  ///< Bounds array for ubrd_affector_list
  std::vector<int> cimp_affector_list;    ///< List of CHARMM impropers affecting a given atom
  std::vector<int> cimp_affector_bounds;  ///< Bounds array for cimp_affector_list
  std::vector<int> cmap_affector_list;    ///< List of all CMAP terms affecting a given atom
  std::vector<int> cmap_affector_bounds;  ///< Bounds array for cmap_affector_list
  std::vector<int> infr_affector_list;    ///< List of all inferred attenuated 1:4 interactions
                                          ///<   affecting a given atom (dihedral-linked attenuated
                                          ///<   1:4 interactions, while relevant to the energy
                                          ///<  surface, do not imply new atoms to consider outside
                                          ///<   of the dihedrals that own them), and new, relevant
                                          ///<   atoms are what these lists are all about).
  std::vector<int> infr_affector_bounds;  ///< Bounds array for cmap_affector_list

  // The virtual site affector arrays will be pre-allcoated to hold virtual sites with up to four
  // frame atoms apiece.
  std::vector<int> vste_affector_list;    ///< List of all virtual sites that name a given atom as
                                          ///<   one of their frame atoms
  std::vector<int> vste_affector_bounds;  ///< Bounds array for vste_affector_list

  // SHAKE and RATTLE groups likewise must have all atoms present in a work group in order to
  // evaluate bond constraints.
  std::vector<int> cnst_affector_list;    ///< List of all constraint groups affecting any atom
                                          ///<   (these must be detected through a special search
                                          ///<   so that only a work unit with all affected atoms
                                          ///<   present will be tasked with evaluating the
                                          ///<   constraints).
  std::vector<int> cnst_affector_bounds;  ///< Bounds array for cnst_affector_list

  // SETTLE groups are special constraint groups for a fast rigid water implementation, kept
  // distinct because of the different treatment they get in the dynamics loop.
  std::vector<int> sett_affector_list;    ///< List of all SETTLE groups affecting any atom
                                          ///<   (only one work group will be tasked to evaluate
                                          ///<   each SETTLE fast water)
  std::vector<int> sett_affector_bounds;  ///< Bounds array for sett_affector_list

  // Positional restraints apply to a single atom, but multiple restraints could still, in theory,
  // restrain a single atom to a ring or other confined region of space.  A bounds array is still
  // necessary.  Other restraints affect multiple atoms and map like their valence term analogs.
  std::vector<int> rposn_affector_list;    ///< List of positional restraints affecting any atom
  std::vector<int> rposn_affector_bounds;  ///< Bounds array for rposn_affector_list
  std::vector<int> rbond_affector_list;    ///< List of bond restraints affecting any atom
  std::vector<int> rbond_affector_bounds;  ///< Bounds array for rbond_affector_list
  std::vector<int> rangl_affector_list;    ///< List of three-point restraints affecting any atom
  std::vector<int> rangl_affector_bounds;  ///< Bounds array for rangl_affector_list
  std::vector<int> rdihe_affector_list;    ///< List of four-point restraints affecting any atom
  std::vector<int> rdihe_affector_bounds;  ///< Bounds array for rdihe_affector_list

  // Individual atoms must leave a record of their whereabouts in the valence work units for rapid
  // retrieval of their locations
  std::vector<int> work_unit_assignment_count;  ///< The numbers of work units in which each atom
                                                ///<   can be found
  std::vector<int> work_unit_presence;          ///< Lists of the work units in which each atom is
                                                ///<   found.  This is a column-format matrix with
                                                ///<   atom_count columns and a number of rows
                                                ///<   expanded as needed to accommodate the
                                                ///<   largest entry in work_unit_assignments.
  std::vector<int> assigned_update_work_units;  ///< List indices for the work units assigned to
                                                ///<   update (and log) each atom's position and
                                                ///<   velocity in the master coordinate set

  // Simple energy terms must leave a record of the work units that will be responsible for
  // accumulating their values in the total reported energy.  Multiple work units may evaluate
  // each energy term, but only one must contribute.
  std::vector<int> assigned_bond_acc_work_units;   ///< Work units assigned to log each bond term
  std::vector<int> assigned_angl_acc_work_units;   ///< Work units assigned to log each angle term
  std::vector<int> assigned_dihe_acc_work_units;   ///< Work units assigned to log each dihedral
  std::vector<int> assigned_ubrd_acc_work_units;   ///< Work units assigned to log each
                                                   ///<   Urey-Bradley term
  std::vector<int> assigned_cimp_acc_work_units;   ///< Work units assigned to log each CHARMM
                                                   ///<   improper dihedral term
  std::vector<int> assigned_cmap_acc_work_units;   ///< Work units assigned to log each CMAP term
  std::vector<int> assigned_infr14_acc_work_units; ///< Work units assigned to log each inferred
                                                   ///<   1:4 attenuated pair interaction
  std::vector<int> assigned_rposn_acc_work_units;  ///< Work units assigned to log each positional
                                                   ///<   restraint penalty energy
  std::vector<int> assigned_rbond_acc_work_units;  ///< Work units assigned to log each distance
                                                   ///<   restraint penalty energy
  std::vector<int> assigned_rangl_acc_work_units;  ///< Work units assigned to log each three-point
                                                   ///<   angle restraint penalty energy
  std::vector<int> assigned_rdihe_acc_work_units;  ///< Work units assigned to log each four-point
                                                   ///<   dihedral angle restraint penalty energy
  
  
  /// Pointers to the original topology and restraint apparatus that formed this object
  const AtomGraph *ag_pointer;
  const RestraintApparatus *ra_pointer;

  /// \brief Make whatever needed space for the arrays indicating which work units any particular
  ///        atom is present in.  This is not a compact array with a bounds list, but rather a
  ///        padded array with space for each atom, due to the frequency with which it might be
  ///        updated.
  ///
  /// \param n_units  The maximum number of work units that any atom might be a part of.  Allocate
  ///                 an array for this many inclusions.
  void resizePresenceArrays(int n_units);
  
  /// \brief Allocate the necessary space for this work unit
  ///
  /// \param vk   Valence term abstract from the original topology
  /// \param vsk  Virtual site abstract from the original topology
  /// \param cnk  Constraint group abstract from the original topology
  /// \param rar  Restraint apparatus abstract
  void allocate();

  /// \brief Fill the arrays describing how different atoms are affected by each potential term
  ///        (including restraints), each virtual site, and each constraint.
  ///
  /// \param vk   Valence term abstract from the original topology
  /// \param vsk  Virtual site abstract from the original topology
  /// \param cnk  Constraint group abstract from the original topology
  /// \param rar  Restraint apparatus abstract
  void fillAffectorArrays(const ValenceKit<double> &vk, const VirtualSiteKit<double> &vsk,
                          const ConstraintKit<double> &cnk,
                          const RestraintKit<double, double2, double4> &rar);
};
  
/// \brief An object to collect the components of a valence work unit (which will also track frozen
///        atoms to implement coordinate updates, velocity updates, and constraints).  While the
///        work unit is encoded in the AtomGraphSynthesis object, the assembly is best done by a
///        dedicated object with plenty of its own methods operating on a single topology
///        (AtomGraph).  All systems in the AtomGraphSynthesis are designed to function
///        independently of one another--the only difference is that they have consensus tables of
///        most parameters and differen atom indexing.  Translating a valence work unit into a
///        list of instructions within an AtomGraphSynthesis is therefore a critical member
///        function of this class.
class ValenceWorkUnit {
public:

  /// \brief The constructor takes a specific input topology (multiple systems using the same
  ///        topology in an AtomGraphSynthesis can thereby take the same valence work unit and
  ///        translate the atom indices as appropriate rather than regenerating the work unit
  ///        for many equivalent systems).  Mapping starts from a specific atom and proceeds until
  ///        a maximum number of atoms has been accumulated in order to process as many related
  ///        valence terms as possible.
  ///
  /// \param vdel_in        Valence delegator managing the creation of this valence work unit
  /// \param tvwu_coverage  Array spanning all atoms in the system to mark whether any one of them
  ///                       has been included in the valence work unit currently under
  ///                       construction.  This is distinct from arrays with similar functionality
  ///                       held by the ValenceDelegator object, which track whether an atom has
  ///                       been included as an update priority in any valence work unit.  Having
  ///                       this work-unit specific resource ensures that atoms are not included
  ///                       multiple times in the import array.  It is cleared after each work
  ///                       unit's construction so as to not require re-allocation.
  /// \param list_index_in  Index of this unit in a larger list (the unit should remember its own
  ///                       index number, for the purposes of coordinating with other work units)
  /// \param seed_atom_in   The first atom to incorporate into the work unit.  Subsequent atoms
  ///                       will be either bonded in some chain to the seed, retracing previous
  ///                       topological indices whereby previous work units left some atoms
  ///                       behind, or jumping forward to the next new molecule.
  /// \param max_atoms_in   The maximum number of atoms to accumulate in the work unit
  ValenceWorkUnit(ValenceDelegator *vdel_in, std::vector<int> *tvwu_coverage, int list_index_in,
                  int seed_atom_in, int max_atoms_in = maximum_valence_work_unit_atoms);

  /// \brief Get the number of atoms currently imported into this work unit.
  int getImportedAtomCount() const;

  /// \brief Get the number of atoms currently set to be moved by this work unit.
  int getMovedAtomCount() const;

  /// \brief Get the number of atoms currently set to be updated by this work unit.
  int getUpdatedAtomCount() const;

  /// \brief Get the list index of this work unit.
  int getListIndex() const;

  /// \brief Get the minimum topological atom index of any used by this work unit.
  int getMinAtomIndex() const;

  /// \brief Get the maximum topological atom index of any used by this work unit.
  int getMaxAtomIndex() const;
  
  /// \brief Get the maximum atom count that this work unit can hold.
  int getMaxAtoms() const;

  /// \brief Get the list of imported atoms.
  ///
  /// \param atom_offset  Offset of the atoms to add to the topological indices, if importing from
  ///                     a synthesis of many systems
  std::vector<int> getAtomImportList(int atom_offset = 0) const;

  /// \brief Get a specific imported atom.
  ///
  /// \param slot         Index of the atom of interest from within the work unit's array (this
  ///                     will return a topological index found in position atom_idx of the local
  ///                     import list)
  /// \param atom_offset  Offset of the atoms to add to the topological indices, if importing from
  ///                     a larger synthesis of topologies / systems
  int getImportedAtomIndex(int slot, int atom_offset = 0) const;
  
  /// \brief Get bitmasks of moving atoms and atoms that this work unit is assigned to update.
  ///        Bits signify 1 for an atom being moving or an update assignment.  The masks are
  ///        encoded in a tuple with the x members making a mask for each segment of moving atoms
  ///        and the y members making a mask for each segment of assigned update atoms.
  std::vector<uint2> getAtomManipulationMasks() const;

  /// \brief Get the padded size of the largest constraint group.  The padding extends the size of
  ///        each constraint group in the work unit to the minimum factor of two, with a maximum
  ///        allowed size of 16.
  int getPaddedConstraintInstructionCount() const;

  /// \brief Get a vector describing the number of each type of item this work unit can be tasked
  ///        to perform.
  std::vector<int> getTaskCounts() const;

  /// \brief Compute and store a vector of the bond instructions.  This function accepts parameter
  ///        interpretation tables in order to produce instructions for a collated topology
  ///        handling many systems (AtomGraphSynthesis).
  ///
  /// \param bond_param_map  Map of the bond parameter sets (optional)
  /// \param ubrd_param_map  Map of the Urey-Bradley parameter sets (optional)
  void storeCompositeBondInstructions(const std::vector<int> &bond_param_map = {},
                                      const std::vector<int> &ubrd_param_map = {});

  /// \brief Store a vector of the harmonic angle instructions.  This function accepts a parameter
  ///        interpretation table in order to produce instructions for a collated topology
  ///        handling many systems (AtomGraphSynthesis).
  ///
  /// \param parameter_map  Map of the angle parameter sets (optional)
  void storeAngleInstructions(const std::vector<int> &parameter_map = {});

  /// \brief Store a vector of the composite (cosine-based) dihedral, associated 1:4 interactions,
  ///        as well as CHARMM improper dihedral instructions.  This function accepts parameter
  ///        interpretation tables in order to produce instructions for a collated topology
  ///        handling many systems (AtomGraphSynthesis).
  ///
  /// \param dihe_param_map  Map of the dihedral parameter sets (optional)
  /// \param dihe_param_map  Map of the 1:4 scaling factor parameter pairs (optional)
  /// \param dihe_param_map  Map of the CHARMM harmonic improper parameter sets (optional)
  void storeCompositeDihedralInstructions(const std::vector<int> &dihe_param_map = {},
                                          const std::vector<int> &dihe14_param_map = {},
                                          const std::vector<int> &cimp_param_map = {});
  
  /// \brief Store a vector of the CMAP instructions.  This function accepts a parameter
  ///        interpretation table in order to produce instructions for a collated topology
  ///        handling many systems (AtomGraphSynthesis).
  ///
  /// \param parameter_map  Map of the one topology's CMAP surface indices onto the synthesis
  ///                       (optional)
  void storeCmapInstructions(const std::vector<int> &parameter_map = {});
  
  /// \brief Store a vector of the inferred 1:4 attenuated pair interaction instructions.  This
  ///        function accepts a parameter interpretation table in order to produce instructions
  ///        for a collated topology handling many systems (AtomGraphSynthesis).
  ///
  /// \param parameter_map  Map of the one topology's attenuated interaction scaling factors onto
  ///                       the synthesis (optional)
  void storeInferred14Instructions(const std::vector<int> &parameter_map = {});

  /// \brief Store a vector of the positional restraint instructions for this work unit.  This
  ///        function accepts parameter interpretation tables showing how the raw list of
  ///        restraint settings (k(2,3), r(1,2,3,4), and x / y / z targets) can be condensed by
  ///        an AtomGraphSynthesis.
  ///
  /// \param kr_param_map   Mapping for k(2,3) and r(1,2,3,4) settings (optional)
  /// \param xyz_param_map  Mapping for Cartesian coordinate targets (optional)
  void storePositionalRestraintInstructions(const std::vector<int> &kr_param_map = {},
                                            const std::vector<int> &xyz_param_map = {});

  /// \brief Store a vector of the distance restraint instructions for this work unit.  This
  ///        function accepts a parameter interpretation table showing how the raw list of
  ///        restraint settings (k(2,3) and r(1,2,3,4)) can be condensed by an AtomGraphSynthesis.
  ///
  /// \param kr_param_map   Mapping for k(2,3) and r(1,2,3,4) settings (optional)
  void storeDistanceRestraintInstructions(const std::vector<int> &kr_param_map = {});

  /// \brief Store a vector of the three-point angle restraint instructions for this work unit.
  ///        This function accepts a parameter interpretation table showing how the raw list of
  ///        restraint settings (k(2,3) and r(1,2,3,4)) can be condensed by an AtomGraphSynthesis.
  ///
  /// \param kr_param_map   Mapping for k(2,3) and r(1,2,3,4) settings (optional)
  void storeAngleRestraintInstructions(const std::vector<int> &kr_param_map = {});

  /// \brief Store a vector of the four-point dihedral restraint instructions for this work unit.
  ///        This function accepts a parameter interpretation table showing how the raw list of
  ///        restraint settings (k(2,3) and r(1,2,3,4)) can be condensed by an AtomGraphSynthesis.
  ///
  /// \param kr_param_map   Mapping for k(2,3) and r(1,2,3,4) settings (optional)
  void storeDihedralRestraintInstructions(const std::vector<int> &kr_param_map = {});

  /// \brief Store a vector of the virtual site instructions for this work unit.  This function
  ///        accepts a parameter interpretation table showing how the raw list of unique virtual
  ///        site frames from one topology maps into a larger selection kept by an
  ///        AtomGraphSynthesis.
  ///
  /// \param parameter_map  Mapping for virtual site frame specifications (optional)
  void storeVirtualSiteInstructions(const std::vector<int> &parameter_map = {});

  /// \brief Store a vector of the SETTLE constraint group instructions for this work unit.  This
  ///        function accepts a parameter interpretation table showing how the raw list of unique
  ///        SETTLE geometries in one topology maps into a possibly more diverse list curated by
  ///        an AtomGraphSynthesis.
  ///
  /// \param parameter_map  Mapping for SETTLE group mass and geometry specifications (optional)
  void storeSettleGroupInstructions(const std::vector<int> &parameter_map = {});

  /// \brief Store a vector of the hub-and-spoke constraint group instructions for this work unit.
  ///        This function accepts a parameter interpretation table showing how the raw list of
  ///        unique constraint groups in one topology maps into a possibly more diverse list
  ///        curated by an AtomGraphSynthesis.
  ///
  /// \param parameter_map       Mapping for SETTLE group mass and geometry specifications
  ///                            (optional)
  /// \param group_param_bounds  Bounds for constraint group parameter sets in an
  ///                            AtomGraphSynthesis (required if parameter_map is supplied, to
  ///                            properly align the parameter indices of each instruction
  void storeConstraintGroupInstructions(const std::vector<int> &parameter_map = {},
                                        const std::vector<int> &group_param_bounds = {});

  /// \brief Get the stored vector of composite bond instructions.
  const std::vector<uint2>& getCompositeBondInstructions() const;

  /// \brief Get the stored vector of angle instructions.
  const std::vector<uint2>& getAngleInstructions() const;

  /// \brief Get the stored vector of composite dihedral instructions.
  const std::vector<uint3>& getCompositeDihedralInstructions() const;

  /// \brief Get the stored vector of CMAP instructions.
  const std::vector<uint2>& getCmapInstructions() const;

  /// \brief Get the stored vector of CMAP instructions.
  const std::vector<uint>& getInferred14Instructions() const;

  /// \brief Get the stored vector of positional restraint instructions.
  const std::vector<uint2>& getPositionalRestraintInstructions() const;

  /// \brief Get the stored vector of distance restraint instructions.
  const std::vector<uint2>& getDistanceRestraintInstructions() const;

  /// \brief Get the stored vector of three-point angle restraint instructions.
  const std::vector<uint2>& getAngleRestraintInstructions() const;

  /// \brief Get the stored vector of four-point dihedral restraint instructions.
  const std::vector<uint2>& getDihedralRestraintInstructions() const;

  /// \brief Get the stored vector of virtual site placement instructions.
  const std::vector<uint2>& getVirtualSiteInstructions() const;  
  
  /// \brief Get the stored vector of SETTLE constraint group instructions.
  const std::vector<uint2>& getSettleGroupInstructions() const;  
  
  /// \brief Get the stored vector of hub-and-spoke constraint group instructions.
  const std::vector<uint2>& getConstraintGroupInstructions() const;  

  /// \brief Get a specific composite bond instruction.
  ///
  /// \param index  Index of the instruction to retrieve
  uint2 getCompositeBondInstruction(int index) const;

  /// \brief Get a specific angle instruction.
  ///
  /// \param index  Index of the instruction to retrieve
  uint2 getAngleInstruction(int index) const;

  /// \brief Get a specific composite dihedral instruction.
  ///
  /// \param index  Index of the instruction to retrieve
  uint3 getCompositeDihedralInstruction(int index) const;

  /// \brief Get a specific CMAP instruction.
  ///
  /// \param index  Index of the instruction to retrieve
  uint2 getCmapInstruction(int index) const;

  /// \brief Get a specific CMAP instruction.
  ///
  /// \param index  Index of the instruction to retrieve
  uint getInferred14Instruction(int index) const;

  /// \brief Get a specific positional restraint instruction.
  ///
  /// \param index  Index of the instruction to retrieve
  uint2 getPositionalRestraintInstruction(int index) const;

  /// \brief Get a specific distance restraint instruction.
  ///
  /// \param index  Index of the instruction to retrieve
  uint2 getDistanceRestraintInstruction(int index) const;

  /// \brief Get a specific three-point angle restraint instruction.
  ///
  /// \param index  Index of the instruction to retrieve
  uint2 getAngleRestraintInstruction(int index) const;

  /// \brief Get a specific four-point dihedral restraint instruction.
  ///
  /// \param index  Index of the instruction to retrieve
  uint2 getDihedralRestraintInstruction(int index) const;

  /// \brief Get a specific virtual site placement instruction.
  ///
  /// \param index  Index of the instruction to retrieve
  uint2 getVirtualSiteInstruction(int index) const;  
  
  /// \brief Get a specific SETTLE constraint group instruction.
  ///
  /// \param index  Index of the instruction to retrieve
  uint2 getSettleGroupInstruction(int index) const;  
  
  /// \brief Get a specific hub-and-spoke constraint group instruction.
  ///
  /// \param index  Index of the instruction to retrieve
  uint2 getConstraintGroupInstruction(int index) const;  

  /// \brief Get the bitstrings indicating which energetic interactions each work unit is
  ///        responsible for accumulating into the official energy outputs.
  ///
  /// \param vtask  The type of task accumulator
  const std::vector<uint>& getAccumulationFlags(VwuTask vtask) const;

  /// \brief Get the bitstrings indicating which of the imported (cached) atoms this work unit is
  ///        responsible for updating in the global postion and velocity arrays.
  const std::vector<uint>& getAtomUpdateFlags() const;
  
  /// \brief Get the topological indices of each task assigned to this work unit.  Assignment of
  ///        an energy / force-producing term or constraint group does not imply that a work unit
  ///        is responsible for updating the final positions of all atoms involved.  Composite
  ///        task lists must be accessed with the special-purpose functions below.
  ///
  /// \param vtask  The type of task, i.e. bonded interactions, or SETTLE constraint groups
  const std::vector<int>& getSimpleTaskList(VwuTask vtask) const;

  /// \brief Get the composite bond tasks assigned to this work unit.  This will return a vector
  ///        of concatenated bond and Urey-Bradley term indices into the original topology.
  ///        Interpreting which is which requires the corresponding vector of composite bond term
  ///        instructions.
  const std::vector<int>& getCompositeBondTaskList() const;

  /// \brief Get the composite dihedral tasks assigned to this work unit.  This will return a
  ///        vector of tuples containing the topological indices of dihedrals or CHARMM impropers
  ///        that the work unit evaluates.  Interpretation of the tuples will depend on knowing
  ///        whether each term index pertains to a standard cosine-based dihedral or a CHARMM
  ///        improper dihedral, for which the corresponding instructions list must be accessed.
  const std::vector<int2>& getCompositeDihedralTaskList() const;
  
  /// \brief Get the pointer to the ValenceDelegator managing the creation of this object.
  ValenceDelegator* getDelegatorPointer();

  /// \brief Get a pointer to the topology for which this work unit applies.
  const AtomGraph* getTopologyPointer() const;

  /// \brief Get a pointer to the restraint collection for which this work unit applies.
  const RestraintApparatus* getRestraintApparatusPointer() const;
  
  /// \brief Set the list index of this work unit, in the event that the list of work units for
  ///        a particular topology needs to be re-ordered.
  ///
  /// \param list_index_in  The new list index for the work unit
  void setListIndex(int list_index_in);

  /// \brief Set the atom limit for a valence work unit.  This can be useful in situations where
  ///        it is desirable to form several work units out of a single molecule, despite there
  ///        being enough room in just one to hold all atoms of the molecule.
  ///
  /// \param new_limit  The new limit on the number of atoms.  This cannot be lower than the
  ///                   number of atoms already in the work unit.
  void setAtomLimit(int new_limit);
  
  /// \brief Add a new atom to a work unit.  This will update the associated ValenceDelegator and
  ///        all assignments therein.
  ///
  /// \param atom_index  Index of the atom of interest
  void addNewAtomImport(int atom_index);

  /// \brief Add a new atom to the list of updates that a work unit shall perform.  The atom must
  ///        already be part of the atom import list.
  void addNewAtomUpdate(const int atom_index);

  /// \brief Create the move list for atoms in the work unit.  Any atom that the work unit is
  ///        responsible for updating must be moved by the work unit, but also any atom sharing a
  ///        constraint group with one of the atoms which is on the official update list.
  void makeAtomMoveList();

  /// \brief Sort the atom lists (import, movement, and update) of this work unit into ascending
  ///        order.  This will optimize memory access when reading the atoms and set the stage
  ///        for mapping valence terms / atom groups to the local list.
  void sortAtomSets();

  /// \brief Create a bit mask spanning the atom imports, marking all of those that the work unit
  ///        is responsible for updating in the global position and velocity arrays.
  void makeAtomUpdateMask();
  
  /// \brief Log all activities of this work unit: valence terms, restraints, virtual sites, and
  ///        constraints.  This will translate the topological indices of atoms into indices of
  ///        the local import list.
  void logActivities();
  
private:
  int imported_atom_count;  ///< Number of atoms imported by the work unit
  int moved_atom_count;     ///< Number of atoms moved by the work unit
  int updated_atom_count;   ///< Number of atoms updated by the work unit
  int bond_term_count;      ///< Number of bond terms in the work unit
  int angl_term_count;      ///< Number of angle terms in the work unit
  int dihe_term_count;      ///< Number of cosine-based dihedral terms in the work unit
  int ubrd_term_count;      ///< Number of Urey-Bradley terms in the work unit
  int cbnd_term_count;      ///< Combined number of bond and Urey-Bradley terms (these share a
                            ///<   form, and thus a code pathway)
  int cimp_term_count;      ///< Number of CHARMM harmonic improper dihedral terms in the work unit
  int cdhe_term_count;      ///< Number of composite dihedral terms, sweeping up dihedrals that
                            ///<   affect the same four atoms, implicit 1:4 interactions, and
                            ///<   CHARMM impropers
  int cmap_term_count;      ///< Number of CMAP terms in the work unit
  int infr14_term_count;    ///< Number of inferred 1:4 interactions (not linked to any dihedral)
  int rposn_term_count;     ///< Number of positional restraints handled by this work unit
  int rbond_term_count;     ///< Number of distance restraints handled by this work unit
  int rangl_term_count;     ///< Number of angle restraints handled by this work unit
  int rdihe_term_count;     ///< Number of dihedral restraints handled by this work unit
  int cnst_group_count;     ///< Number of SHAKE or RATTLE groups managed by this work
                            ///<   unit (excludes SETTLE-constrained waters)
  int sett_group_count;     ///< Number of SETTLE-constrained rigid waters managed by the work unit
  int vste_count;           ///< Number of virtual sites managed by this work unit
  int list_index;           ///< Index of the work unit in a larger list of similar
                            ///<   objects coordinating to cover an entire topology
  int min_atom_index;       ///< Lowest topological index of any imported atom
  int max_atom_index;       ///< Highest topological index of any imported atom
  int atom_limit;           ///< Largest number of atoms that this work unit can hold

  /// The list of imported atoms, indicating indices into the original topology.  The position of
  /// each atom in this list indicates its local index within the work unit, as referenced by
  /// subsequent ????_(i,j,k,...)_atoms arrays.
  std::vector<int> atom_import_list;

  /// The list of atoms that this work unit shall compute complete forces on and thereby move.
  /// All of these atoms (and perhaps more) must, obviously, be contained in atom_import_list.
  std::vector<int> atom_move_list;

  /// The list of atoms that this work unit will be responsible for updating after movement.
  /// All of these atoms must be contained in atom_move_list.  If a work unit moves an atom, it
  /// would be acceptable for that unit to log the updated position and velocity, but only one
  /// work unit will be responsible for doing that.
  std::vector<int> atom_update_list;
  
  /// A mask of local atoms that the work unit is responsible for updating in the global position
  /// arrays.  All local atoms, except virtual sites, will be moved according to the forces acting
  /// upon them (it is faster to move the 1-2% of atoms that do not absolutely need to move than
  /// to try and figure out a detail like that).  Once the particles are moved, however, and
  /// perhaps processed with bond length constraints or virtual site placement, the question of
  /// whether to update the particle positions in the global array is critical.  Traversing the
  /// list of imported atoms while accessing the corresponding bits of this mask answers that
  /// question.
  std::vector<uint> atom_update_mask;
  
  // Valence terms for the work unit, typical force field elements
  std::vector<int> bond_term_list;    ///< List of harmonic bonds for which this work unit is
                                      ///<   responsible (more than one work unit may be tasked
                                      ///<   with computing any of the relevant energy terms.
                                      ///<   One and only one work unit will be tasked with moving
                                      ///<   each (mobile) atom.  If an atom is not mobile, no
                                      ///<   work unit will be tasked with moving it and terms
                                      ///<   pertaining to it may or may not be computed.
  std::vector<int> angl_term_list;    ///< List of harmonic angle terms to be computed by this
                                      ///<   work unit, indexed into the original topology
  std::vector<int> dihe_term_list;    ///< List of cosine-based dihedral terms to be computed by
                                      ///<   this work unit, indexed into the original topology
  std::vector<int> ubrd_term_list;    ///< List of Urey-Bradley harmonic angle terms to be computed
                                      ///<   by this work unit, indexed into the original topology
  std::vector<int> cimp_term_list;    ///< List of CHARMM harmonic improper dihedral terms to be
                                      ///<   computed by this work unit
  std::vector<int> cmap_term_list;    ///< List of CMAP terms to be computed by this work unit
  std::vector<int> infr14_term_list;  ///< List of inferred 1:4 attenuated interaction terms

  // List of local indices for atoms participating in each task (energy term, virtual site,
  // or constraint group) for which the work unit is responsible
  std::vector<int> bond_i_atoms;      ///< List of I atoms in each harmonic bond
  std::vector<int> bond_j_atoms;      ///< List of J atoms in each harmonic bond
  std::vector<int> angl_i_atoms;      ///< List of local indices for I atoms in each angle term
  std::vector<int> angl_j_atoms;      ///< List of local indices for J atoms in each angle term
  std::vector<int> angl_k_atoms;      ///< List of local indices for K atoms in each angle term
  std::vector<int> dihe_i_atoms;      ///< Local indices for I atoms in each dihedral term
  std::vector<int> dihe_j_atoms;      ///< Local indices for J atoms in each dihedral term
  std::vector<int> dihe_k_atoms;      ///< Local indices for K atoms in each dihedral term
  std::vector<int> dihe_l_atoms;      ///< Local indices for L atoms in each dihedral term
  std::vector<int> ubrd_i_atoms;      ///< Local indices for I atoms in each Urey-Bradley term
  std::vector<int> ubrd_k_atoms;      ///< Local indices for K atoms in each Urey-Bradley term
  std::vector<int> cimp_i_atoms;      ///< Local indices for I atoms in each CHARMM improper term
  std::vector<int> cimp_j_atoms;      ///< Local indices for J atoms in each CHARMM improper term
  std::vector<int> cimp_k_atoms;      ///< Local indices for K atoms in each CHARMM improper term
  std::vector<int> cimp_l_atoms;      ///< Local indices for L atoms in each CHARMM improper term
  std::vector<int> cmap_i_atoms;      ///< Local indices for I atoms in each CMAP term
  std::vector<int> cmap_j_atoms;      ///< Local indices for J atoms in each CMAP term
  std::vector<int> cmap_k_atoms;      ///< Local indices for K atoms in each CMAP term
  std::vector<int> cmap_l_atoms;      ///< Local indices for L atoms in each CMAP term
  std::vector<int> cmap_m_atoms;      ///< Local indices for M atoms in each CMAP term
  std::vector<int> infr14_i_atoms;    ///< Local indices for I atoms in each inferred 1:4 term
  std::vector<int> infr14_l_atoms;    ///< Local indices for L atoms in each inferred 1:4 term

  /// Array of composite bond terms, a simple concatenation of tasks with similar forms
  std::vector<int> cbnd_term_list;

  /// Arrays to support the composite bond list, more simple concatenations
  std::vector<bool> cbnd_is_ubrd;  ///< Indicates that the composite bond term serves a CHARMM
                                   ///<   Urey-bradley interaction, not a harmonic bond
  std::vector<int> cbnd_i_atoms;   ///< Local indices for bond or Urey-Bradley I atoms
  std::vector<int> cbnd_jk_atoms;  ///< Local indices for J atoms (K atoms in Urey-Bradley terms)
  
  /// Array of composite dihedral terms (this will be composed from entries in dihe_term_list and
  /// cimp_term_list as well as the atoms those terms affect, and will then replace those lists).
  std::vector<int2> cdhe_term_list;

  // Arrays to support the composite dihedrals
  std::vector<bool> cdhe_is_cimp;     ///< Indicates that the composite dihedral serves a CHARMM
                                      ///<   improper dihedral, not one or two proper torsions
  std::vector<int> cdhe_i_atoms;      ///< Local indices for I atoms in each dihedral term
  std::vector<int> cdhe_j_atoms;      ///< Local indices for J atoms in each dihedral term
  std::vector<int> cdhe_k_atoms;      ///< Local indices for K atoms in each dihedral term
  std::vector<int> cdhe_l_atoms;      ///< Local indices for L atoms in each dihedral term
  
  
  // Restraint terms for this work unit
  std::vector<int> rposn_term_list;  ///< Positional restraint terms, indexing into the original
                                     ///<   restraint apparatus
  std::vector<int> rbond_term_list;  ///< Distance restraint terms, indexing into the original
                                     ///<   restraint apparatus
  std::vector<int> rangl_term_list;  ///< Three-point angle restraint terms, indexing into the
                                     ///<   original restraint apparatus
  std::vector<int> rdihe_term_list;  ///< Four-point dihedral restraint terms, indexing into the
                                     ///<   original restraint apparatus
  std::vector<int> rposn_atoms;      ///< Local indices for atoms subject to each positional
                                     ///<   restraint term in this work unit
  std::vector<int> rbond_i_atoms;    ///< Local indices for I atoms subject to distance restraints
  std::vector<int> rbond_j_atoms;    ///< Local indices for J atoms subject to distance restraints
  std::vector<int> rangl_i_atoms;    ///< Local indices for I atoms subject to angle restraints
  std::vector<int> rangl_j_atoms;    ///< Local indices for J atoms subject to angle restraints
  std::vector<int> rangl_k_atoms;    ///< Local indices for K atoms subject to angle restraints
  std::vector<int> rdihe_i_atoms;    ///< Local indices for I atoms subject to dihedral restraints
  std::vector<int> rdihe_j_atoms;    ///< Local indices for J atoms subject to dihedral restraints
  std::vector<int> rdihe_k_atoms;    ///< Local indices for K atoms subject to dihedral restraints
  std::vector<int> rdihe_l_atoms;    ///< Local indices for L atoms subject to dihedral restraints

  // Bit string directives on whether to contribute energies to the work unit's reported total
  std::vector<uint> acc_bond_energy;    ///< Accumulator flags for bond energy terms
  std::vector<uint> acc_angl_energy;    ///< Accumulator flags for harmonic angle energy terms
  std::vector<uint> acc_dihe_energy;    ///< Accumulator flags for dihedral energy terms
  std::vector<uint> acc_ubrd_energy;    ///< Accumulator flags for Urey-Bradley energy terms
  std::vector<uint> acc_cimp_energy;    ///< Accumulator flags for CHARMM improper energy terms
  std::vector<uint> acc_cmap_energy;    ///< Accumulator flags for CMAP energy terms
  std::vector<uint> acc_infr14_energy;  ///< Accumulator flags for inferred 1:4 energy terms
  std::vector<uint> acc_rposn_energy;   ///< Accumulator flags for positional restraint energies
  std::vector<uint> acc_rbond_energy;   ///< Accumulator flags for distance restraint energies
  std::vector<uint> acc_rangl_energy;   ///< Accumulator flags for angle restraint energies
  std::vector<uint> acc_rdihe_energy;   ///< Accumulator flags for dihedral restraint energies
  std::vector<uint> acc_cbnd_energy;    ///< Accumulator flags for composite bond energies
  std::vector<uint> acc_cdhe_energy;    ///< Accumulator flags for composite dihedral energies

  // Constraint groups for this work unit
  std::vector<int> cnst_group_list;    ///< List of constraint groups, indexing into the group
                                       ///<   enumerated in the original topology, that this
                                       ///<   work unit is responsible for enforcing (the work
                                       ///<   unit must be responsible for moving all atoms in
                                       ///<   any such constraint group)
  std::vector<int> sett_group_list;    ///< List of fast rigid water SETTLE groups, indexing
                                       ///<   into the groups enumerated in the original
                                       ///<   topology, assigned to this work unit
  std::vector<int> cnst_group_atoms;   ///< Local indices of atoms in all constrained groups.
                                       ///<   The bounds of this list, delineating separate groups,
                                       ///<   are found in cnst_group_bounds.
  std::vector<int> cnst_group_bounds;  ///< Bounds array for cnst_group_atoms
  std::vector<int> sett_ox_atoms;      ///< Local indices of oxygen atoms in each of this work
                                       ///<   unit's SETTLE groups
  std::vector<int> sett_h1_atoms;      ///< Local indices of the first hydrogen atoms in each of
                                       ///<   this work unit's SETTLE groups
  std::vector<int> sett_h2_atoms;      ///< Local indices of the second hydrogen atoms in each of
                                       ///<   this work unit's SETTLE groups

  // Virtual sites in this work unit
  std::vector<int> virtual_site_list;   ///< List of virtual sites, indexing into the original
                                        ///<   topology's list of virtual sites (i.e. the 5th
                                        ///<   virtual site, which could be atom index 19 in a box
                                        ///<   of TIP4P-Eq water)
  std::vector<int> vsite_atoms;         ///< Local indices for each of the respective virtual
                                        ///<   particles in this work unit's above list
  std::vector<int> vsite_frame1_atoms;  ///< Local indices of this work unit's virtual site parent
                                        ///<   atoms
  std::vector<int> vsite_frame2_atoms;  ///< Local indices of this work unit's virtual site second
                                        ///<   frame atoms
  std::vector<int> vsite_frame3_atoms;  ///< Local indices of this work unit's virtual site third
                                        ///<   frame atoms
  std::vector<int> vsite_frame4_atoms;  ///< Local indices of this work unit's virtual site fourth
                                        ///<   frame atoms

  // Internal arrays of instructions for each energy term, constraint, or virtual site task,
  // stored for fast retrieval.
  std::vector<uint2> cbnd_instructions;   ///< Composite bond and Urey-Bradley instructions
  std::vector<uint2> angl_instructions;   ///< Harmonic bond angle instructions
  std::vector<uint3> cdhe_instructions;   ///< Composite dihedral and CHARMM improper instructions
  std::vector<uint2> cmap_instructions;   ///< CMAP surface term instructions
  std::vector<uint> infr14_instructions;  ///< Inferred 1:4 attenuated interaction instructions
  std::vector<uint2> rposn_instructions;  ///< Positional restraint instructions
  std::vector<uint2> rbond_instructions;  ///< Distance restraint instructions
  std::vector<uint2> rangl_instructions;  ///< Three-point angle restraint instructions
  std::vector<uint2> rdihe_instructions;  ///< Four-point dihedral restraint instructions
  std::vector<uint2> vste_instructions;   ///< Virtual site placement instructions
  std::vector<uint2> sett_instructions;   ///< Settle group instructions
  std::vector<uint2> cnst_instructions;   ///< Constraint group instructions
  
  // Pointers to important objects
  ValenceDelegator *vdel_pointer;        ///< The delegator managing this object's creation
  const AtomGraph *ag_pointer;           ///< The topology to which this object pertains
  const RestraintApparatus *ra_pointer;  ///< Restraint apparatus to which this object pertains
};

/// \brief Calculate the optimal sizes of valence work units based on a series of system sizes.
///        This will engage some simple heuristics, not examine the actual topologies to optimize
///        a particular number.
///
/// Overloaded:
///   - Accept a C-style array with the system sizes and a trusted length 
///   - Accept a Standard Template Library vector of the system sizes
///   - Accept a Hybrid object with the system sizes
///   - Accept launch parameters to inform the choice of work unit size
///
/// \param atom_counts   The sizes of each system
/// \param system_count  The number of systems (if a C-style array is provided)
/// \param sm_count      The number of blocks that will be launched on the GPU (based on the number
///                      of streaming multiprocessors)
/// \param kwidth        The recommended kernel width in which to launch the valence calculations
///                      (this will be large enough to accommodate the largest work units, but
///                      could be larger if the card is under-filled).  The recommended thread
///                      block size is returned through this pointer.
/// \{
int calculateValenceWorkUnitSize(const int* atom_counts, int system_count);

int calculateValenceWorkUnitSize(const std::vector<int> &atom_counts);

int calculateValenceWorkUnitSize(const Hybrid<int> &atom_counts);

int calculateValenceWorkUnitSize(const int* atom_counts, int system_count, int sm_count,
                                 ValenceKernelSize *kwidth);

int calculateValenceWorkUnitSize(const std::vector<int> &atom_counts, int sm_count,
                                 ValenceKernelSize *kwidth);

int calculateValenceWorkUnitSize(const Hybrid<int> &atom_counts, int sm_count,
                                 ValenceKernelSize *kwidth);
/// \}

/// \brief Build a series of valence work units to cover a topology.
///
/// Overloaded:
///   - Accept a pre-existing ValenceDelegator, so that information used in creating the work units
///     will be available after they are finished
///   - Construct the ValenceDelegator internally (some functionality associated with the
///     ValenceWorkUnits will be unavailable later)
///
/// \param ag                 The topology of interest
/// \param ra                 Restraints linked to the topology
/// \param vdel               Object for managing the creation of the work units
/// \param max_atoms_per_vwu  The maximum number of atoms to permit in each work unit
/// \{
std::vector<ValenceWorkUnit>
buildValenceWorkUnits(const AtomGraph *ag, const RestraintApparatus *ra,
                      int max_atoms_per_vwu = maximum_valence_work_unit_atoms);

std::vector<ValenceWorkUnit>
buildValenceWorkUnits(ValenceDelegator *vdel,
                      int max_atoms_per_vwu = maximum_valence_work_unit_atoms);
/// \}

} // namespace topology
} // namespace stormm

#endif
