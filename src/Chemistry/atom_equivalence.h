// -*-c++-*-
#ifndef STORMM_ATOM_EQUIVALENCE_H
#define STORMM_ATOM_EQUIVALENCE_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "UnitTesting/stopwatch.h"
#include "chemical_features.h"
#include "chemistry_enumerators.h"

namespace stormm {
namespace chemistry {

using topology::AtomGraph;
using testing::StopWatch;

/// \brief Compute ranks for each atom from one or more molecules, assigning a unique (but,
///        essentially arbitrary) integer to each atom based on whether its connections to other
///        atoms in the system can be deemed unique.  If the matchBondingPattern() function below
///        returns TRUE for atoms i and j, then atoms i and j will have the same integer rank in
///        this object.
class AtomRank {
public:

  /// \brief The atom ranks can be constructed with or without temporary arrays to use as a
  ///        workspace.  This optimization can help conserve memory allocations in the context of
  ///        other functions, e.g. constructing an AtomEquivalence object, which use much of the
  ///        same storage space.
  /// \{
  AtomRank(const AtomGraph *ag_in = nullptr);
  
  AtomRank(const AtomGraph *ag_in, const std::vector<double> &formal_charges,
           const std::vector<double> &free_electrons, const std::vector<ullint> &ring_inclusion,
           const std::vector<ChiralOrientation> &chiralities, std::vector<int> *a_idx_tree,
           std::vector<int> *b_idx_tree, std::vector<int> *a_zn_tree, std::vector<int> *b_zn_tree,
           std::vector<double> *a_fc_tree, std::vector<double> *b_fc_tree,
           std::vector<double> *a_fe_tree, std::vector<double> *b_fe_tree,
           std::vector<ullint> *a_ri_tree, std::vector<ullint> *b_ri_tree,
           std::vector<ChiralOrientation> *a_ch_tree, std::vector<ChiralOrientation> *b_ch_tree,
           std::vector<int> *a_coverage, std::vector<int> *b_coverage, int low_molecule_index = 0,
           int high_molecule_index = -1);

  AtomRank(const AtomGraph &ag_in, const std::vector<double> &formal_charges,
           const std::vector<double> &free_electrons, const std::vector<ullint> &ring_inclusion,
           const std::vector<ChiralOrientation> &chiralities, int low_molecule_index = 0,
           int high_molecule_index = -1);

  AtomRank(const AtomGraph *ag_in, const ChemicalFeatures &chemfe, int low_molecule_index = 0,
           int high_molecule_index = -1);

  AtomRank(const AtomGraph &ag_in, const ChemicalFeatures &chemfe, int low_molecule_index = 0,
           int high_molecule_index = -1);
  /// \}

  /// \brief Get the rank of one or more atoms.
  ///
  /// Overload:
  ///   - Get a const reference to the vector of all atoms' ranks
  ///   - Get a copy of the ranks for a segment of the atoms in the molecule
  ///   - Get a particular atom's rank  
  ///
  /// \param low_index   Lower limit of atoms whose ranks to query
  /// \param high_index  Upper limit of atoms whose ranks to query
  /// \param atom_index  Topological index of the atom of interest
  /// \{
  const std::vector<int>& getRanks() const;
  std::vector<int> getRanks(int low_index, int high_index) const;
  int getRank(int atom_index) const;
  /// \}

  /// \brief Get a list of all atoms in the system with a particular rank.
  ///
  /// \param rank_value  The rank of interest
  std::vector<int> getAtomsWithRank(const int rank_value) const;

  /// \brief Get a list of all atoms which have a similar rank to a specific atom.
  ///
  /// \param atom_index  The atom of interest
  std::vector<int> getRankPartners(const int atom_index) const;
  
private:
  int maximum_rank;
  std::vector<int> ranks;
  std::vector<int> rank_degeneracy_bounds;
  std::vector<int> rank_instances;
  AtomGraph *ag_pointer;
};
  
/// \brief Collect the set of symmetry-related and interchangeable groups for one or more molecules
///        within a topology.
class AtomEquivalence {
public:

  /// \brief The constructor relies on a topology pointer plus arrays of formal charges, bond
  ///        orders, and ring inclusion as will have been derived by a ChemicalFeatures object.
  ///
  /// \param ag_in                Topology for the system of interest (a pointer will be retained
  ///                             by the object)
  /// \param formal_charges       Array of formal charges for all atoms in the entire topology
  /// \param free_electrons       Array of free electron content for all atoms in the topology
  /// \param ring_inclusion       Array of ring inclusion marks for all atoms in the topology
  /// \param chiralities          Array of chiral statuses for all atoms in the entire topology
  /// \param chemfe_in            Chemical features of the molecular system in question, from which
  ///                             the previous four arrays as well as the original topology pointer
  ///                             can be obtained
  /// \param timer                Wall time tracker to monitor time spent in various stages of the
  ///                             calculation
  /// \param low_molecule_index   Index of the first molecule in the topology for which to draw
  ///                             atom equivalence groups
  /// \param high_molecule_index  Upper index limit of molecules in the topology for which to draw
  ///                             atom equivalence groups.  The default of -1 triggers the upper
  ///                             limit to be set as one beyond the lower limit.
  /// \{
  AtomEquivalence();
  
  AtomEquivalence(const AtomGraph *ag_in, const std::vector<double> &formal_charges,
                  const std::vector<double> &free_electrons,
                  const std::vector<ullint> &ring_inclusion,
                  const std::vector<ChiralOrientation> &chiralities, StopWatch *timer = nullptr,
                  int low_molecule_index = 0, int high_molecule_index = -1);

  AtomEquivalence(const AtomGraph &ag_in, const std::vector<double> &formal_charges,
                  const std::vector<double> &free_electrons,
                  const std::vector<ullint> &ring_inclusion,
                  const std::vector<ChiralOrientation> &chiralities, StopWatch *timer = nullptr,
                  int low_molecule_index = 0, int high_molecule_index = -1);
  
  AtomEquivalence(const ChemicalFeatures &chemfe_in, StopWatch *timer = nullptr,
                  int low_molecule_index = 0, int high_molecule_index = -1);

  AtomEquivalence(const AtomGraph &ag_in, const CoordinateFrame &cf, StopWatch *timer = nullptr,
                  int low_molecule_index = 0, int high_molecule_index = -1);
  /// \}

  /// The default copy and move constructors, as well as copy and move assignment operators, will
  /// be adequate for this object which contains only Standard Template Library components and no
  /// const member variables.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object to which the current one shall be assigned
  /// \{
  AtomEquivalence(const AtomEquivalence &original) = default;
  AtomEquivalence(AtomEquivalence &&original) = default;
  AtomEquivalence& operator=(const AtomEquivalence &original) = default;
  AtomEquivalence& operator=(AtomEquivalence &&original) = default;
  /// \}
  
  /// \brief Get the number of atom equivalence groups.
  int getGroupCount() const;

  /// \brief Get the number of levels of symmetry-related groups.  This is essential for mapping
  ///        identity swaps that an atom might undergo in a symmetry-corrected RMSD calculation.
  int getSymmetryDepth() const;
  
  /// \brief Get the atoms of one group as an independent Standard Template Library vector.
  ///
  /// \param group_index  The group of interest
  std::vector<int> getGroup(int group_index) const;

  /// \brief Get one atom index from within a symmetry-equivalent atom group.
  ///
  /// \param group_index   The group of interest
  /// \param domain_index  The interchangeable collection of atoms within the group
  /// \param atom_index    Index of the atom within the interchangeable domain (this is not a
  ///                      topological index)
  int getSymmetryRelatedAtom(int group_index, int domain_index, int atom_index) const;

  /// \brief Get a pointer to the atom indices of a particular group.
  ///
  /// \param group_index  The group of interest  
  const int* getGroupPointer(int group_index) const;

  /// \brief Get the size of a particular group.
  ///
  /// \param group_index  The group of interest  
  int getGroupSize(int group_index) const;

  /// \brief Get the order of a particular group
  ///
  /// \param group_index  The group of interest  
  int getGroupOrder(int group_index) const;

  /// \brief Get the depth of the group of interest in terms of its dependence on (containment
  ///        within) other groups of symmetry-related atoms.
  ///
  /// \param group_index  The group of interest  
  int getGroupLevel(int group_index) const;
  
  /// \brief Get the manner by which combinations of the group's interchangeable subunits should be
  ///        laid out in order to search for the best RMSD fit.
  ///
  /// \param group_index  The group of interest  
  EquivalenceSwap getGroupRule(int group_index) const;
  
  /// \brief Get a list of a group's dependencies (other groups that the group in question contains
  ///        completely).  Chained dependencies are handled implicitly by the complete containment
  ///        criterion.
  ///
  /// \param group_index  The group of interest  
  std::vector<int> getGroupDependencies(int group_index) const;

  /// \brief Get the number of asymmetric core atoms
  int getAsymmetricAtomCount() const;
  
  /// \brief Get a const reference to the list of atoms outside any symmetry group.
  const std::vector<int>& getAsymmetricAtoms() const;

  /// \brief Get a const pointer to the topology for which these equivalencies apply.
  const AtomGraph* getTopologyPointer() const;

private:
  int group_count;                     ///< The number of equivalent atom groups
  int symmetry_depth;                  ///< The number of levels in symmetry groups and their
                                       ///<   dependencies 
  std::vector<int> group_atoms;        ///< Indices of atoms in each equivalent group.  For a group
                                       ///<   with K atoms and N-fold symmetry, the equivalent
                                       ///<   atoms of the first group instance are listed
                                       ///<   [ A(1,i), A(1,i+1), ..., A(1,k) ] followed by those
                                       ///<   of the second instance [ A(2,i), ..., A(2,k) ] and
                                       ///<   so on up to [ A(N,i), ..., A(N,k) ].
  std::vector<int> group_sizes;        ///< The number of atoms in each equivalent atom group
  std::vector<int> group_orders;       ///< The number of symmetric instances of each group
  std::vector<int> group_bounds;       ///< Bounds array for group_atoms above
  std::vector<int> group_levels;       ///< The depth at which each group is accessed.  Groups that
                                       ///<   are not the dependents of any other are level 0 (and
                                       ///<   if all groups are at level 0 the object's overall
                                       ///<   symmetry_depth is 1).  Symmetry-related atom groups
                                       ///<   that are the direct dependents of level 0 groups are
                                       ///<   level 1, their dependents are level 2, and so on.
  std::vector<int> dependencies;       ///< Indices of dependent groups for each group in the
                                       ///<   molecule.  If group 4 is completely subsumed within
                                       ///<   group 3, then the dependencies of group 3 will
                                       ///<   include group 4.
  std::vector<int> dependency_bounds;  ///< Bounds array for dependencies above, group_count + 1 in
                                       ///<   length
  std::vector<int> asymmetric_atoms;   ///< List of atoms that participate in no symmetry group

  /// The manner by which each group can exchange the coordinates of its symmetric subunits
  std::vector<EquivalenceSwap> group_rules;

  /// Pointer to the topology this is describing
  AtomGraph *ag_pointer;

  /// \brief Find equivalent atoms from within a subset of the topology, then call the
  ///        drawEquivalentGroups() function below.
  ///
  /// \param subset_atoms      Topological indices for the subset of atoms that this call will look
  ///                          at, i.e. { 4, 5, 6, 7, 21, 22, 23, 24 }.
  /// \param map_to_subset     Mapping of each topological atom its place in the subset, i.e. for
  ///                          the example above, { x, x, x, x, 0, 1, 2, 3, x, x, x, ..., 4, 5,
  ///                          6, 7 }.  Not all indices of this array can be trusted at all times.
  ///                          If the atom is not actually in the subset, then the map will not
  ///                          contain a valid index into the subset_atoms array.
  /// \param formal_charges    Array of formal charges for all atoms in the entire topology
  /// \param free_electrons    Array of free electron content for all atoms in the entire topology
  /// \param ring_inclusion    Array of ring inclusion marks for all atoms in the entire topology
  /// \param chiralities       Array of chiral statuses for all atoms in the entire topology
  /// \param atom_ranks        Ranks of all atoms in the topology.  Ranks are only required for the
  ///                          atoms in the molecule range of interest.  This vector may contain -1
  ///                          entries for atoms outside of the relevant set.
  /// \param domain_coverage   Array (with space for all topological atom indices) for recording
  ///                          the extent to which various symmetry-related groups expanding
  ///                          outwards from each partner have covered the available atoms
  /// \param allowed_atoms     Array (also with space for all topological atom indices) for
  ///                          recording which atoms of the topology are to get consideration for
  ///                          building the symmetry-equivalent partner domains that will form
  ///                          the basis for symmetry groups
  /// \param candidate_hopper  Space reserved for atoms that might be found to contribute to the
  ///                          next layer of each symmetry-related group
  /// \param jumbled_groups    List of non-overlapping atoms making up each symmetry-equivalent
  ///                          group.  This list is filtered by passing through the
  ///                          candidate_hopper array, eliminating atoms that result in
  ///                          double-coverage.
  /// \param aligned_groups    Re-arrangement of jumbled groups that places atoms from each
  ///                          symmetry-equivalent domain, layer by layer working outward from the
  ///                          symmetric partner atom, into a contiguous stretch of memory
  void findEquivalentAtoms(const std::vector<int> &subset_atoms, std::vector<int> *map_to_subset,
                           const std::vector<double> &formal_charges,
                           const std::vector<double> &free_electrons,
                           const std::vector<ullint> &ring_inclusion,
                           const std::vector<ChiralOrientation> &chiralities,
                           const AtomRank &arnks, std::vector<int> *domain_coverage,
                           std::vector<int> *allowed_atoms, std::vector<int> *candidate_hopper,
                           std::vector<int> *domain_assignment, std::vector<int> *jumbled_groups,
                           std::vector<int> *aligned_groups, std::vector<int> *layer_bounds);

  /// \brief Draw the equivalent groups for a molecule, based on a trusted list of partners each
  ///        representing a distinct collection of symmetry-related atoms.  The strategy is to
  ///        begin with a list of atoms that have been found to have equivalent bonding patterns:
  ///        starting from each atom and exploring outwards throughout the molecule, the elements,
  ///        bonds, and electronic structures of each atom are identical.  The search expands
  ///        outwards from these atoms, accumulating symmetry-equivalent "domains" for each partner
  ///        atom, until further expansion would cause the domains to overlap.  If two partners'
  ///        domains both attempt to include the same atom at one time, it will become part of
  ///        neither domain.  Each call to this function may make more than one group of equivalent
  ///        atoms given the set of partners: only those domains which come into contact with one
  ///        another will be considered interchangeable.  If there are symmetric partner atoms A,
  ///        B, C, and D, with A and B at one end of the molecule and C and D at the other, such
  ///        the domains of A and B or those of C and D might collide, and be able to expand no
  ///        further, well before the domains of A and C or those of B and D ever come into
  ///        contact.  In such a case, the domains of A and B would be considered interchangeable,
  ///        and those of C and D as well, but A and B would not be allowed to exchange with
  ///        C or D.
  ///
  /// Argument descriptions follow from findEquivalentAtoms() above (many of the arrays are meant
  /// to be passed on to calls there), with the additions of:
  ///
  /// \param partners          A list of equivalent atoms with one entry spanning each equivalent
  ///                          group
  /// \param subset_atoms      A list of atoms within which the equivalence groups can be drawn
  /// \param layer_bounds      Working array for the bounds of each layer forming the domains of
  ///                          equivalent atom groups
  void drawEquivalentGroups(const std::vector<int> &partners, const std::vector<int> &subset_atoms,
                            const std::vector<double> &formal_charges,
                            const std::vector<double> &free_electrons,
                            const std::vector<ullint> &ring_inclusion,
                            const std::vector<ChiralOrientation> &chiralities,
                            const AtomRank &arnks, std::vector<int> *domain_coverage,
                            std::vector<int> *allowed_atoms, std::vector<int> *candidate_hopper,
                            std::vector<int> *domain_assignment, std::vector<int> *jumbled_groups,
                            std::vector<int> *aligned_groups, std::vector<int> *layer_bounds);

  /// \brief Add a selection of the symmetry-equivalent domains as an interchangeable collection
  ///        of atom groups, with a plan for how to do the swapping.
  ///
  /// \param all_domains       The collection of atom topological indices for all symmetry-related
  ///                          domains of interest.  This may include more domains than just the
  ///                          selected few which can be swapped.  If a dumbell-shaped molecule has
  ///                          three symmtry-equivalent moieties on either end, the all_domains
  ///                          array will contain the atoms of all six moieties, but only three
  ///                          will be listed in selected_domains.
  /// \param domain_size       The size of each domain.  Each domain is presented as a contiguous
  ///                          stream of atm indices in all_domains.
  /// \param selected_domains  The domains which are free to swap
  /// \param plan              Prescribed method for interchanging domains
  void addSymmetryGroup(const std::vector<int> &all_domains, int domain_size,
                        const std::vector<int> &selected_domains, EquivalenceSwap plan);
  
  /// \brief Construct a symmetry group in which free interchanges of all of its subunits are
  ///        permissible.  Descriptions of parameters follow from addSymmetryGroup() above.
  void freeForAllSymmetryGroup(const std::vector<int> &all_domains, int domain_size,
                               const std::vector<int> &selected_domains);

  /// \brief Construct a symmetry group in which interchanges of the subunits ("domains)" are
  ///        permitted in the manner K -> K+1, K+1 -> K+2, K -> K-1, etc., with N total
  ///        permutations for a ring of N subunits.  This function will first verify that there is
  ///        a ring structure defined by the connectivity between subunits, and if not the group
  ///        will be approximated as a "free for all" symmetry group.  Descriptions of parameters
  ///        follow from addSymmetryGroup(), with the addition of:
  ///
  /// \param touch_table  Table indicating whether any pair of domains / subunits touch one another
  void rotarySymmetryGroup(const std::vector<int> &all_domains, int domain_size,
                           const std::vector<int> &selected_domains,
                           const std::vector<bool> &touch_table);

  /// \brief Cull duplicate groups (those containing the same atoms and only the same atoms) from
  ///        the list of symmetry groups to be considered.
  void cullDuplicateGroups();
  
  /// \brief Analyze the detected symmetry groups and determine which ones are subsumed within
  ///        others.  When making plans to create 
  void findDependencies();

  /// \brief Order the symmetry-related atom groups in decreasing order of the ratio of total mass
  ///        in the group's atoms to the combinatorial number of possible arrangements for those
  ///        atoms, given the group's symmetry order and that of all its dependents.
  void orderAllGroups();
};

/// \brief Assign unique ranks to each atom, with an implicit behavior that ties are not broken.
///        Atoms identified to have the same rank in the result will be equivalent in the sense
///        that they belong to distinct groups of atoms which can be interchanged to improve a
///        symmetry-aware RMSD result.  Descriptions of various modifiable vectors (all sized to
///        accommodate the number of atoms in the molecular system) follow from
///        matchBondingPattern() above.
///
/// Overloaded:
///   - Supply a const pointer or const reference to the topology
///   - Supply the arrays of atom properties directly
///   - Supply a chemical features object from which these arrays will be copied
///
/// \param ag                   System topology
/// \param formal_charges       Array of formal charges for all atoms in the entire topology
/// \param bond_orders          Array of bond orders for all bonds in the topology
/// \param free_electrons       Array of ree electron content for all atoms in the entire topology
/// \param ring_inclusion       Array of ring inclusion marks for all atoms in the entire topology
/// \param chiralities          Array of chiral statuses for all atoms in the entire topology
/// \param low_molecule_index   Index of the first molecule in the topology for which to draw
///                             atom equivalence groups
/// \param high_molecule_index  Upper index limit of molecules in the topology for which to draw
///                             atom equivalence groups.  The default of -1 triggers the upper
///                             limit to be set as one beyond the lower limit.
/// \{
std::vector<int> rankAtoms(const AtomGraph *ag, const std::vector<double> &formal_charges,
                           const std::vector<double> &free_electrons,
                           const std::vector<ullint> &ring_inclusion,
                           const std::vector<ChiralOrientation> &chiralities,
                           std::vector<int> *a_idx_tree, std::vector<int> *b_idx_tree,
                           std::vector<int> *a_zn_tree, std::vector<int> *b_zn_tree,
                           std::vector<double> *a_fc_tree, std::vector<double> *b_fc_tree,
                           std::vector<double> *a_fe_tree, std::vector<double> *b_fe_tree,
                           std::vector<ullint> *a_ri_tree, std::vector<ullint> *b_ri_tree,
                           std::vector<ChiralOrientation> *a_ch_tree,
                           std::vector<ChiralOrientation> *b_ch_tree,
                           std::vector<int> *a_coverage, std::vector<int> *b_coverage,
                           int low_molecule_index = 0, int high_molecule_index = -1);

std::vector<int> rankAtoms(const AtomGraph &ag, const std::vector<double> &formal_charges,
                           const std::vector<double> &free_electrons,
                           const std::vector<ullint> &ring_inclusion,
                           const std::vector<ChiralOrientation> &chiralities,
                           int low_molecule_index = 0, int high_molecule_index = -1);

std::vector<int> rankAtoms(const AtomGraph *ag, const ChemicalFeatures &chemfe,
                           int low_molecule_index = 0, int high_molecule_index = -1);

std::vector<int> rankAtoms(const AtomGraph &ag, const ChemicalFeatures &chemfe,
                           int low_molecule_index = 0, int high_molecule_index = -1);
/// \}

} // namespace chemistry
} // namespace stormm

#endif
