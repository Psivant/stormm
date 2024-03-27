// -*-c++-*-
#ifndef STORMM_MATCH_BONDING_PATTERN_H
#define STORMM_MATCH_BONDING_PATTERN_H

#include <vector>
#include "copyright.h"
#include "DataTypes/common_types.h"
#include "Topology/atomgraph.h"
#include "chemistry_enumerators.h"

namespace stormm {
namespace chemistry {

using topology::AtomGraph;
  
/// \brief Beginning with two distinct atoms in a molecule within a topology, proceed throughout
///        the molecular bonding pattern verifying that the atoms one encounters when stepping
///        outward from each atom have identical properties.
///
/// Overloaded:
///   - Pass the topology by pointer or by reference (this is an optimization for when the
///     function gets heavy usage--with its member variables, the topology constitutes an object of
///     significant size)
///   - Accept a collection of pre-allocated vectors to use a scratch space
///   - Accept a ChemicalFeatures object from which to extract critical information
///   - Allocate temporary storage space for all work
///
/// \param ag              System topology
/// \param formal_charges  Array of formal charges for all atoms in the entire topology
/// \param bond_orders     Array of bond orders for all bonds in the topology
/// \param free_electrons  Array of ree electron content for all atoms in the entire topology
/// \param ring_inclusion  Array of ring inclusion marks for all atoms in the entire topology
/// \param chiralities     Array of chiral statuses for all atoms in the entire topology
/// \param atom_a          The first atom to compare   
/// \param atom_b          The second atom to compare
/// \param a_idx_tree      Tree of topological atom indices branching out from the first (A) atom
/// \param b_idx_tree      Tree of topological atom indices branching out from the second (B) atom
/// \param a_zn_tree       Tree of atomic numbers found when branching from the first atom
/// \param b_zn_tree       Tree of atomic numbers found when branching from the second atom
/// \param a_fc_tree       Tree of formal charges branching from the first atom
/// \param b_fc_tree       Tree of formal charges branching from the second atom
/// \param a_fe_tree       Tree of free electron contents branching from the first atom
/// \param b_fe_tree       Tree of free electron contents branching from the second atom
/// \param a_ri_tree       Tree of ring inclusions branching from the first atom
/// \param b_ri_tree       Tree of ring inclusions branching from the second atom
/// \param a_ch_tree       Tree of chiral designations branching from the first atom
/// \param b_ch_tree       Tree of chiral designations branching from the second atom
/// \param a_coverage      Coverage of the first atom's tree
/// \param b_coverage      Coverage of the second atom's tree
/// \{
bool matchBondingPattern(const AtomGraph *ag, const std::vector<double> &formal_charges,
                         const std::vector<double> &free_electrons,
                         const std::vector<ullint> &ring_inclusion,
	                 const std::vector<ChiralOrientation> &chiralities, const int atom_a,
                         const int atom_b, std::vector<int> *a_idx_tree,
                         std::vector<int> *b_idx_tree, std::vector<int> *a_zn_tree,
                         std::vector<int> *b_zn_tree, std::vector<double> *a_fc_tree,
                         std::vector<double> *b_fc_tree, std::vector<double> *a_fe_tree,
                         std::vector<double> *b_fe_tree, std::vector<ullint> *a_ri_tree,
                         std::vector<ullint> *b_ri_tree,
                         std::vector<ChiralOrientation> *a_ch_tree,
                         std::vector<ChiralOrientation> *b_ch_tree,
                         std::vector<int> *a_coverage, std::vector<int> *b_coverage);

bool matchBondingPattern(const AtomGraph &ag, const std::vector<double> &formal_charges,
                         const std::vector<double> &free_electrons,
                         const std::vector<ullint> &ring_inclusion,
                         const std::vector<ChiralOrientation> &chiralities, const int atom_a,
                         const int atom_b, std::vector<int> *a_idx_tree,
                         std::vector<int> *b_idx_tree, std::vector<int> *a_zn_tree,
                         std::vector<int> *b_zn_tree, std::vector<double> *a_fc_tree,
                         std::vector<double> *b_fc_tree, std::vector<double> *a_fe_tree,
                         std::vector<double> *b_fe_tree, std::vector<ullint> *a_ri_tree,
                         std::vector<ullint> *b_ri_tree,
                         std::vector<ChiralOrientation> *a_ch_tree,
                         std::vector<ChiralOrientation> *b_ch_tree,
                         std::vector<int> *a_coverage, std::vector<int> *b_coverage);

bool matchBondingPattern(const AtomGraph *ag, const std::vector<double> &formal_charges,
                         const std::vector<double> &free_electrons,
                         const std::vector<ullint> &ring_inclusion,
                         const std::vector<ChiralOrientation> &chiralities, int atom_a,
                         int atom_b);

bool matchBondingPattern(const AtomGraph &ag, const std::vector<double> &formal_charges,
                         const std::vector<double> &free_electrons,
                         const std::vector<ullint> &ring_inclusion,
                         const std::vector<ChiralOrientation> &chiralities, int atom_a,
                         int atom_b);
/// \}

} // namespace chemistry
} // namespace stormm

#endif
