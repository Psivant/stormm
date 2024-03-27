// -*-c++-*-
#ifndef STORMM_COMBOGRAPH_LJMODEL_H
#define STORMM_COMBOGRAPH_LJMODEL_H

#include <string>
#include <vector>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "Parsing/parse.h"
#include "Potential/energy_enumerators.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "atomgraph.h"
#include "atomgraph_abstracts.h"
#include "atomgraph_analysis.h"
#include "lennard_jones_analysis.h"

namespace stormm {
namespace topology {

using energy::getEnumerationName;
using data_types::operator==;
using synthesis::AtomGraphSynthesis;
using topology::AtomGraph;
  
/// \brief An object to manage the combination of Lennard-Jones parameters from two topologies.
///        This can be used to insert "NBFix" pair-specific (off-rule) parameters where neither
///        topology by itself might contain the requisite information.  It can also be used to
///        build parameter tables for a union of the two underlying topologies and produce arrays
///        indices into that table for a new topology.
class ComboGraphLJModel {
public:

  /// \brief The constructor accepts two topologies or two Lennard-Jones parameter sets.  Edits
  ///        to the Lennard-Jones parameters can be provided to the constructor or entered later.
  ///
  /// \param poly_ag_secondary  A synthesis of topologies form which to harvest unique
  ///                           Lennard-Jones types
  /// \param edits_in           A list of edits to apply to modify pair-specific Lennard-Jones
  //                            interactions
  /// \{
  ComboGraphLJModel(const AtomGraph *primary_ag_in, const AtomGraph *ag_secondary_in = nullptr,
                    VdwCombiningRule default_rule_in = VdwCombiningRule::LORENTZ_BERTHELOT,
                    const std::vector<PairLJInteraction> &edits = {});

  ComboGraphLJModel(const AtomGraph &primary_ag_in, const AtomGraph &ag_secondary_in,
                    VdwCombiningRule default_rule_in = VdwCombiningRule::LORENTZ_BERTHELOT,
                    const std::vector<PairLJInteraction> &edits = {});

  ComboGraphLJModel(const AtomGraph &primary_ag_in, const AtomGraphSynthesis &poly_ag_secondary,
                    VdwCombiningRule default_rule_in = VdwCombiningRule::LORENTZ_BERTHELOT,
                    const std::vector<PairLJInteraction> &edits = {});
  
  ComboGraphLJModel(const AtomGraph *primary_ag_in, const std::vector<double> &lj_a_in,
                    const std::vector<double> &lj_b_in,
                    const std::vector<char4> &lj_type_names = {},
                    VdwCombiningRule default_rule_in = VdwCombiningRule::LORENTZ_BERTHELOT,
                    const std::vector<PairLJInteraction> &edits = {});

  ComboGraphLJModel(const AtomGraph &primary_ag_in, const std::vector<double> &lj_a_in,
                    const std::vector<double> &lj_b_in,
                    const std::vector<char4> &lj_type_names = {},
                    VdwCombiningRule default_rule_in = VdwCombiningRule::LORENTZ_BERTHELOT,
                    const std::vector<PairLJInteraction> &edits = {});
  /// \}

  /// \brief Get the number of Lennard-Jones parameter matrix sets.
  int getSetCount() const;

  /// \brief Get the default combining rule for mixing parameters of different topologies.
  VdwCombiningRule getDefaultCombiningRule() const;

  /// \brief Get the combining rule effective in the primary topology.
  VdwCombiningRule getPrimaryTopologyRule() const;

  /// \brief Get the combining rule effective in one of the secondary topologies.
  ///
  /// \param index  The secondarry topology of interest (the ith secondary topology takes part in
  ///               the ith combination with the primary topology)
  VdwCombiningRule getSecondaryTopologyRule(int index) const;

  /// \brief Get the number of atom types in the primary topology.

  /// \brief Get the size of a particular combination matrix.  The matrix has as many rows as the
  ///        primary topology has atom types and as many columns as the new (other) topology has
  ///        atom types.
  ///
  /// \param index  The combination of interest
  int2 getMatrixSize(int index) const;

  /// \brief Get the matrix of Lennard-Jones A coefficients between the primary topology and
  ///        another topology.  The primary atom types control rows of the matrix, while the other
  ///        topology's atom types control columns.  The matrix has column format, like others in
  ///        STORMM.
  ///
  /// \param index  The combination of interest
  const std::vector<double>& getACoefficients(int index) const;
  
  /// \brief Get the matrix of Lennard-Jones B coefficients between the primary topology and
  ///        another topology.  The primary atom types control rows of the matrix, while the other
  ///        topology's atom types control columns.  The matrix has column format, like others in
  ///        STORMM.
  ///
  /// \param index  The combination of interest
  const std::vector<double>& getBCoefficients(int index) const;

  std::vector<int> computeConsensusParameters() const;
    
  /// \brief Add an interaction matrix to the list of combinations.
  ///
  /// Overloaded:
  ///   - Provide a combination of vectors or C-style arrays for the new parameter set
  ///   - Provide the new topology by const pointer or const reference
  ///
  /// \param lja_in         Array of Lennard-Jones A coefficients, i.e. U = A/r^12 + B/r^6, to be
  ///                       combined per the object's default rule
  /// \param ljb_in         Array of Lennard-Jones B coefficients, to be combined per the object's
  ///                       default rule
  /// \param lj_type_count  The number of Lennard-Jones types, implying a trusted length for lja_in
  ///                       and ljb_in, if provided as C-style arrays
  /// \param ag_new         The topology to incorporate.  A null pointer will abort the inclusion.
  /// \param edits          A list of pair-specific Lennard-Jones parameters to apply
  /// \{
  template <typename T>
  void addCombination(const T* lja_in, const T* ljb_in, int lj_type_count,
                      const char4* lj_type_names = nullptr,
                      const std::vector<PairLJInteraction> &edits = {});

  template <typename T>
  void addCombination(const std::vector<T> &lja_in, const std::vector<T> &ljb_in,
                      const std::vector<char4> &lj_type_names = {},
                      const std::vector<PairLJInteraction> &edits = {});

  void addCombination(const AtomGraph *ag_secondary,
                      const std::vector<PairLJInteraction> &edits = {});

  void addCombination(const AtomGraph &ag_secondary,
                      const std::vector<PairLJInteraction> &edits = {});
  /// \}
  
  /// \brief Compute the number of unique atom types across all topologies other than the primary.
  ///        This is done by finding unique fingerprints among the columns of the A and B
  ///        coefficient matrices for each combination parameter set.
  int getAtomTypeUnionSize();

  /// \brief Count the number of consensus atom types that would be needed in order to span all
  ///        unique Lennard-Jones interactions in the 
  
private:

  /// The number of Lennard-Jones parameter sets being tracked by this object, and the length of
  /// arrays set_sizes, set_lja, set_ljb, set_edits, and new_ag_pointers
  int set_count;                             

  /// The default combining rule to be applied to interactions in each topology.
  VdwCombiningRule default_rule;
  
  /// The rule governing Lennard-Jones parameter mixing in the original topology.  This will be
  /// inferred from the parameter matrices contained therein.
  VdwCombiningRule primary_topology_rule;

  /// The number of Lennard-Jones types in the primary topology
  int primary_atom_type_count;
  
  /// The numbers of Lennard-Jones atom types in each additional topology.
  std::vector<int> secondary_atom_type_counts;

  /// Rules inferred for Lennard-Jones parameter mixing in all other topologies.  The ith index of
  /// this array corresponds to the secondary topology in the ith combination.
  std::vector<VdwCombiningRule> secondary_topology_rules;
  
  /// The rules by which each combination can be taken to occur.
  std::vector<VdwCombiningRule> set_rules;
  
  /// Lennard-Jones A coefficients for each set, given as a column-format matrix with a number of
  /// rows given in the 'x" member of the set_sizes array and a number of columns given in the "y"
  /// member.
  std::vector<std::vector<double>> set_lja;

  /// Lennard-Jones B coefficients for each parameter matrix
  std::vector<std::vector<double>> set_ljb;

  /// The edits that occurred in each set.  The presence of any edits (a nonzero array length) will
  /// change the combining rule for that set to "NBFIX."
  std::vector<std::vector<PairLJInteraction>> set_edits;

  /// Pointer to the original topology, whose atom types (or types derived thereof) feed into each
  /// row of the resulting parameter matrices.
  AtomGraph *primary_ag_pointer;

  /// Pointers to new topologies, whose atom types (or types derived thereof) appear in the columns
  /// of the resulting parameter matrices.
  std::vector<AtomGraph*> secondary_ag_pointers;

  /// \brief Validate the index of some requested matrix set.
  ///
  /// \param index  The index of interest
  void validateSetIndex(int index) const;
};
  
} // namespace topology
} // namespace stormm

#include "combograph_ljmodel.tpp"

#endif
