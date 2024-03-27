#include "copyright.h"
#include "combograph_ljmodel.h"

namespace stormm {
namespace topology {

using topology::NonbondedKit;
  
//-------------------------------------------------------------------------------------------------
ComboGraphLJModel::ComboGraphLJModel(const AtomGraph *primary_ag_in,
                                     const AtomGraph *ag_secondary_in,
                                     const VdwCombiningRule default_rule_in,
                                     const std::vector<PairLJInteraction> &edits) :
    set_count{0},
    default_rule{default_rule_in},
    primary_topology_rule{VdwCombiningRule::LORENTZ_BERTHELOT},
    primary_atom_type_count{0},
    secondary_atom_type_counts{},
    secondary_topology_rules{},
    set_lja{}, set_ljb{}, set_edits{},
    primary_ag_pointer{const_cast<AtomGraph*>(primary_ag_in)},
    secondary_ag_pointers{}
{
  if (primary_ag_pointer == nullptr) {
    rtErr("A primary topology must be specified in order to create parameter combination "
          "matrices.", "ComboGraphLJModel");
  }

  // The default mixing rule cannot be NBFix (this would imply that there is, in effect, no rule).
  switch (default_rule) {
  case VdwCombiningRule::LORENTZ_BERTHELOT:
  case VdwCombiningRule::GEOMETRIC:
    break;    
  case VdwCombiningRule::NBFIX:
    rtErr("A default combining rule of " + getEnumerationName(default_rule) + " would imply that "
          "there is no combining rule and is therefore invalid.", "ComboGraphLJModel");
  }

  // Infer the Lennard-Jones combination rule of the original topology and get its atom count.
  const NonbondedKit<double> primary_nbk = primary_ag_pointer->getDoublePrecisionNonbondedKit();
  primary_topology_rule = inferCombiningRule(primary_nbk.lja_coeff, primary_nbk.ljb_coeff,
                                          primary_nbk.n_lj_types);
  primary_atom_type_count = primary_nbk.n_lj_types;

  // Return if there is no other topology provided (yet).  If such a topology is available,
  // process it.
  if (ag_secondary_in == nullptr) {
    return;
  }
  addCombination(ag_secondary_in, edits);
}

//-------------------------------------------------------------------------------------------------
ComboGraphLJModel::ComboGraphLJModel(const AtomGraph &primary_ag_in,
                                     const AtomGraph &ag_secondary_in,
                                     const VdwCombiningRule default_rule_in,
                                     const std::vector<PairLJInteraction> &edits) :
    ComboGraphLJModel(primary_ag_in.getSelfPointer(), ag_secondary_in.getSelfPointer(),
                      default_rule_in, edits)
{}

//-------------------------------------------------------------------------------------------------
ComboGraphLJModel::ComboGraphLJModel(const AtomGraph &primary_ag_in,
                                     const AtomGraphSynthesis &poly_ag_secondary,
                                     const VdwCombiningRule default_rule_in,
                                     const std::vector<PairLJInteraction> &edits) :
    ComboGraphLJModel(primary_ag_in.getSelfPointer(), nullptr, default_rule_in, edits)
{
  // Add new matrices of Lennard-Jones parameters for each topology in the synthesis.
  const int ntop = poly_ag_secondary.getUniqueTopologyCount();
  const std::vector<AtomGraph*>& unique_tops = poly_ag_secondary.getUniqueTopologies();
  for (int i = 0; i < ntop; i++) {
    addCombination(unique_tops[i], edits);
  }
}

//-------------------------------------------------------------------------------------------------
ComboGraphLJModel::ComboGraphLJModel(const AtomGraph *primary_ag_in,
                                     const std::vector<double> &lj_a_in,
                                     const std::vector<double> &lj_b_in,
                                     const std::vector<char4> &lj_type_names,
                                     const VdwCombiningRule default_rule_in,
                                     const std::vector<PairLJInteraction> &edits) :
    ComboGraphLJModel(primary_ag_in, nullptr, default_rule_in, edits)
{
  addCombination(lj_a_in, lj_b_in, lj_type_names, edits);
}

//-------------------------------------------------------------------------------------------------
ComboGraphLJModel::ComboGraphLJModel(const AtomGraph &primary_ag_in,
                                     const std::vector<double> &lj_a_in,
                                     const std::vector<double> &lj_b_in,
                                     const std::vector<char4> &lj_type_names,
                                     const VdwCombiningRule default_rule_in,
                                     const std::vector<PairLJInteraction> &edits) :
  ComboGraphLJModel(primary_ag_in.getSelfPointer(), lj_a_in, lj_b_in, lj_type_names,
                    default_rule_in, edits)
{}

//-------------------------------------------------------------------------------------------------
int ComboGraphLJModel::getSetCount() const {
  return set_count;
}

//-------------------------------------------------------------------------------------------------
VdwCombiningRule ComboGraphLJModel::getDefaultCombiningRule() const {
  return default_rule;
}

//-------------------------------------------------------------------------------------------------
VdwCombiningRule ComboGraphLJModel::getPrimaryTopologyRule() const {
  return primary_topology_rule;
}

//-------------------------------------------------------------------------------------------------
VdwCombiningRule ComboGraphLJModel::getSecondaryTopologyRule(const int index) const {
  validateSetIndex(index);
  return secondary_topology_rules[index];
}

//-------------------------------------------------------------------------------------------------
int2 ComboGraphLJModel::getMatrixSize(const int index) const {
  validateSetIndex(index);
  return {primary_atom_type_count, secondary_atom_type_counts[index] };
}

//-------------------------------------------------------------------------------------------------
const std::vector<double>& ComboGraphLJModel::getACoefficients(const int index) const {
  validateSetIndex(index);
  return set_lja[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<double>& ComboGraphLJModel::getBCoefficients(const int index) const {
  validateSetIndex(index);
  return set_ljb[index];
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ComboGraphLJModel::computeConsensusParameters() const {

  // There are at least as many types as appear in the original topology.  Loop over all sets,
  // form the super-matrix of the primary topology and one of the others, and see how it can
  // be reduced.  Repeat this process and accumulate a consensus matrix.  This will necessarily
  // involve mixing an two of the secondary topologies, if there is a list of those, and that will
  // be done using the default mixing rule plus any 
  int n_consensus_types = primary_atom_type_count;
  for (int i = 0; i < set_count; i++) {

  }
}

//-------------------------------------------------------------------------------------------------
void ComboGraphLJModel::validateSetIndex(const int index) const {
  if (index < 0 || index >= set_count) {
    rtErr("Combination index " + std::to_string(index) + " is invalid for a collection of " +
          std::to_string(set_count) + " topology combinations.", "ComboGraphLJModel",
          "validateSetIndex");
  }
}

//-------------------------------------------------------------------------------------------------
void ComboGraphLJModel::addCombination(const AtomGraph *ag_secondary,
                                       const std::vector<PairLJInteraction> &edits) {
  if (ag_secondary == nullptr) {
    return;
  }
  const NonbondedKit<double> secondary_nbk = ag_secondary->getDoublePrecisionNonbondedKit();
  addCombination(secondary_nbk.lja_coeff, secondary_nbk.ljb_coeff, secondary_nbk.n_lj_types);

  // Replace the nullptr appended by the topology-free overload of this function.
  secondary_ag_pointers.back() = const_cast<AtomGraph*>(ag_secondary);
}

//-------------------------------------------------------------------------------------------------
void ComboGraphLJModel::addCombination(const AtomGraph &ag_secondary,
                                       const std::vector<PairLJInteraction> &edits) {
  addCombination(ag_secondary.getSelfPointer(), edits);
}

} // namespace topology
} // namespace stormm
