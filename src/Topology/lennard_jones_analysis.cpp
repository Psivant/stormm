#include <string>
#include "copyright.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Reporting/error_format.h"
#include "atomgraph_abstracts.h"
#include "atomgraph_analysis.h"
#include "lennard_jones_analysis.h"
#include "topology_util.h"

namespace stormm {
namespace topology {

using parse::char4ToString;
using parse::minimalRealFormat;
using parse::NumberFormat;
using parse::realToString;
  
//-------------------------------------------------------------------------------------------------
PairLJInteraction::PairLJInteraction(const char4 type_a_in, const char4 type_b_in,
                                     const double lja_in, const double ljb_in,
                                     const double lja_14_in, const double ljb_14_in) :
    type_a{type_a_in}, type_b{type_b_in}, lja{lja_in}, ljb{ljb_in}, lja_14{lja_14_in},
    ljb_14{ljb_14_in}
{}

//-------------------------------------------------------------------------------------------------
LennardJonesAnalysis::
LennardJonesAnalysis(const NonbondedKit<double> &nbk,
                     const std::vector<std::vector<char4>> &atom_type_aliases_in) :
    lj_type_count{nbk.n_lj_types},
    atom_type_count{0}, set_count{0}, prevalent_rule{}, absolute_rule{}, sigma{}, sigma_14{},
    epsilon{}, epsilon_14{}, lja_coeff{}, ljb_coeff{}, ljc_coeff{}, lja_14_coeff{}, ljb_14_coeff{},
    ljc_14_coeff{}, edits{}, atom_type_aliases{atom_type_aliases_in}, atom_type_map{},
    set_to_consensus_map{}, consensus_to_set_map{}
{
  // Check that there are enough alias lists for each atom type
  if (static_cast<int>(atom_type_aliases.size()) != nbk.n_lj_types) {
    rtErr("One atom type name (as part of a list that may include possible aliases) must be "
          "provided for each Lennard-Jones type.  There were " +
          std::to_string(atom_type_aliases.size()) + " lists of atom type names for " +
          std::to_string(nbk.n_lj_types) + " Lennard-Jones types.", "LennardJonesAnalysis");
  }
  
  // Factor out the sigma and epsilon parameters.  If the general non-bonded parameters are zero,
  // the 1:4 non-bonded interaction parameters will be assumed to be negligible as well.
  sigma.resize(nbk.n_lj_types);
  sigma_14.resize(nbk.n_lj_types);
  epsilon.resize(nbk.n_lj_types);
  epsilon_14.resize(nbk.n_lj_types);
  for (int i = 0; i < lj_type_count; i++) {
    const size_t diag_idx = (i * nbk.n_lj_types) + i;
    if (fabs(nbk.ljb_coeff[diag_idx]) > constants::small) {
      sigma[i] = sqrt(cbrt(nbk.lja_coeff[diag_idx] / nbk.ljb_coeff[diag_idx]));
      sigma_14[i] = sqrt(cbrt(nbk.lja_14_coeff[diag_idx] / nbk.ljb_14_coeff[diag_idx]));
      epsilon[i] = 0.25 * nbk.ljb_coeff[diag_idx] / pow(sigma[i], 6.0);
      epsilon_14[i] = 0.25 * nbk.ljb_14_coeff[diag_idx] / pow(sigma_14[i], 6.0);
    }
    else {
      sigma[i] = 0.0;
      sigma_14[i] = 0.0;
      epsilon[i] = 0.0;
      epsilon_14[i] = 0.0;
    }
  }

  // Group the atom types involved in each Lennard-Jones interaction
  for (int i = 0; i < nbk.n_lj_types; i++) {
    const int n_type_names = static_cast<int>(atom_type_aliases[i].size());
    atom_type_count += n_type_names;
    for (int j = 0; j < n_type_names; j++) {
      const uint ij_key = char4ToUint(atom_type_aliases[i][j]);
      std::map<uint, std::vector<int>>::iterator it = atom_type_map.find(ij_key);
      if (it != atom_type_map.end()) {
        rtErr("Atom type " + char4ToString(atom_type_aliases[i][j]) + " controls multiple "
              "Lennard-Jones parameter sets in the consensus tables.", "LennardJonesAnalysis");
      }
      atom_type_map[ij_key] = std::vector<int>(1, i);
    }
  }
  
  // Extract the Lennard-Jones A and B coefficients directly from the first topology.
  prevalent_rule = inferCombiningRule<double>(nbk.lja_coeff, nbk.ljb_coeff, nbk.n_lj_types,
                                              ExceptionResponse::SILENT, true);
  absolute_rule = inferCombiningRule<double>(nbk.lja_coeff, nbk.ljb_coeff, nbk.n_lj_types,
                                              ExceptionResponse::SILENT, false);
  const size_t nlj_squared = nbk.n_lj_types * nbk.n_lj_types;
  lja_coeff.resize(nlj_squared);
  ljb_coeff.resize(nlj_squared);
  ljc_coeff.resize(nlj_squared);
  lja_14_coeff.resize(nlj_squared);
  ljb_14_coeff.resize(nlj_squared);
  ljc_14_coeff.resize(nlj_squared);
  for (size_t i = 0; i < nlj_squared; i++) {
    lja_coeff[i] = nbk.lja_coeff[i];
    ljb_coeff[i] = nbk.ljb_coeff[i];
    ljc_coeff[i] = nbk.ljc_coeff[i];
    lja_14_coeff[i] = nbk.lja_14_coeff[i];
    ljb_14_coeff[i] = nbk.ljb_14_coeff[i];
    ljc_14_coeff[i] = nbk.ljc_14_coeff[i];
  }

  // Identify any edits that would be characteristic of NBFix.
  for (int i = 1; i < nbk.n_lj_types; i++) {
    for (int j = 0; j < i; j++) {
      const size_t ji_idx = (nbk.n_lj_types * i) + j;
      double ij_eps = epsilon[i] * epsilon[j];
      ij_eps = (ij_eps > constants::tiny) ? sqrt(ij_eps) : 0.0;
      double ij_sig, ij_14_sig;
      switch (prevalent_rule) {
      case VdwCombiningRule::GEOMETRIC:
        ij_sig = sigma[i] * sigma[j];
        ij_14_sig = sigma_14[i] * sigma_14[j];
        ij_sig = (ij_sig > constants::tiny) ? sqrt(ij_sig) : 0.0;
        ij_14_sig = (ij_14_sig > constants::tiny) ? sqrt(ij_14_sig) : 0.0;
        break;
      case VdwCombiningRule::LORENTZ_BERTHELOT:
        ij_sig = 0.5 * (sigma[i] + sigma[j]);
        ij_14_sig = 0.5 * (sigma_14[i] + sigma_14[j]);
        break;
      case VdwCombiningRule::NBFIX:
        edits.push_back(PairLJInteraction(atom_type_aliases[i][0], atom_type_aliases[i][j],
                                          nbk.lja_coeff[ji_idx], nbk.ljb_coeff[ji_idx],
                                          nbk.lja_14_coeff[ji_idx], nbk.ljb_14_coeff[ji_idx]));
        continue;
      }
      const double ij_sig_six = pow(ij_sig, 6.0);
      const double ij_sig_14_six = pow(ij_14_sig, 6.0);
      const double pred_ljb = 4.0 * ij_eps * ij_sig_six;
      const double pred_lja = pred_ljb * ij_sig_six;
      if (pred_ljb < constants::small) {
        if (fabs(pred_lja - nbk.lja_coeff[ji_idx]) > 1.0e-5 ||
            fabs(pred_ljb - nbk.ljb_coeff[ji_idx]) > 1.0e-5) {
          edits.push_back(PairLJInteraction(atom_type_aliases[i][0], atom_type_aliases[i][j],
                                            nbk.lja_coeff[ji_idx], nbk.ljb_coeff[ji_idx],
                                            nbk.lja_14_coeff[ji_idx], nbk.ljb_14_coeff[ji_idx]));
        }
      }
      else {
        if (fabs((pred_lja - nbk.lja_coeff[ji_idx]) / nbk.lja_coeff[ji_idx]) > 1.0e-5 ||
            fabs((pred_ljb - nbk.ljb_coeff[ji_idx]) / nbk.ljb_coeff[ji_idx]) > 1.0e-5) {
          edits.push_back(PairLJInteraction(atom_type_aliases[i][0], atom_type_aliases[i][j],
                                            nbk.lja_coeff[ji_idx], nbk.ljb_coeff[ji_idx],
                                            nbk.lja_14_coeff[ji_idx], nbk.ljb_14_coeff[ji_idx]));
        }
      }
    }
  }

  // Set the origins of each Lennard-Jones type index in the original topology.
  set_to_consensus_map.resize(1);
  set_to_consensus_map[0].resize(nbk.n_lj_types);
  for (int i = 0; i < nbk.n_lj_types; i++) {
    set_to_consensus_map[0][i] = i;
  }
  
  // Initialize the consensus map, an explanation of which Lennard-Jones type index in the
  // consensus tables each Lennard-Jones index in one of the included topologies references.
  consensus_to_set_map.resize(nbk.n_lj_types);
  for (int i = 0; i < nbk.n_lj_types; i++) {
    consensus_to_set_map[i].resize(1);
    consensus_to_set_map[i][0] = { 0, i };
  }
  
  // Note that there is one set of Lennard-Jones parameters in the analysis
  set_count = 1;
}

//-------------------------------------------------------------------------------------------------
int LennardJonesAnalysis::getLJTypeCount() const {
  return lj_type_count;
}

//-------------------------------------------------------------------------------------------------
int LennardJonesAnalysis::getAtomTypeCount() const {
  return atom_type_count;
}

//-------------------------------------------------------------------------------------------------
int LennardJonesAnalysis::getSetCount() const {
  return set_count;
}

//-------------------------------------------------------------------------------------------------
VdwCombiningRule LennardJonesAnalysis::getMostPrevalentCombiningRule() const {
  return prevalent_rule;
}

//-------------------------------------------------------------------------------------------------
const std::vector<char4>& LennardJonesAnalysis::getLJAliases(const int consensus_index) const {
  if (consensus_index < 0 || consensus_index >= lj_type_count) {
    rtErr("Index " + std::to_string(consensus_index) + " is invalid for a collection of " +
          std::to_string(lj_type_count) + " unique Lennard-Jones interactions.",
          "LennardJonesAnalysis", "getLJAliases");
  }
  return atom_type_aliases[consensus_index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<char4>& LennardJonesAnalysis::getLJAliases(const double sigma_query,
                                                             const double epsilon_query,
                                                             const double tolerance) const {

  // Find a match and return immediately, if possible.
  int best_match = -1;
  double best_error = 4.0 * tolerance * tolerance;
  for (int i = 0; i < lj_type_count; i++) {
    if (fabs(sigma[i] - sigma_query) < tolerance && fabs(epsilon[i] - epsilon_query) < tolerance) {
      const double dsig = sigma[i] - sigma_query;
      const double deps = epsilon[i] - epsilon_query;
      const double derr = (dsig * dsig) + (deps * deps);
      if (derr < best_error) {
        best_match = i;
        best_error = derr;
      }
    }
  }
  if (best_match == -1) {
    rtErr("No sigma and epsilon parameter combination " +
          realToString(sigma_query, 7, 4, NumberFormat::STANDARD_REAL) + " (sigma) and " +
          realToString(epsilon_query, 7, 4, NumberFormat::STANDARD_REAL) + " (epsilon) was "
          "found in the consensus tables.");
  }
  return atom_type_aliases[best_match];
}

//-------------------------------------------------------------------------------------------------
const std::vector<char4>&
LennardJonesAnalysis::getLJAliases(const char4 atom_type_query,
                                   const ExceptionResponse policy) const {
  const uint query_key = char4ToUint(atom_type_query);
  const std::map<uint, std::vector<int>>::const_iterator it = atom_type_map.find(query_key);
  if (it == atom_type_map.end()) {
    rtErr("No atom type \"" + char4ToString(atom_type_query) + "\" could be located in the "
          "consensus tables.", "LennardJonesAnalysis", "getLJAliases");
  }
  const std::vector<int>& lj_table_indices = it->second;
  if (lj_table_indices.size() > 1) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("More than one parameter definition is available for atom type " +
            char4ToString(atom_type_query) + ".", "LennardJonesAnalysis", "getLJAliases");
    case ExceptionResponse::WARN:
      rtWarn("More than one parameter definition is available for atom type " +
             char4ToString(atom_type_query) + ".  Aliases for the first parameter set will be "
             "returned.", "LennardJonesAnalysis", "getLJAliases");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  return atom_type_aliases[lj_table_indices[0]];
}

//-------------------------------------------------------------------------------------------------
double LennardJonesAnalysis::getLJSigma(const int consensus_index) const {
  return sigma[consensus_index];
}

//-------------------------------------------------------------------------------------------------
double LennardJonesAnalysis::getLJ14Sigma(const int consensus_index) const {
  return sigma_14[consensus_index];
}

//-------------------------------------------------------------------------------------------------
double LennardJonesAnalysis::getLJEpsilon(const int consensus_index) const {
  return epsilon[consensus_index];
}

//-------------------------------------------------------------------------------------------------
double LennardJonesAnalysis::getLJ14Epsilon(const int consensus_index) const {
  return epsilon_14[consensus_index];
}

//-------------------------------------------------------------------------------------------------
double2 LennardJonesAnalysis::getLJParameters(const int consensus_index) const {
  validateConsensusIndexQuery(consensus_index);
  return { sigma[consensus_index], epsilon[consensus_index] };
}

//-------------------------------------------------------------------------------------------------
double2 LennardJonesAnalysis::getLJParameters(const char4 atom_type_query,
                                              const ExceptionResponse policy) const {
  confirmDistinctType(atom_type_query, policy, "getLJParameters");
  const uint query_key = char4ToUint(atom_type_query);
  const std::map<uint, std::vector<int>>::const_iterator it = atom_type_map.find(query_key);
  const size_t on_diag_idx = it->second[0] * (lj_type_count + 1);
  return { lja_coeff[on_diag_idx], ljb_coeff[on_diag_idx] };
}

//-------------------------------------------------------------------------------------------------
double2 LennardJonesAnalysis::getLJ14Parameters(const char4 atom_type_query,
                                                const ExceptionResponse policy) const {
  confirmDistinctType(atom_type_query, policy, "getLJ14Parameters");
  const uint query_key = char4ToUint(atom_type_query);
  const std::map<uint, std::vector<int>>::const_iterator it = atom_type_map.find(query_key);
  const size_t on_diag_idx = it->second[0] * (lj_type_count + 1);
  return { lja_14_coeff[on_diag_idx], ljb_14_coeff[on_diag_idx] };
}

//-------------------------------------------------------------------------------------------------
double3 LennardJonesAnalysis::getLJCoefficients(const int index_i, const int index_j) const {
  validateConsensusIndexQuery(index_i);
  validateConsensusIndexQuery(index_j);
  const size_t tab_idx = (index_j * lj_type_count) + index_i;
  return { lja_coeff[tab_idx], ljb_coeff[tab_idx], ljc_coeff[tab_idx] };
}

//-------------------------------------------------------------------------------------------------
double3 LennardJonesAnalysis::getLJ14Coefficients(const int index_i, const int index_j) const {
  validateConsensusIndexQuery(index_i);
  validateConsensusIndexQuery(index_j);
  const size_t tab_idx = (index_j * lj_type_count) + index_i;
  return { lja_14_coeff[tab_idx], ljb_14_coeff[tab_idx], ljc_14_coeff[tab_idx] };
}

//-------------------------------------------------------------------------------------------------
double3 LennardJonesAnalysis::getLJCoefficients(const char4 atype_i, const char4 atype_j,
                                                const ExceptionResponse policy) const {
  confirmDistinctType(atype_i, policy, "getLJCoefficients");
  confirmDistinctType(atype_j, policy, "getLJCoefficients");
  const uint query_key_i = char4ToUint(atype_i);
  const uint query_key_j = char4ToUint(atype_j);
  const std::map<uint, std::vector<int>>::const_iterator it_i = atom_type_map.find(query_key_i);
  const std::map<uint, std::vector<int>>::const_iterator it_j = atom_type_map.find(query_key_j);
  const std::vector<int>& lji_tab_idx = it_i->second;
  const std::vector<int>& ljj_tab_idx = it_j->second;
  const size_t tab_idx = (it_j->second[0] * lj_type_count) + it_i->second[0];
  return { lja_coeff[tab_idx], ljb_coeff[tab_idx], ljc_coeff[tab_idx] };
}

//-------------------------------------------------------------------------------------------------
double3 LennardJonesAnalysis::getLJ14Coefficients(const char4 atype_i, const char4 atype_j,
                                                const ExceptionResponse policy) const {
  confirmDistinctType(atype_i, policy, "getLJ14Coefficients");
  confirmDistinctType(atype_j, policy, "getLJ14Coefficients");
  const uint query_key_i = char4ToUint(atype_i);
  const uint query_key_j = char4ToUint(atype_j);
  const std::map<uint, std::vector<int>>::const_iterator it_i = atom_type_map.find(query_key_i);
  const std::map<uint, std::vector<int>>::const_iterator it_j = atom_type_map.find(query_key_j);
  const std::vector<int>& lji_tab_idx = it_i->second;
  const std::vector<int>& ljj_tab_idx = it_j->second;
  const size_t tab_idx = (it_j->second[0] * lj_type_count) + it_i->second[0];
  return { lja_14_coeff[tab_idx], ljb_14_coeff[tab_idx], ljc_14_coeff[tab_idx] };
}

//-------------------------------------------------------------------------------------------------
int LennardJonesAnalysis::getCorrespondence(int set_index, int type_index) const {
  if (set_index < 0 || set_index >= set_count) {
    rtErr("Set index " + std::to_string(set_index) + " is invalid for a collection of " +
          std::to_string(set_count) + " sets.", "LennardJonesAnalysis", "getSetCorrespondence");
  }
  if (type_index < 0 || type_index >= set_to_consensus_map[set_index].size()) {
    rtErr("Lennard-Jones tom type index " + std::to_string(type_index) + " is invalid for a "
          "collection of " + std::to_string(set_count) + " atom types in set " +
          std::to_string(set_index) + ".", "LennardJonesAnalysis", "getSetCorrespondence");
  }
  return set_to_consensus_map[set_index][type_index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& LennardJonesAnalysis::getSetCorrespondence(int set_index) const {
  if (set_index < 0 || set_index >= set_count) {
    rtErr("Set index " + std::to_string(set_index) + " is invalid for a collection of " +
          std::to_string(set_count) + " sets.", "LennardJonesAnalysis", "getSetCorrespondence");
  }
  return set_to_consensus_map[set_index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<int2>& LennardJonesAnalysis::getInputInstances(int consensus_index) const {
  if (consensus_index < 0 || consensus_index >= lj_type_count) {
    rtErr("Index " + std::to_string(consensus_index) + " is invalid for a consensus set of " +
          std::to_string(lj_type_count) + " Lennard-Jones atom types.", "LennardJonesAnalysis",
          "getSetCorrespondence");
  }
  return consensus_to_set_map[consensus_index];
}

//-------------------------------------------------------------------------------------------------
void LennardJonesAnalysis::addSet(const NonbondedKit<double> &nbk,
                                  const std::vector<std::vector<char4>> &othr_type_aliases) {

  // Check that there are enough alias lists for each atom type
  if (static_cast<int>(othr_type_aliases.size()) != nbk.n_lj_types) {
    rtErr("One atom type name (as part of a list that may include possible aliases) must be "
          "provided for each Lennard-Jones type.  There were " +
          std::to_string(othr_type_aliases.size()) + " lists of atom type names for " +
          std::to_string(nbk.n_lj_types) + " Lennard-Jones types.", "LennardJonesAnalysis",
          "addSet");
  }

  // Obtain sigma and epsilon parameters for the new set.  Note that, if the general non-bonded
  // Lennard-Jones parameters are zero, the 1:4 Lennard-Jones parameters will be assumed to be
  // negligible as well.
  std::vector<double> othr_sigma(nbk.n_lj_types), othr_14_sigma(nbk.n_lj_types);
  std::vector<double> othr_epsilon(nbk.n_lj_types), othr_14_epsilon(nbk.n_lj_types);
  for (int i = 0; i < nbk.n_lj_types; i++) {
    const size_t diag_idx = (i * nbk.n_lj_types) + i;
    if (fabs(nbk.ljb_coeff[diag_idx]) > constants::small) {
      othr_sigma[i] = sqrt(cbrt(nbk.lja_coeff[diag_idx] / nbk.ljb_coeff[diag_idx]));
      othr_epsilon[i] = 0.25 * nbk.ljb_coeff[diag_idx] / pow(othr_sigma[i], 6.0);
      othr_14_sigma[i] = sqrt(cbrt(nbk.lja_14_coeff[diag_idx] / nbk.ljb_14_coeff[diag_idx]));
      othr_14_epsilon[i] = 0.25 * nbk.ljb_14_coeff[diag_idx] / pow(othr_14_sigma[i], 6.0);
    }
    else {
      othr_sigma[i] = 0.0;
      othr_epsilon[i] = 0.0;
      othr_14_sigma[i] = 0.0;
      othr_14_epsilon[i] = 0.0;
    }
  }

  // Obtain a table of which atoms in the current table have pair-specific combining rules.  If
  // either the general non-bonded or the 1:4 non-bonded terms have pair-specific off-diagonal
  // elements then NBFix will be considered to be in effect.
  std::vector<bool> base_has_nbfix_terms =
    findPairSpecificParticipation(lja_coeff.data(), ljb_coeff.data(), lj_type_count, absolute_rule,
                                  prevalent_rule);
  const std::vector<bool> base_has_14nbf_terms = 
    findPairSpecificParticipation(lja_14_coeff.data(), ljb_14_coeff.data(), lj_type_count,
                                  absolute_rule, prevalent_rule);
  std::vector<bool> othr_has_nbfix_terms =
    findPairSpecificParticipation(nbk.lja_coeff, nbk.ljb_coeff, nbk.n_lj_types, absolute_rule,
                                  prevalent_rule);
  const std::vector<bool> othr_has_14nbf_terms =
    findPairSpecificParticipation(nbk.lja_14_coeff, nbk.ljb_14_coeff, nbk.n_lj_types,
                                  absolute_rule, prevalent_rule);
  for (int i = 0; i < nbk.n_lj_types; i++) {
    base_has_nbfix_terms[i] = (base_has_nbfix_terms[i] || base_has_14nbf_terms[i]);
    othr_has_nbfix_terms[i] = (othr_has_nbfix_terms[i] || othr_has_14nbf_terms[i]);
  }
  
  // Check for overlap in the atom types.  Make a list of all atom type names that are found in
  // both topologies, look up their indices, and check how they interact in each case.
  std::vector<int> ljt_correspondence(nbk.n_lj_types, -1);
  std::vector<bool> found_but_different(nbk.n_lj_types, false);
  for (int i = 0; i < nbk.n_lj_types; i++) {

    // The atom type of the incoming set will be considered unique if it has pair-specific
    // Lennard-Jones interactions with other types in its own set.
    if (othr_has_nbfix_terms[i]) {
      continue;
    }
    const int jlim = othr_type_aliases[i].size();
    int j = 0;
    bool not_found = true;
    while (not_found && j < jlim) {
      std::map<uint, std::vector<int>>::iterator it =
        atom_type_map.find(char4ToUint(othr_type_aliases[i][j]));
      if (it != atom_type_map.end()) {
        const std::vector<int>& current_table_indices = it->second;
        
        // An atom type representing the new topology's ith Lennard-Jones type is represented in
        // the existing tables.  But does it have the same Lennard-Jones properties in the
        // existing tables?  If not, it still needs to be considered a new type.
        const size_t n_defs = current_table_indices.size();
        for (size_t k = 0; k < n_defs; k++) {
          if (fabs(sigma[current_table_indices[k]] - othr_sigma[i]) <= 1.0e-6 &&
              fabs(sigma_14[current_table_indices[k]] - othr_14_sigma[i]) <= 1.0e-6 &&
              fabs(epsilon[current_table_indices[k]] - othr_epsilon[i]) <= 1.0e-6 &&
              fabs(epsilon_14[current_table_indices[k]] - othr_14_epsilon[i]) <= 1.0e-6 &&
              base_has_nbfix_terms[current_table_indices[k]] == false) {
            ljt_correspondence[i] = current_table_indices[k];
          }
          else {
            found_but_different[i] = true;
          }
        }
        not_found = false;
      }
      j++;
    }
    if (ljt_correspondence[i] == -1) {

      // If the atom type was not found, check that the Lennard-Jones parameters are indeed unique
      // and that no existing pair-specific parameters might confound the equivalence.  This
      // condition could be triggered after the ith element of found_but_different[] has been set
      // to TRUE, which will cause the map of atom type names to consensus parameter indices to be
      // expanded beyond a single index.
      for (int j = 0; j < lj_type_count; j++) {
        if (fabs(sigma[j] - othr_sigma[i]) <= 1.0e-6 &&
            fabs(sigma_14[j] - othr_14_sigma[i]) <= 1.0e-6 &&
            fabs(epsilon[j] - othr_epsilon[i]) <= 1.0e-6 &&
            fabs(epsilon_14[j] - othr_14_epsilon[i]) <= 1.0e-6 &&
            base_has_nbfix_terms[j] == false) {
          ljt_correspondence[i] = j;
        }
      }
    }
  }

  // The combination procedure should now isolate the new Lennard-Jones types and expand the
  // tables or atom type maps as appropriate.
  int expanded_cons_lj_types = lj_type_count;
  for (int i = 0; i < nbk.n_lj_types; i++) {
    expanded_cons_lj_types += (ljt_correspondence[i] < 0);
  }
  std::vector<double> tmp_lja_coeff(expanded_cons_lj_types * expanded_cons_lj_types, 0.0);
  std::vector<double> tmp_ljb_coeff(expanded_cons_lj_types * expanded_cons_lj_types, 0.0);
  std::vector<double> tmp_ljc_coeff(expanded_cons_lj_types * expanded_cons_lj_types, 0.0);
  std::vector<double> tmp_lja_14_coeff(expanded_cons_lj_types * expanded_cons_lj_types, 0.0);
  std::vector<double> tmp_ljb_14_coeff(expanded_cons_lj_types * expanded_cons_lj_types, 0.0);
  std::vector<double> tmp_ljc_14_coeff(expanded_cons_lj_types * expanded_cons_lj_types, 0.0);
  for (int i = 0; i < lj_type_count; i++) {
    for (int j = 0; j < lj_type_count; j++) {
      const int current_ij_idx = (i * lj_type_count) + j;
      const int current_ji_idx = (j * lj_type_count) + i;
      const int next_ij_idx = (i * expanded_cons_lj_types) + j;
      const int next_ji_idx = (j * expanded_cons_lj_types) + i;
      tmp_lja_coeff[next_ij_idx] = lja_coeff[current_ij_idx];
      tmp_ljb_coeff[next_ij_idx] = ljb_coeff[current_ij_idx];
      tmp_ljc_coeff[next_ij_idx] = ljc_coeff[current_ij_idx];
      tmp_lja_coeff[next_ji_idx] = lja_coeff[current_ji_idx];
      tmp_ljb_coeff[next_ji_idx] = ljb_coeff[current_ji_idx];
      tmp_ljc_coeff[next_ji_idx] = ljc_coeff[current_ji_idx];
      tmp_lja_14_coeff[next_ij_idx] = lja_14_coeff[current_ij_idx];
      tmp_ljb_14_coeff[next_ij_idx] = ljb_14_coeff[current_ij_idx];
      tmp_ljc_14_coeff[next_ij_idx] = ljc_14_coeff[current_ij_idx];
      tmp_lja_14_coeff[next_ji_idx] = lja_14_coeff[current_ji_idx];
      tmp_ljb_14_coeff[next_ji_idx] = ljb_14_coeff[current_ji_idx];
      tmp_ljc_14_coeff[next_ji_idx] = ljc_14_coeff[current_ji_idx];
    }
  }
  std::vector<int> tmp_set_to_consensus_map(nbk.n_lj_types);
  int next_lj_idx  = lj_type_count;
  consensus_to_set_map.resize(expanded_cons_lj_types);
  for (int i = lj_type_count; i < expanded_cons_lj_types; i++) {
    consensus_to_set_map[i] = std::vector<int2>();
  }
  for (int i = 0; i < nbk.n_lj_types; i++) {
    if (ljt_correspondence[i] == -1) {
      tmp_set_to_consensus_map[i] = next_lj_idx;
      consensus_to_set_map[next_lj_idx].push_back({ set_count, i });
      const size_t jlim = othr_type_aliases[i].size();
      if (found_but_different[i]) {
        for (size_t j = 0; j < jlim; j++) {
          atom_type_map[char4ToUint(othr_type_aliases[i][j])].push_back(next_lj_idx);
        }
      }
      atom_type_aliases.push_back(othr_type_aliases[i]);
      sigma.push_back(othr_sigma[i]);
      sigma_14.push_back(othr_14_sigma[i]);
      epsilon.push_back(othr_epsilon[i]);
      epsilon_14.push_back(othr_14_epsilon[i]);
      ljt_correspondence[i] = next_lj_idx;
      next_lj_idx++;
    }
    else {
      tmp_set_to_consensus_map[i] = ljt_correspondence[i];
      consensus_to_set_map[ljt_correspondence[i]].push_back({ set_count, i });
    }
  }
  set_to_consensus_map.push_back(tmp_set_to_consensus_map);

  // Extend the map.
  for (int i = 0; i < nbk.n_lj_types; i++) {
    const size_t n_alias = othr_type_aliases[i].size();
    for (size_t j = 0; j < n_alias; j++) {
      const uint ij_key = char4ToUint(othr_type_aliases[i][j]);
      const std::map<uint, std::vector<int>>::iterator it = atom_type_map.find(ij_key);
      if (it == atom_type_map.end()) {
        atom_type_map[ij_key] = std::vector<int>(1, ljt_correspondence[i]);
      }
      else {
        const std::vector<int>& at_idx = it->second;
        const size_t n_at = at_idx.size();
        bool found = false;
        for (size_t k = 0; k < n_at; k++) {
          found = (found || at_idx[k] == ljt_correspondence[i]);
        }
        if (found == false) {
          it->second.push_back(ljt_correspondence[i]);
        }
      }
    }
  }
  
  // Fill in the new regions of the consensus table using the appropriate combining rule, or by
  // direct copy of the contents of the incoming set.  Take the prevalent rule from the original
  // parameter set as the means for expanding the table to interactions between the original types
  // and any new types in the incoming set.
  const int prior_lj_type_count = lj_type_count;
  lj_type_count = next_lj_idx;
  atom_type_count = atom_type_map.size();
  lja_coeff = tmp_lja_coeff;
  ljb_coeff = tmp_ljb_coeff;
  ljc_coeff = tmp_ljc_coeff;
  lja_14_coeff = tmp_lja_14_coeff;
  ljb_14_coeff = tmp_ljb_14_coeff;
  ljc_14_coeff = tmp_ljc_14_coeff;
  for (int i = 0; i < prior_lj_type_count; i++) {
    for (int j = prior_lj_type_count; j < lj_type_count; j++) {
      const size_t ij_idx = (i * lj_type_count) + j;
      const size_t ji_idx = (j * lj_type_count) + i;
      const double eps_ij = sqrt(epsilon[i] * epsilon[j]);
      switch (prevalent_rule) {
      case VdwCombiningRule::GEOMETRIC:
        lja_coeff[ij_idx] = 4.0 * eps_ij * pow(sigma[i] * sigma[j], 6.0);
        ljb_coeff[ij_idx] = 4.0 * eps_ij * pow(sigma[i] * sigma[j], 3.0);
        lja_14_coeff[ij_idx] = 4.0 * eps_ij * pow(sigma_14[i] * sigma_14[j], 6.0);
        ljb_14_coeff[ij_idx] = 4.0 * eps_ij * pow(sigma_14[i] * sigma_14[j], 3.0);
        break;
      case VdwCombiningRule::LORENTZ_BERTHELOT:
      case VdwCombiningRule::NBFIX:

        // Assume that Lorentz-Berthelot is the combining rule for adding new types to a topology
        // where the prevalent rule is a free-for-all.
        lja_coeff[ij_idx] = 4.0 * eps_ij * pow(0.5 * (sigma[i] + sigma[j]), 12.0);
        ljb_coeff[ij_idx] = 4.0 * eps_ij * pow(0.5 * (sigma[i] + sigma[j]),  6.0);
        lja_14_coeff[ij_idx] = 4.0 * eps_ij * pow(0.5 * (sigma_14[i] + sigma_14[j]), 12.0);
        ljb_14_coeff[ij_idx] = 4.0 * eps_ij * pow(0.5 * (sigma_14[i] + sigma_14[j]),  6.0);
        break;
      }
      lja_coeff[ji_idx] = lja_coeff[ij_idx];
      ljb_coeff[ji_idx] = ljb_coeff[ij_idx];
      lja_14_coeff[ji_idx] = lja_14_coeff[ij_idx];
      ljb_14_coeff[ji_idx] = ljb_14_coeff[ij_idx];
    }
  }
  std::vector<int> add_row_col;
  add_row_col.reserve(lj_type_count - prior_lj_type_count);
  for (int i = 0; i < nbk.n_lj_types; i++) {
    if (ljt_correspondence[i] >= prior_lj_type_count) {
      add_row_col.push_back(i);
    }
  }
  for (int i = prior_lj_type_count; i < lj_type_count; i++) {
    const int iothr = add_row_col[i - prior_lj_type_count];
    for (int j = prior_lj_type_count; j < lj_type_count; j++) {
      const int jothr = add_row_col[j - prior_lj_type_count];
      lja_coeff[(i * lj_type_count) + j] = nbk.lja_coeff[(iothr * nbk.n_lj_types) + jothr];
      lja_coeff[(j * lj_type_count) + i] = nbk.lja_coeff[(jothr * nbk.n_lj_types) + iothr];
      ljb_coeff[(i * lj_type_count) + j] = nbk.ljb_coeff[(iothr * nbk.n_lj_types) + jothr];
      ljb_coeff[(j * lj_type_count) + i] = nbk.ljb_coeff[(jothr * nbk.n_lj_types) + iothr];
      lja_14_coeff[(i * lj_type_count) + j] = nbk.lja_14_coeff[(iothr * nbk.n_lj_types) + jothr];
      lja_14_coeff[(j * lj_type_count) + i] = nbk.lja_14_coeff[(jothr * nbk.n_lj_types) + iothr];
      ljb_14_coeff[(i * lj_type_count) + j] = nbk.ljb_14_coeff[(iothr * nbk.n_lj_types) + jothr];
      ljb_14_coeff[(j * lj_type_count) + i] = nbk.ljb_14_coeff[(jothr * nbk.n_lj_types) + iothr];
    }
  }  

  // Increment the number of covered sets
  set_count += 1;
}

//-------------------------------------------------------------------------------------------------
void LennardJonesAnalysis::validateConsensusIndexQuery(const int ljt_query) const {
  if (ljt_query >= lj_type_count || ljt_query < 0) {
    rtErr("Type index " + std::to_string(ljt_query) + " is invalid for a consensus set of " +
          std::to_string(lj_type_count) + " types.", "LennardJonesAnalysis",
          "validateConsensusIndexQuery");
  }
}

//-------------------------------------------------------------------------------------------------
void LennardJonesAnalysis::confirmDistinctType(const char4 atom_type_query,
                                               const ExceptionResponse policy,
                                               const char* caller) const {
  const uint query_key = char4ToUint(atom_type_query);
  const std::map<uint, std::vector<int>>::const_iterator it = atom_type_map.find(query_key);
  if (it == atom_type_map.end()) {
    rtErr("No atom type \"" + char4ToString(atom_type_query) + "\" could be located in the "
          "consensus tables.", "LennardJonesAnalysis", caller);
  }
  if (it->second.size() > 1) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("More than one parameter definition is available for atom type " +
            char4ToString(atom_type_query) + ".", "LennardJonesAnalysis", "getLJAliases");
    case ExceptionResponse::WARN:
      rtWarn("More than one parameter definition is available for atom type " +
             char4ToString(atom_type_query) + ".  Parameters for the first parameter set will be "
             "returned.", "LennardJonesAnalysis", "getLJParameters");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
std::string listAtomTypesAsString(const std::vector<char4> &atyp_list) {
  std::string result("[ ");
  const size_t n_alias = atyp_list.size();
  for (size_t i = 0; i < n_alias; i++) {
    result += char4ToString(atyp_list[i]) + " ";
  }
  result += "]";
  return result;
}

//-------------------------------------------------------------------------------------------------
int inferLennardJonesTypeCount(const int length_a, const int length_b, const char* caller) {
  if (length_a != length_b) {
    rtErr("Parameter matrices must have identical sizes (" + std::to_string(length_a) + " and " +
          std::to_string(length_b) + " provided).", caller);
  }
  const int n_lj_types = round(sqrt(length_a));
  if (n_lj_types * n_lj_types != static_cast<int>(length_a)) {
    rtErr("A number of atom types can only be inferred for square parameter matrices.", caller);
  }
  return n_lj_types;
}
  
} // namespace topology
} // namespace stormm
