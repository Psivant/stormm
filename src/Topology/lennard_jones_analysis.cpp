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
using parse::NumberFormat;
using parse::realToString;
  
//-------------------------------------------------------------------------------------------------
PairLJInteraction::PairLJInteraction(const char4 type_a_in, const char4 type_b_in,
                                     const double lja_in, const double ljb_in) :
    type_a{type_a_in}, type_b{type_b_in}, lja{lja_in}, ljb{ljb_in}
{}

//-------------------------------------------------------------------------------------------------
LennardJonesAnalysis::LennardJonesAnalysis(const AtomGraph *ag_in) :
    lj_type_count{0}, atom_type_count{0}, prevalent_rule{VdwCombiningRule::NBFIX},
    sigma{}, epsilon{}, lja_coeff{}, ljb_coeff{}, edits{}, ag_index_origins{},
    ag_index_origins_bounds{}, atom_type_aliases{}, atom_type_interactions{},
    consensus_index_map{}, topology_index_map{}, topology_index_map_bounds{},
    ag_pointers{1, const_cast<AtomGraph*>(ag_in)}
{
  // Reset the object to reflect holding zero topologies if the input was a null pointer
  if (ag_pointers[0] == nullptr) {
    ag_pointers.resize(0);
    return;
  }

  // Group the atom types involved in each Lennard-Jones interaction
  const NonbondedKit<double> nbk = ag_in->getDoublePrecisionNonbondedKit();
  lj_type_count = nbk.n_lj_types;
  atom_type_aliases = ag_in->getAtomTypeNameTable();
  atom_type_count = 0;
  for (int i = 0; i < nbk.n_lj_types; i++) {
    const int n_type_names = static_cast<int>(atom_type_aliases[i].size());
    atom_type_count += n_type_names;
    for (int j = 0; j < n_type_names; j++) {
      std::map<uint, int>::iterator it =
        atom_type_interactions.find(char4ToUint(atom_type_aliases[i][j]));
      if (it != atom_type_interactions.end()) {
        rtErr("Atom type " + char4ToString(atom_type_aliases[i][j]) + " controls multiple "
              "Lennard-Jones parameter sets in the consensus tables.", "LennardJonesAnalysis");
      }
      atom_type_interactions[char4ToUint(atom_type_aliases[i][j])] = i;
    }
  }
  
  // Compute the most prevalent Lennard-Jones rule.
  prevalent_rule = inferCombiningRule(ag_in, ExceptionResponse::SILENT, true);

  // Extract the sigma and epsilon parameters for the first topology.
  sigma   = ag_in->getLennardJonesSigma<double>();
  epsilon = ag_in->getLennardJonesEpsilon<double>();

  // Extract the Lennard-Jones A and B coefficients directly from the first topology.
  const size_t nlj_squared = nbk.n_lj_types * nbk.n_lj_types;
  lja_coeff.resize(nlj_squared);
  ljb_coeff.resize(nlj_squared);
  for (size_t i = 0; i < nlj_squared; i++) {
    lja_coeff[i] = nbk.lja_coeff[i];
    ljb_coeff[i] = nbk.ljb_coeff[i];
  }

  // Identify any edits that would be characteristic of NBFix.
  for (int i = 1; i < nbk.n_lj_types; i++) {
    for (int j = 0; j < i; j++) {
      const size_t ji_idx = (nbk.n_lj_types * i) + j;
      double ij_eps = epsilon[i] * epsilon[j];
      ij_eps = (ij_eps > constants::tiny) ? sqrt(ij_eps) : 0.0;
      double ij_sig;
      switch (prevalent_rule) {
      case VdwCombiningRule::GEOMETRIC:
        ij_sig = sigma[i] * sigma[j];
        ij_sig = (ij_sig > constants::tiny) ? sqrt(ij_sig) : 0.0;
        break;
      case VdwCombiningRule::LORENTZ_BERTHELOT:
        ij_sig = 0.5 * (sigma[i] + sigma[j]);
        break;
      case VdwCombiningRule::NBFIX:
        edits.push_back(PairLJInteraction(atom_type_aliases[i][0], atom_type_aliases[i][j],
                                          nbk.lja_coeff[ji_idx], nbk.ljb_coeff[ji_idx]));
        continue;
      }
      double ij_sig_six = ij_sig * ij_sig * ij_sig;
      ij_sig_six *= ij_sig_six;
      const double pred_ljb = 4.0 * ij_eps * ij_sig_six;
      const double pred_lja = pred_ljb * ij_sig_six;
      if (pred_ljb < constants::small) {
        if (fabs(pred_lja - nbk.lja_coeff[ji_idx]) > 1.0e-5 ||
            fabs(pred_ljb - nbk.ljb_coeff[ji_idx]) > 1.0e-5) {
          edits.push_back(PairLJInteraction(atom_type_aliases[i][0], atom_type_aliases[i][j],
                                            nbk.lja_coeff[ji_idx], nbk.ljb_coeff[ji_idx]));
        }
      }
      else {
        if (fabs((pred_lja - nbk.lja_coeff[ji_idx]) / nbk.lja_coeff[ji_idx]) > 1.0e-5 ||
            fabs((pred_ljb - nbk.ljb_coeff[ji_idx]) / nbk.ljb_coeff[ji_idx]) > 1.0e-5) {
          edits.push_back(PairLJInteraction(atom_type_aliases[i][0], atom_type_aliases[i][j],
                                            nbk.lja_coeff[ji_idx], nbk.ljb_coeff[ji_idx]));
        }
      }
    }
  }

  // Set the origins of each Lennard-Jones type index in the original topology.
  ag_index_origins.resize(nbk.n_lj_types);
  ag_index_origins_bounds.resize(nbk.n_lj_types + 1);
  ag_index_origins_bounds[0] = 0;
  for (int i = 0; i < nbk.n_lj_types; i++) {
    ag_index_origins[i] = { 0, i };
    ag_index_origins_bounds[i + 1] = i + 1;
  }
  
  // Initialize the consensus map
  consensus_index_map.resize(1);
  consensus_index_map[0].resize(nbk.n_lj_types);
  topology_index_map.resize(nbk.n_lj_types);
  topology_index_map_bounds.resize(nbk.n_lj_types + 1);
  topology_index_map_bounds[0] = 0;
  for (int i = 0; i < nbk.n_lj_types; i++) {
    consensus_index_map[0][i] = i;
    topology_index_map[i] = { 0, i };
    topology_index_map_bounds[i + 1] = i + 1;
  }
}

//-------------------------------------------------------------------------------------------------
LennardJonesAnalysis::LennardJonesAnalysis(const AtomGraph &ag_in) :
    LennardJonesAnalysis(ag_in.getSelfPointer())
{}

//-------------------------------------------------------------------------------------------------
int LennardJonesAnalysis::getLJTypeCount() const {
  return lj_type_count;
}

//-------------------------------------------------------------------------------------------------
int LennardJonesAnalysis::getAtomTypeCount() const {
  return atom_type_count;
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
const std::vector<char4>& LennardJonesAnalysis::getLJAliases(const int ag_query_index,
                                                             const int lj_type_index) const {
  if (ag_query_index < 0 || ag_query_index >= static_cast<int>(ag_pointers.size())) {
    rtErr("Topology index " + std::to_string(ag_query_index) + " is invalid for a list of " +
          std::to_string(ag_pointers.size()) + " referenced topologies.", "LennardJonesAnalysis",
          "getLJAliases");
  }
  if (lj_type_index < 0 || lj_type_index >= ag_pointers[ag_query_index]->getLJTypeCount()) {
    rtErr("Lennard-Jones interaction index " + std::to_string(lj_type_index) + " is invalid for "
          "topology " + ag_pointers[ag_query_index]->getFileName() + ", with " +
          std::to_string(ag_pointers[ag_query_index]->getLJTypeCount()) + ".",
          "LennardJonesAnalysis", "getLJAliases");
  }
  return atom_type_aliases[consensus_index_map[ag_query_index][lj_type_index]];
}

//-------------------------------------------------------------------------------------------------
const std::vector<char4>& LennardJonesAnalysis::getLJAliases(const AtomGraph *ag_query,
                                                             const int lj_type_index) const {
  const int ag_query_index = matchTopology(ag_query, ag_pointers);
  if (ag_query_index == static_cast<int>(ag_pointers.size())) {
    rtErr("No topology originating in " + ag_query->getFileName() + " was found in the list of "
          "referenced AtomGraph objects.", "LennardJonesAnalysis", "getLJAliases");
  }
  return getLJAliases(ag_query_index, lj_type_index);
}

//-------------------------------------------------------------------------------------------------
const std::vector<char4>& LennardJonesAnalysis::getLJAliases(const AtomGraph &ag_query,
                                                             const int lj_type_index) const {
  return getLJAliases(ag_query.getSelfPointer(), lj_type_index);
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
const std::vector<char4>& LennardJonesAnalysis::getLJAliases(const char4 atom_type_query) const {
  for (int i = 0; i < lj_type_count; i++) {
    const size_t n_names = atom_type_aliases[i].size();
    for (size_t j = 0; j < n_names; j++) {
      if (atom_type_aliases[i][j] == atom_type_query) {
        return atom_type_aliases[i];
      }
    }
  }
  rtErr("No atom type \"" + char4ToString(atom_type_query) + "\" could be located in the "
        "consensus tables.", "LennardJonesAnalysis", "getLJAliases");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double LennardJonesAnalysis::getLJSigma(const int consensus_index) const {
  return sigma[consensus_index];
}

//-------------------------------------------------------------------------------------------------
double LennardJonesAnalysis::getLJEpsilon(const int consensus_index) const {
  return epsilon[consensus_index];
}

//-------------------------------------------------------------------------------------------------
double2 LennardJonesAnalysis::getLJParameters(const int consensus_index) const {
  return { sigma[consensus_index], epsilon[consensus_index] };
}

//-------------------------------------------------------------------------------------------------
double2 LennardJonesAnalysis::getLJParameters(const char4 atom_type_query) const {
  for (int i = 0; i < lj_type_count; i++) {
    const size_t n_names = atom_type_aliases[i].size();
    for (size_t j = 0; j < n_names; j++) {
      if (atom_type_aliases[i][j] == atom_type_query) {
        return { sigma[i], epsilon[i] };
      }
    }
  }
  rtErr("No atom type \"" + char4ToString(atom_type_query) + "\" could be located in the "
        "consensus tables.", "LennardJonesAnalysis", "getLJSParameters");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void LennardJonesAnalysis::addTopology(const AtomGraph *ag_new_in) {
  const std::vector<std::vector<char4>> othr_atyp_aliases = ag_new_in->getAtomTypeNameTable();

  // Check for overlap in the atom types.  Make a list of all atom type names that are found in
  // both topologies, look up their indices, and check how they interact in each case.
  const std::vector<double> othr_sigma   = ag_new_in->getLennardJonesSigma<double>();
  const std::vector<double> othr_epsilon = ag_new_in->getLennardJonesEpsilon<double>();
  const NonbondedKit<double> nbk = ag_new_in->getDoublePrecisionNonbondedKit();
  std::vector<int> existing_ljt, othr_ljt;
  std::vector<bool> unique_types(nbk.n_lj_types, true);
  for (int i = 0; i < nbk.n_lj_types; i++) {
    const int jlim = othr_atyp_aliases.size();
    int lj_type_footprint = -1;
    for (int j = 0; j < jlim; j++) {
      std::map<uint, int>::iterator it =
        atom_type_interactions.find(char4ToUint(othr_atyp_aliases[i][j]));
      if (it != atom_type_interactions.end()) {
        unique_types[i] = false;

        // An atom type representing the new topology's ith Lennard-Jones type is represented in
        // the existing tables.  But does it have the same Lennard-Jones properties in the
        // existing tables?
        const int current_table_ljidx = it->second;
        othr_ljt.push_back(i);
        existing_ljt.push_back(current_table_ljidx);

        // Ensure that overlapping atom types are confined to specific Lennard-Jones types.
        if (lj_type_footprint == -1) {
          lj_type_footprint = current_table_ljidx;
        }
        else if (current_table_ljidx != lj_type_footprint) {
          std::string i_alias_list, cfoot_alias_list, ccurr_alias_list;
          for (size_t k = 0; k < othr_atyp_aliases[i].size(); k++) {
            i_alias_list += char4ToString(othr_atyp_aliases[i][k]) + " ";
          }
          for (size_t k = 0; k < atom_type_aliases[lj_type_footprint].size(); k++) {
            cfoot_alias_list += char4ToString(atom_type_aliases[lj_type_footprint][k]) + " ";
          }
          for (size_t k = 0; k < atom_type_aliases[current_table_ljidx].size(); k++) {
            ccurr_alias_list += char4ToString(atom_type_aliases[current_table_ljidx][k]) + " ";
          }
          rtErr(char4ToString(othr_atyp_aliases[i][j]) + " should map to the same consensus "
                "Lennard-Jones type as other atom types in its group [ " + i_alias_list +
                "].  Instead it maps to [ " + cfoot_alias_list + "] as well as [ " +
                ccurr_alias_list + "], and possibly the Lennard-Jones interactions of other "
                "groups in the consensus tables.", "LennardJonesAnalysis", "addTopology");
        }
      }
    }
  }
  
  // Check the Lennard-Jones A and B coefficients for all interaction types that are found in
  // common between the two systems.  If there are discrepancies, this is an error.
  const int n_common = othr_ljt.size();
  std::vector<double> common_lja(n_common * n_common), common_ljb(n_common * n_common);
  for (int i = 0; i < n_common; i++) {
    const int exst_ljt_i = existing_ljt[i];
    const int othr_ljt_i = othr_ljt[i];
    for (int j = 0; j < n_common; j++) {
      const int exst_ljt_j = existing_ljt[j];
      const int othr_ljt_j = othr_ljt[j];
      const double lja_exst = lja_coeff[(exst_ljt_j * lj_type_count) + exst_ljt_i];
      const double ljb_exst = ljb_coeff[(exst_ljt_j * lj_type_count) + exst_ljt_i];
      const double lja_othr = nbk.lja_coeff[(othr_ljt_j * nbk.n_lj_types) + othr_ljt_i];
      const double ljb_othr = nbk.ljb_coeff[(othr_ljt_j * nbk.n_lj_types) + othr_ljt_i];
      if (lja_exst < constants::tiny) {
        if (fabs(lja_exst - lja_othr) >= 1.0e-5) {
          std::string exst_type_list_i("[ "), exst_type_list_j("[ ");
          for (size_t k = 0; k < atom_type_aliases[exst_ljt_i].size(); k++) {
            const char4 exst_atyp = atom_type_aliases[exst_ljt_i][k];
            bool found = false;
            for (size_t m = 0; m < othr_atyp_aliases[othr_ljt_i].size(); m++) {
              found = (found || exst_atyp == othr_atyp_aliases[othr_ljt_i][m]);
            }
            if (found == false) {
              exst_type_list_i += char4ToString(exst_atyp) + " ";
            }
          }
          for (size_t k = 0; k < atom_type_aliases[exst_ljt_j].size(); k++) {
            const char4 exst_atyp = atom_type_aliases[exst_ljt_j][k];
            bool found = false;
            for (size_t m = 0; m < othr_atyp_aliases[othr_ljt_j].size(); m++) {
              found = (found || exst_atyp == othr_atyp_aliases[othr_ljt_j][m]);
            }
            if (found == false) {
              exst_type_list_j += char4ToString(exst_atyp) + " ";
            }
          }
          exst_type_list_i += "]";
          exst_type_list_j += "]";
          rtErr("Atom types " + exst_type_list_i + " and " + exst_type_list_j + "interact with "
                "different properties in the consensus tables than they do in the topology found "
                "in file " + ag_new_in->getFileName() + ".", "LennardJonesAnalysis",
                "addTopology");
        }
      }
    }
  }

  // Having survived the type comparisons, the combination procedure should now isolate the
  // new types and create new tables.

}

//-------------------------------------------------------------------------------------------------
void LennardJonesAnalysis::addTopology(const AtomGraph &ag_new_in) {
  addTopology(ag_new_in.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
void LennardJonesAnalysis::addTopology(const std::vector<AtomGraph*> &ag_list_in) {
  for (size_t i = 0; i < ag_list_in.size(); i++) {
    addTopology(ag_list_in[i]);
  }
}

//-------------------------------------------------------------------------------------------------
void LennardJonesAnalysis::addTopology(const AtomGraphSynthesis *poly_ag_in) {
  addTopology(poly_ag_in->getUniqueTopologies());
}

//-------------------------------------------------------------------------------------------------
void LennardJonesAnalysis::addTopology(const AtomGraphSynthesis &poly_ag_in) {
  addTopology(poly_ag_in.getUniqueTopologies());
}

} // namespace topology
} // namespace stormm
