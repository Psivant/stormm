#include <string>
#include "copyright.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "topology_util.h"

namespace stormm {
namespace topology {

using stmath::findBin;
using stmath::prefixSumInPlace;
using stmath::PrefixSumType;
using parse::char4ToString;
  
//-------------------------------------------------------------------------------------------------
std::string writeAtomList(const std::vector<int> &atom_list, const ChemicalDetailsKit &cdk) {
  std::string result;
  const size_t natom = atom_list.size(); 
  for (size_t i = 0; i < natom; i++) {
    const int atom_idx = atom_list[i];
    const int res_idx = findBin(cdk.res_limits, atom_idx, cdk.nres);
    result += std::to_string(cdk.atom_numbers[atom_idx]) + " " +
              char4ToString(cdk.atom_names[atom_idx]) + " " +
              char4ToString(cdk.res_names[res_idx]) + " " +
              std::to_string(cdk.res_numbers[atom_idx]);
    if (i < natom - 1LLU) {
      result += ", ";
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
bool isBonded(const NonbondedKit<double> &nbk, const int atom_i, const int atom_j) {
  const int hlim = nbk.nb12_bounds[atom_i + 1];
  bool result = false;
  for (int i = nbk.nb12_bounds[atom_i]; i < hlim; i++) {
    result = (result || (atom_j == nbk.nb12x[i]));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
bool isBonded(const AtomGraph &ag, const int atom_i, const int atom_j) {
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  return isBonded(nbk, atom_i, atom_j);
}

//-------------------------------------------------------------------------------------------------
bool isBonded(const AtomGraph *ag, const int atom_i, const int atom_j) {
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  return isBonded(nbk, atom_i, atom_j);
}

//-------------------------------------------------------------------------------------------------
int matchTopology(const AtomGraph *query_ag, const std::vector<AtomGraph*> &repo) {
  const int ntop = repo.size();
  for (int i = 0; i < ntop; i++) {
    if (query_ag == repo[i]) {
      return i;
    }
  }
  const int natom = query_ag->getAtomCount();
  const std::string fname = query_ag->getFileName();
  for (int i = 0; i < ntop; i++) {
    if (natom == repo[i]->getAtomCount() && fname == repo[i]->getFileName()) {
      return i;
    }
  }
  return ntop;
}

//-------------------------------------------------------------------------------------------------
int matchTopology(const AtomGraph &query_ag, const std::vector<AtomGraph*> &repo) {
  return matchTopology(query_ag.getSelfPointer(), repo);
}

//-------------------------------------------------------------------------------------------------
int matchTopology(const AtomGraph *query_ag, const std::vector<AtomGraph> &repo) {
  const int ntop = repo.size();
  for (int i = 0; i < ntop; i++) {
    if (query_ag == &repo[i]) {
      return i;
    }
  }
  const int natom = query_ag->getAtomCount();
  const std::string fname = query_ag->getFileName();
  for (int i = 0; i < ntop; i++) {
    if (natom == repo[i].getAtomCount() && fname == repo[i].getFileName()) {
      return i;
    }
  }
  return ntop;
}

//-------------------------------------------------------------------------------------------------
int matchTopology(const AtomGraph &query_ag, const std::vector<AtomGraph> &repo) {
  return matchTopology(query_ag.getSelfPointer(), repo);
}

//-------------------------------------------------------------------------------------------------
bool hasVdwProperties(const AtomGraph *ag, const int atom_index, const VdwCombiningRule lj_rule) {
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  return hasVdwProperties<double>(nbk, atom_index, lj_rule);
}

//-------------------------------------------------------------------------------------------------
bool hasVdwProperties(const AtomGraph &ag, const int atom_index, const VdwCombiningRule lj_rule) {
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  return hasVdwProperties<double>(nbk, atom_index, lj_rule);
}

} // namespace topology
} // namespace stormm
