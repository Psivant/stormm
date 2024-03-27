#include <cmath>
#include <cstdio>
#include <climits>
#include "copyright.h"
#include "atomgraph.h"

namespace stormm {
namespace topology {

//-------------------------------------------------------------------------------------------------
std::vector<char4> AtomGraph::matchOverflowAtomName(const std::string &query) const {
  const std::vector<int> match_indices =
    matchExtendedName(atom_overflow_names.data(), atom_overflow_names.size() / 4, query);
  std::vector<char4> result;
  const int nmatch = match_indices.size();
  for (int i = 0; i < nmatch; i++) {
    result.push_back(atom_overflow_names.readHost((4 * match_indices[i]) + 3));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<char4> AtomGraph::matchOverflowAtomType(const std::string &query) const {
  const std::vector<int> match_indices =
    matchExtendedName(atom_overflow_types.data(),
                      atom_overflow_types.size() / 4, query);
  std::vector<char4> result;
  const int nmatch = match_indices.size();
  for (int i = 0; i < nmatch; i++) {
    result.push_back(atom_overflow_names.readHost((4 * match_indices[i]) + 3));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<char4> AtomGraph::matchOverflowResidueName(const std::string &query) const {
  const std::vector<int> match_indices =
    matchExtendedName(residue_overflow_names.data(), residue_overflow_names.size() / 4, query);
  std::vector<char4> result;
  const int nmatch = match_indices.size();
  for (int i = 0; i < nmatch; i++) {
    result.push_back(residue_overflow_names.readHost((4 * match_indices[i]) + 3));
  }
  return result;
}

} // namespace topology
} // namespace stormm
