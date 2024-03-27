// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace topology {

//-------------------------------------------------------------------------------------------------
template <typename T>
void ComboGraphLJModel::addCombination(const T* lja_in, const T* ljb_in, const int lj_type_count,
                                       const char4* lj_type_names,
                                       const std::vector<PairLJInteraction> &edits) {
  secondary_topology_rules.push_back(inferCombiningRule(lja_in, ljb_in, lj_type_count));

  // Compute the combination of the original (primary) topology's Lennard-Jones parameters with the
  // new topology's Lennard-Jones parameters, according to the default rule.  The default rule for
  // mixing the two topologies could, in principle, be different than the combining rules apparent
  // in either topology itself.  Each atom type of the new topology controls a column of the
  // resulting matrix, each atom type of the primary topology controls a row.
  const NonbondedKit<double> primary_nbk = primary_ag_pointer->getDoublePrecisionNonbondedKit();
  std::vector<double> combi_lja(primary_nbk.n_lj_types * lj_type_count);
  std::vector<double> combi_ljb(primary_nbk.n_lj_types * lj_type_count);
  const std::vector<double> primary_sig = primary_ag_pointer->getLennardJonesSigma<double>();
  const std::vector<double> primary_eps = primary_ag_pointer->getLennardJonesEpsilon<double>();
  for (int j = 0; j < lj_type_count; j++) {
    const double secondary_sig  = sqrt(cbrt(lja_in[j] / ljb_in[j]));
    const double nsig_three = secondary_sig * secondary_sig * secondary_sig;
    const double secondary_eps = 0.25 * ljb_in[j] / (nsig_three * nsig_three);
    for (int i = 0; i < primary_nbk.n_lj_types; i++) {
      const double ij_eps = sqrt(secondary_eps * primary_eps[i]);
      double ij_sig;
      switch (default_rule) {
      case VdwCombiningRule::LORENTZ_BERTHELOT:
        ij_sig = 0.5 * (secondary_sig + primary_sig[i]);
        break;
      case VdwCombiningRule::GEOMETRIC:
        ij_sig = sqrt(secondary_sig * primary_sig[i]);        
        break;
      case VdwCombiningRule::NBFIX:

        // This invalid case is trapped in the constructor.
        break;
      }
      double ij_sig_six = (ij_sig * ij_sig * ij_sig);
      ij_sig_six *= ij_sig_six;
      const size_t ij_idx = (j * primary_nbk.n_lj_types) + i;
      combi_ljb[ij_idx] = 4.0 * ij_eps * ij_sig_six;
      combi_lja[ij_idx] = combi_ljb[ij_idx] * ij_sig_six;
    }      
  }

  // Apply the Lennard-Jones edits
  VdwCombiningRule combi_rule = default_rule;
  if (edits.size() > 0) {
    const std::vector<std::vector<char4>> primary_type_names =
      primary_ag_pointer->getAtomTypeNameTable();
    const size_t nedit = edits.size();
    for (size_t i = 0; i < nedit; i++) {

      // The atom types are not required to be found.
      for (int j = 0; j < primary_nbk.n_lj_types; j++) {
        const int primary_type_aliases = primary_type_names[j].size();
        for (int k = 0; k < lj_type_count; k++) {
          bool type_a_found = false;
          bool type_b_found = false;
          for (int m = 0; m < primary_type_aliases; m++) {
            type_a_found = (type_a_found || (edits[i].type_a == primary_type_names[j][m]));
            type_b_found = (type_b_found || (edits[i].type_b == primary_type_names[j][m]));
          }
          if ((type_a_found && edits[i].type_b == lj_type_names[k]) ||
              (type_b_found && edits[i].type_a == lj_type_names[k])) {

            // The types to be edited are matched.  Assign the new parameters and note that the
            // combining rule must now be assumed to be "NBFIX".
            const size_t jk_idx = (k * primary_nbk.n_lj_types) + j;
            combi_lja[jk_idx] = edits[i].lja;
            combi_ljb[jk_idx] = edits[i].ljb;
            combi_rule = VdwCombiningRule::NBFIX;
          }
        }
      }
    }
  }
  
  // Push the new set sizes, parameter tables, and edits to the object's growing arrays.
  secondary_atom_type_counts.push_back(lj_type_count);
  set_rules.push_back(combi_rule);
  set_lja.push_back(combi_lja);
  set_ljb.push_back(combi_ljb);
  set_edits.push_back(edits);

  // Set the source of the new tables as the null pointer.  If there was an actual topology
  // providing the input Lennard-Jones parameters, a pointer to it will be supplied by one of the
  // calling overloads and replace the nullptr.
  secondary_ag_pointers.push_back(nullptr);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void ComboGraphLJModel::addCombination(const std::vector<T> &lja_in, const std::vector<T> &ljb_in,
                                       const std::vector<char4> &lj_type_names,
                                       const std::vector<PairLJInteraction> &edits) {
  const int lj_type_count = inferLennardJonesTypeCount(lja_in.size(), ljb_in.size(),
                                                       "addCombination");

  // Trap cases where edits are provided but the names of atom types are not.
  if (edits.size() > 0 && static_cast<int>(lj_type_names.size()) != lj_type_count) {
    rtErr("Indicating Lennard-Jones edits requires that the names of atom types be provided "
          "for matching purposes.", "ComboGraphLJModel", "AddCombination");
  }
  addCombination(lja_in.data(), ljb_in.data(), lj_type_count, lj_type_names.data(), edits);
}

} // namespace topology
} // namespace stormm
