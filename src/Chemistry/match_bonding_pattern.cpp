#include "copyright.h"
#include "Topology/atomgraph_abstracts.h"
#include "Math/summation.h"
#include "match_bonding_pattern.h"

namespace stormm {
namespace chemistry {

using stmath::sum;
using topology::ChemicalDetailsKit;
using topology::NonbondedKit;
  
//-------------------------------------------------------------------------------------------------
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
                         std::vector<int> *a_coverage, std::vector<int> *b_coverage) {

  // If the two atoms have the same index, return true.
  if (atom_a == atom_b) {
    return true;
  }
  
  // Check that both atoms are part of the same molecule
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  const int mol_ab = cdk.mol_home[atom_a];
  if (mol_ab != cdk.mol_home[atom_b]) {
    rtErr("Atom indices " + std::to_string(atom_a) + " and " + std::to_string(atom_b) +
          " reside in molecules " + std::to_string(cdk.mol_home[atom_a]) + " and " +
          std::to_string(cdk.mol_home[atom_b]) + ", respectively.  Bonding patterns can only "
          "match for atoms in the same molecule.", "matchBondingPatterns");
  }
  
  // Get the non-bonded abstract to track connectivity in the topology
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();

  // Set pointers for various storage arrays.  They are passed in to avoid needing to re-allocate
  // large amounts of memory for finding equivalencies in smaller and smaller subdivisions of the
  // system.
  int* a_idx_tptr = a_idx_tree->data();
  int* b_idx_tptr = b_idx_tree->data();
  int* a_zn_tptr = a_zn_tree->data();
  int* b_zn_tptr = b_zn_tree->data();
  double* a_fc_tptr = a_fc_tree->data();
  double* b_fc_tptr = b_fc_tree->data();
  double* a_fe_tptr = a_fe_tree->data();
  double* b_fe_tptr = b_fe_tree->data();
  ullint* a_ri_tptr = a_ri_tree->data();
  ullint* b_ri_tptr = b_ri_tree->data();
  ChiralOrientation* a_ch_tptr = a_ch_tree->data();
  ChiralOrientation* b_ch_tptr = b_ch_tree->data();
  int* a_coverage_ptr = a_coverage->data();
  int* b_coverage_ptr = b_coverage->data();
  
  // Initialize trees for each atom
  const int mol_natom = cdk.mol_limits[mol_ab + 1] - cdk.mol_limits[mol_ab];
  for (int k = cdk.mol_limits[mol_ab]; k < cdk.mol_limits[mol_ab + 1]; k++) {
    a_coverage_ptr[cdk.mol_contents[k]] = 0;
    b_coverage_ptr[cdk.mol_contents[k]] = 0;
  }
  a_coverage_ptr[atom_a] = 1;
  b_coverage_ptr[atom_b] = 1;
  a_idx_tptr[0] = atom_a;
  b_idx_tptr[0] = atom_b;
  bool add_to_trees = true;
  bool previous_layer_same = false;
  int last_layer_start = 0;
  int last_layer_end = 1;
  int n_atree_atoms = 1;
  int n_btree_atoms = 1;
  while (add_to_trees) {
    for (int i = last_layer_start; i < last_layer_end; i++) {
      const int a_tree_atom = a_idx_tptr[i];
      for (int j = nbk.nb12_bounds[a_tree_atom]; j < nbk.nb12_bounds[a_tree_atom + 1]; j++) {
        if (a_coverage_ptr[nbk.nb12x[j]] == 0) {
          const size_t next_atom = nbk.nb12x[j];
          a_idx_tptr[n_atree_atoms] = next_atom;
          a_zn_tptr[n_atree_atoms] = cdk.z_numbers[next_atom];
          a_fc_tptr[n_atree_atoms] = formal_charges[next_atom];
          a_fe_tptr[n_atree_atoms] = free_electrons[next_atom];
          a_ri_tptr[n_atree_atoms] = ring_inclusion[next_atom];
          a_ch_tptr[n_atree_atoms] = chiralities[next_atom];
          a_coverage_ptr[next_atom] = 1;
          n_atree_atoms++;
        }
      }
      const int b_tree_atom = b_idx_tptr[i];
      for (int j = nbk.nb12_bounds[b_tree_atom]; j < nbk.nb12_bounds[b_tree_atom + 1]; j++) {
        if (b_coverage_ptr[nbk.nb12x[j]] == 0) {
          const size_t next_atom = nbk.nb12x[j];
          b_idx_tptr[n_btree_atoms] = next_atom;
          b_zn_tptr[n_btree_atoms] = cdk.z_numbers[next_atom];
          b_fc_tptr[n_btree_atoms] = formal_charges[next_atom];
          b_fe_tptr[n_btree_atoms] = free_electrons[next_atom];
          b_ri_tptr[n_btree_atoms] = ring_inclusion[next_atom];
          b_ch_tptr[n_btree_atoms] = chiralities[next_atom];
          b_coverage_ptr[next_atom] = 1;
          n_btree_atoms++;
        }
      }
    }
    if (n_atree_atoms != n_btree_atoms) {
      return false;
    }
    
    // If the previous check passed, both trees are still the same size.  Update the layer bounds,
    // then check sums of the various Z numbers, free electrons, ring inclusions, and formal
    // charges.  Even if there are atoms in 60+ membered rings, the sums of unsigned long long ints
    // should still overflow in the same manner.
    last_layer_start = last_layer_end;
    last_layer_end = n_atree_atoms;
    int zn_sum_a = 0;
    int zn_sum_b = 0;
    int ch_sum_a = 0;
    int ch_sum_b = 0;
    double fc_sum_a = 0.0;
    double fc_sum_b = 0.0;
    double fe_sum_a = 0.0;
    double fe_sum_b = 0.0;
    for (int i = last_layer_start; i < last_layer_end; i++) {
      zn_sum_a += a_zn_tptr[i];
      zn_sum_b += b_zn_tptr[i];
    }
    if (zn_sum_a != zn_sum_b) {
      return false; 
    }
    for (int i = last_layer_start; i < last_layer_end; i++) {
      fc_sum_a += a_fc_tptr[i];
      fc_sum_b += b_fc_tptr[i];
      fe_sum_a += a_fe_tptr[i];
      fe_sum_b += b_fe_tptr[i];
      ch_sum_a += static_cast<int>(a_ch_tptr[i]);
      ch_sum_b += static_cast<int>(b_ch_tptr[i]);
    }
    if (fabs(fc_sum_a - fc_sum_b) > 1.0e-6 || fabs(fe_sum_a - fe_sum_b) > 1.0e-6 ||
        ch_sum_a != ch_sum_b) {
      return false; 
    }

    // Check to see whether the trees are adding atoms in the same order, and can thus be expected
    // to do so forever more.  While the sum of integer atom indices in a very, very large system
    // could, in principle, overflow, the results can still be expected to be equal for two
    // identical series of numbers, and int addition is faster than int => long long int conversion
    // with long long int addition.
    if (sum<int>(&a_idx_tptr[last_layer_start], last_layer_end - last_layer_start) ==
        sum<int>(&b_idx_tptr[last_layer_start], last_layer_end - last_layer_start)) {

      // Because all entries in each array will be unique, just check that each entry of one is
      // present in the other.
      bool all_covered = true;
      int search_start = last_layer_start;
      for (int i = last_layer_start; i < last_layer_end; i++) {
        const int a_tree_atom = a_idx_tptr[i];
        int j = search_start;
        bool found = false;
        while (j < last_layer_end && (! found)) {
          found = (a_tree_atom == b_idx_tptr[j]);
          search_start += (found && j == search_start);
          j++;
        }
        all_covered = (all_covered && found);
      }

      // If two successive layers are identical, then it can be inferred that subsequent layers
      // will only continue to add the same atoms with the same properties and connections.
      if (previous_layer_same) {
        if (all_covered) {
          return true;
        }
        else {
          previous_layer_same = false;
        }
      }
      else {
        previous_layer_same = all_covered;
      }
    }

    // Did the trees grow?
    add_to_trees = (last_layer_end > last_layer_start);
  }

  // If the trees stop growing, the entirety of the molecule is covered.  If nothing has indicated
  // a mismatch by that point, the bonding patterns must be identical.
  return true;
}

//-------------------------------------------------------------------------------------------------
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
                         std::vector<int> *a_coverage, std::vector<int> *b_coverage) {
  return matchBondingPattern(ag.getSelfPointer(), formal_charges, free_electrons, ring_inclusion,
                             chiralities, atom_a, atom_b, a_idx_tree, b_idx_tree, a_zn_tree,
                             b_zn_tree, a_fc_tree, b_fc_tree, a_fe_tree, b_fe_tree, a_ri_tree,
                             b_ri_tree, a_ch_tree, b_ch_tree, a_coverage, b_coverage);
}

//-------------------------------------------------------------------------------------------------
bool matchBondingPattern(const AtomGraph *ag, const std::vector<double> &formal_charges,
                         const std::vector<double> &free_electrons,
                         const std::vector<ullint> &ring_inclusion,
                         const std::vector<ChiralOrientation> &chiralities, const int atom_a,
                         const int atom_b) {
  const int natom = ag->getAtomCount();
  std::vector<int> a_idx_tree(natom), b_idx_tree(natom), a_zn_tree(natom), b_zn_tree(natom);
  std::vector<double> a_fc_tree(natom), b_fc_tree(natom), a_fe_tree(natom), b_fe_tree(natom);
  std::vector<ullint> a_ri_tree(natom), b_ri_tree(natom);
  std::vector<ChiralOrientation> a_ch_tree(natom), b_ch_tree(natom);
  std::vector<int> a_coverage(natom), b_coverage(natom);
  return matchBondingPattern(ag, formal_charges, free_electrons, ring_inclusion, chiralities,
                             atom_a, atom_b, &a_idx_tree, &b_idx_tree, &a_zn_tree, &b_zn_tree,
                             &a_fc_tree, &b_fc_tree, &a_fe_tree, &b_fe_tree, &a_ri_tree,
                             &b_ri_tree, &a_ch_tree, &b_ch_tree, &a_coverage, &b_coverage);
}

//-------------------------------------------------------------------------------------------------
bool matchBondingPattern(const AtomGraph &ag, const std::vector<double> &formal_charges,
                         const std::vector<double> &free_electrons,
                         const std::vector<ullint> &ring_inclusion,
                         const std::vector<ChiralOrientation> &chiralities, const int atom_a,
                         const int atom_b) {
  return matchBondingPattern(ag.getSelfPointer(), formal_charges, free_electrons, ring_inclusion,
                             chiralities, atom_a, atom_b);
}

} // namespace chemistry
} // namespace stormm
