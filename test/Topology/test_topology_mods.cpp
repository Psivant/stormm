#include "copyright.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/ForceField/forcefield_element.h"
#include "../../src/ForceField/forcefield_enumerators.h"
#include "../../src/MolecularMechanics/minimization.h"
#include "../../src/MolecularMechanics/mm_evaluation.h"
#include "../../src/Namelists/nml_minimize.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/static_exclusionmask.h"
#include "../../src/Potential/local_exclusionmask.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/structure_ops.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Topology/atomgraph_refinement.h"
#include "../../src/Topology/atomgraph_stage.h"
#include "../../src/Topology/lennard_jones_analysis.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::mm;
using namespace stormm::modeling;
using namespace stormm::namelist;
using namespace stormm::parse;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::topology;
using namespace stormm::trajectory;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Confirm the match with a bond in another topology.  The atom names, size of the home molecule,
// and relative displacment of the two atom indices in the topological order are all considered.
// The function returns the bond term index if a match could be found, or -1 if no match was found.
//
// Arguments:
//   cdk:              Chemical details of the system in which to attempt to find the bond,
//                     containing atom names and molecule sizes
//   vk:               Valence interaction details of the system in which to attempt to find the
//                     bond
//   target_mol_size:  The target size of molecule in which the bond is to be found
//   target_topl_idx:  The target relative difference between the two atom indices, assuming that
//                     the corresponding molecules are laid out in the same order in each topology
//   target_bi_name:   The name of the bond's i atom
//   target_bj_name:   The name of the bond's j atom
//   bonds_covered:    Mask to indicate which bonds in the input topology have already been covered
//                     and how many times they have been covered
//   max_coverage:     The maximum number of times a given bond can be covered
//-------------------------------------------------------------------------------------------------
int confirmBondMatch(const ChemicalDetailsKit &cdk, const ValenceKit<double> &vk,
                     const int target_mol_size, const int target_topl_idx,
                     const char4 target_bi_name, const char4 target_bj_name,
                     const std::vector<int> &bonds_covered, const int max_coverage) {
  for (int pos = 0; pos < vk.nbond; pos++) {
    const int rel_topl_idx = abs(vk.bond_i_atoms[pos] - vk.bond_j_atoms[pos]);
    if (rel_topl_idx == target_topl_idx) {
      const int mol_idx = cdk.mol_home[vk.bond_i_atoms[pos]];
      const int mol_size = cdk.mol_limits[mol_idx + 1] - cdk.mol_limits[mol_idx];
      if (mol_size == target_mol_size) {
        const char4 bi_name = cdk.atom_names[vk.bond_i_atoms[pos]];
        const char4 bj_name = cdk.atom_names[vk.bond_j_atoms[pos]];
        if (((bi_name == target_bi_name && bj_name == target_bj_name) ||
             (bi_name == target_bj_name && bj_name == target_bi_name)) &&
            bonds_covered[pos] < max_coverage) {
          return pos;
        }
      }
    }
  }
  return -1;
}

//-------------------------------------------------------------------------------------------------
// Confirm the match with an angle in another topology.  Descriptions of input arguments follow
// from confirmBondMatch(), above, in addition to:
//
// Arguments:
//   target_ij_idx:   Target difference between the topological indices of atoms i and j in the
//                    angle term
//   target_jk_idx:   Target difference between the topological indices of atoms j and k in the
//                    angle term
//   target_bk_name:  Expected name of the third atom in the angle term
//-------------------------------------------------------------------------------------------------
int confirmAngleMatch(const ChemicalDetailsKit &cdk, const ValenceKit<double> &vk,
                      const int target_mol_size, const int target_ij_idx, const int target_jk_idx,
                      const char4 target_bi_name, const char4 target_bj_name,
                      const char4 target_bk_name, const std::vector<int> &angls_covered,
                      const int max_coverage) {
  for (int pos = 0; pos < vk.nangl; pos++) {
    const int ij_topl_idx = abs(vk.angl_i_atoms[pos] - vk.angl_j_atoms[pos]);
    const int jk_topl_idx = abs(vk.angl_j_atoms[pos] - vk.angl_k_atoms[pos]);
    if (ij_topl_idx == target_ij_idx && jk_topl_idx == target_jk_idx) {
      const int mol_idx = cdk.mol_home[vk.angl_i_atoms[pos]];
      const int mol_size = cdk.mol_limits[mol_idx + 1] - cdk.mol_limits[mol_idx];
      if (mol_size == target_mol_size) {
        const char4 bi_name = cdk.atom_names[vk.angl_i_atoms[pos]];
        const char4 bj_name = cdk.atom_names[vk.angl_j_atoms[pos]];
        const char4 bk_name = cdk.atom_names[vk.angl_k_atoms[pos]];
        if (bj_name == target_bj_name &&
            ((bi_name == target_bi_name && bk_name == target_bk_name) ||
             (bi_name == target_bk_name && bk_name == target_bi_name)) &&
            angls_covered[pos] < max_coverage) {
          return pos;
        }
      }
    }
  }
  return -1;
}

//-------------------------------------------------------------------------------------------------
// Confirm the match with an angle in another topology.  Descriptions of input arguments follow
// from confirmBondMatch() and confirmAngleMatch(), above, in addition to:
//
// Arguments:
//   target_kl_idx:   Target difference between the topological indices of atoms k and k in the
//                    dihedral term
//   target_bl_name:  Expected name of the fourth atom in the dihedral term
//-------------------------------------------------------------------------------------------------
int confirmDihedralMatch(const ChemicalDetailsKit &cdk, const ValenceKit<double> &vk,
                         const int target_mol_size, const int target_ij_idx,
                         const int target_jk_idx, const int target_kl_idx,
                         const char4 target_bi_name, const char4 target_bj_name,
                         const char4 target_bk_name, const char4 target_bl_name,
                         const std::vector<int> &dihes_covered, const int max_coverage) {
  for (int pos = 0; pos < vk.ndihe; pos++) {
    const int ij_topl_idx = abs(vk.dihe_i_atoms[pos] - vk.dihe_j_atoms[pos]);
    const int jk_topl_idx = abs(vk.dihe_j_atoms[pos] - vk.dihe_k_atoms[pos]);
    const int kl_topl_idx = abs(vk.dihe_k_atoms[pos] - vk.dihe_l_atoms[pos]);
    if (ij_topl_idx == target_ij_idx && jk_topl_idx == target_jk_idx &&
        kl_topl_idx == target_kl_idx) {
      const int mol_idx = cdk.mol_home[vk.dihe_i_atoms[pos]];
      const int mol_size = cdk.mol_limits[mol_idx + 1] - cdk.mol_limits[mol_idx];
      if (mol_size == target_mol_size) {
        const char4 bi_name = cdk.atom_names[vk.dihe_i_atoms[pos]];
        const char4 bj_name = cdk.atom_names[vk.dihe_j_atoms[pos]];
        const char4 bk_name = cdk.atom_names[vk.dihe_k_atoms[pos]];
        const char4 bl_name = cdk.atom_names[vk.dihe_l_atoms[pos]];
        if (((bi_name == target_bi_name && bj_name == target_bj_name &&
              bk_name == target_bk_name && bl_name == target_bl_name) ||
             (bi_name == target_bl_name && bj_name == target_bk_name &&
              bk_name == target_bj_name && bl_name == target_bi_name)) &&
            dihes_covered[pos] < max_coverage) {
          return pos;
        }
      }
    }
  }
  return -1;
}

//-------------------------------------------------------------------------------------------------
// Get the list of attenuated 1:4 pairs along with scaling constants for electrostatic and van-der
// Waals interactions.  An array of tuples is returned.  For each pair, the atom index of the first
// atom is returned in the "x" member of the tuple, the atom index of the second in the "y" member,
// the electrostatic scaling factor in the "z" member, and the van-der Waals scaling factor in the
// "w" member.
//
// Arguments:
//   vk:        Valence parameters for the system in question
//   nbk:       Non-bonded parameters for the system in question
//   desc:      Description of the system from which the list of interactions originated
//   do_tests:  Indication of whether testing is possible
//-------------------------------------------------------------------------------------------------
std::vector<double4_16a> getAttenuated14PairList(const ValenceKit<double> &vk,
                                             const NonbondedKit<double> &nbk,
                                             const std::string &desc,
                                             const TestPriority do_tests) {
  std::vector<double4_16a> result;
  for (int pos = 0; pos < vk.ndihe; pos++) {
    const int attn_idx = vk.dihe14_param_idx[pos];
    if (attn_idx == 0) {
      continue;
    }
    
    result.push_back({ static_cast<double>(vk.dihe_i_atoms[pos]),
                       static_cast<double>(vk.dihe_l_atoms[pos]),
                       vk.attn14_elec[attn_idx], vk.attn14_vdw[attn_idx] });
  }

  // Check the result for duplicates
  std::vector<double4_16a> lcopy = result;
  std::sort(lcopy.begin(), lcopy.end(), [](double4_16a a, double4_16a b) { return (a.x < b.x); });
  const int n = lcopy.size();
  std::vector<int2> duplicates;
  for (int i = 0; i < n; i++) {
    int j = i + 1;
    while (j < n && lcopy[i].x == lcopy[j].x) {
      if (fabs(lcopy[i].y - lcopy[j].y) < 1.0e-4 && fabs(lcopy[i].z - lcopy[j].z) < 1.0e-4 &&
          fabs(lcopy[i].w - lcopy[j].w) < 1.0e-4) {
        duplicates.push_back({ i, j });
      }
      j++;
    }
  }
  check(duplicates.size(), RelationalOperator::EQUAL, 0, "Duplicate entries were found in a list "
        "of attenuated 1-4 interactions pertaining to the " + desc + " system.", do_tests);
  return result;
}

//-------------------------------------------------------------------------------------------------
// Match the details of an attenuated 1:4 interaction to one member of a list.
//
// Arguments:
//   item:
//   list:
//   offset:
//   coverage:
//-------------------------------------------------------------------------------------------------
bool matchAttenuationToList(const double4_16a item, const std::vector<double4_16a> &list,
                            const double offset, std::vector<bool> *coverage) {
  const int llen = list.size();
  for (int i = 0; i < llen; i++) {
    if (((fabs(item.x - (list[i].x - offset)) < 1.0e-4 &&
          fabs(item.y - (list[i].y - offset)) < 1.0e-4) ||
         (fabs(item.y - (list[i].x - offset)) < 1.0e-4 &&
          fabs(item.x - (list[i].y - offset)) < 1.0e-4)) &&
        fabs(item.z - list[i].z) < 1.0e-6 && fabs(item.w - list[i].w) < 1.0e-6) {
      coverage->at(i) = true;
      return true;
    }
  }
  return false;
}

//-------------------------------------------------------------------------------------------------
// Construct a message to reveal the indices and names of some examples of missing pairs out of a
// list of attenuations.
//
// Arguments:
//   pairs:  The list of missing pairs, with the atom indices represented as real numbers in the
//           "x" and "y" members of the tuple.  Electrostatic and van-der Waals scaling factors for
//           the attenuated interaction are listed in the "z" and "w" members, respectively.
//   cdk:    Contains atom names of the system
//-------------------------------------------------------------------------------------------------
std::string listMissingPairs(const std::vector<double4_16a> &pairs,
                             const ChemicalDetailsKit &cdk) {
  std::string result;
  if (pairs.size() > 0) {
    const int nrep = std::min(static_cast<int>(pairs.size()), 8);
    result += "Examples of missing pairs include:\n";
    for (int i = 0; i < nrep; i++) {
      const int i_atom = round(pairs[i].x);
      const int j_atom = round(pairs[i].y);
      result += "  - Indices { " + std::to_string(i_atom) + ", " + std::to_string(j_atom) +
                " } " + char4ToString(cdk.atom_names[i_atom]) + " " +
                char4ToString(cdk.atom_names[j_atom]) + "\n"; 
    }
    result += "\n";
  }
  return result;
}

std::string listMissingPairs(const std::vector<int2> &pairs, const ChemicalDetailsKit &cdk) {
  std::string result;
  if (pairs.size() > 0) {
    const int nrep = std::min(static_cast<int>(pairs.size()), 8);
    result += "Examples of missing pairs include:\n";
    for (int i = 0; i < nrep; i++) {
      result += "  - Indices { " + std::to_string(pairs[i].x) + ", " + std::to_string(pairs[i].y) +
                " } " + char4ToString(cdk.atom_names[pairs[i].x]) + " " +
                char4ToString(cdk.atom_names[pairs[i].y]) + "\n"; 
    }
    result += "\n";
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Check that, for a simple combination of two topologies (each containing one molecule and adding
// them one after the other), all pairs of attentuated 1:4 interactions are accounted for and that
// the interaction parameters are consistent.
//
// Arguments:
//   list_a:      The first topology's list of 1-4 interactions
//   list_b:      The second topology's list of 1-4 interactions
//   list_comp:   The combined list of 1-4 interactions
//   sys_a_file:  Original file name of the first system (for error reporting)
//   sys_b_file:  Original file name of the second system
//   cdk_a:       Details of the first topology, including atom and residue names
//   cdk_b:       Details of the second topology, including atom and residue names
//   cdk_comp:    Details of the combined topology
//   do_tests:    Indication of whether testing is possible
//-------------------------------------------------------------------------------------------------
void checkAttenuated14PairLists(const std::vector<double4_16a> &list_a,
                                const std::vector<double4_16a> &list_b,
                                const std::vector<double4_16a> &list_comp,
                                const std::string &sys_a_file, const std::string &sys_b_file,
                                const ChemicalDetailsKit &cdk_a, const ChemicalDetailsKit &cdk_b,
                                const ChemicalDetailsKit &cdk_comp, const TestPriority do_tests) {
  const int na = list_a.size();
  const int nb = list_b.size();
  const int nc = list_comp.size();
  std::vector<bool> sysa_in_complex(na);
  std::vector<bool> sysb_in_complex(nb);
  std::vector<bool> complex_coverage(nc, false);
  int na_not_represented = 0;
  std::vector<double4_16a> unrepresented_a_pairs, unrepresented_b_pairs;
  for (int pos = 0; pos < na; pos++) {
    sysa_in_complex[pos] = matchAttenuationToList(list_a[pos], list_comp, 0, &complex_coverage);
    if (sysa_in_complex[pos] == false) {
      unrepresented_a_pairs.push_back(list_a[pos]);
      na_not_represented++;
    }
  }
  int nb_not_represented = 0;
  for (int pos = 0; pos < nb; pos++) {
    sysb_in_complex[pos] = matchAttenuationToList(list_b[pos], list_comp, cdk_a.natom,
                                                  &complex_coverage);
    if (sysb_in_complex[pos] == false) {
      unrepresented_b_pairs.push_back(list_b[pos]);
      nb_not_represented++;
    }
  }
  const std::string err_a_msg = listMissingPairs(unrepresented_a_pairs, cdk_a);
  check(na_not_represented, RelationalOperator::EQUAL, 0, "Some attenuated non-bonded pair "
        "interactions found in the first topology (" + sys_a_file + ", " +
        std::to_string(cdk_a.natom) + " atoms) were not represented in the complex (" +
        std::to_string(cdk_comp.natom) + " atoms).  " + err_a_msg, do_tests);
  const std::string err_b_msg = listMissingPairs(unrepresented_b_pairs, cdk_b);
  check(nb_not_represented, RelationalOperator::EQUAL, 0, "Some attenuated non-bonded pair "
        "interactions found in the second topology (" + sys_b_file + ", " +
        std::to_string(cdk_b.natom) + " atoms) were not represented in the complex (" +
        std::to_string(cdk_comp.natom) + " atoms).  " + err_b_msg, do_tests);
}

//-------------------------------------------------------------------------------------------------
// Test the combination of two topologies.  This testing apparatus contains implicit assumptions
// that the topologies contain only one molecule each, and that the first topology, larger than the
// first, will have all of its instances come first in the product.
//
// Arguments:
//   ag_a:      The first topology to combine
//   ag_b:      The second topology to roll into the product
//   ps_a:      Coordinates of the first system
//   ps_b:      Coordinates of the second system
//   n_a:       Optional number of copies of the first topology to include
//   n_b:       Optional number of copies of the second topology to include
//   dr_a:      Displacement vector to be applied to copies of systems described by the first
//              topology.  The coordinates of the first instance will remain "as is" while the
//              coordinates of any subsequent additions will be translated by the multiples of the
//              components of this tuple.
//   dr_b:      Displacement vector applicable to copies of systems described by the second
//              topology.
//   do_tests:  Flag to indicate that testing is feasible
//-------------------------------------------------------------------------------------------------
void testTopologyFusion(const AtomGraph &ag_a, const AtomGraph &ag_b, const PhaseSpace &ps_a,
                        const PhaseSpace &ps_b, const int n_a = 1, const int n_b = 1,
                        const double3 dr_a = { 1000.0, 1000.0, 1000.0 },
                        const double3 dr_b = { 1000.0, 1000.0, -1000.0 },
                        const TestPriority do_tests = TestPriority::CRITICAL) {
  AtomGraph ag_ab(ag_a, n_a, ag_b, n_b);
  const ValenceKit<double> agb_vk = ag_b.getDoublePrecisionValenceKit();
  const ValenceKit<double> aga_vk = ag_a.getDoublePrecisionValenceKit();
  const NonbondedKit<double> aga_nbk = ag_a.getDoublePrecisionNonbondedKit();
  const NonbondedKit<double> agb_nbk = ag_b.getDoublePrecisionNonbondedKit();
  const ValenceKit<double> comp_vk = ag_ab.getDoublePrecisionValenceKit();
  const NonbondedKit<double> comp_nbk = ag_ab.getDoublePrecisionNonbondedKit();
  const ChemicalDetailsKit aga_cdk = ag_a.getChemicalDetailsKit();
  const ChemicalDetailsKit agb_cdk = ag_b.getChemicalDetailsKit();
  const ChemicalDetailsKit comp_cdk = ag_ab.getChemicalDetailsKit();
  std::vector<int> aga_bond_used(aga_vk.nbond, 0), agb_bond_used(agb_vk.nbond, 0);
  std::vector<int> aga_angl_used(aga_vk.nangl, 0), agb_angl_used(agb_vk.nangl, 0);
  std::vector<int> aga_dihe_used(aga_vk.ndihe, 0), agb_dihe_used(agb_vk.ndihe, 0);
  std::vector<int> unmatched_bonds;
  check((n_a * aga_vk.natom) + (n_b * agb_vk.natom), RelationalOperator::EQUAL, comp_vk.natom,
        "The combined number of atoms in two topologies does not meet expectations.", do_tests);
  for (int pos = 0; pos < comp_vk.nbond; pos++) {

    // Find the size of the molecule that the bond is a part of.  Find the relative difference in
    // topological indices of each atom.
    const int mol_idx = comp_cdk.mol_home[comp_vk.bond_i_atoms[pos]];
    const int mol_size = comp_cdk.mol_limits[mol_idx + 1] - comp_cdk.mol_limits[mol_idx];
    const int rel_topl_idx = abs(comp_vk.bond_i_atoms[pos] - comp_vk.bond_j_atoms[pos]);
    const char4 bi_name = comp_cdk.atom_names[comp_vk.bond_i_atoms[pos]];
    const char4 bj_name = comp_cdk.atom_names[comp_vk.bond_j_atoms[pos]];
    
    // Search each topology.  Seek bonds in a molecule of the same size, and if that can be found
    // then bonds with similar atom and residue names and similar relative topological indices.
    const int aga_match = confirmBondMatch(aga_cdk, aga_vk, mol_size, rel_topl_idx, bi_name,
                                           bj_name, aga_bond_used, n_a);
    const int agb_match = confirmBondMatch(agb_cdk, agb_vk, mol_size, rel_topl_idx, bi_name,
                                           bj_name, agb_bond_used, n_b);
    if (aga_match >= 0) {
      aga_bond_used[aga_match] += 1;      
    }
    else if (agb_match >= 0) {
      agb_bond_used[agb_match] += 1;
    }
    else {
      unmatched_bonds.push_back(pos);
    }
  }
  check(unmatched_bonds.size() == 0, "Some bonds (" + std::to_string(unmatched_bonds.size()) +
        ") in a fused topology of " + getBaseName(ag_a.getFileName()) + " and " +
        getBaseName(ag_b.getFileName()) + " could not be traced to their origins in the input "
        "topologies.", do_tests);
  std::vector<int> unmatched_angls;
  for (int pos = 0; pos < comp_vk.nangl; pos++) {

    // Find the size of the molecule that the bond is a part of.  Find the relative difference in
    // topological indices of each atom.
    const int mol_idx = comp_cdk.mol_home[comp_vk.angl_i_atoms[pos]];
    const int mol_size = comp_cdk.mol_limits[mol_idx + 1] - comp_cdk.mol_limits[mol_idx];
    const int ij_topl_idx = abs(comp_vk.angl_i_atoms[pos] - comp_vk.angl_j_atoms[pos]);
    const int jk_topl_idx = abs(comp_vk.angl_j_atoms[pos] - comp_vk.angl_k_atoms[pos]);
    const char4 bi_name = comp_cdk.atom_names[comp_vk.angl_i_atoms[pos]];
    const char4 bj_name = comp_cdk.atom_names[comp_vk.angl_j_atoms[pos]];
    const char4 bk_name = comp_cdk.atom_names[comp_vk.angl_k_atoms[pos]];
    
    // Search each topology.  Seek bonds in a molecule of the same size, and if that can be found
    // then bonds with similar atom and residue names and similar relative topological indices.
    const int aga_match = confirmAngleMatch(aga_cdk, aga_vk, mol_size, ij_topl_idx, jk_topl_idx,
                                            bi_name, bj_name, bk_name, aga_angl_used, n_a);
    const int agb_match = confirmAngleMatch(agb_cdk, agb_vk, mol_size, ij_topl_idx, jk_topl_idx,
                                            bi_name, bj_name, bk_name, agb_angl_used, n_b);
    if (aga_match >= 0) {
      aga_angl_used[aga_match] += 1;      
    }
    else if (agb_match >= 0) {
      agb_angl_used[agb_match] += 1;
    }
    else {
      unmatched_angls.push_back(pos);
    }
  }
  check(unmatched_angls.size() == 0, "A total of " + std::to_string(unmatched_angls.size()) +
        " angles in a fused topology of " + getBaseName(ag_a.getFileName()) + " and " +
        getBaseName(ag_b.getFileName()) + " could not be traced to their origins in the input "
        "topologies.", do_tests);
  std::vector<int> unmatched_dihes;
  for (int pos = 0; pos < comp_vk.ndihe; pos++) {

    // Find the size of the molecule that the bond is a part of.  Find the relative difference in
    // topological indices of each atom.
    const int mol_idx = comp_cdk.mol_home[comp_vk.dihe_i_atoms[pos]];
    const int mol_size = comp_cdk.mol_limits[mol_idx + 1] - comp_cdk.mol_limits[mol_idx];
    const int ij_topl_idx = abs(comp_vk.dihe_i_atoms[pos] - comp_vk.dihe_j_atoms[pos]);
    const int jk_topl_idx = abs(comp_vk.dihe_j_atoms[pos] - comp_vk.dihe_k_atoms[pos]);
    const int kl_topl_idx = abs(comp_vk.dihe_k_atoms[pos] - comp_vk.dihe_l_atoms[pos]);
    const char4 bi_name = comp_cdk.atom_names[comp_vk.dihe_i_atoms[pos]];
    const char4 bj_name = comp_cdk.atom_names[comp_vk.dihe_j_atoms[pos]];
    const char4 bk_name = comp_cdk.atom_names[comp_vk.dihe_k_atoms[pos]];
    const char4 bl_name = comp_cdk.atom_names[comp_vk.dihe_l_atoms[pos]];
    
    // Search each topology.  Seek bonds in a molecule of the same size, and if that can be found
    // then bonds with similar atom and residue names and similar relative topological indices.
    const int aga_match = confirmDihedralMatch(aga_cdk, aga_vk, mol_size, ij_topl_idx,
                                               jk_topl_idx, kl_topl_idx, bi_name, bj_name,
                                               bk_name, bl_name, aga_dihe_used, n_a);
    const int agb_match = confirmDihedralMatch(agb_cdk, agb_vk, mol_size, ij_topl_idx,
                                               jk_topl_idx, kl_topl_idx, bi_name, bj_name,
                                               bk_name, bl_name, agb_dihe_used, n_b);
    if (aga_match >= 0) {
      aga_dihe_used[aga_match] += 1;      
    }
    else if (agb_match >= 0) {
      agb_dihe_used[agb_match] += 1;
    }
    else {
      unmatched_dihes.push_back(pos);
    }
  }
  check(unmatched_dihes.size() == 0, "A total of " + std::to_string(unmatched_dihes.size()) +
        " angles in a fused topology of " + getBaseName(ag_a.getFileName()) + " and " +
        getBaseName(ag_b.getFileName()) + " could not be traced to their origins in the input "
        "topologies.", do_tests);

  // Create a combined PhaseSpace object and test energies.  Test the PhaseSpace combination
  // procedure.
  PhaseSpace comp_ps({ const_cast<PhaseSpace*>(ps_a.getSelfPointer()),
                       const_cast<PhaseSpace*>(ps_b.getSelfPointer()) },
                     { const_cast<AtomGraph*>(ag_a.getSelfPointer()),
                       const_cast<AtomGraph*>(ag_b.getSelfPointer()) }, { n_a, n_b });
  check(comp_ps.getAtomCount(), RelationalOperator::EQUAL,
        (n_a * ag_a.getAtomCount()) + (n_b * ag_b.getAtomCount()), "The number of atoms in a "
        "combined AtomGraph object does not meet expectations.", do_tests);
  check(ag_ab.getBondTermCount(), RelationalOperator::EQUAL,
        (n_a * ag_a.getBondTermCount()) + (n_b * ag_b.getBondTermCount()), "The number of bond "
        "terms in a combined AtomGraph object does not meet expectations.", do_tests);
  check(ag_ab.getAngleTermCount(), RelationalOperator::EQUAL,
        (n_a * ag_a.getAngleTermCount()) + (n_b * ag_b.getAngleTermCount()), "The number of bond "
        "terms in a combined AtomGraph object does not meet expectations.", do_tests);
  check(ag_ab.getDihedralTermCount(), RelationalOperator::EQUAL,
        (n_a * ag_a.getDihedralTermCount()) + (n_b * ag_b.getDihedralTermCount()), "The number of "
        "bond terms in a combined AtomGraph object does not meet expectations.", do_tests);
  PhaseSpace ps_a_copy = ps_a;
  PhaseSpace ps_b_copy = ps_b;
  PhaseSpaceWriter aga_psw = ps_a_copy.data();
  PhaseSpaceWriter agb_psw = ps_b_copy.data();
  PhaseSpaceWriter comp_psw = comp_ps.data();
  const double3 aga_com = centerOfMass(ag_a, ps_a, 0);
  const double3 agb_com = centerOfMass(ag_b, ps_b, 0);
  int atmcon = 0;
  for (int i = 0; i < n_a; i++) {
    const double di = i;
    for (int j = 0; j < aga_psw.natom; j++) {
      comp_psw.xcrd[atmcon] += -500.0 - aga_com.x + (di * dr_a.x);
      comp_psw.ycrd[atmcon] += -500.0 - aga_com.y + (di * dr_a.y);
      comp_psw.zcrd[atmcon] += -500.0 - aga_com.z + (di * dr_a.z);
      comp_psw.xalt[atmcon] += -500.0 - aga_com.x + (di * dr_a.x);
      comp_psw.yalt[atmcon] += -500.0 - aga_com.y + (di * dr_a.y);
      comp_psw.zalt[atmcon] += -500.0 - aga_com.z + (di * dr_a.z);
      atmcon++;
    }
  }
  for (int i = 0; i < n_b; i++) {
    const double di = i;
    for (int j = 0; j < agb_psw.natom; j++) {
      comp_psw.xcrd[atmcon] += 500.0 - agb_com.x + (di * dr_b.x);
      comp_psw.ycrd[atmcon] += 500.0 - agb_com.y + (di * dr_b.y);
      comp_psw.zcrd[atmcon] += 500.0 - agb_com.z + (di * dr_b.z);
      comp_psw.xalt[atmcon] += 500.0 - agb_com.x + (di * dr_b.x);
      comp_psw.yalt[atmcon] += 500.0 - agb_com.y + (di * dr_b.y);
      comp_psw.zalt[atmcon] += 500.0 - agb_com.z + (di * dr_b.z);
      atmcon++;
    }
  }

  // Check valence interactions within the two molecules and their combined structure
  ScoreCard sc(3, 1, 32);
  const double dn_a = n_a;
  const double dn_b = n_b;
  const double bond_nrg_a  = evaluateBondTerms(aga_vk, aga_psw, &sc, EvaluateForce::YES, 0);
  const double bond_nrg_b  = evaluateBondTerms(agb_vk, agb_psw, &sc, EvaluateForce::YES, 1);
  const double bond_nrg_ab = evaluateBondTerms(comp_vk, comp_psw, &sc, EvaluateForce::YES, 2);
  const double angl_nrg_a  = evaluateAngleTerms(aga_vk, aga_psw, &sc, EvaluateForce::YES, 0);
  const double angl_nrg_b  = evaluateAngleTerms(agb_vk, agb_psw, &sc, EvaluateForce::YES, 1);
  const double angl_nrg_ab = evaluateAngleTerms(comp_vk, comp_psw, &sc, EvaluateForce::YES, 2);
  const double2 dihe_nrg_a  = evaluateDihedralTerms(aga_vk, aga_psw, &sc, EvaluateForce::YES, 0);
  const double2 dihe_nrg_b  = evaluateDihedralTerms(agb_vk, agb_psw, &sc, EvaluateForce::YES, 1);
  const double2 dihe_nrg_ab = evaluateDihedralTerms(comp_vk, comp_psw, &sc, EvaluateForce::YES, 2);
  check((dn_a * bond_nrg_a) + (dn_b * bond_nrg_b), RelationalOperator::EQUAL, bond_nrg_ab,
        "The total energy of bonds in a complex of " + std::to_string(n_a) + "  miniproteins "
        "and " + std::to_string(n_b) + " small molecules does not match the sum of energies in "
        "the individual systems.", do_tests);
  check((dn_a * angl_nrg_a) + (dn_b * angl_nrg_b), RelationalOperator::EQUAL, angl_nrg_ab, "The "
        "total energy of bond angles in a complex of " + std::to_string(n_a) + " miniproteins "
        "and " + std::to_string(n_b) + " small molecules does not match the sum of energies in "
        "the individual systems.", do_tests);
  check((dn_a * dihe_nrg_a.x) + (dn_b * dihe_nrg_b.x), RelationalOperator::EQUAL, dihe_nrg_ab.x,
        "The total energy of dihedral terms in a complex of " + std::to_string(n_a) +
        " miniproteins and " + std::to_string(n_b) + " small molecules does not match the sum of "
        "energies in the individual systems.", do_tests);
  check((dn_a * dihe_nrg_a.y) + (dn_b * dihe_nrg_b.y), RelationalOperator::EQUAL, dihe_nrg_ab.y,
        "The total energy of planarity terms in a complex of " + std::to_string(n_a) +
        " miniproteins and " + std::to_string(n_b) + " small molecules does not match the sum of "
        "energies in the individual systems.", do_tests);
  
  // Check 1-4 non-bonded interactions within the two molecules and their combined structure
  const double2 attn_nrg_a  = evaluateAttenuated14Terms(aga_vk, aga_nbk, aga_psw, &sc,
                                                        EvaluateForce::YES, EvaluateForce::YES, 0);
  const double2 attn_nrg_b  = evaluateAttenuated14Terms(agb_vk, agb_nbk, agb_psw, &sc,
                                                        EvaluateForce::YES, EvaluateForce::YES, 1);
  const double2 attn_nrg_ab = evaluateAttenuated14Terms(comp_vk, comp_nbk, comp_psw, &sc,
                                                        EvaluateForce::YES, EvaluateForce::YES, 2);
  check((dn_a * attn_nrg_a.x) + (dn_b * attn_nrg_b.x), RelationalOperator::EQUAL, attn_nrg_ab.x,
        "The total electrostatic energy of attenuated 1:4 non-bonded interactions in a complex "
        "of " + std::to_string(n_a) + " " + std::to_string(ag_a.getAtomCount()) + "-atom and " +
        std::to_string(n_b) + " " + std::to_string(ag_b.getAtomCount()) + "-atom molecules does "
        "not match the sum of energies in the individual systems.");
  check((dn_a * attn_nrg_a.y) + (dn_b * attn_nrg_b.y), RelationalOperator::EQUAL, attn_nrg_ab.y,
        "The total van-der Waals energy of attenuated 1:4 non-bonded interactions in a complex "
        "of " + std::to_string(n_a) + " " + std::to_string(ag_a.getAtomCount()) + "-atom and " +
        std::to_string(n_b) + " " + std::to_string(ag_b.getAtomCount()) + "-atom molecules "
        "does not match the sum of energies in the individual systems.");
  if (n_a == 1 && n_b == 1) {
    const std::vector<double4_16a> attn_pairs_a = getAttenuated14PairList(aga_vk, aga_nbk, "first",
                                                                          do_tests);
    const std::vector<double4_16a> attn_pairs_b = getAttenuated14PairList(agb_vk, agb_nbk,
                                                                          "second", do_tests);
    const std::vector<double4_16a> attn_pairs_ab = getAttenuated14PairList(comp_vk, comp_nbk,
                                                                           "combined", do_tests);
    checkAttenuated14PairLists(attn_pairs_a, attn_pairs_b, attn_pairs_ab,
                               getBaseName(ag_a.getFileName()), getBaseName(ag_b.getFileName()),
                               aga_cdk, agb_cdk, comp_cdk, do_tests);
  }

  // Check the dihedral-related and inferred 1:4 attenuations.
  int indiv_dihe_attn = 0;
  int indiv_infr_attn = 0;
  for (int i = 0; i < aga_vk.ndihe; i++) {
    indiv_dihe_attn += n_a * (aga_vk.dihe14_param_idx[i] == 0);
  }
  for (int i = 0; i < agb_vk.ndihe; i++) {
    indiv_dihe_attn += n_b * (agb_vk.dihe14_param_idx[i] == 0);
  }
  for (int i = 0; i < aga_vk.ninfr14; i++) {
    indiv_infr_attn += n_a * (aga_vk.infr14_param_idx[i] == 0);
  }
  for (int i = 0; i < agb_vk.ninfr14; i++) {
    indiv_infr_attn += n_b * (agb_vk.infr14_param_idx[i] == 0);
  }
  int comp_dihe_attn = 0;
  int comp_infr_attn = 0;
  for (int i = 0; i < comp_vk.ndihe; i++) {
    comp_dihe_attn += (comp_vk.dihe14_param_idx[i] == 0);
  }
  for (int i = 0; i < comp_vk.ninfr14; i++) {
    comp_infr_attn += (comp_vk.infr14_param_idx[i] == 0);
  }
  check(indiv_dihe_attn, RelationalOperator::EQUAL, comp_dihe_attn, "The number of dihedral 1:4 "
        "interactions calculated in a fused topology does not match the sum of 1:4 interaction "
        "counts in the underlying topologies, " + getBaseName(ag_a.getFileName()) + " and " +
        getBaseName(ag_b.getFileName()) + ".", do_tests);
  check(indiv_infr_attn, RelationalOperator::EQUAL, comp_infr_attn, "The number of inferred 1:4 "
        "interactions calculated in a fused topology does not match the sum of 1:4 interaction "
        "counts in the underlying topologies, " + getBaseName(ag_a.getFileName()) + " and " +
        getBaseName(ag_b.getFileName()) + ".", do_tests);

  // Check non-bonded interactions within the two molecules and their combined structure
  const StaticExclusionMask aga_se(ag_a), agb_se(ag_b), comp_se(ag_ab);
  const double2 nonb_nrg_a = evaluateNonbondedEnergy(ag_a, aga_se, &ps_a_copy, &sc,
                                                     EvaluateForce::YES, EvaluateForce::YES, 0);
  const double2 nonb_nrg_b = evaluateNonbondedEnergy(ag_b, agb_se, &ps_b_copy, &sc,
                                                     EvaluateForce::YES, EvaluateForce::YES, 1);
  const double2 nonb_nrg_ab = evaluateNonbondedEnergy(ag_ab, comp_se, &comp_ps, &sc,
                                                      EvaluateForce::YES, EvaluateForce::YES, 2);
  check((dn_a * nonb_nrg_a.x) + (dn_b * nonb_nrg_b.x), RelationalOperator::EQUAL,
        Approx(nonb_nrg_ab.x).margin(2.0), "The electrostatic energy of two combined systems "
        "does not match (within reason) the sum of electrostatic energies of the original "
        "systems.", do_tests);
  check((dn_a * nonb_nrg_a.y) + (dn_b * nonb_nrg_b.y), RelationalOperator::EQUAL,
        Approx(nonb_nrg_ab.y).margin(1.0e-3), "The Lennard-Jones energy of two combined systems "
        "does not match (within reason) the sum of Lennard-Jones energies of the original "
        "systems.", do_tests);

  // Check exclusions within the two molecules and their combined structure
  const LocalExclusionMask aga_lem(ag_a), agb_lem(ag_b), comp_lem(ag_ab);
  std::vector<int2> aga_omitted_exclusions, aga_excess_exclusions;
  std::vector<int2> agb_omitted_exclusions, agb_excess_exclusions;
  std::vector<int2> aga_omitted_local_excl, aga_excess_local_excl;
  std::vector<int2> agb_omitted_local_excl, agb_excess_local_excl;
  int static_a_nexcl  = 0;
  int static_ab_nexcl = 0;
  int local_a_nexcl  = 0;
  int local_ab_nexcl = 0;
  for (int i = 0; i < aga_nbk.natom; i++) {
    for (int j = 0; j <= i; j++) {
      const bool aga_excl = aga_se.testExclusion(i, j);
      const bool comp_excl = comp_se.testExclusion(i, j);
      if (aga_excl == false && comp_excl == true) {
        aga_excess_exclusions.push_back({ i, j });
      }
      if (aga_excl == true && comp_excl == false) {
        aga_omitted_exclusions.push_back({ i, j });
      }
      const bool aga_lcl_excl = aga_lem.testExclusion(i, j);
      const bool comp_lcl_excl = comp_lem.testExclusion(i, j);
      if (aga_lcl_excl == false && comp_lcl_excl == true) {
        aga_excess_local_excl.push_back({ i, j });
      }
      if (aga_lcl_excl == true && comp_lcl_excl == false) {
        aga_omitted_local_excl.push_back({ i, j });
      }
      static_a_nexcl += static_cast<int>(aga_excl);
      local_a_nexcl += static_cast<int>(aga_lcl_excl);
      static_ab_nexcl += static_cast<int>(comp_excl);
      local_ab_nexcl += static_cast<int>(comp_lcl_excl);
    }
  }
  const std::string omit_a_msg = listMissingPairs(aga_omitted_exclusions, aga_cdk);
  const std::string omit_local_a_msg = listMissingPairs(aga_omitted_local_excl, agb_cdk);
  const std::string xs_a_msg = listMissingPairs(aga_excess_exclusions, aga_cdk);
  const std::string xs_local_a_msg = listMissingPairs(aga_excess_local_excl, agb_cdk);
  check(aga_omitted_exclusions.size() == 0, "Some exclusions (total " +
        std::to_string(aga_omitted_exclusions.size()) + " of " + std::to_string(static_a_nexcl) +
        ") in the first molecule (" + std::to_string(ag_a.getAtomCount()) + " atoms in all) were "
        "excluded from the complex topology (" + std::to_string(ag_ab.getAtomCount()) +
        " atoms).  This happened when comparing to a STATIC exclusion mask.  " + omit_a_msg,
        do_tests);
  check(aga_excess_exclusions.size() == 0, "Some pairs (total " +
        std::to_string(aga_excess_exclusions.size()) + " of " + std::to_string(static_a_nexcl) +
        ") in the first molecule (" + std::to_string(ag_a.getAtomCount()) + " atoms in all) "
        "were excluded in the complex topology (" + std::to_string(ag_ab.getAtomCount()) +
        " atoms) by mistake.  This happened when comparing to a STATIC exclusion mask.  " +
        xs_a_msg, do_tests);
  check(aga_omitted_local_excl.size() == 0, "Some exclusions (total " +
        std::to_string(aga_omitted_local_excl.size()) + " of " +
        std::to_string(local_a_nexcl) + ") in the first molecule (" +
        std::to_string(ag_a.getAtomCount()) + " in all) were excluded from the complex "
        "topology (" + std::to_string(ag_ab.getAtomCount()) + " atoms).  This happened when "
        "comparing to a LOCAL exclusion mask.  " + omit_local_a_msg, do_tests);
  check(aga_excess_local_excl.size() == 0, "Some pairs (total " +
        std::to_string(aga_excess_local_excl.size()) + " of " + std::to_string(local_a_nexcl) +
        ") in the first molecule (" + std::to_string(ag_a.getAtomCount()) + " atoms in all) were "
        "excluded in the complex topology (" + std::to_string(ag_ab.getAtomCount()) + " atoms) "
        "by mistake.  This happened when comparing to a LOCAL exclusion mask.  " +
        xs_local_a_msg, do_tests);
  int static_b_nexcl  = 0;
  int local_b_nexcl  = 0;
  for (int i = 0; i < agb_nbk.natom; i++) {
    for (int j = 0; j <= i; j++) {
      const bool agb_excl = agb_se.testExclusion(i, j);
      const bool comp_excl = comp_se.testExclusion(i + (n_a * aga_nbk.natom),
                                                   j + (n_a * aga_nbk.natom));
      if (agb_excl == false && comp_excl == true) {
        agb_excess_exclusions.push_back({ i, j });
      }
      if (agb_excl == true && comp_excl == false) {
        agb_omitted_exclusions.push_back({ i, j });
      }
      const bool agb_lcl_excl = agb_lem.testExclusion(i, j);
      const bool comp_lcl_excl = comp_lem.testExclusion(i + (n_a * aga_nbk.natom),
                                                        j + (n_a * aga_nbk.natom));
      if (agb_lcl_excl == false && comp_lcl_excl == true) {
        agb_excess_local_excl.push_back({ i, j });
      }
      if (agb_lcl_excl == true && comp_lcl_excl == false) {
        agb_omitted_local_excl.push_back({ i, j });
      }
      static_b_nexcl += static_cast<int>(agb_excl);
      local_b_nexcl += static_cast<int>(agb_lcl_excl);
    }
  }
  const std::string omit_b_msg = listMissingPairs(agb_omitted_exclusions, aga_cdk);
  const std::string omit_local_b_msg = listMissingPairs(agb_omitted_local_excl, agb_cdk);
  const std::string xs_b_msg = listMissingPairs(agb_excess_exclusions, aga_cdk);
  const std::string xs_local_b_msg = listMissingPairs(agb_excess_local_excl, agb_cdk);
  check(aga_omitted_exclusions.size() == 0, "Some exclusions (total " +
        std::to_string(agb_omitted_exclusions.size()) + " of " + std::to_string(static_b_nexcl) +
        ") in the second molecule (" + std::to_string(ag_b.getAtomCount()) + " atoms in all) were "
        "excluded from the complex topology (" + std::to_string(ag_ab.getAtomCount()) +
        " atoms).  This happened when comparing to a STATIC exclusion mask.  " + omit_b_msg,
        do_tests);
  check(agb_excess_exclusions.size() == 0, "Some pairs (total " +
        std::to_string(agb_excess_exclusions.size()) + " of " + std::to_string(static_b_nexcl) +
        ") in the second molecule (" + std::to_string(ag_b.getAtomCount()) + " atoms in all) "
        "were excluded in the complex topology (" + std::to_string(ag_ab.getAtomCount()) +
        " atoms) by mistake.  This happened when comparing to a STATIC exclusion mask.  " +
        xs_b_msg, do_tests);
  check(agb_omitted_local_excl.size() == 0, "Some exclusions (total " +
        std::to_string(agb_omitted_local_excl.size()) + " of " +
        std::to_string(local_b_nexcl) + ") in the second molecule (" +
        std::to_string(ag_b.getAtomCount()) + " in all) were excluded from the complex "
        "topology (" + std::to_string(ag_ab.getAtomCount()) + " atoms).  This happened when "
        "comparing to a LOCAL exclusion mask.  " + omit_local_b_msg, do_tests);
  check(agb_excess_local_excl.size() == 0, "Some pairs (total " +
        std::to_string(agb_excess_local_excl.size()) + " of " + std::to_string(local_b_nexcl) +
        ") in the second molecule (" + std::to_string(ag_b.getAtomCount()) + " atoms in all) were "
        "excluded in the complex topology (" + std::to_string(ag_ab.getAtomCount()) + " atoms) "
        "by mistake.  This happened when comparing to a LOCAL exclusion mask.  " +
        xs_local_b_msg, do_tests);
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Initialize the test environment
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1: Test bond modifications
  section("Test bond parameter modifications");

  // Section 2: Test angle modifications
  section("Test angle parameter modifications");

  // Section 3: Test dihedral modifications
  section("Test dihedral parameter modifications");

  // Section 4: Test topology staging
  section("Test topology staging");
  
  // Load a series of topologies, then modify them a bit at a time
  const char osc = osSeparator();
  const std::string all_base_name   = oe.getStormmSourcePath() + osc + "test";
  const std::string topology_home   = "Topology";
  const std::string trajectory_home = "Trajectory";
  const std::string chemistry_home  = "Chemistry";
  const std::vector<std::string> top_names = {
    topology_home + osc + "stereo_L1", topology_home + osc + "stereo_L1_vs",
    topology_home + osc + "symmetry_L1", topology_home + osc + "symmetry_L1_vs",
    topology_home + osc + "bromobenzene", topology_home + osc + "bromobenzene_vs",
    chemistry_home + osc + "lig1_c8h8", chemistry_home + osc + "lig2_c8h8",
    chemistry_home + osc + "lig3_1fsj", chemistry_home + osc + "morphine_like"
  };
  const std::vector<std::string> crd_names = {
    trajectory_home + osc + "stereo_L1", trajectory_home + osc + "stereo_L1_vs",
    trajectory_home + osc + "symmetry_L1", trajectory_home + osc + "symmetry_L1_vs",
    trajectory_home + osc + "bromobenzene", trajectory_home + osc + "bromobenzene_vs",
    chemistry_home + osc + "lig1_c8h8", chemistry_home + osc + "lig2_c8h8",
    chemistry_home + osc + "lig3_1fsj", chemistry_home + osc + "morphine_like",
  };
  TestSystemManager tsm(all_base_name, "top", top_names, all_base_name, "inpcrd", crd_names);
  const TestPriority do_tests = tsm.getTestingStatus();
  std::vector<AtomGraph> ag_list;
  std::vector<PhaseSpace> ps_list;
  std::vector<StaticExclusionMask> se_list;
  const int system_count = tsm.getSystemCount();
  ag_list.reserve(system_count);
  ps_list.reserve(system_count);
  se_list.reserve(system_count);
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    ag_list.push_back(tsm.exportAtomGraph(i));
    ps_list.push_back(tsm.exportPhaseSpace(i));
    se_list.emplace_back(StaticExclusionMask(ag_list[i]));
  }

  // Evaluate energies with the canonical topologies
  ScoreCard sc_orig(system_count);
  MinimizeControls mincon;
  mincon.setSteepestDescentCycles(25);
  mincon.setTotalCycles(50);
  for (size_t i = 0; i < system_count; i++) {
    RestraintApparatus ra(&ag_list[i]);
    ScoreCard min_sc = minimize(&ps_list[i], &ag_list[i], se_list[i], mincon, 30);
    evalNonbValeMM(&ps_list[i], &sc_orig, ag_list[i], se_list[i], EvaluateForce::NO, i);
  }
  const std::vector<double> bond_e0 = sc_orig.reportInstantaneousStates(StateVariable::BOND);
  const std::vector<double> angl_e0 = sc_orig.reportInstantaneousStates(StateVariable::ANGLE);
  const std::vector<double> dihe_e0 =
    sc_orig.reportInstantaneousStates(StateVariable::PROPER_DIHEDRAL);
  const std::vector<double> impr_e0 =
    sc_orig.reportInstantaneousStates(StateVariable::IMPROPER_DIHEDRAL);
  ForceFieldElement ca_ha_bond(ParameterKind::BOND, stringToChar4("ca"), stringToChar4("ha"));
  ForceFieldElement bl_bm_bond(ParameterKind::BOND, stringToChar4("bl"), stringToChar4("bm"));
  ca_ha_bond.setStiffness(100.0);
  ca_ha_bond.setEquilibrium(1.5);
  bl_bm_bond.setStiffness(120.0);
  bl_bm_bond.setEquilibrium(1.7);
  std::vector<AtomGraph> ag_mods;
  std::vector<PhaseSpace> ps_mods;
  for (size_t i = 0; i < system_count; i++) {
    ag_mods.emplace_back(ag_list[i]);
    ps_mods.emplace_back(ps_list[i]);
  }
  ScoreCard sc_mods(system_count);
  for (size_t i = 0; i < system_count; i++) {
    ca_ha_bond.apply(&ag_mods[i]);
    bl_bm_bond.apply(&ag_mods[i]);
    RestraintApparatus ra(&ag_mods[i]);
    ScoreCard min_sc = minimize(&ps_mods[i], &ag_mods[i], se_list[i], mincon, 30);
    evalNonbValeMM(&ps_mods[i], &sc_mods, ag_mods[i], se_list[i], EvaluateForce::NO, i);
  }
  const std::vector<double> bond_em = sc_mods.reportInstantaneousStates(StateVariable::BOND);
  const std::vector<double> angl_em= sc_mods.reportInstantaneousStates(StateVariable::ANGLE);
  const std::vector<double> dihe_em =
    sc_mods.reportInstantaneousStates(StateVariable::PROPER_DIHEDRAL);
  const std::vector<double> impr_em =
    sc_mods.reportInstantaneousStates(StateVariable::IMPROPER_DIHEDRAL);

  // Try deconstructing a topology and putting it back together
  section(4);
  AtomGraphStage ags(ag_list[0]);
  AtomGraph nchiral_sys = ags.exportTopology();
  check(nchiral_sys.getAtomCount(), RelationalOperator::EQUAL, ag_list[0].getAtomCount(),
        "The AtomGraphStage does not convey the atom count correctly when initialized based on an "
        "existing topology.", do_tests);
  AtomGraphStage ags_vs(ag_list[5]);
  AtomGraph nbromo_vs_sys = ags_vs.exportTopology();
  check(nbromo_vs_sys.getAtomCount(), RelationalOperator::EQUAL, ag_list[5].getAtomCount(),
        "The AtomGraphStage does not convey the atom count correctly when initialized based on an "
        "existing topology with virtual sites.", do_tests);
  AtomGraphStage ion_web(100, { 0, 10, 20, 30, 40, 50, 100 });
  AtomGraph nion_sys = ion_web.exportTopology();
  check(nion_sys.getAtomCount(), RelationalOperator::EQUAL, 100, "The AtomGraphStage does not "
        "convey the atom count correctly when initialized based on a number of particles.");

  // Test the Lennard-Jones analysis
  const std::vector<std::string> ffpro = { "trpcage", "tamavidin", "stereo_L1", "stereo_L1_vs" };
  TestSystemManager combi(all_base_name + osc + "Topology", "top", ffpro,
                          all_base_name + osc + "Trajectory", "inpcrd", ffpro);
  const AtomGraph& trpi_ag = combi.getTopologyReference(0);
  const AtomGraph& tama_ag = combi.getTopologyReference(1);
  LennardJonesAnalysis lja(tama_ag.getDoublePrecisionNonbondedKit(),
                           tama_ag.getAtomTypeNameTable());
  lja.addSet(trpi_ag.getDoublePrecisionNonbondedKit(), trpi_ag.getAtomTypeNameTable());
  const char4 ambiguous_hydrogen = stringToChar4("H");
  CHECK_THROWS_SOFT(std::vector<char4> h_names = lja.getLJAliases(ambiguous_hydrogen),
                    "A request for an ambiguous atom type's aliases went through unreported.",
                    combi.getTestingStatus());
  CHECK_THROWS_SOFT(double2 ljab_h = lja.getLJParameters(ambiguous_hydrogen),
                    "A request for an ambiguous atom type's parameters went through unreported.",
                    combi.getTestingStatus());
  std::vector<char4> h_names = lja.getLJAliases(stringToChar4("TP"), ExceptionResponse::WARN);
  check(h_names.size(), RelationalOperator::EQUAL, 6, "The number of atom type aliases found for "
        "atom type \"TP\" does not meet expectations.", combi.getTestingStatus());
  const NonbondedKit<double> tama_nbk = tama_ag.getDoublePrecisionNonbondedKit();
  const NonbondedKit<double> trpi_nbk = trpi_ag.getDoublePrecisionNonbondedKit();
  const size_t nlj_tama_sq = tama_nbk.n_lj_types * tama_nbk.n_lj_types;
  const size_t nlj_trpi_sq = trpi_nbk.n_lj_types * trpi_nbk.n_lj_types;
  std::vector<double> tama_lja(nlj_tama_sq), tama_ljb(nlj_tama_sq);
  std::vector<double> trpi_lja(nlj_tama_sq), trpi_ljb(nlj_tama_sq);
  std::vector<double> ani_lja(nlj_tama_sq), ani_ljb(nlj_tama_sq);
  std::vector<double> anii_lja(nlj_tama_sq), anii_ljb(nlj_tama_sq);
  for (size_t i = 0; i < nlj_tama_sq; i++) {
    tama_lja[i] = tama_nbk.lja_coeff[i];
    tama_ljb[i] = tama_nbk.ljb_coeff[i];
  }
  size_t ijcon = 0;
  for (int i = 0; i < tama_nbk.n_lj_types; i++) {
    for (int j = 0; j < tama_nbk.n_lj_types; j++) {
      const double3 pij = lja.getLJCoefficients(i, j);
      ani_lja[ijcon] = pij.x;
      ani_ljb[ijcon] = pij.y;
      ijcon++;
    }
  }
  check(tama_lja, RelationalOperator::EQUAL, ani_lja, "The Lennard-Jones analysis does not "
        "reflect the A coefficients of the first system (" + std::to_string(tama_nbk.n_lj_types) +
        " Lennard-Jones types) as expected.", combi.getTestingStatus());
  check(tama_ljb, RelationalOperator::EQUAL, ani_ljb, "The Lennard-Jones analysis does not "
        "reflect the B coefficients of the first system (" + std::to_string(tama_nbk.n_lj_types) +
        " Lennard-Jones types) as expected.", combi.getTestingStatus());
  const std::vector<int> trpi_corr_ans = { 0, 14, 2, 3, 9, 7, 8, 4, 11, 5, 6, 10, 15 };
  check(lja.getSetCorrespondence(1), RelationalOperator::EQUAL, trpi_corr_ans, "The Lennard-Jones "
        "type correspondence obtained for an added system based on " +
        getBaseName(trpi_ag.getFileName()) + " did not meet expectations.",
        combi.getTestingStatus()); 

  // Test parameter table combination methods
  Xoroshiro128pGenerator xrs;
  const std::vector<double> phi_a = uniformRand(&xrs, 50, 1.0);
  const std::vector<double> amp_a = uniformRand(&xrs, 50, 1.0);
  const std::vector<double> per_a(50, 2.0);
  std::vector<double> phi_b(25), amp_b(25), per_b(25);
  for (int i = 0; i < 25; i++) {
    phi_b[i] = (i >= 10 && i < 15) ? xrs.uniformRandomNumber() - 1.0 : phi_a[2 * i];
    amp_b[i] = amp_a[2 * i];
    per_b[i] = per_a[2 * i];
  }
  ParameterUnion<double> prmu(phi_a, amp_a, per_a, phi_b, amp_b, per_b);
  check(prmu.getUniqueParameterCount(), RelationalOperator::EQUAL, 55, "The number of unique "
        "parameters in a union of 50 and 25, with 20 overlapping, does not meet expectations.");
  const std::vector<int> cmap_dim_a = { 24, 24, 24, 8, 24, 8, 24 };
  const std::vector<int> cmap_dim_b = { 24, 16, 24, 8, 24 };
  std::vector<int> cmap_bounds_a(cmap_dim_a.size() + 1, 0);
  std::vector<int> cmap_bounds_b(cmap_dim_b.size() + 1, 0);
  int sq_sum_a = 0;
  for (size_t i = 0; i < cmap_dim_a.size(); i++) {
    sq_sum_a += cmap_dim_a[i] * cmap_dim_a[i];
    cmap_bounds_a[i + 1] = sq_sum_a;
  }
  int sq_sum_b = 0;
  for (size_t i = 0; i < cmap_dim_b.size(); i++) {
    sq_sum_b += cmap_dim_b[i] * cmap_dim_b[i];
    cmap_bounds_b[i + 1] = sq_sum_b;
  }
  const std::vector<double> cmap_val_a = uniformRand(&xrs, sq_sum_a, 1.0);
  std::vector<double> cmap_val_b = uniformRand(&xrs, sq_sum_b, 1.0);
  for (size_t i = 1; i < cmap_dim_b.size() - 1; i++) {
    if (cmap_dim_a[i] == cmap_dim_b[i]) {
      for (int j = cmap_bounds_a[i]; j < cmap_bounds_a[i + 1]; j++) {
        cmap_val_b[cmap_bounds_b[i] + j - cmap_bounds_a[i]] = cmap_val_a[j];
      }
    }
  }
  CmapSurfaceUnion cmapu(cmap_val_a, cmap_dim_a, cmap_val_b, cmap_dim_b);
  check(cmapu.getUniqueSurfaceCount(), RelationalOperator::EQUAL, 10, "The number of unique "
        "CMAP surfaces found between two sets of lengths seven and five, with two surfaces in "
        "common, does not meet expectations.");

  // Test topology combinations
  const AtomGraph& ster_ag = combi.getTopologyReference(2);
  const PhaseSpace trpi_ps = combi.exportPhaseSpace(0);
  const PhaseSpace ster_ps = combi.exportPhaseSpace(2);
  testTopologyFusion(ster_ag, ster_ag, ster_ps, ster_ps, 1, 1, { 1000.0, 1000.0, 1000.0 },
                     { 1000.0, 1000.0, 1000.0 }, combi.getTestingStatus());
  testTopologyFusion(trpi_ag, ster_ag, trpi_ps, ster_ps, 1, 1, { 1000.0, 1000.0, 1000.0 },
                     { 1000.0, 1000.0, 1000.0 }, combi.getTestingStatus());
  testTopologyFusion(trpi_ag, ster_ag, trpi_ps, ster_ps, 1, 2, { 1000.0, 1000.0, 1000.0 },
                     { 5000.0, 5000.0, -5000.0 }, combi.getTestingStatus());
  const AtomGraph& stvs_ag = combi.getTopologyReference(3);
  const PhaseSpace& stvs_ps = combi.exportPhaseSpace(3);
  testTopologyFusion(stvs_ag, stvs_ag, stvs_ps, stvs_ps, 1, 1, { 1000.0, -1000.0, 1000.0 },
                     { 5000.0, 5000.0, -5000.0 }, combi.getTestingStatus());
  testTopologyFusion(trpi_ag, stvs_ag, trpi_ps, stvs_ps, 1, 1, { 1000.0, 1000.0, 1000.0 },
                     { 1000.0, 1000.0, 1000.0 }, combi.getTestingStatus());
  testTopologyFusion(trpi_ag, stvs_ag, trpi_ps, stvs_ps, 2, 2, { 1000.0, -1000.0, 1000.0 },
                     { 5000.0, 5000.0, -5000.0 }, combi.getTestingStatus());

  // Test topology extraction
  std::vector<bool> trpi_part_mask(trpi_ag.getAtomCount(), false);
  for (int i = 0; i < 100; i++) {
    trpi_part_mask[i] = true;
  }
  AtomGraph trpi_part_ag(trpi_ag, trpi_part_mask);

  // Test topology printing
  const std::string pruned_top_name = oe.getTemporaryDirectoryPath() + osc + "trpi_part_ag.top";
  trpi_part_ag.printToFile(pruned_top_name);
  oe.logFileCreated(pruned_top_name);

  // Summary evaluation
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
