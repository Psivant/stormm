#include <algorithm>
#include <cmath>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Chemistry/periodic_table.h"
#include "Chemistry/znumber.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Math/matrix_ops.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "UnitTesting/approx.h"
#include "atomgraph_refinement.h"

namespace stormm {
namespace topology {

using chemistry::massToZNumber;
using card::HybridTargetLevel;
using stmath::crossProduct;
using stmath::magnitude;
using stmath::matrixMultiply;
using stmath::maxValue;
using stmath::minValue;
using stmath::PrefixSumType;
using stmath::prefixSumInPlace;
using stmath::sum;
using parse::char4ToString;
using parse::NumberFormat;
using parse::PolyNumeric;
using parse::realToString;
using symbols::tetrahedral_angle;
using testing::Approx;

//-------------------------------------------------------------------------------------------------
BasicValenceTable::BasicValenceTable() :
    total_bonds{0}, total_angls{0}, total_dihes{0},
    bond_i_atoms{}, bond_j_atoms{}, bond_param_idx{}, angl_i_atoms{}, angl_j_atoms{},
    angl_k_atoms{}, angl_param_idx{}, dihe_i_atoms{}, dihe_j_atoms{}, dihe_k_atoms{},
    dihe_l_atoms{}, dihe_param_idx{}, bond_mods{}, angl_mods{}, dihe_mods{}, bond_assigned_atoms{},
    bond_assigned_index{}, bond_assigned_terms{}, bond_assigned_bounds{}, angl_assigned_atoms{},
    angl_assigned_index{}, angl_assigned_terms{}, angl_assigned_bounds{}, dihe_assigned_atoms{},
    dihe_assigned_index{}, dihe_assigned_terms{}, dihe_assigned_bounds{}, bond_assigned_mods{},
    angl_assigned_mods{}, dihe_assigned_mods{}
{}

//-------------------------------------------------------------------------------------------------
BasicValenceTable::BasicValenceTable(const int natom_in, const int nbond_in,
                                     const int nangl_in, const int ndihe_in,
                                     const std::vector<int> &bond_i_atoms_in,
                                     const std::vector<int> &bond_j_atoms_in,
                                     const std::vector<int> &bond_param_idx_in,
                                     const std::vector<int> &angl_i_atoms_in,
                                     const std::vector<int> &angl_j_atoms_in,
                                     const std::vector<int> &angl_k_atoms_in,
                                     const std::vector<int> &angl_param_idx_in,
                                     const std::vector<int> &dihe_i_atoms_in,
                                     const std::vector<int> &dihe_j_atoms_in,
                                     const std::vector<int> &dihe_k_atoms_in,
                                     const std::vector<int> &dihe_l_atoms_in,
                                     const std::vector<int> &dihe_param_idx_in) :
    BasicValenceTable()
{
  total_bonds = nbond_in;
  total_angls = nangl_in;
  total_dihes = ndihe_in;
  bond_i_atoms.resize(nbond_in);
  bond_j_atoms.resize(nbond_in);
  bond_param_idx.resize(nbond_in);
  angl_i_atoms.resize(nangl_in);
  angl_j_atoms.resize(nangl_in);
  angl_k_atoms.resize(nangl_in);
  angl_param_idx.resize(nangl_in);
  dihe_i_atoms.resize(ndihe_in);
  dihe_j_atoms.resize(ndihe_in);
  dihe_k_atoms.resize(ndihe_in);
  dihe_l_atoms.resize(ndihe_in);
  dihe_param_idx.resize(ndihe_in);
  bond_mods.resize(nbond_in);
  angl_mods.resize(nangl_in);
  dihe_mods.resize(ndihe_in);
  bond_assigned_atoms.resize(    nbond_in);
  angl_assigned_atoms.resize(2 * nangl_in);
  dihe_assigned_atoms.resize(3 * ndihe_in);
  bond_assigned_index.resize(nbond_in);
  angl_assigned_index.resize(nangl_in);
  dihe_assigned_index.resize(ndihe_in);
  bond_assigned_terms.resize(nbond_in);
  angl_assigned_terms.resize(nangl_in);
  dihe_assigned_terms.resize(ndihe_in);
  bond_assigned_bounds.resize(natom_in + 1);
  angl_assigned_bounds.resize(natom_in + 1);
  dihe_assigned_bounds.resize(natom_in + 1);
  bond_assigned_mods.resize(nbond_in);
  angl_assigned_mods.resize(nangl_in);
  dihe_assigned_mods.resize(ndihe_in);

  // Fill out the valence term information
  bool bonds_provided = false;
  bool angls_provided = false;
  bool dihes_provided = false;
  if (bond_i_atoms.size() == nbond_in && bond_j_atoms.size() == nbond_in &&
      bond_param_idx_in.size() == nbond_in) {
    for (int pos = 0; pos < nbond_in; pos++) {
      bond_i_atoms[pos] = bond_i_atoms_in[pos];
      bond_j_atoms[pos] = bond_j_atoms_in[pos];
      bond_param_idx[pos] = bond_param_idx_in[pos];
    }
    bonds_provided = true;
  }
  if (angl_i_atoms.size() == nangl_in && angl_j_atoms.size() == nangl_in &&
      angl_k_atoms.size() == nangl_in && angl_param_idx_in.size() == nangl_in) {
    for (int pos = 0; pos < nangl_in; pos++) {
      angl_i_atoms[pos] = angl_i_atoms_in[pos];
      angl_j_atoms[pos] = angl_j_atoms_in[pos];
      angl_k_atoms[pos] = angl_k_atoms_in[pos];
      angl_param_idx[pos] = angl_param_idx_in[pos];
    }
    angls_provided = true;
  }
  if (dihe_i_atoms.size() == ndihe_in && dihe_j_atoms.size() == ndihe_in &&
      dihe_k_atoms.size() == ndihe_in && dihe_l_atoms.size() == ndihe_in &&
      dihe_param_idx_in.size() == ndihe_in) {
    for (int pos = 0; pos < ndihe_in; pos++) {
      dihe_i_atoms[pos] = dihe_i_atoms_in[pos];
      dihe_j_atoms[pos] = dihe_j_atoms_in[pos];
      dihe_k_atoms[pos] = dihe_k_atoms_in[pos];
      dihe_l_atoms[pos] = dihe_l_atoms_in[pos];
      dihe_param_idx[pos] = dihe_param_idx_in[pos];
    }
    dihes_provided = true;
  }
  if (bonds_provided && angls_provided && dihes_provided) {
    makeAtomAssignments();
  }
}

//-------------------------------------------------------------------------------------------------
void BasicValenceTable::makeAtomAssignments() {

  // Bonds go to the first atom in the pair (i), angles go to the center atom (j),
  // and dihedrals are the responsibility of the second atom (j) in the list.
  for (int pos = 0; pos < total_bonds; pos++) {
    bond_assigned_bounds[bond_i_atoms[pos]] += 1;
  }
  for (int pos = 0; pos < total_angls; pos++) {
    angl_assigned_bounds[angl_j_atoms[pos]] += 1;
  }
  for (int pos = 0; pos < total_dihes; pos++) {
    dihe_assigned_bounds[dihe_k_atoms[pos]] += 1;
  }
  const PrefixSumType pfx_excl = PrefixSumType::EXCLUSIVE;
  prefixSumInPlace(&bond_assigned_bounds, pfx_excl, "basicValenceIndexing");
  prefixSumInPlace(&angl_assigned_bounds, pfx_excl, "basicValenceIndexing");
  prefixSumInPlace(&dihe_assigned_bounds, pfx_excl, "basicValenceIndexing");

  // Populate the bond, angle, and dihedral assignements.  Dihedrals are assigned to the third
  // atom in the quartet.  In Amber-format topologies this is the central atom of an improper.
  for (int pos = 0; pos < total_bonds; pos++) {
    const int control_atom = bond_i_atoms[pos];
    const int fill_slot = bond_assigned_bounds[control_atom];
    bond_assigned_atoms[fill_slot]     = bond_j_atoms[pos];
    bond_assigned_index[fill_slot]     = bond_param_idx[pos];
    bond_assigned_terms[fill_slot]     = pos;
    bond_assigned_mods[fill_slot]      = bond_mods[pos];
    bond_assigned_bounds[control_atom] = fill_slot + 1;
  }
  for (int pos = 0; pos < total_angls; pos++) {
    const int control_atom = angl_j_atoms[pos];
    const int fill_slot = angl_assigned_bounds[control_atom];
    angl_assigned_atoms[ 2 * fill_slot     ] = angl_i_atoms[pos];
    angl_assigned_atoms[(2 * fill_slot) + 1] = angl_k_atoms[pos];
    angl_assigned_index[fill_slot]           = angl_param_idx[pos];
    angl_assigned_terms[fill_slot]           = pos;
    angl_assigned_mods[fill_slot]            = angl_mods[pos];
    angl_assigned_bounds[control_atom]       = fill_slot + 1;
  }
  for (int pos = 0; pos < total_dihes; pos++) {
    const int control_atom = dihe_k_atoms[pos];
    const int fill_slot = dihe_assigned_bounds[control_atom];
    dihe_assigned_atoms[ 3 * fill_slot     ] = dihe_i_atoms[pos];
    dihe_assigned_atoms[(3 * fill_slot) + 1] = dihe_j_atoms[pos];
    dihe_assigned_atoms[(3 * fill_slot) + 2] = dihe_l_atoms[pos];
    dihe_assigned_index[fill_slot]           = dihe_param_idx[pos];
    dihe_assigned_terms[fill_slot]           = pos;
    dihe_assigned_mods[fill_slot]            = dihe_mods[pos];
    dihe_assigned_bounds[control_atom]       = fill_slot + 1;
  }

  // Trim the prefix sums--they have been advanced one atom's register
  const int atom_count = static_cast<int>(bond_assigned_bounds.size()) - 1;
  for (int pos = atom_count; pos > 0; pos--) {
    bond_assigned_bounds[pos] = bond_assigned_bounds[pos - 1];
    angl_assigned_bounds[pos] = angl_assigned_bounds[pos - 1];
    dihe_assigned_bounds[pos] = dihe_assigned_bounds[pos - 1];
  }
  bond_assigned_bounds[0] = 0;
  angl_assigned_bounds[0] = 0;
  dihe_assigned_bounds[0] = 0;
}
  
//-------------------------------------------------------------------------------------------------
CharmmValenceTable::CharmmValenceTable() :
    total_ub_angles{0}, total_impropers{0}, total_cmaps{0}, ubrd_i_atoms{}, ubrd_k_atoms{},
    ubrd_param_idx{}, impr_i_atoms{}, impr_j_atoms{}, impr_k_atoms{}, impr_l_atoms{},
    impr_param_idx{}, cmap_i_atoms{}, cmap_j_atoms{}, cmap_k_atoms{}, cmap_l_atoms{},
    cmap_m_atoms{}, cmap_param_idx{}, ubrd_assigned_atoms{}, ubrd_assigned_index{},
    ubrd_assigned_terms{}, ubrd_assigned_bounds{}, impr_assigned_atoms{}, impr_assigned_index{},
    impr_assigned_terms{}, impr_assigned_bounds{}, cmap_assigned_atoms{}, cmap_assigned_index{},
    cmap_assigned_terms{}, cmap_assigned_bounds{}
{}

//-------------------------------------------------------------------------------------------------
CharmmValenceTable::CharmmValenceTable(const int natom_in, const int nubrd_in,
                                       const int nimpr_in, const int ncmap_in,
                                       const std::vector<int> &ubrd_i_atoms_in,
                                       const std::vector<int> &ubrd_k_atoms_in,
                                       const std::vector<int> &ubrd_param_idx_in,
                                       const std::vector<int> &impr_i_atoms_in,
                                       const std::vector<int> &impr_j_atoms_in,
                                       const std::vector<int> &impr_k_atoms_in,
                                       const std::vector<int> &impr_l_atoms_in,
                                       const std::vector<int> &impr_param_idx_in,
                                       const std::vector<int> &cmap_i_atoms_in,
                                       const std::vector<int> &cmap_j_atoms_in,
                                       const std::vector<int> &cmap_k_atoms_in,
                                       const std::vector<int> &cmap_l_atoms_in,
                                       const std::vector<int> &cmap_m_atoms_in,
                                       const std::vector<int> &cmap_param_idx_in) :
    CharmmValenceTable()
{
  total_ub_angles = nubrd_in;
  total_impropers = nimpr_in;
  total_cmaps     = ncmap_in;
  ubrd_i_atoms.resize(nubrd_in);
  ubrd_k_atoms.resize(nubrd_in);
  ubrd_param_idx.resize(nubrd_in);
  impr_i_atoms.resize(nimpr_in);
  impr_j_atoms.resize(nimpr_in);
  impr_k_atoms.resize(nimpr_in);
  impr_l_atoms.resize(nimpr_in);
  impr_param_idx.resize(nimpr_in);
  cmap_i_atoms.resize(ncmap_in);
  cmap_j_atoms.resize(ncmap_in);
  cmap_k_atoms.resize(ncmap_in);
  cmap_l_atoms.resize(ncmap_in);
  cmap_m_atoms.resize(ncmap_in);
  cmap_param_idx.resize(ncmap_in);
  ubrd_assigned_atoms.resize(    nubrd_in);
  impr_assigned_atoms.resize(3 * nimpr_in);
  cmap_assigned_atoms.resize(4 * ncmap_in);
  ubrd_assigned_index.resize(nubrd_in);
  impr_assigned_index.resize(nimpr_in);
  cmap_assigned_index.resize(ncmap_in);
  ubrd_assigned_terms.resize(nubrd_in);
  impr_assigned_terms.resize(nimpr_in);
  cmap_assigned_terms.resize(ncmap_in);
  ubrd_assigned_bounds.resize(natom_in + 1, 0);
  impr_assigned_bounds.resize(natom_in + 1, 0);
  cmap_assigned_bounds.resize(natom_in + 1, 0);

  // Fill out the valecne term information
  bool ubrd_provided = false;
  bool impr_provided = false;
  bool cmap_provided = false;
  if (ubrd_i_atoms_in.size() == nubrd_in && ubrd_k_atoms_in.size() == nubrd_in &&
      ubrd_param_idx_in.size() == nubrd_in) {
    for (int pos = 0; pos < nubrd_in; pos++) {
      ubrd_i_atoms[pos] = ubrd_i_atoms_in[pos];
      ubrd_k_atoms[pos] = ubrd_k_atoms_in[pos];
      ubrd_param_idx[pos] = ubrd_param_idx_in[pos];
    }
    ubrd_provided = true;
  }
  if (impr_i_atoms_in.size() == nimpr_in && impr_j_atoms_in.size() == nimpr_in &&
      impr_k_atoms_in.size() == nimpr_in && impr_l_atoms_in.size() == nimpr_in &&
      impr_param_idx_in.size() == nimpr_in) {
    for (int pos = 0; pos < nimpr_in; pos++) {
      impr_i_atoms[pos] = impr_i_atoms_in[pos]; 
      impr_j_atoms[pos] = impr_j_atoms_in[pos];
      impr_k_atoms[pos] = impr_k_atoms_in[pos];
      impr_l_atoms[pos] = impr_l_atoms_in[pos];
      impr_param_idx[pos] = impr_param_idx_in[pos];
    }
    impr_provided = true;
  }
  if (cmap_i_atoms_in.size() == ncmap_in && cmap_j_atoms_in.size() == ncmap_in &&
      cmap_k_atoms_in.size() == ncmap_in && cmap_l_atoms_in.size() == ncmap_in &&
      cmap_m_atoms_in.size() == ncmap_in && cmap_param_idx_in.size() == ncmap_in) {
    for (int pos = 0; pos < ncmap_in; pos++) {
      cmap_i_atoms[pos] = cmap_i_atoms_in[pos]; 
      cmap_j_atoms[pos] = cmap_j_atoms_in[pos];
      cmap_k_atoms[pos] = cmap_k_atoms_in[pos];
      cmap_l_atoms[pos] = cmap_l_atoms_in[pos];
      cmap_m_atoms[pos] = cmap_m_atoms_in[pos];
      cmap_param_idx[pos] = cmap_param_idx_in[pos];
    }
    cmap_provided = true;
  }
  if (ubrd_provided && impr_provided && cmap_provided) {
    makeAtomAssignments();
  }
}

//-------------------------------------------------------------------------------------------------
void CharmmValenceTable::makeAtomAssignments() {

  // Urey Bradley angles go to the first atom in the pair, CHARMM harmonic impropers and CMAP
  // terms go to the third atom in the quartet / quintet.
  for (int pos = 0; pos < total_ub_angles; pos++) {
    ubrd_assigned_bounds[ubrd_i_atoms[pos]] += 1;
  }
  for (int pos = 0; pos < total_impropers; pos++) {
    impr_assigned_bounds[impr_k_atoms[pos]] += 1;
  }
  for (int pos = 0; pos < total_cmaps; pos++) {
    cmap_assigned_bounds[cmap_k_atoms[pos]] += 1;
  }
  const PrefixSumType pfx_excl = PrefixSumType::EXCLUSIVE;
  prefixSumInPlace(&ubrd_assigned_bounds, pfx_excl, "charmmValenceIndexing");
  prefixSumInPlace(&impr_assigned_bounds, pfx_excl, "charmmValenceIndexing");
  prefixSumInPlace(&cmap_assigned_bounds, pfx_excl, "charmmValenceIndexing");

  // Populate the CHARMM valence assignments.
  for (int pos = 0; pos < total_ub_angles; pos++) {
    const int control_atom = ubrd_i_atoms[pos];
    const int fill_slot = ubrd_assigned_bounds[control_atom];
    ubrd_assigned_atoms[fill_slot]     = ubrd_k_atoms[pos];
    ubrd_assigned_index[fill_slot]     = ubrd_param_idx[pos];
    ubrd_assigned_terms[fill_slot]     = pos;
    ubrd_assigned_bounds[control_atom] = fill_slot + 1;
  }
  for (int pos = 0; pos < total_impropers; pos++) {
    const int control_atom = impr_k_atoms[pos];
    const int fill_slot = impr_assigned_bounds[control_atom];
    impr_assigned_atoms[ 3 * fill_slot     ]  = impr_i_atoms[pos];
    impr_assigned_atoms[(3 * fill_slot) + 1]  = impr_k_atoms[pos];
    impr_assigned_atoms[(3 * fill_slot) + 2]  = impr_l_atoms[pos];
    impr_assigned_index[fill_slot]            = impr_param_idx[pos];
    impr_assigned_terms[fill_slot]            = pos;
    impr_assigned_bounds[control_atom]        = fill_slot + 1;
  }
  for (int pos = 0; pos < total_cmaps; pos++) {
    const int control_atom = cmap_k_atoms[pos];
    const int fill_slot = cmap_assigned_bounds[control_atom];
    cmap_assigned_atoms[ 4 * fill_slot     ]  = cmap_i_atoms[pos];
    cmap_assigned_atoms[(4 * fill_slot) + 1]  = cmap_k_atoms[pos];
    cmap_assigned_atoms[(4 * fill_slot) + 2]  = cmap_l_atoms[pos];
    cmap_assigned_atoms[(4 * fill_slot) + 3]  = cmap_m_atoms[pos];
    cmap_assigned_index[fill_slot]            = cmap_param_idx[pos];
    cmap_assigned_terms[fill_slot]            = pos;
    cmap_assigned_bounds[control_atom]        = fill_slot + 1;
  }

  // Readjust the bounds arrays (trim the prefix sums)
  const int atom_count = static_cast<int>(ubrd_assigned_bounds.size()) - 1;
  for (int pos = atom_count; pos > 0; pos--) {
    ubrd_assigned_bounds[pos] = ubrd_assigned_bounds[pos - 1];
    impr_assigned_bounds[pos] = impr_assigned_bounds[pos - 1];
    cmap_assigned_bounds[pos] = cmap_assigned_bounds[pos - 1];
  }
  ubrd_assigned_bounds[0] = 0;
  impr_assigned_bounds[0] = 0;
  cmap_assigned_bounds[0] = 0;
}

//-------------------------------------------------------------------------------------------------
AttenuationParameterSet::AttenuationParameterSet(const int set_count_in) :
    total_14_sets{set_count_in},
    dihe14_parameter_indices{}, elec_screening_factors{}, vdw_screening_factors{}
{
  dihe14_parameter_indices.resize(set_count_in);
  elec_screening_factors.resize(set_count_in);
  vdw_screening_factors.resize(set_count_in);
}

//-------------------------------------------------------------------------------------------------
CondensedExclusions::CondensedExclusions() :
    total_exclusions{0}, atom_excl_bounds{}, atom_excl_list{}
{}

//-------------------------------------------------------------------------------------------------
CondensedExclusions::CondensedExclusions(const int natom_in, const int total_exclusions_in) :
    total_exclusions{total_exclusions_in},
    atom_excl_bounds{}, atom_excl_list{}
{
  atom_excl_bounds.resize(natom_in + 1, 0);
  atom_excl_list.resize(total_exclusions);
}

//-------------------------------------------------------------------------------------------------
VirtualSiteTable::VirtualSiteTable() :
    vs_count{0}, vs_numbers{}, vs_atoms{}, frame_types{}, frame1_atoms{}, frame2_atoms{},
    frame3_atoms{}, frame4_atoms{}, param_idx{}, frame_dim1{}, frame_dim2{}, frame_dim3{}  
{}

//-------------------------------------------------------------------------------------------------
VirtualSiteTable::VirtualSiteTable(const int natom_in, const int vs_count_in) :
    VirtualSiteTable()
{
  vs_count = vs_count_in;
  vs_numbers.resize(natom_in);
  vs_atoms.resize(vs_count_in);
  frame_types.resize(vs_count_in);
  frame1_atoms.resize(vs_count_in);
  frame2_atoms.resize(vs_count_in);
  frame3_atoms.resize(vs_count_in);
  frame4_atoms.resize(vs_count_in);
  param_idx.resize(vs_count_in);
  frame_dim1.resize(vs_count_in);
  frame_dim2.resize(vs_count_in);
  frame_dim3.resize(vs_count_in);
}

//-------------------------------------------------------------------------------------------------
Map1234::Map1234() :
    nb11_excl_bounds{}, nb11_excl_list{}, nb12_excl_bounds{}, nb12_excl_list{}, nb13_excl_bounds{},
    nb13_excl_list{}, nb14_excl_bounds{}, nb14_excl_list{}
{}

//-------------------------------------------------------------------------------------------------
Map1234::Map1234(const int natom_in, const int nb11_count_in, const int nb12_count_in,
                 const int nb13_count_in, const int nb14_count_in) :
    Map1234()
{
  nb11_excl_bounds.resize(natom_in + 1, 0);
  nb12_excl_bounds.resize(natom_in + 1, 0);
  nb13_excl_bounds.resize(natom_in + 1, 0);
  nb14_excl_bounds.resize(natom_in + 1, 0);
  nb11_excl_list.resize(nb11_count_in);
  nb12_excl_list.resize(nb12_count_in);
  nb13_excl_list.resize(nb13_count_in);
  nb14_excl_list.resize(nb14_count_in);
}

//-------------------------------------------------------------------------------------------------
CmapAccessories::CmapAccessories() :
    phi_derivatives{}, psi_derivatives{}, phi_psi_derivatives{}, patch_matrix_bounds{},
    patch_matrix_form{}
{}

//-------------------------------------------------------------------------------------------------
CmapAccessories::CmapAccessories(const std::vector<int> &cmap_dimensions_in) :
    CmapAccessories()
{
  const int nmaps = cmap_dimensions_in.size();
  int nelem = 0;
  for (int i = 0; i < nmaps; i++) {
    nelem += cmap_dimensions_in[i] * cmap_dimensions_in[i];
  }
  phi_derivatives.resize(nelem);
  psi_derivatives.resize(nelem);
  phi_psi_derivatives.resize(nelem);
  patch_matrix_bounds.resize(nmaps + 1);
  patch_matrix_form.resize(nelem * 16);
}

//-------------------------------------------------------------------------------------------------
ConstraintTable::ConstraintTable(const std::vector<int> &atomic_numbers,
                                 const std::vector<double> &atomic_masses,
                                 const std::vector<int> &mol_limits,
                                 const std::vector<int> &mol_contents,
                                 const std::vector<int> &mol_home,
                                 const BasicValenceTable &bvt, const Map1234 &all_nb_excl,
                                 const std::vector<double> &bond_equilibria,
                                 const std::vector<double> &angl_equilibria) :
    cnst_group_count{0}, settle_group_count{0}, cnst_parameter_count{0},
    settle_parameter_count{0}, cnst_group_list{}, cnst_group_bounds{}, cnst_group_param_idx{},
    settle_ox_atoms{}, settle_h1_atoms{}, settle_h2_atoms{}, settle_param_idx{},
    cnst_parameter_list{}, cnst_parameter_bounds{}, settle_measurements{}
{
  const int natom = atomic_numbers.size();
  const int nmol = static_cast<int>(mol_limits.size()) - 1;
  
  // Loop over all molecules, identifying suitable candidates for SETTLE groups
  std::vector<bool> settle_ready_atoms(natom, false);
  std::vector<bool> settle_ready_mols(nmol, false);
  std::vector<int2> real_populations(3);
  for (int i = 0; i < nmol; i++) {
    int n_real_atoms = 0;
    int n_distinct_z = 0;
    for (int j = mol_limits[i]; j < mol_limits[i + 1]; j++) {
      const int j_znum = atomic_numbers[mol_contents[j]];
      if (j_znum > 0) {
        n_real_atoms++;
        if (n_real_atoms <= 3) {
          bool identity_found = false;
          for (int k = 0; k < n_distinct_z; k++) {
            if (real_populations[k].x == j_znum) {
              real_populations[k].y += 1;
              identity_found = true;
            }
          }
          if (identity_found == false) {
            real_populations[n_distinct_z].x = j_znum;
            real_populations[n_distinct_z].y = 1;
            n_distinct_z++;
          }
        }
      }
    }

    // This is a SETTLE-capable group if there are three real atoms and two distinct Z numbers.
    // Otherwise, mark this molecule as not subject to SETTLE.
    if (n_real_atoms == 3 && n_distinct_z == 2) {
      settle_ready_mols[i] = true;
      for (int j = mol_limits[i]; j < mol_limits[i + 1]; j++) {
        settle_ready_atoms[mol_contents[j]] = true;
      }
    }
    else {
      for (int j = mol_limits[i]; j < mol_limits[i + 1]; j++) {
        settle_ready_atoms[mol_contents[j]] = false;
      }
    }
  }

  // Loop over all bonds, identifying those that bond to hydrogen and within SETTLE candidates.
  // The differential_bonds vector stores, for every molecule, the number of bonds between two
  // atoms of different atomic numbers.  The sum is only relevant for SETTLE-capable molecules as
  // identified in the settle_ready_mols vector, but if settle_ready_mols[k] is TRUE and
  // differential_bonds[k] == 2, then molecule k is suitable for SETTLE.
  std::vector<int> differential_bonds(nmol, 0);
  for (int pos = 0; pos < bvt.total_bonds; pos++) {
    const int iatom = bvt.bond_i_atoms[pos];
    const int jatom = bvt.bond_j_atoms[pos];
    if (settle_ready_atoms[iatom] && settle_ready_atoms[jatom]) {
      if (mol_home[iatom] != mol_home[jatom]) {
        rtErr("Bond " + std::to_string(pos) + " between atoms " + std::to_string(iatom) + " and " +
              std::to_string(jatom) + " was detected to span two molecules, " +
              std::to_string(mol_home[iatom]) + " and " + std::to_string(mol_home[jatom]) + ".",
              "ConstraintTable");
      }
      differential_bonds[mol_home[iatom]] += (atomic_numbers[iatom] != atomic_numbers[jatom]);
    }
    else if (settle_ready_atoms[iatom] && settle_ready_atoms[jatom] == false) {
      rtErr("Inconsistent SETTLE readiness was detected in two bonded atoms.  This should not be "
            "possible, as SETTLE readiness is applied to an antire molecule and bonds only exist "
            "within a molecule.", "ConstraintTable");
    }
  }
  int nsett = 0;
  for (int i = 0; i < nmol; i++) {
    nsett += (settle_ready_mols[i] && differential_bonds[i] == 2);
  }
  settle_group_count = nsett;
  settle_ox_atoms.resize(nsett);
  settle_h1_atoms.resize(nsett);
  settle_h2_atoms.resize(nsett);
  settle_param_idx.resize(nsett);
  nsett = 0;
  for (int i = 0; i < nmol; i++) {
    if (settle_ready_mols[i] && differential_bonds[i] == 2) {
      int real_atom_idx[3];
      int real_atom_znum[3];
      int n_real_found = 0;
      for (int j = mol_limits[i]; j < mol_limits[i + 1]; j++) {
        if (atomic_numbers[mol_contents[j]] > 0) {
          real_atom_idx[n_real_found] = mol_contents[j];
          real_atom_znum[n_real_found] = atomic_numbers[mol_contents[j]];
          n_real_found++;
        }
      }
      if (n_real_found != 3) {
        rtErr("Expected 3 real atoms in molecule " + std::to_string(i) + ", found " +
              std::to_string(n_real_found) + ".", "ConstraintTable");
      }

      // Find the pair of "hydrogen" atoms as a guide
      if (real_atom_znum[0] == real_atom_znum[1]) {
        settle_ox_atoms[nsett] = real_atom_idx[2];
        settle_h1_atoms[nsett] = std::min(real_atom_idx[0], real_atom_idx[1]);
        settle_h2_atoms[nsett] = std::max(real_atom_idx[0], real_atom_idx[1]);
      }
      else if (real_atom_znum[0] == real_atom_znum[2]) {
        settle_ox_atoms[nsett] = real_atom_idx[1];
        settle_h1_atoms[nsett] = std::min(real_atom_idx[0], real_atom_idx[2]);
        settle_h2_atoms[nsett] = std::max(real_atom_idx[0], real_atom_idx[2]);        
      }
      else if (real_atom_znum[1] == real_atom_znum[2]) {
        settle_ox_atoms[nsett] = real_atom_idx[0];
        settle_h1_atoms[nsett] = std::min(real_atom_idx[1], real_atom_idx[2]);
        settle_h2_atoms[nsett] = std::max(real_atom_idx[1], real_atom_idx[2]);        
      }
      else {
        rtErr("Expected two atoms of identical mass and one of a different mass among " +
              std::to_string(real_atom_idx[0]) + ", " + std::to_string(real_atom_idx[1]) +
              ", and " + std::to_string(real_atom_idx[2]) + ".", "ConstraintTable");
      }
      const SettleParm stt_pre = getSettleParameters(settle_ox_atoms[nsett],
                                                     settle_h1_atoms[nsett],
                                                     settle_h2_atoms[nsett], atomic_masses, bvt,
                                                     bond_equilibria, angl_equilibria);
      bool parmset_found = false;
      for (int j = 0; j < settle_parameter_count; j++) {
        if (fabs(settle_measurements[j].mormt - stt_pre.mormt) < constants::tiny &&
            fabs(settle_measurements[j].mhrmt - stt_pre.mhrmt) < constants::tiny &&
            fabs(settle_measurements[j].ra - stt_pre.ra) < constants::tiny &&
            fabs(settle_measurements[j].rc - stt_pre.rc) < constants::tiny) {
          settle_param_idx[nsett] = j;
          parmset_found = true;
        }
      }
      if (parmset_found == false) {
        settle_param_idx[nsett] = settle_parameter_count;
        settle_measurements.push_back(stt_pre);
        settle_parameter_count += 1;
      }
      nsett++;
    }
  }

  // Loop over all molecules, identifying constraint groups with one or more hydrogen atoms bound
  // to a heavier atom.  Expend the settle_ready_atoms array to record other hydrogens as they
  // are incorporated into constraint groups.
  std::vector<int> tmp_cnst_list(8);
  std::vector<double2> tmp_cnst_parm(8);
  std::vector<bool> tmp_cnst_made(8);
  cnst_parameter_bounds.push_back(0);
  cnst_group_bounds.push_back(0);
  int n_available_hydrogen = 0;
  int ngrp = 0;
  for (int i = 0; i < natom; i++) {
    n_available_hydrogen += (settle_ready_atoms[i] == false && atomic_numbers[i] == 1);
  }
  cnst_group_list.reserve(2 * n_available_hydrogen);
  cnst_group_param_idx.reserve(n_available_hydrogen);
  for (int i = 0; i < natom; i++) {
    if (settle_ready_atoms[i] || atomic_numbers[i] != 1) {
      continue;
    }

    // This is a hydrogen with no constraint or SETTLE group.  Find its parent atom.
    settle_ready_atoms[i] = true;
    int n_heavy = 0;
    int heavy_atom_idx;
    for (int j = all_nb_excl.nb12_excl_bounds[i]; j < all_nb_excl.nb12_excl_bounds[i + 1]; j++) {
      if (atomic_numbers[all_nb_excl.nb12_excl_list[j]] >= 1) {
        n_heavy++;
        heavy_atom_idx = all_nb_excl.nb12_excl_list[j];
      }
    }
    
    // There is a possibility of simulations having free hydrogen atoms.  It would be strange,
    // but not against the rules of the simulation.  Just let it be.  Do not attempt to constrain
    // bonds to hydrogen if more than one heavy atom (or more than one other atom) is
    // participating, however.  That is for SETTLE, and in other situations the behavior will
    // remain undefined by this implementation.
    if (n_heavy != 1) {
      continue;
    }

    // Form the constraint group.
    tmp_cnst_list.resize(0);
    tmp_cnst_list.push_back(heavy_atom_idx);
    tmp_cnst_list.push_back(i);
    for (int j = all_nb_excl.nb12_excl_bounds[heavy_atom_idx];
         j < all_nb_excl.nb12_excl_bounds[heavy_atom_idx + 1]; j++) {
      if (settle_ready_atoms[all_nb_excl.nb12_excl_list[j]] == false &&
          atomic_numbers[all_nb_excl.nb12_excl_list[j]] == 1) {
        tmp_cnst_list.push_back(all_nb_excl.nb12_excl_list[j]);
        settle_ready_atoms[all_nb_excl.nb12_excl_list[j]] = true;
      }
    }
    const int group_size = tmp_cnst_list.size();
    cnst_group_bounds.push_back(cnst_group_bounds[ngrp] + group_size);
    for (int j = 0; j < group_size; j++) {
      cnst_group_list.push_back(tmp_cnst_list[j]);
    }

    // Make the parameter set for the constraint group.
    tmp_cnst_made.resize(group_size);
    tmp_cnst_parm.resize(group_size);
    for (int j = 0; j < group_size; j++) {
      tmp_cnst_made[j] = false;
      tmp_cnst_parm[j].x = 1.0 / atomic_masses[tmp_cnst_list[j]];
    }
    tmp_cnst_parm[0].y = 0.0;
    tmp_cnst_made[0] = true;
    const int catom = tmp_cnst_list[0];
    for (int j = bvt.bond_assigned_bounds[catom]; j < bvt.bond_assigned_bounds[catom + 1]; j++) {
      const int other_atom = bvt.bond_assigned_atoms[j];
      if (atomic_numbers[other_atom] == 1) {
        for (int k = 1; k < group_size; k++) {
          if (other_atom == tmp_cnst_list[k]) {
            const double l_eq = bond_equilibria[bvt.bond_param_idx[bvt.bond_assigned_terms[j]]];
            tmp_cnst_parm[k].y = l_eq * l_eq;
            tmp_cnst_made[k] = true;
          }
        }
      }
    }
    for (int j = 1; j < group_size; j++) {
      if (tmp_cnst_made[j]) {
        continue;
      }
      const int jatom = tmp_cnst_list[j];
      for (int k = bvt.bond_assigned_bounds[jatom]; k < bvt.bond_assigned_bounds[jatom + 1]; k++) {
        if (bvt.bond_assigned_atoms[k] == catom) {
          const double l_eq = bond_equilibria[bvt.bond_param_idx[bvt.bond_assigned_terms[k]]];
          tmp_cnst_parm[j].y = l_eq * l_eq;
          tmp_cnst_made[j] = true;
        }
      }
    }

    // Check whether these parameters are unique
    bool parmset_found = false;
    for (int j = 0; j < cnst_parameter_count; j++) {
      if (cnst_parameter_bounds[j + 1] - cnst_parameter_bounds[j] != group_size) {
        continue;
      }
      const int joffset = cnst_parameter_bounds[j];
      bool match = true;
      for (int k = 0; k < group_size; k++) {
        const double xdiff = fabs(tmp_cnst_parm[k].x - cnst_parameter_list[joffset + k].x);
        const double ydiff = fabs(tmp_cnst_parm[k].y - cnst_parameter_list[joffset + k].y);
        match = (match && xdiff < constants::tiny && ydiff < constants::tiny);
      }
      if (match) {
        parmset_found = true;
        cnst_group_param_idx.push_back(j);
      }
    }
    if (parmset_found == false) {
      cnst_group_param_idx.push_back(cnst_parameter_count);
      cnst_parameter_bounds.push_back(cnst_parameter_bounds[cnst_parameter_count] + group_size);
      for (int j = 0; j < group_size; j++) {
        cnst_parameter_list.push_back(tmp_cnst_parm[j]);
      }
      cnst_parameter_count += 1;
    }

    // Increment the number of constraint groups
    ngrp++;
  }
  cnst_group_count = ngrp;
}

//-------------------------------------------------------------------------------------------------
void smoothCharges(std::vector<double> *q, std::vector<double> *tmp_charge_parameters,
                   std::vector<int> *tmp_charge_indices, int *q_param_count,
                   const double rounding_tol, const double charge_discretization,
                   const std::string &filename) {
  const double qsum = sum<double>(*q);
  if (rounding_tol < 0.0) {
    rtErr("A tolerance of " + realToString(rounding_tol, 15, 8, NumberFormat::SCIENTIFIC) +
          " is invalid for rounding total charge to integral values.", "smoothCharges");
  }
  if (charge_discretization > 0.001) {
    rtErr("A charge discretization of 0.001e or better is required for charge smoothing.",
          "smoothCharges");
  }
  if (charge_discretization < 0.0) {
    rtErr("A discretization of " +
          realToString(charge_discretization, 15, 8, NumberFormat::SCIENTIFIC) +
          " is invalid for discretizing partial atomic charges.", "smoothCharges");
  }
  double* qdata = q->data();
  const int natom = q->size();
  if (fabs(qsum - round(qsum)) < rounding_tol) {
    const double dq = (round(qsum) - qsum) / static_cast<double>(natom);
    for (int i = 0; i < natom; i++) {
      qdata[i] += dq;
    }

    // Return if the discretization is extremely small
    if (charge_discretization < constants::tiny) {
      return;
    }

    // There is no more work to do if the charges do not need to add to the target integral
    // value when summed.  Otherwise round all charges to units of the discretization and then
    // sprinkle increments of charge among the largest values.
    double disc_qsum = 0.0;
    for (int i = 0; i < natom; i++) {
      qdata[i] = llround(qdata[i] / charge_discretization);
      qdata[i] *= charge_discretization;
      disc_qsum += qdata[i];
    }
    const double discretized_dq = (round(disc_qsum) - disc_qsum > 0.0) ?  charge_discretization :
                                                                         -charge_discretization;
    const int discretized_dq_increments = round(fabs((round(disc_qsum) - disc_qsum) /
                                                     charge_discretization));

    // Check that no atomic partial charge goes over 21.4.  If it does, we just can't do the
    // final sprinkling of charge to make the sum of discretized charges come out to an exact
    // integral.
    if (maxValue(*q) < 21.4 && minValue(*q) > -21.4 && sizeof(int) >= 4) {
      if (discretized_dq_increments > 0) {

	// No atomic partial charge goes over 21.4 and the standard (long) int is big enough.
        // The discretization is possible.
        int2 tmpi = { 0, 0 };
        std::vector<int2> sortable_q(natom, tmpi);
        for (int i = 0; i < natom; i++) {
          tmpi = { static_cast<int>(qdata[i] * 100000000), i };
          sortable_q[i] = tmpi;
        }
        std::sort(sortable_q.begin(), sortable_q.end(),
                  [](int2 a, int2 b) { return abs(a.x) > abs(b.x); });

        // Within the list, there can be many partial charges with the same values (to within
        // the specified discretization).  Because of the way different machines / compilers
        // execute std::sort when there are ties, this can create different orderings of atoms and
        // thereby sprinkle the charge increments needed to neutralize the system onto different
        // atoms.  To prevent this behavior, step back through the list and make sure to order any
        // subsets that share the same partial charge by increasing order of their atom ID.
        int seci = 0;
        while (seci < natom) {
          int secj = seci + 1;
          while (secj < natom && fabs(qdata[sortable_q[seci].y] -
                                      qdata[sortable_q[secj].y]) < charge_discretization) {
            secj++;
          }
          std::sort(sortable_q.begin() + seci, sortable_q.begin() + secj,
                    [](int2 a, int2 b) { return abs(a.y) < abs(b.y); });
          seci = secj;
        }

        // Add discrete charge increments as necessary to get to a neutral system, or an integer
        // net charge.
        int j = 0;
        for (int i = 0; i < discretized_dq_increments; i++) {
          qdata[sortable_q[j].y] += discretized_dq;
          j++;
          if (j == natom) {
            j = 0;
          }
        }
      }
    }
    else {
      const std::string opener = (filename.size() > 0) ? "Topology " + filename : "Charge set ";
      rtWarn(opener + " contains atomic partial charges which are too large for the " +
             std::to_string(sizeof(int)) + "-byte integer format.  No partial charge may be "
             "larger than 1 / 100,000,000th of the largest representable positive or negative "
             "integers.", "smoothCharges");
    }
  }

  // Create an additional array of charges based on their discretized values, assigning each an
  // index.
  std::vector<double> unique_q(natom, 0);
  std::vector<llint> discretized_charges(natom, 0);
  for (int i = 0; i < natom; i++) {
    discretized_charges[i] = llround(qdata[i] / charge_discretization);
  }
  std::vector<bool> found(natom, false);
  int* tq_idx_ptr = tmp_charge_indices->data();
  double* tq_param_ptr = tmp_charge_parameters->data();
  int n_unique = 0;
  for (int i = 0; i < natom; i++) {
    if (found[i]) {
      continue;
    }
    tq_param_ptr[n_unique] = qdata[i];
    for (int j = i; j < natom; j++) {
      if (found[j] == false && discretized_charges[j] == discretized_charges[i]) {
        found[j] = true;
        tq_idx_ptr[j] = n_unique;
      }
    }
    n_unique++;
  }

  // Relay the number of unique charges back to the calling function and also shrink the array
  // of these parameters to fit that size (this will be critical when uploading the array to the
  // topology's Hybrid objects later, as space is cut precisely based on constants such as the
  // number of atoms or certain types of parameters).
  *q_param_count = n_unique;
  tmp_charge_parameters->resize(n_unique);
}

//-------------------------------------------------------------------------------------------------
void expandLennardJonesTables(std::vector<double> *lj_a_values, std::vector<double> *lj_b_values,
                              std::vector<double> *lj_c_values,
                              std::vector<double> *lj_14_a_values,
                              std::vector<double> *lj_14_b_values,
                              std::vector<double> *lj_14_c_values, 
                              std::vector<double> *hb_a_values, std::vector<double> *hb_b_values,
                              int n_lj_types, const std::vector<int> &nb_param_index) {

  // Resize the Lennard-Jones arrays, then get data pointers
  const int n2_types = n_lj_types * n_lj_types;
  lj_a_values->resize(n2_types);
  lj_b_values->resize(n2_types);
  lj_c_values->resize(n2_types, 0.0);
  lj_14_a_values->resize(n2_types);
  lj_14_b_values->resize(n2_types);
  lj_14_c_values->resize(n2_types, 0.0);
  hb_a_values->resize(n2_types, 0.0);
  hb_b_values->resize(n2_types, 0.0);
  double* lja_ptr = lj_a_values->data();
  double* ljb_ptr = lj_b_values->data();
  double* ljc_ptr = lj_c_values->data();
  double* lja_14_ptr = lj_14_a_values->data();
  double* ljb_14_ptr = lj_14_b_values->data();
  double* ljc_14_ptr = lj_14_c_values->data();
  double* hba_ptr = hb_a_values->data();
  double* hbb_ptr = hb_a_values->data();
  std::vector<double> lj_a_tmp(n2_types);
  std::vector<double> lj_b_tmp(n2_types);
  std::vector<double> lj_c_tmp(n2_types);
  std::vector<double> lj_14_a_tmp(n2_types);
  std::vector<double> lj_14_b_tmp(n2_types);
  std::vector<double> lj_14_c_tmp(n2_types);
  std::vector<double> hb_a_tmp(n2_types);
  std::vector<double> hb_b_tmp(n2_types);
  for (int i = 0; i < n_lj_types; i++) {
    for (int j = 0; j <= i; j++) {
      const int nb_lkp = nb_param_index[(n_lj_types * i) + j];
      if (nb_lkp >= 0) {
        lj_a_tmp[(n_lj_types * i) + j] = lja_ptr[nb_lkp];
        lj_a_tmp[(n_lj_types * j) + i] = lja_ptr[nb_lkp];
        lj_b_tmp[(n_lj_types * i) + j] = ljb_ptr[nb_lkp];
        lj_b_tmp[(n_lj_types * j) + i] = ljb_ptr[nb_lkp];
        lj_c_tmp[(n_lj_types * i) + j] = ljc_ptr[nb_lkp];
        lj_c_tmp[(n_lj_types * j) + i] = ljc_ptr[nb_lkp];
        lj_14_a_tmp[(n_lj_types * i) + j] = lja_14_ptr[nb_lkp];
        lj_14_a_tmp[(n_lj_types * j) + i] = lja_14_ptr[nb_lkp];
        lj_14_b_tmp[(n_lj_types * i) + j] = ljb_14_ptr[nb_lkp];
        lj_14_b_tmp[(n_lj_types * j) + i] = ljb_14_ptr[nb_lkp];
        lj_14_c_tmp[(n_lj_types * i) + j] = ljc_14_ptr[nb_lkp];
        lj_14_c_tmp[(n_lj_types * j) + i] = ljc_14_ptr[nb_lkp];
      }
      else {
        hb_a_tmp[(n_lj_types * i) + j] = hba_ptr[-nb_lkp];
        hb_a_tmp[(n_lj_types * j) + i] = hba_ptr[-nb_lkp];
        hb_b_tmp[(n_lj_types * i) + j] = hbb_ptr[-nb_lkp];
        hb_b_tmp[(n_lj_types * j) + i] = hbb_ptr[-nb_lkp];
      }
    }
  }

  // Refill the Lennard-Jones arrays
  for (int i = 0; i < n2_types; i++) {
    lja_ptr[i] = lj_a_tmp[i];
    ljb_ptr[i] = lj_b_tmp[i];
    ljc_ptr[i] = lj_c_tmp[i];
    lja_14_ptr[i] = lj_14_a_tmp[i];
    ljb_14_ptr[i] = lj_14_b_tmp[i];
    ljc_14_ptr[i] = lj_14_c_tmp[i];
    hba_ptr[i] = hb_a_tmp[i];
    hbb_ptr[i] = hb_b_tmp[i];
  }
}

//-------------------------------------------------------------------------------------------------
CondensedExclusions processExclusions(const std::vector<int> &raw_counts,
                                      const std::vector<int> &raw_exclusions,
                                      const std::string &file_name) {
  CondensedExclusions ce;

  // Check through the counts, making a preliminary prefix sum.  If there is only one exclusion
  // and that atom is 0, count zero actual exclusions.  Also adjust the Fortran indexing to
  // C / C++ indexing.
  const int natom = raw_counts.size();
  std::vector<int> prefix(natom + 1, 0);
  int raw_sum = 0;
  int p_sum = 0;
  for (int i = 0; i < natom; i++) {
    prefix[i] = p_sum;
    p_sum = (raw_counts[i] == 1 &&
             raw_exclusions[raw_sum] == 0) ? p_sum : p_sum + raw_counts[i];

    // Check that there are no zero entries if an atom has more than one exclusions.
    // Check to confirm that no atom names the same excluded atom twice.
    if (raw_counts[i] > 1) {
      const int jlim = raw_sum + raw_counts[i];
      for (int j = raw_sum; j < jlim; j++) {
        if (raw_exclusions[j] == 0) {
          rtErr("Error parsing exclusions list for atom " + std::to_string(i) + " in topology " +
                file_name + ".  An atom with " + std::to_string(raw_counts[i]) +
                " exclusions cannot exclude a blank atom.", "processExclusions");
        }
        for (int k = j + 1; k < jlim; k++) {
          if (raw_exclusions[j] == raw_exclusions[k]) {
            rtErr("Error parsing exclusions for atom " + std::to_string(i) + " in topology " +
                  file_name + ".  The atom excludes atom " + std::to_string(raw_exclusions[j]) +
                  " at least twice.", "processExclusions");
          }
        }
      }
    }
    else if (raw_counts[i] == 0) {
      rtErr("Error parsing exclusions for atom " + std::to_string(i) + " in topology " +
            file_name + ".  The topology cannot explicitly state there are no exclusions.",
            "processExclusions");
    }

    // Increment the raw sum
    raw_sum += raw_counts[i];
  }
  prefix[natom] = p_sum;

  // Construct and return the result
  ce.total_exclusions = p_sum;
  ce.atom_excl_bounds = prefix;
  ce.atom_excl_list.resize(p_sum);
  int real_atoms = 0;
  for (int i = 0; i < raw_sum; i++) {
    if (raw_exclusions[i] > 0) {
      ce.atom_excl_list[real_atoms] = raw_exclusions[i] - 1;
      real_atoms++;
    }
    else if (raw_exclusions[i] < 0) {
      rtErr("Error parsing exclusions for topology " + file_name + ".  Invalid excluded atom " +
            std::to_string(raw_exclusions[i]) + ".", "processExclusions");
    }
  }
  return ce;
}

//-------------------------------------------------------------------------------------------------
CondensedExclusions calculatePrmtopExclusions(const Map1234 nb_excl) {

  // Allocate to hold the reduced number of exclusions.
  CondensedExclusions result;
  const int natom = nb_excl.nb11_excl_bounds.size() - 1;
  result.total_exclusions = (nb_excl.nb11_excl_bounds[natom] +
                             nb_excl.nb12_excl_bounds[natom] +
                             nb_excl.nb13_excl_bounds[natom] +
                             nb_excl.nb14_excl_bounds[natom]) / 2;
  result.atom_excl_list.resize(result.total_exclusions);
  result.atom_excl_bounds.resize(natom + 1, 0);

  // Set up a prefix sum over each atom's forward exclusion complement.  Fill out the exclusions
  // list as the prefix sum accumulates.
  int excl_idx = 0;
  for (int i = 0; i < natom; i++) {
    int ni_excl = 0;
    for (int j = nb_excl.nb11_excl_bounds[i]; j < nb_excl.nb11_excl_bounds[i + 1]; j++) {
      if (nb_excl.nb11_excl_list[j] > i) {
        result.atom_excl_list[excl_idx] = nb_excl.nb11_excl_list[j];
        excl_idx++;
        ni_excl++;
      }
    }
    for (int j = nb_excl.nb12_excl_bounds[i]; j < nb_excl.nb12_excl_bounds[i + 1]; j++) {
      if (nb_excl.nb12_excl_list[j] > i) {
        result.atom_excl_list[excl_idx] = nb_excl.nb12_excl_list[j];
        excl_idx++;
        ni_excl++;
      }
    }
    for (int j = nb_excl.nb13_excl_bounds[i]; j < nb_excl.nb13_excl_bounds[i + 1]; j++) {
      if (nb_excl.nb13_excl_list[j] > i) {
        result.atom_excl_list[excl_idx] = nb_excl.nb13_excl_list[j];
        excl_idx++;
        ni_excl++;
      }
    }
    for (int j = nb_excl.nb14_excl_bounds[i]; j < nb_excl.nb14_excl_bounds[i + 1]; j++) {
      if (nb_excl.nb14_excl_list[j] > i) {
        result.atom_excl_list[excl_idx] = nb_excl.nb14_excl_list[j];
        excl_idx++;
        ni_excl++;
      }
    }
    result.atom_excl_bounds[i] = ni_excl;
  }
  prefixSumInPlace(&result.atom_excl_bounds, PrefixSumType::EXCLUSIVE);
  return result;
}

//-------------------------------------------------------------------------------------------------
int countPrmtopZeroExclusions(const Map1234 nb_excl) {
  int result = 0;
  const int natom = nb_excl.nb11_excl_bounds.size() - 1;
  for (int i = 0; i < natom; i++) {
    bool fw_excl = false;
    for (int j = nb_excl.nb11_excl_bounds[i]; j < nb_excl.nb11_excl_bounds[i + 1]; j++) {
      fw_excl = (fw_excl || (nb_excl.nb11_excl_list[j] > i));
    }
    for (int j = nb_excl.nb12_excl_bounds[i]; j < nb_excl.nb12_excl_bounds[i + 1]; j++) {
      fw_excl = (fw_excl || (nb_excl.nb12_excl_list[j] > i));
    }
    for (int j = nb_excl.nb13_excl_bounds[i]; j < nb_excl.nb13_excl_bounds[i + 1]; j++) {
      fw_excl = (fw_excl || (nb_excl.nb13_excl_list[j] > i));
    }
    for (int j = nb_excl.nb14_excl_bounds[i]; j < nb_excl.nb14_excl_bounds[i + 1]; j++) {
      fw_excl = (fw_excl || (nb_excl.nb14_excl_list[j] > i));
    }
    result += (fw_excl == false);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
BasicValenceTable basicValenceIndexing(const int atom_count,
                                       const std::vector<int> &tmp_bond_atoms_h,
                                       const std::vector<int> &tmp_bond_atoms_noh,
                                       const std::vector<int> &tmp_angl_atoms_h,
                                       const std::vector<int> &tmp_angl_atoms_noh,
                                       const std::vector<int> &tmp_dihe_atoms_h,
                                       const std::vector<int> &tmp_dihe_atoms_noh) {

  // Resize the arrays to accommdate the incoming data
  const int n_bonds_h   = tmp_bond_atoms_h.size() / 3;
  const int n_bonds_noh = tmp_bond_atoms_noh.size() / 3;
  const int n_angls_h   = tmp_angl_atoms_h.size() / 4;
  const int n_angls_noh = tmp_angl_atoms_noh.size() / 4;
  const int n_dihes_h   = tmp_dihe_atoms_h.size() / 5;
  const int n_dihes_noh = tmp_dihe_atoms_noh.size() / 5;
  BasicValenceTable bvt(atom_count, n_bonds_h + n_bonds_noh, n_angls_h + n_angls_noh,
                        n_dihes_h + n_dihes_noh);

  // Pluck the indices and parameters, adjusting everything for C / C++ along the way
  for (int pos = 0; pos < n_bonds_h; pos++) {
    bvt.bond_i_atoms[pos]   = tmp_bond_atoms_h[3*pos    ] / 3;
    bvt.bond_j_atoms[pos]   = tmp_bond_atoms_h[3*pos + 1] / 3;
    bvt.bond_param_idx[pos] = tmp_bond_atoms_h[3*pos + 2] - 1;
    bvt.bond_mods[pos].y    = static_cast<char>(HydrogenContent::HAS_HYDROGEN);
  }
  for (int pos = 0; pos < n_bonds_noh; pos++) {
    bvt.bond_i_atoms[pos + n_bonds_h]   = tmp_bond_atoms_noh[3*pos    ] / 3;
    bvt.bond_j_atoms[pos + n_bonds_h]   = tmp_bond_atoms_noh[3*pos + 1] / 3;
    bvt.bond_param_idx[pos + n_bonds_h] = tmp_bond_atoms_noh[3*pos + 2] - 1;
    bvt.bond_mods[pos + n_bonds_h].y    = static_cast<char>(HydrogenContent::NO_HYDROGEN);
  }
  for (int pos = 0; pos < n_angls_h; pos++) {
    bvt.angl_i_atoms[pos]   = tmp_angl_atoms_h[4*pos    ] / 3;
    bvt.angl_j_atoms[pos]   = tmp_angl_atoms_h[4*pos + 1] / 3;
    bvt.angl_k_atoms[pos]   = tmp_angl_atoms_h[4*pos + 2] / 3;
    bvt.angl_param_idx[pos] = tmp_angl_atoms_h[4*pos + 3] - 1;
    bvt.angl_mods[pos].y    = static_cast<char>(HydrogenContent::HAS_HYDROGEN);
  }
  for (int pos = 0; pos < n_angls_noh; pos++) {
    bvt.angl_i_atoms[pos + n_angls_h]   = tmp_angl_atoms_noh[4*pos    ] / 3;
    bvt.angl_j_atoms[pos + n_angls_h]   = tmp_angl_atoms_noh[4*pos + 1] / 3;
    bvt.angl_k_atoms[pos + n_angls_h]   = tmp_angl_atoms_noh[4*pos + 2] / 3;
    bvt.angl_param_idx[pos + n_angls_h] = tmp_angl_atoms_noh[4*pos + 3] - 1;
    bvt.angl_mods[pos + n_angls_h].y    = static_cast<char>(HydrogenContent::NO_HYDROGEN);
  }
  for (int pos = 0; pos < n_dihes_h; pos++) {
    bvt.dihe_i_atoms[pos]   = tmp_dihe_atoms_h[5*pos    ] / 3;
    bvt.dihe_j_atoms[pos]   = tmp_dihe_atoms_h[5*pos + 1] / 3;
    bvt.dihe_k_atoms[pos]   = tmp_dihe_atoms_h[5*pos + 2] / 3;
    bvt.dihe_l_atoms[pos]   = tmp_dihe_atoms_h[5*pos + 3] / 3;
    bvt.dihe_param_idx[pos] = tmp_dihe_atoms_h[5*pos + 4] - 1;
    bvt.dihe_mods[pos].y    = static_cast<char>(HydrogenContent::HAS_HYDROGEN);
  }
  for (int pos = 0; pos < n_dihes_noh; pos++) {
    bvt.dihe_i_atoms[pos + n_dihes_h]   = tmp_dihe_atoms_noh[5*pos    ] / 3;
    bvt.dihe_j_atoms[pos + n_dihes_h]   = tmp_dihe_atoms_noh[5*pos + 1] / 3;
    bvt.dihe_k_atoms[pos + n_dihes_h]   = tmp_dihe_atoms_noh[5*pos + 2] / 3;
    bvt.dihe_l_atoms[pos + n_dihes_h]   = tmp_dihe_atoms_noh[5*pos + 3] / 3;
    bvt.dihe_param_idx[pos + n_dihes_h] = tmp_dihe_atoms_noh[5*pos + 4] - 1;
    bvt.dihe_mods[pos + n_dihes_h].y    = static_cast<char>(HydrogenContent::NO_HYDROGEN);
  }

  // Mark additional modifiers for each term.  Constraints cannot be assigned without dynamics
  // input, so mark all terms as free.
  for (int pos = 0; pos < bvt.total_bonds; pos++) {
    bvt.bond_mods[pos].x = static_cast<char>(ConstraintStatus::FREE);
    bvt.bond_mods[pos].z = static_cast<char>(ForceFieldFamily::BASIC);
  }
  for (int pos = 0; pos < bvt.total_angls; pos++) {
    bvt.angl_mods[pos].x = static_cast<char>(ConstraintStatus::FREE);
    bvt.angl_mods[pos].z = static_cast<char>(ForceFieldFamily::BASIC);
  }
  for (int pos = 0; pos < bvt.total_dihes; pos++) {
    bvt.dihe_mods[pos].x = static_cast<char>(ConstraintStatus::FREE);
    bvt.dihe_mods[pos].z = static_cast<char>(ForceFieldFamily::BASIC);
    if (bvt.dihe_k_atoms[pos] < 0 && bvt.dihe_l_atoms[pos] >= 0) {
      bvt.dihe_k_atoms[pos] *= -1;
      bvt.dihe_mods[pos].w = static_cast<char>(TorsionKind::PROPER_NO_14);
    }
    else if (bvt.dihe_k_atoms[pos] < 0 && bvt.dihe_l_atoms[pos] < 0) {
      bvt.dihe_k_atoms[pos] *= -1;
      bvt.dihe_l_atoms[pos] *= -1;
      bvt.dihe_mods[pos].w = static_cast<char>(TorsionKind::IMPROPER_NO_14);
    }
    else if (bvt.dihe_k_atoms[pos] >= 0 && bvt.dihe_l_atoms[pos] < 0) {
      bvt.dihe_l_atoms[pos] *= -1;
      bvt.dihe_mods[pos].w = static_cast<char>(TorsionKind::IMPROPER);
    }
    else {
      bvt.dihe_mods[pos].w = static_cast<char>(TorsionKind::PROPER);
    }
  }

  // Assign terms to each atom
  bvt.makeAtomAssignments();
  
  return bvt;
}

//-------------------------------------------------------------------------------------------------
CharmmValenceTable charmmValenceIndexing(const int atom_count,
                                         const std::vector<int> &tmp_ub_atoms,
                                         const std::vector<int> &tmp_charmm_impr_atoms,
                                         const std::vector<int> &tmp_cmap_atoms) {
  CharmmValenceTable mvt(atom_count, tmp_ub_atoms.size() / 3, tmp_charmm_impr_atoms.size() / 5,
                         tmp_cmap_atoms.size() / 6);

  // The indexing for CHARMM force field terms does not contain pre-computations like the basic
  // valence term indexing in an Amber topology.  All elements are Fortran-indexed.
  for (int pos = 0; pos < mvt.total_ub_angles; pos++) {
    mvt.ubrd_i_atoms[pos] =   tmp_ub_atoms[3*pos    ] - 1;
    mvt.ubrd_k_atoms[pos] =   tmp_ub_atoms[3*pos + 1] - 1;
    mvt.ubrd_param_idx[pos] = tmp_ub_atoms[3*pos + 2] - 1;
  }
  for (int pos = 0; pos < mvt.total_impropers; pos++) {
    mvt.impr_i_atoms[pos]   = tmp_charmm_impr_atoms[5*pos    ] - 1;
    mvt.impr_j_atoms[pos]   = tmp_charmm_impr_atoms[5*pos + 1] - 1;
    mvt.impr_k_atoms[pos]   = tmp_charmm_impr_atoms[5*pos + 2] - 1;
    mvt.impr_l_atoms[pos]   = tmp_charmm_impr_atoms[5*pos + 3] - 1;
    mvt.impr_param_idx[pos] = tmp_charmm_impr_atoms[5*pos + 4] - 1;
  }
  for (int pos = 0; pos < mvt.total_cmaps; pos++) {
    mvt.cmap_i_atoms[pos]   = tmp_cmap_atoms[6*pos    ] - 1;
    mvt.cmap_j_atoms[pos]   = tmp_cmap_atoms[6*pos + 1] - 1;
    mvt.cmap_k_atoms[pos]   = tmp_cmap_atoms[6*pos + 2] - 1;
    mvt.cmap_l_atoms[pos]   = tmp_cmap_atoms[6*pos + 3] - 1;
    mvt.cmap_m_atoms[pos]   = tmp_cmap_atoms[6*pos + 4] - 1;
    mvt.cmap_param_idx[pos] = tmp_cmap_atoms[6*pos + 5] - 1;
  }

  // Assign terms to each atom
  mvt.makeAtomAssignments();

  return mvt;
}

//-------------------------------------------------------------------------------------------------
AttenuationParameterSet condenseScreeningFactors(const BasicValenceTable &bvt,
                                                 const std::vector<double> &dihe_elec_screenings,
                                                 const std::vector<double> &dihe_vdw_screenings,
                                                 const double default_elec14_screening,
                                                 const double default_vdw14_screening) {
  AttenuationParameterSet attn_parm;

  // Sanity checks: did the topology file provide equal numbers of electrostatic and van-der Waals
  // screening factors?  If the topology provided screening factors, are there enough of them to
  // cover the dihedrals?
  const int nscreen_param = dihe_elec_screenings.size();
  if (nscreen_param != dihe_vdw_screenings.size()) {
    rtErr("The number of electrostatic and van-der Waals 1:4 screening factors are not equal (" +
          std::to_string(dihe_elec_screenings.size()) + " vs " +
          std::to_string(dihe_vdw_screenings.size()) + ", respectively).  No remedy is available.",
          "condenseScreeningFactors");
  }
  const int ndihe_param = (bvt.dihe_param_idx.size() > 0) ? maxValue(bvt.dihe_param_idx) + 1 : 0;
  if (nscreen_param > 0 && nscreen_param != ndihe_param) {
    rtErr("A total of " + std::to_string(nscreen_param) + " screening factors for 1:4 non-bonded "
          "interactions are provided, but do not appear to cover the range of dihedral parameter "
          "sets (" + std::to_string(ndihe_param) + ").  No remedy is available.",
          "condenseScreeningFactors");
  }

  // Count the number of unique screening factor pairs.  Make the (0.0, 0.0) pair an integral,
  // obligatory part of the set.  If there were no screening factors provided in the topology,
  // add the defaults to the set of parameter pairs and have all parameters point there.
  std::vector<int> set_map;
  attn_parm.elec_screening_factors.push_back(0.0);
  attn_parm.vdw_screening_factors.push_back(0.0);
  attn_parm.total_14_sets = 1;
  if (nscreen_param > 0) {
    set_map.resize(nscreen_param, -1);
    const Approx dzero(0.0);
    for (int i = 0; i < nscreen_param; i++) {
      if (dihe_elec_screenings[i] == dzero && dihe_vdw_screenings[i] == dzero) {
        set_map[i] = 0;
      }
    }
    int nsets = attn_parm.total_14_sets;
    for (int i = 0; i < nscreen_param; i++) {
      if (set_map[i] >= 0) {
        continue;
      }
      attn_parm.elec_screening_factors.push_back(dihe_elec_screenings[i]);
      attn_parm.vdw_screening_factors.push_back(dihe_vdw_screenings[i]);
      const Approx elecfac(dihe_elec_screenings[i]);
      const Approx vdwfac(dihe_vdw_screenings[i]);
      for (int j = i; j < nscreen_param; j++) {
        if (dihe_elec_screenings[j] == elecfac && dihe_vdw_screenings[j] == vdwfac) {
          set_map[j] = nsets;
        }
      }
      nsets++;
    }
    attn_parm.total_14_sets = nsets;
  }
  else {
    attn_parm.elec_screening_factors.push_back(default_elec14_screening);
    attn_parm.vdw_screening_factors.push_back(default_vdw14_screening);
    set_map.resize(ndihe_param, 1);
    attn_parm.total_14_sets += 1;
  }

  // Map all dihedrals to the correct scaling factors, using both their modifiers (to detect
  // dihedrals that do not control a 1:4 interaction between their I and L atoms) and their
  // parameter indices (to reference against the map created above).
  attn_parm.dihe14_parameter_indices.resize(bvt.total_dihes);
  for (int i = 0; i < bvt.total_dihes; i++) {
    const TorsionKind trkind = static_cast<TorsionKind>(bvt.dihe_mods[i].w);
    switch (trkind) {
    case TorsionKind::PROPER:
      attn_parm.dihe14_parameter_indices[i] = set_map[bvt.dihe_param_idx[i]];
      break;
    case TorsionKind::PROPER_NO_14:
    case TorsionKind::IMPROPER:
    case TorsionKind::IMPROPER_NO_14:
      attn_parm.dihe14_parameter_indices[i] = 0;
      break;
    }
  }

  return attn_parm;
}

//-------------------------------------------------------------------------------------------------
VirtualSiteTable listVirtualSites(const int expected_vsite_count,
                                  const std::string &file_name,
                                  const std::vector<double> &tmp_masses,
                                  const BasicValenceTable &bvt,
                                  const std::vector<double> &tmp_bond_equilibria,
                                  const std::vector<double> &tmp_angl_equilibria,
                                  const std::vector<char4> &tmp_atom_types,
                                  const std::vector<char4> &tmp_atom_names,
                                  const std::vector<int> &vsite_custom_frames,
                                  const std::vector<double> &vsite_custom_details) {

  // Count the number of massless particles and allocate the table
  const int natom = tmp_masses.size();
  std::vector<bool> is_vs(natom, false);
  std::vector<bool> is_custom_vs(natom, false);
  std::vector<bool> is_implicit_vs(natom, false);
  int n_virtual_site = 0;
  for (int i = 0; i < natom; i++) {
    if (tmp_masses[i] < constants::tiny) {
      is_vs[i] = true;
      n_virtual_site++;
    }
  }
  if (n_virtual_site != expected_vsite_count) {

    // If identifying virtual sites as particles with no mass did not work, the parm94.dat file
    // lists virtual sites as having a mass of 3.0 Daltons.  Try identifying the virtual sites by
    // atom type, EP or LP.
    for (int i = 0; i < natom; i++) {
      if ((tmp_atom_types[i].x == 'E' || tmp_atom_types[i].x == 'L') &&
          tmp_atom_types[i].y == 'P') {
        is_vs[i] = true;
        n_virtual_site++;
      }
    }
  }
  if (n_virtual_site != expected_vsite_count) {

    // As a last resort, try identifying virtual sites by names beginning with EP or LP.
    for (int i = 0; i < natom; i++) {
      if ((tmp_atom_names[i].x == 'E' || tmp_atom_names[i].x == 'L') &&
          tmp_atom_names[i].y == 'P') {
        is_vs[i] = true;
        n_virtual_site++;
      }
    }
  }

  // If the correct number of virtual sites still have not been found, abort.
  if (n_virtual_site != expected_vsite_count) {
    rtErr("Error.  Topology " + file_name + " states " + std::to_string(expected_vsite_count) +
          " virtual sites but " + std::to_string(n_virtual_site) +
          " massless particles were found.", "listVirtualSites");
  }
  VirtualSiteTable vst(natom, n_virtual_site);

  // If custom virtual site arrays are provided, identify all of the custom sites from amongst the
  // virtual sites in general.  All remaining sites must be implicit, based on bonding patterns.
  const int ncustom = vsite_custom_frames.size() / 6;
  if (static_cast<size_t>(ncustom) != vsite_custom_details.size() / 3) {
    rtErr("Error.  Information for " + std::to_string(ncustom) + " sets of custom frame indices "
          "and " + std::to_string(vsite_custom_details.size() / 3) +
          " sets of custom frame dimensions was provided.", "listVirtualSites");
  }
  for (int i = 0; i < ncustom; i++) {

    // Subtract 1 from the vsite_custom_frames atom indices here and elsewhere--they arrive
    // in this function still in Fortran format from reading the Amber topology file.
    is_custom_vs[vsite_custom_frames[6 * i] - 1] = true;
  }
  for (int i = 0; i < natom; i++) {
    is_implicit_vs[i] = (is_vs[i] && is_custom_vs[i] == false);
  }

  // If there are more virtual sites than custom entries, prepare for a long fight, searching the
  // bonds list over and over.  This process is improved by the facts that an implicit virtual
  // site must be bonded to one and only one other atom (its parent) and that none of the implicit
  // virtual site types can exist with parent atoms bonded to more than two other atoms.  At
  // present, all bonds are assumed to be harmonic bonds and therefore described in the "basic"
  // valence table.  However, if other sorts of specialized interactions, such as a force field
  // family-specific Morse oscillator potential, could be used to bond two atoms together, their
  // indexing tables will need to be included in this calculation as well.
  std::vector<int> bonds_to_real_atoms(natom, 0);
  std::vector<int> bonds_to_impvs(natom, 0);
  std::vector<int> parent_atoms(natom, -1);
  int2 init_int2 = {-1, -1};
  std::vector<int2> real_atom_neighbors(natom, init_int2);
  std::vector<int2> impvs_neighbors(natom, init_int2);
  std::vector<double> impvs_bond_lengths(natom, 0.0);
  for (int i = 0; i < bvt.total_bonds; i++) {

    // Shorthand for the bond's atom indicies
    const int atom_i = bvt.bond_i_atoms[i];
    const int atom_j = bvt.bond_j_atoms[i];

    // Check if a bonded interaction indicates some virtual site's parent atom
    int identified_vs = -1;
    int putative_parent = -1;
    if (is_implicit_vs[atom_i]) {
      identified_vs = atom_i;
      putative_parent = atom_j;
    }
    else if (is_implicit_vs[atom_j]) {
      identified_vs = atom_j;
      putative_parent = atom_i;
    }
    if (identified_vs >= 0) {
      if (parent_atoms[identified_vs] >= 0) {
        rtErr("Virtual site " + std::to_string(identified_vs + 1) + " is bonded to atoms " +
              std::to_string(parent_atoms[identified_vs] + 1) + " and " +
              std::to_string(putative_parent + 1) + " in topology " + file_name,
              "listVirtualSites");
      }
      parent_atoms[identified_vs] = putative_parent;
      impvs_bond_lengths[identified_vs] = tmp_bond_equilibria[bvt.bond_param_idx[i]];
      if (bonds_to_impvs[putative_parent] == 0) {
        impvs_neighbors[putative_parent].x = identified_vs;
      }
      else if (bonds_to_impvs[putative_parent] == 1) {
        impvs_neighbors[putative_parent].y = identified_vs;
      }
      bonds_to_impvs[putative_parent] += 1;
    }

    // Mark each atom as having additional connections to other atoms with mass, or to
    // implicit virtual sites.
    if (is_vs[atom_i] && is_vs[atom_j]) {
      rtErr("Virtual sites " + std::to_string(atom_i + 1) + " and " + std::to_string(atom_j + 1) +
            " are bonded to one another.  This is forbidden.", "listVirtualSites");
    }
    if (is_vs[atom_i] == false && is_vs[atom_j] == false) {

      // Both atoms of this connection have mass.  Mark them and record the actual identities
      // of up to two connecting atoms with mass.
      if (bonds_to_real_atoms[atom_i] == 0) {
        real_atom_neighbors[atom_i].x = atom_j;
      }
      else if (bonds_to_real_atoms[atom_i] == 1) {
        real_atom_neighbors[atom_i].y = atom_j;
      }
      if (bonds_to_real_atoms[atom_j] == 0) {
        real_atom_neighbors[atom_j].x = atom_i;
      }
      else if (bonds_to_real_atoms[atom_j] == 1) {
        real_atom_neighbors[atom_j].y = atom_i;
      }
      bonds_to_real_atoms[atom_i] += 1;
      bonds_to_real_atoms[atom_j] += 1;
    }
  }

  // Check that all implicit virtual sites have parent atoms, and that those parent atoms have
  // appropriate numbers of connections to other atoms with mass.  Also check that no atom with
  // mass has more than two implicit virtual sites, and that all atoms marked parent atoms do
  // in fact have implicit virtual sites attached to them.
  for (int i = 0; i < natom; i++) {
    if (is_implicit_vs[i] && parent_atoms[i] < 0) {
      rtErr("Implicit virtual site " + std::to_string(i + 1) + " in topology " + file_name +
            " has no identifiable parent atom.", "listVirtualSites");
    }
    if (is_implicit_vs[i] && bonds_to_real_atoms[parent_atoms[i]] > 2) {
      rtErr("Implicit virtual site " + std::to_string(i + 1) + " in topology " + file_name +
            " connects to a parent atom with " +
            std::to_string(bonds_to_real_atoms[parent_atoms[i]]) + " atoms.", "listVirtualSites");
    }
    if (is_vs[i] && bonds_to_impvs[i] > 0) {
      rtErr("Atom " + std::to_string(i + 1) + " in topology " + file_name + " is a virtual site, "
            "but has connections to " + std::to_string(bonds_to_impvs[i]) + " implicit virtual "
            "sites of its own.", "listVirtualSites");
    }
    if (bonds_to_impvs[i] > 2) {
      rtErr("Atom " + std::to_string(i + 1) + " is connected to " +
            std::to_string(bonds_to_impvs[i]) + " implicit virtual sites (the maximum is 2).",
            "listVirtualSites");
    }
  }

  // Catalog every virtual site.  The custom sites are already laid out.  Implicit sites can be
  // classified according to the bonded connections already detected.
  int pos = 0;
  int custom_vs_search_start = 0;
  for (int i = 0; i < natom; i++) {
    if (is_vs[i] == false) {
      vst.vs_numbers[i] = -1;
      continue;
    }
    vst.vs_atoms[pos] = i;
    vst.vs_numbers[i] = pos;
    if (is_custom_vs[i]) {

      // Ideally, the list of custom virtual sites will proceed such that the topological atom
      // indices of the virtual sites occur in increasing order.  If that is not the case, the
      // rest of the list must be searched until the right custom virtual site is found.
      int cidx = -1;
      if (vsite_custom_frames[6 * custom_vs_search_start] - 1 == i) {
        cidx = custom_vs_search_start;
        custom_vs_search_start++;
      }
      else {
        for (int j = custom_vs_search_start; j < ncustom; j++) {
          if (vsite_custom_frames[6 * j] - 1 == i) {
            cidx = j;
            break;
          }
        }
      }
      if (cidx < 0) {
        rtErr("Atom " + std::to_string(i) + " is a custom virtual site but no frame indexing "
              "could be identified for it.", "listVirtualSites");
      }
      vst.frame1_atoms[pos] = vsite_custom_frames[6*cidx + 1] - 1;
      vst.frame2_atoms[pos] = vsite_custom_frames[6*cidx + 2] - 1;
      vst.frame3_atoms[pos] = vsite_custom_frames[6*cidx + 3] - 1;
      vst.frame4_atoms[pos] = vsite_custom_frames[6*cidx + 4] - 1;
      vst.frame_types[pos]  = vsite_custom_frames[6*cidx + 5];
      vst.frame_dim1[pos]   = vsite_custom_details[3*cidx    ];
      vst.frame_dim2[pos]   = vsite_custom_details[3*cidx + 1];
      vst.frame_dim3[pos]   = vsite_custom_details[3*cidx + 2];
    }
    else {

      // Implicit virtual sites' parent atoms have been vetted and their parent atoms identified.
      // Determine the frame type, then record the frame atoms, and finally get frame dimensions
      // either from bond lengths or by convention.
      vst.frame1_atoms[pos] = parent_atoms[i];
      if (bonds_to_real_atoms[parent_atoms[i]] == 1) {

        // In the Amber code, this becomes a "type II" (frame index 0 or 2) virtual site,
        // applicable to carbonyl oxygen lone pairs or sigma holes on halo-organic compounds,
        // respectively.  STORMM will approximate those frames with the appropriate GROMACS frame
        // types: fixed-distance, two-atom frames (halogen sigma hole case) and fixed angle and
        // distance, three-atom frames (carbonyl lone pair case).  The GROMACS types are better
        // suited for most situations and more easily understood.
        vst.frame2_atoms[pos] = real_atom_neighbors[parent_atoms[i]].x;
        if (bonds_to_real_atoms[vst.frame2_atoms[pos]] != 2) {
          rtErr("Implicit virtual site " + std::to_string(i + 1) + " bonds to atom " +
                std::to_string(parent_atoms[i] + 1) + ", appearing as a type II frame, but the "
                "parent bonds to " + std::to_string(vst.frame2_atoms[pos] + 1) + ", which in turn "
                "has not 2 but " + std::to_string(bonds_to_real_atoms[vst.frame2_atoms[pos]]) +
                " real atom neighbors.  This is cannot be a type II frame.", "listVirtualSites");
        }
        vst.frame4_atoms[pos] = -1;
        if (bonds_to_impvs[parent_atoms[i]] == 2) {

          // Place the two sites at fixed angles of +/- 120 degrees from the parent-frame atom
          // line, at a fixed distance defined by the bond length between the virtual site and its
          // parent atom.  Take the first atom as the reference point for defining the sign of the
          // angle (the choice is not really consequential)
          vst.frame_types[pos] = static_cast<int>(VirtualSiteKind::FAD_3);
          vst.frame_dim1[pos]  = impvs_bond_lengths[i];
          vst.frame_dim2[pos]  = 2.0 * symbols::pi / 3.0;
          if (i == impvs_neighbors[parent_atoms[i]].y) {
            vst.frame_dim2[pos] *= -1.0;
          }
          else if (i != impvs_neighbors[parent_atoms[i]].x) {
            rtErr("Implicit virtual site " + std::to_string(i + 1) + " was not listed among the "
                  "substituents of atom " + std::to_string(parent_atoms[i] + 1) + ".",
                  "listVirtualSites");
          }
          vst.frame3_atoms[pos] = real_atom_neighbors[vst.frame2_atoms[pos]].x;
        }
        else {

          // Place the one site at a fixed distance from the parent atom, along the line with
          // the other frame atom.
          vst.frame_types[pos] = static_cast<int>(VirtualSiteKind::FIXED_2);
          vst.frame_dim1[pos]  = impvs_bond_lengths[i];
          vst.frame3_atoms[pos] = -1;
        }
      }
      else if (bonds_to_real_atoms[parent_atoms[i]] == 2) {

        // In the Amber code, this becomes a "type I" virtual site, type 0 for TIP5P and type 2
        // for TIP4P.  GROMACS frames offer exact substitutes for both of these important cases.
        vst.frame2_atoms[pos] = real_atom_neighbors[parent_atoms[i]].x;
        vst.frame3_atoms[pos] = real_atom_neighbors[parent_atoms[i]].y;
        vst.frame4_atoms[pos] = -1;

        // The TIP5P arrangement gets type 1 frame designation (TIP4P gets 3... pmemd-inspired
        // convention)
        if (bonds_to_impvs[parent_atoms[i]] == 2) {
          vst.frame_types[pos] = static_cast<int>(VirtualSiteKind::OUT_3);

          // Find the bond lengths between the parent atom and its two substituents.  These
          // will determine the way that the frame is laid out.
          double pf2, pf3, f2f3, f2_p_f3;
          bool pf2_found = false;
          bool pf3_found = false;
          bool f2f3_found = false;
          bool f2_p_f3_found = false;
          const int p_atom = parent_atoms[i];
          const int f2_atom = real_atom_neighbors[p_atom].x;
          const int f3_atom = real_atom_neighbors[p_atom].y;
          for (int j = bvt.bond_assigned_bounds[p_atom];
               j < bvt.bond_assigned_bounds[p_atom + 1]; j++) {
            if (bvt.bond_assigned_atoms[j] == f2_atom) {
              pf2 = tmp_bond_equilibria[bvt.bond_param_idx[bvt.bond_assigned_terms[j]]];
              pf2_found = true;
            }
            if (bvt.bond_assigned_atoms[j] == f3_atom) {
              pf3 = tmp_bond_equilibria[bvt.bond_param_idx[bvt.bond_assigned_terms[j]]];
              pf3_found = true;
            }
          }
          for (int j = bvt.bond_assigned_bounds[f2_atom];
               j < bvt.bond_assigned_bounds[f2_atom + 1]; j++) {
            if (bvt.bond_assigned_atoms[j] == p_atom) {
              pf2 = tmp_bond_equilibria[bvt.bond_param_idx[bvt.bond_assigned_terms[j]]];
              pf2_found = true;
            }
            if (bvt.bond_assigned_atoms[j] == f3_atom) {
              f2f3 = tmp_bond_equilibria[bvt.bond_param_idx[bvt.bond_assigned_terms[j]]];
              f2f3_found = true;
            }
          }
          for (int j = bvt.bond_assigned_bounds[f3_atom];
               j < bvt.bond_assigned_bounds[f3_atom + 1]; j++) {
            if (bvt.bond_assigned_atoms[j] == p_atom) {
              pf3 = tmp_bond_equilibria[bvt.bond_param_idx[bvt.bond_assigned_terms[j]]];
              pf3_found = true;
            }
            if (bvt.bond_assigned_atoms[j] == f2_atom) {
              f2f3 = tmp_bond_equilibria[bvt.bond_param_idx[bvt.bond_assigned_terms[j]]];
              f2f3_found = true;
            }
          }

          // Check that the trivial bond lengths were found
          if (pf2_found == false || pf3_found == false) {
            rtErr("Bond lengths could not be located between parent atom " +
                  char4ToString(tmp_atom_names[p_atom]) + " " + std::to_string(p_atom + 1) +
                  " and frame atoms " + char4ToString(tmp_atom_names[f2_atom]) + " " +
                  std::to_string(f2_atom + 1) + " / " + char4ToString(tmp_atom_names[f3_atom]) +
                  " " + std::to_string(f3_atom + 1) + ".", "listVirtualSites");
          }

          // If the equilibrium distance between frame atoms cannot be found, seek out the angle
          if (f2f3_found == false) {
            for (int j = bvt.angl_assigned_bounds[p_atom];
               j < bvt.angl_assigned_bounds[p_atom + 1]; j++) {
              if ((bvt.angl_assigned_atoms[ 2 * j     ] == f2_atom &&
                   bvt.angl_assigned_atoms[(2 * j) + 1] == f3_atom) ||
                  (bvt.angl_assigned_atoms[ 2 * j     ] == f3_atom &&
                   bvt.angl_assigned_atoms[(2 * j) + 1] == f2_atom)) {
                f2_p_f3 = tmp_angl_equilibria[bvt.angl_param_idx[bvt.angl_assigned_terms[j]]];
                f2_p_f3_found = true;
              }
            }
            if (f2_p_f3_found == false) {
              rtErr("Neither a bond between frame atoms " +
                    char4ToString(tmp_atom_names[f2_atom]) + " " + std::to_string(f2_atom + 1) +
                    " and " + char4ToString(tmp_atom_names[f3_atom]) + " " +
                    std::to_string(f3_atom + 1) + " nor an angle between these atoms and parent "
                    "atom " + char4ToString(tmp_atom_names[p_atom]) + " " +
                    std::to_string(p_atom + 1) + " could be found.", "listVirtualSites");
            }
            const double dx = pf2 - (pf3 * cos(f2_p_f3));
            const double dy = pf3 * sin(f2_p_f3);
            f2f3 = sqrt((dx * dx) + (dy * dy));
          }

          // A balanced ratio of parent->frame2 and parent->frame3 vectors will trace a line down
          // the bisector of angle frame2 - parent - frame3.  The magnitude of the displacement
          // from the parent atom will be determined by a sacling factor and the distance between
          // the frame2 and frame3 atoms at equilibrium.  The total displacement needed is given
          // by the cosine of half the tetrahedral angle times the bond length from the parent
          // atom to each virtual site.
          const double disp_needed = cos(0.5 * tetrahedral_angle) * impvs_bond_lengths[i];
          const double scale_pf2f3 = disp_needed / f2f3;
          vst.frame_dim1[pos] = -scale_pf2f3;
          vst.frame_dim2[pos] = -scale_pf2f3;

          // The cross product of parent->frame2 and parent->frame3 vectors must be scaled by a
          // factor so as to move the virtual site out of the plane by the sine of half the
          // tetrahedral angle times the bond length from the parent atom to each virtual site.
          double pf2v[3], pf3v[3], cr_pf2f3[3];
          pf2v[0] = pf2;
          pf2v[1] = 0.0;
          pf2v[2] = 0.0;
          if (f2_p_f3_found == false) {

            // Law of cosines: a^2 = b^2 + c^2 - 2 * b * c * cos(A).  Here, a^2 is (f2f3^2) and
            // b and c are (pf2) and (pf3), respectively.
            f2_p_f3 = acos(((pf2 * pf2) + (pf3 * pf3) - (f2f3 * f2f3)) / (2.0 * pf2 * pf3));
          }
          pf3v[0] = pf3 * cos(f2_p_f3);
          pf3v[1] = pf3 * sin(f2_p_f3);
          pf3v[2] = 0.0;
          crossProduct(pf2v, pf3v, cr_pf2f3);
          const double mcr = magnitude(cr_pf2f3, 3);
          vst.frame_dim3[pos] = sin(0.5 * tetrahedral_angle) * impvs_bond_lengths[i] / mcr;

          // Flip the second of the two virtual sites in how it is displaced from the plane.
          // Check that both virtual sites are substituents of the same parent atom.
          if (i == impvs_neighbors[parent_atoms[i]].y) {
            vst.frame_dim3[pos] *= -1.0;
          }
          else if (i != impvs_neighbors[parent_atoms[i]].x) {
            rtErr("Implicit virtual site " + std::to_string(i + 1) + " was not listed among the "
                  "substituents of atom " + std::to_string(parent_atoms[i] + 1) + ".",
                  "listVirtualSites");
          }
        }
        else {
          vst.frame_types[pos] = static_cast<int>(VirtualSiteKind::FIXED_3);
          vst.frame_dim1[pos]  = impvs_bond_lengths[i];
          vst.frame_dim2[pos]  = 0.5;
        }
      }
      else {
        rtErr("Atom " + std::to_string(parent_atoms[i] + 1) + ", parent atom of virtual site " +
              std::to_string(i + 1) + ", bonds to an invalid number of other atoms (" +
              std::to_string(bonds_to_real_atoms[parent_atoms[i]]) + ").", "listVirtualSites");
      }
    }
    pos++;
  }

  // Condense the tables of frame type and dimensions.  Assign each virtual site a parameter
  // index.
  int n_unique_frames = 0;
  std::vector<bool> coverage(vst.vs_count, false);
  for (int i = 0; i < vst.vs_count; i++) {
    if (coverage[i]) {
      continue;
    }
    const int frtype = vst.frame_types[i];
    const Approx dim1a(vst.frame_dim1[i], constants::tiny);
    const Approx dim2a(vst.frame_dim2[i], constants::tiny);
    const Approx dim3a(vst.frame_dim3[i], constants::tiny);
    for (int j = i; j < vst.vs_count; j++) {
      if (vst.frame_types[j] == frtype && dim1a.test(vst.frame_dim1[j]) &&
          dim2a.test(vst.frame_dim2[j]) && dim3a.test(vst.frame_dim3[j])) {
        vst.param_idx[j] = n_unique_frames;
        coverage[j] = true;
      }
    }
    vst.frame_types[n_unique_frames] = frtype;
    vst.frame_dim1[n_unique_frames]  = vst.frame_dim1[i];
    vst.frame_dim2[n_unique_frames]  = vst.frame_dim2[i];
    vst.frame_dim3[n_unique_frames]  = vst.frame_dim3[i];
    n_unique_frames++;
  }
  vst.frame_types.resize(n_unique_frames);
  vst.frame_dim1.resize(n_unique_frames);
  vst.frame_dim2.resize(n_unique_frames);
  vst.frame_dim3.resize(n_unique_frames);
  
  return vst;
}

//-------------------------------------------------------------------------------------------------
void accumulateExclusionBounds(std::vector<int> *current_bounds, const int atom_a,
                               const int atom_b, const std::vector<int> &vsite_child_bounds,
                               const std::vector<int> &vsite_child_list) {
  
  // Extend each real atom's bounds to include exclusions with the other's virtual site children
  const int atom_a_cluster = vsite_child_bounds[atom_a + 1] - vsite_child_bounds[atom_a] + 1;
  const int atom_b_cluster = vsite_child_bounds[atom_b + 1] - vsite_child_bounds[atom_b] + 1;
  current_bounds->data()[atom_a] += atom_b_cluster;
  current_bounds->data()[atom_b] += atom_a_cluster;

  // Extend the bounds for each virtual site child of atom_i to include atom_j plus its children
  for (int i = vsite_child_bounds[atom_a]; i < vsite_child_bounds[atom_a + 1]; i++) {
    current_bounds->data()[vsite_child_list[i]] += atom_b_cluster;
  }

  // Extend the bounds for each virtual site child of atom_j to include atom_i plus its children
  for (int i = vsite_child_bounds[atom_b]; i < vsite_child_bounds[atom_b + 1]; i++) {
    current_bounds->data()[vsite_child_list[i]] += atom_a_cluster;
  }
}

//-------------------------------------------------------------------------------------------------
void markPairExclusion(std::vector<int> *excl_list, std::vector<int> *counters,
                       const std::vector<int> &bounds, const int atom_i, const int atom_j) {

  excl_list->data()[bounds[atom_i] + counters->data()[atom_i]] = atom_j;
  excl_list->data()[bounds[atom_j] + counters->data()[atom_j]] = atom_i;
  counters->data()[atom_i] += 1;
  counters->data()[atom_j] += 1;
}

//-------------------------------------------------------------------------------------------------
void markExclusions(std::vector<int> *excl_list, std::vector<int> *excl_counters,
                    const std::vector<int> &excl_bounds, const int atom_a, const int atom_b,
                    const std::vector<int> &vsite_child_bounds,
                    const std::vector<int> &vsite_child_list) {

  // Mark the exclusions between the two real atoms
  markPairExclusion(excl_list, excl_counters, excl_bounds, atom_a, atom_b);

  // Mark exclusions between the first atom and virtual site children of the second
  for (int i = vsite_child_bounds[atom_b]; i < vsite_child_bounds[atom_b + 1]; i++) {
    markPairExclusion(excl_list, excl_counters, excl_bounds, atom_a, vsite_child_list[i]);
  }

  // Loop over the first atom's virtual site children
  for (int i = vsite_child_bounds[atom_a]; i < vsite_child_bounds[atom_a + 1]; i++) {

    // Mark exclusions with the second real atom
    markPairExclusion(excl_list, excl_counters, excl_bounds, vsite_child_list[i], atom_b);

    // Mark exclusions with virtual site children of the second real atom
    for (int j = vsite_child_bounds[atom_b]; j < vsite_child_bounds[atom_b + 1]; j++) {
      markPairExclusion(excl_list, excl_counters, excl_bounds, vsite_child_list[i],
                        vsite_child_list[j]);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void cullDuplicateExclusions(std::vector<int> *excl_list, std::vector<int> *excl_bounds) {

  // Step through the atoms one by one identifying the unique pair partners, then compact
  // the lists within the original list and original bounds.  Each atom will then have a
  // list of unique pair partners followed by garbage in indices of excl_list that are no
  // longer needed.
  const int natom = excl_bounds->size() - 1;
  int *list_data = excl_list->data();
  int *bounds_data = excl_bounds->data();
  std::vector<int> new_counts(natom + 1, 0);
  int largest_group = 0;
  for (int i = 0; i < natom; i++) {
    largest_group = std::max(largest_group, bounds_data[i + 1] - bounds_data[i]);
  }
  std::vector<int> excl_tmp_list(largest_group, 0);
  for (int i = 0; i < natom; i++) {
    int unique_count = 0;
    for (int j = bounds_data[i]; j < bounds_data[i + 1]; j++) {

      // Cull -1 entries: these are blanks where extra space was allocated but not filled.  Also
      // Cull self entries: non-bonded 1:2 must cull 1:1 prior exclusions, 1:3 exclusions must
      // cull 1:2 and 1:1, but this is a 1:0 culling.  Nothing excludes non-bonded interactions
      // with itself explicitly because it's always done implicitly in the non-bonded tiles.
      if (list_data[j] == -1 || list_data[j] == i) {
        continue;
      }
      bool found = false;
      for (int k = 0; k < unique_count; k++) {
        found = (found || list_data[j] == excl_tmp_list[k]);
      }
      if (found) {
        continue;
      }
      excl_tmp_list[unique_count] = list_data[j];
      unique_count++;
    }
    new_counts[i] = unique_count;
    for (int j = 0; j < unique_count; j++) {
      list_data[bounds_data[i] + j] = excl_tmp_list[j];
    }
  }

  // Compact the entries of excl_list based on the new limits
  int total_excl_count = 0;
  for (int i = 0; i < natom; i++) {
    const int i_offset = bounds_data[i];
    for (int j = 0; j < new_counts[i]; j++) {
      list_data[total_excl_count] = list_data[i_offset + j];
      total_excl_count++;
    }
  }
  prefixSumInPlace(&new_counts, PrefixSumType::EXCLUSIVE, "cullDuplicateExclusions");
  for (int i = 0; i < natom + 1; i++) {
    bounds_data[i] = new_counts[i];
  }
}

//-------------------------------------------------------------------------------------------------
void cullPriorExclusions(std::vector<int> *excl_list, std::vector<int> *excl_bounds,
                         const std::vector<int> &prior_list,
                         const std::vector<int> &prior_bounds) {

  // Step through the atoms, again identifying unique pair partners but this time by comparing
  // against some other list (whose elements are assumed to be unique themselves).  Condense and
  // compact the list as was done in cullDuplicateExclusions.
  const int natom = excl_bounds->size() - 1;
  int *list_data = excl_list->data();
  int *bounds_data = excl_bounds->data();
  std::vector<int> new_counts(natom + 1, 0);
  int largest_group = 0;
  for (int i = 0; i < natom; i++) {
    largest_group = std::max(largest_group, bounds_data[i + 1] - bounds_data[i]);    
  }
  std::vector<int> excl_tmp_list(largest_group, 0);
  for (int i = 0; i < natom; i++) {
    int unique_count = 0;
    for (int j = bounds_data[i]; j < bounds_data[i + 1]; j++) {

      // There should be no -1's in the list by this point, but cull them anyway and warn
      if (list_data[j] == -1 || list_data[j] == i) {
        std::string other_atom = (list_data[j] == -1) ? "a blank atom" : "itself";
        rtWarn("Atom " + std::to_string(i + 1) + " excludes " + other_atom + ".",
               "cullPriorExclusions");
        continue;
      }
      if (list_data[j] == i) {
        continue;
      }
      const int excluded_atom = list_data[j];
      bool found = false;
      for (int k = prior_bounds[i]; k < prior_bounds[i + 1]; k++) {
        found = (found || (prior_list[k] == excluded_atom));
      }
      if (found) {
        continue;
      }
      excl_tmp_list[unique_count] = list_data[j];
      unique_count++;
    }
    new_counts[i] = unique_count;
    for (int j = 0; j < unique_count; j++) {
      list_data[bounds_data[i] + j] = excl_tmp_list[j];
    }
  }

  // Compact the entries of excl_list based on the new limits
  int total_excl_count = 0;
  for (int i = 0; i < natom; i++) {
    const int i_offset = bounds_data[i];
    for (int j = 0; j < new_counts[i]; j++) {
      list_data[total_excl_count] = list_data[i_offset + j];
      total_excl_count++;
    }
  }
  prefixSumInPlace(&new_counts, PrefixSumType::EXCLUSIVE, "cullPriorExclusions");
  for (int i = 0; i < natom + 1; i++) {
    bounds_data[i] = new_counts[i];
  }
}

//-------------------------------------------------------------------------------------------------
Map1234 mapExclusions(const int atom_count, const std::vector<int> &vs_particles,
                      const std::vector<int> &vs_parent_atoms,
                      const std::vector<int> &bond_i_atoms, const std::vector<int> &bond_j_atoms) {
  Map1234 result;

  // Allocate arrays for the numbers of exclusions
  result.nb11_excl_bounds.resize(atom_count + 1);
  result.nb12_excl_bounds.resize(atom_count + 1);
  result.nb13_excl_bounds.resize(atom_count + 1);
  result.nb14_excl_bounds.resize(atom_count + 1);
  
  // Allocate an extra array of counters that will serve for each atom in various situations
  std::vector<int> atom_counters(atom_count, 0);
  if (vs_particles.size() != vs_parent_atoms.size()) {
    rtErr("The number of virtual site particle indices and parent atoms (" +
          std::to_string(vs_particles.size()) + ", " + std::to_string(vs_parent_atoms.size()) +
          ") must be the same.", "mapExclusions");
  }
  const int virtual_site_count = vs_particles.size();

  // Make a list of all virtual site children for every parent atom.  Because every virtual site
  // will have one and only one parent atom, this list will tell us how many 1:1 interactions to
  // expect, in addition to 1:1 interactions (exclusions) between virtual sites and their parents.
  std::vector<int> vsite_child_bounds(atom_count + 1, 0);
  for (int i = 0; i < virtual_site_count; i++) {
    vsite_child_bounds[vs_parent_atoms[i]] += 1;
  }
  prefixSumInPlace(&vsite_child_bounds, PrefixSumType::EXCLUSIVE, "mapExclusions");
  std::vector<int> vsite_child_list(vsite_child_bounds[atom_count], -1);
  for (int i = 0; i < virtual_site_count; i++) {
    const int parent = vs_parent_atoms[i];
    vsite_child_list[vsite_child_bounds[parent] + atom_counters[parent]] = vs_particles[i];
    atom_counters[parent] += 1;
  }

  // List all 1:1 interactions as virtual sites interacting with their parents.  Virtual sites
  // that have the same parent atom are also 1:1 to one another.
  for (int i = 0; i < virtual_site_count; i++) {
    const int vsite  = vs_particles[i];
    const int parent = vs_parent_atoms[i];
    result.nb11_excl_bounds[vsite] += vsite_child_bounds[parent + 1] -
                                      vsite_child_bounds[parent];
    result.nb11_excl_bounds[parent] += 1;
  }

  // Compute the capped prefix sum over all atoms, then resize the 1:1 exclusion list as needed
  prefixSumInPlace(&result.nb11_excl_bounds, PrefixSumType::EXCLUSIVE, "mapExclusions");
  result.nb11_excl_list.resize(result.nb11_excl_bounds[atom_count], -1);

  // Loop over all atoms, composing the list of 1:1 exclusions
  for (int i = 0; i < atom_count; i++) {
    atom_counters[i] = 0;
  }
  for (int i = 0; i < atom_count; i++) {
    for (int j = vsite_child_bounds[i]; j < vsite_child_bounds[i + 1]; j++) {
      markPairExclusion(&result.nb11_excl_list, &atom_counters, result.nb11_excl_bounds, i,
                        vsite_child_list[j]);
      for (int k = vsite_child_bounds[i]; k < j; k++) {
        markPairExclusion(&result.nb11_excl_list, &atom_counters, result.nb11_excl_bounds,
                          vsite_child_list[j], vsite_child_list[k]);
      }
    }
  }
  
  // Estimate the 1:2 exclusion bounds, then populate them, then cull duplicates
  if (bond_i_atoms.size() != bond_j_atoms.size()) {
    rtErr("The number of bonded atoms (" + std::to_string(bond_i_atoms.size()) + ", " +
          std::to_string(bond_j_atoms.size()) + ") must be the same in both I and J arrays.",
          "mapExclusions");
  }
  const int bond_count = bond_i_atoms.size();
  for (int i = 0; i < atom_count + 1; i++) {
    result.nb12_excl_bounds[i] = 0;
  }
  for (int i = 0; i < bond_count; i++) {
    accumulateExclusionBounds(&result.nb12_excl_bounds, bond_i_atoms[i], bond_j_atoms[i],
                              vsite_child_bounds, vsite_child_list);
  }
  prefixSumInPlace(&result.nb12_excl_bounds, PrefixSumType::EXCLUSIVE, "mapExclusions");
  result.nb12_excl_list.resize(result.nb12_excl_bounds[atom_count], -1);
  for (int i = 0; i < atom_count; i++) {
    atom_counters[i] = 0;
  }
  for (int i = 0; i < bond_count; i++) {
    markExclusions(&result.nb12_excl_list, &atom_counters, result.nb12_excl_bounds,
                   bond_i_atoms[i], bond_j_atoms[i], vsite_child_bounds, vsite_child_list);
  }
  cullDuplicateExclusions(&result.nb12_excl_list, &result.nb12_excl_bounds);
  cullPriorExclusions(&result.nb12_excl_list, &result.nb12_excl_bounds, result.nb11_excl_list,
                      result.nb11_excl_bounds);
  result.nb12_excl_list.resize(result.nb12_excl_bounds[atom_count]);
  result.nb12_excl_list.shrink_to_fit();
  
  // Estimate the 1:3 exclusion bounds, then populate them, then cull duplicates
  for (int i = 0; i < atom_count + 1; i++) {
    result.nb13_excl_bounds[i] = 0;
  }
  for (int i = 0; i < atom_count; i++) {
    const int atom_i = i;
    for (int j = result.nb12_excl_bounds[atom_i]; j < result.nb12_excl_bounds[atom_i + 1]; j++) {
      const int atom_j = result.nb12_excl_list[j];
      for (int k = result.nb12_excl_bounds[atom_j]; k < result.nb12_excl_bounds[atom_j + 1]; k++) {
        const int atom_k = result.nb12_excl_list[k];
        if (atom_k != atom_i) {
          accumulateExclusionBounds(&result.nb13_excl_bounds, atom_i, atom_k,
                                    vsite_child_bounds, vsite_child_list);
        }
      }
    }
  }
  prefixSumInPlace(&result.nb13_excl_bounds, PrefixSumType::EXCLUSIVE, "mapExclusions");
  result.nb13_excl_list.resize(result.nb13_excl_bounds[atom_count], -1);
  for (int i = 0; i < atom_count; i++) {
    atom_counters[i] = 0;
  }
  for (int i = 0; i < atom_count; i++) {
    const int atom_i = i;
    for (int j = result.nb12_excl_bounds[atom_i]; j < result.nb12_excl_bounds[atom_i + 1]; j++) {
      const int atom_j = result.nb12_excl_list[j];
      for (int k = result.nb12_excl_bounds[atom_j]; k < result.nb12_excl_bounds[atom_j + 1]; k++) {
        const int atom_k = result.nb12_excl_list[k];
        if (atom_k != atom_i) {
          markExclusions(&result.nb13_excl_list, &atom_counters, result.nb13_excl_bounds, atom_i,
                         atom_k, vsite_child_bounds, vsite_child_list);
        }
      }
    }
  }
  cullDuplicateExclusions(&result.nb13_excl_list, &result.nb13_excl_bounds);
  cullPriorExclusions(&result.nb13_excl_list, &result.nb13_excl_bounds, result.nb11_excl_list,
                      result.nb11_excl_bounds);
  cullPriorExclusions(&result.nb13_excl_list, &result.nb13_excl_bounds, result.nb12_excl_list,
                      result.nb12_excl_bounds);
  result.nb13_excl_list.resize(result.nb13_excl_bounds[atom_count]);
  result.nb13_excl_list.shrink_to_fit();
  
  // Repeat the process for the 1:4 exclusions: estimate, populate, cull
  for (int i = 0; i < atom_count + 1; i++) {
    result.nb14_excl_bounds[i] = 0;
  }
  for (int i = 0; i < atom_count; i++) {
    const int atom_i = i;
    for (int j = result.nb12_excl_bounds[atom_i]; j < result.nb12_excl_bounds[atom_i + 1]; j++) {
      const int atom_j = result.nb12_excl_list[j];
      for (int k = result.nb12_excl_bounds[atom_j]; k < result.nb12_excl_bounds[atom_j + 1]; k++) {
        const int atom_k = result.nb12_excl_list[k];
        if (atom_k != atom_i) {
          for (int m = result.nb12_excl_bounds[atom_k]; m < result.nb12_excl_bounds[atom_k + 1];
               m++) {
            const int atom_m = result.nb12_excl_list[m];
            if (atom_m != atom_i && atom_m != atom_j) {
              accumulateExclusionBounds(&result.nb14_excl_bounds, atom_i, atom_m,
                                        vsite_child_bounds, vsite_child_list);
            }
          }
        }
      }
    }
  }
  prefixSumInPlace(&result.nb14_excl_bounds, PrefixSumType::EXCLUSIVE, "mapExclusions");
  result.nb14_excl_list.resize(result.nb14_excl_bounds[atom_count], -1);
  for (int i = 0; i < atom_count; i++) {
    atom_counters[i] = 0;
  }
  for (int i = 0; i < atom_count; i++) {
    const int atom_i = i;
    for (int j = result.nb12_excl_bounds[atom_i]; j < result.nb12_excl_bounds[atom_i + 1]; j++) {
      const int atom_j = result.nb12_excl_list[j];
      for (int k = result.nb12_excl_bounds[atom_j]; k < result.nb12_excl_bounds[atom_j + 1]; k++) {
        const int atom_k = result.nb12_excl_list[k];
        if (atom_k != atom_i) {
          for (int m = result.nb12_excl_bounds[atom_k]; m < result.nb12_excl_bounds[atom_k + 1];
               m++) {
            const int atom_m = result.nb12_excl_list[m];
            if (atom_m != atom_i && atom_m != atom_j) {
              markExclusions(&result.nb14_excl_list, &atom_counters, result.nb14_excl_bounds,
                             atom_i, atom_m, vsite_child_bounds, vsite_child_list);
            }
          }
        }
      }
    }
  }
  cullDuplicateExclusions(&result.nb14_excl_list, &result.nb14_excl_bounds);
  cullPriorExclusions(&result.nb14_excl_list, &result.nb14_excl_bounds, result.nb11_excl_list,
                      result.nb11_excl_bounds);
  cullPriorExclusions(&result.nb14_excl_list, &result.nb14_excl_bounds, result.nb12_excl_list,
                      result.nb12_excl_bounds);
  cullPriorExclusions(&result.nb14_excl_list, &result.nb14_excl_bounds, result.nb13_excl_list,
                      result.nb13_excl_bounds);
  result.nb14_excl_list.resize(result.nb14_excl_bounds[atom_count]);
  result.nb14_excl_list.shrink_to_fit();

  return result;
}

//-------------------------------------------------------------------------------------------------
Map1234 mapExclusions(const int atom_count, const BasicValenceTable &bvt,
                      const VirtualSiteTable &vst) {
  return mapExclusions(atom_count, vst.vs_atoms, vst.frame1_atoms, bvt.bond_i_atoms,
                       bvt.bond_j_atoms);
}

//-------------------------------------------------------------------------------------------------
void exclusionSearchLoop(std::vector<int> &coverage, const std::vector<int> &bounds,
                         const std::vector<int> &excl_list, const CondensedExclusions &ce) {

  // Deduce the number of atoms in the system
  const int natom = ce.atom_excl_bounds.size() - 1;

  // Loop over each exclusion for every atom, incrementing the coverage array
  for (int i = 0; i < natom; i++) {
    for (int j = bounds[i]; j < bounds[i + 1]; j++) {
      const int excluded_atom = excl_list[j];
      bool found = false;
      for (int k = ce.atom_excl_bounds[i]; k < ce.atom_excl_bounds[i + 1]; k++) {
        if (ce.atom_excl_list[k] == excluded_atom) {
          coverage[k] += 1;
          found = true;
          break;
        }
      }
      for (int k = ce.atom_excl_bounds[excluded_atom];
           k < ce.atom_excl_bounds[excluded_atom + 1]; k++) {
        if (ce.atom_excl_list[k] == i) {
          coverage[k] += 1;
          found = true;
          break;
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void checkExclusions(const CondensedExclusions &ce, const Map1234 &all_nb_excl,
                     const std::string &file_name) {

  // Make masks to detect whether an interaction has been accounted for
  std::vector<int> coverage(ce.atom_excl_list.size(), 0);

  // Add coverage for 1:1 interactions
  exclusionSearchLoop(coverage, all_nb_excl.nb11_excl_bounds, all_nb_excl.nb11_excl_list, ce);

  // Add coverage for 1:2 interactions
  exclusionSearchLoop(coverage, all_nb_excl.nb12_excl_bounds, all_nb_excl.nb12_excl_list, ce);

  // Add coverage for 1:3 interactions
  exclusionSearchLoop(coverage, all_nb_excl.nb13_excl_bounds, all_nb_excl.nb13_excl_list, ce);

  // Add coverage for 1:4 interactions
  exclusionSearchLoop(coverage, all_nb_excl.nb14_excl_bounds, all_nb_excl.nb14_excl_list, ce);

  // Check every atom's exclusions to verify that each is covered exactly twice
  const int natom = ce.atom_excl_bounds.size() - 1;
  for (int i = 0; i < natom; i++) {
    for (int j = ce.atom_excl_bounds[i]; j < ce.atom_excl_bounds[i + 1]; j++) {
      if (coverage[j] != 2) {
        rtWarn("Exclusion " + std::to_string(i + 1) + " -> " +
               std::to_string(ce.atom_excl_list[j] + 1) + " was found in topology file " +
               file_name + ", but present " + std::to_string(coverage[j]) +
               " times in the derived exclusion lists.", "checkExclusions");
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<int> atomicNumbersFromMasses(const std::vector<double> &masses,
                                         const std::vector<char4> &atom_names,
                                         const std::vector<int> &nb12_list,
                                         const std::vector<int> &nb12_bounds,
                                         const std::string &source,
                                         const ExceptionResponse policy) {
  const int natom = masses.size();
  if (static_cast<int>(nb12_bounds.size() - 1) != natom) {
    rtErr("The atom masses array and bounds for 1:2 exclusions indicate different numbers of "
          "atoms (" + std::to_string(natom) + " vs. " + std::to_string(nb12_bounds.size() - 1) +
          ").", "atomicNumbersFromMasses");
  }

  // Attempt to identify atom numbers right now, starting with masses and also picking out atoms
  // with extra point (EP) or lone pair (LP) names
  std::vector<int> z_numbers = massToZNumber(masses);
  int n_named_vsite = 0;
  for (int i = 0; i < natom; i++) {
    if (z_numbers[i] == 0) {
      continue;
    }
    if (((atom_names[i].x == 'L' || atom_names[i].x == 'E') && atom_names[i].y == 'P') ||
        (atom_names[i].x == ' ' && (atom_names[i].y == 'L' || atom_names[i].y == 'E') &&
         atom_names[i].z == 'P') ||
        (atom_names[i].x == ' ' && atom_names[i].y == ' ' &&
         (atom_names[i].z == 'L' || atom_names[i].z == 'E') && atom_names[i].w == 'P')) {
      z_numbers[i] = 0;
      n_named_vsite++;
    }
  }
  if (n_named_vsite > 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      rtWarn(std::to_string(n_named_vsite) + " atoms were taken to be virtual sites due to names "
             "beginning with \"LP\" or \"EP\", in spite of having non-zero mass.",
             "atomicNumbersFromMasses");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (minValue(z_numbers) >= 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      rtWarn("Unable to find atomic numbers in topology source file " + source + ".  Z-numbers "
             "were inferred from atomic masses by comparing to averages based on the natural "
             "isotopic abundance for each element.", "atomicNumbersFromMasses");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return z_numbers;
  }

  // Are there any atoms that look like they may be hydrogen which can return repartitioned mass
  // to their parent atoms?
  std::vector<int> bond_counts(natom, 0);
  std::vector<double> reconstituted_masses(masses.begin(), masses.end());
  for (int i = 0; i < natom; i++) {

    // Even when hydrogen mass repartitioning is in effect, many atoms' masses will not have been
    // repartitioned (i.e. water).  Search for the lightest particle that looks like a bound
    // hydrogen but does not have the natural mass of a hydrogen.
    if (z_numbers[i] == 1) {
      continue;
    }
    for (int j = nb12_bounds[i]; j < nb12_bounds[i + 1]; j++) {
      bond_counts[i] += (masses[nb12_list[j]] > constants::tiny);
    }
    
    // Hydrogen mass repartitioning could conceivably go as high as 6.0, which is less than the
    // mass of any natural Lithium atom and would involve a more than 50% repartitioning of mass
    // for one hydrogen atom bound to an oxygen.  Helium does not form any bonds.
    if (bond_counts[i] == 1 && masses[i] < 6.0 && masses[i] < masses[nb12_list[nb12_bounds[i]]]) {
      const double dmass = masses[i] - chemistry::elemental_masses[1];
      reconstituted_masses[i] -= dmass;
      reconstituted_masses[nb12_list[nb12_bounds[i]]] += dmass;
    }
  }
  z_numbers = massToZNumber(reconstituted_masses);
  int n_unknown = 0;
  for (int i = 0; i < natom; i++) {
    n_unknown += z_numbers[i] < 0;
  }
  if (n_unknown == 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      rtWarn("Unable to find atomic numbers in topology source file " + source + ".  Taking all "
             "unidentified atoms with one bond to be hydrogens with mass donated from their heavy "
             "atom parents.  Atomic numbers have been recovered.", "atomicNumbersFromMasses");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  else {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Unable to find atomic numbers in topology source file " + source + ".  No candidate "
            "for hydrogen in a hydrogen mass repartitioning scheme was identified for " +
            std::to_string(n_unknown) + " atoms.", "atomicNumbersFromMasses");
      break;
    case ExceptionResponse::WARN:
      rtWarn("Unable to find atomic numbers in topology source file " + source + ".  No candidate "
             "for hydrogen in a hydrogen mass repartitioning scheme was identified for " +
             std::to_string(n_unknown) + " atoms.", "atomicNumbersFromMasses");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  return z_numbers;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> traceBondedPatterns(const Map1234 &all_nb_excl) {

  // Get the atom count from an array length (the array is padded
  // by one to store the upper bound of the final atom)
  const int natom = all_nb_excl.nb11_excl_bounds.size() - 1;
  std::vector<int> result(natom, -1);
  std::vector<bool> placed(natom, false);
  std::vector<int> current_mol(natom, 0);
  int molid = 0;
  for (int i = 0; i < natom; i++) {
    if (placed[i]) {
      continue;
    }

    // The current number of atoms in the molecule is 1: the ith atom.  Find the rest.
    placed[i] = true;
    current_mol[0] = i;
    int ncurr = 1;
    int pivot = 0;
    while (pivot < ncurr) {
      const int pivot_atom = current_mol[pivot];
      for (int j = all_nb_excl.nb11_excl_bounds[pivot_atom];
           j < all_nb_excl.nb11_excl_bounds[pivot_atom + 1]; j++) {
        if (placed[all_nb_excl.nb11_excl_list[j]]) {
          continue;
        }
        placed[all_nb_excl.nb11_excl_list[j]] = true;
        current_mol[ncurr] = all_nb_excl.nb11_excl_list[j];
        ncurr++;
      }
      for (int j = all_nb_excl.nb12_excl_bounds[pivot_atom];
           j < all_nb_excl.nb12_excl_bounds[pivot_atom + 1]; j++) {
        if (placed[all_nb_excl.nb12_excl_list[j]]) {
          continue;
        }
        placed[all_nb_excl.nb12_excl_list[j]] = true;
        current_mol[ncurr] = all_nb_excl.nb12_excl_list[j];
        ncurr++;
      }
      for (int j = all_nb_excl.nb13_excl_bounds[pivot_atom];
           j < all_nb_excl.nb13_excl_bounds[pivot_atom + 1]; j++) {
        if (placed[all_nb_excl.nb13_excl_list[j]]) {
          continue;
        }
        placed[all_nb_excl.nb13_excl_list[j]] = true;
        current_mol[ncurr] = all_nb_excl.nb13_excl_list[j];
        ncurr++;
      }
      for (int j = all_nb_excl.nb14_excl_bounds[pivot_atom];
           j < all_nb_excl.nb14_excl_bounds[pivot_atom + 1]; j++) {
        if (placed[all_nb_excl.nb14_excl_list[j]]) {
          continue;
        }
        placed[all_nb_excl.nb14_excl_list[j]] = true;
        current_mol[ncurr] = all_nb_excl.nb14_excl_list[j];
        ncurr++;
      }
      pivot++;
    }

    // Log this molecule in the result
    for (int j = 0; j < ncurr; j++) {
      result[current_mol[j]] = molid;
    }
    molid++;
  }

  // Check that all atoms were covered
  for (int i = 0; i < natom; i++) {
    if (result[i] == -1) {
      rtErr("Atom " + std::to_string(i + 1) + " was not placed into any molecule.");
    }
  }
  
  return result;
}

//-------------------------------------------------------------------------------------------------
void mapMolecules(const int atom_count, int *molecule_count, const Map1234 &all_nb_excl,
                  std::vector<int> *molecule_membership, std::vector<int> *molecule_limits,
                  std::vector<int> *molecule_contents) {
  *molecule_membership = traceBondedPatterns(all_nb_excl);
  const int local_molecule_count = stmath::maxValue(*molecule_membership) + 1;
  molecule_limits->resize(local_molecule_count + 1, 0);
  molecule_contents->resize(atom_count, -1);
  int* tml_data = molecule_limits->data();
  int* tmm_data = molecule_membership->data();
  int* tmc_data = molecule_contents->data();
  for (int i = 0; i < atom_count; i++) {
    tml_data[tmm_data[i]] += 1;
  }
  prefixSumInPlace(molecule_limits, PrefixSumType::EXCLUSIVE, "mapMolecules");
  for (int i = 0; i < atom_count; i++) {
    const int mol_idx = tmm_data[i];
    tmc_data[tml_data[mol_idx]] = i;
    tml_data[mol_idx] += 1;
  }
  for (int i = local_molecule_count; i > 0; i--) {
    tml_data[i] = tml_data[i - 1];
  }
  tml_data[0] = 0;
  *molecule_count = local_molecule_count;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> cubicSplineDerivativeStencil(const int npts) {

  // This procedure is best done with row-major matrices to allow partial pivoting Gaussian
  // elimination in the square matrix A (of rank 4 x npts) to manipulate the corresponding rows 
  // of the tracking matrix T.  Construct vectors of double pointers to lay out "maps" into
  // separate vectors of actual data, so that row swapping can proceed accordingly.
  const int rank = 4 * npts;
  std::vector<double> a_matrix(rank * rank, 0.0);
  std::vector<double*> a_map(rank);
  std::vector<double> t_matrix(rank * npts, 0.0);
  std::vector<double*> t_map(rank);
  double* adata = a_matrix.data();
  double* tdata = t_matrix.data();
  for (int i = 0; i < rank; i++) { 
    a_map[i] = &adata[i * rank];
    t_map[i] = &tdata[i * npts];
  }
  const double dnpts = static_cast<double>(npts);
  const double stepsize = 2.0 * symbols::pi / dnpts;
  for (int i = 0; i < npts; i++) {
    const double phi = 0.0;
    const double phip = phi + stepsize;
    int ip = i + 1;
    if (ip == npts) {
      ip = 0;
    }

    /// Energy values for the piecewise cubic spline must meet
    a_map[(4 * i)    ][(4 * i)    ] = phi * phi * phi;
    a_map[(4 * i)    ][(4 * i) + 1] = phi * phi;
    a_map[(4 * i)    ][(4 * i) + 2] = phi;
    a_map[(4 * i)    ][(4 * i) + 3] = 1.0;
    t_map[(4 * i)    ][i          ] = 1.0;
    a_map[(4 * i) + 1][(4 * i)    ] = phip * phip * phip;
    a_map[(4 * i) + 1][(4 * i) + 1] = phip * phip;
    a_map[(4 * i) + 1][(4 * i) + 2] = phip;
    a_map[(4 * i) + 1][(4 * i) + 3] = 1.0;
    t_map[(4 * i) + 1][ip         ] = 1.0;

    // First derivatives of each spline segment must
    // be equal at the boundary between two segments
    a_map[(4 * i) + 2][(4 * i)     ] =  3.0 * phip * phip;
    a_map[(4 * i) + 2][(4 * i)  + 1] =  2.0 * phip;
    a_map[(4 * i) + 2][(4 * i)  + 2] =  1.0;
    a_map[(4 * i) + 2][(4 * ip)    ] = -3.0 * phi * phi;
    a_map[(4 * i) + 2][(4 * ip) + 1] = -2.0 * phi;
    a_map[(4 * i) + 2][(4 * ip) + 2] = -1.0;

    // Second derivatives of each spline segment must
    // be equal at the boundary between two segments,
    // despite the fact that second derivatives don't
    // carry over into the bicubic spline table.
    a_map[(4 * i) + 3][(4 * i)     ] =  6.0 * phip;
    a_map[(4 * i) + 3][(4 * i)  + 1] =  2.0;
    a_map[(4 * i) + 3][(4 * ip)    ] = -6.0 * phi;
    a_map[(4 * i) + 3][(4 * ip) + 1] = -2.0;
  }

  // Proceed with Gaussian elimination, make a lower triangular matrix
  for (int i = 0; i < rank; i++) {

    // Find the row with the largest (absolute) non-zero value in this column
    double maxval = fabs(a_map[i][i]);
    int pivot = i;
    for (int j = i + 1; j < rank; j++) {
      if (fabs(a_map[j][i]) > maxval) {
        maxval = fabs(a_map[j][i]);
        pivot = j;
      }
    }

    // Make this the ith row
    if (pivot != i) {
      std::swap(a_map[pivot], a_map[i]);
      std::swap(t_map[pivot], t_map[i]);
    }

    // Kill everything below the ith row
    for (int j = i + 1; j < rank; j++) {
      if (fabs(a_map[j][i]) > 0.0) {
        double factor = a_map[j][i] / a_map[i][i];
        for (int k = i; k < rank; k++) {
          a_map[j][k] -= factor * a_map[i][k];
        }
        for (int k = 0; k < npts; k++) {
          t_map[j][k] -= factor * t_map[i][k];
        }
      }
    }
  }

  // Make a diagonal matrix
  for (int i = rank - 1; i >= 0; i--) {

    // Normalize each row
    double factor = 1.0 / a_map[i][i];
    for (int j = 0; j < npts; j++) {
      t_map[i][j] *= factor;
    }
    a_map[i][i] = 1.0;

    // Use it to kill everything above
    for (int j = i - 1; j >= 0; j--) {
      factor = a_map[j][i];
      a_map[j][i] -= factor * a_map[i][i];
      for (int k = 0; k < npts; k++) {
        t_map[j][k] -= factor * t_map[i][k];
      }
    }
  }

  // Average the formulas, which are found in rows 2, 6, ..., 4 * i + 2 (indices starting from
  // zero) of the tracking matrix (t_matrix, t_map) and right-shifted with every increment of i.
  // Adjust the formula based on the units of each grid spacing (2 pi / npts). 
  std::vector<double> result(npts, 0.0);
  for (int i = 0; i < npts; i++) {
    for (int rowcon = 0; rowcon < npts; rowcon++) {
      const int colcon = i + rowcon - ((i + rowcon >= npts) * npts);
      result[i] += t_map[(4 * rowcon) + 2][colcon];
    }
    result[i] /= dnpts;
    result[i] *= stepsize;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
CmapAccessories computeCmapDerivatives(const int cmap_surf_count,
                                       const std::vector<int> tmp_cmap_surf_dims,
                                       const std::vector<int> tmp_cmap_surf_bounds,
                                       const std::vector<double> tmp_cmap_surfaces) {
  CmapAccessories result;

  // Allocate data for the derivative arrays
  result.phi_derivatives.resize(tmp_cmap_surfaces.size());
  result.psi_derivatives.resize(tmp_cmap_surfaces.size());
  result.phi_psi_derivatives.resize(tmp_cmap_surfaces.size());

  // Allocate data and set bounds for the patch matrix form of each CMAP
  result.patch_matrix_bounds.resize(cmap_surf_count + 1);
  if (static_cast<size_t>(cmap_surf_count) > tmp_cmap_surf_dims.size()) {
    rtErr("Topology specifies " + std::to_string(cmap_surf_count) + " CMAP surfaces but provides "
          "dimensions for only " + std::to_string(tmp_cmap_surf_dims.size()) + ".",
          "computeCmapDerivatives");
  }
  int boundary = 0;
  for (int i = 0; i < cmap_surf_count; i++) {
    result.patch_matrix_bounds[i] = boundary;
    boundary += tmp_cmap_surf_dims[i] * tmp_cmap_surf_dims[i] * 16;
  }
  result.patch_matrix_bounds[cmap_surf_count] = boundary;
  result.patch_matrix_form.resize(boundary);

  // Transformation matrices for converting potentials and derivatives to the bicubic spline
  // coefficients.  Qt * U * Q = A, where U is the block matrix of CMAP values and derivatives,
  // Q and Qt are given below, and A is the spline coefficients from which the potential and
  // derivative at an arbitrary point can be obtained.
  const std::vector<double> q  = {  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  0.0,
                                   -3.0,  3.0, -2.0, -1.0,  2.0, -2.0,  1.0,  1.0  };
  const std::vector<double> qt = {  1.0,  0.0, -3.0,  2.0,  0.0,  0.0,  3.0, -2.0,
                                    0.0,  1.0, -2.0,  1.0,  0.0,  0.0, -1.0,  1.0  };
  std::vector<double> qt_u(16, 0.0);

  // Loop over all CMAPs and compute phi / psi derivatives, plus the cross derivative.
  int current_stencil_size = -1;
  std::vector<double> stencil;
  std::vector<double> patch_buffer(16);
  for (int i = 0; i < cmap_surf_count; i++) {
    
    // Compute the stencil for a CMAP of this dimension
    if (current_stencil_size != tmp_cmap_surf_dims[i]) {
      stencil = cubicSplineDerivativeStencil(tmp_cmap_surf_dims[i]);
      current_stencil_size = tmp_cmap_surf_dims[i];
    }
    const double pt_spacing = static_cast<double>(current_stencil_size) / symbols::twopi;

    // Use the stencil to determine the first direction (phi) derivatives of the map
    const int offset = tmp_cmap_surf_bounds[i];
    for (int j = 0; j < current_stencil_size; j++) {
      const int offset_j = offset + j;
      for (int k = 0; k < current_stencil_size; k++) {
        const int k_col = k * current_stencil_size;

        // Evaluate the stencil for point (j, k) of the column-major format surface matrix
        double deriv_sum = 0.0;
        for (int m = 0; m < current_stencil_size; m++) {
          const int pos = j + m - ((j + m >= current_stencil_size) * current_stencil_size);
          deriv_sum += stencil[m] * tmp_cmap_surfaces[offset + k_col + pos];
        }
        result.phi_derivatives[offset_j + k_col] = deriv_sum;
      }
    }

    // Use the stencil to determine the second direction (psi) derivatives of the map
    for (int j = 0; j < current_stencil_size; j++) {
      const int offset_j = offset + j;
      for (int k = 0; k < current_stencil_size; k++) {

        // Evaluate the stencil for point (j, k) of the column-major format surface matrix
        double deriv_sum = 0.0;
        for (int m = 0; m < current_stencil_size; m++) {
          const int pos = k + m - ((k + m >= current_stencil_size) * current_stencil_size);
          deriv_sum += stencil[m] * tmp_cmap_surfaces[(pos * current_stencil_size) + offset_j];
        }
        result.psi_derivatives[offset_j + (k * current_stencil_size)] = deriv_sum;
      }
    }

    // Finally, determine the cross derivatives (phi and psi) using the stencil and either of
    // the existing first derivative maps.
    for (int j = 0; j < current_stencil_size; j++) {
      const int offset_j = offset + j;
      for (int k = 0; k < current_stencil_size; k++) {
        const int k_col = k * current_stencil_size;

        // Evaluate the stencil for point (j, k) of the column-major format surface matrix
        double deriv_sum = 0.0;
        for (int m = 0; m < current_stencil_size; m++) {
          const int mpos = j + m - ((j + m >= current_stencil_size) * current_stencil_size);
          for (int n = 0; n < current_stencil_size; n++) {
            const int npos = k + n - ((k + n >= current_stencil_size) * current_stencil_size);
            deriv_sum += stencil[m] * stencil[n] *
                         tmp_cmap_surfaces[offset + (npos * current_stencil_size) + mpos];
          }
        }
        result.phi_psi_derivatives[offset_j + k_col] = deriv_sum;
      }
    }

    // Collate the potential, first derivatives in phi and psi, and cross derivative maps into
    // a patch matrix.  This matrix is 16 times the size of any of the above grids but collects
    // all the necessary information for a contiguous access.
    const int patch_offset = result.patch_matrix_bounds[i];
    for (int j = 0; j < current_stencil_size; j++) {
      for (int k = 0; k < current_stencil_size; k++) {
        const int jk_offset   = offset + (k * current_stencil_size) + j;
        const int jpk_offset  = offset + (k * current_stencil_size) +
                                ((j < current_stencil_size - 1) * (j + 1));
        const int jkp_offset  = offset + j +
                                ((k < current_stencil_size - 1) * (k + 1) * current_stencil_size);
        const int jpkp_offset = offset + ((j < current_stencil_size - 1) * (j + 1)) +
                                ((k < current_stencil_size - 1) * (k + 1) * current_stencil_size);

        // Fill out the patch buffer with the surface and derivative matrix: a concatenated,
        // column-major format rank 4 matrix for each CMAP aggregated's patch representation.
        // Patch (j, k) will become a collection of coefficients for the surface (S), but to get
        // there a matrix must be constructed as follows: the surface (S), first dimension
        // derivative (dS/dphi), second dimension first derivative (dS/dpsi), and cross derivative
        // (d2S/dphi-dpsi).
        //
        // [        S(j  ,k  )         S(j  ,k+1)      (dS/dpsi)(j  ,k  )      (dS/dpsi)(j  ,k+1)
        //          S(j+1,k  )         S(j+1,k+1)      (dS/dpsi)(j+1,k  )      (dS/dpsi)(j+1,k+1)
        //  (dS/dphi)(j  ,k  ) (dS/dphi)(j  ,k+1) (dS/dphi-dpsi)(j  ,k  ) (dS/dphi-dpsi)(j  ,k+1)
        //  (dS/dphi)(j+1,k  ) (dS/dphi)(j+1,k+1) (dS/dphi-dpsi)(j+1,k  ) (dS/dphi-dpsi)(j+1,k+1) ]
        patch_buffer[ 0] = tmp_cmap_surfaces[jk_offset];
        patch_buffer[ 1] = tmp_cmap_surfaces[jpk_offset];
        patch_buffer[ 2] = result.phi_derivatives[jk_offset];
        patch_buffer[ 3] = result.phi_derivatives[jpk_offset];
        patch_buffer[ 4] = tmp_cmap_surfaces[jkp_offset];
        patch_buffer[ 5] = tmp_cmap_surfaces[jpkp_offset];
        patch_buffer[ 6] = result.phi_derivatives[jkp_offset];
        patch_buffer[ 7] = result.phi_derivatives[jpkp_offset];
        patch_buffer[ 8] = result.psi_derivatives[jk_offset];
        patch_buffer[ 9] = result.psi_derivatives[jpk_offset];
        patch_buffer[10] = result.phi_psi_derivatives[jk_offset];
        patch_buffer[11] = result.phi_psi_derivatives[jpk_offset];
        patch_buffer[12] = result.psi_derivatives[jkp_offset];
        patch_buffer[13] = result.psi_derivatives[jpkp_offset];
        patch_buffer[14] = result.phi_psi_derivatives[jkp_offset];
        patch_buffer[15] = result.phi_psi_derivatives[jpkp_offset];

        // Convert the initial matrix into the matrix of coefficients: [ a00 a01 ... a33 ].
        // This does a lot of pre-computation for the same overall amount of memory reads, and has
        // the added nenefit of doing that portion of the CMAP arithmetic in double precision
        // before converting the coefficients matrix [ a00 a01 ... a33 ] to single precision, the
        // form that will ultimately be used to evaluate the CMAPs.
        matrixMultiply(qt.data(), 4, 4, patch_buffer.data(), 4, 4, qt_u.data(), 1.0, 1.0, 0.0);
        matrixMultiply(qt_u.data(), 4, 4, q.data(), 4, 4, patch_buffer.data(), 1.0, 1.0, 0.0);

        // Graft the buffer into the main array
        const int jk16 = ((k * current_stencil_size) + j) * 16;
        for (int m = 0; m < 16; m++) {
          result.patch_matrix_form[patch_offset + jk16 + m] = patch_buffer[m];
        }
      }
    }
  }

  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int3> checkDihedral14Coverage(const int atom_count,
                                          const std::vector<int> &atomic_numbers,
                                          const BasicValenceTable &bvt, const Map1234 &all_nb_excl,
                                          const VirtualSiteTable &vs_tbl,
                                          const AttenuationParameterSet &attn_parm,
                                          const ExceptionResponse policy) {

  // Initialize the result to an empty vector.  If no 1:4 interactions are left without a dihedral
  // to address them, the empty vector will indicate that successful result.
  std::vector<int3> result;

  // Create an array of atoms and bounds which indicates the dihedrals to which each atom is either
  // the I or L atom.
  std::vector<int> dihe_participation_bounds(atom_count + 1, 0);
  std::vector<int> dihe_participation(2 * bvt.total_dihes);
  for (int pos = 0; pos < bvt.total_dihes; pos++) {
    dihe_participation_bounds[bvt.dihe_i_atoms[pos]] += 1;
    dihe_participation_bounds[bvt.dihe_l_atoms[pos]] += 1;
  }
  prefixSumInPlace(&dihe_participation_bounds, PrefixSumType::EXCLUSIVE,
                          "checkDihedral14Coverage");
  for (int pos = 0; pos < bvt.total_dihes; pos++) {
    const int i_atom = bvt.dihe_i_atoms[pos];
    const int l_atom = bvt.dihe_l_atoms[pos];
    dihe_participation[dihe_participation_bounds[i_atom]] = pos;
    dihe_participation[dihe_participation_bounds[l_atom]] = pos;
    dihe_participation_bounds[i_atom] += 1;
    dihe_participation_bounds[l_atom] += 1;
  }
  for (int i = atom_count; i > 0; i--) {
    dihe_participation_bounds[i] = dihe_participation_bounds[i - 1];
  }
  dihe_participation_bounds[0] = 0;

  // Loop over all atoms, find its 1:4 interactions, and, if the 1:4 interaction goes to an
  // atom of greater index (to avoid double-counting), figure out which dihedral might be able
  // to cover it.
  for (int i = 0; i < atom_count; i++) {
    for (int j = all_nb_excl.nb14_excl_bounds[i]; j < all_nb_excl.nb14_excl_bounds[i + 1]; j++) {
      const int jatom = all_nb_excl.nb14_excl_list[j]; 
      if (jatom < i) {
        continue;
      }

      // The 1:4 interaction of atoms i and jatom should be covered by some dihedral, ideally.
      // Seek that dihedral out by looking through the lists of dihedrals for either atom.
      bool found = false;
      for (int k = dihe_participation_bounds[i]; k < dihe_participation_bounds[i + 1]; k++) {
        const int dihe_index = dihe_participation[k];
        if ((bvt.dihe_i_atoms[dihe_index] == i && bvt.dihe_l_atoms[dihe_index] == jatom) ||
            (bvt.dihe_l_atoms[dihe_index] == i && bvt.dihe_i_atoms[dihe_index] == jatom)) {
          const TorsionKind trkind = static_cast<TorsionKind>(bvt.dihe_mods[dihe_index].w);
          if (trkind == TorsionKind::PROPER) {

            // If this 1:4 interaction is already covered by a dihedral, handle the exception
            if (found) {
              switch (policy) {
              case ExceptionResponse::DIE:
                rtErr("A 1:4 attenuated interaction between atoms with indices " +
                      std::to_string(i + 1) + " and " + std::to_string(jatom + 1) + " is covered "
                      "by more than one dihedral.", "checkDihedral14Coverage");
              case ExceptionResponse::WARN:
                rtWarn("A 1:4 attenuated interaction between atoms with indices " +
                       std::to_string(i + 1) + " and " + std::to_string(jatom + 1) + " is covered "
                       "by more than one dihedral.", "checkDihedral14Coverage");
                break;
              case ExceptionResponse::SILENT:
                break;
              }
            }
            found = true;
          }
        }
      }

      // If the 1:4 interaction was not found to be covered by some dihedral, it will not get
      // picked up as part of the loop over dihedrals.  Some additional, explicit 1:4 interaction
      // with accompanying Lennard-Jones and electrostatic scaling will have to be prepared.  To
      // do pair interactions with at least three pieces of input (two parameters and a bit-packed
      // pair of indices) is the height of memory-boundness, which is why parameter tables
      // detailing the various interactions by index will be cached.  These extraneous 1:4
      // attenuated interactions will then be feasible with just a single 4-byte word of input.
      if (found == false) {

        // Initialize the interaction to have zero scaling factors.  If nothing is found to
        // inform the interaction otherwise, these coefficients will remain at zero.
        int3 atp;
        atp.x = i;
        atp.y = jatom;

        // If this interaction involves a virtual site or pair of virtual sites, check for
        // dihedrals affecting the parent atom or atoms
        if (atomic_numbers[i] == 0 || atomic_numbers[jatom] == 0) {
          const int actual_iatom = (atomic_numbers[i] == 0) ?
                                   vs_tbl.frame1_atoms[vs_tbl.vs_numbers[i]] : i;
          const int actual_jatom = (atomic_numbers[jatom] == 0) ?
                                   vs_tbl.frame1_atoms[vs_tbl.vs_numbers[jatom]] : jatom;
          for (int k = dihe_participation_bounds[actual_iatom];
               k < dihe_participation_bounds[actual_iatom + 1]; k++) {
            const int dihe_index = dihe_participation[k];
            if ((bvt.dihe_i_atoms[dihe_index] == actual_iatom &&
                 bvt.dihe_l_atoms[dihe_index] == actual_jatom) ||
                (bvt.dihe_l_atoms[dihe_index] == actual_iatom &&
                 bvt.dihe_i_atoms[dihe_index] == actual_jatom)) {
              const TorsionKind trkind = static_cast<TorsionKind>(bvt.dihe_mods[dihe_index].w);
              if (trkind == TorsionKind::PROPER) {

                // This dihedral becomes the latest candidate to supply the attenuated interaction
                // parameters for this exclusion.  If there are more than one dihedral which
                // could supply parameters, this case would have been trapped above.  If the
                // policy to to merely warn, or even be silent, then this 1:4 term will likewise
                // take parameters from the last encountered dihedral that could do the job.
                found = true;
                const int dihe_param_idx = bvt.dihe_param_idx[dihe_index];
                atp.z = attn_parm.dihe14_parameter_indices[dihe_param_idx];
              }
            }
          }
        }

        // Commit this result to the list
        result.push_back(atp);
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int reviewLargestResidue(const std::vector<int> residue_limits, const int current_lr_size,
                         const ExceptionResponse policy) {
  int max_atoms = 0;
  const size_t nres = residue_limits.size() - 1ULL;
  for (size_t i = 0; i < nres; i++) {
    max_atoms = std::max(max_atoms, residue_limits[i + 1] - residue_limits[i]);
  }
  if (max_atoms != current_lr_size) {
    const std::string msg("The largest residue size recorded in the preamble (" +
                          std::to_string(current_lr_size) + ") does not match the residue "
                          "limits (" + std::to_string(max_atoms) + ").");
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(msg, "reviewLargestResidue");
    case ExceptionResponse::WARN:
      rtWarn(msg, "reviewLargestResidue");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  return max_atoms;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> matchExtendedName(const char4* overflow_names, const int n_overflow,
                                   const std::string &query,
                                   const std::vector<WildCardKind> &wildcards) {
  const int n_wild = wildcards.size();
  std::vector<WildCardKind> local_wildcards =
    (wildcards.size() == query.size()) ? wildcards :
                                         std::vector<WildCardKind>(query.size(),
                                                                   WildCardKind::NONE);
  std::vector<int> result;
  for (int i = 0; i < n_overflow; i++) {
    bool match = true;

    // The entire overflow name must be accommodated by the query and its wildcards
    std::string ofl_name(12, ' ');
    for (int j = 0; j < 3; j++) {
      ofl_name[4*j    ] = overflow_names[(4 * i) + j].x;
      ofl_name[4*j + 1] = overflow_names[(4 * i) + j].y;
      ofl_name[4*j + 2] = overflow_names[(4 * i) + j].z;
      ofl_name[4*j + 3] = overflow_names[(4 * i) + j].w;
    }
    if (strcmpWildCard(ofl_name, query, local_wildcards)) {
      result.push_back(i);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
SettleParm getSettleParameters(const int ox_idx, const int h1_idx, const int h2_idx,
                               const std::vector<double> &atomic_mass,
                               const BasicValenceTable &bvt,
                               const std::vector<double> &bond_equilibria,
                               const std::vector<double> &angl_equilibria) {

  // Find the mass of the oxygen and of the two hydrogens
  const double ox_mass = atomic_mass[ox_idx];
  const double hd_mass = atomic_mass[h1_idx];
  if (fabs(atomic_mass[h1_idx] - atomic_mass[h2_idx]) > constants::small) {
    rtErr("Masses of atoms " + std::to_string(h1_idx) + " and " + std::to_string(h2_idx) +
          "must agree in order to apply the SETTLE algorithm.  Their masses come to " +
          realToString(atomic_mass[h1_idx], 4, NumberFormat::STANDARD_REAL) + " and " +
          realToString(atomic_mass[h2_idx], 4, NumberFormat::STANDARD_REAL) + ".",
          "getSettleParameters");
  }

  // Seek out the H-H and O-H distances based on the known equilibrium bond and angle parameters
  double oh1_dist, oh2_dist, hh_dist;
  bool oh1_found = false;
  bool oh2_found = false;
  bool hh_found = false;
  for (int i = bvt.bond_assigned_bounds[ox_idx]; i < bvt.bond_assigned_bounds[ox_idx + 1]; i++) {
    if (bvt.bond_assigned_atoms[i] == h1_idx) {
      oh1_found = true;
      oh1_dist = bond_equilibria[bvt.bond_param_idx[bvt.bond_assigned_terms[i]]];
    }
    else if (bvt.bond_assigned_atoms[i] == h2_idx) {
      oh2_found = true;
      oh2_dist = bond_equilibria[bvt.bond_param_idx[bvt.bond_assigned_terms[i]]];
    }
  }
  for (int i = bvt.bond_assigned_bounds[h1_idx]; i < bvt.bond_assigned_bounds[h1_idx + 1]; i++) {
    if (bvt.bond_assigned_atoms[i] == ox_idx) {
      oh1_found = true;
      oh1_dist = bond_equilibria[bvt.bond_param_idx[bvt.bond_assigned_terms[i]]];
    }
    else if (bvt.bond_assigned_atoms[i] == h2_idx) {
      hh_found = true;
      hh_dist = bond_equilibria[bvt.bond_param_idx[bvt.bond_assigned_terms[i]]];
    }    
  }
  for (int i = bvt.bond_assigned_bounds[h2_idx]; i < bvt.bond_assigned_bounds[h2_idx + 1]; i++) {
    if (bvt.bond_assigned_atoms[i] == h1_idx) {
      hh_found = true;
      hh_dist = bond_equilibria[bvt.bond_param_idx[bvt.bond_assigned_terms[i]]];
    }
    else if (bvt.bond_assigned_atoms[i] == ox_idx) {
      oh2_found = true;
      oh2_dist = bond_equilibria[bvt.bond_param_idx[bvt.bond_assigned_terms[i]]];
    }
  }
  if (oh1_found == false || oh2_found == false) {
    rtErr("Bond parameters were not found to link oxygen atom " + std::to_string(ox_idx) +
          " with hydrogens " + std::to_string(h1_idx) + " and " + std::to_string(h2_idx) + ".",
          "getSettleParameters");
  }
  if (hh_found == false) {

    // Search for an H-O-H angle in order to determine the H-H equilibrium distance.
    bool hoh_found = false;
    double hoh_theta;
    for (int i = bvt.angl_assigned_bounds[ox_idx]; i < bvt.angl_assigned_bounds[ox_idx + 1]; i++) {
      if ((bvt.angl_assigned_atoms[(2 * i)    ] == h1_idx &&
           bvt.angl_assigned_atoms[(2 * i) + 1] == h2_idx) ||
          (bvt.angl_assigned_atoms[(2 * i)    ] == h2_idx &&
           bvt.angl_assigned_atoms[(2 * i) + 1] == h1_idx)) {
        hoh_found = true;
        hoh_theta = angl_equilibria[bvt.angl_param_idx[bvt.angl_assigned_terms[i]]];
      }
    }
    if (hoh_found) {
      hh_dist = 2.0 * oh1_dist * sin(0.5 * hoh_theta);
    }
    else {
      rtErr("No angle parameter was found between oxygen atom " + std::to_string(ox_idx) +
            " and hydrogen atoms " + std::to_string(h1_idx) + " and " + std::to_string(h2_idx) +
            ".  In lieu of a bond between the hydrogens, this angle is required to determine an "
            "equilibrium geometry.", "getSettleParameters");
    }
  }

  // With O-H1, O-H2, and H-H distances known, construct the equilibrium geometry (with the oxygen
  // positioned at the origin and the H-O-H bisector along the x axis with the molecule in the
  // xy plane in order to determine a center of geometry and then the critical distances of the
  // oxygen and each hydrogen to the center of mass.
  if (fabs(oh1_dist - oh2_dist) > constants::small) {
    rtErr("The distance between oxygen " + std::to_string(ox_idx) + " and hydrogens " +
          std::to_string(h1_idx) + " and " + std::to_string(h2_idx) + " must be identical in "
          "order to apply SETTLE.", "getSettleParameters");
  }
  const double inv_total_mass = 1.0 / (ox_mass + (2.0 * hd_mass));
  const double t1 = 0.5 * ox_mass / hd_mass;
  const double rc = 0.5 * hh_dist;
  const double ra = sqrt((oh1_dist * oh1_dist) - (rc * rc)) / (1.0 + t1);
  return { ox_mass * inv_total_mass, hd_mass * inv_total_mass, ra, t1 * ra, rc, 1.0 / ra };
}

} // namespace topology
} // namespace stormm
