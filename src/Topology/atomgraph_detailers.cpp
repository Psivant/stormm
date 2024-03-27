#include <cmath>
#include <cstdio>
#include <climits>
#include "copyright.h"
#include "Math/rounding.h"
#include "Math/vector_ops.h"
#include "Parsing/ascii_numbers.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "amber_prmtop_util.h"
#include "atomgraph.h"

namespace stormm {
namespace topology {

using stmath::addScalarToVector;
using stmath::roundUp;
using parse::NumberFormat;
using parse::PolyNumeric;
using parse::TextFile;
using parse::TextFileReader;
using parse::TextOrigin;
using parse::verifyNumberFormat;

//-------------------------------------------------------------------------------------------------
void AtomGraph::loadHybridArrays(const std::vector<int> &tmp_desc,
                                 const std::vector<int> &tmp_residue_limits,
                                 const std::vector<int> &tmp_atom_struc_numbers,
                                 const std::vector<int> &tmp_residue_numbers,
                                 const std::vector<int> &tmp_molecule_limits,
                                 const std::vector<int> &tmp_atomic_numbers,
                                 const std::vector<int> &tmp_molecule_membership,
                                 const std::vector<int> &tmp_mobile_atoms,
                                 const std::vector<int> &tmp_molecule_contents,
                                 const std::vector<int> &tmp_cmap_surf_dims,
                                 const std::vector<int> &tmp_cmap_surf_bounds,
                                 const std::vector<int> &tmp_charge_indices,
                                 const std::vector<int> &tmp_lennard_jones_indices,
                                 const std::vector<int> &tmp_inferred_14_i_atoms,
                                 const std::vector<int> &tmp_inferred_14_l_atoms,
                                 const std::vector<int> &tmp_inferred_14_param_idx,
                                 const std::vector<int> &tmp_neck_gb_indices,
                                 const std::vector<int> &tmp_tree_joining_info,
                                 const std::vector<int> &tmp_last_rotator_info,
                                 const std::vector<double> &tmp_charges,
                                 const std::vector<double> &tmp_masses,
                                 const std::vector<double> &tmp_ub_stiffnesses,
                                 const std::vector<double> &tmp_ub_equilibria,
                                 const std::vector<double> &tmp_charmm_impr_stiffnesses,
                                 const std::vector<double> &tmp_charmm_impr_phase_angles,
                                 const std::vector<double> &tmp_cmap_surfaces,
                                 const std::vector<double> &tmp_bond_stiffnesses,
                                 const std::vector<double> &tmp_bond_equilibria,
                                 const std::vector<double> &tmp_angl_stiffnesses,
                                 const std::vector<double> &tmp_angl_equilibria,
                                 const std::vector<double> &tmp_dihe_amplitudes,
                                 const std::vector<double> &tmp_dihe_periodicities,
                                 const std::vector<double> &tmp_dihe_phase_angles,
                                 const std::vector<double> &tmp_charge_parameters,
                                 const std::vector<double> &tmp_lj_a_values,
                                 const std::vector<double> &tmp_lj_b_values,
                                 const std::vector<double> &tmp_lj_c_values,
                                 const std::vector<double> &tmp_lj_14_a_values,
                                 const std::vector<double> &tmp_lj_14_b_values,
                                 const std::vector<double> &tmp_lj_14_c_values,
                                 const std::vector<double> &tmp_atomic_pb_radii,
                                 const std::vector<double> &tmp_gb_screening_factors,
                                 const std::vector<double> &tmp_gb_coef,
                                 const std::vector<double> &tmp_solty_info,
                                 const std::vector<double> &tmp_hbond_a_values,
                                 const std::vector<double> &tmp_hbond_b_values,
                                 const std::vector<double> &tmp_hbond_cutoffs,
                                 const std::vector<char4> &tmp_atom_names,
                                 const std::vector<char4> &tmp_atom_types,
                                 const std::vector<char4> &tmp_residue_names,
                                 const std::vector<char4> &tmp_tree_symbols,
                                 const CmapAccessories &cmap_table,
                                 const CondensedExclusions &cond_excl,
                                 const BasicValenceTable &basic_vtable,
                                 const CharmmValenceTable &charmm_vtable,
                                 const AttenuationParameterSet &attn_parm,
                                 const VirtualSiteTable &vsite_table,
                                 const Map1234 &all_nb_excl, const ConstraintTable &cnst_table) {
  
  // Allocate the Hybrid ARRAY-kind objects based on the compiled data
  size_t int_items = roundUp<int>(tmp_desc.size(), warp_size_int) +
                     roundUp(residue_count + 1, warp_size_int) +
                     roundUp(molecule_count + 1, warp_size_int) +
                     roundUp(tmp_mobile_atoms.size(), warp_size_zu) +
                     13 * roundUp(atom_count, warp_size_int) +
                     10 * roundUp(atom_count + 1, warp_size_int) +
                      7 * roundUp(urey_bradley_term_count, warp_size_int) +
                     10 * roundUp(charmm_impr_term_count, warp_size_int) +
                     12 * roundUp(cmap_term_count, warp_size_int) +
                     roundUp(cmap_surface_count, warp_size_int) +
                      2 * roundUp(cmap_surface_count + 1, warp_size_int) +
                      6 * roundUp(bond_term_count, warp_size_int) +
                      8 * roundUp(angl_term_count, warp_size_int) +
                     11 * roundUp(dihe_term_count, warp_size_int) +
                      7 * roundUp(virtual_site_count, warp_size_int) +
                     roundUp(total_exclusions, warp_size_int) +
                      3 * roundUp(inferred_14_attenuations, warp_size_int) +
                     roundUp<int>(all_nb_excl.nb11_excl_list.size(), warp_size_int) +
                     roundUp<int>(all_nb_excl.nb12_excl_list.size(), warp_size_int) +
                     roundUp<int>(all_nb_excl.nb13_excl_list.size(), warp_size_int) +
                     roundUp<int>(all_nb_excl.nb14_excl_list.size(), warp_size_int) +
                      4 * roundUp(cnst_table.settle_group_count, warp_size_int) +
                     roundUp(cnst_table.cnst_group_list.size(), warp_size_zu) +
                     roundUp(cnst_table.cnst_group_count + 1, warp_size_int) +
                     roundUp(cnst_table.cnst_group_count, warp_size_int) +
                     roundUp(cnst_table.cnst_parameter_count + 1, warp_size_int);
  size_t double_items =  9 * roundUp(atom_count, warp_size_int) +
                         2 * roundUp(urey_bradley_parameter_count, warp_size_int) +
                         2 * roundUp(charmm_impr_parameter_count, warp_size_int) +
                        20 * roundUp<int>(tmp_cmap_surfaces.size(), warp_size_int) +
                         2 * roundUp(bond_parameter_count, warp_size_int) +
                         2 * roundUp(angl_parameter_count, warp_size_int) +
                         3 * roundUp(dihe_parameter_count, warp_size_int) +
                         3 * roundUp(virtual_site_count, warp_size_int) +
                        10 * roundUp(lj_type_count * lj_type_count, warp_size_int) +
                         2 * roundUp(lj_type_count, warp_size_int) +
                         roundUp(charge_type_count,  warp_size_int) +
                         2 * roundUp(attenuated_14_type_count, warp_size_int) +
                         6 * roundUp(cnst_table.settle_parameter_count, warp_size_int) +
                         2 * roundUp(cnst_table.cnst_parameter_list.size(), warp_size_zu);
  size_t float_items = double_items;
  size_t char4_items = 4 * roundUp(atom_count, warp_size_int) +
                       roundUp(residue_count, warp_size_int) +
                       2 * roundUp(bond_term_count, warp_size_int) +
                       2 * roundUp(angl_term_count, warp_size_int) +
                       2 * roundUp(dihe_term_count, warp_size_int);
  int_data.resize(int_items);
  double_data.resize(double_items);
  float_data.resize(float_items);
  char4_data.resize(char4_items);
  
  // Lay out Hybrid POINTER-kind int objects based on the compiled data.  The putHost member
  // function of the Hybrid class has an overloaded version that sets a POINTER-kind object to its
  // target and then fills the object, thus populating the appropriate segment of the target.
  size_t ic = descriptors.putHost(&int_data, tmp_desc, 0, warp_size_zu);
  ic = residue_limits.putHost(&int_data, tmp_residue_limits, ic, warp_size_zu);
  ic = atom_struc_numbers.putHost(&int_data, tmp_atom_struc_numbers, ic, warp_size_zu);
  ic = residue_numbers.putHost(&int_data, tmp_residue_numbers, ic, warp_size_zu);
  ic = molecule_limits.putHost(&int_data, tmp_molecule_limits, ic, warp_size_zu);
  ic = atomic_numbers.putHost(&int_data, tmp_atomic_numbers, ic, warp_size_zu);
  ic = molecule_membership.putHost(&int_data, tmp_molecule_membership, ic, warp_size_zu);
  ic = mobile_atoms.putHost(&int_data, tmp_mobile_atoms, ic, warp_size_zu);
  ic = molecule_contents.putHost(&int_data, tmp_molecule_contents, ic, warp_size_zu);
  ic = urey_bradley_i_atoms.putHost(&int_data, charmm_vtable.ubrd_i_atoms, ic, warp_size_zu);
  ic = urey_bradley_k_atoms.putHost(&int_data, charmm_vtable.ubrd_k_atoms, ic, warp_size_zu);
  ic = urey_bradley_parameter_indices.putHost(&int_data, charmm_vtable.ubrd_param_idx, ic,
                                              warp_size_zu);
  ic = urey_bradley_assigned_atoms.putHost(&int_data, charmm_vtable.ubrd_assigned_atoms, ic,
                                           warp_size_zu);
  ic = urey_bradley_assigned_index.putHost(&int_data, charmm_vtable.ubrd_assigned_index, ic,
                                           warp_size_zu);
  ic = urey_bradley_assigned_terms.putHost(&int_data, charmm_vtable.ubrd_assigned_terms, ic,
                                           warp_size_zu);
  ic = urey_bradley_assigned_bounds.putHost(&int_data, charmm_vtable.ubrd_assigned_bounds, ic,
                                            warp_size_zu);
  ic = charmm_impr_i_atoms.putHost(&int_data, charmm_vtable.impr_i_atoms, ic, warp_size_zu);
  ic = charmm_impr_j_atoms.putHost(&int_data, charmm_vtable.impr_j_atoms, ic, warp_size_zu);
  ic = charmm_impr_k_atoms.putHost(&int_data, charmm_vtable.impr_k_atoms, ic, warp_size_zu);
  ic = charmm_impr_l_atoms.putHost(&int_data, charmm_vtable.impr_l_atoms, ic, warp_size_zu);
  ic = charmm_impr_parameter_indices.putHost(&int_data, charmm_vtable.impr_param_idx, ic,
                                             warp_size_zu);
  ic = charmm_impr_assigned_atoms.putHost(&int_data, charmm_vtable.impr_assigned_atoms, ic,
                                          warp_size_zu);
  ic = charmm_impr_assigned_index.putHost(&int_data, charmm_vtable.impr_assigned_index, ic,
                                          warp_size_zu);
  ic = charmm_impr_assigned_terms.putHost(&int_data, charmm_vtable.impr_assigned_terms, ic,
                                          warp_size_zu);
  ic = charmm_impr_assigned_bounds.putHost(&int_data, charmm_vtable.impr_assigned_bounds, ic,
                                           warp_size_zu);
  ic = cmap_i_atoms.putHost(&int_data, charmm_vtable.cmap_i_atoms, ic, warp_size_zu);
  ic = cmap_j_atoms.putHost(&int_data, charmm_vtable.cmap_j_atoms, ic, warp_size_zu);
  ic = cmap_k_atoms.putHost(&int_data, charmm_vtable.cmap_k_atoms, ic, warp_size_zu);
  ic = cmap_l_atoms.putHost(&int_data, charmm_vtable.cmap_l_atoms, ic, warp_size_zu);
  ic = cmap_m_atoms.putHost(&int_data, charmm_vtable.cmap_m_atoms, ic, warp_size_zu);
  ic = cmap_surface_dimensions.putHost(&int_data, tmp_cmap_surf_dims, ic, warp_size_zu);
  ic = cmap_surface_bounds.putHost(&int_data, tmp_cmap_surf_bounds, ic, warp_size_zu);
  ic = cmap_patch_bounds.putHost(&int_data, cmap_table.patch_matrix_bounds, ic, warp_size_zu);
  ic = cmap_surface_indices.putHost(&int_data, charmm_vtable.cmap_param_idx, ic, warp_size_zu);
  ic = cmap_assigned_atoms.putHost(&int_data, charmm_vtable.cmap_assigned_atoms, ic, warp_size_zu);
  ic = cmap_assigned_index.putHost(&int_data, charmm_vtable.cmap_assigned_index, ic, warp_size_zu);
  ic = cmap_assigned_terms.putHost(&int_data, charmm_vtable.cmap_assigned_terms, ic, warp_size_zu);
  ic = cmap_assigned_bounds.putHost(&int_data, charmm_vtable.cmap_assigned_bounds, ic,
                                    warp_size_zu);
  ic = bond_i_atoms.putHost(&int_data, basic_vtable.bond_i_atoms, ic, warp_size_zu);
  ic = bond_j_atoms.putHost(&int_data, basic_vtable.bond_j_atoms, ic, warp_size_zu);
  ic = bond_parameter_indices.putHost(&int_data, basic_vtable.bond_param_idx, ic, warp_size_zu);
  ic = bond_assigned_atoms.putHost(&int_data, basic_vtable.bond_assigned_atoms, ic, warp_size_zu);
  ic = bond_assigned_index.putHost(&int_data, basic_vtable.bond_assigned_index, ic, warp_size_zu);
  ic = bond_assigned_terms.putHost(&int_data, basic_vtable.bond_assigned_terms, ic, warp_size_zu);
  ic = bond_assigned_bounds.putHost(&int_data, basic_vtable.bond_assigned_bounds, ic,
                                    warp_size_zu);
  ic = angl_i_atoms.putHost(&int_data, basic_vtable.angl_i_atoms, ic, warp_size_zu);
  ic = angl_j_atoms.putHost(&int_data, basic_vtable.angl_j_atoms, ic, warp_size_zu);
  ic = angl_k_atoms.putHost(&int_data, basic_vtable.angl_k_atoms, ic, warp_size_zu);
  ic = angl_parameter_indices.putHost(&int_data, basic_vtable.angl_param_idx, ic, warp_size_zu);
  ic = angl_assigned_atoms.putHost(&int_data, basic_vtable.angl_assigned_atoms, ic, warp_size_zu);
  ic = angl_assigned_index.putHost(&int_data, basic_vtable.angl_assigned_index, ic, warp_size_zu);
  ic = angl_assigned_terms.putHost(&int_data, basic_vtable.angl_assigned_terms, ic, warp_size_zu);
  ic = angl_assigned_bounds.putHost(&int_data, basic_vtable.angl_assigned_bounds, ic,
                                    warp_size_zu);
  ic = dihe_i_atoms.putHost(&int_data, basic_vtable.dihe_i_atoms, ic, warp_size_zu);
  ic = dihe_j_atoms.putHost(&int_data, basic_vtable.dihe_j_atoms, ic, warp_size_zu);
  ic = dihe_k_atoms.putHost(&int_data, basic_vtable.dihe_k_atoms, ic, warp_size_zu);
  ic = dihe_l_atoms.putHost(&int_data, basic_vtable.dihe_l_atoms, ic, warp_size_zu);
  ic = dihe_parameter_indices.putHost(&int_data, basic_vtable.dihe_param_idx, ic, warp_size_zu);
  ic = dihe14_parameter_indices.putHost(&int_data, attn_parm.dihe14_parameter_indices, ic,
                                        warp_size_zu);
  ic = dihe_assigned_atoms.putHost(&int_data, basic_vtable.dihe_assigned_atoms, ic, warp_size_zu);
  ic = dihe_assigned_index.putHost(&int_data, basic_vtable.dihe_assigned_index, ic, warp_size_zu);
  ic = dihe_assigned_terms.putHost(&int_data, basic_vtable.dihe_assigned_terms, ic, warp_size_zu);
  ic = dihe_assigned_bounds.putHost(&int_data, basic_vtable.dihe_assigned_bounds, ic,
                                    warp_size_zu);
  ic = virtual_site_atoms.putHost(&int_data, vsite_table.vs_atoms, ic, warp_size_zu);
  ic = virtual_site_frame_types.putHost(&int_data, vsite_table.frame_types, ic, warp_size_zu);
  ic = virtual_site_frame1_atoms.putHost(&int_data, vsite_table.frame1_atoms, ic, warp_size_zu);
  ic = virtual_site_frame2_atoms.putHost(&int_data, vsite_table.frame2_atoms, ic, warp_size_zu);
  ic = virtual_site_frame3_atoms.putHost(&int_data, vsite_table.frame3_atoms, ic, warp_size_zu);
  ic = virtual_site_frame4_atoms.putHost(&int_data, vsite_table.frame4_atoms, ic, warp_size_zu);
  ic = virtual_site_parameter_indices.putHost(&int_data, vsite_table.param_idx, ic, warp_size_zu);
  ic = charge_indices.putHost(&int_data, tmp_charge_indices, ic, warp_size_zu);
  ic = lennard_jones_indices.putHost(&int_data, tmp_lennard_jones_indices, ic, warp_size_zu);
  ic = atom_exclusion_bounds.putHost(&int_data, cond_excl.atom_excl_bounds, ic, warp_size_zu);
  ic = atom_exclusion_list.putHost(&int_data, cond_excl.atom_excl_list, ic, warp_size_zu);
  ic = nb11_exclusion_bounds.putHost(&int_data, all_nb_excl.nb11_excl_bounds, ic, warp_size_zu);
  ic = nb11_exclusion_list.putHost(&int_data, all_nb_excl.nb11_excl_list, ic, warp_size_zu);
  ic = nb12_exclusion_bounds.putHost(&int_data, all_nb_excl.nb12_excl_bounds, ic, warp_size_zu);
  ic = nb12_exclusion_list.putHost(&int_data, all_nb_excl.nb12_excl_list, ic, warp_size_zu);
  ic = nb13_exclusion_bounds.putHost(&int_data, all_nb_excl.nb13_excl_bounds, ic, warp_size_zu);
  ic = nb13_exclusion_list.putHost(&int_data, all_nb_excl.nb13_excl_list, ic, warp_size_zu);
  ic = nb14_exclusion_bounds.putHost(&int_data, all_nb_excl.nb14_excl_bounds, ic, warp_size_zu);
  ic = nb14_exclusion_list.putHost(&int_data, all_nb_excl.nb14_excl_list, ic, warp_size_zu);
  ic = infr14_i_atoms.putHost(&int_data, tmp_inferred_14_i_atoms, ic, warp_size_zu);
  ic = infr14_l_atoms.putHost(&int_data, tmp_inferred_14_l_atoms, ic, warp_size_zu);
  ic = infr14_parameter_indices.putHost(&int_data, tmp_inferred_14_param_idx, ic, warp_size_zu);
  ic = neck_gb_indices.putHost(&int_data, tmp_neck_gb_indices, ic, warp_size_zu);
  ic = settle_oxygen_atoms.putHost(&int_data, cnst_table.settle_ox_atoms, ic, warp_size_zu);
  ic = settle_hydro1_atoms.putHost(&int_data, cnst_table.settle_h1_atoms, ic, warp_size_zu);
  ic = settle_hydro2_atoms.putHost(&int_data, cnst_table.settle_h2_atoms, ic, warp_size_zu);
  ic = settle_parameter_indices.putHost(&int_data, cnst_table.settle_param_idx, ic, warp_size_zu);
  ic = constraint_group_atoms.putHost(&int_data, cnst_table.cnst_group_list, ic, warp_size_zu);
  ic = constraint_group_bounds.putHost(&int_data, cnst_table.cnst_group_bounds, ic, warp_size_zu);
  ic = constraint_parameter_bounds.putHost(&int_data, cnst_table.cnst_parameter_bounds, ic,
                                           warp_size_zu);
  ic = constraint_parameter_indices.putHost(&int_data, cnst_table.cnst_group_param_idx, ic,
                                            warp_size_zu);
  ic = tree_joining_info.putHost(&int_data, tmp_tree_joining_info, ic, warp_size_zu);
  ic = last_rotator_info.putHost(&int_data, tmp_last_rotator_info, ic, warp_size_zu);

  // Do the same for double Hybrid POINTER-kind objects
  size_t dc = atomic_charges.putHost(&double_data, tmp_charges, 0, warp_size_zu);
  dc = atomic_masses.putHost(&double_data, tmp_masses, dc, warp_size_zu);
  std::vector<double> inv_mass(atom_count, 0.0);
  for (int i = 0; i < atom_count; i++) {
    inv_mass[i] = (tmp_masses[i] > constants::tiny) ? 1.0 / tmp_masses[i] : 0.0;
  }
  const size_t nlj_pairs = tmp_lj_a_values.size();
  std::vector<double> tmp_lj_sigma_values(nlj_pairs), tmp_lj_14_sigma_values(nlj_pairs);
  for (size_t i = 0; i < nlj_pairs; i++) {
    if (tmp_lj_b_values[i] > 1.0e-6) {
      tmp_lj_sigma_values[i] = sqrt(cbrt(tmp_lj_a_values[i] / tmp_lj_b_values[i]));
    }
    else {
      tmp_lj_sigma_values[i] = 0.0;
    }
    if (tmp_lj_14_b_values[i] > 1.0e-6) {
      tmp_lj_14_sigma_values[i] = sqrt(cbrt(tmp_lj_14_a_values[i] / tmp_lj_14_b_values[i]));
    }
    else {
      tmp_lj_14_sigma_values[i] = 0.0;
    }
  }
  dc = inverse_atomic_masses.putHost(&double_data, inv_mass, dc, warp_size_zu);
  dc = urey_bradley_stiffnesses.putHost(&double_data, tmp_ub_stiffnesses, dc, warp_size_zu);
  dc = urey_bradley_equilibria.putHost(&double_data, tmp_ub_equilibria, dc, warp_size_zu);
  dc = charmm_impr_stiffnesses.putHost(&double_data, tmp_charmm_impr_stiffnesses, dc,
                                       warp_size_zu);
  dc = charmm_impr_phase_angles.putHost(&double_data, tmp_charmm_impr_phase_angles, dc,
                                        warp_size_zu);
  dc = cmap_surfaces.putHost(&double_data, tmp_cmap_surfaces, dc, warp_size_zu);
  dc = cmap_phi_derivatives.putHost(&double_data, cmap_table.phi_derivatives, dc, warp_size_zu);
  dc = cmap_psi_derivatives.putHost(&double_data, cmap_table.psi_derivatives, dc, warp_size_zu);
  dc = cmap_phi_psi_derivatives.putHost(&double_data, cmap_table.phi_psi_derivatives, dc,
                                        warp_size_zu);
  dc = cmap_patches.putHost(&double_data, cmap_table.patch_matrix_form, dc, warp_size_zu);
  dc = bond_stiffnesses.putHost(&double_data, tmp_bond_stiffnesses, dc, warp_size_zu);
  dc = bond_equilibria.putHost(&double_data, tmp_bond_equilibria, dc, warp_size_zu);
  dc = angl_stiffnesses.putHost(&double_data, tmp_angl_stiffnesses, dc, warp_size_zu);
  dc = angl_equilibria.putHost(&double_data, tmp_angl_equilibria, dc, warp_size_zu);
  dc = dihe_amplitudes.putHost(&double_data, tmp_dihe_amplitudes, dc, warp_size_zu);
  dc = dihe_periodicities.putHost(&double_data, tmp_dihe_periodicities, dc, warp_size_zu);
  dc = dihe_phase_angles.putHost(&double_data, tmp_dihe_phase_angles, dc, warp_size_zu);
  dc = virtual_site_frame_dim1.putHost(&double_data, vsite_table.frame_dim1, dc, warp_size_zu);
  dc = virtual_site_frame_dim2.putHost(&double_data, vsite_table.frame_dim2, dc, warp_size_zu);
  dc = virtual_site_frame_dim3.putHost(&double_data, vsite_table.frame_dim3, dc, warp_size_zu);
  dc = charge_parameters.putHost(&double_data, tmp_charge_parameters, dc, warp_size_zu);
  dc = lj_a_values.putHost(&double_data, tmp_lj_a_values, dc, warp_size_zu);
  dc = lj_b_values.putHost(&double_data, tmp_lj_b_values, dc, warp_size_zu);
  dc = lj_c_values.putHost(&double_data, tmp_lj_c_values, dc, warp_size_zu);
  dc = lj_14_a_values.putHost(&double_data, tmp_lj_14_a_values, dc, warp_size_zu);
  dc = lj_14_b_values.putHost(&double_data, tmp_lj_14_b_values, dc, warp_size_zu);
  dc = lj_14_c_values.putHost(&double_data, tmp_lj_14_c_values, dc, warp_size_zu);
  dc = lj_sigma_values.putHost(&double_data, tmp_lj_sigma_values, dc, warp_size_zu);
  dc = lj_14_sigma_values.putHost(&double_data, tmp_lj_14_sigma_values, dc, warp_size_zu);
  dc = lj_type_corrections.putHost(&double_data, std::vector<double>(lj_type_count, 0.0), dc,
                                   warp_size_zu);
  dc = attn14_elec_factors.putHost(&double_data, attn_parm.elec_screening_factors, dc,
                                   warp_size_zu);
  dc = attn14_vdw_factors.putHost(&double_data, attn_parm.vdw_screening_factors, dc, warp_size_zu);
  dc = atomic_pb_radii.putHost(&double_data, tmp_atomic_pb_radii, dc, warp_size_zu);
  dc = gb_screening_factors.putHost(&double_data, tmp_gb_screening_factors, dc, warp_size_zu);
  dc = gb_alpha_parameters.putHost(&double_data, tmp_gb_coef, dc, warp_size_zu);
  dc = gb_beta_parameters.putHost(&double_data, tmp_gb_coef, dc, warp_size_zu);
  dc = gb_gamma_parameters.putHost(&double_data, tmp_gb_coef, dc, warp_size_zu);

  // Parse the SETTLE and constraint group parameters into arrays of simpler reals.
  std::vector<double> tmp_settle_mormt(cnst_table.settle_parameter_count);
  std::vector<double> tmp_settle_mhrmt(cnst_table.settle_parameter_count);
  std::vector<double> tmp_settle_ra(cnst_table.settle_parameter_count);
  std::vector<double> tmp_settle_rb(cnst_table.settle_parameter_count);
  std::vector<double> tmp_settle_rc(cnst_table.settle_parameter_count);
  std::vector<double> tmp_settle_invra(cnst_table.settle_parameter_count);
  std::vector<float> sp_tmp_settle_mormt(cnst_table.settle_parameter_count);
  std::vector<float> sp_tmp_settle_mhrmt(cnst_table.settle_parameter_count);
  std::vector<float> sp_tmp_settle_ra(cnst_table.settle_parameter_count);
  std::vector<float> sp_tmp_settle_rb(cnst_table.settle_parameter_count);
  std::vector<float> sp_tmp_settle_rc(cnst_table.settle_parameter_count);
  std::vector<float> sp_tmp_settle_invra(cnst_table.settle_parameter_count);
  for (int i = 0; i < cnst_table.settle_parameter_count; i++) {
    tmp_settle_mormt[i]    = cnst_table.settle_measurements[i].mormt;
    tmp_settle_mhrmt[i]    = cnst_table.settle_measurements[i].mhrmt;
    tmp_settle_ra[i]       = cnst_table.settle_measurements[i].ra;
    tmp_settle_rb[i]       = cnst_table.settle_measurements[i].rb;
    tmp_settle_rc[i]       = cnst_table.settle_measurements[i].rc;
    tmp_settle_invra[i]    = cnst_table.settle_measurements[i].invra;
    sp_tmp_settle_mormt[i] = cnst_table.settle_measurements[i].mormt;
    sp_tmp_settle_mhrmt[i] = cnst_table.settle_measurements[i].mhrmt;
    sp_tmp_settle_ra[i]    = cnst_table.settle_measurements[i].ra;
    sp_tmp_settle_rb[i]    = cnst_table.settle_measurements[i].rb;
    sp_tmp_settle_rc[i]    = cnst_table.settle_measurements[i].rc;
    sp_tmp_settle_invra[i] = cnst_table.settle_measurements[i].invra;
  }
  const int ncnst_val = cnst_table.cnst_parameter_list.size();
  std::vector<double> tmp_cnst_inv_masses(ncnst_val);
  std::vector<double> tmp_cnst_squared_lengths(ncnst_val);
  std::vector<float> sp_tmp_cnst_inv_masses(ncnst_val);
  std::vector<float> sp_tmp_cnst_squared_lengths(ncnst_val);
  for (int i = 0; i < ncnst_val; i++) {
    tmp_cnst_inv_masses[i]        = cnst_table.cnst_parameter_list[i].x;
    tmp_cnst_squared_lengths[i]    = cnst_table.cnst_parameter_list[i].y;
    sp_tmp_cnst_inv_masses[i]     = cnst_table.cnst_parameter_list[i].x;
    sp_tmp_cnst_squared_lengths[i] = cnst_table.cnst_parameter_list[i].y;
  }
  dc = settle_mormt.putHost(&double_data, tmp_settle_mormt, dc, warp_size_zu);
  dc = settle_mhrmt.putHost(&double_data, tmp_settle_mhrmt, dc, warp_size_zu);
  dc = settle_ra.putHost(&double_data, tmp_settle_ra, dc, warp_size_zu);
  dc = settle_rb.putHost(&double_data, tmp_settle_rb, dc, warp_size_zu);
  dc = settle_rc.putHost(&double_data, tmp_settle_rc, dc, warp_size_zu);
  dc = settle_invra.putHost(&double_data, tmp_settle_invra, dc, warp_size_zu);  
  dc = constraint_inverse_masses.putHost(&double_data, tmp_cnst_inv_masses, dc, warp_size_zu);
  dc = constraint_squared_lengths.putHost(&double_data, tmp_cnst_squared_lengths, dc,
                                          warp_size_zu);
  dc = solty_info.putHost(&double_data, tmp_solty_info, dc, warp_size_zu);
  dc = hbond_a_values.putHost(&double_data, tmp_hbond_a_values, dc, warp_size_zu);
  dc = hbond_b_values.putHost(&double_data, tmp_hbond_b_values, dc, warp_size_zu);
  dc = hbond_cutoffs.putHost(&double_data, tmp_hbond_cutoffs, dc, warp_size_zu);

  // Do the same for float Hybrid POINTER-kind objects, using the range constructor to make
  // single-precision vectors out of their dobule-precision counterparts before loading the
  // Hybrid objects.
  const std::vector<float>sp_tmp_charges(tmp_charges.begin(), tmp_charges.end());
  size_t fc = sp_atomic_charges.putHost(&float_data, sp_tmp_charges, 0, warp_size_zu);
  const std::vector<float>sp_tmp_masses(tmp_masses.begin(), tmp_masses.end());
  fc = sp_atomic_masses.putHost(&float_data, sp_tmp_masses, fc, warp_size_zu);
  std::vector<float> sp_inv_mass(atom_count, 0.0);
  for (int i = 0; i < atom_count; i++) {
    sp_inv_mass[i] = (tmp_masses[i] > constants::tiny) ? 1.0 / tmp_masses[i] : 0.0;
  }
  fc = sp_inverse_atomic_masses.putHost(&float_data, sp_inv_mass, fc, warp_size_zu);
  const std::vector<float> sp_tmp_ub_stiffnesses(tmp_ub_stiffnesses.begin(),
                                                 tmp_ub_stiffnesses.end());
  fc = sp_urey_bradley_stiffnesses.putHost(&float_data, sp_tmp_ub_stiffnesses, fc, warp_size_zu);
  const std::vector<float> sp_tmp_ub_equilibria(tmp_ub_equilibria.begin(),
                                                tmp_ub_equilibria.end());
  fc = sp_urey_bradley_equilibria.putHost(&float_data, sp_tmp_ub_equilibria, fc, warp_size_zu);
  const std::vector<float> sp_tmp_charmm_impr_stiffnesses(tmp_charmm_impr_stiffnesses.begin(),
                                                          tmp_charmm_impr_stiffnesses.end());
  fc = sp_charmm_impr_stiffnesses.putHost(&float_data, sp_tmp_charmm_impr_stiffnesses, fc,
                                          warp_size_zu);
  const std::vector<float> sp_tmp_charmm_impr_phase_angles(tmp_charmm_impr_phase_angles.begin(),
                                                           tmp_charmm_impr_phase_angles.end());
  fc = sp_charmm_impr_phase_angles.putHost(&float_data, sp_tmp_charmm_impr_phase_angles, fc,
                                           warp_size_zu);
  const std::vector<float> sp_tmp_cmap_surfaces(tmp_cmap_surfaces.begin(),
                                                tmp_cmap_surfaces.end());
  fc = sp_cmap_surfaces.putHost(&float_data, sp_tmp_cmap_surfaces, fc, warp_size_zu);
  const std::vector<float> sp_tmp_cmap_dphi(cmap_table.phi_derivatives.begin(),
                                            cmap_table.phi_derivatives.end());
  fc = sp_cmap_phi_derivatives.putHost(&float_data, sp_tmp_cmap_dphi, fc, warp_size_zu);
  const std::vector<float> sp_tmp_cmap_dpsi(cmap_table.psi_derivatives.begin(),
                                            cmap_table.psi_derivatives.end());
  fc = sp_cmap_psi_derivatives.putHost(&float_data, sp_tmp_cmap_dpsi, fc, warp_size_zu);
  const std::vector<float> sp_tmp_cmap_dphipsi(cmap_table.phi_psi_derivatives.begin(),
                                               cmap_table.phi_psi_derivatives.end());
  fc = sp_cmap_phi_psi_derivatives.putHost(&float_data, sp_tmp_cmap_dphipsi, fc, warp_size_zu);
  const std::vector<float> sp_tmp_cmap_patches(cmap_table.patch_matrix_form.begin(),
                                               cmap_table.patch_matrix_form.end());
  fc = sp_cmap_patches.putHost(&float_data, sp_tmp_cmap_patches, fc, warp_size_zu);
  const std::vector<float> sp_tmp_bond_stiffnesses(tmp_bond_stiffnesses.begin(),
                                                   tmp_bond_stiffnesses.end());
  fc = sp_bond_stiffnesses.putHost(&float_data, sp_tmp_bond_stiffnesses, fc, warp_size_zu);
  const std::vector<float> sp_tmp_bond_equilibria(tmp_bond_equilibria.begin(),
                                                  tmp_bond_equilibria.end());
  fc = sp_bond_equilibria.putHost(&float_data, sp_tmp_bond_equilibria, fc, warp_size_zu);
  const std::vector<float> sp_tmp_angl_stiffnesses(tmp_angl_stiffnesses.begin(),
                                                   tmp_angl_stiffnesses.end());
  fc = sp_angl_stiffnesses.putHost(&float_data, sp_tmp_angl_stiffnesses, fc, warp_size_zu);
  const std::vector<float> sp_tmp_angl_equilibria(tmp_angl_equilibria.begin(),
                                                  tmp_angl_equilibria.end());
  fc = sp_angl_equilibria.putHost(&float_data, sp_tmp_angl_equilibria, fc, warp_size_zu);
  const std::vector<float> sp_tmp_dihe_amplitudes(tmp_dihe_amplitudes.begin(),
                                                  tmp_dihe_amplitudes.end());
  fc = sp_dihe_amplitudes.putHost(&float_data, sp_tmp_dihe_amplitudes, fc, warp_size_zu);
  const std::vector<float> sp_tmp_dihe_periodicities(tmp_dihe_periodicities.begin(),
                                                     tmp_dihe_periodicities.end());
  fc = sp_dihe_periodicities.putHost(&float_data, sp_tmp_dihe_periodicities, fc, warp_size_zu);
  const std::vector<float> sp_tmp_dihe_phase_angles(tmp_dihe_phase_angles.begin(),
                                                    tmp_dihe_phase_angles.end());
  fc = sp_dihe_phase_angles.putHost(&float_data, sp_tmp_dihe_phase_angles, fc, warp_size_zu);
  const std::vector<float> sp_frame_dim1(vsite_table.frame_dim1.begin(),
                                         vsite_table.frame_dim1.end());
  fc = sp_virtual_site_frame_dim1.putHost(&float_data, sp_frame_dim1, fc, warp_size_zu);
  const std::vector<float> sp_frame_dim2(vsite_table.frame_dim2.begin(),
                                         vsite_table.frame_dim2.end());
  fc = sp_virtual_site_frame_dim2.putHost(&float_data, sp_frame_dim2, fc, warp_size_zu);
  const std::vector<float> sp_frame_dim3(vsite_table.frame_dim3.begin(),
                                         vsite_table.frame_dim3.end());
  fc = sp_virtual_site_frame_dim3.putHost(&float_data, sp_frame_dim3, fc, warp_size_zu);
  const std::vector<float> sp_tmp_charge_parameters(tmp_charge_parameters.begin(),
                                                    tmp_charge_parameters.end());
  fc = sp_charge_parameters.putHost(&float_data, sp_tmp_charge_parameters, fc, warp_size_zu);
  const std::vector<float> sp_tmp_lj_a_values(tmp_lj_a_values.begin(), tmp_lj_a_values.end());
  fc = sp_lj_a_values.putHost(&float_data, sp_tmp_lj_a_values, fc, warp_size_zu);
  const std::vector<float> sp_tmp_lj_b_values(tmp_lj_b_values.begin(), tmp_lj_b_values.end());
  fc = sp_lj_b_values.putHost(&float_data, sp_tmp_lj_b_values, fc, warp_size_zu);
  const std::vector<float> sp_tmp_lj_c_values(tmp_lj_c_values.begin(), tmp_lj_c_values.end());
  fc = sp_lj_c_values.putHost(&float_data, sp_tmp_lj_c_values, fc, warp_size_zu);
  const std::vector<float> sp_tmp_lj_14_a_values(tmp_lj_14_a_values.begin(),
                                                 tmp_lj_14_a_values.end());
  fc = sp_lj_14_a_values.putHost(&float_data, sp_tmp_lj_14_a_values, fc, warp_size_zu);
  const std::vector<float> sp_tmp_lj_14_b_values(tmp_lj_14_b_values.begin(),
                                                 tmp_lj_14_b_values.end());
  fc = sp_lj_14_b_values.putHost(&float_data, sp_tmp_lj_14_b_values, fc, warp_size_zu);
  const std::vector<float> sp_tmp_lj_14_c_values(tmp_lj_14_c_values.begin(),
                                                 tmp_lj_14_c_values.end());
  fc = sp_lj_14_c_values.putHost(&float_data, sp_tmp_lj_14_c_values, fc, warp_size_zu);
  const std::vector<float>sp_tmp_lj_sigma_values(tmp_lj_sigma_values.begin(),
                                                 tmp_lj_sigma_values.end());
  fc = sp_lj_sigma_values.putHost(&float_data, sp_tmp_lj_sigma_values, fc, warp_size_zu);
  const std::vector<float>sp_tmp_lj_14_sigma_values(tmp_lj_14_sigma_values.begin(),
                                                    tmp_lj_14_sigma_values.end());
  fc = sp_lj_14_sigma_values.putHost(&float_data, sp_tmp_lj_14_sigma_values, fc, warp_size_zu);
  fc = sp_lj_type_corrections.putHost(&float_data, std::vector<float>(lj_type_count, 0.0), fc,
                                      warp_size_zu);
  const std::vector<float> sp_tmp_elec_screening_factors(attn_parm.elec_screening_factors.begin(),
                                                         attn_parm.elec_screening_factors.end());
  fc = sp_attn14_elec_factors.putHost(&float_data, sp_tmp_elec_screening_factors, fc,
                                      warp_size_zu);
  const std::vector<float> sp_tmp_vdw_screening_factors(attn_parm.vdw_screening_factors.begin(),
                                                        attn_parm.vdw_screening_factors.end());
  fc = sp_attn14_vdw_factors.putHost(&float_data, sp_tmp_vdw_screening_factors, fc, warp_size_zu);
  const std::vector<float> sp_tmp_atomic_pb_radii(tmp_atomic_pb_radii.begin(),
                                                  tmp_atomic_pb_radii.end());
  fc = sp_atomic_pb_radii.putHost(&float_data, sp_tmp_atomic_pb_radii, fc, warp_size_zu);
  const std::vector<float> sp_tmp_gb_screening_factors(tmp_gb_screening_factors.begin(),
                                                       tmp_gb_screening_factors.end());
  fc = sp_gb_screening_factors.putHost(&float_data, sp_tmp_gb_screening_factors, fc, warp_size_zu);
  const std::vector<float> sp_tmp_gb_coef(tmp_gb_coef.begin(), tmp_gb_coef.end());
  fc = sp_gb_alpha_parameters.putHost(&float_data, sp_tmp_gb_coef, fc, warp_size_zu);
  fc = sp_gb_beta_parameters.putHost(&float_data, sp_tmp_gb_coef, fc, warp_size_zu);
  fc = sp_gb_gamma_parameters.putHost(&float_data, sp_tmp_gb_coef, fc, warp_size_zu);
  fc = sp_settle_mormt.putHost(&float_data, sp_tmp_settle_mormt, fc, warp_size_zu);
  fc = sp_settle_mhrmt.putHost(&float_data, sp_tmp_settle_mhrmt, fc, warp_size_zu);
  fc = sp_settle_ra.putHost(&float_data, sp_tmp_settle_ra, fc, warp_size_zu);
  fc = sp_settle_rb.putHost(&float_data, sp_tmp_settle_rb, fc, warp_size_zu);
  fc = sp_settle_rc.putHost(&float_data, sp_tmp_settle_rc, fc, warp_size_zu);
  fc = sp_settle_invra.putHost(&float_data, sp_tmp_settle_invra, fc, warp_size_zu);
  fc = sp_constraint_inverse_masses.putHost(&float_data, sp_tmp_cnst_inv_masses, fc, warp_size_zu);
  fc = sp_constraint_squared_lengths.putHost(&float_data, sp_tmp_cnst_squared_lengths, fc,
                                            warp_size_zu);
  
  // Do the same for char4 Hybrid POINTER-kind objects
  size_t c4c = atom_names.putHost(&char4_data, tmp_atom_names, 0, warp_size_zu);
  c4c = atom_types.putHost(&char4_data, tmp_atom_types, c4c, warp_size_zu);
  c4c = residue_names.putHost(&char4_data, tmp_residue_names, c4c, warp_size_zu);
  c4c = bond_modifiers.putHost(&char4_data, basic_vtable.bond_mods, c4c, warp_size_zu);
  c4c = angl_modifiers.putHost(&char4_data, basic_vtable.angl_mods, c4c, warp_size_zu);
  c4c = dihe_modifiers.putHost(&char4_data, basic_vtable.dihe_mods, c4c, warp_size_zu);
  c4c = bond_assigned_mods.putHost(&char4_data, basic_vtable.bond_assigned_mods, c4c,
                                   warp_size_zu);
  c4c = angl_assigned_mods.putHost(&char4_data, basic_vtable.angl_assigned_mods, c4c,
                                   warp_size_zu);
  c4c = dihe_assigned_mods.putHost(&char4_data, basic_vtable.dihe_assigned_mods, c4c,
                                   warp_size_zu);
  c4c = tree_symbols.putHost(&char4_data, tmp_tree_symbols, c4c, warp_size_zu);
  
  // The Amber topology read here does not contain overflow names of any sort
  const std::vector<char4> blank_overflow;
  c4c = atom_overflow_names.putHost(&char4_data, blank_overflow, c4c, warp_size_zu);
  c4c = atom_overflow_types.putHost(&char4_data, blank_overflow, c4c, warp_size_zu);
  c4c = residue_overflow_names.putHost(&char4_data, blank_overflow, c4c, warp_size_zu);
}
  
//-------------------------------------------------------------------------------------------------
void AtomGraph::buildFromPrmtop(const std::string &file_name, const ExceptionResponse policy,
                                const double coulomb_constant_in,
                                const double default_elec14_screening,
                                const double default_vdw14_screening,
                                const double charge_rounding_tol,
                                const double charge_discretization) {

  // Log the file that was the source
  source = file_name;

  // Set Coulomb's constant (to the Amber-specific value if the default setting for this function
  // remains in place)
  coulomb_constant = coulomb_constant_in;

  // Get the input file as a big text vector
  const TextFile fmem(file_name, TextOrigin::DISK, std::string(""),
                      "Prmtop-based AtomGraph constructor");

  // Begin parsing, starting with the version stamp and date
  const TextFileReader tfr = fmem.data();
  if (tfr.line_count > 0 && strncmp(tfr.text, "%VERSION", 8) == 0) {
    if (tfr.line_limits[1] >= 34) {
      for (int i = 26; i < 35; i++) {
        version_stamp[i - 26] = tfr.text[i];
      }
      version_stamp[9] = '\0';
    }
    else {
      snprintf(version_stamp, 16, "UNKNOWN");
    }
    if (tfr.line_limits[1] >= 61 && strncmp(&tfr.text[37], "DATE", 4) == 0) {
      char td[4];
      td[0] = tfr.text[44];
      td[1] = tfr.text[45];
      td[2] = '\0';
      date.tm_mon = verifyNumberFormat(td, NumberFormat::INTEGER, 0, 2) ? atoi(td) - 1 : 0;
      td[0] = tfr.text[47];
      td[1] = tfr.text[48];
      td[2] = '\0';
      date.tm_mday = verifyNumberFormat(td, NumberFormat::INTEGER, 0, 2) ? atoi(td) : 0;
      td[0] = tfr.text[50];
      td[1] = tfr.text[51];
      td[2] = '\0';
      date.tm_year = verifyNumberFormat(td, NumberFormat::INTEGER, 0, 2) ? atoi(td) + 100 : 0;
      td[0] = tfr.text[54];
      td[1] = tfr.text[55];
      td[2] = '\0';
      date.tm_hour = verifyNumberFormat(td, NumberFormat::INTEGER, 0, 2) ? atoi(td) : 0;
      td[0] = tfr.text[57];
      td[1] = tfr.text[58];
      td[2] = '\0';
      date.tm_min = verifyNumberFormat(td, NumberFormat::INTEGER, 0, 2) ? atoi(td) : 0;
      td[0] = tfr.text[60];
      td[1] = tfr.text[61];
      td[2] = '\0';
      date.tm_sec = verifyNumberFormat(td, NumberFormat::INTEGER, 0, 2) ? atoi(td) : 0;
    }
    else {
      time_t rawtime;
      tm *timeinfo;
      time(&rawtime);
      timeinfo = localtime(&rawtime);
      date = *timeinfo;
    }
  }

  // Vector to hold the detected format from section to section (this is VERY volatile)
  std::vector<int4> dfmt;

  // Read title
  int lstart = scanToFlag(fmem, "TITLE", &dfmt, TopologyRequirement::OPTIONAL);
  if (lstart == -1) {
    title = "";
  }
  else {
    const int title_char_count = tfr.line_limits[lstart + 1] - tfr.line_limits[lstart];
    title.resize(title_char_count);
    for (int i = 0; i < title_char_count; i++) {
      title[i] = tfr.text[tfr.line_limits[lstart] + i];
    }
  }

  // Get all of the descriptors.  The current Amber topology has one optional descriptor.
  const ulint max_descriptors = static_cast<ulint>(TopologyDescriptor::N_VALUES);
  lstart = scanToFlag(fmem, "POINTERS", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<int> tmp_desc = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                               max_descriptors - 1, max_descriptors);

  // Assign descriptors to counts of atoms, residues, and other parts of the system
  atom_count = tmp_desc[0];
  residue_count = tmp_desc[11];
  largest_residue_size = tmp_desc[28];
  implicit_copy_count = (tmp_desc.size() > 31) ? tmp_desc[31] : 0;

  // Assign descriptors relevant to the bonded calculation
  bond_term_with_hydrogen = tmp_desc[2];
  angl_term_with_hydrogen = tmp_desc[4];
  dihe_term_with_hydrogen = tmp_desc[6];
  bond_term_without_hydrogen = tmp_desc[3];
  angl_term_without_hydrogen = tmp_desc[5];
  dihe_term_without_hydrogen = tmp_desc[7];
  bond_term_count = bond_term_with_hydrogen + bond_term_without_hydrogen;
  angl_term_count = angl_term_with_hydrogen + angl_term_without_hydrogen;
  dihe_term_count = dihe_term_with_hydrogen + dihe_term_without_hydrogen;
  bond_parameter_count = tmp_desc[15];
  angl_parameter_count = tmp_desc[16];
  dihe_parameter_count = tmp_desc[17];
  bond_perturbation_term_count = tmp_desc[21];
  angl_perturbation_term_count = tmp_desc[22];
  dihe_perturbation_term_count = tmp_desc[23];
  bonds_in_perturbed_group = tmp_desc[24];
  angls_in_perturbed_group = tmp_desc[25];
  dihes_in_perturbed_group = tmp_desc[26];

  // Assign descriptors relevant to virtual site placement
  virtual_site_count = tmp_desc[30];

  // Assign descriptors relevant to the non-bonded calculation
  lj_type_count = tmp_desc[1];
  total_exclusions = tmp_desc[10];
  switch (tmp_desc[27]) {
  case 0:
    periodic_box_class = UnitCellType::NONE;
    break;
  case 1:
    periodic_box_class = UnitCellType::ORTHORHOMBIC;
    break;
  case 2:
    periodic_box_class = UnitCellType::TRICLINIC;
    break;
  default:
    rtErr("The IFBOX field has an invalid value of " + std::to_string(tmp_desc[27]), "AtomGraph");
  }

  // Assign descriptors relevant to the MD propagation algorithm
  switch (tmp_desc[20]) {
  case 0:
    use_perturbation_info = PerturbationSetting::OFF;
    break;
  case 1:
    use_perturbation_info = PerturbationSetting::ON;
    break;
  default:
    rtErr("The IFPERT field has an invalid value of " + std::to_string(tmp_desc[20]), "AtomGraph");
  }
  switch (tmp_desc[29]) {
  case 0:
    use_solvent_cap_option = SolventCapSetting::OFF;
    break;
  case 1:
    use_solvent_cap_option = SolventCapSetting::ON;
    break;
  default:
    rtErr("The IFCAP field has an invalid value of " + std::to_string(tmp_desc[29]), "AtomGraph");
  }

  // Assign descriptors pertaining to deprecated or unused information
  unused_nhparm = tmp_desc[8];
  unused_nparm = tmp_desc[9];
  unused_natyp = tmp_desc[18];
  hbond_10_12_parameter_count = tmp_desc[19];
  heavy_bonds_plus_constraints = tmp_desc[12];
  heavy_angls_plus_constraints = tmp_desc[13];
  heavy_dihes_plus_constraints = tmp_desc[14];

  // Read the force field citations
  lstart = scanToFlag(fmem, "FORCE_FIELD_TYPE", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  if (lstart >= 0) {
    force_fields = readForceFieldReferences(fmem, lstart);
  }

  // Read atom names
  lstart = scanToFlag(fmem, "ATOM_NAME", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<char4> tmp_atom_names = c4AmberPrmtopData(fmem, lstart, dfmt[0].x, atom_count);

  // Read atomic numbers
  lstart = scanToFlag(fmem, "ATOMIC_NUMBER", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<int> tmp_atomic_numbers;
  if (lstart >= 0) {
    tmp_atomic_numbers = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, atom_count);
  }

  // Read charges and convert to internal units (atomic units for charges)
  lstart = scanToFlag(fmem, "CHARGE", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<double> tmp_charges = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                     atom_count);
  for (int i = 0; i < atom_count; i++) {
    tmp_charges[i] *= inv_amber_charge_scaling;
  }
  std::vector<double> tmp_charge_parameters(atom_count);
  std::vector<int> tmp_charge_indices(atom_count);
  smoothCharges(&tmp_charges, &tmp_charge_parameters, &tmp_charge_indices,
                &charge_type_count, charge_rounding_tol, charge_discretization, file_name);

  // Read masses
  lstart = scanToFlag(fmem, "MASS", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<double> tmp_masses = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                    atom_count);

  // Read atom type (Lennard-Jones type) indices
  lstart = scanToFlag(fmem, "ATOM_TYPE_INDEX", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<int> tmp_lennard_jones_indices = iAmberPrmtopData(fmem, lstart, dfmt[0].x,
                                                                dfmt[0].z, atom_count);
  addScalarToVector(&tmp_lennard_jones_indices, -1);

  // Read excluded atom counts and immediately convert to a prefix sum, then cap it to create
  // a bounds array indexing into the exclusion list.
  lstart = scanToFlag(fmem, "NUMBER_EXCLUDED_ATOMS", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<int> tmp_atom_exclusion_counts = iAmberPrmtopData(fmem, lstart, dfmt[0].x,
                                                                dfmt[0].z, atom_count);

  // Read non-bonded parameter indices
  lstart = scanToFlag(fmem, "NONBONDED_PARM_INDEX", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<int> tmp_nonbonded_parameter_index =
    iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, lj_type_count * lj_type_count);
  addScalarToVector(&tmp_nonbonded_parameter_index, -1);

  // Read residue labels
  lstart = scanToFlag(fmem, "RESIDUE_LABEL", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<char4> tmp_residue_names = c4AmberPrmtopData(fmem, lstart, dfmt[0].x, residue_count);

  // Read residue limits and cap the result with the total number of atoms
  lstart = scanToFlag(fmem, "RESIDUE_POINTER", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<int> tmp_residue_limits = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                         residue_count);
  addScalarToVector(&tmp_residue_limits, -1);
  tmp_residue_limits.push_back(atom_count);

  // Assemble the array of atom structural numbers
  lstart = scanToFlag(fmem, "ATOMIC_STRUCTURE_NUMBERS", &dfmt, TopologyRequirement::OPTIONAL,
                      lstart);
  std::vector<int> tmp_atom_struc_numbers;
  if (lstart >= 0) {
    tmp_atom_struc_numbers = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, atom_count);
  }
  else {
    tmp_atom_struc_numbers.resize(atom_count);
    for (int i = 0; i < atom_count; i++) {
      tmp_atom_struc_numbers[i] = i + 1;
    }
  }
  
  // Assemble the array of residue numbers for each atom
  lstart = scanToFlag(fmem, "RESIDUE_NUMBERS", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<int> tmp_residue_numbers(atom_count);
  if (lstart >= 0) {
    std::vector<int> stated_residue_numbers = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                               residue_count);
    for (int i = 0; i < residue_count; i++) {
      for (int j = tmp_residue_limits[i]; j < tmp_residue_limits[i + 1]; j++) {
        tmp_residue_numbers[j] = stated_residue_numbers[i];
      }
    }
  }
  else {
    for (int i = 0; i < residue_count; i++) {
      for (int j = tmp_residue_limits[i]; j < tmp_residue_limits[i + 1]; j++) {
        tmp_residue_numbers[j] = i + 1;
      }
    }
  }
  
  // Read bond parameters
  lstart = scanToFlag(fmem, "BOND_FORCE_CONSTANT", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<double> tmp_bond_stiffnesses =
    eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, bond_parameter_count);
  lstart = scanToFlag(fmem, "BOND_EQUIL_VALUE", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<double> tmp_bond_equilibria =
    eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, bond_parameter_count);

  // Read angle parameters
  lstart = scanToFlag(fmem, "ANGLE_FORCE_CONSTANT", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<double> tmp_angl_stiffnesses =
    eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, angl_parameter_count);
  lstart = scanToFlag(fmem, "ANGLE_EQUIL_VALUE", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<double> tmp_angl_equilibria =
    eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, angl_parameter_count);

  // Read CHARMM Urey-Bradley angle terms and parameters
  std::vector<int> tmp_ub_atoms;
  std::vector<double> tmp_ub_stiffnesses;
  std::vector<double> tmp_ub_equilibria;
  lstart = scanToFlag(fmem, "CHARMM_UREY_BRADLEY_COUNT", &dfmt, TopologyRequirement::OPTIONAL,
                      lstart);
  if (lstart >= 0) {
    std::vector<int> tmp_ub_counts = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, 2);
    urey_bradley_term_count = tmp_ub_counts[0];
    urey_bradley_parameter_count = tmp_ub_counts[1];
    if (urey_bradley_term_count > 0) {
      lstart = scanToFlag(fmem, "CHARMM_UREY_BRADLEY", &dfmt, TopologyRequirement::ESSENTIAL,
                          lstart);
      tmp_ub_atoms = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                      urey_bradley_term_count * 3);
      lstart = scanToFlag(fmem, "CHARMM_UREY_BRADLEY_FORCE_CONSTANT", &dfmt,
                          TopologyRequirement::ESSENTIAL, lstart);
      tmp_ub_stiffnesses = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                            urey_bradley_parameter_count);
      lstart = scanToFlag(fmem, "CHARMM_UREY_BRADLEY_EQUIL_VALUE", &dfmt,
                          TopologyRequirement::ESSENTIAL, lstart);
      tmp_ub_equilibria = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                           urey_bradley_parameter_count);
    }
  }

  // Read dihedral parameters
  lstart = scanToFlag(fmem, "DIHEDRAL_FORCE_CONSTANT", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<double> tmp_dihe_amplitudes =
    eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, dihe_parameter_count);
  lstart = scanToFlag(fmem, "DIHEDRAL_PERIODICITY", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<double> tmp_dihe_periodicities =
    eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, dihe_parameter_count);
  lstart = scanToFlag(fmem, "DIHEDRAL_PHASE", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<double> tmp_dihe_phase_angles =
    eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, dihe_parameter_count);
  lstart = scanToFlag(fmem, "SCEE_SCALE_FACTOR", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<double> tmp_dihe_elec_screenings;
  if (lstart >= 0) {
    tmp_dihe_elec_screenings = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                dihe_parameter_count);
  }
  lstart = scanToFlag(fmem, "SCNB_SCALE_FACTOR", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<double> tmp_dihe_vdw_screenings;
  if (lstart >= 0) {
    tmp_dihe_vdw_screenings = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                               dihe_parameter_count);
  }

  // Read CHARMM improper dihedral parameters
  std::vector<int> tmp_charmm_impr_atoms;
  std::vector<double> tmp_charmm_impr_stiffnesses;
  std::vector<double> tmp_charmm_impr_phase_angles;
  lstart = scanToFlag(fmem, "CHARMM_NUM_IMPROPERS", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  if (lstart >= 0) {
    std::vector<int> tmp_cimpr_count = iAmberPrmtopData(fmem, lstart, 1, 8, 1);
    charmm_impr_term_count = tmp_cimpr_count[0];
    if (charmm_impr_term_count > 0) {
      lstart = scanToFlag(fmem, "CHARMM_IMPROPERS", &dfmt, TopologyRequirement::ESSENTIAL,
                          lstart);
      tmp_charmm_impr_atoms = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                               charmm_impr_term_count * 5);
      lstart = scanToFlag(fmem, "CHARMM_NUM_IMPR_TYPES", &dfmt, TopologyRequirement::ESSENTIAL,
                          lstart);
      std::vector<int> tmp_cimpr_types = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, 1);
      charmm_impr_parameter_count = tmp_cimpr_types[0];
      lstart = scanToFlag(fmem, "CHARMM_IMPROPER_FORCE_CONSTANT", &dfmt,
                          TopologyRequirement::ESSENTIAL, lstart);
      tmp_charmm_impr_stiffnesses = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                     charmm_impr_parameter_count);
      lstart = scanToFlag(fmem, "CHARMM_IMPROPER_PHASE", &dfmt, TopologyRequirement::ESSENTIAL,
                          lstart);
      tmp_charmm_impr_phase_angles = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                      charmm_impr_parameter_count);
    }
  }

  // Read hydrogen bonding parameters
  lstart = scanToFlag(fmem, "SOLTY", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<double> tmp_solty_info;
  if (lstart >= 0) {
    tmp_solty_info = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, unused_natyp);
  }

  // Read Lennard-Jones coefficients
  const int nljparm = lj_type_count * (lj_type_count + 1) / 2;
  lstart = scanToFlag(fmem, "LENNARD_JONES_ACOEF", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<double> tmp_lj_a_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                         nljparm);
  lstart = scanToFlag(fmem, "LENNARD_JONES_BCOEF", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<double> tmp_lj_b_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                         nljparm);
  lstart = scanToFlag(fmem, "LENNARD_JONES_CCOEF", &dfmt, TopologyRequirement::OPTIONAL,
                      lstart);
  std::vector<double> tmp_lj_c_values;
  if (lstart >= 0) {
    tmp_lj_c_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, nljparm);
  }
  std::vector<double> tmp_lj_14_a_values;
  std::vector<double> tmp_lj_14_b_values;
  std::vector<double> tmp_lj_14_c_values;
  lstart = scanToFlag(fmem, "LENNARD_JONES_14_ACOEF", &dfmt, TopologyRequirement::OPTIONAL,
                      lstart);
  if (lstart >= 0) {
    tmp_lj_14_a_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, nljparm);
  }
  else {
    tmp_lj_14_a_values = tmp_lj_a_values;
  }
  lstart = scanToFlag(fmem, "LENNARD_JONES_14_BCOEF", &dfmt, TopologyRequirement::OPTIONAL,
                      lstart);
  if (lstart >= 0) {
    tmp_lj_14_b_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, nljparm);
  }
  else {
    tmp_lj_14_b_values = tmp_lj_b_values;
  }
  lstart = scanToFlag(fmem, "LENNARD_JONES_14_CCOEF", &dfmt, TopologyRequirement::OPTIONAL,
                      lstart);
  if (lstart >= 0) {
    tmp_lj_14_c_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, nljparm);
  }
  else {
    tmp_lj_14_c_values = tmp_lj_c_values;
  }

  // Read angle, and dihedral indexing data.  Bond indexing data was read nearer the beginning in
  // order to calculate or confirm the number of separate molecules in the system.
  lstart = scanToFlag(fmem, "BONDS_INC_HYDROGEN", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<int> tmp_bond_atoms_h =
    iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, bond_term_with_hydrogen * 3);
  lstart = scanToFlag(fmem, "BONDS_WITHOUT_HYDROGEN", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<int> tmp_bond_atoms_noh =
    iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, bond_term_without_hydrogen * 3);
  lstart = scanToFlag(fmem, "ANGLES_INC_HYDROGEN", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<int> tmp_angl_atoms_h =
    iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, angl_term_with_hydrogen * 4);
  lstart = scanToFlag(fmem, "ANGLES_WITHOUT_HYDROGEN", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<int> tmp_angl_atoms_noh =
    iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, angl_term_without_hydrogen * 4);
  lstart = scanToFlag(fmem, "DIHEDRALS_INC_HYDROGEN", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<int> tmp_dihe_atoms_h =
    iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, dihe_term_with_hydrogen * 5);
  lstart = scanToFlag(fmem, "DIHEDRALS_WITHOUT_HYDROGEN", &dfmt, TopologyRequirement::ESSENTIAL,
                      lstart);
  std::vector<int> tmp_dihe_atoms_noh =
    iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, dihe_term_without_hydrogen * 5);

  // Read the excluded atoms list (despite counts being read a very long way up in the file).
  lstart = scanToFlag(fmem, "EXCLUDED_ATOMS_LIST", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<int> tmp_exclusion_list = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                         total_exclusions);

  // Read more hydrogen bonding parameters
  std::vector<double> tmp_hbond_a_values;
  std::vector<double> tmp_hbond_b_values;
  std::vector<double> tmp_hbond_cutoffs;
  if (hbond_10_12_parameter_count > 0) {
    lstart = scanToFlag(fmem, "HBOND_ACOEF", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
    tmp_hbond_a_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                          hbond_10_12_parameter_count);
    lstart = scanToFlag(fmem, "HBOND_BCOEF", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
    tmp_hbond_b_values = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                          hbond_10_12_parameter_count);
    lstart = scanToFlag(fmem, "HBCUT", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
    tmp_hbond_cutoffs = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                         hbond_10_12_parameter_count);
  }
  expandLennardJonesTables(&tmp_lj_a_values, &tmp_lj_b_values, &tmp_lj_c_values,
                           &tmp_lj_14_a_values, &tmp_lj_14_b_values, &tmp_lj_14_c_values,
                           &tmp_hbond_a_values, &tmp_hbond_b_values, lj_type_count,
                           tmp_nonbonded_parameter_index);

  // Read the atom type names
  lstart = scanToFlag(fmem, "AMBER_ATOM_TYPE", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
  std::vector<char4> tmp_atom_types = c4AmberPrmtopData(fmem, lstart, dfmt[0].x, atom_count);

  // Read the tree chain information
  lstart = scanToFlag(fmem, "TREE_CHAIN_CLASSIFICATION", &dfmt, TopologyRequirement::OPTIONAL,
                      lstart);
  std::vector<char4> tmp_tree_symbols;
  if (lstart >= 0) {
    tmp_tree_symbols = c4AmberPrmtopData(fmem, lstart, dfmt[0].x, atom_count);
  }
  lstart = scanToFlag(fmem, "JOIN_ARRAY", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<int> tmp_tree_joining_info;
  if (lstart >= 0) {
    tmp_tree_joining_info = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, atom_count);
  }
  lstart = scanToFlag(fmem, "IROTAT", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<int> tmp_last_rotator_info;
  if (lstart >= 0) {
    tmp_last_rotator_info = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, atom_count);
  }

  // Read the GB radius set and other parameters.  These are always allocated, to prepare for
  // applying a radius set as part of the implicit solvent model later.
  std::vector<double> tmp_atomic_pb_radii(atom_count, 0.0);
  std::vector<double> tmp_gb_screening_factors(atom_count, 0.0);
  std::vector<double> tmp_gb_coef(atom_count, 0.0);
  std::vector<int> tmp_neck_gb_indices(atom_count, 0);
  lstart = scanToFlag(fmem, "RADIUS_SET", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  if (lstart >= 0) {
    for (int i = tfr.line_limits[lstart]; i < tfr.line_limits[lstart + 1]; i++) {
      pb_radii_set += tfr.text[i];
    }
    lstart = scanToFlag(fmem, "RADII", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
    tmp_atomic_pb_radii = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, atom_count);
    lstart = scanToFlag(fmem, "SCREEN", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
    tmp_gb_screening_factors = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, atom_count);
  }

  // Read an extra descriptor for polarizability
  lstart = scanToFlag(fmem, "IPOL", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  std::vector<int> tmp_ipol;
  std::vector<double> tmp_atpol;
  if (lstart >= 0) {
    tmp_ipol = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, 1);

    // Proceed to read more information about polarizabilities
    if (tmp_ipol[0] == 1) {
      use_polarization = PolarizationSetting::ON;
      lstart = scanToFlag(fmem, "POLARIZABILITY", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
      tmp_atpol = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, atom_count);
    }
  }

  // Read virtual site frame information
  std::vector<int> vsite_custom_frames;
  std::vector<double> vsite_custom_details;
  lstart = scanToFlag(fmem, "VIRTUAL_SITE_FRAMES", &dfmt, TopologyRequirement::OPTIONAL, lstart);
  if (lstart >= 0) {
    vsite_custom_frames = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, 0,
                                           virtual_site_count * 6);
    
    // Convert the virtual site frame type enumerations in the topology to the internal
    // representation.
    int vs_est_count = static_cast<int>(vsite_custom_frames.size()) / 6;
    vs_est_count -= (vs_est_count * 6 < static_cast<int>(vsite_custom_frames.size()));
    for (int i = 0; i < vs_est_count; i++) {
      const int frame_slot = (i * 6) + 5;
      switch (vsite_custom_frames[frame_slot]) {
      case 0:
      case 1:
      case 2:
      case 3:

        // Cases 0 - 3 are reserved for innate Amber virtual site types and do not show up in
        // Amber topology customized VIRTUAL_SITE_FRAMES format fields.
        rtErr("Virtual site types 0 - 3 are reserved for Amber's inferred virtual site types.  "
              "No such number (" + std::to_string(vsite_custom_frames[i]) + ") should be present "
              "in an Amber topology file's custom VIRTUAL_SITE_FRAMES field, as was found in " +
              "virtual site " + std::to_string((i / 6) + 1) + " of " + source + ".", "AtomGraph",
              "buildAmberPrmtop");
      case 4:
        vsite_custom_frames[frame_slot] = static_cast<int>(VirtualSiteKind::FLEX_2);
        break;
      case 5:
        vsite_custom_frames[frame_slot] = static_cast<int>(VirtualSiteKind::FLEX_3);
        break;
      case 6:
        vsite_custom_frames[frame_slot] = static_cast<int>(VirtualSiteKind::FIXED_3);
        break;
      case 7:
        vsite_custom_frames[frame_slot] = static_cast<int>(VirtualSiteKind::FAD_3);
        break;
      case 8:
        vsite_custom_frames[frame_slot] = static_cast<int>(VirtualSiteKind::OUT_3);
        break;
      case 9:
        vsite_custom_frames[frame_slot] = static_cast<int>(VirtualSiteKind::FIXED_2);
        break;
      case 10:
        vsite_custom_frames[frame_slot] = static_cast<int>(VirtualSiteKind::FIXED_4);
        break;
      }
    }
  }
  lstart = scanToFlag(fmem, "VIRTUAL_SITE_FRAME_DETAILS", &dfmt,
                      TopologyRequirement::OPTIONAL, lstart);
  if (lstart >= 0) {
    vsite_custom_details = eAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, 0,
                                            virtual_site_count * 3);
  }

  // Read CHARMM CMAP terms and parameters
  std::vector<int> tmp_cmap_atoms;
  std::vector<int> tmp_cmap_surf_dims;
  std::vector<int> tmp_cmap_surf_bounds;
  std::vector<double> tmp_cmap_surfaces;
  std::vector<std::string> cmap_aliases = {"CHARMM_CMAP_", "CMAP_"};
  const int n_alias = cmap_aliases.size();
  for (int i = 0; i < n_alias; i++) {
    const std::string count_flag = cmap_aliases[i] + "COUNT";
    lstart = scanToFlag(fmem, count_flag.c_str(), &dfmt, TopologyRequirement::OPTIONAL, lstart);
    if (lstart >= 0) {
      const std::string res_flag = cmap_aliases[i] + "RESOLUTION";
      const std::string parm_flag = cmap_aliases[i] + "PARAMETER_";
      const std::string index_flag = cmap_aliases[i] + "INDEX";
      const std::vector<int> tmp_cmap_count = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                                               2);
      cmap_term_count = tmp_cmap_count[0];
      cmap_surface_count = tmp_cmap_count[1];
      lstart = scanToFlag(fmem, res_flag.c_str(), &dfmt, TopologyRequirement::ESSENTIAL, lstart);
      tmp_cmap_surf_dims = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z,
                                            cmap_surface_count);
      tmp_cmap_surf_bounds.resize(cmap_surface_count + 1);
      int cmap_bound_acc = 0;
      for (int i = 0; i < cmap_surface_count; i++) {
        tmp_cmap_surf_bounds[i] = cmap_bound_acc;
        cmap_bound_acc += tmp_cmap_surf_dims[i] * tmp_cmap_surf_dims[i];
        std::string cmap_name = parm_flag;
        cmap_name += (i + 1 < 10) ? "0" + std::to_string(i + 1) : std::to_string(i + 1);
        lstart = scanToFlag(fmem, cmap_name.c_str(), &dfmt, TopologyRequirement::ESSENTIAL,
                            lstart);
        const int npts = tmp_cmap_surf_dims[i] * tmp_cmap_surf_dims[i];

        // The CMAPs are read in row-major format and must be converted to column-major format.
        std::vector<double> tmp_surf = dAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, npts);
        for (int j = 0; j < tmp_cmap_surf_dims[i]; j++) {
          for (int k = 0; k < j; k++) {
            std::swap(tmp_surf[(j * tmp_cmap_surf_dims[i]) + k],
                      tmp_surf[(k * tmp_cmap_surf_dims[i]) + j]);
          }
        }
        tmp_cmap_surfaces.insert(tmp_cmap_surfaces.end(), tmp_surf.begin(), tmp_surf.end());
      }
      tmp_cmap_surf_bounds[cmap_surface_count] = cmap_bound_acc;
      lstart = scanToFlag(fmem, index_flag.c_str(), &dfmt, TopologyRequirement::ESSENTIAL,
                          lstart);
      tmp_cmap_atoms = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, 6 * cmap_term_count);
    }
  }
  CmapAccessories cmap_table = computeCmapDerivatives(cmap_surface_count, tmp_cmap_surf_dims,
                                                      tmp_cmap_surf_bounds, tmp_cmap_surfaces);

  // Condense the exclusion list to avoid counting atoms with no actual exclusions as having one
  // "blank atom" exclusion.  In the Amber topology, an atom with no exclusions is customarily
  // listed as having one, and the excluded atom is index 0, which in Fortran is not a valid array
  // index.
  CondensedExclusions cond_excl = processExclusions(tmp_atom_exclusion_counts, tmp_exclusion_list,
                                                    source);

  // Create tables of all valence interactions in preparation for organizing other connections
  // and aspects of the topology
  BasicValenceTable basic_vtable = basicValenceIndexing(atom_count, tmp_bond_atoms_h,
                                                        tmp_bond_atoms_noh, tmp_angl_atoms_h,
                                                        tmp_angl_atoms_noh, tmp_dihe_atoms_h,
                                                        tmp_dihe_atoms_noh);
  CharmmValenceTable charmm_vtable = charmmValenceIndexing(atom_count, tmp_ub_atoms,
                                                           tmp_charmm_impr_atoms, tmp_cmap_atoms);

  // Condense the non-bonded 1:4 scaling factors, and flag dihedrals that are responsible for no
  // such interaction.
  elec14_screening_factor = default_elec14_screening;
  vdw14_screening_factor = default_vdw14_screening;
  AttenuationParameterSet attn_parm =
    condenseScreeningFactors(basic_vtable, tmp_dihe_elec_screenings, tmp_dihe_vdw_screenings,
                             default_elec14_screening, default_vdw14_screening);
  attenuated_14_type_count = attn_parm.total_14_sets;

  // Create vectors for virtual site frame indexing and dimensions.  Glean the number of unique
  // virtual site frames from the result.
  VirtualSiteTable vsite_table = listVirtualSites(virtual_site_count, source,
                                                  tmp_masses, basic_vtable, tmp_bond_equilibria,
                                                  tmp_angl_equilibria, tmp_atom_types,
                                                  tmp_atom_names, vsite_custom_frames,
                                                  vsite_custom_details);
  virtual_site_parameter_set_count = vsite_table.frame_dim1.size();

  // Check the largest residue--it may have changed given the virtual site placement
  largest_residue_size = reviewLargestResidue(tmp_residue_limits, largest_residue_size, policy);

  // Elaborate on the bond connections information.  Make arrays of all 1:1 (virtual site to parent
  // atom), 1:2 (bonded atoms), 1:3 (atoms connected by a shortest path of two bonds), and 1:4
  // (atoms connected by a shortest path of three bonds) exclusions, double-counting everything to
  // list every excluded pair interaction, atom by atom.
  Map1234 all_nb_excl = mapExclusions(atom_count, basic_vtable, vsite_table);
  checkExclusions(cond_excl, all_nb_excl, source);

  // If the atomic numbers were not read from the topology file itself, infer them from the masses.
  // This can get tricky if there has been hydrogen mass repartitioning, but try and unroll that.
  if (tmp_atomic_numbers.size() == 0) {
    tmp_atomic_numbers = atomicNumbersFromMasses(tmp_masses, tmp_atom_names,
                                                 all_nb_excl.nb12_excl_list,
                                                 all_nb_excl.nb12_excl_bounds, source, policy);
  }

  // Examine dihedral coverge: are all 1:4 interactions covered by some dihedral, with a scaling
  // factor assoicated with those parameters?
  const std::vector<int3> outstanding_14_pairs =
    checkDihedral14Coverage(atom_count, tmp_atomic_numbers, basic_vtable, all_nb_excl, vsite_table,
                            attn_parm, policy);
  inferred_14_attenuations = outstanding_14_pairs.size();
  std::vector<int> tmp_inferred_14_i_atoms(inferred_14_attenuations);
  std::vector<int> tmp_inferred_14_l_atoms(inferred_14_attenuations);
  std::vector<int> tmp_inferred_14_param_idx(inferred_14_attenuations);
  for (int i = 0; i < inferred_14_attenuations; i++) {
    tmp_inferred_14_i_atoms[i]   = outstanding_14_pairs[i].x;
    tmp_inferred_14_l_atoms[i]   = outstanding_14_pairs[i].y;
    tmp_inferred_14_param_idx[i] = outstanding_14_pairs[i].z;
  }
  
  // With the bond connections information at the ready, calculate the number of molecules.  If
  // periodic box information is present, compare the information generated by the search to the
  // information presented in the topology itself.
  std::vector<int> tmp_molecule_membership, tmp_molecule_limits, tmp_molecule_contents;
  mapMolecules(atom_count, &molecule_count, all_nb_excl, &tmp_molecule_membership,
               &tmp_molecule_limits, &tmp_molecule_contents);
  std::vector<int> tmp_solvent_pointers;
  std::vector<int> tmp_atoms_per_molecule;
  switch(periodic_box_class) {
  case UnitCellType::NONE:
    last_solute_residue = residue_count - 1;
    last_solute_atom = atom_count - 1;
    first_solvent_molecule = molecule_count;
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:

    // Read the demarcations between solute and solvent
    lstart = scanToFlag(fmem, "SOLVENT_POINTERS", &dfmt, TopologyRequirement::ESSENTIAL, lstart);
    tmp_solvent_pointers = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, 3);
    last_solute_residue = tmp_solvent_pointers[0] - 1;
    last_solute_atom = (last_solute_residue == -1) ?
                       -1 : tmp_residue_limits[last_solute_residue + 1] - 1;
    first_solvent_molecule = tmp_solvent_pointers[2] - 1;

    // Read the numbers of atoms per molecule
    lstart = scanToFlag(fmem, "ATOMS_PER_MOLECULE", &dfmt, TopologyRequirement::ESSENTIAL,
                        lstart);
    tmp_atoms_per_molecule = iAmberPrmtopData(fmem, lstart, dfmt[0].x, dfmt[0].z, molecule_count);
    break;
  }
  for (int i = 0; i < molecule_count; i++) {
    largest_molecule_size = std::max(largest_molecule_size,
                                     tmp_molecule_limits[i + 1] - tmp_molecule_limits[i]);
  }
  
  // Mobile atoms are not stated in the topology, but they are an essential part of the information
  // on how the system will move.  Allocate a bitmask for all atoms to be mobile by default.
  const int mobile_atom_mask_size =
    roundUp<int>(std::max(1, atom_count / (static_cast<int>(sizeof(int)) * 8)), warp_size_int);
  const std::vector<int> tmp_mobile_atoms(mobile_atom_mask_size, -1);
  
  // Constraints are not yet known, but allocate space for any foreseeable constraint model (up to
  // and including all bonds to hydrogen, including those on water molecules, and all fast waters)
  // so that the topology does not have to be re-allocated later.
  const ConstraintTable cnst_table(tmp_atomic_numbers, tmp_masses, tmp_molecule_limits,
                                   tmp_molecule_contents, tmp_molecule_membership, basic_vtable,
                                   all_nb_excl, tmp_bond_equilibria, tmp_angl_equilibria);
  settle_group_count = cnst_table.settle_group_count;
  constraint_group_count = cnst_table.cnst_group_count;
  settle_parameter_count = cnst_table.settle_parameter_count;
  constraint_parameter_count = cnst_table.cnst_parameter_count;

  // Count the degrees of freedom, with and without constraints
  unconstrained_dof = (3 * (atom_count - vsite_table.vs_count));
  switch (periodic_box_class) {
  case UnitCellType::NONE:
    unconstrained_dof -= 6;
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    unconstrained_dof -= 3;
    break;
  }
  constrained_dof = unconstrained_dof - (3 * cnst_table.settle_group_count) -
                    (cnst_table.cnst_group_bounds[cnst_table.cnst_group_count] -
                     cnst_table.cnst_group_count);
  
  // Transfer data from the CPU-bound std::vectors and unguarded structs into HPC-capable Hybrid
  // objects.
  loadHybridArrays(tmp_desc, tmp_residue_limits, tmp_atom_struc_numbers, tmp_residue_numbers,
                   tmp_molecule_limits, tmp_atomic_numbers, tmp_molecule_membership,
                   tmp_mobile_atoms, tmp_molecule_contents, tmp_cmap_surf_dims,
                   tmp_cmap_surf_bounds, tmp_charge_indices, tmp_lennard_jones_indices,
                   tmp_inferred_14_i_atoms, tmp_inferred_14_l_atoms, tmp_inferred_14_param_idx,
                   tmp_neck_gb_indices, tmp_tree_joining_info, tmp_last_rotator_info, tmp_charges,
                   tmp_masses, tmp_ub_stiffnesses, tmp_ub_equilibria, tmp_charmm_impr_stiffnesses,
                   tmp_charmm_impr_phase_angles, tmp_cmap_surfaces, tmp_bond_stiffnesses,
                   tmp_bond_equilibria, tmp_angl_stiffnesses, tmp_angl_equilibria,
                   tmp_dihe_amplitudes, tmp_dihe_periodicities, tmp_dihe_phase_angles,
                   tmp_charge_parameters, tmp_lj_a_values, tmp_lj_b_values, tmp_lj_c_values,
                   tmp_lj_14_a_values, tmp_lj_14_b_values, tmp_lj_14_c_values,
                   tmp_atomic_pb_radii, tmp_gb_screening_factors, tmp_gb_coef, tmp_solty_info,
                   tmp_hbond_a_values, tmp_hbond_b_values, tmp_hbond_cutoffs, tmp_atom_names,
                   tmp_atom_types, tmp_residue_names, tmp_tree_symbols, cmap_table,
                   cond_excl, basic_vtable, charmm_vtable, attn_parm, vsite_table, all_nb_excl,
                   cnst_table);
}

} // namespace topology
} // namespace stormm
