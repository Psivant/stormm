#include <ctime>
#include <unistd.h>
#include "copyright.h"
#include "Parsing/parse.h"
#include "Math/vector_ops.h"
#include "atomgraph.h"

namespace stormm {
namespace topology {

using parse::char4ToString;
using stmath::enumerateMask;
using stmath::getSubsetIndexPattern;
using stmath::extractIndexedValues;
using stmath::prefixSumInPlace;
using stmath::PrefixSumType;
  
//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const AtomGraph &original, const std::vector<int> &atom_subset,
                     const ExceptionResponse policy, const double charge_rounding_tol,
                     const double charge_discretization,
                     const ApplyConstraints use_bond_constraints_in,
                     const ApplyConstraints use_settle_in, const HybridFormat format_in) :
    AtomGraph(format_in)
{
  // Add file header information to the object
  snprintf(version_stamp, 16, "STORMM 0.2");
  std::time_t raw_time = std::time(nullptr);
  std::tm* current_time = std::localtime(&raw_time);
  date = *current_time;
  title = "Extraction from an AtomGraph topology in STORMM";
  source = original.source;

  // Sort the subset in ascending order
  std::vector<int> local_subset(atom_subset);
  std::sort(local_subset.begin(), local_subset.end(), [](int a, int b) { return a < b; });
  
  // The number of atoms is the first thing that can be known.  Load all properties of atoms.
  const int nsubset = local_subset.size();
  std::vector<bool> subset_valid(nsubset, true);
  int nskip = 0;
  for (int i = 0; i < nsubset; i++) {
    if (local_subset[i] < 0 || local_subset[i] >= original.atom_count) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Subset index " + std::to_string(i) + ", for atom " +
              std::to_string(local_subset[i]) + ", is invalid for a topology containing " +
              std::to_string(original.atom_count) + " atoms.", "AtomGraph");
      case ExceptionResponse::WARN:
        rtWarn("Subset index " + std::to_string(i) + ", for atom " +
               std::to_string(local_subset[i]) + ", is invalid for a topology containing " +
               std::to_string(original.atom_count) + " atoms.  This entry will be skipped",
               "AtomGraph");
        nskip++;
        break;
      case ExceptionResponse::SILENT:
        nskip++;
        break;
      }
      subset_valid[i] = false;
    }
  }
  int valid_atom_counter = 0;
  for (int i = 0; i < nsubset; i++) {
    if (subset_valid[i]) {
      local_subset[valid_atom_counter] = local_subset[i];
      valid_atom_counter++;
    }
  }
  local_subset.resize(valid_atom_counter);
  atom_count = valid_atom_counter;
  if (atom_count == 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("No valid atoms were found in the subset of " + std::to_string(nsubset) + ", applied "
            "to topology " + original.source + " with " + std::to_string(original.atom_count) +
            " atoms.", "AtomGraph");
    case ExceptionResponse::WARN:
      rtWarn("No valid atoms were found in the subset of " + std::to_string(nsubset) + ", applied "
             "to topology " + original.source + " with " + std::to_string(original.atom_count) +
             " atoms.  An empty object will be returned.", "AtomGraph");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  force_fields = original.force_fields;
  const int* orig_rlim_ptr = original.residue_limits.data();
  std::vector<bool> atom_present(original.atom_count, false);
  for (int i = 0; i < atom_count; i++) {
    atom_present[atom_subset[i]] = true;
  }
  std::vector<int> residue_carryover(original.residue_count, 0);
  for (int i = 0; i < original.residue_count; i++) {
    for (int j = orig_rlim_ptr[i]; j < orig_rlim_ptr[i + 1]; j++) {
      residue_carryover[i] += static_cast<int>(atom_present[j]);
    }
    if (residue_carryover[i] > 0) {
      residue_count += 1;
    }
  }
  std::vector<int> tmp_residue_limits(residue_count + 1, 0);
  int rcon = 0;
  for (int i = 0; i < original.residue_count; i++) {
    if (residue_carryover[i] > 0) {
      tmp_residue_limits[rcon + 1] = tmp_residue_limits[rcon] + residue_carryover[i];
      largest_residue_size = std::max(largest_residue_size, residue_carryover[i]);
      rcon++;
    }
  }
  last_solute_residue = 0;
  for (int i = 0; i <= original.last_solute_residue; i++) {
    if (residue_carryover[i] > 0) {
      last_solute_residue += 1;
    }
  }
  last_solute_residue -= (last_solute_residue > 0);
  last_solute_atom = tmp_residue_limits[last_solute_residue + 1] - 1;
  std::vector<bool> molecule_present(original.molecule_count, false);
  const ChemicalDetailsKit cdk = original.getChemicalDetailsKit();
  for (int i = 0; i < atom_count; i++) {
    molecule_present[cdk.mol_home[local_subset[i]]] = true;
  }
  for (size_t i = 0; i < original.molecule_count; i++) {
    molecule_count += static_cast<int>(molecule_present[i]);
  }

  // Get the relevant abstracts to help with pointers
  const ChemicalDetailsKit orig_cdk         = original.getChemicalDetailsKit();
  const NonbondedKit<double> orig_nbk       = original.getDoublePrecisionNonbondedKit();
  const ValenceKit<double> orig_vk          = original.getDoublePrecisionValenceKit();
  const ImplicitSolventKit<double> orig_isk = original.getDoublePrecisionImplicitSolventKit();
  const VirtualSiteKit orig_vsk             = original.getDoublePrecisionVirtualSiteKit();

  // Check for virtual sites whose frames are not completely represented in the subset.  Die, or
  // reject the virtual sites if their frames are incomplete, depending on the fault tolerance
  // setting.
  bool vs_in_ascending_order  = true;
  bool vs_in_descending_order = true;
  const int* vs_indices = original.virtual_site_atoms.data();
  for (int i = 0; i < original.virtual_site_count - 1; i++) {
    vs_in_ascending_order  = (vs_in_ascending_order  && vs_indices[i] < vs_indices[i + 1]);
    vs_in_descending_order = (vs_in_descending_order && vs_indices[i] > vs_indices[i + 1]);
  }
  std::vector<int> missing_frame_atoms;
  for (int i = 0; i < orig_vsk.nsite; i++) {
    if (atom_present[orig_vsk.vs_atoms[i]] == false) {
      continue;
    }

    // This virtual site is present in the subset.  Make sure that the frame is completely
    // represented.
    bool problem = false;
    switch (static_cast<VirtualSiteKind>(orig_vsk.vs_types[i])) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
      if (atom_present[orig_vsk.frame1_idx[i]] == false ||
          atom_present[orig_vsk.frame2_idx[i]] == false) {
        problem = true;
        if (atom_present[orig_vsk.frame1_idx[i]] == false) {
          missing_frame_atoms.push_back(orig_vsk.frame1_idx[i]);
        }
        if (atom_present[orig_vsk.frame2_idx[i]] == false) {
          missing_frame_atoms.push_back(orig_vsk.frame2_idx[i]);
        }
      }
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      if (atom_present[orig_vsk.frame1_idx[i]] == false ||
          atom_present[orig_vsk.frame2_idx[i]] == false ||
          atom_present[orig_vsk.frame3_idx[i]] == false) {
        problem = true;
        if (atom_present[orig_vsk.frame1_idx[i]] == false) {
          missing_frame_atoms.push_back(orig_vsk.frame1_idx[i]);
        }
        if (atom_present[orig_vsk.frame2_idx[i]] == false) {
          missing_frame_atoms.push_back(orig_vsk.frame2_idx[i]);
        }
        if (atom_present[orig_vsk.frame3_idx[i]] == false) {
          missing_frame_atoms.push_back(orig_vsk.frame3_idx[i]);
        }
      }
      break;
    case VirtualSiteKind::FIXED_4:
      if (atom_present[orig_vsk.frame1_idx[i]] == false ||
          atom_present[orig_vsk.frame2_idx[i]] == false ||
          atom_present[orig_vsk.frame3_idx[i]] == false ||
          atom_present[orig_vsk.frame4_idx[i]] == false) {
        problem = true;
        if (atom_present[orig_vsk.frame1_idx[i]] == false) {
          missing_frame_atoms.push_back(orig_vsk.frame1_idx[i]);
        }
        if (atom_present[orig_vsk.frame2_idx[i]] == false) {
          missing_frame_atoms.push_back(orig_vsk.frame2_idx[i]);
        }
        if (atom_present[orig_vsk.frame3_idx[i]] == false) {
          missing_frame_atoms.push_back(orig_vsk.frame3_idx[i]);
        }
        if (atom_present[orig_vsk.frame4_idx[i]] == false) {
          missing_frame_atoms.push_back(orig_vsk.frame4_idx[i]);
        }
      }
      break;
    case VirtualSiteKind::NONE:
      break;
    }
    if (problem) {
      const int vs_atom_idx = orig_vsk.vs_atoms[i];
      const int res_idx = original.getResidueIndex(vs_atom_idx);
      std::string frame_atom_error_string;
      for (size_t i = 0; i < missing_frame_atoms.size(); i++) {
        frame_atom_error_string += original.getFullAtomName(missing_frame_atoms[i]);
        if (i < missing_frame_atoms.size() - 1) {
          frame_atom_error_string += ", ";
        }
      }
      rtErr("Virtual site " + char4ToString(cdk.res_names[res_idx]) + " " +
            std::to_string(cdk.res_numbers[vs_atom_idx]) + " " +
            char4ToString(cdk.atom_names[vs_atom_idx]) + " (topology index " +
            std::to_string(vs_atom_idx) + ") cannot be included in the topology subset without "
            "its frame atoms: " + frame_atom_error_string + ".", "AtomGraph");
    }
  }
  
  // Prepare a table of residue indices for the original topology
  const std::vector<int> base_residue_indices = original.getResidueIndex();
  std::vector<char4> tmp_residue_names(residue_count);
  std::vector<int> tmp_residue_index(residue_count);
  int res_con = 0;
  for (int i = 0; i < original.residue_count; i++) {
    if (residue_carryover[i] > 0) {
      tmp_residue_names[res_con] = orig_cdk.res_names[i];
      tmp_residue_index[res_con] = i;
      res_con++;
    }
  }

  // Prepare a mask of atoms in the subset, based on the original topology.  The map of atoms in
  // the original topology entering into the new topology is set to -1 to indicate that an atom is
  // not in the subset, and >= 0 to indicate the atom's index in the subset, hence its index in
  // the new topology.  All elements of the reciprocal map (atommap_new_to_orig) are greater than
  // or equal to zero, as every atom in the new topology has its original as some valid index in
  // the original topology.
  std::vector<int> atommap_orig_to_new(original.atom_count, -1);
  std::vector<int> atommap_new_to_orig(atom_count);
  std::vector<int> molmap_orig_to_new(original.molecule_count, -1);
  std::vector<int> molmap_new_to_orig(molecule_count);
  int mol_con = 0;
  for (int i = 0; i < original.molecule_count; i++) {
    if (molecule_present[i]) {
      molmap_orig_to_new[i] = mol_con;
      molmap_new_to_orig[mol_con] = i;
      mol_con++;
    }
  }
  
  // Allocate and tabulate all properties of atoms.
  const int nmask_elem = (atom_count + 31) / 32;
  std::vector<int> tmp_atom_struc_numbers(atom_count);
  std::vector<int> tmp_residue_numbers(atom_count);
  std::vector<int> tmp_atomic_numbers(atom_count);
  std::vector<int> tmp_molecule_membership(atom_count);
  std::vector<int> tmp_molecule_limits(molecule_count + 1, 0);
  std::vector<int> tmp_mobile_atoms(nmask_elem);
  std::vector<int> tmp_molecule_contents(atom_count);
  std::vector<int> tmp_lennard_jones_indices(atom_count);
  std::vector<int> tmp_neck_gb_indices(atom_count);
  std::vector<int> tmp_tree_joining_info(atom_count);
  std::vector<int> tmp_last_rotator_info(atom_count);
  std::vector<double> tmp_charges(atom_count), tmp_atomic_polarizabilities(atom_count);
  std::vector<double> tmp_masses(atom_count), tmp_inv_masses(atom_count);
  std::vector<double> tmp_atomic_pb_radii(atom_count);
  std::vector<double> tmp_gb_screening_factors(atom_count);
  std::vector<char4> tmp_atom_names(atom_count);
  std::vector<char4> tmp_atom_types(atom_count);
  std::vector<char4> tmp_tree_symbols(atom_count);
  for (int i = 0; i < atom_count; i++) {
    const int orig_idx = local_subset[i];
    tmp_atom_struc_numbers[i] = orig_cdk.atom_numbers[orig_idx];
    tmp_residue_numbers[i] = orig_cdk.res_numbers[orig_idx];
    tmp_atomic_numbers[i] = orig_cdk.z_numbers[orig_idx];
    tmp_lennard_jones_indices[i] = orig_nbk.lj_idx[orig_idx];
    tmp_neck_gb_indices[i] = orig_isk.neck_gb_idx[orig_idx];
    tmp_tree_joining_info[i] = original.tree_joining_info.readHost(orig_idx);
    tmp_last_rotator_info[i] = original.last_rotator_info.readHost(orig_idx);
    tmp_charges[i] = orig_nbk.charge[orig_idx];
    tmp_masses[i] = original.atomic_masses.readHost(orig_idx);
    tmp_inv_masses[i] = 1.0 / tmp_masses[i];
    tmp_atomic_pb_radii[i] = orig_isk.pb_radii[orig_idx];
    tmp_gb_screening_factors[i] = orig_isk.gb_screen[orig_idx];
    tmp_atomic_polarizabilities[i] = original.atomic_polarizabilities.readHost(orig_idx);
    tmp_atom_names[i] = orig_cdk.atom_names[orig_idx];
    tmp_atom_types[i] = orig_cdk.atom_types[orig_idx];
    tmp_tree_symbols[i] = original.tree_symbols.readHost(orig_idx);
    atommap_orig_to_new[orig_idx] = i;
    atommap_new_to_orig[i] = orig_idx;
    virtual_site_count += (orig_cdk.z_numbers[orig_idx] == 0);
  }
  for (int i = 0; i < atom_count; i++) {
    const int orig_idx = local_subset[i];
    tmp_molecule_membership[i] = molmap_orig_to_new[orig_cdk.mol_home[orig_idx]];
    tmp_molecule_limits[tmp_molecule_membership[i]] += 1;
  }
  prefixSumInPlace(&tmp_molecule_limits, PrefixSumType::EXCLUSIVE);
  std::vector<int> mol_counters = tmp_molecule_limits;
  for (int i = 0; i < atom_count; i++) {
    const int orig_atom_idx = local_subset[i];

    // Determine which molecule the atom was originally part of, then what molecule it is headed
    // to.  Increment the molecule counter to place subsequent atoms in the same molecule.
    const int orig_mol_idx = orig_cdk.mol_home[orig_atom_idx];
    const int new_mol_idx = molmap_orig_to_new[orig_mol_idx];
    tmp_molecule_contents[mol_counters[tmp_molecule_membership[i]]] = i;
    mol_counters[new_mol_idx] += 1;
  }
  for (int i = 0; i < nmask_elem; i++) {
    uint tmask = 0U;
    for (int j = 0; j < 32; j++) {
      if ((i * 32) + j < atom_count) {
        if (original.getAtomMobility(local_subset[(i * 32) + j])) {
          tmask |= (0x1U << j);
        }
      }
    }
    tmp_mobile_atoms[i] = static_cast<int>(tmask);
  }
  for (int i = 0; i < molecule_count; i++) {
    largest_molecule_size = std::max(largest_molecule_size,
                                     tmp_molecule_limits[i + 1] - tmp_molecule_limits[i]);
  }
  std::vector<double> tmp_charge_parameters;
  std::vector<int> tmp_charge_type_indices;
  smoothCharges(&tmp_charges, &tmp_charge_parameters, &tmp_charge_type_indices,
                &charge_type_count, charge_rounding_tol, charge_discretization, source);
  
  // Allocate and tabulate residue and molecule-level information.  Adjust the residue indices
  // to fit a pattern appropriate to the subset and prepare molecule-specific information.
  const std::vector<int> residue_correspondence = getSubsetIndexPattern(&tmp_residue_index,
                                                                        orig_cdk.nres,
                                                                        &residue_count);
  const std::vector<int> molecule_correspondence = getSubsetIndexPattern(&tmp_molecule_membership,
                                                                         orig_cdk.nmol,
                                                                         &molecule_count);

  // Collect virtual sites.  Any virtual sites from the original topology will have to come with
  // their frame atoms, but a complete set of frame atoms in the new topology does not carry a
  // requirement that the virtual site also be included.  The net charge of the subset topology
  // is the thing most likely to be affected, and it will not be required to be integral.
  std::vector<int> tmp_vs_atoms(virtual_site_count), tmp_vs_frame_types(virtual_site_count);
  std::vector<int> tmp_vs_frame1_idx(virtual_site_count), tmp_vs_frame2_idx(virtual_site_count);
  std::vector<int> tmp_vs_frame3_idx(virtual_site_count), tmp_vs_frame4_idx(virtual_site_count);
  std::vector<int> tmp_vs_parameter_indices(virtual_site_count);
  int vscon = 0;
  for (int i = 0; i < atom_count; i++) {
    if (tmp_atomic_numbers[i] == 0) {
      const int orig_atom_idx = local_subset[i];
      tmp_vs_atoms[vscon] = i;
      tmp_vs_frame_types[vscon] = original.virtual_site_frame_types.readHost(orig_atom_idx);
      const int orig_vsparm_idx = original.virtual_site_parameter_indices.readHost(orig_atom_idx);
      tmp_vs_parameter_indices[vscon] = orig_vsparm_idx;
      switch (static_cast<VirtualSiteKind>(tmp_vs_frame_types[vscon])) {
      case VirtualSiteKind::FLEX_2:
      case VirtualSiteKind::FIXED_2:
        tmp_vs_frame1_idx[vscon] = original.virtual_site_frame1_atoms.readHost(orig_atom_idx);
        tmp_vs_frame2_idx[vscon] = original.virtual_site_frame2_atoms.readHost(orig_atom_idx);
        break;
      case VirtualSiteKind::FLEX_3:
      case VirtualSiteKind::FIXED_3:
      case VirtualSiteKind::FAD_3:
      case VirtualSiteKind::OUT_3:
        tmp_vs_frame1_idx[vscon] = original.virtual_site_frame1_atoms.readHost(orig_atom_idx);
        tmp_vs_frame2_idx[vscon] = original.virtual_site_frame2_atoms.readHost(orig_atom_idx);
        tmp_vs_frame3_idx[vscon] = original.virtual_site_frame3_atoms.readHost(orig_atom_idx);
        break;
      case VirtualSiteKind::FIXED_4:
        tmp_vs_frame1_idx[vscon] = original.virtual_site_frame1_atoms.readHost(orig_atom_idx);
        tmp_vs_frame2_idx[vscon] = original.virtual_site_frame2_atoms.readHost(orig_atom_idx);
        tmp_vs_frame3_idx[vscon] = original.virtual_site_frame3_atoms.readHost(orig_atom_idx);
        tmp_vs_frame4_idx[vscon] = original.virtual_site_frame4_atoms.readHost(orig_atom_idx);
        break;
      case VirtualSiteKind::NONE:
        break;
      }
      vscon++;
    }
  }
  const std::vector<int> vsparm_correspondence =
    getSubsetIndexPattern(&tmp_vs_parameter_indices, orig_vsk.nframe_set,
                          &virtual_site_parameter_set_count);
  std::vector<double> tmp_vs_frame_dim1(virtual_site_parameter_set_count);
  std::vector<double> tmp_vs_frame_dim2(virtual_site_parameter_set_count);
  std::vector<double> tmp_vs_frame_dim3(virtual_site_parameter_set_count);
  for (int i = 0; i < virtual_site_parameter_set_count; i++) {
    const size_t vp_corr = vsparm_correspondence[i];
    tmp_vs_frame_dim1[i] = original.virtual_site_frame_dim1.readHost(vp_corr);
    tmp_vs_frame_dim2[i] = original.virtual_site_frame_dim2.readHost(vp_corr);
    tmp_vs_frame_dim3[i] = original.virtual_site_frame_dim3.readHost(vp_corr);
  }
  VirtualSiteTable vst(virtual_site_count, tmp_vs_atoms, tmp_vs_frame_types, tmp_vs_frame1_idx,
                       tmp_vs_frame2_idx, tmp_vs_frame3_idx, tmp_vs_frame4_idx,
                       tmp_vs_parameter_indices, tmp_vs_frame_dim1, tmp_vs_frame_dim2,
                       tmp_vs_frame_dim3);

  // Scan the original topology for each force field term and add the terms to the new topology
  // if all atoms are present in the subset.
  size_t new_bond_counter = 0;
  for (int pos = 0; pos < orig_vk.nbond; pos++) {
    if (atommap_orig_to_new[orig_vk.bond_i_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.bond_j_atoms[pos]] >= 0) {
      new_bond_counter++;
    }
  }
  size_t new_angl_counter = 0;
  for (int pos = 0; pos < orig_vk.nangl; pos++) {
    if (atommap_orig_to_new[orig_vk.angl_i_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.angl_j_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.angl_k_atoms[pos]] >= 0) {
      new_angl_counter++;
    }
  }
  size_t new_dihe_counter = 0;
  for (int pos = 0; pos < orig_vk.ndihe; pos++) {
    if (atommap_orig_to_new[orig_vk.dihe_i_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.dihe_j_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.dihe_k_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.dihe_l_atoms[pos]] >= 0) {
      new_dihe_counter++;
    }
  }
  bond_term_count = new_bond_counter;
  angl_term_count = new_angl_counter;
  dihe_term_count = new_dihe_counter;
  std::vector<int> tmp_bond_i_atoms(bond_term_count), tmp_bond_j_atoms(bond_term_count);
  std::vector<int> tmp_bond_param_idx(bond_term_count);
  std::vector<int> tmp_angl_i_atoms(angl_term_count), tmp_angl_j_atoms(angl_term_count);
  std::vector<int> tmp_angl_k_atoms(angl_term_count), tmp_angl_param_idx(angl_term_count);
  std::vector<int> tmp_dihe_i_atoms(dihe_term_count), tmp_dihe_j_atoms(dihe_term_count);
  std::vector<int> tmp_dihe_k_atoms(dihe_term_count), tmp_dihe_l_atoms(dihe_term_count);
  std::vector<int> tmp_dihe_param_idx(dihe_term_count), tmp_attn_param_idx(dihe_term_count);
  std::vector<char4> tmp_bond_mods(bond_term_count),  tmp_angl_mods(angl_term_count);
  std::vector<char4> tmp_dihe_mods(dihe_term_count);
  new_bond_counter = 0;
  for (int pos = 0; pos < orig_vk.nbond; pos++) {
    if (atommap_orig_to_new[orig_vk.bond_i_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.bond_j_atoms[pos]] >= 0) {
      tmp_bond_i_atoms[new_bond_counter]    = atommap_orig_to_new[orig_vk.bond_i_atoms[pos]];
      tmp_bond_j_atoms[new_bond_counter]    = atommap_orig_to_new[orig_vk.bond_j_atoms[pos]];
      tmp_bond_param_idx[new_bond_counter]  = orig_vk.bond_param_idx[pos];
      tmp_bond_mods[new_bond_counter] = orig_vk.bond_modifiers[pos];
      new_bond_counter++;
    }
  }
  new_angl_counter = 0;
  for (int pos = 0; pos < orig_vk.nangl; pos++) {
    if (atommap_orig_to_new[orig_vk.angl_i_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.angl_j_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.angl_k_atoms[pos]] >= 0) {
      tmp_angl_i_atoms[new_angl_counter]    = atommap_orig_to_new[orig_vk.angl_i_atoms[pos]];
      tmp_angl_j_atoms[new_angl_counter]    = atommap_orig_to_new[orig_vk.angl_j_atoms[pos]];
      tmp_angl_k_atoms[new_angl_counter]    = atommap_orig_to_new[orig_vk.angl_k_atoms[pos]];
      tmp_angl_param_idx[new_angl_counter]  = orig_vk.angl_param_idx[pos];
      tmp_angl_mods[new_angl_counter] = orig_vk.angl_modifiers[pos];
      new_angl_counter++;
    }
  }
  new_dihe_counter = 0;
  for (int pos = 0; pos < orig_vk.ndihe; pos++) {
    if (atommap_orig_to_new[orig_vk.dihe_i_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.dihe_j_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.dihe_k_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.dihe_l_atoms[pos]] >= 0) {
      tmp_dihe_i_atoms[new_dihe_counter]    = atommap_orig_to_new[orig_vk.dihe_i_atoms[pos]];
      tmp_dihe_j_atoms[new_dihe_counter]    = atommap_orig_to_new[orig_vk.dihe_j_atoms[pos]];
      tmp_dihe_k_atoms[new_dihe_counter]    = atommap_orig_to_new[orig_vk.dihe_k_atoms[pos]];
      tmp_dihe_l_atoms[new_dihe_counter]    = atommap_orig_to_new[orig_vk.dihe_l_atoms[pos]];
      tmp_dihe_param_idx[new_dihe_counter]  = orig_vk.dihe_param_idx[pos];
      tmp_attn_param_idx[new_dihe_counter]  = orig_vk.dihe14_param_idx[pos];
      tmp_dihe_mods[new_dihe_counter] = orig_vk.dihe_modifiers[pos];
      new_dihe_counter++;
    }
  }
  
  // These correspondence arrays have the length of the number of bond, angle, and dihedral
  // parameters in the original topology.  The kth element of each array indicates the
  // which parameter index of the new, smaller topology the kth bond, angle, or dihedral of the
  // original topology inspires.
  const std::vector<int> bond_correspondence = getSubsetIndexPattern(&tmp_bond_param_idx,
                                                                     orig_vk.nbond_param,
                                                                     &bond_parameter_count);
  const std::vector<int> angl_correspondence = getSubsetIndexPattern(&tmp_angl_param_idx,
                                                                     orig_vk.nangl_param,
                                                                     &angl_parameter_count);
  const std::vector<int> dihe_correspondence = getSubsetIndexPattern(&tmp_dihe_param_idx,
                                                                     orig_vk.ndihe_param,
                                                                     &dihe_parameter_count);
  const std::vector<int> attn_correspondence = getSubsetIndexPattern(&tmp_attn_param_idx,
                                                                     orig_vk.nattn14_param,
                                                                     &attenuated_14_type_count);
  
  // Compose the basic valence table and the CHARMM valence table.  Initially, the terms in the
  // new topology's holding arrays will bear the parameter indices of the old topology.  This
  // data will be used, in turn, to create the filtered parameter tables for the new topology
  // before updating the indexing tables.
  BasicValenceTable bvt(atom_count, bond_term_count, angl_term_count, dihe_term_count,
                        tmp_bond_i_atoms, tmp_bond_j_atoms, tmp_bond_param_idx, tmp_angl_i_atoms,
			tmp_angl_j_atoms, tmp_angl_k_atoms, tmp_angl_param_idx, tmp_dihe_i_atoms,
                        tmp_dihe_j_atoms, tmp_dihe_k_atoms, tmp_dihe_l_atoms, tmp_dihe_param_idx,
                        tmp_bond_mods, tmp_angl_mods, tmp_dihe_mods);
  bvt.checkBondAngleRelevance(original.use_shake, original.use_settle, tmp_atomic_numbers);

  // Create parameter tables for the basic valence terms, condensed for the new topology.
  const std::vector<double> tmp_bond_stiffnesses =
    extractIndexedValues(original.bond_stiffnesses, bond_correspondence, bond_parameter_count);
  const std::vector<double> tmp_bond_equilibria =
    extractIndexedValues(original.bond_equilibria, bond_correspondence, bond_parameter_count);
  const std::vector<double> tmp_angl_stiffnesses =
    extractIndexedValues(original.angl_stiffnesses, angl_correspondence, angl_parameter_count);
  const std::vector<double> tmp_angl_equilibria =
    extractIndexedValues(original.angl_equilibria, angl_correspondence, angl_parameter_count);
  const std::vector<double> tmp_dihe_amplitudes =
    extractIndexedValues(original.dihe_amplitudes, dihe_correspondence, dihe_parameter_count);
  const std::vector<double> tmp_dihe_periodicities =
    extractIndexedValues(original.dihe_periodicities, dihe_correspondence, dihe_parameter_count);
  const std::vector<double> tmp_dihe_phase_angles =
    extractIndexedValues(original.dihe_phase_angles, dihe_correspondence, dihe_parameter_count);

  // The 1:4 attenuation factors must be set up as if they came from a file.  Look up the 1:4
  // correspondence and fill the arrays as they would be read from a file containing one pair of
  // 1:4 scaling factors for every unique dihedral parameter set.
  std::vector<double> tmp_attn_14_elec_factors(dihe_parameter_count);
  std::vector<double> tmp_attn_14_vdw_factors(dihe_parameter_count);
  if (dihe_parameter_count > 0) {
    for (int i = 0; i < original.dihe_parameter_count; i++) {
      
      // Find the original topology's attenuation parameter set index.  Based on the dihedral
      // parameter set correspondence, set the values in the new topology.
      const int dihe_parm_idx = dihe_correspondence[i];
      const int attn14_idx = original.dihe14_parameter_indices.readHost(i);
      if (dihe_parm_idx >= 0) {
        tmp_attn_14_elec_factors[dihe_parm_idx] =
          original.attn14_elec_factors.readHost(attn14_idx);
        tmp_attn_14_vdw_factors[dihe_parm_idx] = original.attn14_vdw_factors.readHost(attn14_idx);
      }
    }
  }
  
  // Tabulate 1:4 attenuations.
  const AttenuationParameterSet attn_parm =
    condenseScreeningFactors(bvt, tmp_attn_14_elec_factors, tmp_attn_14_vdw_factors,
                             original.elec14_screening_factor, original.vdw14_screening_factor);
  attenuated_14_type_count = attn_parm.total_14_sets;
  
  // Finish up the basic valence parameter detailing with the atom assignments.  Count the number
  // of terms with and without hydrogen atoms.  Set the parameter indices of the new topology to
  // fit the condensed table.
  for (int pos = 0; pos < bond_term_count; pos++) {
    if (tmp_atomic_numbers[bvt.bond_i_atoms[pos]] == 1 ||
        tmp_atomic_numbers[bvt.bond_j_atoms[pos]] == 1) {
      bond_term_with_hydrogen += 1;
    }
    else {
      bond_term_without_hydrogen += 1;
    }
  }
  for (int pos = 0; pos < angl_term_count; pos++) {
    if (tmp_atomic_numbers[bvt.angl_i_atoms[pos]] == 1 ||
        tmp_atomic_numbers[bvt.angl_j_atoms[pos]] == 1 ||
        tmp_atomic_numbers[bvt.angl_k_atoms[pos]] == 1) {
      angl_term_with_hydrogen += 1;
    }
    else {
      angl_term_without_hydrogen += 1;
    }
  }
  for (int pos = 0; pos < dihe_term_count; pos++) {
    if (tmp_atomic_numbers[bvt.dihe_i_atoms[pos]] == 1 ||
        tmp_atomic_numbers[bvt.dihe_j_atoms[pos]] == 1 ||
        tmp_atomic_numbers[bvt.dihe_k_atoms[pos]] == 1 ||
        tmp_atomic_numbers[bvt.dihe_l_atoms[pos]] == 1) {
      dihe_term_with_hydrogen += 1;
    }
    else {
      dihe_term_without_hydrogen += 1;
    }
  }

  // Compose the CHARMM force field terms in a similar manner.
  size_t new_ubrd_counter = 0;
  for (int pos = 0; pos < orig_vk.nubrd; pos++) {
    if (atommap_orig_to_new[orig_vk.ubrd_i_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.ubrd_k_atoms[pos]] >= 0) {
      new_ubrd_counter++;
    }
  }
  size_t new_impr_counter = 0;
  for (int pos = 0; pos < orig_vk.ncimp; pos++) {
    if (atommap_orig_to_new[orig_vk.cimp_i_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.cimp_j_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.cimp_k_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.cimp_l_atoms[pos]] >= 0) {
      new_impr_counter++;
    }
  }
  size_t new_cmap_counter = 0;
  for (int pos = 0; pos < orig_vk.ncmap; pos++) {
    if (atommap_orig_to_new[orig_vk.cmap_i_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.cmap_j_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.cmap_k_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.cmap_l_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.cmap_m_atoms[pos]] >= 0) {
      new_cmap_counter++;
    }
  }
  urey_bradley_term_count = new_ubrd_counter;
  charmm_impr_term_count = new_impr_counter;
  cmap_term_count = new_cmap_counter;
  std::vector<int> tmp_ubrd_i_atoms(urey_bradley_term_count);
  std::vector<int> tmp_ubrd_k_atoms(urey_bradley_term_count);
  std::vector<int> tmp_ubrd_param_idx(urey_bradley_term_count);
  std::vector<int> tmp_impr_i_atoms(charmm_impr_term_count);
  std::vector<int> tmp_impr_j_atoms(charmm_impr_term_count);
  std::vector<int> tmp_impr_k_atoms(charmm_impr_term_count);
  std::vector<int> tmp_impr_l_atoms(charmm_impr_term_count);
  std::vector<int> tmp_impr_param_idx(charmm_impr_term_count);
  std::vector<int> tmp_cmap_i_atoms(cmap_term_count);
  std::vector<int> tmp_cmap_j_atoms(cmap_term_count);
  std::vector<int> tmp_cmap_k_atoms(cmap_term_count);
  std::vector<int> tmp_cmap_l_atoms(cmap_term_count);
  std::vector<int> tmp_cmap_m_atoms(cmap_term_count);
  std::vector<int> tmp_cmap_param_idx(cmap_term_count);
  new_ubrd_counter = 0;
  for (int pos = 0; pos < orig_vk.nubrd; pos++) {
    if (atommap_orig_to_new[orig_vk.ubrd_i_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.ubrd_k_atoms[pos]] >= 0) {
      tmp_ubrd_i_atoms[new_ubrd_counter] = atommap_orig_to_new[orig_vk.ubrd_i_atoms[pos]];
      tmp_ubrd_k_atoms[new_ubrd_counter] = atommap_orig_to_new[orig_vk.ubrd_k_atoms[pos]];
      tmp_ubrd_param_idx[new_ubrd_counter] = orig_vk.ubrd_param_idx[pos];
      new_ubrd_counter++;
    }
  }
  new_impr_counter = 0;
  for (int pos = 0; pos < orig_vk.ncimp; pos++) {
    if (atommap_orig_to_new[orig_vk.cimp_i_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.cimp_j_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.cimp_k_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.cimp_l_atoms[pos]] >= 0) {
      tmp_impr_i_atoms[new_impr_counter] = atommap_orig_to_new[orig_vk.cimp_i_atoms[pos]];
      tmp_impr_j_atoms[new_impr_counter] = atommap_orig_to_new[orig_vk.cimp_j_atoms[pos]];
      tmp_impr_k_atoms[new_impr_counter] = atommap_orig_to_new[orig_vk.cimp_k_atoms[pos]];
      tmp_impr_l_atoms[new_impr_counter] = atommap_orig_to_new[orig_vk.cimp_l_atoms[pos]];
      tmp_impr_param_idx[new_impr_counter] = orig_vk.cimp_param_idx[pos];
      new_impr_counter++;
    }
  }
  new_cmap_counter = 0;
  for (int pos = 0; pos < orig_vk.ncmap; pos++) {
    if (atommap_orig_to_new[orig_vk.cmap_i_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.cmap_j_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.cmap_k_atoms[pos]] >= 0 &&
        atommap_orig_to_new[orig_vk.cmap_l_atoms[pos]] >= 0) {
      tmp_cmap_i_atoms[new_cmap_counter] = atommap_orig_to_new[orig_vk.cmap_i_atoms[pos]];
      tmp_cmap_j_atoms[new_cmap_counter] = atommap_orig_to_new[orig_vk.cmap_j_atoms[pos]];
      tmp_cmap_k_atoms[new_cmap_counter] = atommap_orig_to_new[orig_vk.cmap_k_atoms[pos]];
      tmp_cmap_l_atoms[new_cmap_counter] = atommap_orig_to_new[orig_vk.cmap_l_atoms[pos]];
      tmp_cmap_m_atoms[new_cmap_counter] = atommap_orig_to_new[orig_vk.cmap_m_atoms[pos]];
      tmp_cmap_param_idx[new_cmap_counter] = orig_vk.cmap_surf_idx[pos];
      new_cmap_counter++;
    }
  }
  const std::vector<int> ubrd_correspondence =
    getSubsetIndexPattern(&tmp_ubrd_param_idx, orig_vk.nubrd_param, &urey_bradley_parameter_count);
  const std::vector<int> cimp_correspondence =
    getSubsetIndexPattern(&tmp_impr_param_idx, orig_vk.ncimp_param, &charmm_impr_parameter_count);
  const std::vector<int> cmap_correspondence =
    getSubsetIndexPattern(&tmp_cmap_param_idx, orig_vk.ncmap_surf, &cmap_surface_count);
  CharmmValenceTable mvt(atom_count, urey_bradley_term_count, charmm_impr_term_count,
                         cmap_term_count, tmp_ubrd_i_atoms, tmp_ubrd_k_atoms, tmp_ubrd_param_idx,
                         tmp_impr_i_atoms, tmp_impr_j_atoms, tmp_impr_k_atoms, tmp_impr_l_atoms,
                         tmp_impr_param_idx, tmp_cmap_i_atoms, tmp_cmap_j_atoms, tmp_cmap_k_atoms,
                         tmp_cmap_l_atoms, tmp_cmap_m_atoms, tmp_cmap_param_idx);

  // Make parameter tables for CHARMM force field terms
  const std::vector<double> tmp_urey_bradley_stiffnesses =
    extractIndexedValues(original.urey_bradley_stiffnesses, ubrd_correspondence,
                         urey_bradley_parameter_count);
  const std::vector<double> tmp_urey_bradley_equilibria =
    extractIndexedValues(original.urey_bradley_equilibria, ubrd_correspondence,
                         urey_bradley_parameter_count);
  const std::vector<double> tmp_charmm_impr_stiffnesses =
    extractIndexedValues(original.charmm_impr_stiffnesses, cimp_correspondence,
                         charmm_impr_parameter_count);
  const std::vector<double> tmp_charmm_impr_phase_angles =
    extractIndexedValues(original.charmm_impr_phase_angles, cimp_correspondence,
                         charmm_impr_parameter_count);
  cmap_term_count = cmap_correspondence.size();
  std::vector<int> tmp_cmap_surface_bounds(cmap_surface_count + 1, 0);
  std::vector<int> tmp_cmap_surface_dims;
  tmp_cmap_surface_dims.reserve(cmap_surface_count);
  for (int i = 0; i < cmap_term_count; i++) {
    if (cmap_correspondence[i] >= 0) {
      tmp_cmap_surface_bounds[cmap_correspondence[i]] = orig_vk.cmap_dim[i] * orig_vk.cmap_dim[i];
      tmp_cmap_surface_dims.push_back(orig_vk.cmap_dim[i]);
    }
  }
  prefixSumInPlace(&tmp_cmap_surface_bounds, PrefixSumType::EXCLUSIVE, "AtomGraph");
  const int cmap_alloc_size = tmp_cmap_surface_bounds[cmap_surface_count];
  std::vector<double> tmp_cmap_surfaces(cmap_alloc_size);
  for (int i = 0; i < cmap_term_count; i++) {
    if (cmap_correspondence[i] >= 0) {
      const int write_pos = tmp_cmap_surface_bounds[cmap_correspondence[i]];
      const int read_pos = orig_vk.cmap_surf_bounds[i];
      const int map_element_count = orig_vk.cmap_dim[i] * orig_vk.cmap_dim[i];
      for (int j = 0; j < map_element_count; j++) {
        tmp_cmap_surfaces[write_pos + j] = orig_vk.cmap_surf[read_pos + j];
      }
    }
  }
  const CmapAccessories cma = computeCmapDerivatives(cmap_term_count, tmp_cmap_surface_dims,
                                                     tmp_cmap_surface_bounds, tmp_cmap_surfaces);
  
  // Populate the non-bonded Lennard-Jones tables and re-number the new topology's indexing.  In
  // this map covering how the original indices map to the new topology's indices, the "x" member
  // indicates the original topology's Lennard-Jones index and the "y" member indicates the new
  // topology's Lennard-Jones index.  Also handle (deprecated) hydrogen bonding parameters.
  std::vector<int> ljmap_orig_to_new(original.lj_type_count, -1);
  std::vector<bool> lj_type_present(original.lj_type_count, false);
  const int* orig_ljt_ptr = original.lennard_jones_indices.data();
  for (int i = 0; i < original.atom_count; i++) {
    if (atom_present[i]) {
      lj_type_present[orig_ljt_ptr[i]] = true;
    }
  }
  int ntcon = 0;
  for (int i = 0; i < original.lj_type_count; i++) {
    if (lj_type_present[i]) {
      ljmap_orig_to_new[i] = ntcon;
      ntcon++;
    }
  }
  int* lj_idx_ptr = tmp_lennard_jones_indices.data();
  for (int i = 0; i < atom_count; i++) {
    lj_idx_ptr[i] = ljmap_orig_to_new[lj_idx_ptr[i]];
  }
  lj_type_count = ntcon;
  std::vector<double> tmp_lja_coef(lj_type_count * lj_type_count);
  std::vector<double> tmp_ljb_coef(lj_type_count * lj_type_count);
  std::vector<double> tmp_ljc_coef(lj_type_count * lj_type_count);
  std::vector<double> tmp_lja_14_coef(lj_type_count * lj_type_count);
  std::vector<double> tmp_ljb_14_coef(lj_type_count * lj_type_count);
  std::vector<double> tmp_ljc_14_coef(lj_type_count * lj_type_count);
  const std::vector<double> tmp_gb_coef(atom_count, 0.0);
  const std::vector<double> tmp_solty_info = original.solty_info.readHost();
  std::vector<double> tmp_hbond_a_values(lj_type_count * lj_type_count);
  std::vector<double> tmp_hbond_b_values(lj_type_count * lj_type_count);
  std::vector<double> tmp_hbond_cutoffs(lj_type_count * lj_type_count);
  for (int i = 0; i < original.lj_type_count; i++) {
    const int iljt = ljmap_orig_to_new[i];
    if (iljt < 0) {
      continue;
    }
    for (int j = 0; j < original.lj_type_count; j++) {
      const int jljt = ljmap_orig_to_new[j];
      if (jljt < 0) {
        continue;
      }
      const size_t ij_ljt = (iljt * lj_type_count) + jljt;
      const size_t orig_ij_ljt = (i * original.lj_type_count) + j;
      tmp_lja_coef[ij_ljt] = orig_nbk.lja_coeff[orig_ij_ljt];
      tmp_ljb_coef[ij_ljt] = orig_nbk.ljb_coeff[orig_ij_ljt];
      tmp_ljc_coef[ij_ljt] = orig_nbk.ljc_coeff[orig_ij_ljt];
      tmp_lja_14_coef[ij_ljt] = orig_nbk.lja_14_coeff[orig_ij_ljt];
      tmp_ljb_14_coef[ij_ljt] = orig_nbk.ljb_14_coeff[orig_ij_ljt];
      tmp_ljc_14_coef[ij_ljt] = orig_nbk.ljc_14_coeff[orig_ij_ljt];
      tmp_hbond_a_values[ij_ljt] = original.hbond_a_values.readHost(orig_ij_ljt);
      tmp_hbond_b_values[ij_ljt] = original.hbond_b_values.readHost(orig_ij_ljt);
      tmp_hbond_cutoffs[ij_ljt] = original.hbond_cutoffs.readHost(orig_ij_ljt);
    }
  }

  // Transfer non-bonded exclusions from the original topology: these might not be counted if they
  // are dervied from bonded interactions, as a 1:3 or 1:4 interaction may be lost without the
  // middle atoms.
  std::vector<int> tmp_raw_excl_counts(atom_count, 0), tmp_excl_bounds(atom_count + 1);
  const int* orig_excl_list_ptr = original.atom_exclusion_list.data();
  const int* orig_excl_bounds_ptr = original.atom_exclusion_bounds.data();
  for (int i = 0; i < atom_count; i++) {
    const int orig_atom_idx = local_subset[i];
    const int j_llim = orig_excl_bounds_ptr[orig_atom_idx];
    const int j_hlim = orig_excl_bounds_ptr[orig_atom_idx + 1];
    for (int j = j_llim; j < j_hlim; j++) {
      const int orig_partner = orig_excl_list_ptr[j];
      if (atommap_orig_to_new[orig_partner] >= 0) {
        tmp_raw_excl_counts[i] += 1;
      }
    }
    tmp_raw_excl_counts[i] = std::max(tmp_raw_excl_counts[i], 1);
    tmp_excl_bounds[i] = total_exclusions;
    total_exclusions += tmp_raw_excl_counts[i];
  }
  tmp_excl_bounds[atom_count] = total_exclusions;
  std::vector<int> tmp_raw_exclusions(total_exclusions, 0);
  std::vector<int> excl_counters = tmp_excl_bounds;
  for (int i = 0; i < atom_count; i++) {
    const int orig_atom_idx = local_subset[i];
    const int j_llim = orig_excl_bounds_ptr[orig_atom_idx];
    const int j_hlim = orig_excl_bounds_ptr[orig_atom_idx + 1];
    for (int j = j_llim; j < j_hlim; j++) {
      const int orig_partner = orig_excl_list_ptr[j];
      if (atommap_orig_to_new[orig_partner] >= 0) {
        tmp_raw_exclusions[excl_counters[i]] = atommap_orig_to_new[orig_partner] + 1;
        excl_counters[i] += 1;
      }
    }
  }
  const CondensedExclusions cond_excl = processExclusions(tmp_raw_excl_counts, tmp_raw_exclusions,
                                                          source);
  Map1234 all_nb_excl = mapExclusions(atom_count, tmp_vs_atoms, tmp_vs_frame1_idx,
                                      tmp_bond_i_atoms, tmp_bond_j_atoms);
  
  // All excluded interactions in the original topology should be reflected in the new topology,
  // by convention.  Additional steps would be required to remove a bond and the associated
  // exclusions if both atoms were otherwise present in the new topology.  Make a list of missing
  // exclusions so that they can be added.  Each exclusion is expected to remain of the same order,
  // as to change the order of an exclusion one would need to delete bonds connecting the atoms and
  // then form a new bond, which does not occur in this routine.
  std::vector<int3> nb11_missing, nb12_missing, nb13_missing, nb14_missing;
  for (int i = 0; i < orig_nbk.natom; i++) {
    if (atom_present[i]) {
      const int ci_atom_idx = atommap_orig_to_new[i];
      const int nb11_llim = all_nb_excl.nb11_excl_bounds[ci_atom_idx];
      const int nb11_hlim = all_nb_excl.nb11_excl_bounds[ci_atom_idx + 1];
      const int nb12_llim = all_nb_excl.nb12_excl_bounds[ci_atom_idx];
      const int nb12_hlim = all_nb_excl.nb12_excl_bounds[ci_atom_idx + 1];
      const int nb13_llim = all_nb_excl.nb13_excl_bounds[ci_atom_idx];
      const int nb13_hlim = all_nb_excl.nb13_excl_bounds[ci_atom_idx + 1];
      const int nb14_llim = all_nb_excl.nb14_excl_bounds[ci_atom_idx];
      const int nb14_hlim = all_nb_excl.nb14_excl_bounds[ci_atom_idx + 1];
      for (int j = orig_nbk.nb11_bounds[i]; j < orig_nbk.nb11_bounds[i + 1]; j++) {
        if (atom_present[orig_nbk.nb11x[j]]) {
          const int cj_atom_idx = atommap_orig_to_new[orig_nbk.nb11x[j]];
          if (ci_atom_idx < cj_atom_idx) {
            if (scanForExclusion(all_nb_excl.nb11_excl_list, nb11_llim,
                                 nb11_hlim, cj_atom_idx) == false &&
                scanForExclusion(all_nb_excl.nb12_excl_list, nb12_llim,
                                 nb12_hlim, cj_atom_idx) == false &&
                scanForExclusion(all_nb_excl.nb13_excl_list, nb13_llim,
                                 nb13_hlim, cj_atom_idx) == false &&
                scanForExclusion(all_nb_excl.nb14_excl_list, nb14_llim,
                                 nb14_hlim, cj_atom_idx) == false) {
              nb11_missing.push_back({ ci_atom_idx, cj_atom_idx, 2 });
            }
          }
        }
      }
      for (int j = orig_nbk.nb12_bounds[i]; j < orig_nbk.nb12_bounds[i + 1]; j++) {
        if (atom_present[orig_nbk.nb12x[j]]) {
          const int cj_atom_idx = atommap_orig_to_new[orig_nbk.nb12x[j]];
          if (ci_atom_idx < cj_atom_idx) {
            if (scanForExclusion(all_nb_excl.nb12_excl_list, nb12_llim,
                                 nb12_hlim, cj_atom_idx) == false &&
                scanForExclusion(all_nb_excl.nb11_excl_list, nb11_llim,
                                 nb11_hlim, cj_atom_idx) == false &&
                scanForExclusion(all_nb_excl.nb13_excl_list, nb13_llim,
                                 nb13_hlim, cj_atom_idx) == false &&
                scanForExclusion(all_nb_excl.nb14_excl_list, nb14_llim,
                                 nb14_hlim, cj_atom_idx) == false) {
              nb12_missing.push_back({ ci_atom_idx, cj_atom_idx, 2 });
            }
          }
        }
      }
      for (int j = orig_nbk.nb13_bounds[i]; j < orig_nbk.nb13_bounds[i + 1]; j++) {
        if (atom_present[orig_nbk.nb13x[j]]) {
          const int cj_atom_idx = atommap_orig_to_new[orig_nbk.nb13x[j]];
          if (ci_atom_idx < cj_atom_idx) {
            if (scanForExclusion(all_nb_excl.nb13_excl_list, nb13_llim,
                                 nb13_hlim, cj_atom_idx) == false &&
                scanForExclusion(all_nb_excl.nb11_excl_list, nb11_llim,
                                 nb11_hlim, cj_atom_idx) == false &&
                scanForExclusion(all_nb_excl.nb12_excl_list, nb12_llim,
                                 nb12_hlim, cj_atom_idx) == false &&
                scanForExclusion(all_nb_excl.nb14_excl_list, nb14_llim,
                                 nb14_hlim, cj_atom_idx) == false) {
              nb13_missing.push_back({ ci_atom_idx, cj_atom_idx, 2 });
            }
          }
        }
      }
      for (int j = orig_nbk.nb14_bounds[i]; j < orig_nbk.nb14_bounds[i + 1]; j++) {
        if (atom_present[orig_nbk.nb14x[j]]) {
          const int cj_atom_idx = atommap_orig_to_new[orig_nbk.nb14x[j]];
          if (ci_atom_idx < cj_atom_idx) {
            if (scanForExclusion(all_nb_excl.nb14_excl_list, nb14_llim,
                                 nb14_hlim, cj_atom_idx) == false &&
                scanForExclusion(all_nb_excl.nb11_excl_list, nb11_llim,
                                 nb11_hlim, cj_atom_idx) == false &&
                scanForExclusion(all_nb_excl.nb12_excl_list, nb12_llim,
                                 nb12_hlim, cj_atom_idx) == false &&
                scanForExclusion(all_nb_excl.nb13_excl_list, nb13_llim,
                                 nb13_hlim, cj_atom_idx) == false) {
              nb14_missing.push_back({ ci_atom_idx, cj_atom_idx, 2 });
            }
          }
        }
      }
    }
  }
  if (nb11_missing.size() > 0) {    
    all_nb_excl.addExclusions(nb11_missing, policy);
  }
  if (nb12_missing.size() > 0) {
    cullPriorExclusions(&nb12_missing, nb11_missing);
    all_nb_excl.addExclusions(nb12_missing, policy);
  }
  if (nb13_missing.size() > 0) {
    cullPriorExclusions(&nb13_missing, nb11_missing);
    cullPriorExclusions(&nb13_missing, nb12_missing);
    all_nb_excl.addExclusions(nb13_missing, policy);
  }
  if (nb14_missing.size() > 0) {
    cullPriorExclusions(&nb14_missing, nb11_missing);
    cullPriorExclusions(&nb14_missing, nb12_missing);
    cullPriorExclusions(&nb14_missing, nb13_missing);
    all_nb_excl.addExclusions(nb14_missing, policy);
  }

  // Assemble the constraints table.
  const ConstraintTable cnst_table(tmp_atomic_numbers, tmp_masses, tmp_molecule_limits,
                                   tmp_molecule_contents, tmp_molecule_membership, bvt,
                                   all_nb_excl, tmp_bond_equilibria, tmp_angl_equilibria);

  // Examine dihedral coverge: are all 1:4 interactions covered by some dihedral, with a scaling
  // factor associated with those parameters?
  const std::vector<int3> outstanding_14_pairs =
    checkDihedral14Coverage(atom_count, tmp_atomic_numbers, bvt, all_nb_excl, vst, attn_parm,
                            policy);
  inferred_14_attenuations = outstanding_14_pairs.size();
  std::vector<int> tmp_inferred_14_i_atoms(inferred_14_attenuations);
  std::vector<int> tmp_inferred_14_l_atoms(inferred_14_attenuations);
  std::vector<int> tmp_inferred_14_param_idx(inferred_14_attenuations);
  for (int i = 0; i < inferred_14_attenuations; i++) {
    tmp_inferred_14_i_atoms[i]   = outstanding_14_pairs[i].x;
    tmp_inferred_14_l_atoms[i]   = outstanding_14_pairs[i].y;
    tmp_inferred_14_param_idx[i] = outstanding_14_pairs[i].z;
  }

  // Carry over settings with rare or deprecated uses.
  implicit_copy_count = original.implicit_copy_count;
  use_solvent_cap_option = original.use_solvent_cap_option;
  unused_nhparm = original.unused_nhparm;
  unused_nparm = original.unused_nparm;
  unused_natyp = original.unused_natyp;
  hbond_10_12_parameter_count = original.hbond_10_12_parameter_count;
  heavy_bonds_plus_constraints = original.heavy_bonds_plus_constraints;
  heavy_angls_plus_constraints = original.heavy_angls_plus_constraints;
  heavy_dihes_plus_constraints = original.heavy_dihes_plus_constraints;

  // Perturbations are not permitted as part of a reduced topology.
  bond_perturbation_term_count = 0;
  angl_perturbation_term_count = 0;
  dihe_perturbation_term_count = 0;
  bonds_in_perturbed_group = 0;
  angls_in_perturbed_group = 0;
  dihes_in_perturbed_group = 0;
  
  // In contrast to topology combinations, the general descriptions of the topology can only take
  // shape once everything is populated.
  const int n_desc = static_cast<int>(TopologyDescriptor::N_VALUES);
  std::vector<int> tmp_desc(n_desc);
  for (int i = 0; i < n_desc; i++) {
    switch (static_cast<TopologyDescriptor>(i)) {
    case TopologyDescriptor::ATOM_COUNT:
      tmp_desc[i] = atom_count;
      break;
    case TopologyDescriptor::ATOM_TYPE_COUNT:
      tmp_desc[i] = lj_type_count;
      break;
    case TopologyDescriptor::BONDS_WITH_HYDROGEN:
      tmp_desc[i] = bond_term_with_hydrogen;
      break;
    case TopologyDescriptor::BONDS_WITHOUT_HYDROGEN:
      tmp_desc[i] = bond_term_without_hydrogen;
      break;
    case TopologyDescriptor::ANGLES_WITH_HYDROGEN:
      tmp_desc[i] = angl_term_with_hydrogen;
      break;
    case TopologyDescriptor::ANGLES_WITHOUT_HYDROGEN:
      tmp_desc[i] = angl_term_without_hydrogen;
      break;
    case TopologyDescriptor::DIHEDRALS_WITH_HYDROGEN:
      tmp_desc[i] = dihe_term_with_hydrogen;
      break;
    case TopologyDescriptor::DIHEDRALS_WITHOUT_HYDROGEN:
      tmp_desc[i] = dihe_term_without_hydrogen;
      break;
    case TopologyDescriptor::NHPARM_UNUSED:
      tmp_desc[i] = unused_nhparm;
      break;
    case TopologyDescriptor::ADDLES_CREATED:
      tmp_desc[i] = unused_nparm;
      break;
    case TopologyDescriptor::TOTAL_EXCLUDED_ATOMS:
      tmp_desc[i] = total_exclusions;
      break;
    case TopologyDescriptor::RESIDUE_COUNT:
      tmp_desc[i] = residue_count;
      break;
    case TopologyDescriptor::NBONA_UNUSED:
      tmp_desc[i] = heavy_bonds_plus_constraints;
      break;
    case TopologyDescriptor::NTHETA_UNUSED:
      tmp_desc[i] = heavy_angls_plus_constraints;
      break;
    case TopologyDescriptor::NPHIA_UNUSED:
      tmp_desc[i] = heavy_dihes_plus_constraints;
      break;
    case TopologyDescriptor::BOND_TYPE_COUNT:
      tmp_desc[i] = bond_parameter_count;
      break;
    case TopologyDescriptor::ANGLE_TYPE_COUNT:
      tmp_desc[i] = angl_parameter_count;
      break;
    case TopologyDescriptor::DIHEDRAL_TYPE_COUNT:
      tmp_desc[i] = dihe_parameter_count;
      break;
    case TopologyDescriptor::NATYP_UNUSED:
      tmp_desc[i] = unused_natyp;
      break;
    case TopologyDescriptor::NPHB_UNUSED:
      tmp_desc[i] = hbond_10_12_parameter_count;
      break;
    case TopologyDescriptor::PERTURBATION:
      break;
    case TopologyDescriptor::BOND_PERTURBATIONS:
      tmp_desc[i] = bond_perturbation_term_count;
      break;
    case TopologyDescriptor::ANGLE_PERTURBATIONS:
      tmp_desc[i] = angl_perturbation_term_count;
      break;
    case TopologyDescriptor::DIHEDRAL_PERTURBATIONS:
      tmp_desc[i] = dihe_perturbation_term_count;
      break;
    case TopologyDescriptor::BONDS_IN_PERTURBED_GROUP:
      tmp_desc[i] = bonds_in_perturbed_group;
      break;
    case TopologyDescriptor::ANGLES_IN_PERTURBED_GROUP:
      tmp_desc[i] = angls_in_perturbed_group;
      break;
    case TopologyDescriptor::DIHEDRALS_IN_PERTURBED_GROUP:
      tmp_desc[i] = dihes_in_perturbed_group;
      break;
    case TopologyDescriptor::BOX_TYPE_INDEX:
      switch (original.periodic_box_class) {
      case UnitCellType::NONE:
        tmp_desc[i] = 0;
      case UnitCellType::ORTHORHOMBIC:
        tmp_desc[i] = 1;
      case UnitCellType::TRICLINIC:
        tmp_desc[i] = 2;
      }
      break;
    case TopologyDescriptor::ATOM_COUNT_LARGEST_RESIDUE:
      tmp_desc[i] = largest_residue_size;
      break;
    case TopologyDescriptor::CAP:
      tmp_desc[i] = (use_solvent_cap_option == SolventCapSetting::ON);
      break;
    case TopologyDescriptor::EXTRA_POINT_COUNT:
      tmp_desc[i] = vst.vs_count;
      break;
    case TopologyDescriptor::PIMD_SLICE_COUNT:
      tmp_desc[i] = implicit_copy_count;
      break;
    case TopologyDescriptor::N_VALUES:
      break;
    }
  }
  
  // Populate the arrays of the new topology in the manner used to read data in from a file
  loadHybridArrays(tmp_desc, tmp_residue_limits, tmp_atom_struc_numbers, tmp_residue_numbers,
                   tmp_molecule_limits, tmp_atomic_numbers, tmp_molecule_membership,
                   tmp_mobile_atoms, tmp_molecule_contents, tmp_cmap_surface_dims,
                   tmp_cmap_surface_bounds, tmp_charge_type_indices, tmp_lennard_jones_indices,
                   tmp_inferred_14_i_atoms, tmp_inferred_14_l_atoms, tmp_inferred_14_param_idx,
                   tmp_neck_gb_indices, tmp_tree_joining_info, tmp_last_rotator_info,
                   tmp_charges, tmp_masses, tmp_urey_bradley_stiffnesses,
                   tmp_urey_bradley_equilibria, tmp_charmm_impr_stiffnesses,
                   tmp_charmm_impr_phase_angles, tmp_cmap_surfaces, tmp_bond_stiffnesses,
                   tmp_bond_equilibria, tmp_angl_stiffnesses, tmp_angl_equilibria,
                   tmp_dihe_amplitudes, tmp_dihe_periodicities, tmp_dihe_phase_angles,
                   tmp_charge_parameters, tmp_lja_coef, tmp_ljb_coef, tmp_ljc_coef,
                   tmp_lja_14_coef, tmp_ljb_14_coef, tmp_ljc_14_coef, tmp_atomic_pb_radii,
                   tmp_gb_screening_factors, tmp_gb_coef, tmp_solty_info, tmp_hbond_a_values,
                   tmp_hbond_b_values, tmp_hbond_cutoffs, tmp_atomic_polarizabilities,
                   tmp_atom_names, tmp_atom_types, tmp_residue_names, tmp_tree_symbols, cma,
                   cond_excl, bvt, mvt, attn_parm, vst, all_nb_excl, cnst_table);
}

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const AtomGraph &original, const std::vector<bool> &mask,
                     const ExceptionResponse policy, const double charge_rounding_tol,
                     const double charge_discretization,
                     const ApplyConstraints use_bond_constraints_in,
                     const ApplyConstraints use_settle_in, const HybridFormat format_in) :
    AtomGraph(original, enumerateMask(mask), policy, charge_rounding_tol, charge_discretization,
              use_bond_constraints_in, use_settle_in, format_in)
{}

//-------------------------------------------------------------------------------------------------
bool scanForExclusion(const std::vector<int> &list, const int llim, const int hlim,
                      const int partner) {
  bool result = false;
  for (int i = llim; i < hlim; i++) {
    result = (result || list[i] == partner);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void cullPriorExclusions(std::vector<int3> *proposals, const std::vector<int3> &priors) {
  int3* prop_ptr = proposals->data();
  const size_t nprop = proposals->size();
  const size_t nprrs = priors.size();
  size_t nunique = 0;
  for (size_t i = 0; i < nprop; i++) {
    bool found = false;
    const int atom_a = prop_ptr[i].x;
    const int atom_b = prop_ptr[i].y;
    for (size_t j = 0; j < nprrs; j++) {
      found = (found || (priors[j].x == atom_a && priors[j].y == atom_b) ||
                        (priors[j].x == atom_b && priors[j].y == atom_a));
    }
    if (! found) {
      prop_ptr[nunique] = prop_ptr[i];
      nunique++;
    }
  }
  proposals->resize(nunique);
}

} // namespace topology
} // namespace stormm
