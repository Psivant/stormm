#include <ctime>
#include <unistd.h>
#include "copyright.h"
#include "FileManagement/file_util.h"
#include "Math/series_ops.h"
#include "Math/vector_ops.h"
#include "atomgraph.h"
#include "lennard_jones_analysis.h"

namespace stormm {
namespace topology {

using diskutil::getBaseName;
using stmath::indexingArray;
using stmath::prefixSumInPlace;
using stmath::PrefixSumType;
  
//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const std::vector<AtomGraph*> &agv, const std::vector<int> &counts,
                     const MoleculeOrdering arrangement, const ExceptionResponse policy,
                     const HybridFormat format_in) :
    AtomGraph(format_in)
{
  snprintf(version_stamp, 16, "V8000.000");
  std::time_t raw_time = std::time(nullptr);
  std::tm* current_time = std::localtime(&raw_time);
  date = *current_time;
  const int nsys = agv.size();
  if (nsys == 0) {
    rtErr("At least one valid topology must be provided in order to compile multiple AtomGraph "
          "objects into a new topology.", "AtomGraph");
  }
  title = "Compilation of " + std::to_string(nsys) + " AtomGraph topologies in STORMM";
  if (nsys >= 2) {
    source = "\"Combination of ";
    for (int i = 0; i < nsys; i++) {
      if (counts[i] > 1) {
        source += std::to_string(counts[i]) + "x " + agv[i]->source;
      }
      else {
        source += agv[i]->source;
      }
      if (nsys == 2 && i == 0) {
        source += " and ";
      }
      else if (nsys > 2 && i < nsys - 2) {
        source += ", ";
      }
      else if (nsys > 2 && i < nsys - 1) {
        source += ", and ";
      }
    }
    source += "\"";
  }
  else {
    source = "Aggregate of " + std::to_string(counts[0]) + " x " + agv[0]->source;
  }
  force_fields = agv[0]->force_fields;
  for (int i = 0; i < nsys; i++) { 
    force_fields.insert(force_fields.end(), agv[i]->force_fields.begin(),
                        agv[i]->force_fields.end());
    atom_count += counts[i] * agv[i]->atom_count;
    residue_count += counts[i] * agv[i]->residue_count;
    molecule_count += counts[i] * agv[i]->molecule_count;
    largest_residue_size = std::max(largest_residue_size, agv[i]->largest_residue_size);
    largest_molecule_size = std::max(largest_molecule_size, agv[i]->largest_molecule_size);
    water_residue_size = std::max(water_residue_size, agv[i]->water_residue_size);
    last_solute_residue += counts[i] * agv[i]->last_solute_residue;
    last_solute_atom += counts[i] * agv[i]->last_solute_atom;
    first_solvent_molecule += counts[i] * agv[i]->first_solvent_molecule;
    last_atom_before_cap += counts[i] * agv[i]->last_atom_before_cap;
    unconstrained_dof += counts[i] * agv[i]->unconstrained_dof;
    constrained_dof += counts[i] * agv[i]->constrained_dof;
    bond_term_with_hydrogen += counts[i] * agv[i]->bond_term_with_hydrogen;
    bond_term_without_hydrogen += counts[i] * agv[i]->bond_term_without_hydrogen;
    angl_term_with_hydrogen += counts[i] * agv[i]->angl_term_with_hydrogen;
    angl_term_without_hydrogen += counts[i] * agv[i]->angl_term_without_hydrogen;
    dihe_term_with_hydrogen += counts[i] * agv[i]->dihe_term_with_hydrogen;
    dihe_term_without_hydrogen += counts[i] * agv[i]->dihe_term_without_hydrogen;
  }

  // Form consensus tables of bond, angle, dihedral, Urey-Bradley, CHARMM improper, and CMAP terms.
  // Try to preserve the order of parameter indices in each topology as much as possible.
  std::vector<ValenceKit<double>> vk_v;
  std::vector<NonbondedKit<double>> nbk_v;
  std::vector<std::vector<std::vector<char4>>> atyp_list_v;
  std::vector<VirtualSiteKit<double>> vsk_v;
  vk_v.reserve(nsys);
  nbk_v.reserve(nsys);
  atyp_list_v.reserve(nsys);
  vsk_v.reserve(nsys);
  for (int i = 0; i < nsys; i++) {
    vk_v.push_back(agv[i]->getDoublePrecisionValenceKit());
    nbk_v.push_back(agv[i]->getDoublePrecisionNonbondedKit());
    atyp_list_v.push_back(agv[i]->getAtomTypeNameTable());
    vsk_v.push_back(agv[i]->getDoublePrecisionVirtualSiteKit());
  }

  // Treat each dihedral parameter set and its associated 1:4 scaling factor as a coherent set of
  // parameters.  This implies that each unique dihedral parameter set will have one and only one
  // set of 1:4 attenuation factors, but that is the convention in Amber topologies and if there
  // are situations where a single set of dihedral parameters is associated with more than one
  // set of 1:4 attenuations then that can be worked out by cloning the dihedral parameter sets.
  std::vector<std::vector<double>> tmp_elec14_attns(nsys), tmp_vdw14_attns(nsys);  
  for (int i = 0; i < nsys; i++) {
    tmp_elec14_attns[i].resize(vk_v[i].ndihe_param, 0.0);
    tmp_vdw14_attns[i].resize(vk_v[i].ndihe_param, 0.0);
    for (int j = 0; j < vk_v[i].ndihe; j++) {
      const int attn14_param_idx = vk_v[i].dihe14_param_idx[j];
      const int dihe_param_idx = vk_v[i].dihe_param_idx[j];
      tmp_elec14_attns[i][dihe_param_idx] = std::max(tmp_elec14_attns[i][dihe_param_idx],
                                                     vk_v[i].attn14_elec[attn14_param_idx]);
      tmp_vdw14_attns[i][dihe_param_idx] = std::max(tmp_vdw14_attns[i][dihe_param_idx],
                                                    vk_v[i].attn14_vdw[attn14_param_idx]);
    }
  }
  ParameterUnion<double> uni_chrgs(nbk_v[0].q_parameter, nbk_v[0].n_q_types);
  ParameterUnion<double> uni_bonds(vk_v[0].bond_keq, vk_v[0].bond_leq, vk_v[0].nbond_param);
  ParameterUnion<double> uni_angls(vk_v[0].angl_keq, vk_v[0].angl_theta, vk_v[0].nangl_param);

  // Some topologies may not have any 1:4 interactions.  In this case, the entries of
  // tmp_elec14_attns and tmp_vdw141_attns will have null pointers.
  std::vector<double> attn_ph(1, 0.0);
  double* tmp_elec14_attn_ptr;
  double* tmp_vdw14_attn_ptr;
  if (vk_v[0].ndihe_param == 0 || vk_v[0].nattn14_param == 0) {
    tmp_elec14_attn_ptr = attn_ph.data();
    tmp_vdw14_attn_ptr = attn_ph.data();
  }
  else {
    tmp_elec14_attn_ptr = tmp_elec14_attns[0].data();
    tmp_vdw14_attn_ptr = tmp_vdw14_attns[0].data();
  }
  ParameterUnion<double> uni_dihes(vk_v[0].dihe_amp, vk_v[0].dihe_freq, vk_v[0].dihe_phi,
                                   tmp_elec14_attn_ptr, tmp_vdw14_attn_ptr, vk_v[0].ndihe_param);
  ParameterUnion<double> uni_hbonds(agv[0]->hbond_a_values.data(), agv[0]->hbond_b_values.data(),
                                    agv[0]->hbond_10_12_parameter_count);
  ParameterUnion<double> uni_attns(vk_v[0].attn14_elec, vk_v[0].attn14_vdw, vk_v[0].nattn14_param);
  ParameterUnion<double> uni_ubrds(vk_v[0].ubrd_keq, vk_v[0].ubrd_leq, vk_v[0].nubrd_param);
  ParameterUnion<double> uni_cimps(vk_v[0].cimp_keq, vk_v[0].cimp_phi, vk_v[0].ncimp_param);
  CmapSurfaceUnion uni_cmaps(vk_v[0].cmap_surf, vk_v[0].cmap_dim, vk_v[0].ncmap_surf);
  LennardJonesAnalysis uni_ljtab(nbk_v[0], atyp_list_v[0]);
  ParameterUnion<double> uni_vsites(vsk_v[0].dim1, vsk_v[0].dim2, vsk_v[0].dim3,
                                    vsk_v[0].nframe_set);
  for (int i = 1; i < nsys; i++) {
    uni_chrgs.addSet(nbk_v[i].q_parameter, nbk_v[i].n_q_types);
    uni_bonds.addSet(vk_v[i].bond_keq, vk_v[i].bond_leq, vk_v[i].nbond_param);
    uni_angls.addSet(vk_v[i].angl_keq, vk_v[i].angl_theta, vk_v[i].nangl_param);
    if (vk_v[i].ndihe_param == 0 || vk_v[i].nattn14_param == 0) {
      tmp_elec14_attn_ptr = attn_ph.data();
      tmp_vdw14_attn_ptr = attn_ph.data();
    }
    else {
      tmp_elec14_attn_ptr = tmp_elec14_attns[i].data();
      tmp_vdw14_attn_ptr = tmp_vdw14_attns[i].data();
    }
    uni_dihes.addSet(vk_v[i].dihe_amp, vk_v[i].dihe_freq, vk_v[i].dihe_phi,
                     tmp_elec14_attn_ptr, tmp_vdw14_attn_ptr, vk_v[i].ndihe_param);
    uni_attns.addSet(vk_v[i].attn14_elec, vk_v[i].attn14_vdw, vk_v[i].nattn14_param);
    uni_hbonds.addSet(agv[i]->hbond_a_values.data(), agv[i]->hbond_b_values.data(),
                      agv[i]->hbond_10_12_parameter_count);
    uni_ubrds.addSet(vk_v[i].ubrd_keq, vk_v[i].ubrd_leq, vk_v[i].nubrd_param);
    uni_cimps.addSet(vk_v[i].cimp_keq, vk_v[i].cimp_phi, vk_v[i].ncimp_param);
    uni_cmaps.addSet(vk_v[i].cmap_surf, vk_v[i].cmap_dim, vk_v[i].ncmap_surf);
    uni_ljtab.addSet(nbk_v[i], atyp_list_v[i]);
    uni_vsites.addSet(vsk_v[i].dim1, vsk_v[i].dim2, vsk_v[i].dim3, vsk_v[i].nframe_set);
  }
  
  // The arrays that enter into the topology must be constructed from scratch in order to feed them
  // to the loadHybridArrays() function used by other constructors.
  const int n_desc = static_cast<int>(TopologyDescriptor::N_VALUES);
  std::vector<int> tmp_desc(n_desc);
  for (int i = 0; i < n_desc; i++) {
    switch (static_cast<TopologyDescriptor>(i)) {
    case TopologyDescriptor::ATOM_COUNT:
      tmp_desc[i] = atom_count;
      break;
    case TopologyDescriptor::ATOM_TYPE_COUNT:

      // Draw upon the consensus Lennard-Jones tables for the two topologies.
      tmp_desc[i] = uni_ljtab.getLJTypeCount();      
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
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        tmp_desc[i] += counts[j] * agv[j]->unused_nhparm;
      }
      break;
    case TopologyDescriptor::ADDLES_CREATED:
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        tmp_desc[i] += counts[j] * agv[j]->unused_nparm;
      }
      break;
    case TopologyDescriptor::TOTAL_EXCLUDED_ATOMS:
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        tmp_desc[i] += counts[j] * agv[j]->total_exclusions;
      }
      break;
    case TopologyDescriptor::RESIDUE_COUNT:
      tmp_desc[i] = residue_count;
      break;
    case TopologyDescriptor::NBONA_UNUSED:
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        tmp_desc[i] += counts[j] * agv[j]->heavy_bonds_plus_constraints;
      }
      break;
    case TopologyDescriptor::NTHETA_UNUSED:
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        tmp_desc[i] += counts[j] * agv[j]->heavy_angls_plus_constraints;
      }
      break;
    case TopologyDescriptor::NPHIA_UNUSED:
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        tmp_desc[i] += counts[j] * agv[j]->heavy_dihes_plus_constraints;
      }
      break;
    case TopologyDescriptor::BOND_TYPE_COUNT:

      // Draw upon the consensus bond parameter tables.
      tmp_desc[i] = uni_bonds.getUniqueParameterCount();
      break;
    case TopologyDescriptor::ANGLE_TYPE_COUNT:

      // Draw upon the consensus angle parameter tables.
      tmp_desc[i] = uni_angls.getUniqueParameterCount();
      break;
    case TopologyDescriptor::DIHEDRAL_TYPE_COUNT:

      // Draw upon the consensus dihedral parameter tables.
      tmp_desc[i] = uni_dihes.getUniqueParameterCount();
      break;
    case TopologyDescriptor::NATYP_UNUSED:

      // This information is unused, and it is not clear where to find the atom types involved in
      // hydrogen bonding.  Assume that the types in each of the various systems are independent.
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        tmp_desc[i] += counts[j] * agv[j]->unused_natyp;
      }
      break;
    case TopologyDescriptor::NPHB_UNUSED:

      // Draw upon the consensus hydrogen bonding parameter tables
      tmp_desc[i] = uni_hbonds.getUniqueParameterCount();
      break;
    case TopologyDescriptor::PERTURBATION:
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        const int active = (agv[j]->use_perturbation_info == PerturbationSetting::ON);
        tmp_desc[i] = std::max(tmp_desc[i], active);
      }
      break;
    case TopologyDescriptor::BOND_PERTURBATIONS:
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        tmp_desc[i] += counts[j] * agv[j]->bond_perturbation_term_count;
      }
      break;
    case TopologyDescriptor::ANGLE_PERTURBATIONS:
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        tmp_desc[i] += counts[j] * agv[j]->angl_perturbation_term_count;
      }
      break;
    case TopologyDescriptor::DIHEDRAL_PERTURBATIONS:
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        tmp_desc[i] += counts[j] * agv[j]->dihe_perturbation_term_count;
      }
      break;
    case TopologyDescriptor::BONDS_IN_PERTURBED_GROUP:
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        tmp_desc[i] += counts[j] * agv[j]->bonds_in_perturbed_group;
      }
      break;
    case TopologyDescriptor::ANGLES_IN_PERTURBED_GROUP:
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        tmp_desc[i] += counts[j] * agv[j]->angls_in_perturbed_group;
      }
      break;
    case TopologyDescriptor::DIHEDRALS_IN_PERTURBED_GROUP:
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        tmp_desc[i] += counts[j] * agv[j]->dihes_in_perturbed_group;
      }
      break;
    case TopologyDescriptor::BOX_TYPE_INDEX:

      // Type promotion: while it could be infeasible to combine certain systems' coordinates,
      // assume that the system with the "higher order" unit cell dictates the unit cell of the
      // combined system.
      tmp_desc[i] = 0;
      int cls_idx;
      for (int j = 0; j < nsys; j++) {
        switch (agv[j]->periodic_box_class) {
        case UnitCellType::NONE:
          cls_idx = 0;
          break;
        case UnitCellType::ORTHORHOMBIC:
          cls_idx = 1;
          break;
        case UnitCellType::TRICLINIC:
          cls_idx = 2;
          break;
        }
        tmp_desc[i] = std::max(tmp_desc[i], cls_idx);
      }
      break;
    case TopologyDescriptor::ATOM_COUNT_LARGEST_RESIDUE:
      tmp_desc[i] = largest_residue_size;
      break;
    case TopologyDescriptor::CAP:
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        const int active = (agv[j]->use_solvent_cap_option == SolventCapSetting::ON);
        tmp_desc[i] = std::max(tmp_desc[i], active);
      }
      break;
    case TopologyDescriptor::EXTRA_POINT_COUNT:
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        tmp_desc[i] = counts[j] * agv[j]->virtual_site_count;
      }
      break;
    case TopologyDescriptor::PIMD_SLICE_COUNT:
      tmp_desc[i] = 0;
      for (int j = 0; j < nsys; j++) {
        tmp_desc[i] = std::max(tmp_desc[i], agv[j]->implicit_copy_count);
      }
      break;
    case TopologyDescriptor::N_VALUES:
      break;
    }
  }
  switch (tmp_desc[static_cast<size_t>(TopologyDescriptor::BOX_TYPE_INDEX)]) {
  case 0:
    periodic_box_class = UnitCellType::NONE;
    break;
  case 1:
    periodic_box_class = UnitCellType::ORTHORHOMBIC;
    break;
  case 2:
    periodic_box_class = UnitCellType::TRICLINIC;
    break;
  }
  
  // The details above are sums over the entire topology that is to be created.  The ordering of
  // atoms and molecules is still up in the air, and to simply tack one topology onto the end of
  // another may deal unexpected twists to a lot of scripts.  Loop over each topology and find its
  // biomolecular, general organic, ions or other small molecules, and finally water molecules.
  // These will be ordered in the new topology in the above order, arranging similar components by
  // the order in which their topologies appear in the list.
  std::vector<int2> molecule_system_instance_origins;
  const std::vector<int2> molecule_order = findMoleculeOrder(agv, counts,
                                                             &molecule_system_instance_origins,
                                                             arrangement);

  // As in other applications, it will be important to have maps of where each atom of the new
  // topology came from, and reverse lookups to understand where in the new topology each atom from
  // one of the original topologies will go.
  std::vector<int2> atom_origins(atom_count);
  std::vector<std::vector<int>> atom_destinations(nsys);
  std::vector<std::vector<int>> destcon(nsys);
  for (int i = 0; i < nsys; i++) {
    const int ni_atom = agv[i]->atom_count;
    atom_destinations[i] = std::vector<int>(counts[i] * ni_atom);
    destcon[i] = std::vector<int>(counts[i] * ni_atom);
    for (int j = 0; j < ni_atom; j++) {
      destcon[i][j] = j * counts[i];
    }
  }

  // Find an order of the atoms and residues in the original topologies upon which to base the
  // properties and residue limits of the new topology.  Lay down atomic, residue, and molecular
  // descriptors.  Structural atom numbers will proceed in the original order of each topology,
  // again incremented by the total number of atoms already laid down.
  std::vector<ChemicalDetailsKit> cdk_v;
  cdk_v.reserve(nsys);
  for (int i = 0; i < nsys; i++) {
    cdk_v.push_back(agv[i]->getChemicalDetailsKit());
  }
  std::vector<std::vector<int>> all_res_idx(nsys);
  for (int i = 0; i < nsys; i++) {
    all_res_idx[i] = agv[i]->getResidueIndex();
  }
  std::vector<int> tmp_residue_limits(residue_count + 1), tmp_molecule_limits(molecule_count + 1);
  std::vector<int> tmp_atom_struc_numbers(atom_count), tmp_residue_numbers(atom_count);
  std::vector<int> tmp_molecule_membership(atom_count), tmp_molecule_contents(atom_count);
  std::vector<int> tmp_atomic_numbers(atom_count);
  std::vector<double> tmp_atomic_charges(atom_count), tmp_atomic_masses(atom_count);
  std::vector<double> tmp_inv_atomic_masses(atom_count), tmp_atomic_polarizabilities(atom_count);
  std::vector<char4> tmp_atom_names(atom_count), tmp_atom_types(atom_count);
  std::vector<char4> tmp_residue_names(residue_count);
  int atmcon = 0;
  int rsdcon = 0;
  int base_atom_idx, base_resi_idx;
  for (int i = 0; i < molecule_count; i++) {
    const int sys_orig = molecule_order[i].x;
    const int mol_orig = molecule_order[i].y;
    const int mol_llim = cdk_v[sys_orig].mol_limits[mol_orig];
    const int mol_hlim = cdk_v[sys_orig].mol_limits[mol_orig + 1];
    int current_topl_res, mol_base_res;
    base_atom_idx = atmcon;
    base_resi_idx = rsdcon;
    tmp_molecule_limits[i] = atmcon;
    const double* polarizability_ptr = agv[sys_orig]->atomic_polarizabilities.data();
    for (int j = mol_llim; j < mol_hlim; j++) {
      const int sys_atom_idx = cdk_v[sys_orig].mol_contents[j];

      // Initialize and track the residues within each molecule, incrementing a counter for the
      // combined topology.  
      if (j == mol_llim) {
        current_topl_res = all_res_idx[sys_orig][sys_atom_idx];
        mol_base_res = current_topl_res;
        tmp_residue_limits[rsdcon] = atmcon;
        tmp_residue_names[rsdcon] = agv[sys_orig]->residue_names.readHost(current_topl_res);
      }
      if (current_topl_res != all_res_idx[sys_orig][sys_atom_idx]) {
        rsdcon++;
        tmp_residue_limits[rsdcon] = atmcon;
        current_topl_res = all_res_idx[sys_orig][sys_atom_idx];
        tmp_residue_names[rsdcon] = agv[sys_orig]->residue_names.readHost(current_topl_res);
      }
      tmp_atom_struc_numbers[atmcon] = cdk_v[sys_orig].atom_numbers[sys_atom_idx] + j - mol_llim +
                                       base_atom_idx;
      tmp_residue_numbers[atmcon] = cdk_v[sys_orig].res_numbers[sys_atom_idx] - mol_base_res +
                                    base_resi_idx;
      tmp_molecule_membership[atmcon] = i;
      tmp_molecule_contents[atmcon] = atmcon;
      tmp_atomic_charges[atmcon] = nbk_v[sys_orig].charge[sys_atom_idx];
      tmp_atomic_masses[atmcon] = cdk_v[sys_orig].masses[sys_atom_idx];
      tmp_inv_atomic_masses[atmcon] = cdk_v[sys_orig].inv_masses[sys_atom_idx];
      tmp_atomic_numbers[atmcon] = cdk_v[sys_orig].z_numbers[sys_atom_idx];
      tmp_atom_names[atmcon] = cdk_v[sys_orig].atom_names[sys_atom_idx];
      tmp_atom_types[atmcon] = cdk_v[sys_orig].atom_types[sys_atom_idx];
      tmp_atomic_polarizabilities[atmcon] = polarizability_ptr[sys_atom_idx];

      // Mark the origins and destinations of atoms in the combined and original topologies
      atom_origins[atmcon] = { sys_orig, sys_atom_idx };
      const int destination_counter = destcon[sys_orig][sys_atom_idx];
      atom_destinations[sys_orig][destination_counter] = atmcon;
      destcon[sys_orig][sys_atom_idx] = destination_counter + 1;
      atmcon++;
    }
    rsdcon++;
  }
  tmp_molecule_limits[molecule_count] = atmcon;
  tmp_residue_limits[residue_count] = atmcon;
  
  // Map the CHARMM valence terms for the combined topology.
  for (int i = 0; i < nsys; i++) {
    urey_bradley_term_count += counts[i] * agv[i]->urey_bradley_term_count;
    charmm_impr_term_count += counts[i] * agv[i]->charmm_impr_term_count;
    cmap_term_count += counts[i] * agv[i]->cmap_term_count;
  }
  urey_bradley_parameter_count = uni_ubrds.getUniqueParameterCount();
  charmm_impr_parameter_count = uni_cimps.getUniqueParameterCount();
  cmap_surface_count = uni_cmaps.getUniqueSurfaceCount();
  std::vector<int> tmp_ubrd_i_atoms(urey_bradley_term_count);
  std::vector<int> tmp_ubrd_k_atoms(urey_bradley_term_count);
  std::vector<int> tmp_ubrd_param_idx(urey_bradley_term_count);
  std::vector<int> tmp_cimp_i_atoms(charmm_impr_term_count);
  std::vector<int> tmp_cimp_j_atoms(charmm_impr_term_count);
  std::vector<int> tmp_cimp_k_atoms(charmm_impr_term_count);
  std::vector<int> tmp_cimp_l_atoms(charmm_impr_term_count);
  std::vector<int> tmp_cimp_param_idx(charmm_impr_term_count);
  std::vector<int> tmp_cmap_i_atoms(cmap_term_count);
  std::vector<int> tmp_cmap_j_atoms(cmap_term_count);
  std::vector<int> tmp_cmap_k_atoms(cmap_term_count);
  std::vector<int> tmp_cmap_l_atoms(cmap_term_count);
  std::vector<int> tmp_cmap_m_atoms(cmap_term_count);
  std::vector<int> tmp_cmap_param_idx(cmap_term_count);
  int ubrd_con = 0;
  int cimp_con = 0;
  int cmap_con = 0;
  for (int i = 0; i < molecule_count; i++) {
    const int sys_instance = molecule_system_instance_origins[i].y;
    for (int j = tmp_molecule_limits[i]; j < tmp_molecule_limits[i + 1]; j++) {
      const int j_atom = tmp_molecule_contents[j];
      
      // Obtain the origin of the atom in the new topology
      const int sys_orig = atom_origins[j_atom].x;
      const int idx_orig = atom_origins[j_atom].y;

      // Loop over all associated Urey-Bradley terms
      const int ubrd_llim = vk_v[sys_orig].ubrd_asgn_bounds[idx_orig];
      const int ubrd_hlim = vk_v[sys_orig].ubrd_asgn_bounds[idx_orig + 1];
      for (int k = ubrd_llim; k < ubrd_hlim; k++) {
        tmp_ubrd_i_atoms[ubrd_con] = j_atom;

        // The destinations of each atom from one of the input topologies are known, but there may
        // be several such destinations.  Choose the one that is within the same molecule.  A
        // specific molecule must come from one and only one instance of a particular topology.
        // That instance was recorded as the molecules for the new topology were being laid out.
        const int dest_ref_idx = (counts[sys_orig] * vk_v[sys_orig].ubrd_asgn_atoms[k]) +
                                 sys_instance;
        tmp_ubrd_k_atoms[ubrd_con] = atom_destinations[sys_orig][dest_ref_idx];
        tmp_ubrd_param_idx[ubrd_con] =
          uni_ubrds.getCorrespondence(sys_orig, vk_v[sys_orig].ubrd_asgn_index[k]);
        ubrd_con++;
      }

      // Loop over all associated CHARMM improper terms
      const int cimp_llim = vk_v[sys_orig].cimp_asgn_bounds[idx_orig];
      const int cimp_hlim = vk_v[sys_orig].cimp_asgn_bounds[idx_orig + 1];
      for (int k = cimp_llim; k < cimp_hlim; k++) {
        tmp_cimp_j_atoms[cimp_con] = j_atom;
        const int idest_ref_idx = (counts[sys_orig] *
                                   vk_v[sys_orig].cimp_asgn_atoms[3 * k]) + sys_instance;
        const int kdest_ref_idx = (counts[sys_orig] *
                                   vk_v[sys_orig].cimp_asgn_atoms[(3 * k) + 1]) + sys_instance;
        const int ldest_ref_idx = (counts[sys_orig] *
                                   vk_v[sys_orig].cimp_asgn_atoms[(3 * k) + 2]) + sys_instance;
        tmp_cimp_i_atoms[cimp_con] = atom_destinations[sys_orig][idest_ref_idx];
        tmp_cimp_k_atoms[cimp_con] = atom_destinations[sys_orig][kdest_ref_idx];
        tmp_cimp_l_atoms[cimp_con] = atom_destinations[sys_orig][ldest_ref_idx];
        tmp_cimp_param_idx[cimp_con] =
          uni_cimps.getCorrespondence(sys_orig, vk_v[sys_orig].cimp_asgn_index[k]);
        cimp_con++;
      }

      // Loop over all associated CMAP terms
      const int cmap_llim = vk_v[sys_orig].cmap_asgn_bounds[idx_orig];
      const int cmap_hlim = vk_v[sys_orig].cmap_asgn_bounds[idx_orig + 1];
      for (int k = cmap_llim; k < cmap_hlim; k++) {
        tmp_cmap_k_atoms[cmap_con] = j_atom;
        const int idest_ref_idx = (counts[sys_orig] *
                                   vk_v[sys_orig].cmap_asgn_atoms[4 * k]) + sys_instance;
        const int jdest_ref_idx = (counts[sys_orig] *
                                   vk_v[sys_orig].cmap_asgn_atoms[(4 * k) + 1]) + sys_instance;
        const int ldest_ref_idx = (counts[sys_orig] *
                                   vk_v[sys_orig].cmap_asgn_atoms[(4 * k) + 2]) + sys_instance;
        const int mdest_ref_idx = (counts[sys_orig] *
                                   vk_v[sys_orig].cmap_asgn_atoms[(4 * k) + 3]) + sys_instance;
        tmp_cmap_i_atoms[cmap_con] = atom_destinations[sys_orig][idest_ref_idx];
        tmp_cmap_j_atoms[cmap_con] = atom_destinations[sys_orig][jdest_ref_idx];
        tmp_cmap_l_atoms[cmap_con] = atom_destinations[sys_orig][ldest_ref_idx];
        tmp_cmap_m_atoms[cmap_con] = atom_destinations[sys_orig][mdest_ref_idx];
        tmp_cmap_param_idx[cmap_con] =
          uni_cmaps.getCorrespondence(sys_orig, vk_v[sys_orig].cmap_asgn_index[k]);
        cmap_con++;
      }
    }
  }
  CharmmValenceTable cvt(atom_count, urey_bradley_term_count, charmm_impr_term_count,
                         cmap_term_count, tmp_ubrd_i_atoms, tmp_ubrd_k_atoms, tmp_ubrd_param_idx,
                         tmp_cimp_i_atoms, tmp_cimp_j_atoms, tmp_cimp_k_atoms, tmp_cimp_l_atoms,
                         tmp_cimp_param_idx, tmp_cmap_i_atoms, tmp_cmap_j_atoms, tmp_cmap_k_atoms,
                         tmp_cmap_l_atoms, tmp_cmap_m_atoms, tmp_cmap_param_idx);

  // Construct the CMAP surface derivatives.
  cmap_term_count = uni_cmaps.getUniqueSurfaceCount();
  std::vector<int> tmp_cmap_surface_bounds(cmap_term_count + 1, 0);
  const std::vector<int> tmp_cmap_surface_dims = uni_cmaps.getSurfaceDimensions();
  for (int i = 0; i < cmap_term_count; i++) {
    tmp_cmap_surface_bounds[i] = tmp_cmap_surface_dims[i] * tmp_cmap_surface_dims[i];
  }
  prefixSumInPlace<int>(&tmp_cmap_surface_bounds, PrefixSumType::EXCLUSIVE);
  const CmapAccessories cma = computeCmapDerivatives(cmap_term_count,
                                                     uni_cmaps.getSurfaceDimensions(),
                                                     tmp_cmap_surface_bounds,
                                                     uni_cmaps.getAllSurfaces());
  
  // Map the basic force field valence terms for the combined topology.  This loop follows the same
  // structure as above, but writing a separate block of code helps to organize the different
  // terms.
  for (int i = 0; i < nsys; i++) {
    bond_term_count += counts[i] * agv[i]->bond_term_count;
    angl_term_count += counts[i] * agv[i]->angl_term_count;
    dihe_term_count += counts[i] * agv[i]->dihe_term_count;
  }
  bond_parameter_count = uni_bonds.getUniqueParameterCount();
  angl_parameter_count = uni_angls.getUniqueParameterCount();
  dihe_parameter_count = uni_dihes.getUniqueParameterCount();
  std::vector<int> tmp_bond_i_atoms(bond_term_count), tmp_bond_j_atoms(bond_term_count);
  std::vector<int> tmp_bond_param_idx(bond_term_count);
  std::vector<int> tmp_angl_i_atoms(angl_term_count), tmp_angl_j_atoms(angl_term_count);
  std::vector<int> tmp_angl_k_atoms(angl_term_count), tmp_angl_param_idx(angl_term_count);
  std::vector<int> tmp_dihe_i_atoms(dihe_term_count), tmp_dihe_j_atoms(dihe_term_count);
  std::vector<int> tmp_dihe_k_atoms(dihe_term_count), tmp_dihe_l_atoms(dihe_term_count);
  std::vector<int> tmp_dihe_param_idx(dihe_term_count);
  std::vector<char4> tmp_bond_mods(bond_term_count),  tmp_angl_mods(angl_term_count);
  std::vector<char4> tmp_dihe_mods(dihe_term_count);
  int bond_con = 0;
  int angl_con = 0;
  int dihe_con = 0;
  for (int i = 0; i < molecule_count; i++) {
    const int sys_instance = molecule_system_instance_origins[i].y;
    for (int j = tmp_molecule_limits[i]; j < tmp_molecule_limits[i + 1]; j++) {
      const int j_atom = tmp_molecule_contents[j];
      
      // Obtain the origin of the atom in the new topology
      const int sys_orig = atom_origins[j_atom].x;
      const int idx_orig = atom_origins[j_atom].y;

      // Loop over all associated bond terms
      const int bond_llim = vk_v[sys_orig].bond_asgn_bounds[idx_orig];
      const int bond_hlim = vk_v[sys_orig].bond_asgn_bounds[idx_orig + 1];
      for (int k = bond_llim; k < bond_hlim; k++) {
        tmp_bond_i_atoms[bond_con] = j_atom;
        const int bond_term_idx = vk_v[sys_orig].bond_asgn_terms[k];
        const int dest_ref_idx = (counts[sys_orig] * vk_v[sys_orig].bond_asgn_atoms[k]) +
                                 sys_instance;
        tmp_bond_j_atoms[bond_con] = atom_destinations[sys_orig][dest_ref_idx];
        tmp_bond_mods[bond_con] = vk_v[sys_orig].bond_modifiers[bond_term_idx];
        tmp_bond_param_idx[bond_con] =
          uni_bonds.getCorrespondence(sys_orig, vk_v[sys_orig].bond_asgn_index[k]);
        bond_con++;
      }

      // Loop over all associated angle terms
      const int angl_llim = vk_v[sys_orig].angl_asgn_bounds[idx_orig];
      const int angl_hlim = vk_v[sys_orig].angl_asgn_bounds[idx_orig + 1];
      for (int k = angl_llim; k < angl_hlim; k++) {
        tmp_angl_j_atoms[angl_con] = j_atom;
        const int angl_term_idx = vk_v[sys_orig].angl_asgn_terms[k];
        const int idest_ref_idx = (counts[sys_orig] *
                                   vk_v[sys_orig].angl_asgn_atoms[2 * k]) + sys_instance;
        const int kdest_ref_idx = (counts[sys_orig] *
                                   vk_v[sys_orig].angl_asgn_atoms[(2 * k) + 1]) + sys_instance;
        tmp_angl_i_atoms[angl_con] = atom_destinations[sys_orig][idest_ref_idx];
        tmp_angl_k_atoms[angl_con] = atom_destinations[sys_orig][kdest_ref_idx];
        tmp_angl_mods[angl_con] = vk_v[sys_orig].angl_modifiers[angl_term_idx];
        tmp_angl_param_idx[angl_con] =
          uni_angls.getCorrespondence(sys_orig, vk_v[sys_orig].angl_asgn_index[k]);
        angl_con++;
      }

      // Loop over all associated dihedral terms
      const int dihe_llim = vk_v[sys_orig].dihe_asgn_bounds[idx_orig];
      const int dihe_hlim = vk_v[sys_orig].dihe_asgn_bounds[idx_orig + 1];
      for (int k = dihe_llim; k < dihe_hlim; k++) {
        tmp_dihe_k_atoms[dihe_con] = j_atom;
        const int dihe_term_idx = vk_v[sys_orig].dihe_asgn_terms[k];
        const int idest_ref_idx = (counts[sys_orig] *
                                   vk_v[sys_orig].dihe_asgn_atoms[3 * k]) + sys_instance;
        const int jdest_ref_idx = (counts[sys_orig] *
                                   vk_v[sys_orig].dihe_asgn_atoms[(3 * k) + 1]) + sys_instance;
        const int ldest_ref_idx = (counts[sys_orig] *
                                   vk_v[sys_orig].dihe_asgn_atoms[(3 * k) + 2]) + sys_instance;
        tmp_dihe_i_atoms[dihe_con] = atom_destinations[sys_orig][idest_ref_idx];
        tmp_dihe_j_atoms[dihe_con] = atom_destinations[sys_orig][jdest_ref_idx];
        tmp_dihe_l_atoms[dihe_con] = atom_destinations[sys_orig][ldest_ref_idx];
        tmp_dihe_mods[dihe_con] = vk_v[sys_orig].dihe_modifiers[dihe_term_idx];
        tmp_dihe_param_idx[dihe_con] =
          uni_dihes.getCorrespondence(sys_orig, vk_v[sys_orig].dihe_asgn_index[k]);
        dihe_con++;
      }
    }
  }
  BasicValenceTable bvt(atom_count, bond_term_count, angl_term_count, dihe_term_count,
                        tmp_bond_i_atoms, tmp_bond_j_atoms, tmp_bond_param_idx, tmp_angl_i_atoms,
                        tmp_angl_j_atoms, tmp_angl_k_atoms, tmp_angl_param_idx, tmp_dihe_i_atoms,
                        tmp_dihe_j_atoms, tmp_dihe_k_atoms, tmp_dihe_l_atoms, tmp_dihe_param_idx,
                        tmp_bond_mods, tmp_angl_mods, tmp_dihe_mods);
  bvt.checkBondAngleRelevance(agv[0]->use_shake, agv[0]->use_settle, tmp_atomic_numbers);
  
  // Map the virtual sites of the new topology.  Unlike valence parameter terms, topologies do not
  // keep their own maps of "assigned" virtual sites for each atom, although the parent atoms can
  // be considered to "own" each virtual site and then reference the associated atoms.  Therefore,
  // make temporary maps of each of the input topologies, noting whether any particular atom is the
  // parent atoms of one or more virtual sites.
  std::vector<std::vector<int>> vs_map(nsys);
  std::vector<std::vector<int>> vs_map_bounds(nsys);
  for (int i = 0; i < nsys; i++) {
    virtual_site_count += counts[i] * agv[i]->virtual_site_count;
    vs_map_bounds[i].resize(agv[i]->atom_count + 1, 0);
    vs_map[i].resize(vsk_v[i].nsite);

    // While it may seem appropriate to use indexingArray() here, the fact that the bounds go over
    // parent atoms while the map needs to list virtual site particle indices is the problem.
    for (int j = 0; j < vsk_v[i].nsite; j++) {
      vs_map_bounds[i][vsk_v[i].frame1_idx[j]] +=1;
    }
    prefixSumInPlace<int>(&vs_map_bounds[i], PrefixSumType::EXCLUSIVE);
    std::vector<int> vs_map_counters = vs_map_bounds[i];
    for (int j = 0; j < vsk_v[i].nsite; j++) {
      const int parent_atom = vsk_v[i].frame1_idx[j];
      const int map_idx = vs_map_counters[parent_atom];
      vs_map[i][map_idx] = vsk_v[i].vs_atoms[j];
      vs_map_counters[parent_atom] = map_idx + 1;
    }
  }
  virtual_site_parameter_set_count = uni_vsites.getUniqueParameterCount();
  std::vector<int> tmp_vs_atoms(virtual_site_count);
  std::vector<int> tmp_vs_frame_types(virtual_site_count);
  std::vector<int> tmp_vs_frame1_idx(virtual_site_count);
  std::vector<int> tmp_vs_frame2_idx(virtual_site_count);
  std::vector<int> tmp_vs_frame3_idx(virtual_site_count);
  std::vector<int> tmp_vs_frame4_idx(virtual_site_count);
  std::vector<int> tmp_vs_parameter_indices(virtual_site_count);
  int vsite_con = 0;
  for (int i = 0; i < molecule_count; i++) {
    const int sys_instance = molecule_system_instance_origins[i].y;
    for (int j = tmp_molecule_limits[i]; j < tmp_molecule_limits[i + 1]; j++) {
      const int j_atom = tmp_molecule_contents[j];
      
      // Obtain the origin of the atom in the new topology
      const int sys_orig = atom_origins[j_atom].x;
      const int idx_orig = atom_origins[j_atom].y;
      const int vs_llim = vs_map_bounds[sys_orig][idx_orig];
      const int vs_hlim = vs_map_bounds[sys_orig][idx_orig + 1];
      for (int k = vs_llim; k < vs_hlim; k++) {

        // Construct the new virtual site frame in the topology.  The actual particles will be
        // present in the topology and their order will not be affected, but this will add the
        // details of the frame.
        const int vs_dest_ref_idx = (counts[sys_orig] * vsk_v[sys_orig].vs_atoms[k]) +
                                    sys_instance;
        const int f1_dest_ref_idx = (counts[sys_orig] * vsk_v[sys_orig].frame1_idx[k]) +
                                    sys_instance;
        const int f2_dest_ref_idx = (counts[sys_orig] * vsk_v[sys_orig].frame2_idx[k]) +
                                    sys_instance;
        const int f3_dest_ref_idx = (counts[sys_orig] * vsk_v[sys_orig].frame3_idx[k]) +
                                    sys_instance;
        const int f4_dest_ref_idx = (counts[sys_orig] * vsk_v[sys_orig].frame4_idx[k]) +
                                    sys_instance;
        tmp_vs_atoms[vsite_con] = atom_destinations[sys_orig][vs_dest_ref_idx];
        tmp_vs_frame_types[vsite_con] = vsk_v[sys_orig].vs_types[k];

        // Map the frame atom indices.  Some virtual sites will not have a frame3 or frame4 atom.
        // In these cases, the proper choice is to take the convention of the input topologies.
        tmp_vs_frame1_idx[vsite_con] = atom_destinations[sys_orig][f1_dest_ref_idx];
        tmp_vs_frame2_idx[vsite_con] = atom_destinations[sys_orig][f2_dest_ref_idx];
        tmp_vs_frame3_idx[vsite_con] = vsk_v[sys_orig].frame3_idx[k];
        tmp_vs_frame4_idx[vsite_con] = vsk_v[sys_orig].frame4_idx[k];
        switch (static_cast<VirtualSiteKind>(tmp_vs_frame_types[vsite_con])) {
        case VirtualSiteKind::FLEX_2:
        case VirtualSiteKind::FIXED_2:
        case VirtualSiteKind::NONE:
          break;
        case VirtualSiteKind::FLEX_3:
        case VirtualSiteKind::FIXED_3:
        case VirtualSiteKind::FAD_3:
        case VirtualSiteKind::OUT_3:
          tmp_vs_frame3_idx[vsite_con] = atom_destinations[sys_orig][f3_dest_ref_idx];
          break;
        case VirtualSiteKind::FIXED_4:
          tmp_vs_frame3_idx[vsite_con] = atom_destinations[sys_orig][f3_dest_ref_idx];
          tmp_vs_frame4_idx[vsite_con] = atom_destinations[sys_orig][f4_dest_ref_idx];
          break;
        }
        tmp_vs_parameter_indices[vsite_con] = vsk_v[sys_orig].vs_param_idx[k];
        vsite_con++;
      }
    }
  }
  VirtualSiteTable vst(atom_count, tmp_vs_atoms, tmp_vs_frame_types, tmp_vs_frame1_idx,
                       tmp_vs_frame2_idx, tmp_vs_frame3_idx, tmp_vs_frame4_idx,
                       tmp_vs_parameter_indices, uni_vsites.getUnion(0), uni_vsites.getUnion(1),
                       uni_vsites.getUnion(2));

  // Map exclusions in the new topology based on those present in the input topologies
  std::vector<int> tmp_raw_counts(atom_count, 0);
  std::vector<int> tmp_excl_bounds(atom_count + 1, 0);
  for (int i = 0; i < atom_count; i++) {
    const int sys_orig = atom_origins[i].x;
    const int idx_orig = atom_origins[i].y;
    const int* excl_ptr = agv[sys_orig]->atom_exclusion_bounds.data();
    const int i_excl = std::max(excl_ptr[idx_orig + 1] - excl_ptr[idx_orig], 1);
    tmp_raw_counts[i] = i_excl;
    tmp_excl_bounds[i] = total_exclusions;
    total_exclusions += i_excl;
  }
  tmp_excl_bounds[atom_count] = total_exclusions;
  std::vector<int> tmp_raw_exclusions(total_exclusions);
  int excl_con = 0;
  for (int i = 0; i < molecule_count; i++) {
    const int sys_instance = molecule_system_instance_origins[i].y;
    for (int j = tmp_molecule_limits[i]; j < tmp_molecule_limits[i + 1]; j++) {
      const int j_atom = tmp_molecule_contents[j];
      const int sys_orig = atom_origins[j_atom].x;
      const int idx_orig = atom_origins[j_atom].y;
      const int* excl_ptr = agv[sys_orig]->atom_exclusion_list.data();
      const int llim = agv[sys_orig]->atom_exclusion_bounds.readHost(idx_orig);
      const int hlim = agv[sys_orig]->atom_exclusion_bounds.readHost(idx_orig + 1);
      if (llim == hlim) {
        tmp_raw_exclusions[excl_con] = 0;
        excl_con++;
      }
      else {
        for (int k = llim; k < hlim; k++) {
          const int dest_ref_idx = (counts[sys_orig] * excl_ptr[k]) + sys_instance;
          tmp_raw_exclusions[excl_con] = atom_destinations[sys_orig][dest_ref_idx];
          excl_con++;
        }
      }
    }
  }
  const CondensedExclusions cond_excl = processExclusions(tmp_raw_counts, tmp_raw_exclusions,
                                                          source);
  const Map1234 all_nb_excl = mapExclusions(atom_count, tmp_vs_atoms, tmp_vs_frame1_idx,
                                            tmp_bond_i_atoms, tmp_bond_j_atoms);
  const AttenuationParameterSet attn_parm =
    condenseScreeningFactors(bvt, uni_dihes.getUnion(3), uni_dihes.getUnion(4),
                             agv[0]->elec14_screening_factor, agv[0]->vdw14_screening_factor);
  attenuated_14_type_count = attn_parm.total_14_sets;
  const ConstraintTable cnst_table(tmp_atomic_numbers, tmp_atomic_masses, tmp_molecule_limits,
                                   tmp_molecule_contents, tmp_molecule_membership, bvt,
                                   all_nb_excl, uni_bonds.getUnion(1), uni_angls.getUnion(1));

  // Examine dihedral coverge: are all 1:4 interactions covered by some dihedral, with a scaling
  // factor associated with those parameters?
  const std::vector<int3> outstanding_14_pairs =
    checkDihedral14Coverage(atom_count, tmp_atomic_numbers, bvt, all_nb_excl, vst, attn_parm,
                            policy);
  inferred_14_attenuations = outstanding_14_pairs.size();
  std::vector<int> tmp_inferred_14_i_atoms(inferred_14_attenuations);
  std::vector<int> tmp_inferred_14_l_atoms(inferred_14_attenuations);
  std::vector<int> tmp_inferred_14_param_idx(inferred_14_attenuations);
  for (int pos = 0; pos < inferred_14_attenuations; pos++) {
    tmp_inferred_14_i_atoms[pos]   = outstanding_14_pairs[pos].x;
    tmp_inferred_14_l_atoms[pos]   = outstanding_14_pairs[pos].y;
    tmp_inferred_14_param_idx[pos] = outstanding_14_pairs[pos].z;
  }

  // Miscellaneous items: atom mobility must be inherited from the input topologies, implicit
  // solvent information must be determined, obscure AtomGraph information on joining information
  // must be filled out.  The results of the Lennard Jones table fusion must be unpacked.
  std::vector<bool> tf_mobile_atoms(atom_count, true);
  const int n_packed_mobile_mask = (atom_count + 31) / 32;
  unused_natyp = tmp_desc[static_cast<size_t>(TopologyDescriptor::NATYP_UNUSED)];
  std::vector<uint> ui_mobile_atoms(n_packed_mobile_mask, 0);
  std::vector<int> tmp_lennard_jones_indices(atom_count);
  std::vector<int> tmp_charge_type_indices(atom_count);
  std::vector<int> tmp_neck_gb_indices(atom_count);
  std::vector<int> tmp_tree_joining_info(atom_count);
  std::vector<int> tmp_last_rotator_info(atom_count);
  std::vector<double> tmp_atomic_pb_radii(atom_count);
  std::vector<double> tmp_gb_screening_factors(atom_count);
  std::vector<char4> tmp_tree_symbols(atom_count);
  for (int i = 1; i < nsys; i++) {
    if (agv[i]->gb_style != agv[0]->gb_style) {
      rtErr("The implicit solvent model in use by topology " + getBaseName(agv[i]->source) +
            "(" + getEnumerationName(agv[i]->gb_style) + ") is inconsistent with that in use by "
            "topology " + getBaseName(agv[0]->source) + "(" +
            getEnumerationName(agv[0]->gb_style) + ").", "AtomGraph");
    }
  }
  charge_type_count = uni_chrgs.getUniqueParameterCount();
  for (int i = 0; i < atom_count; i++) {
    const int sys_orig = atom_origins[i].x;
    const int idx_orig = atom_origins[i].y;
    tf_mobile_atoms[i] = agv[sys_orig]->getAtomMobility(idx_orig);
    if (tf_mobile_atoms[i]) {
      const int i_elem = (i >> 5);
      const int i_shft = i - (i_elem << 5);
      ui_mobile_atoms[i_elem] |= (0x1 << i_shft);
    }

    // Find the Lennard-Jones atom type index in the original topology, then look up the
    // corresponding Lennard-Jones atom type index in the new topology.
    const int ljt_orig = agv[sys_orig]->getLennardJonesIndex(idx_orig);
    tmp_lennard_jones_indices[i] = uni_ljtab.getCorrespondence(sys_orig, ljt_orig);

    // Look up the corresponding charge type index in the new topology.
    const int qt_orig = agv[sys_orig]->getChargeIndex(idx_orig);
    tmp_charge_type_indices[i] = uni_chrgs.getCorrespondence(sys_orig, qt_orig);

    // Copy various implicit solvent properties of the atoms.
    tmp_atomic_pb_radii[i] = agv[sys_orig]->getAtomPBRadius<double>(idx_orig);
    tmp_gb_screening_factors[i] = agv[sys_orig]->getGBScreeningFactor<double>(idx_orig);
    tmp_neck_gb_indices[i] = agv[sys_orig]->neck_gb_indices.readHost(idx_orig);

    // Copy obscure symbology from the input topologies.
    tmp_tree_symbols[i] = agv[sys_orig]->tree_symbols.readHost(idx_orig);
  }
  std::vector<double> tmp_solty_info;
  for (int i = 0; i < nsys; i++) {
    const std::vector<double> isolty = agv[i]->solty_info.readHost();
    tmp_solty_info.insert(tmp_solty_info.end(), isolty.begin(), isolty.end());
  }
  for (int i = 0; i < n_packed_mobile_mask; i++) {
    ui_mobile_atoms[i] = ~ui_mobile_atoms[i]; 
  }
  const std::vector<int> tmp_mobile_atoms(ui_mobile_atoms.begin(), ui_mobile_atoms.end());
  lj_type_count = uni_ljtab.getLJTypeCount();
  const int sq_nlj = lj_type_count * lj_type_count;
  std::vector<double> tmp_lja_values(sq_nlj), tmp_lja_14_values(sq_nlj);
  std::vector<double> tmp_ljb_values(sq_nlj), tmp_ljb_14_values(sq_nlj);
  std::vector<double> tmp_ljc_values(sq_nlj), tmp_ljc_14_values(sq_nlj);
  std::vector<double> tmp_hbond_a_values(sq_nlj), tmp_hbond_b_values(sq_nlj);
  std::vector<double> tmp_hbond_cutoffs(sq_nlj);
  std::vector<double> tmp_gb_coef(atom_count, 0.0);
  for (int i = 0; i < lj_type_count; i++) {
    const std::vector<int2> i_examples = uni_ljtab.getInputInstances(i);
    const int niex = i_examples.size();
    for (int j = 0; j < lj_type_count; j++) {
      const double3 ij_parm = uni_ljtab.getLJCoefficients(i, j);
      const int ij_idx = (j * lj_type_count) + i;
      tmp_lja_values[ij_idx] = ij_parm.x;
      tmp_ljb_values[ij_idx] = ij_parm.y;
      tmp_ljc_values[ij_idx] = ij_parm.z;
      const double3 ij_14_parm = uni_ljtab.getLJ14Coefficients(i, j);
      tmp_lja_14_values[ij_idx] = ij_14_parm.x;
      tmp_ljb_14_values[ij_idx] = ij_14_parm.y;
      tmp_ljc_14_values[ij_idx] = ij_14_parm.z;

      // Hydrogen bonding is difficult to combine because the rules for deriving hydrogen bonding
      // parameters are dated or deprecated.  Hydrogen bonding parameters track Lennard-Jones atom
      // types in the input topologies.  Look up the Lennard-Jones atom types in the original
      // topologies, then find the corresponding hydrogen bonding parameters.  If the two atom
      // types do not coexist in one topology, neglect the hydrogen bonding interaction.
      const std::vector<int2> j_examples = uni_ljtab.getInputInstances(j);
      const int njex = j_examples.size();
      int shared_sys = -1;
      int shared_sys_iljt, shared_sys_jljt;
      for (int k = 0; k < niex; k++) {
        for (int m = 0; m < njex; m++) {
          if (i_examples[k].x == j_examples[m].x) {
            shared_sys = i_examples[k].x;
            shared_sys_iljt = i_examples[k].y;
            shared_sys_jljt = j_examples[m].y;
          }
        }
      }
      if (shared_sys >= 0) {
        const int km_idx = (shared_sys_jljt * agv[shared_sys]->lj_type_count) +
                           shared_sys_iljt;
        tmp_hbond_a_values[ij_idx] = agv[shared_sys]->hbond_a_values.readHost(km_idx);
        tmp_hbond_b_values[ij_idx] = agv[shared_sys]->hbond_b_values.readHost(km_idx);
        tmp_hbond_cutoffs[ij_idx] = agv[shared_sys]->hbond_cutoffs.readHost(km_idx);
      }
    }
  }
  
  // With all of the components collected, load the new topology's arrays
  loadHybridArrays(tmp_desc, tmp_residue_limits, tmp_atom_struc_numbers, tmp_residue_numbers,
                   tmp_molecule_limits, tmp_atomic_numbers, tmp_molecule_membership,
                   tmp_mobile_atoms, tmp_molecule_contents, tmp_cmap_surface_dims,
                   tmp_cmap_surface_bounds, tmp_charge_type_indices, tmp_lennard_jones_indices,
                   tmp_inferred_14_i_atoms, tmp_inferred_14_l_atoms, tmp_inferred_14_param_idx,
                   tmp_neck_gb_indices, tmp_tree_joining_info, tmp_last_rotator_info,
                   tmp_atomic_charges, tmp_atomic_masses, uni_ubrds.getUnion(0),
                   uni_ubrds.getUnion(1), uni_cimps.getUnion(0), uni_cimps.getUnion(1),
                   uni_cmaps.getAllSurfaces(), uni_bonds.getUnion(0),
                   uni_bonds.getUnion(1), uni_angls.getUnion(0), uni_angls.getUnion(1),
                   uni_dihes.getUnion(0), uni_dihes.getUnion(1), uni_dihes.getUnion(2),
                   uni_chrgs.getUnion(0), tmp_lja_values, tmp_ljb_values, tmp_ljc_values,
                   tmp_lja_14_values, tmp_ljb_14_values, tmp_ljc_14_values, tmp_atomic_pb_radii,
                   tmp_gb_screening_factors, tmp_gb_coef, tmp_solty_info, tmp_hbond_a_values,
                   tmp_hbond_b_values, tmp_hbond_cutoffs, tmp_atomic_polarizabilities,
                   tmp_atom_names, tmp_atom_types, tmp_residue_names, tmp_tree_symbols, cma,
                   cond_excl, bvt, cvt, attn_parm, vst, all_nb_excl, cnst_table);
}

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const AtomGraph &ag_a, const int n_a, const AtomGraph &ag_b, const int n_b,
                     const MoleculeOrdering arrangement, const ExceptionResponse policy,
                     const HybridFormat format_in) :
    AtomGraph({ const_cast<AtomGraph*>(ag_a.getSelfPointer()),
                const_cast<AtomGraph*>(ag_b.getSelfPointer()) }, { n_a, n_b }, arrangement, policy)
{}

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const AtomGraph &ag_a, const AtomGraph &ag_b, const int n_b,
                     const MoleculeOrdering arrangement, const ExceptionResponse policy,
                     const HybridFormat format_in) :
  AtomGraph({ const_cast<AtomGraph*>(ag_a.getSelfPointer()),
              const_cast<AtomGraph*>(ag_b.getSelfPointer()) }, { 1, n_b }, arrangement, policy)
{}

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const AtomGraph &ag_a, const AtomGraph &ag_b,
                     const MoleculeOrdering arrangement, const ExceptionResponse policy,
                     const HybridFormat format_in) :
  AtomGraph({ const_cast<AtomGraph*>(ag_a.getSelfPointer()),
              const_cast<AtomGraph*>(ag_b.getSelfPointer()) }, { 1, 1 }, arrangement, policy)
{}

//-------------------------------------------------------------------------------------------------
std::vector<int2> findMoleculeOrder(const std::vector<AtomGraph*> &agv,
                                    const std::vector<int> &counts, std::vector<int2> *origins,
                                    const MoleculeOrdering arrangement) {
  int molecule_count = 0;
  const	size_t nsys = agv.size();
  if (nsys != counts.size()) {
    rtErr("A number of copies must be provided for each topology instance involved in the "
          "recombination (" + std::to_string(agv.size()) + " topologies and " +
          std::to_string(counts.size()) + " counts of each system were provided).",
          "findMoleculeOrder");
  }
  for (size_t i = 0; i < nsys; i++) {
    molecule_count += agv[i]->getMoleculeCount() * counts[i];
  }
  if (origins != nullptr) {
    origins->resize(molecule_count);
  }
  std::vector<int2> result(molecule_count);
  std::vector<std::vector<MoleculeKind>> molecule_kinds(nsys);
  std::vector<std::vector<int>> molecule_usage(nsys);
  for (int i = 0; i < nsys; i++) {
    molecule_kinds[i].resize(agv[i]->getMoleculeCount());
    for (int j = 0; j < agv[i]->getMoleculeCount(); j++) {
      molecule_kinds[i][j] = agv[i]->getMoleculeKind(j);
    }
    if (origins != nullptr) {
      molecule_usage[i] = std::vector<int>(agv[i]->getMoleculeCount(), 0);
    }
  }
  const int n_mol_type = static_cast<int>(MoleculeKind::OTHER);
  int molcon = 0;
  switch (arrangement) {
  case MoleculeOrdering::RETAIN_ORDER:
    {
      for (int i = 0; i < nsys; i++) {
        for (int j = 0; j < counts[j]; j++) {
          for (int k = 0; k < agv[i]->getMoleculeCount(); k++) {
            result[molcon] = { i, k };
            if (origins != nullptr) {
              origins->at(molcon) = { k, molecule_usage[i][k] };
              molecule_usage[i][k] += 1;
            }
            molcon++;
          }
        }
      }
    }
    break;
  case MoleculeOrdering::WATER_LAST:
    for (int i = 0; i < nsys; i++) {
      for (int j = 0; j < agv[i]->getMoleculeCount(); j++) {
        if (molecule_kinds[i][j] != MoleculeKind::WATER) {
          for (int k = 0; k < counts[j]; k++) {
            result[molcon] = { i, j };
            if (origins != nullptr) {
              origins->at(molcon) = { j, molecule_usage[i][j] };
              molecule_usage[i][j] += 1;
            }
            molcon++;
          }
        }
      }
      for (int j = 0; j < agv[i]->getMoleculeCount(); j++) {
        if (molecule_kinds[i][j] == MoleculeKind::WATER) {
          for (int k = 0; k < counts[j]; k++) {
            result[molcon] = { i, j };
            if (origins != nullptr) {
              origins->at(molcon) = { i, molecule_usage[i][j] };
              molecule_usage[i][j] += 1;
            }
            molcon++;
          }
        }
      }
      
    }
    break;
  case MoleculeOrdering::REORDER_ALL:
    for (int i = 0; i <= n_mol_type; i++) {
      const MoleculeKind iseek = static_cast<MoleculeKind>(i);
      for (int j = 0; j < nsys; j++) {
        for (int k = 0; k < agv[j]->getMoleculeCount(); k++) {
          if (molecule_kinds[j][k] == iseek) {
            for (int m = 0; m < counts[j]; m++) {
              result[molcon] = { j, k };
              if (origins != nullptr) {
                origins->at(molcon) = { j, molecule_usage[j][k] };
                molecule_usage[j][k] += 1;
              }
              molcon++;
            }
          }
        }
      }
    }
    break;
  }

  return result;
}

} // namespace topology
} // namespace stormm
