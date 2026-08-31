#include <algorithm>
#include <cmath>
#include <cstdio>
#include <climits>
#include <ctime>
#include <unistd.h>
#include "copyright.h"
#include "FileManagement/file_util.h"
#include "Math/series_ops.h"
#include "Math/summation.h"
#include "Reporting/error_format.h"
#include "atomgraph.h"

namespace stormm {
namespace topology {

using card::HybridTargetLevel;
using card::HybridKind;
using diskutil::detectTopologyKind;
using stmath::prefixSumInPlace;
using stmath::PrefixSumType;
using stmath::sum;

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const HybridFormat format_in) :
    format{format_in},
    version_stamp{""},
    date{0, 0, 0, 0, 0, 0, 0, 0, 0},
    title{""},
    source{""},
    force_fields{},

    // Counts of atoms, residues, and other parts of the system
    atom_count{0}, residue_count{0}, molecule_count{0}, largest_residue_size{0},
    last_solute_residue{0}, last_solute_atom{0}, first_solvent_molecule{0},
    last_atom_before_cap{0}, implicit_copy_count{0}, largest_molecule_size{0},
    water_residue_size{0}, water_residue_count{0}, unconstrained_dof{0}, constrained_dof{0},
    descriptors{HybridKind::POINTER, "tp_desc", format},
    residue_limits{HybridKind::POINTER, "tp_res_limits", format},
    atom_struc_numbers{HybridKind::POINTER, "tp_atom_struc_nums", format},
    residue_numbers{HybridKind::POINTER, "tp_res_struc_nums", format},
    molecule_limits{HybridKind::POINTER, "tp_mol_limits", format},

    // Atom and residue details
    atomic_numbers{HybridKind::POINTER, "tp_znum", format},
    mobile_atoms{HybridKind::POINTER, "tp_belly", format},
    molecule_membership{HybridKind::POINTER, "tp_molnum", format},
    molecule_contents{HybridKind::POINTER, "tp_mol_contents", format},
    atomic_charges{HybridKind::POINTER, "tp_atomq", format},
    atomic_masses{HybridKind::POINTER, "tp_mass", format},
    inverse_atomic_masses{HybridKind::POINTER, "tp_invmass", format},
    sp_atomic_charges{HybridKind::POINTER, "tp_atomq_sp", format},
    sp_atomic_masses{HybridKind::POINTER, "tp_mass_sp", format},
    sp_inverse_atomic_masses{HybridKind::POINTER, "tp_invmass_sp", format},
    atom_names{HybridKind::POINTER, "tp_atom_names", format},
    atom_types{HybridKind::POINTER, "tp_atom_types", format},
    residue_names{HybridKind::POINTER, "tp_res_names", format},

    // CHARMM force field family parameters
    urey_bradley_term_count{0}, charmm_impr_term_count{0}, cmap_term_count{0},
    urey_bradley_parameter_count{0}, charmm_impr_parameter_count{0}, cmap_surface_count{0},
    urey_bradley_pert_term_count{0}, charmm_impr_pert_term_count{0}, cmap_pert_term_count{0},
    urey_bradleys_in_perturbed_group{0}, charmm_imprs_in_perturbed_group{0},
    cmaps_in_perturbed_group{0},
    urey_bradley_i_atoms{HybridKind::POINTER, "tp_urey_i", format},
    urey_bradley_k_atoms{HybridKind::POINTER, "tp_urey_k", format},
    urey_bradley_parameter_indices{HybridKind::POINTER, "tp_urey_parm", format},
    urey_bradley_assigned_atoms{HybridKind::POINTER, "tp_urey_asatoms", format},
    urey_bradley_assigned_index{HybridKind::POINTER, "tp_urey_asindex", format},
    urey_bradley_assigned_terms{HybridKind::POINTER, "tp_urey_asterms", format},
    urey_bradley_assigned_bounds{HybridKind::POINTER, "tp_urey_asbounds", format},
    charmm_impr_i_atoms{HybridKind::POINTER, "tp_charmm_impr_i", format},
    charmm_impr_j_atoms{HybridKind::POINTER, "tp_charmm_impr_j", format},
    charmm_impr_k_atoms{HybridKind::POINTER, "tp_charmm_impr_k", format},
    charmm_impr_l_atoms{HybridKind::POINTER, "tp_charmm_impr_l", format},
    charmm_impr_parameter_indices{HybridKind::POINTER, "tp_charmm_impr_parm", format},
    charmm_impr_assigned_atoms{HybridKind::POINTER, "tp_cimpr_asatoms", format},
    charmm_impr_assigned_index{HybridKind::POINTER, "tp_cimpr_asindex", format},
    charmm_impr_assigned_terms{HybridKind::POINTER, "tp_cimpr_asterms", format},
    charmm_impr_assigned_bounds{HybridKind::POINTER, "tp_cimpr_asbounds", format},
    cmap_i_atoms{HybridKind::POINTER, "tp_cmap_i", format},
    cmap_j_atoms{HybridKind::POINTER, "tp_cmap_j", format},
    cmap_k_atoms{HybridKind::POINTER, "tp_cmap_k", format},
    cmap_l_atoms{HybridKind::POINTER, "tp_cmap_l", format},
    cmap_m_atoms{HybridKind::POINTER, "tp_cmap_m", format},
    cmap_surface_dimensions{HybridKind::POINTER, "tp_cmap_dims", format},
    cmap_surface_bounds{HybridKind::POINTER, "tp_cmap_bounds", format},
    cmap_patch_bounds{HybridKind::POINTER, "tp_cmap_patch_bounds", format},
    cmap_surface_indices{HybridKind::POINTER, "tp_cmap_parm", format},
    cmap_assigned_atoms{HybridKind::POINTER, "tp_cmap_asatoms", format},
    cmap_assigned_index{HybridKind::POINTER, "tp_cmap_asindex", format},
    cmap_assigned_terms{HybridKind::POINTER, "tp_cmap_asterms", format},
    cmap_assigned_bounds{HybridKind::POINTER, "tp_cmap_asbounds", format},
    urey_bradley_stiffnesses{HybridKind::POINTER, "tp_ub_stiffness", format},
    urey_bradley_equilibria{HybridKind::POINTER, "tp_ub_equilibria", format},
    charmm_impr_stiffnesses{HybridKind::POINTER, "tp_cimp_stiffness", format},
    charmm_impr_phase_angles{HybridKind::POINTER, "tp_cimp_equilibria", format},
    cmap_surfaces{HybridKind::POINTER, "tp_cmap_surface", format},
    cmap_phi_derivatives{HybridKind::POINTER, "tp_cmap_dphi", format},
    cmap_psi_derivatives{HybridKind::POINTER, "tp_cmap_dphi", format},
    cmap_phi_psi_derivatives{HybridKind::POINTER, "tp_cmap_dphi_dpsi", format},
    cmap_patches{HybridKind::POINTER, "tp_cmap_patches", format},
    sp_urey_bradley_stiffnesses{HybridKind::POINTER, "tp_ub_stiffness_sp", format},
    sp_urey_bradley_equilibria{HybridKind::POINTER, "tp_ub_equilibria_sp", format},
    sp_charmm_impr_stiffnesses{HybridKind::POINTER, "tp_ub_stiffness_sp", format},
    sp_charmm_impr_phase_angles{HybridKind::POINTER, "tp_ub_equilibria_sp", format},
    sp_cmap_surfaces{HybridKind::POINTER, "tp_cmap_surface_sp", format},
    sp_cmap_phi_derivatives{HybridKind::POINTER, "tp_cmap_dphi", format},
    sp_cmap_psi_derivatives{HybridKind::POINTER, "tp_cmap_dphi", format},
    sp_cmap_phi_psi_derivatives{HybridKind::POINTER, "tp_cmap_dphi_dpsi", format},
    sp_cmap_patches{HybridKind::POINTER, "tp_cmap_patches", format},

    // Relevant information for the bonded calculation
    bond_term_with_hydrogen{0}, angl_term_with_hydrogen{0}, dihe_term_with_hydrogen{0},
    bond_term_without_hydrogen{0}, angl_term_without_hydrogen{0}, dihe_term_without_hydrogen{0},
    bond_term_count{0}, angl_term_count{0}, dihe_term_count{0}, bond_parameter_count{0},
    angl_parameter_count{0}, dihe_parameter_count{0}, bond_perturbation_term_count{0},
    angl_perturbation_term_count{0}, dihe_perturbation_term_count{0}, bonds_in_perturbed_group{0},
    angls_in_perturbed_group{0}, dihes_in_perturbed_group{0}, bonded_group_count{0},
    bond_stiffnesses{HybridKind::POINTER, "tp_bondk", format},
    bond_equilibria{HybridKind::POINTER, "tp_bondl0", format},
    angl_stiffnesses{HybridKind::POINTER, "tp_anglk", format},
    angl_equilibria{HybridKind::POINTER, "tp_anglt0", format},
    dihe_amplitudes{HybridKind::POINTER, "tp_dihek", format},
    dihe_periodicities{HybridKind::POINTER, "tp_dihen", format},
    dihe_phase_angles{HybridKind::POINTER, "tp_dihepsi", format},
    sp_bond_stiffnesses{HybridKind::POINTER, "tp_bondk_sp", format},
    sp_bond_equilibria{HybridKind::POINTER, "tp_bondl0_sp", format},
    sp_angl_stiffnesses{HybridKind::POINTER, "tp_anglk_sp", format},
    sp_angl_equilibria{HybridKind::POINTER, "tp_anglt0_sp", format},
    sp_dihe_amplitudes{HybridKind::POINTER, "tp_dihek_sp", format},
    sp_dihe_periodicities{HybridKind::POINTER, "tp_dihen_sp", format},
    sp_dihe_phase_angles{HybridKind::POINTER, "tp_dihepsi_sp", format},
    bond_i_atoms{HybridKind::POINTER, "tp_bond_i", format},
    bond_j_atoms{HybridKind::POINTER, "tp_bond_j", format},
    bond_parameter_indices{HybridKind::POINTER, "tp_bond_parm", format},
    bond_assigned_atoms{HybridKind::POINTER, "tp_bond_asatoms", format},
    bond_assigned_index{HybridKind::POINTER, "tp_bond_asindex", format},
    bond_assigned_terms{HybridKind::POINTER, "tp_bond_asterms", format},
    bond_assigned_bounds{HybridKind::POINTER, "tp_bond_asbounds", format},
    angl_i_atoms{HybridKind::POINTER, "tp_angl_i", format},
    angl_j_atoms{HybridKind::POINTER, "tp_angl_j", format},
    angl_k_atoms{HybridKind::POINTER, "tp_angl_k", format},
    angl_parameter_indices{HybridKind::POINTER, "tp_angl_parm", format},
    angl_assigned_atoms{HybridKind::POINTER, "tp_angl_asatoms", format},
    angl_assigned_index{HybridKind::POINTER, "tp_angl_asindex", format},
    angl_assigned_terms{HybridKind::POINTER, "tp_angl_asterms", format},
    angl_assigned_bounds{HybridKind::POINTER, "tp_angl_asbounds", format},
    dihe_i_atoms{HybridKind::POINTER, "tp_dihe_i", format},
    dihe_j_atoms{HybridKind::POINTER, "tp_dihe_j", format},
    dihe_k_atoms{HybridKind::POINTER, "tp_dihe_k", format},
    dihe_l_atoms{HybridKind::POINTER, "tp_dihe_l", format},
    dihe_parameter_indices{HybridKind::POINTER, "tp_dihe_parm", format},
    dihe14_parameter_indices{HybridKind::POINTER, "tp_dihe_attn14_parm", format},
    dihe_assigned_atoms{HybridKind::POINTER, "tp_dihe_asatoms", format},
    dihe_assigned_index{HybridKind::POINTER, "tp_dihe_asindex", format},
    dihe_assigned_terms{HybridKind::POINTER, "tp_dihe_asterms", format},
    dihe_assigned_bounds{HybridKind::POINTER, "tp_dihe_asbounds", format},
    bond_evaluation_mask{HybridKind::POINTER, "tp_bond_eval", format},
    angl_evaluation_mask{HybridKind::POINTER, "tp_angl_eval", format},
    bond_modifiers{HybridKind::POINTER, "tp_bond_mods", format},
    angl_modifiers{HybridKind::POINTER, "tp_angl_mods", format},
    dihe_modifiers{HybridKind::POINTER, "tp_dihe_mods", format},
    bond_assigned_mods{HybridKind::POINTER, "tp_bond_asmods", format},
    angl_assigned_mods{HybridKind::POINTER, "tp_angl_asmods", format},
    dihe_assigned_mods{HybridKind::POINTER, "tp_dihe_asmods", format},

    // Information relevant to virtual site placement
    virtual_site_count{0}, virtual_site_parameter_set_count{0},
    virtual_site_atoms{HybridKind::POINTER, "tp_vsidx", format},
    virtual_site_frame_types{HybridKind::POINTER, "tp_vsfrm", format},
    virtual_site_frame1_atoms{HybridKind::POINTER, "tp_vsfr1", format},
    virtual_site_frame2_atoms{HybridKind::POINTER, "tp_vsfr2", format},
    virtual_site_frame3_atoms{HybridKind::POINTER, "tp_vsfr3", format},
    virtual_site_frame4_atoms{HybridKind::POINTER, "tp_vsfr4", format},
    virtual_site_parameter_indices{HybridKind::POINTER, "tp_vs_parm_idx", format},
    virtual_site_frame_dim1{HybridKind::POINTER, "tp_vsdim1", format},
    virtual_site_frame_dim2{HybridKind::POINTER, "tp_vsdim2", format},
    virtual_site_frame_dim3{HybridKind::POINTER, "tp_vsdim3", format},
    sp_virtual_site_frame_dim1{HybridKind::POINTER, "tp_vsdim1_sp", format},
    sp_virtual_site_frame_dim2{HybridKind::POINTER, "tp_vsdim2_sp", format},
    sp_virtual_site_frame_dim3{HybridKind::POINTER, "tp_vsdim3_sp", format},

    // Relevant information for the non-bonded calculation
    charge_type_count{0}, lj_type_count{0}, has_14_lennard_jones_data{false}, total_exclusions{0},
    attenuated_14_type_count{0}, inferred_14_attenuations{0},
    periodic_box_class{UnitCellType::NONE},
    gb_style{ImplicitSolventModel::NONE}, dielectric_constant{1.0}, salt_concentration{0.0},
    coulomb_constant{amber_ancient_bioq}, pb_radii_set{""},
    charge_indices{HybridKind::POINTER, "tp_qidx", format},
    lennard_jones_indices{HybridKind::POINTER, "tp_ljidx", format},
    atom_exclusion_bounds{HybridKind::POINTER, "tp_nexcl", format},
    atom_exclusion_list{HybridKind::POINTER, "tp_excllist", format},
    nb11_exclusion_bounds{HybridKind::POINTER, "tp_nnb11", format},
    nb11_exclusion_list{HybridKind::POINTER, "tp_nb11_excl", format},
    nb12_exclusion_bounds{HybridKind::POINTER, "tp_nnb12", format},
    nb12_exclusion_list{HybridKind::POINTER, "tp_nb12_excl", format},
    nb13_exclusion_bounds{HybridKind::POINTER, "tp_nnb13", format},
    nb13_exclusion_list{HybridKind::POINTER, "tp_nb13_excl", format},
    nb14_exclusion_bounds{HybridKind::POINTER, "tp_nnb14", format},
    nb14_exclusion_list{HybridKind::POINTER, "tp_nb14_excl", format},
    infr14_i_atoms{HybridKind::POINTER, "tp_inferred14_i", format},
    infr14_l_atoms{HybridKind::POINTER, "tp_inferred14_j", format},
    infr14_parameter_indices{HybridKind::POINTER, "tp_inferred14_param", format},
    neck_gb_indices{HybridKind::POINTER, "tp_gbneck_idx", format},
    charge_parameters{HybridKind::POINTER, "tp_qparam", format},
    lj_a_values{HybridKind::POINTER, "tp_lja", format},
    lj_b_values{HybridKind::POINTER, "tp_ljb", format},
    lj_c_values{HybridKind::POINTER, "tp_ljc", format},
    lj_14_a_values{HybridKind::POINTER, "tp_14_lja", format},
    lj_14_b_values{HybridKind::POINTER, "tp_14_ljb", format},
    lj_14_c_values{HybridKind::POINTER, "tp_14_ljc", format},
    lj_sigma_values{HybridKind::POINTER, "tp_sigma", format},
    lj_14_sigma_values{HybridKind::POINTER, "tp_14_sigma", format},
    lj_type_corrections{HybridKind::POINTER, "tp_lj_long", format},
    attn14_elec_factors{HybridKind::POINTER, "tp_scee_param", format},
    attn14_vdw_factors{HybridKind::POINTER, "tp_scnb_param", format},
    atomic_pb_radii{HybridKind::POINTER, "tp_pbradii", format},
    gb_screening_factors{HybridKind::POINTER, "tp_screen", format},
    gb_alpha_parameters{HybridKind::POINTER, "tp_gb_alpha_sp", format},
    gb_beta_parameters{HybridKind::POINTER, "tp_gb_beta_sp", format},
    gb_gamma_parameters{HybridKind::POINTER, "tp_gb_gamma_sp", format},
    sp_charge_parameters{HybridKind::POINTER, "tp_qparam_sp", format},
    sp_lj_a_values{HybridKind::POINTER, "tp_lja_sp", format},
    sp_lj_b_values{HybridKind::POINTER, "tp_ljb_sp", format},
    sp_lj_c_values{HybridKind::POINTER, "tp_ljc_sp", format},
    sp_lj_14_a_values{HybridKind::POINTER, "tp_lja_sp", format},
    sp_lj_14_b_values{HybridKind::POINTER, "tp_ljb_sp", format},
    sp_lj_14_c_values{HybridKind::POINTER, "tp_ljc_sp", format},
    sp_lj_sigma_values{HybridKind::POINTER, "tp_sigma", format},
    sp_lj_14_sigma_values{HybridKind::POINTER, "tp_14_sigma", format},
    sp_lj_type_corrections{HybridKind::POINTER, "tp_lj_long_sp", format},
    sp_attn14_elec_factors{HybridKind::POINTER, "tp_scee_param_sp", format},
    sp_attn14_vdw_factors{HybridKind::POINTER, "tp_scnb_param_sp", format},
    sp_atomic_pb_radii{HybridKind::POINTER, "tp_pbradii_sp", format},
    sp_gb_screening_factors{HybridKind::POINTER, "tp_screen_sp", format},
    sp_gb_alpha_parameters{HybridKind::POINTER, "tp_gb_alpha_sp", format},
    sp_gb_beta_parameters{HybridKind::POINTER, "tp_gb_beta_sp", format},
    sp_gb_gamma_parameters{HybridKind::POINTER, "tp_gb_gamma_sp", format},

    // MD propagation algorithm directives
    use_shake{ApplyConstraints::NO}, use_settle{ApplyConstraints::NO},
    use_perturbation_info{PerturbationSetting::OFF}, use_solvent_cap_option{SolventCapSetting::ON},
    use_polarization{PolarizationSetting::ON}, water_residue_name{' ', ' ', ' ', ' '},
    bond_constraint_mask{""}, bond_constraint_omit_mask{""}, bond_constraint_count{0},
    nonrigid_particle_count{0}, settle_group_count{0},
    settle_parameter_count{0}, constraint_group_count{0}, constraint_parameter_count{0},
    settle_oxygen_atoms{HybridKind::POINTER, "tp_settle_ox", format},
    settle_hydro1_atoms{HybridKind::POINTER, "tp_settle_h1", format},
    settle_hydro2_atoms{HybridKind::POINTER, "tp_settle_h2", format},
    settle_parameter_indices{HybridKind::POINTER, "tp_sett_param_idx", format},
    constraint_group_atoms{HybridKind::POINTER, "tp_cnst_atoms", format},
    constraint_group_bounds{HybridKind::POINTER, "tp_cnst_bounds", format},
    constraint_parameter_indices{HybridKind::POINTER, "tp_cnst_param_idx", format},
    constraint_parameter_bounds{HybridKind::POINTER, "tp_cnst_param_bnds", format},
    settle_mo{HybridKind::POINTER, "tp_settle_mo", format},
    settle_mh{HybridKind::POINTER, "tp_settle_mh", format},
    settle_moh{HybridKind::POINTER, "tp_settle_moh", format},
    settle_mormt{HybridKind::POINTER, "tp_settle_mormt", format},
    settle_mhrmt{HybridKind::POINTER, "tp_settle_mhrmt", format},
    settle_ra{HybridKind::POINTER, "tp_settle_ra", format},
    settle_rb{HybridKind::POINTER, "tp_settle_rb", format},
    settle_rc{HybridKind::POINTER, "tp_settle_rc", format},
    constraint_inverse_masses{HybridKind::POINTER, "tp_cnst_invms", format},
    constraint_squared_lengths{HybridKind::POINTER, "tp_cnst_targets", format},
    sp_settle_mo{HybridKind::POINTER, "tp_settle_mo_sp", format},
    sp_settle_mh{HybridKind::POINTER, "tp_settle_mh_sp", format},
    sp_settle_moh{HybridKind::POINTER, "tp_settle_moh_sp", format},
    sp_settle_mormt{HybridKind::POINTER, "tp_settle_mormt_sp", format},
    sp_settle_mhrmt{HybridKind::POINTER, "tp_settle_mhrmt_sp", format},
    sp_settle_ra{HybridKind::POINTER, "tp_settle_ra_sp", format},
    sp_settle_rb{HybridKind::POINTER, "tp_settle_rb_sp", format},
    sp_settle_rc{HybridKind::POINTER, "tp_settle_rc_sp", format},
    sp_constraint_inverse_masses{HybridKind::POINTER, "tp_cnst_invms", format},
    sp_constraint_squared_lengths{HybridKind::POINTER, "tp_cnst_targets", format},

    // Overflow name keys
    atom_overflow_names{HybridKind::POINTER, "atom_name_xtkey", format},
    atom_overflow_types{HybridKind::POINTER, "atom_type_xtkey", format},
    residue_overflow_names{HybridKind::POINTER, "residue_name_xtkey", format},

    // Information currently unused
    unused_nhparm{0}, unused_nparm{0}, unused_natyp{0}, hbond_10_12_parameter_count{0},
    heavy_bonds_plus_constraints{0}, heavy_angls_plus_constraints{0},
    heavy_dihes_plus_constraints{0},
    tree_joining_info{HybridKind::POINTER, "tp_join", format},
    last_rotator_info{HybridKind::POINTER, "tp_irotat", format},
    solty_info{HybridKind::POINTER, "tp_solty", format},
    hbond_a_values{HybridKind::POINTER, "tp_sola", format},
    hbond_b_values{HybridKind::POINTER, "tp_solb", format},
    hbond_cutoffs{HybridKind::POINTER, "tp_hbcut", format},
    atomic_polarizabilities{HybridKind::POINTER, "tp_polarizability", format},
    tree_symbols{HybridKind::POINTER, "tp_symbl", format},

    // Hybrid data structures (actual arrays)
    int_data{HybridKind::ARRAY, "tp_int", format},
    double_data{HybridKind::ARRAY, "tp_double", format},
    float_data{HybridKind::ARRAY, "tp_float", format},
    char4_data{HybridKind::ARRAY, "tp_char4", format}
{
  // The format of the Hybrid arrays must include a component of memory on the CPU host.
#ifdef STORMM_USE_HPC
  switch (format) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::HOST_ONLY:
    break;
  case HybridFormat::DEVICE_ONLY:
    rtErr("Memory for a basic topology must include data on the CPU host.", "Hybrid");    
  }
#endif
}

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const std::string &file_name, const ExceptionResponse policy,
                     const TopologyKind engine_format, const HybridFormat format_in) :
    AtomGraph(format_in)
{
  TopologyKind detected_engine_format;
  switch (engine_format) {
  case TopologyKind::AMBER:
  case TopologyKind::CHARMM:
  case TopologyKind::GROMACS:
  case TopologyKind::OPENMM:
    detected_engine_format = engine_format;
    break;
  case TopologyKind::UNKNOWN:
    detected_engine_format = detectTopologyKind(file_name, "AtomGraph");
    break;
  }
  switch (detected_engine_format) {
  case TopologyKind::AMBER:
    buildFromPrmtop(file_name, policy);
    break;
  case TopologyKind::CHARMM:
  case TopologyKind::GROMACS:
  case TopologyKind::OPENMM:
    rtErr("Construction from " + getEnumerationName(detected_engine_format) + " format files is "
          "not yet implemented.", "AtomGraph");
  case TopologyKind::UNKNOWN:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const std::string &file_name, const ExceptionResponse policy,
                     const TopologyKind engine_format, const double coulomb_constant_in,
                     const double default_elec14_screening, const double default_vdw14_screening,
                     const double charge_rounding_tol, const double charge_discretization,
                     const ApplyConstraints use_shake_in, const ApplyConstraints use_settle_in,
                     const HybridFormat format_in) :
    AtomGraph(format_in)
{
  TopologyKind detected_engine_format;
  use_shake  = use_shake_in;
  use_settle = use_settle_in;
  switch (engine_format) {
  case TopologyKind::AMBER:
  case TopologyKind::CHARMM:
  case TopologyKind::GROMACS:
  case TopologyKind::OPENMM:
    detected_engine_format = engine_format;
    break;
  case TopologyKind::UNKNOWN:
    detected_engine_format = detectTopologyKind(file_name, "AtomGraph");
    break;
  }
  switch (detected_engine_format) {
  case TopologyKind::AMBER:
    buildFromPrmtop(file_name, policy, coulomb_constant_in, default_elec14_screening,
                    default_vdw14_screening, charge_rounding_tol, charge_discretization);
    break;
  case TopologyKind::CHARMM:
  case TopologyKind::GROMACS:
  case TopologyKind::OPENMM:
  case TopologyKind::UNKNOWN:
    rtErr("Construction from non-Amber format files is not yet implemented.", "AtomGraph");
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::rebasePointers() {

  // Repair counts of atoms, residues, and other parts of the system
  descriptors.swapTarget(&int_data);
  residue_limits.swapTarget(&int_data);
  atom_struc_numbers.swapTarget(&int_data);
  residue_numbers.swapTarget(&int_data);
  molecule_limits.swapTarget(&int_data);

  // Atom and residue details
  atomic_numbers.swapTarget(&int_data);
  mobile_atoms.swapTarget(&int_data);
  molecule_membership.swapTarget(&int_data);
  molecule_contents.swapTarget(&int_data);
  atomic_charges.swapTarget(&double_data);
  atomic_masses.swapTarget(&double_data);
  inverse_atomic_masses.swapTarget(&double_data);
  sp_atomic_charges.swapTarget(&float_data);
  sp_atomic_masses.swapTarget(&float_data);
  sp_inverse_atomic_masses.swapTarget(&float_data);
  atom_names.swapTarget(&char4_data);
  atom_types.swapTarget(&char4_data);
  residue_names.swapTarget(&char4_data);

  // CHARMM force field family parameters
  urey_bradley_i_atoms.swapTarget(&int_data);
  urey_bradley_k_atoms.swapTarget(&int_data);
  urey_bradley_parameter_indices.swapTarget(&int_data);
  urey_bradley_assigned_atoms.swapTarget(&int_data);
  urey_bradley_assigned_index.swapTarget(&int_data);
  urey_bradley_assigned_terms.swapTarget(&int_data);
  urey_bradley_assigned_bounds.swapTarget(&int_data);
  charmm_impr_i_atoms.swapTarget(&int_data);
  charmm_impr_j_atoms.swapTarget(&int_data);
  charmm_impr_k_atoms.swapTarget(&int_data);
  charmm_impr_l_atoms.swapTarget(&int_data);
  charmm_impr_parameter_indices.swapTarget(&int_data);
  charmm_impr_assigned_atoms.swapTarget(&int_data);
  charmm_impr_assigned_index.swapTarget(&int_data);
  charmm_impr_assigned_terms.swapTarget(&int_data);
  charmm_impr_assigned_bounds.swapTarget(&int_data);
  cmap_i_atoms.swapTarget(&int_data);
  cmap_j_atoms.swapTarget(&int_data);
  cmap_k_atoms.swapTarget(&int_data);
  cmap_l_atoms.swapTarget(&int_data);
  cmap_m_atoms.swapTarget(&int_data);
  cmap_surface_dimensions.swapTarget(&int_data);
  cmap_surface_bounds.swapTarget(&int_data);
  cmap_patch_bounds.swapTarget(&int_data);
  cmap_surface_indices.swapTarget(&int_data);
  cmap_assigned_atoms.swapTarget(&int_data);
  cmap_assigned_index.swapTarget(&int_data);
  cmap_assigned_terms.swapTarget(&int_data);
  cmap_assigned_bounds.swapTarget(&int_data);
  urey_bradley_stiffnesses.swapTarget(&double_data);
  urey_bradley_equilibria.swapTarget(&double_data);
  charmm_impr_stiffnesses.swapTarget(&double_data);
  charmm_impr_phase_angles.swapTarget(&double_data);
  cmap_surfaces.swapTarget(&double_data);
  cmap_phi_derivatives.swapTarget(&double_data);
  cmap_psi_derivatives.swapTarget(&double_data);
  cmap_phi_psi_derivatives.swapTarget(&double_data);
  cmap_patches.swapTarget(&double_data);
  sp_urey_bradley_stiffnesses.swapTarget(&float_data);
  sp_urey_bradley_equilibria.swapTarget(&float_data);
  sp_charmm_impr_stiffnesses.swapTarget(&float_data);
  sp_charmm_impr_phase_angles.swapTarget(&float_data);
  sp_cmap_surfaces.swapTarget(&float_data);
  sp_cmap_phi_derivatives.swapTarget(&float_data);
  sp_cmap_psi_derivatives.swapTarget(&float_data);
  sp_cmap_phi_psi_derivatives.swapTarget(&float_data);
  sp_cmap_patches.swapTarget(&float_data);

  // Repair pointers relevant to the basic force field bonded terms
  bond_stiffnesses.swapTarget(&double_data);
  bond_equilibria.swapTarget(&double_data);
  angl_stiffnesses.swapTarget(&double_data);
  angl_equilibria.swapTarget(&double_data);
  dihe_amplitudes.swapTarget(&double_data);
  dihe_periodicities.swapTarget(&double_data);
  dihe_phase_angles.swapTarget(&double_data);
  sp_bond_stiffnesses.swapTarget(&float_data);
  sp_bond_equilibria.swapTarget(&float_data);
  sp_angl_stiffnesses.swapTarget(&float_data);
  sp_angl_equilibria.swapTarget(&float_data);
  sp_dihe_amplitudes.swapTarget(&float_data);
  sp_dihe_periodicities.swapTarget(&float_data);
  sp_dihe_phase_angles.swapTarget(&float_data);
  bond_i_atoms.swapTarget(&int_data);
  bond_j_atoms.swapTarget(&int_data);
  bond_parameter_indices.swapTarget(&int_data);
  bond_assigned_atoms.swapTarget(&int_data);
  bond_assigned_index.swapTarget(&int_data);
  bond_assigned_terms.swapTarget(&int_data);
  bond_assigned_bounds.swapTarget(&int_data);
  angl_i_atoms.swapTarget(&int_data);
  angl_j_atoms.swapTarget(&int_data);
  angl_k_atoms.swapTarget(&int_data);
  angl_parameter_indices.swapTarget(&int_data);
  angl_assigned_atoms.swapTarget(&int_data);
  angl_assigned_index.swapTarget(&int_data);
  angl_assigned_terms.swapTarget(&int_data);
  angl_assigned_bounds.swapTarget(&int_data);
  dihe_i_atoms.swapTarget(&int_data);
  dihe_j_atoms.swapTarget(&int_data);
  dihe_k_atoms.swapTarget(&int_data);
  dihe_l_atoms.swapTarget(&int_data);
  dihe_parameter_indices.swapTarget(&int_data);
  dihe14_parameter_indices.swapTarget(&int_data);
  dihe_assigned_atoms.swapTarget(&int_data);
  dihe_assigned_index.swapTarget(&int_data);
  dihe_assigned_terms.swapTarget(&int_data);
  dihe_assigned_bounds.swapTarget(&int_data);
  bond_evaluation_mask.swapTarget(&int_data);
  angl_evaluation_mask.swapTarget(&int_data);
  bond_modifiers.swapTarget(&char4_data);
  angl_modifiers.swapTarget(&char4_data);
  dihe_modifiers.swapTarget(&char4_data);
  bond_assigned_mods.swapTarget(&char4_data);
  angl_assigned_mods.swapTarget(&char4_data);
  dihe_assigned_mods.swapTarget(&char4_data);

  // Repair pointers relevant to virtual site placement
  virtual_site_atoms.swapTarget(&int_data);
  virtual_site_frame_types.swapTarget(&int_data);
  virtual_site_frame1_atoms.swapTarget(&int_data);
  virtual_site_frame2_atoms.swapTarget(&int_data);
  virtual_site_frame3_atoms.swapTarget(&int_data);
  virtual_site_frame4_atoms.swapTarget(&int_data);
  virtual_site_parameter_indices.swapTarget(&int_data);
  virtual_site_frame_dim1.swapTarget(&double_data);
  virtual_site_frame_dim2.swapTarget(&double_data);
  virtual_site_frame_dim3.swapTarget(&double_data);
  sp_virtual_site_frame_dim1.swapTarget(&float_data);
  sp_virtual_site_frame_dim2.swapTarget(&float_data);
  sp_virtual_site_frame_dim3.swapTarget(&float_data);

  // Repair pointers relevant to the non-bonded calculation
  charge_indices.swapTarget(&int_data);
  lennard_jones_indices.swapTarget(&int_data);
  atom_exclusion_bounds.swapTarget(&int_data);
  atom_exclusion_list.swapTarget(&int_data);
  nb11_exclusion_bounds.swapTarget(&int_data);
  nb11_exclusion_list.swapTarget(&int_data);
  nb12_exclusion_bounds.swapTarget(&int_data);
  nb12_exclusion_list.swapTarget(&int_data);
  nb13_exclusion_bounds.swapTarget(&int_data);
  nb13_exclusion_list.swapTarget(&int_data);
  nb14_exclusion_bounds.swapTarget(&int_data);
  nb14_exclusion_list.swapTarget(&int_data);
  infr14_i_atoms.swapTarget(&int_data);
  infr14_l_atoms.swapTarget(&int_data);
  infr14_parameter_indices.swapTarget(&int_data);
  neck_gb_indices.swapTarget(&int_data);
  charge_parameters.swapTarget(&double_data);
  lj_a_values.swapTarget(&double_data);
  lj_b_values.swapTarget(&double_data);
  lj_c_values.swapTarget(&double_data);
  lj_14_a_values.swapTarget(&double_data);
  lj_14_b_values.swapTarget(&double_data);
  lj_14_c_values.swapTarget(&double_data);
  lj_sigma_values.swapTarget(&double_data);
  lj_14_sigma_values.swapTarget(&double_data);
  lj_type_corrections.swapTarget(&double_data);
  attn14_elec_factors.swapTarget(&double_data);
  attn14_vdw_factors.swapTarget(&double_data);
  atomic_pb_radii.swapTarget(&double_data);
  gb_screening_factors.swapTarget(&double_data);
  gb_alpha_parameters.swapTarget(&double_data);
  gb_beta_parameters.swapTarget(&double_data);
  gb_gamma_parameters.swapTarget(&double_data);
  sp_charge_parameters.swapTarget(&float_data);
  sp_lj_a_values.swapTarget(&float_data);
  sp_lj_b_values.swapTarget(&float_data);
  sp_lj_c_values.swapTarget(&float_data);
  sp_lj_14_a_values.swapTarget(&float_data);
  sp_lj_14_b_values.swapTarget(&float_data);
  sp_lj_14_c_values.swapTarget(&float_data);
  sp_lj_sigma_values.swapTarget(&float_data);
  sp_lj_14_sigma_values.swapTarget(&float_data);
  sp_lj_type_corrections.swapTarget(&float_data);
  sp_attn14_elec_factors.swapTarget(&float_data);
  sp_attn14_vdw_factors.swapTarget(&float_data);
  sp_atomic_pb_radii.swapTarget(&float_data);
  sp_gb_screening_factors.swapTarget(&float_data);
  sp_gb_alpha_parameters.swapTarget(&float_data);
  sp_gb_beta_parameters.swapTarget(&float_data);
  sp_gb_gamma_parameters.swapTarget(&float_data);

  // Repair pointers relevant to the MD propagation algorithm (constraints)
  settle_oxygen_atoms.swapTarget(&int_data);
  settle_hydro1_atoms.swapTarget(&int_data);
  settle_hydro2_atoms.swapTarget(&int_data);
  settle_parameter_indices.swapTarget(&int_data);
  constraint_group_atoms.swapTarget(&int_data);
  constraint_group_bounds.swapTarget(&int_data);
  constraint_parameter_indices.swapTarget(&int_data);
  constraint_parameter_bounds.swapTarget(&int_data);
  settle_mo.swapTarget(&double_data);
  settle_mh.swapTarget(&double_data);
  settle_moh.swapTarget(&double_data);
  settle_mormt.swapTarget(&double_data);
  settle_mhrmt.swapTarget(&double_data);
  settle_ra.swapTarget(&double_data);
  settle_rb.swapTarget(&double_data);
  settle_rc.swapTarget(&double_data);
  constraint_inverse_masses.swapTarget(&double_data);
  constraint_squared_lengths.swapTarget(&double_data);
  sp_settle_mo.swapTarget(&float_data);
  sp_settle_mh.swapTarget(&float_data);
  sp_settle_moh.swapTarget(&float_data);
  sp_settle_mormt.swapTarget(&float_data);
  sp_settle_mhrmt.swapTarget(&float_data);
  sp_settle_ra.swapTarget(&float_data);
  sp_settle_rb.swapTarget(&float_data);
  sp_settle_rc.swapTarget(&float_data);
  sp_constraint_inverse_masses.swapTarget(&float_data);
  sp_constraint_squared_lengths.swapTarget(&float_data);

  // Repair overflow name pointers
  atom_overflow_names.swapTarget(&char4_data);
  atom_overflow_types.swapTarget(&char4_data);
  residue_overflow_names.swapTarget(&char4_data);

  // Repair information to unused pointers
  tree_joining_info.swapTarget(&int_data);
  last_rotator_info.swapTarget(&int_data);
  solty_info.swapTarget(&double_data);
  hbond_a_values.swapTarget(&double_data);
  hbond_b_values.swapTarget(&double_data);
  hbond_cutoffs.swapTarget(&double_data);
  atomic_polarizabilities.swapTarget(&double_data);
  tree_symbols.swapTarget(&char4_data);
}

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(const AtomGraph &original) :

    // General information
    version_stamp{},
    date{original.date},
    title{original.title},
    source{original.source},
    force_fields{original.force_fields},

    // Counts of atoms, residues, and other parts of the system
    atom_count{original.atom_count},
    residue_count{original.residue_count},
    molecule_count{original.molecule_count},
    largest_residue_size{original.largest_residue_size},
    last_solute_residue{original.last_solute_residue},
    last_solute_atom{original.last_solute_atom},
    first_solvent_molecule{original.first_solvent_molecule},
    last_atom_before_cap{original.last_atom_before_cap},
    implicit_copy_count{original.implicit_copy_count},
    largest_molecule_size{original.largest_molecule_size},
    water_residue_size{original.water_residue_size},
    water_residue_count{original.water_residue_count},
    unconstrained_dof{original.unconstrained_dof},
    constrained_dof{original.constrained_dof},
    descriptors{original.descriptors},
    residue_limits{original.residue_limits},
    atom_struc_numbers{original.atom_struc_numbers},
    residue_numbers{original.residue_numbers},
    molecule_limits{original.molecule_limits},

    // Atom and residue details
    atomic_numbers{original.atomic_numbers},
    mobile_atoms{original.mobile_atoms},
    molecule_membership{original.molecule_membership},
    molecule_contents{original.molecule_contents},
    atomic_charges{original.atomic_charges},
    atomic_masses{original.atomic_masses},
    inverse_atomic_masses{original.inverse_atomic_masses},
    sp_atomic_charges{original.sp_atomic_charges},
    sp_atomic_masses{original.sp_atomic_masses},
    sp_inverse_atomic_masses{original.sp_inverse_atomic_masses},
    atom_names{original.atom_names},
    atom_types{original.atom_types},
    residue_names{original.residue_names},

    // CHARMM force field family parameters
    urey_bradley_term_count{original.urey_bradley_term_count},
    charmm_impr_term_count{original.charmm_impr_term_count},
    cmap_term_count{original.cmap_term_count},
    urey_bradley_parameter_count{original.urey_bradley_parameter_count},
    charmm_impr_parameter_count{original.charmm_impr_parameter_count},
    cmap_surface_count{original.cmap_surface_count},
    urey_bradley_pert_term_count{original.urey_bradley_pert_term_count},
    charmm_impr_pert_term_count{original.charmm_impr_pert_term_count},
    cmap_pert_term_count{original.cmap_pert_term_count},
    urey_bradleys_in_perturbed_group{original.urey_bradleys_in_perturbed_group},
    charmm_imprs_in_perturbed_group{original.charmm_imprs_in_perturbed_group},
    cmaps_in_perturbed_group{original.cmaps_in_perturbed_group},
    urey_bradley_i_atoms{original.urey_bradley_i_atoms},
    urey_bradley_k_atoms{original.urey_bradley_k_atoms},
    urey_bradley_parameter_indices{original.urey_bradley_parameter_indices},
    urey_bradley_assigned_atoms{original.urey_bradley_assigned_atoms},
    urey_bradley_assigned_index{original.urey_bradley_assigned_index},
    urey_bradley_assigned_terms{original.urey_bradley_assigned_terms},
    urey_bradley_assigned_bounds{original.urey_bradley_assigned_bounds},
    charmm_impr_i_atoms{original.charmm_impr_i_atoms},
    charmm_impr_j_atoms{original.charmm_impr_j_atoms},
    charmm_impr_k_atoms{original.charmm_impr_k_atoms},
    charmm_impr_l_atoms{original.charmm_impr_l_atoms},
    charmm_impr_parameter_indices{original.charmm_impr_parameter_indices},
    charmm_impr_assigned_atoms{original.charmm_impr_assigned_atoms},
    charmm_impr_assigned_index{original.charmm_impr_assigned_index},
    charmm_impr_assigned_terms{original.charmm_impr_assigned_terms},
    charmm_impr_assigned_bounds{original.charmm_impr_assigned_bounds},
    cmap_i_atoms{original.cmap_i_atoms},
    cmap_j_atoms{original.cmap_j_atoms},
    cmap_k_atoms{original.cmap_k_atoms},
    cmap_l_atoms{original.cmap_l_atoms},
    cmap_m_atoms{original.cmap_m_atoms},
    cmap_surface_dimensions{original.cmap_surface_dimensions},
    cmap_surface_bounds{original.cmap_surface_bounds},
    cmap_patch_bounds{original.cmap_patch_bounds},
    cmap_surface_indices{original.cmap_surface_indices},
    cmap_assigned_atoms{original.cmap_assigned_atoms},
    cmap_assigned_index{original.cmap_assigned_index},
    cmap_assigned_terms{original.cmap_assigned_terms},
    cmap_assigned_bounds{original.cmap_assigned_bounds},
    urey_bradley_stiffnesses{original.urey_bradley_stiffnesses},
    urey_bradley_equilibria{original.urey_bradley_equilibria},
    charmm_impr_stiffnesses{original.charmm_impr_stiffnesses},
    charmm_impr_phase_angles{original.charmm_impr_phase_angles},
    cmap_surfaces{original.cmap_surfaces},
    cmap_phi_derivatives{original.cmap_phi_derivatives},
    cmap_psi_derivatives{original.cmap_psi_derivatives},
    cmap_phi_psi_derivatives{original.cmap_phi_psi_derivatives},
    cmap_patches{original.cmap_patches},
    sp_urey_bradley_stiffnesses{original.sp_urey_bradley_stiffnesses},
    sp_urey_bradley_equilibria{original.sp_urey_bradley_equilibria},
    sp_charmm_impr_stiffnesses{original.sp_charmm_impr_stiffnesses},
    sp_charmm_impr_phase_angles{original.sp_charmm_impr_phase_angles},
    sp_cmap_surfaces{original.sp_cmap_surfaces},
    sp_cmap_phi_derivatives{original.sp_cmap_phi_derivatives},
    sp_cmap_psi_derivatives{original.sp_cmap_psi_derivatives},
    sp_cmap_phi_psi_derivatives{original.sp_cmap_phi_psi_derivatives},
    sp_cmap_patches{original.sp_cmap_patches},

    // Relevant information for the bonded calculation
    bond_term_with_hydrogen{original.bond_term_with_hydrogen},
    angl_term_with_hydrogen{original.angl_term_with_hydrogen},
    dihe_term_with_hydrogen{original.dihe_term_with_hydrogen},
    bond_term_without_hydrogen{original.bond_term_without_hydrogen},
    angl_term_without_hydrogen{original.angl_term_without_hydrogen},
    dihe_term_without_hydrogen{original.dihe_term_without_hydrogen},
    bond_term_count{original.bond_term_count},
    angl_term_count{original.angl_term_count},
    dihe_term_count{original.dihe_term_count},
    bond_parameter_count{original.bond_parameter_count},
    angl_parameter_count{original.angl_parameter_count},
    dihe_parameter_count{original.dihe_parameter_count},
    bond_perturbation_term_count{original.bond_perturbation_term_count},
    angl_perturbation_term_count{original.angl_perturbation_term_count},
    dihe_perturbation_term_count{original.dihe_perturbation_term_count},
    bonds_in_perturbed_group{original.bonds_in_perturbed_group},
    angls_in_perturbed_group{original.angls_in_perturbed_group},
    dihes_in_perturbed_group{original.dihes_in_perturbed_group},
    bonded_group_count{original.bonded_group_count},
    bond_stiffnesses{original.bond_stiffnesses},
    bond_equilibria{original.bond_equilibria},
    angl_stiffnesses{original.angl_stiffnesses},
    angl_equilibria{original.angl_equilibria},
    dihe_amplitudes{original.dihe_amplitudes},
    dihe_periodicities{original.dihe_periodicities},
    dihe_phase_angles{original.dihe_phase_angles},
    sp_bond_stiffnesses{original.sp_bond_stiffnesses},
    sp_bond_equilibria{original.sp_bond_equilibria},
    sp_angl_stiffnesses{original.sp_angl_stiffnesses},
    sp_angl_equilibria{original.sp_angl_equilibria},
    sp_dihe_amplitudes{original.sp_dihe_amplitudes},
    sp_dihe_periodicities{original.sp_dihe_periodicities},
    sp_dihe_phase_angles{original.sp_dihe_phase_angles},
    bond_i_atoms{original.bond_i_atoms},
    bond_j_atoms{original.bond_j_atoms},
    bond_parameter_indices{original.bond_parameter_indices},
    bond_assigned_atoms{original.bond_assigned_atoms},
    bond_assigned_index{original.bond_assigned_index},
    bond_assigned_terms{original.bond_assigned_terms},
    bond_assigned_bounds{original.bond_assigned_bounds},
    angl_i_atoms{original.angl_i_atoms},
    angl_j_atoms{original.angl_j_atoms},
    angl_k_atoms{original.angl_k_atoms},
    angl_parameter_indices{original.angl_parameter_indices},
    angl_assigned_atoms{original.angl_assigned_atoms},
    angl_assigned_index{original.angl_assigned_index},
    angl_assigned_terms{original.angl_assigned_terms},
    angl_assigned_bounds{original.angl_assigned_bounds},
    dihe_i_atoms{original.dihe_i_atoms},
    dihe_j_atoms{original.dihe_j_atoms},
    dihe_k_atoms{original.dihe_k_atoms},
    dihe_l_atoms{original.dihe_l_atoms},
    dihe_parameter_indices{original.dihe_parameter_indices},
    dihe14_parameter_indices{original.dihe14_parameter_indices},
    dihe_assigned_atoms{original.dihe_assigned_atoms},
    dihe_assigned_index{original.dihe_assigned_index},
    dihe_assigned_terms{original.dihe_assigned_terms},
    dihe_assigned_bounds{original.dihe_assigned_bounds},
    bond_evaluation_mask{original.bond_evaluation_mask},
    angl_evaluation_mask{original.angl_evaluation_mask},
    bond_modifiers{original.bond_modifiers},
    angl_modifiers{original.angl_modifiers},
    dihe_modifiers{original.dihe_modifiers},
    bond_assigned_mods{original.bond_assigned_mods},
    angl_assigned_mods{original.angl_assigned_mods},
    dihe_assigned_mods{original.dihe_assigned_mods},

    // Information relevant to virtual site placement
    virtual_site_count{original.virtual_site_count},
    virtual_site_parameter_set_count{original.virtual_site_parameter_set_count},
    virtual_site_atoms{original.virtual_site_atoms},
    virtual_site_frame_types{original.virtual_site_frame_types},
    virtual_site_frame1_atoms{original.virtual_site_frame1_atoms},
    virtual_site_frame2_atoms{original.virtual_site_frame2_atoms},
    virtual_site_frame3_atoms{original.virtual_site_frame3_atoms},
    virtual_site_frame4_atoms{original.virtual_site_frame4_atoms},
    virtual_site_parameter_indices{original.virtual_site_parameter_indices},
    virtual_site_frame_dim1{original.virtual_site_frame_dim1},
    virtual_site_frame_dim2{original.virtual_site_frame_dim2},
    virtual_site_frame_dim3{original.virtual_site_frame_dim3},
    sp_virtual_site_frame_dim1{original.sp_virtual_site_frame_dim1},
    sp_virtual_site_frame_dim2{original.sp_virtual_site_frame_dim2},
    sp_virtual_site_frame_dim3{original.sp_virtual_site_frame_dim3},

    // Relevant information for the non-bonded calculation
    charge_type_count{original.charge_type_count},
    lj_type_count{original.lj_type_count},
    has_14_lennard_jones_data{original.has_14_lennard_jones_data},
    total_exclusions{original.total_exclusions},
    attenuated_14_type_count{original.attenuated_14_type_count},
    inferred_14_attenuations{original.inferred_14_attenuations},
    periodic_box_class{original.periodic_box_class},
    gb_style{original.gb_style},
    dielectric_constant{original.dielectric_constant},
    salt_concentration{original.salt_concentration},
    coulomb_constant{original.coulomb_constant},
    pb_radii_set{original.pb_radii_set},
    charge_indices{original.charge_indices},
    lennard_jones_indices{original.lennard_jones_indices},
    atom_exclusion_bounds{original.atom_exclusion_bounds},
    atom_exclusion_list{original.atom_exclusion_list},
    nb11_exclusion_bounds{original.nb11_exclusion_bounds},
    nb11_exclusion_list{original.nb11_exclusion_list},
    nb12_exclusion_bounds{original.nb12_exclusion_bounds},
    nb12_exclusion_list{original.nb12_exclusion_list},
    nb13_exclusion_bounds{original.nb13_exclusion_bounds},
    nb13_exclusion_list{original.nb13_exclusion_list},
    nb14_exclusion_bounds{original.nb14_exclusion_bounds},
    nb14_exclusion_list{original.nb14_exclusion_list},
    infr14_i_atoms{original.infr14_i_atoms},
    infr14_l_atoms{original.infr14_l_atoms},
    infr14_parameter_indices{original.infr14_parameter_indices},
    neck_gb_indices{original.neck_gb_indices},
    charge_parameters{original.charge_parameters},
    lj_a_values{original.lj_a_values},
    lj_b_values{original.lj_b_values},
    lj_c_values{original.lj_c_values},
    lj_14_a_values{original.lj_14_a_values},
    lj_14_b_values{original.lj_14_b_values},
    lj_14_c_values{original.lj_14_c_values},
    lj_sigma_values{original.lj_sigma_values},
    lj_14_sigma_values{original.lj_14_sigma_values},
    lj_type_corrections{original.lj_type_corrections},
    attn14_elec_factors{original.attn14_elec_factors},
    attn14_vdw_factors{original.attn14_vdw_factors},
    atomic_pb_radii{original.atomic_pb_radii},
    gb_screening_factors{original.gb_screening_factors},
    gb_alpha_parameters{original.gb_alpha_parameters},
    gb_beta_parameters{original.gb_beta_parameters},
    gb_gamma_parameters{original.gb_gamma_parameters},
    sp_charge_parameters{original.sp_charge_parameters},
    sp_lj_a_values{original.sp_lj_a_values},
    sp_lj_b_values{original.sp_lj_b_values},
    sp_lj_c_values{original.sp_lj_c_values},
    sp_lj_14_a_values{original.sp_lj_14_a_values},
    sp_lj_14_b_values{original.sp_lj_14_b_values},
    sp_lj_14_c_values{original.sp_lj_14_c_values},
    sp_lj_sigma_values{original.sp_lj_sigma_values},
    sp_lj_14_sigma_values{original.sp_lj_14_sigma_values},
    sp_lj_type_corrections{original.sp_lj_type_corrections},
    sp_attn14_elec_factors{original.sp_attn14_elec_factors},
    sp_attn14_vdw_factors{original.sp_attn14_vdw_factors},
    sp_atomic_pb_radii{original.sp_atomic_pb_radii},
    sp_gb_screening_factors{original.sp_gb_screening_factors},
    sp_gb_alpha_parameters{original.sp_gb_alpha_parameters},
    sp_gb_beta_parameters{original.sp_gb_beta_parameters},
    sp_gb_gamma_parameters{original.sp_gb_gamma_parameters},

    // MD propagation algorithm directives (including constraints)
    use_shake{original.use_shake},
    use_settle{original.use_settle},
    use_perturbation_info{original.use_perturbation_info},
    use_solvent_cap_option{original.use_solvent_cap_option},
    use_polarization{original.use_polarization},
    water_residue_name{original.water_residue_name},
    bond_constraint_mask{original.bond_constraint_mask},
    bond_constraint_omit_mask{original.bond_constraint_omit_mask},
    bond_constraint_count{original.bond_constraint_count},
    nonrigid_particle_count{original.nonrigid_particle_count},
    settle_group_count{original.settle_group_count},
    settle_parameter_count{original.settle_parameter_count},
    constraint_group_count{original.constraint_group_count},
    constraint_parameter_count{original.constraint_parameter_count},
    settle_oxygen_atoms{original.settle_oxygen_atoms},
    settle_hydro1_atoms{original.settle_hydro1_atoms},
    settle_hydro2_atoms{original.settle_hydro2_atoms},
    settle_parameter_indices{original.settle_parameter_indices},
    constraint_group_atoms{original.constraint_group_atoms},
    constraint_group_bounds{original.constraint_group_bounds},
    constraint_parameter_indices{original.constraint_parameter_indices},
    constraint_parameter_bounds{original.constraint_parameter_bounds},
    settle_mo{original.settle_mo},
    settle_mh{original.settle_mh},
    settle_moh{original.settle_moh},
    settle_mormt{original.settle_mormt},
    settle_mhrmt{original.settle_mhrmt},
    settle_ra{original.settle_ra},
    settle_rb{original.settle_rb},
    settle_rc{original.settle_rc},
    constraint_inverse_masses{original.constraint_inverse_masses},
    constraint_squared_lengths{original.constraint_squared_lengths},
    sp_settle_mo{original.sp_settle_mo},
    sp_settle_mh{original.sp_settle_mh},
    sp_settle_moh{original.sp_settle_moh},
    sp_settle_mormt{original.sp_settle_mormt},
    sp_settle_mhrmt{original.sp_settle_mhrmt},
    sp_settle_ra{original.sp_settle_ra},
    sp_settle_rb{original.sp_settle_rb},
    sp_settle_rc{original.sp_settle_rc},
    sp_constraint_inverse_masses{original.sp_constraint_inverse_masses},
    sp_constraint_squared_lengths{original.sp_constraint_squared_lengths},

    // Overflow name keys
    atom_overflow_names{original.atom_overflow_names},
    atom_overflow_types{original.atom_overflow_types},
    residue_overflow_names{original.residue_overflow_names},

    // Information currently unused
    unused_nhparm{original.unused_nhparm},
    unused_nparm{original.unused_nparm},
    unused_natyp{original.unused_natyp},
    hbond_10_12_parameter_count{original.hbond_10_12_parameter_count},
    heavy_bonds_plus_constraints{original.heavy_bonds_plus_constraints},
    heavy_angls_plus_constraints{original.heavy_angls_plus_constraints},
    heavy_dihes_plus_constraints{original.heavy_dihes_plus_constraints},
    tree_joining_info{original.tree_joining_info},
    last_rotator_info{original.last_rotator_info},
    solty_info{original.solty_info},
    hbond_a_values{original.hbond_a_values},
    hbond_b_values{original.hbond_b_values},
    hbond_cutoffs{original.hbond_cutoffs},
    atomic_polarizabilities{original.atomic_polarizabilities},
    tree_symbols{original.tree_symbols},

    // Copy Hybrid data
    int_data{original.int_data},
    double_data{original.double_data},
    float_data{original.float_data},
    char4_data{original.char4_data}
{
  const int nvst_char = strlen(original.version_stamp);
  for (int i = 0; i < nvst_char; i++) {
    version_stamp[i] = original.version_stamp[i];
  }

  // Repair pointers.  It is tedious to write the initializer list in the way that was done, and
  // it adds insult to injury to now have to list out all of the pointers yet again.  However, this
  // obviates the need for additional depth within the topology object to differentiate its trivial
  // and more complex member variables, so best to accept the tedium here for the sake of
  // simplicity elsewhere.
  rebasePointers();
}

//-------------------------------------------------------------------------------------------------
AtomGraph& AtomGraph::operator=(const AtomGraph &other) {
  
  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // Copy the elements of the original 
  const int nvst_char = strlen(other.version_stamp);
  for (int i = 0; i < nvst_char; i++) {
    version_stamp[i] = other.version_stamp[i];
  }
  date = other.date;
  title = other.title;
  source = other.source;
  force_fields = other.force_fields;

  // Copy counts of atoms, residues, and other parts of the system
  atom_count = other.atom_count;
  residue_count = other.residue_count;
  molecule_count = other.molecule_count;
  largest_residue_size = other.largest_residue_size;
  last_solute_residue = other.last_solute_residue;
  last_solute_atom = other.last_solute_atom;
  first_solvent_molecule = other.first_solvent_molecule;
  last_atom_before_cap = other.last_atom_before_cap;
  implicit_copy_count = other.implicit_copy_count;
  largest_molecule_size = other.largest_molecule_size;
  water_residue_size = other.water_residue_size;
  water_residue_count = other.water_residue_count;
  unconstrained_dof = other.unconstrained_dof;
  constrained_dof = other.constrained_dof;
  descriptors = other.descriptors;
  residue_limits = other.residue_limits;
  atom_struc_numbers = other.atom_struc_numbers;
  residue_numbers = other.residue_numbers;
  molecule_limits = other.molecule_limits;

  // Copy atom and residue details
  atomic_numbers = other.atomic_numbers;
  mobile_atoms = other.mobile_atoms;
  molecule_membership = other.molecule_membership;
  molecule_contents = other.molecule_contents;
  atomic_charges = other.atomic_charges;
  atomic_masses = other.atomic_masses;
  inverse_atomic_masses = other.inverse_atomic_masses;
  sp_atomic_charges = other.sp_atomic_charges;
  sp_atomic_masses = other.sp_atomic_masses;
  sp_inverse_atomic_masses = other.sp_inverse_atomic_masses;
  atom_names = other.atom_names;
  atom_types = other.atom_types;
  residue_names = other.residue_names;

  // Copy CHARMM force field family parameters
  urey_bradley_term_count = other.urey_bradley_term_count;
  charmm_impr_term_count = other.charmm_impr_term_count;
  cmap_term_count = other.cmap_term_count;
  urey_bradley_parameter_count = other.urey_bradley_parameter_count;
  charmm_impr_parameter_count = other.charmm_impr_parameter_count;
  cmap_surface_count = other.cmap_surface_count;
  urey_bradley_pert_term_count = other.urey_bradley_pert_term_count;
  charmm_impr_pert_term_count = other.charmm_impr_pert_term_count;
  cmap_pert_term_count = other.cmap_pert_term_count;
  urey_bradleys_in_perturbed_group = other.urey_bradleys_in_perturbed_group;
  charmm_imprs_in_perturbed_group = other.charmm_imprs_in_perturbed_group;
  cmaps_in_perturbed_group = other.cmaps_in_perturbed_group;
  urey_bradley_i_atoms = other.urey_bradley_i_atoms;
  urey_bradley_k_atoms = other.urey_bradley_k_atoms;
  urey_bradley_parameter_indices = other.urey_bradley_parameter_indices;
  urey_bradley_assigned_atoms = other.urey_bradley_assigned_atoms;
  urey_bradley_assigned_index = other.urey_bradley_assigned_index;
  urey_bradley_assigned_terms = other.urey_bradley_assigned_terms;
  urey_bradley_assigned_bounds = other.urey_bradley_assigned_bounds;
  charmm_impr_i_atoms = other.charmm_impr_i_atoms;
  charmm_impr_j_atoms = other.charmm_impr_j_atoms;
  charmm_impr_k_atoms = other.charmm_impr_k_atoms;
  charmm_impr_l_atoms = other.charmm_impr_l_atoms;
  charmm_impr_parameter_indices = other.charmm_impr_parameter_indices;
  charmm_impr_assigned_atoms = other.charmm_impr_assigned_atoms;
  charmm_impr_assigned_index = other.charmm_impr_assigned_index;
  charmm_impr_assigned_terms = other.charmm_impr_assigned_terms;
  charmm_impr_assigned_bounds = other.charmm_impr_assigned_bounds;
  cmap_i_atoms = other.cmap_i_atoms;
  cmap_j_atoms = other.cmap_j_atoms;
  cmap_k_atoms = other.cmap_k_atoms;
  cmap_l_atoms = other.cmap_l_atoms;
  cmap_m_atoms = other.cmap_m_atoms;
  cmap_surface_dimensions = other.cmap_surface_dimensions;
  cmap_surface_bounds = other.cmap_surface_bounds;
  cmap_patch_bounds = other.cmap_patch_bounds;
  cmap_surface_indices = other.cmap_surface_indices;
  cmap_assigned_atoms = other.cmap_assigned_atoms;
  cmap_assigned_index = other.cmap_assigned_index;
  cmap_assigned_terms = other.cmap_assigned_terms;
  cmap_assigned_bounds = other.cmap_assigned_bounds;
  urey_bradley_stiffnesses = other.urey_bradley_stiffnesses;
  urey_bradley_equilibria = other.urey_bradley_equilibria;
  charmm_impr_stiffnesses = other.charmm_impr_stiffnesses;
  charmm_impr_phase_angles = other.charmm_impr_phase_angles;
  cmap_surfaces = other.cmap_surfaces;
  cmap_phi_derivatives = other.cmap_phi_derivatives;
  cmap_psi_derivatives = other.cmap_psi_derivatives;
  cmap_phi_psi_derivatives = other.cmap_phi_psi_derivatives;
  cmap_patches = other.cmap_patches;
  sp_urey_bradley_stiffnesses = other.sp_urey_bradley_stiffnesses;
  sp_urey_bradley_equilibria = other.sp_urey_bradley_equilibria;
  sp_charmm_impr_stiffnesses = other.sp_charmm_impr_stiffnesses;
  sp_charmm_impr_phase_angles = other.sp_charmm_impr_phase_angles;
  sp_cmap_surfaces = other.sp_cmap_surfaces;
  sp_cmap_phi_derivatives = other.sp_cmap_phi_derivatives;
  sp_cmap_psi_derivatives = other.sp_cmap_psi_derivatives;
  sp_cmap_phi_psi_derivatives = other.sp_cmap_phi_psi_derivatives;
  sp_cmap_patches = other.sp_cmap_patches;

  // Copy relevant information for the bonded calculation
  bond_term_with_hydrogen = other.bond_term_with_hydrogen;
  angl_term_with_hydrogen = other.angl_term_with_hydrogen;
  dihe_term_with_hydrogen = other.dihe_term_with_hydrogen;
  bond_term_without_hydrogen = other.bond_term_without_hydrogen;
  angl_term_without_hydrogen = other.angl_term_without_hydrogen;
  dihe_term_without_hydrogen = other.dihe_term_without_hydrogen;
  bond_term_count = other.bond_term_count;
  angl_term_count = other.angl_term_count;
  dihe_term_count = other.dihe_term_count;
  bond_parameter_count = other.bond_parameter_count;
  angl_parameter_count = other.angl_parameter_count;
  dihe_parameter_count = other.dihe_parameter_count;
  bond_perturbation_term_count = other.bond_perturbation_term_count;
  angl_perturbation_term_count = other.angl_perturbation_term_count;
  dihe_perturbation_term_count = other.dihe_perturbation_term_count;
  bonds_in_perturbed_group = other.bonds_in_perturbed_group;
  angls_in_perturbed_group = other.angls_in_perturbed_group;
  dihes_in_perturbed_group = other.dihes_in_perturbed_group;
  bonded_group_count = other.bonded_group_count;
  bond_stiffnesses = other.bond_stiffnesses;
  bond_equilibria = other.bond_equilibria;
  angl_stiffnesses = other.angl_stiffnesses;
  angl_equilibria = other.angl_equilibria;
  dihe_amplitudes = other.dihe_amplitudes;
  dihe_periodicities = other.dihe_periodicities;
  dihe_phase_angles = other.dihe_phase_angles;
  sp_bond_stiffnesses = other.sp_bond_stiffnesses;
  sp_bond_equilibria = other.sp_bond_equilibria;
  sp_angl_stiffnesses = other.sp_angl_stiffnesses;
  sp_angl_equilibria = other.sp_angl_equilibria;
  sp_dihe_amplitudes = other.sp_dihe_amplitudes;
  sp_dihe_periodicities = other.sp_dihe_periodicities;
  sp_dihe_phase_angles = other.sp_dihe_phase_angles;
  bond_i_atoms = other.bond_i_atoms;
  bond_j_atoms = other.bond_j_atoms;
  bond_parameter_indices = other.bond_parameter_indices;
  bond_assigned_atoms = other.bond_assigned_atoms;
  bond_assigned_index = other.bond_assigned_index;
  bond_assigned_terms = other.bond_assigned_terms;
  bond_assigned_bounds = other.bond_assigned_bounds;
  angl_i_atoms = other.angl_i_atoms;
  angl_j_atoms = other.angl_j_atoms;
  angl_k_atoms = other.angl_k_atoms;
  angl_parameter_indices = other.angl_parameter_indices;
  angl_assigned_atoms = other.angl_assigned_atoms;
  angl_assigned_index = other.angl_assigned_index;
  angl_assigned_terms = other.angl_assigned_terms;
  angl_assigned_bounds = other.angl_assigned_bounds;
  dihe_i_atoms = other.dihe_i_atoms;
  dihe_j_atoms = other.dihe_j_atoms;
  dihe_k_atoms = other.dihe_k_atoms;
  dihe_l_atoms = other.dihe_l_atoms;
  dihe_parameter_indices = other.dihe_parameter_indices;
  dihe14_parameter_indices = other.dihe14_parameter_indices;
  dihe_assigned_atoms = other.dihe_assigned_atoms;
  dihe_assigned_index = other.dihe_assigned_index;
  dihe_assigned_terms = other.dihe_assigned_terms;
  dihe_assigned_bounds = other.dihe_assigned_bounds;
  bond_evaluation_mask = other.bond_evaluation_mask;
  angl_evaluation_mask = other.angl_evaluation_mask;
  bond_modifiers = other.bond_modifiers;
  angl_modifiers = other.angl_modifiers;
  dihe_modifiers = other.dihe_modifiers;
  bond_assigned_mods = other.bond_assigned_mods;
  angl_assigned_mods = other.angl_assigned_mods;
  dihe_assigned_mods = other.dihe_assigned_mods;

  // Copy information relevant to virtual site placement
  virtual_site_count = other.virtual_site_count;
  virtual_site_parameter_set_count = other.virtual_site_parameter_set_count;
  virtual_site_atoms = other.virtual_site_atoms;
  virtual_site_frame_types = other.virtual_site_frame_types;
  virtual_site_frame1_atoms = other.virtual_site_frame1_atoms;
  virtual_site_frame2_atoms = other.virtual_site_frame2_atoms;
  virtual_site_frame3_atoms = other.virtual_site_frame3_atoms;
  virtual_site_frame4_atoms = other.virtual_site_frame4_atoms;
  virtual_site_parameter_indices = other.virtual_site_parameter_indices;
  virtual_site_frame_dim1 = other.virtual_site_frame_dim1;
  virtual_site_frame_dim2 = other.virtual_site_frame_dim2;
  virtual_site_frame_dim3 = other.virtual_site_frame_dim3;
  sp_virtual_site_frame_dim1 = other.sp_virtual_site_frame_dim1;
  sp_virtual_site_frame_dim2 = other.sp_virtual_site_frame_dim2;
  sp_virtual_site_frame_dim3 = other.sp_virtual_site_frame_dim3;

  // Copy relevant information for the non-bonded calculation
  charge_type_count = other.charge_type_count;
  lj_type_count = other.lj_type_count;
  has_14_lennard_jones_data = other.has_14_lennard_jones_data;
  total_exclusions = other.total_exclusions;
  attenuated_14_type_count = other.attenuated_14_type_count;
  inferred_14_attenuations = other.inferred_14_attenuations;
  periodic_box_class = other.periodic_box_class;
  gb_style = other.gb_style;
  dielectric_constant = other.dielectric_constant;
  salt_concentration = other.salt_concentration;
  coulomb_constant = other.coulomb_constant;
  pb_radii_set = other.pb_radii_set;
  charge_indices = other.charge_indices;
  lennard_jones_indices = other.lennard_jones_indices;
  atom_exclusion_bounds = other.atom_exclusion_bounds;
  atom_exclusion_list = other.atom_exclusion_list;
  nb11_exclusion_bounds = other.nb11_exclusion_bounds;
  nb11_exclusion_list = other.nb11_exclusion_list;
  nb12_exclusion_bounds = other.nb12_exclusion_bounds;
  nb12_exclusion_list = other.nb12_exclusion_list;
  nb13_exclusion_bounds = other.nb13_exclusion_bounds;
  nb13_exclusion_list = other.nb13_exclusion_list;
  nb14_exclusion_bounds = other.nb14_exclusion_bounds;
  nb14_exclusion_list = other.nb14_exclusion_list;
  infr14_i_atoms = other.infr14_i_atoms;
  infr14_l_atoms = other.infr14_l_atoms;
  infr14_parameter_indices = other.infr14_parameter_indices;
  neck_gb_indices = other.neck_gb_indices;
  charge_parameters = other.charge_parameters;
  lj_a_values = other.lj_a_values;
  lj_b_values = other.lj_b_values;
  lj_c_values = other.lj_c_values;
  lj_14_a_values = other.lj_14_a_values;
  lj_14_b_values = other.lj_14_b_values;
  lj_14_c_values = other.lj_14_c_values;
  lj_sigma_values = other.lj_sigma_values;
  lj_14_sigma_values = other.lj_14_sigma_values;
  lj_type_corrections = other.lj_type_corrections;
  attn14_elec_factors = other.attn14_elec_factors;
  attn14_vdw_factors = other.attn14_vdw_factors;
  atomic_pb_radii = other.atomic_pb_radii;
  gb_screening_factors = other.gb_screening_factors;
  gb_alpha_parameters = other.gb_alpha_parameters;
  gb_beta_parameters = other.gb_beta_parameters;
  gb_gamma_parameters = other.gb_gamma_parameters;
  sp_charge_parameters = other.sp_charge_parameters;
  sp_lj_a_values = other.sp_lj_a_values;
  sp_lj_b_values = other.sp_lj_b_values;
  sp_lj_c_values = other.sp_lj_c_values;
  sp_lj_14_a_values = other.sp_lj_14_a_values;
  sp_lj_14_b_values = other.sp_lj_14_b_values;
  sp_lj_14_c_values = other.sp_lj_14_c_values;
  sp_lj_sigma_values = other.sp_lj_sigma_values;
  sp_lj_14_sigma_values = other.sp_lj_14_sigma_values;
  sp_lj_type_corrections = other.sp_lj_type_corrections;
  sp_attn14_elec_factors = other.sp_attn14_elec_factors;
  sp_attn14_vdw_factors = other.sp_attn14_vdw_factors;
  sp_atomic_pb_radii = other.sp_atomic_pb_radii;
  sp_gb_screening_factors = other.sp_gb_screening_factors;
  sp_gb_alpha_parameters = other.sp_gb_alpha_parameters;
  sp_gb_beta_parameters = other.sp_gb_beta_parameters;
  sp_gb_gamma_parameters = other.sp_gb_gamma_parameters;

  // Copy MD propagation algorithm directives
  use_shake = other.use_shake;
  use_settle = other.use_settle;
  use_perturbation_info = other.use_perturbation_info;
  use_solvent_cap_option = other.use_solvent_cap_option;
  use_polarization = other.use_polarization;
  water_residue_name = other.water_residue_name;
  bond_constraint_mask = other.bond_constraint_mask;
  bond_constraint_omit_mask = other.bond_constraint_omit_mask;
  bond_constraint_count = other.bond_constraint_count;
  nonrigid_particle_count = other.nonrigid_particle_count;
  settle_group_count = other.settle_group_count;
  settle_parameter_count = other.settle_parameter_count;
  constraint_group_count = other.constraint_group_count;
  constraint_parameter_count = other.constraint_parameter_count;
  settle_oxygen_atoms = other.settle_oxygen_atoms;
  settle_hydro1_atoms = other.settle_hydro1_atoms;
  settle_hydro2_atoms = other.settle_hydro2_atoms;
  settle_parameter_indices = other.settle_parameter_indices;
  constraint_group_atoms = other.constraint_group_atoms;
  constraint_group_bounds = other.constraint_group_bounds;
  constraint_parameter_indices = other.constraint_parameter_indices;
  constraint_parameter_bounds = other.constraint_parameter_bounds;
  settle_mo = other.settle_mo;
  settle_mh = other.settle_mh;
  settle_moh = other.settle_moh;
  settle_mormt = other.settle_mormt;
  settle_mhrmt = other.settle_mhrmt;
  settle_ra = other.settle_ra;
  settle_rb = other.settle_rb;
  settle_rc = other.settle_rc;
  constraint_inverse_masses = other.constraint_inverse_masses;
  constraint_squared_lengths = other.constraint_squared_lengths;
  sp_settle_mo = other.sp_settle_mo;
  sp_settle_mh = other.sp_settle_mh;
  sp_settle_moh = other.sp_settle_moh;
  sp_settle_mormt = other.sp_settle_mormt;
  sp_settle_mhrmt = other.sp_settle_mhrmt;
  sp_settle_ra = other.sp_settle_ra;
  sp_settle_rb = other.sp_settle_rb;
  sp_settle_rc = other.sp_settle_rc;
  sp_constraint_inverse_masses = other.sp_constraint_inverse_masses;
  sp_constraint_squared_lengths = other.sp_constraint_squared_lengths;

  // Copy overflow name keys
  atom_overflow_names = other.atom_overflow_names;
  atom_overflow_types = other.atom_overflow_types;
  residue_overflow_names = other.residue_overflow_names;

  // Copy information currently unused
  unused_nhparm = other.unused_nhparm;
  unused_nparm = other.unused_nparm;
  unused_natyp = other.unused_natyp;
  hbond_10_12_parameter_count = other.hbond_10_12_parameter_count;
  heavy_bonds_plus_constraints = other.heavy_bonds_plus_constraints;
  heavy_angls_plus_constraints = other.heavy_angls_plus_constraints;
  heavy_dihes_plus_constraints = other.heavy_dihes_plus_constraints;
  tree_joining_info = other.tree_joining_info;
  last_rotator_info = other.last_rotator_info;
  solty_info = other.solty_info;
  hbond_a_values = other.hbond_a_values;
  hbond_b_values = other.hbond_b_values;
  hbond_cutoffs = other.hbond_cutoffs;
  atomic_polarizabilities = other.atomic_polarizabilities;
  tree_symbols = other.tree_symbols;

  // Copy Hybrid data structures
  int_data = other.int_data;
  double_data = other.double_data;
  float_data = other.float_data;
  char4_data = other.char4_data;
  
  // Repair pointers and return the result, analogous to what happens in the PhaseSpace object
  rebasePointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
AtomGraph::AtomGraph(AtomGraph &&original) :

    // General information
    version_stamp{},
    date{original.date},
    title{std::move(original.title)},
    source{std::move(original.source)},
    force_fields{std::move(original.force_fields)},

    // Move or copy counts of atoms, residues, and other parts of the system
    atom_count{original.atom_count},
    residue_count{original.residue_count},
    molecule_count{original.molecule_count},
    largest_residue_size{original.largest_residue_size},
    last_solute_residue{original.last_solute_residue},
    last_solute_atom{original.last_solute_atom},
    first_solvent_molecule{original.first_solvent_molecule},
    last_atom_before_cap{original.last_atom_before_cap},
    implicit_copy_count{original.implicit_copy_count},
    largest_molecule_size{original.largest_molecule_size},
    water_residue_size{original.water_residue_size},
    water_residue_count{original.water_residue_count},
    unconstrained_dof{original.unconstrained_dof},
    constrained_dof{original.constrained_dof},
    descriptors{std::move(original.descriptors)},
    residue_limits{std::move(original.residue_limits)},
    atom_struc_numbers{std::move(original.atom_struc_numbers)},
    residue_numbers{std::move(original.residue_numbers)},
    molecule_limits{std::move(original.molecule_limits)},

    // Move atom and residue details
    atomic_numbers{std::move(original.atomic_numbers)},
    mobile_atoms{std::move(original.mobile_atoms)},
    molecule_membership{std::move(original.molecule_membership)},
    molecule_contents{std::move(original.molecule_contents)},
    atomic_charges{std::move(original.atomic_charges)},
    atomic_masses{std::move(original.atomic_masses)},
    inverse_atomic_masses{std::move(original.inverse_atomic_masses)},
    sp_atomic_charges{std::move(original.sp_atomic_charges)},
    sp_atomic_masses{std::move(original.sp_atomic_masses)},
    sp_inverse_atomic_masses{std::move(original.sp_inverse_atomic_masses)},
    atom_names{std::move(original.atom_names)},
    atom_types{std::move(original.atom_types)},
    residue_names{std::move(original.residue_names)},

    // Move CHARMM force field family parameters
    urey_bradley_term_count{original.urey_bradley_term_count},
    charmm_impr_term_count{original.charmm_impr_term_count},
    cmap_term_count{original.cmap_term_count},
    urey_bradley_parameter_count{original.urey_bradley_parameter_count},
    charmm_impr_parameter_count{original.charmm_impr_parameter_count},
    cmap_surface_count{original.cmap_surface_count},
    urey_bradley_pert_term_count{original.urey_bradley_pert_term_count},
    charmm_impr_pert_term_count{original.charmm_impr_pert_term_count},
    cmap_pert_term_count{original.cmap_pert_term_count},
    urey_bradleys_in_perturbed_group{original.urey_bradleys_in_perturbed_group},
    charmm_imprs_in_perturbed_group{original.charmm_imprs_in_perturbed_group},
    cmaps_in_perturbed_group{original.cmaps_in_perturbed_group},
    urey_bradley_i_atoms{std::move(original.urey_bradley_i_atoms)},
    urey_bradley_k_atoms{std::move(original.urey_bradley_k_atoms)},
    urey_bradley_parameter_indices{std::move(original.urey_bradley_parameter_indices)},
    urey_bradley_assigned_atoms{std::move(original.urey_bradley_assigned_atoms)},
    urey_bradley_assigned_index{std::move(original.urey_bradley_assigned_index)},
    urey_bradley_assigned_terms{std::move(original.urey_bradley_assigned_terms)},
    urey_bradley_assigned_bounds{std::move(original.urey_bradley_assigned_bounds)},
    charmm_impr_i_atoms{std::move(original.charmm_impr_i_atoms)},
    charmm_impr_j_atoms{std::move(original.charmm_impr_j_atoms)},
    charmm_impr_k_atoms{std::move(original.charmm_impr_k_atoms)},
    charmm_impr_l_atoms{std::move(original.charmm_impr_l_atoms)},
    charmm_impr_parameter_indices{std::move(original.charmm_impr_parameter_indices)},
    charmm_impr_assigned_atoms{std::move(original.charmm_impr_assigned_atoms)},
    charmm_impr_assigned_index{std::move(original.charmm_impr_assigned_index)},
    charmm_impr_assigned_terms{std::move(original.charmm_impr_assigned_terms)},
    charmm_impr_assigned_bounds{std::move(original.charmm_impr_assigned_bounds)},
    cmap_i_atoms{std::move(original.cmap_i_atoms)},
    cmap_j_atoms{std::move(original.cmap_j_atoms)},
    cmap_k_atoms{std::move(original.cmap_k_atoms)},
    cmap_l_atoms{std::move(original.cmap_l_atoms)},
    cmap_m_atoms{std::move(original.cmap_m_atoms)},
    cmap_surface_dimensions{std::move(original.cmap_surface_dimensions)},
    cmap_surface_bounds{std::move(original.cmap_surface_bounds)},
    cmap_patch_bounds{std::move(original.cmap_patch_bounds)},
    cmap_surface_indices{std::move(original.cmap_surface_indices)},
    cmap_assigned_atoms{std::move(original.cmap_assigned_atoms)},
    cmap_assigned_index{std::move(original.cmap_assigned_index)},
    cmap_assigned_terms{std::move(original.cmap_assigned_terms)},
    cmap_assigned_bounds{std::move(original.cmap_assigned_bounds)},
    urey_bradley_stiffnesses{std::move(original.urey_bradley_stiffnesses)},
    urey_bradley_equilibria{std::move(original.urey_bradley_equilibria)},
    charmm_impr_stiffnesses{std::move(original.charmm_impr_stiffnesses)},
    charmm_impr_phase_angles{std::move(original.charmm_impr_phase_angles)},
    cmap_surfaces{std::move(original.cmap_surfaces)},
    cmap_phi_derivatives{std::move(original.cmap_phi_derivatives)},
    cmap_psi_derivatives{std::move(original.cmap_psi_derivatives)},
    cmap_phi_psi_derivatives{std::move(original.cmap_phi_psi_derivatives)},
    cmap_patches{std::move(original.cmap_patches)},
    sp_urey_bradley_stiffnesses{std::move(original.sp_urey_bradley_stiffnesses)},
    sp_urey_bradley_equilibria{std::move(original.sp_urey_bradley_equilibria)},
    sp_charmm_impr_stiffnesses{std::move(original.sp_charmm_impr_stiffnesses)},
    sp_charmm_impr_phase_angles{std::move(original.sp_charmm_impr_phase_angles)},
    sp_cmap_surfaces{std::move(original.sp_cmap_surfaces)},
    sp_cmap_phi_derivatives{std::move(original.sp_cmap_phi_derivatives)},
    sp_cmap_psi_derivatives{std::move(original.sp_cmap_psi_derivatives)},
    sp_cmap_phi_psi_derivatives{std::move(original.sp_cmap_phi_psi_derivatives)},
    sp_cmap_patches{std::move(original.sp_cmap_patches)},

    // Move information relevant to the basic force field valence terms
    bond_term_with_hydrogen{original.bond_term_with_hydrogen},
    angl_term_with_hydrogen{original.angl_term_with_hydrogen},
    dihe_term_with_hydrogen{original.dihe_term_with_hydrogen},
    bond_term_without_hydrogen{original.bond_term_without_hydrogen},
    angl_term_without_hydrogen{original.angl_term_without_hydrogen},
    dihe_term_without_hydrogen{original.dihe_term_without_hydrogen},
    bond_term_count{original.bond_term_count},
    angl_term_count{original.angl_term_count},
    dihe_term_count{original.dihe_term_count},
    bond_parameter_count{original.bond_parameter_count},
    angl_parameter_count{original.angl_parameter_count},
    dihe_parameter_count{original.dihe_parameter_count},
    bond_perturbation_term_count{original.bond_perturbation_term_count},
    angl_perturbation_term_count{original.angl_perturbation_term_count},
    dihe_perturbation_term_count{original.dihe_perturbation_term_count},
    bonds_in_perturbed_group{original.bonds_in_perturbed_group},
    angls_in_perturbed_group{original.angls_in_perturbed_group},
    dihes_in_perturbed_group{original.dihes_in_perturbed_group},
    bonded_group_count{original.bonded_group_count},
    bond_stiffnesses{std::move(original.bond_stiffnesses)},
    bond_equilibria{std::move(original.bond_equilibria)},
    angl_stiffnesses{std::move(original.angl_stiffnesses)},
    angl_equilibria{std::move(original.angl_equilibria)},
    dihe_amplitudes{std::move(original.dihe_amplitudes)},
    dihe_periodicities{std::move(original.dihe_periodicities)},
    dihe_phase_angles{std::move(original.dihe_phase_angles)},
    sp_bond_stiffnesses{std::move(original.sp_bond_stiffnesses)},
    sp_bond_equilibria{std::move(original.sp_bond_equilibria)},
    sp_angl_stiffnesses{std::move(original.sp_angl_stiffnesses)},
    sp_angl_equilibria{std::move(original.sp_angl_equilibria)},
    sp_dihe_amplitudes{std::move(original.sp_dihe_amplitudes)},
    sp_dihe_periodicities{std::move(original.sp_dihe_periodicities)},
    sp_dihe_phase_angles{std::move(original.sp_dihe_phase_angles)},
    bond_i_atoms{std::move(original.bond_i_atoms)},
    bond_j_atoms{std::move(original.bond_j_atoms)},
    bond_parameter_indices{std::move(original.bond_parameter_indices)},
    bond_assigned_atoms{std::move(original.bond_assigned_atoms)},
    bond_assigned_index{std::move(original.bond_assigned_index)},
    bond_assigned_terms{std::move(original.bond_assigned_terms)},
    bond_assigned_bounds{std::move(original.bond_assigned_bounds)},
    angl_i_atoms{std::move(original.angl_i_atoms)},
    angl_j_atoms{std::move(original.angl_j_atoms)},
    angl_k_atoms{std::move(original.angl_k_atoms)},
    angl_parameter_indices{std::move(original.angl_parameter_indices)},
    angl_assigned_atoms{std::move(original.angl_assigned_atoms)},
    angl_assigned_index{std::move(original.angl_assigned_index)},
    angl_assigned_terms{std::move(original.angl_assigned_terms)},
    angl_assigned_bounds{std::move(original.angl_assigned_bounds)},
    dihe_i_atoms{std::move(original.dihe_i_atoms)},
    dihe_j_atoms{std::move(original.dihe_j_atoms)},
    dihe_k_atoms{std::move(original.dihe_k_atoms)},
    dihe_l_atoms{std::move(original.dihe_l_atoms)},
    dihe_parameter_indices{std::move(original.dihe_parameter_indices)},
    dihe14_parameter_indices{std::move(original.dihe14_parameter_indices)},
    dihe_assigned_atoms{std::move(original.dihe_assigned_atoms)},
    dihe_assigned_index{std::move(original.dihe_assigned_index)},
    dihe_assigned_terms{std::move(original.dihe_assigned_terms)},
    dihe_assigned_bounds{std::move(original.dihe_assigned_bounds)},
    bond_evaluation_mask{std::move(original.bond_evaluation_mask)},
    angl_evaluation_mask{std::move(original.angl_evaluation_mask)},
    bond_modifiers{std::move(original.bond_modifiers)},
    angl_modifiers{std::move(original.angl_modifiers)},
    dihe_modifiers{std::move(original.dihe_modifiers)},
    bond_assigned_mods{std::move(original.bond_assigned_mods)},
    angl_assigned_mods{std::move(original.angl_assigned_mods)},
    dihe_assigned_mods{std::move(original.dihe_assigned_mods)},

    // Move information relevant to virtual site handling
    virtual_site_count{original.virtual_site_count},
    virtual_site_parameter_set_count{original.virtual_site_parameter_set_count},
    virtual_site_atoms{std::move(original.virtual_site_atoms)},
    virtual_site_frame_types{std::move(original.virtual_site_frame_types)},
    virtual_site_frame1_atoms{std::move(original.virtual_site_frame1_atoms)},
    virtual_site_frame2_atoms{std::move(original.virtual_site_frame2_atoms)},
    virtual_site_frame3_atoms{std::move(original.virtual_site_frame3_atoms)},
    virtual_site_frame4_atoms{std::move(original.virtual_site_frame4_atoms)},
    virtual_site_parameter_indices{std::move(original.virtual_site_parameter_indices)},
    virtual_site_frame_dim1{std::move(original.virtual_site_frame_dim1)},
    virtual_site_frame_dim2{std::move(original.virtual_site_frame_dim2)},
    virtual_site_frame_dim3{std::move(original.virtual_site_frame_dim3)},
    sp_virtual_site_frame_dim1{std::move(original.sp_virtual_site_frame_dim1)},
    sp_virtual_site_frame_dim2{std::move(original.sp_virtual_site_frame_dim2)},
    sp_virtual_site_frame_dim3{std::move(original.sp_virtual_site_frame_dim3)},

    // Move information relevant to the non-bonded calculation
    charge_type_count{original.charge_type_count},
    lj_type_count{original.lj_type_count},
    has_14_lennard_jones_data{original.has_14_lennard_jones_data},
    total_exclusions{original.total_exclusions},
    attenuated_14_type_count{original.attenuated_14_type_count},
    inferred_14_attenuations{original.inferred_14_attenuations},
    periodic_box_class{original.periodic_box_class},
    gb_style{original.gb_style},
    dielectric_constant{original.dielectric_constant},
    salt_concentration{original.salt_concentration},
    coulomb_constant{original.coulomb_constant},
    pb_radii_set{std::move(original.pb_radii_set)},
    charge_indices{std::move(original.charge_indices)},
    lennard_jones_indices{std::move(original.lennard_jones_indices)},
    atom_exclusion_bounds{std::move(original.atom_exclusion_bounds)},
    atom_exclusion_list{std::move(original.atom_exclusion_list)},
    nb11_exclusion_bounds{std::move(original.nb11_exclusion_bounds)},
    nb11_exclusion_list{std::move(original.nb11_exclusion_list)},
    nb12_exclusion_bounds{std::move(original.nb12_exclusion_bounds)},
    nb12_exclusion_list{std::move(original.nb12_exclusion_list)},
    nb13_exclusion_bounds{std::move(original.nb13_exclusion_bounds)},
    nb13_exclusion_list{std::move(original.nb13_exclusion_list)},
    nb14_exclusion_bounds{std::move(original.nb14_exclusion_bounds)},
    nb14_exclusion_list{std::move(original.nb14_exclusion_list)},
    infr14_i_atoms{std::move(original.infr14_i_atoms)},
    infr14_l_atoms{std::move(original.infr14_l_atoms)},
    infr14_parameter_indices{std::move(original.infr14_parameter_indices)},
    neck_gb_indices{std::move(original.neck_gb_indices)},
    charge_parameters{std::move(original.charge_parameters)},
    lj_a_values{std::move(original.lj_a_values)},
    lj_b_values{std::move(original.lj_b_values)},
    lj_c_values{std::move(original.lj_c_values)},
    lj_14_a_values{std::move(original.lj_14_a_values)},
    lj_14_b_values{std::move(original.lj_14_b_values)},
    lj_14_c_values{std::move(original.lj_14_c_values)},
    lj_sigma_values{std::move(original.lj_sigma_values)},
    lj_14_sigma_values{std::move(original.lj_14_sigma_values)},
    lj_type_corrections{std::move(original.lj_type_corrections)},
    attn14_elec_factors{std::move(original.attn14_elec_factors)},
    attn14_vdw_factors{std::move(original.attn14_vdw_factors)},
    atomic_pb_radii{std::move(original.atomic_pb_radii)},
    gb_screening_factors{std::move(original.gb_screening_factors)},
    gb_alpha_parameters{std::move(original.gb_alpha_parameters)},
    gb_beta_parameters{std::move(original.gb_beta_parameters)},
    gb_gamma_parameters{std::move(original.gb_gamma_parameters)},
    sp_charge_parameters{std::move(original.sp_charge_parameters)},
    sp_lj_a_values{std::move(original.sp_lj_a_values)},
    sp_lj_b_values{std::move(original.sp_lj_b_values)},
    sp_lj_c_values{std::move(original.sp_lj_c_values)},
    sp_lj_14_a_values{std::move(original.sp_lj_14_a_values)},
    sp_lj_14_b_values{std::move(original.sp_lj_14_b_values)},
    sp_lj_14_c_values{std::move(original.sp_lj_14_c_values)},
    sp_lj_sigma_values{std::move(original.sp_lj_sigma_values)},
    sp_lj_14_sigma_values{std::move(original.sp_lj_14_sigma_values)},
    sp_lj_type_corrections{std::move(original.sp_lj_type_corrections)},
    sp_attn14_elec_factors{std::move(original.sp_attn14_elec_factors)},
    sp_attn14_vdw_factors{std::move(original.sp_attn14_vdw_factors)},
    sp_atomic_pb_radii{std::move(original.sp_atomic_pb_radii)},
    sp_gb_screening_factors{std::move(original.sp_gb_screening_factors)},
    sp_gb_alpha_parameters{std::move(original.sp_gb_alpha_parameters)},
    sp_gb_beta_parameters{std::move(original.sp_gb_beta_parameters)},
    sp_gb_gamma_parameters{std::move(original.sp_gb_gamma_parameters)},

    // Move information relevant to constraints and propagation of the MD algorithm
    use_shake{original.use_shake},
    use_settle{original.use_settle},
    use_perturbation_info{original.use_perturbation_info},
    use_solvent_cap_option{original.use_solvent_cap_option},
    use_polarization{original.use_polarization},
    water_residue_name{original.water_residue_name},
    bond_constraint_mask{std::move(original.bond_constraint_mask)},
    bond_constraint_omit_mask{std::move(original.bond_constraint_omit_mask)},
    bond_constraint_count{original.bond_constraint_count},
    nonrigid_particle_count{original.nonrigid_particle_count},
    settle_group_count{original.settle_group_count},
    settle_parameter_count{original.settle_parameter_count},
    constraint_group_count{original.constraint_group_count},
    constraint_parameter_count{original.constraint_parameter_count},
    settle_oxygen_atoms{std::move(original.settle_oxygen_atoms)},
    settle_hydro1_atoms{std::move(original.settle_hydro1_atoms)},
    settle_hydro2_atoms{std::move(original.settle_hydro2_atoms)},
    settle_parameter_indices{std::move(original.settle_parameter_indices)},
    constraint_group_atoms{std::move(original.constraint_group_atoms)},
    constraint_group_bounds{std::move(original.constraint_group_bounds)},
    constraint_parameter_indices{std::move(original.constraint_parameter_indices)},
    constraint_parameter_bounds{std::move(original.constraint_parameter_bounds)},
    settle_mo{std::move(original.settle_mo)},
    settle_mh{std::move(original.settle_mh)},
    settle_moh{std::move(original.settle_moh)},
    settle_mormt{std::move(original.settle_mormt)},
    settle_mhrmt{std::move(original.settle_mhrmt)},
    settle_ra{std::move(original.settle_ra)},
    settle_rb{std::move(original.settle_rb)},
    settle_rc{std::move(original.settle_rc)},
    constraint_inverse_masses{std::move(original.constraint_inverse_masses)},
    constraint_squared_lengths{std::move(original.constraint_squared_lengths)},
    sp_settle_mo{std::move(original.sp_settle_mo)},
    sp_settle_mh{std::move(original.sp_settle_mh)},
    sp_settle_moh{std::move(original.sp_settle_moh)},
    sp_settle_mormt{std::move(original.sp_settle_mormt)},
    sp_settle_mhrmt{std::move(original.sp_settle_mhrmt)},
    sp_settle_ra{std::move(original.sp_settle_ra)},
    sp_settle_rb{std::move(original.sp_settle_rb)},
    sp_settle_rc{std::move(original.sp_settle_rc)},
    sp_constraint_inverse_masses{std::move(original.sp_constraint_inverse_masses)},
    sp_constraint_squared_lengths{std::move(original.sp_constraint_squared_lengths)},

    // Move overflow name keys
    atom_overflow_names{std::move(original.atom_overflow_names)},
    atom_overflow_types{std::move(original.atom_overflow_types)},
    residue_overflow_names{std::move(original.residue_overflow_names)},

    // Move unused information
    unused_nhparm{original.unused_nhparm},
    unused_nparm{original.unused_nparm},
    unused_natyp{original.unused_natyp},
    hbond_10_12_parameter_count{std::move(original.hbond_10_12_parameter_count)},
    heavy_bonds_plus_constraints{std::move(original.heavy_bonds_plus_constraints)},
    heavy_angls_plus_constraints{std::move(original.heavy_angls_plus_constraints)},
    heavy_dihes_plus_constraints{std::move(original.heavy_dihes_plus_constraints)},
    tree_joining_info{std::move(original.tree_joining_info)},
    last_rotator_info{std::move(original.last_rotator_info)},
    solty_info{std::move(original.solty_info)},
    hbond_a_values{std::move(original.hbond_a_values)},
    hbond_b_values{std::move(original.hbond_b_values)},
    hbond_cutoffs{std::move(original.hbond_cutoffs)},
    atomic_polarizabilities{std::move(original.atomic_polarizabilities)},
    tree_symbols{std::move(original.tree_symbols)},

    // Move the Hybrid data storage arrays
    int_data{std::move(original.int_data)},
    double_data{std::move(original.double_data)},
    float_data{std::move(original.float_data)},
    char4_data{std::move(original.char4_data)}
{
  // Copy the last item, the explicit character version stamp
  const int nvst_char = strlen(original.version_stamp);
  for (int i = 0; i < nvst_char; i++) {
    version_stamp[i] = original.version_stamp[i];
  }
}

//-------------------------------------------------------------------------------------------------
AtomGraph& AtomGraph::operator=(AtomGraph &&other) {
  
  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // Copy or move the general elements of the original, as appropriate
  const int nvst_char = strlen(other.version_stamp);
  for (int i = 0; i < nvst_char; i++) {
    version_stamp[i] = other.version_stamp[i];
  }
  date = other.date;
  title = std::move(other.title);
  source = std::move(other.source);
  force_fields = std::move(other.force_fields);

  // Copy or move counts of atoms, residues, and other parts of the system
  atom_count = other.atom_count;
  residue_count = other.residue_count;
  molecule_count = other.molecule_count;
  largest_residue_size = other.largest_residue_size;
  last_solute_residue = other.last_solute_residue;
  last_solute_atom = other.last_solute_atom;
  first_solvent_molecule = other.first_solvent_molecule;
  last_atom_before_cap = other.last_atom_before_cap;
  implicit_copy_count = other.implicit_copy_count;
  largest_molecule_size = other.largest_molecule_size;
  water_residue_size = other.water_residue_size;
  water_residue_count = other.water_residue_count;
  unconstrained_dof = other.unconstrained_dof;
  constrained_dof = other.constrained_dof;
  descriptors = std::move(other.descriptors);
  residue_limits = std::move(other.residue_limits);
  atom_struc_numbers = std::move(other.atom_struc_numbers);
  residue_numbers = std::move(other.residue_numbers);
  molecule_limits = std::move(other.molecule_limits);

  // Copy or move atom and residue details
  atomic_numbers = std::move(other.atomic_numbers);
  mobile_atoms = std::move(other.mobile_atoms);
  molecule_membership = std::move(other.molecule_membership);
  molecule_contents = std::move(other.molecule_contents);
  atomic_charges = std::move(other.atomic_charges);
  atomic_masses = std::move(other.atomic_masses);
  inverse_atomic_masses = std::move(other.inverse_atomic_masses);
  sp_atomic_charges = std::move(other.sp_atomic_charges);
  sp_atomic_masses = std::move(other.sp_atomic_masses);
  sp_inverse_atomic_masses = std::move(other.sp_inverse_atomic_masses);
  atom_names = std::move(other.atom_names);
  atom_types = std::move(other.atom_types);
  residue_names = std::move(other.residue_names);

  // Copy or move CHARMM force field family parameters
  urey_bradley_term_count = other.urey_bradley_term_count;
  charmm_impr_term_count = other.charmm_impr_term_count;
  cmap_term_count = other.cmap_term_count;
  urey_bradley_parameter_count = other.urey_bradley_parameter_count;
  charmm_impr_parameter_count = other.charmm_impr_parameter_count;
  cmap_surface_count = other.cmap_surface_count;
  urey_bradley_pert_term_count = other.urey_bradley_pert_term_count;
  charmm_impr_pert_term_count = other.charmm_impr_pert_term_count;
  cmap_pert_term_count = other.cmap_pert_term_count;
  urey_bradleys_in_perturbed_group = other.urey_bradleys_in_perturbed_group;
  charmm_imprs_in_perturbed_group = other.charmm_imprs_in_perturbed_group;
  cmaps_in_perturbed_group = other.cmaps_in_perturbed_group;
  urey_bradley_i_atoms = std::move(other.urey_bradley_i_atoms);
  urey_bradley_k_atoms = std::move(other.urey_bradley_k_atoms);
  urey_bradley_parameter_indices = std::move(other.urey_bradley_parameter_indices);
  urey_bradley_assigned_atoms = std::move(other.urey_bradley_assigned_atoms);
  urey_bradley_assigned_index = std::move(other.urey_bradley_assigned_index);
  urey_bradley_assigned_terms = std::move(other.urey_bradley_assigned_terms);
  urey_bradley_assigned_bounds = std::move(other.urey_bradley_assigned_bounds);
  charmm_impr_i_atoms = std::move(other.charmm_impr_i_atoms);
  charmm_impr_j_atoms = std::move(other.charmm_impr_j_atoms);
  charmm_impr_k_atoms = std::move(other.charmm_impr_k_atoms);
  charmm_impr_l_atoms = std::move(other.charmm_impr_l_atoms);
  charmm_impr_parameter_indices = std::move(other.charmm_impr_parameter_indices);
  charmm_impr_assigned_atoms = std::move(other.charmm_impr_assigned_atoms);
  charmm_impr_assigned_index = std::move(other.charmm_impr_assigned_index);
  charmm_impr_assigned_terms = std::move(other.charmm_impr_assigned_terms);
  charmm_impr_assigned_bounds = std::move(other.charmm_impr_assigned_bounds);
  cmap_i_atoms = std::move(other.cmap_i_atoms);
  cmap_j_atoms = std::move(other.cmap_j_atoms);
  cmap_k_atoms = std::move(other.cmap_k_atoms);
  cmap_l_atoms = std::move(other.cmap_l_atoms);
  cmap_m_atoms = std::move(other.cmap_m_atoms);
  cmap_surface_dimensions = std::move(other.cmap_surface_dimensions);
  cmap_surface_bounds = std::move(other.cmap_surface_bounds);
  cmap_patch_bounds = std::move(other.cmap_patch_bounds);
  cmap_surface_indices = std::move(other.cmap_surface_indices);
  cmap_assigned_atoms = std::move(other.cmap_assigned_atoms);
  cmap_assigned_index = std::move(other.cmap_assigned_index);
  cmap_assigned_terms = std::move(other.cmap_assigned_terms);
  cmap_assigned_bounds = std::move(other.cmap_assigned_bounds);
  urey_bradley_stiffnesses = std::move(other.urey_bradley_stiffnesses);
  urey_bradley_equilibria = std::move(other.urey_bradley_equilibria);
  charmm_impr_stiffnesses = std::move(other.charmm_impr_stiffnesses);
  charmm_impr_phase_angles = std::move(other.charmm_impr_phase_angles);
  cmap_surfaces = std::move(other.cmap_surfaces);
  cmap_phi_derivatives = std::move(other.cmap_phi_derivatives);
  cmap_psi_derivatives = std::move(other.cmap_psi_derivatives);
  cmap_phi_psi_derivatives = std::move(other.cmap_phi_psi_derivatives);
  cmap_patches = std::move(other.cmap_patches);
  sp_urey_bradley_stiffnesses = std::move(other.sp_urey_bradley_stiffnesses);
  sp_urey_bradley_equilibria = std::move(other.sp_urey_bradley_equilibria);
  sp_charmm_impr_stiffnesses = std::move(other.sp_charmm_impr_stiffnesses);
  sp_charmm_impr_phase_angles = std::move(other.sp_charmm_impr_phase_angles);
  sp_cmap_surfaces = std::move(other.sp_cmap_surfaces);
  sp_cmap_phi_derivatives = std::move(other.sp_cmap_phi_derivatives);
  sp_cmap_psi_derivatives = std::move(other.sp_cmap_psi_derivatives);
  sp_cmap_phi_psi_derivatives = std::move(other.sp_cmap_phi_psi_derivatives);
  sp_cmap_patches = std::move(other.sp_cmap_patches);

  // Copy or move information relevant to basic force field valence terms
  bond_term_with_hydrogen = other.bond_term_with_hydrogen;
  angl_term_with_hydrogen = other.angl_term_with_hydrogen;
  dihe_term_with_hydrogen = other.dihe_term_with_hydrogen;
  bond_term_without_hydrogen = other.bond_term_without_hydrogen;
  angl_term_without_hydrogen = other.angl_term_without_hydrogen;
  dihe_term_without_hydrogen = other.dihe_term_without_hydrogen;
  bond_term_count = other.bond_term_count;
  angl_term_count = other.angl_term_count;
  dihe_term_count = other.dihe_term_count;
  bond_parameter_count = other.bond_parameter_count;
  angl_parameter_count = other.angl_parameter_count;
  dihe_parameter_count = other.dihe_parameter_count;
  bond_perturbation_term_count = other.bond_perturbation_term_count;
  angl_perturbation_term_count = other.angl_perturbation_term_count;
  dihe_perturbation_term_count = other.dihe_perturbation_term_count;
  bonds_in_perturbed_group = other.bonds_in_perturbed_group;
  angls_in_perturbed_group = other.angls_in_perturbed_group;
  dihes_in_perturbed_group = other.dihes_in_perturbed_group;
  bonded_group_count = other.bonded_group_count;
  bond_stiffnesses = std::move(other.bond_stiffnesses);
  bond_equilibria = std::move(other.bond_equilibria);
  angl_stiffnesses = std::move(other.angl_stiffnesses);
  angl_equilibria = std::move(other.angl_equilibria);
  dihe_amplitudes = std::move(other.dihe_amplitudes);
  dihe_periodicities = std::move(other.dihe_periodicities);
  dihe_phase_angles = std::move(other.dihe_phase_angles);
  sp_bond_stiffnesses = std::move(other.sp_bond_stiffnesses);
  sp_bond_equilibria = std::move(other.sp_bond_equilibria);
  sp_angl_stiffnesses = std::move(other.sp_angl_stiffnesses);
  sp_angl_equilibria = std::move(other.sp_angl_equilibria);
  sp_dihe_amplitudes = std::move(other.sp_dihe_amplitudes);
  sp_dihe_periodicities = std::move(other.sp_dihe_periodicities);
  sp_dihe_phase_angles = std::move(other.sp_dihe_phase_angles);
  bond_i_atoms = std::move(other.bond_i_atoms);
  bond_j_atoms = std::move(other.bond_j_atoms);
  bond_parameter_indices = std::move(other.bond_parameter_indices);
  bond_assigned_atoms = std::move(other.bond_assigned_atoms);
  bond_assigned_index = std::move(other.bond_assigned_index);
  bond_assigned_terms = std::move(other.bond_assigned_terms);
  bond_assigned_bounds = std::move(other.bond_assigned_bounds);
  angl_i_atoms = std::move(other.angl_i_atoms);
  angl_j_atoms = std::move(other.angl_j_atoms);
  angl_k_atoms = std::move(other.angl_k_atoms);
  angl_parameter_indices = std::move(other.angl_parameter_indices);
  angl_assigned_atoms = std::move(other.angl_assigned_atoms);
  angl_assigned_index = std::move(other.angl_assigned_index);
  angl_assigned_terms = std::move(other.angl_assigned_terms);
  angl_assigned_bounds = std::move(other.angl_assigned_bounds);
  dihe_i_atoms = std::move(other.dihe_i_atoms);
  dihe_j_atoms = std::move(other.dihe_j_atoms);
  dihe_k_atoms = std::move(other.dihe_k_atoms);
  dihe_l_atoms = std::move(other.dihe_l_atoms);
  dihe_parameter_indices = std::move(other.dihe_parameter_indices);
  dihe14_parameter_indices = std::move(other.dihe14_parameter_indices);
  dihe_assigned_atoms = std::move(other.dihe_assigned_atoms);
  dihe_assigned_index = std::move(other.dihe_assigned_index);
  dihe_assigned_terms = std::move(other.dihe_assigned_terms);
  dihe_assigned_bounds = std::move(other.dihe_assigned_bounds);
  bond_evaluation_mask = std::move(other.bond_evaluation_mask);
  angl_evaluation_mask = std::move(other.angl_evaluation_mask);
  bond_modifiers = std::move(other.bond_modifiers);
  angl_modifiers = std::move(other.angl_modifiers);
  dihe_modifiers = std::move(other.dihe_modifiers);
  bond_assigned_mods = std::move(other.bond_assigned_mods);
  angl_assigned_mods = std::move(other.angl_assigned_mods);
  dihe_assigned_mods = std::move(other.dihe_assigned_mods);

  // Copy or move information relevant to virtual site handling
  virtual_site_count = other.virtual_site_count;
  virtual_site_parameter_set_count = other.virtual_site_parameter_set_count;
  virtual_site_atoms = std::move(other.virtual_site_atoms);
  virtual_site_frame_types = std::move(other.virtual_site_frame_types);
  virtual_site_frame1_atoms = std::move(other.virtual_site_frame1_atoms);
  virtual_site_frame2_atoms = std::move(other.virtual_site_frame2_atoms);
  virtual_site_frame3_atoms = std::move(other.virtual_site_frame3_atoms);
  virtual_site_frame4_atoms = std::move(other.virtual_site_frame4_atoms);
  virtual_site_parameter_indices = std::move(other.virtual_site_parameter_indices);
  virtual_site_frame_dim1 = std::move(other.virtual_site_frame_dim1);
  virtual_site_frame_dim2 = std::move(other.virtual_site_frame_dim2);
  virtual_site_frame_dim3 = std::move(other.virtual_site_frame_dim3);
  sp_virtual_site_frame_dim1 = std::move(other.sp_virtual_site_frame_dim1);
  sp_virtual_site_frame_dim2 = std::move(other.sp_virtual_site_frame_dim2);
  sp_virtual_site_frame_dim3 = std::move(other.sp_virtual_site_frame_dim3);

  // Copy or move information relevant to the non-bonded calculation
  charge_type_count = other.charge_type_count;
  lj_type_count = other.lj_type_count;
  has_14_lennard_jones_data = other.has_14_lennard_jones_data;
  total_exclusions = other.total_exclusions;
  attenuated_14_type_count = other.attenuated_14_type_count;
  inferred_14_attenuations = other.inferred_14_attenuations;
  periodic_box_class = other.periodic_box_class;
  gb_style = other.gb_style;
  dielectric_constant = other.dielectric_constant;
  salt_concentration = other.salt_concentration;
  coulomb_constant = other.coulomb_constant;
  pb_radii_set = std::move(other.pb_radii_set);
  charge_indices = std::move(other.charge_indices);
  lennard_jones_indices = std::move(other.lennard_jones_indices);
  atom_exclusion_bounds = std::move(other.atom_exclusion_bounds);
  atom_exclusion_list = std::move(other.atom_exclusion_list);
  nb11_exclusion_bounds = std::move(other.nb11_exclusion_bounds);
  nb11_exclusion_list = std::move(other.nb11_exclusion_list);
  nb12_exclusion_bounds = std::move(other.nb12_exclusion_bounds);
  nb12_exclusion_list = std::move(other.nb12_exclusion_list);
  nb13_exclusion_bounds = std::move(other.nb13_exclusion_bounds);
  nb13_exclusion_list = std::move(other.nb13_exclusion_list);
  nb14_exclusion_bounds = std::move(other.nb14_exclusion_bounds);
  nb14_exclusion_list = std::move(other.nb14_exclusion_list);
  infr14_i_atoms = std::move(other.infr14_i_atoms);
  infr14_l_atoms = std::move(other.infr14_l_atoms);
  infr14_parameter_indices = std::move(other.infr14_parameter_indices);
  neck_gb_indices = std::move(other.neck_gb_indices);
  charge_parameters = std::move(other.charge_parameters);
  lj_a_values = std::move(other.lj_a_values);
  lj_b_values = std::move(other.lj_b_values);
  lj_c_values = std::move(other.lj_c_values);
  lj_14_a_values = std::move(other.lj_14_a_values);
  lj_14_b_values = std::move(other.lj_14_b_values);
  lj_14_c_values = std::move(other.lj_14_c_values);
  lj_sigma_values = std::move(other.lj_sigma_values);
  lj_14_sigma_values = std::move(other.lj_14_sigma_values);
  lj_type_corrections = std::move(other.lj_type_corrections);
  attn14_elec_factors = std::move(other.attn14_elec_factors);
  attn14_vdw_factors = std::move(other.attn14_vdw_factors);
  atomic_pb_radii = std::move(other.atomic_pb_radii);
  gb_screening_factors = std::move(other.gb_screening_factors);
  gb_alpha_parameters = std::move(other.gb_alpha_parameters);
  gb_beta_parameters = std::move(other.gb_beta_parameters);
  gb_gamma_parameters = std::move(other.gb_gamma_parameters);
  sp_charge_parameters = std::move(other.sp_charge_parameters);
  sp_lj_a_values = std::move(other.sp_lj_a_values);
  sp_lj_b_values = std::move(other.sp_lj_b_values);
  sp_lj_c_values = std::move(other.sp_lj_c_values);
  sp_lj_14_a_values = std::move(other.sp_lj_14_a_values);
  sp_lj_14_b_values = std::move(other.sp_lj_14_b_values);
  sp_lj_14_c_values = std::move(other.sp_lj_14_c_values);
  sp_lj_sigma_values = std::move(other.sp_lj_sigma_values);
  sp_lj_14_sigma_values = std::move(other.sp_lj_14_sigma_values);
  sp_lj_type_corrections = std::move(other.sp_lj_type_corrections);
  sp_attn14_elec_factors = std::move(other.sp_attn14_elec_factors);
  sp_attn14_vdw_factors = std::move(other.sp_attn14_vdw_factors);
  sp_atomic_pb_radii = std::move(other.sp_atomic_pb_radii);
  sp_gb_screening_factors = std::move(other.sp_gb_screening_factors);
  sp_gb_alpha_parameters = std::move(other.sp_gb_alpha_parameters);
  sp_gb_beta_parameters = std::move(other.sp_gb_beta_parameters);
  sp_gb_gamma_parameters = std::move(other.sp_gb_gamma_parameters);

  // Copy or move information relevant to constraints and propagation of the MD algorithm
  use_shake = other.use_shake;
  use_settle = other.use_settle;
  use_perturbation_info = other.use_perturbation_info;
  use_solvent_cap_option = other.use_solvent_cap_option;
  use_polarization = other.use_polarization;
  water_residue_name = other.water_residue_name;
  bond_constraint_mask = std::move(other.bond_constraint_mask);
  bond_constraint_omit_mask = std::move(other.bond_constraint_omit_mask);
  bond_constraint_count = other.bond_constraint_count;
  nonrigid_particle_count = other.nonrigid_particle_count;
  settle_group_count = other.settle_group_count;
  settle_parameter_count = other.settle_parameter_count;
  constraint_group_count = other.constraint_group_count;
  constraint_parameter_count = other.constraint_parameter_count;
  settle_oxygen_atoms = std::move(other.settle_oxygen_atoms);
  settle_hydro1_atoms = std::move(other.settle_hydro1_atoms);
  settle_hydro2_atoms = std::move(other.settle_hydro2_atoms);
  settle_parameter_indices = std::move(other.settle_parameter_indices);
  constraint_group_atoms = std::move(other.constraint_group_atoms);
  constraint_group_bounds = std::move(other.constraint_group_bounds);
  constraint_parameter_indices = std::move(other.constraint_parameter_indices);
  constraint_parameter_bounds = std::move(other.constraint_parameter_bounds);
  settle_mo = std::move(other.settle_mo);
  settle_mh = std::move(other.settle_mh);
  settle_moh = std::move(other.settle_moh);
  settle_mormt = std::move(other.settle_mormt);
  settle_mhrmt = std::move(other.settle_mhrmt);
  settle_ra = std::move(other.settle_ra);
  settle_rb = std::move(other.settle_rb);
  settle_rc = std::move(other.settle_rc);
  constraint_inverse_masses = std::move(other.constraint_inverse_masses);
  constraint_squared_lengths = std::move(other.constraint_squared_lengths);
  sp_settle_mo = std::move(other.sp_settle_mo);
  sp_settle_mh = std::move(other.sp_settle_mh);
  sp_settle_moh = std::move(other.sp_settle_moh);
  sp_settle_mormt = std::move(other.sp_settle_mormt);
  sp_settle_mhrmt = std::move(other.sp_settle_mhrmt);
  sp_settle_ra = std::move(other.sp_settle_ra);
  sp_settle_rb = std::move(other.sp_settle_rb);
  sp_settle_rc = std::move(other.sp_settle_rc);
  sp_constraint_inverse_masses = std::move(other.sp_constraint_inverse_masses);
  sp_constraint_squared_lengths = std::move(other.sp_constraint_squared_lengths);

  // Move overflow name keys
  atom_overflow_names = std::move(other.atom_overflow_names);
  atom_overflow_types = std::move(other.atom_overflow_types);
  residue_overflow_names = std::move(other.residue_overflow_names);

  // Copy or move unused information
  unused_nhparm = other.unused_nhparm;
  unused_nparm = other.unused_nparm;
  unused_natyp = other.unused_natyp;
  hbond_10_12_parameter_count = other.hbond_10_12_parameter_count;
  heavy_bonds_plus_constraints = std::move(other.heavy_bonds_plus_constraints);
  heavy_angls_plus_constraints = std::move(other.heavy_angls_plus_constraints);
  heavy_dihes_plus_constraints = std::move(other.heavy_dihes_plus_constraints);
  tree_joining_info = std::move(other.tree_joining_info);
  last_rotator_info = std::move(other.last_rotator_info);
  solty_info = std::move(other.solty_info);
  hbond_a_values = std::move(other.hbond_a_values);
  hbond_b_values = std::move(other.hbond_b_values);
  hbond_cutoffs = std::move(other.hbond_cutoffs);
  atomic_polarizabilities = std::move(other.atomic_polarizabilities);
  tree_symbols = std::move(other.tree_symbols);

  // Move the Hybrid data storage arrays
  int_data = std::move(other.int_data);
  double_data = std::move(other.double_data);
  float_data = std::move(other.float_data);
  char4_data = std::move(other.char4_data);
  
  return *this;
}

} // namespace topology
} // namespace stormm
