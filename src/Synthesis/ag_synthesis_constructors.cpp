#include "copyright.h"
#include "atomgraph_synthesis.h"

namespace stormm {
namespace synthesis {

using card::HybridKind;
using topology::amber_ancient_bioq;

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis::AtomGraphSynthesis(const std::vector<AtomGraph*> &topologies_in,
                                       const std::vector<RestraintApparatus*> &restraints_in,
                                       const std::vector<int> &topology_indices_in,
                                       const std::vector<int> &restraint_indices_in,
                                       const ExceptionResponse policy_in, const GpuDetails &gpu,
                                       StopWatch *timer_in) :

    // Counts spanning all topologies
    policy{policy_in},
    topology_count{static_cast<int>(topologies_in.size())},
    restraint_network_count{static_cast<int>(restraints_in.size())},
    system_count{static_cast<int>(topology_indices_in.size())},
    total_atoms{0}, total_virtual_sites{0}, total_bond_terms{0}, total_angl_terms{0},
    total_dihe_terms{0}, total_ubrd_terms{0}, total_cimp_terms{0}, total_cmap_terms{0},
    total_lj_types{0}, total_charge_types{0}, total_bond_params{0}, total_angl_params{0},
    total_dihe_params{0}, total_ubrd_params{0}, total_cimp_params{0}, total_cmap_surfaces{0},
    total_attn14_params{0}, total_vste_params{0}, total_position_restraints{0},
    total_distance_restraints{0}, total_angle_restraints{0}, total_dihedral_restraints{0},

    // Descriptors spanning all systems
    periodic_box_class{UnitCellType::NONE}, gb_style{ImplicitSolventModel::NONE},
    dielectric_constant{1.0}, is_kappa{0.0}, salt_concentration{0.0},
    gb_offset{default_gb_radii_offset}, gb_neckscale{default_gb_neck_scale},
    gb_neckcut{default_gb_neck_cut}, coulomb_constant{amber_ancient_bioq},
    use_bond_constraints{ShakeSetting::OFF}, use_settle{SettleSetting::OFF},
    largest_constraint_group{0}, water_residue_name{' ', ' ', ' ', ' '},
    pb_radii_sets{std::vector<AtomicRadiusSet>(system_count, AtomicRadiusSet::NONE)},

    // Individual topologies and restraint networks (RestraintApparatus objects), each of them
    // describing one or more of the individual systems within the synthesis
    topologies{topologies_in},
    restraint_networks{restraints_in},
    restraint_dummies{},
    topology_indices{HybridKind::POINTER, "tpsyn_top_indices"},
    restraint_indices{HybridKind::POINTER, "tpsyn_rst_indices"},
    atom_counts{HybridKind::POINTER, "tpsyn_atom_counts"},
    residue_counts{HybridKind::POINTER, "tpsyn_res_counts"},
    molecule_counts{HybridKind::POINTER, "tpsyn_mol_counts"},
    largest_residue_sizes{HybridKind::POINTER, "tpsyn_max_res"},
    last_solute_residues{HybridKind::POINTER, "tpsyn_last_sol_res"},
    last_solute_atoms{HybridKind::POINTER, "tpsyn_last_sol_atm"},
    first_solvent_molecules{HybridKind::POINTER, "tpsyn_1st_solv_mol"},

    // Counts of energy terms and particle-related quantities
    ubrd_term_counts{HybridKind::POINTER, "tpsyn_ubrd_counts"},
    cimp_term_counts{HybridKind::POINTER, "tpsyn_cimp_counts"},
    cmap_term_counts{HybridKind::POINTER, "tpsyn_cmap_counts"},
    bond_term_counts{HybridKind::POINTER, "tpsyn_bond_counts"},
    angl_term_counts{HybridKind::POINTER, "tpsyn_angl_counts"},
    dihe_term_counts{HybridKind::POINTER, "tpsyn_dihe_counts"},
    virtual_site_counts{HybridKind::POINTER, "tpsyn_vsite_counts"},
    posn_restraint_counts{HybridKind::POINTER, "tpsyn_rposn_counts"},
    bond_restraint_counts{HybridKind::POINTER, "tpsyn_rbond_counts"},
    angl_restraint_counts{HybridKind::POINTER, "tpsyn_rangl_counts"},
    dihe_restraint_counts{HybridKind::POINTER, "tpsyn_rdihe_counts"},
    lj_type_counts{HybridKind::POINTER, "tpsyn_atype_counts"},
    total_exclusion_counts{HybridKind::POINTER, "tpsyn_excl_counts"},
    rigid_water_counts{HybridKind::POINTER, "tpsyn_rwat_counts"},
    bond_constraint_counts{HybridKind::POINTER, "tpsyn_bcnst_counts"},
    degrees_of_freedom{HybridKind::POINTER, "tpsyn_deg_freedom"},
    cnst_degrees_of_freedom{HybridKind::POINTER, "tpsyn_cdeg_freedom"},
    nonrigid_particle_counts{HybridKind::POINTER, "tpsyn_n_nonrigid"},

    // Offsets for each system's term lists
    atom_offsets{HybridKind::POINTER, "tpsyn_atom_offsets"},
    atom_bit_offsets{HybridKind::POINTER, "tpsyn_abit_offsets"},
    residue_offsets{HybridKind::POINTER, "tpsyn_res_offsets"},
    molecule_offsets{HybridKind::POINTER, "tpsyn_res_offsets"},
    ubrd_term_offsets{HybridKind::POINTER, "tpsyn_ubrd_offset"},
    cimp_term_offsets{HybridKind::POINTER, "tpsyn_cimp_offset"},
    cmap_term_offsets{HybridKind::POINTER, "tpsyn_cmap_offset"},
    bond_term_offsets{HybridKind::POINTER, "tpsyn_bond_offset"},
    angl_term_offsets{HybridKind::POINTER, "tpsyn_angl_offset"},
    dihe_term_offsets{HybridKind::POINTER, "tpsyn_dihe_offset"},
    virtual_site_offsets{HybridKind::POINTER, "tpsyn_vsite_offset"},
    posn_restraint_offsets{HybridKind::POINTER, "tpsyn_rposn_offset"},
    bond_restraint_offsets{HybridKind::POINTER, "tpsyn_rbond_offset"},
    angl_restraint_offsets{HybridKind::POINTER, "tpsyn_rangl_offset"},
    dihe_restraint_offsets{HybridKind::POINTER, "tpsyn_rdihe_offset"},
    sett_group_offsets{HybridKind::POINTER, "tpsyn_sett_offset"},
    cnst_group_offsets{HybridKind::POINTER, "tpsyn_cnst_offset"},
    nb_exclusion_offsets{HybridKind::POINTER, "tpsyn_nbexcl_offset"},
    lennard_jones_abc_offsets{HybridKind::POINTER, "tpsyn_ljtable_offset"},
    int_system_data{HybridKind::ARRAY, "tpsyn_int_data"},

    // Atom and residue details
    residue_limits{HybridKind::POINTER, "tpsyn_res_lims"},
    atom_struc_numbers{HybridKind::POINTER, "tpsyn_atom_struc_nums"},
    residue_numbers{HybridKind::POINTER, "tpsyn_res_numbers"},
    molecule_limits{HybridKind::POINTER, "tpsyn_mol_limits"},
    atomic_numbers{HybridKind::POINTER, "tpsyn_znum"},
    mobile_atoms{HybridKind::POINTER, "tpsyn_belly"},
    molecule_membership{HybridKind::POINTER, "tpsyn_molnum"},
    molecule_contents{HybridKind::POINTER, "tpsyn_mol_contents"},
    atomic_charges{HybridKind::POINTER, "tpsyn_atomq"},
    atomic_masses{HybridKind::POINTER, "tpsyn_mass"},
    inverse_atomic_masses{HybridKind::POINTER, "tpsyn_invmass"},
    sp_atomic_charges{HybridKind::POINTER, "tpsyn_atomq_sp"},
    sp_atomic_masses{HybridKind::POINTER, "tpsyn_mass_sp"},
    sp_inverse_atomic_masses{HybridKind::POINTER, "tpsyn_invmass_sp"},
    atom_names{HybridKind::POINTER, "tpsyn_atom_names"},
    atom_types{HybridKind::POINTER, "tpsyn_atom_types"},
    residue_names{HybridKind::POINTER, "tpsyn_res_names"},
    chem_int_data{HybridKind::ARRAY, "tpsyn_chem_ints"},
    chem_int2_data{HybridKind::ARRAY, "tpsyn_chem_int2s"},
    chem_double_data{HybridKind::ARRAY, "tpsyn_chem_doubles"},
    chem_float_data{HybridKind::ARRAY, "tpsyn_chem_floats"},
    chem_char4_data{HybridKind::ARRAY, "tpsyn_chem_char4s"},

    // Force field parameter maps from individual systems to the consensus tables
    ubrd_param_map{HybridKind::POINTER, "tpsyn_ubrd_map"},
    ubrd_param_map_bounds{HybridKind::POINTER, "tpsyn_ubrd_map_bnd"},
    cimp_param_map{HybridKind::POINTER, "tpsyn_cimp_map"},
    cimp_param_map_bounds{HybridKind::POINTER, "tpsyn_cimp_map_bnd"},
    cmap_param_map{HybridKind::POINTER, "tpsyn_cmap_map"},
    cmap_param_map_bounds{HybridKind::POINTER, "tpsyn_cmap_map_bnd"},
    bond_param_map{HybridKind::POINTER, "tpsyn_bond_map"},
    bond_param_map_bounds{HybridKind::POINTER, "tpsyn_bond_map_bnd"},
    angl_param_map{HybridKind::POINTER, "tpsyn_angl_map"},
    angl_param_map_bounds{HybridKind::POINTER, "tpsyn_angl_map_bnd"},
    dihe_param_map{HybridKind::POINTER, "tpsyn_dihe_map"},
    dihe_param_map_bounds{HybridKind::POINTER, "tpsyn_dihe_map_bnd"},
    attn14_param_map{HybridKind::POINTER, "tpsyn_attn_map"},
    attn14_param_map_bounds{HybridKind::POINTER, "tpsyn_attn_map_bnd"},
    vste_param_map{HybridKind::POINTER, "tpsyn_vsite_map"},
    vste_param_map_bounds{HybridKind::POINTER, "tpsyn_vste_map_bnd"},
    sett_param_map{HybridKind::POINTER, "tpsyn_settle_map"},
    sett_param_map_bounds{HybridKind::POINTER, "tpsyn_sett_map_bnd"},
    cnst_param_map{HybridKind::POINTER, "tpsyn_constraint_map"},
    cnst_param_map_bounds{HybridKind::POINTER, "tpsyn_cnst_map_bnd"},

    // Force field consensus table condensed parameters
    ubrd_stiffnesses{HybridKind::POINTER, "tpsyn_ub_stiff"},
    ubrd_equilibria{HybridKind::POINTER, "tpsyn_ub_equil"},
    cimp_stiffnesses{HybridKind::POINTER, "tpsyn_cimp_stiff"},
    cimp_phase_angles{HybridKind::POINTER, "tpsyn_cimp_equil"},
    cmap_surface_dimensions{HybridKind::POINTER, "tpsyn_cmap_dims"},
    cmap_surface_bounds{HybridKind::POINTER, "tpsyn_cmap_bnd"},
    cmap_patch_bounds{HybridKind::POINTER, "tpsyn_cmpatch_bnd"},
    cmap_surfaces{HybridKind::POINTER, "tpsyn_cmap_surf"},
    cmap_patches{HybridKind::POINTER, "tpsyn_cmap_patch"},
    sp_ubrd_stiffnesses{HybridKind::POINTER, "tpsyn_ub_stiff_sp"},
    sp_ubrd_equilibria{HybridKind::POINTER, "tpsyn_ub_equil_sp"},
    sp_cimp_stiffnesses{HybridKind::POINTER, "tpsyn_cimp_stiff_sp"},
    sp_cimp_phase_angles{HybridKind::POINTER, "tpsyn_cimp_equil_sp"},
    sp_cmap_surfaces{HybridKind::POINTER, "tpsyn_cmap_surf_sp"},
    sp_cmap_patches{HybridKind::POINTER, "tpsyn_cmap_patch_sp"},
    bond_stiffnesses{HybridKind::POINTER, "tpsyn_bondk"},
    bond_equilibria{HybridKind::POINTER, "tpsyn_bondl0"},
    angl_stiffnesses{HybridKind::POINTER, "tpsyn_anglk"},
    angl_equilibria{HybridKind::POINTER, "tpsyn_anglt0"},
    dihe_amplitudes{HybridKind::POINTER, "tpsyn_dihek"},
    dihe_periodicities{HybridKind::POINTER, "tpsyn_dihen"},
    dihe_phase_angles{HybridKind::POINTER, "tpsyn_dihepsi"},
    attn14_elec_factors{HybridKind::POINTER, "tpsyn_elec14_scale"},
    attn14_vdw_factors{HybridKind::POINTER, "tpsyn_vdw14_scale"},
    sp_bond_stiffnesses{HybridKind::POINTER, "tpsyn_bondk_sp"},
    sp_bond_equilibria{HybridKind::POINTER, "tpsyn_bondl0_sp"},
    sp_angl_stiffnesses{HybridKind::POINTER, "tpsyn_anglk_sp"},
    sp_angl_equilibria{HybridKind::POINTER, "tpsyn_anglt0_sp"},
    sp_dihe_amplitudes{HybridKind::POINTER, "tpsyn_dihek_sp"},
    sp_dihe_periodicities{HybridKind::POINTER, "tpsyn_dihen_sp"},
    sp_dihe_phase_angles{HybridKind::POINTER, "tpsyn_dihepsi_sp"},
    sp_attn14_elec_factors{HybridKind::POINTER, "tpsynf_elec14_scale"},
    sp_attn14_vdw_factors{HybridKind::POINTER, "tpsynf_vdw14_scale"},
    valparam_double_data{HybridKind::ARRAY, "tpsyn_vparm_dbl"},
    valparam_float_data{HybridKind::ARRAY, "tpsyn_vparm_flt"},
    valparam_int_data{HybridKind::ARRAY, "tpsyn_vparm_int"},
    valparam_int2_data{HybridKind::ARRAY, "tpsyn_vparm_int2"},

    // Atoms participating in each valence term, plus the parameter indices of each term.  Atom
    // indices refer to concatenated lists of all systems.  Parameter indices refer to the
    // consensus tables in this object.
    ubrd_i_atoms{HybridKind::POINTER, "tpsyn_ubrd_i"},
    ubrd_k_atoms{HybridKind::POINTER, "tpsyn_ubrd_k"},
    ubrd_param_idx{HybridKind::POINTER, "tpsyn_ubrd_idx"},
    cimp_i_atoms{HybridKind::POINTER, "tpsyn_cimp_i"},
    cimp_j_atoms{HybridKind::POINTER, "tpsyn_cimp_j"},
    cimp_k_atoms{HybridKind::POINTER, "tpsyn_cimp_k"},
    cimp_l_atoms{HybridKind::POINTER, "tpsyn_cimp_l"},
    cimp_param_idx{HybridKind::POINTER, "tpsyn_cimp_idx"},
    cmap_i_atoms{HybridKind::POINTER, "tpsyn_cmap_i"},
    cmap_j_atoms{HybridKind::POINTER, "tpsyn_cmap_j"},
    cmap_k_atoms{HybridKind::POINTER, "tpsyn_cmap_k"},
    cmap_l_atoms{HybridKind::POINTER, "tpsyn_cmap_l"},
    cmap_m_atoms{HybridKind::POINTER, "tpsyn_cmap_m"},
    cmap_param_idx{HybridKind::POINTER, "tpsyn_cmap_idx"},
    bond_i_atoms{HybridKind::POINTER, "tpsyn_bond_i"},
    bond_j_atoms{HybridKind::POINTER, "tpsyn_bond_j"},
    bond_param_idx{HybridKind::POINTER, "tpsyn_bond_idx"},
    angl_i_atoms{HybridKind::POINTER, "tpsyn_angl_i"},
    angl_j_atoms{HybridKind::POINTER, "tpsyn_angl_j"},
    angl_k_atoms{HybridKind::POINTER, "tpsyn_angl_k"},
    angl_param_idx{HybridKind::POINTER, "tpsyn_angl_idx"},
    dihe_i_atoms{HybridKind::POINTER, "tpsyn_dihe_i"},
    dihe_j_atoms{HybridKind::POINTER, "tpsyn_dihe_j"},
    dihe_k_atoms{HybridKind::POINTER, "tpsyn_dihe_k"},
    dihe_l_atoms{HybridKind::POINTER, "tpsyn_dihe_l"},
    dihe_param_idx{HybridKind::POINTER, "tpsyn_dihe_idx"},
    valence_int_data{HybridKind::ARRAY, "tpsyn_val_ints"},

    // Non-bonded parameters for all systems
    charge_indices{HybridKind::POINTER, "tpsyn_q_idx"},
    lennard_jones_indices{HybridKind::POINTER, "tpsyn_lj_idx"},
    charge_parameters{HybridKind::ARRAY, "tpsyn_q_parm"},
    lennard_jones_ab_coeff{HybridKind::ARRAY, "tpsyn_lj_ab"},
    lennard_jones_c_coeff{HybridKind::ARRAY, "tpsyn_lj_c"},
    lennard_jones_14_a_coeff{HybridKind::ARRAY, "tpsyn_lj_14_a"},
    lennard_jones_14_b_coeff{HybridKind::ARRAY, "tpsyn_lj_14_b"},
    lennard_jones_14_c_coeff{HybridKind::ARRAY, "tpsyn_lj_14_c"},
    lennard_jones_sigma{HybridKind::ARRAY, "tpsyn_lj_sigma"},
    lennard_jones_14_sigma{HybridKind::ARRAY, "tpsyn_lj_14_sigma"},
    sp_charge_parameters{HybridKind::ARRAY, "tpsyn_q_parm_sp"},
    sp_lennard_jones_ab_coeff{HybridKind::ARRAY, "tpsyn_lj_ab_sp"},
    sp_lennard_jones_c_coeff{HybridKind::ARRAY, "tpsyn_lj_c_sp"},
    sp_lennard_jones_14_a_coeff{HybridKind::ARRAY, "tpsyn_lj_14_a_sp"},
    sp_lennard_jones_14_b_coeff{HybridKind::ARRAY, "tpsyn_lj_14_b_sp"},
    sp_lennard_jones_14_c_coeff{HybridKind::ARRAY, "tpsyn_lj_14_c_sp"},
    sp_lennard_jones_sigma{HybridKind::ARRAY, "tpsyn_lj_sigma_sp"},
    sp_lennard_jones_14_sigma{HybridKind::ARRAY, "tpsyn_lj_14_sigma_sp"},

    // Implicit solvent model parameters
    neck_table_size{0},
    neck_gb_indices{HybridKind::POINTER, "tpsyn_gb_idx"},
    atomic_pb_radii{HybridKind::POINTER, "tpsyn_atom_pb_rad"},
    gb_screening_factors{HybridKind::POINTER, "tpsyn_gbscreen"},
    gb_alpha_parameters{HybridKind::POINTER, "tpsyn_gbalpha"},
    gb_beta_parameters{HybridKind::POINTER, "tpsyn_gbbeta"},
    gb_gamma_parameters{HybridKind::POINTER, "tpsyn_gbgamma"},
    neck_limit_tables{HybridKind::ARRAY, "tpsyn_neck_limits"},
    sp_atomic_pb_radii{HybridKind::POINTER, "tpsyn_atom_pb_rad_sp"},
    sp_gb_screening_factors{HybridKind::POINTER, "tpsyn_gbscreen_sp"},
    sp_gb_alpha_parameters{HybridKind::POINTER, "tpsyn_gbalpha_sp"},
    sp_gb_beta_parameters{HybridKind::POINTER, "tpsyn_gbbeta_sp"},
    sp_gb_gamma_parameters{HybridKind::POINTER, "tpsyn_gbgamma_sp"},    
    sp_neck_limit_tables{HybridKind::ARRAY, "tpsyn_neck_limits_sp"},
    
    // Restraint parameters for all systems
    rposn_step_bounds{HybridKind::POINTER, "tpsyn_rposn_steps"},
    rbond_step_bounds{HybridKind::POINTER, "tpsyn_rbond_steps"},
    rangl_step_bounds{HybridKind::POINTER, "tpsyn_rangl_steps"},
    rdihe_step_bounds{HybridKind::POINTER, "tpsyn_rdihe_steps"},
    rposn_init_k{HybridKind::POINTER, "tpsyn_rposn_init_k"},
    rposn_final_k{HybridKind::POINTER, "tpsyn_rposn_finl_k"},
    rposn_init_r{HybridKind::POINTER, "tpsyn_rposn_init_k"},
    rposn_final_r{HybridKind::POINTER, "tpsyn_rposn_finl_k"},
    rposn_init_xy{HybridKind::POINTER, "tpsyn_rposn_init_xy"},
    rposn_init_z{HybridKind::POINTER, "tpsyn_rposn_init_z"},
    rposn_final_xy{HybridKind::POINTER, "tpsyn_rposn_final_xy"},
    rposn_final_z{HybridKind::POINTER, "tpsyn_rposn_final_z"},
    rbond_init_k{HybridKind::POINTER, "tpsyn_rbond_init_k"},
    rbond_final_k{HybridKind::POINTER, "tpsyn_rbond_finl_k"},
    rbond_init_r{HybridKind::POINTER, "tpsyn_rbond_init_k"},
    rbond_final_r{HybridKind::POINTER, "tpsyn_rbond_finl_k"},
    rangl_init_k{HybridKind::POINTER, "tpsyn_rangl_init_k"},
    rangl_final_k{HybridKind::POINTER, "tpsyn_rangl_finl_k"},
    rangl_init_r{HybridKind::POINTER, "tpsyn_rangl_init_k"},
    rangl_final_r{HybridKind::POINTER, "tpsyn_rangl_finl_k"},
    rdihe_init_k{HybridKind::POINTER, "tpsyn_rdihe_init_k"},
    rdihe_final_k{HybridKind::POINTER, "tpsyn_rdihe_finl_k"},
    rdihe_init_r{HybridKind::POINTER, "tpsyn_rdihe_init_k"},
    rdihe_final_r{HybridKind::POINTER, "tpsyn_rdihe_finl_k"},
    sp_rposn_init_k{HybridKind::POINTER, "tpsynf_rposn_init_k"},
    sp_rposn_final_k{HybridKind::POINTER, "tpsynf_rposn_finl_k"},
    sp_rposn_init_r{HybridKind::POINTER, "tpsynf_rposn_init_k"},
    sp_rposn_final_r{HybridKind::POINTER, "tpsynf_rposn_finl_k"},
    sp_rposn_init_xy{HybridKind::POINTER, "tpsynf_rposn_init_xy"},
    sp_rposn_init_z{HybridKind::POINTER, "tpsynf_rposn_init_z"},
    sp_rposn_final_xy{HybridKind::POINTER, "tpsynf_rposn_final_xy"},
    sp_rposn_final_z{HybridKind::POINTER, "tpsynf_rposn_final_z"},
    sp_rbond_init_k{HybridKind::POINTER, "tpsynf_rbond_init_k"},
    sp_rbond_final_k{HybridKind::POINTER, "tpsynf_rbond_finl_k"},
    sp_rbond_init_r{HybridKind::POINTER, "tpsynf_rbond_init_k"},
    sp_rbond_final_r{HybridKind::POINTER, "tpsynf_rbond_finl_k"},
    sp_rangl_init_k{HybridKind::POINTER, "tpsynf_rangl_init_k"},
    sp_rangl_final_k{HybridKind::POINTER, "tpsynf_rangl_finl_k"},
    sp_rangl_init_r{HybridKind::POINTER, "tpsynf_rangl_init_k"},
    sp_rangl_final_r{HybridKind::POINTER, "tpsynf_rangl_finl_k"},
    sp_rdihe_init_k{HybridKind::POINTER, "tpsynf_rdihe_init_k"},
    sp_rdihe_final_k{HybridKind::POINTER, "tpsynf_rdihe_finl_k"},
    sp_rdihe_init_r{HybridKind::POINTER, "tpsynf_rdihe_init_k"},
    sp_rdihe_final_r{HybridKind::POINTER, "tpsynf_rdihe_finl_k"},
    nmr_double_data{HybridKind::ARRAY, "tpsyn_nmr_dbl_data"},
    nmr_double2_data{HybridKind::ARRAY, "tpsyn_nmr_dbl2_data"},
    nmr_double4_data{HybridKind::ARRAY, "tpsyn_nmr_dbl4_data"},
    nmr_float_data{HybridKind::ARRAY, "tpsyn_nmr_flt_data"},
    nmr_float2_data{HybridKind::ARRAY, "tpsyn_nmr_flt2_data"},
    nmr_float4_data{HybridKind::ARRAY, "tpsyn_nmr_flt4_data"},

    // Restraint term atom indices and parameter indices.  Like standard force field terms, the
    // atom indices refer to the concatenated list of all systems and the parameter indices refer
    // to the consensus tables.
    rposn_atoms{HybridKind::POINTER, "tpsyn_rposn_at"},
    rposn_kr_param_idx{HybridKind::POINTER, "tpsyn_rposn_kr_idx"},
    rposn_xyz_param_idx{HybridKind::POINTER, "tpsyn_rposn_xyz_idx"},
    rbond_i_atoms{HybridKind::POINTER, "tpsyn_rbond_iat"},
    rbond_j_atoms{HybridKind::POINTER, "tpsyn_rbond_jat"},
    rbond_param_idx{HybridKind::POINTER, "tpsyn_rbond_param"},
    rangl_i_atoms{HybridKind::POINTER, "tpsyn_rangl_iat"},
    rangl_j_atoms{HybridKind::POINTER, "tpsyn_rangl_jat"},
    rangl_k_atoms{HybridKind::POINTER, "tpsyn_rangl_kat"},
    rangl_param_idx{HybridKind::POINTER, "tpsyn_rangl_param"},
    rdihe_i_atoms{HybridKind::POINTER, "tpsyn_rdihe_iat"},
    rdihe_j_atoms{HybridKind::POINTER, "tpsyn_rdihe_jat"},
    rdihe_k_atoms{HybridKind::POINTER, "tpsyn_rdihe_kat"},
    rdihe_l_atoms{HybridKind::POINTER, "tpsyn_rdihe_lat"},
    rdihe_param_idx{HybridKind::POINTER, "tpsyn_rdihe_param"},
    rposn_kr_param_map{HybridKind::POINTER, "tpsyn_rposn_kr_map"},
    rposn_xyz_param_map{HybridKind::POINTER, "tpsyn_rposn_xyz_map"},
    rposn_param_map_bounds{HybridKind::POINTER, "tpsyn_rposn_map_bnd"},
    rbond_param_map{HybridKind::POINTER, "tpsyn_rbond_map"},
    rbond_param_map_bounds{HybridKind::POINTER, "tpsyn_rbond_map_bnd"},
    rangl_param_map{HybridKind::POINTER, "tpsyn_rangl_map"},
    rangl_param_map_bounds{HybridKind::POINTER, "tpsyn_rangl_map_bnd"},
    rdihe_param_map{HybridKind::POINTER, "tpsyn_rdihe_map"},
    rdihe_param_map_bounds{HybridKind::POINTER, "tpsyn_rdihe_map_bnd"},
    nmr_int_data{HybridKind::ARRAY, "tpsyn_nmr_ints"},
    nmr_int2_data{HybridKind::ARRAY, "tpsyn_nmr_int2_data"},

    // Virtual site parameter arrays, term atom indices, and parameter indices
    virtual_site_parameters{HybridKind::ARRAY, "tpsyn_vs_params"},
    sp_virtual_site_parameters{HybridKind::ARRAY, "tpsynf_vs_params"},
    virtual_site_atoms{HybridKind::POINTER, "tpsyn_vs_atoms"},
    virtual_site_frame1_atoms{HybridKind::POINTER, "tpsyn_vs_frame1"},
    virtual_site_frame2_atoms{HybridKind::POINTER, "tpsyn_vs_frame2"},
    virtual_site_frame3_atoms{HybridKind::POINTER, "tpsyn_vs_frame3"},
    virtual_site_frame4_atoms{HybridKind::POINTER, "tpsyn_vs_frame4"},
    virtual_site_parameter_indices{HybridKind::POINTER, "tpsyn_vs_param_idx"},
    vsite_int_data{HybridKind::ARRAY, "tpsyn_vsite_ints"},

    // SETTLE (analytic rigid water constraints) parameter arrays, plus fused atom and parameter
    // indices.
    settle_group_geometry{HybridKind::ARRAY, "tpsyn_settle_geom"},
    settle_group_masses{HybridKind::ARRAY, "tpsyn_settle_mass"},
    sp_settle_group_geometry{HybridKind::ARRAY, "tpsynf_settle_geom"},
    sp_settle_group_masses{HybridKind::ARRAY, "tpsynf_settle_mass"},
    settle_group_indexing{HybridKind::ARRAY, "tpsyn_settle_idx"},

    // Constraint group parameter indices and bounds (referencing consensus tables from the
    // synthesis) plus length / inverse mass parameters (also fused)
    constraint_group_indices{HybridKind::ARRAY, "tpsyn_cnst_atom_idx"},
    constraint_group_bounds{HybridKind::ARRAY, "tpsyn_cnst_atom_bnd"},
    constraint_group_param_idx{HybridKind::ARRAY, "tpsyn_cnst_parm_idx"},
    constraint_param_bounds{HybridKind::ARRAY, "tpsyn_cnst_parm_bnd"},
    constraint_group_params{HybridKind::ARRAY, "tpsyn_cnst_lm"},
    sp_constraint_group_params{HybridKind::ARRAY, "tpsynf_cnst_lm"},

    // Valence work unit instruction sets and energy accumulation masks
    total_valence_work_units{0}, valence_work_unit_size{maximum_valence_work_unit_atoms},
    valence_thread_block_size{ValenceKernelSize::XL},
    vwu_instruction_sets{HybridKind::ARRAY, "tpsyn_vwu_insr_sets"},
    vwu_import_lists{HybridKind::ARRAY, "tpsyn_vwu_imports"},
    vwu_manipulation_masks{HybridKind::POINTER, "tpsyn_vwu_manip"},
    cbnd_instructions{HybridKind::POINTER, "tpsyn_cbnd_insr"},
    angl_instructions{HybridKind::POINTER, "tpsyn_angl_insr"},
    cdhe_instructions{HybridKind::POINTER, "tpsyn_cdhe_insr"},
    cdhe_overtones{HybridKind::POINTER, "tpsyn_ovrt_insr"},
    cmap_instructions{HybridKind::POINTER, "tpsyn_cmap_insr"},
    infr14_instructions{HybridKind::POINTER, "tpsyn_infr14_insr"},
    rposn_instructions{HybridKind::POINTER, "tpsyn_nmr1_insr"},
    rbond_instructions{HybridKind::POINTER, "tpsyn_nmr2_insr"},
    rangl_instructions{HybridKind::POINTER, "tpsyn_nmr3_insr"},
    rdihe_instructions{HybridKind::POINTER, "tpsyn_nmr4_insr"},
    vste_instructions{HybridKind::POINTER, "tpsyn_vste_insr"},
    sett_instructions{HybridKind::POINTER, "tpsyn_sett_insr"},
    cnst_instructions{HybridKind::POINTER, "tpsyn_cnst_insr"},
    accumulate_cbnd_energy{HybridKind::POINTER, "tpsyn_cbnd_edir"},
    accumulate_angl_energy{HybridKind::POINTER, "tpsyn_angl_edir"},
    accumulate_cdhe_energy{HybridKind::POINTER, "tpsyn_cdhe_edir"},
    accumulate_cmap_energy{HybridKind::POINTER, "tpsyn_cmap_edir"},
    accumulate_infr14_energy{HybridKind::POINTER, "tpsyn_infr14_edir"},
    accumulate_rposn_energy{HybridKind::POINTER, "tpsyn_rposn_edir"},
    accumulate_rbond_energy{HybridKind::POINTER, "tpsyn_rbond_edir"},
    accumulate_rangl_energy{HybridKind::POINTER, "tpsyn_rangl_edir"},
    accumulate_rdihe_energy{HybridKind::POINTER, "tpsyn_rdihe_edir"},
    insr_uint_data{HybridKind::ARRAY, "tpsyn_insr_data1"},
    insr_uint2_data{HybridKind::ARRAY, "tpsyn_insr_data2"},

    // Non-bonded work units and their instructions
    total_nonbonded_work_units{0},
    nonbonded_work_type{NbwuKind::UNKNOWN},
    nonbonded_abstracts{HybridKind::ARRAY, "tpsyn_nbwu_abstract"},
    nbwu_instructions{HybridKind::ARRAY, "tpsyn_nbwu_insr"},

    // Reduction work units and their abstracts
    total_reduction_work_units{0},
    rdwu_per_system{RdwuPerSystem::ONE},
    reduction_abstracts{HybridKind::ARRAY, "tpsyn_rdwu_abstract"},
    
    // Pointer to the timer used to track time needed for assembling this object.  The timer can
    // also be used to track access and other usage of the object.
    timer{timer_in}
{
  // Set up dummy restaint apparatuses for systems that do not have a restraint apparatus.  These
  // restraint apparatuses will exist for the duration of the object.
  const std::vector<int> new_restraint_indices = createDummyRestraints(restraint_indices_in,
                                                                       topology_indices_in);

  // Setup and memory layout
  const std::vector<int> topology_index_rebase = checkTopologyList(topology_indices_in);
  const std::vector<int> restraint_index_rebase = checkRestraintList(new_restraint_indices,
                                                                     topology_indices_in,
                                                                     topology_index_rebase);
  checkCommonSettings();
  int term_array_timings, param_cond_timings, vwu_creation_timings;
  if (timer != nullptr) {
    term_array_timings   = timer->addCategory("[AtomGraphSynthesis] Build term arrays");
    param_cond_timings   = timer->addCategory("[AtomGraphSynthesis] Condense parameters");
    vwu_creation_timings = timer->addCategory("[AtomGraphSynthesis] VWU creation");
    timer->assignTime(0);
  }
  buildAtomAndTermArrays(topology_indices_in, topology_index_rebase, new_restraint_indices,
                         restraint_index_rebase);
  if (timer != nullptr) timer->assignTime(term_array_timings);

  // Condense valence and charge parameters into compact tables of unique values
  condenseParameterTables();
  
  // The unique Lennard-Jones parameters are hard to map out.  To the degree that there are
  // unique parameters in each system, the Lennard-Jones tables might as well be block matrices.
  // Keep a list of the A, B, and possibly C coefficients of each Lennard-Jones type interacting
  // with all other types.
  extendLJMatrices();

  // Restraint data can now be incorporated.
  condenseRestraintNetworks();
  if (timer != nullptr) timer->assignTime(param_cond_timings);

  // Create valence work units for all topologies, then load them into the synthesis
  valence_work_unit_size = calculateValenceWorkUnitSize(atom_counts, gpu.getSMPCount(),
                                                        &valence_thread_block_size);
  loadValenceWorkUnits(valence_work_unit_size);
  if (timer != nullptr) timer->assignTime(vwu_creation_timings);

  // Create reduction work units for all topologies, then load them into the synthesis.  Like the
  // process of creating the valence work units, the original array of C++ objects is destroyed
  // upon completion of this constructor.  The critical information is retained in an array of
  // abstracts.
  loadReductionWorkUnits();

  // Apply an implicit solvent model if one is already present in the underlying topologies.
  setImplicitSolventModel();
}

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis::AtomGraphSynthesis(const std::vector<AtomGraph*> &topologies_in,
                                       const std::vector<RestraintApparatus*> &restraints_in,
                                       const ExceptionResponse policy_in, const GpuDetails &gpu,
                                       StopWatch *timer_in) :
  AtomGraphSynthesis(topologies_in, restraints_in,
                     incrementingSeries<int>(0, topologies_in.size()),
                     incrementingSeries<int>(0, topologies_in.size()), policy_in, gpu, timer_in)
{}

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis::AtomGraphSynthesis(const std::vector<AtomGraph*> &topologies_in,
                                       const std::vector<int> &topology_indices_in,
                                       const ExceptionResponse policy_in, const GpuDetails &gpu,
                                       StopWatch *timer_in) :
    AtomGraphSynthesis(topologies_in, std::vector<RestraintApparatus*>(1, nullptr),
                       topology_indices_in, std::vector<int>(topology_indices_in.size(), 0),
                       policy_in, gpu, timer_in)
{}

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis::AtomGraphSynthesis(const std::vector<AtomGraph*> &topologies_in,
                                       const ExceptionResponse policy_in, const GpuDetails &gpu,
                                       StopWatch *timer_in) :
    AtomGraphSynthesis(topologies_in, std::vector<RestraintApparatus*>(1, nullptr),
                       incrementingSeries<int>(0, topologies_in.size()),
                       std::vector<int>(topologies_in.size(), 0), policy_in, gpu, timer_in)
{}

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis::AtomGraphSynthesis(const AtomGraphSynthesis &original) :

    // Counts spanning all topologies
    policy{original.policy},
    topology_count{original.topology_count},
    restraint_network_count{original.restraint_network_count},
    system_count{original.system_count},
    total_atoms{original.total_atoms},
    total_virtual_sites{original.total_virtual_sites},
    total_bond_terms{original.total_bond_terms},
    total_angl_terms{original.total_angl_terms},
    total_dihe_terms{original.total_dihe_terms},
    total_ubrd_terms{original.total_ubrd_terms},
    total_cimp_terms{original.total_cimp_terms},
    total_cmap_terms{original.total_cmap_terms},
    total_lj_types{original.total_lj_types},
    total_charge_types{original.total_charge_types},
    total_bond_params{original.total_bond_params},
    total_angl_params{original.total_angl_params},
    total_dihe_params{original.total_dihe_params},
    total_ubrd_params{original.total_ubrd_params},
    total_cimp_params{original.total_cimp_params},
    total_cmap_surfaces{original.total_cmap_surfaces},
    total_attn14_params{original.total_attn14_params},
    total_vste_params{original.total_vste_params},
    total_position_restraints{original.total_position_restraints},
    total_distance_restraints{original.total_distance_restraints},
    total_angle_restraints{original.total_angle_restraints},
    total_dihedral_restraints{original.total_dihedral_restraints},

    // Descriptors spanning all systems
    periodic_box_class{original.periodic_box_class},
    gb_style{original.gb_style},
    dielectric_constant{original.dielectric_constant},
    is_kappa{original.is_kappa},
    salt_concentration{original.salt_concentration},
    gb_offset{original.gb_offset},
    gb_neckscale{original.gb_neckscale},
    gb_neckcut{original.gb_neckcut},
    coulomb_constant{original.coulomb_constant},
    use_bond_constraints{original.use_bond_constraints},
    use_settle{original.use_settle},
    largest_constraint_group{original.largest_constraint_group},
    water_residue_name{original.water_residue_name},
    pb_radii_sets{original.pb_radii_sets},

    // Individual topologies and restraint networks (RestraintApparatus objects), each of them
    // describing one or more of the individual systems within the synthesis
    topologies{original.topologies},
    restraint_networks{original.restraint_networks},
    restraint_dummies{original.restraint_dummies},
    topology_indices{original.topology_indices},
    restraint_indices{original.restraint_indices},
    atom_counts{original.atom_counts},
    residue_counts{original.residue_counts},
    molecule_counts{original.molecule_counts},
    largest_residue_sizes{original.largest_residue_sizes},
    last_solute_residues{original.last_solute_residues},
    last_solute_atoms{original.last_solute_atoms},
    first_solvent_molecules{original.first_solvent_molecules},

    // Counts of energy terms and particle-related quantities
    ubrd_term_counts{original.ubrd_term_counts},
    cimp_term_counts{original.cimp_term_counts},
    cmap_term_counts{original.cmap_term_counts},
    bond_term_counts{original.bond_term_counts},
    angl_term_counts{original.angl_term_counts},
    dihe_term_counts{original.dihe_term_counts},
    virtual_site_counts{original.virtual_site_counts},
    posn_restraint_counts{original.posn_restraint_counts},
    bond_restraint_counts{original.bond_restraint_counts},
    angl_restraint_counts{original.angl_restraint_counts},
    dihe_restraint_counts{original.dihe_restraint_counts},
    lj_type_counts{original.lj_type_counts},
    total_exclusion_counts{original.total_exclusion_counts},
    rigid_water_counts{original.rigid_water_counts},
    bond_constraint_counts{original.bond_constraint_counts},
    degrees_of_freedom{original.degrees_of_freedom},
    cnst_degrees_of_freedom{original.cnst_degrees_of_freedom},
    nonrigid_particle_counts{original.nonrigid_particle_counts},

    // Offsets for each system's term lists
    atom_offsets{original.atom_offsets},
    atom_bit_offsets{original.atom_bit_offsets},
    residue_offsets{original.residue_offsets},
    molecule_offsets{original.molecule_offsets},
    ubrd_term_offsets{original.ubrd_term_offsets},
    cimp_term_offsets{original.cimp_term_offsets},
    cmap_term_offsets{original.cmap_term_offsets},
    bond_term_offsets{original.bond_term_offsets},
    angl_term_offsets{original.angl_term_offsets},
    dihe_term_offsets{original.dihe_term_offsets},
    virtual_site_offsets{original.virtual_site_offsets},
    posn_restraint_offsets{original.posn_restraint_offsets},
    bond_restraint_offsets{original.bond_restraint_offsets},
    angl_restraint_offsets{original.angl_restraint_offsets},
    dihe_restraint_offsets{original.dihe_restraint_offsets},
    sett_group_offsets{original.sett_group_offsets},
    cnst_group_offsets{original.cnst_group_offsets},
    nb_exclusion_offsets{original.nb_exclusion_offsets},
    lennard_jones_abc_offsets{original.lennard_jones_abc_offsets},
    int_system_data{original.int_system_data},

    // Atom and residue details
    residue_limits{original.residue_limits},
    atom_struc_numbers{original.atom_struc_numbers},
    residue_numbers{original.residue_numbers},
    molecule_limits{original.molecule_limits},
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
    chem_int_data{original.chem_int_data},
    chem_int2_data{original.chem_int2_data},
    chem_double_data{original.chem_double_data},
    chem_float_data{original.chem_float_data},
    chem_char4_data{original.chem_char4_data},

    // Force field parameter maps from individual systems to the consensus tables
    ubrd_param_map{original.ubrd_param_map},
    ubrd_param_map_bounds{original.ubrd_param_map_bounds},
    cimp_param_map{original.cimp_param_map},
    cimp_param_map_bounds{original.cimp_param_map_bounds},
    cmap_param_map{original.cmap_param_map},
    cmap_param_map_bounds{original.cmap_param_map_bounds},
    bond_param_map{original.bond_param_map},
    bond_param_map_bounds{original.bond_param_map_bounds},
    angl_param_map{original.angl_param_map},
    angl_param_map_bounds{original.angl_param_map_bounds},
    dihe_param_map{original.dihe_param_map},
    dihe_param_map_bounds{original.dihe_param_map_bounds},
    attn14_param_map{original.attn14_param_map},
    attn14_param_map_bounds{original.attn14_param_map_bounds},
    vste_param_map{original.vste_param_map},
    vste_param_map_bounds{original.vste_param_map_bounds},
    sett_param_map{original.sett_param_map},
    sett_param_map_bounds{original.sett_param_map_bounds},
    cnst_param_map{original.cnst_param_map},
    cnst_param_map_bounds{original.cnst_param_map_bounds},

    // Force field consensus table condensed parameters
    ubrd_stiffnesses{original.ubrd_stiffnesses},
    ubrd_equilibria{original.ubrd_equilibria},
    cimp_stiffnesses{original.cimp_stiffnesses},
    cimp_phase_angles{original.cimp_phase_angles},
    cmap_surface_dimensions{original.cmap_surface_dimensions},
    cmap_surface_bounds{original.cmap_surface_bounds},
    cmap_patch_bounds{original.cmap_patch_bounds},
    cmap_surfaces{original.cmap_surfaces},
    cmap_patches{original.cmap_patches},
    sp_ubrd_stiffnesses{original.sp_ubrd_stiffnesses},
    sp_ubrd_equilibria{original.sp_ubrd_equilibria},
    sp_cimp_stiffnesses{original.sp_cimp_stiffnesses},
    sp_cimp_phase_angles{original.sp_cimp_phase_angles},
    sp_cmap_surfaces{original.sp_cmap_surfaces},
    sp_cmap_patches{original.sp_cmap_patches},
    bond_stiffnesses{original.bond_stiffnesses},
    bond_equilibria{original.bond_equilibria},
    angl_stiffnesses{original.angl_stiffnesses},
    angl_equilibria{original.angl_equilibria},
    dihe_amplitudes{original.dihe_amplitudes},
    dihe_periodicities{original.dihe_periodicities},
    dihe_phase_angles{original.dihe_phase_angles},
    attn14_elec_factors{original.attn14_elec_factors},
    attn14_vdw_factors{original.attn14_vdw_factors},
    sp_bond_stiffnesses{original.sp_bond_stiffnesses},
    sp_bond_equilibria{original.sp_bond_equilibria},
    sp_angl_stiffnesses{original.sp_angl_stiffnesses},
    sp_angl_equilibria{original.sp_angl_equilibria},
    sp_dihe_amplitudes{original.sp_dihe_amplitudes},
    sp_dihe_periodicities{original.sp_dihe_periodicities},
    sp_dihe_phase_angles{original.sp_dihe_phase_angles},
    sp_attn14_elec_factors{original.sp_attn14_elec_factors},
    sp_attn14_vdw_factors{original.sp_attn14_vdw_factors},
    valparam_double_data{original.valparam_double_data},
    valparam_float_data{original.valparam_float_data},
    valparam_int_data{original.valparam_int_data},
    valparam_int2_data{original.valparam_int2_data},

    // Atoms participating in each valence term, plus the parameter indices of each term
    ubrd_i_atoms{original.ubrd_i_atoms},
    ubrd_k_atoms{original.ubrd_k_atoms},
    ubrd_param_idx{original.ubrd_param_idx},
    cimp_i_atoms{original.cimp_i_atoms},
    cimp_j_atoms{original.cimp_j_atoms},
    cimp_k_atoms{original.cimp_k_atoms},
    cimp_l_atoms{original.cimp_l_atoms},
    cimp_param_idx{original.cimp_param_idx},
    cmap_i_atoms{original.cmap_i_atoms},
    cmap_j_atoms{original.cmap_j_atoms},
    cmap_k_atoms{original.cmap_k_atoms},
    cmap_l_atoms{original.cmap_l_atoms},
    cmap_m_atoms{original.cmap_m_atoms},
    cmap_param_idx{original.cmap_param_idx},
    bond_i_atoms{original.bond_i_atoms},
    bond_j_atoms{original.bond_j_atoms},
    bond_param_idx{original.bond_param_idx},
    angl_i_atoms{original.angl_i_atoms},
    angl_j_atoms{original.angl_j_atoms},
    angl_k_atoms{original.angl_k_atoms},
    angl_param_idx{original.angl_param_idx},
    dihe_i_atoms{original.dihe_i_atoms},
    dihe_j_atoms{original.dihe_j_atoms},
    dihe_k_atoms{original.dihe_k_atoms},
    dihe_l_atoms{original.dihe_l_atoms},
    dihe_param_idx{original.dihe_param_idx},
    valence_int_data{original.valence_int_data},

    // Non-bonded parameters for all systems
    charge_indices{original.charge_indices},
    lennard_jones_indices{original.lennard_jones_indices},
    charge_parameters{original.charge_parameters},
    lennard_jones_ab_coeff{original.lennard_jones_ab_coeff},
    lennard_jones_c_coeff{original.lennard_jones_c_coeff},
    lennard_jones_14_a_coeff{original.lennard_jones_14_a_coeff},
    lennard_jones_14_b_coeff{original.lennard_jones_14_b_coeff},
    lennard_jones_14_c_coeff{original.lennard_jones_14_c_coeff},
    lennard_jones_sigma{original.lennard_jones_sigma},
    lennard_jones_14_sigma{original.lennard_jones_14_sigma},
    sp_charge_parameters{original.sp_charge_parameters},
    sp_lennard_jones_ab_coeff{original.sp_lennard_jones_ab_coeff},
    sp_lennard_jones_c_coeff{original.sp_lennard_jones_c_coeff},
    sp_lennard_jones_14_a_coeff{original.sp_lennard_jones_14_a_coeff},
    sp_lennard_jones_14_b_coeff{original.sp_lennard_jones_14_b_coeff},
    sp_lennard_jones_14_c_coeff{original.sp_lennard_jones_14_c_coeff},
    sp_lennard_jones_sigma{original.sp_lennard_jones_sigma},
    sp_lennard_jones_14_sigma{original.sp_lennard_jones_14_sigma},

    // Implicit solvent model parameters
    neck_table_size{original.neck_table_size},
    neck_gb_indices{original.neck_gb_indices},
    atomic_pb_radii{original.atomic_pb_radii},
    gb_screening_factors{original.gb_screening_factors},
    gb_alpha_parameters{original.gb_alpha_parameters},
    gb_beta_parameters{original.gb_beta_parameters},
    gb_gamma_parameters{original.gb_gamma_parameters},
    neck_limit_tables{original.neck_limit_tables},
    sp_atomic_pb_radii{original.sp_atomic_pb_radii},
    sp_gb_screening_factors{original.sp_gb_screening_factors},
    sp_gb_alpha_parameters{original.sp_gb_alpha_parameters},
    sp_gb_beta_parameters{original.sp_gb_beta_parameters},
    sp_gb_gamma_parameters{original.sp_gb_gamma_parameters},
    sp_neck_limit_tables{original.sp_neck_limit_tables},
    
    // Restraint parameters for all systems
    rposn_step_bounds{original.rposn_step_bounds},
    rbond_step_bounds{original.rbond_step_bounds},
    rangl_step_bounds{original.rangl_step_bounds},
    rdihe_step_bounds{original.rdihe_step_bounds},
    rposn_init_k{original.rposn_init_k},
    rposn_final_k{original.rposn_final_k},
    rposn_init_r{original.rposn_init_r},
    rposn_final_r{original.rposn_final_r},
    rposn_init_xy{original.rposn_init_xy},
    rposn_init_z{original.rposn_init_z},
    rposn_final_xy{original.rposn_final_xy},
    rposn_final_z{original.rposn_final_z},
    rbond_init_k{original.rbond_init_k},
    rbond_final_k{original.rbond_final_k},
    rbond_init_r{original.rbond_init_r},
    rbond_final_r{original.rbond_final_r},
    rangl_init_k{original.rangl_init_k},
    rangl_final_k{original.rangl_final_k},
    rangl_init_r{original.rangl_init_r},
    rangl_final_r{original.rangl_final_r},
    rdihe_init_k{original.rdihe_init_k},
    rdihe_final_k{original.rdihe_final_k},
    rdihe_init_r{original.rdihe_init_r},
    rdihe_final_r{original.rdihe_final_r},
    sp_rposn_init_k{original.sp_rposn_init_k},
    sp_rposn_final_k{original.sp_rposn_final_k},
    sp_rposn_init_r{original.sp_rposn_init_r},
    sp_rposn_final_r{original.sp_rposn_final_r},
    sp_rposn_init_xy{original.sp_rposn_init_xy},
    sp_rposn_init_z{original.sp_rposn_init_z},
    sp_rposn_final_xy{original.sp_rposn_final_xy},
    sp_rposn_final_z{original.sp_rposn_final_z},
    sp_rbond_init_k{original.sp_rbond_init_k},
    sp_rbond_final_k{original.sp_rbond_final_k},
    sp_rbond_init_r{original.sp_rbond_init_r},
    sp_rbond_final_r{original.sp_rbond_final_r},
    sp_rangl_init_k{original.sp_rangl_init_k},
    sp_rangl_final_k{original.sp_rangl_final_k},
    sp_rangl_init_r{original.sp_rangl_init_r},
    sp_rangl_final_r{original.sp_rangl_final_r},
    sp_rdihe_init_k{original.sp_rdihe_init_k},
    sp_rdihe_final_k{original.sp_rdihe_final_k},
    sp_rdihe_init_r{original.sp_rdihe_init_r},
    sp_rdihe_final_r{original.sp_rdihe_final_r},
    nmr_double_data{original.nmr_double_data},
    nmr_double2_data{original.nmr_double2_data},
    nmr_double4_data{original.nmr_double4_data},
    nmr_float_data{original.nmr_float_data},
    nmr_float2_data{original.nmr_float2_data},
    nmr_float4_data{original.nmr_float4_data},

    // Restraint term atom indices and parameter indices
    rposn_atoms{original.rposn_atoms},
    rposn_kr_param_idx{original.rposn_kr_param_idx},
    rposn_xyz_param_idx{original.rposn_xyz_param_idx},
    rbond_i_atoms{original.rbond_i_atoms},
    rbond_j_atoms{original.rbond_j_atoms},
    rbond_param_idx{original.rbond_param_idx},
    rangl_i_atoms{original.rangl_i_atoms},
    rangl_j_atoms{original.rangl_j_atoms},
    rangl_k_atoms{original.rangl_k_atoms},
    rangl_param_idx{original.rangl_param_idx},
    rdihe_i_atoms{original.rdihe_i_atoms},
    rdihe_j_atoms{original.rdihe_j_atoms},
    rdihe_k_atoms{original.rdihe_k_atoms},
    rdihe_l_atoms{original.rdihe_l_atoms},
    rdihe_param_idx{original.rdihe_param_idx},
    rposn_kr_param_map{original.rposn_kr_param_map},
    rposn_xyz_param_map{original.rposn_xyz_param_map},
    rposn_param_map_bounds{original.rposn_param_map_bounds},
    rbond_param_map{original.rbond_param_map},
    rbond_param_map_bounds{original.rbond_param_map_bounds},
    rangl_param_map{original.rangl_param_map},
    rangl_param_map_bounds{original.rangl_param_map_bounds},
    rdihe_param_map{original.rdihe_param_map},
    rdihe_param_map_bounds{original.rdihe_param_map_bounds},
    nmr_int_data{original.nmr_int_data},
    nmr_int2_data{original.nmr_int2_data},

    // Virtual site parameter arrays, term atom indices, and parameter indices
    virtual_site_parameters{original.virtual_site_parameters},
    sp_virtual_site_parameters{original.sp_virtual_site_parameters},
    virtual_site_atoms{original.virtual_site_atoms},
    virtual_site_frame1_atoms{original.virtual_site_frame1_atoms},
    virtual_site_frame2_atoms{original.virtual_site_frame2_atoms},
    virtual_site_frame3_atoms{original.virtual_site_frame3_atoms},
    virtual_site_frame4_atoms{original.virtual_site_frame4_atoms},
    virtual_site_parameter_indices{original.virtual_site_parameter_indices},
    vsite_int_data{original.vsite_int_data},

    // SETTLE (analytic rigid water constraints) parameter arrays, plus fused atom and parameter
    // indices.
    settle_group_geometry{original.settle_group_geometry},
    settle_group_masses{original.settle_group_masses},
    sp_settle_group_geometry{original.sp_settle_group_geometry},
    sp_settle_group_masses{original.sp_settle_group_masses},
    settle_group_indexing{original.settle_group_indexing},

    // Constraint group parameter indices and bounds (referencing consensus tables from the
    // synthesis) plus length / inverse mass parameters (also fused)
    constraint_group_indices{original.constraint_group_indices},
    constraint_group_bounds{original.constraint_group_bounds},
    constraint_group_param_idx{original.constraint_group_param_idx},
    constraint_param_bounds{original.constraint_param_bounds},
    constraint_group_params{original.constraint_group_params},
    sp_constraint_group_params{original.sp_constraint_group_params},

    // Valence work unit instruction sets and energy accumulation masks
    total_valence_work_units{original.total_valence_work_units},
    valence_work_unit_size{original.valence_work_unit_size},
    valence_thread_block_size{original.valence_thread_block_size},
    vwu_instruction_sets{original.vwu_instruction_sets},
    vwu_import_lists{original.vwu_import_lists},
    vwu_manipulation_masks{original.vwu_manipulation_masks},
    cbnd_instructions{original.cbnd_instructions},
    angl_instructions{original.angl_instructions},
    cdhe_instructions{original.cdhe_instructions},
    cdhe_overtones{original.cdhe_overtones},
    cmap_instructions{original.cmap_instructions},
    infr14_instructions{original.infr14_instructions},
    rposn_instructions{original.rposn_instructions},
    rbond_instructions{original.rbond_instructions},
    rangl_instructions{original.rangl_instructions},
    rdihe_instructions{original.rdihe_instructions},
    vste_instructions{original.vste_instructions},
    sett_instructions{original.sett_instructions},
    cnst_instructions{original.cnst_instructions},
    accumulate_cbnd_energy{original.accumulate_cbnd_energy},
    accumulate_angl_energy{original.accumulate_angl_energy},
    accumulate_cdhe_energy{original.accumulate_cdhe_energy},
    accumulate_cmap_energy{original.accumulate_cmap_energy},
    accumulate_infr14_energy{original.accumulate_infr14_energy},
    accumulate_rposn_energy{original.accumulate_rposn_energy},
    accumulate_rbond_energy{original.accumulate_rbond_energy},
    accumulate_rangl_energy{original.accumulate_rangl_energy},
    accumulate_rdihe_energy{original.accumulate_rdihe_energy},
    insr_uint_data{original.insr_uint_data},
    insr_uint2_data{original.insr_uint2_data},

    // Non-bonded work units and their instructions
    total_nonbonded_work_units{original.total_nonbonded_work_units},
    nonbonded_work_type{original.nonbonded_work_type},
    nonbonded_abstracts{original.nonbonded_abstracts},
    nbwu_instructions{original.nbwu_instructions},

    // Reduction work units and their abstracts
    total_reduction_work_units{original.total_reduction_work_units},
    rdwu_per_system{original.rdwu_per_system},
    reduction_abstracts{original.reduction_abstracts},
    
    // Pointer to the timer used to track time needed for assembling this object.  The timer can
    // also be used to track access and other usage of the object.
    timer{original.timer}
{
  rebasePointers();
}

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis::AtomGraphSynthesis(AtomGraphSynthesis &&original) :

    // Counts spanning all topologies
    policy{original.policy},
    topology_count{original.topology_count},
    restraint_network_count{original.restraint_network_count},
    system_count{original.system_count},
    total_atoms{original.total_atoms},
    total_virtual_sites{original.total_virtual_sites},
    total_bond_terms{original.total_bond_terms},
    total_angl_terms{original.total_angl_terms},
    total_dihe_terms{original.total_dihe_terms},
    total_ubrd_terms{original.total_ubrd_terms},
    total_cimp_terms{original.total_cimp_terms},
    total_cmap_terms{original.total_cmap_terms},
    total_lj_types{original.total_lj_types},
    total_charge_types{original.total_charge_types},
    total_bond_params{original.total_bond_params},
    total_angl_params{original.total_angl_params},
    total_dihe_params{original.total_dihe_params},
    total_ubrd_params{original.total_ubrd_params},
    total_cimp_params{original.total_cimp_params},
    total_cmap_surfaces{original.total_cmap_surfaces},
    total_attn14_params{original.total_attn14_params},
    total_vste_params{original.total_vste_params},
    total_position_restraints{original.total_position_restraints},
    total_distance_restraints{original.total_distance_restraints},
    total_angle_restraints{original.total_angle_restraints},
    total_dihedral_restraints{original.total_dihedral_restraints},

    // Descriptors spanning all systems
    periodic_box_class{original.periodic_box_class},
    gb_style{original.gb_style},
    dielectric_constant{original.dielectric_constant},
    is_kappa{original.is_kappa},
    salt_concentration{original.salt_concentration},
    gb_offset{original.gb_offset},
    gb_neckscale{original.gb_neckscale},
    gb_neckcut{original.gb_neckcut},
    coulomb_constant{original.coulomb_constant},
    use_bond_constraints{original.use_bond_constraints},
    use_settle{original.use_settle},
    largest_constraint_group{original.largest_constraint_group},
    water_residue_name{original.water_residue_name},
    pb_radii_sets{std::move(original.pb_radii_sets)},

    // Individual topologies and restraint networks (RestraintApparatus objects), each of them
    // describing one or more of the individual systems within the synthesis
    topologies{std::move(original.topologies)},
    restraint_networks{std::move(original.restraint_networks)},
    restraint_dummies{std::move(original.restraint_dummies)},
    topology_indices{std::move(original.topology_indices)},
    restraint_indices{std::move(original.restraint_indices)},
    atom_counts{std::move(original.atom_counts)},
    residue_counts{std::move(original.residue_counts)},
    molecule_counts{std::move(original.molecule_counts)},
    largest_residue_sizes{std::move(original.largest_residue_sizes)},
    last_solute_residues{std::move(original.last_solute_residues)},
    last_solute_atoms{std::move(original.last_solute_atoms)},
    first_solvent_molecules{std::move(original.first_solvent_molecules)},

    // Counts of energy terms and particle-related quantities
    ubrd_term_counts{std::move(original.ubrd_term_counts)},
    cimp_term_counts{std::move(original.cimp_term_counts)},
    cmap_term_counts{std::move(original.cmap_term_counts)},
    bond_term_counts{std::move(original.bond_term_counts)},
    angl_term_counts{std::move(original.angl_term_counts)},
    dihe_term_counts{std::move(original.dihe_term_counts)},
    virtual_site_counts{std::move(original.virtual_site_counts)},
    posn_restraint_counts{std::move(original.posn_restraint_counts)},
    bond_restraint_counts{std::move(original.bond_restraint_counts)},
    angl_restraint_counts{std::move(original.angl_restraint_counts)},
    dihe_restraint_counts{std::move(original.dihe_restraint_counts)},
    lj_type_counts{std::move(original.lj_type_counts)},
    total_exclusion_counts{std::move(original.total_exclusion_counts)},
    rigid_water_counts{std::move(original.rigid_water_counts)},
    bond_constraint_counts{std::move(original.bond_constraint_counts)},
    degrees_of_freedom{std::move(original.degrees_of_freedom)},
    cnst_degrees_of_freedom{std::move(original.cnst_degrees_of_freedom)},
    nonrigid_particle_counts{std::move(original.nonrigid_particle_counts)},

    // Offsets for each system's term lists
    atom_offsets{std::move(original.atom_offsets)},
    atom_bit_offsets{std::move(original.atom_bit_offsets)},
    residue_offsets{std::move(original.residue_offsets)},
    molecule_offsets{std::move(original.molecule_offsets)},
    ubrd_term_offsets{std::move(original.ubrd_term_offsets)},
    cimp_term_offsets{std::move(original.cimp_term_offsets)},
    cmap_term_offsets{std::move(original.cmap_term_offsets)},
    bond_term_offsets{std::move(original.bond_term_offsets)},
    angl_term_offsets{std::move(original.angl_term_offsets)},
    dihe_term_offsets{std::move(original.dihe_term_offsets)},
    virtual_site_offsets{std::move(original.virtual_site_offsets)},
    posn_restraint_offsets{std::move(original.posn_restraint_offsets)},
    bond_restraint_offsets{std::move(original.bond_restraint_offsets)},
    angl_restraint_offsets{std::move(original.angl_restraint_offsets)},
    dihe_restraint_offsets{std::move(original.dihe_restraint_offsets)},
    sett_group_offsets{std::move(original.sett_group_offsets)},
    cnst_group_offsets{std::move(original.cnst_group_offsets)},
    nb_exclusion_offsets{std::move(original.nb_exclusion_offsets)},
    lennard_jones_abc_offsets{std::move(original.lennard_jones_abc_offsets)},
    int_system_data{std::move(original.int_system_data)},

    // Atom and residue details
    residue_limits{std::move(original.residue_limits)},
    atom_struc_numbers{std::move(original.atom_struc_numbers)},
    residue_numbers{std::move(original.residue_numbers)},
    molecule_limits{std::move(original.molecule_limits)},
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
    chem_int_data{std::move(original.chem_int_data)},
    chem_int2_data{std::move(original.chem_int2_data)},
    chem_double_data{std::move(original.chem_double_data)},
    chem_float_data{std::move(original.chem_float_data)},
    chem_char4_data{std::move(original.chem_char4_data)},

    // Force field parameter maps from individual systems to the consensus tables
    ubrd_param_map{std::move(original.ubrd_param_map)},
    ubrd_param_map_bounds{std::move(original.ubrd_param_map_bounds)},
    cimp_param_map{std::move(original.cimp_param_map)},
    cimp_param_map_bounds{std::move(original.cimp_param_map_bounds)},
    cmap_param_map{std::move(original.cmap_param_map)},
    cmap_param_map_bounds{std::move(original.cmap_param_map_bounds)},
    bond_param_map{std::move(original.bond_param_map)},
    bond_param_map_bounds{std::move(original.bond_param_map_bounds)},
    angl_param_map{std::move(original.angl_param_map)},
    angl_param_map_bounds{std::move(original.angl_param_map_bounds)},
    dihe_param_map{std::move(original.dihe_param_map)},
    dihe_param_map_bounds{std::move(original.dihe_param_map_bounds)},
    attn14_param_map{std::move(original.attn14_param_map)},
    attn14_param_map_bounds{std::move(original.attn14_param_map_bounds)},
    vste_param_map{std::move(original.vste_param_map)},
    vste_param_map_bounds{std::move(original.vste_param_map_bounds)},
    sett_param_map{std::move(original.sett_param_map)},
    sett_param_map_bounds{std::move(original.sett_param_map_bounds)},
    cnst_param_map{std::move(original.cnst_param_map)},
    cnst_param_map_bounds{std::move(original.cnst_param_map_bounds)},

    // Force field consensus table condensed parameters
    ubrd_stiffnesses{std::move(original.ubrd_stiffnesses)},
    ubrd_equilibria{std::move(original.ubrd_equilibria)},
    cimp_stiffnesses{std::move(original.cimp_stiffnesses)},
    cimp_phase_angles{std::move(original.cimp_phase_angles)},
    cmap_surface_dimensions{std::move(original.cmap_surface_dimensions)},
    cmap_surface_bounds{std::move(original.cmap_surface_bounds)},
    cmap_patch_bounds{std::move(original.cmap_patch_bounds)},
    cmap_surfaces{std::move(original.cmap_surfaces)},
    cmap_patches{std::move(original.cmap_patches)},
    sp_ubrd_stiffnesses{std::move(original.sp_ubrd_stiffnesses)},
    sp_ubrd_equilibria{std::move(original.sp_ubrd_equilibria)},
    sp_cimp_stiffnesses{std::move(original.sp_cimp_stiffnesses)},
    sp_cimp_phase_angles{std::move(original.sp_cimp_phase_angles)},
    sp_cmap_surfaces{std::move(original.sp_cmap_surfaces)},
    sp_cmap_patches{std::move(original.sp_cmap_patches)},
    bond_stiffnesses{std::move(original.bond_stiffnesses)},
    bond_equilibria{std::move(original.bond_equilibria)},
    angl_stiffnesses{std::move(original.angl_stiffnesses)},
    angl_equilibria{std::move(original.angl_equilibria)},
    dihe_amplitudes{std::move(original.dihe_amplitudes)},
    dihe_periodicities{std::move(original.dihe_periodicities)},
    dihe_phase_angles{std::move(original.dihe_phase_angles)},
    attn14_elec_factors{std::move(original.attn14_elec_factors)},
    attn14_vdw_factors{std::move(original.attn14_vdw_factors)},
    sp_bond_stiffnesses{std::move(original.sp_bond_stiffnesses)},
    sp_bond_equilibria{std::move(original.sp_bond_equilibria)},
    sp_angl_stiffnesses{std::move(original.sp_angl_stiffnesses)},
    sp_angl_equilibria{std::move(original.sp_angl_equilibria)},
    sp_dihe_amplitudes{std::move(original.sp_dihe_amplitudes)},
    sp_dihe_periodicities{std::move(original.sp_dihe_periodicities)},
    sp_dihe_phase_angles{std::move(original.sp_dihe_phase_angles)},
    sp_attn14_elec_factors{std::move(original.sp_attn14_elec_factors)},
    sp_attn14_vdw_factors{std::move(original.sp_attn14_vdw_factors)},
    valparam_double_data{std::move(original.valparam_double_data)},
    valparam_float_data{std::move(original.valparam_float_data)},
    valparam_int_data{std::move(original.valparam_int_data)},
    valparam_int2_data{std::move(original.valparam_int2_data)},

    // Atoms participating in each valence term, plus the parameter indices of each term
    ubrd_i_atoms{std::move(original.ubrd_i_atoms)},
    ubrd_k_atoms{std::move(original.ubrd_k_atoms)},
    ubrd_param_idx{std::move(original.ubrd_param_idx)},
    cimp_i_atoms{std::move(original.cimp_i_atoms)},
    cimp_j_atoms{std::move(original.cimp_j_atoms)},
    cimp_k_atoms{std::move(original.cimp_k_atoms)},
    cimp_l_atoms{std::move(original.cimp_l_atoms)},
    cimp_param_idx{std::move(original.cimp_param_idx)},
    cmap_i_atoms{std::move(original.cmap_i_atoms)},
    cmap_j_atoms{std::move(original.cmap_j_atoms)},
    cmap_k_atoms{std::move(original.cmap_k_atoms)},
    cmap_l_atoms{std::move(original.cmap_l_atoms)},
    cmap_m_atoms{std::move(original.cmap_m_atoms)},
    cmap_param_idx{std::move(original.cmap_param_idx)},
    bond_i_atoms{std::move(original.bond_i_atoms)},
    bond_j_atoms{std::move(original.bond_j_atoms)},
    bond_param_idx{std::move(original.bond_param_idx)},
    angl_i_atoms{std::move(original.angl_i_atoms)},
    angl_j_atoms{std::move(original.angl_j_atoms)},
    angl_k_atoms{std::move(original.angl_k_atoms)},
    angl_param_idx{std::move(original.angl_param_idx)},
    dihe_i_atoms{std::move(original.dihe_i_atoms)},
    dihe_j_atoms{std::move(original.dihe_j_atoms)},
    dihe_k_atoms{std::move(original.dihe_k_atoms)},
    dihe_l_atoms{std::move(original.dihe_l_atoms)},
    dihe_param_idx{std::move(original.dihe_param_idx)},
    valence_int_data{std::move(original.valence_int_data)},

    // Non-bonded parameters for all systems
    charge_indices{std::move(original.charge_indices)},
    lennard_jones_indices{std::move(original.lennard_jones_indices)},
    charge_parameters{std::move(original.charge_parameters)},
    lennard_jones_ab_coeff{std::move(original.lennard_jones_ab_coeff)},
    lennard_jones_c_coeff{std::move(original.lennard_jones_c_coeff)},
    lennard_jones_14_a_coeff{std::move(original.lennard_jones_14_a_coeff)},
    lennard_jones_14_b_coeff{std::move(original.lennard_jones_14_b_coeff)},
    lennard_jones_14_c_coeff{std::move(original.lennard_jones_14_c_coeff)},
    lennard_jones_sigma{std::move(original.lennard_jones_sigma)},
    lennard_jones_14_sigma{std::move(original.lennard_jones_14_sigma)},
    sp_charge_parameters{std::move(original.sp_charge_parameters)},
    sp_lennard_jones_ab_coeff{std::move(original.sp_lennard_jones_ab_coeff)},
    sp_lennard_jones_c_coeff{std::move(original.sp_lennard_jones_c_coeff)},
    sp_lennard_jones_14_a_coeff{std::move(original.sp_lennard_jones_14_a_coeff)},
    sp_lennard_jones_14_b_coeff{std::move(original.sp_lennard_jones_14_b_coeff)},
    sp_lennard_jones_14_c_coeff{std::move(original.sp_lennard_jones_14_c_coeff)},
    sp_lennard_jones_sigma{std::move(original.sp_lennard_jones_sigma)},
    sp_lennard_jones_14_sigma{std::move(original.sp_lennard_jones_14_sigma)},

    // Implicit solvent model parameters
    neck_table_size{original.neck_table_size},
    neck_gb_indices{std::move(original.neck_gb_indices)},
    atomic_pb_radii{std::move(original.atomic_pb_radii)},
    gb_screening_factors{std::move(original.gb_screening_factors)},
    gb_alpha_parameters{std::move(original.gb_alpha_parameters)},
    gb_beta_parameters{std::move(original.gb_beta_parameters)},
    gb_gamma_parameters{std::move(original.gb_gamma_parameters)},
    neck_limit_tables{std::move(original.neck_limit_tables)},
    sp_atomic_pb_radii{std::move(original.sp_atomic_pb_radii)},
    sp_gb_screening_factors{std::move(original.sp_gb_screening_factors)},
    sp_gb_alpha_parameters{std::move(original.sp_gb_alpha_parameters)},
    sp_gb_beta_parameters{std::move(original.sp_gb_beta_parameters)},
    sp_gb_gamma_parameters{std::move(original.sp_gb_gamma_parameters)},
    sp_neck_limit_tables{std::move(original.sp_neck_limit_tables)},
    
    // Restraint parameters for all systems
    rposn_step_bounds{std::move(original.rposn_step_bounds)},
    rbond_step_bounds{std::move(original.rbond_step_bounds)},
    rangl_step_bounds{std::move(original.rangl_step_bounds)},
    rdihe_step_bounds{std::move(original.rdihe_step_bounds)},
    rposn_init_k{std::move(original.rposn_init_k)},
    rposn_final_k{std::move(original.rposn_final_k)},
    rposn_init_r{std::move(original.rposn_init_r)},
    rposn_final_r{std::move(original.rposn_final_r)},
    rposn_init_xy{std::move(original.rposn_init_xy)},
    rposn_init_z{std::move(original.rposn_init_z)},
    rposn_final_xy{std::move(original.rposn_final_xy)},
    rposn_final_z{std::move(original.rposn_final_z)},
    rbond_init_k{std::move(original.rbond_init_k)},
    rbond_final_k{std::move(original.rbond_final_k)},
    rbond_init_r{std::move(original.rbond_init_r)},
    rbond_final_r{std::move(original.rbond_final_r)},
    rangl_init_k{std::move(original.rangl_init_k)},
    rangl_final_k{std::move(original.rangl_final_k)},
    rangl_init_r{std::move(original.rangl_init_r)},
    rangl_final_r{std::move(original.rangl_final_r)},
    rdihe_init_k{std::move(original.rdihe_init_k)},
    rdihe_final_k{std::move(original.rdihe_final_k)},
    rdihe_init_r{std::move(original.rdihe_init_r)},
    rdihe_final_r{std::move(original.rdihe_final_r)},
    sp_rposn_init_k{std::move(original.sp_rposn_init_k)},
    sp_rposn_final_k{std::move(original.sp_rposn_final_k)},
    sp_rposn_init_r{std::move(original.sp_rposn_init_r)},
    sp_rposn_final_r{std::move(original.sp_rposn_final_r)},
    sp_rposn_init_xy{std::move(original.sp_rposn_init_xy)},
    sp_rposn_init_z{std::move(original.sp_rposn_init_z)},
    sp_rposn_final_xy{std::move(original.sp_rposn_final_xy)},
    sp_rposn_final_z{std::move(original.sp_rposn_final_z)},
    sp_rbond_init_k{std::move(original.sp_rbond_init_k)},
    sp_rbond_final_k{std::move(original.sp_rbond_final_k)},
    sp_rbond_init_r{std::move(original.sp_rbond_init_r)},
    sp_rbond_final_r{std::move(original.sp_rbond_final_r)},
    sp_rangl_init_k{std::move(original.sp_rangl_init_k)},
    sp_rangl_final_k{std::move(original.sp_rangl_final_k)},
    sp_rangl_init_r{std::move(original.sp_rangl_init_r)},
    sp_rangl_final_r{std::move(original.sp_rangl_final_r)},
    sp_rdihe_init_k{std::move(original.sp_rdihe_init_k)},
    sp_rdihe_final_k{std::move(original.sp_rdihe_final_k)},
    sp_rdihe_init_r{std::move(original.sp_rdihe_init_r)},
    sp_rdihe_final_r{std::move(original.sp_rdihe_final_r)},
    nmr_double_data{std::move(original.nmr_double_data)},
    nmr_double2_data{std::move(original.nmr_double2_data)},
    nmr_double4_data{std::move(original.nmr_double4_data)},
    nmr_float_data{std::move(original.nmr_float_data)},
    nmr_float2_data{std::move(original.nmr_float2_data)},
    nmr_float4_data{std::move(original.nmr_float4_data)},

    // Restraint term atom indices and parameter indices
    rposn_atoms{std::move(original.rposn_atoms)},
    rposn_kr_param_idx{std::move(original.rposn_kr_param_idx)},
    rposn_xyz_param_idx{std::move(original.rposn_xyz_param_idx)},
    rbond_i_atoms{std::move(original.rbond_i_atoms)},
    rbond_j_atoms{std::move(original.rbond_j_atoms)},
    rbond_param_idx{std::move(original.rbond_param_idx)},
    rangl_i_atoms{std::move(original.rangl_i_atoms)},
    rangl_j_atoms{std::move(original.rangl_j_atoms)},
    rangl_k_atoms{std::move(original.rangl_k_atoms)},
    rangl_param_idx{std::move(original.rangl_param_idx)},
    rdihe_i_atoms{std::move(original.rdihe_i_atoms)},
    rdihe_j_atoms{std::move(original.rdihe_j_atoms)},
    rdihe_k_atoms{std::move(original.rdihe_k_atoms)},
    rdihe_l_atoms{std::move(original.rdihe_l_atoms)},
    rdihe_param_idx{std::move(original.rdihe_param_idx)},
    rposn_kr_param_map{std::move(original.rposn_kr_param_map)},
    rposn_xyz_param_map{std::move(original.rposn_xyz_param_map)},
    rposn_param_map_bounds{std::move(original.rposn_param_map_bounds)},
    rbond_param_map{std::move(original.rbond_param_map)},
    rbond_param_map_bounds{std::move(original.rbond_param_map_bounds)},
    rangl_param_map{std::move(original.rangl_param_map)},
    rangl_param_map_bounds{std::move(original.rangl_param_map_bounds)},
    rdihe_param_map{std::move(original.rdihe_param_map)},
    rdihe_param_map_bounds{std::move(original.rdihe_param_map_bounds)},
    nmr_int_data{std::move(original.nmr_int_data)},
    nmr_int2_data{std::move(original.nmr_int2_data)},

    // Virtual site parameter arrays, term atom indices, and parameter indices
    virtual_site_parameters{std::move(original.virtual_site_parameters)},
    sp_virtual_site_parameters{std::move(original.sp_virtual_site_parameters)},
    virtual_site_atoms{std::move(original.virtual_site_atoms)},
    virtual_site_frame1_atoms{std::move(original.virtual_site_frame1_atoms)},
    virtual_site_frame2_atoms{std::move(original.virtual_site_frame2_atoms)},
    virtual_site_frame3_atoms{std::move(original.virtual_site_frame3_atoms)},
    virtual_site_frame4_atoms{std::move(original.virtual_site_frame4_atoms)},
    virtual_site_parameter_indices{std::move(original.virtual_site_parameter_indices)},
    vsite_int_data{std::move(original.vsite_int_data)},

    // SETTLE (analytic rigid water constraints) parameter arrays, plus fused atom and parameter
    // indices.
    settle_group_geometry{std::move(original.settle_group_geometry)},
    settle_group_masses{std::move(original.settle_group_masses)},
    sp_settle_group_geometry{std::move(original.sp_settle_group_geometry)},
    sp_settle_group_masses{std::move(original.sp_settle_group_masses)},
    settle_group_indexing{std::move(original.settle_group_indexing)},

    // Constraint group parameter indices and bounds (referencing consensus tables from the
    // synthesis) plus length / inverse mass parameters (also fused)
    constraint_group_indices{std::move(original.constraint_group_indices)},
    constraint_group_bounds{std::move(original.constraint_group_bounds)},
    constraint_group_param_idx{std::move(original.constraint_group_param_idx)},
    constraint_param_bounds{std::move(original.constraint_param_bounds)},
    constraint_group_params{std::move(original.constraint_group_params)},
    sp_constraint_group_params{std::move(original.sp_constraint_group_params)},

    // Valence work unit instruction sets and energy accumulation masks
    total_valence_work_units{original.total_valence_work_units},
    valence_work_unit_size{original.valence_work_unit_size},
    valence_thread_block_size{original.valence_thread_block_size},
    vwu_instruction_sets{std::move(original.vwu_instruction_sets)},
    vwu_import_lists{std::move(original.vwu_import_lists)},
    vwu_manipulation_masks{std::move(original.vwu_manipulation_masks)},
    cbnd_instructions{std::move(original.cbnd_instructions)},
    angl_instructions{std::move(original.angl_instructions)},
    cdhe_instructions{std::move(original.cdhe_instructions)},
    cdhe_overtones{std::move(original.cdhe_overtones)},
    cmap_instructions{std::move(original.cmap_instructions)},
    infr14_instructions{std::move(original.infr14_instructions)},
    rposn_instructions{std::move(original.rposn_instructions)},
    rbond_instructions{std::move(original.rbond_instructions)},
    rangl_instructions{std::move(original.rangl_instructions)},
    rdihe_instructions{std::move(original.rdihe_instructions)},
    vste_instructions{std::move(original.vste_instructions)},
    sett_instructions{std::move(original.sett_instructions)},
    cnst_instructions{std::move(original.cnst_instructions)},
    accumulate_cbnd_energy{std::move(original.accumulate_cbnd_energy)},
    accumulate_angl_energy{std::move(original.accumulate_angl_energy)},
    accumulate_cdhe_energy{std::move(original.accumulate_cdhe_energy)},
    accumulate_cmap_energy{std::move(original.accumulate_cmap_energy)},
    accumulate_infr14_energy{std::move(original.accumulate_infr14_energy)},
    accumulate_rposn_energy{std::move(original.accumulate_rposn_energy)},
    accumulate_rbond_energy{std::move(original.accumulate_rbond_energy)},
    accumulate_rangl_energy{std::move(original.accumulate_rangl_energy)},
    accumulate_rdihe_energy{std::move(original.accumulate_rdihe_energy)},
    insr_uint_data{std::move(original.insr_uint_data)},
    insr_uint2_data{std::move(original.insr_uint2_data)},

    // Non-bonded work units and their instructions
    total_nonbonded_work_units{original.total_nonbonded_work_units},
    nonbonded_work_type{std::move(original.nonbonded_work_type)},
    nonbonded_abstracts{std::move(original.nonbonded_abstracts)},
    nbwu_instructions{std::move(original.nbwu_instructions)},

    // Reduction work units and their abstracts
    total_reduction_work_units{original.total_reduction_work_units},
    rdwu_per_system{original.rdwu_per_system},
    reduction_abstracts{std::move(original.reduction_abstracts)},
    
    // Pointer to the timer used to track time needed for assembling this object.  The timer can
    // also be used to track access and other usage of the object.
    timer{std::move(original.timer)}
{}

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis& AtomGraphSynthesis::operator=(const AtomGraphSynthesis &other) {

  // Gaurd against self-assignment
  if (this == &other) {
    return *this;
  }
  
  // Counts spanning all topologies
  policy = other.policy;
  topology_count = other.topology_count;
  restraint_network_count = other.restraint_network_count;
  system_count = other.system_count;
  total_atoms = other.total_atoms;
  total_virtual_sites = other.total_virtual_sites;
  total_bond_terms = other.total_bond_terms;
  total_angl_terms = other.total_angl_terms;
  total_dihe_terms = other.total_dihe_terms;
  total_ubrd_terms = other.total_ubrd_terms;
  total_cimp_terms = other.total_cimp_terms;
  total_cmap_terms = other.total_cmap_terms;
  total_lj_types = other.total_lj_types;
  total_charge_types = other.total_charge_types;
  total_bond_params = other.total_bond_params;
  total_angl_params = other.total_angl_params;
  total_dihe_params = other.total_dihe_params;
  total_ubrd_params = other.total_ubrd_params;
  total_cimp_params = other.total_cimp_params;
  total_cmap_surfaces = other.total_cmap_surfaces;
  total_attn14_params = other.total_attn14_params;
  total_vste_params = other.total_vste_params;
  total_position_restraints = other.total_position_restraints;
  total_distance_restraints = other.total_distance_restraints;
  total_angle_restraints = other.total_angle_restraints;
  total_dihedral_restraints = other.total_dihedral_restraints;

  // Descriptors spanning all systems
  periodic_box_class = other.periodic_box_class;
  gb_style = other.gb_style;
  dielectric_constant = other.dielectric_constant;
  is_kappa = other.is_kappa;
  salt_concentration = other.salt_concentration;
  gb_offset = other.gb_offset;
  gb_neckscale = other.gb_neckscale;
  gb_neckcut = other.gb_neckcut;
  coulomb_constant = other.coulomb_constant;
  use_bond_constraints = other.use_bond_constraints;
  use_settle = other.use_settle;
  largest_constraint_group = other.largest_constraint_group;
  water_residue_name = other.water_residue_name;
  pb_radii_sets = other.pb_radii_sets;

  // Individual topologies and restraint networks (RestraintApparatus objects), each of them
  // describing one or more of the individual systems within the synthesis
  topologies = other.topologies;
  restraint_networks = other.restraint_networks;
  restraint_dummies = other.restraint_dummies;
  topology_indices = other.topology_indices;
  restraint_indices = other.restraint_indices;
  atom_counts = other.atom_counts;
  residue_counts = other.residue_counts;
  molecule_counts = other.molecule_counts;
  largest_residue_sizes = other.largest_residue_sizes;
  last_solute_residues = other.last_solute_residues;
  last_solute_atoms = other.last_solute_atoms;
  first_solvent_molecules = other.first_solvent_molecules;

  // Counts of energy terms and particle-related quantities
  ubrd_term_counts = other.ubrd_term_counts;
  cimp_term_counts = other.cimp_term_counts;
  cmap_term_counts = other.cmap_term_counts;
  bond_term_counts = other.bond_term_counts;
  angl_term_counts = other.angl_term_counts;
  dihe_term_counts = other.dihe_term_counts;
  virtual_site_counts = other.virtual_site_counts;
  posn_restraint_counts = other.posn_restraint_counts;
  bond_restraint_counts = other.bond_restraint_counts;
  angl_restraint_counts = other.angl_restraint_counts;
  dihe_restraint_counts = other.dihe_restraint_counts;
  lj_type_counts = other.lj_type_counts;
  total_exclusion_counts = other.total_exclusion_counts;
  rigid_water_counts = other.rigid_water_counts;
  bond_constraint_counts = other.bond_constraint_counts;
  degrees_of_freedom = other.degrees_of_freedom;
  cnst_degrees_of_freedom = other.cnst_degrees_of_freedom;
  nonrigid_particle_counts = other.nonrigid_particle_counts;

  // Offsets for each system's term lists
  atom_offsets = other.atom_offsets;
  atom_bit_offsets = other.atom_bit_offsets;
  residue_offsets = other.residue_offsets;
  molecule_offsets = other.molecule_offsets;
  ubrd_term_offsets = other.ubrd_term_offsets;
  cimp_term_offsets = other.cimp_term_offsets;
  cmap_term_offsets = other.cmap_term_offsets;
  bond_term_offsets = other.bond_term_offsets;
  angl_term_offsets = other.angl_term_offsets;
  dihe_term_offsets = other.dihe_term_offsets;
  virtual_site_offsets = other.virtual_site_offsets;
  posn_restraint_offsets = other.posn_restraint_offsets;
  bond_restraint_offsets = other.bond_restraint_offsets;
  angl_restraint_offsets = other.angl_restraint_offsets;
  dihe_restraint_offsets = other.dihe_restraint_offsets;
  sett_group_offsets = other.sett_group_offsets;
  cnst_group_offsets = other.cnst_group_offsets;
  nb_exclusion_offsets = other.nb_exclusion_offsets;
  lennard_jones_abc_offsets = other.lennard_jones_abc_offsets;
  int_system_data = other.int_system_data;

  // Atom and residue details
  residue_limits = other.residue_limits;
  atom_struc_numbers = other.atom_struc_numbers;
  residue_numbers = other.residue_numbers;
  molecule_limits = other.molecule_limits;
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
  chem_int_data = other.chem_int_data;
  chem_int2_data = other.chem_int2_data;
  chem_double_data = other.chem_double_data;
  chem_float_data = other.chem_float_data;
  chem_char4_data = other.chem_char4_data;

  // Force field parameter maps from individual systems to the consensus tables
  ubrd_param_map = other.ubrd_param_map;
  ubrd_param_map_bounds = other.ubrd_param_map_bounds;
  cimp_param_map = other.cimp_param_map;
  cimp_param_map_bounds = other.cimp_param_map_bounds;
  cmap_param_map = other.cmap_param_map;
  cmap_param_map_bounds = other.cmap_param_map_bounds;
  bond_param_map = other.bond_param_map;
  bond_param_map_bounds = other.bond_param_map_bounds;
  angl_param_map = other.angl_param_map;
  angl_param_map_bounds = other.angl_param_map_bounds;
  dihe_param_map = other.dihe_param_map;
  dihe_param_map_bounds = other.dihe_param_map_bounds;
  attn14_param_map = other.attn14_param_map;
  attn14_param_map_bounds = other.attn14_param_map_bounds;
  vste_param_map = other.vste_param_map;
  vste_param_map_bounds = other.vste_param_map_bounds;
  sett_param_map = other.sett_param_map;
  sett_param_map_bounds = other.sett_param_map_bounds;
  cnst_param_map = other.cnst_param_map;
  cnst_param_map_bounds = other.cnst_param_map_bounds;

  // Force field consensus table condensed parameters
  ubrd_stiffnesses = other.ubrd_stiffnesses;
  ubrd_equilibria = other.ubrd_equilibria;
  cimp_stiffnesses = other.cimp_stiffnesses;
  cimp_phase_angles = other.cimp_phase_angles;
  cmap_surface_dimensions = other.cmap_surface_dimensions;
  cmap_surface_bounds = other.cmap_surface_bounds;
  cmap_patch_bounds = other.cmap_patch_bounds;
  cmap_surfaces = other.cmap_surfaces;
  cmap_patches = other.cmap_patches;
  sp_ubrd_stiffnesses = other.sp_ubrd_stiffnesses;
  sp_ubrd_equilibria = other.sp_ubrd_equilibria;
  sp_cimp_stiffnesses = other.sp_cimp_stiffnesses;
  sp_cimp_phase_angles = other.sp_cimp_phase_angles;
  sp_cmap_surfaces = other.sp_cmap_surfaces;
  sp_cmap_patches = other.sp_cmap_patches;
  bond_stiffnesses = other.bond_stiffnesses;
  bond_equilibria = other.bond_equilibria;
  angl_stiffnesses = other.angl_stiffnesses;
  angl_equilibria = other.angl_equilibria;
  dihe_amplitudes = other.dihe_amplitudes;
  dihe_periodicities = other.dihe_periodicities;
  dihe_phase_angles = other.dihe_phase_angles;
  attn14_elec_factors = other.attn14_elec_factors;
  attn14_vdw_factors = other.attn14_vdw_factors;
  sp_bond_stiffnesses = other.sp_bond_stiffnesses;
  sp_bond_equilibria = other.sp_bond_equilibria;
  sp_angl_stiffnesses = other.sp_angl_stiffnesses;
  sp_angl_equilibria = other.sp_angl_equilibria;
  sp_dihe_amplitudes = other.sp_dihe_amplitudes;
  sp_dihe_periodicities = other.sp_dihe_periodicities;
  sp_dihe_phase_angles = other.sp_dihe_phase_angles;
  sp_attn14_elec_factors = other.sp_attn14_elec_factors;
  sp_attn14_vdw_factors = other.sp_attn14_vdw_factors;
  valparam_double_data = other.valparam_double_data;
  valparam_float_data = other.valparam_float_data;
  valparam_int_data = other.valparam_int_data;
  valparam_int2_data = other.valparam_int2_data;

  // Atoms participating in each valence term, plus the parameter indices of each term
  ubrd_i_atoms = other.ubrd_i_atoms;
  ubrd_k_atoms = other.ubrd_k_atoms;
  ubrd_param_idx = other.ubrd_param_idx;
  cimp_i_atoms = other.cimp_i_atoms;
  cimp_j_atoms = other.cimp_j_atoms;
  cimp_k_atoms = other.cimp_k_atoms;
  cimp_l_atoms = other.cimp_l_atoms;
  cimp_param_idx = other.cimp_param_idx;
  cmap_i_atoms = other.cmap_i_atoms;
  cmap_j_atoms = other.cmap_j_atoms;
  cmap_k_atoms = other.cmap_k_atoms;
  cmap_l_atoms = other.cmap_l_atoms;
  cmap_m_atoms = other.cmap_m_atoms;
  cmap_param_idx = other.cmap_param_idx;
  bond_i_atoms = other.bond_i_atoms;
  bond_j_atoms = other.bond_j_atoms;
  bond_param_idx = other.bond_param_idx;
  angl_i_atoms = other.angl_i_atoms;
  angl_j_atoms = other.angl_j_atoms;
  angl_k_atoms = other.angl_k_atoms;
  angl_param_idx = other.angl_param_idx;
  dihe_i_atoms = other.dihe_i_atoms;
  dihe_j_atoms = other.dihe_j_atoms;
  dihe_k_atoms = other.dihe_k_atoms;
  dihe_l_atoms = other.dihe_l_atoms;
  dihe_param_idx = other.dihe_param_idx;
  valence_int_data = other.valence_int_data;

  // Non-bonded parameters for all systems
  charge_indices = other.charge_indices;
  lennard_jones_indices = other.lennard_jones_indices;
  charge_parameters = other.charge_parameters;
  lennard_jones_ab_coeff = other.lennard_jones_ab_coeff;
  lennard_jones_c_coeff = other.lennard_jones_c_coeff;
  lennard_jones_14_a_coeff = other.lennard_jones_14_a_coeff;
  lennard_jones_14_b_coeff = other.lennard_jones_14_b_coeff;
  lennard_jones_14_c_coeff = other.lennard_jones_14_c_coeff;
  lennard_jones_sigma = other.lennard_jones_sigma;
  lennard_jones_14_sigma = other.lennard_jones_14_sigma;
  sp_charge_parameters = other.sp_charge_parameters;
  sp_lennard_jones_ab_coeff = other.sp_lennard_jones_ab_coeff;
  sp_lennard_jones_c_coeff = other.sp_lennard_jones_c_coeff;
  sp_lennard_jones_14_a_coeff = other.sp_lennard_jones_14_a_coeff;
  sp_lennard_jones_14_b_coeff = other.sp_lennard_jones_14_b_coeff;
  sp_lennard_jones_14_c_coeff = other.sp_lennard_jones_14_c_coeff;
  sp_lennard_jones_sigma = other.sp_lennard_jones_sigma;
  sp_lennard_jones_14_sigma = other.sp_lennard_jones_14_sigma;

  // Implicit solvent model parameters
  neck_table_size = other.neck_table_size;
  neck_gb_indices = other.neck_gb_indices;
  atomic_pb_radii = other.atomic_pb_radii;
  gb_screening_factors = other.gb_screening_factors;
  gb_alpha_parameters = other.gb_alpha_parameters;
  gb_beta_parameters = other.gb_beta_parameters;
  gb_gamma_parameters = other.gb_gamma_parameters;
  neck_limit_tables = other.neck_limit_tables;
  sp_atomic_pb_radii = other.sp_atomic_pb_radii;
  sp_gb_screening_factors = other.sp_gb_screening_factors;
  sp_gb_alpha_parameters = other.sp_gb_alpha_parameters;
  sp_gb_beta_parameters = other.sp_gb_beta_parameters;
  sp_gb_gamma_parameters = other.sp_gb_gamma_parameters;
  sp_neck_limit_tables = other.sp_neck_limit_tables;
    
  // Restraint parameters for all systems
  rposn_step_bounds = other.rposn_step_bounds;
  rbond_step_bounds = other.rbond_step_bounds;
  rangl_step_bounds = other.rangl_step_bounds;
  rdihe_step_bounds = other.rdihe_step_bounds;
  rposn_init_k = other.rposn_init_k;
  rposn_final_k = other.rposn_final_k;
  rposn_init_r = other.rposn_init_r;
  rposn_final_r = other.rposn_final_r;
  rposn_init_xy = other.rposn_init_xy;
  rposn_init_z = other.rposn_init_z;
  rposn_final_xy = other.rposn_final_xy;
  rposn_final_z = other.rposn_final_z;
  rbond_init_k = other.rbond_init_k;
  rbond_final_k = other.rbond_final_k;
  rbond_init_r = other.rbond_init_r;
  rbond_final_r = other.rbond_final_r;
  rangl_init_k = other.rangl_init_k;
  rangl_final_k = other.rangl_final_k;
  rangl_init_r = other.rangl_init_r;
  rangl_final_r = other.rangl_final_r;
  rdihe_init_k = other.rdihe_init_k;
  rdihe_final_k = other.rdihe_final_k;
  rdihe_init_r = other.rdihe_init_r;
  rdihe_final_r = other.rdihe_final_r;
  sp_rposn_init_k = other.sp_rposn_init_k;
  sp_rposn_final_k = other.sp_rposn_final_k;
  sp_rposn_init_r = other.sp_rposn_init_r;
  sp_rposn_final_r = other.sp_rposn_final_r;
  sp_rposn_init_xy = other.sp_rposn_init_xy;
  sp_rposn_init_z = other.sp_rposn_init_z;
  sp_rposn_final_xy = other.sp_rposn_final_xy;
  sp_rposn_final_z = other.sp_rposn_final_z;
  sp_rbond_init_k = other.sp_rbond_init_k;
  sp_rbond_final_k = other.sp_rbond_final_k;
  sp_rbond_init_r = other.sp_rbond_init_r;
  sp_rbond_final_r = other.sp_rbond_final_r;
  sp_rangl_init_k = other.sp_rangl_init_k;
  sp_rangl_final_k = other.sp_rangl_final_k;
  sp_rangl_init_r = other.sp_rangl_init_r;
  sp_rangl_final_r = other.sp_rangl_final_r;
  sp_rdihe_init_k = other.sp_rdihe_init_k;
  sp_rdihe_final_k = other.sp_rdihe_final_k;
  sp_rdihe_init_r = other.sp_rdihe_init_r;
  sp_rdihe_final_r = other.sp_rdihe_final_r;
  nmr_double_data = other.nmr_double_data;
  nmr_double2_data = other.nmr_double2_data;
  nmr_double4_data = other.nmr_double4_data;
  nmr_float_data = other.nmr_float_data;
  nmr_float2_data = other.nmr_float2_data;
  nmr_float4_data = other.nmr_float4_data;

  // Restraint term atom indices and parameter indices.  Like standard force field terms, the
  // atom indices refer to the concatenated list of all systems and the parameter indices refer
  // to the consensus tables.
  rposn_atoms = other.rposn_atoms;
  rposn_kr_param_idx = other.rposn_kr_param_idx;
  rposn_xyz_param_idx = other.rposn_xyz_param_idx;
  rbond_i_atoms = other.rbond_i_atoms;
  rbond_j_atoms = other.rbond_j_atoms;
  rbond_param_idx = other.rbond_param_idx;
  rangl_i_atoms = other.rangl_i_atoms;
  rangl_j_atoms = other.rangl_j_atoms;
  rangl_k_atoms = other.rangl_k_atoms;
  rangl_param_idx = other.rangl_param_idx;
  rdihe_i_atoms = other.rdihe_i_atoms;
  rdihe_j_atoms = other.rdihe_j_atoms;
  rdihe_k_atoms = other.rdihe_k_atoms;
  rdihe_l_atoms = other.rdihe_l_atoms;
  rdihe_param_idx = other.rdihe_param_idx;
  rposn_kr_param_map = other.rposn_kr_param_map;
  rposn_xyz_param_map = other.rposn_xyz_param_map;
  rposn_param_map_bounds = other.rposn_param_map_bounds;
  rbond_param_map = other.rbond_param_map;
  rbond_param_map_bounds = other.rbond_param_map_bounds;
  rangl_param_map = other.rangl_param_map;
  rangl_param_map_bounds = other.rangl_param_map_bounds;
  rdihe_param_map = other.rdihe_param_map;
  rdihe_param_map_bounds = other.rdihe_param_map_bounds;
  nmr_int_data = other.nmr_int_data;
  nmr_int2_data = other.nmr_int2_data;

  // Virtual site parameter arrays, term atom indices, and parameter indices
  virtual_site_parameters = other.virtual_site_parameters;
  sp_virtual_site_parameters = other.sp_virtual_site_parameters;
  virtual_site_atoms = other.virtual_site_atoms;
  virtual_site_frame1_atoms = other.virtual_site_frame1_atoms;
  virtual_site_frame2_atoms = other.virtual_site_frame2_atoms;
  virtual_site_frame3_atoms = other.virtual_site_frame3_atoms;
  virtual_site_frame4_atoms = other.virtual_site_frame4_atoms;
  virtual_site_parameter_indices = other.virtual_site_parameter_indices;
  vsite_int_data = other.vsite_int_data;

  // SETTLE (analytic rigid water constraints) parameter arrays, plus fused atom and parameter
  // indices.
  settle_group_geometry = other.settle_group_geometry;
  settle_group_masses = other.settle_group_masses;
  sp_settle_group_geometry = other.sp_settle_group_geometry;
  sp_settle_group_masses = other.sp_settle_group_masses;
  settle_group_indexing = other.settle_group_indexing;

  // Constraint group parameter indices and bounds (referencing consensus tables from the
  // synthesis) plus length / inverse mass parameters (also fused)
  constraint_group_indices = other.constraint_group_indices;
  constraint_group_bounds = other.constraint_group_bounds;
  constraint_group_param_idx = other.constraint_group_param_idx;
  constraint_param_bounds = other.constraint_param_bounds;
  constraint_group_params = other.constraint_group_params;
  sp_constraint_group_params = other.sp_constraint_group_params;

  // Valence work unit instruction sets and energy accumulation masks
  total_valence_work_units = other.total_valence_work_units;
  valence_work_unit_size = other.valence_work_unit_size;
  valence_thread_block_size = other.valence_thread_block_size;
  vwu_instruction_sets = other.vwu_instruction_sets;
  vwu_import_lists = other.vwu_import_lists;
  vwu_manipulation_masks = other.vwu_manipulation_masks;
  cbnd_instructions = other.cbnd_instructions;
  angl_instructions = other.angl_instructions;
  cdhe_instructions = other.cdhe_instructions;
  cdhe_overtones = other.cdhe_overtones;
  cmap_instructions = other.cmap_instructions;
  infr14_instructions = other.infr14_instructions;
  rposn_instructions = other.rposn_instructions;
  rbond_instructions = other.rbond_instructions;
  rangl_instructions = other.rangl_instructions;
  rdihe_instructions = other.rdihe_instructions;
  vste_instructions = other.vste_instructions;
  sett_instructions = other.sett_instructions;
  cnst_instructions = other.cnst_instructions;
  accumulate_cbnd_energy = other.accumulate_cbnd_energy;
  accumulate_angl_energy = other.accumulate_angl_energy;
  accumulate_cdhe_energy = other.accumulate_cdhe_energy;
  accumulate_cmap_energy = other.accumulate_cmap_energy;
  accumulate_infr14_energy = other.accumulate_infr14_energy;
  accumulate_rposn_energy = other.accumulate_rposn_energy;
  accumulate_rbond_energy = other.accumulate_rbond_energy;
  accumulate_rangl_energy = other.accumulate_rangl_energy;
  accumulate_rdihe_energy = other.accumulate_rdihe_energy;
  insr_uint_data = other.insr_uint_data;
  insr_uint2_data = other.insr_uint2_data;

  // Non-bonded work units and their instructions
  total_nonbonded_work_units = other.total_nonbonded_work_units;
  nonbonded_work_type = other.nonbonded_work_type;
  nonbonded_abstracts = other.nonbonded_abstracts;
  nbwu_instructions = other.nbwu_instructions;

  // Reduction work units and their abstracts
  total_reduction_work_units = other.total_reduction_work_units;
  rdwu_per_system = other.rdwu_per_system;
  reduction_abstracts = other.reduction_abstracts;
  
  // Pointer to the timer used to track time needed for assembling this object.  The timer can
  // also be used to track access and other usage of the object.
  timer = other.timer;

  // Repair all pointers
  rebasePointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis& AtomGraphSynthesis::operator=(AtomGraphSynthesis &&other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }
  
  // Counts spanning all topologies
  policy = other.policy;
  topology_count = other.topology_count;
  restraint_network_count = other.restraint_network_count;
  system_count = other.system_count;
  total_atoms = other.total_atoms;
  total_virtual_sites = other.total_virtual_sites;
  total_bond_terms = other.total_bond_terms;
  total_angl_terms = other.total_angl_terms;
  total_dihe_terms = other.total_dihe_terms;
  total_ubrd_terms = other.total_ubrd_terms;
  total_cimp_terms = other.total_cimp_terms;
  total_cmap_terms = other.total_cmap_terms;
  total_lj_types = other.total_lj_types;
  total_charge_types = other.total_charge_types;
  total_bond_params = other.total_bond_params;
  total_angl_params = other.total_angl_params;
  total_dihe_params = other.total_dihe_params;
  total_ubrd_params = other.total_ubrd_params;
  total_cimp_params = other.total_cimp_params;
  total_cmap_surfaces = other.total_cmap_surfaces;
  total_attn14_params = other.total_attn14_params;
  total_vste_params = other.total_vste_params;
  total_position_restraints = other.total_position_restraints;
  total_distance_restraints = other.total_distance_restraints;
  total_angle_restraints = other.total_angle_restraints;
  total_dihedral_restraints = other.total_dihedral_restraints;

  // Descriptors spanning all systems
  periodic_box_class = other.periodic_box_class;
  gb_style = other.gb_style;
  dielectric_constant = other.dielectric_constant;
  is_kappa = other.is_kappa;
  salt_concentration = other.salt_concentration;
  gb_offset = other.gb_offset;
  gb_neckscale = other.gb_neckscale;
  gb_neckcut = other.gb_neckcut;
  coulomb_constant = other.coulomb_constant;
  use_bond_constraints = other.use_bond_constraints;
  use_settle = other.use_settle;
  largest_constraint_group = other.largest_constraint_group;
  water_residue_name = other.water_residue_name;
  pb_radii_sets = std::move(other.pb_radii_sets);

  // Individual topologies and restraint networks (RestraintApparatus objects), each of them
  // describing one or more of the individual systems within the synthesis
  topologies = std::move(other.topologies);
  restraint_networks = std::move(other.restraint_networks);
  restraint_dummies = std::move(other.restraint_dummies);
  topology_indices = std::move(other.topology_indices);
  restraint_indices = std::move(other.restraint_indices);
  atom_counts = std::move(other.atom_counts);
  residue_counts = std::move(other.residue_counts);
  molecule_counts = std::move(other.molecule_counts);
  largest_residue_sizes = std::move(other.largest_residue_sizes);
  last_solute_residues = std::move(other.last_solute_residues);
  last_solute_atoms = std::move(other.last_solute_atoms);
  first_solvent_molecules = std::move(other.first_solvent_molecules);

  // Counts of energy terms and particle-related quantities
  ubrd_term_counts = std::move(other.ubrd_term_counts);
  cimp_term_counts = std::move(other.cimp_term_counts);
  cmap_term_counts = std::move(other.cmap_term_counts);
  bond_term_counts = std::move(other.bond_term_counts);
  angl_term_counts = std::move(other.angl_term_counts);
  dihe_term_counts = std::move(other.dihe_term_counts);
  virtual_site_counts = std::move(other.virtual_site_counts);
  posn_restraint_counts = std::move(other.posn_restraint_counts);
  bond_restraint_counts = std::move(other.bond_restraint_counts);
  angl_restraint_counts = std::move(other.angl_restraint_counts);
  dihe_restraint_counts = std::move(other.dihe_restraint_counts);
  lj_type_counts = std::move(other.lj_type_counts);
  total_exclusion_counts = std::move(other.total_exclusion_counts);
  rigid_water_counts = std::move(other.rigid_water_counts);
  bond_constraint_counts = std::move(other.bond_constraint_counts);
  degrees_of_freedom = std::move(other.degrees_of_freedom);
  cnst_degrees_of_freedom = std::move(other.cnst_degrees_of_freedom);
  nonrigid_particle_counts = std::move(other.nonrigid_particle_counts);

  // Offsets for each system's term lists
  atom_offsets = std::move(other.atom_offsets);
  atom_bit_offsets = std::move(other.atom_bit_offsets);
  residue_offsets = std::move(other.residue_offsets);
  molecule_offsets = std::move(other.molecule_offsets);
  ubrd_term_offsets = std::move(other.ubrd_term_offsets);
  cimp_term_offsets = std::move(other.cimp_term_offsets);
  cmap_term_offsets = std::move(other.cmap_term_offsets);
  bond_term_offsets = std::move(other.bond_term_offsets);
  angl_term_offsets = std::move(other.angl_term_offsets);
  dihe_term_offsets = std::move(other.dihe_term_offsets);
  virtual_site_offsets = std::move(other.virtual_site_offsets);
  posn_restraint_offsets = std::move(other.posn_restraint_offsets);
  bond_restraint_offsets = std::move(other.bond_restraint_offsets);
  angl_restraint_offsets = std::move(other.angl_restraint_offsets);
  dihe_restraint_offsets = std::move(other.dihe_restraint_offsets);
  sett_group_offsets = std::move(other.sett_group_offsets);
  cnst_group_offsets = std::move(other.cnst_group_offsets);
  nb_exclusion_offsets = std::move(other.nb_exclusion_offsets);
  lennard_jones_abc_offsets = std::move(other.lennard_jones_abc_offsets);
  int_system_data = std::move(other.int_system_data);

  // Atom and residue details
  residue_limits = std::move(other.residue_limits);
  atom_struc_numbers = std::move(other.atom_struc_numbers);
  residue_numbers = std::move(other.residue_numbers);
  molecule_limits = std::move(other.molecule_limits);
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
  chem_int_data = std::move(other.chem_int_data);
  chem_int2_data = std::move(other.chem_int2_data);
  chem_double_data = std::move(other.chem_double_data);
  chem_float_data = std::move(other.chem_float_data);
  chem_char4_data = std::move(other.chem_char4_data);

  // Force field parameter maps from individual systems to the consensus tables
  ubrd_param_map = std::move(other.ubrd_param_map);
  ubrd_param_map_bounds = std::move(other.ubrd_param_map_bounds);
  cimp_param_map = std::move(other.cimp_param_map);
  cimp_param_map_bounds = std::move(other.cimp_param_map_bounds);
  cmap_param_map = std::move(other.cmap_param_map);
  cmap_param_map_bounds = std::move(other.cmap_param_map_bounds);
  bond_param_map = std::move(other.bond_param_map);
  bond_param_map_bounds = std::move(other.bond_param_map_bounds);
  angl_param_map = std::move(other.angl_param_map);
  angl_param_map_bounds = std::move(other.angl_param_map_bounds);
  dihe_param_map = std::move(other.dihe_param_map);
  dihe_param_map_bounds = std::move(other.dihe_param_map_bounds);
  attn14_param_map = std::move(other.attn14_param_map);
  attn14_param_map_bounds = std::move(other.attn14_param_map_bounds);
  vste_param_map = std::move(other.vste_param_map);
  vste_param_map_bounds = std::move(other.vste_param_map_bounds);
  sett_param_map = std::move(other.sett_param_map);
  sett_param_map_bounds = std::move(other.sett_param_map_bounds);
  cnst_param_map = std::move(other.cnst_param_map);
  cnst_param_map_bounds = std::move(other.cnst_param_map_bounds);

  // Force field consensus table condensed parameters
  ubrd_stiffnesses = std::move(other.ubrd_stiffnesses);
  ubrd_equilibria = std::move(other.ubrd_equilibria);
  cimp_stiffnesses = std::move(other.cimp_stiffnesses);
  cimp_phase_angles = std::move(other.cimp_phase_angles);
  cmap_surface_dimensions = std::move(other.cmap_surface_dimensions);
  cmap_surface_bounds = std::move(other.cmap_surface_bounds);
  cmap_patch_bounds = std::move(other.cmap_patch_bounds);
  cmap_surfaces = std::move(other.cmap_surfaces);
  cmap_patches = std::move(other.cmap_patches);
  sp_ubrd_stiffnesses = std::move(other.sp_ubrd_stiffnesses);
  sp_ubrd_equilibria = std::move(other.sp_ubrd_equilibria);
  sp_cimp_stiffnesses = std::move(other.sp_cimp_stiffnesses);
  sp_cimp_phase_angles = std::move(other.sp_cimp_phase_angles);
  sp_cmap_surfaces = std::move(other.sp_cmap_surfaces);
  sp_cmap_patches = std::move(other.sp_cmap_patches);
  bond_stiffnesses = std::move(other.bond_stiffnesses);
  bond_equilibria = std::move(other.bond_equilibria);
  angl_stiffnesses = std::move(other.angl_stiffnesses);
  angl_equilibria = std::move(other.angl_equilibria);
  dihe_amplitudes = std::move(other.dihe_amplitudes);
  dihe_periodicities = std::move(other.dihe_periodicities);
  dihe_phase_angles = std::move(other.dihe_phase_angles);
  attn14_elec_factors = std::move(other.attn14_elec_factors);
  attn14_vdw_factors = std::move(other.attn14_vdw_factors);
  sp_bond_stiffnesses = std::move(other.sp_bond_stiffnesses);
  sp_bond_equilibria = std::move(other.sp_bond_equilibria);
  sp_angl_stiffnesses = std::move(other.sp_angl_stiffnesses);
  sp_angl_equilibria = std::move(other.sp_angl_equilibria);
  sp_dihe_amplitudes = std::move(other.sp_dihe_amplitudes);
  sp_dihe_periodicities = std::move(other.sp_dihe_periodicities);
  sp_dihe_phase_angles = std::move(other.sp_dihe_phase_angles);
  sp_attn14_elec_factors = std::move(other.sp_attn14_elec_factors);
  sp_attn14_vdw_factors = std::move(other.sp_attn14_vdw_factors);
  valparam_double_data = std::move(other.valparam_double_data);
  valparam_float_data = std::move(other.valparam_float_data);
  valparam_int_data = std::move(other.valparam_int_data);
  valparam_int2_data = std::move(other.valparam_int2_data);

  // Atoms participating in each valence term, plus the parameter indices of each term
  ubrd_i_atoms = std::move(other.ubrd_i_atoms);
  ubrd_k_atoms = std::move(other.ubrd_k_atoms);
  ubrd_param_idx = std::move(other.ubrd_param_idx);
  cimp_i_atoms = std::move(other.cimp_i_atoms);
  cimp_j_atoms = std::move(other.cimp_j_atoms);
  cimp_k_atoms = std::move(other.cimp_k_atoms);
  cimp_l_atoms = std::move(other.cimp_l_atoms);
  cimp_param_idx = std::move(other.cimp_param_idx);
  cmap_i_atoms = std::move(other.cmap_i_atoms);
  cmap_j_atoms = std::move(other.cmap_j_atoms);
  cmap_k_atoms = std::move(other.cmap_k_atoms);
  cmap_l_atoms = std::move(other.cmap_l_atoms);
  cmap_m_atoms = std::move(other.cmap_m_atoms);
  cmap_param_idx = std::move(other.cmap_param_idx);
  bond_i_atoms = std::move(other.bond_i_atoms);
  bond_j_atoms = std::move(other.bond_j_atoms);
  bond_param_idx = std::move(other.bond_param_idx);
  angl_i_atoms = std::move(other.angl_i_atoms);
  angl_j_atoms = std::move(other.angl_j_atoms);
  angl_k_atoms = std::move(other.angl_k_atoms);
  angl_param_idx = std::move(other.angl_param_idx);
  dihe_i_atoms = std::move(other.dihe_i_atoms);
  dihe_j_atoms = std::move(other.dihe_j_atoms);
  dihe_k_atoms = std::move(other.dihe_k_atoms);
  dihe_l_atoms = std::move(other.dihe_l_atoms);
  dihe_param_idx = std::move(other.dihe_param_idx);
  valence_int_data = std::move(other.valence_int_data);

  // Non-bonded parameters for all systems
  charge_indices = std::move(other.charge_indices);
  lennard_jones_indices = std::move(other.lennard_jones_indices);
  charge_parameters = std::move(other.charge_parameters);
  lennard_jones_ab_coeff = std::move(other.lennard_jones_ab_coeff);
  lennard_jones_c_coeff = std::move(other.lennard_jones_c_coeff);
  lennard_jones_14_a_coeff = std::move(other.lennard_jones_14_a_coeff);
  lennard_jones_14_b_coeff = std::move(other.lennard_jones_14_b_coeff);
  lennard_jones_14_c_coeff = std::move(other.lennard_jones_14_c_coeff);
  lennard_jones_sigma = std::move(other.lennard_jones_sigma);
  lennard_jones_14_sigma = std::move(other.lennard_jones_14_sigma);
  sp_charge_parameters = std::move(other.sp_charge_parameters);
  sp_lennard_jones_ab_coeff = std::move(other.sp_lennard_jones_ab_coeff);
  sp_lennard_jones_c_coeff = std::move(other.sp_lennard_jones_c_coeff);
  sp_lennard_jones_14_a_coeff = std::move(other.sp_lennard_jones_14_a_coeff);
  sp_lennard_jones_14_b_coeff = std::move(other.sp_lennard_jones_14_b_coeff);
  sp_lennard_jones_14_c_coeff = std::move(other.sp_lennard_jones_14_c_coeff);
  sp_lennard_jones_sigma = std::move(other.sp_lennard_jones_sigma);
  sp_lennard_jones_14_sigma = std::move(other.sp_lennard_jones_14_sigma);

  // Implicit solvent model parameters
  neck_table_size = other.neck_table_size;
  neck_gb_indices = std::move(other.neck_gb_indices);
  atomic_pb_radii = std::move(other.atomic_pb_radii);
  gb_screening_factors = std::move(other.gb_screening_factors);
  gb_alpha_parameters = std::move(other.gb_alpha_parameters);
  gb_beta_parameters = std::move(other.gb_beta_parameters);
  gb_gamma_parameters = std::move(other.gb_gamma_parameters);
  neck_limit_tables = std::move(other.neck_limit_tables);
  sp_atomic_pb_radii = std::move(other.sp_atomic_pb_radii);
  sp_gb_screening_factors = std::move(other.sp_gb_screening_factors);
  sp_gb_alpha_parameters = std::move(other.sp_gb_alpha_parameters);
  sp_gb_beta_parameters = std::move(other.sp_gb_beta_parameters);
  sp_gb_gamma_parameters = std::move(other.sp_gb_gamma_parameters);
  sp_neck_limit_tables = std::move(other.sp_neck_limit_tables);
  
  // Restraint parameters for all systems
  rposn_step_bounds = std::move(other.rposn_step_bounds);
  rbond_step_bounds = std::move(other.rbond_step_bounds);
  rangl_step_bounds = std::move(other.rangl_step_bounds);
  rdihe_step_bounds = std::move(other.rdihe_step_bounds);
  rposn_init_k = std::move(other.rposn_init_k);
  rposn_final_k = std::move(other.rposn_final_k);
  rposn_init_r = std::move(other.rposn_init_r);
  rposn_final_r = std::move(other.rposn_final_r);
  rposn_init_xy = std::move(other.rposn_init_xy);
  rposn_init_z = std::move(other.rposn_init_z);
  rposn_final_xy = std::move(other.rposn_final_xy);
  rposn_final_z = std::move(other.rposn_final_z);
  rbond_init_k = std::move(other.rbond_init_k);
  rbond_final_k = std::move(other.rbond_final_k);
  rbond_init_r = std::move(other.rbond_init_r);
  rbond_final_r = std::move(other.rbond_final_r);
  rangl_init_k = std::move(other.rangl_init_k);
  rangl_final_k = std::move(other.rangl_final_k);
  rangl_init_r = std::move(other.rangl_init_r);
  rangl_final_r = std::move(other.rangl_final_r);
  rdihe_init_k = std::move(other.rdihe_init_k);
  rdihe_final_k = std::move(other.rdihe_final_k);
  rdihe_init_r = std::move(other.rdihe_init_r);
  rdihe_final_r = std::move(other.rdihe_final_r);
  sp_rposn_init_k = std::move(other.sp_rposn_init_k);
  sp_rposn_final_k = std::move(other.sp_rposn_final_k);
  sp_rposn_init_r = std::move(other.sp_rposn_init_r);
  sp_rposn_final_r = std::move(other.sp_rposn_final_r);
  sp_rposn_init_xy = std::move(other.sp_rposn_init_xy);
  sp_rposn_init_z = std::move(other.sp_rposn_init_z);
  sp_rposn_final_xy = std::move(other.sp_rposn_final_xy);
  sp_rposn_final_z = std::move(other.sp_rposn_final_z);
  sp_rbond_init_k = std::move(other.sp_rbond_init_k);
  sp_rbond_final_k = std::move(other.sp_rbond_final_k);
  sp_rbond_init_r = std::move(other.sp_rbond_init_r);
  sp_rbond_final_r = std::move(other.sp_rbond_final_r);
  sp_rangl_init_k = std::move(other.sp_rangl_init_k);
  sp_rangl_final_k = std::move(other.sp_rangl_final_k);
  sp_rangl_init_r = std::move(other.sp_rangl_init_r);
  sp_rangl_final_r = std::move(other.sp_rangl_final_r);
  sp_rdihe_init_k = std::move(other.sp_rdihe_init_k);
  sp_rdihe_final_k = std::move(other.sp_rdihe_final_k);
  sp_rdihe_init_r = std::move(other.sp_rdihe_init_r);
  sp_rdihe_final_r = std::move(other.sp_rdihe_final_r);
  nmr_double_data = std::move(other.nmr_double_data);
  nmr_double2_data = std::move(other.nmr_double2_data);
  nmr_double4_data = std::move(other.nmr_double4_data);
  nmr_float_data = std::move(other.nmr_float_data);
  nmr_float2_data = std::move(other.nmr_float2_data);
  nmr_float4_data = std::move(other.nmr_float4_data);

  // Restraint term atom indices and parameter indices
  rposn_atoms = std::move(other.rposn_atoms);
  rposn_kr_param_idx = std::move(other.rposn_kr_param_idx);
  rposn_xyz_param_idx = std::move(other.rposn_xyz_param_idx);
  rbond_i_atoms = std::move(other.rbond_i_atoms);
  rbond_j_atoms = std::move(other.rbond_j_atoms);
  rbond_param_idx = std::move(other.rbond_param_idx);
  rangl_i_atoms = std::move(other.rangl_i_atoms);
  rangl_j_atoms = std::move(other.rangl_j_atoms);
  rangl_k_atoms = std::move(other.rangl_k_atoms);
  rangl_param_idx = std::move(other.rangl_param_idx);
  rdihe_i_atoms = std::move(other.rdihe_i_atoms);
  rdihe_j_atoms = std::move(other.rdihe_j_atoms);
  rdihe_k_atoms = std::move(other.rdihe_k_atoms);
  rdihe_l_atoms = std::move(other.rdihe_l_atoms);
  rdihe_param_idx = std::move(other.rdihe_param_idx);
  rposn_kr_param_map = std::move(other.rposn_kr_param_map);
  rposn_xyz_param_map = std::move(other.rposn_xyz_param_map);
  rposn_param_map_bounds = std::move(other.rposn_param_map_bounds);
  rbond_param_map = std::move(other.rbond_param_map);
  rbond_param_map_bounds = std::move(other.rbond_param_map_bounds);
  rangl_param_map = std::move(other.rangl_param_map);
  rangl_param_map_bounds = std::move(other.rangl_param_map_bounds);
  rdihe_param_map = std::move(other.rdihe_param_map);
  rdihe_param_map_bounds = std::move(other.rdihe_param_map_bounds);
  nmr_int_data = std::move(other.nmr_int_data);
  nmr_int2_data = std::move(other.nmr_int2_data);

  // Virtual site parameter arrays, term atom indices, and parameter indices
  virtual_site_parameters = std::move(other.virtual_site_parameters);
  sp_virtual_site_parameters = std::move(other.sp_virtual_site_parameters);
  virtual_site_atoms = std::move(other.virtual_site_atoms);
  virtual_site_frame1_atoms = std::move(other.virtual_site_frame1_atoms);
  virtual_site_frame2_atoms = std::move(other.virtual_site_frame2_atoms);
  virtual_site_frame3_atoms = std::move(other.virtual_site_frame3_atoms);
  virtual_site_frame4_atoms = std::move(other.virtual_site_frame4_atoms);
  virtual_site_parameter_indices = std::move(other.virtual_site_parameter_indices);
  vsite_int_data = std::move(other.vsite_int_data);

  // SETTLE (analytic rigid water constraints) parameter arrays, plus fused atom and parameter
  // indices.
  settle_group_geometry = std::move(other.settle_group_geometry);
  settle_group_masses = std::move(other.settle_group_masses);
  sp_settle_group_geometry = std::move(other.sp_settle_group_geometry);
  sp_settle_group_masses = std::move(other.sp_settle_group_masses);
  settle_group_indexing = std::move(other.settle_group_indexing);

  // Constraint group parameter indices and bounds (referencing consensus tables from the
  // synthesis) plus length / inverse mass parameters (also fused)
  constraint_group_indices = std::move(other.constraint_group_indices);
  constraint_group_bounds = std::move(other.constraint_group_bounds);
  constraint_group_param_idx = std::move(other.constraint_group_param_idx);
  constraint_param_bounds = std::move(other.constraint_param_bounds);
  constraint_group_params = std::move(other.constraint_group_params);
  sp_constraint_group_params = std::move(other.sp_constraint_group_params);

  // Valence work unit instruction sets and energy accumulation masks
  total_valence_work_units = other.total_valence_work_units;
  valence_work_unit_size = other.valence_work_unit_size;
  valence_thread_block_size = other.valence_thread_block_size;
  vwu_instruction_sets = std::move(other.vwu_instruction_sets);
  vwu_import_lists = std::move(other.vwu_import_lists);
  vwu_manipulation_masks = std::move(other.vwu_manipulation_masks);
  cbnd_instructions = std::move(other.cbnd_instructions);
  angl_instructions = std::move(other.angl_instructions);
  cdhe_instructions = std::move(other.cdhe_instructions);
  cdhe_overtones = std::move(other.cdhe_overtones);
  cmap_instructions = std::move(other.cmap_instructions);
  infr14_instructions = std::move(other.infr14_instructions);
  rposn_instructions = std::move(other.rposn_instructions);
  rbond_instructions = std::move(other.rbond_instructions);
  rangl_instructions = std::move(other.rangl_instructions);
  rdihe_instructions = std::move(other.rdihe_instructions);
  vste_instructions = std::move(other.vste_instructions);
  sett_instructions = std::move(other.sett_instructions);
  cnst_instructions = std::move(other.cnst_instructions);
  accumulate_cbnd_energy = std::move(other.accumulate_cbnd_energy);
  accumulate_angl_energy = std::move(other.accumulate_angl_energy);
  accumulate_cdhe_energy = std::move(other.accumulate_cdhe_energy);
  accumulate_cmap_energy = std::move(other.accumulate_cmap_energy);
  accumulate_infr14_energy = std::move(other.accumulate_infr14_energy);
  accumulate_rposn_energy = std::move(other.accumulate_rposn_energy);
  accumulate_rbond_energy = std::move(other.accumulate_rbond_energy);
  accumulate_rangl_energy = std::move(other.accumulate_rangl_energy);
  accumulate_rdihe_energy = std::move(other.accumulate_rdihe_energy);
  insr_uint_data = std::move(other.insr_uint_data);
  insr_uint2_data = std::move(other.insr_uint2_data);

  // Non-bonded work units and their instructions
  total_nonbonded_work_units = other.total_nonbonded_work_units;
  nonbonded_work_type = std::move(other.nonbonded_work_type);
  nonbonded_abstracts = std::move(other.nonbonded_abstracts);
  nbwu_instructions = std::move(other.nbwu_instructions);

  // Reduction work units and their abstracts
  total_reduction_work_units = other.total_reduction_work_units;
  rdwu_per_system = other.rdwu_per_system;
  reduction_abstracts = std::move(other.reduction_abstracts);
  
  // Pointer to the timer used to track time needed for assembling this object.  The timer can
  // also be used to track access and other usage of the object.
  timer = std::move(other.timer);
  return *this;
}

//-------------------------------------------------------------------------------------------------
void AtomGraphSynthesis::rebasePointers() {

  // Individual topologies and restraint networks (RestraintApparatus objects)
  topology_indices.swapTarget(&int_system_data);
  restraint_indices.swapTarget(&int_system_data);
  atom_counts.swapTarget(&int_system_data);
  residue_counts.swapTarget(&int_system_data);
  molecule_counts.swapTarget(&int_system_data);
  largest_residue_sizes.swapTarget(&int_system_data);
  last_solute_residues.swapTarget(&int_system_data);
  last_solute_atoms.swapTarget(&int_system_data);
  first_solvent_molecules.swapTarget(&int_system_data);

  // Counts of energy terms and particle-related quantities
  ubrd_term_counts.swapTarget(&int_system_data);
  cimp_term_counts.swapTarget(&int_system_data);
  cmap_term_counts.swapTarget(&int_system_data);
  bond_term_counts.swapTarget(&int_system_data);
  angl_term_counts.swapTarget(&int_system_data);
  dihe_term_counts.swapTarget(&int_system_data);
  virtual_site_counts.swapTarget(&int_system_data);
  posn_restraint_counts.swapTarget(&int_system_data);
  bond_restraint_counts.swapTarget(&int_system_data);
  angl_restraint_counts.swapTarget(&int_system_data);
  dihe_restraint_counts.swapTarget(&int_system_data);
  lj_type_counts.swapTarget(&int_system_data);
  total_exclusion_counts.swapTarget(&int_system_data);
  rigid_water_counts.swapTarget(&int_system_data);
  bond_constraint_counts.swapTarget(&int_system_data);
  degrees_of_freedom.swapTarget(&int_system_data);
  cnst_degrees_of_freedom.swapTarget(&int_system_data);
  nonrigid_particle_counts.swapTarget(&int_system_data);

  // Offsets for each system's term lists
  atom_offsets.swapTarget(&int_system_data);
  atom_bit_offsets.swapTarget(&int_system_data);
  residue_offsets.swapTarget(&int_system_data);
  molecule_offsets.swapTarget(&int_system_data);
  ubrd_term_offsets.swapTarget(&int_system_data);
  cimp_term_offsets.swapTarget(&int_system_data);
  cmap_term_offsets.swapTarget(&int_system_data);
  bond_term_offsets.swapTarget(&int_system_data);
  angl_term_offsets.swapTarget(&int_system_data);
  dihe_term_offsets.swapTarget(&int_system_data);
  virtual_site_offsets.swapTarget(&int_system_data);
  posn_restraint_offsets.swapTarget(&int_system_data);
  bond_restraint_offsets.swapTarget(&int_system_data);
  angl_restraint_offsets.swapTarget(&int_system_data);
  dihe_restraint_offsets.swapTarget(&int_system_data);
  sett_group_offsets.swapTarget(&int_system_data);
  cnst_group_offsets.swapTarget(&int_system_data);
  nb_exclusion_offsets.swapTarget(&int_system_data);
  lennard_jones_abc_offsets.swapTarget(&int_system_data);

  // Atom and residue details
  residue_limits.swapTarget(&chem_int2_data);
  atom_struc_numbers.swapTarget(&chem_int_data);
  residue_numbers.swapTarget(&chem_int_data);
  molecule_limits.swapTarget(&chem_int2_data);
  atomic_numbers.swapTarget(&chem_int_data);
  mobile_atoms.swapTarget(&chem_int_data);
  molecule_membership.swapTarget(&chem_int_data);
  molecule_contents.swapTarget(&chem_int_data);
  atomic_charges.swapTarget(&chem_double_data);
  atomic_masses.swapTarget(&chem_double_data);
  inverse_atomic_masses.swapTarget(&chem_double_data);
  sp_atomic_charges.swapTarget(&chem_float_data);
  sp_atomic_masses.swapTarget(&chem_float_data);
  sp_inverse_atomic_masses.swapTarget(&chem_float_data);
  atom_names.swapTarget(&chem_char4_data);
  atom_types.swapTarget(&chem_char4_data);
  residue_names.swapTarget(&chem_char4_data);

  // Force field parameter maps from individual systems to the consensus tables
  ubrd_param_map.swapTarget(&valparam_int_data);
  ubrd_param_map_bounds.swapTarget(&valparam_int2_data);
  cimp_param_map.swapTarget(&valparam_int_data);
  cimp_param_map_bounds.swapTarget(&valparam_int2_data);
  cmap_param_map.swapTarget(&valparam_int_data);
  cmap_param_map_bounds.swapTarget(&valparam_int2_data);
  bond_param_map.swapTarget(&valparam_int_data);
  bond_param_map_bounds.swapTarget(&valparam_int2_data);
  angl_param_map.swapTarget(&valparam_int_data);
  angl_param_map_bounds.swapTarget(&valparam_int2_data);
  dihe_param_map.swapTarget(&valparam_int_data);
  dihe_param_map_bounds.swapTarget(&valparam_int2_data);
  attn14_param_map.swapTarget(&valparam_int_data);
  attn14_param_map_bounds.swapTarget(&valparam_int2_data);
  vste_param_map.swapTarget(&valparam_int_data);
  vste_param_map_bounds.swapTarget(&valparam_int2_data);
  sett_param_map.swapTarget(&valparam_int_data);
  sett_param_map_bounds.swapTarget(&valparam_int2_data);
  cnst_param_map.swapTarget(&valparam_int_data);
  cnst_param_map_bounds.swapTarget(&valparam_int2_data);

  // Force field consensus table condensed parameters
  ubrd_stiffnesses.swapTarget(&valparam_double_data);
  ubrd_equilibria.swapTarget(&valparam_double_data);
  cimp_stiffnesses.swapTarget(&valparam_double_data);
  cimp_phase_angles.swapTarget(&valparam_double_data);
  cmap_surface_dimensions.swapTarget(&valparam_int_data);
  cmap_surface_bounds.swapTarget(&valparam_int_data);
  cmap_patch_bounds.swapTarget(&valparam_int_data);
  cmap_surfaces.swapTarget(&valparam_double_data);
  cmap_patches.swapTarget(&valparam_double_data);
  sp_ubrd_stiffnesses.swapTarget(&valparam_float_data);
  sp_ubrd_equilibria.swapTarget(&valparam_float_data);
  sp_cimp_stiffnesses.swapTarget(&valparam_float_data);
  sp_cimp_phase_angles.swapTarget(&valparam_float_data);
  sp_cmap_surfaces.swapTarget(&valparam_float_data);
  sp_cmap_patches.swapTarget(&valparam_float_data);
  bond_stiffnesses.swapTarget(&valparam_double_data);
  bond_equilibria.swapTarget(&valparam_double_data);
  angl_stiffnesses.swapTarget(&valparam_double_data);
  angl_equilibria.swapTarget(&valparam_double_data);
  dihe_amplitudes.swapTarget(&valparam_double_data);
  dihe_periodicities.swapTarget(&valparam_double_data);
  dihe_phase_angles.swapTarget(&valparam_double_data);
  attn14_elec_factors.swapTarget(&valparam_double_data);
  attn14_vdw_factors.swapTarget(&valparam_double_data);
  sp_bond_stiffnesses.swapTarget(&valparam_float_data);
  sp_bond_equilibria.swapTarget(&valparam_float_data);
  sp_angl_stiffnesses.swapTarget(&valparam_float_data);
  sp_angl_equilibria.swapTarget(&valparam_float_data);
  sp_dihe_amplitudes.swapTarget(&valparam_float_data);
  sp_dihe_periodicities.swapTarget(&valparam_float_data);
  sp_dihe_phase_angles.swapTarget(&valparam_float_data);
  sp_attn14_elec_factors.swapTarget(&valparam_float_data);
  sp_attn14_vdw_factors.swapTarget(&valparam_float_data);
  
  // Atoms participating in each valence term, plus the parameter indices of each term
  ubrd_i_atoms.swapTarget(&valence_int_data);
  ubrd_k_atoms.swapTarget(&valence_int_data);
  ubrd_param_idx.swapTarget(&valence_int_data);
  cimp_i_atoms.swapTarget(&valence_int_data);
  cimp_j_atoms.swapTarget(&valence_int_data);
  cimp_k_atoms.swapTarget(&valence_int_data);
  cimp_l_atoms.swapTarget(&valence_int_data);
  cimp_param_idx.swapTarget(&valence_int_data);
  cmap_i_atoms.swapTarget(&valence_int_data);
  cmap_j_atoms.swapTarget(&valence_int_data);
  cmap_k_atoms.swapTarget(&valence_int_data);
  cmap_l_atoms.swapTarget(&valence_int_data);
  cmap_m_atoms.swapTarget(&valence_int_data);
  cmap_param_idx.swapTarget(&valence_int_data);
  bond_i_atoms.swapTarget(&valence_int_data);
  bond_j_atoms.swapTarget(&valence_int_data);
  bond_param_idx.swapTarget(&valence_int_data);
  angl_i_atoms.swapTarget(&valence_int_data);
  angl_j_atoms.swapTarget(&valence_int_data);
  angl_k_atoms.swapTarget(&valence_int_data);
  angl_param_idx.swapTarget(&valence_int_data);
  dihe_i_atoms.swapTarget(&valence_int_data);
  dihe_j_atoms.swapTarget(&valence_int_data);
  dihe_k_atoms.swapTarget(&valence_int_data);
  dihe_l_atoms.swapTarget(&valence_int_data);
  dihe_param_idx.swapTarget(&valence_int_data);

  // Non-bonded parameter indices for all systems
  charge_indices.swapTarget(&chem_int_data);
  lennard_jones_indices.swapTarget(&chem_int_data);

  // Implicit solvent model parameters
  neck_gb_indices.swapTarget(&chem_int_data);
  atomic_pb_radii.swapTarget(&chem_double_data);
  gb_screening_factors.swapTarget(&chem_double_data);
  gb_alpha_parameters.swapTarget(&chem_double_data);
  gb_beta_parameters.swapTarget(&chem_double_data);
  gb_gamma_parameters.swapTarget(&chem_double_data);
  sp_atomic_pb_radii.swapTarget(&chem_float_data);
  sp_gb_screening_factors.swapTarget(&chem_float_data);
  sp_gb_alpha_parameters.swapTarget(&chem_float_data);
  sp_gb_beta_parameters.swapTarget(&chem_float_data);
  sp_gb_gamma_parameters.swapTarget(&chem_float_data);

  // Restraint parameters for all systems
  rposn_step_bounds.swapTarget(&nmr_int2_data);
  rbond_step_bounds.swapTarget(&nmr_int2_data);
  rangl_step_bounds.swapTarget(&nmr_int2_data);
  rdihe_step_bounds.swapTarget(&nmr_int2_data);
  rposn_init_k.swapTarget(&nmr_double2_data);
  rposn_final_k.swapTarget(&nmr_double2_data);
  rposn_init_r.swapTarget(&nmr_double4_data);
  rposn_final_r.swapTarget(&nmr_double4_data);
  rposn_init_xy.swapTarget(&nmr_double2_data);
  rposn_init_z.swapTarget(&nmr_double_data);
  rposn_final_xy.swapTarget(&nmr_double2_data);
  rposn_final_z.swapTarget(&nmr_double_data);
  rbond_init_k.swapTarget(&nmr_double2_data);
  rbond_final_k.swapTarget(&nmr_double2_data);
  rbond_init_r.swapTarget(&nmr_double4_data);
  rbond_final_r.swapTarget(&nmr_double4_data);
  rangl_init_k.swapTarget(&nmr_double2_data);
  rangl_final_k.swapTarget(&nmr_double2_data);
  rangl_init_r.swapTarget(&nmr_double4_data);
  rangl_final_r.swapTarget(&nmr_double4_data);
  rdihe_init_k.swapTarget(&nmr_double2_data);
  rdihe_final_k.swapTarget(&nmr_double2_data);
  rdihe_init_r.swapTarget(&nmr_double4_data);
  rdihe_final_r.swapTarget(&nmr_double4_data);
  sp_rposn_init_k.swapTarget(&nmr_float2_data);
  sp_rposn_final_k.swapTarget(&nmr_float2_data);
  sp_rposn_init_r.swapTarget(&nmr_float4_data);
  sp_rposn_final_r.swapTarget(&nmr_float4_data);
  sp_rposn_init_xy.swapTarget(&nmr_float2_data);
  sp_rposn_init_z.swapTarget(&nmr_float_data);
  sp_rposn_final_xy.swapTarget(&nmr_float2_data);
  sp_rposn_final_z.swapTarget(&nmr_float_data);
  sp_rbond_init_k.swapTarget(&nmr_float2_data);
  sp_rbond_final_k.swapTarget(&nmr_float2_data);
  sp_rbond_init_r.swapTarget(&nmr_float4_data);
  sp_rbond_final_r.swapTarget(&nmr_float4_data);
  sp_rangl_init_k.swapTarget(&nmr_float2_data);
  sp_rangl_final_k.swapTarget(&nmr_float2_data);
  sp_rangl_init_r.swapTarget(&nmr_float4_data);
  sp_rangl_final_r.swapTarget(&nmr_float4_data);
  sp_rdihe_init_k.swapTarget(&nmr_float2_data);
  sp_rdihe_final_k.swapTarget(&nmr_float2_data);
  sp_rdihe_init_r.swapTarget(&nmr_float4_data);
  sp_rdihe_final_r.swapTarget(&nmr_float4_data);

  // Restraint term atom indices and parameter indices
  rposn_atoms.swapTarget(&nmr_int_data);
  rposn_kr_param_idx.swapTarget(&nmr_int_data);
  rposn_xyz_param_idx.swapTarget(&nmr_int_data);
  rbond_i_atoms.swapTarget(&nmr_int_data);
  rbond_j_atoms.swapTarget(&nmr_int_data);
  rbond_param_idx.swapTarget(&nmr_int_data);
  rangl_i_atoms.swapTarget(&nmr_int_data);
  rangl_j_atoms.swapTarget(&nmr_int_data);
  rangl_k_atoms.swapTarget(&nmr_int_data);
  rangl_param_idx.swapTarget(&nmr_int_data);
  rdihe_i_atoms.swapTarget(&nmr_int_data);
  rdihe_j_atoms.swapTarget(&nmr_int_data);
  rdihe_k_atoms.swapTarget(&nmr_int_data);
  rdihe_l_atoms.swapTarget(&nmr_int_data);
  rdihe_param_idx.swapTarget(&nmr_int_data);
  rposn_kr_param_map.swapTarget(&nmr_int_data);
  rposn_xyz_param_map.swapTarget(&nmr_int_data);
  rposn_param_map_bounds.swapTarget(&nmr_int2_data);
  rbond_param_map.swapTarget(&nmr_int_data);
  rbond_param_map_bounds.swapTarget(&nmr_int2_data);
  rangl_param_map.swapTarget(&nmr_int_data);
  rangl_param_map_bounds.swapTarget(&nmr_int2_data);
  rdihe_param_map.swapTarget(&nmr_int_data);
  rdihe_param_map_bounds.swapTarget(&nmr_int2_data);

  // Virtual site parameter arrays, term atom indices, and parameter indices
  virtual_site_atoms.swapTarget(&vsite_int_data);
  virtual_site_frame1_atoms.swapTarget(&vsite_int_data);
  virtual_site_frame2_atoms.swapTarget(&vsite_int_data);
  virtual_site_frame3_atoms.swapTarget(&vsite_int_data);
  virtual_site_frame4_atoms.swapTarget(&vsite_int_data);
  virtual_site_parameter_indices.swapTarget(&vsite_int_data);

  // Valence work unit instruction sets and energy accumulation masks
  vwu_manipulation_masks.swapTarget(&insr_uint2_data);
  cbnd_instructions.swapTarget(&insr_uint2_data);
  angl_instructions.swapTarget(&insr_uint2_data);
  cdhe_instructions.swapTarget(&insr_uint2_data);
  cdhe_overtones.swapTarget(&insr_uint_data);
  cmap_instructions.swapTarget(&insr_uint2_data);
  infr14_instructions.swapTarget(&insr_uint_data);
  rposn_instructions.swapTarget(&insr_uint2_data);
  rbond_instructions.swapTarget(&insr_uint2_data);
  rangl_instructions.swapTarget(&insr_uint2_data);
  rdihe_instructions.swapTarget(&insr_uint2_data);
  vste_instructions.swapTarget(&insr_uint2_data);
  sett_instructions.swapTarget(&insr_uint2_data);
  cnst_instructions.swapTarget(&insr_uint2_data);
  accumulate_cbnd_energy.swapTarget(&insr_uint_data);
  accumulate_angl_energy.swapTarget(&insr_uint_data);
  accumulate_cdhe_energy.swapTarget(&insr_uint_data);
  accumulate_cmap_energy.swapTarget(&insr_uint_data);
  accumulate_infr14_energy.swapTarget(&insr_uint_data);
  accumulate_rposn_energy.swapTarget(&insr_uint_data);
  accumulate_rbond_energy.swapTarget(&insr_uint_data);
  accumulate_rangl_energy.swapTarget(&insr_uint_data);
  accumulate_rdihe_energy.swapTarget(&insr_uint_data);
}

} // namespace synthesis
} // namespace stormm
