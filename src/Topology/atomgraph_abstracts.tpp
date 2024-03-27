// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace topology {

//-------------------------------------------------------------------------------------------------
template <typename T>
ValenceKit<T>::ValenceKit(const int natom_in, const int nbond_in, const int nangl_in,
                          const int ndihe_in, const int nbond_param_in, const int nangl_param_in,
                          const int ndihe_param_in, const int ninfr14_in,
                          const int nattn14_param_in, const int nubrd_in, const int ncimp_in,
                          const int ncmap_in, const int nubrd_param_in, const int ncimp_param_in,
                          const int ncmap_surf_in, const T* bond_keq_in,
                          const T* bond_leq_in, const T* angl_keq_in,
                          const T* angl_theta_in, const T* dihe_amp_in, const T* dihe_freq_in,
                          const T* dihe_phi_in, const T* attn14_elec_in, const T* attn14_vdw_in,
                          const int* bond_i_atoms_in, const int* bond_j_atoms_in,
                          const int* bond_param_idx_in, const char4* bond_modifiers_in,
                          const int* angl_i_atoms_in, const int* angl_j_atoms_in,
                          const int* angl_k_atoms_in, const int* angl_param_idx_in,
                          const char4* angl_modifiers_in, const int* dihe_i_atoms_in,
                          const int* dihe_j_atoms_in, const int* dihe_k_atoms_in,
                          const int* dihe_l_atoms_in, const int* dihe_param_idx_in,
                          const int* dihe14_param_idx_in, const char4* dihe_modifiers_in,
                          const int* infr14_i_atoms_in, const int* infr14_l_atoms_in,
                          const int* infr14_param_idx_in, const int* ubrd_i_atoms_in,
                          const int* ubrd_k_atoms_in, const int* ubrd_param_idx_in,
                          const int* cimp_i_atoms_in, const int* cimp_j_atoms_in,
                          const int* cimp_k_atoms_in, const int* cimp_l_atoms_in,
                          const int* cimp_param_idx_in, const int* cmap_i_atoms_in,
                          const int* cmap_j_atoms_in, const int* cmap_k_atoms_in,
                          const int* cmap_l_atoms_in, const int* cmap_m_atoms_in,
                          const int* cmap_dim_in, const int* cmap_surf_bounds_in,
                          const int* cmap_patch_bounds_in, const int* cmap_surf_idx_in,
                          const T* ubrd_keq_in, const T* ubrd_leq_in,
                          const T* cimp_keq_in, const T* cimp_phi_in, const T* cmap_surf_in,
                          const T* cmap_dphi_in, const T* cmap_dpsi_in, const T* cmap_dphi_dpsi_in,
                          const T* cmap_patches_in, const int* bond_asgn_atoms_in,
                          const int* bond_asgn_index_in, const int* bond_asgn_terms_in,
                          const int* bond_asgn_bounds_in, const int* angl_asgn_atoms_in,
                          const int* angl_asgn_index_in, const int* angl_asgn_terms_in,
                          const int* angl_asgn_bounds_in, const int* dihe_asgn_atoms_in,
                          const int* dihe_asgn_index_in, const int* dihe_asgn_terms_in,
                          const int* dihe_asgn_bounds_in, const int* ubrd_asgn_atoms_in,
                          const int* ubrd_asgn_index_in, const int* ubrd_asgn_terms_in,
                          const int* ubrd_asgn_bounds_in, const int* cimp_asgn_atoms_in,
                          const int* cimp_asgn_index_in, const int* cimp_asgn_terms_in,
                          const int* cimp_asgn_bounds_in, const int* cmap_asgn_atoms_in,
                          const int* cmap_asgn_index_in, const int* cmap_asgn_terms_in,
                          const int* cmap_asgn_bounds_in) :
    natom{natom_in}, nbond{nbond_in}, nangl{nangl_in}, ndihe{ndihe_in},
    nbond_param{nbond_param_in}, nangl_param{nangl_param_in}, ndihe_param{ndihe_param_in},
    ninfr14{ninfr14_in}, nattn14_param{nattn14_param_in}, nubrd{nubrd_in}, ncimp{ncimp_in},
    ncmap{ncmap_in}, nubrd_param{nubrd_param_in}, ncimp_param{ncimp_param_in},
    ncmap_surf{ncmap_surf_in}, bond_keq{bond_keq_in}, bond_leq{bond_leq_in}, angl_keq{angl_keq_in},
    angl_theta{angl_theta_in}, dihe_amp{dihe_amp_in}, dihe_freq{dihe_freq_in},
    dihe_phi{dihe_phi_in}, attn14_elec{attn14_elec_in}, attn14_vdw{attn14_vdw_in},
    bond_i_atoms{bond_i_atoms_in}, bond_j_atoms{bond_j_atoms_in},
    bond_param_idx{bond_param_idx_in}, bond_modifiers{bond_modifiers_in},
    angl_i_atoms{angl_i_atoms_in}, angl_j_atoms{angl_j_atoms_in}, angl_k_atoms{angl_k_atoms_in},
    angl_param_idx{angl_param_idx_in}, angl_modifiers{angl_modifiers_in},
    dihe_i_atoms{dihe_i_atoms_in}, dihe_j_atoms{dihe_j_atoms_in}, dihe_k_atoms{dihe_k_atoms_in},
    dihe_l_atoms{dihe_l_atoms_in}, dihe_param_idx{dihe_param_idx_in},
    dihe14_param_idx{dihe14_param_idx_in}, dihe_modifiers{dihe_modifiers_in},
    infr14_i_atoms{infr14_i_atoms_in}, infr14_l_atoms{infr14_l_atoms_in},
    infr14_param_idx{infr14_param_idx_in}, ubrd_i_atoms{ubrd_i_atoms_in},
    ubrd_k_atoms{ubrd_k_atoms_in}, ubrd_param_idx{ubrd_param_idx_in},
    cimp_i_atoms{cimp_i_atoms_in}, cimp_j_atoms{cimp_j_atoms_in}, cimp_k_atoms{cimp_k_atoms_in},
    cimp_l_atoms{cimp_l_atoms_in}, cimp_param_idx{cimp_param_idx_in},
    cmap_i_atoms{cmap_i_atoms_in}, cmap_j_atoms{cmap_j_atoms_in}, cmap_k_atoms{cmap_k_atoms_in},
    cmap_l_atoms{cmap_l_atoms_in}, cmap_m_atoms{cmap_m_atoms_in}, cmap_dim{cmap_dim_in},
    cmap_surf_bounds{cmap_surf_bounds_in}, cmap_patch_bounds{cmap_patch_bounds_in},
    cmap_surf_idx{cmap_surf_idx_in}, ubrd_keq{ubrd_keq_in}, ubrd_leq{ubrd_leq_in},
    cimp_keq{cimp_keq_in}, cimp_phi{cimp_phi_in}, cmap_surf{cmap_surf_in}, cmap_dphi{cmap_dphi_in},
    cmap_dpsi{cmap_dpsi_in}, cmap_dphi_dpsi{cmap_dphi_dpsi_in}, cmap_patches{cmap_patches_in},
    bond_asgn_atoms{bond_asgn_atoms_in}, bond_asgn_index{bond_asgn_index_in},
    bond_asgn_terms{bond_asgn_terms_in}, bond_asgn_bounds{bond_asgn_bounds_in},
    angl_asgn_atoms{angl_asgn_atoms_in}, angl_asgn_index{angl_asgn_index_in},
    angl_asgn_terms{angl_asgn_terms_in}, angl_asgn_bounds{angl_asgn_bounds_in},
    dihe_asgn_atoms{dihe_asgn_atoms_in}, dihe_asgn_index{dihe_asgn_index_in},
    dihe_asgn_terms{dihe_asgn_terms_in}, dihe_asgn_bounds{dihe_asgn_bounds_in},
    ubrd_asgn_atoms{ubrd_asgn_atoms_in}, ubrd_asgn_index{ubrd_asgn_index_in},
    ubrd_asgn_terms{ubrd_asgn_terms_in}, ubrd_asgn_bounds{ubrd_asgn_bounds_in},
    cimp_asgn_atoms{cimp_asgn_atoms_in}, cimp_asgn_index{cimp_asgn_index_in},
    cimp_asgn_terms{cimp_asgn_terms_in}, cimp_asgn_bounds{cimp_asgn_bounds_in},
    cmap_asgn_atoms{cmap_asgn_atoms_in}, cmap_asgn_index{cmap_asgn_index_in},
    cmap_asgn_terms{cmap_asgn_terms_in}, cmap_asgn_bounds{cmap_asgn_bounds_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
NonbondedKit<T>::NonbondedKit(const int natom_in, const int n_q_types_in, const int n_lj_types_in,
                              const T coulomb_constant_in, const T* charge_in, const int* q_idx_in,
                              const int* lj_idx_in, const T* q_parameter_in, const T* lja_coeff_in,
                              const T* ljb_coeff_in, const T* ljc_coeff_in,
                              const T* lja_14_coeff_in, const T* ljb_14_coeff_in,
                              const T* ljc_14_coeff_in, const T* lj_sigma_in,
                              const T* lj_14_sigma_in, const int* nb11x_in,
                              const int* nb11_bounds_in, const int* nb12x_in,
                              const int* nb12_bounds_in, const int* nb13x_in,
                              const int* nb13_bounds_in, const int* nb14x_in,
                              const int* nb14_bounds_in, const T* lj_type_corr_in) :
    natom{natom_in}, n_q_types{n_q_types_in}, n_lj_types{n_lj_types_in},
    coulomb_constant{coulomb_constant_in}, charge{charge_in}, q_idx{q_idx_in}, lj_idx{lj_idx_in},
    q_parameter{q_parameter_in}, lja_coeff{lja_coeff_in}, ljb_coeff{ljb_coeff_in},
    ljc_coeff{ljc_coeff_in}, lja_14_coeff{lja_14_coeff_in}, ljb_14_coeff{ljb_14_coeff_in},
    ljc_14_coeff{ljc_14_coeff_in}, lj_sigma{lj_sigma_in}, lj_14_sigma{lj_14_sigma_in},
    nb11x{nb11x_in}, nb11_bounds{nb11_bounds_in}, nb12x{nb12x_in}, nb12_bounds{nb12_bounds_in},
    nb13x{nb13x_in}, nb13_bounds{nb13_bounds_in}, nb14x{nb14x_in}, nb14_bounds{nb14_bounds_in},
    lj_type_corr{lj_type_corr_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
ImplicitSolventKit<T>::ImplicitSolventKit(int natom_in, ImplicitSolventModel igb_in,
                                          const T dielectric_in, const T saltcon_in,
                                          const int* neck_gb_idx_in, const T* pb_radii_in,
                                          const T* gb_screen_in, const T* gb_alpha_in,
                                          const T* gb_beta_in, const T* gb_gamma_in) :
    natom{natom_in}, igb{igb_in}, dielectric{dielectric_in}, saltcon{saltcon_in},
    neck_gb_idx{neck_gb_idx_in}, pb_radii{pb_radii_in}, gb_screen{gb_screen_in},
    gb_alpha{gb_alpha_in}, gb_beta{gb_beta_in}, gb_gamma{gb_gamma_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
VirtualSiteKit<T>::VirtualSiteKit(const int nsite_in, const int nframe_set_in,
                                  const int* vs_atoms_in, const int* frame1_idx_in,
                                  const int* frame2_idx_in, const int* frame3_idx_in,
                                  const int* frame4_idx_in, const int* vs_param_idx_in,
                                  const int* vs_types_in, const T* dim1_in, const T* dim2_in,
                                  const T* dim3_in) :
    nsite{nsite_in}, nframe_set{nframe_set_in}, vs_atoms{vs_atoms_in}, frame1_idx{frame1_idx_in},
    frame2_idx{frame2_idx_in}, frame3_idx{frame3_idx_in}, frame4_idx{frame4_idx_in},
    vs_param_idx{vs_param_idx_in}, vs_types{vs_types_in}, dim1{dim1_in}, dim2{dim2_in},
    dim3{dim3_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
ConstraintKit<T>::ConstraintKit(const int nsettle_in, const int nsett_param_in,
                                const int ngroup_in, const int ncnst_param_in,
                                const int* settle_ox_atoms_in, const int* settle_h1_atoms_in,
                                const int* settle_h2_atoms_in, const int* settle_param_idx_in,
                                const int* group_list_in, const int* group_bounds_in,
                                const int* group_param_idx_in, const int* group_param_bounds_in,
                                const T* settle_mormt_in, const T* settle_mhrmt_in,
                                const T* settle_ra_in, const T* settle_rb_in,
                                const T* settle_rc_in, const T* settle_invra_in,
                                const T* group_sq_lengths_in, const T* group_inv_masses_in) :
    nsettle{nsettle_in}, nsett_param{nsett_param_in}, ngroup{ngroup_in},
    ncnst_param{ncnst_param_in}, settle_ox_atoms{settle_ox_atoms_in},
    settle_h1_atoms{settle_h1_atoms_in}, settle_h2_atoms{settle_h2_atoms_in},
    settle_param_idx{settle_param_idx_in}, group_list{group_list_in},
    group_bounds{group_bounds_in}, group_param_idx{group_param_idx_in},
    group_param_bounds{group_param_bounds_in}, settle_mormt{settle_mormt_in},
    settle_mhrmt{settle_mhrmt_in}, settle_ra{settle_ra_in}, settle_rb{settle_rb_in},
    settle_rc{settle_rc_in}, settle_invra{settle_invra_in}, group_sq_lengths{group_sq_lengths_in},
    group_inv_masses{group_inv_masses_in}
{}

} // namespace topology
} // namespace stormm
