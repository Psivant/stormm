// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
template <typename T>
SyValenceKit<T>::SyValenceKit(const int nvwu_in, const T coulomb_in, const T* bond_keq_in,
                              const T* bond_leq_in, const T* angl_keq_in, const T* angl_theta_in,
                              const T* dihe_amp_in, const T* dihe_freq_in, const T* dihe_phi_in,
                              const T* attn14_elec_in, const T* attn14_vdw_in, const T* charges_in,
                              const T* lja_14_coeff_in, const T* ljb_14_coeff_in,
                              const T* ljc_14_coeff_in, const T* lj_14_sigma_in,
                              const int* lj_idx_in, const int* n_lj_types_in,
                              const int* ljabc_offsets_in, const T* ubrd_keq_in,
                              const T* ubrd_leq_in, const T* cimp_keq_in, const T* cimp_phi_in,
                              const int* cmap_dim_in, const T* cmap_patches_in,
                              const int* cmap_patch_bounds_in, const int2* vwu_abstracts_in,
                              const int* vwu_imports_in, const uint2* cbnd_insr_in,
                              const uint2* angl_insr_in, const uint2* cdhe_insr_in,
                              const uint* cdhe_ovrt_insr_in, const uint2* cmap_insr_in,
                              const uint* infr14_insr_in, const uint* cbnd_acc_in,
                              const uint* angl_acc_in, const uint* cdhe_acc_in,
                              const uint* cmap_acc_in, const uint* infr14_acc_in) :
    nvwu{nvwu_in}, coulomb{coulomb_in}, bond_keq{bond_keq_in}, bond_leq{bond_leq_in},
    angl_keq{angl_keq_in}, angl_theta{angl_theta_in}, dihe_amp{dihe_amp_in},
    dihe_freq{dihe_freq_in}, dihe_phi{dihe_phi_in}, attn14_elec{attn14_elec_in},
    attn14_vdw{attn14_vdw_in}, charges{charges_in}, lja_14_coeff{lja_14_coeff_in},
    ljb_14_coeff{ljb_14_coeff_in}, ljc_14_coeff{ljc_14_coeff_in}, lj_14_sigma{lj_14_sigma_in},
    lj_idx{lj_idx_in}, n_lj_types{n_lj_types_in}, ljabc_offsets{ljabc_offsets_in},
    ubrd_keq{ubrd_keq_in}, ubrd_leq{ubrd_leq_in}, cimp_keq{cimp_keq_in}, cimp_phi{cimp_phi_in},
    cmap_dim{cmap_dim_in}, cmap_patches{cmap_patches_in}, cmap_patch_bounds{cmap_patch_bounds_in},
    vwu_abstracts{vwu_abstracts_in}, vwu_imports{vwu_imports_in}, cbnd_insr{cbnd_insr_in},
    angl_insr{angl_insr_in}, cdhe_insr{cdhe_insr_in}, cdhe_ovrt_insr{cdhe_ovrt_insr_in},
    cmap_insr{cmap_insr_in}, infr14_insr{infr14_insr_in}, cbnd_acc{cbnd_acc_in},
    angl_acc{angl_acc_in}, cdhe_acc{cdhe_acc_in}, cmap_acc{cmap_acc_in}, infr14_acc{infr14_acc_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2, typename T4>
SyRestraintKit<T, T2, T4>::SyRestraintKit(const int2* rposn_step_bounds_in,
                                          const int2* rbond_step_bounds_in,
                                          const int2* rangl_step_bounds_in,
                                          const int2* rdihe_step_bounds_in,
                                          const T2* rposn_init_k_in, const T2* rposn_finl_k_in,
                                          const T4* rposn_init_r_in, const T4* rposn_finl_r_in,
                                          const T2* rposn_init_xy_in, const T* rposn_init_z_in,
                                          const T2* rposn_finl_xy_in, const T* rposn_finl_z_in,
                                          const T2* rbond_init_k_in, const T2* rbond_finl_k_in,
                                          const T4* rbond_init_r_in, const T4* rbond_finl_r_in,
                                          const T2* rangl_init_k_in, const T2* rangl_finl_k_in,
                                          const T4* rangl_init_r_in, const T4* rangl_finl_r_in,
                                          const T2* rdihe_init_k_in, const T2* rdihe_finl_k_in,
                                          const T4* rdihe_init_r_in, const T4* rdihe_finl_r_in,
                                          const uint2* rposn_insr_in, const uint2* rbond_insr_in,
                                          const uint2* rangl_insr_in, const uint2* rdihe_insr_in,
                                          const uint* rposn_acc_in, const uint* rbond_acc_in,
                                          const uint* rangl_acc_in, const uint* rdihe_acc_in) :
    rposn_step_bounds{rposn_step_bounds_in}, rbond_step_bounds{rbond_step_bounds_in},
    rangl_step_bounds{rangl_step_bounds_in}, rdihe_step_bounds{rdihe_step_bounds_in},
    rposn_init_k{rposn_init_k_in}, rposn_finl_k{rposn_finl_k_in}, rposn_init_r{rposn_init_r_in},
    rposn_finl_r{rposn_finl_r_in}, rposn_init_xy{rposn_init_xy_in}, rposn_init_z{rposn_init_z_in},
    rposn_finl_xy{rposn_finl_xy_in}, rposn_finl_z{rposn_finl_z_in}, rbond_init_k{rbond_init_k_in},
    rbond_finl_k{rbond_finl_k_in}, rbond_init_r{rbond_init_r_in}, rbond_finl_r{rbond_finl_r_in},
    rangl_init_k{rangl_init_k_in}, rangl_finl_k{rangl_finl_k_in}, rangl_init_r{rangl_init_r_in},
    rangl_finl_r{rangl_finl_r_in}, rdihe_init_k{rdihe_init_k_in}, rdihe_finl_k{rdihe_finl_k_in},
    rdihe_init_r{rdihe_init_r_in}, rdihe_finl_r{rdihe_finl_r_in}, rposn_insr{rposn_insr_in},
    rbond_insr{rbond_insr_in}, rangl_insr{rangl_insr_in}, rdihe_insr{rdihe_insr_in},
    rposn_acc{rposn_acc_in}, rbond_acc{rbond_acc_in}, rangl_acc{rangl_acc_in},
    rdihe_acc{rdihe_acc_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2>
SyNonbondedKit<T, T2>::SyNonbondedKit(const int nsys_in, const UnitCellType unit_cell_in,
                                      const int nnbwu_in, const int* nbwu_abstracts_in,
                                      const uint2* nbwu_insr_in, const int* atom_offsets_in,
                                      const int* atom_counts_in, const T coulomb_in,
                                      const ImplicitSolventModel igb_in,
                                      const int neck_table_size_in, const T dielectric_in,
                                      const T kappa_in, const T saltcon_in, const T gb_offset_in,
                                      const T gb_neckscale_in, const T gb_neckcut_in,
                                      const T* charge_in, const int* q_idx_in,
                                      const T* q_params_in, const int* lj_idx_in,
                                      const int* n_lj_types_in, const int* ljabc_offsets_in,
                                      const T2* ljab_coeff_in, const T* ljc_coeff_in,
                                      const T* lj_sigma_in, const int* neck_gb_idx_in,
                                      const T* pb_radii_in, const T* gb_screen_in,
                                      const T* gb_alpha_in, const T* gb_beta_in,
                                      const T* gb_gamma_in, const T2* neck_limits_in) :
    nsys{nsys_in}, unit_cell{unit_cell_in}, nnbwu{nnbwu_in}, nbwu_abstracts{nbwu_abstracts_in},
    nbwu_insr{nbwu_insr_in}, atom_offsets{atom_offsets_in}, atom_counts{atom_counts_in},
    coulomb{coulomb_in}, igb{igb_in}, neck_table_size{neck_table_size_in},
    dielectric{dielectric_in}, kappa{kappa_in}, saltcon{saltcon_in}, gb_offset{gb_offset_in},
    gb_neckscale{gb_neckscale_in}, gb_neckcut{gb_neckcut_in}, charge{charge_in}, q_idx{q_idx_in},
    q_params{q_params_in}, lj_idx{lj_idx_in}, n_lj_types{n_lj_types_in},
    ljabc_offsets{ljabc_offsets_in}, ljab_coeff{ljab_coeff_in}, ljc_coeff{ljc_coeff_in},
    lj_sigma{lj_sigma_in}, neck_gb_idx{neck_gb_idx_in}, pb_radii{pb_radii_in},
    gb_screen{gb_screen_in}, gb_alpha{gb_alpha_in}, gb_beta{gb_beta_in}, gb_gamma{gb_gamma_in},
    neck_limits{neck_limits_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2, typename T4>
SyAtomUpdateKit<T, T2, T4>::SyAtomUpdateKit(const T* masses_in, const T* inv_masses_in,
                                            const int largest_group_in, const T4* vs_params_in,
                                            const T4* settle_geom_in, const T2* settle_mass_in,
                                            const T2* cnst_grp_params_in,
                                            const uint2* vste_insr_in, const uint2* sett_insr_in,
                                            const uint2* cnst_insr_in, const uint2* vwu_manip_in) :
    masses{masses_in}, inv_masses{inv_masses_in}, largest_group{largest_group_in},
    vs_params{vs_params_in}, settle_geom{settle_geom_in}, settle_mass{settle_mass_in},
    cnst_grp_params{cnst_grp_params_in}, vste_insr{vste_insr_in}, sett_insr{sett_insr_in},
    cnst_insr{cnst_insr_in}, vwu_manip{vwu_manip_in}
{}

} // namespace synthesis
} // namespace stormm
