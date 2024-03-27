// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace restraints {

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2, typename T4>
RestraintKit<T, T2, T4>::RestraintKit(const int total_rst_in, const int nposn_in,
                                      const int nbond_in, const int nangl_in, const int ndihe_in,
                                      const bool time_dependence_in, const int* rposn_atoms_in,
                                      const int* rbond_i_atoms_in, const int* rbond_j_atoms_in,
                                      const int* rangl_i_atoms_in, const int* rangl_j_atoms_in,
                                      const int* rangl_k_atoms_in, const int* rdihe_i_atoms_in,
                                      const int* rdihe_j_atoms_in, const int* rdihe_k_atoms_in,
                                      const int* rdihe_l_atoms_in, const int* rposn_init_step_in,
                                      const int* rposn_finl_step_in, const int* rbond_init_step_in,
                                      const int* rbond_finl_step_in, const int* rangl_init_step_in,
                                      const int* rangl_finl_step_in, const int* rdihe_init_step_in,
                                      const int* rdihe_finl_step_in, const T2* rposn_init_keq_in,
                                      const T2* rposn_finl_keq_in, const T2* rposn_init_xy_in,
                                      const T2* rposn_finl_xy_in, const T* rposn_init_z_in,
                                      const T* rposn_finl_z_in, const T2* rbond_init_keq_in,
                                      const T2* rbond_finl_keq_in, const T2* rangl_init_keq_in,
                                      const T2* rangl_finl_keq_in, const T2* rdihe_init_keq_in,
                                      const T2* rdihe_finl_keq_in, const T4* rposn_init_r_in,
                                      const T4* rposn_finl_r_in, const T4* rbond_init_r_in,
                                      const T4* rbond_finl_r_in, const T4* rangl_init_r_in,
                                      const T4* rangl_finl_r_in, const T4* rdihe_init_r_in,
                                      const T4* rdihe_finl_r_in, const AtomGraph *ag_pointer_in) :
    total_rst{total_rst_in}, nposn{nposn_in}, nbond{nbond_in}, nangl{nangl_in}, ndihe{ndihe_in},
    time_dependence{time_dependence_in}, rposn_atoms{rposn_atoms_in},
    rbond_i_atoms{rbond_i_atoms_in}, rbond_j_atoms{rbond_j_atoms_in},
    rangl_i_atoms{rangl_i_atoms_in}, rangl_j_atoms{rangl_j_atoms_in},
    rangl_k_atoms{rangl_k_atoms_in}, rdihe_i_atoms{rdihe_i_atoms_in},
    rdihe_j_atoms{rdihe_j_atoms_in}, rdihe_k_atoms{rdihe_k_atoms_in},
    rdihe_l_atoms{rdihe_l_atoms_in}, rposn_init_step{rposn_init_step_in},
    rposn_finl_step{rposn_finl_step_in}, rbond_init_step{rbond_init_step_in},
    rbond_finl_step{rbond_finl_step_in}, rangl_init_step{rangl_init_step_in},
    rangl_finl_step{rangl_finl_step_in}, rdihe_init_step{rdihe_init_step_in},
    rdihe_finl_step{rdihe_finl_step_in}, rposn_init_keq{rposn_init_keq_in},
    rposn_finl_keq{rposn_finl_keq_in}, rposn_init_xy{rposn_init_xy_in},
    rposn_finl_xy{rposn_finl_xy_in}, rposn_init_z{rposn_init_z_in},
    rposn_finl_z{rposn_finl_z_in}, rbond_init_keq{rbond_init_keq_in},
    rbond_finl_keq{rbond_finl_keq_in}, rangl_init_keq{rangl_init_keq_in},
    rangl_finl_keq{rangl_finl_keq_in}, rdihe_init_keq{rdihe_init_keq_in},
    rdihe_finl_keq{rdihe_finl_keq_in}, rposn_init_r{rposn_init_r_in},
    rposn_finl_r{rposn_finl_r_in}, rbond_init_r{rbond_init_r_in},
    rbond_finl_r{rbond_finl_r_in}, rangl_init_r{rangl_init_r_in},
    rangl_finl_r{rangl_finl_r_in}, rdihe_init_r{rdihe_init_r_in},
    rdihe_finl_r{rdihe_finl_r_in}, ag_pointer{ag_pointer_in}
{}

} // namespace restraints
} // namespace stormm
