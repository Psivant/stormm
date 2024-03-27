// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename T>
ImplicitSolventRecipe<T>::ImplicitSolventRecipe(const ImplicitSolventKit<T> &isk,
                                                const NeckGeneralizedBornKit<T> &ngbk) :
    natom{isk.natom},
    igb{isk.igb},
    table_size{ngbk.table_size},
    dielectric{isk.dielectric},
    kappa{(isk.saltcon > static_cast<T>(constants::tiny)) ?
          static_cast<T>(sqrt(default_salt_kappa_dependence * isk.saltcon)) : static_cast<T>(0.0)},
    gb_offset{(igb == ImplicitSolventModel::NECK_GB_II) ?
              static_cast<T>(default_neck_ii_gb_radii_offset) :
              static_cast<T>(default_gb_radii_offset)},
    gb_neckscale{(igb == ImplicitSolventModel::NECK_GB_II) ?
                 static_cast<T>(default_gb_neck_ii_scale) : static_cast<T>(default_gb_neck_scale)},
    gb_neckcut{ngbk.neck_cut},
    neck_gb_idx{isk.neck_gb_idx},
    pb_radii{isk.pb_radii},
    gb_screen{isk.gb_screen},
    gb_alpha{isk.gb_alpha},
    gb_beta{isk.gb_beta},
    gb_gamma{isk.gb_gamma},
    neck_max_sep{ngbk.max_separation},
    neck_max_val{ngbk.max_value}
{}

} // namespace energy
} // namespace stormm
