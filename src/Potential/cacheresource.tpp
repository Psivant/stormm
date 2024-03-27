// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename T>
CacheResourceKit<T>::CacheResourceKit(const int max_blocks_in, const int max_atoms_in,
                                      llint* xcrd_in, llint* ycrd_in, llint* zcrd_in,
                                      llint* xvel_in, llint* yvel_in, llint* zvel_in,
                                      int* xcrd_ovrf_in, int* ycrd_ovrf_in, int* zcrd_ovrf_in,
                                      int* xvel_ovrf_in, int* yvel_ovrf_in, int* zvel_ovrf_in,
                                      int* xfrc_ovrf_in, int* yfrc_ovrf_in, int* zfrc_ovrf_in,
                                      T* charges_in, int* lj_idx_in) :
    max_blocks{max_blocks_in}, max_atoms{max_atoms_in}, xcrd{xcrd_in}, ycrd{ycrd_in},
    zcrd{zcrd_in}, xvel{xvel_in}, yvel{yvel_in}, zvel{zvel_in}, xcrd_ovrf{xcrd_ovrf_in},
    ycrd_ovrf{ycrd_ovrf_in}, zcrd_ovrf{zcrd_ovrf_in}, xvel_ovrf{xvel_ovrf_in},
    yvel_ovrf{yvel_ovrf_in}, zvel_ovrf{zvel_ovrf_in}, xfrc_ovrf{xfrc_ovrf_in},
    yfrc_ovrf{yfrc_ovrf_in}, zfrc_ovrf{zfrc_ovrf_in}, charges{charges_in}, lj_idx{lj_idx_in}
{}

} // namespace energy
} // namespace stormm
