// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename T>
ConvolutionKit<T>::ConvolutionKit(const int system_count_in, const T ew_coeff_in,
                                  const int* sys_offsets_in, const T* self_ecorr_in,
                                  const T* bmesh_a_in, const T* bmesh_b_in, const T* bmesh_c_in,
                                  const T* mval_a_in, const T* mval_b_in, const T* mval_c_in,
                                  const T* msval_a_in, const T* msval_b_in, const T* msval_c_in,
                                  const T* cmesh_a_in, const T* cmesh_b_in, const T* cmesh_c_in) :
    system_count{system_count_in}, ew_coeff{ew_coeff_in}, sys_offsets{sys_offsets_in},
    self_ecorr{self_ecorr_in}, bmesh_a{bmesh_a_in}, bmesh_b{bmesh_b_in}, bmesh_c{bmesh_c_in},
    mval_a{mval_a_in}, mval_b{mval_b_in}, mval_c{mval_c_in}, msval_a{msval_a_in},
    msval_b{msval_b_in}, msval_c{msval_c_in}, cmesh_a{cmesh_a_in}, cmesh_b{cmesh_b_in},
    cmesh_c{cmesh_c_in}
{}

} // namespace energy
} // namespace stormm
