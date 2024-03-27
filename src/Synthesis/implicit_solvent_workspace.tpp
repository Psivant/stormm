// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
template <typename T>
ISWorkspaceKit<T>::ISWorkspaceKit(const int fp_bits_in, llint* psi_in, int* psi_ovrf_in,
                                  llint* sum_deijda_in, int* sum_deijda_ovrf_in, llint* alt_psi_in,
                                  int* alt_psi_ovrf_in, llint* alt_sum_deijda_in,
                                  int* alt_sum_deijda_ovrf_in) :
    fp_bits{fp_bits_in}, fp_scale{static_cast<T>(pow(2.0, fp_bits_in))},
    inv_fp_scale{static_cast<T>(1.0 / fp_scale)}, psi{psi_in}, psi_ovrf{psi_ovrf_in},
    sum_deijda{sum_deijda_in}, sum_deijda_ovrf{sum_deijda_ovrf_in}, alt_psi{alt_psi_in},
    alt_psi_ovrf{alt_psi_ovrf_in}, alt_sum_deijda{alt_sum_deijda_in},
    alt_sum_deijda_ovrf{alt_sum_deijda_ovrf_in}
{}

} // namespace synthesis
} // namespace stormm
