// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace generalized_born_defaults {

//-------------------------------------------------------------------------------------------------
template <typename T>
NeckGeneralizedBornKit<T>::NeckGeneralizedBornKit(const int table_size_in, const T neck_cut_in,
                                                  const T kscale_in, const T* max_separation_in,
                                                  const T* max_value_in) :
    table_size{table_size_in}, neck_cut{neck_cut_in}, kscale{kscale_in},
    max_separation{max_separation_in}, max_value{max_value_in}
{}

} // namespace generalized_born_defaults
} // namespace stormm
