// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace mm {

//-------------------------------------------------------------------------------------------------
template <typename T>
MMControlKit<T>::MMControlKit(const int step_in, const int sd_cycles_in, const int max_cycles_in,
                              const T dt_in, const T rattle_tol_in, const T initial_step_in,
                              int* vwu_progress_in, int* vupt_progress_in, int* vcns_progress_in,
                              int* pupt_progress_in, int* gcns_progress_in, int* nbwu_progress_in,
                              int* pmewu_progress_in, int* gbrwu_progress_in,
                              int* gbdwu_progress_in, int* gtwu_progress_in, int* scwu_progress_in,
                              int* rdwu_progress_in) :
    step{step_in}, sd_cycles{sd_cycles_in}, max_cycles{max_cycles_in}, dt{dt_in},
    rattle_tol{rattle_tol_in}, initial_step{initial_step_in}, vwu_progress{vwu_progress_in},
    vupt_progress{vupt_progress_in}, vcns_progress{vcns_progress_in},
    pupt_progress{pupt_progress_in}, gcns_progress{gcns_progress_in},
    nbwu_progress{nbwu_progress_in}, pmewu_progress{pmewu_progress_in},
    gbrwu_progress{gbrwu_progress_in}, gbdwu_progress{gbdwu_progress_in},
    gtwu_progress{gtwu_progress_in}, scwu_progress{scwu_progress_in},
    rdwu_progress{rdwu_progress_in}
{}

} // namespace mm
} // namespace stormm
