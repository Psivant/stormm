// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace mm {

//-------------------------------------------------------------------------------------------------
template <typename T>
MMControlKit<T>::MMControlKit(const int step_in, const int sd_cycles_in, const int max_cycles_in,
                              const T initial_step_in, const int nt_warp_mult_in,
                              const T elec_cut_in, const T vdw_cut_in, int* vwu_progress_in,
                              int* vupt_progress_in, int* vcns_progress_in, int* pupt_progress_in,
                              int* gcns_progress_in, int* nbwu_progress_in, int* pmewu_progress_in,
                              int* gbrwu_progress_in, int* gbdwu_progress_in,
                              int* gtwu_progress_in, int* scwu_progress_in,
                              int* rdwu_progress_in) :
    step{step_in}, sd_cycles{sd_cycles_in}, max_cycles{max_cycles_in},
    initial_step{initial_step_in}, nt_warp_mult{nt_warp_mult_in}, elec_cut{elec_cut_in},
    vdw_cut{vdw_cut_in}, elec_cut_sq{elec_cut_in * elec_cut_in},
    vdw_cut_sq{vdw_cut_in * vdw_cut_in}, vwu_progress{vwu_progress_in},
    vupt_progress{vupt_progress_in}, vcns_progress{vcns_progress_in},
    pupt_progress{pupt_progress_in}, gcns_progress{gcns_progress_in},
    nbwu_progress{nbwu_progress_in}, pmewu_progress{pmewu_progress_in},
    gbrwu_progress{gbrwu_progress_in}, gbdwu_progress{gbdwu_progress_in},
    gtwu_progress{gtwu_progress_in}, scwu_progress{scwu_progress_in},
    rdwu_progress{rdwu_progress_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void MolecularMechanicsControls::setNTWarpMultiplicity(const CellGrid<T, Tacc, Tcalc, T4> *cg_a,
                                                       const CellGrid<T, Tacc, Tcalc, T4> *cg_b,
                                                       const GpuDetails &gpu) {
  const int total_cells = (cg_b == nullptr) ? cg_a->getTotalCellCount() :
                                              cg_a->getTotalCellCount() +
                                              cg_b->getTotalCellCount();
  const int gpu_warps = 40 * gpu.getSMPCount();
  if (total_cells > 3 * gpu_warps) {
    nt_warp_multiplicity = 1;
  }
  else {
    if (static_cast<double>(gpu_warps * 2) <= total_cells) {
      nt_warp_multiplicity = std::min(8, total_cells / gpu_warps);
    }
    else if (static_cast<double>(gpu_warps * 3) <= total_cells * 2) {
      nt_warp_multiplicity = 3;
    }
    else {
      nt_warp_multiplicity = 1;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
void MolecularMechanicsControls::setNTWarpMultiplicity(const CellGrid<T, Tacc, Tcalc, T4> *cg_a,
                                                       const GpuDetails &gpu) {
  setNTWarpMultiplicity(cg_a, nullptr, gpu);
}

} // namespace mm
} // namespace stormm
