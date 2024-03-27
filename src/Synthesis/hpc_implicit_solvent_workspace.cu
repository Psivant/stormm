// -*-c++-*-
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "implicit_solvent_workspace.h"

namespace stormm {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
// Set accumulators in an ImplicitSolventWorkspace to zero to prepare for a new calculation.
//
// Arguments:
//   psi:              Accumulators for effective GB radii
//   psi_ovrf:         Overflow accumulators for effective GB radii
//   sum_deijda:       Accumulators for the GB radii derivatives
//   sum_deijda_ovrf:  Overflow accumulators for the GB radii derivatives
//   natom:            The total (padded) number of atoms across all systems in the synthesis
//                     served by this kernel
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kInitISWorkspace(llint* psi, int* psi_ovrf, llint* sum_deijda, int* sum_deijda_ovrf,
                 const int natom) {
  for (int pos = threadIdx.x + (blockIdx.x * blockDim.x); pos < natom;
       pos += blockDim.x * gridDim.x) {
    psi[pos] = 0LL;
    sum_deijda[pos] = 0LL;
    psi_ovrf[pos] = 0;
    sum_deijda_ovrf[pos] = 0;
  }
}

//-------------------------------------------------------------------------------------------------
void ImplicitSolventWorkspace::launchInitialization(const GpuDetails &gpu,
                                                    const CoordinateCycle orientation) {
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  const int nsmp = gpu.getSMPCount();
  switch (orientation) {
  case CoordinateCycle::BLACK:
    kInitISWorkspace<<<nsmp, large_block_size>>>(alt_psi.data(devc_tier),
                                                 alt_psi_overflow.data(devc_tier),
                                                 alt_sum_deijda.data(devc_tier),
                                                 alt_sum_deijda_overflow.data(devc_tier),
                                                 padded_atom_count);
    break;
  case CoordinateCycle::WHITE:
    kInitISWorkspace<<<nsmp, large_block_size>>>(psi.data(devc_tier), psi_overflow.data(devc_tier),
                                                 sum_deijda.data(devc_tier),
                                                 sum_deijda_overflow.data(devc_tier),
                                                 padded_atom_count);
    break;
  }
}

} // namespace synthesis
} // namespace stormm
