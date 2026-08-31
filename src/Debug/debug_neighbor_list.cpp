#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "MolecularMechanics/dynamics_intervention.h"

namespace stormm {
namespace mm {

using card::GpuDetails;
using card::HybridTargetLevel;
  
//-------------------------------------------------------------------------------------------------
void execNLCompositionDebug(const int step, const int index, const HybridTargetLevel tier) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      // Find the neighbor list data types
      const size_t nl_tcoord  = dyna_tk.getNeighborListCoordinateType();
      const size_t nl_tacc    = dyna_tk.getNeighborListAccumulationType();
      const size_t nl_tcalc   = dyna_tk.getNeighborListCalculationType();
      const size_t nl_tcoord4 = dyna_tk.getNeighborListTupleType();

      // Branch based on the types, calling the appropriate templated variant of the CPU-based
      // checking function.
      if (nl_tcoord == float_type_index) {
        if (nl_tacc == int_type_index) {
          if (nl_tcalc == float_type_index) {
            debugNLComposition<float, int, float, float4>(step);
          }
          else if (nl_tcalc == double_type_index) {
            debugNLComposition<float, int, double, float4>(step);
          }
        }
        else if (nl_tacc == llint_type_index) {
          if (nl_tcalc == float_type_index) {
            debugNLComposition<float, llint, float, float4>(step);
          }
          else if (nl_tcalc == double_type_index) {
            debugNLComposition<float, llint, double, float4>(step);
          }
        }
      }
      else if (nl_tcoord == double_type_index) {
        if (nl_tacc == int_type_index) {
          if (nl_tcalc == float_type_index) {
            debugNLComposition<double, int, float, double4_16a>(step);
          }
          else if (nl_tcalc == double_type_index) {
            debugNLComposition<double, int, double, double4_16a>(step);
          }
        }
        else if (nl_tacc == llint_type_index) {
          if (nl_tcalc == float_type_index) {
            debugNLComposition<double, llint, float, double4_16a>(step);
          }
          else if (nl_tcalc == double_type_index) {
            debugNLComposition<double, llint, double, double4_16a>(step);
          }
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const GpuDetails &gpu = dyna_tk.getGpuDetails();
      debugNLComposition(step, gpu);
    }
    break;
#endif
  }
}

} // namespace mm
} // namespace stormm
