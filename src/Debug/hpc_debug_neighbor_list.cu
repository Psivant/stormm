// -*-c++-*-
#include "copyright.h"
#include "hpc_debug_neighbor_list.cuh"

namespace stormm {
namespace mm {

using card::GpuDetails;

//-------------------------------------------------------------------------------------------------
void debugNLComposition(const int step, const GpuDetails &gpu) {

  // Find the neighbor list data types
  const size_t nl_tcoord  = dyna_tk.getNeighborListCoordinateType();
  const size_t nl_tacc    = dyna_tk.getNeighborListAccumulationType();
  const size_t nl_tcalc   = dyna_tk.getNeighborListCalculationType();
  const size_t nl_tcoord4 = dyna_tk.getNeighborListTupleType();

  // Based on the neighbor list layout, prepare the downloads
  std::vector<NonbondedTheme> nl_themes;
  switch (dyna_tk.getNeighborListLayout()) {
  case NeighborListKind::MONO:
    nl_themes = { NonbondedTheme::ALL };
    break;
  case NeighborListKind::DUAL:
    nl_themes = { NonbondedTheme::ELECTROSTATIC, NonbondedTheme::VAN_DER_WAALS };
    break;
  }

  // Loop over all applicable themes and copy the neighbor list contents to workspaces
  for (size_t i = 0; i < nl_themes.size(); i++) {
    if (nl_tcoord == float_type_index) {
      if (nl_tacc == int_type_index) {
        if (nl_tcalc == float_type_index) {
          debugNLComposition<float, int, float, float4>(step, nl_themes[i], gpu);
        }
        else if (nl_tcalc == double_type_index) {
          debugNLComposition<float, int, double, float4>(step, nl_themes[i], gpu);
        }
      }
      else if (nl_tacc == llint_type_index) {
        if (nl_tcalc == float_type_index) {
          debugNLComposition<float, llint, float, float4>(step, nl_themes[i], gpu);
        }
        else if (nl_tcalc == double_type_index) {
          debugNLComposition<float, llint, double, float4>(step, nl_themes[i], gpu);
        }
      }
    }
    else if (nl_tcoord == double_type_index) {
      if (nl_tacc == int_type_index) {
        if (nl_tcalc == float_type_index) {
          debugNLComposition<double, int, float, double4_16a>(step, nl_themes[i], gpu);
        }
        else if (nl_tcalc == double_type_index) {
          debugNLComposition<double, int, double, double4_16a>(step, nl_themes[i], gpu);
        }
      }
      else if (nl_tacc == llint_type_index) {
        if (nl_tcalc == float_type_index) {
          debugNLComposition<double, llint, float, double4_16a>(step, nl_themes[i], gpu);
        }
        else if (nl_tcalc == double_type_index) {
          debugNLComposition<double, llint, double, double4_16a>(step, nl_themes[i], gpu);
        }
      }
    }
  }
}

} // namespace mm
} // namespace stormm
