// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace mm {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void debugNLComposition(const int step) {
  switch (dyna_tk.getNeighborListLayout()) {
  case NeighborListKind::MONO:
    {
      CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> *test_cg;
      test_cg = dyna_tk.getNLWorkspacePointer<Tcoord, Tacc, Tcalc, Tcoord4>(NonbondedTheme::ALL);
    }
    break;
  case NeighborListKind::DUAL:
    {
      CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> *test_qqcg, *test_ljcg;
      test_qqcg = dyna_tk.getNLWorkspacePointer<Tcoord, Tacc,
                                                Tcalc, Tcoord4>(NonbondedTheme::ELECTROSTATIC);
      test_ljcg = dyna_tk.getNLWorkspacePointer<Tcoord, Tacc,
                                                Tcalc, Tcoord4>(NonbondedTheme::VAN_DER_WAALS);
    }
    break;
  }
}
  
} // namespace mm 
} // namespace stormm
