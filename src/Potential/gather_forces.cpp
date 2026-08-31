#include "copyright.h"
#include "gather_forces.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
void gatherForces(PMIGrid *pm, const AtomGraphSynthesis *poly_ag) {

  // In order to unpack the CellGrid underlying the PMIGrid, its templated parameters must be
  // understood.
  const size_t cg_tmat  = pm->getCellGridMatrixTypeID();
  const size_t cg_tacc  = pm->getCellGridAccumulatorTypeID();
  const size_t cg_tcalc = pm->getCellGridCalculationTypeID();
  if (cg_tmat == double_type_index) { 
    unrollGatherForcesCall<double, double4_16a>(pm, cg_tacc, cg_tcalc, poly_ag);
  }
  else if (cg_tmat == float_type_index) {
    unrollGatherForcesCall<float, float4>(pm, cg_tacc, cg_tcalc, poly_ag);
  }
  else if (cg_tmat == llint_type_index) {
    unrollGatherForcesCall<llint, llint4>(pm, cg_tacc, cg_tcalc, poly_ag);
  }
  else if (cg_tmat == int_type_index) {
    unrollGatherForcesCall<int, int4>(pm, cg_tacc, cg_tcalc, poly_ag);
  }
}

//-------------------------------------------------------------------------------------------------
void gatherForces(PMIGrid *pm, const AtomGraphSynthesis &poly_ag) {
  gatherForces(pm, poly_ag.getSelfPointer());
}

} // namespace energy
} // namespace stormm
