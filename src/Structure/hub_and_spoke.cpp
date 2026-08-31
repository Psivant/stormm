#include "copyright.h"
#include "hub_and_spoke.h"

namespace stormm {
namespace structure {

using trajectory::getNextCyclePosition;

//-------------------------------------------------------------------------------------------------
void shakePositions(PhaseSpace *ps, const AtomGraph *ag, const PrecisionModel prec,
                     const double dt, const double tol, const int max_iter,
                     const RattleMethod style) {

  // The writeable abstract is obtained with respect to whichever point in the time cycle the
  // original coordinate object has as "current."
  PhaseSpaceWriter psw = ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    shakePositions<double>(&psw, ag->getDoublePrecisionConstraintKit(), dt, tol, max_iter, style);
    break;
  case PrecisionModel::SINGLE:
    shakePositions<float>(&psw, ag->getSinglePrecisionConstraintKit(), dt, tol, max_iter, style);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void shakePositions(PhaseSpace *ps, const AtomGraph &ag, const PrecisionModel prec,
                     const double dt, const double tol, const int max_iter,
                     const RattleMethod style) {
  shakePositions(ps, ag.getSelfPointer(), prec, dt, tol, max_iter, style);
}

//-------------------------------------------------------------------------------------------------
void shakePositions(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                     const PrecisionModel prec, const double dt, const double tol,
                     const int max_iter) {

  // The writeable abstract is obtained with respect to whichever point in the time cycle the
  // original coordinate object has as "current."
  PsSynthesisWriter poly_psw = poly_ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    shakePositions<double,
                   double2, double4_16a>(&poly_psw, poly_ag->getDoublePrecisionValenceKit(),
                                     poly_ag->getDoublePrecisionAtomUpdateKit(), dt, tol,
                                     max_iter);
    break;
  case PrecisionModel::SINGLE:
    shakePositions<float,
                   float2, float4>(&poly_psw, poly_ag->getSinglePrecisionValenceKit(),
                                   poly_ag->getSinglePrecisionAtomUpdateKit(), dt, tol, max_iter);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void shakePositions(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                     const PrecisionModel prec, const double dt, const double tol,
                     const int max_iter) {
  shakePositions(poly_ps, poly_ag.getSelfPointer(), prec, dt, tol, max_iter);
}
  
//-------------------------------------------------------------------------------------------------
void rattleVelocities(PhaseSpace *ps, const AtomGraph *ag, const PrecisionModel prec,
                      const double dt, const double tol, const int max_iter,
                      const RattleMethod style) {

  // The writeable abstract is obtained with respect to whichever point in the time cycle the
  // original coordinate object has as "current."
  PhaseSpaceWriter psw = ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    rattleVelocities<double>(&psw, ag->getDoublePrecisionConstraintKit(), dt, tol, max_iter,
                             style);
    break;
  case PrecisionModel::SINGLE:
    rattleVelocities<float>(&psw, ag->getSinglePrecisionConstraintKit(), dt, tol, max_iter, style);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void rattleVelocities(PhaseSpace *ps, const AtomGraph &ag, const PrecisionModel prec,
                      const double dt, const double tol, const int max_iter,
                      const RattleMethod style) {
  rattleVelocities(ps, ag.getSelfPointer(), prec, dt, tol, max_iter, style);
}

//-------------------------------------------------------------------------------------------------
void rattleVelocities(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                      const PrecisionModel prec, const double dt, const double tol,
                      const int max_iter) {
  PsSynthesisWriter poly_psw = poly_ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag->getDoublePrecisionValenceKit();
      const SyAtomUpdateKit<double,
                            double2,
                            double4_16a> poly_auk = poly_ag->getDoublePrecisionAtomUpdateKit();
      rattleVelocities<double, double2, double4_16a>(&poly_psw, poly_vk, poly_auk, dt, tol,
                                                     max_iter);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag->getSinglePrecisionValenceKit();
      const SyAtomUpdateKit<float,
                            float2, float4> poly_auk = poly_ag->getSinglePrecisionAtomUpdateKit();
      rattleVelocities<float, float2, float4>(&poly_psw, poly_vk, poly_auk, dt, tol, max_iter);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void rattleVelocities(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                      const PrecisionModel prec, const double dt, const double tol,
                      const int max_iter) {
  rattleVelocities(poly_ps, poly_ag.getSelfPointer(), prec, dt, tol, max_iter);
}

} // namespace structure
} // namespace stormm
