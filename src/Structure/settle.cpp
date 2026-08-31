#include "copyright.h"
#include "settle.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
void settleVelocities(PhaseSpace *ps, const AtomGraph *ag, const PrecisionModel prec) {
  PhaseSpaceWriter psw = ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const ConstraintKit<double> cnst = ag->getDoublePrecisionConstraintKit();
      settleVelocities<double>(&psw, cnst);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const ConstraintKit<float> cnst = ag->getSinglePrecisionConstraintKit();
      settleVelocities<float>(&psw, cnst);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void settleVelocities(PhaseSpace *ps, const AtomGraph &ag, const PrecisionModel prec) {
  settleVelocities(ps, ag.getSelfPointer(), prec);
}

//-------------------------------------------------------------------------------------------------
void settleVelocities(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                      const PrecisionModel prec) {
  PsSynthesisWriter poly_psw = poly_ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag->getDoublePrecisionValenceKit();
      const SyAtomUpdateKit<double,
                            double2,
                            double4_16a> poly_auk = poly_ag->getDoublePrecisionAtomUpdateKit();
      settleVelocities(&poly_psw, poly_vk, poly_auk);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag->getSinglePrecisionValenceKit();
      const SyAtomUpdateKit<float,
                            float2, float4> poly_auk = poly_ag->getSinglePrecisionAtomUpdateKit();
      settleVelocities(&poly_psw, poly_vk, poly_auk);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void settleVelocities(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                     const PrecisionModel prec) {
  settleVelocities(poly_ps, poly_ag.getSelfPointer(), prec);
}

//-------------------------------------------------------------------------------------------------
void settlePositions(PhaseSpace *ps, const AtomGraph *ag, const double dt,
                     const PrecisionModel prec) {
  PhaseSpaceWriter psw = ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const ConstraintKit<double> cnst = ag->getDoublePrecisionConstraintKit();
      settlePositions<double, double3>(&psw, cnst, dt);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const ConstraintKit<float> cnst = ag->getSinglePrecisionConstraintKit();
      settlePositions<float, float3>(&psw, cnst, dt);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void settlePositions(PhaseSpace *ps, const AtomGraph &ag, const double dt,
                     const PrecisionModel prec) {
  settlePositions(ps, ag.getSelfPointer(), dt, prec);
}

//-------------------------------------------------------------------------------------------------
void settlePositions(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                     const double dt, const PrecisionModel prec) {
  PsSynthesisWriter poly_psw = poly_ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag->getDoublePrecisionValenceKit();
      const SyAtomUpdateKit<double,
                            double2,
                            double4_16a> poly_auk = poly_ag->getDoublePrecisionAtomUpdateKit();
      settlePositions<double, double2, double3, double4_16a>(&poly_psw, poly_vk, poly_auk, dt);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag->getSinglePrecisionValenceKit();
      const SyAtomUpdateKit<float,
                            float2, float4> poly_auk = poly_ag->getSinglePrecisionAtomUpdateKit();
      settlePositions<float, float2, float3, float4>(&poly_psw, poly_vk, poly_auk, dt);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void settlePositions(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                     const double dt, const PrecisionModel prec) {
  settlePositions(poly_ps, poly_ag.getSelfPointer(), dt, prec);
}

//-------------------------------------------------------------------------------------------------
void idealizeSettleReference(PhaseSpace *ps, const AtomGraph *ag, const PrecisionModel prec) {
  PhaseSpaceWriter psw = ps->data();
  switch (prec)	{
  case PrecisionModel::DOUBLE:
    idealizeSettleReference<double>(&psw, ag->getDoublePrecisionConstraintKit());
    break;
  case PrecisionModel::SINGLE:
    idealizeSettleReference<float>(&psw, ag->getSinglePrecisionConstraintKit());
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void idealizeSettleReference(PhaseSpace *ps, const AtomGraph &ag, const PrecisionModel prec) {
  idealizeSettleReference(ps, ag.getSelfPointer(), prec);
}
  
} // namespace structure
} // namespace stormm

