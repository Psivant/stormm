#include "copyright.h"
#include "kinetic.h"

namespace stormm {
namespace mm {

using trajectory::PhaseSpaceReader;
using symbols::boltzmann_constant;
  
//-------------------------------------------------------------------------------------------------
void evalKineticEnergy(const PhaseSpaceWriter psw, ScoreCard *sc, const ChemicalDetailsKit &cdk,
                       const int system_index) {
  evalKineticEnergy<double, double>(psw.vxalt, psw.vyalt, psw.vzalt, sc, cdk, system_index, 1.0);
}

//-------------------------------------------------------------------------------------------------
void evalKineticEnergy(const PhaseSpace *ps, CoordinateCycle orientation, ScoreCard *sc,
                       const AtomGraph *ag, int system_index) {
  const PhaseSpaceReader psr = ps->data(orientation);
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  evalKineticEnergy<double, double>(psr.vxalt, psr.vyalt, psr.vzalt, sc, cdk, system_index, 1.0);
}

//-------------------------------------------------------------------------------------------------
void evalKineticEnergy(const PhaseSpace &ps, CoordinateCycle orientation, ScoreCard *sc,
                       const AtomGraph &ag, int system_index) {
  const PhaseSpaceReader psr = ps.data(orientation);
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  evalKineticEnergy<double, double>(psr.vxalt, psr.vyalt, psr.vzalt, sc, cdk, system_index, 1.0);
}

//-------------------------------------------------------------------------------------------------
void evalKineticEnergy(const PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                       const int system_index) {
  evalKineticEnergy(ps, ps->getCyclePosition(), sc, ag, system_index);
}

//-------------------------------------------------------------------------------------------------
void evalKineticEnergy(const PhaseSpace &ps, ScoreCard *sc, const AtomGraph &ag,
                       const int system_index) {
  evalKineticEnergy(ps, ps.getCyclePosition(), sc, ag, system_index);
}

//-------------------------------------------------------------------------------------------------
void computeTemperature(ScoreCard *sc, const ChemicalDetailsKit &cdk, const bool cnst,
                        const int index) {
  const double nrg_scale = sc->getEnergyScalingFactor<double>();
  const double ke = sc->reportInstantaneousStates(StateVariable::KINETIC, index);
  const double dof = (cnst) ? cdk.cnst_dof : cdk.free_dof;
  const llint acc = llround(2.0 * ke * nrg_scale / (dof * boltzmann_constant));
  sc->contribute(StateVariable::TEMPERATURE_ALL, acc, index);
}

//-------------------------------------------------------------------------------------------------
void computeTemperature(ScoreCard *sc, const AtomGraph *ag, const bool cnst, const int index) {
  const double nrg_scale = sc->getEnergyScalingFactor<double>();
  const double ke = sc->reportInstantaneousStates(StateVariable::KINETIC, index);
  const double dof = (cnst) ? ag->getConstrainedDegreesOfFreedom() : ag->getDegreesOfFreedom();
  const llint acc = llround(2.0 * ke * nrg_scale / (dof * boltzmann_constant));
  sc->contribute(StateVariable::TEMPERATURE_ALL, acc, index);
}

//-------------------------------------------------------------------------------------------------
void computeTemperature(ScoreCard *sc, const AtomGraph &ag, const bool cnst, const int index) {
  computeTemperature(sc, ag.getSelfPointer(), cnst, index);
}

//-------------------------------------------------------------------------------------------------
double computeTemperature(const PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                          const ApplyConstraints use_cnst, const int system_index,
                          const PrecisionModel prec) {
  if (system_index < 0 || system_index >= poly_ps->getSystemCount()) {
    rtErr("System index " + std::to_string(system_index) + " is invalid for a synthesis of " +
          std::to_string(poly_ps->getSystemCount()) + " systems.", "computeTemperature");
  }
  PsSynthesisReader poly_psr = poly_ps->data();
  const AtomGraph *sys_ag = poly_ps->getSystemTopologyPointer(system_index);
  int ndof;
  switch (use_cnst) {
  case ApplyConstraints::YES:
    ndof = sys_ag->getConstrainedDegreesOfFreedom();
    break;
  case ApplyConstraints::NO:
    ndof = sys_ag->getDegreesOfFreedom();
    break;
  }
  const size_t start_idx = poly_psr.atom_starts[system_index];
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyAtomUpdateKit<double,
                            double2,
                            double4> poly_auk = poly_ag->getDoublePrecisionAtomUpdateKit();
      return computeTemperature<llint,
                                double>(&poly_psr.xvel[start_idx], &poly_psr.yvel[start_idx],
                                        &poly_psr.zvel[start_idx], &poly_psr.xvel_ovrf[start_idx],
                                        &poly_psr.yvel_ovrf[start_idx],
                                        &poly_psr.zvel_ovrf[start_idx], poly_auk.masses,
                                        poly_psr.atom_counts[system_index], ndof, pow(2.0, 32),
                                        poly_psr.inv_vel_scale);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyAtomUpdateKit<float,
                            float2, float4> poly_auk = poly_ag->getSinglePrecisionAtomUpdateKit();
      return computeTemperature<llint,
                                float>(&poly_psr.xvel[start_idx], &poly_psr.yvel[start_idx],
                                       &poly_psr.zvel[start_idx], nullptr, nullptr, nullptr,
                                       poly_auk.masses, poly_psr.atom_counts[system_index], ndof,
                                       pow(2.0, 32), poly_psr.inv_vel_scale);
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double computeTemperature(const PhaseSpaceSynthesis &poly_ps, const AtomGraphSynthesis &poly_ag,
                          const ApplyConstraints use_cnst, const int system_index,
                          const PrecisionModel prec) {
  return computeTemperature(poly_ps.getSelfPointer(), poly_ag.getSelfPointer(), use_cnst,
                            system_index, prec);
}

} // namespace mm
} // namespace stormm
