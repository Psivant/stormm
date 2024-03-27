// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace mm {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tmass, typename Tcalc>
llint evalKineticEnergy(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                        const int* xvel_ovrf, const int* yvel_ovrf, const int* zvel_ovrf,
                        const Tmass* masses, const int natom, const Tcalc nrg_scale_factor,
                        const Tcalc inv_vel_scale) {
  llint result = 0LL;
  const Tcalc my_gafs_to_kcal = gafs_to_kcal;
  if (isSignedIntegralScalarType<Tcoord>()) {
    if (xvel_ovrf == nullptr && yvel_ovrf == nullptr && zvel_ovrf == nullptr) {
      for (int i = 0; i < natom; i++) {
        const Tcalc mss_i = static_cast<Tcalc>(masses[i]) * static_cast<Tcalc>(0.5);
        const Tcalc vx = static_cast<Tcalc>(xvel[i]) * inv_vel_scale;
        const Tcalc vy = static_cast<Tcalc>(yvel[i]) * inv_vel_scale;
        const Tcalc vz = static_cast<Tcalc>(zvel[i]) * inv_vel_scale;
        const Tcalc contrib_i = mss_i * ((vx * vx) + (vy * vy) + (vz * vz)) * my_gafs_to_kcal;
        result += llround(contrib_i * nrg_scale_factor);
      }
    }
    else {
      for (int i = 0; i < natom; i++) {
        const Tcalc mss_i = static_cast<Tcalc>(masses[i]) * static_cast<Tcalc>(0.5);
        const Tcalc vx = hostInt95ToDouble(xvel[i], xvel_ovrf[i]) * inv_vel_scale;
        const Tcalc vy = hostInt95ToDouble(yvel[i], yvel_ovrf[i]) * inv_vel_scale;
        const Tcalc vz = hostInt95ToDouble(zvel[i], zvel_ovrf[i]) * inv_vel_scale;
        const Tcalc contrib_i = mss_i * ((vx * vx) + (vy * vy) + (vz * vz)) * my_gafs_to_kcal;
        result += llround(contrib_i * nrg_scale_factor);
      }
    }
  }
  else {
    for (int i = 0; i < natom; i++) {
      const Tcalc mss_i = static_cast<Tcalc>(masses[i]) * static_cast<Tcalc>(0.5);
      const Tcalc contrib_i = mss_i * my_gafs_to_kcal *
                              ((xvel[i] * xvel[i]) + (yvel[i] * yvel[i]) + (zvel[i] * zvel[i]));
      result += llround(contrib_i * nrg_scale_factor);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void evalKineticEnergy(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                       ScoreCard *sc, const ChemicalDetailsKit &cdk, const int system_index,
                       const Tcalc inv_vel_scale) {
  const Tcalc nrg_scale_factor = sc->getEnergyScalingFactor<Tcalc>();
  const llint acc = evalKineticEnergy<Tcoord, double, Tcalc>(xvel, yvel, zvel, nullptr, nullptr,
                                                             nullptr, cdk.masses, cdk.natom,
                                                             nrg_scale_factor, inv_vel_scale);
  sc->contribute(StateVariable::KINETIC, acc, system_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tmass, typename Tcalc>
Tcalc computeTemperature(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                         const int* xvel_ovrf, const int* yvel_ovrf, const int* zvel_ovrf,
                         const Tmass* masses, const int natom, const int ndof,
                         const Tcalc nrg_scale_factor, const Tcalc inv_vel_scale) {

  // Begin by computing the kinetic energy
  const llint ke = evalKineticEnergy<Tcoord, Tmass, Tcalc>(xvel, yvel, zvel, xvel_ovrf, yvel_ovrf,
                                                           zvel_ovrf, masses, natom,
                                                           nrg_scale_factor, inv_vel_scale);
  const Tcalc tpre = static_cast<Tcalc>(ke) / (nrg_scale_factor * static_cast<Tcalc>(ndof));
  const Tcalc result = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) ?
                       2.0 * tpre / boltzmann_constant : 2.0f * tpre / boltzmann_constant_f;
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc computeTemperature(const PhaseSpace *ps, const AtomGraph *ag,
                         const ApplyConstraints use_cnst) {
  PhaseSpaceReader psr = ps->data();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  int ndof;
  switch (use_cnst) {
  case ApplyConstraints::YES:
    ndof = cdk.cnst_dof;
    break;
  case ApplyConstraints::NO:
    ndof = cdk.free_dof;
    break;
  }
  return computeTemperature<double, double, Tcalc>(psr.xvel, psr.yvel, psr.zvel, nullptr,
                                                   nullptr, nullptr, cdk.masses, cdk.natom, ndof,
                                                   pow(2.0, 32), 1.0);
}
  
//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc computeTemperature(const PhaseSpace &ps, const AtomGraph &ag,
                         const ApplyConstraints use_cnst) {
  return computeTemperature<Tcalc>(ps.getSelfPointer(), ag.getSelfPointer(), use_cnst);
}
                       
} // namespace mm
} // namespace stormm
