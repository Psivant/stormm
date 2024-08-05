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
llint evalRescaledKineticEnergy(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                                const int* xvel_ovrf, const int* yvel_ovrf, const int* zvel_ovrf,
                                const Tmass* masses, const int natom, const Tcalc nrg_scale_factor,
                                const Tcalc inv_vel_scale) {
  llint result = 0LL;
  const Tcalc my_gafs_to_kcal = gafs_to_kcal;
  if (isSignedIntegralScalarType<Tcoord>()){
    if (xvel_ovrf == nullptr && yvel_ovrf == nullptr && zvel_ovrf == nullptr) {
      for (int i = 0; i < natom; i++){
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
template <typename Tcoord, typename Tmass, typename Tcalc>
llint evalMomenta(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                  const int* xvel_ovrf, const int* yvel_ovrf, const int* zvel_ovrf,
                  const Tmass* masses, const int natom, const Tcalc inv_vel_scale) {
  llint result = 0LL;
  if (isSignedIntegralScalarType<Tcoord>()){
    if (xvel_ovrf == nullptr && yvel_ovrf == nullptr && zvel_ovrf == nullptr) {
      for (int i = 0; i < natom; i++){
        const Tcalc mss_i = static_cast<Tcalc>(masses[i]) * static_cast<Tcalc>(0.5);
        const Tcalc vx = static_cast<Tcalc>(xvel[i]) * inv_vel_scale;
        const Tcalc vy = static_cast<Tcalc>(yvel[i]) * inv_vel_scale;
        const Tcalc vz = static_cast<Tcalc>(zvel[i]) * inv_vel_scale;
        result += mss_i * ((vx) + (vy) + (vz));
      }   
    }
    else { 
      for (int i = 0; i < natom; i++) {
        const Tcalc mss_i = static_cast<Tcalc>(masses[i]) * static_cast<Tcalc>(0.5);
        const Tcalc vx = hostInt95ToDouble(xvel[i], xvel_ovrf[i]) * inv_vel_scale;
        const Tcalc vy = hostInt95ToDouble(yvel[i], yvel_ovrf[i]) * inv_vel_scale;
        const Tcalc vz = hostInt95ToDouble(zvel[i], zvel_ovrf[i]) * inv_vel_scale;
        result += mss_i * ((vx) + (vy) + (vz));
      }
    }
  }
  else {
    for (int i = 0; i < natom; i++) {
      const Tcalc mss_i = static_cast<Tcalc>(masses[i]) * static_cast<Tcalc>(0.5);
      const Tcalc contrib_i = mss_i * ((xvel[i]) + (yvel[i]) + (zvel[i]));
      result += llround(contrib_i);
    }
  }
  return result;
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
//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
std::vector<llint> initiateMomentaRescale(const PhaseSpaceSynthesis *PsSynthesis, 
                                          const AtomGraphSynthesis *AgSynthesis, 
                                          const ScoreCard *sc, std::vector<int> swap_indices, 
                                          std::vector<double> temp) {
  std::vector<llint> rescaledMomenta;
  const PsSynthesisReader pssr = PsSynthesis->data();
  const std::vector<AtomGraph*> topologies = AgSynthesis->getSystemTopologyPointer();
  const int system_count = pssr.system_count;
  const int* atom_starts = pssr.atom_starts;
  for(int i = 0; i < system_count; i++){
    const llint* xvel = &pssr.xvel[atom_starts[i]];
    const llint* yvel = &pssr.yvel[atom_starts[i]];
    const llint* zvel = &pssr.zvel[atom_starts[i]];
    const int* xvel_ovrf = nullptr;
    const int* yvel_ovrf = nullptr;
    const int* zvel_ovrf = nullptr;
    if (pssr.vel_bits >= velocity_scale_nonoverflow_bits){
    xvel_ovrf = &pssr.xvel_ovrf[atom_starts[i]];
    yvel_ovrf = &pssr.yvel_ovrf[atom_starts[i]];
    zvel_ovrf = &pssr.zvel_ovrf[atom_starts[i]];
    }
    const AtomGraph* top = topologies[pssr.unique_ag_idx[atom_starts[i]]];
    const ChemicalDetailsKit cdk = top->getChemicalDetailsKit();
    const int natom = pssr.atom_counts[atom_starts[i]];
    const int inv_vel_scale = pssr.inv_vel_scale;
    rescaledMomenta.push_back(evalMomenta(xvel, yvel, zvel, xvel_ovrf, yvel_ovrf, 
                                                zvel_ovrf, cdk.masses, natom, inv_vel_scale));
  }
  
  for(int i = 0; i < static_cast<int>(swap_indices.size()); i++) {
    const llint p1 = rescaledMomenta[swap_indices[i]];
    const llint p2 = rescaledMomenta[swap_indices[i + 1]];
    const double t1 = temp[swap_indices[i]];
    const double t2 = temp[swap_indices[i + 1]];
    rescaledMomenta[swap_indices[i]] = (std::sqrt(t2/t1) * p1);
    rescaledMomenta[swap_indices[i + 1]] = (std::sqrt(t1/t2) * p2);
  }
  return rescaledMomenta;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
std::vector<llint> initiateKineticEnergyRescale(const PhaseSpaceSynthesis *PsSynthesis, 
                                                const AtomGraphSynthesis *AgSynthesis, 
                                                const ScoreCard *sc) {
  std::vector<llint> results;
  const PsSynthesisReader pssr = PsSynthesis->data();
  const std::vector<AtomGraph*> topologies = AgSynthesis->getSystemTopologyPointer();
  const int system_count = pssr.system_count;
  const int* atom_starts = pssr.atom_starts;
  for(int i = 0; i < system_count; i++){
    const llint* xvel = &pssr.xvel[atom_starts[i]];
    const llint* yvel = &pssr.yvel[atom_starts[i]];
    const llint* zvel = &pssr.zvel[atom_starts[i]];
    const int* xvel_ovrf = nullptr;
    const int* yvel_ovrf = nullptr;
    const int* zvel_ovrf = nullptr;
    if (pssr.vel_bits >= velocity_scale_nonoverflow_bits){
    xvel_ovrf = &pssr.xvel_ovrf[atom_starts[i]];
    yvel_ovrf = &pssr.yvel_ovrf[atom_starts[i]];
    zvel_ovrf = &pssr.zvel_ovrf[atom_starts[i]];
    }
    const AtomGraph* top = topologies[pssr.unique_ag_idx[atom_starts[i]]];
    const ChemicalDetailsKit cdk = top->getChemicalDetailsKit();
    const int natom = pssr.atom_counts[atom_starts[i]];
    const double inv_vel_scale = pssr.inv_vel_scale;
    // Need to discuss this nrg_scale_factor (and how it is taken from scorecard)
    // This implementation is copied from line 167 of the same file
    llint result = evalRescaledKineticEnergy(xvel, yvel, zvel, xvel_ovrf, yvel_ovrf, zvel_ovrf,
                                            cdk.masses, natom, pow(2.0, 32), inv_vel_scale); 
    results.push_back(result);
  }
  return results;
}
                      
} // namespace mm
} // namespace stormm
