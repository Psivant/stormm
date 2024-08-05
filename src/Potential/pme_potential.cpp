#include "copyright.h"
#include "pme_potential.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
double2 evaluateParticleParticleEnergy(PhaseSpace *ps, const AtomGraph *ag,
                                       const LocalExclusionMask &lema, const PrecisionModel prec,
                                       const double elec_cutoff, const double vdw_cutoff,
                                       const double qqew_coeff, const double ljew_coeff,
                                       const VdwSumMethod vdw_sum, const EvaluateForce eval_frc,
                                       const NonbondedTheme theme) {
  PhaseSpaceWriter psw = ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    return evaluateParticleParticleEnergy<double>(&psw, ag->getDoublePrecisionNonbondedKit(),
                                                  lema.data(), elec_cutoff, vdw_cutoff, qqew_coeff,
                                                  ljew_coeff, vdw_sum, eval_frc, theme);
    break;
  case PrecisionModel::SINGLE:
    return evaluateParticleParticleEnergy<float>(&psw, ag->getSinglePrecisionNonbondedKit(),
                                                 lema.data(), elec_cutoff, vdw_cutoff, qqew_coeff,
                                                 ljew_coeff, vdw_sum, eval_frc, theme);
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double2 evaluateParticleParticleEnergy(PhaseSpace *ps, const AtomGraph &ag,
                                       const LocalExclusionMask &lema, const PrecisionModel prec,
                                       const double elec_cutoff, const double vdw_cutoff,
                                       const double qqew_coeff, const double ljew_coeff,
                                       const VdwSumMethod vdw_sum, const EvaluateForce eval_frc,
                                       const NonbondedTheme theme) {
  return evaluateParticleParticleEnergy(ps, ag.getSelfPointer(), lema, prec, elec_cutoff,
                                        vdw_cutoff, qqew_coeff, ljew_coeff, vdw_sum, eval_frc,
                                        theme);
}

} // namespace energy
} // namespace stormm
