#include "copyright.h"
#include "mm_evaluation.h"

namespace stormm {
namespace mm {

//-------------------------------------------------------------------------------------------------
void evalValeMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                const NonbondedKit<double> &nbk, const EvaluateForce eval_force,
                const int system_index, const double clash_distance, const double clash_ratio) {
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, vk, nbk,
                                     eval_force, system_index, 1.0, 1.0, clash_distance,
                                     clash_ratio);
}

//-------------------------------------------------------------------------------------------------
void evalValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                const EvaluateForce eval_force, const int system_index,
                const double clash_distance, const double clash_ratio) {
  PhaseSpaceWriter psw = ps->data();
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc,
                                     ag.getDoublePrecisionValenceKit(),
                                     ag.getDoublePrecisionNonbondedKit(), eval_force,
                                     system_index, 1.0, 1.0, clash_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
void evalValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                const EvaluateForce eval_force, const int system_index,
                const double clash_distance, const double clash_ratio) {
  PhaseSpaceWriter psw = ps->data();
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc,
                                     ag->getDoublePrecisionValenceKit(),
                                     ag->getDoublePrecisionNonbondedKit(), eval_force,
                                     system_index, 1.0, 1.0, clash_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
void evalValeMM(PhaseSpaceSynthesis* poly_ps, ScoreCard* sc, const AtomGraphSynthesis* poly_ag,
		const EvaluateForce eval_force, const PrecisionModel prec,
                const double clash_distance, const double clash_ratio) {
  PsSynthesisWriter poly_psw = poly_ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag->getDoublePrecisionValenceKit();
      const SyAtomUpdateKit<double,
                            double2,
                            double4_16a> poly_auk = poly_ag->getDoublePrecisionAtomUpdateKit();
      const SyRestraintKit<double, double2, double4_16a> empty_rk;
      evalValeMM<double, double2, double4_16a>(&poly_psw, sc, poly_vk, poly_auk, eval_force,
                                           VwuTask::ALL_TASKS, clash_distance, clash_ratio, 0,
                                           empty_rk);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag->getSinglePrecisionValenceKit();
      const SyAtomUpdateKit<float,
                            float2, float4> poly_auk = poly_ag->getSinglePrecisionAtomUpdateKit();
      const SyRestraintKit<float, float2, float4> empty_rk;
      evalValeMM<float, float2, float4>(&poly_psw, sc, poly_vk, poly_auk, eval_force,
                                        VwuTask::ALL_TASKS, clash_distance, clash_ratio, 0,
                                        empty_rk);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void evalValeRestMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                    const NonbondedKit<double> &nbk,
                    const RestraintKit<double, double2, double4_16a> &rar,
                    const EvaluateForce eval_force, const int step, const int system_index,
                    const double clash_distance, const double clash_ratio) {
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, vk, nbk,
                                     eval_force, system_index, clash_distance,
                                     clash_ratio);
  evaluateRestraints<double, double,
                     double, double2, double4_16a>(rar, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                               psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                               psw.zfrc, sc, eval_force, system_index, step, 1.0,
                                               1.0);
}

//-------------------------------------------------------------------------------------------------
void evalValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                    const RestraintApparatus &ra, const EvaluateForce eval_force,
                    const int system_index, const int step, const double clash_distance,
                    const double clash_ratio) {
  PhaseSpaceWriter psw = ps->data();
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc,
                                     ag.getDoublePrecisionValenceKit(),
                                     ag.getDoublePrecisionNonbondedKit(), eval_force,
                                     system_index, 1.0, 1.0, clash_distance, clash_ratio);
  evaluateRestraints<double, double,
                     double, double2, double4_16a>(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd,
                                                   psw.umat, psw.invu, psw.unit_cell, psw.xfrc,
                                                   psw.yfrc, psw.zfrc, sc, eval_force,
                                                   system_index, step, 1.0, 1.0);
}

//-------------------------------------------------------------------------------------------------
void evalValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                    const RestraintApparatus &ra, const EvaluateForce eval_force,
                    const int system_index, const int step, const double clash_distance,
                    const double clash_ratio) {
  PhaseSpaceWriter psw = ps->data();
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc,
                                     psw.zfrc, sc, ag->getDoublePrecisionValenceKit(),
                                     ag->getDoublePrecisionNonbondedKit(), eval_force,
                                     system_index, 1.0, 1.0, clash_distance, clash_ratio);
  evaluateRestraints<double, double,
                     double, double2, double4_16a>(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd,
                                                   psw.umat, psw.invu, psw.unit_cell, psw.xfrc,
                                                   psw.yfrc, psw.zfrc, sc, eval_force,
                                                   system_index, step, 1.0, 1.0);
}

//-------------------------------------------------------------------------------------------------
void evalValeRestMM(PhaseSpaceSynthesis* poly_ps, ScoreCard* sc, const AtomGraphSynthesis* poly_ag,
                    const int step_number, const EvaluateForce eval_force,
                    const PrecisionModel prec, const double clash_distance,
                    const double clash_ratio) {
  PsSynthesisWriter poly_psw = poly_ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag->getDoublePrecisionValenceKit();
      const SyAtomUpdateKit<double,
                            double2,
                            double4_16a> poly_auk = poly_ag->getDoublePrecisionAtomUpdateKit();
      const SyRestraintKit<double,
                           double2,
                           double4_16a> poly_rk = poly_ag->getDoublePrecisionRestraintKit();
      evalValeMM<double, double2, double4_16a>(&poly_psw, sc, poly_vk, poly_auk, eval_force,
                                           VwuTask::ALL_TASKS, clash_distance, clash_ratio,
                                           step_number, poly_rk);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag->getSinglePrecisionValenceKit();
      const SyAtomUpdateKit<float,
                            float2, float4> poly_auk = poly_ag->getSinglePrecisionAtomUpdateKit();
      const SyRestraintKit<float,
                           float2, float4> poly_rk = poly_ag->getSinglePrecisionRestraintKit();
      evalValeMM<float, float2, float4>(&poly_psw, sc, poly_vk, poly_auk, eval_force,
                                        VwuTask::ALL_TASKS, clash_distance, clash_ratio,
                                        step_number, poly_rk);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                    const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                    const EvaluateForce eval_force, const int system_index,
                    const double clash_distance, const double clash_ratio) {
  evaluateNonbondedEnergy<double, double, double>(nbk, ser, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                                  psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                                  psw.zfrc, sc, eval_force, eval_force,
                                                  system_index, 1.0, 1.0, clash_distance,
                                                  clash_ratio);
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, vk, nbk,
                                     eval_force, system_index, 1.0, 1.0, clash_distance,
                                     clash_ratio);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                    const StaticExclusionMask &se, const EvaluateForce eval_force,
                    const int system_index, const double clash_distance,
                    const double clash_ratio) {
  PhaseSpaceWriter psw = ps->data();
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  evaluateNonbondedEnergy<double, double, double>(nbk, se.data(), psw.xcrd, psw.ycrd, psw.zcrd,
                                                  psw.umat, psw.invu, psw.unit_cell, psw.xfrc,
                                                  psw.yfrc, psw.zfrc, sc, eval_force, eval_force,
                                                  system_index, 1.0, 1.0, clash_distance,
                                                  clash_ratio);
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc,
                                     ag.getDoublePrecisionValenceKit(), nbk, eval_force,
                                     system_index, 1.0, 1.0, clash_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                    const StaticExclusionMask &se, const EvaluateForce eval_force,
                    const int system_index, const double clash_distance,
                    const double clash_ratio) {
  PhaseSpaceWriter psw = ps->data();
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  evaluateNonbondedEnergy<double, double, double>(nbk, se.data(), psw.xcrd, psw.ycrd, psw.zcrd,
                                                  psw.umat, psw.invu, psw.unit_cell, psw.xfrc,
                                                  psw.yfrc, psw.zfrc, sc, eval_force, eval_force,
                                                  system_index, 1.0, 1.0, clash_distance,
                                                  clash_ratio);
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc,
                                     ag->getDoublePrecisionValenceKit(), nbk, eval_force,
                                     system_index, 1.0, 1.0, clash_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeRestMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                        const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                        const RestraintKit<double, double2, double4_16a> &rar,
                        const EvaluateForce eval_force, const int system_index, const int step,
                        const double clash_distance, const double clash_ratio) {
  evaluateNonbondedEnergy<double, double, double>(nbk, ser, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                                  psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                                  psw.zfrc, sc, eval_force, eval_force,
                                                  system_index, 1.0, 1.0, clash_distance,
                                                  clash_ratio);
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, vk, nbk,
                                     eval_force, system_index, 1.0, 1.0, clash_distance,
                                     clash_ratio);
  evaluateRestraints<double, double,
                     double, double2, double4_16a>(rar, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                               psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                               psw.zfrc, sc, eval_force, system_index, step, 1.0,
                                               1.0);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                        const StaticExclusionMask &se, const RestraintApparatus &ra,
                        const EvaluateForce eval_force, const int system_index, const int step,
                        const double clash_distance, const double clash_ratio) {
  PhaseSpaceWriter psw = ps->data();
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  evaluateNonbondedEnergy<double, double, double>(nbk, se.data(), psw.xcrd, psw.ycrd, psw.zcrd,
                                                  psw.umat, psw.invu, psw.unit_cell, psw.xfrc,
                                                  psw.yfrc, psw.zfrc, sc, eval_force, eval_force,
                                                  system_index, 1.0, 1.0, clash_distance,
                                                  clash_ratio);
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc,
                                     ag.getDoublePrecisionValenceKit(), nbk, eval_force,
                                     system_index, 1.0, 1.0, clash_distance, clash_ratio);
  evaluateRestraints<double, double,
                     double, double2, double4_16a>(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd,
                                                   psw.umat, psw.invu, psw.unit_cell, psw.xfrc,
                                                   psw.yfrc, psw.zfrc, sc, eval_force,
                                                   system_index, step, 1.0, 1.0);
}

//-------------------------------------------------------------------------------------------------
void evalNonbValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                        const StaticExclusionMask &se, const RestraintApparatus *ra,
                        const EvaluateForce eval_force, const int system_index, const int step,
                        const double clash_distance, const double clash_ratio) {
  PhaseSpaceWriter psw = ps->data();
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  evaluateNonbondedEnergy<double, double, double>(nbk, se.data(), psw.xcrd, psw.ycrd, psw.zcrd,
                                                  psw.umat, psw.invu, psw.unit_cell, psw.xfrc,
                                                  psw.yfrc, psw.zfrc, sc, eval_force, eval_force,
                                                  system_index, 1.0, 1.0, clash_distance,
                                                  clash_ratio);
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc,
                                     ag->getDoublePrecisionValenceKit(), nbk, eval_force,
                                     system_index, 1.0, 1.0, clash_distance, clash_ratio);
  evaluateRestraints<double, double,
                     double, double2, double4_16a>(ra->dpData(), psw.xcrd, psw.ycrd, psw.zcrd,
                                               psw.umat, psw.invu, psw.unit_cell, psw.xfrc,
                                               psw.yfrc, psw.zfrc, sc, eval_force, system_index,
                                               step, 1.0, 1.0);
}

//-------------------------------------------------------------------------------------------------
void evalRestrainedMMGB(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                        const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                        const ImplicitSolventKit<double> &isk,
                        const NeckGeneralizedBornKit<double> &neck_gbk,
                        double* effective_gb_radii, double *psi, double *sumdeijda,
                        const RestraintKit<double, double2, double4_16a> &rar,
                        const EvaluateForce eval_force, const int system_index, const int step,
                        const double clash_distance, const double clash_ratio) {
  evaluateGeneralizedBornEnergy(nbk, ser, isk, neck_gbk, psw, sc, eval_force, system_index);
  evalNonbValeRestMM(psw, sc, vk, nbk, ser, rar, eval_force, system_index, step, clash_distance,
                     clash_ratio);
}

//-------------------------------------------------------------------------------------------------
void evalRestrainedMMGB(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                        const NeckGeneralizedBornTable &neck_gbtab,
			const StaticExclusionMask &se, const RestraintApparatus &ra,
                        const EvaluateForce eval_force, const int system_index, const int step,
                        const double clash_distance, const double clash_ratio) {
  const ImplicitSolventKit<double> isk = ag.getDoublePrecisionImplicitSolventKit();
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const StaticExclusionMaskReader ser = se.data();
  evaluateGeneralizedBornEnergy(nbk, ser, isk, neck_gbtab.dpData(), ps->data(), sc, eval_force,
                                system_index);
  evalNonbValeRestMM(ps->data(), sc, vk, nbk, ser, ra.dpData(), eval_force, system_index, step,
                     clash_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
void evalRestrainedMMGB(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                        const NeckGeneralizedBornTable &neck_gbtab,
			const StaticExclusionMask &se, const RestraintApparatus *ra,
                        const EvaluateForce eval_force, const int system_index, const int step,
                        const double clash_distance, const double clash_ratio) {
  const ImplicitSolventKit<double> isk = ag->getDoublePrecisionImplicitSolventKit();
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const ValenceKit<double> vk = ag->getDoublePrecisionValenceKit();
  const StaticExclusionMaskReader ser = se.data();
  evaluateGeneralizedBornEnergy(nbk, ser, isk, neck_gbtab.dpData(), ps->data(), sc, eval_force,
                                system_index);
  evalNonbValeRestMM(ps->data(), sc, vk, nbk, ser, ra->dpData(), eval_force, system_index, step,
                     clash_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
void evalVwuInitEnergy(ScoreCard *ecard, const VwuTask activity, const int sysid) {
  switch (activity) {
  case VwuTask::BOND:
    ecard->initialize(StateVariable::BOND, sysid);
    break;
  case VwuTask::ANGL:
    ecard->initialize(StateVariable::ANGLE, sysid);
    break;
  case VwuTask::DIHE:
    ecard->initialize(StateVariable::PROPER_DIHEDRAL, sysid);
    ecard->initialize(StateVariable::IMPROPER_DIHEDRAL, sysid);
    break;
  case VwuTask::UBRD:
    ecard->initialize(StateVariable::UREY_BRADLEY, sysid);
    break;
  case VwuTask::CIMP:
    ecard->initialize(StateVariable::CHARMM_IMPROPER, sysid);
    break;
  case VwuTask::CMAP:
    ecard->initialize(StateVariable::CMAP, sysid);
    break;
  case VwuTask::INFR14:
    ecard->initialize(StateVariable::ELEC_ONE_FOUR, sysid);
    ecard->initialize(StateVariable::VDW_ONE_FOUR, sysid);
    break;
  case VwuTask::RPOSN:
  case VwuTask::RBOND:
  case VwuTask::RANGL:
  case VwuTask::RDIHE:
    ecard->initialize(StateVariable::RESTRAINT, sysid);
    break;
  case VwuTask::SETTLE:
  case VwuTask::CGROUP:
  case VwuTask::VSITE:
  case VwuTask::CDHE:
  case VwuTask::CBND:
    break;
  case VwuTask::ALL_TASKS:
    ecard->initialize({ StateVariable::BOND, StateVariable::ANGLE, StateVariable::PROPER_DIHEDRAL,
                        StateVariable::IMPROPER_DIHEDRAL, StateVariable::UREY_BRADLEY,
                        StateVariable::CHARMM_IMPROPER, StateVariable::CMAP,
                        StateVariable::ELEC_ONE_FOUR, StateVariable::VDW_ONE_FOUR,
                        StateVariable::RESTRAINT }, sysid);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void commitVwuEnergies(const llint bond_acc, const llint angl_acc, const llint dihe_acc,
                       const llint impr_acc, const llint ubrd_acc, const llint cimp_acc,
                       const llint cmap_acc, const llint qq14_acc, const llint lj14_acc,
                       const llint rest_acc, const int sysid, const VwuTask activity,
                       ScoreCard *ecard) {
  switch (activity) {
  case VwuTask::BOND:
    ecard->add(StateVariable::BOND, bond_acc, sysid);
    break;
  case VwuTask::ANGL:
    ecard->add(StateVariable::ANGLE, angl_acc, sysid);
    break;
  case VwuTask::DIHE:
    ecard->add(StateVariable::PROPER_DIHEDRAL, dihe_acc, sysid);
    ecard->add(StateVariable::IMPROPER_DIHEDRAL, impr_acc, sysid);
    break;
  case VwuTask::UBRD:
    ecard->add(StateVariable::UREY_BRADLEY, ubrd_acc, sysid);
    break;
  case VwuTask::CIMP:
    ecard->add(StateVariable::CHARMM_IMPROPER, cimp_acc, sysid);
    break;
  case VwuTask::CMAP:
    ecard->add(StateVariable::CMAP, cmap_acc, sysid);
    break;
  case VwuTask::INFR14:
    ecard->add(StateVariable::ELEC_ONE_FOUR, qq14_acc, sysid);
    ecard->add(StateVariable::VDW_ONE_FOUR, lj14_acc, sysid);
    break;
  case VwuTask::RPOSN:
  case VwuTask::RBOND:
  case VwuTask::RANGL:
  case VwuTask::RDIHE:
    ecard->add(StateVariable::RESTRAINT, rest_acc, sysid);
    break;
  case VwuTask::SETTLE:
  case VwuTask::CGROUP:
  case VwuTask::VSITE:
  case VwuTask::CDHE:
  case VwuTask::CBND:
    break;
  case VwuTask::ALL_TASKS:
    ecard->add(StateVariable::BOND, bond_acc, sysid);
    ecard->add(StateVariable::ANGLE, angl_acc, sysid);
    ecard->add(StateVariable::PROPER_DIHEDRAL, dihe_acc, sysid);
    ecard->add(StateVariable::IMPROPER_DIHEDRAL, impr_acc, sysid);
    ecard->add(StateVariable::UREY_BRADLEY, ubrd_acc, sysid);
    ecard->add(StateVariable::CHARMM_IMPROPER, cimp_acc, sysid);
    ecard->add(StateVariable::CMAP, cimp_acc, sysid);
    ecard->add(StateVariable::ELEC_ONE_FOUR, qq14_acc, sysid);
    ecard->add(StateVariable::VDW_ONE_FOUR, lj14_acc, sysid);
    ecard->add(StateVariable::RESTRAINT, cimp_acc, sysid);
    break;
  }
}

} // namespace mm
} // namespace stormm
