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
void evalValeRestMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                    const NonbondedKit<double> &nbk,
                    const RestraintKit<double, double2, double4> &rar,
                    const EvaluateForce eval_force, const int step, const int system_index,
                    const double clash_distance, const double clash_ratio) {
  evalValeMM<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc, sc, vk, nbk,
                                     eval_force, system_index, clash_distance,
                                     clash_ratio);
  evaluateRestraints<double, double,
                     double, double2, double4>(rar, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
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
                     double, double2, double4>(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                               psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                               psw.zfrc, sc, eval_force, system_index, step, 1.0,
                                               1.0);
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
                     double, double2, double4>(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                               psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                               psw.zfrc, sc, eval_force, system_index, step, 1.0,
                                               1.0);
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
                        const RestraintKit<double, double2, double4> &rar,
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
                     double, double2, double4>(rar, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
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
                     double, double2, double4>(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                               psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                               psw.zfrc, sc, eval_force, system_index, step, 1.0,
                                               1.0);
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
                     double, double2, double4>(ra->dpData(), psw.xcrd, psw.ycrd, psw.zcrd,
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
                        const RestraintKit<double, double2, double4> &rar,
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

} // namespace mm
} // namespace stormm
