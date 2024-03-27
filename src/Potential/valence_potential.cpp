#include "copyright.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "Topology/atomgraph.h"
#include "valence_potential.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {
  return evaluateBondTerms<double, double, double>(vk, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                                   psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                                   psw.zfrc, ecard, eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {

  // Unpack the topology and initialize the energy result
  return evaluateBondTerms(ag.getDoublePrecisionValenceKit(), ps->data(), ecard, eval_force,
                           system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {

  // Unpack the topology and initialize the energy result
  return evaluateBondTerms(ag->getDoublePrecisionValenceKit(), ps->data(), ecard, eval_force,
                           system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                         ScoreCard *ecard, const int system_index) {
  return evaluateBondTerms<double, double, double>(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat,
                                                   cfr.invu, cfr.unit_cell, nullptr, nullptr,
                                                   nullptr, ecard, EvaluateForce::NO,
                                                   system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const ValenceKit<double> &vk, const CoordinateFrameWriter &cfw,
                         ScoreCard *ecard, const int system_index) {
  return evaluateBondTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const AtomGraph &ag, const CoordinateFrame &cf, ScoreCard *ecard,
                         const int system_index) {
  return evaluateBondTerms(ag.getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateBondTerms(const AtomGraph *ag, const CoordinateFrame &cf, ScoreCard *ecard,
                         const int system_index) {
  return evaluateBondTerms(ag->getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                          const EvaluateForce eval_force, const int system_index) {
  return evaluateAngleTerms<double, double, double>(vk, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                                    psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                                    psw.zfrc, ecard, eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                          const EvaluateForce eval_force, const int system_index) {
  return evaluateAngleTerms(ag.getDoublePrecisionValenceKit(), ps->data(), ecard,
                            eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                          const EvaluateForce eval_force, const int system_index) {
  return evaluateAngleTerms(ag->getDoublePrecisionValenceKit(), ps->data(), ecard,
                            eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                          ScoreCard *ecard, const int system_index) {
  return evaluateAngleTerms<double, double, double>(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat,
                                                    cfr.invu, cfr.unit_cell, nullptr, nullptr,
                                                    nullptr, ecard, EvaluateForce::NO,
                                                    system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                          ScoreCard *ecard, const int system_index) {
  return evaluateAngleTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const AtomGraph &ag, const CoordinateFrame &cf, ScoreCard *ecard,
                          const int system_index) {
  return evaluateAngleTerms(ag.getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateAngleTerms(const AtomGraph *ag, const CoordinateFrame &cf, ScoreCard *ecard,
                          const int system_index) {
  return evaluateAngleTerms(ag->getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                              const EvaluateForce eval_force, const int system_index) {
  return evaluateDihedralTerms<double, double, double>(vk, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                                       psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                                       psw.zfrc, ecard, eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                              const EvaluateForce eval_force, const int system_index) {
  return evaluateDihedralTerms(ag.getDoublePrecisionValenceKit(), ps->data(), ecard, eval_force,
                               system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                              const EvaluateForce eval_force, const int system_index) {
  return evaluateDihedralTerms(ag->getDoublePrecisionValenceKit(), ps->data(), ecard, eval_force,
                               system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                              ScoreCard *ecard, const int system_index) {
  return evaluateDihedralTerms<double, double, double>(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat,
                                                       cfr.invu, cfr.unit_cell, nullptr, nullptr,
                                                       nullptr, ecard, EvaluateForce::NO,
                                                       system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                              ScoreCard *ecard, const int system_index) {
  return evaluateDihedralTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const AtomGraph &ag, const CoordinateFrame &cf, ScoreCard *ecard,
                              const int system_index) {
  return evaluateDihedralTerms(ag.getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateDihedralTerms(const AtomGraph *ag, const CoordinateFrame &cf, ScoreCard *ecard,
                              const int system_index) {
  return evaluateDihedralTerms(ag->getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw,
                                ScoreCard *ecard, const EvaluateForce eval_force,
                                const int system_index) {
  return evaluateUreyBradleyTerms<double, double, double>(vk, psw.xcrd, psw.ycrd, psw.zcrd,
                                                          psw.umat, psw.invu, psw.unit_cell,
                                                          psw.xfrc, psw.yfrc, psw.zfrc, ecard,
                                                          eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                const EvaluateForce eval_force, const int system_index) {
  return evaluateUreyBradleyTerms(ag.getDoublePrecisionValenceKit(), ps->data(), ecard, eval_force,
                                  system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                                const EvaluateForce eval_force, const int system_index) {
  return evaluateUreyBradleyTerms(ag->getDoublePrecisionValenceKit(), ps->data(), ecard,
                                  eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                                ScoreCard *ecard, const int system_index) {
  return evaluateUreyBradleyTerms<double, double, double>(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd,
                                                          cfr.umat, cfr.invu, cfr.unit_cell,
                                                          nullptr, nullptr, nullptr, ecard,
                                                          EvaluateForce::NO, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                                ScoreCard *ecard, const int system_index) {
  return evaluateUreyBradleyTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const AtomGraph &ag, const CoordinateFrame &cf, ScoreCard *ecard,
                                const int system_index) {
  return evaluateUreyBradleyTerms(ag.getDoublePrecisionValenceKit(), cf.data(), ecard,
                                  system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateUreyBradleyTerms(const AtomGraph *ag, const CoordinateFrame &cf, ScoreCard *ecard,
                                const int system_index) {
  return evaluateUreyBradleyTerms(ag->getDoublePrecisionValenceKit(), cf.data(), ecard,
                                  system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw,
                                   ScoreCard *ecard, const EvaluateForce eval_force,
                                   const int system_index) {
  return evaluateCharmmImproperTerms<double, double, double>(vk, psw.xcrd, psw.ycrd, psw.zcrd,
                                                             psw.umat, psw.invu, psw.unit_cell,
                                                             psw.xfrc, psw.yfrc, psw.zfrc, ecard,
                                                             eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                   const EvaluateForce eval_force, const int system_index) {
  return evaluateCharmmImproperTerms(ag.getDoublePrecisionValenceKit(), ps->data(), ecard,
                                     eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                                   const EvaluateForce eval_force, const int system_index) {
  return evaluateCharmmImproperTerms(ag->getDoublePrecisionValenceKit(), ps->data(), ecard,
                                     eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                                   ScoreCard *ecard, const int system_index) {
  return evaluateCharmmImproperTerms<double, double, double>(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd,
                                                             cfr.umat, cfr.invu, cfr.unit_cell,
                                                             nullptr, nullptr, nullptr, ecard,
                                                             EvaluateForce::NO, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                                   ScoreCard *ecard, const int system_index) {
  return evaluateCharmmImproperTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const AtomGraph &ag, const CoordinateFrame &cf,
                                   ScoreCard *ecard, const int system_index) {
  return evaluateCharmmImproperTerms(ag.getDoublePrecisionValenceKit(), cf.data(), ecard,
                                     system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCharmmImproperTerms(const AtomGraph *ag, const CoordinateFrame &cf,
                                   ScoreCard *ecard, const int system_index) {
  return evaluateCharmmImproperTerms(ag->getDoublePrecisionValenceKit(), cf.data(), ecard,
                                     system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const ValenceKit<double> vk, PhaseSpaceWriter psw, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {
  return evaluateCmapTerms<double, double, double>(vk, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                                   psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                                   psw.zfrc, ecard, eval_force, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {
  return evaluateCmapTerms(ag.getDoublePrecisionValenceKit(), ps->data(), ecard, eval_force,
                           system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                         const EvaluateForce eval_force, const int system_index) {
  return evaluateCmapTerms(ag->getDoublePrecisionValenceKit(), ps->data(), ecard, eval_force,
                           system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const ValenceKit<double> vk, const CoordinateFrameReader cfr,
                         ScoreCard *ecard, const int system_index) {
  return evaluateCmapTerms<double, double, double>(vk, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat,
                                                   cfr.invu, cfr.unit_cell, nullptr, nullptr,
                                                   nullptr, ecard, EvaluateForce::NO,
                                                   system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const ValenceKit<double> vk, const CoordinateFrameWriter &cfw,
                         ScoreCard *ecard, const int system_index) {
  return evaluateCmapTerms(vk, CoordinateFrameReader(cfw), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const AtomGraph &ag, const CoordinateFrame &cf, ScoreCard *ecard,
                         const int system_index) {
  return evaluateCmapTerms(ag.getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double evaluateCmapTerms(const AtomGraph *ag, const CoordinateFrame &cf, ScoreCard *ecard,
                         const int system_index) {
  return evaluateCmapTerms(ag->getDoublePrecisionValenceKit(), cf.data(), ecard, system_index);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  PhaseSpaceWriter psw, ScoreCard *ecard,
                                  const EvaluateForce eval_elec_force,
                                  const EvaluateForce eval_vdw_force, const int system_index,
                                  const double clash_minimum_distance, const double clash_ratio) {
  return evaluateAttenuated14Terms<double, double, double>(vk, nbk, psw.xcrd, psw.ycrd, psw.zcrd,
                                                           psw.umat, psw.invu, psw.unit_cell,
                                                           psw.xfrc, psw.yfrc, psw.zfrc, ecard,
                                                           eval_elec_force, eval_vdw_force,
                                                           system_index, clash_minimum_distance,
                                                           clash_ratio);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const AtomGraph &ag, PhaseSpace *ps, ScoreCard *ecard,
                                  const EvaluateForce eval_elec_force,
                                  const EvaluateForce eval_vdw_force, const int system_index,
                                  const double clash_minimum_distance, const double clash_ratio) {
  return evaluateAttenuated14Terms(ag.getDoublePrecisionValenceKit(),
                                   ag.getDoublePrecisionNonbondedKit(), ps->data(), ecard,
                                   eval_elec_force, eval_vdw_force, system_index,
                                   clash_minimum_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const AtomGraph *ag, PhaseSpace *ps, ScoreCard *ecard,
                                  const EvaluateForce eval_elec_force,
                                  const EvaluateForce eval_vdw_force, const int system_index,
                                  const double clash_minimum_distance, const double clash_ratio) {
  return evaluateAttenuated14Terms(ag->getDoublePrecisionValenceKit(),
                                   ag->getDoublePrecisionNonbondedKit(), ps->data(), ecard,
                                   eval_elec_force, eval_vdw_force, system_index,
                                   clash_minimum_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  const CoordinateFrameReader cfr, ScoreCard *ecard,
                                  const int system_index, const double clash_minimum_distance,
                                  const double clash_ratio) {
  return evaluateAttenuated14Terms<double, double, double>(vk, nbk, cfr.xcrd, cfr.ycrd, cfr.zcrd,
                                                           cfr.umat, cfr.invu, cfr.unit_cell,
                                                           nullptr, nullptr, nullptr, ecard,
                                                           EvaluateForce::NO, EvaluateForce::NO,
                                                           system_index, clash_minimum_distance,
                                                           clash_ratio);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const ValenceKit<double> vk, const NonbondedKit<double> nbk,
                                  const CoordinateFrameWriter &cfw, ScoreCard *ecard,
                                  const int system_index, const double clash_minimum_distance,
                                  const double clash_ratio) {
  return evaluateAttenuated14Terms(vk, nbk, CoordinateFrameReader(cfw), ecard, system_index,
                                   clash_minimum_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const AtomGraph &ag, const CoordinateFrame &cf,
                                  ScoreCard *ecard, const int system_index,
                                  const double clash_minimum_distance, const double clash_ratio) {
  return evaluateAttenuated14Terms(ag.getDoublePrecisionValenceKit(),
                                   ag.getDoublePrecisionNonbondedKit(), cf.data(), ecard,
                                   system_index, clash_minimum_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
double2 evaluateAttenuated14Terms(const AtomGraph *ag, const CoordinateFrame &cf,
                                  ScoreCard *ecard, const int system_index,
                                  const double clash_minimum_distance, const double clash_ratio) {
  return evaluateAttenuated14Terms(ag->getDoublePrecisionValenceKit(),
                                   ag->getDoublePrecisionNonbondedKit(), cf.data(), ecard,
                                   system_index, clash_minimum_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
double evaluateRestraints(const RestraintKit<double, double2, double4> rar, PhaseSpaceWriter psw,
                          ScoreCard *ecard, const EvaluateForce eval_force, const int system_index,
                          const int step_number) {
  return evaluateRestraints<double, double,
                            double, double2, double4>(rar, psw.xcrd, psw.ycrd, psw.zcrd, psw.umat,
                                                      psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc,
                                                      psw.zfrc, ecard, eval_force, system_index,
                                                      step_number);
}

//-------------------------------------------------------------------------------------------------
double evaluateRestraints(const RestraintApparatus &ra, PhaseSpace *ps, ScoreCard *ecard,
                          const EvaluateForce eval_force, const int system_index,
                          const int step_number) {
  PhaseSpaceWriter psw = ps->data();
  return evaluateRestraints<double, double,
                            double, double2, double4>(ra.dpData(), psw.xcrd, psw.ycrd, psw.zcrd,
                                                      psw.umat, psw.invu, psw.unit_cell, psw.xfrc,\
                                                      psw.yfrc, psw.zfrc, ecard, eval_force,
                                                      system_index, step_number);
}

//-------------------------------------------------------------------------------------------------
double evaluateRestraints(const RestraintApparatus *ra, PhaseSpace *ps, ScoreCard *ecard,
                          const EvaluateForce eval_force, const int system_index,
                          const int step_number) {
  PhaseSpaceWriter psw = ps->data();
  return evaluateRestraints<double, double,
                            double, double2, double4>(ra->dpData(), psw.xcrd, psw.ycrd, psw.zcrd,
                                                      psw.umat, psw.invu, psw.unit_cell, psw.xfrc,
                                                      psw.yfrc, psw.zfrc, ecard, eval_force,
                                                      system_index, step_number);
}

//-------------------------------------------------------------------------------------------------
double evaluateRestraints(const RestraintKit<double, double2, double4> rar,
                          const CoordinateFrameReader cfr, ScoreCard *ecard,
                          const int system_index, const int step_number) {
  return evaluateRestraints<double, double,
                            double, double2, double4>(rar, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat,
                                                      cfr.invu, cfr.unit_cell, nullptr, nullptr,
                                                      nullptr, ecard, EvaluateForce::NO,
                                                      system_index, step_number);
}

//-------------------------------------------------------------------------------------------------
double evaluateRestraints(const RestraintKit<double, double2, double4> rar,
                          const CoordinateFrameWriter &cfw, ScoreCard *ecard,
                          const int system_index, const int step_number) {
  const CoordinateFrameReader cfr(cfw);
  return evaluateRestraints<double, double,
                            double, double2, double4>(rar, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat,
                                                      cfr.invu, cfr.unit_cell, nullptr, nullptr,
                                                      nullptr, ecard, EvaluateForce::NO,
                                                      system_index, step_number);
}

//-------------------------------------------------------------------------------------------------
double evaluateRestraints(const RestraintApparatus &ra, const CoordinateFrameReader cfr,
                          ScoreCard *ecard, const int system_index, const int step_number) {
  return evaluateRestraints<double, double,
                            double, double2, double4>(ra.dpData(), cfr.xcrd, cfr.ycrd, cfr.zcrd,
                                                      cfr.umat, cfr.invu, cfr.unit_cell, nullptr,
                                                      nullptr, nullptr, ecard, EvaluateForce::NO,
                                                      system_index, step_number);
}

//-------------------------------------------------------------------------------------------------
double evaluateRestraints(const RestraintApparatus *ra, const CoordinateFrameReader cfr,
                          ScoreCard *ecard, const int system_index, const int step_number) {
  return evaluateRestraints<double, double,
                            double, double2, double4>(ra->dpData(), cfr.xcrd, cfr.ycrd, cfr.zcrd,
                                                      cfr.umat, cfr.invu, cfr.unit_cell, nullptr,
                                                      nullptr, nullptr, ecard, EvaluateForce::NO,
                                                      system_index, step_number);
}

} // namespace energy
} // namespace stormm
