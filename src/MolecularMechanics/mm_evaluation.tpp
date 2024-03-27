// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace mm {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
void evalValeMM(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                const double* invu, const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                Tforce* zfrc, ScoreCard *sc, const ValenceKit<Tcalc> &vk,
                const NonbondedKit<Tcalc> &nbk, const EvaluateForce eval_force,
                const int system_index, const Tcalc inv_gpos_factor, const Tcalc force_factor,
                const Tcalc clash_distance, const Tcalc clash_ratio) {
  evaluateBondTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc,
                                           zfrc, sc, eval_force, system_index, inv_gpos_factor,
                                           force_factor);
  evaluateAngleTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                            yfrc, zfrc, sc, eval_force, system_index,
                                            inv_gpos_factor, force_factor);
  evaluateDihedralTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                               yfrc, zfrc, sc, eval_force, system_index,
                                               inv_gpos_factor, force_factor);
  evaluateUreyBradleyTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell,
                                                  xfrc, yfrc, zfrc, sc, eval_force, system_index,
                                                  inv_gpos_factor, force_factor);
  evaluateCharmmImproperTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell,
                                                     xfrc, yfrc, zfrc, sc, eval_force,
                                                     system_index, inv_gpos_factor, force_factor);
  evaluateCmapTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc,
                                           zfrc, sc, eval_force, system_index, inv_gpos_factor,
                                           force_factor);
  evaluateAttenuated14Terms<Tcoord, Tforce, Tcalc>(vk, nbk, xcrd, ycrd, zcrd, umat, invu,
                                                   unit_cell, xfrc, yfrc, zfrc, sc, eval_force,
                                                   eval_force, system_index, inv_gpos_factor,
                                                   force_factor, clash_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
void evalValeRestMM(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                    const double* invu, const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                    Tforce* zfrc, ScoreCard *sc, const ValenceKit<Tcalc> &vk,
                    const NonbondedKit<Tcalc> &nbk, const RestraintKit<Tcalc, Tcalc2, Tcalc4> &rar,
                    const EvaluateForce eval_force, const int step, const int system_index,
                    const Tcalc inv_gpos_factor, const Tcalc force_factor,
                    const Tcalc clash_distance, const Tcalc clash_ratio) {
  evalValeMM<Tcoord, Tforce, Tcalc>(xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                                    vk, nbk, eval_force, system_index, inv_gpos_factor,
                                    force_factor, clash_distance, clash_ratio);
  evaluateRestraints<Tcoord, Tforce,
                     Tcalc, Tcalc2, Tcalc4>(rar, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                            yfrc, zfrc, sc, eval_force, system_index, step,
                                            inv_gpos_factor, force_factor);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
void evalNonbValeMM(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                    const double* invu, const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                    Tforce* zfrc, ScoreCard *sc, const ValenceKit<Tcalc> &vk,
                    const NonbondedKit<Tcalc> &nbk, const StaticExclusionMaskReader &ser,
                    const EvaluateForce eval_force, const int system_index,
                    const Tcalc inv_gpos_factor, const Tcalc force_factor,
                    const Tcalc clash_distance, const Tcalc clash_ratio) {
  evaluateNonbondedEnergy<Tcoord, Tforce, Tcalc>(nbk, ser, xcrd, ycrd, zcrd, umat, invu, unit_cell,
                                                 xfrc, yfrc, zfrc, sc, eval_force, eval_force,
                                                 system_index, inv_gpos_factor, force_factor,
                                                 clash_distance, clash_ratio);
  evalValeMM<Tcoord, Tforce, Tcalc>(xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                                    vk, nbk, eval_force, system_index, inv_gpos_factor,
                                    force_factor, clash_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
void evalNonbValeRestMM(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                        const double* umat, const double* invu, const UnitCellType unit_cell,
                        Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, ScoreCard *sc,
                        const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                        const StaticExclusionMaskReader &ser,
                        const RestraintKit<Tcalc, Tcalc2, Tcalc4> &rar,
                        const EvaluateForce eval_force, const int system_index, const int step,
                        const Tcalc inv_gpos_factor, const Tcalc force_factor,
                        const Tcalc clash_distance, const Tcalc clash_ratio) {
  evaluateNonbondedEnergy<Tcoord, Tforce, Tcalc>(nbk, ser, xcrd, ycrd, zcrd, umat, invu, unit_cell,
                                                 xfrc, yfrc, zfrc, sc, eval_force, eval_force,
                                                 system_index, inv_gpos_factor, force_factor,
                                                 clash_distance, clash_ratio);
  evalValeMM<Tcoord, Tforce, Tcalc>(xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                                    vk, nbk, eval_force, system_index, inv_gpos_factor,
                                    force_factor, clash_distance, clash_ratio);
  evaluateRestraints<Tcoord, Tforce,
                     Tcalc, Tcalc2, Tcalc4>(rar, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                            yfrc, zfrc, sc, eval_force, system_index, step,
                                            inv_gpos_factor, force_factor);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
void evalRestrainedMMGB(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                        const double* umat, const double* invu, const UnitCellType unit_cell,
                        Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, ScoreCard *sc,
                        const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                        const StaticExclusionMaskReader &ser, const ImplicitSolventKit<Tcalc> &isk,
                        const NeckGeneralizedBornKit<Tcalc> &neck_gbk,
                        Tforce* effective_gb_radii, Tforce *psi, Tforce *sumdeijda,
                        const RestraintKit<Tcalc, Tcalc2, Tcalc4> &rar,
                        const EvaluateForce eval_force, const int system_index, const int step,
                        const Tcalc inv_gpos_factor, const Tcalc force_factor,
                        const Tcalc clash_distance, const Tcalc clash_ratio) {
  evalNonbValeRestMM(xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc, vk, nbk, ser,
                     rar, eval_force, system_index, step, inv_gpos_factor, force_factor,
                     clash_distance, clash_ratio);
  evaluateGeneralizedBornEnergy<Tcoord, Tforce, Tcalc>(nbk, ser, isk, neck_gbk, xcrd, ycrd, zcrd,
                                                       xfrc, yfrc, zfrc, effective_gb_radii, psi,
                                                       sumdeijda, sc, eval_force, system_index,
                                                       inv_gpos_factor, force_factor);
}

} // namespace mm
} // namespace stormm
