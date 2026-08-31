// -*-c++-*-
#include "copyright.h"
#include "eval_synthesis.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
int hostGetTileSideAtomCount(const std::vector<int> &nbwu_map, const int pos) {
  const int key_idx  = pos / 4;
  const int key_slot = pos - (key_idx * 4);
  return ((nbwu_map[small_block_max_imports + 1 + key_idx] >> (8 * key_slot)) & 0xff);
}

//-------------------------------------------------------------------------------------------------
void evalSyNonbondedEnergy(const AtomGraphSynthesis &poly_ag,
                           const StaticExclusionMaskSynthesis &poly_se,
                           PhaseSpaceSynthesis *poly_ps, ImplicitSolventWorkspace *isw,
                           ScoreCard *ecard,
                           const NonbondedTask task, const PrecisionModel prec,
                           const EvaluateForce eval_elec_force, const EvaluateForce eval_vdw_force,
                           const double clash_minimum_distance, const double clash_ratio) {
  PsSynthesisWriter psyw = poly_ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    if (isw == nullptr) {
      evalSyNonbondedTileGroups<double>(poly_ag.getDoublePrecisionNonbondedKit(), poly_se.data(),
                                        &psyw, ecard, task, eval_elec_force, eval_vdw_force,
                                        clash_minimum_distance, clash_ratio, nullptr);
    }
    else {
      ISWorkspaceKit<double> iswk = isw->dpData();
      evalSyNonbondedTileGroups<double>(poly_ag.getDoublePrecisionNonbondedKit(), poly_se.data(),
                                        &psyw, ecard, task, eval_elec_force, eval_vdw_force,
                                        clash_minimum_distance, clash_ratio, &iswk);
    }
    break;
  case PrecisionModel::SINGLE:
    if (isw == nullptr) {
      evalSyNonbondedTileGroups<float>(poly_ag.getSinglePrecisionNonbondedKit(), poly_se.data(),
                                       &psyw, ecard, task, eval_elec_force, eval_vdw_force,
                                       clash_minimum_distance, clash_ratio, nullptr);
    }
    else {
      ISWorkspaceKit<float> iswk = isw->spData();
      evalSyNonbondedTileGroups<float>(poly_ag.getSinglePrecisionNonbondedKit(), poly_se.data(),
                                       &psyw, ecard, task, eval_elec_force, eval_vdw_force,
                                       clash_minimum_distance, clash_ratio, &iswk);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void evalSyNonbondedEnergy(const AtomGraphSynthesis &poly_ag,
                           const StaticExclusionMaskSynthesis &poly_se,
                           PhaseSpaceSynthesis *poly_ps, ScoreCard *ecard,
                           const NonbondedTask task, const PrecisionModel prec,
                           const EvaluateForce eval_elec_force, const EvaluateForce eval_vdw_force,
                           const double clash_minimum_distance, const double clash_ratio) {
  switch (poly_ag.getImplicitSolventModel()) {
  case ImplicitSolventModel::NONE:
    evalSyNonbondedEnergy(poly_ag, poly_se, poly_ps, nullptr, ecard, task, prec, eval_elec_force,
                          eval_vdw_force, clash_minimum_distance, clash_ratio);
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    {
      ImplicitSolventWorkspace isw(poly_ag.getSystemAtomOffsets(), poly_ag.getSystemAtomCounts(),
                                   prec);
      evalSyNonbondedEnergy(poly_ag, poly_se, poly_ps, &isw, ecard, task, prec, eval_elec_force,
                            eval_vdw_force, clash_minimum_distance, clash_ratio);
    }
    break;
  }
}

} // namespace energy
} // namespace stormm
