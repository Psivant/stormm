// -*-c++-*-
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "MolecularMechanics/mm_controls.h"
#include "Synthesis/synthesis_abstracts.h"
#include "pmigrid.h"
#include "hpc_pme_potential.h"
#include "hpc_pme_potential.cuh"

namespace stormm {
namespace energy {

#ifdef STORMM_USE_CUDA
//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes queryPMEPairsKernelRequirements(const PrecisionModel coord_prec,
                                                          const PrecisionModel calc_prec,
                                                          const NeighborListKind neighbor_layout,
                                                          const EvaluateForce eval_frc,
                                                          const EvaluateEnergy eval_nrg,
                                                          const TinyBoxPresence has_tiny_box,
                                                          const ClashResponse clash_handling) {

  // As with other kernel querying functions, the kernel manager calling this function will have
  // specifications of the GPU in use.  It is the overall thread occupancy and multiplicity of each
  // kernel that this function must return.
  switch (coord_prec) {
  case PrecisionModel::DOUBLE:
    switch (calc_prec) {
    case PrecisionModel::DOUBLE:
      return queryDDPMEPairsKernelRequirements(neighbor_layout, eval_frc, eval_nrg, has_tiny_box,
                                               clash_handling);
    case PrecisionModel::SINGLE:
      return queryDFPMEPairsKernelRequirements(neighbor_layout, eval_frc, eval_nrg, has_tiny_box,
                                               clash_handling);
    }
    break;
  case PrecisionModel::SINGLE:
    switch (calc_prec) {
    case PrecisionModel::DOUBLE:
      return queryFDPMEPairsKernelRequirements(neighbor_layout, eval_frc, eval_nrg, has_tiny_box,
                                               clash_handling);
    case PrecisionModel::SINGLE:
      return queryFFPMEPairsKernelRequirements(neighbor_layout, eval_frc, eval_nrg, has_tiny_box,
                                               clash_handling);
    }
    break;
  }
    __builtin_unreachable();
}
#endif

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const PrecisionModel prec, const LocalExclusionMask &lem,
                           const PPITable &pairs_tbl,
                           CellGrid<double, llint, double, double4_16a> *cg, TileManager *tlmn,
                           ScoreCard *sc, MolecularMechanicsControls *mmctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const CoreKlManager &launcher, const double clash_distance,
                           const double clash_ratio) {

  // Obtain the correct kernel launch parameters
  const ClashResponse mitigate_clash = (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) ?
                                       ClashResponse::FORGIVE : ClashResponse::NONE;
  const AtomGraphSynthesis *poly_ag = cg->getTopologySynthesisPointer();
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  const LocalExclusionMaskReader lemr = lem.data(devc);
  TilePlan tlpn = tlmn->data(devc);
  ScoreCardWriter scw = sc->data(devc);
  CellGridWriter<double, llint, double, double4_16a> cgw = cg->data(devc);
  const TinyBoxPresence has_tiny_box = cg->getTinyBoxPresence();
  const int2 bt_tp = launcher.getPMEPairsKernelDims(PrecisionModel::DOUBLE, prec, 
                                                    NeighborListKind::MONO, has_tiny_box,
                                                    eval_frc, eval_nrg, mitigate_clash);
  const PsSynthesisBorders sysbrd = cg->getUnitCellTransforms(devc);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag->getDoublePrecisionNonbondedKit(devc);
      const PPIKit<double, double4_16a> nrg_tab = pairs_tbl.dpData();
      MMControlKit<double> ctrl = mmctrl->dpData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw, &tlpn, &scw, &ctrl, eval_frc, eval_nrg,
                       bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag->getSinglePrecisionNonbondedKit(devc);
      const PPIKit<float, float4> nrg_tab = pairs_tbl.spData();
      MMControlKit<float> ctrl = mmctrl->spData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw, &tlpn, &scw, &ctrl, eval_frc, eval_nrg,
                       bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const PrecisionModel prec, const LocalExclusionMask &lem,
                           const PPITable &pairs_tbl,
                           CellGrid<double, llint, double, double4_16a> *cg_qq,
                           CellGrid<double, llint, double, double4_16a> *cg_lj, TileManager *tlmn,
                           ScoreCard *sc, MolecularMechanicsControls *mmctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const CoreKlManager &launcher, const double clash_distance,
                           const double clash_ratio) {

  // Obtain the correct kernel launch parameters
  const ClashResponse mitigate_clash = (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) ?
                                       ClashResponse::FORGIVE : ClashResponse::NONE;
  const AtomGraphSynthesis *poly_ag = cg_qq->getTopologySynthesisPointer();
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  const LocalExclusionMaskReader lemr = lem.data(devc);
  TilePlan tlpn = tlmn->data(devc);
  ScoreCardWriter scw = sc->data(devc);
  CellGridWriter<double, llint, double, double4_16a> cgw_qq = cg_qq->data(devc);
  CellGridWriter<double, llint, double, double4_16a> cgw_lj = cg_lj->data(devc);
  const TinyBoxPresence has_tiny_box = (cg_qq->getTinyBoxPresence() == TinyBoxPresence::YES ||
                                        cg_lj->getTinyBoxPresence() == TinyBoxPresence::YES) ?
                                       TinyBoxPresence::YES : TinyBoxPresence::NO;
  const int2 bt_tp = launcher.getPMEPairsKernelDims(PrecisionModel::DOUBLE, prec, 
                                                    NeighborListKind::DUAL, has_tiny_box,
                                                    eval_frc, eval_nrg, mitigate_clash);
  const PsSynthesisBorders sysbrd = cg_qq->getUnitCellTransforms(devc);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag->getDoublePrecisionNonbondedKit(devc);
      const PPIKit<double, double4_16a> nrg_tab = pairs_tbl.dpData();
      MMControlKit<double> ctrl = mmctrl->dpData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl,
                       eval_frc, eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag->getSinglePrecisionNonbondedKit(devc);
      const PPIKit<float, float4> nrg_tab = pairs_tbl.spData();
      MMControlKit<float> ctrl = mmctrl->spData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl,
                       eval_frc, eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const PrecisionModel prec, const LocalExclusionMask &lem,
                           const PPITable &pairs_tbl, CellGrid<float, int, float, float4> *cg,
                           TileManager *tlmn, ScoreCard *sc, MolecularMechanicsControls *mmctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const CoreKlManager &launcher, const double clash_distance,
                           const double clash_ratio) {

  // Obtain the correct kernel launch parameters
  const ClashResponse mitigate_clash = (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) ?
                                       ClashResponse::FORGIVE : ClashResponse::NONE;
  const AtomGraphSynthesis *poly_ag = cg->getTopologySynthesisPointer();
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  const LocalExclusionMaskReader lemr = lem.data(devc);
  TilePlan tlpn = tlmn->data(devc);
  ScoreCardWriter scw = sc->data(devc);
  CellGridWriter<float, int, float, float4> cgw = cg->data(devc);
  const TinyBoxPresence has_tiny_box = cg->getTinyBoxPresence();
  const int2 bt_tp = launcher.getPMEPairsKernelDims(PrecisionModel::SINGLE, prec,
                                                    NeighborListKind::MONO,
                                                    cg->getTinyBoxPresence(), eval_frc, eval_nrg,
                                                    mitigate_clash);
  const PsSynthesisBorders sysbrd = cg->getUnitCellTransforms(devc);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag->getDoublePrecisionNonbondedKit(devc);
      const PPIKit<double, double4_16a> nrg_tab = pairs_tbl.dpData();
      MMControlKit<double> ctrl = mmctrl->dpData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw, &tlpn, &scw, &ctrl, eval_frc, eval_nrg,
                       bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag->getSinglePrecisionNonbondedKit(devc);
      const PPIKit<float, float4> nrg_tab = pairs_tbl.spData();
      MMControlKit<float> ctrl = mmctrl->spData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw, &tlpn, &scw, &ctrl, eval_frc, eval_nrg,
                       bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchPMEPairs(const PrecisionModel prec, const LocalExclusionMask &lem,
                           const PPITable &pairs_tbl, CellGrid<float, int, float, float4> *cg_qq,
                           CellGrid<float, int, float, float4> *cg_lj, TileManager *tlmn,
                           ScoreCard *sc, MolecularMechanicsControls *mmctrl,
                           const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                           const CoreKlManager &launcher, const double clash_distance,
                           const double clash_ratio) {

  // Obtain the correct kernel launch parameters
  const ClashResponse mitigate_clash = (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) ?
                                       ClashResponse::FORGIVE : ClashResponse::NONE;
  const AtomGraphSynthesis *poly_ag = cg_qq->getTopologySynthesisPointer();
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  const LocalExclusionMaskReader lemr = lem.data(devc);
  TilePlan tlpn = tlmn->data(devc);
  ScoreCardWriter scw = sc->data(devc);
  CellGridWriter<float, int, float, float4> cgw_qq = cg_qq->data(devc);
  CellGridWriter<float, int, float, float4> cgw_lj = cg_lj->data(devc);
  const TinyBoxPresence has_tiny_box = (cg_qq->getTinyBoxPresence() == TinyBoxPresence::YES ||
                                        cg_lj->getTinyBoxPresence() == TinyBoxPresence::YES) ?
                                       TinyBoxPresence::YES : TinyBoxPresence::NO;
  const int2 bt_tp = launcher.getPMEPairsKernelDims(PrecisionModel::SINGLE, prec,
                                                    NeighborListKind::DUAL, has_tiny_box,
                                                    eval_frc, eval_nrg, mitigate_clash);
  const PsSynthesisBorders sysbrd = cg_qq->getUnitCellTransforms(devc);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag->getDoublePrecisionNonbondedKit(devc);
      const PPIKit<double, double4_16a> nrg_tab = pairs_tbl.dpData();
      MMControlKit<double> ctrl = mmctrl->dpData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl,
                       eval_frc, eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag->getSinglePrecisionNonbondedKit(devc);
      const PPIKit<float, float4> nrg_tab = pairs_tbl.spData();
      MMControlKit<float> ctrl = mmctrl->spData(devc);
      switch (has_tiny_box) {
      case TinyBoxPresence::NO:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl, eval_frc,
                       eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      case TinyBoxPresence::YES:
        launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw_qq, &cgw_lj, &tlpn, &scw, &ctrl,
                       eval_frc, eval_nrg, bt_tp, clash_distance, clash_ratio);
        break;
      }
    }
    break;
  }
}

} // namespace energy
} // namespace stormm
