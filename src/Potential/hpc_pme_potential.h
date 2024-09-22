// -*-c++-*-
#ifndef STORMM_HPC_PME_POTENTIAL_H
#define STORMM_HPC_PME_POTENTIAL_H

#ifdef STORMM_USE_CUDA
#  include <cuda_runtime.h>
#endif
#include "copyright.h"
#include "Accelerator/core_kernel_manager.h"
#include "Constants/behavior.h"
#include "MolecularMechanics/mm_controls.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "energy_enumerators.h"
#include "local_exclusionmask.h"
#include "pme_potential.h"
#include "ppitable.h"
#include "scorecard.h"
#include "tile_manager.h"

namespace stormm {
namespace energy {

using card::CoreKlManager;
using constants::PrecisionModel;
using mm::MMControlKit;
using mm::MolecularMechanicsControls;
using synthesis::SyNonbondedKit;
using synthesis::PsSynthesisBorders;

#ifdef STORMM_USE_CUDA
/// \brief Return critical attributes of a selected particle-particle pair interactions kernel
///        based on selected features.  The kernel must serve tower-plate 
///
/// \param coord_prec
/// \param calc_prec
/// \param eval_frc
/// \param eval_nrg
/// \param neighbor_list
/// \param has_tiny_box        
/// \param clash_handling
cudaFuncAttributes queryPMEPairsKernelRequirements(PrecisionModel coord_prec,
                                                   PrecisionModel calc_prec,
                                                   NeighborListKind neighbor_list,
                                                   EvaluateForce eval_frc, EvaluateEnergy eval_nrg,
                                                   TinyBoxPresence has_tiny_box,
                                                   ClashResponse clash_handling);
#endif

/// \brief Launch the appropriate kernel to evaluate particle-particle pair interactions in a
///        neighbor list for periodic simulations.
///
/// Overloaded:
///   - Supply abstracts, differentiating single- and double-precision calculations
///   - Supply the original objects, differentiating a unified neighbor list from separated
///     neighbor lists
///
/// \param poly_nbk
/// \param lemr
/// \param tlpn
/// \param nrg_tab
/// \param sysbrd
/// \param scw
/// \param cgw
/// \param cgw_qq
/// \param cgw_lj
/// \param ctrl
/// \param eval_frc
/// \param eval_nrg
/// \param has_tiny_box
/// \param bt
/// \param clash_distance
/// \param clash_ratio
/// \{
void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<double, double4> &nrg_tab,
                    CellGridWriter<double, llint, double, double4> *cgw, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<double> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<double, double4> &nrg_tab,
                    const PsSynthesisBorders &sysbrd,
                    CellGridWriter<double, llint, double, double4> *cgw, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<double> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<double, double4> &nrg_tab,
                    CellGridWriter<float, int, float, float4> *cgw, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<double> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<double, double4> &nrg_tab,
                    const PsSynthesisBorders &sysbrd,
                    CellGridWriter<float, int, float, float4> *cgw, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<double> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<double, double4> &nrg_tab,
                    CellGridWriter<double, llint, double, double4> *cgw_qq,
                    CellGridWriter<double, llint, double, double4> *cgw_lj, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<double> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<double, double4> &nrg_tab,
                    const PsSynthesisBorders &sysbrd,
                    CellGridWriter<double, llint, double, double4> *cgw_qq,
                    CellGridWriter<double, llint, double, double4> *cgw_lj, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<double> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<double, double4> &nrg_tab,
                    CellGridWriter<float, int, float, float4> *cgw_qq,
                    CellGridWriter<float, int, float, float4> *cgw_lj, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<double> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<double, double4> &nrg_tab,
                    const PsSynthesisBorders &sysbrd,
                    CellGridWriter<float, int, float, float4> *cgw_qq,
                    CellGridWriter<float, int, float, float4> *cgw_lj, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<double> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);
  
void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<float, float4> &nrg_tab,
                    CellGridWriter<float, int, float, float4> *cgw, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<float> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<float, float4> &nrg_tab,
                    const PsSynthesisBorders &sysbrd,
                    CellGridWriter<float, int, float, float4> *cgw, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<float> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);
  
void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<float, float4> &nrg_tab,
                    CellGridWriter<double, llint, double, double4> *cgw, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<float> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<float, float4> &nrg_tab,
                    const PsSynthesisBorders &sysbrd,
                    CellGridWriter<double, llint, double, double4> *cgw, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<float> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);
  
void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<float, float4> &nrg_tab,
                    CellGridWriter<float, int, float, float4> *cgw_qq,
                    CellGridWriter<float, int, float, float4> *cgw_lj, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<float> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<float, float4> &nrg_tab,
                    const PsSynthesisBorders &sysbrd,
                    CellGridWriter<float, int, float, float4> *cgw_qq,
                    CellGridWriter<float, int, float, float4> *cgw_lj, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<float> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<float, float4> &nrg_tab,
                    CellGridWriter<double, llint, double, double4> *cgw_qq,
                    CellGridWriter<double, llint, double, double4> *cgw_lj, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<float> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<float, float4> &nrg_tab,
                    const PsSynthesisBorders &sysbrd,
                    CellGridWriter<double, llint, double, double4> *cgw_qq,
                    CellGridWriter<double, llint, double, double4> *cgw_lj, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<float> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(PrecisionModel prec, const LocalExclusionMask &lem,
                    const PPITable &pairs_tbl, CellGrid<double, llint, double, double4> *cg,
                    TileManager *tlmn, ScoreCard *sc, MolecularMechanicsControls *mmctrl,
                    EvaluateForce eval_frc, EvaluateEnergy eval_nrg,
                    const CoreKlManager &launcher, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(PrecisionModel prec, const LocalExclusionMask &lem,
                    const PPITable &pairs_tbl, CellGrid<double, llint, double, double4> *cg_qq,
                    CellGrid<double, llint, double, double4> *cg_lj, TileManager *tlmn,
                    ScoreCard *sc, MolecularMechanicsControls *mmctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, const CoreKlManager &launcher,
                    double clash_distance = 0.0, double clash_ratio = 0.0);

void launchPMEPairs(PrecisionModel prec, const LocalExclusionMask &lem,
                    const PPITable &pairs_tbl, CellGrid<float, int, float, float4> *cg,
                    TileManager *tlmn, ScoreCard *sc, MolecularMechanicsControls *mmctrl,
                    EvaluateForce eval_frc, EvaluateEnergy eval_nrg,
                    const CoreKlManager &launcher, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(PrecisionModel prec, const LocalExclusionMask &lem,
                    const PPITable &pairs_tbl, CellGrid<float, int, float, float4> *cg_qq,
                    CellGrid<float, int, float, float4> *cg_lj, TileManager *tlmn, ScoreCard *sc,
                    MolecularMechanicsControls *mmctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, const CoreKlManager &launcher,
                    double clash_distance = 0.0, double clash_ratio = 0.0);
/// \}

} // namespace energy
} // namespace stormm

#endif
