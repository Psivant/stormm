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
/// Overloaded:
///   - Find the attributes of a kernel taking single- or double-precision coordinates
///   - Find the attributes of a kernel performing single- or double-precision calculations
///
/// \param coord_prec       The precision of the positions representaton in the neighbor list
/// \param calc_prec        Precision of the arithmetic
/// \param eval_frc         Indicate whether to evaluate forces in the kernel
/// \param eval_nrg         Indicate whether to evaluate energies in the kernel
/// \param neighbor_layout  The layout of the non-bonded neighbor list: one grid, or two?
/// \param has_tiny_box     Indicate whether one or more of the systems occur in a tiny box, with
///                         dimensions smaller than 2.5 times the respective cutoff between any two
///                         parallel faces
/// \param clash_handling   Indicate whether to forgive clashes between particles
/// \{
cudaFuncAttributes
queryDDPMEPairsKernelRequirements(NeighborListKind neighbor_layout, EvaluateForce eval_frc,
                                  EvaluateEnergy eval_nrg, TinyBoxPresence has_tiny_box,
                                  ClashResponse clash_handling);

cudaFuncAttributes
queryDFPMEPairsKernelRequirements(NeighborListKind neighbor_layout, EvaluateForce eval_frc,
                                  EvaluateEnergy eval_nrg, TinyBoxPresence has_tiny_box,
                                  ClashResponse clash_handling);

cudaFuncAttributes
queryFDPMEPairsKernelRequirements(NeighborListKind neighbor_layout, EvaluateForce eval_frc,
                                  EvaluateEnergy eval_nrg, TinyBoxPresence has_tiny_box,
                                  ClashResponse clash_handling);

cudaFuncAttributes
queryFFPMEPairsKernelRequirements(NeighborListKind neighbor_layout, EvaluateForce eval_frc,
                                  EvaluateEnergy eval_nrg, TinyBoxPresence has_tiny_box,
                                  ClashResponse clash_handling);

cudaFuncAttributes
queryPMEPairsKernelRequirements(PrecisionModel coord_prec, PrecisionModel calc_prec,
                                NeighborListKind neighbor_layout, EvaluateForce eval_frc,
                                EvaluateEnergy eval_nrg, TinyBoxPresence has_tiny_box,
                                ClashResponse clash_handling);
/// \}
#endif

/// \brief Launch the appropriate kernel to evaluate particle-particle pair interactions in a
///        neighbor list for periodic simulations.
///
/// Overloaded:
///   - Supply abstracts, differentiating single- and double-precision calculations
///   - Supply the original objects, differentiating a unified neighbor list from separated
///     neighbor lists
///
/// \param poly_nbk        Non-bonded parameters for all particles in the synthesis
/// \param lemr            Local exclusion mask abstract containing profile indices for all atoms
///                        and mask profiles for all systems
/// \param tlpn            Tile management object for directing thread communication while
///                        evaluating tiles
/// \param nrg_tab         Compiled table of splined interactions
/// \param sysbrd          System transformation matrices, needed when systems may exist in tiny
///                        unit cells.  This information comes from the coordinate synthesis, but
///                        is a much more compact abstract than the full PsSynthesisReader.
/// \param scw             Mutable abstract of the energy tracker, ready to accept contributions
/// \param cgw             Abstract of the neighbor list for all particle types
/// \param cgw_qq          Abstract of the neighbor list for electrostatic interactions
/// \param cgw_lj          Abstract of the neighbor list for van-der Waals interactions
/// \param ctrl            Holds counters for asynchronous advancement of work units
/// \param eval_frc        Indicate whether to launch a kernel evaluating forces between particles
/// \param eval_nrg        Indicate whether to launch a kernel that evaluates total system energy
/// \param has_tiny_box    Indicate whether one or more systems in the synthesis has a small box
///                        which requires special imaging considerations
/// \param bt              Block and thread counts for the kernel launch, obtained from a
///                        CoreKLManager (core kernel manager) object
/// \param clash_distance  The distance at which electrostatic interactions will be deemed in need
///                        of softening due to extreme proximity
/// \param clash_ratio     The ratio of Lennard-Jones pairwise sigma values at which a non-bonded
///                        Lennard-Jones interaction will be deemed in need of softening
/// \{
void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr,
                    const PPIKit<double, double4_16a> &nrg_tab,
                    CellGridWriter<double, llint, double, double4_16a> *cgw, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<double> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr,
                    const PPIKit<double, double4_16a> &nrg_tab, const PsSynthesisBorders &sysbrd,
                    CellGridWriter<double, llint, double, double4_16a> *cgw, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<double> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr,
                    const PPIKit<double, double4_16a> &nrg_tab,
                    CellGridWriter<float, int, float, float4> *cgw, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<double> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr,
                    const PPIKit<double, double4_16a> &nrg_tab, const PsSynthesisBorders &sysbrd,
                    CellGridWriter<float, int, float, float4> *cgw, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<double> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr,
                    const PPIKit<double, double4_16a> &nrg_tab,
                    CellGridWriter<double, llint, double, double4_16a> *cgw_qq,
                    CellGridWriter<double, llint, double, double4_16a> *cgw_lj, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<double> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr,
                    const PPIKit<double, double4_16a> &nrg_tab, const PsSynthesisBorders &sysbrd,
                    CellGridWriter<double, llint, double, double4_16a> *cgw_qq,
                    CellGridWriter<double, llint, double, double4_16a> *cgw_lj, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<double> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr,
                    const PPIKit<double, double4_16a> &nrg_tab,
                    CellGridWriter<float, int, float, float4> *cgw_qq,
                    CellGridWriter<float, int, float, float4> *cgw_lj, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<double> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<double, double2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr,
                    const PPIKit<double, double4_16a> &nrg_tab, const PsSynthesisBorders &sysbrd,
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
                    CellGridWriter<double, llint, double, double4_16a> *cgw, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<float> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<float, float4> &nrg_tab,
                    const PsSynthesisBorders &sysbrd,
                    CellGridWriter<double, llint, double, double4_16a> *cgw, TilePlan *tlpn,
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
                    CellGridWriter<double, llint, double, double4_16a> *cgw_qq,
                    CellGridWriter<double, llint, double, double4_16a> *cgw_lj, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<float> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(const SyNonbondedKit<float, float2> &poly_nbk,
                    const LocalExclusionMaskReader &lemr, const PPIKit<float, float4> &nrg_tab,
                    const PsSynthesisBorders &sysbrd,
                    CellGridWriter<double, llint, double, double4_16a> *cgw_qq,
                    CellGridWriter<double, llint, double, double4_16a> *cgw_lj, TilePlan *tlpn,
                    ScoreCardWriter *scw, MMControlKit<float> *ctrl, EvaluateForce eval_frc,
                    EvaluateEnergy eval_nrg, int2 bt_tp, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(PrecisionModel prec, const LocalExclusionMask &lem,
                    const PPITable &pairs_tbl, CellGrid<double, llint, double, double4_16a> *cg,
                    TileManager *tlmn, ScoreCard *sc, MolecularMechanicsControls *mmctrl,
                    EvaluateForce eval_frc, EvaluateEnergy eval_nrg,
                    const CoreKlManager &launcher, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void launchPMEPairs(PrecisionModel prec, const LocalExclusionMask &lem,
                    const PPITable &pairs_tbl, CellGrid<double, llint, double, double4_16a> *cg_qq,
                    CellGrid<double, llint, double, double4_16a> *cg_lj, TileManager *tlmn,
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
