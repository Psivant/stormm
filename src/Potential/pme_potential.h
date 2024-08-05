// -*-c++-*-
#ifndef STORMM_PME_POTENTIAL_H
#define STORMM_PME_POTENTIAL_H

#include "copyright.h"
#include "Constants/symbol_values.h"
#include "Math/series_ops.h"
#include "Math/vector_ops.h"
#include "Structure/local_arrangement.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/phasespace.h"
#include "local_exclusionmask.h"
#include "cellgrid.h"
#include "energy_enumerators.h"
#include "pme_util.h"
#include "scorecard.h"

namespace stormm {
namespace energy {

using stmath::findBin;
using stmath::indexingArray;
using symbols::amber_ancient_bioq;
using synthesis::AtomGraphSynthesis;
using synthesis::SyNonbondedKit;
using topology::AtomGraph;
using topology::NonbondedKit;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;

/// \brief Staggering constants for warps in the non-bonded particle-particle interactions.  Each
///        value is a 32-bit string.  Its high eight bits encode the baseline depth of the tile
///        (assuming that there is a complete batch of receiving atoms and therefore no further
///        subdivision of the tile occurs).  Bits 1-8, 9-16, and 17-24 encode the baseline stagger
///        for the second, third, and fourth quarters of the warp, respectively. (The stagger of
///        the first quarter of the warp is always zero.) If the tile is further subdivided by a
///        short batch of receiving atoms, the baseline staggers will advance and the baseline
///        depth of the tile evaluation will diminish as appropriate.
/// \{
#ifdef STORMM_USE_CUDA
constexpr uint no_stagger = (0x20 << 24);
constexpr uint half_stagger = ((0x10 << 24) | (0x10 << 16) | (0x10 << 8));
constexpr uint quarter_stagger = ((0x8 << 24) | (0x18 << 16) | (0x10 << 8) | 0x8);
#endif
/// \}
  
/// \brief Evaluate all interacting pairs between two cells of a simple neighbor list.  This is
///        only effective for a single system and is bound to perform double-precision
///        accumulation regardless of the precision of the calculations themselves.
///
/// \param psw          Abstract of the coordinates object containing positions, velocities,
///                     and forces of all particles in the system
/// \param cell_list    List of indices identifying the particle contents of each cell
/// \param cell_bounds  Bounds array for cell_list 
/// \param ncell_a      The number of neighbor list cells along the unit cell A axis
/// \param ncell_b      The number of neighbor list cells along the unit cell B axis
/// \param ncell_c      The number of neighbor list cells along the unit cell C axis
/// \param cell_idx_a   Position of the home cell along the system's A axis
/// \param cell_idx_b   Position of the home cell along the system's B axis
/// \param cell_idx_c   Position of the home cell along the system's C axis
/// \param i_pair       Relative position of the partner cell along the unit cell A axis
/// \param j_pair       Relative position of the partner cell along the unit cell B axis
/// \param k_pair       Relative position of the partner cell along the unit cell C axis
/// \param nbk          Non-bonded parameters for the system
/// \param lemr         Map of all exclusions based on relative particle indices
/// \param elec_cutoff  Particle-particle interaction cutoff for electrostatic interactions
/// \param vdw_cutoff   Particle-particle interaction cutoff for van-der Waals interactions
/// \param qqew_coeff   Ewald coefficient for electrostatic interactions
/// \param ljew_coeff   Ewald coefficient for van-der Waals interactions
/// \param vdw_sum      Method for computing the van-der Waals sum pairwise
/// \param eval_frc     Indicate whether to evaluate forces (energy computations are obligatory,
///                     as they are with other CPU-bound routines)
/// \param theme        The nature of particle-particle interactions to compute (electrostatic,
///                     van-der Waals, or both).
template <typename Tcalc>
double2 cellToCellInteractions(PhaseSpaceWriter *psw, const std::vector<int> &cell_list,
                               const std::vector<int> &cell_bounds, int ncell_a, int ncell_b,
                               int ncell_c, int cell_idx_a, int cell_idx_b, int cell_idx_c,
                               int i_pair, int j_pair, int k_pair,
                               const NonbondedKit<Tcalc> &nbk,
                               const LocalExclusionMaskReader &lemr, Tcalc elec_cutoff,
                               Tcalc vdw_cutoff, Tcalc qqew_coeff, Tcalc ljew_coeff,
                               VdwSumMethod vdw_sum = VdwSumMethod::CUTOFF,
                               EvaluateForce eval_frc = EvaluateForce::YES,
                               NonbondedTheme theme = NonbondedTheme::ALL);
  
/// \brief Evaluate the short-ranged, particle-particle energy of a system or group of systems.
///
/// Overloaded:
///   - Evaluate the energy of a single system
///   - Evaluate the energy of a group of systems
///   - Provide topology information via const pointer or const reference
///   - Provide abstracts to the relevant objects, or the objects themselves
///        
/// Descriptions of input parameters follow from cellToCellInteractions(), above, in addition to:
///
/// \param ps     Coordinates of a single system, including force accumulators suitable for
///               single-threaded accumulation
/// \param ag     Topology for a single system
/// \param cg     The cell grid, containing particle positions and force accumulators (in the
///               pre-arranged cell decomposition)
/// \param lema   Local exclusion masks for one structure or the synthesis of structures
/// \param prec   The precision model with which to carry out calculations
/// \param theme  The nature of particle-particle interactions to compute (electrostatic, van-der
///               Waals, or both).  While this information is present in the CellGrid object, an
///               indeoendent value is provided to this function which must be compatible with
///               the contents of the CellGrid object but can be used to test components of the
///               force and energy computations.
/// \{
template <typename Tcalc>
double2 evaluateParticleParticleEnergy(PhaseSpaceWriter *psw, const NonbondedKit<Tcalc> &nbk,
                                       const LocalExclusionMaskReader &lemr,
                                       Tcalc elec_cutoff = static_cast<Tcalc>(default_pme_cutoff),
                                       Tcalc vdw_cutoff = static_cast<Tcalc>(default_pme_cutoff),
                                       Tcalc qqew_coeff = static_cast<Tcalc>(0.0),
                                       Tcalc ljew_coeff = static_cast<Tcalc>(0.0),
                                       VdwSumMethod vdw_sum = VdwSumMethod::CUTOFF,
                                       EvaluateForce eval_frc = EvaluateForce::YES,
                                       NonbondedTheme theme = NonbondedTheme::ALL);

double2 evaluateParticleParticleEnergy(PhaseSpace *ps, const AtomGraph *ag,
                                       const LocalExclusionMask &lema,
                                       const PrecisionModel prec = PrecisionModel::SINGLE,
                                       double elec_cutoff = default_pme_cutoff,
                                       double vdw_cutoff = default_pme_cutoff,
                                       double qqew_coeff = 0.0, double ljew_coeff = 0.0,
                                       VdwSumMethod vdw_sum = VdwSumMethod::CUTOFF,
                                       EvaluateForce eval_frc = EvaluateForce::YES,
                                       NonbondedTheme theme = NonbondedTheme::ALL);

double2 evaluateParticleParticleEnergy(PhaseSpace *ps, const AtomGraph &ag,
                                       const LocalExclusionMask &lema,
                                       const PrecisionModel prec = PrecisionModel::SINGLE,
                                       double elec_cutoff = default_pme_cutoff,
                                       double vdw_cutoff = default_pme_cutoff,
                                       double qqew_coeff = 0.0, double ljew_coeff = 0.0,
                                       VdwSumMethod vdw_sum = VdwSumMethod::CUTOFF,
                                       EvaluateForce eval_frc = EvaluateForce::YES,
                                       NonbondedTheme theme = NonbondedTheme::ALL);

template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcalc2, typename Tcalc4>
void evaluateParticleParticleEnergy(const PhaseSpaceSynthesis &poly_ps,
                                    CellGridWriter<void, void, void, void> *cgw_v,
                                    ScoreCardWriter *scw, const PsSynthesisReader &poly_psr,
                                    const SyNonbondedKit<Tcalc, Tcalc2> &poly_nbk,
                                    const LocalExclusionMaskReader &lemr,
                                    Tcalc elec_cutoff = static_cast<Tcalc>(default_pme_cutoff),
                                    Tcalc vdw_cutoff = static_cast<Tcalc>(default_pme_cutoff),
                                    Tcalc qqew_coeff = static_cast<Tcalc>(0.0),
                                    Tcalc ljew_coeff = static_cast<Tcalc>(0.0),
                                    VdwSumMethod vdw_sum = VdwSumMethod::CUTOFF,
                                    EvaluateForce eval_frc = EvaluateForce::YES,
                                    NonbondedTheme theme = NonbondedTheme::ALL);
  
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcalc4>
void evaluateParticleParticleEnergy(CellGrid<Tcoord, Tacc, Tcalc, Tcalc4> *cg, ScoreCard *sc,
                                    const LocalExclusionMask &lema,
                                    Tcalc elec_cutoff = static_cast<Tcalc>(default_pme_cutoff),
                                    Tcalc vdw_cutoff = static_cast<Tcalc>(default_pme_cutoff),
                                    Tcalc qqew_coeff = static_cast<Tcalc>(0.0),
                                    Tcalc ljew_coeff = static_cast<Tcalc>(0.0),
                                    VdwSumMethod vdw_sum = VdwSumMethod::CUTOFF,
                                    EvaluateForce eval_frc = EvaluateForce::YES,
                                    NonbondedTheme theme = NonbondedTheme::ALL);
/// \}

/// \brief Compute all interactions in a tile.  This is helpful for encapsulating repetitive and
///        detailed code appearing in the tower-plate interaction protocol.  Descriptions of input
///        parameters follow from cellToCellInteractions(), above, in addition to:
///
/// \param a_xpos            Cartesian X positions of "a" atoms, shifted as required by their home
///                          cell positions
/// \param a_ypos            Cartesian Y positions of "a" atoms, shifted as required by home cells
/// \param a_zpos            Cartesian Z positions of "a" atoms, shifted as required by home cells
/// \param b_xpos            Cartesian X positions of "b" atoms, shifted as required by their home
///                          cell positions
/// \param b_ypos            Cartesian Y positions of "b" atoms, shifted as required by home cells
/// \param b_zpos            Cartesian Z positions of "b" atoms, shifted as required by home cells
/// \param scl_aq            Charges of the "a" atoms, scaled by the Coulomb constant
/// \param bq                Charges of the "b" atoms
/// \param ofs_aljidx        Lennard-Jones type indices of the "a" atoms, multiplied by the number
///                          of Lennard-Jones atom types in the system and then shifted according
///                          to the system offset within the concatenated tables of the synthesis
///                          such that the index of the "b" atom will step along the rows of each
///                          Lennard-Jones matrix to deliver the correct parameters for the a->b
///                          interaction
/// \param bljidx            Lennard-Jones type indices of the "b" atoms
/// \param top_aidx          Indices of the "a" atoms within the topology synthesis
/// \param top_bidx          Indices of the "b" atoms within the topology synthesis
/// \param img_aidx          Indices of the "a" atoms within the cell grid image
/// \param img_bidx          Indices of the "b" atoms within the cell grid image
/// \param na                Number of atoms in the "a" arrays
/// \param nb                Number of atoms in the "b" arrays
/// \param self_interaction  Set to TRUE if the list of atoms described by the "a" arrays is
///                          identical to the list of atoms described by the "b" arrays
/// \param elec_cutsq        Squared electrostatic interaction cutoff
/// \param vdw_cutsq         Squared van-der Waals interaction cutoff
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcalc2, typename Tcalc4>
double2 basicTileInteractions(const std::vector<Tcalc> &a_xpos, const std::vector<Tcalc> &a_ypos,
                              const std::vector<Tcalc> &a_zpos, const std::vector<Tcalc> &b_xpos,
                              const std::vector<Tcalc> &b_ypos, const std::vector<Tcalc> &b_zpos,
                              const std::vector<Tcalc> &scl_aq, const std::vector<Tcalc> &bq,
                              const std::vector<int> &ofs_aljidx, const std::vector<int> &bljidx,
                              const std::vector<int> &top_aidx, const std::vector<int> &top_bidx,
                              const std::vector<uint> &img_aidx, const std::vector<uint> &img_bidx,
                              int system_index, int na, int nb, bool self_interaction,
                              CellGridWriter<Tcoord, Tacc, Tcalc, Tcalc4> *cgw,
                              const PsSynthesisReader &poly_psr,
                              const SyNonbondedKit<Tcalc, Tcalc2> &poly_nbk,
                              const LocalExclusionMaskReader &lemr, Tcalc elec_cutsq,
                              Tcalc vdw_cutsq, Tcalc qqew_coeff, Tcalc ljew_coeff,
                              VdwSumMethod vdw_sum, EvaluateForce eval_frc, NonbondedTheme theme);

/// \brief Evaluate all interacting pairs in cells comprised by a neutral-territory tower/plate
///        decomposition.  Descriptions of input parameters follow from cellToCellInteractions(),
///        above, in addition to:
///
/// \param cgw         The cell grid for all systems
/// \param poly_psr    Global coordinates for all systems (needed in case small boxes require a
///                    re-imaging using the system's complete unit cell transformation matrices)
/// \param system_idx  The system within the synthesis to analyze
/// \param cell_idx    The home cell within the system, owning the Neutral Territory stencil.  The
///                    location within the cell grid will be deduced internally.
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcalc2, typename Tcoord4>
double2 towerPlatePairInteractions(CellGridWriter<Tcoord, Tacc, Tcalc, Tcoord4> *cgw,
                                   const PsSynthesisReader &poly_psr,
                                   int system_idx, int cell_idx,
                                   const SyNonbondedKit<Tcalc, Tcalc2> &poly_nbk,
                                   const LocalExclusionMaskReader &lemr, Tcalc elec_cutoff,
                                   Tcalc vdw_cutoff, Tcalc qqew_coeff, Tcalc ljew_coeff,
                                   VdwSumMethod vdw_sum, EvaluateForce eval_frc,
                                   NonbondedTheme theme);

} // namespace energy
} // namespace stormm

#include "pme_potential.tpp"

#endif
