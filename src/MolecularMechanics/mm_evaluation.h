// -*-c++-*-
#ifndef STORMM_MM_EVALUATION_H
#define STORMM_MM_EVALUATION_H

#include "copyright.h"
#include "Constants/generalized_born.h"
#include "Potential/energy_enumerators.h"
#include "Potential/nonbonded_potential.h"
#include "Potential/scorecard.h"
#include "Potential/static_exclusionmask.h"
#include "Potential/valence_potential.h"
#include "Restraints/restraint_apparatus.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"

namespace stormm {
namespace mm {

using energy::EvaluateForce;
using energy::evaluateAngleTerms;
using energy::evaluateAttenuated14Terms;
using energy::evaluateBondTerms;
using energy::evaluateCharmmImproperTerms;
using energy::evaluateCmapTerms;
using energy::evaluateDihedralTerms;
using energy::evaluateGeneralizedBornEnergy;
using energy::evaluateNonbondedEnergy;
using energy::evaluateRestraints;
using energy::evaluateUreyBradleyTerms;
using energy::ScoreCard;
using energy::StaticExclusionMask;
using energy::StaticExclusionMaskReader;
using restraints::RestraintApparatus;
using restraints::RestraintKit;
using topology::AtomGraph;
using topology::ImplicitSolventKit;
using topology::NonbondedKit;
using topology::UnitCellType;
using topology::ValenceKit;
using topology::VirtualSiteKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;
using namespace generalized_born_defaults;
  
/// \brief Evaluate molecular mechanics valence energies and forces in a single system.  Energy
///        results are collected in fixed precision in the ScoreCard object fed to this routine.
///
/// Overloaded:
///   - Take raw pointers to coordinates, box transformations, and forces, along with abstracts
///     to the required parameters
///   - Take a PhaseSpace object and the relevant topology, or abstracts thereof
///
/// \param xcrd             Cartesian X positions of all particles
/// \param ycrd             Cartesian Y positions of all particles
/// \param zcrd             Cartesian Z positions of all particles
/// \param umat             Transformation matrix taking coordinates into fractional (unit cell)
///                         space
/// \param invu             Transformation matrix to take fractional coordinates back into real
///                         space
/// \param unit_cell        The type of simulation cell, to guide re-imaging considerations
/// \param xfrc             Forces acting on all particles in the Cartesian X direction
///                         (accumulated and returned)
/// \param yfrc             Forces acting on all particles in the Cartesian Y direction
/// \param zfrc             Forces acting on all particles in the Cartesian Z direction
/// \param ps               PhaseSpace object with all coordinates, box information, and forces
/// \param psw              PhaseSpace object abstract with coordinates, box information, and
///                         forces
/// \param sc               Object to hold the resulting energies (also keeps running sums)
/// \param vk               Valence parameters abstract from the original topology
/// \param nbk              Non-bonded parameter abstract from the original topology
/// \param ag               System topology from which relevant parameter abstracts can be obtained
/// \param eval_force       Directive to have forces evaluated or not
/// \param system_index     Index of the system within the score card
/// \param inv_gpos_factor  De-scaling factor to apply to fixed-precision coordinates to get them
///                         back in units of Angstroms
/// \param force factor     Scaling factor to apply to real-values forces for fixed-precision
///                         accumulation
/// \param clash_distance   Minimum absolute separation, in Angstroms, below which two particles
///                         can be considered to be clashing
/// \param clash_ratio      Minimum ratio of the separation distance to the particle pair's van-der
///                         Waals sigma parameter, below which the two particles can be consdiered
///                         to be clashing
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc>
void evalValeMM(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                const double* invu, UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                Tforce* zfrc, ScoreCard *sc, const ValenceKit<Tcalc> &vk,
                const NonbondedKit<Tcalc> &nbk, EvaluateForce eval_force, int system_index = 0,
                Tcalc inv_gpos_factor = 1.0, Tcalc force_factor = 1.0, Tcalc clash_distance = 0.0,
                Tcalc clash_ratio = 0.0);

void evalValeMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                const NonbondedKit<double> &nbk, EvaluateForce eval_force, int system_index = 0,
                double clash_distance = 0.0, double clash_ratio = 0.0);

void evalValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag, EvaluateForce eval_force,
                int system_index = 0, double clash_distance = 0.0, double clash_ratio = 0.0);

void evalValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag, EvaluateForce eval_force,
                int system_index = 0, double clash_distance = 0.0, double clash_ratio = 0.0);
/// \}

/// \brief Evaluate molecular mechanics energies and forces due to valence interactions and
///        NMR restraints in a single system.  Energy results are collected in fixed precision in
///        the ScoreCard object fed to this routine.
///
/// Overloaded:
///   - Take raw pointers to coordinates, box transformations, and forces, along with abstracts
///     to the required parameters
///   - Take a PhaseSpace object and the relevant topology plus restraint collection, or abstracts
///     thereof
///
/// Descriptions of input variables follow from evalValeMM() above, with the addition of:
///
/// \param rar   Restraint apparatus abstract guiding NMR restraint potentials, supplemental to
///              the system's topology
/// \param ra    Restraint apparatus from which an appropriate abstract can be derived
/// \param step  The step number, which may affect restraint activation
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
void evalValeRestMM(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                    const double* invu, UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                    Tforce* zfrc, ScoreCard *sc, const ValenceKit<Tcalc> &vk,
                    const NonbondedKit<Tcalc> &nbk, const RestraintKit<Tcalc, Tcalc2, Tcalc4> &rar,
                    EvaluateForce eval_force, int system_index = 0, int step = 0,
                    Tcalc inv_gpos_factor = 1.0, Tcalc force_factor = 1.0,
                    Tcalc clash_distance = 0.0, Tcalc clash_ratio = 0.0);

void evalValeRestMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                    const NonbondedKit<double> &nbk,
                    const RestraintKit<double, double2, double4> &rar, EvaluateForce eval_force,
                    int system_index = 0, int step = 0, double clash_distance = 0.0,
                    double clash_ratio = 0.0);

void evalValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                    const RestraintApparatus &ra, EvaluateForce eval_force, int system_index = 0,
                    int step = 0, double clash_distance = 0.0, double clash_ratio = 0.0);

void evalValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                    const RestraintApparatus &ra, EvaluateForce eval_force, int system_index = 0,
                    int step = 0, double clash_distance = 0.0, double clash_ratio = 0.0);
/// \}
  
/// \brief Evaluate the molecular mechanics energies and forces due to valence and non-bonded
///        pair interactions.  Energy results are collected in fixed precision in the ScoreCard
///        object fed to this routine.
///
/// Overloaded:
///   - Take raw pointers to coordinates, box transformations, and forces, along with abstracts
///     to the required parameters
///   - Take a PhaseSpace object and the relevant topology plus restraint collection, or abstracts
///     thereof
///
/// Descriptions of input variables follow from evalValeMM() and evalValeRestMM() above.
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc>
void evalNonbValeMM(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                    const double* invu, UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                    Tforce* zfrc, ScoreCard *sc, const ValenceKit<Tcalc> &vk,
                    const NonbondedKit<Tcalc> &nbk, const StaticExclusionMaskReader &ser,
                    EvaluateForce eval_force, int system_index = 0, Tcalc inv_gpos_factor = 1.0,
                    Tcalc force_factor = 1.0, Tcalc clash_distance = 0.0, Tcalc clash_ratio = 0.0);

void evalNonbValeMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                    const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                    EvaluateForce eval_force, int system_index = 0,
                    double clash_distance = 0.0, double clash_ratio = 0.0);

void evalNonbValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                    const StaticExclusionMask &se, EvaluateForce eval_force, int system_index = 0,
                    double clash_distance = 0.0, double clash_ratio = 0.0);

void evalNonbValeMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                    const StaticExclusionMask &se, EvaluateForce eval_force, int system_index = 0,
                    double clash_distance = 0.0, double clash_ratio = 0.0);
/// \}
  
/// \brief Evaluate molecular mechanics energies and forces due to valence interactions, NMR
///        restraints, and non-bonded interactions in a single system.  Energy results are
///        collected in fixed precision in the ScoreCard object fed to this routine.
///
/// Overloaded:
///   - Take raw pointers to coordinates, box transformations, and forces, along with abstracts
///     to the required parameters
///   - Take a PhaseSpace object and the relevant topology plus restraint collection, or abstracts
///     thereof
///
/// Descriptions of input variables follow from evalValeMM() and evalValeRestMM() above, with the
/// addition of:
///
/// \param ser  Abstract for a map of non-bonded exclusions for all pairs of atoms in the system
/// \param se   Permanent map of non-bonded exclusions for all pairs of atoms in the system
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
void evalNonbValeRestMM(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                        const double* umat, const double* invu, UnitCellType unit_cell,
                        Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, ScoreCard *sc,
                        const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                        const StaticExclusionMaskReader &ser,
                        const RestraintKit<Tcalc, Tcalc2, Tcalc4> &rar, EvaluateForce eval_force,
                        int system_index = 0, int step = 0, Tcalc inv_gpos_factor = 1.0,
                        Tcalc force_factor = 1.0, Tcalc clash_distance = 0.0,
                        Tcalc clash_ratio = 0.0);

void evalNonbValeRestMM(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                        const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                        const RestraintKit<double, double2, double4> &rar,
                        EvaluateForce eval_force, int system_index = 0, int step = 0,
                        double clash_distance = 0.0, double clash_ratio = 0.0);

void evalNonbValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                        const StaticExclusionMask &se, const RestraintApparatus &ra,
                        EvaluateForce eval_force, int system_index = 0, int step = 0,
                        double clash_distance = 0.0, double clash_ratio = 0.0);

void evalNonbValeRestMM(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                        const StaticExclusionMask &se, const RestraintApparatus *ra,
                        EvaluateForce eval_force, int system_index = 0, int step = 0,
                        double clash_distance = 0.0, double clash_ratio = 0.0);
/// \}

/// \brief Evaluate the complete molecular mechanics energies and forces for a system in isolated
///        boundary conditions, including valence terms, restraints, non-bonded interactions, and
///        Generalized Born implicit solvent contributions.
///
/// Overloaded:
///   - Take raw pointers to coordinates, box transformations, and forces, along with abstracts
///     to the required parameters
///   - Take a PhaseSpace object and the relevant topology plus restraint collection, or abstracts
///     thereof
///
/// Descriptions of input variables follow from evalValeMM(), evalValeRestMM(), and
/// evalNonbValeRestMM() above, with the addition of:
///
/// \param isk         A recipe for evaluating implicit solvent contributions constructed from the
///                    topology (containing the implicit solvent model) and a Generalized Born
///                    parameter table
/// \param neck_gbk    Abstract of a neck GB table at the appropriate (double) precision level
/// \param neck_gbtab  Table of "neck" Generalized Born parameters, if applicable
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
void evalRestrainedMMGB(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                        const double* umat, const double* invu, UnitCellType unit_cell,
                        Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, ScoreCard *sc,
                        const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                        const StaticExclusionMaskReader &ser, const ImplicitSolventKit<Tcalc> &isk,
                        const NeckGeneralizedBornKit<Tcalc> &neck_gbk,
                        Tforce* effective_gb_radii, Tforce* psi, Tforce* sumdeijda,
                        const RestraintKit<Tcalc, Tcalc2, Tcalc4> &rar, EvaluateForce eval_force,
                        int system_index = 0, int step = 0, Tcalc inv_gpos_factor = 1.0,
                        Tcalc force_factor = 1.0, Tcalc clash_distance = 0.0,
                        Tcalc clash_ratio = 0.0);

void evalRestrainedMMGB(PhaseSpaceWriter psw, ScoreCard *sc, const ValenceKit<double> &vk,
                        const NonbondedKit<double> &nbk, const StaticExclusionMaskReader &ser,
                        const ImplicitSolventKit<double> &isk,
                        const NeckGeneralizedBornKit<double> &neck_gbk,
                        const RestraintKit<double, double2, double4> &rar,
                        EvaluateForce eval_force, int system_index = 0, int step = 0,
                        double clash_distance = 0.0, double clash_ratio = 0.0);

void evalRestrainedMMGB(PhaseSpace *ps, ScoreCard *sc, const AtomGraph &ag,
                        const NeckGeneralizedBornTable &neck_gbtab,
                        const StaticExclusionMask &se, const RestraintApparatus &ra,
                        EvaluateForce eval_force, int system_index = 0, int step = 0,
                        double clash_distance = 0.0, double clash_ratio = 0.0);
                        
void evalRestrainedMMGB(PhaseSpace *ps, ScoreCard *sc, const AtomGraph *ag,
                        const NeckGeneralizedBornTable &neck_gbtab,
                        const StaticExclusionMask &se, const RestraintApparatus *ra,
                        EvaluateForce eval_force, int system_index = 0, int step = 0,
                        double clash_distance = 0.0, double clash_ratio = 0.0);
/// \}
  
} // namespace mm
} // namespace stormm

#include "mm_evaluation.tpp"

#endif
