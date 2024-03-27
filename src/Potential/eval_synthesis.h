// -*-c++-*-
#ifndef STORMM_EVAL_SYNTHESIS_H
#define STORMM_EVAL_SYNTHESIS_H

#include "copyright.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/vector_ops.h"
#include "Potential/eval_valence_workunit.h"
#include "Potential/scorecard.h"
#include "Potential/static_exclusionmask.h"
#include "Potential/valence_potential.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/nonbonded_workunit.h"
#include "Synthesis/static_mask_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Synthesis/valence_workunit.h"

namespace stormm {
namespace energy {
  
using stmath::readBitFromMask;
using synthesis::AtomGraphSynthesis;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::small_block_max_imports;
using synthesis::supertile_wu_abstract_length;
using synthesis::tile_groups_wu_abstract_length;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using synthesis::SeMaskSynthesisReader;
using synthesis::StaticExclusionMaskSynthesis;
using synthesis::SyAtomUpdateKit;
using synthesis::SyNonbondedKit;
using synthesis::SyValenceKit;
using synthesis::SyRestraintKit;
using synthesis::vwu_abstract_length;
using synthesis::VwuAbstractMap;
using synthesis::VwuGoal;
using synthesis::VwuTask;

/// \brief Carry out the instructions in a single valence work unit, as presented in the topology
///        synthesis.  This routine is called by various incarnations of the evalSyValenceEnergy()
///        function.  While any of the local position, velocity, or force data passed into it
///        could be modified, the calling functions may impose const-ness on the global positions,
///        velocities, and forces of various particles.
///
/// \param syvk         Consensus tables of valence parameters and instructions
/// \param syrk         Consensus tables of restraint parameters and instructions  
/// \param sh_xcrd      Mock data array for locally cached particle X coordinates
/// \param sh_ycrd      Mock data array for locally cached particle Y coordinates
/// \param sh_zcrd      Mock data array for locally cached particle Z coordinates
/// \param sh_xvel      Mock data array for locally cached particle X velocities
/// \param sh_yvel      Mock data array for locally cached particle Y velocities
/// \param sh_zvel      Mock data array for locally cached particle Z velocities
/// \param sh_xfrc      Mock data array for locally accumulating particle X forces
/// \param sh_yfrc      Mock data array for locally accumulating particle Y forces
/// \param sh_zfrc      Mock data array for locally accumulating particle Z forces
/// \param ecard        Energy tracking object
/// \param vwu_idx      Index of the valence work unit to evaluate, as stored in the topology
///                     synthesis
/// \param eval_force   Flag to have forces evaluated
/// \param activity     Evaluate a particular energy component, or all components
/// \param purpose      Purpose of the evaluation: to accumulate forces and / or energies, or to
///                     move particles
/// \param step_number  Number of the step in the simulation (relevant to restraint applications)
template <typename Tcalc, typename Tcalc2, typename Tcalc4>
void synthesisVwuEvaluation(const SyValenceKit<Tcalc> syvk,
                            const SyRestraintKit<Tcalc, Tcalc2, Tcalc4> syrk,
                            const Tcalc* sh_charges, const int* sh_lj_idx, llint* sh_xcrd,
                            llint* sh_ycrd, llint* sh_zcrd, llint* sh_xvel, llint* sh_yvel,
                            llint* sh_zvel, llint* sh_xfrc, llint* sh_yfrc, llint* sh_zfrc,
                            double inv_gpos_scale, double force_scale, ScoreCard *ecard,
                            int vwu_idx, EvaluateForce eval_force, VwuTask activity,
                            VwuGoal purpose, int step_number);

/// \brief Evaluate all work units in an AtomGraphSynthesis (synthesis of topologies).  This
///        function will allocate mock data arrays to drive synthesisVwuEvaluation() above.
///
/// \param syvk         Topology synthesis-based abstract for parameters on valence interactions
/// \param syrk         Topology synthesis-based abstract for parameters on restraints
/// \param psyw         Writeable abstract for the coordinate synthesis
/// \param ecard        Energy tracker object
/// \param eval_force   Flag to also carry out force evaluation (energy is always evaluated in CPU
///                     functions)
/// \param activity     Evaluate a particular energy component, or all components
/// \param purpose      Purpose of the evaluation: to accumulate forces and / or energies, or to
///                     move particles
/// \param step_number  Number of the step in the simulation (relevant to restraint applications)
template <typename Tcalc, typename Tcalc2, typename Tcalc4>
void evalSyValenceEnergy(const SyValenceKit<Tcalc> syvk,
                         const SyAtomUpdateKit<Tcalc, Tcalc2, Tcalc4> syauk,
                         const SyRestraintKit<Tcalc, Tcalc2, Tcalc4> syrk, PsSynthesisWriter psyw,
                         ScoreCard *ecard, EvaluateForce eval_force, VwuTask activity,
                         VwuGoal purpose, int step_number);

/// \brief Evaluate the non-bonded energies (and possibly forces) of a synthesis of systems in
///        isolated boundary conditions using non-bonded work units composed of tile groups.  These
///        smaller forms of the all-to-all non-bonded work units enumerate all of the tiles they
///        process.
///
/// \param synbk            Non-bonded parameters for all atoms in the compilation of systems
/// \param psyr             Abstract for the coordinate synthesis
/// \param ecardw           Energy tracker object writer
/// \param eval_elec_force  Flag to also carry out force evaluation of electrostatic interactions
///                         (energy is always evaluated in CPU functions)
/// \param eval_vdw_force   Flag to also carry out force evaluation of van-der Waals interactions
template <typename Tcalc, typename Tcalc2>
void evalSyNonbondedTileGroups(const SyNonbondedKit<Tcalc, Tcalc2> synbk,
                               const SeMaskSynthesisReader syse, PsSynthesisWriter *psyw,
                               ScoreCard *ecard, EvaluateForce eval_elec_force,
                               EvaluateForce eval_vdw_force);

/// \brief Evaluate the non-bonded energy with a particular precision level.  This will invoke
///        the proper C++ function.
template <typename Tcalc>
void evalSyNonbondedEnergy(const AtomGraphSynthesis &poly_ag,
                           const StaticExclusionMaskSynthesis &poly_se,
                           PhaseSpaceSynthesis *poly_ps, ScoreCard *ecard,
                           EvaluateForce eval_elec_force, EvaluateForce eval_vdw_force);

} // namespace energy
} // namespace stormm

#include "eval_synthesis.tpp"

#endif
