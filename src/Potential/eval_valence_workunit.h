// -*-c++-*-
#ifndef STORMM_EVAL_VALENCE_WORKUNIT_H
#define STORMM_EVAL_VALENCE_WORKUNIT_H

#include "copyright.h"
#include "Restraints/restraint_apparatus.h"
#include "Synthesis/valence_workunit.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/coordinateframe.h"
#include "energy_enumerators.h"
#include "scorecard.h"

namespace stormm {
namespace energy {

using restraints::RestraintApparatus;
using restraints::RestraintKit;
using synthesis::ValenceWorkUnit;
using synthesis::VwuTask;
using synthesis::VwuGoal;
using topology::AtomGraph;
using topology::UnitCellType;
using topology::NonbondedKit;
using topology::ValenceKit;
using topology::VirtualSiteKit;
using trajectory::CoordinateFrame;
using trajectory::PhaseSpace;

/// \brief Carry out the local evaluations need for each valence work unit using locally cached
///        data.  This routine makes no assumptions about const-ness of the coordinates to permit
///        maximum flexibility.
///
/// \param vk           Valence parameters for the system as a whole
/// \param vsk          Virtual site frame information
/// \param nbk          Non-bonded parameters for the system (for 1:4 interactions)
/// \param rar          Restraint apparatus pointers
/// \param sh_charges   Charges of cached particles
/// \param sh_lj_idx    Lennard Jones types indices of cached particles
/// \param sh_xcrd      Cached Cartesian X positions of all particles, imaged such that unit cell
///                     transformations are irrelevant
/// \param sh_ycrd      Cached Cartesian Y positions
/// \param sh_zcrd      Cached Cartesian Z positions
/// \param sh_xfrc      Cached Cartesian X forces for all particles (allocated and set to 0.0
///                     before input, ready to accumulate and pass back to the calling function)
/// \param sh_yfrc      Cached Y forces (allocated and set to 0.0 before input)
/// \param sh_zfrc      Cached Z forces (allocated and set to 0.0 before input)
/// \param ecard        Report comprising all energy terms
/// \param sysid        System ID number for catalogging purposes in the ScoreCard object
/// \param my_vwu       The specific work unit to evaluate (atomic coordinates will have been
///                     pre-cached, and forces initialized)
/// \param eval_force   Flag to indicate that forces are desired
/// \param activity     The specific tasks within each work unit to carry out (or, ALL_TASKS)
/// \param purpose      Reason for evaluatin the work units--one can compute global forces or move
///                     particles, but not both at once
/// \param step_number  MD or minimization step number (for calculating restraint activation)
void localVwuEvaluation(const ValenceKit<double> vk, const VirtualSiteKit<double> vsk,
                        const NonbondedKit<double> nbk,
                        const RestraintKit<double, double2, double4_16a> rar,
                        const double* sh_charges, const int* sh_lj_idx, double* sh_xcrd,
                        double* sh_ycrd, double* sh_zcrd, double* sh_xfrc, double* sh_yfrc,
                        double* sh_zfrc, ScoreCard *ecard, int sysid,
                        const ValenceWorkUnit &my_vwu, EvaluateForce eval_force, VwuTask activity,
                        VwuGoal purpose, int step_number);

/// \brief Evaluate force and energy-related tasks in a list of valence work units given a system
///        or synthesis of systems with one or more topologies and coordinate sets.
///
/// Overloaded:
///   - Evaluate for a single topology and coordinate set
///   - Evaluate for a synthesis of topologies and coordinate sets
///   - Accept pointers or references to the above objects, or in the case of a single system
///     offer the option of passing critical topology abstracts by value with pointers to C-style
///     arrays for coordinates and forces
///   - Evaluate forces, or energies alone
///
/// Parameters for these function follow from localVwuEvaluation() above, except for:
///
/// \param xcrd  Cartesian X coordinates of all particles in the system (the correct import list
///              will be cached and passed to the local evaluation)
/// \param ycrd  Cartesian Y coordinates of all particles in the system
/// \param ycrd  Cartesian Z coordinates of all particles in the system
/// \param xfrc  Global storage array of Cartesian X forces on all particles
/// \param yfrc  Global storage array of Cartesian Y forces on all particles
/// \param zfrc  Global storage array of Cartesian Z forces on all particles
/// \{
void evalValenceWorkUnits(const ValenceKit<double> vk, const VirtualSiteKit<double> vsk,
                          const NonbondedKit<double> nbk,
                          const RestraintKit<double, double2, double4_16a> rar, double* xcrd,
                          double* ycrd, double* zcrd, const double* umat, const double* invu,
                          UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                          ScoreCard *ecard, int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          EvaluateForce eval_force = EvaluateForce::NO,
                          VwuTask activity = VwuTask::ALL_TASKS,
                          VwuGoal goal = VwuGoal::ACCUMULATE, int step_number = 0);

void evalValenceWorkUnits(const AtomGraph *ag, PhaseSpace *ps, const RestraintApparatus *ra,
                          ScoreCard *ecard, int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          EvaluateForce eval_force = EvaluateForce::NO,
                          VwuTask activity = VwuTask::ALL_TASKS,
                          VwuGoal goal = VwuGoal::ACCUMULATE, int step_number = 0);

void evalValenceWorkUnits(const ValenceKit<double> vk, const VirtualSiteKit<double> vsk,
                          const NonbondedKit<double> nbk,
                          const RestraintKit<double, double2, double4_16a> rar, const double* xcrd,
                          const double* ycrd, const double* zcrd, const double* umat,
                          const double* invu, UnitCellType unit_cell, double* xfrc, double* yfrc,
                          double* zfrc, ScoreCard *ecard, const int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list, EvaluateForce eval_force,
                          VwuTask activity, int step_number);

void evalValenceWorkUnits(const AtomGraph &ag, const PhaseSpace &ps, const RestraintApparatus &ra,
                          ScoreCard *ecard, int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          VwuTask activity = VwuTask::ALL_TASKS, int step_number = 0);
/// \}

} // namespace energy
} // namespace stormm

#endif
