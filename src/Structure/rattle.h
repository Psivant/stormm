// -*-c++-*-
#ifndef STORMM_RATTLE_H
#define STORMM_RATTLE_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/common_types.h"
#include "Namelists/nml_dynamics.h"
#include "Numerics/split_fixed_precision.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/trajectory_enumerators.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

using constants::PrecisionModel;
using data_types::isSignedIntegralScalarType;
using namelist::default_rattle_max_iter;
using synthesis::AtomGraphSynthesis;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using synthesis::SyAtomUpdateKit;
using synthesis::SyValenceKit;
using synthesis::VwuAbstractMap;
using synthesis::vwu_abstract_length;
using topology::AtomGraph;
using topology::ConstraintKit;
using trajectory::CoordinateCycle;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;

/// \brief Apply RATTLE constraints to bonds in one or more topologies, in fulfillment of the
///        coordinate update constraints.  This is applied after the Velocity-Verlet positions
///        and second half-step velocity update.  The corrections are applied to positions "under
///        construction," e.g. the arrays xalt, yalt, and zalt in a PhaseSpace or
///        PhaseSpaceSynthesis abstract obtained for the original coordinate object's current time
///        frame.
///
/// Overloaded:
///   - Provide a single system's PhaseSpace object, with the appropriate topological constraints
///   - Provide a synthesis of coordinates (PhaseSpaceSynthesis), with the corresponding
///     topological synthesis
///   - Provide the original objects or their abstracts
///
/// \param xdev  The developing Cartesian X coordinates of all particles.  These are taken from
///              arrays such as xalt in a PhaseSpaceWriter.
/// \param ydev  The developing Cartesian Y coordinates of all particles
/// \param zdev  The developing Cartesian Z coordinates of all particles
/// \param xref  The reference Cartesian X coordinates of all particles
/// \param yref  The reference Cartesian Y coordinates of all particles
/// \param zref  The reference Cartesian Z coordinates of all particles
/// \param ps    The coordinates and velocities of particles in the system
/// \param psw   Mutable abstract of the coordinates and velocities
/// \param ag    The topology guiding motion in the system
/// \param cnk   Constraints abstract from the topology
/// \param dt    The time step, necessary for scaling the postional adjustments against the length
///              of the time step to make a post-hoc velocity correction
/// \param prec  The precision model in which to perform calculations
/// \{
template <typename Tcoord, typename Tcalc>
void rattlePositions(Tcoord* xdev, Tcoord* ydev, Tcoord* zdev, Tcoord* xvel_dev, Tcoord* yvel_dev,
                     Tcoord* zvel_dev, const Tcoord* xref, const Tcoord* yref, const Tcoord* zref,
                     const ConstraintKit<Tcalc> &cnk, Tcalc dt, Tcalc tol, int max_iter,
                     RattleMethod style, Tcalc gpos_scale_factor = 1.0,
                     Tcalc vel_scale_factor = 1.0);

template <typename T>
void rattlePositions(PhaseSpaceWriter *psw, const ConstraintKit<T> &cnk, T dt, T tol,
                     int max_iter = default_rattle_max_iter,
                     RattleMethod style = RattleMethod::SEQUENTIAL);

void rattlePositions(PhaseSpace *ps, const AtomGraph *ag, PrecisionModel prec, double dt,
                     double tol, int max_iter = default_rattle_max_iter,
                     RattleMethod style = RattleMethod::SEQUENTIAL);

void rattlePositions(PhaseSpace *ps, const AtomGraph &ag, PrecisionModel prec, double dt,
                     double tol, int max_iter = default_rattle_max_iter,
                     RattleMethod style = RattleMethod::SEQUENTIAL);

template <typename T, typename T2, typename T4>
void rattlePositions(PsSynthesisWriter *poly_psw, const SyValenceKit<T> &poly_vk,
                     const SyAtomUpdateKit<T, T2, T4> &poly_auk, T dt, T tol,
                     int max_iter = default_rattle_max_iter);

void rattlePositions(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                     PrecisionModel prec, double dt, double tol,
                     int max_iter = default_rattle_max_iter);

void rattlePositions(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                     PrecisionModel prec, double dt, double tol,
                     int max_iter = default_rattle_max_iter);
/// \}

/// \brief Apply the RATTLE constraints to bonds in one or more topologies, in fulfillment of the
///        velocity update constraints.  This is applied after the Velocity-Verlet force
///        calculation and subsequent half-step velocity update.  The corrections are applied to
///        the velocities "under construction", e.g. the arrays vxalt, vyalt, and vzalt in a
///        PhaseSpace or PhaseSpaceSynthesis abstract obtained for the original coordinate object's
///        current time frame.  Overloads and descriptions of input arguments follow from
///        rattlePositions(), above.
/// \{
template <typename Tcoord, typename Tcalc>
void rattleVelocities(Tcoord* xvel_dev, Tcoord* yvel_dev, Tcoord *zvel_dev, const Tcoord* xcrd_ref,
                      const Tcoord* ycrd_ref, const Tcoord* zcrd_ref,
                      const ConstraintKit<Tcalc> &cnk, Tcalc dt, Tcalc tol,
                      int max_iter = default_rattle_max_iter,
                      RattleMethod style = RattleMethod::SEQUENTIAL, Tcalc gpos_scale_factor = 1.0,
                      Tcalc vel_scale_factor = 1.0);

template <typename Tcalc>
void rattleVelocities(PhaseSpaceWriter *psw, const ConstraintKit<Tcalc> &cnk, Tcalc dt, Tcalc tol,
                      int max_iter = default_rattle_max_iter,
                      RattleMethod style = RattleMethod::SEQUENTIAL);

void rattleVelocities(PhaseSpace *ps, const AtomGraph *ag, PrecisionModel prec, double dt,
                      double tol, int max_iter = default_rattle_max_iter,
                      RattleMethod style = RattleMethod::SEQUENTIAL);

void rattleVelocities(PhaseSpace *ps, const AtomGraph &ag, PrecisionModel prec, double dt,
                      double tol, int max_iter = default_rattle_max_iter,
                      RattleMethod style = RattleMethod::SEQUENTIAL);

template <typename T, typename T2, typename T4>
void rattleVelocities(PsSynthesisWriter *poly_psw, const SyAtomUpdateKit<T, T2, T4> &poly_auk,
                      T dt, T tol, int max_iter = default_rattle_max_iter);

void rattleVelocities(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                      PrecisionModel prec, double dt, double tol,
                      int max_iter = default_rattle_max_iter);

void rattleVelocities(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                      PrecisionModel prec, double dt, double tol,
                      int max_iter = default_rattle_max_iter);
/// \}
  
} // namespace structure
} // namespace stormm

#include "rattle.tpp"

#endif
