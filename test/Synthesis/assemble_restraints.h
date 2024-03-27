// -*-c++-*-
#include "../../src/Restraints/restraint_apparatus.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/phasespace.h"

using stormm::restraints::RestraintApparatus;
using stormm::topology::AtomGraph;
using stormm::trajectory::PhaseSpace;

/// \brief Create a collection of positional, distance, angle, and dihedral restraints acting on a
///        system which will exert mild to moderate forces.  All of the restraints, together, will
///        be checked for total energy and forces.  This feeds into tests that are not designed to
///        weed out individual problems but instead take a "pooled sample" and verify that there
///        are no problems in it.
///
/// \param ag  System topology
/// \param ps  System coordinates
RestraintApparatus assembleRestraints(const AtomGraph *ag, const PhaseSpace &ps);
