// -*-c++-*-
#ifndef STORMM_TEST_AMBER_PRMTOP_H
#define STORMM_TEST_AMBER_PRMTOP_H

#include <vector>
#include "../../src/Constants/scaling.h"
#include "../../src/Math/summation.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Topology/atomgraph_analysis.h"
#include "../../src/Topology/atomgraph_enumerators.h"

using namespace stormm::topology;

//-------------------------------------------------------------------------------------------------
// Enumerated real-valued analyses of every topology in a list
//-------------------------------------------------------------------------------------------------
enum class RealInfoCode {
  FIRST_MOLECULE_MASS,     // Compute the average mass of all molecules
  TOTAL_CHARGE,            // Compute the total charge
  MAX_ABS_CHARGE,          // Find the maximum absolute charge for each system
  AVERAGE_LJ_SIGMA,        // Compute the average Lennard-Jones self-interaction sigma value for
	                   //   all atoms
  AVERAGE_LJ_EPS,          // Compute the average Lennard-Jones self-interaction sigma value for
                           //   all atoms
  AVERAGE_BOND_STIFFNESS,  // Compute the average of all bond stiffnesses, weighted by the number
                           //   of such bonds in the system as a whole
  AVERAGE_ANGL_THETA,      // Compute the average of all equilibrium angles (in radians), weighted
                           //   by the number of such bonds in the system as a whole
  AVERAGE_DIHE_COMP,       // Compute the average of all dihedral amplitudes multiplied by their
                           //   respective periodicities
  AVERAGE_UBRD_COMP,       // Compute the average of all Urey-Bradley stiffnesses, times their
                           //   equilibria
  AVERAGE_CIMP_COMP,       // Compute the average of all CHARM improper stiffnesses, plus their
	                   //   phase angles (in radians)
  AVERAGE_CMAP_VALUE       // Compute the average of all grid point values for all CMAP surfaces
};

//-------------------------------------------------------------------------------------------------
// Create real-valued vectors for a particular property based on inputs from a series of
// topologies.
//
// Arguments:
//   topols:         The list of topologies
//   all_top_exist:  Flag to indicate that the topologies do indeed exist and should therefore be
//                   analyzed
//   analysis_code:  Codification of the analysis to perform (see the RealInfoCode enumerator)
//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<double> realAGProp(const std::vector<AtomGraph*> &topols, const bool all_top_exist,
                               const RealInfoCode analysis_code);

#include "test_amber_prmtop.tpp"

#endif

