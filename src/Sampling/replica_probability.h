// -*-c++-*-
#ifndef STORMM_REPLICA_PROBABILITY_H
#define STORMM_REPLICA_PROBABILITY_H

#include <cmath>
#include <vector>
#include "copyright.h"
#include "Constants/symbol_values.h"

namespace stormm {
namespace sampling {

using symbols::boltzmann_constant;

/// \brief  Function to calculate the Hamiltonian of all particles in the current system
///
/// \param hamiltonian  Vector containing all the particle hamiltoninas in the REMD simulation
double getTotalHamiltonian(std::vector<double> hamiltonian);

/// \brief  The following function returns the probabilities of a MD at every replica during a
///         REMD simulation. Currently the Math has been coded for temperature REMD.
///
/// \param temperatures Vector containing all the temperatures in the current REMD simulation
/// \param hamiltonian  Vector containing all the particle hamiltonians in the REMD simulation
std::vector<double> getProbabilityAtState(std::vector<double> &temperatures, 
                                            std::vector<double> &hamiltonian);

} // namespace sampling
} // namespace stormm

#endif
