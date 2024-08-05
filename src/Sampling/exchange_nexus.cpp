// -*-c++-*-
#ifndef STORMM_EXCHANGE_NEXUS_H
#define STORMM_EXCHANGE_NEXUS_H

#include <cmath>
#include "copyright.h"
#include "exchange_nexus.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Namelists/nml_remd.h"
#include "Namelists/nml_dynamics.h"
#include "Potential/scorecard.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/systemcache.h"
#include "Synthesis/synthesis_cache_map.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinate_series.h"

namespace stormm {
namespace sampling {

using symbols::boltzmann_constant;
using synthesis::Condensate;
using synthesis::PhaseSpaceSynthesis;
using synthesis::StructureSource;
using synthesis::SynthesisCacheMap;
using synthesis::SystemGrouping;
using namelist::RemdControls;
using namelist::DynamicsControls;
using topology::AtomGraph;
using trajectory::CoordinateSeries;

//-------------------------------------------------------------------------------------------------
RemdConstraints::RemdConstraints(int system_count_in, int total_swap_count,
                           std::string remd_type_in, int frequency_swaps_in, 
                           std::string swap_store_in, std::string temperature_dist_in,
                           double exchange_probability_in, int max_replicas_in, 
                           double initial_temperature_in, double equilibrium_temperature_in,
                           const PhaseSpaceSynthesis *ps_in, const AtomGraphSynthesis *ag_in,
                           const ScoreCard *sc_in) :
  system_count{system_count_in}, total_swaps{total_swap_count},
  remd_type{remd_type_in}, frequency_swaps{frequency_swaps_in},
  swap_store{swap_store_in}, temperature_dist{temperature_dist_in},
  exchange_probability{exchange_probability_in}, max_replicas{max_replicas_in},
  initial_temperature{initial_temperature_in}, 
  equilibrium_temperature{equilibrium_temperature_in},
  PsSynthesis{*ps_in}, AgSynthesis{*ag_in}, sc{*sc_in}, ps{*ps_in}, ag{*ag_in}
{}

//-------------------------------------------------------------------------------------------------
std::vector<double> RemdConstraints::getPotentialEnergies(){
  return sc.reportInstantaneousStates(StateVariable::POTENTIAL_ENERGY);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> RemdConstraints::getKineticEnergies(){
  return sc.reportInstantaneousStates(StateVariable::KINETIC);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> RemdConstraints::getHamiltonian(std::vector<double> kinetic_energy, 
                                                 std::vector<double> potential_energy){
    std::vector<double> hamiltonian;
    if (kinetic_energy.size() >= potential_energy.size()) {
        for (size_t i = 0; i < kinetic_energy.size(); ++i) {
            hamiltonian.push_back(kinetic_energy[i] + potential_energy[i]);
        }
    }
    // TODO: Else condition, as well as making sure vector indices refer to same particle
    return hamiltonian;
}

//-------------------------------------------------------------------------------------------------
double RemdConstraints::getTotalHamiltonian(std::vector<double> hamiltonian){
    double totalHamiltonian = 0;
    for (double h : hamiltonian) {
        totalHamiltonian += h;
    }
    return totalHamiltonian;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> RemdConstraints::getTempDistribution(double initial_temperature, 
																												 double equilibrium_temperature,
                                     std::string temperature_dist, double exchange_probability,
                                     const AtomGraphSynthesis cur_ag){
    // TODO: Implement different temperature algorithms and return a vector with temp distribution
    
    std::vector<double> temperatures;
    return temperatures;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> RemdConstraints::getProbabilityAtState(std::vector<double> temperatures, 
                                                        std::vector<double> hamiltonian){
    std::vector<double> probabilities;
    for (size_t i = 0; i < temperatures.size(); ++i) {
        double totalHamiltonian = getTotalHamiltonian(hamiltonian);
        double beta = 1 / (temperatures[i] * boltzmann_constant);
        probabilities.push_back(exp((-beta) * totalHamiltonian));
    }
    return probabilities;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> RemdConstraints::RemdMain(){ 
    std::vector<double> temperatures = getTempDistribution(initial_temperature, 
                                                equilibrium_temperature, temperature_dist,
                                                exchange_probability, ag);
    // According to my knowledge, we have the temperature distribution at this point
    // (please ignore the incomplete code for now)
    // For each replica in this distribution, we have to create a "Synthesis" instance
    // for AtomGraph and PhaseSpace
    // Then in each replica we just go odd/even, calculate probabilities
    // If they're higher for a given temperature, we do the swap. 
    // I have attached a slack message here to explain my roadblock.
    return temperatures;
}
}// namespace sampling
} // namespace stormm

#endif