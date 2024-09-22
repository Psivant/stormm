// -*-c++-*-
#include <cmath>
#include "copyright.h"
#include "exchange_nexus.h"
#include "replica_probability.h"
#include "temp_distributions.h"
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
ExchangeNexus::ExchangeNexus( const int system_count_in, const int total_swap_count,
                              const std::string &remd_type_in, const int frequency_swaps_in, 
                              const std::string &swap_store_in,
                              const std::string &temperature_dist_in,
                              const double exchange_probability_in, const double tolerance_in,
                              const int max_replicas_in, 
                              const double initial_temperature_in,
                              const double equilibrium_temperature_in,
                              const PhaseSpaceSynthesis *ps_in, const AtomGraphSynthesis *ag_in,
                              const ScoreCard *sc_in) :
    system_count{system_count_in}, total_swaps{total_swap_count},
    remd_type{remd_type_in}, frequency_swaps{frequency_swaps_in},
    swap_store{swap_store_in}, temperature_dist{temperature_dist_in},
    exchange_probability{exchange_probability_in}, tolerance{tolerance_in}, 
    max_replicas{max_replicas_in},
    initial_temperature{initial_temperature_in}, 
    equilibrium_temperature{equilibrium_temperature_in},
    ps{const_cast<PhaseSpaceSynthesis*>(ps_in)},
    ag{const_cast<AtomGraphSynthesis*>(ag_in)},
    sc{const_cast<ScoreCard*>(sc_in)}
{}

//-------------------------------------------------------------------------------------------------
std::vector<double> ExchangeNexus::getPotentialEnergies(){
  return sc->reportInstantaneousStates(StateVariable::POTENTIAL_ENERGY);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ExchangeNexus::getKineticEnergies(){
  return sc->reportInstantaneousStates(StateVariable::KINETIC);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ExchangeNexus::getHamiltonian(std::vector<double> kinetic_energy, 
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
std::vector<double> ExchangeNexus::getTempDistribution(const double initial_temperature, 
                                                         const double equilibrium_temperature,
                                                         const std::string &temperature_dist,
                                                         const double exchange_probability,
                                                         const AtomGraphSynthesis &cur_ag) {

  // TODO: Implement different temperature algorithms and return a vector with temp distribution
  std::vector<double> temperatures;
  return temperatures;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ExchangeNexus::initiateRemd() {
    std::vector<double> temperatures = getTempDistribution(initial_temperature, 
                                                           equilibrium_temperature,
                                                           temperature_dist, exchange_probability,
                                                           *ag);
    int num_replicas = temperatures.size();
}

} // namespace sampling
} // namespace stormm
