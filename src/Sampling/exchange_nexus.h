// -*-c++-*-
#ifndef STORMM_EXCHANGE_NEXUS_H
#define STORMM_EXCHANGE_NEXUS_H

#include "copyright.h"
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

using chemistry::PhaseSpaceSynthesis;
using energy::ScoreCard;
using energy::StateVariable;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using topology::AtomGraph;
using topology::UnitCellType;
using trajectory::CoordinateCycle;
using trajectory::CoordinateFileKind;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateFrameWriter;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
using trajectory::CoordinateSeriesWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceReader;
using trajectory::PhaseSpaceWriter;
using trajectory::TrajectoryKind;

/// \brief
class ExchangeNexus {
public:

  /// \brief 
  ExchangeNexus(int system_count_in, int total_swap_count,
                  const std::string &remd_type_in, int frequency_swaps_in, 
                  const std::string &swap_store_in, const std::string &temperature_dist_in,
                  double exchange_probability_in, double tolerance_in, int max_replicas_in, 
                  double initial_temperature_in, double equilibrium_temperature_in,
                  const PhaseSpaceSynthesis *ps_in, const AtomGraphSynthesis *ag_in,
                  const ScoreCard *sc_in);
    
  /// \brief Copy and Move Constructors
  ExchangeNexus(const ExchangeNexus &original) = default;
  ExchangeNexus(ExchangeNexus &&other) = default;
  
  /// \brief Function to return the potential energies of all particles in the system
  ///         in a vector.
  std::vector<double> getPotentialEnergies();
  
  /// \brief Function to return the kinetic energies of all particles in the system
  ///         in a vector.
  std::vector<double> getKineticEnergies();
  
  /// \brief Function to calculate the Hamiltonian of each particle in the current system
  ///
  /// \param  kinetic_energy    Kinetic energies of all particles in the simulation.
  /// \param  potential_energy  Potential energies of all particles in the simulation
  std::vector<double> getHamiltonian(std::vector<double> kinetic_energy, 
                                     std::vector<double> potential_energy);
  
  
  /// \brief Function to return the temperature distributions by invoking one of
  ///         temp_distributions.h and passing the correct algorithm in.
  ///
  ///  \param  initial_temperature      The temperature at the bottom of the range in which a REMD
  ///                                   simulation has to be conducted in.
  ///  \param  equilibrium_temperature  The temperature at the top of the range in which a REMD
  ///                                   simulation has to be conducted in.
  ///  \param  temperature_dist         The algorithm to follow when calculating a temperature
  ///                                   distribution.
  ///  \param  exchange_probability     The user input probability at which to perform a swap.
  ///  \param  cur_ag                   The AtomGraph system for the current REMD.
  std::vector<double> getTempDistribution(const double initial_temperature, 
                                          const double equilibrium_temperature,
                                          const std::string &temperature_dist,
                                          const double exchange_probability,
                                          const AtomGraphSynthesis &cur_ag);
  
  /// \brief The main function to get a REMD input from namelist and invoke the process
  std::vector<double> initiateRemd();
  
private:  
  int system_count;
  int total_swaps;
  std::string remd_type;
  int frequency_swaps;
  int max_replicas;
  std::string swap_store;
  std::string temperature_dist;
  double initial_temperature;
  double equilibrium_temperature;
  double exchange_probability;
  double tolerance;
  PhaseSpaceSynthesis *ps;
  AtomGraphSynthesis *ag;
  ScoreCard *sc;
};

/// \brief A class to conduct Replica Exchange Molecular Dynamics; taking in constraints
///        across multiple simulations. Takes in data from the current running synthesis
///        and user inputs regarding temperatures, probabilities, and temperature distribution
///        algorithms.

} // namespace sampling
} // namespace stormm

#endif
