// -*-c++-*-

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

struct RemdConstraints {
	RemdConstraints(int system_count_in, int total_swap_count,
                  std::string remd_type_in, int frequency_swaps_in, 
                  std::string swap_store_in, std::string temperature_dist_in,
                  double exchange_probability_in, int max_replicas_in, 
                  double initial_temperature_in, double equilibrium_temperature_in,
                  const PhaseSpaceSynthesis *ps_in, const AtomGraphSynthesis *ag_in,
                  const ScoreCard *sc_in);
    
  /// \brief Copy and Move Constructors
  RemdConstraints(const RemdConstraints &original) = default;
  RemdConstraints(RemdConstraints &&other) = default;
  
  /// \brief Function to return the potential energies of all particles in the system
  ///				 in a vector.
  std::vector<double> getPotentialEnergies();
  
  /// \brief Function to return the kinetic energies of all particles in the system
  ///				 in a vector.
  std::vector<double> getKineticEnergies();
  
  /// \brief Function to calculate the Hamiltonian of each particle in the current system
  ///
  ///	\param	kinetic_energy: 	Kinetic energies of all particles in the simulation.
  ///	\param	potential_energy: Potential energies of all particles in the simulation
  std::vector<double> getHamiltonian(std::vector<double> kinetic_energy, 
                                     std::vector<double> potential_energy);
  
  /// \brief Function to calculate the Hamiltonian of all particles in the current system
  ///
  /// \param hamiltonian:	Vector of the hamiltonians of each particle in the system.
  double getTotalHamiltonian(std::vector<double> hamiltonian);
  
  /// \brief Function to return the temperature distributions by invoking one of
  ///				 temp_distributions.h and passing the correct algorithm in.
  ///
  ///	\param	initial_temperature:			The temperature at the bottom of the range in which a REMD
  ///																		simulation has to be conducted in.
  ///	\param	equilibrium_temperature:	The temperature at the top of the range in which a REMD
  ///																		simulation has to be conducted in.
  ///	\param	temperature_dist:					The algorithm to follow when calculating a temperature
  ///																		distribution.
  ///	\param	exchange_probability:			The user input probability at which to perform a swap.
  ///	\param	cur_ag:										The AtomGraph system for the current REMD.
  std::vector<double> getTempDistribution(double initial_temperature, 
																					double equilibrium_temperature,
                                     			std::string temperature_dist, double exchange_probability,
                                     			const AtomGraphSynthesis cur_ag);
                      										
  ///	\brief:	Function to check the probability at a given state between two adjacent replicas.
  ///					If this is greater than the exchange_probability, a swap occurs.
  ///
  ///	\param	temperatures:	The temperatures of all the adjacent replicas.
  ///	\param	hamiltonian:	The hamiltonians of all the particles in the current REMD simulation.
  std::vector<double> getProbabilityAtState(std::vector<double> temperatures, 
                                            std::vector<double> hamiltonian);
  
  /// \brief The main function to get a REMD input from namelist and invoke the process
  std::vector<double> RemdMain();
  
  private:  
    const PhaseSpaceSynthesis PsSynthesis;
    const AtomGraphSynthesis AgSynthesis;
    const ScoreCard sc;
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
    const PhaseSpaceSynthesis ps;
    const AtomGraphSynthesis ag;
};

/// \brief A class to conduct Replica Exchange Molecular Dynamics; taking in constraints
///        across multiple simulations. Takes in data from the current running synthesis
//         and user inputs regarding temperatures, probabilities, and temperature distribution
//         algorithms.

}
}
