// -*-c++-*-
#ifndef STORMM_NML_REMD_H
#define STORMM_NML_REMD_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/symbol_values.h"
#include "Parsing/textfile.h"
#include "Potential/energy_enumerators.h"
#include "Structure/structure_enumerators.h"
#include "Trajectory/trajectory_enumerators.h"
#include "input.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using constants:: CartesianDimension;
using constants::PrecisionModel;
using constants::ExceptionResponse;
using parse::TextFile;
using parse::WrapTextSearch;
using structure::ApplyConstraints;
using namelist::NamelistEmulator;

constexpr int default_total_swaps = 0;
constexpr char default_remd_type[] = "Temperature";
constexpr int default_freq_swaps = 500;
constexpr char default_swap_store[] = "Successful";
constexpr char default_temp_distribution[] = "Van Der Spoel";
constexpr double default_exchange_probability = 0.2;
constexpr double default_tolerance = 0.00001;
constexpr int default_max_replicas = 1000;
constexpr double default_low_temperature = 298.15;
constexpr double default_high_temperature = 398.15;
 
/// \brief Encapsulating the data extracted from the REMD namelist
class RemdControls {
public:

  /// \brief The Constructors can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &REMD namelist
  /// \param found_nml   Indication of whether the namelist was found in the input file
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &conformer namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  RemdControls(ExceptionResponse policy_in = ExceptionResponse::DIE);

  RemdControls(const TextFile &tf, int *start_line, bool *found_nml,
               ExceptionResponse policy_in = ExceptionResponse::DIE,
               WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief Copy and Move Constructors, inspired from nml_mesh.h
  /// \{
  RemdControls(const RemdControls &original) = default;
  RemdControls(RemdControls &&original) = default;
  RemdControls& operator=(const RemdControls &original) = default;
  RemdControls& operator=(RemdControls &&original) = default;
  /// \}
    
  /// \brief Getters for the REMD Namelist
  /// \brief Get the total number of swaps to be attempted in this simulation
  int getTotalSwapCount() const;

  /// \brief Get the type of REMD to be performed in this simulation
  std::string getRemdType() const;

  /// \brief Get the frequency of steps to perform before a swap is attempted
  int getFrequencyOfSwaps() const;

  /// \brief Get the string encoding for the method of storage to use for swaps
  /// TODO: Storage methods need to be revisited upon creating ExchangeNexus
  std::string getSwapStore() const;
  
  /// \brief Get the string encoding for the type of temperature distribution in use
  /// TODO: This needs to be revisited upon coding out ExchangeNexus
  std::string getTemperatureDistributionMethod() const;

  /// \brief Get the user input Exchange Probability between 2 adjacent replicas
  double getExchangeProbability();

  /// \brief Get the tolerance set between desired exchange probability and predicted probability
  double getTolerance();

  /// \brief Get the user input for Maximum Replicas allowed
  int getMaxReplicas() const;
  
  /// \brief Get the vector of initial temperature targets for all groups of atoms and systems.
  // TODO: Change this to just a double
  const double getInitialTemperature();

  /// \brief Get the vector of final temperature targets for all groups of atoms and systems.
  const double getEquilibriumTemperature();
  
  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator& getTranscript() const;

  /// \brief Set the total number of swaps to be performed in a REMD simulation
  /// 
  /// \param total_swap_count_in The number of total swaps to take
  void setTotalSwapCount(int total_swap_count_in);
  
  /// \brief Set the REMD Type: 1) Temperature, or 2) Hamiltonian
  ///
  /// \param remd_type_in String encoding for the type of REMD to perform
  void setRemdType(std::string remd_type_in);

  /// \brief Set the frequency of steps to take before performing a swap
  ///
  /// \param frequency_swaps_in The number of steps before a swap
  void setFrequencyOfSwaps(int frequency_swaps_in);

  /// \brief Set the amount of data to be stored for swap history
  ///
  /// \param swap_store_in The amount of data to be stored
  void setSwapStore(std::string swap_store_in);

  /// \brief Set the algorithm to be used to calculate the temperature distribution
  ///
  /// \param temperature_dist_in String encoding for the type of algorithm to use to calculate
  /// the temperature distribution
  void setTemperatureDistributionMethod(std::string temperature_dist_in);

  /// \brief Set the exchange probability between two adjacent replicas in REMD
  ///
  /// \param exchange_probability_in New exchange probability for adjacent replicas
  void setExchangeProbability(double exchange_probability_in);

  /// \brief Set the tolerance between desired probability and poredicted probability
  ///
  /// \param tolerance_in Tolerance
  void setTolerance(double tolerance_in);

  /// \brief Set the maximum number of replicas the REMD algorithm is allowed to create
  ///
  /// \param max_replicas_in Number of replicas
  void setMaxReplicas(int max_replicas_in);

  /// \brief Set the initial temperature for the REMD simulation
  /// 
  /// \param low_temperature_in Initial temperature
  void setInitialTemperature(double low_temperature_in);

  /// \brief Set the equilibrium temperature for the REMD simulation
  ///
  /// \param high_temperature_in The equilibrium temperature
  void setEquilibriumTemperature(double high_temperature_in);

private:
  ExceptionResponse policy;      ///< The course to take when encountering bad input
  int total_swap_count;          ///< The number of total times to attempt for swaps
  std::string remd_type;         ///< The type of REMD to perform: Temperature or Hamiltonian
  int frequency_swaps;           ///< The number of steps to perform before a swap is attempted
  std::string swap_store;        ///< String encoding for how much data the user wants to store
                                 ///<   for a swap: [None, Successful, State]
  std::string temperature_dist;  ///< String encoding for the type of distribution to use
                                 ///<   for setting up adjacent replicas
  double exchange_probability;   ///< The probability for a successful swap between two adjacent
                                 ///<   replicas in a REMD simulation
  double tolerance;              ///< Tolerance for convergence criterion between
                                 ///  exchange_probability and the predicted probability.
  int max_replicas;              ///< The maximum number of replicas in a given REMD simulation

  /// A series of the initial target temperatures for various thermostat partitions of the 
  /// synthesis of systems
  double low_temperature;

  /// A series of the final target temperatures for the various thermostat partitions of the
  /// synthesis of systems
  double high_temperature;

  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;
    
  /// Validator functions for user input.
  /// \brief Validate the total number of swaps
  void validateSwapCount();
    
  /// \brief Validate the Kind of REMD (Hamiltonian or Temperature)
  void validateRemdKind();

  /// \brief  Validate which of either to use: total swaps or frequencxy of swaps
  void validateFreqOrTotalSwaps();

  /// \brief Validate the frequency of swaps in the MD
  void validateSwapFrequency();

  /// \brief Validate the amount of data being stored for swaps
  void validateSwapStore();
 
  /// \brief Validate a temperature requested for thermoregulation.
  ///
  /// \param t  The temperature in question
  void validateTemperature();

  /// \brief Validate the Temperature Distribution Algorithm in use
  void validateTemperatureDistribution();

  /// \brief Validate the Exchange Probability = [0, 1]
  void validateExchangeProbability();

  /// \brief Validate the tolerance between desired exchange probability and predicted probability
  void validateTolerance();

  /// \brief Validate Max Replicas > 0
  void validateMaxReplicas();
};

/// \brief Free function to read the &remd namelist
///
/// \param tf          Text of file containing the inputs, read into RAM
/// \param start_line  Line of the input file at which to begin the scan
/// \param found       Indicator that the namelist was found in the input file
/// \param policy      Response to bad inputs
/// \param wrap        Indicate that the search for a &conformer namelist should carry on from the
///                    beginning of an input file if no such namelist is found from the original 
///                    starting point
NamelistEmulator remdInput(const TextFile &tf, int *start_line, bool *found,
                           ExceptionResponse policy = ExceptionResponse::DIE,
                           WrapTextSearch wrap = WrapTextSearch::NO);

} // namespace namelist
} // namespace stormm

#endif
