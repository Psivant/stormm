#include "copyright.h"
#include "Constants/symbol_values.h"
#include "Math/matrix_ops.h"
#include "Math/math_enumerators.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Reporting/error_format.h"
#include "nml_remd.h"

namespace stormm {
namespace namelist {

using constants::getEnumerationName;
using parse::minimalRealFormat;
using parse::NumberFormat;
using parse::realToString;
using parse::strcmpCased;
using stmath::computeBoxTransform;

//-------------------------------------------------------------------------------------------------
RemdControls::RemdControls(const ExceptionResponse policy_in) :
    policy{policy_in},
    total_swap_count{default_total_swaps},
    remd_type{std::string(default_remd_type)},
    frequency_swaps{default_freq_swaps},
    swap_store{std::string(default_swap_store)},
    temperature_dist{std::string(default_temp_distribution)},
    exchange_probability{default_exchange_probability},
    max_replicas{default_max_replicas},
    low_temperature{default_low_temperature},
    high_temperature{default_high_temperature},
    nml_transcript{"remd"}
{}

//-------------------------------------------------------------------------------------------------
RemdControls::RemdControls(const TextFile&tf, int *start_line, bool *found_nml,
                           const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    RemdControls(policy_in)
{
  const NamelistEmulator t_nml = remdInput(tf, start_line, found_nml, policy, wrap);
  nml_transcript = t_nml;
  t_nml.assignVariable(&total_swap_count, "total_swaps");
  t_nml.assignVariable(&remd_type, "remd_type");
  t_nml.assignVariable(&frequency_swaps, "freq_swaps");
  t_nml.assignVariable(&swap_store, "swap_store");
  t_nml.assignVariable(&temperature_dist, "temp_distribution");
  t_nml.assignVariable(&exchange_probability, "exchange_probability");
  t_nml.assignVariable(&tolerance, "tolerance");
  t_nml.assignVariable(&max_replicas, "max_replicas");
  t_nml.assignVariable(&low_temperature, "low_temperature");
  t_nml.assignVariable(&high_temperature, "high_temperature");

  // Validate Input
  validateSwapCount();
  validateRemdKind();
  validateFreqOrTotalSwaps();
  validateSwapFrequency();
  validateSwapStore();
  validateTemperature();
  validateTemperatureDistribution();
  validateExchangeProbability();
  validateTolerance();
  validateMaxReplicas();
}

//-------------------------------------------------------------------------------------------------
int RemdControls::getTotalSwapCount() const{
  return total_swap_count;
}

//-------------------------------------------------------------------------------------------------
std::string RemdControls::getRemdType() const{
  return remd_type;
}

//-------------------------------------------------------------------------------------------------
int RemdControls::getFrequencyOfSwaps() const {
  return frequency_swaps;
}

//-------------------------------------------------------------------------------------------------
std::string RemdControls::getSwapStore() const {
  return swap_store;
}

//-------------------------------------------------------------------------------------------------
std::string RemdControls::getTemperatureDistributionMethod() const{
  return temperature_dist;
}

//-------------------------------------------------------------------------------------------------
double RemdControls::getExchangeProbability(){
  return exchange_probability;
}

//-------------------------------------------------------------------------------------------------
double RemdControls::getTolerance() {
  return tolerance;
}

//-------------------------------------------------------------------------------------------------
int RemdControls::getMaxReplicas() const{
  return max_replicas;
}
//-------------------------------------------------------------------------------------------------
const double RemdControls::getInitialTemperature() {
  return low_temperature;
}

//-------------------------------------------------------------------------------------------------
const double RemdControls::getEquilibriumTemperature() {
  return high_temperature;
}

//-------------------------------------------------------------------------------------------------
const NamelistEmulator& RemdControls::getTranscript() const {
  return nml_transcript;
}

//-------------------------------------------------------------------------------------------------
void RemdControls::setTotalSwapCount(int total_swap_count_in) {
  total_swap_count = total_swap_count_in;
}

//-------------------------------------------------------------------------------------------------
void RemdControls::setRemdType(std::string remd_type_in) {
  remd_type = remd_type_in;
}

//-------------------------------------------------------------------------------------------------
void RemdControls::setFrequencyOfSwaps(int frequency_swaps_in) {
  frequency_swaps = frequency_swaps_in;
}

//-------------------------------------------------------------------------------------------------
void RemdControls::setSwapStore(std::string swap_store_in) {
  swap_store = swap_store_in;
}

//-------------------------------------------------------------------------------------------------
void RemdControls::setTemperatureDistributionMethod(std::string temperature_dist_in) {
  temperature_dist = temperature_dist_in;
}

//-------------------------------------------------------------------------------------------------
void RemdControls::setExchangeProbability(double exchange_probability_in){
  exchange_probability = exchange_probability_in;
}

//-------------------------------------------------------------------------------------------------
void RemdControls::setTolerance(double tolerance_in) {
  tolerance = tolerance_in;
}

//-------------------------------------------------------------------------------------------------
void RemdControls::setMaxReplicas(int max_replicas_in) {
  max_replicas = max_replicas_in;
}

//-------------------------------------------------------------------------------------------------
void RemdControls::setInitialTemperature(double low_temperature_in) {
  low_temperature = low_temperature_in;
}

//-------------------------------------------------------------------------------------------------
void RemdControls::setEquilibriumTemperature(double high_temperature_in) {
  high_temperature = high_temperature_in;
}

//-------------------------------------------------------------------------------------------------
void RemdControls::validateSwapCount(){
  if (total_swap_count < 0) {
    switch(policy) {
    case ExceptionResponse::DIE:
      rtErr("The total number of swaps to attempt cannot be less than 0", "RemdControls",
            "validateSwapCount");
    case ExceptionResponse::WARN:
      rtWarn("The total number of swaps cannot be less than 0.", "RemdControls",
             "validateSwapCount");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    total_swap_count = default_total_swaps;
  }
}

//-------------------------------------------------------------------------------------------------
void RemdControls::validateRemdKind() {
  if (strcmpCased(remd_type, "temperature", CaseSensitivity::NO) == false &&
      strcmpCased(remd_type, "hamiltonian", CaseSensitivity::NO) == false) {
    switch(policy) {
    case ExceptionResponse::DIE:
      rtErr("REMD Type has to be either Temperature or Hamiltonian", "RemdControls",
            "validateRemdKind");
      break;
    case ExceptionResponse::WARN:
      rtWarn("REMD Type has to be either Temperature or Hamiltonian. Default value " + 
             static_cast<std::string>(default_remd_type) + " will be applied.", "RemdControls",
             "validateRemdKind");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    remd_type = default_remd_type;
  }
}
//-------------------------------------------------------------------------------------------------
void RemdControls::validateFreqOrTotalSwaps() {
 if(nml_transcript.getKeywordStatus("freq_swaps") == InputStatus::DEFAULT && 
    total_swap_count > 0) {
  switch(policy) {
  case ExceptionResponse::DIE:
    rtErr("Using total swap count and frequency of steps for a swap at the same time leads to a "
          "conflict. Please specify either value.", "RemdControls", "validateFreqOrTotalSwaps");
    break;
  case ExceptionResponse::WARN:
    rtWarn("Using total swap count and frequency of steps for a swap at the same time leads to a "
           "conflict. Total swap count takes priority, and will be considered.", "RemdControls", 
           "validateFreqOrTotalSwaps");
    break;
  case ExceptionResponse::SILENT:
    break;
  }
  frequency_swaps = -1;
 }  
}


//-------------------------------------------------------------------------------------------------
void RemdControls::validateSwapFrequency(){
  if (frequency_swaps <= 0) {
    switch(policy) {
      case ExceptionResponse::DIE:
        rtErr("The total steps before a swap cannot be less than 0", "RemdControls",
            "validateSwapFrequency");
      case ExceptionResponse::WARN:
        rtWarn("The total steps before a swap cannot be less than 0.", "RemdControls",
               "validateSwapFrequency");
        break;
      case ExceptionResponse::SILENT:
        break;
    }
    frequency_swaps = default_freq_swaps;
  }
}

//-------------------------------------------------------------------------------------------------
void RemdControls::validateSwapStore() {
  //TODO: Change this based on design decision
  if (strcmpCased(swap_store, "successful", CaseSensitivity::NO) == false &&
      strcmpCased(swap_store, "swaps", CaseSensitivity::NO) == false) {
    switch(policy) {
      case ExceptionResponse::DIE:
        rtErr("Swap Store has to be either Successful or Swaps", "RemdControls",
            "validateSwapStore");
      case ExceptionResponse::WARN:
        rtWarn("Swap Store has to be either Successful or Swaps. Default value " + 
            static_cast<std::string>(default_swap_store) + " will be applied.", 
            "RemdControls", "validateSwapStore");
        break;
      case ExceptionResponse::SILENT:
        rtWarn("Swap Store has to be either Successful or Swaps. Default value " + 
               static_cast<std::string>(default_swap_store) + " will be applied.", 
               "RemdControls", "validateSwapStore");
        break;
    }
    swap_store = default_swap_store;
  }
}

//-------------------------------------------------------------------------------------------------
void RemdControls::validateTemperature(){
  if (low_temperature < 0 || high_temperature < 0 || 
      high_temperature < low_temperature) {
    switch(policy) {
      case ExceptionResponse::DIE:
        rtErr("The temperatures cannot be less than 0", "RemdControls",
            "validateTemperature");
      case ExceptionResponse::WARN:
        rtWarn("The temperatures cannot be less than 0, default values will be applied.", 
               "RemdControls", "validateExchangeProbability");
        break;
      case ExceptionResponse::SILENT:
        break;
    }
    low_temperature = default_low_temperature;
    high_temperature = default_high_temperature;
  }
}

//-------------------------------------------------------------------------------------------------
void RemdControls::validateTemperatureDistribution() {
  //TODO: Change this based on design decision
  if (strcmpCased(temperature_dist, "van der spoel", CaseSensitivity::NO) == false &&
      strcmpCased(temperature_dist, "tiger2", CaseSensitivity::NO) == false) {
    switch(policy) {
    case ExceptionResponse::DIE:
      rtErr("Invalid temperature distribution algorithm detected", "RemdControls",
            "validateTemperatureDistribution");
    case ExceptionResponse::WARN:
      rtWarn("Invalid temperature distribution algorithm detected. Default value " + 
             static_cast<std::string>(default_temp_distribution) + " will be applied.", 
             "RemdControls", "validateTemperatureDistribution");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    temperature_dist = default_temp_distribution;
  }
}

//-------------------------------------------------------------------------------------------------
void RemdControls::validateExchangeProbability() {
  if (exchange_probability < 0 || exchange_probability > 1) {
    switch(policy) {
    case ExceptionResponse::DIE:
      rtErr("The exchange probability should be in range [0, 1]", "RemdControls",
            "validateExchangeProbability");
    case ExceptionResponse::WARN:
      rtWarn("The exchange probability should be in range [0, 1]", "RemdControls",
             "validateExchangeProbability");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    exchange_probability = default_exchange_probability;
  }
}

//-------------------------------------------------------------------------------------------------
void RemdControls::validateTolerance() {
  if(!std::isfinite(tolerance)) {
    switch(policy) {
    case ExceptionResponse::DIE:
      rtErr("Tolerance between desired exchange probability and predicted probability is not a "
            "valid numerical value", "RemdControls", "validateTolerance");
      break;
    case ExceptionResponse::WARN:
      rtWarn("The tolerance between desired exchange probability and predicted probability is "
             "not a valid numerical value", "RemdControls", "validateTolerance");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void RemdControls::validateMaxReplicas() {
  if (max_replicas <= 0) {
    switch(policy) {
    case ExceptionResponse::DIE:
      rtErr("The Maximum Replicas should be greater than 0", "RemdControls",
            "validateMaxReplicas");
    case ExceptionResponse::WARN:
      rtWarn("The Maximum Replicas should be greater than 0", "RemdControls",
             "validateMaxReplicas");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    max_replicas = default_max_replicas;
  }
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator remdInput(const TextFile &tf, int *start_line, bool *found,
                           const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("remd", CaseSensitivity::AUTOMATIC, policy, "Collects details of a "
                         "REMD simulation input (particles, energies, coords, temperature)");

  // User input keywords
  t_nml.addKeyword("total_swaps", NamelistType::INTEGER, std::to_string(default_total_swaps));
  t_nml.addHelp("total_swaps", "The total number of times to attempt for swaps across the "
                "simulation.");
  t_nml.addKeyword("remd_type", NamelistType::STRING, std::string(default_remd_type));
  t_nml.addHelp("remd_type", "String encoding for the type of REMD to be performed.  Possible "
                "values include TEMPERATURE or HAMILTONIAN.");
  t_nml.addKeyword("freq_swaps", NamelistType::INTEGER, std::to_string(default_freq_swaps));
  t_nml.addHelp("freq_swaps", "The number (frequency) of steps to perform before a swap is "
                "attempted.  Specifying a value here will override any value given for "
                "total_swaps, or produce a conflict and an exception depending on the error "
                "tolerance.");
  t_nml.addKeyword("swap_store", NamelistType::STRING, std::string(default_swap_store));
  t_nml.addHelp("swap_store", "The amount of information to be stored regarding checking and "
                "success of swaps.  Options include NONE, SUCCESSFUL, and STATE.");
  t_nml.addKeyword("temp_distribution", NamelistType::STRING,
                   std::string(default_temp_distribution));
  t_nml.addHelp("temp_distribution", "String encoding for the specific algorithm to use for "
                "creating a temperature distribution between an initial temperature and a final "
                "temperature");
  t_nml.addKeyword("exchange_probability", NamelistType::REAL,
                   realToString(default_exchange_probability));
  t_nml.addHelp("exchange_probability", "the likelihood that a swap between adjacent replicas "
                "will be successful. To achieve effective statistical sampling, this probability "
                "should be high enough to maintain strong coupling between replicas "
                "(typically 0.2 - 0.4)");
  t_nml.addKeyword("tolerance", NamelistType::REAL,
                   realToString(default_tolerance));
  t_nml.addHelp("tolerance", "The tolerance between desired exchange probability and actual "
                "calculated probability. A typical value of 10^-4 is good.");
  t_nml.addKeyword("max_replicas", NamelistType::INTEGER, std::to_string(default_max_replicas));
  t_nml.addHelp("max_replicas", "The number of maximum total replicas to be made. Higher numbers "
                "are better to ensure that temperature differences between adjacent replicas are "
                "not too high.");
  t_nml.addKeyword("low_temperature", NamelistType::REAL,
                   realToString(default_low_temperature));
  t_nml.addHelp("low_temperature", "The temperature of the coldest replica in a temperature "
                "REMD simulation, in Kelvin");
  t_nml.addKeyword("high_temperature", NamelistType::REAL,
                   realToString(default_high_temperature));
  t_nml.addHelp("high_temperature", "The temperature of the hottest replica in a temperature REMD "
                "simulation, in Kelvin");
  
  // Search for the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  return t_nml;
}

} // namespace namelist
} // namespace stormm
