// -*-c++-*-
#ifndef STORMM_TICKCOUNTER_H
#define STORMM_TICKCOUNTER_H

#include <string>
#include <vector>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "Reporting/error_format.h"
#include "Structure/local_arrangement.h"
#include "Structure/structure_enumerators.h"
#include "multiplication.h"
#include "summation.h"

namespace stormm {
namespace stmath {

using structure::imageValue;
using structure::ImagingMethod;

/// \brief Make an array of integral settings and tick it forward (or backward).  The analogy is a
///        series of wheels A - N, each with n_A, n_B, ... , n_N settings.  With each tick, the
///        lowest wheel (position 0 of the array) advances one tick.  When each wheel completes one
///        revolution, it triggers an advance of the next wheel.  This collapses to a simple base-K
///        counting system for a series of wheels with equal numbers of K settings apiece.
template <typename T> class TickCounter {
public:

  /// \brief The constructor takes a series of state counts for each "wheel" in the counter.  The
  ///        initial state can also be provided.
  /// \{
  explicit TickCounter(const std::vector<int> &state_limits_in,
                       const std::vector<int> &settings_in = std::vector<int>());

  explicit TickCounter(const std::vector<std::vector<T>> &state_values_in,
                       const std::vector<int> &settings_in = std::vector<int>());
  /// \}

  /// \brief The copy and move constructors as well as assignment operators are all valid.
  /// \{
  TickCounter(const TickCounter &original) = default;
  TickCounter(TickCounter &&original) = default;
  TickCounter& operator=(const TickCounter &original) = default;
  TickCounter& operator=(TickCounter &&original) = default;
  /// \}

  /// \brief Get the number of independent settings, the number of variables.
  int getVariableCount() const;

  /// \brief Get the logarithm of the number of overall permutations.
  double getLogPermutations() const;
  
  /// \brief Get a const reference to the vector of current settings.
  const std::vector<int>& getSettings() const;

  /// \brief Get the number of available states.
  ///
  /// Overloaded:
  ///   - Get the degeneracy of states for a particular variable
  ///   - Get a const reference to the vector of numbers of options for each state.
  ///
  /// \param var_index  The variable of interest
  /// \{
  int getStateLimits(int var_index) const;
  const std::vector<int>& getStateLimits() const;
  /// \}

  /// \brief Get the state represented by a particular combination of settings, drawing on the
  ///        internal representations of each state.
  ///
  /// Overloaded:
  ///   - Get the current state of one of the independent variables.
  ///   - Get the current state of all variables.
  ///
  /// \param var_index  Index of the variable of interest
  /// \{
  T getState(int var_index) const;
  std::vector<T> getState() const;
  /// \}

  /// \brief Get a vector of all possible states for one of the variables.
  ///
  /// \param var_index  Index of the variable of interest
  /// \{
  std::vector<T> getPossibleStates(int var_index) const;
  /// \}

  /// \brief Get the exact number of permutations for this object.
  llint getExactPermutationCount() const;

  /// \brief Get the number of permutations count for this object as a floating-point number.
  double getApproximatePermutationCount() const;

  /// \brief Get the natural logarithm of the number of permutations of this object.
  double getLogPermutationCount() const;
  
  /// \brief Advance the counter by a specified number of steps.
  void advance(int steps = 1);

  /// \brief Reverse the counter by a specified number of steps.
  void reverse(int steps = 1);
  
  /// \brief Randomize some or all variables of the object.
  ///
  /// Overloaded:
  ///   - Randomize all counter variables
  ///   - Randomize selected counter variables
  ///
  /// \param xrs        The random number generator (the type of this argument is the basis of the
  ///                   templating)
  /// \param apply_rng  A vector containing entries for all counter variables and marked TRUE to
  ///                   cause the variable to be randomly reassigned, or a vector containing a
  ///                   limited number of counter indices
  /// \{
  template <typename Trng> void randomize(Trng *xrs);
  template <typename Trng> void randomize(Trng *xrs, const std::vector<bool> &apply_rng);
  template <typename Trng> void randomize(Trng *xrs, const std::vector<int> &apply_rng);
  /// \}

  /// \brief Set the object's state to an arbitrary series of values.
  ///
  /// Overloaded:
  ///   - Set all variables
  ///   - Set a specific variable
  ///
  /// \param new_state  The value or values to set
  /// \param var_index  Index of a particular counter variable to set
  /// \{
  void set(const std::vector<int> &new_state);

  void set(int new_state, int var_index);
  /// \}
  
  /// \brief Reset all states of the counter to zero.
  void reset();

  /// \brief Load the available states for one or more independent variables.
  ///
  /// Overloaded:
  ///   - Load the available states for a particular variable.
  ///   - Load the available states for all variables.
  ///
  /// \param state_values_in  The enumerated states that one or all independent variables can take
  /// \param var_index        The specific variable of interest
  /// \{
  void loadStateValues(const std::vector<T> &state_values_in, int var_index);
  void loadStateValues(const std::vector<std::vector<T>> &state_values_in);
  /// \}

private:
  int variable_count;             ///< The number of independent variables in the counter
  double log_permutations;        ///< Logarithm of the total number of permutations
  std::vector<int> settings;      ///< "Fingerprint" of the counter, i.e. a sequence such as
                                  ///<   { 0, 1, 0, 6, 0, 4, 0, 2 } for a counter with
                                  ///<   { 3, 5, 1, 7, 4, 5, 5, 8 } unique settings for each state
  std::vector<int> state_limits;  ///< The number of unique settings for each state, thus the
                                  ///<   maximum values of entries in the settings vector.
  std::vector<T> state_values;    ///< List of values for each state.  By default, the values will
                                  ///<   be uninitialized. (It is very difficult to write code to
                                  ///<   automatically assign values to data that could be of any
                                  ///<   type, including structs or tuples with attributes that
                                  ///<   cannot be generalized.)

  /// Bounds aray for state_values above
  std::vector<int> state_value_bounds;

  /// \brief Check that the index of a given variable is within the number of independent variables
  ///        for the object.
  ///
  /// \param var_index  The specific variable of interest
  void variableBoundsCheck(const int var_index) const;
  
  /// \brief Validate the settings to ensure that they are within the stated limits and above zero
  ///        in all counter variables.
  void validateSettings() const;
};

/// \brief Load states appropriate for an even sampling of a specific range in each variable of a
///        TickCounter object.  These explicitly typed functions exist to serve the most common
///        types of TickCounter objects.
///
/// Overloaded:
///   - Accept double- or single-precision floating point limits the ranges
///   - Impose periodic boundary conditions on the ranges
///
/// \param tc                   The object to load with state values
/// \param ranges               Ranges for the state values (the values will be loaded with an
///                             even spread over the stated range).  The lower limit of each range
///                             is held in the "x" member of the tuple, while the upper limit is
///                             held in the "y" member.
/// \param periodic_boundaries  Limits of periodic boundaries applicable to each state
/// \{
template <typename T, typename T2>
void loadScalarStateValues(TickCounter<T> *tc, const std::vector<T2> &ranges,
                           const std::vector<T> &periodic_boundaries);

template <typename T, typename T2>
void loadScalarStateValues(TickCounter<T> *tc, const std::vector<T2> &ranges);
/// \}

} // namespace stmath
} // namespace stormm

#include "tickcounter.tpp"

#endif
