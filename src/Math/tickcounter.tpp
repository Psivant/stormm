// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename T>
TickCounter<T>::TickCounter(const std::vector<int> &state_limits_in,
                            const std::vector<int> &settings_in) :
    variable_count{static_cast<int>(state_limits_in.size())},
    log_permutations{logProduct(state_limits_in)},
    settings{settings_in.size() == 0 ? std::vector<int>(state_limits_in.size(), 0) : settings_in},
    state_limits{state_limits_in},
    state_values{}, state_value_bounds{}
{
  if (variable_count == 0) {
    rtErr("No variables were indicated.", "TickCounter");
  }
  int minloc = 0;
  int minval = state_limits[0];
  for (int i = 1; i < variable_count; i++) {
    if (state_limits[i] < minval) {
      minval = state_limits[i];
      minloc = i;
    }
  }
  if (minval <= 0) {
    rtErr("A variable with " + std::to_string(minval) + " possible states is invalid at "
          "position " + std::to_string(minloc) + ".", "TickCounter");
  }
  state_value_bounds.resize(variable_count + 1, 0);
  for (int i = 0; i < variable_count; i++) { 
    state_value_bounds[i] = state_limits[i];
  }
  prefixSumInPlace(&state_value_bounds, PrefixSumType::EXCLUSIVE);
  state_values.resize(state_value_bounds[variable_count]);
  validateSettings();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
TickCounter<T>::TickCounter(const std::vector<std::vector<T>> &state_values_in,
                            const std::vector<int> &settings_in) :
    variable_count{static_cast<int>(state_values_in.size())},
    log_permutations{0.0},
    settings{settings_in.size() == 0 ? std::vector<int>(state_values_in.size(), 0) : settings_in},
    state_limits{}, state_values{}, state_value_bounds{}
{
  // Adjust the state limits to conform to the input values
  int nval = 0;
  state_limits.resize(state_values_in.size());
  state_value_bounds.resize(state_values_in.size() + 1);
  for (int i = 0; i < variable_count; i++) {
    state_limits[i] = state_values_in[i].size();
    state_value_bounds[i] = nval;
    nval += state_values_in[i].size();
  }
  log_permutations = logProduct(state_limits);
  state_value_bounds[variable_count] = nval;
  state_values.resize(nval);
  nval = 0;
  for (int i = 0; i < variable_count; i++) {
    const size_t nval_i = state_values_in[i].size();
    for (size_t j = 0; j < nval_i; j++) {
      state_values[nval] = state_values_in[i][j];
      nval++;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
int TickCounter<T>::getVariableCount() const {
  return variable_count;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
double TickCounter<T>::getLogPermutations() const {
  return log_permutations;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const std::vector<int>& TickCounter<T>::getSettings() const {
  return settings;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
int TickCounter<T>::getStateLimits(const int var_index) const {
  return state_limits[var_index];
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const std::vector<int>& TickCounter<T>::getStateLimits() const {
  return state_limits;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
T TickCounter<T>::getState(const int var_index) const {
  return state_values[state_value_bounds[var_index] + settings[var_index]];
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> TickCounter<T>::getState() const {
  std::vector<T> result(variable_count);
  for (int i = 0; i < variable_count; i++) {
    result[i] = state_values[state_value_bounds[i] + settings[i]];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> TickCounter<T>::getPossibleStates(const int var_index) const {
  std::vector<T> result(state_limits[var_index]);
  const int llim = state_value_bounds[var_index];
  const int hlim = state_value_bounds[var_index + 1];
  for (int i = llim; i < hlim; i++) {
    result[i - llim] = state_values[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
llint TickCounter<T>::getExactPermutationCount() const {
  if (getLogPermutationCount() >= static_cast<double>((sizeof(llint) * 8LLU) - 1LLU) * log(2.0)) {
    rtErr("There are too many permutations to represent as a long long integer.", "TickCounter",
          "getExactPermutationCount");
  }
  if (variable_count == 0) {
    return 0LL;
  }
  return seriesProduct<llint>(state_limits);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
double TickCounter<T>::getApproximatePermutationCount() const {
  if (variable_count == 0) {
    return 0.0;
  }
  return seriesProduct<double>(state_limits);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
double TickCounter<T>::getLogPermutationCount() const {
  if (variable_count == 0) {
    return 0.0;
  }
  return logProduct(state_limits);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TickCounter<T>::advance(const int steps) {

  // Return immediately if there are no variables to increment.
  if (variable_count == 0) {
    return;
  }
  if (steps == 1) {
    int wheel = 0;
    settings[wheel] += 1;
    while (wheel < variable_count && settings[wheel] == state_limits[wheel]) {
      settings[wheel] = 0;
      wheel++;
      if (wheel < variable_count) {
        settings[wheel] += 1;
      }
    }
  }
  else if (steps > 0) {
    int wheel = 0;
    int balance = steps;
    while (wheel < variable_count && balance > 0) {

      // The current wheel ticks forward by the modulo of the remaining balance with the wheel's
      // maximum number of ticks.
      settings[wheel] += balance % state_limits[wheel];

      // Clean up overflows resulting from the current wheel ticking forward.
      int wheel_b = wheel;
      while (wheel_b < variable_count && settings[wheel_b] >= state_limits[wheel_b]) {
        settings[wheel_b] -= state_limits[wheel_b];
        wheel_b++;
        if (wheel_b < variable_count) {
          settings[wheel_b] += 1;
        }
      }

      // Update the balance that later wheels will need to incorporate.
      balance /= state_limits[wheel];
      wheel++;
    }
  }
  else if (steps == 0) {
    return;
  }
  else {
    rtErr("A forward move of " + std::to_string(steps) + " is invalid.", "TickCounter", "advance");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TickCounter<T>::reverse(const int steps) {

  // Return immediately if there are no variables to increment.
  if (variable_count == 0) {
    return;
  }
  if (steps == 1) {
    int wheel = 0;
    settings[wheel] -= 1;
    while (wheel < variable_count && settings[wheel] < 0) {
      settings[wheel] = state_limits[wheel] - 1;
      wheel++;
      if (wheel < variable_count) {
        settings[wheel] -= 1;
      }
    }
  }
  else if (steps > 0) {
    int wheel = 0;
    int balance = steps;
    while (wheel < variable_count && balance > 0) {

      // The current wheel ticks forward by the modulo of the remaining balance with the wheel's
      // maximum number of ticks.
      settings[wheel] -= balance % state_limits[wheel];

      // Clean up overflows resulting from the current wheel ticking forward.
      int wheel_b = wheel;
      while (wheel_b < variable_count && settings[wheel_b] < 0) {
        settings[wheel_b] += state_limits[wheel_b];
        wheel_b++;
        if (wheel_b < variable_count) {
          settings[wheel_b] -= 1;
        }
      }

      // Update the balance that later wheels will need to incorporate.
      balance /= state_limits[wheel];
      wheel++;
    }
  }
  else if (steps == 0) {
    return;
  }
  else {
    rtErr("A reverse move of " + std::to_string(steps) + " is invalid.", "TickCounter", "reverse");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TickCounter<T>::set(const std::vector<int> &new_state) {
  if (static_cast<int>(new_state.size()) != variable_count) {
    rtErr("The length of the input state (" + std::to_string(new_state.size()) + ") must equal "
          "the number of variables (" + std::to_string(variable_count) + ").", "TickCounter",
          "set");
  }
  for (int i = 0; i < variable_count; i++) {
    settings[i] = new_state[i];
  }
  validateSettings();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TickCounter<T>::set(const int new_state, const int var_index) {
  variableBoundsCheck(var_index);
  settings[var_index] = new_state;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TickCounter<T>::reset() {
  for (int i = 0; i < variable_count; i++) {
    settings[i] = 0;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TickCounter<T>::loadStateValues(const std::vector<T> &state_values_in, int var_index) {
  variableBoundsCheck(var_index);
  if (state_values_in.size() != static_cast<size_t>(state_limits[var_index])) {
    rtErr("There are " + std::to_string(state_limits[var_index]) + " possible settings for "
          "variable " + std::to_string(var_index) + ", but " +
          std::to_string(state_values_in.size()) + "elements in the input array.",
          "loadStateValues");
  }
  const int imin = state_value_bounds[var_index];
  const int imax = state_value_bounds[var_index + 1];
  for (int i = imin; i < imax; i++) {
    state_values[i] = state_values_in[i - imin];
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TickCounter<T>::loadStateValues(const std::vector<std::vector<T>> &state_values_in) {


}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Trng>
void TickCounter<T>::randomize(Trng *xrs) {
  randomize(xrs, std::vector<bool>(variable_count, true));
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Trng>
void TickCounter<T>::randomize(Trng *xrs, const std::vector<bool> &apply_rng) {
  for (int i = 0; i < variable_count; i++) {
    if (apply_rng[i]) {
      settings[i] = static_cast<int>(xrs->uniformRandomNumber() * state_limits[i]);
      settings[i] -= ((settings[i] < 0) * settings[i]) +
                     ((settings[i] >= state_limits[i]) * (settings[i] - state_limits[i] + 1));
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Trng>
void TickCounter<T>::randomize(Trng *xrs, const std::vector<int> &apply_rng) {
  std::vector<bool> bapply_rng(variable_count, false);
  const	int nr	= apply_rng.size();
  for (int i = 0; i < nr; i++) {
    bapply_rng[apply_rng[i]] = true;
  }
  randomize(xrs, bapply_rng);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TickCounter<T>::validateSettings() const {
  if (settings.size() != state_limits.size()) {
    rtErr("The number of settings (" + std::to_string(settings.size()) + ") does not equal the "
          "number of mutable variables (" + std::to_string(state_limits.size()) + ").",
          "TickCounter", "validateSettings");
  }
  for (int i = 0; i < variable_count; i++) {
    if (settings[i] >= state_limits[i] || settings[i] < 0) {
      rtErr("Counter variable " + std::to_string(i) + " has a maximum of " +
            std::to_string(state_limits[i]) + " possible values.  A setting of " +
            std::to_string(settings[i]) + " is invalid.", "TickCounter", "validateSettings");
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TickCounter<T>::variableBoundsCheck(const int var_index) const {
  if (var_index >= variable_count || var_index < 0) {
    rtErr("Variable index " + std::to_string(var_index) + " is invalid for " +
          std::to_string(variable_count) + " overall counters.", "TickCounter", "set");
  }  
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2>
void loadScalarStateValues(TickCounter<T> *tc, const std::vector<T2> &ranges,
                           const std::vector<T> &periodic_boundaries) {
  const size_t nvar = tc->getVariableCount();
  if (ranges.size() != nvar) {
    rtErr("There must be a range provided for each variable in the TickCounter (" +
          std::to_string(ranges.size()) + " ranges and " + std::to_string(nvar) + " TickCounter "
          "variables are present).", "loadStateValues");
  }
  for (int i = 0; i < nvar; i++) {
    const int ni_states = tc->getStateLimits(i);
    const double dni_states = ni_states;
    std::vector<T> states(ni_states);
    for (int j = 0; j < ni_states; j++) {
      states[j] = ranges[i].x + (static_cast<double>(j) * (ranges[i].y - ranges[i].x) /
                                 dni_states);
      imageValue<double, double>(states[j], periodic_boundaries[j],
                                 ImagingMethod::PRIMARY_UNIT_CELL);
    }
    tc->loadStateValues(states, i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2>
void loadScalarStateValues(TickCounter<T> *tc, const std::vector<T2> &ranges) {

  // It would be preferable to call the previous overload, but it cannot be guaranteed that the
  // stated ranges would be preserved under all circumstances if any periodic boundary conditions
  // are designed.
  const size_t nvar = tc->getVariableCount();
  if (ranges.size() != nvar) {
    rtErr("There must be a range provided for each variable in the TickCounter (" +
          std::to_string(ranges.size()) + " ranges and " + std::to_string(nvar) + " TickCounter "
          "variables are present).", "loadStateValues");
  }
  for (int i = 0; i < nvar; i++) {
    const int ni_states = tc->getStateLimits(i);
    const double dni_states = ni_states;
    std::vector<T> states(ni_states);
    for (int j = 0; j < ni_states; j++) {
      states[j] = ranges[i].x + (static_cast<double>(j) * (ranges[i].y - ranges[i].x) /
                                 dni_states);
    }
    tc->loadStateValues(states, i);
  }
}

} // namespace stmath
} // namespace stormm
