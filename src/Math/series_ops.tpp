// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> incrementingSeries(const T start_value, const T end_value,
                                                        const T increment) {
  const T zero = 0;
  const T actual_increment = (end_value - start_value > zero) ?  std::abs(increment) :
                                                                -std::abs(increment);
  std::vector<T> result(std::abs((end_value - start_value) / increment));
  T v = start_value;
  const size_t nval = result.size();
  for (size_t i = 0; i < nval; i++) {
    result[i] = v;
    v += actual_increment;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> extractIndexedValues(const T* original_values,
                                                          const size_t values_length,
                                                          const std::vector<int> reduced_indices,
                                                          const int reduced_length) {
  const int actual_length = (reduced_length == 0) ? maxValue(reduced_indices) : reduced_length;
  const size_t ridx_size = reduced_indices.size();
  if (values_length < ridx_size) {
    rtErr("The array of original values is smaller than the number of unique indices.",
          "extractIndexedValues");    
  }
  std::vector<T> result;
  for (size_t i = 0; i < ridx_size; i++) {
    if (reduced_indices[i] >= 0) {
      result[reduced_indices[i]] = original_values[i];
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> extractIndexedValues(const std::vector<T> &original_values,
                                                          const std::vector<int> reduced_indices,
                                                          const int reduced_length) {
  return extractIndexedValues(original_values.data(), original_values.size(), reduced_indices,
                              reduced_length);
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> extractIndexedValues(const Hybrid<T> &original_values,
                                                          const std::vector<int> reduced_indices,
                                                          const int reduced_length) {
  return extractIndexedValues(original_values.data(), original_values.size(), reduced_indices,
                              reduced_length);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tloc>
void indexingArray(const Tdata* raw_values, Tloc* value_locations, Tloc* value_bounds,
                   const size_t value_count, const size_t value_limit) {

  // Zero the bounds array
  for (size_t i = 0; i <= value_limit; i++) {
    value_bounds[i] = 0;
  }

  // Loop over the raw values and tally their quantities.  Cast the raw_values elements to size_t
  // in order to perform the safety check and prepare for indexing the proper element of the bounds
  // array.
  for (size_t i = 0; i < value_count; i++) {
    const size_t rvi = raw_values[i];
    if (rvi < value_limit) {
      value_bounds[rvi] += 1;
    }
  }
  prefixSumInPlace(value_bounds, value_limit + 1, PrefixSumType::EXCLUSIVE, "indexingArray");

  // Loop back over all values and sort them into the result.  Use the bounds array to guide the
  // indexing, then rewind it after the locations have been catalogged.
  for (size_t i = 0; i < value_count; i++) {
    const size_t rvi = raw_values[i];
    if (rvi < value_limit) {
      const Tloc vbi = value_bounds[rvi];
      value_locations[vbi] = i;
      value_bounds[rvi] = vbi + 1;
    }
  }
  for (size_t i = value_limit; i > 0; i--) {
    value_bounds[i] = value_bounds[i - 1];
  }
  value_bounds[0] = 0;
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tloc>
void indexingArray(const std::vector<Tdata> &raw_values, std::vector<Tloc> *value_locations,
                   std::vector<Tloc> *value_bounds, size_t value_limit) {
  const size_t actual_value_limit = (value_limit == 0) ? value_bounds->size() - 1 : value_limit;
  indexingArray(raw_values.data(), value_locations->data(), value_bounds->data(),
                raw_values.size(), actual_value_limit);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tloc>
void indexingArray(const Hybrid<Tdata> &raw_values, Hybrid<Tloc> *value_locations,
                   Hybrid<Tloc> *value_bounds, size_t value_limit) {
  const size_t actual_value_limit = (value_limit == 0) ? value_bounds->size() - 1 : value_limit;
  indexingArray(raw_values.data(), value_locations->data(), value_bounds->data(),
                raw_values.size(), actual_value_limit);
}

} // namespace stmath
} // namespace stormm
