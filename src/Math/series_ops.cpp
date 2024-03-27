#include <algorithm>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "series_ops.h"
#include "vector_ops.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
std::vector<uint> numberSeriesToBitMask(const int* number_series, const size_t length,
                                        const int output_size) {
  const int actual_size = (output_size == 0LLU) ? std::max(maxValue(number_series, length), 1) :
                                                  output_size;
  const int n_bits = sizeof(uint) * 8;
  const int n_uint = (actual_size + n_bits - 1) / n_bits;
  std::vector<uint> result(n_uint, 0);

  // Loop over all polar hydrogens and check the appropriate bits
  for (size_t i = 0; i < length; i++) {
    const int mask_idx = number_series[i] / n_bits;
    const int bit_idx = number_series[i] - (mask_idx * n_bits);
    result[mask_idx] |= (0x1 << bit_idx);
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
std::vector<uint> numberSeriesToBitMask(const std::vector<int> &number_series,
                                        const int output_size) {
  return numberSeriesToBitMask(number_series.data(), number_series.size(), output_size);
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> numberSeriesToBitMask(const Hybrid<int> &number_series, const int output_size) {
  return numberSeriesToBitMask(number_series.data(), number_series.size(), output_size);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> getSubsetIndexPattern(const std::vector<int> &x_subset, int *n_subset_indices) {

  // Obtain the maximum index
  const size_t nidx = x_subset.size();
  const int max_index = maxValue(x_subset);
  std::vector<bool> index_found(max_index + 1, false);
  for (size_t i = 0; i < nidx; i++) {
    index_found[x_subset[i]] = true;
  }
  int new_index_counter = 0;
  std::vector<int> result(max_index + 1, -1);
  for (size_t i = 0; i <= max_index; i++) {
    if (index_found[i]) {
      result[i] = new_index_counter;     
      new_index_counter++;
    }
  }
  if (n_subset_indices != nullptr) {
    *n_subset_indices = new_index_counter;
  }
  return result;
}
  
} // namespace stmath
} // namespace stormm
