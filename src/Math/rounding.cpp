#include "copyright.h"
#include "DataTypes/common_types.h"
#include "Constants/scaling.h"
#include "rounding.h"
#include "tickcounter.h"

namespace stormm {
namespace stmath {

using constants::mega;

//-------------------------------------------------------------------------------------------------
std::vector<uint> primeFactors(const ullint number, const std::vector<uint> &primes,
                               const int n_primes) {
  std::vector<uint> factors;
  ullint residual = number;
  const int actual_prime_count = (n_primes == 0) ? primes.size() : n_primes;
  for (int i = 0; i < actual_prime_count; i++) {
    const ullint ul_pi = primes[i];
    int nfac = 0;
    while ((residual / ul_pi) * ul_pi == residual) {
      nfac++;
      residual /= ul_pi;
    }
    if (nfac > 0) {
      factors.push_back(primes[i]);
    }
  }
  return factors;
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> primeFactorCounts(const ullint number, const std::vector<uint> &primes,
                                    const int n_primes) {
  std::vector<uint> factors((n_primes <= 0 || n_primes > primes.size()) ? primes.size() : n_primes,
                            0);
  ullint residual = number;
  const int actual_prime_count = (n_primes == 0) ? factors.size() : n_primes;
  for (int i = 0; i < actual_prime_count; i++) {
    const ullint ul_pi = primes[i];
    while ((residual / ul_pi) * ul_pi == residual) {
      factors[i] += 1;
      residual /= ul_pi;
    }
  }
  return factors;
}

//-------------------------------------------------------------------------------------------------
ulint getSmallestLot(const int element_size, const int increment, const int n_primes) {

  // Small primes
  std::vector<uint> primes = { 2, 3, 5, 7, 11, 13, 17, 19 };
  const int useable_primes = (n_primes > 8) ? 8 : n_primes;

  // Naive factorization of each number
  std::vector<uint> efac = primeFactorCounts(element_size, primes, useable_primes);
  std::vector<uint> ifac = primeFactorCounts(increment, primes, useable_primes);

  // Find the lowest common multiple
  ulint common_multiple = 1;
  for (int i = 0; i < useable_primes; i++) {
    for (int j = 0; j < std::min(efac[i], ifac[i]); j++) {
      common_multiple *= primes[i];
    }
  }

  // Divide the increment by the common multiple, then multiply the element size by the result
  return increment / common_multiple;
}

//-------------------------------------------------------------------------------------------------
ullint nearestFactor(const ullint number, const ullint target, const std::vector<uint> &primes,
                     const LimitApproach approach, const int n_primes) {
  const std::vector<uint> factors = primeFactorCounts(number, primes, n_primes);
  const int actual_pc = (n_primes <= 0 || n_primes >= primes.size()) ? primes.size() : n_primes;
  std::vector<int> maximum_prime_inclusions(actual_pc);
  int unique_prime_factors = 0;
  for (int i = 0; i < actual_pc; i++) {
    uint ilog_calc;
    switch (approach) {
    case LimitApproach::BELOW:
      ilog_calc = floor(log2(target) / log2(primes[i]));
      break;
    case LimitApproach::ABOVE:
      ilog_calc = ceil(log2(target) / log2(primes[i]));
      break;
    }
    maximum_prime_inclusions[i] = std::min(factors[i], ilog_calc);
    unique_prime_factors += (maximum_prime_inclusions[i] > 0);
  }
  if (unique_prime_factors > 0) {
    std::vector<std::vector<int>> pf_settings(unique_prime_factors);
    unique_prime_factors = 0;
    for (int i = 0; i < actual_pc; i++) {
      if (maximum_prime_inclusions[i] > 0) {
        std::vector<int> tmp_powers(maximum_prime_inclusions[i] + 1);
        int k = 1;
        for (int j = 0; j <= maximum_prime_inclusions[i]; j++) {
          tmp_powers[j] = k;
          k *= primes[i];
        }
        pf_settings[unique_prime_factors] = tmp_powers;
        unique_prime_factors++;
      }
    }
    
    // The problem of finding the best combination of prime factors to get as close as possible
    // to the target number of entries per line is hard, perhaps NP-hard, but the combinatorics
    // is not excessively large for a small target number.
    TickCounter<int> dials(pf_settings);
    const std::vector<int>& dial_settings = dials.getSettings();
    int ncombo = (dials.getLogPermutations() > 3.0) ? 1000 : dials.getExactPermutationCount();
    int best_approx = 1;
    for (int i = 0; i < ncombo; i++) {
      int trial_approx = 1;
      for (int j = 0; j < unique_prime_factors; j++) {
        trial_approx *= dials.getState(j);
      }
      switch (approach) {
      case LimitApproach::BELOW:
        if (trial_approx <= target && target - trial_approx < target - best_approx) {
          best_approx = trial_approx;
        }
        break;
      case LimitApproach::ABOVE:
        if (trial_approx >= target && trial_approx - target < best_approx - target) {
          best_approx = trial_approx;
        }
        break;
      }
      dials.advance();
    }
    return best_approx;
  }
  else {
    return number;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
size_t getPaddedMemorySize(const size_t length, const size_t growth_increment,
                           const size_t element_size) {

  // The actual growth increment
  const size_t incr = growth_increment * element_size;

  // Handle the case of zero
  if (length == 0) {
    return incr;
  }
  else if (length < 4 * incr) {
    return roundUp(length, incr);
  }
  else if (length < 8 * incr) {
    return roundUp(length, 2 * incr);
  }
  else if (length < 4 * constants::mega) {
    return roundUp(2 * length, incr);
  }
  else if (length < 16 * constants::mega) {
    return roundUp((3 * length) / 2, incr);
  }
  else if (length < 64 * constants::mega) {
    return roundUp((5 * length) / 4, incr);
  }
  else {
    return roundUp((9 * length) / 8, incr);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void partition(const int n, const int pmax, const std::vector<int> &pref,
               const std::vector<int> &dscr, std::vector<int> *result) {

  // Return immediately if there is no partitioning needed
  if (pmax >= n) {
    result->resize(1);
    result->at(0) = n;
    return;
  }

  // Adjust the maximum length downwards to meet the highest possible preferred length.
  // Determine a suitable combination of lengths that are all in the "preferred" category, or
  // otherwise avoid the "discouraged" category.
  int pbase = pmax;
  int nearest_match = 0;
  const int npref = pref.size();
  for (int i = 0; i < npref; i++) {
    if (pref[i] < pmax && (pmax - pref[i] < pmax - nearest_match)) {
      nearest_match = pref[i];
    }
  }
  if (nearest_match > 0) {
    pbase = nearest_match;
  }

  // Compute the number of base partitions.
  int base_count = n / pbase;
  
  // Compute the length of the final partition.  
  int last_partition = n - (base_count * pbase);
  if (last_partition == 0) {

    // If no further partitioning is needed, fill out the result.
    result->resize(base_count);
    int* res_ptr = result->data();
    for (int i = 0; i < base_count; i++) {
      res_ptr[i] = pbase;
    }
  }
  else {

    // Determine whether the last partition is already a preferred size.
    bool last_is_good = false;
    for (int i = 0; i < npref; i++) {
      last_is_good = (last_is_good || (pref[i] == last_partition));
    }
    if (last_is_good) {
      result->resize(base_count + 1);
      int* res_ptr = result->data();
      for (int i = 0; i < base_count; i++) {
        res_ptr[i] = pbase;
      }
      res_ptr[base_count] = last_partition;
    }
    else {

      // Determine whether the final partition can be subdivided into roughly equal portions of
      // preferred sizes.
      int2 best_pair = { 0, last_partition };
      for (int i = 0; i < npref; i++) {
        for (int j = 0; j < npref; j++) {
          if (pref[i] + pref[j] == last_partition &&
              abs(pref[i] - pref[j]) < abs(best_pair.x - best_pair.y)) {
            best_pair.x = pref[i];
            best_pair.y = pref[j];
          }
        }
      }
      if (best_pair.x != 0) {
        result->resize(base_count + 2);
        int* res_ptr = result->data();
        for (int i = 0; i < base_count; i++) {
          res_ptr[i] = pbase;
        }
        res_ptr[base_count    ] = best_pair.x;
        res_ptr[base_count + 1] = best_pair.y;
      }
      else {

        // Determine whether the final partition, which so far is at least smaller than pmax and
        // therefore permitted, is truly bad.
        const int ndscr = dscr.size();
        bool last_is_bad = false;
        for (int i = 0; i < ndscr; i++) {
          last_is_bad = (last_is_bad || dscr[i] == last_partition);
        }
        if (last_is_bad == false) {
          result->resize(base_count + 1);
          int* res_ptr = result->data();
          for (int i = 0; i < base_count; i++) {
            res_ptr[i] = pbase;
          }
          res_ptr[base_count] = last_partition;
        }
        else {

          // Find the largest preferred partition less than the current remainder, take this out
          // as its own partition, and then fold in the rest as one final partition.
          nearest_match = 0;
          for (int i = 0; i < npref; i++) {
            if (pref[i] < last_partition &&
                (last_partition - pref[i] < last_partition - nearest_match)) {
              nearest_match = pref[i];
            }
          }
          if (nearest_match > 0) {
            result->resize(base_count + 2);
            int* res_ptr = result->data();
            for (int i = 0; i < base_count; i++) {
              res_ptr[i] = pbase;
            }
            res_ptr[base_count    ] = nearest_match;
            res_ptr[base_count + 1] = last_partition - nearest_match;
          }
          else {

            // The last partition is discouraged, and nothing preferable was found.  It is probably
            // still better to have one work unit of a discouraged size than two of any other size.
            result->resize(base_count + 1);
            int* res_ptr = result->data();
            for (int i = 0; i < base_count; i++) {
              res_ptr[i] = pbase;
            }
            res_ptr[base_count] = last_partition;
          }
        }
      }
    }
  }

  // Check the result.
  const int nres = result->size();
  const int* chk_ptr = result->data();
  for (size_t i = 0; i < nres; i++) {
    if (chk_ptr[i] > pmax) {
      rtErr("Element " + std::to_string(i) + " (" + std::to_string(chk_ptr[i]) + ") exceeds the "
            "maximum partition of " + std::to_string(pmax) + ".", "partition");
    }
  }
  const int tsum = sum<int>(chk_ptr, nres);
  if (tsum != n) {
    rtErr("The sum of elements (" + std::to_string(tsum) + ") did not meet the target of " +
          std::to_string(n) + ".", "partition");
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<int> partition(const int n, const int pmax, const std::vector<int> &pref,
                           const std::vector<int> &dscr) {
  std::vector<int> result;
  partition(n, pmax, pref, dscr, &result);
  return result;
}

//-------------------------------------------------------------------------------------------------
int sanctionedDensityGridDimension(const ullint big_product, const double* invu,
                                   const UnitCellAxis edge, const std::vector<uint> &primes,
                                   const double max_spacing) {
  int ne;
  switch (edge) {
  case UnitCellAxis::A:
    ne = ceil(invu[0] / max_spacing);
    break;
  case UnitCellAxis::B:
    ne = ceil(sqrt((invu[3] * invu[3]) + (invu[4] * invu[4])) / max_spacing);
    break;
  case UnitCellAxis::C: 
    ne = ceil(sqrt((invu[6] * invu[6]) + (invu[7] * invu[7]) + (invu[8] * invu[8])) / max_spacing);
    break;
  }
  int result = nearestFactor(big_product, ne, primes, LimitApproach::ABOVE);

  // The combination of radices 7 and 11 is very costly for FFTs, but either radix by itself along
  // one dimension is fine.  Check both to see which provides the minimal result.
  if (result % 77 == 0) {
    const int r_seven  = nearestFactor(big_product / 11LLU, ne, primes, LimitApproach::ABOVE);
    const int r_eleven = nearestFactor(big_product / 7LLU,  ne, primes, LimitApproach::ABOVE);
    result = std::min(r_seven, r_eleven);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int sanctionedDensityGridDimension(const ullint big_product, const double* invu,
                                   const CartesianDimension edge,
                                   const std::vector<uint> &primes, const double max_spacing) {
  switch (edge) {
  case CartesianDimension::X:
    return sanctionedDensityGridDimension(big_product, invu, UnitCellAxis::A, primes, max_spacing);
  case CartesianDimension::Y:
    return sanctionedDensityGridDimension(big_product, invu, UnitCellAxis::B, primes, max_spacing);
  case CartesianDimension::Z:
    return sanctionedDensityGridDimension(big_product, invu, UnitCellAxis::C, primes, max_spacing);
  }
  __builtin_unreachable();
}
} // namespace stmath
} // namespace stormm
