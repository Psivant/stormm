// -*-c++-*-
#ifndef STORMM_HPC_RANDOM_H
#define STORMM_HPC_RANDOM_H

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/stormm_vector_types.h"
#include "Random/random.h"

namespace stormm {
namespace random {

using card::GpuDetails;
using card::Hybrid;
  
/// \brief Launch the eponymous kernel to seed Xoroshiro128+ random number generator states.  The
///        algorithm's long jump is executed using CPU code as this work must be serialized.
///
/// \param state_vector  Blank vector allocated to hold the requested number of random number
///                      generator states.  Filled and returned.
/// \param igseed        Random number seed for the first of the generator states (all others
///                      follow from this)
/// \param scrub_cycles  Number of times to run the random number generator and discard results in
///                      order to raise the quality of random numbers coming out of it
/// \param gpu           Details of the GPU on which to launch this operation
void initXoroshiro128pArray(Hybrid<ullint2> *state_vector, int igseed, int scrub_cycles,
                            const GpuDetails &gpu);

/// \brief Launch the eponymous kernel to seed Xoshiro256++ random number generator states.  The
///        algorithm's long jump is executed using CPU code as this work must be serialized.
///
/// \param state_xy      Blank vector allocated to hold the first half of all random number
///                      generator states, in the requested number.  Filled and returned.
/// \param state_zw      Blank vector allocated to hold the second half of all random number
///                      generator states, in the requested number.  Filled and returned.
/// \param igseed        Random number seed for the first of the generator states (all others
///                      follow from this)
/// \param scrub_cycles  Number of times to run the random number generator and discard results in
///                      order to raise the quality of random numbers coming out of it
/// \param gpu           Details of the GPU on which to launch this operation
void initXoshiro256ppArray(Hybrid<ullint2> *state_xy, Hybrid<ullint2> *state_zw, int igseed,
                           int scrub_cycles, const GpuDetails &gpu);

/// \brief Fill a cache of random numbers using a series (an array) of state vectors (generators).
///        The second of two ullint2 state vectors may be supplied as nullptr if only 128-bit
///        states are in use.
///
/// Overloaded:
///   - Work with 128-bit or 256-bit generator states.
///   - Operate on C-style arrays or Hybrid objects
///   - Produce single- or double-precision random number results on the requested distribution
///     (this is not templated like the eponymous CPU function, to avoid 
///
/// \param state_xy     First halves of each 256-bit generator state vector, or the array of
///                     128-bit state vectors (state_zw should be the null pointer in this case)
/// \param state_zw     Second havles of each 256-bit generator state vector, if required
/// \param cache        Array of real-valued random number results to fill out (the template
///                     parameter indicates the data type of this array)
/// \param length       Trusted length of state_xy and state_zw
/// \param depth        Quantity of random numbers for each generator to produce during checkout
///                     from the state vector arrays.  The (warp size-padded) length parameter
///                     times depth gives the size of the cache array.
/// \param method       Method for generating random numbers (some XOR-shift technique)
/// \param product      Shape of the random distribution over which to take results
/// \param index_start  Starting index at which to begin drawing upon random number generators
/// \param index_end    Upper bound of random number generators from which to draw results
/// \{
void fillRandomCache(ullint2* state_xy, ullint2* state_zw, double* cache, size_t length,
                     size_t depth, RandomAlgorithm method, RandomNumberKind product,
                     size_t index_start, size_t index_end, const GpuDetails &gpu);

void fillRandomCache(Hybrid<ullint2> *state_xy, Hybrid<ullint2> *state_zw, Hybrid<double> *cache,
                     size_t length, size_t depth, RandomAlgorithm method, RandomNumberKind product,
                     size_t index_start, size_t index_end, const GpuDetails &gpu);

void fillRandomCache(ullint2* state_xy, ullint2* state_zw, float* cache, size_t length,
                     size_t depth, RandomAlgorithm method, RandomNumberKind product,
                     size_t index_start, size_t index_end, const GpuDetails &gpu);

void fillRandomCache(Hybrid<ullint2> *state_xy, Hybrid<ullint2> *state_zw, Hybrid<float> *cache,
                     size_t length, size_t depth, RandomAlgorithm method, RandomNumberKind product,
                     size_t index_start, size_t index_end, const GpuDetails &gpu);
/// \}

} // namespace random
} // namespace stormm

#endif
