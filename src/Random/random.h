// -*-c++-*-
#ifndef STORMM_RANDOM_H
#define STORMM_RANDOM_H

#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/scaling.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/rounding.h"
#include "Numerics/split_fixed_precision.h"
#include "random_enumerators.h"

namespace stormm {
namespace random {

using card::Hybrid;
using card::HybridKind;
using card::HybridTargetLevel;
using constants::giga_zu;
using data_types::isFloatingPointScalarType;
using stmath::roundUp;

/// \brief Constants for the "ran2" pseudo-random number generator
/// \{
constexpr int im1 = 2147483563;
constexpr int im2 = 2147483399;
constexpr double am = 1.0 / im1;
constexpr int im1_minus_one = im1 - 1;
constexpr int ia1 = 40014;
constexpr int ia2 = 40692;
constexpr int iq1 = 53668;
constexpr int iq2 = 52774;
constexpr int ir1 = 12211;
constexpr int ir2 = 3791;
constexpr int ntab = 32;
constexpr int ndiv = 1 + im1_minus_one / ntab;
constexpr double double_increment = 2.2205e-16;
constexpr double ran2_max = 1.0 - double_increment;
/// \}

/// \brief The default random seed for STORMM
constexpr int default_random_seed = 827493;

/// \brief Numbers of "scrub cycles" expected to sufficiently randomize non-zero seeds in any of
///        the various XOR-shift generators.
/// \{
constexpr int default_xoroshiro128p_scrub = 25;
constexpr int default_xoshiro256pp_scrub  = 50;
/// \}
  
/// \brief The maximum number of long jumps to take with any of the scrambled linear pseudo-random
///        number generators.  The CPU will take the long jumps, while individual GPU threads make
///        standard jumps ahead in the pseudo-random sequence.  The CPU is 10-30 times faster
///        than a single GPU thread, so this number could in fact scale with the number of state
///        vectors being requested, but the optimization would be on the order of milliseconds
///        and making this a constant makes it easier to reproduce the sequence at a given point.
constexpr int max_xo_long_jumps = 1024;

/// \brief When using the XOR-shift random number generators, uniform random numbers are produced
///        on the range [0.0, 1.0).  The Box-Muller transformation to produce random numbers on a
///        normal distribution is valid for the range (0.0, 1.0), otherwise it produces nan when
///        the log() or sqrt() functions encounter a zero, respectively.  In order to solve this,
///        XOR-shift random number generators will commit to producing uniform random numbers with
///        less random bits (see the implementation), but make use of the 24th or 53rd bits in the
///        form of an offset one-half the width of the smallest representable increment on the
///        range [0.5, 1.0).
/// \{
constexpr double rng_unit_bin_offset  = 1.1102230246251565404236316680908203125e-16;
constexpr float rng_unit_bin_offset_f = 5.960464477539062500e-08f;
/// \}
  
/// \brief Stores the state of a Ran2 pseudo-random number generator.  Member functions produce
///        random numbers along various distributions, as required.  While it is not as performant
///        to have a member function, these random number generators are intended for convenience
///        and unit testing purposes.  Developers that wish to use higher-performance random number
///        generators should use the Xoroshiro128Generator (below) or link other libraries for the
///        C++ layer (i.e. Intel MKL) or libraries that run on the GPU layer (i.e. cuRAND). 
class Ran2Generator {
public:

  /// \brief Initialize the unit testing on-board pseudo-random number generator.
  ///
  /// \param igseed  The pseudo-randome number seed
  Ran2Generator(int igseed = default_random_seed);  

  /// \brief Return a single random number distributed over a uniform distribution [0, 1).  This is
  ///        an internal generator based on an integer vector state which will provide reproducible
  ///        results across platforms.
  double uniformRandomNumber();

  /// \brief Return a single random number distributed over a uniform distribution [0, 1), in
  ///        single-precision.  While the XOR-shift random generators have a distinct method for
  ///        producing single-precision random numbers (selecting the highest-quality bits), this
  ///        merely calls uniformRandomNumber() and exists for compatibility with the other
  ///        generators in templated xxxxRand() functions.
  float spUniformRandomNumber();

  /// \brief Return a normally distributed random number with standard deviation 1.0.  This works
  ///        off of the uniform random number generator and will thus advance the state vector of
  ///        the random number generator that produces the result.
  double gaussianRandomNumber();

  /// \brief Produce a single-precision random number on the normal distribution with a standard
  ///        deviation of 1.0.  Like spUniformRandomNumber(), this function merely calls
  ///        gaussianRandomNumber() and exists for compatibility with other random number
  ///        generators in templated functions.
  float spGaussianRandomNumber();

private:
  int iunstable_a;
  int iunstable_b;
  int state_sample;
  int state_vector[ntab + 8];
};

/// \brief Constants for Xoroshiro128+ generator.
/// \{
constexpr ullint xrs128p_jump_i      = 0xdf900294d8f554a5LLU;
constexpr ullint xrs128p_jump_ii     = 0x170865df4b3201fcLLU;
constexpr ullint xrs128p_longjump_i  = 0xd2a98b26625eee7bLLU;
constexpr ullint xrs128p_longjump_ii = 0xdddf9b1090aa7ac1LLU;
/// \}

/// \brief Constants for Xoshiro256++ generator.
/// \{
constexpr ullint xrs256pp_jump_i       = 0x180ec6d33cfd0abaLLU;
constexpr ullint xrs256pp_jump_ii      = 0xd5a61266f0c9392cLLU;
constexpr ullint xrs256pp_jump_iii     = 0xa9582618e03fc9aaLLU;
constexpr ullint xrs256pp_jump_iv      = 0x39abdc4529b1661cLLU;
constexpr ullint xrs256pp_longjump_i   = 0x76e15d3efefdcbbfLLU;
constexpr ullint xrs256pp_longjump_ii  = 0xc5004e441c522fb3LLU;
constexpr ullint xrs256pp_longjump_iii = 0x77710069854ee241LLU;
constexpr ullint xrs256pp_longjump_iv  = 0x39109bb02acbe635LLU;
/// \}

/// \brief The "Xoroshiro128+" random number generator.  It's decent, but not recommended for
///        situations where a quarter million or more streams are producing random number sequences
///        in unison.  Those streams can end up containing correlated patterns.  This generator has
///        been shown to fail BigCrush.
class Xoroshiro128pGenerator {
public:

  /// \brief The constructor can start or restart the generator.
  ///
  /// Overloaded:
  ///   - Initialize the generator upon construction.
  ///   - Take the given state.
  ///
  /// \param igseed    The pseudo-randome number seed
  /// \param niter     Number of iterations to run through while equilibrating the generator
  /// \param state_in  State to accept
  /// \{
  Xoroshiro128pGenerator(int igseed = default_random_seed, int niter = 25);
  Xoroshiro128pGenerator(const ullint2 state_in);
  /// \}
  
  /// \brief Return a single random number distributed over a uniform distribution [0, 1).  This is
  ///        an internal generator based on an integer vector state which will provide reproducible
  ///        results across platforms.
  double uniformRandomNumber();

  /// \brief Return a single random number distributed over a uniform distribution [0, 1), in
  ///        single-precision.  While the XOR-shift random generators have a distinct method for
  ///        producing single-precision random numbers (selecting the highest-quality bits), this
  ///        merely calls uniformRandomNumber() and exists for compatibility with the other
  ///        generators in templated xxxxRand() functions.
  float spUniformRandomNumber();

  /// \brief Return a normally distributed random number with standard deviation 1.0.  This works
  ///        off of the uniform random number generator and will thus advance the state vector of
  ///        the random number generator that produces the result.
  double gaussianRandomNumber();

  /// \brief Produce a single-precision random number on the normal distribution with a standard
  ///        deviation of 1.0.  Like spUniformRandomNumber(), this function merely calls
  ///        gaussianRandomNumber() and exists for compatibility with other random number
  ///        generators in templated functions.
  float spGaussianRandomNumber();

  /// \brief Jump forward 2^64 iterations in the sequence.
  void jump();
  
  /// \brief Jump forward 2^96 iterations in the sequence.
  void longJump();

  /// \brief Reveal the current state of the generator.
  ullint2 revealState() const;

  /// \brief Reveal the random bit string.
  ullint revealBitString() const;

  /// \brief Set the current state of the generator.
  void setState(const ullint2 state_in);
  
private:
  ullint2 state;  ///< 128-bit state vector for the generator
  
  /// \brief Transform a 32-bit int into 128 well-spaced bits for the seed
  ///
  /// \param igseed  Random seed passed down from the constructor
  ullint2 seed128(int igseed);
  
  /// \brief Iterate to the next random number
  ullint next();

  /// \brief Jump forward in the sequence according to the jump or long-jump perturbations
  ///
  /// \param stride  The perturbation to use.  This will be either xrs128p_jump_(i,ii) or
  ///                xrs128p_longjump_(i,ii).
  void fastForward(const ullint2 stride);
};

/// \brief The "Xoshiro256++" random number generator.  While not cryptographically useful, it is
///        a rock-solid random number generator for both floating-point and 64-bit integer results.
class Xoshiro256ppGenerator {
public:

  /// \brief The constructor can start or restart the generator.
  ///
  /// Overloaded:
  ///   - Initialize the generator upon construction.
  ///   - Take the given state.
  ///
  /// \param igseed    The pseudo-randome number seed
  /// \param niter     Number of iterations to run through while equilibrating the generator
  /// \param state_in  State to accept
  Xoshiro256ppGenerator(int igseed = default_random_seed, int niter = 50);
  Xoshiro256ppGenerator(const ullint4 state_in);

  /// \brief Return a single random number distributed over a uniform distribution [0, 1).  This is
  ///        an internal generator based on an integer vector state which will provide reproducible
  ///        results across platforms.
  double uniformRandomNumber();

  /// \brief Return a single random number distributed over a uniform distribution [0, 1), in
  ///        single-precision.  While the XOR-shift random generators have a distinct method for
  ///        producing single-precision random numbers (selecting the highest-quality bits), this
  ///        merely calls uniformRandomNumber() and exists for compatibility with the other
  ///        generators in templated xxxxRand() functions.
  float spUniformRandomNumber();

  /// \brief Return a normally distributed random number with standard deviation 1.0.  This works
  ///        off of the uniform random number generator and will thus advance the state vector of
  ///        the random number generator that produces the result.
  double gaussianRandomNumber();

  /// \brief Produce a single-precision random number on the normal distribution with a standard
  ///        deviation of 1.0.  Like spUniformRandomNumber(), this function merely calls
  ///        gaussianRandomNumber() and exists for compatibility with other random number
  ///        generators in templated functions.
  float spGaussianRandomNumber();

  /// \brief Jump forward 2^128 iterations in the sequence.
  void jump();
  
  /// \brief Jump forward 2^192 iterations in the sequence.
  void longJump();

  /// \brief Reveal the current state of the generator.
  ullint4 revealState() const;

  /// \brief Reveal the random bit string.
  ullint revealBitString() const;

  /// \brief Set the current state of the generator.
  void setState(const ullint4 state_in);
  
private:
  ullint4 state;  ///< 128-bit state vector for the generator
  
  /// \brief Transform a 32-bit int into 128 well-spaced bits for the seed
  ///
  /// \param igseed  Random seed passed down from the constructor
  ullint4 seed256(int igseed);
  
  /// \brief Iterate to the next random number
  ullint next();

  /// \brief Jump forward in the sequence according to the jump or long-jump perturbations
  ///
  /// \param stride  The perturbation to use.  This will be either xrs256pp_jump_(i,ii,iii,iv) or
  ///                xrs256pp_longjump_(i,ii,iii,iv).
  void fastForward(const ullint4 stride);
};

/// \brief An series of "Xorshiro128+" generators, with state vectors for all of them and the
///        means for seeding the series based on long jumps from a single state vector.
template <typename T> class RandomNumberMill {
public:

  /// \brief The constructor can work off of a simple random initialization seed integer or a
  ///        specific state vector for the first generator in the series.
  ///
  /// \param generators_in   The count of generators in the series
  /// \param depth_in        Quantity of random numbers from each generator to store in the bank
  /// \param init_kind       Style in which to initialize the random numbers in the bank,
  ///                        i.e. from a uniform or a normal distribution
  /// \param igseed_in       The seed for the first generator in the series
  /// \param niter           The number of iterations to use in initializing each generator
  /// \param bank_limit      The maximum length of the random number cache
  /// \param state_in        The state to apply to generator zero, thus determining the initial
  ///                        states of all other generators
  /// \{
  RandomNumberMill(const ullint2 state_in, size_t generators_in = 1LLU, size_t depth_in = 2LLU,
                   RandomNumberKind init_kind = RandomNumberKind::GAUSSIAN,
                   size_t bank_limit = constants::giga_zu);

  RandomNumberMill(const ullint4 state_in, size_t generators_in = 1LLU, size_t depth_in = 2LLU,
                   RandomNumberKind init_kind = RandomNumberKind::GAUSSIAN,
                   size_t bank_limit = constants::giga_zu);

  RandomNumberMill(size_t generators_in = 1LLU, size_t depth_in = 2LLU,
                   RandomAlgorithm style_in = RandomAlgorithm::XOSHIRO_256PP,
                   RandomNumberKind init_kind = RandomNumberKind::GAUSSIAN,
                   int igseed_in = default_random_seed, int niter = default_xoshiro256pp_scrub,
                   size_t bank_limit = constants::giga_zu);
  /// \}

  /// \brief Get the number of (forward-jumped) generators.
  size_t getGeneratorCount() const;

  /// \brief Get the depth of the bank for each generator.
  size_t getDepth() const;
  
  /// \brief Get the number of generators for which to refesh all banked values at one time.
  size_t getRefreshStride() const;

  /// \brief Get a random number out of the bank.
  ///
  /// \param generator_index  The generator series from which to obtain the value
  /// \param layer_index      Layer from which to obtain the value
  T getBankValue(size_t generator_index, size_t layer_index) const;
  
  /// \brief Populate a portion of this object's bank with random numbers from each of the
  ///        respective generators.
  ///
  /// \param first_gen  Index of the first generator to draw random numbers from
  /// \param last_gen   Index of the generator before which to stop drawing new random numbers
  void uniformRandomNumbers(size_t first_gen, size_t last_gen);

  /// \brief Populate a portion of this object's bank with normally distributed random numbers
  ///        from each of the respective generators.
  ///
  /// \param first_gen  Index of the first generator to draw random numbers from
  /// \param last_gen   Index of the generator before which to stop drawing new random numbers
  void gaussianRandomNumbers(size_t first_gen, size_t last_gen);

private:
  RandomAlgorithm style;     ///< The kind of random number generator that this mill uses
  size_t generators;         ///< The quantity of (pseudo-) random number generators in the series
  size_t depth;              ///< The depth of each generator's random numbers in the memory bank
  size_t refresh_stride;     ///< The number of segments in which the generators series can be
                             ///<   refreshed, calling a subset of the generators to recalculate
                             ///<   the banked random values for that generator / lane.
  Hybrid<ullint2> state_xy;  ///< State vectors for all random number generators in the series,
                             ///<   sufficient for a Xoroshiro128+ generator or any other requiring
                             ///<   a 128-bit state
  Hybrid<ullint2> state_zw;  ///< State vectors for all random number generators in the series,
                             ///<   extended for Xoshiro256++ and any other generator requiring a
                             ///<   256-bit state vector
  Hybrid<T> bank;            ///< Bank of random numbers pre-computed and saved for later use.  The
                             ///<   numbers are stored only in the specifictype format of the
                             ///<   series, float or double, to save space and avoid ambiguities
                             ///<   as to when the numbers should be refreshed if multiple levels
                             ///<   of precision in the random numbers are required.

  /// \brief Check the inputs for validity and safety, to avoid excessive allocations.
  ///
  /// \param bank_limit  The maximum length of the random number cache
  void checkDimensions(size_t bank_limit);

  /// \brief Prime the bank with its first complement of random numbers.
  ///
  /// \param init_kind  Style in which to initialize the random numbers in the bank, i.e. from a
  ///                   uniform or a normal distribution
  void initializeBank(RandomNumberKind init_kind);
};

/// \brief Return double-precision random number(s) distributed over a uniform distribution.
///
/// Overloaded:
///   - Issue a single random number
///   - Issue multiple random numbers as a std::vector<double>
///   - Issue random numbers to fill a Standard Template Library vector of arbitrary type (these
///     random numbers will still be issued as double-precision, using a procedure that draws on
///     53 bits of the 64-bit string)
///
/// \param count     The quantity of random numbers to produce
/// \param rows      The total rows of a matrix to fill up with random numbers
/// \param columns   The total columns of a matrix to fill up with random numbers
/// \param xv        Pre-allocated array to fill with random numbers
/// \param fp_scale  Fixed precision scaling factor to apply when filling arrays of integer type
/// \param scale     Scaling factor to apply to all random numbers
/// \{
template <typename Tprng> double uniformRand(Tprng *rng, double scale = 1.0);

template <typename Tprng> std::vector<double> uniformRand(Tprng *rng, size_t count,
                                                          double scale);

template <typename Tprng>
std::vector<double> uniformRand(Tprng *rng, size_t rows, size_t columns, double scale,
                                RngFillMode mode = RngFillMode::COLUMNS);

template <typename Tprng, typename Tprod>
void uniformRand(Tprng *rng, std::vector<Tprod> *xv, double scale = 1.0, double fp_scale = 1.0);
/// \}

/// \brief Return single-precision random number(s) distributed over a uniform distribution.
///        Descriptions of parameters and overloading follows from uniformRand() above.
/// \{
template <typename Tprng> float spUniformRand(Tprng *rng, float scale = 1.0);

template <typename Tprng> std::vector<float> spUniformRand(Tprng *rng, size_t count, float scale);

template <typename Tprng>
std::vector<float> spUniformRand(Tprng *rng, size_t rows, size_t columns, float scale,
                                 RngFillMode mode = RngFillMode::COLUMNS);

template <typename Tprng, typename Tprod>
void spUniformRand(Tprng *rng, std::vector<Tprod> *xv, float scale = 1.0, float fp_scale = 1.0);
/// \}

/// \brief Return double-precision random number(s) distributed over a normal distribution.
///        Descriptions of parameters and overloading follows from uniformRand() above.
/// \{
template <typename Tprng> double guassianRand(Tprng *rng, double scale = 1.0);

template <typename Tprng> std::vector<double> guassianRand(Tprng *rng, size_t count, double scale);

template <typename Tprng>
std::vector<double> guassianRand(Tprng *rng, size_t rows, size_t columns, double scale,
                                 RngFillMode mode = RngFillMode::COLUMNS);

template <typename Tprng, typename Tprod>
void guassianRand(Tprng *rng, std::vector<Tprod> *xv, double fp_scale = 1.0, double scale = 1.0);
/// \}

/// \brief Return single-precision random number(s) distributed over a normal distribution.  
///        Descriptions of parameters and overloading follows from uniformRand() above.
/// \{
template <typename Tprng> float spGuassianRand(Tprng *rng, float scale = 1.0);

template <typename Tprng> std::vector<float> spGuassianRand(Tprng *rng, size_t count, float scale);

template <typename Tprng>
std::vector<float> spGuassianRand(Tprng *rng, size_t rows, size_t columns, float scale,
                                  RngFillMode mode = RngFillMode::COLUMNS);

template <typename Tprng, typename Tprod>
void spGuassianRand(Tprng *rng, std::vector<Tprod> *xv, float fp_scale = 1.0, float scale = 1.0);
/// \}

/// \brief Initialize an array of Xoroshiro128p generators using CPU-based code.
///
/// Overloaded:
///   - Operate on a C-style array
///   - Operate on a Standard Template Library Vector of RNG state vectors
///   - Operate on a Hybrid object (HOST only--see the overloaded form in hpc_random.h for GPU
///     functionality)
///
/// \brief state_vector  The array of RNG state vectors
/// \brief n_generators  Trusted length of the state_vector array, the number of RN generators
/// \brief igseed        Seed to apply in creating the first random number generator, from which
///                      all others are long jumps plus short jumps 
/// \param scrub_cycles  Number of times to run the random number generator and discard results in
///                      order to raise the quality of random numbers coming out of it
/// \{
void initXoroshiro128pArray(ullint2* state_vector, int n_generators, int igseed,
                            int scrub_cycles = default_xoroshiro128p_scrub);
void initXoroshiro128pArray(std::vector<ullint2> *state_vector, int igseed,
                            int scrub_cycles = default_xoroshiro128p_scrub);
void initXoroshiro128pArray(Hybrid<ullint2> *state_vector, int igseed,
                            int scrub_cycles = default_xoroshiro128p_scrub);
/// \}

/// \brief Initialize an array of Xoroshiro128p generators using CPU-based code.
///
/// Overloaded:
///   - Operate on a C-style array
///   - Operate on a Standard Template Library Vector of RNG state vectors
///   - Operate on a Hybrid object (HOST only--see the overloaded form in hpc_random.h for GPU
///     functionality)
///
/// \brief state_xy      
/// \brief n_generators  Trusted length of the state_vector array, the number of RN generators
/// \brief igseed        Seed to apply in creating the first random number generator, from which
///                      all others are long jumps plus short jumps
/// \param scrub_cycles  Number of times to run the random number generator and discard results in
///                      order to raise the quality of random numbers coming out of it
/// \{
void initXoshiro256ppArray(ullint2* state_xy, ullint2* state_zw, int n_generators, int igseed,
                           int scrub_cycles = default_xoshiro256pp_scrub);
void initXoshiro256ppArray(std::vector<ullint2> *state_xy, std::vector<ullint2> *state_zw,
                           int igseed, int scrub_cycles = default_xoshiro256pp_scrub);
void initXoshiro256ppArray(Hybrid<ullint2> *state_xy, Hybrid<ullint2> *state_zw, int igseed,
                           int scrub_cycles = default_xoshiro256pp_scrub);
/// \}

/// \brief Fill a cache of random numbers using a series (an array) of state vectors (generators).
///        The second of two ullint2 state vectors may be supplied as nullptr if only 128-bit
///        states are in use.
///
/// Overloaded:
///   - Work with 128-bit or 256-bit generator states.
///   - Operate on C-style arrays, Standard Template Library Vectors, or Hybrid objects
///   - An alternative form in the hpc_random library makes use of the eponymous GPU kernels.
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
template <typename T>
void fillRandomCache(ullint2* state_xy, ullint2* state_zw, T* cache, size_t length, size_t depth,
                     RandomAlgorithm method, RandomNumberKind product, size_t index_start,
                     size_t index_end);

template <typename T>
void fillRandomCache(std::vector<ullint2> *state_xy, std::vector<ullint2> *state_zw,
                     std::vector<T> *cache, size_t length, size_t depth, RandomAlgorithm method,
                     RandomNumberKind product, size_t index_start, size_t index_end);

template <typename T>
void fillRandomCache(Hybrid<ullint2> *state_xy, Hybrid<ullint2> *state_zw, std::vector<T> *cache,
                     size_t length, size_t depth, RandomAlgorithm method, RandomNumberKind product,
                     size_t index_start, size_t index_end);
/// \}

/// \brief Add random noise to a series of numbers.  If the number series is represented in a
///        floating point type, its scaling factor is assumed to be 1.0
///
/// Overloaded:
///   - Operate on C-style arrays with a trusted length, Standard Template Library vectors, or
///     Hybrid objects
///   - Operate on the 95-bit fixed-precision fused data type
///   - Operate on one array or three arrays of numbers, anticipating applications to
///     three-dimensional coordinate sets
/// 
/// \param prng    The random number generator to use
/// \param x       The (first) array of numbers to add noise onto
/// \param y       The second array of numbers to add noise onto
/// \param z       The third array of numbers to add noise onto
/// \param x_ovrf  Overflow bits for the (first) array of fixed precision numbers
/// \param y_ovrf  Overflow bits for the third array of fixed precision numbers
/// \param z_ovrf  Overflow bits for the second array of fixed precision numbers
/// \param length  The length of the number series (if provided as C-style arrays)
/// \param mult    Strength of the random noise
/// \param scale   Scaling factor applied to fixed-precision number representations in x, y, and z
/// \param kind    The type of random noise to apply (default Gaussian, normal distribution)
/// \{
template <typename Trng, typename Tvar>
void addRandomNoise(Trng *prng, Tvar* x, size_t length, double mult, double scale = 1.0,
                    RandomNumberKind kind = RandomNumberKind::GAUSSIAN);

template <typename Trng, typename Tvar>
void addRandomNoise(Trng *prng, std::vector<Tvar> *x, double mult, double scale = 1.0,
                    RandomNumberKind kind = RandomNumberKind::GAUSSIAN);

template <typename Trng, typename Tvar>
void addRandomNoise(Trng *prng, Hybrid<Tvar> *x, double mult, double scale = 1.0,
                    RandomNumberKind kind = RandomNumberKind::GAUSSIAN);

template <typename Trng, typename Tvar>
void addRandomNoise(Trng *prng, Tvar* x, Tvar* y, Tvar *z, size_t length, double mult,
                    double scale = 1.0, RandomNumberKind kind = RandomNumberKind::GAUSSIAN);

template <typename Trng, typename Tvar>
void addRandomNoise(Trng *prng, std::vector<Tvar> *x, std::vector<Tvar> *y, std::vector<Tvar> *z,
                    double mult, double scale = 1.0,
                    RandomNumberKind kind = RandomNumberKind::GAUSSIAN);

template <typename Trng, typename Tvar>
void addRandomNoise(Trng *prng, Hybrid<Tvar> *x, Hybrid<Tvar> *y, Hybrid<Tvar> *z, double mult,
                    double scale = 1.0, RandomNumberKind kind = RandomNumberKind::GAUSSIAN);

template <typename Trng>
void addRandomNoise(Trng *prng, llint* x, int* x_ovrf, size_t length, double mult, double scale,
                    RandomNumberKind kind = RandomNumberKind::GAUSSIAN);

template <typename Trng>
void addRandomNoise(Trng *prng, std::vector<llint> *x, std::vector<int> *x_ovrf, size_t length,
                    double mult, double scale, RandomNumberKind kind = RandomNumberKind::GAUSSIAN);

template <typename Trng>
void addRandomNoise(Trng *prng, Hybrid<llint> *x, Hybrid<int> *x_ovrf, size_t length,
                    double mult, double scale, RandomNumberKind kind = RandomNumberKind::GAUSSIAN);

template <typename Trng>
void addRandomNoise(Trng *prng, llint* x, int* x_ovrf, llint* y, int* y_ovrf, llint* z,
                    int* z_ovrf, size_t length, double mult, double scale,
                    RandomNumberKind kind = RandomNumberKind::GAUSSIAN);

template <typename Trng>
void addRandomNoise(Trng *prng, std::vector<llint> *x, std::vector<int> *x_ovrf,
                    std::vector<llint> *y, std::vector<int> *y_ovrf, std::vector<llint> *z,
                    std::vector<int> *z_ovrf, size_t length, double mult, double scale,
                    RandomNumberKind kind = RandomNumberKind::GAUSSIAN);

template <typename Trng>
void addRandomNoise(Trng *prng, Hybrid<llint> *x, Hybrid<int> *x_ovrf, Hybrid<llint> *y,
                    Hybrid<int> *y_ovrf, Hybrid<llint> *z, Hybrid<int> *z_ovrf, size_t length,
                    double mult, double scale, RandomNumberKind kind = RandomNumberKind::GAUSSIAN);
/// \}
  
} // namespace random
} // namespace stormm

#include "random.tpp"

#endif
