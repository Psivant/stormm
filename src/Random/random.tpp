// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace random {

//-------------------------------------------------------------------------------------------------
template <typename Tprng> double uniformRand(Tprng *rng, double scale) {
  return rng->uniformRandomNumber() * scale;
}

//-------------------------------------------------------------------------------------------------
template <typename Tprng>
std::vector<double> uniformRand(Tprng *rng, const size_t count, const double scale) {
  std::vector<double> result(count);
  for (size_t i = 0; i < count; i++) {
    result[i] = rng->uniformRandomNumber() * scale;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tprng>
std::vector<double> uniformRand(Tprng *rng, const size_t rows, const size_t columns,
                                const double scale, const RngFillMode mode) {
  const size_t count = rows * columns;
  std::vector<double> result(rows * columns);
  switch (mode) {
  case RngFillMode::COLUMNS:
    for (size_t i = 0; i < count; i++) {
      result[i] = rng->uniformRandomNumber() * scale;
    }
    break;
  case RngFillMode::ROWS:
    for (size_t i = 0; i < rows; i++) {
      for (size_t j = 0; j < columns; j++) {
        result[(j * rows) + i] = rng->uniformRandomNumber() * scale;
      }
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tprng, typename Tprod>
void uniformRand(Tprng *rng, std::vector<Tprod> *xv, const double scale, const double fp_scale) {
  const double eff_scale = fp_scale * scale;
  const size_t count = xv->size();
  Tprod* xv_ptr = xv->data();
  for (size_t i = 0; i < count; i++) {
    xv_ptr[i] = rng->uniformRandomNumber() * eff_scale;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tprng> float spUniformRand(Tprng *rng, float scale) {
  return rng->spUniformRandomNumber() * scale;
}

//-------------------------------------------------------------------------------------------------
template <typename Tprng>
std::vector<float> spUniformRand(Tprng *rng, const size_t count, const float scale) {
  std::vector<float> result(count);
  for (size_t i = 0; i < count; i++) {
    result[i] = rng->spUniformRandomNumber() * scale;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tprng>
std::vector<float> spUniformRand(Tprng *rng, const size_t rows, const size_t columns,
                                 const float scale, const RngFillMode mode) {
  const size_t count = rows * columns;
  std::vector<float> result(rows * columns);
  switch (mode) {
  case RngFillMode::COLUMNS:
    for (size_t i = 0; i < count; i++) {
      result[i] = rng->spUniformRandomNumber() * scale;
    }
    break;
  case RngFillMode::ROWS:
    for (size_t i = 0; i < rows; i++) {
      for (size_t j = 0; j < columns; j++) {
        result[(j * rows) + i] = rng->spUniformRandomNumber() * scale;
      }
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tprng, typename Tprod>
void spUniformRand(Tprng *rng, std::vector<Tprod> *xv, const float scale, const float fp_scale) {
  const float eff_scale = fp_scale * scale;
  const size_t count = xv->size();
  Tprod* xv_ptr = xv->data();
  for (size_t i = 0; i < count; i++) {
    xv_ptr[i] = rng->spUniformRandomNumber() * eff_scale;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tprng> double gaussianRand(Tprng *rng, double scale) {
  return rng->gaussianRandomNumber() * scale;
}

//-------------------------------------------------------------------------------------------------
template <typename Tprng>
std::vector<double> gaussianRand(Tprng *rng, const size_t count, const double scale) {
  std::vector<double> result(count);
  for (size_t i = 0; i < count; i++) {
    result[i] = rng->gaussianRandomNumber() * scale;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tprng>
std::vector<double> gaussianRand(Tprng *rng, const size_t rows, const size_t columns,
                                 const double scale, const RngFillMode mode) {
  const size_t count = rows * columns;
  std::vector<double> result(rows * columns);
  switch (mode) {
  case RngFillMode::COLUMNS:
    for (size_t i = 0; i < count; i++) {
      result[i] = rng->gaussianRandomNumber() * scale;
    }
    break;
  case RngFillMode::ROWS:
    for (size_t i = 0; i < rows; i++) {
      for (size_t j = 0; j < columns; j++) {
        result[(j * rows) + i] = rng->gaussianRandomNumber() * scale;
      }
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tprng, typename Tprod>
void gaussianRand(Tprng *rng, std::vector<Tprod> *xv, const double fp_scale, const double scale) {
  const double eff_scale = fp_scale * scale;
  const size_t count = xv->size();
  Tprod* xv_ptr = xv->data();
  for (size_t i = 0; i < count; i++) {
    xv_ptr[i] = rng->gaussianRandomNumber() * eff_scale;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tprng> float spGaussianRand(Tprng *rng, float scale) {
  return rng->spGaussianRandomNumber() * scale;
}

//-------------------------------------------------------------------------------------------------
template <typename Tprng>
std::vector<float> spGaussianRand(Tprng *rng, const size_t count, const float scale) {
  std::vector<float> result(count);
  for (size_t i = 0; i < count; i++) {
    result[i] = rng->spGaussianRandomNumber() * scale;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tprng>
std::vector<float> spGaussianRand(Tprng *rng, const size_t rows, const size_t columns,
                                  const float scale, const RngFillMode mode) {
  const size_t count = rows * columns;
  std::vector<float> result(rows * columns);
  switch (mode) {
  case RngFillMode::COLUMNS:
    for (size_t i = 0; i < count; i++) {
      result[i] = rng->spGaussianRandomNumber() * scale;
    }
    break;
  case RngFillMode::ROWS:
    for (size_t i = 0; i < rows; i++) {
      for (size_t j = 0; j < columns; j++) {
        result[(j * rows) + i] = rng->spGaussianRandomNumber() * scale;
      }
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tprng, typename Tprod>
void spGaussianRand(Tprng *rng, std::vector<Tprod> *xv, const float fp_scale, const float scale) {
  const float eff_scale = fp_scale * scale;
  const size_t count = xv->size();
  Tprod* xv_ptr = xv->data();
  for (size_t i = 0; i < count; i++) {
    xv_ptr[i] = rng->spGaussianRandomNumber() * eff_scale;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
RandomNumberMill<T>::RandomNumberMill(const ullint2 state_in, const size_t generators_in,
                                      const size_t depth_in, const RandomNumberKind init_kind,
                                      const size_t bank_limit) :
    style{RandomAlgorithm::XOROSHIRO_128P}, generators{generators_in}, depth{depth_in},
    refresh_stride{(generators + depth - 1LLU) / depth},
    state_xy{generators, "state_xy"},
    state_zw{0, "state_zw"},
    bank{HybridKind::ARRAY, "rng_bank"}
{
  checkDimensions(bank_limit);
  bank.resize(generators * depth);

  // Seed the first random number generator, and from there all the rest.
  Xoroshiro128pGenerator xrs;
  xrs.setState(state_in);
  for (size_t i = 0LLU; i < generators; i++) {
    state_xy.putHost(xrs.revealState(), i);
    xrs.longJump();
  }

  // Initialize the random number bank based on the states of all generators.
  initializeBank(init_kind);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
RandomNumberMill<T>::RandomNumberMill(const ullint4 state_in, const size_t generators_in,
                                      const size_t depth_in, const RandomNumberKind init_kind,
                                      const size_t bank_limit) :
    style{RandomAlgorithm::XOSHIRO_256PP}, generators{generators_in}, depth{depth_in},
    refresh_stride{(generators + depth - 1LLU) / depth},
    state_xy{generators, "state_xy"},
    state_zw{generators, "state_zw"},
    bank{HybridKind::ARRAY, "rng_bank"}
{
  checkDimensions(bank_limit);
  bank.resize(generators * depth);

  // Seed the first random number generator, and from there all the rest.
  Xoshiro256ppGenerator xrs;
  xrs.setState(state_in);
  for (size_t i = 0LLU; i < generators; i++) {
    const ullint4 curr_state = xrs.revealState();
    const ullint2 xy_part = { curr_state.x, curr_state.y };
    const ullint2 zw_part = { curr_state.z, curr_state.w };
    state_xy.putHost(xy_part, i);
    state_zw.putHost(zw_part, i);
    xrs.longJump();
  }

  // Initialize the random number bank based on the states of all generators.
  initializeBank(init_kind);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
RandomNumberMill<T>::RandomNumberMill(const size_t generators_in, const size_t depth_in,
                                      const RandomAlgorithm style_in,
                                      const RandomNumberKind init_kind, const int igseed_in,
                                      const int niter, const size_t bank_limit) :
    style{style_in}, generators{generators_in}, depth{depth_in},
    refresh_stride{(generators + depth - 1LLU) / depth},
    state_xy{generators, "state_xy"},
    state_zw{generators, "state_zw"},
    bank{HybridKind::ARRAY, "rng_bank"}
{
  checkDimensions(bank_limit);
  bank.resize(generators * depth);

  // Seed the first random number generator, and from there all the rest.
  switch (style) {
  case RandomAlgorithm::XOROSHIRO_128P:
    {
      Xoroshiro128pGenerator xrs(igseed_in, niter);
      for (size_t i = 0LLU; i < generators; i++) {
        state_xy.putHost(xrs.revealState(), i);
        xrs.longJump();
      }
    }
    break;
  case RandomAlgorithm::XOSHIRO_256PP:
    {
      Xoshiro256ppGenerator xrs(igseed_in, niter);
      for (size_t i = 0LLU; i < generators; i++) {
        const ullint4 curr_state = xrs.revealState();
        const ullint2 xy_part = { curr_state.x, curr_state.y };
        const ullint2 zw_part = { curr_state.z, curr_state.w };
        state_xy.putHost(xy_part, i);
        state_zw.putHost(zw_part, i);
        xrs.longJump();
      }
    }
    break;
  }

  // Initialize the random number bank based on the states of all generators.
  initializeBank(init_kind);

}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t RandomNumberMill<T>::getGeneratorCount() const {
  return generators;
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t RandomNumberMill<T>::getDepth() const {
  return depth;
}

//-------------------------------------------------------------------------------------------------
template <typename T> size_t RandomNumberMill<T>::getRefreshStride() const {
  return refresh_stride;
}

//-------------------------------------------------------------------------------------------------
template <typename T> T RandomNumberMill<T>::getBankValue(const size_t generator_index,
                                                          const size_t layer_index) const {
  if (generator_index >= generators || layer_index >= depth) {
    rtErr("Value index (" + std::to_string(generator_index) + ", " + std::to_string(layer_index) +
          ") is not valid in a series with " + std::to_string(generators) + " pseudo-random "
          "number generators and " + std::to_string(depth) + " bank depth.", "RandomNumberMill",
          "getBankValue");
  }
  return bank.readHost((layer_index * generators) + generator_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void RandomNumberMill<T>::uniformRandomNumbers(const size_t first_gen, const size_t last_gen) {
  Xoroshiro128pGenerator xrs;
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool t_is_double = (ct == double_type_index);
  T* bank_ptr = bank.data();
  for (size_t i = first_gen; i < last_gen; i++) {
    xrs.setState(state_xy.readHost(i));
    if (t_is_double) {
      for (size_t j = 0LLU; j < depth; j++) {
        bank_ptr[(j * generators) + i] = xrs.uniformRandomNumber();
      }
    }
    else {
      for (size_t j = 0LLU; j < depth; j++) {
        bank_ptr[(j * generators) + i] = xrs.spUniformRandomNumber();
      }
    }
    state_xy.putHost(xrs.revealState(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void RandomNumberMill<T>::gaussianRandomNumbers(const size_t first_gen, const size_t last_gen) {
  Xoroshiro128pGenerator xrs;
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool t_is_double = (ct == double_type_index);
  T* bank_ptr = bank.data();
  for (size_t i = first_gen; i < last_gen; i++) {
    xrs.setState(state_xy.readHost(i));
    if (t_is_double) {
      for (size_t j = 0LLU; j < depth; j++) {
        bank_ptr[(j * generators) + i] = xrs.gaussianRandomNumber();
      }
    }
    else {
      for (size_t j = 0LLU; j < depth; j++) {
        bank_ptr[(j * generators) + i] = xrs.spGaussianRandomNumber();
      }
    }
    state_xy.putHost(xrs.revealState(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void RandomNumberMill<T>::checkDimensions(const size_t bank_limit) {

  // Check that the data type is scalar
  if (isFloatingPointScalarType<T>() == false) {
    rtErr("Random number generator series can only store results in floating point scalar data "
          "types.\n", "Xoroshiro128pSeries");
  }

  // Check the bank's dimensions before allocating it, as it is a product of two numbers
  if (generators * depth == 0LLU) {
    rtErr("A random number generator series cannot be constructed with " +
          std::to_string(generators) + " generators and " + std::to_string(depth) + " depth.",
          "Xoroshiro128pSeries");
  }
  if (generators * depth * sizeof(T) > bank_limit) {
    rtErr("It is not permitted to allocate " + std::to_string(generators * depth) + " random "
          "numbers.  A maximum of " + std::to_string(bank_limit) + " numbers may be stored within "
          "the machine's stated limits.  If no keyword is available to refine this quantity, it "
          "may be necessary to recompile the code with a more permissive built-in limit.",
          "Xoroshiro128pSeries");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void RandomNumberMill<T>::initializeBank(const RandomNumberKind init_kind) {
  T* bank_ptr = bank.data();
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool t_is_double = (ct == double_type_index);
  Xoroshiro128pGenerator xrs128p;
  Xoshiro256ppGenerator xrs256pp;
  for (size_t i = 0LLU; i < generators; i++) {
    switch (style) {
    case RandomAlgorithm::XOROSHIRO_128P:
      {
        xrs128p.setState(state_xy.readHost(i));
        switch (init_kind) {
        case RandomNumberKind::UNIFORM:
          if (t_is_double) {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs128p.uniformRandomNumber();
            }
          }
          else {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs128p.spUniformRandomNumber();
            }
          }
          break;
        case RandomNumberKind::GAUSSIAN:
          if (t_is_double) {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs128p.gaussianRandomNumber();
            }
          }
          else {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs128p.spGaussianRandomNumber();
            }
          }
          break;
        }
        state_xy.putHost(xrs128p.revealState(), i);
      }
      break;
    case RandomAlgorithm::XOSHIRO_256PP:
      {
        ullint4 init_state;
        ullint2 tmp_state = state_xy.readHost(i);
        init_state.x = tmp_state.x;
        init_state.y = tmp_state.y;
        tmp_state = state_zw.readHost(i);
        init_state.z = tmp_state.x;
        init_state.w = tmp_state.y;
        xrs256pp.setState(init_state);
        switch (init_kind) {
        case RandomNumberKind::UNIFORM:
          if (t_is_double) {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs256pp.uniformRandomNumber();
            }
          }
          else {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs256pp.spUniformRandomNumber();
            }
          }
          break;
        case RandomNumberKind::GAUSSIAN:
          if (t_is_double) {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs256pp.gaussianRandomNumber();
            }
          }
          else {
            for (size_t j = 0LLU; j < depth; j++) {
              bank_ptr[(j * generators) + i] = xrs256pp.spGaussianRandomNumber();
            }
          }
          break;
        }
        init_state = xrs256pp.revealState();
        tmp_state.x = init_state.x;
        tmp_state.y = init_state.y;
        state_xy.putHost(tmp_state, i);
        tmp_state.x = init_state.z;
        tmp_state.y = init_state.w;
        state_zw.putHost(tmp_state, i);
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void fillRandomCache(ullint2* state_xy, ullint2* state_zw, T* cache, const size_t length,
                     const size_t depth, const RandomAlgorithm method,
                     const RandomNumberKind product, const size_t index_start,
                     const size_t index_end) {
  const size_t padded_atom_count = roundUp(length, warp_size_zu);
  Xoroshiro128pGenerator xrs128p;
  Xoshiro256ppGenerator xrs256pp;
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool t_is_double = (ct == double_type_index);
  for (size_t i = index_start; i < index_end; i++) {

    // Set the state of the generator object to the current value in the array of generators.
    // Next, generate a series of random numbers with this seed, then store the resulting state.
    // This is the most economical way to go about the process in terms of overall memory
    // transactions--the state is 256 bits loaded and stored, and each random number is 32 or 64
    // bits (depending on the mode), so best to store at least a handful of random numbers with
    // each checkout of the state.  On the CPU, this process is complicated by the fact that every
    // write of a random number result is going to be a cache miss.  On the GPU, memory coalescence
    // will be ideal.  Store the evolved state vector back in its original place.
    size_t pos = i;
    switch (method) {
    case RandomAlgorithm::XOROSHIRO_128P:
      {
        xrs128p.setState(state_xy[i]);
        if (t_is_double) {
          switch (product) {
          case RandomNumberKind::UNIFORM:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs128p.uniformRandomNumber();
              pos += padded_atom_count;
            }
            break;
          case RandomNumberKind::GAUSSIAN:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs128p.gaussianRandomNumber();
              pos += padded_atom_count;
            }
            break;
          }
        }
        else {
          switch (product) {
          case RandomNumberKind::UNIFORM:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs128p.spUniformRandomNumber();
              pos += padded_atom_count;
            }
            break;
          case RandomNumberKind::GAUSSIAN:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs128p.spGaussianRandomNumber();
              pos += padded_atom_count;
            }
            break;
          }
        }
        state_xy[i] = xrs128p.revealState();
      }
      break;
    case RandomAlgorithm::XOSHIRO_256PP:
      {
        xrs256pp.setState({ state_xy[i].x, state_xy[i].y, state_zw[i].x, state_zw[i].y });
        if (t_is_double) {
          switch (product) {
          case RandomNumberKind::UNIFORM:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs256pp.uniformRandomNumber();
              pos += padded_atom_count;
            }
            break;
          case RandomNumberKind::GAUSSIAN:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs256pp.gaussianRandomNumber();
              pos += padded_atom_count;
            }
            break;
          }
        }
        else {
          switch (product) {
          case RandomNumberKind::UNIFORM:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs256pp.spUniformRandomNumber();
              pos += padded_atom_count;
            }
            break;
          case RandomNumberKind::GAUSSIAN:
            for (size_t j = 0; j < depth; j++) {
              cache[pos] = xrs256pp.spGaussianRandomNumber();
              pos += padded_atom_count;
            }
            break;
          }
        }
        const ullint4 current_state = xrs256pp.revealState();
        state_xy[i] = { current_state.x, current_state.y };
        state_zw[i] = { current_state.x, current_state.y };
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void fillRandomCache(std::vector<ullint2> *state_xy, std::vector<ullint2> *state_zw,
                     std::vector<T> *cache, const size_t length, const size_t depth,
                     const RandomAlgorithm method, const RandomNumberKind product,
                     const size_t index_start, const size_t index_end) {
  fillRandomCache(state_xy->data(), state_zw->data(), cache->data(), length, depth, method,
                  product, index_start, index_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void fillRandomCache(Hybrid<ullint2> *state_xy, Hybrid<ullint2> *state_zw, Hybrid<T> *cache,
                     const size_t length, const size_t depth, const RandomAlgorithm method,
                     const RandomNumberKind product, const size_t index_start,
                     const size_t index_end) {
  fillRandomCache(state_xy->data(), state_zw->data(), cache->data(), length, depth, method,
                  product, index_start, index_end);
}

//-------------------------------------------------------------------------------------------------
template <typename Trng, typename Tvar>
void addRandomNoise(Trng *prng, Tvar* x, const size_t length, const double mult,
                    const double scale, const RandomNumberKind kind) {
  if (isFloatingPointScalarType<Tvar>()) {
    if (fabs(scale - 1.0) > 1.0e-6) {
      rtErr("Floating point data types assume that there is no scaling factor on the numerical "
            "representation.  Fold the scaling factor into the multiplier if required.",
            "addRandomNoise");
    }
    switch (kind) {
    case RandomNumberKind::UNIFORM:
      for (size_t i = 0; i < length; i++) {
        x[i] += (0.5 - prng->uniformRandomNumber()) * mult;
      }
      break;
    case RandomNumberKind::GAUSSIAN:
      for (size_t i = 0; i < length; i++) {
        x[i] += prng->gaussianRandomNumber() * mult;
      }
      break;
    }
  }
  else {
    const double relevant_mult = mult * scale;
    switch (kind) {
    case RandomNumberKind::UNIFORM:
      if (std::type_index(typeid(Tvar)).hash_code() == llint_type_index) {
        for (size_t i = 0; i < length; i++) {
          x[i] += llround((0.5 - prng->uniformRandomNumber()) * relevant_mult);
        }
      }
      else {
        for (size_t i = 0; i < length; i++) {
          x[i] += round((0.5 - prng->uniformRandomNumber()) * relevant_mult);
        }
      }
      break;
    case RandomNumberKind::GAUSSIAN:
      if (std::type_index(typeid(Tvar)).hash_code() == llint_type_index) {
        for (size_t i = 0; i < length; i++) {
          x[i] += llround((0.5 - prng->gaussianRandomNumber()) * relevant_mult);
        }
      }
      else {
        for (size_t i = 0; i < length; i++) {
          x[i] += round((0.5 - prng->gaussianRandomNumber()) * relevant_mult);
        }
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Trng, typename Tvar>
void addRandomNoise(Trng *prng, std::vector<Tvar> *x, const double mult, const double scale,
                    const RandomNumberKind kind) {
  addRandomNoise(prng, x->data(), x->size(), mult, scale, kind);
}

//-------------------------------------------------------------------------------------------------
template <typename Trng, typename Tvar>
void addRandomNoise(Trng *prng, Hybrid<Tvar> *x, const double mult, const double scale,
                    const RandomNumberKind kind) {
  addRandomNoise(prng, x->data(), x->size(), mult, scale, kind);
}

//-------------------------------------------------------------------------------------------------
template <typename Trng, typename Tvar>
void addRandomNoise(Trng *prng, Tvar* x, Tvar* y, Tvar *z, const size_t length, const double mult,
                    const double scale, const RandomNumberKind kind) {
  if (isFloatingPointScalarType<Tvar>()) {
    if (fabs(scale - 1.0) > 1.0e-6) {
      rtErr("Floating point data types assume that there is no scaling factor on the numerical "
            "representation.  Fold the scaling factor into the multiplier if required.",
            "addRandomNoise");
    }
    switch (kind) {
    case RandomNumberKind::UNIFORM:
      for (size_t i = 0; i < length; i++) {
        x[i] += (0.5 - prng->uniformRandomNumber()) * mult;
        y[i] += (0.5 - prng->uniformRandomNumber()) * mult;
        z[i] += (0.5 - prng->uniformRandomNumber()) * mult;
      }
      break;
    case RandomNumberKind::GAUSSIAN:
      for (size_t i = 0; i < length; i++) {
        x[i] += prng->gaussianRandomNumber() * mult;
        y[i] += prng->gaussianRandomNumber() * mult;
        z[i] += prng->gaussianRandomNumber() * mult;
      }
      break;
    }
  }
  else {
    const double relevant_mult = mult * scale;
    switch (kind) {
    case RandomNumberKind::UNIFORM:
      if (std::type_index(typeid(Tvar)).hash_code() == llint_type_index) {
        for (size_t i = 0; i < length; i++) {
          x[i] += llround((0.5 - prng->uniformRandomNumber()) * relevant_mult);
          y[i] += llround((0.5 - prng->uniformRandomNumber()) * relevant_mult);
          z[i] += llround((0.5 - prng->uniformRandomNumber()) * relevant_mult);
        }
      }
      else {
        for (size_t i = 0; i < length; i++) {
          x[i] += round((0.5 - prng->uniformRandomNumber()) * relevant_mult);
          y[i] += round((0.5 - prng->uniformRandomNumber()) * relevant_mult);
          z[i] += round((0.5 - prng->uniformRandomNumber()) * relevant_mult);
        }
      }
      break;
    case RandomNumberKind::GAUSSIAN:
      if (std::type_index(typeid(Tvar)).hash_code() == llint_type_index) {
        for (size_t i = 0; i < length; i++) {
          x[i] += llround((0.5 - prng->gaussianRandomNumber()) * relevant_mult);
          y[i] += llround((0.5 - prng->gaussianRandomNumber()) * relevant_mult);
          z[i] += llround((0.5 - prng->gaussianRandomNumber()) * relevant_mult);
        }
      }
      else {
        for (size_t i = 0; i < length; i++) {
          x[i] += round((0.5 - prng->gaussianRandomNumber()) * relevant_mult);
          y[i] += round((0.5 - prng->gaussianRandomNumber()) * relevant_mult);
          z[i] += round((0.5 - prng->gaussianRandomNumber()) * relevant_mult);
        }
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Trng, typename Tvar>
void addRandomNoise(Trng *prng, std::vector<Tvar> *x, std::vector<Tvar> *y, std::vector<Tvar> *z,
                    const double mult, const double scale, const RandomNumberKind kind) {
  const size_t length = x->size();
  if (length != y->size() || length != z->size()) {
    rtErr("The arrays must all be the same lengths (currently " + std::to_string(length) + ", " +
          std::to_string(y->size()) + ", " + std::to_string(z->size()) + ").", "addRandomNoise");
  }
  addRandomNoise(prng, x->data(), y->data(), z->data(), x->size(), mult, kind);
}

//-------------------------------------------------------------------------------------------------
template <typename Trng, typename Tvar>
void addRandomNoise(Trng *prng, Hybrid<Tvar> *x, Hybrid<Tvar> *y, Hybrid<Tvar> *z,
                    const double mult, const double scale, const RandomNumberKind kind) {
  const size_t length = x->size();
  if (length != y->size() || length != z->size()) {
    rtErr("The arrays must all be the same lengths (currently " + std::to_string(length) + ", " +
          std::to_string(y->size()) + ", " + std::to_string(z->size()) + ").", "addRandomNoise");
  }
  addRandomNoise(prng, x->data(), y->data(), z->data(), x->size(), mult, kind);
}

//-------------------------------------------------------------------------------------------------
template <typename Trng>
void addRandomNoise(Trng *prng, llint* x, int* x_ovrf, const size_t length, const double mult,
                    const double scale, const RandomNumberKind kind) {
  const double relevant_mult = mult * scale;
  switch (kind) {
  case RandomNumberKind::UNIFORM:
    for (size_t i = 0; i < length; i++) {
      const double noise = (0.5 - prng->uniformRandomNumber()) * relevant_mult;
      const int95_t nxi = hostInt95Sum(x[i], x_ovrf[i], noise);
      x[i] = nxi.x;
      x_ovrf[i] = nxi.y;
    }
    break;
  case RandomNumberKind::GAUSSIAN:
    for (size_t i = 0; i < length; i++) {
      const double noise = (0.5 - prng->gaussianRandomNumber()) * relevant_mult;
      const int95_t nxi = hostInt95Sum(x[i], x_ovrf[i], noise);
      x[i] = nxi.x;
      x_ovrf[i] = nxi.y;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Trng>
void addRandomNoise(Trng *prng, std::vector<llint> *x, std::vector<int> *x_ovrf, const double mult,
                    const double scale, const RandomNumberKind kind) {
  if (x->size() != x_ovrf->size()) {
    rtErr("The array of overflow bits must be the same length as the primary array of "
          "fixed-precision numbers.", "addRandomNoise");
  }
  addRandomNoise(prng, x->data(), x_ovrf->data(), x->size(), mult, scale, kind);
}

//-------------------------------------------------------------------------------------------------
template <typename Trng>
void addRandomNoise(Trng *prng, Hybrid<llint> *x, Hybrid<int> *x_ovrf, const double mult,
                    const double scale, const RandomNumberKind kind) {
  if (x->size() != x_ovrf->size()) {
    rtErr("The array of overflow bits must be the same length as the primary array of "
          "fixed-precision numbers.", "addRandomNoise");
  }
  addRandomNoise(prng, x->data(), x_ovrf->data(), x->size(), mult, scale, kind);
}

//-------------------------------------------------------------------------------------------------
template <typename Trng>
void addRandomNoise(Trng *prng, llint* x, int* x_ovrf, llint* y, int* y_ovrf, llint* z,
                    int* z_ovrf, const size_t length, const double mult, const double scale,
                    const RandomNumberKind kind) {
  const double relevant_mult = mult * scale;
  switch (kind) {
  case RandomNumberKind::UNIFORM:
    for (size_t i = 0; i < length; i++) {
      const double noise_x = (0.5 - prng->uniformRandomNumber()) * relevant_mult;
      const double noise_y = (0.5 - prng->uniformRandomNumber()) * relevant_mult;
      const double noise_z = (0.5 - prng->uniformRandomNumber()) * relevant_mult;
      const int95_t nxi = hostInt95Sum(x[i], x_ovrf[i], noise_x);
      const int95_t nyi = hostInt95Sum(y[i], y_ovrf[i], noise_y);
      const int95_t nzi = hostInt95Sum(z[i], z_ovrf[i], noise_z);
      x[i] = nxi.x;
      y[i] = nyi.x;
      z[i] = nzi.x;
      x_ovrf[i] = nxi.y;
      y_ovrf[i] = nyi.y;
      z_ovrf[i] = nzi.y;
    }
    break;
  case RandomNumberKind::GAUSSIAN:
    for (size_t i = 0; i < length; i++) {
      const double noise_x = (0.5 - prng->gaussianRandomNumber()) * relevant_mult;
      const double noise_y = (0.5 - prng->gaussianRandomNumber()) * relevant_mult;
      const double noise_z = (0.5 - prng->gaussianRandomNumber()) * relevant_mult;
      const int95_t nxi = hostInt95Sum(x[i], x_ovrf[i], noise_x);
      const int95_t nyi = hostInt95Sum(y[i], y_ovrf[i], noise_y);
      const int95_t nzi = hostInt95Sum(z[i], z_ovrf[i], noise_z);
      x[i] = nxi.x;
      y[i] = nyi.x;
      z[i] = nzi.x;
      x_ovrf[i] = nxi.y;
      y_ovrf[i] = nyi.y;
      z_ovrf[i] = nzi.y;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Trng>
void addRandomNoise(Trng *prng, std::vector<llint> *x, std::vector<int> *x_ovrf,
                    std::vector<llint> *y, std::vector<int> *y_ovrf, std::vector<llint> *z,
                    std::vector<int> *z_ovrf, const double mult, const double scale,
                    const RandomNumberKind kind) {
  const size_t length = x->size();
  if (length != y->size() || length != z->size()) {
    rtErr("The arrays must all be the same lengths (currently " + std::to_string(length) + ", " +
          std::to_string(y->size()) + ", " + std::to_string(z->size()) + ").", "addRandomNoise");
  }
  if (length != x_ovrf->size() || length != y_ovrf->size() || length != z_ovrf->size()) {
    rtErr("The arrays of overflow bits must be the same length as the primary array of "
          "fixed-precision numbers.", "addRandomNoise");
  }
  addRandomNoise(prng, x->data(), x_ovrf->data(), x->size(), mult, scale, kind);
}

//-------------------------------------------------------------------------------------------------
template <typename Trng>
void addRandomNoise(Trng *prng, Hybrid<llint> *x, Hybrid<int> *x_ovrf, Hybrid<llint> *y,
                    Hybrid<int> *y_ovrf, Hybrid<llint> *z, Hybrid<int> *z_ovrf, const double mult,
                    const double scale, const RandomNumberKind kind) {
  const size_t length = x->size();
  if (length != y->size() || length != z->size()) {
    rtErr("The arrays must all be the same lengths (currently " + std::to_string(length) + ", " +
          std::to_string(y->size()) + ", " + std::to_string(z->size()) + ").", "addRandomNoise");
  }
  if (length != x_ovrf->size() || length != y_ovrf->size() || length != z_ovrf->size()) {
    rtErr("The arrays of overflow bits must be the same length as the primary array of "
          "fixed-precision numbers.", "addRandomNoise");
  }
  addRandomNoise(prng, x->data(), x_ovrf->data(), x->size(), mult, scale, kind);
}

} // namespace random
} // namespace stormm
