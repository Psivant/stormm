#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "Math/rounding.h"
#include "Reporting/error_format.h"
#include "line_minimization.h"

namespace stormm {
namespace mm {

using card::HybridKind;
using stmath::roundUp;

//-------------------------------------------------------------------------------------------------
LinMinWriter::LinMinWriter(const int nsys_in, double* l_move_in, double* s_move_in,
                           double* mfac_a_in, double* mfac_b_in, double* mfac_c_in,
                           double* nrg_a_in, double* nrg_b_in, double* nrg_c_in,
                           double* nrg_d_in) :
    nsys{nsys_in}, l_move{l_move_in}, s_move{s_move_in}, mfac_a{mfac_a_in}, mfac_b{mfac_b_in},
    mfac_c{mfac_c_in}, nrg_a{nrg_a_in}, nrg_b{nrg_b_in}, nrg_c{nrg_c_in}, nrg_d{nrg_d_in}
{}

//-------------------------------------------------------------------------------------------------
LinMinReader::LinMinReader(const int nsys_in, const double* l_move_in, const double* s_move_in,
                           const double* mfac_a_in, const double* mfac_b_in,
                           const double* mfac_c_in, const double* nrg_a_in, const double* nrg_b_in,
                           const double* nrg_c_in, const double* nrg_d_in) :
    nsys{nsys_in}, l_move{l_move_in}, s_move{s_move_in}, mfac_a{mfac_a_in}, mfac_b{mfac_b_in},
    mfac_c{mfac_c_in}, nrg_a{nrg_a_in}, nrg_b{nrg_b_in}, nrg_c{nrg_c_in}, nrg_d{nrg_d_in}
{}

//-------------------------------------------------------------------------------------------------
LineMinimization::LineMinimization(const int system_count_in, const double dx0) :
    system_count{system_count_in},
    move_length{HybridKind::POINTER, "linmin_mlen"},
    save_length{HybridKind::POINTER, "linmin_slen"},
    move_factor_a{HybridKind::POINTER, "linmin_mfac_a"},
    move_factor_b{HybridKind::POINTER, "linmin_mfac_b"},
    move_factor_c{HybridKind::POINTER, "linmin_mfac_c"},
    energy_a{HybridKind::POINTER, "linmin_nrg_a"},
    energy_b{HybridKind::POINTER, "linmin_nrg_b"},
    energy_c{HybridKind::POINTER, "linmin_nrg_c"},
    energy_d{HybridKind::POINTER, "linmin_nrg_d"},
    storage{static_cast<size_t>(roundUp(system_count, warp_size_int) * 9), "linmin_storage"}
{
  const size_t padded_nsys = roundUp(system_count, warp_size_int);
  move_length.setPointer(&storage,                 0LLU, system_count);
  save_length.setPointer(&storage,          padded_nsys, system_count);
  move_factor_a.setPointer(&storage, 2LLU * padded_nsys, system_count);
  move_factor_b.setPointer(&storage, 3LLU * padded_nsys, system_count);
  move_factor_c.setPointer(&storage, 4LLU * padded_nsys, system_count);
  energy_a.setPointer(&storage, 5LLU * padded_nsys, system_count);
  energy_b.setPointer(&storage, 6LLU * padded_nsys, system_count);
  energy_c.setPointer(&storage, 7LLU * padded_nsys, system_count);
  energy_d.setPointer(&storage, 8LLU * padded_nsys, system_count);
  primeMoveLengths(dx0);
}

//-------------------------------------------------------------------------------------------------
LineMinimization::LineMinimization(const LineMinimization &original) :
    system_count{original.system_count},
    move_length{original.move_length},
    save_length{original.save_length},
    move_factor_a{original.move_factor_a},
    move_factor_b{original.move_factor_b},
    move_factor_c{original.move_factor_c},
    energy_a{original.energy_a},
    energy_b{original.energy_b},
    energy_c{original.energy_c},
    energy_d{original.energy_d},
    storage{original.storage}
{
  // Repair pointers
  move_length.swapTarget(&storage);
  save_length.swapTarget(&storage);
  move_factor_a.swapTarget(&storage);
  move_factor_b.swapTarget(&storage);
  move_factor_c.swapTarget(&storage);
  energy_a.swapTarget(&storage);
  energy_b.swapTarget(&storage);
  energy_c.swapTarget(&storage);
  energy_d.swapTarget(&storage);
}

//-------------------------------------------------------------------------------------------------
LineMinimization& LineMinimization::operator=(const LineMinimization &other) {
  if (this == &other) {
    return *this;
  }
  system_count = other.system_count;
  move_length = other.move_length;
  save_length = other.save_length;
  move_factor_a = other.move_factor_a;
  move_factor_b = other.move_factor_b;
  move_factor_c = other.move_factor_c;
  energy_a = other.energy_a;
  energy_b = other.energy_b;
  energy_c = other.energy_c;
  energy_d = other.energy_d;
  storage = other.storage;

  // Repair pointers
  move_length.swapTarget(&storage);
  save_length.swapTarget(&storage);
  move_factor_a.swapTarget(&storage);
  move_factor_b.swapTarget(&storage);
  move_factor_c.swapTarget(&storage);
  energy_a.swapTarget(&storage);
  energy_b.swapTarget(&storage);
  energy_c.swapTarget(&storage);
  energy_d.swapTarget(&storage);

  return *this;
}

//-------------------------------------------------------------------------------------------------
int LineMinimization::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> LineMinimization::getMoveLength(const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return move_length.readHost();
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return move_length.readDevice();
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double LineMinimization::getMoveLength(const int system_index,
                                       const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return move_length.readHost(system_index);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return move_length.readDevice(system_index);
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> LineMinimization::getMoveFactor(const int move_index,
                                                    const HybridTargetLevel tier) const {
  bool problem = false;
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (move_index) {
    case 0:
      return move_factor_a.readHost();
    case 1:
      return move_factor_b.readHost();
    case 2:
      return move_factor_c.readHost();
    default:
      problem = true;
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (move_index) {
    case 0:
      return move_factor_a.readDevice();
    case 1:
      return move_factor_b.readDevice();
    case 2:
      return move_factor_c.readDevice();
    default:
      problem = true;
      break;
    }
    break;
#endif
  }
  if (problem) {    
    rtErr("Up to three moves are taken along any computed gradient during a line minimization, "
          "prior to solving a cubic polynomial to determine the optimal move.  A move index of " +
          std::to_string(move_index) + " is invalid.", "LineMinimization", "getMoveLength");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double LineMinimization::getMoveFactor(const int move_index, const int system_index,
                                       const HybridTargetLevel tier) const {
  bool problem = false;
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (move_index) {
    case 0:
      return move_factor_a.readHost(system_index);
    case 1:
      return move_factor_b.readHost(system_index);
    case 2:
      return move_factor_c.readHost(system_index);
    default:
      problem = true;
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (move_index) {
    case 0:
      return move_factor_a.readDevice(system_index);
    case 1:
      return move_factor_b.readDevice(system_index);
    case 2:
      return move_factor_c.readDevice(system_index);
    default:
      problem = true;
      break;
    }
    break;
#endif
  }
  if (problem) {    
    rtErr("Up to three moves are taken along any computed gradient during a line minimization, "
          "prior to solving a cubic polynomial to determine the optimal move.  A move index of " +
          std::to_string(move_index) + " is invalid.", "LineMinimization", "getMoveLength");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> LineMinimization::getEnergy(const int move_index,
                                                const HybridTargetLevel tier) const {
  bool problem = false;
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (move_index) {
    case 0:
      return energy_a.readHost();
    case 1:
      return energy_b.readHost();
    case 2:
      return energy_c.readHost();
    case 3:
      return energy_d.readHost();
    default:
      problem = true;
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (move_index) {
    case 0:
      return energy_a.readDevice();
    case 1:
      return energy_b.readDevice();
    case 2:
      return energy_c.readDevice();
    case 3:
      return energy_d.readDevice();
    default:
      problem = true;
      break;
    }
    break;
#endif
  }
  if (problem) {    
    rtErr("Up to three moves are taken along any computed gradient during a line minimization, "
          "prior to solving a cubic polynomial to determine the optimal move.  A move index of " +
          std::to_string(move_index) + " is invalid.", "LineMinimization", "getEnergy");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double LineMinimization::getEnergy(const int move_index,
                                   const int system_index, const HybridTargetLevel tier) const {
  bool problem = false;
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (move_index) {
    case 0:
      return energy_a.readHost(system_index);
    case 1:
      return energy_b.readHost(system_index);
    case 2:
      return energy_c.readHost(system_index);
    case 3:
      return energy_d.readHost(system_index);
    default:
      problem = true;
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (move_index) {
    case 0:
      return energy_a.readDevice(system_index);
    case 1:
      return energy_b.readDevice(system_index);
    case 2:
      return energy_c.readDevice(system_index);
    case 3:
      return energy_d.readDevice(system_index);
    default:
      problem = true;
      break;
    }
    break;
#endif
  }
  if (problem) {
    rtErr("Up to three moves are taken along any computed gradient during a line minimization, "
          "prior to solving a cubic polynomial to determine the optimal move.  A move index of " +
          std::to_string(move_index) + " is invalid.", "LineMinimization", "getEnergy");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
LinMinReader LineMinimization::data(const HybridTargetLevel tier) const {
  return LinMinReader(system_count, move_length.data(tier), save_length.data(tier),
                      move_factor_a.data(tier), move_factor_b.data(tier), move_factor_c.data(tier),
                      energy_a.data(tier), energy_b.data(tier), energy_c.data(tier),
                      energy_d.data(tier));
}

//-------------------------------------------------------------------------------------------------
LinMinWriter LineMinimization::data(const HybridTargetLevel tier) {
  return LinMinWriter(system_count, move_length.data(tier), save_length.data(tier),
                      move_factor_a.data(tier), move_factor_b.data(tier), move_factor_c.data(tier),
                      energy_a.data(tier), energy_b.data(tier), energy_c.data(tier),
                      energy_d.data(tier));
}

//-------------------------------------------------------------------------------------------------
void LineMinimization::primeMoveLengths(const double dx0) {
  double* l_move_ptr = move_length.data();
  for (int i = 0; i < system_count; i++) {
    l_move_ptr[i] = dx0;
  }
#ifdef STORMM_USE_HPC
  move_length.upload();
#endif
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void LineMinimization::upload() {
  storage.upload();
}

//-------------------------------------------------------------------------------------------------
void LineMinimization::download() {
  storage.download();
}
#endif

} // namespace mm
} // namespace stormm
