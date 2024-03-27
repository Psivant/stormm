#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "Math/rounding.h"
#include "implicit_solvent_workspace.h"

namespace stormm {
namespace synthesis {

using card::HybridKind;
using stmath::roundUp;

//-------------------------------------------------------------------------------------------------
ImplicitSolventWorkspace::ImplicitSolventWorkspace(const Hybrid<int> &atom_starts,
                                                   const Hybrid<int> &atom_counts,
                                                   const int bit_count) :
    fp_bits{bit_count},
    padded_atom_count{atom_starts.readHost(atom_starts.size() - 1LLU) +
                      roundUp(atom_counts.readHost(atom_starts.size() - 1LLU), warp_size_int)},
    cycle_position{CoordinateCycle::WHITE},
    psi{HybridKind::POINTER, "isw_eff_gbrad"},
    psi_overflow{HybridKind::POINTER, "isw_eff_gbrad_ovrf"},
    alt_psi{HybridKind::POINTER, "alt_eff_gbrad"},
    alt_psi_overflow{HybridKind::POINTER, "alt_eff_gbrad_ovrf"},
    sum_deijda{HybridKind::POINTER, "isw_sumdeijda"},
    sum_deijda_overflow{HybridKind::POINTER, "isw_sumdeijda_ovrf"},
    alt_sum_deijda{HybridKind::POINTER, "alt_sumdeijda"},
    alt_sum_deijda_overflow{HybridKind::POINTER, "alt_sumdeijda_ovrf"},
    llint_data{4LLU * static_cast<size_t>(padded_atom_count), "isw_llint_data"},
    int_data{4LLU * static_cast<size_t>(padded_atom_count), "isw_int_data"}
{
  psi.setPointer(&llint_data,                                       0LLU, padded_atom_count);
  alt_psi.setPointer(&llint_data,                      padded_atom_count, padded_atom_count);
  sum_deijda.setPointer(&llint_data,            2LLU * padded_atom_count, padded_atom_count);
  alt_sum_deijda.setPointer(&llint_data,        3LLU * padded_atom_count, padded_atom_count);
  psi_overflow.setPointer(&int_data,                                0LLU, padded_atom_count);
  alt_psi_overflow.setPointer(&int_data,               padded_atom_count, padded_atom_count);
  sum_deijda_overflow.setPointer(&int_data,     2LLU * padded_atom_count, padded_atom_count);
  alt_sum_deijda_overflow.setPointer(&int_data, 3LLU * padded_atom_count, padded_atom_count);
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventWorkspace::ImplicitSolventWorkspace(const Hybrid<int> &atom_starts,
                                                   const Hybrid<int> &atom_counts,
                                                   const PrecisionModel prec) :
  ImplicitSolventWorkspace(atom_starts, atom_counts, (prec == PrecisionModel::SINGLE) ? 24 : 55)
{}

//-------------------------------------------------------------------------------------------------
ImplicitSolventWorkspace::ImplicitSolventWorkspace(const ImplicitSolventWorkspace &original) :
    fp_bits{original.fp_bits},
    padded_atom_count{original.padded_atom_count},
    cycle_position{original.cycle_position},
    psi{original.psi},
    psi_overflow{original.psi_overflow},
    alt_psi{original.alt_psi},
    alt_psi_overflow{original.alt_psi_overflow},
    sum_deijda{original.sum_deijda},
    sum_deijda_overflow{original.sum_deijda_overflow},
    alt_sum_deijda{original.alt_sum_deijda},
    alt_sum_deijda_overflow{original.alt_sum_deijda_overflow},
    llint_data{original.llint_data},
    int_data{original.int_data}
{
  rebasePointers();
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventWorkspace::ImplicitSolventWorkspace(ImplicitSolventWorkspace &&original) :
    fp_bits{original.fp_bits},
    padded_atom_count{original.padded_atom_count},
    cycle_position{original.cycle_position},
    psi{std::move(original.psi)},
    psi_overflow{std::move(original.psi_overflow)},
    alt_psi{std::move(original.alt_psi)},
    alt_psi_overflow{std::move(original.alt_psi_overflow)},
    sum_deijda{std::move(original.sum_deijda)},
    sum_deijda_overflow{std::move(original.sum_deijda_overflow)},
    alt_sum_deijda{std::move(original.alt_sum_deijda)},
    alt_sum_deijda_overflow{std::move(original.alt_sum_deijda_overflow)},
    llint_data{original.llint_data},
    int_data{original.int_data}
{}

//-------------------------------------------------------------------------------------------------
ImplicitSolventWorkspace&
ImplicitSolventWorkspace::operator=(const ImplicitSolventWorkspace &other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // Perform the customary copying operations
  fp_bits = other.fp_bits;
  padded_atom_count = other.padded_atom_count;
  cycle_position = other.cycle_position;
  psi = other.psi;
  psi_overflow = other.psi_overflow;
  alt_psi = other.alt_psi;
  alt_psi_overflow = other.alt_psi_overflow;
  sum_deijda = other.sum_deijda;
  sum_deijda_overflow = other.sum_deijda_overflow;
  alt_sum_deijda = other.alt_sum_deijda;
  alt_sum_deijda_overflow = other.alt_sum_deijda_overflow;
  llint_data = other.llint_data;
  int_data = other.int_data;

  // Repair the pointers and return the result
  rebasePointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventWorkspace&
ImplicitSolventWorkspace::operator=(ImplicitSolventWorkspace &&other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // Perform the customary move operations
  fp_bits = other.fp_bits;
  padded_atom_count = other.padded_atom_count;
  cycle_position = other.cycle_position;
  psi = std::move(other.psi);
  psi_overflow = std::move(other.psi_overflow);
  alt_psi = std::move(other.alt_psi);
  alt_psi_overflow = std::move(other.alt_psi_overflow);
  sum_deijda = std::move(other.sum_deijda);
  sum_deijda_overflow = std::move(other.sum_deijda_overflow);
  alt_sum_deijda = std::move(other.alt_sum_deijda);
  alt_sum_deijda_overflow = std::move(other.alt_sum_deijda_overflow);
  llint_data = std::move(other.llint_data);
  int_data = std::move(other.int_data);
  return *this;
}

//-------------------------------------------------------------------------------------------------
void ImplicitSolventWorkspace::rebasePointers() {
  psi.swapTarget(&llint_data);
  psi_overflow.swapTarget(&int_data);
  alt_psi.swapTarget(&llint_data);
  alt_psi_overflow.swapTarget(&int_data);
  sum_deijda.swapTarget(&llint_data);
  sum_deijda_overflow.swapTarget(&int_data);
  alt_sum_deijda.swapTarget(&llint_data);
  alt_sum_deijda_overflow.swapTarget(&int_data);
}

//-------------------------------------------------------------------------------------------------
int ImplicitSolventWorkspace::getFixedPrecisionBits() const {
  return fp_bits;
}

//-------------------------------------------------------------------------------------------------
CoordinateCycle ImplicitSolventWorkspace::getCyclePosition() const {
  return cycle_position;
}

//-------------------------------------------------------------------------------------------------
ISWorkspaceKit<double> ImplicitSolventWorkspace::dpData(const CoordinateCycle orientation,
                                                        const HybridTargetLevel tier) {
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return ISWorkspaceKit<double>(fp_bits, alt_psi.data(tier), alt_psi_overflow.data(tier),
                                  alt_sum_deijda.data(tier), alt_sum_deijda_overflow.data(tier),
                                  psi.data(tier), psi_overflow.data(tier),
                                  sum_deijda.data(tier), sum_deijda_overflow.data(tier));
  case CoordinateCycle::WHITE:
    return ISWorkspaceKit<double>(fp_bits, psi.data(tier), psi_overflow.data(tier),
                                  sum_deijda.data(tier), sum_deijda_overflow.data(tier),
                                  alt_psi.data(tier), alt_psi_overflow.data(tier),
                                  alt_sum_deijda.data(tier), alt_sum_deijda_overflow.data(tier));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
ISWorkspaceKit<double> ImplicitSolventWorkspace::dpData(const HybridTargetLevel tier) {
  return dpData(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
ISWorkspaceKit<float> ImplicitSolventWorkspace::spData(const CoordinateCycle orientation,
                                                       const HybridTargetLevel tier) {
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return ISWorkspaceKit<float>(fp_bits, alt_psi.data(tier), alt_psi_overflow.data(tier),
                                 alt_sum_deijda.data(tier), alt_sum_deijda_overflow.data(tier),
                                 psi.data(tier), psi_overflow.data(tier),
                                 sum_deijda.data(tier), sum_deijda_overflow.data(tier));
  case CoordinateCycle::WHITE:
    return ISWorkspaceKit<float>(fp_bits, psi.data(tier), psi_overflow.data(tier),
                                 sum_deijda.data(tier), sum_deijda_overflow.data(tier),
                                 alt_psi.data(tier), alt_psi_overflow.data(tier),
                                 alt_sum_deijda.data(tier), alt_sum_deijda_overflow.data(tier));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
ISWorkspaceKit<float> ImplicitSolventWorkspace::spData(const HybridTargetLevel tier) {
  return spData(cycle_position, tier);
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void ImplicitSolventWorkspace::upload() {
  llint_data.upload();
  int_data.upload();
}

//-------------------------------------------------------------------------------------------------
void ImplicitSolventWorkspace::download() {
  llint_data.download();
  int_data.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void ImplicitSolventWorkspace::updateCyclePosition() {
  switch (cycle_position) {
  case CoordinateCycle::WHITE:
    cycle_position = CoordinateCycle::BLACK;
    break;
  case CoordinateCycle::BLACK:
    cycle_position = CoordinateCycle::WHITE;
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void ImplicitSolventWorkspace::initialize(const HybridTargetLevel tier,
                                          const CoordinateCycle orientation,
                                          const GpuDetails &gpu) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      llint* psi_ptr;
      llint* sum_deijda_ptr;
      int* psi_ovrf_ptr;
      int* sum_deijda_ovrf_ptr;
      switch (orientation) {
      case CoordinateCycle::BLACK:
        psi_ptr = alt_psi.data();
        psi_ovrf_ptr = alt_psi_overflow.data();
        sum_deijda_ptr = alt_sum_deijda.data();
        sum_deijda_ovrf_ptr = alt_sum_deijda_overflow.data();
        break;
      case CoordinateCycle::WHITE:
        psi_ptr = psi.data();
        psi_ovrf_ptr = psi_overflow.data();
        sum_deijda_ptr = sum_deijda.data();
        sum_deijda_ovrf_ptr = sum_deijda_overflow.data();
        break;
      }
      for (int i = 0; i < padded_atom_count; i++) {
        psi_ptr[i] = 0LL;
        psi_ovrf_ptr[i] = 0;
        sum_deijda_ptr[i] = 0LL;
        sum_deijda_ovrf_ptr[i] = 0;
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    launchInitialization(gpu, orientation);
    break;
#endif
  }
}

} // namespace synthesis
} // namespace stormm
