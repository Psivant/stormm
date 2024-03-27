// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace synthesis {
  
//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> AtomGraphSynthesis::getPartialCharges(const HybridTargetLevel tier) const {
  const int high_pos = atom_offsets.readHost(system_count - 1) +
                       atom_counts.readHost(system_count - 1);
  return getRealParameters<T>(atomic_charges, sp_atomic_charges, tier, 0, high_pos,
                              "getPartialCharges");
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> AtomGraphSynthesis::getAtomicMasses(const HybridTargetLevel tier) const {
  const int high_pos = atom_offsets.readHost(system_count - 1) +
                       atom_counts.readHost(system_count - 1);
  return getRealParameters<T>(atomic_masses, sp_atomic_masses, tier, 0, high_pos,
                              "getAtomicMasses");
}

} // namespace synthesis
} // namespace stormm
