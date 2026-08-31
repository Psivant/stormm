#include "copyright.h"
#include "Accelerator/gpu_enumerators.h"
#include "Trajectory/coordinate_copy.h"
#include "Trajectory/trajectory_enumerators.h"
#include "replay.h"

namespace stormm {
namespace debug {

using card::HybridFormat;
using card::HybridKind;
using card::HybridTargetLevel;
using trajectory::coordCopy;
using trajectory::CoordinateCycle;
using trajectory::TrajectoryKind;
  
//-------------------------------------------------------------------------------------------------
Replay::Replay(const PhaseSpaceSynthesis *poly_ps_in, const int depth_in, const int current_step) :
    depth{depth_in}, next_holdings_index{0},
    backtrace_order{}, holdings{},
    system_pairs{HybridKind::ARRAY, "replay_map"},
    poly_ps_ptr{const_cast<PhaseSpaceSynthesis*>(poly_ps_in)}
{
  // Allocate the holdings as copies of the original coordinate synthesis.
  holdings.reserve(depth);
  backtrace_order.resize(depth);
  for (int i = 0; i < depth; i++) {
    holdings.push_back(*poly_ps_in);
    backtrace_order[i] = 0;
  }
  incrementNextHoldingsIndex();

  // Prepare the system correspondence
  const int nsys = poly_ps_ptr->getSystemCount();
  system_pairs.resize(nsys);
  for (int i = 0; i < nsys; i++) {
    system_pairs.putHost({ i, i }, i);
  }
}

//-------------------------------------------------------------------------------------------------
int Replay::getDepth() const {
  return depth;
}

//-------------------------------------------------------------------------------------------------
int Replay::getSystemCount() const {
  return poly_ps_ptr->getSystemCount();
}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceSynthesis* Replay::getCoordinateSynthesisPointer() const {
  return poly_ps_ptr;
}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceSynthesis* Replay::getSnapshot(int relative_index) const {
  if (relative_index > 0) {
    rtErr("Unable to return a state of the system in the future (" +
          std::to_string(relative_index) + " steps ahead).", "Replay", "getSnapshot");
  }
  else {
    const int actual_index = -relative_index;
    if (actual_index < depth) {
      return (&holdings[actual_index]);
    }
    else {
      rtErr("Unable to return a state of the system " + std::to_string(actual_index) + " steps in "
            "the past.  At most " + std::to_string(depth) + " previous steps are recorded.",
            "Replay", "getSnapshot");
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void Replay::takeSnapshot(const GpuDetails &gpu) {

  // Begin by setting the coordinate time cycle
  const CoordinateCycle cpos = poly_ps_ptr->getCyclePosition();
  while (holdings[next_holdings_index].getCyclePosition() != cpos) {
    holdings[next_holdings_index].updateCyclePosition();
  }
  std::vector<HybridTargetLevel> all_tiers;
  switch (poly_ps_ptr->getFormat()) {
  case HybridFormat::HOST_ONLY:
    all_tiers.push_back(HybridTargetLevel::HOST);
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::UNIFIED:
    all_tiers.push_back(HybridTargetLevel::HOST);
    break;
  case HybridFormat::DECOUPLED:
  case HybridFormat::EXPEDITED:
    all_tiers.push_back(HybridTargetLevel::HOST);
    all_tiers.push_back(HybridTargetLevel::DEVICE);
    break;
  case HybridFormat::DEVICE_ONLY:
    all_tiers.push_back(HybridTargetLevel::DEVICE);
    break;
#endif
  }
  std::vector<TrajectoryKind> xvf = { TrajectoryKind::POSITIONS, TrajectoryKind::VELOCITIES,
                                      TrajectoryKind::FORCES };
  std::vector<CoordinateCycle> stages = { CoordinateCycle::WHITE, CoordinateCycle::BLACK };
  for (size_t i = 0; i < all_tiers.size(); i++) {
    coordCopy(&holdings[next_holdings_index], *poly_ps_ptr, system_pairs, all_tiers[i],
              all_tiers[i], gpu);
  }
}
  
//-------------------------------------------------------------------------------------------------
void Replay::incrementNextHoldingsIndex() {
  next_holdings_index += 1;
  if (next_holdings_index == depth) {
    next_holdings_index = 0;
  }
}
  
} // namespace debug
} // namespace stormm
