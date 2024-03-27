#include "copyright.h"
#include "Math/rounding.h"
#include "Parsing/parse.h"
#include "watcher.h"

namespace stormm {
namespace review {

using card::HybridKind;
using parse::realToString;
using stmath::roundUp;

//-------------------------------------------------------------------------------------------------
Watcher::Watcher(const PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                 const float force_threshold_in, const float speed_threshold_in,
                 const bool track_momentum_purge_in, int max_reports_in,
                 const ExceptionResponse policy_in) :
    policy{policy_in},
    max_reports{max_reports_in},
    force_threshold{force_threshold_in},
    large_force_count{HybridKind::POINTER, "watcher_lf_cnt"},
    large_forces{HybridKind::ARRAY, "watcher_lf"},
    large_force_steps{HybridKind::POINTER, "watcher_lf_steps"},
    large_force_stages{HybridKind::POINTER, "watcher_lv_stages"},
    speed_threshold{speed_threshold_in},
    high_velocity_count{HybridKind::POINTER, "watcher_lv_cnt"},
    high_velocities{HybridKind::ARRAY, "watcher_lv"},
    high_velocity_steps{HybridKind::POINTER, "watcher_lv_steps"},
    high_velocity_stages{HybridKind::POINTER, "watcher_lv_stages"},
    failed_rattle_count{HybridKind::POINTER, "watcher_rattle_cnt"},
    failed_shake_count{HybridKind::POINTER, "watcher_shake_cnt"},
    rattle_group_failures{HybridKind::ARRAY, "watcher_shake"},
    shake_group_failures{HybridKind::ARRAY, "watcher_shake"},
    rattle_violations{HybridKind::POINTER, "watcher_rattle_viol"},
    shake_violations{HybridKind::POINTER, "watcher_shake_viol"},
    track_momentum_purge{track_momentum_purge_in},
    x_velocity_purge{HybridKind::POINTER, "watcher_xvel"},
    y_velocity_purge{HybridKind::POINTER, "watcher_yvel"},
    z_velocity_purge{HybridKind::POINTER, "watcher_zvel"},
    x_angular_purge{HybridKind::POINTER, "watcher_xvel"},
    y_angular_purge{HybridKind::POINTER, "watcher_yvel"},
    z_angular_purge{HybridKind::POINTER, "watcher_zvel"},
    int_data{HybridKind::ARRAY, "watcher_int"},
    float_data{HybridKind::ARRAY, "watcher_float"},
    poly_ps_ptr{const_cast<PhaseSpaceSynthesis*>(poly_ps)},
    poly_ag_ptr{const_cast<AtomGraphSynthesis*>(poly_ag.getSelfPointer())}
{
  allocate();
}

//-------------------------------------------------------------------------------------------------
Watcher::Watcher(const PhaseSpaceSynthesis &poly_ps, const AtomGraphSynthesis &poly_ag,
                 const float force_threshold_in, const float speed_threshold_in,
                 const bool track_momentum_purge_in, int max_reports_in,
                 const ExceptionResponse policy_in) :
    Watcher(poly_ps.getSelfPointer(), poly_ag, force_threshold_in, speed_threshold_in,
            track_momentum_purge_in, max_reports_in, policy_in)
{}

//-------------------------------------------------------------------------------------------------
Watcher::Watcher(const Watcher &original) :
    policy{original.policy},
    max_reports{original.max_reports},
    force_threshold{original.force_threshold},
    large_force_count{original.large_force_count},
    large_forces{original.large_forces},
    large_force_steps{original.large_force_steps},
    large_force_stages{original.large_force_stages},
    speed_threshold{original.speed_threshold},
    high_velocity_count{original.high_velocity_count},
    high_velocities{original.high_velocities},
    high_velocity_steps{original.high_velocity_steps},
    high_velocity_stages{original.high_velocity_stages},
    failed_rattle_count{original.failed_rattle_count},
    failed_shake_count{original.failed_shake_count},
    rattle_group_failures{original.rattle_group_failures},
    shake_group_failures{original.shake_group_failures},
    rattle_violations{original.rattle_violations},
    shake_violations{original.shake_violations},
    track_momentum_purge{original.track_momentum_purge},
    x_velocity_purge{original.x_velocity_purge},
    y_velocity_purge{original.y_velocity_purge},
    z_velocity_purge{original.z_velocity_purge},
    x_angular_purge{original.x_angular_purge},
    y_angular_purge{original.y_angular_purge},
    z_angular_purge{original.z_angular_purge},
    int_data{original.int_data},
    float_data{original.float_data},
    poly_ps_ptr{original.poly_ps_ptr},
    poly_ag_ptr{original.poly_ag_ptr}
{
  allocate();
}

//-------------------------------------------------------------------------------------------------
Watcher::Watcher(Watcher &&original) :
    policy{original.policy},
    max_reports{original.max_reports},
    force_threshold{original.force_threshold},
    large_force_count{std::move(original.large_force_count)},
    large_forces{std::move(original.large_forces)},
    large_force_steps{std::move(original.large_force_steps)},
    large_force_stages{std::move(original.large_force_stages)},
    speed_threshold{original.speed_threshold},
    high_velocity_count{std::move(original.high_velocity_count)},
    high_velocities{std::move(original.high_velocities)},
    high_velocity_steps{std::move(original.high_velocity_steps)},
    high_velocity_stages{std::move(original.high_velocity_stages)},
    failed_rattle_count{std::move(original.failed_rattle_count)},
    failed_shake_count{std::move(original.failed_shake_count)},
    rattle_group_failures{std::move(original.rattle_group_failures)},
    shake_group_failures{std::move(original.shake_group_failures)},
    rattle_violations{std::move(original.rattle_violations)},
    shake_violations{std::move(original.shake_violations)},
    track_momentum_purge{original.track_momentum_purge},
    x_velocity_purge{std::move(original.x_velocity_purge)},
    y_velocity_purge{std::move(original.y_velocity_purge)},
    z_velocity_purge{std::move(original.z_velocity_purge)},
    x_angular_purge{std::move(original.x_angular_purge)},
    y_angular_purge{std::move(original.y_angular_purge)},
    z_angular_purge{std::move(original.z_angular_purge)},
    int_data{std::move(original.int_data)},
    float_data{std::move(original.float_data)},
    poly_ps_ptr{original.poly_ps_ptr},
    poly_ag_ptr{original.poly_ag_ptr}
{}

//-------------------------------------------------------------------------------------------------
Watcher& Watcher::operator=(const Watcher &other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }
  policy = other.policy;
  max_reports = other.max_reports;
  force_threshold = other.force_threshold;
  large_force_count = other.large_force_count;
  large_forces = other.large_forces;
  large_force_steps = other.large_force_steps;
  large_force_stages = other.large_force_stages;
  speed_threshold = other.speed_threshold;
  high_velocity_count = other.high_velocity_count;
  high_velocities = other.high_velocities;
  high_velocity_steps = other.high_velocity_steps;
  high_velocity_stages = other.high_velocity_stages;
  failed_rattle_count = other.failed_rattle_count;
  failed_shake_count = other.failed_shake_count;
  rattle_group_failures = other.rattle_group_failures;
  shake_group_failures = other.shake_group_failures;
  rattle_violations = other.rattle_violations;
  shake_violations = other.shake_violations;
  track_momentum_purge = other.track_momentum_purge;
  x_velocity_purge = other.x_velocity_purge;
  y_velocity_purge = other.y_velocity_purge;
  z_velocity_purge = other.z_velocity_purge;
  x_angular_purge = other.x_angular_purge;
  y_angular_purge = other.y_angular_purge;
  z_angular_purge = other.z_angular_purge;
  int_data = other.int_data;
  float_data = other.float_data;
  poly_ps_ptr = other.poly_ps_ptr;
  poly_ag_ptr = other.poly_ag_ptr;

  // Rebase pointers and return the result
  allocate();
  return *this;
}

//-------------------------------------------------------------------------------------------------
Watcher& Watcher::operator=(Watcher &&other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }
  policy = other.policy;
  max_reports = other.max_reports;
  force_threshold = other.force_threshold;
  large_force_count = std::move(other.large_force_count);
  large_forces = std::move(other.large_forces);
  large_force_steps = std::move(other.large_force_steps);
  large_force_stages = std::move(other.large_force_stages);
  speed_threshold = other.speed_threshold;
  high_velocity_count = std::move(other.high_velocity_count);
  high_velocities = std::move(other.high_velocities);
  high_velocity_steps = std::move(other.high_velocity_steps);
  high_velocity_stages = std::move(other.high_velocity_stages);
  failed_rattle_count = std::move(other.failed_rattle_count);
  failed_shake_count = std::move(other.failed_shake_count);
  rattle_group_failures = std::move(other.rattle_group_failures);
  shake_group_failures = std::move(other.shake_group_failures);
  rattle_violations = std::move(other.rattle_violations);
  shake_violations = std::move(other.shake_violations);
  track_momentum_purge = other.track_momentum_purge;
  x_velocity_purge = std::move(other.x_velocity_purge);
  y_velocity_purge = std::move(other.y_velocity_purge);
  z_velocity_purge = std::move(other.z_velocity_purge);
  x_angular_purge = std::move(other.x_angular_purge);
  y_angular_purge = std::move(other.y_angular_purge);
  z_angular_purge = std::move(other.z_angular_purge);
  int_data = std::move(other.int_data);
  float_data = std::move(other.float_data);
  poly_ps_ptr = other.poly_ps_ptr;
  poly_ag_ptr = other.poly_ag_ptr;
  return *this;
}

//-------------------------------------------------------------------------------------------------
int Watcher::getReportCount() const {
  return max_reports;
}

//-------------------------------------------------------------------------------------------------
int Watcher::getSystemCount() const {
  return poly_ps_ptr->getSystemCount();
}

//-------------------------------------------------------------------------------------------------
float Watcher::getForceThreshold() const {
  return force_threshold;
}

//-------------------------------------------------------------------------------------------------
int Watcher::getFailedRattleCount(const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return failed_rattle_count.readHost(0);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return failed_rattle_count.readDevice(0);
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int Watcher::getFailedShakeCount(const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return failed_shake_count.readHost(0);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return failed_shake_count.readDevice(0);
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<uint2> Watcher::getRattleFailures(const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return rattle_group_failures.readHost();
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return rattle_group_failures.readDevice();
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<uint2> Watcher::getShakeFailures(const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return shake_group_failures.readHost();
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return shake_group_failures.readDevice();
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<float> Watcher::getRattleViolations(const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return rattle_violations.readHost();
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return rattle_violations.readDevice();
#endif
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
std::vector<float> Watcher::getShakeViolations(const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return shake_violations.readHost();
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return shake_violations.readDevice();
#endif
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
void Watcher::setForceThreshold(const float force_threshold_in) {
  force_threshold = force_threshold_in;
  validateForceThreshold();
}

//-------------------------------------------------------------------------------------------------
void Watcher::setSpeedThreshold(const float speed_threshold_in) {
  speed_threshold = speed_threshold_in;
  validateSpeedThreshold();
}

//-------------------------------------------------------------------------------------------------
void Watcher::allocate() {
  
  // Allocate integer data
  const size_t padded_max_reports = roundUp(max_reports, warp_size_int);
  int_data.resize(4 + (4 * padded_max_reports));
  large_force_count.setPointer(&int_data, 0, 1);
  high_velocity_count.setPointer(&int_data, 1, 1);
  failed_rattle_count.setPointer(&int_data, 2, 1);
  failed_shake_count.setPointer(&int_data, 3, 1);
  large_force_steps.setPointer(&int_data,                              4, max_reports);
  large_force_stages.setPointer(&int_data,        4 + padded_max_reports, max_reports);
  high_velocity_steps.setPointer(&int_data, 4 + (2 * padded_max_reports), max_reports);
  high_velocity_stages.setPointer(&int_data, 4 + (3 * padded_max_reports), max_reports);

  // Allocate float data
  if (track_momentum_purge) {
    float_data.resize(8 * padded_max_reports);
    x_velocity_purge.setPointer(&float_data,                      0, max_reports);
    y_velocity_purge.setPointer(&float_data,     padded_max_reports, max_reports);
    z_velocity_purge.setPointer(&float_data, 2 * padded_max_reports, max_reports);
    x_angular_purge.setPointer(&float_data,  3 * padded_max_reports, max_reports);
    y_angular_purge.setPointer(&float_data,  4 * padded_max_reports, max_reports);
    z_angular_purge.setPointer(&float_data,  5 * padded_max_reports, max_reports);
  }
  else {
    float_data.resize(2 * padded_max_reports);
  }
  rattle_violations.setPointer(&float_data,                  0, max_reports);
  rattle_violations.setPointer(&float_data, padded_max_reports, max_reports);

  // Allocate specific arrays
  large_forces.resize(max_reports);
  high_velocities.resize(max_reports);
  rattle_group_failures.resize(max_reports);
  shake_group_failures.resize(max_reports);
}

//-------------------------------------------------------------------------------------------------
void Watcher::validateForceThreshold() const {
  if (force_threshold < minimum_force_threshold) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A force reporting threshold of " + realToString(force_threshold) + " will result in "
            "far too many forces being reported.", "Watcher", "validateForceThreshold");
    case ExceptionResponse::WARN:
      rtWarn("A force reporting threshold of " + realToString(force_threshold) + " may result in "
             "far too many forces being reported.", "Watcher", "validateForceThreshold");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void Watcher::validateSpeedThreshold() const {
  if (speed_threshold < minimum_speed_threshold) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A particle speed reporting threshold of " + realToString(speed_threshold) + " will "
            "result in far too many speeds being reported.", "Watcher", "validateSpeedThreshold");
    case ExceptionResponse::WARN:
      rtWarn("A particle speed reporting threshold of " + realToString(speed_threshold) + " may "
             "result in far too many speeds being reported.", "Watcher", "validateSpeedThreshold");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

} // namespace review
} // namespace stormm
