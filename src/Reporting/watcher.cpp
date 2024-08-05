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
WatcherWriter::WatcherWriter(const int nsystem_in, const int max_reports_in,
                             const float force_limit_in, const float speed_limit_in,
                             const bool track_purge_in, int* nforce_in, int* nspeed_in,
                             int* nrattle_in, int* nshake_in, float4* forces_in,
                             int* force_steps_in, int* force_stages_in, float4* speeds_in,
                             int* speed_steps_in, int* speed_stages_in, uint2* rattle_fails_in,
                             uint2* shake_fails_in, float* rattle_ext_in, float* shake_ext_in,
                             float* xvel_purge_in, float* yvel_purge_in, float* zvel_purge_in,
                             float* xang_purge_in, float* yang_purge_in, float* zang_purge_in) :
    nsystem{nsystem_in}, max_reports{max_reports_in}, force_limit{force_limit_in},
    speed_limit{speed_limit_in}, track_purge{track_purge_in}, nforce{nforce_in},
    nspeed{nspeed_in}, nrattle{nrattle_in}, nshake{nshake_in}, forces{forces_in},
    force_steps{force_steps_in}, force_stages{force_stages_in}, speeds{speeds_in},
    speed_steps{speed_steps_in}, speed_stages{speed_stages_in}, rattle_fails{rattle_fails_in},
    shake_fails{shake_fails_in}, rattle_ext{rattle_ext_in}, shake_ext{shake_ext_in},
    xvel_purge{xvel_purge_in}, yvel_purge{yvel_purge_in}, zvel_purge{zvel_purge_in}
{}

//-------------------------------------------------------------------------------------------------
WatcherReader::WatcherReader(const int nsystem_in, const int max_reports_in,
                             const float force_limit_in, const float speed_limit_in,
                             const bool track_purge_in, const int* nforce_in, const int* nspeed_in,
                             const int* nrattle_in, const int* nshake_in, const float4* forces_in,
                             const int* force_steps_in, const int* force_stages_in,
                             const float4* speeds_in, const int* speed_steps_in,
                             const int* speed_stages_in, const uint2* rattle_fails_in,
                             const uint2* shake_fails_in, const float* rattle_ext_in,
                             const float* shake_ext_in, const float* xvel_purge_in,
                             const float* yvel_purge_in, const float* zvel_purge_in,
                             const float* xang_purge_in, const float* yang_purge_in,
                             const float* zang_purge_in) :
    nsystem{nsystem_in}, max_reports{max_reports_in}, force_limit{force_limit_in},
    speed_limit{speed_limit_in}, track_purge{track_purge_in}, nforce{nforce_in},
    nspeed{nspeed_in}, nrattle{nrattle_in}, nshake{nshake_in}, forces{forces_in},
    force_steps{force_steps_in}, force_stages{force_stages_in}, speeds{speeds_in},
    speed_steps{speed_steps_in}, speed_stages{speed_stages_in}, rattle_fails{rattle_fails_in},
    shake_fails{shake_fails_in}, rattle_ext{rattle_ext_in}, shake_ext{shake_ext_in},
    xvel_purge{xvel_purge_in}, yvel_purge{yvel_purge_in}, zvel_purge{zvel_purge_in}
{}

//-------------------------------------------------------------------------------------------------
WatcherReader::WatcherReader(const WatcherWriter &w) :
    nsystem{w.nsystem}, max_reports{w.max_reports}, force_limit{w.force_limit},
    speed_limit{w.speed_limit}, track_purge{w.track_purge}, nforce{w.nforce},
    nspeed{w.nspeed}, nrattle{w.nrattle}, nshake{w.nshake}, forces{w.forces},
    force_steps{w.force_steps}, force_stages{w.force_stages}, speeds{w.speeds},
    speed_steps{w.speed_steps}, speed_stages{w.speed_stages}, rattle_fails{w.rattle_fails},
    shake_fails{w.shake_fails}, rattle_ext{w.rattle_ext}, shake_ext{w.shake_ext},
    xvel_purge{w.xvel_purge}, yvel_purge{w.yvel_purge}, zvel_purge{w.zvel_purge}
{}

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
    high_speed_count{HybridKind::POINTER, "watcher_lv_cnt"},
    high_speeds{HybridKind::ARRAY, "watcher_lv"},
    high_speed_steps{HybridKind::POINTER, "watcher_lv_steps"},
    high_speed_stages{HybridKind::POINTER, "watcher_lv_stages"},
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
    high_speed_count{original.high_speed_count},
    high_speeds{original.high_speeds},
    high_speed_steps{original.high_speed_steps},
    high_speed_stages{original.high_speed_stages},
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
    high_speed_count{std::move(original.high_speed_count)},
    high_speeds{std::move(original.high_speeds)},
    high_speed_steps{std::move(original.high_speed_steps)},
    high_speed_stages{std::move(original.high_speed_stages)},
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
  high_speed_count = other.high_speed_count;
  high_speeds = other.high_speeds;
  high_speed_steps = other.high_speed_steps;
  high_speed_stages = other.high_speed_stages;
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
  high_speed_count = std::move(other.high_speed_count);
  high_speeds = std::move(other.high_speeds);
  high_speed_steps = std::move(other.high_speed_steps);
  high_speed_stages = std::move(other.high_speed_stages);
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
int Watcher::getLargeForceCount(const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return large_force_count.readHost(0);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return large_force_count.readDevice(0);
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<float4> Watcher::getLargeForces(const HybridTargetLevel tier) const {
  const int count = getLargeForceCount(tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return large_forces.readHost(0, count);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return large_forces.readDevice(0, count);
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int> Watcher::getLargeForceSteps(const HybridTargetLevel tier) const {
  const int count = getLargeForceCount(tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return large_force_steps.readHost(0, count);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return large_force_steps.readDevice(0, count);
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<IntegrationStage> Watcher::getLargeForceStages(const HybridTargetLevel tier) const {
  const int count = getLargeForceCount(tier);
  std::vector<int> ires;
  switch (tier) {
  case HybridTargetLevel::HOST:
    ires = large_force_steps.readHost(0, count);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    ires = large_force_steps.readDevice(0, count);
    break;
#endif
  }
  std::vector<IntegrationStage> result(count);
  for (int i = 0; i < count; i++) {
    result[i] = static_cast<IntegrationStage>(ires[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
float Watcher::getSpeedThreshold() const {
  return speed_threshold;
}

//-------------------------------------------------------------------------------------------------
int Watcher::getHighSpeedCount(const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return high_speed_count.readHost(0);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return high_speed_count.readDevice(0);
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<float4> Watcher::getHighSpeeds(const HybridTargetLevel tier) const {
  const int count = getHighSpeedCount(tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return high_speeds.readHost(0, count);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return high_speeds.readDevice(0, count);
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int> Watcher::getHighSpeedSteps(const HybridTargetLevel tier) const {
  const int count = getHighSpeedCount(tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return high_speed_steps.readHost(0, count);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return high_speed_steps.readDevice(0, count);
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<IntegrationStage> Watcher::getHighSpeedStages(const HybridTargetLevel tier) const {
  const int count = getHighSpeedCount(tier);
  std::vector<int> ires;
  switch (tier) {
  case HybridTargetLevel::HOST:
    ires = high_speed_stages.readHost(0, count);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    ires = high_speed_stages.readDevice(0, count);
    break;
#endif
  }
  std::vector<IntegrationStage> result(count);
  for (int i = 0; i < count; i++) {
    result[i] = static_cast<IntegrationStage>(ires[i]);
  }
  return result;
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
  const int count = getFailedRattleCount(tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return rattle_group_failures.readHost(0, count);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return rattle_group_failures.readDevice(0, count);
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<uint2> Watcher::getShakeFailures(const HybridTargetLevel tier) const {
  const int count = getFailedShakeCount(tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return shake_group_failures.readHost(0, count);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return shake_group_failures.readDevice(0, count);
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<float> Watcher::getRattleViolations(const HybridTargetLevel tier) const {
  const int count = getFailedRattleCount(tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return rattle_violations.readHost(0, count);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return rattle_violations.readDevice(0, count);
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<float> Watcher::getShakeViolations(const HybridTargetLevel tier) const {
  const int count = getFailedShakeCount(tier);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return shake_violations.readHost(0, count);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return shake_violations.readDevice(0, count);
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const WatcherReader Watcher::data(const HybridTargetLevel tier) const {
  return WatcherReader(poly_ps_ptr->getSystemCount(), max_reports, force_threshold,
                       speed_threshold, track_momentum_purge, large_force_count.data(tier),
                       high_speed_count.data(tier), failed_rattle_count.data(tier),
                       failed_shake_count.data(tier), large_forces.data(tier),
                       large_force_steps.data(tier), large_force_stages.data(tier),
                       high_speeds.data(tier), high_speed_steps.data(tier),
                       high_speed_stages.data(tier), rattle_group_failures.data(tier),
                       shake_group_failures.data(tier), rattle_violations.data(tier),
                       shake_violations.data(tier), x_velocity_purge.data(tier),
                       y_velocity_purge.data(tier), z_velocity_purge.data(tier),
                       x_angular_purge.data(tier), y_angular_purge.data(tier),
                       z_angular_purge.data(tier));
}

//-------------------------------------------------------------------------------------------------
WatcherWriter Watcher::data(const HybridTargetLevel tier) {
  return WatcherWriter(poly_ps_ptr->getSystemCount(), max_reports, force_threshold,
                       speed_threshold, track_momentum_purge, large_force_count.data(tier),
                       high_speed_count.data(tier), failed_rattle_count.data(tier),
                       failed_shake_count.data(tier), large_forces.data(tier),
                       large_force_steps.data(tier), large_force_stages.data(tier),
                       high_speeds.data(tier), high_speed_steps.data(tier),
                       high_speed_stages.data(tier), rattle_group_failures.data(tier),
                       shake_group_failures.data(tier), rattle_violations.data(tier),
                       shake_violations.data(tier), x_velocity_purge.data(tier),
                       y_velocity_purge.data(tier), z_velocity_purge.data(tier),
                       x_angular_purge.data(tier), y_angular_purge.data(tier),
                       z_angular_purge.data(tier));
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
  high_speed_count.setPointer(&int_data, 1, 1);
  failed_rattle_count.setPointer(&int_data, 2, 1);
  failed_shake_count.setPointer(&int_data, 3, 1);
  large_force_steps.setPointer(&int_data,                             4, max_reports);
  large_force_stages.setPointer(&int_data,       4 + padded_max_reports, max_reports);
  high_speed_steps.setPointer(&int_data,   4 + (2 * padded_max_reports), max_reports);
  high_speed_stages.setPointer(&int_data,  4 + (3 * padded_max_reports), max_reports);

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
  high_speeds.resize(max_reports);
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
