#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "debug.h"

namespace stormm {
namespace debug {

using numerics::force_scale_nonoverflow_bits;
using numerics::velocity_scale_nonoverflow_bits;

//-------------------------------------------------------------------------------------------------
void checkAtomicForces(WatcherWriter *bugw, const PsSynthesisReader &poly_psr, const int step,
                       const IntegrationStage when) {
  for (int i = 0; i < poly_psr.system_count; i++) {
    const int j_llim = poly_psr.atom_starts[i];
    const int j_hlim = j_llim + poly_psr.atom_counts[i];
    for (int j = j_llim; j < j_hlim; j++) {
      float fx, fy, fz;
      if (poly_psr.frc_bits > force_scale_nonoverflow_bits) {
        fx = hostInt95ToDouble(poly_psr.xfrc[j], poly_psr.xfrc_ovrf[j]) * poly_psr.inv_frc_scale;
        fy = hostInt95ToDouble(poly_psr.yfrc[j], poly_psr.yfrc_ovrf[j]) * poly_psr.inv_frc_scale;
        fz = hostInt95ToDouble(poly_psr.zfrc[j], poly_psr.zfrc_ovrf[j]) * poly_psr.inv_frc_scale;
      }
      else {
        fx = (float)(poly_psr.xfrc[j]) * poly_psr.inv_frc_scale;
        fy = (float)(poly_psr.yfrc[j]) * poly_psr.inv_frc_scale;
        fz = (float)(poly_psr.zfrc[j]) * poly_psr.inv_frc_scale;
      }
      const float fmag = sqrtf((fx * fx) + (fy * fy) + (fz * fz));
      if (fmag >= bugw->force_limit) {
        const int event_idx = bugw->nforce[0];
        if (event_idx < bugw->max_reports) {
          Ecumenical4 ec_idx = { .i = j };
          bugw->forces[event_idx] = { fx, fy, fz, ec_idx.f };
          bugw->force_steps[event_idx] = step;
          bugw->force_stages[event_idx] = static_cast<int>(when);
          bugw->force_contexts[event_idx] = static_cast<int>(AnomalyContext::PHASESPACE_SYNTHESIS);
          bugw->nforce[0] = event_idx + 1;
        }
        else {
          j = j_hlim;
          i = poly_psr.system_count;
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void checkAtomicForces(Watcher *anom, const PhaseSpaceSynthesis *poly_ps, const int step,
		       const IntegrationStage when) {
  WatcherWriter bugw = anom->data();
  const PsSynthesisReader poly_psr = poly_ps->data();
  checkAtomicForces(&bugw, poly_psr, step, when);
}

//-------------------------------------------------------------------------------------------------
void checkAtomicForces(Watcher *anom, const PhaseSpaceSynthesis *poly_ps,
		       const CoordinateCycle orientation, const int step,
                       const IntegrationStage when) {
  WatcherWriter bugw = anom->data();
  const PsSynthesisReader poly_psr = poly_ps->data(orientation);
  checkAtomicForces(&bugw, poly_psr, step, when);
}

//-------------------------------------------------------------------------------------------------
void checkAtomicSpeeds(WatcherWriter *bugw, const PsSynthesisReader &poly_psr, const int step,
                       const IntegrationStage when) {
  for (int i = 0; i < poly_psr.system_count; i++) {
    const int j_llim = poly_psr.atom_starts[i];
    const int j_hlim = j_llim + poly_psr.atom_counts[i];
    for (int j = j_llim; j < j_hlim; j++) {
      float vx, vy, vz;
      if (poly_psr.vel_bits > velocity_scale_nonoverflow_bits) {
        vx = hostInt95ToDouble(poly_psr.vxalt[j], poly_psr.vxalt_ovrf[j]) * poly_psr.inv_vel_scale;
        vy = hostInt95ToDouble(poly_psr.vyalt[j], poly_psr.vyalt_ovrf[j]) * poly_psr.inv_vel_scale;
        vz = hostInt95ToDouble(poly_psr.vzalt[j], poly_psr.vzalt_ovrf[j]) * poly_psr.inv_vel_scale;
      }
      else {
        vx = (float)(poly_psr.vxalt[j]) * poly_psr.inv_vel_scale;
        vy = (float)(poly_psr.vyalt[j]) * poly_psr.inv_vel_scale;
        vz = (float)(poly_psr.vzalt[j]) * poly_psr.inv_vel_scale;
      }
      const float vmag = sqrtf((vx * vx) + (vy * vy) + (vz * vz));
      if (vmag >= bugw->speed_limit) {
        const int event_idx = bugw->nspeed[0];
        if (event_idx < bugw->max_reports) {
          Ecumenical4 ec_idx = { .i = j };
          bugw->speeds[event_idx] = { vx, vy, vz, ec_idx.f };
          bugw->speed_steps[event_idx] = step;
          bugw->speed_stages[event_idx] = static_cast<int>(when);
          bugw->nspeed[0] = event_idx + 1;
        }
        else {
          j = j_hlim;
          i = poly_psr.system_count;
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void checkAtomicSpeeds(Watcher *anom, const PhaseSpaceSynthesis *poly_ps, const int step,
		       const IntegrationStage when) {
  WatcherWriter bugw = anom->data();
  const PsSynthesisReader poly_psr = poly_ps->data();
  checkAtomicSpeeds(&bugw, poly_psr, step, when);
}

//-------------------------------------------------------------------------------------------------
void checkAtomicSpeeds(Watcher *anom, const PhaseSpaceSynthesis *poly_ps,
		       const CoordinateCycle orientation, const int step,
                       const IntegrationStage when) {
  WatcherWriter bugw = anom->data();
  const PsSynthesisReader poly_psr = poly_ps->data(orientation);
  checkAtomicSpeeds(&bugw, poly_psr, step, when);
}

} // namespace debug
} // namespace stormm
