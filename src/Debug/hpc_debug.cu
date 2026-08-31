// -*-c++-*-
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Constants/fixed_precision.h"
#include "Constants/hpc_bounds.h"
#include "MolecularMechanics/dynamics_intervention.h"
#include "Debug/watcher.h"
#include "Synthesis/phasespace_synthesis.h"
#include "debugging_enumerators.h"
#include "hpc_debug.h"

namespace stormm {
namespace debug {

using card::HybridTargetLevel;
using numerics::force_scale_nonoverflow_bits;
using numerics::velocity_scale_nonoverflow_bits;
using mm::dyna_tk;

#include "Math/rounding.cui"
#include "Numerics/accumulation.cui"
  
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kCheckAtomicForces(WatcherWriter bugw, const PsSynthesisReader poly_psr, const int step,
                   const IntegrationStage when) {
  __shared__ int current_report_count;
  
  // As a guard against overflow in the total number of reports, each thread block will begin by
  // checking the number of reports submitted thus far.  If the number of reports already
  // exceeds the maximum allowed number, no further work will be done.
  if (threadIdx.x == 0) {
    current_report_count = bugw.nforce[0];
  }
  __syncthreads();
  const int last_sys_idx = poly_psr.system_count - 1;
  const int max_atoms = poly_psr.atom_starts[last_sys_idx] + poly_psr.atom_counts[last_sys_idx];
  const int max_padded_atoms = devcRoundUp(max_atoms, warp_size_int);
  const int kstride = gridDim.x * blockDim.x;
  const int ilim = (current_report_count >= bugw.max_reports) ? 0 : max_padded_atoms;
  for (int i = (blockDim.x * blockIdx.x) + threadIdx.x; i < max_padded_atoms; i += kstride) {
    float fx, fy, fz;
    if (poly_psr.frc_bits > force_scale_nonoverflow_bits) {
      fx = int95ToDouble(poly_psr.xfrc[i], poly_psr.xfrc_ovrf[i]) * poly_psr.inv_frc_scale;
      fy = int95ToDouble(poly_psr.yfrc[i], poly_psr.yfrc_ovrf[i]) * poly_psr.inv_frc_scale;
      fz = int95ToDouble(poly_psr.zfrc[i], poly_psr.zfrc_ovrf[i]) * poly_psr.inv_frc_scale;
    }
    else {
      fx = (float)(poly_psr.xfrc[i]) * poly_psr.inv_frc_scale;
      fy = (float)(poly_psr.yfrc[i]) * poly_psr.inv_frc_scale;
      fz = (float)(poly_psr.zfrc[i]) * poly_psr.inv_frc_scale;
    }
    const float fmag = sqrtf((fx * fx) + (fy * fy) + (fz * fz));
    if (fmag >= bugw.force_limit) {

      // No work units are needed because the likelihood of a threshold-exceeding force is
      // expected to be very small.  Only in such cases will a lengthier binary search be done
      // to identify the system in which the force occurred and its exact atom index.
      int llim = 0;
      int hlim = last_sys_idx + 1;
      int mid = llim + ((hlim - llim) / 2);
      bool seek_system;
      do {
        if (mid == last_sys_idx) {
          seek_system = false;
        }
        else {
          if (i < poly_psr.atom_starts[mid]) {
            seek_system = true;
            hlim = mid;
            mid = llim + ((hlim - llim) / 2);
          }
          else if (i >= poly_psr.atom_starts[mid] && i < poly_psr.atom_starts[mid + 1]) {
            seek_system = false;
          }
          else {
            seek_system = true;
            llim = mid + 1;
            mid = llim + ((hlim - llim) / 2);
          }
        }
      } while (seek_system);

      // Guard against the case that the "large force" is found in the padded part of the array.
      // This is not a problem per se, but it may indicate that something is dumping data where it
      // should not.  Such behavior will be checked by a separate kernel.
      if (i < poly_psr.atom_starts[mid] + poly_psr.atom_counts[mid]) {
        const int event_idx = atomicAdd(bugw.nforce, 1);
        if (event_idx < bugw.max_reports) {
          bugw.forces[event_idx] = { fx, fy, fz, __int_as_float(i) };
          bugw.force_steps[event_idx] = step;
          bugw.force_stages[event_idx] = (int)(when);
          bugw.force_contexts[event_idx] = (int)(AnomalyContext::PHASESPACE_SYNTHESIS);
        }
        else {

          // Set the loop control variable such that no further iterations of the loop over all
          // atoms will occur.  In this way, the total number of reports will not grow to be
          // greater than the maximum allowed number of reports plus the number of threads in the
          // kernel, and only as many reports as have been allocated will be transcribed.
          i = ilim;
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void checkAtomicForces(WatcherWriter *bugw, const PsSynthesisReader &poly_psr, const int step,
                       const IntegrationStage when, const GpuDetails &gpu) {
  kCheckAtomicForces<<<gpu.getSMPCount(), large_block_size>>>(*bugw, poly_psr, step, when);
}

//-------------------------------------------------------------------------------------------------
void checkAtomicForces(Watcher *anom, const PhaseSpaceSynthesis *poly_ps, const int step,
                       const IntegrationStage when, const GpuDetails &gpu) {
  WatcherWriter bugw = anom->data(HybridTargetLevel::DEVICE);
  const PsSynthesisReader poly_psr = poly_ps->data(HybridTargetLevel::DEVICE);
  checkAtomicForces(&bugw, poly_psr, step, when, gpu);
}

//-------------------------------------------------------------------------------------------------
void checkAtomicForces(Watcher *anom, const PhaseSpaceSynthesis *poly_ps,
                       const CoordinateCycle orientation, const int step,
                       const IntegrationStage when, const GpuDetails &gpu) {
  WatcherWriter bugw = anom->data(HybridTargetLevel::DEVICE);
  const PsSynthesisReader poly_psr = poly_ps->data(orientation, HybridTargetLevel::DEVICE);
  checkAtomicForces(&bugw, poly_psr, step, when, gpu);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kCheckAtomicSpeeds(WatcherWriter bugw, const PsSynthesisReader poly_psr, const int step,
                   const IntegrationStage when) {
  __shared__ int current_report_count;
  
  // As a guard against overflow in the total number of reports, each thread block will begin by
  // checking the number of reports submitted thus far.  If the number of reports already
  // exceeds the maximum allowed number, no further work will be done.
  if (threadIdx.x == 0) {
    current_report_count = bugw.nspeed[0];
  }
  __syncthreads();
  
  const int last_sys_idx = poly_psr.system_count - 1;
  const int max_atoms = poly_psr.atom_starts[last_sys_idx] + poly_psr.atom_counts[last_sys_idx];
  const int max_padded_atoms = devcRoundUp(max_atoms, warp_size_int);
  const int kstride = gridDim.x * blockDim.x;
  const int ilim = (current_report_count >= bugw.max_reports) ? 0 : max_padded_atoms;
  for (int i = (blockDim.x * blockIdx.x) + threadIdx.x; i < max_padded_atoms; i += kstride) {
    float vx, vy, vz;
    if (poly_psr.vel_bits > velocity_scale_nonoverflow_bits) {
      vx = int95ToDouble(poly_psr.vxalt[i], poly_psr.vxalt_ovrf[i]) * poly_psr.inv_vel_scale;
      vy = int95ToDouble(poly_psr.vyalt[i], poly_psr.vyalt_ovrf[i]) * poly_psr.inv_vel_scale;
      vz = int95ToDouble(poly_psr.vzalt[i], poly_psr.vzalt_ovrf[i]) * poly_psr.inv_vel_scale;
    }
    else {
      vx = (float)(poly_psr.vxalt[i]) * poly_psr.inv_vel_scale;
      vy = (float)(poly_psr.vyalt[i]) * poly_psr.inv_vel_scale;
      vz = (float)(poly_psr.vzalt[i]) * poly_psr.inv_vel_scale;
    }
    const float vmag = sqrtf((vx * vx) + (vy * vy) + (vz * vz));
    if (vmag >= bugw.speed_limit) {

      // No work units are needed because the likelihood of a threshold-exceeding force is
      // expected to be very small.  Only in such cases will a lengthier binary search be done
      // to identify the system in which the force occurred and its exact atom index.
      int llim = 0;
      int hlim = last_sys_idx + 1;
      int mid = llim + ((hlim - llim) / 2);
      bool seek_system;
      do {
        if (mid == last_sys_idx) {
          seek_system = false;
        }
        else {
          if (i < poly_psr.atom_starts[mid]) {
            seek_system = true;
            hlim = mid;
            mid = llim + ((hlim - llim) / 2);
          }
          else if (i >= poly_psr.atom_starts[mid] && i < poly_psr.atom_starts[mid + 1]) {
            seek_system = false;
          }
          else {
            seek_system = true;
            llim = mid + 1;
            mid = llim + ((hlim - llim) / 2);
          }
        }
      } while (seek_system);

      // Guard against the case that the "large force" is found in the padded part of the array.
      // This is not a problem per se, but it may indicate that something is dumping data where it
      // should not.  Such behavior will be checked by a separate kernel.
      if (i < poly_psr.atom_starts[mid] + poly_psr.atom_counts[mid]) {
        const int event_idx = atomicAdd(bugw.nspeed, 1);
        if (event_idx < bugw.max_reports) {
          bugw.speeds[event_idx] = { vx, vy, vz, __int_as_float(i - poly_psr.atom_starts[mid]) };
          bugw.speed_steps[event_idx] = step;
          bugw.speed_stages[event_idx] = (int)(when);
        }
        else {

          // Set the loop control variable such that no further iterations of the loop over all
          // atoms will occur.  In this way, the total number of reports will not grow to be
          // greater than the maximum allowed number of reports plus the number of threads in the
          // kernel, and only as many reports as have been allocated will be transcribed.
          i = ilim;
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void checkAtomicSpeeds(WatcherWriter *bugw, const PsSynthesisReader &poly_psr, const int step,
                       const IntegrationStage when, const GpuDetails &gpu) {
  kCheckAtomicSpeeds<<<gpu.getSMPCount(), large_block_size>>>(*bugw, poly_psr, step, when);
}

//-------------------------------------------------------------------------------------------------
void checkAtomicSpeeds(Watcher *anom, const PhaseSpaceSynthesis *poly_ps, const int step,
                       const IntegrationStage when, const GpuDetails &gpu) {
  WatcherWriter bugw = anom->data(HybridTargetLevel::DEVICE);
  const PsSynthesisReader poly_psr = poly_ps->data(HybridTargetLevel::DEVICE);
  checkAtomicSpeeds(&bugw, poly_psr, step, when, gpu);
}

//-------------------------------------------------------------------------------------------------
void checkAtomicSpeeds(Watcher *anom, const PhaseSpaceSynthesis *poly_ps,
                       const CoordinateCycle orientation, const int step,
                       const IntegrationStage when, const GpuDetails &gpu) {
  WatcherWriter bugw = anom->data(HybridTargetLevel::DEVICE);
  const PsSynthesisReader poly_psr = poly_ps->data(orientation, HybridTargetLevel::DEVICE);
  checkAtomicSpeeds(&bugw, poly_psr, step, when, gpu);
}

//-------------------------------------------------------------------------------------------------
void checkPeriodicDynamics(int step_number, IntegrationStage when) {

  // Make a local copy of the GPU specs, for kernel launch control.  This and other intervention
  // functions will not be able to store launch grid information in the various kernel managers,
  // but it is feasible to design kernels of reasonable efficiency with launch grids built around
  // a fixed size and the number of GPU streaming multiprocessors.
  const GpuDetails gpu = dyna_tk.getGpuDetails();
  
  // Unpack the anomaly reporting abstract for the GPU.  The Watcher class object it points to
  // will have been established for the synthesis of the calculations.
  WatcherWriter bugw = dyna_tk.getAnomalyReportingData(HybridTargetLevel::DEVICE);
}

} // namespace debug
} // namespase stormm
