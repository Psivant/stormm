// -*-c++-*-
#ifndef STORMM_TRAJECTORY_UTIL_H
#define STORMM_TRAJECTORY_UTIL_H

#include <vector>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Reporting/reporting_enumerators.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_cache_map.h"
#include "coordinateframe.h"
#include "phasespace.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

using review::BrokenAsciiCode;
using synthesis::Condensate;
using synthesis::PhaseSpaceSynthesis;
using synthesis::SynthesisCacheMap;

/// \brief Write trajectory frames for all systems within a synthesis, accessing memory from the
///        coordinate synthesis on the CPU host or GPU device and in terms of the positions,
///        velocities, or forces on particles as required.  The current simulation time,
///        expected to be consistent for all systems in the synthesis, is also written out.
///
/// \param poly_ps       The coordinate synthesis
/// \param staging_zone  A staging area in which to render molecules for export to files.  The
///                      condensate can be the destination for copy operations of coordinates held
///                      in CPU or GPU memory, and these copy operations will implicitly perform
///                      the fixed-precision to real-valued conversion needed to use writeFrame().
///                      The memory resources of this object may reside entirely on the CPU host,
///                      but it must have memory available on the CPU host.
/// \param scmap         A map relating systems in the synthesis back to the original collection
///                      of input files
/// \param crd_format    Format of teh resulting trajectory files
/// \param traj_kind     The type of trajectory to write: positions, velocities, or forces
/// \param current_time  Current time, a synthesis-wide parameter
/// \param tier          Indicate whether to access data on the CPU host or GPU device
/// \param gpu           Details of the GPU wherein the memory of interest resides, and the device
///                      which will carry out the extraction to CPU temporary memory if coordinates
///                      are to be pulled from the device side
/// \param recovery      Indicate what to do if a fixed-column file format cannot accommodate a
///                      very large number
void writeSynthesisFrames(const PhaseSpaceSynthesis &poly_ps, Condensate *staging_zone,
                          const Hybrid<int2> &system_pairs, const SynthesisCacheMap &scmap,
                          const CoordinateFileKind crd_format = CoordinateFileKind::AMBER_CRD,
                          double current_time = 0.0,
                          const TrajectoryKind traj_kind = TrajectoryKind::POSITIONS,
                          HybridTargetLevel tier = HybridTargetLevel::HOST,
                          const GpuDetails &gpu = null_gpu,
                          const BrokenAsciiCode recovery = BrokenAsciiCode::NONE);

} // namespace trajectory
} // namesapce stormm

#endif
