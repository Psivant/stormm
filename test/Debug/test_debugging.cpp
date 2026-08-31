#include <algorithm>
#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/Debug/debug.h"
#include "../../src/Debug/hpc_debug.h"
#include "../../src/FileManagement/file_util.h"
#include "../../src/Math/series_ops.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/section_contents.h"
#include "../../src/Debug/watcher.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/synthesis_abstracts.h"
#include "../../src/Synthesis/nonbonded_workunit.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/systemcache.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test_enumerators.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::card;
using namespace stormm::data_types;
using namespace stormm::debug;
using namespace stormm::diskutil;
using namespace stormm::numerics;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::synthesis;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Add forces to a synthesis of systems, then make some of the forces arbitrarily large to test
// whether the debugging function, and associated GPU kernel if running in HPC mode, will pick
// them up.
// 
// Arguments:
//   poly_ps:   The synthesis of coordinates, into which bogus forces will be written
//   anom:      Tracking object for recording anomalous events, such as large forces
//   xrs:       Soure of random numbers for noise
//   poly_ag:   The synthesis of topologies
//   do_tests:  
//   gpu:       Specifications of the GPU that will carry out calculations, if GPU mode is enabled
//-------------------------------------------------------------------------------------------------
void challengeWithLargeForces(PhaseSpaceSynthesis *poly_ps, Watcher *anom,
                              Xoshiro256ppGenerator *xrs, const AtomGraphSynthesis &poly_ag,
                              const TestPriority do_tests, const GpuDetails &gpu) {

  // Seed the system with random forces, on occasion with a large force
  PsSynthesisWriter host_poly_psw = poly_ps->data();
  const double inv_scl = host_poly_psw.inv_frc_scale;
  const double limiter = anom->getForceThreshold() * host_poly_psw.frc_scale;
  const int95_t posi_limit = hostDoubleToInt95(limiter * 0.5);
  const int95_t nega_limit = hostDoubleToInt95(-limiter * 0.5);
  for (int i = 0; i < host_poly_psw.system_count; i++) {
    const size_t sys_start = host_poly_psw.atom_starts[i];
    const size_t natom = host_poly_psw.atom_counts[i];
    addRandomNoise(xrs, &host_poly_psw.xfrc[sys_start], &host_poly_psw.xfrc_ovrf[sys_start],
                   &host_poly_psw.yfrc[sys_start], &host_poly_psw.yfrc_ovrf[sys_start],
                   &host_poly_psw.zfrc[sys_start], &host_poly_psw.zfrc_ovrf[sys_start], natom,
                   0.25, host_poly_psw.frc_scale);

    // Clamp the noise to ensure that nothing is above the threshold (this would be very unlikely
    // for the Gaussian spread applied).
    for (size_t j = sys_start; j < sys_start + natom; j++) {
      const double xc = hostInt95ToDouble(host_poly_psw.xfrc[j], host_poly_psw.xfrc_ovrf[j]);
      const double yc = hostInt95ToDouble(host_poly_psw.yfrc[j], host_poly_psw.yfrc_ovrf[j]);
      const double zc = hostInt95ToDouble(host_poly_psw.zfrc[j], host_poly_psw.zfrc_ovrf[j]);
      const double rc = sqrt((xc * xc) + (yc * yc) + (zc * zc)) * inv_scl;
      if (fabs(rc) >= limiter) {
        if (xc < 0.0) {
          host_poly_psw.xfrc[j] = nega_limit.x;
          host_poly_psw.xfrc_ovrf[j] = nega_limit.y;
        }
        else {
          host_poly_psw.xfrc[j] = posi_limit.x;
          host_poly_psw.xfrc_ovrf[j] = posi_limit.y;
        }
        if (yc < 0.0) {
          host_poly_psw.yfrc[j] = nega_limit.x;
          host_poly_psw.yfrc_ovrf[j] = nega_limit.y;
        }
        else {
          host_poly_psw.yfrc[j] = posi_limit.x;
          host_poly_psw.yfrc_ovrf[j] = posi_limit.y;
        }
        if (zc < 0.0) {
          host_poly_psw.zfrc[j] = nega_limit.x;
          host_poly_psw.zfrc_ovrf[j] = nega_limit.y;
        }
        else {
          host_poly_psw.zfrc[j] = posi_limit.x;
          host_poly_psw.zfrc_ovrf[j] = posi_limit.y;
        }
      }
    }

    // Create three arbitrarily large forces in each system.
    const double dnatom = natom;
    std::vector<int> used_indices(3, -1);
    for (int j = 0; j < 3; j++) {
      const int95_t breaker_x = hostDoubleToInt95(limiter * (1.1 + xrs->uniformRandomNumber()));
      const int95_t breaker_y = hostDoubleToInt95(limiter * (1.1 + xrs->uniformRandomNumber()));
      const int95_t breaker_z = hostDoubleToInt95(limiter * (1.1 + xrs->uniformRandomNumber()));
      int jpos;
      bool seek_new_index;
      do {
        jpos = static_cast<int>(dnatom * xrs->uniformRandomNumber());
        jpos = std::max(0, jpos);
        jpos = std::min(jpos, static_cast<int>(natom) - 1);
        seek_new_index = (locateValue(used_indices, jpos) < 3);
      } while (seek_new_index);
      used_indices[j] = jpos;
      const size_t jpos_zu = sys_start + static_cast<size_t>(jpos);
      host_poly_psw.xfrc[jpos_zu]      = breaker_x.x;
      host_poly_psw.xfrc_ovrf[jpos_zu] = breaker_x.y;
      host_poly_psw.yfrc[jpos_zu]      = breaker_y.x;
      host_poly_psw.yfrc_ovrf[jpos_zu] = breaker_y.y;
      host_poly_psw.zfrc[jpos_zu]      = breaker_z.x;
      host_poly_psw.zfrc_ovrf[jpos_zu] = breaker_z.y;
    }
  }
#ifdef STORMM_USE_HPC
  HybridTargetLevel tier = HybridTargetLevel::DEVICE;
#else
  HybridTargetLevel tier = HybridTargetLevel::HOST;
#endif
  WatcherWriter bugw = anom->data(tier);
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
#ifdef STORMM_USE_HPC
  poly_ps->upload();
  checkAtomicForces(&bugw, poly_psw, 15, IntegrationStage::GEOMETRY_CONSTRAINT, gpu);
  anom->downloadEventCounts();
#else
  checkAtomicForces(&bugw, poly_psw, 15, IntegrationStage::GEOMETRY_CONSTRAINT);
#endif
  
  // If calculations were run on the GPU, the event counts have been downloaded and should be
  // available on the host.  If calculations were run on the host, the event count will be ready.
  check(anom->getLargeForceCount(), RelationalOperator::EQUAL, 3 * host_poly_psw.system_count,
        "A total of 3 large forces were seeded in each of " +
        std::to_string(host_poly_psw.system_count) + " systems, but " +
        std::to_string(anom->getLargeForceCount()) + " events were recorded in all.", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Add small, random velocities to a synthesis of systems, then make some of the velocities
// arbitrarily large to test whether the debugging function, and associated GPU kernel if running
// in HPC mode, will pick them up.
//-------------------------------------------------------------------------------------------------
void challengeWithHighSpeeds(PhaseSpaceSynthesis *poly_ps, Watcher *anom,
                             Xoshiro256ppGenerator *xrs, const AtomGraphSynthesis &poly_ag,
                             const TestPriority do_tests, const GpuDetails &gpu) {

  // Seed the system with random forces, on occasion with a large force
  PsSynthesisWriter host_poly_psw = poly_ps->data();
  const double inv_scl = host_poly_psw.inv_vel_scale;
  const double limiter = anom->getSpeedThreshold() * host_poly_psw.vel_scale;
  const int95_t posi_limit = hostDoubleToInt95(limiter * 0.5);
  const int95_t nega_limit = hostDoubleToInt95(-limiter * 0.5);
  for (int i = 0; i < host_poly_psw.system_count; i++) {
    const size_t sys_start = host_poly_psw.atom_starts[i];
    const size_t natom = host_poly_psw.atom_counts[i];
    addRandomNoise(xrs, &host_poly_psw.vxalt[sys_start], &host_poly_psw.vxalt_ovrf[sys_start],
                   &host_poly_psw.vyalt[sys_start], &host_poly_psw.vyalt_ovrf[sys_start],
                   &host_poly_psw.vzalt[sys_start], &host_poly_psw.vzalt_ovrf[sys_start], natom,
                   0.01, host_poly_psw.vel_scale);

    // Clamp the noise to ensure that nothing is above the threshold (this would be very unlikely
    // for the Gaussian spread applied).
    for (size_t j = sys_start; j < sys_start + natom; j++) {
      const double xv = hostInt95ToDouble(host_poly_psw.vxalt[j], host_poly_psw.vxalt_ovrf[j]);
      const double yv = hostInt95ToDouble(host_poly_psw.vyalt[j], host_poly_psw.vyalt_ovrf[j]);
      const double zv = hostInt95ToDouble(host_poly_psw.vzalt[j], host_poly_psw.vzalt_ovrf[j]);
      const double rv = sqrt((xv * xv) + (yv * yv) + (zv * zv)) * inv_scl;
      if (fabs(rv) >= limiter) {
        if (xv < 0.0) {
          host_poly_psw.vxalt[j] = nega_limit.x;
          host_poly_psw.vxalt_ovrf[j] = nega_limit.y;
        }
        else {
          host_poly_psw.vxalt[j] = posi_limit.x;
          host_poly_psw.vxalt_ovrf[j] = posi_limit.y;
        }
        if (yv < 0.0) {
          host_poly_psw.vyalt[j] = nega_limit.x;
          host_poly_psw.vyalt_ovrf[j] = nega_limit.y;
        }
        else {
          host_poly_psw.vyalt[j] = posi_limit.x;
          host_poly_psw.vyalt_ovrf[j] = posi_limit.y;
        }
        if (zv < 0.0) {
          host_poly_psw.vzalt[j] = nega_limit.x;
          host_poly_psw.vzalt_ovrf[j] = nega_limit.y;
        }
        else {
          host_poly_psw.vzalt[j] = posi_limit.x;
          host_poly_psw.vzalt_ovrf[j] = posi_limit.y;
        }
      }
    }

    // Create three arbitrarily large particle velocities in each system.
    const double dnatom = natom;
    std::vector<int> used_indices(3, -1);
    for (int j = 0; j < 3; j++) {
      const int95_t breaker_x = hostDoubleToInt95(limiter * (1.1 + xrs->uniformRandomNumber()));
      const int95_t breaker_y = hostDoubleToInt95(limiter * (1.1 + xrs->uniformRandomNumber()));
      const int95_t breaker_z = hostDoubleToInt95(limiter * (1.1 + xrs->uniformRandomNumber()));
      int jpos;
      bool seek_new_index;
      do {
        jpos = static_cast<int>(dnatom * xrs->uniformRandomNumber());
        jpos = std::max(0, jpos);
        jpos = std::min(jpos, static_cast<int>(natom) - 1);
        seek_new_index = (locateValue(used_indices, jpos) < 3);
      } while (seek_new_index);
      used_indices[j] = jpos;
      const size_t jpos_zu = sys_start + static_cast<size_t>(jpos);
      host_poly_psw.vxalt[jpos_zu]      = breaker_x.x;
      host_poly_psw.vxalt_ovrf[jpos_zu] = breaker_x.y;
      host_poly_psw.vyalt[jpos_zu]      = breaker_y.x;
      host_poly_psw.vyalt_ovrf[jpos_zu] = breaker_y.y;
      host_poly_psw.vzalt[jpos_zu]      = breaker_z.x;
      host_poly_psw.vzalt_ovrf[jpos_zu] = breaker_z.y;
    }
  }
#ifdef STORMM_USE_HPC
  HybridTargetLevel tier = HybridTargetLevel::DEVICE;
#else
  HybridTargetLevel tier = HybridTargetLevel::HOST;
#endif
  WatcherWriter bugw = anom->data(tier);
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
#ifdef STORMM_USE_HPC
  poly_ps->upload();
  checkAtomicSpeeds(&bugw, poly_psw, 15, IntegrationStage::GEOMETRY_CONSTRAINT, gpu);
  anom->downloadEventCounts();
#else
  checkAtomicSpeeds(&bugw, poly_psw, 15, IntegrationStage::GEOMETRY_CONSTRAINT);
#endif

  // If calculations were run on the GPU, the event counts have been downloaded and should be
  // available on the host.  If calculations were run on the host, the event count will be ready.
  check(anom->getHighSpeedCount(), RelationalOperator::EQUAL, 3 * host_poly_psw.system_count,
        "A total of 3 high speeds were seeded in each of " +
        std::to_string(host_poly_psw.system_count) + " systems, but " +
        std::to_string(anom->getHighSpeedCount()) + " events were recorded in all.", do_tests);
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  StopWatch timer("Debugging tests");
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Get the GPU specs
#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
#else
  const GpuDetails gpu = null_gpu;
#endif
  
  // Section 1
  section("Test PhaseSpaceSynthesis for events");

  // Collect coordinates and topologies
  const char osc = osSeparator();
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";

  // Compile systems
  const std::vector<std::string> pbc_mols_ph = { "bromobenzene", "tip3p", "tip4p",
                                                 "trpcage_in_water", "trpcage_in_water",
                                                 "ubiquitin", "ubiquitin", "drug_example" };
  TestSystemManager pbc_tsm_ph(base_top_name, "top", pbc_mols_ph, base_crd_name, "inpcrd",
                               pbc_mols_ph);
  const int nsys = pbc_mols_ph.size();
  const std::vector<int> all_systems = incrementingSeries<int>(0, nsys);
  AtomGraphSynthesis poly_ag = pbc_tsm_ph.exportAtomGraphSynthesis(all_systems);
  PhaseSpaceSynthesis poly_ps = pbc_tsm_ph.exportPhaseSpaceSynthesis(all_systems);

  // Detect large forces, speeds, and positional changes
  section(1);
  Xoshiro256ppGenerator xrs(oe.getRandomSeed());
  Watcher anom(poly_ps, poly_ag);
  challengeWithLargeForces(&poly_ps, &anom, &xrs, poly_ag, pbc_tsm_ph.getTestingStatus(), gpu);
  challengeWithHighSpeeds(&poly_ps, &anom, &xrs, poly_ag, pbc_tsm_ph.getTestingStatus(), gpu);
  
  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
