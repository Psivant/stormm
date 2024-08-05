#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Analysis/comparison_guide.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/rmsd.h"
#include "../../src/Structure/rmsd_plan.h"
#include "../../src/Structure/structure_enumerators.h"
#include "../../src/Synthesis/condensate.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/synthesis_cache_map.h"
#include "../../src/Synthesis/synthesis_enumerators.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinate_copy.h"
#include "../../src/Trajectory/trajectory_enumerators.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::analysis;
using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::stmath;
using namespace stormm::numerics;
using namespace stormm::parse;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::structure;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Quantify the difference in Cartesian X, Y, and Z coordinates.  Various overloads serve different
// coordinate formats.
//
// Arguments:
//   xcrd:      Experimental Cartesian X coordinates
//   ycrd:      Experimental Cartesian Y coordinates
//   zcrd:      Experimental Cartesian Z coordinates
//   xref:      Reference Cartesian X coordinates
//   yref:      Reference Cartesian Y coordinates
//   zref:      Reference Cartesian Z coordinates
//   natom:     The number of atoms to analyze
//   copy_msg:  Message to display, describing the type of copying that led to the experimental
//              coordinates, in the event that any inconsistencies are found
//   do_test:   Indication of whether the test is feasible, based on prior success in reading
//              critical coordiante files
//   overall:   The expected overall difference that should be found
//   bounds:    The expected overall difference that should be found
//-------------------------------------------------------------------------------------------------
template <typename Tcrd, typename Tref>
void diffCoordinates(const Tcrd* xcrd, const Tcrd* ycrd, const Tcrd* zcrd, const int* xcrd_ovrf,
                     const int* ycrd_ovrf, const int* zcrd_ovrf, const double crd_scale,
                     const Tref* xref, const Tref* yref, const Tref* zref, const int* xref_ovrf,
                     const int* yref_ovrf, const int* zref_ovrf, const double ref_scale,
                     const int natom, const double* umat, const double* umat_ref,
                     const std::string &copy_msg, const TestPriority do_test,
                     const double overall = tiny, const double bounds = tiny) {
  bool xoff = false;
  bool yoff = false;
  bool zoff = false;
  double rmsd = 0.0;
  std::vector<double> x_errors, y_errors, z_errors;
  for (int i = 0; i < natom; i++) {
    const double tx_crd = (xcrd_ovrf == nullptr) ?
                          static_cast<double>(xcrd[i]) / crd_scale :
                          hostInt95ToDouble(xcrd[i], xcrd_ovrf[i]) / crd_scale;
    const double ty_crd = (ycrd_ovrf == nullptr) ?
                          static_cast<double>(ycrd[i]) / crd_scale :
                          hostInt95ToDouble(ycrd[i], ycrd_ovrf[i]) / crd_scale;
    const double tz_crd = (zcrd_ovrf == nullptr) ?
                          static_cast<double>(zcrd[i]) / crd_scale :
                          hostInt95ToDouble(zcrd[i], zcrd_ovrf[i]) / crd_scale;
    const double tx_ref = (xref_ovrf == nullptr) ?
                          static_cast<double>(xref[i]) / ref_scale :
                          hostInt95ToDouble(xref[i], xref_ovrf[i]) / ref_scale;
    const double ty_ref = (yref_ovrf == nullptr) ?
                          static_cast<double>(yref[i]) / ref_scale :
                          hostInt95ToDouble(yref[i], yref_ovrf[i]) / ref_scale;
    const double tz_ref = (zref_ovrf == nullptr) ?
                          static_cast<double>(zref[i]) / ref_scale :
                          hostInt95ToDouble(zref[i], zref_ovrf[i]) / ref_scale;
    xoff = (xoff || fabs(tx_crd - tx_ref) > bounds);
    yoff = (yoff || fabs(ty_crd - ty_ref) > bounds);
    zoff = (zoff || fabs(tz_crd - tz_ref) > bounds);
    if (fabs(tx_crd - tx_ref) > bounds) {
      xoff = true;
      x_errors.push_back(fabs(tx_crd - tx_ref));
    }
    if (fabs(ty_crd - ty_ref) > bounds) {
      yoff = true;
      y_errors.push_back(fabs(ty_crd - ty_ref));
    }
    if (fabs(tz_crd - tz_ref) > bounds) {
      zoff = true;
      z_errors.push_back(fabs(tz_crd - tz_ref));
    }
    const double dx = tx_crd - tx_ref;
    const double dy = ty_crd - ty_ref;
    const double dz = tz_crd - tz_ref;
    rmsd += (dx * dx) + (dy * dy) + (dz * dz);
  }
  if (umat != nullptr && umat_ref != nullptr) {
    std::vector<double> umat_stl(9), umat_ref_stl(9);
    for (int i = 0; i < 9; i++) {
      umat_stl[i] = umat[i];
      umat_ref_stl[i] = umat_ref[i];
    }
    check(umat_stl, RelationalOperator::EQUAL, umat_ref_stl, "Transformation matrices for two "
          "coordinate sets disagree after copying " + copy_msg + ".", do_test);
  }
  rmsd = sqrt(rmsd / static_cast<double>(natom));
  std::string fail_msg;
  if (xoff) {
    fail_msg = "X (" + realToString(mean(x_errors), 11, 4, NumberFormat::SCIENTIFIC) + " +/- " +
               realToString(variance(x_errors, VarianceMethod::STANDARD_DEVIATION), 11, 4,
                            NumberFormat::SCIENTIFIC) + ")";
  }
  if (yoff) {
    fail_msg += (xoff) ? ", Y" : "Y";
    fail_msg += " (" + realToString(mean(y_errors), 11, 4, NumberFormat::SCIENTIFIC) + " +/- " +
                realToString(variance(y_errors, VarianceMethod::STANDARD_DEVIATION), 11, 4,
                             NumberFormat::SCIENTIFIC) + ")";
  }
  if (zoff) {
    fail_msg += (xoff || yoff) ? ", Z" : "Z";
    fail_msg += " (" + realToString(mean(z_errors), 11, 4, NumberFormat::SCIENTIFIC) + " +/- " +
                realToString(variance(z_errors, VarianceMethod::STANDARD_DEVIATION), 11, 4,
                             NumberFormat::SCIENTIFIC) + ")";
  }
  check((xoff || yoff || zoff) == false, "Coordinates exceeded the specified deviations after "
        "copying " + copy_msg + ".  Deviations occur in " + fail_msg + ".", do_test);
  if (overall > 1.5 * tiny) {
    check(rmsd, RelationalOperator::EQUAL, Approx(0.0).margin(overall), "The root mean squared "
          "coordinate deviation obtained after copying " + copy_msg + " does not meet "
          "expectations.", do_test);
  }
}

template <typename Tcrd, typename Tref>
void diffCoordinates(const Tcrd* xcrd, const Tcrd* ycrd, const Tcrd* zcrd, const int* xcrd_ovrf,
                     const int* ycrd_ovrf, const int* zcrd_ovrf, const double crd_scale,
                     const Tref* xref, const Tref* yref, const Tref* zref, const int* xref_ovrf,
                     const int* yref_ovrf, const int* zref_ovrf, const double ref_scale,
                     const int natom, const std::string &copy_msg, const TestPriority do_test,
                     const double overall = tiny, const double bounds = tiny) {
  diffCoordinates(xcrd, ycrd, zcrd, xcrd_ovrf, ycrd_ovrf, zcrd_ovrf, crd_scale, xref, yref, zref,
                  xref_ovrf, yref_ovrf, zref_ovrf, ref_scale, natom, nullptr, nullptr, copy_msg,
                  do_test, overall, bounds);
}

template <typename Tcrd, typename Tref>
void diffCoordinates(const Tcrd* xcrd, const Tcrd* ycrd, const Tcrd* zcrd, const double crd_scale,
                     const Tref* xref, const Tref* yref, const Tref* zref, const double ref_scale,
                     const int natom, const std::string &copy_msg, const TestPriority do_test,
                     const double overall = 0.0, const double bounds = tiny) {
  diffCoordinates(xcrd, ycrd, zcrd, nullptr, nullptr, nullptr, crd_scale, xref, yref, zref,
                  nullptr, nullptr, nullptr, ref_scale, natom, nullptr, nullptr, copy_msg, do_test,
                  overall, bounds);
}

template <typename Tcrd, typename Tref>
void diffCoordinates(const Tcrd* xcrd, const Tcrd* ycrd, const Tcrd* zcrd, const Tref* xref,
                     const Tref* yref, const Tref* zref, const int natom,
                     const std::string &copy_msg, const TestPriority do_test,
                     const double overall = 0.0, const double bounds = tiny) {
  diffCoordinates(xcrd, ycrd, zcrd, nullptr, nullptr, nullptr, 1.0, xref, yref, zref, nullptr,
                  nullptr, nullptr, 1.0, natom, nullptr, nullptr, copy_msg, do_test, overall,
                  bounds);
}

template <typename Tcrd, typename Tref>
void diffCoordinates(const Tcrd* xcrd, const Tcrd* ycrd, const Tcrd* zcrd, const double crd_scale,
                     const Tref* xref, const Tref* yref, const Tref* zref, const double ref_scale,
                     const int natom, const double* umat, const double* umat_ref,
                     const std::string &copy_msg, const TestPriority do_test,
                     const double overall = 0.0, const double bounds = tiny) {
  diffCoordinates(xcrd, ycrd, zcrd, nullptr, nullptr, nullptr, crd_scale, xref, yref, zref,
                  nullptr, nullptr, nullptr, ref_scale, natom, umat, umat_ref,
                  copy_msg, do_test, overall, bounds);
}

template <typename Tcrd, typename Tref>
void diffCoordinates(const Tcrd* xcrd, const Tcrd* ycrd, const Tcrd* zcrd, const Tref* xref,
                     const Tref* yref, const Tref* zref, const int natom, const double* umat,
                     const double* umat_ref, const std::string &copy_msg,
                     const TestPriority do_test, const double overall = 0.0,
                     const double bounds = tiny) {
  diffCoordinates(xcrd, ycrd, zcrd, nullptr, nullptr, nullptr, 1.0, xref, yref, zref, nullptr,
                  nullptr, nullptr, 1.0, natom, umat, umat_ref, copy_msg, do_test, overall,
                  bounds);
}

//-------------------------------------------------------------------------------------------------
// Replicate a series of integers.
//
// Arguments:
//   length:  The length of the original series { 0, 1, 2, ..., length - 1 }
//   nrep:    The number of times to replicate the series
//-------------------------------------------------------------------------------------------------
std::vector<int> replicateSeries(const int length, const int nrep) {
  std::vector<int> orig(length);
  for (int i = 0; i < length; i++) {
    orig[i] = i;
  }
  return tileVector(orig, nrep);
}

//-------------------------------------------------------------------------------------------------
// Create arrays of various single-system coordinate objects from the entries in a test system
// manager.
//
// Arguments:
//   tsm:          The basis test system manager, containing the coordinates to use as templates
//   tsm_cf:       Array of coordinate frames for each system, built and returned
//   tsm_cfr:      Array of abstracts for the coordinate frames, built and returned
//   tsm_ps:       Array of PhaseSpace objects for each system, built and returned
//   tsm_psr:      Array of abstracts for the phase space objects, built and returned
//   tsm_cs:       Array of coordinate series for each system, built and returned.  After the first
//                 frame, subsequent frames' atoms will be randomly perturbed.
//   tsm_csr:      Array of abstracts for the coordinate series, built and returned
//   random_seed:  Random number generator seed for the coordinate perturbations
//-------------------------------------------------------------------------------------------------
template <typename T>
PhaseSpaceSynthesis spawnMutableCoordinateObjects(const TestSystemManager &tsm,
                                                  std::vector<CoordinateFrame> *tsm_cf,
                                                  std::vector<CoordinateFrameReader> *tsm_cfr,
                                                  std::vector<PhaseSpace> *tsm_ps,
                                                  std::vector<PhaseSpaceReader> *tsm_psr,
                                                  std::vector<CoordinateSeries<T>> *tsm_cs,
                                                  std::vector<CoordinateSeriesReader<T>> *tsm_csr,
                                                  const int random_seed = 57983,
                                                  const int gpos_bits = 28,
                                                  const int vel_bits = 36,
                                                  const int frc_bits = 21) {
  Xoroshiro128pGenerator xrs(random_seed);
  tsm_cf->reserve(tsm.getSystemCount());
  tsm_ps->reserve(tsm.getSystemCount());
  tsm_cs->reserve(tsm.getSystemCount());
  const int small_mol_frames = 4;
  const int nbits = (isFloatingPointScalarType<T>()) ? 0 : 54;
  const double random_factor = pow(2.0, nbits - 4); 
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    tsm_cf->emplace_back(tsm.exportCoordinateFrame(i));
    tsm_ps->emplace_back(tsm.exportPhaseSpace(i));
    tsm_cs->emplace_back(tsm.exportCoordinateFrame(i), small_mol_frames, nbits);
  }
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    CoordinateSeriesWriter<T> icsw = tsm_cs->data()[i].data();
    for (int j = 1; j < small_mol_frames; j++) {
      int ij_start = j * roundUp(icsw.natom, warp_size_int);
      for (int k = ij_start; k < ij_start + icsw.natom; k++) {
        icsw.xcrd[k] += xrs.gaussianRandomNumber() * random_factor;
        icsw.ycrd[k] += xrs.gaussianRandomNumber() * random_factor;
        icsw.zcrd[k] += xrs.gaussianRandomNumber() * random_factor;
      }
      if (icsw.unit_cell != UnitCellType::NONE) {
        const int bdim_offset = j * roundUp(6, warp_size_int);
        const int xfrm_offset = j * roundUp(9, warp_size_int);
        for (int k = 0; k < 6; k++) {
          icsw.boxdim[bdim_offset + k] += (0.01 * xrs.uniformRandomNumber());
        }
        computeBoxTransform(&icsw.boxdim[bdim_offset], &icsw.umat[xfrm_offset],
                            &icsw.invu[xfrm_offset]);
      }
    }
  }

  // Create a vector of readers for the original, double-precision coordinates
  tsm_cfr->reserve(tsm.getSystemCount());
  tsm_psr->reserve(tsm.getSystemCount());
  tsm_csr->reserve(tsm.getSystemCount());
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    tsm_cfr->emplace_back(tsm_cf->data()[i].data());
    tsm_psr->emplace_back(tsm_ps->data()[i].data());
    tsm_csr->emplace_back(tsm_cs->data()[i].data());
  }
  
  // Return a synthesis of everything, times two
  const std::vector<int> replicator = replicateSeries(tsm.getSystemCount(), 2);
  return PhaseSpaceSynthesis(*tsm_ps, tsm.getTopologyPointer(), replicator, gpos_bits, 24,
                             vel_bits, frc_bits);
}

//-------------------------------------------------------------------------------------------------
// Check the footprint of all-to-all calculations guided by work instructions against the footprint
// generated by a straightforward application of system groupings.
//
// Arguments:
//   cg:            Instructions for comparing structures in the Condensate object
//   cdns:          The coordinates of the synthesis, distilled into floating point format
//   scmap:         Map tracing each member system of the synthesis back to the user-specified
//                  systems cache
//   organization:  The manner by which to group systems in the synthesis
//-------------------------------------------------------------------------------------------------
void checkATAFootprint(const ComparisonGuide &cg, const Condensate &cdns,
                       const SynthesisCacheMap &scmap, const SystemGrouping organization) {
  const PsSynthesisReader poly_psr = cdns.getSynthesisPointer()->data();
  const CondensateReader cdw = cdns.data();
  const SynthesisMapReader scmapr = scmap.data();
  std::vector<int> ata_occupancy(cdw.system_count * cdw.system_count, 0);
  const int* bounds_list = (organization == SystemGrouping::SOURCE) ?
                           scmapr.csystem_bounds : (organization == SystemGrouping::TOPOLOGY) ?
                           scmapr.ctopol_bounds : scmapr.clabel_bounds;
  const int* system_list = (organization == SystemGrouping::SOURCE) ?
                           scmapr.csystem_proj : (organization == SystemGrouping::TOPOLOGY) ?
                           scmapr.ctopol_proj : scmapr.clabel_proj;
  const std::vector<int4> ata_insr = cg.getATAInstructionMembers(organization);
  const std::vector<int> ata_group = cg.getATAInstructionGroups(organization);
  const int instruction_count = cg.getATAInstructionCount(organization);
  const int partition_count = scmap.getPartitionCount(organization);
  const int half_uint_bits = sizeof(uint) * 4;
  for (int i = 0; i < partition_count; i++) {
    for (int j = bounds_list[i]; j < bounds_list[i + 1]; j++) {
      const int sysj = system_list[j];
      for (int k = bounds_list[i]; k < j; k++) {
        const int sysk = system_list[k];
        ata_occupancy[(sysk * scmapr.nsynth) + sysj] = 1;
      }
    }
  }
  int half_mask = 0;
  for (int i = 0; i < half_uint_bits; i++) {
    half_mask |= (0x1 << i);
  }
  for (int i = 0; i < instruction_count; i++) {
    const int4 tinsr = ata_insr[i];
    const int tgroup = ata_group[i];
    const int x_ext = (tinsr.w & half_mask);
    const int y_ext = ((tinsr.w >> half_uint_bits) & half_mask);
    for (int j = 0; j < x_ext; j++) {
      const int jsys_idx = system_list[bounds_list[tgroup] + tinsr.x + j];
      for (int k = 0; k < y_ext; k++) {
        const int ksys_idx = system_list[bounds_list[tgroup] + tinsr.y + k];
        if (ksys_idx < jsys_idx) {
          ata_occupancy[(ksys_idx * scmapr.nsynth) + jsys_idx] -= 1;
        }
      }
    }
  }
  check(ata_occupancy, RelationalOperator::EQUAL, std::vector<int>(ata_occupancy.size(), 0),
        "All-to-all comparison instructions do not cover all systems (instances of +1), or do not "
        "cover all systems percisely once (instances of -1) when systems are gouped by " +
        getEnumerationName(organization) + ".");
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
// Prepare coordinate objects for comparison by downloading their device-side data to the host if
// the copy originated from or went to the device side.
//
// Arguments:
//   crd_a:   The first of two coordinate objects to check
//   crd_b:   The second of two coordinate objects to check
//   tier_a:  The tier containing relevant data on crd_a
//   tier_b:  The tier containing relevant data on crd_b
//-------------------------------------------------------------------------------------------------
template <typename Tcrda, typename Tcrdb>
void criticalCoordTransfer(Tcrda *crd_a, Tcrdb *crd_b, const HybridTargetLevel tier_a,
                           const HybridTargetLevel tier_b) {
  switch (tier_a) {
  case HybridTargetLevel::HOST:
    break;
  case HybridTargetLevel::DEVICE:
    crd_a->download();
  }
  switch (tier_b) {
  case HybridTargetLevel::HOST:
    break;
  case HybridTargetLevel::DEVICE:
    crd_b->download();
  }
}

//-------------------------------------------------------------------------------------------------
// Copy a coordinate set from one object to another, from a particular tier in the origin to a
// particular tier in the destination.  This routine begins by using the copy assignment operator
// to duplicate the original data so as not to corrupt it, performs the copy operation on the
// duplicates, then compares the destination of the copy procedure to the source of the original.
//
// Overloaded:
//   - Read various types of coordinate objects
//   - Copy different types of coordinates (positions, velocities, forces, alternate or primary
//     versions of each) depending on the objects
//
// Arguments:
//   crd_a:       Basis for the origin in the copy operation
//   crd_b:       Basis for the destination in the copy operation
//   dest_tier:   Specify whether to copy coordinates into the CPU host or GPU device layer
//   orig_tier:   Specify whether to copy coordinates from the CPU host or GPU device layer
//   dest_type:   The type of coordinates to which the copied information will go
//   orig_type:   The type of coordinates from which the copied information will come
//   dest_cycle:  The coordinate cycle position to which the coordinates will go
//   orig_cycle:  The coordinate cycle position from which the coordinates emerge
//   gpu:         Details of the available GPU
//   do_test:     Give the go-ahead to performing the test
//-------------------------------------------------------------------------------------------------
void copySetup(const CoordinateFrame &dest, const CoordinateFrame &orig,
               const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
               const GpuDetails &gpu, const TestPriority do_test) {

  // Copy between the two objects on the specified tiers
  CoordinateFrame active_dest = dest;
  CoordinateFrame active_orig = orig;
  CoordinateFrameWriter chkr = active_dest.data();
  CoordinateFrameWriter refr = active_orig.data();
  coordCopy(&active_dest, active_orig, dest_tier, orig_tier, gpu);

  // Prepare for comparisons in the HOST tier.
  criticalCoordTransfer(&active_dest, &active_orig, dest_tier, orig_tier);
  diffCoordinates(chkr.xcrd, chkr.ycrd, chkr.zcrd, refr.xcrd, refr.ycrd, refr.zcrd, refr.natom,
                  chkr.umat, refr.umat, "one CoordinateFrame to another, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);
}

void copySetup(const CoordinateFrame &cf, const PhaseSpace &ps,
               const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
               const GpuDetails &gpu, const TestPriority do_test) {

  // Copy from the PhaseSpace to the CoordinateFrame
  CoordinateFrame active_cf = cf;
  PhaseSpace active_ps = ps;
  CoordinateFrameWriter cfw = active_cf.data();
  PhaseSpaceWriter psw = active_ps.data();
  coordCopy(&active_cf, active_ps, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cf, &active_ps, dest_tier, orig_tier);
  diffCoordinates(cfw.xcrd, cfw.ycrd, cfw.zcrd, psw.xcrd, psw.ycrd, psw.zcrd, psw.natom,
                  cfw.umat, psw.umat, "a PhaseSpace object to a CoordinateFrame, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);

  // Reload by host-to-host copy, then upload to make the device data consistent
  Xoroshiro128pGenerator xrs(671843927);
  coordCopy(&active_cf, cf);
  std::vector<double> umat_hold(9);
  for (int i = 0; i < 9; i++) {
    umat_hold[i] = 5 * i;
    cfw.umat[i] = umat_hold[i];
  }
  active_cf.upload();
  coordCopy(&active_ps, ps);
  addRandomNoise(&xrs, psw.xfrc, psw.yfrc, psw.zfrc, psw.natom, 8.4);
  active_ps.upload();

  // Copy a different aspect of the PhaseSpace object to the CoordinateFrame.
  coordCopy(&active_cf, active_ps, TrajectoryKind::FORCES, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cf, &active_ps, dest_tier, orig_tier);
  diffCoordinates(cfw.xcrd, cfw.ycrd, cfw.zcrd, psw.xfrc, psw.yfrc, psw.zfrc, psw.natom,
                  cfw.umat, umat_hold.data(), "a PhaseSpace object (forces) to a "
                  "CoordinateFrame, " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);  
  
  // Copy from the CoordinateFrame to the PhaseSpace  
  coordCopy(&active_cf, cf);
  active_cf.upload();
  coordCopy(&active_ps, ps);
  active_ps.upload();
  coordCopy(&active_ps, active_cf, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_ps, &active_cf, dest_tier, orig_tier);
  diffCoordinates(psw.xcrd, psw.ycrd, psw.zcrd, cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.natom,
                  psw.umat, cfw.umat, "a CoordinateFrame to a PhaseSpace, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);

  // Try copying the CoordinateFrame's coordinates to a separate aspect of the PhaseSpace object,
  // at a separate point in the time cycle.  Mark up the box space transformation matrix to verify
  // that it doesn't get copied in this instance.
  coordCopy(&active_cf, cf);
  for (int i = 0; i < 9; i++) {
    umat_hold[i] = psw.umat[i];
    cfw.umat[i] = 3 * i;
  }
  active_cf.upload();
  coordCopy(&active_ps, ps);
  active_ps.upload();
  coordCopy(&active_ps, TrajectoryKind::VELOCITIES, CoordinateCycle::BLACK, active_cf,
            dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_ps, &active_cf, dest_tier, orig_tier);
  diffCoordinates(psw.vxalt, psw.vyalt, psw.vzalt, cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.natom,
                  psw.umat, umat_hold.data(), "a CoordinateFrame to a PhaseSpace (alt. "
                  "velocities), " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);

  // Try copying another aspect of the PhaseSpace object into the CoordinateFrame.  Make sure that
  // the (invalid, but distinct) transformation matrix placed in the CoordinateFrame stays intact.
  coordCopy(&active_cf, cf);
  for (int i = 0; i < 9; i++) {
    umat_hold[i] = -2 * i;
    cfw.umat[i] = umat_hold[i];
  }
  active_cf.upload();
  coordCopy(&active_ps, ps);
  addRandomNoise(&xrs, psw.fxalt, psw.fyalt, psw.fzalt, psw.natom, 7.2);
  active_ps.upload();
  coordCopy(&active_cf, active_ps, TrajectoryKind::FORCES, CoordinateCycle::BLACK,
            dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cf, &active_ps, dest_tier, orig_tier);
  diffCoordinates(cfw.xcrd, cfw.ycrd, cfw.zcrd, psw.fxalt, psw.fyalt, psw.fzalt, psw.natom,
                  cfw.umat, umat_hold.data(), "a PhaseSpace (alt. forces) to a CoordinateFrame " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);
}

template <typename T>
void copySetup(const CoordinateFrame &cf, const CoordinateSeries<T> &cs,
               const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
               const GpuDetails &gpu, const TestPriority do_test) {

  // Copy the CoordinateSeries to the CoordinateFrame
  const int frm_idx = cs.getFrameCount() / 2;
  CoordinateFrame active_cf = cf;
  CoordinateSeries<T> active_cs = cs;
  CoordinateFrameWriter cfw = active_cf.data();
  CoordinateSeriesWriter<T> csw = active_cs.data();
  coordCopy<T>(&active_cf, active_cs, frm_idx, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cf, &active_cs, dest_tier, orig_tier);
  const int atom_offset = frm_idx * roundUp(csw.natom, warp_size_int);
  const int xfrm_offset = frm_idx * roundUp(9, warp_size_int);
  diffCoordinates(cfw.xcrd, cfw.ycrd, cfw.zcrd, 1.0, &csw.xcrd[atom_offset],
                  &csw.ycrd[atom_offset], &csw.zcrd[atom_offset], csw.gpos_scale, csw.natom,
                  cfw.umat, &csw.umat[xfrm_offset], "a CoordinateSeries snapshot to a "
                  "CoordinateFrame, " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);

  // Refresh the coordinates
  coordCopy(&active_cf, cf);
  active_cf.upload();
  for (int i = 0; i < cs.getFrameCount(); i++) {
    coordCopy<T, T>(&active_cs, i, cs, i);
  }
  active_cs.upload();
  
  // Copy the CoordinateFrame into the CoordinateSeries
  coordCopy<T>(&active_cs, frm_idx, active_cf, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cs, &active_cf, dest_tier, orig_tier);
  diffCoordinates(&csw.xcrd[atom_offset], &csw.ycrd[atom_offset], &csw.zcrd[atom_offset],
                  csw.gpos_scale, cfw.xcrd, cfw.ycrd, cfw.zcrd, 1.0, cfw.natom,
                  &csw.umat[xfrm_offset], cfw.umat, "a CoordinateFrame to a snapshot in a "
                  "CoordinateSeries, " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);
}

void copySetup(const CoordinateFrame &cf, const Condensate &cdns, const int sys_idx,
               const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
               const GpuDetails &gpu, const TestPriority do_test) {

  // Copy one of the condensate's systems to the CoordinateFrame
  CoordinateFrame active_cf = cf;
  Condensate active_cdns = cdns;
  CoordinateFrameWriter cfw = active_cf.data();
  CondensateWriter cdnsw = active_cdns.data();
  coordCopy(&active_cf, active_cdns, sys_idx, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cf, &active_cdns, dest_tier, orig_tier);
  const int atom_offset = active_cdns.getAtomOffset(sys_idx);
  const int xfrm_offset = sys_idx * roundUp(9, warp_size_int);
  diffCoordinates(cfw.xcrd, cfw.ycrd, cfw.zcrd, 1.0, &cdnsw.xcrd_sp[atom_offset],
                  &cdnsw.ycrd_sp[atom_offset], &cdnsw.zcrd_sp[atom_offset], 1.0, cfw.natom,
                  cfw.umat, &cdnsw.umat[xfrm_offset], "a Condensate system to a "
                  "CoordinateFrame, " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);

  // Refresh the coordinates
  coordCopy(&active_cf, cf);
  active_cf.upload();
  for (int i = 0; i < cdns.getSystemCount(); i++) {
    coordCopy(&active_cdns, i, cdns, i);    
  }
  active_cdns.upload();

  // Copy the CoordinateFrame into one of the Condensate's systems
  coordCopy(&active_cdns, sys_idx, active_cf, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cdns, &active_cf, dest_tier, orig_tier);
  diffCoordinates(&cdnsw.xcrd_sp[atom_offset], &cdnsw.ycrd_sp[atom_offset],
                  &cdnsw.zcrd_sp[atom_offset], 1.0, cfw.xcrd, cfw.ycrd, cfw.zcrd, 1.0, cfw.natom,
                  &cdnsw.umat[xfrm_offset], cfw.umat, "a CoordinateFrame to a Condensate "
                  "system, " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test, 7.5e-7, 1.5e-6);
}

void copySetup(const CoordinateFrame &cf, const PhaseSpaceSynthesis &poly_ps, const int sys_idx,
               const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
               const GpuDetails &gpu, const TestPriority do_test) {

  // Copy one of the synthesis's systems to the CoordinateFrame
  CoordinateFrame active_cf = cf;
  PhaseSpaceSynthesis active_pps = poly_ps;
  CoordinateFrameWriter cfw = active_cf.data();
  PsSynthesisWriter poly_psw = active_pps.data();
  coordCopy(&active_cf, active_pps, sys_idx, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cf, &active_pps, dest_tier, orig_tier);
  const int atom_offset = active_pps.getAtomOffset(sys_idx);
  const int xfrm_offset = sys_idx * roundUp(9, warp_size_int);
  diffCoordinates(cfw.xcrd, cfw.ycrd, cfw.zcrd, nullptr, nullptr, nullptr, 1.0,
                  &poly_psw.xcrd[atom_offset], &poly_psw.ycrd[atom_offset],
                  &poly_psw.zcrd[atom_offset], &poly_psw.xcrd_ovrf[atom_offset],
                  &poly_psw.ycrd_ovrf[atom_offset], &poly_psw.zcrd_ovrf[atom_offset],
                  poly_psw.gpos_scale, cfw.natom, cfw.umat, &poly_psw.umat[xfrm_offset],
                  "a PhaseSpaceSynthesis system to a CoordinateFrame, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);

  // Reload and try loading a separate aspect of the PhaseSpaceSynthesis into the CoordinateFrame
  coordCopy(&active_cf, cf);
  std::vector<double> umat_hold(9);
  for (int i = 0; i < 9; i++) {
    umat_hold[i] = -2 * i;
    cfw.umat[i] = umat_hold[i];
  }
  active_cf.upload();
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    coordCopy(&active_pps, i, poly_ps, i);    
  }
  active_pps.upload();
  coordCopy(&active_cf, active_pps, sys_idx, TrajectoryKind::VELOCITIES,
            CoordinateCycle::BLACK, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cf, &active_pps, dest_tier, orig_tier);
  diffCoordinates(cfw.xcrd, cfw.ycrd, cfw.zcrd, nullptr, nullptr, nullptr, 1.0,
                  &poly_psw.vxalt[atom_offset], &poly_psw.vyalt[atom_offset],
                  &poly_psw.vzalt[atom_offset], &poly_psw.vxalt_ovrf[atom_offset],
                  &poly_psw.vyalt_ovrf[atom_offset], &poly_psw.vzalt_ovrf[atom_offset],
                  poly_psw.vel_scale, cfw.natom, cfw.umat, umat_hold.data(),
                  "a PhaseSpaceSynthesis system (alt. velocities) to a CoordinateFrame, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);

  // Reload and try copying the CoordinateFrame into the PhaseSpaceSynthesis
  coordCopy(&active_cf, cf);
  active_cf.upload();
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    coordCopy(&active_pps, i, poly_ps, i);    
  }
  active_pps.upload();
  coordCopy(&active_pps, sys_idx, active_cf, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_pps, &active_cf, dest_tier, orig_tier);
  diffCoordinates(&poly_psw.xcrd[atom_offset], &poly_psw.ycrd[atom_offset],
                  &poly_psw.zcrd[atom_offset], &poly_psw.xcrd_ovrf[atom_offset],
                  &poly_psw.ycrd_ovrf[atom_offset], &poly_psw.zcrd_ovrf[atom_offset],
                  poly_psw.gpos_scale, cfw.xcrd, cfw.ycrd, cfw.zcrd, nullptr, nullptr, nullptr,
                  1.0, cfw.natom, cfw.umat, &poly_psw.umat[xfrm_offset],
                  "a CoordinateFrame to a PhaseSpaceSynthesis system, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test,
                  1.2e-8, 2.5e-8);

  // Try copying the CoordinateFrame into another aspect of the PhaseSpaceSynthesis
  coordCopy(&active_cf, cf);
  for (int i = 0; i < 9; i++) {
    umat_hold[i] = -i;
    cfw.umat[i] = umat_hold[i];
  }
  active_cf.upload();
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    coordCopy(&active_pps, i, poly_ps, i);    
  }
  active_pps.upload();
  coordCopy(&active_pps, sys_idx, TrajectoryKind::VELOCITIES, CoordinateCycle::BLACK,
            active_cf, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_pps, &active_cf, dest_tier, orig_tier);
  diffCoordinates(&poly_psw.vxalt[atom_offset], &poly_psw.vyalt[atom_offset],
                  &poly_psw.vzalt[atom_offset], &poly_psw.vxalt_ovrf[atom_offset],
                  &poly_psw.vyalt_ovrf[atom_offset], &poly_psw.vzalt_ovrf[atom_offset],
                  poly_psw.vel_scale, cfw.xcrd, cfw.ycrd, cfw.zcrd, nullptr, nullptr, nullptr,
                  1.0, cfw.natom, cfw.umat, umat_hold.data(),
                  "a CoordinateFrame to a PhaseSpaceSynthesis system (alt. velocities), " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test,
                  3.0e-9, 6.2e-9);
}

void copySetup(const PhaseSpace &dest, const PhaseSpace &orig,
               const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
               const GpuDetails &gpu, const TestPriority do_test) {
  
  // Duplicate the inputs and add random forces and velocities to each.
  Xoroshiro128pGenerator xrs(817252330);
  PhaseSpace active_dest = dest;
  PhaseSpace active_orig = orig;
  PhaseSpaceWriter chkw = active_dest.data();
  PhaseSpaceWriter refw = active_orig.data();
  addRandomNoise(&xrs, chkw.xvel, chkw.yvel, chkw.zvel, chkw.natom, 5.0);
  addRandomNoise(&xrs, chkw.xfrc, chkw.yfrc, chkw.zfrc, chkw.natom, 9.0);
  addRandomNoise(&xrs, refw.xvel, refw.yvel, refw.zvel, refw.natom, 2.0);
  addRandomNoise(&xrs, refw.xfrc, refw.yfrc, refw.zfrc, refw.natom, 4.0);
  active_dest.upload();
  active_orig.upload();
  
  // Copy between the two objects on the specified tiers
  coordCopy(&active_dest, active_orig, dest_tier, orig_tier, gpu);

  // Prepare for comparisons in the HOST tier.
  criticalCoordTransfer(&active_dest, &active_orig, dest_tier, orig_tier);
  diffCoordinates(chkw.xcrd, chkw.ycrd, chkw.zcrd, refw.xcrd, refw.ycrd, refw.zcrd, refw.natom,
                  chkw.umat, refw.umat, "one PhaseSpace object to another, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);
  diffCoordinates(chkw.xvel, chkw.yvel, chkw.zvel, refw.xvel, refw.yvel, refw.zvel, refw.natom,
                  chkw.umat, refw.umat, "one PhaseSpace object to another, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);
  diffCoordinates(chkw.xfrc, chkw.yfrc, chkw.zfrc, refw.xfrc, refw.yfrc, refw.zfrc, refw.natom,
                  chkw.umat, refw.umat, "one PhaseSpace object to another, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);
}

template <typename T>
void copySetup(const PhaseSpace &ps, const CoordinateSeries<T> &cs,
               const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
               const GpuDetails &gpu, const TestPriority do_test) {

  // Duplicate the inputs.  Add random forces and velocities to the PhaseSpace object.
  Xoroshiro128pGenerator xrs(30184725);
  PhaseSpace active_ps = ps;
  CoordinateSeries<T> active_cs = cs;
  PhaseSpaceWriter psw = active_ps.data();
  CoordinateSeriesWriter<T> csw = active_cs.data();
  active_ps.upload();
  active_cs.upload();

  // Copy one frame of the coordinate series into the PhaseSpace's coordinate holdings.
  const int frm_idx = cs.getFrameCount() - 1;
  coordCopy(&active_ps, active_cs, frm_idx, dest_tier, orig_tier, gpu);
  const int atom_offset = roundUp(active_cs.getAtomCount(), warp_size_int) * frm_idx;
  const int xfrm_offset = roundUp(9, warp_size_int) * frm_idx;
  criticalCoordTransfer(&active_ps, &active_cs, dest_tier, orig_tier);
  diffCoordinates(psw.xcrd, psw.ycrd, psw.zcrd, &csw.xcrd[atom_offset], &csw.ycrd[atom_offset],
                  &csw.zcrd[atom_offset], csw.natom, psw.umat, &csw.umat[xfrm_offset],
                  "a CoordinateSeries snapshot to a PhaseSpace object, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);

  // Transfer the coordinate series data to an alternate aspect of the PhaseSpace object.
  addRandomNoise(&xrs, csw.xcrd, csw.ycrd, csw.zcrd, csw.natom, 4.1);
  active_cs.upload();
  std::vector<double> umat_hold(9);
  for (int i = 0; i < 9; i++) {
    umat_hold[i] = 5 - (3 * i);
    psw.umat[i] = umat_hold[i];
  }
  active_ps.upload();
  coordCopy(&active_ps, TrajectoryKind::FORCES, CoordinateCycle::BLACK, active_cs, 0,
           dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_ps, &active_cs, dest_tier, orig_tier);
  diffCoordinates(psw.fxalt, psw.fyalt, psw.fzalt, csw.xcrd, csw.ycrd, csw.zcrd, csw.natom,
                  psw.umat, umat_hold.data(), "a CoordinateSeries snapshot to a PhaseSpace "
                  "object, " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);

  // Reload both objects and try copying the PhaseSpace data into the CoordinateSeries.
  coordCopy(&active_ps, ps);
  addRandomNoise(&xrs, psw.xvel, psw.yvel, psw.zvel, psw.natom, 6.7);
  addRandomNoise(&xrs, psw.xfrc, psw.yfrc, psw.zfrc, psw.natom, 6.7);
  active_ps.upload();
  for (int i = 0; i < cs.getFrameCount(); i++) {
    coordCopy<T, T>(&active_cs, i, cs, i);
  }
  active_cs.upload();
  const int nfrm_idx = cs.getFrameCount() / 2;
  coordCopy(&active_cs, nfrm_idx, active_ps, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cs, &active_ps, dest_tier, orig_tier);
  const int natom_offset = roundUp(active_cs.getAtomCount(), warp_size_int) * nfrm_idx;
  const int nxfrm_offset = roundUp(9, warp_size_int) * nfrm_idx;
  diffCoordinates(&csw.xcrd[natom_offset], &csw.ycrd[natom_offset], &csw.zcrd[natom_offset],
                  psw.xcrd, psw.ycrd, psw.zcrd, csw.natom, &csw.umat[nxfrm_offset], psw.umat,
                  "a PhaseSpace object's positions into a CoordinateSeries snapshot, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);
  for (int i = 0; i < 9; i++) {
    umat_hold[i] = 7 - (2 * i);
    csw.umat[i] = umat_hold[i];
  }
  active_cs.upload();
  coordCopy(&active_cs, 0, active_ps, TrajectoryKind::VELOCITIES, CoordinateCycle::BLACK,
            dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cs, &active_ps, dest_tier, orig_tier);
  diffCoordinates(csw.xcrd, csw.ycrd, csw.zcrd, psw.vxalt, psw.vyalt, psw.vzalt, csw.natom,
                  csw.umat, umat_hold.data(), "a PhaseSpace object into a CoordinateSeries "
                  "snapshot (alt. velocities), " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);
}

void copySetup(const PhaseSpace &ps, const Condensate &cdns, const int sys_idx,
               const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
               const GpuDetails &gpu, const TestPriority do_test) {
  
  // Duplicate the inputs. 
  Xoroshiro128pGenerator xrs(7719382);
  PhaseSpace active_ps = ps;
  Condensate active_cdns = cdns;
  PhaseSpaceWriter psw = active_ps.data();
  CondensateWriter cdnsw = active_cdns.data();
  active_ps.upload();
  active_cdns.upload();

  // Copy one system from the condensate into the PhaseSpace object.
  coordCopy(&active_ps, active_cdns, sys_idx, dest_tier, orig_tier, gpu);
  const size_t atom_offset = active_cdns.getAtomOffset(sys_idx);
  const int xfrm_offset = roundUp(9, warp_size_int) * sys_idx;
  criticalCoordTransfer(&active_ps, &active_cdns, dest_tier, orig_tier);
  diffCoordinates(psw.xcrd, psw.ycrd, psw.zcrd, &cdnsw.xcrd_sp[atom_offset],
                  &cdnsw.ycrd_sp[atom_offset], &cdnsw.zcrd_sp[atom_offset], psw.natom, psw.umat,
                  &cdnsw.umat[xfrm_offset], "a Condensate into a PhaseSpace object, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);

  // Copy the system from the condensate into an alternate aspect of the PhaseSpace object.
  addRandomNoise(&xrs, &cdnsw.xcrd_sp[atom_offset], &cdnsw.ycrd_sp[atom_offset],
                 &cdnsw.zcrd_sp[atom_offset], cdns.getAtomCount(sys_idx), 1.7);
  active_cdns.upload();
  std::vector<double> umat_hold(9);
  for (int i = 0; i < 9; i++) {
    umat_hold[i] = xrs.gaussianRandomNumber();
    psw.umat[i] = umat_hold[i];
  }
  active_ps.upload();
  coordCopy(&active_ps, TrajectoryKind::VELOCITIES, CoordinateCycle::BLACK, active_cdns,
            sys_idx, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_ps, &active_cdns, dest_tier, orig_tier);
  diffCoordinates(psw.vxalt, psw.vyalt, psw.vzalt, &cdnsw.xcrd_sp[atom_offset],
                  &cdnsw.ycrd_sp[atom_offset], &cdnsw.zcrd_sp[atom_offset], psw.natom, psw.umat,
                  umat_hold.data(), "a Condensate into a PhaseSpace object, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);

  // Reload and copy the PhaseSpace object to one system in the Condensate
  coordCopy(&active_ps, ps);
  active_ps.upload();
  for (int i = 0; i < cdns.getSystemCount(); i++) {
    coordCopy(&active_cdns, i, cdns, i);
  }
  active_cdns.upload();
  coordCopy(&active_cdns, sys_idx, active_ps, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cdns, &active_ps, dest_tier, orig_tier);
  diffCoordinates(&cdnsw.xcrd_sp[atom_offset], &cdnsw.ycrd_sp[atom_offset],
                  &cdnsw.zcrd_sp[atom_offset], psw.xcrd, psw.ycrd, psw.zcrd, psw.natom,
                  &cdnsw.umat[xfrm_offset], psw.umat, "a PhaseSpace object into a Condensate, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test,
                  7.5e-7, 1.5e-6);

  // Copy information from another aspect of the PhaseSpace object into the Condensate
  addRandomNoise(&xrs, psw.fxalt, psw.fyalt, psw.fzalt, psw.natom, 5.4);
  active_ps.upload();
  for (int i = 0; i < 9; i++) {
    umat_hold[i] = xrs.gaussianRandomNumber();
    cdnsw.umat[xfrm_offset + i] = umat_hold[i];
  }
  active_cdns.upload();
  coordCopy(&active_cdns, sys_idx, active_ps, TrajectoryKind::FORCES, CoordinateCycle::BLACK,
            dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cdns, &active_ps, dest_tier, orig_tier);
  diffCoordinates(&cdnsw.xcrd_sp[atom_offset], &cdnsw.ycrd_sp[atom_offset],
                  &cdnsw.zcrd_sp[atom_offset], psw.fxalt, psw.fyalt, psw.fzalt, psw.natom,
                  &cdnsw.umat[xfrm_offset], umat_hold.data(), "a PhaseSpace object (alt. forces) "
                  "into a Condensate, " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test, 7.5e-7, 1.5e-6);  
}

void copySetup(const PhaseSpace &ps, const PhaseSpaceSynthesis &poly_ps, const int sys_idx,
               const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
               const GpuDetails &gpu, const TestPriority do_test) {
  
  // Duplicate the inputs and copy one system from the synthesis into the PhaseSpace object.
  Xoroshiro128pGenerator xrs(7719382);
  PhaseSpace active_ps = ps;
  PhaseSpaceSynthesis active_pps = poly_ps;
  PhaseSpaceWriter psw = active_ps.data();
  PsSynthesisWriter poly_psw = active_pps.data();
  coordCopy(&active_ps, active_pps, sys_idx, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_ps, &active_pps, dest_tier, orig_tier);
  const int atom_offset = active_pps.getAtomOffset(sys_idx);
  const int xfrm_offset = roundUp(9, warp_size_int) * sys_idx;
  diffCoordinates(psw.xcrd, psw.ycrd, psw.zcrd, nullptr, nullptr, nullptr, 1.0,
                  &poly_psw.xcrd[atom_offset], &poly_psw.ycrd[atom_offset],
                  &poly_psw.zcrd[atom_offset], &poly_psw.ycrd_ovrf[atom_offset],
                  &poly_psw.xcrd_ovrf[atom_offset], &poly_psw.zcrd_ovrf[atom_offset],
                  poly_psw.gpos_scale, psw.natom, psw.umat, &poly_psw.umat[xfrm_offset],
                  "a PhaseSpaceSynthesis (primary positions) into a PhaseSpace object, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);
  diffCoordinates(psw.xvel, psw.yvel, psw.zvel, nullptr, nullptr, nullptr, 1.0,
                  &poly_psw.xvel[atom_offset], &poly_psw.yvel[atom_offset],
                  &poly_psw.zvel[atom_offset], &poly_psw.yvel_ovrf[atom_offset],
                  &poly_psw.xvel_ovrf[atom_offset], &poly_psw.zvel_ovrf[atom_offset],
                  poly_psw.vel_scale, psw.natom, psw.umat, &poly_psw.umat[xfrm_offset],
                  "a PhaseSpaceSynthesis (primary velocities) into a PhaseSpace object, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);
  diffCoordinates(psw.xfrc, psw.yfrc, psw.zfrc, nullptr, nullptr, nullptr, 1.0,
                  &poly_psw.xfrc[atom_offset], &poly_psw.yfrc[atom_offset],
                  &poly_psw.zfrc[atom_offset], &poly_psw.yfrc_ovrf[atom_offset],
                  &poly_psw.xfrc_ovrf[atom_offset], &poly_psw.zfrc_ovrf[atom_offset],
                  poly_psw.frc_scale, psw.natom, psw.umat, &poly_psw.umat[xfrm_offset],
                  "a PhaseSpaceSynthesis (primary forces) into a PhaseSpace object, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);
  diffCoordinates(psw.xalt, psw.yalt, psw.zalt, nullptr, nullptr, nullptr, 1.0,
                  &poly_psw.xalt[atom_offset], &poly_psw.yalt[atom_offset],
                  &poly_psw.zalt[atom_offset], &poly_psw.yalt_ovrf[atom_offset],
                  &poly_psw.xalt_ovrf[atom_offset], &poly_psw.zalt_ovrf[atom_offset],
                  poly_psw.gpos_scale, psw.natom, psw.umat, &poly_psw.umat[xfrm_offset],
                  "a PhaseSpaceSynthesis (alternate positions) into a PhaseSpace object, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);
  diffCoordinates(psw.vxalt, psw.vyalt, psw.vzalt, nullptr, nullptr, nullptr, 1.0,
                  &poly_psw.vxalt[atom_offset], &poly_psw.vyalt[atom_offset],
                  &poly_psw.vzalt[atom_offset], &poly_psw.vyalt_ovrf[atom_offset],
                  &poly_psw.vxalt_ovrf[atom_offset], &poly_psw.vzalt_ovrf[atom_offset],
                  poly_psw.vel_scale, psw.natom, psw.umat, &poly_psw.umat[xfrm_offset],
                  "a PhaseSpaceSynthesis (alternate positions) into a PhaseSpace object, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);
  diffCoordinates(psw.fxalt, psw.fyalt, psw.fzalt, nullptr, nullptr, nullptr, 1.0,
                  &poly_psw.fxalt[atom_offset], &poly_psw.fyalt[atom_offset],
                  &poly_psw.fzalt[atom_offset], &poly_psw.fyalt_ovrf[atom_offset],
                  &poly_psw.fxalt_ovrf[atom_offset], &poly_psw.fzalt_ovrf[atom_offset],
                  poly_psw.frc_scale, psw.natom, psw.umat, &poly_psw.umat[xfrm_offset],
                  "a PhaseSpaceSynthesis (alternate positions) into a PhaseSpace object, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);

  // Copy a system from the synthesis into the single-structure object.
  coordCopy(&active_ps, ps);
  active_ps.upload();
  coordCopy(&active_pps, sys_idx, poly_ps, sys_idx);
  active_pps.upload();
  coordCopy(&active_pps, sys_idx, active_ps, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_pps, &active_ps, dest_tier, orig_tier);
  diffCoordinates(&poly_psw.xcrd[atom_offset], &poly_psw.ycrd[atom_offset],
                  &poly_psw.zcrd[atom_offset], &poly_psw.ycrd_ovrf[atom_offset],
                  &poly_psw.xcrd_ovrf[atom_offset], &poly_psw.zcrd_ovrf[atom_offset],
                  poly_psw.gpos_scale, psw.xcrd, psw.ycrd, psw.zcrd, nullptr, nullptr, nullptr,
                  1.0, psw.natom, &poly_psw.umat[xfrm_offset], psw.umat,
                  "a PhaseSpaceSynthesis (primary positions) into a PhaseSpace object, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test,
                  poly_psw.inv_gpos_scale, poly_psw.inv_gpos_scale);
  diffCoordinates(&poly_psw.xvel[atom_offset], &poly_psw.yvel[atom_offset],
                  &poly_psw.zvel[atom_offset], &poly_psw.yvel_ovrf[atom_offset],
                  &poly_psw.xvel_ovrf[atom_offset], &poly_psw.zvel_ovrf[atom_offset],
                  poly_psw.vel_scale, psw.xvel, psw.yvel, psw.zvel, nullptr, nullptr, nullptr, 1.0,
                  psw.natom, &poly_psw.umat[xfrm_offset], psw.umat,
                  "a PhaseSpaceSynthesis (primary velocities) into a PhaseSpace object, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test,
                  poly_psw.inv_vel_scale, poly_psw.inv_vel_scale);
  diffCoordinates(&poly_psw.xfrc[atom_offset], &poly_psw.yfrc[atom_offset],
                  &poly_psw.zfrc[atom_offset], &poly_psw.yfrc_ovrf[atom_offset],
                  &poly_psw.xfrc_ovrf[atom_offset], &poly_psw.zfrc_ovrf[atom_offset],
                  poly_psw.frc_scale, psw.xfrc, psw.yfrc, psw.zfrc, nullptr, nullptr, nullptr, 1.0,
                  psw.natom, &poly_psw.umat[xfrm_offset], psw.umat,
                  "a PhaseSpaceSynthesis (primary forces) into a PhaseSpace object, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test,
                  poly_psw.inv_frc_scale, poly_psw.inv_frc_scale);
  diffCoordinates(&poly_psw.xalt[atom_offset], &poly_psw.yalt[atom_offset],
                  &poly_psw.zalt[atom_offset], &poly_psw.yalt_ovrf[atom_offset],
                  &poly_psw.xalt_ovrf[atom_offset], &poly_psw.zalt_ovrf[atom_offset],
                  poly_psw.gpos_scale, psw.xalt, psw.yalt, psw.zalt, nullptr, nullptr, nullptr,
                  1.0, psw.natom, &poly_psw.umat[xfrm_offset], psw.umat,
                  "a PhaseSpaceSynthesis (alternate positions) into a PhaseSpace object, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test,
                  poly_psw.inv_gpos_scale, poly_psw.inv_gpos_scale);
  diffCoordinates(&poly_psw.vxalt[atom_offset], &poly_psw.vyalt[atom_offset],
                  &poly_psw.vzalt[atom_offset], &poly_psw.vyalt_ovrf[atom_offset],
                  &poly_psw.vxalt_ovrf[atom_offset], &poly_psw.vzalt_ovrf[atom_offset],
                  poly_psw.vel_scale, psw.vxalt, psw.vyalt, psw.vzalt, nullptr, nullptr, nullptr,
                  1.0, psw.natom, &poly_psw.umat[xfrm_offset], psw.umat,
                  "a PhaseSpaceSynthesis (alternate positions) into a PhaseSpace object, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test,
                  poly_psw.inv_vel_scale, poly_psw.inv_vel_scale);
  diffCoordinates(&poly_psw.fxalt[atom_offset], &poly_psw.fyalt[atom_offset],
                  &poly_psw.fzalt[atom_offset], &poly_psw.fyalt_ovrf[atom_offset],
                  &poly_psw.fxalt_ovrf[atom_offset], &poly_psw.fzalt_ovrf[atom_offset],
                  poly_psw.frc_scale, psw.fxalt, psw.fyalt, psw.fzalt, nullptr, nullptr, nullptr,
                  1.0, psw.natom, &poly_psw.umat[xfrm_offset], psw.umat,
                  "a PhaseSpaceSynthesis (alternate positions) into a PhaseSpace object, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test,
                  poly_psw.inv_frc_scale, poly_psw.inv_frc_scale);
}

template <typename Tdest, typename Torig>
void copySetup(const CoordinateSeries<Tdest> &dest, const CoordinateSeries<Torig> &orig,
               const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
               const GpuDetails &gpu, const TestPriority do_test) {

  // Duplicate the inputs and copy one into the other.
  CoordinateSeries<Tdest> active_dest = dest;
  CoordinateSeries<Torig> active_orig = orig;
  CoordinateSeriesWriter<Tdest> destw = active_dest.data();
  CoordinateSeriesWriter<Torig> origw = active_orig.data();
  const int frame_dest = active_dest.getFrameCount() - 1;
  const int frame_orig = std::max(0, active_orig.getFrameCount() - 2);
  const int natom = active_orig.getAtomCount();
  const int padded_natom = roundUp(natom, warp_size_int);
  const int padded_xfrm  = roundUp(9, warp_size_int);
  const int orig_atom_offset = padded_natom * frame_orig;
  const int dest_atom_offset = padded_natom * frame_dest;
  const int orig_xfrm_offset = padded_xfrm  * frame_orig;
  const int dest_xfrm_offset = padded_xfrm  * frame_dest;
  coordCopy<Tdest, Torig>(&active_dest, frame_dest, active_orig, frame_orig, dest_tier, orig_tier,
                          gpu);
  criticalCoordTransfer(&active_dest, &active_orig, dest_tier, orig_tier);
  diffCoordinates(&destw.xcrd[dest_atom_offset], &destw.ycrd[dest_atom_offset],
                  &destw.zcrd[dest_atom_offset], destw.gpos_scale, &origw.xcrd[orig_atom_offset],
                  &origw.ycrd[orig_atom_offset], &origw.zcrd[orig_atom_offset], origw.gpos_scale,
                  natom, &destw.umat[dest_xfrm_offset], &origw.umat[orig_xfrm_offset],
                  "one CoordinateSeries snapshot into another, " + getEnumerationName(orig_tier) +
                  " to " + getEnumerationName(dest_tier), do_test);
}

template <typename T>
void copySetup(const CoordinateSeries<T> &cs, const PhaseSpaceSynthesis &poly_ps,
               const int sys_idx, const HybridTargetLevel dest_tier,
               const HybridTargetLevel orig_tier, const GpuDetails &gpu,
               const TestPriority do_test) {

  // Duplicate the inputs and copy the PhaseSpaceSynthesis positions into the CoordinateSeries at
  // a selected frame.
  CoordinateSeries<T> active_cs = cs;
  PhaseSpaceSynthesis active_pps = poly_ps;
  CoordinateSeriesWriter<T> csw = active_cs.data();
  PsSynthesisWriter poly_psw = active_pps.data();
  const int frm_idx = std::max(0, cs.getFrameCount() - 2);
  const int natom = cs.getAtomCount();
  const int padded_xfrm = roundUp(9, warp_size_int);
  const int cs_atom_offset = roundUp(natom, warp_size_int) * frm_idx;
  const int cs_xfrm_offset = padded_xfrm * frm_idx;
  const int pps_atom_offset = poly_ps.getAtomOffset(sys_idx);
  const int pps_xfrm_offset = padded_xfrm * sys_idx;
  coordCopy(&active_cs, frm_idx, active_pps, sys_idx, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cs, &active_pps, dest_tier, orig_tier);
  diffCoordinates(&csw.xcrd[cs_atom_offset], &csw.ycrd[cs_atom_offset],
                  &csw.zcrd[cs_atom_offset], nullptr, nullptr, nullptr, csw.gpos_scale,
                  &poly_psw.xcrd[pps_atom_offset], &poly_psw.ycrd[pps_atom_offset],
                  &poly_psw.zcrd[pps_atom_offset], &poly_psw.xcrd_ovrf[pps_atom_offset],
                  &poly_psw.ycrd_ovrf[pps_atom_offset], &poly_psw.zcrd_ovrf[pps_atom_offset],
                  poly_psw.gpos_scale, natom, &csw.umat[cs_xfrm_offset],
                  &poly_psw.umat[pps_xfrm_offset], "a PhaseSpaceSynthesis system into a "
                  "CoordinateSeries snapshot, " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);

  // Copy a different aspect of the PhaseSpaceSynthesis
  Xoroshiro128pGenerator xrs(7719382);
  addRandomNoise(&xrs, &poly_psw.xalt[pps_atom_offset], &poly_psw.xalt_ovrf[pps_atom_offset],
                 &poly_psw.yalt[pps_atom_offset], &poly_psw.yalt_ovrf[pps_atom_offset],
                 &poly_psw.zalt[pps_atom_offset], &poly_psw.zalt_ovrf[pps_atom_offset], natom,
                 4.9, poly_psw.gpos_scale);
  for (int i = 0; i < 9; i++) {
    poly_psw.umat_alt[pps_xfrm_offset + i] = xrs.gaussianRandomNumber();
  }
  active_pps.upload();
  coordCopy(&active_cs, frm_idx, active_pps, sys_idx, TrajectoryKind::POSITIONS,
            CoordinateCycle::BLACK, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cs, &active_pps, dest_tier, orig_tier);
  diffCoordinates(&csw.xcrd[cs_atom_offset], &csw.ycrd[cs_atom_offset],
                  &csw.zcrd[cs_atom_offset], nullptr, nullptr, nullptr, csw.gpos_scale,
                  &poly_psw.xalt[pps_atom_offset], &poly_psw.yalt[pps_atom_offset],
                  &poly_psw.zalt[pps_atom_offset], &poly_psw.xalt_ovrf[pps_atom_offset],
                  &poly_psw.yalt_ovrf[pps_atom_offset], &poly_psw.zalt_ovrf[pps_atom_offset],
                  poly_psw.gpos_scale, natom, &csw.umat[cs_xfrm_offset],
                  &poly_psw.umat_alt[pps_xfrm_offset], "a PhaseSpaceSynthesis system into a "
                  "CoordinateSeries snapshot, " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);

  // Reload and copy form the CoordinateSeries to the PhaseSpaceSynthesis
  for (int i = 0; i < cs.getFrameCount(); i++) {
    coordCopy<T, T>(&active_cs, i, cs, i);
  }
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    coordCopy(&active_pps, i, poly_ps, i);
  }
  active_cs.upload();
  active_pps.upload();
  coordCopy(&active_pps, sys_idx, active_cs, frm_idx, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_pps, &active_cs, dest_tier, orig_tier);
  diffCoordinates(&poly_psw.xcrd[pps_atom_offset], &poly_psw.ycrd[pps_atom_offset],
                  &poly_psw.zcrd[pps_atom_offset], &poly_psw.xcrd_ovrf[pps_atom_offset],
                  &poly_psw.ycrd_ovrf[pps_atom_offset], &poly_psw.zcrd_ovrf[pps_atom_offset],
                  poly_psw.gpos_scale, &csw.xcrd[cs_atom_offset], &csw.ycrd[cs_atom_offset],
                  &csw.zcrd[cs_atom_offset], nullptr, nullptr, nullptr, csw.gpos_scale,
                  natom, &poly_psw.umat[pps_xfrm_offset], &csw.umat[cs_xfrm_offset],
                  "a CoordinateSeries snapshot into a PhaseSpaceSynthesis system, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test,
                  7.5e-7, 1.5e-6);
  addRandomNoise(&xrs, &csw.xcrd[cs_atom_offset], &csw.ycrd[cs_atom_offset],
                 &csw.zcrd[cs_atom_offset], natom, 3.1);
  active_cs.upload();
  std::vector<double> umat_hold(9);
  for (int i = 0; i < 9; i++) {
    umat_hold[i] = xrs.gaussianRandomNumber();
    poly_psw.umat_alt[pps_xfrm_offset + i] = umat_hold[i];
  }
  active_pps.upload();
  coordCopy(&active_pps, sys_idx, TrajectoryKind::FORCES, CoordinateCycle::BLACK, active_cs,
            frm_idx, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_pps, &active_cs, dest_tier, orig_tier);
  diffCoordinates(&poly_psw.fxalt[pps_atom_offset], &poly_psw.fyalt[pps_atom_offset],
                  &poly_psw.fzalt[pps_atom_offset], &poly_psw.fxalt_ovrf[pps_atom_offset],
                  &poly_psw.fyalt_ovrf[pps_atom_offset], &poly_psw.fzalt_ovrf[pps_atom_offset],
                  poly_psw.frc_scale, &csw.xcrd[cs_atom_offset], &csw.ycrd[cs_atom_offset],
                  &csw.zcrd[cs_atom_offset], nullptr, nullptr, nullptr, csw.gpos_scale, natom,
                  &poly_psw.umat_alt[pps_xfrm_offset], umat_hold.data(),
                  "a CoordinateSeries snapshot into a PhaseSpaceSynthesis system (alt. forces), " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test,
                  3.2e-7, 1.5e-6);
}

template <typename T>
void copySetup(const CoordinateSeries<T> &cs, const Condensate &cdns, const int sys_idx,
               const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
               const GpuDetails &gpu, const TestPriority do_test) {

  // Duplicate the inputs and copy the PhaseSpaceSynthesis positions into the CoordinateSeries at
  // a selected frame.
  CoordinateSeries<T> active_cs = cs;
  Condensate active_cdns = cdns;
  CoordinateSeriesWriter<T> csw = active_cs.data();
  CondensateWriter cdnsw = active_cdns.data();
  const int frm_idx = std::max(0, cs.getFrameCount() - 3);
  const int natom = cs.getAtomCount();
  const int padded_xfrm = roundUp(9, warp_size_int);
  const int cs_atom_offset = roundUp(natom, warp_size_int) * frm_idx;
  const int cs_xfrm_offset = padded_xfrm * frm_idx;
  const int cdns_atom_offset = cdns.getAtomOffset(sys_idx);
  const int cdns_xfrm_offset = padded_xfrm * sys_idx;
  coordCopy(&active_cs, frm_idx, active_cdns, sys_idx, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cs, &active_cdns, dest_tier, orig_tier);
  diffCoordinates(&csw.xcrd[cs_atom_offset], &csw.ycrd[cs_atom_offset], &csw.zcrd[cs_atom_offset],
                  csw.gpos_scale, &cdnsw.xcrd_sp[cdns_atom_offset],
                  &cdnsw.ycrd_sp[cdns_atom_offset], &cdnsw.zcrd_sp[cdns_atom_offset], 1.0, natom,
                  &csw.umat[cs_xfrm_offset], &cdnsw.umat[cdns_xfrm_offset], "a Condensate system "
                  "into a CoordinateSeries snapshot, " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);

  // Reload and perform the reverse copy operation.
  for (int i = 0; i < cs.getFrameCount(); i++) {
    coordCopy<T, T>(&active_cs, i, cs, i);
  }
  for (int i = 0; i < cdns.getSystemCount(); i++) {
    coordCopy(&active_cdns, i, cdns, i);
  }
  active_cs.upload();
  active_cdns.upload();
  coordCopy(&active_cdns, sys_idx, active_cs, frm_idx, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cdns, &active_cs, dest_tier, orig_tier);  
  diffCoordinates(&cdnsw.xcrd_sp[cdns_atom_offset], &cdnsw.ycrd_sp[cdns_atom_offset],
                  &cdnsw.zcrd_sp[cdns_atom_offset], 1.0, &csw.xcrd[cs_atom_offset],
                  &csw.ycrd[cs_atom_offset], &csw.zcrd[cs_atom_offset], csw.gpos_scale, natom,
                  &cdnsw.umat[cdns_xfrm_offset], &csw.umat[cs_xfrm_offset], "a CoordinateSeries "
                  "snapshot into a Condensate system, " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test, 7.5e-7, 1.0e-6);
}

void copySetup(const PhaseSpaceSynthesis &poly_ps, const int psys_idx, const Condensate &cdns,
               const int csys_idx, const HybridTargetLevel dest_tier,
               const HybridTargetLevel orig_tier, const GpuDetails &gpu,
               const TestPriority do_test) {

  // Duplicate the inputs and copy the PhaseSpaceSynthesis positions into the CoordinateSeries at
  // a selected frame.
  PhaseSpaceSynthesis active_pps = poly_ps;
  Condensate active_cdns = cdns;
  PsSynthesisWriter poly_psw = active_pps.data();
  CondensateWriter cdnsw = active_cdns.data();
  const int natom = cdns.getAtomCount(csys_idx);
  const int padded_xfrm = roundUp(9, warp_size_int);
  const int pps_atom_offset = poly_ps.getAtomOffset(psys_idx);
  const int pps_xfrm_offset = padded_xfrm * psys_idx;
  const int cdns_atom_offset = cdns.getAtomOffset(csys_idx);
  const int cdns_xfrm_offset = padded_xfrm * csys_idx;
  coordCopy(&active_pps, psys_idx, active_cdns, csys_idx, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_pps, &active_cdns, dest_tier, orig_tier);  
  diffCoordinates(&poly_psw.xcrd[pps_atom_offset], &poly_psw.ycrd[pps_atom_offset],
                  &poly_psw.zcrd[pps_atom_offset], &poly_psw.xcrd_ovrf[pps_atom_offset],
                  &poly_psw.ycrd_ovrf[pps_atom_offset], &poly_psw.zcrd_ovrf[pps_atom_offset],
                  poly_psw.gpos_scale, &cdnsw.xcrd_sp[cdns_atom_offset],
                  &cdnsw.ycrd_sp[cdns_atom_offset], &cdnsw.zcrd_sp[cdns_atom_offset], nullptr,
                  nullptr, nullptr, 1.0, natom, &poly_psw.umat[pps_xfrm_offset],
                  &cdnsw.umat[cdns_xfrm_offset], "a Condensate system into a "
                  "PhaseSpaceSynthesis, " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);

  // Load the Condensate data into different aspects of the PhaseSpaceSynthesis.
  coordCopy(&active_pps, psys_idx, TrajectoryKind::POSITIONS, CoordinateCycle::BLACK,
            active_cdns, csys_idx, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_pps, &active_cdns, dest_tier, orig_tier);  
  diffCoordinates(&poly_psw.xalt[pps_atom_offset], &poly_psw.yalt[pps_atom_offset],
                  &poly_psw.zalt[pps_atom_offset], &poly_psw.xalt_ovrf[pps_atom_offset],
                  &poly_psw.yalt_ovrf[pps_atom_offset], &poly_psw.zalt_ovrf[pps_atom_offset],
                  poly_psw.gpos_scale, &cdnsw.xcrd_sp[cdns_atom_offset],
                  &cdnsw.ycrd_sp[cdns_atom_offset], &cdnsw.zcrd_sp[cdns_atom_offset], nullptr,
                  nullptr, nullptr, 1.0, natom, &poly_psw.umat_alt[pps_xfrm_offset],
                  &cdnsw.umat[cdns_xfrm_offset], "a Condensate system into a "
                  "PhaseSpaceSynthesis (alt. positions), " + getEnumerationName(orig_tier) +
                  " to " + getEnumerationName(dest_tier), do_test);
  Xoroshiro128pGenerator xrs(260148296);
  std::vector<double> umat_hold(9);
  for (int i = 0; i < 9; i++) {
    umat_hold[i] = xrs.gaussianRandomNumber();
    poly_psw.umat[pps_xfrm_offset + i] = umat_hold[i];
  }
  active_pps.upload();
  coordCopy(&active_pps, psys_idx, TrajectoryKind::VELOCITIES, active_cdns, csys_idx, dest_tier,
            orig_tier, gpu);
  criticalCoordTransfer(&active_pps, &active_cdns, dest_tier, orig_tier);  
  diffCoordinates(&poly_psw.xvel[pps_atom_offset], &poly_psw.yvel[pps_atom_offset],
                  &poly_psw.zvel[pps_atom_offset], &poly_psw.xvel_ovrf[pps_atom_offset],
                  &poly_psw.yvel_ovrf[pps_atom_offset], &poly_psw.zvel_ovrf[pps_atom_offset],
                  poly_psw.vel_scale, &cdnsw.xcrd_sp[cdns_atom_offset],
                  &cdnsw.ycrd_sp[cdns_atom_offset], &cdnsw.zcrd_sp[cdns_atom_offset], nullptr,
                  nullptr, nullptr, 1.0, natom, &poly_psw.umat[pps_xfrm_offset], umat_hold.data(),
                  "a PhaseSpaceSynthesis system into a Condensate, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);

  // Reload and transfer information from the PhaseSpaceSynthesis to the Condensate.
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    coordCopy(&active_pps, i, poly_ps, i);
  }
  for (int i = 0; i < cdns.getSystemCount(); i++) {
    coordCopy(&active_cdns, i, cdns, i);    
  }
  active_pps.upload();
  active_cdns.upload();
  coordCopy(&active_cdns, csys_idx, active_pps, psys_idx, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_cdns, &active_pps, dest_tier, orig_tier);  
  diffCoordinates(&cdnsw.xcrd_sp[cdns_atom_offset], &cdnsw.ycrd_sp[cdns_atom_offset],
                  &cdnsw.zcrd_sp[cdns_atom_offset], nullptr, nullptr, nullptr, 1.0,
                  &poly_psw.xcrd[pps_atom_offset], &poly_psw.ycrd[pps_atom_offset],
                  &poly_psw.zcrd[pps_atom_offset], &poly_psw.xcrd_ovrf[pps_atom_offset],
                  &poly_psw.ycrd_ovrf[pps_atom_offset], &poly_psw.zcrd_ovrf[pps_atom_offset],
                  poly_psw.gpos_scale, natom, &cdnsw.umat[cdns_xfrm_offset],
                  &poly_psw.umat[pps_xfrm_offset], "a PhaseSpaceSynthesis system into a "
                  "Condensate, " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test, 7.5e-7, 1.0e-6);
  for (int i = 0; i < 9; i++) {
    umat_hold[i] = xrs.gaussianRandomNumber();
    cdnsw.umat[cdns_xfrm_offset + i] = umat_hold[i];
  }
  active_cdns.upload();
  coordCopy(&active_cdns, csys_idx, active_pps, psys_idx, TrajectoryKind::FORCES, dest_tier,
            orig_tier, gpu);
  criticalCoordTransfer(&active_cdns, &active_pps, dest_tier, orig_tier);  
  diffCoordinates(&cdnsw.xcrd_sp[cdns_atom_offset], &cdnsw.ycrd_sp[cdns_atom_offset],
                  &cdnsw.zcrd_sp[cdns_atom_offset], nullptr, nullptr, nullptr, 1.0,
                  &poly_psw.xfrc[pps_atom_offset], &poly_psw.yfrc[pps_atom_offset],
                  &poly_psw.zfrc[pps_atom_offset], &poly_psw.xfrc_ovrf[pps_atom_offset],
                  &poly_psw.yfrc_ovrf[pps_atom_offset], &poly_psw.zfrc_ovrf[pps_atom_offset],
                  poly_psw.frc_scale, natom, &cdnsw.umat[cdns_xfrm_offset], umat_hold.data(),
                  "a PhaseSpaceSynthesis (primary forces) system into a Condensate, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test,
                  7.5e-7, 1.0e-6);
}

void copySetup(const PhaseSpaceSynthesis &dest, const int dest_idx,
               const PhaseSpaceSynthesis &orig, const int orig_idx,
               const HybridTargetLevel dest_tier, const HybridTargetLevel orig_tier,
               const GpuDetails &gpu, const TestPriority do_test) {

  // Duplicate the inputs and copy the PhaseSpaceSynthesis positions into the CoordinateSeries at
  // a selected frame.
  PhaseSpaceSynthesis active_dest = dest;
  PhaseSpaceSynthesis active_orig = orig;
  PsSynthesisWriter destw = active_dest.data();
  PsSynthesisWriter origw = active_orig.data();
  coordCopy(&active_dest, dest_idx, active_orig, orig_idx, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_dest, &active_orig, dest_tier, orig_tier);
  const int natom = orig.getAtomCount(orig_idx);
  const int dest_atom_start = dest.getAtomOffset(dest_idx);
  const int orig_atom_start = orig.getAtomOffset(orig_idx);
  const int padded_xfrm = roundUp(9, warp_size_int);
  const int dest_xfrm_start = padded_xfrm * dest_idx;
  const int orig_xfrm_start = padded_xfrm * orig_idx;
  diffCoordinates(&destw.xcrd[dest_atom_start], &destw.ycrd[dest_atom_start],
                  &destw.zcrd[dest_atom_start], &destw.xcrd_ovrf[dest_atom_start],
                  &destw.ycrd_ovrf[dest_atom_start], &destw.zcrd_ovrf[dest_atom_start],
                  destw.gpos_scale, &origw.xcrd[orig_atom_start], &origw.ycrd[orig_atom_start],
                  &origw.zcrd[orig_atom_start], &origw.xcrd_ovrf[orig_atom_start],
                  &origw.ycrd_ovrf[orig_atom_start], &origw.zcrd_ovrf[orig_atom_start],
                  origw.gpos_scale, natom, &destw.umat[dest_xfrm_start],
                  &origw.umat[orig_xfrm_start], "a system in one PhaseSpaceSynthesis to "
                  "another, " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);
  diffCoordinates(&destw.xvel[dest_atom_start], &destw.yvel[dest_atom_start],
                  &destw.zvel[dest_atom_start], &destw.xvel_ovrf[dest_atom_start],
                  &destw.yvel_ovrf[dest_atom_start], &destw.zvel_ovrf[dest_atom_start],
                  destw.vel_scale, &origw.xvel[orig_atom_start], &origw.yvel[orig_atom_start],
                  &origw.zvel[orig_atom_start], &origw.xvel_ovrf[orig_atom_start],
                  &origw.yvel_ovrf[orig_atom_start], &origw.zvel_ovrf[orig_atom_start],
                  origw.vel_scale, natom, &destw.umat[dest_xfrm_start],
                  &origw.umat[orig_xfrm_start], "a system in one PhaseSpaceSynthesis to "
                  "another (velocities), " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);
  diffCoordinates(&destw.xfrc[dest_atom_start], &destw.yfrc[dest_atom_start],
                  &destw.zfrc[dest_atom_start], &destw.xfrc_ovrf[dest_atom_start],
                  &destw.yfrc_ovrf[dest_atom_start], &destw.zfrc_ovrf[dest_atom_start],
                  destw.frc_scale, &origw.xfrc[orig_atom_start], &origw.yfrc[orig_atom_start],
                  &origw.zfrc[orig_atom_start], &origw.xfrc_ovrf[orig_atom_start],
                  &origw.yfrc_ovrf[orig_atom_start], &origw.zfrc_ovrf[orig_atom_start],
                  origw.frc_scale, natom, &destw.umat[dest_xfrm_start],
                  &origw.umat[orig_xfrm_start], "a system in one PhaseSpaceSynthesis to "
                  "another (forces), " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);
  diffCoordinates(&destw.xalt[dest_atom_start], &destw.yalt[dest_atom_start],
                  &destw.zalt[dest_atom_start], &destw.xalt_ovrf[dest_atom_start],
                  &destw.yalt_ovrf[dest_atom_start], &destw.zalt_ovrf[dest_atom_start],
                  destw.gpos_scale, &origw.xalt[orig_atom_start], &origw.yalt[orig_atom_start],
                  &origw.zalt[orig_atom_start], &origw.xalt_ovrf[orig_atom_start],
                  &origw.yalt_ovrf[orig_atom_start], &origw.zalt_ovrf[orig_atom_start],
                  origw.gpos_scale, natom, &destw.umat_alt[dest_xfrm_start],
                  &origw.umat_alt[orig_xfrm_start], "a system in one PhaseSpaceSynthesis to "
                  "another (alt. positions), " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);
  diffCoordinates(&destw.vxalt[dest_atom_start], &destw.vyalt[dest_atom_start],
                  &destw.vzalt[dest_atom_start], &destw.vxalt_ovrf[dest_atom_start],
                  &destw.vyalt_ovrf[dest_atom_start], &destw.vzalt_ovrf[dest_atom_start],
                  destw.vel_scale, &origw.vxalt[orig_atom_start], &origw.vyalt[orig_atom_start],
                  &origw.vzalt[orig_atom_start], &origw.vxalt_ovrf[orig_atom_start],
                  &origw.vyalt_ovrf[orig_atom_start], &origw.vzalt_ovrf[orig_atom_start],
                  origw.vel_scale, natom, &destw.umat_alt[dest_xfrm_start],
                  &origw.umat_alt[orig_xfrm_start], "a system in one PhaseSpaceSynthesis to "
                  "another (alt. velocities), " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);
  diffCoordinates(&destw.fxalt[dest_atom_start], &destw.fyalt[dest_atom_start],
                  &destw.fzalt[dest_atom_start], &destw.fxalt_ovrf[dest_atom_start],
                  &destw.fyalt_ovrf[dest_atom_start], &destw.fzalt_ovrf[dest_atom_start],
                  destw.frc_scale, &origw.fxalt[orig_atom_start], &origw.fyalt[orig_atom_start],
                  &origw.fzalt[orig_atom_start], &origw.fxalt_ovrf[orig_atom_start],
                  &origw.fyalt_ovrf[orig_atom_start], &origw.fzalt_ovrf[orig_atom_start],
                  origw.frc_scale, natom, &destw.umat_alt[dest_xfrm_start],
                  &origw.umat_alt[orig_xfrm_start], "a system in one PhaseSpaceSynthesis to "
                  "another (alt. forces), " + getEnumerationName(orig_tier) + " to " +
                  getEnumerationName(dest_tier), do_test);
}

void copySetup(const Condensate &dest, const int dest_idx, const Condensate &orig,
               const int orig_idx, const HybridTargetLevel dest_tier,
               const HybridTargetLevel orig_tier, const GpuDetails &gpu,
               const TestPriority do_test) {

  // Duplicate the inputs and copy one to the other.
  Condensate active_dest = dest;
  Condensate active_orig = orig;
  CondensateWriter destw = active_dest.data();
  CondensateWriter origw = active_orig.data();
  coordCopy(&active_dest, dest_idx, active_orig, orig_idx, dest_tier, orig_tier, gpu);
  criticalCoordTransfer(&active_dest, &active_orig, dest_tier, orig_tier);
  const int dest_atom_offset = dest.getAtomOffset(dest_idx);
  const int orig_atom_offset = orig.getAtomOffset(orig_idx);
  const int padded_xfrm = roundUp(9, warp_size_int);
  const int dest_xfrm_offset = padded_xfrm * dest_idx;
  const int orig_xfrm_offset = padded_xfrm * orig_idx;
  diffCoordinates(&destw.xcrd_sp[dest_atom_offset], &destw.ycrd_sp[dest_atom_offset],
                  &destw.zcrd_sp[dest_atom_offset], 1.0, &origw.xcrd_sp[orig_atom_offset],
                  &origw.ycrd_sp[orig_atom_offset], &origw.zcrd_sp[orig_atom_offset], 1.0,
                  orig.getAtomCount(orig_idx), &destw.umat[dest_xfrm_offset],
                  &origw.umat[orig_xfrm_offset],"one system of a Condensate into another, " +
                  getEnumerationName(orig_tier) + " to " + getEnumerationName(dest_tier), do_test);
}
#endif

//-------------------------------------------------------------------------------------------------
// Add noise to a PhaseSpace object, all aspects in both coordinate cycles.  The abstract must be
// valid on the CPU host.
//
// Arguments:
//   xrs:  Source of random numbers
//   psw:  Abstract of the coordinate object
//-------------------------------------------------------------------------------------------------
void perturbPhaseSpace(Xoroshiro128pGenerator *xrs, PhaseSpaceWriter *psw) {
  addRandomNoise(xrs, psw->xcrd, psw->ycrd, psw->zcrd, psw->natom, 0.1);
  addRandomNoise(xrs, psw->xvel, psw->yvel, psw->zfrc, psw->natom, 0.1);
  addRandomNoise(xrs, psw->xfrc, psw->yfrc, psw->zvel, psw->natom, 0.1);
  addRandomNoise(xrs, psw->xalt, psw->yalt, psw->zalt, psw->natom, 0.1);
  addRandomNoise(xrs, psw->vxalt, psw->vyalt, psw->vzalt, psw->natom, 0.1);
  addRandomNoise(xrs, psw->fxalt, psw->fyalt, psw->fzalt, psw->natom, 0.1);
}

//-------------------------------------------------------------------------------------------------
// Add noise to one system of a PhaseSpaceSynthesis based on the abstract taken in one orientation.
// The abstract must be valid on the CPU host.
//
// Arguments:
//   xrs:       Source of random numbers
//   poly_psw:  Abstract of the coordinate synthesis
//   istart:    Staring index of atoms to perturb
//   natom:     The number of atoms in the system of interest
//-------------------------------------------------------------------------------------------------
void perturbSynthesis(Xoroshiro128pGenerator *xrs, PsSynthesisWriter *poly_psw, const int istart,
                      const int natom) {
  addRandomNoise(xrs, &poly_psw->xcrd[istart], &poly_psw->xcrd_ovrf[istart],
                 &poly_psw->ycrd[istart], &poly_psw->ycrd_ovrf[istart], &poly_psw->zcrd[istart],
                 &poly_psw->zcrd_ovrf[istart], natom, 0.1, poly_psw->gpos_scale);
  addRandomNoise(xrs, &poly_psw->xvel[istart], &poly_psw->xvel_ovrf[istart],
                 &poly_psw->yvel[istart], &poly_psw->yvel_ovrf[istart], &poly_psw->zvel[istart],
                 &poly_psw->zvel_ovrf[istart], natom, 0.1, poly_psw->vel_scale);
  addRandomNoise(xrs, &poly_psw->xfrc[istart], &poly_psw->xfrc_ovrf[istart],
                 &poly_psw->yfrc[istart], &poly_psw->yfrc_ovrf[istart], &poly_psw->zfrc[istart],
                 &poly_psw->zfrc_ovrf[istart], natom, 0.1, poly_psw->frc_scale);
  addRandomNoise(xrs, &poly_psw->xalt[istart], &poly_psw->xalt_ovrf[istart],
                 &poly_psw->yalt[istart], &poly_psw->yalt_ovrf[istart], &poly_psw->zalt[istart],
                 &poly_psw->zalt_ovrf[istart], natom, 0.1, poly_psw->gpos_scale);
  addRandomNoise(xrs, &poly_psw->vxalt[istart], &poly_psw->vxalt_ovrf[istart],
                 &poly_psw->vyalt[istart], &poly_psw->vyalt_ovrf[istart], &poly_psw->vzalt[istart],
                 &poly_psw->vzalt_ovrf[istart], natom, 0.1, poly_psw->vel_scale);
  addRandomNoise(xrs, &poly_psw->fxalt[istart], &poly_psw->fxalt_ovrf[istart],
                 &poly_psw->fyalt[istart], &poly_psw->fyalt_ovrf[istart], &poly_psw->fzalt[istart],
                 &poly_psw->fzalt_ovrf[istart], natom, 0.1, poly_psw->frc_scale);
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Baseline initialization
  const TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1: test the veracity of the Condensate object
  section("Verify Condensate contents");

  // Section 2: test RMSD calculations with the Condensate
  section("Test RMSD calcuations");
  
  // Read all systems
  const char osc = osSeparator();
  const std::string base_top_dir = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_dir = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::vector<std::string> system_names = { "drug_example", "drug_example_vs",
                                                  "bromobenzene", "bromobenzene_vs", "biotin" };
  TestSystemManager tsm(base_top_dir, "top", system_names, base_crd_dir, "inpcrd", system_names);

  // Form a synthesis of the coordinates.  Perturb each structure slightly in the synthesis.
  // Simultaneously construct a systems cache to map the synthesis back to its origins.
  Xoroshiro128pGenerator xrs(517389207);
  std::vector<AtomGraph*> pbc_ag_list;
  std::vector<PhaseSpace> pbc_ps_list;
  const int n_structures = 100;
  std::vector<int> index_memory(n_structures);
  pbc_ag_list.reserve(n_structures);
  pbc_ps_list.reserve(n_structures);
  int pbc_con = 0;
  while (pbc_con < n_structures) {
    const int ag_id = static_cast<double>(tsm.getSystemCount()) * xrs.uniformRandomNumber();
    if (tsm.getTopologyPointer(ag_id)->getUnitCellType() == UnitCellType::NONE) {
      continue;
    }
    index_memory[pbc_con] = ag_id;
    pbc_ag_list.push_back(tsm.getTopologyPointer(ag_id));
    pbc_ps_list.push_back(tsm.exportPhaseSpace(ag_id));
    PhaseSpaceWriter psw = pbc_ps_list.back().data();
    for (int i = 0; i < psw.natom; i++) {
      psw.xcrd[i] += 0.5 * (xrs.uniformRandomNumber() - 0.5);
      psw.ycrd[i] += 0.5 * (xrs.uniformRandomNumber() - 0.5);
      psw.zcrd[i] += 0.5 * (xrs.uniformRandomNumber() - 0.5);
    }
    pbc_con++;
  }
  PhaseSpaceSynthesis poly_ps(pbc_ps_list, pbc_ag_list, 48);
  std::string fcon_string = "&files\n";
  for (int i = 0; i < 3; i++) {
    std::string ilabel;
    if (i == 0) {
      ilabel = "first_";
    }
    else if (i == 1) {
      ilabel = "second_";
    }
    else {
      ilabel = "third_";
    }
    ilabel += "group";
    for (int j = 0; j < tsm.getSystemCount(); j++) {
      fcon_string += "  -sys { -p " + tsm.getTopologyFile(j) + " -c " + tsm.getCoordinateFile(j) +
                     " -label " + ilabel + " }\n";
    }
  }
  fcon_string += "&end\n";
  const TextFile fcon_tf(fcon_string, TextOrigin::RAM);
  int start_line = 0;
  const FilesControls fcon(fcon_tf, &start_line);
  const SystemCache sc(fcon);
  for (int i = 0; i < n_structures; i++) {
    index_memory[i] += tsm.getSystemCount() * 2.0 * round(xrs.uniformRandomNumber());
  }
  const SynthesisCacheMap scmap(index_memory, sc, poly_ps);
  
  // Create a condensate from the synthesis.
  Condensate cdns_dbl(poly_ps, PrecisionModel::DOUBLE);
  Condensate cdns_flt(poly_ps, PrecisionModel::SINGLE);
  section(1);
  const int n_test_pt = 100;
  std::vector<double> xdbl_dev(n_test_pt), ydbl_dev(n_test_pt), zdbl_dev(n_test_pt);
  std::vector<double> xflt_dev(n_test_pt), yflt_dev(n_test_pt), zflt_dev(n_test_pt);
  CondensateWriter cdw_dbl = cdns_dbl.data();
  CondensateWriter cdw_flt = cdns_flt.data();
  PsSynthesisWriter poly_psw = poly_ps.data();
  const ComparisonGuide cg_dbl(cdns_dbl, scmap);
  const CompGuideKit cgr_dbl = cg_dbl.data();
  for (int i = 0; i < n_test_pt; i++) {
    const int sys_idx = poly_psw.system_count * xrs.uniformRandomNumber();
    const int atom_idx = static_cast<int>(static_cast<double>(poly_psw.atom_counts[sys_idx]) *
                                          xrs.uniformRandomNumber()) +
                         poly_psw.atom_starts[sys_idx];
    xdbl_dev[i] = (hostInt95ToDouble(poly_psw.xcrd[atom_idx], poly_psw.xcrd_ovrf[atom_idx]) *
                   poly_psw.inv_gpos_scale) - cdw_dbl.xcrd[atom_idx];
    ydbl_dev[i] = (hostInt95ToDouble(poly_psw.ycrd[atom_idx], poly_psw.ycrd_ovrf[atom_idx]) *
                   poly_psw.inv_gpos_scale) - cdw_dbl.ycrd[atom_idx];
    zdbl_dev[i] = (hostInt95ToDouble(poly_psw.zcrd[atom_idx], poly_psw.zcrd_ovrf[atom_idx]) *
                   poly_psw.inv_gpos_scale) - cdw_dbl.zcrd[atom_idx];
    xflt_dev[i] = (hostInt95ToDouble(poly_psw.xcrd[atom_idx], poly_psw.xcrd_ovrf[atom_idx]) *
                   poly_psw.inv_gpos_scale) - cdw_flt.xcrd_sp[atom_idx];
    yflt_dev[i] = (hostInt95ToDouble(poly_psw.ycrd[atom_idx], poly_psw.ycrd_ovrf[atom_idx]) *
                   poly_psw.inv_gpos_scale) - cdw_flt.ycrd_sp[atom_idx];
    zflt_dev[i] = (hostInt95ToDouble(poly_psw.zcrd[atom_idx], poly_psw.zcrd_ovrf[atom_idx]) *
                   poly_psw.inv_gpos_scale) - cdw_flt.zcrd_sp[atom_idx];
  }
  check(xdbl_dev, RelationalOperator::EQUAL, std::vector<double>(n_test_pt, 0.0), "Cartesian X "
        "coordinates kept in the Condensate differ from the original PhaseSpaceSynthesis, even "
        "when recorded in double-precision.", tsm.getTestingStatus());
  check(ydbl_dev, RelationalOperator::EQUAL, std::vector<double>(n_test_pt, 0.0), "Cartesian Y "
        "coordinates kept in the Condensate differ from the original PhaseSpaceSynthesis, even "
        "when recorded in double-precision.", tsm.getTestingStatus());
  check(zdbl_dev, RelationalOperator::EQUAL, std::vector<double>(n_test_pt, 0.0), "Cartesian Z "
        "coordinates kept in the Condensate differ from the original PhaseSpaceSynthesis, even "
        "when recorded in double-precision.", tsm.getTestingStatus());
  check(xflt_dev, RelationalOperator::EQUAL,
        Approx(std::vector<double>(n_test_pt, 0.0)).margin(6.0e-6), "Cartesian X coordinates kept "
        "in the Condensate differ from the original PhaseSpaceSynthesis, even when recorded in "
        "double-precision.", tsm.getTestingStatus());
  check(yflt_dev, RelationalOperator::EQUAL,
        Approx(std::vector<double>(n_test_pt, 0.0)).margin(6.0e-6), "Cartesian Y coordinates kept "
        "in the Condensate differ from the original PhaseSpaceSynthesis, even when recorded in "
        "double-precision.", tsm.getTestingStatus());
  check(zflt_dev, RelationalOperator::EQUAL,
        Approx(std::vector<double>(n_test_pt, 0.0)).margin(6.0e-6), "Cartesian Z coordinates kept "
        "in the Condensate differ from the original PhaseSpaceSynthesis, even when recorded in "
        "double-precision.", tsm.getTestingStatus());
  check(cgr_dbl.natr_insr_top, RelationalOperator::EQUAL, 27, "The number of all-to-reference "
        "instructions for systems grouped by topology does not agree with the number of systems "
        "for a CPU-bound Condensate.", tsm.getTestingStatus());
  check(cgr_dbl.nata_insr_top, RelationalOperator::EQUAL, 27, "The number of all-to-all "
        "instructions for systems grouped by topology does not agree with the number of systems "
        "for a CPU-bound Condensate.", tsm.getTestingStatus());

  // Check that all-to-one instructions for various system groupings cover all systems.
  std::vector<int> atr_occupancy_src(cdw_dbl.system_count, 0);
  std::vector<int> atr_occupancy_top(cdw_dbl.system_count, 0);
  std::vector<int> atr_occupancy_lbl(cdw_dbl.system_count, 0);
  const SynthesisMapReader scmapr = scmap.data();
  const int uint_bits = sizeof(uint) * 8;
  for (int i = 0; i < cgr_dbl.natr_insr_src; i++) {
    const int4 tinsr = cgr_dbl.atr_member_src[i];
    const int tgroup = cgr_dbl.atr_group_src[i];
    for (int j = 0; j < tinsr.w; j++) {
      atr_occupancy_src[scmapr.csystem_proj[scmapr.csystem_bounds[tgroup] + tinsr.x + j]] = 1;
    }
  }
  for (int i = 0; i < cgr_dbl.natr_insr_top; i++) {
    const int4 tinsr = cgr_dbl.atr_member_top[i];
    const int tgroup = cgr_dbl.atr_group_top[i];
    for (int j = 0; j < tinsr.w; j++) {
      atr_occupancy_top[scmapr.ctopol_proj[scmapr.ctopol_bounds[tgroup] + tinsr.x + j]] = 1;
    }
  }
  for (int i = 0; i < cgr_dbl.natr_insr_lbl; i++) {
    const int4 tinsr = cgr_dbl.atr_member_lbl[i];
    const int tgroup = cgr_dbl.atr_group_lbl[i];
    for (int j = 0; j < tinsr.w; j++) {
      atr_occupancy_lbl[scmapr.clabel_proj[scmapr.clabel_bounds[tgroup] + tinsr.x + j]] = 1;
    }
  }
  check(atr_occupancy_src, RelationalOperator::EQUAL, std::vector<int>(cdw_dbl.system_count, 1),
        "Not all systems are covered in all-to-reference calculations for source groups.");
  check(atr_occupancy_top, RelationalOperator::EQUAL, std::vector<int>(cdw_dbl.system_count, 1),
        "Not all systems are covered in all-to-reference calculations for topology groups.");
  check(atr_occupancy_lbl, RelationalOperator::EQUAL, std::vector<int>(cdw_dbl.system_count, 1),
        "Not all systems are covered in all-to-reference calculations for label groups.");
  
  // Check that all-to-all instructions for various sytem groupings match the coverage of a
  // direct listing.
  checkATAFootprint(cg_dbl, cdns_dbl, scmap, SystemGrouping::SOURCE);
  checkATAFootprint(cg_dbl, cdns_dbl, scmap, SystemGrouping::TOPOLOGY);
  checkATAFootprint(cg_dbl, cdns_dbl, scmap, SystemGrouping::LABEL);
  
  // Try an RMSD calculation with the Condensate guiding the process on the CPU.  Compare it to
  // an RMSD calculation using the synthesis directly.
  RMSDPlan rplan(poly_ps, sc, scmap);
  const ComparisonGuide synth_cg(poly_ps, scmap);
  const SystemGrouping top_grp = SystemGrouping::TOPOLOGY;
  Hybrid<double> atr_rmsd_a(synth_cg.getAllToReferenceOutputSize());
  Hybrid<double> atr_rmsd_b(synth_cg.getAllToReferenceOutputSize());
  Hybrid<double> atr_rmsd_c(synth_cg.getAllToReferenceOutputSize());
  Hybrid<double> ata_rmsd_a(synth_cg.getSymmetryEquivalentPairOutputSize(top_grp));
  Hybrid<double> ata_rmsd_b(synth_cg.getSymmetryEquivalentPairOutputSize(top_grp));
  Hybrid<double> ata_rmsd_c(synth_cg.getSymmetryEquivalentPairOutputSize(top_grp));
  Hybrid<int> examples(sc.getTopologyCount(), "top_grouping");
  std::vector<int> example_indices = poly_ps.getUniqueTopologyExampleIndices();
  for (int i = 0; i < poly_ps.getUniqueTopologyCount(); i++) {
    examples.putHost(example_indices[i], scmap.getTopologyCacheIndex(example_indices[i]));
  }
  rmsd(cg_dbl, rplan, poly_ps, examples, &atr_rmsd_a);
  rmsd(cg_dbl, rplan, poly_ps, &ata_rmsd_a);
  rmsd(cg_dbl, rplan, poly_ps, cdns_dbl, examples, &atr_rmsd_b);
  rmsd(cg_dbl, rplan, poly_ps, cdns_dbl, &ata_rmsd_b);
  rmsd(cg_dbl, rplan, poly_ps, cdns_flt, examples, &atr_rmsd_c);
  rmsd(cg_dbl, rplan, poly_ps, cdns_flt, &ata_rmsd_c);
  const std::vector<double> atr_rmsd_av = atr_rmsd_a.readHost();
  const std::vector<double> atr_rmsd_bv = atr_rmsd_b.readHost();
  const std::vector<double> atr_rmsd_cv = atr_rmsd_c.readHost();
  const std::vector<double> ata_rmsd_av = ata_rmsd_a.readHost();
  const std::vector<double> ata_rmsd_bv = ata_rmsd_b.readHost();
  const std::vector<double> ata_rmsd_cv = ata_rmsd_c.readHost();
  section(2);
  check(atr_rmsd_bv, RelationalOperator::EQUAL, atr_rmsd_av, "RMSD (all to reference) "
        "calculations guided by a Condensate object's instructions on the CPU do not match RMSD "
        "calculations performed on the original PhaseSpaceSynthesis.", tsm.getTestingStatus());
  check(ata_rmsd_bv, RelationalOperator::EQUAL, ata_rmsd_av, "RMSD (all to all) calculations "
        "guided by a Condensate object's instructions on the CPU do not match RMSD calculations "
        "performed on the original PhaseSpaceSynthesis.", tsm.getTestingStatus());
  check(ata_rmsd_cv, RelationalOperator::EQUAL, Approx(ata_rmsd_av).margin(1.5e-5), "RMSD (all to "
        "all) calculations guided by a Condensate object's instructions on the CPU do not match "
        "RMSD calculations performed on the original PhaseSpaceSynthesis when the Condensate is "
        "cast in 32 bit floating-point precision.", tsm.getTestingStatus());

  // Try additional RMSD calculations based on the source.
  const SystemGrouping src_grp = SystemGrouping::SOURCE;
  Hybrid<double> atr_rmsd_src_a(cg_dbl.getAllToReferenceOutputSize());
  Hybrid<double> atr_rmsd_src_b(cg_dbl.getAllToReferenceOutputSize());
  Hybrid<double> ata_rmsd_src_a(cg_dbl.getSymmetryEquivalentPairOutputSize(src_grp));
  Hybrid<double> ata_rmsd_src_b(cg_dbl.getSymmetryEquivalentPairOutputSize(src_grp));
  Hybrid<int> src_examples(sc.getSystemCount());
  std::vector<double> atr_rmsd_src_ans(cg_dbl.getAllToReferenceOutputSize());
  int vresult_idx = 0;
  std::vector<std::vector<int>> system_examples(sc.getSystemCount(), std::vector<int>());
  for (int i = 0 ; i < sc.getSystemCount(); i++) {

    // Construct the all-to-reference RMSD results and construct a list of all systems spawned by
    // each source.
    bool iseek = true;
    for (int j = 0; j < scmapr.nsynth; j++) {
      if (iseek && scmapr.cache_origins[j] == i) {
        src_examples.putHost(j, i);
        iseek = false;
        const CoordinateFrame cf_ref = poly_ps.exportCoordinates(j);
        const CoordinateFrameReader cf_refr = cf_ref.data();
        const ChemicalDetailsKit cdk =
          poly_ps.getSystemTopologyPointer(j)->getChemicalDetailsKit();
        for (int k = j; k < scmapr.nsynth; k++) {
          if (scmapr.cache_origins[k] == i) {
            const CoordinateFrame cfk = poly_ps.exportCoordinates(k);
            atr_rmsd_src_ans[vresult_idx] = rmsd(cf_refr, cfk.data(), cdk, RMSDMethod::ALIGN_MASS);
            vresult_idx++;
            system_examples[i].push_back(k);
          }
        }
      }
    }
  }

  // Construct the all-to-all RMSD results
  std::vector<double> ata_rmsd_src_ans(cg_dbl.getSymmetryEquivalentPairOutputSize(src_grp));
  int isq_offset = 0;
  for (int i = 0 ; i < sc.getSystemCount(); i++) {
    const int ni = system_examples[i].size();
    if (ni == 0) {
      continue;
    }
    const AtomGraph *iag = poly_ps.getSystemTopologyPointer(system_examples[i][0]);
    const ChemicalDetailsKit icdk = iag->getChemicalDetailsKit();
    for (int j = 0; j < ni; j++) {
      const CoordinateFrame cf_j = poly_ps.exportCoordinates(system_examples[i][j]);
      const CoordinateFrameReader cf_jr = cf_j.data();
      size_t jkr_idx = isq_offset + (j * (j - 1) / 2);
      for (int k = 0; k < j; k++) {
        const CoordinateFrame cf_k = poly_ps.exportCoordinates(system_examples[i][k]);
        const CoordinateFrameReader cf_kr = cf_k.data();
        ata_rmsd_src_ans[jkr_idx + k] = rmsd(cf_jr, cf_kr, icdk, RMSDMethod::ALIGN_MASS);
      }
    }
    isq_offset += roundUp(ni * (ni - 1) / 2, warp_size_int);
  }
  rmsd(synth_cg, rplan, poly_ps, src_examples, &atr_rmsd_src_a, SystemGrouping::SOURCE);
  rmsd(synth_cg, rplan, poly_ps, &ata_rmsd_src_a, SystemGrouping::SOURCE);
  rmsd(cg_dbl, rplan, poly_ps, cdns_dbl, src_examples, &atr_rmsd_src_b, SystemGrouping::SOURCE);
  rmsd(cg_dbl, rplan, poly_ps, cdns_dbl, &ata_rmsd_src_b, SystemGrouping::SOURCE);
  const std::vector<double> atr_rmsd_src_va = atr_rmsd_src_a.readHost();
  const std::vector<double> ata_rmsd_src_va = ata_rmsd_src_a.readHost();
  const std::vector<double> atr_rmsd_src_vb = atr_rmsd_src_a.readHost();
  const std::vector<double> ata_rmsd_src_vb = ata_rmsd_src_a.readHost();
  check(atr_rmsd_src_va, RelationalOperator::EQUAL, atr_rmsd_src_ans, "RMSD (all to reference) "
        "calculations guide by a plan did not meet expectations when grouping systems by cache "
        "system (source).", tsm.getTestingStatus());
  check(ata_rmsd_src_va, RelationalOperator::EQUAL, ata_rmsd_src_ans, "RMSD (all to all) "
        "calculations guide by a plan did not meet expectations when grouping systems by cache "
        "system (source).", tsm.getTestingStatus());
  check(atr_rmsd_src_vb, RelationalOperator::EQUAL, atr_rmsd_src_ans, "RMSD (all to reference) "
        "calculations guide by a plan did not meet expectations when grouping systems by cache "
        "system (source) and using a Condensate object with work instructions.",
        tsm.getTestingStatus());
  check(ata_rmsd_src_vb, RelationalOperator::EQUAL, ata_rmsd_src_ans, "RMSD (all to all) "
        "calculations guide by a plan did not meet expectations when grouping systems by cache "
        "system (source) and using a Condensate object with work instructions.",
        tsm.getTestingStatus());
  
  // With the most complex coordinate objects now established, check various coordinate copying
  // overloads.
  std::vector<CoordinateFrame> tsm_cf;
  std::vector<PhaseSpace> tsm_ps;
  std::vector<CoordinateSeries<llint>> tsm_cs;
  std::vector<CoordinateFrameReader> tsm_cfr;
  std::vector<PhaseSpaceReader> tsm_psr;
  std::vector<CoordinateSeriesReader<llint>> tsm_csr;
  PhaseSpaceSynthesis tsm_psyn = spawnMutableCoordinateObjects(tsm, &tsm_cf, &tsm_cfr, &tsm_ps,
                                                               &tsm_psr, &tsm_cs, &tsm_csr);
  Condensate tsm_cdns(tsm_psyn, PrecisionModel::SINGLE);

  // Test coordinate copying into a PhaseSpace object
  PhaseSpace recv(tsm_cf[0].getAtomCount());
  const PhaseSpaceReader recv_r = recv.data();
  coordCopy(&recv, tsm_cf[0]);
  coordCopy(&recv, TrajectoryKind::VELOCITIES, tsm_cs[0].exportFrame(1));
  coordCopy(&recv, TrajectoryKind::FORCES, tsm_cs[0].exportFrame(2));
  coordCopy(&recv, TrajectoryKind::POSITIONS, CoordinateCycle::BLACK, tsm_cf[0]);
  coordCopy(&recv, TrajectoryKind::FORCES, CoordinateCycle::BLACK, tsm_cdns,
            tsm.getSystemCount());

  // The X, Y, and Z coordinates of the PhaseSpace object recv's WHITE coordinates should match
  // those of the test system manager's first coordinate frame, exactly.
  diffCoordinates<double, double>(recv_r.xcrd, recv_r.ycrd, recv_r.zcrd, tsm_cfr[0].xcrd,
                                  tsm_cfr[0].ycrd, tsm_cfr[0].zcrd, recv_r.natom,
                                  "CoordinateFrame => PhaseSpace(POSITIONS, WHITE)",
                                  tsm.getTestingStatus());

  // The PhaseSpace object recv's WHITE velocities and forces should match those of the second
  // and third frames of the series created for the same system, to a precision of about one part
  // in sixteen quadrillion (just over 64-bit floating point precision).
  const int second_frame_offset =     roundUp(tsm_csr[0].natom, warp_size_int);
  const int third_frame_offset  = 2 * roundUp(tsm_csr[0].natom, warp_size_int);
  diffCoordinates<double, llint>(recv_r.xvel, recv_r.yvel, recv_r.zvel, 1.0,
                                 &tsm_csr[0].xcrd[second_frame_offset],
                                 &tsm_csr[0].ycrd[second_frame_offset],
                                 &tsm_csr[0].zcrd[second_frame_offset], tsm_csr[0].gpos_scale,
                                 recv_r.natom,
                                 "CoordinateSeries(1) => PhaseSpace(VELOCITIES, WHITE)",
                                 tsm.getTestingStatus());
  diffCoordinates<double, llint>(recv_r.xfrc, recv_r.yfrc, recv_r.zfrc, 1.0,
                                 &tsm_csr[0].xcrd[third_frame_offset],
                                 &tsm_csr[0].ycrd[third_frame_offset],
                                 &tsm_csr[0].zcrd[third_frame_offset], tsm_csr[0].gpos_scale,
                                 recv_r.natom,
                                 "CoordinateSeries(2) => PhaseSpace(FORCES, WHITE)",
                                 tsm.getTestingStatus());
  const CondensateWriter tsm_cdnsw = tsm_cdns.data();
  const int second_set_offset = tsm_cdnsw.atom_starts[tsm.getSystemCount()];
  diffCoordinates<double, float>(recv_r.fxalt, recv_r.fyalt, recv_r.fzalt,
                                 &tsm_cdnsw.xcrd_sp[second_set_offset],
                                 &tsm_cdnsw.ycrd_sp[second_set_offset],
                                 &tsm_cdnsw.zcrd_sp[second_set_offset], recv_r.natom,
                                 "Condensate(" + std::to_string(tsm_cdnsw.system_count) +
                                 ") => PhaseSpace(FORCES, BLACK)", tsm.getTestingStatus());

  // Bring synthesis objects to the fore, filling out more aspects of the recv PhaseSpace object.
  CoordinateFrame buffer_a_cf(tsm_cf[0].getAtomCount()), buffer_b_cf(tsm_cf[0].getAtomCount());
  CoordinateFrameWriter buffer_b_cfw = buffer_b_cf.data();
  coordCopy(&buffer_a_cf, tsm_psyn, tsm.getSystemCount());
  coordCopy(&tsm_psyn, 0, TrajectoryKind::FORCES, CoordinateCycle::BLACK, buffer_a_cf);
  coordCopy(&buffer_b_cf, tsm_psyn, 0, TrajectoryKind::FORCES, CoordinateCycle::BLACK);
  coordCopy(&recv, TrajectoryKind::FORCES, CoordinateCycle::BLACK, buffer_a_cf);
  diffCoordinates<double, double>(recv_r.fxalt, recv_r.fyalt, recv_r.fzalt, buffer_b_cfw.xcrd,
                                  buffer_b_cfw.ycrd, buffer_b_cfw.zcrd, recv_r.natom,
                                  "CoordinateFrame => PhaseSpace(FORCES, BLACK)",
                                  tsm.getTestingStatus(), 3.0e-7, 3.0e-7);
  CHECK_THROWS(coordCopy(&recv, TrajectoryKind::VELOCITIES, CoordinateCycle::BLACK, tsm_cdns,
                         1), "Systems with mismatching numbers of atoms were submitted to a copy "
               "operation between PhaseSpace and Condensate objects.");
  coordCopy(&recv, TrajectoryKind::VELOCITIES, CoordinateCycle::BLACK, tsm_cdns, 0);
  diffCoordinates<double, double>(recv_r.vxalt, recv_r.vyalt, recv_r.vzalt, tsm_cfr[0].xcrd,
                                  tsm_cfr[0].ycrd, tsm_cfr[0].zcrd, recv_r.natom,
                                  "CoordinateFrame => PhaseSpace(FORCES, BLACK)",
                                  tsm.getTestingStatus(), 3.0e-6, 2.2e-6);
  
  // Push coordinates into a new PhaseSpaceSynthesis with higher precision
  const std::vector<int> replicator = replicateSeries(tsm.getSystemCount(), 2);
  PhaseSpaceSynthesis test_psyn(tsm_ps, replicator, tsm.getTopologyPointer(), replicator, 62, 24,
                                61, 60);
  PsSynthesisWriter test_psynw = test_psyn.data();
  coordCopy(&test_psyn, 1, TrajectoryKind::VELOCITIES, tsm_cf[1]);
  const int sys1_start = test_psynw.atom_starts[1];
  diffCoordinates(&test_psynw.xvel[sys1_start], &test_psynw.yvel[sys1_start],
                  &test_psynw.zvel[sys1_start], &test_psynw.xvel_ovrf[sys1_start],
                  &test_psynw.yvel_ovrf[sys1_start], &test_psynw.zvel_ovrf[sys1_start],
                  test_psynw.vel_scale, tsm_cfr[1].xcrd, tsm_cfr[1].ycrd, tsm_cfr[1].zcrd,
                  nullptr, nullptr, nullptr, 1.0, tsm_cfr[1].natom,
                  "CoordinateFrame => PhaseSpaceSynthesis(VELOCITIES, WHITE)",
                  tsm.getTestingStatus());
  coordCopy(&test_psyn, 1, TrajectoryKind::FORCES, CoordinateCycle::BLACK, tsm_cs[1], 3);
  const int fourth_frame_offset = 3 * roundUp(tsm_csr[1].natom, warp_size_int);
  diffCoordinates(&test_psynw.fxalt[sys1_start], &test_psynw.fyalt[sys1_start],
                  &test_psynw.fzalt[sys1_start], &test_psynw.fxalt_ovrf[sys1_start],
                  &test_psynw.fyalt_ovrf[sys1_start], &test_psynw.fzalt_ovrf[sys1_start],
                  test_psynw.frc_scale, &tsm_csr[1].xcrd[fourth_frame_offset],
                  &tsm_csr[1].ycrd[fourth_frame_offset], &tsm_csr[1].zcrd[fourth_frame_offset],
                  nullptr, nullptr, nullptr, tsm_csr[1].gpos_scale, tsm_cfr[1].natom,
                  "CoordinateFrame => PhaseSpaceSynthesis(FORCES, BLACK)",
                  tsm.getTestingStatus());
  coordCopy(&test_psyn, 2, TrajectoryKind::FORCES, CoordinateCycle::BLACK, tsm_cdns, 2);
  const int sys2_start = test_psynw.atom_starts[2];
  diffCoordinates(&test_psynw.fxalt[sys2_start], &test_psynw.fyalt[sys2_start],
                  &test_psynw.fzalt[sys2_start], &test_psynw.fxalt_ovrf[sys2_start],
                  &test_psynw.fyalt_ovrf[sys2_start], &test_psynw.fzalt_ovrf[sys2_start],
                  test_psynw.frc_scale, &tsm_cdnsw.xcrd_sp[sys2_start],
                  &tsm_cdnsw.ycrd_sp[sys2_start], &tsm_cdnsw.zcrd_sp[sys2_start], nullptr, nullptr,
                  nullptr, 1.0, tsm_cdnsw.atom_counts[2],
                  "Condensate => PhaseSpaceSynthesis(FORCES, BLACK)", tsm.getTestingStatus());
  coordCopy(&test_psyn, tsm.getSystemCount(), recv);
  CoordinateFrameWriter buffer_a_cfw = buffer_a_cf.data();
  diffCoordinates(&test_psynw.fxalt[second_set_offset], &test_psynw.fyalt[second_set_offset],
                  &test_psynw.fzalt[second_set_offset], &test_psynw.fxalt_ovrf[second_set_offset],
                  &test_psynw.fyalt_ovrf[second_set_offset],
                  &test_psynw.fzalt_ovrf[second_set_offset], test_psynw.frc_scale,
                  buffer_a_cfw.xcrd, buffer_a_cfw.ycrd, buffer_a_cfw.zcrd, nullptr, nullptr,
                  nullptr, 1.0, buffer_a_cfw.natom,
                  "CoordinateFrame => PhaseSpaceSynthesis(FORCES, BLACK)",
                  tsm.getTestingStatus(), 1.5e-7, 2.9e-7);

  // Create more test systems with periodic boundary conditions.  Test that the transformation
  // matrices are transferred properly, while also checking other combinations of origin and
  // destination for coordCopy() overloads.
  const std::vector<std::string> pbc_fi_names = { "symmetry_C1_in_water", "symmetry_C2_in_water",
                                                  "symmetry_C3_in_water", "symmetry_C4_in_water" };
  TestSystemManager tsm_pbc(base_top_dir, "top", pbc_fi_names, base_crd_dir, "inpcrd",
                            pbc_fi_names);
  const TestPriority do_pbc_tests = tsm_pbc.getTestingStatus();
  std::vector<CoordinateFrame> tsm_pbc_cf;
  std::vector<PhaseSpace> tsm_pbc_ps;
  std::vector<CoordinateSeries<double>> tsm_pbc_cs;
  std::vector<CoordinateFrameReader> tsm_pbc_cfr;
  std::vector<PhaseSpaceReader> tsm_pbc_psr;
  std::vector<CoordinateSeriesReader<double>> tsm_pbc_csr;
  PhaseSpaceSynthesis tsm_pbc_psyn = spawnMutableCoordinateObjects(tsm_pbc, &tsm_pbc_cf,
                                                                   &tsm_pbc_cfr, &tsm_pbc_ps,
                                                                   &tsm_pbc_psr, &tsm_pbc_cs,
                                                                   &tsm_pbc_csr);
  Condensate tsm_pbc_cdns(tsm_pbc_cs[2], PrecisionModel::DOUBLE);
  PhaseSpace recv_pbc(tsm_pbc_cfr[2].natom);
  PhaseSpaceReader recv_pbc_r = recv_pbc.data();
  check(tsm_pbc_cdns.getBasis() == StructureSource::SERIES, "A condensate based on a coordinate "
        "series does not properly report its basis.", do_pbc_tests);
  check(tsm_pbc_cdns.ownsCoordinates() == false, "A condensate based on a coordinate series of "
        "double-precision reals should rely on the source to store its data, but does not.",
        do_pbc_tests);
  coordCopy(&recv_pbc, TrajectoryKind::POSITIONS, CoordinateCycle::BLACK, tsm_pbc_cdns, 3);
  coordCopy(&recv_pbc, TrajectoryKind::FORCES, CoordinateCycle::BLACK, tsm_pbc_cdns, 2);
  const CoordinateFrame sm2_f3 = tsm_pbc_cs[2].exportFrame(3);
  const CoordinateFrameReader sm2_f3r = sm2_f3.data();
  diffCoordinates(recv_pbc_r.xalt, recv_pbc_r.yalt, recv_pbc_r.zalt, sm2_f3r.xcrd, sm2_f3r.ycrd,
                  sm2_f3r.zcrd, sm2_f3r.natom, recv_pbc_r.umat_alt, sm2_f3r.umat,
                  "CoordinateSeries => PhaseSpace(POSITIONS, BLACK)", do_pbc_tests);
  const CoordinateFrame sm2_f2 = tsm_pbc_cs[2].exportFrame(2);
  const CoordinateFrameReader sm2_f2r = sm2_f2.data();
  diffCoordinates(recv_pbc_r.fxalt, recv_pbc_r.fyalt, recv_pbc_r.fzalt, sm2_f2r.xcrd, sm2_f2r.ycrd,
                  sm2_f2r.zcrd, sm2_f2r.natom,
                  "CoordinateSeries => PhaseSpace(FORCES, BLACK)", do_pbc_tests);
  
#ifdef STORMM_USE_HPC
  // The following unit tests use HPC resources, copying between the host and device layers.
  // Start by getting details of the available GPU.
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);  

  // Make additional copies of the periodic coordinate objects, distinguished by perturbations.
  std::vector<CoordinateFrame> clone_pbc_cf;
  std::vector<PhaseSpace> clone_pbc_ps;
  std::vector<CoordinateSeries<double>> clone_pbc_cs;
  std::vector<CoordinateFrameReader> clone_pbc_cfr;
  std::vector<PhaseSpaceReader> clone_pbc_psr;
  std::vector<CoordinateSeriesReader<double>> clone_pbc_csr;
  PhaseSpaceSynthesis clone_pbc_psyn = spawnMutableCoordinateObjects(tsm_pbc, &clone_pbc_cf,
                                                                     &clone_pbc_cfr, &clone_pbc_ps,
                                                                     &clone_pbc_psr, &clone_pbc_cs,
                                                                     &clone_pbc_csr, 187693308);
  Condensate clone_pbc_cdns(clone_pbc_psyn, PrecisionModel::SINGLE);
  std::vector<CoordinateFrameWriter> host_clone_cfw;
  std::vector<PhaseSpaceWriter> host_clone_psw;
  std::vector<CoordinateSeriesWriter<double>> host_clone_csw;
  PsSynthesisWriter host_clone_psynw = clone_pbc_psyn.data();
  PsSynthesisWriter host_tsm_psynw = tsm_pbc_psyn.data();
  CondensateWriter host_clone_cdnsw = clone_pbc_cdns.data();
  for (int i = 0; i < tsm_pbc.getSystemCount(); i++) {
    host_clone_cfw.push_back(clone_pbc_cf[i].data());
    host_clone_psw.push_back(clone_pbc_ps[i].data());
    host_clone_csw.push_back(clone_pbc_cs[i].data());
    const int natom = host_clone_cfw[i].natom;
    const int padded_natom = roundUp(natom, warp_size_int);
    addRandomNoise(&xrs, host_clone_cfw[i].xcrd, host_clone_cfw[i].ycrd, host_clone_cfw[i].zcrd,
                   natom, 0.1);
    perturbPhaseSpace(&xrs, &host_clone_psw[i]);
    for (int j = 0; j < host_clone_csw[i].nframe; j++) {
      const int jpos = j * padded_natom;
      addRandomNoise(&xrs, &host_clone_csw[i].xcrd[jpos], &host_clone_csw[i].ycrd[jpos],
                     &host_clone_csw[i].zcrd[jpos], natom, 0.1);
    }
    const int istart = host_clone_psynw.atom_starts[i];
    perturbSynthesis(&xrs, &host_clone_psynw, istart, natom);
    perturbSynthesis(&xrs, &host_tsm_psynw, istart, natom);
    addRandomNoise(&xrs, &host_clone_cdnsw.xcrd_sp[istart], &host_clone_cdnsw.ycrd_sp[istart],
                   &host_clone_cdnsw.zcrd_sp[istart], natom, 0.1);

    // After these uploads, all structures will hold the same coordinates on both the host and
    // the device, but each structure will have unique perturbations.
    clone_pbc_cf[i].upload();
    clone_pbc_ps[i].upload();
    clone_pbc_cs[i].upload();
  }
  clone_pbc_psyn.upload();
  tsm_pbc_psyn.upload();
  clone_pbc_cdns.upload();

  // The various coordinate objects now contain unique structures, but the host and device memory
  // for each particular object contains identical copies of one structure.  Run tests of various
  // copy operations to verify their fidelity.
  const HybridTargetLevel hl_host = HybridTargetLevel::HOST;
  const HybridTargetLevel hl_devc = HybridTargetLevel::DEVICE;
  const std::vector<HybridTargetLevel> dest_tiers = { hl_devc, hl_host, hl_devc };
  const std::vector<HybridTargetLevel> orig_tiers = { hl_host, hl_devc, hl_devc };
  for (int i = 0; i < 3; i++) {
    copySetup(clone_pbc_cf[2], tsm_pbc_cf[2], dest_tiers[i], orig_tiers[i], gpu, do_pbc_tests);
    copySetup(clone_pbc_cf[1], tsm_pbc_ps[1], dest_tiers[i], orig_tiers[i], gpu, do_pbc_tests);
    copySetup(clone_pbc_cf[0], tsm_pbc_cs[0], dest_tiers[i], orig_tiers[i], gpu, do_pbc_tests);
    copySetup(clone_pbc_cf[3], clone_pbc_cdns, 3, dest_tiers[i], orig_tiers[i], gpu, do_pbc_tests);
    copySetup(clone_pbc_cf[1], clone_pbc_psyn, 1, dest_tiers[i], orig_tiers[i], gpu, do_pbc_tests);
    copySetup(clone_pbc_ps[3], tsm_pbc_ps[3], dest_tiers[i], orig_tiers[i], gpu, do_pbc_tests);
    copySetup(clone_pbc_ps[2], tsm_pbc_cs[2], dest_tiers[i], orig_tiers[i], gpu, do_pbc_tests);
    copySetup(clone_pbc_ps[1], clone_pbc_cdns, 1, dest_tiers[i], orig_tiers[i], gpu, do_pbc_tests);
    copySetup(clone_pbc_ps[0], clone_pbc_psyn, 0, dest_tiers[i], orig_tiers[i], gpu, do_pbc_tests);
    copySetup(clone_pbc_cs[3], tsm_pbc_cs[3], dest_tiers[i], orig_tiers[i], gpu, do_pbc_tests);
    copySetup(clone_pbc_cs[2], clone_pbc_psyn, 2, dest_tiers[i], orig_tiers[i], gpu, do_pbc_tests);
    copySetup(clone_pbc_cs[1], clone_pbc_cdns, 1, dest_tiers[i], orig_tiers[i], gpu, do_pbc_tests);
    copySetup(clone_pbc_psyn, 0, clone_pbc_cdns, 4, dest_tiers[i], orig_tiers[i], gpu,
              do_pbc_tests);
    copySetup(clone_pbc_psyn, 5, tsm_pbc_psyn, 1, dest_tiers[i], orig_tiers[i], gpu,
              do_pbc_tests);
    copySetup(clone_pbc_cdns, 2, clone_pbc_cdns, 6, dest_tiers[i], orig_tiers[i], gpu,
              do_pbc_tests);
  }
#endif
  
  // Print results
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
