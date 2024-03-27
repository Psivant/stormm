#include <string>
#include <vector>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hpc_config.h"
#include "Accelerator/hybrid.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "FileManagement/file_listing.h"
#include "Math/rounding.h"
#include "Math/series_ops.h"
#include "Numerics/split_fixed_precision.h"
#include "Random/random.h"
#include "Reporting/summary_file.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Trajectory/coordinate_copy.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/coordinate_swap.h"
#include "Trajectory/trajectory_enumerators.h"
#include "UnitTesting/approx.h"
#include "UnitTesting/test_environment.h"
#include "UnitTesting/test_system_manager.h"
#include "UnitTesting/unit_test.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::numerics;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Save the box dimensions of a multi-system coordinate object to temporary buffers.
//
// Arguments:
//   umat:    Pointer to the array of transformation matrices taking Cartesian coordinates into
//            unit cell fractional space
//   invu:    Pointer to the array of inverse transformation matrices taking unit cell fractional
//            coordinates into Cartesian space
//   bdim:    Pointer to the array of box dimensions (these, as well as the transformation
//            matrices, are arranged in separate cache lines padded by the warp size)
//   nfrm:    The number of frames or system indices in the multi-coordinate object
//   utmp:    Temporary array of transformation matrices, holding the packed contents of umat.
//            Allocated, filled and returned.
//   itmp:    Temporary array of inverse transformation matrices, holding the packed contents of
//            invu.  Allocated, filled and returned.
//   btmp:    Temporary array of box dimensions, holding the packed contents of bdim.  Allocated,
//            filled and returned.
//-------------------------------------------------------------------------------------------------
void saveBoxDimensions(const double* umat, const double* invu, const double* bdim,
                       const size_t nfrm, std::vector<double> *utmp, std::vector<double> *itmp,
                       std::vector<double> *btmp) {
  const size_t xfrm_offset = roundUp(9, warp_size_int); 
  const size_t bdim_offset = roundUp(6, warp_size_int); 
  utmp->resize(9 * nfrm);
  itmp->resize(9 * nfrm);
  btmp->resize(6 * nfrm);
  double* u_ptr = utmp->data();
  double* i_ptr = itmp->data();
  double* b_ptr = btmp->data();
  for (size_t i = 0; i < nfrm; i++) {
    for (size_t j = 0; j < 9; j++) {
      u_ptr[(i * 9) + j] = umat[(i * xfrm_offset) + j];
      i_ptr[(i * 9) + j] = invu[(i * xfrm_offset) + j];
    }
    for (size_t j = 0; j < 6; j++) {
      b_ptr[(i * 6) + j] = bdim[(i * bdim_offset) + j];
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Restore the box dimensions of a multi-system coordinate object from temporary buffers.
// Descriptions of input parameters follow from saveBoxDimensions(), above.
//-------------------------------------------------------------------------------------------------
void restoreBoxDimensions(double* umat, double* invu, double* bdim, const size_t nfrm,
                          const std::vector<double> &utmp, const std::vector<double> &itmp,
                          const std::vector<double> &btmp) {
  const size_t xfrm_offset = roundUp(9, warp_size_int); 
  const size_t bdim_offset = roundUp(6, warp_size_int); 
  for (size_t i = 0; i < nfrm; i++) {
    for (size_t j = 0; j < 9; j++) {
      umat[(i * xfrm_offset) + j] = utmp[(i * 9) + j];
      invu[(i * xfrm_offset) + j] = itmp[(i * 9) + j];
    }
    for (size_t j = 0; j < 6; j++) {
      bdim[(i * bdim_offset) + j] = btmp[(i * 6) + j];
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Blur the box dimensions by random resizings, maintaining the original unit cell angles.
//
// Arguments:
//   umat:    Pointer to the array of transformation matrices taking Cartesian coordinates into
//            unit cell fractional space
//   invu:    Pointer to the array of inverse transformation matrices taking unit cell fractional
//            coordinates into Cartesian space
//   bdim:    Pointer to the array of box dimensions (these, as well as the transformation
//            matrices, are arranged in separate cache lines padded by the warp size)
//   nfrm:    The number of frames or system indices in the multi-coordinate object
//   xrs:     Random number generator to introduce instability
//-------------------------------------------------------------------------------------------------
void stretchBoxes(double* umat, double* invu, double* bdim, const size_t nfrm,
                  Xoroshiro128pGenerator *xrs) {
  const size_t bstride = roundUp(6, warp_size_int);
  const size_t xstride = roundUp(9, warp_size_int);
  for (size_t i = 0; i < nfrm; i++) {
    for (size_t j = 0; j < 3; j++) {
      bdim[(i * bstride) + j] *= (1.0 + xrs->gaussianRandomNumber());
    }
    computeBoxTransform(bdim[(i * bstride)    ], bdim[(i * bstride) + 1], bdim[(i * bstride) + 2],
                        bdim[(i * bstride) + 3], bdim[(i * bstride) + 4], bdim[(i * bstride) + 5],
                        &umat[(i * xstride)], &invu[(i * xstride)]);
  }
}

//-------------------------------------------------------------------------------------------------
// Perturb the coordinates of a coordinate object.
//
// Overloaded:
//   - Operate on two CoordinateSeries objects
//   - Operate on two PhaseSpaceSynthesis objects
//   - Operate on two Condensate objects
//
// Arguments:
//   c:       The coordinate object to modify
//   xrs:     Random number generator to introduce instability
//-------------------------------------------------------------------------------------------------
template <typename T>
void blurCoordinates(CoordinateSeriesWriter<T> *cw, Xoroshiro128pGenerator *xrs) {
  const size_t natom_zu = static_cast<size_t>(cw->natom);
  const size_t nfrm_zu  = static_cast<size_t>(cw->nframe);
  for (size_t i = 0; i < nfrm_zu; i++) {
    const size_t offset = static_cast<size_t>(i) * natom_zu;
    addRandomNoise(xrs, &cw->xcrd[offset], &cw->ycrd[offset], &cw->zcrd[offset], natom_zu, 0.5,
                   cw->gpos_scale);
  }
  stretchBoxes(cw->umat, cw->invu, cw->boxdim, cw->nframe, xrs);
}

void blurCoordinates(PsSynthesisWriter *cw, Xoroshiro128pGenerator *xrs) {
  for (int i = 0; i < cw->system_count; i++) {
    const size_t offset = cw->atom_starts[i];
    addRandomNoise(xrs, &cw->xcrd[offset], &cw->xcrd_ovrf[offset], &cw->ycrd[offset],
                   &cw->ycrd_ovrf[offset], &cw->zcrd[offset], &cw->zcrd_ovrf[offset],
                   cw->atom_counts[i], 0.5, cw->gpos_scale);
    addRandomNoise(xrs, &cw->xvel[offset], &cw->xvel_ovrf[offset], &cw->yvel[offset],
                   &cw->yvel_ovrf[offset], &cw->zvel[offset], &cw->zvel_ovrf[offset],
                   cw->atom_counts[i], 0.5, cw->gpos_scale);
    addRandomNoise(xrs, &cw->xfrc[offset], &cw->xfrc_ovrf[offset], &cw->yfrc[offset],
                   &cw->yfrc_ovrf[offset], &cw->zfrc[offset], &cw->zfrc_ovrf[offset],
                   cw->atom_counts[i], 0.5, cw->gpos_scale);
  }
  stretchBoxes(cw->umat, cw->invu, cw->boxdims, cw->system_count, xrs);
}

void blurCoordinates(CondensateWriter *cw, Xoroshiro128pGenerator *xrs) {
  for (int i = 0; i < cw->system_count; i++) {
    const size_t offset = cw->atom_starts[i];
    switch (cw->mode) {
    case PrecisionModel::DOUBLE:
      addRandomNoise(xrs, &cw->xcrd[offset], &cw->ycrd[offset], &cw->zcrd[offset],
                     cw->atom_counts[i], 0.5);
      break;
    case PrecisionModel::SINGLE:
      addRandomNoise(xrs, &cw->xcrd_sp[offset], &cw->ycrd_sp[offset], &cw->zcrd_sp[offset],
                     cw->atom_counts[i], 0.5);
      break;
    }
  }
  stretchBoxes(cw->umat, cw->invu, cw->boxdims, cw->system_count, xrs);
}

//-------------------------------------------------------------------------------------------------
// Introduce perturbations to the coordinates in some coordinate object.  Upload the results if
// applicable and further perturb coordinates on the host side to differentiate all coordinate
// sets anywhere in the object.  Overloading and descriptions of parameters follow from
// blurCoordinates() above.
//-------------------------------------------------------------------------------------------------
template <typename T>
void differentiateCoordinates(CoordinateSeries<T> *c, Xoroshiro128pGenerator *xrs) {
  CoordinateSeriesWriter<T> cw = c->data();
  const size_t natom_zu = static_cast<size_t>(cw.natom);
  const size_t nfrm_zu  = static_cast<size_t>(cw.nframe);
#ifdef STORMM_USE_HPC
  const size_t padded_natom = roundUp(natom_zu, warp_size_zu);
  std::vector<double> xtmp(natom_zu * nfrm_zu);
  std::vector<double> ytmp(natom_zu * nfrm_zu);
  std::vector<double> ztmp(natom_zu * nfrm_zu);
  std::vector<double> utmp, itmp, btmp;
  for (size_t i = 0; i < nfrm_zu; i++) {
    for (size_t j = 0; j < natom_zu; j++) {
      xtmp[(i * natom_zu) + j] = cw.xcrd[(i * padded_natom) + j];
      ytmp[(i * natom_zu) + j] = cw.ycrd[(i * padded_natom) + j];
      ztmp[(i * natom_zu) + j] = cw.zcrd[(i * padded_natom) + j];
    }
  }
  saveBoxDimensions(cw.umat, cw.invu, cw.boxdim, cw.nframe, &utmp, &itmp, &btmp);
#endif
  blurCoordinates(&cw, xrs);
#ifdef STORMM_USE_HPC
  c->upload();
  for (size_t i = 0; i < nfrm_zu; i++) {
    for (size_t j = 0; j < natom_zu; j++) {
      cw.xcrd[(i * padded_natom) + j] = xtmp[(i * natom_zu) + j];
      cw.ycrd[(i * padded_natom) + j] = ytmp[(i * natom_zu) + j];
      cw.zcrd[(i * padded_natom) + j] = ztmp[(i * natom_zu) + j];
    }
  }
  restoreBoxDimensions(cw.umat, cw.invu, cw.boxdim, cw.nframe, utmp, itmp, btmp);
  blurCoordinates(&cw, xrs);
#endif
}

void differentiateCoordinates(PhaseSpaceSynthesis *c, Xoroshiro128pGenerator *xrs) {
  PsSynthesisWriter cw = c->data();
#ifdef STORMM_USE_HPC
  const size_t last_sys = cw.system_count - 1;
  const size_t total_atoms = cw.atom_starts[last_sys] + cw.atom_counts[last_sys];
  std::vector<double> xtmp(total_atoms), ytmp(total_atoms), ztmp(total_atoms);
  std::vector<double> utmp, itmp, btmp;
  for (size_t i = 0; i < total_atoms; i++) {
    xtmp[i] = hostInt95ToDouble(cw.xcrd[i], cw.xcrd_ovrf[i]);
    ytmp[i] = hostInt95ToDouble(cw.ycrd[i], cw.ycrd_ovrf[i]);
    ztmp[i] = hostInt95ToDouble(cw.zcrd[i], cw.zcrd_ovrf[i]);
  }
  saveBoxDimensions(cw.umat, cw.invu, cw.boxdims, cw.system_count, &utmp, &itmp, &btmp);
#endif
  blurCoordinates(&cw, xrs);
#ifdef STORMM_USE_HPC
  c->upload();

  // Replacement of the original coordinates can happen, to within any reasonable precision, by
  // simply swapping back the values "compressed" into the double-precision array.
  for (size_t i = 0; i < total_atoms; i++) {
    const int95_t xrep = hostDoubleToInt95(xtmp[i]);
    const int95_t yrep = hostDoubleToInt95(ytmp[i]);
    const int95_t zrep = hostDoubleToInt95(ztmp[i]);
    cw.xcrd[i] = xrep.x;
    cw.ycrd[i] = yrep.x;
    cw.zcrd[i] = zrep.x;
    cw.xcrd_ovrf[i] = xrep.y;
    cw.ycrd_ovrf[i] = yrep.y;
    cw.zcrd_ovrf[i] = zrep.y;
    cw.xvel[i] = 0LL;
    cw.yvel[i] = 0LL;
    cw.zvel[i] = 0LL;
    cw.xvel_ovrf[i] = 0;
    cw.yvel_ovrf[i] = 0;
    cw.zvel_ovrf[i] = 0;
    cw.xfrc[i] = 0LL;
    cw.yfrc[i] = 0LL;
    cw.zfrc[i] = 0LL;
    cw.xfrc_ovrf[i] = 0;
    cw.yfrc_ovrf[i] = 0;
    cw.zfrc_ovrf[i] = 0;
  }
  restoreBoxDimensions(cw.umat, cw.invu, cw.boxdims, cw.system_count, utmp, itmp, btmp);
  blurCoordinates(&cw, xrs);
#endif
}

void differentiateCoordinates(Condensate *c, Xoroshiro128pGenerator *xrs) {
  CondensateWriter cw = c->data();
#ifdef STORMM_USE_HPC
  const size_t last_sys = cw.system_count - 1;
  const size_t total_atoms = cw.atom_starts[last_sys] + cw.atom_counts[last_sys];
  std::vector<double> xtmp(total_atoms), ytmp(total_atoms), ztmp(total_atoms);
  std::vector<double> utmp, itmp, btmp;
  switch (cw.mode) {
  case PrecisionModel::DOUBLE:
    for (size_t i = 0; i < total_atoms; i++) {
      xtmp[i] = cw.xcrd[i];
      ytmp[i] = cw.ycrd[i];
      ztmp[i] = cw.zcrd[i];
    }
    break;
  case PrecisionModel::SINGLE:
    for (size_t i = 0; i < total_atoms; i++) {
      xtmp[i] = cw.xcrd_sp[i];
      ytmp[i] = cw.ycrd_sp[i];
      ztmp[i] = cw.zcrd_sp[i];
    }
    break;
  }
  saveBoxDimensions(cw.umat, cw.invu, cw.boxdims, cw.system_count, &utmp, &itmp, &btmp);
#endif
  blurCoordinates(&cw, xrs);
#ifdef STORMM_USE_HPC
  c->upload();

  // Double-precision arrays will have stored in the original coordinates with no approximation.
  switch (cw.mode) {
  case PrecisionModel::DOUBLE:
    for (size_t i = 0; i < total_atoms; i++) {
      cw.xcrd[i] = xtmp[i];
      cw.ycrd[i] = ytmp[i];
      cw.zcrd[i] = ztmp[i];
    }
    break;
  case PrecisionModel::SINGLE:
    for (size_t i = 0; i < total_atoms; i++) {
      cw.xcrd_sp[i] = xtmp[i];
      cw.ycrd_sp[i] = ytmp[i];
      cw.zcrd_sp[i] = ztmp[i];
    }
    break;
  }
  restoreBoxDimensions(cw.umat, cw.invu, cw.boxdims, cw.system_count, utmp, itmp, btmp);
  blurCoordinates(&cw, xrs);
#endif
}

//-------------------------------------------------------------------------------------------------
// Extract the box information from a CoordinateFrame as a compact vector of 24 values.
//
// Arguments:
//   cf:    The coordinate frame to examine
//-------------------------------------------------------------------------------------------------
std::vector<double> getBoxData(const CoordinateFrame &cf) {
  const CoordinateFrameReader cfr = cf.data();
  std::vector<double> result(24);
  for (int i = 0; i < 9; i++) {
    result[i    ] = cfr.umat[i];
    result[i + 9] = cfr.invu[i];
  }
  for (int i = 0; i < 6; i++) {
    result[i + 18] = cfr.boxdim[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Convert a string of failed system comparisons into a parseable error message, with proper
// grammar.
//
// Arguments:
//   failures:  The list of failed comparisons
//-------------------------------------------------------------------------------------------------
std::string listErrorIndices(const std::vector<int> &failures) {
  std::string result = (failures.size() == 1) ? "index " : "indices ";
  const int nfail = failures.size();
  for (size_t i = 0; i < nfail; i++) {
    result += std::to_string(failures[i]);
    if (i < nfail - 1 && nfail > 2) {
      result += ", ";
      if (i == nfail - 2) {
        result += "and ";
      }
    }
    else if (i == 0 && nfail == 2) {
      result += " and ";
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Check that two coordinate systems (including velocities and forces, if applicable) are
// identical.  Extend a vector of failed comparisons if not.
//
// Overloaded:
//   - Operate on two CoordinateSeries objects
//   - Operate on two PhaseSpaceSynthesis objects
//   - Operate on two Condensate objects
//
// Arguments:
//   init:          The initial state of the originating object (this could be the first or second
//                  coordinate object involved in the swap)
//   finl:          The final state of the destination coordinate object (this would be the second
//                  or first coordinate object involved in the swap)
//   init_idx:      Index of the structure in the originating object
//   finl_idx:      Index of the structure in the destination object
//   init_tier:     Tier of the originating object where the structure comes from
//   finl_tier:     Tier of the destination object where the structure goes to
//   failures:      A list of failed comparisons for positions, velocities, or forces, modified
//                  and returned if any failures occur
//   box_failures:  A list of failed comparisons for box geometries, modified and returned if any
//                  failures occur
//   gpu:           Details of the GPU that carried out the swap(s)
//-------------------------------------------------------------------------------------------------
template <typename T>
void verifyTransfer(const CoordinateSeries<T> &init, const CoordinateSeries<T> &finl,
                    CoordinateFrame *workspace, const int init_idx, const int finl_idx,
                    const HybridTargetLevel init_tier, const HybridTargetLevel finl_tier,
                    std::vector<int> *failures, std::vector<int> *box_failures,
                    const GpuDetails &gpu = null_gpu) {
  const HybridTargetLevel host_tier = HybridTargetLevel::HOST;
  coordCopy(workspace, init, init_idx, host_tier, init_tier, gpu);
  const std::vector<double> ref_xyz = workspace->getInterlacedCoordinates();
  const std::vector<double> ref_box = getBoxData(*workspace);
  coordCopy(workspace, finl, finl_idx, host_tier, finl_tier, gpu);
  const std::vector<double> end_xyz = workspace->getInterlacedCoordinates();
  const std::vector<double> end_box = getBoxData(*workspace);
  const Approx chk_first(ref_xyz);
  if (chk_first.test(end_xyz) == false) {
    failures->push_back(init_idx);
  }
  const Approx chk_first_box(ref_box);
  if (chk_first_box.test(end_box) == false) {
    box_failures->push_back(init_idx);
  }
}

void verifyTransfer(const PhaseSpaceSynthesis &init, const PhaseSpaceSynthesis &finl,
                    CoordinateFrame *workspace, const int init_idx, const int finl_idx,
                    const TrajectoryKind kind, const HybridTargetLevel init_tier,
                    const HybridTargetLevel finl_tier, std::vector<int> *failures,
                    std::vector<int> *box_failures = nullptr, const GpuDetails &gpu = null_gpu) {
  const HybridTargetLevel host_tier = HybridTargetLevel::HOST;
  coordCopy(workspace, init, init_idx, kind, host_tier, init_tier, gpu);
  const std::vector<double> ref_xyz = workspace->getInterlacedCoordinates();
  const std::vector<double> ref_box = getBoxData(*workspace);
  coordCopy(workspace, finl, finl_idx, kind, host_tier, finl_tier, gpu);
  const std::vector<double> end_xyz = workspace->getInterlacedCoordinates();
  const std::vector<double> end_box = getBoxData(*workspace);
  const Approx chk_first(ref_xyz);
  if (chk_first.test(end_xyz) == false) {
    failures->push_back(init_idx);
  }
  if (box_failures != nullptr) {
    const Approx chk_first_box(ref_box);
    if (chk_first_box.test(end_box) == false) {
      box_failures->push_back(init_idx);
    }
  }
}

void verifyTransfer(const Condensate &init, const Condensate &finl, CoordinateFrame *workspace,
                    const int init_idx, const int finl_idx, const HybridTargetLevel init_tier,
                    const HybridTargetLevel finl_tier, std::vector<int> *failures,
                    std::vector<int> *box_failures, const GpuDetails &gpu = null_gpu) {
  const HybridTargetLevel host_tier = HybridTargetLevel::HOST;
  coordCopy(workspace, init, init_idx, host_tier, init_tier, gpu);
  const std::vector<double> ref_xyz = workspace->getInterlacedCoordinates();
  const std::vector<double> ref_box = getBoxData(*workspace);
  coordCopy(workspace, finl, finl_idx, host_tier, finl_tier, gpu);
  const std::vector<double> end_xyz = workspace->getInterlacedCoordinates();
  const std::vector<double> end_box = getBoxData(*workspace);
  const Approx chk_first(ref_xyz);
  if (chk_first.test(end_xyz) == false) {
    failures->push_back(init_idx);
  }
  const Approx chk_first_box(ref_box);
  if (chk_first_box.test(end_box) == false) {
    box_failures->push_back(init_idx);
  }
}

//-------------------------------------------------------------------------------------------------
// Check that information was properly conveyed over course of a swap.  This processes lists of
// failures accumulated by various functions above.
//
// Arguments:
//   forward:   List of failures taking systems in the first object to the second
//   backward:  List of failures taking systems in the second object to the first
//   desc:      Description of the type of information being conveyed
//   tier_a:    Indicate whether the first object's data resides on the CPU host or GPU device
//   tier_b:    Indicate whether the first object's data resides on the GPU device or CPU host
//   do_tests:  Indicate whether tests are viable (based on the availability of files needed to
//              set up the tests)
//-------------------------------------------------------------------------------------------------
void testSwapFidelity(const std::vector<int> &forward, const std::vector<int> &backward,
                      const std::string &desc, const HybridTargetLevel tier_a,
                      const HybridTargetLevel tier_b, const TestPriority do_tests) {
  const std::string forward_list = listErrorIndices(forward);
  const std::string backward_list = listErrorIndices(backward);
  check(forward.size() == 0, "Some " + desc + " swaps were unsuccessful.  System " +
        forward_list + " residing in " + getEnumerationName(tier_a) + " memory of the first "
        "object were not transferred to " + getEnumerationName(tier_b) + " memory of the second "
        "object.", do_tests);
  check(backward.size() == 0, "Some " + desc + " swaps were unsuccessful.  System " +
        backward_list + " residing in " + getEnumerationName(tier_a) + " memory of the second "
        "object were not transferred to " + getEnumerationName(tier_b) + " memory of the first "
        "object.", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Take snapshots of two coordinate objects, perform a coordinate swap, and compare the results.
//
// Overloaded:
//   - Operate on two CoordinateSeries objects
//   - Operate on two PhaseSpaceSynthesis objects
//   - Operate on two Condensate objects
//
// Arguments:
//   a:         The first object participating in the swap
//   b:         The second object participating in the swap
//   swaps:     A list of pairs of systems to swap.  If the list has only one element, the
//              appropriate single-system variant of coordSwap will be called.
//   do_tests:  Indicate whether tests are viable
//   tier_a:    The tier of a to target in the swap
//   tier_b:    The tier of b to target in the swap
//-------------------------------------------------------------------------------------------------
template <typename T>
void checkSwap(CoordinateSeries<T> *a, CoordinateSeries<T> *b, const Hybrid<int2> &swaps,
               const TestPriority do_tests,
               const HybridTargetLevel tier_a = HybridTargetLevel::HOST,
               const HybridTargetLevel tier_b = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu) {
  const int natom = a->getAtomCount();
  const CoordinateSeries<T> init_a(*a);
  const CoordinateSeries<T> init_b(*b);
  const int nswaps = swaps.size();
  if (nswaps == 1) {
    const int2 pair = swaps.readHost(0);
    coordSwap(a, pair.x, b, pair.y, tier_a, tier_b, gpu);
  }
  else {
    coordSwap(a, b, swaps, tier_a, tier_b, gpu);
  }
  
  // Loop through the swaps and check that the current state of the input coordinate series is
  // correct in terms of their original contents.
  CoordinateFrame exam(natom);
  std::vector<int> first_failures, second_failures, first_box_failures, second_box_failures;
  for (int i = 0; i < nswaps; i++) {
    const int2 pair = swaps.readHost(i);
    verifyTransfer(init_a, *b, &exam, pair.x, pair.y, tier_a, tier_b, &first_failures,
                   &first_box_failures);
    verifyTransfer(init_b, *a, &exam, pair.y, pair.x, tier_b, tier_a, &second_failures,
                   &second_box_failures);
  }
  testSwapFidelity(first_failures, second_failures, "position", tier_a, tier_b, do_tests);
  testSwapFidelity(first_box_failures, second_box_failures, "box geometry", tier_a, tier_b,
                   do_tests);
}

void checkSwap(PhaseSpaceSynthesis *a, PhaseSpaceSynthesis *b, const Hybrid<int2> &swaps,
               const TestPriority do_tests,
               const HybridTargetLevel tier_a = HybridTargetLevel::HOST,
               const HybridTargetLevel tier_b = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu) {
  const PhaseSpaceSynthesis init_a(*a);
  const PhaseSpaceSynthesis init_b(*b);
  const int nswaps = swaps.size();
  if (nswaps == 1) {
    const int2 pair = swaps.readHost(0);
    coordSwap(a, pair.x, b, pair.y, tier_a, tier_b, gpu);
  }
  else {
    coordSwap(a, b, swaps, tier_a, tier_b, gpu);
  }
  
  // Loop through the swaps and check that the current state of the input coordinate series is
  // correct in terms of their original contents.
  std::vector<int> fpos_failures, spos_failures, fbox_failures, sbox_failures;
  std::vector<int> fvel_failures, svel_failures, ffrc_failures, sfrc_failures;
  for (int i = 0; i < nswaps; i++) {
    const int2 pair = swaps.readHost(i);
    CoordinateFrame exam(init_a.getAtomCount(pair.x));
    verifyTransfer(init_a, *b, &exam, pair.x, pair.y, TrajectoryKind::POSITIONS, tier_a, tier_b,
                   &fpos_failures, &fbox_failures);
    verifyTransfer(init_b, *a, &exam, pair.y, pair.x, TrajectoryKind::POSITIONS, tier_b, tier_a,
                   &spos_failures, &sbox_failures);
    verifyTransfer(init_a, *b, &exam, pair.x, pair.y, TrajectoryKind::VELOCITIES, tier_a, tier_b,
                   &fvel_failures);
    verifyTransfer(init_b, *a, &exam, pair.y, pair.x, TrajectoryKind::VELOCITIES, tier_b, tier_a,
                   &svel_failures);
    verifyTransfer(init_a, *b, &exam, pair.x, pair.y, TrajectoryKind::FORCES, tier_a, tier_b,
                   &ffrc_failures);
    verifyTransfer(init_b, *a, &exam, pair.y, pair.x, TrajectoryKind::FORCES, tier_b, tier_a,
                   &sfrc_failures);
  }
  testSwapFidelity(fpos_failures, spos_failures, "position", tier_a, tier_b, do_tests);
  testSwapFidelity(fbox_failures, sbox_failures, "box geometry", tier_a, tier_b, do_tests);
  testSwapFidelity(fvel_failures, svel_failures, "velocity", tier_a, tier_b, do_tests);
  testSwapFidelity(ffrc_failures, sfrc_failures, "force", tier_a, tier_b, do_tests);
}

void checkSwap(Condensate *a, Condensate *b, const Hybrid<int2> &swaps,
               const TestPriority do_tests,
               const HybridTargetLevel tier_a = HybridTargetLevel::HOST,
               const HybridTargetLevel tier_b = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu) {
  const Condensate init_a(*a);
  const Condensate init_b(*b);
  const int nswaps = swaps.size();
  if (nswaps == 1) {
    const int2 pair = swaps.readHost(0);
    coordSwap(a, pair.x, b, pair.y, tier_a, tier_b, gpu);
  }
  else {
    coordSwap(a, b, swaps, tier_a, tier_b, gpu);
  }
  
  // Loop through the swaps and check that the current state of the input coordinate series is
  // correct in terms of their original contents.
  std::vector<int> fpos_failures, spos_failures, fbox_failures, sbox_failures;
  for (int i = 0; i < nswaps; i++) {
    const int2 pair = swaps.readHost(i);
    CoordinateFrame exam(init_a.getAtomCount(pair.x));
    verifyTransfer(init_a, *b, &exam, pair.x, pair.y, tier_a, tier_b, &fpos_failures,
                   &fbox_failures);
    verifyTransfer(init_b, *a, &exam, pair.y, pair.x, tier_b, tier_a, &spos_failures,
                   &sbox_failures);
  }
  testSwapFidelity(fpos_failures, spos_failures, "position", tier_a, tier_b, do_tests);
  testSwapFidelity(fbox_failures, sbox_failures, "box geometry", tier_a, tier_b, do_tests);
}

//-------------------------------------------------------------------------------------------------
// Reset the scheduled swaps for a particular array.
//
// Overloaded:
//   - Assume that all frames have the same atom count (valid for CoordinateSeries objects)
//   - Accept lists of the atom counts for systems in each object (valid for either type of
//     coordinate synthesis object)
//
// Arguments:
//   swaps:    The array of swaps
//   nsys_a:   The number of systems in the first object that will be involved in the swap
//   nsys_b:   The number of systems in the second object that will be involved in the swap (if
//             this is less than zero, it indicates that the two objects are one and the same)
//   natom_a:  The number of atoms in each system of the first synthesis
//   natom_b:  The number of atoms in each system of the second synthesis.  If there is no data,
//             it is assumed that the two syntheses are one and the same.
//   xrs:      Source of random numbers to program the swaps
//-------------------------------------------------------------------------------------------------
void resetSwaps(Hybrid<int2> *swaps, const double nsys_a, const double nsys_b,
                Xoroshiro128pGenerator *xrs) {
  const int ns = swaps->size();
  int2* sw_ptr = swaps->data();
  int i = 0;
  const bool a_is_b = (nsys_b < 0.0);
  const double eff_nsys_b = (a_is_b) ? nsys_a : nsys_b;
  while (i < ns) {
    sw_ptr[i] = { static_cast<int>(xrs->uniformRandomNumber() * nsys_a),
                  static_cast<int>(xrs->uniformRandomNumber() * eff_nsys_b) };

    // Check for problematic selections.
    bool x_blocked, y_blocked;
    int iter = 0;
    do {
      x_blocked = false;
      y_blocked = false;
      for (int j = 0; j < i; j++) {
        x_blocked = (x_blocked || sw_ptr[j].x == sw_ptr[i].x);
        y_blocked = (y_blocked || sw_ptr[j].y == sw_ptr[i].y);
        if (a_is_b) {
          x_blocked = (x_blocked || sw_ptr[j].x == sw_ptr[i].y);
          y_blocked = (y_blocked || sw_ptr[j].y == sw_ptr[i].x);
        }
      }
      if (x_blocked) {
        sw_ptr[i].x = static_cast<int>(xrs->uniformRandomNumber() * nsys_a);
      }
      if (y_blocked) {
        sw_ptr[i].y = static_cast<int>(xrs->uniformRandomNumber() * eff_nsys_b);
      }
      iter++;
    } while (x_blocked || y_blocked && iter < 100);
    if (x_blocked == false && y_blocked == false) {
      i++;
    }
  }
  if (i < ns) {
    swaps->resize(i);
  }
}

void resetSwaps(Hybrid<int2> *swaps, const std::vector<int> &natom_a,
                const std::vector<int> &natom_b, Xoroshiro128pGenerator *xrs) {
  const int nsys_a = natom_a.size();
  const bool a_is_b = (natom_b.size() == 0);
  const std::vector<int>& eff_natom_b = (a_is_b) ? natom_a : natom_b;
  const int nsys_b = eff_natom_b.size();
  const int ns = swaps->size();
  int2* sw_ptr = swaps->data();
  int i = 0;
  int iter = 0;
  std::vector<int> matches(16);
  while (i < ns && iter < 1000) {
    int iter_b = 0;
    do {
      sw_ptr[i].x = xrs->uniformRandomNumber() * static_cast<double>(nsys_a);
      matches.resize(0);

      // Count the number of pair candidates
      for (int j = 0; j < nsys_b; j++) {
        if (natom_a[sw_ptr[i].x] == eff_natom_b[j]) {;
          matches.push_back(j);
        }
      }
      iter_b++;
    } while (matches.size() == 0 && iter_b < 10);

    // Cycle if no matches were found.
    if (matches.size() == 0) {
      iter++;
      continue;
    }
    else {
      const int mj = xrs->uniformRandomNumber() * static_cast<double>(matches.size());
      sw_ptr[i].y = matches[mj];
    }
    
    // Check for problematic selections.  Only advance the counter if a suitable pair has been
    // found, but more trials are available.
    bool x_blocked = false;
    bool y_blocked = false;
    for (int j = 0; j < i; j++) {
      x_blocked = (x_blocked || sw_ptr[j].x == sw_ptr[i].x);
      y_blocked = (y_blocked || sw_ptr[j].y == sw_ptr[i].y);
      if (a_is_b) {
        x_blocked = (x_blocked || sw_ptr[j].x == sw_ptr[i].y);
        y_blocked = (y_blocked || sw_ptr[j].y == sw_ptr[i].x);
      }
    }
    if (x_blocked == false && y_blocked == false) {
      i++;
    }
    iter++;
  }
  if (i < ns) {
    swaps->resize(i);
  }
}

//-------------------------------------------------------------------------------------------------
// Create coordinate series and attempt swaps.
//
// Arguments:
//   tsm:      A collection of test systems to draw from
//   sys_idx:  Index of the system within the collection of test systems
//   na:       The number of replicas of the system in the first series
//   nb:       The number of replicas of the system in the second series
//   nswap:    The number of swaps to attempt
//   xrs:      Source of random numbers to differentiate coordinates and program the swaps
//   gpu:      Specifications of the GPU availabel to perform the swap(s)
//-------------------------------------------------------------------------------------------------
template <typename T>
void csSwapTests(const TestSystemManager &tsm, const int sys_idx, const int na, const int nb,
                 const int nswap, Xoroshiro128pGenerator *xrs, const GpuDetails &gpu = null_gpu) {
  CoordinateSeries<T> cs_a(tsm.exportCoordinateFrame(sys_idx), 16);
  CoordinateSeries<T> cs_b(tsm.exportCoordinateFrame(sys_idx), 16);
  differentiateCoordinates(&cs_a, xrs);
  differentiateCoordinates(&cs_b, xrs);
  Hybrid<int2> cs_swaps(nswap, "cs_swaps");
#ifdef STORMM_USE_HPC
  const std::vector<HybridTargetLevel> levels = { HybridTargetLevel::HOST,
                                                  HybridTargetLevel::DEVICE };
#else
  const std::vector<HybridTargetLevel> levels = { HybridTargetLevel::HOST };
#endif
  for (size_t i = 0; i < levels.size(); i++) {
    for (size_t j = 0; j < levels.size(); j++) {
      resetSwaps(&cs_swaps, cs_a.getFrameCount(), cs_b.getFrameCount(), xrs);
      checkSwap(&cs_a, &cs_b, cs_swaps, tsm.getTestingStatus(), levels[i], levels[j], gpu);
      if (na >= nb) {
        resetSwaps(&cs_swaps, cs_a.getFrameCount(), -1, xrs);
        checkSwap(&cs_a, &cs_a, cs_swaps, tsm.getTestingStatus());
      }
      else {
        resetSwaps(&cs_swaps, cs_b.getFrameCount(), -1, xrs);
        checkSwap(&cs_b, &cs_b, cs_swaps, tsm.getTestingStatus());
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Create PhaseSpace syntheses and attempt swaps.  Descriptions of imput arguments follow from
// csSwapTests(), above, in addition to:
//
// Arguments:
//   fp_bits:   The number of fixed precision bits with which to create the coordinate synthesis
//-------------------------------------------------------------------------------------------------
void psSwapTests(const TestSystemManager &tsm, const int na, const int nb, const int nswap,
                 Xoroshiro128pGenerator *xrs, const GpuDetails &gpu = null_gpu,
                 const int fp_bits = 48) {
  const double dnsys = tsm.getSystemCount();
  std::vector<int> index_a = incrementingSeries(0, na);
  std::vector<int> index_b = incrementingSeries(0, nb);
  for (int i = tsm.getSystemCount(); i < na; i++) {
    index_a[i] = dnsys * xrs->uniformRandomNumber();
  }
  for (int i = tsm.getSystemCount(); i < nb; i++) {
    index_b[i] = dnsys * xrs->uniformRandomNumber();
  }
  PhaseSpaceSynthesis poly_psa = tsm.exportPhaseSpaceSynthesis(index_a, 0.0, 1980764, fp_bits);
  PhaseSpaceSynthesis poly_psb = tsm.exportPhaseSpaceSynthesis(index_b, 0.0, 7710835, fp_bits);
  differentiateCoordinates(&poly_psa, xrs);
  differentiateCoordinates(&poly_psb, xrs);
  std::vector<int> atom_counts_a(na), atom_counts_b(nb);
  for (int i = 0; i < na; i++) {
    atom_counts_a[i] = poly_psa.getAtomCount(i);
  }
  for (int i = 0; i < nb; i++) {
    atom_counts_b[i] = poly_psb.getAtomCount(i);
  }
  Hybrid<int2> poly_ps_swaps(nswap, "poly_ps_swaps");
#ifdef STORMM_USE_HPC
  const std::vector<HybridTargetLevel> levels = { HybridTargetLevel::HOST,
                                                  HybridTargetLevel::DEVICE };
#else
  const std::vector<HybridTargetLevel> levels = { HybridTargetLevel::HOST };
#endif
  for (size_t i = 0; i < levels.size(); i++) {
    for (size_t j = 0; j < levels.size(); j++) {
      resetSwaps(&poly_ps_swaps, atom_counts_a, atom_counts_b, xrs);
      checkSwap(&poly_psa, &poly_psb, poly_ps_swaps, tsm.getTestingStatus(), levels[i], levels[j],
                gpu);
      if (na > nb) {
        resetSwaps(&poly_ps_swaps, atom_counts_a, std::vector<int>(), xrs);
        checkSwap(&poly_psa, &poly_psa, poly_ps_swaps, tsm.getTestingStatus(), levels[i],
                  levels[j], gpu);
      }
      else {
        resetSwaps(&poly_ps_swaps, atom_counts_b, std::vector<int>(), xrs);
        checkSwap(&poly_psb, &poly_psb, poly_ps_swaps, tsm.getTestingStatus(), levels[i],
                  levels[j], gpu);
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Create Condensate objects and attempt swaps.  Descriptions of input arguments follow from
// psSwapTests(), above.
//-------------------------------------------------------------------------------------------------
void cdnsSwapTests(const TestSystemManager &tsm, const int na, const int nb, const int nswap,
                   Xoroshiro128pGenerator *xrs, const GpuDetails &gpu = null_gpu) {
  const double dnsys = tsm.getSystemCount();
  std::vector<int> index_a = incrementingSeries(0, na);
  std::vector<int> index_b = incrementingSeries(0, nb);
  for (int i = tsm.getSystemCount(); i < na; i++) {
    index_a[i] = dnsys * xrs->uniformRandomNumber();
  }
  for (int i = tsm.getSystemCount(); i < nb; i++) {
    index_b[i] = dnsys * xrs->uniformRandomNumber();
  }
  PhaseSpaceSynthesis poly_psa = tsm.exportPhaseSpaceSynthesis(index_a, 0.0, 1980764, 30);
  PhaseSpaceSynthesis poly_psb = tsm.exportPhaseSpaceSynthesis(index_b, 0.0, 7710835, 30);
  Condensate cdns_a(poly_psa, PrecisionModel::DOUBLE);
  Condensate cdns_b(poly_psb, PrecisionModel::DOUBLE);
  differentiateCoordinates(&cdns_a, xrs);
  differentiateCoordinates(&cdns_b, xrs);
  std::vector<int> atom_counts_a(na), atom_counts_b(nb);
  for (int i = 0; i < na; i++) {
    atom_counts_a[i] = cdns_a.getAtomCount(i);
  }
  for (int i = 0; i < nb; i++) {
    atom_counts_b[i] = cdns_b.getAtomCount(i);
  }
  Hybrid<int2> cdns_swaps(nswap, "poly_ps_swaps");
#ifdef STORMM_USE_HPC
  const std::vector<HybridTargetLevel> levels = { HybridTargetLevel::HOST,
                                                  HybridTargetLevel::DEVICE };
#else
  const std::vector<HybridTargetLevel> levels = { HybridTargetLevel::HOST };
#endif
  for (size_t i = 0; i < levels.size(); i++) {
    for (size_t j = 0; j < levels.size(); j++) {
      resetSwaps(&cdns_swaps, atom_counts_a, atom_counts_b, xrs);
      checkSwap(&cdns_a, &cdns_b, cdns_swaps, tsm.getTestingStatus(), levels[i], levels[j], gpu);
      if (na >= nb) {
        resetSwaps(&cdns_swaps, atom_counts_a, std::vector<int>(), xrs);
        checkSwap(&cdns_a, &cdns_a, cdns_swaps, tsm.getTestingStatus(), levels[i], levels[j], gpu);
      }
      else {
        resetSwaps(&cdns_swaps, atom_counts_b, std::vector<int>(), xrs);
        checkSwap(&cdns_b, &cdns_b, cdns_swaps, tsm.getTestingStatus(), levels[i], levels[j], gpu);
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Test input traps associated with CoordinateSeries swaps.
//
// Arguments:
//   tsm:         Collection of test ssytems from which to draw examples
//   system_idx:  Index of the system to work with
//-------------------------------------------------------------------------------------------------
void testBadCoordinateSeriesInputs(const TestSystemManager &tsm) {
  const std::vector<int> indices = incrementingSeries(0, tsm.getSystemCount());
  PhaseSpaceSynthesis poly_psa = tsm.exportPhaseSpaceSynthesis(indices);
  PhaseSpaceSynthesis poly_psb = tsm.exportPhaseSpaceSynthesis(indices);
  Hybrid<int2> swaps(4, "bad_swaps");
  int scon = 0;
  for (int i = 0; i < poly_psa.getSystemCount(); i++) {
    for (int j = 0; j < poly_psb.getSystemCount(); j++) {
      if (poly_psa.getAtomCount(i) != poly_psb.getAtomCount(j) && scon == 0) {
        swaps.putHost({ i, j }, scon);
      }
    }
  }
  swaps.putHost({ 0, 1 }, 1);
  swaps.putHost({ 1, 2 }, 2);
  swaps.putHost({ 2, 0 }, 3);
  CHECK_THROWS_SOFT(coordSwap(&poly_psa, &poly_psb, swaps), "An impossible swap between systems "
                    "of different atom counts was attempted.", tsm.getTestingStatus());
  PhaseSpaceSynthesis poly_psc = tsm.exportPhaseSpaceSynthesis(indices, 0.0, 67267483, 39);
  swaps.putHost({ 1, 1 }, 1);
  swaps.putHost({ 2, 2 }, 2);
  swaps.putHost({ 3, 3 }, 3);
  CHECK_THROWS_SOFT(coordSwap(&poly_psa, &poly_psc, swaps), "A swap between two "
                    "PhaseSpaceSyntheses of different bit counts was carried out.",
                    tsm.getTestingStatus());
  CoordinateSeries<float> csa(tsm.exportCoordinateFrame(4), 12);
  CoordinateSeries<float> csb(tsm.exportCoordinateFrame(4), 12);
  swaps.putHost({ 2, 3 }, 2);
  CHECK_THROWS_SOFT(coordSwap(&csa, &csb, swaps), "An impossible swap involving the same system "
                    "in multiple positions was attempted.", tsm.getTestingStatus());
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
#else
  const GpuDetails gpu = null_gpu;
#endif

  // Section 1
  section("Swaps in CoordinateSeries objects");

  // Section 2
  section("Swaps in PhaseSpaceSynthesis objects");

  // Section 3
  section("Swaps in Condensate objects");
  
  // Section 4
  section("Test input traps");

  // Read in systems.  Make perturbations and, if HPC mode is engaged, upload the coordinates.
  const std::vector<std::string> iso_systems = { "drug_example_iso", "med_1", "med_2", "med_3",
                                                 "med_4", "med_5", "symmetry_C1", "symmetry_C2",
                                                 "bromobenzene_vs_iso", "symmetry_L1_vs" };
  const std::vector<std::string> pbc_systems = { "bromobenzene", "drug_example_vs", "tamavidin",
                                                 "symmetry_C1_in_water", "symmetry_C6_in_water",
                                                 "tip3p", "tip4p" };
  const char osc = osSeparator();
  const std::string top_dir = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string crd_dir = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  TestSystemManager tsm_iso(top_dir, "top", iso_systems, crd_dir, "inpcrd", iso_systems);
  TestSystemManager tsm_pbc(top_dir, "top", pbc_systems, crd_dir, "inpcrd", pbc_systems);
  Xoroshiro128pGenerator xrs;
  
  // Create various coordinate series and test swaps
  section(1);
  csSwapTests<double>(tsm_iso, 3, 16, 16, 4, &xrs, gpu);
  csSwapTests<float>(tsm_iso, 7, 16, 16, 4, &xrs, gpu);
  csSwapTests<llint>(tsm_iso, 8, 16, 12, 4, &xrs, gpu);
  csSwapTests<double>(tsm_iso, 9, 16, 16, 4, &xrs, gpu);
  csSwapTests<float>(tsm_iso, 0, 16, 16, 1, &xrs, gpu);
  csSwapTests<llint>(tsm_iso, 5, 16, 16, 1, &xrs, gpu);

  // Create various coordinate syntheses and test swaps
  section(2);
  psSwapTests(tsm_iso, 48, 24, 1, &xrs, gpu);
  psSwapTests(tsm_pbc, 24, 72, 1, &xrs, gpu);
  psSwapTests(tsm_iso, 40, 40, 8, &xrs, gpu);
  psSwapTests(tsm_pbc, 40, 20, 8, &xrs, gpu);
  psSwapTests(tsm_iso, 36, 18, 1, &xrs, gpu, 30);
  psSwapTests(tsm_pbc, 20, 60, 1, &xrs, gpu, 30);
  psSwapTests(tsm_iso, 40, 40, 8, &xrs, gpu, 30);
  psSwapTests(tsm_pbc, 20, 60, 8, &xrs, gpu, 30);
  section(3);
  cdnsSwapTests(tsm_iso, 30, 36,  1, &xrs, gpu);
  cdnsSwapTests(tsm_pbc, 48, 30,  1, &xrs, gpu);
  cdnsSwapTests(tsm_iso, 30, 36,  8, &xrs, gpu);
  cdnsSwapTests(tsm_pbc, 48, 30, 10, &xrs, gpu);

  // Test bad inputs
  section(4);
  testBadCoordinateSeriesInputs(tsm_iso);
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
