#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/coordinate_series.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::symbols::pi;

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Extract data from a coordinate object, whether on the CPU host or GPU device
//
// Overloaded:
//   - Extract data from a CoordinateFrame
//   - Extract data from a PhaseSpace object
//   - Extract data from one frame of a CoordinateSeries
//   - Extract data from one system of a PhaseSpaceSynthesis
//
// Arguments:
//   v:     The coordinate object to extract data from
//   tier:  Obtain coordinate data from the CPU host or the GPU device
//   x:     Cartesian X coordinates of particles, modified and returned
//   y:     Cartesian Y coordinates of particles, modified and returned
//   z:     Cartesian Z coordinates of particles, modified and returned
//   umat:  Fractional space transformation to extract from the object, modified and returned
//   invu:  Inverse transformation to extract from the object, modified and returned
//   bdim:  Box dimensions to extract from the object, modified and returned
//-------------------------------------------------------------------------------------------------
void extractCoordinates(const CoordinateFrame &v, const HybridTargetLevel tier,
                        std::vector<double> *x, std::vector<double> *y, std::vector<double> *z,
                        std::vector<double> *umat, std::vector<double> *invu,
                        std::vector<double> *bdim) {
  const size_t natom = v.getAtomCount();
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (v.getFormat()) {
    case HybridFormat::HOST_ONLY:
#ifdef STORMM_USE_HPC
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_MOUNTED:
#endif
      {
        x->resize(natom);
        y->resize(natom);
        z->resize(natom);
        umat->resize(9);
        invu->resize(9);
        bdim->resize(6);
        double* x_ptr = x->data();
        double* y_ptr = y->data();
        double* z_ptr = z->data();
        double* umat_ptr = umat->data();
        double* invu_ptr = invu->data();
        double* bdim_ptr = bdim->data();
        const CoordinateFrameReader vr(v.data(tier));
        for (size_t i = 0; i < natom; i++) {
          x_ptr[i] = vr.xcrd[i];
          y_ptr[i] = vr.ycrd[i];
          z_ptr[i] = vr.zcrd[i];
        }
        for (size_t i = 0; i < 9; i++) {
          umat_ptr[i] = vr.umat[i];
          invu_ptr[i] = vr.invu[i];
        }
        for (size_t i = 0; i < 6; i++) {
          bdim_ptr[i] = vr.boxdim[i];
        }
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridFormat::DEVICE_ONLY:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (v.getFormat()) {
    case HybridFormat::HOST_ONLY:
    case HybridFormat::HOST_MOUNTED:
      break;
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
    case HybridFormat::DEVICE_ONLY:
      {
        const Hybrid<double>* vx = v.getCoordinateHandle(CartesianDimension::X);
        const Hybrid<double>* vy = v.getCoordinateHandle(CartesianDimension::Y);
        const Hybrid<double>* vz = v.getCoordinateHandle(CartesianDimension::Z);
        *x = vx->readDevice();
        *y = vy->readDevice();
        *z = vz->readDevice();
        const Hybrid<double>* vumat = v.getBoxTransformHandle();
        const Hybrid<double>* vinvu = v.getInverseTransformHandle();
        const Hybrid<double>* vbdim = v.getBoxDimensionsHandle();
        *umat = vumat->readDevice();
        *invu = vinvu->readDevice();
        *bdim = vbdim->readDevice();
      }
      break;
    }
    break;
#endif
  }
}

void extractCoordinates(const PhaseSpace &v, const TrajectoryKind kind,
                        const CoordinateCycle stage, const HybridTargetLevel tier,
                        std::vector<double> *x, std::vector<double> *y, std::vector<double> *z,
                        std::vector<double> *umat, std::vector<double> *invu,
                        std::vector<double> *bdim) {
  const size_t natom = v.getAtomCount();
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (v.getFormat()) {
    case HybridFormat::HOST_ONLY:
#ifdef STORMM_USE_HPC
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_MOUNTED:
#endif
      {
        x->resize(natom);
        y->resize(natom);
        z->resize(natom);
        umat->resize(9);
        invu->resize(9);
        bdim->resize(6);
        double* x_ptr = x->data();
        double* y_ptr = y->data();
        double* z_ptr = z->data();
        double* umat_ptr = umat->data();
        double* invu_ptr = invu->data();
        double* bdim_ptr = bdim->data();
        const PhaseSpaceReader vr = v.data(stage, tier);
        switch (kind) {
        case TrajectoryKind::POSITIONS:
          for (size_t i = 0; i < natom; i++) {
            x_ptr[i] = vr.xcrd[i];
            y_ptr[i] = vr.ycrd[i];
            z_ptr[i] = vr.zcrd[i];
          }
          for (size_t i = 0; i < 9; i++) {
            umat_ptr[i] = vr.umat[i];
            invu_ptr[i] = vr.invu[i];
          }
          for (size_t i = 0; i < 6; i++) {
            bdim_ptr[i] = vr.boxdim[i];
          }
          break;
        case TrajectoryKind::VELOCITIES:
          for (size_t i = 0; i < natom; i++) {
            x_ptr[i] = vr.xvel[i];
            y_ptr[i] = vr.yvel[i];
            z_ptr[i] = vr.zvel[i];
          }
          break;
        case TrajectoryKind::FORCES:
          for (size_t i = 0; i < natom; i++) {
            x_ptr[i] = vr.xfrc[i];
            y_ptr[i] = vr.yfrc[i];
            z_ptr[i] = vr.zfrc[i];
          }
          break;
        }
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridFormat::DEVICE_ONLY:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (v.getFormat()) {
    case HybridFormat::HOST_ONLY:
    case HybridFormat::HOST_MOUNTED:
      break;
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
    case HybridFormat::DEVICE_ONLY:
      {
        const Hybrid<double>* vx = v.getCoordinateHandle(CartesianDimension::X, kind, stage);
        const Hybrid<double>* vy = v.getCoordinateHandle(CartesianDimension::Y, kind, stage);
        const Hybrid<double>* vz = v.getCoordinateHandle(CartesianDimension::Z, kind, stage);
        *x = vx->readDevice();
        *y = vy->readDevice();
        *z = vz->readDevice();
        switch (kind) {
        case TrajectoryKind::POSITIONS:
          *umat = v.getBoxTransformHandle()->readDevice();
          *invu = v.getInverseTransformHandle()->readDevice();
          *bdim = v.getBoxDimensionsHandle()->readDevice();
          break;
        case TrajectoryKind::VELOCITIES:
        case TrajectoryKind::FORCES:
          break;
        }
      }
      break;
    }
    break;
#endif
  }
}

template <typename T>
void extractCoordinates(const CoordinateSeries<T> &v, const size_t frame_idx,
                        const HybridTargetLevel tier, std::vector<double> *x,
                        std::vector<double> *y, std::vector<double> *z, std::vector<double> *umat,
                        std::vector<double> *invu, std::vector<double> *bdim) {
  const size_t natom = v.getAtomCount();
  const size_t padded_natom = roundUp(natom, warp_size_zu);
  const size_t padded_xfrm = roundUp(9, warp_size_int);
  const size_t padded_bdim = roundUp(6, warp_size_int);
  const double descaling_factor = (isSignedIntegralScalarType<T>()) ?
                                  pow(2.0, -(v.getFixedPrecisionBits())) : 1.0;
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (v.getFormat()) {
    case HybridFormat::HOST_ONLY:
#ifdef STORMM_USE_HPC
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_MOUNTED:
#endif
      {
        x->resize(natom);
        y->resize(natom);
        z->resize(natom);
        umat->resize(9);
        invu->resize(9);
        bdim->resize(6);
        double* x_ptr = x->data();
        double* y_ptr = y->data();
        double* z_ptr = z->data();
        double* umat_ptr = umat->data();
        double* invu_ptr = invu->data();
        double* bdim_ptr = bdim->data();
        const CoordinateSeriesReader<T> vr(v.data(tier));
        for (size_t i = 0; i < natom; i++) {
          x_ptr[i] = vr.xcrd[(padded_natom * frame_idx) + i] * descaling_factor;
          y_ptr[i] = vr.ycrd[(padded_natom * frame_idx) + i] * descaling_factor;
          z_ptr[i] = vr.zcrd[(padded_natom * frame_idx) + i] * descaling_factor;
        }
        for (size_t i = 0; i < 9; i++) {
          umat_ptr[i] = vr.umat[(padded_xfrm * frame_idx) + i];
          invu_ptr[i] = vr.invu[(padded_xfrm * frame_idx) + i];
        }
        for (size_t i = 0; i < 6; i++) {
          bdim_ptr[i] = vr.boxdim[(padded_bdim * frame_idx) + i];
        }
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridFormat::DEVICE_ONLY:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (v.getFormat()) {
    case HybridFormat::HOST_ONLY:
    case HybridFormat::HOST_MOUNTED:
      break;
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
    case HybridFormat::DEVICE_ONLY:
      {
        const Hybrid<T>* vx = v.getFramesHandle(CartesianDimension::X);
        const Hybrid<T>* vy = v.getFramesHandle(CartesianDimension::Y);
        const Hybrid<T>* vz = v.getFramesHandle(CartesianDimension::Z);
        const std::vector<T> x_tmp = vx->readDevice();
        const std::vector<T> y_tmp = vy->readDevice();
        const std::vector<T> z_tmp = vz->readDevice();
        double* x_ptr = x->data();
        double* y_ptr = y->data();
        double* z_ptr = z->data();
        for (size_t i = 0; i < natom; i++) {
          x_ptr[i] = static_cast<double>(x_tmp[(padded_natom * frame_idx) + i]) * descaling_factor;
          y_ptr[i] = static_cast<double>(y_tmp[(padded_natom * frame_idx) + i]) * descaling_factor;
          z_ptr[i] = static_cast<double>(z_tmp[(padded_natom * frame_idx) + i]) * descaling_factor;
        }
        const Hybrid<double>* vumat = v.getBoxTransformsHandle();
        const Hybrid<double>* vinvu = v.getInverseTransformsHandle();
        const Hybrid<double>* vbdim = v.getBoxDimensionsHandle();
        *umat = vumat->readDevice(padded_xfrm * frame_idx, 9);
        *invu = vinvu->readDevice(padded_xfrm * frame_idx, 9);
        *bdim = vbdim->readDevice(padded_bdim * frame_idx, 6);
      }
      break;
    }
    break;
#endif
  }
}

void extractCoordinates(const PhaseSpaceSynthesis &v, const size_t system_idx,
                        const TrajectoryKind kind, const CoordinateCycle stage,
                        const HybridTargetLevel tier, std::vector<double> *x,
                        std::vector<double> *y, std::vector<double> *z, std::vector<double> *umat,
                        std::vector<double> *invu, std::vector<double> *bdim) {
  const size_t natom = v.getAtomCount(system_idx);
  const size_t offset = v.getAtomOffset(system_idx);
  double descaling_factor;
  x->resize(natom);
  y->resize(natom);
  z->resize(natom);
  umat->resize(9);
  invu->resize(9);
  bdim->resize(6);
  double* x_ptr = x->data();
  double* y_ptr = y->data();
  double* z_ptr = z->data();
  double* umat_ptr = umat->data();
  double* invu_ptr = invu->data();
  double* bdim_ptr = bdim->data();
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    descaling_factor = pow(2.0, -(v.getGlobalPositionBits()));
    break;
  case TrajectoryKind::VELOCITIES:
    descaling_factor = pow(2.0, -(v.getVelocityBits()));
    break;
  case TrajectoryKind::FORCES:
    descaling_factor = pow(2.0, -(v.getForceAccumulationBits()));
    break;
  }
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (v.getFormat()) {
    case HybridFormat::HOST_ONLY:
#ifdef STORMM_USE_HPC
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_MOUNTED:
#endif
      {
        const PsSynthesisReader vr = v.data(stage, tier);
        switch (kind) {
        case TrajectoryKind::POSITIONS:
          for (size_t i = 0; i < natom; i++) {
            x_ptr[i] = hostInt95ToDouble(vr.xcrd[i + offset], vr.xcrd_ovrf[i + offset]);
            y_ptr[i] = hostInt95ToDouble(vr.ycrd[i + offset], vr.ycrd_ovrf[i + offset]);
            z_ptr[i] = hostInt95ToDouble(vr.zcrd[i + offset], vr.zcrd_ovrf[i + offset]);
          }
          for (size_t i = 0; i < 9; i++) {
            umat_ptr[i] = vr.umat[i];
            invu_ptr[i] = vr.invu[i];
          }
          for (size_t i = 0; i < 6; i++) {
            bdim_ptr[i] = vr.boxdims[i];
          }
          break;
        case TrajectoryKind::VELOCITIES:
          for (size_t i = 0; i < natom; i++) {
            x_ptr[i] = hostInt95ToDouble(vr.xvel[i + offset], vr.xvel_ovrf[i + offset]);
            y_ptr[i] = hostInt95ToDouble(vr.yvel[i + offset], vr.yvel_ovrf[i + offset]);
            z_ptr[i] = hostInt95ToDouble(vr.zvel[i + offset], vr.zvel_ovrf[i + offset]);
          }
          break;
        case TrajectoryKind::FORCES:
          for (size_t i = 0; i < natom; i++) {
            x_ptr[i] = hostInt95ToDouble(vr.xfrc[i + offset], vr.xfrc_ovrf[i + offset]);
            y_ptr[i] = hostInt95ToDouble(vr.yfrc[i + offset], vr.yfrc_ovrf[i + offset]);
            z_ptr[i] = hostInt95ToDouble(vr.zfrc[i + offset], vr.zfrc_ovrf[i + offset]);
          }
          break;
        }
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridFormat::DEVICE_ONLY:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (v.getFormat()) {
    case HybridFormat::HOST_ONLY:
    case HybridFormat::HOST_MOUNTED:
      break;
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
    case HybridFormat::DEVICE_ONLY:
      {
        const Hybrid<llint>* vx = v.getCoordinateHandle(CartesianDimension::X, kind, stage);
        const Hybrid<llint>* vy = v.getCoordinateHandle(CartesianDimension::Y, kind, stage);
        const Hybrid<llint>* vz = v.getCoordinateHandle(CartesianDimension::Z, kind, stage);
        const Hybrid<int>* vx_ovrf = v.getCoordinateOverflowHandle(CartesianDimension::X, kind,
                                                                   stage);
        const Hybrid<int>* vy_ovrf = v.getCoordinateOverflowHandle(CartesianDimension::Y, kind,
                                                                   stage);
        const Hybrid<int>* vz_ovrf = v.getCoordinateOverflowHandle(CartesianDimension::Z, kind,
                                                                   stage);

        // Set up a three-vector conversion with the host-bound split fixed-precision converter,
        // but save the descaling for the final segment of the function, below.
        const std::vector<llint> xprim_tmp = vx->readDevice(offset, natom);
        const std::vector<llint> yprim_tmp = vy->readDevice(offset, natom);
        const std::vector<llint> zprim_tmp = vz->readDevice(offset, natom);
        const std::vector<int> xovrf_tmp = vx_ovrf->readDevice(offset, natom);
        const std::vector<int> yovrf_tmp = vy_ovrf->readDevice(offset, natom);
        const std::vector<int> zovrf_tmp = vz_ovrf->readDevice(offset, natom);
        hostInt95ToDouble(x, y, z, xprim_tmp, xovrf_tmp, yprim_tmp, yovrf_tmp, zprim_tmp,
                          zovrf_tmp);
        switch (kind) {
        case TrajectoryKind::POSITIONS:
          *umat = v.getBoxTransformsHandle()->readDevice();
          *invu = v.getInverseTransformsHandle()->readDevice();
          *bdim = v.getBoxDimensionsHandle()->readDevice();
          break;
        case TrajectoryKind::VELOCITIES:
        case TrajectoryKind::FORCES:
          break;
        }
      }
      break;
    }
    break;
#endif
  }
  for (size_t i = 0; i < natom; i++) {
    x_ptr[i] *= descaling_factor;
    y_ptr[i] *= descaling_factor;
    z_ptr[i] *= descaling_factor;
  }
}

//-------------------------------------------------------------------------------------------------
// Develop box dimensions for a system.  Multiple calls to this function can produce different
// unit cells in the CPU host and GPU device memory.
//
// Arguments:
//   xrs:   Source of random numbers to create a triclinic unit cell
//   umat:  Transformation taking coordinates in Cartesian space into unit cell fractional space
//   invu:  The inverse box transform, taking fractional coordinates into real space
//   bdim:  Summary of the box dimensions, modified and returned to be consistent with umat and
//          invu
//-------------------------------------------------------------------------------------------------
void createBoxDimensions(Xoroshiro128pGenerator *xrs, std::vector<double> *umat,
                         std::vector<double> *invu, std::vector<double> *bdim) {
  const double box_x = 45.0 + (10.0 * xrs->uniformRandomNumber());
  const double box_y = 45.0 + (10.0 * xrs->uniformRandomNumber());
  const double box_z = 45.0 + (10.0 * xrs->uniformRandomNumber());
  const double alpha = (0.45 + (0.1 * xrs->uniformRandomNumber())) * pi;
  const double beta  = (0.45 + (0.1 * xrs->uniformRandomNumber())) * pi;
  const double gamma = (0.45 + (0.1 * xrs->uniformRandomNumber())) * pi;
  umat->resize(9);
  invu->resize(9);
  bdim->resize(6);
  bdim->at(0) = box_x;
  bdim->at(1) = box_y;
  bdim->at(2) = box_z;
  bdim->at(3) = alpha;
  bdim->at(4) = beta;
  bdim->at(5) = gamma;
  computeBoxTransform(bdim->data(), umat->data(), invu->data());
}

//-------------------------------------------------------------------------------------------------
// Execute a specific comparison between two coordinate objects.
//
// Overloaded:
//   - Operate on two CoordinateFrame objects
//   - Operate on two PhaseSpace objects
//
// Arguments:
//   a:         The original object
//   b:         The copied object
//   tier_a:    Obtain memory from the original object on the CPU host or GPU device
//   tier_b:    Obtain memory from the copied object on the CPU host or GPU device
//   do_tests:  Indicate whether tests are possible (relevant in the case of PhaseSpaceSynthesis
//              comparisons, when a topology file was critical to creating the synthesis)
//-------------------------------------------------------------------------------------------------
void makeCoordinateComparison(const CoordinateFrame &a, const CoordinateFrame &b,
                              const HybridTargetLevel tier_a, const HybridTargetLevel tier_b,
                              const TestPriority do_tests = TestPriority::CRITICAL) {

  // Lay out vectors for any scenario.  The extractCoordinates() function will fill them as
  // appropriate.
  const size_t atom_count = a.getAtomCount();
  std::vector<double> xa_snap(atom_count), ya_snap(atom_count), za_snap(atom_count);
  std::vector<double> xb_snap(atom_count), yb_snap(atom_count), zb_snap(atom_count);
  std::vector<double> a_umat(9), b_umat(9), a_invu(9), b_invu(9), a_bdim(6), b_bdim(6);
  extractCoordinates(a, tier_a, &xa_snap, &ya_snap, &za_snap, &a_umat, &a_invu, &a_bdim);
  extractCoordinates(b, tier_b, &xb_snap, &yb_snap, &zb_snap, &b_umat, &b_invu, &b_bdim);
  check(ya_snap, RelationalOperator::EQUAL, yb_snap, "Cartesian Y coordinates were not properly "
        "copied between CoordinateFrame objects of formats " + getEnumerationName(a.getFormat()) +
        " and " + getEnumerationName(b.getFormat()) + ". (" + getEnumerationName(tier_a) + " -> " +
        getEnumerationName(tier_b) + ".)");
  check(za_snap, RelationalOperator::EQUAL, zb_snap, "Cartesian Z coordinates were not properly "
        "copied between CoordinateFrame objects of formats " + getEnumerationName(a.getFormat()) +
        " and " + getEnumerationName(b.getFormat()) + ". (" + getEnumerationName(tier_a) + " -> " +
        getEnumerationName(tier_b) + ".)");
  std::vector<double> xdiff(atom_count);
  for (int i = 0; i < atom_count; i++) {
    xdiff[i] = xb_snap[i] - xa_snap[i];
  }
  check(mean(xdiff), RelationalOperator::EQUAL, 1.05, "Cartesian X coordinates were not shifted "
        "as expected after copying a CoordinateFrame object.");
}

void makeCoordinateComparison(const PhaseSpace &a, const PhaseSpace &b,
                              const HybridTargetLevel tier_a, const HybridTargetLevel tier_b,
                              const TestPriority do_tests = TestPriority::CRITICAL) {

  // Lay out vectors for any scenario.  The extractCoordinates() function will fill them as
  // appropriate.
  const size_t atom_count = a.getAtomCount();
  std::vector<double> xa_snap(atom_count), ya_snap(atom_count), za_snap(atom_count);
  std::vector<double> xb_snap(atom_count), yb_snap(atom_count), zb_snap(atom_count);
  std::vector<double> a_umat(9), b_umat(9), a_invu(9), b_invu(9), a_bdim(6), b_bdim(6);
  const std::vector<TrajectoryKind> rvf = { TrajectoryKind::POSITIONS, TrajectoryKind::VELOCITIES,
                                            TrajectoryKind::FORCES };
  const std::vector<CoordinateCycle> tps = { CoordinateCycle::WHITE, CoordinateCycle::BLACK }; 
  for (size_t i = 0; i < rvf.size(); i++) {
    for (size_t j = 0; j < tps.size(); j++) {
      extractCoordinates(a, rvf[i], tps[j], tier_a, &xa_snap, &ya_snap, &za_snap, &a_umat,
                         &a_invu, &a_bdim);
      extractCoordinates(b, rvf[i], tps[j], tier_b, &xb_snap, &yb_snap, &zb_snap, &b_umat,
                         &b_invu, &b_bdim);
      check(ya_snap, RelationalOperator::EQUAL, yb_snap, "Cartesian Y coordinates were not "
            "properly copied between PhaseSpace objects of formats " +
            getEnumerationName(a.getFormat()) + " and " + getEnumerationName(b.getFormat()) +
            ". (" + getEnumerationName(tier_a) + " -> " + getEnumerationName(tier_b) + ".)");
      check(za_snap, RelationalOperator::EQUAL, zb_snap, "Cartesian Z coordinates were not "
            "properly copied between PhaseSpace objects of formats " +
            getEnumerationName(a.getFormat()) + " and " + getEnumerationName(b.getFormat()) +
            ". (" + getEnumerationName(tier_a) + " -> " + getEnumerationName(tier_b) + ".)");
      std::vector<double> xdiff(atom_count);
      for (int k = 0; k < atom_count; k++) {
        xdiff[k] = xb_snap[k] - xa_snap[k];
      }
      if (rvf[i] == TrajectoryKind::VELOCITIES && tps[j] == CoordinateCycle::BLACK) {
        check(mean(xdiff), RelationalOperator::EQUAL, 1.05, "Cartesian X velocities were not "
              "modified as expected after copying a PhaseSpace object.");
      }
      else {
        check(mean(xdiff), RelationalOperator::EQUAL, 0.0, "Cartesian X coordinates were not "
              "left unchanged as expected after copying a PhaseSpace object.");
      }
    }
  }
}

template <typename T>
void makeCoordinateComparison(const CoordinateSeries<T> &a, const CoordinateSeries<T> &b,
                              const HybridTargetLevel tier_a, const HybridTargetLevel tier_b,
                              const TestPriority do_tests = TestPriority::CRITICAL) {

  // Lay out vectors for any scenario.  The extractCoordinates() function will fill them as
  // appropriate.
  const size_t atom_count = a.getAtomCount();
  const size_t frame_count = a.getFrameCount();
  std::vector<double> xa_snap(atom_count), ya_snap(atom_count), za_snap(atom_count);
  std::vector<double> xb_snap(atom_count), yb_snap(atom_count), zb_snap(atom_count);
  std::vector<double> a_umat(9), b_umat(9), a_invu(9), b_invu(9), a_bdim(6), b_bdim(6);
  extractCoordinates(a, frame_count / 2, tier_a, &xa_snap, &ya_snap, &za_snap, &a_umat, &a_invu,
                     &a_bdim);
  extractCoordinates(b, frame_count / 2, tier_b, &xb_snap, &yb_snap, &zb_snap, &b_umat, &b_invu,
                     &b_bdim);
  check(ya_snap, RelationalOperator::EQUAL, yb_snap, "Cartesian Y coordinates were not "
          "properly copied between CoordinateSeries objects of formats " +
          getEnumerationName(a.getFormat()) + " and " + getEnumerationName(b.getFormat()) +
          ". (" + getEnumerationName(tier_a) + " -> " + getEnumerationName(tier_b) + ".)");
  check(za_snap, RelationalOperator::EQUAL, zb_snap, "Cartesian Z coordinates were not "
        "properly copied between CoordinateSeries objects of formats " +
        getEnumerationName(a.getFormat()) + " and " + getEnumerationName(b.getFormat()) +
        ". (" + getEnumerationName(tier_a) + " -> " + getEnumerationName(tier_b) + ".)");
  std::vector<double> xdiff(atom_count);
  for (int i = 0; i < atom_count; i++) {
    xdiff[i] = xb_snap[i] - xa_snap[i];
  }

  // Tune the tolerance based on the coordinate series data type
  const size_t tc = std::type_index(typeid(T)).hash_code();
  Approx target(1.05);
  double mg = 1.0e-6;
  if (tc == float_type_index) {
    mg = 1.0e-5;
  }
  check(mean(xdiff), RelationalOperator::EQUAL, target.margin(mg), "Cartesian X coordinates were "
        "not shifted as expected after copying a CoordinateSeries object.");
}

void makeCoordinateComparison(const PhaseSpaceSynthesis &a, const PhaseSpaceSynthesis &b,
                              const HybridTargetLevel tier_a, const HybridTargetLevel tier_b,
                              const TestPriority do_tests) {

  // Lay out vectors for any scenario.  The extractCoordinates() function will fill them as
  // appropriate.
  const int system_idx = a.getSystemCount() / 2;
  const size_t atom_count = a.getAtomCount(system_idx);
  std::vector<double> xa_snap(atom_count), ya_snap(atom_count), za_snap(atom_count);
  std::vector<double> xb_snap(atom_count), yb_snap(atom_count), zb_snap(atom_count);
  std::vector<double> a_umat(9), b_umat(9), a_invu(9), b_invu(9), a_bdim(6), b_bdim(6);
  const std::vector<TrajectoryKind> rvf = { TrajectoryKind::POSITIONS, TrajectoryKind::VELOCITIES,
                                            TrajectoryKind::FORCES };
  const std::vector<CoordinateCycle> tps = { CoordinateCycle::WHITE, CoordinateCycle::BLACK }; 
  for (size_t i = 0; i < rvf.size(); i++) {
    for (size_t j = 0; j < tps.size(); j++) {
      extractCoordinates(a, system_idx, rvf[i], tps[j], tier_a, &xa_snap, &ya_snap, &za_snap,
                         &a_umat, &a_invu, &a_bdim);
      extractCoordinates(b, system_idx, rvf[i], tps[j], tier_b, &xb_snap, &yb_snap, &zb_snap,
                         &b_umat, &b_invu, &b_bdim);
      check(ya_snap, RelationalOperator::EQUAL, yb_snap, "Cartesian Y coordinates were not "
            "properly copied between PhaseSpaceSynthesis objects of formats " +
            getEnumerationName(a.getFormat()) + " and " + getEnumerationName(b.getFormat()) +
            ". (" + getEnumerationName(tier_a) + " -> " + getEnumerationName(tier_b) + ".)");
      check(za_snap, RelationalOperator::EQUAL, zb_snap, "Cartesian Z coordinates were not "
            "properly copied between PhaseSpaceSynthesis objects of formats " +
            getEnumerationName(a.getFormat()) + " and " + getEnumerationName(b.getFormat()) +
            ". (" + getEnumerationName(tier_a) + " -> " + getEnumerationName(tier_b) + ".)");
      std::vector<double> xdiff(atom_count);
      for (int k = 0; k < atom_count; k++) {
        xdiff[k] = xb_snap[k] - xa_snap[k];
      }
      if (rvf[i] == TrajectoryKind::FORCES && tps[j] == CoordinateCycle::WHITE) {
        check(mean(xdiff), RelationalOperator::EQUAL, Approx(1.05).margin(1.0e-6), "Cartesian X "
              "forces were not shifted as expected after copying a PhaseSpaceSynthesis object.",
              do_tests);
      }
      else {
        check(mean(xdiff), RelationalOperator::EQUAL, 0.0, "Cartesian X coordinates were not "
              "left unchanged as expected after copying a PhaseSpaceSynthesis object.", do_tests);
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Run the appropriate comparisons for two coordinate objects.  Different formats will require (and
// allow) particular comparisons across (or perhaps between) the CPU host and GPU device data.
//
// Arguments:
//   crd_a:  The first coordinate object
//   crd_b:  The second coordinate object
//-------------------------------------------------------------------------------------------------
template <typename T>
void runLegalComparisons(const T &crd_a, const T &crd_b,
                         const TestPriority do_tests = TestPriority::CRITICAL) {
  const HybridFormat fmt_a = crd_a.getFormat();
  const HybridFormat fmt_b = crd_b.getFormat();
  switch (fmt_a) {
  case HybridFormat::HOST_ONLY:
    switch (fmt_b) {
    case HybridFormat::HOST_ONLY:
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::HOST, HybridTargetLevel::HOST,
                               do_tests);
      break;
#ifdef STORMM_USE_HPC
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::HOST, HybridTargetLevel::HOST,
                               do_tests);
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::HOST, HybridTargetLevel::DEVICE,
                               do_tests);
      break;
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_MOUNTED:
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::HOST, HybridTargetLevel::HOST,
                               do_tests);
      break;
    case HybridFormat::DEVICE_ONLY:
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::HOST, HybridTargetLevel::DEVICE,
                               do_tests);
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    switch (fmt_b) {
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_MOUNTED:
    case HybridFormat::HOST_ONLY:
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::HOST, HybridTargetLevel::HOST,
                               do_tests);
      break;
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::HOST, HybridTargetLevel::HOST,
                               do_tests);
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::DEVICE, HybridTargetLevel::DEVICE,
                               do_tests);
      break;
    case HybridFormat::DEVICE_ONLY:
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::DEVICE, HybridTargetLevel::DEVICE,
                               do_tests);
      break;
    }
    break;
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
    switch (fmt_b) {
    case HybridFormat::HOST_ONLY:
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::HOST, HybridTargetLevel::HOST,
                               do_tests);
      break;
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::HOST, HybridTargetLevel::HOST,
                               do_tests);
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::HOST, HybridTargetLevel::DEVICE,
                               do_tests);
      break;
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_MOUNTED:
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::HOST, HybridTargetLevel::HOST,
                               do_tests);
      break;
    case HybridFormat::DEVICE_ONLY:
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::HOST, HybridTargetLevel::DEVICE,
                               do_tests);
      break;
    }
    break;
  case HybridFormat::DEVICE_ONLY:
    switch (fmt_b) {
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_MOUNTED:
    case HybridFormat::HOST_ONLY:
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::DEVICE, HybridTargetLevel::HOST,
                               do_tests);
      break;
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::DEVICE, HybridTargetLevel::HOST,
                               do_tests);
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::DEVICE, HybridTargetLevel::DEVICE,
                               do_tests);
      break;
    case HybridFormat::DEVICE_ONLY:
      makeCoordinateComparison(crd_a, crd_b, HybridTargetLevel::DEVICE, HybridTargetLevel::DEVICE,
                               do_tests);
      break;
    }
    break;
#endif
  }
}


//-------------------------------------------------------------------------------------------------
// Test the copy constructor for a CoordinateFrame object, with possible memory format changes.
//
// Arguments:
//   fmt_a:       The memory layout in which the original object shall be created
//   fmt_b:       The memory format of the newly copied object
//   atom_count:  The number of atoms to include in the object
//-------------------------------------------------------------------------------------------------
void testCoordinateFrameCopy(Xoroshiro128pGenerator *xrs, const HybridFormat fmt_a,
                             const HybridFormat fmt_b, const int atom_count) {
  CoordinateFrame cf_a(atom_count, UnitCellType::TRICLINIC, fmt_a);
  Hybrid<double>* x_crd = cf_a.getCoordinateHandle(CartesianDimension::X);
  Hybrid<double>* y_crd = cf_a.getCoordinateHandle(CartesianDimension::Y);
  Hybrid<double>* z_crd = cf_a.getCoordinateHandle(CartesianDimension::Z);
  Hybrid<double>* umat_handle = cf_a.getBoxTransformHandle();
  Hybrid<double>* invu_handle = cf_a.getInverseTransformHandle();
  Hybrid<double>* bdim_handle = cf_a.getBoxDimensionsHandle();
  const std::vector<double> hx_loc = uniformRand(xrs, atom_count, 50.0);
  const std::vector<double> hy_loc = uniformRand(xrs, atom_count, 50.0);
  const std::vector<double> hz_loc = uniformRand(xrs, atom_count, 50.0);
  const std::vector<double> dx_loc = uniformRand(xrs, atom_count, 50.0);
  const std::vector<double> dy_loc = uniformRand(xrs, atom_count, 50.0);
  const std::vector<double> dz_loc = uniformRand(xrs, atom_count, 50.0);
  std::vector<double> h_umat, h_invu, h_bdim, d_umat, d_invu, d_bdim;
  createBoxDimensions(xrs, &h_umat, &h_invu, &h_bdim);
  createBoxDimensions(xrs, &d_umat, &d_invu, &d_bdim);
  switch (fmt_a) {
  case HybridFormat::HOST_ONLY:
    x_crd->putHost(hx_loc);
    y_crd->putHost(hy_loc);
    z_crd->putHost(hz_loc);
    umat_handle->putHost(h_umat);
    invu_handle->putHost(h_invu);
    bdim_handle->putHost(h_bdim);
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::UNIFIED:
    x_crd->putHost(hx_loc);
    y_crd->putHost(hy_loc);
    z_crd->putHost(hz_loc);
    umat_handle->putHost(h_umat);
    invu_handle->putHost(h_invu);
    bdim_handle->putHost(h_bdim);
    break;
  case HybridFormat::DEVICE_ONLY:
    x_crd->putDevice(dx_loc);
    y_crd->putDevice(dy_loc);
    z_crd->putDevice(dz_loc);
    umat_handle->putDevice(d_umat);
    invu_handle->putDevice(d_invu);
    bdim_handle->putDevice(d_bdim);
    break;
  case HybridFormat::DECOUPLED:
  case HybridFormat::EXPEDITED:
    x_crd->putHost(hx_loc);
    y_crd->putHost(hy_loc);
    z_crd->putHost(hz_loc);
    x_crd->putDevice(dx_loc);
    y_crd->putDevice(dy_loc);
    z_crd->putDevice(dz_loc);
    umat_handle->putHost(h_umat);
    invu_handle->putHost(h_invu);
    bdim_handle->putHost(h_bdim);
    umat_handle->putDevice(d_umat);
    invu_handle->putDevice(d_invu);
    bdim_handle->putDevice(d_bdim);
    break;
#endif
  }
  CoordinateFrame cf_b(cf_a, fmt_b);
  x_crd = cf_b.getCoordinateHandle(CartesianDimension::X);
  std::vector<double> xh_tmp, xd_tmp;
  switch (fmt_b) {
  case HybridFormat::HOST_ONLY:
    xh_tmp = x_crd->readHost();
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    xh_tmp = x_crd->readHost();
    xd_tmp = x_crd->readDevice();
    break;
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::UNIFIED:
    xh_tmp = x_crd->readHost();
    break;
  case HybridFormat::DEVICE_ONLY:
    xd_tmp = x_crd->readDevice();
    break;
#endif
  }
  const bool has_host = (xh_tmp.size() > 0);
  const bool has_devc = (xd_tmp.size() > 0);
  for (int i = 0; i < atom_count; i++) {
    if (has_host) {
      xh_tmp[i] += 1.05;
    }
    if (has_devc) {
      xd_tmp[i] += 1.05;
    }
  }
  switch (fmt_b) {
  case HybridFormat::HOST_ONLY:
    x_crd->putHost(xh_tmp);
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    x_crd->putHost(xh_tmp);
    x_crd->putDevice(xd_tmp);
    break;
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::UNIFIED:
    x_crd->putHost(xh_tmp);
    break;
  case HybridFormat::DEVICE_ONLY:
    x_crd->putDevice(xd_tmp);
    break;
#endif
  }

  // Depending on the formats of the original and copied objects, run specific tests.
  runLegalComparisons(cf_a, cf_b);
}

//-------------------------------------------------------------------------------------------------
// Fill out a PhaseSpace object with random positions, velocities, and forces.  This encapsulates
// code to avoid redundancy.
//
// Arguments:
//   atom_count:  The number of atoms to place in the PhaseSpace object, which is then returned
//   fmt:         The memory layout in which the object is to be created
//   xrs:         Source of random numbers
//-------------------------------------------------------------------------------------------------
PhaseSpace detailPhaseSpace(const int atom_count, const HybridFormat fmt,
                            Xoroshiro128pGenerator *xrs) {
  PhaseSpace result(atom_count, UnitCellType::TRICLINIC, fmt);
  const std::vector<TrajectoryKind> rvf = { TrajectoryKind::POSITIONS, TrajectoryKind::VELOCITIES,
                                            TrajectoryKind::FORCES };
  const std::vector<CoordinateCycle> tps = { CoordinateCycle::WHITE, CoordinateCycle::BLACK };
  for (size_t i = 0; i < rvf.size(); i++) {
    for (size_t j = 0; j < tps.size(); j++) {
      Hybrid<double>* x_crd = result.getCoordinateHandle(CartesianDimension::X, rvf[i], tps[j]);
      Hybrid<double>* y_crd = result.getCoordinateHandle(CartesianDimension::Y, rvf[i], tps[j]);
      Hybrid<double>* z_crd = result.getCoordinateHandle(CartesianDimension::Z, rvf[i], tps[j]);
      Hybrid<double>* umat_handle = result.getBoxTransformHandle();
      Hybrid<double>* invu_handle = result.getInverseTransformHandle();
      Hybrid<double>* bdim_handle = result.getBoxDimensionsHandle();
      const std::vector<double> hx_loc = uniformRand(xrs, atom_count, 50.0);
      const std::vector<double> hy_loc = uniformRand(xrs, atom_count, 50.0);
      const std::vector<double> hz_loc = uniformRand(xrs, atom_count, 50.0);
      const std::vector<double> dx_loc = uniformRand(xrs, atom_count, 50.0);
      const std::vector<double> dy_loc = uniformRand(xrs, atom_count, 50.0);
      const std::vector<double> dz_loc = uniformRand(xrs, atom_count, 50.0);
      std::vector<double> h_umat, h_invu, h_bdim, d_umat, d_invu, d_bdim;
      createBoxDimensions(xrs, &h_umat, &h_invu, &h_bdim);
      createBoxDimensions(xrs, &d_umat, &d_invu, &d_bdim);
      switch (fmt) {
      case HybridFormat::HOST_ONLY:
        x_crd->putHost(hx_loc);
        y_crd->putHost(hy_loc);
        z_crd->putHost(hz_loc);
        switch (rvf[i]) {
        case TrajectoryKind::POSITIONS:
          umat_handle->putHost(h_umat);
          invu_handle->putHost(h_invu);
          bdim_handle->putHost(h_bdim);
          break;
        case TrajectoryKind::VELOCITIES:
        case TrajectoryKind::FORCES:
          break;
        }
        break;
#ifdef STORMM_USE_HPC
      case HybridFormat::HOST_MOUNTED:
      case HybridFormat::UNIFIED:
        x_crd->putHost(hx_loc);
        y_crd->putHost(hy_loc);
        z_crd->putHost(hz_loc);
        switch (rvf[i]) {
        case TrajectoryKind::POSITIONS:
          umat_handle->putHost(h_umat);
          invu_handle->putHost(h_invu);
          bdim_handle->putHost(h_bdim);
          break;
        case TrajectoryKind::VELOCITIES:
        case TrajectoryKind::FORCES:
          break;
        }
        break;
      case HybridFormat::DEVICE_ONLY:
        x_crd->putDevice(dx_loc);
        y_crd->putDevice(dy_loc);
        z_crd->putDevice(dz_loc);
        switch (rvf[i]) {
        case TrajectoryKind::POSITIONS:
          umat_handle->putDevice(d_umat);
          invu_handle->putDevice(d_invu);
          bdim_handle->putDevice(d_bdim);
          break;
        case TrajectoryKind::VELOCITIES:
        case TrajectoryKind::FORCES:
          break;
        }
        break;
      case HybridFormat::DECOUPLED:
      case HybridFormat::EXPEDITED:
        x_crd->putHost(hx_loc);
        y_crd->putHost(hy_loc);
        z_crd->putHost(hz_loc);
        x_crd->putDevice(dx_loc);
        y_crd->putDevice(dy_loc);
        z_crd->putDevice(dz_loc);
        switch (rvf[i]) {
        case TrajectoryKind::POSITIONS:
          umat_handle->putHost(h_umat);
          invu_handle->putHost(h_invu);
          bdim_handle->putHost(h_bdim);
          umat_handle->putDevice(d_umat);
          invu_handle->putDevice(d_invu);
          bdim_handle->putDevice(d_bdim);
          break;
        case TrajectoryKind::VELOCITIES:
        case TrajectoryKind::FORCES:
          break;
        }
        break;
#endif
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Test the copy constructor for a PhaseSpace object, with possible memory format changes.
//
// Arguments:
//   fmt_a:       The memory layout in which the original object shall be created
//   fmt_b:       The memory format of the newly copied object
//   atom_count:  The number of atoms to include in the object
//-------------------------------------------------------------------------------------------------
void testPhaseSpaceCopy(Xoroshiro128pGenerator *xrs, const HybridFormat fmt_a,
                        const HybridFormat fmt_b, const int atom_count) {
  PhaseSpace ps_a = detailPhaseSpace(atom_count, fmt_a, xrs);
  PhaseSpace ps_b(ps_a, fmt_b);
  Hybrid<double>* x_vel = ps_b.getCoordinateHandle(CartesianDimension::X,
                                                   TrajectoryKind::VELOCITIES,
                                                   CoordinateCycle::BLACK);
  std::vector<double> xh_tmp, xd_tmp;
  switch (fmt_b) {
  case HybridFormat::HOST_ONLY:
    xh_tmp = x_vel->readHost();
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    xh_tmp = x_vel->readHost();
    xd_tmp = x_vel->readDevice();
    break;
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::UNIFIED:
    xh_tmp = x_vel->readHost();
    break;
  case HybridFormat::DEVICE_ONLY:
    xd_tmp = x_vel->readDevice();
    break;
#endif
  }
  const bool has_host = (xh_tmp.size() > 0);
  const bool has_devc = (xd_tmp.size() > 0);
  for (int i = 0; i < atom_count; i++) {
    if (has_host) {
      xh_tmp[i] += 1.05;
    }
    if (has_devc) {
      xd_tmp[i] += 1.05;
    }
  }
  switch (fmt_b) {
  case HybridFormat::HOST_ONLY:
    x_vel->putHost(xh_tmp);
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    x_vel->putHost(xh_tmp);
    x_vel->putDevice(xd_tmp);
    break;
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::UNIFIED:
    x_vel->putHost(xh_tmp);
    break;
  case HybridFormat::DEVICE_ONLY:
    x_vel->putDevice(xd_tmp);
    break;
#endif
  }

  // Depending on the formats of the original and copied objects, run specific tests.
  runLegalComparisons(ps_a, ps_b);
}

//-------------------------------------------------------------------------------------------------
// Test the copy constructor for a CoordinateFrame object, with possible memory format changes.
//
// Arguments:
//   fmt_a:        The memory layout in which the original object shall be created
//   fmt_b:        The memory format of the newly copied object
//   atom_count:   The number of atoms to include in the object
//   frame_count:  The number of frames to include in the coordinate series
//   fp_bits:      Extent of the fixed-precision model's fractional component
//-------------------------------------------------------------------------------------------------
template <typename T>
void testCoordinateSeriesCopy(Xoroshiro128pGenerator *xrs, const HybridFormat fmt_a,
                              const HybridFormat fmt_b, const int atom_count,
                              const int frame_count, const int fp_bits = 0) {
  CoordinateSeries<T> cs_a(atom_count, frame_count, UnitCellType::TRICLINIC, fp_bits, fmt_a);
  Hybrid<T>* x_crd = cs_a.getFramesHandle(CartesianDimension::X);
  Hybrid<T>* y_crd = cs_a.getFramesHandle(CartesianDimension::Y);
  Hybrid<T>* z_crd = cs_a.getFramesHandle(CartesianDimension::Z);
  Hybrid<double>* umat_handle = cs_a.getBoxTransformsHandle();
  Hybrid<double>* invu_handle = cs_a.getInverseTransformsHandle();
  Hybrid<double>* bdim_handle = cs_a.getBoxDimensionsHandle();
  const double ascl = pow(2.0, fp_bits);
  const int padded_atom_count = roundUp(atom_count, warp_size_int);
  const int xfrm_stride = roundUp(9, warp_size_int);
  const int bdim_stride = roundUp(6, warp_size_int);
  for (int i = 0; i < frame_count; i++) {
    const std::vector<double> hx_pre = uniformRand(xrs, atom_count, 50.0 * ascl);
    const std::vector<double> hy_pre = uniformRand(xrs, atom_count, 50.0 * ascl);
    const std::vector<double> hz_pre = uniformRand(xrs, atom_count, 50.0 * ascl);
    const std::vector<double> dx_pre = uniformRand(xrs, atom_count, 50.0 * ascl);
    const std::vector<double> dy_pre = uniformRand(xrs, atom_count, 50.0 * ascl);
    const std::vector<double> dz_pre = uniformRand(xrs, atom_count, 50.0 * ascl);
    const std::vector<T> hx_loc(hx_pre.begin(), hx_pre.end());
    const std::vector<T> hy_loc(hy_pre.begin(), hy_pre.end());
    const std::vector<T> hz_loc(hz_pre.begin(), hz_pre.end());
    const std::vector<T> dx_loc(dx_pre.begin(), dx_pre.end());
    const std::vector<T> dy_loc(dy_pre.begin(), dy_pre.end());
    const std::vector<T> dz_loc(dz_pre.begin(), dz_pre.end());
    std::vector<double> h_umat, h_invu, h_bdim, d_umat, d_invu, d_bdim;
    createBoxDimensions(xrs, &h_umat, &h_invu, &h_bdim);
    createBoxDimensions(xrs, &d_umat, &d_invu, &d_bdim);
    switch (fmt_a) {
    case HybridFormat::HOST_ONLY:
      x_crd->putHost(hx_loc, i * padded_atom_count, atom_count);
      y_crd->putHost(hy_loc, i * padded_atom_count, atom_count);
      z_crd->putHost(hz_loc, i * padded_atom_count, atom_count);
      umat_handle->putHost(h_umat, i * xfrm_stride, 9);
      invu_handle->putHost(h_invu, i * xfrm_stride, 9);
      bdim_handle->putHost(h_bdim, i * bdim_stride, 6);
      break;
#ifdef STORMM_USE_HPC
    case HybridFormat::HOST_MOUNTED:
    case HybridFormat::UNIFIED:
      x_crd->putHost(hx_loc, i * padded_atom_count, atom_count);
      y_crd->putHost(hy_loc, i * padded_atom_count, atom_count);
      z_crd->putHost(hz_loc, i * padded_atom_count, atom_count);
      umat_handle->putHost(h_umat, i * xfrm_stride, 9);
      invu_handle->putHost(h_invu, i * xfrm_stride, 9);
      bdim_handle->putHost(h_bdim, i * bdim_stride, 6);
      break;
    case HybridFormat::DEVICE_ONLY:
      x_crd->putDevice(dx_loc, i * padded_atom_count, atom_count);
      y_crd->putDevice(dy_loc, i * padded_atom_count, atom_count);
      z_crd->putDevice(dz_loc, i * padded_atom_count, atom_count);
      umat_handle->putDevice(d_umat, i * xfrm_stride, 9);
      invu_handle->putDevice(d_invu, i * xfrm_stride, 9);
      bdim_handle->putDevice(d_bdim, i * bdim_stride, 6);
      break;
    case HybridFormat::DECOUPLED:
    case HybridFormat::EXPEDITED:
      x_crd->putHost(hx_loc, i * padded_atom_count, atom_count);
      y_crd->putHost(hy_loc, i * padded_atom_count, atom_count);
      z_crd->putHost(hz_loc, i * padded_atom_count, atom_count);
      x_crd->putDevice(dx_loc, i * padded_atom_count, atom_count);
      y_crd->putDevice(dy_loc, i * padded_atom_count, atom_count);
      z_crd->putDevice(dz_loc, i * padded_atom_count, atom_count);
      umat_handle->putHost(h_umat, i * xfrm_stride, 9);
      invu_handle->putHost(h_invu, i * xfrm_stride, 9);
      bdim_handle->putHost(h_bdim, i * bdim_stride, 6);
      umat_handle->putDevice(d_umat, i * xfrm_stride, 9);
      invu_handle->putDevice(d_invu, i * xfrm_stride, 9);
      bdim_handle->putDevice(d_bdim, i * bdim_stride, 6);
      break;
#endif
    }
  }
  CoordinateSeries<T> cs_b(cs_a, fmt_b);
  x_crd = cs_b.getFramesHandle(CartesianDimension::X);
  std::vector<T> xh_tmp, xd_tmp;
  switch (fmt_b) {
  case HybridFormat::HOST_ONLY:
    xh_tmp = x_crd->readHost();
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    xh_tmp = x_crd->readHost();
    xd_tmp = x_crd->readDevice();
    break;
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::UNIFIED:
    xh_tmp = x_crd->readHost();
    break;
  case HybridFormat::DEVICE_ONLY:
    xd_tmp = x_crd->readDevice();
    break;
#endif
  }
  const bool has_host = (xh_tmp.size() > 0);
  const bool has_devc = (xd_tmp.size() > 0);
  for (int i = 0; i < frame_count; i++) {
    for (int j = 0; j < atom_count; j++) {
      if (has_host) {
        xh_tmp[(i * padded_atom_count) + j] += static_cast<T>(1.05 * ascl);
      }
      if (has_devc) {
        xd_tmp[(i * padded_atom_count) + j] += static_cast<T>(1.05 * ascl);
      }
    }
  }
  switch (fmt_b) {
  case HybridFormat::HOST_ONLY:
    x_crd->putHost(xh_tmp);
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    x_crd->putHost(xh_tmp);
    x_crd->putDevice(xd_tmp);
    break;
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::UNIFIED:
    x_crd->putHost(xh_tmp);
    break;
  case HybridFormat::DEVICE_ONLY:
    x_crd->putDevice(xd_tmp);
    break;
#endif
  }

  // Depending on the formats of the original and copied objects, run specific tests.
  runLegalComparisons(cs_a, cs_b);
}

//-------------------------------------------------------------------------------------------------
// Test coordinate copying and format modification in the PhaseSpaceSynthesis object.  Descriptions
// of input variables follow from testCoordinateSeriesCopy(), above, with the exception of:
//
// Arguments:
//   top_name:      The topology file (needed to generate a valid system for the synthesis)
//   system_count:  The number of systems in the synthesis.  For these tests, all systems will be
//                  taken to have identical sizes.
//-------------------------------------------------------------------------------------------------
void testSynthesisCopy(Xoroshiro128pGenerator *xrs, const HybridFormat fmt_a,
                       const HybridFormat fmt_b, const std::string &top_name,
                       const int system_count, const int fp_bits = 36,
                       const GpuDetails &gpu = null_gpu) {
  std::vector<PhaseSpace> psv;
  const TestPriority do_tests = (getDrivePathType(top_name) == DrivePathType::FILE) ?
                                TestPriority::CRITICAL : TestPriority::ABORT;
  AtomGraph ag(top_name);
  std::vector<AtomGraph*> agv;
  const int atom_count = ag.getAtomCount();
  const int padded_atom_count = roundUp(atom_count, warp_size_int);
  for (int i = 0; i < system_count; i++) {
    psv.push_back(detailPhaseSpace(atom_count, fmt_a, xrs));
    agv.push_back(&ag);
  }
  const PhaseSpaceSynthesis poly_ps_a(psv, agv, fp_bits, 24, fp_bits, fp_bits, fmt_a, gpu);
  check(poly_ps_a.getFormat() == fmt_a, "The memory format of a PhaseSpaceSynthesis (" +
        getEnumerationName(poly_ps_a.getFormat()) + ") does not follow from that of its component "
        "PhaseSpace objects (" + getEnumerationName(fmt_a) + ").", do_tests);
  PhaseSpaceSynthesis poly_ps_b(poly_ps_a, fmt_b);
  Hybrid<llint>* xfrc = poly_ps_b.getCoordinateHandle(CartesianDimension::X,
                                                      TrajectoryKind::FORCES,
                                                      CoordinateCycle::WHITE);  
  Hybrid<int>* xfrc_ovrf = poly_ps_b.getCoordinateOverflowHandle(CartesianDimension::X,
                                                                 TrajectoryKind::FORCES,
                                                                 CoordinateCycle::WHITE);
  std::vector<llint> xh_tmp, xd_tmp;
  std::vector<int> xhovrf_tmp, xdovrf_tmp;
  switch (fmt_b) {
  case HybridFormat::HOST_ONLY:
    xh_tmp = xfrc->readHost();
    xhovrf_tmp = xfrc_ovrf->readHost();
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    xh_tmp = xfrc->readHost();
    xd_tmp = xfrc->readDevice();
    xhovrf_tmp = xfrc_ovrf->readHost();
    xdovrf_tmp = xfrc_ovrf->readDevice();
    break;
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::UNIFIED:
    xh_tmp = xfrc->readHost();
    xhovrf_tmp = xfrc_ovrf->readHost();
    break;
  case HybridFormat::DEVICE_ONLY:
    xd_tmp = xfrc->readDevice();
    xdovrf_tmp = xfrc_ovrf->readDevice();
    break;
#endif
  }
  const int sys_idx = poly_ps_b.getSystemCount() / 2;
  const bool has_host = (xh_tmp.size() > 0);
  const bool has_devc = (xd_tmp.size() > 0);
  const double incr = 1.05 * pow(2.0, fp_bits);
  for (int i = 0; i < atom_count; i++) {
    const size_t ilkp = (sys_idx * padded_atom_count) + i;
    if (has_host) {
      const double cfrc = hostInt95ToDouble(xh_tmp[ilkp], xhovrf_tmp[ilkp]);
      const int95_t icfrc = hostDoubleToInt95(cfrc + incr);
      xh_tmp[ilkp] = icfrc.x;
      xhovrf_tmp[ilkp] = icfrc.y;
    }
    if (has_devc) {
      const double cfrc = hostInt95ToDouble(xd_tmp[ilkp], xdovrf_tmp[ilkp]);
      const int95_t icfrc = hostDoubleToInt95(cfrc + incr);
      xd_tmp[ilkp] = icfrc.x;
      xdovrf_tmp[ilkp] = icfrc.y;
    }
  }
  switch (fmt_b) {
  case HybridFormat::HOST_ONLY:
    xfrc->putHost(xh_tmp);
    xfrc_ovrf->putHost(xhovrf_tmp);
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
    xfrc->putHost(xh_tmp);
    xfrc->putDevice(xd_tmp);
    xfrc_ovrf->putHost(xhovrf_tmp);
    xfrc_ovrf->putDevice(xdovrf_tmp);
    break;
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::UNIFIED:
    xfrc->putHost(xh_tmp);
    xfrc_ovrf->putHost(xhovrf_tmp);
    break;
  case HybridFormat::DEVICE_ONLY:
    xfrc->putDevice(xd_tmp);
    xfrc_ovrf->putDevice(xdovrf_tmp);
    break;
#endif
  }

  // Run tests of the correspondence, now expecting Cartesian X forces in one of the systems to
  // be inflated by 1.05 kcal/mol-A.
  runLegalComparisons(poly_ps_a, poly_ps_b);
}
    
//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main (const int argc, const char* argv[]) {

  // Some baseline initializations
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1
  section("Host-exclusive memory");

  // Section 2
  section("Device-exclusive memory");

  // Create some coordinate objects with random positions (in concept, an ideal gas)
  section(1);
  Xoroshiro128pGenerator xrs;
  testCoordinateFrameCopy(&xrs, HybridFormat::HOST_ONLY, HybridFormat::HOST_ONLY, 79);
  testPhaseSpaceCopy(&xrs, HybridFormat::HOST_ONLY, HybridFormat::HOST_ONLY, 181);
  testCoordinateSeriesCopy<float>(&xrs, HybridFormat::HOST_ONLY, HybridFormat::HOST_ONLY, 132, 7);
  testCoordinateSeriesCopy<int>(&xrs, HybridFormat::HOST_ONLY, HybridFormat::HOST_ONLY, 87, 9, 24);
  const char osc = osSeparator();
  const std::string top_file = oe.getStormmSourcePath() + osc + "test" + osc + "Topology" + osc +
                               "drug_example_iso.top";
  testSynthesisCopy(&xrs, HybridFormat::HOST_ONLY, HybridFormat::HOST_ONLY, top_file, 4, 36);

  // The following tests are only valid in HPC mode.
#ifdef STORMM_USE_HPC
  section(2);
  const std::vector<HybridFormat> all_fmt = { HybridFormat::HOST_ONLY, HybridFormat::HOST_MOUNTED,
                                              HybridFormat::DECOUPLED, HybridFormat::DEVICE_ONLY,
                                              HybridFormat::EXPEDITED, HybridFormat::UNIFIED };
  for (size_t i = 0; i < all_fmt.size(); i++) {
    for (size_t j = 0; j < all_fmt.size(); j++) {
      if (all_fmt[i] == HybridFormat::HOST_ONLY && all_fmt[j] == HybridFormat::HOST_ONLY) {
        continue;
      }
      testCoordinateFrameCopy(&xrs, all_fmt[i], all_fmt[j], 95);
      testPhaseSpaceCopy(&xrs, all_fmt[i], all_fmt[j], 181);
      testCoordinateSeriesCopy<float>(&xrs, all_fmt[i], all_fmt[j], 87, 9);
      testCoordinateSeriesCopy<int>(&xrs, all_fmt[i], all_fmt[j], 87, 9, 24);
      //testSynthesisCopy(&xrs, all_fmt[i], all_fmt[j], top_file, 4, 36);
    }
  }
#endif
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
