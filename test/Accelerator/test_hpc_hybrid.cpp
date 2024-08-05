#include <vector>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Accelerator/hybrid_util.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::constants::tiny;
using stormm::constants::ExceptionResponse;
using stormm::random::Ran2Generator;
using stormm::random::Xoroshiro128pGenerator;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using namespace stormm::card;
using namespace stormm::data_types;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Create a Hybrid object containing some double-precision numbers.  In addition to a process
// similar to what the corresponding function in the CPU-only test case performs, upload that
// result to the device and the further modify the host data.  Return the object with both data
// arrays.
//
// Arguments:
//   length:   the number of values to place in the object
//-------------------------------------------------------------------------------------------------
Hybrid<double> makeDoubleHybrid(int length) {
  Hybrid<double> result(length, "mdh_result");
  int mult = 2;
  int j = 0;
  double spacer = 1.375;
  for (int i = 0; i < length; i++) {
    result.putHost(static_cast<double>(mult * i) + spacer, i);
    j++;
    if (j == 2) {
      j = 0;
      mult *= -1;
      spacer += spacer * (mult * 0.65);
    }
  }
  result.upload();
  for (int i = 0; i < length; i++) {
    result.putHost(exp(sqrt(abs(1.0 / (result.readHost(i) * 7.872)))), i);
  }
  
  return result;
}

//-------------------------------------------------------------------------------------------------
// Populate a std::vector of scalars or tuples.  This routine overcomes compiler issues when the
// inputs could be scalars or tuples but the intended result may not be.  It is incumbent on the
// developer to use it properly.
//
// Arguments:
//   prep:  The array to be populated.  Its true type must be of type Tinterp but it will come from
//          a templated function where its data type is given as T.
//   x:     The first array of random number data with which to populate prep
//   y:     A second array of random number data with which to populate prep
//   z:     A third array of random number data with which to populate prep
//   w:     A fourth array of random number data with which to populate prep
//-------------------------------------------------------------------------------------------------
template <typename Tinterp, typename T>
void populateVector(std::vector<T> *prep, const std::vector<double> &x) {
  Tinterp* ptr = reinterpret_cast<Tinterp*>(prep->data());
  const size_t npts = x.size();
  for (size_t i = 0; i < npts; i++) {
    ptr[i] = x[i];
  }
}

template <typename Tinterp, typename T>
void populateVector(std::vector<T> *prep, const std::vector<double> &x,
                    const std::vector<double> &y) {
  if (x.size() != y.size()) {
    rtErr("The sizes of the input arrays (" + std::to_string(x.size()) + " and " +
          std::to_string(y.size()) + ") must match.", "populateVector");
  }
  Tinterp* ptr = reinterpret_cast<Tinterp*>(prep->data());
  const size_t npts = x.size();
  for (size_t i = 0; i < npts; i++) {
    ptr[i].x = x[i];
    ptr[i].y = y[i];
  }
}

template <typename Tinterp, typename T>
void populateVector(std::vector<T> *prep, const std::vector<double> &x,
                    const std::vector<double> &y, const std::vector<double> &z) {
  if (x.size() != y.size() || x.size() != z.size()) {
    rtErr("The sizes of the input arrays (" + std::to_string(x.size()) + ", " +
          std::to_string(y.size()) + ", and " + std::to_string(z.size()) + ") must match.",
          "populateVector");
  }
  Tinterp* ptr = reinterpret_cast<Tinterp*>(prep->data());
  const size_t npts = x.size();
  for (size_t i = 0; i < npts; i++) {
    ptr[i].x = x[i];
    ptr[i].y = y[i];
    ptr[i].z = z[i];
  }
}

template <typename Tinterp, typename T>
void populateVector(std::vector<T> *prep, const std::vector<double> &x,
                    const std::vector<double> &y, const std::vector<double> &z,
                    const std::vector<double> &w) {
  if (x.size() != y.size() || x.size() != z.size() || x.size() != w.size()) {
    rtErr("The sizes of the input arrays (" + std::to_string(x.size()) + ", " +
          std::to_string(y.size()) + ", " + std::to_string(z.size()) + ", and " +
          std::to_string(w.size()) + ") must match.", "populateVector");
  }
  Tinterp* ptr = reinterpret_cast<Tinterp*>(prep->data());
  const size_t npts = x.size();
  for (size_t i = 0; i < npts; i++) {
    ptr[i].x = x[i];
    ptr[i].y = y[i];
    ptr[i].z = z[i];
    ptr[i].w = w[i];
  }
}

//-------------------------------------------------------------------------------------------------
// Create a Hybrid array with randomized data in CPU host and / or GPU device memory.
//
// Arguments:
//   xrs:      Source of random data for the original array
//   hfmt:     Memory format in which to cast the original array object
//   npts:     The number of elements in the array
//   fp_bits:  Number of bits after the point for a fixed-precision model in the original array
//-------------------------------------------------------------------------------------------------
template <typename T>
Hybrid<T> buildHybridArray(Xoroshiro128pGenerator *xrs, const HybridFormat hfmt, const int npts,
                           const int fp_bits) {
  std::vector<T> host_prep(npts);
  std::vector<T> devc_prep(npts);
  const double ascl = pow(2.0, fp_bits);
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (isHpcVectorType<T>()) {
    std::vector<double> host_rand_x, host_rand_y, host_rand_z, host_rand_w;
    std::vector<double> devc_rand_x, devc_rand_y, devc_rand_z, devc_rand_w;
    if (isFloatingPointHpcVectorType<T>() || isSignedIntegralHpcVectorType<T>()) {
      host_rand_x = gaussianRand(xrs, npts, 5.0 * ascl);
      host_rand_y = gaussianRand(xrs, npts, 5.0 * ascl);
      devc_rand_x = gaussianRand(xrs, npts, 5.0 * ascl);
      devc_rand_y = gaussianRand(xrs, npts, 5.0 * ascl);
    }
    else {
      host_rand_x = uniformRand(xrs, npts, 5.0 * ascl);
      host_rand_y = uniformRand(xrs, npts, 5.0 * ascl);
      devc_rand_x = uniformRand(xrs, npts, 5.0 * ascl);
      devc_rand_y = uniformRand(xrs, npts, 5.0 * ascl);
    }
    if (ct == float2_type_index) {
      populateVector<float2, T>(&host_prep, host_rand_x, host_rand_y);
      populateVector<float2, T>(&devc_prep, devc_rand_x, devc_rand_y);
    }
    else if (ct == double2_type_index) {
      populateVector<double2, T>(&host_prep, host_rand_x, host_rand_y);
      populateVector<double2, T>(&devc_prep, devc_rand_x, devc_rand_y);
    }
    else if (ct == int2_type_index) {
      populateVector<int2, T>(&host_prep, host_rand_x, host_rand_y);
      populateVector<int2, T>(&devc_prep, devc_rand_x, devc_rand_y);
    }
    else if (ct == longlong2_type_index) {
      populateVector<longlong2, T>(&host_prep, host_rand_x, host_rand_y);
      populateVector<longlong2, T>(&devc_prep, devc_rand_x, devc_rand_y);
    }
    else if (ct == uint2_type_index) {
      populateVector<uint2, T>(&host_prep, host_rand_x, host_rand_y);
      populateVector<uint2, T>(&devc_prep, devc_rand_x, devc_rand_y);
    }
    else if (ct == ulonglong2_type_index) {
      populateVector<ulonglong2, T>(&host_prep, host_rand_x, host_rand_y);
      populateVector<ulonglong2, T>(&devc_prep, devc_rand_x, devc_rand_y);
    }
    else if (ct == float3_type_index || ct == double3_type_index || ct == int3_type_index ||
             ct == longlong3_type_index || ct == uint3_type_index || ct == ulonglong3_type_index) {
      if (isFloatingPointHpcVectorType<T>() || isSignedIntegralHpcVectorType<T>()) {
        host_rand_z = gaussianRand(xrs, npts, 5.0 * ascl);
        devc_rand_z = gaussianRand(xrs, npts, 5.0 * ascl);
      }
      else {
        host_rand_z = uniformRand(xrs, npts, 5.0 * ascl);
        devc_rand_z = uniformRand(xrs, npts, 5.0 * ascl);
      }
    }
    else if (ct == float4_type_index || ct == double4_type_index || ct == int4_type_index ||
             ct == longlong4_type_index || ct == uint4_type_index || ct == ulonglong4_type_index) {
      if (isFloatingPointHpcVectorType<T>() || isSignedIntegralHpcVectorType<T>()) {
        host_rand_z = gaussianRand(xrs, npts, 5.0 * ascl);
        host_rand_w = gaussianRand(xrs, npts, 5.0 * ascl);
        devc_rand_z = gaussianRand(xrs, npts, 5.0 * ascl);
        devc_rand_w = gaussianRand(xrs, npts, 5.0 * ascl);
      }
      else {
        host_rand_z = uniformRand(xrs, npts, 5.0 * ascl);
        host_rand_w = uniformRand(xrs, npts, 5.0 * ascl);
        devc_rand_z = uniformRand(xrs, npts, 5.0 * ascl);
        devc_rand_w = uniformRand(xrs, npts, 5.0 * ascl);
      }
    }
    else {
      rtErr("HPC vector type not recognized.", "testDeepCopyInner");
    }
  }
  else {
    std::vector<double> host_rand, devc_rand;
    if (isFloatingPointScalarType<T>() || isSignedIntegralScalarType<T>()) {
      host_rand = gaussianRand(xrs, npts, 5.0 * ascl);
      devc_rand = gaussianRand(xrs, npts, 5.0 * ascl);
    }
    else {
      host_rand = uniformRand(xrs, npts, 5.0 * ascl);
      devc_rand = uniformRand(xrs, npts, 5.0 * ascl);
    }
    if (ct == float_type_index) {
      populateVector<float, T>(&host_prep, host_rand);
      populateVector<float, T>(&devc_prep, devc_rand);
    }
    else if (ct == double_type_index) {
      populateVector<double, T>(&host_prep, host_rand);
      populateVector<double, T>(&devc_prep, devc_rand);
    }
    else if (ct == short_type_index) {
      populateVector<short int, T>(&host_prep, host_rand);
      populateVector<short int, T>(&devc_prep, devc_rand);
    }
    else if (ct == int_type_index) {
      populateVector<int, T>(&host_prep, host_rand);
      populateVector<int, T>(&devc_prep, devc_rand);
    }
    else if (ct == llint_type_index) {
      populateVector<llint, T>(&host_prep, host_rand);
      populateVector<llint, T>(&devc_prep, devc_rand);
    }
    else if (ct == ushort_type_index) {
      populateVector<ushort, T>(&host_prep, host_rand);
      populateVector<ushort, T>(&devc_prep, devc_rand);
    }
    else if (ct == uint_type_index) {
      populateVector<uint, T>(&host_prep, host_rand);
      populateVector<uint, T>(&devc_prep, devc_rand);
    }
    else if (ct == ullint_type_index) {
      populateVector<ullint, T>(&host_prep, host_rand);
      populateVector<ullint, T>(&devc_prep, devc_rand);
    }
    else {
      rtErr("Scalar data type not recognized.", "testDeepCopyInner");
    }
  }

  // Allocate the Hybrid object
  Hybrid<T> result(npts, "orig_array", hfmt, HybridKind::ARRAY);

  // Impart data to the original array
  switch (hfmt) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::HOST_ONLY:
    result.putHost(host_prep);
    break;
  case HybridFormat::DEVICE_ONLY:
    CHECK_THROWS(result.putHost(host_prep), "Data was placed in the host-side array of a " +
                 getEnumerationName(hfmt) + " format Hybrid object.");
    break;
  }
  switch (hfmt) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::DEVICE_ONLY:
    result.putDevice(devc_prep);
    break;
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::HOST_ONLY:
    CHECK_THROWS(result.putDevice(devc_prep), "Data was placed in the device-side array of a " +
                 getEnumerationName(hfmt) + " format Hybrid object.");
    break;
  case HybridFormat::UNIFIED:

    // Use the host-side random data for populating unified virtual memory
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Unpack the original and copied Hybrid arrays in terms of their host and device components.
// 
// Arguments:
//   orig_ar:
//   copy_ar:
//   orig_host:
//   orig_devc:
//   copy_host:
//   copy_devc:
//-------------------------------------------------------------------------------------------------
template <typename Torig, typename Tcopy>
void unpackOriginalAndCopy(const Hybrid<Torig> &orig_ar, const Hybrid<Tcopy> &copy_ar,
                           std::vector<Torig> *orig_host, std::vector<Torig> *orig_devc,
                           std::vector<Tcopy> *copy_host, std::vector<Tcopy> *copy_devc) {
  switch (orig_ar.getFormat()) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
    *orig_host = orig_ar.readHost();
    *orig_devc = orig_ar.readDevice();
    break;
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::HOST_ONLY:
    *orig_host = orig_ar.readHost();
    break;
  case HybridFormat::DEVICE_ONLY:
    *orig_devc = orig_ar.readDevice();
    break;
  }
  switch (copy_ar.getFormat()) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
    *copy_host = copy_ar.readHost();
    *copy_devc = copy_ar.readDevice();
    break;
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::HOST_ONLY:
    *copy_host = copy_ar.readHost();
    break;
  case HybridFormat::DEVICE_ONLY:
    *copy_devc = copy_ar.readDevice();
    break;
  }
}

//-------------------------------------------------------------------------------------------------
// Check the correspondence of two arrays on the host or device after a deep copying or recasting
// operation.
//
// Arguments:
//   copy_vec:
//   orig_vec:
//   copy_fmt:
//   orig_fmt:
//   copy_src:
//   orig_src:
//-------------------------------------------------------------------------------------------------
void checkArrayCorrespondence(const std::vector<double> &copy_vec,
                              const std::vector<double> &orig_vec,
                              const HybridFormat copy_fmt, const HybridFormat orig_fmt,
                              const HybridTargetLevel copy_src, const HybridTargetLevel orig_src,
                              const std::string &member_name = std::string("")) {
  std::string member_msg;
  if (member_name.size() > 0) {
    member_msg = "  Inconsistencies found in the \"" + member_name + "\" member of the tuple.";
  }
  check(copy_vec, RelationalOperator::EQUAL, orig_vec, "Copying a Hybrid object of format " +
        getEnumerationName(orig_fmt) + " to a Hybrid object of format " +
        getEnumerationName(copy_fmt) + " should result in the " + getEnumerationName(orig_src) +
        "-side data of the original object carrying over to the " + getEnumerationName(copy_src) +
        "-side data of the destination object." + member_msg);
}

//-------------------------------------------------------------------------------------------------
// Test whether an array of scalar data was copies as intended.
//
// Arguments:
//   orig_ar:       The original Hybrid array
//   copy_ar:       The copied Hybrid array (possibly cast into a different data type)
//   orig_fp_bits:  The fixed-precision bit count for the original array (provide as zero if no
//                  conversion between the original and copied arrays was to take place)
//   copy_fp_bits:  The fixed-precision bit count for the copied array
//-------------------------------------------------------------------------------------------------
template <typename Tointrp, typename Tcintrp, typename Torig, typename Tcopy>
void testScalarCopy(const Hybrid<Torig> &orig_ar, const Hybrid<Tcopy> &copy_ar,
                    const int orig_fp_bits = 0, const int copy_fp_bits = 0) {
  std::vector<Torig> orig_host, orig_devc;
  std::vector<Tcopy> copy_host, copy_devc;
  unpackOriginalAndCopy(orig_ar, copy_ar, &orig_host, &orig_devc, &copy_host, &copy_devc);
  const int npts = orig_ar.size();
  std::vector<double> orig_dh, orig_dd, copy_dh, copy_dd;
  const Tointrp* orig_phost = reinterpret_cast<Tointrp*>(orig_host.data());
  const Tointrp* orig_pdevc = reinterpret_cast<Tointrp*>(orig_devc.data());
  const Tcintrp* copy_phost = reinterpret_cast<Tcintrp*>(copy_host.data());
  const Tcintrp* copy_pdevc = reinterpret_cast<Tcintrp*>(copy_devc.data());
  if (orig_host.size() == npts) {
    orig_dh.resize(npts);
    for (size_t i = 0; i < npts; i++) {
      orig_dh[i] = orig_phost[i];
    }
  }
  if (orig_devc.size() == npts) {
    orig_dd.resize(npts);
    for (size_t i = 0; i < npts; i++) {
      orig_dd[i] = orig_pdevc[i];
    }
  }
  if (copy_host.size() == npts) {
    copy_dh.resize(npts);
    for (size_t i = 0; i < npts; i++) {
      copy_dh[i] = copy_phost[i];
    }
  }
  if (copy_devc.size() == npts) {
    copy_dd.resize(npts);
    for (size_t i = 0; i < npts; i++) {
      copy_dd[i] = copy_pdevc[i];
    }
  }

  // Check the host-side data in the copied Hybrid object
  switch (copy_ar.getFormat()) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
    switch (orig_ar.getFormat()) {
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_ONLY:
    case HybridFormat::HOST_MOUNTED:
      checkArrayCorrespondence(copy_dh, orig_dh, copy_ar.getFormat(), orig_ar.getFormat(),
                               HybridTargetLevel::HOST, HybridTargetLevel::HOST);
      break;
    case HybridFormat::DEVICE_ONLY:
      checkArrayCorrespondence(copy_dh, orig_dd, copy_ar.getFormat(), orig_ar.getFormat(),
                               HybridTargetLevel::HOST, HybridTargetLevel::DEVICE);
      break;
    }
    break;
  case HybridFormat::DEVICE_ONLY:
    break;
  }

  // Check the device-side data in the copied Hybrid object
  switch (copy_ar.getFormat()) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::DEVICE_ONLY:
    switch (orig_ar.getFormat()) {
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::DEVICE_ONLY:
      checkArrayCorrespondence(copy_dd, orig_dd, copy_ar.getFormat(), orig_ar.getFormat(),
                               HybridTargetLevel::DEVICE, HybridTargetLevel::DEVICE);
      break;
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_ONLY:
    case HybridFormat::HOST_MOUNTED:
      checkArrayCorrespondence(copy_dd, orig_dh, copy_ar.getFormat(), orig_ar.getFormat(),
                               HybridTargetLevel::DEVICE, HybridTargetLevel::HOST);
      break;
    }
    break;
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
// Test whether an array of two-part tuples was copied as intended.  Descriptions of input
// parameters follow from testScalarCopy(), above.  The
//-------------------------------------------------------------------------------------------------
template <typename Tinterp, typename Torig>
void testTuple2Copy(const Hybrid<Torig> &orig_ar, const Hybrid<Torig> &copy_ar,
                    const int orig_fp_bits = 0, const int copy_fp_bits = 0) {
  std::vector<Torig> orig_host, orig_devc;
  std::vector<Torig> copy_host, copy_devc;
  unpackOriginalAndCopy(orig_ar, copy_ar, &orig_host, &orig_devc, &copy_host, &copy_devc);
  const int npts = orig_ar.size();
  std::vector<double> orig_dhx, orig_ddx, copy_dhx, copy_ddx;
  std::vector<double> orig_dhy, orig_ddy, copy_dhy, copy_ddy;
  const Tinterp* orig_phost = reinterpret_cast<Tinterp*>(orig_host.data());
  const Tinterp* orig_pdevc = reinterpret_cast<Tinterp*>(orig_devc.data());
  const Tinterp* copy_phost = reinterpret_cast<Tinterp*>(copy_host.data());
  const Tinterp* copy_pdevc = reinterpret_cast<Tinterp*>(copy_devc.data());
  if (orig_host.size() == npts) {
    orig_dhx.resize(npts);
    orig_dhy.resize(npts);
    for (size_t i = 0; i < npts; i++) {
      orig_dhx[i] = orig_phost[i].x;
      orig_dhy[i] = orig_phost[i].y;
    }
  }
  if (orig_devc.size() == npts) {
    orig_ddx.resize(npts);
    orig_ddy.resize(npts);
    for (size_t i = 0; i < npts; i++) {
      orig_ddx[i] = orig_pdevc[i].x;
      orig_ddy[i] = orig_pdevc[i].y;
    }
  }
  if (copy_host.size() == npts) {
    copy_dhx.resize(npts);
    copy_dhy.resize(npts);
    for (size_t i = 0; i < npts; i++) {
      copy_dhx[i] = copy_phost[i].x;
      copy_dhy[i] = copy_phost[i].y;
    }
  }
  if (copy_devc.size() == npts) {
    copy_ddx.resize(npts);
    copy_ddy.resize(npts);
    for (size_t i = 0; i < npts; i++) {
      copy_ddx[i] = copy_pdevc[i].x;
      copy_ddy[i] = copy_pdevc[i].y;
    }
  }

  // Check the host-side data in the copied Hybrid object
  switch (copy_ar.getFormat()) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
    switch (orig_ar.getFormat()) {
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_ONLY:
    case HybridFormat::HOST_MOUNTED:
      checkArrayCorrespondence(copy_dhx, orig_dhx, copy_ar.getFormat(), orig_ar.getFormat(),
                               HybridTargetLevel::HOST, HybridTargetLevel::HOST, "x");
      checkArrayCorrespondence(copy_dhy, orig_dhy, copy_ar.getFormat(), orig_ar.getFormat(),
                               HybridTargetLevel::HOST, HybridTargetLevel::HOST, "y");
      break;
    case HybridFormat::DEVICE_ONLY:
      checkArrayCorrespondence(copy_dhx, orig_ddx, copy_ar.getFormat(), orig_ar.getFormat(),
                               HybridTargetLevel::HOST, HybridTargetLevel::DEVICE, "x");
      checkArrayCorrespondence(copy_dhy, orig_ddy, copy_ar.getFormat(), orig_ar.getFormat(),
                               HybridTargetLevel::HOST, HybridTargetLevel::DEVICE, "y");
      break;
    }
    break;
  case HybridFormat::DEVICE_ONLY:
    break;
  }

  // Check the device-side data in the copied Hybrid object
  switch (copy_ar.getFormat()) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::DEVICE_ONLY:
    switch (orig_ar.getFormat()) {
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::DEVICE_ONLY:
      checkArrayCorrespondence(copy_ddx, orig_ddx, copy_ar.getFormat(), orig_ar.getFormat(),
                               HybridTargetLevel::DEVICE, HybridTargetLevel::DEVICE, "x");
      checkArrayCorrespondence(copy_ddy, orig_ddy, copy_ar.getFormat(), orig_ar.getFormat(),
                               HybridTargetLevel::DEVICE, HybridTargetLevel::DEVICE, "y");
      break;
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_ONLY:
    case HybridFormat::HOST_MOUNTED:
      checkArrayCorrespondence(copy_ddx, orig_dhx, copy_ar.getFormat(), orig_ar.getFormat(),
                               HybridTargetLevel::DEVICE, HybridTargetLevel::HOST, "x");
      checkArrayCorrespondence(copy_ddy, orig_dhy, copy_ar.getFormat(), orig_ar.getFormat(),
                               HybridTargetLevel::DEVICE, HybridTargetLevel::HOST, "y");
      break;
    }
    break;
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
// Confirm that the contents of one Hybrid object have been appropriately transferred to another,
// given data copy priorities.  Descriptions of input parameters follow from testScalarCopy(),
// above.
//-------------------------------------------------------------------------------------------------
template <typename Torig, typename Tcopy>
void confirmHybridCopying(const Hybrid<Torig> &orig_ar, const Hybrid<Tcopy> &copy_ar,
                          const int orig_fp_bits = 0, const int copy_fp_bits = 0) {

  // Create void-casted pointers to the data so that tuples can also be examined with the same
  // routine.
  const size_t ct_orig = std::type_index(typeid(Torig)).hash_code();
  const size_t ct_copy = std::type_index(typeid(Tcopy)).hash_code();
  if (isScalarType<Torig>() && isScalarType<Tcopy>()) {
    if (ct_orig == float_type_index && ct_copy == float_type_index) {
      testScalarCopy<float, float, Torig, Tcopy>(orig_ar, copy_ar, orig_fp_bits, copy_fp_bits);
    }
    else if (ct_orig == double_type_index && ct_copy == double_type_index) {
      testScalarCopy<double, double, Torig, Tcopy>(orig_ar, copy_ar, orig_fp_bits, copy_fp_bits);
    }
    else if (ct_orig == int_type_index && ct_copy == int_type_index) {
      testScalarCopy<int, int, Torig, Tcopy>(orig_ar, copy_ar, orig_fp_bits, copy_fp_bits);
    }
    else if (ct_orig == ullint_type_index && ct_copy == ullint_type_index) {
      testScalarCopy<ullint, ullint, Torig, Tcopy>(orig_ar, copy_ar, orig_fp_bits, copy_fp_bits);
    }
  }
  else if (ct_orig == float2_type_index) {
    testTuple2Copy<float2, Torig>(orig_ar, copy_ar, orig_fp_bits, copy_fp_bits);
  }
  else if (ct_orig == double2_type_index) {
    testTuple2Copy<double2, Torig>(orig_ar, copy_ar, orig_fp_bits, copy_fp_bits);
  }
  else if (ct_orig == int2_type_index) {
    testTuple2Copy<int2, Torig>(orig_ar, copy_ar, orig_fp_bits, copy_fp_bits);
  }
  else if (ct_orig == uint2_type_index) {
    testTuple2Copy<uint2, Torig>(orig_ar, copy_ar, orig_fp_bits, copy_fp_bits);
  }
  else {
    if (isHpcVectorType<Torig>()) {
      rtErr("HPC vector data type " + getStormmHpcVectorTypeName<Torig>() + " is not handled by "
            "the testing routines.", "confirmHybridCopying");
    }
    else {
      rtErr("Copy operation between unknown data types is not supported by the testing "
            "routines.", "confirmHybridCopying");
    }
  }
}

//-------------------------------------------------------------------------------------------------
// A templated function to handle various instances of the testDeepCopy() function.
//
// Arguments:
//   xrs:       Source of random data for the original array
//   format_a:  Memory format in which to cast the original array object
//   fmt_copy:  Memory format of the copied object
//   gpu:       Details of the GPU that will perform the copy
//   fp_bits:   Number of bits after the point for a fixed-precision model in the original array
//-------------------------------------------------------------------------------------------------
template <typename T>
void testDeepCopyInner(Xoroshiro128pGenerator *xrs, const HybridFormat fmt_orig,
                       const HybridFormat fmt_copy, const GpuDetails &gpu, const int npts = 76,
                       const int fp_bits = 0) {

  // Allocate the Hybrid objects
  std::vector<T> host_prep, devc_prep;
  const Hybrid<T> o_ar = buildHybridArray<T>(xrs, fmt_orig, npts, fp_bits);
  Hybrid<T> c_ar(npts, "copy_array", fmt_copy, HybridKind::ARRAY);

  // Impart data to the new array and check the result.
  const std::string data_type_name = (isHpcVectorType<T>()) ? getStormmHpcVectorTypeName<T>() :
                                                              getStormmScalarTypeName<T>();
  deepCopy(&c_ar, o_ar, gpu);

  // Check the correspondence of the results
  confirmHybridCopying(o_ar, c_ar);
}

//-------------------------------------------------------------------------------------------------
// Test the Hybrid array deep copy and deep recasting functionality.
//
// Arguments:
//   xrs:       Source of random data for the original array
//   fmt_orig:  Memory format in which to cast the original array object
//   fmt_copy:  Memory format of the copied object
//   gpu:       Details of the GPU that will perform the copy
//-------------------------------------------------------------------------------------------------
void testDeepCopy(Xoroshiro128pGenerator *xrs, const HybridFormat fmt_orig,
                  const HybridFormat fmt_copy, const GpuDetails &gpu) {
  testDeepCopyInner<float>(xrs, fmt_orig, fmt_copy, gpu);
  testDeepCopyInner<double>(xrs, fmt_orig, fmt_copy, gpu);
  testDeepCopyInner<int>(xrs, fmt_orig, fmt_copy, gpu, 18);
  testDeepCopyInner<ullint>(xrs, fmt_orig, fmt_copy, gpu, 45);
  testDeepCopyInner<int2>(xrs, fmt_orig, fmt_copy, gpu, 21);
  testDeepCopyInner<float2>(xrs, fmt_orig, fmt_copy, gpu);
  testDeepCopyInner<uint2>(xrs, fmt_orig, fmt_copy, gpu, 20);
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Obtain environment variables or command-line input, if available
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> engage_the_gpu(1);
  
  // Section 1
  section("Basic hybrid functionality, with an emphasis on the HPC components");

  // Fill two Hybrid objects, one with integers and the other with near integer values
  const int n_pts = 100;
  Hybrid<int>    ihb_a(n_pts, "ihb_a", HybridFormat::EXPEDITED, HybridKind::ARRAY);
  Hybrid<double> dhb_a(n_pts, "dhb_a", HybridFormat::EXPEDITED, HybridKind::ARRAY);
  std::vector<int> ihb_host_contents(n_pts);
  std::vector<double> dhb_host_contents(n_pts);
  double perturb = 0.05;
  double delta = 0.01;
  for (int i = 0; i < n_pts; i++) {
    ihb_a.putHost(8 + i, i);
    dhb_a.putHost(8.0 + perturb + i, i);
    ihb_host_contents[i] = 8 + i;
    dhb_host_contents[i] = 8.0 + perturb + i;
    if (perturb > 0.045) {
      delta = -0.1;
    }
    else if (perturb < -0.045) {
      delta = 0.1;
    }
    perturb += delta;
  }
  ihb_a.upload();
  dhb_a.upload();
  const std::vector<int>    ihb_devc_contents = ihb_a.readDevice();
  const std::vector<double> dhb_devc_contents = dhb_a.readDevice();
  check(ihb_devc_contents, RelationalOperator::EQUAL, ihb_host_contents, "An ARRAY-kind "
        "Hybrid<int> does not properly upload its contents, or the readDevice() member function "
        "is not working as it should.");
  dhb_a.resize(2 * n_pts);
  Ran2Generator prng(8726159);
  std::vector<double> random_pts(n_pts);
  for (int i = 0; i < n_pts; i++) {
    random_pts[i] = prng.gaussianRandomNumber();
  }
  dhb_host_contents.insert(dhb_host_contents.end(), random_pts.begin(), random_pts.end());
  dhb_a.putDevice(random_pts, n_pts, n_pts);
  check(std::vector<double>(n_pts, 0.0), RelationalOperator::EQUAL, dhb_a.readHost(n_pts, n_pts),
        "Hybrid<double> should resize to contain zeros in its new elements.");
  dhb_a.download();
  check(random_pts, RelationalOperator::EQUAL, dhb_a.readHost(n_pts, n_pts), "Hybrid<double> did "
        "not properly download new data added on the device side.");
  ihb_a.resize(n_pts / 2);
  check(ihb_a.size(), RelationalOperator::EQUAL, n_pts / 2, "Hybrid<int> did not properly "
        "resize.");
  CHECK_THROWS(ihb_a.putDevice(5, n_pts), "A Hybrid<int> object allows data to be placed on the "
               "device in a position beyond its range.");
  CHECK_THROWS(ihb_a.readDevice(n_pts, n_pts), "A Hybrid<int> object allows data to be read from "
               "the device in positions beyond its range.");
  const Hybrid<double> dhb_a_copy = dhb_a;
  check(dhb_a_copy.readDevice(n_pts, n_pts), RelationalOperator::EQUAL, random_pts, "A copy of a "
        "Hybrid<double> object does not contain the right data on the device.");
  check(dhb_a_copy.readHost(), RelationalOperator::EQUAL, dhb_host_contents, "A copy of a "
        "Hybrid<double> object does not contain the right data on the host.");
  dhb_a.putHost(random_pts, 0, n_pts);
  check(dhb_a_copy.readHost(), RelationalOperator::EQUAL, dhb_host_contents, "A copy of a "
        "Hybrid<double> object does not contain data independent of the original.");
  Hybrid<int4> forces(768, "water_forces", HybridFormat::DEVICE_ONLY, HybridKind::ARRAY);
  Hybrid<int4> frc_ptr(HybridKind::POINTER, "water_force_pointer", HybridFormat::DEVICE_ONLY);
  std::vector<int4> ifrc(768);
  for (int i = 0; i < 768; i++) {
    ifrc[i].x = 1024.0 * prng.gaussianRandomNumber();
    ifrc[i].y = 1024.0 * prng.gaussianRandomNumber();
    ifrc[i].z = 1024.0 * prng.gaussianRandomNumber();
    ifrc[i].w =  512.0 * prng.gaussianRandomNumber();
  }
  CHECK_THROWS(forces.putHost(ifrc), "A Hybrid object with DEVICE_ONLY data accepted data for its "
               "non-existent host array.");
  forces.putDevice(ifrc);
  std::vector<int4> force_buffer = forces.readDevice(256, 256);
  check(force_buffer.size(), RelationalOperator::EQUAL, 256LLU, "A DEVICE_ONLY Hybrid<int4> "
        "object did not return the correct amount of data upon request.");
  check(force_buffer[76].z, RelationalOperator::EQUAL, ifrc[256 + 76].z, "Data collected from a "
        "DEVICE_ONLY Hybrid<int4> object does not meet expectations.");
  frc_ptr.setPointer(&forces, 128, 384);
  std::vector<int4> force_buffer_ii = frc_ptr.readDevice(128, 256);
  check(force_buffer[91].y, RelationalOperator::EQUAL, force_buffer_ii[91].y, "Data collected "
        "from a POINTER-kind Hybrid does not match data collected from its target.");
  check(force_buffer.size(), RelationalOperator::EQUAL, 256LLU, "Data collected "
        "from a POINTER-kind Hybrid does not have the expeceted size.");
  Hybrid<double> dhb_blank;
  dhb_blank.resize(n_pts);
  dhb_blank.putHost(random_pts);
  dhb_blank.upload();
  std::vector<double> blank_returns = dhb_blank.readDevice();
  check(blank_returns, RelationalOperator::EQUAL, random_pts, "Resizing an empty Hybrid<double> "
        "object to hold data did not work.");

  // Test the copy assignment operator.  In this modified version of the corresponding CPU-only
  // test, we upload the data to the scope-sequestered Hybrid object's GPU array, then further
  // modify the host data, to verify that the copy grabs everything.
  Hybrid<double> td_copied;
  if (true) {
    const int hidden_pts = 16;
    Hybrid<double> new_hb(hidden_pts, "Not dead");
    for (int i = 0; i < hidden_pts; i++) {
      new_hb.putHost((3.0 * static_cast<double>(i) * prng.gaussianRandomNumber()) + 1.7, i);
    }
    new_hb.upload();
    for (int i = 0; i < hidden_pts; i++) {
      new_hb.putHost((3.0 * static_cast<double>(i) * prng.gaussianRandomNumber()) - 3.8, i);
    }
    td_copied = new_hb;
  }
  check(td_copied.size(), RelationalOperator::EQUAL, 16, "Assignment of a Hybrid object does not "
        "make the object the correct size.");
  check(td_copied.readHost(4), RelationalOperator::EQUAL, Approx(-22.76199814).margin(1.0e-6),
        "The Hybrid assignment operator does not do a deep copy of the values on the host.");
  check(td_copied.readDevice(7), RelationalOperator::EQUAL, Approx(-21.48105323).margin(1.0e-6),
        "The Hybrid assignment operator does not do a deep copy of the values on the device.");
  check(std::string(td_copied.getLabel().name), RelationalOperator::EQUAL, "Not dead",
        "The Hybrid assignment operator does not carry over the proper label.");

  // Test the move constructor
  std::vector<Hybrid<int>> stl_vec_of_hybrids;
  std::vector<Hybrid<double>> stl_vec_of_patterns;
  std::vector<int> host_sum_answer;
  std::vector<int> devc_sum_answer;
  for (int i = 1; i < 10; i++) {
    Hybrid<int> tmp_hyb(i, std::string("tmp_hyb_" + std::to_string(i)).c_str());
    for (int j = 0; j < i; j++) {
      tmp_hyb.putHost(i + 2, j);
    }
    tmp_hyb.upload();
    for (int j = 0; j < i; j++) {
      tmp_hyb.putHost((9 * i) - 7, j);
    }
    const std::vector<int> tmp_host_contents = tmp_hyb.readHost();
    const std::vector<int> tmp_devc_contents = tmp_hyb.readDevice();
    host_sum_answer.push_back(sum<int>(tmp_host_contents));
    devc_sum_answer.push_back(sum<int>(tmp_devc_contents));
    stl_vec_of_hybrids.push_back(tmp_hyb);
    stl_vec_of_patterns.push_back(makeDoubleHybrid((2 * i) - 1));
  }
  std::vector<int> hyb_host_sums;
  std::vector<int> hyb_devc_sums;
  for (size_t i = 0; i < stl_vec_of_hybrids.size(); i++) {
    const std::vector<int> tmp_host_contents = stl_vec_of_hybrids[i].readHost();
    const std::vector<int> tmp_devc_contents = stl_vec_of_hybrids[i].readDevice();
    hyb_host_sums.push_back(sum<int>(tmp_host_contents));
    hyb_devc_sums.push_back(sum<int>(tmp_devc_contents));
  }
  int npatterns = 0;
  for (size_t i = 0; i < stl_vec_of_patterns.size(); i++) {
    npatterns += stl_vec_of_patterns[i].size();
  }
  check(host_sum_answer, RelationalOperator::EQUAL, hyb_host_sums, "Hybrid arrays committed to a "
        "Standard Template Library vector via its push_back() member function and the Hybrid copy "
        "constructor were not maintained correctly on the host side.");
  check(devc_sum_answer, RelationalOperator::EQUAL, hyb_devc_sums, "Hybrid arrays committed to a "
        "Standard Template Library vector via its push_back() member function and the Hybrid copy "
        "constructor were not maintained correctly on the device side.");
  const Hybrid<double> result_of_function = makeDoubleHybrid(10);
  check(sum<double>(result_of_function.readHost()), RelationalOperator::EQUAL,
        Approx(11.56283584).margin(1.0e-6), "A Hybrid array returned by a function does not "
        "contain the expected host data.");
  check(sum<double>(result_of_function.readDevice()), RelationalOperator::EQUAL,
        Approx(19.90602500).margin(1.0e-6), "A Hybrid array returned by a function does not "
        "contain the expected device data.");
  check(npatterns, RelationalOperator::EQUAL, 81, "Hybrid arrays transferred by the object's move "
        "constructor into a Standard Template Library vector via the push_back() method do not "
        "maintain their proper sizes.");
  std::vector<double> host_pattern_sums(stl_vec_of_patterns.size());
  std::vector<double> devc_pattern_sums(stl_vec_of_patterns.size());
  for (size_t i = 0; i < stl_vec_of_patterns.size(); i++) {

    // Take the sum of the current host data, then download the device data to the host
    // to get its sum.  The sum() function operates on the host data of a Hybrid object
    // by default.
    host_pattern_sums[i] = sum<double>(stl_vec_of_patterns[i]);
    stl_vec_of_patterns[i].download();
    devc_pattern_sums[i] = sum<double>(stl_vec_of_patterns[i]);
  }
  const std::vector<double> host_pattern_sums_ans = {   1.35520505,   3.75422996,   6.04901199,
                                                        8.28452913,  10.47681445,  12.64537424,
                                                       14.80017277,  16.94308412,  19.07518675 };
  const std::vector<double> devc_pattern_sums_ans = {   1.37500000,   0.33750000,   0.97625000,
                                                       -1.68787500,   1.25138750,  -0.29036625,
                                                        1.06154263,  -1.25464728,   1.19253559 };
  check(host_pattern_sums, RelationalOperator::EQUAL, Approx(host_pattern_sums_ans).margin(1.0e-6),
        "Hybrid arrays transferred by the object's move constructor into a Standard Template "
        "Library vector via the push_back() method do not retain the correct host data.");
  check(devc_pattern_sums, RelationalOperator::EQUAL, Approx(devc_pattern_sums_ans).margin(1.0e-6),
        "Hybrid arrays transferred by the object's move constructor into a Standard Template "
        "Library vector via the push_back() method do not retain the correct device data.");

  // Check HPC_enabled deep copy and memory format reassignments
  Xoroshiro128pGenerator xrs;
  const std::vector<HybridFormat> all_fmt = { HybridFormat::HOST_ONLY, HybridFormat::DEVICE_ONLY,
                                              HybridFormat::EXPEDITED, HybridFormat::HOST_MOUNTED,
                                              HybridFormat::DECOUPLED, HybridFormat::UNIFIED };
  for (size_t i = 0; i < all_fmt.size(); i++) {
    for (size_t j = 0; j < all_fmt.size(); j++) {
      testDeepCopy(&xrs, all_fmt[i], all_fmt[j], gpu);
    }
  }
  
  // Print a summary of tests run
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
