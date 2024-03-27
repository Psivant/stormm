#include <vector>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::constants::tiny;
using stormm::constants::ExceptionResponse;
using stormm::random::Ran2Generator;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using namespace stormm::card;
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
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Obtain environment variables or command-line input, if available
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

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
  
  // Print a summary of tests run
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
