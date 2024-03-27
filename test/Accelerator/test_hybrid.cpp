#include "copyright.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/Math/summation.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/UnitTesting/unit_test.h"

#ifndef STORMM_USE_HPC
using stormm::data_types::int4;
using stormm::data_types::ushort2;
#endif
using stormm::stmath::sum;
using stormm::random::Ran2Generator;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using namespace stormm::card;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Check to see if two vectors of int4 tuples match in every element for a specified number of
// points.
//
// Arguments:
//   a:     The first vector
//   b:     The secdond vector
//   npts:  The number of points to check, starting from the beginning of each vector
//-------------------------------------------------------------------------------------------------
bool int4_series_match(const std::vector<int4> a, const std::vector<int4> b, size_t npts) {
  if (npts > a.size() || npts > b.size()) {
    printf("Bad inputs to int4_series_match.\n");
  }
  bool agreement = true;
  for (int i = 0; i < npts; i++) {
    agreement = (agreement && a[i].x == b[i].x);
    agreement = (agreement && a[i].y == b[i].y);
    agreement = (agreement && a[i].z == b[i].z);
    agreement = (agreement && a[i].w == b[i].w);
  }
  return agreement;
}

//-------------------------------------------------------------------------------------------------
// Create a Hybrid object containing some double-precision numbers.
//
// Arguments:
//   length:   the number of values to place in the object
//-------------------------------------------------------------------------------------------------
Hybrid<double> makeDoubleHybrid(int length) {
  Hybrid<double> result(length, "mdh_result");
  int mult = 2;
  int j = 0;
  double spacer = 0.75;
  for (int i = 0; i < length; i++) {    
    result.putHost(static_cast<double>(mult * i) + spacer, i);
    j++;
    if (j == 3) {
      j = 0;
      mult *= -1;
      spacer += spacer * (mult * 0.65);
    }
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

  // "Forward declaration" of test sections: section 1
  section("Basic Hybrid object functionality, spanning getter functions as well as data entry "
	  "and retrieval");

  // Section 2
  section("Hybrid object failsafes, exception handling for definite code errors");

  // Section 3
  section("Global memory ledger features");

  // Section 4
  section("Hybrid POINTER-kind object veracity");

  // Section 5
  section("Hybrid copy and move constructors");

  // Fill two Hybrid objects, one with integers and the other with near integer values
  const int n_pts = 100;
  Hybrid<int>    ihybrid_a(n_pts, "ihybrid_a", HybridFormat::HOST_ONLY, HybridKind::ARRAY);
  Hybrid<double> dhybrid_a(n_pts, "dhybrid_a", HybridFormat::HOST_ONLY, HybridKind::ARRAY);
  double perturb = 0.05;
  double delta = 0.01;
  for (int i = 0; i < n_pts; i++) {
    ihybrid_a.putHost(8 + i, i);
    dhybrid_a.putHost(8.0 + perturb + i, i);
    if (perturb > 0.045) {
      delta = -0.1;
    }
    else if (perturb < -0.045) {
      delta = 0.1;
    }
    perturb += delta;
  }

  // Basic comparisons of the two vectors
  section(1);
  check(ihybrid_a.readHost() ==
        Approx(dhybrid_a.readHost(), ComparisonType::ABSOLUTE).margin(0.06),
        "Failure in comparing integer and double-precision Hybrid objects that should contain "
        "roughly similar values");
  check(dhybrid_a.readHost() !=
        Approx(ihybrid_a.readHost(), ComparisonType::ABSOLUTE).margin(0.04),
        "Failure to detect variation in integer and double-precision Hybrid objects once the "
        "tolerance was tightened");

  // Make a pointer to one of the hybrid objects
  Hybrid<double> dhybrid_ptr(n_pts, "dhybrid_ptr", HybridFormat::HOST_ONLY, HybridKind::POINTER);
  dhybrid_ptr.setPointer(&dhybrid_a);

  // Test each of the Hybrid getter member functions
  check(ihybrid_a.getKind() == HybridKind::ARRAY, "Hybrid::getKind() fails to return correct "
        "\"ARRAY\" result.");
  check(dhybrid_ptr.getKind() == HybridKind::POINTER, "Hybrid::getKind() fails to return "
        "correct \"POINTER\" result.");
  const std::string fmt_err_msg("Hybrid::getFormat() gives incorrect reports of the properties "
                                "of Hybrid object ");
  check(ihybrid_a.getFormat() == HybridFormat::HOST_ONLY, fmt_err_msg + "\"ihybrid_a\".");
  check(dhybrid_ptr.getFormat() == HybridFormat::HOST_ONLY, fmt_err_msg + "\"dhybrid_ptr\".");
  HybridLabel hlbl = dhybrid_a.getLabel();
  check(std::string(hlbl.name) == "dhybrid_a", "Hybrid object \"dhybrid_a\" carries label " +
        std::string(hlbl.name));
  check(hlbl.serial_number, RelationalOperator::EQUAL, 1, "The wrong serial number was assigned "
        "to Hybrid object \"dhybrid_a\".");
  check(dhybrid_ptr.readHost(3), RelationalOperator::EQUAL, Approx(10.95).margin(1.0e-6),
        "Hybrid readHost() member function fails to report the proper value of a double-precision "
        "POINTER kind object.");
  check(ihybrid_a.getElementSize(), RelationalOperator::EQUAL, sizeof(int), "Hybrid object "
        "\"ihybrid_a\" records the wrong element size.");
  check(dhybrid_a.getAllocations(), RelationalOperator::EQUAL, 1, "Hybrid object \"dhybrid_a\" "
	"records the wrong number of allocations.");
  
  // Put something into the POINTER object, then check the target ARRAY
  dhybrid_ptr.putHost(85.4, 7);
  check(dhybrid_a.readHost(7), RelationalOperator::EQUAL, Approx(85.4).margin(1.0e-6),
        "Hybrid putHost() member function acting on a POINTER object does not change the correct "
        "value in the target ARRAY object.");

  // Try putting something in the POINTER object that exceeds the bounds of the underlying ARRAY
  section(2);
  CHECK_THROWS(dhybrid_ptr.putHost(93.1, 10 * n_pts), "Hybrid::putHost() failed to implement its "
               "bounds check.");

  // Reset the POINTER-kind Hybrid object and verify that the length changes as appropriate
  section(3);
  std::vector<LedgerEntry> known_hybrids = gbl_mem_balance_sheet.getEntry("dhybrid_ptr");
  check(known_hybrids.size(), RelationalOperator::EQUAL, 1, "Ledger::getEntry() fails to "
	"identify the correct number of hybrid objects named \"dhybrid_ptr\".");
  check(known_hybrids[0].length, RelationalOperator::EQUAL, n_pts, "The global memory ledger "
        "fails to record the correct length of \"dhybrid_ptr\" when targeting \"dhybrid_a\" in "
        "the most basic fashion.");
  const int third_offset = n_pts / 3;
  const int small_index = n_pts / 14;
  dhybrid_ptr.setPointer(&dhybrid_a, third_offset);
  section(1);
  check(dhybrid_ptr.readHost(small_index), RelationalOperator::EQUAL,
        dhybrid_a.readHost(third_offset + small_index), "The POINTER-kind Hybrid object "
        "\"dhybrid_ptr\" fails to produce the correct value from its target array after being "
        "given an offset.");
  known_hybrids = gbl_mem_balance_sheet.getEntry("dhybrid_ptr");
  section(3);
  check(known_hybrids[0].length, RelationalOperator::EQUAL, n_pts - third_offset, "The global "
        "memory ledger fails to record the correct length of \"dhybrid_ptr\" when targeting "
        "\"dhybrid_a\" in the most basic fashion.");
  const int half_offset = n_pts / 2;
  const int short_span = n_pts / 6;
  dhybrid_ptr.setPointer(&dhybrid_a, half_offset, short_span);
  known_hybrids = gbl_mem_balance_sheet.getEntry("dhybrid_ptr");
  check(known_hybrids[0].length, RelationalOperator::EQUAL, short_span, "The global memory "
        "ledger fails to record the correct length of \"dhybrid_ptr\" when targeting "
        "\"dhybrid_a\" with an offset and a specified extent (setting the POINTER-kind object in "
        "the most rigorous fashion).");  

  // Try to put things in the POINTER object that would exceed the specified bounds, even though
  // remaining within the actual bounds of the underlying array.
  section(2);
  CHECK_THROWS(dhybrid_ptr.putHost(34.5, short_span), "Hybrid::putHost() was able to add data to "
               "a POINTER kind object that exceeds its specified extent within the underlying "
               "ARRAY.");

  // Try to make an ARRAY-kind Hybrid object target another as if it is a POINTER-kind object
  Hybrid<double> dhybrid_b(n_pts);
  Hybrid<double> dhybrid_ptr_ii(HybridKind::POINTER);
  CHECK_THROWS(dhybrid_a.setPointer(&dhybrid_b), "An ARRAY-kind Hybrid object was set to target "
               "another as if it were a POINTER-kind object.");
  CHECK_THROWS(dhybrid_ptr_ii.setPointer(&dhybrid_ptr, half_offset, short_span), "A POINTER-kind "
               "Hybrid object was set to target another POINTER-kind object.  This would make "
               "tracking difficult and is therefore illegal.");

  // Try to set a POINTER-kind Hybrid in a fashion that would exceed the bounds of the target
  // ARRAY-kind Hybrid.
  CHECK_THROWS(dhybrid_ptr.setPointer(&dhybrid_a, 0, n_pts + half_offset), "Hybrid::setPointer() "
               "was allowed to set a POINTER-kind object with more elements than the target "
               "ARRAY.");

  // Try to create a Hybrid object with something that shouldn't be in a Hybrid object.
  struct MyMess {
    int x;
    double y;
    char z[4];
  };
  CHECK_THROWS(Hybrid<MyMess> zmess_a(n_pts), "A Hybrid object was created from a non-sanctioned "
               "data type.");

  // Attempt some more advanced styles of the putHost method
  const int natom = 305;
  Hybrid<int4> ifrc(natom, "FP forces");
  std::vector<int4> tmp_forces;
  Ran2Generator prng(7183005);
  tmp_forces.resize(natom);
  for (int i = 0; i < natom; i++) {
    int4 iftmp;
    iftmp.x = static_cast<int>(1000000 * prng.gaussianRandomNumber());
    iftmp.y = static_cast<int>(1000000 * prng.gaussianRandomNumber());
    iftmp.z = static_cast<int>(1000000 * prng.gaussianRandomNumber());
    iftmp.w = ((static_cast<int>(1000 * prng.uniformRandomNumber()) << 20) | 
               (static_cast<int>(1000 * prng.uniformRandomNumber()) << 10) | 
               static_cast<int>(1000 * prng.uniformRandomNumber()));
    tmp_forces[i] = iftmp;
  }
  ifrc.putHost(tmp_forces[5], 3);
  int4 irtmp = ifrc.readHost(3);
  check(irtmp.x == tmp_forces[5].x && irtmp.y == tmp_forces[5].y && irtmp.z == tmp_forces[5].z &&
        irtmp.w == tmp_forces[5].w, "Putting an int4 into a Hybrid object and then reading it "
        "back fails to produce the correct number.");
  ifrc.putHost(tmp_forces, 30, 35);
  std::vector<int4> ntmp_forces = ifrc.readHost(30, 35);
  check(int4_series_match(tmp_forces, ntmp_forces, 35), "Putting a series of int4 tuples into a "
        "Hybrid object's host data at a specific offset and then reading it back fails to produce "
        "the correct result.");
  ifrc.putHost(tmp_forces);
  ntmp_forces = ifrc.readHost();
  check(ntmp_forces.size(), RelationalOperator::EQUAL, natom, "The number of elements read from "
        "an int4 Hybrid's entire array is incorrect.");
  check(int4_series_match(tmp_forces, ntmp_forces, natom), "Copying an entire vector of int4 "
        "tuples into a Hybrid object's host data and then reading it back fails to produce the "
        "correct result.");
  irtmp.x = 0;
  irtmp.y = 1;
  irtmp.z = 2;
  irtmp.w = -340593;
  ifrc.putHost(irtmp, 71);
  ntmp_forces = ifrc.readHost();
  check(int4_series_match(tmp_forces, ntmp_forces, natom) == false, "A spot edit of the host data "
        "in a Hybrid object of int4 tuples does not appear to have changed the contents.");
  Hybrid<int4> ifrc_ptr(HybridKind::POINTER, "FP force pointer");
  for (int i = 0; i < natom; i++) {
    std::swap(tmp_forces[i].x, tmp_forces[i].y);
    std::swap(tmp_forces[i].z, tmp_forces[i].w);
  }
  ifrc_ptr.putHost(&ifrc, tmp_forces);
  CHECK_THROWS(ntmp_forces = ifrc_ptr.readHost(10, natom), "The Hybrid object's readHost() "
               "function permits reading a series of elements that would overrun the end of the "
               "array.");
  ntmp_forces = ifrc.readHost();
  check(int4_series_match(tmp_forces, ntmp_forces, natom), "Posting data to a Hybrid object via "
        "the pointer set-and-populate overload of the putHost() method appears to have failed.");

  // Try verifying the POINTER-kind object's target
  section(4);
  check(dhybrid_ptr.verifyTarget(), "POINTER-kind Hybrid object's target did not pass "
        "verification.");

  // Update the length of the target ARRAY-kind object, which should invalidate any POINTER-kind
  // object set to it.
  section(1);
  check(dhybrid_a.size(), RelationalOperator::EQUAL, n_pts, "ARRAY-kind Hybrid object "
        "\"dhybrid_ptr\" displays the wrong length.");
  dhybrid_a.resize(2 * n_pts);
  check(dhybrid_a.getAllocations(), RelationalOperator::EQUAL, 2, "Hybrid object \"dhybrid_a\" "
	"records the wrong number of allocations after resizing.");
  check(dhybrid_a.size(), RelationalOperator::EQUAL, 2 * n_pts, "ARRAY-kind Hybrid object "
        "\"dhybrid_ptr\" displays the wrong length.");
  section(4);
  check(dhybrid_ptr.verifyTarget() == false, "POINTER-kind Hybrid object verified against a stale "
        "target.");

  // Reset the POINTER-kind object with a size commensurate with the new target ARRAY
  dhybrid_ptr.setPointer(&dhybrid_a, half_offset, (2 * n_pts) - half_offset);
  perturb = 0.55;
  delta = 0.01;
  for (int i = n_pts; i < 2 * n_pts; i++) {
    dhybrid_a.putHost(8.0 + perturb + i, i);
    if (perturb > 1.045) {
      delta = -0.02;
    }
    else if (perturb < -1.045) {
      delta = 0.02;
    }
    perturb += delta;
  }
  section(1);
  check(dhybrid_ptr.readHost(n_pts + short_span), RelationalOperator::EQUAL,
        dhybrid_a.readHost(half_offset + n_pts + short_span), "The reset POINTER-kind object "
        "fails to track with the values in its extended target array.");

  // Keep doing pushBack (analog of std::vector push_back() on a Hybrid object) to add some
  // more room to the integer Hybrid object.
  const int orig_allocations = ihybrid_a.getAllocations();
  for (int i = 0; i < 400; i++) {
    ihybrid_a.pushBack(4*i + 7);
  }
  check(ihybrid_a.size(), RelationalOperator::EQUAL, 500, "Hybrid object \"" +
        std::string(ihybrid_a.getLabel().name) + "\" has the wrong size.");
  const int iha_allocations = ihybrid_a.getAllocations();
  check(iha_allocations - orig_allocations <= 5, "Hybrid::pushBack() is generating too "
        "many reallocations (" + std::to_string(ihybrid_a.getAllocations() - orig_allocations) +
        " to grow from 100 to 500 elements).");
  Hybrid<double> grow_by_strides;
  const std::vector<double> real_numbers = { 0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 0.9, -0.8, -1.7 };
  grow_by_strides.pushBack(9.3);
  grow_by_strides.pushBack(real_numbers);
  check(grow_by_strides.size(), RelationalOperator::EQUAL, 10, "The vectorized pushBack() member "
        "function of the Hybrid class does not function as intended.");
  check(grow_by_strides.readHost(3), RelationalOperator::EQUAL, real_numbers[2], "The vectorized "
        "pushBack() member function does not convey the correct data to a Hybrid object.");
  grow_by_strides.pushBack(grow_by_strides);
  check(grow_by_strides.size(), RelationalOperator::EQUAL, 20, "Pushing a Hybrid object to the "
        "back of itself does not work as intended.");
  check(grow_by_strides.readHost(14), RelationalOperator::EQUAL, real_numbers[3], "Pushing a "
        "Hybrid object to the back of itself does not convey data correctly.");
  grow_by_strides.resize(50);
  grow_by_strides.resize(20);
  grow_by_strides.pushBack(grow_by_strides);
  check(grow_by_strides.size(), RelationalOperator::EQUAL, 40, "Resizing a Hybrid object in "
        "advance, such that there is ample room to double its size, results in unexpected "
        "behavior when pushing the object to the back of itself.");
  check(grow_by_strides.readHost(35), RelationalOperator::EQUAL, real_numbers[4], "Pushing a "
        "Hybrid object to the back of itself, when the object does not get reallocated, does not "
        "convey data correctly.");
  
  // Test the copy constructor
  section(5);
  Hybrid<int> icpy_a = ihybrid_a;
  check(icpy_a.getKind() == ihybrid_a.getKind(), "Copying an ARRAY-kind Hybrid object does not "
        "convey the ARRAY-kind characteristic.");
  check(icpy_a.getFormat() == ihybrid_a.getFormat(), "Copying a Hybrid object does not convey the "
        "correct memory format.");
  check(icpy_a.capacity() == ihybrid_a.capacity(), "The copy of a Hybrid object does not possess "
        "the original's maximum capacity.");
  check(icpy_a.getElementSize() == ihybrid_a.getElementSize(), "The copy of a Hybrid object does "
        "not obtain the correct element size from the original.");
  check(std::string(icpy_a.getLabel().name), RelationalOperator::EQUAL,
        std::string(ihybrid_a.getLabel().name), "The copy of a Hybrid object does not inherit the "
        "original object's label as it should.");
  check(icpy_a.getAllocations(), RelationalOperator::EQUAL, 1, "A copy of a Hybrid object should "
        "register one and only one allocation just after creation.");
  check(ihybrid_a.getAllocations(), RelationalOperator::EQUAL, iha_allocations, "The number of "
        "allocations accrued by a Hybrid object should not be changed by making a copy of it.");
  check(icpy_a.size(), RelationalOperator::EQUAL, 500, "Copied Hybrid object \"" +
        std::string(icpy_a.getLabel().name) + " has the wrong size.");
  check(icpy_a.readHost(), RelationalOperator::EQUAL, ihybrid_a.readHost(), "Copying a "
        "Hybrid<int> object fails to duplicate the original data.");
  const int isample_value = icpy_a.readHost(5);
  ihybrid_a.putHost(isample_value + 500123, 5);
  check(icpy_a.readHost(5), RelationalOperator::EQUAL, isample_value, "Copying a Hybrid<int> "
        "object does not decouple the memory from the original.");
  check(icpy_a.getSerialNumber(), RelationalOperator::NOT_EQUAL, ihybrid_a.getSerialNumber(),
        "The copy of a Hybrid<int> object retains the original's serial number.");
  ihybrid_a.resize(50);
  check(icpy_a.size(), RelationalOperator::EQUAL, 500, "The copied Hybrid<int> object is "
        "affected by resizing the original.");
  icpy_a.resize(65);
  check(ihybrid_a.size(), RelationalOperator::EQUAL, 50, "Resizing a copy of Hybrid<int> object " +
        std::string(ihybrid_a.getLabel().name) + " reciprocally affects the original itself.");
  const int jsample_value = ihybrid_a.readHost(2);
  icpy_a.putHost(icpy_a.readHost(2) + 7021893, 2);
  check(ihybrid_a.readHost(2), RelationalOperator::EQUAL, jsample_value, "Manipulating data in "
        "the copy of a Hybrid object affects the original.");

  // Test the copy assignment operator.  Hide the original object that forms the basis for three
  // copies behind a scope to make sure that they collect the right data and come back out of the
  // scope with it.
  Hybrid<double> td_copied, td_copied2, td_copied3;
  if (true) {
    const int hidden_pts = 100;
    Hybrid<double> new_hb(hidden_pts, "Not dead");
    for (int i = 0; i < hidden_pts; i++) {
      new_hb.putHost((3.0 * static_cast<double>(i)) + 1.7, i);
    }
    td_copied = new_hb;
    td_copied2 = new_hb;
    td_copied3 = new_hb;
    new_hb = new_hb;
  }
  check(td_copied.size(), RelationalOperator::EQUAL, 100, "Assignment of a Hybrid object does not "
        "make the object the correct size.");
  check(td_copied.readHost(4), RelationalOperator::EQUAL, Approx(13.7).margin(1.0e-8),
        "The Hybrid assignment operator does not do a deep copy of the values.");
  check(std::string(td_copied.getLabel().name), RelationalOperator::EQUAL, "Not dead",
        "The Hybrid assignment operator does not carry over the proper label.");
  
  // Test the move constructor
  std::vector<Hybrid<int>> stl_vec_of_hybrids;
  std::vector<Hybrid<double>> stl_vec_of_patterns;
  std::vector<int> sum_answer;
  for (int i = 1; i < 10; i++) {
    Hybrid<int> tmp_hyb(i, std::string("tmp_hyb_" + std::to_string(i)).c_str());
    for (int j = 0; j < i; j++) {
      tmp_hyb.putHost(i + 2, j);
    }
    const std::vector<int> tmp_contents = tmp_hyb.readHost();
    sum_answer.push_back(sum<int>(tmp_contents));
    stl_vec_of_hybrids.push_back(tmp_hyb);
    stl_vec_of_patterns.push_back(makeDoubleHybrid((2 * i) - 1));
  }
  std::vector<int> hyb_sums;
  for (size_t i = 0; i < stl_vec_of_hybrids.size(); i++) {
    const std::vector<int> tmp_contents = stl_vec_of_hybrids[i].readHost();
    hyb_sums.push_back(sum<int>(tmp_contents));
  }
  int npatterns = 0;
  for (size_t i = 0; i < stl_vec_of_patterns.size(); i++) {
    npatterns += stl_vec_of_patterns[i].size();
  }
  check(sum_answer, RelationalOperator::EQUAL, hyb_sums, "Hybrid arrays committed to a Standard "
        "Template Library vector via its push_back() member function and the Hybrid copy "
        "constructor were not maintained correctly.");
  const Hybrid<double> result_of_function = makeDoubleHybrid(10);
  check(sum<double>(result_of_function.readHost()), RelationalOperator::EQUAL,
        Approx(6.177750).margin(1.0e-6), "A Hybrid array returned by a function does not contain "
        "the expected data.");
  check(npatterns, RelationalOperator::EQUAL, 81, "Hybrid arrays transferred by the object's move "
        "constructor into a Standard Template Library vector via the push_back() method do not "
        "maintain their proper sizes.");
  std::vector<double> pattern_sums(stl_vec_of_patterns.size());
  for (size_t i = 0; i < stl_vec_of_patterns.size(); i++) {
    pattern_sums[i] = sum<double>(stl_vec_of_patterns[i]);
  }
  const std::vector<double> pattern_sums_answer = {   0.75000000,   8.25000000,  -6.20000000,
                                                     -4.94250000,  24.02250000, -13.66700000,
                                                    -11.15467500,  43.55947500, -18.65477000 };
  check(pattern_sums, RelationalOperator::EQUAL, Approx(pattern_sums_answer).margin(1.0e-6),
        "Hybrid arrays transferred by the object's move constructor into a Standard Template "
        "Library vector via the push_back() method do not retain the correct data.");

  // Test the move assignment operator.
  Hybrid<double> replicator(4);
  for (int i = 0; i < 4; i++) {
    replicator.putHost(static_cast<double>(i) + 0.5, i);
  }
  std::vector<Hybrid<double>> rep_four = { replicator, replicator, replicator, replicator };
  rep_four.erase(rep_four.begin() + 1);
  rep_four.erase(rep_four.begin() + 1);
  check(rep_four.size(), RelationalOperator::EQUAL, 2, "The Standard Template Library vector "
        "erase() method does not work as intended.  This indicates a problem with the Hybrid "
        "object's move assignment operator.");
  
  // Print a summary of tests run
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
