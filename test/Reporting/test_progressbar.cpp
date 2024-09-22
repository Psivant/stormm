#include <iostream>
#include <thread>
#include <chrono>
#include "../../src/Reporting/progress_bar.h"
#include "../../src/Reporting/present_field.h"
#include "../../src/Reporting/render_options.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::reporting;
using namespace stormm::review;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  
  // Create a ProgressBar object
  ProgressBar progressBar;

  // Section 1
  section("Test default progress bar functionality");

  // Section 2
  section("Test a modified ASCII progress bar");

  // Section 3
  section("Test the minimal progress bar, percentage only");

  // Section 4
  section("Test the progress bar outside a loop");

  section(1);
  
  // Initialize progress bar with 10000 iterations
  const int initial_iterations = 10000;
  progressBar.initialize(initial_iterations);

  // Simulate work with sleep
  for (int i = 0; i < initial_iterations; ++i) {
    progressBar.update();
    std::this_thread::sleep_for(std::chrono::microseconds(150));
  }
  std::cerr << std::endl;
  
  section(2);

  // Reset and modify progress bar settings
  const int new_iterations = 5000;
  progressBar.setIterations(new_iterations);
  progressBar.reset();
  progressBar.setTodoChar(" ");
  progressBar.setDoneChar("█");
  progressBar.setOpeningBracketChar("{");
  progressBar.setClosingBracketChar("}");

  // Simulate work with sleep
  for (int i = 0; i < new_iterations; ++i) {
    progressBar.update();
    std::this_thread::sleep_for(std::chrono::microseconds(150));
  }

  // TODO: Changing this error reporting to STORMM format
  std::cerr << std::endl;

  section(3);

  // Disable bar display
  progressBar.reset();
  progressBar.showBar(false);
  for (int i = 0; i < new_iterations; ++i) {
    progressBar.update();
    std::this_thread::sleep_for(std::chrono::microseconds(150));
  }

  // TODO: Changing this error reporting to STORMM format
  std::cerr << std::endl;
  section(4);

  // Re-enable bar display and change output stream to std::cout
  // As well as conduct this outside the constraints of a loop
  progressBar.setIterations(4);
  progressBar.reset();
  progressBar.showBar(true);
  progressBar.update();

  // Simulate work
  std::this_thread::sleep_for(std::chrono::milliseconds(200));
  progressBar.update();

  // Simulate work
  std::this_thread::sleep_for(std::chrono::milliseconds(200));
  progressBar.update();

  // Simulate work
  std::this_thread::sleep_for(std::chrono::milliseconds(200));
  progressBar.update();

  // Simulate work
  std::this_thread::sleep_for(std::chrono::milliseconds(200));

  // Return an accounting of all errors and test results
  printf("\n");
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
