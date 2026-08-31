#include <iostream>
#include <thread>
#include <chrono>
#include "copyright.h"
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

  // Section 1
  section("Test default progress bar functionality");

  // Section 2
  section("Test a modified ASCII progress bar");

  // Section 3
  section("Test the minimal progress bar, percentage only");

  // Initialize progress bar with 10000 iterations
  section(1);
  const int initial_iterations = 5000;
  ProgressBar pbar(initial_iterations);
  pbar.setTerminalWidth(80);
  for (int i = 0; i < initial_iterations; ++i) {
    pbar.update();

    // Simulate work with sleep
    std::this_thread::sleep_for(std::chrono::microseconds(120));
  }
  std::cout << std::endl;
  section(2);

  // Reset and modify progress bar settings
  const int new_iterations = 2500;
  pbar.setCycleCount(new_iterations);
  pbar.reset();
  pbar.setTodoChar(' ');
  pbar.setDoneChar('|');
  pbar.setOpeningBracket("{");
  pbar.setClosingBracket("}");

  for (int i = 0; i < new_iterations; ++i) {
    pbar.update();
    std::this_thread::sleep_for(std::chrono::microseconds(120));
  }
  
  // TODO: Changing this error reporting to STORMM format
  std::cout << std::endl;

  section(3);
  // Disable bar display
  pbar.reset();
  pbar.setStyle(ProgBarStyle::NONE);
  for (int i = 0; i < new_iterations; ++i) {
    pbar.update();
    std::this_thread::sleep_for(std::chrono::microseconds(150));
  }
  // TODO: Changing this error reporting to STORMM format
  std::cout << std::endl;

  // Simulate work with sleep
  for (int i = 0; i < new_iterations; ++i) {
    pbar.update();
    std::this_thread::sleep_for(std::chrono::microseconds(150));
  }
  pbar.update();
  
  // Simulate work
  std::this_thread::sleep_for(std::chrono::milliseconds(800));
  pbar.update();
  
  // Simulate work
  std::this_thread::sleep_for(std::chrono::milliseconds(700));
  std::cout << std::endl;

  // Disable bar display
  section(3);
  pbar.reset();
  pbar.setStyle(ProgBarStyle::NONE);
  for (int i = 0; i < new_iterations; ++i) {
    pbar.update();
    std::this_thread::sleep_for(std::chrono::microseconds(150));
  }

  // TODO: Changing this error reporting to STORMM format
  std::cerr << std::endl;

  // Re-enable bar display and change output stream to std::cout
  // As well as conduct this outside the constraints of a loop
  pbar.setCycleCount(4);
  pbar.reset();
  pbar.setStyle(ProgBarStyle::FULL);
  pbar.update();

  // Simulate work
  std::this_thread::sleep_for(std::chrono::milliseconds(80));
  pbar.update();

  // Simulate work
  std::this_thread::sleep_for(std::chrono::milliseconds(80));
  pbar.update();

  // Simulate work
  std::this_thread::sleep_for(std::chrono::milliseconds(80));
  pbar.update();

  // Simulate work
  std::this_thread::sleep_for(std::chrono::milliseconds(80));

  // Return an accounting of all errors and test results
  printf("\n");
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
