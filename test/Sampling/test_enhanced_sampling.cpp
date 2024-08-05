#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Sampling/temp_distributions.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Topology/atomgraph_analysis.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::constants;
using namespace stormm::diskutil;
using namespace stormm::errors;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::sampling;
using namespace stormm::stmath;
using namespace stormm::symbols;
using namespace stormm::topology;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Van Der Spoel 
//-------------------------------------------------------------------------------------------------
void testVanDerSpoelMethods(){

  // Test cases for different Van Der Spoel algorithm inputs
  // Testing independent variations
  std::vector<double> results = vanDerSpoel(0.2, 293, 393, 15000, 10000, pow(10, -4), 0, 0, 0, 0, 0,
  																					ExceptionResponse::SILENT);
	std::vector<double> actual = {293.00, 294.28, 295.56, 296.85, 298.14, 299.44, 300.74, 302.05,
																303.35, 304.67, 305.99, 307.31, 308.64, 309.97, 311.30, 312.64,
																313.99, 315.33, 316.69, 318.05, 319.41, 320.78, 322.15, 323.52,
																324.91, 326.29, 327.69, 329.09, 330.49, 331.89, 333.30, 334.72,
																336.13, 337.56, 338.99, 340.42, 341.86, 343.31, 344.76, 346.22,
																347.68, 349.14, 350.61, 352.08, 353.56, 355.04, 356.53, 358.02,
																359.52, 361.02, 362.53, 364.04, 365.56, 367.09, 368.62, 370.15,
																371.69, 373.23, 374.78, 376.34, 377.90, 379.46, 381.03, 382.61,
																384.19, 385.78, 387.37, 388.97, 390.58, 392.18, 393.00};
																
  for(const auto& result : results) {
    std::cout << result << " ";
  }
  // std::cout << std::endl;
//   
// 	check(results, RelationalOperator::EQUAL, actual, "The distribution calculation does not match"
// 				" the original");
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
  StopWatch timer;
  Xoshiro256ppGenerator xrs(oe.getRandomSeed());

  // Section 1
  section("Test REMD Distribution Methods");

  // Create a synthesis of systems and the associated particle-mesh interaction grid
  section(1);
  //TODO: Create a test case here
  testVanDerSpoelMethods();
  
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
