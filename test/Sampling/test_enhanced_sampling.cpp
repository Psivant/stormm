#include <string>
#include <vector>
#include <iostream>
#include <iomanip> 
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
#include "../../src/Sampling/replica_probability.h"
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
void testVanDerSpoelMethods() {

  // Test cases for different Van Der Spoel algorithm inputs
  std::vector<double> results = vanDerSpoel(0.2, 293, 393, 15000, 10000, pow(10, -4), 0, 0, 0,
                                            0, 0, ExceptionResponse::SILENT);
  std::vector<double> actual = {293.00, 294.28, 295.56, 296.85, 298.14, 299.44, 300.74, 302.05,
                                303.35, 304.67, 305.99, 307.31, 308.64, 309.97, 311.30, 312.64,
                                313.99, 315.33, 316.69, 318.05, 319.41, 320.78, 322.15, 323.52,
                                324.91, 326.29, 327.69, 329.09, 330.49, 331.89, 333.30, 334.72,
                                336.13, 337.56, 338.99, 340.42, 341.86, 343.31, 344.76, 346.22,
                                347.68, 349.14, 350.61, 352.08, 353.56, 355.04, 356.53, 358.02,
                                359.52, 361.02, 362.53, 364.04, 365.56, 367.09, 368.62, 370.15,
                                371.69, 373.23, 374.78, 376.34, 377.90, 379.46, 381.03, 382.61,
                                384.19, 385.78, 387.37, 388.97, 390.58, 392.18, 393.00};
                                
  check(results, RelationalOperator::EQUAL, Approx(actual).margin(1.0e-1), "Test 1: Error (Base "
        "case with all free atoms, between 293 to 393 Kelvin.)");

  results = vanDerSpoel(0.4, 273, 373, 15000, 10000, 0.0001, 0, 0, 0, 0, 0,
                        ExceptionResponse::SILENT);
  
  actual = {273.00, 273.82, 274.63, 275.45, 276.27, 277.10, 277.92, 278.76, 279.59, 280.42, 281.26,
            282.09, 282.93, 283.78, 284.60, 285.45, 286.30, 287.15, 288.00, 288.86, 289.71, 290.58,
            291.44, 292.31, 293.17, 294.04, 294.91, 295.79, 296.66, 297.54, 298.42, 299.30, 300.18,
            301.07, 301.96, 302.85, 303.75, 304.64, 305.54, 306.44, 307.35, 308.25, 309.16, 310.07,
            310.99, 311.90, 312.82, 313.74, 314.67, 315.59, 316.52, 317.45, 318.38, 319.32, 320.26,
            321.20, 322.14, 323.08, 324.03, 324.98, 325.93, 326.89, 327.85, 328.81, 329.77, 330.73,
            331.70, 332.67, 333.64, 334.62, 335.60, 336.58, 337.56, 338.55, 339.54, 340.53, 341.52,
            342.52, 343.51, 344.52, 345.52, 346.53, 347.54, 348.55, 349.56, 350.58, 351.60, 352.62,
            353.65, 354.68, 355.69, 356.73, 357.76, 358.81, 359.85, 360.89, 361.94, 362.99, 364.04,
            365.10, 366.16, 367.22, 368.28, 369.35, 370.42, 371.49, 372.57, 373.00 };

  check(results, RelationalOperator::EQUAL, Approx(actual).margin(1.0e-1), "Test 2: Error (Higher "
        "probability, all other conditions remain the same.)");

  results = vanDerSpoel(0.3, 273.00, 373.00, 15000, 10000, pow(10, -4), 1, 0, 0, 0, 0,
                        ExceptionResponse::SILENT);
  actual = {273.00, 274.00, 275.00, 276.00, 277.01, 278.02, 279.03, 280.05, 281.07, 282.09,
            283.11, 284.14, 285.19, 286.22, 287.26, 288.30, 289.34, 290.40, 291.45, 292.50,
            293.56, 294.62, 295.68, 296.75, 297.82, 298.89, 299.96, 301.04, 302.13, 303.21,
            304.30, 305.39, 306.49, 307.59, 308.69, 309.80, 310.91, 312.02, 313.13, 314.25,
            315.37, 316.56, 317.69, 318.82, 319.96, 321.09, 322.24, 323.38, 324.53, 325.68,
            326.83, 327.98, 329.15, 330.31, 331.48, 332.66, 333.84, 335.02, 336.20, 337.39,
            338.58, 339.77, 340.97, 342.17, 343.38, 344.59, 345.80, 347.02, 348.24, 349.46,
            350.69, 351.92, 353.16, 354.39, 355.64, 356.88, 358.01, 359.26, 360.52, 361.78,
            363.05, 364.31, 365.59, 366.86, 368.14, 369.43, 370.71, 372.01, 373.00  };
  
  check(results, RelationalOperator::EQUAL, Approx(actual).margin(1.0e-1), "Test 3: Error (Proteins "
        "bond to hydrogen only, everything else is the same.)", TestPriority::NON_CRITICAL);
  for(int i = 0; i < 87; ++i){
    std::cout << "Expected: " << actual[i] << "\t Result: " << results[i] << std::endl;
  }
  std::cout << "Extra element: " << results[87] << std::endl;
  
  results = vanDerSpoel(0.3, 273.00, 373.00, 15000, 10000, pow(10, -4), 0, 1, 0, 0, 0,
                        ExceptionResponse::SILENT);
  actual = {273.00, 274.03, 275.06, 276.10, 277.14, 278.18, 279.22, 280.27, 281.32, 282.38,
            283.43, 284.50, 285.56, 286.63, 287.70, 288.78, 289.86, 290.94, 292.03, 293.12,
            294.21, 295.30, 296.40, 297.50, 298.60, 299.71, 300.83, 301.94, 303.06, 304.19,
            305.31, 306.44, 307.57, 308.71, 309.85, 311.00, 312.15, 313.30, 314.45, 315.61,
            316.77, 317.94, 319.11, 320.28, 321.45, 322.64, 323.82, 325.04, 326.23, 327.42,
            328.62, 329.82, 331.03, 332.24, 333.45, 334.67, 335.89, 337.12, 338.35, 339.58,
            340.82, 342.06, 343.30, 344.55, 345.80, 347.06, 348.32, 349.58, 350.85, 352.12,
            353.40, 354.68, 355.96, 357.25, 358.54, 359.84, 361.14, 362.44, 363.75, 365.03,
            366.35, 367.67, 369.00, 370.32, 371.66, 372.99, 373.00};
  
  check(results, RelationalOperator::EQUAL, Approx(actual).margin(1.0e-1), "Test 4: Error (Water "
        "molecules have flexible angle, proteins are fully flexible.)",
        TestPriority::NON_CRITICAL);

}

//-------------------------------------------------------------------------------------------------
// Test cases for probability checking functions at a replica
//-------------------------------------------------------------------------------------------------
void testReplicaProbabilityMethods() {

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

  section(1);
  //TODO: Create a test case here
  testVanDerSpoelMethods();

  // Section 2
  section("Test REMD Probability functions at a given replica");

  section(2);
  testReplicaProbabilityMethods();  
  
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
