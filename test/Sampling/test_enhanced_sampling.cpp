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
void testVanDerSpoelMethods(StopWatch *timer) {
  const int vds_timings = timer->getCategoryIndex("van der Spoel T distribution");
  timer->assignTime(0);
  
  // Test cases for different Van Der Spoel algorithm inputs
  std::vector<double> results = vanDerSpoel(0.2, 293, 393, 15000, 10000, 1.0e-4, 0, 0, 0,
                                            0, 0, ExceptionResponse::SILENT);
  timer->assignTime(vds_timings);
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
  timer->assignTime(0);

  // Test with exchange probability 0.4
  results = vanDerSpoel(0.4, 273.0, 373.0, 15000, 10000, 1.0e-4, 0, 0, 0, 0, 0,
                        ExceptionResponse::SILENT);
  timer->assignTime(vds_timings);
  actual = { 273.00, 273.82, 274.63, 275.45, 276.27, 277.10, 277.92, 278.76, 279.59, 280.42,
             281.26, 282.09, 282.93, 283.78, 284.60, 285.45, 286.30, 287.15, 288.00, 288.86,
             289.71, 290.58, 291.44, 292.31, 293.17, 294.04, 294.91, 295.79, 296.66, 297.54,
             298.42, 299.30, 300.18, 301.07, 301.96, 302.85, 303.75, 304.64, 305.54, 306.44,
             307.35, 308.25, 309.16, 310.07, 310.99, 311.90, 312.82, 313.74, 314.67, 315.59,
             316.52, 317.45, 318.38, 319.32, 320.26, 321.20, 322.14, 323.08, 324.03, 324.98,
             325.93, 326.89, 327.85, 328.81, 329.77, 330.73, 331.70, 332.67, 333.64, 334.62,
             335.60, 336.58, 337.56, 338.55, 339.54, 340.53, 341.52, 342.52, 343.51, 344.52,
             345.52, 346.53, 347.54, 348.55, 349.56, 350.58, 351.60, 352.62, 353.65, 354.68,
             355.69, 356.73, 357.76, 358.81, 359.85, 360.89, 361.94, 362.99, 364.04, 365.10,
             366.16, 367.22, 368.28, 369.35, 370.42, 371.49, 372.57, 373.00 };
  check(results, RelationalOperator::EQUAL, Approx(actual).margin(1.0e-1), "Test 2: Error (Higher "
        "exchange probability, all other conditions remain the same.)");
  timer->assignTime(0);

  // Test with exchange probability 0.3
  results = vanDerSpoel(0.3, 273.00, 373.00, 15000, 10000, 1.0e-4, 1, 0, 0, 0, 0,
                        ExceptionResponse::SILENT);
  timer->assignTime(vds_timings);
  actual = { 273.00, 273.98, 274.97, 275.96, 276.96, 277.96, 278.96, 279.97, 280.97, 281.98,
             282.99, 284.01, 285.03, 286.05, 287.08, 288.11, 289.14, 290.18, 291.21, 292.25,
             293.30, 294.35, 295.40, 296.45, 297.51, 298.57, 299.63, 300.70, 301.77, 302.84,
             303.91, 304.99, 306.08, 307.16, 308.25, 309.34, 310.44, 311.53, 312.63, 313.74,
             314.84, 315.96, 317.07, 318.19, 319.31, 320.44, 321.56, 322.70, 323.83, 324.97,
             326.11, 327.25, 328.40, 329.55, 330.71, 331.86, 333.03, 334.19, 335.36, 336.53,
             337.71, 338.88, 340.07, 341.25, 342.44, 343.64, 344.83, 346.04, 347.24, 348.45,
             349.66, 350.87, 352.11, 353.33, 354.56, 355.79, 357.02, 358.26, 359.50, 360.74,
             361.99, 363.24, 364.50, 365.76, 367.02, 368.29, 369.56, 370.83, 372.11, 373.00 };
  check(results, RelationalOperator::EQUAL, Approx(actual).margin(1.0e-1), "Test 3: Error "
        "(Proteins bond to hydrogen only, everything else is the same.)");
  timer->assignTime(0);

  // Test with exchange probability 0.3 and ??? TODO ???
  results = vanDerSpoel(0.3, 273.00, 373.00, 15000, 10000, 1.0e-4, 2, 0, 0, 0, 0,
                        ExceptionResponse::SILENT);
  timer->assignTime(vds_timings);
  actual = { 273.00, 273.98, 274.96, 275.95, 276.94, 277.93, 278.93, 279.43, 280.43, 281.43,
             282.44, 283.45, 284.47, 285.48, 286.50, 287.52, 288.55, 289.58, 290.61, 291.64,
             292.68, 293.72, 294.76, 295.82, 296.87, 297.92, 298.98, 300.04, 301.13, 302.20,
             303.27, 304.34, 305.41, 306.49, 307.57, 308.66, 309.75, 310.84, 311.93, 313.03,
             314.13, 315.24, 316.34, 317.45, 318.57, 319.69, 320.81, 321.93, 323.06, 324.19,
             325.32, 326.46, 327.60, 328.75, 329.89, 331.05, 332.20, 333.36, 334.52, 335.68,
             336.85, 338.03, 339.21, 340.39, 341.57, 342.76, 343.95, 345.14, 346.34, 347.54,
             348.74, 349.95, 351.16, 352.38, 353.60, 354.82, 356.04, 357.27, 358.50, 359.74,
             360.98, 362.22, 363.47, 364.72, 365.98, 367.24, 368.50, 369.77, 371.04, 372.31,
             373.00 };
  check(results, RelationalOperator::EQUAL, Approx(actual).margin(1.0e-1), "Test 4: Error "
        "(Proteins bond to all, everything else is the same.)");
  timer->assignTime(0);
  
  results = vanDerSpoel(0.3, 273.00, 373.00, 15000, 10000, 1.0e-4, 0, 1, 0, 0, 0,
                        ExceptionResponse::SILENT);
  timer->assignTime(vds_timings);
  actual = { 273.00, 274.00, 275.02, 276.03, 277.05, 278.08, 279.10, 280.13, 281.16, 282.20,
             283.23, 284.28, 285.32, 286.37, 287.42, 288.47, 289.53, 290.59, 291.65, 292.72,
             293.79, 294.86, 295.94, 297.02, 298.10, 299.19, 300.28, 301.37, 302.47, 303.57,
             304.67, 305.78, 306.89, 308.00, 309.11, 310.23, 311.36, 312.48, 313.61, 314.75,
             315.88, 317.03, 318.18, 319.32, 320.47, 321.63, 322.79, 323.95, 325.12, 326.28,
             327.46, 328.63, 329.81, 330.99, 332.18, 333.37, 334.56, 335.76, 336.96, 338.17,
             339.38, 340.59, 341.80, 343.02, 344.25, 345.47, 346.71, 347.94, 349.18, 350.43,
             351.67, 352.80, 354.05, 355.31, 356.57, 357.83, 359.10, 360.38, 361.65, 362.93,
             364.22, 365.51, 366.80, 368.10, 369.40, 370.70, 372.01, 373.00 };
  check(results, RelationalOperator::EQUAL, Approx(actual).margin(1.0e-1), "Test 5: Error (Water "
        "molecules have flexible angle, proteins are fully flexible.)");
  timer->assignTime(0);
  
  results = vanDerSpoel(0.3, 273.00, 373.00, 15000, 10000, 1.0e-4, 0, 2, 0, 0, 0,
                        ExceptionResponse::SILENT);
  timer->assignTime("van der Spoel T distribution");
  actual = { 273.00, 274.02, 275.06, 276.09, 277.13, 278.17, 279.22, 280.27, 281.32, 282.37,
             283.43, 284.49, 285.56, 286.62, 287.70, 288.77, 289.86, 290.94, 292.02, 293.11,
             294.20, 295.30, 296.39, 297.50, 298.60, 299.71, 300.82, 301.94, 303.06, 304.18,
             305.31, 306.44, 307.57, 308.71, 309.85, 310.99, 312.14, 313.29, 314.45, 315.60,
             316.77, 317.93, 319.10, 320.27, 321.45, 322.63, 323.81, 325.03, 326.22, 327.42,
             328.62, 329.82, 331.03, 332.24, 333.45, 334.67, 335.89, 337.11, 338.34, 339.58,
             340.81, 342.05, 343.30, 344.55, 345.80, 347.05, 348.31, 349.58, 350.85, 352.12,
             353.39, 354.67, 355.96, 357.25, 358.54, 359.83, 361.13, 362.44, 363.75, 365.03,
             366.34, 367.67, 368.99, 370.32, 371.65, 372.99, 373.00 };
  check(results, RelationalOperator::EQUAL, Approx(actual).margin(1.0e-1), "Test 6: Error (Water "
        "molecules have rigid angle, proteins are fully flexible.)");
  timer->assignTime(0);
}

//-------------------------------------------------------------------------------------------------
// Test cases for probability checking functions at a replica
//-------------------------------------------------------------------------------------------------
void testReplicaProbabilityMethods() {

  //TODO: Create a test case here
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
  timer.addCategory("van der Spoel T distribution");

  // Section 1
  section("Test REMD Distribution Methods");

  // Section 2
  section("Test REMD Probability functions at a given replica");

  section(1);
  testVanDerSpoelMethods(&timer);

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
