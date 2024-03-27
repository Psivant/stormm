#include <algorithm>
#include "copyright.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/bspline.h"
#include "../../src/Math/formulas.h"
#include "../../src/Math/math_enumerators.h"
#include "../../src/Math/matrix.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Math/multiplication.h"
#include "../../src/Math/radial_derivatives.h"
#include "../../src/Math/rounding.h"
#include "../../src/Math/series_ops.h"
#include "../../src/Math/sorting_enumerators.h"
#include "../../src/Math/summation.h"
#include "../../src/Math/tickcounter.h"
#include "../../src/Math/tricubic_cell.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/local_arrangement.h"
#include "../../src/Structure/structure_enumerators.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/file_snapshot.h"

using stormm::llint;
using stormm::ulint;
using stormm::ullint;
#ifndef STORMM_USE_HPC
using stormm::int2;
using stormm::double2;
using stormm::double3;
using stormm::double4;
using stormm::float4;
#endif
using stormm::ullint2;
using stormm::constants::tiny;
using stormm::card::Hybrid;
using stormm::card::HybridTargetLevel;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::errors::rtWarn;
using stormm::parse::TextFile;
using stormm::parse::polyNumericVector;
using stormm::random::Ran2Generator;
using stormm::random::Xoroshiro128pGenerator;
using stormm::random::RandomNumberMill;
using stormm::random::Xoshiro256ppGenerator;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using stormm::structure::BoundaryCondition;
using stormm::structure::imageValue;
using stormm::structure::ImagingMethod;
using stormm::symbols::pi;
using stormm::symbols::twopi;
using stormm::topology::getEnumerationName;
using stormm::topology::UnitCellType;
using namespace stormm::stmath;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Check the matrix multiplication routine with a hand-written, much simpler protocol.
//
// Arguments:
//   row_a:  Row count of matrix A
//   col_a:  Column count of matrix A
//   row_b:  Row count of matrix B
//   col_b:  Column count of matrix B
//   axps:   Transposition state for matrix A
//   bxps:   Transposition state for matrix B
//-------------------------------------------------------------------------------------------------
void checkMatrixMultiply(const int row_a, const int col_a, const int row_b, const int col_b,
                         const TransposeState axps = TransposeState::AS_IS,
                         const TransposeState bxps = TransposeState::AS_IS) {

  // Try a matrix multiply between 5 x 5 and 5 x 2 matrices
  std::vector<double> amat(row_a * col_a), bmat(row_b * col_b);
  Xoroshiro128pGenerator xrsm((row_a * col_b) + (row_b * col_a) +
                              ((row_a < row_b) * (row_a * row_b)) + (col_a * col_b)); 
  for (int i = 0; i < row_a * col_a; i++) {
    amat[i] = xrsm.uniformRandomNumber();
  }
  for (int i = 0; i < row_b * col_b; i++) {
    bmat[i] = xrsm.uniformRandomNumber();
  }
  int actual_row_a, actual_col_a, actual_row_b, actual_col_b;
  switch (axps) {
  case TransposeState::AS_IS:
    actual_row_a = row_a;
    actual_col_a = col_a;
    break;
  case TransposeState::TRANSPOSE:
    actual_row_a = col_a;
    actual_col_a = row_a;
    break;
  }
  switch (bxps) {
  case TransposeState::AS_IS:
    actual_row_b = row_b;
    actual_col_b = col_b;
    break;
  case TransposeState::TRANSPOSE:
    actual_row_b = col_b;
    actual_col_b = row_b;
    break;
  }
  std::vector<double> result(actual_row_a * actual_col_b, 0.0);
  std::vector<double> byhand(actual_row_a * actual_col_b, 0.0);
  matrixMultiply(amat.data(), row_a, col_a, bmat.data(), row_b, col_b, result.data(), 1.0, 1.0,
                 0.0, axps, bxps);
  if (axps == TransposeState::TRANSPOSE) {
    transpose(&amat, row_a, col_a);
  }
  if (bxps == TransposeState::TRANSPOSE) {
    transpose(&bmat, row_b, col_b);
  }
  for (int i = 0; i < actual_row_a; i++) {
    for (int j = 0; j < actual_col_a; j++) {
      for (int k = 0; k < actual_col_b; k++) {
        byhand[(k * actual_row_a) + i] += amat[(j * actual_row_a) + i] *
                                          bmat[(k * actual_row_b) + j];
      }
    }
  }
  const std::string a_state = (axps == TransposeState::TRANSPOSE) ? "(t)" : "";
  const std::string b_state = (bxps == TransposeState::TRANSPOSE) ? "(t)" : "";
  check(result, RelationalOperator::EQUAL, byhand, "A(" + std::to_string(row_a) + " x " +
        std::to_string(col_a) + ")" + a_state + " * B(" + std::to_string(row_b) + " x " +
        std::to_string(col_b) + ")" + b_state + " was not multiplied correctly.");
}

//-------------------------------------------------------------------------------------------------
// Take the (Moore-Penrose) pseudo-inverse of a matrix.
//
// Arguments:
//   nrow:  Number of rows to impart to the matrix
//   ncol:  Number of columns to impart to the matrix
//-------------------------------------------------------------------------------------------------
void takePseudoInverse(const int nrow, const int ncol, const std::string &matrices_snp,
                       const std::string &var_name, const TestEnvironment &oe,
                       const TestPriority do_test) {
  std::vector<double> amat(nrow * ncol);
  Xoroshiro128pGenerator xrsm(nrow * ncol);
  for (int i = 0; i < nrow * ncol; i++) {
    amat[i] = xrsm.gaussianRandomNumber() + 1.5;
  }
  pseudoInverse(&amat, nrow, ncol);
  snapshot(matrices_snp, polyNumericVector(amat), var_name, 1.0e-5, "The matrix pseudo-inversion "
           "result for a " + std::to_string(nrow) + " x " + std::to_string(ncol) + " matrix is "
           "incorrect.", oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL,
           PrintSituation::APPEND, do_test);
}

//-------------------------------------------------------------------------------------------------
// Fashion an over-determined system of linear equations based on random coefficients, a
// pre-determined solution vector, and random perturbations.  The perturbations guarantee that,
// while the pre-determined solution is still a good one, it is not the optimal answer.  The
// optimal answer can then be determined using either the Moore-Penrose left pseudo-inverse or a
// QR decomposition.  Both methods should give very similar answers.
//
// Arguments:
//   n_equations:  The number of equations in the system
//   n_unknowns:   The number of independent coefficients describing the cost function
//-------------------------------------------------------------------------------------------------
void solveLinearSystem(const size_t n_equations, const size_t n_unknowns) {
  Xoroshiro128pGenerator xrsm(n_equations * n_unknowns);
  const std::vector<double> base_coefficients = gaussianRand(&xrsm, n_unknowns, 1.0);
  const std::vector<double> amat = gaussianRand(&xrsm, n_equations * n_unknowns, 1.0);
  std::vector<double> bvec = gaussianRand(&xrsm, n_equations, 1.0);
  elementwiseMultiply(&bvec, 0.5);
  matrixVectorMultiply(amat, base_coefficients, &bvec, n_equations, n_unknowns, 1.0, 1.0, 1.0);

  // Compute the solution using the pseudo-inverse
  const std::vector<double> amat_pinv = pseudoInverse(amat, n_equations, n_unknowns);
  std::vector<double> xvec_pinv(n_unknowns, 0.0);
  matrixVectorMultiply(amat_pinv, bvec, &xvec_pinv, n_unknowns, n_equations);
  std::vector<double> base_result(n_equations, 0.0);
  matrixVectorMultiply(amat, base_coefficients, &base_result, n_equations, n_unknowns);
  std::vector<double> pinv_result(n_equations, 0.0);
  matrixVectorMultiply(amat, xvec_pinv, &pinv_result, n_equations, n_unknowns);
  const double base_error = relativeRmsError(base_result, bvec);
  const double pinv_error = relativeRmsError(pinv_result, bvec);
  check(base_error > pinv_error, "The pseudo-inverse does not produce a better answer than the "
        "original solution before the introduction of noise.");

  // Compute the solution using the QR decomposition
  std::vector<double> xvec_qr(n_unknowns, 0.0);
  qrSolver(amat, &xvec_qr, bvec);

  // Check that the two solutions are similar.
  check(xvec_qr, RelationalOperator::EQUAL, xvec_pinv, "The coefficients computed by "
        "Moore-Penrose pseudo-inverse and QR decomposition do not agree.");
}

//-------------------------------------------------------------------------------------------------
// Produce a series of random numbers from a generator based on an initial state and some
// directions about how to cycle the generator after each sample.
//
// Arguments:
//   xrs:         Random number generator (the templating revolves around this type and its state)
//   init_state:  Initial state of the first generator
//   gen_idx:     Generator index
//   ngen:        Number of random generators
//   depth:       Depth of the pre-computed random numbers cache
//   maxcyc:      Number of cycles to run the generator through after each sample is taken
//-------------------------------------------------------------------------------------------------
template <typename Tgen, typename Tstate>
std::vector<double> generateExpectedSeries(Tgen *xrs, const Tstate init_state, const int gen_idx,
                                           const int ngen, const int depth, const int maxcyc) {

  // Initialize the pseudo-random number generator
  xrs->setState(init_state);
  for (int i = 0; i < gen_idx; i++) {
    xrs->longJump();
  }
  
  // Create a buffer to represent the series that should be held within the generator series
  std::vector<double> buffer(depth);
  for (int i = 0; i < depth; i++) {
    buffer[i] = xrs->gaussianRandomNumber();
  }

  // Allocate the result
  std::vector<double> result(maxcyc * depth);
  
  // Compute the refresh schedule
  const int rstride = (ngen + depth - 1) / depth;
  const int rsched  = gen_idx / rstride;
  for (int i = 0; i < maxcyc * depth; i++) {

    // Pull a value from the table
    const int table_idx = i - ((i / depth) * depth);
    result[i] = buffer[table_idx];
    
    // Refresh the table as appropriate
    if (table_idx == rsched) {
      for (int j = 0; j < depth; j++) {
        buffer[j] = xrs->gaussianRandomNumber();
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Test a TickCounter object to ensure that its incrementation is valid for a variety of series
// lengths.
//
// Arguments:
//   dials:  The tick counter to test, "a collection of dials"
//   ncyc:   The number of cycles to increment forward or backward
//-------------------------------------------------------------------------------------------------
void testDialCounter(TickCounter<int> *dials, const int ncyc) {
  const int ndials = dials->getVariableCount();
  const std::vector<int>& dial_settings = dials->getSettings();
  const std::vector<int>& dial_limits = dials->getStateLimits();
  for (int i = 0; i < ncyc; i++) {
    dials->advance();
  }
  int nperm = 1;
  for (int i = 0; i < ndials; i++) {
    nperm *= dial_limits[i];
  }
  std::vector<int> expected_settings(ndials);
  int balance = (ncyc % nperm);
  int tick_value = nperm / dial_limits[ndials - 1];
  for (int i = ndials - 1; i >= 0; i--) {
    expected_settings[i] = (balance / tick_value);
    balance -= expected_settings[i] * tick_value;
    if (i > 0) {
      tick_value /= dial_limits[i - 1];
    }
  }
  check(dial_settings, RelationalOperator::EQUAL, expected_settings, "A TickCounter with " +
        std::to_string(ndials) + " variables did not increment as expected.");
  const std::vector<int> inc_result = dials->getSettings();
  dials->reset();
  dials->advance(ncyc);
  check(dial_settings, RelationalOperator::EQUAL, expected_settings, "A TickCounter with " +
        std::to_string(ndials) + " variables did not jump as expected.");
  for (int i = 0; i < 50; i++) {
    dials->reverse();
  }
  dials->reverse(ncyc / 7);
  dials->advance((ncyc / 7) + 50);
  check(dial_settings, RelationalOperator::EQUAL, expected_settings, "A TickCounter with " +
        std::to_string(ndials) + " variables did not jump as expected.");
}

//-------------------------------------------------------------------------------------------------
// Evaluate different polynomial functions of three dimensions at a specific point in space.  The
// functions in question include a separable example, evaluated by default, which can be stated:
//
//      F(x, y, z) = (5x^3 + 2x^2 + x + 3) * (y^3 - 4y^2 + 7y - 8) * (6z^3 - 9z + 1)
//
// There is also the option to include a list of coefficients for all combinations of powers of
// the variables x, y, and z.  If supplied, this form will be evaluated instead.
//
// Arguments:
//   x:     The Cartesian X coordinate
//   y:     The Cartesian Y coordinate
//   z:     The Cartesian Z coordinate
//   kind:  The type of function derivative (or the primary value) to return
//-------------------------------------------------------------------------------------------------
double triPolynomial(const double x, const double y, const double z, const FunctionLevel kind,
                     const std::vector<double> coeffs = {}) {
  const bool do_separable_case = (coeffs.size() != 64LLU);
  int imin = 0;
  int jmin = 0;
  int kmin = 0;
  switch (kind) {
  case FunctionLevel::VALUE:
    if (do_separable_case) {

      // Horner's rule would make these computations more efficient, but write expressions in a
      // format that uses fewer parentheses for clarity.
      return ((5.0 * x * x * x) + ( 2.0 * x * x) + ( 1.0 * x) +  3.0) *
             ((1.0 * y * y * y) - ( 4.0 * y * y) + ( 7.0 * y) -  8.0) *
             ((6.0 * z * z * z)                  - ( 9.0 * z) +  1.0);
    }
    break;
  case FunctionLevel::DX:
    if (do_separable_case) {
      return (                    (15.0 * x * x) + ( 4.0 * x) +  1.0) *
             ((1.0 * y * y * y) - ( 4.0 * y * y) + ( 7.0 * y) -  8.0) *
             ((6.0 * z * z * z)                  - ( 9.0 * z) +  1.0);
    }
    else {
      imin = 1;
    }
    break;
  case FunctionLevel::DY:
    if (do_separable_case) {
      return ((5.0 * x * x * x) + ( 2.0 * x * x) + ( 1.0 * x) +  3.0) *
             (                    ( 3.0 * y * y) - ( 8.0 * y) +  7.0) *
             ((6.0 * z * z * z)                  - ( 9.0 * z) +  1.0);
    }
    else {
      jmin = 1;
    }
    break;
  case FunctionLevel::DZ:
    if (do_separable_case) {
      return ((5.0 * x * x * x) + ( 2.0 * x * x) + ( 1.0 * x) +  3.0) *
             ((1.0 * y * y * y) - ( 4.0 * y * y) + ( 7.0 * y) -  8.0) *
             (                    (18.0 * z * z)              -  9.0);
    }
    else {
      kmin = 1;
    }
    break;
  case FunctionLevel::DXX:
    if (do_separable_case) {
      return (                                     (30.0 * x) +  4.0) *
             ((1.0 * y * y * y) - ( 4.0 * y * y) + ( 7.0 * y) -  8.0) *
             ((6.0 * z * z * z)                  - ( 9.0 * z) +  1.0);
    }
    else {
      imin = 2;
    }
    break;
  case FunctionLevel::DXY:
    if (do_separable_case) {
      return (                    (15.0 * x * x) + ( 4.0 * x) +  1.0) *
             (                    ( 3.0 * y * y) - ( 8.0 * y) +  7.0) *
             ((6.0 * z * z * z)                  - ( 9.0 * z) +  1.0);
    }
    else {
      imin = 1;
      jmin = 1;
    }
    break;
  case FunctionLevel::DXZ:
    if (do_separable_case) {
      return (                    (15.0 * x * x) + ( 4.0 * x) +  1.0) *
             ((1.0 * y * y * y) - ( 4.0 * y * y) + ( 7.0 * y) -  8.0) *
             (                    (18.0 * z * z)              -  9.0);
    }
    else {
      imin = 1;
      kmin = 1;
    }
    break;
  case FunctionLevel::DYY:
    if (do_separable_case) {
      return ((5.0 * x * x * x) + ( 2.0 * x * x) + ( 1.0 * x) +  3.0) *
             (                                     ( 6.0 * y) -  8.0) *
             ((6.0 * z * z * z)                  - ( 9.0 * z) +  1.0);
    }
    else {
      jmin = 2;
    }
    break;
  case FunctionLevel::DYZ:
    if (do_separable_case) {
      return ((5.0 * x * x * x) + ( 2.0 * x * x) + ( 1.0 * x) +  3.0) *
             (                    ( 3.0 * y * y) - ( 8.0 * y) +  7.0) *
             (                    (18.0 * z * z)              -  9.0);
    }
    else {
      jmin = 1;
      kmin = 1;
    }
    break;
  case FunctionLevel::DZZ:
    if (do_separable_case) {
      return ((5.0 * x * x * x) + ( 2.0 * x * x) + ( 1.0 * x) +  3.0) *
             ((1.0 * y * y * y) - ( 4.0 * y * y) + ( 7.0 * y) -  8.0) *
             (                                     (36.0 * z)       );
    }
    else {
      kmin = 2;
    }
    break;
  case FunctionLevel::DXXX:
    if (do_separable_case) {
      return (                                                  30.0) *
             ((1.0 * y * y * y) - ( 4.0 * y * y) + ( 7.0 * y) -  8.0) *
             ((6.0 * z * z * z)                  - ( 9.0 * z) +  1.0);
    }
    else {
      imin = 3;
    }
    break;
  case FunctionLevel::DXXY:
    if (do_separable_case) {
      return (                                     (30.0 * x) +  4.0) *
             (                    ( 3.0 * y * y) - ( 8.0 * y) +  7.0) *
             ((6.0 * z * z * z)                  - ( 9.0 * z) +  1.0);
    }
    else {
      imin = 2;
      jmin = 1;
    }
    break;
  case FunctionLevel::DXXZ:
    if (do_separable_case) {
      return (                                     (30.0 * x) +  4.0) *
             ((1.0 * y * y * y) - ( 4.0 * y * y) + ( 7.0 * y) -  8.0) *
             (                    (18.0 * z * z)              -  9.0);
    }
    else {
      imin = 2;
      kmin = 1;
    }
    break;
  case FunctionLevel::DXYY:
    if (do_separable_case) {
      return (                    (15.0 * x * x) + ( 4.0 * x) +  1.0) *
             (                                     ( 6.0 * y) -  8.0) *
             ((6.0 * z * z * z)                  - ( 9.0 * z) +  1.0);
    }
    else {
      imin = 1;
      jmin = 2;
    }
    break;
  case FunctionLevel::DXYZ:
    if (do_separable_case) {
      return (                    (15.0 * x * x) + ( 4.0 * x) +  1.0) *
             (                    ( 3.0 * y * y) - ( 8.0 * y) +  7.0) *
             (                    (18.0 * z * z)              -  9.0);
    }
    else {
      imin = 1;
      jmin = 1;
      kmin = 1;
    }
    break;
  case FunctionLevel::DXZZ:
    if (do_separable_case) {
      return (                    (15.0 * x * x) + ( 4.0 * x) +  1.0) *
             ((1.0 * y * y * y) - ( 4.0 * y * y) + ( 7.0 * y) -  8.0) *
             (                                     (36.0 * z)       );
    }
    else {
      imin = 1;
      kmin = 2;
    }
    break;
  case FunctionLevel::DYYY:
    if (do_separable_case) {
      return ((5.0 * x * x * x) + ( 2.0 * x * x) + ( 1.0 * x) +  3.0) *
             (                                                   6.0) *
             ((6.0 * z * z * z)                  - ( 9.0 * z) +  1.0);
    }
    else {
      jmin = 3;
    }
    break;
  case FunctionLevel::DYYZ:
    if (do_separable_case) {
      return ((5.0 * x * x * x) + ( 2.0 * x * x) + ( 1.0 * x) +  3.0) *
             (                                     ( 6.0 * y) -  8.0) *
             (                    (18.0 * z * z)              -  9.0);
    }
    else {
      jmin = 2;
      kmin = 1;
    }
    break;
  case FunctionLevel::DYZZ:
    if (do_separable_case) {
      return ((5.0 * x * x * x) + ( 2.0 * x * x) + ( 1.0 * x) +  3.0) *
             (                    ( 3.0 * y * y) - ( 8.0 * y) +  7.0) *
             (                                     (36.0 * z)       );
    }
    else {
      jmin = 1;
      kmin = 2;
    }
    break;
  case FunctionLevel::DZZZ:
    if (do_separable_case) {
      return ((5.0 * x * x * x) + ( 2.0 * x * x) + ( 1.0 * x) +  3.0) *
             ((1.0 * y * y * y) - ( 4.0 * y * y) + ( 7.0 * y) -  8.0) *
             (                                                  36.0);
    }
    else {
      kmin = 3;
    }
    break;
  }

  // If this point is reached, the function is evaluating the general case with a vector of
  // input coefficients.
  double dfactors[4][4] = { { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 1.0, 2.0, 3.0 },
                            { 0.0, 0.0, 2.0, 6.0 }, { 0.0, 0.0, 0.0, 6.0 } };
  double result = 0.0;
  double xv = 1.0;
  for (int i = imin; i < 4; i++) {
    double yv = 1.0;
    for (int j = jmin; j < 4; j++) {
      double zv = 1.0;
      for (int k = kmin; k < 4; k++) {
        result += coeffs[(4 * ((4 * k) + j)) + i] *
                  dfactors[imin][i] * xv * dfactors[jmin][j] * yv * dfactors[kmin][k] * zv;
        zv *= z;
      }
      yv *= y;
    }
    xv *= x;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Use the finite difference method to check that the derivatives for the triPolynomial function
// are being computed correctly.  Descriptions of the formal arguments follow from triPolynomial()
// above.
//-------------------------------------------------------------------------------------------------
void finiteDifferenceDisparity(const double x, const double y, const double z,
                               const std::vector<double> coeffs = {}) {
  const double incr = 0.000005;
  const double step = 2.0 * incr;
  const double v = triPolynomial(x, y, z, FunctionLevel::VALUE, coeffs);
  const double vpx = triPolynomial(x + incr, y, z, FunctionLevel::VALUE, coeffs);
  const double vpy = triPolynomial(x, y + incr, z, FunctionLevel::VALUE, coeffs);
  const double vpz = triPolynomial(x, y, z + incr, FunctionLevel::VALUE, coeffs);
  const double vpxy = triPolynomial(x, y + incr, z, FunctionLevel::DX, coeffs);
  const double vpxz = triPolynomial(x, y, z + incr, FunctionLevel::DX, coeffs);
  const double vpyz = triPolynomial(x, y, z + incr, FunctionLevel::DY, coeffs);
  const double vpxyz = triPolynomial(x, y, z + incr, FunctionLevel::DXY, coeffs);
  const double vnx = triPolynomial(x - incr, y, z, FunctionLevel::VALUE, coeffs);
  const double vny = triPolynomial(x, y - incr, z, FunctionLevel::VALUE, coeffs);
  const double vnz = triPolynomial(x, y, z - incr, FunctionLevel::VALUE, coeffs);
  const double vnxy = triPolynomial(x, y - incr, z, FunctionLevel::DX, coeffs);
  const double vnxz = triPolynomial(x, y, z - incr, FunctionLevel::DX, coeffs);
  const double vnyz = triPolynomial(x, y, z - incr, FunctionLevel::DY, coeffs);
  const double vnxyz = triPolynomial(x, y, z - incr, FunctionLevel::DXY, coeffs);
  const std::vector<double> fdv = { v, (vpx - vnx) / step, (vpy - vny) / step, (vpz - vnz) / step,
                                    (vpxy - vnxy) / step, (vpxz - vnxz) / step,
                                    (vpyz - vnyz) / step, (vpxyz - vnxyz) / step };
  const std::vector<double> dv = { triPolynomial(x, y, z, FunctionLevel::VALUE, coeffs),
                                   triPolynomial(x, y, z, FunctionLevel::DX, coeffs),
                                   triPolynomial(x, y, z, FunctionLevel::DY, coeffs),
                                   triPolynomial(x, y, z, FunctionLevel::DZ, coeffs),
                                   triPolynomial(x, y, z, FunctionLevel::DXY, coeffs),
                                   triPolynomial(x, y, z, FunctionLevel::DXZ, coeffs),
                                   triPolynomial(x, y, z, FunctionLevel::DYZ, coeffs),
                                   triPolynomial(x, y, z, FunctionLevel::DXYZ, coeffs) };
  check(dv, RelationalOperator::EQUAL, Approx(fdv).margin(1.0e-6), "Analytic derivatives of the "
        "tricubic polynomial do not agree with results computed by the finite difference "
        "approximation.  Subsequent tests of the TricubicCell object may be unreliable.");
  const std::vector<double> fdv2 = { v, (vpx - vnx) / step, (vpy - vny) / step, (vpz - vnz) / step,
                                     (triPolynomial(x + incr, y, z, FunctionLevel::DY, coeffs) -
                                      triPolynomial(x - incr, y, z, FunctionLevel::DY, coeffs)) /
                                     step,
                                     (triPolynomial(x + incr, y, z, FunctionLevel::DZ, coeffs) -
                                      triPolynomial(x - incr, y, z, FunctionLevel::DZ, coeffs)) /
                                     step,
                                     (triPolynomial(x, y + incr, z, FunctionLevel::DZ, coeffs) -
                                      triPolynomial(x, y - incr, z, FunctionLevel::DZ, coeffs)) /
                                     step,
                                     (triPolynomial(x, y + incr, z, FunctionLevel::DXZ, coeffs) -
                                      triPolynomial(x, y - incr, z, FunctionLevel::DXZ, coeffs)) /
                                     step };
  check(fdv, RelationalOperator::EQUAL, fdv2, "Mixed partial derivatives of the tricubic "
        "polynomial do not show equivalence.  Subsequent tests of the TricubicCell object may be "
        "unreliable.");
}

//-------------------------------------------------------------------------------------------------
// Run a series of tests on a grid of TricubicCell objects.
//
// Arguments:
//   tcmat:               Standard Template Library copy of the tricubic spline weights matrix
//   coeffs:              [Optional] Vector of coefficients for individual terms in the polynomial
//   cell_dimensions_in:  [Optional] Dimensions of the interpolating volume element
//-------------------------------------------------------------------------------------------------
void tricubicTestBundle(const TricubicStencil &tcstencil,
                        const std::vector<double> coeffs = {},
                        const std::vector<double> &cell_dimensions_in = {}) {

  // Prepare the cell dimensions
  std::vector<double> cell_dimensions(9, 0.0), inv_cell_dimensions(9);
  cell_dimensions[0] = 1.0;
  cell_dimensions[4] = 1.0;
  cell_dimensions[8] = 1.0;
  bool orthorhombic = true;
  if (cell_dimensions_in.size() == 3LLU) {
    cell_dimensions[0] = cell_dimensions_in[0];
    cell_dimensions[4] = cell_dimensions_in[1];
    cell_dimensions[8] = cell_dimensions_in[2];
  }
  else if (cell_dimensions_in.size() == 9LLU) {
    cell_dimensions = cell_dimensions_in;
    orthorhombic = false;
  }
  invertSquareMatrix(cell_dimensions, &inv_cell_dimensions);
  
  // Begin by testing the polynomial's analytic evaluation, to ensure subsequent tests are valid.
  finiteDifferenceDisparity( 0.5,  0.9,  0.4, coeffs);
  finiteDifferenceDisparity(-1.7,  2.6, -1.8, coeffs);
  finiteDifferenceDisparity( 2.3,  1.1, -2.0, coeffs);
  std::vector<TricubicCell<double>> grid;
  grid.reserve(1331);
  for (int k = -5; k < 6; k++) {
    const double d_k = k;
    for (int j = -5; j < 6; j++) {
      const double d_j = j;
      for (int i = -5; i < 6; i++) {
        const double d_i = i;

        // Map the function U = (5x^3 + 2x^2 + x + 3) * (y^3 - 4y^2 + 7y - 8) * (6z^3 - 9z + 1)
        // to a grid.  The function is encoded in triPolynomial() above.
        std::vector<double> f(8), dx(8), dy(8), dz(8);
        std::vector<double> dxx(8), dxy(8), dxz(8), dyy(8), dyz(8);
        std::vector<double> dxxx(8), dxxy(8), dxxz(8), dxyy(8), dxyz(8);
        for (int pi = 0; pi < 2; pi++) {
          const double d_ipi = d_i + static_cast<double>(pi);
          for (int pj = 0; pj < 2; pj++) {
            const double d_jpj = d_j + static_cast<double>(pj);
            for (int pk = 0; pk < 2; pk++) {
              const double d_kpk = d_k + static_cast<double>(pk);
              const double xloc = (cell_dimensions[0] * d_ipi) + (cell_dimensions[3] * d_jpj) +
                                  (cell_dimensions[6] * d_kpk);
              const double yloc = (cell_dimensions[4] * d_jpj) + (cell_dimensions[7] * d_kpk);
              const double zloc = (cell_dimensions[8] * d_kpk);
              const int m = (2 * ((2 * pk) + pj)) + pi;
              f[m]    = triPolynomial(xloc, yloc, zloc, FunctionLevel::VALUE, coeffs);
              dx[m]   = triPolynomial(xloc, yloc, zloc, FunctionLevel::DX, coeffs);
              dy[m]   = triPolynomial(xloc, yloc, zloc, FunctionLevel::DY, coeffs);
              dz[m]   = triPolynomial(xloc, yloc, zloc, FunctionLevel::DZ, coeffs);
              switch (tcstencil.getKind()) {
              case Interpolant::SMOOTHNESS:
                dxx[m]  = triPolynomial(xloc, yloc, zloc, FunctionLevel::DXX, coeffs);
                dxy[m]  = triPolynomial(xloc, yloc, zloc, FunctionLevel::DXY, coeffs);
                dxz[m]  = triPolynomial(xloc, yloc, zloc, FunctionLevel::DXZ, coeffs);
                dyy[m]  = triPolynomial(xloc, yloc, zloc, FunctionLevel::DYY, coeffs);
                dyz[m]  = triPolynomial(xloc, yloc, zloc, FunctionLevel::DYZ, coeffs);
                dxxx[m] = triPolynomial(xloc, yloc, zloc, FunctionLevel::DXXX, coeffs);
                dxxy[m] = triPolynomial(xloc, yloc, zloc, FunctionLevel::DXXY, coeffs);
                dxxz[m] = triPolynomial(xloc, yloc, zloc, FunctionLevel::DXXZ, coeffs);
                dxyy[m] = triPolynomial(xloc, yloc, zloc, FunctionLevel::DXYY, coeffs);
                dxyz[m] = triPolynomial(xloc, yloc, zloc, FunctionLevel::DXYZ, coeffs);
                break;
              case Interpolant::FUNCTION_VALUE:
                {
                  const double eorig_x = (cell_dimensions[0] * d_i) + (cell_dimensions[3] * d_j) +
                                         (cell_dimensions[6] * d_k);
                  const double eorig_y = (cell_dimensions[4] * d_j) + (cell_dimensions[7] * d_k);
                  const double eorig_z = (cell_dimensions[8] * d_k);
                  double abloc_x, abloc_y, abloc_z, acloc_x, acloc_y, acloc_z;
                  double bcloc_x, bcloc_y, bcloc_z, cnloc_x, cnloc_y, cnloc_z;
                  fvStencilCoordinates(eorig_x, eorig_y, eorig_z, UnitCellAxis::C, pi, pj, pk,
                                       cell_dimensions.data(), &abloc_x, &abloc_y, &abloc_z);
                  fvStencilCoordinates(eorig_x, eorig_y, eorig_z, UnitCellAxis::B, pi, pj, pk,
                                       cell_dimensions.data(), &acloc_x, &acloc_y, &acloc_z);
                  fvStencilCoordinates(eorig_x, eorig_y, eorig_z, UnitCellAxis::A, pi, pj, pk,
                                       cell_dimensions.data(), &bcloc_x, &bcloc_y, &bcloc_z);
                  fvStencilCoordinates(eorig_x, eorig_y, eorig_z, pi, pj, pk,
                                       cell_dimensions.data(), &cnloc_x, &cnloc_y, &cnloc_z);
                  dxy[m] = triPolynomial(abloc_x, abloc_y, abloc_z, FunctionLevel::VALUE, coeffs);
                  dxz[m] = triPolynomial(acloc_x, acloc_y, acloc_z, FunctionLevel::VALUE, coeffs);
                  dyz[m] = triPolynomial(bcloc_x, bcloc_y, bcloc_z, FunctionLevel::VALUE, coeffs);
                  dxyz[m] = triPolynomial(cnloc_x, cnloc_y, cnloc_z, FunctionLevel::VALUE, coeffs);
                }
                break;
              }
            }
          }
        }
        const double xe_orig = (cell_dimensions[0] * d_i) + (cell_dimensions[3] * d_j) +
                               (cell_dimensions[6] * d_k);
        const double ye_orig = (cell_dimensions[4] * d_j) + (cell_dimensions[7] * d_k);
        const double ze_orig = (cell_dimensions[8] * d_k);
        const std::vector<double> anchor = { xe_orig, ye_orig, ze_orig, cell_dimensions[0],
                                             cell_dimensions[1], cell_dimensions[2],
                                             cell_dimensions[3], cell_dimensions[4],
                                             cell_dimensions[5], cell_dimensions[6],
                                             cell_dimensions[7], cell_dimensions[8] };
        switch (tcstencil.getKind()) {
        case Interpolant::SMOOTHNESS:
          grid.emplace_back(tcstencil, anchor, f, dx, dy, dz, dxx, dxy, dxz, dyy, dyz, dxxx, dxxy,
                            dxxz, dxyy, dxyz);
          break;
        case Interpolant::FUNCTION_VALUE:
          grid.emplace_back(tcstencil, anchor, f, dx, dy, dz, dxy, dxz, dyz, dxyz);
          break;
        }
      }
    }
  }
  const std::vector<double> test_points = {  0.10,  0.50, -0.40,
                                            -0.80,  1.20, -1.30,
                                             1.80,  2.20, -0.70,
                                             0.01,  0.01,  0.01,
                                             1.00,  1.00,  1.00,
                                             1.01,  1.01,  1.01,
                                             1.50,  1.20,  1.40,
                                            -2.00,  1.00,  2.00,
                                            -0.01,  1.01,  0.41 };
  const int npts = static_cast<int>(test_points.size()) / 3;
  std::vector<double> analytic_tc(npts), spline_tc(npts);
  std::vector<double> analytic_tc_dx(npts), spline_tc_dx(npts);
  std::vector<double> analytic_tc_dy(npts), spline_tc_dy(npts);
  std::vector<double> analytic_tc_dz(npts), spline_tc_dz(npts);
  std::vector<double> analytic_tc_dxx(npts), spline_tc_dxx(npts);
  std::vector<double> analytic_tc_dxy(npts), spline_tc_dxy(npts);
  std::vector<double> analytic_tc_dxz(npts), spline_tc_dxz(npts);
  std::vector<double> analytic_tc_dyy(npts), spline_tc_dyy(npts);
  std::vector<double> analytic_tc_dyz(npts), spline_tc_dyz(npts);
  std::vector<double> analytic_tc_dzz(npts), spline_tc_dzz(npts);
  std::vector<double> analytic_tc_dxxx(npts), spline_tc_dxxx(npts);
  std::vector<double> analytic_tc_dxxy(npts), spline_tc_dxxy(npts);
  std::vector<double> analytic_tc_dxxz(npts), spline_tc_dxxz(npts);
  std::vector<double> analytic_tc_dxyy(npts), spline_tc_dxyy(npts);
  std::vector<double> analytic_tc_dxyz(npts), spline_tc_dxyz(npts);
  std::vector<double> analytic_tc_dxzz(npts), spline_tc_dxzz(npts);
  std::vector<double> analytic_tc_dyyy(npts), spline_tc_dyyy(npts);
  std::vector<double> analytic_tc_dyyz(npts), spline_tc_dyyz(npts);
  std::vector<double> analytic_tc_dyzz(npts), spline_tc_dyzz(npts);
  std::vector<double> analytic_tc_dzzz(npts), spline_tc_dzzz(npts);
  const double grid_orig_x = grid[0].getCellOrigin(CartesianDimension::X);
  const double grid_orig_y = grid[0].getCellOrigin(CartesianDimension::Y);
  const double grid_orig_z = grid[0].getCellOrigin(CartesianDimension::Z);
  for (int i = 0; i < npts; i++) {
    const double ptx = test_points[(3 * i)    ];
    const double pty = test_points[(3 * i) + 1];
    const double ptz = test_points[(3 * i) + 2];
    analytic_tc[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::VALUE, coeffs);
    const double rel_ptx = (ptx - grid_orig_x);
    const double rel_pty = (pty - grid_orig_y);
    const double rel_ptz = (ptz - grid_orig_z);
    const int gcell_x = (rel_ptx * inv_cell_dimensions[0]) + (rel_pty * inv_cell_dimensions[3]) +
                        (rel_ptz * inv_cell_dimensions[6]);
    const int gcell_y = (rel_pty * inv_cell_dimensions[4]) + (rel_ptz * inv_cell_dimensions[7]);
    const int gcell_z = rel_ptz * inv_cell_dimensions[8];
    const int gcell = (11 * ((11 * gcell_z) + gcell_y)) + gcell_x;
    spline_tc[i] = grid[gcell].evaluate(ptx, pty, ptz);
    analytic_tc_dx[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DX, coeffs);
    analytic_tc_dy[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DY, coeffs);
    analytic_tc_dz[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DZ, coeffs);
    analytic_tc_dxx[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DXX, coeffs);
    analytic_tc_dxy[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DXY, coeffs);
    analytic_tc_dxz[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DXZ, coeffs);
    analytic_tc_dyy[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DYY, coeffs);
    analytic_tc_dyz[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DYZ, coeffs);
    analytic_tc_dzz[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DZZ, coeffs);
    analytic_tc_dxxx[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DXXX, coeffs);
    analytic_tc_dxxy[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DXXY, coeffs);
    analytic_tc_dxxz[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DXXZ, coeffs);
    analytic_tc_dxyy[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DXYY, coeffs);
    analytic_tc_dxyz[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DXYZ, coeffs);
    analytic_tc_dxzz[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DXZZ, coeffs);
    analytic_tc_dyyy[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DYYY, coeffs);
    analytic_tc_dyyz[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DYYZ, coeffs);
    analytic_tc_dyzz[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DYZZ, coeffs);
    analytic_tc_dzzz[i] = triPolynomial(ptx, pty, ptz, FunctionLevel::DZZZ, coeffs);
    const double3 first_derivative = grid[gcell].derivative<double3>(ptx, pty, ptz);
    spline_tc_dx[i] = first_derivative.x; 
    spline_tc_dy[i] = first_derivative.y;
    spline_tc_dz[i] = first_derivative.z;
    const std::vector<double> second_derivative = grid[gcell].secondDerivative(ptx, pty, ptz);
    spline_tc_dxx[i] = second_derivative[0];
    spline_tc_dxy[i] = second_derivative[1];
    spline_tc_dxz[i] = second_derivative[2];
    spline_tc_dyy[i] = second_derivative[4];
    spline_tc_dyz[i] = second_derivative[5];
    spline_tc_dzz[i] = second_derivative[8];
    const std::vector<double> third_derivative = grid[gcell].thirdDerivative(ptx, pty, ptz);
    spline_tc_dxxx[i] = third_derivative[0];
    spline_tc_dxxy[i] = third_derivative[1];
    spline_tc_dxxz[i] = third_derivative[2];
    spline_tc_dxyy[i] = third_derivative[4];
    spline_tc_dxyz[i] = third_derivative[5];
    spline_tc_dxzz[i] = third_derivative[8];
    spline_tc_dyyy[i] = third_derivative[13];
    spline_tc_dyyz[i] = third_derivative[14];
    spline_tc_dyzz[i] = third_derivative[17];
    spline_tc_dzzz[i] = third_derivative[26];
  }
  double etol, etol_d1, etol_d2, etol_d3;
  switch (tcstencil.getKind()) {
  case Interpolant::SMOOTHNESS:
    etol    = 1.0e-8;
    etol_d1 = 1.0e-8;
    etol_d2 = 1.0e-8;
    etol_d3 = 1.0e-8;
    break;
  case Interpolant::FUNCTION_VALUE:
    etol    = 1.0e-8;
    etol_d1 = 1.0e-8;
    etol_d2 = 5.0e-8;
    etol_d3 = 1.0e-7;
    break;
  }
  check(spline_tc, RelationalOperator::EQUAL, Approx(analytic_tc).margin(etol), "Spline-based "
        "computations for a tricubic function do not match the analytic results.");
  check(spline_tc_dx, RelationalOperator::EQUAL, Approx(analytic_tc_dx).margin(etol_d1),
        "Spline-based computations for d/dx partial derivatives of a tricubic function do not "
        "match the analytic results.");
  check(spline_tc_dy, RelationalOperator::EQUAL, Approx(analytic_tc_dy).margin(etol_d1),
        "Spline-based computations for d/dy partial derivatives of a tricubic function do not "
        "match the analytic results.");
  check(spline_tc_dz, RelationalOperator::EQUAL, Approx(analytic_tc_dz).margin(etol_d1),
        "Spline-based computations for d/dz partial derivatives of a tricubic function do not "
        "match the analytic results.");
  check(spline_tc_dxx, RelationalOperator::EQUAL, Approx(analytic_tc_dxx).margin(etol_d2),
        "Spline-based computations for d2/dx2 mixed partial derivatives of a tricubic function "
        "do not match the analytic results.");
  check(spline_tc_dxy, RelationalOperator::EQUAL, Approx(analytic_tc_dxy).margin(etol_d2),
        "Spline-based computations for d2/dxdy mixed partial derivatives of a tricubic function "
        "do not match the analytic results.");
  check(spline_tc_dxz, RelationalOperator::EQUAL, Approx(analytic_tc_dxz).margin(etol_d2),
        "Spline-based computations for d2/dxdz mixed partial derivatives of a tricubic function "
        "do not match the analytic results.");
  check(spline_tc_dyy, RelationalOperator::EQUAL, Approx(analytic_tc_dyy).margin(etol_d2),
        "Spline-based computations for d2/dy2 mixed partial derivatives of a tricubic function "
        "do not match the analytic results.");
  check(spline_tc_dyz, RelationalOperator::EQUAL, Approx(analytic_tc_dyz).margin(etol_d2),
        "Spline-based computations for d2/dydz mixed partial derivatives of a tricubic function "
        "do not match the analytic results.");
  check(spline_tc_dzz, RelationalOperator::EQUAL, Approx(analytic_tc_dzz).margin(etol_d2),
        "Spline-based computations for d2/dz2 mixed partial derivatives of a tricubic function "
        "do not match the analytic results.");
  check(spline_tc_dxxx, RelationalOperator::EQUAL, Approx(analytic_tc_dxxx).margin(etol_d3),
        "Spline-based computations for d3/dx3 derivatives of a tricubic function do not match the "
        "analytic results.");
  check(spline_tc_dxxy, RelationalOperator::EQUAL, Approx(analytic_tc_dxxy).margin(etol_d3),
        "Spline-based computations for d3/dx2dy mixed partial derivatives of a tricubic function "
        "do not match the analytic results.");
  check(spline_tc_dxxz, RelationalOperator::EQUAL, Approx(analytic_tc_dxxz).margin(etol_d3),
        "Spline-based computations for d3/dx2dz mixed partial derivatives of a tricubic function "
        "do not match the analytic results.");
  check(spline_tc_dxyy, RelationalOperator::EQUAL, Approx(analytic_tc_dxyy).margin(etol_d3),
        "Spline-based computations for d3/dxdy2 mixed partial derivatives of a tricubic function "
        "do not match the analytic results.");
  check(spline_tc_dxyz, RelationalOperator::EQUAL, Approx(analytic_tc_dxyz).margin(etol_d3),
        "Spline-based computations for d3/dxdydz mixed partial derivatives of a tricubic function "
        "do not match the analytic results.");
  check(spline_tc_dxzz, RelationalOperator::EQUAL, Approx(analytic_tc_dxzz).margin(etol_d3),
        "Spline-based computations for d3/dxdz2 mixed partial derivatives of a tricubic function "
        "do not match the analytic results.");
  check(spline_tc_dyyy, RelationalOperator::EQUAL, Approx(analytic_tc_dyyy).margin(etol_d3),
        "Spline-based computations for d3/dy3 mixed partial derivatives of a tricubic function "
        "do not match the analytic results.");
  check(spline_tc_dyyz, RelationalOperator::EQUAL, Approx(analytic_tc_dyyz).margin(etol_d3),
        "Spline-based computations for d3/dy2dz mixed partial derivatives of a tricubic function "
        "do not match the analytic results.");
  check(spline_tc_dyzz, RelationalOperator::EQUAL, Approx(analytic_tc_dyzz).margin(etol_d3),
        "Spline-based computations for d3/dydz2 mixed partial derivatives of a tricubic function "
        "do not match the analytic results.");
  check(spline_tc_dzzz, RelationalOperator::EQUAL, Approx(analytic_tc_dzzz).margin(etol_d3),
        "Spline-based computations for d3/dz3 mixed partial derivatives of a tricubic function "
        "do not match the analytic results.");
}

//-------------------------------------------------------------------------------------------------
// Compute derivatives for the electrostatic potential due to a test charge in the space of a
// tricubic mesh element using finite difference approximations.  These derivaatives can be used to
// test the derivatives that the mesh element works from, to ensure that various chain rules are
// being properly computed.
//
// Arguments:
//   q_x:    Cartesian X location of the charge source
//   q_y:    Cartesian Y location of the charge source
//   q_z:    Cartesian Z location of the charge source
//   a:      The position of the corner of inerest along the cell's a axis
//   b:      The position of the corner of inerest along the cell's b axis
//   c:      The position of the corner of inerest along the cell's c axis
//   invu:   Transformation matrix taking fractional coordinates in the mesh axes into real space
//   order:  Order of the derivatives.  Selected derivatives from each tensor will be computed.
//-------------------------------------------------------------------------------------------------
std::vector<double> cellSpaceFD(const double q_x, const double q_y, const double q_z,
                                const double a, const double b, const double c,
                                const std::vector<double> &invu, const int order) {
  
  // Fill out the potential, then each derivative.
  const double inc = pow(2.0, -12.0);
  const double inv_inc = 0.5 / inc;
  std::vector<double> u(27), du_da(27, 0.0), du_db(27, 0.0), du_dc(27, 0.0);
  for (int i = 0; i < 3; i++) {
    const double da = a + (static_cast<double>(i - 1) * inc);
    for (int j = 0; j < 3; j++) {
      const double db = b + (static_cast<double>(j - 1) * inc);
      for (int k = 0; k < 3; k++) {
        const double dc = c + (static_cast<double>(k - 1) * inc);
        const double x = (invu[0] * da) + (invu[3] * db) + (invu[6] * dc);
        const double y =                  (invu[4] * db) + (invu[7] * dc);
        const double z =                                   (invu[8] * dc);
        const double dx = x - q_x;
        const double dy = y - q_y;
        const double dz = z - q_z;
        u[(((k * 3) + j) * 3) + i] = 1.0 / sqrt((dx * dx) + (dy * dy) + (dz * dz));
      }
    }
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        const int ijk_idx = (((k * 3) + j) * 3) + i;
        if (i == 1) {
          du_da[ijk_idx] = (u[ijk_idx + 1] - u[ijk_idx - 1]) * inv_inc;
        }
        if (j == 1) {
          du_db[ijk_idx] = (u[ijk_idx + 3] - u[ijk_idx - 3]) * inv_inc;
        }
        if (k == 1) {
          du_dc[ijk_idx] = (u[ijk_idx + 9] - u[ijk_idx - 9]) * inv_inc;
        }
      }
    }
  }
  std::vector<double> result(3);
  if (order == 1) {
    result[0] = du_da[13];
    result[1] = du_db[13];
    result[2] = du_dc[13];
    return result;
  }
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        int ijk_idx = (((k * 3) + j) * 3) + i;
        
        // Turn du_da into du_dab
        if (i == 1 && j == 1) {
          du_da[ijk_idx] = (du_da[ijk_idx + 3] - du_da[ijk_idx - 3]) * inv_inc;
        }

        // Turn du_db into du_dbc
        if (j == 1 && k == 1) {
          du_db[ijk_idx] = (du_db[ijk_idx + 9] - du_db[ijk_idx - 9]) * inv_inc;          
        }

        // Turn du_dc into du_dac
        if (k == 1 && i == 1) {
          du_dc[ijk_idx] = (du_dc[ijk_idx + 1] - du_dc[ijk_idx - 1]) * inv_inc;          
        }
      }
    }
  }
  if (order == 2) {
    result[0] = du_da[13];
    result[1] = du_dc[13];
    result[2] = du_db[13];
    return result;
  }

  // Turn du_dab into du_dabc
  du_da[13] = (du_da[22] - du_da[4]) * inv_inc;
  result.resize(1);
  result[0] = du_da[13];
  return result;
}

//-------------------------------------------------------------------------------------------------
// Test the construction of a tricubic cell by computing derivatives at each of its corners using
// finite difference calculations in the {a, b, c} coordinate axes, then comparing the results to
// transformations of Cartesian derivatives.  Use a test charge of 1.0 (unitless) placed outside
// the cell to generate the potential function.
//
// Arguments:
//   invu:  Transformation matrix for taking fractional coordinates in the {a, b, c} axes of the
//          tricubic interpolation cell into real space
//-------------------------------------------------------------------------------------------------
void testTricubicCorners(const std::vector<double> &invu) {

  // Create the transformation matrix take Cartesian coordinates into the cell space
  std::vector<double> umat(9);
  invertSquareMatrix(invu, &umat);

  // Compute the function values and derivatives at the cell corners.  The test charge is placed
  // in the (-, -, -) octant to always be outside the cell.
  std::vector<double> u(8), du_dx(8), du_dy(8), du_dz(8), du_dxx(8), du_dxy(8), du_dxz(8);
  std::vector<double> du_dyy(8), du_dyz(8), du_dxxx(8), du_dxxy(8), du_dxxz(8), du_dxyy(8);
  std::vector<double> du_dxyz(8);
  std::vector<double> du_da(8), du_db(8), du_dc(8), du_dab(8), du_dac(8), du_dbc(8), du_dabc(8);
  std::vector<double> fd_da(8), fd_db(8), fd_dc(8), fd_dab(8), fd_dac(8), fd_dbc(8), fd_dabc(8);
  const double q_x = -5.3;
  const double q_y = -4.7;
  const double q_z = -2.4;
  for (int i = 0; i < 2; i++) {
    const double di = i;
    for (int j = 0; j < 2; j++) {
      const double dj = j;
      for (int k = 0; k < 2; k++) {
        const double dk = k;
        const int ijk_idx = (((k * 2) + j) * 2) + i;
        const double x = (invu[0] * di) + (invu[3] * dj) + (invu[6] * dk);
        const double y =                  (invu[4] * dj) + (invu[7] * dk);
        const double z =                                   (invu[8] * dk);
        const double dx = x - q_x;
        const double dy = y - q_y;
        const double dz = z - q_z;
        const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
        const double r = sqrt(r2);
        const double invr = 1.0 / r;
        const double invr2 = invr * invr;

        // Compute the potential and its derivatives in Cartesian space
        u[ijk_idx] = invr;
        const double du  = -invr2;
        const double d2u =  2.0 * invr2 * invr;
        const double d3u = -6.0 * invr2 * invr2;
        du_dx[ijk_idx] = radialFirstDerivative<double>(du, dx, r);
        du_dy[ijk_idx] = radialFirstDerivative<double>(du, dy, r);
        du_dz[ijk_idx] = radialFirstDerivative<double>(du, dz, r);
        du_dxx[ijk_idx] = radialSecondDerivative<double>(du, d2u, dx, r);
        du_dxy[ijk_idx] = radialSecondDerivative<double>(du, d2u, dx, dy, r, r2);
        du_dxz[ijk_idx] = radialSecondDerivative<double>(du, d2u, dx, dz, r, r2);
        du_dyy[ijk_idx] = radialSecondDerivative<double>(du, d2u, dy, r);
        du_dyz[ijk_idx] = radialSecondDerivative<double>(du, d2u, dy, dz, r, r2);
        du_dxxx[ijk_idx] = radialThirdDerivative<double>(du, d2u, d3u, dx, r, r2);
        du_dxxy[ijk_idx] = radialThirdDerivative<double>(du, d2u, d3u, dx, dy, r, r2);
        du_dxxz[ijk_idx] = radialThirdDerivative<double>(du, d2u, d3u, dx, dz, r, r2);
        du_dxyy[ijk_idx] = radialThirdDerivative<double>(du, d2u, d3u, dy, dx, r, r2);
        du_dxyz[ijk_idx] = radialThirdDerivative<double>(du, d2u, d3u, dx, dy, dz, r, r2);

        // Compute the derivatives in the cell space by analytic methods
        du_da[ijk_idx] = (du_dx[ijk_idx] * invu[0]);
        du_db[ijk_idx] = (du_dx[ijk_idx] * invu[3]) + (du_dy[ijk_idx] * invu[4]);
        du_dc[ijk_idx] = (du_dx[ijk_idx] * invu[6]) + (du_dy[ijk_idx] * invu[7]) +
                         (du_dz[ijk_idx] * invu[8]);
        du_dab[ijk_idx] = invu[0] * ((invu[3] * du_dxx[ijk_idx]) + (invu[4] * du_dxy[ijk_idx]));
        const double tcol_dx_sum = (invu[6] * du_dxx[ijk_idx]) + (invu[7] * du_dxy[ijk_idx]) +
                                   (invu[8] * du_dxz[ijk_idx]);
        du_dac[ijk_idx] = invu[0] * tcol_dx_sum;
        du_dbc[ijk_idx] = (invu[3] * tcol_dx_sum) + (invu[4] * ((invu[6] * du_dxy[ijk_idx]) +
                                                                (invu[7] * du_dyy[ijk_idx]) +
                                                                (invu[8] * du_dyz[ijk_idx])));
        du_dabc[ijk_idx] = invu[0] * ((invu[3] * ((invu[6] * du_dxxx[ijk_idx]) +
                                                  (invu[7] * du_dxxy[ijk_idx]) +
                                                  (invu[8] * du_dxxz[ijk_idx]))) +
                                      (invu[4] * ((invu[6] * du_dxxy[ijk_idx]) +
                                                  (invu[7] * du_dxyy[ijk_idx]) +
                                                  (invu[8] * du_dxyz[ijk_idx]))));
        
        // Compute the derivatives in the cell space by finite difference methods
        const std::vector<double> deriv_one   = cellSpaceFD(q_x, q_y, q_z, di, dj, dk, invu, 1);
        const std::vector<double> deriv_two   = cellSpaceFD(q_x, q_y, q_z, di, dj, dk, invu, 2);
        const std::vector<double> deriv_three = cellSpaceFD(q_x, q_y, q_z, di, dj, dk, invu, 3);
        fd_da[ijk_idx] = deriv_one[0];
        fd_db[ijk_idx] = deriv_one[1];
        fd_dc[ijk_idx] = deriv_one[2];
        fd_dab[ijk_idx] = deriv_two[0];
        fd_dac[ijk_idx] = deriv_two[1];
        fd_dbc[ijk_idx] = deriv_two[2];
        fd_dabc[ijk_idx] = deriv_three[0];
      }
    }
  }
  
  // Run tests
  check(fd_da, RelationalOperator::EQUAL, du_da, "Calculations of df / da made analytically do "
        "not match those computed by finite difference approximations.");
  check(fd_db, RelationalOperator::EQUAL, du_db, "Calculations of df / db made analytically do "
        "not match those computed by finite difference approximations.");
  check(fd_dc, RelationalOperator::EQUAL, du_dc, "Calculations of df / dc made analytically do "
        "not match those computed by finite difference approximations.");
  check(fd_dab, RelationalOperator::EQUAL, du_dab, "Calculations of d2f / dab made analytically "
        "do not match those computed by finite difference approximations.");
  check(fd_dac, RelationalOperator::EQUAL, du_dac, "Calculations of d2f / dac made analytically "
        "do not match those computed by finite difference approximations.");
  check(fd_dbc, RelationalOperator::EQUAL, du_dbc, "Calculations of d2f / dbc made analytically "
        "do not match those computed by finite difference approximations.");
  check(fd_dabc, RelationalOperator::EQUAL, Approx(du_dabc).margin(1.0e-5), "Calculations of "
        "d3f / dabc made analytically do not match those computed by finite difference "
        "approximations.");
}

//-------------------------------------------------------------------------------------------------
// Test the efficacy of the FUNCTION_VALUE stencil mode for tricubic interpolation of typical
// non-bonded potentials on a triclinic mesh element.
//-------------------------------------------------------------------------------------------------
void testCoulombInterpolation(const std::vector<double> &invu, const TricubicStencil &tcs,
                              Xoshiro256ppGenerator *xrs) {

  // Design a set of charges to influence the mesh element
  const int nq = 8;
  std::vector<double> q_x(nq), q_y(nq), q_z(nq), q_v(nq);
  for (int i = 0; i < nq; i++) {
    q_x[i] = -6.0 + (5.0 * xrs->uniformRandomNumber());
    q_y[i] = -6.0 + (5.0 * xrs->uniformRandomNumber());
    q_z[i] = -6.0 + (5.0 * xrs->uniformRandomNumber());
    q_v[i] = -1.0 + (2.0 * xrs->uniformRandomNumber());
  }

  // Compute the locations of midpoints for extra potential calculations in the FUNCTION_VALUE
  // interpolant case.
  std::vector<double> xyblock_x(8), xzblock_x(8), yzblock_x(8), cnblock_x(8);
  std::vector<double> xyblock_y(8), xzblock_y(8), yzblock_y(8), cnblock_y(8);
  std::vector<double> xyblock_z(8), xzblock_z(8), yzblock_z(8), cnblock_z(8);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        const int idx = (((k * 2) + j) * 2) + i;
        double mpt_x, mpt_y, mpt_z;
        fvStencilCoordinates(0.0, 0.0, 0.0, UnitCellAxis::C, i, j, k, invu.data(),
                                   &mpt_x, &mpt_y, &mpt_z);
        xyblock_x[idx] = mpt_x;
        xyblock_y[idx] = mpt_y;
        xyblock_z[idx] = mpt_z;
        fvStencilCoordinates(0.0, 0.0, 0.0, UnitCellAxis::B, i, j, k, invu.data(),
                                   &mpt_x, &mpt_y, &mpt_z);
        xzblock_x[idx] = mpt_x;
        xzblock_y[idx] = mpt_y;
        xzblock_z[idx] = mpt_z;
        fvStencilCoordinates(0.0, 0.0, 0.0, UnitCellAxis::A, i, j, k, invu.data(),
                                   &mpt_x, &mpt_y, &mpt_z);
        yzblock_x[idx] = mpt_x;
        yzblock_y[idx] = mpt_y;
        yzblock_z[idx] = mpt_z;
        fvStencilCoordinates(0.0, 0.0, 0.0, i, j, k, invu.data(), &mpt_x, &mpt_y, &mpt_z);
        cnblock_x[idx] = mpt_x;
        cnblock_y[idx] = mpt_y;
        cnblock_z[idx] = mpt_z;
      }
    }
  }
  
  // Compute the potential and derivatives for the stencil on an element with origin
  // { 0.0, 0.0, 0.0 }.
  std::vector<double> u(8, 0.0), dudx(8, 0.0), dudy(8, 0.0), dudz(8, 0.0), dudxy(8, 0.0);
  std::vector<double> dudxz(8, 0.0), dudyz(8, 0.0), dudxyz(8, 0.0), dudxx(8, 0.0), dudyy(8, 0.0);
  std::vector<double> dudxxx(8, 0.0), dudxxy(8, 0.0), dudxxz(8, 0.0), dudxyy(8, 0.0);
  for (int i = 0; i < 2; i++) {
    const double d_i = i;
    for (int j = 0; j < 2; j++) {
      const double d_j = j;
      for (int k = 0; k < 2; k++) {
        const double d_k = k;
        const int idx = (((k * 2) + j) * 2) + i;
        const double vert_x = (invu[0] * d_i) + (invu[3] * d_j) + (invu[6] * d_k);
        const double vert_y =                   (invu[4] * d_j) + (invu[7] * d_k);
        const double vert_z =                                     (invu[8] * d_k);
        for (int m = 0; m < nq; m++) {
          const double dx = vert_x - q_x[m];
          const double dy = vert_y - q_y[m];
          const double dz = vert_z - q_z[m];
          const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
          const double r = sqrt(r2);
          const double invr = 1.0 / r;
          const double invr2 = invr * invr;
          const double dqqe  = -invr2;
          const double d2qqe =  2.0 * invr2 * invr;
          const double d3qqe = -6.0 * invr2 * invr2;
          u[idx]    += q_v[m] * invr;
          dudx[idx] += q_v[m] * radialFirstDerivative(dqqe, dx, r);
          dudy[idx] += q_v[m] * radialFirstDerivative(dqqe, dy, r);
          dudz[idx] += q_v[m] * radialFirstDerivative(dqqe, dz, r);
          switch (tcs.getKind()) {
          case Interpolant::SMOOTHNESS:
            dudxx[idx] += q_v[m] * radialSecondDerivative(dqqe, d2qqe, dx, r);
            dudxy[idx] += q_v[m] * radialSecondDerivative(dqqe, d2qqe, dx, dy, r, r2);
            dudxz[idx] += q_v[m] * radialSecondDerivative(dqqe, d2qqe, dx, dz, r, r2);
            dudyy[idx] += q_v[m] * radialSecondDerivative(dqqe, d2qqe, dy, r);
            dudyz[idx] += q_v[m] * radialSecondDerivative(dqqe, d2qqe, dy, dz, r, r2);
            dudxxx[idx] += q_v[m] * radialThirdDerivative(dqqe, d2qqe, d3qqe, dx, r, r2);
            dudxxy[idx] += q_v[m] * radialThirdDerivative(dqqe, d2qqe, d3qqe, dx, dy, r, r2);
            dudxxz[idx] += q_v[m] * radialThirdDerivative(dqqe, d2qqe, d3qqe, dx, dz, r, r2);
            dudxyy[idx] += q_v[m] * radialThirdDerivative(dqqe, d2qqe, d3qqe, dy, dx, r, r2);
            dudxyz[idx] += q_v[m] * radialThirdDerivative(dqqe, d2qqe, d3qqe, dx, dy, dz, r, r2);
            break;
          case Interpolant::FUNCTION_VALUE:
            {
              double mdx = xyblock_x[idx] - q_x[m];
              double mdy = xyblock_y[idx] - q_y[m];
              double mdz = xyblock_z[idx] - q_z[m];
              dudxy[idx] += q_v[m] / sqrt((mdx * mdx) + (mdy * mdy) + (mdz * mdz));
              mdx = xzblock_x[idx] - q_x[m];
              mdy = xzblock_y[idx] - q_y[m];
              mdz = xzblock_z[idx] - q_z[m];
              dudxz[idx] += q_v[m] / sqrt((mdx * mdx) + (mdy * mdy) + (mdz * mdz));
              mdx = yzblock_x[idx] - q_x[m];
              mdy = yzblock_y[idx] - q_y[m];
              mdz = yzblock_z[idx] - q_z[m];
              dudyz[idx] += q_v[m] / sqrt((mdx * mdx) + (mdy * mdy) + (mdz * mdz));
              mdx = cnblock_x[idx] - q_x[m];
              mdy = cnblock_y[idx] - q_y[m];
              mdz = cnblock_z[idx] - q_z[m];
              dudxyz[idx] += q_v[m] / sqrt((mdx * mdx) + (mdy * mdy) + (mdz * mdz));
            }
            break;
          }
        }
      }
    }
  }

  // Some arrays may not have been initialized, but these will not be used.
  std::vector<double> tc_bounds(12, 0.0);
  for (int i = 0; i < 9; i++) {
    tc_bounds[i + 3] = invu[i];
  }
  TricubicCell<double> tc(tcs, tc_bounds, u, dudx, dudy, dudz, dudxx, dudxy, dudxz, dudyy, dudyz,
                          dudxxx, dudxxy, dudxxz, dudxyy, dudxyz);

  // Determine a spread of points in the element
  const int npts = 1000;
  std::vector<double> u_ans(npts, 0.0);
  std::vector<double> dudx_ans(npts, 0.0), dudy_ans(npts, 0.0), dudz_ans(npts, 0.0);
  std::vector<double> u_itp(npts), dudx_itp(npts), dudy_itp(npts), dudz_itp(npts);
  for (int i = 0; i < npts; i++) {
    const double loc_a = xrs->uniformRandomNumber();
    const double loc_b = xrs->uniformRandomNumber();
    const double loc_c = xrs->uniformRandomNumber();
    const double loc_x = (invu[0] * loc_a) + (invu[3] * loc_b) + (invu[6] * loc_c);
    const double loc_y =                     (invu[4] * loc_b) + (invu[7] * loc_c);
    const double loc_z =                                         (invu[8] * loc_c);
    for (int j = 0; j < nq; j++) {
      const double dx = loc_x - q_x[j];
      const double dy = loc_y - q_y[j];
      const double dz = loc_z - q_z[j];
      const double r = sqrt((dx * dx) + (dy * dy) + (dz * dz));
      const double invr = 1.0 / r;
      const double invr2 = invr * invr;
      u_ans[i] += q_v[j] * invr;
      dudx_ans[i] += q_v[j] * radialFirstDerivative(-invr2, dx, r);
      dudy_ans[i] += q_v[j] * radialFirstDerivative(-invr2, dy, r);
      dudz_ans[i] += q_v[j] * radialFirstDerivative(-invr2, dz, r);
    }
    u_itp[i] = tc.evaluate(loc_x, loc_y, loc_z);
    const double3 itp_dxyz = tc.derivative<double3>(loc_x, loc_y, loc_z);
    dudx_itp[i] = itp_dxyz.x;
    dudy_itp[i] = itp_dxyz.y;
    dudz_itp[i] = itp_dxyz.z;
  }

  // Check that the errors are reasonable
  std::string space_group;
  double err_margin;
  if (fabs(invu[3]) > tiny || fabs(invu[6]) > tiny || fabs(invu[7]) > tiny) {
    space_group = getEnumerationName(UnitCellType::TRICLINIC);
    err_margin = 1.5e-4;
  }
  else {
    space_group = getEnumerationName(UnitCellType::ORTHORHOMBIC);
    err_margin = 1.0e-4;
  }
  check(u_ans, RelationalOperator::EQUAL, Approx(u_itp).margin(err_margin), "Energies of "
        "Coulombic interactions interpolated from a " + getEnumerationName(tcs.getKind()) +
        " interpolant do not match the analytic results.  The element geometry is " + space_group +
        ".");
  check(dudx_ans, RelationalOperator::EQUAL, Approx(dudx_itp).margin(err_margin), "X-axis forces "
        "due to Coulombic interactions interpolated from a " + getEnumerationName(tcs.getKind()) +
        " interpolant do not match the analytic results.  The element geometry is " + space_group +
        ".");
  check(dudy_ans, RelationalOperator::EQUAL, Approx(dudy_itp).margin(err_margin), "Y-axis forces "
        "due to Coulombic interactions interpolated from a " + getEnumerationName(tcs.getKind()) +
        " interpolant do not match the analytic results.  The element geometry is " + space_group +
        ".");
  check(dudz_ans, RelationalOperator::EQUAL, Approx(dudz_itp).margin(err_margin), "Z-axis forces "
        "due to Coulombic interactions interpolated from a " + getEnumerationName(tcs.getKind()) +
        " interpolant do not match the analytic results.  The element geometry is " + space_group +
        ".");
}

//-------------------------------------------------------------------------------------------------
// Test the sigmoid and sigmoidf functions to verify their returned values and derivatives.
//-------------------------------------------------------------------------------------------------
void testSigmoid(const double crossover, const double intensity) {
  const int npts = 40;
  std::vector<double4> sigm(npts);
  std::vector<double> sigm_x(npts), sigm_y(npts), sigm_z(npts), sigm_w(npts), vchk(npts);
  std::vector<double> sole_x(npts), sole_y(npts), sole_z(npts), sole_w(npts);
  std::vector<double> fd_one(npts), fd_two(npts), fd_thr(npts);
  std::vector<float> fsigm_x(npts), fsigm_y(npts), fsigm_z(npts), fsigm_w(npts);
  std::vector<float> fsole_x(npts), fsole_y(npts), fsole_z(npts), fsole_w(npts);
  const double delta = pow(2.0, -16.0);
  double r = -10.0;
  for (int i = 0; i < npts; i++) {
    sigm[i] = sigmoid(r, crossover, intensity);
    const double4 sigm_p = sigmoid(r + delta, crossover, intensity);
    const double4 sigm_n = sigmoid(r - delta, crossover, intensity);
    vchk[i] = 1.0 / (exp(intensity * (r - crossover)) + 1.0);
    sigm_x[i] = sigm[i].x;
    sigm_y[i] = sigm[i].y;
    sigm_z[i] = sigm[i].z;
    sigm_w[i] = sigm[i].w;
    fd_one[i] = (sigm_p.x - sigm_n.x) / (2.0 * delta);
    fd_two[i] = (sigm_p.y - sigm_n.y) / (2.0 * delta);
    fd_thr[i] = (sigm_p.z - sigm_n.z) / (2.0 * delta);
    sole_x[i] = sigmoid(r, crossover, intensity, 0);
    sole_y[i] = sigmoid(r, crossover, intensity, 1);
    sole_z[i] = sigmoid(r, crossover, intensity, 2);
    sole_w[i] = sigmoid(r, crossover, intensity, 3);
    const float4 fsigm = sigmoidf(r, crossover, intensity);
    fsigm_x[i] = fsigm.x;
    fsigm_y[i] = fsigm.y;
    fsigm_z[i] = fsigm.z;
    fsigm_w[i] = fsigm.w;
    fsole_x[i] = sigmoidf(r, crossover, intensity, 0);
    fsole_y[i] = sigmoidf(r, crossover, intensity, 1);
    fsole_z[i] = sigmoidf(r, crossover, intensity, 2);
    fsole_w[i] = sigmoidf(r, crossover, intensity, 3);
    r += 0.5;
  }
  check(sigm_x, RelationalOperator::EQUAL, Approx(vchk).margin(stormm::constants::tiny),
        "The sigmoid function does not return the expected function values.");
  check(sigm_y, RelationalOperator::EQUAL, fd_one, "The sigmoid function does not agree with the "
        "first derivative computed by a finite difference approximation.");
  check(sigm_z, RelationalOperator::EQUAL, fd_two, "The sigmoid function does not agree with the "
        "second derivative computed by a finite difference approximation.");
  check(sigm_w, RelationalOperator::EQUAL, fd_thr, "The sigmoid function does not agree with "
        "the third derivative computed by a finite difference approximation.");
  check(sole_x, RelationalOperator::EQUAL, sigm_x, "The sigmoid function does not return the "
        "expected result when called explicitly for the function value.");
  check(sole_y, RelationalOperator::EQUAL, sigm_y, "The sigmoid function does not return the "
        "expected result when called explicitly for the first derivative.");
  check(sole_z, RelationalOperator::EQUAL, sigm_z, "The sigmoid function does not return the "
        "expected result when called explicitly for the second derivative.");
  check(sole_w, RelationalOperator::EQUAL, sigm_w, "The sigmoid function does not return the "
        "expected result when called explicitly for the third derivative.");
  check(fsigm_x, RelationalOperator::EQUAL, sigm_x, "The sigmoidf function does not return the "
        "correct result in the function value.");
  check(fsigm_y, RelationalOperator::EQUAL, sigm_y, "The sigmoidf function does not return the "
        "correct result in the first derivative.");
  check(fsigm_z, RelationalOperator::EQUAL, sigm_z, "The sigmoidf function does not return the "
        "correct result in the second derivative.");
  check(fsigm_w, RelationalOperator::EQUAL, Approx(sigm_w).margin(1.65e-6), "The sigmoidf "
        "function does not return the correct result in the third derivative.");
  check(fsole_x, RelationalOperator::EQUAL, sigm_x, "The sigmoidf function does not return the "
        "correct result when queried for the function value.");
  check(fsole_y, RelationalOperator::EQUAL, sigm_y, "The sigmoidf function does not return the "
        "correct result when queried for the first derivative.");
  check(fsole_z, RelationalOperator::EQUAL, sigm_z, "The sigmoidf function does not return the "
        "correct result when queried for the second derivative.");
  check(fsole_w, RelationalOperator::EQUAL, Approx(sigm_w).margin(1.65e-6), "The sigmoidf "
        "function does not return the correct result when queried for the third derivative.");
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
  
  // Section 1
  section("Vector processing capabilities");

  // Section 2
  section("Random number generation");

  // Section 3
  section("Rounding and prime factorization");

  // Section 4
  section("Lightweight matrix math from the hand-coded routines");

  // Section 5
  section("Test the TickCounter object");

  // Section 6
  section("Test tricubic interpolation methods");

  // Section 7
  section("Customized formulas");
  
  // Check vector processing capabilities
  section(1);
  const std::vector<double> dv_i = { 0.1, 0.5, 0.9, 0.7, 0.8 };
  check(mean(dv_i), RelationalOperator::EQUAL, Approx(0.6).margin(1.0e-12), "The mean value of a "
        "simple real number vector is incorrect.");
  const std::vector<int> iv_i = { 1, 5, 9, 7, 38 };
  check(mean(iv_i), RelationalOperator::EQUAL, 12, "The mean value of a "
        "simple integer vector is incorrect.");
  const std::vector<double> dv_i_long = { -0.3, -0.6, 0.1, 0.5, 0.9, 0.7, 0.8 };
  CHECK_THROWS(maxAbsoluteDifference(dv_i_long, dv_i), "Vectors of different lengths were "
               "compared to obtain a maximum difference.");
  const std::vector<double> dv_ii = { 2.1, 3.5, -9.9, 4.7, 7.8 };
  check(maxAbsoluteDifference(dv_i, dv_ii), RelationalOperator::EQUAL, Approx(10.8).margin(1.0e-8),
        "Incorrect evaluation of the maximum absolute difference between two short vectors.");
  const std::vector<int> iv_ii = { -1, 15, -9, 72, -3 };
  check(maxAbsoluteDifference(iv_i, iv_ii), RelationalOperator::EQUAL, Approx(65).margin(1.0e-10),
        "Incorrect evaluation of the maximum absolute difference between two short vectors.");
  const std::vector<int> uv_i = { 92, 47, 81, 36, 29 };
  const std::vector<int> uv_ii = { 82, 37, 101, 19, 147 };
  check(maxAbsoluteDifference(uv_i, uv_ii), RelationalOperator::EQUAL,
        maxAbsoluteDifference(uv_ii, uv_i), "Non-reflexive evaluation of the maximum absolute "
        "difference between two vectors of unsigned integers.");
  check(maxRelativeDifference(dv_ii, dv_i), RelationalOperator::EQUAL, Approx(20.0).margin(1.0e-8),
        "Incorrect evaluation of maximum relative error between two short vectors.");
  check(maxRelativeDifference(uv_ii, uv_i), RelationalOperator::EQUAL,
        Approx(4.06896552).margin(1.0e-7), "Incorrect evaluation of maximum relative error "
        "between two short unsigned integer vectors.");
  check(mean(uv_i), RelationalOperator::EQUAL, 57.0, "Incorrect evaluation of the mean of a "
        "vector of unsigned integers.");
  const std::vector<double> cp_a = {  1.4,  5.8, -8.5 };
  const std::vector<double> cp_b = { -3.4,  2.7, 10.1 };
  std::vector<double> cp_c(3);
  crossProduct(cp_a, cp_b, &cp_c);
  const std::vector<double> cp_ab = { 81.53, 14.76, 23.50 };
  check(cp_c, RelationalOperator::EQUAL, cp_ab, "Vector cross product does not function as "
        "expected when provided Standard Template Library vectors.");
  const std::vector<float> cp_af = {  1.4, 5.8, -8.5 };
  const std::vector<float> cp_bf = { -3.4, 2.7, 10.1 };
  std::vector<float> cp_cf(3);
  crossProduct(cp_af.data(), cp_bf.data(), cp_cf.data());
  check(cp_cf, RelationalOperator::EQUAL, Approx(cp_ab).margin(1.0e-5), "Vector cross product "
        "does not function as expected when provided C-style pointers.");
  const double3 d3_a = { cp_a[0],  cp_a[1],  cp_a[2] };
  const double3 d3_b = { cp_b[0],  cp_b[1],  cp_b[2] };
  const double3 d3_c = crossProduct(d3_a, d3_b);
  check(d3_c.x == Approx(cp_ab[0]).margin(tiny) && d3_c.y == Approx(cp_ab[1]).margin(tiny) &&
        d3_c.z == Approx(cp_ab[2]).margin(tiny), "Vector cross product does not function as "
        "expected when provided double-precision HPC tuples.");
  const std::vector<int> i_searchable_asc = {  0,  1,  2,  3,  4,  5, 16, 17, 28, 99 };
  const std::vector<int> i_searchable_des = {  9,  8,  6,  5,  4,  2,  1, -3, -4, -6 };
  const std::vector<int> i_searchable_non = {  9, 18, 46,  5, -4,  2,  1, -3, 15, -8 };
  std::vector<int> search_order = { 0, 4, 8, 2, 6, 1, 5, 3, 7, 9 };
  const int nsearch = i_searchable_asc.size();
  std::vector<int> asc_result(nsearch), des_result(nsearch), non_result(nsearch);
  for (int i = 0; i < nsearch; i++) {
    asc_result[i] = locateValue(i_searchable_asc, i_searchable_asc[search_order[i]],
                                DataOrder::ASCENDING);
    des_result[i] = locateValue(i_searchable_des, i_searchable_des[search_order[i]],
                                DataOrder::DESCENDING);
    non_result[i] = locateValue(i_searchable_non, i_searchable_non[search_order[i]],
                                DataOrder::NONE);
  }
  check(asc_result, RelationalOperator::EQUAL, search_order, "Searching a data set arranged in "
        "ascending order yielded incorrect results.");
  check(des_result, RelationalOperator::EQUAL, search_order, "Searching a data set arranged in "
        "descending order yielded incorrect results.");
  check(non_result, RelationalOperator::EQUAL, search_order, "Searching a data set arranged in "
        "no discernible order yielded incorrect results.");
  for (int i = 0; i < nsearch - 1; i++) {
    asc_result[i] = locateValue(i_searchable_asc.data(), i_searchable_asc[search_order[i]],
                                nsearch - 1, DataOrder::ASCENDING);
    des_result[i] = locateValue(i_searchable_des.data(), i_searchable_des[search_order[i]],
                                nsearch - 1, DataOrder::DESCENDING);
    non_result[i] = locateValue(i_searchable_non.data(), i_searchable_non[search_order[i]],
                                nsearch - 1, DataOrder::NONE);
  }
  asc_result.resize(nsearch - 1);
  des_result.resize(nsearch - 1);
  non_result.resize(nsearch - 1);
  search_order.resize(nsearch - 1);
  check(asc_result, RelationalOperator::EQUAL, search_order, "Searching a data set of size " +
        std::to_string(nsearch - 1) + ", arranged in ascending order, yielded incorrect results.");
  check(des_result, RelationalOperator::EQUAL, search_order, "Searching a data set of size " +
        std::to_string(nsearch - 1) +
        ", arranged in descending order, yielded incorrect results.");
  check(non_result, RelationalOperator::EQUAL, search_order, "Searching a data set of size " +
        std::to_string(nsearch - 1) +
        ", arranged in no discernible order, yielded incorrect results.");
  const std::vector<double> d_searchable_asc = { 0.1, 0.7, 2.3, 6.8, 7.1, 7.2, 9.3, 9.5, 9.7 };
  const std::vector<double> d_searchable_des = { 6.1, 5.7, 4.3, 3.8, 3.1, 2.2, 1.3, 0.5, 0.2 };
  const std::vector<double> d_searchable_non = { 5.1, 1.7, 4.3, 9.8, 3.1, 0.2, 0.3, 7.5, 1.8 };
  check(locateValue(Approx(2.35).margin(0.06), d_searchable_asc, DataOrder::ASCENDING),
        RelationalOperator::EQUAL, 2, "Binary search in ascending order with an approximate "
        "target value that should hit its mark from above fails.");
  check(locateValue(Approx(7.05).margin(0.06), d_searchable_asc, DataOrder::ASCENDING),
        RelationalOperator::EQUAL, 4, "Binary search in ascending order with an approximate "
        "target value that should hit its mark from below fails.");
  check(locateValue(Approx(0.55).margin(0.06), d_searchable_des, DataOrder::DESCENDING),
        RelationalOperator::EQUAL, 7, "Binary search in descending order with an approximate "
        "target value that should hit its mark from above fails.");
  check(locateValue(Approx(3.75).margin(0.06), d_searchable_des, DataOrder::DESCENDING),
        RelationalOperator::EQUAL, 3, "Binary search in descending order with an approximate "
        "target value that should hit its mark from below fails.");
  const std::vector<double> productize_this = { 1.5, 1.8, 2.9, 3.1, -2.0 };
  check(seriesProduct<double>(productize_this), RelationalOperator::EQUAL,
        Approx(1.5 * 1.8 * 2.9 * 3.1 * -2.0).margin(tiny), "The product of a short series of "
        "real numbers was not computed properly.");
  const std::vector<double> poly_product = tileVector(productize_this, 6);
  check(seriesProduct<double>(poly_product), RelationalOperator::EQUAL,
        Approx(13089429228.941328).margin(1.0e-4), "The product of a long series of real "
        "numbers was not computed correctly.");
  const std::vector<double> poly_product_ii = tileVector(productize_this, 400.0);
  CHECK_THROWS(seriesProduct<double>(poly_product_ii), "The product of a very long series of "
               "numbers was computed, despite the fact that this would break the double-precision "
               "format.");
  const std::vector<int> iproductize_this = { 2, 1, 3, -3, 4, -5, 1, 7 };
  check(seriesProduct<int>(iproductize_this), RelationalOperator::EQUAL, 2520, "The product of a "
        "short series of integers was not computed correctly.");
  const std::vector<int> ipoly_product = tileVector(iproductize_this, 3);
  check(seriesProduct<llint>(ipoly_product), RelationalOperator::EQUAL, 16003008000LL,
        "The product of a triple-length series of integers was not computed correctly.");
  CHECK_THROWS(seriesProduct<int>(ipoly_product), "The product of a triple-length series of "
               "integers was computed, despite the fact that it would overflow the integer "
               "format.");

  // Check series formation capabilities
  const std::vector<int> inc_eight = incrementingSeries(0, 8, 1);
  const std::vector<int> inc_eight_ans = { 0, 1, 2, 3, 4, 5, 6, 7 };
  check(inc_eight, RelationalOperator::EQUAL, inc_eight_ans, "An incrementing series of eight "
        "values was not formed correctly.");
  const std::vector<int> inc_fourteen = incrementingSeries(0, 14, 2);
  const std::vector<int> inc_fourteen_ans = { 0, 2, 4, 6, 8, 10, 12 };
  check(inc_fourteen, RelationalOperator::EQUAL, inc_fourteen_ans, "An incrementing series of "
        "seven values and intervals of two was not formed correctly.");
  const std::vector<int> dec_eight = incrementingSeries(8, 0, 1);
  const std::vector<int> dec_eight_ans = { 8, 7, 6, 5, 4, 3, 2, 1 };
  check(dec_eight, RelationalOperator::EQUAL, dec_eight_ans, "A decrementing series of eight "
        "values was not formed correctly.");
  const std::vector<int> mixed_ints = { 0, 1, 2, 3, 4, 4, 2, 2, 3, 1, 0, 0, 1, 2, 3, 1, 0 };
  const std::vector<int> locations_ans = { 0, 10, 11, 16, 1, 9, 12, 15, 2, 6, 7, 13, 3, 8, 14, 4,
                                           5 };
  const std::vector<int> location_bounds_ans = { 0, 4, 8, 12, 15, 17 };
  std::vector<int> locations(17);
  std::vector<int> location_bounds(6);
  indexingArray(mixed_ints, &locations, &location_bounds, 5);
  check(locations, RelationalOperator::EQUAL, locations_ans, "A series of mixed integers was not "
        "processed correctly to reveal the locations of various values.");
  check(location_bounds, RelationalOperator::EQUAL, location_bounds_ans, "The bounds array for "
        "locations of a series of mixed integers was not formed correctly.");

  // Verify that the internal random number generation is consistent with expectations.  Create
  // three generators, the first two initialized in the same way (which thus should track one
  // another precisely) and the third initialized to a different value (which should decorrelate
  // from the first two immediately thanks to some priming that is done in the initialization).
  section(2);
  check(oe.getRandomSeed(), RelationalOperator::EQUAL, 827493, "Pseudo-random number seed is not "
        "set to the expected value.  While permitted as part of the baseline test environment "
        "initialization, changing the random seed will, in this case, result in several bad "
        "results in later tests.");
  const int n_pts = 2000;
  Ran2Generator prng_a(oe.getRandomSeed());
  Ran2Generator prng_b(oe.getRandomSeed());
  Ran2Generator prng_c(71277);
  std::vector<double> result_a(n_pts, 0.0);
  std::vector<double> result_b(n_pts, 0.0);
  std::vector<double> result_c(n_pts, 0.0);
  for (int i = 0; i < n_pts; i++) {
    result_a[i] = prng_a.uniformRandomNumber();
    result_b[i] = prng_b.uniformRandomNumber();
    result_c[i] = prng_c.gaussianRandomNumber();
  }
  check(result_a, RelationalOperator::EQUAL,
        Approx(result_b, ComparisonType::MEAN_UNSIGNED_ERROR).margin(1.0e-12), "Differences were "
        "detected between random number series created starting from the same seed.");
  check(mean(result_a), RelationalOperator::EQUAL, Approx(0.49674889).margin(1.0e-7), "The mean "
	"of a set of " + std::to_string(n_pts) + " random numbers is incorrect.");
  check(mean(result_c), RelationalOperator::EQUAL, Approx(-0.022698028552).margin(1.0e-7),
        "The mean value of a normal distribution of random numbers is incorrect.");

  // Additional checks, using the file reference system
  const std::string randoms_snp = oe.getStormmSourcePath() + osSeparator() + "test" +
                                  osSeparator() + "Math" + osSeparator() + "randoms.m";
  TestPriority snp_found = (getDrivePathType(randoms_snp) == DrivePathType::FILE) ?
                           TestPriority::CRITICAL : TestPriority::ABORT;
  if (snp_found == TestPriority::ABORT && oe.takeSnapshot() != SnapshotOperation::SNAPSHOT) {
    rtWarn("Snapshot file " + randoms_snp + " was not found.  Check the $STORMM_SOURCE "
           "environment variable and make sure it indicates the root source directory where src/ "
           "and test/ can be found.  Some subsequent tests will be skipped.", "test_math");
  }
  snapshot(oe.getStormmSourcePath() + osSeparator() + "test" + osSeparator() + "Math" +
           osSeparator() + "randoms.m", polyNumericVector(result_c), "rngvec", 1.0e-4,
           "Series of random numbers created by the ran2 method does not conform to expectations.",
           oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::OVERWRITE,
           snp_found);

  // Check scrambled linear random number generators
  Xoroshiro128pGenerator xrs128p(798031);
  Xoshiro256ppGenerator xrs256pp(901835);
  Xoroshiro128pGenerator xrs128pj(798031);
  Xoshiro256ppGenerator xrs256ppj(901835);
  xrs128pj.jump();
  xrs256ppj.jump();
  const int nrand_trial = 16;
  std::vector<double> xrs128p_result(nrand_trial);
  std::vector<double> xrs128p_jump_result(nrand_trial);
  std::vector<double> xrs256pp_result(nrand_trial);
  std::vector<double> xrs256pp_jump_result(nrand_trial);
  for (int i = 0; i < nrand_trial; i++) {
    xrs128p_result[i] = xrs128p.uniformRandomNumber();
    xrs256pp_result[i] = xrs256pp.uniformRandomNumber();
    xrs128pj.uniformRandomNumber();
    xrs256ppj.uniformRandomNumber();
  }
  xrs128p.jump();
  xrs256pp.jump();
  const std::vector<double> ans_128p = { 0.7453648164, 0.4049923254, 0.8584963726, 0.9833355535,
                                         0.0066062865, 0.6311114017, 0.9820136114, 0.2413733841,
                                         0.2753459418, 0.4993040685, 0.0806123499, 0.3691566725,
                                         0.4001401073, 0.0590209187, 0.5804605659, 0.5293153466 };
  const std::vector<double> ans_256pp = { 0.9570185700, 0.9443435021, 0.5241518529, 0.1428242166,
                                          0.5687186981, 0.9839369524, 0.9751010737, 0.8695307120,
                                          0.4293373290, 0.1555431764, 0.0005355056, 0.0820197480,
                                          0.8509827801, 0.7310430980, 0.7770864125, 0.9021266541 };
  check(xrs128p_result, RelationalOperator::EQUAL, Approx(ans_128p).margin(1.0e-8),
        "Random numbers generated by the xoroshift128+ method do not meet expectations.");
  check(xrs256pp_result, RelationalOperator::EQUAL, Approx(ans_256pp).margin(1.0e-8),
        "Random numbers generated by the xoshift256++ method do not meet expectations.");
  for (int i = 0; i < nrand_trial; i++) {
    xrs128p_result[i] = xrs128p.uniformRandomNumber();
    xrs128p_jump_result[i] = xrs128pj.uniformRandomNumber();
    xrs256pp_result[i] = xrs256pp.uniformRandomNumber();
    xrs256pp_jump_result[i] = xrs256ppj.uniformRandomNumber();
  }
  check(xrs128p_result, RelationalOperator::EQUAL, xrs128p_jump_result, "Two xoroshiro128+ "
        "generators do not re-synchronize after different combinations of random bit string "
        "generation and a jump.");
  check(xrs256pp_result, RelationalOperator::EQUAL, xrs256pp_jump_result, "Two xoroshiro256++ "
        "generators do not re-synchronize after different combinations of random bit string "
        "generation and a jump.");
  Xoshiro256ppGenerator xrs256ppb({ 0x7743a154e17a5e9bLLU, 0x7823a1cd9453899bLLU,
                                    0x976589eefbb1c7f5LLU, 0x702cf168260fa29eLLU });
  const std::vector<ullint> orbit = { 0xd5c766557e6e16e4LLU, 0xf8eb8f8747a8cc67LLU,
                                      0xfc18365710a653eeLLU, 0xc698a193593f232LLU,
                                      0xa44ddeaac93b1be7LLU, 0x678dcd0e0516c741LLU,
                                      0x351d94668d7c35eLLU,  0x73160e24fc8768daLLU,
                                      0x562ca11220c31698LLU, 0x3dba336235c48913LLU };
  std::vector<ullint> xrsb_result(orbit.size());
  for (size_t i = 0; i < orbit.size(); i++) {
    xrsb_result[i] = xrs256ppb.revealBitString();
    xrs256ppb.uniformRandomNumber();
  }
  std::vector<int> xrsb_diff(orbit.size());
  for (size_t i = 0; i < orbit.size(); i++) {
    xrsb_diff[i] = xrsb_result[i] - orbit[i];
  }
  check(xrsb_diff, RelationalOperator::EQUAL, std::vector<ullint>(orbit.size(), 0LLU),
        "The Xoshiro256++ random number generator did not produce the expected bit strings when "
        "initialized with a specific state.");
  Xoshiro256ppGenerator xrs256ppc({ 0x7743a154e17a5e9bLLU, 0x7823a1cd9453899bLLU,
                                    0x976589eefbb1c7f5LLU, 0x702cf168260fa29eLLU });
  xrs256ppc.jump();
  const std::vector<ullint> orbit_ii = { 0x39c4396d8759c874LLU, 0x4b948d9de69752ecLLU,
                                         0x871591604b03d9a6LLU, 0x444d6d471322d17bLLU,
                                         0xb0a9eb9383bf0803LLU, 0x481f6c796c1d0ecaLLU,
                                         0xb89a346b480341bfLLU, 0x1494bad1d1b19126LLU,
                                         0xa2f5ca0a0ab0805LLU,  0x75a4de1da308cc8fLLU };
  std::vector<ullint> xrsc_result(orbit_ii.size());
  for (size_t i = 0; i < orbit_ii.size(); i++) {
    xrsc_result[i] = xrs256ppc.revealBitString();
    xrs256ppc.uniformRandomNumber();
  }
  std::vector<int> xrsc_diff(orbit.size());
  for (size_t i = 0; i < orbit.size(); i++) {
    xrsc_diff[i] = xrsc_result[i] - orbit_ii[i];
  }
  check(xrsc_diff, RelationalOperator::EQUAL, std::vector<ullint>(orbit_ii.size(), 0LLU),
        "The Xoshiro256++ random number generator did not produce the expected bit strings after "
        "traversing a jump.");
  
  // Verify rounding results
  section(3);
  const size_t szt_a = 159283;
  const size_t szt_a_round = roundUp<size_t>(szt_a, 32);
  check(szt_a_round, RelationalOperator::EQUAL, Approx(159296).margin(1.0e-6), "Rounding upwards "
        "to the nearest increment of 32 failed.");
  const int szt_b_round = roundUp<size_t>(szt_a, 128);
  check(szt_b_round, RelationalOperator::EQUAL, Approx(159360).margin(1.0e-6), "Rounding upwards "
        "to the nearest increment of 128 failed.");
  const ulint prime_composite = 2 * 2 * 5 * 3 * 7 * 19;
  const std::vector<uint> primes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 };
  const std::vector<uint> factors_i = primeFactorCounts(prime_composite, primes, 7);
  const std::vector<uint> factors_ii = primeFactorCounts(prime_composite, primes, 9);
  const std::vector<uint> answer_i = { 2, 1, 1, 1, 0, 0, 0 };
  const std::vector<uint> answer_ii = { 2, 1, 1, 1, 0, 0, 0, 1, 0 };
  check(factors_i, RelationalOperator::EQUAL, Approx(answer_i).margin(1.0e-8), "Prime "
        "factorization with numbers up to 7 fails.");
  const int near_factor_a = nearestFactor(128, 64, primes, LimitApproach::BELOW);
  const int near_factor_b = nearestFactor(128, 63, primes, LimitApproach::BELOW);
  const int near_factor_c = nearestFactor(128, 63, primes, LimitApproach::ABOVE);
  const int near_factor_d = nearestFactor(128, 69, primes, LimitApproach::ABOVE);
  const int near_factor_e = nearestFactor(12612600, 69, primes, LimitApproach::ABOVE);
  check(near_factor_a, RelationalOperator::EQUAL, 64, "The nearest factor of 128 to 64 was not "
        "found correctly when approaching from " + getEnumerationName(LimitApproach::BELOW) + ".");
  check(near_factor_b, RelationalOperator::EQUAL, 32, "The nearest factor of 128 to 63 was not "
        "found correctly when approaching from " + getEnumerationName(LimitApproach::BELOW) + ".");
  check(near_factor_c, RelationalOperator::EQUAL, 64, "The nearest factor of 128 to 63 was not "
        "found correctly when approaching from " + getEnumerationName(LimitApproach::ABOVE) + ".");
  check(near_factor_d, RelationalOperator::EQUAL, 128, "The nearest factor of 128 to 69 was not "
        "found correctly when approaching from " + getEnumerationName(LimitApproach::ABOVE) + ".");
  check(near_factor_e, RelationalOperator::EQUAL, 70, "The nearest factor of 12612600 to 69 was "
        "not found correctly when approaching from " + getEnumerationName(LimitApproach::ABOVE) +
        ".");
  const ullint big_product = ipowl(2, 14) * ipowl(3, 10) * ipowl(5, 6) * ipowl(7, 4);
  const int near_factor_f = nearestFactor(big_product, 77, primes, LimitApproach::BELOW, 4);
  check(near_factor_f, RelationalOperator::EQUAL, 75, "The nearest factor of a massive product "
        "of 2, 3, 5, and 7 to the target of 77 was not computed correctly when approaching "
        "from " + getEnumerationName(LimitApproach::ABOVE) + ".");
  const std::vector<int> fctr_comps = { factorial(0), factorial(1), factorial(4), factorial(6) };
  const std::vector<int> fctr_ans = { 1, 1, 24, 720 };
  check(fctr_comps, RelationalOperator::EQUAL, fctr_ans, "Factorials of small integers were not "
        "computed correctly.");
  CHECK_THROWS(factorial(-1), "A factorial of a negative number was computed.");
  CHECK_THROWS(factorial(13), "A factorial of an unwieldy number was computed.");
  
  // Check matrix math from the lightweight, onboard library
  section(4);
  const int rows_a = 7;
  const int cols_a = 5;
  const int rows_b = 5;
  const int cols_b = 7;
  Hybrid<double> tr_mat_a(rows_a * cols_a, "random matrix");
  Hybrid<double> tr_mat_b(rows_b * cols_b, "random matrix");
  Hybrid<double> tr_mat_c(rows_a * cols_b, "posdef matrix");
  Hybrid<double> tr_mat_d(cols_a * rows_b, "posdef matrix");
  double* tr_a_ptr = tr_mat_a.data();
  double* tr_b_ptr = tr_mat_b.data();
  for (int i = 0; i < rows_a * cols_a; i++) {
    tr_a_ptr[i] = prng_c.gaussianRandomNumber();
  }
  for (int i = 0; i < rows_b * cols_b; i++) {
    tr_b_ptr[i] = prng_c.gaussianRandomNumber();
  }

  // Tranpose the matrices in-place and out-of-place.
  std::vector<double> tr_mat_a_copy = tr_mat_a.readHost();
  std::vector<double> tr_mat_b_copy = tr_mat_b.readHost();
  std::vector<double> tr_mat_c_copy = tr_mat_c.readHost();
  std::vector<double> tr_mat_a_xpse(tr_mat_a_copy.size());
  std::vector<double> tr_mat_b_xpse(tr_mat_b_copy.size());
  std::vector<double> tr_mat_c_xpse(tr_mat_c_copy.size());
  transpose(tr_mat_a_copy, &tr_mat_a_xpse, rows_a, cols_a);
  transpose(tr_mat_b_copy, &tr_mat_b_xpse, rows_b, cols_b);
  transpose(tr_mat_c_copy, &tr_mat_c_xpse, rows_a, cols_b);
  transpose(&tr_mat_a_copy, rows_a, cols_a);
  transpose(&tr_mat_b_copy, rows_b, cols_b);
  transpose(&tr_mat_c_copy, rows_a, cols_b);
  check(tr_mat_a_copy, RelationalOperator::EQUAL, tr_mat_a_xpse, "In-place transposition of a "
        "non-square matrix does not yield the same result as out-of-place transposition.");
  check(tr_mat_b_copy, RelationalOperator::EQUAL, tr_mat_b_xpse, "In-place transposition of a "
        "non-square matrix does not yield the same result as out-of-place transposition.");
  check(tr_mat_c_copy, RelationalOperator::EQUAL, tr_mat_c_xpse, "In-place transposition of a "
        "square matrix does not yield the same result as out-of-place transposition.");
  transpose(&tr_mat_a_copy, cols_a, rows_a);
  transpose(&tr_mat_b_copy, cols_b, rows_b);
  transpose(&tr_mat_c_copy, cols_b, rows_a);

  // Check that the GEMM matrix multiplication snapshot file exists
  const std::string matrices_snp = oe.getStormmSourcePath() + osSeparator() + "test" +
                                   osSeparator() + "Math" + osSeparator() + "matrices.m";
  snp_found = (getDrivePathType(matrices_snp) == DrivePathType::FILE) ? TestPriority::CRITICAL :
                                                                        TestPriority::ABORT;
  if (snp_found == TestPriority::ABORT && oe.takeSnapshot() != SnapshotOperation::SNAPSHOT) {
    rtWarn("The snapshot file " + matrices_snp + "was not found.  Make sure that the "
           "$STORMM_SOURCE environment variable is set to the root soruce directory, where src/ "
           "and test/ subdirectories can be found.  A number of subsequent tests will be skipped.",
           "test_math");
  }

  // Try a direct matrix-matrix multiplication, producing a non-symmetric [7 by 7] matrix result
  matrixMultiply(tr_mat_a.data(), rows_a, cols_a, tr_mat_b.data(), rows_b, cols_b,
                 tr_mat_c.data());
  snapshot(matrices_snp, polyNumericVector(tr_mat_c.readHost()), "tr_ab", 1.0e-5, "Matrix "
           "multiply result for [7 by 5] x [5 by 7] matrices is incorrect.", oe.takeSnapshot(),
           1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::OVERWRITE, snp_found);
  
  // Overwrite the result with a positive definite [7 by 7] matrix result
  matrixMultiply(tr_mat_a.data(), rows_a, cols_a, tr_mat_a.data(), rows_a, cols_a, tr_mat_c.data(),
                 1.0, 1.0, 0.0, TransposeState::AS_IS, TransposeState::TRANSPOSE);
  snapshot(matrices_snp, polyNumericVector(tr_mat_c.readHost()), "tr_aat", 1.0e-5, "Matrix "
           "multiply result for [7 by 5] x [7 by 5](T) matrices is incorrect.", oe.takeSnapshot(),
           1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, snp_found);

  // Obtain a new, positive definite [5 by 5] matrix result
  matrixMultiply(tr_mat_a.data(), rows_a, cols_a, tr_mat_a.data(), rows_a, cols_a, tr_mat_d.data(),
                 1.0, 1.0, 0.0, TransposeState::TRANSPOSE, TransposeState::AS_IS);
  snapshot(matrices_snp, polyNumericVector(tr_mat_d.readHost()), "tr_ata", 1.0e-5, "Matrix "
           "multiply result for [7 by 5](T) x [7 by 5] matrices is incorrect.", oe.takeSnapshot(),
           1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, snp_found);

  // Try additional GEMM operations
  checkMatrixMultiply( 5,  5,  5,  3);
  checkMatrixMultiply( 5,  5,  5,  3, TransposeState::TRANSPOSE, TransposeState::AS_IS);
  checkMatrixMultiply( 6,  3,  3,  3);
  checkMatrixMultiply( 6,  3,  3,  3, TransposeState::AS_IS, TransposeState::TRANSPOSE);
  checkMatrixMultiply( 6,  3,  3,  3);
  checkMatrixMultiply( 6,  7,  3,  6, TransposeState::TRANSPOSE, TransposeState::TRANSPOSE);
  checkMatrixMultiply( 5,  9,  5,  3, TransposeState::TRANSPOSE, TransposeState::AS_IS);
  checkMatrixMultiply( 5,  9, 37,  5, TransposeState::TRANSPOSE, TransposeState::TRANSPOSE);
  checkMatrixMultiply(17,  5, 37,  5, TransposeState::AS_IS, TransposeState::TRANSPOSE);

  // Compute the Moore-Penrose matrix pseudo-inverse of one of the matrices obtained earlier.
  takePseudoInverse( 7,  5, matrices_snp, "pinv_a", oe, snp_found);
  takePseudoInverse( 8,  6, matrices_snp, "pinv_b", oe, snp_found);
  takePseudoInverse( 6,  8, matrices_snp, "pinv_c", oe, snp_found);
  takePseudoInverse( 8,  6, matrices_snp, "pinv_d", oe, snp_found);
  takePseudoInverse(15, 11, matrices_snp, "pinv_e", oe, snp_found);
  takePseudoInverse(35, 49, matrices_snp, "pinv_f", oe, snp_found);

  // Try solving some overdetermined systems of linear equations using the left pseudo-inverse
  // method, then the QR decomposition.
  solveLinearSystem( 7,  5);
  solveLinearSystem(17,  7);
  solveLinearSystem(12, 12);
  solveLinearSystem(39, 22);
  
  // Try inverting a positive definite matrix
  Hybrid<double> tr_mat_e(cols_a * rows_b, "inverse matrix");
  invertSquareMatrix(tr_mat_d.data(), tr_mat_e.data(), cols_a);
  snapshot(matrices_snp, polyNumericVector(tr_mat_d.readHost()), "inv_ata", 1.0e-5, "Matrix "
           "inversion result for a [5 by 5] matrix is incorrect.", oe.takeSnapshot(),
           1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, snp_found);

  // Overwrite the positive definite matrix with a non-symmetric [5 by 5] matrix result
  matrixMultiply(tr_mat_a.data(), rows_a, cols_a, tr_mat_b.data(), rows_b, cols_b, tr_mat_d.data(),
                 1.0, 1.0, 0.0, TransposeState::TRANSPOSE, TransposeState::TRANSPOSE);
  snapshot(matrices_snp, polyNumericVector(tr_mat_d.readHost()), "tr_atbt", 1.0e-5, "Matrix "
           "multiply result for [7 by 5](T) x [5 by 7](T) matrices is incorrect.",
           oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::APPEND,
           snp_found);
  
  // Create a positive-definite matrix, then compute its eigenvalues and eigenvectors (this is
  // better accomplished by a routine like BLAS dsyevd, as it is more than just real and symmetric,
  // but the slower jacobi routine is all STORMM has got without real BLAS compiled).
  const size_t rank = 8;
  Hybrid<double> base_matrix(rank * rank);
  Hybrid<double> posdef_matrix(rank * rank);
  Hybrid<double> eigenvectors(rank * rank);
  Hybrid<double> eigenvalues(rank);
  double* dbase = base_matrix.data();
  double* dposdef = posdef_matrix.data();
  for (size_t i = 0; i < rank * rank; i++) {
    dbase[i] = prng_c.gaussianRandomNumber();
  }
  matrixMultiply(dbase, rank, rank, dbase, rank, rank, dposdef, 1.0, 1.0, 0.0,
                 TransposeState::TRANSPOSE, TransposeState::AS_IS);
  std::vector<double> posdef_mat_b = posdef_matrix.readHost();
  std::vector<double> posdef_mat_c = posdef_matrix.readHost();
  jacobiEigensolver(&posdef_matrix, &eigenvectors, &eigenvalues, rank);
  snapshot(matrices_snp, polyNumericVector(eigenvectors.readHost()), "eigvec", 1.0e-5,
           "Eigenvectors for a rank-8 positive definite matrix are incorrect.", oe.takeSnapshot(),
           1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, snp_found);
  snapshot(matrices_snp, polyNumericVector(eigenvalues.readHost()), "eigval", 1.0e-5,
           "Eigenvalues for a rank-8 positive definite matrix are incorrect.", oe.takeSnapshot(),
           1.0e-8, NumberFormat::STANDARD_REAL, PrintSituation::APPEND, snp_found);
  std::vector<double> rsym_diag_b(8, 0.0);
  std::vector<double> rsym_eigv_b(8, 0.0);
  realSymmEigensolver(posdef_mat_b.data(), rank, rsym_eigv_b.data(), rsym_diag_b.data());
  bool eigval_match = true;
  bool eigvec_match = true;
  std::vector<bool> eig_taken(rank, false);
  std::vector<double> eigvec_mults(rank);
  const double* evec_ptr = eigenvectors.data();
  for (int i = 0; i < rank; i++) {
    bool val_found = false;
    bool vec_found = false;
    for (int j = 0; j < rank; j++) {
      if (eig_taken[j]) {
        continue;
      }
      if (fabs(rsym_eigv_b[i] - eigenvalues.readHost(j)) < 1.0e-7) {
        eig_taken[j] = true;
        val_found = true;
        for (int k = 0; k < rank; k++) {
          eigvec_mults[k] = eigenvectors.readHost((rank * j) + k) / posdef_mat_b[(rank * i) + k];
        }
        vec_found = (variance(eigvec_mults,
                              VarianceMethod::ROOT_MEAN_SQUARED_DEVIATION) < 1.0e-6);
      }
    }
    eigval_match = (eigval_match && val_found);
    eigvec_match = (eigvec_match && vec_found);
  }
  check(eigval_match, "Eigenvalues computed by the symmetric real eigensolver do not agree with "
        "those computed using the Jacobi solver.");
  check(eigvec_match, "Eigenvectors computed by the symmetric real eigensolver do not agree with "
        "those computed using the Jacobi solver.");
  std::vector<double> rsym_diag_c(8, 0.0);
  std::vector<double> rsym_eigv_c(8, 0.0);
  realSymmEigensolver(posdef_mat_c.data(), rank, rsym_eigv_c.data(), rsym_diag_c.data(),
                      EigenWork::EIGENVALUES);
  check(rsym_eigv_b, RelationalOperator::EQUAL, Approx(rsym_eigv_c).margin(1.0e-8), "Eigenvalues "
        "obtained from the real symmetric solver are not the same when computed in the absence of "
        "eigenvectors.h");
  
  // Try a much bigger eigenvalue problem and check its results
  const size_t big_rank = 95;
  HpcMatrix<double> mtrx_base(big_rank, big_rank, MatrixFillMode::RANDN, &prng_c);
  HpcMatrix<double> mtrx_posdef(big_rank, big_rank);
  HpcMatrix<double> mtrx_eigvec(big_rank, big_rank);
  eigenvalues.resize(big_rank);
  matrixMultiply(mtrx_base, mtrx_base, &mtrx_posdef, 1.0, 1.0, 0.0, TransposeState::TRANSPOSE,
                 TransposeState::AS_IS);
  dposdef = mtrx_posdef.memptr();
  double bad_value = mtrx_posdef(big_rank / 2, 0) * 1.098;
  std::swap(bad_value, dposdef[big_rank / 2]);
  CHECK_THROWS(jacobiEigensolver(&mtrx_posdef, &mtrx_eigvec, &eigenvalues, ExceptionResponse::DIE),
               "The jacobiEigensolver() routine attempted to work on a non-symmetric matrix.");
  std::swap(bad_value, dposdef[big_rank / 2]);
  HpcMatrix<double> copy_posdef = mtrx_posdef;
  jacobiEigensolver(&mtrx_posdef, &mtrx_eigvec, &eigenvalues, ExceptionResponse::DIE);
  std::vector<double> eigvsums(big_rank, 0.0);
  Hybrid<double> mtrx_eigtest(big_rank, "test_eig");
  for (size_t i = 0; i < big_rank; i++) {
    Hybrid<double> icol_ptr = mtrx_eigvec.colptr(i);
    matrixVectorMultiply(copy_posdef, icol_ptr, &mtrx_eigtest, 1.0, 1.0, 0.0,
                         TransposeState::AS_IS);
    eigvsums[i] = sum<double>(mtrx_eigtest) - sum<double>(icol_ptr) * eigenvalues.readHost(i);
  }
  check(eigvsums, RelationalOperator::EQUAL,
        Approx(std::vector<double>(big_rank, 0.0)).margin(1.0e-8), "Eigenvalues and eigenvectors "
        "produced by jacobiEigensolver() are incorrect.");

  // Try computing box transformation matrices
  const double lx = 64.1;
  const double ly = 60.3;
  const double lz = 48.9;
  const double alpha =  98.7 * pi / 180.0;
  const double beta  = 103.2 * pi / 180.0;
  const double gamma =  85.4 * pi / 180.0;
  HpcMatrix<double> umat(3, 3);
  HpcMatrix<double> invu(3, 3);
  HpcMatrix<double> xfrm_prod(3, 3);
  computeBoxTransform(lx, ly, lz, alpha, beta, gamma, umat.memptr(), invu.memptr());
  matrixMultiply(umat, invu, &xfrm_prod);
  HpcMatrix<double> ident(3, 3, MatrixFillMode::EYE);
  std::vector<double> linear_umat(9, 0.0);
  std::vector<double> linear_invu(9, 0.0);
  std::vector<double> linear_xfrm_prod(9, 0.0);
  std::vector<double> linear_ident(9, 0.0);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      linear_umat[(j * 3) + i] = umat(i, j);
      linear_invu[(j * 3) + i] = invu(i, j);
      linear_xfrm_prod[(j * 3) + i] = xfrm_prod(i, j);
      linear_ident[(j * 3) + i] = ident(i, j);
    }
  }

  // Data entered in column-major format may look like a transpose at first.  Check
  // the box transformation matrices against results from an external implementation.
  const std::vector<double> umat_answer = {
    1.5600624025e-02,  0.0000000000e+00,  0.0000000000e+00,
   -1.2551964058e-03,  1.6637338818e-02,  0.0000000000e+00,
    3.5203270719e-03,  2.3009524689e-03,  2.1204798362e-02 };
  const std::vector<double> invu_answer = {
    6.4100000000e+01,  0.0000000000e+00,  0.0000000000e+00,
    4.8359951370e+00,  6.0105766371e+01,  0.0000000000e+00,
   -1.1166357548e+01, -6.5221328286e+00,  4.7159137423e+01 };
  check(linear_umat, RelationalOperator::EQUAL, Approx(umat_answer).margin(1.0e-8),
        "The box transformation matrix is incorrect.");
  check(linear_invu, RelationalOperator::EQUAL, Approx(invu_answer).margin(1.0e-8),
        "The inverse transformation matrix is incorrect.");
  check(linear_xfrm_prod, RelationalOperator::EQUAL, linear_ident, "The product of box "
        "transformation and inverse transformation matrices is not the identity matrix.");

  // Check the computation of unit cell widths: computing the Hessian normal form
  const std::vector<double> invu_s1 = { 2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0 };
  std::vector<double> s1_widths(3);
  const std::vector<double> s1_widths_ans = { 2.0, 3.0, 4.0 };
  double* s1_ptr = s1_widths.data();
  hessianNormalWidths(invu_s1.data(), &s1_ptr[0], &s1_ptr[1], &s1_ptr[2]);
  check(s1_widths, RelationalOperator::EQUAL, s1_widths_ans, "The Hessian normal form is not "
        "correctly computed for an orthorhombic unit cell.");
  const std::vector<double> invu_s2 = { 100.000000000000,   0.000000000000,   0.000000000000,
                                        -33.333333276578,  94.280904178272,   0.000000000000,
                                        -33.333333276578, -47.140451968741,  81.649658179660 };
  std::vector<double> s2_widths(3);
  const std::vector<double> s2_widths_ans = { invu_s2[8], invu_s2[8], invu_s2[8] };
  double* s2_ptr = s2_widths.data();
  hessianNormalWidths(invu_s2.data(), &s2_ptr[0], &s2_ptr[1], &s2_ptr[2]);
  check(s2_widths, RelationalOperator::EQUAL, Approx(s2_widths_ans).margin(1.0e-7), "The Hessian "
        "normal form is not correctly computed for an orthorhombic unit cell.");  

  // Check Leibniz' formula for computing determinants
  const std::vector<double> two_mat   = { -1.0, 0.5, 7.4, 8.3 };
  const std::vector<double> three_mat = { 1.1, 0.8, 3.9, 5.6, 7.7, 9.4, 2.0, -4.3, -0.7 };
  const std::vector<double> four_mat   = {  6.1,  7.9, -1.2, -5.5, -0.2,  5.4,  3.3,  9.1,
                                            3.7,  2.0, -5.0, -8.3,  8.6,  1.9,  9.5,  0.4 };
  std::vector<double> singular_five(25);
  for (int i = 0; i < 20; i++) {
    singular_five[i] = (xrs128p.uniformRandomNumber() * 10.0) - 5.0;
  }
  for (int i = 20; i < 25; i++) {
    singular_five[i] = singular_five[i - 20];
  }
  check(leibnizDeterminant(two_mat), RelationalOperator::EQUAL, Approx(-12.0).margin(tiny),
        "The determinant of a 2 x 2 matrix was not computed correctly.");
  check(leibnizDeterminant(three_mat), RelationalOperator::EQUAL, Approx(-97.263).margin(tiny),
        "The determinant of a 3 x 3 matrix was not computed correctly by the Leibniz formula.");
  check(leibnizDeterminant(four_mat), RelationalOperator::EQUAL, Approx(-1767.2249).margin(tiny),
        "The determinant of a 4 x 4 matrix was not computed correctly by the Leibniz formula.");
  check(leibnizDeterminant(singular_five), RelationalOperator::EQUAL, Approx(0.0).margin(tiny),
        "The determinant of a singular 5 x 5 matrix was not computed correctly by the Leibniz "
        "formula.");

  // Further checks on the mean, standard deviation, variance, correlation, dot product,
  // magnitude, and projection operations
  section(1);
  check(mean(invu_answer), RelationalOperator::EQUAL, Approx(17.6124898394).margin(tiny),
        "The mean value of a double-precision vector is not correct.");
  check(variance(invu_answer, VarianceMethod::VARIANCE), RelationalOperator::EQUAL,
        Approx(816.0346461040).margin(tiny), "The variance of a double-precision vector is not "
        "correct.");
  check(variance(invu_answer, VarianceMethod::STANDARD_DEVIATION), RelationalOperator::EQUAL,
        Approx(30.2991580224).margin(tiny), "The standard deviation of a double-precision vector "
        "is not correct.");
  check(variance(invu_answer, VarianceMethod::ROOT_MEAN_SQUARED_DEVIATION),
        RelationalOperator::EQUAL, Approx(28.5663201359).margin(tiny), "The root mean squared "
        "deviation of a double-precision vector is not correct.");
  check(variance(invu_answer, VarianceMethod::COEFFICIENT_OF_VARIATION),
        RelationalOperator::EQUAL, Approx(1.7203222428).margin(tiny), "The coefficient of "
        "variation for a double-precision vector is not correct.");
  check(variance(invu_answer, VarianceMethod::NORMALIZED_RMSD),
        RelationalOperator::EQUAL, Approx(1.6219353650).margin(tiny), "The normalized root mean "
        "squared deviation for a double-precision vector is not correct.");
  check(magnitude(umat_answer), RelationalOperator::EQUAL, Approx(3.1449747049e-02).margin(tiny),
        "The magnitude of a double-precision vector is not computed correctly.");
  check(dot(invu_answer, umat_answer), RelationalOperator::EQUAL,
        Approx(2.9396135279).margin(tiny), "The dot product of two double-precision vectors is "
        "not computed correctly.");
  std::vector<double> invu_perturb(invu_answer);
  std::vector<double> invu_partial(invu_answer);
  for (size_t i = 6; i < invu_partial.size(); i++) {
    invu_partial[i] = 0.0;
  }
  for (size_t i = 0; i < invu_perturb.size(); i++) {
    invu_perturb[i] += xrs256pp.uniformRandomNumber() - 0.5;
    invu_partial[i] += xrs256pp.uniformRandomNumber() - 0.5;
  }
  std::vector<double> invu_project(invu_answer.size());
  project(invu_answer, invu_partial, &invu_project);
  const std::vector<double> invu_project_answer = { 63.5808430141, -0.2998731521,  0.4369816989,
                                                     5.2640142928, 60.3620556749,  0.0677375859,
                                                    -0.2593021730, -0.1475558457, -0.3946802908 };
  check(invu_project, RelationalOperator::EQUAL, Approx(invu_project_answer).margin(tiny),
        "The projection of one vector onto another does not meet expectations.");
  const std::vector<double> corr_x1 = {  0.5, -1.0, -5.4,  6.7,  7.2,  3.8, -4.1, 9.3 };
  const std::vector<double> corr_x2 = {  0.7, -0.8, -4.9,  9.2,  5.3,  4.0, -5.9, 7.9 };
  check(pearson(corr_x1, corr_x2), RelationalOperator::EQUAL, Approx(0.9651506029).margin(tiny),
        "The correlation coefficient computed for two vectors is incorrect.");

  // Check single-precision random number generation
  section(2);
  Xoroshiro128pGenerator test_128;
  Xoroshiro128pGenerator test_128f;
  Xoshiro256ppGenerator test_256;
  Xoshiro256ppGenerator test_256f;
  const size_t npts = 64;
  std::vector<double> t128_uni(npts), t128f_uni(npts), t128_gss(npts), t128f_gss(npts);
  std::vector<double> t256_uni(npts), t256f_uni(npts), t256_gss(npts), t256f_gss(npts);
  for (size_t i = 0; i < npts; i++) {
    t128_uni[i]  = test_128.uniformRandomNumber();
    t128f_uni[i] = test_128f.spUniformRandomNumber();
    t128_gss[i]  = test_128.gaussianRandomNumber();
    t128f_gss[i] = test_128f.spGaussianRandomNumber();
    t256_uni[i]  = test_256.uniformRandomNumber();
    t256f_uni[i] = test_256f.spUniformRandomNumber();
    t256_gss[i]  = test_256.gaussianRandomNumber();
    t256f_gss[i] = test_256f.spGaussianRandomNumber();
  }
  check(t128f_uni, RelationalOperator::EQUAL, Approx(t128_uni).margin(3.0e-7), "Random numbers "
        "pulled from a uniform distribution with the Xoroshiro128+ generator do not have the same "
        "values when produced in single- versus double-precision.");
  check(t128f_gss, RelationalOperator::EQUAL, Approx(t128_gss).margin(6.0e-6), "Random numbers "
        "pulled from a normal distribution with the Xoroshiro128+ generator do not have the same "
        "values when produced in single- versus double-precision.");
  check(t256f_uni, RelationalOperator::EQUAL, Approx(t256_uni).margin(3.0e-7), "Random numbers "
        "pulled from a uniform distribution with the Xoshiro256++ generator do not have the same "
        "values when produced in single- versus double-precision.");
  check(t256f_gss, RelationalOperator::EQUAL, Approx(t256_gss).margin(6.0e-6), "Random numbers "
        "pulled from a normal distribution with the Xoshiro256++ generator do not have the same "
        "values when produced in single- versus double-precision.");

  // Check the behavior of the Xoroshiro128+ series object
  Xoroshiro128pGenerator leader_xrs128p(9183084, 75);
  const ullint2 init_state = leader_xrs128p.revealState();
  Xoroshiro128pGenerator test_xrs128p(105892, 95);
  test_xrs128p.setState(init_state);
  std::vector<double> orig_uni_output(npts), rstr_uni_output(npts);
  for (int i = 0; i < npts; i++) {
    orig_uni_output[i] = leader_xrs128p.uniformRandomNumber();
    rstr_uni_output[i] = test_xrs128p.uniformRandomNumber();
  }
  check(rstr_uni_output, RelationalOperator::EQUAL, orig_uni_output, "Resetting the state of a "
        "Xoroshiro128+ generator failed to restart the sequence as expected.");
  const int ngen = 1024;
  const int bank_depth = 8;
  RandomNumberMill<double> my_128p_series(init_state, ngen, bank_depth);
  CHECK_THROWS(my_128p_series.getBankValue(ngen + 3, 7), "Invalid random number bank access was "
               "permitted in a Xoroshiro128+ generator series object.");
  int mseries_state = 0;
  const int stride_length = my_128p_series.getRefreshStride();
  check(stride_length, RelationalOperator::EQUAL, ngen / bank_depth, "The length of a refresh "
        "stride in the Xoroshiro128+ generator series should be such that running through a "
        "number of refresh cycles equal to the depth of the object's banks should cover all "
        "generators.");
  const int sample_a = 0;
  const int sample_b = ngen - 53;
  const int sample_c = (ngen / 2) - 37;
  const int sample_d = ngen / 4;
  const int maxcyc = 4;
  std::vector<double> a_path(bank_depth * maxcyc), b_path(bank_depth * maxcyc);
  std::vector<double> c_path(bank_depth * maxcyc), d_path(bank_depth * maxcyc);
  for (int cyc = 0; cyc < maxcyc; cyc++) {
    for (int i = 0; i < bank_depth; i++) {
      a_path[(cyc * bank_depth) + i] = my_128p_series.getBankValue(sample_a, i);
      b_path[(cyc * bank_depth) + i] = my_128p_series.getBankValue(sample_b, i);
      c_path[(cyc * bank_depth) + i] = my_128p_series.getBankValue(sample_c, i);
      d_path[(cyc * bank_depth) + i] = my_128p_series.getBankValue(sample_d, i);
      my_128p_series.gaussianRandomNumbers(i * stride_length, (i + 1) * stride_length);
    }
  }
  const std::vector<double> a_expect = generateExpectedSeries(&test_xrs128p, init_state, sample_a,
                                                              ngen, bank_depth, maxcyc);
  const std::vector<double> b_expect = generateExpectedSeries(&test_xrs128p, init_state, sample_b,
                                                              ngen, bank_depth, maxcyc);
  const std::vector<double> c_expect = generateExpectedSeries(&test_xrs128p, init_state, sample_c,
                                                              ngen, bank_depth, maxcyc);
  const std::vector<double> d_expect = generateExpectedSeries(&test_xrs128p, init_state, sample_d,
                                                              ngen, bank_depth, maxcyc);
  check(a_path, RelationalOperator::EQUAL, a_expect, "A random number sequence obtained from "
        "generator " + std::to_string(sample_a) + " of a Xoroshiro128+ generator series did not "
        "meet expectations.");
  check(b_path, RelationalOperator::EQUAL, b_expect, "A random number sequence obtained from "
        "generator " + std::to_string(sample_b) + " of a Xoroshiro128+ generator series did not "
        "meet expectations.");
  check(c_path, RelationalOperator::EQUAL, c_expect, "A random number sequence obtained from "
        "generator " + std::to_string(sample_c) + " of a Xoroshiro128+ generator series did not "
        "meet expectations.");
  check(d_path, RelationalOperator::EQUAL, d_expect, "A random number sequence obtained from "
        "generator " + std::to_string(sample_d) + " of a Xoroshiro128+ generator series did not "
        "meet expectations.");

  // Create a tick counter with a series of dials with different lengths.
  section(5);
  TickCounter<int> dials({ 2, 2, 3, 5, 7, 4, 1, 2 });
  testDialCounter(&dials, 5000);
  dials.reset();
  testDialCounter(&dials, 2134);
  dials.reset();
  testDialCounter(&dials,  187);
  CHECK_THROWS(TickCounter<int> test_dials({ 2, 0, 7}), "A TickCounter object was created with "
              "zero increments possible in one of the variables.");
  CHECK_THROWS(TickCounter<int> test_dials({ 9, 5, -7}), "A TickCounter object was created with "
               "negative increments in one of the variables.");
  CHECK_THROWS(dials.set({ 0, 0, 1, 1, 5, 3, 0}), "A TickCounter allowed a misaligned vector "
               "to be used in setting its state.");
  CHECK_THROWS(dials.set({ 0, 0, -1, 1, 5, 3, 0, 0}), "A TickCounter allowed a negative number to "
               "enter its state settings.");
  dials.loadStateValues({ 0, 1 }, 0);
  dials.loadStateValues({ 7, 5 }, 1);
  dials.loadStateValues({ 1, 2, 4 }, 2);
  dials.loadStateValues({ 3, 2, 4, 7, 7 }, 3);
  dials.loadStateValues({ 6, 8, 3, 2, 4, 7, 7 }, 4);
  dials.loadStateValues({ 4, 5, 7, 8 }, 5);
  CHECK_THROWS(dials.loadStateValues({ 4, 5, 7, 8 }, 6), "A TickCounter accepted a series of four "
               "values for a variable with only one setting.");
  dials.loadStateValues({ 3 }, 6);
  dials.loadStateValues({ 3, 3 }, 7);
  dials.reset();
  const std::vector<int> zero_state = dials.getState();
  const std::vector<int> zero_state_ans = { 0, 7, 1, 3, 6, 4, 3, 3 };
  check(zero_state, RelationalOperator::EQUAL, zero_state_ans, "A TickCounter does not return the "
        "expected state after a reset operation.");
  dials.advance();
  const std::vector<int> one_state = dials.getState();
  const std::vector<int> one_state_ans = { 1, 7, 1, 3, 6, 4, 3, 3 };
  check(one_state, RelationalOperator::EQUAL, one_state_ans, "A TickCounter does not return the "
        "expected state after a advancing one tick.");
  dials.advance(8);
  const std::vector<int> nine_state = dials.getState();
  const std::vector<int> nine_state_ans = { 1, 7, 4, 3, 6, 4, 3, 3 };
  check(nine_state, RelationalOperator::EQUAL, nine_state_ans, "A TickCounter does not return the "
        "expected state after a advancing one tick.");
  TickCounter<double> k_dials({ 3, 3, 6 });
  loadScalarStateValues(&k_dials, std::vector<double2>({ {0.0, 4.5}, {5.7, 7.8}, {4.5, 12.6} }));
  const std::vector<double> dzero_state = k_dials.getState();
  const std::vector<double> dzero_state_ans = { 0.0, 5.7, 4.5 };
  check(dzero_state, RelationalOperator::EQUAL, dzero_state_ans, "A TickCounter does not return "
        "the expected state just after construction.");
  k_dials.advance(10);
  const std::vector<double> dten_state = k_dials.getState();
  const std::vector<double> dten_state_ans = { 1.5, 5.7, 5.85 };
  check(dten_state, RelationalOperator::EQUAL, dten_state_ans, "A TIckCounter does not return "
        "the expected state after some advancement.");

  // Check tricubic interpolation mechanics.  Begin with the high-continuity stencil mode.
  section(6);
  const TricubicStencil tcstencil_cont(Interpolant::SMOOTHNESS);
  const TricubicStencil tcstencil_accr(Interpolant::FUNCTION_VALUE);
  bool all_integer = true;
  const double* tccont_ptr = tcstencil_cont.data();
  for (int i = 0; i < 4096; i++) {
    all_integer = (all_integer && Approx(tccont_ptr[i]).test(round(tccont_ptr[i])));
  }
  check(all_integer, "The tricubic spline coefficients matrix contains non-integral elements.");
  std::vector<double> tccol_sums(64, 0.0);
  std::vector<double> tccol_sums_ans(64, 0.0);
  tccol_sums_ans[7] = 1.0;
  for (int i = 0; i < 64; i++) {
    for (int j = 0; j < 64; j++) {
      tccol_sums[i] += tccont_ptr[(64 * i) + j];
    }
  }
  check(tccol_sums, RelationalOperator::EQUAL, tccol_sums_ans, "The column sums of the tricubic "
        "spline coefficients matrix do not meet expectations.");
  tricubicTestBundle(tcstencil_cont);
  std::vector<double> random_coefficients(64);
  for (int i = 0; i < 64; i++) {
    random_coefficients[i] = xrs256pp.gaussianRandomNumber();
  }
  tricubicTestBundle(tcstencil_cont, random_coefficients);
  std::vector<double> x_only_coefficients(64, 0.0);
  for (int i = 0; i < 4; i++) {
    x_only_coefficients[i] = xrs256pp.gaussianRandomNumber();
  }
  tricubicTestBundle(tcstencil_cont, x_only_coefficients, { 0.5, 1.0, 1.0 });
  const std::vector<double> equivariant = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
  const std::vector<double> equivariant_ii = { 0.9, 0.0, 0.0, 0.0, 0.8, 0.0, 0.0, 0.0, 0.7 };
  const std::vector<double> triclinic_i = {  0.50000000000,  0.00000000000,  0.00000000000,
                                            -0.33333333277,  0.94280904178,  0.00000000000,
                                            -0.29999999949, -0.42426406772,  0.73484692362 };
  const std::vector<double> triclinic_ii = {  0.50000000000,  0.00000000000,  0.00000000000, 
                                             -0.23754420664,  0.97137672913,  0.00000000000, 
                                             -0.25721634619, -0.37174068675,  0.77823429189 };
  tricubicTestBundle(tcstencil_cont, x_only_coefficients, triclinic_i);
  tricubicTestBundle(tcstencil_cont, random_coefficients,
                     { 0.7, 0.9, 1.1 });
  testTricubicCorners(equivariant_ii);
  testTricubicCorners(triclinic_i);

  // Check mesh mechanics with the high-accuracy stencil mode.
  tricubicTestBundle(tcstencil_accr);
  tricubicTestBundle(tcstencil_accr, random_coefficients);
  tricubicTestBundle(tcstencil_accr, x_only_coefficients, { 0.5, 1.0, 1.0 });
  tricubicTestBundle(tcstencil_accr, x_only_coefficients, triclinic_i);
  tricubicTestBundle(tcstencil_accr, random_coefficients, { 0.7, 0.9, 1.1 });

  // Compate the two modes for the tricubic interpolant.
  testCoulombInterpolation(equivariant, tcstencil_cont, &xrs256pp);
  testCoulombInterpolation(equivariant, tcstencil_accr, &xrs256pp);
  testCoulombInterpolation(triclinic_i, tcstencil_cont, &xrs256pp);
  testCoulombInterpolation(triclinic_i, tcstencil_accr, &xrs256pp);

  // Check customized formulas in STORMM's math library
  testSigmoid(4.3, 2.6);
  check(ipow(5, 5), RelationalOperator::EQUAL, 3125, "The integer power function does not produce "
        "the correct result.");
  check(ipow(1, 20853835), RelationalOperator::EQUAL, 1, "The integer power function does not "
        "produce the correct result.");
  check(ipow(-1, 9816837), RelationalOperator::EQUAL, -1, "The integer power function does not "
        "produce the correct result.");
  check(ipow(-2, 14), RelationalOperator::EQUAL, 16384, "The integer power function does not "
        "produce the correct result.");
  check(ipowl(5, 5), RelationalOperator::EQUAL, 3125, "The integer power function does not "
        "produce the correct result.");
  check(ipowl(-7, 13), RelationalOperator::EQUAL, -96889010407, "The long long integer power "
        "function does not produce the correct result.");
  check(ipowl(2, 37), RelationalOperator::EQUAL, 137438953472, "The long long integer power "
        "function does not produce the correct result.");
  check(ipowl(-1, 30683237), RelationalOperator::EQUAL, -1, "The long long integer power "
        "function does not produce the correct result.");
  CHECK_THROWS(ipow(6, -7), "The integer power function accepted a negative exponent.");
  CHECK_THROWS(ipowl(3, -29), "The long long integer power function accepted a negative "
               "exponent.");
  CHECK_THROWS(bSpline<float>(0.5, 1), "B-spline computations were permitted with order 1.");
  CHECK_THROWS(bSpline<double>(1.01, 3), "B-spline computations were permitted using a coordinate "
               "outside the allowed range of [0.0, 1.0].");
  const std::vector<double> bspl2 = bSpline<double>(0.5, 2);
  check(bspl2, RelationalOperator::EQUAL, std::vector<double>(2, 0.5), "B-spline computation "
        "for the hat function with coordinate 0.5 did not yield the expected results.");
  const std::vector<double> bspl3 = bSpline<double>(0.5, 3);
  const std::vector<double> bspl3_ans = { 0.125, 0.750, 0.125 };
  check(bspl3, RelationalOperator::EQUAL, bspl3_ans, "B-spline computation for the parabolic "
        "function with coordinate 0.5 did not yield the expected results.");
  const std::vector<double> bspl3_asm = bSpline<double>(0.15, 3);
  const std::vector<double> bspl3_asm_ans = { 0.361250000, 0.627500000, 0.011250000 };
  check(bspl3_asm, RelationalOperator::EQUAL, bspl3_asm_ans, "B-spline computation for the "
        "parabolic function with coordinate 0.15 did not yield the expected results.");
  const std::vector<double> bspl4 = bSpline<double>(0.5, 4);
  const std::vector<double> bspl4_ans = { 0.020833333, 0.479166667, 0.479166667, 0.020833333 };
  check(bspl4, RelationalOperator::EQUAL, bspl4_ans, "B-spline computation for the cubic function "
        "with coordinate 0.5 did not yield the expected results.");
  const std::vector<double> bspl4_asm = bSpline<double>(0.38, 4);
  const std::vector<double> bspl4_asm_ans = { 0.039721333, 0.549702667, 0.401430667, 0.009145333 };
  check(bspl4_asm, RelationalOperator::EQUAL, bspl4_asm_ans, "B-spline computation for the cubic "
        "function with coordinate 0.38 did not yield the expected results.");
  const std::vector<double> bspl5 = bSpline<double>(0.137, 5);
  const std::vector<double> bspl5_ans = { 0.023111703, 0.386368047, 0.520943476, 0.069562096,
                                          0.000014678 };
  check(bspl5, RelationalOperator::EQUAL, bspl5_ans, "B-spline computation for the quartic "
        "function with coordinate 0.137 did not yield the expected results.");
  const std::vector<double> dbspl5 = dBSpline<double>(0.137, 5);
  const std::vector<double> dbspl5_ans = { -0.107122608, -0.542060735, 0.405917853, 0.242836931,
                                            0.000428559 };
  check(dbspl5, RelationalOperator::EQUAL, dbspl5_ans, "B-spline derivative computation for the "
        "quartic function with coordinate 0.137 did not yield the expected results.");
  const std::vector<double> xb_crd = {  0.54,  1.50,  2.09,  -0.51, -10.88 };
  const std::vector<double> yb_crd = {  0.47,  1.31,  2.14,  -6.56, -10.05 };
  const std::vector<double> zb_crd = {  8.11,  9.22, 10.58, -10.73,   4.43 };
  const std::vector<double> b_umat = {  0.1, 0.0, 0.0, 0.0,  0.1, 0.0, 0.0, 0.0,  0.1 };
  const std::vector<double> b_invu = { 10.0, 0.0, 0.0, 0.0, 10.0, 0.0, 0.0, 0.0, 10.0 };
  const std::vector<double> e_invu = {  1.0, 0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0,  1.0 };
  std::vector<double> xbspl_coef(20), ybspl_coef(20), zbspl_coef(20);
  std::vector<double> dxbspl_coef(20), dybspl_coef(20), dzbspl_coef(20);
  std::vector<int> a_init(5), b_init(5), c_init(5);
  check(bSplineNoUnity(0.02, 4), RelationalOperator::EQUAL, bSpline(0.02, 4), "B-spline "
        "computations that do not apply the smooth partition of unity do not produce similar "
        "answers to B-splines that do when the order is 4.");
  check(bSplineNoUnity(0.71, 5), RelationalOperator::EQUAL, bSpline(0.71, 5), "B-spline "
        "computations that do not apply the smooth partition of unity do not produce similar "
        "answers to B-splines that do when the order is 5.");
  check(bSplineNoUnity(0.63, 8), RelationalOperator::EQUAL, bSpline(0.63, 8), "B-spline "
        "computations that do not apply the smooth partition of unity do not produce similar "
        "answers to B-splines that do when the order is 8.");
  bSpline<double, double>(xb_crd, yb_crd, zb_crd, 4, b_umat.data(), b_invu.data(), 10, 10, 10,
                          e_invu.data(), &xbspl_coef, &ybspl_coef, &zbspl_coef, &a_init, &b_init,
                          &c_init, &dxbspl_coef, &dybspl_coef, &dzbspl_coef);
  std::vector<double> xbspl_ans, ybspl_ans, zbspl_ans;
  const ImagingMethod img_meth = ImagingMethod::PRIMARY_UNIT_CELL;
  std::vector<int> xbins(5), ybins(5), zbins(5);
  for (int i = 0; i < 5; i++) {
    const double reim_x = imageValue(xb_crd[i], 10.0, img_meth);
    const double reim_y = imageValue(yb_crd[i], 10.0, img_meth);
    const double reim_z = imageValue(zb_crd[i], 10.0, img_meth);
    xbins[i] = floor(reim_x);
    ybins[i] = floor(reim_y);
    zbins[i] = floor(reim_z);
    const std::vector<double> tmp_xb = bSpline(reim_x - floor(reim_x), 4);
    const std::vector<double> tmp_yb = bSpline(reim_y - floor(reim_y), 4);
    const std::vector<double> tmp_zb = bSpline(reim_z - floor(reim_z), 4);
    xbspl_ans.insert(xbspl_ans.end(), tmp_xb.begin(), tmp_xb.end());
    ybspl_ans.insert(ybspl_ans.end(), tmp_yb.begin(), tmp_yb.end());
    zbspl_ans.insert(zbspl_ans.end(), tmp_zb.begin(), tmp_zb.end());
  }
  const std::vector<double> dx_ans = { -0.105800000, -0.642600000,  0.602600000,  0.145800000, 
                                       -0.125000000, -0.625000000,  0.625000000,  0.125000000,
                                       -0.414050000, -0.167850000,  0.577850000,  0.004050000, 
                                       -0.130050000, -0.619850000,  0.629850000,  0.120050000, 
                                       -0.387200000, -0.218400000,  0.598400000,  0.007200000 };
  const std::vector<double> dy_ans = { -0.140450000, -0.608650000,  0.638650000,  0.110450000, 
                                       -0.238050000, -0.475850000,  0.665850000,  0.048050000, 
                                       -0.369800000, -0.250600000,  0.610600000,  0.009800000, 
                                       -0.156800000, -0.589600000,  0.649600000,  0.096800000, 
                                       -0.001250000, -0.546250000,  0.096250000,  0.451250000 };
  const std::vector<double> dz_ans = { -0.396050000, -0.201850000,  0.591850000,  0.006050000, 
                                       -0.304200000, -0.367400000,  0.647400000,  0.024200000, 
                                       -0.088200000, -0.655400000,  0.575400000,  0.168200000, 
                                       -0.266450000, -0.430650000,  0.660650000,  0.036450000, 
                                       -0.162450000, -0.582650000,  0.652650000,  0.092450000 };
  check(dxbspl_coef, RelationalOperator::EQUAL, dx_ans, "B-spline derivatives along the unit "
        "cell A axis were not computed as expected.");
  check(dybspl_coef, RelationalOperator::EQUAL, dy_ans, "B-spline derivatives along the unit "
        "cell B axis were not computed as expected.");
  check(dzbspl_coef, RelationalOperator::EQUAL, dz_ans, "B-spline derivatives along the unit "
        "cell C axis were not computed as expected.");
  check(xbspl_coef, RelationalOperator::EQUAL, xbspl_ans, "B-spline coefficients along the unit "
        "cell A axis were not computed properly based on real coordinates in double-precision "
        "mode.");
  check(ybspl_coef, RelationalOperator::EQUAL, ybspl_ans, "B-spline coefficients along the unit "
        "cell B axis were not computed properly based on real coordinates in double-precision "
        "mode.");
  check(zbspl_coef, RelationalOperator::EQUAL, zbspl_ans, "B-spline coefficients along the unit "
        "cell C axis were not computed properly based on real coordinates in double-precision "
        "mode.");
  check(a_init, RelationalOperator::EQUAL, xbins, "Initial seedings of particles along the unit "
        "cell A axis were not computed as expected.");
  check(b_init, RelationalOperator::EQUAL, ybins, "Initial seedings of particles along the unit "
        "cell B axis were not computed as expected.");
  check(c_init, RelationalOperator::EQUAL, zbins, "Initial seedings of particles along the unit "
        "cell C axis were not computed as expected.");

  // Check the methods for backing out angles relative to the abscissa based on two-dimensional
  // coordinates.
  std::vector<double> recovered_angles(200), recovered_angles_f(200), angles_ans(200);
  for (int i = 0; i < 200; i++) {
    const double angle_i = (static_cast<double>(i - 100) + 0.5) * pi / 100.0;
    const double range_i = (0.9 * xrs128p.uniformRandomNumber()) + 0.1;
    const double x_point = range_i * cos(angle_i);
    const double y_point = range_i * sin(angle_i);
    recovered_angles[i] = angleOnAxes(x_point, y_point);
    recovered_angles_f[i] = angleOnAxesf(x_point, y_point);
    angles_ans[i] = angle_i;
  }
  check(recovered_angles, RelationalOperator::EQUAL, angles_ans, "Angles back-calculated from a "
        "series of points arranged about a set of two-dimensional axes do not obtain the original "
        "values of the angles that generated them.");
  check(recovered_angles_f, RelationalOperator::EQUAL, Approx(angles_ans).margin(4.5e-6),
        "Angles back-calculated from a series of points arranged about a set of two-dimensional "
        "axes do not obtain the original values of the angles that generated them.  The "
        "single-precision mode of the function was used in this test.");
  
  // Check the application of a sorted template
  const int nwild = 10;
  std::vector<int2> wildcards(nwild);
  std::vector<double> associated(nwild);
  for (int i = 0; i < nwild; i++) {
    wildcards[i].x = 100.0 * xrs128p.uniformRandomNumber();
    wildcards[i].y = i;
    associated[i] = xrs128p.gaussianRandomNumber();
  }
  const std::vector<int2> wildcards_cpy = wildcards;
  const std::vector<double> associated_cpy = associated;
  std::sort(wildcards.begin(), wildcards.end(), [](int2 a, int2 b) { return a.x < b.x; });
  const std::vector<double> sorted_association = applyAssociatedSort(associated, wildcards);
  std::vector<size_t> extracted_order(nwild);
  for (int i = 0; i < nwild; i++) {
    extracted_order[i] = wildcards[i].y;
  }
  const std::vector<double> sorted_assoc_ii = applyAssociatedSort(associated, extracted_order);
  bool correspondence = true;
  for (int i = 0; i < nwild; i++) {
    bool match_found = false;
    for (int j = 0; j < nwild; j++) {
      match_found = (match_found ||
                     (wildcards[i].x == wildcards_cpy[j].x &&
                      fabs(sorted_association[i] - associated_cpy[j]) < stormm::constants::tiny));
    }
    correspondence = (correspondence && match_found);
  }
  check(correspondence, "Sorting an array based on an associated, tagged sort did not work as "
        "expected.");
  check(sorted_association, RelationalOperator::EQUAL, sorted_assoc_ii, "Ordering an array based "
        "on the extracted ordering from an assoicated sort does not produce the expected result.");

  // Check the partitioning of some numbers
  const std::vector<int> preferred_lengths = { 15, 13, 11, 9, 7, 5 };
  const std::vector<int> discouraged_lengths = { 12, 10, 8, 6, 4, 3 };
  std::vector<int> segments = partition(23, 16, preferred_lengths, discouraged_lengths);
  const std::vector<int> segments_ans_a = { 15, 7, 1 };
  check(segments, RelationalOperator::EQUAL, segments_ans_a, "The partitioning of 23 did not "
        "return the expected result.");
  segments = partition(27, 11, preferred_lengths, discouraged_lengths);
  const std::vector<int> segments_ans_b = { 9, 9, 9 };
  check(segments, RelationalOperator::EQUAL, segments_ans_b, "The partitioning of 27 did not "
        "return the expected result.");
  segments = partition(23, 23, preferred_lengths, discouraged_lengths);
  const std::vector<int> segments_ans_c = { 23 };
  check(segments, RelationalOperator::EQUAL, segments_ans_c, "The partitioning of 23 with a "
        "maximum value of 23 did not return the expected result.");
  segments = partition(29, 14, preferred_lengths, discouraged_lengths);
  const std::vector<int> segments_ans_d = { 13, 13, 3 };
  check(segments, RelationalOperator::EQUAL, segments_ans_d, "The partitioning of 29 did not "
        "return the expected result.");
  
  // Print results
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
