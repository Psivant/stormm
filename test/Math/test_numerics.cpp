#include <climits>
#include "copyright.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/formulas.h"
#include "../../src/Math/juffa.h"
#include "../../src/Math/log_scale_spline.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Math/radial_derivatives.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Potential/pme_util.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/UnitTesting/file_snapshot.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::constants::PrecisionModel;
using stormm::constants::tiny;
using stormm::constants::small;
using stormm::data_types::int95_t;
using stormm::data_types::llint;
#ifndef STORMM_USE_HPC
using stormm::data_types::int2;
using stormm::data_types::double4;
using stormm::data_types::float2;
using stormm::data_types::float4;
#endif
using stormm::data_types::double_type_index;
using stormm::data_types::double4_type_index;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::energy::ewaldCoefficient;
using stormm::energy::recoverDirectSumTolerance;
using stormm::errors::rtWarn;
using stormm::parse::NumberFormat;
using stormm::parse::polyNumericVector;
using stormm::random::Xoshiro256ppGenerator;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using stormm::symbols::pi;
using namespace stormm::stmath;
using namespace stormm::numerics;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Test the split accumulation method over a segment of the number line.
//
// Arguments:
//   llim:        Low limit for sampling
//   hlim:        High limit for sampling
//   incr:        Incrementation value
//   scale_bits:  The number of bits after the decimal
//-------------------------------------------------------------------------------------------------
void testSplitAccumulation(const double llim, const double hlim, const double incr,
                           const int scale_bits, const PrecisionModel lvl) {
  const double scale_factor = pow(2.0, scale_bits);
  const double inv_scale_factor = 1.0 / scale_factor;
  const double inv_overflow_factor = max_llint_accumulation / scale_factor;
  const float scale_factorf = scale_factor;
  int n_basic_fail = 0;
  for (double r = llim; r < hlim; r += incr) {
    const double scaled_r = r * scale_factor;
    switch (lvl) {
    case PrecisionModel::SINGLE:
      {
        const llint dp_result = scaled_r;
        float fr = r;
        const float scaled_fr = fr * scale_factor;
        const llint fp_result = scaled_fr;
        int overflow, workbin;
        hostFloatToInt63(scaled_fr, &workbin, &overflow);
        const llint lloverflow = overflow;
        const llint llworkbin  = workbin;
        const llint fp_reconst = (lloverflow * max_int_accumulation_ll) + llworkbin;
        if (fp_reconst != fp_result) {
          n_basic_fail++;
        }
      }
      break;
    case PrecisionModel::DOUBLE:
      {
        llint workbin;
        int overflow;
        hostDoubleToInt95(r * scale_factor, &workbin, &overflow);
        double dp_reconst = static_cast<double>(workbin) * inv_scale_factor;
        dp_reconst += static_cast<double>(overflow) * inv_overflow_factor;
        if (fabs(dp_reconst - r) > tiny) {
          n_basic_fail++;
        }
      }
      break;
    }
  }
  const double ellim = llim * 0.125;
  const double ehlim = hlim * 0.125;
  const double eincr = incr * 0.125;
  int n_inc_fail = 0;
  for (double r = ellim; r < ehlim; r += eincr) {
    switch (lvl) {
    case PrecisionModel::SINGLE:
      {
        llint dp_result = 0LL;
        llint fp_result = 0LL;
        int overflow = 0;
        int workbin = 0;
        for (int i = 0; i < 9; i++) {
          float fr = r;
          double scaled_r = r * scale_factor;
          float scaled_fr = fr * scale_factorf;
          hostSplitAccumulation(scaled_fr, &workbin, &overflow);
          dp_result += static_cast<llint>(scaled_r);
          fp_result += static_cast<llint>(scaled_fr);
        }
        const llint lloverflow = overflow;
        const llint llworkbin  = workbin;
        const llint fp_reconst = (lloverflow * max_int_accumulation_ll) + llworkbin;
        if (fp_reconst != fp_result) {
          n_inc_fail++;
        }
      }
      break;
    case PrecisionModel::DOUBLE:
      {
        double dp_tracker = 0.0;
        int overflow = 0;
        llint workbin = 0LL;
        for (int i = 0; i < 9; i++) {
          double scaled_r = r * scale_factor;
          hostSplitAccumulation(scaled_r, &workbin, &overflow);
          dp_tracker += r;
        }
        double dp_reconst = static_cast<double>(workbin) * inv_scale_factor;
        dp_reconst += static_cast<double>(overflow) * inv_overflow_factor;
        if (fabs(dp_reconst - dp_tracker) > tiny) {
          n_inc_fail++;
        }
      }
      break;
    }
  }
  check(n_basic_fail == 0, "A total of " + std::to_string(n_basic_fail) + " failures were "
        "recorded converting the range " + realToString(llim, 8, 4) + " : " +
        realToString(hlim, 8, 4) + " to split integer fixed-precision values when sampled with "
        "a " + realToString(incr, 9, 2, NumberFormat::SCIENTIFIC) + " increment.  Precision "
        "model: " + getEnumerationName(lvl) + ".");
  check(n_inc_fail == 0, "A total of " + std::to_string(n_inc_fail) + " failures were recorded "
        "converting the range " + realToString(ellim, 8, 4) + " : " + realToString(ehlim, 8, 4) +
        " to split integer fixed-precision values by incrementing 9x when sampled with a " +
        realToString(incr, 9, 2, NumberFormat::SCIENTIFIC) + " increment.  Precision model: " +
        getEnumerationName(lvl) + ".");
}

//-------------------------------------------------------------------------------------------------
// Test a change in the fixed-precision bit count using the 95-bit accumulators.
//
// Arguments:
//   x:         The real number to represent in fixed precision
//   fp_orig:   The number of bits after the decimal in the first representation
//   fp_new:    The number of bits after the decimal in the second representation
//-------------------------------------------------------------------------------------------------
void testFixedPrecisionChange(const double x, const int fp_orig, const int fp_new) {
  const double oscale = pow(2.0, fp_orig);
  const double nscale = pow(2.0, fp_new);
  const double etol = std::max(pow(2.0, -std::min(fp_new, fp_orig) + 1), 1.0e-12);
  const int95_t ix = hostDoubleToInt95(x * oscale);
  const int95_t in_x = hostChangeFPBits(ix, fp_orig, fp_new);
  const double nx = hostInt95ToDouble(in_x) / nscale;
  check(nx, RelationalOperator::EQUAL, Approx(x).margin(etol), "A fixed-precision bit change "
        "does not work going from " + std::to_string(fp_orig) + " to " + std::to_string(fp_new) +
        " bits in 95-bit format.");

  // Try the 63-bit fixed-precision format
  if (fp_orig < 48 && fp_new < 48) {
    const int2 ix = hostDoubleToInt63(x * oscale);
    const int2 in_x = hostChangeFPBits(ix, fp_orig, fp_new);
    const double nx = hostInt63ToDouble(in_x) / nscale;
    check(nx, RelationalOperator::EQUAL, Approx(x).margin(etol), "A fixed-precision bit change "
          "does not work going from " + std::to_string(fp_orig) + " to " + std::to_string(fp_new) +
          " bits in 63-bit format.");
  }
}

//-------------------------------------------------------------------------------------------------
// Test differentiation in a radially symmetric polynomial function.
//
//
// Arguments:
//   rmax:   The maximum radius in spherical coordinates at which to evaluate the function and its
//           derivatives
//   level:  The level at which to compute finite difference derivatives
//   xrs:    Random number generator used to spread points throughout space
//   acoef:  The A coefficient for the polynomial Ar^3 + Br^2 + Cr + D
//   bcoef:  The B coefficient for the polynomial
//   ccoef:  The C coefficient for the polynomial
//   dcoef:  The D coefficient for the polynomial
//-------------------------------------------------------------------------------------------------
template <typename T>
void testPolynomialDerivatives(const long double rmax, const int level, Xoshiro256ppGenerator *xrs,
                               const long double acoef = 4.7, const long double bcoef = 7.6,
                               const long double ccoef = -1.9, const long double dcoef = -8.2) {

  // Arrays to hold the errors at each derivative level
  const int npts = 64;
  std::vector<double> dx_err(npts), dy_err(npts), dz_err(npts);
  std::vector<double> dxx_err(npts), dxy_err(npts), dxz_err(npts), dyy_err(npts), dzx_err(npts);
  std::vector<double> dxxx_err(npts), dxxy_err(npts), dxxz_err(npts);
  std::vector<double> dxyy_err(npts), dxyz_err(npts), dzxx_err(npts);
  
  // Test the polynomial at 64 points throughout space
  for (int n = 0; n < npts; n++) {

    // Provide the displacements to the highest possible precision
    long double r, xval, yval, zval;
    do {
      r = rmax * xrs->uniformRandomNumber();
      const long double theta = (stormm::symbols::twopi * xrs->uniformRandomNumber()) -
                                stormm::symbols::pi;
      const long double psi   = (stormm::symbols::twopi * xrs->uniformRandomNumber()) -
                                stormm::symbols::pi;
      xval = r * sinl(theta) * cosl(psi);
      yval = r * sinl(theta) * sinl(psi);
      zval = r * cosl(theta);
    } while (fabs(xval) < 0.25 || fabs(yval) < 0.25 || fabs(zval < 0.25));
    const T tr     = r;
    const T tx     = xval;
    const T ty     = yval;
    const T tz     = zval;
    const T tacoef = acoef;
    const T tbcoef = bcoef;
    const T tccoef = ccoef;
    const T tdcoef = dcoef;
    const T value_three = 3.0;
    const T value_two   = 2.0;
    const T value_six   = 6.0;
    const T dfunc   = (((value_three * tacoef * tr) + (value_two * tbcoef)) * tr) + tccoef;
    const T ddfunc  = (value_six * tacoef * tr) + (value_two * tbcoef);
    const T dddfunc = value_six * tacoef;
    const T tdx = radialFirstDerivative<T>(dfunc, tx, tr);
    const T tdy = radialFirstDerivative<T>(dfunc, ty, tr);
    const T tdz = radialFirstDerivative<T>(dfunc, tz, tr);
    int disc_bits;
    switch (level) {
    case 1:
      disc_bits = (sizeof(long double) > 8) ? 22 : 20;
      break;
    case 2:
      disc_bits = (sizeof(long double) > 8) ? 16 : 14;
      break;
    case 3:
      disc_bits = (sizeof(long double) > 8) ? 12 : 10;
      break;
    default:
      rtErr("Invalid derivative level " + std::to_string(level) + ".",
            "testPolynomialDerivatives");
    }
    std::vector<long double> eval_xyz(343);
    const long double xyz_inc = powl(2.0, -disc_bits);
    const long double inv_inc = static_cast<long double>(0.5) / xyz_inc;
    for (int i = 0; i < 7; i++) {
      const long double di = static_cast<long double>(i - 3) * xyz_inc;
      for (int j = 0; j < 7; j++) {
        const long double dj = static_cast<long double>(j - 3) * xyz_inc;
        for (int k = 0; k < 7; k++) {
          const long double dk = static_cast<long double>(k - 3) * xyz_inc;
          const long double xijk = xval + di;
          const long double yijk = yval + dj;
          const long double zijk = zval + dk;
          const long double rijk = sqrtl((xijk * xijk) + (yijk * yijk) + (zijk * zijk));
          eval_xyz[(((k * 7) + j) * 7) + i] = (((((acoef * rijk) + bcoef) * rijk) +
                                                ccoef) * rijk) + dcoef;
        }
      }
    }
    std::vector<long double> eval_dxyz(125), eval_xdyz(125), eval_xydz(125);
    for (int i = 0; i < 5; i++) {
      const int iu = i + 1;
      for (int j = 0; j < 5; j++) {
        const int ju = j + 1;
        for (int k = 0; k < 5; k++) {
          const int ku = k + 1;
          const long double xp = eval_xyz[(((ku * 7) + ju) * 7) + iu + 1];
          const long double xm = eval_xyz[(((ku * 7) + ju) * 7) + iu - 1];
          const long double yp = eval_xyz[(((ku * 7) + ju + 1) * 7) + iu];
          const long double ym = eval_xyz[(((ku * 7) + ju - 1) * 7) + iu];
          const long double zp = eval_xyz[((((ku + 1) * 7) + ju) * 7) + iu];
          const long double zm = eval_xyz[((((ku - 1) * 7) + ju) * 7) + iu];
          eval_dxyz[(((k * 5) + j) * 5) + i] = (xp - xm) * inv_inc;
          eval_xdyz[(((k * 5) + j) * 5) + i] = (yp - ym) * inv_inc;
          eval_xydz[(((k * 5) + j) * 5) + i] = (zp - zm) * inv_inc;
        }
      }
    }
    if (level == 1) {
      dx_err[n] = fabs(eval_dxyz[62] - tdx);
      dy_err[n] = fabs(eval_xdyz[62] - tdy);
      dz_err[n] = fabs(eval_xydz[62] - tdz);
      continue;
    }
    std::vector<long double> eval_ddxyz(27), eval_dxdyz(27), eval_dxydz(27);
    std::vector<long double> eval_xddyz(27), eval_xdydz(27), eval_xyddz(27);
    std::vector<long double> eval_dydxz(27), eval_xdzdy(27), eval_dzydx(27);
    for (int i = 0; i < 3; i++) {
      const int id = i + 1;
      for (int j = 0; j < 3; j++) {
        const int jd = j + 1;
        for (int k = 0; k < 3; k++) {
          const int kd = k + 1;
          const int xp_idx = (((kd * 5) + jd) * 5) + id + 1;
          const int xm_idx = (((kd * 5) + jd) * 5) + id - 1;
          const int yp_idx = (((kd * 5) + jd + 1) * 5) + id;
          const int ym_idx = (((kd * 5) + jd - 1) * 5) + id;
          const int zp_idx = ((((kd + 1) * 5) + jd) * 5) + id;
          const int zm_idx = ((((kd - 1) * 5) + jd) * 5) + id;
          const int ijk_idx = (((k * 3) + j) * 3) + i;
          eval_ddxyz[ijk_idx] = (eval_dxyz[xp_idx] - eval_dxyz[xm_idx]) * inv_inc;
          eval_dxdyz[ijk_idx] = (eval_xdyz[xp_idx] - eval_xdyz[xm_idx]) * inv_inc;
          eval_dxydz[ijk_idx] = (eval_xydz[xp_idx] - eval_xydz[xm_idx]) * inv_inc;
          eval_xddyz[ijk_idx] = (eval_xdyz[yp_idx] - eval_xdyz[ym_idx]) * inv_inc;
          eval_xdydz[ijk_idx] = (eval_xydz[yp_idx] - eval_xydz[ym_idx]) * inv_inc;
          eval_xyddz[ijk_idx] = (eval_xydz[zp_idx] - eval_xydz[zm_idx]) * inv_inc;
          eval_dydxz[ijk_idx] = (eval_dxyz[yp_idx] - eval_dxyz[ym_idx]) * inv_inc;
          eval_xdzdy[ijk_idx] = (eval_xdyz[zp_idx] - eval_xdyz[zm_idx]) * inv_inc;
          eval_dzydx[ijk_idx] = (eval_dxyz[zp_idx] - eval_dxyz[zm_idx]) * inv_inc;
        }
      }
    }
    const T tr2 = tr * tr;
    if (level == 2) {
      dxx_err[n] = fabs(eval_ddxyz[13] - radialSecondDerivative<T>(dfunc, ddfunc, xval, tr));
      dxy_err[n] = fabs(eval_dxdyz[13] - radialSecondDerivative<T>(dfunc, ddfunc, xval, yval, tr,
                                                                   tr2));
      dxz_err[n] = fabs(eval_dxydz[13] - radialSecondDerivative<T>(dfunc, ddfunc, xval, zval, tr,
                                                                   tr2));
      dyy_err[n] = fabs(eval_xddyz[13] - radialSecondDerivative<T>(dfunc, ddfunc, yval, r));
      dzx_err[n] = fabs(eval_dzydx[13] - radialSecondDerivative<T>(dfunc, ddfunc, zval, xval, tr,
                                                                   tr2));
      continue;
    }
    const int xp_shift = 1;
    const int yp_shift = 3;
    const int zp_shift = 9;
    const long double eval_dddxyz = (eval_ddxyz[13 + xp_shift] - eval_ddxyz[13 - xp_shift]) *
                                    inv_inc;
    const long double eval_ddxdyz = (eval_dxdyz[13 + xp_shift] - eval_dxdyz[13 - xp_shift]) *
                                    inv_inc;
    const long double eval_ddxydz = (eval_dxydz[13 + xp_shift] - eval_dxydz[13 - xp_shift]) *
                                    inv_inc;
    const long double eval_dxddyz = (eval_dxdyz[13 + yp_shift] - eval_dxdyz[13 - yp_shift]) *
                                    inv_inc;
    const long double eval_dxdydz = (eval_dxdyz[13 + zp_shift] - eval_dxdyz[13 - zp_shift]) *
                                    inv_inc;
    const long double eval_dzyddx = (eval_ddxyz[13 + zp_shift] - eval_ddxyz[13 - zp_shift]) *
                                    inv_inc;
    dxxx_err[n] = fabs(eval_dddxyz - radialThirdDerivative<T>(dfunc, ddfunc, dddfunc, tx, tr,
                                                              tr2));
    dxxy_err[n] = fabs(eval_ddxdyz - radialThirdDerivative<T>(dfunc, ddfunc, dddfunc, tx, ty, tr,
                                                              tr2));
    dxxz_err[n] = fabs(eval_ddxydz - radialThirdDerivative<T>(dfunc, ddfunc, dddfunc, tx, tz, tr,
                                                              tr2));
    dxyy_err[n] = fabs(eval_dxddyz - radialThirdDerivative<T>(dfunc, ddfunc, dddfunc, ty, tx, tr,
                                                              tr2));
    dxyz_err[n] = fabs(eval_dxdydz -
                       radialThirdDerivative<T>(dfunc, ddfunc, dddfunc, tx, ty, tz, tr, tr2));
    dzxx_err[n] = fabs(eval_dzyddx - radialThirdDerivative<T>(dfunc, ddfunc, dddfunc, tx, tz, tr,
                                                              tr2));
  }

  // Perform the appropriate checks
  const std::vector<double> zero(64, 0.0);
  const size_t ct = std::type_index(typeid(T)).hash_code();
  double etol;
  switch (level) {
  case 1:
    etol = (ct == double_type_index) ? 1.0e-9 : 1.0e-5;
    check(dx_err, RelationalOperator::EQUAL, Approx(zero).margin(etol), "First derivatives of a "
          "radially symmetric polynomial function computed analytically do not agree with those "
          "computed by finite differences.  Direction: X.  Precision model: " +
          getStormmScalarTypeName<T>() + ".");
    check(dy_err, RelationalOperator::EQUAL, Approx(zero).margin(etol), "First derivatives of a "
          "radially symmetric polynomial function computed analytically do not agree with those "
          "computed by finite differences.  Direction: Y.  Precision model: " +
          getStormmScalarTypeName<T>() + ".");
    check(dz_err, RelationalOperator::EQUAL, Approx(zero).margin(etol), "First derivatives of a "
          "radially symmetric polynomial function computed analytically do not agree with those "
          "computed by finite differences.  Direction: Z.  Precision model: " +
          getStormmScalarTypeName<T>() + ".");
    break;
  case 2:
    etol = (ct == double_type_index) ? 1.0e-7 : 1.0e-4;
    check(dxx_err, RelationalOperator::EQUAL, Approx(zero).margin(etol), "Second derivatives of "
          "a radially symmetric polynomial function computed analytically do not agree with those "
          "computed by finite differences.  Direction: X / X.  Precision model: " +
          getStormmScalarTypeName<T>() + ".");
    check(dxy_err, RelationalOperator::EQUAL, Approx(zero).margin(etol), "Second derivatives of "
          "a radially symmetric polynomial function computed analytically do not agree with those "
          "computed by finite differences.  Direction: X / Y.  Precision model: " +
          getStormmScalarTypeName<T>() + ".");
    check(dxz_err, RelationalOperator::EQUAL, Approx(zero).margin(etol), "Second derivatives of "
          "a radially symmetric polynomial function computed analytically do not agree with those "
          "computed by finite differences.  Direction: X / Z.  Precision model: " +
          getStormmScalarTypeName<T>() + ".");
    check(dyy_err, RelationalOperator::EQUAL, Approx(zero).margin(etol), "Second derivatives of "
          "a radially symmetric polynomial function computed analytically do not agree with those "
          "computed by finite differences.  Direction: Y / Y.  Precision model: " +
          getStormmScalarTypeName<T>() + ".");
    check(dzx_err, RelationalOperator::EQUAL, Approx(zero).margin(etol), "Second derivatives of "
          "a radially symmetric polynomial function computed analytically do not agree with those "
          "computed by finite differences.  Direction: Z / X.  Precision model: " +
          getStormmScalarTypeName<T>() + ".");
    break;
  case 3:
    etol = (ct == double_type_index) ? 1.0e-5 : 2.0e-5;
    check(dxxx_err, RelationalOperator::EQUAL, Approx(zero).margin(etol), "Third derivatives of "
          "a radially symmetric polynomial function computed analytically do not agree with those "
          "computed by finite differences.  Direction: X / X / X.  Precision model: " +
          getStormmScalarTypeName<T>() + ".");
    check(dxxy_err, RelationalOperator::EQUAL, Approx(zero).margin(etol), "Third derivatives of "
          "a radially symmetric polynomial function computed analytically do not agree with those "
          "computed by finite differences.  Direction: X / X / Y.  Precision model: " +
          getStormmScalarTypeName<T>() + ".");
    check(dxxz_err, RelationalOperator::EQUAL, Approx(zero).margin(etol), "Third derivatives of "
          "a radially symmetric polynomial function computed analytically do not agree with those "
          "computed by finite differences.  Direction: X / X / Z.  Precision model: " +
          getStormmScalarTypeName<T>() + ".");
    check(dxyy_err, RelationalOperator::EQUAL, Approx(zero).margin(etol), "Third derivatives of "
          "a radially symmetric polynomial function computed analytically do not agree with those "
          "computed by finite differences.  Direction: X / Y / Y.  Precision model: " +
          getStormmScalarTypeName<T>() + ".");
    check(dxyz_err, RelationalOperator::EQUAL, Approx(zero).margin(etol), "Third derivatives of "
          "a radially symmetric polynomial function computed analytically do not agree with those "
          "computed by finite differences.  Direction: X / Y / Z.  Precision model: " +
          getStormmScalarTypeName<T>() + ".");
    check(dzxx_err, RelationalOperator::EQUAL, Approx(zero).margin(etol), "Third derivatives of "
          "a radially symmetric polynomial function computed analytically do not agree with those "
          "computed by finite differences.  Direction: Z / X / X.  Precision model: " +
          getStormmScalarTypeName<T>() + ".");
    break;
  default:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
// Test the radial derivatives using a Coulomb-like interaction.
//
// Arguments:
//   rmax:   The maximum radius in spherical coordinates at which to evaluate the function and its
//           derivatives
//   xrs:    Random number generator used to spread points throughout space
//-------------------------------------------------------------------------------------------------
template <typename T>
void testCoulombDerivatives(const double rmax, Xoshiro256ppGenerator *xrs) {
  const long double xyz_inc = pow(2.0, -15.0);
  const long double inv_inc = 0.5 / xyz_inc;
  const int npts = 64;
  std::vector<double> anlt_dxx(npts), anlt_dxy(npts), webs_dxx(npts), webs_dxy(npts);
  std::vector<double> anlt_dxxx(npts), anlt_dxxy(npts), anlt_dxyy(npts), anlt_dxyz(npts);
  std::vector<double> webs_dxxx(npts), webs_dxxy(npts), webs_dxyy(npts), webs_dxyz(npts);
  for (int i = 0; i < npts; i++) {

    // Provide the displacements to the highest possible precision
    long double r, x, y, z;
    do {
      r = 0.5 + ((rmax - 0.5) * xrs->uniformRandomNumber());
      const long double theta = (stormm::symbols::twopi * xrs->uniformRandomNumber()) -
                                stormm::symbols::pi;
      const long double psi   = (stormm::symbols::twopi * xrs->uniformRandomNumber()) -
                                stormm::symbols::pi;
      x = r * sinl(theta) * cosl(psi);
      y = r * sinl(theta) * sinl(psi);
      z = r * cosl(theta);
    } while (fabs(x) < 0.25 || fabs(y) < 0.25 || fabs(z < 0.25));

    // Compute radial derivatives
    const long double invr = 1.0 / r;
    const long double invr2 = invr * invr;
    long double u[4] = { invr, -invr2, 2.0 * invr2 * invr, -6.0 * invr2 * invr2 };

    // Compute d2/dxx and then d2/dxy
    const long double r2 = r * r;
    anlt_dxx[i] = radialSecondDerivative<T>(u[1], u[2], x, r);
    anlt_dxy[i] = radialSecondDerivative<T>(u[1], u[2], x, y, r, r2);
    webs_dxx[i] = ((3.0 * x * x) - r2) * invr2 * invr2 * invr;
    webs_dxy[i] = 3.0 * x * y * invr2 * invr2 * invr;

    // Compute d3/dxxx, d3/dxxy, d3/dxyy, and d3/dxyz.  Check against their web-based solutions.
    anlt_dxxx[i] = radialThirdDerivative<T>(u[1], u[2], u[3], x, r, r2);
    anlt_dxxy[i] = radialThirdDerivative<T>(u[1], u[2], u[3], x, y, r, r2);
    anlt_dxyy[i] = radialThirdDerivative<T>(u[1], u[2], u[3], y, x, r, r2);
    anlt_dxyz[i] = radialThirdDerivative<T>(u[1], u[2], u[3], x, y, z, r, r2);
    const long double invr7 = invr2 * invr2 * invr2 * invr;
    webs_dxxx[i] = 3.0 * x * ((3.0 * r2) - (5.0 * x * x)) * invr7;
    webs_dxxy[i] = 3.0 * y * (r2 - (5.0 * x * x)) * invr7;
    webs_dxyy[i] = 3.0 * x * (r2 - (5.0 * y * y)) * invr7;
    webs_dxyz[i] = -15.0 * x * y * z * invr7;
  }
  check(anlt_dxx, RelationalOperator::EQUAL, Approx(webs_dxx).margin(stormm::constants::small),
        "The analytic radial function derivative d2/dxx does not agree with a web-based solution "
        "to the same problem.");
  check(anlt_dxy, RelationalOperator::EQUAL, Approx(webs_dxy).margin(stormm::constants::small),
        "The analytic radial function derivative d2/dxy does not agree with a web-based solution "
        "to the same problem.");
  check(anlt_dxxx, RelationalOperator::EQUAL, Approx(webs_dxxx).margin(stormm::constants::small),
        "The analytic radial function derivative d2/dxxx does not agree with a web-based solution "
        "to the same problem.");
  check(anlt_dxxy, RelationalOperator::EQUAL, Approx(webs_dxxy).margin(stormm::constants::small),
        "The analytic radial function derivative d2/dxxy does not agree with a web-based solution "
        "to the same problem.");
  check(anlt_dxyy, RelationalOperator::EQUAL, Approx(webs_dxyy).margin(stormm::constants::small),
        "The analytic radial function derivative d2/dxxy does not agree with a web-based solution "
        "to the same problem.");
  check(anlt_dxyz, RelationalOperator::EQUAL, Approx(webs_dxyz).margin(stormm::constants::small),
        "The analytic radial function derivative d2/dxyz does not agree with a web-based solution "
        "to the same problem.");
}

//-------------------------------------------------------------------------------------------------
// Test the logarithmic spline production.
//
// Arguments:
//   style:            The type of spline to produce and test
//   idx_meth:         Indexing method for the spline table
//   basis_set:        The basis functions to be used to contruct the spline
//   ew_coeff:         The Ewald coefficient governing the splitting factor
//   manissa_bits:     The number of bits of the mantissa with which to build the table index
//   tol:              Tolerance for errors in the spline tables generated under the specified
//                     precision model
//   ulp_opt:          The number of units of least place with which to attempt optimizing the
//                     spline coefficients
//   indexing_offset:  Argument offset for the spline index when using ARG_OFFSET or SQ_ARG_OFFSET
//                     indexing mode
//   xrs:              Random number generator for creating inputs to spline evaluation
//   timer:            Timings object to help track how long tables take to build on CPU or GPU
//-------------------------------------------------------------------------------------------------
template <typename T4>
void testLogSplineRendering(const LogSplineForm style, const TableIndexing idx_meth,
                            const BasisFunctions basis_set, const double ew_coeff,
                            const int mantissa_bits, const double tol, const int ulp_opt,
                            const float indexing_offset, Xoshiro256ppGenerator *xrs,
                            StopWatch *timer) {
  const std::string desc = (std::type_index(typeid(T4)).hash_code() == double4_type_index) ?
                           "dp" : "sp";
  const int spl_timings = timer->addCategory("Compute " + desc + " PME energy by " +
                                             getEnumerationName(idx_meth) + " (" +
                                             std::to_string(mantissa_bits) + " bits, " +
                                             std::to_string(ulp_opt) + " ulp)");
  const int npts = 1024;
  std::vector<double> rpts = uniformRand(xrs, npts, 10.0);
  addScalarToVector(&rpts, 0.9);
  std::vector<double> dbl_eval(npts), spl_eval(npts), flt_eval(npts);
  timer->assignTime(0);
  const double kcoul = stormm::symbols::charmm_gromacs_bioq;
  const float kcoulf = stormm::symbols::charmm_gromacs_bioq;
  const LogScaleSpline<T4> lgsp(style, ew_coeff, kcoul, mantissa_bits, 64.0, 0.015625, idx_meth,
                                basis_set, ulp_opt, indexing_offset);
  timer->assignTime(spl_timings);
  const float ew_coeff_f = ew_coeff;
  const double beta = 2.0 * ew_coeff / sqrt(pi);
  const float betaf = static_cast<float>(2.0 / sqrt(pi)) * ew_coeff_f;
  for (int i = 0; i < npts; i++) {
    const double r = rpts[i];
    const float rf = r;
    switch (style) {
    case LogSplineForm::ELEC_PME_DIRECT:
      dbl_eval[i] = kcoul * erfc(ew_coeff * r) / r;
      flt_eval[i] = kcoulf * erfcf(ew_coeff_f * rf) / rf;
      break;
    case LogSplineForm::ELEC_PME_DIRECT_EXCL:
      dbl_eval[i] = kcoul * (erfc(ew_coeff * r) - 1.0) / r;
      flt_eval[i] = kcoulf * (erfcf(ew_coeff_f * rf) - 1.0f) / rf;
      break;
    case LogSplineForm::DELEC_PME_DIRECT:
      {
        const double ew_coeff_r = ew_coeff * r;
        const float ew_coeff_rf = ew_coeff_f * rf;
        dbl_eval[i] = -kcoul * ((beta * exp(-ew_coeff_r * ew_coeff_r)) +
                                (erfc(ew_coeff_r) / r)) / (r * r);
        flt_eval[i] = -kcoulf * ((betaf * expf(-ew_coeff_rf * ew_coeff_rf)) +
                                 (erfcf(ew_coeff_rf) / rf)) / (rf * rf);
      }
      break;
    case LogSplineForm::DELEC_PME_DIRECT_EXCL:
      {
        const double ew_coeff_r = ew_coeff * r;
        const double invr = 1.0 / r;
        const float ew_coeff_rf = ew_coeff_f * rf;
        const double invrf = 1.0f / rf;
        dbl_eval[i] = -kcoul * ((beta * exp(-ew_coeff_r * ew_coeff_r)) +
                                ((erfc(ew_coeff_r) - 1.0) * invr)) * invr * invr;
        flt_eval[i] = -kcoulf * ((betaf * expf(-ew_coeff_rf * ew_coeff_rf)) +
                                 ((erfcf(ew_coeff_rf) - 1.0f) * invrf)) * invrf * invrf;
      }
      break;
    case LogSplineForm::CUSTOM:
      break;
    }
    spl_eval[i] = lgsp.evaluateByRealArg(r);
  }
  check(spl_eval, RelationalOperator::EQUAL, Approx(dbl_eval).margin(tol), "The spline table "
        "for " + getEnumerationName(style) + ", indexed by " + getEnumerationName(idx_meth) +
        " and constructed with " + getEnumerationName(basis_set) + " (" +
        std::to_string(mantissa_bits) + " bits of the mantissa for indexing), does not meet "
        "expectations for accuracy.");
  timer->assignTime(0);
}

//-------------------------------------------------------------------------------------------------
// Test various implementations from Norbert Juffa's math library for calculating important
// functions in 32-bit arithmetic.
//-------------------------------------------------------------------------------------------------
void testJuffaImplementations() {

  // Check the expf() function
  const int npts = 401;
  std::vector<float> fpts(npts);
  std::vector<double> exp_ref(npts), exp_jff(npts), exp_isc(npts);
  for (int i = 0; i < npts; i++) {
    fpts[i] = static_cast<double>(i - 200) * 0.05;
    exp_ref[i] = exp(fpts[i]);
    exp_jff[i] = expJf(fpts[i]);
    exp_isc[i] = expf(fpts[i]);
  }
  check(exp_jff, RelationalOperator::EQUAL, Approx(exp_ref).margin(1.5e-3), "The erfc() "
        "approximation for a given Ewald coefficient does not match the reference for values in "
        "a range [0, 12.9).");

  // Check the erfcf() function
  JfErfc ewld(0.31);
  std::vector<double> erfc_ref(npts), erfc_jff(npts), erfc_isc(npts);
  for (int i = 0; i < npts; i++) {
    fpts[i] = static_cast<double>(i + 30) * 0.03;
    erfc_ref[i] = erfc(0.31f * fpts[i]);
    erfc_jff[i] = ewld.erfcJf(fpts[i]);
    erfc_isc[i] = erfcf(0.31f * fpts[i]);
  }
  check(erfc_jff, RelationalOperator::EQUAL, Approx(erfc_ref).margin(8.0e-8), "The erfc() "
        "approximation for a given Ewald coefficient does not match the reference for values in "
        "a range [0, 12.9).");
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

  // Section 1
  section("Long long int <=> float2 conversion");

  // Section 2
  section("Real to int2, 63-bit conversion");

  // Section 3
  section("Real to int95_t, 95-bit conversion");

  // Section 4
  section("Overloads of various conversion functions");

  // Section 5
  section("Compute partial derivatives of radially-symmetric functions");

  // Section 6
  section("Compute logarithmically indexed spline tables");

  // Section 7
  section("Test adaptations of Norbert Juffa's 32-bit math functions");

  // Make a series of prime numbers
  std::vector<llint> primes(1, 2);
  int p = 3;
  int g = 3;
  while (g < 10000) {
    const int q = g / p;
    if (g - (p * q) == 0) {
      g += 2;
      p = 3;
      continue;
    }
    if (p > q) {
      primes.push_back(g);
      g += 2;
      p = 3;
    }
    else {
      p += 2;
    }
  }

  // Make a series of integers covering a large portion of the long long int range
  const int n_pts = 16277216;
  const int half_pts = n_pts / 2;
  std::vector<llint> exact_test_range(n_pts);
  llint acc = 1LL;
  llint inc = 2LL;
  int pcon = 0;
  const int n_primes = primes.size();
  for (int i = half_pts; i < n_pts; i++) {
    exact_test_range[i] = acc;
    exact_test_range[n_pts - i - 1] = -acc + 1LL;
    acc += inc;
    inc = (inc > 134217728LL) ? primes[pcon] : inc + primes[pcon];
    pcon++;
    if (pcon == n_primes) {
      pcon = 0;
    }
  }

  // Check the condition of the number series
  section(1);
  const llint etr_expected_min = -545794316914904LL;  
  const llint etr_expected_max =  545794316914905LL;
  const bool exact_passes = (minValue(exact_test_range) == etr_expected_min &&
                             maxValue(exact_test_range) == etr_expected_max);
  check(exact_passes, "A series of long long integers needed for many tests was not constructed "
        "properly.  Its contents should span a range [" + std::to_string(etr_expected_min) + ", " +
        std::to_string(etr_expected_max) + "], but instead it spans " +
        std::to_string(minValue(exact_test_range)) + " to " +
        std::to_string(maxValue(exact_test_range)) + ".  Dependent tests will be skipped.");
  const TestPriority do_exact_tests = (exact_passes) ? TestPriority::CRITICAL :
                                                       TestPriority::ABORT;
  
  // Compute for a slightly larger series of numbers
  std::vector<llint> large_test_range(n_pts);
  acc = 1LL;
  inc = 2LL;
  pcon = 0;
  for (int i = half_pts; i < n_pts; i++) {
    large_test_range[i] = acc;
    large_test_range[n_pts - i - 1] = -acc + 1LL;
    acc += inc;
    inc = (inc > 260000000LL) ? primes[pcon] : inc + primes[pcon];
    pcon++;
    if (pcon == n_primes) {
      pcon = 0;
    }
  }
  const llint ltr_expected_max =  1056755308863635LL;
  const llint ltr_expected_min = -1056755308863634LL;
  const bool large_passes = (minValue(large_test_range) == ltr_expected_min &&
                             maxValue(large_test_range) == ltr_expected_max);
  check(large_passes, "A series of long long integers needed for many tests was not constructed "
        "properly.  Its contents should span a range [" + std::to_string(ltr_expected_min) + ", " +
        std::to_string(ltr_expected_max) + "], but instead it spans " +
        std::to_string(minValue(large_test_range)) + " to " +
        std::to_string(maxValue(large_test_range)) + ".  Dependent tests will be skipped.");
  const TestPriority do_large_tests = (exact_passes) ? TestPriority::CRITICAL :
                                                       TestPriority::ABORT;

  // Test direct conversion
  section(2);
  testSplitAccumulation(-48.5, -47.5, 1.0e-5, 26, PrecisionModel::SINGLE);
  testSplitAccumulation(-32.5, -31.5, 1.0e-5, 26, PrecisionModel::SINGLE);
  testSplitAccumulation( -0.5,   0.5, 1.0e-5, 26, PrecisionModel::SINGLE);
  testSplitAccumulation( 31.5,  32.5, 1.0e-5, 26, PrecisionModel::SINGLE);
  testSplitAccumulation( 47.5,  48.5, 1.0e-5, 26, PrecisionModel::SINGLE);
  testSplitAccumulation(-64.5, -63.5, 1.0e-5, 26, PrecisionModel::SINGLE);
  testSplitAccumulation( 63.5,  64.5, 1.0e-5, 26, PrecisionModel::SINGLE);
  section(3);
  testSplitAccumulation(-48.5, -47.5, 1.0e-5, 58, PrecisionModel::DOUBLE);
  testSplitAccumulation(-32.5, -31.5, 1.0e-5, 58, PrecisionModel::DOUBLE);
  testSplitAccumulation( -0.5,   0.5, 1.0e-5, 58, PrecisionModel::DOUBLE);
  testSplitAccumulation( 31.5,  32.5, 1.0e-5, 58, PrecisionModel::DOUBLE);
  testSplitAccumulation( 47.5,  48.5, 1.0e-5, 58, PrecisionModel::DOUBLE);
  testSplitAccumulation(-64.5, -63.5, 1.0e-5, 58, PrecisionModel::DOUBLE);
  testSplitAccumulation( 63.5,  64.5, 1.0e-5, 58, PrecisionModel::DOUBLE);

  // Test various overloads of the fixed-precision functions
  const int npts = 16;
  std::vector<float> fx(npts), fy(npts), fz(npts), rec_fx(npts), rec_fy(npts), rec_fz(npts);
  std::vector<int> ifx(npts), ify(npts), ifz(npts), ifx_ovrf(npts), ify_ovrf(npts), ifz_ovrf(npts);
  std::vector<double> dx(npts), dy(npts), dz(npts), rec_dx(npts), rec_dy(npts), rec_dz(npts);
  std::vector<llint> idx(npts), idy(npts), idz(npts);
  std::vector<int> idx_ovrf(npts), idy_ovrf(npts), idz_ovrf(npts);
  Xoshiro256ppGenerator xrs;
  const double dscale = pow(2.0, 68.0);
  const float fscale  = pow(2.0, 48.0);
  for (int i = 0; i < npts; i++) {
    fx[i] = (xrs.uniformRandomNumber() - 0.5) * fscale;
    fy[i] = (xrs.uniformRandomNumber() - 0.5) * fscale;
    fz[i] = (xrs.uniformRandomNumber() - 0.5) * fscale;
    dx[i] = (xrs.uniformRandomNumber() - 0.5) * dscale;
    dy[i] = (xrs.uniformRandomNumber() - 0.5) * dscale;
    dz[i] = (xrs.uniformRandomNumber() - 0.5) * dscale;
  }
  hostFloatToInt63(fx, fy, fz, &ifx, &ifx_ovrf, &ify, &ify_ovrf, &ifz, &ifz_ovrf);
  hostDoubleToInt95(dx, dy, dz, &idx, &idx_ovrf, &idy, &idy_ovrf, &idz, &idz_ovrf);
  hostInt63ToFloat(&rec_fx, &rec_fy, &rec_fz, ifx, ifx_ovrf, ify, ify_ovrf, ifz, ifz_ovrf);
  hostInt95ToDouble(&rec_dx, &rec_dy, &rec_dz, idx, idx_ovrf, idy, idy_ovrf, idz, idz_ovrf);
  check(fx, RelationalOperator::EQUAL, rec_fx, "Single-precision floating point numbers are not "
        "recovered when passed through the int63_t (int2) fixed-precision representation (X).");
  check(dx, RelationalOperator::EQUAL, rec_dx, "Double-precision floating point numbers are not "
        "recovered when passed through the int95_t fixed-precision representation (X).");
  check(fy, RelationalOperator::EQUAL, rec_fy, "Single-precision floating point numbers are not "
        "recovered when passed through the int63_t (int2) fixed-precision representation (Y).");
  check(dy, RelationalOperator::EQUAL, rec_dy, "Double-precision floating point numbers are not "
        "recovered when passed through the int95_t fixed-precision representation (Y).");
  check(fz, RelationalOperator::EQUAL, rec_fz, "Single-precision floating point numbers are not "
        "recovered when passed through the int63_t (int2) fixed-precision representation (Z).");
  check(dz, RelationalOperator::EQUAL, rec_dz, "Double-precision floating point numbers are not "
        "recovered when passed through the int95_t fixed-precision representation (Z).");
  Hybrid<float> hfx(npts), rec_hfx(npts);
  Hybrid<double> hdx(npts), rec_hdx(npts);
  Hybrid<int> ihfx(npts), ihfx_ovrf(npts), ihdx_ovrf(npts);
  Hybrid<llint> ihdx(npts);
  for (int i = 0; i < npts; i++) {
    hfx.putHost((xrs.uniformRandomNumber() - 0.5) * fscale, i);
    hdx.putHost((xrs.uniformRandomNumber() - 0.5) * dscale, i);
  }
  hostFloatToInt63(hfx, &ihfx, &ihfx_ovrf);
  hostDoubleToInt95(hdx, &ihdx, &ihdx_ovrf);
  hostInt63ToFloat(&rec_hfx, ihfx, ihfx_ovrf);
  hostInt95ToDouble(&rec_hdx, ihdx, ihdx_ovrf);
  check(hfx.readHost(), RelationalOperator::EQUAL, rec_hfx.readHost(), "Single-precision floating "
        "point numbers are not recovered when passed through the int63_t (int2) fixed-precision "
        "representation (X), when using a Hybrid object.");
  check(hdx.readHost(), RelationalOperator::EQUAL, rec_hdx.readHost(), "Double-precision floating "
        "point numbers are not recovered when passed through the int95_t fixed-precision "
        "representation (X), when using a Hybrid object.");
  const int nscatter_tests = 2048;
  const int nincrement_tests = 256;
  std::vector<llint> large_z_a(nscatter_tests);
  std::vector<llint> large_z_b(nscatter_tests);
  std::vector<llint> large_z_c(nscatter_tests);
  std::vector<llint> large_z_d(nincrement_tests);
  const llint powertwo_a = ipowl(2, 51);
  const llint powertwo_b = ipowl(2, 26);
  const llint powertwo_c = ipowl(2, 20);
  const llint scatter_llim = -nscatter_tests / 4;
  const llint scatter_hlim =  nscatter_tests / 4;
  for (llint i = scatter_llim; i < scatter_hlim; i++) {
    const llint base_a = i * powertwo_a;
    const llint base_b = i * powertwo_b;
    const llint base_c = i * powertwo_c;
    for (int j = 0; j < 2; j++) {
      const int jidx = 2 * (i - scatter_llim) + j;
      large_z_a[jidx] = base_a + static_cast<int>(xrs.uniformRandomNumber() * 1024.0);
      large_z_b[jidx] = base_b + static_cast<int>(xrs.uniformRandomNumber() * 1024.0);
      large_z_c[jidx] = base_c + static_cast<int>(xrs.uniformRandomNumber() * 1024.0);
    }
  }
  const llint incr_llim = -nincrement_tests / 4;
  const llint incr_hlim =  nincrement_tests / 4;
  for (llint i = incr_llim; i < incr_hlim; i++) {
    large_z_d[2 * (i + 64)    ] = static_cast<llint>(INT_MAX) + i;
    large_z_d[2 * (i + 64) + 1] = static_cast<llint>(INT_MIN) + i;
  }
  std::vector<int2> split_short_a(nscatter_tests);
  std::vector<int2> split_short_b(nscatter_tests);
  std::vector<int2> split_short_c(nscatter_tests);
  std::vector<int2> split_short_d(nincrement_tests);
  const std::vector<double> dblref_a(large_z_a.begin(), large_z_a.end());
  const std::vector<double> dblref_b(large_z_b.begin(), large_z_b.end());
  const std::vector<double> dblref_c(large_z_c.begin(), large_z_c.end());
  const std::vector<double> dblref_d(large_z_d.begin(), large_z_d.end());
  std::vector<double> dblshort_a(nscatter_tests), dblshort_b(nscatter_tests);
  std::vector<double> dblshort_c(nscatter_tests), dblshort_d(nincrement_tests);
  for (int i = 0; i < nscatter_tests; i++) {
    split_short_a[i] = hostLongLongToInt63(large_z_a[i]);
    split_short_b[i] = hostLongLongToInt63(large_z_b[i]);
    split_short_c[i] = hostLongLongToInt63(large_z_c[i]);
    dblshort_a[i] = hostInt63ToDouble(split_short_a[i]);
    dblshort_b[i] = hostInt63ToDouble(split_short_b[i]);
    dblshort_c[i] = hostInt63ToDouble(split_short_c[i]);
  }
  for (int i = 0; i < nincrement_tests; i++) {
    split_short_d[i] = hostLongLongToInt63(large_z_d[i]);
    dblshort_d[i] = hostInt63ToDouble(split_short_d[i]);
  }
  check(dblshort_a, RelationalOperator::EQUAL, dblref_a, "Conversion of long long integers to "
        "int2 did not conerve information in the broadest-spectrum scan.");
  check(dblshort_b, RelationalOperator::EQUAL, dblref_b, "Conversion of long long integers to "
        "int2 did not conerve information in a scan that focuses on the seam between the two "
        "32-bit values.");
  check(dblshort_c, RelationalOperator::EQUAL, dblref_c, "Conversion of long long integers to "
        "int2 did not conerve information in a fine-grained scan that never exceeds the primary "
        "value's range.");
  check(dblshort_d, RelationalOperator::EQUAL, dblref_d, "Conversion of long long integers to "
        "int2 did not conerve information in an incremental scan that never exceeds the primary "
        "value's range.");

  // Test functions for changing the precision model
  section(4);
  testFixedPrecisionChange(67.029, 24, 36);
  testFixedPrecisionChange(67.029, 36, 24);
  testFixedPrecisionChange(-23.981, 24, 36);
  testFixedPrecisionChange(-76.329, 63, 70);
  testFixedPrecisionChange(-76.329, 70, 63);

  // Test the numerical precision of analytic partial derivatives for radially symmetric functions
  section(5);
  testPolynomialDerivatives<double>(1.7, 1, &xrs);
  testPolynomialDerivatives<double>(1.7, 2, &xrs);
  testPolynomialDerivatives<double>(1.7, 3, &xrs);
  testPolynomialDerivatives<float>(1.7, 1, &xrs);
  testPolynomialDerivatives<float>(1.7, 2, &xrs);
  testPolynomialDerivatives<float>(1.7, 3, &xrs);

  // Test Coulomb-like functions with radial differentiation
  testCoulombDerivatives<double>(2.6, &xrs);

  // Test logarithmic spline tables and associated computations
  section(6);
  const double ew_loose = ewaldCoefficient(10.0, 1.0e-5);
  const double ew_tight = ewaldCoefficient(12.0, 2.0e-7);
  check(ew_loose, RelationalOperator::EQUAL, Approx(0.2751063906).margin(1.0e-8), "The Ewald "
        "coefficient was not computed as expected for a cutoff of 10.0 Angstroms and direct sum "
        "tolerance 1.0e-5.");
  check(ew_tight, RelationalOperator::EQUAL, Approx(0.2779192472).margin(1.0e-8), "The Ewald "
        "coefficient was not computed as expected for a cutoff of 12.0 Angstroms and direct sum "
        "tolerance 1.0e-7.");
  const double dsum_tol_loose = recoverDirectSumTolerance(10.5, 0.31);
  check(dsum_tol_loose, RelationalOperator::EQUAL, Approx(3.961125e-7).margin(1.0e-10),
        "The inverse method of computing the direct sum tolerance for a given Ewald coefficient "
        "does not produce the expected result.");  
  testLogSplineRendering<double4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::ARG,
                                  BasisFunctions::POLYNOMIAL, ew_loose, 8, 1.0e-9, 0, 0.0, &xrs,
                                  &timer);
  testLogSplineRendering<float4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::ARG,
                                 BasisFunctions::POLYNOMIAL, ew_loose, 5, 2.5e-5, 0, 0.0, &xrs,
                                 &timer);
  testLogSplineRendering<double4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::SQUARED_ARG,
                                  BasisFunctions::POLYNOMIAL, ew_loose, 8, 1.0e-9, 0, 0.0, &xrs,
                                  &timer);
  testLogSplineRendering<float4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::SQUARED_ARG,
                                 BasisFunctions::POLYNOMIAL, ew_loose, 5, 3.5e-5, 0, 0.0, &xrs,
                                 &timer);
  testLogSplineRendering<double4>(LogSplineForm::ELEC_PME_DIRECT_EXCL, TableIndexing::ARG,
                                  BasisFunctions::POLYNOMIAL, ew_loose, 8, 1.0e-9, 0, 0.0, &xrs,
                                  &timer);
  testLogSplineRendering<float4>(LogSplineForm::ELEC_PME_DIRECT_EXCL, TableIndexing::ARG,
                                 BasisFunctions::POLYNOMIAL, ew_loose, 5, 1.1e-5, 0, 0.0, &xrs,
                                 &timer);
  testLogSplineRendering<double4>(LogSplineForm::ELEC_PME_DIRECT_EXCL, TableIndexing::SQUARED_ARG,
                                  BasisFunctions::POLYNOMIAL, ew_loose, 8, 1.0e-9, 0, 0.0, &xrs,
                                  &timer);
  testLogSplineRendering<float4>(LogSplineForm::ELEC_PME_DIRECT_EXCL, TableIndexing::SQUARED_ARG,
                                 BasisFunctions::POLYNOMIAL, ew_loose, 5, 1.5e-5, 0, 0.0, &xrs,
                                 &timer);
  testLogSplineRendering<double4>(LogSplineForm::DELEC_PME_DIRECT, TableIndexing::ARG,
                                  BasisFunctions::POLYNOMIAL, ew_loose, 8, 1.5e-8, 0, 0.0, &xrs,
                                  &timer);
  testLogSplineRendering<float4>(LogSplineForm::DELEC_PME_DIRECT, TableIndexing::ARG,
                                 BasisFunctions::POLYNOMIAL, ew_loose, 5, 9.6e-5, 0, 0.0, &xrs,
                                 &timer);
  testLogSplineRendering<double4>(LogSplineForm::DELEC_PME_DIRECT, TableIndexing::SQUARED_ARG,
                                  BasisFunctions::POLYNOMIAL, ew_loose, 8, 3.0e-9, 0, 0.0, &xrs,
                                  &timer);
  testLogSplineRendering<float4>(LogSplineForm::DELEC_PME_DIRECT, TableIndexing::SQUARED_ARG,
                                 BasisFunctions::POLYNOMIAL, ew_loose, 5, 1.2e-4, 0, 0.0,
                                 &xrs, &timer);
  testLogSplineRendering<double4>(LogSplineForm::DELEC_PME_DIRECT_EXCL, TableIndexing::ARG,
                                  BasisFunctions::POLYNOMIAL, ew_loose, 8, 1.0e-9, 0, 0.0,
                                  &xrs, &timer);
  testLogSplineRendering<float4>(LogSplineForm::DELEC_PME_DIRECT_EXCL, TableIndexing::ARG,
                                 BasisFunctions::POLYNOMIAL, ew_loose, 5, 3.0e-5, 0, 0.0,
                                 &xrs, &timer);
  testLogSplineRendering<double4>(LogSplineForm::DELEC_PME_DIRECT_EXCL, TableIndexing::SQUARED_ARG,
                                  BasisFunctions::POLYNOMIAL, ew_loose, 8, 1.0e-9, 0, 0.0,
                                  &xrs, &timer);
  testLogSplineRendering<float4>(LogSplineForm::DELEC_PME_DIRECT_EXCL, TableIndexing::SQUARED_ARG,
                                 BasisFunctions::POLYNOMIAL, ew_loose, 5, 3.0e-5, 0, 0.0,
                                 &xrs, &timer);
  testLogSplineRendering<double4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::ARG,
                                  BasisFunctions::MIXED_FRACTIONS, ew_loose, 8, 1.0e-9, 0, 0.0,
                                  &xrs, &timer);
  testLogSplineRendering<float4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::ARG,
                                 BasisFunctions::MIXED_FRACTIONS, ew_loose, 5, 5.0e-5, 0, 0.0,
                                 &xrs, &timer);
  testLogSplineRendering<double4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::SQUARED_ARG,
                                  BasisFunctions::MIXED_FRACTIONS, ew_loose, 8, 1.0e-9, 0, 0.0,
                                  &xrs, &timer);
  testLogSplineRendering<float4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::SQUARED_ARG,
                                 BasisFunctions::MIXED_FRACTIONS, ew_loose, 5, 3.5e-5, 0, 0.0,
                                 &xrs, &timer);

  // Introduce the offset argument indexing methods
  testLogSplineRendering<float4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::ARG_OFFSET,
                                 BasisFunctions::MIXED_FRACTIONS, ew_loose, 5, 6.0e-5, 0, 0.015625,
                                 &xrs, &timer);
  testLogSplineRendering<float4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::ARG_OFFSET,
                                 BasisFunctions::MIXED_FRACTIONS, ew_loose, 5, 6.5e-5, 2, 0.015625,
                                 &xrs, &timer);
  testLogSplineRendering<float4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::SQ_ARG_OFFSET,
                                 BasisFunctions::MIXED_FRACTIONS, ew_loose, 5, 6.5e-5, 0, 0.015625,
                                 &xrs, &timer);
  testLogSplineRendering<float4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::SQ_ARG_OFFSET,
                                 BasisFunctions::MIXED_FRACTIONS, ew_loose, 5, 6.5e-5, 2, 0.015625,
                                 &xrs, &timer);
  testLogSplineRendering<float4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::ARG_OFFSET,
                                 BasisFunctions::POLYNOMIAL, ew_loose, 5, 4.0e-5, 0, 0.015625,
                                 &xrs, &timer);
  testLogSplineRendering<float4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::ARG_OFFSET,
                                 BasisFunctions::POLYNOMIAL, ew_loose, 5, 3.8e-5, 3, 0.015625,
                                 &xrs, &timer);
  testLogSplineRendering<float4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::SQ_ARG_OFFSET,
                                 BasisFunctions::POLYNOMIAL, ew_loose, 5, 5.0e-5, 0, 0.015625,
                                 &xrs, &timer);
  testLogSplineRendering<float4>(LogSplineForm::ELEC_PME_DIRECT, TableIndexing::SQ_ARG_OFFSET,
                                 BasisFunctions::POLYNOMIAL, ew_loose, 5, 2.8e-5, 3, 0.015625,
                                 &xrs, &timer);

  // Test Norbert Juffa's various 32-bit floating point functions
  section(7);
  testJuffaImplementations();
  
  // Print results
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
