#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Reporting/error_format.h"
#include "Math/bspline.h"
#include "Math/vector_ops.h"
#include "pme_util.h"

namespace stormm {
namespace energy {

using constants::ExceptionResponse;
using parse::NumberFormat;
using parse::realToString;
using stmath::bSpline;
using stmath::findBin;

//-------------------------------------------------------------------------------------------------
double ewaldCoefficient(const double cutoff, const double direct_sum_tol) {
  if (cutoff < constants::tiny) {
    rtErr("A cutoff of " + realToString(cutoff, 9, 4, NumberFormat::STANDARD_REAL) +
          " is invalid.", "ewaldCoefficient");
  }
  const double effective_dsum_tol = (direct_sum_tol < 0.0) ? 0.0 : direct_sum_tol;
  double result = 1.0;
  double best_error = fabs((erfc(result * cutoff) / cutoff) - effective_dsum_tol);
  double amin = 0.0;
  double amax = 10.0;
  for (int i = 0; i < 12; i++) {
    const double delta_a = 0.005 * (amax - amin);
    for (int j = 0; j < 200; j++) {
      const double a_update = amin + (delta_a * static_cast<double>(j));
      const double test_error = fabs((erfc(a_update * cutoff) / cutoff) - effective_dsum_tol);
      if (test_error < best_error) {
        best_error = test_error;
        result = a_update;
      }
    }
    amin = result - (5.0 * delta_a);
    amax = result + (5.0 * delta_a);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double pmeGaussianSpread(const double cutoff, const double direct_sum_tol) {
  return 0.5 / ewaldCoefficient(cutoff, direct_sum_tol);
}

//-------------------------------------------------------------------------------------------------
double recoverDirectSumTolerance(const double cutoff, const double ewald_coefficient) {

  // Store a table of the Gaussian spreads for various direct sum tolerances and a cutoff of 9.0
  // Angstroms.  Begin by assuming that the answer probably lies in the range 10^(-1.5) to
  // 10^(-13), then iterating down to a solution.
  static const double known_grms[24] = { 0.890688120, 0.911073483, 0.932906254, 0.956364761,
                                         0.981659487, 1.009040902, 1.038809725, 1.071330613,
                                         1.107050649, 1.146524810, 1.190451675, 1.239724668,
                                         1.295507403, 1.359347738, 1.433356308, 1.520497200,
                                         1.625084093, 1.753677124, 1.916824487, 2.132770200,
                                         2.436381747, 2.904751468, 3.753668530, 5.947450305 };
  const double inp_grms = 0.5 / ewald_coefficient;
  const int bin_idx = findBin(known_grms, inp_grms, 23, ExceptionResponse::SILENT);
  double dtol_low, dtol_hgh, grms_low, grms_hgh;
  if (bin_idx >= 24) {
    rtErr("The direct sum tolerance for this arrangement is far too loose (too large) to be "
          "useful.  Cutoff: " + realToString(cutoff, 9, 4, NumberFormat::STANDARD_REAL) +
          ", gaussian RMS width: " + realToString(inp_grms, 9, 4, NumberFormat::STANDARD_REAL) +
          ".", "recoverDirectSumTolerance");
  }
  else if (bin_idx < 0) {
    grms_low = 0.0;
    grms_hgh = known_grms[0];
    dtol_low = 1.0e-18;
    dtol_hgh = 1.0e-13;
  }
  else {
    grms_low = known_grms[bin_idx];
    grms_hgh = known_grms[bin_idx + 1];
    dtol_low = pow(10.0, -0.5 * static_cast<double>(26 - bin_idx));
    dtol_hgh = pow(10.0, -0.5 * static_cast<double>(26 - (bin_idx + 1)));
  }

  // Estimate the point along the Gaussian spread and scale appropriately
  double dg = (inp_grms - grms_low) / (grms_hgh - grms_low);
  double dtol_est = exp(log(dtol_low) + (dg * log(dtol_hgh / dtol_low)));
  dg *= dtol_est;
  bool converged = false;
  std::vector<bool> overshot;
  overshot.reserve(128);
  int iter = 0;
  do {
    const double grms_est = 0.5 / ewaldCoefficient(cutoff, dtol_est);
    if (grms_est < inp_grms) {

      // The Gaussian spread is too small, meaing that the estimated tolerance is too tight.
      // Loosen (increase) the direct sum tolerance.
      dtol_est += 0.25 * dg;
      overshot.push_back(true);
    }
    else if (grms_est > inp_grms) {

      // The Gaussian spread is too large, meaning that the estimated tolerance is too loose.
      // Tighten (decrease) the direct sum tolerance.
      if (dg < dtol_est) {
        dtol_est -= 0.25 * dg;
      }
      else {
        dtol_est *= 0.5;
      }
      overshot.push_back(false);
    }
    else {
      converged = true;
    }
    if (iter > 0) {
      if (overshot[iter] == overshot[iter - 1]) {
        dg *= 1.25;
      }
      else {
        dg *= 0.5;
      }
    }
    if (converged == false) {
      converged = (dg < 1.0e-15);
    }
    iter++;
  } while (converged == false);
  return dtol_est;
}

//-------------------------------------------------------------------------------------------------
double pmeGammaSum(const int m, const int mesh_length, const int ordr) {
  double result = 1.0;
  if (m != 0) {
    const double dm = m;
    const double dlen = mesh_length;
    const double x = symbols::pi * dm / dlen;
    const double d_ordr = ordr;
    for (int i = 1; i <= 50; i++) {
      const double pi_i = symbols::pi * static_cast<double>(i);
      result += pow(x / (x + pi_i), d_ordr) + pow(x / (x - pi_i), d_ordr);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> pmeLoadBPrefactor(int ordr, int mesh_length) {
  std::vector<double> result(mesh_length, 0.0);

  // Load coefficients for an on-point B-Spline.  Most elements of the array will be zeroes.
  std::vector<double> bspln_arr(mesh_length, 0.0);
  bSpline<double>(1.0, ordr, bspln_arr.data(), nullptr);

  // Compute the moduli of the discrete Fourier transform (DFT).
  std::vector<double> bspln_mod(mesh_length, 0.0);
  const double dlength = mesh_length;
  for (int i = 0; i < mesh_length; i++) {
    double sum_cos = 0.0;
    double sum_sin = 0.0;
    const double difac = twopi * static_cast<double>(i) / dlength;
    for (int j = 0; j < mesh_length; j++) {
      const double arg = difac * static_cast<double>(j);
      sum_cos += bspln_arr[j] * cos(arg);
      sum_sin += bspln_arr[j] * sin(arg);
    }
    bspln_mod[i] = (sum_cos * sum_cos) + (sum_sin * sum_sin);
  }

  // Fix the case where Euler exponential spline interpolation fails.
  for (int i = 0; i < mesh_length; i++) {
    if (bspln_mod[i] < constants::small) {
      if (i > 0 && i < mesh_length) {
        bspln_mod[i] = 0.5 * (bspln_mod[i - 1] + bspln_mod[i + 1]);
      }
      else {

        // This situation should be impossible, but if some numerics have been broken then set
        // the modified B-Spline array to a very high value so that the inverse will have a very
        // small value, effectively zero.
        bspln_mod[i] = 1.0e15;
        rtWarn("DFT modulus at index " + std::to_string(i) + " out of " +
               std::to_string(mesh_length) + " is near zero.", "pmeLoadBPrefactor");
      }
    }
  }

  // Optimize the lambda coefficients
  const int half_length = mesh_length / 2;
  const int twice_ordr = ordr * 2;
  for (int i = 0; i < mesh_length; i++) {
    const int j = (i > half_length) ? i - mesh_length : i;
    const double gsum  = pmeGammaSum(j, mesh_length, ordr);
    const double gsum2 = pmeGammaSum(j, mesh_length, twice_ordr);
    const double lambda = gsum / gsum2;
    result[i] = lambda * lambda / bspln_mod[i];
  }
  
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> pmeLoadMVec(const int mesh_length) {
  std::vector<double> result(mesh_length);
  for (int i = 0; i < mesh_length; i++) {
    result[i] =  (i <= mesh_length / 2) ? i : i - mesh_length;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> pmeLoadMVecShift(const int mesh_length) {
  if (mesh_length <= 0) {
    rtErr("A mesh length of " + std::to_string(mesh_length) + " is invalid.", "pmeLoadMVecShift");
  }
  std::vector<double> result(mesh_length);
  const int half_mesh_length = (mesh_length / 2) + (mesh_length & 0x1);
  for (int i = 0; i < mesh_length; i++) {
    const int idx = (mesh_length - i) % mesh_length;
    result[i] = (idx < (mesh_length >> 1)) ? idx : idx - mesh_length;
  }
  return result;
}

} // namespace energy
} // namespace stormm
