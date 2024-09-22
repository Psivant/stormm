// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void linearCoreElectrostatics(const Tcalc r, const Tcalc clash_distance, const Tcalc qiqj,
                              Tcalc *ele_contrib, Tcalc *fmag) {
  const Tcalc value_one = 1.0;
  if (r < clash_distance) {
    const Tcalc aparm = -value_one / (clash_distance * clash_distance);
    const Tcalc bparm = (aparm + value_one) / clash_distance;
    *ele_contrib += qiqj * ((aparm * r) + bparm);
    if (fmag != nullptr && r > static_cast<Tcalc>(constants::tiny)) {
      *fmag += qiqj * aparm / r;
    }
  }
  else {
    const Tcalc invr = value_one / r;
    const Tcalc invr2 = invr * invr;
    *ele_contrib += qiqj * invr;
    if (fmag != nullptr) {
      *fmag += -(qiqj * invr * invr2);
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void quadraticCoreElectrostatics(const Tcalc r, const Tcalc clash_distance, const Tcalc qiqj,
                                 Tcalc *ele_contrib, Tcalc *fmag) {
  const Tcalc value_one = 1.0;
  if (r < clash_distance) {
    const Tcalc value_half = 0.5;
    const Tcalc value_two  = 2.0;
    const Tcalc aparm = -value_half /
                        (clash_distance * clash_distance * (clash_distance + value_one));
    const Tcalc bparm = (value_one / clash_distance) -
                        (aparm * (clash_distance + value_one) * (clash_distance + value_one));
    *ele_contrib += qiqj * ((aparm * (r + value_one) * (r + value_one)) + bparm);
    if (fmag != nullptr && r > static_cast<Tcalc>(constants::tiny)) {
      *fmag += qiqj * ((value_two * aparm) + (value_two * aparm / r));
    }
  }
  else {
    const Tcalc invr = value_one / r;
    const Tcalc invr2 = invr * invr;
    *ele_contrib += qiqj * invr;
    if (fmag != nullptr) {
      *fmag += -(qiqj * invr * invr2);
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void cubicCoreElectrostatics(const Tcalc r, const Tcalc clash_distance, const Tcalc qiqj,
                             Tcalc *ele_contrib, Tcalc *fmag) {
  const Tcalc value_one = 1.0;
  if (r < clash_distance) {
    const Tcalc sq_clash_distance = clash_distance * clash_distance;
    const Tcalc value_four = 4.0;
    const Tcalc aparm = -value_one / (sq_clash_distance * sq_clash_distance);
    const Tcalc bparm = value_four / (sq_clash_distance * clash_distance);
    const Tcalc cparm = static_cast<Tcalc>(-6.0) / sq_clash_distance;
    const Tcalc dparm = value_four / clash_distance;
    *ele_contrib += qiqj * ((((((aparm * r) + bparm) * r) + cparm) * r) + dparm);
    if (fmag != nullptr && r > static_cast<Tcalc>(constants::tiny)) {
      *fmag += qiqj * ((((static_cast<Tcalc>(3.0) * aparm * r) +
                         (static_cast<Tcalc>(2.0) * bparm)) * r) + cparm) / r;
    }
  }
  else {
    const Tcalc invr = value_one / r;
    const Tcalc invr2 = invr * invr;
    *ele_contrib += qiqj * invr;
    if (fmag != nullptr) {
      *fmag += -(qiqj * invr * invr2);
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void cubicCoreLennardJones(const Tcalc r, const Tcalc clash_ratio, const Tcalc lja,
                           const Tcalc ljb, const Tcalc ij_sigma, Tcalc *vdw_contrib,
                           Tcalc *fmag) {
  const Tcalc value_one = 1.0;
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const Tcalc rlimit = clash_ratio * ij_sigma;
  const Tcalc invr  = value_one / r;
  const Tcalc invr2 = invr * invr;
  const Tcalc invr4 = invr2 * invr2;
  const Tcalc invr8 = invr4 * invr4;
  if (r < rlimit) {
    const Tcalc invr6  = invr4 * invr2;
    const Tcalc invr7  = invr6 * invr;
    const Tcalc invr8  = invr4 * invr4;
    const Tcalc invr9  = invr8 * invr;
    const Tcalc invr12 = invr6 * invr6;
    const Tcalc invr13 = invr7 * invr6;
    const Tcalc invr14 = invr8 * invr6;
    const Tcalc invr15 = invr9 * invr6;
    Tcalc aparm, bparm, cparm, dparm;
    if (tcalc_is_double) {
      aparm = ( -364.0 * lja * invr15) + ( 56.0 * ljb * invr9);
      bparm = ( 2028.0 * lja * invr14) - (378.0 * ljb * invr8);
      cparm = (-2976.0 * lja * invr13) + (594.0 * ljb * invr7);
      dparm = ( 1313.0 * lja * invr12) - (273.0 * ljb * invr6);
    }
    else {
      aparm = ( -364.0f * lja * invr15) + ( 56.0f * ljb * invr9);
      bparm = ( 2028.0f * lja * invr14) - (378.0f * ljb * invr8);
      cparm = (-2976.0f * lja * invr13) + (594.0f * ljb * invr7);
      dparm = ( 1313.0f * lja * invr12) - (273.0f * ljb * invr6);
    }
    *vdw_contrib += (((((aparm * r) + bparm) * r) + cparm) * r) + dparm;
    if (fmag != nullptr) {
      if (tcalc_is_double) {
        *fmag += ((((3.0 * aparm * r) + (2.0 * bparm)) * r) + cparm) * invr;
      }
      else {
        *fmag += ((((3.0f * aparm * r) + (2.0f * bparm)) * r) + cparm) * invr;
      }
    }
  }
  else {
    const Tcalc invr2 = value_one / (r * r);
    *vdw_contrib += (lja * invr4 * invr8) - (ljb * invr4 * invr2);
    if (fmag != nullptr) {
      if (tcalc_is_double) {
        *fmag += ((6.0 * ljb) - (12.0 * lja * invr4 * invr2)) * invr8;
      }
      else {
        *fmag += ((6.0f * ljb) - (12.0f * lja * invr4 * invr2)) * invr8;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void quarticCoreLennardJones(const Tcalc r, const Tcalc	clash_ratio, const Tcalc lja,
                             const Tcalc ljb, const Tcalc ij_sigma, Tcalc *vdw_contrib,
                             Tcalc *fmag) {
  const Tcalc value_one  = 1.0;
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const Tcalc rlimit = clash_ratio * ij_sigma;
  if (r < rlimit) {
    const Tcalc invrlim = value_one / rlimit;
    const Tcalc invrlim2 = invrlim * invrlim;
    const Tcalc invrlim6 = invrlim2 * invrlim2 * invrlim2;
    Tcalc aparm;
    if (tcalc_is_double) {
      aparm = (((6.0 * ljb) - (12.0 * lja * invrlim6)) * invrlim * invrlim6) /
              ((((((4.0 * rlimit) + 12.0) * rlimit) + 12.0) * rlimit) + 4.0);
    }
    else {
      aparm = (((6.0f * ljb) - (12.0f * lja * invrlim6)) * invrlim * invrlim6) /
              ((((((4.0f * rlimit) + 12.0f) * rlimit) + 12.0f) * rlimit) + 4.0f);
    }
    const Tcalc rlimit_plus_one = rlimit + value_one;
    const Tcalc arlimit_po_four = aparm * (rlimit_plus_one * rlimit_plus_one *
                                           rlimit_plus_one * rlimit_plus_one);
    const Tcalc bparm = (((lja * invrlim6) - ljb) * invrlim6) - (arlimit_po_four);
    const Tcalc r_plus_one = r + value_one;
    const Tcalc arpo_three = aparm * r_plus_one * r_plus_one * r_plus_one;
    *vdw_contrib += (arpo_three * r_plus_one) + bparm;
    if (fmag != nullptr && r > static_cast<Tcalc>(constants::tiny)) {
      *fmag += (tcalc_is_double) ? (4.0 * arpo_three / r) : (4.0f * arpo_three / r);
    }
  }
  else {
    const Tcalc invr2 = value_one / (r * r);
    const Tcalc invr4 = invr2 * invr2;
    *vdw_contrib += (lja * invr4 * invr4 * invr4) - (ljb * invr4 * invr2);
    if (fmag != nullptr) {
      if (tcalc_is_double) {
        *fmag += ((6.0 * ljb) - (12.0 * lja * invr4 * invr2)) * invr4 * invr4;
      }
      else {
        *fmag += ((6.0f * ljb) - (12.0f * lja * invr4 * invr2)) * invr4 * invr4;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4>
void cubicSoftCore(Tcalc4* abcd_coefficients, const double rswitch, const double f_rswitch,
                   const double df_rswitch, const double target_zero) {

  // Compute powers of the interval for the softcore potential
  const double rswitch2 = rswitch  * rswitch;
  const double rswitch3 = rswitch2 * rswitch;

  // Solve the coefficients for each spline.  Use double-precision to solve the coefficients for
  // numerical stability, then convert to the chosen type after all coefficients have been settled.
  std::vector<double> amat(16, 0.0), bvec(4), xvec(4);

  // The value of the spline at the right hand limit is the value of the function or next
  // higher spline at that limit.
  amat[ 0] = rswitch3;
  amat[ 4] = rswitch2;
  amat[ 8] = rswitch;
  amat[12] = 1.0;
  bvec[0] = f_rswitch;

  // The first derivative of the spline at the right hand limit is the first derivative of the
  // function or next higher spline at that limit.
  amat[ 1] = 3.0 * rswitch2;
  amat[ 5] = 2.0 * rswitch;
  amat[ 9] = 1.0;
  bvec[1] = df_rswitch;

  // The slope of the spline at the left hand side limit is set to the requested value.
  amat[10] = 1.0;
  bvec[2] = target_zero;
      
  // The second derivative of the spline at the left hand side will be set to zero.
  amat[7] = 2.0;
  bvec[3] = 0.0;

  // Solve the linear least-squares equation using the faster QR decomposition
  qrSolver(&amat, &xvec, &bvec);

  // Transfer the results to the output
  abcd_coefficients->x = xvec[0];
  abcd_coefficients->y = xvec[1];
  abcd_coefficients->z = xvec[2];
  abcd_coefficients->w = xvec[3];
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc>
void quarticSoftCore(Tcalc4* abcd_coefficients, Tcalc* e_coefficient, const double rswitch,
                     const double f_rswitch, const double df_rswitch, const double d2f_rswitch,
                     const double target_zero) {

  // Compute powers of the interval for the softcore potential
  const double rswitch2 = rswitch  * rswitch;
  const double rswitch3 = rswitch2 * rswitch;
  const double rswitch4 = rswitch2 * rswitch2;
  
  // Solve the coefficients for each spline.  Use double-precision to solve the coefficients for
  // numerical stability, then convert to the chosen type after all coefficients have been settled.
  std::vector<double> amat(25, 0.0), bvec(5), xvec(5);

  // The value of the spline at the right hand limit is the value of the function or next
  // higher spline at that limit.
  amat[ 0] = rswitch4;
  amat[ 5] = rswitch3;
  amat[10] = rswitch2;
  amat[15] = rswitch;
  amat[20] = 1.0;
  bvec[0] = f_rswitch;

  // The first derivative of the spline at the right hand limit is the first derivative of the
  // function or next higher spline at that limit.
  amat[ 1] = 4.0 * rswitch3;
  amat[ 6] = 3.0 * rswitch2;
  amat[11] = 2.0 * rswitch;
  amat[16] = 1.0;
  bvec[1] = df_rswitch;

  // The second derivative of the spline at the right hand limit is the second derivative of
  // the function or next higher spline at that limit.
  amat[ 2] = 12.0 * rswitch2;
  amat[ 7] =  6.0 * rswitch;
  amat[12] =  2.0;
  bvec[2] = d2f_rswitch;

  // The slope of the spline at the left hand side limit is set to the requested value.
  amat[18] = 1.0;
  bvec[3] = target_zero;
      
  // The second derivative of the spline at the left hand side will be set to zero.
  amat[14] = 2.0;
  bvec[4] = 0.0;

  // Solve the linear least-squares equation using the faster QR decomposition
  qrSolver(&amat, &xvec, &bvec);

  // Transfer the results to the output
  abcd_coefficients->x = xvec[0];
  abcd_coefficients->y = xvec[1];
  abcd_coefficients->z = xvec[2];
  abcd_coefficients->w = xvec[3];
  *e_coefficient       = xvec[4];
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2>
void quinticSoftCore(Tcalc4* abcd_coefficients, Tcalc2* ef_coefficients, const double rswitch,
                     const double f_rswitch, const double df_rswitch, const double d2f_rswitch,
                     const double d3f_rswitch, const double target_zero) {

  // Compute powers of the interval for the softcore potential
  const double rswitch2 = rswitch  * rswitch;
  const double rswitch3 = rswitch2 * rswitch;
  const double rswitch4 = rswitch2 * rswitch2;
  const double rswitch5 = rswitch3 * rswitch2;
  
  // Solve the coefficients for each spline.  Use double-precision to solve the coefficients for
  // numerical stability, then convert to the chosen type after all coefficients have been settled.
  std::vector<double> amat(36, 0.0), bvec(6), xvec(6);

  // The value of the spline at the right hand limit is the value of the function or next
  // higher spline at that limit.
  amat[ 0] = rswitch5;
  amat[ 6] = rswitch4;
  amat[12] = rswitch3;
  amat[18] = rswitch2;
  amat[24] = rswitch;
  amat[30] = 1.0;
  bvec[0] = f_rswitch;

  // The first derivative of the spline at the right hand limit is the first derivative of the
  // function or next higher spline at that limit.
  amat[ 1] = 5.0 * rswitch4;
  amat[ 7] = 4.0 * rswitch3;
  amat[13] = 3.0 * rswitch2;
  amat[19] = 2.0 * rswitch;
  amat[25] = 1.0;
  bvec[1] = df_rswitch;

  // The second derivative of the spline at the right hand limit is the second derivative of
  // the function or next higher spline at that limit.
  amat[ 2] = 20.0 * rswitch3;
  amat[ 8] = 12.0 * rswitch2;
  amat[14] =  6.0 * rswitch;
  amat[20] =  2.0;
  bvec[2] = d2f_rswitch;

  // The third derivative of the spline at the right hand limit is the third derivative of the
  // function or next higher spline at that limit.
  amat[ 3] = 60.0 * rswitch2;
  amat[ 9] = 24.0 * rswitch;
  amat[15] =  6.0;
  bvec[3] = d3f_rswitch;

  // The slope of the spline at the left hand side limit is set to the requested value.
  amat[28] = 1.0;
  bvec[4] = target_zero;
      
  // The second derivative of the spline at the left hand side will be set to zero.
  amat[23] = 2.0;
  bvec[5] = 0.0;

  // Solve the linear least-squares equation using QR decomposition
  qrSolver(&amat, &xvec, &bvec);

  // Transfer the results to the prepared arrays
  abcd_coefficients->x = xvec[0];
  abcd_coefficients->y = xvec[1];
  abcd_coefficients->z = xvec[2];
  abcd_coefficients->w = xvec[3];
  ef_coefficients->x   = xvec[4];
  ef_coefficients->y   = xvec[5];
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4>
void exponentialSoftCore(Tcalc4 *abcd_coefficients, const double v, const double rswitch,
                         const double f_rswitch, const double df_rswitch, const double d2f_rswitch,
                         const double d3f_rswitch) {
  std::vector<double> amat(16), inv_amat(16);
  amat[ 0] = 1.0;
  amat[ 4] = 1.0;
  amat[ 8] = 1.0;
  amat[12] = 1.0;
  amat[ 1] =  2.0 * v;
  amat[ 5] =  3.0 * v;
  amat[ 9] =  4.0 * v;
  amat[13] = 0.0;
  amat[ 2] =  4.0 * v * v;
  amat[ 6] =  9.0 * v * v;
  amat[10] = 16.0 * v * v;
  amat[14] = 0.0;

  // Use pow to get the cube to higher precision than straight multiplication
  const double v3 = pow(v, 3.0);
  amat[ 3] =  4.0 * v3;
  amat[ 7] =  9.0 * v3;
  amat[11] = 16.0 * v3;
  amat[15] = 0.0;

  // Form the result vector and solve the linear system of equations
  const std::vector<double> bvec = { f_rswitch, df_rswitch, d2f_rswitch, d3f_rswitch };
  std::vector<double> xvec(4, 0.0);
  invertSquareMatrix(amat.data(), inv_amat.data(), 4);
  matrixVectorMultiply(inv_amat, bvec, &xvec, 4, 4, 1.0, 1.0, 0.0);
  abcd_coefficients->x = xvec[0];
  abcd_coefficients->y = xvec[1];
  abcd_coefficients->z = xvec[2];
  abcd_coefficients->w = xvec[3];
}

} // namespace energy
} // namespace stormm
