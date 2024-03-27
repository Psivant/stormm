// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void checkSplineRange(const int spline_idx, const int spline_count, const Tcalc interval,
                      const Tcalc r, const char* caller) {
  if (spline_idx < 0 || spline_idx >= spline_count) {
    rtErr("An inter-particle distance of " + realToString(r, 9, 4, NumberFormat::STANDARD_REAL) +
          " Angstroms is outside the range of a tabulated spline with " +
          std::to_string(spline_count) + " segments and an interval of " +
          realToString(interval, 9, 4, NumberFormat::STANDARD_REAL) + ".", caller);
  }
}
  
//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc>
void evaluateCubicSpline(const Tcalc4 abcd_coefficients, const Tcalc r, Tcalc *u_contrib,
                         Tcalc *fmag) {
  *u_contrib += (((((abcd_coefficients.x * r) + abcd_coefficients.y) * r) +
                  abcd_coefficients.z) * r) + abcd_coefficients.w;
  if (fmag != nullptr) {
    if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
      *fmag += (((3.0 * abcd_coefficients.x * r) + (2.0 * abcd_coefficients.y)) * r) +
               abcd_coefficients.z;
    }
    else {
      *fmag += (((3.0f * abcd_coefficients.x * r) + (2.0f * abcd_coefficients.y)) * r) +
               abcd_coefficients.z;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateCubicSpline(const Tcalc4 abcd_coefficients, const Tcalc r) {
  Tcalc2 result;
  result.x = (((((abcd_coefficients.x * r) + abcd_coefficients.y) * r) +
               abcd_coefficients.z) * r) + abcd_coefficients.w;
  if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
    result.y = (((3.0 * abcd_coefficients.x * r) + (2.0 * abcd_coefficients.y)) * r) +
               abcd_coefficients.z;
  }
  else {
    result.y = (((3.0f * abcd_coefficients.x * r) + (2.0f * abcd_coefficients.y)) * r) +
               abcd_coefficients.z;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc>
void evaluateCubicSpline(const Tcalc4* abcd_coefficients, const int spline_count,
                         const Tcalc interval, const Tcalc r, Tcalc *u_contrib, Tcalc *fmag) {
  const Tcalc r_itvl = r / interval;
  const int spline_idx = r_itvl;
  checkSplineRange(spline_idx, spline_count, interval, r, "evaluateCubicSpline");
  const Tcalc dr = r - static_cast<Tcalc>(spline_idx);
  evaluateCubicSpline(abcd_coefficients[spline_idx], dr, u_contrib, fmag);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc>
void evaluateCubicSpline(const std::vector<Tcalc4> &abcd_coefficients, const Tcalc interval,
                         const Tcalc r, Tcalc *u_contrib, Tcalc *fmag) {
  evaluateCubicSpline(abcd_coefficients.data(), abcd_coefficients.size(), interval, r, u_contrib,
                      fmag);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc>
void evaluateCubicSpline(const Hybrid<Tcalc4> &abcd_coefficients, const Tcalc interval,
                         const Tcalc r, Tcalc *u_contrib, Tcalc *fmag) {
  evaluateCubicSpline(abcd_coefficients.data(), abcd_coefficients.size(), interval, r, u_contrib,
                      fmag);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateCubicSpline(const Tcalc4* abcd_coefficients, const int spline_count,
                           const Tcalc interval, const Tcalc r) {
  const Tcalc r_itvl = r / interval;
  const int spline_idx = r_itvl;
  checkSplineRange(spline_idx, spline_count, interval, r, "evaluateCubicSpline");
  const Tcalc dr = r_itvl - static_cast<Tcalc>(spline_idx);
  return evaluateCubicSpline(abcd_coefficients[spline_idx], dr);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateCubicSpline(const std::vector<Tcalc4> &abcd_coefficients, const int spline_count,
                           const Tcalc interval, const Tcalc r) {
  return evaluateCubicSpline(abcd_coefficients.data(), spline_count, interval, r);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateCubicSpline(const Hybrid<Tcalc4> &abcd_coefficients, const int spline_count,
                           const Tcalc interval, const Tcalc r) {
  return evaluateCubicSpline(abcd_coefficients.data(), spline_count, interval, r);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc>
void evaluateQuarticSpline(const Tcalc4 abcd_coefficients, const Tcalc e_coefficient,
                           const Tcalc r, Tcalc *u_contrib, Tcalc *fmag) {
  *u_contrib += (((((((abcd_coefficients.x * r) + abcd_coefficients.y) * r) +
                    abcd_coefficients.z) * r) + abcd_coefficients.w) * r) + e_coefficient;
  if (fmag != nullptr) {
    if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
      *fmag += ((((((4.0 * abcd_coefficients.x * r) + (3.0 * abcd_coefficients.y)) * r) +
                  (2.0 * abcd_coefficients.z)) * r) + abcd_coefficients.w) / r;
    }
    else {
      *fmag += ((((((4.0f * abcd_coefficients.x * r) + (3.0f * abcd_coefficients.y)) * r) +
                  (2.0f * abcd_coefficients.z)) * r) + abcd_coefficients.w) / r;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateQuarticSpline(const Tcalc4 abcd_coefficients, const Tcalc e_coefficient,
                             const Tcalc r) {
  Tcalc2 result;
  result.x = (((((((abcd_coefficients.x * r) + abcd_coefficients.y) * r) +
                 abcd_coefficients.z) * r) + abcd_coefficients.w) * r) + e_coefficient;
  if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
    result.y = (((((4.0 * abcd_coefficients.x * r) + (3.0 * abcd_coefficients.y)) * r) +
                 (2.0 * abcd_coefficients.z)) * r) + abcd_coefficients.w;
  }
  else {
    result.y = (((((4.0f * abcd_coefficients.x * r) + (3.0f * abcd_coefficients.y)) * r) +
                 (2.0f * abcd_coefficients.z)) * r) + abcd_coefficients.w;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc>
void evaluateQuarticSpline(const Tcalc4* abcd_coefficients, const Tcalc* e_coefficient,
                           const int spline_count, const Tcalc interval, const Tcalc r,
                           Tcalc *u_contrib, Tcalc *fmag) {
  const Tcalc r_itvl = r / interval;
  const int spline_idx = r_itvl;
  checkSplineRange(spline_idx, spline_count, interval, r, "evaluateQuarticSpline");
  const Tcalc dr = r - static_cast<Tcalc>(spline_idx);
  evaluateQuarticSpline(abcd_coefficients[spline_idx], e_coefficient[spline_idx], dr, u_contrib,
                        fmag);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc>
void evaluateQuarticSpline(const std::vector<Tcalc4> &abcd_coefficients,
                           const std::vector<Tcalc> &e_coefficient, const Tcalc interval,
                           const Tcalc r, Tcalc *u_contrib, Tcalc *fmag) {
  evaluateQuarticSpline(abcd_coefficients.data(), e_coefficient.data(), e_coefficient.size(),
                        interval, r, u_contrib, fmag);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc>
void evaluateQuarticSpline(const Hybrid<Tcalc4> &abcd_coefficients,
                           const Hybrid<Tcalc> &e_coefficient, const Tcalc interval,
                           const Tcalc r, Tcalc *u_contrib, Tcalc *fmag) {
  evaluateQuarticSpline(abcd_coefficients.data(), e_coefficient.data(), e_coefficient.size(),
                        interval, r, u_contrib, fmag);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateQuarticSpline(const Tcalc4* abcd_coefficients, const Tcalc* e_coefficient,
                             const int spline_count, const Tcalc interval, const Tcalc r) {
  const Tcalc r_itvl = r / interval;
  const int spline_idx = r_itvl;
  checkSplineRange(spline_idx, spline_count, interval, r, "evaluateQuarticSpline");
  const Tcalc dr = r_itvl - static_cast<Tcalc>(spline_idx);
  return evaluateQuarticSpline(abcd_coefficients[spline_idx], e_coefficient[spline_idx], dr);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateQuarticSpline(const std::vector<Tcalc4> &abcd_coefficients,
                             const std::vector<Tcalc> &e_coefficient, const int spline_count,
                             const Tcalc interval, const Tcalc r) {
  return evaluateQuarticSpline(abcd_coefficients.data(), e_coefficient.data(), spline_count,
                               interval, r);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateQuarticSpline(const Hybrid<Tcalc4> &abcd_coefficients,
                             const Hybrid<Tcalc> &e_coefficient, const int spline_count,
                             const Tcalc interval, const Tcalc r) {
  return evaluateWuarticSpline(abcd_coefficients.data(), e_coefficient.data(), spline_count,
                               interval, r);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
void evaluateQuinticSpline(const Tcalc4 abcd_coefficients, const Tcalc2 ef_coefficients,
                           const Tcalc r, Tcalc *u_contrib, Tcalc *fmag) {
  *u_contrib += (((((((((abcd_coefficients.x * r) + abcd_coefficients.y) * r) +
                      abcd_coefficients.z) * r) + abcd_coefficients.w) * r) +
                  ef_coefficients.x) * r) + ef_coefficients.y;
  if (fmag != nullptr) {
    if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
      *fmag += ((((((((5.0 * abcd_coefficients.x * r) + (4.0 * abcd_coefficients.y)) * r) +
                    (3.0 * abcd_coefficients.z)) * r) + (2.0 * abcd_coefficients.w)) * r) +
                ef_coefficients.x) / r;
    }
    else {
      *fmag += ((((((((5.0f * abcd_coefficients.x * r) + (4.0f * abcd_coefficients.y)) * r) +
                    (3.0f * abcd_coefficients.z)) * r) + (2.0f * abcd_coefficients.w)) * r) +
                ef_coefficients.x) / r;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc3, typename Tcalc2, typename Tcalc>
Tcalc3 evaluateQuinticSpline(const Tcalc4 abcd_coefficients, const Tcalc2 ef_coefficients,
                             const Tcalc r) {
  Tcalc3 result;
  result.x = (((((((((abcd_coefficients.x * r) + abcd_coefficients.y) * r) +
                   abcd_coefficients.z) * r) + abcd_coefficients.w) * r) +
               ef_coefficients.x) * r) + ef_coefficients.y;
  if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
    result.y = (((((((5.0 * abcd_coefficients.x * r) + (4.0 * abcd_coefficients.y)) * r) +
                   (3.0 * abcd_coefficients.z)) * r) + (2.0 * abcd_coefficients.w)) * r) +
               ef_coefficients.x;
    result.z = (((((20.0 * abcd_coefficients.x * r) + (12.0 * abcd_coefficients.y)) * r) +
                 (6.0 * abcd_coefficients.z)) * r) + (2.0 * abcd_coefficients.w);
  }
  else {
    result.y = (((((((5.0f * abcd_coefficients.x * r) + (4.0f * abcd_coefficients.y)) * r) +
                   (3.0 * abcd_coefficients.z)) * r) + (2.0f * abcd_coefficients.w)) * r) +
               ef_coefficients.x;
    result.z = (((((20.0f * abcd_coefficients.x * r) + (12.0f * abcd_coefficients.y)) * r) +
                 (6.0f * abcd_coefficients.z)) * r) + (2.0f * abcd_coefficients.w);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
void evaluateQuinticSpline(const Tcalc4* abcd_coefficients, const Tcalc2* ef_coefficients,
                           const int spline_count, const Tcalc interval, const Tcalc r,
                           Tcalc *u_contrib, Tcalc *fmag) {
  const Tcalc r_itvl = r / interval;
  const int spline_idx = r_itvl;
  checkSplineRange(spline_idx, spline_count, interval, r, "evaluateQuinticSpline");
  const Tcalc dr = r - static_cast<Tcalc>(spline_idx);
  evaluateQuinticSpline(abcd_coefficients[spline_idx], ef_coefficients[spline_idx], dr,
                        u_contrib, fmag);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
void evaluateQuinticSpline(const std::vector<Tcalc4> &abcd_coefficients,
                           const std::vector<Tcalc2> &ef_coefficients, const Tcalc interval,
                           const Tcalc r, Tcalc *u_contrib, Tcalc *fmag) {
  evaluateQuinticSpline(abcd_coefficients.data(), ef_coefficients.data(), ef_coefficients.size(),
                        interval, r, u_contrib, fmag);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
void evaluateQuinticSpline(const Hybrid<Tcalc4> &abcd_coefficients,
                           const Hybrid<Tcalc2> &ef_coefficients, const Tcalc interval,
                           const Tcalc r, Tcalc *u_contrib, Tcalc *fmag) {
  evaluateQuinticSpline(abcd_coefficients.data(), ef_coefficients.data(), ef_coefficients.size(),
                        interval, r, u_contrib, fmag);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc3, typename Tcalc2, typename Tcalc>
Tcalc3 evaluateQuinticSpline(const Tcalc4* abcd_coefficients, const Tcalc* ef_coefficients,
                             const int spline_count, const Tcalc interval, const Tcalc r) {
  const Tcalc r_itvl = r / interval;
  const int spline_idx = r_itvl;
  checkSplineRange(spline_idx, spline_count, interval, r, "evaluateQuinticSpline");
  const Tcalc dr = r_itvl - static_cast<Tcalc>(spline_idx);
  return evaluateQuinticSpline(abcd_coefficients[spline_idx], ef_coefficients[spline_idx], dr);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc3, typename Tcalc2, typename Tcalc>
Tcalc3 evaluateQuinticSpline(const std::vector<Tcalc4> &abcd_coefficients,
                             const std::vector<Tcalc> &ef_coefficients, const int spline_count,
                             const Tcalc interval, const Tcalc r) {
  return evaluateQuinticSpline(abcd_coefficients.data(), ef_coefficients.data(), spline_count,
                               interval, r);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc3, typename Tcalc2, typename Tcalc>
Tcalc3 evaluateQuinticSpline(const Hybrid<Tcalc4> &abcd_coefficients,
                             const Hybrid<Tcalc> &ef_coefficients, const int spline_count,
                             const Tcalc interval, const Tcalc r) {
  return evaluateQuinticSpline(abcd_coefficients.data(), ef_coefficients.data(), spline_count,
                               interval, r);
}

} // namespace stmath
} // namespace stormm
