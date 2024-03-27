// -*-c++-*-
#ifndef STORMM_MATH_ONE_DIMENSIONAL_SPLINES_H
#define STORMM_MATH_ONE_DIMENSIONAL_SPLINES_H

#include "copyright.h"

namespace stormm {
namespace stmath {

/// \brief Evaluate the value and, if requested, the derivative due to a tabulated cubic spline.
///        The spline has the form U = A(dr)^3 + B(dr)^2 + C(dr) + D, where dr is the relevant
///        distance less the left-hand limit of the relevant spline interval (the spline works in
///        local frames of reference to keep the numerics as well-conditioned as possible).  The
///        supplied inter-particle distance will be checked to ensure that it lies within the
///        applicable range of the spline.  The quartic functions will have continuous values, and
///        first derivatives at the boundaries between each segment.
///
/// Overloaded:
///   - Return the value and derivative as the "x" and "y" members of a tuple, respectively.
///   - Find the relevant spline interval and verify the spline's range, or assume that the
///     interval begins at zero and that the provided distance lies within its range.
///   - Operate on C-style arrays with a trusted number of spline intervals.
///   - Operate on Standard Template Library vectors (the length of the vectors indicates the
///     size of the spline table).
///   - Operate on Hybrid objects.
///
/// \param abcd_coefficients  Array of coefficient sets A, B, C, and D for the defining equation
///                           above, in members "x", "y", "z", and "w" of the tuple
/// \param spline_count       The number of spline intervals--this, combined with interval, will be
///                           used to determine the spline interval spacing.
/// \param interval           The switching distance at which the spline table becomes relevant
/// \param r                  Distance, e.g. the inter-paricle distance in units of Angstroms
/// \param u_contrib          Energy (or, more generally, function value) tally, modified and
///                           returned
/// \param fmag               Magnitude of the resulting derivative (accumulated and returned).
///                           The contribution of the spline will be divided by the distance r to
///                           prepare for the customary chain rule application that produces the
///                           force's Cartesian components.
/// \{
template <typename Tcalc4, typename Tcalc>
void evaluateCubicSpline(const Tcalc4 abcd_coefficients, Tcalc r, Tcalc *u_contrib,
                         Tcalc *fmag = nullptr);

template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateCubicSpline(const Tcalc4 abcd_coefficients, Tcalc r);

template <typename Tcalc4, typename Tcalc>
void evaluateCubicSpline(const Tcalc4* abcd_coefficients, int spline_count, Tcalc interval,
                         Tcalc r, Tcalc *u_contrib, Tcalc *fmag = nullptr);

template <typename Tcalc4, typename Tcalc>
void evaluateCubicSpline(const std::vector<Tcalc4> &abcd_coefficients, Tcalc interval, Tcalc r,
                         Tcalc *u_contrib, Tcalc *fmag = nullptr);

template <typename Tcalc4, typename Tcalc>
void evaluateCubicSpline(const Hybrid<Tcalc4> &abcd_coefficients, Tcalc interval, Tcalc r,
                         Tcalc *u_contrib, Tcalc *fmag = nullptr);

template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateCubicSpline(const Tcalc4* abcd_coefficients, int spline_count, Tcalc interval,
                           Tcalc r);

template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateCubicSpline(const std::vector<Tcalc4> &abcd_coefficients, Tcalc interval,
                           Tcalc r);

template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateCubicSpline(const Hybrid<Tcalc4> &abcd_coefficients, Tcalc interval,
                           Tcalc r);
/// \}

/// \brief Evaluate the potential and, if requested, the force due to a tabulated quartic spline.
///        The spline has the form U = A(dr)^4 + B(dr)^3 + C(dr)^2 + D(dr) + E, where dr is the
///        inter-particle distance less the left-hand limit of the relevant spline interval (the
///        spline works in local frames of reference to keep the numerics as well-conditioned as
///        possible).  The supplied inter-particle distance will be checked to ensure that it lies
///        within the applicable range of the spline.  The quartic functions will have continuous
///        values, and first derivatives at the boundaries between each segment.  Overloading and
///        descriptions of input arguments follow from evaluateCubicSpline() above, in addition to:
/// 
/// \param e_coefficient      Array of coefficients E from the defining equation above
/// \{
template <typename Tcalc4, typename Tcalc>
void evaluateQuarticSpline(const Tcalc4 abcd_coefficients, const Tcalc e_coefficient,
                           Tcalc r, Tcalc *u_contrib, Tcalc *fmag = nullptr);

template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateQuarticSpline(const Tcalc4 abcd_coefficients, const Tcalc e_coefficient, Tcalc r);

template <typename Tcalc4, typename Tcalc>
void evaluateQuarticSpline(const Tcalc4* abcd_coefficients, const Tcalc* e_coefficient,
                           int spline_count, Tcalc interval, Tcalc r, Tcalc *u_contrib,
                           Tcalc *fmag = nullptr);

template <typename Tcalc4, typename Tcalc>
void evaluateQuarticSpline(const std::vector<Tcalc4> &abcd_coefficients,
                           const std::vector<Tcalc> &e_coefficient, Tcalc interval, Tcalc r,
                           Tcalc *u_contrib, Tcalc *fmag = nullptr);

template <typename Tcalc4, typename Tcalc>
void evaluateQuarticSpline(const Hybrid<Tcalc4> &abcd_coefficients,
                           const Hybrid<Tcalc> &e_coefficient, Tcalc interval, Tcalc r,
                           Tcalc *u_contrib, Tcalc *fmag = nullptr);

template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateQuarticSpline(const Tcalc4* abcd_coefficients, const Tcalc* e_coefficient,
                             int spline_count, Tcalc interval, Tcalc r);

template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateQuarticSpline(const std::vector<Tcalc4> &abcd_coefficients,
                             const std::vector<Tcalc> &e_coefficient, Tcalc interval, Tcalc r);

template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evaluateQuarticSpline(const Hybrid<Tcalc4> &abcd_coefficients,
                             const Hybrid<Tcalc> &e_coefficient, Tcalc interval, Tcalc r);
/// \}

/// \brief Evaluate the potential and, if requested, the force due to a tabulated quartic spline.
///        The spline has the form U = A(dr)^4 + B(dr)^3 + C(dr)^2 + D(dr) + E, where dr is the
///        inter-particle distance less the left-hand limit of the relevant spline interval (the
///        spline works in local frames of reference to keep the numerics as well-conditioned as
///        possible).  The supplied inter-particle distance will be checked to ensure that it lies
///        within the applicable range of the spline.  The quartic functions will have continuous
///        values, and first derivatives at the boundaries between each segment.  Overloading and
///        descriptions of input arguments follow from evaluateCubicSpline() and
///        evaluateQuarticSpline() above.  For quintic splines, the value, first, and second
///        derivatives of the splined function can be held consistent at each knot and are thus
///        returned.
/// \{
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
void evaluateQuinticSpline(const Tcalc4 abcd_coefficients, const Tcalc2 ef_coefficients,
                           Tcalc r, Tcalc *u_contrib, Tcalc *fmag = nullptr);

template <typename Tcalc4, typename Tcalc3, typename Tcalc2, typename Tcalc>
Tcalc3 evaluateQuinticSpline(const Tcalc4 abcd_coefficients, const Tcalc2 ef_coefficients,
                             Tcalc r);

template <typename Tcalc4, typename Tcalc2, typename Tcalc>
void evaluateQuinticSpline(const Tcalc4* abcd_coefficients, const Tcalc2* ef_coefficients,
                           int spline_count, Tcalc interval, Tcalc r, Tcalc *u_contrib,
                           Tcalc *fmag = nullptr);

template <typename Tcalc4, typename Tcalc2, typename Tcalc>
void evaluateQuinticSpline(const std::vector<Tcalc4> &abcd_coefficients,
                           const std::vector<Tcalc2> &ef_coefficients, Tcalc interval, Tcalc r,
                           Tcalc *u_contrib, Tcalc *fmag = nullptr);

template <typename Tcalc4, typename Tcalc2, typename Tcalc>
void evaluateQuinticSpline(const Hybrid<Tcalc4> &abcd_coefficients,
                           const Hybrid<Tcalc2> &ef_coefficients, Tcalc interval, Tcalc r,
                           Tcalc *u_contrib, Tcalc *fmag = nullptr);

template <typename Tcalc4, typename Tcalc3, typename Tcalc2, typename Tcalc>
Tcalc3 evaluateQuinticSpline(const Tcalc4* abcd_coefficients, const Tcalc* ef_coefficients,
                             int spline_count, Tcalc interval, Tcalc r);

template <typename Tcalc4, typename Tcalc3, typename Tcalc2, typename Tcalc>
Tcalc3 evaluateQuinticSpline(const std::vector<Tcalc4> &abcd_coefficients,
                             const std::vector<Tcalc> &ef_coefficients, Tcalc interval, Tcalc r);

template <typename Tcalc4, typename Tcalc3, typename Tcalc2, typename Tcalc>
Tcalc3 evaluateQuinticSpline(const Hybrid<Tcalc4> &abcd_coefficients,
                             const Hybrid<Tcalc> &ef_coefficients, Tcalc interval, Tcalc r);
/// \}

} // namespace stmath 
} // namespace stormm

#include "one_dimensional_splines.tpp"

#endif
