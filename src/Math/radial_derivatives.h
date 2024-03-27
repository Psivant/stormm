// -*-c++-*-
#ifndef STORMM_RADIAL_DERIVATIVES_H
#define STORMM_RADIAL_DERIVATIVES_H

#include "copyright.h"
#include "Constants/scaling.h"
#include "math_enumerators.h"

namespace stormm {
namespace stmath {

/// \brief Compute the first partial derivative of a radially symmetric function based on the
///        first derivative of the function at the radius of interest and the displacement along
///        the direction of interest.
///
/// \param dfunc  The first derivative of the function of interest, (d/dr) [ u(r) ]
/// \param disp   Displacement along the direction of interest (e.g. i ~ x, y, or z in the
///               Catesian coordinate system)
/// \param r      The radius at which the function is being evaluated
template <typename T> T radialFirstDerivative(T dfunc, T disp, T r);

/// \brief Compute the second partial derivative of a radially symmetric function based on its
///        first and second derivatives at the radius of interest and the displacements along the
///        directions of interest.
///
/// Overloaded:
///   - Compute the second derivative along a single axis
///   - Compute the mixed partial derivative along two axes
///
/// \param dfunc   The first derivative of the function of interest, (d/dr) [ u(r) ]
/// \param ddfunc  The second derivative of the function of interest, (d2/dr2) [ u(r) ]
/// \param disp    Displacement along the one direction of interest (e.g. i ~ x, y, or z in the
///                Catesian coordinate system)
/// \param disp_x  Displacement along the first direction of interest (e.g. i ~ x, y, or z in the
///                Catesian coordinate system)
/// \param disp_y  Displacement along the second direction of interest, must be different from the
///                direction of disp_x
/// \param r       The radius at which the function is being evaluated
/// \param r2      The squared radius at which the function is being evaluated
/// \{
template <typename T> T radialSecondDerivative(T dfunc, T ddfunc, T disp, T r);

template <typename T> T radialSecondDerivative(T dfunc, T ddfunc, T disp_x, T disp_y, T r,
                                               T r2);
/// \}

/// \brief Compute the thrid partial derivative of a radially symmetric function based on its
///        first, second, and third derivatives at the radius of interest and the displacements
///        along the directions of interest.
///
/// Overloaded:
///   - Compute the third derivative along a single axis
///   - Compute the mixed partial derivative along two axes (two differentiations will take place
///     along the first axis, one along the second)
///   - Compute the mixed partial derivative along all three axes
///
/// Descriptions of input parameters follow from radialSecondDerivative() above, in addition to:
///
/// \param dddfunc  The third derivative of the function of interest, (d3/dr3) [ u(r) ]
/// \param disp_z   Displacement along the second direction of interest, must be different from the
///                 directions of the previous two displacements
/// \{
template <typename T> T radialThirdDerivative(T dfunc, T ddfunc, T dddfunc, T disp, T r,
                                              T r2);

template <typename T> T radialThirdDerivative(T dfunc, T ddfunc, T dddfunc, T disp_x, T disp_y,
                                              T r, T r2);

template <typename T> T radialThirdDerivative(T dfunc, T ddfunc, T dddfunc, T disp_x, T disp_y,
                                              T disp_z, T r, T r2);
/// \}

/// \brief Evaluate partial derivatives of a radially symmetric function.
///
/// \param du     The first derivative of the function in question
/// \param d2u    Second derivative of the function in question (may be set to zero if only first
///               derivatives are needed)
/// \param d3u    Third derivative of the function in question (may be set to zero if only the
//                first or second derivatives are needed)
/// \param dx     Displacement along the Cartesian X direction
/// \param dy     Displacement along the Cartesian Y direction
/// \param dz     Displacement along the Cartesian Z direction
/// \param r      Overall displacement
/// \param r2     Square of the overall displacement
/// \param order  The particular derivative to evaluate (choices for all partial derivatives up
///               to and including the third derivative tensor)
template <typename T> T radialPartialDerivative(T du, T d2u, T d3u, T dx, T dy, T dz, T r, T r2,
                                                FunctionLevel order);
  
} // namespace stmath
} // namespace stormm

#include "radial_derivatives.tpp"

#endif
