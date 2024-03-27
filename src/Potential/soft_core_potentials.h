// -*-c++-*-
#ifndef STORMM_SOFT_CORE_POTENTIALS_H
#define STORMM_SOFT_CORE_POTENTIALS_H

#include <cmath>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "Math/matrix_ops.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "energy_enumerators.h"

namespace stormm {
namespace energy {

using card::Hybrid;
using parse::NumberFormat;
using parse::realToString;
using stmath::invertSquareMatrix;
using stmath::matrixVectorMultiply;
using stmath::qrSolver;

/// \brief The maximum number of quartic splines that can define any tabulated softcore function
constexpr int max_quartic_spline_count = 64;

/// \brief Compute the electrostatic contribution for an interaction forced to go linear at a
///         specified distance.  This function will have a continuous potential and first
///         derivative at the crossover point.
///
/// \param r               Distance between the two particles
/// \param clash_distance  Absolute distance between any two particles, below which they would be
///                        declared to clash
/// \param qiqj            Product of the two particles' charges and Coulomb's constant, attenuated
///                        by any other prefactors such as 1:4 scaling terms.
/// \param ele_contrib     The energy contribution from this specific interaction.  This will be
///                        accumulated and returned.
/// \param fmag            The magnitude of the force between two particles.  This will be
///                        accumulated and returned.  Providing a nullptr will cause the force
///                        calculation to be skipped.
template <typename Tcalc>
void linearCoreElectrostatics(const Tcalc r, const Tcalc clash_distance, const Tcalc qiqj,
                              Tcalc *ele_contrib, Tcalc *fmag = nullptr);

/// \brief Compute the electrostatic contribution for an interaction damped to become a quadratic
///        function of the inter-particle distance with a relative maximum at -1.0 Angstroms (the
///        location of the maximum ensures that the function derivative has a finite, nonzero value
///        even as the inter-particle distance approaches zero).  If the inter-particle distance is
///        too close to zero, no force will be logged.  Descriptions of parameters follow from
///        linearCoreElectrostatics() above.  This function will have a continuous potential and
///        first derivative at the crossover point, but the second derivative will bear a
///        discontinuity.
template <typename Tcalc>
void quadraticCoreElectrostatics(const Tcalc r, const Tcalc clash_distance, const Tcalc qiqj,
                                 Tcalc *ele_contrib, Tcalc *fmag = nullptr);

/// \brief Compute the electrostatic contributions for an interaction damped to become a cubic
///        function of the inter-particle distance.  The potential and first three derivatives are
///        all maintained during the switch to the softcore form.  Descriptions of parameters
///        follow from linearCoreElectrostatics() above.
template <typename Tcalc>
void cubicCoreElectrostatics(const Tcalc r, const Tcalc clash_distance, const Tcalc qiqj,
                             Tcalc *ele_contrib, Tcalc *fmag = nullptr);

/// \brief Compute the Lennard-Jones contribution for an interaction damped to become a cubic
///        function of the inter-particle distance.  The potential and first three derivatives are
///        all maintained during the switch to the softcore form.
///
/// \param r            Distance between the two particles
/// \param clash_ratio  The minimum ratio of the distance between any two particles and the pair's
///                     Lennard-Jones sigma parameter.  If the inter-particle distance is too low,
///                     the softcore potential will engage.
/// \param lja          The Lennard-Jones B parameter, pre-scaled by any attenuation constants such
///                     as 1:4 prefactors
/// \param ljb          The Lennard-Jones B parameter, pre-scaled by any attenuation constants
/// \param vdw_contrib  The energy contribution from this specific interaction.  This will be
///                     accumulated and returned.
/// \param fmag         The magnitude of the force between two particles.  This will be accumulated
///                     and returned.  Providing a nullptr will skip the force calculation.
template <typename Tcalc>
void cubicCoreLennardJones(const Tcalc r, const Tcalc clash_ratio, const Tcalc lja,
                           const Tcalc ljb, Tcalc *vdw_contrib, Tcalc *fmag = nullptr);
  
/// \brief Compute the Lennard-Jones contribution for an interaction damped to become a quartic
///        function with a relative maximum at -1.0 Angstroms (the location of the maximum ensures
///        that the function derivative has a finite, nonzero value even as the inter-particle
///        distance approaches zero).  Descriptions of parameters follow from
///        cubicCoreLennardJones() above.
template <typename Tcalc>
void quarticCoreLennardJones(const Tcalc r, const Tcalc clash_ratio, const Tcalc lja,
                             const Tcalc ljb, Tcalc *vdw_contrib, Tcalc *fmag = nullptr);

/// \brief Adjust the target value of the splined function for the left-hand side of a quartic
///        spline segment.  This supports the quarticCoreSplineTable() function below.
///
/// \param interval_targets  The list of spline knots, modified and returned
/// \param index             The spline knot to be adjusted
/// \param f_rswitch         The value of the original function at the point where the handoff to
///                          the softcore potential occurs
/// \param target_zero       The original value targeted for the softcore potential to reach as
///                          the inter-particle distance goes to zero
/// \param move_increment    The amount by which to adjust the left-hand target for the interval in
///                          question, if the first derivative fails to fall in line with the
///                          overall intended trend for the spline.
void adjustIntervalTarget(std::vector<double> *interval_targets, int index, double f_rswitch,
                          double target_zero, double move_increment);

/// \brief Evaluate the first derivative of a quartic polynomial.
///
/// \param coeffs  The first four coefficients for the function Ar^4 + Br^3 + Cr^2 + Dr + E. The
///                "x", "y", "z", and "w" members of the tuple are A, B, C, and D, respectively.
/// \param r       The point at whcih to evaluate the function
double evaluateQuarticFirstDerivative(const double4 coeffs, double r);

/// \brief Compute a softcore potential for an inter-particle potential function based on a
///        cubic spline.  The cubic spline will target a slope for the softcore function as the
///        inter-particle distance approaches zero, respecting the value of the first derivative
///        at the point where the original function takes over.  Furthermore, the softcore
///        potential will hold the second derivative to the smallest possible value as the
///        inter-particle distance approaches zero.  The output of this function is a tuple,
///        holding four coefficients for the spline.
///
///        SC(r) = A(r)^3 + B(r)^2 + C(r) + D
///
/// \param abcd_coefficients  The first four coefficients for each quartic spline in the table,
///                           Ar^4 + Br^3 + Cr^2 + Dr + E.  Filled and returned.
/// \param rswitch            The interparticle distance for switching to the splined softcore
///                           potential
/// \param f_rswitch          Value of the original function at the switching distance rswitch
/// \param df_rswitch         Value of the original function's first derivative at rswitch
/// \param target_zero        Target value for the softcore function's derivative at r = 0
template <typename Tcalc4>
void cubicSoftCore(Tcalc4* abcd_coefficients, double rswitch, double f_rswitch, double df_rswitch,
                   double target_zero = 0.0);

/// \brief Compute a softcore potential for an inter-particle potential function based on a
///        quartic spline.  The quartic spline will target a slope for the softcore function
///        as the inter-particle distance approaches zero, respecting the value of the first and
///        second derivatives at the point where the original function takes over.  Furthermore,
///        the softcore potential will hold the second derivative to the smallest possible value
///        as the inter-particle distance approaches zero.  The output of this function is a tuple
///        and one additional constant, holding five coefficients for the spline.
///
///        SC(r) = A(r)^4 + B(r)^3 + C(r)^2 + D(r) + E
///
/// Descriptions of parameters follow from cubicSoftCore() above, with the addition of:
///
/// \param e_coefficient      The fifth coefficient for the quartic splines, filled and returned
/// \param d2f_rswitch        Value of the original function's second derivative at rswitch
template <typename Tcalc4, typename Tcalc>
void quarticSoftCore(Tcalc4* abcd_coefficients, Tcalc* e_coefficient, double rswitch,
                     double f_rswitch, double df_rswitch, double d2f_rswitch,
                     double target_zero = 0.0);

/// \brief Compute a softcore potential for an inter-particle potential function based on a
///        quintic spline.  The quintic spline will target a slope for the softcore function
///        as the inter-particle distance approaches zero and respect the values of the first,
///        second, and third derivatives at the point where the original function takes over.
///        Furthermore, the potential will have a minimal second derivative as r goes to zero.
///        The output of this function is a pair of tuples holding five coefficients
///        for the spline:
///
///        SC(r) = A(r)^5 + B(r)^4 + C(r)^3 + D(r)^2 + E(r) + F
///
/// Descriptions of parameters follow from cubicSoftCore() and quarticSoftCore above, with the
/// addition of:
///
/// \param ef_coefficients    The fifth and sixth coefficients for the quartic splines, filled and
///                           returned
/// \param d3f_rswitch        Value of the original function's third derivative at rswitch
template <typename Tcalc4, typename Tcalc2>
void quinticSoftCore(Tcalc4* abcd_coefficients, Tcalc2* ef_coefficients, double rswitch,
                     double f_rswitch, double df_rswitch, double d2f_rswitch, double d3f_rswitch,
                     double target_zero = 0.0);

/// \brief Compute a softcore potential based on exponentials.  The function has the form
///
///            U = Ae^(2v(x - rswitch)) + Be^(3v(x - rswitch)) + C e^(4v(x - rswitch)) + D
///
///        where k is a constant defined in the input.  The function will seek to fit the value
///        and the first three derivatives of the function at the handoff point.
///
/// \param abcd_coefficients  The coefficients A, B, C, and D from the equation above, computed
///                           and returned in the "x", "y", "z", and "w" coefficients of the tuple
/// \param v                  The exponential factor v in the equation above
/// \param rswitch            The point below which the standard potential transitions to the
///                           softcore form
/// \param f_rswitch          Value of the original potential function at the switching distance
/// \param df_rswitch         Derivative of the original potential at the switching distance
/// \param d2f_rswitch        Second derivative of the original potential at the switching distance
/// \param d3f_rswitch        Third derivative of the original potential at the switching distance
template <typename Tcalc4>
void exponentialSoftCore(Tcalc4 *abcd_coefficients, double v, double rswitch, double f_rswitch,
                         double df_rswitch, double d2f_rswitch, double d3f_rswitch);

} // namespace energy
} // namespace stormm

#include "soft_core_potentials.tpp"

#endif
