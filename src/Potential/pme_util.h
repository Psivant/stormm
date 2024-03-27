// -*-c++-*-
#ifndef STORMM_PME_UTIL_H
#define STORMM_PME_UTIL_H

#include <cmath>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "Math/math_enumerators.h"
#include "Math/radial_derivatives.h"

namespace stormm {
namespace energy {

using card::Hybrid;
using data_types::isSignedIntegralScalarType;
using stmath::FunctionLevel;
using stmath::radialFirstDerivative;
using stmath::radialSecondDerivative;
using stmath::radialThirdDerivative;
using stmath::radialPartialDerivative;
  
/// \brief Default settings for PME calculations
/// \{
constexpr int default_charge_mapping_order = 5;
constexpr double default_dsum_tol = 1.0e-5;
constexpr double max_dsum_tol = 1.0e-4;
constexpr double default_pme_cutoff = 8.0;
constexpr double default_pme_grid_spacing_target = 1.0;
/// \}

/// \brief A recommended maximum for any Ewald coefficient, obtained for a cutoff of 40.0A and
///        direct sum tolerance 1.0e-5.  The Gaussian density (charge, dispersion) smearing is
///        exceptionally broad.
constexpr double minimum_ewald_coefficient = 0.0625;

/// \brief A recommended maximum for any Ewald coefficient, obtained for cutoffs shorter than
///        4.0A and direct sum tolerances of 1.0e-8 or lower.  Density is hardly spread at all.
constexpr double maximum_ewald_coefficient = 1.0;

/// \brief Compute the "Ewald coefficient" for a given particle-particle ("direct-space")
///        interaction model.
///
/// \param cutoff          The cutoff distance for particle-particle interactions
/// \param direct_sum_tol  The direct sum tolerance
double ewaldCoefficient(double cutoff, double direct_sum_tol);

/// \brief Compute the Gaussian spread of density (charge or otherwise) associated with a specific
///        particle-particle interaction cutoff and direct sum tolerance.  Descriptions of input
///        parameters follow from ewaldCoefficient() above, and in fact this function merely calls
///        that and returns 0.5 over the result.
double pmeGaussianSpread(double cutoff, double direct_sum_tol);

/// \brief Back out the direct sum tolerance implied by a given Ewald coefficient and
///        particle-particle interaction cutoff.
///
/// \param cutoff             The particle-particle interaction cutoff
/// \param ewald_coefficient  The multiplier to the range argument in the error function when
///                           applied to the potential, e.g. "a" in erf(a * r) / r for a Coulomb
///                           potential
double recoverDirectSumTolerance(double cutoff, double ewald_coefficient);

/// \brief For a given Ewald coefficient and Coulomb's constant, compute the derivative of the
///        particle-particle interaction of two protons.
///
/// Overloaded:
///   - Work along one dimension, or with two particles at arbitrary locations
///   - Return a single value
///   - Return a triplet of values based on three distinct orders
///
/// \param pa          Location of the first particle
/// \param pb          Location of the second particle
/// \param ew_coeff    The Ewald coefficient to apply.  It is fed in as a double-precision number
///                    so that a constant factor involved in the derivative calculations can be
///                    computed in as high a precision as possible, but for other purposes the
///                    Ewald coefficient will be cast to the calculation data type.
/// \param kcoul       Coulomb's constant
/// \param r           The distance between particles
/// \param r2          The squared distance between particles.  This is fed in alongside the
///                    distance because both will be available in most "real world" settings where
///                    the PME derivative is being computed, and moreover r will have been derived
///                    from r2, meaning that additional error would likely be incurred by
///                    back-calculating r2 from r.
/// \param derivative  Indicate whether to solve for the function value or some derivative, up to
///                    the third derivative
/// \param order       Indicate whether to compute an energy (integer 0, or the VALUE enumeration),
///                    a first, second or third order radial derivative (integer 1, 2, or 3), or a
///                    specific partial derivative (whatever enumeration)
/// \{
template <typename Tcalc>
Tcalc elecPMEDirectSpace(double ew_coeff, Tcalc kcoul, Tcalc r, Tcalc r2, int order);

template <typename Tcalc, typename T3>
Tcalc elecPMEDirectSpace(const T3 pa, const T3 pb, double ew_coeff, Tcalc kcoul,
                         FunctionLevel order = FunctionLevel::VALUE);

template <typename Tcalc, typename T3>
std::vector<Tcalc> elecPMEDirectSpace(const T3 pa, const T3 pb, double ew_coeff, Tcalc kcoul,
                                      const std::vector<FunctionLevel> &orders);
/// \}

} // namespace energy
} // namespace stormm

#include "pme_util.tpp"

#endif