// -*-c++-*-
#ifndef STORMM_FORMULAS_H
#define STORMM_FORMULAS_H

#include "copyright.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"

namespace stormm {
namespace stmath {

/// \brief Compute the factorial of an integer.  Strict limits are placed on this function, which
///        operates in integer math, to ensure that it does not overflow.
///
/// \param x  Integer for which to compute the factorial (limit 12)
int factorial(int x);

/// \brief Compute the factorial of an integer.  This can handle larger limits than the factorial()
///        function but is still somewhat limited.
///
/// \param x  Integer for which to compute the factorial (limit 20)
llint factoriall(llint x);

/// \brief Compute the factorial of an integer.  This can handle the largest integers using the
///        incomplete gamma function.
///
/// \param x  Integer for which to compute the factorial (limit 170)
double factoriald(int x);

/// \brief Raise a 32-bit signed integer to an integer power.  This routine is not protected
///        against numerical overflow, and uses recursion.
///
/// \param x  The integer to exponentiate
/// \param p  The power to which the integer shall be raised
int ipow(int x, int p);

/// \brief Raise a 64-bit signed integer to an integer power.  This routine is not protected
///        against numerical overflow, and uses recursion.
///
/// \param x  The integer to exponentiate
/// \param p  The power to which the integer shall be raised
llint ipowl(llint x, int p);

/// \brief Compute the value and derivatives of a sigmoidal function of the form
///
///        S(r) = 1 / (exp(p * (r - r0)) + 1)
///
/// Overloaded:
///   - Compute the value plus first, second, and third derivatives, returning the value as a tuple
///     with the value in the "x" member and first through third derivatives in the "y", "z", and
///     "w" members, respectively
///   - Compute one of the above values alone
///
/// \param r           Argument to the sigmoidal function, as shown in the defining equation above
/// \param crossover   The value of r at which the sigmoidal form S(r) crosses 0.5.  This is the
///                    point about which S(r) is odd, represented as r0 in the equation above.
/// \param intensity   The steepness parameter controlling the rate at which the sigmoidal function
///                    makes the transition.  This is the represented as p in the equation above.
/// \param order       Order of the derivative to compute (specify zero for the function value)
/// \{
double4 sigmoid(double r, double crossover, double intensity);
double sigmoid(double r, double crossover, double intensity, int order);
/// \}

/// \brief Single-precision variant of sigmoid, above.  Overloading and descriptions of parameters
///        follows from the double-precision variant.
/// \{
float4 sigmoidf(float r, float crossover, float intensity);
float sigmoidf(float r, float crossover, float intensity, int order);
/// \}

/// \brief Compute the value of an angle between the ray connecting the origin to a point in a
///        two-dimensional plane with the ray along the first principal axis of that plane's
///        coordinate system.  The angle is returned in radians.
///
/// \param x  The coordinate along the first dimension
/// \param y  The coordinate along the second dimension
double angleOnAxes(double x, double y);

/// \brief Single-precision variant of angleOnAxes, above.  Overloading and descriptions of
///        parameters follows from the double-precision variant.
float angleOnAxesf(float x, float y);
  
} // namespace stmath
} // namespace stormm

#endif
