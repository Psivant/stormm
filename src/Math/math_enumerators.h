// -*-c++-*-
#ifndef STORMM_MATH_ENUMERATORS_H
#define STORMM_MATH_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace stmath {

/// \brief A list of the different boundary components that determine a tricubic spline, covering
///        the values and all derivatives at the boundaries of the grid element.
enum class FunctionLevel {
  VALUE,   ///< Value of the function at a particular grid point
  DX,      ///< Derivative of the function along the Cartesian X direction
  DY,      ///< Derivative of the function along the Cartesian Y direction
  DZ,      ///< Derivative of the function along the Cartesian Z direction
  DXX,     ///< Double derivative of the function along the Cartesian X direction
  DXY,     ///< Cross-derivative of the function along the Cartesian X and Y directions
  DXZ,     ///< Cross-derivative of the function along the Cartesian X and Z directions
  DYY,     ///< Double derivative of the function along the Cartesian Y direction
  DYZ,     ///< Cross-derivative of the function along the Cartesian Y and Z directions
  DZZ,     ///< Double derivative of the function along the Cartesian Z direction
  DXXX,    ///< Triple derivative of the function along the Cartesian X direction
  DXXY,    ///< Triple partial derivative of the function along the Cartesian X and Y directions
  DXXZ,    ///< Triple partial derivative of the function along the Cartesian X and Z directions
  DXYY,    ///< Triple partial derivative of the function along the Cartesian X and Y directions
  DXZZ,    ///< Triple partial derivative of the function along the Cartesian X and Z directions
  DXYZ,    ///< Triple partial derivative of the function along all three Cartesian directions
  DYYY,    ///< Triple derivative of the function along the Cartesian Y direction
  DYYZ,    ///< Triple partial derivative of the function along the Cartesian Y and Z directions
  DYZZ,    ///< Triple partial derivative of the function along the Cartesian Y and Z directions
  DZZZ     ///< Triple derivative of the function along the Cartesian Z direction
};

/// \brief For tricubic interpolation, there are many choices as to the 64 pieces of information
///        that will define the polynomial coefficients for any particular mesh element.  The
///        traditional approach taken by Francois Lekien and Jerrold Marsden favors maximal
///        continuity at the corners of the interpolated region of space, taking the values of the
///        function f(x,y,z), its first derivatives df/dx, df/dy, and df/dz, and its mixed partial
///        second and third derivatives d2f/dxy, d2f/dxz, d2f/dyz, and d3f/dxyz.  Despite employing
///        third derivatives, the resulting interpolated function is only C1 continuous--the first
///        derivative alone is guaranteed to be smooth throughout space.
///
///        As Lekien and Marsden showed, fitting local polynomial splines to the function values
///        and first derivatives at each mesh vertex is sufficient to guarantee C1 continuity
///        across mesh element boundaries.  Fitting to higher-order mixed partial derivatives
///        grants only partial C2 and C3 continuity.  d2f/dx2, d2f/dy2, and d2f/dz2, as well as
///        various third derivatives, are not guaranteed to be continuous if there is any degree of
///        approximation in the interpolant.  Furthermore, imposing these restrictions requires the
///        interpolated function to have C2 and C3 continuity of its own.  Finally, results suggest
///        the numerical accuracy of this interpolant is only strong in orthorhombic meshes.
///        Outside of this case, the higher-order derivatives appear to lead to significant
///        over-corrections in the interiors of mesh elements.
///
///        To avoid the limitations and improve the accuracy of tricubic interpolation in the
///        interiors of mesh elements, a second style of interpolant is provided which evaluates
///        additional values of the function f(x,y,z) at selected points in the interior of each
///        mesh element.
/// \{
enum class Interpolant {
  SMOOTHNESS,     ///< Evaluate f(x,y,z), its first derivatives, and various mixed partial
                  ///<   derivatives at mesh corners until sufficient criteria are found to fit the
                  ///<   necessary coefficients for the interpolant.
  FUNCTION_VALUE  ///< Evaluate f(x,y,z) and its first derivatives at each mesh vertex, then make
                  ///<   additional evaluations of f(x,y,z) at points in the interior of the
                  ///<   element to get the most accurate reproduction of f(x,y,z) throughout
                  ///<   space.  For tricubic interpolation, these additional f(x,y,z) evaluations
                  ///<   begin once C1 continuity is guaranteed.
};
/// \}

/// \brief Enumerate the ways to approach a limit.
enum class LimitApproach {
  BELOW,  ///< Approach the limit from below, with the function argument ascending
  ABOVE   ///< Approach the limit from above, with the function argument descending
};

/// \brief Enumerate various functions that can be mapped to a logarthimic spline table by
///        hard-coded methods.
enum class LogSplineForm {
  ELEC_PME_DIRECT,        ///< Electrostatic PME direct space functional form (this will include
                          ///<   Coulomb's constant, in the program's internal units)
  ELEC_PME_DIRECT_EXCL,   ///< Electrostatic PME direct space functional form, plus a Coulombic
                          ///<   exclusion
  DELEC_PME_DIRECT,       ///< Electrostatic PME direct space functional form.  This will include
                          ///<   Coulomb's constant, in the program's internal units, as well as an
                          ///<   additional factor of the inverse interparticle distance, in the
                          ///<   interest of rolling as many factors as possible together in
                          ///<   double-precision before creating the spline.
  DELEC_PME_DIRECT_EXCL,  ///< Electrostatic PME direct space functional form, plus a Coulombic
                          ///<   exclusion.  See above for additional factors in this quantity.
  CUSTOM                  ///< Custom potential drawn from developer-supplied data tables
};

/// \brief Enumerate various methods by which a function can be indexed.
enum class TableIndexing {
  ARG,           ///< Indexing into the spline tables is by the value of the function argument
  SQUARED_ARG,   ///< Indexing is based on the square of the value of the function argument
  ARG_OFFSET,    ///< Indexing into the spline tables is by the value of the function argument plus
                 ///<   an offset (some power of 2), which eliminates all of the table elements
                 ///<   with applicable ranges smaller than the offset
  SQ_ARG_OFFSET  ///< Indexing into the spline tables is base on the square of the value of the
                 ///<   function argument plus an offset (some power of 2)
};

/// \brief List different combinations of basis functions that might be used to construct splines
///        for interpolation.
enum class BasisFunctions {
  MIXED_FRACTIONS,  ///< The basis functions include a fractional series, e.g. x, 1, 1/x, 1/x^2
  POLYNOMIAL        ///< The basis functions are a polynomial series, e.g. 1, x, x^2, x^3
};

/// \brief B-spline coefficients are, by definition, a smooth partition of unity.  However, that
///        result is only guaranteed to the precision of the calculating machine, and even then
///        there are many ways of achieving or enforcing it.  List those methods.
enum class BSplineUnity {
  CENTER_FILL,  ///< One of the central coefficients, the largest for odd orders of interpolation
                ///<   and likely to be the largest for even orders, will be left uncalculated at
                ///<   each step of the recursive calculation.  This coefficient will be computed
                ///<   subtracting all others from 1.0 in the final step.
  NONE          ///< The standard recursive procedure is applied to all coefficients, allowing the
                ///<   sum to diverge from unity as roundoff error seeps in over each iteration.
};

/// \brief Real-to-complex and complex-to-real FFTs can be staged as in-place or out-of-place
///        transformations, with advntages in terms of speed, memory conservation, or both
///        depending on the card and the size of the problem.
enum class FFTMode {
  IN_PLACE,     ///< Perform the FFTs "in place." Code that implements such functionality must
                ///<   allocate the arrays with additional padding along the "x" (or, unit cell
                ///<   A) axis, as ((nx / 2) + 1).
  OUT_OF_PLACE  ///< Perform the FFTs "out of place", dumping complex results in an R2C transform
                ///<   into an array with nx * ny * nz elements for real data of dimensions nx, ny,
                ///<   and nz.
};
  
/// \brief Get a human-readable string describing an enumeration of the provided type.  Various
///        overloads of this function serve enumerators across many libraries.
///
/// \param input  The enumeration to translate
/// \{
std::string getEnumerationName(FunctionLevel input);
std::string getEnumerationName(Interpolant input);
std::string getEnumerationName(LimitApproach input);
std::string getEnumerationName(LogSplineForm input);
std::string getEnumerationName(TableIndexing input);
std::string getEnumerationName(BasisFunctions input);
std::string getEnumerationName(BSplineUnity input);
std::string getEnumerationName(FFTMode input);
/// \}

/// \brief Translate a human-readable string (likely from user input) into one of the modes for
///        real-to-complex and complex-to-real FFTs.
FFTMode translateFFTMode(const std::string &input);

/// \brief Translate a human-readable string (likely from user input) into one of the known
///        functions for logarithmically indexed spline tables.
LogSplineForm translateLogSplineForm(const std::string &input);

/// \brief Translate a human-readable string (likely from user input) into one of the known
///        table indexing methods for logarithmically indexed spline tables.
TableIndexing translateTableIndexing(const std::string &input);
  
} // namespace stmath
} // namespace stormm

#endif
