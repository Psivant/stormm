// -*-c++-*-
#ifndef STORMM_APPROX_H
#define STORMM_APPROX_H

#include <cmath>
#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Reporting/error_format.h"
#include "unit_test_enumerators.h"

namespace stormm {
namespace testing {

using data_types::isScalarType;

/// \brief Class for handling comparisons of floating-point results versus expected values in unit
///        tests.  The general form is to take a vector of real numbers and apply the comparison
///        across all entries.
class Approx {
public:

  /// \brief Constructors for real number comparisons based on a vector of multiple input values.
  ///        Other floating point and integer vectors can be converted to double-precision vectors.
  ///
  /// \param values_in  The value for comparison
  /// \param style_in   Type of comparison, relative or absolute
  /// \param tol_in     Tolerance for deviations that will still yield a correct comparison
  template <typename T> Approx(const std::vector<T> &values_in,
                               ComparisonType style_in = ComparisonType::ABSOLUTE,
                               double tol_in = 1.0e-6);

  /// \brief Constructor for real number comparisons based on a single input value
  ///
  /// Overloaded:
  ///   - Take the comparison style ahead of the tolerance (convenient use of default tolerance)
  ///   - Take the tolerance ahead of the comparison style (tolerance must be specified)
  ///
  /// \param value_in  The value for comparison
  /// \param style_in  The type of comparison, relative or absolute
  /// \param tol_in    Tolerance for deviations that will still yield a correct comparison
  /// \{
  Approx(double value_in, ComparisonType style_in = ComparisonType::ABSOLUTE,
         double tol_in = 1.0e-6);
  Approx(double value_in, double tol_in, ComparisonType style_in = ComparisonType::ABSOLUTE);
  /// \}

  /// \brief Default destructor
  ~Approx() = default;

  /// \brief Get the size of the collection of numbers in a real number comparison
  int size() const;

  /// \brief Get the expectation value, sans any tolerance
  double getValue() const;

  /// \brief Get all values of an approximate comparison, sans any tolerance
  std::vector<double> getValues() const;

  /// \brief Get the style used in this approximate comparison
  ComparisonType getStyle() const;

  /// \brief Get the margin or tolerance associated with this approximate comparison
  /// \{
  double getMargin() const;
  double getTolerance() const;
  double getTol() const;
  /// \}

  /// \brief Set the values that form the basis of the approximate comparison.
  ///
  /// Overloaded:
  ///   - Set the approximation to use a single value
  ///   - Set the target value of a particular index
  ///   - Set the approximation to use multiple values 
  ///
  /// \param value_in   The single value to set (this will change the length of the comparison
  ///                   vector to one)
  /// \param values_in  The vector of values to set
  /// \param index      Index of the value to change (this will be checked against the length of
  ///                   comparison values in the existing object)
  /// \{
  void setValue(double value_in);
  void setValues(const std::vector<double> &values_in);
  void setValue(double value_in, size_t index);
  /// \}

  /// \brief Set the tolerance of the existing object.  This will not emit a new object.
  ///
  /// Overloaded:
  ///   - Three different names that individual developers may find most intuitive
  ///
  /// \param dtol_in  The tolerance to use
  /// \{
  void setMargin(double dtol_in);
  void setTolerance(double dtol_in);
  void setTol(double dtol_in);
  /// \}

  /// \brief Set the tolerance.   This is written in such a way as to mimic the Catch2 unit
  ///        testing framework in some circumstances.  It will not change the original object's
  ///        tolerance, rather it will emit an object with the same comparison values and the
  ///        desired tolerance.
  ///
  /// Overloaded:
  ///   - Three different names that individual developers may find most intuitive
  ///
  /// \param dtol_in  The tolerance to use
  /// \{
  Approx margin(double dtol_in) const;
  Approx tolerance(double dtol_in) const;
  Approx tol(double dtol_in) const;
  /// \}

  /// \brief Test whether a real-valued scalar is the same as the value held for comparison.
  ///
  /// \param test_value  The scalar value to test
  bool test(const double test_value) const;

  /// \brief Test whether a vector of real-valued scalars is the same as as vector held for
  ///        comparison.
  ///
  /// \param test_values  The vector of values to test
  template <typename T> bool test(const std::vector<T> &test_values) const;

private:
  std::vector<double> values;  ///< Reference values stored for later testing
  ComparisonType style;        ///< The type of comparison to make, i.e. ABSOLUTE or RELATIVE
  double dtol;                 ///< Tolerance for a succcessful test
};

/// \brief Check the type and size of vectors for approximate comparisons.  This encapsulates
///        checks that would have to be done in several contexts.
///
/// \param test_values  The vector of values to test
/// \param cr           Approximate object to test against
template <typename T> bool verifyVectorApproxCompatibility(const std::vector<T> &test_values,
                                                           const Approx &cr);

/// \brief Overload the == operator to accommodate a Approx object and a float or double scalar
///        or vector
///
/// \param d   The scalar to compare (everything becomes double precision, so integers to 53 bits)
/// \param cr  Reference data with an associated tolerance
/// \{
bool operator==(const double d, const Approx &cr);
bool operator==(const Approx &cr, const double d);
bool operator!=(const double d, const Approx &cr);
bool operator!=(const Approx &cr, const double d);
bool operator>(const double d, const Approx &cr);
bool operator>(const Approx &cr, const double d);
bool operator<(const double d, const Approx &cr);
bool operator<(const Approx &cr, const double d);
bool operator>=(const double d, const Approx &cr);
bool operator>=(const Approx &cr, const double d);
bool operator<=(const double d, const Approx &cr);
bool operator<=(const Approx &cr, const double d);
/// \}

/// \brief Overloads for equal and non equal relational operators using templated vectors
///
/// \param tvec  The vector to compare
/// \param cr    Reference data with an associated tolerance
/// \{
template <typename T> bool operator==(const std::vector<T> &tvec, const Approx &cr);
template <typename T> bool operator==(const Approx &cr, const std::vector<T> &tvec);
template <typename T> bool operator!=(const std::vector<T> &tvec, const Approx &cr);
template <typename T> bool operator!=(const Approx &cr, const std::vector<T> &tvec);
template <typename T> bool operator>(const std::vector<T> &tvec, const Approx &cr);
template <typename T> bool operator>(const Approx &cr, const std::vector<T> &tvec);
template <typename T> bool operator<(const std::vector<T> &tvec, const Approx &cr);
template <typename T> bool operator<(const Approx &cr, const std::vector<T> &tvec);
template <typename T> bool operator>=(const std::vector<T> &tvec, const Approx &cr);
template <typename T> bool operator>=(const Approx &cr, const std::vector<T> &tvec);
template <typename T> bool operator<=(const std::vector<T> &tvec, const Approx &cr);
template <typename T> bool operator<=(const Approx &cr, const std::vector<T> &tvec);
/// \}


} // namespace testing
} // namespace stormm

#include "approx.tpp"

#endif
