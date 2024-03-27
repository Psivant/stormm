// -*-c++-*-
#ifndef STORMM_MULTIPLICATION_H
#define STORMM_MULTIPLICATION_H

#include "copyright.h"
#include <vector>
#include "Accelerator/hybrid.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Reporting/error_format.h"

namespace stormm {
namespace stmath {

using card::Hybrid;
using data_types::getStormmScalarTypeName;
using data_types::isScalarType;
using data_types::isSignedIntegralScalarType;
using data_types::isUnsignedIntegralScalarType;
using data_types::isFloatingPointScalarType;
using parse::NumberFormat;
using parse::realToString;
  
/// \brief Compute the (natural) logarithm of a product series.
///
/// Overloaded:
///   - Operate on a C-style array of trusted length
///   - Operate on a std::vector of original values
///   - Operate on a Hybrid object of original values
///
/// \param values  The array of values for which to take the product
/// \param length  Trusted length of the array values (if a C-style array is used)
/// \{
template <typename T> double logProduct(const T* values, const size_t length);
template <typename T> double logProduct(const std::vector<T> &values);
template <typename T> double logProduct(const Hybrid<T> &values);
/// \}

/// \brief Compute the product series (pi notation) result from an array of scalar numbers.  If
///        presented with an integral type, the result will not be checked against overflow.  If
///        presented with a floating point scalar type, the range of the developing product will
///        be tracked to ensure that it remains within the limits of the hardware.  If there is a
///        violation, there will be a last-ditch attempt to compute the result by exponentiating
///        the sum of the natural logarithms of the absolute values of each element in the array.
///
/// Overloaded:
///   - Operate on a C-style array with a trusted length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The array of numbers
/// \param length  Trusted length of va (if va is a C-style array)
/// \{
template <typename Tprod, typename Tbase> Tprod seriesProduct(const Tbase* va,
                                                              const size_t length);
template <typename Tprod, typename Tbase> Tprod seriesProduct(const std::vector<Tbase> &va);
template <typename Tprod, typename Tbase> Tprod seriesProduct(const Hybrid<Tbase> &va);
/// \}

} // namespace stmath
} // namespace stormm

#include "multiplication.tpp"

#endif
