// -*-c++-*-
#ifndef STORMM_POLYNUMERIC_H
#define STORMM_POLYNUMERIC_H

#include <vector>
#include "copyright.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Reporting/error_format.h"
#include "parsing_enumerators.h"

namespace stormm {
namespace parse {

using data_types::getStormmHpcVectorTypeName;
using data_types::isFloatingPointScalarType;
using data_types::isHpcVectorType;
using data_types::isScalarType;
using data_types::isSignedIntegralScalarType;
using data_types::isUnsignedIntegralScalarType;

/// \brief Union for storing numbers or other 4-8 byte pieces of information
union PolyNumeric {
  double d;
  float f;
  float2 f2;
  int i;
  unsigned int ui;
  long long int lli;
  unsigned long long int ulli;
  char4 c4;
};

/// \brief Pull a double-precision real vector out of a PolyNumeric vector.
///
/// \param values  The STL vector of PolyNumeric data
std::vector<double> doubleFromPolyNumeric(const std::vector<PolyNumeric> &values);

/// \brief Pull an integer vector out of a PolyNumeric vector.
///
/// \param values  The STL vector of PolyNumeric data
std::vector<int> intFromPolyNumeric(const std::vector<PolyNumeric> &values);

/// \brief Pull a long long integer vector out of a PolyNumeric vector.
///
/// \param values  The STL vector of PolyNumeric data
std::vector<llint> llintFromPolyNumeric(const std::vector<PolyNumeric> &values);

/// \brief Pull an unsigned integer vector out of a PolyNumeric vector.
///
/// \param values  The STL vector of PolyNumeric data
std::vector<uint> uintFromPolyNumeric(const std::vector<PolyNumeric> &values);

/// \brief Pull an unsigned long long integer vector out of a PolyNumeric vector.
///
/// \param values  The STL vector of PolyNumeric data
std::vector<ullint> ullintFromPolyNumeric(const std::vector<PolyNumeric> &values);

/// \brief Pull a vector of char4 tuples out of a PolyNumeric vector.
///
/// \param values  The STL vector of PolyNumeric data
std::vector<char4> char4FromPolyNumeric(const std::vector<PolyNumeric> &values);

/// \brief Put a vector of char4s into PolyNumeric format.  This is a special case of the templated
///        function of the same name, to avoid collisions between the member variables of the char4
///        type and scalar quantities which have no member variables.
///
/// \param values  The vector of interest
std::vector<PolyNumeric> polyNumericVector(const std::vector<char4> &values);

/// \brief Produce a string describing the data type of each of the NumberFormat enumerations.
///        This is useful for describing the contents of the associated PolyNumeric union type.
///
/// \param cform  The format setting
std::string nameNumericalType(NumberFormat cform);

/// \brief Produce a string corresponding to each of the NumberFormat enumerations.  This is
///        useful for describing the contents of the associated PolyNumeric union type.
///
/// \param cform  The format setting
std::string nameNumberFormat(NumberFormat cform);

/// \brief Put a vector of any scalar type into PolyNumeric format.
///
/// \param values  The vector of interest
template <typename T> std::vector<PolyNumeric> polyNumericVector(const std::vector<T> &values);

/// \brief Extract a formatted number from a character array.  Returns a multi-valent number type
///        (union) which can then be interpreted in several ways.
///
/// \param a           The character string
/// \param read_begin  Point in the string at which to being reading the number (defaults to the
///                    start of the string)
/// \param len         Expected length of the number (default 0, continues until the end of the
///                    string)
PolyNumeric extractFormattedNumber(const char* a, NumberFormat cform, int read_begin = 0,
                                   int len = 0);

} // namespace parse
} // namespace stormm

#include "polynumeric.tpp"

#endif
