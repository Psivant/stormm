// -*-c++-*-
#ifndef STORMM_VECTOR_OPS_H
#define STORMM_VECTOR_OPS_H

#include <algorithm>
#include <cmath>
#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "DataTypes/mixed_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Numerics/split_fixed_precision.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "UnitTesting/approx.h"
#include "sorting_enumerators.h"
#include "statistics.h"

namespace stormm {
namespace stmath {

using card::Hybrid;
using constants::ExceptionResponse;
using data_types::getHpcVectorTypeSize;
using data_types::isHpcVectorType;
using data_types::isFloatingPointScalarType;
using data_types::isScalarType;
using data_types::isSignedIntegralScalarType;
using data_types::isUnsignedIntegralScalarType;
using data_types::ValueWithCounter;
using numerics::hostChangeFPBits;
using numerics::hostDoubleToInt95;
using numerics::hostInt95ToDouble;
using parse::NumberFormat;
using parse::realToString;
using testing::Approx;
  
/// \brief Check that two vectors of the same data type are compatible for various arithmetic
///        comparisons.  Throws and exception if the vectors cannot be compared.
///
/// Overloaded:
///   - Check two Standard Template Library vectors
///   - Check two Hybrid objects
///
/// \param va  The first vector
/// \param vb  The second vector
/// \{
template <typename T> void vectorComparisonCheck(const std::vector<T> &va,
                                                 const std::vector<T> &vb, const char* caller);
template <typename T> void vectorComparisonCheck(const Hybrid<T> &va, const Hybrid<T> &vb,
                                                 const char* caller);
/// \}

/// \brief Maximum absolute deviation between two vectors
///
/// Overloaded:
///   - Operate on two C-style vectors of a stated (and trusted) length
///   - Operate on two Standard Template Library vectors
///   - Operate on to Hybrid objects
///
/// \param va      The first vector
/// \param vb      The second vector, taken as the reference
/// \param length  The length of the C-style arrays
/// \{
template <typename T> T maxAbsoluteDifference(const T* va, const T* vb, size_t length);
template <typename T> T maxAbsoluteDifference(const std::vector<T> &va, const std::vector<T> &vb);
template <typename T> T maxAbsoluteDifference(const Hybrid<T> &va, const Hybrid<T> &vb);
/// \}

/// \brief Maximum relative deviation between two vectors.  The signed deviation with the largest
///        absolute value will be returned.
///
/// Overloaded:
///   - Operate on two C-style vectors of a stated (and trusted) length
///   - Operate on two Standard Template Library vectors
///   - Operate on to Hybrid objects
///
/// \param va      The first vector
/// \param vb      The second vector, taken as the reference
/// \param length  The length of the C-style arrays
/// \{
template <typename T> double maxRelativeDifference(const T* va, const T* vb, size_t length);
template <typename T> double maxRelativeDifference(const std::vector<T> &va,
                                                   const std::vector<T> &vb);
template <typename T> double maxRelativeDifference(const Hybrid<T> &va, const Hybrid<T> &vb);
/// \}

/// \brief Compute the mean unsigned error between two vectors.
///
/// Overloaded:
///   - Operate on two C-style vectors of a stated (and trusted) length
///   - Operate on two Standard Template Library vectors
///   - Operate on to Hybrid objects
///
/// \param va      The first vector
/// \param vb      The second vector, taken as the reference
/// \param length  The length of the C-style arrays
/// \{
template <typename T> double meanUnsignedError(const T* va, const T* vb, size_t length);
template <typename T> double meanUnsignedError(const std::vector<T> &va, const std::vector<T> &vb);
template <typename T> double meanUnsignedError(const Hybrid<T> &va, const Hybrid<T> &vb);
/// \}

/// \brief Compute the root mean squared error between two vectors.
///
/// Overloaded:
///   - Operate on two C-style vectors of a stated (and trusted) length
///   - Operate on two Standard Template Library vectors
///   - Operate on to Hybrid objects
///
/// \param va      The first vector
/// \param vb      The second vector, taken as the reference
/// \param length  The length of the C-style arrays
/// \{
template <typename T> double rmsError(const T* va, const T* vb, size_t length);
template <typename T> double rmsError(const std::vector<T> &va, const std::vector<T> &vb);
template <typename T> double rmsError(const Hybrid<T> &va, const Hybrid<T> &vb);
/// \}

/// \brief Compute the root mean squared error between two vectors, normalized by the mean value
///        of the reference vector.
///
/// Overloaded:
///   - Operate on two C-style vectors of a stated (and trusted) length
///   - Operate on two Standard Template Library vectors
///   - Operate on to Hybrid objects
///
/// \param va      The first vector
/// \param vb      The second vector, taken as the reference
/// \param length  The length of the C-style arrays
/// \{
template <typename T> double relativeRmsError(const T* va, const T* vb, size_t length);
template <typename T> double relativeRmsError(const std::vector<T> &va, const std::vector<T> &vb);
template <typename T> double relativeRmsError(const Hybrid<T> &va, const Hybrid<T> &vb);
/// \}

/// \brief Compute the mean value of a vector of familiar data types.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \{
template <typename T> double mean(const T* va, const size_t length);
template <typename T> double mean(const std::vector<T> &va);
template <typename T> double mean(const Hybrid<T> &va);
/// \}

/// \brief Compute the variance of a vector of familiar data types.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \param method  Statistical method for computing a variance-related quantity
/// \{
template <typename T> double variance(const T* va, size_t length, VarianceMethod method);
template <typename T> double variance(const std::vector<T> &va, VarianceMethod method);
template <typename T> double variance(const Hybrid<T> &va, VarianceMethod method);
/// \}

/// \brief Compute the mean value of a vector of two-tuple familiar data types.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \{
template <typename T> double2 bivariateMean(const T* va, size_t length);
template <typename T> double2 bivariateMean(const std::vector<T> &va);
template <typename T> double2 bivariateMean(const Hybrid<T> &va);
/// \}

/// \brief Compute the mean value of a vector of three-tuple familiar data types.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \{
template <typename T> double3 trivariateMean(const T* va, size_t length);
template <typename T> double3 trivariateMean(const std::vector<T> &va);
template <typename T> double3 trivariateMean(const Hybrid<T> &va);
/// \}

/// \brief Compute the mean value of a vector of three-tuple familiar data types.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \{
template <typename T> double4 quadrivariateMean(const T* va, size_t length);
template <typename T> double4 quadrivariateMean(const std::vector<T> &va);
template <typename T> double4 quadrivariateMean(const Hybrid<T> &va);
/// \}

/// \brief Obtain the maximum value of a vector of scalar numbers.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \{
template <typename T> T maxValue(const T* va, size_t length);
template <typename T> T maxValue(const std::vector<T> &va);
template <typename T> T maxValue(const Hybrid<T> &va);
/// \}

/// \brief Obtain the minimum value of a vector of scalar numbers.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \{
template <typename T> T minValue(const T* va, size_t length);
template <typename T> T minValue(const std::vector<T> &va);
template <typename T> T minValue(const Hybrid<T> &va);
/// \}

/// \brief Obtain the largest absolute value of a vector of scalar numbers
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \{
template <typename T> T maxAbsValue(const T* va, size_t length);
template <typename T> T maxAbsValue(const std::vector<T> &va);
template <typename T> T maxAbsValue(const Hybrid<T> &va);
/// \}
  
/// \brief Compute the Pearson correlation coefficient for two vectors of numbers.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The first vector
/// \param vb      The second vector
/// \param length  The length of the C-style arrays
/// \{
template <typename T> double pearson(const T* va, const T* vb, size_t length);
template <typename T> double pearson(const std::vector<T> &va, const std::vector<T> &vb);
template <typename T> double pearson(const Hybrid<T> &va, const Hybrid<T> &vb);
/// \}

/// \brief Add a constant to all elements of a vector.
///
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \param inc     The amount by which to alter the vector (passed in by value with an implicit
///                conversion to the same type as the vector itself)
/// \{
template <typename T> void addScalarToVector(T* va, size_t length, T inc);
template <typename T> void addScalarToVector(std::vector<T> *va, T inc);
template <typename T> void addScalarToVector(Hybrid<T> *va, T inc);
/// \}

/// \brief Multiply all elements in a vector by a value.
/// 
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \param factor  The amount by which to scale the vector (passed in by value with an implicit
///                conversion to the same type as the vector itself)
/// \{
template <typename T> void elementwiseMultiply(T* va, size_t length, T factor);
template <typename T> void elementwiseMultiply(std::vector<T> *va, T factor);
template <typename T> void elementwiseMultiply(Hybrid<T> *va, T factor);
/// \}

/// \brief Divide all elements in a vector by a value.  This is needed for integer representations
///        which cannot work with non-integral divisors.
/// 
/// Overloaded:
///   - Operate on a C-style vectors of a stated (and trusted) length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  The length of the C-style array
/// \param factor  The amount by which to divide the vector (passed in by value with an implicit
///                conversion to the same type as the vector itself)
/// \{
template <typename T> void elementwiseDivide(T* va, size_t length, T factor);
template <typename T> void elementwiseDivide(std::vector<T> *va, T factor);
template <typename T> void elementwiseDivide(Hybrid<T> *va, T factor);
/// \}

/// \brief Compute the cross product of two vectors.  Each should be of dimension 3, and if any of
///        the vectors have more than 3 elements, only the first three will be used.  If supplying
///        vectors, no check is made to ensure the float or double nature of the type (that's up to
///        the developer, and allows for corner cases in which two vectors of ints could produce a
///        vector of other ints, provided that the format is not overflowed).  If supplying HPC 
///        tuples, the compiler will catch cases of bad input where two-element tuples are
///        supplied, but four-element tuples will act like vectors longer than three elements.
///
/// Overloaded:
///   - Operate on C-style arrays of trusted length
///   - Operate on Standard Template Library vectors
///   - Operate on three-element HPC tuples
///
/// \param va  Compute va x vb = vc.
/// \param vb  Compute va x vb = vc.
/// \param vc  Pre-allocated third vector to hold the result.  Compute va x vb = vc.
/// \{
template <typename T> void crossProduct(const T* va, const T* vb, T* vc);
template <typename T> void crossProduct(const std::vector<T> &va, const std::vector<T> &vb,
                                        std::vector<T> *vc);
template <typename T> T crossProduct(const T va, const T vb);
/// \}

/// \brief Compute the magnitude of a vector.  The result is always returned as a double-precision
///        value, and accumulated as such, but individual elements are not converted to double
///        prior to taking their squares.  This can overflow integer data types if used improperly.
///
/// Overloaded:
///   - Operate on a C-style array of trusted length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  Length of va
/// \{
template <typename T> double magnitude(const T* va, const size_t length);
template <typename T> double magnitude(const std::vector<T> &va);
template <typename T> double magnitude(const Hybrid<T> &va);
/// \}

/// \brief Compute the length of the hypotenuse from two lengths (assumed to be the lengths of
///        vectors at a right angle to one another) in a manner resistan to catastrophic
///        cancellation.
///
/// \param a  Length of one vector
/// \param b  Length of the other vector
template <typename T> T hypotenuse(const T a, const T b);
  
/// \brief Normalize a vector, storing the unit normal vector in its place and returning the
///        original magnitude.
///
/// Overloaded:
///   - Operate on a C-style array of trusted length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///
/// \param va      The vector of interest
/// \param length  Length of va  
/// \{
template <typename T> double normalize(T* va, const size_t length);
template <typename T> double normalize(std::vector<T> *va);
template <typename T> double normalize(Hybrid<T> *va);
/// \}
  
/// \brief Compute the dot product of two vectors.  The result is always returned as a
///        double-precision value, and accumulated as such, but individual elements are not
///        converted to double prior to multiplication.  This can overflow integer data types if
///        used improperly.
///
/// Overloaded:
///   - Operate on a pair of C-style arrays of trusted length
///   - Operate on a pair of Standard Template Library vectors
///   - Operate on a pair of Hybrid objects
///
/// \param va      The first vector of interest
/// \param vb      The second vector of interest
/// \param length  The first vector of interest
/// \{
template <typename T> double dot(const T* va, const T* vb, const size_t length);
template <typename T> double dot(const std::vector<T> &va, const std::vector<T> &vb);
template <typename T> double dot(const Hybrid<T> &va, const Hybrid<T> &vb);
/// \}

/// \brief Compute the angle between two vectors using the cosine formula.
///
/// Overloaded:
///   - Operate on a pair of C-style arrays of trusted length
///   - Operate on a pair of Standard Template Library vectors (the vectors will be checked for
///     consistent lengths)
///   - Operate on a pair of Hybrid objects (consistency checked)
///
/// \param va      The first vector of interest
/// \param vb      The second vector of interest
/// \param length  The length of both vectors
/// \{
template <typename T> double angleBetweenVectors(const T* va, const T* vb, const size_t length);
template <typename T> double angleBetweenVectors(const std::vector<T> &va,
                                                 const std::vector<T> &vb);
template <typename T> double angleBetweenVectors(const Hybrid<T> &va, const Hybrid<T> &vb);
/// \}
  
/// \brief Compute the projection of one vector onto another.
///
/// Overloaded:
///   - Operate on C-style arrays of trusted length
///   - Operate on Standard Template Library vectors (must be of consistent length)
///   - Operate on a pair of Hybrid objects (must be of consistent length)
///   - Return a new array or work with a provided, pre-existing one (only works with vector
///     objects, not raw pointers or Hybrids)
///
/// \param va      The projection of va onto vb will be computed 
/// \param vb      The projection of va onto vb will be computed 
/// \param vc      Store the result in a pre-allocated space vc
/// \param length  The length of all vectors
/// \{
template <typename T> void project(const T* va, const T* vb, T* vc, size_t length);
template <typename T> void project(const std::vector<T> &va, const std::vector<T> &vb,
                                   std::vector<T> *vc);
template <typename T> void project(const Hybrid<T> &va, const Hybrid<T> &vb, Hybrid<T> *vc);
template <typename T> std::vector<T> project(const std::vector<T> &va, const std::vector<T> &vb);
/// \}

/// \brief Compute the perpendicular component of a vector A which has no alignment with vector B.
///
/// Overloaded:
///   - Operate on C-style arrays of trusted length
///   - Operate on Standard Template Library vectors (must be of consistent length)
///   - Return a new array or work with a provided, pre-existing one (only works with vector
///     objects, not raw pointers or Hybrids)
///
/// \param va      The projection of va onto vb will be computed 
/// \param vb      The projection of va onto vb will be computed 
/// \param result  Store the result in this pre-allocated space
/// \param length  The length of all vectors
/// \{
template <typename T> void perpendicularComponent(const T* va, const T* vb, T* result,
                                                  const size_t length);
template <typename T> void perpendicularComponent(const std::vector<T> &va,
                                                  const std::vector<T> &vb,
                                                  std::vector<T> *result);
template <typename T> std::vector<T> perpendicularComponent(const std::vector<T> &va,
                                                            const std::vector<T> &vb);
/// \}

/// \brief Compute the distance between a point and a plane (both in three dimensions, not a
///        hyperplane, as the formula uses a cross-product).  This is done by computing a unit
///        normal vector to the plane and taking the projection of the point positional vector
///        onto the unit normal vector.  Requires vectors defining the plane and a point position
///        all in the same local coordinate frame.
///
/// Overloaded:
///   - Operate on a pair of C-style arrays of length trusted to be 3
///   - Operate on a pair of Standard Template Library vectors
///   - Operate on a pair of Hybrid objects
///
/// \param va      First defininig the plane
/// \param vb      Second vector defining the plane
/// \param pt_pos  Positional vector of the point
/// \{
template <typename T> double pointPlaneDistance(const T* va, const T* vb, const T* pt_pos);
template <typename T> double pointPlaneDistance(const std::vector<T> &va, const std::vector<T> &vb,
                                                const std::vector<T> &pt_pos);
template <typename T> double pointPlaneDistance(const Hybrid<T> &va, const Hybrid<T> &vb,
                                                const Hybrid<T> &pt_pos);
/// \}

/// \brief Report an error for the histogram / bin finding function below, if the requested value
///        exceeds the available range.
///
/// \param desc          Description to provide with an error or warning message
/// \param value         Value that violated the known bin limits
/// \param policy        Indicator of how to respond to this contingency
template <typename T> void reportBinLimitError(const std::string &desc, const T value,
                                               const ExceptionResponse policy);
  
/// \brief Find the bin to which a value belongs, given a monotonically increasing set of limits.
///        Returns -1 if the value is less than the leftmost bin limit, or the number of bins
///        plus 1 if the value is greater than the rightmost bin limit, if the exception handling
///        permits.
///
/// Overloaded:
///   - Operate on a C-style array of trusted length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object (host side only)
///
/// \param limits  Bin limits
/// \param value   Value to search for within one of the bins
/// \param length  Trusted length of a C-style array (the array is assumed to be length + 1)
/// \param policy  Indication of what to do if the value is outside the listed limits
/// \{
template <typename T> int findBin(const T* limits, const T value, const int length,
                                  const ExceptionResponse policy = ExceptionResponse::DIE);
template <typename T> int findBin(const std::vector<T> &limits, const T value,
                                  const ExceptionResponse policy = ExceptionResponse::DIE);
template <typename T> int findBin(const Hybrid<T> &limits, const T value,
                                  const ExceptionResponse policy = ExceptionResponse::DIE);
/// \}

/// \brief Search a data set for a particular value.  The index of the element containing the
///        value of interest is returned.  The data can be stated as occurring in a particular
///        order, which will be trusted (certain orderings will permit a binary search, otherwise
///        a direct linear search will occur).  Returns the length of the array if the value is
///        not found.
///
/// Overloaded:
///   - Operate on a C-style array
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object
///   - Search for an approximate value (return the first match meeting specified tolerance)
///
/// \param vdata   The array to search
/// \param value   The value to seek out in the array (might be approximate)
/// \param length  Length of the C-style array (will be trusted)
/// \param format  Statement about the form of the data to expect (will be trusted)
/// \{
template <typename T> size_t locateValue(const T* vdata, const T value, const size_t length,
                                         const DataOrder format = DataOrder::NONE);

template <typename T> size_t locateValue(const std::vector<T> &vdata, const T value,
                                         const DataOrder format = DataOrder::NONE);

template <typename T> size_t locateValue(const Hybrid<T> &vdata, const T value,
                                         const DataOrder format = DataOrder::NONE);

template <typename T> size_t locateValue(const Approx &value, const T* vdata, const size_t length,
                                         const DataOrder format = DataOrder::NONE);

template <typename T> size_t locateValue(const Approx &value, const std::vector<T> &vdata,
                                         const DataOrder format = DataOrder::NONE);

template <typename T> size_t locateValue(const Approx &value, const Hybrid<T> &vdata,
                                         const DataOrder format = DataOrder::NONE);
/// \}

/// \brief Find the unique entries in a vector.  The result will be returned in ascending order.
///        This will only operate on a Standard Template Library vector.
///
/// Overloaded:
///   - Perform an in-place operation on the vector, shrinking its size and changing its content
///     into the result
///   - Perform an out-of-place operation on the vector, leaving the original data unchanged.
///
/// \param va  The vector to sort
/// \{
template <typename T> std::vector<T> reduceUniqueValues(const std::vector<T> &va);
template <typename T> void reduceUniqueValues(std::vector<T> *va);
/// \}

/// \brief Find the unique elements of one vector that are not present in another.  Options are
///        provided to search for replaceable or separate copies of each value.  The function will
///        proceed by making Standard Template Library vectors of the unique values in each array
///        as it finds them and keeping tallies of the number of unique copies of each value that
///        it finds.  The returned vector will show the unique unmatched values, coupled with a
///        positive count if the unique value is present some number of times in the first vector
///        but not in the second, or a negative count if the unique value is present some number
///        of times in the second vector but not the first.
///
/// Overloaded:
///   - Operate on Standard Template Library vectors
///   - Operate on C-style arrays (copied internally to std::vectors for sorting, the original
///     data will be unchanged)
///   - Operate on a pair of Hybrid objects (copied internally to std::vectors for sorting, the
///     original data will be unchanged)
///   
/// \param va             Search for elements of va that are not present in vb
/// \param vb             Search for elements of vb that are not present in va
/// \param length_a       The length of va, if va is a C-style array
/// \param length_b       The length of vb, if vb is a C-style array
/// \param check_repeats  Indication of whether to check that the same numbers of copies of each
///                       unique value are present in both arrays, or to simply verify that each
///                       unique value is covered at least once by both arrays
/// \{
template <typename T> std::vector<ValueWithCounter<T>>
findUnmatchedValues(const std::vector<T> &va, const std::vector<T> &vb,
                    UniqueValueHandling check_repeats = UniqueValueHandling::UNIQUE_VALUES_ONLY);

template <typename T> std::vector<ValueWithCounter<T>>
findUnmatchedValues(const T* va, const T* vb, size_t length_a, size_t length_b,
                    UniqueValueHandling check_repeats = UniqueValueHandling::UNIQUE_VALUES_ONLY);

template <typename T> std::vector<ValueWithCounter<T>>
findUnmatchedValues(const Hybrid<T> &va, const Hybrid<T> &vb,
                    UniqueValueHandling check_repeats = UniqueValueHandling::UNIQUE_VALUES_ONLY);
/// \}

/// \brief Log a result in an array of unsigned integers to accumulate a bitmask.
///
/// Overloaded:
///   - Use unsigned integers
///   - Use long long unsigned integers
///   - Use short unsigned integers
///   - Operate on a C-style array trusted to be of sufficient length, or a Standard Template
///     Library vector (the size will be checked)
///
/// \param va   The developing bit mask
/// \param pos  The index of the bit mask to color by setting the bit to 1
/// \{
void accumulateBitmask(uint* va, const size_t pos);
void accumulateBitmask(std::vector<uint> *va, const size_t pos);
void accumulateBitmask(ullint* va, const size_t pos);
void accumulateBitmask(std::vector<ullint> *va, const size_t pos);
void accumulateBitmask(ushort* va, const size_t pos);
void accumulateBitmask(std::vector<ushort> *va, const size_t pos);
/// \}

/// \brief Uncheck a result in an array of unsigned integers to modify a bitmask.
///
/// Overloaded:
///   - Use unsigned integers
///   - Use long long unsigned integers
///   - Use short unsigned integers
///   - Operate on a C-style array trusted to be of sufficient length, or a Standard Template
///     Library vector (the size will be checked)
///
/// \param va   The developing bit mask
/// \param pos  The index of the bit mask to color by setting the bit to 0
/// \{
void unsetBitInMask(uint* va, const size_t pos);
void unsetBitInMask(std::vector<uint> *va, const size_t pos);
void unsetBitInMask(ullint* va, const size_t pos);
void unsetBitInMask(std::vector<ullint> *va, const size_t pos);
void unsetBitInMask(ushort* va, const size_t pos);
void unsetBitInMask(std::vector<ushort> *va, const size_t pos);
/// \}

/// \brief Retrieve a numbered bit result from an array of unsigned integers.  This function will
///        return 1 or zero depending on the value of the queried bit.
///
/// Overloaded:
///   - Use unsigned integers
///   - Use long long unsigned integers
///   - Use short unsigned integers
///   - Operate on a C-style array trusted to be of sufficient length, or a Standard Template
///     Library vector (the size will be checked)
///
/// \param va   The bit mask to evaluate
/// \param pos  The positon of the bit to access
/// \{
int readBitFromMask(const uint* va, const size_t pos);
int readBitFromMask(const std::vector<uint> &va, const size_t pos);
int readBitFromMask(const ullint* va, const size_t pos);
int readBitFromMask(const std::vector<ullint> &va, const size_t pos);
int readBitFromMask(const ushort* va, const size_t pos);
int readBitFromMask(const std::vector<ushort> &va, const size_t pos);
/// \}

/// \brief Color a vector of boolean data according to the elements of a mask.
///
/// Overloaded:
///   - Have a new vector produced as output, with all elements not named in mask set to false
///   - Supply the pre-allocated result vector (any bits currently marked TRUE will not be wiped)
///
/// \param mask    The indices of the bool vector to mark as TRUE.
/// \param result  The pre-allocated result vector
/// \{
template <typename T> std::vector<bool> colorVectorMask(const std::vector<T> &mask);

template <typename T> std::vector<bool> colorVectorMask(const std::vector<T> &mask, size_t length);

template <typename T> void colorVectorMask(std::vector<bool> *result, const std::vector<T> &mask);
/// \}

/// \brief Construct a vector based on an existing vector of any object type and an ordered list
///        of indices.  The original vector is limited to INT_MAX unique entries.
///
/// Overloaded:
///   - Tile a vector by replicating it a certain number of times
///   - Pluck particular values from a vector, arrange them into a new vector, and replicate that
///     a certain number of times
///
/// \param va    The original vector of objects
/// \param tidx  Vector of indices into va, from which to construct the result
/// \param nrep  The number of replicas to concatenate into the result
/// \{
template <typename T> std::vector<T> tileVector(const std::vector<T> &va, size_t nrep = 1);

template <typename T>
std::vector<T> tileVector(const std::vector<T> &va, const std::vector<int> &tidx, size_t nrep = 1);

template <typename T>
std::vector<T> tileVector(const std::vector<T> &va, const std::vector<size_t> &tidx,
                          size_t nrep = 1);
/// \}

/// \brief Construct a vector of a new type to replace an existing vector.  The new vector will
///        contain all data from the original, scaled by an optional constant to facilitate
///        conversions between real and scaled integral types.  The original vector will be resized
///        and shrunk to zero after transferring its data to the new vector.
///
/// Overloaded:
///   - Convert from an int95_t vector to an arbitrary type
///   - Convert from an arbitrary type to an int95_t vector
///   - Convert between two arbitrary templated types
///   - Work in Standard Template Library vectors or Hybrid objects (only data on the host will be
///     copied)
///
/// \param input         The original data
/// \param input_scale   Scaling factor already applied to the input vector
/// \param output_scale  Optional scaling factor to apply to the output vector 
/// \{
template <typename Tout>
std::vector<Tout> convertData(std::vector<int95_t> *input, double input_scale = 1.0,
                              double output_scale = 1.0);

template <typename Tout>
std::vector<Tout> convertData(std::vector<llint> *primary, std::vector<int> *overflow,
                              double input_scale = 1.0, double output_scale = 1.0);

template <typename Tin>
std::vector<int95_t> convertData(std::vector<Tin> *input, double input_scale = 1.0,
                                 double output_scale = 1.0);

template <typename Tin, typename Tout>
std::vector<Tout> convertData(std::vector<Tin> *input, double input_scale = 1.0,
                              double output_scale = 1.0);

template <typename Tout>
Hybrid<Tout> convertData(Hybrid<int95_t> *input, double input_scale = 1.0,
                         double output_scale = 1.0);

template <typename Tout>
Hybrid<Tout> convertData(Hybrid<llint> *primary, Hybrid<int> *overflow, double input_scale = 1.0,
                         double output_scale = 1.0);

template <typename Tin>
Hybrid<int95_t> convertData(Hybrid<Tin> *input, double input_scale = 1.0,
                            double output_scale = 1.0);

template <typename Tin, typename Tout>
Hybrid<Tout> convertData(Hybrid<Tin> *input, double input_scale = 1.0, double output_scale = 1.0);
/// \}

/// \brief Perform an out-of-place rearrangement of a vector of data based on a sorted pattern.
///        The newly sorted pattern is returned.  This array returns a vector with length equal
///        to the sorting pattern--any rearrangement, even one with missing elements or replication
///        of some elements is permissible.
///
/// Overloaded:
///   - Provide a vector of int or size_t numbers indicating the elements of the original array to
///     transfer, in order, to the new array
///   - Provide a vector of tuples with "x" and "y" members indicating the sorting pattern
///
/// \param unsorted_data    The data set with some or all elements out of their proper place
/// \param sorting_pattern  A sequence of new array positions for the contents of unsorted data,
///                         with the ith element indicating (in its "y" member if this is a
///                         two-tuple with a second template argument), where in the unsorted data
///                         array to obtain the value for the ith position of the output.  If this
///                         is a tuple, its "x" member is irrelevant, and probably was the data
///                         used to perform the sorting and place the "y" members of each tuple in
///                         the proper places.
/// \{
template <typename T>
std::vector<T> applyAssociatedSort(const std::vector<T> &unsorted_data,
                                   const std::vector<int> &sorting_pattern);

template <typename T>
std::vector<T> applyAssociatedSort(const std::vector<T> &unsorted_data,
                                   const std::vector<size_t> &sorting_pattern);

template <typename T, typename T2>
std::vector<T> applyAssociatedSort(const std::vector<T> &unsorted_data,
                                   const std::vector<T2> &sorting_pattern);
/// \}

} // namespace stmath
} // namespace stormm

#include "vector_ops.tpp"

#endif
