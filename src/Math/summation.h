// -*-c++-*-
#ifndef STORMM_SUMMATION_H
#define STORMM_SUMMATION_H

#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "UnitTesting/approx.h"

namespace stormm {
namespace stmath {

using card::Hybrid;
using card::HybridTargetLevel;
using errors::rtErr;
using testing::Approx;
using testing::ComparisonType;

/// \brief The type of prefix sum to compute.
enum class PrefixSumType {
  EXCLUSIVE, INCLUSIVE
};

// \brief Compute an in-place prefix sum over a vector.  For "capped" prefix sums, which are
///       exclusive sums followed by one extra entry containing the total sum over the entire
///       array, provide an array that is one element larger than the actual number of indices of
///       interest with a zero in the final element.
///
/// Overloaded:
///   - Operate on a C-style array with a trusted length
///   - Operate on a Standard Template Library vector
///   - Operate on a Hybrid object (HOST data only)
///
/// \param v       The vector to accumulate
/// \param style   Method for computing the sum
/// \param caller  Name of the calling function (optional)
/// \{
template <typename TBase>
void prefixSumInPlace(TBase* vdata, const size_t n_elements, const PrefixSumType style,
                      const char* caller = nullptr);

template <typename TBase>
void prefixSumInPlace(std::vector<TBase> *v, const PrefixSumType style,
                      const char* caller = nullptr);

template <typename TBase>
void prefixSumInPlace(Hybrid<TBase> *v, const PrefixSumType style, const char* caller = nullptr);
/// \}

/// \brief Sum the elements of a standard template library vector for elements of known types.
///        Developers must take care that the sum type have a sufficient range, particularly in
///        cases with integral types.
///
/// Overloaded:
///   - Signed integer sums, all returning long long int
///   - Unsigned integer sums, all returning unsigned long long int
///   - Real number vector sums, all returning a double precision floating point number
///   - Summations for HPC vector types following the above conventions for each component
///
/// \param v     The vector to sum
/// \param hb    The Hybrid object to sum (available only for HOST data in C++ code--see the
///              hpc_summation.cuh header file for kernels implementations to sum data on the
///              DEVICE side)
/// \param vlen  The length of the vector
/// \{
template <typename TSum, typename TBase>
TSum sum(const TBase* v, const size_t vlen);

template <typename TSum, typename TBase>
TSum sum(const std::vector<TBase> &v);

template <typename TSum, typename TBase>
TSum sum(const Hybrid<TBase> &hb);

template <typename TSum, typename TBase>
TSum sumTuple2(const TBase* v, const size_t vlen);

template <typename TSum, typename TBase>
TSum sumTuple2(const std::vector<TBase> &v);

template <typename TSum, typename TBase>
TSum sumTuple2(const Hybrid<TBase> &hb);

template <typename TSum, typename TBase>
TSum sumTuple3(const TBase* v, const size_t vlen);

template <typename TSum, typename TBase>
TSum sumTuple3(const std::vector<TBase> &v);

template <typename TSum, typename TBase>
TSum sumTuple3(const Hybrid<TBase> &hb);

template <typename TSum, typename TBase>
TSum sumTuple4(const TBase* v, const size_t vlen);

template <typename TSum, typename TBase>
TSum sumTuple4(const std::vector<TBase> &v);

template <typename TSum, typename TBase>
TSum sumTuple4(const Hybrid<TBase> &hb);
/// \}

/// \brief Sum two vectors in an element-wise manner.  The result will always take the type of the
///        first input.  Supplying a non-const first parameter will cause the result to be
///        acccumulated in the first parameter.
///
/// Overloaded:
///   - Sum C-style arrays with a trusted length, Standard Template Library Vectors, or Hybrid
///     Objects.  Hybrid object operations are only available on HOST data in C++ code.
///   - Sum two, three, or four arrays at once
///   - Accumulate into a new object (available for vector type only) or an existing array (all
///     types)
///
/// \param va    The first input
/// \param vb    The second input
/// \param vc    Optional third input
/// \param vd    Optional fourth input
/// \param afac  Multiplyiing factor for the first array (for daxpy-like behavior)
/// \param bfac  Multiplyiing factor for the second array (for daxpy-like behavior)
/// \param cfac  Multiplyiing factor for the third array
/// \param dfac  Multiplyiing factor for the fourth array
/// \{
template <typename Ta, typename Tb>
std::vector<Ta> sum(const std::vector<Ta> &va, const std::vector<Tb> &vb, double afac = 1.0,
                    double bfac = 1.0);

template <typename Ta, typename Tb, typename Tc>
std::vector<Ta> sum(const std::vector<Ta> &va, const std::vector<Tb> &vb,
                    const std::vector<Tc> &vc, double afac = 1.0, double bfac = 1.0,
                    double cfac = 1.0);

template <typename Ta, typename Tb, typename Tc, typename Td>
std::vector<Ta> sum(const std::vector<Ta> &va, const std::vector<Tb> &vb,
                    const std::vector<Tc> &vc, const std::vector<Td> &vd, double afac = 1.0,
                    double bfac = 1.0, double cfac = 1.0, double dfac = 1.0);

template <typename Ta, typename Tb>
void sum(Ta* va, const Tb* vb, const size_t length, double afac = 1.0, double bfac = 1.0);

template <typename Ta, typename Tb, typename Tc>
void sum(Ta* va, const Tb* vb, const Tc* vc, const size_t length, double afac = 1.0,
         double bfac = 1.0, double cfac = 1.0);

template <typename Ta, typename Tb, typename Tc, typename Td>
void sum(Ta* va, const Tb* vb, const Tc* vc, const Td* vd, const size_t length, double afac = 1.0,
         double bfac = 1.0, double cfac = 1.0, double dfac = 1.0);

template <typename Ta, typename Tb>
void sum(std::vector<Ta> *va, const std::vector<Tb> &vb, double afac = 1.0, double bfac = 1.0);

template <typename Ta, typename Tb, typename Tc>
void sum(std::vector<Ta> *va, const std::vector<Tb> &vb, const std::vector<Tc> &vc,
         double afac = 1.0, double bfac = 1.0, double cfac = 1.0);

template <typename Ta, typename Tb, typename Tc, typename Td>
void sum(std::vector<Ta> *va, const std::vector<Tb> &vb, const std::vector<Tc> &vc,
         const std::vector<Td> &vd, double afac = 1.0, double bfac = 1.0, double cfac = 1.0,
         double dfac = 1.0);

template <typename Ta, typename Tb>
void sum(Hybrid<Ta> *va, const Hybrid<Tb> &vb, double afac = 1.0, double bfac = 1.0);

template <typename Ta, typename Tb, typename Tc>
void sum(Hybrid<Ta> *va, const Hybrid<Tb> &vb, const Hybrid<Tc> &vc, double afac = 1.0,
         double bfac = 1.0, double cfac = 1.0);

template <typename Ta, typename Tb, typename Tc, typename Td>
void sum(Hybrid<Ta> *va, const Hybrid<Tb> &vb, const Hybrid<Tc> &vc, const Hybrid<Td> &vd,
         double afac = 1.0, double bfac = 1.0, double cfac = 1.0, double dfac = 1.0);
/// \}

} // namespace stmath
} // namespace stormm

#include "summation.tpp"

#endif

