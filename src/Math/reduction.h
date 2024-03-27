// -*-c++-*-
#ifndef STORMM_REDUCTION_H
#define STORMM_REDUCTION_H

#include "copyright.h"
#include "DataTypes/common_types.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Numerics/split_fixed_precision.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "reduction_abstracts.h"
#include "reduction_enumerators.h"
#include "reduction_workunit.h"

namespace stormm {
namespace stmath {

using numerics::max_llint_accumulation;

/// \brief Perform gathering operations for a normalization reduction.
///
/// Overloaded:
///   - Accept long long integer data and determine whether there is an extended format
///   - Accept standard data (no extended format will be considered)
///
/// \param rsbs       Data to reduce (passed by value, containing pointers to modifiable values)
/// \param start_pos  Index of the first atom for the work unit at hand to process
/// \param end_pos    Upper limit of atoms for the work unit at hand to process
/// \{
double gatherNormalization(const GenericRdSubstrate<llint> rsbs, const int start_pos,
                           const int end_pos);

template <typename T>
double gatherNormalization(const GenericRdSubstrate<T> rsbs, int start_pos, int end_pos);
/// \}
  
/// \brief Perform gathering operations for a centering-on-zero reduction.
///
/// Overloaded:
///   - Accept long long integer data and determine whether there is an extended format
///   - Accept standard data (no extended format will be considered)
///
/// \param rsbs       Data to reduce (passed by value, containing pointers to modifiable values)
/// \param start_pos  Index of the first atom for the work unit at hand to process
/// \param end_pos    Upper limit of atoms for the work unit at hand to process
/// \{
double3 gatherCenterOnZero(const GenericRdSubstrate<llint> rsbs, int start_pos, int end_pos);

template <typename T>
double3 gatherCenterOnZero(const GenericRdSubstrate<T> rsbs, int start_pos, int end_pos);
/// \}

/// \brief Perform scattering operations for a normalization reduction: divide the values found
///        in the writeable arrays by the magnitude of the vector represented by all values in
///        the readable arrays (computed gatherNormalization).
///
/// Overloaded:
///   - Accept long long integer data and determine whether there is an extended format
///   - Accept standard data (no extended format will be considered)
///
/// \param rsbs       Data to reduce (passed by value, containing pointers to modifiable values)
/// \param tsum       Sum of the vector components' squares (presumably the sum of the squares of
///                   values about to be divided by its square root)
/// \param start_pos  Index of the first atom for the work unit at hand to process
/// \param end_pos    Upper limit of atoms for the work unit at hand to process
/// \{
void scatterNormalization(GenericRdSubstrate<llint> rsbs, double tsum, int start_pos,
                          int end_pos);

template <typename T>
void scatterNormalization(GenericRdSubstrate<T> rsbs, double tsum, int start_pos, int end_pos);
/// \}

/// \brief Perform scattering operations for a center-on-zero reduction of data in up to three
///        dimensions.
///
/// Overloaded:
///   - Accept long long integer data and determine whether there is an extended format
///   - Accept standard data (no extended format will be considered)
///
/// \param rsbs       Data to reduce (passed by value, containing pointers to modifiable values)
/// \param tsum_x     Sum of the vector components in the 1st dimension
/// \param tsum_y     Sum of the vector components in the 2nd dimension
/// \param tsum_z     Sum of the vector components in the 3rd dimension
/// \param natom      Total number of atoms in the system (for taking the system-wide average)
/// \param start_pos  Index of the first atom for the work unit at hand to process
/// \param end_pos    Upper limit of atoms for the work unit at hand to process
/// \{
void scatterCenterOnZero(GenericRdSubstrate<llint> rsbs, double tsum_x, double tsum_y,
                         double tsum_z, int natom, int start_pos, int end_pos);

template <typename T>
void scatterCenterOnZero(GenericRdSubstrate<T> rsbs, double tsum_x, double tsum_y, double tsum_z,
                         int natom, int start_pos, int end_pos);
/// \}

/// \brief Evaluate a reduction.  In CPU code, this takes advantage of the fact that the
///        ReductionKit contains information about whether storage of the gathered results is
///        necessary and will recursively call itself to accomplish the all-reduce.  By default,
///        this function will perform normalizations across each system in all provided data
///        vectors.
///
/// \param rsbs     Data to reduce
/// \param redk     Instructions for the reduction operations to perform
/// \param process  The stage of reduction to perform (if an all-reduce is requested but must be
///                 carried out in stages, the function can call itself recursively)
template <typename T>
void evalReduction(GenericRdSubstrate<T> *rsbs, const ReductionKit &redk,
                   const ReductionStage process = ReductionStage::ALL_REDUCE,
                   const ReductionGoal purpose = ReductionGoal::NORMALIZE);

} // namespace stmath
} // namespace stormm

#include "reduction.tpp"

#endif
