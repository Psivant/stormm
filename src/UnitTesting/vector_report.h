// -*-c++-*-
#ifndef STORMM_VECTOR_REPORT_H
#define STORMM_VECTOR_REPORT_H

#include <vector>
#include "copyright.h"
#include "Parsing/polynumeric.h"

namespace stormm {
namespace testing {

using parse::NumberFormat;
using parse::PolyNumeric;

/// \brief Report on the alignment of two vectors.  This should be called once the vectors have
///        been found to be different by some criterion, to check whether shifting them some
///        number of indices relative to one another would lead to a better match.
///
/// \param va           The first vector
/// \param vb           The second vector
/// \param data_format  The numericla format of each vector
/// \param tol          Comparison tolerance for the two vectors
std::string vectorAlignmentReport(const std::vector<PolyNumeric> &va,
                                  const std::vector<PolyNumeric> &vb, NumberFormat data_format,
                                  double tol);

} // namespace testing
} // namespace stormm

#endif
