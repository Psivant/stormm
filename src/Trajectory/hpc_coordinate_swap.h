// -*-c++-*-
#ifndef STORMM_HPC_COORDINATE_SWAP_H
#define STORMM_HPC_COORDINATE_SWAP_H

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "DataTypes/stormm_vector_types.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "coordinate_series.h"

namespace stormm {
namespace trajectory {

using card::GpuDetails;
using synthesis::CondensateWriter;
using synthesis::PsSynthesisWriter;
  
/// \brief Launch the appropriate kernel to swap coordinates (and, if applicable, velocities and
///        forces) between one or more systems in two separate objects.
///
/// \param v_first     Abstract of the first object involved in the swap
/// \param v_second    Abstract of the second object involved in the swap
/// \param frm_first   Index of the system or frame in the first object involved in the swap (its
///                    atom count is assumed to have been validated against the system or frame
///                    indicated by frm_second)
/// \param frm_second  Index of the system in the second object involved in the swap
/// \param frm_pairs   Array of system index pairs to swap.  Any given launchSwapCoordinates
///                    function will either have this array with its length pair_count, or both of
///                    frm_first and frm_second.
/// \param pair_count  The trusted length of frm_pairs
/// \param gpu         Details of the GPU to use in the swap
/// \{
void launchSwapCoordinates(CoordinateSeriesWriter<void> *v_first, size_t frm_first,
                           CoordinateSeriesWriter<void> *v_second, size_t frm_second, size_t ct,
                           const GpuDetails &gpu);

void launchSwapCoordinates(CoordinateSeriesWriter<void> *v_first,
                           CoordinateSeriesWriter<void> *v_second, const int2* frm_pairs,
                           int pair_count, size_t ct, const GpuDetails &gpu);

void launchSwapCoordinates(PsSynthesisWriter *first, size_t frm_first, PsSynthesisWriter *second,
                           size_t frm_second, const GpuDetails &gpu);

void launchSwapCoordinates(PsSynthesisWriter *first, PsSynthesisWriter *second,
                           const int2* frm_pairs, int pair_count, const GpuDetails &gpu);

void launchSwapCoordinates(CondensateWriter *first, size_t frm_first, CondensateWriter *second,
                           size_t frm_second, const GpuDetails &gpu);

void launchSwapCoordinates(CondensateWriter *first, CondensateWriter *second,
                           const int2* frm_pairs, int pair_count, const GpuDetails &gpu);
/// \}

} // namespace trajectory
} // namespace stormm

#endif
