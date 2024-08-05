// -*-c++-*-
#ifndef STORMM_HPC_HYBRID_UTIL_H
#define STORMM_HPC_HYBRID_UTIL_H

#include "copyright.h"
#include "gpu_details.h"
#include "hybrid.h"

namespace stormm {
namespace card {

/// \brief Launch one of the overloaded, templated kernels to perform deep copies and possibly
///        type conversions.
///
/// \{
void launchDeepCopy(void* vdest_host, void* vdest_devc, const void* vorig_host,
                    const void* vorig_devc, const size_t dest_offset, const size_t orig_offset,
                    const size_t length, const size_t ct, const bool do_hdc, const bool do_dhc,
                    const bool do_ddc, const GpuDetails &gpu, const int dest_bits,
                    const int orig_bits);

void launchDeepCopy(void* vdest_host, void* vdest_devc, const void* vorig_host,
                    const void* vorig_devc, const size_t dest_offset, const size_t orig_offset,
                    const size_t length, const size_t dest_ct, const size_t orig_ct,
                    const bool do_hdc, const bool do_dhc, const bool do_ddc, const GpuDetails &gpu,
                    const int dest_bits, const int orig_bits);
/// \}

} // namespace card
} // namespace stormm

#endif
