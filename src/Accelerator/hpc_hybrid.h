// -*-c++-*-
#ifndef STORMM_HPC_HYBRID_H
#define STORMM_HPC_HYBRID_H

#include "copyright.h"

namespace stormm {
namespace card {

/// \brief
///
void launchDeepCopy(void* vdest_host, void* vdest_devc, const void* vorig_host,
                    const void* vorig_devc, size_t dest_offset, size_t orig_offset, size_t length,
                    size_t ct, bool do_hdc, bool do_dhc, bool do_ddc, int dest_bits,
                    int orig_bits);

void launchDeepCopy(void* vdest_host, void* vdest_devc, const void* vorig_host,
                    const void* vorig_devc, size_t dest_offset, size_t orig_offset, size_t length,
                    size_t dest_ct, size_t orig_ct, bool do_hdc, bool do_dhc, bool do_ddc,
                    int dest_bits, int orig_bits);
 
} // namespace card
} // namespace stormm

#endif
