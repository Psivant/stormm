// -*-c++-*-
#include "copyright.h"
#include "Constants/fixed_precision.h"
#include "Reporting/error_format.h"
#include "condensate.h"
#include "hpc_condensate.cuh"

namespace stormm {
namespace synthesis {

using numerics::globalpos_scale_nonoverflow_bits;
using numerics::max_llint_accumulation;

#include "Math/rounding.cui"

//-------------------------------------------------------------------------------------------------
// Transfer data from a PhaseSpaceSynthesis object to the Condensate's coordinate arrays
//
// Arguments:
//   cdw:       Writeable abstract for the Condensate object
//   poly_psr:  Read-only abstract for the source coordinate synthesis
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kCondensateUpdate(CondensateWriter cdw, const PsSynthesisReader poly_psr) {
  const size_t last_sys = poly_psr.system_count - 1;
  const int padded_atoms = poly_psr.atom_starts[last_sys] +
                           devcRoundUp(poly_psr.atom_counts[last_sys], warp_size_int);
  const int istride = gridDim.x * blockDim.x;
  for (int i = threadIdx.x + (blockIdx.x * blockDim.x); i < padded_atoms; i += istride) {
    double xi = poly_psr.xcrd[i];
    double yi = poly_psr.ycrd[i];
    double zi = poly_psr.zcrd[i];
    if (poly_psr.gpos_bits > globalpos_scale_nonoverflow_bits) {
      xi += (double)(poly_psr.xcrd_ovrf[i]) * max_llint_accumulation;
      yi += (double)(poly_psr.ycrd_ovrf[i]) * max_llint_accumulation;
      zi += (double)(poly_psr.zcrd_ovrf[i]) * max_llint_accumulation;
    }
    switch (cdw.mode) {
    case PrecisionModel::DOUBLE:
      cdw.xcrd[i] = xi * poly_psr.inv_gpos_scale;
      cdw.ycrd[i] = yi * poly_psr.inv_gpos_scale;
      cdw.zcrd[i] = zi * poly_psr.inv_gpos_scale;
      break;
    case PrecisionModel::SINGLE:
      cdw.xcrd_sp[i] = (float)(xi) * poly_psr.inv_gpos_scale_f;
      cdw.ycrd_sp[i] = (float)(yi) * poly_psr.inv_gpos_scale_f;
      cdw.zcrd_sp[i] = (float)(zi) * poly_psr.inv_gpos_scale_f;
      break;
    }
  }
  const size_t matrix_entries = (size_t)(poly_psr.system_count) *
                                (size_t)(devcRoundUp(9, warp_size_int));
  for (size_t i = threadIdx.x + (blockIdx.x * blockDim.x); i < matrix_entries; i += istride) {
    if ((i & warp_bits_mask_zu) < 9) {
      cdw.umat[i] = poly_psr.umat[i];
      cdw.invu[i] = poly_psr.invu[i];
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
void Condensate::launchCondensateUpdate(const GpuDetails &gpu) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  if (pps_ptr != nullptr) {
    kCondensateUpdate<<<gpu.getSMPCount(), large_block_size>>>(this->data(tier),
                                                               pps_ptr->data(tier));
  }
  else if (cs_ptr != nullptr) {
    if (csptr_data_type == double_type_index) {
      const CoordinateSeries<double>* dcs_ptr =
        reinterpret_cast<CoordinateSeries<double>*>(cs_ptr);
      kCondensateUpdate<<<gpu.getSMPCount(),large_block_size>>>(this->data(tier),
                                                                dcs_ptr->data(tier));
    }
    else if (csptr_data_type == float_type_index) {
      const CoordinateSeries<float>* fcs_ptr = reinterpret_cast<CoordinateSeries<float>*>(cs_ptr);
      kCondensateUpdate<<<gpu.getSMPCount(),large_block_size>>>(this->data(tier),
                                                                fcs_ptr->data(tier));
    }
    else if (csptr_data_type == llint_type_index) {
      const CoordinateSeries<llint>* lcs_ptr = reinterpret_cast<CoordinateSeries<llint>*>(cs_ptr);
      kCondensateUpdate<<<gpu.getSMPCount(),large_block_size>>>(this->data(tier),
                                                                lcs_ptr->data(tier));
    }
    else if (csptr_data_type == int_type_index) {
      const CoordinateSeries<int>* ics_ptr = reinterpret_cast<CoordinateSeries<int>*>(cs_ptr);
      kCondensateUpdate<<<gpu.getSMPCount(),large_block_size>>>(this->data(tier),
                                                                ics_ptr->data(tier));
    }
    else if (csptr_data_type == short_type_index) {
      const CoordinateSeries<short>* scs_ptr = reinterpret_cast<CoordinateSeries<short>*>(cs_ptr);
      kCondensateUpdate<<<gpu.getSMPCount(),large_block_size>>>(this->data(tier),
                                                                scs_ptr->data(tier));
    }
    else if (csptr_data_type == char_type_index) {
      const CoordinateSeries<char>* ccs_ptr = reinterpret_cast<CoordinateSeries<char>*>(cs_ptr);
      kCondensateUpdate<<<gpu.getSMPCount(),large_block_size>>>(this->data(tier),
                                                                ccs_ptr->data(tier));
    }
    else {
      rtErr("A CoordinateSeries must be typed as double, float, or some signed integer in "
              "order to submit for analysis.  Check the data type, or call the update() function "
              "by directly supplying a pointer to the original CoordinateSeries.", "Condensate",
              "launchCondensaeUpdate");
    }
  }
  else {
    rtErr("There is no current coordinate synthesis or series to base the object upon.",
            "Condensate", "launchCondensateUpdate");
  }
}

} // namespace synthesis
} // namespace stormm
