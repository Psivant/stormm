// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace testing {

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T> TestSystemManager::exportCoordinateSeries(const int base_system,
                                                              const int frame_count,
                                                              const double perturbation_sigma,
                                                              const int xrs_seed,
                                                              const int scale_bits) const {
  CoordinateSeries<T> result(exportCoordinateFrame(base_system), frame_count, scale_bits);
  CoordinateSeriesWriter<T> resr = result.data();
  Xoshiro256ppGenerator xrs(xrs_seed);
  if (fabs(perturbation_sigma) > constants::tiny) {
    for (int i = 0; i < resr.nframe; i++) {
      for (int j = 0; j < resr.natom; j++) {
        const size_t ij_atom = (static_cast<size_t>(i) *
                                roundUp<size_t>(resr.natom, warp_size_zu)) + j;
        resr.xcrd[ij_atom] = ((static_cast<double>(resr.xcrd[ij_atom]) * resr.inv_gpos_scale) +
                              (perturbation_sigma * xrs.gaussianRandomNumber())) * resr.gpos_scale;
        resr.ycrd[ij_atom] = ((static_cast<double>(resr.ycrd[ij_atom]) * resr.inv_gpos_scale) +
                              (perturbation_sigma * xrs.gaussianRandomNumber())) * resr.gpos_scale;
        resr.zcrd[ij_atom] = ((static_cast<double>(resr.zcrd[ij_atom]) * resr.inv_gpos_scale) +
                              (perturbation_sigma * xrs.gaussianRandomNumber())) * resr.gpos_scale;
      }
    }
  }
  return result;
}

} // namespace testing
} // namespace stormm
