// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename T>
void swapValidateAtomCounts(const int frm_first, const T* natom_first, const int frm_second,
                            const T* natom_second) {
  if (natom_first[frm_first] != natom_second[frm_second]) {
    rtErr("The number of atoms in system " + std::to_string(frm_first) + " (" +
          std::to_string(natom_first[frm_first]) + ") does not match the number of atoms in "
          "system " + std::to_string(frm_second) + " (" +
          std::to_string(natom_second[frm_second]) + ").", "swapValidateAtomCounts");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void swapValidateAtomCounts(const int2* frm_pairs, int pair_count, const T* natom_first,
                            const T* natom_second) {
  for (int i = 0; i < pair_count; i++) {
    swapValidateAtomCounts(frm_pairs[i].x, natom_first, frm_pairs[i].y, natom_second);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void swapCoordinates(T* first_x, T* first_y, T* first_z, T* second_x, T* second_y, T* second_z,
                     const size_t first_ofs, const size_t second_ofs, const size_t natom,
                     const size_t frm_first, const size_t frm_second, int* first_xovrf,
                     int* first_yovrf, int* first_zovrf, int* second_xovrf, int* second_yovrf,
                     int* second_zovrf, double* first_umat, double* first_invu, double* first_bdim,
                     double* second_umat, double* second_invu, double* second_bdim) {
  
  // Swap the coordinates
  if (first_xovrf != nullptr && first_yovrf != nullptr && first_zovrf != nullptr &&
      second_xovrf != nullptr && second_yovrf != nullptr && second_zovrf != nullptr) {
    for (size_t i = 0; i < natom; i++) {
      std::swap(first_x[first_ofs + i], second_x[second_ofs + i]);
      std::swap(first_y[first_ofs + i], second_y[second_ofs + i]);
      std::swap(first_z[first_ofs + i], second_z[second_ofs + i]);
      std::swap(first_xovrf[first_ofs + i], second_xovrf[second_ofs + i]);
      std::swap(first_yovrf[first_ofs + i], second_yovrf[second_ofs + i]);
      std::swap(first_zovrf[first_ofs + i], second_zovrf[second_ofs + i]);
    }
  }
  else {
    for (size_t i = 0; i < natom; i++) {
      std::swap(first_x[first_ofs + i], second_x[second_ofs + i]);
      std::swap(first_y[first_ofs + i], second_y[second_ofs + i]);
      std::swap(first_z[first_ofs + i], second_z[second_ofs + i]);
    }
  }

  // Swap the box dimensions
  if (first_umat  != nullptr && first_invu  != nullptr && first_bdim  != nullptr &&
      second_umat != nullptr && second_invu != nullptr && second_bdim != nullptr) {
    const size_t xfrm_stride = roundUp(9, warp_size_int);
    const size_t xfrm_first_ofs  = xfrm_stride * frm_first;
    const size_t xfrm_second_ofs = xfrm_stride * frm_second;
    for (size_t i = 0; i < 9; i++) {
      std::swap(first_umat[xfrm_first_ofs + i], second_umat[xfrm_second_ofs + i]);
      std::swap(first_invu[xfrm_first_ofs + i], second_invu[xfrm_second_ofs + i]);
    }
    const size_t bdim_stride = roundUp(6, warp_size_int);
    const size_t bdim_first_ofs  = bdim_stride * frm_first;
    const size_t bdim_second_ofs = bdim_stride * frm_second;
    for (size_t i = 0; i < 6; i++) {
      std::swap(first_bdim[bdim_first_ofs + i], second_bdim[bdim_second_ofs + i]);
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
template <typename T>
void swapCoordinates(T* first_x, T* first_y, T* first_z, T* second_x, T* second_y, T* second_z,
                     const size_t first_ofs, const size_t second_ofs, const size_t natom,
                     const size_t frm_first, const size_t frm_second, double* first_umat,
                     double* first_invu, double* first_bdim, double* second_umat,
                     double* second_invu, double* second_bdim) {
  swapCoordinates(first_x, first_y, first_z, second_x, second_y, second_z, first_ofs, second_ofs,
                  natom, frm_first, frm_second, nullptr, nullptr, nullptr, nullptr, nullptr,
                  nullptr, first_umat, first_invu, first_bdim, second_umat, second_invu,
                  second_bdim);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void swapCoordinates(CoordinateSeriesWriter<void> *v_first, const size_t frm_first,
                     CoordinateSeriesWriter<void> *v_second, const size_t frm_second) {
  CoordinateSeriesWriter<T> first = restoreType<T>(v_first);
  CoordinateSeriesWriter<T> second = restoreType<T>(v_second);
  const size_t natom = first.natom;
  const size_t padded_natom = roundUp(first.natom, warp_size_int);
  swapCoordinates(first.xcrd, first.ycrd, first.zcrd, second.xcrd, second.ycrd, second.zcrd,
                  padded_natom * frm_first, padded_natom * frm_second, natom, frm_first,
                  frm_second, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, first.umat,
                  first.invu, first.boxdim, second.umat, second.invu, second.boxdim);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void swapCoordinates(CoordinateSeriesWriter<void> *v_first, CoordinateSeriesWriter<void> *v_second,
                     const int2* frm_pairs, const int pair_count) {
  CoordinateSeriesWriter<T> first = restoreType<T>(v_first);
  CoordinateSeriesWriter<T> second = restoreType<T>(v_second);
  const size_t natom = first.natom;
  const size_t padded_natom = roundUp(first.natom, warp_size_int);
  for (int i = 0; i < pair_count; i++) {
    swapCoordinates(first.xcrd, first.ycrd, first.zcrd, second.xcrd, second.ycrd, second.zcrd,
                    padded_natom * static_cast<size_t>(frm_pairs[i].x),
                    padded_natom * static_cast<size_t>(frm_pairs[i].y), natom, frm_pairs[i].x,
                    frm_pairs[i].y, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                    first.umat, first.invu, first.boxdim, second.umat, second.invu, second.boxdim);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordSwap(CoordinateSeries<T> *first, const size_t frm_first, CoordinateSeries<T> *second,
               const size_t frm_second, const HybridTargetLevel tier_first,
               const HybridTargetLevel tier_second, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  switch (tier_first) {
  case HybridTargetLevel::HOST:
    {
      CoordinateSeriesWriter<void> csw_b = second->templateFreeData(tier_second);
      switch (tier_second) {
      case HybridTargetLevel::HOST:
        {
          CoordinateSeriesWriter<void> csw_a = first->templateFreeData(tier_first);
          coordSwap(&csw_a, frm_first, &csw_b, frm_second, ct, tier_first, tier_second, gpu, sync);
        }
        break;
#ifdef STORMM_USE_HPC
      case HybridTargetLevel::DEVICE:
        {
          CoordinateSeriesWriter<void> csw_a = first->deviceViewToTemplateFreeHostData();
          coordSwap(&csw_a, frm_first, &csw_b, frm_second, ct, tier_first, tier_second, gpu, sync);
        }
        break;
#endif
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CoordinateSeriesWriter<void> csw_a = first->templateFreeData(tier_first);
      switch (tier_second) {
      case HybridTargetLevel::HOST:
        {
          CoordinateSeriesWriter<void> csw_b = second->deviceViewToTemplateFreeHostData();
          coordSwap(&csw_a, frm_first, &csw_b, frm_second, ct, tier_first, tier_second, gpu, sync);
        }
        break;
      case HybridTargetLevel::DEVICE:
        {
          CoordinateSeriesWriter<void> csw_b = second->templateFreeData(tier_second);
          coordSwap(&csw_a, frm_first, &csw_b, frm_second, ct, tier_first, tier_second, gpu, sync);
        }
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordSwap(CoordinateSeries<T> *first, CoordinateSeries<T> *second,
               const Hybrid<int2> &frm_pairs, const HybridTargetLevel tier_first,
               const HybridTargetLevel tier_second, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  switch (tier_first) {
  case HybridTargetLevel::HOST:
    {
      CoordinateSeriesWriter<void> csw_b = second->templateFreeData(tier_second);
      switch (tier_second) {
      case HybridTargetLevel::HOST:
        {
          CoordinateSeriesWriter<void> csw_a = first->templateFreeData(tier_first);
          coordSwap(&csw_a, &csw_b, frm_pairs.data(), frm_pairs.size(), ct, tier_first,
                    tier_second, gpu, sync);
        }
        break;
#ifdef STORMM_USE_HPC
      case HybridTargetLevel::DEVICE:
        {
          CoordinateSeriesWriter<void> csw_a = first->deviceViewToTemplateFreeHostData();
          coordSwap(&csw_a, &csw_b, frm_pairs.getDeviceValidHostPointer(), frm_pairs.size(), ct,
                    tier_first, tier_second, gpu, sync);
        }
        break;
#endif
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CoordinateSeriesWriter<void> csw_a = first->templateFreeData(tier_first);
      switch (tier_second) {
      case HybridTargetLevel::HOST:
        {
          CoordinateSeriesWriter<void> csw_b = second->deviceViewToTemplateFreeHostData();
          coordSwap(&csw_a, &csw_b, frm_pairs.getDeviceValidHostPointer(), frm_pairs.size(), ct,
                    tier_first, tier_second, gpu, sync);
        }
        break;
      case HybridTargetLevel::DEVICE:
        {
          CoordinateSeriesWriter<void> csw_b = second->templateFreeData(tier_second);
          coordSwap(&csw_a, &csw_b, frm_pairs.getDeviceValidHostPointer(), frm_pairs.size(), ct,
                    tier_first, tier_second, gpu, sync);
        }
        break;
      }
    }
    break;
#endif
  }
}

} // namespace trajectory
} // namespace stormm
