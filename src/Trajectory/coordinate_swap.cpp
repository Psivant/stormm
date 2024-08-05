#include "copyright.h"
#include "coordinate_swap.h"
#include "hpc_coordinate_swap.h"

namespace stormm {
namespace trajectory {

using topology::getEnumerationName;
  
//-------------------------------------------------------------------------------------------------
void swapValidatePrecisionModels(const int first_bits, const int second_bits, const char* desc) {
  if (first_bits != second_bits) {
    rtErr("Fixed precision bits for " + std::string(desc) + " must be identical (" +
          std::to_string(first_bits) + " vs " + std::to_string(second_bits) + ").",
          "swapValidatePrecisionModels");
  }
}
  
//-------------------------------------------------------------------------------------------------
void swapValidateFrameIndices(const int frm_first, const int count_first, const int frm_second,
                              const int count_second) {
  if (frm_first < 0 || frm_first >= count_first) {
    rtErr("Index " + std::to_string(frm_first) + " is invalid for a collection of " +
          std::to_string(count_first) + " coordinate sets.", "swapValidateFrameIndices");
  }
  if (frm_second < 0 || frm_second >= count_second) {
    rtErr("Index " + std::to_string(frm_second) + " is invalid for a collection of " +
          std::to_string(count_second) + " coordinate sets.", "swapValidateFrameIndices");
  }
}

//-------------------------------------------------------------------------------------------------
void swapValidateFrameIndices(const int2* frm_pairs, const int pair_count, const int count_first,
                              const int count_second) {
  for (int i = 0; i < pair_count; i++) {
    if (frm_pairs[i].x < 0 || frm_pairs[i].x >= count_first) {
      rtErr("Index " + std::to_string(frm_pairs[i].x) + " (pair " + std::to_string(i) + ") is "
            "invalid for a collection of " + std::to_string(count_first) + " coordinate sets.",
            "swapValidateFrameIndices");
    }
    if (frm_pairs[i].y < 0 || frm_pairs[i].y >= count_second) {
      rtErr("Index " + std::to_string(frm_pairs[i].y) + " (pair " + std::to_string(i) + ") is "
            "invalid for a collection of " + std::to_string(count_second) + " coordinate sets.",
            "swapValidateFrameIndices");
    }
  }

  // In addition to the validity of each individual frame / system index, one batch of swaps cannot
  // contain the same indices twice.
  for (int i = 0; i < pair_count; i++) {
    bool x_blocked = false;
    bool y_blocked = false;
    for (int j = 0; j < i; j++) {
      x_blocked = (x_blocked || frm_pairs[j].x == frm_pairs[i].x);
      y_blocked = (y_blocked || frm_pairs[j].y == frm_pairs[i].y);
    }
    if (x_blocked) {
      rtErr("Index " + std::to_string(frm_pairs[i].x) + " is included twice in a batch of " +
            std::to_string(pair_count) + " coordinate sets from the first object (out of " +
            std::to_string(count_first) + " sets in all).  This will lead to undefined behavior.",
            "swalValidateFrameIndices");
    }
    if (y_blocked) {
      rtErr("Index " + std::to_string(frm_pairs[i].y) + " is included twice in a batch of " +
            std::to_string(pair_count) + " coordinate sets from the second object (out of " +
            std::to_string(count_second) + " sets in all).  This will lead to undefined behavior.",
            "swalValidateFrameIndices");
    }
  }  
}

//-------------------------------------------------------------------------------------------------
void swapValidateUnitCells(const UnitCellType uca, const UnitCellType ucb) {
  bool problem = false;
  switch (uca) {
  case UnitCellType::NONE:
    switch (ucb) {
    case UnitCellType::NONE:
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      problem = true;
      break;
    }
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    switch (ucb) {
    case UnitCellType::NONE:
      problem = true;
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      break;
    }
    break;
  }
  if (problem) {
    rtErr("Unable to swap coordinates between systems with different unit cell types (" +
          getEnumerationName(uca) + " vs " + getEnumerationName(ucb) + ").",
          "swapValidateUnitCells");
  }
}
  
//-------------------------------------------------------------------------------------------------
void swapCheckSameObject(const bool is_same, const int2* frm_pairs, const int pair_count,
                         const HybridTargetLevel tier_first, const HybridTargetLevel tier_second) {

  // Return immediately if the objects are not the same, or if the swap is occurring between
  // different CPU and GPU memory tiers of the same object (this can also separate the coordinate
  // sets).
  if (is_same == false || tier_first != tier_second) {
    return;
  }
  for (int i = 0; i < pair_count; i++) {

    // It is actually OK for the same system to be swapped with itself, as this just does nothing,
    // but if the same system appears more than once in any other context then this will lead to
    // undefined behavior: a race condition on the GPU or a dependence on the order of the swap
    // list on the CPU.
    bool x_blocked = false;
    bool y_blocked = false;
    const int fpix = frm_pairs[i].x;
    const int fpiy = frm_pairs[i].y;
    for (int j = 0; j < i; j++) {
      const int fpjx = frm_pairs[j].x;
      const int fpjy = frm_pairs[j].y;
      x_blocked = (x_blocked || fpix == fpjx || fpix == fpjy);
      y_blocked = (y_blocked || fpiy == fpjx || fpiy == fpjy);
    }
    if (x_blocked || y_blocked) {
      const int problem_idx = (x_blocked) ? fpix : fpiy;
      rtErr("Index " + std::to_string(problem_idx) + " is included twice in a batch of " +
            std::to_string(pair_count) + " coordinate sets which swap within the same object.  "
            "This will lead to undefined behavior.", "swapCheckSameObject");
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
void swapValidateAtomCounts(const int natom_first, const int natom_second) {
  if (natom_first != natom_second) {
    rtErr("Cannot swap coordinate sets with different atom counts (" +
          std::to_string(natom_first) + ", " + std::to_string(natom_second) + ").",
          "swapValidateAtomCounts");
  }
}

//-------------------------------------------------------------------------------------------------
void coordSwap(CoordinateSeriesWriter<void> *v_first, const size_t frm_first,
               CoordinateSeriesWriter<void> *v_second, const size_t frm_second, const size_t ct,
               const HybridTargetLevel tier_first, const HybridTargetLevel tier_second,
               const GpuDetails &gpu, const HpcKernelSync sync) {

  // Check the fixed precision bit counts for consistency.  If either object is of floating point
  // type, its FP bit count will return as zero.
  swapValidatePrecisionModels(v_first->gpos_bits, v_second->gpos_bits, "particle positions");
  swapValidateFrameIndices(frm_first, v_first->nframe, frm_second, v_second->nframe);
  swapValidateAtomCounts(v_first->natom, v_second->natom);
  swapValidateUnitCells(v_first->unit_cell, v_second->unit_cell);
  switch (tier_first) {
  case HybridTargetLevel::HOST:
    switch (tier_second) {
    case HybridTargetLevel::HOST:

      // Restore the type of the abstracts
      if (ct == double_type_index) {
        swapCoordinates<double>(v_first, frm_first, v_second, frm_second);
      }
      else if (ct == float_type_index) {
        swapCoordinates<float>(v_first, frm_first, v_second, frm_second);
      }
      else if (ct == llint_type_index) {
        swapCoordinates<llint>(v_first, frm_first, v_second, frm_second);
      }
      else if (ct == int_type_index) {
        swapCoordinates<int>(v_first, frm_first, v_second, frm_second);
      }
      else if (ct == short_type_index) {
        swapCoordinates<short int>(v_first, frm_first, v_second, frm_second);
      }
      return;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // Any swap operations involving data on the GPU device, whether in the first or second
  // coordinate series, will fall through the above switch and execute the following code.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, tier_first, tier_second);
  launchSwapCoordinates(v_first, frm_first, v_second, frm_second, ct, gpu);
  launchResolution(sync, tier_first, tier_second);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordSwap(CoordinateSeriesWriter<void> *v_first, CoordinateSeriesWriter<void> *v_second,
               const int2* frm_pairs, const int pair_count, const size_t ct,
               const HybridTargetLevel tier_first, const HybridTargetLevel tier_second,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  swapValidatePrecisionModels(v_first->gpos_bits, v_second->gpos_bits, "particle positions");
  swapValidateAtomCounts(v_first->natom, v_second->natom);
  swapValidateUnitCells(v_first->unit_cell, v_second->unit_cell);
  switch (tier_first) {
  case HybridTargetLevel::HOST:
    switch (tier_second) {
    case HybridTargetLevel::HOST:
      swapValidateFrameIndices(frm_pairs, pair_count, v_first->nframe, v_second->nframe);
      swapCheckSameObject(v_first == v_second, frm_pairs, pair_count, tier_first, tier_second);

      // Restore the type of the abstracts
      if (ct == double_type_index) {
        swapCoordinates<double>(v_first, v_second, frm_pairs, pair_count);
      }
      else if (ct == float_type_index) {
        swapCoordinates<float>(v_first, v_second, frm_pairs, pair_count);
      }
      else if (ct == llint_type_index) {
        swapCoordinates<llint>(v_first, v_second, frm_pairs, pair_count);
      }
      else if (ct == int_type_index) {
        swapCoordinates<int>(v_first, v_second, frm_pairs, pair_count);
      }
      else if (ct == short_type_index) {
        swapCoordinates<short int>(v_first, v_second, frm_pairs, pair_count);
      }
      return;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // Any swap operations involving data on the GPU device, whether in the first or second
  // coordinate series, will fall through the above switch and execute the following code.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, tier_first, tier_second);
  launchSwapCoordinates(v_first, v_second, frm_pairs, pair_count, ct, gpu);
  swapValidateFrameIndices(frm_pairs, pair_count, v_first->nframe, v_second->nframe);
  swapCheckSameObject(v_first == v_second, frm_pairs, pair_count, tier_first, tier_second);
  launchResolution(sync, tier_first, tier_second);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordSwap(PsSynthesisWriter *first, const int frm_first, PsSynthesisWriter *second,
               const int frm_second, const PsSynthesisBorders &first_bdrs,
               const PsSynthesisBorders &second_bdrs, const HybridTargetLevel tier_first,
               const HybridTargetLevel tier_second, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  swapValidatePrecisionModels(first->gpos_bits, second->gpos_bits, "particle positions");
  swapValidatePrecisionModels(first->vel_bits,  second->vel_bits, "velocities");
  swapValidatePrecisionModels(first->frc_bits,  second->frc_bits, "forces");
  swapValidateFrameIndices(frm_first, first_bdrs.system_count, frm_second,
                           second_bdrs.system_count);
  swapValidateAtomCounts(first_bdrs.atom_counts[frm_first], second_bdrs.atom_counts[frm_second]);
  swapValidateUnitCells(first->unit_cell, second->unit_cell);
  switch (tier_first) {
  case HybridTargetLevel::HOST:
    switch (tier_second) {
    case HybridTargetLevel::HOST:
      {
        const size_t natom = first->atom_counts[frm_first];
        swapCoordinates(first->xcrd, first->ycrd, first->zcrd, second->xcrd, second->ycrd,
                        second->zcrd, first->atom_starts[frm_first],
                        second->atom_starts[frm_second], first->atom_counts[frm_first], frm_first,
                        frm_second, first->xcrd_ovrf, first->ycrd_ovrf, first->zcrd_ovrf,
                        second->xcrd_ovrf, second->ycrd_ovrf, second->zcrd_ovrf, first->umat,
                        first->invu, first->boxdims, second->umat, second->invu, second->boxdims);
        swapCoordinates(first->xvel, first->yvel, first->zvel, second->xvel, second->yvel,
                        second->zvel, first->atom_starts[frm_first],
                        second->atom_starts[frm_second], first->atom_counts[frm_first], frm_first,
                        frm_second, first->xvel_ovrf, first->yvel_ovrf, first->zvel_ovrf,
                        second->xvel_ovrf, second->yvel_ovrf, second->zvel_ovrf);
        swapCoordinates(first->xfrc, first->yfrc, first->zfrc, second->xfrc, second->yfrc,
                        second->zfrc, first->atom_starts[frm_first],
                        second->atom_starts[frm_second], first->atom_counts[frm_first], frm_first,
                        frm_second, first->xfrc_ovrf, first->yfrc_ovrf, first->zfrc_ovrf,
                        second->xfrc_ovrf, second->yfrc_ovrf, second->zfrc_ovrf);
      }
      return;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // Any swap operations involving data on the GPU device, whether in the first or second
  // coordinate series, will fall through the above switch and execute the following code.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, tier_first, tier_second);
  launchSwapCoordinates(first, frm_first, second, frm_second, gpu);
  launchResolution(sync, tier_first, tier_second);
#endif
}
               
//-------------------------------------------------------------------------------------------------
void coordSwap(PhaseSpaceSynthesis *first, const int frm_first, PhaseSpaceSynthesis *second,
               const int frm_second, const HybridTargetLevel tier_first,
               const HybridTargetLevel tier_second, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  checkCopyValidity(first, *second, tier_first, tier_second);
  const PsSynthesisBorders first_bdrs  = first->borders();
  const PsSynthesisBorders second_bdrs = second->borders();
  switch (tier_first) {
  case HybridTargetLevel::HOST:
    {
      PsSynthesisWriter csw_b = second->data(tier_second);
      switch (tier_second) {
      case HybridTargetLevel::HOST:
        {
          PsSynthesisWriter csw_a = first->data(tier_first);
          coordSwap(&csw_a, frm_first, &csw_b, frm_second, first_bdrs, second_bdrs, tier_first,
                    tier_second, gpu, sync);
        }
        break;
#ifdef STORMM_USE_HPC
      case HybridTargetLevel::DEVICE:
        {
          PsSynthesisWriter csw_a = first->deviceViewToHostData();
          coordSwap(&csw_a, frm_first, &csw_b, frm_second, first_bdrs, second_bdrs, tier_first,
                    tier_second, gpu, sync);
        }
        break;
#endif
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      PsSynthesisWriter csw_a = first->data(tier_first);
      switch (tier_second) {
      case HybridTargetLevel::HOST:
        {
          PsSynthesisWriter csw_b = second->deviceViewToHostData();
          coordSwap(&csw_a, frm_first, &csw_b, frm_second, first_bdrs, second_bdrs, tier_first,
                    tier_second, gpu, sync);
        }
        break;
      case HybridTargetLevel::DEVICE:
        {
          PsSynthesisWriter csw_b = second->data(tier_second);
          coordSwap(&csw_a, frm_first, &csw_b, frm_second, first_bdrs, second_bdrs, tier_first,
                    tier_second, gpu, sync);
        }
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordSwap(PsSynthesisWriter *first, PsSynthesisWriter *second,
               const int2* frm_pairs, const int pair_count, const PsSynthesisBorders &first_bdrs,
               const PsSynthesisBorders &second_bdrs, const HybridTargetLevel tier_first,
               const HybridTargetLevel tier_second, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  swapValidatePrecisionModels(first->gpos_bits, second->gpos_bits, "particle positions");
  swapValidatePrecisionModels(first->vel_bits,  second->vel_bits, "velocities");
  swapValidatePrecisionModels(first->frc_bits,  second->frc_bits, "forces");
  swapValidateUnitCells(first->unit_cell, second->unit_cell);
  switch (tier_first) {
  case HybridTargetLevel::HOST:
    switch (tier_second) {
    case HybridTargetLevel::HOST:
      swapValidateFrameIndices(frm_pairs, pair_count, first_bdrs.system_count,
                               second_bdrs.system_count);
      swapValidateAtomCounts(frm_pairs, pair_count, first_bdrs.atom_counts,
                             second_bdrs.atom_counts);
      swapCheckSameObject(first == second, frm_pairs, pair_count, tier_first, tier_second);
      for (int i = 0; i < pair_count; i++) {
        swapCoordinates(first->xcrd, first->ycrd, first->zcrd, second->xcrd, second->ycrd,
                        second->zcrd, first->atom_starts[frm_pairs[i].x],
                        second->atom_starts[frm_pairs[i].y], first->atom_counts[frm_pairs[i].x],
                        frm_pairs[i].x, frm_pairs[i].y, first->xcrd_ovrf, first->ycrd_ovrf,
                        first->zcrd_ovrf, second->xcrd_ovrf, second->ycrd_ovrf,
                        second->zcrd_ovrf, first->umat, first->invu, first->boxdims,
                        second->umat, second->invu, second->boxdims);
        swapCoordinates(first->xvel, first->yvel, first->zvel, second->xvel, second->yvel,
                        second->zvel, first->atom_starts[frm_pairs[i].x],
                        second->atom_starts[frm_pairs[i].y], first->atom_counts[frm_pairs[i].x],
                        frm_pairs[i].x, frm_pairs[i].y, first->xvel_ovrf, first->yvel_ovrf,
                        first->zvel_ovrf, second->xvel_ovrf, second->yvel_ovrf,
                        second->zvel_ovrf);
        swapCoordinates(first->xfrc, first->yfrc, first->zfrc, second->xfrc, second->yfrc,
                        second->zfrc, first->atom_starts[frm_pairs[i].x],
                        second->atom_starts[frm_pairs[i].y], first->atom_counts[frm_pairs[i].x],
                        frm_pairs[i].x, frm_pairs[i].y, first->xfrc_ovrf, first->yfrc_ovrf,
                        first->zfrc_ovrf, second->xfrc_ovrf, second->yfrc_ovrf,
                        second->zfrc_ovrf);
      }
      return;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // Any swap operations involving data on the GPU device, whether in the first or second
  // coordinate series, will fall through the above switch and execute the following code.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, tier_first, tier_second);
  launchSwapCoordinates(first, second, frm_pairs, pair_count, gpu);
  swapValidateFrameIndices(frm_pairs, pair_count, first_bdrs.system_count,
                           second_bdrs.system_count);
  swapValidateAtomCounts(frm_pairs, pair_count, first_bdrs.atom_counts, second_bdrs.atom_counts);
  swapCheckSameObject(first == second, frm_pairs, pair_count, tier_first, tier_second);
  launchResolution(sync, tier_first, tier_second);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordSwap(PhaseSpaceSynthesis *first, PhaseSpaceSynthesis *second,
               const Hybrid<int2> &frm_pairs, const HybridTargetLevel tier_first,
               const HybridTargetLevel tier_second, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  checkCopyValidity(first, *second, tier_first, tier_second);
  const PsSynthesisBorders first_bdrs  = first->borders();
  const PsSynthesisBorders second_bdrs = second->borders();
  switch (tier_first) {
  case HybridTargetLevel::HOST:
    {
      PsSynthesisWriter csw_b = second->data(tier_second);
      switch (tier_second) {
      case HybridTargetLevel::HOST:
        {
          PsSynthesisWriter csw_a = first->data(tier_first);
          coordSwap(&csw_a, &csw_b, frm_pairs.data(), frm_pairs.size(), first_bdrs, second_bdrs,
                    tier_first, tier_second, gpu, sync);
        }
        break;
#ifdef STORMM_USE_HPC
      case HybridTargetLevel::DEVICE:
        {
          PsSynthesisWriter csw_a = first->deviceViewToHostData();
          coordSwap(&csw_a, &csw_b, frm_pairs.getDeviceValidHostPointer(), frm_pairs.size(),
                    first_bdrs, second_bdrs, tier_first, tier_second, gpu, sync);
        }
        break;
#endif
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      PsSynthesisWriter csw_a = first->data(tier_first);
      switch (tier_second) {
      case HybridTargetLevel::HOST:
        {
          PsSynthesisWriter csw_b = second->deviceViewToHostData();
          coordSwap(&csw_a, &csw_b, frm_pairs.data(), frm_pairs.size(), first_bdrs, second_bdrs,
                    tier_first, tier_second, gpu, sync);
        }
        break;
      case HybridTargetLevel::DEVICE:
        {
          PsSynthesisWriter csw_b = second->data(tier_second);
          coordSwap(&csw_a, &csw_b, frm_pairs.getDeviceValidHostPointer(), frm_pairs.size(),
                    first_bdrs, second_bdrs, tier_first, tier_second, gpu, sync);
        }
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordSwap(CondensateWriter *first, const size_t frm_first, CondensateWriter *second,
               const size_t frm_second, const CondensateBorders &first_bdrs,
               const CondensateBorders &second_bdrs, const HybridTargetLevel tier_first,
               const HybridTargetLevel tier_second, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  if (first->mode != second->mode) {
    rtErr("Coordinates can only be swapped between Condensates with the same precision model (" +
          getEnumerationName(first->mode) + ", " + getEnumerationName(second->mode) + ").",
          "coordSwap");
  }
  swapValidateFrameIndices(frm_first, first_bdrs.system_count, frm_second,
                           second_bdrs.system_count);
  swapValidateAtomCounts(first_bdrs.atom_counts[frm_first], second_bdrs.atom_counts[frm_second]);
  swapValidateUnitCells(first->unit_cell, second->unit_cell);
  switch (tier_first) {
  case HybridTargetLevel::HOST:
    switch (tier_second) {
    case HybridTargetLevel::HOST:
      swapCoordinates(first->xcrd, first->ycrd, first->zcrd, second->xcrd, second->ycrd,
                      second->zcrd, first_bdrs.atom_starts[frm_first],
                      second_bdrs.atom_starts[frm_second], first_bdrs.atom_counts[frm_first],
                      frm_first, frm_second, first->umat, first->invu, first->boxdims,
                      second->umat, second->invu, second->boxdims);
      return;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }
  
  // Any swap operations involving data on the GPU device, whether in the first or second
  // coordinate series, will fall through the above switch and execute the following code.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, tier_first, tier_second);
  launchSwapCoordinates(first, frm_first, second, frm_second, gpu);
  launchResolution(sync, tier_first, tier_second);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordSwap(Condensate *first, const size_t frm_first, Condensate *second,
               const size_t frm_second, const HybridTargetLevel tier_first,
               const HybridTargetLevel tier_second, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  checkCopyValidity(first, *second, tier_first, tier_second);
  const CondensateBorders first_bdrs  = first->borders();
  const CondensateBorders second_bdrs = second->borders();
  switch (tier_first) {
  case HybridTargetLevel::HOST:
    {
      CondensateWriter csw_b = second->data(tier_second);
      switch (tier_second) {
      case HybridTargetLevel::HOST:
        {
          CondensateWriter csw_a = first->data(tier_first);
          coordSwap(&csw_a, frm_first, &csw_b, frm_second, first_bdrs, second_bdrs, tier_first,
                    tier_second, gpu, sync);
        }
        break;
#ifdef STORMM_USE_HPC
      case HybridTargetLevel::DEVICE:
        {
          CondensateWriter csw_a = first->deviceViewToHostData();
          coordSwap(&csw_a, frm_first, &csw_b, frm_second, first_bdrs, second_bdrs, tier_first,
                    tier_second, gpu, sync);
        }
        break;
#endif
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CondensateWriter csw_a = first->data(tier_first);
      switch (tier_second) {
      case HybridTargetLevel::HOST:
        {
          CondensateWriter csw_b = second->deviceViewToHostData();
          coordSwap(&csw_a, frm_first, &csw_b, frm_second, first_bdrs, second_bdrs, tier_first,
                    tier_second, gpu, sync);
        }
        break;
      case HybridTargetLevel::DEVICE:
        {
          CondensateWriter csw_b = second->data(tier_second);
          coordSwap(&csw_a, frm_first, &csw_b, frm_second, first_bdrs, second_bdrs, tier_first,
                    tier_second, gpu, sync);
        }
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void coordSwap(CondensateWriter *first, CondensateWriter *second,
               const int2* frm_pairs, const int pair_count, const CondensateBorders &first_bdrs,
               const CondensateBorders &second_bdrs, const HybridTargetLevel tier_first,
               const HybridTargetLevel tier_second, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  if (first->mode != second->mode) {
    rtErr("Coordinates can only be swapped between Condensates with the same precision model (" +
          getEnumerationName(first->mode) + ", " + getEnumerationName(second->mode) + ").",
          "coordSwap");
  }
  swapValidateUnitCells(first->unit_cell, second->unit_cell);
  switch (tier_first) {
  case HybridTargetLevel::HOST:
    switch (tier_second) {
    case HybridTargetLevel::HOST:
      swapValidateFrameIndices(frm_pairs, pair_count, first_bdrs.system_count,
                               second_bdrs.system_count);
      swapValidateAtomCounts(frm_pairs, pair_count, first_bdrs.atom_counts,
                             second_bdrs.atom_counts);
      swapCheckSameObject(first == second, frm_pairs, pair_count, tier_first, tier_second);
      for (int i = 0; i < pair_count; i++) {
        swapCoordinates(first->xcrd, first->ycrd, first->zcrd, second->xcrd, second->ycrd,
                        second->zcrd, first_bdrs.atom_starts[frm_pairs[i].x],
                        second_bdrs.atom_starts[frm_pairs[i].y],
                        first_bdrs.atom_counts[frm_pairs[i].x], frm_pairs[i].x, frm_pairs[i].y,
                        first->umat, first->invu, first->boxdims, second->umat, second->invu,
                        second->boxdims);
      }
      return;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // Any swap operations involving data on the GPU device, whether in the first or second
  // coordinate series, will fall through the above switch and execute the following code.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, tier_first, tier_second);
  launchSwapCoordinates(first, second, frm_pairs, pair_count, gpu);
  swapValidateFrameIndices(frm_pairs, pair_count, first_bdrs.system_count,
                           second_bdrs.system_count);
  swapValidateAtomCounts(frm_pairs, pair_count, first_bdrs.atom_counts, second_bdrs.atom_counts);
  swapCheckSameObject(first == second, frm_pairs, pair_count, tier_first, tier_second);
  launchResolution(sync, tier_first, tier_second);
#endif
}

//-------------------------------------------------------------------------------------------------
void coordSwap(Condensate *first, Condensate *second,
               const Hybrid<int2> &frm_pairs, const HybridTargetLevel tier_first,
               const HybridTargetLevel tier_second, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  checkCopyValidity(first, *second, tier_first, tier_second);
  const CondensateBorders first_bdrs  = first->borders();
  const CondensateBorders second_bdrs = second->borders();
  switch (tier_first) {
  case HybridTargetLevel::HOST:
    {
      CondensateWriter csw_b = second->data(tier_second);
      switch (tier_second) {
      case HybridTargetLevel::HOST:
        {
          CondensateWriter csw_a = first->data(tier_first);
          coordSwap(&csw_a, &csw_b, frm_pairs.data(), frm_pairs.size(), first_bdrs, second_bdrs,
                    tier_first, tier_second, gpu, sync);
        }
        break;
#ifdef STORMM_USE_HPC
      case HybridTargetLevel::DEVICE:
        {
          CondensateWriter csw_a = first->deviceViewToHostData();
          coordSwap(&csw_a, &csw_b, frm_pairs.getDeviceValidHostPointer(), frm_pairs.size(),
                    first_bdrs, second_bdrs, tier_first, tier_second, gpu, sync);
        }
        break;
#endif
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CondensateWriter csw_a = first->data(tier_first);
      switch (tier_second) {
      case HybridTargetLevel::HOST:
        {
          CondensateWriter csw_b = second->deviceViewToHostData();
          coordSwap(&csw_a, &csw_b, frm_pairs.data(), frm_pairs.size(), first_bdrs, second_bdrs,
                    tier_first, tier_second, gpu, sync);
        }
        break;
      case HybridTargetLevel::DEVICE:
        {
          CondensateWriter csw_b = second->data(tier_second);
          coordSwap(&csw_a, &csw_b, frm_pairs.getDeviceValidHostPointer(), frm_pairs.size(),
                    first_bdrs, second_bdrs, tier_first, tier_second, gpu, sync);
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
