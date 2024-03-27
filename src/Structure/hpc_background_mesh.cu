// -*-c++-*-
#include "copyright.h"
#include "Accelerator/cuda_wrappers.h"
#include "hpc_mesh_support.h"
#include "hpc_background_mesh.cuh"

namespace stormm {
namespace structure {

using card::wrapCudaFuncGetAttributes;
  
#define COLOR_EXCLUSION_THREAD_COUNT  256
#define TCALC double
#  define MESH_BLOCKS_MULTIPLIER 2
#  define PERIODIC_MESH
#    define KERNEL_NAME kdColorOcclusionMeshPbc
#      include "color_occlusion_mesh.cui"
#    undef KERNEL_NAME
#  undef PERIODIC_MESH
#  define KERNEL_NAME kdColorOcclusionMesh
#    include "color_occlusion_mesh.cui"
#  undef KERNEL_NAME
#  undef MESH_BLOCKS_MULTIPLIER
#undef TCALC

#define TCALC float
#  define TCALC_IS_SINGLE
#  define MESH_BLOCKS_MULTIPLIER 4
#  define PERIODIC_MESH
#    define KERNEL_NAME kfColorOcclusionMeshPbc
#      include "color_occlusion_mesh.cui"
#    undef KERNEL_NAME
#  undef PERIODIC_MESH
#  define KERNEL_NAME kfColorOcclusionMesh
#    include "color_occlusion_mesh.cui"
#  undef KERNEL_NAME
#  undef MESH_BLOCKS_MULTIPLIER
#  undef TCALC_IS_SINGLE
#undef TCALC
#undef COLOR_EXCLUSION_THREAD_COUNT

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes queryGridKernelRequirements(const PrecisionModel prec,
                                                      const GridDetail picture,
                                                      const BoundaryCondition bounds,
                                                      const MappingActivity process) {
  cudaFuncAttributes result;
  bool picture_problem = false;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    switch (picture) {
    case GridDetail::OCCLUSION:
      switch (process) {
      case MappingActivity::PARTICLE_TO_MESH:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          if (cudaFuncGetAttributes(&result, kdColorOcclusionMesh) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kdColorOcclusionMesh.",
                  "queryGridKernelRequirements");
          }
          break;
        case BoundaryCondition::PERIODIC:
          if (cudaFuncGetAttributes(&result, kdColorOcclusionMeshPbc) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kdColorOcclusionMeshPbc.",
                  "queryGridKernelRequirements");
          }
          break;
        }
        break;
      case MappingActivity::MESH_TO_PARTICLE:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          break;
        case BoundaryCondition::PERIODIC:
          break;
        }
        break;
      }
      break;
    case GridDetail::OCCLUSION_FIELD:
    case GridDetail::NONBONDED_FIELD:
    case GridDetail::NONBONDED_ATOMIC:
      picture_problem = true;
      break;
    }
  case PrecisionModel::SINGLE:
    switch (picture) {
    case GridDetail::OCCLUSION:
      switch (process) {
      case MappingActivity::PARTICLE_TO_MESH:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          if (cudaFuncGetAttributes(&result, kfColorOcclusionMesh) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kdColorOcclusionMesh.",
                  "queryGridKernelRequirements");
          }
          break;
        case BoundaryCondition::PERIODIC:
          if (cudaFuncGetAttributes(&result, kfColorOcclusionMeshPbc) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel kfColorOcclusionMeshPbc.",
                  "queryGridKernelRequirements");
          }
          break;
        }
        break;
      case MappingActivity::MESH_TO_PARTICLE:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          break;
        case BoundaryCondition::PERIODIC:
          break;
        }
        break;
      }
      break;
    case GridDetail::OCCLUSION_FIELD:
    case GridDetail::NONBONDED_FIELD:
    case GridDetail::NONBONDED_ATOMIC:
      picture_problem = true;
      break;
    }
    break;
  }
  if (picture_problem) {
    rtErr("\"" + getEnumerationName(GridDetail::OCCLUSION) + "\" mesh kernels are agnostic to the "
          "unit cell type, and do not use an interpolant of any kind.  These properties must be "
          "specified to get the requirements for a " + getEnumerationName(picture) + " mesh "
          "kernel.", "queryGridKernelRequirements");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes queryGridKernelRequirements(const PrecisionModel prec,
                                                      const GridDetail picture,
                                                      const BoundaryCondition bounds,
                                                      const UnitCellType unit_cell,
                                                      const Interpolant stencil_kind,
                                                      const MappingActivity process) {
  cudaFuncAttributes result;
  bool picture_problem = false;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    switch (picture) {
    case GridDetail::OCCLUSION:
      picture_problem = true;
      break;
    case GridDetail::OCCLUSION_FIELD:
      switch (process) {
      case MappingActivity::PARTICLE_TO_MESH:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          break;
        case BoundaryCondition::PERIODIC:
          break;
        }
        break;
      case MappingActivity::MESH_TO_PARTICLE:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          break;
        case BoundaryCondition::PERIODIC:
          break;
        }
        break;
      }
      break;
    case GridDetail::NONBONDED_FIELD:
      switch (process) {
      case MappingActivity::PARTICLE_TO_MESH:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          switch (unit_cell) {
          case UnitCellType::ORTHORHOMBIC:
            switch (stencil_kind) {
            case Interpolant::SMOOTHNESS:
              if (wrapCudaFuncGetAttributes(&result, kdColorNBFieldMesh<double>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kdColorNBFieldMesh.",
                      "queryGridKernelRequirements");
              }
              break;
            case Interpolant::FUNCTION_VALUE:
              if (wrapCudaFuncGetAttributes(&result, kdColorNBFieldMeshValue<double>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kdColorNBFieldMeshValue.",
                      "queryGridKernelRequirements");
              }
              break;
            }
            break;
          case UnitCellType::TRICLINIC:
            switch (stencil_kind) {
            case Interpolant::SMOOTHNESS:
              if (wrapCudaFuncGetAttributes(&result, kdColorNBFieldMeshTric<double>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kdColorNBFieldMeshTric.",
                      "queryGridKernelRequirements");
              }
              break;
            case Interpolant::FUNCTION_VALUE:
              if (wrapCudaFuncGetAttributes(&result, kdColorNBFieldMeshTricValue<double>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kdColorNBFieldMeshTricValue.",
                      "queryGridKernelRequirements");
              }
              break;
            }
            break;
          case UnitCellType::NONE:
            break;
          }
          break;
        case BoundaryCondition::PERIODIC:
          switch (unit_cell) {
          case UnitCellType::ORTHORHOMBIC:
            switch (stencil_kind) {
            case Interpolant::SMOOTHNESS:
              if (wrapCudaFuncGetAttributes(&result, kdColorNBFieldMeshPbc<double>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kdColorNBFieldMeshPbc.",
                      "queryGridKernelRequirements");
              }
              break;
            case Interpolant::FUNCTION_VALUE:
              if (wrapCudaFuncGetAttributes(&result, kdColorNBFieldMeshPbcValue<double>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kdColorNBFieldMeshPbcValue.",
                      "queryGridKernelRequirements");
              }
              break;
            }
            break;
          case UnitCellType::TRICLINIC:
            switch (stencil_kind) {
            case Interpolant::SMOOTHNESS:
              if (wrapCudaFuncGetAttributes(&result, kdColorNBFieldMeshPbcTric<double>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kdColorNBFieldMeshPbcTric.",
                      "queryGridKernelRequirements");
              }
              break;
            case Interpolant::FUNCTION_VALUE:
              if (wrapCudaFuncGetAttributes(&result, kdColorNBFieldMeshPbcTricValue<double>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kdColorNBFieldMeshPbcTricValue.",
                      "queryGridKernelRequirements");
              }
              break;
            }
            break;
          case UnitCellType::NONE:
            break;
          }
          break;
        }
        break;
      case MappingActivity::MESH_TO_PARTICLE:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          break;
        case BoundaryCondition::PERIODIC:
          break;
        }
        break;
      }
      break;
    case GridDetail::NONBONDED_ATOMIC:
      switch (process) {
      case MappingActivity::PARTICLE_TO_MESH:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          break;
        case BoundaryCondition::PERIODIC:
          break;
        }
        break;
      case MappingActivity::MESH_TO_PARTICLE:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          break;
        case BoundaryCondition::PERIODIC:
          break;
        }
        break;
      }
      break;
    }
    break;
  case PrecisionModel::SINGLE:
    switch (picture) {
    case GridDetail::OCCLUSION:
      picture_problem = true;
      break;
    case GridDetail::OCCLUSION_FIELD:
      switch (process) {
      case MappingActivity::PARTICLE_TO_MESH:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          break;
        case BoundaryCondition::PERIODIC:
          break;
        }
        break;
      case MappingActivity::MESH_TO_PARTICLE:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          break;
        case BoundaryCondition::PERIODIC:
          break;
        }
        break;
      }
      break;
    case GridDetail::NONBONDED_FIELD:
      switch (process) {
      case MappingActivity::PARTICLE_TO_MESH:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          switch (unit_cell) {
          case UnitCellType::ORTHORHOMBIC:
            switch (stencil_kind) {
            case Interpolant::SMOOTHNESS:
              if (wrapCudaFuncGetAttributes(&result, kfColorNBFieldMesh<float>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kfColorNBFieldMesh.",
                      "queryGridKernelRequirements");
              }
              break;
            case Interpolant::FUNCTION_VALUE:
              if (wrapCudaFuncGetAttributes(&result, kfColorNBFieldMeshValue<float>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kfColorNBFieldMeshValue.",
                      "queryGridKernelRequirements");
              }
              break;
            }
            break;
          case UnitCellType::TRICLINIC:
            switch (stencil_kind) {
            case Interpolant::SMOOTHNESS:
              if (wrapCudaFuncGetAttributes(&result, kfColorNBFieldMeshTric<float>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kfColorNBFieldMeshTric.",
                      "queryGridKernelRequirements");
              }
              break;
            case Interpolant::FUNCTION_VALUE:
              if (wrapCudaFuncGetAttributes(&result, kfColorNBFieldMeshTricValue<float>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kfColorNBFieldMeshTricValue.",
                      "queryGridKernelRequirements");
              }
              break;
            }
            break;
          case UnitCellType::NONE:
            break;
          }
          break;
        case BoundaryCondition::PERIODIC:
          switch (unit_cell) {
          case UnitCellType::ORTHORHOMBIC:
            switch (stencil_kind) {
            case Interpolant::SMOOTHNESS:
              if (wrapCudaFuncGetAttributes(&result, kfColorNBFieldMeshPbc<float>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kfColorNBFieldMeshPbc.",
                      "queryGridKernelRequirements");
              }
              break;
            case Interpolant::FUNCTION_VALUE:
              if (wrapCudaFuncGetAttributes(&result, kfColorNBFieldMeshPbcValue<float>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kfColorNBFieldMeshPbcValue.",
                      "queryGridKernelRequirements");
              }
              break;
            }
            break;
          case UnitCellType::TRICLINIC:
            switch (stencil_kind) {
            case Interpolant::SMOOTHNESS:
              if (wrapCudaFuncGetAttributes(&result, kfColorNBFieldMeshPbcTric<float>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kfColorNBFieldMeshPbcTric.",
                      "queryGridKernelRequirements");
              }
              break;
            case Interpolant::FUNCTION_VALUE:
              if (wrapCudaFuncGetAttributes(&result, kfColorNBFieldMeshPbcTricValue<float>) !=
                  cudaSuccess) {
                rtErr("Error obtaining attributes for kernel kfColorNBFieldMeshPbcTricValue.",
                      "queryGridKernelRequirements");
              }
              break;
            }
            break;
          case UnitCellType::NONE:
            break;
          }
          break;
        }
        break;
      case MappingActivity::MESH_TO_PARTICLE:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          break;
        case BoundaryCondition::PERIODIC:
          break;
        }
        break;
      }
      break;
    case GridDetail::NONBONDED_ATOMIC:
      switch (process) {
      case MappingActivity::PARTICLE_TO_MESH:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          break;
        case BoundaryCondition::PERIODIC:
          break;
        }
        break;
      case MappingActivity::MESH_TO_PARTICLE:
        switch (bounds) {
        case BoundaryCondition::ISOLATED:
          break;
        case BoundaryCondition::PERIODIC:
          break;
        }
        break;
      }
      break;
    }
    break;
  }
  if (picture_problem) {
    rtErr("The unit cell type and interpolant must be specified for any mesh other than an \"" +
          getEnumerationName(GridDetail::OCCLUSION) + "\" mesh, including \"" +
          getEnumerationName(picture) + "\".", "queryGridKernelRequirements");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
extern void launchColorOcclusionMesh(BackgroundMeshWriter<void> *bgmw, const AtomGraph *ag,
                                     const CoordinateFrameReader &cfr, const PrecisionModel prec,
                                     const MeshKlManager &launcher) {
  const int2 nb = launcher.getMeshKernelDims(prec, GridDetail::OCCLUSION, bgmw->dims.bounds,
                                             MappingActivity::PARTICLE_TO_MESH);
  BackgroundMeshWriter<ullint> ll_bgmw = restoreType<ullint>(bgmw);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    kdColorOcclusionMesh<<<nb.x, nb.y>>>(ll_bgmw, cfr, ag->getDoublePrecisionNonbondedKit());
    break;
  case PrecisionModel::SINGLE:
    kfColorOcclusionMesh<<<nb.x, nb.y>>>(ll_bgmw, cfr, ag->getSinglePrecisionNonbondedKit());
    break;
  }
}

//-------------------------------------------------------------------------------------------------
// Accumulate an occlusion field mesh based on the occlusion mesh of a single snapshot.  Abstracts
// of each mesh must be provided in a manner that is accessible to the GPU; in fact, the kernel
// will only be engaged if the relevant data is already present on the GPU device.
//
// Arguments:
//   target:        The mesh being accumulated.  This must be a long long int-type mesh in order
//                  enable fixed-precision accumulation.
//   contribution:  The mesh to contribute.  Dimensions will be checked against those of the   
//                  accumulating target.
//-------------------------------------------------------------------------------------------------
__global__ __launch_bounds__(large_block_size, 1)
void kAccumulateOcclusionMesh(BackgroundMeshWriter<llint> target,
                              const BackgroundMeshReader<llint> &contribution) {
  int period_factor;
  switch (target.dims.bounds) { 
  case BoundaryCondition::ISOLATED:
    period_factor = 1;
    break;
  case BoundaryCondition::PERIODIC:
    period_factor = 0;
    break;
  }
  const int samp_a = 4 * contribution.dims.na / (target.dims.na + period_factor);
  const int samp_b = 4 * contribution.dims.nb / (target.dims.nb + period_factor);
  const int samp_c = 4 * contribution.dims.nc / (target.dims.nc + period_factor);

  // Loop over all elements of the target mesh, dividing them up one per warp.  Each warp will
  // access elements of the contributing  relevant to the
  int wpos = (threadIdx.x >> warp_bits);
  while (wpos < target.dims.na * target.dims.nb * target.dims.nc) {
    const int trg_k = wpos / (target.dims.na * target.dims.nb);
    const int trg_j = (wpos - (trg_k * target.dims.na * target.dims.nb)) / target.dims.na;
    const int trg_i = wpos - (((trg_k * target.dims.nb) + trg_j) * target.dims.na);

    // Loop over the limits of the relevant contribution
    int nbits = 0;
    int tpos = (threadIdx.x & warp_bits_mask_int);
    while (tpos < samp_a * samp_b * samp_c) {      
      const int loc_ci = tpos / (samp_a * samp_b);
      const int loc_cj = (tpos - (loc_ci * samp_a * samp_b)) / samp_a;
      const int loc_ck = tpos - (((loc_ci * samp_b) + loc_cj) * samp_a);
      const int ci = loc_ci + (samp_a * trg_i);
      int ci_cube = ci / 4;
      const int ci_cubelet = ci - (4 * ci_cube);
      if (ci_cube >= contribution.dims.na) {
        ci_cube -= contribution.dims.na;
      }
      const int cj = loc_cj + (samp_b * trg_j);
      int cj_cube = cj / 4;
      const int cj_cubelet = cj - (4 * cj_cube);
      if (cj_cube >= contribution.dims.nb) {
        cj_cube -= contribution.dims.nb;
      }
      const int ck = loc_ck + (samp_c * trg_k);
      int ck_cube = ck / 4;
      const int ck_cubelet = ck - (4 * ck_cube);
      if (ck_cube >= contribution.dims.nc) {
        ck_cube -= contribution.dims.nc;
      }
      const size_t ocidx = ((size_t)(64) *
                            (size_t)((((ck_cube * contribution.dims.nb) +
                                       cj_cube) * contribution.dims.na) + ci_cube)) +
                           (size_t)((((ck_cubelet * 4) + cj_cubelet) * 4) + ci_cubelet);
      nbits += __popcll(contribution.coeffs[ocidx]);
      tpos += warp_size_int;
    }

    // Reduce the bits over all threads.  The first eight threads will log the result in the
    // value slots of up to eight elements in the target mesh.
    WARP_REDUCE_DOWN(nbits);
    nbits = SHFL(nbits, 0);
    tpos = (threadIdx.x & warp_bits_mask_int);
    while (tpos < 8) {
      const int ck = tpos / 4;
      const int cj = (tpos - (4 * ck)) / 2;
      const int ci = (tpos & 1);
      int im = trg_i - ci;
      int jm = trg_j - cj;
      int km = trg_k - ck;
      switch (contribution.dims.bounds) {
      case BoundaryCondition::ISOLATED:
        break;
      case BoundaryCondition::PERIODIC:
        im += (im < 0) * target.dims.na;
        jm += (jm < 0) * target.dims.nb;
        km += (km < 0) * target.dims.nc;
        break;
      }

      // Proceed only if the element is valid
      if (im >= 0 && jm >= 0 && km >= 0) {
        const size_t target_idx = (((km * target.dims.nb) + jm) * target.dims.na) + im;
        const size_t target_subidx = 8 * ((((ck * 2) + cj) * 2) + ci);
        target.coeffs[(64LLU * target_idx) + target_subidx] = nbits;
      }
      tpos += warp_size_int;
    }

    // The warp proceeds to the next mesh element
    wpos += (blockDim.x >> warp_bits) * gridDim.x;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchAccOcclusionMesh(BackgroundMeshWriter<llint> *target,
                                   const BackgroundMeshReader<llint> &contribution,
                                   const GpuDetails &gpu) {
  checkMeshCompatibility<llint, llint>(*target, contribution);
  kAccumulateOcclusionMesh<<<gpu.getSMPCount(), large_block_size>>>(*target, contribution);
}

//-------------------------------------------------------------------------------------------------
extern void launchColorNonbondedFieldMesh(BackgroundMeshWriter<void> *target, const size_t ct_targ,
                                          const MeshFFKit<double> &mnbk,
                                          const PrecisionModel prec, const double* stn_xfrm,
                                          const AtomGraph *ag, const CoordinateFrameReader &cfr,
                                          const MeshKlManager &launcher,
                                          const HybridTargetLevel availability) {

  // Make a switch over the allowed types of the mesh data, then execute the appropriate kernel to
  // compute the non-bonded field.
  const int2 lp = launcher.getMeshKernelDims(prec, GridDetail::NONBONDED_FIELD,
                                             target->dims.bounds, target->dims.unit_cell,
                                             target->dims.stencil_kind,
                                             MappingActivity::PARTICLE_TO_MESH);
  if (ct_targ == double_type_index) {
    BackgroundMeshWriter<double> bgmw = restoreType<double>(target);
    unrollColorNBFMesh<double>(&bgmw, mnbk, prec, stn_xfrm, ag, cfr, lp, availability);
  }
  else if (ct_targ == float_type_index) {
    BackgroundMeshWriter<float> bgmw = restoreType<float>(target);
    unrollColorNBFMesh<float>(&bgmw, mnbk, prec, stn_xfrm, ag, cfr, lp, availability);
  }
  else if (ct_targ == llint_type_index) {
    BackgroundMeshWriter<llint> bgmw = restoreType<llint>(target);
    unrollColorNBFMesh<llint>(&bgmw, mnbk, prec, stn_xfrm, ag, cfr, lp, availability);
  }
  else {
    rtErr("The allowed data types for non-bonded fields are float, double, and long long integer.",
          "launchColorNonbondedFieldMesh");
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchAccNonbondedFieldMesh(BackgroundMeshWriter<llint> *target,
                                        const BackgroundMeshReader<void> &contribution,
                                        const size_t ct_contrib, const GpuDetails &gpu) {
  checkMeshCompatibility<llint, void>(*target, contribution);  
  const int nblocks = gpu.getSMPCount();
  if (ct_contrib == double_type_index) {
    BackgroundMeshReader<double> contrib_rt = restoreType<double>(contribution);
    kAccumulateNonbondedFieldMesh<<<nblocks, large_block_size>>>(*target, contrib_rt);
  }
  else if (ct_contrib == float_type_index) {
    BackgroundMeshReader<float> contrib_rt = restoreType<float>(contribution);
    kAccumulateNonbondedFieldMesh<<<nblocks, large_block_size>>>(*target, contrib_rt);
  }
  else if (ct_contrib == llint_type_index) {
    BackgroundMeshReader<llint> contrib_rt = restoreType<llint>(contribution);
    kAccumulateNonbondedFieldMesh<<<nblocks, large_block_size>>>(*target, contrib_rt);
  }
}

//-------------------------------------------------------------------------------------------------
void launchOccFieldDerivativeCalc(BackgroundMeshWriter<void> *bgmw, const double max_occlusion,
                                  size_t ct_mesh, const GpuDetails &gpu) {
  if (ct_mesh == double_type_index) {
    BackgroundMeshWriter<double> occfield = restoreType<double>(bgmw);
    kOccFieldDerivativeCalc<<<2 * gpu.getSMPCount(), medium_block_size>>>(occfield, max_occlusion);
  }
  else if (ct_mesh == float_type_index) {
    BackgroundMeshWriter<float> occfield = restoreType<float>(bgmw);
    kOccFieldDerivativeCalc<float><<<2 * gpu.getSMPCount(), medium_block_size>>>(occfield,
                                                                                 max_occlusion);
  }
  else {
    rtErr("The only available types for this BackgroundMesh operation are float and double.",
          "launchOccFieldDerivativeCalc");
  }
}
  
} // namespace structure
} // namespace stormm
