// -*-c++-*-
#ifndef STORMM_HPC_BACKGROUND_MESH_CUH
#define STORMM_HPC_BACKGROUND_MESH_CUH

#include "Accelerator/ptx_macros.h"
#include "DataTypes/common_types.h"
#include "Reporting/error_format.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "background_mesh.h"

namespace stormm {
namespace structure {

#include "Accelerator/syncwarp.cui"
#include "Numerics/accumulation.cui"

using topology::AtomGraph;
using topology::getEnumerationName;
using topology::NonbondedKit;
using trajectory::CoordinateFrameReader;

#include "Math/radial_derivatives.cui"
#include "local_arrangement.cui"

/// \brief Accumulate a non-bonded field based on the contribution of one system and the average of
///        many systems.
///
/// \param target        The accumulating average of the non-bonded potential due to many snapshots
///                      of the rigid molecule of interest
/// \param contribution  The non-bonded potential due to one snapshot of the molecule
template <typename T> __global__ __launch_bounds__(large_block_size, 1)
void kAccumulateNonbondedFieldMesh(BackgroundMeshWriter<llint> target,
                                   const BackgroundMeshReader<T> contribution) {
  const size_t ncoeff = (size_t)(target.dims.na * target.dims.nb * target.dims.nc) * (size_t)(64);
  const size_t init_pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t total_threads = gridDim.x * blockDim.x;
  const bool data_is_real = (contribution.coeff_scale < 1.001);
  if (data_is_real) {
    for (size_t pos = init_pos; pos < ncoeff; pos += total_threads) {
      target.coeffs[pos] += __double2ll_rn(contribution.coeffs[pos] * target.coeff_scale);
    }
  }
  else {
    const int target_bits = roundf(log2f(target.coeff_scale));
    const int contrib_bits = roundf(log2f(contribution.coeff_scale));
    if (contrib_bits > target_bits) {
      const int bit_diff = contrib_bits - target_bits;
      for (size_t pos = init_pos; pos < ncoeff; pos += total_threads) {
        const llint tc = contribution.coeffs[pos];
        target.coeffs[pos] += (tc >> bit_diff);
      }
    }
    else if (contrib_bits < target_bits) {
      const int bit_diff = target_bits - contrib_bits;
      for (size_t pos = init_pos; pos < ncoeff; pos += total_threads) {
        const llint tc = contribution.coeffs[pos];
        target.coeffs[pos] += (tc << bit_diff);
      }      
    }
    else {
      for (size_t pos = init_pos; pos < ncoeff; pos += total_threads) {
        target.coeffs[pos] += contribution.coeffs[pos];
      }
    }
  }
}

/// \brief Compute the derivatives due to an occlusion field.  This applies a finite-difference
///        approximation to compute the necessary derivatives that set up a tricubic interpolant.
///
/// \param occfield       The occlusion field mesh to compute
/// \param max_occlusion  The maximum value of occlusion
template <typename T> __global__ __launch_bounds__(medium_block_size, 2)
void kOccFieldDerivativeCalc(BackgroundMeshWriter<T> occfield, const T max_occlusion) {

  // Allocate multiple arrays to store the potential and its derivatives, to get the most out of
  // L1 caching, particularly in terms of coalesced writing of the derivatives.
  __shared__ volatile T u_patch[1000], dxu_patch[1000], dyu_patch[1000], dzu_patch[1000];
  
  // Pre-computations to determine the numbers of patches
  int na_stride, nb_stride, nc_stride;
  switch (occfield.dims.bounds) {
  case BoundaryCondition::ISOLATED:
    na_stride = (occfield.dims.na + 8) / 8 ;
    nb_stride = (occfield.dims.nb + 8) / 8 ;
    nc_stride = (occfield.dims.nc + 8) / 8 ;
    break;
  case BoundaryCondition::PERIODIC:
    na_stride = (occfield.dims.na + 7) / 8 ;
    nb_stride = (occfield.dims.nb + 7) / 8 ;
    nc_stride = (occfield.dims.nc + 7) / 8 ;
    break;
  }

  // Loop over all patches
  int block_pos = blockIdx.x;
  while (block_pos < na_stride * nb_stride * nc_stride) {

    // Fill up the potential mesh.  If necessary, make inferences about values on the periphery of
    // the patch.
    const int bpos_c = block_pos / (na_stride * nb_stride);
    const int bpos_b = (block_pos - (bpos_c * na_stride * nb_stride)) / na_stride;
    const int bpos_a = block_pos - (((bpos_c * nb_stride) + bpos_b) * na_stride);
    int tmpr_a_start, tmpr_b_start, tmpr_c_start;
    int patch_dim_a = occfield.dims.na - (bpos_a * 8);
    int patch_dim_b = occfield.dims.nb - (bpos_b * 8);
    int patch_dim_c = occfield.dims.nc - (bpos_c * 8);
    switch (occfield.dims.bounds) {
    case BoundaryCondition::ISOLATED:
      tmpr_a_start = (patch_dim_a == 1);
      tmpr_b_start = (patch_dim_b == 1);
      tmpr_c_start = (patch_dim_c == 1);
      break;
    case BoundaryCondition::PERIODIC:
      tmpr_a_start = 0;
      tmpr_b_start = 0;
      tmpr_c_start = 0;
      break;
    }
    patch_dim_a += tmpr_a_start + 2;
    patch_dim_b += tmpr_b_start + 2;
    patch_dim_c += tmpr_c_start + 2;
    const int mesh_a_start = (8 * bpos_a) - tmpr_a_start - 1;
    const int mesh_b_start = (8 * bpos_b) - tmpr_b_start - 1;
    const int mesh_c_start = (8 * bpos_c) - tmpr_c_start - 1;

    // Each warp will loop over a single mesh element and read its eight energies, then decide
    // where they can be assigned.  Skip odd-numbered relative reads--exploit the fact that
    // information is present in octuplicate.
    int warp_pos = (threadIdx.x >> warp_bits);
    const int n_warps = (blockDim.x >> warp_bits);
    const int lane_idx = (threadIdx.x & warp_bits_mask_int);
    while (warp_pos < patch_dim_a * patch_dim_b * patch_dim_c) {
      const int rel_read_c = warp_pos / (patch_dim_a * patch_dim_b);
      if (rel_read_c & 1) {
        warp_pos += n_warps;
        continue;
      }
      const int rel_read_b = (warp_pos - (rel_read_c * patch_dim_a * patch_dim_b)) / patch_dim_a;
      if (rel_read_b & 1) {
        warp_pos += n_warps;
        continue;
      }
      const int rel_read_a = warp_pos - (((rel_read_c * patch_dim_b) + rel_read_b) * patch_dim_a);
      if (rel_read_a & 1) {
        warp_pos += n_warps;
        continue;
      }

      // Verify that the mesh element is valid
      int elem_read_a = mesh_a_start + rel_read_a;
      int elem_read_b = mesh_b_start + rel_read_b;
      int elem_read_c = mesh_c_start + rel_read_c;
      switch (occfield.dims.bounds) {
      case BoundaryCondition::ISOLATED:
        if (elem_read_a < 0 || elem_read_a >= occfield.dims.na ||
            elem_read_b < 0 || elem_read_b >= occfield.dims.nb ||
            elem_read_c < 0 || elem_read_c >= occfield.dims.nc) {
          warp_pos += n_warps;
          continue;
        }
        break;
      case BoundaryCondition::PERIODIC:
        elem_read_a += occfield.dims.na * ((elem_read_a < 0) - (elem_read_a >= occfield.dims.na));
        elem_read_b += occfield.dims.nb * ((elem_read_b < 0) - (elem_read_b >= occfield.dims.nb));
        elem_read_c += occfield.dims.nc * ((elem_read_c < 0) - (elem_read_c >= occfield.dims.nc));
        break;
      }

      // Read the function values from the mesh element.
      int pos = lane_idx;
      while (pos < 8) {
        const size_t elem_idx = (size_t)((((elem_read_c * occfield.dims.nb) + elem_read_b) *
                                          occfield.dims.na) + elem_read_a) * (size_t)(64);
        const T f_value = occfield.coeffs[elem_idx + (size_t)(pos)];

        // Take into account the origin of the patch relative to the origin of the mesh
        const int pos_c_bump = pos / 4;
        const int pos_b_bump = (pos - (pos_c_bump * 4)) / 2;
        const int pos_a_bump = (pos & 1);
        const int patch_a_idx = rel_read_a + pos_a_bump;
        const int patch_b_idx = rel_read_b + pos_b_bump;
        const int patch_c_idx = rel_read_c + pos_c_bump;
        if (patch_a_idx < patch_dim_a && patch_b_idx < patch_dim_b && patch_c_idx < patch_dim_c) {
          u_patch[(((patch_c_idx * patch_dim_b) +
                    patch_b_idx) * patch_dim_a) + patch_a_idx] = f_value;
        }
        pos += warp_size_int;
      }
      
      // Increment the warp counter
      warp_pos += n_warps;
    }
    __syncthreads();

    // Infer additional potentials if required
    int thrd_pos = threadIdx.x;
    while (thrd_pos < patch_dim_a * patch_dim_b * patch_dim_c) {
      const int patch_c_idx = thrd_pos / (patch_dim_a * patch_dim_b);
      const int patch_b_idx = (thrd_pos - (patch_c_idx * patch_dim_a * patch_dim_b)) / patch_dim_a;
      const int patch_a_idx = thrd_pos -
                              (((patch_c_idx * patch_dim_b) + patch_b_idx) * patch_dim_a);

      // Detect whether the point corresponds to an invalid mesh point
      const int elem_read_a = mesh_a_start + patch_a_idx;
      const int elem_read_b = mesh_b_start + patch_b_idx;
      const int elem_read_c = mesh_c_start + patch_c_idx;
      switch (occfield.dims.bounds) {
      case BoundaryCondition::ISOLATED:
        if (elem_read_a < 0 || elem_read_a >= occfield.dims.na ||
            elem_read_b < 0 || elem_read_b >= occfield.dims.nb ||
            elem_read_c < 0 || elem_read_c >= occfield.dims.nc) {
          int fulcrum_a, fulcrum_b, fulcrum_c, lever_a, lever_b, lever_c;
          if (patch_a_idx == 0) {
            fulcrum_a = 1;
            lever_a = 2;
          }
          else if (patch_a_idx == patch_dim_a - 1) {
            fulcrum_a = patch_dim_a - 2;
            lever_a = patch_dim_a - 3;
          }
          else {
            fulcrum_a = patch_a_idx;
            lever_a = patch_a_idx;
          }
          if (patch_b_idx == 0) {
            fulcrum_b = 1;
            lever_b = 2;
          }
          else if (patch_b_idx == patch_dim_b - 1) {
            fulcrum_b = patch_dim_b - 2;
            lever_b = patch_dim_b - 3;
          }
          else {
            fulcrum_b = patch_b_idx;
            lever_b = patch_b_idx;
          }
          if (patch_c_idx == 0) {
            fulcrum_c = 1;
            lever_c = 2;
          }
          else if (patch_c_idx == patch_dim_c - 1) {
            fulcrum_c = patch_dim_c - 2;
            lever_c = patch_dim_c - 3;
          }
          else {
            fulcrum_c = patch_c_idx;
            lever_c = patch_c_idx;
          }
          const T u_fulcrum = u_patch[(((fulcrum_c * patch_dim_b) + fulcrum_b) * patch_dim_a) +
                                      fulcrum_a];
          const T u_lever = u_patch[(((lever_c * patch_dim_b) + lever_b) * patch_dim_a) +
                                      lever_a];
          u_patch[(((patch_c_idx * patch_dim_b) + patch_b_idx) * patch_dim_a) + patch_a_idx] =
            ((T)(2.0) * u_fulcrum) - u_lever;
        }
        break;
      case BoundaryCondition::PERIODIC:
        break;
      }

      // Increment to the next element in the patch
      thrd_pos += blockDim.x;
    }
    
    // Evaluate df/dx, df/dy, and df/dz, storing the results in works arrays dxu, dyu, and dzu.
    // This phase of the calculation gives those __shared__ arrays their nomenclature.
    for (int abc_idx = threadIdx.x; abc_idx < patch_dim_a * patch_dim_b * patch_dim_c;
         abc_idx += blockDim.x) {
      const int patch_c_idx = abc_idx / (patch_dim_a * patch_dim_b);
      const int patch_b_idx = (abc_idx - (patch_dim_a * patch_dim_b * patch_c_idx)) / patch_dim_a;
      const int patch_a_idx = abc_idx - (((patch_c_idx * patch_dim_b) +
                                          patch_b_idx) * patch_dim_a);
      if (patch_a_idx > 0 && patch_a_idx < patch_dim_a - 1) {
        const size_t abc_mi = abc_idx - 1;
        const size_t abc_pi = abc_idx + 1;
        dxu_patch[abc_idx] = (T)(0.5) * (u_patch[abc_pi] - u_patch[abc_mi]);
      }
      if (patch_b_idx > 0 && patch_b_idx < patch_dim_b - 1) {
        const size_t abc_mj = (((patch_c_idx * patch_dim_b) + patch_b_idx - 1) *
                               patch_dim_a) + patch_a_idx;
        const size_t abc_pj = (((patch_c_idx * patch_dim_b) + patch_b_idx + 1) *
                               patch_dim_a) + patch_a_idx;
        dyu_patch[abc_idx] = (T)(0.5) * (u_patch[abc_pj] - u_patch[abc_mj]);
      }
      if (patch_c_idx > 0 && patch_c_idx < patch_dim_c - 1) {
        const size_t abc_mk = ((((patch_c_idx - 1) * patch_dim_b) + patch_b_idx) *
                               patch_dim_a) + patch_a_idx;
        const size_t abc_pk = ((((patch_c_idx + 1) * patch_dim_b) + patch_b_idx) *
                               patch_dim_a) + patch_a_idx;
        dzu_patch[abc_idx] = (T)(0.5) * (u_patch[abc_pk] - u_patch[abc_mk]);
      }
    }
    __syncthreads();

    // Store the first derivatives

    
    // Evaluate d2f/dxdy, storing the result in the second derivative work array.
    
    
    // Move on to the next patch.
    block_pos += gridDim.x;    
  }
}

#define TCALC float
#  define TCALC2 float2
#  define TCALC4 float4
#  define TCALC_IS_SINGLE
#  define LLCONV_FUNC __float2ll_rn
#  define SQRT_FUNC sqrtf
#  define FLOOR_FUNC floorf
#  define CEIL_FUNC ceilf
#  define SMOOTHNESS_INTP
#    define TRICLINIC_ELEMENT
#      define COLOR_NBF_THREAD_COUNT 896
#       define KERNEL_NAME kfColorNBFieldMeshTric
#        include "color_nbfield_mesh.cui"
#      undef KERNEL_NAME
#      undef COLOR_NBF_THREAD_COUNT
#      define PERIODIC_MESH
#        define COLOR_NBF_THREAD_COUNT 768
#        define KERNEL_NAME kfColorNBFieldMeshPbcTric
#          include "color_nbfield_mesh.cui"
#        undef KERNEL_NAME
#        undef COLOR_NBF_THREAD_COUNT
#      undef PERIODIC_MESH
#    undef TRICLINIC_ELEMENT
#    define COLOR_NBF_THREAD_COUNT 1024
#    define KERNEL_NAME kfColorNBFieldMesh
#      include "color_nbfield_mesh.cui"
#    undef KERNEL_NAME
#    undef COLOR_NBF_THREAD_COUNT
#    define PERIODIC_MESH
#      define COLOR_NBF_THREAD_COUNT 896
#      define KERNEL_NAME kfColorNBFieldMeshPbc
#        include "color_nbfield_mesh.cui"
#      undef KERNEL_NAME
#      undef COLOR_NBF_THREAD_COUNT
#    undef PERIODIC_MESH
#  undef SMOOTHNESS_INTP
#  define VALUE_INTP
#    define COLOR_NBF_THREAD_COUNT 1024
#    define TRICLINIC_ELEMENT
#       define KERNEL_NAME kfColorNBFieldMeshTricValue
#        include "color_nbfield_mesh.cui"
#      undef KERNEL_NAME
#      define PERIODIC_MESH
#        define KERNEL_NAME kfColorNBFieldMeshPbcTricValue
#          include "color_nbfield_mesh.cui"
#        undef KERNEL_NAME
#      undef PERIODIC_MESH
#    undef TRICLINIC_ELEMENT
#    define KERNEL_NAME kfColorNBFieldMeshValue
#      include "color_nbfield_mesh.cui"
#    undef KERNEL_NAME
#    define PERIODIC_MESH
#      define KERNEL_NAME kfColorNBFieldMeshPbcValue
#        include "color_nbfield_mesh.cui"
#      undef KERNEL_NAME
#    undef PERIODIC_MESH
#    undef COLOR_NBF_THREAD_COUNT
#  undef VALUE_INTP
#  undef CEIL_FUNC
#  undef FLOOR_FUNC
#  undef SQRT_FUNC
#  undef LLCONV_FUNC
#  undef TCALC_IS_SINGLE
#  undef TCALC4
#  undef TCALC2
#undef TCALC

#define TCALC double
#  define TCALC2 double2
#  define TCALC4 double4
#  define LLCONV_FUNC __double2ll_rn
#  define SQRT_FUNC sqrt
#  define FLOOR_FUNC floor
#  define CEIL_FUNC ceil
#  define COLOR_NBF_THREAD_COUNT 512
#  define SMOOTHNESS_INTP
#    define TRICLINIC_ELEMENT
#      define KERNEL_NAME kdColorNBFieldMeshTric
#        include "color_nbfield_mesh.cui"
#      undef KERNEL_NAME
#      define PERIODIC_MESH
#        define KERNEL_NAME kdColorNBFieldMeshPbcTric
#          include "color_nbfield_mesh.cui"
#        undef KERNEL_NAME
#      undef PERIODIC_MESH
#    undef TRICLINIC_ELEMENT
#    define KERNEL_NAME kdColorNBFieldMesh
#      include "color_nbfield_mesh.cui"
#    undef KERNEL_NAME
#    define PERIODIC_MESH
#      define KERNEL_NAME kdColorNBFieldMeshPbc
#        include "color_nbfield_mesh.cui"
#      undef KERNEL_NAME
#    undef PERIODIC_MESH
#  undef SMOOTHNESS_INTP
#  undef COLOR_NBF_THREAD_COUNT
#  define COLOR_NBF_THREAD_COUNT 768
#  define VALUE_INTP
#    define TRICLINIC_ELEMENT
#      define KERNEL_NAME kdColorNBFieldMeshTricValue
#        include "color_nbfield_mesh.cui"
#      undef KERNEL_NAME
#      define PERIODIC_MESH
#        define KERNEL_NAME kdColorNBFieldMeshPbcTricValue
#          include "color_nbfield_mesh.cui"
#        undef KERNEL_NAME
#      undef PERIODIC_MESH
#    undef TRICLINIC_ELEMENT
#    define KERNEL_NAME kdColorNBFieldMeshValue
#      include "color_nbfield_mesh.cui"
#    undef KERNEL_NAME
#    define PERIODIC_MESH
#      define KERNEL_NAME kdColorNBFieldMeshPbcValue
#        include "color_nbfield_mesh.cui"
#      undef KERNEL_NAME
#    undef PERIODIC_MESH
#  undef VALUE_INTP
#  undef COLOR_NBF_THREAD_COUNT
#  undef CEIL_FUNC
#  undef FLOOR_FUNC
#  undef SQRT_FUNC
#  undef LLCONV_FUNC
#  undef TCALC4
#  undef TCALC2
#undef TCALC
  
/// \brief Unroll the innermost switch for launching non-bonded field mesh plotting kernels of
///        different precision models.  The data type of the mesh non-bonded parameter kit carries
///        with it the intended calculation mode.  
///
/// \param bgmw          The background mesh to color with the non-bonded field of interest
/// \param mnbk          Mesh non-bonded softcore potential parameters, other essential constants
/// \param prec          Precision model in which to compute the mesh-based potentials
/// \param stn_xfrm      Stencil transformation matrix for converting observations about function
///                      values and various derivatives at mesh points into set of coefficients for
///                      polynomial approximations within each mesh element
/// \param ag            Topology referenced by the mesh (this pointer is not included within the
///                      mesh abstract, but a copy of the frozen atom content of the topology is)
/// \param cfr           Coordinates of the molecule to be used in coloring
/// \param lp            Launch parameters 
/// \param availability  Indicate whether the topology's information is available on the CPU host
///                      or GPU device
template <typename T>
void unrollColorNBFMesh(BackgroundMeshWriter<T> *bgmw, const MeshFFKit<double> &mnbk,
                        const PrecisionModel prec, const double* stn_xfrm, const AtomGraph *ag,
                        const CoordinateFrameReader &cfr, const int2 lp,
                        const HybridTargetLevel availability) {
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  bool problem = false;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const NonbondedKit<double> nbk = (availability == HybridTargetLevel::DEVICE) ?
                                       ag->getDoublePrecisionNonbondedKit(devc_tier) :
                                       ag->getDeviceViewToHostDPNonbondedKit();

      // The switch over boundary conditions is unrolled for different kernels, as the periodic
      // mesh case may require more registers.
      switch (bgmw->dims.bounds) {
      case BoundaryCondition::ISOLATED:
        switch (bgmw->dims.unit_cell) {
        case UnitCellType::ORTHORHOMBIC:
          switch (bgmw->dims.stencil_kind) {
          case Interpolant::SMOOTHNESS:
            kdColorNBFieldMesh<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          case Interpolant::FUNCTION_VALUE:
            kdColorNBFieldMeshValue<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          }
          break;
        case UnitCellType::TRICLINIC:
          switch (bgmw->dims.stencil_kind) {
          case Interpolant::SMOOTHNESS:
            kdColorNBFieldMeshTric<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          case Interpolant::FUNCTION_VALUE:
            kdColorNBFieldMeshTricValue<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          }
          break;
        case UnitCellType::NONE:
          problem = true;
          break;
        }
        break;
      case BoundaryCondition::PERIODIC:
        switch (bgmw->dims.unit_cell) {
        case UnitCellType::ORTHORHOMBIC:
          switch (bgmw->dims.stencil_kind) {
          case Interpolant::SMOOTHNESS:
            kdColorNBFieldMeshPbc<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          case Interpolant::FUNCTION_VALUE:
            kdColorNBFieldMeshPbcValue<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          }
          break;
        case UnitCellType::TRICLINIC:
          switch (bgmw->dims.stencil_kind) {
          case Interpolant::SMOOTHNESS:
            kdColorNBFieldMeshPbcTric<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          case Interpolant::FUNCTION_VALUE:
            kdColorNBFieldMeshPbcTricValue<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          }
          break;
        case UnitCellType::NONE:
          problem = true;
          break;
        }
        break;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const NonbondedKit<float> nbk = (availability == HybridTargetLevel::DEVICE) ?
                                      ag->getSinglePrecisionNonbondedKit(devc_tier) :
                                      ag->getDeviceViewToHostSPNonbondedKit();
      switch (bgmw->dims.bounds) {
      case BoundaryCondition::ISOLATED:
        switch (bgmw->dims.unit_cell) {
        case UnitCellType::ORTHORHOMBIC:
          switch (bgmw->dims.stencil_kind) {
          case Interpolant::SMOOTHNESS:
            kfColorNBFieldMesh<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          case Interpolant::FUNCTION_VALUE:
            kfColorNBFieldMeshValue<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          }
          break;
        case UnitCellType::TRICLINIC:
          switch (bgmw->dims.stencil_kind) {
          case Interpolant::SMOOTHNESS:
            kfColorNBFieldMeshTric<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          case Interpolant::FUNCTION_VALUE:
            kfColorNBFieldMeshTricValue<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          }
          break;
        case UnitCellType::NONE:
          problem = true;
          break;
        }
        break;
      case BoundaryCondition::PERIODIC:
        switch (bgmw->dims.unit_cell) {
        case UnitCellType::ORTHORHOMBIC:
          switch (bgmw->dims.stencil_kind) {
          case Interpolant::SMOOTHNESS:
            kfColorNBFieldMeshPbc<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          case Interpolant::FUNCTION_VALUE:
            kfColorNBFieldMeshPbcValue<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          }
          break;
        case UnitCellType::TRICLINIC:
          switch (bgmw->dims.stencil_kind) {
          case Interpolant::SMOOTHNESS:
            kfColorNBFieldMeshPbcTric<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          case Interpolant::FUNCTION_VALUE:
            kfColorNBFieldMeshPbcTricValue<T><<<lp.x, lp.y>>>(*bgmw, mnbk, stn_xfrm, nbk, cfr);
            break;
          }
          break;
        case UnitCellType::NONE:
          problem = true;
          break;
        }
        break;
      }
    }
    break;
  }
  if (problem) {
    rtErr("A mesh cannot have UnitCellType " + getEnumerationName(bgmw->dims.unit_cell) + ".",
          "unrollColorNBFMesh");
  }
}

} // namespace structure
} // namespace stormm

#endif
