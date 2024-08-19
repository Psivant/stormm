// -*-c++-*-
#ifndef STORMM_HPC_COORDINATE_COPY_CUH
#define STORMM_HPC_COORDINATE_COPY_CUH

#include "Accelerator/gpu_details.h"
#include "Constants/behavior.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Numerics/split_fixed_precision.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "coordinate_series.h"

namespace stormm {
namespace trajectory {

using card::GpuDetails;
using constants::PrecisionModel;
using synthesis::CondensateWriter;
using synthesis::CondensateReader;
using synthesis::PsSynthesisWriter;
using synthesis::PsSynthesisReader;

#include "Math/rounding.cui"
#include "Numerics/accumulation.cui"

/// \brief Copy the box information between coordinate types that contain it.  This __device__
///        function consolidates a great deal of repetitive code.
///
/// \param umat_dest         Box space transformation(s) for the destinations system(s)
/// \param invu_dest         Transformations back to real space for the destination system(s)
//   bdim_dest:         Box dimensions for the destination system(s)
//   umat_orig:         Box space transformation(s) for the original system(s)
//   invu_orig:         Transformations back to real space for the original system(s)
//   bdim_orig:         Box dimensions for the original system(s)
//   dest_xfrm_offset:  Index at which to start writing matrix values to umat_dest, invu_dest,
//                      and (if applicable) dest_boxvecs and dest_boxvec_ovrf
//   dest_bdim_offset:  Index at which to start writing box dimensions to bdim_dest
//   orig_xfrm_offset:  Index at which to start reading matrix values from umat_orig, invu_orig,
//                      and (if applicable) orig_boxvecs and orig_boxvec_ovrf
//   orig_bdim_offset:  Index at which to start reading box dimensions from bdim_orig
//   warp_offset:       This __device__ function assigns each task (copy box-space transform,
//                      copy back-to-real-space inverse transform, copy box dimensions) to a
//                      separate warp.  This function will engage warps warp_offset,
//                      warp_offset + 1, and warp_offset + 2, which must be valid within the
//                      launch parameters.
//   dest_boxvecs:      Integral representation of box vectors.  This will be set according to the
//                      double-precision inverse transformation matrix if that is all that is
//                      available, or based on other integral representations of the box vectors if
//                      those are provided.  This representation matches the way that the positions
//                      of particles are actually represented in the PhaseSpaceSynthesis, and is
//                      authoritative over the double-precision inverse transformation matrices.
//   dest_boxvec_ovrf:  Overflow bits for the destination's integral box vector representation
//   dest_scale:        Scaling factor for the destination coordinates, used only if the integral
//                      box vector representation is provided
//   orig_boxvecs:      Integral representations of the box vectors for the system(s) in the origin
//   orig_boxvec_ovrf:  Overflow bits for the origin's integral box vector representation
//   dest_bits:         Bit count for scaling the destination positions (and box vectors)
//   orig_bits:         Bit count for scaling the origin positions (and box vectors)
__device__ __forceinline__
void kCopyBoxInformation(double* umat_dest, double* invu_dest, double* bdim_dest,
                         const double* umat_orig, const double* invu_orig, const double* bdim_orig,
                         const int dest_xfrm_offset = 0, const int orig_xfrm_offset = 0,
                         const int dest_bdim_offset = 0, const int orig_bdim_offset = 0,
                         const int warp_offset = 0, llint* dest_boxvecs = nullptr,
                         int* dest_boxvec_ovrf = nullptr, const double dest_scale = 1.0,
                         const llint* orig_boxvecs = nullptr,
                         const int* orig_boxvec_ovrf = nullptr, const int dest_bits = 0,
                         const int orig_bits = 0);
  
/// \brief Copy a set of coordinates from one array into another, inferring the conversion
///        operation based on the origin and destination data types.
///
/// Overloaded:
///   - Convert one general type into another (here, a "general type" is a coordinate format
///     which expresses the X, Y, or Z component of a position, velocity, or force as a scalar
///     type, e.g. float, llint, double)
///   - Convert a general type into 95-bit fixed-precision (int95_t)
///   - Convert 95-bit fixed-precision into a general type
///
/// \param dest_crd         Destination array to hold coordinates
/// \param dest_crd_ovrf    Overflow bits for dest_crd (its presence implies dest_crd is llint)
/// \param dest_start_idx   Starting index of the destination array at which to begin writing
///                         (likely to be zero for single-system coordinate objects)
/// \param dest_scale       Scaling factor for coordinates in the destination array
/// \param dest_scale_bits  Number of bits after the decimal for fixed-precision coordinates in
///                         the destination array (the base-2 logarithm of dest_scale)
/// \param orig_crd         Array holding the original coordinates
/// \param orig_crd_ovrf    Overflow bits for orig_crd (its presence implies dest_crd is llint)
/// \param orig_start_idx   Starting index of the origin array at which to begin reading (likely
///                         to be zero for single-system coordinate objects)
/// \param orig_scale       Scaling factor for coordinates in the origin array
/// \param orig_scale_bits  Number of bits after the decimal for fixed-precision coordinates in
///                         the origin array (the base-2 logarithm of dest_scale)
/// \param count            The number of particle coordinates to copy
/// \param pos              Internal counter passed in from the calling function
/// \param iter             The number of passes through this array, each incrementing the counter
///                         behind pos
/// \param stride           The padded form of count (pre-computed and passed in for convenience
///                         and probably reduced register pressure)
/// \param advance          The amount to increment an internal variant of the pos counter (this
///                         can be the width of a warp, the width of a block, or the width of the
///                         entire kernel launch grid depending on the context in which this copy
///                         routine is called)
/// \{
template <typename Tdest, typename Torig> __device__ __forceinline__
size_t copyCoordinateSet(Tdest* dest_crd, const int dest_start_idx, const double dest_scale,
                         const int dest_scale_bits, const Torig* orig_crd,
                         const int orig_start_idx, const double orig_scale,
                         const int orig_scale_bits, const int count, const size_t pos,
                         const size_t iter, const size_t stride, const size_t advance) {
  const size_t pos_limit = (iter + 1) * stride;
  size_t ipos = pos;
  while (ipos < pos_limit) {
    const int rel_pos = ipos - (iter * stride);
    if (rel_pos < count) {
      const size_t drp_idx = dest_start_idx + rel_pos;
      const size_t orp_idx = orig_start_idx + rel_pos;
      if (orig_scale_bits == 0) {
        if (dest_scale_bits == 0) {
          __stwt(&dest_crd[drp_idx], __ldlu(&orig_crd[orp_idx]));
        }
        else {
          const llint ival = __double2ll_rn(__ldlu(&orig_crd[orp_idx]) * dest_scale);
          __stwt(&dest_crd[drp_idx], ival);
        }
      }
      else {
        if (dest_scale_bits == 0) {
          __stwt(&dest_crd[drp_idx], (double)(__ldlu(&orig_crd[orp_idx])) / orig_scale);
        }
        else {
          if (dest_scale_bits > orig_scale_bits) {
            const llint dconv = __ldlu(&orig_crd[orp_idx]);
            __stwt(&dest_crd[drp_idx], (dconv << (dest_scale_bits - orig_scale_bits)));
          }
          else if (dest_scale_bits == orig_scale_bits) {
            __stwt(&dest_crd[drp_idx], __ldlu(&orig_crd[orp_idx]));
          }
          else {
            const llint dconv = __ldlu(&orig_crd[orp_idx]);
            __stwt(&dest_crd[drp_idx], (dconv >> (orig_scale_bits - dest_scale_bits)));
          }
        }
      }
    }
    ipos += advance;
  }
  return ipos;
}

template <typename Torig> __device__ __forceinline__
size_t copyCoordinateSet(llint* dest_crd, int* dest_crd_ovrf, const int dest_start_idx,
                         const double dest_scale, const int dest_scale_bits, const Torig* orig_crd,
                         const int orig_start_idx, const int orig_scale_bits, const int count,
                         const size_t pos, const size_t iter, const size_t stride,
                         const size_t advance) {
  const size_t pos_limit = (iter + 1) * stride;
  size_t ipos = pos;
  while (ipos < pos_limit) {
    const int rel_pos = ipos - (iter * stride);
    if (rel_pos < count) {
      const size_t drp_idx = dest_start_idx + rel_pos;

      // Detect the presence of a signed integer, fixed-precision representation in the origin
      // by a nonzero number of scaling bits.  Copy that directly into the long long int portion
      // of a 95-bit representation and upgrade (or downgrade) the bit count as appropriate.
      if (orig_scale_bits > 0) {
        const int95_t oval = { (long long int)(__ldlu(&orig_crd[orig_start_idx + rel_pos])), 0 };
        const int95_t dval = changeFPBits(oval, orig_scale_bits, dest_scale_bits);
        __stwt(&dest_crd[drp_idx], dval.x);
        __stwt(&dest_crd_ovrf[drp_idx], dval.y);
      }
      else {

        // Ignore limits on various types of coordinates which could support an assumption that
        // the high 32 bits are unnecessary.  The memory bandwidth will be the limiting factor.
        // Treat every conversion as a double-to-int95_t operation.
        const int95_t dval = doubleToInt95(__ldlu(&orig_crd[orig_start_idx + rel_pos]) *
                                           dest_scale);
        __stwt(&dest_crd[drp_idx], dval.x);
        __stwt(&dest_crd_ovrf[drp_idx], dval.y);
      }
    }
    ipos += advance;
  }
  return ipos;
}

template <typename Tdest> __device__ __forceinline__
size_t copyCoordinateSet(Tdest* dest_crd, const int dest_start_idx, const int dest_scale_bits,
                         const llint* orig_crd, const int* orig_crd_ovrf, const int orig_start_idx,
                         const double orig_scale, const int orig_scale_bits, const int count,
                         const size_t pos, const size_t iter, const size_t stride,
                         const size_t advance) {
  const size_t pos_limit = (iter + 1) * stride;
  size_t ipos = pos;
  while (ipos < pos_limit) {
    const int rel_pos = ipos - (iter * stride);
    if (rel_pos < count) {
      const size_t drp_idx = dest_start_idx + rel_pos;
      const size_t orp_idx = orig_start_idx + rel_pos;

      // Detect the presence of a signed integer, fixed-precision representation in the destination
      // by a nonzero number of scaling bits.  The result will emerge as, at most, the low 64 bits
      // of the rescaled 95-bit quantity.
      if (dest_scale_bits > 0) {
        const int95_t oval = { __ldlu(&orig_crd[orp_idx]), __ldlu(&orig_crd_ovrf[orp_idx]) };
        const int95_t dval = changeFPBits(oval, orig_scale_bits, dest_scale_bits);
        __stwt(&dest_crd[drp_idx], dval.x);
      }
      else {
        const double dval = int95ToDouble(__ldlu(&orig_crd[orp_idx]),
                                          __ldlu(&orig_crd_ovrf[orp_idx])) / orig_scale;
        __stwt(&dest_crd[drp_idx], dval);
      }
    }
    ipos += advance;
  }
  return ipos;  
}
/// \}

/// \brief Kernels for transferring groups of three X, Y, and Z coordinate arrays.  Descriptions of
///        parameters follow from kCopyBoxInformation() above, with the addition of:
///
/// \param xdest             Cartesian X coordinates of the destination object
/// \param ydest             Cartesian Y coordinates of the destination object
/// \param zdest             Cartesian Z coordinates of the destination object
/// \param xorig             Cartesian X coordinates of the origin object
/// \param yorig             Cartesian Y coordinates of the origin object
/// \param zorig             Cartesian Z coordinates of the origin object
/// \param dest_atom_offset  Index at which to start writing coordinates in the destination arrays
/// \param orig_atom_offset  Index at which to start reading coordinates in the origin arrays
/// \{
template <typename Tdest, typename Torig>
__global__ void __launch_bounds__(medium_block_size, 2)
kCopyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, const double dest_scale,
                   const int dest_scale_bits, const Torig* xorig, const Torig* yorig,
                   const Torig* zorig, const double orig_scale, const int orig_scale_bits,
                   const int natom, double* umat_dest, double* invu_dest, double* bdim_dest,
                   const double* umat_orig, const double* invu_orig, const double* bdim_orig,
                   const int dest_atom_offset, const int orig_atom_offset,
                   const int dest_xfrm_offset, const int orig_xfrm_offset,
                   const int dest_bdim_offset, const int orig_bdim_offset) {

  // Copy the box information, if provided.  This assumes that there are at least eight warps in
  // the block, which is true for all known architectures.
  if (umat_dest != nullptr && blockIdx.x == 0) {
    kCopyBoxInformation(umat_dest, invu_dest, bdim_dest, umat_orig, invu_orig, bdim_orig,
                        dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset, 5);
  }

  // Copy each of three coordinate arrays
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t stride = devcRoundUp(natom, warp_size_int);
  const size_t advance = (blockDim.x * gridDim.x);
  pos = copyCoordinateSet<Tdest, Torig>(xdest, dest_atom_offset, dest_scale, dest_scale_bits,
                                        xorig, orig_atom_offset, orig_scale, orig_scale_bits,
                                        natom, pos, 0, stride, advance);
  pos = copyCoordinateSet<Tdest, Torig>(ydest, dest_atom_offset, dest_scale, dest_scale_bits,
                                        yorig, orig_atom_offset, orig_scale, orig_scale_bits,
                                        natom, pos, 1, stride, advance);
  pos = copyCoordinateSet<Tdest, Torig>(zdest, dest_atom_offset, dest_scale, dest_scale_bits,
                                        zorig, orig_atom_offset, orig_scale, orig_scale_bits,
                                        natom, pos, 2, stride, advance);
}

template <typename Tdest>
__global__ void __launch_bounds__(medium_block_size, 2)
kCopyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, const double dest_scale,
                   const int dest_scale_bits, const llint* xorig, const llint* yorig,
                   const llint* zorig, const int* xorig_ovrf, const int* yorig_ovrf,
                   const int* zorig_ovrf, const double orig_scale, const int orig_scale_bits,
                   const int natom, double* umat_dest, double* invu_dest, double* bdim_dest,
                   const double* umat_orig, const double* invu_orig, const double* bdim_orig,
                   const int dest_atom_offset, const int orig_atom_offset,
                   const int dest_xfrm_offset, const int orig_xfrm_offset,
                   const int dest_bdim_offset, const int orig_bdim_offset) {

  // Copy the box information, if provided.  This assumes that there are at least eight warps in
  // the block, which is true for all known architectures.
  if (umat_dest != nullptr && blockIdx.x == 0) {
    kCopyBoxInformation(umat_dest, invu_dest, bdim_dest, umat_orig, invu_orig, bdim_orig,
                        dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset, 5);
  }

  // Copy each of three coordinate arrays
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t stride = devcRoundUp(natom, warp_size_int);
  const size_t advance = (blockDim.x * gridDim.x);
  pos = copyCoordinateSet<Tdest>(xdest, dest_atom_offset, dest_scale_bits, xorig, xorig_ovrf,
                                 orig_atom_offset, orig_scale, orig_scale_bits, natom, pos, 0,
                                 stride, advance);
  pos = copyCoordinateSet<Tdest>(ydest, dest_atom_offset, dest_scale_bits, yorig, yorig_ovrf,
                                 orig_atom_offset, orig_scale, orig_scale_bits, natom, pos, 1,
                                 stride, advance);
  pos = copyCoordinateSet<Tdest>(zdest, dest_atom_offset, dest_scale_bits, zorig, zorig_ovrf,
                                 orig_atom_offset, orig_scale, orig_scale_bits, natom, pos, 2,
                                 stride, advance);
}

template <typename Torig>
__global__ void __launch_bounds__(medium_block_size, 2)
kCopyCoordinateXYZ(llint* xdest, llint* ydest, llint* zdest, int* xdest_ovrf, int* ydest_ovrf,
                   int* zdest_ovrf, const double dest_scale, const int dest_scale_bits,
                   const Torig* xorig, const Torig* yorig, const Torig* zorig,
                   const double orig_scale, const int orig_scale_bits, const int natom,
                   double* umat_dest, double* invu_dest, double* bdim_dest,
                   const double* umat_orig, const double* invu_orig, const double* bdim_orig,
                   llint* boxvecs_dest, int* boxvec_ovrf_dest, const int dest_atom_offset,
                   const int orig_atom_offset, const int dest_xfrm_offset,
                   const int orig_xfrm_offset, const int dest_bdim_offset,
                   const int orig_bdim_offset) {

  // Copy the box information, if provided.  This assumes that there are at least eight warps in
  // the block, which is true for all known architectures.
  if (umat_dest != nullptr && blockIdx.x == 0) {
    kCopyBoxInformation(umat_dest, invu_dest, bdim_dest, umat_orig, invu_orig, bdim_orig,
                        dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset, 5);
  }
  
  // Copy each of three coordinate arrays
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t stride = devcRoundUp(natom, warp_size_int);
  const size_t advance = (blockDim.x * gridDim.x);
  pos = copyCoordinateSet<Torig>(xdest, xdest_ovrf, dest_atom_offset, dest_scale, dest_scale_bits,
                                 xorig, orig_atom_offset, orig_scale_bits, natom, pos, 0, stride,
                                 advance);
  pos = copyCoordinateSet<Torig>(ydest, ydest_ovrf, dest_atom_offset, dest_scale, dest_scale_bits,
                                 yorig, orig_atom_offset, orig_scale_bits, natom, pos, 1, stride,
                                 advance);
  pos = copyCoordinateSet<Torig>(zdest, zdest_ovrf, dest_atom_offset, dest_scale, dest_scale_bits,
                                 zorig, orig_atom_offset, orig_scale_bits, natom, pos, 2, stride,
                                 advance);
}
/// \}

/// \brief Unroll the second layer of the double-templated coordinate copy (X, Y, Z) operation.
///
/// \param xdest             Cartesian X coordinates for all particles in the destination system(s)
/// \param ydest             Cartesian Y coordinates for all particles in the destination system(s)
/// \param zdest             Cartesian Z coordinates for all particles in the destination system(s)
/// \param dest_scale        Scaling factor for the destination coordinates
/// \param ct_dest           Specifier code for the destination data type
/// \param xorig             Cartesian X coordinates for all particles in the original system(s)
/// \param yorig             Cartesian Y coordinates for all particles in the original system(s)
/// \param zorig             Cartesian Z coordinates for all particles in the original system(s)
/// \param orig_scale        Scaling factor for the origin coordinates
/// \param ct_orig           Specifier code for the origin data type
/// \param natom             The number of atoms' coordinates to copy
/// \param umat_dest         Transformation matrix taking real coordinates into fractional space
/// \param invu_dest         Transformation matrix taking fractional coordinates into real space
/// \param bdim_dest         The copy target box dimensions
/// \param umat_orig         Transformation matrix taking real coordinates into fractional space
/// \param invu_orig         Transformation matrix taking fractional coordinates into real space
/// \param bdim_orig         The origin's box dimensions
/// \param dest_atom_offset  Offset at which to begin writing coordinates for each atom to the
///                          destination
/// \param orig_atom_offset  Offset at which to begin reading coordinates for each atom from the
///                          origin
/// \param dest_xfrm_offset  Offset for reading transformation matrices in the destination
/// \param orig_xfrm_offset  Offset for reading transformation matrices in the origin
/// \param dest_bdim_offset  Offset for reading box dimensions in the destination
/// \param orig_bdim_offset  Offset for reading box dimensions in the origin
/// \param gpu               Details of the GPU to use
template <typename Tdest>
void unrollCCXYZOrigin(Tdest* xdest, Tdest* ydest, Tdest* zdest, const double dest_scale,
                       const int dest_scale_bits, const void* xorig, const void* yorig,
                       const void* zorig, const double orig_scale, const size_t ct_orig,
                       const int natom, double* umat_dest, double* invu_dest, double* bdim_dest,
                       const double* umat_orig, const double* invu_orig, const double* bdim_orig,
                       const int dest_atom_offset, const int orig_atom_offset,
                       const int dest_xfrm_offset, const int orig_xfrm_offset,
                       const int dest_bdim_offset, const int orig_bdim_offset,
                       const GpuDetails &gpu) {
  const int nblock = 2 * gpu.getSMPCount();
  if (ct_orig == double_type_index) {
    const double* d_xorig = reinterpret_cast<const double*>(xorig);
    const double* d_yorig = reinterpret_cast<const double*>(yorig);
    const double* d_zorig = reinterpret_cast<const double*>(zorig);
    kCopyCoordinateXYZ<<<nblock, medium_block_size>>>(xdest, ydest, zdest, dest_scale,
                                                      dest_scale_bits, d_xorig, d_yorig, d_zorig,
                                                      orig_scale, 0, natom, umat_dest, invu_dest,
                                                      bdim_dest, umat_orig, invu_orig, bdim_orig,
                                                      dest_atom_offset, orig_atom_offset,
                                                      dest_xfrm_offset, orig_xfrm_offset,
                                                      dest_bdim_offset, orig_bdim_offset);
  }
  else if (ct_orig == float_type_index) {
    const float* f_xorig = reinterpret_cast<const float*>(xorig);
    const float* f_yorig = reinterpret_cast<const float*>(yorig);
    const float* f_zorig = reinterpret_cast<const float*>(zorig);
    kCopyCoordinateXYZ<<<nblock, medium_block_size>>>(xdest, ydest, zdest, dest_scale,
                                                      dest_scale_bits, f_xorig, f_yorig, f_zorig,
                                                      orig_scale, 0, natom, umat_dest, invu_dest,
                                                      bdim_dest, umat_orig, invu_orig, bdim_orig,
                                                      dest_atom_offset, orig_atom_offset,
                                                      dest_xfrm_offset, orig_xfrm_offset,
                                                      dest_bdim_offset, orig_bdim_offset);
  }
  else if (ct_orig == short_type_index) {
    const short int* shi_xorig = reinterpret_cast<const short int*>(xorig);
    const short int* shi_yorig = reinterpret_cast<const short int*>(yorig);
    const short int* shi_zorig = reinterpret_cast<const short int*>(zorig);
    const int nbits = round(log(orig_scale));
    kCopyCoordinateXYZ<<<nblock, medium_block_size>>>(xdest, ydest, zdest, dest_scale,
                                                      dest_scale_bits, shi_xorig, shi_yorig,
                                                      shi_zorig, orig_scale, nbits, natom,
                                                      umat_dest, invu_dest, bdim_dest, umat_orig,
                                                      invu_orig, bdim_orig, dest_atom_offset,
                                                      orig_atom_offset, dest_xfrm_offset,
                                                      orig_xfrm_offset, dest_bdim_offset,
                                                      orig_bdim_offset);
  }
  else if (ct_orig == int_type_index) {
    const int* i_xorig = reinterpret_cast<const int*>(xorig);
    const int* i_yorig = reinterpret_cast<const int*>(yorig);
    const int* i_zorig = reinterpret_cast<const int*>(zorig);
    const int nbits = round(log(orig_scale));
    kCopyCoordinateXYZ<<<nblock, medium_block_size>>>(xdest, ydest, zdest, dest_scale,
                                                      dest_scale_bits, i_xorig, i_yorig, i_zorig,
                                                      orig_scale, nbits, natom, umat_dest,
                                                      invu_dest, bdim_dest, umat_orig, invu_orig,
                                                      bdim_orig, dest_atom_offset,
                                                      orig_atom_offset, dest_xfrm_offset,
                                                      orig_xfrm_offset, dest_bdim_offset,
                                                      orig_bdim_offset);
  }
  else if (ct_orig == llint_type_index) {
    const llint* lli_xorig = reinterpret_cast<const llint*>(xorig);
    const llint* lli_yorig = reinterpret_cast<const llint*>(yorig);
    const llint* lli_zorig = reinterpret_cast<const llint*>(zorig);
    const int nbits = round(log(orig_scale));
    kCopyCoordinateXYZ<<<nblock, medium_block_size>>>(xdest, ydest, zdest, dest_scale,
                                                      dest_scale_bits, lli_xorig, lli_yorig,
                                                      lli_zorig, orig_scale, nbits, natom,
                                                      umat_dest, invu_dest, bdim_dest, umat_orig,
                                                      invu_orig, bdim_orig, dest_atom_offset,
                                                      orig_atom_offset, dest_xfrm_offset,
                                                      orig_xfrm_offset, dest_bdim_offset,
                                                      orig_bdim_offset);
  }
}

/// \brief Copy multiple frames from one coordinate series into another.  Type reassignment for
///        each series will have been done on the CPU host, after passing their void-cast abstracts
///        through the C++ : HPC demarcation.
///
/// \param dest          Coordinate synthesis object into which information will be copied
/// \param orig          The series from which coordinates will be obtained
/// \param system_pairs  List of system pairs to transfer
/// \param copy_count    Trusted length of system_pairs
template <typename Tdest, typename Torig>
__global__ void __launch_bounds__(large_block_size, 1)
kMultiSystemCoordinateCopy(CoordinateSeriesWriter<Tdest> dest,
                           const CoordinateSeriesReader<Torig> orig, const int2* system_pairs,
                           const int copy_count) {

  // Each block will take a separate system and carry out the copying, with any necessary
  // fixed-precision bit conversions.
  int pair_idx = blockIdx.x;
  while (pair_idx < copy_count) {
    const int2 origx_to_desty = __ldca(&system_pairs[pair_idx]);

    // Skip copy requests to or from invalid system indices.  This error will be caught by the
    // CPU while the kernel is running and raise an error.
    if (origx_to_desty.x < 0 || origx_to_desty.x >= orig.nframe ||
        origx_to_desty.y < 0 || origx_to_desty.y >= dest.nframe) {
      pair_idx += gridDim.x;
      continue;
    }
    const int natom = dest.natom;

    // Skip copy requests between systems of different sizes.  These will be caught by the CPU
    // while the kernel is running and raise an error.
    if (natom != orig.natom) {
      pair_idx += gridDim.x;
      continue;
    }
    
    // Pre-compute offsets
    const int xfrm_w = devcRoundUp(9, warp_size_int);
    const int bdim_w = devcRoundUp(6, warp_size_int);
    const int orig_xfrm_offset = origx_to_desty.x * xfrm_w;
    const int dest_xfrm_offset = origx_to_desty.y * xfrm_w;
    const int orig_bdim_offset = origx_to_desty.x * bdim_w;
    const int dest_bdim_offset = origx_to_desty.y * bdim_w;
    const size_t padded_natom = devcRoundUp(natom, warp_size_int);
    const size_t dest_atom_offset = padded_natom * origx_to_desty.y;
    const size_t orig_atom_offset = padded_natom * origx_to_desty.x;

    // Coordinate series have box information, even if it is irrelevant for a coordinate series
    // storing something like a trajectory of velocities.  Copy the available data.
    kCopyBoxInformation(dest.umat, dest.invu, dest.boxdim, orig.umat, orig.invu, orig.boxdim,
                        dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset, 5);
    size_t pos = threadIdx.x;
    pos = copyCoordinateSet<Tdest, Torig>(dest.xcrd, dest_atom_offset, dest.gpos_scale,
                                          dest.gpos_bits, orig.xcrd, orig_atom_offset,
                                          orig.gpos_scale, orig.gpos_bits, natom, pos, 0,
                                          padded_natom, blockDim.x);
    pos = copyCoordinateSet<Tdest, Torig>(dest.ycrd, dest_atom_offset, dest.gpos_scale,
                                          dest.gpos_bits, orig.ycrd, orig_atom_offset,
                                          orig.gpos_scale, orig.gpos_bits, natom, pos, 1,
                                          padded_natom, blockDim.x);
    pos = copyCoordinateSet<Tdest, Torig>(dest.zcrd, dest_atom_offset, dest.gpos_scale,
                                          dest.gpos_bits, orig.zcrd, orig_atom_offset,
                                          orig.gpos_scale, orig.gpos_bits, natom, pos, 2,
                                          padded_natom, blockDim.x);

    // Increment the system copy counter
    pair_idx += gridDim.x;
  }
}

/// \brief Copy multiple systems from one coordinate series to another.  This unrolls the inner
///        branch over possible CoordinateSeries data types.
///
/// \param destination       Abstract of the series to receive coordiante information, with its
///                          type restored from <void> by the calling function
/// \param origin            Void-casted form of the abstract for the series from which information
///                          will be derived
/// \param ct_orig           Codified data type of the origin series 
/// \param system_pairs      List of frames to copy between each series, with the origin frame
///                          given in the "x" member of each tuple and the destination frame given
///                          in the "y" member
/// \param copy_count        Trusted length of system_pairs
/// \param gpu               Details of the GPU to use
template <typename Tdest>
void unrollCCXYZOrigin(CoordinateSeriesWriter<Tdest> *destination,
                       const CoordinateSeriesReader<void> &origin, size_t ct_orig,
                       const int2* system_pairs, int copy_count, const GpuDetails &gpu) {
  const int nblock = gpu.getSMPCount();
  if (ct_orig == double_type_index) {
    const CoordinateSeriesReader<double> tr_orig = restoreType<double>(origin);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(*destination, tr_orig, system_pairs,
                                                             copy_count);
  }
  else if (ct_orig == float_type_index) {
    const CoordinateSeriesReader<float> tr_orig = restoreType<float>(origin);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(*destination, tr_orig, system_pairs,
                                                             copy_count);
  }
  else if (ct_orig == short_type_index) {
    const CoordinateSeriesReader<short int> tr_orig = restoreType<short int>(origin);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(*destination, tr_orig, system_pairs,
                                                             copy_count);
  }
  else if (ct_orig == int_type_index) {
    const CoordinateSeriesReader<int> tr_orig = restoreType<int>(origin);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(*destination, tr_orig, system_pairs,
                                                             copy_count);
  }
  else if (ct_orig == llint_type_index) {
    const CoordinateSeriesReader<llint> tr_orig = restoreType<llint>(origin);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(*destination, tr_orig, system_pairs,
                                                             copy_count);
  }
}

/// \brief Copy multiple systems from a PhaseSpaceSynthesis into frames of a coordinate series.  As
///        with other templated coordinate copy kernels, the type reassignment will be done on the
///        CPU, by the launching function.
///
/// \param dest          Coordinate synthesis object into which information will be copied
/// \param kind          Specify whether to copy into positions, velocities, or forces
/// \param orig          The series from which coordinates will be obtained
/// \param system_pairs  List of system pairs to transfer
/// \param copy_count    Trusted length of system_pairs
template <typename Tdest>
__global__ void __launch_bounds__(large_block_size, 1)
kMultiSystemCoordinateCopy(CoordinateSeriesWriter<Tdest> dest, const PsSynthesisReader orig,
                           const TrajectoryKind kind, const int2* system_pairs,
                           const int copy_count) {

  // Each block will take a separate system and carry out the copying, with any necessary
  // fixed-precision bit conversions.
  int pair_idx = blockIdx.x;
  while (pair_idx < copy_count) {
    const int2 origx_to_desty = __ldca(&system_pairs[pair_idx]);

    // Skip copy requests to or from invalid system indices.  This error will be caught by the
    // CPU while the kernel is running and raise an error.
    if (origx_to_desty.x < 0 || origx_to_desty.x >= orig.system_count ||
        origx_to_desty.y < 0 || origx_to_desty.y >= dest.nframe) {
      pair_idx += gridDim.x;
      continue;
    }
    const int natom = dest.natom;

    // Skip copy requests between systems of different sizes.  These will be caught by the CPU
    // while the kernel is running and raise an error.
    if (natom != __ldca(&orig.atom_counts[origx_to_desty.x])) {
      pair_idx += gridDim.x;
      continue;
    }

    // Execute the switch over coordinate kinds
    const size_t padded_natom = devcRoundUp(natom, warp_size_int);
    const size_t orig_atom_offset = __ldca(&orig.atom_starts[origx_to_desty.x]);
    const size_t dest_atom_offset = padded_natom * origx_to_desty.y;
    size_t pos = threadIdx.x;
    switch (kind) {
    case TrajectoryKind::POSITIONS:

      // Copy the box information
      {
        const int xfrm_w = devcRoundUp(9, warp_size_int);
        const int bdim_w = devcRoundUp(6, warp_size_int);
        const int orig_xfrm_offset = origx_to_desty.x * xfrm_w;
        const int dest_xfrm_offset = origx_to_desty.y * xfrm_w;
        const int orig_bdim_offset = origx_to_desty.x * bdim_w;
        const int dest_bdim_offset = origx_to_desty.y * bdim_w;
        kCopyBoxInformation(dest.umat, dest.invu, dest.boxdim, orig.umat, orig.invu, orig.boxdims,
                            dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset,
                            0);

        // Copy into each of three position arrays
        pos = copyCoordinateSet<Tdest>(dest.xcrd, dest_atom_offset, dest.gpos_bits, orig.xcrd,
                                       orig.xcrd_ovrf, orig_atom_offset, orig.gpos_scale,
                                       orig.gpos_bits, natom, pos, 0, padded_natom, blockDim.x);
        pos = copyCoordinateSet<Tdest>(dest.ycrd, dest_atom_offset, dest.gpos_bits, orig.ycrd,
                                       orig.ycrd_ovrf, orig_atom_offset, orig.gpos_scale,
                                       orig.gpos_bits, natom, pos, 1, padded_natom, blockDim.x);
        pos = copyCoordinateSet<Tdest>(dest.zcrd, dest_atom_offset, dest.gpos_bits, orig.zcrd,
                                       orig.zcrd_ovrf, orig_atom_offset, orig.gpos_scale,
                                       orig.gpos_bits, natom, pos, 2, padded_natom, blockDim.x);
      }
      break;
    case TrajectoryKind::VELOCITIES:

      // Copy into each of the three velocity arrays
      pos = copyCoordinateSet<Tdest>(dest.xcrd, dest_atom_offset, dest.gpos_bits, orig.xvel,
                                     orig.xvel_ovrf, orig_atom_offset, orig.vel_scale,
                                     orig.vel_bits, natom, pos, 0, padded_natom, blockDim.x);
      pos = copyCoordinateSet<Tdest>(dest.ycrd, dest_atom_offset, dest.gpos_bits, orig.yvel,
                                     orig.yvel_ovrf, orig_atom_offset, orig.vel_scale,
                                     orig.vel_bits, natom, pos, 1, padded_natom, blockDim.x);
      pos = copyCoordinateSet<Tdest>(dest.zcrd, dest_atom_offset, dest.gpos_bits, orig.zvel,
                                     orig.zvel_ovrf, orig_atom_offset, orig.vel_scale,
                                     orig.vel_bits, natom, pos, 2, padded_natom, blockDim.x);
      break;
    case TrajectoryKind::FORCES:

      // Copy into each of the three force arrays
      pos = copyCoordinateSet<Tdest>(dest.xcrd, dest_atom_offset, dest.gpos_bits, orig.xfrc,
                                     orig.xfrc_ovrf, orig_atom_offset, orig.frc_scale,
                                     orig.frc_bits, natom, pos, 0, padded_natom, blockDim.x);
      pos = copyCoordinateSet<Tdest>(dest.ycrd, dest_atom_offset, dest.gpos_bits, orig.yfrc,
                                     orig.yfrc_ovrf, orig_atom_offset, orig.frc_scale,
                                     orig.frc_bits, natom, pos, 1, padded_natom, blockDim.x);
      pos = copyCoordinateSet<Tdest>(dest.zcrd, dest_atom_offset, dest.gpos_bits, orig.zfrc,
                                     orig.zfrc_ovrf, orig_atom_offset, orig.frc_scale,
                                     orig.frc_bits, natom, pos, 2, padded_natom, blockDim.x);
      break;
    }

    // Increment the system copy counter
    pair_idx += gridDim.x;
  }
}

/// \brief Copy multiple systems from a PhaseSpaceSynthesis into frames of a coordinate series.  As
///        with other templated coordinate copy kernels, the type reassignment will be done on the
///        CPU, by the launching function.  Descriptions of input parameters follow from a previous
///        overloaded variant of kMultiSystemCoordinateCopy() above.
template <typename Tdest>
__global__ void __launch_bounds__(large_block_size, 1)
kMultiSystemCoordinateCopy(CoordinateSeriesWriter<Tdest> dest, const CondensateReader orig,
                           const int2* system_pairs, const int copy_count) {

  // Each block will take a separate system and carry out the copying, with any necessary
  // fixed-precision bit conversions.
  int pair_idx = blockIdx.x;
  while (pair_idx < copy_count) {
    const int2 origx_to_desty = __ldca(&system_pairs[pair_idx]);

    // Skip copy requests to or from invalid system indices.  This error will be caught by the
    // CPU while the kernel is running and raise an error.
    if (origx_to_desty.x < 0 || origx_to_desty.x >= orig.system_count ||
        origx_to_desty.y < 0 || origx_to_desty.y >= dest.nframe) {
      pair_idx += gridDim.x;
      continue;
    }
    const int natom = dest.natom;

    // Skip copy requests between systems of different sizes.  These will be caught by the CPU
    // while the kernel is running and raise an error.
    if (natom != __ldca(&orig.atom_counts[origx_to_desty.x])) {
      pair_idx += gridDim.x;
      continue;
    }

    // Assume that the box information is relevant (a position coordinate transfer)
    const int xfrm_w = devcRoundUp(9, warp_size_int);
    const int bdim_w = devcRoundUp(6, warp_size_int);
    const int orig_xfrm_offset = origx_to_desty.x * xfrm_w;
    const int dest_xfrm_offset = origx_to_desty.y * xfrm_w;
    const int orig_bdim_offset = origx_to_desty.x * bdim_w;
    const int dest_bdim_offset = origx_to_desty.y * bdim_w;
    kCopyBoxInformation(dest.umat, dest.invu, dest.boxdim, orig.umat, orig.invu, orig.boxdims,
                        dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset,
                        0);

    // Pre-compute array offsets and execute the switch over the Condensate's data types.
    const size_t padded_natom = devcRoundUp(natom, warp_size_int);
    const size_t orig_atom_offset = __ldca(&orig.atom_starts[origx_to_desty.x]);
    const size_t dest_atom_offset = padded_natom * origx_to_desty.y;
    size_t pos = threadIdx.x;
    switch (orig.mode) {
    case PrecisionModel::DOUBLE:
      pos = copyCoordinateSet<Tdest, double>(dest.xcrd, dest_atom_offset, dest.gpos_scale,
                                             dest.gpos_bits, orig.xcrd, orig_atom_offset, 1.0, 0,
                                             natom, pos, 0, padded_natom, blockDim.x);
      pos = copyCoordinateSet<Tdest, double>(dest.ycrd, dest_atom_offset, dest.gpos_scale,
                                             dest.gpos_bits, orig.ycrd, orig_atom_offset, 1.0, 0,
                                             natom, pos, 1, padded_natom, blockDim.x);
      pos = copyCoordinateSet<Tdest, double>(dest.zcrd, dest_atom_offset, dest.gpos_scale,
                                             dest.gpos_bits, orig.zcrd, orig_atom_offset, 1.0, 0,
                                             natom, pos, 2, padded_natom, blockDim.x);
      break;
    case PrecisionModel::SINGLE:
      pos = copyCoordinateSet<Tdest, float>(dest.xcrd, dest_atom_offset, dest.gpos_scale,
                                            dest.gpos_bits, orig.xcrd_sp, orig_atom_offset, 1.0, 0,
                                            natom, pos, 0, padded_natom, blockDim.x);
      pos = copyCoordinateSet<Tdest, float>(dest.ycrd, dest_atom_offset, dest.gpos_scale,
                                            dest.gpos_bits, orig.ycrd_sp, orig_atom_offset, 1.0, 0,
                                            natom, pos, 1, padded_natom, blockDim.x);
      pos = copyCoordinateSet<Tdest, float>(dest.zcrd, dest_atom_offset, dest.gpos_scale,
                                            dest.gpos_bits, orig.zcrd_sp, orig_atom_offset, 1.0, 0,
                                            natom, pos, 2, padded_natom, blockDim.x);
      break;
    }

    // Increment the system copy counter
    pair_idx += gridDim.x;
  }
}
    
/// \brief Copy multiple frames from a coordinate series into systems of a PhaseSpaceSynthesis.  As
///        with other templated coordinate copy kernels, the type reassignment will be done on the
///        CPU, by the launching function.  Descriptions of input parameters follow from a previous
///        overloaded variant of kMultiSystemCoordinateCopy() above, with the exception that the
///        origin and destination have reversed their object types.
template <typename Torig>
__global__ void __launch_bounds__(large_block_size, 1)
kMultiSystemCoordinateCopy(PsSynthesisWriter dest, const TrajectoryKind kind,
                           const CoordinateSeriesReader<Torig> orig, const int2* system_pairs,
                           const int copy_count) {
  int pair_idx = blockIdx.x;
  while (pair_idx < copy_count) {
    const int2 origx_to_desty = __ldca(&system_pairs[pair_idx]);
    if (origx_to_desty.x < 0 || origx_to_desty.x >= orig.nframe ||
        origx_to_desty.y < 0 || origx_to_desty.y >= dest.system_count) {
      pair_idx += gridDim.x;
      continue;
    }
    const int natom = __ldca(&dest.atom_counts[origx_to_desty.y]);
    if (natom != orig.natom) {
      pair_idx += gridDim.x;
      continue;
    }

    // Execute the switch over coordinate kinds
    const int padded_natom = devcRoundUp(natom, warp_size_int);
    const int orig_atom_offset = padded_natom * origx_to_desty.x;
    const int dest_atom_offset = __ldca(&dest.atom_starts[origx_to_desty.y]);
    size_t pos = threadIdx.x;
    switch (kind) {
    case TrajectoryKind::POSITIONS:

      // Copy the box information
      {
        const int xfrm_w = devcRoundUp(9, warp_size_int);
        const int bdim_w = devcRoundUp(6, warp_size_int);
        const int orig_xfrm_offset = origx_to_desty.x * xfrm_w;
        const int dest_xfrm_offset = origx_to_desty.y * xfrm_w;
        const int orig_bdim_offset = origx_to_desty.x * bdim_w;
        const int dest_bdim_offset = origx_to_desty.y * bdim_w;
        kCopyBoxInformation(dest.umat, dest.invu, dest.boxdims, orig.umat, orig.invu, orig.boxdim,
                            dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset,
                            0, dest.boxvecs, dest.boxvec_ovrf, dest.gpos_scale);

        // Copy into each of three position arrays
        pos = copyCoordinateSet<Torig>(dest.xcrd, dest.xcrd_ovrf, dest_atom_offset,
                                       dest.gpos_scale, dest.gpos_bits, orig.xcrd,
                                       orig_atom_offset, orig.gpos_bits, natom, pos, 0,
                                       padded_natom, blockDim.x);
        pos = copyCoordinateSet<Torig>(dest.ycrd, dest.ycrd_ovrf, dest_atom_offset,
                                       dest.gpos_scale, dest.gpos_bits, orig.ycrd,
                                       orig_atom_offset, orig.gpos_bits, natom, pos, 1,
                                       padded_natom, blockDim.x);
        pos = copyCoordinateSet<Torig>(dest.zcrd, dest.zcrd_ovrf, dest_atom_offset,
                                       dest.gpos_scale, dest.gpos_bits, orig.zcrd,
                                       orig_atom_offset, orig.gpos_bits, natom, pos, 2,
                                       padded_natom, blockDim.x);
      }
      break;
    case TrajectoryKind::VELOCITIES:

      // Copy into each of the three velocity arrays
      pos = copyCoordinateSet<Torig>(dest.xvel, dest.xvel_ovrf, dest_atom_offset, dest.vel_scale,
                                     dest.vel_bits, orig.xcrd, orig_atom_offset, orig.gpos_bits,
                                     natom, pos, 0, padded_natom, blockDim.x);
      pos = copyCoordinateSet<Torig>(dest.yvel, dest.yvel_ovrf, dest_atom_offset, dest.vel_scale,
                                     dest.vel_bits, orig.ycrd, orig_atom_offset, orig.gpos_bits,
                                     natom, pos, 1, padded_natom, blockDim.x);
      pos = copyCoordinateSet<Torig>(dest.zvel, dest.zvel_ovrf, dest_atom_offset, dest.vel_scale,
                                     dest.vel_bits, orig.zcrd, orig_atom_offset, orig.gpos_bits,
                                     natom, pos, 2, padded_natom, blockDim.x);
      break;
    case TrajectoryKind::FORCES:

      // Copy into each of the three force arrays
      pos = copyCoordinateSet<Torig>(dest.xfrc, dest.xfrc_ovrf, dest_atom_offset, dest.frc_scale,
                                     dest.frc_bits, orig.xcrd, orig_atom_offset, orig.gpos_bits,
                                     natom, pos, 0, padded_natom, blockDim.x);
      pos = copyCoordinateSet<Torig>(dest.yfrc, dest.yfrc_ovrf, dest_atom_offset, dest.frc_scale,
                                     dest.frc_bits, orig.ycrd, orig_atom_offset, orig.gpos_bits,
                                     natom, pos, 1, padded_natom, blockDim.x);
      pos = copyCoordinateSet<Torig>(dest.zfrc, dest.zfrc_ovrf, dest_atom_offset, dest.frc_scale,
                                     dest.frc_bits, orig.zcrd, orig_atom_offset, orig.gpos_bits,
                                     natom, pos, 2, padded_natom, blockDim.x);
      break;
    }
    pair_idx += gridDim.x;
  }
}

/// \brief Copy multiple frames from a coordinate series into systems of a Condensate object.  As
///        with other templated coordinate copy kernels, the type reassignment will be done on the
///        CPU, by the launching function.  Descriptions of input parameters follow from a previous
///        overloaded variant of kMultiSystemCoordinateCopy() above, with the exception that the
///        origin and destination have reversed their object types.
template <typename Torig>
__global__ void __launch_bounds__(large_block_size, 1)
kMultiSystemCoordinateCopy(CondensateWriter dest, const CoordinateSeriesReader<Torig> orig,
                           const int2* system_pairs, const int copy_count) {
  int pair_idx = blockIdx.x;
  while (pair_idx < copy_count) {
    const int2 origx_to_desty = __ldca(&system_pairs[pair_idx]);
    if (origx_to_desty.x < 0 || origx_to_desty.x >= orig.nframe ||
        origx_to_desty.y < 0 || origx_to_desty.y >= dest.system_count) {
      pair_idx += gridDim.x;
      continue;
    }
    const int natom = __ldca(&dest.atom_counts[origx_to_desty.y]);
    if (natom != orig.natom) {
      pair_idx += gridDim.x;
      continue;
    }

    // Assume that the box information is relevant (a position coordinate transfer)
    const int xfrm_w = devcRoundUp(9, warp_size_int);
    const int bdim_w = devcRoundUp(6, warp_size_int);
    const int orig_xfrm_offset = origx_to_desty.x * xfrm_w;
    const int dest_xfrm_offset = origx_to_desty.y * xfrm_w;
    const int orig_bdim_offset = origx_to_desty.x * bdim_w;
    const int dest_bdim_offset = origx_to_desty.y * bdim_w;
    kCopyBoxInformation(dest.umat, dest.invu, dest.boxdims, orig.umat, orig.invu, orig.boxdim,
                        dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset,
                        0);

    // Pre-compute array offsets and execute the switch over the Condensate's data types.
    const size_t padded_natom = devcRoundUp(natom, warp_size_int);
    const size_t dest_atom_offset = __ldca(&dest.atom_starts[origx_to_desty.y]);
    const size_t orig_atom_offset = padded_natom * (size_t)(origx_to_desty.x);
    size_t pos = threadIdx.x;
    switch (dest.mode) {
    case PrecisionModel::DOUBLE:
      pos = copyCoordinateSet<double, Torig>(dest.xcrd, dest_atom_offset, 1.0, 0, orig.xcrd,
                                             orig_atom_offset, orig.gpos_scale, orig.gpos_bits,
                                             natom, pos, 0, padded_natom, blockDim.x);
      pos = copyCoordinateSet<double, Torig>(dest.ycrd, dest_atom_offset, 1.0, 0, orig.ycrd,
                                             orig_atom_offset, orig.gpos_scale, orig.gpos_bits,
                                             natom, pos, 1, padded_natom, blockDim.x);
      pos = copyCoordinateSet<double, Torig>(dest.zcrd, dest_atom_offset, 1.0, 0, orig.zcrd,
                                             orig_atom_offset, orig.gpos_scale, orig.gpos_bits,
                                             natom, pos, 2, padded_natom, blockDim.x);
      break;
    case PrecisionModel::SINGLE:
      pos = copyCoordinateSet<float, Torig>(dest.xcrd_sp, dest_atom_offset, 1.0, 0, orig.xcrd,
                                            orig_atom_offset, orig.gpos_scale, orig.gpos_bits,
                                            natom, pos, 0, padded_natom, blockDim.x);
      pos = copyCoordinateSet<float, Torig>(dest.ycrd_sp, dest_atom_offset, 1.0, 0, orig.ycrd,
                                            orig_atom_offset, orig.gpos_scale, orig.gpos_bits,
                                            natom, pos, 0, padded_natom, blockDim.x);
      pos = copyCoordinateSet<float, Torig>(dest.zcrd_sp, dest_atom_offset, 1.0, 0, orig.zcrd,
                                            orig_atom_offset, orig.gpos_scale, orig.gpos_bits,
                                            natom, pos, 0, padded_natom, blockDim.x);
      break;
    }
    
    // Increment the system copy count
    pair_idx += gridDim.x;
  }
}

} // namespace trajectory
} // namespace stormm

#endif
