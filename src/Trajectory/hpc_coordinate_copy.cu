// -*-c++-*-
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "hpc_coordinate_copy.cuh"
#include "hpc_coordinate_copy.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
__device__ __forceinline__
void kCopyBoxInformation(double* umat_dest, double* invu_dest, double* bdim_dest,
                         const double* umat_orig, const double* invu_orig, const double* bdim_orig,
                         const int dest_xfrm_offset, const int orig_xfrm_offset,
                         const int dest_bdim_offset, const int orig_bdim_offset,
                         const int warp_offset, llint* dest_boxvecs, int* dest_boxvec_ovrf,
                         const double dest_scale, const llint* orig_boxvecs,
                         const int* orig_boxvec_ovrf, const int dest_bits, const int orig_bits) {
  const int warp_idx = (threadIdx.x >> warp_bits);
  if (warp_idx == warp_offset) {
    int pos = (threadIdx.x & warp_bits_mask_int);
    while (pos < 9) {
      __stwt(&umat_dest[dest_xfrm_offset + pos], __ldlu(&umat_orig[orig_xfrm_offset + pos]));
      pos += warp_size_int;
    }
  }
  else if (warp_idx == warp_offset + 1) {
    int pos = (threadIdx.x & warp_bits_mask_int);
    while (pos < 9) {
      const double dval = __ldlu(&invu_orig[orig_xfrm_offset + pos]);
      __stwt(&invu_dest[dest_xfrm_offset + pos], dval);
      if (dest_boxvecs != nullptr) {
        if (orig_boxvecs != nullptr) {
          const int95_t orig_fpval = { __ldlu(&orig_boxvecs[orig_xfrm_offset + pos]),
                                       __ldlu(&orig_boxvec_ovrf[orig_xfrm_offset + pos]) };
          const int95_t fpval = changeFPBits(orig_fpval, orig_bits, dest_bits);
          __stwt(&dest_boxvecs[dest_xfrm_offset + pos], fpval.x);
          __stwt(&dest_boxvec_ovrf[dest_xfrm_offset + pos], fpval.y);
        }
        else {
          const int95_t fpval = doubleToInt95(dval * dest_scale);
          __stwt(&dest_boxvecs[dest_xfrm_offset + pos], fpval.x);
          __stwt(&dest_boxvec_ovrf[dest_xfrm_offset + pos], fpval.y);
        }
      }
      pos += warp_size_int;
    }
  }
  else if (warp_idx == warp_offset + 2) {
    int pos = (threadIdx.x & warp_bits_mask_int);
    while (pos < 6) {
      __stwt(&bdim_dest[dest_bdim_offset + pos], __ldlu(&bdim_orig[orig_bdim_offset + pos]));
      pos += warp_size_int;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinateXYZBox(void* xdest, void* ydest, void* zdest,
                                       const double dest_scale, const size_t ct_dest,
                                       const void* xorig, const void* yorig, const void* zorig,
                                       const double orig_scale, const size_t ct_orig,
                                       const int natom, double* umat_dest, double* invu_dest,
                                       double* bdim_dest, const double* umat_orig,
                                       const double* invu_orig, const double* bdim_orig,
                                       const int dest_atom_offset, const int orig_atom_offset,
                                       const int dest_xfrm_offset, const int orig_xfrm_offset,
                                       const int dest_bdim_offset, const int orig_bdim_offset,
                                       const GpuDetails &gpu) {
  const int nblock = 2 * gpu.getSMPCount();
  if (ct_dest == double_type_index) {
    double* d_xdest = reinterpret_cast<double*>(xdest);
    double* d_ydest = reinterpret_cast<double*>(ydest);
    double* d_zdest = reinterpret_cast<double*>(zdest);
    unrollCCXYZOrigin(d_xdest, d_ydest, d_zdest, dest_scale, 0, xorig, yorig, zorig, orig_scale,
                      ct_orig, natom, umat_dest, invu_dest, bdim_dest, umat_orig, invu_orig,
                      bdim_orig, dest_atom_offset, orig_atom_offset, dest_xfrm_offset,
                      orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset, gpu);
  }
  else if (ct_dest == float_type_index) {
    float* f_xdest = reinterpret_cast<float*>(xdest);
    float* f_ydest = reinterpret_cast<float*>(ydest);
    float* f_zdest = reinterpret_cast<float*>(zdest);
    unrollCCXYZOrigin(f_xdest, f_ydest, f_zdest, dest_scale, 0, xorig, yorig, zorig, orig_scale,
                      ct_orig, natom, umat_dest, invu_dest, bdim_dest, umat_orig, invu_orig,
                      bdim_orig, dest_atom_offset, orig_atom_offset, dest_xfrm_offset,
                      orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset, gpu);
  }
  else if (ct_dest == short_type_index) {
    short int* shi_xdest = reinterpret_cast<short int*>(xdest);
    short int* shi_ydest = reinterpret_cast<short int*>(ydest);
    short int* shi_zdest = reinterpret_cast<short int*>(zdest);
    const int nbits = round(log(dest_scale));
    unrollCCXYZOrigin(shi_xdest, shi_ydest, shi_zdest, dest_scale, nbits, xorig, yorig, zorig,
                      orig_scale, ct_orig, natom, umat_dest, invu_dest, bdim_dest, umat_orig,
                      invu_orig, bdim_orig, dest_atom_offset, orig_atom_offset, dest_xfrm_offset,
                      orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset, gpu);
  }
  else if (ct_dest == int_type_index) {
    int* i_xdest = reinterpret_cast<int*>(xdest);
    int* i_ydest = reinterpret_cast<int*>(ydest);
    int* i_zdest = reinterpret_cast<int*>(zdest);
    const int nbits = round(log(dest_scale));
    unrollCCXYZOrigin(i_xdest, i_ydest, i_zdest, dest_scale, nbits, xorig, yorig, zorig,
                      orig_scale, ct_orig, natom, umat_dest, invu_dest, bdim_dest, umat_orig,
                      invu_orig, bdim_orig, dest_atom_offset, orig_atom_offset, dest_xfrm_offset,
                      orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset, gpu);
  }
  else if (ct_dest == llint_type_index) {
    llint* lli_xdest = reinterpret_cast<llint*>(xdest);
    llint* lli_ydest = reinterpret_cast<llint*>(ydest);
    llint* lli_zdest = reinterpret_cast<llint*>(zdest);
    const int nbits = round(log(dest_scale));
    unrollCCXYZOrigin(lli_xdest, lli_ydest, lli_zdest, dest_scale, nbits, xorig, yorig, zorig,
                      orig_scale, ct_orig, natom, umat_dest, invu_dest, bdim_dest, umat_orig,
                      invu_orig, bdim_orig, dest_atom_offset, orig_atom_offset, dest_xfrm_offset,
                      orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset, gpu);
  }
}
  
//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinateXYZ(void* xdest, void* ydest, void* zdest, const double dest_scale,
                                    const size_t ct_dest, const void* xorig, const void* yorig,
                                    const void* zorig, const double orig_scale,
                                    const size_t ct_orig, const int natom,
                                    const int dest_atom_offset, const int orig_atom_offset,
                                    const GpuDetails &gpu) {
  launchCopyCoordinateXYZBox(xdest, ydest, zdest, dest_scale, ct_dest, xorig, yorig, zorig,
                             orig_scale, ct_orig, natom, nullptr, nullptr, nullptr, nullptr,
                             nullptr, nullptr, dest_atom_offset, orig_atom_offset, 0, 0, 0, 0,
                             gpu);
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinateXYZBox(void* xdest, void* ydest, void* zdest,
                                       const double dest_scale, const int dest_scale_bits,
                                       const size_t ct_dest, const llint* xorig,
                                       const llint* yorig, const llint* zorig,
                                       const int* xorig_ovrf, const int* yorig_ovrf,
                                       const int* zorig_ovrf, const double orig_scale,
                                       const int orig_scale_bits, const int natom,
                                       double* umat_dest, double* invu_dest, double* bdim_dest,
                                       const double* umat_orig, const double* invu_orig,
                                       const double* bdim_orig, const int dest_atom_offset,
                                       const int orig_atom_offset, const int dest_xfrm_offset,
                                       const int orig_xfrm_offset, const int dest_bdim_offset,
                                       const int orig_bdim_offset, const GpuDetails &gpu) {
  const int nblock = 2 * gpu.getSMPCount();
  const int nthread = medium_block_size;
  if (ct_dest == double_type_index) {
    double* d_xdest = reinterpret_cast<double*>(xdest);
    double* d_ydest = reinterpret_cast<double*>(ydest);
    double* d_zdest = reinterpret_cast<double*>(zdest);
    kCopyCoordinateXYZ<<<nblock, nthread>>>(d_xdest, d_ydest, d_zdest, dest_scale, dest_scale_bits,
                                            xorig, yorig, zorig, xorig_ovrf, yorig_ovrf,
                                            zorig_ovrf, orig_scale, orig_scale_bits, natom,
                                            umat_dest, invu_dest, bdim_dest, umat_orig, invu_orig,
                                            bdim_orig, dest_atom_offset, orig_atom_offset,
                                            dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset,
                                            orig_bdim_offset);
  }
  else if (ct_dest == float_type_index) {
    float* f_xdest = reinterpret_cast<float*>(xdest);
    float* f_ydest = reinterpret_cast<float*>(ydest);
    float* f_zdest = reinterpret_cast<float*>(zdest);
    kCopyCoordinateXYZ<<<nblock, nthread>>>(f_xdest, f_ydest, f_zdest, dest_scale, dest_scale_bits,
                                            xorig, yorig, zorig, xorig_ovrf, yorig_ovrf,
                                            zorig_ovrf, orig_scale, orig_scale_bits, natom,
                                            umat_dest, invu_dest, bdim_dest, umat_orig, invu_orig,
                                            bdim_orig, dest_atom_offset, orig_atom_offset,
                                            dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset,
                                            orig_bdim_offset);
  }
  else if (ct_dest == short_type_index) {
    short int* shi_xdest = reinterpret_cast<short int*>(xdest);
    short int* shi_ydest = reinterpret_cast<short int*>(ydest);
    short int* shi_zdest = reinterpret_cast<short int*>(zdest);
    kCopyCoordinateXYZ<<<nblock, nthread>>>(shi_xdest, shi_ydest, shi_zdest, dest_scale,
                                            dest_scale_bits, xorig, yorig, zorig, xorig_ovrf,
                                            yorig_ovrf, zorig_ovrf, orig_scale, orig_scale_bits,
                                            natom, umat_dest, invu_dest, bdim_dest, umat_orig,
                                            invu_orig, bdim_orig, dest_atom_offset,
                                            orig_atom_offset, dest_xfrm_offset, orig_xfrm_offset,
                                            dest_bdim_offset, orig_bdim_offset);
  }
  else if (ct_dest == int_type_index) {
    int* i_xdest = reinterpret_cast<int*>(xdest);
    int* i_ydest = reinterpret_cast<int*>(ydest);
    int* i_zdest = reinterpret_cast<int*>(zdest);
    kCopyCoordinateXYZ<<<nblock, nthread>>>(i_xdest, i_ydest, i_zdest, dest_scale, dest_scale_bits,
                                            xorig, yorig, zorig, xorig_ovrf, yorig_ovrf,
                                            zorig_ovrf, orig_scale, orig_scale_bits, natom,
                                            umat_dest, invu_dest, bdim_dest, umat_orig, invu_orig,
                                            bdim_orig, dest_atom_offset, orig_atom_offset,
                                            dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset,
                                            orig_bdim_offset);
  }
  else if (ct_dest == llint_type_index) {
    llint* lli_xdest = reinterpret_cast<llint*>(xdest);
    llint* lli_ydest = reinterpret_cast<llint*>(ydest);
    llint* lli_zdest = reinterpret_cast<llint*>(zdest);
    kCopyCoordinateXYZ<<<nblock, nthread>>>(lli_xdest, lli_ydest, lli_zdest, dest_scale,
                                            dest_scale_bits, xorig, yorig, zorig, xorig_ovrf,
                                            yorig_ovrf, zorig_ovrf, orig_scale, orig_scale_bits,
                                            natom, umat_dest, invu_dest, bdim_dest, umat_orig,
                                            invu_orig, bdim_orig, dest_atom_offset,
                                            orig_atom_offset, dest_xfrm_offset, orig_xfrm_offset,
                                            dest_bdim_offset, orig_bdim_offset);
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinateXYZ(void* xdest, void* ydest, void* zdest, const double dest_scale,
                                    const int dest_scale_bits, const size_t ct_dest,
                                    const llint* xorig, const llint* yorig, const llint* zorig,
                                    const int* xorig_ovrf, const int* yorig_ovrf,
                                    const int* zorig_ovrf, const double orig_scale,
                                    const int orig_scale_bits, const int natom,
                                    const int dest_atom_offset, const int orig_atom_offset,
                                    const GpuDetails &gpu) {
  launchCopyCoordinateXYZBox(xdest, ydest, zdest, dest_scale, dest_scale_bits, ct_dest, xorig,
                             yorig, zorig, xorig_ovrf, yorig_ovrf, zorig_ovrf, orig_scale,
                             orig_scale_bits, natom, nullptr, nullptr, nullptr, nullptr, nullptr,
                             nullptr, dest_atom_offset, orig_atom_offset, 0, 0, 0, 0, gpu);
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinateXYZBox(llint* xdest, llint* ydest, llint* zdest, int* xdest_ovrf,
                                       int* ydest_ovrf, int* zdest_ovrf, const double dest_scale,
                                       const int dest_scale_bits, const void* xorig,
                                       const void* yorig, const void* zorig,
                                       const double orig_scale, const int orig_scale_bits,
                                       const size_t ct_orig, const int natom, double* umat_dest,
                                       double* invu_dest, double* bdim_dest,
                                       const double* umat_orig, const double* invu_orig,
                                       const double* bdim_orig, llint* boxvecs_dest,
                                       int* boxvec_ovrf_dest, const int dest_atom_offset,
                                       const int orig_atom_offset, const int dest_xfrm_offset,
                                       const int orig_xfrm_offset, const int dest_bdim_offset,
                                       const int orig_bdim_offset, const GpuDetails &gpu) {
  const int nblock = 2 * gpu.getSMPCount();
  const int nthread = medium_block_size;
  if (ct_orig == double_type_index) {
    const double* d_xorig = reinterpret_cast<const double*>(xorig);
    const double* d_yorig = reinterpret_cast<const double*>(yorig);
    const double* d_zorig = reinterpret_cast<const double*>(zorig);
    kCopyCoordinateXYZ<<<nblock, nthread>>>(xdest, ydest, zdest, xdest_ovrf, ydest_ovrf,
                                            zdest_ovrf, dest_scale, dest_scale_bits, d_xorig,
                                            d_yorig, d_zorig, orig_scale, orig_scale_bits, natom,
                                            umat_dest, invu_dest, bdim_dest, umat_orig, invu_orig,
                                            bdim_orig, boxvecs_dest, boxvec_ovrf_dest,
                                            dest_atom_offset, orig_atom_offset, dest_xfrm_offset,
                                            orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset);
  }
  else if (ct_orig == float_type_index) {
    const float* f_xorig = reinterpret_cast<const float*>(xorig);
    const float* f_yorig = reinterpret_cast<const float*>(yorig);
    const float* f_zorig = reinterpret_cast<const float*>(zorig);
    kCopyCoordinateXYZ<<<nblock, nthread>>>(xdest, ydest, zdest, xdest_ovrf, ydest_ovrf,
                                            zdest_ovrf, dest_scale, dest_scale_bits, f_xorig,
                                            f_yorig, f_zorig, orig_scale, orig_scale_bits, natom,
                                            umat_dest, invu_dest, bdim_dest, umat_orig, invu_orig,
                                            bdim_orig, boxvecs_dest, boxvec_ovrf_dest,
                                            dest_atom_offset, orig_atom_offset, dest_xfrm_offset,
                                            orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset);
  }
  else if (ct_orig == short_type_index) {
    const short int* shi_xorig = reinterpret_cast<const short int*>(xorig);
    const short int* shi_yorig = reinterpret_cast<const short int*>(yorig);
    const short int* shi_zorig = reinterpret_cast<const short int*>(zorig);
    kCopyCoordinateXYZ<<<nblock, nthread>>>(xdest, ydest, zdest, xdest_ovrf, ydest_ovrf,
                                            zdest_ovrf, dest_scale, dest_scale_bits, shi_xorig,
                                            shi_yorig, shi_zorig, orig_scale, orig_scale_bits,
                                            natom, umat_dest, invu_dest, bdim_dest, umat_orig,
                                            invu_orig, bdim_orig, boxvecs_dest, boxvec_ovrf_dest,
                                            dest_atom_offset, orig_atom_offset, dest_xfrm_offset,
                                            orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset);
  }
  else if (ct_orig == int_type_index) {
    const int* i_xorig = reinterpret_cast<const int*>(xorig);
    const int* i_yorig = reinterpret_cast<const int*>(yorig);
    const int* i_zorig = reinterpret_cast<const int*>(zorig);
    kCopyCoordinateXYZ<<<nblock, nthread>>>(xdest, ydest, zdest, xdest_ovrf, ydest_ovrf,
                                            zdest_ovrf, dest_scale, dest_scale_bits, i_xorig,
                                            i_yorig, i_zorig, orig_scale, orig_scale_bits, natom,
                                            umat_dest, invu_dest, bdim_dest, umat_orig, invu_orig,
                                            bdim_orig, boxvecs_dest, boxvec_ovrf_dest,
                                            dest_atom_offset, orig_atom_offset, dest_xfrm_offset,
                                            orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset);
  }
  else if (ct_orig == llint_type_index) {
    const llint* lli_xorig = reinterpret_cast<const llint*>(xorig);
    const llint* lli_yorig = reinterpret_cast<const llint*>(yorig);
    const llint* lli_zorig = reinterpret_cast<const llint*>(zorig);
    kCopyCoordinateXYZ<<<nblock, nthread>>>(xdest, ydest, zdest, xdest_ovrf, ydest_ovrf,
                                            zdest_ovrf, dest_scale, dest_scale_bits, lli_xorig,
                                            lli_yorig, lli_zorig, orig_scale, orig_scale_bits,
                                            natom, umat_dest, invu_dest, bdim_dest, umat_orig,
                                            invu_orig, bdim_orig, boxvecs_dest, boxvec_ovrf_dest,
                                            dest_atom_offset, orig_atom_offset, dest_xfrm_offset,
                                            orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset);
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinateXYZ(llint* xdest, llint* ydest, llint* zdest, int* xdest_ovrf,
                                    int* ydest_ovrf, int* zdest_ovrf, const double dest_scale,
                                    const int dest_scale_bits, const void* xorig,
                                    const void* yorig, const void* zorig, const double orig_scale,
                                    const int orig_scale_bits, const size_t ct_orig,
                                    const int natom, const int dest_atom_offset,
                                    const int orig_atom_offset, const GpuDetails &gpu) {
  launchCopyCoordinateXYZBox(xdest, ydest, zdest, xdest_ovrf, ydest_ovrf, zdest_ovrf, dest_scale,
                             dest_scale_bits, xorig, yorig, zorig, orig_scale, orig_scale_bits,
                             ct_orig, natom, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                             nullptr, nullptr, dest_atom_offset, orig_atom_offset, 0, 0, 0, 0,
                             gpu);
}

//-------------------------------------------------------------------------------------------------
// Copy particle coordinates (positions, velocities, or forces) for Cartesian X, Y, and Z
// dimensions in one system.  The abstracts submitted to this function comprise all inputs needed
// by prior overloaded kernels.
//
// Arguments:
//   dest:  Object into which coordinates shall be written
//   orig:  Object from which coordinates shall be obtained
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(medium_block_size, 2)
kCopyCoordinates(PhaseSpaceWriter dest, const PhaseSpaceReader orig) {

  // Copy the box information
  if (blockIdx.x == 0) {
    kCopyBoxInformation(dest.umat, dest.invu, dest.boxdim, orig.umat, orig.invu, orig.boxdim,
                        0, 0, 0, 0, 0);
    kCopyBoxInformation(dest.umat_alt, dest.invu_alt, dest.boxdim_alt, orig.umat_alt,
                        orig.invu_alt, orig.boxdim_alt, 0, 0, 0, 0, 3);
  }
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t stride = devcRoundUp(dest.natom, warp_size_int);
  const size_t advance = (blockDim.x * gridDim.x);
  pos = copyCoordinateSet(dest.xcrd, 0, 1.0, 0, orig.xcrd, 0, 1.0, 0, orig.natom, pos,  0, stride,
                          advance);
  pos = copyCoordinateSet(dest.ycrd, 0, 1.0, 0, orig.ycrd, 0, 1.0, 0, orig.natom, pos,  1, stride,
                          advance);
  pos = copyCoordinateSet(dest.zcrd, 0, 1.0, 0, orig.zcrd, 0, 1.0, 0, orig.natom, pos,  2, stride,
                          advance);
  pos = copyCoordinateSet(dest.xvel, 0, 1.0, 0, orig.xvel, 0, 1.0, 0, orig.natom, pos,  3, stride,
                          advance);
  pos = copyCoordinateSet(dest.yvel, 0, 1.0, 0, orig.yvel, 0, 1.0, 0, orig.natom, pos,  4, stride,
                          advance);
  pos = copyCoordinateSet(dest.zvel, 0, 1.0, 0, orig.zvel, 0, 1.0, 0, orig.natom, pos,  5, stride,
                          advance);
  pos = copyCoordinateSet(dest.xfrc, 0, 1.0, 0, orig.xfrc, 0, 1.0, 0, orig.natom, pos,  6, stride,
                          advance);
  pos = copyCoordinateSet(dest.yfrc, 0, 1.0, 0, orig.yfrc, 0, 1.0, 0, orig.natom, pos,  7, stride,
                          advance);
  pos = copyCoordinateSet(dest.zfrc, 0, 1.0, 0, orig.zfrc, 0, 1.0, 0, orig.natom, pos,  8, stride,
                          advance);
  pos = copyCoordinateSet(dest.xalt, 0, 1.0, 0, orig.xalt, 0, 1.0, 0, orig.natom, pos,  9, stride,
                          advance);
  pos = copyCoordinateSet(dest.yalt, 0, 1.0, 0, orig.yalt, 0, 1.0, 0, orig.natom, pos, 10, stride,
                          advance);
  pos = copyCoordinateSet(dest.zalt, 0, 1.0, 0, orig.zalt, 0, 1.0, 0, orig.natom, pos, 11, stride,
                          advance);
  pos = copyCoordinateSet(dest.vxalt, 0, 1.0, 0, orig.vxalt, 0, 1.0, 0, orig.natom, pos, 12,
                          stride, advance);
  pos = copyCoordinateSet(dest.vyalt, 0, 1.0, 0, orig.vyalt, 0, 1.0, 0, orig.natom, pos, 13,
                          stride, advance);
  pos = copyCoordinateSet(dest.vzalt, 0, 1.0, 0, orig.vzalt, 0, 1.0, 0, orig.natom, pos, 14,
                          stride, advance);
  pos = copyCoordinateSet(dest.fxalt, 0, 1.0, 0, orig.fxalt, 0, 1.0, 0, orig.natom, pos, 15,
                          stride, advance);
  pos = copyCoordinateSet(dest.fyalt, 0, 1.0, 0, orig.fyalt, 0, 1.0, 0, orig.natom, pos, 16,
                          stride, advance);
  pos = copyCoordinateSet(dest.fzalt, 0, 1.0, 0, orig.fzalt, 0, 1.0, 0, orig.natom, pos, 17,
                          stride, advance);
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinates(PhaseSpaceWriter *destination, const PhaseSpaceReader &origin,
                                  const GpuDetails &gpu) {
  const int nblock = 2 * gpu.getSMPCount();
  kCopyCoordinates<<<nblock, medium_block_size>>>(*destination, origin);
}

//-------------------------------------------------------------------------------------------------
// Copy particle coordinates (positions, velocities, or forces) for Cartesian X, Y, and Z
// dimensions in one system.  The abstracts submitted to this function comprise all inputs needed
// by prior overloaded kernels.  Descriptions of all parameters follow from prior overloaded
// implementations.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(medium_block_size, 2)
kCopyCoordinates(PhaseSpaceWriter dest, const PsSynthesisReader orig,
                 const size_t orig_atom_offset, const size_t orig_xfrm_offset,
                 const size_t orig_bdim_offset) {

  // Copy the box information
  if (blockIdx.x == 0) {
    kCopyBoxInformation(dest.umat, dest.invu, dest.boxdim, orig.umat, orig.invu, orig.boxdims,
                        0, orig_xfrm_offset, 0, orig_bdim_offset, 0);
    kCopyBoxInformation(dest.umat_alt, dest.invu_alt, dest.boxdim_alt, orig.umat_alt,
                        orig.invu_alt, orig.alt_boxdims, 0, orig_xfrm_offset, 0, orig_bdim_offset,
                        3);
  }
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t stride = devcRoundUp(dest.natom, warp_size_int);
  const size_t advance = (blockDim.x * gridDim.x);
  pos = copyCoordinateSet(dest.xcrd, 0, 0, orig.xcrd, orig.xcrd_ovrf, orig_atom_offset,
                          orig.gpos_scale, orig.gpos_bits, dest.natom, pos,  0, stride, advance);
  pos = copyCoordinateSet(dest.ycrd, 0, 0, orig.ycrd, orig.ycrd_ovrf, orig_atom_offset,
                          orig.gpos_scale, orig.gpos_bits, dest.natom, pos,  1, stride, advance);
  pos = copyCoordinateSet(dest.zcrd, 0, 0, orig.zcrd, orig.zcrd_ovrf, orig_atom_offset,
                          orig.gpos_scale, orig.gpos_bits, dest.natom, pos,  2, stride, advance);
  pos = copyCoordinateSet(dest.xvel, 0, 0, orig.xvel, orig.xvel_ovrf, orig_atom_offset,
                          orig.vel_scale, orig.gpos_bits, dest.natom, pos,   3, stride, advance);
  pos = copyCoordinateSet(dest.yvel, 0, 0, orig.yvel, orig.yvel_ovrf, orig_atom_offset,
                          orig.vel_scale, orig.gpos_bits, dest.natom, pos,   4, stride, advance);
  pos = copyCoordinateSet(dest.zvel, 0, 0, orig.zvel, orig.zvel_ovrf, orig_atom_offset,
                          orig.vel_scale, orig.gpos_bits, dest.natom, pos,   5, stride, advance);
  pos = copyCoordinateSet(dest.xfrc, 0, 0, orig.xfrc, orig.xfrc_ovrf, orig_atom_offset,
                          orig.frc_scale, orig.gpos_bits, dest.natom, pos,   6, stride, advance);
  pos = copyCoordinateSet(dest.yfrc, 0, 0, orig.yfrc, orig.yfrc_ovrf, orig_atom_offset,
                          orig.frc_scale, orig.gpos_bits, dest.natom, pos,   7, stride, advance);
  pos = copyCoordinateSet(dest.zfrc, 0, 0, orig.zfrc, orig.zfrc_ovrf, orig_atom_offset,
                          orig.frc_scale, orig.gpos_bits, dest.natom, pos,   8, stride, advance);
  pos = copyCoordinateSet(dest.xalt, 0, 0, orig.xalt, orig.xalt_ovrf, orig_atom_offset,
                          orig.gpos_scale, orig.gpos_bits, dest.natom, pos,  9, stride, advance);
  pos = copyCoordinateSet(dest.yalt, 0, 0, orig.yalt, orig.yalt_ovrf, orig_atom_offset,
                          orig.gpos_scale, orig.gpos_bits, dest.natom, pos, 10, stride, advance);
  pos = copyCoordinateSet(dest.zalt, 0, 0, orig.zalt, orig.zalt_ovrf, orig_atom_offset,
                          orig.gpos_scale, orig.gpos_bits, dest.natom, pos, 11, stride, advance);
  pos = copyCoordinateSet(dest.vxalt, 0, 0, orig.vxalt, orig.vxalt_ovrf, orig_atom_offset,
                          orig.vel_scale, orig.gpos_bits, dest.natom, pos,  12, stride, advance);
  pos = copyCoordinateSet(dest.vyalt, 0, 0, orig.vyalt, orig.vyalt_ovrf, orig_atom_offset,
                          orig.vel_scale, orig.gpos_bits, dest.natom, pos,  13, stride, advance);
  pos = copyCoordinateSet(dest.vzalt, 0, 0, orig.vzalt, orig.vzalt_ovrf, orig_atom_offset,
                          orig.vel_scale, orig.gpos_bits, dest.natom, pos,  14, stride, advance);
  pos = copyCoordinateSet(dest.fxalt, 0, 0, orig.fxalt, orig.fxalt_ovrf, orig_atom_offset,
                          orig.frc_scale, orig.gpos_bits, dest.natom, pos,  15, stride, advance);
  pos = copyCoordinateSet(dest.fyalt, 0, 0, orig.fyalt, orig.fyalt_ovrf, orig_atom_offset,
                          orig.frc_scale, orig.gpos_bits, dest.natom, pos,  16, stride, advance);
  pos = copyCoordinateSet(dest.fzalt, 0, 0, orig.fzalt, orig.fzalt_ovrf, orig_atom_offset,
                          orig.frc_scale, orig.gpos_bits, dest.natom, pos,  17, stride, advance);
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinates(PhaseSpaceWriter *destination, const PsSynthesisReader &origin,
                                  const size_t orig_atom_offset, const int index_orig,
                                  const GpuDetails &gpu) {
  const int nblock = 2 * gpu.getSMPCount();

  // Compute the offsets here and pass them to the kernel.  Host memory will have lower latency
  // device memory, and computations of the offsets would simply have to be computed on every
  // thread, so shift the work to the CPU.  The API call for launchCopyCoordinates() remains as
  // simple as possible and compute time is optimized.
  const size_t orig_xfrm_offset = roundUp(9, warp_size_int) * index_orig;
  const size_t orig_bdim_offset = roundUp(6, warp_size_int) * index_orig;
  kCopyCoordinates<<<nblock, medium_block_size>>>(*destination, origin, orig_atom_offset,
                                                  orig_xfrm_offset, orig_bdim_offset);
}

//-------------------------------------------------------------------------------------------------
// Copy particle coordinates (positions, velocities, or forces) for Cartesian X, Y, and Z
// dimensions in one system.  The abstracts submitted to this function comprise all inputs needed
// by prior overloaded kernels.  Descriptions of all parameters follow from prior overloaded
// implementations.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(medium_block_size, 2)
kCopyCoordinates(PsSynthesisWriter dest, const size_t dest_atom_offset,
                 const size_t dest_xfrm_offset, const size_t dest_bdim_offset,
                 const PhaseSpaceReader orig) {

  // Copy the box information
  if (blockIdx.x == 0) {
    kCopyBoxInformation(dest.umat, dest.invu, dest.boxdims, orig.umat, orig.invu, orig.boxdim,
                        dest_xfrm_offset, 0, dest_bdim_offset, 0, 0, dest.boxvecs,
                        dest.boxvec_ovrf, dest.gpos_scale);
    kCopyBoxInformation(dest.umat_alt, dest.invu_alt, dest.alt_boxdims, orig.umat_alt,
                        orig.invu_alt, orig.boxdim_alt, dest_xfrm_offset, 0, dest_bdim_offset, 0,
                        3, dest.alt_boxvecs, dest.alt_boxvec_ovrf, dest.gpos_scale);
  }
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t stride = devcRoundUp(orig.natom, warp_size_int);
  const size_t advance = (blockDim.x * gridDim.x);
  pos = copyCoordinateSet(dest.xcrd, dest.xcrd_ovrf, dest_atom_offset, dest.gpos_scale,
                          dest.gpos_bits, orig.xcrd, 0, 0, orig.natom, pos,  0, stride, advance);
  pos = copyCoordinateSet(dest.ycrd, dest.ycrd_ovrf, dest_atom_offset, dest.gpos_scale,
                          dest.gpos_bits, orig.ycrd, 0, 0, orig.natom, pos,  1, stride, advance);
  pos = copyCoordinateSet(dest.zcrd, dest.zcrd_ovrf, dest_atom_offset, dest.gpos_scale,
                          dest.gpos_bits, orig.zcrd, 0, 0, orig.natom, pos,  2, stride, advance);
  pos = copyCoordinateSet(dest.xvel, dest.xvel_ovrf, dest_atom_offset, dest.vel_scale,
                          dest.vel_bits, orig.xvel, 0, 0, orig.natom, pos,   3, stride, advance);
  pos = copyCoordinateSet(dest.yvel, dest.yvel_ovrf, dest_atom_offset, dest.vel_scale,
                          dest.vel_bits, orig.yvel, 0, 0, orig.natom, pos,   4, stride, advance);
  pos = copyCoordinateSet(dest.zvel, dest.zvel_ovrf, dest_atom_offset, dest.vel_scale,
                          dest.vel_bits, orig.zvel, 0, 0, orig.natom, pos,   5, stride, advance);
  pos = copyCoordinateSet(dest.xfrc, dest.xfrc_ovrf, dest_atom_offset, dest.frc_scale,
                          dest.frc_bits, orig.xfrc, 0, 0, orig.natom, pos,   6, stride, advance);
  pos = copyCoordinateSet(dest.yfrc, dest.yfrc_ovrf, dest_atom_offset, dest.frc_scale,
                          dest.frc_bits, orig.yfrc, 0, 0, orig.natom, pos,   7, stride, advance);
  pos = copyCoordinateSet(dest.zfrc, dest.zfrc_ovrf, dest_atom_offset, dest.frc_scale,
                          dest.frc_bits, orig.zfrc, 0, 0, orig.natom, pos,   8, stride, advance);
  pos = copyCoordinateSet(dest.xalt, dest.xalt_ovrf, dest_atom_offset, dest.gpos_scale,
                          dest.gpos_bits, orig.xalt, 0, 0, orig.natom, pos,  9, stride, advance);
  pos = copyCoordinateSet(dest.yalt, dest.yalt_ovrf, dest_atom_offset, dest.gpos_scale,
                          dest.gpos_bits, orig.yalt, 0, 0, orig.natom, pos, 10, stride, advance);
  pos = copyCoordinateSet(dest.zalt, dest.zalt_ovrf, dest_atom_offset, dest.gpos_scale,
                          dest.gpos_bits, orig.zalt, 0, 0, orig.natom, pos, 11, stride, advance);
  pos = copyCoordinateSet(dest.vxalt, dest.vxalt_ovrf, dest_atom_offset, dest.vel_scale,
                          dest.vel_bits, orig.vxalt, 0, 0, orig.natom, pos,  12, stride, advance);
  pos = copyCoordinateSet(dest.vyalt, dest.vyalt_ovrf, dest_atom_offset, dest.vel_scale,
                          dest.vel_bits, orig.vyalt, 0, 0, orig.natom, pos,  13, stride, advance);
  pos = copyCoordinateSet(dest.vzalt, dest.vzalt_ovrf, dest_atom_offset, dest.vel_scale,
                          dest.vel_bits, orig.vzalt, 0, 0, orig.natom, pos,  14, stride, advance);
  pos = copyCoordinateSet(dest.fxalt, dest.fxalt_ovrf, dest_atom_offset, dest.frc_scale,
                          dest.frc_bits, orig.fxalt, 0, 0, orig.natom, pos,  15, stride, advance);
  pos = copyCoordinateSet(dest.fyalt, dest.fyalt_ovrf, dest_atom_offset, dest.frc_scale,
                          dest.frc_bits, orig.fyalt, 0, 0, orig.natom, pos,  16, stride, advance);
  pos = copyCoordinateSet(dest.fzalt, dest.fzalt_ovrf, dest_atom_offset, dest.frc_scale,
                          dest.frc_bits, orig.fzalt, 0, 0, orig.natom, pos,  17, stride, advance);
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinates(PsSynthesisWriter *destination, const size_t dest_atom_offset,
                                  const int index_dest, const PhaseSpaceReader &origin,
                                  const GpuDetails &gpu) {
  const int nblock = 2 * gpu.getSMPCount();

  // Compute the offsets here and pass them to the kernel.  Host memory will have lower latency
  // device memory, and computations of the offsets would simply have to be computed on every
  // thread, so shift the work to the CPU.  The API call for launchCopyCoordinates() remains as
  // simple as possible and compute time is optimized.
  const size_t dest_xfrm_offset = roundUp(9, warp_size_int) * index_dest;
  const size_t dest_bdim_offset = roundUp(6, warp_size_int) * index_dest;
  kCopyCoordinates<<<nblock, medium_block_size>>>(*destination, dest_atom_offset, dest_xfrm_offset,
                                                  dest_bdim_offset, origin);
}

//-------------------------------------------------------------------------------------------------
// Copy one coordinate set into another when both coordinate sets are based on 95-bit
// fixed-precision representations.  The arguments to this function follow from additional,
// templated overloads found in hpc_coordinate_copy.cuh. 
//-------------------------------------------------------------------------------------------------
__device__ __forceinline__
size_t copyCoordinateSet(llint* dest_crd, int* dest_crd_ovrf, const int dest_start_idx,
                         const int dest_bits, const llint* orig_crd, const int* orig_crd_ovrf,
                         const int orig_start_idx, const int orig_bits, const int count,
                         const size_t pos, const size_t iter, const size_t stride,
                         const size_t advance) {
  size_t ipos = pos;
  const size_t pos_limit = (iter + 1) * stride;
  while (ipos < pos_limit) {
    const int rel_pos = ipos - (iter * stride);
    if (rel_pos < count) {
      const size_t orp_idx = orig_start_idx + rel_pos;
      const int95_t orig_val = { __ldlu(&orig_crd[orp_idx]), __ldlu(&orig_crd_ovrf[orp_idx]) };
      const int95_t dest_val = changeFPBits(orig_val, orig_bits, dest_bits);
      const size_t drp_idx = dest_start_idx + rel_pos;
      __stwt(&dest_crd[drp_idx], dest_val.x);
      __stwt(&dest_crd_ovrf[drp_idx], dest_val.y);
    }
    ipos += advance;
  }
  return ipos;
}
  
//-------------------------------------------------------------------------------------------------
// Copy particle coordinates (positions, velocities, or forces) for Cartesian X, Y, and Z
// dimensions in one system.  The abstracts submitted to this function comprise all inputs needed
// by prior overloaded kernels.  Descriptions of all parameters follow from prior overloaded
// implementations.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(medium_block_size, 2)
kCopyCoordinates(PsSynthesisWriter dest, const size_t dest_atom_offset,
                 const size_t dest_xfrm_offset, const size_t dest_bdim_offset,
                 const PsSynthesisReader orig, const size_t orig_atom_offset,
                 const size_t orig_xfrm_offset, const size_t orig_bdim_offset, const int natom) {

  // Copy the box information
  if (blockIdx.x == 0) {
    kCopyBoxInformation(dest.umat, dest.invu, dest.boxdims, orig.umat, orig.invu, orig.boxdims,
                        dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset, 0,
                        dest.boxvecs, dest.boxvec_ovrf, dest.gpos_scale, orig.boxvecs,
                        orig.boxvec_ovrf, dest.gpos_bits, orig.gpos_bits);
    kCopyBoxInformation(dest.umat_alt, dest.invu_alt, dest.alt_boxdims, orig.umat_alt,
                        orig.invu_alt, orig.alt_boxdims, dest_xfrm_offset, orig_xfrm_offset,
                        dest_bdim_offset, orig_bdim_offset, 3, dest.alt_boxvecs,
                        dest.alt_boxvec_ovrf, dest.gpos_scale, orig.alt_boxvecs,
                        orig.alt_boxvec_ovrf, dest.gpos_bits, orig.gpos_bits);
  }
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t stride = devcRoundUp(natom, warp_size_int);
  const size_t advance = (blockDim.x * gridDim.x);
  pos = copyCoordinateSet(dest.xcrd, dest.xcrd_ovrf, dest_atom_offset, dest.gpos_bits, orig.xcrd,
                          orig.xcrd_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos,  0,
                          stride, advance);
  pos = copyCoordinateSet(dest.ycrd, dest.ycrd_ovrf, dest_atom_offset, dest.gpos_bits, orig.ycrd,
                          orig.ycrd_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos,  1,
                          stride, advance);
  pos = copyCoordinateSet(dest.zcrd, dest.zcrd_ovrf, dest_atom_offset, dest.gpos_bits, orig.zcrd,
                          orig.zcrd_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos,  2,
                          stride, advance);
  pos = copyCoordinateSet(dest.xvel, dest.xvel_ovrf, dest_atom_offset, dest.vel_bits, orig.xvel,
                          orig.xvel_ovrf, orig_atom_offset, orig.vel_bits, natom, pos,   3,
                          stride, advance);
  pos = copyCoordinateSet(dest.yvel, dest.yvel_ovrf, dest_atom_offset, dest.vel_bits, orig.yvel,
                          orig.yvel_ovrf, orig_atom_offset, orig.vel_bits, natom, pos,   4,
                          stride, advance);
  pos = copyCoordinateSet(dest.zvel, dest.zvel_ovrf, dest_atom_offset, dest.vel_bits, orig.zvel,
                          orig.zvel_ovrf, orig_atom_offset, orig.vel_bits, natom, pos,   5,
                          stride, advance);
  pos = copyCoordinateSet(dest.xfrc, dest.xfrc_ovrf, dest_atom_offset, dest.frc_bits, orig.xfrc,
                          orig.xfrc_ovrf, orig_atom_offset, orig.frc_bits, natom, pos,   6,
                          stride, advance);
  pos = copyCoordinateSet(dest.yfrc, dest.yfrc_ovrf, dest_atom_offset, dest.frc_bits, orig.yfrc,
                          orig.yfrc_ovrf, orig_atom_offset, orig.frc_bits, natom, pos,   7,
                          stride, advance);
  pos = copyCoordinateSet(dest.zfrc, dest.zfrc_ovrf, dest_atom_offset, dest.frc_bits, orig.zfrc,
                          orig.zfrc_ovrf, orig_atom_offset, orig.frc_bits, natom, pos,   8,
                          stride, advance);
  pos = copyCoordinateSet(dest.xalt, dest.xalt_ovrf, dest_atom_offset, dest.gpos_bits, orig.xalt,
                          orig.xalt_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos,  9,
                          stride, advance);
  pos = copyCoordinateSet(dest.yalt, dest.yalt_ovrf, dest_atom_offset, dest.gpos_bits, orig.yalt,
                          orig.yalt_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos, 10,
                          stride, advance);
  pos = copyCoordinateSet(dest.zalt, dest.zalt_ovrf, dest_atom_offset, dest.gpos_bits, orig.zalt,
                          orig.zalt_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos, 11,
                          stride, advance);
  pos = copyCoordinateSet(dest.vxalt, dest.vxalt_ovrf, dest_atom_offset, dest.vel_bits, orig.vxalt,
                          orig.vxalt_ovrf, orig_atom_offset, orig.vel_bits, natom, pos,  12,
                          stride, advance);
  pos = copyCoordinateSet(dest.vyalt, dest.vyalt_ovrf, dest_atom_offset, dest.vel_bits, orig.vyalt,
                          orig.vyalt_ovrf, orig_atom_offset, orig.vel_bits, natom, pos,  13,
                          stride, advance);
  pos = copyCoordinateSet(dest.vzalt, dest.vzalt_ovrf, dest_atom_offset, dest.vel_bits, orig.vzalt,
                          orig.vzalt_ovrf, orig_atom_offset, orig.vel_bits, natom, pos,  14,
                          stride, advance);
  pos = copyCoordinateSet(dest.fxalt, dest.fxalt_ovrf, dest_atom_offset, dest.frc_bits, orig.fxalt,
                          orig.fxalt_ovrf, orig_atom_offset, orig.frc_bits, natom, pos,  15,
                          stride, advance);
  pos = copyCoordinateSet(dest.fyalt, dest.fyalt_ovrf, dest_atom_offset, dest.frc_bits, orig.fyalt,
                          orig.fyalt_ovrf, orig_atom_offset, orig.frc_bits, natom, pos,  16,
                          stride, advance);
  pos = copyCoordinateSet(dest.fzalt, dest.fzalt_ovrf, dest_atom_offset, dest.frc_bits, orig.fzalt,
                          orig.fzalt_ovrf, orig_atom_offset, orig.frc_bits, natom, pos,  17,
                          stride, advance);
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinates(PsSynthesisWriter *destination, const size_t dest_atom_offset,
                                  const int index_dest, const PsSynthesisReader &origin,
                                  const size_t orig_atom_offset, const int index_orig,
                                  const int natom, const GpuDetails &gpu) {
  const int nblock = 2 * gpu.getSMPCount();

  // Compute the offsets here and pass them to the kernel.  Host memory will have lower latency
  // device memory, and computations of the offsets would simply have to be computed on every
  // thread, so shift the work to the CPU.  The API call for launchCopyCoordinates() remains as
  // simple as possible and compute time is optimized.
  const int xfrm_w = roundUp(9, warp_size_int);
  const int bdim_w = roundUp(6, warp_size_int);
  const size_t dest_xfrm_offset = xfrm_w * index_dest;
  const size_t dest_bdim_offset = bdim_w * index_dest;
  const size_t orig_xfrm_offset = xfrm_w * index_orig;
  const size_t orig_bdim_offset = bdim_w * index_orig;
  kCopyCoordinates<<<nblock, medium_block_size>>>(*destination, dest_atom_offset, dest_xfrm_offset,
                                                  dest_bdim_offset, origin, orig_atom_offset,
                                                  orig_xfrm_offset, orig_bdim_offset, natom);
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinates(CoordinateSeriesWriter<void> *destination, size_t ct_dest,
                                  const CoordinateSeriesReader<void> &origin, size_t ct_orig,
                                  const int2* system_pairs, const int copy_count,
                                  const GpuDetails &gpu) {
  if (ct_dest == double_type_index) {
    CoordinateSeriesWriter<double> tr_dest = restoreType<double>(*destination);
    unrollCCXYZOrigin(&tr_dest, origin, ct_orig, system_pairs, copy_count, gpu);
  }
  else if (ct_dest == float_type_index) {
    CoordinateSeriesWriter<float> tr_dest = restoreType<float>(*destination);
    unrollCCXYZOrigin(&tr_dest, origin, ct_orig, system_pairs, copy_count, gpu);
  }
  else if (ct_dest == short_type_index) {
    CoordinateSeriesWriter<short int> tr_dest = restoreType<short int>(*destination);
    unrollCCXYZOrigin(&tr_dest, origin, ct_orig, system_pairs, copy_count, gpu);
  }
  else if (ct_dest == int_type_index) {
    CoordinateSeriesWriter<int> tr_dest = restoreType<int>(*destination);
    unrollCCXYZOrigin(&tr_dest, origin, ct_orig, system_pairs, copy_count, gpu);
  }
  else if (ct_dest == llint_type_index) {
    CoordinateSeriesWriter<llint> tr_dest = restoreType<llint>(*destination);
    unrollCCXYZOrigin(&tr_dest, origin, ct_orig, system_pairs, copy_count, gpu);
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinates(CoordinateSeriesWriter<void> *destination, size_t ct_dest,
                                  const PsSynthesisReader &origin, const TrajectoryKind kind,
                                  const int2* system_pairs, const int copy_count,
                                  const GpuDetails &gpu) {
  const int nblock = gpu.getSMPCount();
  if (ct_dest == double_type_index) {
    CoordinateSeriesWriter<double> tr_dest = restoreType<double>(destination);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(tr_dest, origin, kind, system_pairs,
                                                             copy_count);
  }
  else if (ct_dest == float_type_index) {
    CoordinateSeriesWriter<float> tr_dest = restoreType<float>(destination);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(tr_dest, origin, kind, system_pairs,
                                                             copy_count);
  }
  else if (ct_dest == short_type_index) {
    CoordinateSeriesWriter<short int> tr_dest = restoreType<short int>(destination);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(tr_dest, origin, kind, system_pairs,
                                                             copy_count);
  }
  else if (ct_dest == int_type_index) {
    CoordinateSeriesWriter<int> tr_dest = restoreType<int>(destination);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(tr_dest, origin, kind, system_pairs,
                                                             copy_count);
  }
  else if (ct_dest == llint_type_index) {
    CoordinateSeriesWriter<llint> tr_dest = restoreType<llint>(destination);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(tr_dest, origin, kind, system_pairs,
                                                             copy_count);
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinates(CoordinateSeriesWriter<void> *destination, size_t ct_dest,
                                  const CondensateReader &origin, const int2* system_pairs,
                                  const int copy_count, const GpuDetails &gpu) {
  const int nblock = gpu.getSMPCount();
  if (ct_dest == double_type_index) {
    CoordinateSeriesWriter<double> tr_dest = restoreType<double>(destination);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(tr_dest, origin, system_pairs,
                                                             copy_count);
  }
  else if (ct_dest == float_type_index) {
    CoordinateSeriesWriter<float> tr_dest = restoreType<float>(destination);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(tr_dest, origin, system_pairs,
                                                             copy_count);
  }
  else if (ct_dest == short_type_index) {
    CoordinateSeriesWriter<short int> tr_dest = restoreType<short int>(destination);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(tr_dest, origin, system_pairs,
                                                             copy_count);
  }
  else if (ct_dest == int_type_index) {
    CoordinateSeriesWriter<int> tr_dest = restoreType<int>(destination);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(tr_dest, origin, system_pairs,
                                                             copy_count);
  }
  else if (ct_dest == llint_type_index) {
    CoordinateSeriesWriter<llint> tr_dest = restoreType<llint>(destination);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(tr_dest, origin, system_pairs,
                                                             copy_count);
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinates(PsSynthesisWriter *destination, const TrajectoryKind kind,
                                  const CoordinateSeriesReader<void> &origin, size_t ct_orig,
                                  const int2* system_pairs, const int copy_count,
                                  const GpuDetails &gpu) {
  const int nblock = gpu.getSMPCount();
  if (ct_orig == double_type_index) {
    const CoordinateSeriesReader<double> tr_orig = restoreType<double>(origin);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(*destination, kind, tr_orig,
                                                             system_pairs, copy_count);
  }
  else if (ct_orig == float_type_index) {
    const CoordinateSeriesReader<float> tr_orig = restoreType<float>(origin);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(*destination, kind, tr_orig,
                                                             system_pairs, copy_count);
  }
  else if (ct_orig == short_type_index) {
    const CoordinateSeriesReader<short int> tr_orig = restoreType<short int>(origin);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(*destination, kind, tr_orig,
                                                             system_pairs, copy_count);
  }
  else if (ct_orig == int_type_index) {
    const CoordinateSeriesReader<int> tr_orig = restoreType<int>(origin);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(*destination, kind, tr_orig,
                                                             system_pairs, copy_count);
  }
  else if (ct_orig == llint_type_index) {
    const CoordinateSeriesReader<llint> tr_orig = restoreType<llint>(origin);
    kMultiSystemCoordinateCopy<<<nblock, large_block_size>>>(*destination, kind, tr_orig,
                                                             system_pairs, copy_count);
  }
}

//-------------------------------------------------------------------------------------------------
// Copy particle coordinates (positions, velocities, or forces) for Cartesian X, Y, and Z
// dimensions in one or more systems.  The abstracts submitted to this function comprise all
// inputs needed by prior overloaded kernels.  Descriptions of parameters follow from prior
// overloaded implementations.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kMultiSystemCoordinateCopy(PsSynthesisWriter dest, const PsSynthesisReader orig,
                           const int2* system_pairs, const int copy_count) {
  
  // Each block will take a separate system and carry out the copying, with any necessary
  // fixed-precision bit conversions.
  int pair_idx = blockIdx.x;
  while (pair_idx < copy_count) {
    const int2 origx_to_desty = __ldca(&system_pairs[pair_idx]);

    // Skip copy requests to or from invalid system indices.  This error will be caught by the
    // CPU while the kernel is running and raise an error.
    if (origx_to_desty.x < 0 || origx_to_desty.x >= orig.system_count ||
        origx_to_desty.y < 0 || origx_to_desty.y >= dest.system_count) {
      pair_idx += gridDim.x;
      continue;
    }
    const int natom = __ldca(&dest.atom_counts[origx_to_desty.y]);

    // Skip copy requests between systems of different sizes.  These will be caught by the CPU
    // while the kernel is running and raise an error.
    if (natom != __ldca(&orig.atom_counts[origx_to_desty.x])) {
      pair_idx += gridDim.x;
      continue;
    }

    // Copy the box information
    const int orig_atom_offset = __ldca(&orig.atom_starts[origx_to_desty.x]);
    const int dest_atom_offset = __ldca(&dest.atom_starts[origx_to_desty.y]);
    const int xfrm_w = devcRoundUp(9, warp_size_int);
    const int bdim_w = devcRoundUp(6, warp_size_int);
    const int orig_xfrm_offset = origx_to_desty.x * xfrm_w;
    const int dest_xfrm_offset = origx_to_desty.y * xfrm_w;
    const int orig_bdim_offset = origx_to_desty.x * bdim_w;
    const int dest_bdim_offset = origx_to_desty.y * bdim_w;
    kCopyBoxInformation(dest.umat, dest.invu, dest.boxdims, orig.umat, orig.invu, orig.boxdims,
                        dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset, 0,
                        dest.boxvecs, dest.boxvec_ovrf, dest.gpos_scale);
    kCopyBoxInformation(dest.umat_alt, dest.invu_alt, dest.alt_boxdims, orig.umat_alt,
                        orig.invu_alt, orig.alt_boxdims, dest_xfrm_offset, orig_xfrm_offset,
                        dest_bdim_offset, orig_bdim_offset, 3, dest.alt_boxvecs,
                        dest.alt_boxvec_ovrf, dest.gpos_scale);

    // Copy all types of coordinates
    size_t pos = threadIdx.x;
    const size_t stride = devcRoundUp(natom, warp_size_int);
    pos = copyCoordinateSet(dest.xcrd, dest.xcrd_ovrf, dest_atom_offset, dest.gpos_bits, orig.xcrd,
                            orig.xcrd_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos,  0,
                            stride, blockDim.x);
    pos = copyCoordinateSet(dest.ycrd, dest.ycrd_ovrf, dest_atom_offset, dest.gpos_bits, orig.ycrd,
                            orig.ycrd_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos,  1,
                            stride, blockDim.x);
    pos = copyCoordinateSet(dest.zcrd, dest.zcrd_ovrf, dest_atom_offset, dest.gpos_bits, orig.zcrd,
                            orig.zcrd_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos,  1,
                            stride, blockDim.x);
    pos = copyCoordinateSet(dest.xvel, dest.xvel_ovrf, dest_atom_offset, dest.vel_bits, orig.xvel,
                            orig.xvel_ovrf, orig_atom_offset, orig.vel_bits, natom, pos,   3,
                            stride, blockDim.x);
    pos = copyCoordinateSet(dest.yvel, dest.yvel_ovrf, dest_atom_offset, dest.vel_bits, orig.yvel,
                            orig.yvel_ovrf, orig_atom_offset, orig.vel_bits, natom, pos,   4,
                            stride, blockDim.x);
    pos = copyCoordinateSet(dest.zvel, dest.zvel_ovrf, dest_atom_offset, dest.vel_bits, orig.zvel,
                            orig.zvel_ovrf, orig_atom_offset, orig.vel_bits, natom, pos,   5,
                            stride, blockDim.x);
    pos = copyCoordinateSet(dest.xfrc, dest.xfrc_ovrf, dest_atom_offset, dest.frc_bits, orig.xfrc,
                            orig.xfrc_ovrf, orig_atom_offset, orig.frc_bits, natom, pos,   6,
                            stride, blockDim.x);
    pos = copyCoordinateSet(dest.yfrc, dest.yfrc_ovrf, dest_atom_offset, dest.frc_bits, orig.yfrc,
                            orig.yfrc_ovrf, orig_atom_offset, orig.frc_bits, natom, pos,   7,
                            stride, blockDim.x);
    pos = copyCoordinateSet(dest.zfrc, dest.zfrc_ovrf, dest_atom_offset, dest.frc_bits, orig.zfrc,
                            orig.zfrc_ovrf, orig_atom_offset, orig.frc_bits, natom, pos,   8,
                            stride, blockDim.x);
    pos = copyCoordinateSet(dest.xalt, dest.xalt_ovrf, dest_atom_offset, dest.gpos_bits, orig.xalt,
                            orig.xalt_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos,  9,
                            stride, blockDim.x);
    pos = copyCoordinateSet(dest.yalt, dest.yalt_ovrf, dest_atom_offset, dest.gpos_bits, orig.yalt,
                            orig.yalt_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos, 10,
                            stride, blockDim.x);
    pos = copyCoordinateSet(dest.zalt, dest.zalt_ovrf, dest_atom_offset, dest.gpos_bits, orig.zalt,
                            orig.zalt_ovrf, orig_atom_offset, orig.gpos_bits, natom, pos, 11,
                            stride, blockDim.x);
    pos = copyCoordinateSet(dest.vxalt, dest.vxalt_ovrf, dest_atom_offset, dest.vel_bits,
                            orig.vxalt, orig.vxalt_ovrf, orig_atom_offset, orig.vel_bits, natom,
                            pos,  12, stride, blockDim.x);
    pos = copyCoordinateSet(dest.vyalt, dest.vyalt_ovrf, dest_atom_offset, dest.vel_bits,
                            orig.vyalt, orig.vyalt_ovrf, orig_atom_offset, orig.vel_bits, natom,
                            pos,  13, stride, blockDim.x);
    pos = copyCoordinateSet(dest.vzalt, dest.vzalt_ovrf, dest_atom_offset, dest.vel_bits,
                            orig.vzalt, orig.vzalt_ovrf, orig_atom_offset, orig.vel_bits, natom,
                            pos,  14, stride, blockDim.x);
    pos = copyCoordinateSet(dest.fxalt, dest.fxalt_ovrf, dest_atom_offset, dest.frc_bits,
                            orig.fxalt, orig.fxalt_ovrf, orig_atom_offset, orig.frc_bits, natom,
                            pos,  15, stride, blockDim.x);
    pos = copyCoordinateSet(dest.fyalt, dest.fyalt_ovrf, dest_atom_offset, dest.frc_bits,
                            orig.fyalt, orig.fyalt_ovrf, orig_atom_offset, orig.frc_bits, natom,
                            pos,  16, stride, blockDim.x);
    pos = copyCoordinateSet(dest.fzalt, dest.fzalt_ovrf, dest_atom_offset, dest.frc_bits,
                            orig.fzalt, orig.fzalt_ovrf, orig_atom_offset, orig.frc_bits, natom,
                            pos,  17, stride, blockDim.x);
    pair_idx += gridDim.x;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinates(PsSynthesisWriter *destination, const PsSynthesisReader &origin,
                                  const int2* system_pairs, const int copy_count,
                                  const GpuDetails &gpu) {
  kMultiSystemCoordinateCopy<<<gpu.getSMPCount(), large_block_size>>>(*destination, origin,
                                                                      system_pairs, copy_count);
}

//-------------------------------------------------------------------------------------------------
// Copy particle coordinates (positions, velocities, or forces) for Cartesian X, Y, and Z
// dimensions in one or more systems.  The abstracts submitted to this function comprise all
// inputs needed by prior overloaded kernels.  Descriptions of parameters follow from prior
// overloaded implementations.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kMultiSystemCoordinateCopy(PsSynthesisWriter dest, const TrajectoryKind kind,
                           const CondensateReader orig, const int2* system_pairs,
                           const int copy_count) {
  int pair_idx = blockIdx.x;
  while (pair_idx < copy_count) {
    const int2 origx_to_desty = __ldca(&system_pairs[pair_idx]);
    if (origx_to_desty.x < 0 || origx_to_desty.x >= orig.system_count ||
        origx_to_desty.y < 0 || origx_to_desty.y >= dest.system_count) {
      pair_idx += gridDim.x;
      continue;
    }
    const int natom = __ldca(&dest.atom_counts[origx_to_desty.y]);
    if (natom != __ldca(&orig.atom_counts[origx_to_desty.x])) {
      pair_idx += gridDim.x;
      continue;
    }

    // Nested switches cover the destination's coordinate kinds and the origin's data type.
    const int orig_atom_offset = __ldca(&orig.atom_starts[origx_to_desty.x]);
    const int dest_atom_offset = __ldca(&dest.atom_starts[origx_to_desty.y]);
    size_t pos = threadIdx.x;
    const size_t padded_natom = devcRoundUp(natom, warp_size_int);
    switch (kind) {
    case TrajectoryKind::POSITIONS:
      {
        // Copy the box information
        const int xfrm_w = devcRoundUp(9, warp_size_int);
        const int bdim_w = devcRoundUp(6, warp_size_int);
        const int orig_xfrm_offset = origx_to_desty.x * xfrm_w;
        const int dest_xfrm_offset = origx_to_desty.y * xfrm_w;
        const int orig_bdim_offset = origx_to_desty.x * bdim_w;
        const int dest_bdim_offset = origx_to_desty.y * bdim_w;
        kCopyBoxInformation(dest.umat, dest.invu, dest.boxdims, orig.umat, orig.invu, orig.boxdims,
                            dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset,
                            0, dest.boxvecs, dest.boxvec_ovrf, dest.gpos_scale);

        // Copy atomic positions
        switch (orig.mode) {
        case PrecisionModel::DOUBLE:
          pos = copyCoordinateSet(dest.xcrd, dest.xcrd_ovrf, dest_atom_offset, dest.gpos_scale,
                                  dest.gpos_bits, orig.xcrd, orig_atom_offset, 0, natom, pos, 0,
                                  padded_natom, blockDim.x);
          pos = copyCoordinateSet(dest.ycrd, dest.ycrd_ovrf, dest_atom_offset, dest.gpos_scale,
                                  dest.gpos_bits, orig.ycrd, orig_atom_offset, 0, natom, pos, 1,
                                  padded_natom, blockDim.x);
          pos = copyCoordinateSet(dest.zcrd, dest.zcrd_ovrf, dest_atom_offset, dest.gpos_scale,
                                  dest.gpos_bits, orig.zcrd, orig_atom_offset, 0, natom, pos, 2,
                                  padded_natom, blockDim.x);
          break;
        case PrecisionModel::SINGLE:
          pos = copyCoordinateSet(dest.xcrd, dest.xcrd_ovrf, dest_atom_offset, dest.gpos_scale,
                                  dest.gpos_bits, orig.xcrd_sp, orig_atom_offset, 0, natom, pos, 0,
                                  padded_natom, blockDim.x);
          pos = copyCoordinateSet(dest.ycrd, dest.ycrd_ovrf, dest_atom_offset, dest.gpos_scale,
                                  dest.gpos_bits, orig.ycrd_sp, orig_atom_offset, 0, natom, pos, 1,
                                  padded_natom, blockDim.x);
          pos = copyCoordinateSet(dest.zcrd, dest.zcrd_ovrf, dest_atom_offset, dest.gpos_scale,
                                  dest.gpos_bits, orig.zcrd_sp, orig_atom_offset, 0, natom, pos, 2,
                                  padded_natom, blockDim.x);
          break;
        }
      }
      break;
    case TrajectoryKind::VELOCITIES:
      switch (orig.mode) {
      case PrecisionModel::DOUBLE:
        pos = copyCoordinateSet(dest.xvel, dest.xcrd_ovrf, dest_atom_offset, dest.vel_scale,
                                dest.vel_bits, orig.xcrd, orig_atom_offset, 0, natom, pos, 0,
                                padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.yvel, dest.ycrd_ovrf, dest_atom_offset, dest.vel_scale,
                                dest.vel_bits, orig.ycrd, orig_atom_offset, 0, natom, pos, 1,
                                padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.zvel, dest.zcrd_ovrf, dest_atom_offset, dest.vel_scale,
                                dest.vel_bits, orig.zcrd, orig_atom_offset, 0, natom, pos, 2,
                                padded_natom, blockDim.x);
        break;
      case PrecisionModel::SINGLE:
        pos = copyCoordinateSet(dest.xvel, dest.xcrd_ovrf, dest_atom_offset, dest.vel_scale,
                                dest.vel_bits, orig.xcrd_sp, orig_atom_offset, 0, natom, pos, 0,
                                padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.yvel, dest.ycrd_ovrf, dest_atom_offset, dest.vel_scale,
                                dest.vel_bits, orig.ycrd_sp, orig_atom_offset, 0, natom, pos, 1,
                                padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.zvel, dest.zcrd_ovrf, dest_atom_offset, dest.vel_scale,
                                dest.vel_bits, orig.zcrd_sp, orig_atom_offset, 0, natom, pos, 2,
                                padded_natom, blockDim.x);
        break;
      }
      break;
    case TrajectoryKind::FORCES:
      switch (orig.mode) {
      case PrecisionModel::DOUBLE:
        pos = copyCoordinateSet(dest.xfrc, dest.xcrd_ovrf, dest_atom_offset, dest.frc_scale,
                                dest.frc_bits, orig.xcrd, orig_atom_offset, 0, natom, pos, 0,
                                padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.yfrc, dest.ycrd_ovrf, dest_atom_offset, dest.frc_scale,
                                dest.frc_bits, orig.ycrd, orig_atom_offset, 0, natom, pos, 1,
                                padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.zfrc, dest.zcrd_ovrf, dest_atom_offset, dest.frc_scale,
                                dest.frc_bits, orig.zcrd, orig_atom_offset, 0, natom, pos, 2,
                                padded_natom, blockDim.x);
        break;
      case PrecisionModel::SINGLE:
        pos = copyCoordinateSet(dest.xfrc, dest.xcrd_ovrf, dest_atom_offset, dest.frc_scale,
                                dest.frc_bits, orig.xcrd_sp, orig_atom_offset, 0, natom, pos, 0,
                                padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.yfrc, dest.ycrd_ovrf, dest_atom_offset, dest.frc_scale,
                                dest.frc_bits, orig.ycrd_sp, orig_atom_offset, 0, natom, pos, 1,
                                padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.zfrc, dest.zcrd_ovrf, dest_atom_offset, dest.frc_scale,
                                dest.frc_bits, orig.zcrd_sp, orig_atom_offset, 0, natom, pos, 2,
                                padded_natom, blockDim.x);
        break;
      }
      break;
    }

    // Increment the system copy counter
    pair_idx += gridDim.x;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinates(PsSynthesisWriter *destination, const TrajectoryKind kind,
                                  const CondensateReader &origin, const int2* system_pairs,
                                  const int copy_count, const GpuDetails &gpu) {
  kMultiSystemCoordinateCopy<<<gpu.getSMPCount(), large_block_size>>>(*destination, kind, origin,
                                                                      system_pairs, copy_count);
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinates(CondensateWriter *destination,
                                  const CoordinateSeriesReader<void> &origin, const size_t ct_orig,
                                  const int2* system_pairs, const int copy_count,
                                  const GpuDetails &gpu) {
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

//-------------------------------------------------------------------------------------------------
// Copy multiple systems from a PhaseSpaceSynthesis into a Condensate object.  As with other
// templated coordinate copy kernels, the type reassignment will be done on the CPU, by the
// launching function.  Descriptions of input parameters follow from a previous overloaded variant
// of kMultiSystemCoordinateCopy() above, with changes in origin and destination object types.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kMultiSystemCoordinateCopy(CondensateWriter dest, const PsSynthesisReader orig,
                           const TrajectoryKind kind, const int2* system_pairs,
                           const int copy_count) {
  int pair_idx = blockIdx.x;
  while (pair_idx < copy_count) {
    const int2 origx_to_desty = __ldca(&system_pairs[pair_idx]);
    if (origx_to_desty.x < 0 || origx_to_desty.x >= orig.system_count ||
        origx_to_desty.y < 0 || origx_to_desty.y >= dest.system_count) {
      pair_idx += gridDim.x;
      continue;
    }
    const int natom = __ldca(&dest.atom_counts[origx_to_desty.y]);
    if (natom != __ldca(&orig.atom_counts[origx_to_desty.x])) {
      pair_idx += gridDim.x;
      continue;
    }

    // The box information, if it is to be copied, is independent of the Condensate's data type
    // (its mode).  Evaluate the switch over origin coordinate kinds as the outer nested switch.
    const int dest_atom_offset = __ldca(&dest.atom_starts[origx_to_desty.y]);
    const size_t orig_atom_offset = __ldca(&orig.atom_starts[origx_to_desty.x]);
    size_t pos = threadIdx.x;
    const size_t padded_natom = devcRoundUp(natom, warp_size_int);
    switch (kind) {
    case TrajectoryKind::POSITIONS:
      {
        // Copy the box information
        const int xfrm_w = devcRoundUp(9, warp_size_int);
        const int bdim_w = devcRoundUp(6, warp_size_int);
        const int orig_xfrm_offset = origx_to_desty.x * xfrm_w;
        const int dest_xfrm_offset = origx_to_desty.y * xfrm_w;
        const int orig_bdim_offset = origx_to_desty.x * bdim_w;
        const int dest_bdim_offset = origx_to_desty.y * bdim_w;
        kCopyBoxInformation(dest.umat, dest.invu, dest.boxdims, orig.umat, orig.invu, orig.boxdims,
                            dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset,
                            orig_bdim_offset);

        // Copy the atomic positions
        switch (dest.mode) {
        case PrecisionModel::DOUBLE:
          pos = copyCoordinateSet(dest.xcrd, dest_atom_offset, 0, orig.xcrd, orig.xcrd_ovrf,
                                  orig_atom_offset, orig.gpos_scale, orig.gpos_bits, natom, pos,
                                  0, padded_natom, blockDim.x);
          pos = copyCoordinateSet(dest.ycrd, dest_atom_offset, 0, orig.ycrd, orig.ycrd_ovrf,
                                  orig_atom_offset, orig.gpos_scale, orig.gpos_bits, natom, pos,
                                  1, padded_natom, blockDim.x);
          pos = copyCoordinateSet(dest.zcrd, dest_atom_offset, 0, orig.zcrd, orig.zcrd_ovrf,
                                  orig_atom_offset, orig.gpos_scale, orig.gpos_bits, natom, pos,
                                  2, padded_natom, blockDim.x);
          break;
        case PrecisionModel::SINGLE:
          pos = copyCoordinateSet(dest.xcrd_sp, dest_atom_offset, 0, orig.xcrd, orig.xcrd_ovrf,
                                  orig_atom_offset, orig.gpos_scale, orig.gpos_bits, natom, pos,
                                  0, padded_natom, blockDim.x);
          pos = copyCoordinateSet(dest.ycrd_sp, dest_atom_offset, 0, orig.ycrd, orig.ycrd_ovrf,
                                  orig_atom_offset, orig.gpos_scale, orig.gpos_bits, natom, pos,
                                  1, padded_natom, blockDim.x);
          pos = copyCoordinateSet(dest.zcrd_sp, dest_atom_offset, 0, orig.zcrd, orig.zcrd_ovrf,
                                  orig_atom_offset, orig.gpos_scale, orig.gpos_bits, natom, pos,
                                  2, padded_natom, blockDim.x);
          break;
        }
      }
      break;
    case TrajectoryKind::VELOCITIES:
      switch (dest.mode) {
      case PrecisionModel::DOUBLE:
        pos = copyCoordinateSet(dest.xcrd, dest_atom_offset, 0, orig.xvel, orig.xvel_ovrf,
                                orig_atom_offset, orig.vel_scale, orig.vel_bits, natom, pos,
                                0, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.ycrd, dest_atom_offset, 0, orig.yvel, orig.yvel_ovrf,
                                orig_atom_offset, orig.vel_scale, orig.vel_bits, natom, pos,
                                1, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.zcrd, dest_atom_offset, 0, orig.zvel, orig.zvel_ovrf,
                                orig_atom_offset, orig.vel_scale, orig.vel_bits, natom, pos,
                                2, padded_natom, blockDim.x);
        break;
      case PrecisionModel::SINGLE:
        pos = copyCoordinateSet(dest.xcrd_sp, dest_atom_offset, 0, orig.xvel, orig.xvel_ovrf,
                                orig_atom_offset, orig.vel_scale, orig.vel_bits, natom, pos,
                                0, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.ycrd_sp, dest_atom_offset, 0, orig.yvel, orig.yvel_ovrf,
                                orig_atom_offset, orig.vel_scale, orig.vel_bits, natom, pos,
                                1, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.zcrd_sp, dest_atom_offset, 0, orig.zvel, orig.zvel_ovrf,
                                orig_atom_offset, orig.vel_scale, orig.vel_bits, natom, pos,
                                2, padded_natom, blockDim.x);
        break;
      }
      break;
    case TrajectoryKind::FORCES:
      switch (dest.mode) {
      case PrecisionModel::DOUBLE:
        pos = copyCoordinateSet(dest.xcrd, dest_atom_offset, 0, orig.xfrc, orig.xfrc_ovrf,
                                orig_atom_offset, orig.frc_scale, orig.frc_bits, natom, pos,
                                0, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.ycrd, dest_atom_offset, 0, orig.yfrc, orig.yfrc_ovrf,
                                orig_atom_offset, orig.frc_scale, orig.frc_bits, natom, pos,
                                1, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.zcrd, dest_atom_offset, 0, orig.zfrc, orig.zfrc_ovrf,
                                orig_atom_offset, orig.frc_scale, orig.frc_bits, natom, pos,
                                2, padded_natom, blockDim.x);
        break;
      case PrecisionModel::SINGLE:
        pos = copyCoordinateSet(dest.xcrd_sp, dest_atom_offset, 0, orig.xfrc, orig.xfrc_ovrf,
                                orig_atom_offset, orig.frc_scale, orig.frc_bits, natom, pos,
                                0, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.ycrd_sp, dest_atom_offset, 0, orig.yfrc, orig.yfrc_ovrf,
                                orig_atom_offset, orig.frc_scale, orig.frc_bits, natom, pos,
                                1, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.zcrd_sp, dest_atom_offset, 0, orig.zfrc, orig.zfrc_ovrf,
                                orig_atom_offset, orig.frc_scale, orig.frc_bits, natom, pos,
                                2, padded_natom, blockDim.x);
        break;
      }
      break;
    }
    
    // Increment the system copy count
    pair_idx += gridDim.x;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinates(CondensateWriter *destination, const PsSynthesisReader &origin,
                                  const TrajectoryKind kind, const int2* system_pairs,
                                  const int copy_count, const GpuDetails &gpu) {
  kMultiSystemCoordinateCopy<<<gpu.getSMPCount(), large_block_size>>>(*destination, origin, kind,
                                                                      system_pairs, copy_count);
}

//-------------------------------------------------------------------------------------------------
// Copy multiple systems from a PhaseSpaceSynthesis into a Condensate object.  As with other
// templated coordinate copy kernels, the type reassignment will be done on the CPU, by the
// launching function.  Descriptions of input parameters follow from a previous overloaded variant
// of kMultiSystemCoordinateCopy() above, with changes in origin and destination object types.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kMultiSystemCoordinateCopy(CondensateWriter dest, const CondensateReader orig,
                           const int2* system_pairs, const int copy_count) {
  int pair_idx = blockIdx.x;
  while (pair_idx < copy_count) {
    const int2 origx_to_desty = __ldca(&system_pairs[pair_idx]);
    if (origx_to_desty.x < 0 || origx_to_desty.x >= orig.system_count ||
        origx_to_desty.y < 0 || origx_to_desty.y >= dest.system_count) {
      pair_idx += gridDim.x;
      continue;
    }
    const int natom = __ldca(&dest.atom_counts[origx_to_desty.y]);
    if (natom != __ldca(&orig.atom_counts[origx_to_desty.x])) {
      pair_idx += gridDim.x;
      continue;
    }

    // Assume that box information is relevant and copy it.
    const int xfrm_w = devcRoundUp(9, warp_size_int);
    const int bdim_w = devcRoundUp(6, warp_size_int);
    const int orig_xfrm_offset = origx_to_desty.x * xfrm_w;
    const int dest_xfrm_offset = origx_to_desty.y * xfrm_w;
    const int orig_bdim_offset = origx_to_desty.x * bdim_w;
    const int dest_bdim_offset = origx_to_desty.y * bdim_w;
    kCopyBoxInformation(dest.umat, dest.invu, dest.boxdims, orig.umat, orig.invu, orig.boxdims,
                        dest_xfrm_offset, orig_xfrm_offset, dest_bdim_offset, orig_bdim_offset);

    // Evaluate the nested switch over each object's mode.
    const size_t dest_atom_offset = __ldca(&dest.atom_starts[origx_to_desty.y]);
    const size_t orig_atom_offset = __ldca(&orig.atom_starts[origx_to_desty.x]);
    size_t pos = threadIdx.x;
    const size_t padded_natom = devcRoundUp(natom, warp_size_int);
    switch (dest.mode) {
    case PrecisionModel::DOUBLE:
      switch (orig.mode) {
      case PrecisionModel::DOUBLE:
        pos = copyCoordinateSet(dest.xcrd, dest_atom_offset, 1.0, 0, orig.xcrd, orig_atom_offset,
                                1.0, 0, natom, pos, 0, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.ycrd, dest_atom_offset, 1.0, 0, orig.ycrd, orig_atom_offset,
                                1.0, 0, natom, pos, 1, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.zcrd, dest_atom_offset, 1.0, 0, orig.zcrd, orig_atom_offset,
                                1.0, 0, natom, pos, 2, padded_natom, blockDim.x);
        break;
      case PrecisionModel::SINGLE:
        pos = copyCoordinateSet(dest.xcrd, dest_atom_offset, 1.0, 0, orig.xcrd_sp,
                                orig_atom_offset, 1.0, 0, natom, pos, 0, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.ycrd, dest_atom_offset, 1.0, 0, orig.ycrd_sp,
                                orig_atom_offset, 1.0, 0, natom, pos, 1, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.zcrd, dest_atom_offset, 1.0, 0, orig.zcrd_sp,
                                orig_atom_offset, 1.0, 0, natom, pos, 2, padded_natom, blockDim.x);
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      switch (orig.mode) {
      case PrecisionModel::DOUBLE:
        pos = copyCoordinateSet(dest.xcrd_sp, dest_atom_offset, 1.0, 0, orig.xcrd,
                                orig_atom_offset, 1.0, 0, natom, pos, 0, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.ycrd_sp, dest_atom_offset, 1.0, 0, orig.ycrd,
                                orig_atom_offset, 1.0, 0, natom, pos, 1, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.zcrd_sp, dest_atom_offset, 1.0, 0, orig.zcrd,
                                orig_atom_offset, 1.0, 0, natom, pos, 2, padded_natom, blockDim.x);
        break;
      case PrecisionModel::SINGLE:
        pos = copyCoordinateSet(dest.xcrd_sp, dest_atom_offset, 1.0, 0, orig.xcrd_sp,
                                orig_atom_offset, 1.0, 0, natom, pos, 0, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.ycrd_sp, dest_atom_offset, 1.0, 0, orig.ycrd_sp,
                                orig_atom_offset, 1.0, 0, natom, pos, 1, padded_natom, blockDim.x);
        pos = copyCoordinateSet(dest.zcrd_sp, dest_atom_offset, 1.0, 0, orig.zcrd_sp,
                                orig_atom_offset, 1.0, 0, natom, pos, 2, padded_natom, blockDim.x);
        break;
      }
      break;
    }
    
    // Increment the system copy count
    pair_idx += gridDim.x;    
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchCopyCoordinates(CondensateWriter *destination, const CondensateReader &origin,
                                  const int2* system_pairs, const int copy_count,
                                  const GpuDetails &gpu) {
  kMultiSystemCoordinateCopy<<<gpu.getSMPCount(), large_block_size>>>(*destination, origin,
                                                                      system_pairs, copy_count);
}

} // namespace trajectory
} // namespace stormm
