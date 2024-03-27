// -*-c++-*-
#include <vector>
#ifdef STORMM_USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include "copyright.h"
#include "Constants/scaling.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/stormm_vector_types.h"
#include "Numerics/split_fixed_precision.h"
#include "Reporting/error_format.h"
#include "Trajectory/trajectory_enumerators.h"
#include "hpc_phasespace_synthesis.h"
#include "hpc_phasespace_synthesis.cuh"

namespace stormm {
namespace synthesis {

using numerics::globalpos_scale_nonoverflow_bits;
using numerics::velocity_scale_nonoverflow_bits;
using numerics::force_scale_nonoverflow_bits;
  
#include "Numerics/accumulation.cui"
#include "Math/rounding.cui"
  
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kSystemTransfer(PsSynthesisWriter destination, PsSynthesisWriter source, const int low_index,
                const int high_index, const TrajectoryKind material) {
  __shared__ int sread_start, sread_last_system, sread_last_offset;

  // Pre-compute some constants for transferring positions
  const int mtrx_stride = ((9 + warp_bits_mask_int) / warp_size_int) * warp_size_int;
  const int mtrx_read_start = low_index * mtrx_stride;
  const int mtrx_read_end   = (high_index + 1) * mtrx_stride;
  const int bdim_stride = ((6 + warp_bits_mask_int) / warp_size_int) * warp_size_int;
  const int bdim_read_start = low_index * bdim_stride;
  const int bdim_read_end   = (high_index + 1) * bdim_stride;

  // Assume that the bounds arrays for each system have been copied to the device, where they are
  // visible with the device view pointers.
  if (threadIdx.x == 0) {
    sread_start = source.atom_starts[low_index];
  }
  else if (threadIdx.x == warp_size_int) {
    sread_last_system = source.atom_starts[high_index - 1];
  }
  else if (threadIdx.x == 2 * warp_size_int) {
    sread_last_offset = source.atom_counts[high_index - 1];
  }
  __syncthreads();
  const int atom_read_start = sread_start;
  const int atom_read_end   = sread_last_system + sread_last_offset;

  // The amount of total information transfer is now known.  Assess the number of available warps,
  // assign each an index, and determine a way to subdivide the workload among the available warps.
  int atom_work_units, mtrx_work_units, bdim_work_units, box_work_units, box_work_per_system;
  atom_work_units = (atom_read_end - atom_read_start + warp_bits_mask_int) / warp_size_int;
  switch (material) {
  case TrajectoryKind::POSITIONS:
    {
      mtrx_work_units = mtrx_stride / warp_size_int;
      bdim_work_units = bdim_stride / warp_size_int;
      box_work_per_system = (5 * mtrx_work_units) + (2 * bdim_work_units);
      box_work_units = box_work_per_system * (high_index - low_index);
    }
    break;
  case TrajectoryKind::VELOCITIES:
  case TrajectoryKind::FORCES:
    atom_work_units *= 3;
    box_work_per_system = 0;
    box_work_units = 0;
    break;
  }
  const float total_work_units = (float)(atom_work_units + box_work_units);
  const float total_warps = (float)((blockDim.x / warp_size_int) * gridDim.x);
  int tmp_box_warps = (float)(box_work_units) / total_work_units * (float)total_warps;

  // This setup will not fail on any known GPU, but if the GPU were exceptionally tiny (fewer than
  // seven warps per block and only one block available) it would fail.  But, the GPUs are only
  // getting bigger.
  const int box_warps = (box_work_per_system > 0) ?
                        ((tmp_box_warps + (tmp_box_warps == 0 && box_work_units > 0) +
                          box_work_per_system - 1) / box_work_per_system) * box_work_per_system :
                        0;
  const int atom_warps = total_warps - box_warps;
  
  // Warps now serve atom or box data transfer.
  const int warp_idx = (threadIdx.x >> warp_bits) + ((blockDim.x >> warp_bits) * blockIdx.x);
  const int tgx = (threadIdx.x & warp_bits_mask_int);
  switch (material) {
  case TrajectoryKind::POSITIONS:
    {
      if (warp_idx < atom_warps) {
        int read_pos = atom_read_start + (warp_idx * warp_size_int) + tgx;
        while (read_pos < atom_read_end) {
          destination.xcrd[read_pos] = source.xcrd[read_pos];
          destination.ycrd[read_pos] = source.ycrd[read_pos];
          destination.zcrd[read_pos] = source.zcrd[read_pos];
          destination.xalt[read_pos] = source.xalt[read_pos];
          destination.yalt[read_pos] = source.yalt[read_pos];
          destination.zalt[read_pos] = source.zalt[read_pos];
          if (source.gpos_bits > globalpos_scale_nonoverflow_bits) {
            destination.xcrd_ovrf[read_pos] = source.xcrd_ovrf[read_pos];
            destination.ycrd_ovrf[read_pos] = source.ycrd_ovrf[read_pos];
            destination.zcrd_ovrf[read_pos] = source.zcrd_ovrf[read_pos];
            destination.xalt_ovrf[read_pos] = source.xalt_ovrf[read_pos];
            destination.yalt_ovrf[read_pos] = source.yalt_ovrf[read_pos];
            destination.zalt_ovrf[read_pos] = source.zalt_ovrf[read_pos];
          }
          read_pos += atom_warps * warp_size_int;
        }
      }
      else {
        const int warp_nature = (warp_idx - atom_warps) % box_work_per_system;
        const int box_warp_sets = box_warps / box_work_per_system;
        const int mtrx_warps = mtrx_work_units * box_warp_sets;
        const int bdim_warps = bdim_work_units * box_warp_sets;
        const int box_warp_idx = warp_idx - atom_warps;
        if (warp_nature < mtrx_warps) {
          int read_pos = mtrx_read_start + (box_warp_idx * warp_size_int) + tgx;
          while (read_pos < mtrx_read_end) {
            destination.boxvecs[read_pos] = source.boxvecs[read_pos];
            destination.alt_boxvecs[read_pos] = source.alt_boxvecs[read_pos];
            read_pos += mtrx_warps * warp_size_int;
          }
        }
        else if (warp_nature < 2 * mtrx_warps) {
          int read_pos = mtrx_read_start + ((box_warp_idx - mtrx_warps) * warp_size_int) + tgx;
          while (read_pos < mtrx_read_end) {
            destination.umat[read_pos] = source.umat[read_pos];
            destination.umat_alt[read_pos] = source.umat_alt[read_pos];
            read_pos += mtrx_warps * warp_size_int;
          }
        }
        else if (warp_nature < 3 * mtrx_work_units) {
          int read_pos = mtrx_read_start +
                         ((box_warp_idx - (2 * mtrx_warps)) * warp_size_int) + tgx;
          while (read_pos < mtrx_read_end) {
            destination.invu[read_pos] = source.invu[read_pos];
            destination.invu_alt[read_pos] = source.invu_alt[read_pos];
            read_pos += mtrx_warps * warp_size_int;
          }
        }
        else if (warp_nature < (4 * mtrx_work_units) + bdim_work_units) {
          int read_pos = bdim_read_start +
                         ((box_warp_idx - (3 * mtrx_warps)) * warp_size_int) + tgx;
          while (read_pos < bdim_read_end) {
            destination.boxdims[read_pos] = source.boxdims[read_pos];
            destination.alt_boxdims[read_pos] = source.alt_boxdims[read_pos];
            read_pos += bdim_warps * warp_size_int;
          }
        }
      }
    }
    break;
  case TrajectoryKind::VELOCITIES:
  case TrajectoryKind::FORCES:
    {
      const int pdim_warps = atom_warps / 3;
      if (warp_idx < pdim_warps) {
        int read_pos = atom_read_start + (warp_idx * warp_size_int) + tgx;
        while (read_pos < atom_read_end) {
          if (material == TrajectoryKind::VELOCITIES) {
            destination.xvel[read_pos] = source.xvel[read_pos];
            destination.vxalt[read_pos] = source.vxalt[read_pos];
            if (source.vel_bits > velocity_scale_nonoverflow_bits) {
              destination.xvel_ovrf[read_pos]  = source.xvel_ovrf[read_pos];              
              destination.vxalt_ovrf[read_pos] = source.vxalt_ovrf[read_pos];              
            }
          }
          else {
            destination.xfrc[read_pos] = source.xfrc[read_pos];
            destination.fxalt[read_pos] = source.fxalt[read_pos];
            if (source.frc_bits > force_scale_nonoverflow_bits) {
              destination.xfrc_ovrf[read_pos]  = source.xfrc_ovrf[read_pos];              
              destination.fxalt_ovrf[read_pos] = source.fxalt_ovrf[read_pos];              
            }
          }
          read_pos += pdim_warps * warp_size_int;
        }
      }
      else if (warp_idx < 2 * pdim_warps) {
        int read_pos = atom_read_start + ((warp_idx - pdim_warps) * warp_size_int) + tgx;
        while (read_pos < atom_read_end) {
          if (material == TrajectoryKind::VELOCITIES) {
            destination.yvel[read_pos] = source.yvel[read_pos];
            destination.vyalt[read_pos] = source.vyalt[read_pos];
            if (source.vel_bits > velocity_scale_nonoverflow_bits) {
              destination.yvel_ovrf[read_pos]  = source.yvel_ovrf[read_pos];              
              destination.vyalt_ovrf[read_pos] = source.vyalt_ovrf[read_pos];              
            }
          }
          else {
            destination.yfrc[read_pos] = source.yfrc[read_pos];
            destination.fyalt[read_pos] = source.fyalt[read_pos];
            if (source.frc_bits > force_scale_nonoverflow_bits) {
              destination.yfrc_ovrf[read_pos]  = source.yfrc_ovrf[read_pos];              
              destination.fyalt_ovrf[read_pos] = source.fyalt_ovrf[read_pos];              
            }
          }
          read_pos += pdim_warps * warp_size_int;
        }
      }
      else if (warp_idx < 3 * pdim_warps) {
        int read_pos = atom_read_start + ((warp_idx - (2 * pdim_warps)) * warp_size_int) + tgx;
        while (read_pos < atom_read_end) {
          if (material == TrajectoryKind::VELOCITIES) {
            destination.zvel[read_pos] = source.zvel[read_pos];
            destination.vzalt[read_pos] = source.vzalt[read_pos];
            if (source.vel_bits > velocity_scale_nonoverflow_bits) {
              destination.zvel_ovrf[read_pos]  = source.zvel_ovrf[read_pos];              
              destination.vzalt_ovrf[read_pos] = source.vzalt_ovrf[read_pos];              
            }
          }
          else {
            destination.zfrc[read_pos] = source.zfrc[read_pos];
            destination.fzalt[read_pos] = source.fzalt[read_pos];
            if (source.frc_bits > force_scale_nonoverflow_bits) {
              destination.zfrc_ovrf[read_pos]  = source.zfrc_ovrf[read_pos];              
              destination.fzalt_ovrf[read_pos] = source.fzalt_ovrf[read_pos];              
            }
          }
          read_pos += pdim_warps * warp_size_int;
        }
      }
    }
    break;
  }
}
  
//-------------------------------------------------------------------------------------------------
extern void systemTransfer(PsSynthesisWriter *destination, PsSynthesisWriter *source,
                           const TrajectoryKind kind, const int low_index,
                           const int high_index, const GpuDetails &gpu) {
  kSystemTransfer<<<gpu.getSMPCount(), gpu.getMaxThreadsPerBlock()>>>(*destination, *source,
                                                                      low_index, high_index, kind);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kPsyInitializeForces(PsSynthesisWriter psyw, const int index) {
  int minpos, maxpos;
  if (index < 0) {
    minpos = 0;
    maxpos = psyw.atom_starts[psyw.system_count - 1] + psyw.atom_counts[psyw.system_count - 1];
  }
  else {
    minpos = psyw.atom_starts[index];
    maxpos = minpos + psyw.atom_counts[index];
  }
  for (int pos = minpos + threadIdx.x + (blockDim.x * blockIdx.x);
       pos < maxpos; pos += gridDim.x * blockDim.x) {
    psyw.xfrc[pos] = 0LL;
    psyw.yfrc[pos] = 0LL;
    psyw.zfrc[pos] = 0LL;
    if (psyw.frc_bits > force_scale_nonoverflow_bits) {
      psyw.xfrc_ovrf[pos] = 0;
      psyw.yfrc_ovrf[pos] = 0;
      psyw.zfrc_ovrf[pos] = 0;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void psyInitializeForces(PsSynthesisWriter *psyw, const int index, const GpuDetails &gpu) {
  kPsyInitializeForces<<<gpu.getSMPCount(), gpu.getMaxThreadsPerBlock()>>>(*psyw, index);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kPsyPrimeConjugateGradient(PsSynthesisWriter psyw) {
  const int minpos = threadIdx.x + (blockDim.x * blockIdx.x);
  const int last_system = psyw.system_count - 1;
  const int maxpos = psyw.atom_starts[last_system] + psyw.atom_counts[last_system];
  for (int pos = minpos; pos < maxpos; pos += blockDim.x * gridDim.x) {
    psyw.xalt[pos] = psyw.xfrc[pos];
    psyw.yalt[pos] = psyw.yfrc[pos];
    psyw.zalt[pos] = psyw.zfrc[pos];
    psyw.xvel[pos] = 0LL;
    psyw.yvel[pos] = 0LL;
    psyw.zvel[pos] = 0LL;
    if (psyw.frc_bits > force_scale_nonoverflow_bits) {
      psyw.xalt_ovrf[pos] = psyw.xfrc_ovrf[pos];
      psyw.yalt_ovrf[pos] = psyw.yfrc_ovrf[pos];
      psyw.zalt_ovrf[pos] = psyw.zfrc_ovrf[pos];
      psyw.xvel_ovrf[pos] = 0;
      psyw.yvel_ovrf[pos] = 0;
      psyw.zvel_ovrf[pos] = 0;
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
extern void psyPrimeConjugateGradient(PsSynthesisWriter *psyw, const GpuDetails &gpu) {
  kPsyPrimeConjugateGradient<<<gpu.getSMPCount(), gpu.getMaxThreadsPerBlock()>>>(*psyw);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kPsyCopySystem(PhaseSpaceWriter psw, const PsSynthesisReader poly_psw, const int index) {
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  const size_t atom_offset = poly_psw.atom_starts[index];
  pos = splitFPToReal(psw.xcrd, pos,  0, &poly_psw.xcrd[atom_offset],
                      &poly_psw.xcrd_ovrf[atom_offset], psw.natom, poly_psw.inv_gpos_scale);
  pos = splitFPToReal(psw.ycrd, pos,  1, &poly_psw.ycrd[atom_offset],
                      &poly_psw.ycrd_ovrf[atom_offset], psw.natom, poly_psw.inv_gpos_scale);
  pos = splitFPToReal(psw.zcrd, pos,  2, &poly_psw.zcrd[atom_offset],
                      &poly_psw.zcrd_ovrf[atom_offset], psw.natom, poly_psw.inv_gpos_scale);
  pos = splitFPToReal(psw.xalt, pos,  3, &poly_psw.xalt[atom_offset],
                      &poly_psw.xalt_ovrf[atom_offset], psw.natom, poly_psw.inv_gpos_scale);
  pos = splitFPToReal(psw.yalt, pos,  4, &poly_psw.yalt[atom_offset],
                      &poly_psw.yalt_ovrf[atom_offset], psw.natom, poly_psw.inv_gpos_scale);
  pos = splitFPToReal(psw.zalt, pos,  5, &poly_psw.zalt[atom_offset],
                      &poly_psw.zalt_ovrf[atom_offset], psw.natom, poly_psw.inv_gpos_scale);
  pos = splitFPToReal(psw.xvel, pos,  6, &poly_psw.xvel[atom_offset],
                      &poly_psw.xvel_ovrf[atom_offset], psw.natom, poly_psw.inv_vel_scale);
  pos = splitFPToReal(psw.yvel, pos,  7, &poly_psw.yvel[atom_offset],
                      &poly_psw.yvel_ovrf[atom_offset], psw.natom, poly_psw.inv_vel_scale);
  pos = splitFPToReal(psw.zvel, pos,  8, &poly_psw.zvel[atom_offset],
                      &poly_psw.zvel_ovrf[atom_offset], psw.natom, poly_psw.inv_vel_scale);
  pos = splitFPToReal(psw.vxalt, pos,  9, &poly_psw.vxalt[atom_offset],
                      &poly_psw.vxalt_ovrf[atom_offset], psw.natom, poly_psw.inv_vel_scale);
  pos = splitFPToReal(psw.vyalt, pos, 10, &poly_psw.vyalt[atom_offset],
                      &poly_psw.vyalt_ovrf[atom_offset], psw.natom, poly_psw.inv_vel_scale);
  pos = splitFPToReal(psw.vzalt, pos, 11, &poly_psw.vzalt[atom_offset],
                      &poly_psw.vzalt_ovrf[atom_offset], psw.natom, poly_psw.inv_vel_scale);
  pos = splitFPToReal(psw.xfrc, pos, 12, &poly_psw.xfrc[atom_offset],
                      &poly_psw.xfrc_ovrf[atom_offset], psw.natom, poly_psw.inv_frc_scale);
  pos = splitFPToReal(psw.yfrc, pos, 13, &poly_psw.yfrc[atom_offset],
                      &poly_psw.yfrc_ovrf[atom_offset], psw.natom, poly_psw.inv_frc_scale);
  pos = splitFPToReal(psw.zfrc, pos, 14, &poly_psw.zfrc[atom_offset],
                      &poly_psw.zfrc_ovrf[atom_offset], psw.natom, poly_psw.inv_frc_scale);
  pos = splitFPToReal(psw.fxalt, pos, 15, &poly_psw.fxalt[atom_offset],
                      &poly_psw.fxalt_ovrf[atom_offset], psw.natom, poly_psw.inv_frc_scale);
  pos = splitFPToReal(psw.fyalt, pos, 16, &poly_psw.fyalt[atom_offset],
                      &poly_psw.fyalt_ovrf[atom_offset], psw.natom, poly_psw.inv_frc_scale);
  pos = splitFPToReal(psw.fzalt, pos, 17, &poly_psw.fzalt[atom_offset],
                      &poly_psw.fzalt_ovrf[atom_offset], psw.natom, poly_psw.inv_frc_scale);

  // Copy the transformation matrices and box dimensions
  const int mtrx_offset = devcRoundUp(9, warp_size_int) * index;
  const int bdim_offset = devcRoundUp(6, warp_size_int) * index;
  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (warp_idx == 0) {
    int npos = lane_idx;
    while (npos < 9) {
      psw.umat[npos] = poly_psw.umat[mtrx_offset + npos];
      npos += warp_size_int;
    }
  }
  else if (warp_idx == 1) {
    int npos = lane_idx;
    while (npos < 9) {
      psw.invu[npos] = poly_psw.invu[mtrx_offset + npos];
      npos += warp_size_int;
    }
  }
  else if (warp_idx == 2) {
    int npos = lane_idx;
    while (npos < 6) {
      psw.boxdim[npos] = poly_psw.boxdims[bdim_offset + npos];
      npos += warp_size_int;
    }
  }
  else if (warp_idx == 3) {
    int npos = lane_idx;
    while (npos < 9) {
      psw.umat_alt[npos] = poly_psw.umat_alt[mtrx_offset + npos];
      npos += warp_size_int;
    }
  }
  else if (warp_idx == 4) {
    int npos = lane_idx;
    while (npos < 9) {
      psw.invu_alt[npos] = poly_psw.invu_alt[mtrx_offset + npos];
      npos += warp_size_int;
    }
  }
  else if (warp_idx == 5) {
    int npos = lane_idx;
    while (npos < 6) {
      psw.boxdim_alt[npos] = poly_psw.alt_boxdims[bdim_offset + npos];
      npos += warp_size_int;
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
void PhaseSpaceSynthesis::extractSystem(PhaseSpaceWriter *psw, const int index,
                                        const GpuDetails &gpu,
                                        const HybridTargetLevel origin) const {
  const PsSynthesisReader poly_psr = this->data(origin);
  kPsyCopySystem<<<gpu.getSMPCount(), large_block_size>>>(*psw, poly_psr, index);
  cudaDeviceSynchronize();
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kPsyImportSystemData(llint* x_recv, int* x_recv_ovrf, llint* y_recv, int* y_recv_ovrf,
                     llint* z_recv, int* z_recv_ovrf, double* box_xform, double* inverse_xform,
                     double* box_dimensions, llint* box_vectors, int* box_vector_ovrf,
                     const int* atom_starts, const int* atom_counts, const double* x_import,
                     const double* y_import, const double* z_import, const double* box_xform_in,
                     const double* inverse_xform_in, const double* box_dimensions_in,
                     const int system_index, const TrajectoryKind kind,
                     const double conversion_factor) {

  // Transfer the box information, if it exists
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    if (threadIdx.x < 9) {
      const int box_offset = system_index * devcRoundUp(9, warp_size_int);
      box_xform[box_offset + threadIdx.x] = box_xform_in[threadIdx.x];
      inverse_xform[box_offset + threadIdx.x] = inverse_xform_in[threadIdx.x];
    }
    if (threadIdx.x >= warp_size_int && threadIdx.x < warp_size_int + 9) {
      const int box_offset = system_index * devcRoundUp(9, warp_size_int);
      const int95_t fpbv = doubleToInt95(inverse_xform_in[threadIdx.x] * conversion_factor);
      box_vectors[box_offset + threadIdx.x] = fpbv.x;
      box_vector_ovrf[box_offset + threadIdx.x] = fpbv.y;
    }
    if (threadIdx.x >= twice_warp_size_int && threadIdx.x < twice_warp_size_int + 6) {
      const int dim_offset = system_index * devcRoundUp(9, warp_size_int);
      box_dimensions[dim_offset + threadIdx.x] = box_dimensions_in[threadIdx.x];
    }
    break;
  case TrajectoryKind::VELOCITIES:
  case TrajectoryKind::FORCES:
    break;
  }

  // Transfer atomic data
  const int pos_start = atom_starts[system_index];
  const int pos_end   = pos_start + atom_counts[system_index];
  const int stride = pos_end - pos_start;
  const int padded_stride = devcRoundUp(stride, warp_size_int);
  int pos = threadIdx.x;
  while (pos < padded_stride) {
    if (pos < stride) {
      const int95_t fpx = doubleToInt95((double)(x_import[pos]) * conversion_factor);
      const size_t ip = pos + pos_start;
      x_recv[ip]      = fpx.x;
      x_recv_ovrf[ip] = fpx.y;
    }
    pos += blockDim.x * gridDim.x;
  }
  while (pos < 2 * padded_stride) {
    const int rel_pos = pos - padded_stride;
    if (rel_pos < padded_stride) {
      const int95_t fpy = doubleToInt95((double)(y_import[rel_pos]) * conversion_factor);
      const size_t ip = rel_pos + pos_start;
      y_recv[ip]      = fpy.x;
      y_recv_ovrf[ip] = fpy.y;
    }
    pos += blockDim.x * gridDim.x;
  }
  while (pos < 3 * padded_stride) {
    const int rel_pos = pos - (2 * padded_stride);
    if (rel_pos < padded_stride) {
      const int95_t fpz = doubleToInt95((double)(z_import[rel_pos]) * conversion_factor);
      const size_t ip = rel_pos + pos_start;
      z_recv[ip]      = fpz.x;
      z_recv_ovrf[ip] = fpz.y;
    }
    pos += blockDim.x * gridDim.x;
  }
}
  
//-------------------------------------------------------------------------------------------------
void psyImportSystemData(llint* x_recv, int* x_recv_ovrf, llint* y_recv, int* y_recv_ovrf,
                         llint* z_recv, int* z_recv_ovrf, double* box_xform, double* inverse_xform,
                         double* box_dimensions, llint* box_vectors, int* box_vector_ovrf,
                         const int* atom_starts, const int* atom_counts, const double* x_import,
                         const double* y_import, const double* z_import,
                         const double* box_xform_in, const double* inverse_xform_in,
                         const double* box_dimensions_in, const int system_index,
                         const TrajectoryKind kind, const double conversion_factor,
                         const GpuDetails &gpu) {
  kPsyImportSystemData<<<large_block_size,
                         gpu.getSMPCount()>>>(x_recv, x_recv_ovrf, y_recv, y_recv_ovrf, z_recv,
                                              z_recv_ovrf, box_xform, inverse_xform,
                                              box_dimensions, box_vectors, box_vector_ovrf,
                                              atom_starts, atom_counts, x_import, y_import,
                                              z_import, box_xform_in, inverse_xform_in,
                                              box_dimensions_in, system_index, kind,
                                              conversion_factor);
}

} // namespace synthesis
} // namespace stormm
