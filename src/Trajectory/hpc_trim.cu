// -*-c++-*-
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Constants/hpc_bounds.h"
#include "Constants/fixed_precision.h"
#include "hpc_trim.h"

namespace stormm {
namespace trajectory {

#include "Accelerator/syncwarp.cui"
#include "Numerics/accumulation.cui"
#include "Math/matrix_formulas.cui"
#include "Math/vector_formulas.cui"

using numerics::globalpos_scale_nonoverflow_bits;
using numerics::velocity_scale_nonoverflow_bits;
  
/// \brief Accumulate motion of the center of mass for each system in the synthesis, including its
///        position relative to the origin and total momentum.
///
/// \param mosw      Contains accumulators for center of mass, momentum, and moment of inertia
/// \param poly_auk  Contains masses of all particles in all systems
/// \param poly_psr  Contains coordinates and velocities of all particles in all systems
__global__ void __launch_bounds__(small_block_size, 1)
kAccumulateCenterOfMassMotion(MotionSweepWriter mosw,
                              const SyAtomUpdateKit<double, double2, double4> poly_auk,
                              const PsSynthesisReader poly_psr) {
  __shared__ volatile double sh_xcom[small_block_size / warp_size_int];
  __shared__ volatile double sh_ycom[small_block_size / warp_size_int];
  __shared__ volatile double sh_zcom[small_block_size / warp_size_int];
  __shared__ volatile double sh_xmv[small_block_size / warp_size_int];
  __shared__ volatile double sh_ymv[small_block_size / warp_size_int];
  __shared__ volatile double sh_zmv[small_block_size / warp_size_int];
  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  int wu_idx = blockIdx.x;
  while (wu_idx < mosw.nwu) {
    const int4 wu = mosw.work_units[wu_idx];
    
    // If this work unit is the first for the system, initialize the accumulators for the next
    // iteration.
    if (wu.w == 1 && threadIdx.x == 0) {
      mosw.xcom_nxt[wu.z] = 0LL;
      mosw.ycom_nxt[wu.z] = 0LL;
      mosw.zcom_nxt[wu.z] = 0LL;
      mosw.xcom_nxt_ovrf[wu.z] = 0;
      mosw.ycom_nxt_ovrf[wu.z] = 0;
      mosw.zcom_nxt_ovrf[wu.z] = 0;
      mosw.xmv_nxt[wu.z] = 0LL;
      mosw.ymv_nxt[wu.z] = 0LL;
      mosw.zmv_nxt[wu.z] = 0LL;
      mosw.xmv_nxt_ovrf[wu.z] = 0;
      mosw.ymv_nxt_ovrf[wu.z] = 0;
      mosw.zmv_nxt_ovrf[wu.z] = 0;
    }

    // Center of mass
    double xcntrb = 0.0;
    double ycntrb = 0.0;
    double zcntrb = 0.0;
    for (int i = wu.x + threadIdx.x; i < wu.y; i += blockDim.x) {
      const double cfac = poly_psr.inv_gpos_scale * mosw.com_scale * poly_auk.masses[i];
      if (poly_psr.gpos_bits <= globalpos_scale_nonoverflow_bits) {
        xcntrb += (double)(poly_psr.xcrd[i]) * cfac;
        ycntrb += (double)(poly_psr.ycrd[i]) * cfac;
        zcntrb += (double)(poly_psr.zcrd[i]) * cfac;
      }
      else {
        xcntrb += int95ToDouble(poly_psr.xcrd[i], poly_psr.xcrd_ovrf[i]) * cfac;
        ycntrb += int95ToDouble(poly_psr.ycrd[i], poly_psr.ycrd_ovrf[i]) * cfac;
        zcntrb += int95ToDouble(poly_psr.zcrd[i], poly_psr.zcrd_ovrf[i]) * cfac;
      }
    }
    WARP_REDUCE_DOWN(xcntrb);
    WARP_REDUCE_DOWN(ycntrb);
    WARP_REDUCE_DOWN(zcntrb);
    if (lane_idx == 0) {
      sh_xcom[warp_idx] = xcntrb;
      sh_ycom[warp_idx] = ycntrb;
      sh_zcom[warp_idx] = zcntrb;
    }

    // Total momentum accumulation
    double xtm = 0.0;
    double ytm = 0.0;
    double ztm = 0.0;
    for (int i = wu.x + threadIdx.x; i < wu.y; i += blockDim.x) {
      const double cfac = poly_psr.inv_vel_scale * mosw.mv_scale * poly_auk.masses[i];
      if (poly_psr.gpos_bits <= globalpos_scale_nonoverflow_bits) {
        xtm += (double)(poly_psr.xvel[i]) * cfac;
        ytm += (double)(poly_psr.yvel[i]) * cfac;
        ztm += (double)(poly_psr.zvel[i]) * cfac;
      }
      else {
        xtm += int95ToDouble(poly_psr.xvel[i], poly_psr.xvel_ovrf[i]) * cfac;
        ytm += int95ToDouble(poly_psr.yvel[i], poly_psr.yvel_ovrf[i]) * cfac;
        ztm += int95ToDouble(poly_psr.zvel[i], poly_psr.zvel_ovrf[i]) * cfac;
      }
    }
    WARP_REDUCE_DOWN(xtm);
    WARP_REDUCE_DOWN(ytm);
    WARP_REDUCE_DOWN(ztm);
    if (lane_idx == 0) {
      sh_xmv[warp_idx] = xtm;
      sh_ymv[warp_idx] = ytm;
      sh_zmv[warp_idx] = ztm;
    }

    // Merge the work of each warp in the block, then use atomics to add the work unit's
    // contribution to each global accumulator.
    __syncthreads();
    if (warp_idx == 0) {
      double n_xcntrb, n_xtm;
      if (lane_idx < small_block_size / warp_size_int) {
        n_xcntrb = sh_xcom[lane_idx];
        n_xtm = sh_xmv[lane_idx];
      }
      else {
        n_xcntrb = 0.0;
        n_xtm = 0.0;
      }
      WARP_REDUCE_DOWN(n_xcntrb);
      WARP_REDUCE_DOWN(n_xtm);
      if (lane_idx == 0) {
        atomicSplit(n_xcntrb, wu.z, mosw.xcom, mosw.xcom_ovrf);
        atomicSplit(n_xtm, wu.z, mosw.xmv, mosw.xmv_ovrf);
      }
    }
    else if (warp_idx == 1) {
      double n_ycntrb, n_ytm;
      if (lane_idx < small_block_size / warp_size_int) {
        n_ycntrb = sh_ycom[lane_idx];
        n_ytm = sh_ymv[lane_idx];
      }
      else {
        n_ycntrb = 0.0;
        n_ytm = 0.0;
      }
      WARP_REDUCE_DOWN(n_ycntrb);
      WARP_REDUCE_DOWN(n_ytm);
      if (lane_idx == 0) {
        atomicSplit(n_ycntrb, wu.z, mosw.ycom, mosw.ycom_ovrf);
        atomicSplit(n_ytm, wu.z, mosw.ymv, mosw.ymv_ovrf);
      }
    }
    else if (warp_idx == 2) {
      double n_zcntrb, n_ztm;
      if (lane_idx < small_block_size / warp_size_int) {
        n_zcntrb = sh_zcom[lane_idx];
        n_ztm = sh_zmv[lane_idx];
      }
      else {
        n_zcntrb = 0.0;
        n_ztm = 0.0;
      }
      WARP_REDUCE_DOWN(n_zcntrb);
      WARP_REDUCE_DOWN(n_ztm);
      if (lane_idx == 0) {
        atomicSplit(n_zcntrb, wu.z, mosw.zcom, mosw.zcom_ovrf);
        atomicSplit(n_ztm, wu.z, mosw.zmv, mosw.zmv_ovrf);
      }
    }
    __syncthreads();
    wu_idx += gridDim.x;
  }
}
  
//-------------------------------------------------------------------------------------------------
void launchAccCenterOfMassMotion(MotionSweepWriter *mosw,
                                 const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                                 const PsSynthesisReader &poly_psr, const GpuDetails &gpu) {
  const int nblock = gpu.getSMPCount() * (gpu.getMaxThreadsPerBlock() / small_block_size);
  kAccumulateCenterOfMassMotion<<<nblock, small_block_size>>>(*mosw, poly_auk, poly_psr);
}

/// \brief Remove translational motion of the center of mass.
///
/// \param mosr      Contains accumulators for center of mass, momentum, and moment of inertia
/// \param poly_auk  Contains masses of all particles in all systems
/// \param poly_psr  Contains coordinates and velocities of all particles in all systems
__global__ void __launch_bounds__(small_block_size, 1)
kRemoveCenterOfMassMotion(PsSynthesisWriter poly_psw, const MotionSweepReader mosr) {
  int wu_idx = blockIdx.x;
  while (wu_idx < mosr.nwu) {
    const int4 wu = mosr.work_units[wu_idx];

    // Move the system's center of mass to the origin
    const double com_fac = poly_psw.gpos_scale / (mosr.com_scale * mosr.total_mass[wu.z]);
    const double com_x = int95ToDouble(mosr.xcom[wu.z], mosr.xcom_ovrf[wu.z]) * com_fac;
    const double com_y = int95ToDouble(mosr.ycom[wu.z], mosr.ycom_ovrf[wu.z]) * com_fac;
    const double com_z = int95ToDouble(mosr.zcom[wu.z], mosr.zcom_ovrf[wu.z]) * com_fac;
    if (poly_psw.gpos_bits <= globalpos_scale_nonoverflow_bits) {
      const llint iadj_x = __double2ll_rn(-com_x);
      const llint iadj_y = __double2ll_rn(-com_y);
      const llint iadj_z = __double2ll_rn(-com_z);
      for (int i = wu.x + threadIdx.x; i < wu.y; i += blockDim.x) {
        poly_psw.xcrd[i] += iadj_x;
        poly_psw.ycrd[i] += iadj_y;
        poly_psw.zcrd[i] += iadj_z;
      }
    }
    else {
      const int95_t iadj_x = doubleToInt95(-com_x);
      const int95_t iadj_y = doubleToInt95(-com_y);
      const int95_t iadj_z = doubleToInt95(-com_z);
      for (int i = wu.x + threadIdx.x; i < wu.y; i += blockDim.x) {
        const int95_t inx = splitFPSum(iadj_x, poly_psw.xcrd[i], poly_psw.xcrd_ovrf[i]);
        const int95_t iny = splitFPSum(iadj_y, poly_psw.ycrd[i], poly_psw.ycrd_ovrf[i]);
        const int95_t inz = splitFPSum(iadj_z, poly_psw.zcrd[i], poly_psw.zcrd_ovrf[i]);
        poly_psw.xcrd[i] = inx.x;
        poly_psw.ycrd[i] = iny.x;
        poly_psw.zcrd[i] = inz.x;
        poly_psw.xcrd_ovrf[i] = inx.y;
        poly_psw.ycrd_ovrf[i] = iny.y;
        poly_psw.zcrd_ovrf[i] = inz.y;
      }
    }

    // Remove the system's net velocity
    const double mvs_fac = poly_psw.vel_scale / (mosr.mv_scale * mosr.total_mass[wu.z]);
    const double net_xvel = int95ToDouble(mosr.xmv[wu.z], mosr.xmv_ovrf[wu.z]) * mvs_fac;
    const double net_yvel = int95ToDouble(mosr.ymv[wu.z], mosr.ymv_ovrf[wu.z]) * mvs_fac;
    const double net_zvel = int95ToDouble(mosr.zmv[wu.z], mosr.zmv_ovrf[wu.z]) * mvs_fac;
    if (poly_psw.vel_bits <= velocity_scale_nonoverflow_bits) {
      const llint iadj_xvel = __double2ll_rn(-net_xvel);
      const llint iadj_yvel = __double2ll_rn(-net_yvel);
      const llint iadj_zvel = __double2ll_rn(-net_zvel);
      for (int i = wu.x + threadIdx.x; i < wu.y; i += blockDim.x) {
        poly_psw.xvel[i] += iadj_xvel;
        poly_psw.yvel[i] += iadj_yvel;
        poly_psw.zvel[i] += iadj_zvel;
      }
    }
    else {
      const int95_t iadj_xvel = doubleToInt95(-net_xvel);
      const int95_t iadj_yvel = doubleToInt95(-net_yvel);
      const int95_t iadj_zvel = doubleToInt95(-net_zvel);
      for (int i = wu.x + threadIdx.x; i < wu.y; i += blockDim.x) {
        const int95_t inx = splitFPSum(iadj_xvel, poly_psw.xvel[i], poly_psw.xvel_ovrf[i]);
        const int95_t iny = splitFPSum(iadj_yvel, poly_psw.yvel[i], poly_psw.yvel_ovrf[i]);
        const int95_t inz = splitFPSum(iadj_zvel, poly_psw.zvel[i], poly_psw.zvel_ovrf[i]);
        poly_psw.xvel[i] = inx.x;
        poly_psw.yvel[i] = iny.x;
        poly_psw.zvel[i] = inz.x;
        poly_psw.xvel_ovrf[i] = inx.y;
        poly_psw.yvel_ovrf[i] = iny.y;
        poly_psw.zvel_ovrf[i] = inz.y;
      }
    }
    wu_idx += gridDim.x;
  }
}

//-------------------------------------------------------------------------------------------------
void launchRemoveCenterOfMassMotion(PsSynthesisWriter *poly_psw, const MotionSweepReader &mosr,
                                    const GpuDetails &gpu) {
  const int nblock = gpu.getSMPCount() * (gpu.getMaxThreadsPerBlock() / small_block_size);
  kRemoveCenterOfMassMotion<<<nblock, small_block_size>>>(*poly_psw, mosr);
}

/// \brief Accumulate angular momentum based on recentered systems with zero net momentum.
///
/// \param mosw      Contains accumulators for center of mass, momentum, and moment of inertia
/// \param poly_auk  Contains masses of all particles in all systems
/// \param poly_psr  Contains coordinates and velocities of all particles in all systems
__global__ void __launch_bounds__(small_block_size, 1)
kAccumulateAngularMomentum(MotionSweepWriter mosw,
                           SyAtomUpdateKit<double, double2, double4> poly_auk,
                           const PsSynthesisReader poly_psr) {
  __shared__ volatile double sh_inrt_xx[small_block_size / warp_size_int];
  __shared__ volatile double sh_inrt_xy[small_block_size / warp_size_int];
  __shared__ volatile double sh_inrt_xz[small_block_size / warp_size_int];
  __shared__ volatile double sh_inrt_yy[small_block_size / warp_size_int];
  __shared__ volatile double sh_inrt_yz[small_block_size / warp_size_int];
  __shared__ volatile double sh_inrt_zz[small_block_size / warp_size_int];
  __shared__ volatile double sh_inrt_holdings[6];
  __shared__ volatile double sh_rx_mv[small_block_size / warp_size_int];
  __shared__ volatile double sh_ry_mv[small_block_size / warp_size_int];
  __shared__ volatile double sh_rz_mv[small_block_size / warp_size_int];

  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  int wu_idx = blockIdx.x;
  while (wu_idx < mosw.nwu) {
    const int4 wu = mosw.work_units[wu_idx];

    // If this work unit is the first for the system, initialize the accumulators for the next
    // iteration.
    if (wu.w == 1) {
      if (threadIdx.x == 0) {
        mosw.rxmv_nxt[wu.z] = 0LL;
        mosw.rymv_nxt[wu.z] = 0LL;
        mosw.rzmv_nxt[wu.z] = 0LL;
        mosw.rxmv_nxt_ovrf[wu.z] = 0;
        mosw.rymv_nxt_ovrf[wu.z] = 0;
        mosw.rzmv_nxt_ovrf[wu.z] = 0;
      }
      if (warp_idx == 1 && lane_idx < 6) {
        mosw.inrt_nxt[lane_idx] = 0LL;
        mosw.inrt_nxt_ovrf[lane_idx] = 0;
      }
    }

    // Compute the cross-product of centered positions and velocities
    double inrt_xx = 0.0;
    double inrt_xy = 0.0;
    double inrt_xz = 0.0;
    double inrt_yy = 0.0;
    double inrt_yz = 0.0;
    double inrt_zz = 0.0;
    double rx_mv = 0.0;
    double ry_mv = 0.0;
    double rz_mv = 0.0;
    const double rmv_fac = poly_psr.inv_gpos_scale * poly_psr.inv_vel_scale * mosw.mv_scale;
    const double inrt_fac = poly_psr.inv_gpos_scale * poly_psr.inv_gpos_scale * mosw.inrt_scale;
    for (int i = wu.x + threadIdx.x; i < wu.y; i += blockDim.x) {
      double3 r, v;
      if (poly_psr.gpos_bits <= globalpos_scale_nonoverflow_bits) {
        r.x = (double)(poly_psr.xcrd[i]);
        r.y = (double)(poly_psr.ycrd[i]);
        r.z = (double)(poly_psr.zcrd[i]);
      }
      else {
        r.x = int95ToDouble(poly_psr.xcrd[i], poly_psr.xcrd_ovrf[i]);
        r.y = int95ToDouble(poly_psr.ycrd[i], poly_psr.ycrd_ovrf[i]);
        r.z = int95ToDouble(poly_psr.zcrd[i], poly_psr.zcrd_ovrf[i]);
      }
      if (poly_psr.gpos_bits <= velocity_scale_nonoverflow_bits) {
        v.x = (double)(poly_psr.xvel[i]);
        v.y = (double)(poly_psr.yvel[i]);
        v.z = (double)(poly_psr.zvel[i]);
      }
      else {
        v.x = int95ToDouble(poly_psr.xvel[i], poly_psr.xvel_ovrf[i]);
        v.y = int95ToDouble(poly_psr.yvel[i], poly_psr.yvel_ovrf[i]);
        v.z = int95ToDouble(poly_psr.zvel[i], poly_psr.zvel_ovrf[i]);
      }
      const double m_inrt_fac = inrt_fac * poly_auk.masses[i];
      inrt_xx += r.x * r.x * m_inrt_fac;
      inrt_xy += r.x * r.y * m_inrt_fac;
      inrt_xz += r.x * r.z * m_inrt_fac;
      inrt_yy += r.y * r.y * m_inrt_fac;
      inrt_yz += r.y * r.z * m_inrt_fac;
      inrt_zz += r.z * r.z * m_inrt_fac;
      const double3 rcv = crossProduct(r, v);
      const double m_rmv_fac = rmv_fac * poly_auk.masses[i];
      rx_mv += rcv.x * m_rmv_fac;
      ry_mv += rcv.y * m_rmv_fac;
      rz_mv += rcv.z * m_rmv_fac;
    }
    WARP_REDUCE_DOWN(inrt_xx);
    WARP_REDUCE_DOWN(inrt_xy);
    WARP_REDUCE_DOWN(inrt_xz);
    WARP_REDUCE_DOWN(inrt_yy);
    WARP_REDUCE_DOWN(inrt_yz);
    WARP_REDUCE_DOWN(inrt_zz);
    WARP_REDUCE_DOWN(rx_mv);
    WARP_REDUCE_DOWN(ry_mv);
    WARP_REDUCE_DOWN(rz_mv);
    if (lane_idx == 0) {
      sh_inrt_xx[warp_idx] = inrt_xx;
      sh_inrt_xy[warp_idx] = inrt_xy;
      sh_inrt_xz[warp_idx] = inrt_xz;
      sh_inrt_yy[warp_idx] = inrt_yy;
      sh_inrt_yz[warp_idx] = inrt_yz;
      sh_inrt_zz[warp_idx] = inrt_zz;
      sh_rx_mv[warp_idx] = rx_mv;
      sh_ry_mv[warp_idx] = ry_mv;
      sh_rz_mv[warp_idx] = rz_mv;
    }
    __syncthreads();
    if (warp_idx == 0) {
      if (lane_idx < small_block_size / warp_size_int) {
        rx_mv = sh_rx_mv[lane_idx];
        inrt_xx = sh_inrt_xx[lane_idx];
        inrt_xy = sh_inrt_xy[lane_idx];
      }
      else {
        rx_mv = 0.0;
        inrt_xx = 0.0;
        inrt_xy = 0.0;
      }
      WARP_REDUCE_DOWN(rx_mv);
      WARP_REDUCE_DOWN(inrt_xx);
      WARP_REDUCE_DOWN(inrt_xy);
      if (lane_idx == 0) {
        atomicSplit(rx_mv, wu.z, mosw.rxmv, mosw.rxmv_ovrf);
        sh_inrt_holdings[0] = inrt_xx;
        sh_inrt_holdings[1] = inrt_xy;
      }
    }
    else if (warp_idx == 1) {
      if (lane_idx < small_block_size / warp_size_int) {
        ry_mv = sh_ry_mv[lane_idx];
        inrt_xz = sh_inrt_xz[lane_idx];
        inrt_yy = sh_inrt_yy[lane_idx];
      }
      else {
        ry_mv = 0.0;
        inrt_xz = 0.0;
        inrt_yy = 0.0;
      }
      WARP_REDUCE_DOWN(ry_mv);
      WARP_REDUCE_DOWN(inrt_xz);
      WARP_REDUCE_DOWN(inrt_yy);
      if (lane_idx == 0) {
        atomicSplit(ry_mv, wu.z, mosw.rymv, mosw.rymv_ovrf);
        sh_inrt_holdings[2] = inrt_xz;
        sh_inrt_holdings[3] = inrt_yy;
      }
    }
    else if (warp_idx == 2) {
      if (lane_idx < small_block_size / warp_size_int) {
        rz_mv = sh_rz_mv[lane_idx];
        inrt_yz = sh_inrt_yz[lane_idx];
        inrt_zz = sh_inrt_zz[lane_idx];
      }
      else {
        rz_mv = 0.0;
        inrt_yz = 0.0;
        inrt_zz = 0.0;
      }
      WARP_REDUCE_DOWN(rz_mv);
      WARP_REDUCE_DOWN(inrt_yz);
      WARP_REDUCE_DOWN(inrt_zz);
      if (lane_idx == 0) {
        atomicSplit(rz_mv, wu.z, mosw.rzmv, mosw.rzmv_ovrf);
        sh_inrt_holdings[4] = inrt_yz;
        sh_inrt_holdings[5] = inrt_zz;
      }
    }
    __syncthreads();
    if (warp_idx == 0) {
      if (lane_idx < 6) {
        atomicSplit(sh_inrt_holdings[lane_idx], (6 * wu.z) + lane_idx, mosw.inrt, mosw.inrt_ovrf);
      }
    }
    __syncthreads();
    wu_idx += gridDim.x;
  }
}

//-------------------------------------------------------------------------------------------------
void launchAccAngularMomentum(MotionSweepWriter *mosw,
                              const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                              const PsSynthesisReader &poly_psr, const GpuDetails &gpu) {
  const int nblock = gpu.getSMPCount() * (gpu.getMaxThreadsPerBlock() / small_block_size);
  kAccumulateAngularMomentum<<<nblock, small_block_size>>>(*mosw, poly_auk, poly_psr);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(small_block_size, 1)
kRemoveAngularMomentum(PsSynthesisWriter poly_psw, const MotionSweepReader mosr) {
  __shared__ volatile double sh_inrt_elem[6], sh_inrt[9], sh_cof[9], sh_inv[9]; 

  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  int wu_idx = blockIdx.x;
  while (wu_idx < mosr.nwu) {
    const int4 wu = mosr.work_units[wu_idx];

    // Construct and invert the inertial tensor for this system.
    if (warp_idx == 0) {
      if (lane_idx < 6) {
        sh_inrt_elem[lane_idx] = mosr.inrt[(6 * wu.z) + lane_idx] / mosr.inrt_scale;
      }
      SYNCWARP;
      if (lane_idx < 9) {
        int idx_a, idx_b;
        if (lane_idx == 0) {
          idx_a = 3;
          idx_b = 5;
        }
        else if (lane_idx == 4) {
          idx_a = 0;
          idx_b = 5;
        }
        else if (lane_idx == 8) {
          idx_a = 0;
          idx_b = 3;
        }
        else if (lane_idx == 1 || lane_idx == 3) {
          idx_a = 1;
        }
        else if (lane_idx == 2 || lane_idx == 6) {
          idx_a = 2;
        }
        else if (lane_idx == 5 || lane_idx == 7) {
          idx_a = 4;
        }
        if ((lane_idx & 0x3) == 0) {
          sh_inrt[lane_idx] = sh_inrt_elem[idx_a] + sh_inrt_elem[idx_b];
        }
        else {
          sh_inrt[lane_idx] = -sh_inrt_elem[idx_a];
        }
      }
      SYNCWARP;
      invertRankThreeMatrix((const double*)(sh_inrt), (double*)(sh_cof), (double*)(sh_inv));
    }
    __syncthreads();
    const double t_rxmv = int95ToDouble(mosr.rxmv[wu.z], mosr.rxmv_ovrf[wu.z]) / mosr.mv_scale;
    const double t_rymv = int95ToDouble(mosr.rymv[wu.z], mosr.rymv_ovrf[wu.z]) / mosr.mv_scale;
    const double t_rzmv = int95ToDouble(mosr.rzmv[wu.z], mosr.rzmv_ovrf[wu.z]) / mosr.mv_scale;
    const double3 rot_vel = { (sh_inv[0] * t_rxmv) + (sh_inv[3] * t_rymv) + (sh_inv[6] * t_rzmv),
                              (sh_inv[1] * t_rxmv) + (sh_inv[4] * t_rymv) + (sh_inv[7] * t_rzmv),
                              (sh_inv[2] * t_rxmv) + (sh_inv[5] * t_rymv) + (sh_inv[8] * t_rzmv) };
    for (int i = wu.x + threadIdx.x; i < wu.y; i += blockDim.x) {
      double3 r;
      if (poly_psw.gpos_bits <= globalpos_scale_nonoverflow_bits) {
        r.x = poly_psw.xcrd[i];
        r.y = poly_psw.ycrd[i];
        r.z = poly_psw.zcrd[i];
      }
      else {
        r.x = int95ToDouble(poly_psw.xcrd[i], poly_psw.xcrd_ovrf[i]);
        r.y = int95ToDouble(poly_psw.ycrd[i], poly_psw.ycrd_ovrf[i]);
        r.z = int95ToDouble(poly_psw.zcrd[i], poly_psw.zcrd_ovrf[i]);
      }
      r.x *= poly_psw.inv_gpos_scale;
      r.y *= poly_psw.inv_gpos_scale;
      r.z *= poly_psw.inv_gpos_scale;
      const double3 vcr = crossProduct(rot_vel, r);
      if (poly_psw.vel_bits <= velocity_scale_nonoverflow_bits) {
        poly_psw.xvel[i] -= __double2ll_rn(vcr.x * poly_psw.vel_scale);
        poly_psw.yvel[i] -= __double2ll_rn(vcr.y * poly_psw.vel_scale);
        poly_psw.zvel[i] -= __double2ll_rn(vcr.z * poly_psw.vel_scale);
      }
      else {
        const int95_t inx = int95Sum(poly_psw.xvel[i], poly_psw.xvel_ovrf[i],
                                     -vcr.x * poly_psw.vel_scale);
        const int95_t iny = int95Sum(poly_psw.yvel[i], poly_psw.yvel_ovrf[i],
                                     -vcr.y * poly_psw.vel_scale);
        const int95_t inz = int95Sum(poly_psw.zvel[i], poly_psw.zvel_ovrf[i],
                                     -vcr.z * poly_psw.vel_scale);
        poly_psw.xvel[i] = inx.x;
        poly_psw.yvel[i] = iny.x;
        poly_psw.zvel[i] = inz.x;
        poly_psw.xvel_ovrf[i] = inx.y;
        poly_psw.yvel_ovrf[i] = iny.y;
        poly_psw.zvel_ovrf[i] = inz.y;
      }
    }
    wu_idx += gridDim.x;
  }
}

//-------------------------------------------------------------------------------------------------
void launchRemoveAngularMomentum(PsSynthesisWriter *poly_psw, const MotionSweepReader &mosr,
                                 const GpuDetails &gpu) {
  const int nblock = gpu.getSMPCount() * (gpu.getMaxThreadsPerBlock() / small_block_size);
  kRemoveAngularMomentum<<<nblock, small_block_size>>>(*poly_psw, mosr);
}

} // namespace trajectory
} // namespace stormm
