#include "copyright.h"
#include "Constants/fixed_precision.h"
#include "Topology/atomgraph_abstracts.h"
#ifdef STORMM_USE_HPC
#  include "hpc_trim.h"
#endif
#include "trim.h"

namespace stormm {
namespace trajectory {

using numerics::globalpos_scale_nonoverflow_bits;
using numerics::velocity_scale_nonoverflow_bits;
using synthesis::PsSynthesisWriter;
using synthesis::SyAtomUpdateKit;
using topology::ChemicalDetailsKit;
  
//-------------------------------------------------------------------------------------------------
void removeMomentum(PhaseSpace *ps, const AtomGraph *ag, const ExceptionResponse policy) {
  PhaseSpaceWriter psw = ps->data();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  removeMomentum<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, nullptr, nullptr, nullptr,
                                         psw.xvel, psw.yvel, psw.zvel, nullptr, nullptr, nullptr,
                                         cdk.masses, psw.unit_cell, psw.natom, 1.0, 1.0, policy);
}

//-------------------------------------------------------------------------------------------------
void removeMomentum(PhaseSpace *ps, const AtomGraph &ag, const ExceptionResponse policy) {
  removeMomentum(ps, ag.getSelfPointer(), policy);
}

//-------------------------------------------------------------------------------------------------
void removeMomentum(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                    const PrecisionModel prec, const ExceptionResponse policy) {

  // Check that the system counts are consistent.
  PsSynthesisWriter host_poly_psw = poly_ps->data();
  if (poly_ag->getSystemCount() != host_poly_psw.system_count) {
    rtErr("Inconsistent system counts (" + std::to_string(poly_ag->getSystemCount()) + " and " +
          std::to_string(host_poly_psw.system_count) + ") were found in the topology and "
          "coordinate syntheses.", "removeMomentum");
  }

  // Check that the atom counts, and atom alignments are consistent between the coordinate and
  // topology syntheses.  This check ensures that, if the kernels will refuse to work on misaligned
  // data, the program will halt.  Performing the check here keeps the GPU from having to wait on
  // slow CPU operations.  The CPU can keep working once the kernel is launched.  If CPU-based
  // momentum removal is requested, then these checks will happen ahead of the operations.
  for (int i = 0; i < host_poly_psw.system_count; i++) {
    if (poly_ag->getAtomOffset(i) != host_poly_psw.atom_starts[i] ||
        poly_ag->getAtomCount(i) != host_poly_psw.atom_counts[i]) {
      rtErr("Inconsistent atom counts (" + std::to_string(poly_ag->getAtomOffset(i)) + ", " +
            std::to_string(host_poly_psw.atom_starts[i]) + ") or alignments in the overall "
            "list (" + std::to_string(poly_ag->getAtomCount(i)) + ", " +
            std::to_string(host_poly_psw.atom_counts[i]) + ") were found in the topology and "
            "coordinate syntheses.", "removeMomentum");
    }
  }

  // Loop over all systems and delegate to the templated function.
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyAtomUpdateKit<double,
                            double2,
                            double4> poly_auk = poly_ag->getDoublePrecisionAtomUpdateKit();
      for (int i = 0; i < host_poly_psw.system_count; i++) {
        const size_t sysi_offset = host_poly_psw.atom_starts[i];
        removeMomentum<llint,
                       double, double>(&host_poly_psw.xcrd[sysi_offset],
                                       &host_poly_psw.ycrd[sysi_offset],
                                       &host_poly_psw.zcrd[sysi_offset],
                                       &host_poly_psw.xcrd_ovrf[sysi_offset],
                                       &host_poly_psw.ycrd_ovrf[sysi_offset],
                                       &host_poly_psw.zcrd_ovrf[sysi_offset],
                                       &host_poly_psw.xvel[sysi_offset],
                                       &host_poly_psw.yvel[sysi_offset],
                                       &host_poly_psw.zvel[sysi_offset],
                                       &host_poly_psw.xvel_ovrf[sysi_offset],
                                       &host_poly_psw.yvel_ovrf[sysi_offset],
                                       &host_poly_psw.zvel_ovrf[sysi_offset],
                                       &poly_auk.masses[sysi_offset], poly_ag->getUnitCellType(),
                                       host_poly_psw.atom_counts[i], host_poly_psw.gpos_scale_f,
                                       host_poly_psw.vel_scale_f, policy);
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyAtomUpdateKit<float,
                            float2, float4> poly_auk = poly_ag->getSinglePrecisionAtomUpdateKit();
      for (int i = 0; i < host_poly_psw.system_count; i++) {
        const size_t sysi_offset = host_poly_psw.atom_starts[i];
        removeMomentum<llint,
                       float, float>(&host_poly_psw.xcrd[sysi_offset],
                                     &host_poly_psw.ycrd[sysi_offset],
                                     &host_poly_psw.zcrd[sysi_offset], nullptr, nullptr, nullptr,
                                     &host_poly_psw.xvel[sysi_offset],
                                     &host_poly_psw.yvel[sysi_offset],
                                     &host_poly_psw.zvel[sysi_offset], nullptr, nullptr, nullptr,
                                     &poly_auk.masses[sysi_offset], poly_ag->getUnitCellType(),
                                     host_poly_psw.atom_counts[i], host_poly_psw.gpos_scale_f,
                                     host_poly_psw.vel_scale_f, policy);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void removeMomentum(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                    const PrecisionModel prec, const ExceptionResponse policy) {
  removeMomentum(poly_ps, poly_ag.getSelfPointer(), prec, policy);
}

//-------------------------------------------------------------------------------------------------
void removeMomentum(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                    MotionSweeper *mos, const GpuDetails &gpu, const ExceptionResponse policy) {

  // Check to ensure that the topology and coordinate syntheses match.
  const int nsys = poly_ps->getSystemCount();
  if (poly_ag.getSystemCount() != nsys) {
    rtErr("The number of systems in the coordinate synthesis (" + std::to_string(nsys) + ") must "
          "match that in the topology synthesis (" + std::to_string(poly_ag.getSystemCount()) +
          ").", "MotionSweeper", "accumulateCenterOfMass");
  }
  for (int i = 0; i < nsys; i++) {
    if (poly_ps->getAtomCount(i) != poly_ag.getAtomCount(i)) {
      rtErr("The atom count of system " + std::to_string(i) + " of the coordinate "
            "synthesis (" + std::to_string(poly_ps->getAtomCount(i)) + ") does not equal that "
            "of the topology synthesis (" + std::to_string(poly_ag.getAtomCount(i)) + ").",
            "MotionSweeper", "accumulateCenterOfMass");
    }
  }
  const SyAtomUpdateKit<double,
                        double2, double4> poly_auk = poly_ag.getDoublePrecisionAtomUpdateKit();
  if (gpu == null_gpu) {
    PsSynthesisWriter poly_psw = poly_ps->data();
    MotionSweepWriter mosw = mos->data();
    const PsSynthesisReader poly_psr(poly_psw);
    const MotionSweepReader mosr(mosw);
    accumulateCenterOfMassMotion(&mosw, poly_auk, poly_psr);
    removeCenterOfMassMotion(&poly_psw, mosr);
    accumulateAngularMomentum(&mosw, poly_auk, poly_psr);
    removeAngularMomentum(&poly_psw, mosr, gpu, policy);
  }
#ifdef STORMM_USE_HPC
  else {
    PsSynthesisWriter poly_psw = poly_ps->data(HybridTargetLevel::DEVICE);
    MotionSweepWriter mosw = mos->data(HybridTargetLevel::DEVICE);
    const PsSynthesisReader poly_psr(poly_psw);
    const MotionSweepReader mosr(mosw);
    accumulateCenterOfMassMotion(&mosw, poly_auk, poly_psr, gpu);
    removeCenterOfMassMotion(&poly_psw, mosr, gpu);
    accumulateAngularMomentum(&mosw, poly_auk, poly_psr, gpu);
    removeAngularMomentum(&poly_psw, mosr, gpu, policy);
  }
#endif
}

//-------------------------------------------------------------------------------------------------
void accumulateCenterOfMassMotion(MotionSweepWriter *mosw,
                                  const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                                  const PsSynthesisReader &poly_psr, const GpuDetails &gpu) {
  if (gpu == null_gpu) {

    // Momentum removal, and center of mass calculations, will take place just prior to force
    // calculations and therefore on the current (not the developing) coordinates of the
    // synthesis.
    for (int i = 0; i < poly_psr.system_count; i++) {

      // The serial CPU routine can initialize accumulators directly before filling them.
      mosw->xcom[i] = 0LL;
      mosw->xcom_ovrf[i] = 0;
      mosw->ycom[i] = 0LL;
      mosw->ycom_ovrf[i] = 0;
      mosw->zcom[i] = 0LL;
      mosw->zcom_ovrf[i] = 0;
      mosw->xmv[i] = 0LL;
      mosw->xmv_ovrf[i] = 0;
      mosw->ymv[i] = 0LL;
      mosw->ymv_ovrf[i] = 0;
      mosw->zmv[i] = 0LL;
      mosw->zmv_ovrf[i] = 0;
      const int llim = poly_psr.atom_starts[i];
      const int hlim = llim + poly_psr.atom_counts[i];
      for (int j = llim; j < hlim; j++) {

        // Center of mass
        const double cfac = poly_psr.inv_gpos_scale * mosw->com_scale * poly_auk.masses[j];
        double xcntrb, ycntrb, zcntrb;
        if (poly_psr.gpos_bits <= globalpos_scale_nonoverflow_bits) {
          xcntrb = static_cast<double>(poly_psr.xcrd[j]) * cfac;
          ycntrb = static_cast<double>(poly_psr.ycrd[j]) * cfac;
          zcntrb = static_cast<double>(poly_psr.zcrd[j]) * cfac;
        }
        else {
          xcntrb = hostInt95ToDouble(poly_psr.xcrd[j], poly_psr.xcrd_ovrf[j]) * cfac;
          ycntrb = hostInt95ToDouble(poly_psr.ycrd[j], poly_psr.ycrd_ovrf[j]) * cfac;
          zcntrb = hostInt95ToDouble(poly_psr.zcrd[j], poly_psr.zcrd_ovrf[j]) * cfac;
        }
        hostSplitAccumulation(xcntrb, &mosw->xcom[i], &mosw->xcom_ovrf[i]);
        hostSplitAccumulation(ycntrb, &mosw->ycom[i], &mosw->ycom_ovrf[i]);
        hostSplitAccumulation(zcntrb, &mosw->zcom[i], &mosw->zcom_ovrf[i]);

        // Translational momentum
        const double mvfac = poly_psr.inv_vel_scale * mosw->mv_scale * poly_auk.masses[j];
        double xtm, ytm, ztm;
        if (poly_psr.vel_bits <= velocity_scale_nonoverflow_bits) {
          xtm = static_cast<double>(poly_psr.xvel[j]);
          ytm = static_cast<double>(poly_psr.yvel[j]);
          ztm = static_cast<double>(poly_psr.zvel[j]);
        }
        else {
          xtm = hostInt95ToDouble(poly_psr.xvel[j], poly_psr.xvel_ovrf[j]);
          ytm = hostInt95ToDouble(poly_psr.yvel[j], poly_psr.yvel_ovrf[j]);
          ztm = hostInt95ToDouble(poly_psr.zvel[j], poly_psr.zvel_ovrf[j]);
        }
        hostSplitAccumulation(xtm * mvfac, &mosw->xmv[i], &mosw->xmv_ovrf[i]);
        hostSplitAccumulation(ytm * mvfac, &mosw->ymv[i], &mosw->ymv_ovrf[i]);
        hostSplitAccumulation(ztm * mvfac, &mosw->zmv[i], &mosw->zmv_ovrf[i]);
      }
    }
  }
#ifdef STORMM_USE_HPC
  else {
    launchAccCenterOfMassMotion(mosw, poly_auk, poly_psr, gpu);
  }
#endif
}

//-------------------------------------------------------------------------------------------------
void removeCenterOfMassMotion(PsSynthesisWriter *poly_psw, const MotionSweepReader &mosr,
                              const GpuDetails &gpu) {
  if (gpu == null_gpu) {
    for (int i = 0; i < poly_psw->system_count; i++) {
      const int llim = poly_psw->atom_starts[i];
      const int hlim = llim + poly_psw->atom_counts[i];

      // Move the center of mass to the origin
      const double com_fac = poly_psw->gpos_scale / (mosr.com_scale * mosr.total_mass[i]);
      const double com_x = hostInt95ToDouble(mosr.xcom[i], mosr.xcom_ovrf[i]) * com_fac;
      const double com_y = hostInt95ToDouble(mosr.ycom[i], mosr.ycom_ovrf[i]) * com_fac;
      const double com_z = hostInt95ToDouble(mosr.zcom[i], mosr.zcom_ovrf[i]) * com_fac;
      if (poly_psw->gpos_bits <= globalpos_scale_nonoverflow_bits) {
        const llint iadj_x = llround(-com_x);
        const llint iadj_y = llround(-com_y);
        const llint iadj_z = llround(-com_z);
        for (int j = llim; j < hlim; j++) {
          poly_psw->xcrd[j] += iadj_x;
          poly_psw->ycrd[j] += iadj_y;
          poly_psw->zcrd[j] += iadj_z;
        }
      }
      else {
        const int95_t iadj_x = hostDoubleToInt95(-com_x);
        const int95_t iadj_y = hostDoubleToInt95(-com_y);
        const int95_t iadj_z = hostDoubleToInt95(-com_z);
        for (int j = llim; j < hlim; j++) {
          const int95_t inx =  hostSplitFPSum(iadj_x, poly_psw->xcrd[j], poly_psw->xcrd_ovrf[j]);
          const int95_t iny =  hostSplitFPSum(iadj_y, poly_psw->ycrd[j], poly_psw->ycrd_ovrf[j]);
          const int95_t inz =  hostSplitFPSum(iadj_z, poly_psw->zcrd[j], poly_psw->zcrd_ovrf[j]);
          poly_psw->xcrd[j] = inx.x;
          poly_psw->ycrd[j] = iny.x;
          poly_psw->zcrd[j] = inz.x;
          poly_psw->xcrd_ovrf[j] = inx.y;
          poly_psw->ycrd_ovrf[j] = iny.y;
          poly_psw->zcrd_ovrf[j] = inz.y;
        }
      }
      
      // Remove the system's net velocity
      const double mvs_fac = poly_psw->vel_scale / (mosr.mv_scale * mosr.total_mass[i]);
      const double net_xvel = hostInt95ToDouble(mosr.xmv[i], mosr.xmv_ovrf[i]) * mvs_fac;
      const double net_yvel = hostInt95ToDouble(mosr.ymv[i], mosr.ymv_ovrf[i]) * mvs_fac;
      const double net_zvel = hostInt95ToDouble(mosr.zmv[i], mosr.zmv_ovrf[i]) * mvs_fac;
      if (poly_psw->vel_bits <= velocity_scale_nonoverflow_bits) {
        const llint iadj_xvel = llround(-net_xvel);
        const llint iadj_yvel = llround(-net_yvel);
        const llint iadj_zvel = llround(-net_zvel);
        for (int j = llim; j < hlim; j++) {
          poly_psw->xvel[j] += iadj_xvel;
          poly_psw->yvel[j] += iadj_yvel;
          poly_psw->zvel[j] += iadj_zvel;
        }
      }
      else {
        const int95_t iadj_xvel = hostDoubleToInt95(-net_xvel);
        const int95_t iadj_yvel = hostDoubleToInt95(-net_yvel);
        const int95_t iadj_zvel = hostDoubleToInt95(-net_zvel);
        for (int j = llim; j < hlim; j++) {
          const int95_t inx = hostSplitFPSum(iadj_xvel, poly_psw->xvel[j], poly_psw->xvel_ovrf[j]);
          const int95_t iny = hostSplitFPSum(iadj_yvel, poly_psw->yvel[j], poly_psw->yvel_ovrf[j]);
          const int95_t inz = hostSplitFPSum(iadj_zvel, poly_psw->zvel[j], poly_psw->zvel_ovrf[j]);
          poly_psw->xvel[j] = inx.x;
          poly_psw->yvel[j] = iny.x;
          poly_psw->zvel[j] = inz.x;
          poly_psw->xvel_ovrf[j] = inx.y;
          poly_psw->yvel_ovrf[j] = iny.y;
          poly_psw->zvel_ovrf[j] = inz.y;
        }
      }
    }
  }
#ifdef STORMM_USE_HPC
  else {
    launchRemoveCenterOfMassMotion(poly_psw, mosr, gpu);
  }
#endif
}
  
//-------------------------------------------------------------------------------------------------
void accumulateAngularMomentum(MotionSweepWriter *mosw,
                               const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                               const PsSynthesisReader &poly_psr, const GpuDetails &gpu) {

  // Angular momentum is only applicable for systems without periodic boundary conditions.
  switch (poly_psr.unit_cell) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    return;
  }
  if (gpu == null_gpu) {

    // Momentum removal, and center of mass calculations, will take place just prior to force
    // calculations and therefore on the current (not the developing) coordinates of the
    // synthesis.
    for (int i = 0; i < poly_psr.system_count; i++) {

      // As with other accumulators, the serial CPU routine can initialize at the beginning
      // of the loop over each system.  The GPU must initialize out-of-place.
      mosw->rxmv[i] = 0LL;
      mosw->rymv[i] = 0LL;
      mosw->rzmv[i] = 0LL;
      mosw->rxmv_ovrf[i] = 0;
      mosw->rymv_ovrf[i] = 0;
      mosw->rzmv_ovrf[i] = 0;

      // Some local accumulators for the inertial tensor components will make the code cleaner and
      // faster, analogous to the GPU kernel keeping these accumulators in __shared__ memory.
      llint prim_xx = 0LL;
      llint prim_xy = 0LL;
      llint prim_xz = 0LL;
      llint prim_yy = 0LL;
      llint prim_yz = 0LL;
      llint prim_zz = 0LL;
      int ovrf_xx = 0;
      int ovrf_xy = 0;
      int ovrf_xz = 0;
      int ovrf_yy = 0;
      int ovrf_yz = 0;
      int ovrf_zz = 0;
      const int llim = poly_psr.atom_starts[i];
      const int hlim = llim + poly_psr.atom_counts[i];
      double r[3], v[3], rcv[3];
      for (int j = llim; j < hlim; j++) {
        if (poly_psr.gpos_bits <= globalpos_scale_nonoverflow_bits) {
          r[0] = static_cast<double>(poly_psr.xcrd[j]);
          r[1] = static_cast<double>(poly_psr.ycrd[j]);
          r[2] = static_cast<double>(poly_psr.zcrd[j]);
        }
        else {
          r[0] = hostInt95ToDouble(poly_psr.xcrd[j], poly_psr.xcrd_ovrf[j]);
          r[1] = hostInt95ToDouble(poly_psr.ycrd[j], poly_psr.ycrd_ovrf[j]);
          r[2] = hostInt95ToDouble(poly_psr.zcrd[j], poly_psr.zcrd_ovrf[j]);
        }
        if (poly_psr.vel_bits <= velocity_scale_nonoverflow_bits) {
          v[0] = static_cast<double>(poly_psr.xvel[j]);
          v[1] = static_cast<double>(poly_psr.yvel[j]);
          v[2] = static_cast<double>(poly_psr.zvel[j]);
        }
        else {
          v[0] = hostInt95ToDouble(poly_psr.xvel[j], poly_psr.xvel_ovrf[j]);
          v[1] = hostInt95ToDouble(poly_psr.yvel[j], poly_psr.yvel_ovrf[j]);
          v[2] = hostInt95ToDouble(poly_psr.zvel[j], poly_psr.zvel_ovrf[j]);
        }
        for (int k = 0; k < 3; k++) {
          r[k] *= poly_psr.inv_gpos_scale;
          v[k] *= poly_psr.inv_vel_scale;
        }
        const double inrt_fac = mosw->inrt_scale * poly_auk.masses[j];
        const double mv_fac = mosw->mv_scale * poly_auk.masses[j];
        hostSplitAccumulation(r[0] * r[0] * inrt_fac,
                              &mosw->inrt[ 6 * i     ], &mosw->inrt_ovrf[ 6 * i     ]);
        hostSplitAccumulation(r[0] * r[1] * inrt_fac,
                              &mosw->inrt[(6 * i) + 1], &mosw->inrt_ovrf[(6 * i) + 1]);
        hostSplitAccumulation(r[0] * r[2] * inrt_fac,
                              &mosw->inrt[(6 * i) + 2], &mosw->inrt_ovrf[(6 * i) + 2]);
        hostSplitAccumulation(r[1] * r[1] * inrt_fac,
                              &mosw->inrt[(6 * i) + 3], &mosw->inrt_ovrf[(6 * i) + 3]);
        hostSplitAccumulation(r[1] * r[2] * inrt_fac,
                              &mosw->inrt[(6 * i) + 4], &mosw->inrt_ovrf[(6 * i) + 4]);
        hostSplitAccumulation(r[2] * r[2] * inrt_fac,
                              &mosw->inrt[(6 * i) + 5], &mosw->inrt_ovrf[(6 * i) + 5]);
        crossProduct(r, v, rcv);
        hostSplitAccumulation(rcv[0] * mv_fac, &mosw->rxmv[i], &mosw->rxmv_ovrf[i]);
        hostSplitAccumulation(rcv[1] * mv_fac, &mosw->rymv[i], &mosw->rymv_ovrf[i]);
        hostSplitAccumulation(rcv[2] * mv_fac, &mosw->rzmv[i], &mosw->rzmv_ovrf[i]);
        hostSplitAccumulation(r[0] * r[0] * inrt_fac, &prim_xx, &ovrf_xx);
        hostSplitAccumulation(r[0] * r[1] * inrt_fac, &prim_xy, &ovrf_xy);
        hostSplitAccumulation(r[0] * r[2] * inrt_fac, &prim_xz, &ovrf_xz);
        hostSplitAccumulation(r[1] * r[1] * inrt_fac, &prim_yy, &ovrf_yy);
        hostSplitAccumulation(r[1] * r[2] * inrt_fac, &prim_yz, &ovrf_yz);
        hostSplitAccumulation(r[2] * r[2] * inrt_fac, &prim_zz, &ovrf_zz);
      }
      mosw->inrt[(6 * i)    ] = prim_xx;
      mosw->inrt[(6 * i) + 1] = prim_xy;
      mosw->inrt[(6 * i) + 2] = prim_xz;
      mosw->inrt[(6 * i) + 3] = prim_yy;
      mosw->inrt[(6 * i) + 4] = prim_yz;
      mosw->inrt[(6 * i) + 5] = prim_zz;
      mosw->inrt_ovrf[(6 * i)    ] = ovrf_xx;
      mosw->inrt_ovrf[(6 * i) + 1] = ovrf_xy;
      mosw->inrt_ovrf[(6 * i) + 2] = ovrf_xz;
      mosw->inrt_ovrf[(6 * i) + 3] = ovrf_yy;
      mosw->inrt_ovrf[(6 * i) + 4] = ovrf_yz;
      mosw->inrt_ovrf[(6 * i) + 5] = ovrf_zz;
    }
  }
#ifdef STORMM_USE_HPC
  else {
    launchAccAngularMomentum(mosw, poly_auk, poly_psr, gpu);
  }
#endif
}

//-------------------------------------------------------------------------------------------------
void removeAngularMomentum(PsSynthesisWriter *poly_psw, const MotionSweepReader &mosr,
                           const GpuDetails &gpu, ExceptionResponse policy) {
  switch (poly_psw->unit_cell) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    return;
  }
  if (gpu == null_gpu) {
    std::vector<double> cmps(6), tnsr(9), inv_tnsr(9);
    for (int i = 0; i < poly_psw->system_count; i++) {
      const int llim = poly_psw->atom_starts[i];
      const int hlim = llim + poly_psw->atom_counts[i];

      // Construct and invert the inertial tensor.  On the GPU, this will be handled by a each warp
      // ahead of all other work: a relatively inefficient mechanism, but the best available.
      for (int j = 0; j < 6; j++) {
        const llint prim = mosr.inrt[(6 * i) + j];
        const int ovrf = mosr.inrt_ovrf[(6 * i) + j];
        cmps[j] = hostInt95ToDouble(prim, ovrf) / mosr.inrt_scale;
      }
      tnsr[0] =  cmps[3] + cmps[5];
      tnsr[1] = -cmps[1];
      tnsr[2] = -cmps[2];
      tnsr[3] = -cmps[1];
      tnsr[4] =  cmps[0] + cmps[5];
      tnsr[5] = -cmps[4];
      tnsr[6] = -cmps[2];
      tnsr[7] = -cmps[4];
      tnsr[8] =  cmps[0] + cmps[3];

      // Invert the inertial tensor, if possible.  Compute the net rotational velocity and
      // remove it as well.
      if (checkInertialTensor(tnsr, hlim - llim, policy) == false) {
        continue;
      }
      invertSquareMatrix(tnsr, &inv_tnsr);
      const double t_rxmv = hostInt95ToDouble(mosr.rxmv[i], mosr.rxmv_ovrf[i]) / mosr.mv_scale;
      const double t_rymv = hostInt95ToDouble(mosr.rymv[i], mosr.rymv_ovrf[i]) / mosr.mv_scale;
      const double t_rzmv = hostInt95ToDouble(mosr.rzmv[i], mosr.rzmv_ovrf[i]) / mosr.mv_scale;
      double rot_vel[3], r[3], vcr[3];
      rot_vel[0] = (inv_tnsr[0] * t_rxmv) + (inv_tnsr[3] * t_rymv) + (inv_tnsr[6] * t_rzmv);
      rot_vel[1] = (inv_tnsr[1] * t_rxmv) + (inv_tnsr[4] * t_rymv) + (inv_tnsr[7] * t_rzmv);
      rot_vel[2] = (inv_tnsr[2] * t_rxmv) + (inv_tnsr[5] * t_rymv) + (inv_tnsr[8] * t_rzmv);
      for (int j = llim; j < hlim; j++) {
        if (poly_psw->gpos_bits <= globalpos_scale_nonoverflow_bits) {
          r[0] = poly_psw->xcrd[j];
          r[1] = poly_psw->ycrd[j];
          r[2] = poly_psw->zcrd[j];
        }
        else {
          r[0] = hostInt95ToDouble(poly_psw->xcrd[j], poly_psw->xcrd_ovrf[j]);
          r[1] = hostInt95ToDouble(poly_psw->ycrd[j], poly_psw->ycrd_ovrf[j]);
          r[2] = hostInt95ToDouble(poly_psw->zcrd[j], poly_psw->zcrd_ovrf[j]);
        }
        r[0] *= poly_psw->inv_gpos_scale;
        r[1] *= poly_psw->inv_gpos_scale;
        r[2] *= poly_psw->inv_gpos_scale;
        crossProduct(rot_vel, r, vcr);
        if (poly_psw->vel_bits <= velocity_scale_nonoverflow_bits) {
          poly_psw->xvel[j] -= llround(vcr[0] * poly_psw->vel_scale);
          poly_psw->yvel[j] -= llround(vcr[1] * poly_psw->vel_scale);
          poly_psw->zvel[j] -= llround(vcr[2] * poly_psw->vel_scale);
        }
        else {
          const int95_t inx = hostInt95Sum(poly_psw->xvel[j], poly_psw->xvel_ovrf[j],
                                           -vcr[0] * poly_psw->vel_scale);
          const int95_t iny = hostInt95Sum(poly_psw->yvel[j], poly_psw->yvel_ovrf[j],
                                           -vcr[1] * poly_psw->vel_scale);
          const int95_t inz = hostInt95Sum(poly_psw->zvel[j], poly_psw->zvel_ovrf[j],
                                           -vcr[2] * poly_psw->vel_scale);
          poly_psw->xvel[j] = inx.x;
          poly_psw->yvel[j] = iny.x;
          poly_psw->zvel[j] = inz.x;
          poly_psw->xvel_ovrf[j] = inx.y;
          poly_psw->yvel_ovrf[j] = iny.y;
          poly_psw->zvel_ovrf[j] = inz.y;
        }
      }
    }
  }
#ifdef STORMM_USE_HPC
  else {
    launchRemoveAngularMomentum(poly_psw, mosr, gpu);
  }
#endif
}

} // namespace trajectory
} // namespace stormm
