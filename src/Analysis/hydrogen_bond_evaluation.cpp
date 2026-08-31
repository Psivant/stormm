#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/fixed_precision.h"
#include "Constants/symbol_values.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/rounding.h"
#include "MolecularMechanics/dynamics_intervention.h"
#include "Structure/local_arrangement.h"
#include "Structure/structure_enumerators.h"
#include "Topology/atomgraph_enumerators.h"
#include "hydrogen_bond_analysis.h"
#ifdef STORMM_USE_HPC
#include "hpc_hydrogen_bond_analysis.h"
#endif

namespace stormm {
namespace mm {

using symbols::pi_f;
using symbols::twopi_f;
using numerics::globalpos_scale_nonoverflow_bits;
using stmath::roundUp;
using structure::angle;
using structure::imageCoordinates;
using structure::ImagingMethod;
using synthesis::PsSynthesisReader;
using topology::UnitCellType;
  
//-------------------------------------------------------------------------------------------------
void evalHydrogenBondAnalysis(const int step, const int index, const HybridTargetLevel tier) {
  HBondWriter hbw = dyna_tk.getHydrogenBondAnalysisData(index, tier);
  if (step < hbw.init_step || (step % hbw.eval_intv) != 0) {
    return;
  }
  const PsSynthesisReader poly_psr = dyna_tk.getReadOnlyCoordinateData(tier);
  const int xfrm_stride = roundUp(9, warp_size_int);
  switch (tier) {
  case HybridTargetLevel::HOST:
    for (int i = 0; i < hbw.n_partners; i++) {
      const int4 particle_ids = hbw.partners[i];

      // This function, as other analysis functions, is designed to fire off at the end of the
      // simulation time step but before the coordinate objects have updated their positions in
      // the coordinate cycle.  Evaluate the state of the most recently constructed coordinates,
      // those assembled at the end of the time step.  This function, like other analyses, will
      // do the majority of its work in single-precision (float32_t) and use fixed-precision
      // accumulation, but in order to set up the problem with imaging considerations, the local
      // coordinate frame will be established with double-precision (float64_t) calculations.
      double t_dx, t_dy, t_dz, t_ax, t_ay, t_az;
      const size_t xfrm_offset = xfrm_stride * particle_ids.w;
      if (poly_psr.gpos_bits > globalpos_scale_nonoverflow_bits) {
        const int95_t i_hdx = hostInt95Subtract(poly_psr.xalt[particle_ids.x],
                                                poly_psr.xalt_ovrf[particle_ids.x],
                                                poly_psr.xalt[particle_ids.y],
                                                poly_psr.xalt_ovrf[particle_ids.y]);
        const int95_t i_hdy = hostInt95Subtract(poly_psr.yalt[particle_ids.x],
                                                poly_psr.yalt_ovrf[particle_ids.x],
                                                poly_psr.yalt[particle_ids.y],
                                                poly_psr.yalt_ovrf[particle_ids.y]);
        const int95_t i_hdz = hostInt95Subtract(poly_psr.zalt[particle_ids.x],
                                                poly_psr.zalt_ovrf[particle_ids.x],
                                                poly_psr.zalt[particle_ids.y],
                                                poly_psr.zalt_ovrf[particle_ids.y]);
        t_dx = hostInt95ToDouble(i_hdx) * poly_psr.inv_gpos_scale;
        t_dy = hostInt95ToDouble(i_hdy) * poly_psr.inv_gpos_scale;
        t_dz = hostInt95ToDouble(i_hdz) * poly_psr.inv_gpos_scale;
        const int95_t i_hax = hostInt95Subtract(poly_psr.xalt[particle_ids.z],
                                                poly_psr.xalt_ovrf[particle_ids.z],
                                                poly_psr.xalt[particle_ids.y],
                                                poly_psr.xalt_ovrf[particle_ids.y]);
        const int95_t i_hay = hostInt95Subtract(poly_psr.yalt[particle_ids.z],
                                                poly_psr.yalt_ovrf[particle_ids.z],
                                                poly_psr.yalt[particle_ids.y],
                                                poly_psr.yalt_ovrf[particle_ids.y]);
        const int95_t i_haz = hostInt95Subtract(poly_psr.zalt[particle_ids.z],
                                                poly_psr.zalt_ovrf[particle_ids.z],
                                                poly_psr.zalt[particle_ids.y],
                                                poly_psr.zalt_ovrf[particle_ids.y]);
        t_ax = hostInt95ToDouble(i_hax) * poly_psr.inv_gpos_scale;
        t_ay = hostInt95ToDouble(i_hay) * poly_psr.inv_gpos_scale;
        t_az = hostInt95ToDouble(i_haz) * poly_psr.inv_gpos_scale;
      }
      else {
        const llint i_hdx = poly_psr.xalt[particle_ids.x] - poly_psr.xalt[particle_ids.y];
        const llint i_hdy = poly_psr.yalt[particle_ids.x] - poly_psr.yalt[particle_ids.y];
        const llint i_hdz = poly_psr.zalt[particle_ids.x] - poly_psr.zalt[particle_ids.y];
        t_dx = static_cast<double>(i_hdx) * poly_psr.inv_gpos_scale;
        t_dy = static_cast<double>(i_hdy) * poly_psr.inv_gpos_scale;
        t_dz = static_cast<double>(i_hdz) * poly_psr.inv_gpos_scale;
        const llint i_hax = poly_psr.xalt[particle_ids.z] - poly_psr.xalt[particle_ids.y];
        const llint i_hay = poly_psr.yalt[particle_ids.z] - poly_psr.yalt[particle_ids.y];
        const llint i_haz = poly_psr.zalt[particle_ids.z] - poly_psr.zalt[particle_ids.y];
        t_ax = static_cast<double>(i_hax) * poly_psr.inv_gpos_scale;
        t_ay = static_cast<double>(i_hay) * poly_psr.inv_gpos_scale;
        t_az = static_cast<double>(i_haz) * poly_psr.inv_gpos_scale;
      }
      
      // Image the displacements and compute the donor-acceptor distance
      imageCoordinates<double, double>(&t_dx, &t_dy, &t_dz, &poly_psr.umat[xfrm_offset],
                                       &poly_psr.invu[xfrm_offset], poly_psr.unit_cell,
                                       ImagingMethod::MINIMUM_IMAGE);
      imageCoordinates<double, double>(&t_ax, &t_ay, &t_az, &poly_psr.umat[xfrm_offset],
                                       &poly_psr.invu[xfrm_offset], poly_psr.unit_cell,
                                       ImagingMethod::MINIMUM_IMAGE);
      const float local_x[3] = { static_cast<float>(t_dx), 0.0f, static_cast<float>(t_ax) };
      const float local_y[3] = { static_cast<float>(t_dy), 0.0f, static_cast<float>(t_ay) };
      const float local_z[3] = { static_cast<float>(t_dz), 0.0f, static_cast<float>(t_az) };
      const float disp_x = local_x[2] - local_x[0];
      const float disp_y = local_y[2] - local_y[0];
      const float disp_z = local_z[2] - local_z[0];
      const float t_range = sqrtf((disp_x * disp_x) + (disp_y * disp_y) + (disp_z * disp_z));
      const int statblk_idx = (step - hbw.init_step) / hbw.block_steps;
      const int item_idx = i + (statblk_idx * hbw.n_partners);
      if (t_range < hbw.max_separation) {

        // If the particles are close enough to verify a hydrogen bond, check the angle.
        float t_angle = angle<float, float>(0, 1, 2, local_x, local_y, local_z,
                                            &poly_psr.umat[xfrm_offset],
                                            &poly_psr.invu[xfrm_offset],
                                            UnitCellType::NONE);
        float t_segm = t_angle / twopi_f;
        if (t_angle > pi_f) {
          t_angle -= ceilf(t_segm) * twopi_f;
        }
        else if (t_angle < -pi_f) {
          t_angle += floorf(t_segm) * twopi_f;
        }
        t_angle = fabsf(t_angle);
        if (t_angle >= hbw.min_angle) {

          // The arrangement is verified as a hydrogen bond.  Increment the tally and various
          // accumulators.  On the GPU, one thread will process each pair of possible hydrogen
          // bonding partners, averting problems of asynchronous accumulation.
          hbw.tallies[item_idx] += 1;
          hbw.dist_acc[item_idx] += t_range;
          hbw.dist_sqacc[item_idx] += (t_range * t_range);
          hbw.angl_acc[item_idx] += t_angle;
          hbw.angl_sqacc[item_idx] += (t_angle * t_angle);
        }
      }

      // Always log the distance between the donor and acceptor, so that bonds which do manifest
      // themselves can be seen in context of the variability of the donor-acceptor distance over
      // time.
      hbw.all_dist_acc[item_idx] += t_range;
      hbw.all_dist_sqacc[item_idx] += (t_range * t_range);
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    launchEvalHydrogenBondAnalysis(step, &hbw, poly_psr);
    break;
#endif
  }
}

} // namespace mm
} // namespace stormm
