// -*-c++-*-
#include "copyright.h"
#include "Constants/fixed_precision.h"
#include "MolecularMechanics/dynamics_intervention.h"
#include "Numerics/split_fixed_precision.h"
#include "Topology/atomgraph_enumerators.h"
#include "hpc_hydrogen_bond_analysis.h"

namespace stormm {
namespace analysis {

using mm::dyna_tk;
using numerics::globalpos_scale_nonoverflow_bits;
using topology::UnitCellType;
  
#include "Math/rounding.cui"
#include "Numerics/accumulation.cui"
#include "Structure/local_arrangement.cui"
  
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(medium_block_size, 2)
kEvalHydrogenBondAnalysis(const int step, HBondWriter hbw, const PsSynthesisReader poly_psr) {
  const int grid_dim = blockDim.x * gridDim.x;
  const int xfrm_stride = devcRoundUp(9, warp_size_int);
  for (int pos = threadIdx.x + (blockIdx.x * blockDim.x); pos < hbw.n_partners; pos += grid_dim) {
    const int4 particle_ids = hbw.partners[pos];
    double t_dx, t_dy, t_dz, t_ax, t_ay, t_az;
    if (poly_psr.gpos_bits > globalpos_scale_nonoverflow_bits) {
      const int95_t i_hdx = int95Subtract(poly_psr.xalt[particle_ids.x],
                                          poly_psr.xalt_ovrf[particle_ids.x],
                                          poly_psr.xalt[particle_ids.y],
                                          poly_psr.xalt_ovrf[particle_ids.y]);
      const int95_t i_hdy = int95Subtract(poly_psr.yalt[particle_ids.x],
                                          poly_psr.yalt_ovrf[particle_ids.x],
                                          poly_psr.yalt[particle_ids.y],
                                          poly_psr.yalt_ovrf[particle_ids.y]);
      const int95_t i_hdz = int95Subtract(poly_psr.zalt[particle_ids.x],
                                          poly_psr.zalt_ovrf[particle_ids.x],
                                          poly_psr.zalt[particle_ids.y],
                                          poly_psr.zalt_ovrf[particle_ids.y]);
      t_dx = splitFPToReal(i_hdx) * poly_psr.inv_gpos_scale;
      t_dy = splitFPToReal(i_hdy) * poly_psr.inv_gpos_scale;
      t_dz = splitFPToReal(i_hdz) * poly_psr.inv_gpos_scale;
      const int95_t i_hax = int95Subtract(poly_psr.xalt[particle_ids.z],
                                          poly_psr.xalt_ovrf[particle_ids.z],
                                          poly_psr.xalt[particle_ids.y],
                                          poly_psr.xalt_ovrf[particle_ids.y]);
      const int95_t i_hay = int95Subtract(poly_psr.yalt[particle_ids.z],
                                          poly_psr.yalt_ovrf[particle_ids.z],
                                          poly_psr.yalt[particle_ids.y],
                                          poly_psr.yalt_ovrf[particle_ids.y]);
      const int95_t i_haz = int95Subtract(poly_psr.zalt[particle_ids.z],
                                          poly_psr.zalt_ovrf[particle_ids.z],
                                          poly_psr.zalt[particle_ids.y],
                                          poly_psr.zalt_ovrf[particle_ids.y]);
      t_ax = splitFPToReal(i_hax) * poly_psr.inv_gpos_scale;
      t_ay = splitFPToReal(i_hay) * poly_psr.inv_gpos_scale;
      t_az = splitFPToReal(i_haz) * poly_psr.inv_gpos_scale;
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

    // Image the coordinates, if necessary
    switch (poly_psr.unit_cell) {
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      {
        const size_t xfrm_start = particle_ids.w * xfrm_stride;
        devcMinimumImage(t_dx, t_dy, t_dz, &poly_psr.umat[xfrm_start], &poly_psr.invu[xfrm_start]);
        devcMinimumImage(t_ax, t_ay, t_az, &poly_psr.umat[xfrm_start], &poly_psr.invu[xfrm_start]);
      }
      break;
    case UnitCellType::NONE:
      break;
    }
    const float local_x[3] = { (float)(t_dx), 0.0f, (float)(t_ax) };
    const float local_y[3] = { (float)(t_dy), 0.0f, (float)(t_ay) };
    const float local_z[3] = { (float)(t_dz), 0.0f, (float)(t_az) };
    const float disp_x = local_x[2] - local_x[0];
    const float disp_y = local_y[2] - local_y[0];
    const float disp_z = local_z[2] - local_z[0];
    const float t_range = sqrtf((disp_x * disp_x) + (disp_y * disp_y) + (disp_z * disp_z));
    const int statblk_idx = (step - hbw.init_step) / hbw.block_steps;
    const int item_idx = pos + (statblk_idx * hbw.n_partners);
    if (t_range < hbw.max_separation) {

      // If the particles are close enough to verify a hydrogen bond, check the angle.
      float t_angle = devcAngle(0, 1, 2, local_x, local_y, local_z);
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
        // accumulators.  Because one thread processes each pair of possible hydrogen bonding
        // partners, asynchronous accumulation is not an issue.
        hbw.tallies[item_idx] += 1;
        hbw.dist_acc[item_idx] += t_range;
        hbw.dist_sqacc[item_idx] += (t_range * t_range);
        hbw.angl_acc[item_idx] += t_angle;
        hbw.angl_sqacc[item_idx] += (t_angle * t_angle);
      }
    }
    hbw.all_dist_acc[item_idx] += t_range;
    hbw.all_dist_sqacc[item_idx] += (t_range * t_range);
  }
}
  
//-------------------------------------------------------------------------------------------------
void launchEvalHydrogenBondAnalysis(const int step, HBondWriter *hbw,
                                    const PsSynthesisReader &poly_psr) {
  kEvalHydrogenBondAnalysis<<<dyna_tk.getSMPCount() * 2, medium_block_size>>>(step, *hbw,
                                                                              poly_psr);
}
  
} // namespace analysis
} // namespace stormm
