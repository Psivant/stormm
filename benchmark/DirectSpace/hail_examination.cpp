#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/math_enumerators.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Math/radial_derivatives.h"
#include "../../src/Math/rounding.h"
#include "../../src/Math/tricubic_cell.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/parsing_enumerators.h"
#include "../../src/Potential/cellgrid.h" 
#include "../../src/Potential/energy_enumerators.h" 
#include "../../src/Potential/layered_potential.h" 
#include "../../src/Potential/layered_potential_metrics.h" 
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/local_arrangement.h"
#include "../../src/Structure/structure_enumerators.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/test_system_manager.h"

#ifndef STORMM_USE_HPC
using stormm::data_types::uint2;
using stormm::data_types::double3;
using stormm::data_types::double4;
#endif
using stormm::data_types::double_type_index;
using stormm::data_types::llint;
using stormm::data_types::ullint;
using stormm::data_types::int95_t;
using stormm::data_types::Ecumenical8;
using namespace stormm::constants;
using namespace stormm::diskutil;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::numerics;
using namespace stormm::parse;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::structure;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Build a synthesis of the periodic systems in a test manager object, with a given perturbation
// and a highly accurate fixed-precision representation.
//
// Arguments:
//   tsm:                   The catalog of test systems
//   xdev:                  Width of the Gaussian perturbation to impart to coordinates
//   igseed:                Seed for the random number generator guiding the perturbations
//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis collectPeriodicCases(const TestSystemManager &tsm, const double xdev,
                                         const int igseed) {
  const std::vector<UnitCellType> uc = { UnitCellType::ORTHORHOMBIC, UnitCellType::TRICLINIC };
  const std::vector<int> quals = tsm.getQualifyingSystems(uc);
  return tsm.exportPhaseSpaceSynthesis(quals, xdev, igseed, 48);
}

//-------------------------------------------------------------------------------------------------
// Create a cell grid for a given coordinate synthesis based on particles with one type of
// non-bonded property.
//
// Arguments:
//   poly_ps:               The coordinate synthesis containing periodic systems
//   usurf:                 The potential surface to decompose with HAIL
//   particle_pair_cutoff:  The cutoff for pairwise particle-particle interactions
//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
CellGrid<T, llint, double, T4> buildSpatialDecomposition(const PhaseSpaceSynthesis &poly_ps,
                                                         const AtomGraphSynthesis &poly_ag,
                                                         const DecomposablePotential usurf,
                                                         const double particle_pair_cutoff) {
  NonbondedTheme cg_theme;
  switch (usurf) {
  case DecomposablePotential::ELECTROSTATIC:
  case DecomposablePotential::ELEC_PME_DIRECT:
    cg_theme = NonbondedTheme::ELECTROSTATIC;
    break;
  case DecomposablePotential::DISPERSION:
    cg_theme = NonbondedTheme::VAN_DER_WAALS;
    break;
  }
  return CellGrid<T, llint, double, T4>(poly_ps, poly_ag, particle_pair_cutoff, 0.25, 4, cg_theme);
}

//-------------------------------------------------------------------------------------------------
// Export all systems, which will form a valid PhaseSpaceSynthesis as they all entail periodic
// boundary conditions, with various perturbations as a means of repeatedly testing the quality of
// the mesh-based potentials at various levels.
//
// Arguments:
//   tsm:                   The catalog of test systems
//   xdev:                  Width of the Gaussian perturbation to impart to coordinates
//   igseed:                Seed for the random number generator guiding the perturbations
//   usurf:                 The potential surface to decompose with HAIL
//   particle_pair_cutoff:  The cutoff for pairwise particle-particle interactions
//   layer_id:              Index of the layer to test
//   perform_check:         Request checks of the veracity of various calculations
//-------------------------------------------------------------------------------------------------
void testHAIL(const TestSystemManager &tsm, const double xdev, const int igseed,
              const DecomposablePotential usurf, const double particle_pair_cutoff,
              const int layer_id, const bool perform_check = false) {
  const PhaseSpaceSynthesis poly_ps = collectPeriodicCases(tsm, xdev, igseed);
  const AtomGraphSynthesis poly_ag(poly_ps.getSystemTopologyPointer());
  CellGrid<double, llint,
           double, double4> cg = buildSpatialDecomposition<double,
                                                           double4>(poly_ps, poly_ag, usurf,
                                                                    particle_pair_cutoff);
  CellGridWriter<double, llint, double, double4> cgw = cg.data();
  LayeredPotentialMetrics lpm(usurf, BoundaryCondition::PERIODIC);
  lpm.setCutoff(particle_pair_cutoff);
  LayeredPotential<double, double4> lnrg(lpm);
  TricubicStencil weights(Interpolant::SMOOTHNESS);

  // Determine critical parameters for the test based on the layer of interest
  const int search_layers = ceil(lpm.getCutoff(layer_id) / cg.getEffectiveCutoff());
  const PsSynthesisReader poly_psr = poly_ps.data();
  const int xfrm_stride = roundUp(9, warp_size_int);
  cg.initializeForces();
  for (int pos = 0; pos < poly_psr.system_count; pos++) {

    // Take all cells within the system and evaluate the exact particle-particle contributions for
    // each level of the decomposition.  Compare this to the mesh-based representation of the same
    // interaction.
    const ullint pcell_dims = cgw.system_cell_grids[pos];
    const int pcell_offset = pcell_dims & 0xfffffffLLU;
    const int pcell_na = ((pcell_dims >> 28) & 0xfffLLU);
    const int pcell_nb = ((pcell_dims >> 40) & 0xfffLLU);
    const int pcell_nc = ((pcell_dims >> 52) & 0xfffLLU);
    
    // Compute the energy of the system and forces on all particles using direct interaction
    // computations.  These provided benchmark energies and benchmark forces to test the resulting
    // interactions.
    const int natom = poly_psr.atom_counts[pos];
    std::vector<double> bench_xfrc(natom, 0.0), bench_yfrc(natom, 0.0), bench_zfrc(natom, 0.0);
    double bench_nrg = 0.0;
    const double* pinvu = &cgw.system_cell_invu[pos * xfrm_stride];
    const AtomGraph *ag_ptr = poly_ps.getSystemTopologyPointer(pos);
    const NonbondedKit<double> nbk = ag_ptr->getDoublePrecisionNonbondedKit();
    for (int i = 0; i < pcell_na; i++) {
      for (int j = 0; j < pcell_nb; j++) {
        for (int k = 0; k < pcell_nc; k++) {

          // Determine limits for searching this cell.  A floating-point representation for the
          // image is critical here, as a fixed-precision representation is calibrated to overflow
          // if atoms from one cell are translated by more than 3-4 cell widths in order to compute
          // pairwise computations with a neighboring cell.
          const size_t cell_ijk_idx = pcell_offset + (((k * pcell_nb) + j) * pcell_na) + i;
          const uint2 cell_ijk_lims = cgw.cell_limits[cell_ijk_idx];
          const uint mhlim = cell_ijk_lims.x + (cell_ijk_lims.y >> 16);
          
          // Determine the limits for searching other cells
          const int minsr_i = i - search_layers;
          const int maxsr_i = i + search_layers;
          const int minsr_j = j - search_layers;
          const int maxsr_j = j + search_layers;
          const int minsr_k = k - search_layers;

          // Loop over other cells
          for (int ngb_k = minsr_k; ngb_k <= k; ngb_k++) {
            const int jlim = (ngb_k == k) ? j : maxsr_j;
            const double dngb_k = ngb_k - k;
            const int reim_ngb_k = ngb_k + (((ngb_k < 0) - (ngb_k >= pcell_nc)) * pcell_nc);
            for (int ngb_j = minsr_j; ngb_j <= jlim; ngb_j++) {
              const int ilim = (ngb_k == k && ngb_j == j) ? i : maxsr_i;
              const double dngb_j = ngb_j -j ;
              const int reim_ngb_j = ngb_j + (((ngb_j < 0) - (ngb_j >= pcell_nb)) * pcell_nb);
              for (int ngb_i = minsr_i; ngb_i <= ilim; ngb_i++) {
                const bool same_cell = (ngb_i == i && ngb_j == j && ngb_k == k);
                const double dngb_i = ngb_i - i;
                const int reim_ngb_i = ngb_i + (((ngb_i < 0) - (ngb_i >= pcell_na)) * pcell_na);

                // Determine limits for searching the neighbor cell.
                const size_t cell_ngb_idx = pcell_offset +
                                            (((reim_ngb_k * pcell_nb) + reim_ngb_j) * pcell_na) +
                                            reim_ngb_i;
                const uint2 cell_ngb_lims = cgw.cell_limits[cell_ngb_idx];
                const uint nhlim = cell_ngb_lims.x + (cell_ngb_lims.y >> 16);

                // Determine the offset for the neighbor cell
                const double nx_shift = (pinvu[0] * dngb_i) + (pinvu[3] * dngb_j) +
                                        (pinvu[6] * dngb_k);
                const double ny_shift = (pinvu[4] * dngb_j) + (pinvu[7] * dngb_k);
                const double nz_shift = (pinvu[8] * dngb_k);
                
                // Loop over all atoms in the neighbor cell with a nested loop in the home cell
                for (uint n = cell_ngb_lims.x; n < nhlim; n++) {
                  const double4 atom_n = cgw.image[n];
                  const double4 atom_n_shft = { atom_n.x + nx_shift, atom_n.y + ny_shift,
                                                atom_n.z + nz_shift, atom_n.w };
                  const int mlim = (same_cell) ? n : mhlim;
                  for (uint m = cell_ijk_lims.x; m < mlim; m++) {
                    const double4 atom_m = cgw.image[m];
                    const double dx = atom_n_shft.x - atom_m.x;
                    const double dy = atom_n_shft.y - atom_m.y;
                    const double dz = atom_n_shft.z - atom_m.z;
                    const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
                    const double r = sqrt(r2);
                    double q_mn;
                    switch (usurf) {
                    case DecomposablePotential::ELECTROSTATIC:
                    case DecomposablePotential::ELEC_PME_DIRECT:
                      q_mn = atom_m.w * atom_n_shft.w;
                      break;
                    case DecomposablePotential::DISPERSION:
                      {
                        const Ecumenical8 mlj_conv = { .d = atom_m.w };
                        const Ecumenical8 nlj_conv = { .d = atom_n.w };
                        const int atom_m_ljidx = mlj_conv.lli;
                        const int atom_n_ljidx = nlj_conv.lli;

                        // The Lennard-Jones tables must be configured for geometric combining
                        // rules for the following assignment to be valid.
                        q_mn = nbk.ljb_coeff[(atom_m_ljidx * nbk.n_lj_types) + atom_n_ljidx];
                      }
                      break;
                    }
                    
                    // Compute the energy contribution
                    const double u_mn = q_mn * lnrg.getAnalyticValue(layer_id, r, r2);
                    bench_nrg += u_mn;
                    
                    // Compute the force contribution
                    const double f_mn = q_mn * poly_psr.frc_scale *
                                        lnrg.getAnalyticDerivative(layer_id, r, r2) / r;
                    const int95_t nfx = hostInt95Sum(cgw.xfrc[n], cgw.xfrc_ovrf[n], f_mn * dx);
                    const int95_t nfy = hostInt95Sum(cgw.yfrc[n], cgw.yfrc_ovrf[n], f_mn * dy);
                    const int95_t nfz = hostInt95Sum(cgw.zfrc[n], cgw.zfrc_ovrf[n], f_mn * dz);
                    const int95_t mfx = hostInt95Sum(cgw.xfrc[m], cgw.xfrc_ovrf[m], -f_mn * dx);
                    const int95_t mfy = hostInt95Sum(cgw.yfrc[m], cgw.yfrc_ovrf[m], -f_mn * dy);
                    const int95_t mfz = hostInt95Sum(cgw.zfrc[m], cgw.zfrc_ovrf[m], -f_mn * dz);
                    cgw.xfrc[n] = nfx.x;
                    cgw.yfrc[n] = nfy.x;
                    cgw.zfrc[n] = nfz.x;
                    cgw.xfrc_ovrf[n] = nfx.y;
                    cgw.yfrc_ovrf[n] = nfy.y;
                    cgw.zfrc_ovrf[n] = nfz.y;
                    cgw.xfrc[m] = mfx.x;
                    cgw.yfrc[m] = mfy.x;
                    cgw.zfrc[m] = mfz.x;
                    cgw.xfrc_ovrf[m] = mfx.y;
                    cgw.yfrc_ovrf[m] = mfy.y;
                    cgw.zfrc_ovrf[m] = mfz.y;
                  }
                }
              }
            }
          }
        }
      }
    }

    // Translate the cell grid forces back to real space
    const int atom_start = poly_psr.atom_starts[pos];
    for (int i = 0; i < natom; i++) {
      const uint cg_idx = cgw.img_atom_idx[i + atom_start];
      if (cg_idx != 0xffffffffU) {
        bench_xfrc[i] = hostInt95ToDouble(cgw.xfrc[cg_idx], cgw.xfrc_ovrf[cg_idx]) *
                        poly_psr.inv_frc_scale;
        bench_yfrc[i] = hostInt95ToDouble(cgw.yfrc[cg_idx], cgw.yfrc_ovrf[cg_idx]) *
                        poly_psr.inv_frc_scale;
        bench_zfrc[i] = hostInt95ToDouble(cgw.zfrc[cg_idx], cgw.zfrc_ovrf[cg_idx]) *
                        poly_psr.inv_frc_scale;
      }
    }
    
    // Make a check of the forces using a complete loop over all particles
    if (perform_check) {
      const CoordinateFrame cf = poly_ps.exportCoordinates(pos);
      const CoordinateFrameReader cfr = cf.data();
      std::vector<double> check_xfrc(natom, 0.0), check_yfrc(natom, 0.0), check_zfrc(natom, 0.0);
      double check_nrg = 0.0;
      for (int i = 0; i < cfr.natom; i++) {
        for (int j = 0; j < i; j++) {
          double dx = cfr.xcrd[j] - cfr.xcrd[i];
          double dy = cfr.ycrd[j] - cfr.ycrd[i];
          double dz = cfr.zcrd[j] - cfr.zcrd[i];
          imageCoordinates<double, double>(&dx, &dy, &dz, cfr.umat, cfr.invu, cfr.unit_cell,
                                           ImagingMethod::MINIMUM_IMAGE);
          const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
          const double r = sqrt(r2);

          // As above, determine the density factor (charge or dispersion) in the interaction
          double q_ij;
          switch (usurf) {
          case DecomposablePotential::ELECTROSTATIC:
          case DecomposablePotential::ELEC_PME_DIRECT:
            q_ij = nbk.charge[i] * nbk.charge[j];
            break;
          case DecomposablePotential::DISPERSION:

            // The Lennard-Jones tables must be configured for geometric combining
            // rules for the following assignment to be valid.
            const int ilj_idx = nbk.lj_idx[i];
            const int jlj_idx = nbk.lj_idx[j];
            q_ij = nbk.ljb_coeff[(nbk.lj_idx[i] * nbk.n_lj_types) + nbk.lj_idx[j]];
            break;
          }         

          // Compute the energy contribution
          const double u_ij = q_ij * lnrg.getAnalyticValue(layer_id, r, r2);
          check_nrg += u_ij;

          // Compute the force contribution
          const double f_ij = q_ij * lnrg.getAnalyticDerivative(layer_id, r, r2) / r;
          check_xfrc[j] += f_ij * dx;
          check_yfrc[j] += f_ij * dy;
          check_zfrc[j] += f_ij * dz;
          check_xfrc[i] -= f_ij * dx;
          check_yfrc[i] -= f_ij * dy;
          check_zfrc[i] -= f_ij * dz;
        }
      }

      // Check the forces
      const double xfrc_mue = meanUnsignedError(bench_xfrc, check_xfrc);
      const double yfrc_mue = meanUnsignedError(bench_yfrc, check_yfrc);
      const double zfrc_mue = meanUnsignedError(bench_zfrc, check_zfrc);
      if (xfrc_mue > 1.0e-5 || yfrc_mue > 1.0e-5 || zfrc_mue > 1.0e-5) {
        rtErr("Mean unsigned errors indicate that benchmark forces are in disagreement with "
              "forces computed using a simpler, all-to-all loop: X " +
              realToString(xfrc_mue, 9, 2, NumberFormat::SCIENTIFIC) + ", Y " +
              realToString(yfrc_mue, 9, 2, NumberFormat::SCIENTIFIC) + ", Z " +
              realToString(zfrc_mue, 9, 2, NumberFormat::SCIENTIFIC) + ".", "testHAIL");
      }
    }
    
    // Compute the potential and various derivatives at the vertices of the cell grid
    std::vector<double> u_self(8), dudx_self(8), dudy_self(8), dudz_self(8), dudxx_self(8);
    std::vector<double> dudxy_self(8), dudxz_self(8), dudyy_self(8), dudyz_self(8), dudxxx_self(8);
    std::vector<double> dudxxy_self(8), dudxxz_self(8), dudxyy_self(8), dudxyz_self(8);
    const size_t total_cells = pcell_na * pcell_nb * pcell_nc;
    std::vector<double> u_mesh(total_cells, 0.0), dudx_mesh(total_cells, 0.0);
    std::vector<double> dudy_mesh(total_cells, 0.0), dudz_mesh(total_cells, 0.0);
    std::vector<double> dudxx_mesh(total_cells, 0.0), dudxy_mesh(total_cells, 0.0);
    std::vector<double> dudxz_mesh(total_cells, 0.0), dudyy_mesh(total_cells, 0.0);
    std::vector<double> dudyz_mesh(total_cells, 0.0), dudxxx_mesh(total_cells, 0.0);
    std::vector<double> dudxxy_mesh(total_cells, 0.0), dudxxz_mesh(total_cells, 0.0);
    std::vector<double> dudxyy_mesh(total_cells, 0.0), dudxyz_mesh(total_cells, 0.0);
    std::vector<double> intrp_xfrc(natom, 0.0), intrp_yfrc(natom, 0.0), intrp_zfrc(natom, 0.0);
    double intrp_nrg = 0.0;
    const std::vector<double> interpolation_cell = {      0.0,      0.0,      0.0,
                                                     pinvu[0], pinvu[1], pinvu[2],
                                                     pinvu[3], pinvu[4], pinvu[5],
                                                     pinvu[6], pinvu[7], pinvu[8] };
    for (int i = 0; i < pcell_na; i++) {
      for (int j = 0; j < pcell_nb; j++) {
        for (int k = 0; k < pcell_nc; k++) {
          const size_t cell_ijk_idx = pcell_offset + (((k * pcell_nb) + j) * pcell_na) + i;
          const uint2 cell_ijk_lims = cgw.cell_limits[cell_ijk_idx];
          const uint mhlim = cell_ijk_lims.x + (cell_ijk_lims.y >> 16);

          // Subtract off self interactions of particles on the mesh.  While the quantities
          // computed are legitimate interactions and could be contributed to the cell grid mesh
          // itself, skip this optimization in favor of clarity.
          for (uint m = cell_ijk_lims.x; m < mhlim; m++) {
            const double4 atom_m = cgw.image[m];
            const int topol_idx = cgw.nonimg_atom_idx[m] - poly_psr.atom_starts[pos];
            double q_m;
            switch (usurf) {
            case DecomposablePotential::ELECTROSTATIC:
            case DecomposablePotential::ELEC_PME_DIRECT:
              q_m = nbk.charge[topol_idx];
              break;
            case DecomposablePotential::DISPERSION:
              q_m = sqrt(nbk.ljb_coeff[(nbk.n_lj_types + 1) * nbk.lj_idx[topol_idx]]);
              break;
            }
            for (int ngb_i = 0; ngb_i <= 1; ngb_i++) {
              const double dngb_i = ngb_i;
              for (int ngb_j = 0; ngb_j <= 1; ngb_j++) {
                const double dngb_j = ngb_j;
                for (int ngb_k = 0; ngb_k <= 1; ngb_k++) {
                  const double dngb_k = ngb_k;

                  // Determine the relative location of the cell corner
                  const double cnx = (pinvu[0] * dngb_i) + (pinvu[3] * dngb_j) +
                                     (pinvu[6] * dngb_k);
                  const double cny = (pinvu[4] * dngb_j) + (pinvu[7] * dngb_k);
                  const double cnz = (pinvu[8] * dngb_k);
                  const double dx = cnx - atom_m.x;
                  const double dy = cny - atom_m.y;
                  const double dz = cnz - atom_m.z;
                  const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
                  const double r  = sqrt(r2);
                  
                  // Compute the potential and derivatives
                  const int cn_idx = (((ngb_k * 2) + ngb_j) * 2) + ngb_i;
                  const double u   = q_m * lnrg.getAnalyticValue(layer_id, r, r2);
                  const double du  = q_m * lnrg.getAnalyticDerivative(layer_id, r, r2, 1);
                  const double d2u = q_m * lnrg.getAnalyticDerivative(layer_id, r, r2, 2);
                  const double d3u = q_m * lnrg.getAnalyticDerivative(layer_id, r, r2, 3);
                  u_self[cn_idx] = u;
                  dudx_self[cn_idx]   = du * dx / r;
                  dudy_self[cn_idx]   = du * dy / r;
                  dudz_self[cn_idx]   = du * dz / r;
                  dudxy_self[cn_idx]  = radialSecondDerivative<double>(du, d2u, dx, dy, r, r2);
                  dudxz_self[cn_idx]  = radialSecondDerivative<double>(du, d2u, dx, dz, r, r2);
                  dudyz_self[cn_idx]  = radialSecondDerivative<double>(du, d2u, dy, dz, r, r2);
                  dudxyz_self[cn_idx] = radialThirdDerivative<double>(du, d2u, d3u, dx, dy, dz, r,
                                                                      r2);
                  switch (poly_psr.unit_cell) {
                  case UnitCellType::ORTHORHOMBIC:
                    dudxx_self[cn_idx] = 0.0;
                    dudyy_self[cn_idx] = 0.0;
                    dudxxx_self[cn_idx] = 0.0;
                    dudxxy_self[cn_idx] = 0.0;
                    dudxxz_self[cn_idx] = 0.0;
                    dudxyy_self[cn_idx] = 0.0;
                    break;
                  case UnitCellType::TRICLINIC:
                    dudxx_self[cn_idx]  = radialSecondDerivative<double>(du, d2u, dx, dx, r, r2);
                    dudyy_self[cn_idx]  = radialSecondDerivative<double>(du, d2u, dy, dy, r, r2);
                    dudxxx_self[cn_idx] = radialThirdDerivative<double>(du, d2u, d3u, dx, r, r2);
                    dudxxy_self[cn_idx] = radialThirdDerivative<double>(du, d2u, d3u, dx, dy, r,
                                                                        r2);
                    dudxxz_self[cn_idx] = radialThirdDerivative<double>(du, d2u, d3u, dx, dz, r,
                                                                        r2);
                    dudxyy_self[cn_idx] = radialThirdDerivative<double>(du, d2u, d3u, dy, dx, r,
                                                                        r2);
                    break;
                  case UnitCellType::NONE:
                    rtErr("These benchmarks are designed for periodic simulation systems.",
                          "testHAIL");
                    break;
                  }
                }
              }
            }

            // Form a tricubic cell to carry out the interpolation
            const TricubicCell<double> s_tric(weights, interpolation_cell, u_self, dudx_self,
                                              dudy_self, dudz_self, dudxx_self, dudxy_self,
                                              dudxz_self, dudyy_self, dudyz_self, dudxxx_self,
                                              dudxxy_self, dudxxz_self, dudxyy_self, dudxyz_self);

            // The energy can be immediately multipled by the density factor, but each member of
            // the force tuple that must be multiplied separately.
            const double self_nrg = q_m * s_tric.evaluate(atom_m.x, atom_m.y, atom_m.z);
            const double3 self_frc = s_tric.derivative<double3>(atom_m.x, atom_m.y, atom_m.z);
            intrp_nrg -= self_nrg;
            intrp_xfrc[topol_idx] -= q_m * self_frc.x;
            intrp_yfrc[topol_idx] -= q_m * self_frc.y;
            intrp_zfrc[topol_idx] -= q_m * self_frc.z;
          }
          
          // Calculate for points in the applicable region
          for (int ngb_i = i - search_layers; ngb_i <= i + 1 + search_layers; ngb_i++) {
            const double dngb_i = ngb_i - i;
            const int reim_ngb_i = ngb_i + (((ngb_i < 0) - (ngb_i >= pcell_na)) * pcell_na);
            for (int ngb_j = j - search_layers; ngb_j <= j + 1 + search_layers; ngb_j++) {
              const double dngb_j = ngb_j - j;
              const int reim_ngb_j = ngb_j + (((ngb_j < 0) - (ngb_j >= pcell_nb)) * pcell_nb);
              for (int ngb_k = k - search_layers; ngb_k <= k + 1 + search_layers; ngb_k++) {
                const double dngb_k = ngb_k - k;
                const int reim_ngb_k = ngb_k + (((ngb_k < 0) - (ngb_k >= pcell_nc)) * pcell_nc);
                const size_t mesh_pt_idx = (((reim_ngb_k * pcell_nb) + reim_ngb_j) * pcell_na) +
                                           reim_ngb_i;
                
                // Determine the relative location of the cell grid vertex
                const double mesh_nx = (pinvu[0] * dngb_i) + (pinvu[3] * dngb_j) +
                                       (pinvu[6] * dngb_k);
                const double mesh_ny = (pinvu[4] * dngb_j) + (pinvu[7] * dngb_k);
                const double mesh_nz = (pinvu[8] * dngb_k);
                double u_pt = 0.0;
                double dudx_pt = 0.0;
                double dudy_pt = 0.0;
                double dudz_pt = 0.0;
                double dudxx_pt = 0.0;
                double dudxy_pt = 0.0;
                double dudxz_pt = 0.0;
                double dudyy_pt = 0.0;
                double dudyz_pt = 0.0;
                double dudxxx_pt = 0.0;
                double dudxxy_pt = 0.0;
                double dudxxz_pt = 0.0;
                double dudxyy_pt = 0.0;
                double dudxyz_pt = 0.0;
                for (int m = cell_ijk_lims.x; m < mhlim; m++) {
                  const double4 atom_m = cgw.image[m];
                  const int topol_idx = cgw.nonimg_atom_idx[m] - poly_psr.atom_starts[pos];
                  const double dx = mesh_nx - atom_m.x;
                  const double dy = mesh_ny - atom_m.y;
                  const double dz = mesh_nz - atom_m.z;
                  const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
                  const double r  = sqrt(r2);
                  double q_m;
                  switch (usurf) {
                  case DecomposablePotential::ELECTROSTATIC:
                  case DecomposablePotential::ELEC_PME_DIRECT:
                    q_m = nbk.charge[topol_idx];
                    break;
                  case DecomposablePotential::DISPERSION:
                    q_m = sqrt(nbk.ljb_coeff[(nbk.n_lj_types + 1) * nbk.lj_idx[topol_idx]]);
                    break;
                  }
                  const double u   = q_m * lnrg.getAnalyticValue(layer_id, r, r2);
                  const double du  = q_m * lnrg.getAnalyticDerivative(layer_id, r, r2, 1);
                  const double d2u = q_m * lnrg.getAnalyticDerivative(layer_id, r, r2, 2);
                  const double d3u = q_m * lnrg.getAnalyticDerivative(layer_id, r, r2, 3);
                  u_pt += u;
                  dudx_pt   += du * dx / r;
                  dudy_pt   += du * dy / r;
                  dudz_pt   += du * dz / r;
                  dudxy_pt  += radialSecondDerivative<double>(du, d2u, dx, dy, r, r2);
                  dudxz_pt  += radialSecondDerivative<double>(du, d2u, dx, dz, r, r2);
                  dudyz_pt  += radialSecondDerivative<double>(du, d2u, dy, dz, r, r2);
                  dudxyz_pt += radialThirdDerivative<double>(du, d2u, d3u, dx, dy, dz, r, r2);
                  switch (poly_psr.unit_cell) {
                  case UnitCellType::ORTHORHOMBIC:
                    break;
                  case UnitCellType::TRICLINIC:
                    dudxx_pt  += radialSecondDerivative<double>(du, d2u, dx, dx, r, r2);
                    dudyy_pt  += radialSecondDerivative<double>(du, d2u, dy, dy, r, r2);
                    dudxxx_pt += radialThirdDerivative<double>(du, d2u, d3u, dx, r, r2);
                    dudxxy_pt += radialThirdDerivative<double>(du, d2u, d3u, dx, dy, r, r2);
                    dudxxz_pt += radialThirdDerivative<double>(du, d2u, d3u, dx, dz, r, r2);
                    dudxyy_pt += radialThirdDerivative<double>(du, d2u, d3u, dy, dx, r, r2);
                    break;
                  case UnitCellType::NONE:
                    break;
                  }
                }
                u_mesh[mesh_pt_idx]      += u_pt;
                dudx_mesh[mesh_pt_idx]   += dudx_pt;
                dudy_mesh[mesh_pt_idx]   += dudy_pt;
                dudz_mesh[mesh_pt_idx]   += dudz_pt;
                dudxy_mesh[mesh_pt_idx]  += dudxy_pt;
                dudxz_mesh[mesh_pt_idx]  += dudxz_pt;
                dudyz_mesh[mesh_pt_idx]  += dudyz_pt;
                dudxyz_mesh[mesh_pt_idx] += dudxyz_pt;
                switch (poly_psr.unit_cell) {
                case UnitCellType::ORTHORHOMBIC:
                  break;
                case UnitCellType::TRICLINIC:
                  dudxx_mesh[mesh_pt_idx]  += dudxx_pt;
                  dudyy_mesh[mesh_pt_idx]  += dudyy_pt;
                  dudxxx_mesh[mesh_pt_idx]  += dudxxx_pt;
                  dudxxy_mesh[mesh_pt_idx]  += dudxxy_pt;
                  dudxxz_mesh[mesh_pt_idx]  += dudxxz_pt;
                  dudxyy_mesh[mesh_pt_idx]  += dudxyy_pt;
                  break;
                case UnitCellType::NONE:
                  break;
                }
              }
            }
          }
        }
      }
    }

    // With the mesh constructed, formulate TricubicCell objects for each neighbor list cell.
    // Compute the energy contributions and forces for all atoms within.
    for (int i = 0; i < pcell_na; i++) {
      for (int j = 0; j < pcell_nb; j++) {
        for (int k = 0; k < pcell_nc; k++) {

          // Re-use the arrays allocaed for computing the self interactions to define the corners
          // of this mesh cell.
          for (int ngb_i = i; ngb_i <= i + 1; ngb_i++) {
            const int reim_ngb_i = ngb_i + (((ngb_i < 0) - (ngb_i >= pcell_na)) * pcell_na);
            for (int ngb_j = j; ngb_j <= j + 1; ngb_j++) {
              const int reim_ngb_j = ngb_j + (((ngb_j < 0) - (ngb_j >= pcell_nb)) * pcell_nb);
              for (int ngb_k = k; ngb_k <= k + 1; ngb_k++) {
                const int reim_ngb_k = ngb_k + (((ngb_k < 0) - (ngb_k >= pcell_nc)) * pcell_nc);
                const size_t cn_idx = ((((ngb_k - k) * 2) + ngb_j - j) * 2) + ngb_i - i;
                const size_t mesh_idx = (((reim_ngb_k * pcell_nb) + reim_ngb_j) * pcell_na) +
                                        reim_ngb_i;
                u_self[cn_idx]      = u_mesh[mesh_idx];
                dudx_self[cn_idx]   = dudx_mesh[mesh_idx];
                dudy_self[cn_idx]   = dudy_mesh[mesh_idx];
                dudz_self[cn_idx]   = dudz_mesh[mesh_idx];
                dudxy_self[cn_idx]  = dudxy_mesh[mesh_idx];
                dudxz_self[cn_idx]  = dudxz_mesh[mesh_idx];
                dudyz_self[cn_idx]  = dudyz_mesh[mesh_idx];
                dudxyz_self[cn_idx] = dudxyz_mesh[mesh_idx];
                switch (poly_psr.unit_cell) {
                case UnitCellType::ORTHORHOMBIC:
                  break;
                case UnitCellType::TRICLINIC:
                  dudxx_self[cn_idx]  = dudxx_mesh[mesh_idx];
                  dudyy_self[cn_idx]  = dudyy_mesh[mesh_idx];
                  dudxxx_self[cn_idx] = dudxxx_mesh[mesh_idx];
                  dudxxy_self[cn_idx] = dudxxy_mesh[mesh_idx];
                  dudxxz_self[cn_idx] = dudxxz_mesh[mesh_idx];
                  dudxyy_self[cn_idx] = dudxyy_mesh[mesh_idx];
                  break;
                case UnitCellType::NONE:
                  break;
                }
              }
            }
          }

          // Create the interpolation framework.  Apply it to all atoms in the neighbor list cell.
          const TricubicCell<double> s_tric(weights, interpolation_cell, u_self, dudx_self,
                                            dudy_self, dudz_self, dudxx_self, dudxy_self,
                                            dudxz_self, dudyy_self, dudyz_self, dudxxx_self,
                                            dudxxy_self, dudxxz_self, dudxyy_self, dudxyz_self);
          const size_t cell_ijk_idx = pcell_offset + (((k * pcell_nb) + j) * pcell_na) + i;
          const uint2 cell_ijk_lims = cgw.cell_limits[cell_ijk_idx];
          const uint mhlim = cell_ijk_lims.x + (cell_ijk_lims.y >> 16);
          for (uint m = cell_ijk_lims.x; m < mhlim; m++) {
            const double4 atom_m = cgw.image[m];
            const int topol_idx = cgw.nonimg_atom_idx[m] - poly_psr.atom_starts[pos];
            double q_m;
            switch (usurf) {
            case DecomposablePotential::ELECTROSTATIC:
            case DecomposablePotential::ELEC_PME_DIRECT:
              q_m = nbk.charge[topol_idx];
              break;
            case DecomposablePotential::DISPERSION:
              q_m = sqrt(nbk.ljb_coeff[(nbk.n_lj_types + 1) * nbk.lj_idx[topol_idx]]);
              break;
            }

            // The energy can be immediately multipled by the density factor, but each member of
            // the force tuple that must be multiplied separately.
            const double test_nrg = q_m * s_tric.evaluate(atom_m.x, atom_m.y, atom_m.z);
            const double3 test_frc = s_tric.derivative<double3>(atom_m.x, atom_m.y, atom_m.z);
            intrp_nrg += test_nrg;
            intrp_xfrc[topol_idx] += q_m * test_frc.x;
            intrp_yfrc[topol_idx] += q_m * test_frc.y;
            intrp_zfrc[topol_idx] += q_m * test_frc.z;
          }
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
// A plain struct to hold the results of a particle-point computation.  The default constructor and
// various copy and move constructors / assignment operators are all taken.  Various members track
// all components of the interaction that might be necessary for assembling a tricubic interpolant,
// or a straightforward particle-particle interaction.
//-------------------------------------------------------------------------------------------------
struct ParticlePointInfluence {
  double u;       ///< Function value (energy)
  double dudx;    ///< Partial derivative (force) in x
  double dudy;    ///< Partial derivative in y
  double dudz;    ///< Partial derivative in z
  double dudxx;   ///< Double partial derivative in x
  double dudxy;   ///< Mixed partial derivative in x, y
  double dudxz;   ///< Mixed partial derivative in x, z
  double dudyy;   ///< Double partial derivative in y
  double dudyz;   ///< Mixed partial derivative in y, z
  double dudxxx;  ///< Triple partial derivative in x
  double dudxxy;  ///< Mixed partial derivative in x2, y
  double dudxxz;  ///< Mixed partial derivative in x2, z
  double dudxyy;  ///< Mixed partial derivative in x, y2
  double dudxyz;  ///< Mixed partial derivative in x, y, and z
};

//-------------------------------------------------------------------------------------------------
// Compute the influence of all particles in some outer shell of cell grid elements with a point,
// which will be somewhere in the inner cell.  The point in question could be the location of one
// of the atoms in the inner cell, or one of the corners of the inner cell.
//
// Arguments:
//   cgw:         Abstract of the cell grid, containing a neighbor list spatial decomposition of
//                each system in the coordinate synthesis
//   poly_ps:     The coordinate synthesis, containing pointers to topologies and offsets needed
//                to interpret some cell grid indexing
//   system_idx:  Index of the system of interest within the cell grid and the coordinate synthesis
//   icell_a:     Index of the inner cell within the cell grid along the A axis
//   icell_b:     Index of the inner cell within the cell grid along the B axis
//   icell_c:     Index of the inner cell within the cell grid along the C axis
//   spacing:     Number of spatial decomposition cells between the inner cell and the outer shell
//   ptx:         Cartesian x location of the point, relative to the origin of the inner cell
//   pty:         Cartesian y location of the point, relative to the origin of the inner cell
//   ptz:         Cartesian z location of the point, relative to the origin of the inner cell
//   uform:       Form of the potential to apply
//   ew_coeff:    Value of the Ewald coefficient, if applying a PME-style potential
//   eval_all:    Flag to indicate that mixed partial derivatives necessary to compute an
//                interpolant should all be evaluated.  If set to FALSE, only the enery and force
//                will be evaluated.
//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
ParticlePointInfluence outerShellEffect(const CellGridWriter<T, llint, double, T4> &cgw,
                                        const PhaseSpaceSynthesis &poly_ps, const int system_idx,
                                        const int icell_a, const int icell_b, const int icell_c,
                                        const int spacing, const T ptx, const T pty, const T ptz,
                                        const DecomposablePotential uform, const double ew_coeff,
                                        const bool eval_all) {

  // Loop over all cells in the outer shell.
  const ullint pcell_dims = cgw.system_cell_grids[system_idx];
  const int pcell_offset = pcell_dims & 0xfffffffLLU;
  const int pcell_na = ((pcell_dims >> 28) & 0xfffLLU);
  const int pcell_nb = ((pcell_dims >> 40) & 0xfffLLU);
  const int pcell_nc = ((pcell_dims >> 52) & 0xfffLLU);
  const int ci_llim = icell_a - spacing - 1;
  const int ci_hlim = icell_a + spacing + 1;
  const int cj_llim = icell_b - spacing - 1;
  const int cj_hlim = icell_b + spacing + 1;
  const int ck_llim = icell_c - spacing - 1;
  const int ck_hlim = icell_c + spacing + 1;
  const AtomGraph *ag_ptr = poly_ps.getSystemTopologyPointer(system_idx);
  const NonbondedKit<double> nbk = ag_ptr->getDoublePrecisionNonbondedKit();
  const int synth_atom_offset = poly_ps.getAtomOffset(system_idx);
  const bool tcalc_is_double = (std::type_index(typeid(T)).hash_code() == double_type_index);
  ParticlePointInfluence result = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  const T kcoul = stormm::symbols::amber_ancient_bioq;
  const int xfrm_stride = roundUp(9, warp_size_int);
  const T* invu = &cgw.system_cell_invu[system_idx * xfrm_stride];
  for (int i = ci_llim; i <= ci_hlim; i++) {
    const T tngb_i = i;
    const int reim_i = i + (((i < 0) - (i >= pcell_na)) * pcell_na);
    for (int j = cj_llim; j <= cj_hlim; j++) {
      const T tngb_j = j;
      const int reim_j = j + (((j < 0) - (j >= pcell_nb)) * pcell_nb);
      for (int k = ck_llim; k <= ck_hlim; k++) {
        const T tngb_k = k;

        // Skip this cell if it is not part of the shell
        if (i != ci_llim && i != ci_hlim && j != cj_llim && j != cj_hlim &&
            k != ck_llim && k != ck_hlim) {
          continue;
        }

        // Compute the coordinate frame offset relative to the inner cell
        const T shift_ijk_x = (invu[0] * tngb_i) + (invu[3] * tngb_j) + (invu[6] * tngb_k);
        const T shift_ijk_y =                      (invu[4] * tngb_j) + (invu[7] * tngb_k);
        const T shift_ijk_z =                                           (invu[8] * tngb_k);

        // Loop over all atoms of the present neighbor cell
        const int reim_k = k + (((k < 0) - (k >= pcell_nc)) * pcell_nc);
        const int cell_ijk_idx = pcell_offset + (((reim_k * pcell_nb) + reim_j) * pcell_na) +
                                 reim_i;
        const uint2 cell_ijk_lims = cgw.cell_limits[cell_ijk_idx];
        const uint mlim = cell_ijk_lims.x + (cell_ijk_lims.y >> 16);
        for (int m = cell_ijk_lims.x; m < mlim; m++) {
          const T4 atom_m = cgw.image[m];
          const int topol_idx_m = cgw.nonimg_atom_idx[m] - synth_atom_offset;
          const T dx = ptx - (atom_m.x + shift_ijk_x);
          const T dy = pty - (atom_m.y + shift_ijk_y);
          const T dz = ptz - (atom_m.z + shift_ijk_z);
          const T r2 = (dx * dx) + (dy * dy) + (dz * dz);
          const T r = (tcalc_is_double) ? sqrt(r2) : sqrtf(r2);
          switch (uform) {
          case DecomposablePotential::ELECTROSTATIC:
            break;
          case DecomposablePotential::ELEC_PME_DIRECT:
            {
              const T q_m = nbk.charge[topol_idx_m];

              // Always evaluate energy and force
              const T u   = q_m * elecPMEDirectSpace<T>(ew_coeff, kcoul, r, r2, 0);
              const T du  = q_m * elecPMEDirectSpace<T>(ew_coeff, kcoul, r, r2, 1);
              result.u += u;
              result.dudx += radialFirstDerivative<T>(du, dx, r);
              result.dudy += radialFirstDerivative<T>(du, dy, r);
              result.dudz += radialFirstDerivative<T>(du, dz, r);
              if (eval_all) {
                const T d2u = q_m * elecPMEDirectSpace<T>(ew_coeff, kcoul, r, r2, 2);
                const T d3u = q_m * elecPMEDirectSpace<T>(ew_coeff, kcoul, r, r2, 3);
                result.dudxx  += radialSecondDerivative<T>(du, d2u, dx, r);
                result.dudxy  += radialSecondDerivative<T>(du, d2u, dx, dy, r, r2);
                result.dudxz  += radialSecondDerivative<T>(du, d2u, dx, dz, r, r2);
                result.dudyy  += radialSecondDerivative<T>(du, d2u, dy, r);
                result.dudyz  += radialSecondDerivative<T>(du, d2u, dy, dz, r, r2);
                result.dudxxx += radialThirdDerivative<T>(du, d2u, d3u, dx, r, r2);
                result.dudxxy += radialThirdDerivative<T>(du, d2u, d3u, dx, dy, r, r2);
                result.dudxxz += radialThirdDerivative<T>(du, d2u, d3u, dx, dz, r, r2);
                result.dudxyy += radialThirdDerivative<T>(du, d2u, d3u, dy, dx, r, r2);
                result.dudxyz += radialThirdDerivative<T>(du, d2u, d3u, dx, dy, dz, r, r2);
              }
            }
            break;
          case DecomposablePotential::DISPERSION:
            {
              const double ljb = nbk.ljb_coeff[(nbk.n_lj_types + 1) * nbk.lj_idx[topol_idx_m]];
              const T q_m = (tcalc_is_double) ? sqrt(ljb) : sqrtf(ljb);
            }
            break;
          }
        }
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Compute the particle-particle interactions and compare to an interpolated approximation for a
// given potential form and cell separation on the spatial decomposition lattice.
//
// Arguments:
//   tsm:                   The catalog of test systems
//   xdev:                  Width of the Gaussian perturbation to impart to coordinates
//   igseed:                Seed for the random number generator guiding the perturbations
//   uform:                 The potential to evaluate
//   particle_pair_cutoff:  The cutoff for pairwise particle-particle interactions
//   ew_coeff:              The Ewald coefficient, for applying PME-based potential forms
//   standoff:              The number of cell grid layers to skip before computing interactions
//-------------------------------------------------------------------------------------------------
template <typename T, typename T3, typename T4>
void testStandoff(const TestSystemManager &tsm, const double xdev, const int igseed,
                  const DecomposablePotential uform, const double particle_pair_cutoff,
                  const T ew_coeff, const int standoff) {
  const PhaseSpaceSynthesis poly_ps = collectPeriodicCases(tsm, xdev, igseed);
  const AtomGraphSynthesis poly_ag(poly_ps.getSystemTopologyPointer());
  const PsSynthesisReader poly_psr = poly_ps.data();
  CellGrid<T, llint, double, T4> cg = buildSpatialDecomposition<T, T4>(poly_ps, poly_ag, uform,
                                                                       particle_pair_cutoff);
  CellGridWriter<T, llint, double, T4> cgw = cg.data();
  TricubicStencil weights(Interpolant::SMOOTHNESS);

  // Determine the calculation type based on the base data type of the cell grid
  const bool tcalc_is_double = (std::type_index(typeid(T)).hash_code() == double_type_index);
  
  // Loop over all cells, taking each as the "inner" cell.  First compute the forces that particles
  // in the outer shell would exert by explicit particle-particle interactions.  Next, compute the
  // influence of particles in the outer shell at the corners of the inner cell and interpolate the
  // forces from this influence to approximate the interactions of particles in the inner cell with
  // those in the outer shell.
  for (int pos = 0; pos < poly_psr.system_count; pos++) {
    const ullint pcell_dims = cgw.system_cell_grids[pos];
    const int pcell_offset = pcell_dims & 0xfffffffLLU;
    const int pcell_na = ((pcell_dims >> 28) & 0xfffLLU);
    const int pcell_nb = ((pcell_dims >> 40) & 0xfffLLU);
    const int pcell_nc = ((pcell_dims >> 52) & 0xfffLLU);
    const int xfrm_stride = roundUp(9, warp_size_int);
    const T* invu = &cgw.system_cell_invu[pos * xfrm_stride];
    const std::vector<double> interpolation_cell = {     0.0,     0.0,     0.0,
                                                     invu[0], invu[1], invu[2],
                                                     invu[3], invu[4], invu[5],
                                                     invu[6], invu[7], invu[8] };
    const AtomGraph *ag_ptr = poly_ps.getSystemTopologyPointer(pos);
    const NonbondedKit<double> nbk = ag_ptr->getDoublePrecisionNonbondedKit();
    for (int i = 0; i < pcell_na; i++) {
      for (int j = 0; j < pcell_nb; j++) {
        for (int k = 0; k < pcell_nc; k++) {

          // Compute the explicit force of all particles in the outer shell on a particle with
          // locations at the cell corners.
          std::vector<double> u(8), dudx(8), dudy(8), dudz(8), dudxx(8), dudxy(8), dudxz(8);
          std::vector<double> dudyy(8), dudyz(8), dudxxx(8), dudxxy(8), dudxxz(8), dudxyy(8);
          std::vector<double> dudxyz(8);
          for (int ci = 0; ci < 2; ci++) {
            const T tci = ci;
            for (int cj = 0; cj < 2; cj++) {
              const T tcj = cj;
              for (int ck = 0; ck < 2; ck++) {
                const T tck = ck;
                const T ptx = (invu[0] * tci) + (invu[3] * tcj) + (invu[6] * tck);
                const T pty =                   (invu[4] * tcj) + (invu[7] * tck);
                const T ptz =                                     (invu[8] * tck);
                const ParticlePointInfluence ppi = outerShellEffect(cgw, poly_ps, pos, i, j, k,
                                                                    standoff, ptx, pty, ptz, uform,
                                                                    ew_coeff, true);
                const size_t cn_idx = (((ck * 2) + cj) * 2) + ci;
                u[cn_idx] = ppi.u;
                dudx[cn_idx]   = ppi.dudx;
                dudy[cn_idx]   = ppi.dudy;
                dudz[cn_idx]   = ppi.dudz;
                dudxx[cn_idx]  = ppi.dudxx;
                dudxy[cn_idx]  = ppi.dudxy;
                dudxz[cn_idx]  = ppi.dudxz;
                dudyy[cn_idx]  = ppi.dudyy;
                dudyz[cn_idx]  = ppi.dudyz;
                dudxxx[cn_idx] = ppi.dudxxx;
                dudxxy[cn_idx] = ppi.dudxxy;
                dudxxz[cn_idx] = ppi.dudxxz;
                dudxyy[cn_idx] = ppi.dudxyy;
                dudxyz[cn_idx] = ppi.dudxyz;
              }
            }
          }

          // Create the tricubic interpolant, then approximate the force and energy for each
          // particle in the cell.  Compute the explicit interaction energy for comparison.
          const TricubicCell<T> s_tric(weights, interpolation_cell, u, dudx, dudy, dudz, dudxx,
                                       dudxy, dudxz, dudyy, dudyz, dudxxx, dudxxy, dudxxz, dudxyy,
                                       dudxyz);
          const int cell_ijk_idx = pcell_offset + (((k * pcell_nb) + j) * pcell_na) + i;
          const uint2 cell_ijk_lims = cgw.cell_limits[cell_ijk_idx];
          const uint mlim = cell_ijk_lims.x + (cell_ijk_lims.y >> 16);
          for (int m = cell_ijk_lims.x; m < mlim; m++) {
            const T4 atom_m = cgw.image[m];
            const int topol_idx_m = cgw.nonimg_atom_idx[m] - poly_psr.atom_starts[pos];
            T q_m;
            switch (uform) {
            case DecomposablePotential::ELECTROSTATIC:
            case DecomposablePotential::ELEC_PME_DIRECT:
              q_m = nbk.charge[topol_idx_m];
              break;
            case DecomposablePotential::DISPERSION:
              {
                const double ljb = nbk.ljb_coeff[(nbk.n_lj_types + 1) * nbk.lj_idx[topol_idx_m]];
                q_m = (tcalc_is_double) ? sqrt(ljb) : sqrtf(ljb);
              }
              break;
            }
            const ParticlePointInfluence ppi = outerShellEffect(cgw, poly_ps, pos, i, j, k,
                                                                standoff, atom_m.x, atom_m.y,
                                                                atom_m.z, uform, ew_coeff, false);
            const T approx_u = s_tric.evaluate(atom_m.x, atom_m.y, atom_m.z);
            const T3 approx_f = s_tric.template derivative<T3>(atom_m.x, atom_m.y, atom_m.z);
          }

          // CHECK
          if (i == 0 && j == 0 && k == 0 && pos == 0) {
            std::vector<double> act_u(32768), act_fx(32768), act_fy(32768), act_fz(32768);
            std::vector<double> apr_u(32768), apr_fx(32768), apr_fy(32768), apr_fz(32768);
            for (int ci = 0; ci < 32; ci++) {
              const T tci = 0.015625 + (0.03125 * static_cast<T>(ci));
              for (int cj = 0; cj < 32; cj++) {
                const T tcj = 0.015625 + (0.03125 * static_cast<T>(cj));
                for (int ck = 0; ck < 32; ck++) {
                  const T tck = 0.015625 + (0.03125 * static_cast<T>(ck));
                  const T ptx = (invu[0] * tci) + (invu[3] * tcj) + (invu[6] * tck);
                  const T pty =                   (invu[4] * tcj) + (invu[7] * tck);
                  const T ptz =                                     (invu[8] * tck);
                  const ParticlePointInfluence ppi = outerShellEffect(cgw, poly_ps, pos, i, j, k,
                                                                      standoff, ptx, pty, ptz,
                                                                      uform, ew_coeff, false);
                  const T approx_u = s_tric.evaluate(ptx, pty, ptz);
                  const T3 approx_f = s_tric.template derivative<T3>(ptx, pty, ptz);
                  const size_t ijk_idx = (((ck * 32) + cj) * 32) + ci;
                  act_u[ijk_idx] = ppi.u;
                  act_fx[ijk_idx] = ppi.dudx;
                  act_fy[ijk_idx] = ppi.dudy;
                  act_fz[ijk_idx] = ppi.dudz;
                  apr_u[ijk_idx] = approx_u;
                  apr_fx[ijk_idx] = approx_f.x;
                  apr_fy[ijk_idx] = approx_f.y;
                  apr_fz[ijk_idx] = approx_f.z;
                }
              }
            }
            for (int ck = 0; ck < 32; ck++) {
              printf("act_u(:,:,%d) = [\n", ck + 1);
              for (int ci = 0; ci < 32; ci++) {
                for (int cj = 0; cj < 32; cj++) {
                  const size_t ijk_idx = (((ck * 32) + cj) * 32) + ci;
                  printf(" %9.5lf", act_u[ijk_idx]);
                }
                printf("\n");
              }
              printf("];\n");
            }
            for (int ck = 0; ck < 32; ck++) {
              printf("apr_u(:,:,%d) = [\n", ck + 1);
              for (int ci = 0; ci < 32; ci++) {
                for (int cj = 0; cj < 32; cj++) {
                  const size_t ijk_idx = (((ck * 32) + cj) * 32) + ci;
                  printf(" %9.5lf", apr_u[ijk_idx]);
                }
                printf("\n");
              }
              printf("];\n");
            }
          }
          // END CHECK
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Baseline variables
  TestEnvironment oe(argc, argv, ExceptionResponse::SILENT);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  StopWatch timer;

  // Take in additional inputs
  double cutoff = 10.0;
  for (int i = 0; i < argc; i++) {
    if (i < argc - 1 && strcmpCased(argv[i], "-cut", CaseSensitivity::NO)) {
      if (verifyNumberFormat(argv[i + 1], NumberFormat::STANDARD_REAL) ||
          verifyNumberFormat(argv[i + 1], NumberFormat::SCIENTIFIC)) {
        cutoff = stod(std::string(argv[i + 1]));
      }
      else {
        rtErr("The -cut optional keyword must be followed by a real number.", "cellgrid");
      }
    }
  }
  const double half_cutoff = 0.5 * cutoff;
  const std::vector<std::string> systems = { "trpcage_in_water", "drug_example", "ubiquitin" };
  const char osc = osSeparator();
  const std::string base_top_path = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_path = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  TestSystemManager tsm(base_top_path, "top", systems, base_crd_path, "inpcrd", systems);

  // Test the HAIL approximation with various potentials
  testHAIL(tsm, 0.0, 671285832, DecomposablePotential::ELECTROSTATIC, 5.25, 1);
  testHAIL(tsm, 0.0, 981673147, DecomposablePotential::ELEC_PME_DIRECT, 5.25, 1);
  testHAIL(tsm, 0.0, 410491830, DecomposablePotential::DISPERSION, 10.25, 1);

  // Test the Standoff approximation with various potentials
  testStandoff<double, double3, double4>(tsm, 0.0, 918639418,
                                         DecomposablePotential::ELEC_PME_DIRECT, 5.25, 0.3, 1);
  
  // Print results
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
}
