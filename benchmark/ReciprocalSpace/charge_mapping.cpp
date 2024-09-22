#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Accelerator/core_kernel_manager.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Constants/behavior.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/bspline.h"
#include "../../src/Math/math_enumerators.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Namelists/command_line_parser.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/parsing_enumerators.h"
#include "../../src/Potential/cellgrid.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/map_density.h"
#include "../../src/Potential/pmigrid.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/synthesis_abstracts.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/stopwatch.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::energy;
using namespace stormm::namelist;
using namespace stormm::parse;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Test the accuracy of B-spline computations with various orders of interpolation.  How conserved
// is the density when mapped to the particle-mesh interaction grid?
//
// Arguments:
//   apply_unity:  Apply the fact that B-splines form a smooth partition of unity when computing
//                 them.
//-------------------------------------------------------------------------------------------------
void singleBSplineTests(const bool apply_unity) {
  const int low_order_limit = 4;
  const int high_order_limit = 8;
  const int norders = high_order_limit - low_order_limit + 1;
  std::vector<double> lp_mean_qcons(norders), dp_mean_qcons(norders), sp_mean_qcons(norders);
  std::vector<double> lp_stdv_qcons(norders), dp_stdv_qcons(norders), sp_stdv_qcons(norders);
  std::vector<double> dp_mue_force(norders), sp_mue_force(norders);
  std::vector<double> dp_stdv_force(norders), sp_stdv_force(norders);
  for (int ordr = low_order_limit; ordr <= high_order_limit; ordr++) {
    std::vector<long double> lp_coef(ordr * 32), lp_dcoef(ordr * 32);
    long double* lpc_ptr  = lp_coef.data();
    long double* lpdc_ptr = lp_dcoef.data();
    std::vector<double> dp_coef(ordr * 32), dp_dcoef(ordr * 32);
    double* dpc_ptr  = dp_coef.data();
    double* dpdc_ptr = dp_dcoef.data();
    std::vector<float> sp_coef(ordr * 32), sp_dcoef(ordr * 32);
    float* spc_ptr  = sp_coef.data();
    float* spdc_ptr = sp_dcoef.data();
    long double lx = 0.0;
    double dx      = 0.0;
    float fx       = 0.0;
    const int npts = 32;
    const double l_incr = static_cast<long double>(1.0) / static_cast<long double>(npts);
    const double d_incr = 1.0 / static_cast<double>(npts);
    const float f_incr = static_cast<float>(1.0) / static_cast<float>(npts);
    for (int i = 0; i < npts; i++) {
      if (apply_unity) {
        bSpline<long double>(lx, ordr, &lpc_ptr[ordr * i], &lpdc_ptr[ordr * i]);
        bSpline<double>(dx, ordr, &dpc_ptr[ordr * i], &dpdc_ptr[ordr * i]);
        bSpline<float>(fx, ordr, &spc_ptr[ordr * i], &spdc_ptr[ordr * i]);
      }
      else {
        const std::vector<long double> t_l = bSplineNoUnity<long double>(lx, ordr);
        const std::vector<double> t_d      = bSplineNoUnity<double>(dx, ordr);
        const std::vector<float> t_f       = bSplineNoUnity<float>(fx, ordr);
        const std::vector<long double> dt_l = dBSpline<long double>(lx, ordr, false);
        const std::vector<double> dt_d      = dBSpline<double>(dx, ordr, false);
        const std::vector<float> dt_f       = dBSpline<float>(fx, ordr, false);
        const int joffset = ordr * i;
        for (int j = joffset; j < ordr * (i + 1); j++) {
          lpc_ptr[j] = t_l[j - joffset];
          dpc_ptr[j] = t_d[j - joffset];
          spc_ptr[j] = t_f[j - joffset];
          lpdc_ptr[j] = dt_l[j - joffset];
          dpdc_ptr[j] = dt_d[j - joffset];
          spdc_ptr[j] = dt_f[j - joffset];
        }
      }
      lx += l_incr;
      dx += d_incr;
      fx += f_incr;
    }
    const int npts3 = npts * npts * npts;
    std::vector<long double> lp_sums(npts3);
    std::vector<double> lp_dsums(3 * npts3);
    std::vector<double> dp_sums(npts3), dp_dsums(3 * npts3);
    std::vector<double> sp_sums(npts3), sp_dsums(3 * npts3);
    for (int i = 0; i < npts; i++) {
      const long double *lp_icof_ptr = &lpc_ptr[ordr * i];
      const double      *dp_icof_ptr = &dpc_ptr[ordr * i];
      const float       *sp_icof_ptr = &spc_ptr[ordr * i];
      const long double *lp_dicof_ptr = &lpdc_ptr[ordr * i];
      const double      *dp_dicof_ptr = &dpdc_ptr[ordr * i];
      const float       *sp_dicof_ptr = &spdc_ptr[ordr * i];
      for (int j = 0; j < npts; j++) {
        const long double *lp_jcof_ptr = &lpc_ptr[ordr * j];
        const double      *dp_jcof_ptr = &dpc_ptr[ordr * j];
        const float       *sp_jcof_ptr = &spc_ptr[ordr * j];
        const long double *lp_djcof_ptr = &lpdc_ptr[ordr * j];
        const double      *dp_djcof_ptr = &dpdc_ptr[ordr * j];
        const float       *sp_djcof_ptr = &spdc_ptr[ordr * j];
        for (int k = 0; k < npts; k++) {
          const long double *lp_kcof_ptr = &lpc_ptr[ordr * k];
          const double      *dp_kcof_ptr = &dpc_ptr[ordr * k];
          const float       *sp_kcof_ptr = &spc_ptr[ordr * k];
          const long double *lp_dkcof_ptr = &lpdc_ptr[ordr * k];
          const double      *dp_dkcof_ptr = &dpdc_ptr[ordr * k];
          const float       *sp_dkcof_ptr = &spdc_ptr[ordr * k];
          long double lp_sum = 0.0;
          double dp_sum      = 0.0;
          float sp_sum       = 0.0;
          long double lp_dxsum = 0.0;
          long double lp_dysum = 0.0;
          long double lp_dzsum = 0.0;
          double dp_dxsum      = 0.0;
          double dp_dysum      = 0.0;
          double dp_dzsum      = 0.0;
          float sp_dxsum       = 0.0;
          float sp_dysum       = 0.0;
          float sp_dzsum       = 0.0;
          for (int ni = 0; ni < ordr; ni++) {
            for (int nj = 0; nj < ordr; nj++) {
              for (int nk = 0; nk < ordr; nk++) {
                lp_sum += lp_icof_ptr[ni] * lp_jcof_ptr[nj] * lp_kcof_ptr[nk];
                dp_sum += dp_icof_ptr[ni] * dp_jcof_ptr[nj] * dp_kcof_ptr[nk];
                sp_sum += sp_icof_ptr[ni] * sp_jcof_ptr[nj] * sp_kcof_ptr[nk];
                lp_dxsum += lp_dicof_ptr[ni] * lp_jcof_ptr[nj] * lp_kcof_ptr[nk];
                dp_dxsum += dp_dicof_ptr[ni] * dp_jcof_ptr[nj] * dp_kcof_ptr[nk];
                sp_dxsum += sp_dicof_ptr[ni] * sp_jcof_ptr[nj] * sp_kcof_ptr[nk];
                lp_dysum += lp_icof_ptr[ni] * lp_djcof_ptr[nj] * lp_kcof_ptr[nk];
                dp_dysum += dp_icof_ptr[ni] * dp_djcof_ptr[nj] * dp_kcof_ptr[nk];
                sp_dysum += sp_icof_ptr[ni] * sp_djcof_ptr[nj] * sp_kcof_ptr[nk];
                lp_dzsum += lp_icof_ptr[ni] * lp_jcof_ptr[nj] * lp_dkcof_ptr[nk];
                dp_dzsum += dp_icof_ptr[ni] * dp_jcof_ptr[nj] * dp_dkcof_ptr[nk];
                sp_dzsum += sp_icof_ptr[ni] * sp_jcof_ptr[nj] * sp_dkcof_ptr[nk];
              }
            }
          }
          const size_t dest_idx = (((k * npts) + j) * npts) + i;
          lp_sums[dest_idx] = lp_sum - 1.0;
          dp_sums[dest_idx] = dp_sum - 1.0;
          sp_sums[dest_idx] = sp_sum - 1.0;
          lp_dsums[ 3 * dest_idx     ] = lp_dxsum;
          lp_dsums[(3 * dest_idx) + 1] = lp_dysum;
          lp_dsums[(3 * dest_idx) + 2] = lp_dzsum;
          dp_dsums[ 3 * dest_idx     ] = dp_dxsum;
          dp_dsums[(3 * dest_idx) + 1] = dp_dysum;
          dp_dsums[(3 * dest_idx) + 2] = dp_dzsum;
          sp_dsums[ 3 * dest_idx     ] = sp_dxsum;
          sp_dsums[(3 * dest_idx) + 1] = sp_dysum;
          sp_dsums[(3 * dest_idx) + 2] = sp_dzsum;
        }
      }
    }
    const VarianceMethod vmeth = VarianceMethod::STANDARD_DEVIATION;
    const size_t ordr_idx = ordr - low_order_limit;
    lp_mean_qcons[ordr_idx] = mean(lp_sums);
    dp_mean_qcons[ordr_idx] = mean(dp_sums);
    sp_mean_qcons[ordr_idx] = mean(sp_sums);
    lp_stdv_qcons[ordr_idx] = variance(lp_sums, VarianceMethod::STANDARD_DEVIATION);
    dp_stdv_qcons[ordr_idx] = variance(dp_sums, VarianceMethod::STANDARD_DEVIATION);
    sp_stdv_qcons[ordr_idx] = variance(sp_sums, VarianceMethod::STANDARD_DEVIATION);
    dp_mue_force[ordr_idx]  = meanUnsignedError(lp_dsums, dp_dsums);
    sp_mue_force[ordr_idx]  = meanUnsignedError(lp_dsums, sp_dsums);
    for (int i = 0; i < 3 * npts3; i++) {
      dp_dsums[i] -= lp_dsums[i];
      sp_dsums[i] -= lp_dsums[i];
    }
    dp_stdv_force[ordr_idx] = variance(dp_dsums, VarianceMethod::STANDARD_DEVIATION);
    sp_stdv_force[ordr_idx] = variance(sp_dsums, VarianceMethod::STANDARD_DEVIATION);
  }

  printf(" B-spline charge conservation\n");
  printf(" Order    Long Double (80-bit)            Double (64-bit)               "
         "Float (32-bit)\n");
  for (int ordr = low_order_limit; ordr <= high_order_limit; ordr++) {
    const size_t ordr_idx = ordr - low_order_limit;
    printf(" %2d   %14.6e %12.6e  %14.6e %12.6e  %14.6e %12.6e\n", ordr,
           lp_mean_qcons[ordr_idx], lp_stdv_qcons[ordr_idx], dp_mean_qcons[ordr_idx],
           dp_stdv_qcons[ordr_idx], sp_mean_qcons[ordr_idx], sp_stdv_qcons[ordr_idx]);
  }
  printf("\n B-spline derivative accuracy\n");
  printf(" Order    Long Double (80-bit)            Double (64-bit)               "
         "Float (32-bit)\n");
  for (int ordr = low_order_limit; ordr <= high_order_limit; ordr++) {
    const size_t ordr_idx = ordr - low_order_limit;
    printf("                                   %14.6e %12.6e  %14.6e %12.6e\n",
           dp_mue_force[ordr_idx], dp_stdv_force[ordr_idx], sp_mue_force[ordr_idx],
           sp_stdv_force[ordr_idx]);
  }  
}

//-------------------------------------------------------------------------------------------------
// Compute the charge density for a series of systems, mapping first in double-precision and then
// using various other precision models, whether in the accumulation or in the coordinate
// representation of the particles.
//
// Arguments:
//   poly_ps:            The coordinate synthesis spanning all systems
//   poly_ag:            The topology synthesis spanning all systems
//   order:              The order of interpolation
//   unification:        Method for B-spline evaluation
//   half_cutoff:        Indicate the minimum width of the spatial decomposition cells
//   cutoff_pad:         Minimum padding added to each cell width (in real simulations, this
//                       ensures that the cells are unlikely to shrink beneath the cutoff
//                       requirements)
//   mesh_subdivisions:  The number of times to subdivide each spatial decomposition cell, along
//                       each axis, to obtain the particle-mesh interaction grid
//-------------------------------------------------------------------------------------------------
void benchmarkChargeDensity(const PhaseSpaceSynthesis &poly_ps, const AtomGraphSynthesis &poly_ag,
                            const int order,
                            const BSplineUnity unification = BSplineUnity::CENTER_FILL,
                            const double half_cutoff = 4.5, const double cutoff_pad = 0.15,
                            const int mesh_subdivisions = 4) {
  const int nsys = poly_ps.getSystemCount();
  std::vector<std::vector<double>> bnch_density(nsys), fcrd_dcalc_density(nsys);
  std::vector<std::vector<float>> fcrd_fcalc_density(nsys), dcrd_fcalc_density(nsys);
  std::vector<std::vector<double>> d_fcrd_fcalc_density, d_dcrd_fcalc_density;

  // Create a cell grid with as much double-precision content as possible.
  CellGrid<double, llint, double, double4> cg_dd(poly_ps.getSelfPointer(),
                                                 poly_ag.getSelfPointer(), half_cutoff, cutoff_pad,
                                                 mesh_subdivisions, NonbondedTheme::ELECTROSTATIC);
  const CellGridReader cg_ddr = cg_dd.data();
  
  // Loop over all systems, creating the benchmark densities with double-precision calculations.
  d_fcrd_fcalc_density.reserve(nsys);
  d_dcrd_fcalc_density.reserve(nsys);
  printf("B-Spline Order: %d (unification %s)\n", order, getEnumerationName(unification).c_str());
  printf("   64-bit Calc / 32-bit Coord   32-bit Calc / 64-bit Coord   32-bit Calc / 32-bit Coord"
         "\n");
  for (int i = 0; i < nsys; i++) {
    CoordinateFrame cf = poly_ps.exportCoordinates(i);
    const AtomGraph *ag_ptr = poly_ag.getSystemTopologyPointer(i);
    const ullint gdims = cg_ddr.system_cell_grids[i];
    const int na = ((gdims >> 28) & 0xfff) * cg_ddr.mesh_ticks;
    const int nb = ((gdims >> 40) & 0xfff) * cg_ddr.mesh_ticks;
    const int nc = (gdims >> 52) * cg_ddr.mesh_ticks;
    bnch_density[i] = mapDensity(&cf, ag_ptr,  NonbondedTheme::ELECTROSTATIC,
                                 FFTMode::OUT_OF_PLACE, na, nb, nc, order);

    // Create a 32-bit representation of the coordinates and recalculate the density in single-
    // and double-precision with those coordinates.  Which approximation does more damage to the
    // outcome?
    CoordinateFrame cf_32t = cf;
    CoordinateFrameWriter cfr = cf.data();
    CoordinateFrameWriter cfr_32t = cf_32t.data();
    for (int i = 0; i < cfr.natom; i++) {
      const float fx = cfr.xcrd[i];
      const float fy = cfr.ycrd[i];
      const float fz = cfr.zcrd[i];
      cfr_32t.xcrd[i] = fx;
      cfr_32t.ycrd[i] = fy;
      cfr_32t.zcrd[i] = fz;
    }
    fcrd_dcalc_density[i] = mapDensity(&cf_32t, ag_ptr, NonbondedTheme::ELECTROSTATIC,
                                       FFTMode::OUT_OF_PLACE, na, nb, nc, order);
    const NonbondedKit<float> nbk_f = ag_ptr->getSinglePrecisionNonbondedKit();
    fcrd_fcalc_density[i] = mapDensity<float>(cf_32t.data(), nbk_f, NonbondedTheme::ELECTROSTATIC,
                                              FFTMode::OUT_OF_PLACE, na, nb, nc, order,
                                              unification);
    dcrd_fcalc_density[i] = mapDensity<float>(cf.data(), nbk_f, NonbondedTheme::ELECTROSTATIC,
                                              FFTMode::OUT_OF_PLACE, na, nb, nc, order,
                                              unification);
    d_fcrd_fcalc_density.emplace_back(fcrd_fcalc_density[i].begin(), fcrd_fcalc_density[i].end());
    d_dcrd_fcalc_density.emplace_back(dcrd_fcalc_density[i].begin(), dcrd_fcalc_density[i].end());

    // Compute the difference vectors
    const size_t ngrds = na * nb * nc;
    std::vector<double> fcrd_dcalc_delta(ngrds), dcrd_fcalc_delta(ngrds), fcrd_fcalc_delta(ngrds);
    for (size_t j = 0; j < ngrds; j++) {
      fcrd_dcalc_delta[j] = fcrd_dcalc_density[i][j] - bnch_density[i][j];
      dcrd_fcalc_delta[j] = dcrd_fcalc_density[i][j] - bnch_density[i][j];
      fcrd_fcalc_delta[j] = fcrd_fcalc_density[i][j] - bnch_density[i][j];
    }
    
    // Print the results of single-precision coordinates, calculations, or both
    printf("  %13.6e %13.6e  %13.6e %13.6e  %13.6e %13.6e\n",
           mean(fcrd_dcalc_delta), variance(fcrd_dcalc_delta, VarianceMethod::STANDARD_DEVIATION),
           mean(dcrd_fcalc_delta), variance(dcrd_fcalc_delta, VarianceMethod::STANDARD_DEVIATION),
           mean(fcrd_fcalc_delta), variance(fcrd_fcalc_delta, VarianceMethod::STANDARD_DEVIATION));
  }
  printf("\n");
  
  // Create a single-precision cell grid, map the densities, and compute the errors arising from
  // this level of approximation.  The local coordinate axes of the spatial decomposition cells
  // afford an advantage for positional rounding.
  CellGrid<float, llint, float, float4> cg_ff(poly_ps.getSelfPointer(),
                                              poly_ag.getSelfPointer(), half_cutoff, cutoff_pad,
                                              mesh_subdivisions, NonbondedTheme::ELECTROSTATIC);
  const std::vector<PrecisionModel> models = { PrecisionModel::DOUBLE, PrecisionModel::SINGLE };
  std::vector<std::vector<double>> cg_mean_errors(models.size()), cg_stdev(models.size());
  for (size_t calc = 0; calc < models.size(); calc++) {
    PMIGrid pm_fx(cg_ff, NonbondedTheme::ELECTROSTATIC, order, models[calc]);
    mapDensity(&pm_fx, cg_ff, poly_ag);
    PMIGridReader pm_fxr = pm_fx.data();
    cg_mean_errors[calc] = std::vector<double>(nsys);
    cg_stdev[calc] = std::vector<double>(nsys);
    for (int i = 0; i < nsys; i++) {
      std::vector<double> pmi_extract = pm_fx.getGrid(i);
      const uint4 i_dims = pm_fxr.dims[i];
      const size_t npmg = i_dims.x * i_dims.y * i_dims.z;
      for (size_t j = 0; j < npmg; j++) {
        pmi_extract[j] -= bnch_density[i][j];
      }
      cg_mean_errors[calc][i] = mean(pmi_extract);
      cg_stdev[calc][i] = variance(pmi_extract, VarianceMethod::STANDARD_DEVIATION);
    }
  }

  // Present the final table of results
  printf("CellGrid 32-bit local coordinate representation:\n");
  printf("   64-bit Calc / 32-bit Coord                                32-bit Calc / 32-bit Coord"
         "\n");
  for (int i = 0; i < nsys; i++) {
    printf("  %13.6e %13.6e                               %13.6e %13.6e\n", cg_mean_errors[0][i],
           cg_stdev[0][i], cg_mean_errors[1][i], cg_stdev[1][i]);
  }
  printf("\n");
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
// Run timings benchmarks to determine the cost of launching GPU kernels.
//
// Arguments:
//   poly_ps:             The synthesis of coordinates
//   poly_ag:             The synthesis of topologies
//   timer:               Tracks the wall time
//   gpu:                 Specifications of the GPU that will perform the work
//   category_label:      Basic label to apply to timings from various tests.  This should describe
//                        the nature of the system or collection of systems, to which various
//                        markers for the precision model will be appended.
//   half_cutoff:         Half of the non-bonded cutoff for particle interactions computed in real
//                        space
//   cutoff_pad:          The padding applied to each minimum cell width (half_cutoff) ensuring
//                        that contraction of the box does not cause the cutoff to be violated
//                        with the pre-calculated partitioning
//   mesh_subdivisions:   The number of particle-mesh interaction grid intervals lying in each
//                        spatial decomposition cell along any unit cell axis
//   n_trials:            The number of times to repeat each experiment
//   n_repeats:           The number of times to repeat each kernel launch within a single test,
//                        each of which is terminated by a device synchronization call
//   fp_bit_count:        The number of fixed-precision bits in which density is accumulated
//-------------------------------------------------------------------------------------------------
template <typename T, typename Tcalc, typename T4>
void benchmarkAccumulationKernels(const PhaseSpaceSynthesis &poly_ps,
                                  const AtomGraphSynthesis &poly_ag,
                                  const int order, StopWatch *timer, const GpuDetails &gpu,
                                  const std::string &category_label,
                                  const double half_cutoff = 4.5, const double cutoff_pad = 0.15,
                                  const int mesh_subdivisions = 4, const int n_trials = 4,
                                  const int n_repeats = 100, const int fp_bit_count = 32) {
  CoreKlManager launcher(gpu, poly_ag);
  CellGrid<T, llint, Tcalc, T4> cg(poly_ps.getSelfPointer(), poly_ag.getSelfPointer(), half_cutoff,
                                   cutoff_pad, mesh_subdivisions, NonbondedTheme::ELECTROSTATIC);
  cg.upload();
  const std::vector<PrecisionModel> models = { PrecisionModel::DOUBLE, PrecisionModel::SINGLE };
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const PrecisionModel prec = (tcalc_is_double) ? PrecisionModel::DOUBLE : PrecisionModel::SINGLE;
  const size_t cg_tmat = std::type_index(typeid(T)).hash_code();
  PMIGrid pm_fx(cg, NonbondedTheme::ELECTROSTATIC, order, prec, FFTMode::OUT_OF_PLACE,
                fp_bit_count, fp_bit_count, gpu);
  pm_fx.prepareWorkUnits(QMapMethod::ACC_SHARED, gpu);
  pm_fx.upload();
  timer->assignTime(0);
  
  // Create resources needed by some density mapping kernels.
  MolecularMechanicsControls mm_ctrl;
  mm_ctrl.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                                ClashResponse::NONE, VwuGoal::ACCUMULATE, prec, prec,
                                pm_fx.getWorkUnitConfiguration(), pm_fx.getMode(), cg_tmat, order,
                                NeighborListKind::MONO, TinyBoxPresence::NO, poly_ag);
  mm_ctrl.upload();
  MMControlKit<double> ctrl_d = mm_ctrl.dpData(devc_tier);
  MMControlKit<float>  ctrl_f = mm_ctrl.spData(devc_tier);
    
  // Add categories to the timer, if they are not already present
  const int init_timings = timer->addCategory(category_label + ", Init");
  const int conv_timings = timer->addCategory(category_label + ", Convert");

  // Time the initialization and conversion
  PMIGridAccumulator pm_acc = pm_fx.fpData(devc_tier);
  PMIGridWriter pm_wrt = pm_fx.data(devc_tier);
  for (int trial = 0; trial < n_trials; trial++) {
    for (int i = 0; i < n_repeats; i++) {
      launchPMIGridInitialization(&pm_acc, gpu);
    }
    pm_fx.setRealDataFormat(false);
    cudaDeviceSynchronize();
    timer->assignTime(init_timings);
    for (int i = 0; i < n_repeats; i++) {
      launchPMIGridRealConversion(&pm_wrt, pm_acc, gpu);
    }
    pm_fx.setRealDataFormat();
    cudaDeviceSynchronize();
    timer->assignTime(conv_timings);
  }
  
  // Time the entire process with different kernel implementations.
  const int genp_timings = timer->addCategory(category_label + ", General Purpose");
  const int sacc_timings = timer->addCategory(category_label + ", Acc. __shared__");
  const CellGridReader<void, void, void, void> v_cgr = cg.templateFreeData(devc_tier);
  const int2 lp_genp = launcher.getDensityMappingKernelDims(QMapMethod::GENERAL_PURPOSE, prec,
                                                            cg_tmat, pm_wrt.order);
  const int2 lp_sacc = launcher.getDensityMappingKernelDims(QMapMethod::ACC_SHARED, prec,
                                                            pm_wrt.mode,
                                                            pm_fx.useOverflowAccumulation(),
                                                            cg_tmat, pm_wrt.order);
  const SyNonbondedKit<double,
                       double2> synbk_d = poly_ag.getDoublePrecisionNonbondedKit(devc_tier);
  const SyNonbondedKit<float,
                       float2> synbk_f = poly_ag.getSinglePrecisionNonbondedKit(devc_tier);
  for (int trial = 0; trial < n_trials; trial++) {
    if (tcalc_is_double) {
      timer->assignTime(0);
      for (int i = 0; i < n_repeats; i++) {
        mapDensity(&pm_wrt, &pm_acc, nullptr, v_cgr, cg_tmat, synbk_d, gpu.getSMPCount(),
                   lp_genp, QMapMethod::GENERAL_PURPOSE, &pm_fx);
      }
    }
    else {
      timer->assignTime(0);
      for (int i = 0; i < n_repeats; i++) {
        mapDensity(&pm_wrt, &pm_acc, nullptr, v_cgr, cg_tmat, synbk_f, gpu.getSMPCount(),
                   lp_genp, QMapMethod::GENERAL_PURPOSE, &pm_fx);
      }
    }
    cudaDeviceSynchronize();
    timer->assignTime(genp_timings);
  }  

  // Download the results.  Check to ensure that each grid is identical and that they agree with
  // the results of a CPU-based routine.  Remember the CPU-based results to check other kernels.
  pm_fx.download();
  const PMIGridWriter pm_wrt_host = pm_fx.data();
  std::vector<std::vector<double>> cpu_grids(poly_ps.getSystemCount());
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    const CoordinateFrame cf = poly_ps.exportCoordinates(i);
    cpu_grids[i] = mapDensity(&cf, poly_ps.getSystemTopologyPointer(i),
                              NonbondedTheme::ELECTROSTATIC, FFTMode::OUT_OF_PLACE,
                              pm_wrt_host.dims[i].x, pm_wrt_host.dims[i].y, pm_wrt_host.dims[i].z,
                              order);
    const std::vector<double> gpu_grid = pm_fx.getGrid(i);
    const double mue_i = meanUnsignedError(cpu_grids[i], gpu_grid);
    if (mue_i > 1.0e-4) {
      printf("Mean unsigned error in replica %2d is %9.6lf after implementing the %s kernel in "
             "%s precision.\n", i, mue_i, getEnumerationName(QMapMethod::GENERAL_PURPOSE).c_str(),
             getEnumerationName(prec).c_str());
    }
  }

  // Run the __shared__ density accumulation kernel.
  for (int trial = 0; trial < n_trials; trial++) {
    timer->assignTime(0);
    if (tcalc_is_double) {
      for (int i = 0; i < n_repeats; i++) {
        mapDensity(&pm_wrt, &pm_acc, &ctrl_d, v_cgr, cg_tmat, synbk_d, gpu.getSMPCount(), lp_sacc,
                   QMapMethod::ACC_SHARED, &pm_fx);
        ctrl_d.step += 1;
      }
    }
    else {
      for (int i = 0; i < n_repeats; i++) {
        mapDensity(&pm_wrt, &pm_acc, &ctrl_f, v_cgr, cg_tmat, synbk_f, gpu.getSMPCount(), lp_sacc,
                   QMapMethod::ACC_SHARED, &pm_fx);
        ctrl_f.step += 1;
      }
    }
    cudaDeviceSynchronize();
    timer->assignTime(sacc_timings);    
  }
  
  // Download the results.  Check to ensure that each grid is identical and that they agree with
  // the results of the CPU-based routine.
  pm_fx.download();
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    const CoordinateFrame cf = poly_ps.exportCoordinates(i);
    const std::vector<double> gpu_grid = pm_fx.getGrid(i);
    const double mue_i = meanUnsignedError(cpu_grids[i], gpu_grid);
    if (mue_i > 1.0e-4) {
      printf("Mean unsigned error in replica %2d is %9.6lf after implementing the %s kernel in "
             "%s precision.\n", i, mue_i, getEnumerationName(QMapMethod::ACC_SHARED).c_str(),
             getEnumerationName(prec).c_str());
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Make replicas of a system from the test system manager and perform charge mapping on that
// system.  Descriptions of input arguments follow from benchmarkAccumulationKernels(), above, in
// addition to those noted below.  Most inputs to this function are passed down.
//
// Arguments:
//   tsm:           A collection of test systems with periodic boundary conditions
//   system_index:  Index of the system from the test collection to use in making the synthesis
//   replicas:      The number of replicas of the chosen system to use in making the synthesis
//-------------------------------------------------------------------------------------------------
template <typename T, typename Tcalc, typename T4>
void replicateAndMapCharges(TestSystemManager *tsm, const int system_index,
                            const int replicas, const int order, StopWatch *timer,
                            const GpuDetails &gpu, const std::string &category_label,
                            const double half_cutoff = 4.5, const double cutoff_pad = 0.15,
                            const int mesh_subdivisions = 4, const int n_trials = 4,
                            const int n_repeats = 100, const int fp_bit_count = 32) {
  const std::vector<int> tiles(replicas, system_index);
  PhaseSpaceSynthesis poly_ps = tsm->exportPhaseSpaceSynthesis(tiles);
  AtomGraphSynthesis poly_ag = tsm->exportAtomGraphSynthesis(tiles);
  poly_ps.upload();
  poly_ag.upload();
  benchmarkAccumulationKernels<T, Tcalc, T4>(poly_ps, poly_ag, order, timer, gpu, category_label,
                                             half_cutoff, cutoff_pad, mesh_subdivisions, n_trials,
                                             n_repeats, fp_bit_count);
}
#endif

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Baseline variables
  StopWatch timer;
#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> force_gpu_to_engage(1);
#endif

  // Take in additional command line inputs
  CommandLineParser clip("charge_mapping", "A program to benchmark the time to map the partial "
                         "charges of atoms to the PME mesh.", { "-timings" });
  clip.addStandardAmberInputs("-p", "-c");
  clip.addStandardBenchmarkingInputs({ "-iter", "-trials", "-cutoff", "-pad", "-replicas" });
  clip.activateHelpOnNoArgs();
  clip.activateExitOnHelp();
  NamelistEmulator *t_nml = clip.getNamelistPointer();
  t_nml->addKeyword("-conservation", NamelistType::BOOLEAN);
  t_nml->addHelp("-conservation", "Estimate the conservation of charge obtained when mapping "
                 "charges to the mesh.");
  t_nml->addKeyword("-fp_bits", NamelistType::INTEGER, std::to_string(28));
  t_nml->addHelp("-fp_bits", "The number of fixed-precision bits after the point with which "
                 "charges will be accumulated at each mesh point.");
  t_nml->addKeyword("-order", NamelistType::INTEGER, std::to_string(5));
  t_nml->addHelp("-order", "The order of particle-mesh interpolation.");
  TestEnvironment oe(argc, argv, &clip, TmpdirStatus::NOT_REQUIRED, ExceptionResponse::SILENT);
  const char osc = osSeparator();
  const std::string base_top = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  t_nml->setDefaultValue("-p", base_top + osc + "test" + osc + "Topology" + osc + "kinase.top");
  t_nml->setDefaultValue("-c",
                         base_crd + osc + "test" + osc + "Trajectory" + osc + "kinase.inpcrd");
  clip.parseUserInput(argc, argv);
  const int n_trials = t_nml->getIntValue("-trials");
  const int n_repeats = t_nml->getIntValue("-iter");
  const int n_replicas = t_nml->getIntValue("-replicas");
  const int fp_bits = t_nml->getIntValue("-fp_bits");
  const int ordr = t_nml->getIntValue("-order");
  const double cutoff = t_nml->getRealValue("-cutoff");
  const double cutoff_pad = t_nml->getRealValue("-pad");
  const bool conservation_tests = t_nml->getBoolValue("-conservation");
  const std::vector<std::string> all_topologies = t_nml->getAllStringValues("-p");
  const std::vector<std::string> all_coordinates = t_nml->getAllStringValues("-c");
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  
  // A Hybrid object was created to engage the GPU.  Absorb any bootup time into "miscellaneous."
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
  }

  // Test the accuracy of B-spline charge mapping in various precision models
  if (conservation_tests) {
    printf("Taking B-splines defined as a smooth partition of unity:\n");
    singleBSplineTests(true);
    printf("\nIgnoring the partition:\n");
    singleBSplineTests(false);
  }
  
  // Take in a variety of systems and perform charge mapping tests.  Check the precision at
  // individual points in space as a function of the order, grid density, and fixed-precision
  // model.
  TestSystemManager tsm(all_topologies, all_coordinates);
  if (tsm.getTestingStatus() != TestPriority::CRITICAL) {
    rtErr("Some systems were not located.", "charge_mapping");
  }
  
  // Obtain a benchmark density map for each system at each B-spline order.
  const int nsys = tsm.getSystemCount();
  const std::vector<UnitCellType> uc_pbc = { UnitCellType::ORTHORHOMBIC, UnitCellType::TRICLINIC };
  PhaseSpaceSynthesis poly_ps = tsm.exportPhaseSpaceSynthesis(uc_pbc);
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(uc_pbc);
#ifdef STORMM_USE_HPC
  poly_ps.upload();
  poly_ag.upload();
#endif
  if (conservation_tests) {
    for (int i = 4; i <= 6; i++) {
      benchmarkChargeDensity(poly_ps, poly_ag, i, BSplineUnity::CENTER_FILL);
    }
    for (int i = 4; i <= 6; i++) {
      benchmarkChargeDensity(poly_ps, poly_ag, i, BSplineUnity::NONE);
    }
  }
#ifdef STORMM_USE_HPC
  // While the number of grid points per neighbor list cell is not always nine minus the order, for
  // orders 4, 5, and 6 the expression obtains the correct result of 5, 4, and 3 mesh points.
  const int n_grids = 9 - ordr;
  const std::string cat_lab = "Kinase (" + std::to_string(ordr) + "th order) (" +
                              std::to_string(n_replicas) + ")";
  replicateAndMapCharges<float, float, float4>(&tsm, 0, n_replicas, ordr, &timer, gpu, cat_lab,
                                               0.5 * cutoff, cutoff_pad, n_grids, n_trials,
                                               n_repeats, fp_bits);
#endif

  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  return 0;
}
