#include "copyright.h"
#include "map_density.h"

namespace stormm {
namespace energy {

using stmath::sanctionedDensityGridDimension;

//-------------------------------------------------------------------------------------------------
void matchThemes(const NonbondedTheme pm_theme, const NonbondedTheme cg_theme) {
  bool problem = false;
  switch (pm_theme) {
  case NonbondedTheme::ELECTROSTATIC:
    switch (cg_theme) {
    case NonbondedTheme::ELECTROSTATIC:
    case NonbondedTheme::ALL:
      break;
    case NonbondedTheme::VAN_DER_WAALS:
      problem = true;
    }
    break;
  case NonbondedTheme::VAN_DER_WAALS:
    switch (cg_theme) {
    case NonbondedTheme::ELECTROSTATIC:
      problem = true;
    case NonbondedTheme::VAN_DER_WAALS:
    case NonbondedTheme::ALL:
      break;
    }
    break;
  case NonbondedTheme::ALL:
    rtErr("A particle-mesh interaction grid cannot carry more than one kind of non-bonded "
          "potential.", "matchThemes");
    break;
  }
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void mapDensity(PMIGridWriter *pm_wrt, PMIGridAccumulator *pm_acc, MMControlKit<double> *ctrl,
                const CellGridReader<void, void, void, void> &v_cgr, const size_t cg_tmat,
                const SyNonbondedKit<double, double2> &synbk, const int block_count, const int2 lp,
                const QMapMethod approach, PMIGrid *pm) {
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  switch (approach) {
  case QMapMethod::ACC_SHARED:
    launchShrAccDensityKernel(pm_wrt, pm->useOverflowAccumulation(), ctrl, v_cgr, cg_tmat, synbk,
                              lp);
    break;
  case QMapMethod::GENERAL_PURPOSE:
    launchPMIGridInitialization(pm_acc, block_count);
    launchGenPrpDensityKernel(pm_acc, v_cgr, cg_tmat, synbk, lp);
    launchPMIGridRealConversion(pm_wrt, *pm_acc, block_count);    
    break;
  case QMapMethod::AUTOMATIC:
    break;
  }

  // Mark that the particle-mesh interaction grids are now presented in real format
  pm->setRealDataFormat();
}

//-------------------------------------------------------------------------------------------------
void mapDensity(PMIGridWriter *pm_wrt, PMIGridAccumulator *pm_acc, MMControlKit<float> *ctrl,
                const CellGridReader<void, void, void, void> &v_cgr, const size_t cg_tmat,
                const SyNonbondedKit<float, float2> &synbk, const int block_count, const int2 lp,
                const QMapMethod approach, PMIGrid *pm) {
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  switch (approach) {
  case QMapMethod::ACC_SHARED:
    launchShrAccDensityKernel(pm_wrt, pm->useOverflowAccumulation(), ctrl, v_cgr, cg_tmat, synbk,
                              lp);
    break;
  case QMapMethod::GENERAL_PURPOSE:
    launchPMIGridInitialization(pm_acc, block_count);
    launchGenPrpDensityKernel(pm_acc, v_cgr, cg_tmat, synbk, lp);
    launchPMIGridRealConversion(pm_wrt, *pm_acc, block_count);    
    break;
  case QMapMethod::AUTOMATIC:
    break;
  }

  // Mark that the particle-mesh interaction grids are now presented in real format
  pm->setRealDataFormat();
}
#endif

//-------------------------------------------------------------------------------------------------
void mapDensity(PMIGrid *pm, const AtomGraphSynthesis *poly_ag) {

  // Extract the cell grid pointer and unroll its templating
  const size_t cg_tmat  = pm->getCellGridMatrixTypeID();
  const size_t cg_tacc  = pm->getCellGridAccumulatorTypeID();
  const size_t cg_tcalc = pm->getCellGridCalculationTypeID();
  if (cg_tmat == double_type_index) {

    // The type of the cell dimension matrices implies the format of the coordinate data (a
    // four-tuple of the matrix dimension data type).
    unrollMapDensityCall<double, double4>(pm, cg_tacc, cg_tcalc, poly_ag);
  }
  else if (cg_tmat == float_type_index) {
    unrollMapDensityCall<float, float4>(pm, cg_tacc, cg_tcalc, poly_ag);
  }
  else if (cg_tmat == llint_type_index) {
    unrollMapDensityCall<llint, llint4>(pm, cg_tacc, cg_tcalc, poly_ag);
  }
  else if (cg_tmat == int_type_index) {
    unrollMapDensityCall<int, int4>(pm, cg_tacc, cg_tcalc, poly_ag);
  }
}

//-------------------------------------------------------------------------------------------------
void mapDensity(PMIGrid *pm, const AtomGraphSynthesis &poly_ag) {
  mapDensity(pm, poly_ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
std::vector<double> mapDensity(const CoordinateFrame *cf, const AtomGraph *ag,
                               const NonbondedTheme theme, const FFTMode fft_staging,
                               const int grid_dim_a, const int grid_dim_b, const int grid_dim_c,
                               const int order) {
  int actual_grid_dim_a = grid_dim_a;
  int actual_grid_dim_b = grid_dim_b;
  int actual_grid_dim_c = grid_dim_c;
  const CoordinateFrameReader cfr = cf->data();
  if (grid_dim_a < 0 || grid_dim_b < 0 || grid_dim_c < 0) {
    const ullint big_product = ipowl(2, 14) * ipowl(3, 10) * ipowl(5, 6) * 7LL * 11LL;
    const std::vector<uint> primes = { 2, 3, 5, 7, 11 };
    double max_spacing = 1.00;
    switch (order) {
    case 5:
      max_spacing = 1.25;
      break;
    case 6:
      max_spacing = 1.50;
      break;
    case 4:
    default:
      break;
    }
    if (grid_dim_a < 0) {
      actual_grid_dim_a = sanctionedDensityGridDimension(big_product, cfr.invu, UnitCellAxis::A,
                                                         primes, max_spacing);
    }
    if (grid_dim_b < 0) {
      actual_grid_dim_b = sanctionedDensityGridDimension(big_product, cfr.invu, UnitCellAxis::B,
                                                         primes, max_spacing);
    }
    if (grid_dim_c < 0) {
      actual_grid_dim_c = sanctionedDensityGridDimension(big_product, cfr.invu, UnitCellAxis::C,
                                                         primes, max_spacing);
    }
  }
  return mapDensity<double>(cfr, ag->getDoublePrecisionNonbondedKit(), theme, fft_staging,
                            actual_grid_dim_a, actual_grid_dim_b, actual_grid_dim_c, order);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> mapDensity(const CoordinateFrame &cf, const AtomGraph &ag,
                               const NonbondedTheme theme, const FFTMode fft_staging,
                               const int grid_dim_a, const int grid_dim_b, const int grid_dim_c,
                               const int order) {
  return mapDensity<double>(cf.data(), ag.getDoublePrecisionNonbondedKit(), theme, fft_staging,
                            grid_dim_a, grid_dim_b, grid_dim_c, order);
}

} // namespace energy
} // namespace stormm
