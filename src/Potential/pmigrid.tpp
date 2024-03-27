// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
PMIGrid::PMIGrid(const CellGrid<T, Tacc, Tcalc, T4> *cg_in, const NonbondedTheme theme_in,
                 const int b_spline_order_in, const PrecisionModel mode_in,
                 const FFTMode fft_staging_in, const int fp_accumulation_bits_in,
                 const int shared_fp_accumulation_bits_in, const GpuDetails &gpu,
                 const QMapMethod work_unit_configuration_in) :
    theme{theme_in}, mode{mode_in}, fft_staging{fft_staging_in},
    fp_accumulation_bits{fp_accumulation_bits_in},
    shared_fp_accumulation_bits{shared_fp_accumulation_bits_in},
    data_is_real{(fp_accumulation_bits_in > 0)},
    use_short_format_accumulation{false},
    system_count{0}, b_spline_order{b_spline_order_in},
    shared_acc_buffer_size{findSharedBufferSize(gpu)},
    recommendation{QMapMethod::GENERAL_PURPOSE},
    work_unit_configuration{work_unit_configuration_in},
    grid_dimensions{HybridKind::ARRAY, "pmig_dims"},
    capacity{0}, work_unit_count{0}, largest_work_unit_grid_points{0},
    dgrid_stack{HybridKind::ARRAY, "pmig_ddata"},
    fgrid_stack{HybridKind::ARRAY, "pmig_fdata"},
    overflow_stack{HybridKind::ARRAY, "pmig_ovrf_data"},
    work_units{HybridKind::ARRAY, "pmig_work_units"},
    cg_pointer{nullptr},
    cg_tmat{std::type_index(typeid(T)).hash_code()},
    cg_tacc{std::type_index(typeid(Tacc)).hash_code()},
    cg_tcalc{std::type_index(typeid(Tcalc)).hash_code()},
    cg_tcrd{std::type_index(typeid(T4)).hash_code()},
    poly_ps_pointer{nullptr}
{
  if (cg_in == nullptr) {
    rtErr("A valid CellGrid must be presented for array sizing purposes.", "PMIGrid");
  }
  const CellGrid<double, double, double, double4>* cgp =
    reinterpret_cast<const CellGrid<double, double, double, double4>*>(cg_in);
  cg_pointer = const_cast<CellGrid<double, double, double, double4>*>(cgp);
  system_count = cg_in->getSystemCount();
  poly_ps_pointer = const_cast<PhaseSpaceSynthesis*>(cg_in->getCoordinateSynthesisPointer());
  
  // Load the grid dimensions based on the attached CellGrid object.
  std::vector<uint> na_dims(system_count), nb_dims(system_count), nc_dims(system_count);
  std::vector<ullint> grid_totals(system_count + 1, 0);
  const int grid_mult = cg_in->getMeshSubdivisions();
  grid_dimensions.resize(system_count);
  uint4* grid_dim_ptr = grid_dimensions.data();
  for (int i = 0; i < system_count; i++) {
    na_dims[i] = cg_in->getCellCount(i, UnitCellAxis::A) * grid_mult;
    nb_dims[i] = cg_in->getCellCount(i, UnitCellAxis::B) * grid_mult;
    nc_dims[i] = cg_in->getCellCount(i, UnitCellAxis::C) * grid_mult;
    switch (fft_staging) {
    case FFTMode::IN_PLACE:
      grid_totals[i] = roundUp<ullint>(2LLU * static_cast<ullint>((na_dims[i] / 2) + 1) *
                                       nb_dims[i] * nc_dims[i], warp_size_int);
      break;
    case FFTMode::OUT_OF_PLACE:
      grid_totals[i] = roundUp<ullint>(na_dims[i] * nb_dims[i] * nc_dims[i], warp_size_int);
      break;
    }
  }
  prefixSumInPlace(&grid_totals, PrefixSumType::EXCLUSIVE);
  if (grid_totals[system_count] > UINT_MAX) {
    rtErr("The total grid volume of all systems cannot exceed " + std::to_string(UINT_MAX) +
          " points.  Total: " + std::to_string(grid_totals[system_count]) + ".  Reduce the number "
          "of systems or the density of the particle-mesh interaction grids to conserve memory.",
          "PMIGrid");
  }
  for (int i = 0; i < system_count; i++) {
    grid_dim_ptr[i] = { na_dims[i], nb_dims[i], nc_dims[i], static_cast<uint>(grid_totals[i]) };
  }
  capacity = grid_totals[system_count];

  // With the capacity set, which by design is restricted to the constructor, setting the mode and
  // fixed precision bits then controls which arrays are allocated to the necessary size.
  setMode(mode_in);
  prepareFixedPrecisionModel(fp_accumulation_bits_in, shared_fp_accumulation_bits_in);
  setRecommendedMappingMethod(gpu);
  prepareWorkUnits(gpu);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
PMIGrid::PMIGrid(const CellGrid<T, Tacc, Tcalc, T4> &cg_in, const NonbondedTheme theme_in,
                 const int b_spline_order_in, const PrecisionModel mode_in,
                 const FFTMode fft_staging_in, const int fp_accumulation_bits_in,
                 const int shared_fp_accumulation_bits_in, const GpuDetails &gpu,
                 const QMapMethod work_unit_configuration_in) :
    PMIGrid(cg_in.getSelfPointer(), theme_in, b_spline_order_in, mode_in, fft_staging_in,
            fp_accumulation_bits_in, shared_fp_accumulation_bits_in, gpu,
            work_unit_configuration_in)
{}

//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename Tcalc, typename T4>
const CellGrid<T, Tacc, Tcalc, T4>* PMIGrid::getCellGridPointer() const {
  const CellGrid<T, Tacc,
                 Tcalc, T4>* cgp = reinterpret_cast<CellGrid<T, Tacc, Tcalc, T4>*>(cg_pointer);
  return cgp;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T4> const CellGridReader<void, void, void, void>
PMIGrid::unrollTemplateFreeCGReader(const HybridTargetLevel tier) const {
  bool calc_type_problem = false;
  if (cg_tacc == int_type_index) {
    if (cg_tcalc == float_type_index) {
      const CellGrid<T, int, float, T4>* cgp = getCellGridPointer<T, int, float, T4>();
      return cgp->templateFreeData(tier);
    }
    else if (cg_tcalc == double_type_index) {
      const CellGrid<T, int, double, T4>* cgp = getCellGridPointer<T, int, double, T4>();
      return cgp->templateFreeData(tier);
    }
    else {
      calc_type_problem = true;
    }
  }
  else if (cg_tacc == llint_type_index) {
    if (cg_tcalc == float_type_index) {
      const CellGrid<T, llint, float, T4>* cgp = getCellGridPointer<T, llint, float, T4>();
      return cgp->templateFreeData(tier);
    }
    else if (cg_tcalc == double_type_index) {
      const CellGrid<T, llint, double, T4>* cgp = getCellGridPointer<T, llint, double, T4>();
      return cgp->templateFreeData(tier);
    }
    else {
      calc_type_problem = true;
    }
  }
  else {
    rtErr("The only valid types for the CellGrid's accumulation are int and llint.", "PMIGrid",
          "getTemplateFreeCellGridReader");
  }
  if (calc_type_problem) {
    rtErr("The only valid types for the CellGrid's calculations are float and double.",
          "PMIGrid", "getTemplateFreeCellGridReader");
  }
  __builtin_unreachable();
}
  
} // namespace energy
} // namespace stormm
