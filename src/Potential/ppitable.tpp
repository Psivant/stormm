// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
PPIKit<T, T4>::PPIKit(const NonbondedTheme theme_in, const BasisFunctions basis_in,
                      const TableIndexing lookup_in, const int index_bound_in,
                      const int excl_offset_in, const int index_shift_bits_in,
                      const ullint dp_detail_mask_in, const uint sp_detail_mask_in,
                      const T arg_offset_in, const T4* energy_in, const T4* force_in,
                      const T4* energy_excl_in, const T4* force_excl_in) :
    theme{theme_in}, basis{basis_in}, lookup{lookup_in}, index_bound{index_bound_in},
    excl_offset{excl_offset_in}, index_shift_bits{index_shift_bits_in},
    dp_detail_mask{dp_detail_mask_in}, sp_detail_mask{sp_detail_mask_in},
    arg_offset{arg_offset_in}, energy{energy_in}, force{force_in}, energy_excl{energy_excl_in},
    force_excl{force_excl_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
PPIeKit<T>::PPIeKit(const NonbondedTheme theme_in, const BasisFunctions basis_in,
                    const TableIndexing lookup_in, const int index_bound_in,
                    const int excl_offset_in, const int index_shift_bits_in,
                    const ullint dp_detail_mask_in, const uint sp_detail_mask_in,
                    const T arg_offset_in, const T* energy_x_in, const T* energy_y_in,
                    const T* energy_z_in, const T* energy_w_in, const T* force_x_in,
                    const T* force_y_in, const T* force_z_in, const T* force_w_in,
                    const T* energy_excl_x_in, const T* energy_excl_y_in,
                    const T* energy_excl_z_in, const T* energy_excl_w_in,
                    const T* force_excl_x_in, const T* force_excl_y_in, const T* force_excl_z_in,
                    const T* force_excl_w_in) :
    theme{theme_in}, basis{basis_in}, lookup{lookup_in}, index_bound{index_bound_in},
    excl_offset{excl_offset_in}, index_shift_bits{index_shift_bits_in},
    dp_detail_mask{dp_detail_mask_in}, sp_detail_mask{sp_detail_mask_in},
    arg_offset{arg_offset_in}, energy_x{energy_x_in}, energy_y{energy_y_in}, energy_z{energy_z_in},
    energy_w{energy_w_in}, force_x{force_x_in}, force_y{force_y_in}, force_z{force_z_in},
    force_w{force_w_in}, energy_excl_x{energy_excl_x_in}, energy_excl_y{energy_excl_y_in},
    energy_excl_z{energy_excl_z_in}, energy_excl_w{energy_excl_w_in},
    force_excl_x{force_excl_x_in}, force_excl_y{force_excl_y_in}, force_excl_z{force_excl_z_in},
    force_excl_w{force_excl_w_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T4>
PPITable::PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b,
                   const LogScaleSpline<T4> &spl_c, const LogScaleSpline<T4> &spl_d,
                   const double cutoff_in) :
    theme{findTheme(spl_a)},
    basis_set{spl_a.getBasisSet()},
    indexing_method{spl_a.getIndexingMethod()},
    cutoff{cutoff_in},
    max_range{spl_a.getMaximumRange()},
    min_range{spl_a.getMinimumRange()},
    argument_offset{spl_a.getIndexingOffset()},
    dsum_tol{recoverDirectSumTolerance(spl_a.getEwaldCoefficient(), cutoff_in)},
    ew_coeff{spl_a.getEwaldCoefficient()},
    mantissa_bits{spl_a.getBitStride()},
    coulomb{spl_a.getCoulombConstant()},
    dp_exclusion_offset{0}, sp_exclusion_offset{0},
    dp_detail_bitmask{doublePrecisionSplineDetailMask(mantissa_bits)},
    sp_detail_bitmask{singlePrecisionSplineDetailMask(mantissa_bits)},
    energy{HybridKind::POINTER, "etab_energy"},
    force{HybridKind::POINTER, "etab_force"},
    energy_with_exclusions{HybridKind::POINTER, "etab_energy_excl"},
    force_with_exclusions{HybridKind::POINTER, "etab_force_excl"},
    sp_energy{HybridKind::POINTER, "sp_etab_energy"},
    sp_force{HybridKind::POINTER, "sp_etab_force"},
    sp_energy_with_exclusions{HybridKind::POINTER, "sp_etab_energy_excl"},
    sp_force_with_exclusions{HybridKind::POINTER, "sp_etab_force_excl"},
    coeffs{HybridKind::ARRAY, "etab_coeffs"},
    sp_coeffs{HybridKind::ARRAY, "etab_sp_coeffs"},
    energy_x{HybridKind::POINTER, "etab_energy_x"},
    energy_y{HybridKind::POINTER, "etab_energy_y"},
    energy_z{HybridKind::POINTER, "etab_energy_z"},
    energy_w{HybridKind::POINTER, "etab_energy_w"},
    force_x{HybridKind::POINTER, "etab_force_x"},
    force_y{HybridKind::POINTER, "etab_force_y"},
    force_z{HybridKind::POINTER, "etab_force_z"},
    force_w{HybridKind::POINTER, "etab_force_w"},
    energy_with_excl_x{HybridKind::POINTER, "etab_energy_excl_x"},
    energy_with_excl_y{HybridKind::POINTER, "etab_energy_excl_y"},
    energy_with_excl_z{HybridKind::POINTER, "etab_energy_excl_z"},
    energy_with_excl_w{HybridKind::POINTER, "etab_energy_excl_w"},
    force_with_excl_x{HybridKind::POINTER, "etab_force_excl_x"},
    force_with_excl_y{HybridKind::POINTER, "etab_force_excl_y"},
    force_with_excl_z{HybridKind::POINTER, "etab_force_excl_z"},
    force_with_excl_w{HybridKind::POINTER, "etab_force_excl_w"},
    sp_energy_x{HybridKind::POINTER, "sp_etab_energy_x"},
    sp_energy_y{HybridKind::POINTER, "sp_etab_energy_y"},
    sp_energy_z{HybridKind::POINTER, "sp_etab_energy_z"},
    sp_energy_w{HybridKind::POINTER, "sp_etab_energy_w"},
    sp_force_x{HybridKind::POINTER, "sp_etab_force_x"},
    sp_force_y{HybridKind::POINTER, "sp_etab_force_y"},
    sp_force_z{HybridKind::POINTER, "sp_etab_force_z"},
    sp_force_w{HybridKind::POINTER, "sp_etab_force_w"},
    sp_energy_with_excl_x{HybridKind::POINTER, "sp_etab_energy_excl_x"},
    sp_energy_with_excl_y{HybridKind::POINTER, "sp_etab_energy_excl_y"},
    sp_energy_with_excl_z{HybridKind::POINTER, "sp_etab_energy_excl_z"},
    sp_energy_with_excl_w{HybridKind::POINTER, "sp_etab_energy_excl_w"},
    sp_force_with_excl_x{HybridKind::POINTER, "sp_etab_force_excl_x"},
    sp_force_with_excl_y{HybridKind::POINTER, "sp_etab_force_excl_y"},
    sp_force_with_excl_z{HybridKind::POINTER, "sp_etab_force_excl_z"},
    sp_force_with_excl_w{HybridKind::POINTER, "sp_etab_force_excl_w"},
    elemental_coeffs{HybridKind::ARRAY, "etab_ele_coeffs"},
    sp_elemental_coeffs{HybridKind::ARRAY, "etab_ele_cp_coeffs"}
{
  checkSplineCompatibility(spl_a, spl_b);
  checkSplineCompatibility(spl_a, spl_c);
  checkSplineCompatibility(spl_a, spl_d);
  
  // Determine which potential is which and populate the appropriate coefficients array.
  const LogScaleSpline<T4> &u_ptr = findNonExclPotential(spl_a, spl_b, spl_c, spl_d);
  const LogScaleSpline<T4> &du_ptr = findNonExclForce(spl_a, spl_b, spl_c, spl_d);
  const LogScaleSpline<T4> &u_excl_ptr = findExclPotential(spl_a, spl_b, spl_c, spl_d);
  const LogScaleSpline<T4> &du_excl_ptr = findExclForce(spl_a, spl_b, spl_c, spl_d);
  populateCoefficients<T4>(u_ptr, du_ptr, u_excl_ptr, du_excl_ptr);
  
  // Populate the remaining coefficients array.
  const size_t ct = std::type_index(typeid(T4)).hash_code();
  if (ct == double4_type_index) {
    const std::vector<LogScaleSpline<float4>> splv = buildAllSplineTables<float4>();
    populateCoefficients<float4>(splv[0], splv[1], splv[2], splv[3]);
  }
  else {
    const std::vector<LogScaleSpline<double4>> splv = buildAllSplineTables<double4>();
    populateCoefficients<double4>(splv[0], splv[1], splv[2], splv[3]);
  }

  // Upload the contents immediately if possible
#ifdef STORMM_USE_HPC
  this->upload();
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
PPITable::PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b,
                   const LogScaleSpline<T4> &spl_c) :
    PPITable(spl_a, spl_b, spl_c, getTablePriority(spl_a, spl_b, spl_c))
{}

//-------------------------------------------------------------------------------------------------
template <typename T4>
PPITable::PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b) :
    PPITable(spl_a, spl_b, getTablePriority(spl_a, spl_b))
{}

//-------------------------------------------------------------------------------------------------
template <typename T4>
PPITable::PPITable(const LogScaleSpline<T4> &spl_a) :
    PPITable(spl_a, getTablePriority(spl_a))
{}

//-------------------------------------------------------------------------------------------------
template <typename T4> NonbondedTheme PPITable::findTheme(const LogScaleSpline<T4> &spl) const {
  switch (spl.getForm()) {
  case LogSplineForm::ELEC_PME_DIRECT:
  case LogSplineForm::ELEC_PME_DIRECT_EXCL:
  case LogSplineForm::DELEC_PME_DIRECT:
  case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    return NonbondedTheme::ELECTROSTATIC;
  case LogSplineForm::CUSTOM:
    return NonbondedTheme::ALL;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4> std::vector<LogScaleSpline<T4>> PPITable::buildAllSplineTables() const {
  std::vector<LogScaleSpline<T4>> result;
  result.reserve(4);
  switch (theme) {
  case NonbondedTheme::ELECTROSTATIC:
    {
      const std::vector<LogSplineForm> all_forms = { LogSplineForm::ELEC_PME_DIRECT,
                                                     LogSplineForm::DELEC_PME_DIRECT,
                                                     LogSplineForm::ELEC_PME_DIRECT_EXCL,
                                                     LogSplineForm::DELEC_PME_DIRECT_EXCL };
      for (size_t i = 0; i < all_forms.size(); i++) {
        result.emplace_back(all_forms[i], ew_coeff, coulomb, mantissa_bits, max_range, min_range,
                            indexing_method, basis_set, 2, argument_offset,
                            ExceptionResponse::DIE);
      } 
    }
    break;
  case NonbondedTheme::VAN_DER_WAALS:
  case NonbondedTheme::ALL:
    break;
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
template <typename T4>
const LogScaleSpline<T4>& PPITable::findNonExclPotential(const LogScaleSpline<T4> &spl_a,
                                                         const LogScaleSpline<T4> &spl_b,
                                                         const LogScaleSpline<T4> &spl_c,
                                                         const LogScaleSpline<T4> &spl_d) const {
  const std::vector<const LogScaleSpline<T4>*> lss = { spl_a.getSelfPointer(),
                                                       spl_b.getSelfPointer(),
                                                       spl_c.getSelfPointer(),
                                                       spl_d.getSelfPointer() };
  for (int i = 0; i < 4; i++) {
    switch (lss[i]->getForm()) {
    case LogSplineForm::ELEC_PME_DIRECT:
      return *lss[i];
    case LogSplineForm::ELEC_PME_DIRECT_EXCL:
    case LogSplineForm::DELEC_PME_DIRECT:
    case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    case LogSplineForm::CUSTOM:
      break;
    }
  }
  rtErr("No potential function was provided.", "PPITable", "findNonExclPotential");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
const LogScaleSpline<T4>& PPITable::findExclPotential(const LogScaleSpline<T4> &spl_a,
                                                      const LogScaleSpline<T4> &spl_b,
                                                      const LogScaleSpline<T4> &spl_c,
                                                      const LogScaleSpline<T4> &spl_d) const {
  const std::vector<const LogScaleSpline<T4>*> lss = { spl_a.getSelfPointer(),
                                                       spl_b.getSelfPointer(),
                                                       spl_c.getSelfPointer(),
                                                       spl_d.getSelfPointer() };
  for (int i = 0; i < 4; i++) {
    switch (lss[i]->getForm()) {
    case LogSplineForm::ELEC_PME_DIRECT_EXCL:
      return *lss[i];
    case LogSplineForm::ELEC_PME_DIRECT:
    case LogSplineForm::DELEC_PME_DIRECT:
    case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    case LogSplineForm::CUSTOM:
      break;
    }
  }
  rtErr("No potential function (with exclusions) was provided.", "PPITable",
        "findNonExclPotential");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
const LogScaleSpline<T4>& PPITable::findNonExclForce(const LogScaleSpline<T4> &spl_a,
                                                     const LogScaleSpline<T4> &spl_b,
                                                     const LogScaleSpline<T4> &spl_c,
                                                     const LogScaleSpline<T4> &spl_d) const {
  const std::vector<const LogScaleSpline<T4>*> lss = { spl_a.getSelfPointer(),
                                                       spl_b.getSelfPointer(),
                                                       spl_c.getSelfPointer(),
                                                       spl_d.getSelfPointer() };
  for (int i = 0; i < 4; i++) {
    switch (lss[i]->getForm()) {
    case LogSplineForm::DELEC_PME_DIRECT:
      return *lss[i];
    case LogSplineForm::ELEC_PME_DIRECT:
    case LogSplineForm::ELEC_PME_DIRECT_EXCL:
    case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    case LogSplineForm::CUSTOM:
      break;
    }
  }
  rtErr("No force function was provided.", "PPITable", "findNonExclPotential");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
const LogScaleSpline<T4>& PPITable::findExclForce(const LogScaleSpline<T4> &spl_a,
                                                  const LogScaleSpline<T4> &spl_b,
                                                  const LogScaleSpline<T4> &spl_c,
                                                  const LogScaleSpline<T4> &spl_d) const {
  const std::vector<const LogScaleSpline<T4>*> lss = { spl_a.getSelfPointer(),
                                                       spl_b.getSelfPointer(),
                                                       spl_c.getSelfPointer(),
                                                       spl_d.getSelfPointer() };
  for (int i = 0; i < 4; i++) {
    switch (lss[i]->getForm()) {
    case LogSplineForm::DELEC_PME_DIRECT_EXCL:
      return *lss[i];
    case LogSplineForm::ELEC_PME_DIRECT:
    case LogSplineForm::ELEC_PME_DIRECT_EXCL:
    case LogSplineForm::DELEC_PME_DIRECT:
    case LogSplineForm::CUSTOM:
      break;
    }
  }
  rtErr("No force function (with exclusions) was provided.", "PPITable", "findNonExclPotential");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4> uint PPITable::checkPriority(const LogScaleSpline<T4> &spl_x) const {
  switch (spl_x.getForm()) {
  case LogSplineForm::ELEC_PME_DIRECT:
    return 0x1;
  case LogSplineForm::ELEC_PME_DIRECT_EXCL:
    return 0x2;
  case LogSplineForm::DELEC_PME_DIRECT:
    return 0x4;
  case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    return 0x8;
  case LogSplineForm::CUSTOM:
    rtErr("Custom potentials cannot be compiled into this object.", "PPITable", "checkPriority");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
void PPITable::checkSplineCompatibility(const LogScaleSpline<T4> &spl_a,
                                        const LogScaleSpline<T4> &spl_b) const {
  if (spl_a.getIndexingMethod() != spl_b.getIndexingMethod()) {
    rtErr("Spline tables with different indexing methods (" +
          getEnumerationName(spl_a.getIndexingMethod()) + ", " +
          getEnumerationName(spl_b.getIndexingMethod()) + ") cannot be compiled together.",
          "PPITable", "checkSplineCompatibility");
  }
  if (spl_a.getBasisSet() != spl_b.getBasisSet()) {
    rtErr("Spline tables with different basis sets (" +
          getEnumerationName(spl_a.getBasisSet()) + ", " +
          getEnumerationName(spl_b.getBasisSet()) + ") cannot be compiled together.",
          "PPITable", "checkSplineCompatibility");
  }
  if (spl_a.getBitStride() != spl_b.getBitStride()) {
    rtErr("Spline tables with different bit strides (" + std::to_string(spl_a.getBitStride()) +
          ", " + std::to_string(spl_b.getBitStride()) + ") cannot be compiled together.",
          "PPITable", "checkSplineCompatibility");
  }
  if (fabs(spl_a.getEwaldCoefficient() - spl_b.getEwaldCoefficient()) > constants::tiny) {
    rtErr("Spline tables with different Gaussian particle widths (" +
          realToString(0.5 / spl_a.getEwaldCoefficient(), 9, 6, NumberFormat::STANDARD_REAL) +
          ", " + 
          realToString(0.5 / spl_b.getEwaldCoefficient(), 9, 6, NumberFormat::STANDARD_REAL) +
          ") cannot be compiled together.", "PPITable", "checkSplineCompatibility");
  }
  if (fabs(spl_a.getCoulombConstant() - spl_b.getCoulombConstant()) > constants::tiny) {
    rtErr("Spline tables with different variations of Coulomb's constant (" +
          realToString(spl_a.getCoulombConstant(), 10, 6, NumberFormat::STANDARD_REAL) + ", " + 
          realToString(spl_b.getCoulombConstant(), 10, 6, NumberFormat::STANDARD_REAL) +
          ") cannot be compiled together.", "PPITable", "checkSplineCompatibility");
  }
  if (fabs(spl_a.getMaximumRange() - spl_b.getMaximumRange()) > constants::tiny) {
    rtErr("Spline tables with different ranges (" + 
          realToString(spl_a.getMaximumRange(), 10, 6, NumberFormat::STANDARD_REAL) + ", " + 
          realToString(spl_b.getMaximumRange(), 10, 6, NumberFormat::STANDARD_REAL) +
          ") cannot be compiled together.", "PPITable", "checkSplineCompatibility");
  }
  if (spl_a.getIndexingOffset() != spl_b.getIndexingOffset()) {
    rtErr("Spline tables with different indexing offsets (" +
          std::to_string(spl_a.getIndexingOffset()) + ", " +
          std::to_string(spl_a.getIndexingOffset()) + ") cannot be compiled together.",
          "PPITable", "checkSplineCompatibility");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
LogScaleSpline<T4> PPITable::getTablePriority(const LogScaleSpline<T4> &spl_a,
                                              const LogScaleSpline<T4> &spl_b,
                                              const LogScaleSpline<T4> &spl_c) const {
  checkSplineCompatibility(spl_a, spl_b);
  checkSplineCompatibility(spl_b, spl_c);
  const uint holdings = checkPriority(spl_a) | checkPriority(spl_b) | checkPriority(spl_c);
  const LogSplineForm missing_form = findMissingForm(holdings);
  return LogScaleSpline<T4>(missing_form, spl_a.getEwaldCoefficient(), spl_a.getCoulombConstant(),
                            spl_a.getBitStride(), spl_a.getMaximumRange(), spl_a.getMinimumRange(),
                            spl_a.getIndexingMethod(), spl_a.getBasisSet(),
                            spl_a.getOptimizationDepth(), spl_a.getIndexingOffset());
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
LogScaleSpline<T4> PPITable::getTablePriority(const LogScaleSpline<T4> &spl_a,
                                              const LogScaleSpline<T4> &spl_b) const {
  checkSplineCompatibility(spl_a, spl_b);
  const uint holdings = checkPriority(spl_a) | checkPriority(spl_b);
  const LogSplineForm missing_form = findMissingForm(holdings);
  return LogScaleSpline<T4>(missing_form, spl_a.getEwaldCoefficient(), spl_a.getCoulombConstant(),
                            spl_a.getBitStride(), spl_a.getMaximumRange(), spl_a.getMinimumRange(),
                            spl_a.getIndexingMethod(), spl_a.getBasisSet(),
                            spl_a.getOptimizationDepth(), spl_a.getIndexingOffset());
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
LogScaleSpline<T4> PPITable::getTablePriority(const LogScaleSpline<T4> &spl_a) const {
  const uint holdings = checkPriority(spl_a);
  const LogSplineForm missing_form = findMissingForm(holdings);
  return LogScaleSpline<T4>(missing_form, spl_a.getEwaldCoefficient(), spl_a.getCoulombConstant(),
                            spl_a.getBitStride(), spl_a.getMaximumRange(), spl_a.getMinimumRange(),
                            spl_a.getIndexingMethod(), spl_a.getBasisSet(),
                            spl_a.getOptimizationDepth(), spl_a.getIndexingOffset());
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
void PPITable::populateCoefficients(const LogScaleSpline<T4> &u, const LogScaleSpline<T4> &du,
                                    const LogScaleSpline<T4> &ux, const LogScaleSpline<T4> &dux) {
  const LogSplineTable tbl_u   = u.data();
  const LogSplineTable tbl_ux  = ux.data();
  const LogSplineTable tbl_du  = du.data();
  const LogSplineTable tbl_dux = dux.data();
  
  // Allocate a common space for all four potentials.
  const size_t ct = std::type_index(typeid(T4)).hash_code();
  const size_t tbl_len = u.getSplineIndexByRealArg(u.getMaximumRange());
  const double max_imaginary_dist = sqrt(3.0) * ((5.0 * u.getMaximumRange()) +
                                                 static_cast<double>(warp_size_int));
  const size_t catch_region = pow(2.0, u.getBitStride()) * ceil(log2(max_imaginary_dist));
  if (ct == double4_type_index) {
    dp_exclusion_offset = 2 * tbl_len;
    coeffs.resize((4 * tbl_len) + catch_region);
    elemental_coeffs.resize((16 * tbl_len) + catch_region);
    double4* coef_ptr = coeffs.data();
    double* ele_coef_ptr = elemental_coeffs.data();
    for (int i = 0; i < tbl_len; i++) {
      const double4 tmp_ucoef = { static_cast<double>(tbl_u.table[i].x),
                                  static_cast<double>(tbl_u.table[i].y),
                                  static_cast<double>(tbl_u.table[i].z),
                                  static_cast<double>(tbl_u.table[i].w) };
      coef_ptr[                i] = tmp_ucoef;
      const double4 tmp_ducoef = { static_cast<double>(tbl_du.table[i].x),
                                   static_cast<double>(tbl_du.table[i].y),
                                   static_cast<double>(tbl_du.table[i].z),
                                   static_cast<double>(tbl_du.table[i].w) };
      coef_ptr[     tbl_len  + i] = tmp_ducoef;
      const double4 tmp_uxcoef = { static_cast<double>(tbl_ux.table[i].x),
                                   static_cast<double>(tbl_ux.table[i].y),
                                   static_cast<double>(tbl_ux.table[i].z),
                                   static_cast<double>(tbl_ux.table[i].w) };
      coef_ptr[(2 * tbl_len) + i] = tmp_uxcoef;
      const double4 tmp_duxcoef = { static_cast<double>(tbl_dux.table[i].x),
                                    static_cast<double>(tbl_dux.table[i].y),
                                    static_cast<double>(tbl_dux.table[i].z),
                                    static_cast<double>(tbl_dux.table[i].w) };
      coef_ptr[(3 * tbl_len) + i] = tmp_duxcoef;
      ele_coef_ptr[                 i] = tmp_ucoef.x;
      ele_coef_ptr[      tbl_len  + i] = tmp_ucoef.y;
      ele_coef_ptr[( 2 * tbl_len) + i] = tmp_ucoef.z;
      ele_coef_ptr[( 3 * tbl_len) + i] = tmp_ucoef.w;
      ele_coef_ptr[( 4 * tbl_len) + i] = tmp_ducoef.x;
      ele_coef_ptr[( 5 * tbl_len) + i] = tmp_ducoef.y;
      ele_coef_ptr[( 6 * tbl_len) + i] = tmp_ducoef.z;
      ele_coef_ptr[( 7 * tbl_len) + i] = tmp_ducoef.w;
      ele_coef_ptr[( 8 * tbl_len) + i] = tmp_uxcoef.x;
      ele_coef_ptr[( 9 * tbl_len) + i] = tmp_uxcoef.y;
      ele_coef_ptr[(10 * tbl_len) + i] = tmp_uxcoef.z;
      ele_coef_ptr[(11 * tbl_len) + i] = tmp_uxcoef.w;
      ele_coef_ptr[(12 * tbl_len) + i] = tmp_duxcoef.x;
      ele_coef_ptr[(13 * tbl_len) + i] = tmp_duxcoef.y;
      ele_coef_ptr[(14 * tbl_len) + i] = tmp_duxcoef.z;
      ele_coef_ptr[(15 * tbl_len) + i] = tmp_duxcoef.w;    
    }

    // Set POINTER-kind Hybrid objects for tuple data
    energy.setPointer(&coeffs,                           0, tbl_len);
    force.setPointer(&coeffs,                      tbl_len, tbl_len);
    energy_with_exclusions.setPointer(&coeffs, 2 * tbl_len, tbl_len);
    force_with_exclusions.setPointer(&coeffs,  3 * tbl_len, tbl_len);

    // Set POINTER-kind Hybrid objects for elemental data
    energy_x.setPointer(&elemental_coeffs,                      0, tbl_len);
    energy_y.setPointer(&elemental_coeffs,                      1, tbl_len);
    energy_z.setPointer(&elemental_coeffs,                      2, tbl_len);
    energy_w.setPointer(&elemental_coeffs,                      3, tbl_len);
    force_x.setPointer(&elemental_coeffs,             4 * tbl_len, tbl_len);
    force_y.setPointer(&elemental_coeffs,             5 * tbl_len, tbl_len);
    force_z.setPointer(&elemental_coeffs,             6 * tbl_len, tbl_len);
    force_w.setPointer(&elemental_coeffs,             7 * tbl_len, tbl_len);
    energy_with_excl_x.setPointer(&elemental_coeffs,  8 * tbl_len, tbl_len);
    energy_with_excl_y.setPointer(&elemental_coeffs,  9 * tbl_len, tbl_len);
    energy_with_excl_z.setPointer(&elemental_coeffs, 10 * tbl_len, tbl_len);
    energy_with_excl_w.setPointer(&elemental_coeffs, 11 * tbl_len, tbl_len);
    force_with_excl_x.setPointer(&elemental_coeffs,  12 * tbl_len, tbl_len);
    force_with_excl_y.setPointer(&elemental_coeffs,  13 * tbl_len, tbl_len);
    force_with_excl_z.setPointer(&elemental_coeffs,  14 * tbl_len, tbl_len);
    force_with_excl_w.setPointer(&elemental_coeffs,  15 * tbl_len, tbl_len);
  }
  else if (ct == float4_type_index) {
    sp_exclusion_offset = 2 * tbl_len;
    sp_coeffs.resize((4 * tbl_len) + catch_region);
    sp_elemental_coeffs.resize((16 * tbl_len) + catch_region);
    float4* sp_coef_ptr = sp_coeffs.data();
    float* sp_ele_coef_ptr = sp_elemental_coeffs.data();
    for (int i = 0; i < tbl_len; i++) {
      const float4 tmp_ucoef = { static_cast<float>(tbl_u.table[i].x),
                                 static_cast<float>(tbl_u.table[i].y),
                                 static_cast<float>(tbl_u.table[i].z),
                                 static_cast<float>(tbl_u.table[i].w) };
      sp_coef_ptr[                i] = tmp_ucoef;
      const float4 tmp_ducoef = { static_cast<float>(tbl_du.table[i].x),
                                  static_cast<float>(tbl_du.table[i].y),
                                  static_cast<float>(tbl_du.table[i].z),
                                  static_cast<float>(tbl_du.table[i].w) };
      sp_coef_ptr[     tbl_len  + i] = tmp_ducoef;
      const float4 tmp_uxcoef = { static_cast<float>(tbl_ux.table[i].x),
                                  static_cast<float>(tbl_ux.table[i].y),
                                  static_cast<float>(tbl_ux.table[i].z),
                                  static_cast<float>(tbl_ux.table[i].w) };
      sp_coef_ptr[(2 * tbl_len) + i] = tmp_uxcoef;
      const float4 tmp_duxcoef = { static_cast<float>(tbl_dux.table[i].x),
                                   static_cast<float>(tbl_dux.table[i].y),
                                   static_cast<float>(tbl_dux.table[i].z),
                                   static_cast<float>(tbl_dux.table[i].w) };
      sp_coef_ptr[(3 * tbl_len) + i] = tmp_duxcoef;
      sp_ele_coef_ptr[                 i] = tmp_ucoef.x;
      sp_ele_coef_ptr[      tbl_len  + i] = tmp_ucoef.y;
      sp_ele_coef_ptr[( 2 * tbl_len) + i] = tmp_ucoef.z;
      sp_ele_coef_ptr[( 3 * tbl_len) + i] = tmp_ucoef.w;
      sp_ele_coef_ptr[( 4 * tbl_len) + i] = tmp_ducoef.x;
      sp_ele_coef_ptr[( 5 * tbl_len) + i] = tmp_ducoef.y;
      sp_ele_coef_ptr[( 6 * tbl_len) + i] = tmp_ducoef.z;
      sp_ele_coef_ptr[( 7 * tbl_len) + i] = tmp_ducoef.w;
      sp_ele_coef_ptr[( 8 * tbl_len) + i] = tmp_uxcoef.x;
      sp_ele_coef_ptr[( 9 * tbl_len) + i] = tmp_uxcoef.y;
      sp_ele_coef_ptr[(10 * tbl_len) + i] = tmp_uxcoef.z;
      sp_ele_coef_ptr[(11 * tbl_len) + i] = tmp_uxcoef.w;
      sp_ele_coef_ptr[(12 * tbl_len) + i] = tmp_duxcoef.x;
      sp_ele_coef_ptr[(13 * tbl_len) + i] = tmp_duxcoef.y;
      sp_ele_coef_ptr[(14 * tbl_len) + i] = tmp_duxcoef.z;
      sp_ele_coef_ptr[(15 * tbl_len) + i] = tmp_duxcoef.w;
    }

    // Set POINTER-kind Hybrid objects for tuple data
    sp_energy.setPointer(&sp_coeffs,                           0, tbl_len);
    sp_force.setPointer(&sp_coeffs,                      tbl_len, tbl_len);
    sp_energy_with_exclusions.setPointer(&sp_coeffs, 2 * tbl_len, tbl_len);
    sp_force_with_exclusions.setPointer(&sp_coeffs,  3 * tbl_len, tbl_len);

    // Set POINTER-kind Hybrid objects for elemental data
    sp_energy_x.setPointer(&sp_elemental_coeffs,                      0, tbl_len);
    sp_energy_y.setPointer(&sp_elemental_coeffs,                      1, tbl_len);
    sp_energy_z.setPointer(&sp_elemental_coeffs,                      2, tbl_len);
    sp_energy_w.setPointer(&sp_elemental_coeffs,                      3, tbl_len);
    sp_force_x.setPointer(&sp_elemental_coeffs,             4 * tbl_len, tbl_len);
    sp_force_y.setPointer(&sp_elemental_coeffs,             5 * tbl_len, tbl_len);
    sp_force_z.setPointer(&sp_elemental_coeffs,             6 * tbl_len, tbl_len);
    sp_force_w.setPointer(&sp_elemental_coeffs,             7 * tbl_len, tbl_len);
    sp_energy_with_excl_x.setPointer(&sp_elemental_coeffs,  8 * tbl_len, tbl_len);
    sp_energy_with_excl_y.setPointer(&sp_elemental_coeffs,  9 * tbl_len, tbl_len);
    sp_energy_with_excl_z.setPointer(&sp_elemental_coeffs, 10 * tbl_len, tbl_len);
    sp_energy_with_excl_w.setPointer(&sp_elemental_coeffs, 11 * tbl_len, tbl_len);
    sp_force_with_excl_x.setPointer(&sp_elemental_coeffs,  12 * tbl_len, tbl_len);
    sp_force_with_excl_y.setPointer(&sp_elemental_coeffs,  13 * tbl_len, tbl_len);
    sp_force_with_excl_z.setPointer(&sp_elemental_coeffs,  14 * tbl_len, tbl_len);
    sp_force_with_excl_w.setPointer(&sp_elemental_coeffs,  15 * tbl_len, tbl_len);
  }
}

} // namespace energy
} // namespace stormm
