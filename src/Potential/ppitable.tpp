// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename T4>
PPITable<T4>::PPITable(const NonbondedTheme theme_in, const BasisFunctions basis_set_in,
                       const TableIndexing indexing_method_in, const double cutoff_in,
                       const double dsum_tol_in, const int mantissa_bits_in,
                       const double coulomb_in, const double min_spl_compute_range_in,
                       const double min_offset_in) :
    theme{theme_in}, basis_set{basis_set_in}, indexing_method{indexing_method_in},
    cutoff{cutoff_in}, max_range{cutoff_in}, dsum_tol{dsum_tol_in},
    ew_coeff{ewaldCoefficient(cutoff_in, dsum_tol_in)}, mantissa_bits{mantissa_bits_in},
    coulomb{coulomb_in}, exclusion_offset{0},
    energy{HybridKind::POINTER, "etab_energy"},
    force{HybridKind::POINTER, "etab_force"},
    energy_with_exclusions{HybridKind::POINTER, "etab_energy_excl"},
    force_with_exclusions{HybridKind::POINTER, "etab_force_excl"},
    sp_energy{HybridKind::POINTER, "etab_energy"},
    sp_force{HybridKind::POINTER, "etab_force"},
    sp_energy_with_exclusions{HybridKind::POINTER, "etab_energy_excl"},
    sp_force_with_exclusions{HybridKind::POINTER, "etab_force_excl"},
    coeffs{HybridKind::ARRAY, "etab_coeffs"},
    sp_coeffs{HybridKind::ARRAY, "etab_sp_coeffs"}
{
  // Create the four spline tables, then feed into the same templated functions for re-arranging
  // them as the constructors taking spline table inputs.
  max_range = exp2(ceil(log2(cutoff)));
  double min_spl_compute_range = min_spl_compute_range_in;
  double min_offset = min_offset_in;
  switch (indexing_method) {
  case TableIndexing::ARG:
  case TableIndexing::ARG_OFFSET:
    break;
  case TableIndexing::SQUARED_ARG:
  case TableIndexing::SQ_ARG_OFFSET:
    min_offset *= min_offset;
    min_spl_compute_range *= min_spl_compute_range;
    break;
  }
  std::vector<LogScaleSpline<T4>> splv;
  splv.reserve(4);
  switch (theme) {
  case NonbondedTheme::ELECTROSTATIC:
    {
      const std::vector<LogSplineForm> all_forms = { LogSplineForm::ELEC_PME_DIRECT,
                                                     LogSplineForm::DELEC_PME_DIRECT,
                                                     LogSplineForm::ELEC_PME_DIRECT_EXCL,
                                                     LogSplineForm::DELEC_PME_DIRECT_EXCL };
      for (size_t i = 0; i < all_forms.size(); i++) {
        splv.emplace_back(all_forms[i], ew_coeff, coulomb, mantissa_bits, max_range,
                          min_spl_compute_range, indexing_method, basis_set, 2, min_offset,
                          ExceptionResponse::DIE);
      }
    }
  case NonbondedTheme::VAN_DER_WAALS:
  case NonbondedTheme::ALL:
    break;
  }
  populateCoefficients(findNonExclPotential(spl_a, spl_b, spl_c, spl_d),
                       findNonExclForce(spl_a, spl_b, spl_c, spl_d),
                       findExclPotential(spl_a, spl_b, spl_c, spl_d),
                       findExclForce(spl_a, spl_b, spl_c, spl_d));
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
PPITable<T4>::PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b,
                       const LogScaleSpline<T4> &spl_c, const LogScaleSpline<T4> &spl_d,
                       const double cutoff_in) :
    theme{NonbondedTheme::ELECTROSTATIC},
    precision{(std::type_index(typeid).hash_code() == double_type_index) ?
              PrecisionModel::DOUBLE : PrecisionModel::SINGLE},
    basis_set{spl_a.getBasisSet()},
    indexing_method{spl_a.getIndexingMethod()},
    cutoff{cutoff_in},
    dsum_tol{recoverDirectSumTolerance(spl_a.getEwaldCoefficient(), cutoff_in)},
    mantissa_bits{spl_a.getBitStride()},
    coulomb{spl_a.getCoulombConstant()}
    exclusion_offset{0},
    energy{HybridKind::POINTER, "ppi_u"},
    force{HybridKind::POINTER, "ppi_du"},
    energy_with_exclusions{HybridKind::POINTER, "ppi_ux"},
    force_with_exclusions{HybridKind::POINTER, "ppi_dux"},
    coeffs{HybridKind::ARRAY, "ppi_coeffs"}
{
  checkSplineCompatibility(spl_a, spl_b);
  checkSplineCompatibility(spl_a, spl_c);
  checkSplineCompatibility(spl_a, spl_d);

  // Determine which potential is which and populate the coefficients array.
  populateCoefficients(findNonExclPotential(spl_a, spl_b, spl_c, spl_d),
                       findNonExclForce(spl_a, spl_b, spl_c, spl_d),
                       findExclPotential(spl_a, spl_b, spl_c, spl_d),
                       findExclForce(spl_a, spl_b, spl_c, spl_d));
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
PPITable<T4>::PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b,
                       const LogScaleSpline<T4> &spl_c) :
    PPITable(spl_a, spl_b, spl_c, getTablePriority(spl_a, spl_b, spl_c))
{}

//-------------------------------------------------------------------------------------------------
template <typename T4>
PPITable<T4>::PPITable(const LogScaleSpline<T4> &spl_a, const LogScaleSpline<T4> &spl_b) :
    PPITable(spl_a, spl_b, getTablePriority(spl_a, spl_b))
{}

//-------------------------------------------------------------------------------------------------
template <typename T4>
PPITable<T4>::PPITable(const LogScaleSpline<T4> &spl_a) :
    PPITable(spl_a, getTablePriority(spl_a))
{}

//-------------------------------------------------------------------------------------------------
template <typename T4> NonbondedTheme PPITable<T4>::getTheme() const {
  return theme;
}

//-------------------------------------------------------------------------------------------------
template <typename T4> double PPITable<T4>::getCutoff() const {
  return cutoff;
}

//-------------------------------------------------------------------------------------------------
template <typename T4> double PPITable<T4>::getMaximumRange() const {
  return max_range;
}

//-------------------------------------------------------------------------------------------------
template <typename T4> double PPITable<T4>::getDirectSumTolerance() const {
  return dsum_tol;
}

//-------------------------------------------------------------------------------------------------
template <typename T4> int PPITable<T4>::getBitStride() const {
  return mantissa_bits;
}

//-------------------------------------------------------------------------------------------------
template <typename T4> double PPITable<T4>::getEwaldCoefficient() const {
  return ew_coeff;
}

//-------------------------------------------------------------------------------------------------
template <typename T4> double PPITable<T4>::getGaussianWidth() const {
  return dsum_tol;
}

//-------------------------------------------------------------------------------------------------
template <typename T4> const PPIKit PPITable<T4>::data(const HybridTargetLevel tier) const {
  return PPIKit(theme, basis_set, indexing_method, exclusion_offset / 2, exclusion_offset,
                
}

//-------------------------------------------------------------------------------------------------
template <typename T4> const LogScaleSpline<T4>&
PPITable<T4>::findNonExclPotential(const LogScaleSpline<T4> &spl_a,
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
const LogScaleSpline<T4>& PPITable<T4>::findExclPotential(const LogScaleSpline<T4> &spl_a,
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
const LogScaleSpline<T4>& PPITable<T4>::findNonExclForce(const LogScaleSpline<T4> &spl_a,
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
const LogScaleSpline<T4>& PPITable<T4>::findExclForce(const LogScaleSpline<T4> &spl_a,
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
template <typename T4> uint PPITable<T4>::checkPriority(const LogScaleSpline<T4> &spl_x) const {
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
void PPITable<T4>::checkSplineCompatibility(const LogScaleSpline<T4> &spl_a,
                                            const LogScaleSpline<T4> &spl_b) {
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
template <typename T4> LogSplineForm PPITable<T4>::findMissingForm(const uint holdings) const {
  if ((holdings & 0x1) == 0) {
    return LogSplineForm::ELEC_PME_DIRECT;
  }
  else if ((holdings & 0x2) == 0) {
    return LogSplineForm::ELEC_PME_DIRECT_EXCL;
  }
  else if ((holdings & 0x4) == 0) {
    return LogSplineForm::DELEC_PME_DIRECT;
  }
  else if ((holdings & 0x8) == 0) {
    return LogSplineForm::DELEC_PME_DIRECT_EXCL;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
LogScaleSpline<T4> PPITable<T4>::getTablePriority(const LogScaleSpline<T4> &spl_a,
                                                  const LogScaleSpline<T4> &spl_b,
                                                  const LogScaleSpline<T4> &spl_c) const {
  checkSplineCompatibility(spl_a, spl_b);
  checkSplineCompatibility(spl_b, spl_c);
  const uint holdings = checkPriority(spl_a) | checkPriority(spl_b) | checkPriority(spl_c);
  const LogSplineForm missing_form = findMissingForm(holdings);
  return LogScaleSpline<T4>(missing_form, spl_a.getEwaldCoefficient(), spl_a.getCoulombConstant(),
                            spl_a.getBitStride(), spl_a.getMaximumRange(), spl_a.getMinimumRange(),
                            spl_a.getIndexingMethod(), spl_a.getBasisSet(),
                            spl_a.getOptimization_depth, spl_a.getIndexingOffset());
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
LogScaleSpline<T4> PPITable<T4>::getTablePriority(const LogScaleSpline<T4> &spl_a,
                                                  const LogScaleSpline<T4> &spl_b) {
  checkSplineCompatibility(spl_a, spl_b);
  const uint holdings = checkPriority(spl_a) | checkPriority(spl_b);
  const LogSplineForm missing_form = findMissingForm(holdings);
  return LogScaleSpline<T4>(missing_form, spl_a.getEwaldCoefficient(), spl_a.getCoulombConstant(),
                            spl_a.getBitStride(), spl_a.getMaximumRange(), spl_a.getMinimumRange(),
                            spl_a.getIndexingMethod(), spl_a.getBasisSet(),
                            spl_a.getOptimization_depth, spl_a.getIndexingOffset());
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
LogScaleSpline<T4> PPITable<T4>::getTablePriority(const LogScaleSpline<T4> &spl_a) {
  const uint holdings = checkPriority(spl_a);
  const LogSplineForm missing_form = findMissingForm(holdings);
  return LogScaleSpline<T4>(missing_form, spl_a.getEwaldCoefficient(), spl_a.getCoulombConstant(),
                            spl_a.getBitStride(), spl_a.getMaximumRange(), spl_a.getMinimumRange(),
                            spl_a.getIndexingMethod(), spl_a.getBasisSet(),
                            spl_a.getOptimization_depth, spl_a.getIndexingOffset());
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
void PPITable<T4>:: populateCoefficients(const LogScaleSpline<T4> &u,
                                         const LogScaleSpline<T4> &du,
                                         const LogScaleSpline<T4> &ux,
                                         const LogScaleSpline<T4> &dux) {
  const LogSplineTable tbl_u   = u.data();
  const LogSplineTable tbl_ux  = ux.data();
  const LogSplineTable tbl_du  = du.data();
  const LogSplineTable tbl_dux = dux.data();

  // Allocate a common space for all four potentials.
  const size_t ct = std::type_index(typeid(T4)).hash_code();
  const size_t tbl_len = u.getSplineIndex(u.getMaximumRange());
  coeffs.resize(4 * exclusion_offset);
  T4* coef_ptr = coeffs.data();
  for (int i = 0; i < tbl_len; i++) {
    coef_ptr[                i] = tbl_u.table[i];
    coef_ptr[      tbl_len + i] = tbl_du.table[i];
    coef_ptr[(2 * tbl_len) + i] = tbl_ux.table[i];
    coef_ptr[(3 * tbl_len) + i] = tbl_dux.table[i];
  }
  exclusion_offset = 2 * tbl_len;
  energy.setPointer(&coeffs,                           0, tbl_len);
  force.setPointer(&coeffs,                      tbl_len, tbl_len);
  energy_with_exclusions.setPointer(&coeffs, 2 * tbl_len, tbl_len);
  force_with_exclusions.setPointer(&coeffs,  3 * tbl_len, tbl_len);
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
const PPIKit<T4> restoreType(const PPIKit<void> *rasa) {
  return PPIKit<T4>(rasa->theme, rasa->basis rasa->lookup, rasa->index_bound, rasa->excl_offset,
                    rasa->sp_detail_mask, rasa->dp_detail_mask, rasa->idx_offset,
                    reinterpret_cast<const T4*>(rasa->energy),
                    reinterpret_cast<const T4*>(rasa->force),
                    reinterpret_cast<const T4*>(rasa->energy_excl),
                    reinterpret_cast<const T4*>(rasa->force_excl));
}

//-------------------------------------------------------------------------------------------------
template <typename T4>
const PPIKit<T4> restoreType(const PPIKit<void> &rasa) {
  return PPIKit<T4>(rasa.theme, rasa.basis rasa.lookup, rasa.index_bound, rasa.excl_offset,
                    rasa.sp_detail_mask, rasa.dp_detail_mask, rasa.idx_offset,
                    reinterpret_cast<const T4*>(rasa.energy),
                    reinterpret_cast<const T4*>(rasa.force),
                    reinterpret_cast<const T4*>(rasa.energy_excl),
                    reinterpret_cast<const T4*>(rasa.force_excl));
}

} // namespace energy
} // namespace stormm
