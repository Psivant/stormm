#include "copyright.h"
#include "ppitable.h"

namespace stormm {
namespace energy {

using stmath::doublePrecisionSplineArgument;
using stmath::singlePrecisionSplineArgument;
  
//-------------------------------------------------------------------------------------------------
PPITable::PPITable(const NonbondedTheme theme_in, const BasisFunctions basis_set_in,
                   const TableIndexing indexing_method_in, const double cutoff_in,
                   const double argument_offset_in, const double dsum_tol_in,
                   const int mantissa_bits_in, const double coulomb_in,
                   const double min_range_in) :
    theme{theme_in}, basis_set{basis_set_in}, indexing_method{indexing_method_in},
    cutoff{cutoff_in}, max_range{cutoff_in}, min_range{min_range_in},
    argument_offset{argument_offset_in}, dsum_tol{dsum_tol_in},
    ew_coeff{ewaldCoefficient(cutoff_in, dsum_tol_in)}, mantissa_bits{mantissa_bits_in},
    coulomb{coulomb_in}, dp_exclusion_offset{0}, sp_exclusion_offset{0},
    dp_detail_bitmask{doublePrecisionSplineDetailMask(mantissa_bits)},
    sp_detail_bitmask{singlePrecisionSplineDetailMask(mantissa_bits)},
    energy{HybridKind::POINTER, "etab_energy"},
    force{HybridKind::POINTER, "etab_force"},
    energy_with_exclusions{HybridKind::POINTER, "etab_energy_excl"},
    force_with_exclusions{HybridKind::POINTER, "etab_force_excl"},
    sp_energy{HybridKind::POINTER, "etab_energy"},
    sp_force{HybridKind::POINTER, "etab_force"},
    sp_energy_with_exclusions{HybridKind::POINTER, "etab_energy_excl"},
    sp_force_with_exclusions{HybridKind::POINTER, "etab_force_excl"},
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
  // Create the four spline tables, then feed into the same templated functions for re-arranging
  // them as the constructors taking spline table inputs.
  max_range = exp2(ceil(log2(cutoff)));
  switch (indexing_method) {
  case TableIndexing::ARG:
  case TableIndexing::ARG_OFFSET:
    break;
  case TableIndexing::SQUARED_ARG:
  case TableIndexing::SQ_ARG_OFFSET:
    argument_offset *= argument_offset;
    min_range *= min_range;
    break;
  }
  const std::vector<LogScaleSpline<double4>> splv = buildAllSplineTables<double4>();
  const std::vector<LogScaleSpline<float4>> sp_splv = buildAllSplineTables<float4>();
  populateCoefficients<double4>(splv[0], splv[1], splv[2], splv[3]);
  populateCoefficients<float4>(sp_splv[0], sp_splv[1], sp_splv[2], sp_splv[3]);

  // Upload immediately if there is a GPU involved
#ifdef STORMM_USE_HPC
  this->upload();
#endif
}

//-------------------------------------------------------------------------------------------------
NonbondedTheme PPITable::getTheme() const {
  return theme;
}

//-------------------------------------------------------------------------------------------------
double PPITable::getCutoff() const {
  return cutoff;
}

//-------------------------------------------------------------------------------------------------
double PPITable::getMaximumRange() const {
  return max_range;
}

//-------------------------------------------------------------------------------------------------
double PPITable::getIndexingOffset() const {
  return argument_offset;
}

//-------------------------------------------------------------------------------------------------
double PPITable::getDirectSumTolerance() const {
  return dsum_tol;
}

//-------------------------------------------------------------------------------------------------
int PPITable::getBitStride() const {
  return mantissa_bits;
}

//-------------------------------------------------------------------------------------------------
double PPITable::getEwaldCoefficient() const {
  return ew_coeff;
}

//-------------------------------------------------------------------------------------------------
double PPITable::getGaussianWidth() const {
  return dsum_tol;
}

//-------------------------------------------------------------------------------------------------
double PPITable::evaluate(const double arg, const LogSplineForm kind,
                          const PrecisionModel prec, const bool use_elemental_tables) const {

  // Compute the basic access index.
  int access_idx = getTableIndex(arg, prec);
  
  // Refine the index based on the function type.
  int excl_offset_base;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    excl_offset_base = dp_exclusion_offset;
    break;
  case PrecisionModel::SINGLE:
    excl_offset_base = sp_exclusion_offset;
    break;
  }
  switch (kind) {
  case LogSplineForm::ELEC_PME_DIRECT:
    break;
  case LogSplineForm::DELEC_PME_DIRECT:
    access_idx += excl_offset_base / 2;
    break;
  case LogSplineForm::ELEC_PME_DIRECT_EXCL:
    access_idx += excl_offset_base;
    break;
  case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    access_idx += 3 * excl_offset_base / 2;
    break;
  case LogSplineForm::CUSTOM:
    rtErr("The particle-particle interaction tables do not comprise " + getEnumerationName(kind) +
          " functions.", "PPITable", "evaluate");
  }

  // Access the tables as specified.
  double result;
  if (use_elemental_tables) {
    switch (prec) {
    case PrecisionModel::DOUBLE:
      {
        const double coef_x = elemental_coeffs.readHost(access_idx);
        const double coef_y = elemental_coeffs.readHost(access_idx);
        const double coef_z = elemental_coeffs.readHost(access_idx);
        const double coef_w = elemental_coeffs.readHost(access_idx);
        switch (basis_set) {
        case BasisFunctions::MIXED_FRACTIONS:
          {
            const double inv_arg = 1.0 / arg;
            result = (arg * coef_x) + coef_y + ((coef_z + (coef_w * inv_arg)) * inv_arg);
          }
          break;
        case BasisFunctions::POLYNOMIAL:
          {
            const double spl_arg = doublePrecisionSplineArgument(arg, mantissa_bits,
                                                                 dp_detail_bitmask);
            result = coef_x + ((coef_y + ((coef_z + (coef_w * spl_arg)) * spl_arg)) * spl_arg);
          }
          break;
        }
      }
      break;
    case PrecisionModel::SINGLE:
      {
        const float coef_x = sp_elemental_coeffs.readHost(access_idx);
        const float coef_y = sp_elemental_coeffs.readHost(access_idx);
        const float coef_z = sp_elemental_coeffs.readHost(access_idx);
        const float coef_w = sp_elemental_coeffs.readHost(access_idx);
        switch (basis_set) {
        case BasisFunctions::MIXED_FRACTIONS:
          {
            const float inv_arg = 1.0f / arg;
            result = (arg * coef_x) + coef_y + ((coef_z + (coef_w * inv_arg)) * inv_arg);
          }
          break;
        case BasisFunctions::POLYNOMIAL:
          {
            const float spl_arg = singlePrecisionSplineArgument(arg, mantissa_bits,
                                                                sp_detail_bitmask);
            result = coef_x + ((coef_y + ((coef_z + (coef_w * spl_arg)) * spl_arg)) * spl_arg);
          }
          break;
        }
      }
      break;
    }
  }
  else {
    switch (prec) {
    case PrecisionModel::DOUBLE:
      {
        const double4 coef = coeffs.readHost(access_idx);
        switch (basis_set) {
        case BasisFunctions::MIXED_FRACTIONS:
          {
            const double inv_arg = 1.0 / arg;
            result = (arg * coef.x) + coef.y + ((coef.z + (coef.w * inv_arg)) * inv_arg);
          }
          break;
        case BasisFunctions::POLYNOMIAL:
          {
            const double spl_arg = doublePrecisionSplineArgument(arg, mantissa_bits,
                                                                 dp_detail_bitmask);
            result = coef.x + ((coef.y + ((coef.z + (coef.w * spl_arg)) * spl_arg)) * spl_arg);
          }
          break;
        }
      }
      break;
    case PrecisionModel::SINGLE:
      {
        const float4 coef = sp_coeffs.readHost(access_idx);
        switch (basis_set) {
        case BasisFunctions::MIXED_FRACTIONS:
          {
            const float inv_arg = 1.0f / arg;
            result = (arg * coef.x) + coef.y + ((coef.z + (coef.w * inv_arg)) * inv_arg);
          }
          break;
        case BasisFunctions::POLYNOMIAL:
          {
            const float spl_arg = singlePrecisionSplineArgument(arg, mantissa_bits,
                                                                sp_detail_bitmask);
            result = coef.x + ((coef.y + ((coef.z + (coef.w * spl_arg)) * spl_arg)) * spl_arg);
          }
          break;
        }
      }
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double PPITable::evaluateByRealArg(const double arg, const LogSplineForm kind,
                                   const PrecisionModel prec,
                                   const bool use_elemental_tables) const {
  switch (indexing_method) {
  case TableIndexing::ARG:
    return evaluate(arg, kind, prec, use_elemental_tables);
  case TableIndexing::SQUARED_ARG:
    return evaluate(arg * arg, kind, prec, use_elemental_tables);
  case TableIndexing::ARG_OFFSET:
    return evaluate(arg + argument_offset, kind, prec, use_elemental_tables);
  case TableIndexing::SQ_ARG_OFFSET:
    return evaluate((arg + argument_offset) * (arg + argument_offset), kind, prec,
                    use_elemental_tables);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int PPITable::getTableIndex(const double arg, const PrecisionModel prec) const {
  int result;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const Ecumenical8 workspc = { .d = arg };
      result = (workspc.ulli >> (52 - mantissa_bits));
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const Ecumenical4 workspc = { .f = static_cast<float>(arg) };
      result = (workspc.ui >> (23 - mantissa_bits));
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int PPITable::getTableIndexByRealArg(const double arg, const PrecisionModel prec) const {
  switch (indexing_method) {
  case TableIndexing::ARG:
    return getTableIndex(arg, prec);
  case TableIndexing::SQUARED_ARG:
    return getTableIndex(arg * arg, prec);
  case TableIndexing::ARG_OFFSET:
    return getTableIndex(arg + argument_offset, prec);
  case TableIndexing::SQ_ARG_OFFSET:
    return getTableIndex((arg + argument_offset) * (arg + argument_offset), prec);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const PPIKit<double, double4> PPITable::dpData(const HybridTargetLevel tier) const {
  const int index_bit_count = mantissa_bits + 12;
  const int index_shift = 64 - index_bit_count;
  return PPIKit<double, double4>(theme, basis_set, indexing_method, dp_exclusion_offset / 2,
                                 dp_exclusion_offset, index_shift, dp_detail_bitmask,
                                 sp_detail_bitmask, argument_offset, energy.data(tier),
                                 force.data(tier), energy_with_exclusions.data(tier),
                                 force_with_exclusions.data(tier));
}

//-------------------------------------------------------------------------------------------------
const PPIKit<float, float4> PPITable::spData(const HybridTargetLevel tier) const {
  const int index_bit_count = mantissa_bits + 9;
  const int index_shift = 32 - index_bit_count;
  return PPIKit<float, float4>(theme, basis_set, indexing_method, sp_exclusion_offset / 2,
                               sp_exclusion_offset, index_shift, dp_detail_bitmask,
                               sp_detail_bitmask, argument_offset, sp_energy.data(tier),
                               sp_force.data(tier), sp_energy_with_exclusions.data(tier),
                               sp_force_with_exclusions.data(tier));
}

//-------------------------------------------------------------------------------------------------
const PPIeKit<double> PPITable::dpeData(const HybridTargetLevel tier) const {
  const int index_bit_count = mantissa_bits + 12;
  const int index_shift = 64 - index_bit_count;
  return PPIeKit<double>(theme, basis_set, indexing_method, dp_exclusion_offset / 2,
                         dp_exclusion_offset, index_shift, dp_detail_bitmask, sp_detail_bitmask,
                         argument_offset, energy_x.data(tier), energy_y.data(tier),
                         energy_z.data(tier), energy_w.data(tier), force_x.data(tier),
                         force_y.data(tier), force_z.data(tier), force_w.data(tier),
                         energy_with_excl_x.data(tier), energy_with_excl_y.data(tier),
                         energy_with_excl_z.data(tier), energy_with_excl_w.data(tier),
                         force_with_excl_x.data(tier), force_with_excl_y.data(tier),
                         force_with_excl_z.data(tier), force_with_excl_w.data(tier));
}

//-------------------------------------------------------------------------------------------------
const PPIeKit<float> PPITable::speData(const HybridTargetLevel tier) const {
  const int index_bit_count = mantissa_bits + 9;
  const int index_shift = 32 - index_bit_count;
  return PPIeKit<float>(theme, basis_set, indexing_method, sp_exclusion_offset / 2,
                        sp_exclusion_offset, index_shift, dp_detail_bitmask, sp_detail_bitmask,
                        argument_offset, sp_energy_x.data(tier), sp_energy_y.data(tier),
                        sp_energy_z.data(tier), sp_energy_w.data(tier), sp_force_x.data(tier),
                        sp_force_y.data(tier), sp_force_z.data(tier), sp_force_w.data(tier),
                        sp_energy_with_excl_x.data(tier), sp_energy_with_excl_y.data(tier),
                        sp_energy_with_excl_z.data(tier), sp_energy_with_excl_w.data(tier),
                        sp_force_with_excl_x.data(tier), sp_force_with_excl_y.data(tier),
                        sp_force_with_excl_z.data(tier), sp_force_with_excl_w.data(tier));
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void PPITable::upload() {
  coeffs.upload();
  sp_coeffs.upload();
}
  
//-------------------------------------------------------------------------------------------------
void PPITable::download() {
  coeffs.download();
  sp_coeffs.download();
}
#endif

//-------------------------------------------------------------------------------------------------
LogSplineForm PPITable::findMissingForm(const uint holdings) const {
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

} // namespace energy
} // namespace stormm
