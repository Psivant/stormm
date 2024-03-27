// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord> std::vector<Tcoord> MeshParameters::getMeshOrigin() const {
  if (isFloatingPointScalarType<Tcoord>()) {
    std::vector<Tcoord> result(3);
    result[0] = hostInt95ToDouble(origin_x) * inverse_scale_factor;
    result[1] = hostInt95ToDouble(origin_y) * inverse_scale_factor;
    result[2] = hostInt95ToDouble(origin_z) * inverse_scale_factor;
    return result;
  }
  else if (isFloatingPointHpcVectorType<Tcoord>()) {
    rtErr("The mesh coordinate origin is available as an HPC vector type through the "
          "getMeshOriginAsTuple() function.", "MeshParameter", "getMeshOrigin");
  }
  else {
    rtErr("In order to get the mesh coordinate origin as a vector of fixed-precision numbers, use "
          "the getMeshOriginAsFixedPrecision() function.", "MeshParameter", "getMeshOrigin");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T3> T3 MeshParameters::getMeshOriginAsTuple() const {
  if (isFloatingPointScalarType<T3>()) {
    rtErr("The mesh coordinate origin is available as a 3-element array of the desired floating "
          "point type through the getMeshOrigin() function.", "MeshParameter",
          "getMeshOriginAsTuple");
  }
  else if (isFloatingPointHpcVectorType<T3>()) {
    T3 result;
    result.x = hostInt95ToDouble(origin_x) * inverse_scale_factor;
    result.y = hostInt95ToDouble(origin_y) * inverse_scale_factor;
    result.z = hostInt95ToDouble(origin_z) * inverse_scale_factor;
    return result;
  }
  else {
    rtErr("In order to get the mesh coordinate origin as a vector of fixed-precision numbers, use "
          "the getMeshOriginAsFixedPrecision() function.  The fixed-precision representation is "
          "not available as a tuple.", "MeshParameter", "getMeshOriginAsTuple");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord>
std::vector<Tcoord> MeshParameters::getMeshElementVector(const UnitCellAxis dim) const {
  std::vector<Tcoord> result(3);
  const int icol = static_cast<int>(dim);
  if (isFloatingPointScalarType<Tcoord>()) {
    for (int i = 0; i < 3; i++) {
      result[i] = hostInt95ToDouble(fp_element_invu[i + (3 * icol)]) * inverse_scale_factor;
    }
  }
  else {
    rtErr("Mesh element vectors can only be returned as " + getStormmScalarTypeName<float>() +
          " or " + getStormmScalarTypeName<double>() + ".  For the fixed-precision "
          "representation, use getMeshElementVectorFP().", "MeshParameters",
          "getMeshElementVector");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord>
std::vector<Tcoord> MeshParameters::getMeshElementVector(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return getMeshElementVector<Tcoord>(UnitCellAxis::A);
  case CartesianDimension::Y:
    return getMeshElementVector<Tcoord>(UnitCellAxis::B);
  case CartesianDimension::Z:
    return getMeshElementVector<Tcoord>(UnitCellAxis::C);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T3>
T3 MeshParameters::getMeshElementVectorAsTuple(const UnitCellAxis dim) const {
  T3 result;
  const int icol = static_cast<int>(dim);
  if (isFloatingPointHpcVectorType<T3>()) {
    result.x = hostInt95ToDouble(fp_element_invu[(3 * icol)    ]) * inverse_scale_factor;
    result.y = hostInt95ToDouble(fp_element_invu[(3 * icol) + 1]) * inverse_scale_factor;
    result.z = hostInt95ToDouble(fp_element_invu[(3 * icol) + 2]) * inverse_scale_factor;
  }
  else {
    rtErr("Mesh element vectors can only be returned as " + getStormmScalarTypeName<float3>() +
          " or " + getStormmScalarTypeName<double3>() + " tuples.  To get one of the vectors as "
          "a real-valued vector, use getMeshElementVector().  For the fixed-precision "
          "representation, use getMeshElementVectorFP().", "MeshParameters",
          "getMeshElementVector");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T3>
T3 MeshParameters::getMeshElementVectorAsTuple(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return getMeshElementVectorAsTuple<T3>(UnitCellAxis::A);
  case CartesianDimension::Y:
    return getMeshElementVectorAsTuple<T3>(UnitCellAxis::B);
  case CartesianDimension::Z:
    return getMeshElementVectorAsTuple<T3>(UnitCellAxis::C);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord> std::vector<Tcoord> MeshParameters::getMeshTransform() const {
  std::vector<Tcoord> result(9);
  const size_t ct = std::type_index(typeid(Tcoord)).hash_code();
  if (ct == double_type_index) {
    for (size_t i = 0; i < 9LLU; i++) {
      result[i] = element_umat[i];
    }
  }
  else if (ct == float_type_index) {
    for (size_t i = 0; i < 9LLU; i++) {
      result[i] = sp_element_umat[i];
    }
  }
  else {
    rtErr("The transformation matrix into element space is only available in single- or double-"
          "precision floating point numbers.", "MeshParameters", "getMeshTransform");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord> std::vector<Tcoord> MeshParameters::getMeshInverseTransform() const {
  std::vector<Tcoord> result(9);
  const size_t ct = std::type_index(typeid(Tcoord)).hash_code();
  if (isFloatingPointScalarType<Tcoord>()) {
    for (size_t i = 0; i < 9LLU; i++) {
      result[i] = hostInt95ToDouble(fp_element_invu[i]) * inverse_scale_factor;
    }
  }
  else {
    rtErr("The inverse transformation matrix (the column matrix of element vectors) is only "
          "available in fixed-precision format or single- or double-precision floating point "
          "numbers.  To get the fixed-precision format, use getMeshInverseTransformAsFP().",
          "MeshParameters", "getMeshTransform");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshParameters getMeasurements(const AtomGraph &ag, const CoordinateSeries<T> &cs,
                               const double padding, const std::vector<double> &spacing,
                               const int scale_bits_in) {
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  const CoordinateSeriesReader<T> csr = cs.data();
  T xmin, ymin, zmin, xmax, ymax, zmax;
  const bool data_is_int = isSignedIntegralScalarType<T>();
  const T t_padding = (data_is_int) ? llround(padding * csr.gpos_scale) : padding * csr.gpos_scale;
  bool points_unset = true;
  const int lj_idx_offset = nbk.n_lj_types + 1;
  const size_t natom_zu = csr.natom;
  const size_t padded_natom = roundUp(natom_zu, warp_size_zu);
  for (int i = 0; i < csr.nframe; i++) {
    for (size_t pos = 0; pos < natom_zu; pos++) {
      const size_t i_pos = (i * padded_natom) + pos;
      if (ag.getAtomMobility(pos)) {
        continue;
      }
      const size_t plj_idx = lj_idx_offset * nbk.lj_idx[pos];
      const double atom_radius = 0.5 * csr.gpos_scale *
                                 pow(nbk.lja_coeff[plj_idx] / nbk.ljb_coeff[plj_idx], 1.0 / 6.0);
      const T t_radius = (data_is_int) ? llround(atom_radius) : atom_radius;
      if (points_unset) {
        xmin = csr.xcrd[i_pos] - t_radius;
        xmax = csr.xcrd[i_pos] + t_radius;
        ymin = csr.ycrd[i_pos] - t_radius;
        ymax = csr.ycrd[i_pos] + t_radius;
        zmin = csr.zcrd[i_pos] - t_radius;
        zmax = csr.zcrd[i_pos] + t_radius;
        points_unset = false;
      }
      else {
        xmin = std::min(xmin, csr.xcrd[i_pos] - t_radius);
        xmax = std::max(xmax, csr.xcrd[i_pos] + t_radius);
        ymin = std::min(ymin, csr.ycrd[i_pos] - t_radius);
        ymax = std::max(ymax, csr.ycrd[i_pos] + t_radius);
        zmin = std::min(zmin, csr.zcrd[i_pos] - t_radius);
        zmax = std::max(zmax, csr.zcrd[i_pos] + t_radius);
      }
    }
  }
  xmin -= t_padding;
  xmax += t_padding;
  ymin -= t_padding;
  ymax += t_padding;
  zmin -= t_padding;
  zmax += t_padding;
  const std::vector<double> limits = { static_cast<double>(xmin) * csr.inv_gpos_scale,
                                       static_cast<double>(ymin) * csr.inv_gpos_scale,
                                       static_cast<double>(zmin) * csr.inv_gpos_scale,
                                       static_cast<double>(xmax) * csr.inv_gpos_scale,
                                       static_cast<double>(ymax) * csr.inv_gpos_scale,
                                       static_cast<double>(zmax) * csr.inv_gpos_scale };
  return getMeasurements(limits, spacing, scale_bits_in);
}

} // namespace structure
} // namespace stormm
