// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcoord imageValue(const Tcoord x, const Tcalc range, const ImagingMethod style,
                  const Tcalc gpos_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc half = 0.5;
  Tcalc x_frac;
  if (tcoord_is_sgnint) {
    const Tcalc value_one = 1.0;
    const Tcalc inv_gpos_scale_factor = value_one / gpos_scale_factor;
    x_frac = static_cast<Tcalc>(x) * inv_gpos_scale_factor / range;
  }
  else {
    x_frac = x / range;
  }
  switch (style) {
  case ImagingMethod::PRIMARY_UNIT_CELL:
    if (tcalc_is_double) {
      if (tcoord_is_sgnint) {
        return llround((x_frac - floor(x_frac)) * range * gpos_scale_factor);
      }
      else {
        return ((x_frac - floor(x_frac)) * range);
      }
    }
    else {
      if (tcoord_is_sgnint) {
        return llround((x_frac - floorf(x_frac)) * range * gpos_scale_factor);
      }
      else {
        return ((x_frac - floorf(x_frac)) * range);
      }
    }
  case ImagingMethod::MINIMUM_IMAGE:
    if (tcalc_is_double) {
      x_frac -= ((x_frac >= half) *  ceil(x_frac - half)) +
                ((x_frac < -half) * floor(x_frac + half));
    }
    else {
      x_frac -= ((x_frac >= half) *  ceilf(x_frac - half)) +
                ((x_frac < -half) * floorf(x_frac + half));
    }

    // The final subtraction covers the case of the re-imaged coordinate sitting right on 0.5
    if (tcoord_is_sgnint) {
      return llround((x_frac - (x_frac >= half)) * range * gpos_scale_factor);
    }
    else {
      return (x_frac - (x_frac >= half)) * range;
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void imageCoordinates(Tcoord *x, Tcoord *y, Tcoord *z, const double* umat, const double* invu,
                      const UnitCellType unit_cell, const ImagingMethod style,
                      const Tcalc gpos_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc half = 0.5;
  switch (unit_cell) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
    {
      Tcalc local_x, local_y, local_z;
      if (tcoord_is_sgnint) {
        const Tcalc value_one = 1.0;
        const Tcalc inv_gpos_scale_factor = value_one / gpos_scale_factor;
        local_x = static_cast<Tcalc>((*x) * inv_gpos_scale_factor) * umat[0];
        local_y = static_cast<Tcalc>((*y) * inv_gpos_scale_factor) * umat[4];
        local_z = static_cast<Tcalc>((*z) * inv_gpos_scale_factor) * umat[8];
      }
      else {
        local_x = (*x) * umat[0];
        local_y = (*y) * umat[4];
        local_z = (*z) * umat[8];
      }
      switch (style) {
      case ImagingMethod::PRIMARY_UNIT_CELL:
        if (tcalc_is_double) {
          local_x -= floor(local_x);
          local_y -= floor(local_y);
          local_z -= floor(local_z);
        }
        else {
          local_x -= floorf(local_x);
          local_y -= floorf(local_y);
          local_z -= floorf(local_z);
        }
        break;
      case ImagingMethod::MINIMUM_IMAGE:
        if (tcalc_is_double) {
          local_x -= ((local_x >= half) *  ceil(local_x - half)) +
                     ((local_x < -half) * floor(local_x + half));
          local_y -= ((local_y >= half) *  ceil(local_y - half)) +
                     ((local_y < -half) * floor(local_y + half));
          local_z -= ((local_z >= half) *  ceil(local_z - half)) +
                     ((local_z < -half) * floor(local_z + half));
        }
        else {
          local_x -= ((local_x >= half) *  ceilf(local_x - half)) +
                     ((local_x < -half) * floorf(local_x + half));
          local_y -= ((local_y >= half) *  ceilf(local_y - half)) +
                     ((local_y < -half) * floorf(local_y + half));
          local_z -= ((local_z >= half) *  ceilf(local_z - half)) +
                     ((local_z < -half) * floorf(local_z + half));
        }
        local_x -= (local_x >= half);
        local_y -= (local_y >= half);
        local_z -= (local_z >= half);
        break;
      }
      if (tcoord_is_sgnint) {
        *x = llround(local_x * invu[0] * gpos_scale_factor);
        *y = llround(local_y * invu[4] * gpos_scale_factor);
        *z = llround(local_z * invu[8] * gpos_scale_factor);
      }
      else {
        *x = local_x * invu[0];
        *y = local_y * invu[4];
        *z = local_z * invu[8];
      }
    }
    break;
  case UnitCellType::TRICLINIC:
    {      
      Tcalc local_x, local_y, local_z;
      if (tcoord_is_sgnint) {
        const Tcalc value_one = 1.0;
        const Tcalc inv_gpos_scale_factor = value_one / gpos_scale_factor;
        local_x = static_cast<Tcalc>((*x) * inv_gpos_scale_factor);
        local_y = static_cast<Tcalc>((*y) * inv_gpos_scale_factor);
        local_z = static_cast<Tcalc>((*z) * inv_gpos_scale_factor);
      }
      else {
        local_x = (*x);
        local_y = (*y);
        local_z = (*z);
      }
      Tcalc ndx = (umat[0] * local_x) + (umat[3] * local_y) + (umat[6] * local_z);
      Tcalc ndy =                       (umat[4] * local_y) + (umat[7] * local_z);
      Tcalc ndz =                                             (umat[8] * local_z);
      switch (style) {
      case ImagingMethod::PRIMARY_UNIT_CELL:
        if (tcalc_is_double) {
          ndx -= floor(ndx);
          ndy -= floor(ndy);
          ndz -= floor(ndz);
        }
        else {
          ndx -= floorf(ndx);
          ndy -= floorf(ndy);
          ndz -= floorf(ndz);
        }
        break;
      case ImagingMethod::MINIMUM_IMAGE:
        if (tcalc_is_double) {
          ndx -= ((ndx >= half) * ceil(ndx - half)) + ((ndx < -half) * floor(ndx + half));
          ndy -= ((ndy >= half) * ceil(ndy - half)) + ((ndy < -half) * floor(ndy + half));
          ndz -= ((ndz >= half) * ceil(ndz - half)) + ((ndz < -half) * floor(ndz + half));
        }
        else {
          ndx -= ((ndx >= half) * ceilf(ndx - half)) + ((ndx < -half) * floorf(ndx + half));
          ndy -= ((ndy >= half) * ceilf(ndy - half)) + ((ndy < -half) * floorf(ndy + half));
          ndz -= ((ndz >= half) * ceilf(ndz - half)) + ((ndz < -half) * floorf(ndz + half));
        }
        ndx -= (ndx >= half);
        ndy -= (ndy >= half);
        ndz -= (ndz >= half);
        break;
      }
      if (tcoord_is_sgnint) {
        *x = llround(((invu[0] * ndx) + (invu[3] * ndy) + (invu[6] * ndz)) * gpos_scale_factor);
        *y = llround((                  (invu[4] * ndy) + (invu[7] * ndz)) * gpos_scale_factor);
        *z = llround((                                    (invu[8] * ndz)) * gpos_scale_factor);
      }
      else {
        *x = (invu[0] * ndx) + (invu[3] * ndy) + (invu[6] * ndz);
        *y =                   (invu[4] * ndy) + (invu[7] * ndz);
        *z =                                     (invu[8] * ndz);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void imageCoordinates(Tcoord* x, Tcoord* y, Tcoord* z, const int length, const double* umat,
                      const double* invu, const UnitCellType unit_cell, const ImagingMethod style,
                      const Tcalc gpos_scale_factor, int* x_ovrf, int* y_ovrf, int* z_ovrf) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc half = 0.5;
  const Tcalc value_one = 1.0;
  const Tcalc inv_gpos_scale_factor = value_one / gpos_scale_factor;
  const bool xtra_bits = (x_ovrf == nullptr && y_ovrf == nullptr && z_ovrf == nullptr);
  switch (unit_cell) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
    for (int i = 0; i < length; i++) {
      Tcalc local_x, local_y, local_z;
      if (tcoord_is_sgnint) {
        if (xtra_bits) {
          local_x = static_cast<Tcalc>(hostInt95ToDouble(x[i], x_ovrf[i]) * umat[0] *
                                       inv_gpos_scale_factor);
          local_y = static_cast<Tcalc>(hostInt95ToDouble(y[i], y_ovrf[i]) * umat[4] *
                                       inv_gpos_scale_factor);
          local_z = static_cast<Tcalc>(hostInt95ToDouble(z[i], z_ovrf[i]) * umat[8] *
                                       inv_gpos_scale_factor);
        }
        else {
          local_x = static_cast<Tcalc>(x[i] * umat[0] * inv_gpos_scale_factor);
          local_y = static_cast<Tcalc>(y[i] * umat[4] * inv_gpos_scale_factor);
          local_z = static_cast<Tcalc>(z[i] * umat[8] * inv_gpos_scale_factor);
        }
      }
      else {
        local_x = x[i] * umat[0];
        local_y = y[i] * umat[4];
        local_z = z[i] * umat[8];
      }
      switch (style) {
      case ImagingMethod::PRIMARY_UNIT_CELL:
        if (tcalc_is_double) {
          local_x -= floor(local_x);
          local_y -= floor(local_y);
          local_z -= floor(local_z);
        }
        else {
          local_x -= floorf(local_x);
          local_y -= floorf(local_y);
          local_z -= floorf(local_z);
        }
        break;
      case ImagingMethod::MINIMUM_IMAGE:
        if (tcalc_is_double) {
          local_x -= ((local_x >= half) *  ceil(local_x - half)) +
                     ((local_x < -half) * floor(local_x + half));
          local_y -= ((local_y >= half) *  ceil(local_y - half)) +
                     ((local_y < -half) * floor(local_y + half));
          local_z -= ((local_z >= half) *  ceil(local_z - half)) +
                     ((local_z < -half) * floor(local_z + half));
        }
        else {
          local_x -= ((local_x >= half) *  ceilf(local_x - half)) +
                     ((local_x < -half) * floorf(local_x + half));
          local_y -= ((local_y >= half) *  ceilf(local_y - half)) +
                     ((local_y < -half) * floorf(local_y + half));
          local_z -= ((local_z >= half) *  ceilf(local_z - half)) +
                     ((local_z < -half) * floorf(local_z + half));
        }
        local_x -= (local_x >= half);
        local_y -= (local_y >= half);
        local_z -= (local_z >= half);
        break;
      }
      if (tcoord_is_sgnint) {
        if (xtra_bits) {
          const int95_t txi = hostDoubleToInt95(local_x * invu[0] * gpos_scale_factor);
          const int95_t tyi = hostDoubleToInt95(local_y * invu[4] * gpos_scale_factor);
          const int95_t tzi = hostDoubleToInt95(local_z * invu[8] * gpos_scale_factor);
          x[i] = txi.x;
          y[i] = tyi.x;
          z[i] = tzi.x;
          x_ovrf[i] = txi.y;
          y_ovrf[i] = tyi.y;
          z_ovrf[i] = tzi.y;
        }
        else {
          x[i] = llround(local_x * invu[0] * gpos_scale_factor);
          y[i] = llround(local_y * invu[4] * gpos_scale_factor);
          z[i] = llround(local_z * invu[8] * gpos_scale_factor);
        }
      }
      else {
        x[i] = local_x * invu[0];
        y[i] = local_y * invu[4];
        z[i] = local_z * invu[8];
      }
    }
    break;
  case UnitCellType::TRICLINIC:
    for (int i = 0; i < length; i++) {
      Tcalc local_x, local_y, local_z;
      if (tcoord_is_sgnint) {
        if (xtra_bits) {
          local_x = static_cast<Tcalc>(hostInt95ToDouble(x[i], x_ovrf[i])) * inv_gpos_scale_factor;
          local_y = static_cast<Tcalc>(hostInt95ToDouble(y[i], y_ovrf[i])) * inv_gpos_scale_factor;
          local_z = static_cast<Tcalc>(hostInt95ToDouble(z[i], z_ovrf[i])) * inv_gpos_scale_factor;
        }
        else {
          local_x = static_cast<Tcalc>(x[i]) * inv_gpos_scale_factor;
          local_y = static_cast<Tcalc>(y[i]) * inv_gpos_scale_factor;
          local_z = static_cast<Tcalc>(z[i]) * inv_gpos_scale_factor;
        }
      }
      else {
        local_x = x[i];
        local_y = y[i];
        local_z = z[i];
      }
      Tcalc ndx = (umat[0] * local_x) + (umat[3] * local_y) + (umat[6] * local_z);
      Tcalc ndy =                       (umat[4] * local_y) + (umat[7] * local_z);
      Tcalc ndz =                                             (umat[8] * local_z);
      switch (style) {
      case ImagingMethod::PRIMARY_UNIT_CELL:
        if (tcalc_is_double) {
          ndx -= floor(ndx);
          ndy -= floor(ndy);
          ndz -= floor(ndz);
        }
        else {
          ndx -= floorf(ndx);
          ndy -= floorf(ndy);
          ndz -= floorf(ndz);
        }
        break;
      case ImagingMethod::MINIMUM_IMAGE:
        if (tcalc_is_double) {
          ndx -= ((ndx >= half) * ceil(ndx - half)) + ((ndx < -half) * floor(ndx + half));
          ndy -= ((ndy >= half) * ceil(ndy - half)) + ((ndy < -half) * floor(ndy + half));
          ndz -= ((ndz >= half) * ceil(ndz - half)) + ((ndz < -half) * floor(ndz + half));
        }
        else {
          ndx -= ((ndx >= half) * ceilf(ndx - half)) + ((ndx < -half) * floorf(ndx + half));
          ndy -= ((ndy >= half) * ceilf(ndy - half)) + ((ndy < -half) * floorf(ndy + half));
          ndz -= ((ndz >= half) * ceilf(ndz - half)) + ((ndz < -half) * floorf(ndz + half));
        }
        ndx -= (ndx >= half);
        ndy -= (ndy >= half);
        ndz -= (ndz >= half);
        break;
      }
      if (tcoord_is_sgnint) {
        if (xtra_bits) {
          const double t_rx = ((invu[0] * ndx) + (invu[3] * ndy) + (invu[6] * ndz)) *
                              gpos_scale_factor;
          const double t_ry = ((invu[4] * ndy) + (invu[7] * ndz)) * gpos_scale_factor;
          const double t_rz = ((invu[8] * ndz)) * gpos_scale_factor;
          const int95_t it_rx = hostDoubleToInt95(t_rx);
          const int95_t it_ry = hostDoubleToInt95(t_ry);
          const int95_t it_rz = hostDoubleToInt95(t_rz);
          x[i] = it_rx.x;
          y[i] = it_ry.x;
          z[i] = it_rz.x;
          x_ovrf[i] = it_rx.y;
          y_ovrf[i] = it_ry.y;
          z_ovrf[i] = it_rz.y;
        }
        else {
          x[i] = llround(((invu[0] * ndx) + (invu[3] * ndy) + (invu[6] * ndz)) *
                         gpos_scale_factor);
          y[i] = llround(((invu[4] * ndy) + (invu[7] * ndz)) * gpos_scale_factor);
          z[i] = llround(((invu[8] * ndz)) * gpos_scale_factor);
        }
      }
      else {
        x[i] = (invu[0] * ndx) + (invu[3] * ndy) + (invu[6] * ndz);
        y[i] =                   (invu[4] * ndy) + (invu[7] * ndz);
        z[i] =                                     (invu[8] * ndz);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void imageCoordinates(std::vector<Tcoord> *x, std::vector<Tcoord> *y, std::vector<Tcoord> *z,
                      const double* umat, const double* invu, const UnitCellType unit_cell,
                      const ImagingMethod style, const Tcalc gpos_scale_factor) {
  const size_t length = x->size();
  if (length != y->size() || length != z->size()) {
    rtErr("Vectors for x, y, and z coordinates must be the same length for re-imaging.",
          "imageCoordinates");
  }
  imageCoordinates(x->data(), y->data(), z->data(), length, umat, invu, unit_cell, style,
                   gpos_scale_factor);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void imageCoordinates(Hybrid<Tcoord> *x, Hybrid<Tcoord> *y, Hybrid<Tcoord> *z,
                      const double* umat, const double* invu, const UnitCellType unit_cell,
                      const ImagingMethod style, const Tcalc gpos_scale_factor) {
  const size_t length = x->size();
  if (length != y->size() || length != z->size()) {
    rtErr("Vectors for x, y, and z coordinates must be the same length for re-imaging.",
          "imageCoordinates");
  }
  imageCoordinates(x->data(), y->data(), z->data(), length, umat, invu, unit_cell, style,
                   gpos_scale_factor);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void imageCoordinates(CoordinateSeriesWriter<Tcoord> csw, const size_t frame_index,
                      const ImagingMethod style) {
  const size_t atom_offset = frame_index * roundUp<size_t>(csw.natom, warp_size_zu);
  const size_t xfrm_offset = frame_index * roundUp<size_t>(9, warp_size_zu);
  imageCoordinates<Tcoord, Tcalc>(&csw.xcrd[atom_offset], &csw.ycrd[atom_offset],
                                  &csw.zcrd[atom_offset], csw.natom, &csw.umat[xfrm_offset],
                                  &csw.invu[xfrm_offset], csw.unit_cell, style, csw.gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void imageCoordinates(CoordinateSeries<Tcoord> *cs, const size_t frame_index,
                      const ImagingMethod style) {
  imageCoordinates<Tcoord, Tcalc>(cs->data(), frame_index, style);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void imageCoordinates(CondensateWriter cdnsw, const int system_index, const ImagingMethod style) {
  const size_t atom_offset = cdnsw.atom_starts[system_index];
  const size_t xfrm_offset = system_index * roundUp<size_t>(9, warp_size_zu);
  switch (cdnsw.mode) {
  case PrecisionModel::DOUBLE:
    return imageCoordinates<double, Tcalc>(&cdnsw.xcrd[atom_offset], &cdnsw.ycrd[atom_offset],
                                           &cdnsw.zcrd[atom_offset],
                                           cdnsw.atom_counts[system_index],
                                           &cdnsw.umat[xfrm_offset], &cdnsw.invu[xfrm_offset],
                                           cdnsw.unit_cell, style);
  case PrecisionModel::SINGLE:
    return imageCoordinates<float, Tcalc>(&cdnsw.xcrd[atom_offset], &cdnsw.ycrd[atom_offset],
                                          &cdnsw.zcrd[atom_offset],
                                          cdnsw.atom_counts[system_index],
                                          &cdnsw.umat[xfrm_offset], &cdnsw.invu[xfrm_offset],
                                          cdnsw.unit_cell, style);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void imageCoordinates(Condensate *cdns, const int system_index, const ImagingMethod style) {
  imageCoordinates<Tcalc>(cdns->data(), system_index, style);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void imageCoordinates(PsSynthesisWriter poly_psw, const int system_index,
                      const ImagingMethod style) {
  const size_t atom_offset = poly_psw.atom_starts[system_index];
  const size_t xfrm_offset = system_index * roundUp<size_t>(9, warp_size_zu);
  if (poly_psw.gpos_bits <= globalpos_scale_nonoverflow_bits) {
    imageCoordinates(&poly_psw.xcrd[atom_offset], &poly_psw.ycrd[atom_offset],
                     &poly_psw.zcrd[atom_offset], poly_psw.atom_counts[system_index],
                     &poly_psw.umat[xfrm_offset], &poly_psw.invu[xfrm_offset],
                     poly_psw.unit_cell, style, poly_psw.gpos_scale);
  }
  else {
    imageCoordinates(&poly_psw.xcrd[atom_offset], &poly_psw.ycrd[atom_offset],
                     &poly_psw.zcrd[atom_offset], poly_psw.atom_counts[system_index],
                     &poly_psw.umat[xfrm_offset], &poly_psw.invu[xfrm_offset],
                     poly_psw.unit_cell, style, poly_psw.gpos_scale,
                     &poly_psw.xcrd_ovrf[atom_offset], &poly_psw.ycrd_ovrf[atom_offset],
                     &poly_psw.zcrd_ovrf[atom_offset]);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void imageCoordinates(PhaseSpaceSynthesis *poly_ps, const int system_index,
                      const ImagingMethod style) {
  imageCoordinates<Tcalc>(poly_ps->data(), system_index, style);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc4, typename Tcalc>
Tcalc4 distance(const int95_t pti_x, const int95_t pti_y, const int95_t pti_z, const int95_t ptj_x,
                const int95_t ptj_y, const int95_t ptj_z, const double* umat, const double* invu,
                const UnitCellType unit_cell, const Tcalc gpos_scale_factor) {
  const int95_t fp_disp_x = hostInt95Sum(ptj_x.x, ptj_x.y, -pti_x.x, -pti_x.y);
  const int95_t fp_disp_y = hostInt95Sum(ptj_y.x, ptj_y.y, -pti_y.x, -pti_y.y);
  const int95_t fp_disp_z = hostInt95Sum(ptj_z.x, ptj_z.y, -pti_z.x, -pti_z.y);
  const Tcalc inv_gpos_factor = static_cast<Tcalc>(1.0) / gpos_scale_factor;
  Tcalc dx = hostInt95ToDouble(fp_disp_x) * inv_gpos_factor;
  Tcalc dy = hostInt95ToDouble(fp_disp_y) * inv_gpos_factor;
  Tcalc dz = hostInt95ToDouble(fp_disp_z) * inv_gpos_factor;
  imageCoordinates<Tcalc, Tcalc>(&dx, &dy, &dz, umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);
  const Tcalc r = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) ?
                  sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                  sqrtf((dx * dx) + (dy * dy) + (dz * dz));
  return { r, dx, dy, dz };
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc distance(const int atom_i, const int atom_j, const Tcoord* xcrd, const Tcoord* ycrd,
               const Tcoord* zcrd, const double* umat, const double* invu,
               const UnitCellType unit_cell, const Tcalc gpos_scale_factor, const int* xcrd_ovrf,
               const int* ycrd_ovrf, const int* zcrd_ovrf) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  Tcalc dx, dy, dz;
  if (tcoord_is_sgnint) {
    const Tcalc inv_gpos_factor = static_cast<Tcalc>(1.0) / gpos_scale_factor;
    dx = static_cast<Tcalc>(xcrd[atom_j] - xcrd[atom_i]) * inv_gpos_factor;
    dy = static_cast<Tcalc>(ycrd[atom_j] - ycrd[atom_i]) * inv_gpos_factor;
    dz = static_cast<Tcalc>(zcrd[atom_j] - zcrd[atom_i]) * inv_gpos_factor;
  }
  else {
    dx = xcrd[atom_j] - xcrd[atom_i];
    dy = ycrd[atom_j] - ycrd[atom_i];
    dz = zcrd[atom_j] - zcrd[atom_i];
  }
  imageCoordinates(&dx, &dy, &dz, umat, invu, unit_cell, ImagingMethod::MINIMUM_IMAGE,
                   gpos_scale_factor);
  return (tcalc_is_double) ? sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                             sqrtf((dx * dx) + (dy * dy) + (dz * dz));    
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc distance(const int atom_i, const int atom_j, const CoordinateSeriesReader<Tcoord> &csr,
               const size_t frame_index) {
  const size_t atom_offset = frame_index * roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t xfrm_offset = frame_index * roundUp<size_t>(9, warp_size_zu);
  return distance<Tcoord, Tcalc>(atom_i, atom_j, &csr.xcrd[atom_offset], &csr.ycrd[atom_offset],
                                 &csr.zcrd[atom_offset], &csr.umat[xfrm_offset],
                                 &csr.invu[xfrm_offset], csr.unit_cell, csr.gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc distance(const int atom_i, const int atom_j, const CoordinateSeries<Tcoord> *cs,
               const size_t frame_index) {
  return distance<Tcoord, Tcalc>(atom_i, atom_j, cs->data(), frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc distance(const int atom_i, const int atom_j, const CoordinateSeries<Tcoord> &cs,
               const size_t frame_index) {
  return distance<Tcoord, Tcalc>(atom_i, atom_j, cs.data(), frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc distance(const int atom_i, const int atom_j, const CondensateReader &cdnsr,
               const int system_index) {
  const size_t atom_offset = cdnsr.atom_starts[system_index];
  const size_t xfrm_offset = system_index * roundUp<size_t>(9, warp_size_zu);
  switch (cdnsr.mode) {
  case PrecisionModel::DOUBLE:
    return distance<double, double>(atom_i, atom_j, &cdnsr.xcrd[atom_offset],
                                    &cdnsr.ycrd[atom_offset], &cdnsr.zcrd[atom_offset],
                                    &cdnsr.umat[xfrm_offset], &cdnsr.invu[xfrm_offset],
                                    cdnsr.unit_cell);
  case PrecisionModel::SINGLE:
    return distance<float, float>(atom_i, atom_j, &cdnsr.xcrd_sp[atom_offset],
                                  &cdnsr.ycrd_sp[atom_offset], &cdnsr.zcrd_sp[atom_offset],
                                  &cdnsr.umat[xfrm_offset], &cdnsr.invu[xfrm_offset],
                                  cdnsr.unit_cell);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc distance(const int atom_i, const int atom_j, const Condensate *cdns,
               const int system_index) {
  return distance<Tcalc>(atom_i, atom_j, cdns->data(), system_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc distance(const int atom_i, const int atom_j, const Condensate &cdns,
               const int system_index) {
  return distance<Tcalc>(atom_i, atom_j, cdns.data(), system_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc distance(const int atom_i, const int atom_j, const PsSynthesisReader &poly_psr,
               const int system_index) {
  const size_t atom_offset = poly_psr.atom_starts[system_index];
  const size_t xfrm_offset = system_index * roundUp(9, warp_size_int);
  if (poly_psr.gpos_bits <= globalpos_scale_nonoverflow_bits) {
    return distance<llint, Tcalc>(atom_i, atom_j, &poly_psr.xcrd[atom_offset],
                                  &poly_psr.ycrd[atom_offset], &poly_psr.zcrd[atom_offset],
                                  &poly_psr.umat[xfrm_offset], &poly_psr.invu[xfrm_offset],
                                  poly_psr.unit_cell, poly_psr.gpos_scale);
  }
  else {
    return distance<llint, Tcalc>(atom_i, atom_j, &poly_psr.xcrd[atom_offset],
                                  &poly_psr.ycrd[atom_offset], &poly_psr.zcrd[atom_offset],
                                  &poly_psr.umat[xfrm_offset], &poly_psr.invu[xfrm_offset],
                                  poly_psr.unit_cell, poly_psr.gpos_scale,
                                  &poly_psr.xcrd_ovrf[atom_offset],
                                  &poly_psr.ycrd_ovrf[atom_offset],
                                  &poly_psr.zcrd_ovrf[atom_offset]);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc distance(const int atom_i, const int atom_j, const PhaseSpaceSynthesis *poly_ps,
               const int system_index) {
  return distance<Tcalc>(atom_i, atom_j, poly_ps->data(), system_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc distance(const int atom_i, const int atom_j, const PhaseSpaceSynthesis &poly_ps,
               const int system_index) {
  return distance<Tcalc>(atom_i, atom_j, poly_ps.data(), system_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc angle(const int atom_i, const int atom_j, const int atom_k, const Tcoord* xcrd,
            const Tcoord* ycrd, const Tcoord* zcrd, const double* umat, const double* invu,
            const UnitCellType unit_cell, const Tcalc gpos_scale_factor, const int* xcrd_ovrf,
            const int* ycrd_ovrf, const int* zcrd_ovrf) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc value_one = 1.0;

  // Image all three atoms and put the atom J at the origin
  Tcalc rix, riy, riz, rkx, rky, rkz;
  if (tcoord_is_sgnint) {
    const Tcalc inv_gpos_scale_factor = value_one / gpos_scale_factor;
    rix = static_cast<Tcalc>(xcrd[atom_i] - xcrd[atom_j]) * inv_gpos_scale_factor;
    riy = static_cast<Tcalc>(ycrd[atom_i] - ycrd[atom_j]) * inv_gpos_scale_factor;
    riz = static_cast<Tcalc>(zcrd[atom_i] - zcrd[atom_j]) * inv_gpos_scale_factor;
    rkx = static_cast<Tcalc>(xcrd[atom_k] - xcrd[atom_j]) * inv_gpos_scale_factor;
    rky = static_cast<Tcalc>(ycrd[atom_k] - ycrd[atom_j]) * inv_gpos_scale_factor;
    rkz = static_cast<Tcalc>(zcrd[atom_k] - zcrd[atom_j]) * inv_gpos_scale_factor;
  }
  else {
    rix = xcrd[atom_i] - xcrd[atom_j];
    riy = ycrd[atom_i] - ycrd[atom_j];
    riz = zcrd[atom_i] - zcrd[atom_j];
    rkx = xcrd[atom_k] - xcrd[atom_j];
    rky = ycrd[atom_k] - ycrd[atom_j];
    rkz = zcrd[atom_k] - zcrd[atom_j];
  }
  const ImagingMethod imeth = ImagingMethod::MINIMUM_IMAGE;
  imageCoordinates(&rix, &riy, &riz, umat, invu, unit_cell, imeth, gpos_scale_factor);
  imageCoordinates(&rkx, &rky, &rkz, umat, invu, unit_cell, imeth, gpos_scale_factor);

  // Compute the angle, in radians
  const Tcalc mgba = (rix * rix) + (riy * riy) + (riz * riz);
  const Tcalc mgbc = (rkx * rkx) + (rky * rky) + (rkz * rkz);
  const Tcalc invbabc = (tcalc_is_double) ? value_one / sqrt(mgba * mgbc) :
                                            value_one / sqrtf(mgba * mgbc);
  Tcalc costheta = ((rix * rkx) + (riy * rky) + (riz * rkz)) * invbabc;
  costheta = (costheta < -value_one) ? -value_one : (costheta > value_one) ? value_one : costheta;
  return (tcalc_is_double) ? acos(costheta) : acosf(costheta);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc angle(const int atom_i, const int atom_j, const int atom_k,
            const CoordinateSeriesReader<Tcoord> &csr, const size_t frame_index) {
  const size_t atom_offset = frame_index * roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t xfrm_offset = frame_index * roundUp<size_t>(9, warp_size_zu);
  return angle<Tcoord, Tcalc>(atom_i, atom_j, atom_k, &csr.xcrd[atom_offset],
                              &csr.ycrd[atom_offset], &csr.zcrd[atom_offset],
                              &csr.umat[xfrm_offset], &csr.invu[xfrm_offset], csr.unit_cell,
                              csr.gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc angle(const int atom_i, const int atom_j, const int atom_k,
            const CoordinateSeries<Tcoord> *cs, const size_t frame_index) {
  return angle<Tcoord, Tcalc>(atom_i, atom_j, atom_k, cs->data(), frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc angle(const int atom_i, const int atom_j, const int atom_k,
            const CoordinateSeries<Tcoord> &cs, const size_t frame_index) {
  return angle<Tcoord, Tcalc>(atom_i, atom_j, atom_k, cs.data(), frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc angle(const int atom_i, const int atom_j, const int atom_k, const CondensateReader &cdnsr,
            const int system_index) {
  const size_t atom_offset = cdnsr.atom_starts[system_index];
  const size_t xfrm_offset = system_index * roundUp<size_t>(9, warp_size_zu);
  switch (cdnsr.mode) {
  case PrecisionModel::DOUBLE:
    return angle<double, double>(atom_i, atom_j, atom_k, &cdnsr.xcrd[atom_offset],
                                 &cdnsr.ycrd[atom_offset], &cdnsr.zcrd[atom_offset],
                                 &cdnsr.umat[xfrm_offset], &cdnsr.invu[xfrm_offset],
                                 cdnsr.unit_cell);
  case PrecisionModel::SINGLE:
    return angle<float, float>(atom_i, atom_j, atom_k, &cdnsr.xcrd_sp[atom_offset],
                               &cdnsr.ycrd_sp[atom_offset], &cdnsr.zcrd_sp[atom_offset],
                               &cdnsr.umat[xfrm_offset], &cdnsr.invu[xfrm_offset],
                               cdnsr.unit_cell);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc angle(const int atom_i, const int atom_j, const int atom_k, const Condensate *cdns,
            const int system_index) {
  return angle<Tcalc>(atom_i, atom_j, atom_k, cdns->data(), system_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc angle(const int atom_i, const int atom_j, const int atom_k, const Condensate &cdns,
            const int system_index) {
  return angle<Tcalc>(atom_i, atom_j, atom_k, cdns.data(), system_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc angle(const int atom_i, const int atom_j, const int atom_k,
            const PsSynthesisReader &poly_psr, const int system_index) {
  const size_t atom_offset = poly_psr.atom_starts[system_index];
  const size_t xfrm_offset = system_index * roundUp(9, warp_size_int);
  if (poly_psr.gpos_bits <= globalpos_scale_nonoverflow_bits) {
    return angle<llint, Tcalc>(atom_i, atom_j, atom_k, &poly_psr.xcrd[atom_offset],
                               &poly_psr.ycrd[atom_offset], &poly_psr.zcrd[atom_offset],
                               &poly_psr.umat[xfrm_offset], &poly_psr.invu[xfrm_offset],
                               poly_psr.unit_cell, poly_psr.gpos_scale);
  }
  else {
    return angle<llint, Tcalc>(atom_i, atom_j, atom_k, &poly_psr.xcrd[atom_offset],
                               &poly_psr.ycrd[atom_offset], &poly_psr.zcrd[atom_offset],
                               &poly_psr.umat[xfrm_offset], &poly_psr.invu[xfrm_offset],
                               poly_psr.unit_cell, poly_psr.gpos_scale,
                               &poly_psr.xcrd_ovrf[atom_offset], &poly_psr.ycrd_ovrf[atom_offset],
                               &poly_psr.zcrd_ovrf[atom_offset]);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc angle(const int atom_i, const int atom_j, const int atom_k,
            const PhaseSpaceSynthesis *poly_ps, const int system_index) {
  return angle<Tcalc>(atom_i, atom_j, atom_k, poly_ps->data(), system_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc angle(int atom_i, int atom_j, int atom_k, const PhaseSpaceSynthesis &poly_ps,
            int system_index) {
  return angle<Tcalc>(atom_i, atom_j, atom_k, poly_ps.data(), system_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l, const Tcoord* xcrd,
                    const Tcoord* ycrd, const Tcoord* zcrd, const double* umat, const double* invu,
                    const UnitCellType unit_cell, const Tcalc gpos_scale_factor,
                    const int* xcrd_ovrf, const int* ycrd_ovrf, const int* zcrd_ovrf) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc value_one = 1.0;
  
  // Image all four atoms and put atom K at the origin
  Tcalc rix, riy, riz, rjx, rjy, rjz, rlx, rly, rlz;
  if (tcoord_is_sgnint) {
    const Tcalc inv_gpos_scale_factor = value_one / gpos_scale_factor;
    if (xcrd_ovrf == nullptr && ycrd_ovrf == nullptr && zcrd_ovrf == nullptr) {
      rix = static_cast<Tcalc>(xcrd[atom_i] - xcrd[atom_k]) * inv_gpos_scale_factor;
      riy = static_cast<Tcalc>(ycrd[atom_i] - ycrd[atom_k]) * inv_gpos_scale_factor;
      riz = static_cast<Tcalc>(zcrd[atom_i] - zcrd[atom_k]) * inv_gpos_scale_factor;
      rjx = static_cast<Tcalc>(xcrd[atom_j] - xcrd[atom_k]) * inv_gpos_scale_factor;
      rjy = static_cast<Tcalc>(ycrd[atom_j] - ycrd[atom_k]) * inv_gpos_scale_factor;
      rjz = static_cast<Tcalc>(zcrd[atom_j] - zcrd[atom_k]) * inv_gpos_scale_factor;
      rlx = static_cast<Tcalc>(xcrd[atom_l] - xcrd[atom_k]) * inv_gpos_scale_factor;
      rly = static_cast<Tcalc>(ycrd[atom_l] - ycrd[atom_k]) * inv_gpos_scale_factor;
      rlz = static_cast<Tcalc>(zcrd[atom_l] - zcrd[atom_k]) * inv_gpos_scale_factor;
    }
    else {
      const int95_t it_rix = hostInt95Sum( xcrd[atom_i],  xcrd_ovrf[atom_i],
                                          -xcrd[atom_k], -xcrd_ovrf[atom_k]);
      const int95_t it_riy = hostInt95Sum( ycrd[atom_i],  ycrd_ovrf[atom_i],
                                          -ycrd[atom_k], -ycrd_ovrf[atom_k]);
      const int95_t it_riz = hostInt95Sum( zcrd[atom_i],  zcrd_ovrf[atom_i],
                                          -zcrd[atom_k], -zcrd_ovrf[atom_k]);
      const int95_t it_rjx = hostInt95Sum( xcrd[atom_j],  xcrd_ovrf[atom_j],
                                          -xcrd[atom_k], -xcrd_ovrf[atom_k]);
      const int95_t it_rjy = hostInt95Sum( ycrd[atom_j],  ycrd_ovrf[atom_j],
                                          -ycrd[atom_k], -ycrd_ovrf[atom_k]);
      const int95_t it_rjz = hostInt95Sum( zcrd[atom_j],  zcrd_ovrf[atom_j],
                                          -zcrd[atom_k], -zcrd_ovrf[atom_k]);
      const int95_t it_rlx = hostInt95Sum( xcrd[atom_l],  xcrd_ovrf[atom_l],
                                          -xcrd[atom_k], -xcrd_ovrf[atom_k]);
      const int95_t it_rly = hostInt95Sum( ycrd[atom_l],  ycrd_ovrf[atom_l],
                                          -ycrd[atom_k], -ycrd_ovrf[atom_k]);
      const int95_t it_rlz = hostInt95Sum( zcrd[atom_l],  zcrd_ovrf[atom_l],
                                          -zcrd[atom_k], -zcrd_ovrf[atom_k]);
      rix = hostInt95ToDouble(it_rix) * inv_gpos_scale_factor;
      riy = hostInt95ToDouble(it_riy) * inv_gpos_scale_factor;
      riz = hostInt95ToDouble(it_riz) * inv_gpos_scale_factor;
      rjx = hostInt95ToDouble(it_rjx) * inv_gpos_scale_factor;
      rjy = hostInt95ToDouble(it_rjy) * inv_gpos_scale_factor;
      rjz = hostInt95ToDouble(it_rjz) * inv_gpos_scale_factor;
      rlx = hostInt95ToDouble(it_rlx) * inv_gpos_scale_factor;
      rly = hostInt95ToDouble(it_rly) * inv_gpos_scale_factor;
      rlz = hostInt95ToDouble(it_rlz) * inv_gpos_scale_factor;
    }
  }
  else {
    rix = xcrd[atom_i] - xcrd[atom_k];
    riy = ycrd[atom_i] - ycrd[atom_k];
    riz = zcrd[atom_i] - zcrd[atom_k];
    rjx = xcrd[atom_j] - xcrd[atom_k];
    rjy = ycrd[atom_j] - ycrd[atom_k];
    rjz = zcrd[atom_j] - zcrd[atom_k];
    rlx = xcrd[atom_l] - xcrd[atom_k];
    rly = ycrd[atom_l] - ycrd[atom_k];
    rlz = zcrd[atom_l] - zcrd[atom_k];
  }
  const ImagingMethod imeth = ImagingMethod::MINIMUM_IMAGE;
  imageCoordinates(&rix, &riy, &riz, umat, invu, unit_cell, imeth, gpos_scale_factor);
  imageCoordinates(&rjx, &rjy, &rjz, umat, invu, unit_cell, imeth, gpos_scale_factor);
  imageCoordinates(&rlx, &rly, &rlz, umat, invu, unit_cell, imeth, gpos_scale_factor);

  // Compute the dihedral angle, in radians
  Tcalc ab[3], bc[3], cd[3];
  ab[0] = rjx - rix;
  ab[1] = rjy - riy;
  ab[2] = rjz - riz;
  bc[0] = -rjx;
  bc[1] = -rjy;
  bc[2] = -rjz;
  cd[0] = rlx;
  cd[1] = rly;
  cd[2] = rlz;
  
  // Compute cross products and then the angle between the planes
  Tcalc crabbc[3], crbccd[3], scr[3];
  crossProduct(ab, bc, crabbc);
  crossProduct(bc, cd, crbccd);
  Tcalc costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
  if (tcalc_is_double) {
    costheta /= sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                     (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  }
  else {
    costheta /= sqrtf((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                      (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  }
  crossProduct(crabbc, crbccd, scr);
  costheta = (costheta < -value_one) ? -value_one : (costheta > value_one) ? value_one : costheta;
  return angleVerification(costheta, crabbc, crbccd, bc, scr);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l,
                    const CoordinateSeriesReader<Tcoord> &csr, size_t frame_index) {
  const size_t atom_offset = frame_index * roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t xfrm_offset = frame_index * roundUp<size_t>(9, warp_size_zu);
  return dihedralAngle<Tcoord, Tcalc>(atom_i, atom_j, atom_k, atom_l, &csr.xcrd[atom_offset],
                                      &csr.ycrd[atom_offset], &csr.zcrd[atom_offset],
                                      &csr.umat[xfrm_offset], &csr.invu[xfrm_offset],
                                      csr.unit_cell, csr.gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l,
                    const CoordinateSeries<Tcoord> *cs, size_t frame_index) {
  return dihedralAngle<Tcoord, Tcalc>(atom_i, atom_j, atom_k, atom_l, cs->data(), frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l,
                    const CoordinateSeries<Tcoord> &cs, size_t frame_index) {
  return dihedralAngle<Tcoord, Tcalc>(atom_i, atom_j, atom_k, atom_l, cs.data(), frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l, const CondensateReader &cdnsr,
                    int system_index) {
  const size_t atom_offset = cdnsr.atom_starts[system_index];
  const size_t xfrm_offset = system_index * roundUp<size_t>(9, warp_size_zu);
  switch (cdnsr.mode) {
  case PrecisionModel::DOUBLE:
    return dihedralAngle<double, double>(atom_i, atom_j, atom_k, atom_l, &cdnsr.xcrd[atom_offset],
                                         &cdnsr.ycrd[atom_offset], &cdnsr.zcrd[atom_offset],
                                         &cdnsr.umat[xfrm_offset], &cdnsr.invu[xfrm_offset],
                                         cdnsr.unit_cell);
  case PrecisionModel::SINGLE:
    return dihedralAngle<float, float>(atom_i, atom_j, atom_k, atom_l, &cdnsr.xcrd_sp[atom_offset],
                                       &cdnsr.ycrd_sp[atom_offset], &cdnsr.zcrd_sp[atom_offset],
                                       &cdnsr.umat[xfrm_offset], &cdnsr.invu[xfrm_offset],
                                       cdnsr.unit_cell);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l, const Condensate *cdns,
                    int system_index) {
  return dihedralAngle<Tcalc>(atom_i, atom_j, atom_k, atom_l, cdns->data(), system_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l, const Condensate &cdns,
                    int system_index) {
  return dihedralAngle<Tcalc>(atom_i, atom_j, atom_k, atom_l, cdns.data(), system_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l,
                    const PsSynthesisReader &poly_psr, int system_index) {
  const size_t atom_offset = poly_psr.atom_starts[system_index];
  const size_t xfrm_offset = system_index * roundUp(9, warp_size_int);
  if (poly_psr.gpos_bits <= globalpos_scale_nonoverflow_bits) {
    return dihedralAngle<llint, Tcalc>(atom_i, atom_j, atom_k, atom_l, &poly_psr.xcrd[atom_offset],
                                       &poly_psr.ycrd[atom_offset], &poly_psr.zcrd[atom_offset],
                                       &poly_psr.umat[xfrm_offset], &poly_psr.invu[xfrm_offset],
                                       poly_psr.unit_cell, poly_psr.gpos_scale);
  }
  else {
    return dihedralAngle<llint, Tcalc>(atom_i, atom_j, atom_k, atom_l, &poly_psr.xcrd[atom_offset],
                                       &poly_psr.ycrd[atom_offset], &poly_psr.zcrd[atom_offset],
                                       &poly_psr.umat[xfrm_offset], &poly_psr.invu[xfrm_offset],
                                       poly_psr.unit_cell, poly_psr.gpos_scale,
                                       &poly_psr.xcrd_ovrf[atom_offset],
                                       &poly_psr.ycrd_ovrf[atom_offset],
                                       &poly_psr.zcrd_ovrf[atom_offset]);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l,
                    const PhaseSpaceSynthesis *poly_ps, int system_index) {
  return dihedralAngle<Tcalc>(atom_i, atom_j, atom_k, atom_l, poly_ps->data(), system_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l,
                    const PhaseSpaceSynthesis &poly_ps, int system_index) {
  return dihedralAngle<Tcalc>(atom_i, atom_j, atom_k, atom_l, poly_ps.data(), system_index);
}

} // namespace structure
} // namespace stormm
