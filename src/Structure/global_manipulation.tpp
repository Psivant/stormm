// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

using topology::UnitCellType;
  
//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> beardRotationMatrix(const T om_x, const T om_y, const T om_z) {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool t_is_double = (ct == double_type_index);

  std::vector<T> result(9);

  // Convenient quantities
  const T om2_x = om_x * om_x;
  const T om2_y = om_y * om_y;
  const T om2_z = om_z * om_z;
  const T om2 = om2_x + om2_y + om2_z;
  const T om = (t_is_double) ? sqrt(om2) : sqrtf(om2);
  const T cos_om = (t_is_double) ? cos(om) : cosf(om);
  const T sin_om = (t_is_double) ? sin(om) : sinf(om);
  const T inv_om = 1.0 / om;
  const T inv_om2 = inv_om * inv_om;

  // Compute rotation matrix
  result[0] = (((om2_y + om2_z) * cos_om) + om2_x) * inv_om2;
  result[3] = ((om_x * om_y * inv_om2) * (1.0 - cos_om)) + ((om_z * inv_om) * sin_om);
  result[6] = ((om_x * om_z * inv_om2) * (1.0 - cos_om)) - ((om_y * inv_om) * sin_om);
  result[1] = ((om_x * om_y * inv_om2) * (1.0 - cos_om)) - ((om_z * inv_om) * sin_om);
  result[4] = (((om2_x + om2_z) * cos_om) + om2_y) * inv_om2;
  result[7] = ((om_y * om_z * inv_om2) * (1.0 - cos_om)) + ((om_x * inv_om) * sin_om);
  result[2] = ((om_x * om_z * inv_om2) * (1.0 - cos_om)) + ((om_y * inv_om) * sin_om);
  result[5] = ((om_y * om_z * inv_om2) * (1.0 - cos_om)) - ((om_x * inv_om) * sin_om);
  result[8] = (((om2_x + om2_y) * cos_om) + om2_z) * inv_om2;
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void rotateCoordinates(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const Tcalc alpha,
                       const Tcalc beta, const Tcalc gamma, const int lower_limit,
                       const int upper_limit, const VirtualSiteKit<Tcalc> *vsk,
                       const Tcalc globalpos_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const std::vector<Tcalc> rmat = beardRotationMatrix(alpha, beta, gamma);
  if (isSignedIntegralScalarType<Tcoord>()) {
    const Tcalc value_one = 1.0;
    const Tcalc inv_globalpos_scale_factor = value_one / globalpos_scale_factor;
    for (int i = lower_limit; i < upper_limit; i++) {
      const Tcalc xc = static_cast<Tcalc>(xcrd[i]) * inv_globalpos_scale_factor;
      const Tcalc yc = static_cast<Tcalc>(ycrd[i]) * inv_globalpos_scale_factor;
      const Tcalc zc = static_cast<Tcalc>(zcrd[i]) * inv_globalpos_scale_factor;
      xcrd[i] = llround(((rmat[0] * xc) + (rmat[3] * yc) + (rmat[6] * zc)) *
                        globalpos_scale_factor);
      ycrd[i] = llround(((rmat[1] * xc) + (rmat[4] * yc) + (rmat[7] * zc)) *
                        globalpos_scale_factor);
      zcrd[i] = llround(((rmat[2] * xc) + (rmat[5] * yc) + (rmat[8] * zc)) *
                        globalpos_scale_factor);
    }
  }
  else {
    for (int i = lower_limit; i < upper_limit; i++) {
      const Tcoord xc = xcrd[i];
      const Tcoord yc = ycrd[i];
      const Tcoord zc = zcrd[i];
      xcrd[i] = (rmat[0] * xc) + (rmat[3] * yc) + (rmat[6] * zc);
      ycrd[i] = (rmat[1] * xc) + (rmat[4] * yc) + (rmat[7] * zc);
      zcrd[i] = (rmat[2] * xc) + (rmat[5] * yc) + (rmat[8] * zc);
    }
  }

  // Reposition virtual sites, if information is present to do so.  This assumes no unit cell
  // with periodic boundary conditions, as is approriate for coordinates undergoing global
  // rotation.  Dereferencing the pointer is not costly so long as the abstract contains only
  // constants and pointers.
  if (vsk != nullptr) {
    placeVirtualSites<Tcoord, Tcalc>(xcrd, ycrd, zcrd, nullptr, nullptr, UnitCellType::NONE,
                                     *vsk);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void rotateCoordinates(CoordinateSeriesWriter<Tcoord> *csw, const size_t frame_index,
                       const VirtualSiteKit<Tcalc> &vsk, const Tcalc alpha, const Tcalc beta,
                       const Tcalc gamma, const int lower_limit, const int upper_limit) {
  const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit : csw->natom;
  const size_t atom_start = static_cast<size_t>(roundUp(csw->natom, warp_size_int)) * frame_index;
  rotateCoordinates<Tcoord, Tcalc>(&csw->xcrd[atom_start], &csw->ycrd[atom_start],
                                   &csw->zcrd[atom_start], alpha, beta, gamma, lower_limit,
                                   actual_upper_limit);
  placeVirtualSites<Tcoord, Tcalc>(&csw->xcrd[atom_start], &csw->ycrd[atom_start],
                                   &csw->zcrd[atom_start], nullptr, nullptr, UnitCellType::NONE,
                                   vsk);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void rotateCoordinates(CoordinateSeries<Tcoord> *cs, const size_t frame_index, const AtomGraph &ag,
                       const Tcalc alpha, const Tcalc beta, const Tcalc gamma,
                       const int lower_limit, const int upper_limit) {

  // The calculation type is either single- or double-precision
  CoordinateSeriesWriter<Tcoord> csw = cs->data();
  if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
    const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
    rotateCoordinates<Tcoord, double>(&csw, frame_index, vsk, alpha, beta, gamma, lower_limit,
                                      upper_limit);
  }
  else {
    const VirtualSiteKit<float> vsk = ag.getSinglePrecisionVirtualSiteKit();
    rotateCoordinates<Tcoord, float>(&csw, frame_index, vsk, alpha, beta, gamma, lower_limit,
                                     upper_limit);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void rotateCoordinates(PhaseSpaceSynthesis *poly_ps, const int system_index, const Tcalc alpha,
                       const Tcalc beta, const Tcalc gamma, const int lower_limit,
                       const int upper_limit) {
  const AtomGraph *ag = poly_ps->getSystemTopologyPointer(system_index);
  if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
    const VirtualSiteKit<double> vsk = ag->getDoublePrecisionVirtualSiteKit();
    CoordinateFrame cf = poly_ps->exportCoordinates(system_index);
    CoordinateFrameWriter cfr = cf.data();
    const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit : cfr.natom;
    rotateCoordinates<double, double>(cfr.xcrd, cfr.ycrd, cfr.zcrd, alpha, beta, gamma,
                                      lower_limit, actual_upper_limit, &vsk);
    poly_ps->import(cf, system_index, TrajectoryKind::POSITIONS);
  }
  else {
    const VirtualSiteKit<float> vsk = ag->getSinglePrecisionVirtualSiteKit();
    CoordinateFrame cf = poly_ps->exportCoordinates(system_index);
    CoordinateFrameWriter cfr = cf.data();
    const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit : cfr.natom;
    rotateCoordinates<double, float>(cfr.xcrd, cfr.ycrd, cfr.zcrd, alpha, beta, gamma,
                                      lower_limit, actual_upper_limit, &vsk);
    poly_ps->import(cf, system_index, TrajectoryKind::POSITIONS);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void rotateCoordinates(CondensateWriter *cdnsw, const int system_index,
                       const VirtualSiteKit<Tcalc> &vsk, const Tcalc alpha, const Tcalc beta,
                       const Tcalc gamma, const int lower_limit, const int upper_limit) {
  const size_t atom_start = cdnsw->atom_starts[system_index];
  const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit :
                                                               cdnsw->atom_counts[system_index];
  switch (cdnsw->mode) {
  case PrecisionModel::DOUBLE:
    rotateCoordinates<double, Tcalc>(&cdnsw->xcrd[atom_start], &cdnsw->ycrd[atom_start],
                                     &cdnsw->zcrd[atom_start], alpha, beta, gamma, lower_limit,
                                     actual_upper_limit);
    placeVirtualSites<double, Tcalc>(&cdnsw->xcrd[atom_start], &cdnsw->ycrd[atom_start],
                                     &cdnsw->zcrd[atom_start], nullptr, nullptr,
                                     UnitCellType::NONE, vsk);
    break;
  case PrecisionModel::SINGLE:
    rotateCoordinates<float, Tcalc>(&cdnsw->xcrd_sp[atom_start], &cdnsw->ycrd_sp[atom_start],
                                    &cdnsw->zcrd_sp[atom_start], alpha, beta, gamma,
                                    lower_limit, actual_upper_limit);
    placeVirtualSites<float, Tcalc>(&cdnsw->xcrd_sp[atom_start], &cdnsw->ycrd_sp[atom_start],
                                    &cdnsw->zcrd_sp[atom_start], nullptr, nullptr,
                                    UnitCellType::NONE, vsk);
    break;
  }  
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void translateCoordinates(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const Tcalc xmove,
                          const Tcalc ymove, const Tcalc zmove, const int lower_limit,
                          const int upper_limit, const VirtualSiteKit<Tcalc> *vsk,
                          const Tcalc globalpos_scale_factor) {
  if (isSignedIntegralScalarType<Tcoord>()) {
    for (int i = lower_limit; i < upper_limit; i++) {
      xcrd[i] += llround(xmove * globalpos_scale_factor);
      ycrd[i] += llround(ymove * globalpos_scale_factor);
      zcrd[i] += llround(zmove * globalpos_scale_factor);
    }
  }
  else {
    for (int i = lower_limit; i < upper_limit; i++) {
      xcrd[i] += xmove;
      ycrd[i] += ymove;
      zcrd[i] += zmove;
    }
  }
  
  // Reposition virtual sites, if information is present to do so.  This assumes no unit cell
  // with periodic boundary conditions, as is approriate for coordinates undergoing global
  // rotation.  Dereferencing the pointer is not costly so long as the abstract contains only
  // constants and pointers.
  if (vsk != nullptr) {
    placeVirtualSites<Tcoord, Tcalc>(xcrd, ycrd, zcrd, nullptr, nullptr, UnitCellType::NONE,
                                     *vsk);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void translateCoordinates(PsSynthesisWriter *poly_psw, const int system_index,
                          const VirtualSiteKit<Tcalc> &vsk, const Tcalc xmove,
                          const Tcalc ymove, const Tcalc zmove, const int lower_limit,
                          const int upper_limit) {
  const int llim = poly_psw->atom_starts[system_index] + lower_limit;
  const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit :
                                                               poly_psw->atom_counts[system_index];
  const int hlim = poly_psw->atom_starts[system_index] + actual_upper_limit;

  // Translate the coordinates for the selected atoms
  const int95_t ixmove = hostDoubleToInt95(xmove * poly_psw->gpos_scale);
  const int95_t iymove = hostDoubleToInt95(ymove * poly_psw->gpos_scale);
  const int95_t izmove = hostDoubleToInt95(zmove * poly_psw->gpos_scale);
  for (int i = llim; i < hlim; i++) {
    const int95_t xnew = hostSplitFPSum(ixmove, poly_psw->xcrd[i], poly_psw->xcrd_ovrf[i]);
    const int95_t ynew = hostSplitFPSum(iymove, poly_psw->ycrd[i], poly_psw->ycrd_ovrf[i]);
    const int95_t znew = hostSplitFPSum(izmove, poly_psw->zcrd[i], poly_psw->zcrd_ovrf[i]);
    poly_psw->xcrd[i]      = xnew.x;
    poly_psw->xcrd_ovrf[i] = xnew.y;
    poly_psw->ycrd[i]      = ynew.x;
    poly_psw->ycrd_ovrf[i] = ynew.y;
    poly_psw->zcrd[i]      = znew.x;
    poly_psw->zcrd_ovrf[i] = znew.y;
  }

  // Reposition the virtual sites
  const int natom = poly_psw->atom_counts[system_index];
  CoordinateFrame cf(natom);
  CoordinateFrameWriter cfw = cf.data();
  for (int i = 0; i < natom; i++) {
    cfw.xcrd[i] = hostInt95ToDouble(poly_psw->xcrd[i], poly_psw->xcrd_ovrf[i]) *
                  poly_psw->inv_gpos_scale;
    cfw.ycrd[i] = hostInt95ToDouble(poly_psw->ycrd[i], poly_psw->ycrd_ovrf[i]) *
                  poly_psw->inv_gpos_scale;
    cfw.zcrd[i] = hostInt95ToDouble(poly_psw->zcrd[i], poly_psw->zcrd_ovrf[i]) *
                  poly_psw->inv_gpos_scale;
  }
  placeVirtualSites<double, Tcalc>(cfw.xcrd, cfw.ycrd, cfw.zcrd, nullptr, nullptr,
                                   UnitCellType::NONE, vsk);
  for (int i = 0; i < vsk.nsite; i++) {
    const int vs_idx = vsk.vs_atoms[i];
    const int95_t xvs_new = hostDoubleToInt95(cfw.xcrd[vs_idx] * poly_psw->gpos_scale);
    const int95_t yvs_new = hostDoubleToInt95(cfw.ycrd[vs_idx] * poly_psw->gpos_scale);
    const int95_t zvs_new = hostDoubleToInt95(cfw.zcrd[vs_idx] * poly_psw->gpos_scale);
    poly_psw->xcrd[i]      = xvs_new.x;
    poly_psw->xcrd_ovrf[i] = xvs_new.y;
    poly_psw->ycrd[i]      = yvs_new.x;
    poly_psw->ycrd_ovrf[i] = yvs_new.y;
    poly_psw->zcrd[i]      = zvs_new.x;
    poly_psw->zcrd_ovrf[i] = zvs_new.y;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void translateCoordinates(PhaseSpaceSynthesis *poly_ps, const int system_index, const Tcalc xmove,
                          const Tcalc ymove, const Tcalc zmove, const int lower_limit,
                          const int upper_limit) {
  const AtomGraph *ag = poly_ps->getSystemTopologyPointer(system_index);
  PsSynthesisWriter poly_psw = poly_ps->data();
  if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
    const VirtualSiteKit<double> vsk = ag->getDoublePrecisionVirtualSiteKit();
    translateCoordinates(&poly_psw, system_index, vsk, xmove, ymove, zmove, lower_limit,
                         upper_limit);
  }
  else {
    const VirtualSiteKit<float> vsk = ag->getSinglePrecisionVirtualSiteKit();
    translateCoordinates(&poly_psw, system_index, vsk, xmove, ymove, zmove, lower_limit,
                         upper_limit);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void translateCoordinates(CondensateWriter *cdnsw, const int system_index,
                          const VirtualSiteKit<Tcalc> &vsk, const Tcalc xmove, const Tcalc ymove,
                          const Tcalc zmove, const int lower_limit, const int upper_limit) {
  const size_t atom_start = cdnsw->atom_starts[system_index];
  const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit :
                                                               cdnsw->atom_counts[system_index];
  switch (cdnsw->mode) {
  case PrecisionModel::DOUBLE:
    translateCoordinates<double, Tcalc>(&cdnsw->xcrd[atom_start], &cdnsw->ycrd[atom_start],
                                        &cdnsw->zcrd[atom_start], xmove, ymove, zmove, lower_limit,
                                        actual_upper_limit);
    placeVirtualSites<double, Tcalc>(&cdnsw->xcrd[atom_start], &cdnsw->ycrd[atom_start],
                                     &cdnsw->zcrd[atom_start], nullptr, nullptr,
                                     UnitCellType::NONE, vsk);
    break;
  case PrecisionModel::SINGLE:
    translateCoordinates<float, Tcalc>(&cdnsw->xcrd_sp[atom_start], &cdnsw->ycrd_sp[atom_start],
                                       &cdnsw->zcrd_sp[atom_start], xmove, ymove, zmove,
                                       lower_limit, actual_upper_limit);
    placeVirtualSites<float, Tcalc>(&cdnsw->xcrd_sp[atom_start], &cdnsw->ycrd_sp[atom_start],
                                    &cdnsw->zcrd_sp[atom_start], nullptr, nullptr,
                                    UnitCellType::NONE, vsk);
    break;
  }
}

} // namespace structure
} // namespace stormm
