// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void placeVirtualSites(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const double* umat,
                       const double* invu, const UnitCellType unit_cell,
                       const VirtualSiteKit<Tcalc> &vsk, const Tcalc gpos_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc value_one = 1.0;
  const Tcalc inv_gpos_factor = value_one / gpos_scale_factor;
  for (int i = 0; i < vsk.nsite; i++) {
    const int vsite_atom  = vsk.vs_atoms[i];
    const int parent_atom = vsk.frame1_idx[i];
    const int frame2_atom = vsk.frame2_idx[i];
    const int param_idx = vsk.vs_param_idx[i];
    const ImagingMethod minimum_image = ImagingMethod::MINIMUM_IMAGE;
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[param_idx])) {
    case VirtualSiteKind::FLEX_2:
      {
        Tcalc dx, dy, dz;
        if (tcoord_is_sgnint) {
          dx = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          dy = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          dz = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
        }
        else {
          dx = xcrd[frame2_atom] - xcrd[parent_atom];
          dy = ycrd[frame2_atom] - ycrd[parent_atom];
          dz = zcrd[frame2_atom] - zcrd[parent_atom];
        }
        imageCoordinates<Tcalc, Tcalc>(&dx, &dy, &dz, umat, invu, unit_cell, minimum_image);
        const Tcalc p_f2_factor = vsk.dim1[param_idx];
        if (tcoord_is_sgnint) {
          const Tcalc disp_mult = p_f2_factor * gpos_scale_factor;
          xcrd[vsite_atom] = xcrd[parent_atom] + llround(dx * disp_mult);
          ycrd[vsite_atom] = ycrd[parent_atom] + llround(dy * disp_mult);
          zcrd[vsite_atom] = zcrd[parent_atom] + llround(dz * disp_mult);
        }
        else {
          xcrd[vsite_atom] = xcrd[parent_atom] + (p_f2_factor * dx);
          ycrd[vsite_atom] = ycrd[parent_atom] + (p_f2_factor * dy);
          zcrd[vsite_atom] = zcrd[parent_atom] + (p_f2_factor * dz);
        }
      }
      break;
    case VirtualSiteKind::FIXED_2:
      {
        Tcalc dx, dy, dz;
        if (tcoord_is_sgnint) {
          dx = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          dy = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          dz = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
        }
        else {
          dx = xcrd[frame2_atom] - xcrd[parent_atom];
          dy = ycrd[frame2_atom] - ycrd[parent_atom];
          dz = zcrd[frame2_atom] - zcrd[parent_atom];          
        }
        imageCoordinates<Tcalc, Tcalc>(&dx, &dy, &dz, umat, invu, unit_cell, minimum_image);
        const Tcalc invr = 1.0 / sqrt((dx * dx) + (dy * dy) + (dz * dz));
        const Tcalc p_vs_distance = vsk.dim1[param_idx];
        if (tcoord_is_sgnint) {
          const Tcalc disp_mult = p_vs_distance * invr * gpos_scale_factor;
          xcrd[vsite_atom] = xcrd[parent_atom] + llround(disp_mult * dx);
          ycrd[vsite_atom] = ycrd[parent_atom] + llround(disp_mult * dy);
          zcrd[vsite_atom] = zcrd[parent_atom] + llround(disp_mult * dz);
        }
        else {
          xcrd[vsite_atom] = xcrd[parent_atom] + (p_vs_distance * dx * invr);
          ycrd[vsite_atom] = ycrd[parent_atom] + (p_vs_distance * dy * invr);
          zcrd[vsite_atom] = zcrd[parent_atom] + (p_vs_distance * dz * invr);
        }
      }
      break;
    case VirtualSiteKind::FLEX_3:
      {
        const int frame3_atom = vsk.frame3_idx[i];
        Tcalc dx2, dy2, dz2, dx3, dy3, dz3;
        if (tcoord_is_sgnint) {
          dx2 = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          dy2 = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          dz2 = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          dx3 = static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          dy3 = static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          dz3 = static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[parent_atom]) * inv_gpos_factor;
        }
        else {
          dx2 = xcrd[frame2_atom] - xcrd[parent_atom];
          dy2 = ycrd[frame2_atom] - ycrd[parent_atom];
          dz2 = zcrd[frame2_atom] - zcrd[parent_atom];
          dx3 = xcrd[frame3_atom] - xcrd[parent_atom];
          dy3 = ycrd[frame3_atom] - ycrd[parent_atom];
          dz3 = zcrd[frame3_atom] - zcrd[parent_atom];
        }
        imageCoordinates<Tcalc, Tcalc>(&dx2, &dy2, &dz2, umat, invu, unit_cell, minimum_image);
        imageCoordinates<Tcalc, Tcalc>(&dx3, &dy3, &dz3, umat, invu, unit_cell, minimum_image);
        const Tcalc p_f2_factor = vsk.dim1[param_idx];
        const Tcalc p_f3_factor = vsk.dim2[param_idx];
        if (tcoord_is_sgnint) {
          xcrd[vsite_atom] = xcrd[parent_atom] +
                             llround(((p_f2_factor * dx2) + (p_f3_factor * dx3)) *
                                     gpos_scale_factor);
          ycrd[vsite_atom] = ycrd[parent_atom] +
                             llround(((p_f2_factor * dy2) + (p_f3_factor * dy3)) *
                                     gpos_scale_factor);
          zcrd[vsite_atom] = zcrd[parent_atom] +
                             llround(((p_f2_factor * dz2) + (p_f3_factor * dz3)) *
                                     gpos_scale_factor);
        }
        else {
          xcrd[vsite_atom] = xcrd[parent_atom] + (p_f2_factor * dx2) + (p_f3_factor * dx3);
          ycrd[vsite_atom] = ycrd[parent_atom] + (p_f2_factor * dy2) + (p_f3_factor * dy3);
          zcrd[vsite_atom] = zcrd[parent_atom] + (p_f2_factor * dz2) + (p_f3_factor * dz3);
        }
      }
      break;
    case VirtualSiteKind::FIXED_3:
      {
        const int frame3_atom = vsk.frame3_idx[i];
        Tcalc dx23, dy23, dz23;
        if (tcoord_is_sgnint) {
          dx23 = static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[frame2_atom]) * inv_gpos_factor;
          dy23 = static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[frame2_atom]) * inv_gpos_factor;
          dz23 = static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[frame2_atom]) * inv_gpos_factor;
        }
        else {
          dx23 = xcrd[frame3_atom] - xcrd[frame2_atom];
          dy23 = ycrd[frame3_atom] - ycrd[frame2_atom];
          dz23 = zcrd[frame3_atom] - zcrd[frame2_atom];
        }
        imageCoordinates<Tcalc, Tcalc>(&dx23, &dy23, &dz23, umat, invu, unit_cell, minimum_image);
        const Tcalc f2_f3_factor  = vsk.dim2[param_idx];
        Tcoord x_midpoint, y_midpoint, z_midpoint;
        if (tcoord_is_sgnint) {
          const Tcalc disp_mult = f2_f3_factor * gpos_scale_factor;
          x_midpoint = xcrd[frame2_atom] + llround(disp_mult * dx23);
          y_midpoint = ycrd[frame2_atom] + llround(disp_mult * dy23);
          z_midpoint = zcrd[frame2_atom] + llround(disp_mult * dz23);
        }
        else {
          x_midpoint = xcrd[frame2_atom] + (f2_f3_factor * dx23);
          y_midpoint = ycrd[frame2_atom] + (f2_f3_factor * dy23);
          z_midpoint = zcrd[frame2_atom] + (f2_f3_factor * dz23);
        }
        Tcalc dxm, dym, dzm;
        if (tcoord_is_sgnint) {
          dxm = static_cast<Tcalc>(x_midpoint - xcrd[parent_atom]) * inv_gpos_factor;
          dym = static_cast<Tcalc>(y_midpoint - ycrd[parent_atom]) * inv_gpos_factor;
          dzm = static_cast<Tcalc>(z_midpoint - zcrd[parent_atom]) * inv_gpos_factor;
        }
        else {
          dxm = x_midpoint - xcrd[parent_atom];
          dym = y_midpoint - ycrd[parent_atom];
          dzm = z_midpoint - zcrd[parent_atom];
        }
        imageCoordinates<Tcalc, Tcalc>(&dxm, &dym, &dzm, umat, invu, unit_cell, minimum_image);
        const Tcalc invr = 1.0 / sqrt((dxm * dxm) + (dym * dym) + (dzm * dzm));
        const Tcalc p_vs_distance = vsk.dim1[param_idx];
        if (tcoord_is_sgnint) {
          const Tcalc disp_mult = p_vs_distance * invr * gpos_scale_factor;
          xcrd[vsite_atom] = xcrd[parent_atom] + llround(disp_mult * dxm);
          ycrd[vsite_atom] = ycrd[parent_atom] + llround(disp_mult * dym);
          zcrd[vsite_atom] = zcrd[parent_atom] + llround(disp_mult * dzm);
        }
        else {
          xcrd[vsite_atom] = xcrd[parent_atom] + (p_vs_distance * dxm * invr);
          ycrd[vsite_atom] = ycrd[parent_atom] + (p_vs_distance * dym * invr);
          zcrd[vsite_atom] = zcrd[parent_atom] + (p_vs_distance * dzm * invr);
        }
      }
      break;
    case VirtualSiteKind::FAD_3:
      {
        // Use small vectors to make use of a vector operation
        const int frame3_atom  = vsk.frame3_idx[i];
        Tcalc p_f2[3], f2_f3[3], f23_t_pf2[3];
        if (tcoord_is_sgnint) {
          p_f2[0]  = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f2[1]  = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f2[2]  = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          f2_f3[0] = static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[frame2_atom]) * inv_gpos_factor;
          f2_f3[1] = static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[frame2_atom]) * inv_gpos_factor;
          f2_f3[2] = static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[frame2_atom]) * inv_gpos_factor;
        }
        else {
          p_f2[0]  = xcrd[frame2_atom] - xcrd[parent_atom];
          p_f2[1]  = ycrd[frame2_atom] - ycrd[parent_atom];
          p_f2[2]  = zcrd[frame2_atom] - zcrd[parent_atom];
          f2_f3[0] = xcrd[frame3_atom] - xcrd[frame2_atom];
          f2_f3[1] = ycrd[frame3_atom] - ycrd[frame2_atom];
          f2_f3[2] = zcrd[frame3_atom] - zcrd[frame2_atom];
        }
        imageCoordinates<Tcalc, Tcalc>( &p_f2[0],  &p_f2[1],  &p_f2[2], umat, invu, unit_cell,
                                       minimum_image);
        imageCoordinates<Tcalc, Tcalc>(&f2_f3[0], &f2_f3[1], &f2_f3[2], umat, invu, unit_cell,
                                       minimum_image);
        const Tcalc invr2_p_f2 = 1.0 / ((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                        (p_f2[2] * p_f2[2]));

        // Compute the projection of f2_f3 (the displacement vector between frame atoms 2 and 3),
        // onto p_f2 (the displacement vector between the parent atom and frame atom 2).  Subtract
        // this from f2_f3 to get the part of f2_f3 that is perpendicular to p_f2.  Use the vector
        // that will ultimately store the result as temporary space to hold the projection.  In
        // more optimized code, the inverse squared displacement p_f2 computed for the projection
        // can be re-used, but it makes things more clear for the reference code to recompute that
        // quantity when needed.
        project(f2_f3, p_f2, f23_t_pf2, 3);
        f23_t_pf2[0] = f2_f3[0] - f23_t_pf2[0];
        f23_t_pf2[1] = f2_f3[1] - f23_t_pf2[1];
        f23_t_pf2[2] = f2_f3[2] - f23_t_pf2[2];
        const Tcalc invr_p_f2 = (tcalc_is_double) ?
                                1.0 / sqrt((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                           (p_f2[2] * p_f2[2])) :
                                value_one / sqrtf((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                                  (p_f2[2] * p_f2[2]));
          
        const Tcalc invr_t = (tcalc_is_double) ?
                             1.0 / sqrt((f23_t_pf2[0] * f23_t_pf2[0]) +
                                        (f23_t_pf2[1] * f23_t_pf2[1]) +
                                        (f23_t_pf2[2] * f23_t_pf2[2])) :
                             value_one / sqrtf((f23_t_pf2[0] * f23_t_pf2[0]) +
                                               (f23_t_pf2[1] * f23_t_pf2[1]) +
                                               (f23_t_pf2[2] * f23_t_pf2[2]));          
        const Tcalc p_f2_factor = vsk.dim1[param_idx] * cos(vsk.dim2[param_idx]) * invr_p_f2;
        const Tcalc t_factor    = vsk.dim1[param_idx] * sin(vsk.dim2[param_idx]) * invr_t;
        if (tcoord_is_sgnint) {
          xcrd[vsite_atom] = xcrd[parent_atom] +
                             llround(((p_f2_factor * p_f2[0]) + (t_factor * f23_t_pf2[0])) *
                                     gpos_scale_factor);
          ycrd[vsite_atom] = ycrd[parent_atom] +
                             llround(((p_f2_factor * p_f2[1]) + (t_factor * f23_t_pf2[1])) *
                                     gpos_scale_factor);
          zcrd[vsite_atom] = zcrd[parent_atom] +
                             llround(((p_f2_factor * p_f2[2]) + (t_factor * f23_t_pf2[2])) *
                                     gpos_scale_factor);
        }
        else {
          xcrd[vsite_atom] = xcrd[parent_atom] +
                             (p_f2_factor * p_f2[0]) + (t_factor * f23_t_pf2[0]);
          ycrd[vsite_atom] = ycrd[parent_atom] +
                             (p_f2_factor * p_f2[1]) + (t_factor * f23_t_pf2[1]);
          zcrd[vsite_atom] = zcrd[parent_atom] +
                             (p_f2_factor * p_f2[2]) + (t_factor * f23_t_pf2[2]);
        }
      }
      break;
    case VirtualSiteKind::OUT_3:
      {
        const int frame3_atom  = vsk.frame3_idx[i];
        Tcalc p_f2[3], p_f3[3], f2_f3[3], pf2_x_pf3[3];
        if (tcoord_is_sgnint) {
          p_f2[0] = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f2[1] = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f2[2] = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          p_f3[0] = static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f3[1] = static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f3[2] = static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[parent_atom]) * inv_gpos_factor;
        }
        else {
          p_f2[0] = xcrd[frame2_atom] - xcrd[parent_atom];
          p_f2[1] = ycrd[frame2_atom] - ycrd[parent_atom];
          p_f2[2] = zcrd[frame2_atom] - zcrd[parent_atom];
          p_f3[0] = xcrd[frame3_atom] - xcrd[parent_atom];
          p_f3[1] = ycrd[frame3_atom] - ycrd[parent_atom];
          p_f3[2] = zcrd[frame3_atom] - zcrd[parent_atom];
        }
        imageCoordinates<Tcalc, Tcalc>(&p_f2[0], &p_f2[1], &p_f2[2], umat, invu, unit_cell,
                                       minimum_image);
        imageCoordinates<Tcalc, Tcalc>(&p_f3[0], &p_f3[1], &p_f3[2], umat, invu, unit_cell,
                                       minimum_image);
        crossProduct(p_f2, p_f3, pf2_x_pf3);
        const Tcalc pf2_factor = vsk.dim1[param_idx];
        const Tcalc pf3_factor = vsk.dim2[param_idx];
        const Tcalc cr_factor  = vsk.dim3[param_idx];
        if (tcoord_is_sgnint) {
          xcrd[vsite_atom] = xcrd[parent_atom] +
                             llround(((pf2_factor * p_f2[0]) + (pf3_factor * p_f3[0]) +
                                     (cr_factor * pf2_x_pf3[0])) * gpos_scale_factor);
          ycrd[vsite_atom] = ycrd[parent_atom] +
                             llround(((pf2_factor * p_f2[1]) + (pf3_factor * p_f3[1]) +
                                     (cr_factor * pf2_x_pf3[1])) * gpos_scale_factor);
          zcrd[vsite_atom] = zcrd[parent_atom] +
                             llround(((pf2_factor * p_f2[2]) + (pf3_factor * p_f3[2]) +
                                     (cr_factor * pf2_x_pf3[2])) * gpos_scale_factor);
        }
        else {
          xcrd[vsite_atom] = xcrd[parent_atom] + (pf2_factor * p_f2[0]) + (pf3_factor * p_f3[0]) +
                             (cr_factor * pf2_x_pf3[0]);
          ycrd[vsite_atom] = ycrd[parent_atom] + (pf2_factor * p_f2[1]) + (pf3_factor * p_f3[1]) +
                             (cr_factor * pf2_x_pf3[1]);
          zcrd[vsite_atom] = zcrd[parent_atom] + (pf2_factor * p_f2[2]) + (pf3_factor * p_f3[2]) +
                             (cr_factor * pf2_x_pf3[2]);
        }
      }
      break;
    case VirtualSiteKind::FIXED_4:
      {
        const int frame3_atom  = vsk.frame3_idx[i];
        const int frame4_atom  = vsk.frame4_idx[i];
        Tcalc p_f2[3], p_f3[3], p_f4[3], pf3_m_pf2[3], pf4_m_pf2[3], p_vs[3];
        if (tcoord_is_sgnint) {
          p_f2[0] = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f2[1] = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f2[2] = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          p_f3[0] = static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f3[1] = static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f3[2] = static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          p_f4[0] = static_cast<Tcalc>(xcrd[frame4_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f4[1] = static_cast<Tcalc>(ycrd[frame4_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f4[2] = static_cast<Tcalc>(zcrd[frame4_atom] - zcrd[parent_atom]) * inv_gpos_factor;
        }
        else {
          p_f2[0] = xcrd[frame2_atom] - xcrd[parent_atom];
          p_f2[1] = ycrd[frame2_atom] - ycrd[parent_atom];
          p_f2[2] = zcrd[frame2_atom] - zcrd[parent_atom];
          p_f3[0] = xcrd[frame3_atom] - xcrd[parent_atom];
          p_f3[1] = ycrd[frame3_atom] - ycrd[parent_atom];
          p_f3[2] = zcrd[frame3_atom] - zcrd[parent_atom];
          p_f4[0] = xcrd[frame4_atom] - xcrd[parent_atom];
          p_f4[1] = ycrd[frame4_atom] - ycrd[parent_atom];
          p_f4[2] = zcrd[frame4_atom] - zcrd[parent_atom];
        }
        imageCoordinates<Tcalc, Tcalc>(&p_f2[0], &p_f2[1], &p_f2[2], umat, invu, unit_cell,
                                       minimum_image);
        imageCoordinates<Tcalc, Tcalc>(&p_f3[0], &p_f3[1], &p_f3[2], umat, invu, unit_cell,
                                       minimum_image);
        imageCoordinates<Tcalc, Tcalc>(&p_f4[0], &p_f4[1], &p_f4[2], umat, invu, unit_cell,
                                       minimum_image);
        const Tcalc pf3_factor = vsk.dim1[param_idx];
        const Tcalc pf4_factor = vsk.dim2[param_idx];
        pf3_m_pf2[0] = (pf3_factor * p_f3[0]) - p_f2[0];
        pf3_m_pf2[1] = (pf3_factor * p_f3[1]) - p_f2[1];
        pf3_m_pf2[2] = (pf3_factor * p_f3[2]) - p_f2[2];
        pf4_m_pf2[0] = (pf4_factor * p_f4[0]) - p_f2[0];
        pf4_m_pf2[1] = (pf4_factor * p_f4[1]) - p_f2[1];
        pf4_m_pf2[2] = (pf4_factor * p_f4[2]) - p_f2[2];
        crossProduct(pf3_m_pf2, pf4_m_pf2, p_vs);
        const Tcalc pvs_factor = (tcalc_is_double) ?
                                 vsk.dim3[param_idx] /
                                 sqrt((p_vs[0] * p_vs[0]) + (p_vs[1] * p_vs[1]) +
                                      (p_vs[2] * p_vs[2])) :
                                 vsk.dim3[param_idx] /
                                 sqrtf((p_vs[0] * p_vs[0]) + (p_vs[1] * p_vs[1]) +
                                       (p_vs[2] * p_vs[2]));
        if (tcoord_is_sgnint) {
          const Tcalc disp_mult = pvs_factor * gpos_scale_factor;
          xcrd[vsite_atom] = xcrd[parent_atom] + llround(disp_mult * p_vs[0]);
          ycrd[vsite_atom] = ycrd[parent_atom] + llround(disp_mult * p_vs[1]);
          zcrd[vsite_atom] = zcrd[parent_atom] + llround(disp_mult * p_vs[2]);
        }
        else {
          xcrd[vsite_atom] = xcrd[parent_atom] + (pvs_factor * p_vs[0]);
          ycrd[vsite_atom] = ycrd[parent_atom] + (pvs_factor * p_vs[1]);
          zcrd[vsite_atom] = zcrd[parent_atom] + (pvs_factor * p_vs[2]);
        }
      }
      break;
    case VirtualSiteKind::NONE:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void placeVirtualSites(PhaseSpaceWriter psw, const VirtualSiteKit<Tcalc> vsk) {
  placeVirtualSites(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, vsk);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void placeVirtualSites(CoordinateFrameWriter cfw, const VirtualSiteKit<Tcalc> vsk) {
  placeVirtualSites(cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.umat, cfw.invu, cfw.unit_cell, vsk);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
void transmitVirtualSiteForces(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                               Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, const double* umat,
                               const double* invu, const UnitCellType unit_cell,
                               const VirtualSiteKit<Tcalc> &vsk,
                               const Tcalc gpos_scale_factor, const Tcalc force_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const bool tforce_is_sgnint = isSignedIntegralScalarType<Tforce>();
  const Tcalc value_one = 1.0;
  const Tcalc inv_gpos_factor = value_one / gpos_scale_factor;
  const Tcalc inv_force_factor = value_one / force_scale_factor;
  for (int i = 0; i < vsk.nsite; i++) {
    const int vsite_atom = vsk.vs_atoms[i];
    const int parent_atom  = vsk.frame1_idx[i];
    const int frame2_atom  = vsk.frame2_idx[i];
    const int pidx = vsk.vs_param_idx[i];
    const ImagingMethod minimum_image = ImagingMethod::MINIMUM_IMAGE;

    // Copy the force on the virtual site into its own vector.  This will be needed in most
    // subsequent force transmissions as a real-valued representation of the forces (scaling
    // factor for fixed precision removed).
    Tcalc vs_frc[3];
    if (tforce_is_sgnint) {
      vs_frc[0] = static_cast<Tcalc>(xfrc[vsite_atom]) * inv_force_factor;
      vs_frc[1] = static_cast<Tcalc>(yfrc[vsite_atom]) * inv_force_factor;
      vs_frc[2] = static_cast<Tcalc>(zfrc[vsite_atom]) * inv_force_factor;
    }
    else {
      vs_frc[0] = xfrc[vsite_atom];
      vs_frc[1] = yfrc[vsite_atom];
      vs_frc[2] = zfrc[vsite_atom];
    }
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[pidx])) {
    case VirtualSiteKind::FLEX_2:
      {
        const Tcalc p_f2_factor = vsk.dim1[pidx];
        if (tforce_is_sgnint) {
          const Tcalc part_mult = (value_one - p_f2_factor);
          const Tforce xpart = llround(part_mult * static_cast<Tcalc>(xfrc[vsite_atom]));
          const Tforce ypart = llround(part_mult * static_cast<Tcalc>(yfrc[vsite_atom]));
          const Tforce zpart = llround(part_mult * static_cast<Tcalc>(zfrc[vsite_atom]));
          xfrc[parent_atom] += xpart;
          yfrc[parent_atom] += ypart;
          zfrc[parent_atom] += zpart;
          xfrc[frame2_atom] += xfrc[vsite_atom] - xpart;
          yfrc[frame2_atom] += yfrc[vsite_atom] - ypart;
          zfrc[frame2_atom] += zfrc[vsite_atom] - zpart;
        }
        else {
          xfrc[parent_atom] += (1.0 - p_f2_factor) * xfrc[vsite_atom];
          yfrc[parent_atom] += (1.0 - p_f2_factor) * yfrc[vsite_atom];
          zfrc[parent_atom] += (1.0 - p_f2_factor) * zfrc[vsite_atom];
          xfrc[frame2_atom] += p_f2_factor * xfrc[vsite_atom];
          yfrc[frame2_atom] += p_f2_factor * yfrc[vsite_atom];
          zfrc[frame2_atom] += p_f2_factor * zfrc[vsite_atom];
        }
      }
      break;
    case VirtualSiteKind::FIXED_2:
      {
        Tcalc p_f2[3], p_vs[3], vs_frc_proj[3], force_partition[3];
        if (tcoord_is_sgnint) {
          p_f2[0] = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f2[1] = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f2[2] = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          p_vs[0] = static_cast<Tcalc>(xcrd[vsite_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_vs[1] = static_cast<Tcalc>(ycrd[vsite_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_vs[2] = static_cast<Tcalc>(zcrd[vsite_atom] - zcrd[parent_atom]) * inv_gpos_factor;
        }
        else {
          p_f2[0] = xcrd[frame2_atom] - xcrd[parent_atom];
          p_f2[1] = ycrd[frame2_atom] - ycrd[parent_atom];
          p_f2[2] = zcrd[frame2_atom] - zcrd[parent_atom];
          p_vs[0] = xcrd[vsite_atom] - xcrd[parent_atom];
          p_vs[1] = ycrd[vsite_atom] - ycrd[parent_atom];
          p_vs[2] = zcrd[vsite_atom] - zcrd[parent_atom];
        }
        imageCoordinates<Tcalc, Tcalc>(&p_f2[0], &p_f2[1], &p_f2[2], umat, invu, unit_cell,
                                       minimum_image);
        imageCoordinates<Tcalc, Tcalc>(&p_vs[0], &p_vs[1], &p_vs[2], umat, invu, unit_cell,
                                       minimum_image);
        const Tcalc invr_p_f2 = (tcalc_is_double) ?
                                1.0 / sqrt((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                           (p_f2[2] * p_f2[2])) :
                                value_one / sqrtf((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                                  (p_f2[2] * p_f2[2]));

        // Take the projection of the force on the virtual site onto the parent atom -> virtual
        // site displacement vector p_vs.  The force is partitioned between the two frame atoms
        // according to the distance between the parent atom and the virtual site (normalized by
        // the distance between the two frame atoms).
        project(vs_frc, p_vs, vs_frc_proj, 3);
        const Tcalc p_vs_distance = vsk.dim1[pidx];
        for (int j = 0; j < 3; j++) {
          force_partition[j] = invr_p_f2 * p_vs_distance * (vs_frc[j] - vs_frc_proj[j]);
        }
        if (tforce_is_sgnint) {
          const Tforce xpart = llround(force_partition[0] * force_scale_factor);
          const Tforce ypart = llround(force_partition[1] * force_scale_factor);
          const Tforce zpart = llround(force_partition[2] * force_scale_factor);
          xfrc[parent_atom] += xfrc[vsite_atom] - xpart;
          yfrc[parent_atom] += yfrc[vsite_atom] - ypart;
          zfrc[parent_atom] += zfrc[vsite_atom] - zpart;
          xfrc[frame2_atom] += xpart;
          yfrc[frame2_atom] += ypart;
          zfrc[frame2_atom] += zpart;
        }
        else {
          xfrc[parent_atom] += vs_frc[0] - force_partition[0];
          yfrc[parent_atom] += vs_frc[1] - force_partition[1];
          zfrc[parent_atom] += vs_frc[2] - force_partition[2];
          xfrc[frame2_atom] += force_partition[0];
          yfrc[frame2_atom] += force_partition[1];
          zfrc[frame2_atom] += force_partition[2];
        }
      }
      break;
    case VirtualSiteKind::FLEX_3:
      {
        const int frame3_atom  = vsk.frame3_idx[i];
        if (tforce_is_sgnint) {
          const Tcalc p_f2_factor = vsk.dim1[pidx] * force_scale_factor;
          const Tcalc p_f3_factor = vsk.dim2[pidx] * force_scale_factor;
          const Tforce f2x_part = llround(p_f2_factor * vs_frc[0]);
          const Tforce f2y_part = llround(p_f2_factor * vs_frc[1]);
          const Tforce f2z_part = llround(p_f2_factor * vs_frc[2]);
          const Tforce f3x_part = llround(p_f3_factor * vs_frc[0]);
          const Tforce f3y_part = llround(p_f3_factor * vs_frc[1]);
          const Tforce f3z_part = llround(p_f3_factor * vs_frc[2]);
          xfrc[parent_atom] += xfrc[vsite_atom] - f2x_part - f3x_part;
          yfrc[parent_atom] += yfrc[vsite_atom] - f2y_part - f3y_part;
          zfrc[parent_atom] += zfrc[vsite_atom] - f2z_part - f3z_part;
          xfrc[frame2_atom] += f2x_part;
          yfrc[frame2_atom] += f2y_part;
          zfrc[frame2_atom] += f2z_part;
          xfrc[frame3_atom] += f3x_part;
          yfrc[frame3_atom] += f3y_part;
          zfrc[frame3_atom] += f3z_part;
        }
        else {
          const Tcalc p_f2_factor = vsk.dim1[pidx];
          const Tcalc p_f3_factor = vsk.dim2[pidx];
          xfrc[parent_atom] += (1.0 - p_f2_factor - p_f3_factor) * xfrc[vsite_atom];
          yfrc[parent_atom] += (1.0 - p_f2_factor - p_f3_factor) * yfrc[vsite_atom];
          zfrc[parent_atom] += (1.0 - p_f2_factor - p_f3_factor) * zfrc[vsite_atom];
          xfrc[frame2_atom] += p_f2_factor * xfrc[vsite_atom];
          yfrc[frame2_atom] += p_f2_factor * yfrc[vsite_atom];
          zfrc[frame2_atom] += p_f2_factor * zfrc[vsite_atom];
          xfrc[frame3_atom] += p_f3_factor * xfrc[vsite_atom];
          yfrc[frame3_atom] += p_f3_factor * yfrc[vsite_atom];
          zfrc[frame3_atom] += p_f3_factor * zfrc[vsite_atom];
        }
      }
      break;
    case VirtualSiteKind::FIXED_3:
      {
        Tcalc f2_f3[3], p_vs[3], p_mid[3], vs_frc_proj[3], force_partition[3];
        const int frame3_atom  = vsk.frame3_idx[i];
        const Tcalc f2_f3_factor  = vsk.dim2[pidx];
        if (tcoord_is_sgnint) {
          f2_f3[0] = static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[frame2_atom]) * inv_gpos_factor;
          f2_f3[1] = static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[frame2_atom]) * inv_gpos_factor;
          f2_f3[2] = static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[frame2_atom]) * inv_gpos_factor;
          p_vs[0]  = static_cast<Tcalc>(xcrd[vsite_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_vs[1]  = static_cast<Tcalc>(ycrd[vsite_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_vs[2]  = static_cast<Tcalc>(zcrd[vsite_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          p_mid[0] = (static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) *
                      inv_gpos_factor) + (f2_f3_factor * f2_f3[0]);
          p_mid[1] = (static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) *
                      inv_gpos_factor) + (f2_f3_factor * f2_f3[1]);
          p_mid[2] = (static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) *
                      inv_gpos_factor) + (f2_f3_factor * f2_f3[2]);
        }
        else {
          f2_f3[0] = xcrd[frame3_atom] - xcrd[frame2_atom];
          f2_f3[1] = ycrd[frame3_atom] - ycrd[frame2_atom];
          f2_f3[2] = zcrd[frame3_atom] - zcrd[frame2_atom];
          p_vs[0] = xcrd[vsite_atom] - xcrd[parent_atom];
          p_vs[1] = ycrd[vsite_atom] - ycrd[parent_atom];
          p_vs[2] = zcrd[vsite_atom] - zcrd[parent_atom];
          p_mid[0] = xcrd[frame2_atom] - xcrd[parent_atom] + (f2_f3_factor * f2_f3[0]);
          p_mid[1] = ycrd[frame2_atom] - ycrd[parent_atom] + (f2_f3_factor * f2_f3[1]);
          p_mid[2] = zcrd[frame2_atom] - zcrd[parent_atom] + (f2_f3_factor * f2_f3[2]);
        }
        imageCoordinates<Tcalc, Tcalc>(&f2_f3[0], &f2_f3[1], &f2_f3[2], umat, invu, unit_cell,
                                       minimum_image);
        imageCoordinates<Tcalc, Tcalc>(&p_vs[0], &p_vs[1], &p_vs[2], umat, invu, unit_cell,
                                       minimum_image);
        imageCoordinates<Tcalc, Tcalc>(&p_mid[0], &p_mid[1], &p_mid[2], umat, invu, unit_cell,
                                       minimum_image);

        // As with the fixed distance, two-point frame, compute the part of the force on the
        // virtual site that is perpendicular to the parent atom -> virtual site displacement
        // vector and use this to compute the partitioning between the parent atom and some
        // imaginary midpoint on the line between frame atom 2 and frame atom 3.  That's like
        // the FIXED_2 frame type.  The force on the midpoint is subsequently distributed, akin
        // to the FLEX_2 type, between atoms 2 and 3.
        project(vs_frc, p_vs, vs_frc_proj, 3);
        const Tcalc p_vs_distance = vsk.dim1[pidx];
        const Tcalc invr_p_mid = 1.0 / sqrt((p_mid[0] * p_mid[0]) + (p_mid[1] * p_mid[1]) +
                                            (p_mid[2] * p_mid[2]));
        for (int j = 0; j < 3; j++) {
          force_partition[j] = invr_p_mid * p_vs_distance * (vs_frc[j] - vs_frc_proj[j]);
        }
        if (tforce_is_sgnint) {
          const Tforce xpart = llround(force_partition[0] * force_scale_factor);
          const Tforce ypart = llround(force_partition[1] * force_scale_factor);
          const Tforce zpart = llround(force_partition[2] * force_scale_factor);
          const Tforce f2f3_xprt = llround(f2_f3_factor * force_partition[0] * force_scale_factor);
          const Tforce f2f3_yprt = llround(f2_f3_factor * force_partition[1] * force_scale_factor);
          const Tforce f2f3_zprt = llround(f2_f3_factor * force_partition[2] * force_scale_factor);
          xfrc[parent_atom] += xfrc[vsite_atom] - xpart;
          yfrc[parent_atom] += yfrc[vsite_atom] - ypart;
          zfrc[parent_atom] += zfrc[vsite_atom] - zpart;
          xfrc[frame2_atom] += xpart - f2f3_xprt;
          yfrc[frame2_atom] += ypart - f2f3_yprt;
          zfrc[frame2_atom] += zpart - f2f3_zprt;
          xfrc[frame3_atom] += f2f3_xprt;
          yfrc[frame3_atom] += f2f3_yprt;
          zfrc[frame3_atom] += f2f3_zprt;
        }
        else {
          xfrc[parent_atom] += vs_frc[0] - force_partition[0];
          yfrc[parent_atom] += vs_frc[1] - force_partition[1];
          zfrc[parent_atom] += vs_frc[2] - force_partition[2];
          xfrc[frame2_atom] += (1.0 - f2_f3_factor) * force_partition[0];
          yfrc[frame2_atom] += (1.0 - f2_f3_factor) * force_partition[1];
          zfrc[frame2_atom] += (1.0 - f2_f3_factor) * force_partition[2];
          xfrc[frame3_atom] += f2_f3_factor * force_partition[0];
          yfrc[frame3_atom] += f2_f3_factor * force_partition[1];
          zfrc[frame3_atom] += f2_f3_factor * force_partition[2];
        }
      }
      break;
    case VirtualSiteKind::FAD_3:
      {
        Tcalc p_f2[3], f2_f3[3], f23_t_pf2[3], F1[3], F2[3], F3[3];
        const int frame3_atom  = vsk.frame3_idx[i];
        if (tcoord_is_sgnint) {
          p_f2[0]  = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f2[1]  = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f2[2]  = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          f2_f3[0] = static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[frame2_atom]) * inv_gpos_factor;
          f2_f3[1] = static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[frame2_atom]) * inv_gpos_factor;
          f2_f3[2] = static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[frame2_atom]) * inv_gpos_factor;
        }
        else {
          p_f2[0]  = xcrd[frame2_atom] - xcrd[parent_atom];
          p_f2[1]  = ycrd[frame2_atom] - ycrd[parent_atom];
          p_f2[2]  = zcrd[frame2_atom] - zcrd[parent_atom];
          f2_f3[0] = xcrd[frame3_atom] - xcrd[frame2_atom];
          f2_f3[1] = ycrd[frame3_atom] - ycrd[frame2_atom];
          f2_f3[2] = zcrd[frame3_atom] - zcrd[frame2_atom];
        }
        imageCoordinates<Tcalc, Tcalc>( &p_f2[0],  &p_f2[1],  &p_f2[2], umat, invu, unit_cell,
                                       minimum_image);
        imageCoordinates<Tcalc, Tcalc>(&f2_f3[0], &f2_f3[1], &f2_f3[2], umat, invu, unit_cell,
                                       minimum_image);
        project(f2_f3, p_f2, f23_t_pf2, 3);
        f23_t_pf2[0] = f2_f3[0] - f23_t_pf2[0];
        f23_t_pf2[1] = f2_f3[1] - f23_t_pf2[1];
        f23_t_pf2[2] = f2_f3[2] - f23_t_pf2[2];
        const Tcalc invr2_p_f2 = value_one / ((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                              (p_f2[2] * p_f2[2]));
        const Tcalc invr2_t = value_one / ((f23_t_pf2[0] * f23_t_pf2[0]) +
                                           (f23_t_pf2[1] * f23_t_pf2[1]) +
                                           (f23_t_pf2[2] * f23_t_pf2[2]));
        Tcalc invr_p_f2, invr_t, p_f2_factor, t_factor;
        if (tcalc_is_double) {
          invr_p_f2 = sqrt(invr2_p_f2);
          invr_t = sqrt(invr2_t);
          p_f2_factor = vsk.dim1[pidx] * cos(vsk.dim2[pidx]) * invr_p_f2;
          t_factor    = vsk.dim1[pidx] * sin(vsk.dim2[pidx]) * invr_t;
        }
        else {
          invr_p_f2 = sqrtf(invr2_p_f2);
          invr_t = sqrtf(invr2_t);
          p_f2_factor = vsk.dim1[pidx] * cosf(vsk.dim2[pidx]) * invr_p_f2;
          t_factor    = vsk.dim1[pidx] * sinf(vsk.dim2[pidx]) * invr_t;
        }
        const Tcalc f1fac  = dot(p_f2, vs_frc, 3) * invr2_p_f2;
        const Tcalc f2fac  = dot(f23_t_pf2, vs_frc, 3) * invr2_t;
        const Tcalc abbcOabab = dot(p_f2, f2_f3, 3) * invr2_p_f2;
        for (int j = 0; j < 3; j++) {
          F1[j] = vs_frc[j] - (f1fac * p_f2[j]);
          F2[j] = F1[j] - (f2fac * f23_t_pf2[j]);
          F3[j] = f1fac * f23_t_pf2[j];          
        }
        if (tforce_is_sgnint) {
          const Tforce f1_xpart = llround(((p_f2_factor * F1[0]) -
                                           (t_factor * ((abbcOabab * F2[0]) + F3[0]))) *
                                          force_scale_factor);
          const Tforce f1_ypart = llround(((p_f2_factor * F1[1]) -
                                           (t_factor * ((abbcOabab * F2[1]) + F3[1]))) *
                                          force_scale_factor);
          const Tforce f1_zpart = llround(((p_f2_factor * F1[2]) -
                                           (t_factor * ((abbcOabab * F2[2]) + F3[2]))) *
                                          force_scale_factor);
          const Tforce f3_xpart = llround(t_factor * F2[0] * force_scale_factor);
          const Tforce f3_ypart = llround(t_factor * F2[1] * force_scale_factor);
          const Tforce f3_zpart = llround(t_factor * F2[2] * force_scale_factor);
          xfrc[parent_atom] += xfrc[vsite_atom] - f1_xpart;
          yfrc[parent_atom] += yfrc[vsite_atom] - f1_ypart;
          zfrc[parent_atom] += zfrc[vsite_atom] - f1_zpart;
          xfrc[frame2_atom] += f1_xpart - f3_xpart;
          yfrc[frame2_atom] += f1_ypart - f3_ypart;
          zfrc[frame2_atom] += f1_zpart - f3_zpart;
          xfrc[frame3_atom] += f3_xpart;
          yfrc[frame3_atom] += f3_ypart;
          zfrc[frame3_atom] += f3_zpart;
        }
        else {
          xfrc[parent_atom] += vs_frc[0] - (p_f2_factor * F1[0]) +
                               (t_factor * ((abbcOabab * F2[0]) + F3[0]));
          yfrc[parent_atom] += vs_frc[1] - (p_f2_factor * F1[1]) +
                               (t_factor * ((abbcOabab * F2[1]) + F3[1]));
          zfrc[parent_atom] += vs_frc[2] - (p_f2_factor * F1[2]) +
                               (t_factor * ((abbcOabab * F2[2]) + F3[2]));
          xfrc[frame2_atom] += (p_f2_factor * F1[0]) -
                               (t_factor * (F2[0] + (abbcOabab * F2[0]) + F3[0]));
          yfrc[frame2_atom] += (p_f2_factor * F1[1]) -
                               (t_factor * (F2[1] + (abbcOabab * F2[1]) + F3[1]));
          zfrc[frame2_atom] += (p_f2_factor * F1[2]) -
                               (t_factor * (F2[2] + (abbcOabab * F2[2]) + F3[2]));
          xfrc[frame3_atom] += t_factor * F2[0]; 
          yfrc[frame3_atom] += t_factor * F2[1];
          zfrc[frame3_atom] += t_factor * F2[2];
        }
      }
      break;
    case VirtualSiteKind::OUT_3:
      {
        const int frame3_atom  = vsk.frame3_idx[i];
        Tcalc p_f2[3], p_f3[3], partition_f2[3], partition_f3[3];
        Tcalc mf2_01, mf2_02, mf2_12, mf3_01, mf3_02, mf3_12;
        if (tcoord_is_sgnint) {
          const Tcalc d3_factor = vsk.dim3[pidx] * inv_gpos_factor;
          mf2_01 = d3_factor * static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[parent_atom]);
          mf2_02 = d3_factor * static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[parent_atom]);
          mf2_12 = d3_factor * static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[parent_atom]);
          mf3_01 = d3_factor * static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]);
          mf3_02 = d3_factor * static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]);
          mf3_12 = d3_factor * static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]);
        }
        else {
          mf2_01 = vsk.dim3[pidx] * (zcrd[frame3_atom] - zcrd[parent_atom]);
          mf2_02 = vsk.dim3[pidx] * (ycrd[frame3_atom] - ycrd[parent_atom]);
          mf2_12 = vsk.dim3[pidx] * (xcrd[frame3_atom] - xcrd[parent_atom]);
          mf3_01 = vsk.dim3[pidx] * (zcrd[frame2_atom] - zcrd[parent_atom]);
          mf3_02 = vsk.dim3[pidx] * (ycrd[frame2_atom] - ycrd[parent_atom]);
          mf3_12 = vsk.dim3[pidx] * (xcrd[frame2_atom] - xcrd[parent_atom]);
        }
        partition_f2[0] = ( vsk.dim1[pidx] * vs_frc[0]) - (mf2_01 * vs_frc[1]) +
                          (mf2_02 * vs_frc[2]);
        partition_f2[1] = ( mf2_01 * vs_frc[0]) + (vsk.dim1[pidx] * vs_frc[1]) -
                          (mf2_12 * vs_frc[2]);
        partition_f2[2] = (-mf2_02 * vs_frc[0]) + (mf2_12 * vs_frc[1]) +
                          (vsk.dim1[pidx] * vs_frc[2]);
        partition_f3[0] = ( vsk.dim2[pidx] * vs_frc[0]) + (mf3_01 * vs_frc[1]) -
                          (mf3_02 * vs_frc[2]);
        partition_f3[1] = (-mf3_01 * vs_frc[0]) + (vsk.dim2[pidx] * vs_frc[1]) +
                          (mf3_12 * vs_frc[2]);
        partition_f3[2] = ( mf3_02 * vs_frc[0]) - (mf3_12 * vs_frc[1]) +
                          (vsk.dim2[pidx] * vs_frc[2]);
        if (tforce_is_sgnint) {
          const Tforce f2_xpart = llround(partition_f2[0] * force_scale_factor);
          const Tforce f2_ypart = llround(partition_f2[1] * force_scale_factor);
          const Tforce f2_zpart = llround(partition_f2[2] * force_scale_factor);
          const Tforce f3_xpart = llround(partition_f3[0] * force_scale_factor);
          const Tforce f3_ypart = llround(partition_f3[1] * force_scale_factor);
          const Tforce f3_zpart = llround(partition_f3[2] * force_scale_factor);
          xfrc[parent_atom] += xfrc[vsite_atom] - f2_xpart - f3_xpart;
          yfrc[parent_atom] += yfrc[vsite_atom] - f2_ypart - f3_ypart;
          zfrc[parent_atom] += zfrc[vsite_atom] - f2_zpart - f3_zpart;
          xfrc[frame2_atom] += f2_xpart;
          yfrc[frame2_atom] += f2_ypart;
          zfrc[frame2_atom] += f2_zpart;
          xfrc[frame3_atom] += f3_xpart;
          yfrc[frame3_atom] += f3_ypart;
          zfrc[frame3_atom] += f3_zpart;
        }
        else {
          xfrc[parent_atom] += vs_frc[0] - partition_f2[0] - partition_f3[0];
          yfrc[parent_atom] += vs_frc[1] - partition_f2[1] - partition_f3[1];
          zfrc[parent_atom] += vs_frc[2] - partition_f2[2] - partition_f3[2];
          xfrc[frame2_atom] += partition_f2[0];
          yfrc[frame2_atom] += partition_f2[1];
          zfrc[frame2_atom] += partition_f2[2];
          xfrc[frame3_atom] += partition_f3[0];
          yfrc[frame3_atom] += partition_f3[1];
          zfrc[frame3_atom] += partition_f3[2];
        }
      }
      break;
    case VirtualSiteKind::FIXED_4:
      {
        const int frame3_atom  = vsk.frame3_idx[i];
        const int frame4_atom  = vsk.frame4_idx[i];
        Tcalc p_f2[3], rj_f3[3], rj_f4[3], rj_f34[3], rm[3], rt[3], fb[3], fc[3], fd[3];
        if (tcoord_is_sgnint) {
          p_f2[0] = static_cast<Tcalc>(xcrd[frame2_atom] - xcrd[parent_atom]) * inv_gpos_factor;
          p_f2[1] = static_cast<Tcalc>(ycrd[frame2_atom] - ycrd[parent_atom]) * inv_gpos_factor;
          p_f2[2] = static_cast<Tcalc>(zcrd[frame2_atom] - zcrd[parent_atom]) * inv_gpos_factor;
          const Tcalc d1_factor = vsk.dim1[pidx] * inv_gpos_factor;
          const Tcalc d2_factor = vsk.dim2[pidx] * inv_gpos_factor;
          rj_f3[0] = (d1_factor * static_cast<Tcalc>(xcrd[frame3_atom] - xcrd[parent_atom])) -
                     p_f2[0];
          rj_f3[1] = (d1_factor * static_cast<Tcalc>(ycrd[frame3_atom] - ycrd[parent_atom])) -
                     p_f2[1];
          rj_f3[2] = (d1_factor * static_cast<Tcalc>(zcrd[frame3_atom] - zcrd[parent_atom])) -
                     p_f2[2];
          rj_f4[0] = (d2_factor * static_cast<Tcalc>(xcrd[frame4_atom] - xcrd[parent_atom])) -
                     p_f2[0];
          rj_f4[1] = (d2_factor * static_cast<Tcalc>(ycrd[frame4_atom] - ycrd[parent_atom])) -
                     p_f2[1];
          rj_f4[2] = (d2_factor * static_cast<Tcalc>(zcrd[frame4_atom] - zcrd[parent_atom])) -
                     p_f2[2];
        }
        else {
          p_f2[0] = xcrd[frame2_atom] - xcrd[parent_atom];
          p_f2[1] = ycrd[frame2_atom] - ycrd[parent_atom];
          p_f2[2] = zcrd[frame2_atom] - zcrd[parent_atom];
          rj_f3[0] = (vsk.dim1[pidx] * (xcrd[frame3_atom] - xcrd[parent_atom])) - p_f2[0];
          rj_f3[1] = (vsk.dim1[pidx] * (ycrd[frame3_atom] - ycrd[parent_atom])) - p_f2[1];
          rj_f3[2] = (vsk.dim1[pidx] * (zcrd[frame3_atom] - zcrd[parent_atom])) - p_f2[2];
          rj_f4[0] = (vsk.dim2[pidx] * (xcrd[frame4_atom] - xcrd[parent_atom])) - p_f2[0];
          rj_f4[1] = (vsk.dim2[pidx] * (ycrd[frame4_atom] - ycrd[parent_atom])) - p_f2[1];
          rj_f4[2] = (vsk.dim2[pidx] * (zcrd[frame4_atom] - zcrd[parent_atom])) - p_f2[2];
        }
        rj_f34[0] = rj_f4[0] - rj_f3[0];
        rj_f34[1] = rj_f4[1] - rj_f3[1];
        rj_f34[2] = rj_f4[2] - rj_f3[2];
        crossProduct(rj_f3, rj_f4, rm);
        const Tcalc invr2_rm = (tcalc_is_double) ?
                               1.0 / ((rm[0] * rm[0]) + (rm[1] * rm[1]) + (rm[2] * rm[2])) :
                               value_one / ((rm[0] * rm[0]) + (rm[1] * rm[1]) + (rm[2] * rm[2]));
        const Tcalc invr_rm  = (tcalc_is_double) ? sqrt(invr2_rm) : sqrtf(invr2_rm);
        const Tcalc cfx = vsk.dim3[pidx] * invr_rm * vs_frc[0];
        const Tcalc cfy = vsk.dim3[pidx] * invr_rm * vs_frc[1];
        const Tcalc cfz = vsk.dim3[pidx] * invr_rm * vs_frc[2];
        crossProduct(rm, rj_f34, rt);
        rt[0] *= invr2_rm;
        rt[1] *= invr2_rm;
        rt[2] *= invr2_rm;
        fb[0] = (-rm[0] * rt[0] * cfx) + ((rj_f34[2] - (rm[1] * rt[0])) * cfy) -
                ( (rj_f34[1] + (rm[2] * rt[0])) * cfz);
        fb[1] = (-(rj_f34[2] + (rm[0] * rt[1])) * cfx) - (rm[1] * rt[1] * cfy) +
                ( (rj_f34[0] - (rm[2] * rt[1])) * cfz);
        fb[2] = ((rj_f34[1] - (rm[0] * rt[2])) * cfx) - ((rj_f34[0] + (rm[1] * rt[2])) * cfy) -
                (rm[2] * rt[2] * cfz);
        rt[0] = ((rj_f4[1] * rm[2]) - (rj_f4[2] * rm[1])) * invr2_rm * vsk.dim1[pidx];
        rt[1] = ((rj_f4[2] * rm[0]) - (rj_f4[0] * rm[2])) * invr2_rm * vsk.dim1[pidx];
        rt[2] = ((rj_f4[0] * rm[1]) - (rj_f4[1] * rm[0])) * invr2_rm * vsk.dim1[pidx];
        fc[0] = (-rm[0] * rt[0] * cfx) - (((vsk.dim1[pidx] * rj_f4[2]) + (rm[1] * rt[0])) * cfy) +
                (((vsk.dim1[pidx] * rj_f4[1]) - (rm[2] * rt[0])) * cfz);
        fc[1] = (((vsk.dim1[pidx] * rj_f4[2]) - (rm[0] * rt[1])) * cfx) - (rm[1] * rt[1] * cfy) -
                (((vsk.dim1[pidx] * rj_f4[0]) + (rm[2] * rt[1])) * cfz);
        fc[2] = (-((vsk.dim1[pidx] * rj_f4[1]) + (rm[0] * rt[2])) * cfx) +
                (((vsk.dim1[pidx] * rj_f4[0]) - (rm[1] * rt[2])) * cfy) - (rm[2] * rt[2] * cfz);
        rt[0] = ((rm[1] * rj_f3[2]) - (rm[2] * rj_f3[1])) * invr2_rm * vsk.dim2[pidx];
        rt[1] = ((rm[2] * rj_f3[0]) - (rm[0] * rj_f3[2])) * invr2_rm * vsk.dim2[pidx];
        rt[2] = ((rm[0] * rj_f3[1]) - (rm[1] * rj_f3[0])) * invr2_rm * vsk.dim2[pidx];
        fd[0] = (-rm[0] * rt[0] * cfx) + (((vsk.dim2[pidx] * rj_f3[2]) - (rm[1] * rt[0])) * cfy) -
                (((vsk.dim2[pidx] * rj_f3[1]) + (rm[2] * rt[0])) * cfz);
        fd[1] = (-((vsk.dim2[pidx] * rj_f3[2]) + (rm[0] * rt[1])) * cfx) - (rm[1] * rt[1] * cfy) +
                (((vsk.dim2[pidx] * rj_f3[0]) - (rm[2] * rt[1])) * cfz);
        fd[2] = (((vsk.dim2[pidx] * rj_f3[1]) - (rm[0] * rt[2])) * cfx) +
                (-((vsk.dim2[pidx] * rj_f3[0]) + (rm[1] * rt[2])) * cfy) - (rm[2] * rt[2] * cfz);
        if (tforce_is_sgnint) {
          const Tforce fb_xpart = llround(fb[0]);
          const Tforce fb_ypart = llround(fb[1]);
          const Tforce fb_zpart = llround(fb[2]);
          const Tforce fc_xpart = llround(fc[0]);
          const Tforce fc_ypart = llround(fc[1]);
          const Tforce fc_zpart = llround(fc[2]);
          const Tforce fd_xpart = llround(fd[0]);
          const Tforce fd_ypart = llround(fd[1]);
          const Tforce fd_zpart = llround(fd[2]);
          xfrc[parent_atom] += xfrc[vsite_atom] - fb_xpart - fc_xpart - fd_xpart;
          yfrc[parent_atom] += yfrc[vsite_atom] - fb_ypart - fc_ypart - fd_ypart;
          zfrc[parent_atom] += zfrc[vsite_atom] - fb_zpart - fc_zpart - fd_zpart;
          xfrc[frame2_atom] += fb_xpart;
          yfrc[frame2_atom] += fb_ypart;
          zfrc[frame2_atom] += fb_zpart;
          xfrc[frame3_atom] += fc_xpart;
          yfrc[frame3_atom] += fc_ypart;
          zfrc[frame3_atom] += fc_zpart;
          xfrc[frame4_atom] += fd_xpart;
          yfrc[frame4_atom] += fd_ypart;
          zfrc[frame4_atom] += fd_zpart;
        }
        else {
          xfrc[parent_atom] += vs_frc[0] - fb[0] - fc[0] - fd[0];
          yfrc[parent_atom] += vs_frc[1] - fb[1] - fc[1] - fd[1];
          zfrc[parent_atom] += vs_frc[2] - fb[2] - fc[2] - fd[2];
          xfrc[frame2_atom] += fb[0];
          yfrc[frame2_atom] += fb[1];
          zfrc[frame2_atom] += fb[2];
          xfrc[frame3_atom] += fc[0];
          yfrc[frame3_atom] += fc[1];
          zfrc[frame3_atom] += fc[2];
          xfrc[frame4_atom] += fd[0];
          yfrc[frame4_atom] += fd[1];
          zfrc[frame4_atom] += fd[2];
        }
      }
      break;
    case VirtualSiteKind::NONE:
      break;
    }
    if (tforce_is_sgnint) {
      xfrc[vsite_atom] = 0LL;
      yfrc[vsite_atom] = 0LL;
      zfrc[vsite_atom] = 0LL;
    }
    else {
      xfrc[vsite_atom] = 0.0;
      yfrc[vsite_atom] = 0.0;
      zfrc[vsite_atom] = 0.0;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void transmitVirtualSiteForces(PhaseSpaceWriter psw, const VirtualSiteKit<Tcalc> vsk) {
  transmitVirtualSiteForces(psw.xcrd, psw.ycrd, psw.zcrd, psw.xfrc, psw.yfrc, psw.zfrc,
                            psw.umat, psw.invu, psw.unit_cell, vsk);
}

} // namespace structure
} // namespace stormm
