// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace chemistry {

//-------------------------------------------------------------------------------------------------
template <typename T>
int getChiralOrientation(const T* xcrd, const T* ycrd, const T* zcrd, const int center_atom,
                         const int root_atom, const int branch_a_atom, const int branch_b_atom,
                         const int branch_c_atom) {

  // Collect coordinates
  T ra[3], rb[3], rc[3], r_root[3], prja[3], prjb[3], prjc[3], acrb[3], bcrc[3];
  const T cax = xcrd[center_atom];
  const T cay = ycrd[center_atom];
  const T caz = zcrd[center_atom];
  ra[0] = xcrd[branch_a_atom] - cax;
  ra[1] = ycrd[branch_a_atom] - cay;
  ra[2] = zcrd[branch_a_atom] - caz;
  rb[0] = xcrd[branch_b_atom] - cax;
  rb[1] = ycrd[branch_b_atom] - cay;
  rb[2] = zcrd[branch_b_atom] - caz;
  rc[0] = xcrd[branch_c_atom] - cax;
  rc[1] = ycrd[branch_c_atom] - cay;
  rc[2] = zcrd[branch_c_atom] - caz;
  r_root[0] = xcrd[root_atom] - cax;
  r_root[1] = ycrd[root_atom] - cay;
  r_root[2] = zcrd[root_atom] - caz;

  // Remove the projections of each higher priority branch onto the lowest priority branch.
  project(ra, r_root, prja, 3);
  project(rb, r_root, prjb, 3);
  project(rc, r_root, prjc, 3);
  for (int i = 0; i < 3; i++) {
    ra[i] -= prja[i];
    rb[i] -= prjb[i];
    rc[i] -= prjc[i];
  }
  crossProduct(ra, rb, acrb);
  crossProduct(rb, rc, bcrc);

  // L- (S-) chirality occurs when the cross product of ra and rb, as well as the cross product of
  // rb and rc, following the right hand rule, point towards r_root.  However, L- chirality is
  // given a convention of + (as most amino acids, which have chiral centers, are expected to be
  // L-chiral).
  return ((2 * (dot(acrb, r_root, 3) < 0.0 && dot(bcrc, r_root, 3))) - 1);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
int getChiralOrientation(const T* xcrd, const T* ycrd, const T* zcrd, const int* xcrd_ovrf,
                         const int* ycrd_ovrf, const int* zcrd_ovrf, const int center_atom,
                         const int root_atom, const int branch_a_atom, const int branch_b_atom,
                         const int branch_c_atom, const double inv_scale) {

  // Collect coordinates
  double ra[3], rb[3], rc[3], r_root[3], prja[3], prjb[3], prjc[3], acrb[3], bcrc[3];
  if (xcrd_ovrf != nullptr && ycrd_ovrf != nullptr && zcrd_ovrf != nullptr) {

    // While it is of hardly any consequence whether the full precision of the 95-bit type is
    // preserved, the most expensive operation is hostInt95ToDouble() and this approach, keeping
    // the full precision of the data type intact through subtraction, uses the fewest of them.
    const int95_t i_rax = hostInt95Sum(xcrd[branch_a_atom], xcrd_ovrf[branch_a_atom],
                                       -xcrd[center_atom], -xcrd_ovrf[center_atom]);
    const int95_t i_ray = hostInt95Sum(ycrd[branch_a_atom], ycrd_ovrf[branch_a_atom],
                                       -ycrd[center_atom], -ycrd_ovrf[center_atom]);
    const int95_t i_raz = hostInt95Sum(zcrd[branch_a_atom], zcrd_ovrf[branch_a_atom],
                                       -zcrd[center_atom], -zcrd_ovrf[center_atom]);
    const int95_t i_rbx = hostInt95Sum(xcrd[branch_b_atom], xcrd_ovrf[branch_b_atom],
                                       -xcrd[center_atom], -xcrd_ovrf[center_atom]);
    const int95_t i_rby = hostInt95Sum(ycrd[branch_b_atom], ycrd_ovrf[branch_b_atom],
                                       -ycrd[center_atom], -ycrd_ovrf[center_atom]);
    const int95_t i_rbz = hostInt95Sum(zcrd[branch_b_atom], zcrd_ovrf[branch_b_atom],
                                       -zcrd[center_atom], -zcrd_ovrf[center_atom]);
    const int95_t i_rcx = hostInt95Sum(xcrd[branch_c_atom], xcrd_ovrf[branch_c_atom],
                                       -xcrd[center_atom], -xcrd_ovrf[center_atom]);
    const int95_t i_rcy = hostInt95Sum(ycrd[branch_c_atom], ycrd_ovrf[branch_c_atom],
                                       -ycrd[center_atom], -ycrd_ovrf[center_atom]);
    const int95_t i_rcz = hostInt95Sum(zcrd[branch_c_atom], zcrd_ovrf[branch_c_atom],
                                       -zcrd[center_atom], -zcrd_ovrf[center_atom]);
    const int95_t i_rrx = hostInt95Sum(xcrd[root_atom], xcrd_ovrf[root_atom], -xcrd[center_atom],
                                       -xcrd_ovrf[center_atom]);
    const int95_t i_rry = hostInt95Sum(ycrd[root_atom], ycrd_ovrf[root_atom], -ycrd[center_atom],
                                       -ycrd_ovrf[center_atom]);
    const int95_t i_rrz = hostInt95Sum(zcrd[root_atom], zcrd_ovrf[root_atom], -zcrd[center_atom],
                                       -zcrd_ovrf[center_atom]);
    ra[0] = hostInt95ToDouble(i_rax) * inv_scale;
    ra[1] = hostInt95ToDouble(i_ray) * inv_scale;
    ra[2] = hostInt95ToDouble(i_raz) * inv_scale;
    rb[0] = hostInt95ToDouble(i_rbx) * inv_scale;
    rb[1] = hostInt95ToDouble(i_rby) * inv_scale;
    rb[2] = hostInt95ToDouble(i_rbz) * inv_scale;
    rc[0] = hostInt95ToDouble(i_rcx) * inv_scale;
    rc[1] = hostInt95ToDouble(i_rcy) * inv_scale;
    rc[2] = hostInt95ToDouble(i_rcz) * inv_scale;
    r_root[0] = hostInt95ToDouble(i_rrx) * inv_scale;
    r_root[1] = hostInt95ToDouble(i_rry) * inv_scale;
    r_root[2] = hostInt95ToDouble(i_rrz) * inv_scale;
  }
  else {
    const T cax = xcrd[center_atom];
    const T cay = ycrd[center_atom];
    const T caz = zcrd[center_atom];
    ra[0] = (xcrd[branch_a_atom] - cax) * inv_scale;
    ra[1] = (ycrd[branch_a_atom] - cay) * inv_scale;
    ra[2] = (zcrd[branch_a_atom] - caz) * inv_scale;
    rb[0] = (xcrd[branch_b_atom] - cax) * inv_scale;
    rb[1] = (ycrd[branch_b_atom] - cay) * inv_scale;
    rb[2] = (zcrd[branch_b_atom] - caz) * inv_scale;
    rc[0] = (xcrd[branch_c_atom] - cax) * inv_scale;
    rc[1] = (ycrd[branch_c_atom] - cay) * inv_scale;
    rc[2] = (zcrd[branch_c_atom] - caz) * inv_scale;
    r_root[0] = (xcrd[root_atom] - cax) * inv_scale;
    r_root[1] = (ycrd[root_atom] - cay) * inv_scale;
    r_root[2] = (zcrd[root_atom] - caz) * inv_scale;
  }

  // Remove the projections of each higher priority branch onto the lowest priority branch.
  project(ra, r_root, prja, 3);
  project(rb, r_root, prjb, 3);
  project(rc, r_root, prjc, 3);
  for (int i = 0; i < 3; i++) {
    ra[i] -= prja[i];
    rb[i] -= prjb[i];
    rc[i] -= prjc[i];
  }
  crossProduct(ra, rb, acrb);
  crossProduct(rb, rc, bcrc);

  // L- (S-) chirality occurs when the cross product of ra and rb, as well as the cross product of
  // rb and rc, following the right hand rule, point towards r_root.  However, L- chirality is
  // given a convention of + (as most amino acids, which have chiral centers, are expected to be
  // L-chiral).
  return ((2 * (dot(acrb, r_root, 3) < 0.0 && dot(bcrc, r_root, 3))) - 1);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
int getChiralOrientation(const CoordinateSeriesReader<T> &csr, const size_t frame_index,
                         const int center_atom, const int root_atom, const int branch_a_atom,
                         const int branch_b_atom, const int branch_c_atom) {
  const size_t offset = frame_index * static_cast<size_t>(roundUp(csr.natom, warp_size_int));
  if (isFloatingPointScalarType<T>()) {
    return getChiralOrientation<T>(&csr.xcrd[offset], &csr.ycrd[offset], &csr.zcrd[offset],
                                   center_atom, root_atom, branch_a_atom, branch_b_atom,
                                   branch_c_atom);
  }
  else {
    return getChiralOrientation<T>(&csr.xcrd[offset], &csr.ycrd[offset], &csr.zcrd[offset],
                                   nullptr, nullptr, nullptr, center_atom, root_atom,
                                   branch_a_atom, branch_b_atom, branch_c_atom,
                                   csr.inv_gpos_scale);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
int getChiralOrientation(const CoordinateSeries<T> *cs, const size_t frame_index,
                         const int center_atom, const int root_atom, const int branch_a_atom,
                         const int branch_b_atom, const int branch_c_atom) {
  return getChiralOrientation<T>(cs->data(), frame_index, center_atom, root_atom, branch_a_atom,
                                 branch_b_atom, branch_c_atom);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
int getChiralOrientation(const CoordinateSeries<T> &cs, const size_t frame_index,
                         const int center_atom, const int root_atom, const int branch_a_atom,
                         const int branch_b_atom, const int branch_c_atom) {
  return getChiralOrientation<T>(cs.data(), frame_index, center_atom, root_atom, branch_a_atom,
                                 branch_b_atom, branch_c_atom);
}

} // namespace chemistry
} // namespace stormm
