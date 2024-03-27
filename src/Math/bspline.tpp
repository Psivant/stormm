// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void bSplineInputChecks(const std::vector<Tcoord> &xcrd, const std::vector<Tcoord> &ycrd,
                        const std::vector<Tcoord> &zcrd, const int order,
                        std::vector<Tcalc> *a_coefs, std::vector<Tcalc> *b_coefs,
                        std::vector<Tcalc> *c_coefs, std::vector<int> *a_init,
                        std::vector<int> *b_init, std::vector<int> *c_init) {
  const size_t natom = xcrd.size();
  if (natom != ycrd.size() || natom != zcrd.size()) {
    rtErr("Input coordinate vectors must have similar lengths (" + std::to_string(xcrd.size()) +
          ", " + std::to_string(ycrd.size()) + ", and " + std::to_string(zcrd.size()) +
          " provided).", "bSpline");
  }
  const size_t order_zu = order;
  if (natom * order_zu != a_coefs->size() || natom * order_zu != b_coefs->size() ||
      natom * order_zu != c_coefs->size()) {
    rtErr("Output coefficient arrays (length as provided " + std::to_string(a_coefs->size()) +
          ", " + std::to_string(b_coefs->size()) + ", and " + std::to_string(c_coefs->size()) +
          ") must equal the order (" + std::to_string(order) + ") times the number of atoms (" +
          std::to_string(natom) + ").", "bSpline");
  }
  if (natom != a_init->size() || natom != b_init->size() || natom != c_init->size()) {
    rtErr("Output mesh alignment arrays (lengths as provided " + std::to_string(a_init->size()) +
          ", " + std::to_string(b_init->size()) + ", and " + std::to_string(c_init->size()) +
          ") must have lengths equal to the number of atoms (" + std::to_string(natom) + ").",
          "bSpline");
  }
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> void bSpline(const T x, const int order, T* coefs, T* dcoefs) {
  const T value_zero = 0.0;
  const T value_one = 1.0;
  if (order < 2) {
    rtErr("B-splines are undefined for order less than 2.", "bSpline");
  }
  else if (dcoefs != nullptr && order == 2) {
    rtErr("B-spline derivatives are undefined for orders less than 3.", "bSpline");
  }
  else if (x < value_zero || x > value_one) {
    rtErr("B-splines are computed for x in the range [0, 1].", "bSpline");
  }
  coefs[0] = value_one - x;
  coefs[1] = x;
  if (order == 2) {
    return;
  }

  // One pass to order three.  If the array of derivatives is defined, compute them immediately.
  if (order == 3 && dcoefs != nullptr) {
    dcoefs[0] = -coefs[0];
    dcoefs[1] = coefs[0] - coefs[1];
    dcoefs[2] = coefs[1];
  }
  const T value_half = 0.5;
  coefs[2] = value_half * x * coefs[1];
  coefs[0] = value_half * (value_one - x) * coefs[0];
  coefs[1] = value_one - coefs[2] - coefs[0];
  if (order == 3) {
    return;
  }

  // Advance to the fourth-order B-spline.  If derivatives are needed, compute them immediately.
  if (order == 4 && dcoefs != nullptr) {
    dcoefs[0] = -coefs[0];
    dcoefs[1] = coefs[0] - coefs[1];
    dcoefs[2] = coefs[1] - coefs[2];
    dcoefs[3] = coefs[2];
  }
  const T value_three = 3.0;
  const T value_third = value_one / value_three;
  coefs[3] = value_third * x * coefs[2];
  coefs[2] = value_third * (((x + value_one) * coefs[1]) + ((value_three - x) * coefs[2]));
  coefs[0] = value_third * (value_one - x) * coefs[0];
  coefs[1] = value_one - coefs[3] - coefs[2] - coefs[0];
  if (order == 4) {
    return;
  }

  // Advance to higher-order splines.
  for (int n = 4; n < order - 1; n++) {
    const T value_n = n;
    const T value_invn = value_one / value_n;
    coefs[n] = value_invn * x * coefs[n - 1];
    const int halfway = (n >> 1);
    T running_sum = coefs[n];
    for (int i = n - 1; i > halfway; i--) {
      coefs[i] = value_invn * (((x + static_cast<T>(n - i)) * coefs[i - 1]) +
                                ((static_cast<T>(i + 1) - x) * coefs[i]));
      running_sum += coefs[i];
    }
    for (int i = halfway - 1; i >= 1; i--) {
      coefs[i] = value_invn * (((x + static_cast<T>(n - i)) * coefs[i - 1]) +
                                ((static_cast<T>(i + 1) - x) * coefs[i]));
      running_sum += coefs[i];
    }
    coefs[0] = value_invn * (value_one - x) * coefs[0];
    coefs[halfway] = value_one - running_sum - coefs[0];
  }
  
  // Compute the derivatives, if required.
  const int om = order - 1;
  if (dcoefs != nullptr) {
    dcoefs[0] = -coefs[0];
    for (int i = 1; i < om; i++) {
      dcoefs[i] = coefs[i - 1] - coefs[i];
    }
    dcoefs[om] = coefs[order - 2];
  }
  
  // Carry out the final recursive pass to get the spline cofficients.
  const T value_f = om;
  const T value_invf = value_one / value_f;
  coefs[om] = value_invf * x * coefs[order - 2];
  const int halfway = (order >> 1);
  T running_fsum = coefs[om];
  for (int i = om - 1; i > halfway; i--) {
    coefs[i] = value_invf * (((x + static_cast<T>(om - i)) * coefs[i - 1]) +
                              ((static_cast<T>(i + 1) - x) * coefs[i]));
    running_fsum += coefs[i];
  }
  for (int i = halfway - 1; i >= 1; i--) {
    coefs[i] = value_invf * (((x + static_cast<T>(om - i)) * coefs[i - 1]) +
                              ((static_cast<T>(i + 1) - x) * coefs[i]));
    running_fsum += coefs[i];
  }
  coefs[0] = value_invf * (value_one - x) * coefs[0];
  coefs[halfway] = value_one - running_fsum - coefs[0];
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> bSpline(const T x, const int order) {
  std::vector<T> result(order), workspace(order - 1);
  bSpline(x, order, result.data());
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void bSpline(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, const int natom,
             const int order, const Tcalc* umat_cell, const Tcoord* invu_cell, const int mesh_na,
             const int mesh_nb, const int mesh_nc, const Tcoord* invu_mesh, Tcalc* a_coefs,
             Tcalc* b_coefs, Tcalc* c_coefs, int* a_init, int* b_init, int* c_init,
             Tcalc* da_coefs, Tcalc* db_coefs, Tcalc* dc_coefs, const Tcalc coordinate_scale) {
  const Tcalc* null_pointer = static_cast<Tcalc*>(nullptr);
  const Tcalc value_one = 1.0;
  const bool tcoord_is_integer = isSignedIntegralScalarType<Tcoord>();
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const Tcalc inv_crd_scl = value_one / coordinate_scale;
  const Tcalc tmesh_na = mesh_na;
  const Tcalc tmesh_nb = mesh_nb;
  const Tcalc tmesh_nc = mesh_nc;
  for (int i = 0; i < natom; i++) {

    // Compute the bin offsets for the atom.
    Tcalc frac_txi, frac_tyi, frac_tzi;
    if (tcoord_is_integer) {
      const Tcalc txi = static_cast<Tcalc>(xcrd[i]) * inv_crd_scl;
      const Tcalc tyi = static_cast<Tcalc>(ycrd[i]) * inv_crd_scl;
      const Tcalc tzi = static_cast<Tcalc>(zcrd[i]) * inv_crd_scl;
      frac_txi = (umat_cell[0] * txi) + (umat_cell[3] * tyi) + (umat_cell[6] * tzi);
      frac_tyi =                        (umat_cell[4] * tyi) + (umat_cell[7] * tzi);
      frac_tzi =                                               (umat_cell[8] * tzi);
      Tcalc shft_xi, shft_yi, shft_zi;
      if (tcalc_is_double) {
        shft_xi = floor(frac_txi);
        shft_yi = floor(frac_tyi);
        shft_zi = floor(frac_tzi);
      }
      else {
        shft_xi = floorf(frac_txi);
        shft_yi = floorf(frac_tyi);
        shft_zi = floorf(frac_tzi);
      }
      const Tcoord ishft_xi = -shft_xi;
      const Tcoord ishft_yi = -shft_yi;
      const Tcoord ishft_zi = -shft_zi;

      // Re-calculate the properly imaged coordinates making use of the fixed-precision unit cell
      // inverse transformation matrix.
      const Tcoord nxcrdi = xcrd[i] + (ishft_xi * invu_cell[0]) + (ishft_yi * invu_cell[3]) +
                            (ishft_zi * invu_cell[6]);
      const Tcoord nycrdi = ycrd[i] + (ishft_yi * invu_cell[4]) + (ishft_zi * invu_cell[7]);
      const Tcoord nzcrdi = zcrd[i] + (ishft_zi * invu_cell[8]);
      const Tcalc ntxi = static_cast<Tcalc>(nxcrdi) * inv_crd_scl;
      const Tcalc ntyi = static_cast<Tcalc>(nycrdi) * inv_crd_scl;
      const Tcalc ntzi = static_cast<Tcalc>(nzcrdi) * inv_crd_scl;

      // Re-calculate the fractional coordinates
      frac_txi = (umat_cell[0] * ntxi) + (umat_cell[3] * ntyi) + (umat_cell[6] * ntzi);
      frac_tyi =                         (umat_cell[4] * ntyi) + (umat_cell[7] * ntzi);
      frac_tzi =                                                 (umat_cell[8] * ntzi);
      frac_txi *= tmesh_na;
      frac_tyi *= tmesh_nb;
      frac_tzi *= tmesh_nc;
      a_init[i] = frac_txi;
      b_init[i] = frac_tyi;
      c_init[i] = frac_tzi;

      // Re-calculate the delta
      frac_txi = static_cast<Tcalc>(nxcrdi - (a_init[i] * invu_mesh[0]) -
                                             (b_init[i] * invu_mesh[3]) -
                                             (c_init[i] * invu_mesh[6])) * inv_crd_scl;
      frac_tyi = static_cast<Tcalc>(nycrdi - (b_init[i] * invu_mesh[4]) -
                                             (c_init[i] * invu_mesh[7])) * inv_crd_scl;
      frac_tzi = static_cast<Tcalc>(nzcrdi - (c_init[i] * invu_mesh[8])) * inv_crd_scl;
    }
    else {
      const Tcalc txi = xcrd[i];
      const Tcalc tyi = ycrd[i];
      const Tcalc tzi = zcrd[i];
      frac_txi = (umat_cell[0] * txi) + (umat_cell[3] * tyi) + (umat_cell[6] * tzi);
      frac_tyi =                        (umat_cell[4] * tyi) + (umat_cell[7] * tzi);
      frac_tzi =                                               (umat_cell[8] * tzi);
      if (tcalc_is_double) {
        frac_txi -= floor(frac_txi);
        frac_tyi -= floor(frac_tyi);
        frac_tzi -= floor(frac_tzi);
      }
      else {
        frac_txi -= floorf(frac_txi);
        frac_tyi -= floorf(frac_tyi);
        frac_tzi -= floorf(frac_tzi);
      }
      frac_txi *= tmesh_na;
      frac_tyi *= tmesh_nb;
      frac_tzi *= tmesh_nc;
      a_init[i] = frac_txi;
      b_init[i] = frac_tyi;
      c_init[i] = frac_tzi;
      frac_txi -= static_cast<Tcalc>(a_init[i]);
      frac_tyi -= static_cast<Tcalc>(b_init[i]);
      frac_tzi -= static_cast<Tcalc>(c_init[i]);
    }

    // Compute the B-spline coefficients along each unit cell axis.
    if (da_coefs != null_pointer && db_coefs != null_pointer && dc_coefs != null_pointer) {
      bSpline(frac_txi, order, &a_coefs[order * i], &da_coefs[order * i]);
      bSpline(frac_tyi, order, &b_coefs[order * i], &db_coefs[order * i]);
      bSpline(frac_tzi, order, &c_coefs[order * i], &dc_coefs[order * i]);
    }
    else {
      bSpline(frac_txi, order, &a_coefs[order * i]);
      bSpline(frac_tyi, order, &b_coefs[order * i]);
      bSpline(frac_tzi, order, &c_coefs[order * i]);
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void bSpline(const std::vector<Tcoord> &xcrd, const std::vector<Tcoord> &ycrd,
             const std::vector<Tcoord> &zcrd, const int order, const Tcalc* umat_cell,
             const Tcoord* invu_cell, const int mesh_na, const int mesh_nb, const int mesh_nc,
             const Tcoord* invu_mesh, std::vector<Tcalc> *a_coefs, std::vector<Tcalc> *b_coefs,
             std::vector<Tcalc> *c_coefs, std::vector<int> *a_init, std::vector<int> *b_init,
             std::vector<int> *c_init, std::vector<Tcalc> *da_coefs, std::vector<Tcalc> *db_coefs,
             std::vector<Tcalc> *dc_coefs, const Tcalc coordinate_scale) {
  bSplineInputChecks(xcrd, ycrd, zcrd, order, a_coefs, b_coefs, c_coefs, a_init, b_init, c_init);
  if (da_coefs != nullptr && db_coefs != nullptr && dc_coefs != nullptr) {
    bSpline(xcrd.data(), ycrd.data(), zcrd.data(), xcrd.size(), order, umat_cell, invu_cell,
            mesh_na,  mesh_nb, mesh_nc, invu_mesh, a_coefs->data(), b_coefs->data(),
            c_coefs->data(), a_init->data(), b_init->data(), c_init->data(), da_coefs->data(),
            db_coefs->data(), dc_coefs->data(), coordinate_scale);
  }
  else {
    bSpline(xcrd.data(), ycrd.data(), zcrd.data(), xcrd.size(), order, umat_cell, invu_cell,
            mesh_na,  mesh_nb, mesh_nc, invu_mesh, a_coefs->data(), b_coefs->data(),
            c_coefs->data(), a_init->data(), b_init->data(), c_init->data(),
            static_cast<Tcalc*>(nullptr), static_cast<Tcalc*>(nullptr),
            static_cast<Tcalc*>(nullptr), coordinate_scale);    
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void bSpline(const Hybrid<Tcoord> &xcrd, const Hybrid<Tcoord> &ycrd, const Hybrid<Tcoord> &zcrd,
             const int order, const Tcalc* umat_cell, const Tcoord* invu_cell, const int mesh_na,
             const int mesh_nb, const int mesh_nc, const Tcoord* invu_mesh, Hybrid<Tcalc> *a_coefs,
             Hybrid<Tcalc> *b_coefs, Hybrid<Tcalc> *c_coefs, Hybrid<int> *a_init,
             Hybrid<int> *b_init, Hybrid<int> *c_init, Hybrid<Tcalc> *da_coefs,
             Hybrid<Tcalc> *db_coefs, Hybrid<Tcalc> *dc_coefs, const Tcalc coordinate_scale) {
  bSplineInputChecks(xcrd, ycrd, zcrd, order, a_coefs, b_coefs, c_coefs, a_init, b_init, c_init);
  if (da_coefs != nullptr && db_coefs != nullptr && dc_coefs != nullptr) {
    bSpline(xcrd.data(), ycrd.data(), zcrd.data(), xcrd.size(), order, umat_cell, invu_cell,
            mesh_na,  mesh_nb, mesh_nc, invu_mesh, a_coefs->data(), b_coefs->data(),
            c_coefs->data(), a_init->data(), b_init->data(), c_init->data(), da_coefs->data(),
            db_coefs->data(), dc_coefs->data(), coordinate_scale);
  }
  else {
    bSpline(xcrd.data(), ycrd.data(), zcrd.data(), xcrd.size(), order, umat_cell, invu_cell,
            mesh_na,  mesh_nb, mesh_nc, invu_mesh, a_coefs->data(), b_coefs->data(),
            c_coefs->data(), a_init->data(), b_init->data(), c_init->data(),
            static_cast<Tcalc*>(nullptr), static_cast<Tcalc*>(nullptr),
            static_cast<Tcalc*>(nullptr), coordinate_scale);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> bSplineNoUnity(T x, int order) {

  // Begin by constructing the second-order spline.
  const T value_one = 1.0;
  std::vector<T> result(order);
  if (order < 2) {
    rtErr("B-splines are undefined for order less than 2.", "bSplineNoUnity");
  }
  result[0] = value_one - x;
  result[1] = x;
  if (order == 2) {
    return result;
  }

  // Make one pass to build the third-order spline.
  const T value_half = 0.5;
  const T value_two  = 2.0;
  result[2] = value_half * x * result[1];
  result[1] = value_half * (((x + value_one) * result[0]) + ((value_two - x) * result[1]));
  result[0] = value_half * (value_one - x) * result[0];
  if (order == 3) {
    return result;
  }

  // Make one pass to build the fourth-order spline.
  const T value_three = 3.0;
  const T value_third = value_one / value_three;
  result[3] = value_third * x * result[2];
  result[2] = value_third * (((x + value_one) * result[1]) + ((value_three - x) * result[2]));
  result[1] = value_third * (((x + value_two) * result[0]) + ((value_two - x) * result[1]));
  result[0] = value_third * (value_one - x) * result[0];
  if (order == 4) {
    return result;
  }

  // Make additional passes to build higher order splines
  for (int n = 4; n < order; n++) {
    const T value_n = n;
    const T value_invn = value_one / value_n;
    result[n] = value_invn * x * result[n - 1];
    for (int i = n - 1; i >= 1; i--) {
      result[i] = value_invn * (((x + static_cast<T>(n - i)) * result[i - 1]) +
                                ((static_cast<T>(i + 1) - x) * result[i]));
    }
    result[0] = value_invn * (value_one - x) * result[0];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> dBSpline(const T x, const int order,
                                              const bool exploit_partition) {
  std::vector<T> workspace;
  if (exploit_partition) {
    workspace.resize(order - 1);
    bSpline(x, order - 1, workspace.data());
  }
  else {
    workspace = bSplineNoUnity(x, order - 1);
    
  }
  std::vector<T> result(order);
  result[0] = -workspace[0];
  for (int i = 1; i < order - 1; i++) {
    result[i] = workspace[i - 1] - workspace[i];
  }
  result[order - 1] = workspace[order - 2];
  return result;
}


} // namespace stmath
} // namespace stormm
