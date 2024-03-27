// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateFrame::fill(const T* xcrd, const T* ycrd, const T* zcrd, const int scale_bits,
                           const double* box_dims) {
  double* xptr = x_coordinates.data();
  double* yptr = y_coordinates.data();
  double* zptr = z_coordinates.data();
  if (scale_bits == 0) {
    for (int i = 0; i < atom_count; i++) {
      xptr[i] = xcrd[i];
      yptr[i] = ycrd[i];
      zptr[i] = zcrd[i];
    }
  }
  else {
    const double conv_factor = pow(2.0, -scale_bits);
    for (int i = 0; i < atom_count; i++) {
      xptr[i] = static_cast<double>(xcrd[i]) * conv_factor;
      yptr[i] = static_cast<double>(ycrd[i]) * conv_factor;
      zptr[i] = static_cast<double>(zcrd[i]) * conv_factor;
    }
  }
  if (box_dims != nullptr) {
    double* boxptr = box_dimensions.data();
    for (int i = 0; i < 6; i++) {
      boxptr[i] = box_dims[i];
    }
    computeBoxTransform(boxptr, box_space_transform.data(), inverse_transform.data());
    unit_cell = determineUnitCellTypeByShape(inverse_transform.data());
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateFrame::fill(const std::vector<T> &xcrd, const std::vector<T> &ycrd,
                           const std::vector<T> &zcrd, const int scale_bits,
                           const std::vector<double> &box_dims) {
  fill(xcrd.data(), ycrd.data(), zcrd.data(), scale_bits, box_dims.data());
}

} // namespace trajectory
} // namespace stormm
