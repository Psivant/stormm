// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename T>
void PhaseSpace::fill(const T* xcrd, const T* ycrd, const T* zcrd, const TrajectoryKind kind,
                      const CoordinateCycle cycle_in, const int scale_bits,
                      const double* box_dims) {
  double *xptr, *yptr, *zptr;
  switch (cycle_in) {
  case CoordinateCycle::BLACK:
    switch (kind) {
    case TrajectoryKind::POSITIONS:
      xptr = x_alt_coordinates.data();
      yptr = y_alt_coordinates.data();
      zptr = z_alt_coordinates.data();
      break;
    case TrajectoryKind::VELOCITIES:
      xptr = x_alt_velocities.data();
      yptr = y_alt_velocities.data();
      zptr = z_alt_velocities.data();
      break;
    case TrajectoryKind::FORCES:
      xptr = x_alt_forces.data();
      yptr = y_alt_forces.data();
      zptr = z_alt_forces.data();
      break;
    }
    break;
  case CoordinateCycle::WHITE:
    switch (kind) {
    case TrajectoryKind::POSITIONS:

      // Only in the case of WHITE POSITIONS should the box dimensions be filled.
      if (box_dims != nullptr) {
        double* boxptr = box_dimensions.data();
        for (int i = 0; i < 6; i++) {
          boxptr[i] = box_dims[i];
        }
        computeBoxTransform(boxptr, box_space_transform.data(), inverse_transform.data());
        unit_cell = determineUnitCellTypeByShape(inverse_transform.data());
      }
      xptr = x_coordinates.data();
      yptr = y_coordinates.data();
      zptr = z_coordinates.data();
      break;
    case TrajectoryKind::VELOCITIES:
      xptr = x_velocities.data();
      yptr = y_velocities.data();
      zptr = z_velocities.data();
      break;
    case TrajectoryKind::FORCES:
      xptr = x_forces.data();
      yptr = y_forces.data();
      zptr = z_forces.data();
      break;
    }
    break;
  }
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
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void PhaseSpace::fill(const std::vector<T> &xcrd, const std::vector<T> &ycrd,
                      const std::vector<T> &zcrd, const TrajectoryKind kind,
                      const CoordinateCycle cycle_in, const int scale_bits,
                      const std::vector<double> &box_dims) {
  fill(xcrd.data(), ycrd.data(), zcrd.data(), kind, cycle_in, scale_bits, box_dims.data());
}

} // namespace trajectory
} // namespace stormm
