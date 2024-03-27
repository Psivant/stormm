#include "copyright.h"
#include "matrix_ops.h"
#include "tricubic_cell.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
TricubicStencil::TricubicStencil(const Interpolant kind_in) :
    kind{kind_in},
    transform{4096, "tricubic_transform"}
{
  std::vector<double> tmp(4096);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4 ; k++) {
        const int col_idx = (4 * ((4 * k) + j)) + i;
        for (int xv = 0; xv < 2; xv++) {
          const int xfac = (i == 0) ? 1 : xv;
          const int dxfac = (i == 0) ? 0 : (i == 1) ? 1 : (i == 2) ? 2 * xv : 3 * xv * xv;
          for (int yv = 0; yv < 2; yv++) {
            const int yfac = (j == 0) ? 1 : yv;
            const int dyfac = (j == 0) ? 0 : (j == 1) ? 1 : (j == 2) ? 2 * yv : 3 * yv * yv;
            for (int zv = 0; zv < 2; zv++) {
              const int row_idx = (2 * ((2 * zv) + yv)) + xv;
              const int zfac = (k == 0) ? 1 : zv;
              const int dzfac = (k == 0) ? 0 : (k == 1) ? 1 : (k == 2) ? 2 * zv : 3 * zv * zv;

              // Compute the polynomial coefficient for the function value
              tmp[(col_idx * 64) + row_idx] = xfac * yfac * zfac;

              // Compute the polynomial coefficients for function first derivatives
              tmp[(col_idx * 64) + row_idx +  8] = dxfac * yfac * zfac;
              tmp[(col_idx * 64) + row_idx + 16] = xfac * dyfac * zfac;
              tmp[(col_idx * 64) + row_idx + 24] = xfac * yfac * dzfac;

              switch (kind_in) {
              case Interpolant::SMOOTHNESS:

                // Compute the polynomial coefficients for function cross derivatives
                tmp[(col_idx * 64) + row_idx + 32] = dxfac * dyfac * zfac;
                tmp[(col_idx * 64) + row_idx + 40] = dxfac * yfac * dzfac;
                tmp[(col_idx * 64) + row_idx + 48] = xfac * dyfac * dzfac;

                // Compute the polynomial coefficient for the function triple derivative
                tmp[(col_idx * 64) + row_idx + 56] = dxfac * dyfac * dzfac;
                break;
              case Interpolant::FUNCTION_VALUE:
                break;
              }
            }
          }
        }

        // Compute polynomial coefficients for additional points in the element interior.  For a
        // tricubic cell, this selects 32 points arranged in six sets of four around a central set
        // of eight, like even cubes connected in a three-dimensional jack or star configuration.
        // The two sets of four points lying in planes parallel to xy, xz, and yz take the places
        // of mixed partial derivatives along these directions, while the central eight points take
        // the place of the triple derivative evaluations in the "SMOOTHNESS" case.
        switch (kind_in) {
        case Interpolant::SMOOTHNESS:
          break;
        case Interpolant::FUNCTION_VALUE:
          for (int xv = 0; xv < 2; xv++) {
            const double dbl_xv = xv;
            const double xyface_x = 0.375 + (0.25 * dbl_xv);
            const double xzface_x = 0.375 + (0.25 * dbl_xv);
            const double yzface_x = 0.125 + (0.75 * dbl_xv);
            const double center_x = 0.375 + (0.25 * dbl_xv);
            double xy_xfac = 1.0;
            double xz_xfac = 1.0;
            double yz_xfac = 1.0;
            double cn_xfac = 1.0;
            for (int ic = 0; ic < i; ic++) {
              xy_xfac *= xyface_x;
              xz_xfac *= xzface_x;
              yz_xfac *= yzface_x;
              cn_xfac *= center_x;
            }
            for (int yv = 0; yv < 2; yv++) {
              const double dbl_yv = yv;
              const double xyface_y = 0.375 + (0.25 * dbl_yv);
              const double xzface_y = 0.125 + (0.75 * dbl_yv);
              const double yzface_y = 0.375 + (0.25 * dbl_yv);
              const double center_y = 0.375 + (0.25 * dbl_yv);
              double xy_yfac = 1.0;
              double xz_yfac = 1.0;
              double yz_yfac = 1.0;
              double cn_yfac = 1.0;
              for (int jc = 0; jc < j; jc++) {
                xy_yfac *= xyface_y;
                xz_yfac *= xzface_y;
                yz_yfac *= yzface_y;
                cn_yfac *= center_y;
              }
              for (int zv = 0; zv < 2; zv++) {
                const double dbl_zv = zv;
                const double xyface_z = 0.125 + (0.75 * dbl_zv);
                const double xzface_z = 0.375 + (0.25 * dbl_zv);
                const double yzface_z = 0.375 + (0.25 * dbl_zv);
                const double center_z = 0.375 + (0.25 * dbl_zv);
                const int row_idx = (2 * ((2 * zv) + yv)) + xv;
                double xy_zfac = 1.0;
                double xz_zfac = 1.0;
                double yz_zfac = 1.0;
                double cn_zfac = 1.0;
                for (int kc = 0; kc < k; kc++) {
                  xy_zfac *= xyface_z;
                  xz_zfac *= xzface_z;
                  yz_zfac *= yzface_z;
                  cn_zfac *= center_z;
                }
                tmp[(col_idx * 64) + row_idx + 32] = xy_xfac * xy_yfac * xy_zfac;
                tmp[(col_idx * 64) + row_idx + 40] = xz_xfac * xz_yfac * xz_zfac;
                tmp[(col_idx * 64) + row_idx + 48] = yz_xfac * yz_yfac * yz_zfac;
                tmp[(col_idx * 64) + row_idx + 56] = cn_xfac * cn_yfac * cn_zfac;
              }
            }
          }
          break;
        }
      }
    }
  }
  
  // Invert the matrix to obtain the tricubic spline coefficients.  This is best done on the host,
  // as the matrix is too small to effectively utilize a GPU.
  std::vector<double> itmp(4096);
  const double* tmp_ptr = tmp.data();
  invertSquareMatrix(tmp_ptr, itmp.data(), 64);
  transform.putHost(itmp);
}

//-------------------------------------------------------------------------------------------------
Interpolant TricubicStencil::getKind() const {
  return kind;
}

//-------------------------------------------------------------------------------------------------
const double* TricubicStencil::data(const HybridTargetLevel tier) const {
  return transform.data(tier);
}

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
//-------------------------------------------------------------------------------------------------
const double* TricubicStencil::getDeviceViewToHostData() const {
  return transform.getDeviceValidHostPointer();
}
#  endif

//-------------------------------------------------------------------------------------------------
void TricubicStencil::upload() {
  transform.upload();
}

//-------------------------------------------------------------------------------------------------
void TricubicStencil::download() {
  transform.download();
}
#endif

//-------------------------------------------------------------------------------------------------
std::vector<double> TricubicStencil::exportMatrix() const {
  return transform.readHost();
}

//-------------------------------------------------------------------------------------------------
double3 stencilInternalOffset(const UnitCellAxis face_normal, const int a_index, const int b_index,
                              const int c_index, const double* invu) {
  double a_offset = 0.375;
  double b_offset = 0.375;
  double c_offset = 0.375;
  double a_length = 0.25;
  double b_length = 0.25;
  double c_length = 0.25;
  switch (face_normal) {
  case UnitCellAxis::A:
    a_offset -= 0.25;
    a_length = 0.75;
    break;
  case UnitCellAxis::B:
    b_offset -= 0.25;
    b_length = 0.75;
    break;
  case UnitCellAxis::C: 
    c_offset -= 0.25;
    c_length = 0.75;
    break;
  }
  a_offset += static_cast<double>(a_index) * a_length;
  b_offset += static_cast<double>(b_index) * b_length;
  c_offset += static_cast<double>(c_index) * c_length;
  return { (invu[0] * a_offset) + (invu[3] * b_offset) + (invu[6] * c_offset),
                                  (invu[4] * b_offset) + (invu[7] * c_offset),
                                                         (invu[8] * c_offset) };
}

//-------------------------------------------------------------------------------------------------
double3 stencilInternalOffset(const int a_index, const int b_index, const int c_index,
                              const double* invu) {
  double a_offset = 0.375;
  double b_offset = 0.375;
  double c_offset = 0.375;
  double a_length = 0.25;
  double b_length = 0.25;
  double c_length = 0.25;
  a_offset += static_cast<double>(a_index) * a_length;
  b_offset += static_cast<double>(b_index) * b_length;
  c_offset += static_cast<double>(c_index) * c_length;
  return { (invu[0] * a_offset) + (invu[3] * b_offset) + (invu[6] * c_offset),
                                  (invu[4] * b_offset) + (invu[7] * c_offset),
                                                         (invu[8] * c_offset) };
}

//-------------------------------------------------------------------------------------------------
void incorporateStencilOrigin(const int95_t orig_x, const int95_t orig_y, const int95_t orig_z,
                              const double3 pt_xyz, const double scale_factor, int95_t *point_x,
                              int95_t *point_y, int95_t *point_z) {
  const int95_t x_i = hostDoubleToInt95(pt_xyz.x * scale_factor);
  const int95_t y_i = hostDoubleToInt95(pt_xyz.y * scale_factor);
  const int95_t z_i = hostDoubleToInt95(pt_xyz.z * scale_factor);
  *point_x = hostSplitFPSum(orig_x, x_i);
  *point_y = hostSplitFPSum(orig_y, y_i);
  *point_z = hostSplitFPSum(orig_z, z_i);
}

} // namespace stmath
} // namespace stormm
