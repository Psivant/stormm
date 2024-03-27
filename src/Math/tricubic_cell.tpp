// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename T>
TricubicCell<T>::TricubicCell() :
    coefficients{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    f{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dx{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    dy{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    dxx{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dxy{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    dxz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dyy{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    dyz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dxxx{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    dxxy{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dxxz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    dxyy{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, dxyz{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    origin_x{0.0}, origin_y{0.0}, origin_z{0.0},
    umat{ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 },
    invu{ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 }
{}
    
//-------------------------------------------------------------------------------------------------
template <typename T>
TricubicCell<T>::TricubicCell(const TricubicStencil weights_matrix,
                              const std::vector<double> &bounds, const std::vector<T> &f_in,
                              const std::vector<T> &dx_in, const std::vector<T> &dy_in,
                              const std::vector<T> &dz_in, const std::vector<T> &dxx_in,
                              const std::vector<T> &dxy_in, const std::vector<T> &dxz_in,
                              const std::vector<T> &dyy_in, const std::vector<T> &dyz_in,
                              const std::vector<T> &dxxx_in, const std::vector<T> &dxxy_in,
                              const std::vector<T> &dxxz_in, const std::vector<T> &dxyy_in,
                              const std::vector<T> &dxyz_in) :
    coefficients{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    f{ f_in[0], f_in[1], f_in[2], f_in[3], f_in[4], f_in[5], f_in[6], f_in[7] },
    dx{ dx_in[0], dx_in[1], dx_in[2], dx_in[3], dx_in[4], dx_in[5], dx_in[6], dx_in[7] },
    dy{ dy_in[0], dy_in[1], dy_in[2], dy_in[3], dy_in[4], dy_in[5], dy_in[6], dy_in[7] },
    dz{ dz_in[0], dz_in[1], dz_in[2], dz_in[3], dz_in[4], dz_in[5], dz_in[6], dz_in[7] },
    dxx{ dxx_in[0], dxx_in[1], dxx_in[2], dxx_in[3], dxx_in[4], dxx_in[5], dxx_in[6], dxx_in[7] },
    dyy{ dyy_in[0], dyy_in[1], dyy_in[2], dyy_in[3], dyy_in[4], dyy_in[5], dyy_in[6], dyy_in[7] },
    dxy{ dxy_in[0], dxy_in[1], dxy_in[2], dxy_in[3], dxy_in[4], dxy_in[5], dxy_in[6], dxy_in[7] },
    dxz{ dxz_in[0], dxz_in[1], dxz_in[2], dxz_in[3], dxz_in[4], dxz_in[5], dxz_in[6], dxz_in[7] },
    dyz{ dyz_in[0], dyz_in[1], dyz_in[2], dyz_in[3], dyz_in[4], dyz_in[5], dyz_in[6], dyz_in[7] },
    dxxx{ dxxx_in[0], dxxx_in[1], dxxx_in[2], dxxx_in[3], dxxx_in[4], dxxx_in[5], dxxx_in[6],
          dxxx_in[7] },
    dxxy{ dxxy_in[0], dxxy_in[1], dxxy_in[2], dxxy_in[3], dxxy_in[4], dxxy_in[5], dxxy_in[6],
          dxxy_in[7] },
    dxxz{ dxxz_in[0], dxxz_in[1], dxxz_in[2], dxxz_in[3], dxxz_in[4], dxxz_in[5], dxxz_in[6],
          dxxz_in[7] },
    dxyy{ dxyy_in[0], dxyy_in[1], dxyy_in[2], dxyy_in[3], dxyy_in[4], dxyy_in[5], dxyy_in[6],
          dxyy_in[7] },
    dxyz{ dxyz_in[0], dxyz_in[1], dxyz_in[2], dxyz_in[3], dxyz_in[4], dxyz_in[5], dxyz_in[6],
          dxyz_in[7] }, origin_x{0.0}, origin_y{0.0}, origin_z{0.0},
    umat{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    invu{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
{
  // Fill out the grid cell's dimensions
  if (bounds.size() > 3) {
    origin_x = bounds[0];
    origin_y = bounds[1];
    origin_z = bounds[2];
  }
  if (bounds.size() == 4LLU) {
    invu[0] = bounds[3];
    invu[4] = bounds[3];
    invu[8] = bounds[3];
  }
  else if (bounds.size() == 6LLU) {
    invu[0] = bounds[3];
    invu[4] = bounds[4];
    invu[8] = bounds[5];
  }
  else if (bounds.size() == 12LLU) {
    for (int i = 0; i < 9; i++) {
      invu[i] = bounds[i + 3];
    }
  }
  else {
    rtErr("Invalid dimesion " + std::to_string(bounds.size()) + " on the bounds array.  There "
          "must be four, six, or twelve entries (Cartesian origin X, Y, and Z plus 1-3 edge "
          "lengths or all entries of a complete 3 x 3 inverse unit cell transformation matrix.",
          "TricubicCell");
  }

  // Check the cell lengths
  if (leibnizDeterminant(invu, 3) < 1.0e-8) {
    rtErr("An invalid cell was defined with vectors [ " +
          realToString(invu[0], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[1], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[2], 7, 4, NumberFormat::STANDARD_REAL) + "], [" +
          realToString(invu[3], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[4], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[5], 7, 4, NumberFormat::STANDARD_REAL) + "], and [" +
          realToString(invu[6], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[7], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(invu[8], 7, 4, NumberFormat::STANDARD_REAL) + "].", "TricubicCell");
  }
  else {
    const std::vector<T> invu_copy = { invu[0], invu[1], invu[2], invu[3], invu[4], invu[5],
                                       invu[6], invu[7], invu[8] };
    invertSquareMatrix(invu_copy.data(), umat, 3);
  }
  
  // Fill out the solution vector to the system of 64 equations
  std::vector<double> b(64);
  for (int i = 0; i < 8; i++) {
    b[i     ] = f[i];
    if (bounds.size() <= 6) {
      b[i +  8] = dx[i] * invu[0];
      b[i + 16] = dy[i] * invu[4];
      b[i + 24] = dz[i] * invu[8];
      switch (weights_matrix.getKind()) {
      case Interpolant::SMOOTHNESS:
        b[i + 32] = dxy[i] * invu[0] * invu[4];
        b[i + 40] = dxz[i] * invu[0] * invu[8];
        b[i + 48] = dyz[i] * invu[4] * invu[8];
        b[i + 56] = dxyz[i] * invu[0] * invu[4] * invu[8];
        break;
      case Interpolant::FUNCTION_VALUE:
        b[i + 32] = dxy[i];
        b[i + 40] = dxz[i];
        b[i + 48] = dyz[i];
        b[i + 56] = dxyz[i];
        break;
      }
    }
    else {

      // The general, triclinic case.  Use the chain rule to apply the appropriate transformation
      // matrix.
      //
      // du/da = (du/dx)(dx/da) + (du/dy)(dy/da) + (du/dz)(dz/da)
      // x = (invu[0] * a) + (invu[3] * b) + (invu[6] * c)
      // y =                 (invu[4] * b) + (invu[7] * c)
      // z =                                 (invu[8] * c)
      // dx/da = invu[0], dy/da = 0, and dz/da = 0
      b[i +  8] = (dx[i] * invu[0]);

      // du/db = (du/dx)(dx/db) + (du/dy)(dy/db) + (du/dz)(dz/db)
      // dx/db = inuv[3], dy/db = invu[4], and dz/db = 0
      b[i + 16] = (dx[i] * invu[3]) + (dy[i] * invu[4]);

      // du/dc = (du/dx)(dx/dc) + (du/dy)(dy/dc) + (du/dz)(dz/dc)
      // dx/dc = invu[6], dy/dc = invu[7], and dz/dc = invu[8]      
      b[i + 24] = (dx[i] * invu[6]) + (dy[i] * invu[7]) + (dz[i] * invu[8]);
      switch (weights_matrix.getKind()) {
      case Interpolant::SMOOTHNESS:
        {
          // d2u/dab = (d/db)(du/da) = (d/db) [ invu[0] * (du/dx) ] = invu[0] * (d/db)(du/dx)
          //   = invu[0] * ((d2u/dx2)(dx/db) + (d2u/dxy)(dy/db) + (d2u/dxz)(dz/db)
          //   = invi[0] * ((invu[3] * (d2u/dx2)) + (invu[4] * (d2u/dxy)) 
          b[i + 32] = invu[0] * ((invu[3] * dxx[i]) + (invu[4] * dxy[i]));

          // d2u/dac = (d/dc)(du/da) = (d/dc) [ invu[0] * (du/dx) ] = invu[0] * (d/dc)(du/dx)
          //   = invu[0] * ((d2u/dx2)(dx/dc) + (d2u/dxy)(dy/dc) + (d2u/dxz)(dz/dc))
          //   = invu[0] * ((invu[6] * (d2u/dx2)) + (invu[7] * (d2u/dxy)) + (invu[8] * (d2u/dxz)))
          const T tcol_dx_sum = (invu[6] * dxx[i]) + (invu[7] * dxy[i]) + (invu[8] * dxz[i]);
          b[i + 40] = invu[0] * tcol_dx_sum;

          // d2u/dbc = (d/dc)(du/db) = (d/dc) [ (invu[3] * (du/dx)) + (invu[4] * (du/dy)) ]
          //   = (invu[3] * (d/dc) [ (du/dx) ]) + (invu[4] * (d/dc) [ (du/dy) ])
          //   = (invu[3] * ((invu[6] * (d2u/dx2)) + (invu[7] * (d2u/dxy)) +
          //                 (invu[8] * (d2u/dxz)))) +
          //     (invu[4] * ((d2u/dxy)(dx/dc) + (d2u/dyy)(dy/dc) + (d2u/dyz)(dz/dc))
          //   = (invu[3] * ((invu[6] * (d2u/dx2)) + (invu[7] * (d2u/dxy)) +
          //                 (invu[8] * (d2u/dxz)))) +
          //     (invu[4] * ((invu[6] * (d2u/dxy)) + (invu[7] * (d2u/dyy)) +
          //                 (invu[8] * (d2u/dyz))))
          b[i + 48] = (invu[3] * tcol_dx_sum) + (invu[4] * ((invu[6] * dxy[i]) +
                                                            (invu[7] * dyy[i]) +
                                                            (invu[8] * dyz[i])));

          // d3u/dabc = (d/da)(d2u/dbc)
          b[i + 56] = invu[0] * ((invu[3] * ((invu[6] * dxxx[i]) + (invu[7] * dxxy[i]) +
                                             (invu[8] * dxxz[i]))) +
                                 (invu[4] * ((invu[6] * dxxy[i]) + (invu[7] * dxyy[i]) +
                                             (invu[8] * dxyz[i]))));
        }
        break;
      case Interpolant::FUNCTION_VALUE:
        b[i + 32] = dxy[i];
        b[i + 40] = dxz[i];
        b[i + 48] = dyz[i];
        b[i + 56] = dxyz[i];
        break;
      }
    }
  }
  std::vector<double> dcoeffs(64);
  matrixVectorMultiply(weights_matrix.data(), b.data(), dcoeffs.data(), 64, 64);
  for (int i = 0; i < 64; i++) {
    coefficients[i] = dcoeffs[i];
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
TricubicCell<T>::TricubicCell(const TricubicStencil weights_matrix,
                              const std::vector<double> &bounds, const std::vector<T> &f_in,
                              const std::vector<T> &dx_in, const std::vector<T> &dy_in,
                              const std::vector<T> &dz_in, const std::vector<T> &dxy_in,
                              const std::vector<T> &dxz_in, const std::vector<T> &dyz_in,
                              const std::vector<T> &dxyz_in) :
  TricubicCell(weights_matrix, bounds, f_in, dx_in, dy_in, dz_in, std::vector<T>(8, 0.0),
               dxy_in, dxz_in, std::vector<T>(8, 0.0), dyz_in, std::vector<T>(8, 0.0),
               std::vector<T>(8, 0.0), std::vector<T>(8, 0.0), std::vector<T>(8, 0.0), dxyz_in)
{
  switch (weights_matrix.getKind()) {
  case Interpolant::SMOOTHNESS:

    // Check that the bounds array concerns an orthorhombic unit cell
    if (std::abs(invu[3]) > constants::tiny || std::abs(invu[6]) > constants::tiny ||
        std::abs(invu[7]) > constants::tiny) {
      rtErr("Partial derivatives d2/dx2, d2/dy2, d3/dx3, d3/dx2y, d3/dx2z, and d3/dxy2 cannot be "
            "assumed to be zero when making a non-orthorhombic mesh.", "TricubicCell");
    }
    break;
  case Interpolant::FUNCTION_VALUE:

    // Additional partial derivatives are unnecessary with a stencil involving only first
    // derivatives.
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
T TricubicCell<T>::getCoefficient(const int i, const int j, const int k) const {
  if (i > 3 || j > 3 || k > 3 || i < 0 || j < 0 || k < 0) {
    rtErr("A coefficient for x, y, and z powers " + std::to_string(i) + ", " + std::to_string(j) +
          ", and " + std::to_string(k) + " is not acceptable.", "TricubicCell", "getCoefficient");
  }
  return coefficients[(4 * ((4 * k) + j)) + i];
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> TricubicCell<T>::getCoefficients() const {
  std::vector<T> result(64);
  for (size_t i = 0; i < 64LLU; i++) {
    result[i] = coefficients[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TricubicCell<T>::setCoefficient(const T value, const int i, const int j, const int k) {
  if (i > 3 || j > 3 || k > 3 || i < 0 || j < 0 || k < 0) {
    rtErr("A coefficient for x, y, and z powers " + std::to_string(i) + ", " + std::to_string(j) +
          ", and " + std::to_string(k) + " is not acceptable.", "TricubicCell", "setCoefficient");
  }
  coefficients[(4 * ((4 * k) + j)) + i] = value;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
T TricubicCell<T>::getData(const FunctionLevel kind, const int i, const int j, const int k) const {
  switch (kind) {
  case FunctionLevel::VALUE:
    return f[(2 * ((2 * i) + j)) + k];
  case FunctionLevel::DX:
    return dx[(2 * ((2 * i) + j)) + k];
  case FunctionLevel::DY:
    return dy[(2 * ((2 * i) + j)) + k];
  case FunctionLevel::DZ:
    return dz[(2 * ((2 * i) + j)) + k];
  case FunctionLevel::DXY:
    return dxy[(2 * ((2 * i) + j)) + k];
  case FunctionLevel::DXZ:
    return dxz[(2 * ((2 * i) + j)) + k];
  case FunctionLevel::DYZ:
    return dyz[(2 * ((2 * i) + j)) + k];
  case FunctionLevel::DXYZ:
    return dxyz[(2 * ((2 * i) + j)) + k];
  case FunctionLevel::DXX:
  case FunctionLevel::DYY:
  case FunctionLevel::DZZ:
  case FunctionLevel::DXXX:
  case FunctionLevel::DXXY:
  case FunctionLevel::DXXZ:
  case FunctionLevel::DXYY:
  case FunctionLevel::DXZZ:
  case FunctionLevel::DYYY:
  case FunctionLevel::DYYZ:
  case FunctionLevel::DYZZ:
  case FunctionLevel::DZZZ:
    rtErr("A tricubic mesh does not store data for a function's " + getEnumerationName(kind) + ".",
          "TricubicCell", "getData");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> void TricubicCell<T>::setData(const T value, const FunctionLevel kind,
                                                    const int i, const int j, const int k) {
  switch (kind) {
  case FunctionLevel::VALUE:
    f[(2 * ((2 * i) + j)) + k] = value;
    break;
  case FunctionLevel::DX:
    dx[(2 * ((2 * i) + j)) + k] = value;
    break;
  case FunctionLevel::DY:
    dy[(2 * ((2 * i) + j)) + k] = value;
    break;
  case FunctionLevel::DZ:
    dz[(2 * ((2 * i) + j)) + k] = value;
    break;
  case FunctionLevel::DXY:
    dxy[(2 * ((2 * i) + j)) + k] = value;
    break;
  case FunctionLevel::DXZ:
    dxz[(2 * ((2 * i) + j)) + k] = value;
    break;
  case FunctionLevel::DYZ:
    dyz[(2 * ((2 * i) + j)) + k] = value;
    break;
  case FunctionLevel::DXYZ:
    dxyz[(2 * ((2 * i) + j)) + k] = value;
    break;
  case FunctionLevel::DXX:
  case FunctionLevel::DYY:
  case FunctionLevel::DZZ:
  case FunctionLevel::DXXX:
  case FunctionLevel::DXXY:
  case FunctionLevel::DXXZ:
  case FunctionLevel::DXYY:
  case FunctionLevel::DXZZ:
  case FunctionLevel::DYYY:
  case FunctionLevel::DYYZ:
  case FunctionLevel::DYZZ:
  case FunctionLevel::DZZZ:
    rtErr("A tricubic mesh does not store data for a function's " + getEnumerationName(kind) + ".",
          "TricubicCell", "setData");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T TricubicCell<T>::getCellOrigin(CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return origin_x;
  case CartesianDimension::Y:
    return origin_y;
  case CartesianDimension::Z:
    return origin_z;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T TricubicCell<T>::getCellLength(CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return invu[0];
  case CartesianDimension::Y:
    return invu[4];
  case CartesianDimension::Z:
    return invu[8];
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TricubicCell<T>::fractionalPosition(const T x, const T y, const T z, T *a_frac, T *b_frac,
                                         T *c_frac, const ExceptionResponse policy,
                                         const char* caller) const {
  const T xrel = x - origin_x;
  const T yrel = y - origin_y;
  const T zrel = z - origin_z;
  *a_frac = (umat[0] * xrel) + (umat[3] * yrel) + (umat[6] * zrel);
  *b_frac =                    (umat[4] * yrel) + (umat[7] * zrel);
  *c_frac =                                       (umat[8] * zrel);
  const T value_zero = 0.0;
  const T value_one  = 1.0;
  if (*a_frac < value_zero || *a_frac >= value_one || *b_frac < value_zero ||
      *b_frac >= value_one || *c_frac < value_zero || *c_frac >= value_one) {
    const std::string msg("Point (" + realToString(x, 7, 4, NumberFormat::STANDARD_REAL) +
                          ", " + realToString(y, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
                          realToString(z, 7, 4, NumberFormat::STANDARD_REAL) + ") is out of "
                          "bounds for a grid cell with origin (" +
                          realToString(origin_x, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
                          realToString(origin_y, 7, 4, NumberFormat::STANDARD_REAL) + ", " +
                          realToString(origin_z, 7, 4, NumberFormat::STANDARD_REAL) +
                          ") and dimensions [ " +
                          realToString(invu[0], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
                          realToString(invu[1], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
                          realToString(invu[2], 7, 4, NumberFormat::STANDARD_REAL) + " ] x [ " +
                          realToString(invu[3], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
                          realToString(invu[4], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
                          realToString(invu[5], 7, 4, NumberFormat::STANDARD_REAL) + " ] x [ " +
                          realToString(invu[6], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
                          realToString(invu[7], 7, 4, NumberFormat::STANDARD_REAL) + ", " +
                          realToString(invu[8], 7, 4, NumberFormat::STANDARD_REAL) + " ].");
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(msg, "TricubicCell", caller);
    case ExceptionResponse::WARN:
      rtWarn(msg, "TricubicCell", caller);
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> T TricubicCell<T>::evaluate(const T x, const T y, const T z) const {
  T a_frac, b_frac, c_frac;
  fractionalPosition(x, y, z, &a_frac, &b_frac, &c_frac, ExceptionResponse::DIE, "evaluate");
  const T value_zero = 0.0;
  const T value_one  = 1.0;
  T result = value_zero;
  T av = value_one;
  for (int i = 0; i < 4; i++) {
    T bv = value_one;
    for (int j = 0; j < 4; j++) {
      T cv = value_one;
      for (int k = 0; k < 4; k++) {
        result += coefficients[(4 * ((4 * k) + j)) + i] * av * bv * cv;
        cv *= c_frac;
      }
      bv *= b_frac;
    }
    av *= a_frac;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
template <typename T3> T3 TricubicCell<T>::derivative(const T x, const T y, const T z) const {
  T a_frac, b_frac, c_frac;
  fractionalPosition(x, y, z, &a_frac, &b_frac, &c_frac, ExceptionResponse::DIE, "evaluate");
  const T value_zero = 0.0;
  const T value_one  = 1.0;
  T result_da = value_zero;
  T result_db = value_zero;
  T result_dc = value_zero;

  // Compute the derivative along the one appropriate mesh vector and scale as necessary.
  // In orthorhombic meshes, A is parallel to X, B parallel to Y, and C parallel to Z.
  T av = value_one;
  T difac = value_one;
  for (int i = 1; i < 4; i++) {
    T bv = value_one;
    for (int j = 0; j < 4; j++) {
      T cv = value_one;
      for (int k = 0; k < 4; k++) {
        result_da += coefficients[(4 * ((4 * k) + j)) + i] * difac * av * bv * cv;
        cv *= c_frac;
      }
      bv *= b_frac;
    }
    av *= a_frac;
    difac += value_one;
  }
  av = value_one;
  for (int i = 0; i < 4; i++) {
    T bv = value_one;
    T djfac = value_one;
    for (int j = 1; j < 4; j++) {
      T cv = value_one;
      for (int k = 0; k < 4; k++) {
        result_db += coefficients[(4 * ((4 * k) + j)) + i] * av * djfac * bv * cv;
        cv *= c_frac;
      }
      bv *= b_frac;
      djfac += value_one;
    }
    av *= a_frac;
  }
  av = value_one;
  for (int i = 0; i < 4; i++) {
    T bv = value_one;
    for (int j = 0; j < 4; j++) {
      T cv = value_one;
      T dkfac = 1.0;
      for (int k = 1; k < 4; k++) {
        result_dc += coefficients[(4 * ((4 * k) + j)) + i] * av * bv * dkfac * cv;
        cv *= c_frac;
        dkfac += value_one;
      }
      bv *= b_frac;
    }
    av *= a_frac;
  }

  // Apply the chain rules to reverse the transformation and return to expression in Cartesian
  // space.
  //
  // a = (umat[0] * x) + (umat[3] * y) + (umat[6] * z)
  // b =                 (umat[4] * y) + (umat[7] * z)
  // c =                                 (umat[8] * z)
  //
  // du/dx = (du/da)(da/dx) + (du/db)(db/dx) + (du/dc)(dc/dx)
  //       = umat[0] * (du/da)
  // du/dy = (du/da)(da/dy) + (du/db)(db/dy) + (du/dc)(dc/dy)
  //       = (umat[3] * (du/da)) + (umat[4] * (du/db))
  // du/dz = (du/da)(da/dz) + (du/db)(db/dz) + (du/dc)(dc/dz)
  //       = (umat[6] * (du/da)) + (umat[7] * (du/db)) + (umat[8] * (du/dc))
  const T3 result = { (umat[0] * result_da),
                      (umat[3] * result_da) + (umat[4] * result_db),
                      (umat[6] * result_da) + (umat[7] * result_db) + (umat[8] * result_dc) };
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TricubicCell<T>::secondDerivative(const T x, const T y, const T z, T* result) const {
  T a_frac, b_frac, c_frac;
  fractionalPosition(x, y, z, &a_frac, &b_frac, &c_frac, ExceptionResponse::DIE, "evaluate");
  const T value_zero = 0.0;
  const T value_one  = 1.0;
  const T d2factors[4] = { 0.0, 0.0, 2.0, 6.0 };
  const T dfactors[4]  = { 0.0, 1.0, 2.0, 3.0 };
  const bool a_is_zero = (std::abs(a_frac) < constants::tiny);
  const bool b_is_zero = (std::abs(b_frac) < constants::tiny);
  const bool c_is_zero = (std::abs(c_frac) < constants::tiny);
  
  // Compute the matrix of second-order partial derivatives in the mesh element coordinate space
  T result_daa = value_zero;
  T result_dab = value_zero;
  T result_dac = value_zero;
  T result_dbb = value_zero;
  T result_dbc = value_zero;
  T result_dcc = value_zero;
  T av    = value_one;
  T d_av  = (a_is_zero) ? value_zero : value_one / a_frac;
  T d2_av = d_av * d_av;
  for (int i = 0; i < 4; i++) {
    T bv    = value_one;
    T d_bv  = (b_is_zero) ? value_zero : value_one / b_frac;
    T d2_bv = d_bv * d_bv;
    for (int j = 0; j < 4; j++) {
      T cv    = value_one;
      T d_cv  = (c_is_zero) ? value_zero : value_one / c_frac;
      T d2_cv = d_cv * d_cv;
      for (int k = 0; k < 4; k++) {
        const T coef = coefficients[(4 * ((4 * k) + j)) + i];
        result_daa += coef * d2factors[i] * d2_av * bv * cv;
        result_dab += coef * dfactors[i] * d_av * dfactors[j] * d_bv * cv;
        result_dac += coef * dfactors[i] * d_av * bv * dfactors[k] * d_cv;
        result_dbb += coef * av * d2factors[j] * d2_bv * cv;
        result_dbc += coef * av * dfactors[j] * d_bv * dfactors[k] * d_cv;
        result_dcc += coef * av * bv * d2factors[k] * d2_cv;
        cv *= c_frac;
        if (c_is_zero) {
          d_cv  = (k == 0) ? value_one : value_zero;
          d2_cv = (k == 1) ? value_one : value_zero;
        }
        else {
          d_cv  *= c_frac;
          d2_cv *= c_frac;
        }
      }
      bv *= b_frac;
      if (b_is_zero) {
        d_bv  = (j == 0) ? value_one : value_zero;
        d2_bv = (j == 1) ? value_one : value_zero;
      }
      else {
        d_bv  *= b_frac;
        d2_bv *= b_frac;
      }
    }
    av *= a_frac;
    if (a_is_zero) {
      d_av  = (i == 0) ? value_one : value_zero;
      d2_av = (i == 1) ? value_one : value_zero;
    }
    else {
      d_av  *= a_frac;
      d2_av *= a_frac;
    }
  }
  
  // Follow the chain rules for a(x, y, z), b(x, y, z), and c(x, y, z) when evaluating derivatives
  // of f(a, b, c).  This parallels operations to map derivatives of the function with respect to
  // x, y, and z to the mesh.  Take U as the transformation matrix taking coordinates into the
  // mesh element space (U is the member variable umat):
  //
  // a = (U(1,1) * x) + (U(1,2) * y) + (U(1,3) * z)
  // b =                (U(2,2) * y) + (U(2,3) * z)
  // c =                               (U(3,3) * z)
  //
  // df(a,b,c) / dx = (df / da)(da / dx) + (df / db)(db / dx) + (df / dc)(dc / dx)
  // df(a,b,c) / dy = (df / da)(da / dy) + (df / db)(db / dy) + (df / dc)(dc / dy)
  // df(a,b,c) / dz = (df / da)(da / dz) + (df / db)(db / dz) + (df / dc)(dc / dz)
  // df(a,b,c) / dx = U(1,1)*(df / da)
  // df(a,b,c) / dy = U(1,2)*(df / da) + U(2,2)*(df / db)
  // df(a,b,c) / dz = U(1,3)*(df / da) + U(2,3)*(df / db) + U(3,3)*(df / dc)
  //
  // (d/dx)(df/dx) = (d/dx) [ U(1,1)*(df/da) ]
  //               = U(1,1) * [ (d2f/da2)(da/dx) + (d2f/dab)(db/dx) + (d2f/dac)(dc/dx) ]
  //               = U(1,1)*U(1,1)*(d2f/da2)
  // (d/dy)(df/dx) = (d/dy) [ U(1,1)*(df/da) ]
  //               = U(1,1) * [ (d2f/da2)(da/dy) + (d2f/dab)(db/dy) + (d2f/dac)(dc/dy) ]
  //               = U(1,1) * [ U(1,2)*(d2f / da2) + U(2,2)*(d2f / dab) ]
  // (d/dz)(df/dx) = (d/dz) [ U(1,1)*(df/da) ]
  //               = U(1,1) * [ (d2f/da2)(da/dz) + (d2f/dab)(db/dz) + (d2f/dac)(dc/dz) ]
  //               = U(1,1) * [ U(1,3)*(d2f / da2) + U(2,3)*(d2f / dab) + U(3,3)*(d2f / dac) ]
  // (d/dy)(df/dy) = (d/dy) [ U(1,2)*(df / da) + U(2,2)*(df / db) ]
  //               = U(1,2) * [ (d/dy)(df/da) ] + U(2,2) * [ (d/dy)(df/db) ]
  //               = U(1,2) * [ (d2f/da2)(da/dy) + (d2f/dab)(db/dy) + (d2f/dac)(dc/dy) ] +
  //                 U(2,2) * [ (d2f/dab)(da/dy) + (d2f/db2)(db/dy) + (d2f/dbc)(dc/dy) ]
  //               = U(1,2) * [ U(1,2)*(d2f / da2) + U(2,2)*(d2f / dab) ] +
  //                 U(2,2) * [ U(1,2)*(d2f / dab) + U(2,2)*(d2f / db2) ]
  // (d/dz)(df/dy) = (d/dz) [ U(1,2)*(df / da) + U(2,2)*(df / db) ]
  //               = U(1,2) * [ (d/dz)(df/da) ] + U(2,2) * [ (d/dz)(df/db) ]
  //               = U(1,2) * [ (d2f/da2)(da/dz) + (d2f/dab)(db/dz) + (d2f/dac)(dc/dz) ] +
  //                 U(2,2) * [ (d2f/dab)(da/dz) + (d2f/db2)(db/dz) + (d2f/dbc)(dc/dz) ]
  //               = U(1,2) * [ U(1,3)*(d2f / da2) + U(2,3)*(d2f / dab) + U(3,3)*(d2f / dac) ] +
  //                 U(2,2) * [ U(1,3)*(d2f / dab) + U(2,3)*(d2f / db2) + U(3,3)*(d2f / dbc) ]
  // (d/dz)(df/dz) = (d/dz) [ U(1,3)*(df / da) + U(2,3)*(df / db) + U(3,3)*(df / dc) ]
  //               = U(1,3) * [ (d/dz)(df/da) ] + U(2,3) * [ (d/dz)(df/db) ] +
  //                 U(3,3) * [ (d/dz)(df/dc) ]
  //               = U(1,3) * [ (d2f/da2)(da/dz) + (d2f/dab)(db/dz) + (d2f/dac)(dc/dz) ] +
  //                 U(2,3) * [ (d2f/dab)(da/dz) + (d2f/db2)(db/dz) + (d2f/dbc)(dc/dz) ] +
  //                 U(3,3) * [ (d2f/dac)(da/dz) + (d2f/dbc)(db/dz) + (d2f/dc2)(dc/dz) ]
  //               = U(1,3) * [ U(1,3)*(d2f / da2) + U(2,3)*(d2f / dab) + U(3,3)*(d2f / dac) ] +
  //                 U(2,3) * [ U(1,3)*(d2f / dab) + U(2,3)*(d2f / db2) + U(3,3)*(d2f / dbc) ] +
  //                 U(3,3) * [ U(1,3)*(d2f / dac) + U(2,3)*(d2f / dbc) + U(3,3)*(d2f / dc2) ]
  const T result_dxx = umat[0] * umat[0] * result_daa;
  const T result_dxy = umat[0] * ((umat[3] * result_daa) + (umat[4] * result_dab));
  const T top_row = (umat[6] * result_daa) + (umat[7] * result_dab) + (umat[8] * result_dac);
  const T result_dxz = umat[0] * top_row;
  const T result_dyy = (umat[3] * ((umat[3] * result_daa) + (umat[4] * result_dab))) +
                       (umat[4] * ((umat[3] * result_dab) + (umat[4] * result_dbb)));
  const T mid_row = (umat[6] * result_dab) + (umat[7] * result_dbb) + (umat[8] * result_dbc);
  const T result_dyz = (umat[3] * top_row) + (umat[4] * mid_row);
  const T result_dzz = (umat[6] * top_row) + (umat[7] * mid_row) +
                       (umat[8] * ((umat[6] * result_dac) + (umat[7] * result_dbc) +
                                   (umat[8] * result_dcc)));
  result[0] = result_dxx;
  result[1] = result_dxy;
  result[2] = result_dxz;
  result[3] = result_dxy;
  result[4] = result_dyy;
  result[5] = result_dyz;
  result[6] = result_dxz;
  result[7] = result_dyz;
  result[8] = result_dzz;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TricubicCell<T>::secondDerivative(const T x, const T y, const T z,
                                       std::vector<T> *result) const {
  if (result->size() != 9) {
    rtErr("The result must hold nine elements of the second-order partial derivative matrix "
          "(size of provided vector " + std::to_string(result->size()) + ").", "TricubicCell",
          "secondDerivative");
  }
  secondDerivative(x, y, z, result->data());
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> TricubicCell<T>::secondDerivative(const T x, const T y, const T z) const {
  std::vector<T> result(9);
  secondDerivative(x, y, z, result.data());
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TricubicCell<T>::thirdDerivative(const T x, const T y, const T z, T* result) const {
  T a_frac, b_frac, c_frac;
  fractionalPosition(x, y, z, &a_frac, &b_frac, &c_frac, ExceptionResponse::DIE, "evaluate");
  const T value_zero = 0.0;
  const T value_one  = 1.0;
  const T d3factors[4] = { 0.0, 0.0, 0.0, 6.0 };
  const T d2factors[4] = { 0.0, 0.0, 2.0, 6.0 };
  const T dfactors[4]  = { 0.0, 1.0, 2.0, 3.0 };
  const bool a_is_zero = (std::abs(a_frac) < constants::tiny);
  const bool b_is_zero = (std::abs(b_frac) < constants::tiny);
  const bool c_is_zero = (std::abs(c_frac) < constants::tiny);
  
  // Compute the matrix of second-order partial derivatives in the mesh element coordinate space
  T result_daaa = value_zero;
  T result_daab = value_zero;
  T result_daac = value_zero;
  T result_dabb = value_zero;
  T result_dabc = value_zero;
  T result_dacc = value_zero;
  T result_dbbb = value_zero;
  T result_dbbc = value_zero;
  T result_dbcc = value_zero;
  T result_dccc = value_zero;
  T av    = value_one;
  T d_av  = (a_is_zero) ? value_zero : value_one / a_frac;
  T d2_av = d_av * d_av;
  T d3_av = d2_av * d_av;
  for (int i = 0; i < 4; i++) {
    T bv    = value_one;
    T d_bv  = (b_is_zero) ? value_zero : value_one / b_frac;
    T d2_bv = d_bv * d_bv;
    T d3_bv = d2_bv * d_bv;
    for (int j = 0; j < 4; j++) {
      T cv    = value_one;
      T d_cv  = (c_is_zero) ? value_zero : value_one / c_frac;
      T d2_cv = d_cv * d_cv;
      T d3_cv = d2_cv * d_cv;
      for (int k = 0; k < 4; k++) {
        const T coef = coefficients[(4 * ((4 * k) + j)) + i];
        result_daaa += coef * d3factors[i] * d3_av * bv * cv;
        result_daab += coef * d2factors[i] * d2_av * dfactors[j] * d_bv * cv;
        result_daac += coef * d2factors[i] * d2_av * bv * dfactors[k] * d_cv;
        result_dabb += coef * dfactors[i] * d_av * d2factors[j] * d2_bv * cv;
        result_dabc += coef * dfactors[i] * d_av * dfactors[j] * d_bv * dfactors[k] * d_cv;
        result_dacc += coef * dfactors[i] * d_av * bv * d2factors[k] * d2_cv;
        result_dbbb += coef * av * d3factors[j] * d3_bv * cv;
        result_dbbc += coef * av * d2factors[j] * d2_bv * dfactors[k] * d_cv;
        result_dbcc += coef * av * dfactors[j] * d_bv * d2factors[k] * d2_cv;
        result_dccc += coef * av * bv * d3factors[k] * d3_cv;
        cv *= c_frac;
        if (c_is_zero) {
          d_cv  = (k == 0) ? value_one : value_zero;
          d2_cv = (k == 1) ? value_one : value_zero;
          d3_cv = (k == 2) ? value_one : value_zero;
        }
        else {
          d_cv  *= c_frac;
          d2_cv *= c_frac;
          d3_cv *= c_frac;
        }
      }
      bv *= b_frac;
      if (b_is_zero) {
        d_bv  = (j == 0) ? value_one : value_zero;
        d2_bv = (j == 1) ? value_one : value_zero;
        d3_bv = (j == 2) ? value_one : value_zero;
      }
      else {
        d_bv  *= b_frac;
        d2_bv *= b_frac;
        d3_bv *= b_frac;
      }
    }
    av *= a_frac;
    if (a_is_zero) {
      d_av  = (i == 0) ? value_one : value_zero;
      d2_av = (i == 1) ? value_one : value_zero;
      d3_av = (i == 2) ? value_one : value_zero;
    }
    else {
      d_av  *= a_frac;
      d2_av *= a_frac;
      d3_av *= a_frac;
    }
  }
  
  // Follow the chain rules for a(x, y, z), b(x, y, z), and c(x, y, z) when evaluating derivatives
  // of f(a, b, c).  This parallels operations to map derivatives of the function with respect to
  // x, y, and z to the mesh, and continues from the work in secondDerivative() above.
  //
  // (d/dx)(d2f/dx2) = (d/dx) [ U(1,1)*U(1,1)*(d2f/da2) ]
  //                 = U(1,1)*U(1,1) * [ (d3f/da3)(da/dx) + (d3f/daab)(db/dx) + (d3f/daac)(dc/dx) ]
  //                 = U(1,1)*U(1,1)*U(1,1)*(d3f/da3)
  const T result_dxxx = umat[0] * umat[0] * umat[0] * result_daaa;

  // (d/dx)(d2f/dxy) = (d/dx) [ U(1,1) * [ U(1,2)*(d2f / da2) + U(2,2)*(d2f / dab) ] ]
  //                 = U(1,1) * [ U(1,2) * [ (d3f/da3)(da/dx)  + (d3f/daab)(db/dx) +
  //                                         (d3f/daac)(dc/dx) ] +
  //                              U(2,2) * [ (d3f/daab)(da/dx) + (d3f/dabb)(db/dx) +
  //                                         (d3f/dabc)(dc/dx) ] ]
  //                 = U(1,1) * [ U(1,2)*U(1,1)*(d3f / da3) + U(2,2)*U(1,1)*(d3f / daab) ]
  //                 = U(1,1) * U(1,1) * [ U(1,2)*(d3f / da3) + U(2,2)*(d3f / daab) ]
  const T result_dxxy = umat[0] * umat[0] * ((umat[3] * result_daaa) + (umat[4] * result_daab));
  
  // (d/dx)(d2f/dxz) = (d/dx) [ U(1,1) * [ U(1,3)*(d2f / da2) + U(2,3)*(d2f / dab) +
  //                                       U(3,3)*(d2f / dac) ] ]
  //                 = U(1,1) * [ U(1,3) * [ (d3f/da3)(da/dx)  + (d3f/daab)(db/dx) +
  //                                         (d3f/daac)(dc/dx) ] +
  //                              U(2,3) * [ (d3f/daab)(da/dx) + (d3f/dabb)(db/dx) +
  //                                         (d3f/dabc)(dc/dx) ] +
  //                              U(3,3) * [ (d3f/daac)(da/dx) + (d3f/dabc)(db/dx) +
  //                                         (d3f/dacc)(dc/dx) ] ]
  //                 = U(1,1) * U(1,1) * [ U(1,3)*(d3f / da3) + U(2,3)*(d3f / daab) +
  //                                       U(3,3)*(d3f / daac) ]
  const T result_dxxz = umat[0] * umat[0] * ((umat[6] * result_daaa) + (umat[7] * result_daab) +
                                             (umat[8] * result_daac));

  // (d/dy)(d2f/dxy) = (d/dy) [ U(1,1) * [ U(1,2)*(d2f / da2) + U(2,2)*(d2f / dab) ] ]
  //                 = U(1,1) * [ U(1,2) * [ (d3f/da3)(da/dy)  + (d3f/daab)(db/dy) +
  //                                         (d3f/daac)(dc/dy) ] +
  //                              U(2,2) * [ (d3f/daab)(da/dy) + (d3f/dabb)(db/dy) +
  //                                         (d3f/dabc)(dc/dy) ] ]
  //                 = U(1,1) * [ U(1,2) * [ U(1,2)*(d3f / da3)  + U(2,2)*(d3f / daab) ] +
  //                              U(2,2) * [ U(1,2)*(d3f / daab) + U(2,2)*(d3f / dabb) ] ]
  const T result_dxyy = umat[0] *
                        ((umat[3] * ((umat[3] * result_daaa) + (umat[4] * result_daab))) +
                         (umat[4] * ((umat[3] * result_daab) + (umat[4] * result_dabb))));

  // (d/dy)(d2f/dxz) = (d/dy) [ U(1,1) * [ U(1,3)*(d2f / da2) + U(2,3)*(d2f / dab) +
  //                                       U(3,3)*(d2f / dac) ] ]
  //                 = U(1,1) * [ U(1,3) * [ (d3f/da3)(da/dy)  + (d3f/daab)(db/dy) +
  //                                         (d3f/daac)(dc/dy) ] +
  //                              U(2,3) * [ (d3f/daab)(da/dy) + (d3f/dabb)(db/dy) +
  //                                         (d3f/dabc)(dc/dy) ] +
  //                              U(3,3) * [ (d3f/daac)(da/dy) + (d3f/dabc)(db/dy) +
  //                                         (d3f/dacc)(dc/dy) ] ]
  //                 = U(1,1) * [ U(1,3) * [ U(1,2)*(d3f / da3)  + U(2,2)*(d3f / daab) ] +
  //                              U(2,3) * [ U(1,2)*(d3f / daab) + U(2,2)*(d3f / dabb) ] +
  //                              U(3,3) * [ U(1,2)*(d3f / daac) + U(2,2)*(d3f / dabc) ]
  const T result_dxyz = umat[0] *
                        ((umat[6] * ((umat[3] * result_daaa) + (umat[4] * result_daab))) +
                         (umat[7] * ((umat[3] * result_daab) + (umat[4] * result_dabb))) +
                         (umat[8] * ((umat[3] * result_daac) + (umat[4] * result_dabc))));

  // (d/dz)(d2f/dxz) = (d/dz) [ U(1,1) * [ U(1,3)*(d2f / da2) + U(2,3)*(d2f / dab) +
  //                                       U(3,3)*(d2f / dac) ] ]
  //                 = U(1,1) * [ U(1,3) * [ (d3f/da3)(da/dz)  + (d3f/daab)(db/dz) +
  //                                         (d3f/daac)(dc/dz) ] +
  //                              U(2,3) * [ (d3f/daab)(da/dz) + (d3f/dabb)(db/dz) +
  //                                         (d3f/dabc)(dc/dz) ] +
  //                              U(3,3) * [ (d3f/daac)(da/dz) + (d3f/dabc)(db/dz) +
  //                                         (d3f/dacc)(dc/dz) ] ]
  //                 = U(1,1) * [ U(1,3) * [ U(1,3)*(d3f / da3)  + U(2,3)*(d3f / daab) +
  //                                         U(3,3)*(d3f / daac) ] +
  //                              U(2,3) * [ U(1,3)*(d3f / daab) + U(2,3)*(d3f / dabb) +
  //                                         U(3,3)*(d3f / dabc) ] +
  //                              U(3,3) * [ U(1,3)*(d3f / daac) + U(2,3)*(d3f / dabc) +
  //                                         U(3,3)*(d3f / dacc) ] ]
  const T result_dxzz = umat[0] *
                        ((umat[6] * ((umat[6] * result_daaa) + (umat[7] * result_daab) +
                                     (umat[8] * result_daac))) +
                         (umat[7] * ((umat[6] * result_daab) + (umat[7] * result_dabb) +
                                     (umat[8] * result_dabc))) +
                         (umat[8] * ((umat[6] * result_daac) + (umat[7] * result_dabc) +
                                     (umat[8] * result_dacc))));
  
  // (d/dy)(d2f/dy2) = (d/dy) [ U(1,2) * [ U(1,2)*(d2f / da2) + U(2,2)*(d2f / dab) ] +
  //                            U(2,2) * [ U(1,2)*(d2f / dab) + U(2,2)*(d2f / db2) ] ]
  //                 = U(1,2) * [ U(1,2) * [ (d3f/da3)(da/dy)  + (d3f/daab)(db/dy) +
  //                                         (d3f/daac)(dc/dy) ] +
  //                              U(2,2) * [ (d3f/daab)(da/dy) + (d3f/dabb)(db/dy) +
  //                                         (d3f/dabc)(dc/dy) ] ] +
  //                   U(2,2) * [ U(1,2) * [ (d3f/daab)(da/dy) + (d3f/dabb)(db/dy) +
  //                                         (d3f/dabc)(dc/dy) ] +
  //                              U(2,2) * [ (d3f/dabb)(da/dy) + (d3f/db3)(db/dy) +
  //                                         (d3f/dbbc)(dc/dy) ] ]
  //                 = U(1,2) * [ U(1,2) * [ U(1,2)*(d3f / da3)  + U(2,2)*(d3f / daab) ] +
  //                              U(2,2) * [ U(1,2)*(d3f / daab) + U(2,2)*(d3f / dabb) ] ] +
  //                   U(2,2) * [ U(1,2) * [ U(1,2)*(d3f / daab) + U(2,2)*(d3f / dabb) ] +
  //                              U(2,2) * [ U(1,2)*(d3f / dabb) + U(2,2)*(d3f / db3)  ] ]
  const T result_dyyy = (umat[3] *
                         ((umat[3] * ((umat[3] * result_daaa) + (umat[4] * result_daab))) +
                          (umat[4] * ((umat[3] * result_daab) + (umat[4] * result_dabb))))) +
                        (umat[4] *
                         ((umat[3] * ((umat[3] * result_daab) + (umat[4] * result_dabb))) +
                          (umat[4] * ((umat[3] * result_dabb) + (umat[4] * result_dbbb)))));
                         
  // (d/dz)(d2f/dy2) = (d/dz) [ U(1,2) * [ U(1,2)*(d2f / da2) + U(2,2)*(d2f / dab) ] +
  //                            U(2,2) * [ U(1,2)*(d2f / dab) + U(2,2)*(d2f / db2) ] ]
  //                 = U(1,2) * [ U(1,2) * [ (d3f/da3)(da/dz)  + (d3f/daab)(db/dz) +
  //                                         (d3f/daac)(dc/dz) ] +
  //                              U(2,2) * [ (d3f/daab)(da/dz) + (d3f/dabb)(db/dz) +
  //                                         (d3f/dabc)(dc/dz) ] ] +
  //                   U(2,2) * [ U(1,2) * [ (d3f/daab)(da/dz) + (d3f/dabb)(db/dz) +
  //                                         (d3f/dabc)(dc/dz) ] +
  //                              U(2,2) * [ (d3f/dabb)(da/dz) + (d3f/db3)(db/dz) +
  //                                         (d3f/dbbc)(dc/dz) ] ]
  //                 = U(1,2) * [ U(1,2) * [ U(1,3)*(d3f / da3)  + U(2,3)*(d3f / daab) +
  //                                         U(3,3)*(d3f / daac) ] +
  //                              U(2,2) * [ U(1,3)*(d3f / daab) + U(2,3)*(d3f / dabb) +
  //                                         U(3,3)*(d3f / dabc) ] ] +
  //                   U(2,2) * [ U(1,2) * [ U(1,3)*(d3f / daab) + U(2,3)*(d3f / dabb) +
  //                                         U(3,3)*(d3f / dabc) ] +
  //                              U(2,2) * [ U(1,3)*(d3f / dabb) + U(2,3)*(d3f / db3)  +
  //                                         U(3,3)*(d3f / dbbc) ] ]
  const T result_dyyz = (umat[3] * ((umat[3] * ((umat[6] * result_daaa) + (umat[7] * result_daab) +
                                                (umat[8] * result_daac))) +
                                    (umat[4] * ((umat[6] * result_daab) + (umat[7] * result_dabb) +
                                                (umat[8] * result_dabc))))) +
                        (umat[4] * ((umat[3] * ((umat[6] * result_daab) + (umat[7] * result_dabb) +
                                                (umat[8] * result_dabc))) +
                                    (umat[4] * ((umat[6] * result_dabb) + (umat[7] * result_dbbb) +
                                                (umat[8] * result_dbbc)))));

  // (d/dy)(d2f/dz2) = (d/dy) [ U(1,3) * [ U(1,3)*(d2f / da2) + U(2,3)*(d2f / dab) +
  //                                       U(3,3)*(d2f / dac) ] +    
  //                            U(2,3) * [ U(1,3)*(d2f / dab) + U(2,3)*(d2f / db2) +
  //                                       U(3,3)*(d2f / dbc) ] +    
  //                            U(3,3) * [ U(1,3)*(d2f / dac) + U(2,3)*(d2f / dbc) +
  //                                       U(3,3)*(d2f / dc2) ]
  //                 = U(1,3) * [ U(1,3) * [ (d3f/da3)(da/dy)  + (d3f/daab)(db/dy) +
  //                                         (d3f/daac)(dc/dy) ] +
  //                              U(2,3) * [ (d3f/daab)(da/dy) + (d3f/dabb)(db/dy) +
  //                                         (d3f/dabc)(dc/dy) ] +
  //                              U(3,3) * [ (d3f/daac)(da/dy) + (d3f/dabc)(db/dy) +
  //                                         (d3f/dacc)(dc/dy) ] ] +
  //                   U(2,3) * [ U(1,3) * [ (d3f/daab)(da/dy) + (d3f/dabb)(db/dy) +
  //                                         (d3f/dabc)(dc/dy) ] +
  //                              U(2,3) * [ (d3f/dabb)(da/dy) + (d3f/db3)(db/dy)  +
  //                                         (d3f/dbbc)(dc/dy) ] +
  //                              U(3,3) * [ (d3f/dabc)(da/dy) + (d3f/dbbc)(db/dy) +
  //                                         (d3f/dbcc)(dc/dy) ] ] +
  //                   U(3,3) * [ U(1,3) * [ (d3f/daac)(da/dy) + (d3f/dabc)(db/dy) +
  //                                         (d3f/dacc)(dc/dy) ] +
  //                              U(2,3) * [ (d3f/dabc)(da/dy) + (d3f/dbbc)(db/dy) +
  //                                         (d3f/dbcc)(dc/dy) ] +
  //                              U(3,3) * [ (d3f/dacc)(da/dy) + (d3f/dbcc)(db/dy) +
  //                                         (d3f/dc3)(dc/dy)  ] ]
  //                 = U(1,3) * [ U(1,3) * [ U(1,2)*(d3f / da3)  + U(2,2)*(d3f / daab) ] +
  //                              U(2,3) * [ U(1,2)*(d3f / daab) + U(2,2)*(d3f / dabb) ] +
  //                              U(3,3) * [ U(1,2)*(d3f / daac) + U(2,2)*(d3f / dabc) ] ] +
  //                   U(2,3) * [ U(1,3) * [ U(1,2)*(d3f / daab) + U(2,2)*(d3f / dabb) ] +
  //                              U(2,3) * [ U(1,2)*(d3f / dabb) + U(2,2)*(d3f / db3)  ] +
  //                              U(3,3) * [ U(1,2)*(d3f / dabc) + U(2,2)*(d3f / dbbc) ] ] +
  //                   U(3,3) * [ U(1,3) * [ U(1,2)*(d3f / daac) + U(2,2)*(d3f / dabc) ] +
  //                              U(2,3) * [ U(1,2)*(d3f / dabc) + U(2,2)*(d3f / dbbc) ] +
  //                              U(3,3) * [ U(1,2)*(d3f / dacc) + U(2,2)*(d3f / dbcc) ] ]
  const T result_dyzz = (umat[6] *
                         ((umat[6] * ((umat[3] * result_daaa) + (umat[4] * result_daab))) +
                          (umat[7] * ((umat[3] * result_daab) + (umat[4] * result_dabb))) +
                          (umat[8] * ((umat[3] * result_daac) + (umat[4] * result_dabc))))) +
                        (umat[7] *
                         ((umat[6] * ((umat[3] * result_daab) + (umat[4] * result_dabb))) +
                          (umat[7] * ((umat[3] * result_dabb) + (umat[4] * result_dbbb))) +
                          (umat[8] * ((umat[3] * result_dabc) + (umat[4] * result_dbbc))))) +
                        (umat[8] *
                         ((umat[6] * ((umat[3] * result_daac) + (umat[4] * result_dabc))) +
                          (umat[7] * ((umat[3] * result_dabc) + (umat[4] * result_dbbc))) +
                          (umat[8] * ((umat[3] * result_dacc) + (umat[4] * result_dbcc)))));

  // (d/dz)(d2fdz2) = (d/dz) [ U(1,3) * [ U(1,3)*(d2f / da2) + U(2,3)*(d2f / dab) +
  //                                      U(3,3)*(d2f / dac) ] +    
  //                           U(2,3) * [ U(1,3)*(d2f / dab) + U(2,3)*(d2f / db2) +
  //                                      U(3,3)*(d2f / dbc) ] +    
  //                           U(3,3) * [ U(1,3)*(d2f / dac) + U(2,3)*(d2f / dbc) +
  //                                      U(3,3)*(d2f / dc2) ]
  //                = U(1,3) * [ U(1,3) * [ (d3f/da3)(da/dz)  + (d3f/daab)(db/dz) +
  //                                        (d3f/daac)(dc/dz) ] +
  //                             U(2,3) * [ (d3f/daab)(da/dz) + (d3f/dabb)(db/dz) +
  //                                        (d3f/dabc)(dc/dz) ] +
  //                             U(3,3) * [ (d3f/daac)(da/dz) + (d3f/dabc)(db/dz) +
  //                                        (d3f/dacc)(dc/dz) ] ] +
  //                  U(2,3) * [ U(1,3) * [ (d3f/daab)(da/dz) + (d3f/dabb)(db/dz) +
  //                                        (d3f/dabc)(dc/dz) ] +
  //                             U(2,3) * [ (d3f/dabb)(da/dz) + (d3f/db3)(db/dz)  +
  //                                        (d3f/dbbc)(dc/dz) ] +
  //                             U(3,3) * [ (d3f/dabc)(da/dz) + (d3f/dbbc)(db/dz) +
  //                                        (d3f/dbcc)(dc/dz) ] ] +
  //                  U(3,3) * [ U(1,3) * [ (d3f/daac)(da/dz) + (d3f/dabc)(db/dz) +
  //                                        (d3f/dacc)(dc/dz) ] +
  //                             U(2,3) * [ (d3f/dabc)(da/dz) + (d3f/dbbc)(db/dz) +
  //                                        (d3f/dbcc)(dc/dz) ] +
  //                             U(3,3) * [ (d3f/dacc)(da/dz) + (d3f/dbcc)(db/dz) +
  //                                        (d3f/dc3)(dc/dz)  ] ]
  //                = U(1,3) * [ U(1,3) * [ U(1,3)*(d3f / da3)  + U(2,3)*(d3f / daab) +
  //                                        U(3,3)*(d3f / daac) ] +
  //                             U(2,3) * [ U(1,3)*(d3f / daab) + U(2,3)*(d3f / dabb) ] +
  //                                        U(3,3)*(d3f / dabc) ] +
  //                             U(3,3) * [ U(1,3)*(d3f / daac) + U(2,3)*(d3f / dabc) ] +
  //                                        U(3,3)*(d3f / dacc) ] ] +
  //                  U(2,3) * [ U(1,3) * [ U(1,3)*(d3f / daab) + U(2,3)*(d3f / dabb) ] +
  //                                        U(3,3)*(d3f / dabc) ] +
  //                             U(2,3) * [ U(1,3)*(d3f / dabb) + U(2,3)*(d3f / db3)  ] +
  //                                        U(3,3)*(d3f / dbbc) ] +
  //                             U(3,3) * [ U(1,3)*(d3f / dabc) + U(2,3)*(d3f / dbbc) ] +
  //                                        U(3,3)*(d3f / dbcc) ] ] +
  //                  U(3,3) * [ U(1,3) * [ U(1,3)*(d3f / daac) + U(2,3)*(d3f / dabc) ] +
  //                                        U(3,3)*(d3f / dacc) ] +
  //                             U(2,3) * [ U(1,3)*(d3f / dabc) + U(2,3)*(d3f / dbbc) ] +
  //                                        U(3,3)*(d3f / dbcc) ] +
  //                             U(3,3) * [ U(1,3)*(d3f / dacc) + U(2,3)*(d3f / dbcc) ] +
  //                                        U(3,3)*(d3f / dccc) ] ]
  const T result_dzzz = (umat[6] *
                         ((umat[6] * ((umat[6] * result_daaa) + (umat[7] * result_daab) +
                                      (umat[8] * result_daac))) +
                          (umat[7] * ((umat[6] * result_daab) + (umat[7] * result_dabb) +
                                      (umat[8] * result_dabc))) +
                          (umat[8] * ((umat[6] * result_daac) + (umat[7] * result_dabc) +
                                      (umat[8] * result_dacc))))) +
                        (umat[7] *
                         ((umat[6] * ((umat[6] * result_daab) + (umat[7] * result_dabb) +
                                      (umat[8] * result_dabc))) +
                          (umat[7] * ((umat[6] * result_dabb) + (umat[7] * result_dbbb) +
                                      (umat[8] * result_dbbc))) +
                          (umat[8] * ((umat[6] * result_dabc) + (umat[7] * result_dbbc) +
                                      (umat[8] * result_dbcc))))) +
                        (umat[8] *
                         ((umat[6] * ((umat[6] * result_daac) + (umat[7] * result_dabc) +
                                      (umat[8] * result_dacc))) +
                          (umat[7] * ((umat[6] * result_dabc) + (umat[7] * result_dbbc) +
                                      (umat[8] * result_dbcc))) +
                          (umat[8] * ((umat[6] * result_dacc) + (umat[7] * result_dbcc) +
                                      (umat[8] * result_dccc)))));

  // Transcribe the results into the tensor array.  For the ith, jth, and kth derivatives
  // df / dijk and numbering x, y, and z as 0, 1, and 2, the value of the triple derivative
  // df / dijk is found in the tensor array at element (((k * 3) + j) * 3) + i.
  result[ 0] = result_dxxx;
  result[ 1] = result_dxxy;
  result[ 2] = result_dxxz;
  result[ 3] = result_dxxy;
  result[ 4] = result_dxyy;
  result[ 5] = result_dxyz;
  result[ 6] = result_dxxz;
  result[ 7] = result_dxyz;
  result[ 8] = result_dxzz;
  result[ 9] = result_dxxy;
  result[10] = result_dxyy;
  result[11] = result_dxyz;
  result[12] = result_dxyy;
  result[13] = result_dyyy;
  result[14] = result_dyyz;
  result[15] = result_dxyz;
  result[16] = result_dyyz;
  result[17] = result_dyzz;
  result[18] = result_dxxz;
  result[19] = result_dxyz;
  result[20] = result_dxzz;
  result[21] = result_dxyz;
  result[22] = result_dyyz;
  result[23] = result_dyzz;
  result[24] = result_dxzz;
  result[25] = result_dyzz;
  result[26] = result_dzzz;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void TricubicCell<T>::thirdDerivative(const T x, const T y, const T z,
                                      std::vector<T> *result) const {
  if (result->size() != 27) {
    rtErr("The result must hold nine elements of the second-order partial derivative matrix "
          "(size of provided vector " + std::to_string(result->size()) + ").", "TricubicCell",
          "thirdDerivative");
  }
  thirdDerivative(x, y, z, result->data());
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> TricubicCell<T>::thirdDerivative(const T x, const T y, const T z) const {
  std::vector<T> result(27);
  thirdDerivative(x, y, z, result.data());
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void incorporateStencilOrigin(const T orig_x, const T orig_y, const T orig_z, const double3 pt_xyz,
                              const double scale_factor, T *point_x, T *point_y, T *point_z) {
  if (isFloatingPointScalarType<T>()) {
    if (fabs(scale_factor - 1.0) > constants::tiny) {
      rtErr("A scaling factor of " + realToString(scale_factor, 9, 2, NumberFormat::SCIENTIFIC) +
            " is incompatible with a floating point coordinate representations.");
    }
    *point_x = orig_x + pt_xyz.x;
    *point_y = orig_y + pt_xyz.y;
    *point_z = orig_z + pt_xyz.z;
  }
  else {
    const llint x_i = llround(pt_xyz.x * scale_factor);
    const llint y_i = llround(pt_xyz.y * scale_factor);
    const llint z_i = llround(pt_xyz.z * scale_factor);
    *point_x = static_cast<llint>(orig_x) + x_i;
    *point_y = static_cast<llint>(orig_y) + y_i;
    *point_z = static_cast<llint>(orig_z) + z_i;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void fvStencilCoordinates(const T orig_x, const T orig_y, const T orig_z,
                          const double scale_factor, const UnitCellAxis face_normal,
                          const int a_index, const int b_index, const int c_index,
                          const double* invu, T *point_x, T *point_y, T *point_z) {
  const double3 pt_xyz = stencilInternalOffset(face_normal, a_index, b_index, c_index, invu);
  incorporateStencilOrigin(orig_x, orig_y, orig_z, pt_xyz, scale_factor, point_x, point_y,
                           point_z);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void fvStencilCoordinates(const T orig_x, const T orig_y, const T orig_z,
                          const double scale_factor, const int a_index, const int b_index,
                          const int c_index, const double* invu, T *point_x, T *point_y,
                          T *point_z) {
  const double3 pt_xyz = stencilInternalOffset(a_index, b_index, c_index, invu);
  incorporateStencilOrigin(orig_x, orig_y, orig_z, pt_xyz, scale_factor, point_x, point_y,
                           point_z);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void fvStencilCoordinates(const T orig_x, const T orig_y, const T orig_z,
                          const UnitCellAxis face_normal, const int a_index, const int b_index,
                          const int c_index, const double* invu, T *point_x, T *point_y,
                          T *point_z) {
  fvStencilCoordinates(orig_x, orig_y, orig_z, 1.0, face_normal, a_index, b_index, c_index,
                             invu, point_x, point_y, point_z);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void fvStencilCoordinates(const T orig_x, const T orig_y, const T orig_z, const int a_index,
                          const int b_index, const int c_index, const double* invu, T *point_x,
                          T *point_y, T *point_z) {
  fvStencilCoordinates(orig_x, orig_y, orig_z, 1.0, a_index, b_index, c_index, invu, point_x,
                             point_y, point_z);
}

} // namespace stmath
} // namespace stormm
