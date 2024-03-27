// -*-c++-*-
#ifndef STORMM_TRICUBIC_CELL_H
#define STORMM_TRICUBIC_CELL_H

#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "math_enumerators.h"
#include "matrix_ops.h"

namespace stormm {
namespace stmath {

using card::Hybrid;
using card::HybridTargetLevel;
using constants::CartesianDimension;
using constants::UnitCellAxis;
using data_types::isFloatingPointScalarType;

/// \brief A bundle of a 64 x 64 matrix, encoding the transformation of function values and
///        derivatives at a stencil of points inside a tricubic mesh element into coefficients for
///        that mesh element expressing the interpolated tricubic function.  The type of stencil
///        is also indicated, to keep this information coupled.
class TricubicStencil {
public:

  /// \brief The constructor accepts the stencil type.  Dimensions of the matrix and the geometry
  ///        of each stencil option are hard-wired.
  TricubicStencil(Interpolant kind_in = Interpolant::SMOOTHNESS);

  /// \brief The copy and move constructors as well as assignment operators can all take their
  ///        default forms.  There are no const members or pointers to repair.
  ///
  /// \param original  The original object to copy or move
  /// \param other     A pre-existing object to copy or move
  /// \{
  TricubicStencil(const TricubicStencil &original) = default;
  TricubicStencil(TricubicStencil &&original) = default;
  TricubicStencil& operator=(const TricubicStencil &original) = default;
  TricubicStencil& operator=(TricubicStencil &&original) = default;
  /// \}

  /// \brief Get the type of interpolant.
  Interpolant getKind() const;
  
  /// \brief Get a pointer to the matrix data on the host or device.
  ///
  /// \param tier  Specify whether to obtain data on the CPU host or GPU device
  const double* data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  /// \brief Get a pointer to the matrix data on the host visible from the device.
  const double* getDeviceViewToHostData() const;
#  endif

  /// Upload the matrix data to the device.
  void upload();

  /// Download the matrix data from the device.
  void download();
#endif

  /// \brief Export the matrix as a Standard-Template Library vector.
  std::vector<double> exportMatrix() const;
  
private:
  Interpolant kind;          ///< The type of stencil encoded in the transformation
  Hybrid<double> transform;  ///< Transformation matrix to take function and derivative values into
                             ///<   coefficients of the tricubic interpolant
};

/// \brief One element of a three-dimensional, space-filling tile representing a three-dimensional
///        function to be interpolated by a tricubic spline.
template <typename T> class TricubicCell {
public:

  /// \brief The constructor can take nothing and simply initialize all values to zero, or accept
  ///        the tricubic weights matrix, a dimensions array, and details of the potential
  ///        function.
  ///
  /// \param weights_matrix  The inverse matrix of polynomial weights obtained from
  ///                        getTricubicMatrix() in this library
  /// \param bounds          Array containing the Cartesian x, y, and z coordinates of the grid
  ///                        cell origin plus its lengths along each axis (the grid cell is NOT
  ///                        required to be orthorhombic).  This array can have four entries (in
  ///                        which case the final entry is assumed to be the isotropic length
  ///                        parameter), six (for anisotropic cells), or twelve entries, in which
  ///                        case the last nine represent the column matrix of cell vectors.
  /// \{
  TricubicCell();

  TricubicCell(const TricubicStencil weights_matrix, const std::vector<double> &bounds,
               const std::vector<T> &f_in, const std::vector<T> &dx_in,
               const std::vector<T> &dy_in, const std::vector<T> &dz_in,
               const std::vector<T> &dxx_in, const std::vector<T> &dxy_in,
               const std::vector<T> &dxz_in, const std::vector<T> &dyy_in,
               const std::vector<T> &dyz_in, const std::vector<T> &dxxx_in,
               const std::vector<T> &dxxy_in, const std::vector<T> &dxxz_in,
               const std::vector<T> &dxyy_in, const std::vector<T> &dxyz_in);

  TricubicCell(const TricubicStencil weights_matrix, const std::vector<double> &bounds,
               const std::vector<T> &f_in, const std::vector<T> &dx_in,
               const std::vector<T> &dy_in, const std::vector<T> &dz_in,
               const std::vector<T> &dxy_in, const std::vector<T> &dxz_in,
               const std::vector<T> &dyz_in, const std::vector<T> &dxyz_in);
  /// \}

  /// \brief Retrieve one of the 64 coefficients Aijk for the tricubic spline.
  ///
  /// \param i  The ith coefficient relevant to the notation Aijk
  /// \param j  The jth coefficient relevant to the notation Aijk
  /// \param k  The kth coefficient relevant to the notation Aijk
  T getCoefficient(int i, int j, int k) const;

  /// \brief Get a Standard Template Library vector containing all coefficients in the element.
  ///        The order of coefficients is given in Fortran order: for the term Aijk x^i y^j z^k,
  ///        the ith coefficient varies the most rapidly and the index for Aijk is given by
  ///        (((k * 4) + j) * 4 + i.
  std::vector<T> getCoefficients() const;
  
  /// \brief Set one of the 64 coefficients Aijk for the tricubic spline.  Parameter descriptions
  ///        follow from above, with the addition of:
  ///
  /// \param value  The value to apply
  void setCoefficient(T value, int i, int j, int k);

  /// \brief Get one of the data points from the boundary.  Parameter descriptions follow from
  ///        above, with the addition of:
  ///
  /// \param kind  The classification of the boundary condition as a derivative (or value)
  T getData(FunctionLevel kind, int i, int j, int k) const;

  /// \brief Set one of the data items.  Parameter descriptions follow from above.
  void setData(T value, FunctionLevel kind, int i, int j, int k);

  /// \brief Get the cell origin along one dimension.
  ///
  /// \param dim  The Cartesian dimension along which to return the origin coordinate
  T getCellOrigin(CartesianDimension dim) const;

  /// \brief Get the cell length along one dimension.
  ///
  /// \param dim  The Cartesian dimension along which to return the cell length
  T getCellLength(CartesianDimension dim) const;

  /// \brief Get the fractional position of a point in Cartesian coordinates within a mesh cell.
  ///        This includes a check on whether the point is actually within the mesh cell.
  ///
  /// \param x       Cartesian X coordinate
  /// \param y       Cartesian Y coordinate
  /// \param z       Cartesian Z coordinate
  /// \param a_frac  Fractional coordinate along the mesh element's a axis, evaluated and returned
  /// \param b_frac  Fractional coordinate along the mesh element's b axis, evaluated and returned
  /// \param c_frac  Fractional coordinate along the mesh element's c axis, evaluated and returned
  /// \param policy  Indicate the course of action if the point lies outside the mesh element
  /// \param caller  Name of the calling function
  void fractionalPosition(T x, T y, T z, T *a_frac, T *b_frac, T *c_frac,
                          ExceptionResponse policy = ExceptionResponse::SILENT,
                          const char* caller = nullptr) const;

  /// \brief Evaluate the function at a specific point in space.  This will take into account the
  ///        grid cell's origin and lengths to determine where in the grid cell the point of
  ///        interest lies.  If the point is outside the grid cell, produce an error.
  ///
  /// \param x     Cartesian X location of the point
  /// \param y     Cartesian Y location of the point
  /// \param z     Cartesian Z location of the point
  T evaluate(T x, T y, T z) const;

  /// \brief Evaluate the first derivatives of a tricubic mesh element at a point and return the
  ///        results along Cartesian axes.  Descriptions of input parameters follow from the
  ///        evaluate() member function above.
  template <typename T3> T3 derivative(T x, T y, T z) const;

  /// \brief Evaluate the second derivatives of a tricubic mesh element at a point and return the
  ///        results along Cartesian axes.  Descriptions of input parameters follow from the
  ///        evaluate() member function above, with the addition of:
  ///
  /// Overloaded:
  ///   - Use a pre-allocated space to hold the results (provided as a C-style array or Standard
  ///     Template Library vector)
  ///   - Produce a Standard Template Library vector holding the results (of trusted length 9)
  ///   
  /// \param result  Pre-allocated space to hold the output
  /// \{
  void secondDerivative(T x, T y, T z, T* result) const;
  void secondDerivative(T x, T y, T z, std::vector<T> *result) const;
  std::vector<T> secondDerivative(T x, T y, T z) const;
  /// \}
  
  /// \brief Evaluate the second derivatives of a tricubic mesh element at a point and return the
  ///        results along Cartesian axes.  Descriptions of input parameters follow from the
  ///        evaluate() member function above, with the addition of:
  ///
  /// Overloaded:
  ///   - Use a pre-allocated space to hold the results (provided as a C-style array or Standard
  ///     Template Library vector)
  ///   - Produce a Standard Template Library vector holding the results (of trusted length 27)
  ///   
  /// \param result  Pre-allocated space to hold the output tensor
  /// \{
  void thirdDerivative(T x, T y, T z, T* result) const;
  void thirdDerivative(T x, T y, T z, std::vector<T> *result) const;
  std::vector<T> thirdDerivative(T x, T y, T z) const;
  /// \}

private:
  T coefficients[64];  ///< Solved coefficients of the tricubic spline that satisfies all boundary
                       ///<   conditions.

  // The following arrays store their series of values in "Fortran" order: (X0, Y0, Z0),
  // (X1, Y0, Z0), (X0, Y1, Z0), (X1, Y1, Z0), (X0, Y0, Z1), ...
  T f[8];      ///< Values of the function at the bounding grid points
  T dx[8];     ///< Cartesian X derivatives of the function at the bounding grid points
  T dy[8];     ///< Cartesian Y derivatives of the function at the bounding grid points
  T dz[8];     ///< Cartesian Z derivatives of the function at the bounding grid points
  T dxx[8];    ///< Cartesian X second derivatives of the function at the bounding grid points
               ///<   (not required for orthorhombic meshes)
  T dxy[8];    ///< Cartesian X/Y cross-derivatives of the function at the bounding grid points
  T dxz[8];    ///< Cartesian X/Z cross-derivatives of the function at the bounding grid points
  T dyy[8];    ///< Cartesian Y second derivatives of the function at the bounding grid points
               ///<   (not required for orthorhombic meshes)
  T dyz[8];    ///< Cartesian Y/Z cross-derivatives of the function at the bounding grid points
  T dxxx[8];   ///< Cartesian X/X/X triple-derivatives of the function at the bounding grid points
               ///<   (not required for orthorhombic meshes)
  T dxxy[8];   ///< Cartesian X/X/Y triple-derivatives of the function at the bounding grid points
               ///<   (not required for orthorhombic meshes)
  T dxxz[8];   ///< Cartesian X/X/Z triple-derivatives of the function at the bounding grid points
               ///<   (not required for orthorhombic meshes)
  T dxyy[8];   ///< Cartesian X/X/Z triple-derivatives of the function at the bounding grid points
               ///<   (not required for orthorhombic meshes)
  T dxyz[8];   ///< Cartesian X/Y/Z triple-derivatives of the function at the bounding grid points

  // The grid cell boundaries are stored in double precision for accuracy considerations.
  double origin_x;  ///< Cartesian X location of the grid cell origin
  double origin_y;  ///< Cartesian Y location of the grid cell origin
  double origin_z;  ///< Cartesian Z location of the grid cell origin

  /// Transformation matrix to take Cartesian coordinates, as displacements from the cell origin,
  /// into the fractional space of the cell over which the tricubic spline is applicable.
  double umat[9];
  
  /// Column matrix (given in Fortran order, first column being entries 0, 1, and 2, second column
  /// 3, 4, and 5, ...) of the cell's bounding vectors.  For rectilinear (orthorhombic) cells, the
  /// Cartesian X, Y, and Z lengths are given in entries 0, 4, and 8.
  double invu[9];
};

/// \brief Get the internal offset for the point of interest within the mesh element, as expressed
///        in terms of Cartesian displacements.
///
/// Overloaded:
///   - Provide a unit cell axis to evaluate one of the octets of outer points
///   - Provide no axis to evaluate the inner octet of points near the mesh element center
///
/// \param face_normal  Specify the unit cell axis normal to the face of interest.  Specify no face
///                     for the octet of points in the center of the mesh element.
/// \param a_index      Index of the point along the unit cell A axis
/// \param b_index      Index of the point along the unit cell B axis
/// \param c_index      Index of the point along the unit cell C axis
/// \param invu         Transformation matrix for taking coordinates in unit cell space into
///                     Cartesian space (must be provided in double precision, assumed to hold
///                     nine elements)
/// \{
double3 stencilInternalOffset(UnitCellAxis face_normal, int a_index, int b_index, int c_index,
                              const double* invu);

double3 stencilInternalOffset(int a_index, int b_index, int c_index, const double* invu);
/// \}

/// \brief Combine the Cartesian location of the origin with the internal offset of a point of
///        interest to arrive at a Cartesian representation of the point in an absolute frame of
///        reference.
///
/// Overloaded:
///   - Accept the origin coordinates as int95_t split fixed-precision numbers
///   - Accept the origin coordinates in scalar data types, e.g. those compatible with the
///     CoordinateSeries object
///
/// \param orig_x        Cartesian X position of the mesh element origin
/// \param orig_y        Cartesian Y position of the mesh element origin
/// \param orig_z        Cartesian Z position of the mesh element origin
/// \param pt_xyz        Relative Cartesian location of the point within the element 
/// \param scale_factor  Scaling factor applied to fixed-precision coordinates of the origin, and
///                      to be applied to the point's internal displacement within the element
/// \param point_x       Cartesian X location of the point, evaluated and returned
/// \param point_y       Cartesian Y location of the point, evaluated and returned
/// \param point_z       Cartesian Z location of the point, evaluated and returned
/// \{
void incorporateStencilOrigin(const int95_t orig_x, const int95_t orig_y, const int95_t orig_z,
                              const double3 pt_xyz, double scale_factor, int95_t *point_x,
                              int95_t *point_y, int95_t *point_z);

template <typename T>
void incorporateStencilOrigin(T orig_x, T orig_y, T orig_z, const double3 pt_xyz,
                              double scale_factor, T *point_x, T *point_y, T *point_z);
/// \}

/// \brief Compute the Cartesian coordinates of mesh points involved in a tricubic stencil tailored
///        to maximize accuracy in the result.  This abstract CPU-based code but it otherwise not
///        efficient for high-intensity computations.
///
/// Overloaded:
///   - Provide a unit cell axis to evaluate one of the octets of outer points
///   - Provide no axis to evaluate the inner octet of points near the mesh element center
///   - Provide a scaling factor to indicate fixed-precision integer coordinate representations, or
///     let it be assumed that the scaling factor is one (coordinates must be floating-point type)
///
/// Descriptions of parameters follow from stencilInternalOffset() and incorporateStencilOrigin()
/// above.
/// \{
template <typename T>
void fvStencilCoordinates(T orig_x, T orig_y, T orig_z, double scale_factor,
                          UnitCellAxis face_normal, int a_index, int b_index, int c_index,
                          const double* invu, T *point_x, T *point_y, T *point_z);

template <typename T>
void fvStencilCoordinates(T orig_x, T orig_y, T orig_z, double scale_factor, int a_index,
                          int b_index, int c_index, const double* invu, T *point_x,
                          T *point_y, T *point_z);

template <typename T>
void fvStencilCoordinates(T orig_x, T orig_y, T orig_z, UnitCellAxis face_normal, int a_index,
                          int b_index, int c_index, const double* invu, T *point_x, T *point_y,
                          T *point_z);

template <typename T>
void fvStencilCoordinates(T orig_x, T orig_y, T orig_z, int a_index, int b_index, int c_index,
                          const double* invu, T *point_x, T *point_y, T *point_z);
/// \}

} // namespace stmath
} // namespace stormm

#include "tricubic_cell.tpp"

#endif
