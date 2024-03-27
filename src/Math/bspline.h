// -*-c++-*-
#ifndef STORMM_BSPLINE_H
#define STORMM_BSPLINE_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/common_types.h"

namespace stormm {
namespace stmath {

using card::Hybrid;
using data_types::isSignedIntegralScalarType;

/// \brief Check the lengths of various arrays needed by the more complex B-spline routines
///        serving atomic coordinates.  This is performed when the input arrays each have their
///        own lengths which can be checked.
///
/// \param xcrd       Cartesian X coordinates of particles
/// \param ycrd       Cartesian Y coordinates of particles
/// \param zcrd       Cartesian Z coordinates of particles
/// \param a_coefs    Array of coefficients for B-splines along the unit cell A axis
/// \param b_coefs    Array of coefficients for B-splines along the unit cell B axis
/// \param c_coefs    Array of coefficients for B-splines along the unit cell C axis
/// \param a_init     Initial seeding point for knots along the unit cell A axis
/// \param b_init     Initial seeding point for knots along the unit cell B axis
/// \param c_init     Initial seeding point for knots along the unit cell C axis
template <typename Tcoord, typename Tcalc>
void bSplineInputChecks(const std::vector<Tcoord> &xcrd, const std::vector<Tcoord> &ycrd,
			const std::vector<Tcoord> &zcrd, const int order,
			std::vector<Tcalc> *a_coefs, std::vector<Tcalc> *b_coefs,
			std::vector<Tcalc> *c_coefs, std::vector<int> *a_init,
			std::vector<int> *b_init, std::vector<int> *c_init);
  
/// \brief Compute B-spline cofficients for a given set of coordinates.  The B-spline coefficients
///        are returned in reverse order, due to conveniences in the way they can then be computed.
///        When implementing B-splines in the context of particles on a grid, each B-spline knot
///        pertain to the grid point corresponding to the index of the modified array: if the
///        coordinate is entered as x in the range [0, 1), the first knot will influence the grid
///        point corresponding to grid point n-1, the next grid point n-2, then n-3, ..., 1, 0 for
///        an nth order spline.
///
/// Overloaded:
///   - Compute for a specific coordinate on the grid.  Return a Standard Template Library object
///     with the computed coefficients.
///   - Compute for a specific coordinate on the grid.  Fill a pre-allocated array of coefficients
///     (this is the fastest form for B-splines in a signle dimension).
///   - Compute for three series of Cartesian coordinates with a given unit cell and mesh
///     discretization.  Modify and return pre-allocated arrays of coefficients and index starting
///     points on the mesh in terms of its A, B, and C axes.
///
/// \param x          Coordinate of the particle or point, in the range [0, 1)
/// \param xcrd       Cartesian X coordinates of particles
/// \param ycrd       Cartesian Y coordinates of particles
/// \param zcrd       Cartesian Z coordinates of particles
/// \param natom      The number of particles to map.  This is the trusted length of xcrd, ycrd,
///                   zcrd, a_init, b_init, and c_init.  Furthermore, when multiplied by order it
///                   implies a trusted length for a_coefs, b_coefs, and c_coefs.
/// \param order      Order of the B-spline to compute
/// \param umat_cell  Unit cell (or, spatial decomposition cell) transformation matrix to put
///                   particles onto the mesh
/// \param invu_mesh  Inverse transformation matrix whose columns define an element of the mesh
///                   onto which the B-splines will be mapped.  This will be some factor, e.g. 4,
///                   a quarter of the size, of the inverse of umat_cell.  
/// \param coefs      Array of coefficients (modified and returned)
/// \param workspace  Pre-allocated workspace for computing B-spline coefficients, using the
///                   recursive relationship but without using explicit recursion
/// \param a_coefs    Array of coefficients for B-splines along the unit cell A axis
/// \param b_coefs    Array of coefficients for B-splines along the unit cell B axis
/// \param c_coefs    Array of coefficients for B-splines along the unit cell C axis
/// \param a_init     Initial seeding point for knots along the unit cell A axis
/// \param b_init     Initial seeding point for knots along the unit cell B axis
/// \param c_init     Initial seeding point for knots along the unit cell C axis
/// \param da_coefs   Array of coefficients for B-spline derivatives along the unit cell A axis
/// \param db_coefs   Array of coefficients for B-spline derivatives along the unit cell B axis
/// \param dc_coefs   Array of coefficients for B-spline derivatives along the unit cell C axis
/// \{
template <typename T> void bSpline(T x, int order, T* coefs, T* dcoefs = nullptr);

template <typename T> std::vector<T> bSpline(T x, int order);
  
template <typename Tcoord, typename Tcalc>
void bSpline(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, int order,
             const Tcalc* umat_cell, const Tcoord* invu_cell, int mesh_na, int mesh_nb,
             int mesh_nc, const Tcoord* invu_mesh, Tcalc* a_coefs, Tcalc* b_coefs, Tcalc* c_coefs,
             int* a_init, int* b_init, int* c_init, Tcalc* da_coefs = nullptr,
             Tcalc* db_coefs = nullptr, Tcalc* dc_coefs = nullptr, Tcalc coordinate_scale = 1.0);

template <typename Tcoord, typename Tcalc>
void bSpline(const std::vector<Tcoord> &xcrd, const std::vector<Tcoord> &ycrd,
             const std::vector<Tcoord> &zcrd, int order, const Tcalc* umat_cell,
             const Tcoord* invu_cell, int mesh_na, int mesh_nb, int mesh_nc,
             const Tcoord* invu_mesh, std::vector<Tcalc> *a_coefs, std::vector<Tcalc> *b_coefs,
             std::vector<Tcalc> *c_coefs, std::vector<int> *a_init, std::vector<int> *b_init,
             std::vector<int> *c_init, std::vector<Tcalc> *da_coefs = nullptr,
             std::vector<Tcalc> *db_coefs = nullptr, std::vector<Tcalc> *dc_coefs = nullptr,
             Tcalc coordinate_scale = 1.0);

template <typename Tcoord, typename Tcalc>
void bSpline(const Hybrid<Tcoord> &xcrd, const Hybrid<Tcoord> &ycrd, const Hybrid<Tcoord> &zcrd,
             int order, const Tcalc* umat_cell, const Tcoord* invu_cell, int mesh_na, int mesh_nb,
             int mesh_nc, const Tcoord* invu_mesh, Hybrid<Tcalc> *a_coefs, Hybrid<Tcalc> *b_coefs,
             Hybrid<Tcalc> *c_coefs, Hybrid<int> *a_init, Hybrid<int> *b_init,
             Hybrid<int> *c_init, Hybrid<Tcalc> *da_coefs = nullptr,
             Hybrid<Tcalc> *db_coefs = nullptr, Hybrid<Tcalc> *dc_coefs = nullptr,
             Tcalc coordinate_scale = 1.0);
/// \}

/// \brief Compute B-splines using the recursion relationship but without applying the logic that
///        the splines are a smooth partition of unity.  The coefficients are again returned in
///        reverse order.
/// 
/// \param x      Coordinate of the particle or point, in the range [0, 1)
/// \param order  Order of the B-spline to compute
template <typename T> std::vector<T> bSplineNoUnity(T x, int order);

/// \brief Return the derivatives of a B-spline, only.  This function works by composing the
///        derivatives of a B-spline of order K by taking differences in the knot values of a
///        B-spline of order K - 1.  This function is far from optimal and intended for testing
///        purposes only.
///
/// \param x                  Coordinate of the particle or point, in the range [0, 1)
/// \param order              Order of the B-spline to compute
/// \param exploit_partition  Flag to make use of the fact that B-splines form a smooth partition
///                           of unity.  In this way, the most complicated B-spline coefficients
///                           can be calculated by subtracting other values from one, saving
///                           computations as well as improving the conservation of density on the
///                           mesh.  This will toggle the exact B-spline computation method.
template <typename T> std::vector<T> dBSpline(T x, int order, bool exploit_partition = true);

} // namespace stmath
} // namespace stormm

#include "bspline.tpp"

#endif
