// -*-c++-*-
#ifndef STORMM_MESH_COMPONENTS_H
#define STORMM_MESH_COMPONENTS_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/stormm_vector_types.h"
#include "Numerics/split_fixed_precision.h"
#include "local_arrangement.h"
#include "mesh_parameters.h"

namespace stormm {
namespace structure {

using card::Hybrid;
using card::HybridTargetLevel;
#ifndef STORMM_USE_HPC
using data_types::double3;
#endif

struct MeshRulerKit {

  /// \brief As with all abstracts, the constructor takes a straight list of inputs corresponding
  ///        to member variables in the source class, covering all member variables.  For
  ///        optimization purposes, this abstract does not contain lengths of each array needed
  ///        for safe traversal.  That information can be found in the associated MeshParameters
  ///        object which guided the creation of the axes.
  MeshRulerKit(const llint* avec_x_in, const llint* avec_y_in, const llint* avec_z_in,
               const llint* bvec_x_in, const llint* bvec_y_in, const llint* bvec_z_in,
               const llint* cvec_x_in, const llint* cvec_y_in, const llint* cvec_z_in,
               const llint* avec_abs_x_in, const llint* avec_abs_y_in, const llint* avec_abs_z_in,
               const int* avec_x_ovrf_in, const int* avec_y_ovrf_in, const int* avec_z_ovrf_in,
               const int* bvec_x_ovrf_in, const int* bvec_y_ovrf_in, const int* bvec_z_ovrf_in,
               const int* cvec_x_ovrf_in, const int* cvec_y_ovrf_in, const int* cvec_z_ovrf_in,
               const int* avec_abs_x_ovrf_in, const int* avec_abs_y_ovrf_in,
               const int* avec_abs_z_ovrf_in);

  // As with other abstracts, the default copy and move constructors are acceptable, but the copy
  // and move assignment operators are implicitly deleted because of const member variables.
  /// \{
  MeshRulerKit(const MeshRulerKit &original) = default;
  MeshRulerKit(MeshRulerKit &&original) = default;
  /// \}
  
  // Public, const member variables for the coordinates of tick marks along each axis
  const llint* avec_x;         ///< Relative Cartesian X coordinates of a axis tick marks
  const llint* avec_y;         ///< Relative Cartesian Y coordinates of a axis tick marks
  const llint* avec_z;         ///< Relative Cartesian Z coordinates of a axis tick marks
  const llint* bvec_x;         ///< Relative Cartesian X coordinates of b axis tick marks
  const llint* bvec_y;         ///< Relative Cartesian Y coordinates of b axis tick marks
  const llint* bvec_z;         ///< Relative Cartesian Z coordinates of b axis tick marks
  const llint* cvec_x;         ///< Relative Cartesian X coordinates of c axis tick marks
  const llint* cvec_y;         ///< Relative Cartesian Y coordinates of c axis tick marks
  const llint* cvec_z;         ///< Relative Cartesian Z coordinates of c axis tick marks
  const llint* avec_abs_x;     ///< Absolute Cartesian X coordinates of a axis tick marks
  const llint* avec_abs_y;     ///< Absolute Cartesian Y coordinates of a axis tick marks
  const llint* avec_abs_z;     ///< Absolute Cartesian Z coordinates of a axis tick marks
  const int* avec_x_ovrf;      ///< Overflow bits for a axis relative Cartesian X coordinates
  const int* avec_y_ovrf;      ///< Overflow bits for a axis relative Cartesian Y coordinates
  const int* avec_z_ovrf;      ///< Overflow bits for a axis relative Cartesian Z coordinates
  const int* bvec_x_ovrf;      ///< Overflow bits for b axis relative Cartesian X coordinates
  const int* bvec_y_ovrf;      ///< Overflow bits for b axis relative Cartesian Y coordinates
  const int* bvec_z_ovrf;      ///< Overflow bits for b axis relative Cartesian Z coordinates
  const int* cvec_x_ovrf;      ///< Overflow bits for c axis relative Cartesian X coordinates
  const int* cvec_y_ovrf;      ///< Overflow bits for c axis relative Cartesian Y coordinates
  const int* cvec_z_ovrf;      ///< Overflow bits for c axis relative Cartesian Z coordinates
  const int* avec_abs_x_ovrf;  ///< Overflow bits for Cartesian X absolute a axis coordinates
  const int* avec_abs_y_ovrf;  ///< Overflow bits for Cartesian Y absolute a axis coordinates
  const int* avec_abs_z_ovrf;  ///< Overflow bits for Cartesian Z absolute a axis coordinates
};

/// \brief A collection of coordinate vectors describing each axis of a mesh.  If the mesh axes
///        track Cartesian axes, then all member variables but a_line_x, b_line_y, and c_line_z
///        plus their overflow arrays will be zero.
class MeshRulers {
public:
  
  /// \brief The constructor takes a full collection of mesh parameters.
  MeshRulers(const MeshParameters &mps = MeshParameters());

  /// \brief The copy and move constructors must be explicitly written to repair POINTER-kind
  ///        Hybrid objects.
  ///
  /// \param original  An existing object to copy or move
  /// \param other     An existing object placed on the right hand side of an assignment statement
  /// \{
  MeshRulers(const MeshRulers &original);
  MeshRulers(MeshRulers &&original);
  MeshRulers& operator=(const MeshRulers &original);
  MeshRulers& operator=(MeshRulers &&original);
  /// \}

  /// \brief Get the origin from the rulers themselves
  double3 getMeshOrigin() const;

  /// \brief Get the location of a point in Cartesian space based on its coordinates on the mesh.
  ///        This function is not intended for high-volume computations, merely to provide unit
  ///        tests of the object's validity in a given context.
  ///
  /// \param mesh_loc  The location cast in units of mesh elements
  double3 getRealLocation(const double3 mesh_loc) const;

  /// \brief Get the location of a point on the mesh based on its coordinates in Cartesian space.
  ///        This function is not intended for high-volume computations, merely to provide unit
  ///        tests of the object's validity in a given context.
  ///
  /// \param real_loc  The location cast in units of mesh elements
  double3 getMeshLocation(const double3 real_loc) const;

  /// \brief Get a collection of all relevant pointers for memory on either the CPU host or the
  ///        GPU device.  The data will be read-only.
  ///
  /// \param tier  Specify wheter to obtain pointers on the CPU host or GPU device
  const MeshRulerKit data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

#ifdef STORMM_USE_HPC
  /// \brief Upload all data to the GPU.
  void upload();

  /// \brief Download all data from the GPU.
  void download();

#  ifdef STORMM_USE_CUDA
  /// \brief Get a collection of all relevant pointers for memory on the host, visible from the
  ///        GPU device.  The data will be read-only.
  const MeshRulerKit deviceViewToHostData() const;
#  endif
#endif
  
private:
  
  /// Fixed-precision Cartesian coordinates of the mesh grid points stepping along the lines
  /// [ i = 0...nx, 0, 0 ] (the "a" box vector), [ 0, j = 0...ny, 0 ] (the "b" box vector),
  /// [ 0, 0, k = 0...nz ] (the "c" box vector), and a fourth vector with absolute coordinates of
  /// the "a" bo vector points, [ i = 0...nx, 0, 0 ] + [ orig_x, orig_y, orig_z ].  Storing these
  /// values obviates the need to do expensive and complicated multiplications of fixed-precision
  /// numbers during grid-based particle calculations.
  /// \{
  Hybrid<llint> a_line_x;      //
  Hybrid<llint> a_line_y;      // Positions of "a" vector tick marks relative to the mesh origin
  Hybrid<llint> a_line_z;      //
  Hybrid<llint> b_line_x;      //
  Hybrid<llint> b_line_y;      // Positions of "b" vector tick marks relative to the mesh origin
  Hybrid<llint> b_line_z;      //
  Hybrid<llint> c_line_x;      //
  Hybrid<llint> c_line_y;      // Positions of "c" vector tick marks relative to the mesh origin
  Hybrid<llint> c_line_z;      //
  Hybrid<llint> a_abs_line_x;  //
  Hybrid<llint> a_abs_line_y;  // Absolute positions of the "a" vector tick marks
  Hybrid<llint> a_abs_line_z;  //
  /// \}

  /// Overflow for fixed-precision Cartesian coordinates of the mesh grid points stepping along
  /// the a, b, and c box vectors.
  /// \{
  Hybrid<int> a_line_x_overflow;      // Overflow bits for positions of the "a" vector tick marks
  Hybrid<int> a_line_y_overflow;      //   relative to the mesh origin.  These engage when the mesh
  Hybrid<int> a_line_z_overflow;      //   serves double-precision coefficients and calculatons.
  Hybrid<int> b_line_x_overflow;      //
  Hybrid<int> b_line_y_overflow;      // Overflow bits for "b" vector tick mark relative positions
  Hybrid<int> b_line_z_overflow;      //
  Hybrid<int> c_line_x_overflow;      //
  Hybrid<int> c_line_y_overflow;      // Overflow bits for "c" vector tick mark relative positions
  Hybrid<int> c_line_z_overflow;      //
  Hybrid<int> a_abs_line_x_overflow;  //
  Hybrid<int> a_abs_line_y_overflow;  // Overflow bits for "a" vector tick mark absolute positions
  Hybrid<int> a_abs_line_z_overflow;  //
  /// \}

  /// Data storage for the POINTER-kind Hybrid<int> objects above
  Hybrid<int> int_data;

  /// Data storage for the POINTER-kind Hybrid<llint> objects above
  Hybrid<llint> llint_data;

  /// Store a pointer to the original MeshParameters object that created these rulers.  That
  /// object, and its abstract, contain critical scaling constants for interpreting the
  /// fixed-precision data.
  MeshParameters *mps_pointer;
  
  /// Repair POINTER-kind Hybrid objects after copy and copy assignment operations.
  void rebasePointers();
};

} // namespace structure
} // namespace stormm

#endif
