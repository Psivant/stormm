// -*-c++-*-
#ifndef STORMM_MESH_KERNEL_MANAGER_H
#define STORMM_MESH_KERNEL_MANAGER_H

#include <map>
#include <string>
#include "copyright.h"
#include "Constants/behavior.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/math_enumerators.h"
#include "Structure/structure_enumerators.h"
#include "Topology/atomgraph_enumerators.h"
#include "gpu_details.h"
#include "kernel_format.h"

namespace stormm {
namespace card {

using constants::PrecisionModel;
using stmath::Interpolant;
using structure::BoundaryCondition;
using structure::GridDetail;
using structure::MappingActivity;
using topology::UnitCellType;
  
/// \brief Collect and dispense the launch parameters for various mesh-based kernels.  The launches
///        for these kernels will generally take place in the context of a BackgroundMesh object,
///        whether during its construction or for its application.
class MeshKlManager : public KernelManager {
public:

  /// \brief As the construction of virtually any mesh will provide enough operations for a GPU to
  ///        run at acceptable efficiency.  The constructor will therefore be based on the GPU
  ///        specs as if the workload is for a single mesh.
  MeshKlManager(const GpuDetails &gpu_in = null_gpu);
  
  /// \brief Get the architecture-specific block multiplier.  This will run a minimum number of
  ///        blocks per streaming multiprocessor on some cards, specifically NVIDIA's GTX 1080-Ti,
  ///        when large blocks cannot use more than 32k registers in all.
  int getArchBlockMultiplier() const;

  /// \brief Get the launch parameters for a mesh computation kernel.
  ///
  /// Overloaded:
  ///   - Provide precision model, process, mesh type, and boundary conditions, the minimum
  ///     descriptors for occlusion meshes
  ///   - Provide unit cell and interpolant specifications for nonbonded field meshes
  ///
  /// \param prec          The precision model in which to perform calculations and present results
  /// \param picture       The type of potential or field presented by the mesh
  /// \param bounds        Boundary conditions for the mesh mapper of interest
  /// \param unit_cell     The unit cell type of each mesh element (and, by extension, of the mesh
  ///                      as a whole)
  /// \param stencil_kind  The framework for fitting tricubic polynomial coefficients in each mesh
  ///                      element
  /// \param process       The direction of information flow between the mesh and particles
  /// \{
  int2 getMeshKernelDims(PrecisionModel prec, GridDetail picture, BoundaryCondition bounds,
                         MappingActivity process) const;

  int2 getMeshKernelDims(PrecisionModel prec, GridDetail picture, BoundaryCondition bounds,
                         UnitCellType unit_cell, Interpolant stencil_kind,
                         MappingActivity process) const;
  /// \}

private:

  /// The architecture-specific block multipliers for mesh-based computation kernels.  The thread
  /// count for each kernel is set to 256 (small).  Mapping kernels (Particle to Mesh, ptom) tend
  /// to operate as collections of threads, whereas most interpolation kernels (Mesh to Particle,
  /// mtop) operate as a collection of warps.
  /// \{
  int occ_mesh_block_multiplier_dp;
  int occ_mesh_block_multiplier_sp;
  int nbf_mesh_block_multiplier_dp;
  int nbf_mesh_block_multiplier_sp;
  /// \}

  /// \brief Set the maximum block size and block multipliers for various mesh-related kernels.
  ///
  /// Overloaded:
  ///   - Provide precision model, process, mesh type, and boundary conditions, the minimum
  ///     descriptors for occlusion meshes
  ///   - Provide unit cell and interpolant specifications for nonbonded field meshes
  ///
  /// \param prec          The type of floating point numbers in which the kernel shall work
  /// \param process       Indicate whether particles are being mapped to the mesh or information
  ///                      from the mesh is being mapped back to particles
  /// \param picture       Type of mesh that is being filled or interpolated
  /// \param bounds        Boundary conditions for the mesh mapper of interest
  /// \param unit_cell     The unit cell type of each mesh element (and, by extension, of the mesh
  ///                      as a whole)
  /// \param stencil_kind  The framework for fitting tricubic polynomial coefficients in each mesh
  ///                      element
  /// \param kernel_name   [Optional] Name of the kernel in the actual code
  /// \{
  void catalogMeshKernel(PrecisionModel prec, MappingActivity process, GridDetail picture,
                         BoundaryCondition bounds,
                         const std::string &kernel_name = std::string(""));

  void catalogMeshKernel(PrecisionModel prec, MappingActivity process, GridDetail picture,
                         BoundaryCondition bounds, UnitCellType unit_cell,
                         Interpolant stencil_kind,
                         const std::string &kernel_name = std::string(""));
  /// \}
};

/// \brief Obtain the block multiplier for mesh mapping kernels.  Each block will run on 128 to 256
///        threads and most are memory-bound.
///
/// \param prec     The type of floating point numbers in which the kernel shall work
/// \param picture  Type of mesh that is being filled or interpolated
int meshBlockMultiplier(PrecisionModel prec, GridDetail picture);

/// \brief Obtain a unique string identifier for one of the mesh mapping kernels.  Each identifier
///        begins with "mesh_" and is then appended with letter codes as follows:
///        - { d, f }            Perform calculations in double (d) or float (f) arithmetic
///        - { oc, of, nf, na }  The mesh comprises a stepwise occlusion mask (oc), a
///                              differentiable occlusion field (of), a non-bonded field (nf), or a
///                              non-bonded field with rigid atom neighbor lists (n)
///        - { iso, pbc }        The mesh uses isolated or periodic boundary conditions
///        - { rc, tr }          The mesh element is rectilinear (orthorhombic) or triclinic
///        - { sm, va }          Coefficients for the tricubic interpolant are computed to optimize
///                              smoothness in derivatives (sm) or function values throughout the
///                              element (va)
///        - { pm, mp }          Map particle information to the mesh (pm), or mesh-based effects
///                              to particles (mp)
///
/// Overloaded:
///   - Provide precision model, process, mesh type, and boundary conditions, the minimum
///     descriptors for occlusion meshes
///   - Provide unit cell and interpolant specifications for nonbonded field meshes
///
/// \param prec          The precision model in which to compute the mesh field
/// \param picture       The type of non-bonded potential mapped ot the mesh
/// \param bounds        Boundary conditions for the mesh mapper of interest
/// \param unit_cell     The unit cell type of each mesh element (and, by extension, of the mesh
///                      as a whole)
/// \param stencil_kind  The framework for fitting tricubic polynomial coefficients in each mesh
///                      element
/// \param process       Indicate whether particles are being mapped to the mesh or information
///                      from the mesh is being mapped back to particles
/// \{
std::string meshKernelKey(PrecisionModel prec, GridDetail picture, BoundaryCondition bounds,
                          MappingActivity process);

std::string meshKernelKey(PrecisionModel prec, GridDetail picture, BoundaryCondition bounds,
                          UnitCellType unit_cell, Interpolant stencil_kind,
                          MappingActivity process);
/// \}

} // card
} // stormm

#endif
