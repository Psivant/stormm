// -*-c++-*-
#ifndef STORMM_HPC_BACKGROUND_MESH_H
#define STORMM_HPC_BACKGROUND_MESH_H

#ifdef STORMM_USE_CUDA
#  include <cuda_runtime.h>
#endif
#include "copyright.h"
#include "Constants/behavior.h"
#include "Math/math_enumerators.h"
#include "Topology/atomgraph_enumerators.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

using constants::PrecisionModel;
using stmath::Interpolant;
using topology::UnitCellType;
  
#ifdef STORMM_USE_CUDA
/// \brief Get the thread and memory requirements for various kernels managing information flow
///        between particle positions and mesh-based fields.
///
/// Overloaded:
///   - Provide precision model, process, mesh type, and boundary conditions, the minimum
///     descriptors for occlusion meshes
///   - Provide unit cell and interpolant specifications for nonbonded field meshes
///
/// \param prec          Numerical precision model for real-valued numbers
/// \param picture       The type of field described by the mesh
/// \param bounds        Boundary conditions for the mesh mapper to use
/// \param unit_cell     The unit cell type describing the mesh element
/// \param stencil_kind  The plan under which the coefficients of each mesh element's interpolant
///                      will be fitted
/// \param activity      Indicate whether particle information is being placed on the mesh, or
///                      mesh-based forces and energies are being placed back on particles.
/// \{
cudaFuncAttributes queryGridKernelRequirements(PrecisionModel prec, GridDetail picture,
                                               BoundaryCondition bounds, MappingActivity process);

cudaFuncAttributes queryGridKernelRequirements(PrecisionModel prec, GridDetail picture,
                                               BoundaryCondition bounds, UnitCellType unit_cell,
                                               Interpolant stencil_kind, MappingActivity process);
/// \}
#endif

} // namespace structure
} // namespace stormm

#endif
