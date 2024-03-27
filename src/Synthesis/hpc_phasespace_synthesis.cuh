// -*-c++-*-
#ifndef STORMM_HPC_PHASE_SPACE_SYNTHESIS_CUH
#define STORMM_HPC_PHASE_SPACE_SYNTHESIS_CUH

#include <vector>
#ifdef STORMM_USE_CUDA
#include <cuda.h>
#include <cuda_runtime.h>
#endif
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "phasespace_synthesis.h"

namespace stormm {
namespace synthesis {

using constants::large_block_size;

/// \brief Transfer a subset of system coordinates (and perhaps box dimensions) data from specific
///        systems within a PhaseSpaceSynthesis object betwee the host and the accelerator device.
///        PhaseSpaceSynthesis is a central object within STORMM and can be huge (multiple
///        gigabytes worth of coordinate, velocity, and force data.  As such, it deserves a
///        dedicated kernel for managing data transfer to and from the host.
///
/// \param destination  Collection of pointers to PhaseSpaceSynthesis data (this must be visible
///                     on the device, but could be device-only memory or host-mapped memory)
/// \param source       Collection of pointers to PhaseSpaceSynthesis data (this must be visible
///                     on the device, but could be device-only memory or host-mapped memory)
/// \param low_index    Lower bound of systems to upload
/// \param high_index   Upper bound of systems to upload (the range is [low_index, high_index))
__global__ void __launch_bounds__(large_block_size, 1)
kSystemTransfer(PsSynthesisWriter destination, PsSynthesisWriter source, int low_index,
                int high_index, const TrajectoryKind material);

/// \brief Initialize forces for one or all systems of a PhaseSpaceSynthesis on the GPU.
///
/// \param psyw   Writeable abstract for the PhaseSpaceSynthesis, containing system limits and
///               pointers to all coordinates and forces
/// \param index  Index of the system to initialize; if negative, all systems will be initialized.
__global__ void __launch_bounds__(large_block_size, 1)
kPsyInitializeForces(PsSynthesisWriter psyw, const int index);

/// \brief Initialize critical buffers in the phase space (specifically, prior coordinates and
///        velocities) which would otherwise not be used in energy minimization calculations.
///
/// \param psyw   Writeable abstract for the PhaseSpaceSynthesis
__global__ void __launch_bounds__(large_block_size, 1)
kPsyPrimeConjugateGradient(PsSynthesisWriter psyw);

/// \brief Import the Cartesian X, Y, and Z components of positions, velocities, or forces of one
///        system within the synthesis.
///
/// \param x_recv             Long-long integer component of the fixed-precision arrays that the
///                           imported Cartesian X data shall replace
/// \param x_recv_ovrf        Integer overflow component of the fixed-precision arrays that the
///                           imported Cartesian X data shall replace
/// \param y_recv             Long-long integer component of Y fixed-precision arrays
/// \param y_recv_ovrf        Integer overflow component of Y fixed-precision arrays
/// \param z_recv             Long-long integer component of Z fixed-precision arrays
/// \param z_recv_ovrf        Integer overflow component of Z fixed-precision arrays
/// \param box_xform          Box space transformation matrix series within the synthesis (the
///                           matrices have units of inverse Angstroms)
/// \param inverse_xform      Inverse transformation matrix series within the synthesis (the
///                           matrices have units of Angstroms)
/// \param box_dimension      Box dimension series within the synthesis (units of Angstroms)
/// \param box_vectors        Primary bits for the fixed precision box vectors matrix
/// \param box_vector_ovrf    Overflow bits for the fixed precision box vectors matrix
/// \param atom_starts        Starting positions for the atoms of each system in the synthesis
/// \param atom_counts        Counts of atoms in each system of the synthesis
/// \param x_import           Input Cartesian X coordinates (positions, velocities, or forces)
/// \param y_import           Input Cartesian Y coordinates
/// \param z_import           Input Cartesian Z coordinates
/// \param box_xform_in       Transformation matrix to take coordinates into fractional space  
///                           (for positions only--provide nullptr for velocities or forces)
/// \param inverse_xform_in   Transformation matrix to take coordinates back to real space.  The
///                           units of elements in this matrix are Angstroms.
/// \param box_dimensions_in  Dimensions of the box, in internal units of Angstroms
/// \param system_index       Index of the system within this synthesis that the imported
///                           coordinates shall replace
/// \param coversion_factor   Scaling factor to take the input X, Y, and Z data into the apprpriate
///                           fixed-precision scale for incorporation into the synthesis

  
} // namespace synthesis
} // namespace stormm

#endif
