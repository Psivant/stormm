// -*-c++-*-
#ifndef STORMM_HPC_PHASE_SPACE_SYNTHESIS_H
#define STORMM_HPC_PHASE_SPACE_SYNTHESIS_H

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Trajectory/trajectory_enumerators.h"
#include "phasespace_synthesis.h"

namespace stormm {
namespace synthesis {

using card::GpuDetails;
using trajectory::TrajectoryKind;

/// \brief Launch the eponymous kernel to upload scattered data for specific systems within a
///        PhaseSpaceSynthesis object.
///
/// \param destination  Collection of pointers to PhaseSpaceSynthesis data on the host or device.
///                     If on the host, these poiniters must be to host-mapped data visible by the
///                     device.
/// \param source       Collection of pointers to PhaseSpaceSynthesis data on the host or device
///                     If on the host, these poiniters must be to host-mapped data visible by the
///                     device.
/// \param kind         The type of trajectory to upload (positions, velocities, or forces)
/// \param low_index    Lower bound of systems to upload
/// \param high_index   Upper bound of systems to upload (the range is [low_index, high_index))
/// \param gpu          Details of the GPU in use
void systemTransfer(PsSynthesisWriter *destination, PsSynthesisWriter *source,
                    TrajectoryKind kind, int low_index, int high_index, const GpuDetails &gpu);

/// \brief Launch the eponymous kernel to initialize forces on the GPU in a PhaseSpaceSynthesis
///        object.
///
/// \param psyw   Writeable abstract for the PhaseSpaceSynthesis, containing system limits and
///               pointers to all coordinates and forces
/// \param index  Index of the system to initialize; if negative, all systems will be initialized.
/// \param gpu    Details of the GPU in use
void psyInitializeForces(PsSynthesisWriter *psyw, int index, const GpuDetails &gpu);

/// \brief Prepare certain buffers in the phase space (specifically, prior coordinates and
///        velocities) which would otherwise not be used in energy minimization calculations.
///
/// \param psyw  The phase space synthesis to modify
/// \param gpu   Details of the GPU available  
void psyPrimeConjugateGradient(PsSynthesisWriter *psyw, const GpuDetails &gpu);

/// \brief Import the Cartesian X, Y, and Z components of poositions, velocities, or forces of one
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
void psyImportSystemData(llint* x_recv, int* x_recv_ovrf, llint* y_recv, int* y_recv_ovrf,
                         llint* z_recv, int* z_recv_ovrf, double* box_xform, double* inverse_xform,
                         double* box_dimensions, llint* box_vectors, int* box_vector_ovrf,
                         const int* atom_starts, const int* atom_counts, const double* x_import,
                         const double* y_import, const double* z_import,
                         const double* box_xform_in, const double* inverse_xform_in,
                         const double* box_dimensions_in, const int system_index,
                         const TrajectoryKind kind, const double conversion_factor,
                         const GpuDetails &gpu);

} // namespace synthesis
} // namespace stormm

#endif
