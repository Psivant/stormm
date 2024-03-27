// -*-c++-*-
#ifndef HPC_COORDINATE_COPY_H
#define HPC_COORDINATE_COPY_H

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "coordinateframe.h"
#include "coordinate_series.h"
#include "phasespace.h"

namespace stormm {
namespace trajectory {

using card::GpuDetails;
using synthesis::CondensateReader;
using synthesis::CondensateWriter;
using synthesis::PhaseSpaceReader;
using synthesis::PhaseSpaceWriter;
using synthesis::PsSynthesisReader;
using synthesis::PsSynthesisWriter;
  
/// \brief Extract the coordinates and box information of one object on the host or device to
///        another on the host or device.  This will launch a templated CUDA kernel to accomplish
///        the copying and trusts that the coordinate arrays it is being fed will be accessible to
///        the device.
///
/// Overloaded:
///   - Copy between ordinary (up to 64-bit) coordinate representations
///   - Copy from an ordinary representation to a PhaseSpaceSynthesis 95-bit representation
///   - Copy from a PhaseSpaceSynthesis 95-bit representation to an ordinary representation
///
/// \param xdest             Cartesian X coordinates for all particles in the destination system(s)
/// \param ydest             Cartesian X coordinates for all particles in the destination system(s)
/// \param zdest             Cartesian X coordinates for all particles in the destination system(s)
/// \param xdest_ovrf        Overflow data for particle X coordinates in the destination system(s)
/// \param ydest_ovrf        Overflow data for particle Y coordinates in the destination system(s)
/// \param zdest_ovrf        Overflow data for particle Z coordinates in the destination system(s)
/// \param dest_scale        Scaling factor for the destination coordinates
/// \param ct_dest           Specifier code for the destination data type
/// \param xorig             Cartesian X coordinates for all particles in the original system(s)
/// \param yorig             Cartesian Y coordinates for all particles in the original system(s)
/// \param zorig             Cartesian Z coordinates for all particles in the original system(s)
/// \param xorig_ovrf        Overflow data for particle X coordinates in the original system(s)
/// \param yorig_ovrf        Overflow data for particle Y coordinates in the original system(s)
/// \param zorig_ovrf        Overflow data for particle Z coordinates in the original system(s)
/// \param orig_scale        Scaling factor for the origin coordinates
/// \param ct_orig           Specifier code for the origin data type
/// \param natom             The number of atoms' coordinates to copy
/// \param umat_dest         Transformation matrix taking real coordinates into fractional space
/// \param invu_dest         Transformation matrix taking fractional coordinates into real space
/// \param bdim_dest         The copy target box dimensions
/// \param umat_orig         Transformation matrix taking real coordinates into fractional space
/// \param invu_orig         Transformation matrix taking fractional coordinates into real space
/// \param bdim_orig         The origin's box dimensions
/// \param dest_atom_offset  Offset at which to begin writing coordinates for each atom to the
///                          destination
/// \param orig_atom_offset  Offset at which to begin reading coordinates for each atom from the
///                          origin
/// \param dest_xfrm_offset  Offset for reading transformation matrices in the destination
/// \param orig_xfrm_offset  Offset for reading transformation matrices in the origin
/// \param dest_bdim_offset  Offset for reading box dimensions in the destination
/// \param orig_bdim_offset  Offset for reading box dimensions in the origin
/// \param gpu               Details of the GPU to use
/// \{
void launchCopyCoordinateXYZBox(void* xdest, void* ydest, void* zdest, double dest_scale,
                                size_t ct_dest, const void* xorig, const void* yorig,
                                const void* zorig, double orig_scale, size_t ct_orig, int natom,
                                double* umat_dest, double* invu_dest, double* bdim_dest,
                                const double* umat_orig, const double* invu_orig, 
                                const double* bdim_orig, int dest_atom_offset,
                                int orig_atom_offset, int dest_xfrm_offset, int orig_xfrm_offset,
                                int dest_bdim_offset, int orig_bdim_offset, const GpuDetails &gpu);

void launchCopyCoordinateXYZBox(void* xdest, void* ydest, void* zdest, double dest_scale,
                                int dest_scale_bits, size_t ct_dest, const llint* xorig,
                                const llint* yorig, const llint* zorig, const int* xorig_ovrf,
                                const int* yorig_ovrf, const int* zorig_ovrf, double orig_scale,
                                int orig_scale_bits, int natom, double* umat_dest,
                                double* invu_dest, double* bdim_dest, const double* umat_orig,
                                const double* invu_orig, const double* bdim_orig,
                                int dest_atom_offset, int orig_atom_offset, int dest_xfrm_offset,
                                int orig_xfrm_offset, int dest_bdim_offset, int orig_bdim_offset,
                                const GpuDetails &gpu);

void launchCopyCoordinateXYZBox(llint* xdest, llint* ydest, llint* zdest, int* xdest_ovrf,
                                int* ydest_ovrf, int* zdest_ovrf, double dest_scale,
                                int dest_scale_bits, const void* xorig, const void* yorig,
                                const void* zorig, double orig_scale, int orig_scale_bits,
                                size_t ct_orig, int natom, double* umat_dest, double* invu_dest,
                                double* bdim_dest, const double* umat_orig,
                                const double* invu_orig, const double* bdim_orig,
                                llint* boxvecs_dest, int* boxvec_ovrf_dest, int dest_atom_offset,
                                int orig_atom_offset, int dest_xfrm_offset, int orig_xfrm_offset,
                                int dest_bdim_offset, int orig_bdim_offset, const GpuDetails &gpu);
/// \}

/// \brief Extract the coordinates of one object on the host or device to another on the host or
///        device.  This will launch a templated CUDA kernel to accomplish the copying.  This
///        functon trusts that the coordinate arrays it is being fed will be accessible to the
///        device.  Descriptions of parameters follow from launchCopyCoordinateXYZ() above.
///
/// Overloaded:
///   - Copy between ordinary (up to 64-bit) coordinate representations
///   - Copy from an ordinary representation to a PhaseSpaceSynthesis 95-bit representation
///   - Copy from a PhaseSpaceSynthesis 95-bit representation to an ordinary representation
/// \{
void launchCopyCoordinateXYZ(void* xdest, void* ydest, void* zdest, double dest_scale,
                             size_t ct_dest, const void* xorig, const void* yorig,
                             const void* zorig, double orig_scale, size_t ct_orig, int natom,
                             int dest_atom_offset, int orig_atom_offset, const GpuDetails &gpu);

void launchCopyCoordinateXYZ(void* xdest, void* ydest, void* zdest, double dest_scale,
                             int dest_scale_bits, size_t ct_dest, const llint* xorig,
                             const llint* yorig, const llint* zorig, const int* xorig_ovrf,
                             const int* yorig_ovrf, const int* zorig_ovrf, double orig_scale,
                             int orig_scale_bits, int natom, int dest_atom_offset,
                             int orig_atom_offset, const GpuDetails &gpu);

void launchCopyCoordinateXYZ(llint* xdest, llint* ydest, llint* zdest, int* xdest_ovrf,
                             int* ydest_ovrf, int* zdest_ovrf, double dest_scale,
                             int dest_scale_bits, const void* xorig, const void* yorig,
                             const void* zorig, double orig_scale, int orig_scale_bits,
                             size_t ct_orig, int natom, int dest_atom_offset,
                             int orig_atom_offset, const GpuDetails &gpu);
/// \}

/// \brief Launch special kernels for copying coordinate objects with multiple components
///        (positons, velocities, forces) or multiple coordinate sets (PhaseSpaceSynthesis,
///        Condensate, or CoordinateSeries objects).  All abstracts and any array pointers
///        supplied to this function must be accessible by the GPU device.
///
/// Overloaded:
///   - Copy one system between selected coordinate objects
///   - Copy multiple systems between selected coordinate objects
///
/// \param destination       Abstract of the object into which coordinates shall be placed
/// \param origin            Abstract of the object from which coordinates are taken
/// \param dest_atom_offset  The starting index at which to begin writing atomic coordinates in the
///                          destination object.  This is needed as the abstracts are not
///                          guaranteed to be valid on the host, and best obtained on the host and
///                          passed to the kernel to avoid a blocking, high-latency memory call.
/// \param orig_atom_offset  The starting index at which to begin writing atomic coordinates in the
///                          origin object.
/// \param index_orig        Index of the system in the origin coordinate synthesis
/// \param frame_orig        Index of the system in the origin coordinate series (the distinction
///                          between "index" and "frame" is semantic, but helps to indicate whether
///                          the function is working with a synthesis or series)
/// \param kind              Specify whether to copy information into the position, velocity, or
///                          force arrays of a PhaseSpaceSynthesis (the time cycle is handled
///                          when making the object's abstract)
/// \param system_pairs      List of system indices in either object to copy.  The indices of the
///                          origin and destination systems are given in the "x" and "y" members
///                          of the tuple, respectively.
/// \param copy_count        Trusted length of system_pairs
/// \param gpu               Details of the GPU that will manage the transfer
/// \{
void launchCopyCoordinates(PhaseSpaceWriter *destination, const PhaseSpaceReader &origin,
                           const GpuDetails &gpu);

void launchCopyCoordinates(PhaseSpaceWriter *destination, const PsSynthesisReader &origin,
                           size_t orig_atom_offset, int index_orig, const GpuDetails &gpu);

void launchCopyCoordinates(PsSynthesisWriter *destination, size_t dest_atom_offset, int index_dest,
                           const PhaseSpaceReader &origin, const GpuDetails &gpu);

void launchCopyCoordinates(PsSynthesisWriter *destination, size_t dest_atom_offset, int index_dest,
                           const PsSynthesisReader &origin, size_t orig_atom_offset,
                           int index_orig, int natom, const GpuDetails &gpu);

void launchCopyCoordinates(CoordinateSeriesWriter<void> *destination, size_t ct_dest,
                           const CoordinateSeriesReader<void> &origin, size_t ct_orig,
                           const int2* system_pairs, int copy_count, const GpuDetails &gpu);

void launchCopyCoordinates(CoordinateSeriesWriter<void> *destination, size_t ct_dest,
                           const PsSynthesisReader &origin, TrajectoryKind kind,
                           const int2* system_pairs, int copy_count, const GpuDetails &gpu);

void launchCopyCoordinates(CoordinateSeriesWriter<void> *destination, size_t ct_dest,
                           const CondensateReader &origin, const int2* system_pairs,
                           int copy_count, const GpuDetails &gpu);

void launchCopyCoordinates(PsSynthesisWriter *destination, TrajectoryKind kind,
                           const CoordinateSeriesReader<void> &origin, size_t ct_orig,
                           const int2* system_pairs, int copy_count, const GpuDetails &gpu);

void launchCopyCoordinates(PsSynthesisWriter *destination, const PsSynthesisReader &origin,
                           const int2* system_pairs, int copy_count, const GpuDetails &gpu);

void launchCopyCoordinates(PsSynthesisWriter *destination, TrajectoryKind kind,
                           const CondensateReader &origin, const int2* system_pairs,
                           int copy_count, const GpuDetails &gpu);

void launchCopyCoordinates(CondensateWriter *destination,
                           const CoordinateSeriesReader<void> &origin, size_t ct_orig,
                           const int2* system_pairs, int copy_count, const GpuDetails &gpu);

void launchCopyCoordinates(CondensateWriter *destination, const PsSynthesisReader &origin,
                           TrajectoryKind kind, const int2* system_pairs, int copy_count,
                           const GpuDetails &gpu);

void launchCopyCoordinates(CondensateWriter *destination, const CondensateReader &origin,
                           const int2* system_pairs, int copy_count, const GpuDetails &gpu);
/// \}

} // namespace trajectory
} // namespace stormm

#endif
