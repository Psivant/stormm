// -*-c++-*-
#ifndef STORMM_COORDINATE_SWAP_H
#define STORMM_COORDINATE_SWAP_H

#include <algorithm>
#include "copyright.h"
#include "Accelerator/card_utilities.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Reporting/error_format.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Topology/atomgraph_enumerators.h"
#include "coordinate_series.h"

namespace stormm {
namespace trajectory {

using card::GpuDetails;
using card::HpcKernelSync;
using card::Hybrid;
using card::HybridTargetLevel;
using synthesis::Condensate;
using synthesis::CondensateBorders;
using synthesis::CondensateWriter;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisBorders;
using synthesis::PsSynthesisWriter;
using topology::UnitCellType;

/// \brief Check that the fixed-precision scaling of two objects is the same.  This can be called
///        for more than one aspect of each object if needed.
///
/// \param first_bits   Fixed-precision bits after the decimal in the first object's data
///                     (real-valued data will have a bit count of zero)
/// \param second_bits  Fixed-precision bits after the decimal in the second object's data
/// \param desc         Description of the quantity rendered in fixed precision
void swapValidatePrecisionModels(int first_bits, int second_bits, const char* desc);

/// \brief Check the frame indices within two CoordinateSeries or the system indices within two
///        coordinate syntheses for validity.
///
/// Overloaded:
///   - Check a pair of systems
///   - Check multiple pairs of systems (the host-side data in each Hybrid object holding the
///     indices will be taken as authoritative)
///
/// \param frm_first     The frame or system index of the first object, or the list of frame or
///                      system indices of the first object
/// \param count_first   The number of frames or systems in the first object
/// \param frm_second    The frame or system index of the second object
/// \param frm_pairs     The list of frame or system indices in each object.  Indices in the first
///                      object are given in the "x" member of the tuple and indices of the second
///                      objects are given in the "y" member of the tuple.
/// \param pair_count    The number of frame pairs to compare (the trusted length of frm_pairs)
/// \param count_second  The number of frames or systems in the second object
/// \{
void swapValidateFrameIndices(int frm_first, int count_first, int frm_second, int count_second);

void swapValidateFrameIndices(const int2* frm_pairs, int pair_count, int count_first,
                              int count_second);
/// \}

/// \brief Check that the systems being swapped have similar unit cell types (at least that they
///        have periodic boundary conditions or not--there could be corner cases in which it is
///        desirable to swap between a system that is orthorhombic and one which is triclinic, if
///        for example the box angles are nearly right angles).
///
/// \param uca  The unit cell type of the second collection of structures
/// \param ucb  The unit cell type of the second collection of structures
void swapValidateUnitCells(UnitCellType uca, UnitCellType ucb);
  
/// \brief Apply additional considerations when swapping a batch of systems between the same memory
///        tier of the same object, to avoid possible race conditions on the GPU or dependence on
///        the order of swaps on the CPU.
///
/// \param is_same      Calls to the function can provide some comparison of the pointers to each
///                     object involved in the swap for this input.  It is a boolean to generalize
///                     the functionality, without templating and without the need to handle a
///                     templated type within a templated type in the case of swaps between
///                     CoordinateSeries objects.
/// \param frm_pairs    The list of all system index swaps to perform
/// \param pair_count   The number of swap pairs in the batch
/// \param tier_first   The first tier of memory (CPU host or GPU device) involved in the swap
/// \param tier_second  The second tier of memory involved in the swap
void swapCheckSameObject(bool is_same, const int2* frm_pairs, int pair_count,
                         HybridTargetLevel tier_first, HybridTargetLevel tier_second);
  
/// \brief Check the atom counts of two objects to ensure that the frames or systems being 
///        swapped are valid.  Overloading and descriptions of input variables follow from
///        swapValidateFrameIndices() above, in addition to:
///
/// \param natom_first   Atom count of each frame in the first CoordinateSeries object, or an
///                      array of atom counts for system in some coordinate synthesis (the
///                      indices of systems to swap will be assumed valid)
/// \param natom_second  Atom count of each frame in the second CoordinateSeries object, or an
///                      array of atom counts for system in some coordinate synthesis (the
///                      indices of systems to swap will be assumed valid)
/// \param system_count  The number of systems to swap, and the trusted length of both
///                      frm_first and frm_second
/// \{
void swapValidateAtomCounts(int natom_first, int natom_second);

template <typename T>
void swapValidateAtomCounts(int frm_first, const T* natom_first, int frm_second,
                            const T* natom_second);

template <typename T>
void swapValidateAtomCounts(const int2* frm_pairs, int pair_count, const T* natom_first,
                            const T* natom_second);
/// \}

/// \brief Swap coordinates and box information on the host.  All coordinate objects include
///        box information, but as the PhaseSpaceSynthesis object also includes velocity and force
///        information, the box information will only be copied if supplied a non-null pointer.
///
/// Overloaded:
///   - Accept pointers to the primary Cartesian coordinate X, Y, and Z arrays
///   - Accept additional pointers to the overflow arrays for X, Y, and Z coordinates
///   - Accept pointers to the box information, or leave as nullptr to skip box information swap
///   - Accept a pair of void-casted CoordinateSeries abstracts, restore their types, and delegate
///     to the appropriate variant of this function
///
/// \param first_x       Cartesian X coordinates of the first system (these can be positions,
///                      velocities, or forces)
/// \param first_y       Cartesian Y coordinates of the first system
/// \param first_z       Cartesian Z coordinates of the first system
/// \param first_ofs     Offset index for the first system
/// \param second_ofs    Offset index for the second system
/// \param natom         The number of atoms in the systems to swap
/// \param frm_first     Frame number to swap from the first system
/// \param frm_second    Frame number to swap from the second system
/// \param first_xovrf   Overflow bits for fixed-precision representations of the Cartesian X
///                      coordinates of the first system (used only with the PhaseSpaceSynthesis)
/// \param first_yovrf   Overflow bits for Cartesian Y coordinates of the first system
/// \param first_zovrf   Overflow bits for Cartesian Z coordinates of the first system
/// \param second_xovrf  Overflow bits for Cartesian X coordinates of the second system
/// \param second_yovrf  Overflow bits for Cartesian Y coordinates of the second system
/// \param second_zovrf  Overflow bits for Cartesian Z coordinates of the second system
/// \param first_umat    Transformation matrices to take systems in the first object into unit
///                      cell fractional space.  Matrices for each system are padded by the warp
///                      stride.
/// \param first_invu    Transformation matrices to take systems in the first object back into
///                      Cartesian space
/// \param first_bdim    Box dimensions for systems in the first object.  These sets of six
///                      dimensions are likewise padded by the warp stride.
/// \param second_umat   Transformation matrices to take systems in the second object into unit
///                      cell fractional space
/// \param second_invu   Transformation matrices to take systems in the second object back into
///                      Cartesian space
/// \param second_bdim   Box dimensions for systems in the second object
/// \param v_first       Void-casted abstract of the first coordinate series
/// \param v_second      Void-casted abstract of the second coordinate series
/// \param ct            Codified type identifier for void-casted coordinate series abstracts
/// \{
template <typename T>
void swapCoordinates(T* first_x, T* first_y, T* first_z, T* second_x, T* second_y, T* second_z,
                     size_t first_ofs, size_t second_ofs, size_t natom, size_t frm_first = 0,
                     size_t frm_second = 0, int* first_xovrf = nullptr, int* first_yovrf = nullptr,
                     int* first_zovrf = nullptr, int* second_xovrf = nullptr,
                     int* second_yovrf = nullptr, int* second_zovrf = nullptr,
                     double* first_umat = nullptr, double* first_invu = nullptr,
                     double* first_bdim = nullptr, double* second_umat = nullptr,
                     double* second_invu = nullptr, double* second_bdim = nullptr);

template <typename T>
void swapCoordinates(T* first_x, T* first_y, T* first_z, T* second_x, T* second_y, T* second_z,
                     size_t first_ofs, size_t second_ofs, size_t natom, size_t frm_first,
                     size_t frm_second, double* first_umat, double* first_invu, double* first_bdim,
                     double* second_umat, double* second_invu, double* second_bdim);

template <typename T>
void swapCoordinates(CoordinateSeriesWriter<void> *v_first, size_t frm_first,
                     CoordinateSeriesWriter<void> *v_second, size_t frm_second);

template <typename T>
void swapCoordinates(CoordinateSeriesWriter<void> *v_first, CoordinateSeriesWriter<void> *v_second,
                     const int2* frm_pairs, int pair_count);
/// \}

/// \brief Swap the coordinates between two equivalent objects.  A GPU kernel will be engaged to
///        perform swaps between the CPU host and GPU device, or between two locations in device
///        memory.  Otherwise, host-based C++ functions will be use.  If swapping systems in a
///        PhaseSpaceSynthesis, all components of each system (positions, velocities, and forces)
///        will be swapped in order to ensure the integrity of each system.
///
/// Overloaded:
///   - Copy between similarly sized coordinate sets on any multi-coordinate object
///   - Specify whether the origin and destination coordinates lie on the device or the host
///   - Swap a single system or multiple systems
///
/// \param first        The first object involved in the swap
/// \param frm_first    Frame index of the first coordinate series
/// \param idx_first    System index of the first coordinate synthesis
/// \param second       The second coordinate object involved in the swap
/// \param frm_second   Frame index of the second coordinate series
/// \param idx_second   System index of the second coordinate synthesis
/// \param tier_first   Specify whether to obtain coordinates from the CPU host or GPU device in
///                     the first object
/// \param tier_second  Specify whether to obtain coordinates from the CPU host or GPU device in
///                     the second object
/// \param gpu          Details of the GPU which may perform the swap
/// \param sync         Directive for synchronzing the CPU host and the GPU device before and / or
///                     after the swap
/// \{
void coordSwap(CoordinateSeriesWriter<void> *first, size_t frm_first,
               CoordinateSeriesWriter<void> *second, size_t frm_second, size_t ct,
               HybridTargetLevel tier_first = HybridTargetLevel::HOST,
               HybridTargetLevel tier_second = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordSwap(CoordinateSeriesWriter<void> *first, CoordinateSeriesWriter<void> *second,
               const int2* frm_pairs, int pair_count, size_t ct,
               HybridTargetLevel tier_first = HybridTargetLevel::HOST,
               HybridTargetLevel tier_second = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordSwap(CoordinateSeries<T> *first, size_t frm_first, CoordinateSeries<T> *second,
               size_t frm_second, HybridTargetLevel tier_first = HybridTargetLevel::HOST,
               HybridTargetLevel tier_second = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

template <typename T>
void coordSwap(CoordinateSeries<T> *first, CoordinateSeries<T> *second,
               const Hybrid<int2> &frm_pairs,
               HybridTargetLevel tier_first = HybridTargetLevel::HOST,
               HybridTargetLevel tier_second = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordSwap(PsSynthesisWriter *first, int idx_first, PsSynthesisWriter *second,
               int idx_second, const PsSynthesisBorders &first_bdrs,
               const PsSynthesisBorders &second_bdrs,
               HybridTargetLevel tier_first = HybridTargetLevel::HOST,
               HybridTargetLevel tier_second = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordSwap(PsSynthesisWriter *first, PsSynthesisWriter *second, const int2* idx_pairs,
               int pair_count, const PsSynthesisBorders &first_bdrs,
               const PsSynthesisBorders &second_bdrs,
               HybridTargetLevel tier_first = HybridTargetLevel::HOST,
               HybridTargetLevel tier_second = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordSwap(PhaseSpaceSynthesis *first, int idx_first, PhaseSpaceSynthesis *second,
               int idx_second, HybridTargetLevel tier_first = HybridTargetLevel::HOST,
               HybridTargetLevel tier_second = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordSwap(PhaseSpaceSynthesis *first, PhaseSpaceSynthesis *second,
               const Hybrid<int2> &idx_pairs,
               HybridTargetLevel tier_first = HybridTargetLevel::HOST,
               HybridTargetLevel tier_second = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordSwap(CondensateWriter *first, size_t idx_first, CondensateWriter *second,
               size_t idx_second, const CondensateBorders &first_bdrs,
               const CondensateBorders &second_bdrs,
               HybridTargetLevel tier_first = HybridTargetLevel::HOST,
               HybridTargetLevel tier_second = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordSwap(CondensateWriter *first, CondensateWriter *second, const int2* idx_pairs,
               int pair_count, const CondensateBorders &first_bdrs,
               const CondensateBorders &second_bdrs,
               HybridTargetLevel tier_first = HybridTargetLevel::HOST,
               HybridTargetLevel tier_second = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordSwap(Condensate *first, size_t idx_first, Condensate *second, size_t idx_second,
               HybridTargetLevel tier_first = HybridTargetLevel::HOST,
               HybridTargetLevel tier_second = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);

void coordSwap(Condensate *first, Condensate *second, const Hybrid<int2> &idx_pairs,
               HybridTargetLevel tier_first = HybridTargetLevel::HOST,
               HybridTargetLevel tier_second = HybridTargetLevel::HOST,
               const GpuDetails &gpu = null_gpu, HpcKernelSync sync = HpcKernelSync::MEMORY_AUTO);
/// \}
  
} // namespace trajetory
} // namespace stormm

#include "coordinate_swap.tpp"

#endif
