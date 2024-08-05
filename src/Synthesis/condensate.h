// -*-c++-*-
#ifndef STORMM_CONDENSATE_H
#define STORMM_CONDENSATE_H

#include <typeinfo>
#include <typeindex>
#include <sys/types.h>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/hybrid_util.h"
#include "Constants/behavior.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "Math/rounding.h"
#include "Math/series_ops.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "phasespace_synthesis.h"
#include "synthesis_enumerators.h"

namespace stormm {
namespace synthesis {

using card::default_hpc_format;
using card::GpuDetails;
using card::Hybrid;
using card::HybridFormat;
using card::HybridTargetLevel;
using constants::CartesianDimension;
using constants::PrecisionModel;
using stmath::incrementingSeries;
using stmath::roundUp;
using topology::UnitCellType;
using trajectory::CoordinateFrame;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;

/// \brief A read-only abstract for the system demarcations in the object.  This information is
///        sometimes critical, and needed on the CPU host even as the general abstract is needed on
///        the GPU device.
struct CondensateBorders {

  /// \brief The constructor accepts the total number of systems as well as pointers to the number
  ///        of atoms and starting indices of each system.
  CondensateBorders(int system_count_in, const size_t* atom_starts_in, const int* atom_counts_in);

  /// \brief Copy and move constructors--as with any object containing const members, the move
  ///        assignment operator is implicitly deleted.
  /// \{
  CondensateBorders(const CondensateBorders &original) = default;
  CondensateBorders(CondensateBorders &&other) = default;
  /// \}
  
  const int system_count;     ///< The total number of systems in the object, and the trusted
                              ///<   length of each of the arrays below
  const size_t* atom_starts;  ///< Starting indices for the atoms of each system
  const int* atom_counts;     ///< Atom counts for each system
};

/// \brief Writeable abstract for the Condensate class, wherein coordinates (only) can be modified
///        as a consequence of certain analyses.
struct CondensateWriter {

  /// \brief The constructor takes a straight list of all relevant constants and pointers.
  CondensateWriter(PrecisionModel mode_in, StructureSource basis_in, int system_count_in,
                   UnitCellType unit_cell_in, const size_t* atom_starts_in,
                   const int* atom_counts_in, float* xcrd_sp, float* ycrd_sp, float* zcrd_sp,
                   double* xcrd, double* ycrd, double* zcrd, double* umat_in, double* invu_in,
                   double* boxdim_in);

  /// \brief The usual copy and move constructors aply for an abstract, with copy and move
  ///        assignment operators being implicitly deleted.
  ///
  /// \param original  The original abstract to copy or move
  /// \{
  CondensateWriter(const CondensateWriter &original) = default;
  CondensateWriter(CondensateWriter &&original) = default;
  /// \}

  const PrecisionModel mode;     ///< The compression mode
  const StructureSource basis;   ///< The original material on which the original object is based
  const int system_count;        ///< The number of systems held by the underlying Condensate
  const UnitCellType unit_cell;  ///< The type of unit cell, common to all systems
  const size_t* atom_starts;     ///< Starting indices of each system's atoms in the data arrays
  const int* atom_counts;        ///< Number of atoms in each system
  float* xcrd_sp;                ///< Single-precision Cartesian X coordinates
  float* ycrd_sp;                ///< Single-precision Cartesian Y coordinates
  float* zcrd_sp;                ///< Single-precision Cartesian Z coordinates
  double* xcrd;                  ///< Double-precision Cartesian X coordinates
  double* ycrd;                  ///< Double-precision Cartesian Y coordinates
  double* zcrd;                  ///< Double-precision Cartesian Z coordinates
  double* umat;                  ///< Box transform information for all systems
  double* invu;                  ///< Inverse box transform information for all systems
  double* boxdims;               ///< Box dimensions for all systems
};

/// \brief Read-only abstract for the Condensate class.  In most cases, the read-only abstract will
///        be preferred as modifications to coordinates should happen in the original synthesis
///        object with its fixed precision arrays.
struct CondensateReader {

  /// \brief The constructor takes a straight list of all relevant constants and pointers.
  ///        One overloaded form accepts the corresponding writer to convert it to read-only form.
  /// \{
  CondensateReader(PrecisionModel mode_in, StructureSource basis_in, int system_count_in,
                   UnitCellType unit_cell_in, const size_t* atom_starts_in,
                   const int* atom_counts_in, const float* xcrd_sp, const float* ycrd_sp,
                   const float* zcrd_sp, const double* xcrd, const double* ycrd,
                   const double* zcrd, const double* umat_in, const double* invu_in,
                   const double* boxdim_in);

  CondensateReader(const CondensateWriter &cdw);
  /// \}

  /// \brief The usual copy and move constructors aply for an abstract, with copy and move
  ///        assignment operators being implicitly deleted.
  ///
  /// \param original  The original abstract to copy or move
  /// \{
  CondensateReader(const CondensateReader &original) = default;
  CondensateReader(CondensateReader &&original) = default;
  /// \}

  const PrecisionModel mode;     ///< The compression mode
  const StructureSource basis;   ///< The original material on which the original object is based
  const int system_count;        ///< The number of systems held by the underlying Condensate
  const UnitCellType unit_cell;  ///< The type of unit cell, common to all systems
  const size_t* atom_starts;     ///< Starting indices of each system's atoms in the data arrays
  const int* atom_counts;        ///< Number of atoms in each system
  const float* xcrd_sp;          ///< Single-precision Cartesian X coordinates
  const float* ycrd_sp;          ///< Single-precision Cartesian Y coordinates
  const float* zcrd_sp;          ///< Single-precision Cartesian Z coordinates
  const double* xcrd;            ///< Double-precision Cartesian X coordinates
  const double* ycrd;            ///< Double-precision Cartesian Y coordinates
  const double* zcrd;            ///< Double-precision Cartesian Z coordinates
  const double* umat;            ///< Box transform information for all systems
  const double* invu;            ///< Inverse box transform information for all systems
  const double* boxdims;         ///< Box dimensions for all systems
};

/// \brief Condense the data format, and possibly offer a reduced representation of coordinates,
///        along with work units detailing which systems should be imported by each block on a
///        given GPU.
class Condensate {
public:

  /// The constructor requires only a complete coordinate synthesis to build the entire object.
  /// The object can be rebuilt post-creation by providing a new PhaseSpaceSynthesis, or
  /// essentially erased by providing a nullptr to conserve memory.  The memory format of a
  /// Condensate is determined by the PhaseSpaceSynthesis or CoordinateSeries object it serves.
  ///
  /// \param poly_ps_in  The coordinate synthesis
  /// \param mode_in     Compression mode to use
  /// \param gpu         Details of the available GPU
  /// \{
  Condensate(HybridFormat format_in = default_hpc_format);
  
  Condensate(const PhaseSpaceSynthesis *poly_ps_in, PrecisionModel mode_in,
             const GpuDetails &gpu = null_gpu);

  Condensate(const PhaseSpaceSynthesis &poly_ps_in, PrecisionModel mode_in,
             const GpuDetails &gpu = null_gpu);

  template <typename T>
  Condensate(const CoordinateSeries<T> *cs_in, PrecisionModel mode_in,
             const GpuDetails &gpu = null_gpu);

  template <typename T>
  Condensate(const CoordinateSeries<T> &cs_in, PrecisionModel mode_in,
             const GpuDetails &gpu = null_gpu);
  /// \}

  /// \brief Copy and move constructors as well as copy and move assignment operator overloads are
  ///        defined in the code, but there is no copy constructor to create objects with a
  ///        different memory format as the mechanics of the Condensate are more complex than those
  ///        of other coordinate objects.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object to fill the assignment statement's right-hand side
  /// \{
  Condensate(const Condensate &original);
  Condensate(Condensate &&original);
  Condensate& operator=(const Condensate &original);
  Condensate& operator=(Condensate &&original);
  /// \}
  
  /// \brief Get the memory layout of the object.
  HybridFormat getFormat() const;

  /// \brief Get the compression mode.
  PrecisionModel getMode() const;

  /// \brief Get the basis for the coordinates in the object: synthesis or series
  StructureSource getBasis() const;

  /// \brief Get the number of systems found in the condensate.
  int getSystemCount() const;
  
  /// \brief Get an indication of whether the Condensate keeps its own copy of the coordinates.
  bool ownsCoordinates() const;

  /// \brief Get the starting index of atoms for one of the systems, using its index.
  ///
  /// \param system_index  Index of the system of interest
  size_t getAtomOffset(int system_index) const;

  /// \brief Get the number of atoms in one of the systems, using its index.
  ///
  /// \param system_index  Index of the system of interest
  int getAtomCount(int system_index) const;
  
  /// \brief Get the data type of the CoordinateSeries upon which this object is based.
  size_t getCoordinateSeriesTypeID() const;
  
  /// \brief Get a const pointer to the original coordinate series.  The developer must supply the
  ///        template type for the CoordinateSeries in order to re-interpret the arbitrarily stored
  ///        <int> type for the CoordinateSeries pointer.
  template <typename T>
  const CoordinateSeries<T>* getSeriesPointer() const;

  /// \brief Get a const pointer to the original coordinate synthesis.
  const PhaseSpaceSynthesis* getSynthesisPointer() const;

  /// \brief Export a single frame as double-precision coordinates.
  ///
  /// \param system_index  Index of the system of interest
  /// \param tier          Obtain information for the frame on the CPU host or the GPU device
  CoordinateFrame exportCoordinateFrame(int system_index,
                                        HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Produce the interlaced coordinates of a specific system.  As is the case with other
  ///        coordinate objects, this will not include the box transformation matrices.
  ///
  /// \param system_index  Index of the system of interest
  /// \param tier          Obtain information for the frame on the CPU host or the GPU device
  std::vector<double>
  getInterlacedCoordinates(int system_index,
                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a pointer to the original object, useful if a pointer is needed but the object
  ///        was passed into a function by const reference.
  const Condensate* getSelfPointer() const;
  
  /// \brief Get the appropriate abstract based on the const-ness of the object.
  ///
  /// Overloaded:
  ///   - Get a reader for a const Condensate
  ///   - Get a writer for a non-const Condensate (the writer can be quickly converted to a reader)
  ///
  /// \param tier  Obtain pointers at the level of the CPU host or the GPU device
  /// \{
  const CondensateReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  CondensateWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}  

  /// \brief Get the read-only summary of the system sizes.  The data is never needed in as a
  ///        collection of device viewable pointers to host-side data, as this information is
  ///        constant as of the creation of the object and therefore consistent on both the CPU
  ///        host and GPU device.
  ///
  /// \param tier  The level (host or device) at which to get the set of pointers
  const CondensateBorders borders(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

#ifdef STORMM_USE_HPC
  /// \brief Get an abstract to the condensate's host data that is guaranteed to be accessible by
  ///        the GPU device.  
  ///
  /// Overloaded:
  ///   - Get a read-only abstract for a const condensate object
  ///   - Get a writeable abstract for a non-const condensate object
  /// \{
  const CondensateReader deviceViewToHostData() const;
  CondensateWriter deviceViewToHostData();
  /// \}
  
  /// \brief Upload all data from the host to the GPU device.
  void upload();

  /// \brief Download all data from the host to the GPU device.
  void download();

  /// \brief Perform a data transfer of coordinates from the coordinate synthesis upon which this
  ///        object is based, at the level of the GPU device.  This function launches a kernel
  ///        and is called by the update() member function below.
  ///
  /// Overloaded:
  ///   - Pull coordinates from a PhaseSpaceSynthesis object
  ///   - Pull coordinates from a templated CoordinateSeries
  ///
  /// \param cs_basis  A reinterpreted cast of the internally stored CoordinateSeries pointer
  /// \param gpu       Details of the available GPU
  /// \{
  void launchCondensateUpdate(const GpuDetails &gpu);
  /// \}
#endif
  
  /// \brief Rebuild the object using a different PhaseSpaceSynthesis or CoordinateSeries, and
  ///        possibly a new compression mode.
  ///
  /// \param poly_ps_in  The new coordinate synthesis
  /// \param cs_in       The new coordinate series
  /// \param mode_in     Compression mode to use (data for any other mode in use will be cleared)
  /// \param gpu         Details of the available GPU
  /// \{
  void rebuild(const PhaseSpaceSynthesis *poly_ps_in,
               PrecisionModel mode_in = PrecisionModel::SINGLE, const GpuDetails &gpu = null_gpu);

  void rebuild(const PhaseSpaceSynthesis &poly_ps_in,
               PrecisionModel mode_in = PrecisionModel::SINGLE, const GpuDetails &gpu = null_gpu);

  template <typename T>
  void rebuild(const CoordinateSeries<T> *cs_in,
               const PrecisionModel mode_in = PrecisionModel::SINGLE,
               const GpuDetails &gpu = null_gpu);
  /// \}

  /// \brief Reload the coordinates based on an updated coordinate synthesis or series, without
  ///        re-computing the block instructions or changing the compression level.
  ///
  /// Overloaded::
  ///   - Rely on internally stored pointers to the original coordinate objects
  ///   - Supply a pointer or reference to the CoordinateSeries object, if applicable
  ///
  /// \param cs_basis  The original coordinate series (will be checked for validity)
  /// \param tier      Obtain pointers at the level of the CPU host or the GPU device
  /// \param gpu       Details of the available GPU
  /// \{
  void update(HybridTargetLevel tier = HybridTargetLevel::HOST, const GpuDetails &gpu = null_gpu);
  
  template <typename T>
  void update(const CoordinateSeries<T> *cs_basis,
              HybridTargetLevel tier = HybridTargetLevel::HOST, const GpuDetails &gpu = null_gpu);
  
  template <typename T>
  void update(const CoordinateSeries<T> &cs_basis,
              HybridTargetLevel tier = HybridTargetLevel::HOST, const GpuDetails &gpu = null_gpu);
  /// \}

private:
  HybridFormat format;             ///< The memory layout of the object, governing whether data is
                                   ///<   available on the GPU device, CPU host, or both  
  PrecisionModel mode;             ///< Mode in which the data from the PhaseSpaceSynthesis object
                                   ///<   is compressed
  StructureSource basis;           ///< Indicate whether the condensate's coordinates are based on
                                   ///<   PhaseSpaceSynthesis or a CoordinateSeries.
  int system_count;                ///< The number of systems held by the Condensate
  UnitCellType unit_cell;          ///< The type of unit cell describing all systems
  bool holds_own_data;             ///< An indication of whether this object has allocated arrays
                                   ///<   to hold a separate copy of the coordinates it represents.
                                   ///<   A Condensate produced (or rebuilt) based on a
                                   ///<   PhaseSpaceSynthesis object will hold its own data, but a
                                   ///<   Condensate built upon a CoordinateSeries object may or
                                   ///<   may not (if not, its coordinate arrays will be set to
                                   ///<   the appropriate coordinate arrays of the series).
  size_t csptr_data_type;          ///< The original data type of the CoordinateSeries pointer,
                                   ///<   stored here for internal reference.  The pointer itself
                                   ///<   (see cs_ptr, below) is re-cast to an arbitrary <int>
                                   ///<   type to prevent the Condensate object itself from taking
                                   ///<   on a template type requirement.
  Hybrid<size_t> atom_starts;      ///< Starting points of each system in the coordinate arrays.
                                   ///<   This is a size_t array rather than an integer array due
                                   ///<   the need to accommodate coordinate series, which could
                                   ///<   exceed the practical 2-billion atom limit for syntheses.
  Hybrid<int> atom_counts;         ///< The number of atoms in each system.  While this replicates
                                   ///<   information from the synthesis or series that the object
                                   ///<   is based upon, the amount of information is trivial, and
                                   ///<   it frees the condensate of template considerations when
                                   ///<   trying to determine a map of the contents.
  Hybrid<float> x_coordinates_sp;  ///< Cartesian X coordinates of all particles (single precision)
  Hybrid<float> y_coordinates_sp;  ///< Cartesian Y coordinates of all particles (single precision)
  Hybrid<float> z_coordinates_sp;  ///< Cartesian Z coordinates of all particles (single precision)
  Hybrid<double> x_coordinates;    ///< Cartesian X coordinates of all particles
  Hybrid<double> y_coordinates;    ///< Cartesian Y coordinates of all particles
  Hybrid<double> z_coordinates;    ///< Cartesian Z coordinates of all particles
  Hybrid<double> box_transforms;   ///< Box space transformation matrices
  Hybrid<double> inv_transforms;   ///< Inverse box space (back to real space) transformation
                                   ///<   matrices
  Hybrid<double> box_dimensions;   ///< Box dimensions (a, b, c, alpha, beta, gamma) for each
                                   ///<   system

  // Pointers to the original coordinate objects
  PhaseSpaceSynthesis *pps_ptr;   ///< Pointer to the PhaseSpaceSynthesis object upon which this
                                  ///<   Condensate and its block instructions are based.  Each
                                  ///<   Condensate will be based on either a PhaseSpaceSynthesis
                                  ///<   or a CoordinateSeries (see cs_ptr, below), the options
                                  ///<   being exclusive.
  CoordinateSeries<int> *cs_ptr;  ///< Pointer to the CoordinateSeries object upon which this
                                  ///<   Condensate and its block instructions are based.  This
                                  ///<   is an arbitrary type and re-interpreted so that a pointer
                                  ///<   to the original source object may be held without placing
                                  ///<   a template requirement on all Condensate objects in any
                                  ///<   situation.

  // ARRAY-kind Hybrid objects targeted by POINTER-kind Hybrids of like types above
  Hybrid<float> float_data;       ///< Targeted by the _sp coordinates Hybrids
  Hybrid<double> double_data;     ///< Targeted by the coordinates and box transformation matrices

  /// \brief Repair pointers to the coordinates held by this object in the event that it has been
  ///        copied.
  void repairPointers();

  /// \brief Validate the system index from within the synthesis.
  ///
  /// \param index  Index of the system of interest
  void validateSystemIndex(int index) const;
};

/// \brief Define the type index for the Condensate object.
static const size_t condensate_type_index = std::type_index(typeid(Condensate)).hash_code();

} // namespace synthesis
} // namespace stormm

#include "condensate.tpp"

// As with common types and STORMM vector types, define the type indices for general use in the
// STORMM namespace.
namespace stormm {
using synthesis::condensate_type_index;
} // namespace stormm

#endif
