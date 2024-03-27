// -*-c++-*-
#ifndef STORMM_COORDINATE_SERIES_H
#define STORMM_COORDINATE_SERIES_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "FileManagement/file_enumerators.h"
#include "FileManagement/file_listing.h"
#include "FileManagement/file_util.h"
#include "Math/rounding.h"
#include "amber_ascii.h"
#include "coordinateframe.h"
#include "phasespace.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

using card::Hybrid;
using card::HybridTargetLevel;
using constants::CartesianDimension;
using data_types::isFloatingPointScalarType;
using data_types::isSignedIntegralScalarType;
using diskutil::detectCoordinateFileKind;
using diskutil::DataFormat;
using diskutil::DrivePathType;
using diskutil::getDrivePathType;
using diskutil::getTrajectoryFormat;
using numerics::default_trajpos_scale_bits;
using stmath::roundUp;
  
/// \brief A simple list of all the valid type specifiers for coordinate data.  This will be filled
///        at runtime based on the type IDs found (see DataTypes/common_types.h).
struct ValidCoordinateTypes {

  ValidCoordinateTypes(size_t double_id_in, size_t float_id_in, size_t short_id_in,
                       size_t int_id_in, size_t llint_id_in);

  /// \brief The copy and move constructors as well as assignment operators can all take their
  ///        default forms for this simple object with no arrays or const members.
  /// \{
  ValidCoordinateTypes(const ValidCoordinateTypes &original) = default;
  ValidCoordinateTypes(ValidCoordinateTypes &&original) = default;
  ValidCoordinateTypes& operator=(const ValidCoordinateTypes &original) = default;
  ValidCoordinateTypes& operator=(ValidCoordinateTypes &&original) = default;
  /// \}
  
  size_t double_id;  ///< Double-precision real data type identifier
  size_t float_id;   ///< Single-precision real data type identifier
  size_t short_id;   ///< Short (16-bit integer) data type identifier
  size_t int_id;     ///< Typical / long (32-bit integer) data type identifier
  size_t llint_id;   ///< Long long (64-bit integer) data type identifier
};
  
/// \brief Collect C-style pointers and critical constants for a writeable CoordinateSeries object.
template <typename T> struct CoordinateSeriesWriter {

  /// \brief The constructor feeds all arguments straight to the inline initialization list.
  CoordinateSeriesWriter(int natom_in, int nframe_in, UnitCellType unit_cell_in, int gpos_bits_in,
                         double gpos_scale_in, double inv_gpos_scale_in, T* xcrd_in, T* ycrd_in,
                         T* zcrd_in, double* umat_in, double* invu_in, double* boxdim_in);

  /// \brief Copy and move constructors.  The move assignment operator is implicitly deleted.
  /// \{
  CoordinateSeriesWriter(const CoordinateSeriesWriter &original) = default;
  CoordinateSeriesWriter(CoordinateSeriesWriter &&original) = default;
  CoordinateSeriesWriter& operator=(const CoordinateSeriesWriter &other) = default;
  /// \}

  const int natom;               ///< The number of atoms in the system
  const int nframe;              ///< The number of frames in the series
  const UnitCellType unit_cell;  ///< The type of unit cell (i.e. ORTHORHOMBIC, could also be NONE)
  const int gpos_bits;           ///< Global position coordinate bits after the decimal
  const double gpos_scale;       ///< Global position coordinate scaling factor
  const double inv_gpos_scale;   ///< Inverse global coordinate scaling factor
  T* xcrd;                       ///< Cartesian X coordinates of all atoms
  T* ycrd;                       ///< Cartesian Y coordinates of all atoms
  T* zcrd;                       ///< Cartesian Z coordinates of all atoms
  double* umat;                  ///< Transformation matrix to take coordinates into box
                                 ///<   (fractional) space
  double* invu;                  ///< Inverse transformation matrix out of box space
  double* boxdim;                ///< Box dimensions (these will be consistent with umat and invu)
};
  
/// \brief Collect C-style pointers and critical constants for a read-only CoordinateSeries object.
template <typename T> struct CoordinateSeriesReader {

  /// \brief The basic constructor feeds all arguments straight to the inline initialization list.
  ///        An alternative constructor takes a writeable abstract and converts it to a read-only
  ///        abstract.
  /// \{
  CoordinateSeriesReader(int natom_in, int nframe_in, UnitCellType unit_cell_in, int gpos_bits_in,
                         double gpos_scale_in, double inv_gpos_scale_in, const T* xcrd_in,
                         const T* ycrd_in, const T* zcrd_in, const double* umat_in,
                         const double* invu_in, const double* boxdim_in);
  CoordinateSeriesReader(const CoordinateSeriesWriter<T> &csw);
  /// \}

  /// \brief Copy and move constructors.  The move assignment operator is implicitly deleted.
  /// \{
  CoordinateSeriesReader(const CoordinateSeriesReader &original) = default;
  CoordinateSeriesReader(CoordinateSeriesReader &&original) = default;
  CoordinateSeriesReader& operator=(const CoordinateSeriesReader &other) = default;
  /// \}

  const int natom;               ///< The number of atoms in the system
  const int nframe;              ///< The number of frames in the series
  const UnitCellType unit_cell;  ///< The type of unit cell (i.e. ORTHORHOMBIC, could also be NONE)
  const int gpos_bits;           ///< Global position coordinate bits after the decimal
  const double gpos_scale;       ///< Global position coordinate scaling factor
  const double inv_gpos_scale;   ///< Inverse global coordinate scaling factor
  const T* xcrd;                 ///< Cartesian X coordinates of all atoms
  const T* ycrd;                 ///< Cartesian Y coordinates of all atoms
  const T* zcrd;                 ///< Cartesian Z coordinates of all atoms
  const double* umat;            ///< Transformation matrix to take coordinates into box
                                 ///<   (fractional) space
  const double* invu;            ///< Inverse transformation matrix out of box space
  const double* boxdim;          ///< Box dimensions (these will be consistent with umat and invu)
};

/// \brief Store the coordinates and box information for a series of frames, in one of several
///        levels of precision.  Individual frames can be extracted into CoordinateFrame, and
///        PhaseSpace, and PhaseSpaceSynthesis objects, or new frames can be added from
///        CoordinateFrame objects.  This object is not the CoordinateFrame equivalent of a
///        PhaseSpaceSynthesis object, however: its purpose is to collect only the coordinates of
///        many frames of a single system at the appropriate level of precision.
template <typename T>  class CoordinateSeries {
public:

  /// \brief There are several options for constructing this collection of coordinate frames.
  ///
  /// Overloaded:
  ///   - Allocate to hold a given number of atoms and frames
  ///   - Create from any of the coordinate file formats (restart and input coordinate formats
  ///     will create only one frame, but this limit can be increased later)
  ///   - From an existing PhaseSpace or CoordinateFrame object (as a pointer or a copy, with the
  ///     option to make many copies immediately)
  ///
  /// \param natom_in          The number of atoms expected
  /// \param nframe_in         Initial number of frames to allocate for
  /// \param unit_cell_in      The type of unit cell to prepare for (this can be modified after
  ///                          creating the object)
  /// \param file_name         File to read from
  /// \param frame_numbers_in  Frame numbers of the file to read (default all frames)
  /// \param replica_count_in  The number of times to replicate a series of one of more frames
  ///                          read from a file (this is useful for immediately making many copies
  ///                          of a restrart or input coordinates file, which has only one frame)
  /// \param atom_count_in     The number of atoms to expect (critical if the input file is an
  ///                          Amber .crd format trajectory, otherwise can be left at zero to be
  ///                          filled in when the file is read)
  /// \param ps                Pre-existing object with a complete coordinate set to use as a
  ///                          template.  When constructing from a pre-existing PhaseSpace object,
  ///                          nframe_in indicates a number of copies to allocate for and create.
  /// \param cf                Pre-existing object with a complete coordinate set to use as a
  ///                          template.  When constructing from a pre-existing CoordinateFrame,
  ///                          nframe_in indicates a number of copies to allocate for and create.
  /// \{
  explicit CoordinateSeries(int natom_in = 0, int nframe_in = 0,
                            UnitCellType unit_cell_in = UnitCellType::NONE,
                            int globalpos_scale_bits_in = 0);
  explicit CoordinateSeries(const std::string &file_name, int atom_count_in = 0,
                            CoordinateFileKind file_kind = CoordinateFileKind::UNKNOWN,
                            const std::vector<int> &frame_numbers = {},
                            int replica_count = 1, UnitCellType unit_cell_in = UnitCellType::NONE,
                            int globalpos_scale_bits_in = 0);
  explicit CoordinateSeries(PhaseSpace *ps, int nframe_in, int globalpos_scale_bits_in = 0);
  explicit CoordinateSeries(const PhaseSpace &ps, int nframe_in, int globalpos_scale_bits_in = 0);
  explicit CoordinateSeries(CoordinateFrame *cf, int nframe_in, int globalpos_scale_bits_in = 0);
  explicit CoordinateSeries(const CoordinateFrame &cf, int nframe_in,
                            int globalpos_scale_bits_in = 0);
  /// \}

  /// \brief Use the default copy and move constructors, copy and move assignment operators.  This
  ///        object has no const members to trigger implicit deletions and no POINTER-kind Hybrid
  ///        objects to repair.
  ///
  /// \param original  The original object, being copied or moved
  /// \param other     The other object, being copied or moved by assignment
  /// \{
  CoordinateSeries(const CoordinateSeries &original) = default;
  CoordinateSeries(CoordinateSeries &&original) = default;
  CoordinateSeries& operator=(const CoordinateSeries &other) = default;
  CoordinateSeries& operator=(CoordinateSeries &&other) = default;
  /// \}

  /// \brief Special copy constructor to take a CoordinateSeries of one data type into another.
  ///        The implementation requires a special declaration (template <typename T> template
  ///        <typename Toriginal>...) with two uses of the template keyword because of the nested
  ///        nature of these templates.  In other contexts template <typename T,
  ///        typename Toriginal> would suffice.
  ///
  /// \param original  The original object, being copied or moved
  /// \{
  template <typename Toriginal>
  CoordinateSeries(const CoordinateSeries<Toriginal> &original, int globalpos_scale_bits_in = 0);
  /// \}

  /// \brief Get the number of atoms in each frame of the series.
  int getAtomCount() const;

  /// \brief Get the number of frames in the series.
  int getFrameCount() const;

  /// \brief Get the maximum number of frames that the object can hold.
  int getFrameCapacity() const;

  /// \brief Get the unit cell type of the coordinate system
  UnitCellType getUnitCellType() const;

  /// \brief Get the fixed precision bits after the decimal (if applicable, warn if the data type
  ///        is non-integer).
  int getFixedPrecisionBits() const;
  
  /// \brief Get the interlaced coordinates of one frame.
  ///
  /// Overloaded:
  ///   - Get coordinates for all atoms in a frame
  ///   - Get coordinates for a selected range of atoms in a frame
  ///
  /// \param frame_index  Index of the frame to access
  /// \param low_index    The lower atom index of a range
  /// \param high_index   The upper atom index of a range
  /// \param tier  The level at which to access coordinates
  /// \{
  template <typename Treport> std::vector<Treport>
  getInterlacedCoordinates(int frame_index, int globalpos_bits_out = -1,
                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  template <typename Treport> std::vector<Treport>
  getInterlacedCoordinates(int frame_index, int low_index, int high_index,
                           int globalpos_bits_out = -1,
                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Get the transformation matrix to take coordinates into box (fractional) space.  The
  ///        result should be interpreted as a 3x3 matrix in column-major format.
  ///
  /// \param frame_index  The frame for which to extract the transformation matrix
  /// \param tier         Level at which to retrieve data (if STORMM is compiled to run on a GPU)
  std::vector<double> getBoxSpaceTransform(int frame_index,
                                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the transformation matrix to take coordinates from fractional space back into
  ///        real space.  The result should be interpreted as a 3x3 matrix in column-major format.
  ///
  /// \param frame_index  The frame for which to extract the transformation matrix
  /// \param tier         Level at which to retrieve data (if STORMM is compiled to run on a GPU)
  std::vector<double> getInverseTransform(int frame_index,
                                          HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the box dimensions in their pure form for a particular frame.
  ///
  /// Overloaded:
  ///   - Get the dimensions for a specific system on the CPU host or GPU device
  ///   - Get the dimensions for all systems at both levels
  ///
  /// \param frame_index  The frame for which to extract the transformation matrix
  /// \param tier         Level at which to retrieve data (if STORMM is compiled to run on a GPU)
  /// \{
  std::vector<double> getBoxDimensions(int frame_index,
                                       HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  const Hybrid<double>& getBoxDimensions() const;
  /// \}

  /// \brief Extract coordinates to a pre-existing object.
  ///
  /// Overloaded:
  ///   - Load coordinates into a CoordinateFrame object
  ///   - Load coordinates into a PhaseSpace object using the stated point in the time cycle or
  ///     taking the object's current orientation
  ///
  /// \param cf           Coordinate frame into which coordinates shall be placed
  /// \param ps           Phase space object into which coordinates shall be placed
  /// \param frame_index  Frame to use in coordinate extraction
  /// \{
  void extractFrame(CoordinateFrame *cf, int frame_index,
                    HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  void extractFrame(PhaseSpace *ps, int frame_index,
                    TrajectoryKind kind, CoordinateCycle time_point,
                    HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  void extractFrame(PhaseSpace *ps, int frame_index,
                    TrajectoryKind kind, HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  void extractFrame(PhaseSpace *ps, int frame_index,
                    HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Prepare a CoordinateFrame object based on one frame of the series, accomplishing any
  ///        necessary data conversions to put the coordinates back into the familiar
  ///        double-precision format.
  ///
  /// \param frame_index  The frame to extract
  /// \param tier         Level at which to retrieve data (if STORMM is compiled to run on a GPU)
  CoordinateFrame exportFrame(int frame_index,
                              HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Prepare a PhaseSpace object based on one frame of the series, accomplishing any data
  ///        conversions needed to put the coordinates back into the familiar double-precision
  ///        format.
  ///
  /// \param frame_index  The frame to extract
  /// \param tier         Level at which to retrieve data (if STORMM is compiled to run on a GPU)
  PhaseSpace exportPhaseSpace(int frame_index,
                              HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
  /// \brief Export the contents of this coordinate series to a trajectory file.
  ///
  /// \param file_name     Name of the file to write
  /// \param output_kind   The format of the file to write (checkpoint files print position and
  ///                      velocity data by obligation, but trajectory files can contain either of
  ///                      these as well as forces)
  /// \param expectation   The condition in which the output file is expected to be found
  /// \param low_index     The first frame index to export
  /// \param high_index    The upper limit of frame indices to export
  /// \param tier          Level at which to retrieve data (if STORMM is compiled to run on a GPU)
  void exportToFile(const std::string &file_name,
                    CoordinateFileKind output_kind = CoordinateFileKind::AMBER_CRD,
                    PrintSituation expectation = PrintSituation::UNKNOWN, int low_index = 0,
                    int high_index = 0, HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a const reference to the one of the Cartesian coordinate arrays.
  ///
  /// \param dim  The dimension of interest
  const Hybrid<T>& getCoordinateReference(CartesianDimension dim) const;

  /// \brief Get a const pointer to one of the Cartesian coordinate arrays.
  ///
  /// \param dim  The dimension of interest
  const Hybrid<T>* getCoordinatePointer(CartesianDimension dim) const;

  /// \brief Get a const reference to the box space transformation matrices
  const Hybrid<double>& getBoxTransforms() const;

  /// \brief Get a const pointer to the box space transformation matrices
  const Hybrid<double>* getBoxTransformPointer() const;

  /// \brief Get a const reference to the box space transformation matrices
  const Hybrid<double>& getInverseTransforms() const;

  /// \brief Get a const pointer to the box space transformation matrices
  const Hybrid<double>* getInverseTransformPointer() const;

  /// \brief Get a const pointer to the list of box dimensions for all frames
  const Hybrid<double>* getBoxDimensionPointer() const;

  /// \brief Get a const pointer to the object iself, so that it may be passed to functions by
  ///        const reference and still emit a const poiner.
  const CoordinateSeries<T>* getSelfPointer() const;

  /// \brief Get the abstract for this object, containing C-style pointers for the most rapid
  ///        access to any of its member variables.
  ///
  /// Overloaded:
  ///   - Get a read-only abstract from a const coordinate series
  ///   - Get a writeable abstract from a non-const coordinate series
  ///
  /// \param tier  Get pointers to data at the level of the CPU host or GPU device
  /// \{
  const CoordinateSeriesReader<T> data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  CoordinateSeriesWriter<T> data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Get the object's abstract, with all templated type pointers cast to void pointers.
  ///        The true data type can be encoded in a separate 64-bit unsigned integer code (use
  ///        std::type_index() as seen in DataTypes/common_types.h).  Parameter descriptions and
  ///        overloading match those in data() above.
  /// \{
  const CoordinateSeriesReader<void>
  templateFreeData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  CoordinateSeriesWriter<void> templateFreeData(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}
  
#ifdef STORMM_USE_HPC
  /// \brief Get an abstract for the object's host data that is guaranteed to be accessible by the
  ///        GPU device.
  ///
  /// Overloaded:
  ///   - Get a read-only abstract for a const coordinate series
  ///   - Get a writeable abstract for a non-const coordinate series
  /// \{
  const CoordinateSeriesReader<T> deviceViewToHostData() const;
  CoordinateSeriesWriter<T> deviceViewToHostData();
  /// \}

  /// \brief Get an abstract for the object's host data that is guaranteed to be accessible by the
  ///        GPU device with all of its templated pointers cast to void.
  ///
  /// Overloaded:
  ///   - Get a read-only abstract for a const coordinate series
  ///   - Get a writeable abstract for a non-const coordinate series
  /// \{
  const CoordinateSeriesReader<void> deviceViewToTemplateFreeHostData() const;
  CoordinateSeriesWriter<void> deviceViewToTemplateFreeHostData();
  /// \}
  
  /// \brief Upload all information
  void upload();

  /// \brief Download all information
  void download();
#endif
  
  /// \brief Import coordinates from a CoordinateFrame or PhaseSpace object.  The original object
  ///        must have the same number of atoms as the CoordinateSeries itself, or else a range of
  ///        atoms within the original coordinate object must be specified that fits the
  ///        CoordinateSeries.  The default behavior is to push the new coordinates to the back of
  ///        the list, but any frame index within the bounds of the current list may also be
  ///        specified.
  ///
  /// Overloaded:
  ///   - Accept all types of single-frame coordinate objects
  ///   - Accept all atoms or a subset of the atoms that fits the current atom count of the
  ///     CoordinateSeries
  ///
  /// \param cfr          Coordinates to import.  The CoordinateFrameReader is the most basic
  ///                     object available for importation.  Both CoordinateFrame and PhaseSpace
  ///                     objects can create CoordinateFrameReaders or writers, and the writers
  ///                     can be const-ified into readers.
  /// \param cfw          Coordinates to import
  /// \param cf           Coordinates to import
  /// \param ps           Coordinates to import
  /// \param atom_start   First atom from the coordinate set to add to a frame of the series
  /// \param atom_end     Limit of atoms from the coordinate set to add to a frame of the series
  /// \param frame_index  Index of the frame into which the coordinates should be imported.  The
  ///                     default value of -1 adds the new coordinates to the end of the list.
  /// \{
  void import(const CoordinateFrameReader &cfr, int atom_start, int atom_end,
              int frame_index = -1);
  void import(const CoordinateFrameReader &cfr, int frame_index = -1);
  void import(const CoordinateFrameWriter &cfw, int atom_start, int atom_end,
              int frame_index = -1);
  void import(const CoordinateFrameWriter &cfw, int frame_index = -1);
  void import(const CoordinateFrame &cf, int atom_start, int atom_end, int frame_index = -1);
  void import(const CoordinateFrame &cf, int frame_index = -1);
  void import(const CoordinateFrame *cf, int atom_start, int atom_end, int frame_index = -1);
  void import(const CoordinateFrame *cf, int frame_index = -1);
  void import(const PhaseSpace &ps, int atom_start, int atom_end, int frame_index = -1,
              TrajectoryKind kind = TrajectoryKind::POSITIONS,
              CoordinateCycle orientation = CoordinateCycle::WHITE);
  void import(const PhaseSpace &ps, int frame_index = -1,
              TrajectoryKind kind = TrajectoryKind::POSITIONS,
              CoordinateCycle orientation = CoordinateCycle::WHITE);
  void import(const PhaseSpace *ps, int atom_start, int atom_end, int frame_index = -1,
              TrajectoryKind kind = TrajectoryKind::POSITIONS,
              CoordinateCycle orientation = CoordinateCycle::WHITE);
  void import(const PhaseSpace *ps, int frame_index = -1,
              TrajectoryKind kind = TrajectoryKind::POSITIONS,
              CoordinateCycle orientation = CoordinateCycle::WHITE);
  /// \}

  /// \brief Import coordinates from a file.  This function accepts directives to read a subset of
  ///        the coordinates and integrate them into the series at a specified point.
  ///
  /// \param file_name          Name of the file
  /// \param file_kind          The kind of coordinate file being read, default UNKNOWN which will
  ///                           lead to automatic deduction of the file type
  /// \param frame_numbers      List of frame indices (starting from zero) to import from the file
  ///                           (if the file is too small to hold one or more frames with the
  ///                           requested indices, this is an error)
  /// \param replica_count_in   The number of replicas of the frames from this file (default 1).
  ///                           If additional replicas are requested, the whole series will be
  ///                           added in sequence, i.e. frames 3, 5, 8, 9, 3, 5, 8, 9, 3, ...
  /// \param frame_index_start  The starting index of the series for incorporating frames found in
  ///                           the file.  A negative value in this argument triggers addition to
  ///                           the end of any existing series.
  void importFromFile(const std::string &file_name,                       
                      CoordinateFileKind file_kind = CoordinateFileKind::UNKNOWN,
                      const std::vector<int> &frame_numbers = {}, int replica_count = 1,
                      int frame_index_start = -1);

  /// \brief Reserve capacity in this series.  The new frames will be uninitialized.
  ///
  /// \param new_frame_capacity  The new capacity to prepare for.  If such capacity already exists,
  ///                            the function will do nothing.
  void reserve(const int new_frame_capacity);

  /// \brief Shrink the object's various arrays to accommodate only the current frame and atom
  ///        counts.  Mimics the Standard Template Library vector shrink_to_fit() member function.
  void shrinkToFit();

  /// \brief Resize the series, allocating new capacity if needed, initializing new frames with the
  ///        provided coordinate set. ("New" frames are defined as any frames with indices greater
  ///        than the original maximum index, regardless of whether new capacity was allocated to
  ///        hold them or if the size simply increased within the existing space available.)
  ///
  /// Overloaded:
  ///   - Accept all types of single-frame coordinate objects
  ///   - Accept all atoms or a subset of the atoms that fits the current atom count of the
  ///     CoordinateSeries
  ///
  /// \param new_frame_count  The number of new frames that the list shall report holding.  This is
  ///                         less than or equal to the capacity; if there is insufficient
  ///                         capacity when resize() is called, new capacity will be allocated to
  ///                         hold precisely new_frame_count frames.
  /// \param cfr              Optional coordinate set to use in initializing new frames.  If no
  ///                         coordinates are provided, the series will report holding frames but
  ///                         have undefined values in them.
  /// \param cfw              Optional coordinate set to use in initializing new frames
  /// \param cf               Optional coordinate set to use in initializing new frames
  /// \param ps               Optional coordinate set to use in initializing new frames
  /// \param atom_start       First atom from the coordinate set to add to a frame of the series
  ///                         (the default value is zero, to use all atoms)
  /// \param atom_end         Limit of atoms from the coordinate set to add to a frame of the
  ///                         series (the default value of zero will trigger an access to atom
  ///                         count from the coordinate object to load all of its atoms)
  /// \{
  void resize(int new_frame_count);
  void resize(int new_frame_count, const CoordinateFrameReader &cfr, int atom_start = 0,
              int atom_end = 0);
  void resize(int new_frame_count, const CoordinateFrameWriter &cfw, int atom_start = 0,
              int atom_end = 0);
  void resize(int new_frame_count, const CoordinateFrame &cf, int atom_start = 0,
              int atom_end = 0);
  void resize(int new_frame_count, CoordinateFrame *cf, int atom_start = 0, int atom_end = 0);
  void resize(int new_frame_count, const PhaseSpace &ps, int atom_start = 0, int atom_end = 0);
  void resize(int new_frame_count, PhaseSpace *ps, int atom_start = 0, int atom_end = 0);
  /// \}
  
  /// \brief Push a coordinate set to the back of the list.  This invokes the import member
  ///        function after reallocating the frame series with 25% spare capacity if the original
  ///        capacity is insufficient.
  ///
  /// Overloaded:
  ///   - Accept all types of single-frame coordinate objects
  ///   - Accept all atoms or a subset of the atoms that fits the current atom count of the
  ///     CoordinateSeries
  ///
  /// \param cfr          Coordinates to import.  The CoordinateFrameReader is the most basic
  ///                     object available for importation.  Both CoordinateFrame and PhaseSpace
  ///                     objects can create CoordinateFrameReaders or writers, and the writers
  ///                     can be const-ified into readers.
  /// \param cfw          Coordinates to import
  /// \param cf           Coordinates to import
  /// \param ps           Coordinates to import
  /// \param atom_start   First atom from the coordinate set to add to a frame of the series
  /// \param atom_end     Limit of atoms from the coordinate set to add to a frame of the series
  /// \{
  void pushBack(const CoordinateFrameReader &cfr, int atom_start = 0, int atom_end = 0);
  void pushBack(const CoordinateFrameWriter &cfw, int atom_start = 0, int atom_end = 0);
  void pushBack(const CoordinateFrame &cf, int atom_start = 0, int atom_end = 0);
  void pushBack(CoordinateFrame *cf, int atom_start = 0, int atom_end = 0);
  void pushBack(const PhaseSpace &ps, int atom_start = 0, int atom_end = 0);
  void pushBack(PhaseSpace *ps, int atom_start = 0, int atom_end = 0);
  /// \}
  
private:
  int atom_count;                       ///< Number of atoms in each frame.  Between frames the
                                        ///<   space for atoms is padded by the warp size, but
                                        ///<   there is no concept of an atom capacity in the same
                                        ///<   way that there is a frame capacity.
  int frame_count;                      ///< Total number of frames currently in the object
  int frame_capacity;                   ///< Total frame capacity of the object
  int globalpos_scale_bits;             ///< The number of bits for global position scaling, if
                                        ///<   the coordinate series is given in a fixed-precision
                                        ///<   representation
  UnitCellType unit_cell;               ///< The unit cell type for these coordinate frames
  double globalpos_scale;               ///< Scaling factor for converting real-number
                                        ///<   representations of the coordinates into the serie's
                                        ///<   fixed-precision representation, if applicable
  double inverse_globalpos_scale;       ///< Inverse scaling factor for converting a series's
                                        ///<   fixed-precision representation into real numbers,
                                        ///<   if a fixed-precision representation is in effect
  Hybrid<T> x_coordinates;              ///< Cartesian X coordinates of all particles, with each
                                        ///<   frame's coordinates padded by the warp size
  Hybrid<T> y_coordinates;              ///< Cartesian Y coordinates of all particles 
  Hybrid<T> z_coordinates;              ///< Cartesian Z coordinates of all particles
  Hybrid<double> box_space_transforms;  ///< Matrices to transform each frame into box space
  Hybrid<double> inverse_transforms;    ///< Matrix to transform each frame into real space
  Hybrid<double> box_dimensions;        ///< Lengths and angles defining each frame's unit cell
                                        ///<   (lengths are given in Angstroms, angles in radians)

  /// \brief Allocate space for this series.  This will only allocate in the forward direction,
  ///        never less than already exists, but it does so stepwise through each of the x, y, and
  ///        z coordinate arrays and each of the box transformation arrays.  This process will
  ///        reduce the spike in memory needed to allocated a larger array ands transfer over the
  ///        original data.
  ///
  /// \param new_frame_capacity  The new capacity to allocate.  The frame count and atom count are
  ///                            are both int type, although the total size is measured in size_t
  ///                            to avoid integer overflow.
  void allocate(int new_frame_capacity);
};

/// \brief Reconstruct an abstract of the Coordinate series with a specific templated type.  This
///        free function complements the templateFreeData() member function in the CoordinateSeries
///        object and can be used to complete the handoff of information from C++-compiled
///        libraries.  The object as well as this associated free function will be compiled by
///        both C++ and HPC compilers so that, on each side of the fence, there is a means to cast
///        the templated coordinate object to <void> and then unpack the contents in the other
///        side.
///
/// Overloaded:
///   - Produce a writeable abstract based on a non-const CoordinateSeriesWriter
///   - Produce a read-only abstract based on a CoordinateSeriesReader
///   - Accept the original, void-cast abstract by pointer or by reference
///
/// \param rasa  Abstract of the original coordinate series, free of templated pointers (inspired
///              by the phrase "tabula rasa", a blank slate)
/// \{
template <typename T> CoordinateSeriesWriter<T> restoreType(CoordinateSeriesWriter<void> *rasa);

template <typename T>
CoordinateSeriesWriter<T> restoreType(const CoordinateSeriesWriter<void> &rasa);

template <typename T>
const CoordinateSeriesReader<T> restoreType(const CoordinateSeriesReader<void> *rasa);

template <typename T>
const CoordinateSeriesReader<T> restoreType(const CoordinateSeriesReader<void> &rasa);
/// \}

/// \brief Convert a coordinate series from one data type to another.  There are three major data
///        types for a CoordinateSeries: double, float, and long long int.  This function has
///        numerous branches to make efficient conversions between real and signed integral data
///        and minimize any loss of precision in the process.
///
/// \param cs  The original coordinate series
template <typename Torig, typename Tnew>
CoordinateSeries<Tnew> changeCoordinateSeriesType(const CoordinateSeries<Torig> &cs,
                                                  int globalpos_scale_bits_in);

} // namespace trajectory 
} // namespace stormm

#include "coordinate_series.tpp"

#endif
