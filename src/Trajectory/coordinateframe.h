// -*-c++-*-
#ifndef STORMM_COORDINATE_FRAME_H
#define STORMM_COORDINATE_FRAME_H

#include <typeinfo>
#include <typeindex>
#include <sys/types.h>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Parsing/textfile.h"
#include "Topology/atomgraph.h"
#include "phasespace.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

using card::Hybrid;
using card::HybridKind;
using card::HybridTargetLevel;
using constants::ExceptionResponse;
using parse::TextFile;
using topology::AtomGraph;

/// \brief Collect C-style pointers for the elements of a writable CoordinateFrame object.
struct CoordinateFrameWriter {

  /// \brief The constructor feeds all arguments straight to the inline initialization list.
  ///
  /// Overloaded:
  ///   - Take all arguments piecemeal
  ///   - Take a PhaseSpace object
  /// \{
  CoordinateFrameWriter(int natom_in, UnitCellType unit_cell_in, double* xcrd_in, double* ycrd_in,
                        double* zcrd_in, double* umat_in, double* invu_in, double* boxdim_in);
  CoordinateFrameWriter(PhaseSpace *ps, HybridTargetLevel tier = HybridTargetLevel::HOST);
  CoordinateFrameWriter(PhaseSpace *ps, TrajectoryKind kind,
                        HybridTargetLevel tier = HybridTargetLevel::HOST);
  CoordinateFrameWriter(PhaseSpace *ps, TrajectoryKind kind, CoordinateCycle orientation,
                        HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Copy and move constructors.  The move assignment operator is implicitly deleted.
  /// \{
  CoordinateFrameWriter(const CoordinateFrameWriter &original) = default;
  CoordinateFrameWriter(CoordinateFrameWriter &&original) = default;
  /// \}
  
  const int natom;               ///< The number of atoms in the system
  const UnitCellType unit_cell;  ///< The type of unit cell (i.e. ORTHORHOMBIC, could also be NONE)
  double* xcrd;                  ///< Cartesian X coordinates of all atoms
  double* ycrd;                  ///< Cartesian Y coordinates of all atoms
  double* zcrd;                  ///< Cartesian Z coordinates of all atoms
  double* umat;                  ///< Transformation matrix to take coordinates into box
                                 ///<   (fractional) space
  double* invu;                  ///< Inverse transformation matrix out of box space
  double* boxdim;                ///< Box dimensions (these will be consistent with umat and invu)
};

/// \brief Collect C-style pointers for the elements of a read-only CoordinateFrame object.
struct CoordinateFrameReader {

  /// \brief The constructor feeds all arguments straight to the inline initialization list.
  ///
  /// Overloaded:
  ///   - Take all arguments piecemeal
  ///   - Take a CoordinateFrameWriter
  ///   - Take a PhaseSpace object
  /// \{
  CoordinateFrameReader(int natom_in, UnitCellType unit_cell_in, const double* xcrd_in,
                        const double* ycrd_in, const double* zcrd_in, const double* umat_in,
                        const double* invu_in, const double* boxdim_in);
  CoordinateFrameReader(const CoordinateFrameWriter &cfw);
  CoordinateFrameReader(const PhaseSpace &ps, HybridTargetLevel tier = HybridTargetLevel::HOST);
  CoordinateFrameReader(const PhaseSpace *ps, HybridTargetLevel tier = HybridTargetLevel::HOST);
  CoordinateFrameReader(const PhaseSpace *ps, TrajectoryKind kind, CoordinateCycle orientation,
                        HybridTargetLevel tier = HybridTargetLevel::HOST);
  CoordinateFrameReader(const PhaseSpace &ps, TrajectoryKind kind, CoordinateCycle orientation,
                        HybridTargetLevel tier = HybridTargetLevel::HOST);
  CoordinateFrameReader(const PhaseSpace *ps, TrajectoryKind kind,
                        HybridTargetLevel tier = HybridTargetLevel::HOST);
  CoordinateFrameReader(const PhaseSpace &ps, TrajectoryKind kind,
                        HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Copy and move constructors.  The move assignment operator is implicitly deleted.
  /// \{
  CoordinateFrameReader(const CoordinateFrameReader &original) = default;
  CoordinateFrameReader(CoordinateFrameReader &&original) = default;
  /// \}
  
  const int natom;               ///< The number of atoms in the system
  const UnitCellType unit_cell;  ///< The type of unit cell (i.e. ORTHORHOMBIC, could also be NONE)
  const double* xcrd;            ///< Cartesian X coordinates of all atoms
  const double* ycrd;            ///< Cartesian Y coordinates of all atoms
  const double* zcrd;            ///< Cartesian Z coordinates of all atoms
  const double* umat;            ///< Transformation matrix to take coordinates into box
                                 ///<   (fractional) space
  const double* invu;            ///< Inverse transformation matrix out of box space
  const double* boxdim;          ///< Box dimensions (these will be consistent with umat and invu)
};

/// \brief Store the coordinates and box information for a frame, only.  This abridged struct can
///        serve when the full PhaseSpace object would allocate too much memory.  It also comes
///        with its own POINTER mode, such that it allocates no memory of its own and merely points
///        to another CoordinateFrame object or PhaseSpace object that does have memory allocated.
class CoordinateFrame {
public:

  /// \brief There are several options for construction of this abridged, coordinate-only object.
  ///
  /// Overloaded:
  ///   - Allocate to hold a given number of atoms
  ///   - From an atom count and a set of double-precision C-style arrays, including coordinates
  ///     and box dimensions
  ///   - Create from any of the coordinate file formats (velocities will not be read, even if
  ///     they are available)
  ///   - From an existing PhaseSpace object (as a pointer or a copy if the PhaseSpace object is
  ///     non-const, otherwise only as a copy)
  ///
  /// \param  natom_in         The number of atoms expected
  /// \param  xcrd_in          Cartesian X coordinates of all particles
  /// \param  ycrd_in          Cartesian Y coordinates of all particles
  /// \param  zcrd_in          Cartesian Z coordinates of all particles
  /// \param  umat_in          Matrix to transform coordinates into box space (3 x 3)
  /// \param  invu_in          Matrix to transform coordinates into real space (3 x 3)
  /// \param  file_name_in     File to read from
  /// \param  file_kind        The type of coordinate file to expect
  /// \param  frame_number_in  Frame number of the file to read (default 0, the first frame)
  /// \param  ps               Pre-existing object with complete description of the system state
  /// \{
  CoordinateFrame(int natom_in = 0, UnitCellType unit_cell_in = UnitCellType::NONE);
  CoordinateFrame(int natom_in, UnitCellType unit_cell_in, const double* xcrd_in,
                  const double* ycrd_in, const double* zcrd_in, const double* umat_in = nullptr,
                  const double* invu_in = nullptr, const double* boxdim_in = nullptr);
  CoordinateFrame(const std::string &file_name_in,
                  CoordinateFileKind file_kind = CoordinateFileKind::UNKNOWN,
                  int frame_number_in = 0);
  CoordinateFrame(const TextFile &tf, CoordinateFileKind file_kind  = CoordinateFileKind::UNKNOWN,
                  int frame_number_in = 0);
  CoordinateFrame(PhaseSpace *ps);
  CoordinateFrame(const PhaseSpace &ps);
  /// \}

  /// \brief Copy constructor handles assignment of internal POINTER-kind Hybrid objects
  ///
  /// \param original  The PhaseSpace object from which to make a deep copy
  CoordinateFrame(const CoordinateFrame &original);

  /// \brief Copy assignment operator likewise handles reassignment of internal POINTER-kind
  ///        Hybrid objects
  ///
  /// \param other     Another way to say original, in a different semantic context
  CoordinateFrame& operator=(const CoordinateFrame &other);

  /// \brief Move constructor works in similar fashion to the PhaseSpace move constructor,
  ///        with std::move and no need for reassignment of the underlying POINTER-kind Hybrid
  ///        objects.
  ///
  /// \param original  The PhaseSpace object from which to make a deep copy
  CoordinateFrame(CoordinateFrame &&original);

  /// \brief Move assignment operator
  ///
  /// \param other     Another way to say original, in a different semantic context
  CoordinateFrame& operator=(CoordinateFrame &&other);

  /// \brief Fill the object from information some coordinate, restart, or trajectory file.
  ///
  /// Overloaded:
  ///   - Take a file name
  ///   - Take an ASCII-format text file already processed into RAM
  ///
  /// \param file_name     Name of the file from which to obtain coordinates
  /// \param file_kind     The type of coordinate-containing input file
  /// \param frame_number  The frame number to read (if the file is a trajectory, not a single
  ///                      point from the system's phase space)
  /// \param tf            Text file data already processed into RAM
  /// \{
  void buildFromFile(const std::string &file_name_in, const CoordinateFileKind file_kind,
                     int frame_number = 0);
  void buildFromFile(const TextFile &tf, const CoordinateFileKind file_kind,
                     int frame_number = 0);
  /// \}

  /// \brief Fill the object from information in three arrays.
  ///
  /// Overloaded:
  ///   - Fill from three C-style arrays
  ///   - Fill from three Standard Template Library vector objects
  ///
  /// \param xcrd        Cartesian X coordinates of positions, velocities, or forces
  /// \param ycrd        Cartesian Y coordinates of positions, velocities, or forces
  /// \param zcrd        Cartesian Z coordinates of positions, velocities, or forces
  /// \param kind        Type of coordinates coming in: fill the positions, velocities, or forces
  /// \param cycle_in    The point in the coordinate cycle to fill
  /// \param scale_bits  The number of bits after the decimal, applicable to fixed-precision
  ///                    representations of xcrd, ycrd, and zcrd (the box dimensions are always
  ///                    given as a double-precision array, in units of Angstroms)
  /// \param box_dims    Box dimensions, from which the tranformation matrices will be derived
  /// \{
  template <typename T>
  void fill(const T* xcrd, const T* ycrd, const T* zcrd, int scale_bits = 0,
            const double* box_dims = nullptr);

  template <typename T>
  void fill(const std::vector<T> &xcrd, const std::vector<T> &ycrd, const std::vector<T> &zcrd,
            int scale_bits = 0, const std::vector<double> &box_dims = {});
  /// \}

  /// \brief Get the file name that originated this coordinate set
  std::string getFileName() const;

  /// \brief Get the frame number, if this coordinate set originated in a trajectory.
  int getFrameNumber() const;
  
  /// \brief Get the number of atoms in the frame
  int getAtomCount() const;

  /// \brief Get the unit cell type of the coordinate system
  UnitCellType getUnitCellType() const;

  /// \brief Get the coordinates returned in an X/Y/Z interlaced manner
  ///
  /// Overloaded:
  ///   - Get all coordinates
  ///   - Get coordinates for a range of atoms
  ///
  /// \param low_index   The lower atom index of a range
  /// \param high_index  The upper atom index of a range
  /// \param kind        Specify coordinates, velocities, or forces--anything that could be thought
  ///                    of as a trajectory
  /// \param tier        Level at which to extract the data
  /// \{
  std::vector<double>
  getInterlacedCoordinates(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  std::vector<double>
  getInterlacedCoordinates(int low_index, int high_index,
                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Get the transformation matrix to take coordinates into box (fractional) space.  The
  ///        result should be interpreted as a 3x3 matrix in column-major format.
  ///
  /// \param tier  Level at which to retrieve the data (if STORMM is compiled to run on a GPU)
  std::vector<double> getBoxSpaceTransform(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the transformation matrix to take coordinates from fractional space back into
  ///        real space.  The result should be interpreted as a 3x3 matrix in column-major format.
  ///
  /// \param tier  Level at which to retrieve the data (if STORMM is compiled to run on a GPU)
  std::vector<double> getInverseTransform(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the box dimensions in their pure form.  Holding the box dimensions and the
  ///        transformation matrices they imply in separate member variables represents a
  ///        liability, that they might at some point become inconsistent, but it also prevents
  ///        having to continuously extract box dimensions from transformation matrices each time
  ///        the box needs to be resized (beyond isotropic scaling, this is not a trivial process).
  ///
  /// \param tier  Level at which to retrieve the data (if STORMM is compiled to run on a GPU)  
  std::vector<double> getBoxDimensions(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Export the contents of this coordinate series to a trajectory, restart, or
  ///        input coordinates file.
  ///
  /// \param file_name     Name of the file to write
  /// \param output_kind   The format of the file to write (checkpoint files print position and
  ///                      velocity data by obligation, but trajectory files can contain either of
  ///                      these as well as forces)
  /// \param expectation   The condition in which the output file is expected to be found
  void exportToFile(const std::string &file_name,
                    CoordinateFileKind output_kind = CoordinateFileKind::AMBER_CRD,
                    PrintSituation expectation = PrintSituation::UNKNOWN) const;

  /// \brief Get a pointer to the object itself, useful when working with a const reference.
  const CoordinateFrame* getSelfPointer() const;
  
  /// \brief Get the abstract for this object, containing C-style pointers for the most rapid
  ///        access to any of its member variables.
  ///
  /// Overloaded:
  ///   - Get a read-only abstract from a const object
  ///   - Get a writeable abstract from a mutable object
  /// \{
  const CoordinateFrameReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  CoordinateFrameWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

#ifdef STORMM_USE_HPC
  /// \brief Upload all information
  void upload();

  /// \brief Download all information
  void download();

  /// \brief Get a read-only abstract that is guaranteed to be valid on the device.  While many
  ///        GPU devices will accept the output of data() and the level of the HOST, the pointers
  ///        set herein will work for all devices.
  ///
  /// Overloaded:
  ///   - Get a read-only abstract from a const object
  ///   - Get a writeable abstract from a mutable object
  /// \{
  const CoordinateFrameReader deviceViewToHostData() const;
  CoordinateFrameWriter deviceViewToHostData();
  /// \}
#endif

  /// \brief Set the frame number, for bookkeeping purposes.  This function exists that the frame
  ///        number does not become a necessary argument in one of the overloaded constructors.
  ///
  /// \param frame_number_in  The (arbitrary) number to set the frame as
  void setFrameNumber(int frame_number_in);
  
private:
  std::string file_name;              ///< File from which this frame was derived, if applicable.
  int frame_number;                   ///< The frame number at which this frame appears in the
                                      ///<   trajectory file, if applicable.  Indexing begins at 0.
  int atom_count;                     ///< The number of atoms in the system
  UnitCellType unit_cell;             ///< The type of unit cell (can be inferred from the
                                      ///<   transformation matrices)
  Hybrid<double> x_coordinates;       ///< Cartesian X coordinates of all particles
  Hybrid<double> y_coordinates;       ///< Cartesian Y coordinates of all particles
  Hybrid<double> z_coordinates;       ///< Cartesian Z coordinates of all particles
  Hybrid<double> box_space_transform; ///< Matrix to transform coordinates into box space (3 x 3)
  Hybrid<double> inverse_transform;   ///< Matrix to transform coordinates into real space (3 x 3)
  Hybrid<double> box_dimensions;      ///< Three lengths and three angles defining the box (lengths
                                      ///<   are given in Angstroms, angles in radians)

  /// All of the above Hybrid objects can either be pointers into this single large array, which
  /// will then be segmented to hold each type of information with zero-padding to accommodate the
  /// HPC warp size, or pointers into a PhaseSpace object's own data array, in which case this
  /// array is likewise a pointer into the same PhaseSpace object's data array.  If this is itself
  /// a POINTER-kind Hybrid object, the upload and download functions for this CoordinateFrame
  /// will act on only a subset of the data in the target PhaseSpace object.  Modifications of the
  /// coordinates or box dimensions are possible, and would be reflected in any target PhaseSpace
  /// object, so uploading and downloading is permitted.
  Hybrid<double> storage;

  /// \brief Allocate memory and set POINTER-kind Hybrid object member variables for this
  ///        coordinate frame.  Follows the PhaseSpace object's eponymous member function.
  void allocate();
};

/// \brief Define the type index for the CoordinateFrame object.
static const size_t coordinateframe_type_index =
  std::type_index(typeid(CoordinateFrame)).hash_code();
  
/// \brief Read a series of frames from a trajectory file.
///
/// Overloaded:
///   - Read multiple frames from an ASCII-format trajectory (Amber .crd, or a more efficient way
///     to read multiple copies of a single Amber inpcrd or restart file)
///   - Read multiple frames from a binary trajectory file
///
/// \param tf             ASCII text file pre-loaded into memory
/// \param kind           Type of coordinate or restart file
/// \param atom_count     Expected atom count (will be checked against the file, if possible)
/// \param unit_cell      Expected unit cell type (will be written to each frame)
/// \param frame_numbers  Indices of frames to read from the file
/// \{
std::vector<CoordinateFrame> getSelectedFrames(const TextFile &tf, CoordinateFileKind kind,
                                               int atom_count, UnitCellType unit_cell,
                                               const std::vector<int> &frame_numbers);

std::vector<CoordinateFrame> getSelectedFrames(const std::string &file_name, int atom_count,
                                               UnitCellType unit_cell,
                                               const std::vector<int> &frame_numbers);
/// \}

/// \brief Read all available frames from a trajectory file.
///
/// Overloaded:
///   - Read all frames of an ASCII-format trajectory (the frame count will be predicted by the
///     number of lines in the text file and the expected format)
///   - Read all frames from a binary trajectory file
///
/// \param tf             ASCII text file pre-loaded into memory
/// \param kind           Type of coordinate or restart file
/// \param atom_count     Expected atom count (will be checked against the file, if possible)
/// \param unit_cell      Expected unit cell type (will be written to each frame)
/// \{
std::vector<CoordinateFrame> getAllFrames(const TextFile &tf, const int atom_count,
                                          const UnitCellType unit_cell,
                                          const ExceptionResponse policy);

std::vector<CoordinateFrame> getAllFrames(const std::string &file_name, const int atom_count,
                                          const UnitCellType unit_cell,
                                          const ExceptionResponse policy);
/// \}

} // namespace trajectory
} // namespace stormm

#include "coordinateframe.tpp"

// As with common types and STORMM vector types, define the type indices for general use in the
// STORMM namespace.
namespace stormm {
using trajectory::coordinateframe_type_index;
} // namespace stormm

#endif
