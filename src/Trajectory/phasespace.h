// -*-c++-*-
#ifndef STORMM_PHASE_SPACE_H
#define STORMM_PHASE_SPACE_H

#include <typeinfo>
#include <typeindex>
#include <sys/types.h>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/hybrid_util.h"
#include "Constants/behavior.h"
#include "Math/matrix_ops.h"
#include "Topology/atomgraph.h"
#include "trajectory_enumerators.h"
#include "write_frame.h"

namespace stormm {
namespace trajectory {

using card::default_hpc_format;
using card::Hybrid;
using card::HybridFormat;
using card::HybridTargetLevel;
using constants::CartesianDimension;
using stmath::computeBoxTransform;
using topology::AtomGraph;
using topology::UnitCellType;

/// \brief Collect constants and pointers to the components of a modifiable PhaseSpace object.
struct PhaseSpaceWriter {

  /// Constructor takes a large list of arguments passed in from the original PhaseSpace object
  PhaseSpaceWriter(const int natom_in, const UnitCellType unit_cell_in, double* xcrd_in,
                   double* ycrd_in, double* zcrd_in, double* umat_in, double* invu_in,
                   double* boxdim_in, double* umat_alt_in, double* invu_alt_in,
                   double* boxdim_alt_in, double* xvel_in, double* yvel_in, double* zvel_in,
                   double* xfrc_in, double* yfrc_in, double* zfrc_in, double* xalt_in,
                   double* yalt_in, double* zalt_in, double* vxalt_in, double* vyalt_in,
                   double* vzalt_in, double* fxalt_in, double* fyalt_in, double* fzalt_in);

  /// \brief Copy and move constructors.  The assignment operators are implicitly deleted.
  /// \{
  PhaseSpaceWriter(const PhaseSpaceWriter &original) = default;
  PhaseSpaceWriter(PhaseSpaceWriter &&original) = default;
  /// \}

  const int natom;               ///< Atom count for this system (still a constant)
  const UnitCellType unit_cell;  ///< The type of unit cell
  double* xcrd;                  ///< Cartesian X positions of all particles
  double* ycrd;                  ///< Cartesian Y positions of all particles
  double* zcrd;                  ///< Cartesian Z positions of all particles
  double* umat;                  ///< Transformation matrix to take coordinates into box
                                 ///<   (fractional) space
  double* invu;                  ///< Transformation matrix to take coordinates into real space
  double* boxdim;                ///< Box dimensions (stored for convenience and accuracy if the
                                 ///<   box is resized repeatedly)
  double* umat_alt;              ///< Transformation matrix to take coordinates into box
                                 ///<   (fractional) space
  double* invu_alt;              ///< Transformation matrix to take coordinates into real space
  double* boxdim_alt;            ///< Box dimensions (stored for convenience and accuracy if the
                                 ///<   box is resized repeatedly)
  double* xvel;                  ///< Cartesian X velocities of all particles
  double* yvel;                  ///< Cartesian Y velocities of all particles
  double* zvel;                  ///< Cartesian Z velocities of all particles
  double* xfrc;                  ///< Cartesian X forces acting on all particles
  double* yfrc;                  ///< Cartesian Y forces acting on all particles
  double* zfrc;                  ///< Cartesian Z forces acting on all particles
  double* xalt;                  ///< Alternate Cartesian X positions of all particles
  double* yalt;                  ///< Alternate Cartesian Y positions of all particles
  double* zalt;                  ///< Alternate Cartesian Z positions of all particles
  double* vxalt;                 ///< Alternate Cartesian X velocities for all particles
  double* vyalt;                 ///< Alternate Cartesian Y velocities for all particles
  double* vzalt;                 ///< Alternate Cartesian Z velocities for all particles
  double* fxalt;                 ///< Alternate Cartesian X forces acting on all particles
  double* fyalt;                 ///< Alternate Cartesian Y forces acting on all particles
  double* fzalt;                 ///< Alternate Cartesian Z forces acting on all particles
};

/// \brief Collect constants and pointers to the components of a read-only PhaseSpace object.
struct PhaseSpaceReader {

  /// The constructor takes a large list of arguments passed in from the original PhaseSpace
  /// object, or the cognate writer to make all of the associated pointers const.
  /// \{
  PhaseSpaceReader(const int natom_in, const UnitCellType unit_cell_in, const double* xcrd_in,
                   const double* ycrd_in, const double* zcrd_in, const double* umat_in,
                   const double* invu_in, const double* boxdim_in, const double* umat_alt_in,
                   const double* invu_alt_in, const double* boxdim_alt_in, const double* xvel_in,
                   const double* yvel_in, const double* zvel_in, const double* xfrc_in,
                   const double* yfrc_in, const double* zfrc_in, const double* xalt_in,
                   const double* yalt_in, const double* zalt_in, const double* vxalt_in,
                   const double* vyalt_in, const double* vzalt_in, const double* fxalt_in,
                   const double* fyalt_in, const double* fzalt_in);

  PhaseSpaceReader(const PhaseSpaceWriter &psw);
  /// \}
  
  /// \brief Copy and move constructors.  The assignment operators are implicitly deleted.
  /// \{
  PhaseSpaceReader(const PhaseSpaceReader &original) = default;
  PhaseSpaceReader(PhaseSpaceReader &&original) = default;
  /// \}
  
  const int natom;               ///< Atom count for this system (still a constant)
  const UnitCellType unit_cell;  ///< The type of unit cell
  const double* xcrd;            ///< Cartesian X positions of all particles
  const double* ycrd;            ///< Cartesian Y positions of all particles
  const double* zcrd;            ///< Cartesian Z positions of all particles
  const double* umat;            ///< Transformation matrix to take coordinates into box
                                 ///<   (fractional) space
  const double* invu;            ///< Transformation matrix to take coordinates into real space
  const double* boxdim;          ///< Box dimensions (stored for convenience and accuracy if the
                                 ///<   box is resized repeatedly)
  const double* umat_alt;        ///< Transformation matrix to take coordinates into box
                                 ///<   (fractional) space
  const double* invu_alt;        ///< Transformation matrix to take coordinates into real space
  const double* boxdim_alt;      ///< Box dimensions (stored for convenience and accuracy if the
                                 ///<   box is resized repeatedly)
  const double* xvel;            ///< Cartesian X velocities of all particles
  const double* yvel;            ///< Cartesian Y velocities of all particles
  const double* zvel;            ///< Cartesian Z velocities of all particles
  const double* xfrc;            ///< Cartesian X forces acting on all particles
  const double* yfrc;            ///< Cartesian Y forces acting on all particles
  const double* zfrc;            ///< Cartesian Z forces acting on all particles
  const double* xalt;            ///< Alternate Cartesian X positions of all particles
  const double* yalt;            ///< Alternate Cartesian Y positions of all particles
  const double* zalt;            ///< Alternate Cartesian Z positions of all particles
  const double* vxalt;           ///< Alternate Cartesian X velocities for all particles
  const double* vyalt;           ///< Alternate Cartesian Y velocities for all particles
  const double* vzalt;           ///< Alternate Cartesian Z velocities for all particles
  const double* fxalt;           ///< Alternate Cartesian X forces acting on all particles
  const double* fyalt;           ///< Alternate Cartesian Y forces acting on all particles
  const double* fzalt;           ///< Alternate Cartesian Z forces acting on all particles
};

/// \brief An object to complement a topology and hold positions, velocities, and forces of all
///        particles in a system.  This is not designed to be the most performant representation of
///        the system's structure.  Rather, it serves to hold a high-precision representation of a
///        single system and transport it between CPUs and high-performance accelerators.
class PhaseSpace {
public:

  /// \brief Construction of a phase space object, like a topology, is typically done from a file.
  ///
  /// Overloaded:
  ///   - Constructor fot an object with a number of atoms but no other information (assists in
  ///     delegation of initialization)
  ///   - Constructors for coordinate sets read from trajectory or restart files, with the option
  ///     of a specific frame number (if unspecified, the first frame is read).  A topology may
  ///     also be specified to check the atom count and sanity of the coordinates presented.
  ///
  /// \param file_name     Name of the file from which to obtain coordinates
  /// \param file_kind     The type of coordinate-containing input file
  /// \param ag_reference  Topology with which to check atom count, proximity of bonded atoms, and
  ///                      sanity of bond angles (optional)
  /// \param caller        Name of the calling function (to help backtrace errors)
  /// \{
  PhaseSpace(int atom_count_in = 0, UnitCellType unit_cell_in = UnitCellType::NONE,
             HybridFormat format_in = default_hpc_format);

  PhaseSpace(const std::string &file_name_in,
             CoordinateFileKind file_kind = CoordinateFileKind::UNKNOWN, int frame_number = 0,
             HybridFormat format_in = default_hpc_format);

  PhaseSpace(const std::string &file_name_in, const AtomGraph &ag,
             CoordinateFileKind file_kind = CoordinateFileKind::UNKNOWN,
             int frame_number = 0, HybridFormat format_in = default_hpc_format);
  /// \}

  /// \brief The copy constructor handles assignment of internal POINTER-kind Hybrid objects.
  ///
  /// Overloaded:
  ///   - Take the original object's memory layout
  ///   - Apply an alternate memory layout.  The content will be determined by whatever content is
  ///     available in the original object at each tier of memory: if the new object is to have a
  ///     memory component on the GPU device and the original object also has a memory component on
  ///     the GPU, then this is the state that will be copied over.  Otherwise, the original
  ///     object's data from the CPU host will become the new object's GPU component.  This
  ///     priority of "copy data from what exists at the same level, otherwise take from the other
  ///     level" applies everywhere.
  ///
  /// \param original  The PhaseSpace object from which to make a deep copy
  /// \param format_in  An alternate memory format in which to lay out the new CoordinateFrame
  /// \{
  PhaseSpace(const PhaseSpace &original);
  PhaseSpace(const PhaseSpace &original, HybridFormat format_in);
  /// \}

  /// \brief Copy assignment operator likewise handles assignment of internal POINTER-kind Hybrid
  ///        objects
  ///  
  /// \param other     Another way to say original, in a different semantic context
  PhaseSpace& operator=(const PhaseSpace &other);

  /// \brief The move constructor prepares the original PhaseSpace object for destruction
  ///
  /// \param original  The PhaseSpace object from which to preserve content
  PhaseSpace(PhaseSpace &&original);

  /// \brief The move assignment operator looks much like the copy assignment operator.
  ///  
  /// \param other     Another way to say original, in a different semantic context
  PhaseSpace& operator=(PhaseSpace &&other);

  /// \brief Fill the object from information in some coordinate, restart, or trajectory file.
  ///
  /// \param file_name     Name of the file from which to obtain coordinates
  /// \param file_kind     The type of coordinate-containing input file
  /// \param frame_number  The frame number to read (if the file is a trajectory, not a single
  ///                      point from the system's phase space)
  void buildFromFile(const std::string &file_name_in,
                     const CoordinateFileKind file_kind = CoordinateFileKind::UNKNOWN,
                     int frame_number = 0);

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
  void fill(const T* xcrd, const T* ycrd, const T* zcrd,
            TrajectoryKind kind = TrajectoryKind::POSITIONS,
            CoordinateCycle cycle_in = CoordinateCycle::WHITE, int scale_bits = 0,
            const double* box_dims = nullptr);

  template <typename T>
  void fill(const std::vector<T> &xcrd, const std::vector<T> &ycrd, const std::vector<T> &zcrd,
            TrajectoryKind kind = TrajectoryKind::POSITIONS,
            CoordinateCycle cycle_in = CoordinateCycle::WHITE, int scale_bits = 0,
            const std::vector<double> &box_dims = {});
  /// \}

  /// \brief Get the Hybrid format taken by the object, indicating on which resources its memory
  ///        is present.
  HybridFormat getFormat() const;
  
  /// \brief Get the name of the file associated with this object.
  std::string getFileName() const;

  /// \brief Get the number of atoms (particles, including virtual sites) in the object.
  int getAtomCount() const;

  /// \brief Get the unit cell type of the coordinate system
  UnitCellType getUnitCellType() const;

  /// \brief Get the time cycle stage indicating the arrays holding current coordinates.
  CoordinateCycle getCyclePosition() const;  
  
  /// \brief Get a pointer to the particle X, Y, or Z coordinates, velocities, or forces, on either
  ///        the host or device.  Use this when the entire abstract is unnecessary or would be
  ///        inefficient to retrieve.
  ///
  /// Overloaded:
  ///   - Get a const pointer for a const PhaseSpace object's data
  ///   - Get a non-const pointer to a mutable PhaseSpace object's data
  ///
  /// \param dim   Cartesian dimension of interest
  /// \param kind  Specify coordinates, velocities, or forces--anything that could be thought
  ///              of as a trajectory
  /// \param tier  Level at which to extract the data
  /// \{
  const double* getCoordinatePointer(CartesianDimension dim,
                                     TrajectoryKind kind = TrajectoryKind::POSITIONS,
                                     HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double* getCoordinatePointer(CartesianDimension dim,
                               TrajectoryKind kind = TrajectoryKind::POSITIONS,
                               HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}
  
  /// \brief Get the coordinates returned in an X/Y/Z interlaced manner
  ///
  /// Overloaded:
  ///   - Get all coordinates
  ///   - Get coordinates for a range of atoms
  ///   - Choose the point in the time cycle or accept the object's "current" state
  ///
  /// \param low_index   The lower atom index of a range
  /// \param high_index  The upper atom index of a range
  /// \param kind        Specify coordinates, velocities, or forces--anything that could be thought
  ///                    of as a trajectory
  /// \param tier        Level at which to extract the data
  /// \{
  std::vector<double>
  getInterlacedCoordinates(TrajectoryKind kind = TrajectoryKind::POSITIONS,
                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  std::vector<double>
  getInterlacedCoordinates(CoordinateCycle orientation,
                           TrajectoryKind kind = TrajectoryKind::POSITIONS,
                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  std::vector<double>
  getInterlacedCoordinates(int low_index, int high_index,
                           TrajectoryKind kind = TrajectoryKind::POSITIONS,
                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  std::vector<double>
  getInterlacedCoordinates(int low_index, int high_index, CoordinateCycle orientation,
                           TrajectoryKind kind = TrajectoryKind::POSITIONS,
                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Get a pointer to the box space transformation matrix that can track its evolution in
  ///        the PhaseSpace object's Hybrid data arrays, on either the host or device.
  ///
  /// Overloaded:
  ///   - Get a const pointer for a const PhaseSpace object's data
  ///   - Get a non-const pointer to a mutable PhaseSpace object's data
  ///   - Indicate a specific point in the time cycle, or accept whatever the object sees as its
  ///     own current stage of the cycle
  ///
  /// \param orientation  A specific point in the time cycle to query
  /// \param tier         Level at which to extract the data
  /// \{
  const double*
  getBoxSpaceTransformPointer(CoordinateCycle orientation,
                              HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double* getBoxSpaceTransformPointer(CoordinateCycle orientation,
                                      HybridTargetLevel tier = HybridTargetLevel::HOST);
  const double*
  getBoxSpaceTransformPointer(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double* getBoxSpaceTransformPointer(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Get a pointer to the inverse (back to real space) transformation matrix that can track
  ///        its evolution in the PhaseSpace object's Hybrid data arrays, on either the host or
  ///        device.
  ///
  /// Overloaded:
  ///   - Get a const pointer for a const PhaseSpace object's data
  ///   - Get a non-const pointer to a mutable PhaseSpace object's data
  ///   - Indicate a specific point in the time cycle, or accept whatever the object sees as its
  ///     own current stage of the cycle
  ///
  /// \param orientation  A specific point in the time cycle to query
  /// \param tier         Level at which to extract the data
  /// \{
  const double* getInverseTransformPointer(CoordinateCycle orientation,
                                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double* getInverseTransformPointer(CoordinateCycle orientation,
                                     HybridTargetLevel tier = HybridTargetLevel::HOST);
  const double* getInverseTransformPointer(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double* getInverseTransformPointer(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Get a pointer to the box dimensions that can track its evolution in the PhaseSpace
  ///        object's Hybrid data arrays, on either the host or device.
  ///
  /// Overloaded:
  ///   - Get a const pointer for a const PhaseSpace object's data
  ///   - Get a non-const pointer to a mutable PhaseSpace object's data
  ///   - Indicate a specific point in the time cycle, or accept whatever the object sees as its
  ///     own current stage of the cycle
  ///
  /// \param orientation  A specific point in the time cycle to query
  /// \param tier         Level at which to extract the data
  /// \{
  const double* getBoxSizePointer(CoordinateCycle orientation,
                                  HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double* getBoxSizePointer(CoordinateCycle orientation,
                            HybridTargetLevel tier = HybridTargetLevel::HOST);
  const double* getBoxSizePointer(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double* getBoxSizePointer(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Get the transformation matrix to take coordinates into box (fractional) space.  The
  ///        result should be interpreted as a 3x3 matrix in column-major format.
  ///
  /// Overloaded:
  ///   - Indicate a specific point in the time cycle
  ///   - Accept whatever the object sees as its own current stage of the cycle
  ///
  /// \param orientation  A specific point in the time cycle to query
  /// \param tier         Level at which to retrieve the data (if STORMM is compiled to run on a
  ///                     GPU)
  /// \{
  std::vector<double> getBoxSpaceTransform(CoordinateCycle orientation,
                                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  std::vector<double> getBoxSpaceTransform(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}
  
  /// \brief Get the transformation matrix to take coordinates from fractional space back into
  ///        real space.  The result should be interpreted as a 3x3 matrix in column-major format.
  ///
  /// \param orientation  A specific point in the time cycle to query
  /// \param tier         Level at which to retrieve the data (if STORMM is compiled to run on a
  ///                     GPU)
  /// \{
  std::vector<double> getInverseTransform(CoordinateCycle orientation,
                                          HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  std::vector<double> getInverseTransform(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Get a pointer to one of the coordinate arrays, a pointer to a POINTER-kind Hybrid
  ///        object.
  ///
  /// Overloaded:
  ///   - Specify a stage in the time cycle at which to query coordinates.  This may be necessary,
  ///     step by step, as the object advances in its time cycle.
  ///   - Accept coordinates from the object's current point in the time cycle.  The pointer will
  ///     remain set to this array even if the PhaseSpace object later advances in its time cycle.
  ///   - Get a const pointer to data within a const object, or a non-const pointer to data in a
  ///     mutable object.
  ///
  /// \param dim          Indicate Cartesian X, Y, or Z values
  /// \param kind         Indicate positions, velocities, or forces
  /// \param orientation  The point in the time cycle (WHITE or BLACK) at which to set the pointer
  /// \{
  const Hybrid<double>* getCoordinateHandle(CartesianDimension dim, TrajectoryKind kind,
                                            CoordinateCycle orientation) const;

  const Hybrid<double>* getCoordinateHandle(CartesianDimension dim,
                                            TrajectoryKind kind = TrajectoryKind::POSITIONS) const;

  Hybrid<double>* getCoordinateHandle(CartesianDimension dim, TrajectoryKind kind,
                                      CoordinateCycle orientation);

  Hybrid<double>* getCoordinateHandle(CartesianDimension dim,
                                      TrajectoryKind kind = TrajectoryKind::POSITIONS);
  /// \}

  /// \brief Get a pointer to the box space transform.  Overloading follows from
  ///        getCoordinateHandle(), above.
  /// \{
  const Hybrid<double>* getBoxTransformHandle(CoordinateCycle orientation) const;
  const Hybrid<double>* getBoxTransformHandle() const;
  Hybrid<double>* getBoxTransformHandle(CoordinateCycle orientation);
  Hybrid<double>* getBoxTransformHandle();
  /// \}

  /// \brief Get a pointer to the inverse transform that takes coordinates back into real space.
  ///        Overloading follows from getCoordinateHandle(), above.
  /// \{
  const Hybrid<double>* getInverseTransformHandle(CoordinateCycle orientation) const;
  const Hybrid<double>* getInverseTransformHandle() const;
  Hybrid<double>* getInverseTransformHandle(CoordinateCycle orientation);
  Hybrid<double>* getInverseTransformHandle();
  /// \}

  /// \brief Get a pointer to the  transform that takes coordinates back into real space.
  ///        Overloading follows from getCoordinateHandle(), above.
  /// \{
  const Hybrid<double>* getBoxDimensionsHandle(CoordinateCycle orientation) const;
  const Hybrid<double>* getBoxDimensionsHandle() const;
  Hybrid<double>* getBoxDimensionsHandle(CoordinateCycle orientation);
  Hybrid<double>* getBoxDimensionsHandle();
  /// \}
  
  /// \brief Get a pointer to the ARRAY-kind Hybrid object that holds the actual data.  Needed by
  ///        CoordinateFrame objects which want to be pointers into a PhaseSpace object.
  ///
  /// Overloaded:
  ///   - Get a const pointer to a const PhaseSpace object's data storage
  ///   - Get a non-const pointer to a mutable PhaeSpace object's data storage
  /// \{
  const Hybrid<double>* getStorageHandle() const;
  Hybrid<double>* getStorageHandle();
  /// \}

  /// \brief Initialize the forces (set them to zero)
  ///
  /// Overloaded:
  ///   - Update forces for an arbitrary point in the time cycle
  ///   - Update forces for the object's current position in the cycle
  ///
  /// \param orentation  A selected point in the time cycle (WHITE or BLACK)
  /// \{
  void initializeForces(CoordinateCycle orientation);
  void initializeForces();
  /// \}
  
  /// \brief Update the cycle position.
  ///
  /// Overloaded:
  ///   - Advance the cycle position based on its current setting (no input argument):
  ///     present >> future >> past >> present >> ...
  ///   - Set the cycle position to an arbitrary point
  ///
  /// \param  time_point  The point in the time cycle that shall become the PhaseSpace object's
  ///                     "present" coordinates.
  /// \{
  void updateCyclePosition();
  void updateCyclePosition(CoordinateCycle time_point);
  /// \}
  
  /// \brief Put the phase space data into a trajectory or checkpoint file.
  ///
  /// \param file_name     Name of the file to write
  /// \param current_time  Time progress of the simulation
  /// \param traj_kind     The type of trajectory to print (coordinates, velocities, or forces)
  /// \param output_kind   The format of the file to write (checkpoint files print position and
  ///                      velocity data by obligation, but trajectory files can contain either of
  ///                      these as well as forces)
  /// \param expectation   The condition in which the output file is expected to be found
  void exportToFile(const std::string &file_name, double current_time = 0.0,
                    TrajectoryKind traj_kind = TrajectoryKind::POSITIONS,
                    CoordinateFileKind output_kind = CoordinateFileKind::AMBER_INPCRD,
                    PrintSituation expectation = PrintSituation::UNKNOWN) const;

  /// \brief Get a pointer to the object itself (useful when the object has been passed as a const
  ///        reference and a pointer is needed).
  const PhaseSpace* getSelfPointer() const;

  /// \brief Get the abstract for this object, containing C-style pointers for the most rapid
  ///        access to any of its member variables.
  ///
  /// Overloaded:
  ///   - Get a read-only abstract from a const PhaseSpace object
  ///   - Get a writeable abstract from a mutable PhaseSpace object
  ///   - Get either object oriented with the present, alternate positional arrays set as holding
  ///     the current coordinates (three such abstracts can be rotated over successive cycles of
  ///     dynamics to let the coordinates evolve, protected against race conditions, without
  ///     swapping the actual locations in memory)
  ///
  /// \param tier         Specify pointers on the host or device
  /// \param orientation  Arbitrarily selected point on the time cycle to have the reader or writer
  ///                     take as the current coordinates
  /// \{
  const PhaseSpaceReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  const PhaseSpaceReader data(CoordinateCycle orientation,
                              HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  PhaseSpaceWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  PhaseSpaceWriter data(CoordinateCycle orientation,
                        HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}
  
#ifdef STORMM_USE_HPC
  /// \brief Get an abstract for the object's host data such that all pointers are guaranteed to
  ///        be valid on the GPU.
  ///
  /// Overloaded:
  ///   - Get a read-only abstract for a const PhaseSpace object
  ///   - Get a writeable abstract for a non-const PhaseSpace object
  ///   - Get either object oriented with the present, alternate coordinate arrays
  ///
  /// \param orientation  Arbitrarily selected point on the time cycle to have the reader or writer
  ///                     take as the current coordinates
  /// \{
  const PhaseSpaceReader deviceViewToHostData(CoordinateCycle orientation) const;
  const PhaseSpaceReader deviceViewToHostData() const;
  PhaseSpaceWriter deviceViewToHostData(CoordinateCycle orientation);
  PhaseSpaceWriter deviceViewToHostData();
  /// \}

  /// \brief Upload all information
  void upload();

  /// \brief Upload positional information for one stage of the time cycle.
  ///
  /// Overloaded:
  ///   - Upload particle positions for a specific point in the time cycle
  ///   - Upload particle positions for the object's current position in the cycle
  ///
  /// \param orentation  A selected point in the time cycle (WHITE or BLACK)
  /// \{
  void uploadPositions(CoordinateCycle orientation);
  void uploadPositions();
  /// \}

  /// \brief Upload current transformation matrices.
  void uploadTransformations();

  /// \brief Upload velocity information for one stage of the time cycle.
  ///
  /// Overloaded:
  ///   - Upload velocities for a specific point in the time cycle
  ///   - Upload velocities for the object's current position in the cycle
  ///
  /// \param orentation  A selected point in the time cycle (WHITE or BLACK)
  /// \{
  void uploadVelocities(CoordinateCycle orientation);
  void uploadVelocities();
  /// \}

  /// \brief Upload force information for one stage of the time cycle.
  ///
  /// Overloaded:
  ///   - Upload forces for a specific point in the time cycle
  ///   - Upload forces for the object's current position in the cycle
  ///
  /// \param orentation  A selected point in the time cycle (WHITE or BLACK)
  /// \{
  void uploadForces(CoordinateCycle orientation);
  void uploadForces();
  /// \}

  /// \brief Download all information
  void download();

  /// \brief Download positional information for one stage of the time cycle.
  ///
  /// Overloaded:
  ///   - Download particle positions for a specific point in the time cycle
  ///   - Download particle positions for the object's current position in the cycle
  ///
  /// \param orentation  A selected point in the time cycle (WHITE or BLACK)
  /// \{
  void downloadPositions(CoordinateCycle orientation);
  void downloadPositions();
  /// \}

  /// \brief Download current transformation matrices.
  void downloadTransformations();

  /// \brief Download velocity information for one stage of the time cycle.
  ///
  /// Overloaded:
  ///   - Download velocities for a specific point in the time cycle
  ///   - Download velocities for the object's current position in the cycle
  ///
  /// \param orentation  A selected point in the time cycle (WHITE or BLACK)
  /// \{
  void downloadVelocities(CoordinateCycle orientation);
  void downloadVelocities();
  /// \}

  /// \brief Download force information for one stage of the time cycle.
  ///
  /// Overloaded:
  ///   - Download forces for a specific point in the time cycle
  ///   - Download forces for the object's current position in the cycle
  ///
  /// \param orentation  A selected point in the time cycle (WHITE or BLACK)
  /// \{
  void downloadForces(CoordinateCycle orientation);
  void downloadForces();
  /// \}
#endif
  
private:
  HybridFormat format;             ///< The layout of memory in the object
  std::string file_name;           ///< Name of the file from which these coordinates (and
                                   ///<   perhaps velocities) derived.  Empty string indicates
                                   ///<   no file.
  int atom_count;                  ///< The number of atoms in the system
  UnitCellType unit_cell;          ///< The type of unit cell
  CoordinateCycle cycle_position;  ///< Indicates the place in the past >> present >> future
                                   ///<   cycle where the object is currently storing its
                                   ///<   relevant coordinates.  After a dynamics step, present
                                   ///<   becomes past, future becomes present, and what were
                                   ///<   the arrays holding past coordinates stand ready to
                                   ///<   accept the future configuration as it is assembled.
                                   ///<   With the work unit system, constraint updates might
                                   ///<   not occur atomically.  Three arrays, not just two as
                                   ///<   are used in other codes, are needed to protect
                                   ///<   against race conditions.
  Hybrid<double> x_coordinates;    ///< Cartesian X coordinates of all particles
  Hybrid<double> y_coordinates;    ///< Cartesian Y coordinates of all particles
  Hybrid<double> z_coordinates;    ///< Cartesian Z coordinates of all particles

  // The actual nature of the following arrays, as well as the [x,y,z]_coordinates themselves, can
  // change based on the CoordinateCycle position.
  Hybrid<double> x_alt_coordinates;  ///< Previous step Cartesian X coordinates of all particles
  Hybrid<double> y_alt_coordinates;  ///< Previous step Cartesian Y coordinates of all particles
  Hybrid<double> z_alt_coordinates;  ///< Previous step Cartesian Z coordinates of all particles

  // Transformation matrices are given directly after the coordinates in the overall order of data
  Hybrid<double> box_space_transform;  ///< Matrix to transform coordinates into box space (3 x 3)
  Hybrid<double> inverse_transform;    ///< Matrix to transform coordinates into real space (3 x 3)
  Hybrid<double> box_dimensions;       ///< Three lengths and three angles defining the box
                                       ///<   (lengths are given in Angstroms, angles in radians)

  // Alternate transformation matrices are stored next, for consistency.  This provides a complete
  // record of the previous unit cell dimensions in the event that the trajectory needs to
  // backtrack against a rejected Monte-Carlo move, or a complete space to build the next state of
  // the system.
  Hybrid<double> alt_box_space_transform;  ///< Matrix to transform coordinates into box space
  Hybrid<double> alt_inverse_transform;    ///< Matrix to transform coordinates into real space
  Hybrid<double> alt_box_dimensions;       ///< Three lengths and three angles defining the box
                                           ///<   (lengths are given in Angstroms, angles in
                                           ///<   radians)

  // Like coordinates, velocities and forces appear in WHITE -> BLACK blocks
  Hybrid<double> x_velocities;      ///< Cartesian X velocities of all particles
  Hybrid<double> y_velocities;      ///< Cartesian Y velocities of all particles
  Hybrid<double> z_velocities;      ///< Cartesian Z velocities of all particles
  Hybrid<double> x_alt_velocities;  ///< Alternate Cartesian X velocities of all particles
  Hybrid<double> y_alt_velocities;  ///< Alternate Cartesian Y velocities of all particles
  Hybrid<double> z_alt_velocities;  ///< Alternate Cartesian Z velocities of all particles
  Hybrid<double> x_forces;          ///< Cartesian X forces acting on all particles
  Hybrid<double> y_forces;          ///< Cartesian Y forces acting on all particles
  Hybrid<double> z_forces;          ///< Cartesian Z forces acting on all particles 
  Hybrid<double> x_alt_forces;      ///< Alternate Cartesian X forces acting on all particles
  Hybrid<double> y_alt_forces;      ///< Alternate Cartesian Y forces acting on all particles
  Hybrid<double> z_alt_forces;      ///< Alternate Cartesian Z forces acting on all particles

  /// All of the above Hybrid objects are pointers into this single large array, segmented to hold
  /// each type of information with zero-padding to accommodate the HPC warp size.
  Hybrid<double> storage;

  /// \brief Allocate space for the object, based on a known number of atoms
  void allocate();
};

/// \brief Define the type index for the PhaseSpace object.
static const size_t phasespace_type_index = std::type_index(typeid(PhaseSpace)).hash_code();

/// \brief Interpret the inverse transformation matrix to determine whether a unit cell exists, and
///        if it is rectilinear or not.
///
/// Overloaded:
///   - Accept a C-style pointer to the array of matrix elements
///   - Accept a Standard Template Library Vector or a Hybrid object
///
/// \param inverse_transform  The inverse transformation matrix
/// \{
UnitCellType determineUnitCellTypeByShape(const double* inverse_transform);
UnitCellType determineUnitCellTypeByShape(const std::vector<double> &inverse_transform);
UnitCellType determineUnitCellTypeByShape(const Hybrid<double> &inverse_transform);
/// \}

/// \brief Interlace three arrays of X, Y, and Z coordinates (i.e. positions, velocities, or
///        forces) into a result ordered { X(0), Y(0), Z(0), X(1), Y(1), ..., Y(N), Z(N) }.
///
/// \param xptr        Pointer to Cartesian X data
/// \param yptr        Pointer to Cartesian Y data
/// \param zptr        Pointer to Cartesian Z data
/// \param low_index   Starting index in each array
/// \param high_index  Upper limit of data to take from each array
std::vector<double> interlaceXYZ(const double* xptr, const double* yptr, const double* zptr,
                                 int low_index, int high_index);
  
} // namespace trajectory
} // namespace stormm

#include "phasespace.tpp"

// As with common types and STORMM vector types, define the type indices for general use in the
// STORMM namespace.
namespace stormm {
using trajectory::phasespace_type_index;
} // namespace stormm

#endif
