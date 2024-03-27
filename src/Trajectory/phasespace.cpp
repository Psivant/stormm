#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "FileManagement/file_listing.h"
#include "FileManagement/file_util.h"
#include "Math/rounding.h"
#include "Parsing/textfile.h"
#include "amber_ascii.h"
#include "phasespace.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

using card::HybridKind;
using constants::CartesianDimension;
using diskutil::detectCoordinateFileKind;
using diskutil::DataFormat;
using diskutil::DrivePathType;
using diskutil::getDrivePathType;
using diskutil::getTrajectoryFormat;
using parse::TextFile;
using stmath::roundUp;

//-------------------------------------------------------------------------------------------------
PhaseSpaceWriter::PhaseSpaceWriter(const int natom_in, const UnitCellType unit_cell_in,
                                   double* xcrd_in, double* ycrd_in, double* zcrd_in,
                                   double* umat_in, double* invu_in, double* boxdim_in,
                                   double* umat_alt_in, double* invu_alt_in, double* boxdim_alt_in,
                                   double* xvel_in, double* yvel_in, double* zvel_in,
                                   double* xfrc_in, double* yfrc_in, double* zfrc_in,
                                   double* xalt_in, double* yalt_in, double* zalt_in,
                                   double* vxalt_in, double* vyalt_in, double* vzalt_in,
                                   double* fxalt_in, double* fyalt_in, double* fzalt_in) :
    natom{natom_in}, unit_cell{unit_cell_in}, xcrd{xcrd_in}, ycrd{ycrd_in}, zcrd{zcrd_in},
    umat{umat_in}, invu{invu_in}, boxdim{boxdim_in}, umat_alt{umat_alt_in}, invu_alt{invu_alt_in},
    boxdim_alt{boxdim_alt_in}, xvel{xvel_in}, yvel{yvel_in}, zvel{zvel_in}, xfrc{xfrc_in},
    yfrc{yfrc_in}, zfrc{zfrc_in}, xalt{xalt_in}, yalt{yalt_in}, zalt{zalt_in}, vxalt{vxalt_in},
    vyalt{vyalt_in}, vzalt{vzalt_in}, fxalt{fxalt_in}, fyalt{fyalt_in}, fzalt{fzalt_in}
{}

//-------------------------------------------------------------------------------------------------
PhaseSpaceReader::PhaseSpaceReader(const int natom_in, const UnitCellType unit_cell_in,
                                   const double* xcrd_in, const double* ycrd_in,
                                   const double* zcrd_in, const double* umat_in,
                                   const double* invu_in, const double* boxdim_in,
                                   const double* umat_alt_in, const double* invu_alt_in,
                                   const double* boxdim_alt_in, const double* xvel_in,
                                   const double* yvel_in, const double* zvel_in,
                                   const double* xfrc_in, const double* yfrc_in,
                                   const double* zfrc_in, const double* xalt_in,
                                   const double* yalt_in, const double* zalt_in,
                                   const double* vxalt_in, const double* vyalt_in,
                                   const double* vzalt_in, const double* fxalt_in,
                                   const double* fyalt_in, const double* fzalt_in) :
    natom{natom_in}, unit_cell{unit_cell_in}, xcrd{xcrd_in}, ycrd{ycrd_in}, zcrd{zcrd_in},
    umat{umat_in}, invu{invu_in}, boxdim{boxdim_in}, umat_alt{umat_alt_in}, invu_alt{invu_alt_in},
    boxdim_alt{boxdim_alt_in}, xvel{xvel_in}, yvel{yvel_in}, zvel{zvel_in}, xfrc{xfrc_in},
    yfrc{yfrc_in}, zfrc{zfrc_in}, xalt{xalt_in}, yalt{yalt_in}, zalt{zalt_in}, vxalt{vxalt_in},
    vyalt{vyalt_in}, vzalt{vzalt_in}, fxalt{fxalt_in}, fyalt{fyalt_in}, fzalt{fzalt_in}
{}

//-------------------------------------------------------------------------------------------------
PhaseSpaceReader::PhaseSpaceReader(const PhaseSpaceWriter &psw) :
    natom{psw.natom},
    unit_cell{psw.unit_cell},
    xcrd{psw.xcrd}, ycrd{psw.ycrd}, zcrd{psw.zcrd},
    umat{psw.umat}, invu{psw.invu}, boxdim{psw.boxdim},
    umat_alt{psw.umat_alt}, invu_alt{psw.invu_alt}, boxdim_alt{psw.boxdim_alt},
    xvel{psw.xvel}, yvel{psw.yvel}, zvel{psw.zvel},
    xfrc{psw.xfrc}, yfrc{psw.yfrc}, zfrc{psw.zfrc},
    xalt{psw.xalt}, yalt{psw.yalt}, zalt{psw.zalt},
    vxalt{psw.vxalt}, vyalt{psw.vyalt}, vzalt{psw.vzalt},
    fxalt{psw.fxalt}, fyalt{psw.fyalt}, fzalt{psw.fzalt}
{}

//-------------------------------------------------------------------------------------------------
PhaseSpace::PhaseSpace(const int atom_count_in, const UnitCellType unit_cell_in) :
    file_name{std::string("")},
    atom_count{atom_count_in},
    unit_cell{unit_cell_in},
    cycle_position{CoordinateCycle::WHITE},
    x_coordinates{HybridKind::POINTER, "x_coordinates"},
    y_coordinates{HybridKind::POINTER, "y_coordinates"},
    z_coordinates{HybridKind::POINTER, "z_coordinates"},
    x_alt_coordinates{HybridKind::POINTER, "x_alt_coords"},
    y_alt_coordinates{HybridKind::POINTER, "y_alt_coords"},
    z_alt_coordinates{HybridKind::POINTER, "z_alt_coords"},
    box_space_transform{HybridKind::POINTER, "box_transform"},
    inverse_transform{HybridKind::POINTER, "inv_transform"},
    box_dimensions{HybridKind::POINTER, "box_dimensions"},
    alt_box_space_transform{HybridKind::POINTER, "alt_box_xform"},
    alt_inverse_transform{HybridKind::POINTER, "alt_inv_xform"},
    alt_box_dimensions{HybridKind::POINTER, "alt_box_dims"},
    x_velocities{HybridKind::POINTER, "x_velocities"},
    y_velocities{HybridKind::POINTER, "y_velocities"},
    z_velocities{HybridKind::POINTER, "z_velocities"},
    x_alt_velocities{HybridKind::POINTER, "x_alt_vels"},
    y_alt_velocities{HybridKind::POINTER, "y_alt_vels"},
    z_alt_velocities{HybridKind::POINTER, "z_alt_vels"},
    x_forces{HybridKind::POINTER, "x_forces"},
    y_forces{HybridKind::POINTER, "y_forces"},
    z_forces{HybridKind::POINTER, "z_forces"},
    x_alt_forces{HybridKind::POINTER, "x_alt_forces"},
    y_alt_forces{HybridKind::POINTER, "y_alt_forces"},
    z_alt_forces{HybridKind::POINTER, "z_alt_forces"},
    storage{HybridKind::ARRAY, "phase_space_data"}
{
  allocate();
}

//-------------------------------------------------------------------------------------------------
PhaseSpace::PhaseSpace(const std::string &file_name_in, const CoordinateFileKind file_kind,
                       const int frame_number) :
    PhaseSpace()
{
  buildFromFile(file_name_in, file_kind, frame_number);
}

//-------------------------------------------------------------------------------------------------
PhaseSpace::PhaseSpace(const std::string &file_name_in, const AtomGraph &ag,
                       const CoordinateFileKind file_kind, const int frame_number) :
    PhaseSpace(file_name_in, file_kind, frame_number)
{
  // Check that the coordinates agree with the topology
  if (ag.getAtomCount() != atom_count) {
    rtErr("Coordinates in file " + file_name + " describe " + std::to_string(atom_count) +
          "atoms, whereas topology " + ag.getFileName() + " describes a system of " +
          std::to_string(atom_count) + " atoms.", "PhaseSpace");
  }
}

//-------------------------------------------------------------------------------------------------
PhaseSpace::PhaseSpace(const PhaseSpace &original) :
    file_name{original.file_name},
    atom_count{original.atom_count},
    unit_cell{original.unit_cell},
    cycle_position{original.cycle_position},
    x_coordinates{original.x_coordinates},
    y_coordinates{original.y_coordinates},
    z_coordinates{original.z_coordinates},
    x_alt_coordinates{original.x_alt_coordinates},
    y_alt_coordinates{original.y_alt_coordinates},
    z_alt_coordinates{original.z_alt_coordinates},
    box_space_transform{original.box_space_transform},
    inverse_transform{original.inverse_transform},
    box_dimensions{original.box_dimensions},
    alt_box_space_transform{original.alt_box_space_transform},
    alt_inverse_transform{original.alt_inverse_transform},
    alt_box_dimensions{original.alt_box_dimensions},
    x_velocities{original.x_velocities},
    y_velocities{original.y_velocities},
    z_velocities{original.z_velocities},
    x_alt_velocities{original.x_alt_velocities},
    y_alt_velocities{original.y_alt_velocities},
    z_alt_velocities{original.z_alt_velocities},
    x_forces{original.x_forces},
    y_forces{original.y_forces},
    z_forces{original.z_forces},
    x_alt_forces{original.x_alt_forces},
    y_alt_forces{original.y_alt_forces},
    z_alt_forces{original.z_alt_forces},
    storage{original.storage}
{
  // Set the POINTER-kind Hybrids in the new object appropriately.  Even the resize() operation
  // inherent to "allocate" will pass by with little more than a check that the length of the data
  // storage array is already what it should be.
  allocate();
}

//-------------------------------------------------------------------------------------------------
PhaseSpace& PhaseSpace::operator=(const PhaseSpace &other) {
  
  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // Copy the file name (if applicable) and atom count
  file_name = other.file_name;
  atom_count = other.atom_count;
  unit_cell = other.unit_cell;
  cycle_position = other.cycle_position;

  // Copy the Hybrid objects to preserve tags and the proper kinds.  Copying the data storage
  // object carries over all of the other object's contents.
  x_coordinates = other.x_coordinates;
  y_coordinates = other.y_coordinates;
  z_coordinates = other.z_coordinates;
  x_alt_coordinates = other.x_alt_coordinates;
  y_alt_coordinates = other.y_alt_coordinates;
  z_alt_coordinates = other.z_alt_coordinates;
  box_space_transform = other.box_space_transform;
  inverse_transform = other.inverse_transform;
  box_dimensions = other.box_dimensions;
  alt_box_space_transform = other.alt_box_space_transform;
  alt_inverse_transform = other.alt_inverse_transform;
  alt_box_dimensions = other.alt_box_dimensions;
  x_velocities = other.x_velocities;
  y_velocities = other.y_velocities;
  z_velocities = other.z_velocities;
  x_alt_velocities = other.x_alt_velocities;
  y_alt_velocities = other.y_alt_velocities;
  z_alt_velocities = other.z_alt_velocities;
  x_forces = other.x_forces;
  y_forces = other.y_forces;
  z_forces = other.z_forces;
  x_alt_forces = other.x_alt_forces;
  y_alt_forces = other.y_alt_forces;
  z_alt_forces = other.z_alt_forces;
  storage = other.storage;

  // Allocating again sets the internal POINTER-kind Hybrid objects
  allocate();
  return *this;
}

//-------------------------------------------------------------------------------------------------
PhaseSpace::PhaseSpace(PhaseSpace &&original) :
    file_name{std::move(original.file_name)},
    atom_count{original.atom_count},
    unit_cell{original.unit_cell},
    cycle_position{original.cycle_position},
    x_coordinates{std::move(original.x_coordinates)},
    y_coordinates{std::move(original.y_coordinates)},
    z_coordinates{std::move(original.z_coordinates)},
    x_alt_coordinates{std::move(original.x_alt_coordinates)},
    y_alt_coordinates{std::move(original.y_alt_coordinates)},
    z_alt_coordinates{std::move(original.z_alt_coordinates)},
    box_space_transform{std::move(original.box_space_transform)},
    inverse_transform{std::move(original.inverse_transform)},
    box_dimensions{std::move(original.box_dimensions)},
    alt_box_space_transform{std::move(original.alt_box_space_transform)},
    alt_inverse_transform{std::move(original.alt_inverse_transform)},
    alt_box_dimensions{std::move(original.alt_box_dimensions)},
    x_velocities{std::move(original.x_velocities)},
    y_velocities{std::move(original.y_velocities)},
    z_velocities{std::move(original.z_velocities)},
    x_alt_velocities{std::move(original.x_alt_velocities)},
    y_alt_velocities{std::move(original.y_alt_velocities)},
    z_alt_velocities{std::move(original.z_alt_velocities)},
    x_forces{std::move(original.x_forces)},
    y_forces{std::move(original.y_forces)},
    z_forces{std::move(original.z_forces)},
    x_alt_forces{std::move(original.x_alt_forces)},
    y_alt_forces{std::move(original.y_alt_forces)},
    z_alt_forces{std::move(original.z_alt_forces)},
    storage{std::move(original.storage)}
{
  // No repair of the pointers is necessary, as they still point to valid data
}

//-------------------------------------------------------------------------------------------------
PhaseSpace& PhaseSpace::operator=(PhaseSpace &&other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  file_name = std::move(other.file_name);
  atom_count = other.atom_count;
  unit_cell = other.unit_cell;
  cycle_position = other.cycle_position;
  x_coordinates = std::move(other.x_coordinates);
  y_coordinates = std::move(other.y_coordinates);
  z_coordinates = std::move(other.z_coordinates);
  x_alt_coordinates = std::move(other.x_alt_coordinates);
  y_alt_coordinates = std::move(other.y_alt_coordinates);
  z_alt_coordinates = std::move(other.z_alt_coordinates);
  box_space_transform = std::move(other.box_space_transform);
  inverse_transform = std::move(other.inverse_transform);
  box_dimensions = std::move(other.box_dimensions);
  alt_box_space_transform = std::move(other.alt_box_space_transform);
  alt_inverse_transform = std::move(other.alt_inverse_transform);
  alt_box_dimensions = std::move(other.alt_box_dimensions);
  x_velocities = std::move(other.x_velocities);
  y_velocities = std::move(other.y_velocities);
  z_velocities = std::move(other.z_velocities);
  x_alt_velocities = std::move(other.x_alt_velocities);
  y_alt_velocities = std::move(other.y_alt_velocities);
  z_alt_velocities = std::move(other.z_alt_velocities);
  x_forces = std::move(other.x_forces);
  y_forces = std::move(other.y_forces);
  z_forces = std::move(other.z_forces);
  x_alt_forces = std::move(other.x_alt_forces);
  y_alt_forces = std::move(other.y_alt_forces);
  z_alt_forces = std::move(other.z_alt_forces);
  storage = std::move(other.storage);
  return *this;
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::buildFromFile(const std::string &file_name_in, const CoordinateFileKind file_kind,
                               const int frame_number) {
  file_name = file_name_in;
  
  // Try to detect the file format if it is not already specified.  If it remains UNKNOWN, that
  // will ultimately lead to an error.
  CoordinateFileKind actual_kind = file_kind;
  if (file_kind == CoordinateFileKind::UNKNOWN) {
    actual_kind = detectCoordinateFileKind(file_name);
  }
  switch (actual_kind) {
  case CoordinateFileKind::AMBER_CRD:
    {
      // The number of atoms must be known a-priori in order to read from a .crd trajectory file.
      if (atom_count == 0) {
        rtErr("A number of atoms matching the trajectory must be known prior to reading a .crd "
              "file.", "PhaseSpace", "buildFromFile");
      }
      allocate();
      TextFile tf(file_name);
      readAmberCrdFormat(tf, &x_coordinates, &y_coordinates, &z_coordinates, unit_cell,
                         &box_space_transform, &inverse_transform, &box_dimensions, frame_number);
    }
    break;
  case CoordinateFileKind::AMBER_INPCRD:
    {
      TextFile tf(file_name);
      atom_count = getAmberRestartAtomCount(tf);
      allocate();
      getAmberInputCoordinates(tf, &x_coordinates, &y_coordinates, &z_coordinates,
                               &box_space_transform, &inverse_transform, &box_dimensions);

      // Interpret the box transformation, updating the unit cell type based on the file
      unit_cell = determineUnitCellTypeByShape(inverse_transform.data());
    }
    break;
  case CoordinateFileKind::AMBER_ASCII_RST:
    {
      TextFile tf(file_name);
      atom_count = getAmberRestartAtomCount(tf);
      allocate();
      getAmberInputCoordinates(tf, &x_coordinates, &y_coordinates, &z_coordinates,
                               &box_space_transform, &inverse_transform, &box_dimensions);
      const int data_line_count = (atom_count + 1) / 2;

      // The file may have been marked as an Amber ASCII restart because of the presence of a time
      // stamp on the second line when no velocities are actually present.  Check that here.
      if (tf.getLineCount() >= 2 + (2 * data_line_count)) {
        getAmberRestartVelocities(tf, &x_velocities, &y_velocities, &z_velocities);
      }

      // Interpret the box transformation, updating the unit cell type based on the file
      unit_cell = determineUnitCellTypeByShape(inverse_transform.data());
    }
    break;
  case CoordinateFileKind::SDF:
    rtErr("An SD file (BIOVIA .sdf) cannot be used to create a PhaseSpace object directly.  "
          "Create an MdlMol object and export a PhaseSpace object from it instead.",
          "PhaseSpace", "buildFromFile");
    break;
  case CoordinateFileKind::AMBER_NETCDF:
    break;
  case CoordinateFileKind::AMBER_NETCDF_RST:
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("The coordinate file type of " + file_name + " could not be understood.", "PhaseSpace",
          "buildFromFile");
  }
}

//-------------------------------------------------------------------------------------------------
std::string PhaseSpace::getFileName() const {
  return file_name;
}

//-------------------------------------------------------------------------------------------------
int PhaseSpace::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
UnitCellType PhaseSpace::getUnitCellType() const {
  return unit_cell;
}

//-------------------------------------------------------------------------------------------------
CoordinateCycle PhaseSpace::getCyclePosition() const {
  return cycle_position;
}

//-------------------------------------------------------------------------------------------------
const double* PhaseSpace::getCoordinatePointer(const CartesianDimension dim,
                                               const TrajectoryKind kind,
                                               const HybridTargetLevel tier) const {
  switch (dim) {
  case CartesianDimension::X:
    switch (kind) {
    case TrajectoryKind::POSITIONS:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return x_alt_coordinates.data(tier);
      case CoordinateCycle::WHITE:
        return x_coordinates.data(tier);
      }
      break;
    case TrajectoryKind::VELOCITIES:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return x_alt_velocities.data(tier);
      case CoordinateCycle::WHITE:
        return x_velocities.data(tier);
      }
      break;
    case TrajectoryKind::FORCES:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return x_alt_forces.data(tier);
      case CoordinateCycle::WHITE:
        return x_forces.data(tier);
      }
      break;
    }
    break;
  case CartesianDimension::Y:
    switch (kind) {
    case TrajectoryKind::POSITIONS:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return y_alt_coordinates.data(tier);
      case CoordinateCycle::WHITE:
        return y_coordinates.data(tier);
      }
      break;
    case TrajectoryKind::VELOCITIES:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return y_alt_velocities.data(tier);
      case CoordinateCycle::WHITE:
        return y_velocities.data(tier);
      }
      break;
    case TrajectoryKind::FORCES:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return y_alt_forces.data(tier);
      case CoordinateCycle::WHITE:
        return y_forces.data(tier);
      }
      break;
    }
    break;
  case CartesianDimension::Z:
    switch (kind) {
    case TrajectoryKind::POSITIONS:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return z_alt_coordinates.data(tier);
      case CoordinateCycle::WHITE:
        return z_coordinates.data(tier);
      }
      break;
    case TrajectoryKind::VELOCITIES:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return z_alt_velocities.data(tier);
      case CoordinateCycle::WHITE:
        return z_velocities.data(tier);
      }
      break;
    case TrajectoryKind::FORCES:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return z_alt_forces.data(tier);
      case CoordinateCycle::WHITE:
        return z_forces.data(tier);
      }
      break;
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double* PhaseSpace::getCoordinatePointer(const CartesianDimension dim, const TrajectoryKind kind,
                                         const HybridTargetLevel tier) {
  switch (dim) {
  case CartesianDimension::X:
    switch (kind) {
    case TrajectoryKind::POSITIONS:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return x_alt_coordinates.data(tier);
      case CoordinateCycle::WHITE:
        return x_coordinates.data(tier);
      }
      break;
    case TrajectoryKind::VELOCITIES:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return x_alt_velocities.data(tier);
      case CoordinateCycle::WHITE:
        return x_velocities.data(tier);
      }
      break;
    case TrajectoryKind::FORCES:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return x_alt_forces.data(tier);
      case CoordinateCycle::WHITE:
        return x_forces.data(tier);
      }
      break;
    }
    break;
  case CartesianDimension::Y:
    switch (kind) {
    case TrajectoryKind::POSITIONS:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return y_alt_coordinates.data(tier);
      case CoordinateCycle::WHITE:
        return y_coordinates.data(tier);
      }
      break;
    case TrajectoryKind::VELOCITIES:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return y_alt_velocities.data(tier);
      case CoordinateCycle::WHITE:
        return y_velocities.data(tier);
      }
      break;
    case TrajectoryKind::FORCES:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return y_alt_forces.data(tier);
      case CoordinateCycle::WHITE:
        return y_forces.data(tier);
      }
      break;
    }
    break;
  case CartesianDimension::Z:
    switch (kind) {
    case TrajectoryKind::POSITIONS:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return z_alt_coordinates.data(tier);
      case CoordinateCycle::WHITE:
        return z_coordinates.data(tier);
      }
      break;
    case TrajectoryKind::VELOCITIES:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return z_alt_velocities.data(tier);
      case CoordinateCycle::WHITE:
        return z_velocities.data(tier);
      }
      break;
    case TrajectoryKind::FORCES:
      switch (cycle_position) {
      case CoordinateCycle::BLACK:
        return z_alt_forces.data(tier);
      case CoordinateCycle::WHITE:
        return z_forces.data(tier);
      }
      break;
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> PhaseSpace::getInterlacedCoordinates(const TrajectoryKind kind,
                                                         const HybridTargetLevel tier) const {
  return getInterlacedCoordinates(0, atom_count, cycle_position, kind, tier);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> PhaseSpace::getInterlacedCoordinates(const CoordinateCycle orientation,
                                                         const TrajectoryKind kind,
                                                         const HybridTargetLevel tier) const {
  return getInterlacedCoordinates(0, atom_count, orientation, kind, tier);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> PhaseSpace::getInterlacedCoordinates(const int low_index, const int high_index,
                                                         const TrajectoryKind kind,
                                                         const HybridTargetLevel tier) const {
  return getInterlacedCoordinates(low_index, high_index, cycle_position, kind, tier);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> PhaseSpace::getInterlacedCoordinates(const int low_index, const int high_index,
                                                         const CoordinateCycle orientation,
                                                         const TrajectoryKind kind,
                                                         const HybridTargetLevel tier) const {

  // Range check as this will use pointers
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("Invalid atom range " + std::to_string(low_index) + " to " + std::to_string(high_index) +
          " in an object with " + std::to_string(atom_count) + " atoms.", "PhaseSpace",
          "getInterlacedCoordinates");
  }
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (kind) {
    case TrajectoryKind::POSITIONS:
      switch (orientation) {
      case CoordinateCycle::BLACK:
        return interlaceXYZ(x_alt_coordinates.data(), y_alt_coordinates.data(),
                            z_alt_coordinates.data(), low_index, high_index);
      case CoordinateCycle::WHITE:
        return interlaceXYZ(x_coordinates.data(), y_coordinates.data(), z_coordinates.data(),
                            low_index, high_index);
      }
      break;
    case TrajectoryKind::VELOCITIES:
      switch (orientation) {
      case CoordinateCycle::BLACK:
        return interlaceXYZ(x_alt_velocities.data(), y_alt_velocities.data(),
                            z_alt_velocities.data(), low_index, high_index);
      case CoordinateCycle::WHITE:
        return interlaceXYZ(x_velocities.data(), y_velocities.data(), z_velocities.data(),
                            low_index, high_index);
      }
      break;
    case TrajectoryKind::FORCES:
      switch (orientation) {
      case CoordinateCycle::BLACK:
        return interlaceXYZ(x_alt_forces.data(), y_alt_forces.data(), z_alt_forces.data(),
                            low_index, high_index);
      case CoordinateCycle::WHITE:
        return interlaceXYZ(x_forces.data(), y_forces.data(), z_forces.data(), low_index,
                            high_index);
      }
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      std::vector<double> xval, yval, zval;
      switch (kind) {
      case TrajectoryKind::POSITIONS:
        switch (orientation) {
        case CoordinateCycle::BLACK:
          xval = x_alt_coordinates.readDevice(low_index, high_index);
          yval = y_alt_coordinates.readDevice(low_index, high_index);
          zval = z_alt_coordinates.readDevice(low_index, high_index);
          break;
        case CoordinateCycle::WHITE:
          xval = x_coordinates.readDevice(low_index, high_index);
          yval = y_coordinates.readDevice(low_index, high_index);
          zval = z_coordinates.readDevice(low_index, high_index);
          break;
        }
        break;
      case TrajectoryKind::VELOCITIES:
        switch (orientation) {
        case CoordinateCycle::BLACK:
          xval = x_alt_velocities.readDevice(low_index, high_index);
          yval = y_alt_velocities.readDevice(low_index, high_index);
          zval = z_alt_velocities.readDevice(low_index, high_index);
          break;
        case CoordinateCycle::WHITE:
          xval = x_velocities.readDevice(low_index, high_index);
          yval = y_velocities.readDevice(low_index, high_index);
          zval = z_velocities.readDevice(low_index, high_index);
          break;
        }
        break;
      case TrajectoryKind::FORCES:
        switch (orientation) {
        case CoordinateCycle::BLACK:
          xval = x_alt_forces.readDevice(low_index, high_index);
          yval = y_alt_forces.readDevice(low_index, high_index);
          zval = z_alt_forces.readDevice(low_index, high_index);
          break;
        case CoordinateCycle::WHITE:
          xval = x_forces.readDevice(low_index, high_index);
          yval = y_forces.readDevice(low_index, high_index);
          zval = z_forces.readDevice(low_index, high_index);
          break;
        }
        break;
      }
      return interlaceXYZ(xval.data(), yval.data(), zval.data(), low_index, high_index);
    }
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const double* PhaseSpace::getBoxSpaceTransformPointer(const CoordinateCycle orientation,
                                                      const HybridTargetLevel tier) const {
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return alt_box_space_transform.data(tier);
  case CoordinateCycle::WHITE:
    return box_space_transform.data(tier);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double* PhaseSpace::getBoxSpaceTransformPointer(const CoordinateCycle orientation,
                                                const HybridTargetLevel tier) {
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return alt_box_space_transform.data(tier);
  case CoordinateCycle::WHITE:
    return box_space_transform.data(tier);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const double* PhaseSpace::getBoxSpaceTransformPointer(const HybridTargetLevel tier) const {
  return getBoxSpaceTransformPointer(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
double* PhaseSpace::getBoxSpaceTransformPointer(const HybridTargetLevel tier) {
  return getBoxSpaceTransformPointer(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
const double* PhaseSpace::getInverseTransformPointer(const CoordinateCycle orientation,
                                                     const HybridTargetLevel tier) const {
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return alt_inverse_transform.data(tier);
  case CoordinateCycle::WHITE:
    return inverse_transform.data(tier);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double* PhaseSpace::getInverseTransformPointer(const CoordinateCycle orientation,
                                               const HybridTargetLevel tier) {
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return alt_inverse_transform.data(tier);
  case CoordinateCycle::WHITE:
    return inverse_transform.data(tier);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const double* PhaseSpace::getInverseTransformPointer(const HybridTargetLevel tier) const {
  return getInverseTransformPointer(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
double* PhaseSpace::getInverseTransformPointer(const HybridTargetLevel tier) {
  return getInverseTransformPointer(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
const double* PhaseSpace::getBoxSizePointer(const CoordinateCycle orientation,
                                            const HybridTargetLevel tier) const {
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return alt_box_dimensions.data(tier);
  case CoordinateCycle::WHITE:
    return box_dimensions.data(tier);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double* PhaseSpace::getBoxSizePointer(const CoordinateCycle orientation,
                                      const HybridTargetLevel tier) {
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return alt_box_dimensions.data(tier);
  case CoordinateCycle::WHITE:
    return box_dimensions.data(tier);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const double* PhaseSpace::getBoxSizePointer(const HybridTargetLevel tier) const {
  return getBoxSizePointer(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
double* PhaseSpace::getBoxSizePointer(const HybridTargetLevel tier) {
  return getBoxSizePointer(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> PhaseSpace::getBoxSpaceTransform(const CoordinateCycle orientation,
                                                     const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (orientation) {
    case CoordinateCycle::BLACK:
      return alt_box_space_transform.readHost();
    case CoordinateCycle::WHITE:
      return box_space_transform.readHost();
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (orientation) {
    case CoordinateCycle::BLACK:
      return alt_box_space_transform.readDevice();
    case CoordinateCycle::WHITE:
      return box_space_transform.readDevice();
    }
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> PhaseSpace::getBoxSpaceTransform(const HybridTargetLevel tier) const {
  return getBoxSpaceTransform(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> PhaseSpace::getInverseTransform(const CoordinateCycle orientation,
                                                    const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (orientation) {
    case CoordinateCycle::BLACK:
      return alt_inverse_transform.readHost();
    case CoordinateCycle::WHITE:
      return inverse_transform.readHost();
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (orientation) {
    case CoordinateCycle::BLACK:
      return alt_inverse_transform.readDevice();
    case CoordinateCycle::WHITE:
      return inverse_transform.readDevice();
    }
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> PhaseSpace::getInverseTransform(const HybridTargetLevel tier) const {
  return getInverseTransform(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
const Hybrid<double>* PhaseSpace::getStoragePointer() const {
  return &storage;
}

//-------------------------------------------------------------------------------------------------
Hybrid<double>* PhaseSpace::getStoragePointer() {
  return &storage;
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::initializeForces(const CoordinateCycle orientation) {
  double* x_ptr;
  double* y_ptr;
  double* z_ptr;
  switch (orientation) {
  case CoordinateCycle::BLACK:
    x_ptr = x_alt_forces.data();
    y_ptr = y_alt_forces.data();
    z_ptr = z_alt_forces.data();
    break;
  case CoordinateCycle::WHITE:
    x_ptr = x_forces.data();
    y_ptr = y_forces.data();
    z_ptr = z_forces.data();
    break;
  }
  for (int i = 0; i < atom_count; i++) {
    x_ptr[i] = 0.0;
    y_ptr[i] = 0.0;
    z_ptr[i] = 0.0;
  }
#ifdef STORMM_USE_HPC
  uploadForces(orientation);
#endif
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::initializeForces() {
  initializeForces(cycle_position);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::updateCyclePosition() {
  switch (cycle_position) {
  case CoordinateCycle::BLACK:
    cycle_position = CoordinateCycle::WHITE;
    break;
  case CoordinateCycle::WHITE:
    cycle_position = CoordinateCycle::BLACK;
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::updateCyclePosition(const CoordinateCycle time_point) {
  cycle_position = time_point;
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::exportToFile(const std::string &file_name, const double current_time,
                              const TrajectoryKind traj_kind, const CoordinateFileKind output_kind,
                              const PrintSituation expectation) const {

  // Filter bad inputs on file appending behavior and an unknown file type
  const PrintSituation aexp = adjustTrajectoryOpeningProtocol(expectation, output_kind,
                                                              "PhaseSpace", "exportToFile");

  // Set the output type (ASCII or binary)
  const DataFormat style = getTrajectoryFormat(output_kind);

  // Begin writing the file
  const bool fi_exists = (getDrivePathType(file_name) == DrivePathType::FILE);
  std::ofstream foutp;
  foutp = openOutputFile(file_name, aexp, "Open an output file for writing PhaseSpace contents",
                         style);
  switch (output_kind) {
  case CoordinateFileKind::AMBER_CRD:
  case CoordinateFileKind::SDF:
  case CoordinateFileKind::AMBER_NETCDF:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    if (fi_exists == false ||
        aexp == PrintSituation::OVERWRITE || aexp == PrintSituation::OPEN_NEW) {
      initializeTrajectory(&foutp, output_kind, atom_count, current_time);
    }
    break;
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_ASCII_RST:
    initializeTrajectory(&foutp, output_kind, atom_count, current_time);
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("Coordinate file format unspecified.", "PhaseSpace", "exportToFile");
  }
  switch (output_kind) {
  case CoordinateFileKind::AMBER_CRD:
  case CoordinateFileKind::AMBER_INPCRD:
    writeFrame(&foutp, file_name, output_kind, atom_count,
               getCoordinatePointer(CartesianDimension::X, traj_kind),
               getCoordinatePointer(CartesianDimension::Y, traj_kind),
               getCoordinatePointer(CartesianDimension::Z, traj_kind), nullptr, nullptr, nullptr,
               unit_cell, box_dimensions.data());
    break;
  case CoordinateFileKind::AMBER_ASCII_RST:
    writeFrame(&foutp, file_name, output_kind, atom_count,
               getCoordinatePointer(CartesianDimension::X, TrajectoryKind::POSITIONS),
               getCoordinatePointer(CartesianDimension::Y, TrajectoryKind::POSITIONS),
               getCoordinatePointer(CartesianDimension::Z, TrajectoryKind::POSITIONS),
               getCoordinatePointer(CartesianDimension::X, TrajectoryKind::VELOCITIES),
               getCoordinatePointer(CartesianDimension::Y, TrajectoryKind::VELOCITIES),
               getCoordinatePointer(CartesianDimension::Z, TrajectoryKind::VELOCITIES),
               unit_cell, box_dimensions.data());
    break;
  case CoordinateFileKind::SDF:
    rtErr("The object does not have sufficient information to create an annotated SD file.  The "
          "program must use one of the writeFrame() overloads from the write_annotated_frame "
          "library instead.", "PhaseSpace", "exportToFile");
    break;
  case CoordinateFileKind::AMBER_NETCDF:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    break;
  case CoordinateFileKind::UNKNOWN:
    break;
  }
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
const PhaseSpace* PhaseSpace::getSelfPointer() const {
  return this;
}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceReader PhaseSpace::data(const HybridTargetLevel tier) const {
  return data(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceReader PhaseSpace::data(const CoordinateCycle orientation,
                                        const HybridTargetLevel tier) const {
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return PhaseSpaceReader(atom_count, unit_cell, x_alt_coordinates.data(tier),
                            y_alt_coordinates.data(tier), z_alt_coordinates.data(tier),
                            alt_box_space_transform.data(tier),
                            alt_inverse_transform.data(tier), alt_box_dimensions.data(tier),
                            box_space_transform.data(tier), inverse_transform.data(tier),
                            box_dimensions.data(tier), x_alt_velocities.data(tier),
                            y_alt_velocities.data(tier), z_alt_velocities.data(tier),
                            x_alt_forces.data(tier), y_alt_forces.data(tier),
                            z_alt_forces.data(tier), x_coordinates.data(tier),
                            y_coordinates.data(tier), z_coordinates.data(tier),
                            x_velocities.data(tier), y_velocities.data(tier),
                            z_velocities.data(tier), x_forces.data(tier), y_forces.data(tier),
                            z_forces.data(tier));
  case CoordinateCycle::WHITE:
    return PhaseSpaceReader(atom_count, unit_cell, x_coordinates.data(tier),
                            y_coordinates.data(tier), z_coordinates.data(tier),
                            box_space_transform.data(tier), inverse_transform.data(tier),
                            box_dimensions.data(tier), alt_box_space_transform.data(tier),
                            alt_inverse_transform.data(tier), alt_box_dimensions.data(tier),
                            x_velocities.data(tier), y_velocities.data(tier),
                            z_velocities.data(tier), x_forces.data(tier), y_forces.data(tier),
                            z_forces.data(tier), x_alt_coordinates.data(tier),
                            y_alt_coordinates.data(tier), z_alt_coordinates.data(tier),
                            x_alt_velocities.data(tier), y_alt_velocities.data(tier),
                            z_alt_velocities.data(tier), x_alt_forces.data(tier),
                            y_alt_forces.data(tier), z_alt_forces.data(tier));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceWriter PhaseSpace::data(const HybridTargetLevel tier) {
  return data(cycle_position, tier);
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceWriter PhaseSpace::data(const CoordinateCycle orientation,
                                  const HybridTargetLevel tier) {
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return PhaseSpaceWriter(atom_count, unit_cell, x_alt_coordinates.data(tier),
                            y_alt_coordinates.data(tier), z_alt_coordinates.data(tier),
                            alt_box_space_transform.data(tier),
                            alt_inverse_transform.data(tier), alt_box_dimensions.data(tier),
                            box_space_transform.data(tier), inverse_transform.data(tier),
                            box_dimensions.data(tier), x_alt_velocities.data(tier),
                            y_alt_velocities.data(tier), z_alt_velocities.data(tier),
                            x_alt_forces.data(tier), y_alt_forces.data(tier),
                            z_alt_forces.data(tier), x_coordinates.data(tier),
                            y_coordinates.data(tier), z_coordinates.data(tier),
                            x_velocities.data(tier), y_velocities.data(tier),
                            z_velocities.data(tier), x_forces.data(tier), y_forces.data(tier),
                            z_forces.data(tier));
  case CoordinateCycle::WHITE:
    return PhaseSpaceWriter(atom_count, unit_cell, x_coordinates.data(tier),
                            y_coordinates.data(tier), z_coordinates.data(tier),
                            box_space_transform.data(tier), inverse_transform.data(tier),
                            box_dimensions.data(tier), alt_box_space_transform.data(tier),
                            alt_inverse_transform.data(tier), alt_box_dimensions.data(tier),
                            x_velocities.data(tier), y_velocities.data(tier),
                            z_velocities.data(tier), x_forces.data(tier), y_forces.data(tier),
                            z_forces.data(tier), x_alt_coordinates.data(tier),
                            y_alt_coordinates.data(tier), z_alt_coordinates.data(tier),
                            x_alt_velocities.data(tier), y_alt_velocities.data(tier),
                            z_alt_velocities.data(tier), x_alt_forces.data(tier),
                            y_alt_forces.data(tier), z_alt_forces.data(tier));
  }
  __builtin_unreachable();
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
const PhaseSpaceReader PhaseSpace::deviceViewToHostData(CoordinateCycle orientation) const {
  const double* xcrd = x_coordinates.getDeviceValidHostPointer();
  const double* ycrd = y_coordinates.getDeviceValidHostPointer();
  const double* zcrd = z_coordinates.getDeviceValidHostPointer();
  const double* xcrd_alt = x_alt_coordinates.getDeviceValidHostPointer();
  const double* ycrd_alt = y_alt_coordinates.getDeviceValidHostPointer();
  const double* zcrd_alt = z_alt_coordinates.getDeviceValidHostPointer();
  const double* umat = box_space_transform.getDeviceValidHostPointer();
  const double* invu = inverse_transform.getDeviceValidHostPointer();
  const double* boxdim = box_dimensions.getDeviceValidHostPointer();
  const double* umat_alt = alt_box_space_transform.getDeviceValidHostPointer();
  const double* invu_alt = alt_inverse_transform.getDeviceValidHostPointer();
  const double* boxdim_alt = alt_box_dimensions.getDeviceValidHostPointer();
  const double* xvel = x_velocities.getDeviceValidHostPointer();
  const double* yvel = y_velocities.getDeviceValidHostPointer();
  const double* zvel = z_velocities.getDeviceValidHostPointer();
  const double* xvel_alt = x_alt_velocities.getDeviceValidHostPointer();
  const double* yvel_alt = y_alt_velocities.getDeviceValidHostPointer();
  const double* zvel_alt = z_alt_velocities.getDeviceValidHostPointer();
  const double* xfrc = x_forces.getDeviceValidHostPointer();
  const double* yfrc = y_forces.getDeviceValidHostPointer();
  const double* zfrc = z_forces.getDeviceValidHostPointer();
  const double* xfrc_alt = x_alt_forces.getDeviceValidHostPointer();
  const double* yfrc_alt = y_alt_forces.getDeviceValidHostPointer();
  const double* zfrc_alt = z_alt_forces.getDeviceValidHostPointer();
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return PhaseSpaceReader(atom_count, unit_cell, xcrd_alt, ycrd_alt, zcrd_alt, umat_alt,
                            invu_alt, boxdim_alt, umat, invu, boxdim, xvel_alt, yvel_alt, zvel_alt,
                            xfrc_alt, yfrc_alt, zfrc_alt, xcrd, ycrd, zcrd, xvel, yvel, zvel, xfrc,
                            yfrc, zfrc);
  case CoordinateCycle::WHITE:
    return PhaseSpaceReader(atom_count, unit_cell, xcrd, ycrd, zcrd, umat, invu, boxdim, umat_alt,
                            invu_alt, boxdim_alt, xvel, yvel, zvel, xfrc, yfrc, zfrc, xcrd_alt,
                            ycrd_alt, zcrd_alt, xvel_alt, yvel_alt, zvel_alt, xfrc_alt, yfrc_alt,
                            zfrc_alt);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceReader PhaseSpace::deviceViewToHostData() const {
  return deviceViewToHostData(cycle_position);
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceWriter PhaseSpace::deviceViewToHostData(CoordinateCycle orientation) {
  double* xcrd = x_coordinates.getDeviceValidHostPointer();
  double* ycrd = y_coordinates.getDeviceValidHostPointer();
  double* zcrd = z_coordinates.getDeviceValidHostPointer();
  double* xcrd_alt = x_alt_coordinates.getDeviceValidHostPointer();
  double* ycrd_alt = y_alt_coordinates.getDeviceValidHostPointer();
  double* zcrd_alt = z_alt_coordinates.getDeviceValidHostPointer();
  double* umat = box_space_transform.getDeviceValidHostPointer();
  double* invu = inverse_transform.getDeviceValidHostPointer();
  double* boxdim = box_dimensions.getDeviceValidHostPointer();
  double* umat_alt = alt_box_space_transform.getDeviceValidHostPointer();
  double* invu_alt = alt_inverse_transform.getDeviceValidHostPointer();
  double* boxdim_alt = alt_box_dimensions.getDeviceValidHostPointer();
  double* xvel = x_velocities.getDeviceValidHostPointer();
  double* yvel = y_velocities.getDeviceValidHostPointer();
  double* zvel = z_velocities.getDeviceValidHostPointer();
  double* xvel_alt = x_alt_velocities.getDeviceValidHostPointer();
  double* yvel_alt = y_alt_velocities.getDeviceValidHostPointer();
  double* zvel_alt = z_alt_velocities.getDeviceValidHostPointer();
  double* xfrc = x_forces.getDeviceValidHostPointer();
  double* yfrc = y_forces.getDeviceValidHostPointer();
  double* zfrc = z_forces.getDeviceValidHostPointer();
  double* xfrc_alt = x_alt_forces.getDeviceValidHostPointer();
  double* yfrc_alt = y_alt_forces.getDeviceValidHostPointer();
  double* zfrc_alt = z_alt_forces.getDeviceValidHostPointer();
  switch (orientation) {
  case CoordinateCycle::BLACK:
    return PhaseSpaceWriter(atom_count, unit_cell, xcrd_alt, ycrd_alt, zcrd_alt, umat_alt,
                            invu_alt, boxdim_alt, umat, invu, boxdim, xvel_alt, yvel_alt, zvel_alt,
                            xfrc_alt, yfrc_alt, zfrc_alt, xcrd, ycrd, zcrd, xvel, yvel, zvel, xfrc,
                            yfrc, zfrc);
  case CoordinateCycle::WHITE:
    return PhaseSpaceWriter(atom_count, unit_cell, xcrd, ycrd, zcrd, umat, invu, boxdim, umat_alt,
                            invu_alt, boxdim_alt, xvel, yvel, zvel, xfrc, yfrc, zfrc, xcrd_alt,
                            ycrd_alt, zcrd_alt, xvel_alt, yvel_alt, zvel_alt, xfrc_alt, yfrc_alt,
                            zfrc_alt);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceWriter PhaseSpace::deviceViewToHostData() {
  return deviceViewToHostData(cycle_position);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::upload() {
  storage.upload();
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::uploadPositions(const CoordinateCycle orientation) {
  const int stride = roundUp(atom_count, warp_size_int);
  switch (orientation) {
  case CoordinateCycle::BLACK:

    // Upload based on the storage array, as the individual X, Y, or Z velocity arrays will not
    // permit uploading 3x their stated length, even though the extra memory is there.
    storage.upload(3 * stride, 3 * stride);
    break;
  case CoordinateCycle::WHITE:
    storage.upload(0         , 3 * stride);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::uploadPositions() {
  uploadPositions(cycle_position);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::uploadTransformations() {
  const int stride = roundUp(atom_count, warp_size_int);
  storage.upload(9 * stride, 3 * roundUp(9, warp_size_int));
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::uploadVelocities(const CoordinateCycle orientation) {
  const int stride = roundUp(atom_count, warp_size_int);
  const int offset = 3 * roundUp(9, warp_size_int);
  switch (orientation) {
  case CoordinateCycle::BLACK:
    storage.upload((9 * stride) + offset, 3 * stride);
    break;
  case CoordinateCycle::WHITE:
    storage.upload((6 * stride) + offset, 3 * stride);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::uploadVelocities() {
  uploadVelocities(cycle_position);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::uploadForces(const CoordinateCycle orientation) {
  const int stride = roundUp(atom_count, warp_size_int);
  const int offset = 3 * roundUp(9, warp_size_int);
  switch (orientation) {
  case CoordinateCycle::BLACK:
    storage.upload((15 * stride) + offset, 3 * stride);
    break;
  case CoordinateCycle::WHITE:
    storage.upload((12 * stride) + offset, 3 * stride);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::uploadForces() {
  uploadForces(cycle_position);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::download() {
  storage.download();
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::downloadPositions(const CoordinateCycle orientation) {
  const int stride = roundUp(atom_count, warp_size_int);
  switch (orientation) {
  case CoordinateCycle::BLACK:
    storage.download(3 * stride, 3 * stride);
    break;
  case CoordinateCycle::WHITE:
    storage.download(0         , 3 * stride);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::downloadPositions() {
  downloadPositions(cycle_position);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::downloadTransformations() {
  const int stride = roundUp(atom_count, warp_size_int);
  storage.download(9 * stride, 3 * roundUp(9, warp_size_int));
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::downloadVelocities(const CoordinateCycle orientation) {
  const int stride = roundUp(atom_count, warp_size_int);
  const int offset = 3 * roundUp(9, warp_size_int);
  switch (orientation) {
  case CoordinateCycle::BLACK:
    storage.download(( 9 * stride) + offset, 3 * stride);
    break;
  case CoordinateCycle::WHITE:
    storage.download(( 6 * stride) + offset, 3 * stride);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::downloadVelocities() {
  downloadVelocities(cycle_position);
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::downloadForces(const CoordinateCycle orientation) {
  const int stride = roundUp(atom_count, warp_size_int);
  const int offset = 3 * roundUp(9, warp_size_int);
  switch (orientation) {
  case CoordinateCycle::BLACK:
    storage.download((15 * stride) + offset, 3 * stride);
    break;
  case CoordinateCycle::WHITE:
    storage.download((12 * stride) + offset, 3 * stride);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void PhaseSpace::downloadForces() {
  downloadForces(cycle_position);
}
#endif

//-------------------------------------------------------------------------------------------------
void PhaseSpace::allocate() {
  const int padded_atom_count  = roundUp(atom_count, warp_size_int);
  const int padded_matrix_size = roundUp(9, warp_size_int);
  storage.resize((27 * padded_atom_count) + (6 * padded_matrix_size));
  x_coordinates.setPointer(&storage,                           0, atom_count);
  y_coordinates.setPointer(&storage,           padded_atom_count, atom_count);
  z_coordinates.setPointer(&storage,       2 * padded_atom_count, atom_count);
  x_alt_coordinates.setPointer(&storage,   3 * padded_atom_count, atom_count);
  y_alt_coordinates.setPointer(&storage,   4 * padded_atom_count, atom_count);
  z_alt_coordinates.setPointer(&storage,   5 * padded_atom_count, atom_count);
  box_space_transform.setPointer(&storage, 6 * padded_atom_count, 9);
  inverse_transform.setPointer(&storage,  (6 * padded_atom_count) +      padded_matrix_size,  9);
  box_dimensions.setPointer(&storage,     (6 * padded_atom_count) + (2 * padded_matrix_size), 6);
  alt_box_space_transform.setPointer(&storage,
                                     (6 * padded_atom_count) + (3 * padded_matrix_size), 9);
  alt_inverse_transform.setPointer(&storage,
                                   (6 * padded_atom_count) + (4 * padded_matrix_size), 9);
  alt_box_dimensions.setPointer(&storage, (6 * padded_atom_count) + (5 * padded_matrix_size), 6);
  const int thus_far = (6 * padded_atom_count) + (6 * padded_matrix_size);
  x_velocities.setPointer(&storage,      thus_far,                            atom_count);
  y_velocities.setPointer(&storage,      thus_far +       padded_atom_count,  atom_count);
  z_velocities.setPointer(&storage,      thus_far +  (2 * padded_atom_count), atom_count);
  x_alt_velocities.setPointer(&storage,  thus_far +  (3 * padded_atom_count), atom_count);
  y_alt_velocities.setPointer(&storage,  thus_far +  (4 * padded_atom_count), atom_count);
  z_alt_velocities.setPointer(&storage,  thus_far +  (5 * padded_atom_count), atom_count);
  x_forces.setPointer(&storage,          thus_far +  (6 * padded_atom_count), atom_count);
  y_forces.setPointer(&storage,          thus_far +  (7 * padded_atom_count), atom_count);
  z_forces.setPointer(&storage,          thus_far +  (8 * padded_atom_count), atom_count);
  x_alt_forces.setPointer(&storage,      thus_far +  (9 * padded_atom_count), atom_count);
  y_alt_forces.setPointer(&storage,      thus_far + (10 * padded_atom_count), atom_count);
  z_alt_forces.setPointer(&storage,      thus_far + (11 * padded_atom_count), atom_count);
}

//-------------------------------------------------------------------------------------------------
UnitCellType determineUnitCellTypeByShape(const double* inverse_transform) {
  if (inverse_transform[0] < 1.01 && inverse_transform[4] < 1.01 && inverse_transform[8] < 1.01) {

    // No unit cell could be this small, and the identity matrix is what the transforms are set to
    // when there are no boundary conditions.  It is an inocuous matrix that can be put through
    // various operations without leading to catastrophes, although code should check the unit_cell
    // member variable to avoid meaningless calculations.
    return UnitCellType::NONE;
  }
  else if (fabs(inverse_transform[1]) + fabs(inverse_transform[2]) +
           fabs(inverse_transform[3]) + fabs(inverse_transform[5]) +
           fabs(inverse_transform[6]) + fabs(inverse_transform[7]) < constants::small) {
    return UnitCellType::ORTHORHOMBIC;
  }
  else {
    return UnitCellType::TRICLINIC;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
UnitCellType determineUnitCellTypeByShape(const std::vector<double> &inverse_transform) {
  return determineUnitCellTypeByShape(inverse_transform.data());
}

//-------------------------------------------------------------------------------------------------
UnitCellType determineUnitCellTypeByShape(const Hybrid<double> &inverse_transform) {
  return determineUnitCellTypeByShape(inverse_transform.data());
}

//-------------------------------------------------------------------------------------------------
std::vector<double> interlaceXYZ(const double* xptr, const double* yptr, const double* zptr,
                                 const int low_index, const int high_index) {
  std::vector<double> result(3 * (high_index - low_index));
  for (int i = low_index; i < high_index; i++) {
    const int base_idx = 3 * (i - low_index);
    result[base_idx    ] = xptr[i];
    result[base_idx + 1] = yptr[i];
    result[base_idx + 2] = zptr[i];
  }
  return result;
}

} // namespace trajectory
} // namespace stormm
