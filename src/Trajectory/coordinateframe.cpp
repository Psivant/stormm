#include "copyright.h"
#include "Accelerator/gpu_enumerators.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "FileManagement/file_listing.h"
#include "FileManagement/file_util.h"
#include "Math/matrix_ops.h"
#include "Math/rounding.h"
#include "MoleculeFormat/molecule_file_io.h"
#include "amber_ascii.h"
#include "coordinateframe.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

using card::confirmCpuMemory;
using card::confirmGpuMemory;
using card::getEnumerationName;
using constants::CartesianDimension;
using diskutil::DataFormat;
using diskutil::detectCoordinateFileKind;
using diskutil::getTrajectoryFormat;
using diskutil::DrivePathType;
using diskutil::getDrivePathType;
using stmath::extractBoxDimensions;
using stmath::roundUp;
using parse::TextFileReader;
using structure::extractSdfCoordinates;

//-------------------------------------------------------------------------------------------------
CoordinateFrameWriter::CoordinateFrameWriter(const int natom_in, const UnitCellType unit_cell_in,
                                             double* xcrd_in, double* ycrd_in, double* zcrd_in,
                                             double* umat_in, double* invu_in, double* boxdim_in) :
    natom{natom_in}, unit_cell{unit_cell_in}, xcrd{xcrd_in}, ycrd{ycrd_in}, zcrd{zcrd_in},
    umat{umat_in}, invu{invu_in}, boxdim{boxdim_in}
{}

//-------------------------------------------------------------------------------------------------
CoordinateFrameWriter::CoordinateFrameWriter(PhaseSpace *ps, const HybridTargetLevel tier) :
    natom{ps->getAtomCount()},
    unit_cell{ps->getUnitCellType()},
    xcrd{ps->getCoordinatePointer(CartesianDimension::X, TrajectoryKind::POSITIONS, tier)},
    ycrd{ps->getCoordinatePointer(CartesianDimension::Y, TrajectoryKind::POSITIONS, tier)},
    zcrd{ps->getCoordinatePointer(CartesianDimension::Z, TrajectoryKind::POSITIONS, tier)},
    umat{ps->getBoxSpaceTransformPointer(tier)},
    invu{ps->getInverseTransformPointer(tier)},
    boxdim{ps->getBoxSizePointer(tier)}
{}

//-------------------------------------------------------------------------------------------------
CoordinateFrameWriter::CoordinateFrameWriter(PhaseSpace *ps, const TrajectoryKind kind,
                                             const CoordinateCycle orientation,
                                             const HybridTargetLevel tier) :
    natom{ps->getAtomCount()},
    unit_cell{ps->getUnitCellType()},
    xcrd{}, ycrd{}, zcrd{},
    umat{ps->getBoxSpaceTransformPointer(tier)},
    invu{ps->getInverseTransformPointer(tier)},
    boxdim{ps->getBoxSizePointer(tier)}
{
  PhaseSpaceWriter psw = ps->data(orientation);

  // With the PhaseSpace abstract having been taken at the specified point in the time cycle,
  // the CoordinateFrameWriter can be set based on whatever appears in the current position,
  // velocity, or force arrays.
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    xcrd = psw.xcrd;
    ycrd = psw.ycrd;
    zcrd = psw.zcrd;
    break;
  case TrajectoryKind::VELOCITIES:
    xcrd = psw.xvel;
    ycrd = psw.yvel;
    zcrd = psw.zvel;
    break;
  case TrajectoryKind::FORCES:
    xcrd = psw.xfrc;
    ycrd = psw.yfrc;
    zcrd = psw.zfrc;
    break;
  }
}

//-------------------------------------------------------------------------------------------------
CoordinateFrameWriter::CoordinateFrameWriter(PhaseSpace *ps, const TrajectoryKind kind,
                                             const HybridTargetLevel tier) :
    CoordinateFrameWriter(ps, kind, ps->getCyclePosition(), tier)
{}

//-------------------------------------------------------------------------------------------------
CoordinateFrameReader::CoordinateFrameReader(const int natom_in, const UnitCellType unit_cell_in,
                                             const double* xcrd_in, const double* ycrd_in,
                                             const double* zcrd_in, const double* umat_in,
                                             const double* invu_in, const double* boxdim_in) :
    natom{natom_in}, unit_cell{unit_cell_in}, xcrd{xcrd_in}, ycrd{ycrd_in}, zcrd{zcrd_in},
    umat{umat_in}, invu{invu_in}, boxdim{boxdim_in}
{}

//-------------------------------------------------------------------------------------------------
CoordinateFrameReader::CoordinateFrameReader(const CoordinateFrameWriter &cfw) :
    natom{cfw.natom}, unit_cell{cfw.unit_cell}, xcrd{cfw.xcrd}, ycrd{cfw.ycrd}, zcrd{cfw.zcrd},
    umat{cfw.umat}, invu{cfw.invu}, boxdim{cfw.boxdim}
{}

//-------------------------------------------------------------------------------------------------
CoordinateFrameReader::CoordinateFrameReader(const PhaseSpace *ps, const HybridTargetLevel tier) :
    natom{ps->getAtomCount()},
    unit_cell{ps->getUnitCellType()},
    xcrd{ps->getCoordinatePointer(CartesianDimension::X, TrajectoryKind::POSITIONS, tier)},
    ycrd{ps->getCoordinatePointer(CartesianDimension::Y, TrajectoryKind::POSITIONS, tier)},
    zcrd{ps->getCoordinatePointer(CartesianDimension::Z, TrajectoryKind::POSITIONS, tier)},
    umat{ps->getBoxSpaceTransformPointer(tier)},
    invu{ps->getInverseTransformPointer(tier)},
    boxdim{ps->getBoxSizePointer(tier)}
{}

//-------------------------------------------------------------------------------------------------
CoordinateFrameReader::CoordinateFrameReader(const PhaseSpace &ps, const HybridTargetLevel tier) :
    CoordinateFrameReader(ps.getSelfPointer(), tier)
{}

//-------------------------------------------------------------------------------------------------
CoordinateFrameReader::CoordinateFrameReader(const PhaseSpace *ps, const TrajectoryKind kind,
                                             const CoordinateCycle orientation,
                                             const HybridTargetLevel tier) :
    natom{ps->getAtomCount()},
    unit_cell{ps->getUnitCellType()},
    xcrd{}, ycrd{}, zcrd{},
    umat{ps->getBoxSpaceTransformPointer(tier)},
    invu{ps->getInverseTransformPointer(tier)},
    boxdim{ps->getBoxSizePointer(tier)}
{
  const PhaseSpaceReader psr = ps->data(orientation);

  // With the PhaseSpace abstract having been taken at the specified point in the time cycle,
  // the CoordinateFrameWriter can be set based on whatever appears in the current position,
  // velocity, or force arrays.
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    xcrd = psr.xcrd;
    ycrd = psr.ycrd;
    zcrd = psr.zcrd;
    break;
  case TrajectoryKind::VELOCITIES:
    xcrd = psr.xvel;
    ycrd = psr.yvel;
    zcrd = psr.zvel;
    break;
  case TrajectoryKind::FORCES:
    xcrd = psr.xfrc;
    ycrd = psr.yfrc;
    zcrd = psr.zfrc;
    break;
  }
}

//-------------------------------------------------------------------------------------------------
CoordinateFrameReader::CoordinateFrameReader(const PhaseSpace &ps, const TrajectoryKind kind,
                                             const CoordinateCycle orientation,
                                             const HybridTargetLevel tier) :
    CoordinateFrameReader(ps.getSelfPointer(), kind, ps.getCyclePosition(), tier)
{}

//-------------------------------------------------------------------------------------------------
CoordinateFrameReader::CoordinateFrameReader(const PhaseSpace *ps, const TrajectoryKind kind,
                                             const HybridTargetLevel tier) :
    CoordinateFrameReader(ps, kind, ps->getCyclePosition(), tier)
{}

//-------------------------------------------------------------------------------------------------
CoordinateFrameReader::CoordinateFrameReader(const PhaseSpace &ps, const TrajectoryKind kind,
                                             const HybridTargetLevel tier) :
    CoordinateFrameReader(ps.getSelfPointer(), kind, ps.getCyclePosition(), tier)
{}

//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(const int natom_in, const UnitCellType unit_cell_in,
                                 const HybridFormat format_in) :
    format{format_in},
    file_name{std::string("")},
    frame_number{0},
    atom_count{natom_in},
    unit_cell{unit_cell_in},
    x_coordinates{HybridKind::POINTER, "x_coordinates", format_in},
    y_coordinates{HybridKind::POINTER, "y_coordinates", format_in},
    z_coordinates{HybridKind::POINTER, "z_coordinates", format_in},
    box_space_transform{HybridKind::POINTER, "box_transform", format_in},
    inverse_transform{HybridKind::POINTER, "inv_transform", format_in},
    box_dimensions{HybridKind::POINTER, "box_dimensions", format_in},
    storage{HybridKind::ARRAY, "frame_data", format_in}
{
  allocate();
}
  
//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(const int natom_in, const UnitCellType unit_cell_in,
                                 const double* xcrd_in, const double* ycrd_in,
                                 const double* zcrd_in, const double* umat_in,
                                 const double* invu_in, const double* boxdim_in,
                                 const HybridFormat format_in) :
    CoordinateFrame(natom_in, unit_cell_in, format_in)
{
  switch (format) {
#ifdef STORMM_USE_HPC
  case HybridFormat::DEVICE_ONLY:
    {
      // In much the same manner as files can be transferred into device-exclusive objects (see
      // below), a host-exclusive temporary object is used to stage the data.
      const CoordinateFrame tmp(natom_in, unit_cell_in, xcrd_in, ycrd_in, zcrd_in, umat_in,
                                invu_in, boxdim_in, HybridFormat::HOST_ONLY);
      deepCopy(&storage, tmp.storage);
    }
    break;
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
#endif
  case HybridFormat::HOST_ONLY:
    
    // No method yet for putting data into a Hybrid object based on a pointer to an array
    double* xtmp = x_coordinates.data();
    double* ytmp = y_coordinates.data();
    double* ztmp = z_coordinates.data();
    for (int i = 0; i < natom_in; i++) {
      xtmp[i] = xcrd_in[i];
      ytmp[i] = ycrd_in[i];
      ztmp[i] = zcrd_in[i];
    }

    // Create transformation matrices if available, or initialize to the identity matrix
    double* utmp = box_space_transform.data();
    double* invutmp = inverse_transform.data();
    double* bdimtmp = box_dimensions.data();
    if (umat_in != nullptr) {
      for (int i = 0; i < 9; i++) {
        utmp[i] = umat_in[i];
        invutmp[i] = invu_in[i];
      }

      // Extract the box dimensions directly from the inverse transformation matrix if necessary
      if (boxdim_in == nullptr) {
        extractBoxDimensions(&bdimtmp[0], &bdimtmp[1], &bdimtmp[2], &bdimtmp[3], &bdimtmp[4],
                             &bdimtmp[5], invu_in);
      }
      else {
        for (int i = 0; i < 6; i++) {
          bdimtmp[i] = boxdim_in[i];
        }
      }
    }
    else {
      for (int i = 0; i < 9; i++) {
        utmp[i]    = static_cast<double>((i & 0x3) == 0);
        invutmp[i] = static_cast<double>((i & 0x3) == 0);
      }
      for (int i = 0; i < 3; i++) {
        bdimtmp[i    ] = 0.0;
        bdimtmp[i + 3] = 0.5 * symbols::pi;
      }
    }
    break;
  }
#ifdef STORMM_USE_HPC
  upload();
#endif
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(const std::string &file_name_in,
                                 const CoordinateFileKind file_kind, const int frame_number_in,
                                 const HybridFormat format_in) :
    CoordinateFrame(0, UnitCellType::NONE, format_in)
{
#ifdef STORMM_USE_HPC
  switch (format) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
    buildFromFile(file_name_in, file_kind, frame_number_in);
    break;
  case HybridFormat::DEVICE_ONLY:
    {
      // If a file is to be read but there is no host data to receive the information, make a
      // temporary object with basic host memory to expedite a transfer of the data to the device.
      const CoordinateFrame tmp(file_name_in, file_kind, frame_number_in, HybridFormat::HOST_ONLY);
      deepCopy(&storage, tmp.storage);
    }
    break;
  }
  upload();
#else
  buildFromFile(file_name_in, file_kind, frame_number_in);
#endif
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(const TextFile &tf, const CoordinateFileKind file_kind,
                                 const int frame_number_in, const HybridFormat format_in) :
    CoordinateFrame(0, UnitCellType::NONE, format_in)
{
#ifdef STORMM_USE_HPC
  switch (format) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
    buildFromFile(tf, file_kind, frame_number_in);
    break;
  case HybridFormat::DEVICE_ONLY:
    {
      // If a file is to be read but there is no host data to receive the information, make a
      // temporary object with basic host memory to expedite a transfer of the data to the device.
      const CoordinateFrame tmp(tf, file_kind, frame_number_in, HybridFormat::HOST_ONLY);
      deepCopy(&storage, tmp.storage);
    }
    break;
  }
  upload();
#else
  buildFromFile(tf, file_kind, frame_number_in);
#endif
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(const PhaseSpace *ps) :
    CoordinateFrame(ps->getAtomCount(), ps->getUnitCellType(), ps->getFormat())
{
  // Transfer host-side data
  switch (format) {
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
#endif
  case HybridFormat::HOST_ONLY:
    {
      file_name = ps->getFileName();
      allocate();
      const PhaseSpaceReader psr = ps->data(HybridTargetLevel::HOST);
      double* xptr = x_coordinates.data();
      double* yptr = y_coordinates.data();
      double* zptr = z_coordinates.data();
      for (int i = 0; i < atom_count; i++) {
        xptr[i] = psr.xcrd[i];
        yptr[i] = psr.ycrd[i];
        zptr[i] = psr.zcrd[i];
      }
      double* utmp = box_space_transform.data();
      double* invutmp = inverse_transform.data();
      double* bdimtmp = box_dimensions.data();
      for (int i = 0; i < 9; i++) {
        utmp[i]    = psr.umat[i];
        invutmp[i] = psr.invu[i];
      }
      for (int i = 0; i < 6; i++) {
        bdimtmp[i] = psr.boxdim[i];
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::DEVICE_ONLY:
    break;
#endif
  }

  // Transfer device-side data
#ifdef STORMM_USE_HPC
  switch (format) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::DEVICE_ONLY:
    {
      const PhaseSpaceReader psr = ps->data(HybridTargetLevel::DEVICE);
      double* xptr = x_coordinates.data(HybridTargetLevel::DEVICE);
      double* yptr = y_coordinates.data(HybridTargetLevel::DEVICE);
      double* zptr = z_coordinates.data(HybridTargetLevel::DEVICE);
      cudaMemcpyKind d_to_d = cudaMemcpyDeviceToDevice;
      if (cudaMemcpy(xptr, psr.xcrd, psr.natom * sizeof(double), d_to_d) != cudaSuccess) {
        rtErr("Error copying X particle positions from a PhaseSpace object into a new "
              "CoordinateFrame object with " + std::to_string(psr.natom) + " atoms.",
              "CoordinateFrame");
      }
      if (cudaMemcpy(yptr, psr.ycrd, psr.natom * sizeof(double), d_to_d) != cudaSuccess) {
        rtErr("Error copying Y particle positions from a PhaseSpace object into a new "
              "CoordinateFrame object with " + std::to_string(psr.natom) + " atoms.",
              "CoordinateFrame");
      }
      if (cudaMemcpy(zptr, psr.zcrd, psr.natom * sizeof(double), d_to_d) != cudaSuccess) {
        rtErr("Error copying Z particle positions from a PhaseSpace object into a new "
              "CoordinateFrame object with " + std::to_string(psr.natom) + " atoms.",
              "CoordinateFrame");
      }
      double* utmp = box_space_transform.data(HybridTargetLevel::DEVICE);
      double* invutmp = inverse_transform.data(HybridTargetLevel::DEVICE);
      double* bdimtmp = box_dimensions.data(HybridTargetLevel::DEVICE);
      if (cudaMemcpy(utmp, psr.umat, 9 * sizeof(double), d_to_d) != cudaSuccess) {
        rtErr("Error copying box transform matrix from a PhaseSpace object into a new "
              "CoordinateFrame object.", "CoordinateFrame");
      }
      if (cudaMemcpy(invutmp, psr.invu, 9 * sizeof(double), d_to_d) != cudaSuccess) {
        rtErr("Error copying inverse transform matrix from a PhaseSpace object into a new "
              "CoordinateFrame object.", "CoordinateFrame");
      }
      if (cudaMemcpy(bdimtmp, psr.boxdim, 6 * sizeof(double), d_to_d) != cudaSuccess) {
        rtErr("Error copying box dimensions from a PhaseSpace object into a new "
              "CoordinateFrame object.", "CoordinateFrame");
      }
    }
    break;
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_MOUNTED:
  case HybridFormat::HOST_ONLY:
    break;
  }
#endif
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(const PhaseSpace &ps) :
    CoordinateFrame(ps.getSelfPointer())
{}

//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(const CoordinateFrame &original) :
    format{original.format},
    file_name{original.file_name},
    atom_count{original.atom_count},
    unit_cell{original.unit_cell},
    x_coordinates{original.x_coordinates},
    y_coordinates{original.y_coordinates},
    z_coordinates{original.z_coordinates},
    box_space_transform{original.box_space_transform},
    inverse_transform{original.inverse_transform},
    box_dimensions{original.box_dimensions},
    storage{original.storage}
{
  allocate();
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(const CoordinateFrame &original, const HybridFormat format_in) :
    format{format_in},
    file_name{original.file_name},
    atom_count{original.atom_count},
    unit_cell{original.unit_cell},
    x_coordinates{HybridKind::POINTER, original.x_coordinates.getLabel().name, format_in},
    y_coordinates{HybridKind::POINTER, original.y_coordinates.getLabel().name, format_in},
    z_coordinates{HybridKind::POINTER, original.z_coordinates.getLabel().name, format_in},
    box_space_transform{HybridKind::POINTER, original.box_space_transform.getLabel().name,
                        format_in},
    inverse_transform{HybridKind::POINTER, original.inverse_transform.getLabel().name, format_in},
    box_dimensions{HybridKind::POINTER, original.box_dimensions.getLabel().name, format_in},
    storage{HybridKind::ARRAY, original.storage.getLabel().name, format_in}
{
  allocate();
  deepCopy(&storage, original.storage);
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame& CoordinateFrame::operator=(const CoordinateFrame &other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // Copy the format, file name (if applicable), and atom count
  format = other.format;
  file_name = other.file_name;
  atom_count = other.atom_count;
  unit_cell = other.unit_cell;

  // Copy the Hybrid objects to preserve tags and the proper kinds.  Copying the data storage
  // object carries over all of the other object's contents.
  x_coordinates = other.x_coordinates;
  y_coordinates = other.y_coordinates;
  z_coordinates = other.z_coordinates;
  box_space_transform = other.box_space_transform;
  inverse_transform = other.inverse_transform;
  box_dimensions = other.box_dimensions;
  storage = other.storage;

  // Allocating again sets the internal POINTER-kind Hybrid objects
  allocate();
  return *this;
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame::CoordinateFrame(CoordinateFrame &&original) :
    format{original.format},
    file_name{std::move(original.file_name)},
    atom_count{original.atom_count},
    unit_cell{original.unit_cell},
    x_coordinates{std::move(original.x_coordinates)},
    y_coordinates{std::move(original.y_coordinates)},
    z_coordinates{std::move(original.z_coordinates)},
    box_space_transform{std::move(original.box_space_transform)},
    inverse_transform{std::move(original.inverse_transform)},
    box_dimensions{std::move(original.box_dimensions)},
    storage{std::move(original.storage)}
{
  // No repair of the pointers is necessary, as they still point to valid data
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame& CoordinateFrame::operator=(CoordinateFrame &&other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  format = other.format;
  file_name = std::move(other.file_name);
  atom_count = other.atom_count;
  unit_cell = other.unit_cell;
  x_coordinates = std::move(other.x_coordinates);
  y_coordinates = std::move(other.y_coordinates);
  z_coordinates = std::move(other.z_coordinates);
  box_space_transform = std::move(other.box_space_transform);
  inverse_transform = std::move(other.inverse_transform);
  box_dimensions = std::move(other.box_dimensions);
  storage = std::move(other.storage);
  return *this;
}

//-------------------------------------------------------------------------------------------------
void CoordinateFrame::buildFromFile(const std::string &file_name_in,
                                    const CoordinateFileKind file_kind, const int frame_number) {
  confirmCpuMemory(format, "CPU memory must be present in the object in order to load data from "
                   "a file.  Current format: " + getEnumerationName(format) + ".",
                   "CoordinateFrame", "buildFromFile");
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
              "file.", "CoordinateFrame", "buildFromFile");
      }
      allocate();
      TextFile tf(file_name);
      readAmberCrdFormat(tf, &x_coordinates, &y_coordinates, &z_coordinates, unit_cell,
                         &box_space_transform, &inverse_transform, &box_dimensions, frame_number);
    }
    break;
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_ASCII_RST:
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
  case CoordinateFileKind::SDF:
    break;
  case CoordinateFileKind::AMBER_NETCDF:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("The file type of " + file_name + " could not be understood.", "CoordinateFrame",
          "buildFromFile");
  }
#ifdef STORMM_USE_HPC
  upload();
#endif
}

//-------------------------------------------------------------------------------------------------
void CoordinateFrame::buildFromFile(const TextFile &tf, const CoordinateFileKind file_kind,
                                    const int frame_number) {
  confirmCpuMemory(format, "CPU memory must be present in the object in order to load data from "
                   "a TextFile object.  Current format: " + getEnumerationName(format) + ".",
                   "CoordinateFrame", "buildFromFile");
  file_name = tf.getFileName();

  // Try to detect the file format if it is not already specified.  If it remains UNKNOWN, that
  // will ultimately lead to an error.
  CoordinateFileKind actual_kind = file_kind;
  if (file_kind == CoordinateFileKind::UNKNOWN) {
    actual_kind = detectCoordinateFileKind(tf);
  }  
  switch (actual_kind) {
  case CoordinateFileKind::AMBER_CRD:
    {
      // The number of atoms must be known a-priori in order to read from a .crd trajectory file.
      if (atom_count == 0) {
        rtErr("A number of atoms matching the trajectory must be known prior to reading a .crd "
              "file.", "CoordinateFrame", "buildFromFile");
      }
      allocate();
      readAmberCrdFormat(tf, &x_coordinates, &y_coordinates, &z_coordinates, unit_cell,
                         &box_space_transform, &inverse_transform, &box_dimensions, frame_number);
    }
    break;
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_ASCII_RST:
    {
      atom_count = getAmberRestartAtomCount(tf);
      allocate();
      getAmberInputCoordinates(tf, &x_coordinates, &y_coordinates, &z_coordinates,
                               &box_space_transform, &inverse_transform, &box_dimensions);

      // Interpret the box transformation, updating the unit cell type based on the file
      unit_cell = determineUnitCellTypeByShape(inverse_transform.data());
    }
    break;
  case CoordinateFileKind::SDF:
    {
      const std::vector<double3> mdl_crd = extractSdfCoordinates(tf, frame_number);

      // The SDF format is assumed to encode no box information
      unit_cell = UnitCellType::NONE;
      double* box_ptr = box_dimensions.data();
      double* umat_ptr = box_space_transform.data();
      double* invu_ptr = inverse_transform.data();
      for (int i = 0; i < 3; i++) {
        box_ptr[i] = 0.0;
      }
      for (int i = 0; i < 9; i++) {
        umat_ptr[i] = static_cast<double>((i & 0x3) == 0);
        invu_ptr[i] = static_cast<double>((i & 0x3) == 0);
      }
    }
    break;
  case CoordinateFileKind::AMBER_NETCDF:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("The file type of " + file_name + " could not be understood.", "CoordinateFrame",
          "buildFromFile");
  }
#ifdef STORMM_USE_HPC
  upload();
#endif
}

//-------------------------------------------------------------------------------------------------
HybridFormat CoordinateFrame::getFormat() const {
  return format;
}

//-------------------------------------------------------------------------------------------------
std::string CoordinateFrame::getFileName() const {
  return file_name;
}

//-------------------------------------------------------------------------------------------------
int CoordinateFrame::getFrameNumber() const {
  return frame_number;
}

//-------------------------------------------------------------------------------------------------
int CoordinateFrame::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
UnitCellType CoordinateFrame::getUnitCellType() const {
  return unit_cell;
}

//-------------------------------------------------------------------------------------------------
std::vector<double>
CoordinateFrame::getInterlacedCoordinates(const HybridTargetLevel tier) const {
  return getInterlacedCoordinates(0, atom_count, tier);
}

//-------------------------------------------------------------------------------------------------
std::vector<double>
CoordinateFrame::getInterlacedCoordinates(const int low_index, const int high_index,
                                          const HybridTargetLevel tier) const {
  
  // Range check as this will use pointers
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("Invalid atom range " + std::to_string(low_index) + " to " + std::to_string(high_index) +
          " in an object with " + std::to_string(atom_count) + " atoms.", "CoordinateFrame",
          "getInterlacedCoordinates");
  }
  checkFormatCompatibility(tier, format, "CoordinateFrame", "getInterlacedCoordinates");
  std::vector<double> result(3 * (high_index - low_index));
  switch (tier) {
  case HybridTargetLevel::HOST:
    return interlaceXYZ(x_coordinates.data(), y_coordinates.data(), z_coordinates.data(),
                        low_index, high_index);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<double> xval = x_coordinates.readDevice(low_index, high_index);
      const std::vector<double> yval = y_coordinates.readDevice(low_index, high_index);
      const std::vector<double> zval = z_coordinates.readDevice(low_index, high_index);
      return interlaceXYZ(xval.data(), yval.data(), zval.data(), low_index, high_index);
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> CoordinateFrame::getBoxSpaceTransform(const HybridTargetLevel tier) const {
  checkFormatCompatibility(tier, format, "CoordinateFrame", "getBoxSpaceTransform");
  switch (tier) {
  case HybridTargetLevel::HOST:
    return box_space_transform.readHost();
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return box_space_transform.readDevice();
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> CoordinateFrame::getInverseTransform(const HybridTargetLevel tier) const {
  checkFormatCompatibility(tier, format, "CoordinateFrame", "getInverseTransform");
  switch (tier) {
  case HybridTargetLevel::HOST:
    return inverse_transform.readHost();
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return inverse_transform.readDevice();
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> CoordinateFrame::getBoxDimensions(const HybridTargetLevel tier) const {
  checkFormatCompatibility(tier, format, "CoordinateFrame", "getBoxDimensions");
  switch (tier) {
  case HybridTargetLevel::HOST:
    return box_dimensions.readHost();
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return box_dimensions.readDevice();
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const Hybrid<double>* CoordinateFrame::getCoordinateHandle(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return &x_coordinates;
  case CartesianDimension::Y:
    return &y_coordinates;
  case CartesianDimension::Z:
    return &z_coordinates;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
Hybrid<double>* CoordinateFrame::getCoordinateHandle(const CartesianDimension dim) {
  switch (dim) {
  case CartesianDimension::X:
    return &x_coordinates;
  case CartesianDimension::Y:
    return &y_coordinates;
  case CartesianDimension::Z:
    return &z_coordinates;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const Hybrid<double>* CoordinateFrame::getBoxTransformHandle() const {
  return &box_space_transform;
}

//-------------------------------------------------------------------------------------------------
Hybrid<double>* CoordinateFrame::getBoxTransformHandle() {
  return &box_space_transform;
}

//-------------------------------------------------------------------------------------------------
const Hybrid<double>* CoordinateFrame::getInverseTransformHandle() const {
  return &inverse_transform;
}

//-------------------------------------------------------------------------------------------------
Hybrid<double>* CoordinateFrame::getInverseTransformHandle() {
  return &inverse_transform;
}

//-------------------------------------------------------------------------------------------------
const Hybrid<double>* CoordinateFrame::getBoxDimensionsHandle() const {
  return &box_dimensions;
}

//-------------------------------------------------------------------------------------------------
Hybrid<double>* CoordinateFrame::getBoxDimensionsHandle() {
  return &box_dimensions;
}

//-------------------------------------------------------------------------------------------------
const Hybrid<double>* CoordinateFrame::getStorageHandle() const {
  return &storage;
}

//-------------------------------------------------------------------------------------------------
Hybrid<double>* CoordinateFrame::getStorageHandle() {
  return &storage;
}

//-------------------------------------------------------------------------------------------------
const CoordinateFrame* CoordinateFrame::getSelfPointer() const {
  return this;
}

//-------------------------------------------------------------------------------------------------
void CoordinateFrame::exportToFile(const std::string &file_name, const CoordinateFileKind kind,
                                   const PrintSituation expectation,
                                   const HybridTargetLevel tier) const {
  const PrintSituation aexp = adjustTrajectoryOpeningProtocol(expectation, kind,
                                                              "CoordinateFrame", "exportToFile");
  const DataFormat style = getTrajectoryFormat(kind);
  const bool fi_exists = (getDrivePathType(file_name) == DrivePathType::FILE);
  std::ofstream foutp;
  foutp = openOutputFile(file_name, aexp, "Open an output file for writing CoordinateFrame "
                         "contents", style);
  switch (kind) {
  case CoordinateFileKind::AMBER_CRD:
  case CoordinateFileKind::SDF:
  case CoordinateFileKind::AMBER_NETCDF:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    if (fi_exists == false ||
        aexp == PrintSituation::OVERWRITE || aexp == PrintSituation::OPEN_NEW) {
      initializeTrajectory(&foutp, kind, atom_count);
    }
    break;
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_ASCII_RST:
    initializeTrajectory(&foutp, kind, atom_count);
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("Coordinate file format unspecified.", "CoordinateFrame", "exportToFile");
  }
  switch (kind) {
  case CoordinateFileKind::AMBER_CRD:
  case CoordinateFileKind::AMBER_INPCRD:
    switch (tier) {
    case HybridTargetLevel::HOST:
      confirmCpuMemory(format, "No memory is allocated on the CPU host (format " +
                       getEnumerationName(format) + ").", "CoordinateFrame", "exportToFile");
      writeFrame(&foutp, file_name, kind, atom_count, x_coordinates.data(), y_coordinates.data(),
                 z_coordinates.data(), nullptr, nullptr, nullptr, unit_cell,
                 box_dimensions.data());
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        confirmGpuMemory(format, "No memory is allocated on the GPU device (format " +
                         getEnumerationName(format) + ").", "CoordinateFrame", "exportToFile");
        const std::vector<double> xtmp = x_coordinates.readDevice();
        const std::vector<double> ytmp = y_coordinates.readDevice();
        const std::vector<double> ztmp = z_coordinates.readDevice();
        writeFrame(&foutp, file_name, kind, atom_count, xtmp.data(), ytmp.data(), ztmp.data(),
                   nullptr, nullptr, nullptr, unit_cell, box_dimensions.data());
      }
      break;
#endif
    }
    break;
  case CoordinateFileKind::SDF:
    rtErr("The object does not have sufficient information to create an annotated SD file.  The "
          "program must use one of the writeFrame() overloads from the write_annotated_frame "
          "library instead.", "CoordinateSeries", "exportToFile");
    break;
  case CoordinateFileKind::AMBER_NETCDF:
    break;
  case CoordinateFileKind::AMBER_ASCII_RST:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    rtErr("A restart file cannot be written based on a CoordinateFrame.  The object will not be "
          "able to store both coordinates and velocities needed for checkpointing.",
          "CoordinateFrame", "exportToFile");
    break;
  case CoordinateFileKind::UNKNOWN:
    break;
  }
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
const CoordinateFrameReader CoordinateFrame::data(HybridTargetLevel tier) const {
  checkFormatCompatibility(tier, format, "CoordinateFrame", "data");
  return CoordinateFrameReader(atom_count, unit_cell, x_coordinates.data(tier),
                               y_coordinates.data(tier), z_coordinates.data(tier),
                               box_space_transform.data(tier), inverse_transform.data(tier),
                               box_dimensions.data(tier));
}

//-------------------------------------------------------------------------------------------------
CoordinateFrameWriter CoordinateFrame::data(HybridTargetLevel tier) {
  checkFormatCompatibility(tier, format, "CoordinateFrame", "data");
  return CoordinateFrameWriter(atom_count, unit_cell, x_coordinates.data(tier),
                               y_coordinates.data(tier), z_coordinates.data(tier),
                               box_space_transform.data(tier), inverse_transform.data(tier),
                               box_dimensions.data(tier));
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void CoordinateFrame::upload() {
  storage.upload();
}

//-------------------------------------------------------------------------------------------------
void CoordinateFrame::download() {
  storage.download();
}

//-------------------------------------------------------------------------------------------------
const CoordinateFrameReader CoordinateFrame::deviceViewToHostData() const {
  confirmHostVisibleToGpu(format, "Host memory is not visible to the GPU in this object",
                          "CoordinateFrame", "deviceViewToHostData");
  const double* devc_xcrd = x_coordinates.getDeviceValidHostPointer();
  const double* devc_ycrd = y_coordinates.getDeviceValidHostPointer();
  const double* devc_zcrd = z_coordinates.getDeviceValidHostPointer();
  const double* devc_umat = box_space_transform.getDeviceValidHostPointer();
  const double* devc_invu = inverse_transform.getDeviceValidHostPointer();
  const double* devc_bdim = box_dimensions.getDeviceValidHostPointer();
  return CoordinateFrameReader(atom_count, unit_cell, devc_xcrd, devc_ycrd, devc_zcrd, devc_umat,
                               devc_invu, devc_bdim);
}

//-------------------------------------------------------------------------------------------------
CoordinateFrameWriter CoordinateFrame::deviceViewToHostData() {
  confirmHostVisibleToGpu(format, "Host memory is not visible to the GPU in this object",
                          "CoordinateFrame", "deviceViewToHostData");
  double* devc_xcrd = x_coordinates.getDeviceValidHostPointer();
  double* devc_ycrd = y_coordinates.getDeviceValidHostPointer();
  double* devc_zcrd = z_coordinates.getDeviceValidHostPointer();
  double* devc_umat = box_space_transform.getDeviceValidHostPointer();
  double* devc_invu = inverse_transform.getDeviceValidHostPointer();
  double* devc_bdim = box_dimensions.getDeviceValidHostPointer();
  return CoordinateFrameWriter(atom_count, unit_cell, devc_xcrd, devc_ycrd, devc_zcrd, devc_umat,
                               devc_invu, devc_bdim);
}
#endif

//-------------------------------------------------------------------------------------------------
void CoordinateFrame::setFrameNumber(int frame_number_in) {
  frame_number = frame_number_in;
}

//-------------------------------------------------------------------------------------------------
std::vector<CoordinateFrame> getSelectedFrames(const TextFile &tf, const CoordinateFileKind kind,
                                               const int atom_count, const UnitCellType unit_cell,
                                               const std::vector<int> &frame_numbers) {
  CoordinateFileKind actual_kind = kind;
  if (kind == CoordinateFileKind::UNKNOWN) {
    actual_kind = detectCoordinateFileKind(tf);
  }
  const int frame_count = frame_numbers.size();
  std::vector<CoordinateFrame> result(frame_count, CoordinateFrame(atom_count, unit_cell));
  switch (actual_kind) {
  case CoordinateFileKind::AMBER_CRD:
    for (int i = 0; i < frame_count; i++) {
      CoordinateFrameWriter cfw = result[i].data();
      readAmberCrdFormat(tf, &cfw, frame_numbers[i]);
    }
    break;
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_ASCII_RST:
    {
      // Check that the series of frames is one, the first frame
      for (int i = 0; i < frame_count; i++) {
        if (frame_numbers[i] != 0) {
          rtErr("An Amber inpcrd or restart file has only one frame.  Request for frame index " +
                std::to_string(frame_numbers[i]) + " (indexing starts from zero) is invalid.",
                "getSelectedFrames");
        }
      }

      // Since the atom count is available (and this is an unusual application of this function),
      // check that it matches the expected system size.
      int test_atom_count = getAmberRestartAtomCount(tf);
      if (test_atom_count != atom_count) {
        rtErr("When reading multiple frames from a trajectory, the atom count must be known in "
              "advance.  Reading a single frame from an input coordinates file is best done with "
              "one of the CoordinateFrame constructors.  The requested atom count of " +
              std::to_string(atom_count) + " does not agree with the detected atom count of " +
              std::to_string(test_atom_count) + ".", "getSelectedFrames");
      }
      const CoordinateFrame cftmp(tf, kind);
      for (int i = 0; i < frame_count; i++) {
        result[i] = cftmp;
      }
    }
    break;
  case CoordinateFileKind::SDF:
    break;
  case CoordinateFileKind::AMBER_NETCDF:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("The file type of " + tf.getFileName() + " could not be understood.",
          "getSelectedFrames");
  }
  return result;
}
                                               
//-------------------------------------------------------------------------------------------------
std::vector<CoordinateFrame> getSelectedFrames(const std::string &file_name, const int atom_count,
                                               const UnitCellType unit_cell,
                                               const std::vector<int> &frame_numbers) {
  const TextFile tf(file_name);
  const CoordinateFileKind kind = detectCoordinateFileKind(tf);
  return getSelectedFrames(tf, kind, atom_count, unit_cell, frame_numbers);
}

//-------------------------------------------------------------------------------------------------
std::vector<CoordinateFrame> getAllFrames(const TextFile &tf, const int atom_count,
                                          const UnitCellType unit_cell,
                                          const ExceptionResponse policy) {
  const CoordinateFileKind kind = detectCoordinateFileKind(tf);
  std::vector<int> frame_numbers;
  switch (kind) {
  case CoordinateFileKind::AMBER_CRD:
    {
      // Calculate the number of text lines per frame
      int nline_per_frame = ((atom_count * 3) + 9) / 10;

      // Test whether there is box information
      const int box_line = 1 + nline_per_frame;
      const TextFileReader tfr = tf.data();
      bool box_info_detected = false;
      if (atom_count == 1) {

        // There would be no way to tell a trajectory of a single atom in empty space from a
        // trajectory of a single atom in a box, so use the provided unit cell type.
        switch (unit_cell) {
        case UnitCellType::ORTHORHOMBIC:
        case UnitCellType::TRICLINIC:
          box_info_detected = true;
          break;
        case UnitCellType::NONE:
          break;
        }
      }
      else if (tfr.line_count > box_line) {
        const int llim = tfr.line_limits[box_line];
        const int hlim = tfr.line_limits[box_line + 1];
        int ndots = 0;
        for (int i = llim; i < hlim; i++) {
          ndots += (tfr.text[i] == '.');
        }
        box_info_detected = (ndots == 3);
      }
      nline_per_frame += (box_info_detected);

      // Calculate the number of frames in the file
      const int nframe = (tfr.line_count - 1) / nline_per_frame;
      if (tfr.line_count - (nframe * nline_per_frame) != 0) {
        const std::string boxmsg = (box_info_detected) ? "Box dimensions were detected.  " : "";
        const std::string errmsg = "In file " + tf.getFileName() + ", an atom count of " +
                                   std::to_string(atom_count) + "implies a total of " +
                                   std::to_string(nframe) + " frames.  " + boxmsg +
                                   "However, a total of " +
                                   std::to_string(tfr.line_count - (nframe * nline_per_frame)) +
                                   " lines remain, which cannot be construed as a whole frame.";
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr(errmsg, "getAllFrames");
        case ExceptionResponse::WARN:
          rtWarn(errmsg, "getAllFrames");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }
      frame_numbers.resize(nframe);
      for (int i = 0; i < nframe; i++) {
        frame_numbers[i] = i;
      }
    }
    break;
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_ASCII_RST:
    frame_numbers.resize(1, 0);
    break;
  case CoordinateFileKind::SDF:
    break;
  case CoordinateFileKind::AMBER_NETCDF:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("The format of " + tf.getFileName() + " could not be understood.", "getSelectedFrames");
  }
  return getSelectedFrames(tf, kind, atom_count, unit_cell, frame_numbers);
}

//-------------------------------------------------------------------------------------------------
std::vector<CoordinateFrame> getAllFrames(const std::string &file_name, const int atom_count,
                                          const UnitCellType unit_cell,
                                          const ExceptionResponse policy) {
  const TextFile tf(file_name);
  return getAllFrames(tf, atom_count, unit_cell, policy);
}

//-------------------------------------------------------------------------------------------------
void CoordinateFrame::allocate() {
  const int padded_atom_count  = (atom_count > 0) ? roundUp(atom_count, warp_size_int) : 1;
  const int padded_matrix_size = roundUp(9, warp_size_int);
  storage.resize((3 * padded_atom_count) + (3 * padded_matrix_size));
  x_coordinates.setPointer(&storage,                           0, atom_count);
  y_coordinates.setPointer(&storage,           padded_atom_count, atom_count);
  z_coordinates.setPointer(&storage,       2 * padded_atom_count, atom_count);
  box_space_transform.setPointer(&storage, 3 * padded_atom_count,                             9);
  inverse_transform.setPointer(&storage,  (3 * padded_atom_count) +      padded_matrix_size,  9);
  box_dimensions.setPointer(&storage,     (3 * padded_atom_count) + (2 * padded_matrix_size), 6);
}

} // namespace trajectory
} // namespace stormm
