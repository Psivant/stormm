// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename T>
void writeFrame(std::ofstream *foutp, const std::string &filename, const CoordinateFileKind kind,
                int natom, const T* x_crd, const T* y_crd, const T* z_crd, const T* x_vel,
                const T* y_vel, const T* z_vel, const UnitCellType unit_cell,
                const double* box_dimensions, const BrokenAsciiCode recovery) {
  
  // Enforce real-valued data types
  if (isFloatingPointScalarType<T>() == false) {
    rtErr("This routine is restricted to work in floating point data types.  " +
          getStormmScalarTypeName<T>() + " is invalid.", "writeFrame");
  }

  // Declare arrays that will be filled out later, as necessary
  std::vector<PolyNumeric> pn_allcrd;
  std::vector<PolyNumeric> pn_allvel;
  
  // Lay out the appropriate coordinates array: PolyNumeric for ASCII text speed printing,
  // double-precision real vector for NetCDF and other formats
  switch (kind) {
  case CoordinateFileKind::AMBER_CRD:
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_ASCII_RST:
    pn_allcrd.resize(3 * natom);
    for (int i = 0; i < natom; i++) {
      pn_allcrd[(3 * i)    ].d = x_crd[i];
      pn_allcrd[(3 * i) + 1].d = y_crd[i];
      pn_allcrd[(3 * i) + 2].d = z_crd[i];
    }
    break;
  case CoordinateFileKind::SDF:
  case CoordinateFileKind::PDB:
    rtErr("A separate overload of writeFrame is required to print PDB entries.", "writeFrame");
  case CoordinateFileKind::AMBER_NETCDF:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("The coordinate file type of " + filename + " could not be understood.", "writeFrame");
  }
  switch (kind) {
  case CoordinateFileKind::AMBER_CRD:
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_NETCDF:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    break;
  case CoordinateFileKind::PDB:

    // The PDB file is treated as a special kind of trajectory, written by its own overload of
    // writeFrame().
    rtErr("Protein Data Bank (PDB) format files require information about atoms that can be found "
          "in an AtomGraph class object.  Use a different overload to access this functionality.",
          "writeFrame");
    break;
  case CoordinateFileKind::SDF:

    // The SD file is treated as a special kind of trajectory, written by a specific overload of
    // writeFrame().  However, this function is still used to print some of the content.
    break;
  case CoordinateFileKind::AMBER_ASCII_RST:
    pn_allvel.resize(3 * natom);
    for (int i = 0; i < natom; i++) {
      pn_allvel[(3 * i)    ].d = x_vel[i];
      pn_allvel[(3 * i) + 1].d = y_vel[i];
      pn_allvel[(3 * i) + 2].d = z_vel[i];
    }
    break;
  case CoordinateFileKind::UNKNOWN:

    // The case of an unknown format is handled in the switch above
    break;
  }
  
  // Print the coordinates
  switch (kind) {
  case CoordinateFileKind::AMBER_CRD:
    printNumberSeries(foutp, pn_allcrd, 10, 8, 3, NumberFormat::STANDARD_REAL, "writeFrame",
                      "Write a frame to an Amber-format .crd trajectory file, " + filename,
                      recovery);
    break;
  case CoordinateFileKind::AMBER_INPCRD:
    printNumberSeries(foutp, pn_allcrd, 6, 12, 7, NumberFormat::STANDARD_REAL, "writeFrame",
                      "Write a frame to an Amber-format input coordinates file, " + filename,
                      recovery);
    break;
  case CoordinateFileKind::AMBER_ASCII_RST:
    printNumberSeries(foutp, pn_allcrd, 6, 12, 7, NumberFormat::STANDARD_REAL, "writeFrame",
                      "Write a frame to an Amber-format input coordinates file, " + filename,
                      recovery);
    printNumberSeries(foutp, pn_allvel, 6, 12, 7, NumberFormat::STANDARD_REAL, "writeFrame",
                      "Write a frame to an Amber-format input coordinates file, " + filename,
                      recovery);
    break;
  case CoordinateFileKind::SDF:
  case CoordinateFileKind::PDB:
    break;
  case CoordinateFileKind::AMBER_NETCDF:
  case CoordinateFileKind::AMBER_NETCDF_RST:
  case CoordinateFileKind::UNKNOWN:
    break;
  }

  // Print the unit cell information  
  switch (unit_cell) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    switch (kind) {
    case CoordinateFileKind::AMBER_CRD:
      {
        std::vector<PolyNumeric> pn_boxlen(3);
        for (int i = 0; i < 3; i++) {
          pn_boxlen[i].d = box_dimensions[i];
        }        
        printNumberSeries(foutp, pn_boxlen, 3, 8, 3, NumberFormat::STANDARD_REAL, "writeFrame",
                          "Write box dimensions to an Amber-format coordinate trajectory, " +
                          filename + ".", recovery);
      }
      break;
    case CoordinateFileKind::AMBER_INPCRD:
    case CoordinateFileKind::AMBER_ASCII_RST:
      {
        std::vector<PolyNumeric> pn_boxdim(6);
        for (int i = 0; i < 3; i++) {
          pn_boxdim[i].d = box_dimensions[i];
        }
        for (int i = 3; i < 6; i++) {
          pn_boxdim[i].d = box_dimensions[i] * 180.0 / symbols::pi;
        }
        printNumberSeries(foutp, pn_boxdim, 6, 12, 7, NumberFormat::STANDARD_REAL, "writeFrame",
                          "Write box dimensions to an Amber-format input coordinates file, " +
                          filename + ".", recovery);
      }
      break;
    case CoordinateFileKind::PDB:
    case CoordinateFileKind::SDF:
    case CoordinateFileKind::AMBER_NETCDF:
    case CoordinateFileKind::AMBER_NETCDF_RST:
      break;
    case CoordinateFileKind::UNKNOWN:
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void writeFrame(std::ofstream *foutp, const std::string &filename, const CoordinateFileKind kind,
                const std::vector<T> &x_crd, const std::vector<T> &y_crd,
                const std::vector<T> &z_crd, const std::vector<T> &x_vel,
                const std::vector<T> &y_vel, const std::vector<T> &z_vel,
                const UnitCellType unit_cell, const std::vector<double> &box_dimensions,
                const BrokenAsciiCode recovery) {

  // Check that all arrays are of the same size
  if (x_crd.size() != y_crd.size() || x_crd.size() != z_crd.size()) {
    rtErr("Coordinates cannot be written for x, y, and z vectors of different lengths (" +
          std::to_string(x_crd.size()) + ", " + std::to_string(y_crd.size()) + ", " + 
          std::to_string(x_crd.size()) + ").", "writeAmberCrd");
  }

  // Check that all required data is present
  switch (unit_cell) {
  case UnitCellType::NONE:
    if (box_dimensions.size() > 0 &&
        (box_dimensions[0] > 1.01 || box_dimensions[1] > 1.01 || box_dimensions[2] > 1.01)) {
      rtErr("Box dimensions of [ " +
            realToString(box_dimensions[0], 8, 4, NumberFormat::STANDARD_REAL) + " x " +
            realToString(box_dimensions[1], 8, 4, NumberFormat::STANDARD_REAL) + " x " +
            realToString(box_dimensions[2], 8, 4, NumberFormat::STANDARD_REAL) + " ] were "
            "supplied for a system with no boundary conditions.", "writeFrame");
    }
    break;
  case UnitCellType::ORTHORHOMBIC:
    if (Approx(box_dimensions[3]).test(0.5 * symbols::pi) == false ||
        Approx(box_dimensions[4]).test(0.5 * symbols::pi) == false ||
        Approx(box_dimensions[5]).test(0.5 * symbols::pi) == false) {
      rtErr("An orthorhombic system must have all box angles set to right angles.  Current "
            "[ alpha, beta, gamma ] = [ " +
            realToString(box_dimensions[3], 8, 4, NumberFormat::STANDARD_REAL) + ", " +
            realToString(box_dimensions[4], 8, 4, NumberFormat::STANDARD_REAL) + ", " +
            realToString(box_dimensions[5], 8, 4, NumberFormat::STANDARD_REAL) + " ].",
            "writeFrame");
    }
    break;
  case UnitCellType::TRICLINIC:

    // Rebuild the transformation matrices to ensure that the box angles are sane
    std::vector<double> umat(9), invu(9);
    computeBoxTransform(box_dimensions, &umat, &invu);
    for (int i = 0; i < 9; i++) {
      if (std::isnan(invu[i]) || std::isinf(invu[i])) {
        rtErr("The triclinic system evaluates to a nonsensical transformation matrix: [ " +
              realToString(box_dimensions[0], 8, 4, NumberFormat::STANDARD_REAL) + " x " +
              realToString(box_dimensions[1], 8, 4, NumberFormat::STANDARD_REAL) + " x " +
              realToString(box_dimensions[2], 8, 4, NumberFormat::STANDARD_REAL) + ", " +
              realToString(box_dimensions[3], 8, 4, NumberFormat::STANDARD_REAL) + ", " +
              realToString(box_dimensions[4], 8, 4, NumberFormat::STANDARD_REAL) + ", " +
              realToString(box_dimensions[5], 8, 4, NumberFormat::STANDARD_REAL) + " ].",
              "writeFrame");
      }
    }
    break;
  }

  // Check that the box dimensions array is of the correct size
  if (box_dimensions.size() != 6) {
    rtErr("Invalid vector for box dimensions (" + std::to_string(box_dimensions.size()) + ").",
          "writeAmberCrd");
  }
  writeFrame(foutp, filename, kind, x_crd.size(), x_crd.data(), y_crd.data(), z_crd.data(),
             x_vel.data(), y_vel.data(), z_vel.data(), unit_cell, box_dimensions.data(),
             recovery);
}

} // namespace trajectory
} // namespace stormm
