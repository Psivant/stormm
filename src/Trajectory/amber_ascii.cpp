#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "Math/matrix_ops.h"
#include "Math/vector_ops.h"
#include "Parsing/ascii_numbers.h"
#include "amber_ascii.h"

namespace stormm {
namespace trajectory {

using stmath::computeBoxTransform;
using stmath::extractBoxDimensions;
using stmath::maxValue;
using parse::NumberFormat;
using parse::readNumberSeries;
using parse::separateText;
using parse::TextFileReader;
using parse::verifyNumberFormat;

//-------------------------------------------------------------------------------------------------
int getAmberRestartAtomCount(const TextFile &tf) {
  const int line_idx = tf.getLineLimits(1);
  const int n_char = tf.getLineLimits(2) - line_idx;
  std::vector<std::string> headline_numbers = separateText(tf.getTextPointer(line_idx), n_char);
  if (headline_numbers.size() == 0LLU) {
    rtErr("No atom count found in file " + tf.getFileName() + ".", "getAmberRestartAtomCount");
  }
  else if (verifyNumberFormat(headline_numbers[0].c_str(), NumberFormat::INTEGER) == false) {
    rtErr("Corrupted atom count in file " + tf.getFileName() + ": " + headline_numbers[0] + ".",
          "getAmberRestartAtomCount");
  }
  return stol(headline_numbers[0]);
}

//-------------------------------------------------------------------------------------------------
void getAmberRestartAtomCountAndTime(const TextFile &tf, int *atom_count, double *start_time,
                                     const ExceptionResponse policy) {
  const int line_idx = tf.getLineLimits(1);
  const int n_char = tf.getLineLimits(2) - line_idx;
  std::vector<std::string> headline_numbers = separateText(tf.getTextPointer(line_idx), n_char);
  if (verifyNumberFormat(headline_numbers[0].c_str(), NumberFormat::INTEGER)) {
    *atom_count = stol(headline_numbers[0]);
  }
  else {
    rtErr("Corrupted atom count in file " + tf.getFileName() + ": " + headline_numbers[0] + ".",
          "getAmberRestartAtomCount");
  }
  if (headline_numbers.size() < 2LLU) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("No time value was found in file " + tf.getFileName() + ".",
            "getAmberRestartAtomCountAndTime");
    case ExceptionResponse::WARN:
      rtWarn("No time value was found in file " + tf.getFileName() + ".  The initial time will be "
             "set to 0.0ps.", "getAmberRestartAtomCountAndTime");
      *start_time = 0.0;
      break;
    case ExceptionResponse::SILENT:
      *start_time = 0.0;
      break;
    }
  }
  else if (verifyNumberFormat(headline_numbers[1].c_str(), NumberFormat::SCIENTIFIC) ||
           verifyNumberFormat(headline_numbers[1].c_str(), NumberFormat::SCIENTIFIC)) {
    *start_time = stod(headline_numbers[1]);
  }
  else {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Corrupted start time in file " + tf.getFileName() + ": " + headline_numbers[1] + ".",
            "getAmberRestartAtomCount");
    case ExceptionResponse::WARN:
      rtErr("Corrupted start time in file " + tf.getFileName() + ": " + headline_numbers[1] + ".  "
            "The initial time will be set to 0.0ps.", "getAmberRestartAtomCount");
      *start_time = 0.0;
      break;
    case ExceptionResponse::SILENT:
      *start_time = 0.0;
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
int checkXYZDimensions(const int natom, const int natomy, const int natomz, const char* caller) {
  if (natom == 0 || natom != natomy || natom != natomz) {
    rtErr("A definite, nonzero number of atoms must be allocated for reading in the X, Y, and Z "
          "arrays.  Current allocations: " + std::to_string(natom) + ", " +
          std::to_string(natomy) + ", " + std::to_string(natomz) + ".", caller);
  }
  return natom;
}

//-------------------------------------------------------------------------------------------------
void splitInterlacedCoordinates(const std::vector<PolyNumeric> &allcrd,
                                double* x_ptr, double* y_ptr, double* z_ptr, const int natom) {
  for (int i = 0; i < natom; i++) {
    x_ptr[i] = allcrd[(3 * i)    ].d;
    y_ptr[i] = allcrd[(3 * i) + 1].d;
    z_ptr[i] = allcrd[(3 * i) + 2].d;
  }
}

//-------------------------------------------------------------------------------------------------
void getAmberInputCoordinates(const TextFile &tf, double* x_coordinates, double* y_coordinates,
                              double* z_coordinates, const int natom, double* box_space_transform,
                              double* inverse_transform, double* box_dimensions) {
  std::vector<PolyNumeric> allcrd = readNumberSeries(tf, 2, 3 * natom, 6, 12, 7,
                                                     NumberFormat::STANDARD_REAL,
                                                     "getAmberInputCoordinates", "Read data in "
                                                     "%12.7f format from " + tf.getFileName() +
                                                     ".");
  splitInterlacedCoordinates(allcrd, x_coordinates, y_coordinates, z_coordinates, natom);
  const int lines_per_set = ((natom * 3) + 5) / 6;
  const bool has_velocities = (tf.getLineCount() >= 2 + (2 * lines_per_set));
  const int box_line = 2 + ((1 + has_velocities) * lines_per_set);
  const TextFileReader tfr = tf.data();
  if (tfr.line_count > box_line) {

    // Test this line carefully.  It may contain nothing.
    const int nchar = tfr.line_limits[box_line + 1] - tfr.line_limits[box_line];
    if (nchar < 6 * 12 ||
        separateText(&tfr.text[tfr.line_limits[box_line]], nchar).size() < 6LLU) {
      for (int i = 0; i < 6; i++) {
        box_dimensions[i] = 0.0;
      }
      for (int i = 0; i < 9; i++) {
        box_space_transform[i] = static_cast<double>((i & 0x3) == 0);
        inverse_transform[i]   = static_cast<double>((i & 0x3) == 0);
      }
    }
    else {
      std::vector<double> box_dims(6);
      std::vector<PolyNumeric> boxinfo = readNumberSeries(tf, box_line, 6, 6, 12, 7,
                                                          NumberFormat::STANDARD_REAL,
                                                          "getAmberInputCoordinates", "Read box "
                                                          "information in %12.7f format from " +
                                                          tf.getFileName() + ".");
      for (int i = 0; i < 3; i++) {
        box_dims[i    ] = boxinfo[i    ].d;
        box_dims[i + 3] = boxinfo[i + 3].d * symbols::pi / 180.0;
      }
      for (int i = 0; i < 6; i++) {
        box_dimensions[i] = box_dims[i];
      }
      computeBoxTransform(box_dims.data(), box_space_transform, inverse_transform);
    }
  }
  else {

    // Fill the transformation matrices with the identity matrix.  i & 0x3 makes use of the bitwise
    // AND operator to compute whether i is a multiple of 4 (AND an integer with 2^(N-1) to get a
    // no or yes on whether the integer is a multiple of 2^N, in this case 4).  Fill in matrix
    // slots 0, 4, and 8 (that is, (0,0), (1,1), and (2,2)) with 1.0. 
    for (int i = 0; i < 9; i++) {
      box_space_transform[i] = static_cast<double>((i & 0x3) == 0);
      inverse_transform[i]   = static_cast<double>((i & 0x3) == 0);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void getAmberInputCoordinates(const TextFile &tf, Hybrid<double> *x_coordinates,
                              Hybrid<double> *y_coordinates,  Hybrid<double> *z_coordinates,
                              Hybrid<double> *box_space_transform,
                              Hybrid<double> *inverse_transform, Hybrid<double> *box_dimensions) {
  const int natom = checkXYZDimensions(x_coordinates->size(), y_coordinates->size(),
                                       z_coordinates->size(), "getAmberInputCoordinates");
  getAmberInputCoordinates(tf, x_coordinates->data(), y_coordinates->data(), z_coordinates->data(),
                           natom, box_space_transform->data(), inverse_transform->data(),
                           box_dimensions->data());
}

//-------------------------------------------------------------------------------------------------
void getAmberRestartVelocities(const TextFile &tf, Hybrid<double> *x_velocities,
                               Hybrid<double> *y_velocities,  Hybrid<double> *z_velocities) {
  const int natom = checkXYZDimensions(x_velocities->size(), y_velocities->size(),
                                       z_velocities->size(), "getAmberRestartVelocities");
  std::vector<PolyNumeric> allvel = readNumberSeries(tf, 2 + (((3 * natom) + 5) / 6), 3 * natom,
                                                     6, 12, 7, NumberFormat::STANDARD_REAL,
                                                     "getAmberRestartVelocities", "Read data in "
                                                     "%12.7f format from " + tf.getFileName() +
                                                     ".");
  splitInterlacedCoordinates(allvel, x_velocities->data(), y_velocities->data(),
                             z_velocities->data(), natom);
}

//-------------------------------------------------------------------------------------------------
int countAmberCrdFormatFrames(const TextFile &tf, const int natom, const UnitCellType unit_cell) {
  const int lines_per_frame = (((3 * natom) + 9) / 10) + (unit_cell != UnitCellType::NONE);
  if (tf.getLineCount() < lines_per_frame + 1) {
    rtErr("File " + tf.getFileName() + " has only " + std::to_string(tf.getLineCount()) +
          " lines, not enough to accommodate " + std::to_string(natom) + " atoms.",
          "readAmberCrdFormat");
  }
  return (tf.getLineCount() - 1) / lines_per_frame;
}
  
//-------------------------------------------------------------------------------------------------
void readAmberCrdFormat(const TextFile &tf, double* x_coordinates, double* y_coordinates,
                        double* z_coordinates, const int natom, const UnitCellType unit_cell,
                        double* box_space_transform, double* inverse_transform,
                        double* box_dimensions, const int frame_number)
{
  // Check for a frame dimension line
  const TextFileReader tfr = tf.data();
  if (countAmberCrdFormatFrames(tf, natom, unit_cell) <= frame_number) {
    rtErr("File " + tf.getFileName() + " has only " + std::to_string(tfr.line_count) + " lines, "
          "not enough to accommodate " + std::to_string(natom) + " atoms and " +
          std::to_string(frame_number + 1) + " frames.", "readAmberCrdFormat");    
  }
  bool has_box = false;
  const int lines_for_coords = ((3 * natom) + 9) / 10;
  if (tfr.line_count >= lines_for_coords + 2) {
    const int llim = tfr.line_limits[lines_for_coords + 1];
    const int hlim = tfr.line_limits[lines_for_coords + 2];
    std::vector<std::string> frame_test = separateText(&tfr.text[llim], hlim - llim);
    has_box = (natom != 1 && frame_test.size() == 3LLU);
    if (has_box && unit_cell == UnitCellType::NONE) {
      rtWarn("Trajectory " + tf.getFileName() + " was requested to be read as a non-periodic "
             "series of coordinates, but appears to contain box information.",
             "readAmberCrdFormat");
    }
  }
  const int ihas_box = has_box;
  const int start_line = 1 + frame_number * (lines_for_coords + ihas_box);
  std::vector<PolyNumeric> allcrd = readNumberSeries(tf, start_line, 3 * natom,
                                                     10, 8, 3, NumberFormat::STANDARD_REAL,
                                                     "readAmberCrdFormat", "Read data in "
                                                     "%8.3f format from " + tf.getFileName() +
                                                     ".");
  splitInterlacedCoordinates(allcrd, x_coordinates, y_coordinates, z_coordinates, natom);
  if (has_box) {
    std::vector<PolyNumeric> boxlen = readNumberSeries(tf, start_line + lines_for_coords, 3,
                                                       3, 8, 3, NumberFormat::STANDARD_REAL,
                                                       "readAmberCrdFormat", "Read data in %8.3f "
                                                       "format from " + tf.getFileName() + ".");

    // The box has angles as well, and these all affect its transformation matrix.  The matrix
    // must already be known in some form, in all likelihood from reading an Amber inpcrd file.
    // Unpack the matrix, replace the box lengths, and reassemble the matrix.
    if (maxValue(inverse_transform, 9) < constants::tiny) {
      rtWarn("A pair of valid transformation matrices is required for reading a frame with box "
             "length information out of an Amber ASCII-format (.crd) trajectory file.  An "
             "orthorhombic unit cell will be assumed.", "readAmberCrdFormat");
      for (int i = 0; i < 9; i++) {
        inverse_transform[i] = static_cast<double>((i & 0x3) == 0);
      }
    }
    double lx, ly, lz, alpha, beta, gamma;
    extractBoxDimensions(&lx, &ly, &lz, &alpha, &beta, &gamma, inverse_transform);
    lx = boxlen[0].d;
    ly = boxlen[1].d;
    lz = boxlen[2].d;
    computeBoxTransform(lx, ly, lz, alpha, beta, gamma, box_space_transform, inverse_transform);
    box_dimensions[0] = lx;
    box_dimensions[1] = ly;
    box_dimensions[2] = lz;
    box_dimensions[3] = alpha;
    box_dimensions[4] = beta;
    box_dimensions[5] = gamma;
  }
}

//-------------------------------------------------------------------------------------------------
void readAmberCrdFormat(const TextFile &tf, Hybrid<double> *x_coordinates,
                        Hybrid<double> *y_coordinates,  Hybrid<double> *z_coordinates,
                        const UnitCellType unit_cell, Hybrid<double> *box_space_transform,
                        Hybrid<double> *inverse_transform, Hybrid<double> *box_dimensions,
                        const int frame_number) {
  const int natom = checkXYZDimensions(x_coordinates->size(), y_coordinates->size(),
                                       z_coordinates->size(), "readAmberCrdFormat");
  readAmberCrdFormat(tf, x_coordinates->data(), y_coordinates->data(), z_coordinates->data(),
                     natom, unit_cell, box_space_transform->data(), inverse_transform->data(),
                     box_dimensions->data(), frame_number);
}

//-------------------------------------------------------------------------------------------------
void readAmberCrdFormat(const TextFile &tf, CoordinateFrameWriter *cfw, const int frame_number) {
  readAmberCrdFormat(tf, cfw->xcrd, cfw->ycrd, cfw->zcrd, cfw->natom, cfw->unit_cell, cfw->umat,
                     cfw->invu, cfw->boxdim, frame_number);
}

} // namespace trajectory
} // namespace stormm
