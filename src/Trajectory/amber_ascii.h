// -*-c++-*-
#ifndef STORMM_AMBER_ASCII_H
#define STORMM_AMBER_ASCII_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Accelerator/hybrid.h"
#include "Parsing/polynumeric.h"
#include "Parsing/textfile.h"
#include "coordinateframe.h"

namespace stormm {
namespace trajectory {

using constants::ExceptionResponse;
using card::Hybrid;
using parse::PolyNumeric;
using parse::TextFile;

/// \brief Obtain just the atom count from an Amber ASCII-format restart (or inpcrd) file.
///
/// \param tf  The text file to read from (must already be parsed and in RAM, as other routines
///            will go back and read other parts of it--reading the file twice from disk would
///            be wasteful)
int getAmberRestartAtomCount(const TextFile &tf);

/// \brief Obtain the atom count and time from an Amber ASCII-format restart file.
///
/// \param tf          The text file to read from (must already be parsed and in RAM, as other
///                    routines will go back and read other parts of it--reading the file twice
///                    from disk would be wasteful)
/// \param atom_count  The atom count (returned)
/// \param start_time  The initial time (returned)
void getAmberRestartAtomCountAndTime(const TextFile &tf, int *atom_count, double *start_time,
                                     ExceptionResponse policy = ExceptionResponse::WARN);

/// \brief Split the interlaced stream of XYZ coordinates from an Amber file into three arrays of
///        contiguous X, Y, and Z coordinates.
///
/// Overloaded:
///   - Take three Hybrid objects for each coordinate dimension
///   - Take a CoordinateFrameWriter object
///
/// \param allcrd  Stream of numbers read from some ASCII-format trajectory
/// \param x_ptr   Array to hold X coordinates (returned)
/// \param y_ptr   Array to hold Y coordinates (returned)
/// \param z_ptr   Array to hold Z coordinates (returned)
/// \param natom   Number of atoms in the system
void splitInterlacedCoordinates(const std::vector<PolyNumeric> &allcrd, double* x, double* y,
                                double* z, int natom);

/// \brief Read coordinates from an Amber inpcrd or ASCII-format restart file.
///
/// Overloaded:
///   - Accept pointers to Hybrid objects and burrow into their host-side data
///   - Accept pointers to C-style arrays
///
/// \param tf                   The text file to read from
/// \param x_coordinates        Cartesian X coordinates of all particles (returned)
/// \param y_coordinates        Cartesian Y coordinates of all particles (returned)
/// \param z_coordinates        Cartesian Z coordinates of all particles (returned)
/// \param natom                Number of atoms (if C-style arrays are provided)
/// \param box_space_transform  Transformation to take coordinates into fractional (box) space
/// \param inverse_transform    Transformation to take fractional coordinates into real space
/// \param box_dimensions       Box dimensions (taken directly from the file)
/// \{
void getAmberInputCoordinates(const TextFile &tf, double* x_coordinates, double* y_coordinates,
                              double* z_coordinates, int natom, double* box_space_transform,
                              double* inverse_transform, double* box_dimensions);
  
void getAmberInputCoordinates(const TextFile &tf, Hybrid<double> *x_coordinates,
                              Hybrid<double> *y_coordinates,  Hybrid<double> *z_coordinates,
                              Hybrid<double> *box_space_transform,
                              Hybrid<double> *inverse_transform, Hybrid<double> *box_dimensions);
/// \}
  
/// \brief Read velocities from an Amber ASCII-format restart file.
///
/// \param tf            The text file to read from
/// \param x_velocities  Cartesian X velocities of all particles (returned)
/// \param y_velocities  Cartesian Y velocities of all particles (returned)
/// \param z_velocities  Cartesian Z velocities of all particles (returned)
void getAmberRestartVelocities(const TextFile &tf, Hybrid<double> *x_velocities,
                               Hybrid<double> *y_velocities,  Hybrid<double> *z_velocities);

/// \brief Get the number of frames in an Amber .crd format file, based on the expected length of
///        each frame and the number of lines in the file.
///
/// \param tf         Contents of the Amber .crd format trajectory file, read into RAM
/// \param natom      The expected number of atoms (must be known in advance)
/// \param unit_cell  The expected unit cell type (indicates whether there will there be a line
///                   of box lengths in the file)
int countAmberCrdFormatFrames(const TextFile &tf, int natom, UnitCellType unit_cell);

/// \brief Read coordinates from an Amber ASCII-format .crd trajectory.
///
/// Overloaded:
///   - Read into C-style arrays for X, Y, Z, and box transformations
///   - Read into Hybrid objects, i.e. components of CoordinateFrame or PhaseSpace objects
///   - Build coordinates into a CoordinateFrameWriter
///
/// \param tf                   The text file to read from
/// \param x_coordinates        Cartesian X coordinates of all particles (returned)
/// \param y_coordinates        Cartesian Y coordinates of all particles (returned)
/// \param z_coordinates        Cartesian Z coordinates of all particles (returned)
/// \param natom                Number of atoms
/// \param unit_cell            The type of unit cell to expect (this will be checked against the
///                             actual content of the trajectory file)
/// \param box_space_transform  Transformation matrix taking coordinates into fractional space
/// \param inverse_transform    Transformation matrix taking coordinates back to real space
/// \param box_dimensions       The six critical box measurements: three lengths and three angles
/// \param frame_number         The frame number to seek in the file.  Raises an exception if the
///                             end of the file is reached before finding this frame.
/// \param cfw                  Host-side coordinate frame writer object with pointers
/// \{
void readAmberCrdFormat(const TextFile &tf, double *x_coordinates, double *y_coordinates,
                        double *z_coordinates, int natom, UnitCellType unit_cell,
                        double *box_space_transform, double *inverse_transform,
                        double* box_dimensions, int frame_number = 0);

void readAmberCrdFormat(const TextFile &tf, Hybrid<double> *x_coordinates,
                        Hybrid<double> *y_coordinates,  Hybrid<double> *z_coordinates,
                        UnitCellType unit_cell, Hybrid<double> *box_space_transform,
                        Hybrid<double> *inverse_transform, Hybrid<double> *box_dimensions,
                        int frame_number = 0);

void readAmberCrdFormat(const TextFile &tf, CoordinateFrameWriter *cfw, int frame_number = 0);
/// \}

} // namespace trajectory
} // namespace stormm

#endif
