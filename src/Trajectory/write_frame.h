// -*-c++-*-
#ifndef STORMM_WRITE_FRAME_H
#define STORMM_WRITE_FRAME_H

#include <fstream>
#include <string>
#include <vector>
#include "copyright.h"
#include "DataTypes/common_types.h"
#include "FileManagement/file_util.h"
#include "Math/matrix_ops.h"
#include "Parsing/ascii_numbers.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/polynumeric.h"
#include "Parsing/textfile.h"
#include "Reporting/reporting_enumerators.h"
#include "Topology/atomgraph_enumerators.h"
#include "UnitTesting/approx.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

using data_types::isFloatingPointScalarType;
using data_types::getStormmScalarTypeName;
using diskutil::PrintSituation;
using parse::NumberFormat;
using parse::PolyNumeric;
using parse::printNumberSeries;
using parse::realToString;
using review::BrokenAsciiCode;
using stmath::computeBoxTransform;
using testing::Approx;
using topology::UnitCellType;
using parse::TextFile;

/// \brief Modify the expected file opening approach in light of the nature of the output.  Check
///        for impossible expectations.
///
/// \param expectation
/// \param output_kind  The type of output being written
/// \param caller       Name of the calling object (for error reporting purposes)
/// \param method       Name of the calling member function (for error reporting purposes)
PrintSituation adjustTrajectoryOpeningProtocol(const PrintSituation expectation,
                                               const CoordinateFileKind output_kind,
                                               const char* caller = nullptr,
                                               const char* method = nullptr);

/// \brief Print the opening lines or basic details of one of the trajectories.
///
/// Overloaded:
///   - Accept a pointer to the open output stream
///   - Accept a file name and output format, open a file stream to call the first overload,
///     then close the file stream
///
/// \param foutp         The just-opened, new trajectory file to begin writing
/// \param output_kind   The trajectory file kind, i.e. AMBER_INPCRD
/// \param atom_count    Number of atoms in the system (if applicable)
/// \param current_time  Current time in the simulation (if applicable)
/// \{
void initializeTrajectory(std::ofstream *foutp, const CoordinateFileKind output_kind,
                          int atom_count = 0, double current_time = 0.0);

void initializeTrajectory(const std::string &traj_name, const PrintSituation expectation,
                          const CoordinateFileKind output_kind,
                          const int atom_count, const double current_time);
/// \}

/// \brief Write a coordinate trajectory file or input coordinates file.  Each of the overloaded
///        versions of the function feeds into the base case, which uses pointers to the data
///        at hand, with a verified dimension.
///
/// Overloaded:
///   - Takes double pointers and the number of atoms
///   - Takes std::vector<double> objects
///   - Takes a Coordinates object based on Hybrid objects
///
/// \param foutp           Pointer to the open file handle for writing output.  This stream will
///                        be created by something such as openOutputFile() (see file_util.h in
///                        the FileManagement directory) to ensure that the write mode is
///                        appropriate for the use case.
/// \param filename        Name of the file to write
/// \param file_kind       The type of coordinate (or restart) file to write
/// \param expectation     Dictates writing behavior based on the presence or absence of any
///                        existing file with the same name
/// \param title           Title of the coordinate file, to place on the first line
/// \param natom           The number of atoms in the system
/// \param x_crd           Vector / array of Cartesian x coordinates for all atoms
/// \param y_crd           Vector / array of Cartesian y coordinates for all atoms
/// \param z_crd           Vector / array of Cartesian z coordinates for all atoms
/// \param unit_cell       The unit cell type (to avoid needing to infer the lack of periodic
///                        boundary conditions from some special settings of the box dimensions)
/// \param box_dimensions  Six-element vector of box dimensions (can be obtained from a
///                        transformation matrix)
/// \param recovery        Indicate a string to write in place of very large numbers that would
///                        break a fixed-column output format
/// \{
template <typename T>
void writeFrame(std::ofstream *foutp, const std::string &filename, CoordinateFileKind kind,
                int natom, const T* x_crd, const T* y_crd, const T* z_crd, const T* x_vel,
                const T* y_vel, const T* z_vel, UnitCellType unit_cell,
                const double* box_dimensions, BrokenAsciiCode recovery = BrokenAsciiCode::NONE);

template <typename T>
void writeFrame(std::ofstream *foutp, const std::string &filename, CoordinateFileKind kind,
                const std::vector<T> &x_crd, const std::vector<T> &y_crd,
                const std::vector<T> &z_crd, const std::vector<T> &x_vel,
                const std::vector<T> &y_vel, const std::vector<T> &z_vel,
                UnitCellType unit_cell, const std::vector<double> &box_dimensions,
                BrokenAsciiCode recovery = BrokenAsciiCode::NONE);
/// \}
  
} // namespace trajectory
} // namespace stormm

#include "write_frame.tpp"

#endif
