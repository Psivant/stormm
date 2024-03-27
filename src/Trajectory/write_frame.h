// -*-c++-*-
#ifndef STORMM_WRITE_FRAME_H
#define STORMM_WRITE_FRAME_H

#include <fstream>
#include <string>
#include <vector>
#include "copyright.h"
#include "FileManagement/file_util.h"
#include "Parsing/textfile.h"
#include "Topology/atomgraph_enumerators.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

using diskutil::PrintSituation;
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
/// \param foutp         The just-opened, new trajectory file to begin writing
/// \param output_kind   The trajectory file kind, i.e. AMBER_INPCRD
/// \param atom_count    Number of atoms in the system (if applicable)
/// \param current_time  Current time in the simulation (if applicable)
void initializeTrajectory(std::ofstream *foutp, const CoordinateFileKind output_kind,
                          int atom_count = 0, double current_time = 0.0);

/// \brief Write a coordinate trajectory file or input coordinates file.  Each of the overloaded
///        versions of the function feeds into the base case, which uses pointers to the data
///        at hand, with a verified dimension.
///
/// Overloaded:
///   - Takes double pointers and the number of atoms
///   - Takes std::vector<double> objects
///   - Takes a Coordinates object based on Hybrid objects
///
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
/// \param time_point      The time point to write at the top of a restart file
/// \param 
/// \{
void writeFrame(std::ofstream *foutp, const std::string &filename, CoordinateFileKind kind,
                int natom, const double* x_crd, const double* y_crd, const double* z_crd,
                const double* x_vel, const double* y_vel, const double* z_vel,
                UnitCellType unit_cell, const double* box_dimensions);

void writeFrame(std::ofstream *foutp, const std::string &filename, CoordinateFileKind kind,
                const std::vector<double> &x_crd, const std::vector<double> &y_crd,
                const std::vector<double> &z_crd, const std::vector<double> &x_vel,
                const std::vector<double> &y_vel, const std::vector<double> &z_vel,
                UnitCellType unit_cell, const std::vector<double> &box_dimensions);

void writeFrame(std::ofstream *foutp, const std::string &filename, const TextFile &tf);
/// \}
  
} // namespace trajectory
} // namespace stormm

#include "write_frame.h"

#endif
