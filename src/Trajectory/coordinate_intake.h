// -*-c++-*-
#ifndef STORMM_COORDINATE_INTAKE_H
#define STORMM_COORDINATE_INTAKE_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "FileManagement/file_listing.h"
#include "coordinateframe.h"
#include "phasespace.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

/// \brief Read and return one or more PhaseSpace objects based on a list of file names.  If any of
///        the files are not found, indicate the problem.
///
/// Overloaded:
///   - Try to read in a single file and return the resulting PhaseSpace
///   - Try to read in multiple files and return a list of the resulting PhaseSpace objects
///
/// \param file_name    Name of the file to read as the coordinates
/// \param file_names   List of names to read as coordinates
/// \param priority     Indicate the action to take if files are not found
/// \param files_found  Indicator of whether the files were found (modified and returned if not
///                     the null pointer)
/// \{
PhaseSpace loadPhaseSpace(const std::string &file_name, bool *files_found = nullptr,
                          ExceptionResponse priority = ExceptionResponse::WARN,
                          CoordinateFileKind crd_format = CoordinateFileKind::UNKNOWN);

std::vector<PhaseSpace>
loadPhaseSpace(const std::vector<std::string> &file_names, bool *files_found = nullptr,
               ExceptionResponse priority = ExceptionResponse::WARN,
               CoordinateFileKind crd_format = CoordinateFileKind::UNKNOWN);
/// \}
  
} // namespace trajectory
} // namespace stormm

#endif
