#include "copyright.h"
#include "Reporting/error_format.h"
#include "coordinate_intake.h"

namespace stormm {
namespace trajectory {

using diskutil::DrivePathType;
using diskutil::getDrivePathType;
  
//-------------------------------------------------------------------------------------------------
PhaseSpace loadPhaseSpace(const std::string &file_name, bool *files_found,
                          const ExceptionResponse priority, const CoordinateFileKind crd_format) {
  if (getDrivePathType(file_name) == DrivePathType::FILE) {
    PhaseSpace result(file_name, crd_format);
    if (files_found != nullptr) {
      *files_found = true;
    }
    return result;
  }
  else {
    switch (priority) {
    case ExceptionResponse::DIE:
      rtErr("Coordinate file " + file_name + " was not found.", "loadPhaseSpace");
      break;
    case ExceptionResponse::WARN: 
      rtErr("Coordinate file " + file_name + " was not found.  A blank topology will be returned "
            "in its place.", "loadPhaseSpace");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    if (files_found != nullptr) {
      *files_found = false;
    }
    return PhaseSpace();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<PhaseSpace> loadPhaseSpace(const std::vector<std::string> &file_names,
                                       bool *files_found, const ExceptionResponse priority,
                                       const CoordinateFileKind crd_format) {
  std::vector<PhaseSpace> result;
  result.reserve(file_names.size());
  if (files_found != nullptr) {
    *files_found = true;
  }
  for (size_t i = 0LLU; i < file_names.size(); i++) {
    if (getDrivePathType(file_names[i]) == DrivePathType::FILE) {
      result.emplace_back(file_names[i], crd_format);
    }
    else {
      switch (priority) {
      case ExceptionResponse::DIE:
        rtErr("Coordinate file " + file_names[i] + " was not found.", "loadPhaseSpace");
        break;
      case ExceptionResponse::WARN:
        rtErr("Coordinate file " + file_names[i] + " was not found.  A blank coordiante object "
              "will be returned in its place.", "loadPhaseSpace");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
      result.emplace_back();
      if (files_found != nullptr) {
        *files_found = false;
      }
    }
  }
  return result;
}

} // namespace topology
} // namespace stormm
