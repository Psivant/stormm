#include "copyright.h"
#include "Reporting/error_format.h"
#include "atomgraph_intake.h"

namespace stormm {
namespace topology {

using diskutil::DrivePathType;
using diskutil::getDrivePathType;

//-------------------------------------------------------------------------------------------------
AtomGraph loadTopology(const std::string &file_name, bool *files_found,
                       const ExceptionResponse priority, const TopologyKind engine_format) {
  if (getDrivePathType(file_name) == DrivePathType::FILE) {
    AtomGraph result(file_name, priority, engine_format);
    if (files_found != nullptr) {
      *files_found = true;
    }
    return result;
  }
  else {
    switch (priority) {
    case ExceptionResponse::DIE:
      rtErr("Topology file " + file_name + " was not found.", "loadTopology");
      break;
    case ExceptionResponse::WARN: 
      rtErr("Topology file " + file_name + " was not found.  A blank topology will be returned in "
            "its place.", "loadTopology");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    if (files_found != nullptr) {
      *files_found = false;
    }
    return AtomGraph();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<AtomGraph> loadTopology(const std::vector<std::string> &file_names, bool *files_found,
                                    const TopologyKind engine_format,
                                    const ExceptionResponse priority) {
  std::vector<AtomGraph> result;
  result.reserve(file_names.size());
  if (files_found != nullptr) {
    *files_found = true;
  }
  for (size_t i = 0LLU; i < file_names.size(); i++) {
    if (getDrivePathType(file_names[i]) == DrivePathType::FILE) {
      result.emplace_back(file_names[i], priority, engine_format);
    }
    else {
      switch (priority) {
      case ExceptionResponse::DIE:
        rtErr("Topology file " + file_names[i] + " was not found.", "loadTopology");
        break;
      case ExceptionResponse::WARN:
        rtErr("Topology file " + file_names[i] + " was not found.  A blank topology will be "
              "returned in its place.", "loadTopology");
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
