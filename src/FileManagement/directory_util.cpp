#include <algorithm>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "Reporting/error_format.h"
#include "directory_util.h"
#include "file_listing.h"

namespace stormm {
namespace diskutil {

using errors::rtErr;

//-------------------------------------------------------------------------------------------------
int localMkdir(const char* path, const mode_t mode) {
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
  return _mkdir(path, mode);
#else
  return mkdir(path, mode);
#endif
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> stormmMkdir(const std::string &path, const mode_t mode) {

  // Keep a list of directories created on the path to (and including) the one of interest
  std::vector<std::string> dirs_created;

  // Detect a leading slash to see if the regular expression contains an absolute path
  const char sep_char = osSeparator();
  const bool absolute_path = (path[0] == sep_char);
  std::vector<std::string> levels = separatePath(path);
  const int n_levels = levels.size();
  std::string partial_path = (absolute_path) ? "" : ".";
  for (int i = 0; i < n_levels; i++) {
    partial_path += sep_char + levels[i];
    const DrivePathType pth_type = getDrivePathType(partial_path);
    switch(getDrivePathType(partial_path)) {
    case DrivePathType::FILE:
      rtErr("The directory " + path + " could not be created because " + partial_path +
            " is a file.", "stormmMkdir");
    case DrivePathType::DIRECTORY:
      break;
    case DrivePathType::REGEXP:
      if (localMkdir(partial_path.c_str(), mode) != 0) {
        rtErr("The directory " + path + " could not be created because " + partial_path +
              " could not be created.", "stormmMkdir");
      }
      dirs_created.push_back(partial_path);
      break;
    }
  }
  return dirs_created;
}

//-------------------------------------------------------------------------------------------------
int localRmdir(const char* path) {
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
  return RemoveDirectory(path);
#else
  return rmdir(path);
#endif
}

//-------------------------------------------------------------------------------------------------
void stormmRmdir(const std::string &path, const ExceptionResponse protocol) {

  // Detect a leading slash to see if the regular expression contains an absolute path
  const char sep_char = osSeparator();
  DrivePathType pth_type = getDrivePathType(path);
  switch (getDrivePathType(path)) {
  case DrivePathType::DIRECTORY:
    if (localRmdir(path.c_str()) != 0) {
      switch (protocol) {
      case ExceptionResponse::DIE:
        rtErr("The directory " + path + " could not be removed.", "stormmRmdir");
      case ExceptionResponse::WARN:
        rtWarn("The directory " + path + " could not be removed.", "stormmRmdir");
	break;
      case ExceptionResponse::SILENT:
	break;
      }
    }
    break;
  case DrivePathType::FILE:
    switch (protocol) {
    case ExceptionResponse::DIE:
      rtErr("The directory " + path + " is, in fact, a file.", "stormmRmdir");
    case ExceptionResponse::WARN:
      rtWarn("The directory " + path + " is a file and will not be removed.", "stormmRmdir");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    break;
  case DrivePathType::REGEXP:
    switch (protocol) {
    case ExceptionResponse::DIE:
      rtErr("Directory " + path + " cannot be removed because it does not exist.", "stormmRmdir");
    case ExceptionResponse::WARN:
      rtWarn("Directory " + path + " cannot be removed because it does not exist.", "stormmRmdir");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void stormmBatchRmdir(const std::vector<std::string> &paths, const ExceptionResponse protocol) {

  // Organize the list if possible.  The longest, deepest paths should go first, followed by
  // shorter ones.
  const int n_paths = paths.size();
  std::vector<std::string> ordered_paths(n_paths, std::string(""));
  std::vector<int2> branch_lengths;
  for (int i = 0; i < n_paths; i++) {
    int2 tmp;
    tmp.x = separatePath(paths[i]).size();
    tmp.y = i;
    branch_lengths.push_back(tmp);
  }
  std::sort(branch_lengths.begin(), branch_lengths.end(),
            [](int2 a, int2 b) { return a.x > b.x; });
  for (int i = 0; i < n_paths; i++) {
    ordered_paths[branch_lengths[i].y] = paths[i];
  }

  // Now start removing directories
  for (int i = 0; i < n_paths; i++) {
    stormmRmdir(ordered_paths[i], protocol);
  }
}

} // namespace diskutil
} // namespace stormm
