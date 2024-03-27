// -*-c++-*-
#ifndef STORMM_DIRECTORIES_H
#define STORMM_DIRECTORIES_H

#include <sstream>
#include <sys/stat.h>
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
#  include <direct.h>
#else
#  include <unistd.h>
#endif
#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"

namespace stormm {
namespace diskutil {

using constants::ExceptionResponse;

// \brief Wrapper for various mkdir() functions on different operating systems
///
/// \param path  The full path of the directory to create
/// \param mode  Mode for creation (default 0755, octal representation giving the user write
///              permissions and all others read and execute permissions)
int localMkdir(const char* path, const mode_t mode);

/// \brief Recursive wrapper for localMkdir(), to form directories many layers deep
///
/// \param path  The full path of the directory to create
/// \param mode  Mode for creation (default 0755, octal representation giving the user write
///              permissions and all others read and execute permissions)
std::vector<std::string> stormmMkdir(const std::string &path, const mode_t mode = 0755);

/// \brief Wrapper for various rmdir() functions on different operating systems
///
/// \param path  The full path of the directory to create
int localRmdir(const char* path);

/// \brief Wrapper for localRmdir(), to check sanity and report problems in the appropriate degree
///
/// \param path      The full path of the directory to create
/// \param protocol  Response that will be elicited by a failure of the wrapped rmdir command
void stormmRmdir(const std::string &path, ExceptionResponse protocol = ExceptionResponse::WARN);

/// \brief Clear a set of directories by recursivelyissuing rmdir commands.  This routine begins by
///        trying to order the results in such a way as, so long as each directory is only filled
///        with other directories that are also in the list, each directory will be empty when the
///        program tries to remove it.
///
/// \param paths     The list of directory paths to (attempt to) remove
/// \param protocol  Indication of what to do in the event of an exception
void stormmBatchRmdir(const std::vector<std::string> &paths,
                      ExceptionResponse protocol = ExceptionResponse::WARN);

}
}

#endif
