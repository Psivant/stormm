// -*-c++-*-
#ifndef STORMM_FILE_LISTING_H
#define STORMM_FILE_LISTING_H

#include <vector>
#include <string>
#include "copyright.h"

namespace stormm {
namespace diskutil {

/// \brief Enumerate the possible types of directory search (this is equivalent to whether there is
///        a -r modifier after ls).
enum class SearchStyle {
  RECURSIVE, NON_RECURSIVE
};

/// \brief Enumerate whether a path can be determined to be a file, directory, or, barring any
///        other explanation, a regular expression.
enum class DrivePathType {
  FILE, DIRECTORY, REGEXP
};

/// \brief An enumerator to make note of the operating systems
enum class OperatingSystem {
  LINUX, UNIX, WINDOWS, MAC_OS
};

/// Pre-processor directives determine the separators to use in path and file concatenation
/// \{
#if defined(__linux__)
constexpr OperatingSystem detected_os = OperatingSystem::LINUX;
#elif defined(unix) || defined(__unix) || defined(__unix__)
constexpr OperatingSystem detected_os = OperatingSystem::UNIX;
#elif defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
constexpr OperatingSystem detected_os = OperatingSystem::WINDOWS;
#elif defined(__APPLE__) || defined(__MACH__)
constexpr OperatingSystem detected_os = OperatingSystem::MAC_OS;
#endif
/// \}

/// \brief Produce the correct file path separator.
char osSeparator();

/// \brief Test whether a path is a file, directory, or must instead be a regular expression.
///
/// \param path  The path of interest
DrivePathType getDrivePathType(const std::string &path);

/// \brief Given a path that has been established to be a directory, list all files within it.
///        Recursively descend into subdirectories, if requested.
///
/// \param dir_path  The path to search
/// \param r_option  Option to use recursion or not
std::vector<std::string> listDirectory(const std::string &path,
                                       SearchStyle r_option = SearchStyle::NON_RECURSIVE,
                                       DrivePathType entity_kind = DrivePathType::FILE);

/// \brief Detect whether a path appears to be absolute, based on the existence of a leading slash
///        on Unix, Linux, and Mac OS systems or one of the drive letter codes on Windows.
///
/// \param path  The path of the file or directory
bool pathIsAbsolute(const std::string &path);

/// \brief Make a path into an absolute path, modified as appropriate for the system architecture.
///
/// \param path  The name of the file or directory
std::string makePathAbsolute(const std::string &path);

/// \brief Get the normative path by removing any trailing slashes.
///
/// \param path  The path to normalize
std::string getNormPath(const std::string &path);

/// \brief Get the base name of a path by taking everything after the last OS-specific separator.
///
/// \param path  The path from which to extract a base name
std::string getBaseName(const std::string &path);

/// \brief Substitute the extension of a file by chopping off anything past a final '.' character,
///        then adding the supplied extension.  If no such '.' character exists, the extension will
///        be added to the name provided.
///
/// \param path     Path to the file of interest
/// \param new_ext  The new extension to apply to the file
std::string substituteNameExtension(const std::string &path, const std::string &new_ext);

/// \brief Extract the various directory levels in a given path
///
/// \param path  The path to analyze
std::vector<std::string> separatePath(const std::string &path);

/// \brief Split a path into parts before and after the final dot (".")
///
/// \param path    The original path
/// \param before  The part of the path before the final dot, might be an empty string (returned)
/// \param after   The part of the path after the final dot, might be an empty string (returned)
void splitPath(const std::string &path, std::string *before, std::string *after);
  
/// \brief Given a path that has been established to be a regular expression, list all files it
///        could describe, recursively descend into subdirectories if requested.  The recursion
///        only kicks in once the regular expression has been interpreted into a path that is,
///        itself, a directory.  Otherwise, "A/*/B/[a-z]*" will find all subdirectories of A/ such
///        that they contain their own subdirectory B/ and within each such B/ all files beginning
///        with a lowercase letter or, if recursion is activated, descend recursively into any and
///        all subdirectories of B/ whose name begins with a lowercase letter listing all files.
///        The search string "A/*/B/[a-z]*" would not find a file A/foo.txt, but it would find
///        A/bar/B/foo.txt, and if recursion were activated it would find A/bar/B/C/D/E/foo.txt.
///
/// \param regexp_path  The regular expression to evaluate
/// \param r_option     Option to use recursion or not
std::vector<std::string> listFilesInPath(const std::string &regexp_path, SearchStyle r_option);

/// \brief Try to find the STORMM home or source directories based on a hypothetical path.  A known
///        file within the desired STORMM directory is provided to help identify the correct root.
///
/// \param path       The path that may contain clues as to the STORMM home directory
/// \param extension  Name of the file that must be found in some partial path when that partial
///                   path is the correct root
std::string findStormmPath(const std::string &path, const std::string &extension);

/// \brief Form the variable name for common path prefixes.  This encapsulates a convention and
///        applies to output files meant for human comprehension, not actual code or functional
///        scripts produced by STORMM.
///
/// \param index  Index of the path variable to encode
std::string commonPathVariableName(int index);
  
/// \brief Find common paths among a list of file names.  Return the common paths in a separate
///        vector of strings and modify the input to reflect variables set to each common path,
///        which can shorten the overall reporting.  Common paths will be restricted to similar
///        directory tree structures and must begin at the base of each path they modify.  The
///        input will be modified as "${PATH_1}", "${PATH_2}", etc. to substitute the first,
///        second, and additional common, shortcut paths.
///
/// \param input_paths  Input paths to parse for common directory trees.  Modified and returned.
/// \param n_shortcut   The maximum number of shortcut paths to find.  Each shortcut must be
///                     common to at least n_instance paths.
/// \param n_instance   The number of paths that must be feasible to shorten with each shortcut
///                     declaration in order to warrant making a new shortcut
std::vector<std::string> extractCommonPaths(std::vector<std::string> *input_paths,
                                            int n_shortcut = 2, int n_instance = 2);

/// \brief State the list of common path names as a concatenated string.
///
/// \param shortcuts    List of path shortcuts
/// \param indentation  The number of spaces to indent the list of path tokens
std::string listCommonPaths(const std::vector<std::string> &shortcuts, int indentation = 2);
  
} // namespace parse
} // namespace stormm

#endif
