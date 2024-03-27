#include <algorithm>
#include <dirent.h>
#include <unistd.h>
#include <sys/stat.h>
#include <regex>
#include <string.h>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/series_ops.h"
#include "Reporting/error_format.h"
#include "file_listing.h"

namespace stormm {
namespace diskutil {

using errors::rtErr;
using stmath::incrementingSeries;

//-------------------------------------------------------------------------------------------------
char osSeparator() {
  switch (detected_os) {
  case (OperatingSystem::LINUX):
  case (OperatingSystem::UNIX):
  case (OperatingSystem::MAC_OS):
    return '/';
  case (OperatingSystem::WINDOWS):
    return '\\';
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
DrivePathType getDrivePathType(const std::string &path) {

  // Try getting the status.  If this fails, the path is probably a regular expression.
  struct stat path_stat;
  if (stat(path.c_str(), &path_stat) == 0) {
    if (S_ISREG(path_stat.st_mode) == 1) {
      return DrivePathType::FILE;
    }
    else if (S_ISDIR(path_stat.st_mode) == 1) {
      return DrivePathType::DIRECTORY;
    }
  }

  // Interpret anything else as a regular expression
  return DrivePathType::REGEXP;
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> listDirectory(const std::string &dir_path, const SearchStyle r_option,
                                       const DrivePathType entity_kind) {
  
  // Attempt to open the directory
  std::vector<std::string> ls_result;
  DIR *dir;
  if ((dir = opendir(dir_path.c_str())) == NULL) {
    return ls_result;
  }
  struct dirent *ent;
  while ((ent = readdir(dir)) != NULL) {
    if (strcmp(ent->d_name, ".") == 0 || strcmp(ent->d_name, "..") == 0) {
      continue;
    }
    
    // Test whether the result at hand is, itself, a directory.
    std::string nested_dir_path = dir_path + osSeparator() + std::string(ent->d_name);
    DrivePathType item_type = getDrivePathType(nested_dir_path);
    if (item_type == DrivePathType::DIRECTORY) {
      switch (r_option) {
      case SearchStyle::RECURSIVE:
        {
          std::vector<std::string> nested_result = listDirectory(nested_dir_path, r_option,
                                                                 entity_kind);
          ls_result.insert(ls_result.end(), nested_result.begin(), nested_result.end());
          break;
        }
      case SearchStyle::NON_RECURSIVE:
        if (entity_kind == DrivePathType::DIRECTORY) {
          ls_result.push_back(nested_dir_path);
        }
        break;
      }
    }
    else if (item_type == DrivePathType::FILE) {
      if (entity_kind == DrivePathType::FILE) {
        std::string new_name(ent->d_name);
        ls_result.push_back(dir_path + osSeparator() + new_name);
      }
    }
  }
  closedir(dir);

  return ls_result;
}

//-------------------------------------------------------------------------------------------------
bool pathIsAbsolute(const std::string &path) {

  // Assume that a zero-length path is a relative path to the current directory
  if (path.size() == 0) {
    return false;
  }
  const char sep_char = osSeparator();
  switch (detected_os) {
  case (OperatingSystem::LINUX):
  case (OperatingSystem::UNIX):
  case (OperatingSystem::MAC_OS):
    return (path[0] == sep_char);
  case (OperatingSystem::WINDOWS):
    return (path.size() > 2 && path[0] >= 'A' && path[0] <= 'Z' && path[1] == ':' &&
            path[2] == sep_char);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string makePathAbsolute(const std::string &path) {

  // Return immediately if there is nothing to be done
  if (pathIsAbsolute(path)) {
    return path;
  }

  // Detect the current directory
  char tmp_c[512];
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
  if (GetCurrentDirectory(512, tmp_c) == nullptr) {
    rtErr("Encountered getcwd() error.", "makePathAbsolute");
  }
#else
  if (getcwd(tmp_c, 512) == nullptr) {
    rtErr("Encountered getcwd() error.", "makePathAbsolute");
  }
#endif
  const std::string program_cwd(tmp_c);
  return program_cwd + osSeparator() + path;
}

//-------------------------------------------------------------------------------------------------
std::string getNormPath(const std::string &path) {
  const char sep_char = osSeparator();
  int i = path.size() - 1;
  while (i >= 0 && path[i] == sep_char) {
    i--;
  }
  return path.substr(0, i + 1);
}

//-------------------------------------------------------------------------------------------------
std::string getBaseName(const std::string &path) {
  const char sep_char = osSeparator();
  const int plength = path.size();
  int last_separator = 0;
  for (int i = 0; i < plength; i++) {
    if (path[i] == sep_char) {
      last_separator = i + 1;
    }
  }
  return path.substr(last_separator, plength - last_separator);
}

//-------------------------------------------------------------------------------------------------
std::string substituteNameExtension(const std::string &path, const std::string &new_ext) {
  bool dot_sought = true;
  int pos = path.size() - 1;
  while (dot_sought && pos) {
    dot_sought = (dot_sought && (path[pos] != '.'));
    pos -= dot_sought;
  }
  if (dot_sought) {
    return path + '.' + new_ext;
  }
  else {
    return path.substr(0, pos) + '.' + new_ext;
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
std::vector<std::string> separatePath(const std::string &path) {

  // Detect a leading slash to see if the regular expression contains an absolute path
  const char sep_char = osSeparator();
  const bool absolute_path = (path[0] == sep_char);

  // Break the regular expression into a series of directory levels
  std::vector<std::string> levels;
  int i = 0;
  const int p_length = path.size();
  while (i < p_length) {

    // Take in characters until the next separator
    std::string new_level;
    while (i < p_length) {
      char t_char = path[i];
      i++;
      if (t_char != sep_char) {
        new_level += t_char;
      }
      else {
        break;
      }
    }
    if (new_level.size() > 0) {
      levels.push_back(new_level);
    }
  }

  return levels;
}

//-------------------------------------------------------------------------------------------------
void splitPath(const std::string &path, std::string *before, std::string *after) {
  const int pathsize = path.size();
  int dotpos = pathsize;
  for (int i = pathsize - 1; i >= 0; i--) {
    if (path[i] == '.') {
      dotpos = i;
      break;
    }
  }
  *after = (dotpos < pathsize) ? path.substr(dotpos + 1, pathsize - dotpos - 1) : "";
  *before = (dotpos > 0) ? path.substr(0, dotpos) : "";
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> listFilesInPath(const std::string &regexp_path, SearchStyle r_option) {
  std::vector<std::string> ls_result;

  // Detect a regular file or directory
  DrivePathType td = getDrivePathType(regexp_path);
  switch (td) {
  case DrivePathType::FILE:
    ls_result.push_back(regexp_path);
    return ls_result;
  case DrivePathType::DIRECTORY:
    {
      std::vector<std::string> nested_result = listDirectory(regexp_path, r_option);
      ls_result.insert(ls_result.end(), nested_result.begin(), nested_result.end());
      return ls_result;
    }
  case DrivePathType::REGEXP:
    break;
  }
    
  // Detect a leading slash to see if the regular expression contains an absolute path
  const char sep_char = osSeparator();
  const bool absolute_path = (regexp_path[0] == sep_char);
  
  // Break the regular expression into a series of directory levels
  const std::vector<std::string> levels = separatePath(regexp_path);

  // For each level, look at the appropriate directory and search for matching subdirectories
  const int n_levels = levels.size();  
  std::string partial_path = (absolute_path) ? "" : ".";
  std::regex complete_expr(regexp_path);
  for (int i = 0; i < n_levels; i++) {
    std::string test_path = partial_path + sep_char + levels[i];
    DrivePathType pd = getDrivePathType(test_path);
    switch (pd) {
    case DrivePathType::FILE:
      if (i == n_levels - 1) {
        if (regex_match(test_path, complete_expr)) {
          ls_result.push_back(test_path);
        }
      }
      break;
    case DrivePathType::DIRECTORY:
      partial_path = test_path;
      break;
    case DrivePathType::REGEXP:
      {
        std::regex level_expr(levels[i]);
        std::vector<std::string> test_paths = listDirectory(partial_path,
                                                            SearchStyle::NON_RECURSIVE,
                                                            DrivePathType::FILE);
        std::vector<std::string> additional_paths = listDirectory(partial_path,
                                                                  SearchStyle::NON_RECURSIVE,
                                                                  DrivePathType::DIRECTORY);
        test_paths.insert(test_paths.end(), additional_paths.begin(), additional_paths.end());
        const int n_paths = test_paths.size();
        for (int j = 0 ; j < n_paths; j++) {
          
          // Test that the deepest level of this path conforms to the regular expression at hand.
          if (regex_match(getBaseName(getNormPath(test_paths[j])), level_expr)) {

            // This path is valid.  If there are more levels to explore, continue adding to the
            // partial path.  Otherwise, if the path is a file, add it to the list, or if it is
            // a directory, list its files, using recursion as requested.
            std::string refined_path = test_paths[j];
            for (int k = i + 1; k < n_levels; k++) {
              refined_path += sep_char + levels[k];
            }
            std::vector<std::string> nested_result = listFilesInPath(refined_path, r_option);
            ls_result.insert(ls_result.end(), nested_result.begin(), nested_result.end());
          }
        }

        // Continue incrementing the partial path
        partial_path = test_path;
      }
      break;
    }
  }
  
  return ls_result;
}

//-------------------------------------------------------------------------------------------------
std::string findStormmPath(const std::string &path, const std::string &extension) {
  const std::string abs_path = pathIsAbsolute(path) ? path : makePathAbsolute(path);
  const std::vector<std::string> pieces = separatePath(abs_path);
  const int n_piece = pieces.size();
  const char sep_char = osSeparator();
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
  std::string partial_path("");
#else
  std::string partial_path(1, sep_char);
#endif
  for (int i = 0; i < n_piece; i++) {
    const std::string guess_filename = partial_path + sep_char + extension;
    if (getDrivePathType(guess_filename) == DrivePathType::FILE) {
      return partial_path;
    }
    if (i > 0) {
      partial_path += sep_char;
    }
    partial_path += pieces[i];
  }

  // If this point is reached, a good guess as to ${STORMM_HOME} could not be made 
  return std::string("");
}

//-------------------------------------------------------------------------------------------------
std::string commonPathVariableName(const int index) {
  return std::string("${PATH_") + std::to_string(index + 1) + std::string("}");
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> extractCommonPaths(std::vector<std::string> *input_paths,
                                            const int n_shortcut, const int n_instance) {
  const int npath = input_paths->size();
  std::vector<std::string> result;
  std::vector<int> matches(npath);
  std::vector<int> last_ossep(npath, 0);
  const char osc = osSeparator();
  std::string* input_path_ptr = input_paths->data();
  for (int i = 0; i < npath; i++) {
    const int ilen = input_path_ptr[i].size();
    for (int j = 0; j < ilen; j++) {
      if (input_path_ptr[i][j] == osc) {
        last_ossep[i] = j;
      }
    }
  }
  std::vector<llint3> best_results(npath);
  for (int i = 0; i < npath; i++) {
    for (int j = 0; j < npath; j++) {
      if (i == j) {
        matches[j] = last_ossep[i] + 1;
        continue;
      }
      int nchar = 0;
      while (nchar <= last_ossep[i] && nchar <= last_ossep[j] &&
             input_path_ptr[i][nchar] == input_path_ptr[j][nchar]) {
        nchar++;
      }
      while (nchar >= 0 && input_path_ptr[i][nchar] != osc) {
        nchar--;
      }
      matches[j] = nchar + 1;
    }

    // With the number of matching characters to all other paths now known, sort the list of
    // matches to determine the most effective shortcut that could be generated based on a
    // substring of the ith path.  Weight shortcuts that reduce the length by greater amounts
    // more heavily and have paths that condense into common roots. 
    std::sort(matches.begin(), matches.end(), [](int a, int b) { return a < b; });
    llint best_score = 0;
    int best_length;
    for (int j = 0; j < npath; j++) {
      int k = j + 1;
      while (k < npath && matches[k] == matches[j]) {
        k++;
      }
      llint test_score = static_cast<llint>(k - j) * static_cast<llint>(matches[j] - 11);

      // Assign negative scores to any shortcut path that mathces fewer than the number of
      // instances needed to qualify.
      if (test_score > 0 && (k - j < n_instance)) {
        test_score *= -1;
      }
      if (test_score >= best_score) {
        best_score = test_score;
        best_length = matches[j];
      }

      // Advance the outer loop counter to reflect similar scores found in the inner "while" loop.
      j = k - 1;
    }
    best_results[i] = { best_score, static_cast<llint>(i), static_cast<llint>(best_length) };
  }
  
  // Sort the vector of best results by their scores.  Each tuple contains the score, the index
  // of one of the paths containing the correct base, and the length of the base.
  std::sort(best_results.begin(), best_results.end(),
            [](llint3 a, llint3 b) { return a.x > b.x; });
  int good_options_limit = 0;
  while (good_options_limit < npath && best_results[good_options_limit].x > 0) {
    good_options_limit++;
  }
  
  // The goal is to find not one but up to n_shortcut paths that can reduce the lengths of the
  // various file paths.  Ignore cases with scores of zero or less--these would not reduce the
  // path lengths sufficently to be worth it.
  llint best_overall_score = 0;
  bool score_improving = true;
  std::vector<bool> shortened(npath, false);
  std::vector<int2> best_shortcuts;
  int nsh_skip = 0;
  while (score_improving) {

    // Clear all markings for shortened paths
    for (int i = 0; i < npath; i++) {
      shortened[i] = false;
    }

    // The best solution may be to forego one of the top-scoring shortcuts because it, in fact,
    // precludes other shortucts from applying to the same paths.  While not an exhuastive search
    // method, eliminating such solutions from the top down will help to find the best set of
    // specific and mututally exclusive path shortcuts.
    int skips_to_process = nsh_skip;
    int sh_counter = 0;
    while (skips_to_process > 0) {
      const int path_basis = best_results[sh_counter].y;
      const int shortcut_length = best_results[sh_counter].z;
      for (int i = 0; i < npath; i++) {
        if (shortened[i] == false &&
            strncmp(input_path_ptr[i].data(), input_path_ptr[path_basis].data(),
                    shortcut_length) == 0) {
          shortened[i] = true;
        }
      }

      // Increment the shortcut counter until an unaffected path is found.
      while (sh_counter < npath && shortened[best_results[sh_counter].y]) {
        sh_counter++;
      }

      // Mark one possible, high-scoring shortcut as skipped.
      skips_to_process--;
    }

    // Reset the tracking for whether each path has been affected by some candidate shortcut.
    for (int i = 0; i < npath; i++) {
      shortened[i] = false;
    }
    
    // Initialize the test score, the number of available path shortcuts to use, and the counter
    // for traversing the ordered list of { score, shortcut ID, length } tuples.  Many of those
    // tuples which have the same score will in fact describe the same shortcut and apply to the
    // same systems.  Therefore, with each new shortcut added to the list, keep stepping through
    // the ordered tuples until the path to which the shortcut ID corresponds has not already been
    // affected by one of the previous shortcuts.
    llint test_score = 0;
    int nsh_avail = n_shortcut;
    std::vector<int2> test_shortcuts;
    test_shortcuts.reserve(n_shortcut);
    while (nsh_avail > 0 && sh_counter < good_options_limit) {
      const int path_basis = best_results[sh_counter].y;
      const int shortcut_length = best_results[sh_counter].z;
      const llint score_incr = shortcut_length - 11;
      bool sh_used = false;
      for (int i = 0; i < npath; i++) {
        if (shortened[i] == false &&
            strncmp(input_path_ptr[i].data(), input_path_ptr[path_basis].data(),
                    shortcut_length) == 0) {
          shortened[i] = true;
          test_score += score_incr;
          sh_used = true;
        }
      }
      if (sh_used) {
        test_shortcuts.push_back({ path_basis, shortcut_length });
        nsh_avail--;
      }

      // Continue to increment the shortcut counter until an unaffected path is found.
      sh_counter++;
      while (sh_counter < npath && shortened[best_results[sh_counter].y]) {
	sh_counter++;
      }
    }
    if (test_score > best_overall_score) {
      best_overall_score = test_score;
      best_shortcuts = test_shortcuts;
      score_improving = true;
    }
    else {
      score_improving = false;
    }
  }
  
  // Take the best shortcuts and apply them to all possible paths.
  const int nfound = best_shortcuts.size();
  for (int i = 0; i < npath; i++) {
    shortened[i] = false;
  }
  for (int i = 0; i < nfound; i++) {
    const std::string substitute = commonPathVariableName(i);

    // Omit the trailing directory delimiter in each shortcut path.  Add it back when
    // reconstructing the path with the appropriate shortcut.
    const int sh_eff_len = best_shortcuts[i].y - 1;
    const std::string sh_path = input_path_ptr[best_shortcuts[i].x].substr(0, sh_eff_len);
    result.push_back(sh_path);
    for (int j = 0; j < npath; j++) {
      if (shortened[j] == false &&
          strncmp(sh_path.data(), input_path_ptr[j].data(), sh_eff_len) == 0) {
        input_path_ptr[j] = substitute + osc + input_path_ptr[j].substr(best_shortcuts[i].y);
        shortened[j] = true;
      }
    }
  }
  
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string listCommonPaths(const std::vector<std::string> &shortcuts, const int indentation) {
  std::string result;
  const std::string spacer(indentation, ' ');
  for (size_t i = 0; i < shortcuts.size(); i++) {
    result += spacer + commonPathVariableName(i) + " = " + shortcuts[i] + "\n";
  }
  return result;
}

} // namespace diskutil
} // namespace stormm
