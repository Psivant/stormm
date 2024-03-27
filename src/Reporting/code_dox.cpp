#include <vector>
#include <string>
#include "copyright.h"
#include "Parsing/parse.h"
#include "FileManagement/file_listing.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "code_dox.h"

namespace stormm {
namespace docs {

using diskutil::osSeparator;
using diskutil::DrivePathType;
using diskutil::SearchStyle;
using diskutil::listFilesInPath;
using parse::LineSpan;
using parse::TextGuard;
using parse::TextFile;

//-------------------------------------------------------------------------------------------------
ObjectIdentifier::ObjectIdentifier(const std::string &filename_in) :
    file_name{filename_in},
    return_types{"NONE"},
    n_instances{0}
{}

//-------------------------------------------------------------------------------------------------
int ObjectIdentifier::getInstanceCount() const {
  return n_instances;
}

//-------------------------------------------------------------------------------------------------
const std::string ObjectIdentifier::getFileName() const {
  return file_name;
}

//-------------------------------------------------------------------------------------------------
int ObjectIdentifier::getReturnTypeCount() const {
  return static_cast<int>(return_types.size());
}

//-------------------------------------------------------------------------------------------------
const std::string ObjectIdentifier::getReturnType(const int nt) const {
  return return_types[nt];
}

//-------------------------------------------------------------------------------------------------
void ObjectIdentifier::addInstance(const int n_add) {
  n_instances += n_add;
}

//-------------------------------------------------------------------------------------------------
void ObjectIdentifier::addReturnType(const std::string rtype) {
  return_types.push_back(rtype);
}

//-------------------------------------------------------------------------------------------------
CppScope::CppScope(const std::string &filename_in) :
  start_line{-1},
  start_pos{-1},
  end_line{-1},
  end_pos{-1},
  scope_name{""},
  file_name{filename_in},
  namespaces{}
{}

//-------------------------------------------------------------------------------------------------
bool lineIsPreProcessor(const char* line, const int nchar) {
  int i = 0;
  while (i < nchar) {
    if (line[i] == '#') {
      return true;
    }
    else if (line[i] != ' ') {
      return false;
    }
    i++;
  }
  return false;
}

//-------------------------------------------------------------------------------------------------
PreProcessorScopeModifier testScopeModifier(const char* line, const int nchar) {
  int i = 0;
  while (i < nchar) {
    if (i < nchar - 1 && (line[i] == 'i' || line[i] == 'I') &&
        (line[i + 1] == 'f' || line[i + 1] == 'F')) {
      return PreProcessorScopeModifier::IF;
    }
    else if (i < nchar - 3 && (line[i] == 'e' || line[i] == 'E') &&
             (line[i + 1] == 'l' || line[i + 1] == 'L') &&
             (line[i + 2] == 'i' || line[i + 2] == 'I') &&
             (line[i + 3] == 'f' || line[i + 3] == 'F')) {
      return PreProcessorScopeModifier::ELIF;
    }
    else if (i < nchar - 3 && (line[i] == 'e' || line[i] == 'E') &&
             (line[i + 1] == 'l' || line[i + 1] == 'L') &&
             (line[i + 2] == 's' || line[i + 2] == 'S') &&
             (line[i + 3] == 'e' || line[i + 3] == 'E')) {
      return PreProcessorScopeModifier::ELSE;
    }
    else if (i < nchar - 4 && (line[i] == 'e' || line[i] == 'E') &&
             (line[i + 1] == 'n' || line[i + 1] == 'N') &&
             (line[i + 2] == 'd' || line[i + 2] == 'D') &&
             (line[i + 3] == 'i' || line[i + 3] == 'I') &&
             (line[i + 4] == 'f' || line[i + 4] == 'F')) {
      return PreProcessorScopeModifier::ENDIF;
    }
    i++;
  }
  return PreProcessorScopeModifier::NONE;
}

//-------------------------------------------------------------------------------------------------
std::vector<int3> findPreProcessorScopes(const TextFileReader &tfr) {

  // Make a list of lines with pre-processor directives
  std::vector<int> pp_lines;
  for (int i = 0; i < tfr.line_count; i++) {
    const int nchar = tfr.line_limits[i + 1] - tfr.line_limits[i];
    if (lineIsPreProcessor(&tfr.text[tfr.line_limits[i]], nchar)) {
      pp_lines.push_back(i);
    }
  }
  const int n_pp_line = pp_lines.size();
  std::vector<PreProcessorScopeModifier> sc_type(n_pp_line, PreProcessorScopeModifier::NONE);
  for (int i = 0; i < n_pp_line; i++) {
    const int nchar = tfr.line_limits[pp_lines[i] + 1] - tfr.line_limits[pp_lines[i]];
    sc_type[i] = testScopeModifier(&tfr.text[tfr.line_limits[pp_lines[i]]], nchar);
  }

  // Sort the lines with pre-processor directives into a list of scopes
  std::vector<int3> scopes;
  bool do_search = (n_pp_line > 0);
  int search_depth = 1;
  while (do_search) {
    int current_depth = 0;
    ulint n_found = scopes.size();
    int3 tsm;
    for (int i = 0; i < n_pp_line; i++) {
      if (sc_type[i] == PreProcessorScopeModifier::IF) {
        current_depth++;
        if (current_depth == search_depth) {
          tsm.x = search_depth;
          tsm.y = i + 1;
        }
      }
      else if ((sc_type[i] == PreProcessorScopeModifier::ELIF ||
                sc_type[i] == PreProcessorScopeModifier::ELSE) && current_depth == search_depth) {
        tsm.z = i;
        scopes.push_back(tsm);
        n_found++;
        tsm.x = search_depth;
        tsm.y = i + 1;
      }
      else if (sc_type[i] == PreProcessorScopeModifier::ENDIF) {
        if (current_depth == search_depth) {
          tsm.z = i;
          scopes.push_back(tsm);
        }
        current_depth--;
      }
    }
    do_search = (scopes.size() > n_found);
  }

  // Return the result as a vector of int3 objects, x = level, y = first line of scope, and
  // z = last line of scope.  The list will have to be searched exhaustively to determine what
  // scope any particular line of code belongs to.  Some scopes will be 'linked' in that only
  // one of a chain of scopes can be part of the code at once.  This will be detectable by
  // the first line of one scope being only one more than the last line of the previous scope.
  return scopes;
}

//-------------------------------------------------------------------------------------------------
std::vector<CppScope> findCppScopes(const TextFileReader &tfr) {

  // Get text that is in a quote or a comment
  std::vector<TextGuard> quote_guards = { TextGuard("\"", "\"") };
  std::vector<TextGuard> comment_guards = { TextGuard("//"),
                                            TextGuard("/*", "*/", LineSpan::MULTIPLE) };
  std::vector<bool> in_quotes  = markGuardedText(tfr, quote_guards, comment_guards);
  std::vector<bool> in_comment = markGuardedText(tfr, comment_guards, quote_guards);

  // Using text that is not in a quote or comment, find all of the { and } pairs
  std::vector<CppScope> result;

  return result;
}

//-------------------------------------------------------------------------------------------------
ObjectIdentifier searchFileForObject(const std::string &object_name, const std::string &filename) {

  // Read the text file
  const TextFile tf(filename);
  const TextFileReader tfr = tf.data();

  // Get text that is in a quote or a comment
  std::vector<TextGuard> quote_guards = { TextGuard("\"", "\"") };
  std::vector<TextGuard> comment_guards = { TextGuard("//"),
                                            TextGuard("/*", "*/", LineSpan::MULTIPLE) };
  std::vector<bool> in_quotes  = markGuardedText(tfr, quote_guards, comment_guards);
  std::vector<bool> in_comment = markGuardedText(tfr, comment_guards, quote_guards);

  // Search the file
  ObjectIdentifier tobj(filename);
  const int objn_length = object_name.size();
  for (int i = 0; i < tfr.line_count; i++) {
    for (int j = tfr.line_limits[i]; j < tfr.line_limits[i + 1] - objn_length; j++) {
      bool found = true;
      for (int k = 0; k < objn_length; k++) {
        if (tfr.text[j + k] != object_name[k] || in_quotes[j + k] || in_comment[j + k]) {
          found = false;
          break;
        }
      }
      if (found) {
        tobj.addInstance();

        // Does the object define a scope in this instance?

      }
    }
  }

  return tobj;
}

//-------------------------------------------------------------------------------------------------
void searchObject(const std::string &object_name, const std::string &member_name,
                  const ObjectReportType report_format) {

  // Get the STORMM home directory
  const std::string stormm_home = std::getenv("STORMM_HOME");

  // Get a list of directories in STORMM's src/ directory, then make lists of .cpp and .h files
  const std::vector<std::string> src_dirs = listDirectory(stormm_home + osSeparator() + "src",
                                                          SearchStyle::RECURSIVE,
                                                          DrivePathType::DIRECTORY);
  std::vector<std::string> stormm_cpp;
  std::vector<std::string> stormm_hdr;
  const int n_src_dir = src_dirs.size();
  for (int i = 0; i < n_src_dir; i++) {
    std::vector<std::string> tmp_cpp = listFilesInPath(src_dirs[i] + osSeparator() + ".*.cpp",
                                                       SearchStyle::NON_RECURSIVE);
    stormm_cpp.insert(stormm_cpp.end(), tmp_cpp.begin(), tmp_cpp.end());
    std::vector<std::string> tmp_hdr = listFilesInPath(src_dirs[i] + osSeparator() + ".*.h",
                                                       SearchStyle::NON_RECURSIVE);
    stormm_hdr.insert(stormm_hdr.end(), tmp_hdr.begin(), tmp_hdr.end());
  }

  // Go through all files, looking for the object.
  const int n_cpp = stormm_cpp.size();
  const int n_hdr = stormm_hdr.size();
  std::vector<ObjectIdentifier> obj_docs;
  for (int i = 0; i < n_cpp; i++) {
    ObjectIdentifier tmp_id = searchFileForObject(object_name, stormm_cpp[i]);
    obj_docs.push_back(tmp_id);
  }
  for (int i = 0; i < n_hdr; i++) {
    ObjectIdentifier tmp_id = searchFileForObject(object_name, stormm_hdr[i]);
    obj_docs.push_back(tmp_id);
  }
}

} // namespace docs
} // namespace stormm
