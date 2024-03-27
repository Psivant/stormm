// -*-c++-*-
#ifndef STORMM_SEARCH_DOX_H
#define STORMM_SEARCH_DOX_H

#include <vector>
#include <string>
#include "copyright.h"

namespace stormm {
namespace docs {

using parse::TextFileReader;

/// \brief Stores data on instances of an object (i.e. a function or struct found in a given file.
class ObjectIdentifier {
public:

  /// \brief Constructor for the ObjectIdentifier struct
  ///
  /// \param filename_in  Name of the file being searched
  ObjectIdentifier(const std::string &filename_in);

  /// \brief Default destructor
  ~ObjectIdentifier() = default;

  /// \brief Return the number of instances for an object, as recorded in some ObjectIdentifier
  int getInstanceCount() const;

  /// \brief Return the name of the file searched to produce an ObjectIdentifier
  const std::string getFileName() const;

  /// \brief Report the number of return types for an object from some ObjectIdentifier
  int getReturnTypeCount() const;

  /// \brief Report the return type found for an object, as recorded in some ObjectIdentifier
  ///
  /// \param nt  The index of the return type to return
  const std::string getReturnType(const int nt = 0) const;

  /// \brief Tally one more instance of an object in an ObjectIdentifier
  ///
  /// \param n_add  The number of new instances to tally
  void addInstance(const int n_add = 1);

  /// \brief Set the return type of an object in an ObjectIdentifier
  void addReturnType(const std::string rtype);

private:
  const std::string file_name;           ///< Name of a file searched for this ObjectIdenifier; new
                                         ///<   instances will be generated for each file searched.
  std::vector<std::string> return_types; ///< Return type, if the object is a function (default
                                         ///<   empty string)
  int n_instances;                       ///< Number of variable declarations, if the object is a
                                         ///<   struct, or number of distinct calls and
                                         ///<   declarations if the object is a function
};

/// \brief Enumerator to describe how much information to deliver about an object.
enum class ObjectReportType {
  FULL,             ///< Print both annotation and usage statistics
  ANNOTATION_ONLY,  ///< Print annotation on members or function arguments
  STATS_ONLY        ///< Print usage statistics on the instances of this object across all files
};

/// \brief Enumerator for pre-processor lines of interest
enum class PreProcessorScopeModifier {
  NONE,  ///< The pre-processor line does not modify the scope, i.e. a #define statement
  IF,    ///< The line is a pre-processor if statement
  ELIF,  ///< The line is a pre-processor elif (else-if) statement
  ELSE,  ///< The line is a pre-processor else statement
  ENDIF, ///< The line is a pre-processor endif statement
};

/// \brief Describe a named C++ scope in terms of its place in a file, the namespaces it occupies,
///        and any annotation that appears to be associated with it.
struct CppScope {
  int start_line;                      ///< Line at which this scope starts
  int start_pos;                       ///< Position on the line at which this scope starts
  int end_line;                        ///< Line at which this scope ends
  int end_pos;                         ///< Position on the line at which this scope ends
  std::string scope_name;              ///< Name of the scope (only structs, classes, and functions
                                       ///<   will be named, and thus recorded)
  std::string file_name;               ///< File name in which this scope is found
  std::vector<std::string> namespaces; ///< Namespaces in which this scope exists (namespaces are
                                       ///<   not recorded as named scopes in their own right)

  /// \brief Constructor for a C++ scope struct
  ///
  /// \param filename_in  Name of the file being searched
  CppScope(const std::string &filename_in);

  /// \brief Default destructor
  ~CppScope() = default;
};

/// \brief Determine whether a line is a pre-processor directive, based on whether the line begins
///        with a hash (#) after any amount of whitespace.
///
/// \param line   The line of interest
/// \param nchar  Number of characters on the line (the line will not terminate in a carriage
///               return if it is part of a TextFile, as these characters are removed for brevity)
bool lineIsPreProcessor(const char* line, int nchar);

/// \brief Test a pre-processor line to determine whether it is part of a scope-modifying
///        if/then/else statement.
///
/// \param line   The line of interest
/// \param nchar  Number of characters on the line
PreProcessorScopeModifier testScopeModifier(const char* line, int nchar);

/// \brief Find all #if / #elif / #else / #endif pre-processor scopes within a given file.
///
/// \param tfr  An abstract of pointers into a TextFile object, the test of the file read into RAM
std::vector<int3> findPreProcessorScopes(const TextFileReader &tfr);

/// \brief Find all C++ scopes within { } braces in a given file.
///
/// \param tfr  Pointers an line limits for the text file now stored in RAM
std::vector<CppScope> findCppScopes(const TextFileReader &tfr);

/// \brief Search a particular file for instances of a given object.
///
/// \param object_name  Name of the object (a function or struct) to search for
/// \param filename     File to search  
ObjectIdentifier searchFileForObject(const std::string &object_name, const std::string &filename);

/// \brief Search for an object in the entire STORMM source code tree.
///
/// \param object_name    Name of the object (a function or struct) to search for
/// \param member_name    Name of a struct's member variable or a function's argument
/// \param report_format  Format of the report to write
void searchObject(const std::string &object_name,  const std::string &member_name = std::string(),
                  ObjectReportType report_format = ObjectReportType::FULL);

} // namespace docs
} // namespace stormm


#endif
