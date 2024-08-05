// -*-c++-*-
#ifndef ERROR_FORMAT_H
#define ERROR_FORMAT_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace errors {

/// \brief Enumerate the kinds of errors and tabular information that might be printed to the
///        terminal or a file at runtime.
enum class RTMessageKind {
  ERROR,    ///< Method used for runtime errors, warnings, and alerts
  TABULAR   ///< Method used for printing tabular data
};

/// \brief Parse output for the terminal screen based on a perceived side and desired indentation.
///
/// \param message          The message to print
/// \param class_caller     The class (or free function) calling this report
/// \param method_caller    A specific class method calling this report
/// \param implicit_indent  Indentation implied by some message printed by another function (i.e.
///                         runtime_error())
/// \param first_indent     First line indentation (after implicit indentation)
/// \param subsq_indent     The degree of indentation for all lines after the first
/// \param width_in         The width to print (defaults to terminal width if <= 0)  
std::string terminalFormat(const std::string &message, const char* class_caller = nullptr,
                           const char* method_caller = nullptr, int implicit_indent = 0,
                           int first_indent = 0, int subsq_indent = 0, int width_in = 0,
			   RTMessageKind style = RTMessageKind::ERROR);

/// \brief The name of this function (abbreviated "runtime error") is deliberately short to permit
///        ease of writing the strings that constitute its arguments within the character limit
///        of each LOC.  This will throw a Standard Template Library runtime_error with a
///        formatted error message designed for a 96-character terminal window.
///
/// \param message        The message to print
/// \param class_caller   The class (or free function) calling this error report
/// \param method_caller  A specific class method calling this error report
void rtErr(const std::string &message, const char* class_caller = nullptr,
           const char* method_caller = nullptr);

/// \brief Warnings are something that the developer needs to hack into C++, so here it is, a
///        runtime warning.  Text is printed to stdout (by printf) but no exception is thrown.
///        The program will continue to run.
///
/// \param message        The message to print
/// \param class_caller   The class (or free function) calling this error report
/// \param method_caller  A specific class method calling this error report
void rtWarn(const std::string &message, const char* class_caller = nullptr,
            const char* method_caller = nullptr);

/// \brief A second type of runtime warning, very similar to rtWarn above.
///
/// \param message        The message to print
/// \param class_caller   The class (or free function) calling this error report
/// \param method_caller  A specific class method calling this error report
void rtAlert(const std::string &message, const char* class_caller = nullptr,
             const char* method_caller = nullptr);

/// \brief Extend a string as if it were a verbal list of items, with the Oxford comma convention.
///
/// \param current_item  The index number of the current item in the list
/// \param item_count    The total number of items that the list is expected to contain
std::string listSeparator(int current_item, int item_count);

} // namespace errors
} // namespace stormm

namespace stormm {
  using errors::rtErr;
  using errors::rtWarn;
  using errors::rtAlert;
} // namespace stormm

#endif
