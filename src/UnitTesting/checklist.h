// -*-c++-*-
#ifndef STORMM_CHECKLIST_H
#define STORMM_CHECKLIST_H

#include <string>
#include <vector>
#include "copyright.h"
#include "unit_test_enumerators.h"

namespace stormm {
namespace testing {

/// \brief Object for storing a series of test results, with labels
class CheckList {
public:

  /// \brief Constructor for the CheckList object prepares a blank slate
  CheckList();

  // Overload the += operator to act similarly to the logResult member function
  CheckList& operator+=(CheckResult rhs);

  /// \brief Log a result in a checklist based on the current section setting
  ///
  /// \param result  The result to log
  void logResult(CheckResult result);

  /// \brief Switch to a new (or previously started) section of the test case CheckList.
  ///
  /// Overloaded:
  ///   - Jump to a section by name
  ///   - Jump to a section by index
  ///
  /// \param section_name   The name of the section to add or jump to
  /// \param section_index  The index of the section to jump to
  /// \{
  void changeSection(const std::string &section_name);
  void changeSection(const int section_index);
  /// \}

  /// \brief Get the index of a test section based on its name
  ///
  /// \param name  The name of the section of interest
  int getSectionIndex(const std::string &section_name) const;

  /// \brief Get the name of a test section based on its index
  ///
  /// \param index  The index of the section of interest
  std::string getSectionName(const int section_index) const;

  /// \brief Get the current section number
  int getCurrentSection() const;

  /// \brief Get the current number of successes in a particular section.  The default value of -1
  ///        returns the count for the current section.
  ///
  /// \param section_index  Index of the section of interest
  int getSuccessCount(int section_index = -1) const;

  /// \brief Get the current number of failures in a particular section.  The default value of -1
  ///        returns the count for the current section.
  ///
  /// \param section_index  Index of the section of interest
  int getFailureCount(int section_index = -1) const;

  /// \brief Get the current number of skipped tests in a particular section.  The default value
  ///        of -1 returns the count for the current section.
  ///
  /// \param section_index  Index of the section of interest
  int getSkipCount(int section_index = -1) const;

  /// \brief Get the current number of ignored test failures in a particular section.  The default
  ///        value of -1 returns the count for the current section.
  ///
  /// \param section_index  Index of the section of interest
  int getIgnoredFailureCount(int section_index = -1) const;

  /// \brief Print a summary of test results from this checklist
  ///
  /// \param verbosity  The level of verboseness at which to report
  void printSummary(TestVerbosity verbosity = TestVerbosity::COMPACT) const;

  /// \brief Get the total number of failures across all sections.
  int getOverallFailureCount() const;

  /// \brief Get the total number of skipped tests across all sections.
  int getOverallSkipCount() const;
  
private:
  int current_section;                ///< Index of the current section, for which successes and
                                      ///<   failures can be tabulated
  std::vector<std::string> sections;  ///< List of named sections (default "General")
  std::vector<int> successes;         ///< Successes recorded for tests in each section
  std::vector<int> skips;             ///< Successes recorded for tests in each section
  std::vector<int> ignored_failures;  ///< Tests that failed but were ignored in each section
                                      ///<   (tests that succeeded will be counted among successes,
                                      ///<   even if they were marked to be ignored)
  std::vector<int> failures;          ///< Failures recorded for tests in each section
};

} // namespace testing
} // namespace stormm

#endif
