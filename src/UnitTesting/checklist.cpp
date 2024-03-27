#include <cstring>
#include "copyright.h"
#include "Math/summation.h"
#include "Parsing/parse.h"
#include "Parsing/tabulation.h"
#include "Reporting/error_format.h"
#include "checklist.h"

namespace stormm {
namespace testing {

using errors::RTMessageKind;
using errors::terminalFormat;
using stmath::sum;
using parse::BorderFormat;
using parse::findStringInVector;
using parse::NumberFormat;
using parse::polyNumericVector;
using parse::printTable;
using parse::realDecimalPlaces;
using parse::realToString;

//-------------------------------------------------------------------------------------------------
CheckList::CheckList() :
    current_section{0},
    sections{"General"},
    successes{0},
    skips{0},
    ignored_failures{0},
    failures{0}
{}

//-------------------------------------------------------------------------------------------------
CheckList& CheckList::operator+=(CheckResult rhs) {
  switch (rhs) {
  case CheckResult::SUCCESS:
    this->successes[this->current_section] += 1;
    break;
  case CheckResult::SKIPPED:
    this->skips[this->current_section] += 1;
    break;
  case CheckResult::IGNORED:
    this->ignored_failures[this->current_section] += 1;
    break;
  case CheckResult::FAILURE:
    this->failures[this->current_section] += 1;
    break;
  }
  return *this;
}

//-------------------------------------------------------------------------------------------------
void CheckList::logResult(const CheckResult result) {
  switch (result) {
  case CheckResult::SUCCESS:
    successes[current_section] += 1;
    break;
  case CheckResult::SKIPPED:
    skips[current_section] += 1;
    break;
  case CheckResult::IGNORED:
    ignored_failures[current_section] += 1;
    break;
  case CheckResult::FAILURE:
    failures[current_section] += 1;
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void CheckList::changeSection(const std::string &section_name) {

  // Try to find the requested section in the list
  const int n_section = sections.size();
  for (int i = 0; i < n_section; i++) {
    if (sections[i] == section_name) {
      current_section = i;
      return;
    }
  }

  // If the requested was not found, search for elipses (...) and use them as wildcards
  const int lname = section_name.size();
  std::vector<std::string> segments;
  std::string current_segment("");
  for (int i = 0; i < lname; i++) {
   
    // Identify an elipsis
    if (i < lname - 2 && section_name[i] == '.' && section_name[i + 1] == '.' &&
        section_name[i + 2] == '.') {
      segments.push_back(current_segment);
      i += 2;
      segments.push_back("...");
      current_segment = "";
    }
    else {
      current_segment += section_name[i];
    }
  }
  if (current_segment.size() > 0) {
    segments.push_back(current_segment);
  }
  const int nseg = segments.size();
  int whi = 0;
  while (whi < n_section) {
    const int lsection = sections[whi].size();
    int index = 0;
    bool found = true;
    bool slide = false;
    for (int i = 0; i < nseg; i++) {
      const int lsegment = segments[i].size();
      if (lsegment == 3 && segments[i] == "...") {
        slide = true;
        continue;
      }
      else {
        if (slide) {
          while (index < lsection - lsegment &&
                 strncmp(segments[i].c_str(), &sections[whi][index], lsegment) != 0) {
            index++;
          }
          if (index >= lsection - lsegment) {
            found = false;
          }
          else {
            index += lsegment;
          }
        }
        else {
          if (index >= lsection - lsegment ||
              strncmp(segments[i].c_str(), &sections[whi][index], lsegment) != 0) {
            found = false;
          }
        }
      }
    }
    if (found) {
      current_section = whi;
      return;
    }
    whi++;
  }

  // If the requested section still was not found, make a new one
  sections.push_back(section_name);
  successes.push_back(0);
  skips.push_back(0);
  ignored_failures.push_back(0);
  failures.push_back(0);
  current_section = n_section;
}

//-------------------------------------------------------------------------------------------------
void CheckList::changeSection(const int section_index) {
  if (section_index < sections.size()) {
    current_section = section_index;
  }
  else {
    rtWarn("Section index " + std::to_string(section_index) + " was requested from a list of "
           "only " + std::to_string(sections.size()) + " sections.  The current section will "
           "remain unchanged.", "CheckList", "changeSection");
  }
}

//-------------------------------------------------------------------------------------------------
int CheckList::getSectionIndex(const std::string &section_name) const {
  const int index = findStringInVector(sections, section_name);
  if (index < 0 || index >= static_cast<int>(sections.size())) {
    rtErr("There is no section " + section_name + ".", "CheckList", "getSectionIndex");
  }
  return index;
}

//-------------------------------------------------------------------------------------------------
std::string CheckList::getSectionName(const int index) const {
  if (index < 0 || index >= static_cast<int>(sections.size())) {
    rtErr("Index " + std::to_string(index) + " dos not exist (there are " +
          std::to_string(sections.size()) + " sections in all).", "CheckList", "getSectionName");
  }
  return sections[index];
}

//-------------------------------------------------------------------------------------------------
int CheckList::getCurrentSection() const {
  return current_section;
}

//-------------------------------------------------------------------------------------------------
int CheckList::getSuccessCount(const int section_index) const {
  if (section_index == -1) {
    return successes[current_section];
  }
  else {
    return successes[section_index];
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int CheckList::getFailureCount(const int section_index) const {
  if (section_index == -1) {
    return failures[current_section];
  }
  else {
    return failures[section_index];
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int CheckList::getSkipCount(const int section_index) const {
  if (section_index == -1) {
    return skips[current_section];
  }
  else {
    return skips[section_index];
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int CheckList::getIgnoredFailureCount(const int section_index) const {
  if (section_index == -1) {
    return ignored_failures[current_section];
  }
  else {
    return ignored_failures[section_index];
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int CheckList::getOverallFailureCount() const {
  const size_t n_sections = sections.size();
  int nfail = 0;
  for (size_t i = 0; i < n_sections; i++) {
    nfail += failures[i];
  }
  return nfail;
}

//-------------------------------------------------------------------------------------------------
int CheckList::getOverallSkipCount() const {
  const size_t n_sections = sections.size();
  int nskip = 0;
  for (size_t i = 0; i < n_sections; i++) {
    nskip += skips[i];
  }
  return nskip;
}

//-------------------------------------------------------------------------------------------------
void CheckList::printSummary(const TestVerbosity verbosity) const {

  // Pre-determine the number of sections containing actual tests
  const int n_section = sections.size(); 
  std::vector<int> p_succ;
  std::vector<int> p_fail;
  std::vector<int> p_ignr;
  std::vector<int> p_skip;
  std::vector<std::string> p_sect;
  for (int i = 0; i < n_section; i++) {
    if (successes[i] + failures[i] + ignored_failures[i] + skips[i] > 0) {
      p_succ.push_back(successes[i]);
      p_fail.push_back(failures[i]);
      p_ignr.push_back(ignored_failures[i]);
      p_skip.push_back(skips[i]);
      p_sect.push_back(" - " + sections[i]);
    }
  }
  const int n_filled_sections = p_succ.size();

  // Display based on requested verbosity
  const int n_succ = sum<int>(successes);
  const int n_fail = sum<int>(failures);
  const int n_skip = sum<int>(skips);
  const int n_ignr = sum<int>(ignored_failures);
  const int n_test = n_succ + n_fail + n_skip + n_ignr;
  switch (verbosity) {
  case TestVerbosity::FULL:
    {
      printf("Summary of test results (%2d sections, %4d tests):\n", n_filled_sections,
             n_test);
      if (n_succ > 0 && n_fail == 0 && n_skip == 0) {
        if (n_ignr == 0) {
          printf("  COMPLETE SUCCESS\n\n");
        }
        else {
          printf("  QUALIFIED SUCCESS (all tests ran, but %4d failures were ignored)\n\n", n_ignr);
        }
      }
      else if (n_succ > 0) {
        printf("  MIXED RESULTS (%4d tests passed, %4d failed, %4d ignored, %4d skipped)\n\n",
               n_succ, n_fail, n_ignr, n_skip);
      }
      else if (n_fail > 0 || n_ignr > 0) {
        printf("  COMPLETE FAILURE (%4d failures can be ignored, %4d tests were skipped)\n\n",
               n_ignr, n_skip);
      }
      else {
        printf("  NO TESTS WERE RUN\n\n");
      }
      int max_section_name = 0;
      for (int i = 0; i < n_filled_sections; i++) {
	max_section_name = std::max(max_section_name, static_cast<int>(p_sect[i].size()));
      }
      max_section_name = std::min(max_section_name, 48);
      for (int i = 0; i < n_filled_sections; i++) {
        p_sect[i] = terminalFormat(p_sect[i], "", "", 0, 0, 3, max_section_name + 1,
                                   RTMessageKind::TABULAR);
      }
      printTable({"Test Section Name", "Pass", "Fail", "Ignored", "Skipped"}, p_sect,
                 {polyNumericVector(p_succ), polyNumericVector(p_fail), polyNumericVector(p_ignr),
                  polyNumericVector(p_skip)}, {NumberFormat::INTEGER, NumberFormat::INTEGER,
                 NumberFormat::INTEGER, NumberFormat::INTEGER}, "Fail", BorderFormat::LIGHT, 0, 1);
      printf("\n");
    }
    break;
  case TestVerbosity::COMPACT:
    if (n_succ > 0 && n_fail == 0 && n_skip == 0) {
      printf("Success.  %4d tests passed", n_succ);
      if (n_filled_sections > 1) {
        printf(" in %2d sections", n_filled_sections);
      }
      if (n_ignr > 0) {
	printf(" (%4d failures ignored).\n", n_ignr);
      }
      else {
	printf(".\n");
      }
    }
    else if (n_succ > 0) {
      printf("Mixed results (%4d tests passed, %4d failed, %4d ignored, %4d skipped)\n",
             n_succ, n_fail, n_ignr, n_skip);
    }
    else if (n_fail > 0 || n_ignr > 0) {
      printf("Total failure (%4d failures can be ignored, %4d tests were skipped)",
             n_ignr, n_skip);
      if (n_filled_sections > 1) {
        printf(" in %2d sections:\n", n_filled_sections);
      }
      else {
        printf(".\n");
      }
    }
    else {
      printf("No tests were run.\n");
    }
    break;
  case TestVerbosity::FAILURE_ONLY:
    if (n_fail > 0) {
      printf("%4d tests failed out of %4d (%4d ignored, %4d skipped)\n", n_fail, n_test, n_ignr,
             n_skip);
    }
    break;
  }
}

} // namespace testing
} // namespace stormm
