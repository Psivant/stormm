#include "copyright.h"
#include "FileManagement/file_util.h"
#include "Parsing/parse.h"
#include "display.h"
#include "section_contents.h"
#include "summary_file.h"

namespace stormm {
namespace review {

using display::horizontalRule;
using parse::alphabetNumber;
using parse::TextFileReader;

//-------------------------------------------------------------------------------------------------
SectionContents::SectionContents(const std::string &title_in, const std::string &file_name_in,
                                 const int width_in, const bool subsection_in,
                                 const int section_id_in, const int subsection_id_in,
                                 const ListEnumeration section_marking_in,
                                 const ListEnumeration subsection_marking_in,
                                 const OutputSyntax style_in, const double table_precision_in) :
    title{title_in},
    component_count{0},
    width{width_in},
    subsection{subsection_in},
    section_id{section_id_in},
    subsection_id{subsection_id_in},
    section_marking{section_marking_in},
    subsection_marking{subsection_marking_in},
    style{style_in},
    narrations{}, lists{}, tables{}, sources{}, scripts{},
    file_name{file_name_in}
{}

//-------------------------------------------------------------------------------------------------
SectionContents::SectionContents(const std::string &title_in, const std::string &file_name_in,
                                 const int width_in, const OutputSyntax style_in,
                                 const double table_precision_in) :
  SectionContents(title_in, file_name_in, width_in, false, 0, 0, ListEnumeration::NUMBERED,
                  ListEnumeration::NUMBERED, style_in, table_precision_in)
{}

//-------------------------------------------------------------------------------------------------
const std::string& SectionContents::getOutputFileName() {
  return file_name;
}

//-------------------------------------------------------------------------------------------------
int SectionContents::getWidth() const {
  return width;
}

//-------------------------------------------------------------------------------------------------
OutputSyntax SectionContents::getSyntax() const {
  return style;
}

//-------------------------------------------------------------------------------------------------
int SectionContents::getComponentCount() const {
  return component_count;
}

//-------------------------------------------------------------------------------------------------
int SectionContents::getComponentCount(const SectionComponent kind) const {
  switch (kind) {
  case SectionComponent::NARRATIVE:
    return narrations.size();
  case SectionComponent::LIST:
    return lists.size();
  case SectionComponent::TABLE:
    return tables.size();
  case SectionComponent::SCRIPT:
    return scripts.size();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
SectionComponent SectionContents::getComponentType(const int index) const {
  validateComponentIndex(index);
  return sources[index];
}

//-------------------------------------------------------------------------------------------------
bool SectionContents::isSubsection() const {
  return subsection;
}
  
//-------------------------------------------------------------------------------------------------
std::string SectionContents::sectionAsString(const OutputSyntax alt_style, const int alt_width,
                                             const int alt_section_id, const int alt_subsection_id,
                                             const int list_indent, const int nested_list_indent,
                                             const ListEnumeration list_numbering,
                                             const ListEnumeration nested_list_numbering) const {
  const char hide = commentSymbol(alt_style);
  const std::string strhide(1, hide);
  const int actual_width = (alt_width >= 0) ? alt_width : width;
  const int actual_section_id = (alt_section_id >= 0) ? alt_section_id : section_id;
  const int actual_subsection_id = (alt_subsection_id >= 0) ? alt_subsection_id : subsection_id;
  int narrative_con = 0;
  int list_con = 0;
  int table_con = 0;
  int script_con = 0;

  // Print the section heading
  std::string result = strhide + " ";
  if (subsection) {
    result += horizontalRule("//", "-", actual_width - 2, '-');
  }
  else {
    result += horizontalRule("[", "]", actual_width - 2, '=');
  }
  std::string title_card = "Section ";
  switch (section_marking) {
  case ListEnumeration::BULLET:
  case ListEnumeration::NONE:
    rtErr("Sections must be numbered with alphanumeric characters or Roman numerals.",
          "SectionContents", "sectionAsString");
  case ListEnumeration::NUMBERED:
    title_card += std::to_string(actual_section_id + 1);
    break;
  case ListEnumeration::ALPHABETIC:
    title_card += alphabetNumber(actual_section_id);
    break;
  case ListEnumeration::ROMAN:
    title_card += toRoman(actual_section_id + 1);
    break;
  }
  if (subsection) {
    switch (subsection_marking) {
    case ListEnumeration::BULLET:
    case ListEnumeration::NONE:
      rtErr("Subsections must be numbered with alphanumeric characters or Roman numerals.",
            "SectionContents", "sectionAsString");
    case ListEnumeration::NUMBERED:
      title_card += "." + std::to_string(actual_subsection_id + 1);
      break;
    case ListEnumeration::ALPHABETIC:
      title_card += "-" + alphabetNumber(actual_subsection_id);
      break;
    case ListEnumeration::ROMAN:
      title_card += "-" + toRoman(actual_subsection_id + 1);
      break;
    }
  }
  title_card += " :";
  title_card = indentText(title, title_card, 0, actual_width - 2, false);
  result += protectText(title_card, hide, actual_width) + strhide + "\n";
  
  // Print the section contents
  for (int i = 0; i < component_count; i++) {
    switch (sources[i]) {
    case SectionComponent::NARRATIVE:
      result += protectText(narrations[narrative_con], hide, actual_width);
      narrative_con++;
      break;
    case SectionComponent::LIST:
      result += protectText(lists[list_con].printList(actual_width - 2, list_indent,
                                                      nested_list_indent, list_numbering,
                                                      nested_list_numbering),
                            hide, actual_width, TextEnds::NEWLINE);
      list_con++;
      break;
    case SectionComponent::TABLE:
      result += tables[table_con].printTable(alt_style, actual_width);
      table_con++;
      break;
    case SectionComponent::SCRIPT:
      result += indentText(scripts[script_con], 0, actual_width, TextEnds::NEWLINE);
      script_con++;
      break;
    }
    
    // Determine how to handle blank lines between components
    if (i < component_count - 1) {
      switch (sources[i]) {
      case SectionComponent::NARRATIVE:
      case SectionComponent::LIST:
        result += strhide + "\n";
        break;
      case SectionComponent::TABLE:
      case SectionComponent::SCRIPT:
        switch (sources[i + 1]) {
        case SectionComponent::NARRATIVE:
        case SectionComponent::LIST:
        case SectionComponent::TABLE:
          result += "\n";
          break;
        case SectionComponent::SCRIPT:
          break;
        }
        break;
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string SectionContents::sectionAsString(const int alt_section_id, const int alt_subsection_id,
                                             const int list_indent, const int nested_list_indent,
                                             const ListEnumeration list_numbering,
                                             const ListEnumeration nested_list_numbering) const {
  return sectionAsString(style, -1, alt_section_id, alt_subsection_id, list_indent,
                         nested_list_indent, list_numbering, nested_list_numbering);
}

//-------------------------------------------------------------------------------------------------
std::string SectionContents::sectionAsString(const int list_indent, const int nested_list_indent,
                                             const ListEnumeration list_numbering,
                                             const ListEnumeration nested_list_numbering) const {
  return sectionAsString(style, -1, -1, -1, list_indent, nested_list_indent, list_numbering,
                         nested_list_numbering);
}

//-------------------------------------------------------------------------------------------------
void SectionContents::printSection(std::ofstream *foutp, const int list_indent,
                                   const int nested_list_indent,
                                   const ListEnumeration list_numbering,
                                   const ListEnumeration nested_list_numbering) const {
  const std::string sct = sectionAsString(list_indent, nested_list_indent, list_numbering,
                                          nested_list_numbering);
  foutp->write(sct.data(), sct.size());
}

//-------------------------------------------------------------------------------------------------
void SectionContents::printSection(std::ostream *foutp, const int list_indent,
                                   const int nested_list_indent,
                                   const ListEnumeration list_numbering,
                                   const ListEnumeration nested_list_numbering) const {
  const std::string sct = sectionAsString(list_indent, nested_list_indent, list_numbering,
                                          nested_list_numbering);
  foutp->write(sct.data(), sct.size());
}

//-------------------------------------------------------------------------------------------------
void SectionContents::printSection(const std::string &alt_file_name,
                                   const PrintSituation alt_expectation, const int list_indent,
                                   const int nested_list_indent,
                                   const ListEnumeration list_numbering,
                                   const ListEnumeration nested_list_numbering) const {
  std::ofstream foutp = openOutputFile(alt_file_name, alt_expectation, "write a result section to "
                                       "a named file");
  printSection(&foutp, list_indent, nested_list_indent, list_numbering, nested_list_numbering);
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
void SectionContents::printSection(const int list_indent, const int nested_list_indent,
                                   const ListEnumeration list_numbering,
                                   const ListEnumeration nested_list_numbering) const {
  printSection(file_name, expectation, list_indent, nested_list_indent, list_numbering,
               nested_list_numbering);
}

//-------------------------------------------------------------------------------------------------
void SectionContents::setTitle(const std::string &title_in) {
  title = title_in;
}

//-------------------------------------------------------------------------------------------------
void SectionContents::designateSubsection(const bool subsection_in) {
  subsection = subsection_in;
}

//-------------------------------------------------------------------------------------------------
void SectionContents::setSectionDetails(const int section_id_in) {
  section_id = section_id_in;
}

//-------------------------------------------------------------------------------------------------
void SectionContents::setSectionDetails(const ListEnumeration section_marking_in) {
  section_marking = section_marking_in;
}

//-------------------------------------------------------------------------------------------------
void SectionContents::setSectionDetails(const int section_id_in,
                                       const ListEnumeration section_marking_in) {
  section_id = section_id_in;
  section_marking = section_marking_in;
}

//-------------------------------------------------------------------------------------------------
void SectionContents::reserve(const SectionComponent item_kind, const int item_count) {
  switch (item_kind) {
  case SectionComponent::NARRATIVE:
    narrations.reserve(item_count);
    break;
  case SectionComponent::LIST:
    lists.reserve(item_count);
    break;
  case SectionComponent::TABLE:
    tables.reserve(item_count);
    break;
  case SectionComponent::SCRIPT:
    scripts.reserve(item_count);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void SectionContents::addNarration(const std::string &narration_in) {
  narrations.push_back(narration_in);
  sources.push_back(SectionComponent::NARRATIVE);
  component_count += 1;
}

//-------------------------------------------------------------------------------------------------
void SectionContents::addNarration(const TextFile &narration_in, const TextEnds line_endings) {
  addNarration(narration_in.toString());
}

//-------------------------------------------------------------------------------------------------
void SectionContents::addList(const OrderedList &list_in) {
  lists.push_back(list_in);
  sources.push_back(SectionComponent::LIST);
  component_count += 1;
}

//-------------------------------------------------------------------------------------------------
void SectionContents::addTable(const ReportTable &table_in) {
  tables.push_back(table_in);
  sources.push_back(SectionComponent::TABLE);
  component_count += 1;
}

//-------------------------------------------------------------------------------------------------
void SectionContents::addScript(const std::string &script_in) {
  scripts.push_back(script_in);
  sources.push_back(SectionComponent::SCRIPT);
  component_count += 1;
}

//-------------------------------------------------------------------------------------------------
void SectionContents::addScript(const TextFile &script_in, const TextEnds line_endings) {
  addScript(script_in.toString(line_endings));
}

//-------------------------------------------------------------------------------------------------
void SectionContents::validateComponentIndex(const int index) const {
  if (index < 0 || index >= component_count) {
    rtErr("Index " + std::to_string(index) + " is invalid for a collection of " +
          std::to_string(component_count) + " components.", "SectionContents",
          "validateComponentIndex");
  }
}

//-------------------------------------------------------------------------------------------------
std::string printAllSections(const std::vector<SectionContents> &info, const OutputSyntax style,
                             const ListEnumeration numbering,
                             const ListEnumeration nested_numbering, const int width,
                             const int list_indent, const int nested_list_indent,
                             const ListEnumeration list_numbering,
                             const ListEnumeration nested_list_numbering) {
  const char hide = commentSymbol(style);
  const int splash_width = (width >= 0) ? width : (info.size() > 0) ? info[0].getWidth() :
                                                                      default_output_file_width;
  std::string result = protectText(stormmSplash(splash_width - 2), hide, splash_width);
  int curr_section = -1;
  int curr_subsection = 0;
  bool on_subs = false;
  const std::string section_spacer("\n");
  result += section_spacer;
  for (size_t i = 0; i < info.size(); i++) {
    const int i_width = (width >= 0) ? width : info[i].getWidth();

    // Shall the section counter or the subsection counter advance?  Set the section information.
    if (info[i].isSubsection()) {
      if (on_subs) {
        curr_subsection++;
      }
      else {
        curr_subsection = 0;
      }
      on_subs = true;
    }
    else {
      on_subs = false;
      curr_section++;
    }

    // Print the section as appropriate.
    result += info[i].sectionAsString(style, i_width, curr_section, curr_subsection,
                                      list_indent, nested_list_indent, list_numbering,
                                      nested_list_numbering);
    result += section_spacer;
  }
  result += protectText(stormmWatermark(splash_width - 2), hide, splash_width, TextEnds::AS_IS);
  return result;
}

//-------------------------------------------------------------------------------------------------
void printAllSections(std::ofstream *foutp, const std::vector<SectionContents> &info,
                      const OutputSyntax style, const ListEnumeration numbering,
                      const ListEnumeration nested_numbering, const int width,
                      const int list_indent, const int nested_list_indent,
                      const ListEnumeration list_numbering,
                      const ListEnumeration nested_list_numbering) {
  
  // The output file may be voluminous.  Print it in stages, avoiding excessive duplication of
  // data in memory as strings.
  const char hide = commentSymbol(style);
  const int splash_width = (width >= 0) ? width : (info.size() > 0) ? info[0].getWidth() :
                                                                      default_output_file_width;
  printProtectedText(stormmSplash(splash_width - 2), hide, foutp, splash_width, TextEnds::NEWLINE);
  int curr_section = -1;
  int curr_subsection = 0;
  bool on_subs = false;
  const std::string section_spacer("\n");
  foutp->write(section_spacer.data(), section_spacer.size());
  for (size_t i = 0; i < info.size(); i++) {
    const int i_width = (width >= 0) ? width : info[i].getWidth();
    
    // Shall the section counter or the subsection counter advance?  Set the section information.
    if (info[i].isSubsection()) {
      if (on_subs) {
        curr_subsection++;
      }
      else {
        curr_subsection = 0;
      }
      on_subs = true;
    }
    else {
      on_subs = false;
      curr_section++;
    }

    // Print the section as appropriate.
    const std::string isect = info[i].sectionAsString(style, i_width, curr_section,
                                                      curr_subsection, list_indent,
                                                      nested_list_indent, list_numbering,
                                                      nested_list_numbering);
    foutp->write(isect.data(), isect.size());
    foutp->write(section_spacer.data(), section_spacer.size());
  }
  printProtectedText(stormmWatermark(splash_width - 2), hide, foutp, splash_width,
                     TextEnds::AS_IS);
}

//-------------------------------------------------------------------------------------------------
void printAllSections(const std::string &file_name, const PrintSituation expectation,
                      const std::vector<SectionContents> &info, const OutputSyntax style,
                      const ListEnumeration numbering, const ListEnumeration nested_numbering,
                      const int width, const int list_indent, const int nested_list_indent,
                      const ListEnumeration list_numbering,
                      const ListEnumeration nested_list_numbering) {
  std::ofstream foutp = openOutputFile(file_name, expectation, "open a file for printing a series "
                                       "of output sections");
  printAllSections(&foutp, info, style, numbering, nested_numbering, width, list_indent,
                   nested_list_indent, list_numbering, nested_list_numbering);
  foutp.close();
}

} // namespace review
} // namespace stormm
