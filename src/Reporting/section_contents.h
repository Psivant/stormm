// -*-c++-*-
#ifndef STORMM_SECTION_CONTENTS_H
#define STORMM_SECTION_CONTENTS_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "copyright.h"
#include "FileManagement/file_enumerators.h"
#include "ordered_list.h"
#include "reporting_enumerators.h"
#include "report_table.h"
#include "summary_file.h"

namespace stormm {
namespace review {

using diskutil::PrintSituation;

/// \brief Any output file is expected to break into sections, which should be numbered in order
///        and formatted consistently.  This object will collect the elements of a section and
///        provide methods for automatically writing them to an arbitrary output file.
class SectionContents {
public:
  
  /// \brief The constructor takes critical components of the output format: target, format width,
  ///        syntax.  Specifying whether the information constitutes a subsection or not is
  ///        optional, and as with other such details this can be specified after creation of the
  ///        object.
  /// \{
  SectionContents(const std::string &title_in, const std::string &file_name_in, int width_in,
                  bool subsection_in, int section_id_in, int subsection_id_in,
                  ListEnumeration section_marking_in, ListEnumeration subsection_marking_in,
                  OutputSyntax style_in, double table_precision_in);

  SectionContents(const std::string &title_in = std::string(""),
                  const std::string &file_name_in = std::string(""),
                  int width_in = default_output_file_width,
                  OutputSyntax style_in = OutputSyntax::STANDALONE,
                  double table_precision_in = 1.0e-6);
  /// \}

  /// \brief With no const members or pointers to repair and all Standard Template Library
  ///        components, the SectionContents can be copied or moved using default constructors and
  ///        assignement operators.
  /// \{
  SectionContents(const SectionContents &original) = default;
  SectionContents(SectionContents &&original) = default;
  SectionContents& operator=(const SectionContents &original) = default;
  SectionContents& operator=(SectionContents &&original) = default;
  /// \}

  /// \brief Get the name of the output file.
  const std::string& getOutputFileName();

  /// \brief Get the target file width, the maximum number of characters on a typical line.
  int getWidth() const;

  /// \brief Get the output file syntax.
  OutputSyntax getSyntax() const;

  /// \brief Get the number of components in the object
  ///
  /// Overloaded:
  ///   - Get the total number of components
  ///   - Get the number of a specific type of component
  ///
  /// \param kind  The type of component to seek
  /// \{
  int getComponentCount() const;
  int getComponentCount(SectionComponent kind) const;
  /// \}

  /// \brief Get the type of one of the components of the section, based on the index in which
  ///        the components will be listed out.
  ///
  /// \param index  The component index of interest
  SectionComponent getComponentType(int index) const;
  
  /// \brief Get an indication of whether the contents constitute a subsection.
  bool isSubsection() const;
  
  /// \brief Produce a string containing the section contents.
  ///
  /// Overloaded:
  ///   - Provide alternate section and subsection identification numbers
  ///   - Rely on the object's internal section and subsection numbers
  ///
  /// \param alt_style              Alternate formatting style to use (this can help enforce a
  ///                               consistent style if printing from a central output management
  ///                               function)
  /// \param alt_section_id         Alternate section index (to be translated into the enumeration
  ///                               scheme used by the object).  Any value < 0 here will defer to
  ///                               the object's internal subsection index.  Again, this can be
  ///                               useful when printing from central output management functions.
  /// \param alt_subsection_id      Alternate subsection index (to be translated into the
  ///                               enumeration scheme used by the object). Any value < 0 here will
  ///                               defer to the object's internal subsection index.
  /// \param list_indent            Number of white space characters to indent main list elements
  ///                               (a value < 0 here will defer to the indentation indicated by
  ///                               the OrderedList object itself)
  /// \param nested_list_indent     Number of white space characters to indent nested list elements
  ///                               (a value < 0 here will defer to the indentation indicated by
  ///                               the OrderedList object itself)
  /// \param list_numbering         The style in which to increment main list enumerations (a NONE
  ///                               value here will defer to the enumerative scheme within the
  ///                               OrderedList object itself)
  /// \param nested_list_numbering  The style in which to increment nested list enumerations (a
  ///                               NONE value here will defer to the enumerative scheme within the
  ///                               OrderedList object itself)
  /// \{
  std::string sectionAsString(OutputSyntax alt_style, int alt_width = -1, int alt_section_id = -1,
                              int alt_subsection_id = -1, int list_indent = -1,
                              int nested_list_indent = -1,
                              ListEnumeration list_marking = ListEnumeration::NONE,
                              ListEnumeration nested_list_marking = ListEnumeration::NONE) const;

  std::string sectionAsString(int alt_section_id, int alt_subsection_id,
                              int list_indent, int nested_list_indent,
                              ListEnumeration list_numbering = ListEnumeration::NONE,
                              ListEnumeration nested_list_numbering = ListEnumeration::NONE) const;

  std::string sectionAsString(int list_indent = -1, int nested_list_indent = -1,
                              ListEnumeration list_numbering = ListEnumeration::NONE,
                              ListEnumeration nested_list_numbering = ListEnumeration::NONE) const;
  /// \}

  /// \brief Print the section contents to a file.
  ///
  /// Overloaded:
  ///   - Print to a file based on a pointer provided directly to the member function
  ///   - Print to a file based on the pointer stored internally
  ///
  /// Parameter descriptions for this function follow from sectionAsString() above, in addition to:
  /// 
  /// \param foutp                  Alternate file pointer (overrides the file pointer stored
  ///                               internally)
  /// \param alt_file_name          Name of a file to write (the file will be opened and then
  ///                               closed)
  /// \param alt_expectation        State in which the alternate file to write is expected to be
  ///                               found, if the file is indicated by name
  /// \{
  void printSection(std::ofstream *foutp, int list_indent = -1, int nested_list_indent = -1,
                    ListEnumeration list_numbering = ListEnumeration::NONE,
                    ListEnumeration nested_list_numbering = ListEnumeration::NONE) const;
  
  void printSection(std::ostream *foutp, int list_indent = -1, int nested_list_indent = -1,
                    ListEnumeration list_numbering = ListEnumeration::NONE,
                    ListEnumeration nested_list_numbering = ListEnumeration::NONE) const;

  void printSection(const std::string &alt_file_name,
                    PrintSituation alt_expectation = PrintSituation::APPEND,
                    int list_indent = -1, int nested_list_indent = -1,
                    ListEnumeration list_numbering = ListEnumeration::NONE,
                    ListEnumeration nested_list_numbering = ListEnumeration::NONE) const;

  void printSection(int list_indent = -1, int nested_list_indent = -1,
                    ListEnumeration list_numbering = ListEnumeration::NONE,
                    ListEnumeration nested_list_numbering = ListEnumeration::NONE) const;
  /// \}

  /// \brief Set the title of the section (this can be used if the title was not specified in the
  ///        constructor, and will overwritte any existing title).
  ///
  /// \param title_in  The title to set
  void setTitle(const std::string &title_in);
  
  /// \brief Set whether the object contains a subsection or a full section.
  ///
  /// \param subsection_in  Indicate that the information contained in the object is a subsection
  ///                       by setting this to TRUE.  The default setting allows the function to
  ///                       be called with no argument to achieve the obvious effect.
  void designateSubsection(bool subsection_in = true);

  /// \brief Set the section identifier (and, optionally, the numbering style).
  ///
  /// Overloaded:
  ///   - Set the section number only
  ///   - Set the section marking style
  ///   - Set the section number and marking style
  ///
  /// \param section_id_in       The number to be assigned to the section
  /// \param section_marking_in  The manner in which to mark and increment section numbers
  /// \{
  void setSectionDetails(int section_id_in);
  void setSectionDetails(ListEnumeration section_marking_in);
  void setSectionDetails(int section_id_in, ListEnumeration section_marking_in);
  /// \}
  
  /// \brief Reserve space for a number of blocks of narrative, ordered lists, or tables.  This
  ///        can be used to limit memory movement in cases when output sections may contain
  ///        voluminous data.
  ///
  /// \param item_kind   The type of output section item for which to allocate space
  /// \param item_count  The number of such items to prepare the section to hold
  void reserve(SectionComponent item_kind, int item_count);
  
  /// \brief Add a block of narrative to the section.
  ///
  /// Overloaded:
  ///   - Provide the narrative as a string, containing carriage returns where explicit line breaks
  ///     are intended.
  ///   - Provide the narrative as a TextFile object, the various lines of which will be delimited
  ///     by carriage returns or spaces.
  ///
  /// \param narration_in  The new block of narrative
  /// \param line_endings  The manner in which to terminate lines of a TextFile object
  /// \{
  void addNarration(const std::string &narration_in);
  void addNarration(const TextFile &narration_in, TextEnds line_endings = TextEnds::AS_IS);
  /// \}

  /// \brief Add a new ordered list to the section.
  ///
  /// \param list_in  The list to add
  void addList(const OrderedList &list_in);

  /// \brief Add a new table to the section.
  ///
  /// \param table_in  The table to add
  void addTable(const ReportTable &table_in);

  /// \brief Add new commands to the section as script, unprotected by comment characters.
  ///
  /// Overloaded:
  ///   - Provide the script commands as a string, containing carriage returns where explicit
  ///     line breaks are intended.
  ///   - Provide the script as a TextFile object, the various lines of which will be delimited
  ///     by carriage returns or spaces.
  ///
  /// \param script_in  The new block of narrative
  /// \param line_endings  The manner in which to terminate lines of a TextFile object
  /// \{
  void addScript(const std::string &script_in);
  void addScript(const TextFile &script_in, TextEnds line_endings = TextEnds::AS_IS);
  /// \}
  
private:
  std::string title;                      ///< The title of the section (may contain multiple words
                                          ///<   and will be formatted along with the section
                                          ///<   number to the specified output width)
  int component_count;                    ///< The number of components in the section, and the
                                          ///<   length of the sources member array
  int width;                              ///< Target width for the output file.  Lines will be
                                          ///<   held to this width unless other formatting needs
                                          ///<   override it.
  bool subsection;                        ///< Indicate whether the contents represent a section of
                                          ///<   their own (or the beginning of a new section which
                                          ///<   may contain multiple subsections), or if the
                                          ///<   contents are a subsection of some larger entry.
                                          ///<   Unlike nested OrderedList entries, subsection
                                          ///<   contents are not further indented.
  int section_id;                         ///< The identification number of the section for the
                                          ///<   information contained by the object
  int subsection_id;                      ///< The identification number of the subsection for the
                                          ///<   information contained by the object
  ListEnumeration section_marking;        ///< The manner by which section identification numbers
                                          ///<   are to be expressed and incremented
  ListEnumeration subsection_marking;     ///< The manner by which subsection identification
                                          ///<   numbers are to be expressed and incremented
  OutputSyntax style;                     ///< The style in which to print the output, based on
                                          ///<   various plotting programs that the user might
                                          ///<   prefer.  The output generated is intended to be
                                          ///<   fed as a script to such a program.
  std::vector<std::string> narrations;    ///< A collection of narratives, each string a paragraph
                                          ///<   or series of paragraphs awaiting final formatting
                                          ///<   in the width specified.  Elements of this array
                                          ///<   will be selected in order, based on directives
                                          ///<   found in the sources array.
  std::vector<OrderedList> lists;         ///< A collection of lists, typically found after blocks
                                          ///<   of narration.  Elements of this array will be
                                          ///<   selected in order, based on directives found in
                                          ///<   the sources array.
  std::vector<ReportTable> tables;        ///< A collection of tables, each formatted for column
                                          ///<   alignment.  If tables stretch beyond the stated
                                          ///<   output file width, groups of columns will be
                                          ///<   printed in succession to keep within the stated
                                          ///<   width.  Elements of this array will be selected
                                          ///<   in order, based on directives found ina the
                                          ///<   sources array.
  std::vector<std::string> scripts;       ///< A collection of script commands, each string being
                                          ///<   potentially a series of script statements to be
                                          ///<   executed by the plotting program implied by the
                                          ///<   main output syntax.
  std::vector<SectionComponent> sources;  ///< Directives of whether to take from narrations,
                                          ///<   tables, or bulleted list elements as the section
                                          ///<   is printed.
  std::string file_name;                  ///< Name of the file to which the section should be
                                          ///<   written (it is possible to override this setting
                                          ///<   and print to an alternative file).  An empty
                                          ///<   string, or providing an empty string in a call to
                                          ///<   one of the printSection() overloads, triggers
                                          ///<   printing to the terminal window.
  PrintSituation expectation;             ///< State in which the file named by file_name is to be
                                          ///<   found in order for printing to proceed

  /// \brief Validate the index of a section component.  This number may be greater than any of
  ///        the numbers of narrations, lists, or tables individually, but not the total size of
  ///        the sources array.
  ///
  /// \param index  The component index to validate
  void validateComponentIndex(int index) const;
};

/// \brief Print the contents of a list of sections to a specific output file.
///
/// Overloaded:
///   - Specify the file by an open stream pointer
///   - Specify the file by name, with the state in which it is expected to be found
///
/// \param foutp                  Pointer to the file sgtream for output
/// \param info                   An array of the various sections to print.  This function will
///                               check each element of the list for whether it is a section or
///                               subsection, and number the output as appropriate.
/// \param style                  The style in which to print output (this will cascade down and
///                               override the styles of any objects in the list)
/// \param numbering              The manner in which to number successive sections
/// \param nested_numbering       The manner in which to number successive subsections (the
///                               combination of numbering and nested_numbering will determine any
///                               punctuation between section and subsection identifiers, e.g. 2a
///                               or 2.1).  Only numeric, alphabetic, or Roman numeral section
///                               markers are permitted.
/// \param width                  The overall format width to use for all content in the output
/// \param list_indent            The number of indented spaces to use for main list items (this
///                               will override the internal settings of content in the info
///                               parameter if set to a positive, nonzero integer)
/// \param nested_list_indent     The number of indented spaces for list nested items (the
///                               overriding behavior of this parameter is similar to list_indent)
/// \param list_numbering         The manner in which to number main list elements in all content.
///                               This will override the specifications in OrderedList objects
///                               contained in elements of the info parameter if set to anything
///                               other than NONE.
/// \param nested_list_numbering  The manner in which to do nested list numbering (the overriding
///                               behavior of this parameter is similar to list_numbering)
/// \{
std::string printAllSections(const std::vector<SectionContents> &info,
                             OutputSyntax style = OutputSyntax::STANDALONE,
                             ListEnumeration numbering = ListEnumeration::NUMBERED,
                             ListEnumeration nested_numbering = ListEnumeration::NUMBERED,
                             int width = -1, int list_indent = -1, int nested_list_indent = -1,
                             ListEnumeration list_numbering = ListEnumeration::NONE,
                             ListEnumeration nested_list_numbering = ListEnumeration::NONE);

void printAllSections(std::ofstream *foutp, const std::vector<SectionContents> &info,
                      OutputSyntax style = OutputSyntax::STANDALONE,
                      ListEnumeration numbering = ListEnumeration::NUMBERED,
                      ListEnumeration nested_numbering = ListEnumeration::NUMBERED,
                      int width = -1, int list_indent = -1, int nested_list_indent = -1,
                      ListEnumeration list_numbering = ListEnumeration::NONE,
                      ListEnumeration nested_list_numbering = ListEnumeration::NONE);

void printAllSections(const std::string &file_name, const PrintSituation expectation,
                      const std::vector<SectionContents> &info,
                      OutputSyntax style = OutputSyntax::STANDALONE,
                      ListEnumeration numbering = ListEnumeration::NUMBERED,
                      ListEnumeration nested_numbering = ListEnumeration::NUMBERED,
                      int width = -1, int list_indent = -1, int nested_list_indent = -1,
                      ListEnumeration list_numbering = ListEnumeration::NONE,
                      ListEnumeration nested_list_numbering = ListEnumeration::NONE);
/// \}

} // namespace review
} // namespace stormm

#endif
