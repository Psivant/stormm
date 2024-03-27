// -*-c++-*-
#ifndef STORMM_ORDERED_LIST_H
#define STORMM_ORDERED_LIST_H

#include <string>
#include <vector>
#include "copyright.h"
#include "FileManagement/file_enumerators.h"
#include "Parsing/textfile.h"
#include "reporting_enumerators.h"
#include "summary_file.h"

namespace stormm {
namespace review {

using diskutil::PrintSituation;
using parse::TextFile;

/// \brief The ability to express multiple items as an ordered list in formatted output is a great
///        luxury in a command-line code.  This object will accept multiple text strings and
///        manage their presentation as a list with appropriate identifiers.
class OrderedList {
public:

  /// \brief The constructor takes an indication of the enumerative style (e.g. bullet points with
  ///        a specific symbol, numbers, letters), a proposed length to the list, indentation, and
  ///        an optional nesting style (with nested identation).  Items must be subsequently added
  ///        to the list in order to fill out primary and nested items.
  ///
  /// \param reserved_length         Number of list items to prepare for catalogging
  /// \param nested_reserved_length  Number of nested items to prepare for catalogging
  /// \{
  OrderedList(ListEnumeration style_in, int reserved_length = 0, int indentation_in = 2,
              char bullet_in = '*', ListEnumeration nested_style_in = ListEnumeration::BULLET,
              int nested_reserved_length = 0, int nested_indentation_in = 4,
              char nested_bullet_in = '-', int format_width_in = default_output_file_width);

  OrderedList(ListEnumeration style_in, ListEnumeration nested_style_in);  
  /// \}

  /// \brief With no const members or pointers to repair and all Standard Template Library
  ///        components, the OrderedList can be copied or moved using default constructors and
  ///        assignement operators.
  /// \{
  OrderedList(const OrderedList &original) = default;
  OrderedList(OrderedList &&original) = default;
  OrderedList& operator=(const OrderedList &original) = default;
  OrderedList& operator=(OrderedList &&original) = default;
  /// \}

  /// \brief Get the number of items in the list
  int getItemCount() const;

  /// \brief Get the number of nested items.
  ///
  /// Overloaded:
  ///   - Get the total number of nested items in the entire list
  ///   - Get the number of nested items for an item of a particular index
  ///
  /// \param item_index  Index of the item of interest
  /// \{
  int getNestedItemCount() const;
  int getNestedItemCount(int item_index) const;
  /// \}

  /// \brief Get one of the items in the list, formatted with indentation and the marker symbol.
  ///
  /// \param item_index  Index of the item of interest
  /// \param width       The width at which to print the formatted text, including indentation and
  ///                    marker characters.  The default of -1 triggers use of the object's
  ///                    format_width member variable.
  /// \param alt_indent         An alternate indentation for the nested list item (a value of < 0
  ///                           will defer to the object's internal setting)
  /// \param alt_style   An alternate style to apply to the item numbering (a value of NONE will
  ///                    defer to the object's internal setting)
  std::string getItem(int item_index, int width = -1, int alt_indent = -1,
                      ListEnumeration alt_style = ListEnumeration::NONE) const;

  /// \brief Get one of the nested items associates with a main list item, formated with
  ///        indentation and the nested marker symbol.
  ///
  /// \param item_index         Index of the main list item
  /// \param nested_item_index  Index of the nested item from within the sub-list
  /// \param width              The width at which to print the formatted text, including
  ///                           indentation and marker characters.  The default of -1 triggers use
  ///                           of the object's format_width member variable.
  /// \param alt_indent         An alternate indentation for the nested list item (a value of < 0
  ///                           will defer to the object's internal setting)
  /// \param alt_style          An alternate style to apply to the item numbering (a value of NONE
  ///                           will defer to the object's internal setting)
  std::string getNestedItem(int item_index, int nested_item_index, int width = -1,
                            int alt_indent = -1,
                            ListEnumeration alt_style = ListEnumeration::NONE) const;

  /// \brief Print all items and nested items in the list, in order.  Each overloads of this
  ///        function delegates work to one of the prior overloads.
  ///
  /// Overloaded:
  ///   - Print the output to a string (this output will not be protected by any comment character
  ///     to prevent interpretation by a program reading it)
  ///   - Print the output to a TextFile object as it would be output to an open file stream,
  ///     using character protection for various output formats.
  ///   - Print the output to a file stream using the indicated style (this dertermines the comment
  ///     protection character used to shield the list and its items from interpretation by the
  ///     chosen plotting program)
  ///   - Print the output to a named file so long as it is found in the state indicated
  ///
  /// \param width        The width at which to print the formatted text, including indentation and
  ///                     marker characters.  The default of -1 triggers use of the object's
  ///                     format_width member variable.
  /// \param foutp        File stream pointer for the output
  /// \param file_name    Name of the file to open and write the table into
  /// \param expectation  The state in which the output file is to be found to permit writing
  /// \param style        One of several options for formatting the non-protected portion of the
  ///                     table to be amenable to certain plotting programs
  /// \{
  std::string printList(int width = default_output_file_width, int alt_indent = -1,
                        int alt_nested_indent = -1,
                        ListEnumeration alt_style = ListEnumeration::NONE,
                        ListEnumeration alt_nested_style = ListEnumeration::NONE) const;

  TextFile printList(OutputSyntax style, int width, int alt_indent = -1,
                     int alt_nested_indent = -1,
                     ListEnumeration alt_marking = ListEnumeration::NONE,
                     ListEnumeration alt_nested_marking = ListEnumeration::NONE) const;
  
  void printList(std::ofstream *foutp, OutputSyntax style, int width = default_output_file_width,
                 int alt_indent = -1, int alt_nested_indent = -1,
                 ListEnumeration alt_marking = ListEnumeration::NONE,
                 ListEnumeration alt_nested_marking = ListEnumeration::NONE) const;

  void printList(const std::string &file_name, PrintSituation expectation, OutputSyntax style,
                 int width = default_output_file_width, int alt_indent = -1,
                 int alt_nested_indent = -1, ListEnumeration alt_marking = ListEnumeration::NONE,
                 ListEnumeration alt_nested_marking = ListEnumeration::NONE) const;
  /// \}

  /// \brief Set the main list itemization style.  Recalculate the maximum marker lengths,
  ///        anticipating format changes.
  ///
  /// \param style_in  The style to apply
  void setStyle(ListEnumeration style_in);

  /// \brief Set the nested itemization style.  Recalculate the maximum marker lengths,
  ///        anticipating format changes.
  ///
  /// \param nested_style_in  The style to apply
  void setNestedStyle(ListEnumeration nested_style_in);

  /// \brief Set the bullet character.  This produces an error if the LineEnumeration type of main
  ///        list items is not BULLET.
  ///
  /// \param bullet_in  The bullet character to set
  void setBullet(char bullet_in);

  /// \brief Set the nested bullet character.  This produces an error if the LineEnumeration type
  ///        of main list items is not BULLET.
  ///
  /// \param nested_bullet_in  The bullet character to set
  void setNestedBullet(char nested_bullet_in);

  /// \brief Set the format width for all getItem() and getNestedItem() calls.  The width set here
  ///        can be overridden by providing an explicit parameter to either call.
  ///
  /// \param format_width_in  The format width to set
  void setFormatWidth(int format_width_in);
  
  /// \brief Add an item to the list.
  ///
  /// \param item_in     The item to add
  /// \param item_index  The index point at which to add the item (the item currently occupying
  ///                    this index, and all others behind it, will be shifted back by one--it is
  ///                    std::vector insert)
  void addItem(const std::string &item_in, int item_index = -1);

  /// \brief Add an item before a specific member of the current list.
  ///
  /// Overloaded:
  ///   - Identify the main list item by index
  ///   - Identify the main list item by matching the leading part of the string
  ///
  /// \param item_in       The item to add
  /// \param item_index    Index of the main list item to add a new entry in front of
  /// \param current_item  A sequence of characters matching an item in the current list.  The
  ///                      first item matching this sequence in its totality will be taken as the
  ///                      point before which to add the new item.  Failing to find this sequence
  ///                      at the front of any list item produces a runtime error.
  /// \{
  void addItemBefore(const std::string &item_in, int item_index);
  void addItemBefore(const std::string &item_in, const std::string &current_item);
  /// \}

  /// \brief Add an item after a specific member of the current list.
  ///
  /// Overloaded:
  ///   - Identify the main list item by index
  ///   - Identify the main list item by matching the leading part of the string
  ///
  /// \param item_in       The item to add
  /// \param item_index    Index of the main list item to add a new entry in front of
  /// \param current_item  A sequence of characters matching an item in the current list.  The
  ///                      first item matching this sequence in its totality will be taken as the
  ///                      point after which to add the new item.  Failing to find this sequence
  ///                      at the front of any list item produces a runtime error.
  /// \{
  void addItemAfter(const std::string &item_in, int item_index);
  void addItemAfter(const std::string &item_in, const std::string &current_item);
  /// \}

  /// \brief Add a nested item to one of the main list items.
  ///
  /// Overloaded:
  ///   - Identify the main list item by index
  ///   - Identify the main list item by matching the leading part of the string
  ///
  /// \param nested_item_in  The nested item to add
  /// \param current_item    Item in the current main item list (a positive value must be valid).
  ///                        Negative values (including the default) provided here will cause the
  ///                        nested items to be assigned to the most recently added list item.
  /// \{
  void addNestedItem(const std::string &nested_item_in, int current_item = -1);
  void addNestedItem(const std::string &nested_item_in, const std::string &current_item);
  /// \}

private:
  int item_count;                         ///< The number of items in the list
  ListEnumeration style;                  ///< Style by which list items are to be enumerated
  int indentation;                        ///< Number of spaces from the leftmost point on the line
                                          ///<   by which the markers for each list item are to be
                                          ///<   indented (list items will be further justified
                                          ///<   behind the longest possible marker item)
  ListEnumeration nested_style;           ///< Style by which nested items are to be enumerated.
                                          ///<   If the style is NUMERIC or ALPHABETIC, the nested
                                          ///<   items will be enumerated starting from 1.) or a.)
                                          ///<   with each subsequent main list item.
  int nested_indentation;                 ///< Number of spaces from the leftmost point on the line
                                          ///<   by which the markers for each nested item are to
                                          ///<   be indented (list items will be further justified
                                          ///<   behind the longest possible marker item)
  char bullet;                            ///< Character symbol to use for the marker, if the list
                                          ///<   is not enumerated in alpha-numeric symbols
  char nested_bullet;                     ///< Character symbol to use for the nested marker, if
                                          ///<   the list is not enumerated in alpha-numeric
                                          ///<   symbols 
  int max_marker_length;                  ///< Number of characters involved in the longest marker
                                          ///<   symbol
  int max_nested_marker_length;           ///< Number of characters involved in the longest nested
                                          ///<   marker symbol
  int format_width;                       ///< Total width available to the formatted text in the
                                          ///<   ordered list.  This provides a backup in case the
                                          ///<   getItem() and getNestedItem() member functions are
                                          ///<   called without a width parameter.
  std::vector<std::string> items;         ///< Ordered contents of the list, each a word, phrase,
                                          ///<   or even a paragraph
  std::vector<int> nested_item_limits;    ///< The ranges of nested items owned by each main list
                                          ///<   item, a bounds array on nested_item contents
  std::vector<int> nested_item_contents;  ///< Members of the nested_items array owned by each
                                          ///<   main list item
  std::vector<std::string> nested_items;  ///< All nested items for any main list item
  std::vector<int> item_homes;            ///< Indices of items to which each nested item belongs

  /// \brief Validate the index of an item to ensure safe memory access.
  ///
  /// \param item_index  The item index of interest
  /// \param caller      Name of the calling function
  void validateItemIndex(int item_index, const char* caller = nullptr) const;  

  /// \brief Validate the index of a nested item to ensure safe memory access.
  ///
  /// \param nested_item_index  Index of the nested item within the main list item's ownership
  /// \param caller             Name of the calling function
  void validateNestedItemIndex(int item_index, int nested_item_index,
                               const char* caller = nullptr) const;
  
  /// \brief Locate an item in the list (not a nested item) containing a specific sequence of
  ///        characters.  The first such match will be taken as the item of interest and its index
  ///        returned.  If no match is found, this function returns the number of items currently
  ///        in the main item list.
  ///
  /// \param current_item  A sequence of characters matching an item in the current list.  Items
  ///                      with additional contents will also mathc provided their characters in
  ///                      the front match this string.
  int locateItem(const std::string &current_item) const;

  /// \brief Compute the maximum length of markers for the main list items and nested items.
  ///        These values are recalculated every time an item or nested item is added to the list,
  ///        but are inexpensive.
  void calculateMarkerLengths();
};
  
} // namespace review
} // namespace stormm

#endif
