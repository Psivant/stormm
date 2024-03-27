#include <cstring>
#include "copyright.h"
#include "FileManagement/file_util.h"
#include "Math/series_ops.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "error_format.h"
#include "ordered_list.h"

namespace stormm {
namespace review {

using diskutil::openOutputFile;
using errors::terminalFormat;
using parse::addLeadingWhiteSpace;
using parse::addTailingWhiteSpace;
using parse::alphabetNumber;
using parse::TextOrigin;
using stmath::indexingArray;

//-------------------------------------------------------------------------------------------------
OrderedList::OrderedList(const ListEnumeration style_in, const int reserved_length,
                         const int indentation_in, const char bullet_in,
                         const ListEnumeration nested_style_in, const int nested_reserved_length,
                         const int nested_indentation_in, const char nested_bullet_in,
                         const int format_width_in) :
    item_count{0}, style{style_in}, indentation{indentation_in}, nested_style{nested_style_in},
    nested_indentation{nested_indentation_in}, bullet{bullet_in}, nested_bullet{nested_bullet_in},
    max_marker_length{0}, max_nested_marker_length{0}, format_width{format_width_in},
    items{}, nested_item_limits{}, nested_item_contents{}, nested_items{}, item_homes{}
{
  // Reserve space for items and nested items, as requested (this does not yet increment the
  // item counter, which remains the official tracker of the numbers of items and nested items)
  items.reserve(reserved_length);
  nested_items.reserve(nested_reserved_length);
  nested_item_limits.reserve(reserved_length + 1);
  nested_item_contents.reserve(nested_reserved_length);
  item_homes.reserve(nested_reserved_length);

  // Initiate the nested item limits with the obligatory zero
  nested_item_limits.push_back(0);
}

//-------------------------------------------------------------------------------------------------
OrderedList::OrderedList(const ListEnumeration style_in, const ListEnumeration nested_style_in) :
    OrderedList(style_in, 0, 2, '*', nested_style_in)
{}
  
//-------------------------------------------------------------------------------------------------
int OrderedList::getItemCount() const {
  return item_count;
}

//-------------------------------------------------------------------------------------------------
int OrderedList::getNestedItemCount() const {
  return nested_items.size();
}

//-------------------------------------------------------------------------------------------------
int OrderedList::getNestedItemCount(const int item_index) const {
  if (item_index < 0 || item_index >= item_count) {
    rtErr("Item index " + std::to_string(item_index) + " is invalid for a list of " +
          std::to_string(item_count) + " items.", "OrderedList", "getNestedItemCount");
  }
  return nested_item_limits[item_index + 1] - nested_item_limits[item_index];
}

//-------------------------------------------------------------------------------------------------
std::string OrderedList::getItem(const int item_index, const int width, const int alt_indent,
                                 const ListEnumeration alt_style) const {
  validateItemIndex(item_index, "getItem");
  std::string marker;
  const ListEnumeration style_in_use = (alt_style == ListEnumeration::NONE) ? style : alt_style;
  switch (style_in_use) {
  case ListEnumeration::BULLET:
    marker.resize(1);
    marker[0] = bullet;
    break;
  case ListEnumeration::NUMBERED:
    marker = std::to_string(item_index + 1) + ")";
    break;
  case ListEnumeration::ALPHABETIC:
    marker = alphabetNumber(item_index) + ")";
    break;
  case ListEnumeration::ROMAN:
    marker = toRoman(item_index + 1) + ")";
    break;
  case ListEnumeration::NONE:
    break;
  }
  marker = addTailingWhiteSpace(marker, max_marker_length);
  const int indent_in_use = (alt_indent < 0) ? indentation : alt_indent;
  const int width_in_use  = (width <= 0) ? format_width : width;
  return indentText(items[item_index], marker, indent_in_use, width_in_use, false,
                    TextEnds::NEWLINE);
}

//-------------------------------------------------------------------------------------------------
std::string OrderedList::getNestedItem(const int item_index, const int nested_item_index,
                                       const int width, const int alt_indent,
                                       const ListEnumeration alt_style) const {
  validateNestedItemIndex(item_index, nested_item_index, "getNestedItem");
  const int n_idx = nested_item_contents[nested_item_limits[item_index] + nested_item_index];
  std::string marker;
  const ListEnumeration style_in_use = (alt_style == ListEnumeration::NONE) ? nested_style :
                                                                              alt_style;
  switch (style_in_use) {
  case ListEnumeration::BULLET:
    marker.resize(1);
    marker[0] = nested_bullet;
    break;
  case ListEnumeration::NUMBERED:
    marker = std::to_string(nested_item_index + 1) + ")";
    break;
  case ListEnumeration::ALPHABETIC:
    marker = alphabetNumber(nested_item_index) + ")";
    break;
  case ListEnumeration::ROMAN:
    marker = toRoman(nested_item_index + 1) + ")";
    break;
  case ListEnumeration::NONE:
    break;
  }
  marker = addTailingWhiteSpace(marker, max_nested_marker_length);
  const int indent_in_use = (alt_indent < 0) ? nested_indentation : alt_indent;
  const int width_in_use  = (width <= 0) ? format_width : width;
  return indentText(nested_items[n_idx], marker, indent_in_use, width_in_use, false,
                    TextEnds::NEWLINE);
}

//-------------------------------------------------------------------------------------------------
std::string OrderedList::printList(const int width, const int alt_indent,
                                   const int alt_nested_indent, const ListEnumeration alt_style,
                                   const ListEnumeration alt_nested_style) const {
  std::string result;
  for (int i = 0; i < item_count; i++) {
    result += getItem(i, width, alt_indent, alt_style);
    const int j_llim = nested_item_limits[i];
    const int j_hlim = nested_item_limits[i + 1];
    for (int j = j_llim; j < j_hlim; j++) {
      result += getNestedItem(i, j - j_llim, width, alt_nested_indent, alt_nested_style);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
TextFile OrderedList::printList(const OutputSyntax style, const int width, const int alt_indent,
                                const int alt_nested_indent, const ListEnumeration alt_marking,
                                const ListEnumeration alt_nested_marking) const {
  if (width <= 2) {
    rtErr("A line width of " + std::to_string(width) + " characters is not enough space to print "
          "a protected list.", "OrderedList", "printList");
  }
  const char protector = commentSymbol(style);
  const std::string content = protectText(printList(width - 2, alt_indent, alt_nested_indent,
                                                    alt_marking, alt_nested_marking), protector,
                                          width, TextEnds::AS_IS);
  return TextFile(content, TextOrigin::RAM);
}

//-------------------------------------------------------------------------------------------------
void OrderedList::printList(std::ofstream *foutp, const OutputSyntax style, const int width,
                            const int alt_indent, const int alt_nested_indent,
                            const ListEnumeration alt_marking,
                            const ListEnumeration alt_nested_marking) const {
  if (width <= 2) {
    rtErr("A line width of " + std::to_string(width) + " characters is not enough space to print "
          "a protected list.", "OrderedList", "printList");
  }
  const char protector = commentSymbol(style);
  const std::string result = protectText(printList(width - 2, alt_indent, alt_nested_indent,
                                                   alt_marking, alt_nested_marking), protector,
                                         width, TextEnds::AS_IS);
  foutp->write(result.data(), result.size());
}

//-------------------------------------------------------------------------------------------------
void OrderedList::printList(const std::string &file_name, const PrintSituation expectation,
                            const OutputSyntax style, const int width,
                            const int alt_indent, const int alt_nested_indent,
                            const ListEnumeration alt_marking,
                            const ListEnumeration alt_nested_marking) const {
  std::ofstream foutp = openOutputFile(file_name, expectation, "print an ordered list");
  printList(&foutp, style, width, alt_indent, alt_nested_indent, alt_marking, alt_nested_marking);
  foutp.close();
}
  
//-------------------------------------------------------------------------------------------------
void OrderedList::setStyle(const ListEnumeration style_in) {
  style = style_in;
  calculateMarkerLengths();  
}

//-------------------------------------------------------------------------------------------------
void OrderedList::setNestedStyle(const ListEnumeration nested_style_in) {
  nested_style = nested_style_in;
  calculateMarkerLengths();  
}

//-------------------------------------------------------------------------------------------------
void OrderedList::setBullet(const char bullet_in) {
  style = ListEnumeration::BULLET;
  bullet = bullet_in;
}

//-------------------------------------------------------------------------------------------------
void OrderedList::setNestedBullet(const char nested_bullet_in) {
  nested_style = ListEnumeration::BULLET;
  nested_bullet = nested_bullet_in;
}

//-------------------------------------------------------------------------------------------------
void OrderedList::setFormatWidth(const int format_width_in) {
  format_width = format_width_in;
}
  
//-------------------------------------------------------------------------------------------------
void OrderedList::addItem(const std::string &item_in, const int item_index) {
  if (item_index < 0 || item_index == item_count) {
    items.push_back(item_in);
    item_count++;
    nested_item_limits.push_back(nested_item_limits.back());
  }
  else if (item_index < item_count) {

    // Perform the vector insertion
    items.insert(items.begin() + item_index, item_in);
    item_count++;
    
    // Repair any nested item indexing
    const size_t nested_item_count = nested_items.size();
    for (size_t i = 0; i < nested_item_count; i++) {
      if (item_homes[i] >= item_index) {
        item_homes[i] += 1;
      }
    }
    nested_item_limits.insert(nested_item_limits.begin() + item_index,
                              nested_item_limits[item_index]);
  }
  else {
    rtErr("Item index " + std::to_string(item_index) + " is invalid for an ordered list of "
          "current length " + std::to_string(item_count) + " items.", "OrderedList", "addItem");
  }
  calculateMarkerLengths();
}

//-------------------------------------------------------------------------------------------------
void OrderedList::addItemBefore(const std::string &item_in, const int item_index) {
  validateItemIndex(item_index, "addItemBefore");
  addItem(item_in, item_index);
}

//-------------------------------------------------------------------------------------------------
void OrderedList::addItemBefore(const std::string &item_in, const std::string &current_item) {

  // Seek out the item matching the given string
  const int loc = locateItem(current_item);
  if (loc < item_count) {
    addItem(item_in, loc);
  }
  else {
    rtErr("No item matching the initial characters \"" + current_item + "\" was found.",
          "OrderedList", "addItemBefore");
  }
}

//-------------------------------------------------------------------------------------------------
void OrderedList::addItemAfter(const std::string &item_in, const int item_index) {
  validateItemIndex(item_index, "addItemAfter");
  addItem(item_in, item_index + 1);
}

//-------------------------------------------------------------------------------------------------
void OrderedList::addItemAfter(const std::string &item_in, const std::string &current_item) {

  // Seek out the item matching the given string
  const int loc = locateItem(current_item);
  if (loc < item_count) {
    addItem(item_in, loc + 1);
  }
  else {
    rtErr("No item matching the initial characters \"" + current_item + "\" was found.",
          "OrderedList", "addItemAfter");
  }
}

//-------------------------------------------------------------------------------------------------
void OrderedList::addNestedItem(const std::string &nested_item_in, const int current_item) {
  if (current_item >= item_count) {
    rtErr("Index " + std::to_string(current_item) + " is invalid for a collection of " +
          std::to_string(item_count) + " items.", "OrderedList", "addNestedItem");    
  }
  nested_items.push_back(nested_item_in);
  nested_item_contents.resize(nested_items.size());
  if (current_item >= 0) {
    item_homes.push_back(current_item);
  }
  else {
    item_homes.push_back(item_count - 1);
  }
  indexingArray<int, int>(item_homes, &nested_item_contents, &nested_item_limits, item_count);
  calculateMarkerLengths();
}    

//-------------------------------------------------------------------------------------------------
void OrderedList::addNestedItem(const std::string &nested_item_in,
                                const std::string &current_item) {
  const int loc = locateItem(current_item);
  if (loc == item_count) {
    rtErr("No item matching the initial characters \"" + current_item + "\" was found.",
          "OrderedList", "addNestedItem");    
  }
  addNestedItem(nested_item_in, loc);
}    

//-------------------------------------------------------------------------------------------------
  void OrderedList::validateItemIndex(const int item_index, const char* caller) const {
  if (item_index < 0 || item_index >= item_count) {
    rtErr("Index " + std::to_string(item_index) + " is out of bounds for a list of " +
          std::to_string(item_count) + " items.", "OrderedList", caller);
  }
}

//-------------------------------------------------------------------------------------------------
void OrderedList::validateNestedItemIndex(const int item_index,
                                          const int nested_item_index, const char* caller) const {
  validateItemIndex(item_index, "validateNestedItemIndex");
  const int nested_item_count = nested_item_limits[item_index + 1] -
                                nested_item_limits[item_index];
  if (nested_item_index < 0 || nested_item_index >= nested_item_count) {
    rtErr("Nested item index " + std::to_string(nested_item_index) + " is out of bounds for a "
          "main list item member with " + std::to_string(nested_item_count) + " nested items.",
          "OrderedList", caller);
  }
}

//-------------------------------------------------------------------------------------------------
int OrderedList::locateItem(const std::string &current_item) const {
  const size_t cis = current_item.size();
  for (int i = 0; i < item_count; i++) {
    if (current_item.size() == 0 && items[i].size() == 0) {
      return i;
    }
    else {
        if (cis > 0 && cis <= items[i].size() && current_item[0] == items[i][0] &&
          current_item[cis - 1] == items[i][cis - 1] &&
          strncmp(current_item.data(), items[i].data(), cis) == 0) {
          return i;
      }
    }
  }
  return item_count;
}

//-------------------------------------------------------------------------------------------------
void OrderedList::calculateMarkerLengths() {

  // Calculate the longest marker in the main list items
  switch (style) {
  case ListEnumeration::BULLET:
    max_marker_length = 1;
    break;
  case ListEnumeration::NUMBERED:
    max_marker_length = static_cast<int>(ceil(log10(item_count + 1))) + 1;
    break;
  case ListEnumeration::ALPHABETIC:
    {
      int combo_inc = 26;
      int nchar = 1;
      int combo_count = combo_inc;
      while (item_count > combo_count) {
        combo_inc *= 26;
        nchar++;
        combo_count += combo_inc;
      }
      max_marker_length = nchar + 1;
    }
    break;
  case ListEnumeration::ROMAN:
    {
      max_marker_length = 0;
      for (int i = 1; i <= item_count; i++) {
        max_marker_length = std::max(max_marker_length, static_cast<int>(toRoman(i).size()) + 1);
      }
    }
    break;
  case ListEnumeration::NONE:
    max_marker_length = 0;
    break;
  }

  // Calculate the longest marker in the nested list items
  switch (nested_style) {
  case ListEnumeration::BULLET:
    max_nested_marker_length = 1;
    break;
  case ListEnumeration::NUMBERED:
    for (int i = 0; i < item_count; i++) {
      const int nsub = nested_item_limits[i + 1] - nested_item_limits[i];
      if (nsub > 0) {
        max_nested_marker_length = std::max(static_cast<int>(ceil(log10(nsub + 1))) + 1,
                                            max_nested_marker_length);
      }
    }
    break;
  case ListEnumeration::ALPHABETIC:
    for (int i = 0; i < item_count; i++) {
      const int nsub = nested_item_limits[i + 1] - nested_item_limits[i];
      int combo_inc = 26;
      int nchar = 1;
      int combo_count = combo_inc;
      while (nsub > combo_count) {
        combo_inc *= 26;
        nchar++;
        combo_count += combo_inc;
      }
      max_nested_marker_length = std::max(nchar + 1, max_nested_marker_length);
    }
    break;
  case ListEnumeration::ROMAN:
    {
      max_nested_marker_length = 0;
      for (int i = 0; i < item_count; i++) {
        const int nsub = nested_item_limits[i + 1] - nested_item_limits[i];
        for (int j = 1; j <= nsub; j++) {
          max_nested_marker_length = std::max(max_nested_marker_length,
                                              static_cast<int>(toRoman(j).size()) + 1);
        }
      }
    }
    break;
  case ListEnumeration::NONE:
    max_nested_marker_length = 0;
    break;
  }
}

} // namespace review
} // namespace stormm
