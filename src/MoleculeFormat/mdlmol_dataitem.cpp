#include "copyright.h"
#include "FileManagement/file_listing.h"
#include "Parsing/parse.h"
#include "mdlmol_dataitem.h"

namespace stormm {
namespace structure {

using constants::CaseSensitivity;
using diskutil::getBaseName;
using parse::NumberFormat;
using parse::strcmpCased;
using parse::stringToChar4;
using parse::strncmpCased;
using parse::separateText;
using parse::verifyContents;

//-------------------------------------------------------------------------------------------------
MdlMolDataItem::MdlMolDataItem(const std::string &item_name_in,
                               const std::string &external_regno_in,
                               const int internal_regno_in, const int maccs_ii_number_in,
                               const uint header_info, const std::vector<std::string> &body_in,
                               const ModificationPolicy mpol, const ExceptionResponse notify) :
    kind{MdlMolDataItemKind::NONE}, tracked_state{StateVariable::ALL_STATES},
    item_name{item_name_in}, output_item_name{item_name_in}, external_regno{external_regno_in},
    internal_regno{internal_regno_in}, maccs_ii_number{maccs_ii_number_in},
    use_internal_regno{((header_info & 0x1) == 1U)},
    use_external_regno{((header_info & 0x2) == 2U)}, use_item_name{((header_info & 0x4) == 4U)},
    use_maccs_ii_number{((header_info & 0x8) == 8U)}, note_archives{((header_info & 0x10) == 16U)},
    body{body_in}
{
  // Check the sanity of the new object and update the output name, if needed
  validateItemName(mpol, notify);
}

//-------------------------------------------------------------------------------------------------
MdlMolDataItem::MdlMolDataItem(const TextFile &tf, const int line_number, int *line_advance,
                               const int compound_line_end, const std::string &title,
                               const ModificationPolicy mpol, const ExceptionResponse notify) :
    MdlMolDataItem()
{
  // Find an item name in the header line
  const char* line_ptr = tf.getLinePointer(line_number);
  const int lnlength = tf.getLineLength(line_number);
  bool on_item_name = false;
  bool on_regno = false;
  bool name_read = false;
  bool regno_read = false;
  for (int i = 1; i < lnlength; i++) {
    if (line_ptr[i] == '<') {
      if (name_read) {
        switch (notify) {
        case ExceptionResponse::DIE:
        case ExceptionResponse::WARN:
          rtWarn("Two item names, given between < >, were found in a data item "
                 "of SD file " + tf.getFileName() + " (line " + std::to_string(line_number + 1) +
                 ").  Only the first such name will be printed to subsequent output files.",
                 "MdlMolDataItem");
          break;
        case ExceptionResponse::SILENT:
          break;
        }          
      }
      if (on_item_name) {
        rtErr("The reserved character '<' appears twice on line " + std::to_string(line_number) +
              " of file " + getBaseName(tf.getFileName()) + ", compound " + title + ".",
              "MdlMolDataItem");
      }
      on_item_name = true;
      use_item_name = true;
    }
    else if (line_ptr[i] == '>') {
      if (on_item_name == false) {
        rtErr("The reserved character '>' appears before its counterpart '<' on line " +
              std::to_string(line_number) + " of file " + getBaseName(tf.getFileName()) + ", "
              "failing to define a data item name in compound " + title + ".", "MdlMolDataItem");
      }
      on_item_name = false;
    }
    else if (line_ptr[i] == '(') {
      if (regno_read) {
        switch (notify) {
        case ExceptionResponse::DIE:
        case ExceptionResponse::WARN:
          rtWarn("Two external registry numbers, given between ( ), were found in a data item "
                 "of SD file " + tf.getFileName() + " (line " + std::to_string(line_number + 1) +
                 ").  Only the first such registry number will be printed to subsequent output "
                 "files.", "MdlMolDataItem");
          break;
        case ExceptionResponse::SILENT:
          break;
        }          
      }
      if (on_regno) {
        rtErr("The reserved character '(' appears twice on line " + std::to_string(line_number) +
              " of file " + getBaseName(tf.getFileName()) + ", in compound " + title + ".",
              "MdlMolDataItem");
      }
      on_regno = true;
      use_external_regno = true;
    }
    else if (line_ptr[i] == ')') {
      if (on_regno == false) {
        rtErr("The reserved character ')' appears before its counterpart '(' on line " +
              std::to_string(line_number) + " of file " + getBaseName(tf.getFileName()) + ", "
              "failing to define an external registry number in compound " + title + ".",
              "MdlMolDataItem");
      }
      on_regno = false;
      regno_read = true;
    }
    else if (on_item_name == false && on_regno == false) {
      if (line_ptr[i] == 'D' && i < lnlength - 2 && line_ptr[i + 1] == 'T') {
        maccs_ii_number = 0;
        int j = i + 2;
        while (j < lnlength && line_ptr[j] >= '0' && line_ptr[j] <= '9') {
          maccs_ii_number *= 10;
          maccs_ii_number += static_cast<int>(line_ptr[j]) - static_cast<int>('0');
          j++;
        }
        i = j - 1;
        use_maccs_ii_number = true;
      }
      else if (line_ptr[i] == 'F' && i <= lnlength - 13 &&
               strncmpCased(std::string("FROM ARCHIVES"), &line_ptr[i])) {
        note_archives = true;
        i += 12;
      }
      else if (line_ptr[i] >= '0' && line_ptr[i] <= '9') {
        internal_regno = 0;
        int j = i;
        while (j < lnlength && line_ptr[j] >= '0' && line_ptr[j] <= '9') {
          internal_regno *= 10;
          internal_regno += static_cast<int>(line_ptr[j]) - static_cast<int>('0');
          j++;
        }
        i = j - 1;
        use_internal_regno = true;
      }
    }
    else {
      if (on_item_name) {
        item_name += line_ptr[i];
      }
      else if (on_regno) {
        external_regno += line_ptr[i];
      }
    }
  }

  // The data item is native to the file.
  kind = MdlMolDataItemKind::NATIVE;

  // Validate the header line information
  validateItemName(mpol, notify);
  
  // Read the data lines
  int tmp_advance = line_number + 1;
  bool search_on = true;
  const int actual_compound_line_end = (compound_line_end == -1) ? tf.getLineCount() :
                                                                   compound_line_end;
  while (tmp_advance < actual_compound_line_end && search_on) {
    const int dl_length = tf.getLineLength(tmp_advance);
    if (dl_length == 0) {
      search_on = false;
    }
    else {
      search_on = false;
      const char* dl_ptr = tf.getLinePointer(tmp_advance);
      for (int i = 0; i < dl_length; i++) {
        search_on = (search_on || dl_ptr[i] != ' ');
      }
    }
    if (search_on) {
      tmp_advance++;
    }
  }
  tmp_advance--;
  body.reserve(tmp_advance - line_number);
  *line_advance = tmp_advance;
  for (int pos = line_number + 1; pos <= tmp_advance; pos++) {
    body.push_back(tf.extractString(pos));
  }
}

//-------------------------------------------------------------------------------------------------
MdlMolDataItem::MdlMolDataItem(const MdlMolDataRequest &ask,
                               const std::vector<std::string> &body_in,
                               const ModificationPolicy mpol, const ExceptionResponse notify) :
    MdlMolDataItem(ask.getTitle(), ask.getExternalRegistryNumber(), -1, ask.getMaccsFieldNumber(),
                   getDataItemHeaderCode(ask), body_in, mpol, notify)
{
  // The request may incorporate one or more values that are not yet known, e.g. details of the
  // energy evaluation.  Place indications in the data item so that a post-analysis can add the
  // details later.
  switch (ask.getKind()) {
  case DataRequestKind::STATE_VARIABLE:
    kind = MdlMolDataItemKind::STATE_VARIABLE;
    tracked_state = ask.getEnergyComponent();
    break;
  case DataRequestKind::ATOM_INFLUENCES:
    kind = MdlMolDataItemKind::ATOM_INFLUENCES;
    break;
  case DataRequestKind::TOPOLOGY_PARAMETER:
    kind = MdlMolDataItemKind::TOPOLOGY_PARAMETER;
    break;
  case DataRequestKind::STRING:
    kind = MdlMolDataItemKind::STRING;
    break;
  case DataRequestKind::ALL_KINDS:
    rtErr("A nonsensical request type of " + getEnumerationName(ask.getKind()) + " was conveyed.",
          "MdlMolDataItem");
  }
}

//-------------------------------------------------------------------------------------------------
bool MdlMolDataItem::placeInternalRegnoInHeader() const {
  return use_internal_regno;
}

//-------------------------------------------------------------------------------------------------
bool MdlMolDataItem::placeExternalRegnoInHeader() const {
  return use_external_regno;
}

//-------------------------------------------------------------------------------------------------
bool MdlMolDataItem::placeTitleInHeader() const {
  return use_item_name;
}

//-------------------------------------------------------------------------------------------------
bool MdlMolDataItem::placeMaccsIIFieldInHeader() const {
  return use_maccs_ii_number;
}

//-------------------------------------------------------------------------------------------------
bool MdlMolDataItem::noteArchivesInHeader() const {
  return note_archives;
}

//-------------------------------------------------------------------------------------------------
const std::string& MdlMolDataItem::getItemName() const {
  return item_name;
}

//-------------------------------------------------------------------------------------------------
const std::string& MdlMolDataItem::getOutputItemName() const {
  return output_item_name;
}

//-------------------------------------------------------------------------------------------------
const std::string& MdlMolDataItem::getExternalRegistryNumber() const {
  return external_regno;
}

//-------------------------------------------------------------------------------------------------
int MdlMolDataItem::getInternalRegistryNumber() const {
  return internal_regno;
}

//-------------------------------------------------------------------------------------------------
int MdlMolDataItem::getMaccsFieldNumber() const {
  return maccs_ii_number;
}

//-------------------------------------------------------------------------------------------------
StateVariable MdlMolDataItem::getTrackedState() const {
  return tracked_state;
}

//-------------------------------------------------------------------------------------------------
int MdlMolDataItem::getDataLineCount() const {
  return body.size();
}

//-------------------------------------------------------------------------------------------------
const std::string& MdlMolDataItem::getDataLine(const int line_index) const {
  if (line_index < 0 || line_index >= static_cast<int>(body.size())) {
    rtErr("Request for data line " + std::to_string(line_index) + " is invalid for a data item "
          "with " + std::to_string(body.size()) + " lines.", "MdlMolDataItem", "getDataLine");
  }
  return body[line_index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<std::string>& MdlMolDataItem::getBody() const {
  return body;
}

//-------------------------------------------------------------------------------------------------
MdlMolDataItemKind MdlMolDataItem::getKind() const {
  return kind;
}

//-------------------------------------------------------------------------------------------------
std::string MdlMolDataItem::parseString(const int element_number, const int line_number) const {
  if (line_number >= static_cast<int>(body.size())) {
    const std::string identifier = (item_name.size() > 0LLU) ? "Item name : <" + item_name + ">" :
                                                               std::string("");
    rtErr("Line index " + std::to_string(line_number) + " is inaccessible in a data item with " +
          std::to_string(body.size()) + " data lines.  " + identifier + ".", "MdlMolDataItem",
          "parseString");
  }
  const std::vector<std::string> line_components = separateText(body[line_number]);
  if (element_number < static_cast<int>(line_components.size()) || element_number < 0) {
    const std::string identifier = (item_name.size() > 0LLU) ? "Item name : <" + item_name + ">" :
                                                               std::string("");
    rtErr("Line index " + std::to_string(line_number) + " has " +
          std::to_string(line_components.size()) + " components.  Index " +
          std::to_string(element_number) + " is invalid.", "MdlMolDataItem", "parseString");
  }
  return line_components[element_number];
}

//-------------------------------------------------------------------------------------------------
std::string MdlMolDataItem::parseString(const int start_pos, const int length,
                                        const int line_number) const {
  if (line_number >= static_cast<int>(body.size())) {
    const std::string identifier = (item_name.size() > 0LLU) ? "Item name : <" + item_name + ">" :
                                                               std::string("");
    rtErr("Line index " + std::to_string(line_number) + " is inaccessible in a data item with " +
          std::to_string(body.size()) + " data lines.  " + identifier + ".", "MdlMolDataItem",
          "parseString");
  }
  if (start_pos < 0) {
    rtErr("Starting position " + std::to_string(start_pos) + " is invalid.", "MdlMolDataItem",
          "parseString");
  }
  const int lnlength = body[line_number].size();
  const int actual_length = (length < 0) ? lnlength - start_pos : length;
  if (start_pos + length > lnlength) {
    rtErr("Starting position " + std::to_string(start_pos) + " and length " +
          std::to_string(actual_length) + " combine to make an invalid read of a string with "
          "length " + std::to_string(body[line_number].size()) + " characters.", "MdlMolDataItem",
          "parseString");
  }
  return body[line_number].substr(start_pos, actual_length);
}

//-------------------------------------------------------------------------------------------------
llint MdlMolDataItem::parseInteger(const int element_number, const int line_number) const {
  const std::string proto = parseString(element_number, line_number);
  if (verifyContents(proto, NumberFormat::LONG_LONG_INTEGER)) {
    return stoll(proto);
  }
  else {
    rtErr("Element " + std::to_string(element_number) + " of line " + std::to_string(line_number) +
          ", " + proto + ", could not be parsed as an integer.", "MdlMolDataItem", "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
llint MdlMolDataItem::parseInteger(const int start_pos, const int length,
                                   const int line_number) const {
  const std::string proto = parseString(start_pos, length, line_number);
  if (verifyContents(proto, NumberFormat::LONG_LONG_INTEGER)) {
    return stoll(proto);
  }
  else {
    rtErr("Column-formatted characters \"" + proto + "\" at position " +
          std::to_string(start_pos) + " of line " + std::to_string(line_number) + " could not be "
          "parsed as an integer.", "MdlMolDataItem", "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
ullint MdlMolDataItem::parseUnsigned(const int element_number, const int line_number) const {
  const std::string proto = parseString(element_number, line_number);
  if (verifyContents(proto, NumberFormat::UNSIGNED_LONG_LONG_INTEGER)) {
    return stoull(proto);
  }
  else {
    rtErr("Element " + std::to_string(element_number) + " of line " + std::to_string(line_number) +
          ", " + proto + ", could not be parsed as an unsigned integer.", "MdlMolDataItem",
          "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
ullint MdlMolDataItem::parseUnsigned(const int start_pos, const int length,
                                     const int line_number) const {
  const std::string proto = parseString(start_pos, length, line_number);
  if (verifyContents(proto, NumberFormat::UNSIGNED_LONG_LONG_INTEGER)) {
    return stoll(proto);
  }
  else {
    rtErr("Column-formatted characters \"" + proto + "\" at position " +
          std::to_string(start_pos) + " of line " + std::to_string(line_number) + " could not be "
          "parsed as an unsigned integer.", "MdlMolDataItem", "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double MdlMolDataItem::parseReal(const int element_number, const int line_number) const {
  const std::string proto = parseString(element_number, line_number);
  if (verifyContents(proto, NumberFormat::STANDARD_REAL) ||
      verifyContents(proto, NumberFormat::SCIENTIFIC)) {
    return stod(proto);
  }
  else {
    rtErr("Element " + std::to_string(element_number) + " of line " + std::to_string(line_number) +
          ", " + proto + ", could not be parsed as a real number.", "MdlMolDataItem",
          "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double MdlMolDataItem::parseReal(const int start_pos, const int length,
                                 const int line_number) const {
  const std::string proto = parseString(start_pos, length, line_number);
  if (verifyContents(proto, NumberFormat::STANDARD_REAL) ||
      verifyContents(proto, NumberFormat::SCIENTIFIC)) {
    return stod(proto);
  }
  else {
    rtErr("Column-formatted characters \"" + proto + "\" at position " +
          std::to_string(start_pos) + " of line " + std::to_string(line_number) + " could not be "
          "parsed as a real number.", "MdlMolDataItem", "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
char4 MdlMolDataItem::parseChar4(const int element_number, const int line_number) const {
  const std::string proto = parseString(element_number, line_number);
  if (verifyContents(proto, NumberFormat::CHAR4)) {
    return stringToChar4(proto);
  }
  else {
    rtErr("Element " + std::to_string(element_number) + " of line " + std::to_string(line_number) +
          ", " + proto + ", could not be parsed as a four-character tuple, having " +
          std::to_string(proto.size()) + " characters in all.", "MdlMolDataItem",
          "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
char4 MdlMolDataItem::parseChar4(const int start_pos, const int length,
                                 const int line_number) const {
  const std::string proto = parseString(start_pos, length, line_number);
  if (verifyContents(proto, NumberFormat::CHAR4)) {
    return stringToChar4(proto);
  }
  else {
    rtErr("Column-formatted characters \"" + proto + "\" at position " +
          std::to_string(start_pos) + " of line " + std::to_string(line_number) + " could not be "
          "parsed as a four-character tuple.", "MdlMolDataItem", "parseString");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool MdlMolDataItem::matchItemName(const std::string &item_name_comp,
                                   const std::string &ext_regno_comp,
                                   const int maccs_ii_no_comp) const {
  return ((use_item_name && item_name_comp == item_name) &&
          (use_external_regno && ext_regno_comp == external_regno) &&
          (use_maccs_ii_number && maccs_ii_no_comp == maccs_ii_number));
}

//-------------------------------------------------------------------------------------------------
bool MdlMolDataItem::matchItemName(const std::string &item_name_comp,
                                   const int maccs_ii_no_comp) const {
  return ((use_item_name && item_name_comp == item_name) &&
          (use_maccs_ii_number && maccs_ii_no_comp == maccs_ii_number));
}

//-------------------------------------------------------------------------------------------------
bool MdlMolDataItem::matchRegistryNumber(const std::string &ext_regno_comp,
                                         const int maccs_ii_no_comp) const {
  return ((use_external_regno && ext_regno_comp == external_regno) &&
          (use_maccs_ii_number && maccs_ii_no_comp == maccs_ii_number));
}

//-------------------------------------------------------------------------------------------------
bool MdlMolDataItem::matchMaccsField(const int maccs_ii_no_comp) const {
  return (use_maccs_ii_number && maccs_ii_no_comp == maccs_ii_number);
}

//-------------------------------------------------------------------------------------------------
void MdlMolDataItem::setItemName(const std::string &item_name_in,
                                 const ModificationPolicy mpol, const ExceptionResponse notify) {
  item_name = item_name_in;
  validateItemName(mpol, notify);
}

//-------------------------------------------------------------------------------------------------
void MdlMolDataItem::addDataLine(const std::string &text) {
  body.push_back(text);
}

//-------------------------------------------------------------------------------------------------
void MdlMolDataItem::validateItemName(const ModificationPolicy mpol,
                                      const ExceptionResponse notify) {
  bool problem = false;
  bool modified = false;
  output_item_name = item_name;
  const int nchar_in = item_name.size();
  if (nchar_in > 0) {
    switch (mpol) {
    case ModificationPolicy::MODIFY:
      {
        if (output_item_name[0] == '_' ||
            (output_item_name[0] >= '0' && output_item_name[0] <= '9')) {
          output_item_name = 'x' + output_item_name;
          modified = true;
        }
        const int nchar_oin = output_item_name.size();
        for (int i = 0; i < nchar_oin; i++) {
          if (output_item_name[i] == '-' || output_item_name[i] == '.' ||
              output_item_name[i] == '=' || output_item_name[i] == '%' ||
              output_item_name[i] == ' ') {
            output_item_name[i] = '_';
            modified = true;
          }
        }
      }
      break;
    case ModificationPolicy::DO_NOT_MODIFY:
      break;
    }
    for (int i = 0; i < nchar_in; i++) {
      problem = (problem || item_name[i] == '>' || item_name[i] == '<');
    }
  }
  if (modified) {
    switch (notify) {
    case ExceptionResponse::DIE:
    case ExceptionResponse::WARN:
      rtWarn("Data item <" + item_name + "> will be interpreted as <" + output_item_name +
             "> in any output to remain in compliance with Biovia SD file formatting "
             "requirements.", "MdlMolDataItem", "validateItemName");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  if (problem) {
    rtErr("Data item name " + item_name + " is invalid.  An item name must begin with an "
          "alphabetical character and thereafter contain alphanumeric characters and "
          "underscores, with no white space.", "MdlMolDataItem", "validateItemName");
  }
}

//-------------------------------------------------------------------------------------------------
uint getDataItemHeaderCode(const bool use_internal_regno, const bool use_external_regno,
                           const bool use_item_name, const bool use_maccs_ii_field,
                           const bool state_from_archives) {
  return static_cast<int>(use_internal_regno) + (static_cast<int>(use_external_regno) * 2) +
         (static_cast<int>(use_item_name) * 4) + (static_cast<int>(use_maccs_ii_field) * 8) +
         (static_cast<int>(state_from_archives) * 16);
}

//-------------------------------------------------------------------------------------------------
uint getDataItemHeaderCode(const MdlMolDataRequest &ask) {

  // It is obligatory to display the item name for requested data items, and "FROM ARCHIVES" is
  // not displayed in such data items.
  uint result = 4U;
  if (ask.placeInternalRegistryInHeader()) {
    result += 1U;
  }
  if (ask.getExternalRegistryNumber().size() > 0LLU) { 
    result += 2U;
  }
  if (ask.placeMaccsFieldInHeader()) {
    result += 8U;
  }
  return result;
}

} // namespace structure
} // namespace stormm
