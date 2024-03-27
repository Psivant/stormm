#include "copyright.h"
#include "Constants/behavior.h"
#include "FileManagement/file_listing.h"
#include "Parsing/ascii_numbers.h"
#include "Parsing/parse.h"
#include "mdlmol_property.h"

namespace stormm {
namespace structure {

using constants::CaseSensitivity;
using constants::ExceptionResponse;
using diskutil::getBaseName;
using parse::NumberFormat;
using parse::readIntegerValue;
using parse::readRealValue;
using parse::strncmpCased;
using parse::verifyContents;
  
//-------------------------------------------------------------------------------------------------
MdlMolProperty::MdlMolProperty(const char4 code_in, const int substrate_in,
                               const int entry_count_in, const int entry_depth_in,
                               const bool exclusions_in, const int substrate_line_pos_in,
                               const int entry_count_line_pos_in,
                               const int entry_read_start_pos_in,
                               const std::vector<int> &entry_data_bounds_in,
                               const std::vector<MolObjPropField> &entry_detail_in,
                               const std::vector<MolObjIndexKind> &entry_adjustment_in,
                               const std::vector<int> &int_data_in,
                               const std::vector<double> &real_data_in,
                               const std::vector<std::string> &str_data_in,
                               const std::vector<std::string> &data_lines_in) :
    code{code_in},
    kind{translateMdlMolPropertyKind(code)},
    substrate{substrate_in}, entry_count{entry_count_in}, exclusions{exclusions_in},
    substrate_line_pos{substrate_line_pos_in}, entry_count_line_pos{entry_count_line_pos_in},
    entry_read_start_pos{entry_read_start_pos_in}, entry_data_bounds{entry_data_bounds_in},
    entry_depth{entry_depth_in}, entry_detail{entry_detail_in},
    entry_adjustment{entry_adjustment_in}, int_data{int_data_in}, real_data{real_data_in},
    str_data{str_data_in}, data_lines{data_lines_in}
{}

//-------------------------------------------------------------------------------------------------
MdlMolProperty::MdlMolProperty(const TextFile &tf, const int line_number, int *line_advance,
                               const std::string &title) :
    MdlMolProperty()
{
  const char* line_ptr = tf.getLinePointer(line_number);
  const int lnlength = tf.getLineLength(line_number);
  if (lnlength < 6) {
    rtErr("Line " + std::to_string(line_number) + " of " + getBaseName(tf.getFileName()) +
          " cannot contain MDL MOL property information due to its length being only " +
          std::to_string(lnlength) + ".");
  }

  // Read the initial letter and three-letter code, then translate the property kind
  code.w = line_ptr[0];
  code.x = line_ptr[3];
  code.y = line_ptr[4];
  code.z = line_ptr[5];
  kind = translateMdlMolPropertyKind(code);

  // Interpret the code according to known options, starting with the maximum number of entries.
  int max_entries = 10000;
  switch (kind) {
  case MdlMolPropertyKind::ATOM_ALIAS:
  case MdlMolPropertyKind::ATOM_VALUE:
  case MdlMolPropertyKind::GROUP_ABBREVIATION:
  case MdlMolPropertyKind::SGROUP_SUBSCRIPT:
  case MdlMolPropertyKind::SGROUP_BOND_VECTOR:
  case MdlMolPropertyKind::SGROUP_FIELD:
  case MdlMolPropertyKind::SGROUP_DISPLAY:
  case MdlMolPropertyKind::SGROUP_DATA:
  case MdlMolPropertyKind::SPATIAL_FEATURE:
  case MdlMolPropertyKind::PHANTOM_ATOM:
  case MdlMolPropertyKind::SGROUP_CLASS:
  case MdlMolPropertyKind::LARGE_REGNO:
    break;
  case MdlMolPropertyKind::ATOM_LIST:
  case MdlMolPropertyKind::SGROUP_EXPANSION:
  case MdlMolPropertyKind::SGROUP_ATOM_LIST:
  case MdlMolPropertyKind::SGROUP_BOND_LIST:
  case MdlMolPropertyKind::MG_PARENT_ATOM_LIST:
    max_entries = 15;
    break;
  case MdlMolPropertyKind::CHARGE:
  case MdlMolPropertyKind::RADICAL:
  case MdlMolPropertyKind::ISOTOPE:
  case MdlMolPropertyKind::RING_BOND_COUNT:
  case MdlMolPropertyKind::SUBSTITUTION_COUNT:
  case MdlMolPropertyKind::UNSATURATED_COUNT:
  case MdlMolPropertyKind::RGROUP_LABEL_LOCATION:
  case MdlMolPropertyKind::SGROUP_TYPE:
  case MdlMolPropertyKind::SGROUP_SUBTYPE:
  case MdlMolPropertyKind::SGROUP_LABELS:
  case MdlMolPropertyKind::SGROUP_CONNECTIVITY:
  case MdlMolPropertyKind::SGROUP_HIERARCHY:
  case MdlMolPropertyKind::SGROUP_COMP_NUMBER:
  case MdlMolPropertyKind::SGROUP_BRACKET_STYLE:
    max_entries = 8;
    break;
  case MdlMolPropertyKind::SGROUP_CORRESPONDENCE:
  case MdlMolPropertyKind::SGROUP_ATTACH_POINT:
    max_entries = 6;
    break;
  case MdlMolPropertyKind::LINK_ATOM:
  case MdlMolPropertyKind::SGROUP_DISPLAY_INFO:
    max_entries = 4;
    break;
  case MdlMolPropertyKind::ATTACHMENT_POINT:
  case MdlMolPropertyKind::ATTACHMENT_ORDER:
    max_entries = 2;
    break;
  case MdlMolPropertyKind::RGROUP_LOGIC:
    max_entries = 1;
    break;
  case MdlMolPropertyKind::SKIP:
  case MdlMolPropertyKind::NONE:
    break;
  }

  // Determine the layout and allocate the necessary space for each property.  Set the line
  // advancement, if appropriate.
  int tmp_advance = 0;
  bool entry_count_unrecognized = false;
  bool substrate_unrecognized = false;
  switch (kind) {
  case MdlMolPropertyKind::ATOM_ALIAS:
  case MdlMolPropertyKind::ATOM_VALUE:
    substrate_unrecognized = readSubstrateIndex(line_ptr, 3);
    entry_count = 1;
    tmp_advance = 1;
    entry_read_start_pos = 3;
    entry_data_bounds = { 0, 3 };
    break;
  case MdlMolPropertyKind::GROUP_ABBREVIATION:
    substrate_unrecognized = readSubstrateIndex(line_ptr, 3);
    entry_count = 1;
    entry_detail = { MolObjPropField::INTEGER };
    entry_adjustment = { MolObjIndexKind::ATOM };
    tmp_advance = 1;
    entry_read_start_pos = 3;
    entry_data_bounds = { 0, 3, 6 };
    break;
  case MdlMolPropertyKind::CHARGE:
  case MdlMolPropertyKind::RADICAL:
  case MdlMolPropertyKind::ISOTOPE:
  case MdlMolPropertyKind::RING_BOND_COUNT:
  case MdlMolPropertyKind::SUBSTITUTION_COUNT:
  case MdlMolPropertyKind::UNSATURATED_COUNT:
  case MdlMolPropertyKind::ATTACHMENT_POINT:
  case MdlMolPropertyKind::RGROUP_LABEL_LOCATION:
  case MdlMolPropertyKind::SGROUP_TYPE:
  case MdlMolPropertyKind::SGROUP_SUBTYPE:
  case MdlMolPropertyKind::SGROUP_LABELS:
  case MdlMolPropertyKind::SGROUP_HIERARCHY:
  case MdlMolPropertyKind::SGROUP_COMP_NUMBER:
  case MdlMolPropertyKind::SGROUP_BRACKET_STYLE:
    entry_count_unrecognized = readEntryCount(line_ptr);
    entry_detail = { MolObjPropField::INTEGER, MolObjPropField::INTEGER };
    entry_adjustment = { MolObjIndexKind::ATOM, MolObjIndexKind::OTHER };
    entry_read_start_pos = 9;
    entry_data_bounds = { 0, 4, 8 };
    break;
  case MdlMolPropertyKind::SGROUP_CONNECTIVITY:
    entry_count_unrecognized = readEntryCount(line_ptr);
    entry_detail = { MolObjPropField::INTEGER, MolObjPropField::CHAR4 };
    entry_adjustment = { MolObjIndexKind::ATOM, MolObjIndexKind::OTHER };
    entry_read_start_pos = 9;
    entry_data_bounds = { 0, 4, 8 };
    break;
  case MdlMolPropertyKind::LINK_ATOM:
  case MdlMolPropertyKind::RGROUP_LOGIC:
    entry_count_unrecognized = readEntryCount(line_ptr);
    entry_detail = std::vector<MolObjPropField>(4, MolObjPropField::INTEGER);
    entry_adjustment = std::vector<MolObjIndexKind>(4, MolObjIndexKind::OTHER);
    entry_read_start_pos = 9;
    entry_data_bounds = { 0, 4, 8, 12, 16 };
    break;
  case MdlMolPropertyKind::ATOM_LIST:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    exclusions = (lnlength >= 15 && line_ptr[14] == 'T');
    entry_detail = std::vector<MolObjPropField>(4, MolObjPropField::CHAR4);
    entry_adjustment = std::vector<MolObjIndexKind>(4, MolObjIndexKind::OTHER);
    entry_read_start_pos = 16;
    entry_data_bounds = { 0, 4 };
    break;
  case MdlMolPropertyKind::ATTACHMENT_ORDER:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_detail = std::vector<MolObjPropField>(4, MolObjPropField::INTEGER);
    entry_adjustment = { MolObjIndexKind::ATOM, MolObjIndexKind::OTHER, MolObjIndexKind::ATOM,
                         MolObjIndexKind::OTHER };
    entry_read_start_pos = 13;
    entry_data_bounds = { 0, 4, 8, 12, 16 };
    break;
  case MdlMolPropertyKind::SGROUP_EXPANSION:
    if (strncmpCased("EXP", tf.extractString(line_number, 7, 3), CaseSensitivity::YES) == false) {
      rtErr("A malformed S-Group expansion entry was encountered on line " +
            std::to_string(line_number) + " of file " + getBaseName(tf.getFileName()) + ".",
            "MdlMolProperty");
    }
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_detail = { MolObjPropField::INTEGER };
    entry_adjustment = { MolObjIndexKind::BOND };
    entry_read_start_pos = 13;
    entry_data_bounds = { 0, 4 };
    break;
  case MdlMolPropertyKind::SGROUP_ATOM_LIST:
  case MdlMolPropertyKind::SGROUP_BOND_LIST:
  case MdlMolPropertyKind::MG_PARENT_ATOM_LIST:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_detail = { MolObjPropField::INTEGER };
    if (kind == MdlMolPropertyKind::SGROUP_BOND_LIST) {
      entry_adjustment = { MolObjIndexKind::BOND };
    }
    else {
      entry_adjustment = { MolObjIndexKind::ATOM };
    }
    entry_read_start_pos = 13;
    entry_data_bounds = { 0, 4 };
    break;
  case MdlMolPropertyKind::SGROUP_SUBSCRIPT:
    entry_count = 1;
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_detail = { MolObjPropField::STRING };
    entry_adjustment = { MolObjIndexKind::OTHER };
    entry_read_start_pos = 10;
    entry_data_bounds = { 0, tf.getLineLength(line_number) - 11 };
    break;
  case MdlMolPropertyKind::SGROUP_CORRESPONDENCE:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_detail = std::vector<MolObjPropField>(3, MolObjPropField::INTEGER);
    entry_adjustment = std::vector<MolObjIndexKind>(3, MolObjIndexKind::BOND);
    entry_read_start_pos = 13;
    entry_data_bounds = { 0, 4, 8, 12 };
    break;
  case MdlMolPropertyKind::SGROUP_DISPLAY_INFO:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_detail = std::vector<MolObjPropField>(4, MolObjPropField::INTEGER);
    entry_adjustment = std::vector<MolObjIndexKind>(4, MolObjIndexKind::BOND);
    entry_read_start_pos = 13;
    entry_data_bounds = { 0, 3, 6, 9, 12 };
    break;
  case MdlMolPropertyKind::SGROUP_BOND_VECTOR:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count = 1;
    entry_detail = std::vector<MolObjPropField>(3, MolObjPropField::INTEGER);
    entry_adjustment = { MolObjIndexKind::BOND, MolObjIndexKind::OTHER, MolObjIndexKind::OTHER };
    entry_read_start_pos = 10;
    entry_data_bounds = { 0, 4, 7, 10 };
    break;
  case MdlMolPropertyKind::SGROUP_FIELD:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count = 1;
    entry_detail = { MolObjPropField::STRING, MolObjPropField::CHAR4, MolObjPropField::STRING,
                     MolObjPropField::CHAR4, MolObjPropField::STRING };
    entry_adjustment = std::vector<MolObjIndexKind>(5, MolObjIndexKind::OTHER);
    entry_read_start_pos = 10;
    entry_data_bounds = { 0, 30, 32, 52, 54, tf.getLineLength(line_number) - 65 };
    break;
  case MdlMolPropertyKind::SGROUP_DISPLAY:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count = 1;
    entry_detail = { MolObjPropField::REAL, MolObjPropField::REAL, MolObjPropField::INTEGER,
                     MolObjPropField::CHAR4, MolObjPropField::CHAR4, MolObjPropField::CHAR4,
                     MolObjPropField::INTEGER, MolObjPropField::INTEGER, MolObjPropField::INTEGER,
                     MolObjPropField::INTEGER, MolObjPropField::CHAR4, MolObjPropField::INTEGER,
                     MolObjPropField::INTEGER };
    entry_adjustment = std::vector<MolObjIndexKind>(13, MolObjIndexKind::OTHER);
    entry_read_start_pos = 10;
    entry_data_bounds = { 0, 10, 20, 24, 25, 26, 27, 29, 33, 36, 39, 41, 43, 45 };

    // Each "M  SDD" entry will be followed by zero or more "M  SCD" entries and a final "M  SED"
    // entry.  Scan ahead to determine the number of subsequent lines to skip.
    tmp_advance = line_number + 1;
    while (tmp_advance < tf.getLineCount() && tf.getLineLength(tmp_advance) >= 6 &&
           strncmpCased(tf.extractString(tmp_advance, 0, 6), "M  SCD")) {
      tmp_advance++;
    }
    if (strncmpCased(tf.extractString(tmp_advance, 0, 6), "M  SED")) {
      tmp_advance++;
    }
    else {
      rtErr("An MDL MOL property S-group display (\"M  SDD\") entry beginning on line " +
            std::to_string(line_number) + " of file " + getBaseName(tf.getFileName()) + " must "
            "be followed by S-group data line (optional \"M  SCD\" lines, and a terminating "
            "\"M  SED\" line).", "MdlMolProperty");
    }
    tmp_advance -= line_number;
    break;
  case MdlMolPropertyKind::SGROUP_DATA:

    // The S-group data will be read as part of the parent S-group display entry.
    break;
  case MdlMolPropertyKind::SPATIAL_FEATURE:

    // The spatial features properties constitute a distinct block and are read by a separate
    // object.
    break;
  case MdlMolPropertyKind::PHANTOM_ATOM:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count = 1;
    entry_detail = { MolObjPropField::REAL, MolObjPropField::REAL, MolObjPropField::REAL,
                     MolObjPropField::CHAR4, MolObjPropField::STRING };
    entry_adjustment = std::vector<MolObjIndexKind>(5, MolObjIndexKind::OTHER);
    entry_read_start_pos = 9;
    entry_data_bounds = { 0, 10, 20, 31, 35, tf.getLineLength(line_number) - 45 };
    break;
  case MdlMolPropertyKind::SGROUP_ATTACH_POINT:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count_unrecognized = readEntryCount(line_ptr, 10);
    entry_detail = { MolObjPropField::INTEGER, MolObjPropField::INTEGER, MolObjPropField::CHAR4 };
    entry_adjustment = { MolObjIndexKind::ATOM, MolObjIndexKind::ATOM, MolObjIndexKind::OTHER };
    entry_read_start_pos = 13;
    entry_data_bounds = { 0, 4, 8, 10 };
    break;
  case MdlMolPropertyKind::SGROUP_CLASS:
    substrate_unrecognized = readSubstrateIndex(line_ptr);
    entry_count = 1;
    entry_detail = { MolObjPropField::STRING };
    entry_adjustment = { MolObjIndexKind::OTHER };
    entry_read_start_pos = 10;
    entry_data_bounds = { 0, tf.getLineLength(line_number) - 11 };
    break;
  case MdlMolPropertyKind::LARGE_REGNO:
    entry_count = 1;
    entry_detail = { MolObjPropField::INTEGER };
    entry_adjustment = { MolObjIndexKind::OTHER };
    entry_read_start_pos = 6;
    entry_data_bounds = { 0, tf.getLineLength(line_number) - 7 };
    break;
  case MdlMolPropertyKind::SKIP:
    entry_count = 0;
    if (verifyContents(line_ptr, 6, 3, NumberFormat::INTEGER)) {
      tmp_advance = readIntegerValue(line_ptr, 6, 3);
    }
    else {
      rtErr("The S  SKP property cannot skip " + tf.extractString(line_number, 6, 3) + " lines.",
            "MdlMolProperty");
    }
    break;
  case MdlMolPropertyKind::NONE:
    break;
  }
  entry_depth = entry_adjustment.size();
  *line_advance = line_number + tmp_advance;
  if (*line_advance >= tf.getLineCount()) {
    rtErr("The advancement due to property " + std::to_string(code.w) + "  " +
          std::to_string(code.x) + std::to_string(code.y) + std::to_string(code.z) +
          " overruns the length of file " + getBaseName(tf.getFileName()) + ".", "MdlMolProperty");
  }
  if (entry_count_unrecognized) {
    rtErr("The entry count on line " + std::to_string(line_number) + " of file " +
          getBaseName(tf.getFileName()) + " was not recognizable as an integer for property \"" +
          tf.extractString(line_number, 0, 6) + "\".", "MdlMolProperty");
  }
  if (substrate_unrecognized) {
    rtErr("The substrate index on line " + std::to_string(line_number) + " of file " +
          getBaseName(tf.getFileName()) + " was not recognizable as an integer for property \"" +
          tf.extractString(line_number, 0, 6) + "\".", "MdlMolProperty");
  }
  if (entry_count < 0 || entry_count > max_entries) {
    rtErr("Property \"" + tf.extractString(line_number, 0, 6) + "\" on line " +
          std::to_string(line_number) + " of file " + getBaseName(tf.getFileName()) + " has an "
          "invalid number of entries (" + std::to_string(entry_count) + ") (max " +
          std::to_string(max_entries) + ").", "MdlMolProperty");
  }

  // Read the line for each type of property.  Advance the auxiliary counter line_advance if there
  // are additional data lines associated with the property, to avoid later trying to read such
  // information as new properties.
  switch (kind) {
  case MdlMolPropertyKind::ATOM_ALIAS:
  case MdlMolPropertyKind::ATOM_VALUE:
  case MdlMolPropertyKind::GROUP_ABBREVIATION:
    data_lines.push_back(tf.extractString(tmp_advance, 0, tf.getLineLength(line_number)));
    break;
  case MdlMolPropertyKind::SGROUP_DISPLAY:
    parseEntries(tf, line_number, entry_read_start_pos, entry_data_bounds);
    for (int i = line_number + 1; i < tmp_advance; i++) {
      data_lines.push_back(tf.extractString(i, 12, std::min(tf.getLineLength(i) - 12, 69)));
    }
  case MdlMolPropertyKind::SKIP:
    for (int i = line_number + 1; i < tmp_advance; i++) {
      data_lines.push_back(tf.extractString(i, 12, std::min(tf.getLineLength(i) - 12, 69)));
    }
    break;
  case MdlMolPropertyKind::CHARGE:
  case MdlMolPropertyKind::RADICAL:
  case MdlMolPropertyKind::ISOTOPE:
  case MdlMolPropertyKind::RING_BOND_COUNT:
  case MdlMolPropertyKind::SUBSTITUTION_COUNT:
  case MdlMolPropertyKind::UNSATURATED_COUNT:
  case MdlMolPropertyKind::LINK_ATOM:
  case MdlMolPropertyKind::ATOM_LIST:
  case MdlMolPropertyKind::ATTACHMENT_POINT:
  case MdlMolPropertyKind::ATTACHMENT_ORDER:
  case MdlMolPropertyKind::RGROUP_LABEL_LOCATION:
  case MdlMolPropertyKind::RGROUP_LOGIC:
  case MdlMolPropertyKind::SGROUP_TYPE:
  case MdlMolPropertyKind::SGROUP_SUBTYPE:
  case MdlMolPropertyKind::SGROUP_LABELS:
  case MdlMolPropertyKind::SGROUP_CONNECTIVITY:
  case MdlMolPropertyKind::SGROUP_EXPANSION:
  case MdlMolPropertyKind::SGROUP_ATOM_LIST:
  case MdlMolPropertyKind::SGROUP_BOND_LIST:
  case MdlMolPropertyKind::MG_PARENT_ATOM_LIST:
  case MdlMolPropertyKind::SGROUP_CORRESPONDENCE:
  case MdlMolPropertyKind::SGROUP_DISPLAY_INFO:
  case MdlMolPropertyKind::SGROUP_BOND_VECTOR:
  case MdlMolPropertyKind::SGROUP_FIELD:
  case MdlMolPropertyKind::SGROUP_SUBSCRIPT:
  case MdlMolPropertyKind::SGROUP_HIERARCHY:
  case MdlMolPropertyKind::SGROUP_COMP_NUMBER:
  case MdlMolPropertyKind::SGROUP_BRACKET_STYLE:
  case MdlMolPropertyKind::PHANTOM_ATOM:
  case MdlMolPropertyKind::SGROUP_ATTACH_POINT:
  case MdlMolPropertyKind::SGROUP_CLASS:
  case MdlMolPropertyKind::LARGE_REGNO:
    parseEntries(tf, line_number, entry_read_start_pos, entry_data_bounds);
    break;
  case MdlMolPropertyKind::SGROUP_DATA:
  case MdlMolPropertyKind::SPATIAL_FEATURE:
  case MdlMolPropertyKind::NONE:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
MdlMolPropertyKind MdlMolProperty::getKind() const {
  return kind;
}

//-------------------------------------------------------------------------------------------------
char4 MdlMolProperty::getCode() const {
  return code;
}

//-------------------------------------------------------------------------------------------------
int MdlMolProperty::getSubstrate() const {
  return substrate;
}

//-------------------------------------------------------------------------------------------------
int MdlMolProperty::getPrintedSubstrate() const {
  return substrate + 1;
}

//-------------------------------------------------------------------------------------------------
bool MdlMolProperty::applyToExclusions() const {
  if (kind != MdlMolPropertyKind::ATOM_LIST) {
    rtErr("Only the atom list property (\"M  ALS\") has defined exclusion behavior.  \"" +
          std::to_string(code.w) + "  " + std::to_string(code.x) + std::to_string(code.y) +
          std::to_string(code.z) + "\" does not.", "MdlMolProperty", "applyToExclusions");    
  }
  return exclusions;
}

//-------------------------------------------------------------------------------------------------
char MdlMolProperty::getExclusionCode() const {

  // Invoke the public member function to engage its error-checking behavior.
  return (applyToExclusions()) ? 'T' : 'F';
}

//-------------------------------------------------------------------------------------------------
int MdlMolProperty::getEntryCount() const {
  return entry_count;
}

//-------------------------------------------------------------------------------------------------
int MdlMolProperty::getDataLineCount() const {
  return data_lines.size();
}

//-------------------------------------------------------------------------------------------------
int MdlMolProperty::getIntegerValue(const int entry_index, const int attribute_index) const {
  checkAttributeValidity(entry_index, attribute_index, MolObjPropField::INTEGER);
  return int_data[(entry_depth * entry_index) + attribute_index];
}

//-------------------------------------------------------------------------------------------------
int MdlMolProperty::getPrintedIntegerValue(const int entry_index,
                                           const int attribute_index) const {
  checkAttributeValidity(entry_index, attribute_index, MolObjPropField::INTEGER);
  switch (entry_adjustment[attribute_index]) {
  case MolObjIndexKind::ATOM:
  case MolObjIndexKind::BOND:
    return int_data[(entry_depth * entry_index) + attribute_index] + 1;
  case MolObjIndexKind::OTHER:
    return int_data[(entry_depth * entry_index) + attribute_index];
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double MdlMolProperty::getRealValue(const int entry_index, const int attribute_index) const {
  checkAttributeValidity(entry_index, attribute_index, MolObjPropField::REAL);
  return real_data[int_data[(entry_depth * entry_index) + attribute_index]];
}

//-------------------------------------------------------------------------------------------------
char4 MdlMolProperty::getChar4Value(const int entry_index, const int attribute_index) const {
  checkAttributeValidity(entry_index, attribute_index, MolObjPropField::CHAR4);
  const uint ui_result = int_data[(entry_depth * entry_index) + attribute_index];
  return { static_cast<char>(ui_result & 0xff),
           static_cast<char>((ui_result >>  8) & 0xff),
           static_cast<char>((ui_result >> 16) & 0xff),
           static_cast<char>((ui_result >> 24) & 0xff) };
}

//-------------------------------------------------------------------------------------------------
std::string MdlMolProperty::getStringValue(const int entry_index,
                                           const int attribute_index) const {
  checkAttributeValidity(entry_index, attribute_index, MolObjPropField::STRING);
  return str_data[int_data[(entry_depth * entry_index) + attribute_index]];
}

//-------------------------------------------------------------------------------------------------
const std::string& MdlMolProperty::getDataLine(const int index) const {
  if (index < 0 || index > static_cast<int>(data_lines.size())) {
    rtErr("A property (code \"" + std::to_string(code.w) + "  " + std::to_string(code.x) +
          std::to_string(code.y) + std::to_string(code.z) + "\") with " +
          std::to_string(data_lines.size()) + " cannot produce a data line with index " +
          std::to_string(index) + ".", "MdlMolProperty", "getDataLine");
  }
  return data_lines[index];
}

//-------------------------------------------------------------------------------------------------
std::string MdlMolProperty::getMdlText() const {

  // Initialize an empty string and reserve an amount of space that will cover most properties.
  std::string result(128, ' ');
  char* result_data = result.data();

  // Add the property code
  result_data[0] = code.w;
  if (kind != MdlMolPropertyKind::ATOM_ALIAS && kind != MdlMolPropertyKind::ATOM_VALUE &&
      kind != MdlMolPropertyKind::GROUP_ABBREVIATION) {
    result_data[3] = code.x;
    result_data[4] = code.y;
    result_data[5] = code.z;
  }
  
  // Write the substrate index and the number of entries, each time erasing the string terminating
  // character '\0' placed by snprintf().
  if (substrate_line_pos + 3 > 128 || entry_count_line_pos + 3 > 128) {
    result.resize(std::max(substrate_line_pos, entry_count_line_pos) + 3, ' ');
    result_data = result.data();
  }
  if (substrate_line_pos > 0) {
    snprintf(&result_data[substrate_line_pos], 128 - substrate_line_pos, "%3d",
             getPrintedSubstrate());
    result_data[substrate_line_pos + 3] = ' ';
  }
  if (entry_count_line_pos > 0) {
    snprintf(&result_data[entry_count_line_pos], 128 - entry_count_line_pos, "%3d", entry_count);
    result_data[entry_count_line_pos + 3] = ' ';
  }
  if (kind == MdlMolPropertyKind::ATOM_LIST) {
    result_data[14] = getExclusionCode();
  }
  int pos = entry_read_start_pos;
  for (int i = 0; i < entry_count; i++) {
    for (int j = 0; j < entry_depth; j++) {
      const int span = entry_data_bounds[j + 1] - entry_data_bounds[j];
      if (pos + span > static_cast<int>(result.size())) {
        result.resize(result.size() + 128LLU, ' ');
        result_data = result.data();
      }
      switch (entry_detail[j]) {
      case MolObjPropField::INTEGER:
        snprintf(&result_data[pos], 128 - pos, "%*d", span, getPrintedIntegerValue(i, j));
        break;
      case MolObjPropField::REAL:
        snprintf(&result_data[pos], 128 - pos, "%*.4lf", span, getRealValue(i, j));
        break;
      case MolObjPropField::CHAR4:
        {
          const char4 cvv = getChar4Value(i, j);
          if (span >= 1) result_data[pos    ] = cvv.x;
          if (span >= 2) result_data[pos + 1] = cvv.y;
          if (span >= 3) result_data[pos + 2] = cvv.z;
          if (span >= 4) result_data[pos + 3] = cvv.w;
        }
        break;
      case MolObjPropField::STRING:
        snprintf(&result_data[pos], 128 - pos, "%*.*s", span, span, getStringValue(i, j).c_str());
      }
      pos += span;
    }
  }
  if (pos > static_cast<int>(result.size()) - 2) {
    result.resize(pos + 2, ' ');
    result_data = result.data();
  }
  result_data[pos] = '\n';
  pos++;
  result_data[pos] = '\0';
  result.resize(pos);
  const int n_dl = data_lines.size();
  for (int i = 0; i < n_dl; i++) {
    result.append(data_lines[i]);
    result.push_back('\n');
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void MdlMolProperty::setCode(const char x, const char y, const char z, const char major) {
  code.x = x;
  code.y = y;
  code.z = z;
  code.w = major;
}

//-------------------------------------------------------------------------------------------------
void MdlMolProperty::setCode(const char4 code_in) {
  code = code_in;
}

//-------------------------------------------------------------------------------------------------
void MdlMolProperty::setSubstrate(const int index) {
  substrate = index;
}

//-------------------------------------------------------------------------------------------------
void MdlMolProperty::setEntryFormat(const std::vector<MolObjPropField> &entry_detail_in,
                                    const std::vector <MolObjIndexKind> &entry_adjustment_in) {
  if (entry_depth != 0) {
    rtErr("A property with defined fields cannot be redesigned.", "MdlMolProperty",
          "setEntryFormat");
  }
  if (entry_detail_in.size() != entry_adjustment_in.size()) {
    rtErr("Details and adjustment instructions for each property field must have a one-to-one "
          "correspondence.", "MdlMolProperty", "setEntryFormat");
  }
  entry_detail = entry_detail_in;
  entry_adjustment = entry_adjustment_in;
  entry_depth = entry_detail.size();
  for (int i = 0; i < entry_depth; i++) {
    if (entry_detail[i] != MolObjPropField::INTEGER &&
        (entry_adjustment[i] == MolObjIndexKind::ATOM ||
         entry_adjustment[i] == MolObjIndexKind::BOND)) {
      rtErr("Atom and bond indices must have integer type.  Property field " +
            std::to_string(i) + " violates convention by combining a " +
            getEnumerationName(entry_detail[i]) + " data type with " +
            getEnumerationName(entry_adjustment[i]) + " index adjustment.", "MdlMolProperty",
            "setEntryFormat");
    }
  }
}

//-------------------------------------------------------------------------------------------------
bool MdlMolProperty::readEntryCount(const char* line_ptr, const int start_pos, const int length) {
  if (verifyContents(line_ptr, start_pos, length, NumberFormat::INTEGER)) {
    entry_count_line_pos = start_pos;
    entry_count = readIntegerValue(line_ptr, start_pos, length);
    return false;
  }
  else {
    return true;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool MdlMolProperty::readSubstrateIndex(const char* line_ptr, const int start_pos,
                                        const int length) {
  if (verifyContents(line_ptr, start_pos, length, NumberFormat::INTEGER)) {
    substrate_line_pos = start_pos;
    substrate = readIntegerValue(line_ptr, start_pos, length) - 1;
    return false;
  }
  else {
    return true;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void MdlMolProperty::parseEntries(const TextFile &tf, const int line_number, const int start_pos,
                                  const std::vector<int> &limits) {

  // There are integers for every piece of information: char4 data gets converted to a bit-packed
  // integer, pieces of string and double-precision data log integers for their index in the
  // corresponding string or double-precision lists.
  int_data.resize(entry_depth * entry_count);
  int n_real   = 0;
  int n_string = 0;
  for (int i = 0; i < entry_depth; i++) {
    n_real   += (entry_detail[i] == MolObjPropField::REAL);
    n_string += (entry_detail[i] == MolObjPropField::STRING);
  }
  real_data.reserve(n_real * entry_count);
  str_data.reserve(n_string * entry_count);
  const char* line_ptr = tf.getLinePointer(line_number);
  const int lnlength = tf.getLineLength(line_number);
  n_real = 0;
  n_string = 0;
  for (int i = 0; i < entry_count; i++) {
    for (int j = 0; j < entry_depth; j++) {
      const int llim = start_pos + (i * limits[entry_depth]) + limits[j];
      int slen = limits[j + 1] - limits[j];
      slen = (entry_detail[j] == MolObjPropField::CHAR4) ? std::min(slen, 4) : slen;
      if (llim + slen > lnlength) {
        rtErr("Reading entry " + std::to_string(i) + " field " + std::to_string(j) +
              " of property \"" + tf.extractString(line_number, 0, 6) + "\" at line " +
              std::to_string(line_number) + " of file " + getBaseName(tf.getFileName()) +
              "would overrun the line length (" + std::to_string(lnlength) + ").",
              "MdlMolProperty", "parseEntries");
      }
      switch (entry_detail[j]) {
      case MolObjPropField::INTEGER:
        if (verifyContents(line_ptr, llim, slen, NumberFormat::INTEGER)) {
          int_data[(i * entry_depth) + j] = readIntegerValue(line_ptr, llim, slen);
          switch (entry_adjustment[j]) {
          case MolObjIndexKind::ATOM:
          case MolObjIndexKind::BOND:
            int_data[(i * entry_depth) + j] -= 1;
            break;
          case MolObjIndexKind::OTHER:
            break;
          }
        }
        else {
          rtErr("Failed to parse an integer from characters " + std::to_string(llim) + " - " +
                std::to_string(llim + slen - 1) + " of line " + std::to_string(line_number) +
                "(an \"" + tf.extractString(line_number, 0, 6) + "\" property) of file " +
                getBaseName(tf.getFileName()) + ".", "MdlMolProperty", "parseEntries");
        }
        break;
      case MolObjPropField::CHAR4:
        if (verifyContents(line_ptr, llim, slen, NumberFormat::CHAR4)) {
          char4 result;
          result.x = (slen > 0) ? line_ptr[llim    ] : ' ';
          result.y = (slen > 1) ? line_ptr[llim + 1] : ' ';
          result.z = (slen > 2) ? line_ptr[llim + 2] : ' ';
          result.w = (slen > 3) ? line_ptr[llim + 3] : ' ';
          uint uiresult = ((static_cast<uint>(result.w) << 24) |
                           (static_cast<uint>(result.z) << 16) |
                           (static_cast<uint>(result.y) <<  8) |
                           (static_cast<uint>(result.x)));
          int_data[(i * entry_depth) + j] = uiresult;
        }
        break;
      case MolObjPropField::REAL:
        if (verifyContents(line_ptr, llim, slen, NumberFormat::STANDARD_REAL) ||
            verifyContents(line_ptr, llim, slen, NumberFormat::SCIENTIFIC)) {
          int_data[(i * entry_depth) + j] = n_real;
          real_data.push_back(readRealValue(line_ptr, llim, slen));
        }
        break;
      case MolObjPropField::STRING:
        int_data[(i * entry_depth) + j] = n_string;
        str_data.push_back(tf.extractString(line_number, llim, slen));
        break;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MdlMolProperty::checkAttributeValidity(const int entry_index, const int attribute_index,
                                            const MolObjPropField expectation) const {
  if (entry_index < 0 || entry_index >= entry_count) {
    rtErr("Entry " + std::to_string(entry_index) + " does not exist for an \"" +
          std::to_string(code.w) + "  " + std::to_string(code.x) + std::to_string(code.y) +
          std::to_string(code.z) + "\" property with " + std::to_string(entry_count) + " entries.",
          "MdlMolProperty", "checkAttributeValidity");
  }
  if (attribute_index < 0 || attribute_index >= entry_depth) {
    rtErr("Attribute index " + std::to_string(attribute_index) + " is invalid for an \"" +
          std::to_string(code.w) + "  " + std::to_string(code.x) + std::to_string(code.y) +
          std::to_string(code.z) + "\" property.", "MdlMolProperty", "checkAttributeValidity");
  }
  if (entry_detail[attribute_index] != expectation) {
    rtErr("Attribute " + std::to_string(attribute_index) + " is of type " +
          getEnumerationName(entry_detail[attribute_index]).c_str() + ", not " +
          getEnumerationName(expectation).c_str() + ".", "MdlMolProperty",
          "checkAttributeValidity");
  }
}

} // namespace structure
} // namespace stormm
