#include <string>
#include <vector>
#include "copyright.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "namelist_element.h"

namespace stormm {
namespace namelist {

using parse::countDelimiters;
using parse::findStringInVector;
using parse::NumberFormat;
using parse::realToString;
using parse::vectorStrtod;
using parse::vectorStrtol;
using parse::verifyNumberFormat;
using parse::strncmpCased;

//-------------------------------------------------------------------------------------------------
NamelistElement::NamelistElement(const std::string &keyword_in, const NamelistType kind_in,
                                 const std::string &default_in, const DefaultIsObligatory obligate,
                                 const InputRepeats rep_in, const std::string &help_in) :
  label{keyword_in},
  kind{kind_in},
  policy{ExceptionResponse::WARN},
  next_entry_index{(obligate == DefaultIsObligatory::YES)},
  entry_count{1},
  int_values{vectorStrtol(std::vector<std::string>(2, default_in), ExceptionResponse::SILENT)},
  real_values{vectorStrtod(std::vector<std::string>(2, default_in), ExceptionResponse::SILENT)},
  string_values{std::vector<std::string>(2, default_in)},
  accept_multiple_values{(rep_in == InputRepeats::YES)},
  help_message{help_in},
  sub_keys{},
  sub_kinds{},
  sub_int_values{},
  sub_real_values{},
  sub_string_values{},
  sub_help_messages{},
  template_size{0},
  template_ints{},
  template_reals{},
  template_strings{},
  default_int_values{int_values[0]},
  default_real_values{real_values[0]},
  default_string_values{string_values[0]},
  establishment{setEstablishment(std::vector<std::string>(1, default_in))[0]},
  imperative{KeyRequirement::REQUIRED},
  template_establishment{},
  template_requirements{},
  instance_establishment{},
  sub_key_found{}
{
  // Check that the keyword is valid
  if (countDelimiters(keyword_in, {' ', '\n', ',', '\'', '\"', '#', '!', '/'}) > 0) {
    rtErr("Keyword \"" + label + "\" separates into multiple words or may contain comments "
          "or quotations.  The keyword must not contain whitespace or carriage returns and avoid "
          "the namelist restricted characters [ ', \", #, !, and / ].", "NamelistElement");
  }

  // Check that this namelist element is NOT a struct.  The data buffers were resized during
  // initialization in this non-STRUCT case.
  switch (kind) {
  case NamelistType::INTEGER:
  case NamelistType::REAL:
  case NamelistType::STRING:
    break;
  case NamelistType::STRUCT:
    rtErr("Specifying a struct to be contained in a namelist requires separate lists of keywords "
          "and input kinds.  See the overloaded constructor documentation.", "NamelistElement");
  }

  // Check other input: obligatory defaults that don't make sense need to be caught.  The checks
  // were also done when converting the default_in string to vectors of ints and doubles in the
  // initializer list, but the exception response was SILENT because it would nearly always fail
  // for one type or the other, and then be replaced by zero.  Knowing the exact type that is
  // needed, the checks can now be performed with more care.
  std::string param_obl;
  switch (obligate) {
  case DefaultIsObligatory::YES:
    if (accept_multiple_values == false) {
      rtErr("Keyword \"" + label + "\" comes with an obligatory default but does not accept "
            "additional values from user input.  This should not be a keyword but rather a "
            "hard-coded quantity", "NamelistElement");
    }
    param_obl = "n obligatory";
    break;
  case DefaultIsObligatory::NO:

    // Keep entry_count set at 0, per its initialization
    break;
  }
  if (obligate == DefaultIsObligatory::YES) {
    if (default_in.size() == 0) {
      rtErr("A" + param_obl + " default value cannot be obtained from a blank string for "
            "keyword \"" + label + "\".", "NamelistElement");
    }
    switch (kind) {
    case NamelistType::INTEGER:
      if (verifyNumberFormat(default_in.c_str(), NumberFormat::INTEGER)) {
        int_values[0] = stol(default_in);
      }
      else {
        rtErr("A" + param_obl + " default INTEGER cannot be obtained from \"" + default_in +
              "\" for keyword \"" + label + "\".", "NamelistElement");
      }
      break;
    case NamelistType::REAL:
      if (verifyNumberFormat(default_in.c_str(), NumberFormat::SCIENTIFIC) ||
          verifyNumberFormat(default_in.c_str(), NumberFormat::STANDARD_REAL)) {
        real_values[0] = stod(default_in);
      }
      else {
        rtErr("A" + param_obl + " default REAL number cannot be obtained from \"" + default_in +
              "\" for keyword \"" + label + "\".", "NamelistElement");
      }
      break;
    case NamelistType::STRING:
      string_values[0] = default_in;
      break;
    case NamelistType::STRUCT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
NamelistElement::NamelistElement(const std::string keyword_in,
                                 const std::vector<std::string> &sub_keys_in,
                                 const std::vector<NamelistType> &sub_kinds_in,
                                 const std::vector<std::string> &default_list,
                                 const DefaultIsObligatory obligate, const InputRepeats rep_in,
                                 const std::string &help_in,
                                 const std::vector<std::string> &sub_help_in,
                                 const std::vector<KeyRequirement> &template_requirements_in) :
  label{keyword_in},
  kind{NamelistType::STRUCT},
  policy{ExceptionResponse::WARN},
  next_entry_index{(obligate == DefaultIsObligatory::YES)},
  entry_count{1},
  int_values{},
  real_values{},
  string_values{},
  accept_multiple_values{(rep_in == InputRepeats::YES)},
  help_message{help_in},
  sub_keys{sub_keys_in},
  sub_kinds{sub_kinds_in},
  sub_int_values{},
  sub_real_values{},
  sub_string_values{},
  sub_help_messages{sub_help_in},
  template_size{static_cast<int>(sub_keys_in.size())},
  template_ints{vectorStrtol(default_list, ExceptionResponse::SILENT)},
  template_reals{vectorStrtod(default_list, ExceptionResponse::SILENT)},
  template_strings{default_list},
  default_int_values{},
  default_real_values{},
  default_string_values{},
  establishment{InputStatus::DEFAULT},
  imperative{KeyRequirement::REQUIRED},
  template_establishment{setEstablishment(default_list)},
  template_requirements{template_requirements_in},
  instance_establishment{static_cast<size_t>(next_entry_index + 1) * sub_keys_in.size(),
                         InputStatus::MISSING},
  sub_key_found{std::vector<bool>(sub_keys_in.size(), false)}
{
  // Check that this namelist element is a struct.
  switch (kind) {
  case NamelistType::INTEGER:
  case NamelistType::REAL:
  case NamelistType::STRING:
    rtErr("Specifying separate lists of keywords and input kinds implies a STRUCT namelist "
          "variable.", "NamelistElement");
  case NamelistType::STRUCT:
    break;
  }

  // Check that there is a member keyword for every member kind.
  if (sub_keys.size() != sub_kinds.size()) {
    rtErr("The STRUCT keywords \"" + label + "\" has " + std::to_string(sub_keys.size()) +
          " sub-keys associated with it and " + std::to_string(sub_kinds.size()) +
          " data types for them.", "NamelistElement");
  }

  // Fill out any remaining STRUCT template requirements--assume that all fields are required.
  for (size_t i = template_requirements.size(); i < sub_kinds.size(); i++) {
    template_requirements.push_back(KeyRequirement::REQUIRED);
  }
  
  // If the right number of help messages for sub-key was not provided, fill in that array with
  // default values.
  for (size_t i = sub_help_messages.size(); i < sub_keys.size(); i++) {
    sub_help_messages.push_back("No description provided");
  }

  // Resize vector storage and set default values.
  sub_int_values.resize(2 * template_size);
  sub_real_values.resize(2 * template_size);
  sub_string_values.resize(2 * template_size);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < template_size; j++) {
      switch (sub_kinds[j]) {
      case NamelistType::INTEGER:
        sub_int_values[(template_size * i) + j] = template_ints[j];
        break;
      case NamelistType::REAL:
        sub_real_values[(template_size * i) + j] = template_reals[j];
        break;
      case NamelistType::STRING:
        sub_string_values[(template_size * i) + j] = template_strings[j];
        break;
      case NamelistType::STRUCT:
        break;
      }
    }
  }

  // Correct the establishment based on whether the entire template has defaults.  If this is
  // not true, then the establishment of the entire STRUCT will be set to MISSING.
  for (int i = 0; i < template_size; i++) {
    if (template_establishment[i] == InputStatus::MISSING) {
      establishment = InputStatus::MISSING;
    }
  }
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistElement::getLabel() const {
  return label;
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistElement::getSubLabel(const size_t index) const {
  if (kind == NamelistType::STRUCT && index < sub_keys.size()) {
    return sub_keys[index];
  }
  else {
    rtErr("Keyword \"" + label + "\" carries \"" + getEnumerationName(kind) + "\" data with " +
          std::to_string(sub_keys.size()) + " elements.  Access to sub-key " +
          std::to_string(index) + " is invalid.", "NamelistElement", "getSubLabel");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
NamelistType NamelistElement::getKind() const {
  return kind;
}

//-------------------------------------------------------------------------------------------------
NamelistType NamelistElement::getKind(const std::string &sub_key_query) const {

  // Filter out bad use cases
  switch (kind) {
  case NamelistType::INTEGER:
  case NamelistType::REAL:
  case NamelistType::STRING:
    rtErr("Keyword \"" + label + "\" is a " + getEnumerationName(kind) + " with none of "
          "its own member variables.", "NamelistElement", "getKind");
  case NamelistType::STRUCT:
    break;
  }

  // Run through the sub-keys, reporting errors if the query is not found
  const int qry_index = findStringInVector(sub_keys, sub_key_query);
  if (qry_index == static_cast<int>(sub_keys.size())) {
    rtErr("STRUCT keyword \"" + label + "\" has no sub-key \"" + sub_key_query + "\".",
          "NamelistElement", "getKind");
  }
  return sub_kinds[qry_index];
}

//-------------------------------------------------------------------------------------------------
void NamelistElement::setPolicy(ExceptionResponse policy_in) {
  policy = policy_in;
}
  
//-------------------------------------------------------------------------------------------------
void NamelistElement::reportNamelistTypeProblem(const std::string &caller,
                                                const std::string &data_type) const {
  std::string nl_typestr = getEnumerationName(kind);
  std::string term_str;
  switch (kind) {
  case NamelistType::INTEGER:
  case NamelistType::REAL:
  case NamelistType::STRING:
    term_str = ".";
    break;
  case NamelistType::STRUCT:
    term_str = ", apart from one of its members.  Access data by including one of the STRUCT's ";
    term_str += "sub-keys.";
    break;
  }
  rtErr("The keyword \"" + label + "\" is of type " + nl_typestr + " and therefore has no "
        "meaningful " + data_type + " value" + term_str, caller.c_str());
}

//-------------------------------------------------------------------------------------------------
int NamelistElement::getEntryCount() const {
  return entry_count;
}

//-------------------------------------------------------------------------------------------------
int NamelistElement::getIntValue(const std::string &sub_key_query, const int index) const {

  // Check the entry count
  if (index > entry_count) {
    rtErr("Request for INTEGER data in entry " + std::to_string(index) + " out of " +
          std::to_string(entry_count) + " actual values assigned to keyword \"" + label + "\".",
          "getIntValue");
  }

  // If there is a STRUCT member to access, check for INTEGER type and return as appropriate.
  // Otherwise, assume that this is a request for an INTEGER namelist element.
  if (sub_key_query.size() > 0) {
    if (kind != NamelistType::STRUCT) {
      rtErr("Request for INTEGER data assigned to sub-key \"" + sub_key_query + "\" in non-STRUCT "
            "keyword \"" + label + "\".", "getIntValue");
    }
    const size_t member_index = validateSubKey(sub_key_query, NamelistType::INTEGER,
                                               "getIntValue");
    return sub_int_values[(template_size * index) + member_index];
  }
  else {
    switch (kind) {
    case NamelistType::INTEGER:
      return int_values[index];
    case NamelistType::REAL:
    case NamelistType::STRING:
    case NamelistType::STRUCT:
      reportNamelistTypeProblem("getIntValue", "int");
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int NamelistElement::getIntValue(const int index) const {
  return getIntValue(std::string(""), index);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> NamelistElement::getIntValue(const std::string &sub_key_query) const {
  std::vector<int> result(entry_count, 0);

  // Pluck values one by one in order to engage the checks on a keyword's associated type.  This
  // is not efficient but it produces the clearest code and is not a significant computation.
  for (int i = 0; i < entry_count; i++) {
    result[i] = getIntValue(sub_key_query, i);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double NamelistElement::getRealValue(const std::string &sub_key_query, const int index) const {

  // Check the entry count
  if (index > entry_count) {
    rtErr("Request for REAL data in entry " + std::to_string(index) + " out of " +
          std::to_string(entry_count) + " actual values in keyword \"" + label + "\".",
          "getRealValue");
  }

  // If there is a STRUCT member to access, check for REAL type and return as appropriate.
  // Otherwise, assume that this is a request for an REAL namelist element.
  if (sub_key_query.size() > 0) {
    if (kind != NamelistType::STRUCT) {
      rtErr("Request for data associated with REAL sub-key \"" + sub_key_query +
            "\" in non-STRUCT keyword \"" + label + "\".", "getRealValue");
    }
    const size_t member_index = validateSubKey(sub_key_query, NamelistType::REAL, "getRealValue");
    return sub_real_values[(template_size * index) + member_index];
  }
  else {
    switch (kind) {
    case NamelistType::REAL:
      return real_values[index];
    case NamelistType::INTEGER:
    case NamelistType::STRING:
    case NamelistType::STRUCT:
      reportNamelistTypeProblem("getRealValue", "real");
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double NamelistElement::getRealValue(const int index) const {
  return getRealValue(std::string(""), index);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> NamelistElement::getRealValue(const std::string &sub_key_query) const {
  std::vector<double> result(entry_count, 0.0);

  // Pluck values one by one in order to engage the checks on a keyword's associated type.  This
  // is not efficient but it produces the clearest code and is not a significant computation.
  for (int i = 0; i < entry_count; i++) {
    result[i] = getRealValue(sub_key_query, i);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistElement::getStringValue(const std::string &sub_key_query,
                                                   const int index) const {

  // Check the entry count
  if (index > entry_count) {
    rtErr("Request for STRING data in entry " + std::to_string(index) + " out of " +
          std::to_string(entry_count) + " actual values associated with keyword \"" + label +
          "\".", "getStringValue");
  }

  // If there is a STRUCT member to access, check for INTEGER type and return as appropriate.
  // Otherwise, assume that this is a request for an INTEGER namelist element.
  if (sub_key_query.size() > 0) {
    if (kind != NamelistType::STRUCT) {
      rtErr("Request for data associated with STRING sub-key \"" + sub_key_query + "\" in "
            "non-STRUCT keyword \"" + label + "\".", "getStringValue");
    }
    const size_t member_index = validateSubKey(sub_key_query, NamelistType::STRING,
                                               "getStringValue");
    return sub_string_values[(template_size * index) + member_index];
  }
  else {
    switch (kind) {
    case NamelistType::STRING:
      return string_values[index];
    case NamelistType::REAL:
    case NamelistType::INTEGER:
    case NamelistType::STRUCT:
      reportNamelistTypeProblem("getStringValue", "string");
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistElement::getStringValue(const int index) const {
  return getStringValue(std::string(""), index);
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> NamelistElement::getStringValue(const std::string &sub_key_query) const {
  std::vector<std::string> result(entry_count, std::string(""));

  // Pluck values one by one in order to engage the checks on a keyword's associated type.  This
  // is not efficient but it produces the clearest code and is not a significant computation.
  for (int i = 0; i < entry_count; i++) {
    result[i] = getStringValue(sub_key_query, i);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
InputStatus NamelistElement::getEstablishment(const std::string &sub_key_query,
                                              const int repeat_no) const {
  if (sub_key_query.size() > 0) {
    
    // Need to use the sub-key names for the search, as the sub-key type is not known
    const size_t member_index = findStringInVector(sub_keys, sub_key_query);
    if (member_index < sub_keys.size()) {
      if (repeat_no < entry_count) {
        return instance_establishment[(repeat_no * template_size) + member_index];
      }
      else {
        rtErr("The sub-key \"" + sub_key_query + "\" has " + std::to_string(entry_count) +
              " entries, not enough to serve a request for entry index " +
              std::to_string(repeat_no) + ".", "NamelistElement", "getEstablishment");
      }
    }
    else {
      rtErr("No sub-key named \"" + sub_key_query + "\" is associated with STRUCT keyword \"" +
            label + "\".", "NamelistElement", "getEstablishment");
    }
  }
  return establishment;
}

//-------------------------------------------------------------------------------------------------
KeyRequirement NamelistElement::getImperative() const {
  return imperative;
}

//-------------------------------------------------------------------------------------------------
KeyRequirement NamelistElement::getImperative(const std::string &sub_key_query) const {
  if (kind != NamelistType::STRUCT) {
    rtErr("Request for data associated with STRING sub-key \"" + sub_key_query + "\" in "
          "non-STRUCT keyword \"" + label + "\".", "getImperative");
  }
  const size_t member_index = validateSubKey(sub_key_query, NamelistType::STRING, "getImperative");
  return template_requirements[member_index];
}

//-------------------------------------------------------------------------------------------------
bool NamelistElement::getRepeatableValueState() const {
  return accept_multiple_values;
}

//-------------------------------------------------------------------------------------------------
int NamelistElement::getTemplateSize() const {
  return template_size;
}

//-------------------------------------------------------------------------------------------------
void NamelistElement::badInputResponse(const std::string &errmsg, const char* caller) {
  switch (policy) {
  case ExceptionResponse::DIE:
    rtErr(errmsg, "NamelistElement", caller);
  case ExceptionResponse::WARN:
    rtWarn(errmsg, "NamelistElement", caller);
    break;
  case ExceptionResponse::SILENT:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistElement::addDefaultValue(const std::string &next_default) {

  // Check that user-specified values have not already been entered
  if (next_entry_index > 0) {
    rtErr("Additional default values may not be specified after user input has already been "
          "accepted in the \"" + label + "\" keyword.  This indicates a problem with the program, "
          "not the usage.", "NamelistElement", "addDefaultValue");
  }
  if (accept_multiple_values == false) {
    rtErr("The keyword \"" + label + "\" does not accept multiple values, and thus cannot accept "
          "multiple default values.", "NamelistElement", "addDefaultValue");
  }
  switch (kind) {
  case NamelistType::INTEGER:
    int_values.resize(entry_count + 1);
    if (verifyContents(next_default, NumberFormat::INTEGER)) {      
      int_values[entry_count] = stoi(next_default);
      default_int_values.push_back(int_values[entry_count]);
      entry_count++;
    }
    else {
      rtErr("Invalid default input \"" + next_default + "\" to " + getEnumerationName(kind) +
            " keyword \"" + label + "\".", "NamelistElement", "addDefaultValue");
    }
    break;
  case NamelistType::REAL:
    real_values.resize(entry_count + 1);
    if (verifyContents(next_default, NumberFormat::STANDARD_REAL) ||
        verifyContents(next_default, NumberFormat::SCIENTIFIC)) {
      real_values[entry_count] = stod(next_default);
      default_real_values.push_back(real_values[entry_count]);
      entry_count++;
    }
    else {
      rtErr("Invalid default input \"" + next_default + "\" to " + getEnumerationName(kind) +
            " keyword \"" + label + "\".", "NamelistElement", "addDefaultValue");
    }
    break;
  case NamelistType::STRING:
    string_values.resize(entry_count + 1);
    string_values[entry_count] = next_default;
    default_string_values.push_back(string_values[entry_count]);
    entry_count++;
  case NamelistType::STRUCT:

    // STRUCT-type keywords' default values are stored in the template_(ints / reals / strings)
    // arrays.
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistElement::setIntValue(const int value) {
  if (kind != NamelistType::INTEGER) {
    reportNamelistTypeProblem("setIntValue", "int");
  }
  if (next_entry_index > 0 && accept_multiple_values == false) {
    const std::string errmsg = "An extra entry was found for keyword \"" + label + "\".  " +
                               std::to_string(value) + " will replace " +
                               std::to_string(int_values[0]) + ".";
    int_values[0] = value;
    entry_count = 1;
    next_entry_index = 1;
    badInputResponse(errmsg, "setIntValue");
  }
  else {
    int_values[next_entry_index] = value;
    resizeBuffer();
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistElement::setIntValue(const std::string &sub_key_query, const int value) {
  if (kind != NamelistType::STRUCT) {
    reportNamelistTypeProblem("setIntValue", "int");
  }
  const int member_index = validateSubKey(sub_key_query, NamelistType::INTEGER, "setIntValue");
  if ((next_entry_index > 0 && accept_multiple_values == false) || sub_key_found[member_index]) {
    const std::string errmsg = "An extra entry was found for keyword \"" + label +
                               "\", sub-key \"" + sub_key_query + "\".  " +
                               std::to_string(value) + " will replace " +
                               std::to_string(sub_int_values[member_index]) + ".";
    sub_int_values[member_index] = value;
    badInputResponse(errmsg, "setIntValue");
  }
  else {
    sub_int_values[(next_entry_index * template_size) +  member_index] = value;
    sub_key_found[member_index] = true;
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistElement::setRealValue(const double value) {
  if (kind != NamelistType::REAL) {
    reportNamelistTypeProblem("setRealValue", "double");
  }
  if (next_entry_index > 0 && accept_multiple_values == false) {
    const std::string errmsg = "An extra entry was found for keyword \"" + label + "\".  " +
                               realToString(value) + " will replace " +
                               realToString(real_values[0]) + ".";
    real_values[0] = value;
    entry_count = 1;
    next_entry_index = 1;
    badInputResponse(errmsg, "setRealValue");
  }
  else {
    real_values[next_entry_index] = value;
    resizeBuffer();
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistElement::setRealValue(const std::string &sub_key_query, const double value) {
  if (kind != NamelistType::STRUCT) {
    reportNamelistTypeProblem("setRealValue", "double");
  }
  const int member_index = validateSubKey(sub_key_query, NamelistType::REAL, "setRealValue");
  if ((next_entry_index > 0 && accept_multiple_values == false) || sub_key_found[member_index]) {
    const std::string errmsg = "An extra entry was found for keyword \"" + label +
                               "\", sub-key \"" + sub_key_query + "\".  " + realToString(value) +
                               " will replace " + realToString(sub_real_values[member_index]) +
                               ".";
    sub_real_values[member_index] = value;
    badInputResponse(errmsg, "setRealValue");
  }
  else {
    sub_real_values[(next_entry_index * template_size) + member_index] = value;
    sub_key_found[member_index] = true;
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistElement::setStringValue(const std::string &value) {
  if (kind != NamelistType::STRING) {
    reportNamelistTypeProblem("setStringValue", "string");
  }
  if (next_entry_index > 0 && accept_multiple_values == false) {
    const std::string errmsg = "An extra entry was found for keyword \"" + label + "\".  \"" +
                               value + "\" will replace \"" + string_values[0] + "\".";
    string_values[0] = value;
    entry_count = 1;
    next_entry_index = 1;
    badInputResponse(errmsg, "setStringValue");
  }
  else {
    string_values[next_entry_index] = value;
    resizeBuffer();
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistElement::setStringValue(const std::string &sub_key_query, const std::string &value) {
  if (kind != NamelistType::STRUCT) {
    reportNamelistTypeProblem("setStringValue", "string");
  }
  const int member_index = validateSubKey(sub_key_query, NamelistType::STRING, "setStringValue");
  if ((next_entry_index > 0 && accept_multiple_values == false) || sub_key_found[member_index]) {
    const std::string errmsg = "An extra entry was found for keyword \"" + label +
                               "\", sub-key \"" + sub_key_query + "\".  \"" + value +
                               "\" will replace \"" + sub_string_values[member_index] + "\".";
    sub_string_values[member_index] = value;
    badInputResponse(errmsg, "setStringValue");
  }
  else {
    sub_string_values[(next_entry_index * template_size) +  member_index] = value;
    sub_key_found[member_index] = true;
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistElement::setImperative(const KeyRequirement imperative_in) {
  imperative = imperative_in;
}

//-------------------------------------------------------------------------------------------------
void NamelistElement::setImperative(const std::string &sub_key_query,
                                    const KeyRequirement imperative_in) {
  const int member_index = validateSubKey(sub_key_query, NamelistType::REAL, "setRealValue");
  template_requirements[member_index] = imperative_in;
}

//-------------------------------------------------------------------------------------------------
bool NamelistElement::defaultIsMissing(const std::string &default_input) const {

  // A zero-length string means missing data, as do two reserved words
  return (default_input.size() == 0 || strncmpCased(default_input, "MISSING") ||
          strncmpCased(default_input, "BLANK"));
}

//-------------------------------------------------------------------------------------------------
std::vector<InputStatus>
NamelistElement::setEstablishment(const std::vector<std::string> &list_of_defaults) {
  std::vector<InputStatus> result;
  const size_t n_defaults = list_of_defaults.size();
  switch (kind) {
  case NamelistType::INTEGER:
  case NamelistType::REAL:
  case NamelistType::STRING:
    if (n_defaults != 1) {
      rtErr("\"" + label + "\", a non-STRUCT namelist element, can only have one point of input "
            "establishment.", "NamelistElement", "setEstablishment");
    }
    break;
  case NamelistType::STRUCT:
    if (n_defaults > sub_keys.size()) {
      rtErr("The STRUCT namelist element \"" + label + "\" has " +
            std::to_string(sub_keys.size()) + " member variables but " +
            std::to_string(n_defaults) + " default values.", "NamelistElement",
            "setEstablishment");
    }
    break;
  }
  for (size_t i = 0; i < n_defaults; i++) {

    // All elements but the STRUCT are looking at the list of defaults through their main
    // keyword-level kind, but the STRUCT is looking at it each default through its corresponding
    // entry in a list of sub-kinds
    const NamelistType relevant_kind = (kind == NamelistType::STRUCT) ? sub_kinds[i] : kind;
    switch (relevant_kind) {
    case NamelistType::INTEGER:
      if (verifyNumberFormat(list_of_defaults[i].c_str(), NumberFormat::INTEGER)) {
        result.push_back(InputStatus::DEFAULT);
      }
      else {
        if (defaultIsMissing(list_of_defaults[i])) {
          result.push_back(InputStatus::MISSING);
        }
        else {
          rtErr("A default INTEGER cannot be obtained from " + list_of_defaults[i] +
                " for namelist element \"" + label + "\".", "NamelistElement",
                "setEstablishment");
        }
      }
      break;
    case NamelistType::REAL:
      if (verifyNumberFormat(list_of_defaults[i].c_str(), NumberFormat::STANDARD_REAL) ||
          verifyNumberFormat(list_of_defaults[i].c_str(), NumberFormat::SCIENTIFIC)) {
        result.push_back(InputStatus::DEFAULT);
      }
      else {
        if (defaultIsMissing(list_of_defaults[i])) {
          result.push_back(InputStatus::MISSING);
        }
        else {
          rtErr("A default REAL number cannot be obtained from " + list_of_defaults[i] +
                " for namelist element \"" + label + "\".", "NamelistElement",
                "setEstablishment");
        }
      }
      break;
    case NamelistType::STRING:
      if (defaultIsMissing(list_of_defaults[i])) {
        result.push_back(InputStatus::MISSING);
      }
      else {
        result.push_back(InputStatus::DEFAULT);
      }
      break;
    case NamelistType::STRUCT:

      // There are no STRUCT elements nested within other namelist STRUCTs, and such errors
      // will have been caught already.
      break;
    }
  }
  for (size_t i = n_defaults; i < sub_keys.size(); i++) {
    result.push_back(InputStatus::MISSING);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int NamelistElement::validateSubKey(const std::string &sub_key_query, const NamelistType nmlt,
                                    const std::string &caller) const {
  if (kind != NamelistType::STRUCT) {
    rtErr("A member keyword such as \"" + sub_key_query + "\" can only be specified for a "
          "STRUCT-type element.", "NamelistElement", caller.c_str());
  }
  const int member_index = findStringInVector(sub_keys, sub_key_query);
  if (member_index == static_cast<int>(sub_keys.size())) {
    rtErr("The namelist element \"" + label + "\" contains no sub-element \"" + sub_key_query +
          "\".", "NameListElement", caller.c_str());
  }
  std::string data_type = getEnumerationName(nmlt);
  switch (nmlt) {
  case NamelistType::INTEGER:
  case NamelistType::REAL:
  case NamelistType::STRING:
    break;
  case NamelistType::STRUCT:
    rtErr("A STRUCT namelist element cannot contain other STRUCTs, only INTEGER, REAL, and STRING "
          "variables.", "NamelistElement", "validateSubKey");
  }
  if (sub_kinds[member_index] != nmlt) {
    rtErr("The namelist element \"" + label + "\", sub-element \"" + sub_key_query +
          "\", is not of type " + data_type + ".", "NameListElement", caller.c_str());
  }
  return member_index;
}

//-------------------------------------------------------------------------------------------------
void NamelistElement::resizeBuffer() {
  next_entry_index++;
  entry_count = next_entry_index;
  switch (kind) {
  case NamelistType::INTEGER:
    int_values.resize(next_entry_index + 1);
    break;
  case NamelistType::REAL:
    real_values.resize(next_entry_index + 1);
    break;
  case NamelistType::STRING:
    string_values.resize(next_entry_index + 1);
    break;
  case NamelistType::STRUCT:
    {
      const size_t current_size = template_size * next_entry_index;
      const size_t next_size = template_size * (next_entry_index + 1);
      sub_int_values.resize(next_size);
      sub_real_values.resize(next_size);
      sub_string_values.resize(next_size);
      instance_establishment.resize(next_size);

      // STRUCT sub-keys must be initialized to their default values, as they may not all get set
      // when a new STRUCT is added to the list of entries.  Also reset the read counters for the
      // STRUCT at hand.
      for (int i = 0; i < template_size; i++) {
        sub_int_values[current_size + i] = template_ints[i];
        sub_real_values[current_size + i] = template_reals[i];
        sub_string_values[current_size + i] = template_strings[i];
        instance_establishment[current_size + i] = InputStatus::MISSING;
        sub_key_found[i] = false;
      }
    }
    break;
  }
}

} // namespace namelist
} // namespace stormm
