#include <string>
#include <vector>
#include "copyright.h"
#include "FileManagement/file_util.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/tabulation.h"
#include "Reporting/display.h"
#include "Reporting/error_format.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using errors::RTMessageKind;
using display::horizontalRule;
using display::terminalHorizontalRule;
using parse::addLeadingWhiteSpace;
using parse::findStringInVector;
using parse::JustifyText;
using parse::justifyStrings;
using parse::lowercase;
using parse::minimalRealFormat;
using parse::NumberFormat;
using parse::removeLeadingWhiteSpace;
using parse::strcmpCased;
  
//-------------------------------------------------------------------------------------------------
NamelistEmulator::NamelistEmulator(const std::string &title_in, const CaseSensitivity casing_in,
                                   const ExceptionResponse unknown_keyword_policy,
                                   const std::string &help_in, const bool cli_content_in) :
    cli_content{cli_content_in},
    title{title_in},
    keywords{},
    casing{casing_in},
    policy{unknown_keyword_policy},
    help_message{help_in},
    category_names{},
    categories{}
{}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistEmulator::getTitle() const {
  return title;
}

//-------------------------------------------------------------------------------------------------
bool NamelistEmulator::isCommandLineContent() const {
  return cli_content;
}
  
//-------------------------------------------------------------------------------------------------
int NamelistEmulator::getKeywordCount() const {
  return keywords.size();
}

//-------------------------------------------------------------------------------------------------
CaseSensitivity NamelistEmulator::getCaseSensitivity() const {
  return casing;
}

//-------------------------------------------------------------------------------------------------
ExceptionResponse NamelistEmulator::getPolicy() const {
  return policy;
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistEmulator::getKeyword(const size_t index) const {
  if (index > keywords.size()) {
    rtErr("Namelist \"" + title + "\" has " + std::to_string(keywords.size()) + " keywords, " +
          "not enough to return index " + std::to_string(index) + ".", "getKeyword");
  }
  return keywords[index].getLabel();
}

//-------------------------------------------------------------------------------------------------
NamelistType NamelistEmulator::getKeywordKind(const std::string &keyword_query) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getKeywordKind");
  }
  return keywords[p_index].getKind();
}

//-------------------------------------------------------------------------------------------------
int NamelistEmulator::getKeywordEntries(const std::string &keyword_query) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getKeywordEntries");
  }
  switch (keywords[p_index].kind) {
  case NamelistType::BOOLEAN:
  case NamelistType::INTEGER:
  case NamelistType::REAL:
  case NamelistType::STRING:
    return keywords[p_index].getEntryCount();
  case NamelistType::STRUCT:
    return keywords[p_index].getEntryCount() *
           (keywords[p_index].getEstablishment() != InputStatus::MISSING);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int NamelistEmulator::getSubKeyCount(const std::string &keyword_query) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getSubKeyCount");
  }
  return keywords[p_index].getTemplateSize();
}

//-------------------------------------------------------------------------------------------------
InputStatus NamelistEmulator::getKeywordStatus(const std::string &keyword_query) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getKeywordEntries");
  }
  return keywords[p_index].getEstablishment();
}

//-------------------------------------------------------------------------------------------------
InputStatus NamelistEmulator::getKeywordStatus(const std::string &keyword_query,
                                               const std::string &sub_key,
                                               const int repeat_no) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getKeywordEntries");
  }
  return keywords[p_index].getEstablishment(sub_key, repeat_no);
}

//-------------------------------------------------------------------------------------------------
bool NamelistEmulator::hasKeyword(const std::string &query) const {
  const size_t p_index = findIndexByKeyword(query);
  return (p_index < keywords.size());
}

//-------------------------------------------------------------------------------------------------
bool NamelistEmulator::hasKeyword(const std::string &query, const NamelistType query_kind) const {
  const size_t p_index = findIndexByKeyword(query);
  return (p_index < keywords.size() && keywords[p_index].getKind() == query_kind);
}

//-------------------------------------------------------------------------------------------------
bool NamelistEmulator::getBoolValue(const std::string &keyword_query) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getBoolValue");
  }
  verifyEstablishment(keyword_query, p_index, "getBoolValue");
  return keywords[p_index].getBoolValue();
}

//-------------------------------------------------------------------------------------------------
bool NamelistEmulator::getBoolValue(const std::string &keyword_query, const std::string &sub_key,
                                    const int index) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getBoolValue");
  }
  verifyEstablishment(keyword_query, p_index, "getBoolValue");
  return keywords[p_index].getBoolValue(sub_key, index);
}

//-------------------------------------------------------------------------------------------------
int NamelistEmulator::getIntValue(const std::string &keyword_query, const int index) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getIntValue");
  }
  verifyEstablishment(keyword_query, p_index, "getIntValue");
  return keywords[p_index].getIntValue(index);
}

//-------------------------------------------------------------------------------------------------
int NamelistEmulator::getIntValue(const std::string &keyword_query, const std::string &sub_key,
                                  const int index) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getIntValue");
  }
  verifyEstablishment(keyword_query, p_index, "getIntValue");
  return keywords[p_index].getIntValue(sub_key, index);
}

//-------------------------------------------------------------------------------------------------
double NamelistEmulator::getRealValue(const std::string &keyword_query, const int index) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getRealValue");
  }
  verifyEstablishment(keyword_query, p_index, "getRealValue");
  return keywords[p_index].getRealValue(index);
}

//-------------------------------------------------------------------------------------------------
double NamelistEmulator::getRealValue(const std::string &keyword_query, const std::string &sub_key,
                                      const int index) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getRealValue");
  }
  verifyEstablishment(keyword_query, p_index, "getRealValue");
  return keywords[p_index].getRealValue(sub_key, index);
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistEmulator::getStringValue(const std::string &keyword_query,
                                                    const int index) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getStringValue");
  }
  verifyEstablishment(keyword_query, p_index, "getStringValue");
  return keywords[p_index].getStringValue(index);
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistEmulator::getStringValue(const std::string &keyword_query,
                                                    const std::string &sub_key,
                                                    const int index) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getStringValue");
  }
  verifyEstablishment(keyword_query, p_index, "getStringValue");
  return keywords[p_index].getStringValue(sub_key, index);
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> NamelistEmulator::getAllBoolValues(const std::string &keyword_query,
                                                     const std::string &sub_key) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getAllBoolValues");
  }
  verifyEstablishment(keyword_query, p_index, "getAllBoolValues");
  return keywords[p_index].getBoolValue(sub_key);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> NamelistEmulator::getAllIntValues(const std::string &keyword_query,
                                                   const std::string &sub_key) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getAllIntValues");
  }
  verifyEstablishment(keyword_query, p_index, "getAllIntValues");
  return keywords[p_index].getIntValue(sub_key);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> NamelistEmulator::getAllRealValues(const std::string &keyword_query,
                                                       const std::string &sub_key) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getAllRealValues");
  }
  verifyEstablishment(keyword_query, p_index, "getAllRealValues");
  return keywords[p_index].getRealValue(sub_key);
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> NamelistEmulator::getAllStringValues(const std::string &keyword_query,
                                                              const std::string &sub_key) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getAllStringValues");
  }
  verifyEstablishment(keyword_query, p_index, "getAllStringValues");
  return keywords[p_index].getStringValue(sub_key);
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistEmulator::getHelp() const {
  return help_message;
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistEmulator::getHelp(const std::string &keyword_query) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getHelp");
  }
  return keywords[p_index].help_message;
}

//-------------------------------------------------------------------------------------------------
const std::string& NamelistEmulator::getHelp(const std::string &keyword_query,
                                      const std::string &sub_key) const {
  const size_t p_index = findIndexByKeyword(keyword_query);
  if (p_index >= keywords.size()) {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + keyword_query + "\".",
          "NamelistEmulator", "getHelp");
  }
  if (keywords[p_index].kind != NamelistType::STRUCT) {
    rtErr("The keyword \"" + keyword_query + "\" is a " +
          getEnumerationName(keywords[p_index].kind) + " keyword.  Only STRUCTs can produce user "
          "documentation on sub-keys such as \"" + sub_key + "\".", "NamelistEmulator", "getHelp");
  }
  const size_t member_index = findStringInVector(keywords[p_index].sub_keys, sub_key);
  if (member_index >= keywords[p_index].sub_help_messages.size()) {
    rtErr("No sub-key \"" + sub_key + "\" is present in namelist \"" + title + "\" keyword \"" +
          keyword_query + "\".", "NamelistEmulator", "getHelp");
  }
  return keywords[p_index].sub_help_messages[member_index];
}

//-------------------------------------------------------------------------------------------------
const NamelistEmulator* NamelistEmulator::getSelfPointer() const {
  return this;
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::assignVariable(std::string *var, const std::string &keyword_query,
                                      const int index) const {
  if (getKeywordStatus(keyword_query) != InputStatus::MISSING) {
    *var = getStringValue(keyword_query, index);
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::assignVariable(std::string *var, const std::string &keyword_query,
                                      const std::string &sub_key, const int index) const {
  if (getKeywordStatus(keyword_query, sub_key, index) != InputStatus::MISSING) {
    *var = getStringValue(keyword_query, sub_key, index);
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::setTitle(const std::string &title_in) {
  title = title_in;
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::setCommandLineContent(const bool cli_content_in) {
  cli_content = cli_content_in;
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addKeyword(const std::vector<NamelistElement> &new_keys) {
  const int n_new_key = new_keys.size();
  for (int i = 0; i < n_new_key; i++) {

    // Check that the namelist does not already contain this keyword.  Such would be an outright
    // error, as no developer should be doing this when making an actual program.
    const int n_stored = keywords.size();
    for (int j = 0; j < n_stored; j++) {
      if (new_keys[i].label == keywords[j].label) {
        rtErr("Namelist \"" + title + "\" already has a " +
              getEnumerationName(keywords[j].kind) + " keyword " + keywords[j].label +
              ".  A " + getEnumerationName(new_keys[i].kind) + " keyword of the same name cannot "
              "be added.", "NamelistEmulator", "addKeyword");
      }
    }
    keywords.push_back(new_keys[i]);
    keywords.back().setPolicy(policy);
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addKeywords(const std::vector<NamelistElement> &new_keys) {
  addKeyword(new_keys);
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addKeyword(const NamelistElement &new_key) {
  std::vector<NamelistElement> new_keys(1, new_key);
  addKeyword(new_keys);
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addKeyword(const std::string &keyword_in, const NamelistType kind_in,
                                  const std::string &default_in,
                                  const DefaultIsObligatory obligate, const InputRepeats rep_in,
                                  const std::string &help_in) {
  addKeyword(NamelistElement(keyword_in, kind_in, default_in, obligate, rep_in, help_in));
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addKeyword(const std::string keyword_in,
                                  const std::vector<std::string> &sub_keys_in,
                                  const std::vector<NamelistType> &sub_kinds_in,
                                  const std::vector<std::string> &default_list,
                                  const DefaultIsObligatory obligate_list,
                                  const InputRepeats rep_in, const std::string &help_in,
                                  const std::vector<std::string> &sub_help_in,
                                  const std::vector<KeyRequirement> &template_requirements_in) {
  addKeyword(NamelistElement(keyword_in, sub_keys_in, sub_kinds_in, default_list, obligate_list,
                             rep_in, help_in, sub_help_in, template_requirements_in));
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addKeyword(const NamelistEmulator *other, const std::string &query) {
  if (this->hasKeyword(query)) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Namelist \"" + other->getTitle() + "\" already contains keyword \"" + query + "\".",
            "NamelistEmulator", "addKeyword");
    case ExceptionResponse::WARN:
      rtWarn("Namelist \"" + other->getTitle() + "\" already contains keyword \"" + query + "\".  "
             "The keyword will be neither imported nor modified.", "NamelistEmulator",
             "addKeyword");
      return;
    case ExceptionResponse::SILENT:
      return;
    }
  }
  const int noth_keys = other->getKeywordCount();
  if (noth_keys == 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Namelist \"" + other->getTitle() + "\" has no keywords.", "NamelistEmulator",
            "addKeyword");
    case ExceptionResponse::WARN:
      rtWarn("Namelist \"" + other->getTitle() + "\" has no keywords.", "NamelistEmulator",
             "addKeyword");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  const int oth_idx = other->findIndexByKeyword(query);
  if (oth_idx < noth_keys) {

    // Copy the key from the other namelist and set its error response policy to that of this
    // namelist object.
    keywords.push_back(other->keywords[oth_idx]);
    keywords.back().setPolicy(policy);
  }
  else {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Namelist \"" + other->getTitle() + "\" has no keyword \"" + query + "\".",
            "NamelistEmulator", "addKeyword");
    case ExceptionResponse::WARN:
      rtWarn("Namelist \"" + other->getTitle() + "\" has no keyword \"" + query + "\".",
             "NamelistEmulator", "addKeyword");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::setDefaultValue(const std::string &key,
                                       const std::string &modified_default,
                                       const int default_idx) {
  const size_t param_index = findIndexByKeyword(key);
  if (param_index >= keywords.size()) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
            "setDefaultValue");
    case ExceptionResponse::WARN:
      rtWarn("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
             "setDefaultValue");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  keywords[param_index].setDefaultValue(modified_default, default_idx);
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::setDefaultValue(const std::string &key,
                                       const std::vector<std::string> &modified_defaults,
                                       const std::vector<std::string> &sub_key_specs) {
  const size_t param_index = findIndexByKeyword(key);
  if (param_index >= keywords.size()) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
            "setDefaultValue");
    case ExceptionResponse::WARN:
      rtWarn("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
             "setDefaultValue");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  keywords[param_index].setDefaultValue(modified_defaults, sub_key_specs);
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addDefaultValue(const std::string &key, const std::string &next_default) {
  const size_t param_index = findIndexByKeyword(key);
  if (param_index >= keywords.size()) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
            "addDefaultValue");
    case ExceptionResponse::WARN:
      rtWarn("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
             "addDefaultValue");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  keywords[param_index].addDefaultValue(next_default);
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::activateBool(const std::string &key) {
  
  // The presence of a keyword and no paired value is only valid in the case of a boolean input.
  const size_t param_index = findIndexByKeyword(key);
  if (param_index >= keywords.size()) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
            "activateBool");
    case ExceptionResponse::WARN:
      rtWarn("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
             "activateBool");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return;
  }
  const NamelistType param_type = keywords[param_index].kind;
  bool problem = false;
  switch (param_type) {
  case NamelistType::BOOLEAN:
    keywords[param_index].activateBool();
    keywords[param_index].establishment = InputStatus::USER_SPECIFIED;
    break;
  case NamelistType::INTEGER:
  case NamelistType::REAL:
  case NamelistType::STRING:
  case NamelistType::STRUCT:
    break;
  }
  return;
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::activateBool(const std::string &key, const std::string &sub_key) {

  // The inputs contain a keyword, a sub-key, and a value.  This is only valid for STRUCTs.
  const size_t param_index = findIndexByKeyword(key);
  if (param_index >= keywords.size()) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
            "activateBool");
    case ExceptionResponse::WARN:
      rtWarn("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
             "activateBool");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return;
  }
  const NamelistType param_type = keywords[param_index].kind;
  bool problem = false;
  switch (param_type) {
  case NamelistType::STRUCT:
    {
      // Search within this keyword.  Use the getter function for encapsulated functionality,
      // even though friendship lets this NamelistEmulator go right in.
      NamelistType subtype = keywords[param_index].getKind(sub_key);
      bool problem = false;
      switch (subtype) {
      case NamelistType::STRUCT:
        break;
      case NamelistType::BOOLEAN:
        keywords[param_index].activateBool(sub_key);
        break;
      case NamelistType::INTEGER:
      case NamelistType::REAL:
      case NamelistType::STRING:
        problem = true;
        break;
      }

      // Respond to input errors
      if (problem) {
        switch(policy) {
        case ExceptionResponse::DIE:
          rtErr("In namelist \"" + title + "\", keyword \"" + key + "\", sub-key \"" +
                sub_key + "\" accepts " + getEnumerationName(subtype) + " values.",
                "NamelistEmulator", "activateBool");
        case ExceptionResponse::WARN:
          rtWarn("In namelist \"" + title + "\", keyword \"" + key + "\", sub-key \"" +
                 sub_key + "\" accepts " + getEnumerationName(subtype) + " values.",
                 "NamelistEmulator", "activateBool");
          return;
        case ExceptionResponse::SILENT:
          return;
        }
      }
    }
    break;
  case NamelistType::BOOLEAN:
  case NamelistType::INTEGER:
  case NamelistType::REAL:
  case NamelistType::STRING:
    break;
  }
  return;
}

//-------------------------------------------------------------------------------------------------
int NamelistEmulator::assignElement(const std::string &key, const std::string &value) {

  // The inputs presumably contain a keyword and value pair.  The question is what to do with the
  // value, and the answer is to find out what sort of input the namelist element expects.
  const size_t param_index = findIndexByKeyword(key);
  if (param_index >= keywords.size()) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
            "assignElement");
    case ExceptionResponse::WARN:
      rtWarn("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
             "assignElement");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return 0;
  }
  const NamelistType param_type = keywords[param_index].kind;
  bool problem = false;
  switch (param_type) {
  case NamelistType::BOOLEAN:
    return 0;
  case NamelistType::STRUCT:
    rtErr(getEnumerationName(param_type) + " keyword \"" + key + "\" should never be handled in "
          "this context.", "NamelistEmulator", "assignElement");
  case NamelistType::INTEGER:
    if (verifyNumberFormat(value.c_str(), NumberFormat::INTEGER)) {
      keywords[param_index].setIntValue(stol(value));
      keywords[param_index].establishment = InputStatus::USER_SPECIFIED;
    }
    else {
      problem = true;
    }
    break;
  case NamelistType::REAL:
    if (verifyNumberFormat(value.c_str(), NumberFormat::STANDARD_REAL) ||
        verifyNumberFormat(value.c_str(), NumberFormat::SCIENTIFIC)) {
      keywords[param_index].setRealValue(stod(value));
      keywords[param_index].establishment = InputStatus::USER_SPECIFIED;
    }
    else {
      problem = true;
    }
    break;
  case NamelistType::STRING:
    keywords[param_index].setStringValue(value);
    keywords[param_index].establishment = InputStatus::USER_SPECIFIED;
    break;
  }

  // Respond to input errors
  if (problem) {
    switch(policy) {
    case ExceptionResponse::DIE:
      rtErr("Keyword \"" + key + "\" in namelist \"" + title + "\" accepts " +
            getEnumerationName(param_type) + " values.  " + value + " is invalid.",
            "NamelistEmulator", "assignElement");
    case ExceptionResponse::WARN:
      rtWarn("Keyword \"" + key + "\" in namelist \"" + title + "\" accepts " +
             getEnumerationName(param_type) + " values.  " + value + " is invalid and no new "
             "value will be assigned.", "NamelistEmulator", "assignElement");
      return 0;
    case ExceptionResponse::SILENT:
      return 0;
    }
  }
  return 1;
}

//-------------------------------------------------------------------------------------------------
int NamelistEmulator::assignElement(const std::string &key, const std::string &sub_key,
                                    const std::string &value) {

  // The inputs contain a keyword, a sub-key, and a value.  This is only valid for STRUCTs.
  const size_t param_index = findIndexByKeyword(key);
  if (param_index >= keywords.size()) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
            "assignElement");
    case ExceptionResponse::WARN:
      rtWarn("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
             "assignElement");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return 0;
  }
  const NamelistType param_type = keywords[param_index].kind;
  bool problem = false;
  switch (param_type) {
  case NamelistType::STRUCT:
    {
      // Search within this keyword.  Use the getter function for encapsulated functionality,
      // even though friendship lets this NamelistEmulator go right in.
      NamelistType subtype = keywords[param_index].getKind(sub_key);
      bool problem = false;
      switch (subtype) {
      case NamelistType::STRUCT:
        break;
      case NamelistType::BOOLEAN:
        return 0;
      case NamelistType::INTEGER:
        if (verifyNumberFormat(value.c_str(), NumberFormat::INTEGER)) {
          keywords[param_index].setIntValue(sub_key, stol(value));
        }
        else {
          problem = true;
        }
        break;
      case NamelistType::REAL:
        if (verifyNumberFormat(value.c_str(), NumberFormat::STANDARD_REAL) ||
            verifyNumberFormat(value.c_str(), NumberFormat::SCIENTIFIC)) {
          keywords[param_index].setRealValue(sub_key, stod(value));
        }
        else {
          problem = true;
        }
        break;
      case NamelistType::STRING:
        keywords[param_index].setStringValue(sub_key, value);
      }

      // Respond to input errors
      if (problem) {
        switch(policy) {
        case ExceptionResponse::DIE:
          rtErr("In namelist \"" + title + "\", keyword \"" + key + "\", sub-key \"" +
                sub_key + "\" accepts " + getEnumerationName(subtype) + " values.  " +
                value + " is invalid.", "NamelistEmulator", "assignElement");
        case ExceptionResponse::WARN:
          rtWarn("In namelist \"" + title + "\", keyword \"" + key + "\", sub-key \"" +
                 sub_key + "\" accepts " + getEnumerationName(subtype) + " values.  " +
                 value + " is invalid and no new value will be assigned.", "NamelistEmulator",
                 "assignElement");
          return 0;
        case ExceptionResponse::SILENT:
          return 0;
        }
      }
    }
    break;
  case NamelistType::BOOLEAN:
    return 0;
  case NamelistType::INTEGER:
  case NamelistType::REAL:
  case NamelistType::STRING:
    rtErr(getEnumerationName(param_type) + " keyword \"" + key +
          "\" should never be handled in this context.", "NamelistEmulator", "assignElement");
  }
  return 1;
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::triggerResizeBuffer(const std::string &key) {
  const size_t param_index = findIndexByKeyword(key);

  // No exception handling is needed here.  Problems will be caught earlier by assignElement().
  // Check that the keyword is part of the namelist, however.  If the keyword is associated with
  // a STRUCT, which is the purpose of providing this function for NamelistEmulators to control
  // advancement of their unerlying NamelistElement entry counts, ensure that all sub-keys in the
  // STRUCT either have default values or have been specified by the user.
  if (param_index < keywords.size()) {
    if (keywords[param_index].kind == NamelistType::STRUCT) {
      bool advance = false;
      bool problem = false;
      const int est_offset = keywords[param_index].next_entry_index *
                             keywords[param_index].template_size;
      for (int i = 0; i < keywords[param_index].template_size; i++) {
        const size_t eopi = est_offset + i;
        if (keywords[param_index].sub_key_found[i] == false) {

          // The subkey was not specified by the user.  Apply the template default, if available,
          // or issue an alert.  Continue issuing alerts until the entire STRUCT has been examined,
          // then trigger a runtime error.
          if (keywords[param_index].template_establishment[i] == InputStatus::MISSING) {
            keywords[param_index].instance_establishment[eopi] = InputStatus::MISSING;
            if (keywords[param_index].template_requirements[i] == KeyRequirement::REQUIRED) {
              rtAlert("Incomplete STRUCT entry for keyword \"" + key + "\" of namelist \"" +
                      title + "\": sub-key \"" + keywords[param_index].sub_keys[i] +
                      "\" does not have a user specification or a default value.",
                      "NamelistEmulator", "triggerResizeBuffer");              
              problem = true;
            }
          }
          else {
            keywords[param_index].instance_establishment[eopi] = InputStatus::DEFAULT;

            // The STRUCT keyword, as a whole, will be described as USER_SPECIFIED if any subkey
            // for at least one instance of the keyword has been specified by the user.
            if (keywords[param_index].establishment == InputStatus::MISSING) {
              keywords[param_index].establishment = InputStatus::DEFAULT;
            }
            
            // Without checking which type of variable this subkey pertains to, apply all of the
            // template values for the ith position to their respective arrays.
            const int tmp_int_default         = keywords[param_index].template_ints[i];
            const double tmp_real_default     = keywords[param_index].template_reals[i];
            const std::string tmp_str_default = keywords[param_index].template_strings[i];
            keywords[param_index].sub_int_values[eopi]    = tmp_int_default;
            keywords[param_index].sub_real_values[eopi]   = tmp_real_default;
            keywords[param_index].sub_string_values[eopi] = tmp_str_default;
            advance = true;
          }
        }
        else {

          // Check that the subkey is valid, then mark the user-specified input
          if (keywords[param_index].template_requirements[i] == KeyRequirement::BOGUS) {
            rtAlert("Situational violation in STRUCT keyword \"" + key + "\" of namelist \"" +
                    title + "\": sub-key \"" + keywords[param_index].sub_keys[i] +
                    "\" is bogus in this context.", "NamelistEmulator", "triggerResizeBuffer");
            problem = true;
          }
          else {
            keywords[param_index].establishment = InputStatus::USER_SPECIFIED;
            keywords[param_index].instance_establishment[eopi] = InputStatus::USER_SPECIFIED;
            advance = true;
          }
        }
      }
      if (problem) {
        rtErr("Incomplete specification of data associated with a STRUCT keyword.",
              "NamelistEmulator", "triggerResizeBuffer");
      }
      else if (advance) {
        keywords[param_index].resizeBuffer();
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addHelp(const std::string &blurb) {
  help_message = blurb;
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addHelp(const std::string &key, const std::string &blurb) {
  const size_t p_index = findIndexByKeyword(key);
  if (p_index >= keywords.size()) {
    rtErr("No keyword \"" + key + "\" exists in namelist \"" + title + "\".",
          "NamelistEmulator", "addHelp");
  }
  keywords[p_index].help_message = blurb;
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addHelp(const std::string &key, const std::string &sub_key,
                               const std::string &blurb) {
  const size_t p_index = findIndexByKeyword(key);
  if (p_index >= keywords.size()) {
    rtErr("No keyword \"" + key + "\" is present in namelist \"" + title + "\".",
          "NamelistEmulator", "addHelp");
  }
  if (keywords[p_index].kind != NamelistType::STRUCT) {
    rtErr("The keyword \"" + key + "\" accepts an " +
          getEnumerationName(keywords[p_index].kind) + " input.  It must be a STRUCT in "
          "order to incorporate documentation for a member variable's documentation.");
  }
  const size_t member_index = findStringInVector(keywords[p_index].sub_keys, sub_key);
  if (member_index >= keywords[p_index].sub_help_messages.size()) {
    rtErr("No sub-key \"" + sub_key + "\" is present in namelist \"" + title + "\" keyword \"" +
          key + "\".", "NamelistEmulator", "addHelp");
  }
  keywords[p_index].sub_help_messages[member_index] = blurb;
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::addCategory(const std::string &new_category) {
  category_names.push_back(new_category);
}
  
//-------------------------------------------------------------------------------------------------
void NamelistEmulator::categorizeKeyword(const std::string &key,
                                         const std::string &category_label) {

  // Ensure that the keyword and category labels are valid
  const size_t category_index = findStringInVector(category_names, category_label);
  if (category_index == category_names.size()) {
    rtErr("Namelist \"" + title + "\" has no category \"" + category_label + "\".",
          "NamelistEmulator", "categorizeKeyword");
  }
  const size_t p_index = findIndexByKeyword(key);
  if (p_index < keywords.size()) {
    categories[category_index].push_back(key);
  }
  else {
    rtErr("Namelist \"" + title + "\" has no keyword \"" + key + "\".", "NamelistEmulator",
          "categorizeKeyword");
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::setImperative(const std::string &key, const KeyRequirement req) {
  const size_t p_index = findIndexByKeyword(key);
  if (p_index >= keywords.size()) {
    rtErr("No keyword \"" + key + "\" is present in namelist \"" + title + "\".",
          "NamelistEmulator", "setImperative");
  }
  keywords[p_index].setImperative(req);  
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::setImperative(const std::vector<std::string> &directives) {
  const int ndir = directives.size();
  for (int i = 0; i < ndir; i++) {
    const size_t ndchar = directives[i].size() - 1;
    if (directives[i][ndchar] == 'r' || directives[i][ndchar] == 'o' ||
        directives[i][ndchar] == 'g') {

      // Try interpreting the directive as a "keyword / requirement" fusion, where the name of a
      // valid keyword is followed by 'r' (REQUIRED), 'o' (OPTIONAL), or 'g' (BOGUS).  If this
      // interpretation succeeds, let the directive counter advance with the loop control.
      try {
        KeyRequirement req;
        switch (directives[i][ndchar]) {
        case 'r':
          req = KeyRequirement::REQUIRED;
          break;
        case 'o':
          req = KeyRequirement::OPTIONAL;
          break;
        case 'g':
          req = KeyRequirement::BOGUS;
          break;
        default:

          // This case cannot be reached.
          break;
        }
        setImperative(directives[i].substr(0, ndchar), req);
        continue;
      }
      catch (std::runtime_error) {}
    }
    if (i < ndir - 1) {

      // Try interpreting the directive as a keyword and the one directly behind it as the
      // requirement code.  If this succeeds, advance the loop counter one extra time.
      try {
        setImperative(directives[i], translateKeyRequirement(directives[i + 1]));
        i++;
        continue;
      }
      catch (std::runtime_error) {}
    }

    // This situation indicates an error in parsing the directives.
    rtErr("The directive " + directives[i] + " is invalid for setting input requirements.",
          "NamelistEmulator", "setImperative");
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::setImperative(const std::string &key,
                                     const std::vector<std::string> &directives) {
  const size_t p_index = findIndexByKeyword(key);
  if (p_index >= keywords.size()) {
    rtErr("No keyword \"" + key + "\" is present in namelist \"" + title + "\".",
          "NamelistEmulator", "setImperative");
  }
  const int ndir = directives.size();
  for (int i = 0; i < ndir; i++) {
    const size_t ndchar = directives[i].size() - 1;
    if (directives[i][ndchar] == 'r' || directives[i][ndchar] == 'o' ||
        directives[i][ndchar] == 'g') {

      // Try interpreting the directive as a "keyword / requirement" fusion, where the name of a
      // valid keyword is followed by 'r' (REQUIRED), 'o' (OPTIONAL), or 'g' (BOGUS).  If this
      // interpretation succeeds, let the directive counter advance with the loop control.
      try {
        KeyRequirement req;
        switch (directives[i][ndchar]) {
        case 'r':
          req = KeyRequirement::REQUIRED;
          break;
        case 'o':
          req = KeyRequirement::OPTIONAL;
          break;
        case 'g':
          req = KeyRequirement::BOGUS;
          break;
        default:

          // This case cannot be reached.
          break;
        }
        keywords[p_index].setImperative(directives[i].substr(0, ndchar), req);
        continue;
      }
      catch (std::runtime_error) {}
    }
    if (i < ndir - 1) {

      // Try interpreting the directive as a keyword and the one directly behind it as the
      // requirement code.  If this succeeds, advance the loop counter one extra time.
      try {
        keywords[p_index].setImperative(directives[i], translateKeyRequirement(directives[i + 1]));
        i++;
        continue;
      }
      catch (std::runtime_error) {}
    }

    // This situation indicates an error in parsing the directives.
    rtErr("The directive " + directives[i] + " is invalid for setting input requirements in "
          "keyword \"" + key + "\".", "NamelistEmulator", "setImperative");
  }
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::printKeywordDocumentation(const int p_idx, const int name_width,
                                                 const int kw_kind_width,
                                                 const std::string &kw_dflt) const {
  std::string kw_kind = getEnumerationName(keywords[p_idx].kind);
  kw_kind = addLeadingWhiteSpace(kw_kind, kw_kind_width);
  printf(" + %-*.*s : %s\n", name_width, name_width, keywords[p_idx].label.c_str(),
         terminalFormat("[" + kw_kind + ", " + kw_dflt + "] " + keywords[p_idx].help_message, "",
                        "", 6 + name_width, 0, 6 + name_width, 0, RTMessageKind::TABULAR).c_str());
  if (keywords[p_idx].kind == NamelistType::STRUCT) {
    int dflt_width = 0;
    int member_width = 0;
    int kind_width = 0;
    std::vector<std::string> kind_pres(keywords[p_idx].template_size);
    for (int i = 0; i < keywords[p_idx].template_size; i++) {
      member_width = std::max(member_width, static_cast<int>(keywords[p_idx].sub_keys[i].size()));
      kind_pres[i] = getEnumerationName(keywords[p_idx].sub_kinds[i]);
      kind_width = std::max(kind_width, static_cast<int>(kind_pres[i].size()));
    }
    std::vector<std::string> dflt_pres(keywords[p_idx].template_size);
    for (int i = 0; i < keywords[p_idx].template_size; i++) {
      switch (keywords[p_idx].template_establishment[i]) {
      case InputStatus::MISSING:
        dflt_pres[i] = "None";
        break;
      case InputStatus::DEFAULT:
        switch (keywords[p_idx].sub_kinds[i]) {
        case NamelistType::BOOLEAN:
          dflt_pres[i] = true;
          break;
        case NamelistType::INTEGER:
          dflt_pres[i] = std::to_string(keywords[p_idx].template_ints[i]);
          break;
        case NamelistType::REAL:
          dflt_pres[i] = minimalRealFormat(keywords[p_idx].template_reals[i], 1.0e-4, true);
          break;
        case NamelistType::STRING:
          dflt_pres[i] = "'" + keywords[p_idx].template_strings[i] + "'";
          break;
        case NamelistType::STRUCT:
          break;
        }
        break;
      case InputStatus::USER_SPECIFIED:

        // Keywords either have default values or are missing default values.
        break;
      }
      dflt_width = std::max(dflt_width, static_cast<int>(dflt_pres[i].size()));
    }
    for (int i = 0; i < keywords[p_idx].template_size; i++) {
      dflt_pres[i] = addLeadingWhiteSpace(dflt_pres[i], dflt_width);
      kind_pres[i] = addLeadingWhiteSpace(kind_pres[i], kind_width);
      printf("   - %-*.*s : %s\n", member_width, member_width,
             keywords[p_idx].sub_keys[i].c_str(),
             terminalFormat("[" + kind_pres[i] + ", " + dflt_pres[i] + "] " +
                            keywords[p_idx].sub_help_messages[i], "", "", 8 + member_width, 0,
                            8 + member_width, 0, RTMessageKind::TABULAR).c_str());
    }
  }
}

//-------------------------------------------------------------------------------------------------
std::string NamelistEmulator::convertDefaultToString(const NamelistElement &tkw) const {
  std::string result;
  int ndval;
  switch (tkw.establishment) {
  case InputStatus::MISSING:
    return std::string("None");
  case InputStatus::DEFAULT:
    switch (tkw.kind) {
    case NamelistType::BOOLEAN:
      result += "FALSE";
      break;
    case NamelistType::INTEGER:
      ndval = tkw.default_int_values.size();
      for (int j = 0; j < ndval; j++) {
        result += std::to_string(tkw.default_int_values[j]);
        if (j < ndval - 1) {
          result += ", ";
        }
      }
      break;
    case NamelistType::REAL:
      ndval = tkw.default_real_values.size();
      for (int j = 0; j < ndval; j++) {
        result += minimalRealFormat(tkw.default_real_values[j], 1.0e-4, true);
        if (j < ndval - 1) {
          result += ", ";
        }
      }
      break;
    case NamelistType::STRING:
      ndval = tkw.default_string_values.size();
      for (int j = 0; j < ndval; j++) {
        result += "'" + tkw.default_string_values[j] + "'";
        if (j < ndval - 1) {
          result += ", ";
        }
      }
      break;
    case NamelistType::STRUCT:

      // STRUCT-type keywords' default values are stored in the template_(ints / reals / strings)
      // arrays.
      break;
    }
    break;
  case InputStatus::USER_SPECIFIED:

    // This case should never be reached
    break;
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
void NamelistEmulator::printHelp() const {
  terminalHorizontalRule();
  const std::string prefix_str = (cli_content) ? "" : "&";
  printf("%s\n", terminalFormat(prefix_str + title + ": " + help_message, "", "", 0, 0,
                                title.size() + 3 - static_cast<int>(cli_content), 0,
                                RTMessageKind::TABULAR).c_str());
  int param_width = 0;
  int param_kind_width = 0;
  std::vector<std::string> param_defaults;
  if (categories.size() == 0) {
    if (cli_content) {
      printf("\n Command line inputs");
    }
    else {
      printf("\n Keywords");
    }
    printf(" [ type, default value ]:\n");
    terminalHorizontalRule();
    const int n_params = keywords.size();

    // Make a list of the parameter defaults and determine formatting
    param_defaults.resize(n_params);
    for (int i = 0; i < n_params; i++) {
      param_defaults[i] = convertDefaultToString(keywords[i]);
    }
    justifyStrings(&param_defaults, JustifyText::RIGHT, 7, 3, 2);

    // Determine the maximum width of any keyword
    for (int i = 0; i < n_params; i++) {
      param_width = std::max(param_width, static_cast<int>(keywords[i].label.size()));
      param_kind_width = std::max(param_kind_width,
                                  static_cast<int>(getEnumerationName(keywords[i].kind).size()));
    }
    for (int i = 0; i < n_params; i++) {
      printKeywordDocumentation(i, param_width, param_kind_width, param_defaults[i]);
      if (i < n_params - 1) {
        printf("\n");
      }
    }
  }
  else {
    const int n_categories = categories.size();
    for (int i = 0; i < n_categories; i++) {
      printf("\n %s Keywords [ type, default value ]:\n", category_names[i].c_str());
      const int n_cati = categories[i].size();

      // Make a list of the parameter defaults and determine formatting
      param_defaults.resize(n_cati);
      for (int j = 0; j < n_cati; j++) {
        param_defaults[j] = convertDefaultToString(keywords[findIndexByKeyword(categories[i][j])]);
      }
      justifyStrings(&param_defaults, JustifyText::RIGHT, 7);

      // Determine the maximum width of any keyword
      int param_width = 0;
      for (int j = 0; j < n_cati; j++) {
        param_width = std::max(param_width, static_cast<int>(categories[i][j].size()));
      }
      for (int j = 0; j < n_cati; j++) {
        // Safe to use the output of findIndexByKeyword() as an array index because the contents of
        // every category have already been checked as valid members of the namelist
        printKeywordDocumentation(findIndexByKeyword(categories[i][j]), param_width,
                                  param_kind_width, param_defaults[j]);
        if (i < n_cati - 1) {
          printf("\n");
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
std::string NamelistEmulator::printContents(const int file_width, const int max_entry_reports,
                                            const NamelistIntroduction print_decor) const {
  int key_width = 0;
  int kind_width = 0;
  int source_width = 0;
  int value_width = 20;
  const std::string needs_ui("[ Needs user input ]");
  const std::string blank_mark;
  const std::string skip_mark("  (...)");
  int nkw = keywords.size();
  int nkw_elements = 1;
  for (int i = 0; i < nkw; i++) {
    int nrep;
    if (keywords[i].getEstablishment() == InputStatus::MISSING) {
      nrep = 1;      
    }
    else {
      nrep = std::min(max_entry_reports, keywords[i].getEntryCount());
      nkw_elements += (max_entry_reports < keywords[i].getEntryCount());
    }
    switch (keywords[i].kind) {
    case NamelistType::BOOLEAN:
    case NamelistType::INTEGER:
    case NamelistType::REAL:
    case NamelistType::STRING:
      nkw_elements += nrep;
      break;
    case NamelistType::STRUCT:
      nkw_elements += 1 + (nrep * keywords[i].getTemplateSize());
      break;
    }
  }
  std::vector<std::string> kw_names, kw_kinds, kw_sources, kw_repeats, kw_values;
  kw_names.reserve(nkw_elements);
  kw_kinds.reserve(nkw_elements);
  kw_sources.reserve(nkw_elements);
  kw_values.reserve(nkw_elements);
  switch (print_decor) {
  case NamelistIntroduction::HEADER:
  case NamelistIntroduction::COMPACT_HEADER:
    kw_names.push_back(std::string("Keyword"));
    kw_kinds.push_back(std::string("Kind"));
    kw_sources.push_back(std::string("Source"));
    kw_repeats.push_back(std::string("Rep"));
    kw_values.push_back(std::string("Value"));
    break;
  case NamelistIntroduction::BLANK_LINE:
  case NamelistIntroduction::NONE:
    break;
  }
  bool multiple_entries = false;
  for (int i = 0; i < nkw; i++) {

    // Skip bogus keywords
    if (keywords[i].getImperative() == KeyRequirement::BOGUS) {
      continue;
    }
    
    // The first line of each entry is devoted to the keyword and, unless it is a STRUCT, the
    // first value of that keyword.
    kw_names.push_back(keywords[i].getLabel());
    kw_kinds.push_back(getEnumerationName(keywords[i].getKind()));
    kw_sources.push_back(getEnumerationName(keywords[i].getEstablishment()));
    if (keywords[i].getEstablishment() == InputStatus::MISSING) {
      kw_values.push_back(needs_ui);
      kw_repeats.push_back(blank_mark);
    }
    else {

      // Detail the first instance of the keyword
      switch (keywords[i].getKind()) {
      case NamelistType::BOOLEAN:
        kw_values.push_back(keywords[i].getBoolValue() ? std::string("true") :
                                                         std::string("false"));
        kw_repeats.push_back(std::to_string(1));
        break;
      case NamelistType::INTEGER:
        kw_values.push_back(std::to_string(keywords[i].getIntValue(0)));
        kw_repeats.push_back(std::to_string(1));
        break;
      case NamelistType::REAL:
        kw_values.push_back(minimalRealFormat(keywords[i].getRealValue(0), 1.0e-4));
        kw_repeats.push_back(std::to_string(1));
        break;
      case NamelistType::STRING:
        kw_values.push_back(keywords[i].getStringValue(0));
        kw_repeats.push_back(std::to_string(1));
        break;
      case NamelistType::STRUCT:
        kw_values.push_back(blank_mark);
        kw_repeats.push_back(blank_mark);
        for (int j = 0; j < keywords[i].getTemplateSize(); j++) {

          // Skip bogus subkeys
          if (keywords[i].template_requirements[j] == KeyRequirement::BOGUS) {
            continue;
          }
          const std::string &jsl = keywords[i].getSubLabel(j);
          const NamelistType jkind = keywords[i].getKind(jsl);
          kw_names.push_back("  + " + jsl);
          kw_kinds.push_back(getEnumerationName(jkind));
          kw_sources.push_back(getEnumerationName(keywords[i].getEstablishment(jsl, 0)));
          kw_repeats.push_back(std::to_string(1));
          switch (jkind) {
          case NamelistType::BOOLEAN:
            kw_values.push_back(std::to_string(keywords[i].getBoolValue(jsl, 0)));
            break;
          case NamelistType::INTEGER:
            kw_values.push_back(std::to_string(keywords[i].getIntValue(jsl, 0)));
            break;
          case NamelistType::REAL:
            kw_values.push_back(minimalRealFormat(keywords[i].getRealValue(jsl, 0), 1.0e-4));
            break;
          case NamelistType::STRING:
            kw_values.push_back(keywords[i].getStringValue(jsl, 0));
            break;
          case NamelistType::STRUCT:
            break;
          }
        }
        break;
      }

      // Detail subsequent instances of the keyword
      const bool skip_reports = (keywords[i].getEntryCount() > max_entry_reports);
      bool skip_placed = false;
      const int n_entry = keywords[i].getEntryCount();
      for (int j = 1; j < n_entry; j++) {
        multiple_entries = true;
        if (skip_reports && j > max_entry_reports - 2 && j < n_entry - 1) {
          if (skip_placed == false) {
            kw_names.push_back(skip_mark);
            kw_kinds.push_back(skip_mark);
            kw_sources.push_back(skip_mark);
            kw_repeats.push_back(blank_mark);
            kw_values.push_back(skip_mark);
            skip_placed = true;
          }
        }
        else {
          switch (keywords[i].getKind()) {
          case NamelistType::BOOLEAN:
            break;
          case NamelistType::INTEGER:
          case NamelistType::REAL:
          case NamelistType::STRING:
            kw_names.push_back(blank_mark);
            kw_kinds.push_back(blank_mark);
            kw_sources.push_back(blank_mark);
            kw_repeats.push_back(std::to_string(j + 1));
            break;
          case NamelistType::STRUCT:
            break;
          }
          switch (keywords[i].getKind()) {
          case NamelistType::BOOLEAN:
            break;
          case NamelistType::INTEGER:
            kw_values.push_back(std::to_string(keywords[i].getIntValue(j)));
            break;
          case NamelistType::REAL:
            kw_values.push_back(minimalRealFormat(keywords[i].getRealValue(j), 1.0e-4));
            break;
          case NamelistType::STRING:
            kw_values.push_back(keywords[i].getStringValue(j));
            break;
          case NamelistType::STRUCT:
            for (int k = 0; k < keywords[i].getTemplateSize(); k++) {

              // Skip bogus subkeys
              if (keywords[i].template_requirements[k] == KeyRequirement::BOGUS) {
                continue;
              }
              const std::string &ksl = keywords[i].getSubLabel(k);
              const NamelistType kkind = keywords[i].getKind(ksl);
              kw_names.push_back("  + " + ksl);
              kw_kinds.push_back(getEnumerationName(kkind));
              kw_sources.push_back(getEnumerationName(keywords[i].getEstablishment(ksl, j)));
              kw_repeats.push_back(std::to_string(j + 1));
              switch (kkind) {
              case NamelistType::BOOLEAN:
                kw_values.push_back(std::to_string(keywords[i].getBoolValue(ksl, j)));
                break;
              case NamelistType::INTEGER:
                kw_values.push_back(std::to_string(keywords[i].getIntValue(ksl, j)));
                break;
              case NamelistType::REAL:
                kw_values.push_back(minimalRealFormat(keywords[i].getRealValue(ksl, j), 1.0e-4));
                break;
              case NamelistType::STRING:
                kw_values.push_back(keywords[i].getStringValue(ksl, j));
                break;
              case NamelistType::STRUCT:
                break;
              }
            }
            break;
          }
        }
      }
    }
  }

  // Improve the presentation by de-emphasizing variable tyes and abridging other terms
  const std::string user_source = getEnumerationName(InputStatus::USER_SPECIFIED);
  for (size_t i = 1; i < kw_sources.size(); i++) {
    if (strcmpCased(kw_sources[i], user_source)) {
      kw_sources[i] = "USER";
    }
    kw_sources[i] = lowercase(kw_sources[i]);
  }
  for (size_t i = 0; i < kw_kinds.size(); i++) {
    kw_kinds[i] = lowercase(kw_kinds[i]);
  }
  
  // Justify the non-value components of each input
  const size_t ntab_ent = kw_names.size();
  justifyStrings(&kw_names, 0, ntab_ent, JustifyText::LEFT, 20, 0, 20);
  justifyStrings(&kw_kinds, 0, ntab_ent, JustifyText::LEFT, 20, 0, 20);
  justifyStrings(&kw_sources, 0, ntab_ent, JustifyText::LEFT, 20, 0, 20);
  justifyStrings(&kw_repeats, JustifyText::RIGHT);

  // Print the table headings
  std::string result;
  switch (print_decor) {
  case NamelistIntroduction::HEADER:
    result += horizontalRule("+", "+", file_width);
    result += "Contents of &" + title + "\n\n";
    break;
  case NamelistIntroduction::COMPACT_HEADER:
    result += "Contents of &" + title + "\n";
    break;
  case NamelistIntroduction::BLANK_LINE:
    result += '\n';
    break;
  case NamelistIntroduction::NONE:
    break;
  }
  for (size_t i = 0; i < ntab_ent; i++) {
    result += " " + kw_names[i] + " | " + kw_kinds[i] + " | " + kw_sources[i] + " | ";
    if (multiple_entries) {
      result += kw_repeats[i] + " | ";
    }
    result += kw_values[i] + "\n";
    if (i == 0) {
      switch (print_decor) {
      case NamelistIntroduction::HEADER:
      case NamelistIntroduction::COMPACT_HEADER:
        {
          result += " " + std::string(kw_names[0].size(), '-') + "-+-" +
                    std::string(kw_kinds[0].size(), '-') + "-+-" +
                    std::string(kw_sources[0].size(), '-') + "-+-";
          if (multiple_entries) {
            result += std::string(kw_repeats[0].size(), '-') + "-+-";
          }
          size_t mv_len = 0;
          for (size_t j = 0; j < ntab_ent; j++) {
            mv_len = std::max(mv_len, kw_values[j].size());
          }
          result += std::string(mv_len, '-') + "\n";
        }
        break;
      case NamelistIntroduction::BLANK_LINE:
      case NamelistIntroduction::NONE:
        break;
      }
    }
  }
  result += "\n";
  return result;
}
  
//-------------------------------------------------------------------------------------------------
void NamelistEmulator::printContents(std::ostream *foutp, const int file_width,
                                     const int max_entry_reports,
                                     const NamelistIntroduction print_decor) const {
  const std::string result = printContents(file_width, max_entry_reports, print_decor);
  foutp->write(result.data(), result.size());
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::printContents(std::ofstream *foutp, const int file_width,
                                     const int max_entry_reports,
                                     const NamelistIntroduction print_decor) const {
  const std::string result = printContents(file_width, max_entry_reports, print_decor);
  foutp->write(result.data(), result.size());
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::printContents(const std::string &file_name,
                                     const PrintSituation expectation, const int file_width,
                                     const int max_entry_reports,
                                     const NamelistIntroduction print_decor) const {
  std::ofstream foutp = openOutputFile(file_name, expectation, "print &namelist information to a "
                                       "file detailing up to " +
                                       std::to_string(max_entry_reports) + " entry repetitions");
  const std::string result = printContents(file_width, max_entry_reports, print_decor);
  foutp.write(result.data(), result.size());
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
int NamelistEmulator::findIndexByKeyword(const std::string &query) const {
  const int n_element = keywords.size();
  if (n_element == 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Namelist \"" + title + "\" contains no keywords with which to match keyword " +
            query + ".", "NamelistEmulator", "findIndexByKeyword");
    case ExceptionResponse::WARN:
      rtWarn("Namelist \"" + title + "\" contains no keywords with which to match keyword " +
             query + ".", "NamelistEmulator", "findIndexByKeyword");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  const int n_letters = query.size();
  if (n_letters == 0) {
    rtErr("A blank namelist (\"" + title + "\") is not searchable.", "NamelistEmulator",
          "findIndexByKeyword");
  }
  for (int i = 0; i < n_element; i++) {
    if (keywords[i].label.size() == n_letters && keywords[i].label[0] == query[0] &&
        keywords[i].label == query) {
      return i;
    }
  }
  
  // If this point is reached, return the end of the list (code that accepts this value will need
  // to check for this possibility and not act on it--there are cases where the code must always
  // throw an exception lest it cause a segmentation fault, or cases where it is OK to warn the
  // user (or even remain silent) and go on.
  return n_element;
}

//-------------------------------------------------------------------------------------------------
void NamelistEmulator::verifyEstablishment(const std::string &keyword_query, const size_t p_index,
                                           const char* caller) const {
  if (keywords[p_index].getEstablishment() == InputStatus::MISSING) {
    rtErr("Namelist \"" + title + "\" keyword \"" + keyword_query + "\" has not been set by "
          "default or by the user.", "NamelistEmulator", caller);
  }
}
  
} // namespace namelist
} // namespace stormm
