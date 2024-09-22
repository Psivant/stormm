// -*-c++-*-
#ifndef STORMM_NAMELIST_ELEMENT_H
#define STORMM_NAMELIST_ELEMENT_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "namelist_enumerators.h"

namespace stormm {
namespace namelist {

using constants::ExceptionResponse;

/// \brief One keyword found in a namelist, ready to store the namelist variable moniker, the
///        type, and the value read from the input file.
class NamelistElement {
public:

  /// \brief Constructor for a non-STRUCT NamelistElement
  ///
  /// Overloaded:
  ///   - Constructor for INTEGER, REAL, and STRING keywords (not STRUCTs)
  ///   - Constructor for STRUCT keywords
  ///
  /// \param keyword_in     The name of the namelist variable
  /// \param sub_keys_in    The name of STRUCT members within the namelist variable
  /// \param kind_in        The kind of the namelist element (it will be checked to ensure that it
  ///                       is not STRUCT)
  /// \param sub_kinds_in   The kinds of STRUCT members (the namelist variable itself is a STRUCT
  ///                       if this input argument is a vector)
  /// \param rep_in         Flag to indicate whether the namelist input can repeat for additional
  ///                       values
  /// \param help_in        User documentation for the keyword
  /// \param sub_help_in    User documentation for sub-keys
  /// \{
  NamelistElement(const std::string &keyword_in, NamelistType kind_in,
                  const std::string &default_in = std::string(""),
                  DefaultIsObligatory obligate = DefaultIsObligatory::NO,
                  InputRepeats rep_in = InputRepeats::NO,
                  const std::string &help_in = std::string("No description provided"));

  NamelistElement(const std::string keyword_in, const std::vector<std::string> &sub_keys_in,
                  const std::vector<NamelistType> &sub_kinds_in,
                  const std::vector<std::string> &default_list,
                  DefaultIsObligatory obligate_list = DefaultIsObligatory::NO,
                  InputRepeats rep_in = InputRepeats::NO,
                  const std::string &help_in = std::string("No description provided"),
                  const std::vector<std::string> &sub_help_in =
                  std::vector<std::string>(1, "No description provided"),
                  const std::vector<KeyRequirement> &template_requirements_in = {});
  /// \}

  /// \brief With no const members and only Standard Template Library components, the default copy
  ///        and move constructors, as well as the copy and move assignment operators, are valid.
  ///
  /// \param original  The pre-existing object to copy or move
  /// \param other     A pre-existing object to fulfill the right hand side of an assignment
  ///                  statement
  /// \{
  NamelistElement(const NamelistElement &original) = default;
  NamelistElement(NamelistElement &&original) = default;
  NamelistElement& operator=(const NamelistElement &original) = default;
  NamelistElement& operator=(NamelistElement &&original) = default;
  /// \}
  
  /// \brief Get the keyword for a namelist element, i.e. nstlim in Amber &ctrl
  const std::string& getLabel() const;

  /// \brief Get a sub-key for this element
  const std::string& getSubLabel(size_t index) const;

  /// \brief Get the kind associated with a namelist element, i.e. nstlim is an INTEGER in Amber
  ///        &ctrl.
  ///
  /// Overloaded:
  ///   - Get the kind of the element itself
  ///   - Take the name of a sub-key to retrieve the kind of a member variable from a STRUCT
  /// \{
  NamelistType getKind() const;
  NamelistType getKind(const std::string &sub_key_query) const;
  /// \}

  /// \brief Set the way that this Namelist element will respond if it encounters bad input.  The
  ///        member variable associated with this function is not passed down as part of the
  ///        constructor argument list, and would be tedious for the developer to enter over and
  ///        over, so it is passed down based on the NamelistEmulator's own setting by calling
  ///        this function.
  ///
  /// \param policy_in  The policy to set
  void setPolicy(ExceptionResponse policy_in);

  /// \brief Report an error based on an incorrect namelist element data type request.  This is an
  ///        assertion that debugs a program when a developer makes a mistake, not something that
  ///        an end user should encounter, but it's helpful to have a more detailed description
  ///        than assert() can provide.  This always ends by throwing a runtime error and so is an
  ///        acceptable way to fill a switch case.
  ///
  /// \param caller     The name of the member function of NamelistElement calling this error
  /// \param data_type  The name of the (erroneously) requested data type
  void reportNamelistTypeProblem(const std::string &caller, const std::string &data_type) const;

  /// \brief Get the depth of the namelist element, the number of values that it stores.
  int getEntryCount() const;

  /// \brief Get the boolean value to which a namelist element has been set (this value is set to
  ///        TRUE if the keyword was found in the namelist or collection of STRUCT sub-keys, FALSE
  ///        if not).
  ///
  /// Overloaded:
  ///   - Take a STRUCT member variable name and the index of the particular STRUCT entry
  ///   - Return the boolean value of the one and only valid entry of a non-STRUCT keyword
  ///   - Return all BOOL values associated with a boolean member of a STRUCT keyword, which
  ///     could have multiple entries
  ///
  /// \param member_key    The (optional) name of a STRUCT member to access (an empty string
  ///                      indicates that the keyword is associated with INTEGER, REAL, or STRING
  ///                      data)
  /// \param index         The entry to access
  /// \{
  bool getBoolValue(const std::string &member_key, int index) const;
  bool getBoolValue() const;
  std::vector<bool> getBoolValue(const std::string &member_key) const;
  /// \}
  
  /// \brief Get the integer value to which a namelist element has been set (this value is read
  ///        from the namelist input file, i.e. mdin).  Descriptions of input parameters follow
  ///        from getBoolValue(), above.
  ///
  /// Overloaded:
  ///   - Take a STRUCT member variable name and the index of the particular STRUCT entry
  ///   - Take an index only and return the integer value (the namelist element must be a
  ///     non-STRUCT)
  ///   - Return all INTEGER values associated with this keyword
  /// \{
  int getIntValue(const std::string &member_key, int index) const;
  int getIntValue(int index) const;
  std::vector<int> getIntValue(const std::string &member_key = std::string("")) const;
  /// \}

  /// \brief Get the real value to which a namelist element has been set (this value is read from
  ///        the namelist input file, i.e. mdin).  Overloading and descriptions of input parameters
  ///        follow from getIntegerValue()
  /// \{
  double getRealValue(const std::string &member_key, int index) const;
  double getRealValue(int index) const;
  std::vector<double> getRealValue(const std::string &member_key = std::string("")) const;
  /// \}

  /// \brief Get the string value to which a namelist element has been set (this value is read from
  ///        the namelist input file, i.e. mdin)
  ///
  /// Overloaded:
  ///   - Take a STRUCT member variable name and the index of the particular STRUCT entry
  ///   - Take an index only and return the string value (the namelist element must be a
  ///     non-STRUCT)
  ///   - Return all STRING values associated with this keyword
  ///
  /// \param member_key    The (optional) name of a STRUCT member to access
  /// \param index         The entry to access
  /// \{
  const std::string& getStringValue(const std::string &member_key, int index) const;
  const std::string& getStringValue(int index) const;
  std::vector<std::string> getStringValue(const std::string &member_key = std::string("")) const;
  /// \}

  /// \brief Report whether a value has been assigned, read from input, to this namelist element.
  ///
  /// \param member_key  Name of the sub-key to search if the namelist element is a STRUCT
  /// \param repeat_no   Index number of the repeated keyword application, if applicable.  The
  ///                    establishment of a member keyword in the fourth application of a
  ///                    repeatable STRUCT keyword will be found in element 4 * (template_size) +
  ///                    member_idx of the input status array.
  InputStatus getEstablishment(const std::string &member_key = std::string(""),
                               int repeat_no = 0) const;

  /// \brief Report whether a keyword has been deemed essential, optional, or bogus in a particular
  ///        context.  This enables the same namelist to take on different profiles in different
  ///        applications, and to report only that data which is relevant in the program output.
  ///
  /// Overloaded:
  ///   - Get the criticality of the keyword
  ///   - Get the criticality of a subkey within a STRUCT keyword
  ///
  /// \param sub_key_query  Name of the sub-key to search if the namelist element is a STRUCT
  /// \{
  KeyRequirement getImperative() const;
  KeyRequirement getImperative(const std::string &sub_key_query) const;
  /// \}

  /// \brief Report whether an element accepts multiple values
  bool getRepeatableValueState() const;

  /// \brief Obtain the size of this namelist element's template (if it is a STRUCT)
  int getTemplateSize() const;

  /// \brief Respond to bad user input at the NamelistElement level, based on a policy handed down
  ///        from a NamelistEmulator.  This does not apply to developer-level mistakes, which
  ///        always throw runtime errors.
  void badInputResponse(const std::string &errmsg, const char* caller);

  /// \brief Set the default value for a particular keyword, overriding anything that might have
  ///        been set beforehand.  This can be useful when applying common keyword loading
  ///        functions, to adapt the namelist elements for a specific program or context.
  ///
  /// Overloaded:
  ///   - Set the default value for a non-STRUCT keyword
  ///   - Set the default values for various sub-keys of a STRUCT keyword
  ///
  /// \param modified_default   The new default setting to apply to the keyword
  /// \param default_idx        The index of the default to set, in the event that there are
  ///                           already multiple default values
  /// \param modified_defaults  The list of modified default values
  /// \param sub_key_specs      The list of sub-keys to which each default corresponds
  /// \{
  void setDefaultValue(const std::string &modified_default, int default_idx = 0);
  void setDefaultValue(const std::vector<std::string> &modified_defaults,
                       const std::vector<std::string> &sub_key_specs);
  /// \}
  
  /// \brief Include an additional value as a default setting for a particular keyword.  This
  ///        enables a single keyword to have a default series of values.  The supplied value will
  ///        be interpreted according to the type of the element.
  ///
  /// \param next_default  The extra default value to include
  void addDefaultValue(const std::string &next_default);

  /// \brief Activate a boolean value based on the mere presence of a keyword in the namelist.
  ///        This function will check to verify that the keyword is a boolean.
  ///
  /// Overloaded:
  ///   - Take a the keyword name
  ///   - Take a sub-keyword and a value to fill in a BOOL member of a STRUCT kind element
  /// \{
  void activateBool();
  void activateBool(const std::string &su_key_query);
  /// \}
  
  /// \brief Set an integer value within this namelist element.  This function will check to ensure
  ///        that the element expects an int.
  ///
  /// Overloaded:
  ///   - Take a sub-keyword and a value to fill in an INTEGER member of a STRUCT kind element
  ///   - Take a value only
  ///
  /// \param member_key    If given, implies that the namelist element is a struct, and triggers a
  ///                      search through its members to find a match and finally place the int.
  /// \param value         The value of the keyword or sub-key
  /// \{
  void setIntValue(int value);
  void setIntValue(const std::string &member_key, int value);
  /// \}

  /// \brief Set a real value within this namelist element.  This function will check to ensure
  ///        that the element expects a real number.
  ///
  /// Overloaded:
  ///   - Take a sub-keyword and a value to fill in a REAL member of a STRUCT kind element
  ///   - Take a value only
  ///
  /// \param member_key    If given, implies that the namelist element is a struct, and triggers a
  ///                      search through its members to find a match and finally place the real
  ///                      number.
  /// \param value         The value of the keyword or sub-key
  /// \{
  void setRealValue(double value);
  void setRealValue(const std::string &member_key, double value);
  /// \}

  /// \brief Set a string value within this namelist element.  This function will check to ensure
  ///        that the element expects a string value.
  ///
  /// Overloaded:
  ///   - Take a sub-keyword and a value to fill in a STRING member of a STRUCT kind element
  ///   - Take a value only
  ///
  /// \param member_key  If given, implies that the namelist element is a struct, and triggers a
  ///                      search through its members to find a match and finally place the string.
  /// \param value         The value of the keyword or sub-key
  /// \{
  void setStringValue(const std::string &value);
  void setStringValue(const std::string &member_key, const std::string &value);
  /// \}

  /// \brief Set the requirement associated with a keyword.
  ///
  /// Overloaded:
  ///   - Set the criticality of a keyword
  ///   - Set the criticality of a subkey within a STRUCT-type keyword
  ///
  /// \param imperative_in  The criticality level of the keyword to set.  Is it required, optional,
  ///                       or bogus in some context?
  /// \{
  void setImperative(const KeyRequirement imperative_in);
  void setImperative(const std::string &sub_key_query, const KeyRequirement imperative_in);
  /// \}

private:
  std::string label;                          ///< The keyword sought in the input file
  NamelistType kind;                          ///< The kind of namelist variable
  ExceptionResponse policy;                   ///< Indicator of what to do if this namelist element
                                              ///<   encounters an error (passed down from the
                                              ///<   NamelistEmulator containing this
                                              ///<   NamelistElement keyword)
  int next_entry_index;                       ///< The index in any of the value arrays to which
                                              ///<   the next user input shall go.  This is set to
                                              ///<   0 upon initialization but can be set to 1 or
                                              ///<   even higher values in order to protect one or
                                              ///<   more default values from being overwitten by
                                              ///<   further input from the user (e.g. make these
                                              ///<   default values obligatory).
   int entry_count;                           ///< The number of entries that this element is
                                              ///<   officially tracking.  This is separated from
                                              ///<   next_entry_index for the reason that there may
                                              ///<   be a default entry or entries that will be
                                              ///<   superceded by any user-supplied input.  Upon
                                              ///<   receiving any user-supplied input, entry_count
                                              ///<   will update to next_entry_index.
  std::vector<bool> bool_values;              ///< Possible boolean values
  std::vector<int> int_values;                ///< Possible integer values
  std::vector<double> real_values;            ///< Possible real values
  std::vector<std::string> string_values;     ///< Possible string values
  bool accept_multiple_values;                ///< Indicates whether this keyword may be given
                                              ///<   multiple times for unique inputs
  std::string help_message;                   ///< User documentation for this namelist element
  std::vector<std::string> sub_keys;          ///< Labels for members of a STRUCT kind variable
  std::vector<NamelistType> sub_kinds;        ///< List of subkinds if the kind is "STRUCT".  It is
                                              ///<   not feasible to nest STRUCTs, however.
  std::vector<bool> sub_bool_values;          ///< Possible boolean values of STRUCT keywords
  std::vector<int> sub_int_values;            ///< Possible integer values of STRUCT subkinds
  std::vector<double> sub_real_values;        ///< Possible real values of STRUCT subkinds
  std::vector<std::string> sub_string_values; ///< Possible string values of STRUCT subkinds
  std::vector<std::string> sub_help_messages; ///< Help messages for individual member variables
                                              ///<   if this is a STRUCT kind of namelist element
  int template_size;                          ///< The number of fields, if this is a STRUCT kind
                                              ///<   namelist element, for each distinct entry.
                                              ///<   This is taken from sub_keys.size() for
                                              ///<   convenience and clarity.
  std::vector<int> template_ints;             ///< List of ints providing template default
                                              ///<   values for STRUCT kind namelist elements
  std::vector<double> template_reals;         ///< List of real numbers providing template default
                                              ///<   values for STRUCT kind namelist elements
  std::vector<std::string> template_strings;  ///< List of strings providing template default
                                              ///<   values for STRUCT kind namelist elements

  // Just as STRUCT-type keywords have the template_(...) arrays to store the default values of
  // their member subkeys, STRING-, INTEGER-, and REAL-type keywords have special storage for their
  // default(s).  The default value for any BOOL member of a struct, and any BOOL keyword in
  // general, is FALSE.
  std::vector<int> default_int_values;             ///< A record of default integer number values
                                                   ///<   passed to the keyword
  std::vector<double> default_real_values;         ///< A record of default real number values
                                                   ///<   passed to the keyword
  std::vector<std::string> default_string_values;  ///< A record of default string values passed to
                                                   ///<   the keyword
  
  /// Indicates whether a value for the appropriate keyword has been found in the input, taken from
  /// a default value, or is missing
  InputStatus establishment;     

  /// Indicate the criticality of the keyword, enabling contextual dependence for various
  /// &namelists.
  KeyRequirement imperative;
  
  /// Default statuses for each member of a STRUCT kind namelist element.  For STRUCTs accepting
  /// multiple input values, each new instance of the struct will have multiple int, real, or
  /// string variables which need to be specified or assigned default values from the template.
  std::vector<InputStatus> template_establishment;

  /// Requirement settings for each subkey in a STRUCT kind namelist element.  STRUCTS can have
  /// optional keywords with no defaults, and for reuseable code STRUCT elements can accept
  /// different subkey requirements in different contexts.  Subkeys can be flagged as "bogus" if,
  /// in some context, they should not be present.
  std::vector<KeyRequirement> template_requirements;
  
  /// Current status of each member of every application of a STRUCT kind namelist element.
  /// Check this array to see if a particular instance of a STRUCT keyword supplied by the user
  /// has a default, user-specified, or missing subkey value.
  std::vector<InputStatus> instance_establishment;

  /// Status of sub-key searches for the next instance of a STRUCT keyword.  This is kept on hand
  /// to catch double-specification of sub-keys within one set of STRUCT input.
  std::vector<bool> sub_key_found;

  /// \brief Authoritative function to determine whether a NamelistElement is missing a default
  ///        value for its keyword or a sub-key.
  ///
  /// \param default_input
  bool defaultIsMissing(const std::string &default_input) const;

  /// \brief Set the degree of establishment for the namelist element (if it contains INTEGER,
  ///        REAL, or STRING inputs) or the establishments of all member variables if it is a
  ///        STRUCT.
  ///
  /// \param list_of_defaults  The list of default keywords, which will be compared against the
  ///                          list of expected data types in the namelist element
  std::vector<InputStatus> setEstablishment(const std::vector<std::string> &list_of_defaults);

  /// \brief Confirm that a particular sub-key is present in the namelist element (only applicable
  ///        to STRUCT-type elements).  This function will return an error if the sub-key is not
  ///        found, obviating the need for further checks against the index it returns.
  ///
  /// \param member_key    The sub-key in question
  /// \param nmlt          The type of data associated with the sub-key
  /// \param caller        The function calling this check
  int validateSubKey(const std::string &sub_key_query, NamelistType nmlt,
                     const std::string &caller) const;

  /// \brief Resize a namelist buffer in preparation to take the next entry.  The NamelistElement
  ///        will have just had a new value set.  This function will also increment the next entry
  ///        counter and the total entry count.
  void resizeBuffer();

  // The namelist emulator is a friend, so that excessive copying doesn't have to happen just
  // to search and access data
  friend struct NamelistEmulator;
};

} // namespace namelist
} // namespace stormm

#endif
