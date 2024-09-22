// -*-c++-*-
#ifndef STORMM_NAMELIST_EMULATOR_H
#define STORMM_NAMELIST_EMULATOR_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "DataTypes/common_types.h"
#include "FileManagement/file_enumerators.h"
#include "Reporting/error_format.h"
#include "Reporting/summary_file.h"
#include "namelist_element.h"
#include "namelist_enumerators.h"

namespace stormm {
namespace namelist {

using constants::CaseSensitivity;
using constants::ExceptionResponse;
using data_types::isFloatingPointScalarType;
using data_types::isSignedIntegralScalarType;
using data_types::isUnsignedIntegralScalarType;
using diskutil::PrintSituation;
using review::default_output_file_width;
  
/// \brief Collection of variables to transcribe information contained within a namelist
class NamelistEmulator {
public:

  /// \brief Construct an object to emulate Fortran namelist functionality, with improvements.
  ///
  /// \param title_in   The title of the namelist
  /// \param casing_in  Case sensitivity to abide (default "AUTOMATIC", which in this context means
  ///                   that namelist titles and keywords are case insensitive but values are
  ///                   case sensitive)
  NamelistEmulator(const std::string &title_in,
                   CaseSensitivity casing_in = CaseSensitivity::AUTOMATIC,
		   ExceptionResponse unknown_keyword_policy = ExceptionResponse::WARN,
                   const std::string &help_in = std::string("No description provided"),
                   bool cli_content_in = false);

  /// \brief Obtain the title of this namelist (i.e. &cntrl or &dock)
  const std::string& getTitle() const;

  /// \brief Detect whether the namelist emulator serves command-line inputs for a program.
  bool isCommandLineContent() const;
  
  /// \brief Obtain the number of parameters catalogged within this namelist emulator.
  int getKeywordCount() const;

  /// \brief Obtain the case sensitivity setting for this namelist.
  CaseSensitivity getCaseSensitivity() const;

  /// \brief Relay the exception handling policy for this namelist.
  ExceptionResponse getPolicy() const;

  /// \brief Get a keyword from this namelist based on an index.  This is for retrieving the
  ///        keyword itself, not a value associated with a keyword.
  ///
  /// \param index  Index of the keyword in the list held by this namelist.  In fact, this is
  ///               best used to step through the list of keywords in the order they were added,
  ///               not much more.
  const std::string& getKeyword(size_t index) const;

  /// \brief Get the type of a specific keyword within this namelist.
  ///
  /// \param keyword_query  The keyword of interest
  NamelistType getKeywordKind(const std::string &keyword_query) const;

  /// \brief Get the number of entries associated with a specific keyword.
  ///
  /// \param keyword_query  The keyword of interest
  int getKeywordEntries(const std::string &keyword_query) const;

  /// \brief Get the template size of a keyword, the number of sub-keys it contains.  For
  ///        non-struct keywords the template size is reported as 0.
  int getSubKeyCount(const std::string &keyword_query) const;
  
  /// \brief Test whether a keyword has been set, be that by default or user input.
  ///
  /// Overloaded:
  ///   - Get the status of the first instance of a keyword
  ///   - Get the status of a sub-key within a keyword
  ///   - Get the status of a specific, repeated specification of a keyword
  ///
  /// \param keyword_query  The keyword of interest
  /// \param sub_key        The keyword of interest
  /// \param repeat_no      The number of the repetition to check for its status.  This is needed
  ///                       only in the case of STRUCT-type keywords, as a plain INTEGER, REAL,
  ///                       or STRING keyword will be "established" if there is but one default
  ///                       or user-specified value given--subsequent applications of such a
  ///                       keyword are, by construction, established (status cannot be missing).
  ///                       In contrast, a STRUCT keyword can be specified many times, but the
  ///                       status of its individual members may still be inquestion if not all
  ///                       components of the STRUCT are given in each application.
  /// \{
  InputStatus getKeywordStatus(const std::string &keyword_query) const;
  InputStatus getKeywordStatus(const std::string &keyword_query, const std::string &sub_key,
                               int repeat_no = 0) const;
  /// \}

  /// \brief Test whether a namelist contains a particular keyword at all.
  ///
  /// Overloaded:
  ///   - Test for a keyword by name only
  ///   - Test for a named keyword of a specific type
  ///
  /// \param query       The keyword to search for
  /// \param query_kind  The type that the keyword must have if it is present
  /// \{
  bool hasKeyword(const std::string &query) const;
  bool hasKeyword(const std::string &query, NamelistType query_kind) const;
  /// \}
  
  /// \brief Get the value of a boolean keyword from the within the namelist.
  ///
  /// Overloaded:
  ///   - Get a BOOL value for a non-STRUCT keyword
  ///   - Get a BOOL value from within a STRUCT keyword
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  /// \param sub_key        Identifier for the member variable within the STRUCT of interest
  /// \param index          For keywords that store multiple values, retrieve this value
  /// \{
  bool getBoolValue(const std::string &keyword_query) const;
  bool getBoolValue(const std::string &keyword_query, const std::string &sub_key,
                   int index = 0) const;
  /// \}
  
  /// \brief Get a labeled integer value from within the namelist.
  ///
  /// Overloaded:
  ///   - Get an INTEGER value for a non-STRUCT keyword
  ///   - Get an INTEGER value from within a STRUCT keyword
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  /// \param sub_key        Identifier for the member variable within the STRUCT of interest
  /// \param index          For keywords that store multiple values, retrieve this value
  /// \{
  int getIntValue(const std::string &keyword_query, int index = 0) const;
  int getIntValue(const std::string &keyword_query, const std::string &sub_key,
                  int index = 0) const;
  /// \}

  /// \brief Get a labeled real number value from within the namelist.
  ///
  /// Overloaded:
  ///   - Get a REAL value associated with a non-STRUCT keyword
  ///   - Get a REAL value from within a STRUCT keyword
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  /// \param sub_key        Identifier for the member variable within the STRUCT of interest
  /// \param index          For keywords that store multiple values, retrieve this value
  /// \{
  double getRealValue(const std::string &keyword_query, int index = 0) const;
  double getRealValue(const std::string &keyword_query, const std::string &sub_key,
                      int index = 0) const;
  /// \}

  /// \brief Get a labeled string value from within the namelist.
  ///
  /// Overloaded:
  ///   - Get a STRING value associated with a non-STRUCT keyword
  ///   - Get a STRING value from within a STRUCT keyword
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  /// \param sub_key        Identifier for the member variable within the STRUCT of interest
  /// \param index          For keywords that store multiple values, retrieve this value
  /// \{
  const std::string& getStringValue(const std::string &keyword_query, int index = 0) const;
  const std::string& getStringValue(const std::string &keyword_query, const std::string &sub_key,
                                    int index = 0) const;
  /// \}

  /// \brief Get all boolen values assigned to a particular keyword.  The keyword must have STRUCT
  ///        type, as the only other option for a BOOL keyword is a single value.
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  std::vector<bool> getAllBoolValues(const std::string &keyword_query,
                                    const std::string &sub_key) const;
  
  /// \brief Get all integer values assigned to a particular keyword
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  std::vector<int> getAllIntValues(const std::string &keyword_query,
                                   const std::string &sub_key = std::string("")) const;

  /// \brief Get all real values assigned to a particular keyword
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  std::vector<double> getAllRealValues(const std::string &keyword_query,
                                       const std::string &sub_key = std::string("")) const;

  /// \brief Get all string values assigned to a particular keyword
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  std::vector<std::string>
  getAllStringValues(const std::string &keyword_query,
                     const std::string &sub_key = std::string("")) const;

  /// \brief Report the help message associated with a keyword or sub-key.  This can be useful for
  ///        developers who wish to alert users to erroneous input.
  ///
  /// Overloaded:
  ///   - Report user documentation for the namelist as a whole
  ///   - Report user documentation for any keyword (including STRUCTs) within the namelist
  ///   - Report user documentation for a sub-key within a STRUCT keyword in the namelist
  ///
  /// \param keyword_query  Identifier of the keyword of interest
  /// \param sub_key        Identifier for the member variable within the STRUCT of interest
  const std::string& getHelp() const;
  const std::string& getHelp(const std::string &keyword_query) const;
  const std::string& getHelp(const std::string &keyword_query, const std::string &sub_key) const;
  /// \}

  /// \brief Get a pointer to the object itself.
  const NamelistEmulator* getSelfPointer() const;
  
  /// \brief Assign a value, external to the object, based on the content inside of it.  This will
  ///        first check whether the appropriate keyword is not missing (that is has a default
  ///        value, or has been specified by the user).  If the associated keyword is indeed
  ///        missing, there will be no effect on the external variable.
  ///
  /// Overloaded:
  ///   - Assign other integers
  ///   - Assign other real values
  ///   - Assign other strings
  ///   - Assign triplets of other numerical variables (useful for cases where a generic keyword
  ///     can assign settings for all three dimensions of an object)
  ///   - Provide a multiplier to scale the units of user input into internal units (a common case
  ///     is degrees to radians)
  ///
  /// \param var            The integer, real, or string ariable to assign
  /// \param mult           Multiplication factor to apply to any value extracted from the
  ///                       &namelist.  This factor is only available for scalar results (or
  ///                       triplicate input extractions) and will be cast to the data type of var
  ///                       for integral types.
  /// \param keyword_query  The keyword associated with the input data of interest
  /// \param sub_key        The sub-key within a STRUCT associated with the input data of interest
  /// \param index          Index of the keyword repeat to retrieve
  /// \{
  template <typename T>
  void assignVariable(T *var, double mult, const std::string &keyword_query, int index = 0) const;

  template <typename T>
  void assignVariable(T *var, const std::string &keyword_query, int index = 0) const;

  void assignVariable(std::string *var, double mult, const std::string &keyword_query,
                      int index = 0) const;

  void assignVariable(std::string *var, const std::string &keyword_query, int index = 0) const;
  
  template <typename T>
  void assignVariable(T *var_x, T *var_y, T *var_z, double mult, const std::string &keyword_query,
                      int index = 0) const;

  template <typename T>
  void assignVariable(T *var_x, T *var_y, T *var_z, const std::string &keyword_query,
                      int index = 0) const;

  template <typename T>
  void assignVariable(T *var, const std::string &keyword_query, const std::string &sub_key,
                      int index = 0) const;

  template <typename T>
  void assignVariable(T *var, double mult, const std::string &keyword_query,
                      const std::string &sub_key, int index = 0) const;

  void assignVariable(std::string *var, double mult, const std::string &keyword_query,
                      const std::string &sub_key, int index = 0) const;

  void assignVariable(std::string *var, const std::string &keyword_query,
                      const std::string &sub_key, int index = 0) const;

  template <typename T>
  void assignVariable(T *var_x, T *var_y, T *var_z, double mult, const std::string &keyword_query,
                      const std::string &sub_key, int index = 0) const;

  template <typename T>
  void assignVariable(T *var_x, T *var_y, T *var_z, const std::string &keyword_query,
                      const std::string &sub_key, int index = 0) const;
  /// \}

  /// \brief Set the title of the namelist.
  ///
  /// \param title_in  The title to set
  void setTitle(const std::string &title_in);

  /// \brief Set whether the namelist actually serves command line input for a whole program.
  ///
  /// \param cli_content_in  Set to TRUE to indicate that the namelist does command line parsing
  void setCommandLineContent(bool cli_content_in = true);
  
  /// \brief Add a keyword to the namelist.
  ///
  /// Overloaded:
  ///   - Add a single keyword (be it a INTEGER, REAL, STRING, or STRUCT namelist element)
  ///   - Add multiple keywords
  ///   - Provide a NamelistElement object
  ///   - Provide input parameters with one-to-one correspondence to NameListElement objects (this
  ///     can save space and make the API cleaner)
  ///   - Import a named keyword from another NamelistEmulator object (the other object will be
  ///     checked to see that it has the keyword of interest)
  ///
  /// \param new_keys  The keywords to add
  /// \param new_key   The keyword to add
  /// \param other     Another NamelistEmulator from which to import a keyword
  /// \param query     Name of the keyword to copy from another NamelistEmulator
  /// \{
  void addKeyword(const std::vector<NamelistElement> &new_keys);

  void addKeywords(const std::vector<NamelistElement> &new_keys);

  void addKeyword(const NamelistElement &new_key);

  void addKeyword(const std::string &keyword_in, NamelistType kind_in,
                  const std::string &default_in = std::string(""),
                  DefaultIsObligatory obligate = DefaultIsObligatory::NO,
                  InputRepeats rep_in = InputRepeats::NO,
                  const std::string &help_in = std::string("No description provided"));

  void addKeyword(const std::string keyword_in, const std::vector<std::string> &sub_keys_in,
                  const std::vector<NamelistType> &sub_kinds_in,
                  const std::vector<std::string> &default_list,
                  DefaultIsObligatory obligate_list = DefaultIsObligatory::NO,
                  InputRepeats rep_in = InputRepeats::NO,
                  const std::string &help_in = std::string("No description provided"),
                  const std::vector<std::string> &sub_help_in =
                  std::vector<std::string>(1, "No description provided"),
                  const std::vector<KeyRequirement> &template_requirements_in = {});

  void addKeyword(const NamelistEmulator *other, const std::string &query);

  void addKeyword(const NamelistEmulator &other, const std::string &query);
  /// \}

  /// \brief Set a default value for one of the namelist's keywords.  Like addDefaultValue() and
  ///        other functions below, the effect will be to pass through the NamelistEmulator to
  //         call the eponymous member function in an individual NamelistElement.
  ///
  /// \param key                The name of the keyword within the namelist
  /// \param modified_default   The new default setting to apply to the keyword
  /// \param default_idx        The index of the default to set, in the event that there are
  ///                           already multiple default values
  /// \param modified_defaults  The list of modified default values
  /// \param sub_key_specs      The list of sub-keys to which each default corresponds
  /// \{
  void setDefaultValue(const std::string &key, const std::string &modified_default,
                       int default_idx = 0);
  void setDefaultValue(const std::string &key, const std::vector<std::string> &modified_defaults,
                       const std::vector<std::string> &sub_key_specs);
  /// \}
  
  /// \brief Add a value to one keyword's default settings.  This enables a single keyword to have
  ///        a collection of default values.  The keyword will be checked to ensure that it permits
  ///        multiple values.
  ///
  /// \param key           The keyword of interest
  /// \param next_default  The new value to include, as a string, to be re-interpreted as necessary
  void addDefaultValue(const std::string &key, const std::string &next_default);

  /// \brief Activate a BOOL-type keyword, or a BOOL-type member of a STRUCT-type keyword.  These
  ///        functions are analogous to assignElement() below, but because BOOL-type keywords do
  ///        not take distinct values the functions that set BOOL-type variables to ON are named
  ///        differently.
  ///
  /// Overloaded:
  ///   - Activate the named BOOL keyword
  ///   - Activate a BOOL sub-key within the named STRUCT
  ///
  /// \param key           The keyword of interest
  /// \param sub_key       The sub-key of interest
  /// \{
  void activateBool(const std::string &key);
  void activateBool(const std::string &key, const std::string &sub_key);
  /// \}

  /// \brief Assign values to elements of each particular NamelistType.  These overloaded functions
  ///        can be called from anywhere, but constructors making control objects for programs
  ///        using the STORMM libraries are the ideal place to use them.  They in turn call the
  ///        set(...)Value member functions of the target NamelistElement object in the
  ///        NamelistEmulator.  Returns 1 if the given value was successfully assigned to the
  ///        label or 0 if not.
  ///
  /// Overloaded:
  ///   - Assign a single integer, real, or string value to an INTEGER, REAL, or STRING namelist
  ///     element, respectively
  ///   - Assign a single integer, real, or string value to the INTEGER, REAL, or STRING member of
  ///     a STRUCT namelist element
  ///
  /// \param key           The keyword of interest
  /// \param sub_key       The sub-key of interest
  /// \param sub_key       The value to assign to the keyword or STRUCT sub-key
  /// \{
  int assignElement(const std::string &key, const std::string &value);
  int assignElement(const std::string &key, const std::string &sub_key, const std::string &value);
  /// \}

  /// \brief When loading data for STRUCT-type keywords, the decision to increment the number of
  ///        entries on file cannot be made with the first sub-key assignment.  Instead, the entire
  ///        struct must be read from input before the number of entries can be incremented.
  ///        Because a function (readNamelist(), see input.h) that manages a NamelistEmulator
  ///        loops over the input that can be associated with each list of subk-eys for a given
  ///        STRUCT, that function must go through the NamelistEmulator in order to increment
  ///        the keyword's entry count.
  ///
  /// \param key  The STRUCT-type keyword of interest (its membership in the namelist will be
  ///             verified)
  void triggerResizeBuffer(const std::string &key);

  /// \brief Attach a help message to the namelist itself, to a keyword within the namelist, or
  ///        even to a member variable of a STRUCT-associated keyword in the namelist.  This is
  ///        provided so that developers do not have to include help messages at the initialization
  ///        of each namelist keyword.  This will overwrite existing help messages.
  ///
  /// Overloaded:
  ///   - Attach a help message to the namelist as a whole
  ///   - Attach a help message to any keyword (including STRUCTs) within the namelist
  ///   - Attach a help message to a sub-key within a STRUCT keyword in the namelist
  ///
  /// \param blurb      The help message to attach
  /// \param key        Label of the namelist keyword to which the help message gets attached
  /// \param sub_key    Label of the member variable of the STRUCT namelist keyword to which the
  ///                   help message gets attached
  /// \{
  void addHelp(const std::string &blurb);
  void addHelp(const std::string &key, const std::string &blurb);
  void addHelp(const std::string &key, const std::string &sub_key, const std::string &blurb);
  /// \}

  /// \brief Add a category to a namelist to group its keywords for user documentation
  ///
  /// \param new_category  The name of the new keyword category
  void addCategory(const std::string &new_category);

  /// \brief Place a namelist keyword into one of a list of arbitrary categories defined by the
  ///        developer.  This is for organizing the user documentation.
  ///
  /// \param key             The namelist keyword to find and categorize
  /// \param category_label  The category to put it in
  void categorizeKeyword(const std::string &key, const std::string &category_label);

  /// \brief Change the requirement associated with a keyword.  By default, all keywords are set to
  ///        "REQUIRED" just as all STRUCT-type keyword subkeys are required unless stated
  ///        otherwise.  Even if required, many keywords, like STRUCT subkeys, will have default
  ///        values that satisfy the requirements.
  ///
  /// Overloaded:
  ///   - Set the criticality of a keyword
  ///   - Set the criticality of a subkey within a STRUCT keyword
  ///
  /// \param key            The namelist keyword to find and alter
  /// \param sub_key_query  Name of the sub-key to search if the namelist element is a STRUCT
  /// \param req            The requirement level to impose
  /// \{
  void setImperative(const std::string &key, KeyRequirement req);
  void setImperative(const std::vector<std::string> &directives);
  void setImperative(const std::string &key, const std::vector<std::string> &directives);
  /// \}

  /// \brief Print the documentation for a specific keyword.  The format is fixed in the sense
  ///        that it will have a set indentation, a dash for a bullet point, and the keyword
  ///        printed in a space large enough for a series of related keywords in a column.
  ///
  /// \param p_idx          Index of the keyword within the namelist emulator
  /// \param name_width     Width at which to print keyword names
  /// \param kw_kind_width  Width at which to print keyword kinds
  /// \param kw_dflt        Default value of the parameter, pre-converted to a string and
  ///                       pre-pended with white space for alignment
  void printKeywordDocumentation(int p_idx, int name_width, int kw_kind_width,
                                 const std::string &kw_dflt) const;

  /// \brief Convert the default value of a keyword to a string for output in a formatted table.
  ///
  /// \param tkw  The keyword of interest
  std::string convertDefaultToString(const NamelistElement &tkw) const;

  /// \brief Print a detailed message concerning the user documentation for keywords in this
  ///        namelist.  This function can be called from the main program, for any namelists that
  ///        accept its inputs, and is itself called by the printProgramDocumentation() function
  ///        in the docs namespace (see Reporting/custom_help.h)
  void printHelp() const;

  /// \brief Print a complete table of the values for all parameters in this namelist, starting
  ///        including their sources (input statuses, i.e. DEFAULT, MISSING, or USER-SPECIFIED).
  ///
  /// Overloaded:
  ///   - Write to a string for further processing
  ///   - Write directly to an output file stream
  ///   - Write to a named output file, given a state that it must be found in
  ///
  /// \param foutp              File to which information should be printed, defaulting to the
  ///                           terminal
  /// \param file_width         The width at which to print results (this parameter will be ignored
  ///                           and replaced with the terminal width if the output goes to the
  ///                           terminal)
  /// \param max_entry_reports  The maximum number of entries to report from any given keyword
  /// \param print_decor        Indicate whether to introduce the namelist or rely on some previous
  ///                           iteration writing to the same file
  /// \{
  std::string printContents(int file_width = default_output_file_width, int max_entry_counts = 4,
                            NamelistIntroduction print_decor = NamelistIntroduction::HEADER) const;
  
  void printContents(std::ostream *foutp, int file_width = default_output_file_width,
                     int max_entry_counts = 4,
                     NamelistIntroduction print_decor = NamelistIntroduction::HEADER) const;

  void printContents(std::ofstream *foutp, int file_width = default_output_file_width,
                     int max_entry_counts = 4,
                     NamelistIntroduction print_decor = NamelistIntroduction::HEADER) const;

  void printContents(const std::string &file_name, PrintSituation expectation,
                     int file_width = default_output_file_width, int max_entry_counts = 4,
                     NamelistIntroduction print_decor = NamelistIntroduction::HEADER) const;
  /// \}

private:
  bool cli_content;                        ///< Flag to indicate that the namelist actually
                                           ///<   contains information from a command-line
                                           ///<   interface
  std::string title;                       ///< Title of this namelist, i.e. &cntrl
  std::vector<NamelistElement> keywords;   ///< List of all keywords stored in this namelist
  CaseSensitivity casing;                  ///< Case sensitivity of the namelist
  ExceptionResponse policy;                ///< Reaction to unusual inputs: more than one instance
                                           ///<   of a non-repeating keyword, an unknown keyword,
                                           ///<   an invalid value
  std::string help_message;                ///< Help message for this namelist, as a whole
  std::vector<std::string> category_names; ///< Names of categories for keywords / labels in this
                                           ///<   namelist (this is for end-user documentation)

  /// Categories containing lists of labels (for sorting the end-user documentation output)
  std::vector<std::vector<std::string>> categories;

  /// \brief Find the index of a keyword based on its label.  This encapsulates throwing an
  ///        exception if the label cannot be found.  While this is similar to findStringInVector()
  ///        in the nearby parse library (see parse.h), the necessity of making some response if
  ///        the label cannot be found and the unusual shape of the vector elements makes it
  ///        reasonable to keep a separate function for this.
  ///
  /// \param query   The label to seek inside the namelist
  int findIndexByKeyword(const std::string &query) const;

  /// \brief Verify the establishment of a keyword, either by user specification or default.
  ///        This is called before returning any information collected from a namelist and will
  ///        trigger an error if the requested information does not exist.
  ///
  /// \param keyword_query  The label that was sought inside the namelist
  /// \param p_index        The keyword (parameter) index within the larger namelist emulator
  /// \param caller         Calling function that is seeking some keyword's information, i.e.
  ///                       getIntValue()
  void verifyEstablishment(const std::string &keyword_query, const size_t p_index,
                           const char* caller) const;
};

} // namespace namelist
} // namespace stormm

#include "namelist_emulator.tpp"

#endif
