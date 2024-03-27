// -*-c++-*-
#ifndef STORMM_MOLOBJ_DATAITEM_H
#define STORMM_MOLOBJ_DATAITEM_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "DataTypes/stormm_vector_types.h"
#include "Parsing/textfile.h"
#include "Potential/energy_enumerators.h"
#include "molecule_format_enumerators.h"
#include "mdlmol_request.h"

namespace stormm {
namespace structure {

using constants::ExceptionResponse;
using constants::ModificationPolicy;
using parse::TextFile;
  
/// \brief Store a data item from within an SD file.  Data items begin with a line of the form
///        "> <ITEM_NAME>" (additional specifications are possible), end with a single blank line,
///        and can contain customized format information.  The information will be read as a series
///        of raw strings, one per line and each representing a whole line.  Once read, member
///        functions of the class can extract specific integer, real, char4, or string coded data
///        from within the item.
class MdlMolDataItem {
public:

  /// \brief The constructor takes the original text file and the number of the first line of the
  ///        data item, or the information pertinent to the header line and a series of strings
  ///        that will form the body of data in the item.  Any information submitted as part of a
  ///        string such as an atom or bond index, which could then go to an output file, must be
  ///        in the Fortran array indexing (starting at 1).  Integer information read from the
  ///        data body can be adjusted to meet internal array indexing.
  ///
  /// \param item_name_in        Name of the data item (an item_name of "ABCD" will appear as
  ///                            "<ABCD>" on the data item's first line)
  /// \param external_regno_in   An identification number for the compound referencing an external
  ///                            database.  An external identification number of "PD-3647" will be
  ///                            represented as "(PD-3647)" on the data item's header line.
  /// \param maccs_ii_number_in  An integer representing the corresponding field in a MACCS-II
  ///                            database.  A field number such as "6519" will be represented as
  ///                            "DT6519" in the data item header line.
  /// \param header_info         A bit-packed unsigned integer containing, in its low to high bits,
  ///                            whether to display the internal identification number (based on
  ///                            the order in the SD file), the external identification number,
  ///                            the item name, the field number, and finally a "FROM ARCHIVES"
  ///                            declaration.
  /// \param tf                  The original text of the SD file, committed to RAM
  /// \param line_number         Line of the file at which to begin reading the data item
  /// \param line_advance        Line of the file that the construction of the data item leads to
  ///                            (after taking in lines of the item)
  /// \param compound_line_end   Index of the last line of the compound within the SD file (the
  ///                            line contains $$$$)
  /// \param title               The title of the structure, if known, for error tracing purposes
  /// \param mpol                Policy about modifying details of the data items that may not
  ///                            conform to the Biovia standard
  /// \param notify              Policy on whether to notify the user if modifications to certain
  ///                            details of the data item take place
  /// \{
  MdlMolDataItem(const std::string &item_name_in = std::string(""),
                 const std::string &external_regno_in = std::string(""),
                 int internal_regno_in = -1, int maccs_ii_number_in = -1, uint header_info = 0U,
                 const std::vector<std::string> &body_in = {},
                 ModificationPolicy mpol = ModificationPolicy::DO_NOT_MODIFY,
                 ExceptionResponse notify = ExceptionResponse::WARN);

  MdlMolDataItem(const TextFile &tf, int line_number, int *line_advance,
                 int compound_line_end = -1, const std::string &title = std::string(""),
                 ModificationPolicy mpol = ModificationPolicy::DO_NOT_MODIFY,
                 ExceptionResponse notify = ExceptionResponse::WARN);

  MdlMolDataItem(const MdlMolDataRequest &ask, const std::vector<std::string> &body_in,
                 ModificationPolicy mpol = ModificationPolicy::DO_NOT_MODIFY,
                 ExceptionResponse notify = ExceptionResponse::WARN);
  /// \}

  /// \brief The default copy and move constructors, as well as copy and move assignment operators,
  ///        are applicable to this object which has no const members or pointers to repair.
  /// \{
  MdlMolDataItem(const MdlMolDataItem &original) = default;
  MdlMolDataItem(MdlMolDataItem &&original) = default;
  MdlMolDataItem& operator=(const MdlMolDataItem &other) = default;
  MdlMolDataItem& operator=(MdlMolDataItem &&other) = default;
  /// \}

  /// \brief Indicate whether to display the internal identification number (based on the position
  ///        of the molecule in the SD file) in the header line.
  bool placeInternalRegnoInHeader() const;

  /// \brief Indicate whether to display the (user-supplied) external identification number in the
  ///        header line.
  bool placeExternalRegnoInHeader() const;
  
  /// \brief Indicate whether to display the data item name (title) in the header line, between
  ///        angular brackets.
  bool placeTitleInHeader() const;

  /// \brief Indicate whether to display the MACCS-II field number in the header.
  bool placeMaccsIIFieldInHeader() const;

  /// \brief Indicate that the "FROM ARCHIVES" tag should go in the header line.
  bool noteArchivesInHeader() const;
  
  /// \brief Get a const reference to the item name, if it exists.  If there is no item name, a
  ///        const reference to a blank string will be returned.
  const std::string& getItemName() const;

  /// \brief Get a const reference to the output item name, if it exists.  If there is no item
  ///        name, a const reference to a blank string will be returned.
  const std::string& getOutputItemName() const;

  /// \brief Get the external registry number.
  const std::string& getExternalRegistryNumber() const;

  /// \brief Get the internal registry number.  This will return -1 if no such registry number
  ///        exists.
  int getInternalRegistryNumber() const;

  /// \brief Get the MACCS-II database field number.
  int getMaccsFieldNumber() const;

  /// \brief Get the state variable enumeration for this data item.  This is used in the case of a
  ///        STATE_VARIABLE kind of data item, to signal that a particular component of the
  ///        molecular mechanics energy, or some intensive or expansive system property, should be
  ///        retrieved.
  StateVariable getTrackedState() const;

  /// \brief Get the data line count.
  int getDataLineCount() const;
  
  /// \brief Get a particular data lines from the item body.
  ///
  /// \param line_index  The data line of interest
  const std::string& getDataLine(const int line_index) const;
  
  /// \brief Get all data lines.
  const std::vector<std::string>& getBody() const;

  /// \brief Get the custom data item type.  Here, the "NATIVE" code indicates that no custom
  ///        enumeration is present.
  MdlMolDataItemKind getKind() const;
  
  /// \brief Retrieve a string from the data lines of a data item, assuming that individual words
  ///        on each line are separated by one or more white space characters and that the
  ///        quotation marks "" and '' collect everything between them into a single word.
  ///
  /// Overloaded:
  ///   - Extract a free-format string based on words separated by white space, assuming that
  ///     individual words on each line are separated by one or more white space characters and
  ///     that the quotation marks "" and '' collect everything between them into a single word.
  ///   - Extract a column-formatted string based on strict character indices
  ///
  /// \param element_number  Number of the word on the line (words are separated by white space)
  /// \param start_pos       The starting position within the line at which to begin reading
  /// \param start_pos       Length of the fixed-column reading to perform
  /// \param line_number     Number of the data line on which to find the string
  /// \{
  std::string parseString(int element_number, int line_number) const;
  std::string parseString(int start_pos, int length, int line_number) const;
  /// \}

  /// \brief Retrieve a signed integer from the data lines of a data item.  The assumptions of
  ///        parseString() above apply, as do the overloads and descriptions of formal arguments.
  /// \{
  llint parseInteger(int element_number, int line_number) const;
  llint parseInteger(int start_pos, int length, int line_number) const;
  /// \}

  /// \brief Retrieve an unsigned integer from the data lines of a data item.  The assumptions of
  ///        parseString() above apply, as do the overloads and descriptions of formal arguments.
  /// \{
  ullint parseUnsigned(int element_number, int line_number) const;
  ullint parseUnsigned(int start_pos, int length, int line_number) const;
  /// \}

  /// \brief Retrieve a real number from the data lines of a data item.  The assumptions of
  ///        parseString() above apply, as do the overloads and descriptions of formal arguments.
  /// \{
  double parseReal(int element_number, int line_number) const;
  double parseReal(int start_pos, int length, int line_number) const;
  /// \}

  /// \brief Retrieve a tuple of four characters from the data lines of a data item.  Assumptions
  ///        from parseString() above apply, as do overloads and descriptions of formal arguments.
  /// \{
  char4 parseChar4(int element_number, int line_number) const;
  char4 parseChar4(int start_pos, int length, int line_number) const;
  /// \}

  /// \brief Match this data item with a series of identification tags.
  ///
  /// \param item_name_comp  Item name for comparison
  /// \param ext_regno_comp  External registry number for comparison
  /// \param maccs_ii_no_comp  MACCS-II database field number for comparison (omit leading "DT")
  /// \{
  bool matchItemName(const std::string &item_name_comp, const std::string &ext_regno_comp,
                     int maccs_ii_no_comp = -1) const;

  bool matchItemName(const std::string &item_name_comp, int maccs_ii_no_comp = -1) const;

  bool matchRegistryNumber(const std::string &ext_regno_comp, int maccs_ii_no_comp = -1) const;

  bool matchMaccsField(int maccs_ii_no_comp) const;
  /// \}

  /// \brief Set the item name and apply checks to the result.
  ///
  /// \param item_name_in  The item name to assign
  /// \param mpol          Indicate what to do if the proposed name does not meet strict Biovia
  ///                      SD file standards but is otherwise salvageable.
  /// \param notify        Indicate what to do if modifications are made to the proposed name.
  void setItemName(const std::string &item_name_in,
                   ModificationPolicy mpol = ModificationPolicy::DO_NOT_MODIFY,
                   ExceptionResponse notify = ExceptionResponse::WARN);

  /// \brief Add a line to the data item.
  ///
  /// \param text  The text of the line to add.  This should not include carriage returns, but if
  ///              it does then it will all be taken as one "data line" to be printed.
  void addDataLine(const std::string &text);
  
private:
  MdlMolDataItemKind kind;       ///< Classifies the data item.  All data items found in the
                                 ///<   original SD file are considered "NATIVE", whereas each
                                 ///<   custom data item falls under one of the other
                                 ///<   classifications.
  StateVariable tracked_state;   ///< The type of energy term to track if the kind is
                                 ///<   STATE_VARIABLE.  If the kind is anything else, this member
                                 ///<   variable is set to ALL_STATES and ignored.
  std::string item_name;         ///< Name of the data item (optional, but an effective means of
                                 ///<   distiction)
  std::string output_item_name;  ///< Name of the data item (optional, but an effective means of
                                 ///<   distiction)
  std::string external_regno;    ///< External registry number
  int internal_regno;            ///< Internal registry number, based on the structure index within
                                 ///<   the SD file and a pure integer
  int maccs_ii_number;           ///< MACCS-II database field number
  bool use_internal_regno;       ///< Flag to have the header line use the SD file's internal
                                 ///<   structure numbering in the data item header line
  bool use_external_regno;       ///< Flag to have the header line use the external registry number
  bool use_item_name;            ///< Flag to use the item name in the header between < > marks
  bool use_maccs_ii_number;      ///< Flag to use the MACCS-II field number in the header after
                                 ///<   the prefix "DT"
  bool note_archives;            ///< Flag to have the header line note "FROM ARCHIVES"

  /// Data lines from the item, one per array element
  std::vector<std::string> body;
  
  /// \brief Validate the choice of item name.
  void validateItemName(ModificationPolicy mpol = ModificationPolicy::DO_NOT_MODIFY,
                        ExceptionResponse notify = ExceptionResponse::WARN);
};

/// \brief Get the unsigned integer bit-packed code expressing various aspects of the data item
///        header line.
///
/// Overloaded:
///   - Accept boolean arguments for all five aspects
///   - Accept a data item request and glean the details from that
///
/// \param use_internal_regno   Flag to display the internal registry number
/// \param use_external_regno   Flag to display the external registry number (in parentheses)
/// \param use_item_name        Flag to display the item name (in angular brackets)
/// \param use_maccs_ii_field   Flag to display the MACCS-II database field number (after "DT")
/// \param state_from_archives  Flag to display "FROM ARCHIVES"
/// \param ask                  The data item request
/// \{
uint getDataItemHeaderCode(bool use_internal_regno = false, bool use_external_regno = false,
                           bool use_item_name = true, bool use_maccs_ii_field = false,
                           bool state_from_archives = false);

uint getDataItemHeaderCode(const MdlMolDataRequest &ask);
/// \}

} // namespace structure
} // namespace stormm

#endif
