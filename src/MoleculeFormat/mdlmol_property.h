// -*-c++-*-
#ifndef STORMM_MOLOBJ_PROPERTY_H
#define STORMM_MOLOBJ_PROPERTY_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Parsing/textfile.h"
#include "molecule_format_enumerators.h"

namespace stormm {
namespace structure {

using parse::TextFile;
  
/// \brief A molecular or atomic property read from an MDL .mol or SDF file
class MdlMolProperty {
public:

  /// \brief The constructor takes all member variable inputs, or the original text and a line
  ///        number within it.
  ///
  /// \param tf            Text of the original .sdf or .mol file, read into RAM
  /// \param line_number   Number of the line on which to read the data
  /// \param line_advance  The line advancement position, which will account for data lines in
  ///                      some properties
  /// \param title         The title of the structure, if known, for error tracing purposes
  /// \{
  MdlMolProperty(const char4 code_in = { ' ', ' ', ' ', ' ' }, int substrate_in = -1,
                 int entry_count_in = 0, int entry_depth_in = 0, bool exclusions_in = false,
                 int substrate_line_pos_in = -1, int entry_count_line_pos_in = -1,
                 int entry_read_start_pos_in = 0,
                 const std::vector<int> &entry_data_bounds_in = {},
                 const std::vector<MolObjPropField> &entry_detail_in = {},
                 const std::vector<MolObjIndexKind> &entry_adjustment_in = {},
                 const std::vector<int> &int_data_in = {},
                 const std::vector<double> &real_data_in = {},
                 const std::vector<std::string> &str_data_in = {},
                 const std::vector<std::string> &data_lines_in = {});

  MdlMolProperty(const TextFile &tf, int line_number, int *line_advance,
                 const std::string &title = std::string(""));
  /// \}

  /// \brief The default copy and move constructors, as well as copy and move assignment operators,
  ///        are applicable to this object which has no const members or pointers to repair.
  /// \{
  MdlMolProperty(const MdlMolProperty &original) = default;
  MdlMolProperty(MdlMolProperty &&original) = default;
  MdlMolProperty& operator=(const MdlMolProperty &other) = default;
  MdlMolProperty& operator=(MdlMolProperty &&other) = default;
  /// \}
  
  /// \brief Get the kind of property according to the internal enumerator.
  MdlMolPropertyKind getKind() const;

  /// \brief Get the property code.  This will indicate whether to obliterate certain types of
  ///        information from the atom block of the V2000 MDL MOL format entry.
  char4 getCode() const;

  /// \brief Get the substrate atom or S-group.  The nature of the property will indicate which
  ///        type of substrate this is.
  int getSubstrate() const;

  /// \brief Get the substrate atom or S-group for the purposes of printing.  The nature of the
  ///        property will indicate which type of substrate this is, but all substrates' indices
  ///        shift by -1 when being read and therefore shift by +1 when being written to the MDL
  ///        MOL format.
  int getPrintedSubstrate() const;

  /// \brief Indicate whether the atom list encoded by this property is based on excluded elements
  ///        (return TRUE) or included ones (return FALSE).  If the property is not an "M  ALS"
  ///        entry, raise a runtime error.
  bool applyToExclusions() const;

  /// \brief Get a code ('T' if TRUE, 'F' if FALSE) to indicate whether the atom list describes
  ///        exclusions.
  char getExclusionCode() const;
  
  /// \brief Get the number of entries for this property.
  int getEntryCount() const;

  /// \brief Get the number of data lines for this property.
  int getDataLineCount() const;
  
  /// \brief Get a value from the property for a specific entry.  Inputs to these functions
  ///        will be checked by checkAttributeValidity().
  ///
  /// \param entry_index      The index of the entry to access
  /// \param attribute_index  The index of the attribute to retrieve
  /// \{
  int getIntegerValue(int entry_index, int attribute_index) const;
  int getPrintedIntegerValue(int entry_index, int attribute_index) const;
  double getRealValue(int entry_index, int attribute_index) const;
  char4 getChar4Value(int entry_index, int attribute_index) const;
  std::string getStringValue(int entry_index, int attribute_index) const;
  /// \}

  /// \brief Get a const pointer to one of the data lines.
  const std::string& getDataLine(int index) const;

  /// \brief Recompose the text string for this property as it would appear in an MDL MOL file.
  std::string getMdlText() const;
  
  /// \brief Define the property code.
  ///
  /// Overloaded:
  ///   - Provide the first letter and the final three
  ///   - Provide all four letters as a pre-made char tuple
  ///
  /// \param x        The first letter of the three-letter tuple in the property code
  /// \param y        The second letter of the three-letter tuple defining the property code
  /// \param z        The last letter of the three-letter tuple defining the property code
  /// \param major    The first letter that will appear on the line if the property is written to
  ///                 a file
  /// \param code_in  A pre-made code to apply, ordered 3-4-5-0 in terms of the indices that each
  ///                 character would occupy on a property line of an MDL MOL format entry
  /// \{
  void setCode(char x, char y, char z, char major = 'M');

  void setCode(char4 code_in);
  /// \}

  /// \brief Define the substrate to be used, whether an atom or an S-group.
  ///
  /// \param index  Atom or S-group index of the substrate.  This should be provided for the
  ///               C / C++ array element, not the file format (+1 will be added if and when the
  ///               information is committed to a file).
  void setSubstrate(int index);

  /// \brief Define the entry format for the property.  If applied to an existing MdlMolProperty
  ///        with nonzero depth, this will create an error.
  ///
  /// \param entry_detail_in  List of data types of the entry elements
  /// \param entry_adjustment_in  List of index adjustments for each element of the entry.  Atom
  ///        and bond elements will have their indices adjusted by -1, but only if the information
  ///        is read from a file--programmatic input is still expected to occur in the C / C++
  ///        array numbering.
  void setEntryFormat(const std::vector<MolObjPropField> &entry_detail_in,
                      const std::vector<MolObjIndexKind> &entry_adjustment_in);

  /// \brief Add an entry to the MDL MOL V2000 format property.  The layout must match that
  ///        established in the detail array.  Adjustments will not be applied to the indexing for
  ///        input originating within the program.  Other methods will adjust the indexing of
  ///        atoms, bonds, and S-groups named in auxiliary user input when making properties out
  ///        of such information.
  ///
  /// \brief int_data_in     All integer data, given in the order that INTEGER MolObjPropField
  ///                        types appear in the detail array.  No char4, real or string data
  ///                        indices are included in this array.
  /// \param char4_data_in   List of char4 data elements, given back-to-back in the order that
  ///                        CHAR4 MolObjPropField types appear in the detail array.
  /// \param real_data_in    List of real data, ordered as the other arrays in this input.
  /// \param string_data_in  List of string data, ordered as the other arrays in this input.
  void addEntry(const std::vector <int> &int_data_in, const std::vector <char4> &char4_data_in,
                const std::vector<double> &real_data_in,
                const std::vector<std::string> &str_data_in);

private:
  char4 code;               ///< A three-letter code indicating what the property is.  The "w"
                            ///<   member stores the first letter on the line, which is usually
                            ///<   but not always 'M'.
  MdlMolPropertyKind kind;  ///< The type of property
  int substrate;            ///< One atom or S-group that is central to all entries, relevant to
                            ///<   some properties
  int entry_count;          ///< Number of entries (some properties have maximum numbers of entries
                            ///<   hard-wired into the format, and thus into the constructor)
  int entry_depth;          ///< The number of fields in each entry
  bool exclusions;          ///< Specific to the ATOM_LIST property enumeration, stores the value
                            ///<   of the exclusions flag for elements named in the list ('T' if
                            ///<   the elements are excluded, 'F' if they are not)
  int substrate_line_pos;   ///< The substrate of this property should be displayed on the MDL MOL
                            ///<   file line at this position.  A negative value indicates that the
                            ///<   substrate is not present on the line.  The substrate is always
                            ///<   displayed in %3d format.
  int entry_count_line_pos; ///< The entry count of this property should be displayed on the MDL
                            ///<   MOL file line at this position.  A negative value indicates that
                            ///<   the entry count is not displayed.
  int entry_read_start_pos; ///< The position on the text line at which entries of the property
                            ///<   begin to appear.

  /// Pattern of character lengths for accessing entries of the property, beginning at the
  /// entry_read_start_pos position.
  std::vector<int> entry_data_bounds;
  
  /// Nature of each field in each entry.  The most important classifications are INTEGER and
  ///   CHAR4, although some properties contain longer strings or real numbers.
  std::vector<MolObjPropField> entry_detail;

  /// Indications of adjustments that must be made to convert from the MDL MOL format conventions
  /// into C / C++ array indexing.
  std::vector<MolObjIndexKind> entry_adjustment;

  std::vector<int> int_data;                 ///< Data for all entries, ordered by entries and
                                             ///<   by depth.  For entries A, B, and C with depth
                                             ///<   3: { A1, A2, A3, B1, B2, B3, C1, C2, C3 }.
  std::vector<double> real_data;             ///< Real-valued information for the property.  If a
                                             ///<   component of an entry is REAL, the
                                             ///<   corresponding value in int_data will refer to
                                             ///<   the index of data_str at which to find the
                                             ///<   information.
  std::vector<std::string> str_data;         ///< String data for all entries, covering the
                                             ///<   STRING enumerations of the MolObjPropField
                                             ///<   entry details.  If a component of an entry is a
                                             ///<   STRING, the corresponding value in int_data
                                             ///<   will refer to the index of data_str at which to
                                             ///<   find the information.
  std::vector<std::string> data_lines;       ///< Lines of data, stored one line per string, that
                                             ///<   complement the property

  /// \brief Extract the number of entries for the property.  Returns FALSE if there is no error
  ///        encountered, TRUE if there is a problem.  This will also mark the position on the
  ///        text line at which the entry count occurs.
  ///
  /// \param line_ptr   The line containing the property text
  /// \param start_pos  The starting position at which to read the value (default 6)
  /// \param length     The expected length of the entry count (default 3)
  bool readEntryCount(const char* line_ptr, int start_pos = 6, int length = 3);

  /// \brief Extract the index of the substrate atom or group for the property.  Returns FALSE if
  ///        there is no error encountered, TRUE if there is a problem.  This will also mark the
  ///        position on the text line at which the entry count occurs.
  ///
  /// \param line_ptr   The line containing the property text
  /// \param start_pos  The starting position at which to read the value (default 7)
  /// \param length     The expected length of the entry count (default 3)
  bool readSubstrateIndex(const char* line_ptr, int start_pos = 7, int length = 3);

  /// \brief Parse the list of entries for the property, filling out the arrays of integer, real,
  ///        and string data in the process.
  ///
  /// \param tf           The original text file (passed in for error-tracing purposes)
  /// \param line_number  Number of the line containing the property text
  /// \param start_pos    Starting position for data pertaining to the first entry
  /// \param limits       Bounds array for reading each entry.  Formatting present in entry_detail
  ///                     will impose some limits on what is read (white space can be implied by
  ///                     using a char4 to read between limits 5 characters apart, for example),
  ///                     but the critical feature is that the character at which to begin reading
  ///                     the nth data element is given in the nth position of limits.
  void parseEntries(const TextFile &tf, int line_number, int start_pos,
                    const std::vector<int> &limits);

  /// \brief Check the validity of an attribute request.  Getting information from properties can
  ///        be a tedious and error-prone process.  These guardrails will help developers find
  ///        information safely.
  ///
  /// \param entry_index      The index of the entry to access (receives a bounds check)
  /// \param attribute_index  The index of the attribute to retrieve (receives a bounds check)
  /// \param expectation      The expected nature of the attribute (checked for sanity)
  void checkAttributeValidity(int entry_index, int attribute_index,
                              MolObjPropField expectation) const;
};

} // namespace structure
} // namespace stormm

#endif
