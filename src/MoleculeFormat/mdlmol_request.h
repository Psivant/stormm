// -*-c++-*-
#ifndef STORMM_MDLMOL_REQUEST_H
#define STORMM_MDLMOL_REQUEST_H

#include <string>
#include <vector>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "Potential/energy_enumerators.h"
#include "molecule_format_enumerators.h"

namespace stormm {
namespace structure {

using energy::StateVariable;

/// \brief A request for one of a sanctioned list of information types to be included in a data
///        item of an SD file.  When printing an SD file, MdlMolDataItem objects will be created
///        based on these requests.
class MdlMolDataRequest {
public:

  /// \brief The constructor takes arguments corresponding to each member variable.  Various
  ///        overloads allow the object ot be constructed for a specific type of request, ignoring
  ///        member variables that are not relevant.
  /// \{
  MdlMolDataRequest(const std::string &title_in = std::string(""),
                    const std::string &label_in = std::string(""));

  MdlMolDataRequest(const std::string &title_in, StateVariable energy_component_in,
                    const std::string &label_in);

  MdlMolDataRequest(DataRequestKind kind_in, const std::string &title_in,
                    const std::string &message_in, const std::string &label_in);

  MdlMolDataRequest(const std::string &title_in, StateVariable valence_kind_in,
                    const std::vector<char4> &atom_types_in, const std::string &label_in);
  /// \}

  /// \brief The default copy and move constructors, as well as copy and move assignment operators,
  ///        are applicable to this object which has no const members or pointers to repair.
  /// \{
  MdlMolDataRequest(const MdlMolDataRequest &original) = default;
  MdlMolDataRequest(MdlMolDataRequest &&original) = default;
  MdlMolDataRequest& operator=(const MdlMolDataRequest &other) = default;
  MdlMolDataRequest& operator=(MdlMolDataRequest &&other) = default;
  /// \}

  /// \brief Get the kind of request.
  DataRequestKind getKind() const;

  /// \brief Get the title for the printed SD file's data item.
  const std::string& getTitle() const;

  /// \brief Get the type of energy requested, or produce a runtime error if this request is not
  ///        of the appropriate kind.
  StateVariable getEnergyComponent() const;

  /// \brief Get the atom mask string.  Raise a runtime error if this request is not of the
  ///        appropriate kind.
  const std::string& getAtomMask() const;

  /// \brief Get the valence parameter type.  Raise a runtime error if this request is not of the
  ///        appropriate kind.
  StateVariable getValenceParameter() const;

  /// \brief Get the custom string, again raising a runtime error if the request is not of the
  ///        appropriate kind.
  const std::string& getMessage() const;

  /// \brief Get the vector of atom types defining a valence parameter of interest, again raising
  ///        a runtime error if the request is not of the appropriate kind.
  const std::vector<char4>& getAtomTypes() const;

  /// \brief Get the system label to which the data item should be applied.
  const std::string& getSystemLabel() const;

  /// \brief Get the external registry number for the compound, as transcribed for this request
  const std::string& getExternalRegistryNumber() const;

  /// \brief Indicate that the MACCS-II field number should go in the resulting data item's header
  ///        line.
  bool placeMaccsFieldInHeader() const;

  /// \brief Indicate that the internal registry numbre (known only once the SD file is ready for
  ///        assembly) should be placed in the resulting data item's header.
  bool placeInternalRegistryInHeader() const;

  /// \brief Get the MACCS-II field number.
  int getMaccsFieldNumber() const;

  /// \brief Set the external registry number.
  ///
  /// \param regno_in  The external registry number for the compound described by this data item
  void setExternalRegistryNumber(const std::string &regno_in);

  /// \brief Set the MACCS-II field number.
  ///
  /// \param maccs_in  The MACCS-II database field number for the requested data item
  void setMaccsFieldNumber(int maccs_in);

  /// \brief Set the data item to make use of the SD file archive's internal registry number (the
  ///        number of the molecule or conformation within the file) in its header line.
  ///
  /// \param input  The directive can be set to "OFF" or "ON / ACTIVE"
  void setInternalRegistryUsage(const std::string &input);

private:
  DataRequestKind kind;           ///< Define the type of data request
  std::string title;              ///< Title for the data item, to be displayed as "> <title>" in
                                  ///<   the resulting SD file.
  StateVariable energy_component; ///< The type of energy quantity requested, i.e. PROPER_DIHEDRAL
  std::string atom_mask;          ///< Mask defining atoms for which valence interactions are to be
                                  ///<   listed in the SD file's data item, if the request is for
                                  ///<   atom influences
  StateVariable valence_kind;     ///< A type of valence parameter (including restraint terms)
                                  ///<   which, in conjunction with the appropriate number of named
                                  ///<   atom types, will cause the program to output all instances
                                  ///<   of the parameter in the model, listing the parameter's
                                  ///<   settings, atom indices to which it applies and the
                                  ///<   energies of each interaction.
  std::string message;            ///< A custom string provided by the user
  std::vector<char4> atom_types;  ///< Vector of atom types defining a valence
  std::string system_label;       ///< When applying the request to create data items in SD file
                                  ///<   output, this label will be checked against the molecular
                                  ///<   system's label to ensure that the information goes only
                                  ///<   into the intended places.

  // The following information is intended to cover the breadth of the SD file format's data item
  // header lines.
  bool use_maccs_ii_number;       ///< Flag to have the MACCS-II number used in the header line
  int maccs_ii_number;            ///< The field number of the property described by this data
                                  ///<   item in the MACCS-II database
  bool use_internal_registry;     ///< Flag to have the data item make use of the internal order
                                  ///<   of molecules present in the SD file archive in the
                                  ///<   header identification
  std::string external_regno;     ///< An external registry number for the compound to include
                                  ///<   in data items created in response to this request

  /// \brief Check that the data request is of the right kind to return a particular member
  ///        variable, or possibly to perform other functions.
  ///
  /// \param accepted_kind  The kind of request that this object must be making in order for it to
  ///                       function as the developer wants
  void checkKind(DataRequestKind accepted_kind) const;
};
  
} // namespace structure
} // namespace stormm

#endif
