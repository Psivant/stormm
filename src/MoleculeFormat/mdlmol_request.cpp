#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "mdlmol_request.h"

namespace stormm {
namespace structure {

using constants::CaseSensitivity;
using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
MdlMolDataRequest::MdlMolDataRequest(const std::string &title_in, const std::string &label_in) :
    kind{DataRequestKind::STRING}, title{title_in}, energy_component{StateVariable::BOND},
    atom_mask{std::string("")}, valence_kind{StateVariable::BOND}, message{std::string("")},
    atom_types{}, system_label{label_in}, use_maccs_ii_number{false}, maccs_ii_number{0},
    use_internal_registry{false}, external_regno{std::string("")}
{}

//-------------------------------------------------------------------------------------------------
MdlMolDataRequest::MdlMolDataRequest(const std::string &title_in,
                                     const StateVariable energy_component_in,
                                     const std::string &label_in) :
    MdlMolDataRequest(title_in, label_in)
{
  kind = DataRequestKind::STATE_VARIABLE;
  energy_component = energy_component_in;
}

//-------------------------------------------------------------------------------------------------
MdlMolDataRequest::MdlMolDataRequest(const DataRequestKind kind_in, const std::string &title_in,
                                     const std::string &message_in, const std::string &label_in) :
    MdlMolDataRequest(title_in, label_in)
{
  kind = kind_in;
  switch (kind) {
  case DataRequestKind::STATE_VARIABLE:
    rtErr("In order to construct a data request for a particular energy component, use "
          "MdlMolDataRequest(<title>, <energy component>, <label>).", "MdlMolDataRequest");
    break;
  case DataRequestKind::TOPOLOGY_PARAMETER:
    rtErr("In order to construct a data request for the actions of a topological energy "
          "parameter, use MdlMolDataRequest(<title>, <interaction type>, <list of atom types>, "
          "<label>).", "MdlMolDataRequest");
    break;
  case DataRequestKind::ATOM_INFLUENCES:
    atom_mask = message_in;
    break;
  case DataRequestKind::STRING:
    message = message_in;
    break;
  case DataRequestKind::ALL_KINDS:
    rtErr("The ALL_KINDS enumeration does not indicate an actual data request and no program "
          "should invoke it here.", "MdlMolDataRequest");
    break;
  }
}

//-------------------------------------------------------------------------------------------------
MdlMolDataRequest::MdlMolDataRequest(const std::string &title_in,
                                     const StateVariable valence_kind_in,
                                     const std::vector<char4> &atom_types_in,
                                     const std::string &label_in) :
    MdlMolDataRequest(title_in, label_in)
{
  kind = DataRequestKind::TOPOLOGY_PARAMETER;
  valence_kind = valence_kind_in;
  atom_types = atom_types_in;

  // Check the type of valence interaction and the corresponding number of atom types supplied.
  const int ntypes = atom_types.size();
  int nreq;
  switch (valence_kind) {
  case StateVariable::BOND:
  case StateVariable::UREY_BRADLEY:
    nreq = 2;
    break;
  case StateVariable::ANGLE:
    nreq = 3;
    break;
  case StateVariable::PROPER_DIHEDRAL:
  case StateVariable::IMPROPER_DIHEDRAL:
  case StateVariable::CHARMM_IMPROPER:
    nreq = 4;
    break;
  case StateVariable::CMAP:
    nreq = 5;
    break;
  case StateVariable::VDW:
  case StateVariable::VDW_ONE_FOUR:
  case StateVariable::ELECTROSTATIC:
  case StateVariable::ELEC_ONE_FOUR:
  case StateVariable::GENERALIZED_BORN:
  case StateVariable::RESTRAINT:
  case StateVariable::KINETIC:
  case StateVariable::PRESSURE:
  case StateVariable::VIRIAL_11:
  case StateVariable::VIRIAL_12:
  case StateVariable::VIRIAL_22:
  case StateVariable::VIRIAL_13:
  case StateVariable::VIRIAL_23:
  case StateVariable::VIRIAL_33:
  case StateVariable::VOLUME:
  case StateVariable::TEMPERATURE_ALL:
  case StateVariable::TEMPERATURE_PROTEIN:
  case StateVariable::TEMPERATURE_LIGAND:
  case StateVariable::TEMPERATURE_SOLVENT:
  case StateVariable::DU_DLAMBDA:
  case StateVariable::POTENTIAL_ENERGY:
  case StateVariable::TOTAL_ENERGY:
  case StateVariable::ALL_STATES:
    rtErr("The accepted topology parameter types for printing in data items of an SD file "
          "include BOND, ANGLE, PROPER_DIHEDRAL, IMPROPER_DIHEDRAL, UREY_BRADLEY, "
          "CHARMM_IMPROPER, and CMAP.", "MdlMolDataRequest");
    break;
  }
  if (ntypes != nreq) {
    rtErr("A " + getEnumerationName(valence_kind) + " requires specification of " +
          std::to_string(nreq) + " atom types, but " + std::to_string(ntypes) + " were provided.",
          "MdlMolDataRequest");
  }
}

//-------------------------------------------------------------------------------------------------
DataRequestKind MdlMolDataRequest::getKind() const {
  return kind;
}

//-------------------------------------------------------------------------------------------------
const std::string& MdlMolDataRequest::getTitle() const {
  return title;
}

//-------------------------------------------------------------------------------------------------
StateVariable MdlMolDataRequest::getEnergyComponent() const {
  checkKind(DataRequestKind::STATE_VARIABLE);
  return energy_component;
}

//-------------------------------------------------------------------------------------------------
const std::string& MdlMolDataRequest::getAtomMask() const {
  checkKind(DataRequestKind::ATOM_INFLUENCES);
  return atom_mask;
}

//-------------------------------------------------------------------------------------------------
StateVariable MdlMolDataRequest::getValenceParameter() const {
  checkKind(DataRequestKind::TOPOLOGY_PARAMETER);
  return valence_kind;
}

//-------------------------------------------------------------------------------------------------
const std::string& MdlMolDataRequest::getMessage() const {
  checkKind(DataRequestKind::STRING);
  return message;
}

//-------------------------------------------------------------------------------------------------
const std::vector<char4>& MdlMolDataRequest::getAtomTypes() const {
  return atom_types;
}

//-------------------------------------------------------------------------------------------------
const std::string& MdlMolDataRequest::getSystemLabel() const {
  return system_label;
}

//-------------------------------------------------------------------------------------------------
const std::string& MdlMolDataRequest::getExternalRegistryNumber() const {
  return external_regno;
}

//-------------------------------------------------------------------------------------------------
bool MdlMolDataRequest::placeMaccsFieldInHeader() const {
  return use_maccs_ii_number;
}

//-------------------------------------------------------------------------------------------------
bool MdlMolDataRequest::placeInternalRegistryInHeader() const {
  return use_internal_registry;
}

//-------------------------------------------------------------------------------------------------
int MdlMolDataRequest::getMaccsFieldNumber() const {
  return maccs_ii_number;
}

//-------------------------------------------------------------------------------------------------
void MdlMolDataRequest::setExternalRegistryNumber(const std::string &regno_in) {
  external_regno = regno_in;
}

//-------------------------------------------------------------------------------------------------
void MdlMolDataRequest::setMaccsFieldNumber(int maccs_in) {
  maccs_ii_number = maccs_in;
  use_maccs_ii_number = true;
}

//-------------------------------------------------------------------------------------------------
void MdlMolDataRequest::setInternalRegistryUsage(const std::string &input) {
  if (strcmpCased(input, "on", CaseSensitivity::NO) ||
      strcmpCased(input, "active", CaseSensitivity::NO)) {
    use_internal_registry = true;
  }
  else if (strcmpCased(input, "off", CaseSensitivity::NO)) {
    use_internal_registry = false;
  }
  else {
    rtErr("Invalid directive " + input + " to an SD file data item request, pertaining to "
          "internal registry number usage.  Use ON or OFF.", "MdlMolDataRequest",
          "setInternalRegistryUsage");
  }
}

//-------------------------------------------------------------------------------------------------
void MdlMolDataRequest::checkKind(const DataRequestKind accepted_kind) const {
  if (kind != accepted_kind) {
    rtErr("A data request of type " + getEnumerationName(kind) + " cannot function as a request "
          "for " + getEnumerationName(accepted_kind) + " data.", "MdlMolDataRequest",
          "checkKind");
  }
}

} // namespace structure
} // namespace stormm
