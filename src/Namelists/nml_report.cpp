#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
#  include <windows.h>
#  include <Lmcons.h>
#else
#  include <unistd.h>
#endif
#include "copyright.h"
#include "Constants/behavior.h"
#include "MoleculeFormat/molecule_format_enumerators.h"
#include "Parsing/parse.h"
#include "Reporting/summary_file.h"
#include "namelist_enumerators.h"
#include "nml_report.h"

namespace stormm {
namespace namelist {

using constants::CaseSensitivity;
using energy::translateEnergySample;
using parse::minimalRealFormat;
using parse::strcmpCased;
using parse::stringToChar4;
using review::default_output_file_width;
using review::translateOutputScope;
using structure::DataRequestKind;
  
//-------------------------------------------------------------------------------------------------
ReportControls::ReportControls(const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    policy{policy_in}, report_layout{OutputSyntax::STANDALONE},
    report_scope{OutputScope::AVERAGES}, state_sampling{EnergySample::TIME_SERIES},
    username{std::string("")},
    report_variable{std::string(default_report_variable_name)},
    start_date{}, print_walltime_data{true},
    report_file_width{default_output_file_width},
    common_path_limit{default_common_path_limit},
    common_path_threshold{default_common_path_threshold},
    energy_decimal_places{default_energy_decimal_places},
    outlier_sigma_factor{default_energy_outlier_sigmas},
    outlier_count{default_outlier_limit},
    reported_quantities{}, sdf_addons{},
    nml_transcript{"report"}
{
  // Detect the time of day
  gettimeofday(&start_date, nullptr);

  // Get the username
#if defined(_WIN32) || defined(_WIN64) || defined(__CYGWIN__)
  std::string buffer(UNLEN+1, ' ');
  DWORD buffer_len = UNLEN+1;
  GetUserName(buffer.data(), &buffer_len);
#else
  std::string buffer(64, ' ');
  getlogin_r(buffer.data(), buffer.size());
  buffer.resize(strlen(buffer.data()));
#endif
  username = buffer;
}

//-------------------------------------------------------------------------------------------------
ReportControls::ReportControls(const TextFile &tf, int *start_line, bool *found_nml,
                               const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    ReportControls(policy_in, wrap)
{
  NamelistEmulator t_nml = reportInput(tf, start_line, found_nml, policy, wrap);
  nml_transcript = t_nml;
  if (t_nml.getKeywordStatus("syntax") != InputStatus::MISSING) {
    setOutputSyntax(t_nml.getStringValue("syntax"));
  }
  if (t_nml.getKeywordStatus("scope") != InputStatus::MISSING) {
    setOutputScope(t_nml.getStringValue("scope"));
  }
  if (t_nml.getKeywordStatus("nrgsample") != InputStatus::MISSING) {
    setStateSampling(t_nml.getStringValue("nrgsample"));
  }
  if (t_nml.getKeywordStatus("username") != InputStatus::MISSING) {
    setUsername(t_nml.getStringValue("username"));
  }
  if (t_nml.getKeywordStatus("varname") != InputStatus::MISSING) {
    setReportVariable(t_nml.getStringValue("varname"));
  }
  if (t_nml.getKeywordStatus("timings") != InputStatus::MISSING) {
    setWallTimeData(t_nml.getStringValue("timings"));
  }
  report_file_width = t_nml.getIntValue("report_width");
  if (t_nml.getKeywordStatus("energy") != InputStatus::MISSING) {
    const int ndetail = t_nml.getKeywordEntries("energy");
    for (int i = 0; i < ndetail; i++) {
      setReportedQuantities(t_nml.getStringValue("energy", i));
    }
  }
  if (t_nml.getKeywordStatus("state") != InputStatus::MISSING) {
    const int ndetail = t_nml.getKeywordEntries("state");
    for (int i = 0; i < ndetail; i++) {
      setReportedQuantities(translateStateQuantity(t_nml.getStringValue("state", i)));
    }
  }
  if (t_nml.getKeywordStatus("sdf_item") != InputStatus::MISSING) {
    const int ndata = t_nml.getKeywordEntries("sdf_item");
    for (int i = 0; i < ndata; i++) {
      const std::vector<MdlMolDataRequest> new_items = translateSdfKeywordInput(t_nml, i);
      for (size_t j = 0; j < new_items.size(); j++) {
        addDataItem(new_items[j]);
      }
    }
  }
  common_path_limit = t_nml.getIntValue("common_path_limit");
  common_path_threshold = t_nml.getIntValue("common_path_threshold");
  energy_decimal_places = t_nml.getIntValue("e_precision");
  outlier_sigma_factor = t_nml.getRealValue("outlier_sigmas");
  outlier_count = t_nml.getIntValue("outlier_count");

  // Validate inputs found thus far
  validateEnergyDecimalPlaces();
  validateOutlierMetrics();
}

//-------------------------------------------------------------------------------------------------
OutputSyntax ReportControls::getOutputSyntax() const {
  return report_layout;
}

//-------------------------------------------------------------------------------------------------
OutputScope ReportControls::getOutputScope() const {
  return report_scope;
}

//-------------------------------------------------------------------------------------------------
EnergySample ReportControls::getEnergySamplingMethod() const {
  return state_sampling;
}

//-------------------------------------------------------------------------------------------------
const std::string& ReportControls::getReportVariable() const {
  return report_variable;
}

//-------------------------------------------------------------------------------------------------
const std::string& ReportControls::getUsername() const {
  return username;
}

//-------------------------------------------------------------------------------------------------
timeval ReportControls::getStartDate() const {
  return start_date;
}

//-------------------------------------------------------------------------------------------------
const NamelistEmulator& ReportControls::getTranscript() const {
  return nml_transcript;
}

//-------------------------------------------------------------------------------------------------
bool ReportControls::printWallTimeData() const {
  return print_walltime_data;
}

//-------------------------------------------------------------------------------------------------
int ReportControls::getReportFileWidth() const {
  return report_file_width;
}

//-------------------------------------------------------------------------------------------------
int ReportControls::getReportedQuantityCount() const {
  return reported_quantities.size();
}
  
//-------------------------------------------------------------------------------------------------
const std::vector<StateVariable>& ReportControls::getReportedQuantities() const {
  return reported_quantities;
}

//-------------------------------------------------------------------------------------------------
int ReportControls::getSDFileDataRequestCount() const {
  return sdf_addons.size();
}

//-------------------------------------------------------------------------------------------------
const std::vector<MdlMolDataRequest>& ReportControls::getSDFileDataRequests() const {
  return sdf_addons;
}

//-------------------------------------------------------------------------------------------------
MdlMolDataRequest ReportControls::getSDFileDataRequest(const int index) const {
  if (index < 0 || index > static_cast<int>(sdf_addons.size())) {
    rtErr("A &report namelist with " + std::to_string(sdf_addons.size()) + " requests cannot "
          "produce request index " + std::to_string(index) + ".", "ReportControls",
          "getSDFileDataRequest");
  }
  return sdf_addons[index];
}

//-------------------------------------------------------------------------------------------------
int ReportControls::getCommonPathLimit() const {
  return common_path_limit;
}

//-------------------------------------------------------------------------------------------------
int ReportControls::getCommonPathThreshold() const {
  return common_path_threshold;
}

//-------------------------------------------------------------------------------------------------
int ReportControls::getEnergyDecimalPlaces() const {
  return energy_decimal_places;
}

//-------------------------------------------------------------------------------------------------
double ReportControls::getOutlierSigmaFactor() const {
  return outlier_sigma_factor;
}

//-------------------------------------------------------------------------------------------------
int ReportControls::getOutlierCount() const {
  return outlier_count;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setOutputSyntax(const OutputSyntax report_layout_in) {
  report_layout = report_layout_in;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setOutputSyntax(const std::string &report_layout_in) {
  if (strcmpCased(report_layout_in, "matplotlib", CaseSensitivity::NO) ||
      strcmpCased(report_layout_in, "mat_plot_lib", CaseSensitivity::NO)) {
    report_layout = OutputSyntax::MATPLOTLIB;
  }
  else if (strcmpCased(report_layout_in, "matlab", CaseSensitivity::NO) ||
           strcmpCased(report_layout_in, "octave", CaseSensitivity::NO) ||
           strcmpCased(report_layout_in, "matrix_pkg", CaseSensitivity::NO)) {
    report_layout = OutputSyntax::MATRIX_PKG;
  }
  else if (strcmpCased(report_layout_in, "standalone", CaseSensitivity::NO) ||
           strcmpCased(report_layout_in, "stormm", CaseSensitivity::NO)) {
    report_layout = OutputSyntax::STANDALONE;
  }
  else {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(report_layout_in + " was not recognized as a valid output syntax.", "ReportControls",
            "setOutputSyntax");
    case ExceptionResponse::WARN:
      rtWarn(report_layout_in + " was not recognized as a valid output syntax.  The " +
             getEnumerationName(OutputSyntax::STANDALONE) + " format will be taken instead.",
             "ReportControls", "setOutputSyntax");
      report_layout = OutputSyntax::STANDALONE;
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setOutputScope(const OutputScope report_scope_in) {
  report_scope = report_scope_in;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setOutputScope(const std::string &report_scope_in) {
  try {
    report_scope = translateOutputScope(report_scope_in);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(report_scope_in + " was not recognized as a valid scope for reporting energies and "
            "other system diagnostics.", "ReportControls", "setOutputSyntax");
    case ExceptionResponse::WARN:
      rtWarn(report_scope_in + " was not recognized as a valid output syntax.  The " +
             getEnumerationName(OutputScope::AVERAGES) + " format will be taken instead.",
             "ReportControls", "setOutputScope");
      report_scope = OutputScope::AVERAGES;
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setStateSampling(EnergySample state_sampling_in) {
  state_sampling = state_sampling_in;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setStateSampling(const std::string &state_sampling_in) {
  try {
    state_sampling = translateEnergySample(state_sampling_in);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(state_sampling_in + " was not recognized as a valid scope for sampling energies and "
            "other system diagnostics.", "ReportControls", "setStateSampling");
    case ExceptionResponse::WARN:
      rtWarn(state_sampling_in + " was not recognized as a valid output syntax.  The " +
             getEnumerationName(EnergySample::TIME_SERIES) + " format will be taken instead.",
             "ReportControls", "setStateSampling");
      state_sampling = EnergySample::TIME_SERIES;
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setUsername(const std::string &username_in) {
  username = username_in;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setReportVariable(const std::string &report_variable_in) {
  report_variable = report_variable_in;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setWallTimeData(const bool preference) {
  print_walltime_data = preference;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setWallTimeData(const std::string &preference) {
  if (strcmpCased(preference, "on") || strcmpCased(preference, "active")) {
    print_walltime_data = true;
  }
  else if (strcmpCased(preference, "off")) {
    print_walltime_data = false;
  }
  else {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(preference + " is not a recognized setting for the \"timings\" keyword.  Use ON / "
            "ACTIVE or OFF.", "ReportControls", "setWallTimeData");
    case ExceptionResponse::WARN:
      rtWarn(preference + " is not a recognized setting for the \"timings\" keyword.  Timings "
             "display will remain ON / ACTIVE unless set to OFF.", "ReportControls",
             "setWallTimeData");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setReportFileWidth(const int report_file_width_in) {
  report_file_width = report_file_width_in;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setReportedQuantities(const std::string &quantity_in) {
  setReportedQuantities(translateEnergyComponent(quantity_in));
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setReportedQuantities(const StateVariable quantities_in) {
  setReportedQuantities(std::vector<StateVariable>(1, quantities_in));
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setReportedQuantities(const std::vector<StateVariable> &quantities_in) {
  const int ns = static_cast<int>(StateVariable::ALL_STATES);
  std::vector<bool> activated(ns, false);
  const int nq = quantities_in.size();
  const int ncurr = reported_quantities.size();
  int n_monitored = 0;
  for (int i = 0; i < ncurr; i++) {
    const int qno = static_cast<int>(reported_quantities[i]);
    if (reported_quantities[i] != StateVariable::ALL_STATES && activated[qno] == false) {
      activated[qno] = true;
      n_monitored++;
    }    
  }
  for (int i = 0; i < nq; i++) {
    const int qno = static_cast<int>(quantities_in[i]);
    if (quantities_in[i] != StateVariable::ALL_STATES && activated[qno] == false) {
      activated[qno] = true;
      n_monitored++;
    }
  }
  reported_quantities.resize(n_monitored);
  int j = 0;
  for (int i = 0LLU; i < ns; i++) {
    if (activated[i]) {
      reported_quantities[j] = static_cast<StateVariable>(i);
      j++;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ReportControls::addDataItem(const MdlMolDataRequest &ask) {
  
  // Check that there is no other data item with a conflicting title.
  bool problem = false;
  const int n_items = sdf_addons.size();
  for (int i = 0; i < n_items; i++) {
    problem = (problem || (sdf_addons[i].getTitle() == ask.getTitle()));
  }
  if (problem == false) {
    sdf_addons.push_back(ask);
  }
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setCommonPathLimit(const int common_path_limit_in) {
  common_path_limit = common_path_limit_in;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setCommonPathThreshold(const int common_path_threshold_in) {
  common_path_threshold = common_path_threshold_in;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::setEnergyDecimalPlaces(const int energy_decimal_places_in) {
  energy_decimal_places = energy_decimal_places_in;
}

//-------------------------------------------------------------------------------------------------
std::vector<StateVariable> ReportControls::translateEnergyComponent(const std::string &inpstr) {
  std::vector<StateVariable> result;
  const CaseSensitivity noc = CaseSensitivity::NO;
  if (strcmpCased(inpstr, "bond", noc) || strcmpCased(inpstr, "bonds", noc)) {
    result.push_back(StateVariable::BOND);
  }
  else if (strcmpCased(inpstr, "angle", noc) || strcmpCased(inpstr, "angles", noc)) {
    result.push_back(StateVariable::ANGLE);
    result.push_back(StateVariable::UREY_BRADLEY);
  }
  else if (strcmpCased(inpstr, "harmonic_angle", noc) ||
           strcmpCased(inpstr, "harmonicangle", noc)) {
    result.push_back(StateVariable::ANGLE);
  }
  else if (strcmpCased(inpstr, "urey_bradley", noc) || strcmpCased(inpstr, "ureybradley", noc)) {
    result.push_back(StateVariable::UREY_BRADLEY);
  }
  else if (strcmpCased(inpstr, "dihedral", noc) || strcmpCased(inpstr, "torsion", noc) ||
           strcmpCased(inpstr, "dihedrals", noc) || strcmpCased(inpstr, "torsions", noc)) {
    result.push_back(StateVariable::PROPER_DIHEDRAL);
    result.push_back(StateVariable::IMPROPER_DIHEDRAL);
    result.push_back(StateVariable::CHARMM_IMPROPER);
  }
  else if (strcmpCased(inpstr, "proper_dihedral", noc) ||
           strcmpCased(inpstr, "proper_dihedrals", noc) ||
           strcmpCased(inpstr, "properdihedral", noc) ||
           strcmpCased(inpstr, "properdihedrals", noc)) {
    result.push_back(StateVariable::PROPER_DIHEDRAL);    
  }
  else if (strcmpCased(inpstr, "improper_dihedral", noc) ||
           strcmpCased(inpstr, "improper_dihedrals", noc) ||
           strcmpCased(inpstr, "improperdihedral", noc) ||
           strcmpCased(inpstr, "improperdihedrals", noc)) {
    result.push_back(StateVariable::IMPROPER_DIHEDRAL);
  }
  else if (strcmpCased(inpstr, "charmm_improper", noc) ||
           strcmpCased(inpstr, "charmmimproper", noc)) {
    result.push_back(StateVariable::CHARMM_IMPROPER);
  }
  else if (strcmpCased(inpstr, "cmap", noc) || strcmpCased(inpstr, "cmaps", noc)) {
    result.push_back(StateVariable::CMAP);
  }
  else if (strcmpCased(inpstr, "vdw", noc) || strcmpCased(inpstr, "van_der_waals", noc) ||
           strcmpCased(inpstr, "vanderwaals", noc) || strcmpCased(inpstr, "lennardjones", noc) ||
           strcmpCased(inpstr, "lj", noc)) {
    result.push_back(StateVariable::VDW);
    result.push_back(StateVariable::VDW_ONE_FOUR);
  }
  else if (strcmpCased(inpstr, "vdw_nonbonded", noc) || strcmpCased(inpstr, "vdwnonbonded", noc) ||
           strcmpCased(inpstr, "nonbonded_lj", noc) || strcmpCased(inpstr, "nonbondedlj", noc) ||
           strcmpCased(inpstr, "lj_nonbonded", noc) || strcmpCased(inpstr, "ljnonbonded", noc)) {
    result.push_back(StateVariable::VDW);
  }
  else if (strcmpCased(inpstr, "vdw_14", noc) || strcmpCased(inpstr, "vdw14", noc) ||
           strcmpCased(inpstr, "vdw_near", noc) || strcmpCased(inpstr, "vdwnear", noc) ||
           strcmpCased(inpstr, "lj_14", noc) || strcmpCased(inpstr, "lj14", noc) ||
           strcmpCased(inpstr, "lj_near", noc) || strcmpCased(inpstr, "ljnear", noc) ||
           strcmpCased(inpstr, "near_lj", noc) || strcmpCased(inpstr, "nearlj", noc)) {
    result.push_back(StateVariable::VDW_ONE_FOUR);
  }
  else if (strcmpCased(inpstr, "elec", noc) || strcmpCased(inpstr, "electrostatic", noc)) {
    result.push_back(StateVariable::ELECTROSTATIC);
    result.push_back(StateVariable::ELEC_ONE_FOUR);
  }
  else if (strcmpCased(inpstr, "elec_nonbonded", noc) ||
           strcmpCased(inpstr, "elecnonbonded", noc) ||
           strcmpCased(inpstr, "electrostatic_nonbonded", noc) ||
           strcmpCased(inpstr, "electrostaticnonbonded", noc) ||
           strcmpCased(inpstr, "nonbonded_electrostatic", noc) ||
           strcmpCased(inpstr, "nonbondedelectrostatic", noc)) {
    result.push_back(StateVariable::ELECTROSTATIC);
  }
  else if (strcmpCased(inpstr, "elec_14", noc) || strcmpCased(inpstr, "elec14", noc) ||
           strcmpCased(inpstr, "elec_near", noc) || strcmpCased(inpstr, "elecnear", noc) ||
           strcmpCased(inpstr, "electrostatic_near", noc) ||
           strcmpCased(inpstr, "near_electrostatic", noc) ||
           strcmpCased(inpstr, "electrostaticnear", noc) ||
           strcmpCased(inpstr, "nearelectrostatic", noc)) {
    result.push_back(StateVariable::ELEC_ONE_FOUR);
  }
  else if (strcmpCased(inpstr, "gb", noc) || strcmpCased(inpstr, "generalized_born", noc) ||
           strcmpCased(inpstr, "solvent", noc) || strcmpCased(inpstr, "generalizedborn", noc)) {
    result.push_back(StateVariable::GENERALIZED_BORN);
  }
  else if (strcmpCased(inpstr, "restraint", noc) || strcmpCased(inpstr, "nmr", noc)) {
    result.push_back(StateVariable::RESTRAINT);
  }
  else if (strcmpCased(inpstr, "total_energy", noc) || strcmpCased(inpstr, "totalenergy", noc) ||
           strcmpCased(inpstr, "etot", noc) || strcmpCased(inpstr, "total", noc)) {
    result.push_back(StateVariable::TOTAL_ENERGY);
  }
  else if (strcmpCased(inpstr, "kinetic", noc)) {
    result.push_back(StateVariable::KINETIC);
  }
  else if (strcmpCased(inpstr, "potential_energy", noc) || strcmpCased(inpstr, "potential", noc) ||
           strcmpCased(inpstr, "potentialenergy", noc) || strcmpCased(inpstr, "pe", noc)) {
    result.push_back(StateVariable::POTENTIAL_ENERGY);
  }
  else {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(inpstr + " was not recognized as a molecular mechanics energy component valid "
            "for focused reporting.  Providing no input in this section will have all molecular "
            "mechanics terms be reported.", "ReportControls", "setReportedQuantities");
    case ExceptionResponse::WARN:
      rtWarn(inpstr + " was not recognized as a molecular mechanics energy component valid "
             "for focused reporting.  Providing no input in this section will have all molecular "
             "mechanics terms be reported.  The present input will be ignored.", "ReportControls",
             "setReportedQuantities");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<StateVariable> ReportControls::translateStateQuantity(const std::string &inpstr) {
  std::vector<StateVariable> result;
  const CaseSensitivity noc = CaseSensitivity::NO;
  if (strcmpCased(inpstr, "temperature", noc) || strcmpCased(inpstr, "temp", noc)) {
    result.push_back(StateVariable::TEMPERATURE_ALL);
  }
  else if (strcmpCased(inpstr, "local_temperature", noc)) {
    result.push_back(StateVariable::TEMPERATURE_PROTEIN);
    result.push_back(StateVariable::TEMPERATURE_LIGAND);
    result.push_back(StateVariable::TEMPERATURE_SOLVENT);
  }
  else if (strcmpCased(inpstr, "pressure", noc)) {
    result.push_back(StateVariable::PRESSURE);
  }
  else if (strcmpCased(inpstr, "virial_components", noc)) {
    result.push_back(StateVariable::VIRIAL_11);
    result.push_back(StateVariable::VIRIAL_12);
    result.push_back(StateVariable::VIRIAL_22);
    result.push_back(StateVariable::VIRIAL_13);
    result.push_back(StateVariable::VIRIAL_23);
    result.push_back(StateVariable::VIRIAL_33);
  }
  else if (strcmpCased(inpstr, "volume", noc) || strcmpCased(inpstr, "vol", noc)) {
    result.push_back(StateVariable::VOLUME);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<MdlMolDataRequest>
ReportControls::translateSdfKeywordInput(const NamelistEmulator &t_nml, const int index) {

  // Get the title and system label.
  const std::string ttl = t_nml.getStringValue("sdf_item", "-title", index);
  const std::string lbl = t_nml.getStringValue("sdf_item", "-label", index);

  // Search for clues as to what the data item could be about.
  const InputStatus missing = InputStatus::MISSING;
  std::vector<bool> rqt(static_cast<size_t>(DataRequestKind::ALL_KINDS), false);
  rqt[static_cast<size_t>(DataRequestKind::STATE_VARIABLE)] =
    (t_nml.getKeywordStatus("sdf_item", "-energy", index) != missing);
  rqt[static_cast<size_t>(DataRequestKind::ATOM_INFLUENCES)] =
    (t_nml.getKeywordStatus("sdf_item", "-mask", index) != missing);
  rqt[static_cast<size_t>(DataRequestKind::TOPOLOGY_PARAMETER)] =
    (t_nml.getKeywordStatus("sdf_item", "-parameter", index) != missing);
  rqt[static_cast<size_t>(DataRequestKind::STRING)] =
    (t_nml.getKeywordStatus("sdf_item", "-message", index) != missing);
  
  // The data item can only be one thing.  Request types will be prioritized in the order in which
  // they appear in the dedicated enum class.
  DataRequestKind result_kind = DataRequestKind::ALL_KINDS;
  for (int i = 0; i < static_cast<int>(DataRequestKind::ALL_KINDS); i++) {
    if (rqt[i] == false) {
      continue;
    }
    result_kind = static_cast<DataRequestKind>(i);
    for (int j = i + 1; j < static_cast<int>(DataRequestKind::ALL_KINDS); j++) {
      if (rqt[j]) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("An SD file data item can encode only one type of information.  Both " +
                getEnumerationName(static_cast<DataRequestKind>(i)) + " and " +
                getEnumerationName(static_cast<DataRequestKind>(j)) + " were indicated.",
                "ReportControls", "translateSdfKeywordInput");
        case ExceptionResponse::WARN:
          rtWarn("An SD file data item can encode only one type of information.  Both " +
                 getEnumerationName(static_cast<DataRequestKind>(i)) + " and " +
                 getEnumerationName(static_cast<DataRequestKind>(j)) + " were indicated.  " +
                 getEnumerationName(static_cast<DataRequestKind>(i)) + " will take precedence.",
                 "ReportControls", "translateSdfKeywordInput");
          rqt[j] = false;
          break;
        case ExceptionResponse::SILENT:
          rqt[j] = false;
          break;
        }
      }
    }
  }

  // Create the appropriate result.  Track the additions so that further edits may take place
  // later in this function.
  std::vector<MdlMolDataRequest> result;
  switch (result_kind) {
  case DataRequestKind::STATE_VARIABLE:
    {
      // Do not count composite terms in the SD file energy recordings.
      const std::string& detail_name = t_nml.getStringValue("sdf_item", "-energy", index);
      const std::vector<StateVariable> tmp_sv = translateEnergyComponent(detail_name);
      if (tmp_sv.size() > 1LLU) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr(detail_name + " implies a composite of " + std::to_string(tmp_sv.size()) +
                " molecular mechanics energy components.  Use a term corresponding to only one "
                "energy component.", "ReportControls", "translateSdfKeywordInput");
        case ExceptionResponse::WARN:
          rtWarn(detail_name + " implies a composite of " + std::to_string(tmp_sv.size()) +
                 " molecular mechanics energy components.  Use a term corresponding to only one "
                 "energy component.  This item will be omitted from the SD file results.",
                 "ReportControls", "translateSdfKeywordInput");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }
      else {
        for (size_t i = 0; i < tmp_sv.size(); i++) {
          result.emplace_back(ttl, tmp_sv[i], lbl);
        }
      }
    }
    break;
  case DataRequestKind::ATOM_INFLUENCES:
    result.emplace_back(DataRequestKind::ATOM_INFLUENCES, ttl,
                        t_nml.getStringValue("sdf_item", "-mask", index), lbl);
    break;
  case DataRequestKind::TOPOLOGY_PARAMETER:
    {
      std::vector<char4> atom_types;
      if (t_nml.getKeywordStatus("sdf_item", "-typeI", index) != missing) {
        atom_types.push_back(stringToChar4(t_nml.getStringValue("sdf_item", "-typeI", index)));
      }
      if (t_nml.getKeywordStatus("sdf_item", "-typeJ", index) != missing) {
        atom_types.push_back(stringToChar4(t_nml.getStringValue("sdf_item", "-typeJ", index)));
      }
      if (t_nml.getKeywordStatus("sdf_item", "-typeK", index) != missing) {
        atom_types.push_back(stringToChar4(t_nml.getStringValue("sdf_item", "-typeK", index)));
      }
      if (t_nml.getKeywordStatus("sdf_item", "-typeL", index) != missing) {
        atom_types.push_back(stringToChar4(t_nml.getStringValue("sdf_item", "-typeL", index)));
      }
      if (t_nml.getKeywordStatus("sdf_item", "-typeM", index) != missing) {
        atom_types.push_back(stringToChar4(t_nml.getStringValue("sdf_item", "-typeM", index)));
      }
      const std::string& sval = t_nml.getStringValue("sdf_item", "-parameter", index);
      StateVariable term_type = StateVariable::ALL_STATES;
      int ntypes_req;
      const CaseSensitivity noc = CaseSensitivity::NO;
      if (strcmpCased(sval, "bond", noc)) {
        ntypes_req = 2;
        term_type = StateVariable::BOND;
      }
      else if (strcmpCased(sval, "angle", noc) || strcmpCased(sval, "harmonic_angle", noc) ||
               strcmpCased(sval, "harmonicangle", noc)) {
        ntypes_req = 3;
        term_type = StateVariable::ANGLE;
      }
      else if (strcmpCased(sval, "proper_dihedral", noc) ||
               strcmpCased(sval, "properdihedral", noc) ||
               strcmpCased(sval, "proper_torsion", noc) ||
               strcmpCased(sval, "propertorsion", noc)) {
        ntypes_req = 4;
        term_type = StateVariable::PROPER_DIHEDRAL;
      }
      else if (strcmpCased(sval, "improper_dihedral", noc) ||
               strcmpCased(sval, "improperdihedral", noc) ||
               strcmpCased(sval, "improper_torsion", noc) ||
               strcmpCased(sval, "impropertorsion", noc)) {
        ntypes_req = 4;
        term_type = StateVariable::IMPROPER_DIHEDRAL;
      }
      else if (strcmpCased(sval, "urey_bradley", noc) || strcmpCased(sval, "ureybradley", noc)) {
        ntypes_req = 2;
        term_type = StateVariable::UREY_BRADLEY;
      }
      else if (strcmpCased(sval, "charmm_improper", noc) ||
               strcmpCased(sval, "charmmimproper", noc)) {
        ntypes_req = 4;
        term_type = StateVariable::CHARMM_IMPROPER;
      }
      else if (strcmpCased(sval, "cmap", noc) || strcmpCased(sval, "correction_map", noc) ||
               strcmpCased(sval, "correctionmap", noc)) {
        ntypes_req = 5;
        term_type = StateVariable::CMAP;
      }

      // Report or die if any errors were encountered.
      if (term_type == StateVariable::ALL_STATES) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("No valid energy term could be translated from \"" + sval + "\" in sdf_item " +
                std::to_string(index + 1) + ".  Name a valid valence term in this field.",
                "ReportControls", "translateSdfKeywordInput");
        case ExceptionResponse::WARN:
          rtWarn("No valid energy term could be translated from \"" + sval + "\" in sdf_item " +
                 std::to_string(index + 1) + ".  No such output will be produced.",
                 "ReportControls", "translateSdfKeywordInput");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }
      else if (static_cast<int>(atom_types.size()) != ntypes_req) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("The energy term \"" + sval + "\" is defined by " + std::to_string(ntypes_req) +
                " atom types, but " + std::to_string(atom_types.size()) + " types were provided.",
                "ReportControls", "translateSdfKeywordInput");
        case ExceptionResponse::WARN:
          rtWarn("The energy term \"" + sval + "\" is defined by " + std::to_string(ntypes_req) +
                 " atom types, but " + std::to_string(atom_types.size()) + " types were "
                 "provided.  This item will be skipped.", "ReportControls",
                 "translateSdfKeywordInput");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }
      else {
        result.emplace_back(ttl, term_type, atom_types, lbl);
      }
    }
    break;
  case DataRequestKind::STRING:
    result.emplace_back(DataRequestKind::STRING, ttl,
                        t_nml.getStringValue("sdf_item", "-message", index), lbl);
    break;
  case DataRequestKind::ALL_KINDS:
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("No kind was discernible for sdf_item " + std::to_string(index) + " in the &report "
            "namelist.", "ReportControls", "translateSdfKeywordInput");
    case ExceptionResponse::WARN:
      rtWarn("No kind was discernible for sdf_item " + std::to_string(index) + " in the &report "
             "namelist.", "ReportControls", "translateSdfKeywordInput");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }

  // Further adjustments based on additional user input
  for (size_t i = 0LLU; i < result.size(); i++) {
    if (t_nml.getKeywordStatus("sdf_item", "-internal_regno", index) != missing) {
      result[i].setInternalRegistryUsage(t_nml.getStringValue("sdf_item", "-internal_regno",
                                                              index));
    }
    if (t_nml.getKeywordStatus("sdf_item", "-maccsid", index) != missing) {
      result[i].setMaccsFieldNumber(t_nml.getIntValue("sdf_item", "-maccsid", index));    
    }
  }  
  return result;
}

//-------------------------------------------------------------------------------------------------
void ReportControls::validateEnergyDecimalPlaces() {
  if (energy_decimal_places < 1 || energy_decimal_places > 12) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Invalid number of digits after the decimal for reporting energies: " +
            std::to_string(energy_decimal_places) + ".", "ReportControls",
            "validateEnergyDecimalPlaces");
    case ExceptionResponse::WARN:
      rtWarn("An invalid number of digits after the decimal was specified for reporting "
             "energies (" + std::to_string(energy_decimal_places) + ").  The default value of " +
             std::to_string(default_energy_decimal_places) + " will be restored.",
             "ReportControls", "validateEnergyDecimalPlaces");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    energy_decimal_places = default_energy_decimal_places;
  }
}

//-------------------------------------------------------------------------------------------------
void ReportControls::validateOutlierMetrics() {
  if (outlier_sigma_factor <= 0.0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Invalid multiplier for detecting outliers: " +
            minimalRealFormat(outlier_sigma_factor, 1.0e-4) + " standard deviations from the "
            "mean.", "ReportControls", "validateOutlierMetrics");
    case ExceptionResponse::WARN:
      rtWarn("An invalid multiplier (" + minimalRealFormat(outlier_sigma_factor, 1.0e-4) +
             " standard deviations from the mean, was provided.  The default value of " +
             minimalRealFormat(default_energy_outlier_sigmas, 1.0e-4) + " will be restored.",
             "ReportControls", "validateOutlierMetrics");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    outlier_sigma_factor = default_energy_outlier_sigmas;
  }
  if (outlier_count < 0) {
    outlier_count = 0;
  }
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator reportInput(const TextFile &tf, int *start_line, bool *found,
                             const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("report", CaseSensitivity::AUTOMATIC, policy, "Collects directives "
                         "pertaining to the content and layout of the diagnostics output file.  "
                         "This degree of control, while optional, can be important for managing "
                         "the voluminous output that can come from a program which runs multiple "
                         "simulations in a single runtime process.");
  const std::string sdf_help("Transfer one of a recognized assortment of energy or other "
                             "structural quantities computed by STORMM to a data item in an SD "
                             "file (.sdf format).  The results will appear after the MDL MOL "
                             "section detailing atoms, properties, and bonding patterns.");
  const std::vector<std::string> sdf_keys_help = {
    "Custom title of the data item, which will appear as \"> <TITLE>\" in the resulting SD file",
    "System label linking the data item to particular systems in the calculation.  Only those "
    "systems matching this label will have the data item in their output.",
    "One of a number of recognized energy terms for singular molecular mechanics energy "
    "components.  Accepted arguments include BOND, HARMONIC_ANGLE, PROPER_DIHEDRAL, "
    "IMPROPER_DIHEDRAL, UREY_BRADLEY, CHARMM_IMPROPER, CMAP, ELECTROSTATIC (NEAR or NONBONDED), "
    "VDW (NEAR OR NONBONDED), GENERALIZED_BORN, RESTRAINT, TOTAL_ENERGY, and POTENTIAL_ENERGY.",
    "A customizable message from the user conveyed into the SD file output",
    "An atom mask indicating a subset of the atoms for which to enumerate relevant valence "
    "interactions",
    "A type of parameter for which to search the relevant topology, in conjunction with the "
    "following atom types.  If found, the settings of the parameter will be printed to the SD "
    "file output.",
    "An atom type that helps define the topology parameter being sought",
    "An atom type that helps define the topology parameter being sought",
    "An atom type that helps define the topology parameter being sought",
    "An atom type that helps define the topology parameter being sought",
    "An atom type that helps define the topology parameter being sought",
    "The external registry number for the compound (will be displayed on the headline in "
    "parentheses)",
    "Set to ON to have this data item state the internal registry number within the SDF on its "
    "headline",
    "The field number of this type of data in a MACCS-II database"
  };
  t_nml.addKeyword(NamelistElement("syntax", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("scope", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("nrgsample", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("username", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("varname", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("timings", NamelistType::STRING, "MISSING"));
  t_nml.addKeyword(NamelistElement("report_width", NamelistType::INTEGER,
                                   std::to_string(default_output_file_width)));
  t_nml.addKeyword(NamelistElement("energy", NamelistType::STRING, "MISSING",
                                   DefaultIsObligatory::NO, InputRepeats::YES));
  t_nml.addKeyword(NamelistElement("state", NamelistType::STRING, "MISSING",
                                   DefaultIsObligatory::NO, InputRepeats::YES));
  t_nml.addKeyword(NamelistElement("e_precision", NamelistType::INTEGER,
                                   std::to_string(default_energy_decimal_places)));
  t_nml.addKeyword(NamelistElement("sdf_item", { "-title", "-label", "-energy", "-message",
                                                 "-mask", "-parameter", "-typeI", "-typeJ",
                                                 "-typeK", "-typeL", "-typeM", "-exregno",
                                                 "-internal_regno", "-maccsid" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::INTEGER },
                                   { std::string(""), std::string("ALL"), std::string(""),
                                     std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string(""),
                                     std::string("OFF"), std::string("") },
                                   DefaultIsObligatory::NO, InputRepeats::YES, sdf_help,
                                   sdf_keys_help,
                                   { KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL }));
  t_nml.addKeyword(NamelistElement("common_path_limit", NamelistType::INTEGER,
                                   std::to_string(default_common_path_limit)));
  t_nml.addKeyword(NamelistElement("common_path_threshold", NamelistType::INTEGER,
                                   std::to_string(default_common_path_threshold)));
  t_nml.addKeyword(NamelistElement("outlier_sigmas", NamelistType::REAL,
                                   std::to_string(default_energy_outlier_sigmas)));
  t_nml.addKeyword(NamelistElement("outlier_count", NamelistType::INTEGER,
                                   std::to_string(default_outlier_limit)));
  t_nml.addHelp("syntax", "Layout of the diagnostics report file, intended to make it amenable "
                "to one of a variety of plotting programs for further analysis.  Options include "
                "MATLAB, OCTAVE, MATRIX_PKG (various matrix algebra programs), MATPLOTLIB or "
                "MAT_PLOT_LIB (python-based plotting library), or STANDALONE or STORMM (generic "
                "format).");
  t_nml.addHelp("scope", "The extent of reporting that shall take place for the energies and "
                "other properties of individual systems.");
  t_nml.addHelp("nrgsample", "The depth of reporting for energy and state variable diagnostics.  "
                "Choose TIME_SERIES, SERIES, or ALL for complete reporting of all measurements "
                "every ntpr steps.  Choose MEAN or AVERAGE to report just the averages and "
                "standard deviations of such quantities.  Choose FINAL or LAST to report only the "
                "values for the trajectory's final frame.");
  t_nml.addHelp("username", "Name of the user driving the run (if different from that which would "
                "be detected automatically).");
  t_nml.addHelp("timings", "By default, the wall time devoted to various aspects of a calculation "
                "will be displayed at the end of the run.  Set to ON or ACTIVE to ensure this "
                "behavior, or OFF to decline printed timings.");
  t_nml.addHelp("report_width", "Indicate the desired width of the report file.  For example, "
                "sander's mdout has a width of 80 characters.  This width will be respected "
                "except as required by aspects of line formatting and unbreakable variable or "
                "file names.");
  t_nml.addHelp("energy", "If unspecified, STORMM will print all relevant molecular mechanics "
                "energy components.  This keyword allows a user to select specific components for "
                "printing, to help focus the output and reduce file sizes.  Acceptable arguments "
                "include: BOND, ANGLE, TORSION / DIHEDRAL, PROPER_DIHEDRAL, IMPROPER_DIHEDRAL, "
                "CMAP, VDW / LJ, ELEC / ELECTROSTATIC, GB / SOLVENT, and NMR / RESTRAINT.  "
                "Various plurals or omissions of underscores may also be recognized, and the "
                "arguments are not case sensitive.");
  t_nml.addHelp("state", "If unspecified, the overall temperature will be printed.  If applicable "
                "and \"state\" is unspecified, volume and overall pressure will be printed.  A "
                "user may also opt to print temperatures of specific regions of the simulation by "
                "supplying the LOCAL_TEMPERATURE argument, and components of the virial using the "
                "VIRIAL_COMPONENTS argument.");
  t_nml.addHelp("e_precision", "The number of decimal places with which to report all energy "
                "quantities in the output tables.  Energies are reported in units of kcal/mol.  "
                "The default setting matches Amber output and is appropriate for most molecular "
                "simulations and geometry optimizations.");
  t_nml.addHelp("sdf_item", "Detail a data item to be included in an SD file archive.  These "
                "items can be attached to specific systems and display particular aspects of the "
                "energy or the model that calculated it.");
  t_nml.addHelp("common_path_limit", "The maximum number of common paths that will be used to "
                "condense output tables detailing the origins of each system in various files.");
  t_nml.addHelp("common_path_threshold", "The number of files that must use a common root path in "
                "order for it to be declared a common path and abstracted into a ${token} to "
                "condense output.");
  t_nml.addHelp("outlier_sigmas", "Multiplier for the number of standard deviations from the mean "
                "energy at which point a given eneregy (and the system displaying it) will be "
                "deemed an outlier.");
  t_nml.addHelp("outlier_count", "The maximum number of outliers that will be reported, per group "
                "if there is a grouping method in place or over the whole system otherwise.");

  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  
  return t_nml;
}

} // namespace namelist
} // namespace namelist
