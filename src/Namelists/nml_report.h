// -*-c++-*-
#ifndef STORMM_NML_OUTPUT_H
#define STORMM_NML_OUTPUT_H

#include <cstring>
#include <string>
#include <vector>
#include <cstdio>
#include <sys/time.h>
#include "copyright.h"
#include "Constants/behavior.h"
#include "DataTypes/common_types.h"
#include "MoleculeFormat/mdlmol_dataitem.h"
#include "Parsing/textfile.h"
#include "Potential/energy_enumerators.h"
#include "Reporting/reporting_enumerators.h"
#include "input.h"
#include "namelist_element.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using energy::EnergySample;
using energy::StateVariable;
using parse::TextFile;
using parse::WrapTextSearch;
using review::OutputScope;
using review::OutputSyntax;
using structure::MdlMolDataRequest;

/// \brief Default values for the &report namelist
/// \{
constexpr int default_common_path_limit = 8;
constexpr int default_common_path_threshold = 4;
constexpr int default_energy_decimal_places = 4;
constexpr double default_energy_outlier_sigmas = 3.0;
constexpr int default_outlier_limit = 4;
constexpr char default_report_variable_name[] = "stormm";
/// \}

/// \brief Collect output directives relating to the diagnostics file.  While the output frequency
///        is controlled by another namelist such as &minimize or &dynamics, STORMM is designed to
///        handle many trajectories and thus the user is faced with a choice of what information
///        to present and in what format.
class ReportControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &minimize namelist
  /// \param found_nml   Indicator of whether namelist input was found
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &minimize namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  ReportControls(ExceptionResponse policy_in = ExceptionResponse::DIE,
                 WrapTextSearch wrap = WrapTextSearch::NO);
  ReportControls(const TextFile &tf, int *start_line, bool *found_nml,
                 ExceptionResponse policy_in = ExceptionResponse::DIE,
                 WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  ReportControls(const ReportControls &original) = default;
  ReportControls(ReportControls &&original) = default;
  ReportControls& operator=(const ReportControls &original) = default;
  ReportControls& operator=(ReportControls &&original) = default;
  /// \}

  /// \brief Get the output format--what program will be able to read the results?
  OutputSyntax getOutputSyntax() const;

  /// \brief Get the scope of results: will all systems' results be displayed, or some statistical
  ///        reduction?
  OutputScope getOutputScope() const;

  /// \brief Get the type of energy or state variable sampling (mean, time series, final value).
  EnergySample getEnergySamplingMethod() const;

  /// \brief Get the method for grouping systems.
  
  /// \brief Get the name of the user responsible for generating the results.
  const std::string& getUsername() const;

  /// \brief Get the (base) name of the matrix variable under which to store results.
  const std::string& getReportVariable() const;
  
  /// \brief Get the date on which the program began to run, taken as the time that this
  ///        ReportControls object was constructed.
  timeval getStartDate() const;
  
  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator& getTranscript() const;
  
  /// \brief Indicate whether to print timings data from the calculations.
  bool printWallTimeData() const;

  /// \brief Get the requested width for the output file.
  int getReportFileWidth() const;
  
  /// \brief Get the number of reported quantities.
  int getReportedQuantityCount() const;
  
  /// \brief Get a list of the energy components to display.  The actual results will depend on
  ///        whether each energy term is applicable, i.e. kinetic energy will not be displayed in
  ///        a conjugate gradient energy minimization run and the pressure will not be displayed
  ///        for a system run in isolated boundary conditions.
  const std::vector<StateVariable>& getReportedQuantities() const;

  /// \brief Get the number of data item requests for printing to an SD file.
  int getSDFileDataRequestCount() const;
  
  /// \brief Get a list of data items to be displayed in an SD file.  Reported quantities can
  ///        include topological parameters, components of the molecular mechanics energy for the
  ///        structure, the username, and notes about external factors affecting the final
  ///        conformation.
  const std::vector<MdlMolDataRequest>& getSDFileDataRequests() const;

  /// \brief Get a specific data request to be displayed in an SD file.
  ///
  /// \param index  Index of the request to produce
  MdlMolDataRequest getSDFileDataRequest(int index) const;

  /// \brief Get the maximum number of common paths to be permitted in system output.
  int getCommonPathLimit() const;

  /// \brief Get the required number of times a root path must be repeated in order to become one
  ///        of the common paths used to condense output tables.
  int getCommonPathThreshold() const;

  /// \brief Get the number of decimal places with which to express energetic quantities in output
  ///        tables.  This precision will be consistent across all energy quantities.
  int getEnergyDecimalPlaces() const;

  /// \brief Get the multiplier for the number of standard deviations from the mean that defines
  ///        an outlier among the energetics of a group of systems.
  double getOutlierSigmaFactor() const;

  /// \brief Get the maximum number of outliers to report (per group, or per the whole synthesis
  ///        if there is no grouping) in terms of systems deviating from the expected (mean)
  ///        energy.
  int getOutlierCount() const;
  
  /// \brief Set the output format.
  ///
  /// Overloaded:
  ///   - Supply the enumeration directly
  ///   - Supply a string that can be translated into the enumeration
  ///
  /// \param report_layout_in  The chosen layout for the output report, setting up for parsing by
  ///                          one of several plotting programs that can be used in analysis
  /// \{
  void setOutputSyntax(OutputSyntax report_layout_in);
  void setOutputSyntax(const std::string &report_layout_in);
  /// \}

  /// \brief Set the scope of the output.
  ///
  /// Overloaded:
  ///   - Supply the enumeration directly
  ///   - Supply a string that can be translated into the enumeration
  ///
  /// \param report_scope_in  The breadth of reporting for energies and the outcomes of other state
  ///                         variables across all systems in the calculation
  /// \{
  void setOutputScope(OutputScope report_scope_in);
  void setOutputScope(const std::string &report_scope_in);
  /// \}

  /// \brief Set the method for reporting samples of system state variables.
  ///
  /// Overloaded:
  ///   - Supply the enumeration directly.
  ///   - Supply a string that can be translate into the enumeration.
  ///
  /// \param state_sampling_in  The sampling method to set
  /// \{
  void setStateSampling(EnergySample state_sampling_in);
  void setStateSampling(const std::string &state_sampling_in);
  /// \}
  
  /// \brief Set the username.
  ///
  /// \param username_in  A username reported in the input file, not detected automatically from
  ///                     the environment (as is the default behavior)
  void setUsername(const std::string &username_in);

  /// \brief Set the name of the matrix variable for reporting state quantities (volume,
  ///        temperature, energy components).
  ///
  /// \param report_variable_in  The base name of the matrix variable to use
  void setReportVariable(const std::string &report_variable_in);
  
  /// \brief Set the preference for printing wall clock timings.
  ///
  /// Overloaded:
  ///   - Accept a single state variable to report
  ///   - Accept a vector of state variables to report
  ///
  /// \param preference  Set to TRUE to activate printing of wall time results
  /// \{
  void setWallTimeData(bool preference);
  void setWallTimeData(const std::string &preference);
  /// \}

  /// \brief Set the report file width.
  ///
  /// \param report_file_width_in  The general width (in ASCII characters) of report output
  void setReportFileWidth(int report_file_width_in);
  
  /// \brief Set the reported quantities that will be reported.  Some quantities are obligatory:
  ///        Total, total potential, and total kinetic energy; thermodynamic integration energy
  ///        derivatives (if applicable); pressure and volume (if applicable); and overall system
  ///        temperature.  By default, all components of the molecular mechanics energy are also
  ///        recorded, but requesting specific components through this function will disable all
  ///        molecular mechanics energy components not specifically requested.
  ///
  /// Overloaded:
  ///   - Accept a string to translate into some state variable code
  ///   - Accept a single state variable to report
  ///   - Accept a vector of state variables to report
  ///
  /// \param quantity_in    String indicating some state variable
  /// \param quantities_in  A state variable quantity, or list of quantities, to include in the
  ///                       diagnostics report file
  /// \{
  void setReportedQuantities(const std::string &quantity_in);
  void setReportedQuantities(StateVariable quantities_in);
  void setReportedQuantities(const std::vector<StateVariable> &quantities_in);
  /// \}
  
  /// \brief Add a request for a data item in an output SD file.
  ///
  /// \param ask  The requested data item
  void addDataItem(const MdlMolDataRequest &ask);

  /// \brief Set the upper limit on the number of common paths used to condense output files.
  ///
  /// \param common_path_limit_in  The maximum number of path shortcuts to create
  void setCommonPathLimit(int common_path_limit_in);
  
  /// \brief Set the number of instances that a common root path must appear in order for it to be
  ///        declared a common path with a token variable to condense output tables.
  ///
  /// \param common_path_threshold  The number of path tree instances needed to create a shortcut
  void setCommonPathThreshold(int common_path_threshold);

  /// \brief Set the number of digits after the decimal with which to express energetic quantities.
  ///
  /// \param energy_decimal_places_in  The number of digits after each energy quantity decimal
  void setEnergyDecimalPlaces(int energy_decimal_places_in);
  
private:
  ExceptionResponse policy;     ///< Course of action if bad input is encountered
  OutputSyntax report_layout;   ///< Layout of the diagnostics report file, making it amenable to
                                ///<   one or more plotting or matrix analysis programs
  OutputScope report_scope;     ///< Detail in which to report each system's results
  EnergySample state_sampling;  ///< The method for sampling and reporting state variables in the
                                ///<   system.
  std::string username;         ///< Name of the user (detected automatically)
  std::string report_variable;  ///< Name of the variable under which to file results for
                                ///<   simulation diagnostics (this is appended with various state
                                ///<   variable names, as appropriate)
  timeval start_date;           ///< Date and time on which this object was constructed, which will
                                ///<   be taken as a good approximation of the date and time on
                                ///<   which the calling program was executed.
  bool print_walltime_data;     ///< Flag to activate printing of wall time data (collected from
                                ///<   functions that issue results) to profile the calculation
  int report_file_width;        ///< Indicate the width that a report file should keep to, within
                                ///<   the limits of required formatting an variable or file names
                                ///<   that cannot be broken up.
  int common_path_limit;        ///< The maximum number of distinct common paths to find when
                                ///<   seeking ways to condense tables of file names for output
  int common_path_threshold;    ///< Threshold for the number of repeated instances of a common
                                ///<   path needed to declare a new path shortcut when seeking ways
                                ///<   to condense tables of file names for output
  int energy_decimal_places;    ///< The number of places with which to express energy quantities
                                ///<   in output tables
  double outlier_sigma_factor;  ///< The multiple of sigmas from the mean energy at which a system
                                ///<   might qualify as an outlier
  int outlier_count;            ///< The maxmum number of outliers to report, whether across the
                                ///<   entire set or within each partition of the set
  
  /// List the energetic quantities to report
  std::vector<StateVariable> reported_quantities;

  /// List of data item requests for information to be included in an SD file output
  std::vector<MdlMolDataRequest> sdf_addons;

  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;
  
  /// \brief Translate a user-supplied string into one or more energy quantities classifiable with
  ///        the state variable.
  ///
  /// \param inpstr  The user input
  std::vector<StateVariable> translateEnergyComponent(const std::string &inpstr);

  /// \brief Translate a user-supplied string into one or more state variables (pressure,
  ///        temperature, volume, or subclassefs of each).
  ///
  /// \param inpstr  The user input
  std::vector<StateVariable> translateStateQuantity(const std::string &inpstr);
  
  /// \brief Translate the namelist STRUCT input for an SD file data item into the internal object
  ///        encoding such a request.
  ///
  /// \param t_nml  The original namelist object, complete with keywords and values read from the
  ///               input file
  /// \param index  Index of the STRUCT keyword value to translate into a data item request
  std::vector<MdlMolDataRequest> translateSdfKeywordInput(const NamelistEmulator &t_nml,
                                                          int index);

  /// \brief Validate the number of decimal places chosen for energetic quantities.  The number
  ///        must be at least 1 (to ensure results are useful) and at most twelve (the practical
  ///        limit of double-precision calculations under even the most favorable circumstances).
  void validateEnergyDecimalPlaces();

  /// \brief Validate the outlier specifications, both in number and in degree.
  void validateOutlierMetrics();
};

/// \brief Produce a namelist for specifying content and format of the output diagnostics report,
///        similar in purpose to sander and pmemd's mdout.  This level of control is a unique
///        feature that STORMM apps can use to focus the reporting on results that are most
///        important to a user.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will wrap back to the beginning of the TextFile object, if needed,
///                    to find a &minimize namelist)
/// \param found       Indicate that the namelist was found
/// \param policy      Reaction to exceptions encountered during namelist reading
/// \param wrap        Indicate that the search for an &minimize namelist should carry on from the
///                    beginning of an input file if no such namelist is found starting from the
///                    original starting point
NamelistEmulator reportInput(const TextFile &tf, int *start_line, bool *found,
                             ExceptionResponse policy = ExceptionResponse::DIE,
                             WrapTextSearch wrap = WrapTextSearch::NO);  

} // namespace namelist
} // namespace stormm

#endif
