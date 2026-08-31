// -*-c++-*-
#ifndef STORMM_NML_DEBUG_H
#define STORMM_NML_DEBUG_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/textfile.h"
#include "input.h"
#include "namelist_common.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

/// \brief Default values for debugging
/// \{
constexpr int default_inspection_interval = 1000;
constexpr int default_max_anomaly_reports = 1024;
constexpr double default_large_force_threshold = 512.0;
constexpr double default_high_speed_threshold = 0.5;
constexpr double default_bond_max_strain = 160.0;
constexpr double default_angle_max_strain = 32.0;
constexpr char default_debug_variable_prefix[] = "dbg";
/// \}
  
using parse::WrapTextSearch;

/// \brief Object to encapsulate debugging control information.  Like other namelist encapsualtors,
///        this object can take input file data as part of its construction, or by a series of
///        setters.  Validation of each piece of data is handled as it appears either in the
///        contructor or via setters.  Getter functions dispense the internal information to any
///        application using STORMM libraries.
class DebugControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.  In general, the debugging namelist keywords will be
  ///        seeded with an input status of MISSING.  Debugging activities will be part of any
  ///        program that utilizes a global object of some intervention class (e.g.,
  ///        DynamicsIntervention), but not invoked unless some keyword comes with a USER_SPECIFIED
  ///        status.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &minimize namelist
  /// \param found_nml   Indicator of whether namelist input was found
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &minimize namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  DebugControls(ExceptionResponse policy_in = ExceptionResponse::DIE,
                WrapTextSearch wrap = WrapTextSearch::NO);

  DebugControls(const TextFile &tf, int *start_line, bool *found_nml,
                ExceptionResponse policy_in = ExceptionResponse::DIE,
                WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another class object, on the right hand side of an assignment statement
  /// \{
  DebugControls(const DebugControls &original) = default;
  DebugControls(DebugControls &&original) = default;
  DebugControls& operator=(const DebugControls &other) = default;
  DebugControls& operator=(DebugControls &&other) = default;
  /// \}

  /// \brief Check the combined forces acting on each particle.  This is applicable in the context
  ///        of non-periodic simulations, where all forces are tallied in the same buffers before
  ///        the debugging can intervene at a point accessible to other interventions.
  bool checkForces() const;
  
  /// \brief Get a directive on whether to report unusually large forces.  This check will apply
  ///        to forces due to interactions in the neighbor list, if applicable, from non-bonded
  ///        tiles, and from local interactions.  Each source of forces will be checked separately.
  bool reportLargeForces() const;

  /// \brief Get a directive on whether to perform checks on the neighbor list composition,
  ///        including atom placement and cell population counts.
  bool checkNeighborListComp() const;

  /// \brief Get a directive on whether to keep a record of momentum purges.
  bool trackMomentumPurge() const;
  
  /// \brief Get the interval at which all active checks (forces, neighbor list composition, and
  ///        more) are obligatory.  If large forces or unconverged constraint adjustments are to
  ///        be reported, active checks will also be performed when they are encountered.
  int getInspectionInterval() const;

  /// \brief Get the maximum number of reports (whether of large forces, high particle speeds,
  ///        failed convergence) to log.  Events beyond this limit will no longer be recorded for
  ///        making the final report, although they could still trigger force checks.
  int getMaximumReports() const;

  /// \brief Get the threshold at which forces acting on any given atom can be considered "large."
  ///        This value is not scaled according to the mass of the atom.
  double getLargeForceThreshold() const;

  /// \brief Get the threshold at which speeds of any particular atom can be considered "high."
  double getHighSpeedThreshold() const;

  /// \brief Indicate whether to run (and act upon) sanity checks during calculations.  The nature
  ///        of the checks may depend on the program implementing them.  For typical molecular
  ///        dynamics, this will mean inspecting systems after energy minimization is complete.
  ///        The rebuilding the synthesis after removing any which have very high bond or angle
  ///        strain in one or more particular terms.
  bool runSanityChecks() const;

  /// \brief Get the maximum strain that any particular bond may experience before a structure is
  ///        deemed to fail a sanity check.
  double getBondStrainTolerance() const;

  /// \brief Get the maximum strain that any particular bond angle may experience before a
  ///        structure is deemed to fail a sanity check.
  double getAngleStrainTolerance() const;

  /// \brief Get the variable name to prefix tables reporting debugging data.
  const std::string& getVariableBase() const;

  /// \brief Toggle the setting on whether to implement a CPU check of total forces acting on each
  ///        particle.  Descriptions of input variables follow from setNeighborListForceCheck(),
  ///        above.
  void setForceCheck(bool setting_in);

  /// \brief Toggle the setting on whether to implement a check for large magnitude forces
  ///        appearing in the simulation.  To set the force_threshold keyword, whether via user
  ///        input or setLargeForceThreshold, will also toggle the check to TRUE (reporting for
  ///        large forces will be ON).  Descriptions of input variables follow from
  ///        setNeighborListForceCheck(), above.
  void setLargeForceReport(bool setting_in);

  /// \brief Toggle the setting on whether to implement a check on the neighbor list composition.
  ///
  /// \param setting_in  Indicate whether to change the setting to TRUE or FALSE
  void setNeighborListCompCheck(bool setting_in);

  /// \brief Toggle the setting on whether to record net translational momentum purges (and
  ///        rotational purges, if applicable).
  ///
  /// \param setting_in  Indicate whether to change the setting to TRUE or FALSE
  void setMomentumPurgeTracking(bool setting_in);

  /// \brief Set the maximum number of anomalous incident reports.
  ///
  /// \param maximum_reports_in  The maximum number of reports to allocate for
  void setMaximumReports(int maximum_reports_in);
  
  /// \brief Set the step interval (e.g. number of minimization cycles or number of dynamics time
  ///        steps) at which all active checks become obligatory, even if no violations have been
  ///        detected.  A value of zero or a negative value will indicate that there is no
  ///        obligatory, periodic enforcement of active checks.
  ///
  /// \param interval_in  The requested interval
  void setInspectionInterval(int interval_in);

  /// \brief Set the threshold for determining the magnitude of the force acting on some atom to
  ///        be "large."
  ///
  /// \param threshold_in  The requested threshold, in units of kcal/mol-Angstrom
  void setLargeForceThreshold(double threshold_in);

  /// \brief Set the threshold for determining that the speed of some particle is "large."
  ///
  /// \param threshold_in  The requested threshold, in units of Angstrom / fs
  void setHighSpeedThreshold(double threshold_in);

  /// \brief Toggle the setting on whether to run sanity checks on systems over the course of a
  ///        calculation.  Actual implementations of the checks, if the setting is TRUE, will vary
  ///        based on the program.
  ///
  /// \param setting_in  Indicate whether to change the setting to TRUE or FALSE
  void setRunSanityChecks(bool setting_in);

  /// \brief Set the prefix for debugging-related report variables.
  ///
  /// \param varname_in  The preferred variable prefix
  void setVariableBase(const std::string &varname_in);
  
private:

  /// Set the behavior when bad inputs are encountered.  DIE = abort program, WARN = warn the user
  /// and likely reset to the default value if one is available, SILENT = do not warn the user,
  /// but also likely reset to the default value if one is available.
  ExceptionResponse policy;

  // Triggers for most debugging activities are controlled by boolean member variables.
  bool check_forces;              ///< Inspect general forces computed in non-periodic simulations
  bool report_large_forces;       ///< Report forces on any particle exceeding the magnitude
                                  ///<   specified in large_force_threshold
  bool check_neighbor_list;       ///< Check the neighbor list placement of particles.  This is
                                  ///<   independent of the check on neighbor list forces.
  bool track_momentum_purge;      ///< Keep records of the amount of momentum purged from the
                                  ///<   system at each interval triggered by the &dynamics
                                  ///<   namelist nscm keyword
  bool enforce_sanity;            ///< Run sanity checks, the nature of which depend on the program

  /// The maximum number of reports to store for anomalous forces, speeds, and other events
  /// recorded over the course of a calculation.  The lists compiled in this manner, up to the
  /// stated maximum number of events, may be used in a post-mortem report.
  int maximum_reports;
  
  /// All checks (forces, neighbor list placement) may be implemented at periodic intervals as
  /// opposed to every step.
  int inspection_interval;
  
  /// The threshold at which forces will be reported.  This can also serve as the threshold at
  /// which a more laborious CPU-based check on the forces will be triggered.  Any force with an
  /// overall magnitude greater than this value will trigger reporting or whatever applicable
  /// checks have been specified.
  double large_force_threshold;

  /// The threshold at which particle speeds will be reported. (Velocity is a vector, speed is its
  /// overall magnitude.)
  double high_speed_threshold;

  // Sanity metrics for various checks
  double bond_strain_sanity;   ///< The tolerated bond strain in any single term
  double angle_strain_sanity;  ///< The tolerated angle strain in any single term

  /// Base name of the variable to which debugging results are printed.  These variables (this
  /// prefix plus descriptive extensions) make it possible to load debugging results into third
  /// party programs such as matrix algebra packages.  See the &report namelist for other
  /// controls on the output.
  std::string varname;
  
  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;
};

/// \brief Produce a namelist for specifying debugging directives.
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
NamelistEmulator debugInput(const TextFile &tf, int *start_line, bool *found,
                            ExceptionResponse policy = ExceptionResponse::DIE,
                            WrapTextSearch wrap = WrapTextSearch::NO);

} // namespace namelist
} // namespace stormm

#endif
