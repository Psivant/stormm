// -*-c++-*-
#ifndef STORMM_NML_ANALYSIS_H
#define STORMM_NML_ANALYSIS_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/symbol_values.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/textfile.h"
#include "input.h"
#include "namelist_common.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using symbols::pi;
  
/// \brief Default values for analytic operations
/// \{
constexpr int default_analysis_stat_blocks = 1;           ///< Number of blocks into which results
                                                          ///<   will be divided for various
                                                          ///<   analyses.  The default behavior is
                                                          ///<   not to do block averaging.
constexpr int default_analysis_sample_interval = 1000;    ///< The default frequency at which to
                                                          ///<   perform general analysis actions
constexpr int default_analysis_init_step = 0;             ///< Initial step at which to begin
                                                          ///<   analysis actions, in general
constexpr double default_hbond_da_max_separation = 3.3;   ///< Maximum D - A distance to verify a
                                                          ///<  hydrogen bond
constexpr double default_hbond_dha_min_angle = pi * 0.5;  ///< Minimum D-H :: A angle to verify a
                                                          ///<   hydrogen bond
constexpr double default_hbond_initial_proximity = 20.0;  ///< Maximum initial distance between
                                                          ///<   donors and acceptors in order to
                                                          ///<   consider the pair as a possible
                                                          ///<   hydrogen bond in simulations
constexpr double default_hbond_occ_threshold = 0.05;      ///< The minimum threshold for which to
                                                          ///<   report a specific hydrogen bond
                                                          ///<   arrangement.  Hydrogen bonds with
                                                          ///<   lower occupancies will still be
                                                          ///<   tallied and enter into other
                                                          ///<   staistics collected for a given
                                                          ///<   system, but will not show up in a
                                                          ///<   listing of notable arrangements.
constexpr double default_hbond_act_threshold = 0.25;      ///< The minimum threshold for which to
                                                          ///<   report a system as one with
                                                          ///<   noteworthy hydrogen bonding
constexpr char default_analysis_identifier[] = "an";      ///< A simple string compatible with
                                                          ///<   integer extensions to give unique
                                                          ///<   tags to each of the analyses
                                                          ///<   enclosed within the namelist block
constexpr char default_hbond_variable_base[] = "hbana";   ///< Hydrogen bonding analysis report
                                                          ///<   variable base name
/// \}

using parse::WrapTextSearch;

/// \brief During molecular dynamics or other sequential calculations, it can be useful to take
///        customized measurements of the system in particular.  This can include, but is not
///        limited to, measurements of the distance between two particles or local angles in a
///        particular arrangement.
class AnalysisControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.  In general, the analysis namelist keywords will be
  ///        seeded with an input status of MISSING if they call for a specific analysis, or set to
  ///        default values if they provide specific parameters for carrying out a requested
  ///        activity.  Analysis of a simulation will be invoked by the presence of one of more
  ///        &analysis namelists
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &minimize namelist
  /// \param found_nml   Indicator of whether namelist input was found
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &minimize namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  AnalysisControls(ExceptionResponse policy_in = ExceptionResponse::DIE,
                   WrapTextSearch wrap = WrapTextSearch::NO);

  AnalysisControls(const TextFile &tf, int *start_line, bool *found_nml,
                   ExceptionResponse policy_in = ExceptionResponse::DIE,
                   WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another class object, on the right hand side of an assignment statement
  /// \{
  AnalysisControls(const AnalysisControls &original) = default;
  AnalysisControls(AnalysisControls &&original) = default;
  AnalysisControls& operator=(const AnalysisControls &other) = default;
  AnalysisControls& operator=(AnalysisControls &&other) = default;
  /// \}

  /// \brief Get the identifier attached to this analysis as a whole.  This is not a system label
  ///        matching something in the &files namelist (although it could overlap with something
  ///        from that set of names, if the user sees convenience).  This identifier will be
  ///        carried through to the output and placed in the preambles to results sections
  ///        displaying each analysis as well as to extensions of variable names within that
  ///        result if (and only if) multiple analyses are present in the input.  If a given
  ///        analysis namelist control block does not contain an identifier, a unique identifier
  ///        will be assigned at runtime.  If the analysis block contains multiple activities,
  ///        each will fall under the same general identifier, but with unique default names for
  ///        each type of analysis there will not be collisions in the space of variable names
  ///        found in the output.
  const std::string& getIdentifier() const;

  /// Indicate whether the set of activities prescribed in this object is the only set in the
  /// calculation (in the user's input), in which case the function returns FALSE, or if the
  /// object exists in the context of others (e.g. multiple &analysis namelist control blocks
  /// were present in the input).
  bool multipleAnalysesPresent() const;

  /// \brief Get the general block averaging directive, applicable to any statistical analysis
  ///        which doesn't have its own user-specified block averaging subdivision count.  To get
  ///        an authoritative count on a specific analysis, use its dedicated accessor.
  int getGeneralStatBlocks() const;

  /// \brief Get the sampling frequency (in terms of simulation time steps) for general analyses.
  ///        For the sampling frequency of a particular analysis, use its dedicated accessor.
  int getGeneralSamplingFrequency() const;

  /// \brief Get the simulation step at which to initiate general analyses.  For an authoritative
  ///        initiation point on a particular analysis, use its dedicated accessor.
  int getGeneralInitiationStep() const;
  
  /// \brief Get the total number of mask / system label pairs for identifying possible hydrogen
  ///        bond donors and acceptors among all systems in the calculation.
  int getHBondMaskCount() const;

  /// \brief Get the maximum separation at which a (proton) donor and acceptor can be considered
  ///        to form a hydrogen bond, units of in Angstroms.
  double getHBondMaxSeparation() const;

  /// \brief Get the minimum angle made by the arrangement of the proton donor, bonded hydrogen,
  ///        and proton acceptor in order for a hydrogen bond to be considered valid.
  double getHBondMinAngle() const;

  /// \brief Get the maximum separation between candidate donors and acceptors in the initial
  ///        configuration.  Donors and acceptors initially separated by more than this distance
  ///        will not be scanned for hydrogen bonding throughout the simulation.
  double getHBondCandidacyProximity() const;

  /// \brief Get the number of blocks to use in block averaging analysis results.
  int getHBondStatBlocks() const;

  /// \brief Get the sampling frequency, expressed in simulation time steps, at which to search for
  ///        hydrogen bonding.
  int getHBondSamplingFrequency() const;

  /// \brief Get the simulation step at which to begin searching for hydrogen bonds.
  int getHBondInitiationStep() const;
  
  /// \brief Get the minimum occupancy for reporting a specific hydrogen bond.
  double getHBondOccupancyThreshold() const;
  
  /// \brief Get the minimum hydrogen bonding propensity for reporting a specific hydrogen bond.
  double getHBondActivityThreshold() const;
  
  /// \brief Get an atom mask string specifying collections of potential proton donors and
  ///        acceptors.  The mask merely selects a group of atoms, from which actual donors and
  ///        acceptors will be identified.
  ///
  /// \param mask_index  Index of the mask, out of many that the user input may contain
  /// \param array_id    Index of the mask array, 0 for hbond_mask_i or 1 for hbond_mask_ii.  Any
  ///                    other value is invalid and will raise an exception.  This can be thought
  ///                    as the "dimension" of the hydrogen bonding matrix specified by a pair of
  ///                    atom masks.
  const std::string& getHBondMask(int mask_index, int array_id) const;

  /// \brief Get the system labels which put bounds on the systems that a particular atom mask can
  ///        be applied to when selecting a group of hydrogen bond donors and acceptors.
  ///
  /// \param label_index  Index of the label, out of many that the user input may contain
  const std::string& getHBondMaskLabel(int label_index) const;

  /// \brief Get the variable base name to be associated with hydrogen bond analysis.  This will
  ///        create multiple variables with different suffixes in the final report file, each a
  ///        matrix interpretable by the user's choice of third-party software package.  This will
  ///        return the variable base name specified by the user if there are no other &analysis
  ///        namelist control blocks, or a composite of the variable base name and the block's
  ///        (unique, unless duplicated by the user) identifier tag if multiple blocks are
  ///        involved.
  const std::string& getHBondVariableBase() const;

  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator& getTranscript() const;

  /// \brief Set the identifier for this analysis.
  ///
  /// \brief identifier_in  The new indetification label to set
  void setIdentifier(const std::string &identifier_in);

  /// \brief Set the presence of multiple analyses.  The default value allows the function to be
  ///        called without arguments in order to set the condition to TRUE.  The associated member
  ///        variable to this setter cannot be set through a namelist keyword.  Rather, the program
  ///        (the UserSettings class) will use this setter to mark the coexistence of multiple
  ///        &analysis namelists.
  ///
  /// \param found  Indicate whether multiple &analysis namelists are present
  void multipleAnalysesFound(bool found = true);
  
  /// \brief Set the general-purpose statistical block averaging subdivision count.  This will
  ///        apply to all analyses without their own user-specified block averaging.
  ///
  /// \param blocks_in  The number of blocks, taken from contiguous segments of the simulation,
  ///                   into which to subdivide the data set for block averaging
  void setGeneralStatBlocks(int blocks_in);

  /// \brief Set the general-purpose sampling frequency.  This will apply to all analyses without
  ///        their own user-specified smpling frequency.
  ///
  /// \param frequency_in  The frequency at which to sample the simulation trajectory, in units of
  ///                      time steps
  void setGeneralSamplingFrequency(int frequency_in);

  /// \brief Set the general-purpose initiation step for all analyses.  Like other general-purpose
  ///        settings, this can be overridden by specific directives in a particular analysis.
  ///
  /// \param step_in  The simulation time step at which to begin analysis
  void setGeneralInitiationStep(int step_in);
  
  /// \brief Set the masks and applicable system labels for a hydrogen bonding matrix.  This is the
  ///        equivalent setter for accessors getHBondMask() and getHBondMaskLabel(), in one.
  void setHBondMatrix(const std::string &mask_i, const std::string &mask_ii,
                      const std::string	&label);  
  
  /// \brief Set the maximum separation between hydrogen bond donor and acceptor atoms at which to
  ///        delcare a hydrogen bond to exist.
  ///
  /// \param separation_in  The maximum allowed separation to set
  void setHBondMaxSeparation(double separation_in);

  /// \brief Set the maximum separation between hydrogen bond donor and acceptor atoms at which to
  ///        delcare a hydrogen bond to exist.
  ///
  /// \param angle_in  The minimum angle to set
  void setHBondMinAngle(double angle_in);

  /// \brief Set the maximum separation between hydrogen bond donor and acceptor atoms at which to
  ///        delcare a hydrogen bond to exist.
  ///
  /// \param proximity_in  The requested maximum initial proximity
  void setHBondInitialProximity(double proximity_in);

  /// \brief Set the number of bins, collected from contiguous segments of the simulation, into
  ///        which to divide hydrogen bonding data for block averaging.
  ///
  /// \param blocks_in  The number of statistical blocks to set
  void setHBondStatBlocks(int blocks_in);

  /// \brief Set the sampling frequency for hydrogen bond analysis.
  ///
  /// \param frequency_in  The time step interval at which to conduct hydrogen bond analysis
  void setHBondSamplingFrequency(int frequency_in);

  /// \brief Set the first step of the simulation at which hydrogen bond analysis will begin.
  ///
  /// \param step_in  The requested initial step
  void setHBondInitiationStep(int step_in);

  /// \brief Set the minimum occupancy for reporting a specific hydrogen bond.
  ///
  /// \param occupancy_in  Occupancy required for a particular hydrogen bonding arrangement to be
  ///                      reported in the output
  void setHBondOccupancyThreshold(double occupancy_in);

  /// \brief Set the minimum hydrogen bonding activity for reporting a noteworthy system in the
  ///        output.
  ///
  /// \param activity_in  The required level of hydrogen bonding (sum of individual occupancies
  ///                     across all possible donor-acceptor pairs) in a given system
  void setHBondActivityThreshold(double activity_in);
  
  /// \brief Set the hydrogen bond variable base name for the subsequent report.
  ///
  /// \param varname_in  The requested variable base name
  void setHBondVariableBase(const std::string &varname_in);
  
private:

  /// Set the behavior when bad inputs are encountered.  DIE = abort program, WARN = warn the user
  /// and likely reset to the default value if one is available, SILENT = do not warn the user,
  /// but also likely reset to the default value if one is available.
  ExceptionResponse policy;

  /// A label to identify the analysis as a whole, including any activities it triggers.  If one
  /// &analysis namelist control block specifies hydrogen bonding as well as RMSD measurement,
  /// both will fall under the same identifier and the necessary distinctions in the results will
  /// be (must be) assured by unique variable base names.
  std::string identifier;

  // Counts of various analysis directives
  bool multiple_analyses_present;  ///< Flag to indicate the presence of multiple &analysis
                                   ///<   namelist blocks in the input, multiple collections of
                                   ///<   analysis procedures.  Each of the analysis blocks is
                                   ///<   expected to have a unique identifier tag, but some user
                                   ///<   input patterns could violate such an assumption.
  int general_statistical_blocks;  ///< The number of blocks to use in block averaging of any
                                   ///<   applicable analysis, barring a specific block count
                                   ///<   setting for that particular analysis
  int general_sample_frequency;    ///< The general sampling frequency applied to all analyses
                                   ///<   unless superceded by a specific frequency
  int general_initiation_step;     ///< The simulation step at which to initiate general analysis
                                   ///<   actions, superceded by specific directives in the manner
                                   ///<   that other general settings can be overridden
  int hbond_mask_count;            ///< The number of hydrogen bonding masks
  int hbond_statistical_blocks;    ///< Number of blocks to use in block averaging of hydrogen
                                   ///<   bonding results
  int hbond_sample_frequency;      ///< Sampling frequency applied to hydrogen bond analysis
  int hbond_initiation_step;       ///< Simulation step at which to initiate hydrogen bond analysis

  // Lists of mask strings critical to analysis directives
  std::vector<std::string> hbond_mask_i;       ///< Masks for selecting groups of atoms which may
                                               ///<   contain hydrogen bond donors and acceptors
  std::vector<std::string> hbond_mask_ii;      ///< Opposing masks for selecting groups of atoms
                                               ///<   which may contain hydrogen bond donors and
                                               ///<   acceptors.  Corresponding elements of
                                               ///<   hbond_mask_i and hbond_mask_ii will be
                                               ///<   applied to select atoms from the same system.
  std::vector<std::string> hbond_mask_labels;  ///< System labels indicating which systems in the
                                               ///<   calculation each corresponding mask in
                                               ///<   hbond_masks should apply to
  std::string hbond_variable_base;             ///< Base name of report / output tables associated
                                               ///<   with hydrogen bond analysis results
  std::string hbond_variable_composite;        ///< Base name of report / output tables associated
                                               ///<   with hydrogen bond analysis results when
                                               ///<   multiple &analysis namelist control blocks
                                               ///<   are in use

  // Other parameters for the analysis operations
  double hbond_max_separation;       ///< The maximum distance between donor and acceptor heavy
                                     ///<   atoms at which a hydrogen bond may be considered to
                                     ///<   exist
  double hbond_min_angle;            ///< The minimum angle formed by the donor, attached polar
                                     ///<   hydrogen, and acceptor atoms at which a hydrogen bond
                                     ///<   may be considered to exist
  double hbond_candidacy_proximity;  ///< Maximum distance between donor and acceptor atoms in the
                                     ///<   initial configuration by which it will be considered
                                     ///<   worthwhile to search for hydrogen bonds during the
                                     ///<   simulation
  double hbond_occupancy_threshold;  ///< Minimum occupancy for reporting a specific hydrogen
                                     ///<   bonding arrangement in the output
  double hbond_activity_threshold;   ///< Minimum sum of hydrogen bond occupancies across all
                                     ///<   possible partners within a system required for the
                                     ///<   system to be reported as significant in the analysis
  
  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;  

  /// \brief Validate a block averaging setting.
  ///
  /// \param setting  Pointer to the setting to validate.  This may be modified by the validator
  ///                 if found to contain spurious input.
  /// \param desc     Description of the particular setting to validate, for backtracing purposes
  void validateStatBlocks(int *setting, const char* desc);

  /// \brief Validate the sampling frequency for a named analysis.  Descriptions of input
  ///        parameters follow from validateStatBlocks(), above.
  void validateSamplingFrequency(int *setting, const char* desc);

  /// \brief Validate the initiation step for a named analysis.  Descriptions of input parameters
  ///        follow from validateStatBlocks(), above.
  void validateInitiationStep(int *setting, const char* desc);

  /// \brief Validate the distance between hydrogen bond donors and acceptors used to verify a
  ///        hydrogen bond.  Based on the policy, the a bogus distance may be replaced with the
  ///        default value.
  void validateHBondMaxSeparation();

  /// \brief Validate the minimum angle made by the donor, attached polar hydrogen, and acceptor
  ///        atoms needed to confirm the existence of a hydrogen bond.  Based on the policy, a
  ///        bogus angle may be replaced with the default value.
  void validateHBondMinAngle();

  /// \brief Validate the maximum distnace between donor and acceptor heavy atoms in the initial
  ///        configuration at which it will be deemed worthwhile to check the pair for hydrogen
  ///        bonding during a simulation.  Based on the policy, a bogus value may be replaced with
  ///        the default.
  void validateHBondCandidacyProximity();

  /// \brief Apply the general-purpose block averaging directive to any analysis which does not
  ///        already have its own setting.
  void applyGeneralStatBlocks();

  /// \brief Apply the general-purpose sampling frequency to any analysys which does not already
  ///        have its own setting.
  void applyGeneralSampleFrequency();

  /// \brief Apply the general-purpose initiation step to any analysis which does not already have
  ///        a user-specified initiation step.
  void applyGeneralInitiationStep();
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
NamelistEmulator analysisInput(const TextFile &tf, int *start_line, bool *found,
                               ExceptionResponse policy = ExceptionResponse::DIE,
                               WrapTextSearch wrap = WrapTextSearch::NO);

} // namespace namelist
} // namespace stormm

#endif
