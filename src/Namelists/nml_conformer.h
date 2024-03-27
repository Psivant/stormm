// -*-c++-*-
#ifndef STORMM_NML_CONFORMER_H
#define STORMM_NML_CONFORMER_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/symbol_values.h"
#include "Structure/structure_enumerators.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Parsing/textfile.h"
#include "namelist_element.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using constants::ExceptionResponse;
using parse::TextFile;
using parse::WrapTextSearch;
using structure::SamplingIntensity;
using synthesis::SystemGrouping;
using synthesis::VariableTorsionAdjustment;

/// \brief Default input settings for the &conformer namelist
/// \{
constexpr int default_conf_rotation_samples       = 3;
constexpr int default_conf_cis_trans_samples      = 2;
constexpr int default_conf_max_rotatable_bonds    = 4;
constexpr int default_conf_max_seeding_attempts   = 3;
constexpr int default_conf_clash_pairs            = 0;
constexpr int default_conf_max_system_trials      = 16384;
constexpr int default_conf_sample_trials          = 512;
constexpr int default_conf_running_states         = 16;
constexpr int default_conf_final_states           = 100;
constexpr int default_conf_reshuffle_iterations   = 0;
constexpr int active_states_limit                 = 524288;
constexpr int seeding_attempts_limit              = 100;
constexpr double default_conf_rmsd_tolerance      = 1.5;
constexpr double default_conf_core_restraint      = 16.0;
constexpr char default_conf_chirality[]           = "false";
constexpr char default_conf_cis_trans[]           = "false";
constexpr char default_conf_stop_hbonds[]         = "false";
constexpr char default_conf_output_grouping[]     = "system";
constexpr char default_conf_sampling_effort[]     = "light";
constexpr char default_conf_adjustment_method[]   = "none";
constexpr char default_conf_rotation_set_zero[]   = "60.0";
constexpr char default_conf_rotation_set_one[]    = "180.0";
constexpr char default_conf_rotation_set_two[]    = "-60.0";
constexpr char default_conf_cis_trans_set_zero[]  = "0.0";
constexpr char default_conf_cis_trans_set_one[]   = "180.0";
constexpr char default_conf_rotation_snap[]       = "10.0";
constexpr char default_conf_cis_trans_snap[]      = "5.0";
/// \}

/// \brief Object to encapsulate the data that can be extracted from the &conformer namelist.
class ConformerControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \param found_nml   Indication of whether the namelist was found in the input file
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &conformer namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  ConformerControls(ExceptionResponse policy_in = ExceptionResponse::DIE);
  
  ConformerControls(const TextFile &tf, int *start_line, bool *found_nml,
                    ExceptionResponse policy_in = ExceptionResponse::DIE,
                    WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  ConformerControls(const ConformerControls &original) = default;
  ConformerControls(ConformerControls &&original) = default;
  ConformerControls& operator=(const ConformerControls &original) = default;
  ConformerControls& operator=(ConformerControls &&original) = default;
  /// \}
  
  /// \brief Get the core atom mask string.
  const std::string& getCoreAtomMask() const;

  /// \brief Get the name of the data item expected to contain core atoms.
  const std::string& getCoreDataItemName() const;

  /// \brief Get the core restraint r2 penalty (near-repulsive potential) stiffness.
  double getCoreRK2Value() const;

  /// \brief Get the core restraint r3 penalty (position-keeping potential) stiffness.
  double getCoreRK3Value() const;

  /// \brief Get the core restraint r2 displacement.
  double getCoreR2Value() const;

  /// \brief Get the core restraint r3 displacement.
  double getCoreR3Value() const;

  /// \brief Get an indicator of whether to sample chirality.
  bool sampleChirality() const;

  /// \brief Get an indicator of whether to sample cis- and trans- isomers.
  bool sampleCisTrans() const;

  /// \brief Get an indicator as to whether to apply restraints that will prevent hydrogen bond
  ///        formation in the resulting conformers.
  bool preventHydrogenBonding() const;

  /// \brief Get the total number of states to attempt minimizing at one time.  This will put a
  ///        limit on the expanded population of conformer systems that the program will attempt
  ///        to model and minimize on the GPU, which takes a significant amount of memory.  If
  ///        there are more systems, the program will expand each and attempt minimizations as
  ///        this limit (a safeguard against overrunning available resources) permits.
  int getRunningStateCount() const;

  /// \brief Get the number of final states to produce for each initial system.
  int getFinalStateCount() const;

  /// \brief Get the number of samples to apply to each explicitly sampled rotatable bond.
  int getRotationSampleCount() const;

  /// \brief Get the number of samples to apply to each explicitly sampled cis-trans isomeric bond.
  int getCisTransSampleCount() const;

  /// \brief Get the list of rotational angle values to sample about each rotatable bond.
  const std::vector<double>& getRotationSampleValues() const;

  /// \brief Get the list of rotational angle values to sample about each cis-trans isomeric bond.
  const std::vector<double>& getCisTransSampleValues() const;

  /// \brief Get the maximum number of rotatable bonds to sample.
  int getRotatableBondLimit() const;

  /// \brief Get the maximum number of conformer seeding attempts
  int getMaxSeedingAttempts() const;

  /// \brief Get the number of clashing pairs of atoms that can be seen in a structure before it is
  ///        declared to contain too much internal conflict and a clash is reported.
  int getClashPairTolerance() const;
  
  /// \brief Get the maximum number of sampling trials that will be permitted across the various
  ///        systems when attempting to sample the viable rotameric states of each bond and other
  ///        possible isomerizations of each compound.  This applies on a per-system basis.
  int getSamplingTrialLimit() const;

  /// \brief Get the maximum number of minimizations to attempt with any one molecule.  Each
  ///        initial state provided by the user will be subject to this limit, so if the limit
  ///        is 5000 and one molecule has two initial states listed in the input deck, the total
  ///        number of conformations sampled will be no greater than 10000.
  int getMaximumTrialLimit() const;
  
  /// \brief Get the positional root mean squared deviation that will distinguish each reported
  ///        confomer.
  double getRMSDTolerance() const;

  /// \brief Get the output grouping strategy, indicating whether to take the specified number of
  ///        final states for each system entered by the user, for each label group spanning one
  ///        or more systems, or for each unique topology spanning one or more systems.
  SystemGrouping getGroupingMethod() const;

  /// \brief Get a general sampling strategy from the user.
  SamplingIntensity getSamplingIntensity() const;

  /// \brief Get the tolerance for snapping rotatable bond settings to consensus values.
  double getRotatableBondSnapThreshold() const;

  /// \brief Get the tolerance for snapping cis-trans isomeric bond settings to consensus values.
  double getCisTransBondSnapThreshold() const;

  /// \brief Get the torsion angle adjustment strategy (applies to rotatable and cis-trans
  ///        isomeric bonds) used in reconciling general settings to the known values emerging from
  ///        a collection of structures for each molecule.
  VariableTorsionAdjustment getTorsionAdjustmentProtocol() const;
  
  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator& getTranscript() const;
  
private:
  ExceptionResponse policy;         ///< Set the behavior when bad inputs are encountered.  DIE =
                                    ///<   abort program, WARN = warn the user, and likely reset to
                                    ///<   the default value if one is available, SILENT = do not
                                    ///<   warn the user, but also likely reset to the default
                                    ///<   value if one is available.
  std::string core_atom_mask;       ///< Mask of common core atoms, applied to all systems
  std::string core_data_item;       ///< Name of a data item from a Biovia SD file indicating the
                                    ///<   core atom mask.  Various formats of the input in the
                                    ///<   SD file will be automatically detected and accepted.
  double core_rk2;                  ///< Core atom repulsive positional restraint stiffness
  double core_rk3;                  ///< Core atom attractive positional restraint stiffness
  double core_r2;                   ///< Core atom maximum repulsive potential range
  double core_r3;                   ///< Range for onset of core atom attractive potential
  std::string anchor_conformation;  ///< The other half of the common core mask, the thing to align
                                    ///<   core atoms against.  This anchor conformation will be
                                    ///<   docked or placed within the receptor of interest and
                                    ///<   must contain atoms that register in the common atom
                                    ///<   mask.  All conformations will be aligned against this
                                    ///<   reference in order to bias the search towards results
                                    ///<   that fit inside of the receptor.
  bool sample_chirality;            ///< Flag to have chiral enantiomers sampled
  bool sample_cis_trans;            ///< Flag to have cis-trans isomers sampled
  bool prevent_hbonds;              ///< Flag to apply restraints that will prevent hydrogen bonds
                                    ///<   from forming during energy minimizations
  int running_states;               ///< Number of states to try minimizing at one time
  int final_states;                 ///< Number of final states to collect
  int rotation_sample_count;        ///< Number of times to sample about a rotatable bond
  int rotatable_bond_limit;         ///< Maximum number of rotatable bonds to explicitly sample
  int cis_trans_sample_count;       ///< Number of times to sample about a cis-trans isomeric bond
  int max_seeding_attempts;         ///< Maximum number of attempts to make in seeding each
                                    ///<   conformer.  If the seeding fails, the conformer will be
                                    ///<   left in its original state from the user input.
  int clash_pair_tolerance;         ///< The number of clashing pairs that will be permitted in any
                                    ///<   given structure before the internal conflicts force the
                                    ///<   program to randomize whatever components it can and
                                    ///<   reassess the structure.
  int sampling_trial_limit;         ///< The maximum number of trials per system to run for rotamer
                                    ///<   and cis-trans bond settings exploration, which will be
                                    ///<   interpreted in light of the total number of systems in
                                    ///<   the calculation.
  int maximum_trial_limit;          ///< Maximum number of distinct minimizations to attempt with
                                    ///<   one molecule
  double rmsd_tolerance;            ///< Minimum mass-weighted root-mean squared deviation between
                                    ///<   unique conformers
  std::string group_method;         ///< String indicating the way in which to group trials from
                                    ///<   one or more user-specified system inputs when selecting
                                    ///<   output conformations
  std::string sample_effort;        ///< String indicating the degree of effort to exert on any
                                    ///<   particular compound based on the number or combinations
                                    ///<   of rotatable bonds, chiral centers, or isomeriable
                                    ///<   cis-trans bonds
  double rotation_snap_threshold;   ///< The threshold below which rotatable bond angle settings
                                    ///<   are to be snapped to nearby consensus values found in a
                                    ///<   collection of conformations.
  double cis_trans_snap_threshold;  ///< The threshold below which cis-trans isomeric bond angle
                                    ///<   settings are to be snapped to nearby consensus values
                                    ///<   found in a collection of conformations.
  std::string adjustment_method;    ///< The method for adjusting torsion angle settings (whether
                                    ///<   for rotatable or cis-trans isomeric bond settings)

  /// The specific values at which to set each rotatable bond for initial poses.
  std::vector<double> rotation_sample_values;

  /// The specific values at which to set each cis-trans isomeric bond for initial poses.
  std::vector<double> cis_trans_sample_values;

  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;
  
  /// \brief Validate the restraining potential that will define the common core.
  void validateCoreRestraint() const;
  
  /// \brief Validate input pertaining to chiral sampling.
  ///
  /// \param directive  The keyword setting for chirality sampling (must be 'true' or 'false',
  ///                   without case sensitivity)
  void validateSampleChirality(const std::string &directive) const;

  /// \brief Validate input pertaining to sampling cis- and trans- states of molecules.
  ///
  /// \param directive  The keyword setting for cis- and trans- sampling (must be 'true' or
  ///                   'false', without case sensitivity)
  void validateSampleCisTrans(const std::string &directive) const;

  /// \brief Validate input pertaining to hydrogen bonding prevention.
  ///
  /// \param directive  The keyword setting for cis- and trans- sampling (must be 'true' or
  ///                   'false', without case sensitivity)
  void validatePreventHBonds(const std::string &directive) const;  
  
  /// \brief Validate the replica counts and criteria for distinguishing unique conformers.
  void validateStateCounts();

  /// \brief Validate the manner in which to group systems for analysis purposes.
  void validateGroupingMethod();

  /// \brief Validate the number of seeding attempts to be made with each conformer.
  void validateSeedingAttempts();

  /// \brief Validate the number of clashing atom pairs that will be tolerated before a seeded
  ///        confomation is deemed unfixable and must be re-seeded before further refinement can
  ///        take place.
  void validateClashCounts();
  
  /// \brief Validate the level of effort to be invested in sampling any given molecule.
  void validateSamplingIntensity();

  /// \brief Validate the snapping threshold for either rotatable bond or cis-trans isomeric bond
  ///        settings.
  ///
  /// \param snap_setting  The snapping threshold to validate, ensuring that it lies in the range
  ///                      [ 0, 2pi ).  May be modified and returned if the setting is invalid and
  ///                      the exception handling permits. 
  /// \param desc          A description of the threshold to validate
  void validateTorsionSnapping(double *snap_setting, const std::string &desc) const;

  /// \brief Validate the string used to specify the torsion adjustment protocol.
  void validateTorsionAdjustmentProtocol();
  
  /// \brief Process values for setting the angle about each rotatable bond or cis-trans isomeric
  ///        bond from the namelist object.  Store the values in the appropriate vector and
  ///        length counter.  Length counters are provided separately to give the user more options
  ///        as to how to specify the array of angles.
  ///
  /// \param sample_count   The number of rotatory samples to attempt.  Set and returned.
  /// \param count_keyword  Keyword associated with the sample_count integer
  /// \param sample_values  Array of values for all rotation angles to be sampled.  Filled and
  ///                       returned.
  /// \param value_keyword  Keyword associated with the sample_values array
  /// \param value_stride   The amount to stride forward when sampling values around the rotation
  ///                       circle with automated sample delineation, given in units of degrees.
  ///                       120 for rotatable bonds and 180 for cis-trans isomeric bonds.
  /// \param value_notch    The amount by which to perturb the circular orbit with each successive
  ///                       rotation, to prevent successive orbits from creating the same set of
  ///                       samples.  10 for rotatable bonds and 5 for cis-trans isomeric bonds.  
  /// \param t_nml          Contains user input from a &conformer namelist
  void processSamplingValues(int *sample_count, const std::string &count_keyword,
                             std::vector<double> *sample_values, const std::string &value_keyword,
                             double value_stride, double value_notch,
                             const NamelistEmulator &t_nml);
};

/// \brief Free function to read the &conformer namelist.  
///
/// \param tf          Text of file containing the input deck, read into RAM
/// \param start_line  Line of the input file at which to begin the scan
/// \param found       Indicator that the namelist was found in the input file
/// \param policy      Response to bad inputs
/// \param wrap        Indicate that the search for a &conformer namelist should carry on from the
///                    beginning of an input file if no such namelist is found starting from the
///                    original starting point
NamelistEmulator conformerInput(const TextFile &tf, int *start_line, bool *found,
                                ExceptionResponse policy = ExceptionResponse::DIE,
                                WrapTextSearch wrap = WrapTextSearch::NO);
  
} // namespace namelist
} // namespace stormm

#endif
