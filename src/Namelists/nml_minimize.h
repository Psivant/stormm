// -*-c++-*-
#ifndef STORMM_NML_MINIMIZE_H
#define STORMM_NML_MINIMIZE_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/textfile.h"
#include "input.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using parse::WrapTextSearch;

/// \brief Default values for energy minimization
/// \{
constexpr int default_minimize_maxcyc         = 200;
constexpr int default_minimize_ncyc           = 50;
constexpr int default_minimize_cdcyc          = 25;
constexpr int default_minimize_ntpr           = 50;
constexpr char default_minimize_checkpoint[]  = "true";
constexpr double default_minimize_cut         = 8.0;
constexpr double default_minimize_dx0         = 0.01;
constexpr double default_minimize_drms        = 0.0001;
constexpr double default_minimize_clash_r0    = 0.2;
constexpr double default_minimize_clash_ratio = 0.5;
/// \}
  
/// \brief Object to encapsulate energy minimization control information.  Like other namelist
///        encapsualtors, this object can take input file data as part of its construction, or
///        by a series of setters.  Validation of each piece of data is handled as it appears
///        either in the contructor or via setters.  Getter functions dispense the internal
///        information to any application using STORMM libraries.
class MinimizeControls {
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
  MinimizeControls(ExceptionResponse policy_in = ExceptionResponse::DIE,
                   WrapTextSearch wrap = WrapTextSearch::NO);
  MinimizeControls(const TextFile &tf, int *start_line, bool *found_nml,
                   ExceptionResponse policy_in = ExceptionResponse::DIE,
                   WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  MinimizeControls(const MinimizeControls &original) = default;
  MinimizeControls(MinimizeControls &&original) = default;
  MinimizeControls& operator=(const MinimizeControls &original) = default;
  MinimizeControls& operator=(MinimizeControls &&original) = default;
  /// \}

  /// \brief Get the total number of minimization cycles.
  int getTotalCycles() const;

  /// \brief Get the number of steepest descent cycles.
  int getSteepestDescentCycles() const;

  /// \brief Get the number of clash-relaxation cycles.
  int getClashDampingCycles() const;
  
  /// \brief Get the diagnostic output printing frequency, akin to the major contribution to pmemd
  ///        and sander mdout files.
  int getDiagnosticPrintFrequency() const;

  /// \brief Get the directive on whether to produce a checkpoint file for the final state of each
  ///        energy minimization run.
  bool getCheckpointProduction() const;
  
  /// \brief Get the electrostatic cutoff.
  double getElectrostaticCutoff() const;
  
  /// \brief Get the Lennard-Jones cutoff.
  double getLennardJonesCutoff() const;
  
  /// \brief Get the initial step length.
  double getInitialStep() const;

  /// \brief Get the convergence criterion.
  double getConvergenceTarget() const;

  /// \brief Get the absolute clash distance, below which two particles will be deemed colliding.
  double getAbsoluteClashDistance() const;

  /// \brief Get the minimum ratio of inter-particle distance to the pairwise van-der Waals sigma
  ///        parameter, below which two particles will be deemed colliding.
  double getVdwClashRatio() const;

  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator& getTranscript() const;
  
  /// \brief Set the total number of minimization cycles.
  ///
  /// \param cycles_in  The requested number of minimization cycles
  void setTotalCycles(int cycles_in);
  
  /// \brief Set the number of steepest descent cycles.
  ///
  /// \param cycles_in  The requested number of steepest descent cycles
  void setSteepestDescentCycles(int cycles_in);

  /// \brief Set the number of clash damping cycles.
  ///
  /// \param cycles_in  The requested number of clash damping cycles
  void setClashDampingCycles(int cycles_in);

  /// \brief Set the minimum inter-particle separation defining an electrostatic clash.
  ///
  /// \param clash_minimum_distance_in
  void setAbsoluteClashDistance(double clash_minimum_distance_in);

  /// \brief Set the minimum ratio of inter-particle separation to Lennard-Jones sigma defining
  ///        a van-der Waals clash.
  ///
  /// \param clash_vdw_ratio_in  The ratio to apply
  void setVdwClashRatio(double clash_vdw_ratio_in);
  
  /// \brief Set the diagnostic printing frequency.
  ///
  /// \param frequency_in  The chosen printing interval
  void setDiagnosticPrintFrequency(int frequency_in);

  /// \brief Set the checkpoint production flag.
  ///
  /// \param produce_in  Whether to produce a checkpoint at the end of the run
  void setCheckpointProduction(bool produce_in);
  
  /// \brief Set the electrostatic cutoff
  void setElectrostaticCutoff(double cutoff_in);
  
  /// \brief Set the Lennard-Jones cutoff
  void setLennardJonesCutoff(double cutoff_in);
  
  /// \brief Set the initial step length.
  ///
  /// \param step_size_in  The requested initial step length
  void setInitialStep(double step_size_in);

  /// \brief Set the convergence criterion, the target for the root mean squared value of all
  ///        gradients obtained after the minimization, in kcal/mol.
  ///
  /// \param target_in  The requested convergence target
  void setConvergenceTarget(double target_in);
  
private:
  ExceptionResponse policy;       ///< Set the behavior when bad inputs are encountered.  DIE =
                                  ///<   abort program, WARN = warn the user, and likely reset to
                                  ///<   the default value if one is available, SILENT = do not
                                  ///<   warn the user, but also likely reset to the default value
                                  ///<   if one is available.
  int total_cycles;               ///< Maximum number of minimization steps to attempt (equivalent
                                  ///<   to maxcyc in sander)
  int steepest_descent_cycles;    ///< Number of steepest descent steps to perform prior to
                                  ///<   beginning conjugate gradient moves (equivalent to ncyc in
                                  ///<   sander)
  int clash_damping_cycles;       ///< Number of clash-damping moves to include--this applies on
                                  ///<   top of steepest descent and conjugate gradient cycles. 
                                  ///<   The first clash_damping_cycles of the minimization will
                                  ///<   take place with special kernels that identify clashes
                                  ///<   between atoms and clamp the inter-particle distances at
                                  ///<   certain minimum values.  If there are no clashes, the
                                  ///<   provisions for clash-damping has no effect.
  int print_frequency;            ///< Print results at step 0 and, thereafter, after each interval
                                  ///<   of this many line minimizations.  The default of 0
                                  ///<   suppresses output except at the outset of the run.
  bool produce_checkpoint;        ///< Indicate that a checkpoint file should be produced at the
                                  ///<   end of the energy minimization run (default TRUE), with
                                  ///<   the name of the checkpoint file for each system found in
                                  ///<   the &files namelist.
  double electrostatic_cutoff;    ///< Cutoff for (short-ranged) electrostatic interactions, or for
                                  ///<   all gas-phase Coulombic electrostatics in a non-periodic
                                  ///<   system.  Units of Angstroms (A).
  double lennard_jones_cutoff;    ///< Cutoff for van-der Waals interactions, in Angstroms (A).
  double initial_step;            ///< Magnitude of the initial displacement along the gradient
                                  ///<   vector.  The size of subsequent moves will grow or shrink
                                  ///<   based on the history of success in previous optimizations.
                                  ///<   Units of Angstroms (A).
  double convergence_target;      ///< Convergence target for root mean squared value of all
                                  ///<   gradients obtained after the minimization, in kcal/mol-A.
  double clash_minimum_distance;  ///< The minimum separation between two particles, below which a
                                  ///<   clash will be declared
  double clash_vdw_ratio;         ///< The minimum ratio of inter-particle distance to the pairwise
                                  ///<   Lennard-Jones (van-der Waals) sigma parameter, below which
                                  ///<   a clash will be declared

  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;
  
  /// \brief Validate the total number of minimization cycles.
  void validateTotalCycles();

  /// \brief Validate the number of steepest descent or clash-damping cycles.
  ///
  /// \param cycles_in       The number of cycles for either auxiliary process.  May be modified
  ///                        and returned.
  /// \param default_cycles  The default number of cycles to impose if there is an error and
  ///                        automatic corrections are allowed (a behavior which depends on the
  ///                        runtime exception response setting)
  void validateAuxiliaryCycles(int *cycles_in, int default_cycles);
  
  /// \brief Validate the diagnostic printing frequency
  void validatePrintFrequency();

  /// \brief Validate the checkpointing behavior.
  ///
  /// \param directive  The command given for whether to write a checkpoint file
  void validateCheckpointProduction(const std::string &directive) const;
  
  /// \brief Validate the electrostatic cutoff.
  void validateElectrostaticCutoff();

  /// \brief Validate the Lennard-Jones cutoff.
  void validateLennardJonesCutoff();
  
  /// \brief Validate the initial step size.
  void validateInitialStep();

  /// \brief Validate the convergence target value.
  void validateConvergenceTarget();
};
  
/// \brief Produce a namelist for specifying energy minimization directives, similar to those found
///        in the &cntrl namelist of sander or pmemd.  This is a separate namelist from the
///        molecular dynamics input, obviating the need for the imin setting found in the general
///        &cntrl namelist of sander and pmemd.
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
NamelistEmulator minimizeInput(const TextFile &tf, int *start_line, bool *found,
                               ExceptionResponse policy = ExceptionResponse::DIE,
                               WrapTextSearch wrap = WrapTextSearch::NO);

} // namespace namelist
} // namespace stormm

#endif
