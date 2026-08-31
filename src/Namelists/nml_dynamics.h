// -*-c++-*-
#ifndef STORMM_NML_DYNAMICS_H
#define STORMM_NML_DYNAMICS_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/textfile.h"
#include "Potential/energy_enumerators.h"
#include "Structure/structure_enumerators.h"
#include "Trajectory/trajectory_enumerators.h"
#include "input.h"
#include "namelist_common.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using constants::CartesianDimension;
using constants::PrecisionModel;
using constants::UnitCellAxis;
using energy::VdwSumMethod;
using parse::WrapTextSearch;
using structure::ApplyConstraints;
using structure::RattleMethod;
using trajectory::BarostatKind;
using trajectory::ThermostatKind;
  
/// \brief Default number of molecular dynamics cycles
constexpr int default_dynamics_nstlim = 100;

/// \brief Default printing frequencies for diagnostics and trajectory snapshots
/// \{
constexpr int default_dynamics_ntpr = 0;
constexpr int default_dynamics_ntwx = 0;
/// \}

/// \brief The default number of steps at which to purge the system's center of mass motion.  This
///        includes rotation about the center of mass, if the system has no boundary conditions.
constexpr int default_dynamics_nscm = 10000;
  
/// \brief Default time step for molecular dynamics, in femtoseconds
constexpr double default_dynamics_time_step = 1.0;

/// \brief Default tolerance for RATTLE bond constraints
constexpr double default_rattle_tolerance = 1.0e-6;

/// \brief The default maximum number of RATTLE iterations that should be attempted
constexpr int default_rattle_max_iter = 20;

/// \brief The default behavior for geometry constraints is "OFF"
constexpr char default_geometry_constraint_behavior[] = "off";
  
/// \brief The default RATTLE approach (to be used on CPU resources only--GPU resources will use
///        the "CENTER_SUM" method)
constexpr char default_rattle_protocol[] = "sequential";
  
/// \brief The minimum molecular dynamics time step, in units of femtoseconds
constexpr double minimum_dynamics_time_step = 0.0009765625;

/// \brief The tightest possible RATTLE tolerance
constexpr double minimum_rattle_tolerance = 1.0e-9;

/// \brief The maximum number of iterations that RATTLE may run when attempting to reach the 
///        constrained geometry target
constexpr int maximum_rattle_iterations = 100;

/// \brief The default thermostat is no thermostat, anticipating an isoenergetic run
constexpr char default_thermostat_kind[] = "none";

/// \brief The default thermostat cahce configuration is 4-byte (single-precision, float32_t).  The
///        exact precision of random numbers is not so consequential, although the manner in which
///        they are generated could be.  This is not so much a precision model as a memory
///        optimization, and therefore included in this namelist rather than &precision.
constexpr char default_thermostat_cache_config[] = "single";
  
/// \brief The maximum depth that the random number cache can take, to guard against excessive
///        memory use.
constexpr int maximum_thermostat_cache_depth = 15;

/// \brief The default cache depth for a thermostat will give high efficiency in the number of
///        write transactions issued for the 256-bit read required: eight bytes written for every
///        byte read.
constexpr int default_thermostat_cache_depth = 1;

/// \brief The default Andersen thermostat velocity reassignment frequency, expressed as a number
///        of time steps.
constexpr int default_andersen_frequency = 10000;

/// \brief The default Langevin thermostat collision frequency, in events per femtosecond.
constexpr double default_langevin_frequency = 0.001;

/// \brief The default simulation temperature, chemistry's standard temperature and pressure (in
///        units of Kelvin).
constexpr double default_simulation_temperature = 298.15;

/// \brief The default random seed for thermostats.  This need not be the same random seed used
///        for any other application--in fact, having it be distinct guards against generators in
///        other parts of the program re-using the same sequence.
constexpr int default_thermostat_random_seed = 1329440765;

/// \brief The limits of the temperature evolution window default to no evolution--the initial
///        temperature is used to kick-start dynamics, if needed, and the equilibrium temperature
///        is immediately enforced.
/// \{
constexpr int default_tstat_evo_window_start = 0;
constexpr int default_tstat_evo_window_end   = 0;
/// \}

/// \brief The default barostat is no barostat, a constant volume condition
constexpr char default_barostat_kind[] = "none";

/// \brief The default external pressure is given in bar
constexpr double default_external_pressure = 1.0;

/// \brief The default MC barostat move rate is one per hundred ordinary steps.  An MC barostat
///        move is attempted by moving particles forward according to the typical integration
///        sequence (from WHITE to BLACK or BLACK to WHITE points in the coordinate time cycle),
///        then computing the energy (hereafter, E_i).  The unit cell then expands or contracts,
///        taking with it all coordinates therein (with consideration to constrained geometries)
///        by a random scaling factor that has been chosen for that specific MC barsotat attempt.
///        The energy is re-computed (with consideration to the volume differential against the
///        constant external pressure), hereafter E_f.  If a move succeeds according to the
///        Metropolis criterion after comparing E_i and E_f at the simulation target temperature,
///        the new coordinates and velocities (velocities may change due to constraints) will then
///        be retained.  If the move is rejected, the integration step will proceed again but
///        without the subsequent unit cell volume change, to recover the original phase space
///        exactly without using additional memory.
constexpr int default_mcbarostat_frequency = 100;

/// \brief The default MC barostat rescaling move will rescale the system by a factor between
///        0.9995 and 1.0005.
constexpr double default_mcbarostat_rescale = 0.0005;
  
/// \brief The warp multiplicity of non-bonded pairwise calculations can have an effect on the
///        overall speed of the computation.  While the unit cell subdivision is aggressive, it
///        is still a lot of work for one warp to do all tiles bewteen perhaps 50+ atoms in the
///        "tower" and 120+ atoms in the "plate" regions of a single neutral-territory region
///        controlled by one of the subdivisions.  Multiple warps can cooperate to get all of the
///        work done.  Allowing the user to control this parameter in certain situations may
///        permit further optimizations.
/// \{
constexpr int maximum_nt_warp_multiplicity = 8;
constexpr int default_nt_warp_multiplicity = 1;
/// \}

/// \brief Object to encapsulate molecular dynamics control information.  Like other namelist
///        encapsualtors, this object can take input file data as part of its construction, or
///        by a series of setters.  Validation of each piece of data is handled as it appears
///        either in the contructor or via setters.  Getter functions dispense the internal
///        information to any application using STORMM libraries.
class DynamicsControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \param found_nml   Indicator that the namelist was found in the input file
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &dynamics namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  DynamicsControls(ExceptionResponse policy_in = ExceptionResponse::DIE);
  
  DynamicsControls(const TextFile &tf, int *start_line, bool *found_nml,
                   ExceptionResponse policy_in = ExceptionResponse::DIE,
                   WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  DynamicsControls(const DynamicsControls &original) = default;
  DynamicsControls(DynamicsControls &&original) = default;
  DynamicsControls& operator=(const DynamicsControls &original) = default;
  DynamicsControls& operator=(DynamicsControls &&original) = default;
  /// \}
  
  /// \brief Get the total number of dynamics steps.
  int getStepCount() const;

  /// \brief Get the printing frequency for energies and other state variables of the system or
  ///        collection of systems.
  int getDiagnosticPrintFrequency() const;
  
  /// \brief Get the printing frequency for trajectory snapshots.
  int getTrajectoryPrintFrequency() const;

  /// \brief Get the frequency for purging motion of the system's center of mass.
  int getCenterOfMassMotionPurgeFrequency() const;
  
  /// \brief Get the simulation time step.
  double getTimeStep() const;

  /// \brief Get the electrostatic pairwise cutoff in periodic simulations.
  double getElectrostaticCutoff() const;
    
  /// \brief Get the van-der Waals pairwise cutoff in periodic simulations.
  double getVanDerWaalsCutoff() const;
    
  /// \brief Get the van-der Waals pairwise cutoff in periodic simulations.
  double getLennardJonesCutoff() const;

  /// \brief Get the method for computing van-der Waals interactions between particles.
  VdwSumMethod getVdwSummation() const;

  /// \brief Get the value of Coulomb's constant.
  double getCoulombConstant() const;

  /// \brief Get the electrostatic 1:4 screening value.
  double getElec14Screening() const;
  
  /// \brief Get the van-der Waals 1:4 screening value.
  double getVdw14Screening() const;

  /// \brief Report whether the value of Coulomb's constant was set explicitly or left at the
  ///        default value.
  bool coulombSetByUser() const;

  /// \brief Report whether the electrostatic 1:4 scaling was set explicitly or left at its
  ///        default value.
  bool elec14SetByUser() const;
  
  /// \brief Report whether the van-der Waals 1:4 scaling was set explicitly or left at its
  ///        default value.
  bool vdw14SetByUser() const;
  
  /// \brief Indicate whether general geometric constraints (RATTLE for hub-and-spoke bonded
  ///        groups, SETTLE for rigid water molecules) are to be implemented.
  ApplyConstraints constrainGeometry() const;

  /// \brief Indicate whether SHAKE (and, by extension, RATTLE in the velocity-Verlet integrator)
  ///        should be implemented.
  ApplyConstraints useShake() const;

  /// \brief Indicate whether SHAKE (and, by extension, RATTLE in the velocity-Verlet integrator)
  ///        should be implemented.
  ApplyConstraints useSettle() const;

  /// \brief Get the RATTLE tolerance.
  double getRattleTolerance() const;

  /// \brief Get the maximum number of RATTLE iterations.
  int getRattleIterations() const;

  /// \brief Get the CPU-based method of implementing RATTLE for multiple constrained bonds
  ///        connected to one central atom.  This determines the order in which the position of the
  ///        central atom is updated with each iteration.
  RattleMethod getCpuRattleMethod() const;

  /// \brief Get the thermostat kind.
  ThermostatKind getThermostatKind() const;

  /// \brief Get the starting step for temperature evolution from the initial temperature profile
  ///        to the final temperature profile across all atoms and all systems.
  int getThermostatEvolutionStart() const;

  /// \brief Get the final step for temperature evolution from the initial temperature profile
  ///        to the final temperature profile across all atoms and all systems.
  int getThermostatEvolutionEnd() const;
  
  /// \brief Get the depth of the thermostat's random number cache.
  int getThermostatCacheDepth() const;

  /// \brief Get the random seed to be used in the thermostat's random number generators.
  int getThermostatSeed() const;

  /// \brief Get the number of distinct layers to the thermostat, the number of different
  ///        temperature and system groups that must be parsed (in order) that the temperature
  ///        profiles at which to maintain all systems of the synthesis are to be maintained.
  ///        This function also serves as a check on the consistency of array sizes: for each
  ///        initial temeprature specification, there must be a final temperature specification,
  ///        a label, a label index, and an atom mask (even if the latter three are "all", -1
  ///        (integer code for all), and "@=" (all atoms).
  int getThermostatLayerCount() const;

  /// \brief Get the frequncy of velocity reassignments for a "massive" Andersen thermostat, one
  ///        in which all velocities are reassigned after a particular number of steps.
  int getAndersenFrequency() const;

  /// \brief Get the Langevin collison frequency, in units of inverse femtoseconds.
  double getLangevinFrequency() const;

  /// \brief Get the type of barostat to be put into effect.
  BarostatKind getBarostatKind() const;

  /// \brief Get the frequency of Monte-Carlo barostat moves, in time steps between attempts.
  ///
  /// Overloaded:
  ///   - Specify no dimension to produce the isotropic move frequency (this will raise an
  ///     exception if moves are anisotropic)
  ///   - Specify the Cartesian dimension most like the unit cell axis of interest
  ///   - Specify the unit cell axis of interest
  /// \{
  int getMCBarostatFrequency() const;
  int getMCBarostatFrequency(const CartesianDimension dim) const;
  int getMCBarostatFrequency(const UnitCellAxis dim) const;
  /// \}
  
  /// \brief Get the rescaling factor for Monte-Carlo barostat moves.  Overloading and descriptions
  ///        of input parameters follow from getMCBarostatFrequency(), above.
  /// \{
  double getMCBarostatRescaling() const;
  double getMCBarostatRescaling(const CartesianDimension dim) const;
  double getMCBarostatRescaling(const UnitCellAxis dim) const;
  /// \}
  
  /// \brief Get the thermostat cache configuration.
  PrecisionModel getThermostatCacheConfig() const;
  
  /// \brief Get the vector of initial temperature targets for all groups of atoms and systems.
  const std::vector<double>& getInitialTemperatureTargets() const;

  /// \brief Get the vector of final temperature targets for all groups of atoms and systems.
  const std::vector<double>& getFinalTemperatureTargets() const;

  /// \brief Get the list of label groups within the systems cache to which each thermostated
  ///        group of atoms belongs.
  const std::vector<std::string>& getThermostatLabels() const;

  /// \brief Get the indices of systems within the label groups to which each thermostated group of
  ///        atoms belongs.
  const std::vector<int>& getThermostatLabelIndices() const;

  /// \brief Get the atom mask strings for each thermostated group of atoms.
  const std::vector<std::string>& getThermostatMasks() const;

  /// \brief Get the full array of external pressures specified in the series of pressure keyword
  ///        declarations.
  const std::vector<double>& getExternalPressures() const;

  /// \brief Get the full array of system cache labels specified in the series of pressure keyword
  ///        declarations.
  const std::vector<std::string>& getBarostatLabels() const;

  /// \brief Get the full array of system cache group indices in the series of pressure keyword
  ///        declarations.
  const std::vector<int>& getBarostatLabelIndices() const;
  
  /// \brief Get the warp multiplicity for non-bonded pairwise calculations.
  int getNTWarpMultiplicity() const;
  
  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator& getTranscript() const;
  
  /// \brief Set the total number of dynamics steps.
  ///
  /// \param total_step_count_in  The number of steps to take
  void setStepCount(int total_step_count_in);

  /// \brief Set the diagnostic printing frequency.
  ///
  /// \param diagnostic_frequency_in  The frequency with which to store energy and state variables
  void setDiagnosticPrintFrequency(int diagnostic_frequency_in);
  
  /// \brief Set the trajectory snapshot printing frequency.
  ///
  /// \param trajectory_frequency_in  The frequency with which to store energy and state variables
  void setTrajectoryPrintFrequency(int trajectory_frequency_in);

  /// \brief Set the center of mass motion purge frequency.
  ///
  /// \param com_motion_purge_ferquency_in  The frequency with which to purge center of mass motion
  void setCenterOfMassMotionPurgeFrequency(int com_motion_purge_ferquency_in);
  
  /// \brief Set the simulation time step
  ///
  /// \param time_step_in  The requested time step, in femtoseconds
  void setTimeStep(double time_step_in);

  /// \brief Set the short-ranged electrostatic cutoff for the simulation.
  ///
  /// \param cutoff_in  The requested particle-particle electrostatic cutoff, in Angstroms
  void setElectrostaticCutoff(double cutoff_in);

  /// \brief Set the short-ranged van-der Waals cutoff for the simulation.
  ///
  /// \param cutoff_in  The requested particle-particle van-der Waals cutoff, in Angstroms
  void setVanDerWaalsCutoff(double cutoff_in);

  /// \brief Set the cutoffs for short-ranged van-der Waals as well as electrostatic interactions
  ///        in a periodic simulation.
  ///
  /// \param cutoff_in  The requested cutoff, in Angstroms
  void setCutoff(double cutoff_in);

  /// \brief Set the strategy for evaluating the tails of van-der Waals (Lennard-Jones)
  ///        interactions.
  ///
  /// Overloaded:
  ///   - Provide a keyword that indicates some value of the enumerator
  ///   - Provide a value of the enumerator itself
  ///
  /// \param vdw_method_in  The input determining how van-der Waals interactions will vanish
  /// \{
  void setVdwSummation(const std::string &vdw_method_in);
  void setVdwSummation(const VdwSumMethod vdw_method_in);
  /// \}
  
  /// \brief Set the value of Coulomb's constant.
  ///
  /// \param coulomb_in  Define Coulomb's constant for the simulation
  void setCoulombConstant(double coulomb_in);

  /// \brief Set the screening factor for 1:4 electrostatic interactions.
  ///
  /// \param screening_in  The screening factor to set
  void setElec14Screening(double screening_in);
  
  /// \brief Set the screening factor for 1:4 van-der Waals interactions.
  ///
  /// \param screening_in  The screening factor to set
  void setVdw14Screening(double screening_in);
  
  /// \brief Stipulate whether geometric constraints will be implemented.
  ///
  /// Overloaded:
  ///   - Provide a parseable string to indicate the constraint instruction
  ///   - Provide the explicit enumeration
  ///
  /// \param constrain_geometry_in  Indicate whether to activate RATTLE and SETTLE
  /// \{
  void setGeometricConstraints(const std::string &constrain_geometry_in);
  void setGeometricConstraints(ApplyConstraints constrain_geometry_in = ApplyConstraints::YES);
  /// \}

  /// \brief Stipulate whether bonds to hydrogen atoms will be held rigid at their equilibrium
  ///        length.  Overloading follows from setGeometricConstraints(), above.
  ///
  /// \param use_shake_in
  /// \{
  void setUseShake(const std::string &use_shake_in);
  void setUseShake(ApplyConstraints use_shake_in = ApplyConstraints::YES);
  /// \}
  
  /// \brief Stipulate whether three-atom isosceles arrangements with two light atoms (e.g. water)
  ///        will be held rigid at their equilibrium geometries.  Overloading follows from
  ///        setGeometricConstraints(), above.
  ///
  /// \param use_settle_in
  /// \{
  void setUseSettle(const std::string &use_settle_in);
  void setUseSettle(ApplyConstraints use_settle_in = ApplyConstraints::YES);
  /// \}
  
  /// \brief Set the RATTLE protocol to be used in CPU operations (the GPU will always use the
  ///        "CENTER_SUM" approach).
  ///
  /// Overloaded:
  ///   - Provide a string to be translated into the protocol enumeration
  ///   - Provide the enumeration itself
  ///
  /// \param rattle_protocol_in  The RATTLE protocol to set
  /// \{
  void setCpuRattleMethod(const std::string &rattle_protocol_in);
  
  /// \brief Set the simulation RATTLE tolerance.
  ///
  /// \param tol_in  The requested tolerance
  void setRattleTolerance(double rattle_tolerance_in);

  /// \brief Get the thermostat kind.
  ///
  /// Overloaded:
  ///   - Provide the type by parseable string
  ///   - Provide the explicit enumeration
  ///
  /// \param thermostat_kind_in  The type of thermostat requested
  /// \{
  void setThermostatKind(const std::string &thermostat_kind_in);
  void setThermostatKind(ThermostatKind thermostat_kind_in);
  /// \}

  /// \brief Set the starting step for temperature evolution from the initial temperature profile
  ///        to the final temperature profile across all atoms and all systems.
  ///
  /// \param step_number  The requested step number
  void setThermostatEvolutionStart(int thermo_evolution_start_in);

  /// \brief Set the final step for temperature evolution from the initial temperature profile
  ///        to the final temperature profile across all atoms and all systems.
  ///
  /// \param step_number  The requested step number
  void setThermostatEvolutionEnd(int step_number);

  /// \brief Set the depth of the thermostat's random number cache.
  ///
  /// \param depth_in  The depth that the cache shall take
  void setThermostatCacheDepth(int depth_in);

  /// \brief Set the random seed to be used in the thermostat's random number generators.
  ///
  /// \param igseed  The seed to set
  void setThermostatSeed(int igseed);

  /// \brief Set the velocity restart frequency for a "massive" Andersen thermostat.
  ///
  /// \param frequency_in  The step frequency to set
  void setAndersenFrequency(int frequency_in);

  /// \brief Set the stochastic collision frequency for a Langevin thermostat.
  ///
  /// \param frequency_in  The collision frequency to set, in units of inverse femtoseconds
  void setLangevinFrequency(double frequency_in);

  /// \brief Set the type of barostat that will be applied to all systems.
  ///
  /// Overloaded:
  ///   - Provide the type by parseable string
  ///   - Provide the explicit enumeration
  ///
  /// \param barostat_kind_in  The type of barostat requested
  /// \{
  void setBarostatKind(const std::string &barostat_kind_in);
  void setBarostatKind(BarostatKind barostat_kind_in);
  /// \}

  /// \brief Set the number of steps between unit cell rescaling moves made by the Monte-Carlo
  ///        barostat.
  ///
  /// Overloaded:
  ///   - Set the rate along all unit cell axes (isotropic rescaling)
  ///   - Set the rate along a specific unit cell axis (anisotropic rescaling)
  ///
  /// \param frequency_in  The requested rate of MC barostat moves
  /// \param dim           Unit cell axis along which the specified move frequency will apply
  /// \{
  void setMCBarostatFrequency(int frequency_in);
  void setMCBarostatFrequency(int frequency_in, CartesianDimension dim);
  void setMCBarostatFrequency(int frequency_in, UnitCellAxis dim);
  /// \}

  /// \brief Set the range of possible Monte-Carlo unit cell volume rescaling moves.  Overloading
  ///        follows from setMCBarostatFrequency(), above.
  ///
  /// \param mcbarostat_rescale_in  The isotropic rescaling factor to set
  /// \param dim                    The unit cell axis along which the specified rescaling factor
  ///                               will be applied
  /// \{
  void setMCBarostatRescaling(double mcbarostat_rescale_in);
  void setMCBarostatRescaling(double mcbarostat_rescale_in, CartesianDimension dim);
  void setMCBarostatRescaling(double mcbarostat_rescale_in, UnitCellAxis dim);
  /// \}
  
  /// \brief Set the thermostat cache configuration.
  ///
  /// Overloaded:
  ///   - Provide a parseable string to indicate the configuration
  ///   - Provide the explicit enumeration
  ///
  /// \param cache_config_in  The cache configuration to set
  /// \{
  void setThermostatCacheConfig(const std::string &cache_config_in);
  void setThermostatCacheConfig(PrecisionModel cache_config_in);
  /// \}

  /// \brief Set the temperature target and start-up schedule for a thermostat group affecting one
  ///        or more systems within the synthesis.
  ///
  /// \param initial_target  The initial target temperature
  /// \param final_target    The equilibrium target temperature
  /// \param label           The label group from within the cache to which the stated temperatures
  ///                        will be applied
  /// \param label_index     Index of the system within the label group to which the stated
  ///                        temperatures apply, e.g. 4 would indicate the fifth system in the
  ///                        named label group.  A value of -1 implies that the temperature
  ///                        regulation will apply to all systems in the label group.
  /// \param mask            String encoding an atom mask for the atoms in the label group and
  ///                        index to which the stated temperatures apply
  void setThermostatGroup(double initial_target = default_simulation_temperature,
                          double final_target = default_simulation_temperature,
                          const std::string &label = std::string("all"), int label_index = -1,
                          const std::string &mask = std::string("@="));

  /// \brief Set the pressure target for a barostat group affecting one or more systems within the
  ///        synthesis.
  ///
  /// \param pressure_target  The target external pressure of the barostat
  /// \param label            Label group from within the systems cache to which the stated
  ///                         pressure will be applied
  /// \param label_index      Index of the system within the label group to which the stated
  ///                         pressure shall apply to all systems in the label group.  The default
  ///                         value of -1 implies all systems in the label group.
  void setBarostatGroup(double pressure_target = default_external_pressure,
                        const std::string &label = std::string("all"), int label_index = -1);

  /// \brief Get the warp multiplicity for non-bonded pairwise calculations.  This applies only to
  ///        periodic dynamics.
  ///
  /// \param mult_in  The warp multiplicity to set
  void setNTWarpMultiplicity(int mult_in);

private:
  ExceptionResponse policy;        ///< Set the behavior when bad inputs are encountered.  DIE =
                                   ///<   abort program, WARN = warn the user, and likely reset to
                                   ///<   the default value if one is available, SILENT = do not
                                   ///<   warn the user, but also likely reset to the default value
                                   ///<   if one is available.
  int total_step_count;            ///< Total number of dynamics steps to perform (equivalent to
                                   ///<   maxcyc in sander)
  int diagnostic_frequency;        ///< The frequency (step interval) with which to print energy
                                   ///<   and state variable diagnostics
  int trajectory_frequency;        ///< The frequency (step interval) with which to store the
                                   ///<   coordinates of all systems
  int com_motion_purge_frequency;  ///< The frequency (step interval) with which to purge the
                                   ///<   motion of the center of mass (including, if appropriate,
                                   ///<   rotation about the center of mass)
  double time_step;                ///< Time step to take after each force evaluation
  double electrostatic_cutoff;     ///< Cutoff applied to electrostatic short-ranged interactions
                                   ///<   in periodic simulations
  double van_der_waals_cutoff;     ///< Cutoff applied to van-der Waals short-ranged interactions
                                   ///<   in periodic simulations
  double coulomb;                  ///< Coulomb's constant may be defined by the user in the input
  VdwSumMethod vdw_style;          ///< The manner in which van-der Waals interactions will vanish
                                   ///<   at the cutoff
  double elec_14_screening;        ///< The screening factor on 1:4 non-bonded interactions.  To
                                   ///<   specify that such interactions take place at 5/6 nominal
                                   ///<   strength, provide the value 6/5, 1.2.
  double vdw_14_screening;         ///< The screening factor on 1:4 non-bonded van-der Waals
                                   ///<   (Lennard-Jones) interactions
  bool coulomb_set_by_user;        ///< Flag to indicate that Coulomb's constant has been set by
                                   ///<   the user.  This will also read TRUE if a value is fed in
                                   ///<   through the setCoulombConstant() member function.
  bool elec_14_set_by_user;        ///< Flag to indicate that the electrostatic 1:4 screening
                                   ///<   factor was set by the user.  This will also read true if
                                   ///<   a value is fed in through the setElec14Screeningt()
                                   ///<   member function.
  bool vdw_14_set_by_user;         ///< Flag to indicate that the van-der Waals 1:4 screening
      	      	      	      	   ///<   factor was set by the user.  This will also read true if
      	      	      	      	   ///<   a value is fed in through the setVdw14Screening() member
                                   ///<   function.
  std::string constrain_geometry;  ///< Indicate whether SHAKE bond length constraints and SETTLE
                                   ///<   rigid water constraints should be implemented.  This will
                                   ///<   enforce both types of constraints.
  std::string use_shake;           ///< Indicate whether SHAKE (and, by extension, RATTLE in the
                                   ///<   velocity-Verlet integrator) should be implemented.
  std::string use_settle;          ///< Indicate whether SETTLE rigid water constraints should be
                                   ///<   implemented.
  double rattle_tolerance;         ///< The tolerance to apply to bond constraint calculations
  int rattle_iterations;           ///< Maximum number of RATTLE iterations to attempt in order to
                                   ///<   achieve convergence at or below the specified tolerance
  std::string rattle_protocol;     ///< The protocol for implementing RATTLE.  This affects the way
                                   ///<   in which the position of the central atom is updated in
                                   ///<   hub-and-spoke constraint groups.
  std::string thermostat_kind;     ///< String encoding the type of thermostat to be used for all
                                   ///<   systems
  int thermo_evolution_start;      ///< The step at which to begin evolution of the thermostat
                                   ///<   target temperature from the initial to the final state
  int thermo_evolution_end;        ///< The step at which to conclude evolution of the thermostat
                                   ///<   target temperature from the initial to the final state
  int thermostat_cache_depth;      ///< The quantity of random numbers that will be held in reserve
                                   ///<   by the thermostat.  This is a matter of optimization and
                                   ///<   the actual quantity of random results held will be three
                                   ///<   times the user-specified value (x, y, and z components
                                   ///<   for every atom).
  int thermostat_seed;             ///< Random number generator seed for the thermostat
  int andersen_frequency;          ///< Frequency (step count) at which all particle velocities are
                                   ///<   reset in the Andersen thermostating scheme
  double langevin_frequency;       ///< The frequency of stochastic Langevin collisions, in units
                                   ///<   of inverse femtoseconds
  std::string barostat_kind;       ///< String encoding the type of barostat to be used to for all
                                   ///<   systems
  int mcbarostat_frequency_a;      ///< Number of moves between attempted MC barostat volumetric
                                   ///<   changes along the unit cell "a" axis
  int mcbarostat_frequency_b;      ///< Number of moves between attempted MC barostat volumetric
                                   ///<   changes along the unit cell "a" axis
  int mcbarostat_frequency_c;      ///< Number of moves between attempted MC barostat volumetric
                                   ///<   changes along the unit cell "a" axis
  double mcbarostat_rescale_a;     ///< Monte-Carlo barostat moves will be chosen from a uniform
                                   ///<   distribution in the range +/- mcbarostat_rescale_a along
                                   ///<   the unit cell "a" axis
  double mcbarostat_rescale_b;     ///< Monte-Carlo barostat moves will be chosen from a uniform
                                   ///<   distribution in the range +/- mcbarostat_rescale_b along
                                   ///<   the unit cell "b" axis
  double mcbarostat_rescale_c;     ///< Monte-Carlo barostat moves will be chosen from a uniform
                                   ///<   distribution in the range +/- mcbarostat_rescale_c along
                                   ///<   the unit cell "c" axis
  int nt_warp_multiplicity;        ///< The number of warps that will cooperate to complete work in
                                   ///<   each neutral-territory subdivision of periodic dynamics
                                   ///<   unit cells
  
  /// Configuration for the thermostat's random number cache (if applicable)
  std::string thermostat_cache_config;    
  
  /// A series of initial target temperatures for various thermostat partitions of the synthesis
  /// of systems
  std::vector<double> initial_temperature_targets;

  /// A series of final target temperatures for various thermostat partitions of the synthesis of
  /// systems
  std::vector<double> final_temperature_targets;

  /// A series of label strings describing which groups from the systems cache each temperature
  /// profile will apply to.
  std::vector<std::string> thermostat_labels;

  /// A series of indices within each label group to which each temperature profile applies.
  /// Each declaration of the "temperature" keyword may specify a specific index within the chosen
  /// label group (valid index >= 0), or all indices within that group (any negative index, the
  /// deault case).
  std::vector<int> thermostat_label_indices;

  /// A series of strings expressing atom masks for subgroups of atoms within each identified
  /// system that will be subject to each temperature profile
  std::vector<std::string> thermostat_masks;

  /// A series of external pressures controlling barostats acting on each system in the entire
  /// synthesis.  Different simulations may run at different pressures, although the mechanics of
  /// all barostats (e.g. the rate of Monte-Carlo attempts to change the volume, and the size of
  /// each change) must be the same.
  std::vector<double> external_pressures;

  /// A series of label strings describing which groups from the systems cache each pressure
  /// setting will apply to.  Analogous to thermostat_labels.
  std::vector<std::string> barostat_labels;

  /// A series of indices within each label group to which each pressure setting applies.  
  std::vector<int> barostat_label_indices;
  
  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;
  
  /// \brief Validate the total number of steps.
  void validateStepCount();

  /// \brief Validate the diagonstic print frequency.
  void validateDiagnosticPrintFrequency();
  
  /// \brief Validate the trajectory print frequency.
  void validateTrajectoryPrintFrequency();

  /// \brief Validate the requested rate at which to purge the motion of the center of mass.  If
  ///        negative values are applied, a frequency of zero (no removal, apart from the initial
  ///        purge with a velocity kick-start) will be applied.
  void validateCenterOfMassMotionPurgeFrequency();
  
  /// \brief Validate the time step.
  void validateTimeStep();

  /// \brief Validate the RATTLE tolerance.
  void validateRattleTolerance();

  /// \brief Validate the number of RATTLE iterations.
  void validateRattleIterations();

  /// \brief Validate a temperature requested for thermoregulation.
  ///
  /// \param t  The temperature in question
  void validateTemperature(double t) const;
  
  /// \brief Validate the type of thermostat to be used across all systems.
  void validateThermostatKind();

  /// \brief Validate the amount of random number caching requested.
  void validateCacheDepth();

  /// \brief Validate the type of barostat to be used across all systems.
  void validateBarostatKind();

  /// \brief Validate the requested pressure.
  ///
  /// \param p  The value of external pressure in question
  void validatePressure(const double p) const;

  /// \brief Validate the Monte-Carlo barostat move frequency.
  void validateMCBarostatFrequency();

  /// \brief Validate the Monte-Carlo barostat rescaling factor.
  void validateMCBarostatRescaling();
  
  /// \brief Validate the precision model for caching random numbers.
  void validateCacheConfiguration();

  /// \brief Validate the warp multiplicity for non-bonded pairwise calculations.
  void validateNTWarpMultiplicity();
};
  
/// \brief Produce a namelist for specifying molecular dynamics directives, similar to those found
///        in the &cntrl namelist of sander or pmemd.  This is a separate namelist from the energy
///        minimization input.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will not wrap back to the beginning of the TextFile object, as the
///                    &rst namelist is intended to be repeatable)
/// \param found       Indicator that the namelist was present in the input file
/// \param policy      Reaction to exceptions encountered during namelist reading
/// \param wrap        Indicate that the search for an &dynamics namelist should carry on from the
///                    beginning of an input file if no such namelist is found starting from the
///                    original starting point
NamelistEmulator dynamicsInput(const TextFile &tf, int *start_line, bool *found,
                               ExceptionResponse policy = ExceptionResponse::DIE,
                               WrapTextSearch wrap = WrapTextSearch::NO);

} // namespace namelist
} // namespace stormm

#endif
