// -*-c++-*-
#ifndef STORMM_THERMOSTAT_H
#define STORMM_THERMOSTAT_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "Math/vector_ops.h"
#include "Namelists/nml_dynamics.h"
#include "Reporting/error_format.h"
#include "Structure/structure_enumerators.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/synthesis_cache_map.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

using card::GpuDetails;
using card::Hybrid;
using card::HybridTargetLevel;
using constants::PrecisionModel;
using data_types::isSignedIntegralScalarType;
using namelist::default_andersen_frequency;
using namelist::default_langevin_frequency;
using namelist::default_rattle_tolerance;
using namelist::default_simulation_temperature;
using namelist::default_thermostat_cache_depth;
using namelist::default_thermostat_random_seed;
using namelist::default_tstat_evo_window_start;
using namelist::default_tstat_evo_window_end;
using namelist::DynamicsControls;
using namelist::maximum_thermostat_cache_depth;
using stmath::findBin;
using structure::ApplyConstraints;
using symbols::boltzmann_constant_gafs;
using symbols::boltzmann_constant_gafs_f;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using synthesis::SyAtomUpdateKit;
using synthesis::SynthesisCacheMap;
using synthesis::SystemCache;
using topology::ChemicalDetailsKit;
  
/// \brief Partially writeable abstract for the Thermostat object.  As with the MMControlKit struct
///        (see the library header mm_controls.h), the step member variable is modifiable so that
///        it can be updated without setting it in the parent object and regenerating the entire
///        abstract.
template <typename T> struct ThermostatWriter {

  /// \brief The constructor takes the usual list of parent object attributes with pointers
  ///        applicable on either the CPU host or GPU device.
  ThermostatWriter(ThermostatKind kind_in, ThermostatPartition layout_in, int natom_in,
                   size_t padded_natom_in, int npart_in, int step_in, int depth_in,
                   int init_evolution_in, int end_evolution_in, T init_temperature_in,
                   T final_temperature_in, T dt_in, bool csnt_geom_in, T rattle_tol_in,
                   int rattle_iter_in, int andr_cyc_in, T gamma_ln_in, T ln_implicit_in,
                   T ln_explicit_in, const int4* partition_bounds_in,
                   const T* init_temperatures_in, const T* final_temperatures_in,
                   ullint2* state_xy_in, ullint2* state_zw_in, PrecisionModel rng_mode_in,
                   double* cache_in, float* sp_cache_in);

  /// \brief The copy and move assignment operators will be implicitly deleted due to the presence
  ///        of const member variables.
  /// \{
  ThermostatWriter(const ThermostatWriter &original) = default;
  ThermostatWriter(ThermostatWriter &&original) = default;
  /// \}
  
  const ThermostatKind kind;         ///< The type of this thermostat, i.e. Andersen
  const ThermostatPartition layout;  ///< The manner in which different atoms or systems are held
                                     ///<   to one or more heat baths
  const int natom;                   ///< Number of atoms controlled by this thermostat
  const size_t padded_natom;         ///< The number of atoms padded by the warp size
  const int npart;                   ///< Number of unique partitions in the thermostat, each
                                     ///<   being the atom index bounds on a distinct heat bath
                                     ///<   at a particular temperature
  int step;                          ///< Current step number of the simulation
  const int depth;                   ///< Depth of the cache (the true depth will be three times
                                     ///<   higher, to serve each atom with Cartesian X, Y, and Z
                                     ///<   influences)
  const int init_evolution;          ///< Final step at which the thermostat target temperatures
                                     ///<   come from init_temperature or the eponymous array, and
                                     ///<   the first step at which temperature targets begin to
                                     ///<   shift linearly towards final_temperature(s)
  const int end_evolution;           ///< The first step at which the thermostat begins to target
                                     ///<   final_temperature(s)
  const T init_temperature;          ///< Initial target temperature, common to all atoms (if it is
                                     ///<   applicable)
  const T final_temperature;         ///< Final target temperature, common to all atoms
  const T dt;                        ///< The time step, in units of femtoseconds
  const bool cnst_geom;              ///< Indicate whether to apply geometry constraints
  const T rattle_tol;                ///< Convergence criterion for RATTLE bond length constraints
  const int rattle_iter;             ///< Number of iterations to allow when attempting to converge
                                     ///<   RATTLE
  const int andr_cyc;                ///< The frequency of Andersen velocity resets, in terms of
                                     ///<   simulation steps.  Multiply by the time step dt to get
                                     ///<   the inverse strength of the Andersen thermostat.
  const T gamma_ln;                  ///< The "Gamma" factor for the Langevin thermostat, giving
                                     ///<   the number of collisions per femtosecond scaled by the
                                     ///<   conversion factor for taking forces in kcal/mol-A into
                                     ///<   Angstroms per femtosecond.
  const T ln_implicit;               ///< Implicit Langevin factor, applied in the first velocity
                                     ///<   Verlet update after new forces have been computed
  const T ln_explicit;               ///< Explicit Langevin factor, applied in the second velocity
                                     ///<   Verlet update as new coordinates are determined
  const int4* partition_bounds;      ///< Bounds array on the atoms that this thermostat regulates.
                                     ///<   Also includes the number of degrees of freedom in each
                                     ///<   partition.
  const T* init_temperatures;        ///< Array of initial target temperatures for each atom
  const T* final_temperatures;       ///< Array of final target gemperatures for each atom
  ullint2* state_xy;                 ///< First 128 bits of each 256-bit state vector
  ullint2* state_zw;                 ///< Second 128 bits of each 256-bit state vector
  const PrecisionModel rng_mode;     ///< The mode in which random numbers are available for this
                                     ///<   thermostat.  One of cache or sp_cache, below, will be
                                     ///<   valid as indicated by this member variable.
  double* cache;                     ///< Cache of double-precision pre-computed random values
                                     ///<   (depth gives the number of such values for each atom)
  float* sp_cache;                   ///< Cache of single-precision pre-computed random values
                                     ///<   (depth gives the number of such values for each atom)
};

/// \brief Read-only abstract for the Thermostat object.
template <typename T> struct ThermostatReader {

  /// \brief The constructor takes the usual list of parent object attributes with pointers
  ///        applicable on either the CPU host or GPU device.  As with some other readers, there
  ///        is a means to construct this one by "const-ifying" the writer.
  /// \{
  ThermostatReader(ThermostatKind kind_in, ThermostatPartition layout_in, int natom_in,
                   size_t padded_natom_in, int npart_in, int step_in, int depth_in,
                   int init_evolution_in, int end_evolution_in, T init_temperature_in,
                   T final_temperature_in, T dt_in, bool csnt_geom_in, T rattle_tol_in,
                   int rattle_iter_in, int andr_cyc_in, T gamma_ln_in, T ln_implicit_in,
                   T ln_explicit_in, const int4* partition_bounds_in,
                   const T* init_temperatures_in, const T* final_temperatures_in,
                   const ullint2* state_xy_in, const ullint2* state_zw_in,
                   PrecisionModel rng_mode_in, const double* cache_in, const float* sp_cache_in);

  ThermostatReader(const ThermostatWriter<T> &w);
  /// \}
  
  /// \brief The copy and move assignment operators will be implicitly deleted due to the presence
  ///        of const member variables.
  /// \{
  ThermostatReader(const ThermostatReader &original) = default;
  ThermostatReader(ThermostatReader &&original) = default;
  /// \}
  
  const ThermostatKind kind;         ///< The type of this thermostat, i.e. Andersen
  const ThermostatPartition layout;  ///< The manner in which different atoms or systems are held
                                     ///<   to one or more heat baths
  const int natom;                   ///< Number of atoms controlled by this thermostat
  const size_t padded_natom;         ///< The number of atoms padded by the warp size
  const int npart;                   ///< Number of unique partitions in the thermostat, each
                                     ///<   being the atom index bounds on a distinct heat bath
                                     ///<   at a particular temperature
  const int step;                    ///< Current step number of the simulation
  const int depth;                   ///< Depth of the cache (the true depth will be three times
                                     ///<   higher, to serve each atom with Cartesian X, Y, and Z
                                     ///<   influences)
  const int init_evolution;          ///< Final step at which the thermostat target temperatures
                                     ///<   come from init_temperature or the eponymous array, and
                                     ///<   the first step at which temperature targets begin to
                                     ///<   shift linearly towards final_temperature(s)
  const int end_evolution;           ///< The first step at which the thermostat begins to target
                                     ///<   final_temperature(s)
  const T init_temperature;          ///< Initial target temperature, common to all atoms (if it is
                                     ///<   applicable)
  const T final_temperature;         ///< Final target temperature, common to all atoms (if it is
                                     ///<   applicable)
  const T dt;                        ///< The time step, in units of femtoseconds
  const bool cnst_geom;              ///< Indicate whether to apply geometry constraints
  const T rattle_tol;                ///< Convergence criterion for RATTLE bond length constraints
  const int rattle_iter;             ///< Number of iterations to allow when attempting to converge
                                     ///<   RATTLE
  const int andr_cyc;                ///< The frequency of Andersen velocity resets, in terms of
                                     ///<   simulation steps.  Multiply by the time step dt to get
                                     ///<   the inverse strength of the Andersen thermostat.
  const T gamma_ln;                  ///< The "Gamma" factor for the Langevin thermostat, giving
                                     ///<   the number of collisions per femtosecond
  const T ln_implicit;               ///< Implicit Langevin factor, applied in the first velocity
                                     ///<   Verlet update after new forces have been computed
  const T ln_explicit;               ///< Explicit Langevin factor, applied in the second velocity
                                     ///<   Verlet update as new coordinates are determined
  const int4* partition_bounds;      ///< Bounds array on the atoms that this thermostat regulates.
                                     ///<   Also includes the number of degrees of freedom in each
                                     ///<   partition.
  const T* init_temperatures;        ///< Array of initial target temperatures for each atom
  const T* final_temperatures;       ///< Array of final target gemperatures for each atom
  const ullint2* state_xy;           ///< First 128 bits of each 256-bit state vector
  const ullint2* state_zw;           ///< Second 128 bits of each 256-bit state vector
  const PrecisionModel rng_mode;     ///< The mode in which random numbers are available for this
                                     ///<   thermostat.  One of cache or sp_cache, below, will be
                                     ///<   valid as indicated by this member variable.
  const double* cache;               ///< Cache of double-precision pre-computed random values
                                     ///<   (depth gives the number of such values for each atom)
  const float* sp_cache;             ///< Cache of single-precision pre-computed random values
                                     ///<   (depth gives the number of such values for each atom)
};

/// \brief Store the parameters for a simulation thermostat.  Includes Berendsen, Andersen, and
///        Langevin methods.  This class can be assembled like many of the control objects, i.e.
///        MinimizationControls, based on namelists.
class Thermostat {
public:

  /// \brief Constructors for the Thermostat object.  Any information can be added via setter
  ///        functions after constructing the object.
  ///
  /// Overloaded:
  ///   - Construct a blank thermostat that applies no regulation
  ///   - Construct a specific type of thermostat with default settings
  ///   - Accept a type of thermostat and a temperature to maintain
  ///   - Accept a temperature evolution profile
  ///   - Accept a total number of atoms and a compartmentalization plan
  ///
  /// \param temperature_in    A flat temperature to apply at all times
  /// \{
  Thermostat(ThermostatKind kind_in = ThermostatKind::NONE,
             PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);
  
  Thermostat(int atom_count_in, ThermostatKind kind_in, double temperature_in,
             PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);

  Thermostat(int atom_count_in, ThermostatKind kind_in, double init_temperature_in,
             double final_temperature_in, int initial_evolution_step_in,
             int final_evolution_step_in, PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);

  Thermostat(int atom_count_in, ThermostatKind kind_in,
             const std::vector<double> &initial_temperatures_in,
             const std::vector<double> &final_temperatures_in,
             const std::vector<int2> &paritions_in,
             int initial_evolution_step_in = default_tstat_evo_window_start,
             int final_evolution_step_in = default_tstat_evo_window_end,
             PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);

  Thermostat(const AtomGraph *ag, ThermostatKind kind_in,
             double initial_temperature_in, double final_temperature_in,
             int initial_evolution_step_in, int final_evolution_step_in,
             PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);

  Thermostat(const AtomGraph &ag, ThermostatKind kind_in,
             double initial_temperature_in,
             double final_temperature_in,
             int initial_evolution_step_in, int final_evolution_step_in,
             PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);

  Thermostat(const AtomGraph *ag, ThermostatKind kind_in,
             double temperature_in = default_simulation_temperature,
             PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);
  
  Thermostat(const AtomGraph &ag, ThermostatKind kind_in,
             double temperature_in = default_simulation_temperature,
             PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);
  
  Thermostat(const AtomGraph *ag, ThermostatKind kind_in,
             const std::vector<double> &initial_temperatures_in,
             const std::vector<double> &final_temperatures_in,
             const std::vector<int2> &paritions_in,
             int initial_evolution_step_in = default_tstat_evo_window_start,
             int final_evolution_step_in = default_tstat_evo_window_end,
             PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);

  Thermostat(const AtomGraph &ag, ThermostatKind kind_in,
             const std::vector<double> &initial_temperatures_in,
             const std::vector<double> &final_temperatures_in,
             const std::vector<int2> &paritions_in,
             int initial_evolution_step_in = default_tstat_evo_window_start,
             int final_evolution_step_in = default_tstat_evo_window_end,
             PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);

  Thermostat(const AtomGraphSynthesis *poly_ag, ThermostatKind kind_in,
             double initial_temperature_in, double final_temperature_in,
             int initial_evolution_step_in, int final_evolution_step_in,
             PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);

  Thermostat(const AtomGraphSynthesis &poly_ag, ThermostatKind kind_in,
             double initial_temperature_in, double final_temperature_in,
             int initial_evolution_step_in, int final_evolution_step_in,
             PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);

  Thermostat(const AtomGraphSynthesis *poly_ag, ThermostatKind kind_in,
             double temperature_in = default_simulation_temperature,
             PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);
  
  Thermostat(const AtomGraphSynthesis &poly_ag, ThermostatKind kind_in,
             double temperature_in = default_simulation_temperature,
             PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);
  
  Thermostat(const AtomGraphSynthesis *poly_ag, ThermostatKind kind_in,
             const std::vector<double> &initial_temperatures_in,
             const std::vector<double> &final_temperatures_in,
             const std::vector<int2> &paritions_in,
             int initial_evolution_step_in = default_tstat_evo_window_start,
             int final_evolution_step_in = default_tstat_evo_window_end, 
             PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);

  Thermostat(const AtomGraphSynthesis &poly_ag, ThermostatKind kind_in,
             const std::vector<double> &initial_temperatures_in,
             const std::vector<double> &final_temperatures_in,
             const std::vector<int2> &paritions_in,
             int initial_evolution_step_in = default_tstat_evo_window_start,
             int final_evolution_step_in = default_tstat_evo_window_end,
             PrecisionModel cache_config_in = PrecisionModel::SINGLE,
             int random_seed_in = default_thermostat_random_seed,
             const GpuDetails &gpu = null_gpu);

  Thermostat(const AtomGraphSynthesis *poly_ag, const DynamicsControls &dyncon,
             const SystemCache &sc, const std::vector<int> &sc_origins,
             const GpuDetails &gpu = null_gpu);

  Thermostat(const AtomGraphSynthesis &poly_ag, const DynamicsControls &dyncon,
             const SystemCache &sc, const std::vector<int> &sc_origins,
             const GpuDetails &gpu = null_gpu);
  /// \}

  /// \brief The default copy and move constructors and assignment operators are best for this
  ///        object with no POINTER-kind Hybrid objects.
  ///
  /// \param original  The object to copy or move
  /// \param other     An object to clone or replace the present object with
  /// \{
  Thermostat(const Thermostat& original) = default;
  Thermostat(Thermostat&& original) = default;
  Thermostat& operator=(const Thermostat& original) = default;
  Thermostat& operator=(Thermostat&& original) = default;
  /// \}
  
  /// \brief Get the kind of thermostat
  ThermostatKind getKind() const;

  /// \brief Indicate whether this thermostat applies a common target temperature to all atoms, to
  ///        individual systems, or to some finer compartmentalization based on individual atoms.
  ThermostatPartition getBathPartitions() const;

  /// \brief Get the random number caching precision
  PrecisionModel getCacheConfiguration() const;
  
  /// \brief Get the total number of atoms that this thermostat can serve.
  int getAtomCount() const;

  /// \brief Get the step number according to this thermostat.  The thermostat will serve as an
  ///        official reference for the official simulation step number.
  int getStepNumber() const;

  /// \brief Get the random seed used to initialize Xoshiro256++ state vectors.
  int getRandomSeed() const;

  /// \brief Get the random number cache depth.  The quantity of random numbers that will be
  ///        produced and cached for each atom is three times this setting, as it implies values
  ///        for random perturbations in the Cartesian X, Y, and Z directions.
  int getRandomCacheDepth() const;

  /// \brief Get the step at which to begin changing the temperature from its initial value T(init)
  ///        to its final value T(final).  This is the last step that the thermostat will pull the
  ///        system towards a value of T(init).  Afterwards, the target will begin to change.
  int getInitialEvolutionStep() const;

  /// \brief Get the step at which the temperature evolution is expected to be complete.  This is
  ///        the step when the thermostat will begin pulling the system towards T(final), although
  ///        the system itself may or may not be there by this time.
  int getFinalEvolutionStep() const;

  /// \brief Get the generator state for a particular atom.  For retrieving large numbers of
  ///        states, one should used the abstract.
  ullint4 getGeneratorState(int atom_index,
                            HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
  /// \brief Get one of the cached random numbers.  This function returns the result as a
  ///        double-precision real, but if a single-precision result is queried it will be
  ///        converted to double and can then be put back in a float without changing its value.
  ///        For retrieving large numbers of results, one should used the abstract.
  ///
  /// \param atom_index  Index of the atom for which to get a result
  /// \param cache_row   Row of the cache to query (the 0th row offers a result for Cartesian X
  ///                    influence on the atom, the 4th row offers a result for the Carteisan Y
  ///                    influence on the atom one cycle in the future)
  /// \param tier        Indicate whether to take information from the CPU host or GPU device
  double getCachedRandomResult(int atom_index, int cache_row,
                               HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
  /// \brief Get the initial target temperature for this thermostat, in units of Kelvin.
  ///
  /// \param atom_index  Index of the atom of interest (only provided if there are different
  ///                    regulated temperatures for unique subsets of atoms)
  double getInitialTemperature(int atom_index = 0) const;

  /// \brief Get the end-point target temperature for this thermostat, in units of Kelvin.
  ///
  /// \param atom_index  Index of the atom of interest
  double getFinalTemperature(int atom_index = 0) const;

  /// \brief Get the spread of temperature targets for all segments of the thermostat at a specific
  ///        time step.
  ///
  /// \param step_number  The step at which to assess the temperature targets
  template <typename T> std::vector<T> getTemperatureSpread(int step_number = 0) const;

  /// \brief Get the partition boundaries between different segments of the thermostat as a vector
  ///        of integer tuples, with the partition lower and upper boundaries in the "x" and "y"
  ///        members, the number of constrained degrees of freedom in the "z" member, and the
  ///        number of unconstrained degrees of freedom in the "w" member of each tuple.
  std::vector<int4> getPartitionMap() const;
  
  /// \brief Get the current target temperature for this thermostat, in units of Kelvin, for a
  ///        specific atom.
  ///
  /// \param atom_index  Index of the atom of interest
  template <typename Tcalc> Tcalc getAtomTarget(int atom_index = 0) const;

  /// \brief Get the current target temperature for a specific partition within this thermostat,
  ///        in units of Kelvin.
  ///
  /// \param index  Index of the partition of interest
  template <typename Tcalc>
  Tcalc getPartitionTarget(int index = 0) const;

  /// \brief Get the number of constrained degrees of freedom for a particular partition of this
  ///        thermostat.
  ///
  /// \param index  Index of the partition of interest
  int getPartitionConstrainedDoF(int index = 0) const;

  /// \brief Get the number of degrees of freedom for a particular partition of this thermostat,
  ///        in absence of any constraints.
  ///
  /// \param index  Index of the partition of interest
  int getPartitionFreeDoF(int index = 0) const;

  /// \brief Get the Andersen reset frequency.
  int getAndersenResetFrequency() const;
  
  /// \brief Get the collision frequency for the Langevin thermostat (in units of femtoseconds).
  double getLangevinCollisionFrequency() const;

  /// \brief Get the "implicit" Langevin factor, the portion of the velocity update to apply
  ///        immediately after forces have been computed.
  double getLangevinImplicitFactor() const;

  /// \brief Get the "explicit" Langevin factor, the portion of the velocity update to apply
  ///        on the second half of the update, as new particle positions are being computed.
  double getLangevinExplicitFactor() const;
  
  /// \brief Get the time step, in units of femtoseconds.
  double getTimeStep() const;

  /// \brief Get the order whether to apply geometric constraints.
  ApplyConstraints constrainGeometry() const;

  /// \brief Get the convergence criterion for RATTLE bond length constraints.
  double getRattleTolerance() const;

  /// \brief Get the number of RATTLE iterations to attempt.
  int getRattleIterations() const;

  /// \brief Get a pointer to the object itself, useful if it has been passed to a function by
  //         const reference.
  const Thermostat* getSelfPointer() const;
  
  /// \brief Get the abstract.
  ///
  /// Overloaded:
  ///   - Get a reader for a const object
  ///   - Get a (partially) writeable abstract for a non-const object
  ///
  /// \param tier  Indicate whether to take pointers to memory on the CPU host or on the GPU device
  /// \{
  const ThermostatReader<double> dpData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  ThermostatWriter<double> dpData(HybridTargetLevel tier = HybridTargetLevel::HOST);
  const ThermostatReader<float> spData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  ThermostatWriter<float> spData(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}
  
  /// \brief Set the number of atoms for which the thermostat is responsible.
  ///
  /// \param atom_count_in  The total number of atoms to assign to this thermostat, whether all
  ///                       part of one system or a padded number collected within a synthesis
  /// \param new_seed       New random number generator seed for initializing the random cache (if
  ///                       applicable).  If set to -1, the object's original seed will be taken.
  ///                       Otherwise, the object will later record having been initialized with
  ///                       new_seed.
  /// \param gpu            Details of any available GPU.  If a real GPU is available, the random
  ///                       state vectors and cache will be initialized on the CPU host as well as
  ///                       the GPU device.
  void setAtomCount(int atom_count_in, int new_seed = -1, const GpuDetails &gpu = null_gpu);

  /// \brief Set the initial target temperature.
  ///
  /// Overloaded:
  ///   - Set the temperature to a single value (if the final temperatures are also consistent,
  ///     this will set common_temperature to TRUE and de-allocate any existing Hybrid data
  ///     related to unique temperature compartments).
  ///   - Set the temperatures of unique subsets of the atoms to a series of values based on a
  ///     compartmentalization
  ///
  /// \param init_temp_in   Temperature or temperatures that the thermostat shall start with, and
  ///                       continue to apply until the initial step in its evolution
  /// \param partitions_in  A new compartmentalization scheme (the same interpretation described
  ///                       in setCompartments() above applies)
  /// \{
  void setTemperature(double initial_temperature_in, double final_temperature_in = -1.0);
  void setTemperature(const std::vector<double> &initial_temperatures_in,
                      const std::vector<double> &final_temperatures_in,
                      const std::vector<int2> &partitions_in);
  /// \}

  /// \brief Check the compartmentalization if the thermostat regulates a synthesis of systems.
  ///
  /// \param poly_ag  The topology synthesis describing systems to regulate
  void checkCompartments(const AtomGraphSynthesis *poly_ag);

  /// \brief Set the step at which temperature evolution, away from T(init), shall begin.
  ///
  /// \param init_step_in  Step number at which temperature evolution begins.  Prior to this step,
  ///                      each bath's target temperature is its initial temperature.  Following
  ///                      this step, the target temperature of each bath becomes a linear mixture
  ///                      of the initial and final temperatures.
  /// \param init_step_in  Step number at which temperature evolution ends.  After this point, the
  ///                      final temperature of each bath is the target. 
  void setEvolutionWindow(int init_step_in, int final_step_in);
  
  /// \brief Set the random number cache configuration.
  ///
  /// \param cache_config_in  The manner in which the cache is to be configured, storing random
  ///                         numbers in single- or double-precision
  /// \param gpu              Details of the GPU that will handle future calculations with this
  ///                         thermostat.  This is critical for getting random state vectors
  ///                         initialized on the GPU.
  void setCacheConfiguration(PrecisionModel cache_config_in, const GpuDetails &gpu = null_gpu);
  
  /// \brief Set the random number cache depth.  This value is kept on a tight leash, as it can
  ///        easily lead to allocating too much memory.  A hundred-thousand atom system with cache
  ///        depth 8 leads to over 13.7MB of memory allocation.
  ///
  /// \param depth_in  The depth to take
  /// \param gpu       Details of the GPU that will handle future calculations with this
  ///                  thermostat.  This is critical for getting random state vectors initialized
  ///                  on the GPU.
  void setRandomCacheDepth(int depth_in, const GpuDetails &gpu = null_gpu);

  /// \brief Set the collision frequency for the Langevin thermostat.  This will also set the
  ///        scaled form of the value.
  ///
  /// \param langevin_frequency_in  The collision frequency to take, in events per femtosecond
  void setLangevinCollisionFrequency(double langevin_frequency_in);

  /// \brief Set the Andersen velocity reassignment frequency.  This will be taken from user input.
  ///
  /// \param andersen_frequency_in  The number of time steps between Andersen velocity resets
  void setAndersenResetFrequency(int andersen_frequency_in);
  
  /// \brief Set the time step.  This will also set the scaled form of the value.
  ///
  /// \param time_step_in  The time step to take, in units of femtoseconds
  void setTimeStep(double time_step_in);
  
  /// \brief Increase the step count by one.  This will advance the simulation's official step
  ///        counter.
  void incrementStep();

  /// \brief Decrease the step count by one.  This will cause the simulation's official step
  ///        counter to backtrack.
  void decrementStep();

  /// \brief Set the instruction for carrying out geometric constraints during dynamics.
  ///
  /// \param cnst_geometry_in  The order for constraints to be applied
  void setGeometryConstraints(ApplyConstraints cnst_geometry_in);

  /// \brief Set the convergence criterion for RATTLE.
  ///
  /// \param rattle_tolerance_in  The convergence criterion to set
  void setRattleTolerance(double rattle_tolerance_in);

  /// \brief Set the maximum number of iterations to take when attempting to converge RATTLE bond
  ///        length constraints.
  ///
  /// \param rattle_iterations_in  The number of iterations to permit
  void setRattleIterations(int rattle_iterations_in);
  
  /// \brief Validate the initial or final target temperature.
  void validateTemperature(double temperature_in);

  /// \brief Validate the initial and final steps for applying the temperature evolution.
  ///
  /// \param step_in  The step number to validate (it will be checked for a positive value)
  void validateEvolutionWindow();

  /// \brief Validate the random number cache depth.  A maxmimum cache depth of 15 (45 total
  ///        random numbers per atom) is enforced, and tightened in cases of exceptionally large
  ///        systems.
  void validateRandomCacheDepth();
  
  /// \brief Initialize the random state vectors of a thermostat, whether for Langevin or Andersen
  ///        thermostating.  This passes a call on to a general-purpose function for seeding the
  ///        generators and will recruit the GPU if available to expedite the process.
  ///
  /// \param new_seed      A new random number seed to apply, if different from the one used
  ///                      during construction of this Thermostat object
  /// \param scrub_cycles  Number of cycles of Xoshiro256++ generation to run in order to ensure
  ///                      that the newly seeded generators produce high-quality results
  /// \param tier          Indicate whether to seed random states and cached numbers on the CPU
  ///                      host alone, or on the GPU device.  Random number state vectors will be
  ///                      mirrored at both levels of memory if the GPU device layer is
  ///                      initialized, but the generators will not be advanced on the host.
  /// \param gpu           Specifications of the GPU in use (if applicable)
  void initializeRandomStates(int new_seed = default_thermostat_random_seed,
                              int scrub_cycles = 25,
                              HybridTargetLevel tier = HybridTargetLevel::HOST,
                              const GpuDetails &gpu = null_gpu);
  
  /// \brief Fill or refill the random number cache for a subset of the atoms, advancing the
  ///        thermostat's random number generators in the process.  This CPU-based function is
  ///        intended to mimic GPU functionality, but the GPU operations are expected to be fused
  ///        with other kernels.
  ///
  /// \param index_start    Starting index in the list of all atoms
  /// \param index_end      Upper bound of atoms for which to cache random numbers
  /// \param refresh_depth  The depth to which the cache is to be refreshed
  void refresh(size_t index_start, size_t index_end, int refresh_depth = 0);

#ifdef STORMM_USE_HPC
  /// \brief Upload the thermostat's data from the host to the HPC device.  Because the GPU is
  ///        used to initialize the vector of random number states (random_sv), this array is
  ///        uploaded and downloaded for synchronization upon initialization.  However, the initial
  ///        and final temperature arrays covering each atom require explicit upload.
  void upload();

  /// \brief Download the thermostat's data from the HPC device to the host.  This is needed to
  ///        synchronize progress made on the GPU if any processes are to be carried out by the
  ///        CPU host.
  void download();

  /// \brief Upload the temperature arrays and paritions only.
  void uploadPartitions();

  /// \brief Download the temperature arrays and partitions only.
  void downloadPartitions();
#endif
  
private:
  ThermostatKind kind;              ///< The type of thermostat
  ThermostatPartition bath_layout;  ///< Indicate whether the thermostat holds a single pair of
                                    ///<   initial and final temperatures for all atoms, a series
                                    ///<   of initial and final temperature pairs for all systems,
                                    ///<   or a series of initial and final temperature pairs for
                                    ///<   all atoms.
  PrecisionModel cache_config;      ///< Configuration for random numbers in the cache.  Because it
                                    ///<   would take a great deal of memory to allocate double- as
                                    ///<   well as single-precision arrays for any given quantity
                                    ///<   of cached random numbers, and a potential source of
                                    ///<   error for new pseudo-random numbers to be deposited in
                                    ///<   one array but mined from another, this mode setting
                                    ///<   locks the object into a particular configuration.
                                    ///<   Either precision model is compatible with arithmetic in
                                    ///<   single- or double-precision.  To store random numbers in
                                    ///<   single-precision and feed them into double-precision
                                    ///<   calculations is hardly a loss of information.
  int atom_count;                   ///< The total number of atoms, and thus the number of random
                                    ///<   state vectors, for which this thermostat is responsible
                                    ///<   (this times random_cache_depth gives the length of
                                    ///<   random_sv_bank)
  size_t padded_atom_count;         ///< The total number of atoms padded by the warp size
  int partition_count;              ///< The number of partitions, essentially the number of unique
                                    ///<   temperature baths that this thermostat comprises
  int step_number;                  ///< Number of the dynamics step, which should be consistent
                                    ///<   with any other record of the current simulation step.
                                    ///<   The thermostat, which will be included in any simulation
                                    ///<   even if set to NONE kind, is a good place to keep the
                                    ///<   official step number.
  int random_seed;                  ///< Seed for the first random state vector (position 0 in the
                                    ///<   random_sv_bank array)
  int random_cache_depth;           ///< Depth of the random cache
  int initial_evolution_step;       ///< The first step at which to initiate temperature evolution
  int final_evolution_step;         ///< Final step at which to temperature evolution is complete
  double initial_temperature;       ///< Initial temperature to apply to all particles
  double final_temperature;         ///< Final temperature to apply to all particles
  int andersen_frequency;           ///< The rate of Andersen velocity reassignments (valid only if
                                    ///<   an Andersen thermostat is in effect)
  double langevin_frequency;        ///< Frequency of collisions between all particles managed by
                                    ///<   this thermostat and a Langevin bath, in femtoseconds
  double time_step;                 ///< Simulation time step (the thermostat is also the official
                                    ///<   reference of this information in various C++ functions
                                    ///<   and GPU kernels), in femtoseconds
  ApplyConstraints cnst_geometry;   ///< Flag to indicate whether geometric constraints are in
                                    ///<   effect
  double rattle_tolerance;          ///< Convergence criterion for RATTLE bond length constraints
  int rattle_iterations;            ///< The maximum number of iterations to use when attempting to
                                    ///<   converge RATTLE constraint groups

  /// Temperatures to apply from step 0 to the initiation of any requested evolution, across
  ///   various compartments of the simulation.  Different compartments, i.e. systems, or
  ///   different components of each system, may have different starting temperatures, but the
  ///   evolution must proceed along the same schedule for all systems.
  /// \{
  Hybrid<double> initial_temperatures;
  Hybrid<float> sp_initial_temperatures;
  /// \}
  
  /// Temperatures to apply from the end of any requested evolution until the end of the simulation
  /// \{
  Hybrid<double> final_temperatures;
  Hybrid<float> sp_final_temperatures;
  /// \}

  /// Bounds array for various compartments of the particles affected by this thermostat.  The name
  /// is more general than "atom_starts" or "system_bounds" because this object is designed to
  /// serve one or many systems, and could, in theory, be used to thermostat different parts of one
  /// or all systems at different values.  Its "x" and "y" members hold the lower and upper bounds
  /// of atoms in the compartment, respectively, while its "z" and "w" members hold the constrained
  /// and unconstrained degrees of freedom for atoms in the compartment, respectively.
  Hybrid<int4> partitions;
  
  /// Xoshiro256++ state vectors for creating random numbers, one per atom of the simulation.
  /// These are stored in two ullint2 vectors to match the mechanics of templated random number
  /// generator arrays, all based on the fact that most GPUs' largest loads and stores are 128-bit.
  /// \{
  Hybrid<ullint2> random_state_vector_xy;
  Hybrid<ullint2> random_state_vector_zw;
  /// \}
  
  /// Space for the random state vector and random number cache
  /// \{
  Hybrid<double> random_cache;
  Hybrid<float> sp_random_cache;
  /// \}

  /// Pointers to the underlying topology or topology synthesis.  At most one of these pointers,
  /// perhaps neither, will be valid.
  /// \{
  AtomGraph *ag_pointer;
  AtomGraphSynthesis *poly_ag_pointer;
  /// \}
  
  /// \brief Allocate storage space for random number generator state vectors and cached random
  ///        numbers.  This will occur when the thermostat is set to Langevin or Andersen, or when
  ///        the atom count changes and the thermostat type is already Langevin or Andersen.
  ///
  /// \param new_seed  The random number cache, if it has any significance, will be initialized
  ///                  as it is reallocated.  This input provides the option of supplying a new
  ///                  random seed.  If set to -1, the object's original random seed will be
  ///                  taken.
  /// \param gpu       Details of the GPU to be used in the calculations.  If an actual GPU is
  ///                  available, initialization will be mirrored on the CPU host and GPU device.
  void allocateRandomStorage(int new_seed = -1, const GpuDetails &gpu = null_gpu);

  /// \brief Resize the temperature arrays as appropriate to support the stated partitioning.
  void resizeTemperatureArrays();
};

/// \brief Return the name of the thermostat choice (an enumerator string conversion function)
///
/// \param kind  The type of thermostat
std::string getThermostatName(ThermostatKind kind);

/// \brief Get the instantaneous target temperature from a thermostat for a particular atom in the
///        system or collection of systems regulated by the thermostat.
///
/// Overloaded:
///   - Provide the original object
///   - Provide the read-only or writeable thermostat abstract
///
/// \param tstw        The writeable thermostat abstract
/// \param tstr        The read-only thermostat abstract
/// \param tst         The original thermostat object
/// \param atom_index  Index of the atom of interest.  The partition in which the atom resides
///                    will be determined inside the function.
/// \{
template <typename Tcalc>
Tcalc getAtomTemperatureTarget(const ThermostatReader<Tcalc> &tstr, int atom_index = 0);

template <typename Tcalc>
Tcalc getAtomTemperatureTarget(const ThermostatWriter<Tcalc> &tstw, int atom_index = 0);

double getAtomTemperatureTarget(const Thermostat &tst, int atom_index = 0);
/// \}

/// \brief Get the instantaneous target temperature from a thermostat for a specific partition of
///        the atoms within the system or collection of systems regulated by the thermostat.
///
/// Overloaded:
///   - Provide the original object
///   - Provide the read-only or writeable thermostat abstract
///
/// \param tstw   The writeable thermostat abstract
/// \param tstr   The read-only thermostat abstract
/// \param tst    The original thermostat object
/// \param index  Index of the partition of interest.  This could be a specific system within a
///               synthesis, if the thermostat is regulating systems by specific temperatures.
/// \{
template <typename Tcalc>
Tcalc getPartitionTemperatureTarget(const ThermostatReader<Tcalc> &tstr, int index = 0);

template <typename Tcalc>
Tcalc getPartitionTemperatureTarget(const ThermostatWriter<Tcalc> &tstw, int index = 0);

double getPartitionTemperatureTarget(const Thermostat &tst, int index = 0);
/// \}

/// \brief Apply Andersen thermostating to reset (or initialize) all velocities in a system.
///
/// Overloaded:
///   - Set velocities in one system or many systems
///   - Provide abstracts or the original coordinate objects
///
/// \param xvel              Particle velocities along the Cartesian X axis
/// \param yvel              Particle velocities along the Cartesian Y axis
/// \param zvel              Particle velocities along the Cartesian Z axis
/// \param xvel_ovrf         Overflow bits for int95_t representations of Cartesian X velocities
/// \param yvel_ovrf         Overflow bits for int95_t representations of Cartesian Y velocities
/// \param zvel_ovrf         Overflow bits for int95_t representations of Cartesian Z velocities
/// \param psw               Abstract for a single system's coordinate object
/// \param ps                Coordinates (positions, velocities, and forces) for all particles in a
///                          single system
/// \param poly_psw          Abstract for a synthesis of coordinates (velocities will be assigned
///                          to the WHITE xvel, yvel, and zvel arrays of this abstract).
/// \param poly_ps           The synthesis of coordinates for many systems
/// \param inv_masses        Inverse masses of all particles, in units of Daltons or g/mol
/// \param cdk               Chemical details of a single system, including inverse masses
/// \param poly_auk          Atom update abstract for a coordinate synthesis, holding inverse
///                          masses
/// \param ag                Topology for a single system, containing atomic masses
/// \param poly_ag           Collated topologies for a synthesis of systems, holding inverse masses
/// \param tstr              Abstract of the thermostat, holding temperatures and a cache of random
///                          numbers.  The RNG cache will not be updated herein.  It is expected to
///                          be refilled on its own schedule, also held by settings in the
///                          thermostat.
/// \param tst               The thermostat controlling the velocity reassignment (see description
///                          above)
/// \param prec              Precision of calculations used in determining particle velocities.
///                          This will affect the representation of particle inverse masses and
///                          parameters (including the temperature) obtained from the thermostat,
///                          but the representation of random numbers in the thermostat is
///                          indicated by the thermostat itself.
/// \param vel_scale_factor  Velocity scaling factor for fixed-precision representations, obtained
///                          from the underlying coordinate synthesis
/// \{
template <typename Tcoord, typename Tcalc>
void andersenVelocityReset(Tcoord* xvel, Tcoord* yvel, Tcoord* zvel, int* xvel_ovrf,
                           int* yvel_ovrf, int* zvel_ovrf, const Tcalc* inv_masses,
                           const ThermostatReader<Tcalc> &tstr, Tcalc vel_scale_factor = 1.0);

template <typename Tcoord, typename Tcalc>
void andersenVelocityReset(Tcoord* xvel, Tcoord* yvel, Tcoord* zvel, const Tcalc* masses,
                           const ThermostatReader<Tcalc> &tstr, Tcalc vel_scale_factor = 1.0);

void andersenVelocityReset(PhaseSpaceWriter *psw, const ChemicalDetailsKit &cdk,
                           const ThermostatReader<double> &tstr);

void andersenVelocityReset(PhaseSpace *ps, const AtomGraph &ag, const Thermostat &tst);

template <typename Tcalc, typename Tcalc2, typename Tcalc4>
void andersenVelocityReset(PsSynthesisWriter *poly_psw,
                           const SyAtomUpdateKit<Tcalc, Tcalc2, Tcalc4> &poly_auk,
                           const ThermostatReader<Tcalc> &tstr);

#ifdef STORMM_USE_HPC
void andersenVelocityReset(PsSynthesisWriter *poly_psw,
                           const SyAtomUpdateKit<void, void, void> &poly_auk,
                           const ThermostatReader<void> &tstr,
                           PrecisionModel prec = PrecisionModel::SINGLE,
                           const GpuDetails &gpu = null_gpu);
#endif

void andersenVelocityReset(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                           const Thermostat *tst, PrecisionModel prec = PrecisionModel::SINGLE);

void andersenVelocityReset(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                           const Thermostat &tst, PrecisionModel prec = PrecisionModel::SINGLE);
/// \}

/// \brief Initialize velocities based on a particular thermostat configuration and temperature or
///        set of temperatures.  If the thermostat provided is configured for Andersen or Langevin
///        regulation, the random numbers expended for this initialization will be replaced by new
///        issues from their respective generators.  If the thermostat does not involve its own
///        random number streams, a new thermostat of similar heat bath configuration will be
///        produced internally, then destroyed after velocities are initialized.
///
/// Overloaded:
///   - Accept a single system.
///   - Accept a synthesis of multiple systems.
///
/// Descriptions of parameters follow from andersenVelocityReset(), above, in addition to:
///
/// \param new_seed  New random number seed guiding the generators that will provide the initial
///                  velocities, if the provided thermostat does not contain its own random number
///                  cache
/// \{
void velocityKickStart(PhaseSpace *ps, const AtomGraph *ag, Thermostat *tst,
                       const DynamicsControls &dyncon,
                       EnforceExactTemperature scale_temp = EnforceExactTemperature::YES);

void velocityKickStart(PhaseSpace *ps, const AtomGraph &ag, Thermostat *tst,
                       const DynamicsControls &dyncon,
                       EnforceExactTemperature scale_temp = EnforceExactTemperature::YES);

void velocityKickStart(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                       Thermostat *tst, const DynamicsControls &dyncon,
                       PrecisionModel prec = PrecisionModel::SINGLE,
                       EnforceExactTemperature scale_temp = EnforceExactTemperature::YES,
                       const GpuDetails &gpu = null_gpu);

void velocityKickStart(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                       Thermostat *tst, const DynamicsControls &dyncon,
                       PrecisionModel prec = PrecisionModel::SINGLE,
                       EnforceExactTemperature scale_temp = EnforceExactTemperature::YES,
                       const GpuDetails &gpu = null_gpu);
/// \}
  
} // namespace trajectory
} // namespace stormm

#include "thermostat.tpp"

#endif
