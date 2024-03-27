#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "Math/rounding.h"
#include "Math/series_ops.h"
#include "MolecularMechanics/kinetic.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Random/random.h"
#ifdef STORMM_USE_HPC
#include "Random/hpc_random.h"
#endif
#include "Reporting/error_format.h"
#include "Structure/rattle.h"
#include "Topology/atomgraph_analysis.h"
#include "UnitTesting/approx.h"
#include "thermostat.h"
#include "trim.h"

namespace stormm {
namespace trajectory {

using card::HybridKind;
using card::HybridTargetLevel;
using constants::time_step_factor;
using constants::CaseSensitivity;
using mm::computeTemperature;
using namelist::default_dynamics_time_step;
using namelist::default_geometry_constraint_behavior;
using namelist::default_rattle_max_iter;
using namelist::default_rattle_tolerance;
using parse::NumberFormat;
using parse::realToString;
using parse::strcmpCased;
using random::fillRandomCache;
using random::initXoshiro256ppArray;
using random::RandomAlgorithm;
using random::RandomNumberKind;
using random::Xoshiro256ppGenerator;
using stmath::incrementingSeries;
using stmath::roundUp;
using structure::rattlePositions;
using structure::rattleVelocities;
using structure::translateApplyConstraints;
using synthesis::SyNonbondedKit;
using topology::ConstraintKit;
using topology::getConstrainedDegreesOfFreedom;
  
//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const ThermostatKind kind_in, const PrecisionModel cache_config_in,
                       const int random_seed_in, const GpuDetails &gpu) :
    kind{kind_in},
    bath_layout{ThermostatPartition::COMMON},
    cache_config{cache_config_in},
    atom_count{0}, padded_atom_count{0}, partition_count{0}, step_number{0},
    random_seed{random_seed_in},
    random_cache_depth{default_thermostat_cache_depth},
    initial_evolution_step{default_tstat_evo_window_start},
    final_evolution_step{default_tstat_evo_window_end}, 
    initial_temperature{default_simulation_temperature},
    final_temperature{default_simulation_temperature},
    andersen_frequency{default_andersen_frequency},
    langevin_frequency{default_langevin_frequency},
    time_step{default_dynamics_time_step},
    cnst_geometry{translateApplyConstraints(std::string(default_geometry_constraint_behavior))},
    rattle_tolerance{default_rattle_tolerance},
    rattle_iterations{default_rattle_max_iter},
    initial_temperatures{HybridKind::ARRAY, "tstat_init_temp"},
    sp_initial_temperatures{HybridKind::ARRAY, "tstat_init_tempf"},
    final_temperatures{HybridKind::ARRAY, "tstat_final_temp"},
    sp_final_temperatures{HybridKind::ARRAY, "tstat_final_tempf"},
    partitions{HybridKind::ARRAY, "tstat_partitions"},
    random_state_vector_xy{HybridKind::ARRAY, "tstat_rng_xycache"},
    random_state_vector_zw{HybridKind::ARRAY, "tstat_rng_zwcache"},
    random_cache{HybridKind::ARRAY, "tstat_cache"},
    sp_random_cache{HybridKind::ARRAY, "tstat_cache_sp"},
    ag_pointer{nullptr}, poly_ag_pointer{nullptr}
{}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const int atom_count_in, const ThermostatKind kind_in,
                       const double temperature_in, const PrecisionModel cache_config_in,
                       const int random_seed_in, const GpuDetails &gpu) :
    Thermostat(kind_in, cache_config_in, random_seed_in, gpu)
{
  setAtomCount(atom_count_in);
  setTemperature(temperature_in, temperature_in);
  setEvolutionWindow(0, 0);
}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const int atom_count_in, const ThermostatKind kind_in,
                       const double initial_temperature_in, const double final_temperature_in,
                       const int initial_evolution_step_in, const int final_evolution_step_in,
                       const PrecisionModel cache_config_in, const int random_seed_in,
                       const GpuDetails &gpu) :
    Thermostat(kind_in, cache_config_in, random_seed_in, gpu)
{
  setAtomCount(atom_count_in);
  setTemperature(initial_temperature_in, final_temperature_in);
  setEvolutionWindow(initial_evolution_step_in, final_evolution_step_in);
}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const int atom_count_in, const ThermostatKind kind_in,
                       const std::vector<double> &initial_temperatures_in,
                       const std::vector<double> &final_temperatures_in,
                       const std::vector<int2> &partitions_in,
                       const int initial_evolution_step_in, const int final_evolution_step_in,
                       const PrecisionModel cache_config_in, const int random_seed_in,
                       const GpuDetails &gpu) :
    Thermostat(kind_in, cache_config_in, random_seed_in, gpu)
{
  setAtomCount(atom_count_in);
  setTemperature(initial_temperatures_in, final_temperatures_in, partitions_in);
  setEvolutionWindow(initial_evolution_step_in, final_evolution_step_in);
}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const AtomGraph *ag, const ThermostatKind kind_in,
                       const double initial_temperature_in, const double final_temperature_in,
                       const int initial_evolution_step_in, const int final_evolution_step_in,
                       const PrecisionModel cache_config_in, const int random_seed_in,
                       const GpuDetails &gpu) :
    Thermostat(kind_in, cache_config_in, random_seed_in, gpu)
{
  ag_pointer = const_cast<AtomGraph*>(ag);
  setAtomCount(ag->getAtomCount());
  setTemperature(initial_temperature_in, final_temperature_in);
  setEvolutionWindow(initial_evolution_step_in, final_evolution_step_in);
}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const AtomGraph &ag, const ThermostatKind kind_in,
                       const double initial_temperature_in, const double final_temperature_in,
                       const int initial_evolution_step_in, const int final_evolution_step_in,
                       const PrecisionModel cache_config_in, const int random_seed_in,
                       const GpuDetails &gpu) :
    Thermostat(ag.getSelfPointer(), kind_in, initial_temperature_in, final_temperature_in,
               initial_evolution_step_in, final_evolution_step_in, cache_config_in, random_seed_in,
               gpu)
{}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const AtomGraph *ag, const ThermostatKind kind_in,
                       const double temperature_in, const PrecisionModel cache_config_in,
                       const int random_seed_in, const GpuDetails &gpu) :
    Thermostat(ag, kind_in, temperature_in, temperature_in, 0, 0, cache_config_in, random_seed_in,
               gpu)
{}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const AtomGraph &ag, const ThermostatKind kind_in,
                       const double temperature_in, const PrecisionModel cache_config_in,
                       const int random_seed_in, const GpuDetails &gpu) :
    Thermostat(ag.getSelfPointer(), kind_in, temperature_in, temperature_in, 0, 0,
               cache_config_in, random_seed_in, gpu)
{}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const AtomGraph *ag, const ThermostatKind kind_in,
                       const std::vector<double> &initial_temperatures_in,
                       const std::vector<double> &final_temperatures_in,
                       const std::vector<int2> &partitions_in,
                       const int initial_evolution_step_in, const int final_evolution_step_in,
                       const PrecisionModel cache_config_in, const int random_seed_in,
                       const GpuDetails &gpu) :
    Thermostat(kind_in, cache_config_in, random_seed_in, gpu)
{
  ag_pointer = const_cast<AtomGraph*>(ag);
  setAtomCount(ag->getAtomCount());
  setTemperature(initial_temperatures_in, final_temperatures_in, partitions_in);
  setEvolutionWindow(initial_evolution_step_in, final_evolution_step_in);
}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const AtomGraph &ag, const ThermostatKind kind_in,
                       const std::vector<double> &initial_temperatures_in,
                       const std::vector<double> &final_temperatures_in,
                       const std::vector<int2> &partitions_in,
                       const int initial_evolution_step_in, const int final_evolution_step_in,
                       const PrecisionModel cache_config_in, const int random_seed_in,
                       const GpuDetails &gpu) :
    Thermostat(ag.getSelfPointer(), kind_in, initial_temperatures_in, final_temperatures_in,
               partitions_in, initial_evolution_step_in, final_evolution_step_in,
               cache_config_in, random_seed_in, gpu)
{}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const AtomGraphSynthesis *poly_ag, const ThermostatKind kind_in,
                       const double initial_temperature_in, const double final_temperature_in,
                       const int initial_evolution_step_in, const int final_evolution_step_in,
                       const PrecisionModel cache_config_in, const int random_seed_in,
                       const GpuDetails &gpu) :
    Thermostat(kind_in, cache_config_in, random_seed_in, gpu)
{
  poly_ag_pointer = const_cast<AtomGraphSynthesis*>(poly_ag);
  setAtomCount(poly_ag_pointer->getPaddedAtomCount());

  // Assume that all systems are to maintain the same temperature.  While it would seem acceptable
  // to leave the partitioning at "COMMON," each system has a unique number of degrees of freedom,
  // which the thermostat is also tasked to track.  Develop vectors of initial and final
  // temperatures as well as partitions for each system.
  const int nsys = poly_ag->getSystemCount();
  const std::vector<double> initial_temperatures_in(nsys, initial_temperature_in);
  const std::vector<double> final_temperatures_in(nsys, final_temperature_in);
  std::vector<int2> partitions_in(nsys);
  for (int i = 0; i < nsys; i++) {
    partitions_in[i].x = poly_ag->getAtomOffset(i);
    partitions_in[i].y = poly_ag->getAtomOffset(i) + poly_ag->getAtomCount(i);
  }
  setTemperature(initial_temperatures_in, final_temperatures_in, partitions_in);
  setEvolutionWindow(initial_evolution_step_in, final_evolution_step_in);
}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const AtomGraphSynthesis &poly_ag, const ThermostatKind kind_in,
                       const double initial_temperature_in, const double final_temperature_in,
                       const int initial_evolution_step_in, const int final_evolution_step_in,
                       const PrecisionModel cache_config_in, const int random_seed_in,
                       const GpuDetails &gpu) :
    Thermostat(poly_ag.getSelfPointer(), kind_in, initial_temperature_in, final_temperature_in,
               initial_evolution_step_in, final_evolution_step_in, cache_config_in, random_seed_in,
               gpu)
{}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const AtomGraphSynthesis *poly_ag, const ThermostatKind kind_in,
                       const double temperature_in, const PrecisionModel cache_config_in,
                       const int random_seed_in, const GpuDetails &gpu) :
    Thermostat(poly_ag, kind_in, temperature_in, temperature_in, 0, 0, cache_config_in,
               random_seed_in, gpu)
{}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const AtomGraphSynthesis &poly_ag, const ThermostatKind kind_in,
                       const double temperature_in, const PrecisionModel cache_config_in,
                       const int random_seed_in, const GpuDetails &gpu) :
    Thermostat(poly_ag.getSelfPointer(), kind_in, temperature_in, temperature_in, 0, 0,
               cache_config_in, random_seed_in, gpu)
{}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const AtomGraphSynthesis *poly_ag, const ThermostatKind kind_in,
                       const std::vector<double> &initial_temperatures_in,
                       const std::vector<double> &final_temperatures_in,
                       const std::vector<int2> &partitions_in,
                       const int initial_evolution_step_in, const int final_evolution_step_in,
                       const PrecisionModel cache_config_in, const int random_seed_in,
                       const GpuDetails &gpu) :
    Thermostat(poly_ag->getPaddedAtomCount(), kind_in, initial_temperatures_in,
               final_temperatures_in, partitions_in, initial_evolution_step_in,
               final_evolution_step_in, cache_config_in, random_seed_in, gpu)
{
  poly_ag_pointer = const_cast<AtomGraphSynthesis*>(poly_ag);
  setAtomCount(poly_ag_pointer->getPaddedAtomCount());
  setTemperature(initial_temperatures_in, final_temperatures_in, partitions_in);
  setEvolutionWindow(initial_evolution_step_in, final_evolution_step_in);
}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const AtomGraphSynthesis &poly_ag, const ThermostatKind kind_in,
                       const std::vector<double> &initial_temperatures_in,
                       const std::vector<double> &final_temperatures_in,
                       const std::vector<int2> &partitions_in,
                       const int initial_evolution_step_in, const int final_evolution_step_in,
                       const PrecisionModel cache_config_in, const int random_seed_in,
                       const GpuDetails &gpu) :
    Thermostat(poly_ag.getSelfPointer(), kind_in, initial_temperatures_in, final_temperatures_in,
               partitions_in, initial_evolution_step_in, final_evolution_step_in,
               cache_config_in, random_seed_in, gpu)
{}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const AtomGraphSynthesis *poly_ag, const DynamicsControls &dyncon,
                       const SystemCache &sc, const std::vector<int> &sc_origins,
                       const GpuDetails &gpu) :
  Thermostat(poly_ag, dyncon.getThermostatKind(), 0.0, 0.0, dyncon.getThermostatEvolutionStart(),
             dyncon.getThermostatEvolutionEnd(), dyncon.getThermostatCacheConfig(),
             dyncon.getThermostatSeed(), gpu)
{
  // Create a synthesis cache map and use it to build the thermostat partitions.
  SynthesisCacheMap scmap(sc_origins, sc.getSelfPointer(), poly_ag);

  // Create a baseline list of partitions based on the number of systems in the synthesis.  Make a
  // bounds array describing the extent of distinct partitions for each system in the synthesis.
  const int nsys = poly_ag->getSystemCount();
  int nprt = nsys;
  std::vector<int> system_partition_list_bounds = incrementingSeries(0, nsys + 1);
  std::vector<int2> system_partitions(nsys);
  std::vector<double> system_init_temps(nsys);
  std::vector<double> system_finl_temps(nsys);
  for (int i = 0; i < nsys; i++) {
    system_partitions[i].x = poly_ag->getAtomOffset(i);
    system_partitions[i].y = poly_ag->getAtomOffset(i) + poly_ag->getAtomCount(i);
  }
  
  // Work backwards from the label groups in the synthesis cache map to identify the system indices
  // to which each temperature applies.
  const std::vector<double>& init_user_temps = dyncon.getInitialTemperatureTargets();
  const std::vector<double>& finl_user_temps = dyncon.getFinalTemperatureTargets();
  const std::vector<std::string>& user_temp_labels = dyncon.getThermostatLabels();
  const std::vector<int>& user_temp_label_idx = dyncon.getThermostatLabelIndices();
  const std::vector<std::string>& user_temp_atom_masks = dyncon.getThermostatMasks();
  const int nt_user = dyncon.getThermostatLayerCount();
  for (int i = 0; i < nt_user; i++) {

    // Get the footprint of this directive across the synthesis.  The first entry should be a flat
    // directive to apply the default simulation temperature to all systems.  Future directives may
    // modify this condition, but in case they do not cover all systems in the synthesis the rest
    // will be thermostated at the default temperature.
    const double i_ut = init_user_temps[i];
    const double f_ut = finl_user_temps[i];
    if (strcmpCased(user_temp_labels[i], "all", CaseSensitivity::NO)) {
      
      // Apply the directive to all existing partitions
      for (int j = 0; j < nprt; j++) {

        // There may yet be an atom mask that further subdivides the systems.
        
        system_init_temps[j] = i_ut;
        system_finl_temps[j] = f_ut;
      }
    }
    else {
      const std::vector<int> ilbl_grp = scmap.getLabelGroup(user_temp_labels[i]);
      if (user_temp_label_idx[i] < 0) {
        
        // Apply the effect to all systems in the label group
        const int nlbl_sys = ilbl_grp.size();
        for (int j = 0; j < nlbl_sys; j++) {
          const int jsys = ilbl_grp[j];
          for (int k = system_partition_list_bounds[jsys];
               k < system_partition_list_bounds[jsys + 1]; k++) {
            system_init_temps[k] = i_ut;
            system_finl_temps[k] = f_ut;
          }
        }
      }
      else if (user_temp_label_idx[i] >= ilbl_grp.size()) {
        rtErr("The label group \"" + user_temp_labels[i] + "\" comprises " +
              std::to_string(ilbl_grp.size()) + ", not enough that index " +
              std::to_string(user_temp_label_idx[i]) + " would be valid.", "Thermostat");
      }
      else {
        const int jsys = ilbl_grp[user_temp_label_idx[i]];
        for (int k = system_partition_list_bounds[jsys];
             k < system_partition_list_bounds[jsys + 1]; k++) {
          system_init_temps[k] = i_ut;
          system_finl_temps[k] = f_ut;
        }
      }
    }
  }
  
  // The purpose of the above was to create a list of partitions and temperature profiles for them.
  // Fill the object's member arrays with this information.  Other details, including the pointer
  // to the topology synthesis, the evolution end points, and the overall atom count will have
  // already been set by the simpler delegated constructor used in the initialization phase.
  setTemperature(system_init_temps, system_finl_temps, system_partitions);
  
  // Set critical integration-related parameters based on the availabe dynamics controls object.
  time_step = dyncon.getTimeStep();
  cnst_geometry = dyncon.constrainGeometry();
  rattle_tolerance = dyncon.getRattleTolerance();
  rattle_iterations = dyncon.getRattleIterations();
  random_seed = dyncon.getThermostatSeed();
  setRandomCacheDepth(dyncon.getThermostatCacheDepth(), gpu);
  andersen_frequency = dyncon.getAndersenFrequency();
  langevin_frequency = dyncon.getLangevinFrequency();
}

//-------------------------------------------------------------------------------------------------
Thermostat::Thermostat(const AtomGraphSynthesis &poly_ag, const DynamicsControls &dyncon,
                       const SystemCache &sc, const std::vector<int> &sc_origins,
                       const GpuDetails &gpu) :
    Thermostat(poly_ag.getSelfPointer(), dyncon, sc, sc_origins, gpu)
{}

//-------------------------------------------------------------------------------------------------
void Thermostat::allocateRandomStorage(const int new_seed, const GpuDetails &gpu) {

  // Do not allocate extra memory which will not be used by certain thermostats.  When seeding
  // particle velocities, the process looks very much like applying an Andersen thermostat once at
  // the outset of the simulation.  In that case, it is faster to have a single generator, driven
  // by a single CPU thread, follow a similar strategy of creating up to 1024 generators and then
  // seeding a series of atoms' velocities with each of them.  That velocity initialization process
  // will get its own function, deferring to the array of generators for each atom in the case of
  // Andersen and Langevin thermostats, or using a temporary array of generators, outside of the
  // Thermostat object, in the cases of Berendesen or no thermostat.
  switch (kind) {
  case ThermostatKind::NONE:
  case ThermostatKind::BERENDSEN:
    random_state_vector_xy.resize(0);
    random_state_vector_zw.resize(0);
    random_cache.resize(0);
    sp_random_cache.resize(0);
    break;
  case ThermostatKind::ANDERSEN:
  case ThermostatKind::LANGEVIN:
    {
      const size_t padded_atom_count_zu = padded_atom_count;
      const size_t rc_depth_zu = random_cache_depth * 3;
      random_state_vector_xy.resize(padded_atom_count_zu);
      random_state_vector_zw.resize(padded_atom_count_zu);
      switch (cache_config) {
      case PrecisionModel::DOUBLE:
        random_cache.resize(padded_atom_count_zu * rc_depth_zu);
        sp_random_cache.resize(0);
        break;
      case PrecisionModel::SINGLE:
        random_cache.resize(0);
        sp_random_cache.resize(padded_atom_count_zu * rc_depth_zu);
        break;
      }

      // Automatically initialize the random state vectors after resizing them.
      const int actual_seed = (new_seed == -1) ? random_seed : new_seed;
      if (gpu == null_gpu) {
        initializeRandomStates(actual_seed, 25, HybridTargetLevel::HOST, gpu);
      }
      else {
#ifdef STORMM_USE_HPC
        initializeRandomStates(actual_seed, 25, HybridTargetLevel::DEVICE, gpu);
#else
        initializeRandomStates(actual_seed, 25, HybridTargetLevel::HOST, gpu);
#endif
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void Thermostat::resizeTemperatureArrays() {
  switch (bath_layout) {
  case ThermostatPartition::COMMON:
    initial_temperatures.resize(0);
    sp_initial_temperatures.resize(0);
    final_temperatures.resize(0);
    sp_final_temperatures.resize(0);
    break;
  case ThermostatPartition::SYSTEMS:
    initial_temperatures.resize(partition_count);
    sp_initial_temperatures.resize(partition_count);
    final_temperatures.resize(partition_count);
    sp_final_temperatures.resize(partition_count);
    break;
  case ThermostatPartition::ATOMS:
    initial_temperatures.resize(padded_atom_count);
    sp_initial_temperatures.resize(padded_atom_count);
    final_temperatures.resize(padded_atom_count);
    sp_final_temperatures.resize(padded_atom_count);
    break;
  }
}
  
//-------------------------------------------------------------------------------------------------
ThermostatKind Thermostat::getKind() const {
  return kind;
}

//-------------------------------------------------------------------------------------------------
ThermostatPartition Thermostat::getBathPartitions() const {
  return bath_layout;
}

//-------------------------------------------------------------------------------------------------
PrecisionModel Thermostat::getCacheConfiguration() const {
  return cache_config;
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getStepNumber() const {
  return step_number;
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getRandomSeed() const {
  return random_seed;
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getRandomCacheDepth() const {
  return random_cache_depth;
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getInitialEvolutionStep() const {
  return initial_evolution_step;
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getFinalEvolutionStep() const {
  return final_evolution_step;
}

//-------------------------------------------------------------------------------------------------
ullint4 Thermostat::getGeneratorState(const int atom_index, const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const ullint2 xy_state = random_state_vector_xy.readHost(atom_index);
      const ullint2 zw_state = random_state_vector_zw.readHost(atom_index);
      return { xy_state.x, xy_state.y, zw_state.x, zw_state.y };
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const ullint2 xy_state = random_state_vector_xy.readDevice(atom_index);
      const ullint2 zw_state = random_state_vector_zw.readDevice(atom_index);
      return { xy_state.x, xy_state.y, zw_state.x, zw_state.y };
    }
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double Thermostat::getCachedRandomResult(const int atom_index, const int cache_row,
                                         const HybridTargetLevel tier) const {
  const size_t pos = (static_cast<size_t>(cache_row) * static_cast<size_t>(padded_atom_count)) +
                     static_cast<size_t>(atom_index);
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (cache_config) {
    case PrecisionModel::DOUBLE:
      return random_cache.readHost(pos);
    case PrecisionModel::SINGLE:
      return sp_random_cache.readHost(pos);
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (cache_config) {
    case PrecisionModel::DOUBLE:
      return random_cache.readDevice(pos);
    case PrecisionModel::SINGLE:
      return sp_random_cache.readDevice(pos);
    }
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double Thermostat::getInitialTemperature(const int index) const {
  switch (bath_layout) {
  case ThermostatPartition::COMMON:
    return initial_temperature;
  case ThermostatPartition::SYSTEMS:
  case ThermostatPartition::ATOMS:
    return initial_temperatures.readHost(index);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double Thermostat::getFinalTemperature(const int index) const {
  switch (bath_layout) {
  case ThermostatPartition::COMMON:
    return final_temperature;
  case ThermostatPartition::SYSTEMS:
  case ThermostatPartition::ATOMS:
    return final_temperatures.readHost(index);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int4> Thermostat::getPartitionMap() const {
  return partitions.readHost();
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getPartitionConstrainedDoF(const int index) const {
  const int4 prtn = partitions.readHost(index);
  return prtn.z;
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getPartitionFreeDoF(const int index) const {
  const int4 prtn = partitions.readHost(index);
  return prtn.w;
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getAndersenResetFrequency() const {
  return andersen_frequency;
}

//-------------------------------------------------------------------------------------------------
double Thermostat::getLangevinCollisionFrequency() const {
  return langevin_frequency;
}

//-------------------------------------------------------------------------------------------------
double Thermostat::getLangevinImplicitFactor() const {
  return 1.0 / (1.0 + (langevin_frequency * 0.5 * time_step));
}

//-------------------------------------------------------------------------------------------------
double Thermostat::getLangevinExplicitFactor() const {
  return 1.0 - (langevin_frequency * 0.5 * time_step);
}

//-------------------------------------------------------------------------------------------------
double Thermostat::getTimeStep() const {
  return time_step;
}

//-------------------------------------------------------------------------------------------------
ApplyConstraints Thermostat::constrainGeometry() const {
  return cnst_geometry;
}

//-------------------------------------------------------------------------------------------------
double Thermostat::getRattleTolerance() const {
  return rattle_tolerance;
}

//-------------------------------------------------------------------------------------------------
int Thermostat::getRattleIterations() const {
  return rattle_iterations;
}

//-------------------------------------------------------------------------------------------------
const Thermostat* Thermostat::getSelfPointer() const {
  return this;
}

//-------------------------------------------------------------------------------------------------
const ThermostatReader<double> Thermostat::dpData(const HybridTargetLevel tier) const {
  return ThermostatReader<double>(kind, bath_layout, atom_count, padded_atom_count,
                                  partition_count, step_number, random_cache_depth,
                                  initial_evolution_step, final_evolution_step,
                                  initial_temperature, final_temperature, time_step,
                                  (cnst_geometry == ApplyConstraints::YES), rattle_tolerance,
                                  rattle_iterations, andersen_frequency, langevin_frequency,
                                  getLangevinImplicitFactor(), getLangevinExplicitFactor(),
                                  partitions.data(tier), initial_temperatures.data(tier),
                                  final_temperatures.data(tier), random_state_vector_xy.data(tier),
                                  random_state_vector_zw.data(tier), cache_config,
                                  random_cache.data(tier), sp_random_cache.data(tier));
}

//-------------------------------------------------------------------------------------------------
ThermostatWriter<double> Thermostat::dpData(const HybridTargetLevel tier) {
  return ThermostatWriter<double>(kind, bath_layout, atom_count, padded_atom_count,
                                  partition_count, step_number, random_cache_depth,
                                  initial_evolution_step, final_evolution_step,
                                  initial_temperature, final_temperature, time_step,
                                  (cnst_geometry == ApplyConstraints::YES), rattle_tolerance,
                                  rattle_iterations, andersen_frequency, langevin_frequency,
                                  getLangevinImplicitFactor(), getLangevinExplicitFactor(),
                                  partitions.data(tier), initial_temperatures.data(tier),
                                  final_temperatures.data(tier), random_state_vector_xy.data(tier),
                                  random_state_vector_zw.data(tier), cache_config,
                                  random_cache.data(tier), sp_random_cache.data(tier));
}

//-------------------------------------------------------------------------------------------------
const ThermostatReader<float> Thermostat::spData(const HybridTargetLevel tier) const {
  return ThermostatReader<float>(kind, bath_layout, atom_count, padded_atom_count,
                                 partition_count, step_number, random_cache_depth,
                                 initial_evolution_step, final_evolution_step,
                                 initial_temperature, final_temperature, time_step,
                                 (cnst_geometry == ApplyConstraints::YES), rattle_tolerance,
                                 rattle_iterations, andersen_frequency, langevin_frequency,
                                 getLangevinImplicitFactor(), getLangevinExplicitFactor(),
                                 partitions.data(tier), sp_initial_temperatures.data(tier),
                                 sp_final_temperatures.data(tier),
                                 random_state_vector_xy.data(tier),
                                 random_state_vector_zw.data(tier), cache_config,
                                 random_cache.data(tier), sp_random_cache.data(tier));
}

//-------------------------------------------------------------------------------------------------
ThermostatWriter<float> Thermostat::spData(const HybridTargetLevel tier) {
  return ThermostatWriter<float>(kind, bath_layout, atom_count, padded_atom_count,
                                 partition_count, step_number, random_cache_depth,
                                 initial_evolution_step, final_evolution_step,
                                 initial_temperature, final_temperature, time_step,
                                 (cnst_geometry == ApplyConstraints::YES), rattle_tolerance,
                                 rattle_iterations, andersen_frequency, langevin_frequency,
                                 getLangevinImplicitFactor(), getLangevinExplicitFactor(),
                                 partitions.data(tier), sp_initial_temperatures.data(tier),
                                 sp_final_temperatures.data(tier),
                                 random_state_vector_xy.data(tier),
                                 random_state_vector_zw.data(tier), cache_config,
                                 random_cache.data(tier), sp_random_cache.data(tier));
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setAtomCount(const int atom_count_in, const int new_seed, const GpuDetails &gpu) {
  atom_count = atom_count_in;
  padded_atom_count = roundUp(atom_count, warp_size_int);
  allocateRandomStorage(new_seed, gpu);
}
  
//-------------------------------------------------------------------------------------------------
void Thermostat::setTemperature(const double initial_temperature_in,
                                const double final_temperature_in) {
  bath_layout = ThermostatPartition::COMMON;
  partition_count = 1;
  partitions.resize(1);

  // Determine the number of constraints
  int ncnst_dof, nfree_dof;
  if (ag_pointer != nullptr) {
    const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
    ncnst_dof = cdk.cnst_dof;
    nfree_dof = cdk.free_dof;
  }
  else if (poly_ag_pointer != nullptr) {

    // The proper temperature setting for a synthesis of systems is below, even if the
    // vectors of temperatures and partitions have been set automatically.
    rtErr("A synthesis must have defined partitions, even if all systems adhere to the same "
          "temperatures.", "Thermostat", "setTemperature");
  }
  else {

    // Assume that there are no constraints in effect
    ncnst_dof = (3 * atom_count) - 6;
    nfree_dof = (3 * atom_count) - 6;
  }
  partitions.putHost({ 0, atom_count, ncnst_dof, nfree_dof }, 0);
  resizeTemperatureArrays();
  initial_temperature = initial_temperature_in;
  if (final_temperature_in < 0.0) {
    final_temperature = initial_temperature;
  }
  validateTemperature(initial_temperature);
  validateTemperature(final_temperature);
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setTemperature(const std::vector<double> &initial_temperatures_in,
                                const std::vector<double> &final_temperatures_in,
                                const std::vector<int2> &partitions_in) {

  // Check that a synthesis of systems if properly supported with at least one partition for each
  // system, and that no partition spans more that one system.
  std::vector<int> system_homes(padded_atom_count, 0);
  std::vector<bool> critical_atoms(padded_atom_count, false);
  partition_count = partitions_in.size();
  if (poly_ag_pointer != nullptr) {
    if (partitions_in.size() < poly_ag_pointer->getSystemCount()) {
      rtErr("There are not enough thermostat partitions (" + std::to_string(partitions_in.size()) +
            ") to account for each system of the synthesis (" +
            std::to_string(poly_ag_pointer->getSystemCount()) + ").", "Thermostat",
            "setTemperature");
    }
    
    // Use the non-bonded kit, which contains the individual system offsets and atom counts.  Fix
    // any errors in the compartmentalization such that all compartments apply to at most one
    // system out of the synthesis but all systems are covered according to the list of
    // temperatures laid out iin the original instructions.
    const SyNonbondedKit<double,
                         double2> synbk = poly_ag_pointer->getDoublePrecisionNonbondedKit();
    const size_t last_sys = synbk.nsys - 1;
    std::vector<bool> covered_atoms(padded_atom_count, false);
    for (int i = 0; i < synbk.nsys; i++) {
      const int jlim = synbk.atom_offsets[i] + synbk.atom_counts[i];
      for (int j = synbk.atom_offsets[i]; j < jlim; j++) {
        critical_atoms[j] = true;
        system_homes[j] = i;
      }
    }
    bool partitions_ok = true;
    for (int i = 0; i < partition_count; i++) {
      const int llim = partitions_in[i].x;
      const int hlim = partitions_in[i].y;
      int current_system;
      for (int j = llim; j < hlim; j++) {
        if (j == llim) {
          current_system = system_homes[j];
        }
        covered_atoms[j] = true;
        partitions_ok = (partitions_ok && critical_atoms[j] && system_homes[j] == current_system);
      }
      if (partitions_ok == false) {
        rtErr("Partition index " + std::to_string(i) + " covers more than one system in the "
              "underlying synthesis, or an invalid atom.");
      }
    }
    int atoms_missing = 0;
    for (int i = 0; i < padded_atom_count; i++) {
      atoms_missing += (critical_atoms[i] && (covered_atoms[i] == false));
    }
    if (atoms_missing > 0) {
      rtErr("A total of " + std::to_string(atoms_missing) + " atoms are not covered by any "
            "thermostat partition.", "Thermostat", "setTemperature");
    }
  }
  else if (ag_pointer != nullptr) {
    for (int i = atom_count; i < padded_atom_count; i++) {
      critical_atoms[i] = false;      
    }
  }
  
  // Contingency: no temperatures were actually provided.  This can only be valid if there is a
  // number of atoms provided with no topology pointer, or the thermostat is based on a pointer to
  // a single topology.  The other case, that the thermostat is based on a synthesis of topologies,
  // is trapped above.  Transform into a COMMON bath layout with the default simulation
  // temperature.
  if (initial_temperatures_in.size() == 0 && final_temperatures_in.size() == 0) {
    setTemperature(default_simulation_temperature, default_simulation_temperature);
    return;
  }
  else if (poly_ag_pointer == nullptr &&
           initial_temperatures_in.size() == 1 && final_temperatures_in.size() == 1) {
    setTemperature(initial_temperatures_in[0], final_temperatures_in[0]);
    return;
  }
  else if (initial_temperatures_in.size() > 0 && final_temperatures_in.size() == 0) {
    setTemperature(initial_temperatures_in, initial_temperatures_in, partitions_in);
    return;
  }
  else if (initial_temperatures_in.size() != final_temperatures_in.size()) {

    // Fail immediately if the initial and final temperature vectors do not align.
    rtErr("Different numbers of initial and final temperatures (" +
          std::to_string(initial_temperatures_in.size()) + ", " +
          std::to_string(final_temperatures_in.size()) + ") indicate an impossible layout.",
          "Thermostat");
  }
  if (initial_temperatures_in.size() != partitions_in.size()) {
    rtErr("Partitions for heat baths (" + std::to_string(partitions_in.size()) +
          ") must be provided in accordance with the number of temperatures provided (" +
          std::to_string(initial_temperatures_in.size()) + ").", "Thermostat", "setTemperature");
  }
  
  // Reaching this point in the function indicates that there are two vectors of temperatures, of
  // equal size, determining how to lay out the heat baths.  There is a valid array of
  // compartments, partitions sanning the atom content, but it is not known at this level whether
  // those partitions refer to demarcations within one or more systems or strictly to the
  // boundaries between systems.  If they are boundaries between systems (each system has its own
  // distinct heat bath), then a degree of optimization with respect to memory and bandwidth is
  // possible.  However, this must be conveyed prior to the call, by setting the bath_layout to
  // the special "SYSTEMS" enumeration.  Otherwise, the more expensive but completely general
  // "ATOMS" enumeration will be taken.
  partitions.resize(partitions_in.size());
  int4* part_ptr = partitions.data();
  if (ag_pointer != nullptr) {
    for (int i = 0; i < partition_count; i++) {
      const int cnst_dof = getConstrainedDegreesOfFreedom(ag_pointer, partitions_in[i].x,
                                                          partitions_in[i].y);
      const int free_dof = (3 * (partitions_in[i].y - partitions_in[i].x)) - 6;
      part_ptr[i] = { partitions_in[i].x, partitions_in[i].y, cnst_dof, free_dof };
    }
    bath_layout = (partition_count > 1) ? ThermostatPartition::ATOMS : ThermostatPartition::COMMON;
  }
  else if (poly_ag_pointer != nullptr) {
    const SyNonbondedKit<double,
                         double2> synbk = poly_ag_pointer->getDoublePrecisionNonbondedKit();
    for (int i = 0; i < partition_count; i++) {
      const int sysid = system_homes[partitions_in[i].x];
      const AtomGraph *iag = poly_ag_pointer->getSystemTopologyPointer(sysid);
      const int sys_offset = synbk.atom_offsets[sysid];
      const int cnst_dof = getConstrainedDegreesOfFreedom(iag, partitions_in[i].x - sys_offset,
                                                          partitions_in[i].y - sys_offset);
      const int free_dof = (3 * (partitions_in[i].y - partitions_in[i].x)) - 6;
      part_ptr[i] = { partitions_in[i].x, partitions_in[i].y, cnst_dof, free_dof };
    }
    bath_layout = (partition_count > poly_ag_pointer->getSystemCount()) ?
                  ThermostatPartition::ATOMS : ThermostatPartition::SYSTEMS;
  }
  else {

    // For a collection of atoms with no defined topology, no constraints can be assumed.
    for (int i = 0; i < partition_count; i++) {
      const int free_dof = (3 * (partitions_in[i].y - partitions_in[i].x)) - 6;
      part_ptr[i] = { partitions_in[i].x, partitions_in[i].y, free_dof, free_dof };
    }
  }
  resizeTemperatureArrays();
  switch (bath_layout) {
  case ThermostatPartition::COMMON:
    break;
  case ThermostatPartition::SYSTEMS:
    {
      initial_temperatures.putHost(initial_temperatures_in);
      final_temperatures.putHost(final_temperatures_in);
      const std::vector<float> sp_init_temp(initial_temperatures_in.begin(),
                                            initial_temperatures_in.end());
      const std::vector<float> sp_final_temp(final_temperatures_in.begin(),
                                             final_temperatures_in.end());
      sp_initial_temperatures.putHost(sp_init_temp);
      sp_final_temperatures.putHost(sp_final_temp);
    }
    break;
  case ThermostatPartition::ATOMS:
    {
      double* init_temp_ptr = initial_temperatures.data();
      double* finl_temp_ptr = final_temperatures.data();
      float* sp_init_temp_ptr = sp_initial_temperatures.data();
      float* sp_finl_temp_ptr = sp_final_temperatures.data();
      for (int i = 0; i < partition_count; i++) {
        for (int j = partitions_in[i].x; j < partitions_in[i].y; j++) {
          init_temp_ptr[j] = initial_temperatures_in[i];
          finl_temp_ptr[j] = final_temperatures_in[i];
          sp_init_temp_ptr[j] = initial_temperatures_in[i];
          sp_finl_temp_ptr[j] = final_temperatures_in[i];
        }
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setEvolutionWindow(const int init_step_in, const int final_step_in) {
  initial_evolution_step = init_step_in;
  final_evolution_step = final_step_in;
  validateEvolutionWindow();
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setCacheConfiguration(const PrecisionModel cache_config_in,
                                       const GpuDetails &gpu) {
  cache_config = cache_config_in;
  allocateRandomStorage(random_seed, gpu);
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setRandomCacheDepth(const int depth_in, const GpuDetails &gpu) {
  random_cache_depth = depth_in;
  validateRandomCacheDepth();
  allocateRandomStorage(random_seed, gpu);
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setAndersenResetFrequency(const int andersen_frequency_in) {
  andersen_frequency = std::max(andersen_frequency_in, 0);
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setLangevinCollisionFrequency(const double langevin_frequency_in) {
  langevin_frequency = langevin_frequency_in;
  if (langevin_frequency > 0.020) {

    // The Langvein collision frequency is expressed in units of inverse femtoseconds, in contrast
    // to a code like Amber where the frequency is expressed in inverse picoseconds.  As a result,
    // Amber inputs will give collision frequencies 1000x larger than intended.
    rtErr("A Langevin frequency of " +
          realToString(langevin_frequency, 7, 4, NumberFormat::STANDARD_REAL) + "/fs is far too "
          "high for realistic thermocoupling.  The input may have been mistakenly assumed to have "
          "units of inverse picoseconds, but its units of inverse femtoseconds demand that such "
          "inputs be divided by 1000.", "Thermostat", "setLangevinCollisionFrequency");
  }
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setTimeStep(const double time_step_in) {
  time_step = time_step_in;
}

//-------------------------------------------------------------------------------------------------
void Thermostat::incrementStep() {
  step_number += 1;
}

//-------------------------------------------------------------------------------------------------
void Thermostat::decrementStep() {
  step_number -= 1;
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setGeometryConstraints(const ApplyConstraints cnst_geometry_in) {
  cnst_geometry = cnst_geometry_in;
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setRattleTolerance(const double rattle_tolerance_in) {
  rattle_tolerance = rattle_tolerance_in;
}

//-------------------------------------------------------------------------------------------------
void Thermostat::setRattleIterations(const int rattle_iterations_in) {
  rattle_iterations = rattle_iterations_in;
}

//-------------------------------------------------------------------------------------------------
void Thermostat::validateTemperature(const double temperature_in) {
  if (temperature_in < 0.0 || temperature_in >= 1.0e5) {
    rtErr("A temperature of " + realToString(temperature_in, 11, 4, NumberFormat::STANDARD_REAL) +
          " is not a sensible choice.", "Thermostat", "validateTemperature");
  }
}

//-------------------------------------------------------------------------------------------------
void Thermostat::validateEvolutionWindow() {
  if (initial_evolution_step < 0) {
    initial_evolution_step = 0;
  }
  if (final_evolution_step <= initial_evolution_step) {
    final_evolution_step = initial_evolution_step;
  }
}

//-------------------------------------------------------------------------------------------------
void Thermostat::validateRandomCacheDepth() {
  if (random_cache_depth > maximum_thermostat_cache_depth) {

    // A silent change occurs here.  The only effect that a bad cache depth could have is that the
    // simulation runs at an imperceptibly different rate than the one that might be expected.
    random_cache_depth = maximum_thermostat_cache_depth;
  }
  switch (kind) {
  case ThermostatKind::NONE:
  case ThermostatKind::BERENDSEN:
    random_cache_depth = 0;
    break;
  case ThermostatKind::ANDERSEN:
  case ThermostatKind::LANGEVIN:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void Thermostat::initializeRandomStates(const int new_seed, const int scrub_cycles,
                                        const HybridTargetLevel tier, const GpuDetails &gpu) {

  // Always replace the random seed with the new seed, but return immediately if there is no random
  // number cache to fill.
  random_seed = new_seed;
  switch (kind) {
  case ThermostatKind::NONE:
  case ThermostatKind::BERENDSEN:
    return;
  case ThermostatKind::ANDERSEN:
  case ThermostatKind::LANGEVIN:
    break;
  }

  // Initialize random states on the CPU host or GPU device.
  switch (tier) {
  case HybridTargetLevel::HOST:
    initXoshiro256ppArray(&random_state_vector_xy, &random_state_vector_zw, new_seed,
                          scrub_cycles);
    switch (cache_config) {
    case PrecisionModel::DOUBLE:
      fillRandomCache(&random_state_vector_xy, &random_state_vector_zw, &random_cache,
                      padded_atom_count, random_cache_depth * 3, RandomAlgorithm::XOSHIRO_256PP,
                      RandomNumberKind::GAUSSIAN, 0, atom_count);
      break;
    case PrecisionModel::SINGLE:
      fillRandomCache(&random_state_vector_xy, &random_state_vector_zw, &sp_random_cache,
                      padded_atom_count, random_cache_depth * 3, RandomAlgorithm::XOSHIRO_256PP,
                      RandomNumberKind::GAUSSIAN, 0, atom_count);
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    initXoshiro256ppArray(&random_state_vector_xy, &random_state_vector_zw, new_seed, scrub_cycles,
                          gpu);
    switch (cache_config) {
    case PrecisionModel::DOUBLE:
      fillRandomCache(&random_state_vector_xy, &random_state_vector_zw, &random_cache,
                      padded_atom_count, random_cache_depth * 3, RandomAlgorithm::XOSHIRO_256PP,
                      RandomNumberKind::GAUSSIAN, 0, atom_count, gpu);
      break;
    case PrecisionModel::SINGLE:
      fillRandomCache(&random_state_vector_xy, &random_state_vector_zw, &sp_random_cache,
                      padded_atom_count, random_cache_depth * 3, RandomAlgorithm::XOSHIRO_256PP,
                      RandomNumberKind::GAUSSIAN, 0, atom_count, gpu);
      break;
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void Thermostat::refresh(const size_t index_start, const size_t index_end,
                         const int refresh_depth) {
  const int actual_depth = (refresh_depth > 0) ? std::min(refresh_depth, random_cache_depth) :
                                                 random_cache_depth;
  
  // Return immediately with a warning if a developer tries to issue random numbers for thermostats
  // that do not use them.
  switch (kind) {
  case ThermostatKind::NONE:
  case ThermostatKind::BERENDSEN:
    return;
  case ThermostatKind::ANDERSEN:
  case ThermostatKind::LANGEVIN:
    break;
  }
  switch (cache_config) {
  case PrecisionModel::DOUBLE:
    fillRandomCache(random_state_vector_xy.data(), random_state_vector_zw.data(),
                    random_cache.data(), atom_count, actual_depth * 3,
                    RandomAlgorithm::XOSHIRO_256PP, RandomNumberKind::GAUSSIAN, 0, atom_count);
    break;
  case PrecisionModel::SINGLE:
    fillRandomCache(random_state_vector_xy.data(), random_state_vector_zw.data(),
                    sp_random_cache.data(), atom_count, actual_depth * 3,
                    RandomAlgorithm::XOSHIRO_256PP, RandomNumberKind::GAUSSIAN, 0, atom_count);
    break;
  }
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void Thermostat::upload() {
  partitions.upload();
  initial_temperatures.upload();
  sp_initial_temperatures.upload();
  final_temperatures.upload();
  sp_final_temperatures.upload();
  random_state_vector_xy.upload();
  random_state_vector_zw.upload();
}

//-------------------------------------------------------------------------------------------------
void Thermostat::download() {
  partitions.download();
  initial_temperatures.download();
  sp_initial_temperatures.download();
  final_temperatures.download();
  sp_final_temperatures.download();
  random_state_vector_xy.download();
  random_state_vector_zw.download();
}

//-------------------------------------------------------------------------------------------------
void Thermostat::uploadPartitions() {
  partitions.upload();
  initial_temperatures.upload();
  sp_initial_temperatures.upload();
  final_temperatures.upload();
  sp_final_temperatures.upload();
}

//-------------------------------------------------------------------------------------------------
void Thermostat::downloadPartitions() {
  partitions.download();
  initial_temperatures.download();
  sp_initial_temperatures.download();
  final_temperatures.download();
  sp_final_temperatures.download();
}
#endif

//-------------------------------------------------------------------------------------------------
std::string getThermostatName(const ThermostatKind kind) {
  switch (kind) {
  case ThermostatKind::NONE:
    return std::string("NONE");
  case ThermostatKind::ANDERSEN:
    return std::string("ANDERSEN");
  case ThermostatKind::LANGEVIN:
    return std::string("LANGEVIN");
  case ThermostatKind::BERENDSEN:
    return std::string("BERENDSEN");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double getAtomTemperatureTarget(const Thermostat &tst, const int atom_index) {
  return tst.getAtomTarget<double>(atom_index);
}

//-------------------------------------------------------------------------------------------------
double getPartitionTemperatureTarget(const Thermostat &tst, const int index) {
  return tst.getPartitionTarget<double>(index);
}

//-------------------------------------------------------------------------------------------------
void andersenVelocityReset(PhaseSpaceWriter *psw, const ChemicalDetailsKit &cdk,
                           const ThermostatReader<double> &tstr) {

  // Check that the thermostat matches the coordinate object for size
  if (tstr.natom != psw->natom) {
    rtErr("The atom content of the thermostat (" + std::to_string(tstr.natom) + ") and the "
          "coordinate set (" + std::to_string(psw->natom) + ") do not match.",
          "andersenVelocityReset");
  }
  andersenVelocityReset<double, double>(psw->xvel, psw->yvel, psw->zvel, nullptr, nullptr, nullptr,
                                        cdk.inv_masses, tstr);
}

//-------------------------------------------------------------------------------------------------
void andersenVelocityReset(PhaseSpace *ps, const AtomGraph *ag, const Thermostat *tst) {
  PhaseSpaceWriter psw = ps->data();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  const ThermostatReader<double> tstr = tst->dpData();
  andersenVelocityReset(&psw, cdk, tstr);
}

//-------------------------------------------------------------------------------------------------
void andersenVelocityReset(PhaseSpace *ps, const AtomGraph &ag, const Thermostat &tst) {
  PhaseSpaceWriter psw = ps->data();
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  const ThermostatReader<double> tstr = tst.dpData();
  andersenVelocityReset(&psw, cdk, tstr);
}

//-------------------------------------------------------------------------------------------------
void andersenVelocityReset(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                           const Thermostat *tst, const PrecisionModel prec) {

  // Check that the thermostat matches the coordinate object for size
  const int pps_natom = poly_ps->getPaddedAtomCount();
  if (tst->getAtomCount() != pps_natom) {
    rtErr("The atom content of the thermostat (" + std::to_string(tst->getAtomCount()) +
          ") and the coordinate set (" + std::to_string(pps_natom) + ") do not match.",
          "andersenVelocityReset");
  }

  // Implement thermostating in the requested precision model
  PsSynthesisWriter poly_psw = poly_ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyAtomUpdateKit<double,
                            double2,
                            double4> poly_auk = poly_ag->getDoublePrecisionAtomUpdateKit();
      const ThermostatReader<double> tstr = tst->dpData();
      andersenVelocityReset<double, double2, double4>(&poly_psw, poly_auk, tstr);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyAtomUpdateKit<float,
                            float2, float4> poly_auk = poly_ag->getSinglePrecisionAtomUpdateKit();
      const ThermostatReader<float> tstr = tst->spData();
      andersenVelocityReset<float, float2, float4>(&poly_psw, poly_auk, tstr);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void andersenVelocityReset(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                           const Thermostat &tst, const PrecisionModel prec) {
  andersenVelocityReset(poly_ps, poly_ag.getSelfPointer(), tst.getSelfPointer(), prec);
}

//-------------------------------------------------------------------------------------------------
void velocityKickStart(PhaseSpace *ps, const AtomGraph *ag, Thermostat *tst,
                       const DynamicsControls &dyncon, const EnforceExactTemperature scale_temp) {

  // If constraints are in effect, it is only possible to apply valid velocities to the constrained
  // groups if they are at their constrained lengths.  This will generate a degree of velocity on
  // each of the particles, but this effect is ignored as the velocities are about to be reset.
  // Use the non-current point in the coordinate time cycle to develop the constrained positions,
  // then copy the results back to the current point.
  PhaseSpaceWriter psw = ps->data();
  switch (tst->constrainGeometry()) {
  case ApplyConstraints::YES:
    {
      const ConstraintKit<double> cnk = ag->getDoublePrecisionConstraintKit();
      for (int i = 0; i < psw.natom; i++) {
        psw.xalt[i] = psw.xcrd[i];
        psw.yalt[i] = psw.ycrd[i];
        psw.zalt[i] = psw.zcrd[i];
      }
      rattlePositions<double,
                      double>(psw.xalt, psw.yalt, psw.zalt, psw.vxalt, psw.vyalt, psw.vzalt,
                              psw.xcrd, psw.ycrd, psw.zcrd, cnk, 1.0, dyncon.getRattleTolerance(),
                              dyncon.getRattleIterations(), dyncon.getCpuRattleMethod());
      for (int i = 0; i < psw.natom; i++) {
        psw.xcrd[i] = psw.xalt[i];
        psw.ycrd[i] = psw.yalt[i];
        psw.zcrd[i] = psw.zalt[i];
      }
    }
    break;
  case ApplyConstraints::NO:
    break;
  }

  // Apply the basic velocity reset
  switch (tst->getKind()) {
  case ThermostatKind::ANDERSEN:
  case ThermostatKind::LANGEVIN:

    // Set the velocities of particles using the set of random numbers at the top of the cache.
    // Replace those numbers immediately.
    andersenVelocityReset(ps, ag, tst);
    tst->refresh(0, tst->getAtomCount(), 1);
    break;
  case ThermostatKind::NONE:
  case ThermostatKind::BERENDSEN:
    {
      // Create a new thermostat with random numbers to populate the velocity arrays.  Destroy that
      // thermostat once the velocities are set.  The dynamics control object's main thermostat
      // random seed can be used, as the thermostat is not otherwise based on any random sequence.
      const std::vector<double> temperatures = tst->getTemperatureSpread<double>(0);
      const std::vector<int4> tmp_partitions = tst->getPartitionMap();
      const size_t nprt = tmp_partitions.size();
      std::vector<int2> partition_limits(nprt);
      for (size_t i = 0; i < nprt; i++) {
        partition_limits[i].x = tmp_partitions[i].x;
        partition_limits[i].y = tmp_partitions[i].y;
      }
      Thermostat disposable(tst->getAtomCount(), ThermostatKind::ANDERSEN, temperatures,
                            temperatures, partition_limits);
      disposable.setRandomCacheDepth(1);
      disposable.initializeRandomStates(dyncon.getThermostatSeed());
      andersenVelocityReset(ps, ag, &disposable);
    }
    break;
  }

  // Remove net translational and, if applicable, net rotational momentum
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  removeMomentum<double, double, double>(psw.xcrd, psw.ycrd, psw.zcrd, nullptr, nullptr, nullptr,
                                         psw.xvel, psw.yvel, psw.zvel, nullptr, nullptr, nullptr,
                                         cdk.masses, psw.unit_cell, psw.natom);
  
  // Stash the unconstrained velocities in the BLACK time point.  As before with positions in
  // the presence of constraints, it is assumed that these velocities are irrelevant and will be
  // overwritten in the subsequent dynamics step.
  switch (scale_temp) {
  case EnforceExactTemperature::YES:
    for (int i = 0; i < psw.natom; i++) {
      psw.vxalt[i] = psw.xvel[i];
      psw.vyalt[i] = psw.yvel[i];
      psw.vzalt[i] = psw.zvel[i];
    }
    break;
  case EnforceExactTemperature::NO:
    break;
  }
  
  // Apply constraints, working with the coordinate object's WHITE time cycle point as the
  // developmental velocities, as well as the reference positions.
  switch (tst->constrainGeometry()) {
  case ApplyConstraints::YES:
    {
      const ConstraintKit<double> cnk = ag->getDoublePrecisionConstraintKit();
      rattleVelocities<double,
                       double>(psw.xvel, psw.yvel, psw.zvel, psw.xcrd, psw.ycrd, psw.zcrd,
                               cnk, tst->getTimeStep(), dyncon.getRattleTolerance(),
                               dyncon.getRattleIterations(), dyncon.getCpuRattleMethod());
    }
    break;
  case ApplyConstraints::NO:
    break;
  }
  
  // Enforce the exact temperature, if requested.
  switch (scale_temp) {
  case EnforceExactTemperature::YES:
    {
      const ThermostatWriter<double> tstw = tst->dpData();
      switch (tstw.layout) {
      case ThermostatPartition::COMMON:
        {
          const double tcr = computeTemperature<double,
                                                double,
                                                double>(psw.xvel, psw.yvel, psw.zvel, nullptr,
                                                        nullptr, nullptr, cdk.masses, psw.natom,
                                                        (tstw.cnst_geom) ? cdk.cnst_dof :
                                                                           cdk.free_dof,
                                                        pow(2.0, 32), 1.0);
          const double vscl = sqrt(tstw.init_temperature / tcr);
          for (int i = 0; i < psw.natom; i++) {
            psw.vxalt[i] *= vscl;
            psw.vyalt[i] *= vscl;
            psw.vzalt[i] *= vscl;
          }
        }
        break;
      case ThermostatPartition::SYSTEMS:
      case ThermostatPartition::ATOMS:
        {
          //  Compute the number of degrees of freedom in any given partition of the system.
          for (int i = 0; i < tstw.npart; i++) {
            const int4 pinfo = tstw.partition_bounds[i];
            const size_t llim = pinfo.x;
            const size_t hlim = pinfo.y;
            const int prt_natom = hlim - llim;
            const double tcr = computeTemperature<double,
                                                  double,
                                                  double>(&psw.xvel[llim], &psw.yvel[llim],
                                                          &psw.zvel[llim], nullptr, nullptr,
                                                          nullptr, &cdk.masses[llim], prt_natom,
                                                          (tstw.cnst_geom) ? pinfo.z : pinfo.w,
                                                          pow(2.0, 32), 1.0);
            const double vscl = sqrt(tstw.init_temperatures[i] / tcr);
            for (size_t j = llim; j < hlim; j++) {
              psw.vxalt[j] *= vscl;
              psw.vyalt[j] *= vscl;
              psw.vzalt[j] *= vscl;
            }
          }
        }
        break;
      }

      // Copy the velocities back and repeat the constraints.  The result will have a quantity of
      // kinetic energy sufficient to give it the exact temperature requested.
      for (int i = 0; i < psw.natom; i++) {
        psw.xvel[i] = psw.vxalt[i];
        psw.yvel[i] = psw.vyalt[i];
        psw.zvel[i] = psw.vzalt[i];
      }
      if (tstw.cnst_geom) {
        const ConstraintKit<double> cnk = ag->getDoublePrecisionConstraintKit();
        rattleVelocities<double,
                         double>(psw.xvel, psw.yvel, psw.zvel, psw.xcrd, psw.ycrd, psw.zcrd,
                                 cnk, tst->getTimeStep(), dyncon.getRattleTolerance(),
                                 dyncon.getRattleIterations(), dyncon.getCpuRattleMethod());
      }
    }    
    break;
  case EnforceExactTemperature::NO:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void velocityKickStart(PhaseSpace *ps, const AtomGraph &ag, Thermostat *tst,
                       const DynamicsControls &dyncon, const EnforceExactTemperature scale_temp) {
  velocityKickStart(ps, ag.getSelfPointer(), tst, dyncon, scale_temp);
}

//-------------------------------------------------------------------------------------------------
void velocityKickStart(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                       Thermostat *tst, const DynamicsControls &dyncon, const PrecisionModel prec,
                       const EnforceExactTemperature scale_temp, const GpuDetails &gpu) {
  if (gpu == null_gpu) {
    PsSynthesisWriter poly_psw = poly_ps->data();
    const SyAtomUpdateKit<double,
                          double2, double4> auk_d = poly_ag->getDoublePrecisionAtomUpdateKit();
    const SyAtomUpdateKit<float,
                          float2, float4> auk_f = poly_ag->getSinglePrecisionAtomUpdateKit();
    const ThermostatWriter<double> tstw_d = tst->dpData();
    const ThermostatWriter<float> tstw_f = tst->spData();

    // Apply positional constraints
    if (tstw_d.cnst_geom) {
      rattlePositions(poly_ps, poly_ag, prec, tst->getTimeStep(), tst->getRattleTolerance(),
                      tst->getRattleIterations());
    }

    // Perform the basic velocity reset.
    switch (tst->getKind()) {
    case ThermostatKind::ANDERSEN:
    case ThermostatKind::LANGEVIN:

      // Set the velocities of particles using the set of random numbers at the top of the cache.
      // Replace those numbers immediately.
      andersenVelocityReset(poly_ps, poly_ag, tst, prec);
      tst->refresh(0, tst->getAtomCount(), 1);
      break;
    case ThermostatKind::NONE:
    case ThermostatKind::BERENDSEN:
      {
        // Create a new thermostat with random numbers to populate the velocity arrays.  Destroy
        // that thermostat once the velocities are set.  The dynamics control object's main
        // thermostat random seed can be used, as the thermostat is not otherwise based on any
        // random sequence.
        const std::vector<double> temperatures = tst->getTemperatureSpread<double>(0);
        const std::vector<int4> tmp_partitions = tst->getPartitionMap();
        const size_t nprt = tmp_partitions.size();
        std::vector<int2> partition_limits(nprt);
        for (size_t i = 0; i < nprt; i++) {
          partition_limits[i].x = tmp_partitions[i].x;
          partition_limits[i].y = tmp_partitions[i].y;
        }
        Thermostat disposable(tst->getAtomCount(), ThermostatKind::ANDERSEN, temperatures,
                              temperatures, partition_limits);
        disposable.setRandomCacheDepth(1);
        disposable.initializeRandomStates(dyncon.getThermostatSeed());
        andersenVelocityReset(poly_ps, poly_ag, &disposable);
      }
      break;
    }

    // Remove net translational and, if applicable, rotational momentum from all systems.  This
    // involves multiple reductions over all atoms in each system and therefore multiple kernel
    // calls in the GPU code.
    removeMomentum(poly_ps, poly_ag, prec);

    // As in the single-system case, stash the unconstrained velocities in the BLACK time
    // point.  Again, it is assumed that these velocities are irrelevant and will be overwritten
    // in the subsequent dynamics step.
    switch (scale_temp) {
    case EnforceExactTemperature::YES:
      for (int i = 0; i < poly_psw.system_count; i++) {
        const int llim = poly_psw.atom_starts[i];
        const int hlim = llim + poly_psw.atom_counts[i];
        for (int j = llim; j < hlim; j++) {
          poly_psw.vxalt[j] = poly_psw.xvel[j];
          poly_psw.vyalt[j] = poly_psw.yvel[j];
          poly_psw.vzalt[j] = poly_psw.zvel[j];
        }
        switch (prec) {
        case PrecisionModel::DOUBLE:
          for (int j = llim; j < hlim; j++) {
            poly_psw.vxalt_ovrf[j] = poly_psw.xvel_ovrf[j];
            poly_psw.vyalt_ovrf[j] = poly_psw.yvel_ovrf[j];
            poly_psw.vzalt_ovrf[j] = poly_psw.zvel_ovrf[j];
          }
          break;
        case PrecisionModel::SINGLE:
          break;
        }
      }
      break;
    case EnforceExactTemperature::NO:
      break;
    }
    
    // Apply velocity constraints, working with the coordinate object's WHITE time cycle point
    // as the developmental velocities, as well as the reference positions.
    if (tstw_d.cnst_geom) {
      rattleVelocities(poly_ps, poly_ag, prec, tst->getTimeStep(), tst->getRattleTolerance(),
                       tst->getRattleIterations());
    }

    // Enforce the exact temperature, if requested.  Reset the velocities once more if temperature
    // rescaling was performed.
    switch (scale_temp) {
    case EnforceExactTemperature::YES:
      {
        switch (tstw_d.layout) {
        case ThermostatPartition::COMMON:

          // A thermostat serving a coordinate synthesis cannot pool all of its systems into a
          // single partition.
          break;
        case ThermostatPartition::SYSTEMS:
        case ThermostatPartition::ATOMS:
          {
            // Compute the temperature, based on the (perhaps constrained) WHITE velocities
            // and the number of degrees of freedom in any given partition of the system.
            for (int i = 0; i < tstw_d.npart; i++) {

              // Both the single- and double-precision thermostat abstracts will contain the same
              // partition bounds and counts of the degrees of freedom.
              const int4 pinfo = tstw_d.partition_bounds[i];
              const size_t llim = pinfo.x;
              const size_t hlim = pinfo.y;
              const int prt_natom = hlim - llim;
              switch (prec) {
              case PrecisionModel::DOUBLE:
                {
                  const double tcr =
                    computeTemperature<llint,
                                       double,
                                       double>(&poly_psw.xvel[llim], &poly_psw.yvel[llim],
                                               &poly_psw.zvel[llim], &poly_psw.xvel_ovrf[llim],
                                               &poly_psw.yvel_ovrf[llim],
                                               &poly_psw.zvel_ovrf[llim], &auk_d.masses[llim],
                                               prt_natom, (tstw_d.cnst_geom) ? pinfo.z : pinfo.w,
                                               pow(2.0, 32), poly_psw.inv_vel_scale_f);
                  const double vscl = sqrt(tstw_d.init_temperatures[i] / tcr);
                  for (size_t j = llim; j < hlim; j++) {
                    const double vxa = hostInt95ToDouble(poly_psw.vxalt[j],
                                                         poly_psw.vxalt_ovrf[j]) * vscl;
                    const double vya = hostInt95ToDouble(poly_psw.vyalt[j],
                                                         poly_psw.vyalt_ovrf[j]) * vscl;
                    const double vza = hostInt95ToDouble(poly_psw.vzalt[j],
                                                         poly_psw.vzalt_ovrf[j]) * vscl;
                    const int95_t ivxa = hostDoubleToInt95(vxa);
                    const int95_t ivya = hostDoubleToInt95(vya);
                    const int95_t ivza = hostDoubleToInt95(vza);
                    poly_psw.vxalt[j] = ivxa.x;
                    poly_psw.vyalt[j] = ivya.x;
                    poly_psw.vzalt[j] = ivza.x;
                    poly_psw.vxalt_ovrf[j] = ivxa.y;
                    poly_psw.vyalt_ovrf[j] = ivya.y;
                    poly_psw.vzalt_ovrf[j] = ivza.y;
                  }
                }
                break;
              case PrecisionModel::SINGLE:
                {
                  const float tcr =
                    computeTemperature<llint,
                                       float,
                                       float>(&poly_psw.xvel[llim], &poly_psw.yvel[llim],
                                              &poly_psw.zvel[llim], nullptr, nullptr, nullptr,
                                              &auk_f.masses[llim], prt_natom,
                                              (tstw_d.cnst_geom) ? pinfo.z : pinfo.w,
                                              powf(2.0, 32), poly_psw.inv_vel_scale_f);
                  const float vscl = sqrtf(tstw_f.init_temperatures[i] / tcr);
                  for (size_t j = llim; j < hlim; j++) {
                    poly_psw.vxalt[j] = llround(static_cast<float>(poly_psw.vxalt[j]) * vscl);
                    poly_psw.vyalt[j] = llround(static_cast<float>(poly_psw.vyalt[j]) * vscl);
                    poly_psw.vzalt[j] = llround(static_cast<float>(poly_psw.vzalt[j]) * vscl);
                  }
                }
                break;
              }
            }
          }
          break;
        }

        // Copy the velocities back and repeat the constraints.  The result will have a quantity of
        // kinetic energy sufficient to give it the exact temperature requested.
        for (int i = 0; i < poly_psw.system_count; i++) {
          const int llim = poly_psw.atom_starts[i];
          const int hlim = llim + poly_psw.atom_counts[i];
          for (int j = llim; j < hlim; j++) {
            poly_psw.xvel[j] = poly_psw.vxalt[j];
            poly_psw.yvel[j] = poly_psw.vyalt[j];
            poly_psw.zvel[j] = poly_psw.vzalt[j];
          }
          switch (prec) {
          case PrecisionModel::DOUBLE:
            for (int j = llim; j < hlim; j++) {
              poly_psw.xvel_ovrf[j] = poly_psw.vxalt_ovrf[j];
              poly_psw.yvel_ovrf[j] = poly_psw.vyalt_ovrf[j];
              poly_psw.zvel_ovrf[j] = poly_psw.vzalt_ovrf[j];
            }
            break;
          case PrecisionModel::SINGLE:
            break;
          }
        }

        // Constrain the velocities once more.  The strategy was to subtract off the system's net
        // momentum.  If constraints were to be applied, these velocities for the non-migrating,
        // non-rotating system were stashed in the BLACK arrays of the coordinate synthesis.
        // Whether constraints were then applied or not, the system's actual temperature was
        // computed and a scaling factor needed to bring the system to the target temperature was
        // computed.  This scaling was then applied to the stashed velocities, which were then
        // copied back to the current velocities in the WHITE stage of the time cycle.  If
        // constraints were applied before testing to obtain the correct scaling factor, they must
        // now be applied again to the rescaled, raw velocities that now stand in the WHITE stage
        // of the time cycle.
        if (tstw_d.cnst_geom) {
          rattleVelocities(poly_ps, poly_ag, prec, tst->getTimeStep(), tst->getRattleTolerance(),
                           tst->getRattleIterations());
        }
      }    
      break;
    case EnforceExactTemperature::NO:
      break;
    }
  }
#ifdef STORMM_USE_HPC
  else {
    PsSynthesisWriter poly_psw = poly_ps->data(HybridTargetLevel::DEVICE);
  }
#endif

}

//-------------------------------------------------------------------------------------------------
void velocityKickStart(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                       Thermostat *tst, const DynamicsControls &dyncon, const PrecisionModel prec,
                       const EnforceExactTemperature scale_temp, const GpuDetails &gpu) {
  velocityKickStart(poly_ps, poly_ag.getSelfPointer(), tst, dyncon, prec, scale_temp, gpu);
}

} // namespace trajectory
} // namespace stormm
