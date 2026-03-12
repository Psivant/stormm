// -*-c++-*-
#include "stormm_pybridge.h"

#include <algorithm>
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

#include "Constants/generalized_born.h"
#include "MolecularMechanics/dynamics.h"
#include "Namelists/nml_dynamics.h"
#include "Potential/energy_enumerators.h"
#include "Potential/scorecard.h"
#include "Potential/static_exclusionmask.h"
#include "Restraints/restraint_apparatus.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Topology/atomgraph.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/trajectory_enumerators.h"

namespace {

// Import the STORMM types wrapped by this C ABI layer.  Keeping aliases local avoids leaking
// STORMM namespaces into extern "C" declarations while reducing verbosity in wrapper code.
using stormm::energy::ScoreCard;
using stormm::energy::StateVariable;
using stormm::energy::StaticExclusionMask;
using stormm::generalized_born_defaults::NeckGeneralizedBornTable;
using stormm::mm::dynamics;
using stormm::namelist::DynamicsControls;
using stormm::restraints::RestraintApparatus;
using stormm::synthesis::AtomGraphSynthesis;
using stormm::synthesis::PhaseSpaceSynthesis;
using stormm::topology::AtomGraph;
using stormm::trajectory::PhaseSpace;
using stormm::trajectory::Thermostat;
using stormm::trajectory::ThermostatKind;
using stormm::trajectory::TrajectoryKind;

// Bridge errors are stored per-thread so concurrent Python use can report diagnostics correctly.
thread_local std::string stormm_last_error;

void setError(const std::string &msg) {
  stormm_last_error = msg;
}

void clearError() {
  stormm_last_error.clear();
}

// Validate opaque pointers originating from Python before reinterpreting as C++ handles.
template <typename T>
T* requireHandle(void* handle, const char* name) {
  if (handle == nullptr) {
    throw std::runtime_error(std::string("Null handle for ") + name + ".");
  }
  return reinterpret_cast<T*>(handle);
}

template <typename T>
const T* requireHandle(const void* handle, const char* name) {
  if (handle == nullptr) {
    throw std::runtime_error(std::string("Null handle for ") + name + ".");
  }
  return reinterpret_cast<const T*>(handle);
}

// Convert C++ exceptions to integer return codes and record a string message for Python.
template <typename Fn>
int wrapVoid(Fn &&fn) {
  try {
    fn();
    clearError();
    return 0;
  }
  catch (const std::exception &e) {
    setError(e.what());
  }
  catch (...) {
    setError("Unknown C++ exception.");
  }
  return -1;
}

template <typename Fn, typename RetType>
RetType wrapValue(Fn &&fn, RetType error_value) {
  try {
    const RetType result = fn();
    clearError();
    return result;
  }
  catch (const std::exception &e) {
    setError(e.what());
  }
  catch (...) {
    setError("Unknown C++ exception.");
  }
  return error_value;
}

// Translate compact Python thermostat codes into STORMM enum values.
ThermostatKind toThermostatKind(const int kind_code) {
  switch (kind_code) {
  case 0:
    return ThermostatKind::NONE;
  case 1:
    return ThermostatKind::ANDERSEN;
  case 2:
    return ThermostatKind::LANGEVIN;
  case 3:
    return ThermostatKind::BERENDSEN;
  default:
    throw std::runtime_error("Unknown thermostat code " + std::to_string(kind_code) +
                             ". Valid values are 0 (none), 1 (andersen), 2 (langevin), "
                             "3 (berendsen).");
  }
}

// Translate compact Python trajectory codes into STORMM enum values.
TrajectoryKind toTrajectoryKind(const int kind_code) {
  switch (kind_code) {
  case 0:
    return TrajectoryKind::POSITIONS;
  case 1:
    return TrajectoryKind::VELOCITIES;
  case 2:
    return TrajectoryKind::FORCES;
  default:
    throw std::runtime_error("Unknown trajectory kind " + std::to_string(kind_code) +
                             ". Valid values are 0 (positions), 1 (velocities), 2 (forces).");
  }
}

// Convert nullable C strings from Python into C++ strings with safe defaults.
std::string toString(const char* text) {
  return (text == nullptr) ? std::string() : std::string(text);
}

// Standard two-call vector export helper: first call with out = null to get required length.
int copyVectorToBuffer(const std::vector<double> &data, double* out, const int out_len) {
  const int need = static_cast<int>(data.size());
  if (out == nullptr) {
    return need;
  }
  if (out_len < need) {
    throw std::runtime_error("Output buffer length (" + std::to_string(out_len) +
                             ") is smaller than required size (" + std::to_string(need) + ").");
  }
  std::copy(data.begin(), data.end(), out);
  return need;
}

} // namespace

extern "C" {

// -----------------------------------------------------------------------------
// Global diagnostics
// -----------------------------------------------------------------------------
const char* stormm_get_last_error() {
  return stormm_last_error.c_str();
}

int stormm_state_variable_count() {
  return static_cast<int>(StateVariable::ALL_STATES);
}

const char* stormm_state_variable_name(const int state_index) {
  static thread_local std::string name;
  return wrapValue([&]() -> const char* {
    const int nstate = static_cast<int>(StateVariable::ALL_STATES);
    if (state_index < 0 || state_index >= nstate) {
      throw std::runtime_error("State variable index " + std::to_string(state_index) +
                               " is out of range [0, " + std::to_string(nstate - 1) + "].");
    }
    name = stormm::energy::getEnumerationName(static_cast<StateVariable>(state_index));
    return name.c_str();
  }, static_cast<const char*>(nullptr));
}

// -----------------------------------------------------------------------------
// AtomGraph wrappers
// -----------------------------------------------------------------------------
void* stormm_atomgraph_create(const char* topology_file) {
  return wrapValue([&]() -> void* {
    if (topology_file == nullptr) {
      throw std::runtime_error("Topology file path cannot be null.");
    }
    AtomGraph* ag = new AtomGraph(std::string(topology_file));
    return reinterpret_cast<void*>(ag);
  }, static_cast<void*>(nullptr));
}

int stormm_atomgraph_destroy(void* atomgraph_handle) {
  return wrapVoid([&]() {
    AtomGraph* ag = requireHandle<AtomGraph>(atomgraph_handle, "AtomGraph");
    delete ag;
  });
}

int stormm_atomgraph_get_atom_count(const void* atomgraph_handle) {
  return wrapValue([&]() -> int {
    const AtomGraph* ag = requireHandle<AtomGraph>(atomgraph_handle, "AtomGraph");
    return ag->getAtomCount();
  }, -1);
}

int stormm_atomgraph_get_residue_count(const void* atomgraph_handle) {
  return wrapValue([&]() -> int {
    const AtomGraph* ag = requireHandle<AtomGraph>(atomgraph_handle, "AtomGraph");
    return ag->getResidueCount();
  }, -1);
}

double stormm_atomgraph_get_total_mass(const void* atomgraph_handle) {
  return wrapValue([&]() -> double {
    const AtomGraph* ag = requireHandle<AtomGraph>(atomgraph_handle, "AtomGraph");
    return ag->getTotalMass();
  }, -1.0);
}

int stormm_atomgraph_get_atomic_masses(const void* atomgraph_handle, double* out_masses,
                                       const int out_len) {
  return wrapValue([&]() -> int {
    const AtomGraph* ag = requireHandle<AtomGraph>(atomgraph_handle, "AtomGraph");
    const std::vector<double> masses = ag->getAtomicMass<double>();
    return copyVectorToBuffer(masses, out_masses, out_len);
  }, -1);
}

// -----------------------------------------------------------------------------
// PhaseSpace wrappers
// -----------------------------------------------------------------------------
void* stormm_phasespace_create(const char* coordinate_file, const void* atomgraph_handle) {
  return wrapValue([&]() -> void* {
    if (coordinate_file == nullptr) {
      throw std::runtime_error("Coordinate file path cannot be null.");
    }
    const AtomGraph* ag = requireHandle<AtomGraph>(atomgraph_handle, "AtomGraph");
    PhaseSpace* ps = new PhaseSpace(std::string(coordinate_file), *ag);
    return reinterpret_cast<void*>(ps);
  }, static_cast<void*>(nullptr));
}

int stormm_phasespace_destroy(void* phasespace_handle) {
  return wrapVoid([&]() {
    PhaseSpace* ps = requireHandle<PhaseSpace>(phasespace_handle, "PhaseSpace");
    delete ps;
  });
}

int stormm_phasespace_get_atom_count(const void* phasespace_handle) {
  return wrapValue([&]() -> int {
    const PhaseSpace* ps = requireHandle<PhaseSpace>(phasespace_handle, "PhaseSpace");
    return ps->getAtomCount();
  }, -1);
}

int stormm_phasespace_get_interlaced_coordinates(const void* phasespace_handle,
                                                 const int trajectory_kind_code, double* out_xyz,
                                                 const int out_len) {
  return wrapValue([&]() -> int {
    const PhaseSpace* ps = requireHandle<PhaseSpace>(phasespace_handle, "PhaseSpace");
    const std::vector<double> xyz = ps->getInterlacedCoordinates(toTrajectoryKind(trajectory_kind_code));
    return copyVectorToBuffer(xyz, out_xyz, out_len);
  }, -1);
}

int stormm_phasespace_export(const void* phasespace_handle, const char* output_file,
                             const double current_time) {
  return wrapVoid([&]() {
    if (output_file == nullptr) {
      throw std::runtime_error("Output file path cannot be null.");
    }
    const PhaseSpace* ps = requireHandle<PhaseSpace>(phasespace_handle, "PhaseSpace");
    ps->exportToFile(std::string(output_file), current_time);
  });
}

// -----------------------------------------------------------------------------
// DynamicsControls wrappers
// -----------------------------------------------------------------------------
void* stormm_dynamics_controls_create() {
  return wrapValue([&]() -> void* {
    DynamicsControls* dyncon = new DynamicsControls();
    return reinterpret_cast<void*>(dyncon);
  }, static_cast<void*>(nullptr));
}

int stormm_dynamics_controls_destroy(void* controls_handle) {
  return wrapVoid([&]() {
    DynamicsControls* dyncon = requireHandle<DynamicsControls>(controls_handle, "DynamicsControls");
    delete dyncon;
  });
}

int stormm_dynamics_controls_get_step_count(const void* controls_handle) {
  return wrapValue([&]() -> int {
    const DynamicsControls* dyncon = requireHandle<DynamicsControls>(controls_handle,
                                                                     "DynamicsControls");
    return dyncon->getStepCount();
  }, -1);
}

int stormm_dynamics_controls_set_step_count(void* controls_handle, const int step_count) {
  return wrapVoid([&]() {
    DynamicsControls* dyncon = requireHandle<DynamicsControls>(controls_handle, "DynamicsControls");
    dyncon->setStepCount(step_count);
  });
}

int stormm_dynamics_controls_set_time_step(void* controls_handle, const double time_step) {
  return wrapVoid([&]() {
    DynamicsControls* dyncon = requireHandle<DynamicsControls>(controls_handle, "DynamicsControls");
    dyncon->setTimeStep(time_step);
  });
}

int stormm_dynamics_controls_set_diagnostic_frequency(void* controls_handle,
                                                      const int diagnostic_frequency) {
  return wrapVoid([&]() {
    DynamicsControls* dyncon = requireHandle<DynamicsControls>(controls_handle, "DynamicsControls");
    dyncon->setDiagnosticPrintFrequency(diagnostic_frequency);
  });
}

int stormm_dynamics_controls_set_trajectory_frequency(void* controls_handle,
                                                      const int trajectory_frequency) {
  return wrapVoid([&]() {
    DynamicsControls* dyncon = requireHandle<DynamicsControls>(controls_handle, "DynamicsControls");
    dyncon->setTrajectoryPrintFrequency(trajectory_frequency);
  });
}

int stormm_dynamics_controls_set_com_purge_frequency(void* controls_handle,
                                                     const int purge_frequency) {
  return wrapVoid([&]() {
    DynamicsControls* dyncon = requireHandle<DynamicsControls>(controls_handle, "DynamicsControls");
    dyncon->setCenterOfMassMotionPurgeFrequency(purge_frequency);
  });
}

int stormm_dynamics_controls_set_thermostat_kind(void* controls_handle,
                                                 const int thermostat_kind_code) {
  return wrapVoid([&]() {
    DynamicsControls* dyncon = requireHandle<DynamicsControls>(controls_handle, "DynamicsControls");
    dyncon->setThermostatKind(toThermostatKind(thermostat_kind_code));
  });
}

int stormm_dynamics_controls_set_temperature(void* controls_handle, const double initial_temperature,
                                             const double final_temperature) {
  return wrapVoid([&]() {
    DynamicsControls* dyncon = requireHandle<DynamicsControls>(controls_handle, "DynamicsControls");
    dyncon->setThermostatGroup(initial_temperature, final_temperature);
  });
}

int stormm_dynamics_controls_set_langevin_frequency(void* controls_handle,
                                                    const double collision_frequency) {
  return wrapVoid([&]() {
    DynamicsControls* dyncon = requireHandle<DynamicsControls>(controls_handle, "DynamicsControls");
    dyncon->setLangevinFrequency(collision_frequency);
  });
}

int stormm_dynamics_controls_set_andersen_frequency(void* controls_handle,
                                                    const int reassignment_frequency) {
  return wrapVoid([&]() {
    DynamicsControls* dyncon = requireHandle<DynamicsControls>(controls_handle, "DynamicsControls");
    dyncon->setAndersenFrequency(reassignment_frequency);
  });
}

int stormm_dynamics_controls_set_seed(void* controls_handle, const int random_seed) {
  return wrapVoid([&]() {
    DynamicsControls* dyncon = requireHandle<DynamicsControls>(controls_handle, "DynamicsControls");
    dyncon->setThermostatSeed(random_seed);
  });
}

// -----------------------------------------------------------------------------
// One-system dynamics execution wrapper
// -----------------------------------------------------------------------------
int stormm_dynamics_run(void* phasespace_handle, const void* atomgraph_handle,
                        const void* controls_handle, const int thermostat_kind_code,
                        const double thermostat_temperature, const char* trajectory_file,
                        const char* restart_file, double* potential_energy_out,
                        double* total_energy_out, int* sample_count_out,
                        double* instantaneous_states_out, const int state_buffer_len) {
  return wrapValue([&]() -> int {
    PhaseSpace* ps = requireHandle<PhaseSpace>(phasespace_handle, "PhaseSpace");
    const AtomGraph* ag = requireHandle<AtomGraph>(atomgraph_handle, "AtomGraph");

    DynamicsControls dyncon = (controls_handle == nullptr) ? DynamicsControls() :
      *requireHandle<DynamicsControls>(controls_handle, "DynamicsControls");
    if (thermostat_kind_code >= 0) {
      dyncon.setThermostatKind(toThermostatKind(thermostat_kind_code));
    }
    if (thermostat_temperature > 0.0) {
      dyncon.setThermostatGroup(thermostat_temperature, thermostat_temperature);
    }

    const ThermostatKind thermostat_kind = dyncon.getThermostatKind();
    const double target_temperature = (thermostat_temperature > 0.0) ? thermostat_temperature :
                                      stormm::namelist::default_simulation_temperature;
    Thermostat heat_bath(ag, thermostat_kind, target_temperature);
    ScoreCard score(1);
    const NeckGeneralizedBornTable gb_table;
    const StaticExclusionMask semask(ag);
    const RestraintApparatus restraints(ag);
    dynamics(ps, &heat_bath, &score, *ag, gb_table, semask, restraints, dyncon, 0,
             toString(trajectory_file), toString(restart_file));

    if (potential_energy_out != nullptr) {
      *potential_energy_out = score.reportPotentialEnergy(0);
    }
    if (total_energy_out != nullptr) {
      *total_energy_out = score.reportTotalEnergy(0);
    }
    if (sample_count_out != nullptr) {
      *sample_count_out = score.getSampleSize();
    }

    const std::vector<double> inst = score.reportInstantaneousStates(0);
    const int state_count = static_cast<int>(inst.size());
    if (instantaneous_states_out != nullptr) {
      if (state_buffer_len < state_count) {
        throw std::runtime_error("State output buffer length (" +
                                 std::to_string(state_buffer_len) + ") is smaller than required "
                                 "size (" + std::to_string(state_count) + ").");
      }
      std::copy(inst.begin(), inst.end(), instantaneous_states_out);
    }
    return state_count;
  }, -1);
}

// -----------------------------------------------------------------------------
// AtomGraphSynthesis wrappers
// -----------------------------------------------------------------------------
void* stormm_atomgraph_synthesis_create(void* const* atomgraph_handles, const int count) {
  return wrapValue([&]() -> void* {
    if (count <= 0) {
      throw std::runtime_error("AtomGraphSynthesis requires at least one AtomGraph.");
    }
    if (atomgraph_handles == nullptr) {
      throw std::runtime_error("Array of AtomGraph handles cannot be null.");
    }
    std::vector<AtomGraph*> ag_list(count);
    for (int i = 0; i < count; i++) {
      ag_list[i] = requireHandle<AtomGraph>(atomgraph_handles[i], "AtomGraph");
    }
    AtomGraphSynthesis* ags = new AtomGraphSynthesis(ag_list);
    return reinterpret_cast<void*>(ags);
  }, static_cast<void*>(nullptr));
}

int stormm_atomgraph_synthesis_destroy(void* ags_handle) {
  return wrapVoid([&]() {
    AtomGraphSynthesis* ags = requireHandle<AtomGraphSynthesis>(ags_handle, "AtomGraphSynthesis");
    delete ags;
  });
}

int stormm_atomgraph_synthesis_get_system_count(const void* ags_handle) {
  return wrapValue([&]() -> int {
    const AtomGraphSynthesis* ags = requireHandle<AtomGraphSynthesis>(ags_handle,
                                                                      "AtomGraphSynthesis");
    return ags->getSystemCount();
  }, -1);
}

int stormm_atomgraph_synthesis_get_unique_topology_count(const void* ags_handle) {
  return wrapValue([&]() -> int {
    const AtomGraphSynthesis* ags = requireHandle<AtomGraphSynthesis>(ags_handle,
                                                                      "AtomGraphSynthesis");
    return ags->getUniqueTopologyCount();
  }, -1);
}

int stormm_atomgraph_synthesis_get_atom_count(const void* ags_handle) {
  return wrapValue([&]() -> int {
    const AtomGraphSynthesis* ags = requireHandle<AtomGraphSynthesis>(ags_handle,
                                                                      "AtomGraphSynthesis");
    return ags->getAtomCount();
  }, -1);
}

int stormm_atomgraph_synthesis_get_atom_count_system(const void* ags_handle, const int system_index) {
  return wrapValue([&]() -> int {
    const AtomGraphSynthesis* ags = requireHandle<AtomGraphSynthesis>(ags_handle,
                                                                      "AtomGraphSynthesis");
    return ags->getAtomCount(system_index);
  }, -1);
}

// -----------------------------------------------------------------------------
// PhaseSpaceSynthesis wrappers
// -----------------------------------------------------------------------------
void* stormm_phasespace_synthesis_create(void* const* phasespace_handles,
                                         void* const* atomgraph_handles, const int count) {
  return wrapValue([&]() -> void* {
    if (count <= 0) {
      throw std::runtime_error("PhaseSpaceSynthesis requires at least one system.");
    }
    if (phasespace_handles == nullptr || atomgraph_handles == nullptr) {
      throw std::runtime_error("PhaseSpace and AtomGraph handle arrays cannot be null.");
    }

    std::vector<PhaseSpace> ps_list;
    ps_list.reserve(count);
    std::vector<AtomGraph*> ag_list(count);
    for (int i = 0; i < count; i++) {
      const PhaseSpace* ps = requireHandle<PhaseSpace>(phasespace_handles[i], "PhaseSpace");
      ag_list[i] = requireHandle<AtomGraph>(atomgraph_handles[i], "AtomGraph");
      ps_list.push_back(*ps);
    }
    PhaseSpaceSynthesis* pss = new PhaseSpaceSynthesis(ps_list, ag_list);
    return reinterpret_cast<void*>(pss);
  }, static_cast<void*>(nullptr));
}

int stormm_phasespace_synthesis_destroy(void* pss_handle) {
  return wrapVoid([&]() {
    PhaseSpaceSynthesis* pss = requireHandle<PhaseSpaceSynthesis>(pss_handle,
                                                                  "PhaseSpaceSynthesis");
    delete pss;
  });
}

int stormm_phasespace_synthesis_get_system_count(const void* pss_handle) {
  return wrapValue([&]() -> int {
    const PhaseSpaceSynthesis* pss = requireHandle<PhaseSpaceSynthesis>(pss_handle,
                                                                        "PhaseSpaceSynthesis");
    return pss->getSystemCount();
  }, -1);
}

int stormm_phasespace_synthesis_get_unique_topology_count(const void* pss_handle) {
  return wrapValue([&]() -> int {
    const PhaseSpaceSynthesis* pss = requireHandle<PhaseSpaceSynthesis>(pss_handle,
                                                                        "PhaseSpaceSynthesis");
    return pss->getUniqueTopologyCount();
  }, -1);
}

int stormm_phasespace_synthesis_get_padded_atom_count(const void* pss_handle) {
  return wrapValue([&]() -> int {
    const PhaseSpaceSynthesis* pss = requireHandle<PhaseSpaceSynthesis>(pss_handle,
                                                                        "PhaseSpaceSynthesis");
    return pss->getPaddedAtomCount();
  }, -1);
}

int stormm_phasespace_synthesis_get_atom_count(const void* pss_handle, const int system_index) {
  return wrapValue([&]() -> int {
    const PhaseSpaceSynthesis* pss = requireHandle<PhaseSpaceSynthesis>(pss_handle,
                                                                        "PhaseSpaceSynthesis");
    return pss->getAtomCount(system_index);
  }, -1);
}

int stormm_phasespace_synthesis_get_interlaced_coordinates(const void* pss_handle,
                                                           const int system_index,
                                                           const int trajectory_kind_code,
                                                           double* out_xyz, const int out_len) {
  return wrapValue([&]() -> int {
    const PhaseSpaceSynthesis* pss = requireHandle<PhaseSpaceSynthesis>(pss_handle,
                                                                        "PhaseSpaceSynthesis");
    const std::vector<double> xyz = pss->getInterlacedCoordinates(system_index,
                                                                   toTrajectoryKind(trajectory_kind_code));
    return copyVectorToBuffer(xyz, out_xyz, out_len);
  }, -1);
}

} // extern "C"
