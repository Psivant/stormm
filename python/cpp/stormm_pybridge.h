// -*-c++-*-
#ifndef STORMM_PYBRIDGE_H
#define STORMM_PYBRIDGE_H

#include "copyright.h"

/// \file
/// \brief C ABI bridge declarations for exposing selected STORMM C++ features to Python.
///
/// This header declares a narrow ABI that Python can load with ctypes.  All object handles are
/// opaque pointers managed by creation / destruction functions below.

extern "C" {

/// \brief Get a thread-local error message from the last failing bridge call.
const char* stormm_get_last_error();

/// \brief Get the number of entries in the StateVariable enum.
int stormm_state_variable_count();

/// \brief Get the canonical name of one StateVariable enum entry.
///
/// \param state_index  Zero-based index into StateVariable values
const char* stormm_state_variable_name(int state_index);

/// \brief Construct and destroy an AtomGraph from an input topology.
///
/// \param topology_file  Path to an Amber topology (prmtop-compatible) file
/// \{
void* stormm_atomgraph_create(const char* topology_file);
int stormm_atomgraph_destroy(void* atomgraph_handle);
/// \}

/// \brief Query basic AtomGraph metadata and per-atom masses.
///
/// For vector data, pass out pointer null to get the required length, then call again with a
/// sufficiently large buffer.
///
/// \{
int stormm_atomgraph_get_atom_count(const void* atomgraph_handle);
int stormm_atomgraph_get_residue_count(const void* atomgraph_handle);
double stormm_atomgraph_get_total_mass(const void* atomgraph_handle);
int stormm_atomgraph_get_atomic_masses(const void* atomgraph_handle, double* out_masses,
                                       int out_len);
/// \}

/// \brief Construct and destroy a PhaseSpace from coordinates plus an AtomGraph.
///
/// \param coordinate_file   Path to coordinate / restart input
/// \param atomgraph_handle  Opaque AtomGraph pointer
/// \{
void* stormm_phasespace_create(const char* coordinate_file, const void* atomgraph_handle);
int stormm_phasespace_destroy(void* phasespace_handle);
/// \}

/// \brief Query basic PhaseSpace properties and interlaced coordinates.
///
/// \param trajectory_kind_code  0 positions, 1 velocities, 2 forces
/// \{
int stormm_phasespace_get_atom_count(const void* phasespace_handle);
int stormm_phasespace_get_interlaced_coordinates(const void* phasespace_handle,
                                                 int trajectory_kind_code, double* out_xyz,
                                                 int out_len);
int stormm_phasespace_export(const void* phasespace_handle, const char* output_file,
                             double current_time);
/// \}

/// \brief Construct and destroy DynamicsControls objects.
/// \{
void* stormm_dynamics_controls_create();
int stormm_dynamics_controls_destroy(void* controls_handle);
/// \}

/// \brief Query and mutate dynamics control parameters from Python.
/// \{
int stormm_dynamics_controls_get_step_count(const void* controls_handle);
int stormm_dynamics_controls_set_step_count(void* controls_handle, int step_count);
int stormm_dynamics_controls_set_time_step(void* controls_handle, double time_step);
int stormm_dynamics_controls_set_diagnostic_frequency(void* controls_handle,
                                                      int diagnostic_frequency);
int stormm_dynamics_controls_set_trajectory_frequency(void* controls_handle,
                                                      int trajectory_frequency);
int stormm_dynamics_controls_set_com_purge_frequency(void* controls_handle, int purge_frequency);
int stormm_dynamics_controls_set_thermostat_kind(void* controls_handle, int thermostat_kind_code);
int stormm_dynamics_controls_set_temperature(void* controls_handle, double initial_temperature,
                                             double final_temperature);
int stormm_dynamics_controls_set_langevin_frequency(void* controls_handle,
                                                    double collision_frequency);
int stormm_dynamics_controls_set_andersen_frequency(void* controls_handle,
                                                    int reassignment_frequency);
int stormm_dynamics_controls_set_seed(void* controls_handle, int random_seed);
/// \}

/// \brief Run CPU dynamics on one PhaseSpace + AtomGraph system.
///
/// Pass `controls_handle` as null to use default DynamicsControls.  Pass output pointers as null
/// for values that are not needed.
///
/// \param thermostat_kind_code  -1 keep controls setting, otherwise 0..3 override thermostat type
/// \param thermostat_temperature  If > 0, applies a flat thermostat target in Kelvin
/// \return Number of instantaneous state values written (or required if buffer is null)
int stormm_dynamics_run(void* phasespace_handle, const void* atomgraph_handle,
                        const void* controls_handle, int thermostat_kind_code,
                        double thermostat_temperature, const char* trajectory_file,
                        const char* restart_file, double* potential_energy_out,
                        double* total_energy_out, int* sample_count_out,
                        double* instantaneous_states_out, int state_buffer_len);

/// \brief Construct and destroy an AtomGraphSynthesis from multiple AtomGraph handles.
/// \{
void* stormm_atomgraph_synthesis_create(void* const* atomgraph_handles, int count);
int stormm_atomgraph_synthesis_destroy(void* ags_handle);
/// \}

/// \brief Query AtomGraphSynthesis metadata.
/// \{
int stormm_atomgraph_synthesis_get_system_count(const void* ags_handle);
int stormm_atomgraph_synthesis_get_unique_topology_count(const void* ags_handle);
int stormm_atomgraph_synthesis_get_atom_count(const void* ags_handle);
int stormm_atomgraph_synthesis_get_atom_count_system(const void* ags_handle, int system_index);
/// \}

/// \brief Construct and destroy a PhaseSpaceSynthesis from PhaseSpace + AtomGraph handles.
/// \{
void* stormm_phasespace_synthesis_create(void* const* phasespace_handles,
                                         void* const* atomgraph_handles, int count);
int stormm_phasespace_synthesis_destroy(void* pss_handle);
/// \}

/// \brief Query PhaseSpaceSynthesis metadata and interlaced coordinates by system index.
/// \{
int stormm_phasespace_synthesis_get_system_count(const void* pss_handle);
int stormm_phasespace_synthesis_get_unique_topology_count(const void* pss_handle);
int stormm_phasespace_synthesis_get_padded_atom_count(const void* pss_handle);
int stormm_phasespace_synthesis_get_atom_count(const void* pss_handle, int system_index);
int stormm_phasespace_synthesis_get_interlaced_coordinates(const void* pss_handle,
                                                           int system_index,
                                                           int trajectory_kind_code,
                                                           double* out_xyz, int out_len);
/// \}

} // extern "C"

#endif
