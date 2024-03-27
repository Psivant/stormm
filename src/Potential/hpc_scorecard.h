// -*-c++-*-
#ifndef STORMM_HPC_SCORECARD_H
#define STORMM_HPC_SCORECARD_H

#include <vector>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "energy_enumerators.h"

namespace stormm {
namespace energy {

using card::GpuDetails;

/// \brief Launch the appropriate kernel to initialize some or all of a scorecard on the device.
///
/// Overloaded:
///   - Accept a single state variable
///   - Accept a list of state variables
///
/// \param var           One or more state variables to initialize
/// \param system_index  Index of the system to initialize (-1 will initialize all systems)
/// \param accumulators  Pointer to the energy tracking object's instantaneous accumulators on
///                      the HPC device
/// \param system_count  Total number of systems traked by the ScoreCard (energy tracking object)
/// \{
void launchScoreCardInitialization(StateVariable var, int system_index, llint* accumulators,
                                   int system_count, const GpuDetails &gpu);

void launchScoreCardInitialization(const std::vector<StateVariable> &var, int system_index,
                                   llint* accumulators, int system_count, const GpuDetails &gpu);
/// \}

/// \brief Launch the kernel to sum potential and total (potential plus kinetic) energies for all
///        systems and all time points on the GPU.
///
/// \param system_count   Total number of systems being tracked
/// \param sample_count   Current number of samples being stored in the energy time series array
/// \param inst_acc       Pointer to the energy tracking object's instantaneous accumulators on
///                       the HPC device
/// \param run_acc        Pointer to the energy tracking object's running average accumulators on
///                       the HPC device
/// \param sqd_acc        Pointer to the energy tracking object's variance accumulators on the HPC
///                       device
/// \param time_ser_acc   Pointer to the energy tracking object's time series data storage on the
///                       HPC device
/// \param inv_nrg_scale  Inverse energy sclaing factor to convert out of fixed-precision format
/// \param gpu            Details of the GPU to use
void launchScoreCardEnergySum(const int system_count, const int sample_count, llint* inst_acc,
                              llint* time_ser_acc, double* run_acc, double* sqd_acc,
                              const double inv_nrg_scale, const GpuDetails &gpu);

/// \brief Launch the kernel to transfer the requested energies to running first and second moment
///        accumulators, as well as commit them to saved time series data on the HPC device.
///
/// Overloaded:
///   - Accept a single state variable
///   - Accept a list of state variables
///
/// Parameters follow from launchScoreCardEnergySum() above, with the additions of:
///
/// \param var           One or more state variables to initialize
/// \param system_index  Index of the system to initialize (-1 will initialize all systems)
/// \{
void launchScoreCardCommit(const std::vector<StateVariable> &var, const int system_index,
                           const int system_count, const size_t sample_count,
                           const llint* inst_acc, double* run_acc, double* sqd_acc,
                           llint* time_ser_acc, const GpuDetails &gpu);

void launchScoreCardCommit(StateVariable var, const int system_index, const int system_count,
                           const size_t sample_count, const llint* inst_acc, double* run_acc,
                           double* sqd_acc, llint* time_ser_acc, const GpuDetails &gpu);
/// \}

} // namespace energy
} // namespace stormm

#endif
