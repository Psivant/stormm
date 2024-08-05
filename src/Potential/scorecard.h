// -*-c++-*-
#ifndef STORMM_SCORECARD_H
#define STORMM_SCORECARD_H

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/stormm_vector_types.h"
#include "energy_enumerators.h"

namespace stormm {
namespace energy {

using card::GpuDetails;
using card::Hybrid;
using card::HybridTargetLevel;
using numerics::default_energy_scale_bits;

/// \brief Read-only abstract for the ScoreCard object.  This is needed more for completeness than
///        anything else, but could be useful in cases where virials are to be translated into
///        rescaling coefficients on the GPU and then immediately used to rescale positions in a
///        Berendsen barostatting situation.
struct ScoreCardReader {

  /// \brief The constructor is, as usual, a collection of the relevant constants and pointers.
  ScoreCardReader(int system_count_in, int data_stride_in, int sampled_step_count_in,
                  float nrg_scale_f_in, double nrg_scale_lf_in, float inverse_nrg_scale_f_in,
                  double inverse_nrg_scale_lf_in, const int* time_steps_in,
                  const llint* instantaneous_accumulators_in,
                  const double* running_accumulators_in, const double* squared_accumulators_in,
                  const llint* time_series_in);

  /// \brief Take only the default copy constructor and copy assignemnt operator for this
  ///        struct with const elements.
  /// \{
  ScoreCardReader(const ScoreCardReader &original) = default;
  ScoreCardReader(ScoreCardReader &&original) = default;
  /// \}
  
  const int system_count;                  ///< Number of independent systems tracked
  const int data_stride;                   ///< Size of the StateVariable enumerator rounded up
                                           ///<   to the nearest multiple of the HPC warp size
  const int sampled_step_count;            ///< The number of steps in the sample (a nice result
                                           ///<   is that updating this member variable in the
                                           ///<   parent ScoreCard right before making the
                                           ///<   abstract can send this data to the GPU)
  const float nrg_scale_f;                 ///< Conversion factor for fixed-precision accumulation
  const double nrg_scale_lf;               ///< Conversion factor for fixed-precision accumulation
  const float inverse_nrg_scale_f;         ///< Conversion for fixed-precision interpretation
  const double inverse_nrg_scale_lf;       ///< Conversion for fixed-precision interpretation
  const int* time_steps;                   ///< Time steps in the simulation at which each energy
                                           ///<   sample was taken
  const llint* instantaneous_accumulators; ///< State variables for each system
  const double* running_accumulators;      ///< Running sums of state variables for each system
  const double* squared_accumulators;      ///< Running squared sums of state variables for each
                                           ///<   system
  const llint* time_series;                ///< Details of energy values in each component's
                                           ///<   accumulator at each time step
};

/// \brief Writeable abstract for the ScoreCard object, useful for accumulating energies in many
///        kernels.
struct ScoreCardWriter {

  /// \brief The constructor is, as usual, a collection of the relevant constants and pointers.
  ScoreCardWriter(int system_count_in, int data_stride_in, int sampled_step_count_in,
                  float nrg_scale_f_in, double nrg_scale_lf_in, float inverse_nrg_scale_f_in,
                  double inverse_nrg_scale_lf_in, int* time_steps_in,
                  llint* instantaneous_accumulators_in, double* running_accumulators_in,
                  double* squared_accumulators_in, llint* time_series_in);

  /// \brief Take only the default copy constructor and copy assignemnt operator for this
  ///        struct with const elements.
  /// \{
  ScoreCardWriter(const ScoreCardWriter &original) = default;
  ScoreCardWriter(ScoreCardWriter &&original) = default;
  /// \}

  const int system_count;             ///< Number of independent systems tracked
  const int data_stride;              ///< Size of the StateVariable enumerator rounded up to the
                                      ///<   nearest multiple of the HPC warp size
  const int sampled_step_count;       ///< The number of steps in the sample
  const float nrg_scale_f;            ///< Conversion factor for fixed-precision accumulation
  const double nrg_scale_lf;          ///< Conversion factor for fixed-precision accumulation
  const float inverse_nrg_scale_f;    ///< Conversion for fixed-precision interpretation
  const double inverse_nrg_scale_lf;  ///< Conversion for fixed-precision interpretation
  int* time_steps;                    ///< Time steps in the simulation at which each energy
                                      ///<   sample was taken
  llint* instantaneous_accumulators;  ///< State variables for each system
  double* running_accumulators;       ///< Running sums of state variables for each system
  double* squared_accumulators;       ///< Running squared sums of state variables for each system
  llint* time_series;                 ///< Details of energy values in each component's
                                      ///<   accumulator at each time step
};

/// \brief Track the energy components of a collection of systems in an HPC-capable array.  This
///        object uses the familiar trick of defining an enumerator (StateVariables) with a final
///        entry to indicate its total length, so that if more energy components need to be tracked
///        in the future the storage and indexing can automatically adjust with new entries.
class ScoreCard {
public:

  /// \brief The constructor requires only the number of systems.
  ///
  /// \param system_count_in    The number of systems to track
  /// \param capacity_in        The capacity (samples per energy component in each system) to
  ///                           initially allocate for
  /// \param nrg_scale_bits_in  Number of bits after the decimal with which to store energy values
  ScoreCard(int system_count_in, int capacity_in = 16,
            int nrg_scale_bits_in = default_energy_scale_bits);

  /// \brief Copy and move constructors help make the energy tracker a first-class C++ object.
  ///        The defaults apply as there are no const member variables and no POINTER-kind Hybrid
  ///        objects.
  ///
  /// \param original  The original energy tracking object to build from
  /// \{
  ScoreCard(const ScoreCard &original) = default;
  ScoreCard(ScoreCard &&original) = default;
  /// \}

  /// \brief Copy and move assignment operators complete the support for integrating this energy
  ///        tracker with Standard Template Library mechanics.  The defaults apply as there are no
  ///        const member variables and no POINTER-kind Hybrid objects.
  ///
  /// \param other  Object to which the present ScoreCard shall be equated
  /// \{
  ScoreCard& operator=(const ScoreCard &other) = default;
  ScoreCard& operator=(ScoreCard &&other) = default;
  /// \}
  
  /// \brief Get the number of systems that this object is tracking
  int getSystemCount() const;

  /// \brief Get the number of steps that have been sampled
  int getSampleSize() const;

  /// \brief Get the stride used to store separate results for each system in the instantaneous
  ///        accumulator arrays.
  int getDataStride() const;

  /// \brief Get the sample capacity (the number of individual energy measurements that the object
  ///        is prepared to store).
  int getSampleCapacity() const;
  
  /// \brief Get the number of bits of fixed precision to which results are stored
  int getEnergyScaleBits() const;
  
  /// \brief Get the energy scaling factors in single- or double-precision floating point format
  /// \{
  template <typename T> T getEnergyScalingFactor() const;
  template <typename T> T getInverseEnergyScalingFactor() const;
  /// \}

  /// \brief Get the time step at which a particular energy value was taken.
  ///
  /// \param time_index  The index of the energies and time step number in question
  int getTimeStep(int time_index) const;

#ifdef STORMM_USE_HPC
  /// \brief Upload data to the device (this could be useful in situations where the CPU is
  ///        required to compute additional energy terms, or other energy quantities are being
  ///        brought in from some module outside the typical STORMM molecular mechanics routines).
  void upload();

  /// \brief Download all data from the device (each of the report(...) functions below will call
  ///        the download method for the appropriate array in an HPC setting, as the data on the
  ///        device is assumed to be the most important).
  void download();
#endif

  /// \brief Get the appropriate abstract based on the const-ness of the abstract
  ///
  /// \param tier  Get pointers to data on the host or on the HPC accelerator device
  /// \{
  const ScoreCardReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  ScoreCardWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Reserve space for storing sampled energy component values.
  ///
  /// \param new_capacity  The new capacity to allocate for
  void reserve(int new_capacity);
  
  /// \brief Place a result into one of the instantaneous state variable accumulators.  This is
  ///        for CPU activity; on the GPU, the contributions will occur as part of each energy
  ///        kernel using pointers.  This will automatically update running accumulators for
  ///        statistical tracking.  It is add + commit, below.
  ///
  /// \param var           The state variable to which this contribution belongs, i.e. bond energy
  /// \param amount        Amount to contribute to said state variable (in fixed precision format)
  /// \param system_index  Index of the system (among a list of those being tracked) that the
  ///                      contrbution describes
  void contribute(StateVariable var, llint amount, int system_index = 0);

  /// \brief Initialize some or all instantaneous state variable accumulators.  If all state
  ///        variable accumulators are initialized in all systems, the sample counter will also be
  ///        reset.
  ///
  /// Overloaded:
  ///   - Initialize a single state variable accumulator
  ///   - Initialize a list of state variable accumulators
  ///   - Initialize all state variable accumulators in one or more systems
  ///
  /// \param var           The state variable to initialize
  /// \param system_index  Index of the system (among a list of those being tracked) that the
  ///                      contrbution describes
  /// \param tier          Indicate whether to initialize accumulators on the host or device
  /// \param gpu           Details of the HPC device in use
  /// \{
  void initialize(StateVariable var, int system_index = 0,
                  HybridTargetLevel tier = HybridTargetLevel::HOST,
                  const GpuDetails &gpu = null_gpu);
  void initialize(const std::vector<StateVariable> &var, int system_index = 0,
                  HybridTargetLevel tier = HybridTargetLevel::HOST,
                  const GpuDetails &gpu = null_gpu);
  void initialize(int system_index, HybridTargetLevel tier = HybridTargetLevel::HOST,
                  const GpuDetails &gpu = null_gpu);
  void initialize(HybridTargetLevel tier = HybridTargetLevel::HOST,
                  const GpuDetails &gpu = null_gpu);
  /// \}
  
  /// \brief Add a result to a growing total in one of the instantaneous state variable
  ///        accumulators.  This is for CPU activity; on the GPU, the contributions will occur as
  ///        part of each energy kernel using pointers.
  ///
  /// \param var           The state variable to which this contribution belongs, i.e. bond energy
  /// \param amount        Amount to contribute to said state variable (in fixed precision format)
  /// \param system_index  Index of the system (among a list of those being tracked) that the
  ///                      contrbution describes
  void add(StateVariable var, llint amount, int system_index = 0);  

  /// \brief Commit the (assumed complete) accumulated results in one or more state variables to
  ///        the running tallies and time series kept for statistical purposes.
  ///
  /// Overloaded:
  ///   - Commit results for a single state variable
  ///   - Commit results for a list of specific state variables
  ///   - Commit results for all state variables
  ///   - Commit results for a single system, or for all systems
  ///
  /// \param var           The state variable to which this contribution belongs, i.e. bond energy
  /// \param system_index  Index of the system (among a list of those being tracked) that the
  ///                      contrbution describes
  /// \{
  void commit(StateVariable var, int system_index,
              HybridTargetLevel tier = HybridTargetLevel::HOST, const GpuDetails &gpu = null_gpu);

  void commit(const std::vector<StateVariable> &var, int system_index,
              HybridTargetLevel tier = HybridTargetLevel::HOST, const GpuDetails &gpu = null_gpu);

  void commit(StateVariable var,
              HybridTargetLevel tier = HybridTargetLevel::HOST, const GpuDetails &gpu = null_gpu);

  void commit(const std::vector<StateVariable> &var,
              HybridTargetLevel tier = HybridTargetLevel::HOST, const GpuDetails &gpu = null_gpu);

  void commit(int system_index,
              HybridTargetLevel tier = HybridTargetLevel::HOST, const GpuDetails &gpu = null_gpu);
  
  void commit(HybridTargetLevel tier = HybridTargetLevel::HOST, const GpuDetails &gpu = null_gpu);
  /// \}

  /// \brief Compute total potential energies for all systems in the instantaneous accumulators as
  ///        well as all time series accumulators.  This total is not automatically computed by
  ///        various interaction evaluations, which only compute individual "components" of this
  ///        quantity.
  ///
  /// \param tier  Perform this operation on the CPU host or GPU device
  void computePotentialEnergy(HybridTargetLevel tier = HybridTargetLevel::HOST,
                              const GpuDetails &gpu = null_gpu);

  /// \brief Compute total energies for all systems in the instantaneous accumulators as well as
  ///        all time series accumulators.  This total is not automatically computed by various
  ///        interaction evaluations, which only compute individual "components" of this quantity.
  ///
  /// \param tier  Perform this operation on the CPU host or GPU device
  void computeTotalEnergy(HybridTargetLevel tier = HybridTargetLevel::HOST,
                          const GpuDetails &gpu = null_gpu);

  /// \brief Increment the number of sampled steps.  This will automatically allocate additional
  ///        capacity if the sampled step count reaches the object's capacity.
  void incrementSampleCount();

  /// \brief Reset the sample counter (this is implicitly done by initialize() if that function is
  ///        called to operate on all start variables and all systems, on either the host or HPC
  ///        device).
  ///
  /// \param count_in  The number of samples to set the object as having (default 0, full reset)
  void resetSampleCount(int count_in = 0);
  
  /// \brief Report the total energy, or total potential energy, for one system or for all systems.
  ///        Each result will be summed in the internal fixed-point accumulation before conversion
  ///        to real values in units of kcal/mol.
  ///
  /// Overloaded:
  ///   - Sum the total or potential energies of a selected system
  ///   - Sum the total or potential energies of all systems
  ///
  /// \param system_index  Index of the system of interest
  /// \param tier          Level from which to extract the data
  /// \{
  std::vector<double>
  reportPotentialEnergies(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  double reportPotentialEnergy(int system_index = 0,
                               HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  std::vector<double> reportTotalEnergies(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  double reportTotalEnergy(int system_index = 0,
                           HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Report instantaneous results in kcal/mol, as a double-precision vector.
  ///
  /// Overloaded:
  ///   - Report results for all systems (the vector will be concatenated, with padding removed)
  ///   - Report results for a single system
  ///   - Report a specific result for all systems
  ///   - Report a specific result for a single system
  ///
  /// \param system_index  Index of the system of interest within all of those being tracked
  /// \param aspect        The type of energy or virial quantity of interest
  /// \param tier          Level from which to extract the data
  /// \{
  std::vector<double>
  reportInstantaneousStates(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  std::vector<double>
  reportInstantaneousStates(int system_index,
                            HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  std::vector<double>
  reportInstantaneousStates(StateVariable aspect,
                            HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  double reportInstantaneousStates(StateVariable aspect, int system_index,
                                   HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Report averaged results in kcal/mol, as a double-precision vector.
  ///
  /// Overloaded:
  ///   - Report results for all systems (the vector will be concatenated, with padding removed)
  ///   - Report results for a single system
  ///   - Report a specific result for all systems
  ///   - Report a specific result for a single system
  ///
  /// \param system_index  Index of the system of interest within all of those being tracked
  /// \param aspect        The type of energy or virial quantity of interest
  /// \param tier          Level from which to extract the data
  /// \{
  std::vector<double> reportAverageStates(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  std::vector<double> reportAverageStates(int system_index,
                                          HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  std::vector<double> reportAverageStates(StateVariable aspect,
                                          HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  double reportAverageStates(StateVariable aspect, int system_index,
                             HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Report standard deviations in kcal/mol, as a double-precision vector.
  ///
  /// Overloaded:
  ///   - Report results for all systems (the vector will be concatenated, with padding removed)
  ///   - Report results for a single system
  ///   - Report a specific result for all systems
  ///   - Report a specific result for a single system
  ///
  /// \param system_index  Index of the system of interest within all of those being tracked
  /// \param aspect        The type of energy or virial quantity of interest
  /// \param tier          Level from which to extract the data
  /// \{
  std::vector<double>
  reportVarianceOfStates(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  std::vector<double>
  reportVarianceOfStates(int system_index, HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  std::vector<double>
  reportVarianceOfStates(StateVariable aspect,
                         HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  double reportVarianceOfStates(StateVariable aspect, int system_index,
                                HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Report the energy history for one or more systems.
  ///
  /// Overloaded:
  ///   - Report the total energy history, as recorded on either the host or GPU device
  ///   - Report the history of a specific energy component
  ///   - Report results for a specific system or as the mean and standard deviation (in the x and
  ///     y components of a double2 tuple) of all systems
  ///
  /// \param system_index  Index of the system of interest within all of those being tracked
  /// \param aspect        The type of energy or virial quantity of interest (if none is specified
  ///                      the history of systems' total energy will be reported, taking elements
  ///                      from all potential and kinetic energy sources)
  /// \param tier          Level from which to extract the data
  /// \{
  std::vector<double> reportHistory(int system_index,
                                    HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  std::vector<double> reportHistory(StateVariable aspect, int system_index,
                                    HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  std::vector<double2> reportHistory(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  std::vector<double2> reportHistory(StateVariable aspect,
                                     HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}

  /// \brief Get a const pointer to this object.
  const ScoreCard* getSelfPointer() const;

  /// \brief Set the time index for the most recent stored set of energies.
  ///
  /// \param time_index  Set the time index
  /// \param tier        The memory level at which to set the time step
  void setLastTimeStep(int time_index, HybridTargetLevel tier = HybridTargetLevel::HOST);
  
  /// \brief Import the results of another ScoreCard into this one, including all components and
  ///        the associated energy history.  If the current object does not have sufficient space
  ///        either in terms of systems or sample capacity (depth of history), it will be
  ///        re-allocated.  This functionality is available on the GPU as a free function
  ///        energyCopy(), which allows the reader from the source and the writer from the
  ///        destination to be used rather than a long series of pointers and array size constants.
  ///
  /// Overloaded:
  ///   - Supply a const pointer or const reference to the other energy tracking object
  ///   - Import a single system into a specific index.
  ///   - Import multiple systems into multiple indices.
  ///
  /// \param other           The other energy tracking object
  /// \param fill_index      System index of the current object to import results into
  /// \param fill_indices    System indices of the current object to import results into
  /// \param source_index    System index of the original object to import results from
  /// \param source_indices  System indices of the original object to import results from
  /// \{
  void import(const ScoreCard *other, size_t fill_index, size_t source_index);
  void import(const ScoreCard &other, size_t fill_index, size_t source_index);
  void import(const ScoreCard *other, const std::vector<int> &fill_indices,
              const std::vector<int> &source_indices);
  void import(const ScoreCard &other, const std::vector<int> &fill_indices,
              const std::vector<int> &source_indices);
  /// \}
  
private:
  int system_count;                          ///< The number of systems in the collection (each
                                             ///<   system will get a separate set of accumulators
                                             ///<   and averages, although combined statistics can
                                             ///<   be collected).
  int data_stride;                           ///< The number of accumulators (i.e. bond, angle, ...
                                             ///<   dU/dLambda) per system.  This indicates the
                                             ///<   length of each system's subset of the data in
                                             ///<   the following arrays.
  int sampled_step_count;                    ///< The number of samples going into running averages
                                             ///<   and standard deviations
  int sample_capacity;                       ///< The maximum number of samples available to this
                                             ///<   ScoreCard (the time_series_accumulators array
                                             ///<   will be resized if the sampled_step_count
                                             ///<   reaches this capacity).
  
  /// Scaling factors for fixed-precision accumulation and interpretation.  The defaults from
  /// Constants/fixed_precision.h are good suggestions, but the precision that the code works in
  /// may offer a tradeoff between speed and stability.  For energy accumulation, the question is
  /// not what happens to the future of the dynamics, merely what values are reported and how
  /// accurate they might be.  However, in the interest of knowing quantities to very high
  /// precision these values can be altered at the user or developer's discretion.
  /// \{
  int nrg_scale_bits;
  float nrg_scale_f;
  double nrg_scale_lf;
  float inverse_nrg_scale_f;
  double inverse_nrg_scale_lf;
  /// \}

  Hybrid<int> time_steps;                    ///< The numbers of the time steps in the simulation
                                             ///<   at which each energy sample was taken.
                                             ///<   Multiply each step number by the time step
                                             ///<   increment (available from a DynamicsControls
                                             ///<   object, in units of femtoseconds) to determine
                                             ///<   the time at which each energy sample was taken.
                                             ///<   In other contexts, these values may be step
                                             ///<   numbers in a line minimization calculation.
  Hybrid<llint> instantaneous_accumulators;  ///< Instantaneous accumulators for reporting the
                                             ///<   energy at one time step, accumulated in fixed
                                             ///<   precision after scaling by
                                             ///<   default_energy_scale_lf (see fixed_precision.h
                                             ///<   in Constants/ ).
  Hybrid<double> running_accumulators;       ///< Running sums of the energy and other state
                                             ///<   variables, collected over each sampled time
                                             ///<   step.  These are stored in double precision to
                                             ///<   guard against overflow in a long simulation.
                                             ///<   In an HPC setting, these will be computed on
                                             ///<   the device after preparing the instantaneous
                                             ///<   accumulators.
  Hybrid<double> squared_accumulators;       ///< Running sums of the squared energy components and
                                             ///<   other stat variables collected over each
                                             ///<   sampled time step.  In an HPC setting, these
                                             ///<   will be computed on the device after preparing
                                             ///<   the instantaneous accumulators.
  Hybrid<llint> time_series_accumulators;    ///< A history of values for each energy accumulator
                                             ///<   at each sampled time step.  This detailed
                                             ///<   array is resized as needed, or can be reserved.

  /// \brief Sum the potential energy contributions for a single system arranged in some array of
  ///        long long integers.
  ///
  /// Overloaded:
  ///   - Sum the energy as a long long integer (native fixed-precision format)
  ///   - Sum the energy and return a double-precision real number, scaled to units of kcal/mol
  ///
  /// \param nrg_data  Energies for the system.  The first element of the array will coorespond to
  ///                  the first element of the StateVariable enumerator, i.e. BOND.
  /// \{
  llint sumPotentialEnergyAsLlint(const llint* nrg_data) const;
  double sumPotentialEnergy(const llint* nrg_data) const;
  /// \}

  /// \brief Sum the total energy contributions for a single sustem arranged in some array of
  ///        long long integers.  This adds kinetic energy to the usual potential contributions.
  ///
  /// Overloaded:
  ///   - Sum the energy as a long long integer (native fixed-precision format)
  ///   - Sum the energy and return a double-precision real number, scaled to units of kcal/mol
  ///
  /// \param nrg_data  Energies for the system.  The first element of the array will coorespond to
  ///                  the first element of the StateVariable enumerator, i.e. BOND.
  /// \{
  llint sumTotalEnergyAsLlint(const llint* nrg_data) const;
  double sumTotalEnergy(const llint* nrg_data) const;
  /// \}
};

/// \brief Add a quantity to a ScoreCard using its abstract.  This will have the same effect as the
///        eponymous member function.
///
/// \param scw           The abstract of the energy tracking object
/// \param var           The energy quantity into which the sum goes
/// \param amount        The amount of energy to contribute, scaled to some fixed-precision model
/// \param system_index  Index of the system to which the energy pertains
void add(ScoreCardWriter *csw, StateVariable var, llint amount, int system_index);

} // namespace energy
} // namespace stormm

#include "scorecard.tpp"

#endif
