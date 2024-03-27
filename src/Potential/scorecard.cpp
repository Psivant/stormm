#include <cmath>
#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "DataTypes/common_types.h"
#include "Math/rounding.h"
#include "Math/vector_ops.h"
#include "Math/statistics.h"
#include "scorecard.h"
#ifdef STORMM_USE_HPC
#include "hpc_scorecard.h"
#endif

namespace stormm {
namespace energy {

using stmath::roundUp;
using stmath::mean;
using stmath::variance;
using stmath::VarianceMethod;
using numerics::default_energy_scale_f;
using numerics::default_energy_scale_lf;
using numerics::default_inverse_energy_scale_f;
using numerics::default_inverse_energy_scale_lf;

//-------------------------------------------------------------------------------------------------
ScoreCardReader::ScoreCardReader(const int system_count_in, const int data_stride_in,
                                 const int sampled_step_count_in, const float nrg_scale_f_in,
                                 const double nrg_scale_lf_in, const float inverse_nrg_scale_f_in,
                                 const double inverse_nrg_scale_lf_in, const int* time_steps_in,
                                 const llint* instantaneous_accumulators_in,
                                 const double* running_accumulators_in,
                                 const double* squared_accumulators_in,
                                 const llint* time_series_in) :
    system_count{system_count_in}, data_stride{data_stride_in},
    sampled_step_count{sampled_step_count_in}, nrg_scale_f{nrg_scale_f_in},
    nrg_scale_lf{nrg_scale_lf_in}, inverse_nrg_scale_f{inverse_nrg_scale_f_in},
    inverse_nrg_scale_lf{inverse_nrg_scale_lf_in}, time_steps{time_steps_in},
    instantaneous_accumulators{instantaneous_accumulators_in},
    running_accumulators{running_accumulators_in},
    squared_accumulators{squared_accumulators_in},
    time_series{time_series_in}
{}

//-------------------------------------------------------------------------------------------------
ScoreCardWriter::ScoreCardWriter(const int system_count_in, const int data_stride_in,
                                 const int sampled_step_count_in, const float nrg_scale_f_in,
                                 const double nrg_scale_lf_in, const float inverse_nrg_scale_f_in,
                                 const double inverse_nrg_scale_lf_in, int* time_steps_in,
                                 llint* instantaneous_accumulators_in,
                                 double* running_accumulators_in, double* squared_accumulators_in,
                                 llint* time_series_in) :
    system_count{system_count_in}, data_stride{data_stride_in},
    sampled_step_count{sampled_step_count_in}, nrg_scale_f{nrg_scale_f_in},
    nrg_scale_lf{nrg_scale_lf_in}, inverse_nrg_scale_f{inverse_nrg_scale_f_in},
    inverse_nrg_scale_lf{inverse_nrg_scale_lf_in}, time_steps{time_steps_in},
    instantaneous_accumulators{instantaneous_accumulators_in},
    running_accumulators{running_accumulators_in},
    squared_accumulators{squared_accumulators_in},
    time_series{time_series_in}
{}

//-------------------------------------------------------------------------------------------------
ScoreCard::ScoreCard(const int system_count_in, const int capacity_in,
                     const int nrg_scale_bits_in) :
    system_count{system_count_in},
    data_stride{roundUp(static_cast<int>(StateVariable::ALL_STATES), warp_size_int)},
    sampled_step_count{0},
    sample_capacity{capacity_in},
    nrg_scale_bits{nrg_scale_bits_in},
    nrg_scale_f{0.0},
    nrg_scale_lf{0.0},
    inverse_nrg_scale_f{0.0},
    inverse_nrg_scale_lf{0.0},
    time_steps{static_cast<size_t>(capacity_in), "scorecard_time_steps"},
    instantaneous_accumulators{static_cast<size_t>(data_stride * system_count), "scorecard_acc"},
    running_accumulators{static_cast<size_t>(data_stride * system_count), "score_running_acc"},
    squared_accumulators{static_cast<size_t>(data_stride * system_count), "score_squared_acc"},
    time_series_accumulators{static_cast<size_t>(data_stride * system_count * sample_capacity),
                             "score_card_ts"}
{
  // Check the fixed precision bit count: there are limits that the user or developer will be
  // allowed, for the sake of energy estimates with some degree of accuracy and remaining within
  // the bounds of the long long integer accumulator format.
  if (nrg_scale_bits < 11) {
    rtErr("Energy and virial accumulation must take place in a precision of at least one part in "
          "2048 of one kcal/mol.  A precision of " + std::to_string(nrg_scale_bits) +
          " bits after the decimal is unacceptably low.", "ScoreCard");
  }
  else if (nrg_scale_bits > 40) {
    rtErr("Energy and virial accumulation can overflow the 64-bit integer accumulators.  Storing "
          "the numbers to a precision of more than one part in 1.0995 trillionths of a kcal/mol, "
          "a precision of " + std::to_string(nrg_scale_bits) + " after the decimal, is "
          "unnecessary and risks producing overflow in large systems.", "ScoreCard");
  }

  // Determine the energy scaling factors
  if (nrg_scale_bits == default_energy_scale_bits) {
    nrg_scale_lf = default_energy_scale_lf;
    nrg_scale_f = default_energy_scale_f;
    inverse_nrg_scale_lf = default_inverse_energy_scale_lf;
    inverse_nrg_scale_f = default_inverse_energy_scale_f;
  }
  else {
    int ib = 0;
    nrg_scale_lf = 1.0;
    while (ib + 10 <= nrg_scale_bits) {
      nrg_scale_lf *= 1024.0;
      ib += 10;
    }
    while (ib < nrg_scale_bits) {
      nrg_scale_lf *= 2.0;
      ib++;
    }
    nrg_scale_f = (float)nrg_scale_lf;
    inverse_nrg_scale_lf = 1.0 / nrg_scale_lf;
    inverse_nrg_scale_f = (float)1.0 / nrg_scale_f;
  }
}

//-------------------------------------------------------------------------------------------------
llint ScoreCard::sumPotentialEnergyAsLlint(const llint* nrg_data) const {
  llint lacc = nrg_data[static_cast<size_t>(StateVariable::BOND)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::ANGLE)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::PROPER_DIHEDRAL)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::IMPROPER_DIHEDRAL)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::UREY_BRADLEY)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::CHARMM_IMPROPER)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::CMAP)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::VDW)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::VDW_ONE_FOUR)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::ELECTROSTATIC)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::ELEC_ONE_FOUR)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::GENERALIZED_BORN)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::RESTRAINT)];
  return lacc;
}

//-------------------------------------------------------------------------------------------------
double ScoreCard::sumPotentialEnergy(const llint* nrg_data) const {
  return static_cast<double>(sumPotentialEnergyAsLlint(nrg_data)) * inverse_nrg_scale_lf;
}

//-------------------------------------------------------------------------------------------------
llint ScoreCard::sumTotalEnergyAsLlint(const llint* nrg_data) const {
  llint lacc = nrg_data[static_cast<size_t>(StateVariable::BOND)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::ANGLE)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::PROPER_DIHEDRAL)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::IMPROPER_DIHEDRAL)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::UREY_BRADLEY)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::CHARMM_IMPROPER)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::CMAP)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::VDW)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::VDW_ONE_FOUR)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::ELECTROSTATIC)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::ELEC_ONE_FOUR)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::GENERALIZED_BORN)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::RESTRAINT)];
  lacc += nrg_data[static_cast<size_t>(StateVariable::KINETIC)];
  return lacc;
}

//-------------------------------------------------------------------------------------------------
double ScoreCard::sumTotalEnergy(const llint* nrg_data) const {
  const llint lacc = nrg_data[static_cast<size_t>(StateVariable::KINETIC)];
  return static_cast<double>(sumTotalEnergyAsLlint(nrg_data)) * inverse_nrg_scale_lf;
}

//-------------------------------------------------------------------------------------------------
int ScoreCard::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
int ScoreCard::getSampleSize() const {
  return sampled_step_count;
}

//-------------------------------------------------------------------------------------------------
int ScoreCard::getDataStride() const {
  return data_stride;
}

//-------------------------------------------------------------------------------------------------
int ScoreCard::getSampleCapacity() const {
  return sample_capacity;
}

//-------------------------------------------------------------------------------------------------
int ScoreCard::getEnergyScaleBits() const {
  return nrg_scale_bits;
}

//-------------------------------------------------------------------------------------------------
int ScoreCard::getTimeStep(const int time_index) const {
  return time_steps.readHost(time_index);
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void ScoreCard::upload() {
  time_steps.upload();
  instantaneous_accumulators.upload();
  running_accumulators.upload();
  squared_accumulators.upload();
  time_series_accumulators.upload();
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::download() {
  time_steps.download();
  instantaneous_accumulators.download();
  running_accumulators.download();
  squared_accumulators.download();
  time_series_accumulators.download();
}
#endif

//-------------------------------------------------------------------------------------------------
const ScoreCardReader ScoreCard::data(const HybridTargetLevel tier) const {
  return ScoreCardReader(system_count, data_stride, sampled_step_count, nrg_scale_f, nrg_scale_lf,
                         inverse_nrg_scale_f, inverse_nrg_scale_lf, time_steps.data(tier),
                         instantaneous_accumulators.data(tier), running_accumulators.data(tier),
                         squared_accumulators.data(tier), time_series_accumulators.data(tier));
}

//-------------------------------------------------------------------------------------------------
ScoreCardWriter ScoreCard::data(const HybridTargetLevel tier) {
  return ScoreCardWriter(system_count, data_stride, sampled_step_count, nrg_scale_f, nrg_scale_lf,
                         inverse_nrg_scale_f, inverse_nrg_scale_lf, time_steps.data(tier),
                         instantaneous_accumulators.data(tier), running_accumulators.data(tier),
                         squared_accumulators.data(tier), time_series_accumulators.data(tier));
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::reserve(const int new_capacity) {
  if (new_capacity > sample_capacity) {

    // Allocate one space beyond the requested capacity, so that reserve() can be called at the
    // outset of a molecular dynamics simulation for the projected number of sampled steps and
    // not trigger a resize when the final step gets logged.
    sample_capacity = new_capacity + 1;
    time_series_accumulators.resize(static_cast<size_t>(system_count * data_stride) *
                                    static_cast<size_t>(sample_capacity));
    time_steps.resize(sample_capacity);
  }
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::contribute(const StateVariable var, const llint amount, const int system_index) {
  const size_t slot = static_cast<int>(var) + (data_stride * system_index);
  instantaneous_accumulators.putHost(amount, slot);
  const double dbl_amount = static_cast<double>(amount) * inverse_nrg_scale_lf;
  running_accumulators.putHost(running_accumulators.readHost(slot) + dbl_amount, slot);
  squared_accumulators.putHost(squared_accumulators.readHost(slot) + dbl_amount * dbl_amount,
                               slot);
  const size_t ts_slot = slot + (static_cast<size_t>(data_stride * system_count) *
                                 static_cast<size_t>(sampled_step_count));
  if (sampled_step_count >= sample_capacity) {
    if (sample_capacity < 100) {
      reserve(2LLU * static_cast<size_t>(sampled_step_count));
    }
    else {
      reserve(5LLU * static_cast<size_t>(sampled_step_count) / 4LLU);
    }
  }
  time_series_accumulators.putHost(amount, ts_slot);
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::initialize(const StateVariable var, const int system_index,
                           const HybridTargetLevel tier, const GpuDetails &gpu) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const size_t slot = static_cast<int>(var) + (data_stride * system_index);
      instantaneous_accumulators.putHost(0LL, slot);
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    launchScoreCardInitialization(var, system_index,
                                  instantaneous_accumulators.data(HybridTargetLevel::DEVICE),
                                  system_count, gpu); 
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::initialize(const std::vector<StateVariable> &var, const int system_index,
                           const HybridTargetLevel tier, const GpuDetails &gpu) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const size_t nvar = var.size();
      for (size_t i = 0LLU; i < nvar; i++) {
        initialize(var[i], system_index);
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    launchScoreCardInitialization(var, system_index,
                                  instantaneous_accumulators.data(HybridTargetLevel::DEVICE),
                                  system_count, gpu); 
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::initialize(const int system_index, const HybridTargetLevel tier,
                           const GpuDetails &gpu) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const int last = static_cast<int>(StateVariable::ALL_STATES);
      llint* data_ptr = instantaneous_accumulators.data();
      for (int i = 0; i < last; i++) {
        data_ptr[(system_index * data_stride) + i] = 0LL;    
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    launchScoreCardInitialization(StateVariable::ALL_STATES, system_index,
                                  instantaneous_accumulators.data(HybridTargetLevel::DEVICE),
                                  system_count, gpu); 
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::initialize(const HybridTargetLevel tier, const GpuDetails &gpu) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const int last = static_cast<int>(StateVariable::ALL_STATES);
      llint* data_ptr = instantaneous_accumulators.data();
      for (int i = 0; i < system_count; i++) {
        for (int j = 0; j < last; j++) {
          data_ptr[(i * data_stride) + j] = 0LL;
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    launchScoreCardInitialization(StateVariable::ALL_STATES, -1,
                                  instantaneous_accumulators.data(HybridTargetLevel::DEVICE),
                                  system_count, gpu); 
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::add(const StateVariable var, const llint amount, const int system_index) {
  const size_t slot = static_cast<int>(var) + (data_stride * system_index);
  instantaneous_accumulators.putHost(instantaneous_accumulators.readHost(slot) + amount, slot);
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::commit(const StateVariable var, const int system_index,
                       const HybridTargetLevel tier, const GpuDetails &gpu) {
  if (sampled_step_count == sample_capacity) {
    sample_capacity *= 2;
    time_series_accumulators.resize(static_cast<size_t>(system_count * data_stride) *
                                    static_cast<size_t>(sample_capacity));
  }
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const size_t slot = static_cast<int>(var) + (data_stride * system_index);
      const llint amount = instantaneous_accumulators.readHost(slot);
      const double dbl_amount = static_cast<double>(amount) * inverse_nrg_scale_lf;
      running_accumulators.putHost(running_accumulators.readHost(slot) + dbl_amount, slot);
      squared_accumulators.putHost(squared_accumulators.readHost(slot) + (dbl_amount * dbl_amount),
                                   slot);
      const size_t ts_slot = slot + (static_cast<size_t>(data_stride * system_count) *
                                     static_cast<size_t>(sampled_step_count));
      time_series_accumulators.putHost(amount, ts_slot);
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    launchScoreCardCommit(var, -1, system_count, sampled_step_count,
                          instantaneous_accumulators.data(tier), running_accumulators.data(tier),
                          squared_accumulators.data(tier), time_series_accumulators.data(tier),
                          gpu);
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::commit(const std::vector<StateVariable> &var, const int system_index,
                       const HybridTargetLevel tier, const GpuDetails &gpu) {
  if (sampled_step_count == sample_capacity) {
    sample_capacity *= 2;
    time_series_accumulators.resize(static_cast<size_t>(system_count * data_stride) *
                                    static_cast<size_t>(sample_capacity));
  }
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const llint* inst_acc_ptr = instantaneous_accumulators.data();
      double* run_acc_ptr = running_accumulators.data();
      double* sqd_acc_ptr = squared_accumulators.data();
      llint* time_ser_ptr = time_series_accumulators.data();
      const size_t nvar = var.size();
      for (size_t i = 0LLU; i < nvar; i++) {
        const size_t slot = static_cast<int>(var[i]) + (data_stride * system_index);
        const llint amount = inst_acc_ptr[slot];
        const double dbl_amount = static_cast<double>(amount) * inverse_nrg_scale_lf;
        run_acc_ptr[slot] += dbl_amount;
        sqd_acc_ptr[slot] += dbl_amount * dbl_amount;
        const size_t ts_slot = slot + (static_cast<size_t>(data_stride * system_count) *
                                       static_cast<size_t>(sampled_step_count));
        time_ser_ptr[ts_slot] = amount;
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    launchScoreCardCommit(var, -1, system_count, sampled_step_count,
                          instantaneous_accumulators.data(tier), running_accumulators.data(tier),
                          squared_accumulators.data(tier), time_series_accumulators.data(tier),
                          gpu);
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::commit(const StateVariable var, const HybridTargetLevel tier,
                       const GpuDetails &gpu) {

  // The ALL_STATES variable triggers the commit for all details
  if (var == StateVariable::ALL_STATES) {
    commit(tier, gpu);
    return;
  }
  if (sampled_step_count == sample_capacity) {
    sample_capacity *= 2;
    time_series_accumulators.resize(static_cast<size_t>(system_count * data_stride) *
                                    static_cast<size_t>(sample_capacity));
  }
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const llint* inst_acc_ptr = instantaneous_accumulators.data();
      double* run_acc_ptr = running_accumulators.data();
      double* sqd_acc_ptr = squared_accumulators.data();
      llint* time_ser_ptr = time_series_accumulators.data();
      for (int i = 0; i < system_count; i++) {
        const size_t slot = static_cast<int>(var) + (data_stride * i);
        const llint amount = inst_acc_ptr[slot];
        const double dbl_amount = static_cast<double>(amount) * inverse_nrg_scale_lf;
        run_acc_ptr[slot] += dbl_amount;
        sqd_acc_ptr[slot] += dbl_amount * dbl_amount;
        const size_t ts_slot = slot + (static_cast<size_t>(data_stride * system_count) *
                                       static_cast<size_t>(sampled_step_count));
        time_ser_ptr[ts_slot] = amount;
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    launchScoreCardCommit(var, -1, system_count, sampled_step_count,
                          instantaneous_accumulators.data(tier), running_accumulators.data(tier),
                          squared_accumulators.data(tier), time_series_accumulators.data(tier),
                          gpu);
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::commit(const std::vector<StateVariable> &var, const HybridTargetLevel tier,
                       const GpuDetails &gpu) {
  if (sampled_step_count == sample_capacity) {
    sample_capacity *= 2;
    time_series_accumulators.resize(static_cast<size_t>(system_count * data_stride) *
                                    static_cast<size_t>(sample_capacity));
  }
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const llint* inst_acc_ptr = instantaneous_accumulators.data();
      double* run_acc_ptr = running_accumulators.data();
      double* sqd_acc_ptr = squared_accumulators.data();
      llint* time_ser_ptr = time_series_accumulators.data();
      const int nvar = var.size();
      for (int i = 0; i < system_count; i++) {
        for (int j = 0; j < nvar; j++) {
          const size_t slot = static_cast<int>(var[j]) + (data_stride * i);
          const llint amount = inst_acc_ptr[slot];
          const double dbl_amount = static_cast<double>(amount) * inverse_nrg_scale_lf;
          run_acc_ptr[slot] += dbl_amount;
          sqd_acc_ptr[slot] += dbl_amount * dbl_amount;
          const size_t ts_slot = slot + (static_cast<size_t>(data_stride * system_count) *
                                         static_cast<size_t>(sampled_step_count));
          time_ser_ptr[ts_slot] = amount;
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    launchScoreCardCommit(var, -1, system_count, sampled_step_count,
                          instantaneous_accumulators.data(tier), running_accumulators.data(tier),
                          squared_accumulators.data(tier), time_series_accumulators.data(tier),
                          gpu);
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::commit(int system_index, const HybridTargetLevel tier, const GpuDetails &gpu) {
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  std::vector<StateVariable> var_vec(nvar);
  for (int i = 0; i < nvar; i++) {
    var_vec[i] = static_cast<StateVariable>(i);
  }
  commit(var_vec, system_index, tier, gpu);
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::commit(const HybridTargetLevel tier, const GpuDetails &gpu) {
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  std::vector<StateVariable> var_vec(nvar);
  for (int i = 0; i < nvar; i++) {
    var_vec[i] = static_cast<StateVariable>(i);
  }
  commit(var_vec, tier, gpu);
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::computePotentialEnergy(const HybridTargetLevel tier, const GpuDetails &gpu) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const llint* inst_acc_ptr = instantaneous_accumulators.data();
      const llint* tser_acc_ptr = time_series_accumulators.data();
      const int pe_int = static_cast<int>(StateVariable::POTENTIAL_ENERGY);
      const size_t pe_zu = static_cast<size_t>(pe_int);
      for (int i = 0; i < system_count; i++) {
        const llint* iptr = &inst_acc_ptr[i * data_stride];
        instantaneous_accumulators.putHost(sumPotentialEnergyAsLlint(iptr),
                                           (i * data_stride) + pe_int);
        for (int j = 0; j < sampled_step_count; j++) {
          const size_t sys_offset = static_cast<size_t>(data_stride * system_count) *
                                    static_cast<size_t>(j) + static_cast<size_t>(data_stride * i);
          time_series_accumulators.putHost(sumPotentialEnergyAsLlint(&tser_acc_ptr[sys_offset]),
                                           sys_offset + pe_zu);
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    launchScoreCardEnergySum(system_count, sampled_step_count,
                             instantaneous_accumulators.data(tier),
                             time_series_accumulators.data(tier), running_accumulators.data(tier),
                             squared_accumulators.data(tier), inverse_nrg_scale_lf, gpu);
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::computeTotalEnergy(const HybridTargetLevel tier, const GpuDetails &gpu) {
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const llint* inst_acc_ptr = instantaneous_accumulators.data();
      const llint* tser_acc_ptr = time_series_accumulators.data();
      const int te_int = static_cast<int>(StateVariable::TOTAL_ENERGY);
      const size_t te_zu = static_cast<size_t>(te_int);
      for (int i = 0; i < system_count; i++) {
        instantaneous_accumulators.putHost(sumTotalEnergyAsLlint(&inst_acc_ptr[i * data_stride]),
                                           (i * data_stride) + te_int);
        for (int j = 0; j < sampled_step_count; j++) {
          const size_t sys_offset = static_cast<size_t>(data_stride * system_count) *
                                    static_cast<size_t>(j) + static_cast<size_t>(data_stride * i);
          time_series_accumulators.putHost(sumTotalEnergyAsLlint(&tser_acc_ptr[sys_offset]),
                                           sys_offset + te_zu);
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    launchScoreCardEnergySum(system_count, sampled_step_count,
                             instantaneous_accumulators.data(tier),
                             time_series_accumulators.data(tier), running_accumulators.data(tier),
                             squared_accumulators.data(tier), inverse_nrg_scale_lf, gpu);
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::incrementSampleCount() {
  sampled_step_count += 1;
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::resetSampleCount(const int count_in) {
  sampled_step_count = count_in;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportPotentialEnergies(const HybridTargetLevel tier) const {
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  std::vector<double> result(system_count);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const llint* inst_acc_ptr = instantaneous_accumulators.data();
      for (int i = 0; i < system_count; i++) {
        const size_t info_start = (i * padded_nvar);
        result[i] = sumPotentialEnergy(&inst_acc_ptr[info_start]);
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<llint> devc_acc = instantaneous_accumulators.readDevice();
      const llint* inst_acc_ptr = devc_acc.data();
      for (int i = 0; i < system_count; i++) {
        const size_t info_start = (i * padded_nvar);
        result[i] = sumPotentialEnergy(&inst_acc_ptr[info_start]);
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double ScoreCard::reportPotentialEnergy(const int system_index,
                                        const HybridTargetLevel tier) const {
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  const int info_start = system_index * padded_nvar;
  switch (tier) {
  case HybridTargetLevel::HOST:
    return sumPotentialEnergy(&instantaneous_accumulators.data()[info_start]);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<llint> devc_acc = instantaneous_accumulators.readDevice(info_start,
                                                                                padded_nvar);
      return sumPotentialEnergy(devc_acc.data());
    }
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportTotalEnergies(const HybridTargetLevel tier) const {
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  std::vector<double> result(system_count);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const llint* inst_acc_ptr = instantaneous_accumulators.data();
      for (int i = 0; i < system_count; i++) {
        const size_t info_start = (i * padded_nvar);
        result[i] = sumTotalEnergy(&inst_acc_ptr[info_start]);
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<llint> devc_acc = instantaneous_accumulators.readDevice();
      const llint* inst_acc_ptr = devc_acc.data();
      for (int i = 0; i < system_count; i++) {
        const size_t info_start = (i * padded_nvar);
        result[i] = sumTotalEnergy(&inst_acc_ptr[info_start]);
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double ScoreCard::reportTotalEnergy(const int system_index, const HybridTargetLevel tier) const {
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  const int info_start = system_index * padded_nvar;
  switch (tier) {
  case HybridTargetLevel::HOST:
    return sumTotalEnergy(&instantaneous_accumulators.data()[info_start]);
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<llint> devc_acc = instantaneous_accumulators.readDevice(info_start,
                                                                                padded_nvar);
      return sumTotalEnergy(devc_acc.data());
    }
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportInstantaneousStates(const HybridTargetLevel tier) const {
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  std::vector<double> result(nvar * system_count);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const llint* inst_acc_ptr = instantaneous_accumulators.data();
      for (int i = 0; i < system_count; i++) {
        for (int j = 0; j < nvar; j++) {
          const llint lacc = inst_acc_ptr[(i * padded_nvar) + j];
          result[(i * nvar) + j] = inverse_nrg_scale_lf * static_cast<double>(lacc);
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<llint> devc_acc = instantaneous_accumulators.readDevice();
      for (int i = 0; i < system_count; i++) {
        for (int j = 0; j < nvar; j++) {
          const llint lacc = devc_acc[(i * padded_nvar) + j];
          result[(i * nvar) + j] = inverse_nrg_scale_lf * static_cast<double>(lacc);
        }
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportInstantaneousStates(const int system_index,
                                                         const HybridTargetLevel tier) const {
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int offset = system_index * roundUp(nvar, warp_size_int);
  std::vector<double> result(nvar);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const llint* inst_acc_ptr = instantaneous_accumulators.data();
      for (int i = 0; i < nvar; i++) {
        result[i] = inverse_nrg_scale_lf * static_cast<double>(inst_acc_ptr[offset + i]);
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<llint> devc_acc = instantaneous_accumulators.readDevice();
      for (int i = 0; i < nvar; i++) {
        result[i] = inverse_nrg_scale_lf * static_cast<double>(devc_acc[offset + i]);
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportInstantaneousStates(const StateVariable aspect,
                                                         const HybridTargetLevel tier) const {
  const int aspect_no = static_cast<int>(aspect);
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  std::vector<double> result(system_count);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const llint* inst_acc_ptr = instantaneous_accumulators.data();
      for (int i = 0; i < system_count; i++) {
        const llint lacc = inst_acc_ptr[(i * padded_nvar) + aspect_no];
        result[i] = inverse_nrg_scale_lf * static_cast<double>(lacc);
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<llint> devc_acc = instantaneous_accumulators.readDevice();
      for (int i = 0; i < system_count; i++) {
        const llint lacc = devc_acc[(i * padded_nvar) + aspect_no];
        result[i] = inverse_nrg_scale_lf * static_cast<double>(lacc);
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double ScoreCard::reportInstantaneousStates(const StateVariable aspect, const int system_index,
                                            const HybridTargetLevel tier) const {
  const int aspect_no = static_cast<int>(aspect);
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  const int didx = (system_index * padded_nvar) + aspect_no;
  switch (tier) {
  case HybridTargetLevel::HOST:
    return inverse_nrg_scale_lf * static_cast<double>(instantaneous_accumulators.readHost(didx));
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return inverse_nrg_scale_lf * static_cast<double>(instantaneous_accumulators.readDevice(didx));
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportAverageStates(const HybridTargetLevel tier) const {
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  const double nsamp = static_cast<double>(sampled_step_count);
  std::vector<double> result(nvar * system_count);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const double* run_acc_ptr = running_accumulators.data();
      for (int i = 0; i < system_count; i++) {
        for (int j = 0; j < nvar; j++) {
          result[(i * nvar) + j] = run_acc_ptr[(i * padded_nvar) + j] / nsamp;
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<double> devc_acc = running_accumulators.readDevice();
      for (int i = 0; i < system_count; i++) {
        for (int j = 0; j < nvar; j++) {
          result[(i * nvar) + j] = devc_acc[(i * padded_nvar) + j] / nsamp;
        }
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportAverageStates(const int system_index,
                                                   const HybridTargetLevel tier) const {
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int offset = system_index * roundUp(nvar, warp_size_int);
  const double nsamp = static_cast<double>(sampled_step_count);
  std::vector<double> result(nvar);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const double* run_acc_ptr = running_accumulators.data();
      for (int i = 0; i < nvar; i++) {
        result[i] = run_acc_ptr[offset + i] / nsamp;
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<double> devc_acc = running_accumulators.readDevice();
      for (int i = 0; i < nvar; i++) {
        result[i] = devc_acc[offset + i] / nsamp;
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportAverageStates(const StateVariable aspect,
                                                   const HybridTargetLevel tier) const {
  const int aspect_no = static_cast<int>(aspect);
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  const double nsamp = static_cast<double>(sampled_step_count);
  std::vector<double> result(system_count);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const double* run_acc_ptr = running_accumulators.data();
      for (int i = 0; i < system_count; i++) {
        result[i] = run_acc_ptr[(i * padded_nvar) + aspect_no] / nsamp;
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<double> devc_acc = running_accumulators.readDevice();
      const double* run_acc_ptr = devc_acc.data();
      for (int i = 0; i < system_count; i++) {
        result[i] = run_acc_ptr[(i * padded_nvar) + aspect_no] / nsamp;
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double ScoreCard::reportAverageStates(const StateVariable aspect, const int system_index,
                                      const HybridTargetLevel tier) const {
  const int aspect_no = static_cast<int>(aspect);
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int offset = system_index * roundUp(nvar, warp_size_int);
  const double nsamp = static_cast<double>(sampled_step_count);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return running_accumulators.readHost(offset + aspect_no) / nsamp;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return running_accumulators.readDevice(offset + aspect_no) / nsamp;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportVarianceOfStates(const HybridTargetLevel tier) const {
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  const double nsamp = static_cast<double>(sampled_step_count);
  std::vector<double> result(nvar * system_count, 0.0);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const double* run_acc_ptr = running_accumulators.data();
      const double* sqr_acc_ptr = squared_accumulators.data();
      for (int i = 0; i < system_count; i++) {
        for (int j = 0; j < nvar; j++) {
          const double s1 = run_acc_ptr[(i * padded_nvar) + j];
          const double s2 = sqr_acc_ptr[(i * padded_nvar) + j];
          result[(i * nvar) + j] = sqrt((nsamp * s2) - (s1 * s1)) / sqrt(nsamp * (nsamp - 1.0));
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<double> devc_run_acc = running_accumulators.readDevice();
      const std::vector<double> devc_sqr_acc = running_accumulators.readDevice();
      for (int i = 0; i < system_count; i++) {
        for (int j = 0; j < nvar; j++) {
          const double s1 = devc_run_acc[(i * padded_nvar) + j];
          const double s2 = devc_sqr_acc[(i * padded_nvar) + j];
          result[(i * nvar) + j] = sqrt((nsamp * s2) - (s1 * s1)) / sqrt(nsamp * (nsamp - 1.0));
        }
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportVarianceOfStates(const int system_index,
                                                      const HybridTargetLevel tier) const {
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int offset = system_index * roundUp(nvar, warp_size_int);
  const double nsamp = static_cast<double>(sampled_step_count);
  std::vector<double> result(nvar);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const double* run_acc_ptr = running_accumulators.data();
      const double* sqr_acc_ptr = squared_accumulators.data();
      for (int i = 0; i < nvar; i++) {
        const double s1 = run_acc_ptr[offset + i];
        const double s2 = sqr_acc_ptr[offset + i];
        result[i] = sqrt((nsamp * s2) - (s1 * s1)) / sqrt(nsamp * (nsamp - 1.0));
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<double> devc_run_acc = running_accumulators.readDevice();
      const std::vector<double> devc_sqr_acc = running_accumulators.readDevice();
      for (int i = 0; i < nvar; i++) {
        const double s1 = devc_run_acc[offset + i];
        const double s2 = devc_sqr_acc[offset + i];
        result[i] = sqrt((nsamp * s2) - (s1 * s1)) / sqrt(nsamp * (nsamp - 1.0));
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportVarianceOfStates(const StateVariable aspect,
                                                      const HybridTargetLevel tier) const {
  const int aspect_no = static_cast<int>(aspect);
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int padded_nvar = roundUp(nvar, warp_size_int);
  const double nsamp = static_cast<double>(sampled_step_count);
  std::vector<double> result(system_count, 0.0);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const double* run_acc_ptr = running_accumulators.data();
      const double* sqr_acc_ptr = squared_accumulators.data();
      for (int i = 0; i < system_count; i++) {
        const double s1 = run_acc_ptr[(i * padded_nvar) + aspect_no];
        const double s2 = sqr_acc_ptr[(i * padded_nvar) + aspect_no];
        result[i] = sqrt((nsamp * s2) - (s1 * s1)) / sqrt(nsamp * (nsamp - 1.0));
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const std::vector<double> devc_run_acc = running_accumulators.readDevice();
      const std::vector<double> devc_sqr_acc = running_accumulators.readDevice();
      for (int i = 0; i < system_count; i++) {
        const double s1 = devc_run_acc[(i * padded_nvar) + aspect_no];
        const double s2 = devc_sqr_acc[(i * padded_nvar) + aspect_no];
        result[i] = sqrt((nsamp * s2) - (s1 * s1)) / sqrt(nsamp * (nsamp - 1.0));
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double ScoreCard::reportVarianceOfStates(const StateVariable aspect, const int system_index,
                                         const HybridTargetLevel tier) const {
  const int aspect_no = static_cast<int>(aspect);
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int offset = system_index * roundUp(nvar, warp_size_int);
  const double nsamp = static_cast<double>(sampled_step_count);
  double s1, s2;
  switch (tier) {
  case HybridTargetLevel::HOST:
    s1 = running_accumulators.readHost(offset + aspect_no);
    s2 = squared_accumulators.readHost(offset + aspect_no);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    s1 = running_accumulators.readDevice(offset + aspect_no);
    s2 = squared_accumulators.readDevice(offset + aspect_no);    
    break;
#endif
  }
  return sqrt((nsamp * s2) - (s1 * s1)) / sqrt(nsamp * (nsamp - 1.0));
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportHistory(const int system_index,
                                             const HybridTargetLevel tier) const {
  return reportHistory(StateVariable::TOTAL_ENERGY, system_index, tier);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ScoreCard::reportHistory(const StateVariable aspect, const int system_index,
                                             const HybridTargetLevel tier) const {
  std::vector<double> result(sampled_step_count);
  const llint* ts_ptr = time_series_accumulators.data();
  const size_t data_stride_zu = static_cast<size_t>(data_stride);
  for (int i = 0; i < sampled_step_count; i++) {
    const size_t info_start = static_cast<size_t>(system_index + (i * system_count)) *
                              data_stride_zu;
    result[i] = static_cast<double>(ts_ptr[info_start + static_cast<size_t>(aspect)]) *
                inverse_nrg_scale_lf;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double2> ScoreCard::reportHistory(const HybridTargetLevel tier) const {
  return reportHistory(StateVariable::TOTAL_ENERGY, tier);
}

//-------------------------------------------------------------------------------------------------
std::vector<double2> ScoreCard::reportHistory(const StateVariable aspect,
                                              const HybridTargetLevel tier) const {
  std::vector<double2> result(sampled_step_count);
  const llint* ts_ptr = time_series_accumulators.data();
  const size_t data_stride_zu = static_cast<size_t>(data_stride);
  std::vector<double> step_buffer(system_count);
  for (int i = 0; i < sampled_step_count; i++) {
    for (int j = 0; j < system_count; j++) {
      const size_t info_start = static_cast<size_t>(j + (i * system_count)) * data_stride_zu;
      if (aspect == StateVariable::TOTAL_ENERGY) {
        step_buffer[j] = sumTotalEnergy(&ts_ptr[info_start]);
      }
      else if (aspect == StateVariable::POTENTIAL_ENERGY) {
        step_buffer[j] = sumPotentialEnergy(&ts_ptr[info_start]);
      }
      else {
        step_buffer[j] = static_cast<double>(ts_ptr[info_start + static_cast<size_t>(aspect)]) *
                         inverse_nrg_scale_lf;
      }
    }
    result[i].x = mean(step_buffer);
    result[i].y = variance(step_buffer, VarianceMethod::STANDARD_DEVIATION);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const ScoreCard* ScoreCard::getSelfPointer() const {
  return this;
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::setLastTimeStep(const int time_index, const HybridTargetLevel tier) {
  if (sampled_step_count <= 0) {
    rtErr("No energy values have yet been committed to the sampled history.", "ScoreCard",
          "setLastTimeStep");
  }
  switch (tier) {
  case HybridTargetLevel::HOST:
    time_steps.putHost(time_index, sampled_step_count - 1);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    time_steps.putDevice(time_index, sampled_step_count - 1);
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::import(const ScoreCard *other, const size_t fill_index,
                       const size_t source_index) {
  ScoreCardReader src_r = other->data();
  if (source_index >= src_r.system_count) {
    rtErr("System index " + std::to_string(source_index) + " is invalid for a source ScoreCard "
          "with " + std::to_string(src_r.system_count) + " systems.", "ScoreCard", "import");
  }
  if (fill_index >= system_count) {
    rtErr("System index " + std::to_string(fill_index) + " is invalid for a ScoreCard with " +
          std::to_string(system_count) + " systems.", "ScoreCard", "import");
  }

  // Take the other object's time steps if there are no settings
  bool all_zero_ts = true;
  int tscon = 0;
  int* time_step_ptr = time_steps.data();
  while (tscon < sampled_step_count && all_zero_ts) {
    all_zero_ts = (all_zero_ts && time_step_ptr[tscon] == 0);
    tscon++;
  }
  bool other_zero_ts = true;
  tscon = 0;
  while (tscon < src_r.sampled_step_count && other_zero_ts) {
    other_zero_ts = (other_zero_ts && src_r.time_steps[tscon] == 0);
    tscon++;
  }

  // Reserve space as needed.  Apply the other object's time series if this object has no time
  // series of its own.  Otherwise, check that this object's time series matches the other
  // object's, or that at least the other object has no time series of its own.  Warn if there is
  // an explicit mismatch.
  if (src_r.sampled_step_count > sample_capacity) {
    reserve(src_r.sampled_step_count);
  }
  if (all_zero_ts == false && other_zero_ts == false) {
    bool ts_consistent = true;
    const int minpts = std::min(sampled_step_count, src_r.sampled_step_count);
    for (int i = 0; i < minpts; i++) {
      ts_consistent = (ts_consistent && (time_step_ptr[i] == src_r.time_steps[i]));
    }
    if (ts_consistent == false) {
      rtWarn("Importing energy tracking with inconsistent time step numbering.  The first " +
             std::to_string(src_r.sampled_step_count) + " time steps will be altered to match the "
             "imported energies.", "ScoreCard", "import");
    }
    const int istart = (ts_consistent) ? minpts : 0;
    for (int i = istart; i < src_r.sampled_step_count; i++) {
      time_step_ptr[i] = src_r.time_steps[i];
    }
  }
  else if (all_zero_ts && other_zero_ts == false) {
    for (int i = 0; i < src_r.sampled_step_count; i++) {
      time_step_ptr[i] = src_r.time_steps[i];
    }
  }
  sampled_step_count = std::max(sampled_step_count, src_r.sampled_step_count);
  const size_t nvar = data_stride;
  const size_t padded_nvar = data_stride;
  size_t src_llim = padded_nvar * source_index;
  size_t fill_llim = padded_nvar * fill_index;
  llint* inst_acc_ptr = instantaneous_accumulators.data();
  double* run_acc_ptr = running_accumulators.data();
  double* sqd_acc_ptr = squared_accumulators.data();
  llint* ts_acc_ptr = time_series_accumulators.data();
  const int bit_shift = nrg_scale_bits - other->getEnergyScaleBits();
  const llint bit_fac = llround(pow(2.0, fabs(bit_shift)));
  for (size_t i = 0; i < nvar; i++) {
    if (bit_shift < 0) {
      inst_acc_ptr[fill_llim + i] = (src_r.instantaneous_accumulators[src_llim + i] / bit_fac);
    }
    else if (bit_shift == 0) {
      inst_acc_ptr[fill_llim + i] = src_r.instantaneous_accumulators[src_llim + i];
    }
    else {
      inst_acc_ptr[fill_llim + i] = (src_r.instantaneous_accumulators[src_llim + i] * bit_fac);
    }
    run_acc_ptr[fill_llim + i] = src_r.running_accumulators[src_llim + i];
    sqd_acc_ptr[fill_llim + i] = src_r.squared_accumulators[src_llim + i];
  }
  const size_t nsys = system_count;
  for (size_t i = 0; i < src_r.sampled_step_count; i++) {
    src_llim  = padded_nvar * ((src_r.system_count * i) + source_index);
    fill_llim = padded_nvar * ((nsys * i) + fill_index);
    for (size_t j = 0; j < nvar; j++) {
      if (bit_shift < 0) {
        ts_acc_ptr[fill_llim + j] = (src_r.time_series[src_llim + j] / bit_fac);
      }
      else if (bit_shift == 0) {
        ts_acc_ptr[fill_llim + j] = src_r.time_series[src_llim + j];
      }
      else {
        ts_acc_ptr[fill_llim + j] = (src_r.time_series[src_llim + j] * bit_fac);
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::import(const ScoreCard &other, const size_t fill_index,
                       const size_t source_index) {
  import(other.getSelfPointer(), fill_index, source_index);
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::import(const ScoreCard *other, const std::vector<int> &fill_indices,
                       const std::vector<int> &source_indices) {
  const size_t nfill = fill_indices.size();
  if (nfill != source_indices.size()) {
    rtErr("A series of " + std::to_string(source_indices.size()) + " source systems cannot fill " +
          std::to_string(nfill) + ".", "ScoreCard", "import"); 
  }
  for (size_t i = 0; i < nfill; i++) {
    import(other, fill_indices[i], source_indices[i]);
  }
}

//-------------------------------------------------------------------------------------------------
void ScoreCard::import(const ScoreCard &other, const std::vector<int> &fill_indices,
                       const std::vector<int> &source_indices) {
  import(other.getSelfPointer(), fill_indices, source_indices);
}

} // namespace energy
} // namespace stormm
