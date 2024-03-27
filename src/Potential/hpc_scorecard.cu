// -*-c++-*-
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "Namelists/nml_minimize.h"
#include "hpc_scorecard.h"

namespace stormm {
namespace energy {

using namelist::MinimizeControls;
  
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kInitializeScoreCard(const ullint state_mask, llint* accumulators, const int system_index,
                     const int system_count) {

  // Get the extent of the StateVariable enumerator
  const int max_states = (int)(StateVariable::ALL_STATES);
  const int padded_max_states = (((max_states + warp_size_int - 1) >> warp_bits) << warp_bits);
  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (system_index >= 0) {
    if (warp_idx == 0 && blockIdx.x == 0) {
      for (int i = lane_idx; i < max_states; i += warp_size_int) {
        if ((state_mask >> i) & 0x1LLU) {
          accumulators[(padded_max_states * system_index) + i] = 0LL;
        }
      }
    }
  }
  else {
    int syspos = warp_idx + (blockIdx.x * (blockDim.x >> warp_bits));
    while (syspos < system_count) {
      for (int i = lane_idx; i < max_states; i += warp_size_int) {
        if ((state_mask >> i) & 0x1LLU) {
          accumulators[(padded_max_states * syspos) + i] = 0LL;
        }
      }
      syspos += gridDim.x * (blockDim.x >> warp_bits);
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchScoreCardInitialization(const StateVariable var, const int system_index,
                                          llint* accumulators, const int system_count,
                                          const GpuDetails &gpu) {
  const ullint state_mask = (var == StateVariable::ALL_STATES) ? 0xffffffffffffffffLLU :
                                                                 (0x1 << static_cast<int>(var));
  kInitializeScoreCard<<<gpu.getSMPCount(), large_block_size>>>(state_mask, accumulators,
                                                                system_index, system_count);
}

//-------------------------------------------------------------------------------------------------
extern void launchScoreCardInitialization(const std::vector<StateVariable> &var,
                                          const int system_index, llint* accumulators,
                                          const int system_count, const GpuDetails &gpu) {
  const size_t nvar = var.size();
  ullint state_mask = 0LLU;
  for (size_t i = 0; i < nvar; i++) {
    state_mask |= (0x1 << static_cast<int>(var[i]));
  }
  kInitializeScoreCard<<<gpu.getSMPCount(), large_block_size>>>(state_mask, accumulators,
                                                                system_index, system_count);
}

//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kScoreCardSumEnergy(const int system_count, const int sample_count, llint* inst_acc,
                    llint* time_ser_acc, double* run_acc, double* sqd_acc,
                    const double inv_nrg_scale) {
  const int max_states = (int)(StateVariable::ALL_STATES);
  const int padded_max_states = (((max_states + warp_size_int - 1) >> warp_bits) << warp_bits);
  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);

  // Sum the potential or total energy in the instantaneous accumulators.
  int syspos = warp_idx + (blockIdx.x * (blockDim.x >> warp_bits));
  while (syspos < system_count) {
    llint tot_nrg = 0LL;
    llint ptn_nrg = 0LL;
    for (int i = lane_idx; i < max_states; i += warp_size_int) {
      const int slot = (padded_max_states * syspos) + i;
      if (i == (int)(StateVariable::PROPER_DIHEDRAL) || i == (int)(StateVariable::ANGLE) ||
          i == (int)(StateVariable::IMPROPER_DIHEDRAL) || i == (int)(StateVariable::BOND) || 
          i == (int)(StateVariable::UREY_BRADLEY) || i == (int)(StateVariable::ELECTROSTATIC) ||
          i == (int)(StateVariable::CHARMM_IMPROPER) || i == (int)(StateVariable::CMAP) ||
          i == (int)(StateVariable::ELEC_ONE_FOUR) || i == (int)(StateVariable::VDW) ||
          i == (int)(StateVariable::VDW_ONE_FOUR) || i == (int)(StateVariable::RESTRAINT) ||
          i == (int)(StateVariable::GENERALIZED_BORN)) {
        tot_nrg += inst_acc[slot];
        ptn_nrg += inst_acc[slot];
      }
      if (i == (int)(StateVariable::KINETIC)) {
        tot_nrg += inst_acc[slot];
      }
    }
    WARP_REDUCE_DOWN(tot_nrg);
    WARP_REDUCE_DOWN(ptn_nrg);
    if (lane_idx == 0) {
      inst_acc[(padded_max_states * syspos) + (int)(StateVariable::POTENTIAL_ENERGY)] = ptn_nrg;
      inst_acc[(padded_max_states * syspos) + (int)(StateVariable::TOTAL_ENERGY)]     = tot_nrg;
    }
    syspos += gridDim.x * (blockDim.x >> warp_bits);
  }

  // Sum the potential or total energy in each time series, and contribute to the running sums.
  while (syspos < 2 * system_count) {
    const int rel_pos = syspos - system_count;
    double dt_nrg = 0.0;
    double d2t_nrg = 0.0;
    double dp_nrg = 0.0;
    double d2p_nrg = 0.0;
    const size_t sys_offset = (padded_max_states * rel_pos);
    for (int i = 0; i < sample_count; i++) {
      llint tot_nrg = 0LL;
      llint ptn_nrg = 0LL;
      const size_t batch_offset = (size_t)(i) * (size_t)(padded_max_states * system_count);
      for (int j = lane_idx; j < max_states; j += warp_size_int) {
        const size_t batch_pos = (padded_max_states * rel_pos) + j;
        if (j == (int)(StateVariable::PROPER_DIHEDRAL) || j == (int)(StateVariable::ANGLE) ||
            j == (int)(StateVariable::IMPROPER_DIHEDRAL) || j == (int)(StateVariable::BOND) ||
            j == (int)(StateVariable::UREY_BRADLEY) || j == (int)(StateVariable::ELECTROSTATIC) ||
            j == (int)(StateVariable::CHARMM_IMPROPER) || j == (int)(StateVariable::CMAP) ||
            j == (int)(StateVariable::ELEC_ONE_FOUR) || j == (int)(StateVariable::VDW) ||
            j == (int)(StateVariable::VDW_ONE_FOUR) || j == (int)(StateVariable::RESTRAINT) ||
            j == (int)(StateVariable::GENERALIZED_BORN)) {
          tot_nrg += time_ser_acc[batch_offset + batch_pos];
          ptn_nrg += time_ser_acc[batch_offset + batch_pos];
        }
        if (j == (int)(StateVariable::KINETIC)) {
          tot_nrg += time_ser_acc[batch_offset + batch_pos];
        }
      }
      WARP_REDUCE_DOWN(tot_nrg);
      WARP_REDUCE_DOWN(ptn_nrg);
      if (lane_idx == 0) {
        time_ser_acc[batch_offset + sys_offset +
                     (size_t)(StateVariable::POTENTIAL_ENERGY)] = ptn_nrg;
        time_ser_acc[batch_offset + sys_offset +
                     (size_t)(StateVariable::TOTAL_ENERGY)]     = tot_nrg;
        const double ldt_nrg = (double)(tot_nrg) * inv_nrg_scale;
        const double ldp_nrg = (double)(ptn_nrg) * inv_nrg_scale;
        dt_nrg += ldt_nrg;
        d2t_nrg += ldt_nrg * ldt_nrg;
        dp_nrg += ldp_nrg;
        d2p_nrg += ldp_nrg * ldp_nrg;
      }
    }
    run_acc[sys_offset + (size_t)(StateVariable::POTENTIAL_ENERGY)] = dp_nrg;
    sqd_acc[sys_offset + (size_t)(StateVariable::POTENTIAL_ENERGY)] = d2p_nrg;
    syspos += gridDim.x * (blockDim.x >> warp_bits);
  }
}
  
//-------------------------------------------------------------------------------------------------
extern void launchScoreCardEnergySum(const int system_count, const int sample_count,
                                     llint* inst_acc, llint* time_ser_acc, double* run_acc,
                                     double* sqd_acc, const double inv_nrg_scale,
                                     const GpuDetails &gpu) {
  kScoreCardSumEnergy<<<gpu.getSMPCount(), large_block_size>>>(system_count, sample_count,
                                                               inst_acc, time_ser_acc, run_acc,
                                                               sqd_acc, inv_nrg_scale);
}
  
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(large_block_size, 1)
kScoreCardCommit(const ullint state_mask, const int system_index, const int system_count,
                 const size_t sample_count, const llint* inst_acc, double* run_acc,
                 double* sqd_acc, llint* time_ser_acc) {
  
  // Get the extent of the StateVariable enumerator
  const int max_states = (int)(StateVariable::ALL_STATES);
  const int padded_max_states = (((max_states + warp_size_int - 1) >> warp_bits) << warp_bits);
  const int warp_idx = (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);
  if (system_index >= 0) {
    if (warp_idx == 0 && blockIdx.x == 0) {
      for (int i = lane_idx; i < max_states; i += warp_size_int) {
        if ((state_mask >> i) & 0x1LLU) {
          const size_t slot = (padded_max_states * system_index) + i;
          const llint ll_amt = inst_acc[slot];
          const double d_amt = (double)(ll_amt);
          run_acc[slot] += d_amt;
          sqd_acc[slot] += d_amt * d_amt;
          const size_t ts_slot = slot + (static_cast<size_t>(padded_max_states * system_count) *
                                         sample_count);
          time_ser_acc[ts_slot] = ll_amt;
        }
      }
    }    
  }
  else {
    int syspos = warp_idx + (blockIdx.x * (blockDim.x >> warp_bits));
    while (syspos < system_count) {
      for (int i = lane_idx; i < max_states; i += warp_size_int) {
        if ((state_mask >> i) & 0x1LLU) {
          const size_t slot = (padded_max_states * syspos) + i;
          const llint ll_amt = inst_acc[slot];
          const double d_amt = (double)(ll_amt);
          run_acc[slot] += d_amt;
          sqd_acc[slot] += d_amt * d_amt;
          const size_t ts_slot = slot + (static_cast<size_t>(padded_max_states * system_count) *
                                         sample_count);
          time_ser_acc[ts_slot] = ll_amt;
        }
      }
      syspos += gridDim.x * (blockDim.x >> warp_bits);
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchScoreCardCommit(const std::vector<StateVariable> &var, const int system_index,
                                  const int system_count, const size_t sample_count,
                                  const llint* inst_acc, double* run_acc, double* sqd_acc,
                                  llint* time_ser_acc, const GpuDetails &gpu) {
  const size_t nvar = var.size();
  ullint state_mask = 0LLU;
  for (size_t i = 0; i < nvar; i++) {
    state_mask |= (0x1 << static_cast<int>(var[i]));
  }
  kScoreCardCommit<<<gpu.getSMPCount(), large_block_size>>>(state_mask, system_index, system_count,
                                                            sample_count, inst_acc, run_acc,
                                                            sqd_acc, time_ser_acc);
}

//-------------------------------------------------------------------------------------------------
extern void launchScoreCardCommit(StateVariable var, const int system_index,
                                  const int system_count, const size_t sample_count,
                                  const llint* inst_acc, double* run_acc, double* sqd_acc,
                                  llint* time_ser_acc, const GpuDetails &gpu) {
  const ullint state_mask = (var == StateVariable::ALL_STATES) ? 0xffffffffffffffffLLU :
                                                                 (0x1 << static_cast<int>(var));
  kScoreCardCommit<<<gpu.getSMPCount(), large_block_size>>>(state_mask, system_index, system_count,
                                                            sample_count, inst_acc, run_acc,
                                                            sqd_acc, time_ser_acc);
}

} // namespace energy
} // namespace stormm
