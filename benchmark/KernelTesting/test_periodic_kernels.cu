// -*-c++-*-
#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Accelerator/core_kernel_manager.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Constants/fixed_precision.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Namelists/command_line_parser.h"
#include "../../src/Potential/cellgrid.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Random/hpc_random.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/synthesis_abstracts.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/trajectory_enumerators.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::energy;
using namespace stormm::namelist;
using namespace stormm::numerics;
using namespace stormm::random;
using namespace stormm::symbols;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::trajectory;
using stormm::symbols::twopi;
using stormm::symbols::twopi_f;

#include "../../src/Math/rounding.cui"
#include "../../src/Random/xor_shift_rng.cui"
#include "../../src/Numerics/accumulation.cui"

//-------------------------------------------------------------------------------------------------
// Scramble the coordinates of each system in a synthesis based on the original positions in a
// solitary coordinate object.
//
// Arguments:
//   cfr:        Pristine coordinates of the original system, in double-precision
//   poly_psw:   Coordinates of the synthesis at either stage of the coordinate cycle.  The
//               synthesis is assumed to consist of only one system.
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(small_block_size, 4)
kScrambleCoordinates(PsSynthesisWriter poly_psw, const CoordinateFrameReader cfr,
                     ullint2* rng_states) {

  // Each thread pulls its private random number generator from the global array.
  const int global_tidx = (blockIdx.x * blockDim.x) + threadIdx.x;
  const int total_natom = cfr.natom * poly_psw.system_count;
  const int padded_total_natom = devcRoundUp(total_natom, warp_size_int);
  ullint2 state_v;
  if (global_tidx < 3 * padded_total_natom) {
    state_v = rng_states[global_tidx];
  }
  else {
    state_v = { 0LLU, 0LLU };
  }
  int i = global_tidx;
  while (i < padded_total_natom) {
    if (i < total_natom) {
      const int sys_idx = i / cfr.natom;
      const int atom_idx = i - (cfr.natom * sys_idx);
      const size_t local_idx = atom_idx;
      const size_t synth_idx = __ldca(&poly_psw.atom_starts[sys_idx]) + atom_idx;
      const float pert_x = (float)(0.1) * xoroshiro128p_normalf(&state_v) +
                           __ldca(&cfr.xcrd[local_idx]);
      if (poly_psw.gpos_bits <= globalpos_scale_nonoverflow_bits) {
        const llint ipert_x = __float2ll_rn(pert_x * poly_psw.gpos_scale_f);
        __stwt(&poly_psw.xcrd[synth_idx], ipert_x);
      }
      else {
        const int95_t ipert_x = doubleToInt95(pert_x * poly_psw.gpos_scale_f);
        __stwt(&poly_psw.xcrd[synth_idx], ipert_x.x);
        __stwt(&poly_psw.xcrd_ovrf[synth_idx], ipert_x.y);
      }
    }
    i += blockDim.x * gridDim.x;
  }
  while (i < 2 * padded_total_natom) {
    const int rel_i = i - padded_total_natom;
    if (rel_i < total_natom) {
      const int sys_idx = rel_i / cfr.natom;
      const int atom_idx = rel_i - (cfr.natom * sys_idx);
      const size_t local_idx = atom_idx;
      const size_t synth_idx = __ldca(&poly_psw.atom_starts[sys_idx]) + atom_idx;
      const float pert_y = (float)(0.1) * xoroshiro128p_normalf(&state_v) +
                           __ldca(&cfr.ycrd[local_idx]);
      if (poly_psw.gpos_bits <= globalpos_scale_nonoverflow_bits) {
        const llint ipert_y = __float2ll_rn(pert_y * poly_psw.gpos_scale_f);
        __stwt(&poly_psw.ycrd[synth_idx], ipert_y);
      }
      else {
        const int95_t ipert_y = doubleToInt95(pert_y * poly_psw.gpos_scale_f);
        __stwt(&poly_psw.ycrd[synth_idx], ipert_y.x);
        __stwt(&poly_psw.ycrd_ovrf[synth_idx], ipert_y.y);
      }
    }
    i += blockDim.x * gridDim.x;
  }
  while (i < 3 * padded_total_natom) {
    const int rel_i = i - (2 * padded_total_natom);
    if (rel_i < total_natom) {
      const int sys_idx = rel_i / cfr.natom;
      const int atom_idx = rel_i - (cfr.natom * sys_idx);
      const size_t local_idx = atom_idx;
      const size_t synth_idx = __ldca(&poly_psw.atom_starts[sys_idx]) + atom_idx;
      const float pert_z = (float)(0.1) * xoroshiro128p_normalf(&state_v) +
                           __ldca(&cfr.zcrd[local_idx]);
      if (poly_psw.gpos_bits <= globalpos_scale_nonoverflow_bits) {
        const llint ipert_z = __float2ll_rn(pert_z * poly_psw.gpos_scale_f);
        __stwt(&poly_psw.zcrd[synth_idx], ipert_z);
      }
      else {
        const int95_t ipert_z = doubleToInt95(pert_z * poly_psw.gpos_scale_f);
        __stwt(&poly_psw.zcrd[synth_idx], ipert_z.x);
        __stwt(&poly_psw.zcrd_ovrf[synth_idx], ipert_z.y);
      }
    }
    i += blockDim.x * gridDim.x;
  }

  // Store the random number generator state, having advanced it
  if (global_tidx < 3 * padded_total_natom) {
    rng_states[global_tidx] = state_v;
  }
}

//-------------------------------------------------------------------------------------------------
// Test the migration kernels as applied to a specific system from the collection of test cases.
//
// Arguments:
//   tsm:         The collection of test systems
//   system_idx:  The system index within the collection to examine
//   t_nml:       Namelist object holding user input from the command line
//   timer:       Object for timing particle migration kernels and filtering out other processes
//   gpu:         Details of the GPU on which the resorting will occur
//-------------------------------------------------------------------------------------------------
template <typename T, typename Tacc, typename T4>
void testSystemMigration(const TestSystemManager &tsm, const int system_idx,
                         const NamelistEmulator &t_nml,
                         StopWatch *timer, const GpuDetails &gpu) {
  const int replicas = t_nml.getIntValue("-replicas");
  const int iter = t_nml.getIntValue("-iter");
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  const double elec_cutoff = (t_nml.getKeywordStatus("-elec_cutoff") != InputStatus::MISSING) ?
                             t_nml.getRealValue("-elec_cutoff") : t_nml.getRealValue("-cutoff");
  const double vdw_cutoff = (t_nml.getKeywordStatus("-vdw_cutoff") != InputStatus::MISSING) ?
                            t_nml.getRealValue("-vdw_cutoff") : t_nml.getRealValue("-cutoff");
  const double cell_padding = t_nml.getRealValue("-pad");
  const bool use_dual_cg = (fabs(elec_cutoff - vdw_cutoff) > 1.0e-6 ||
                            t_nml.getBoolValue("-dual_cg"));
  const std::vector<int> index_key(replicas, system_idx);
  CoordinateFrame cf = tsm.exportCoordinateFrame(system_idx);
  cf.upload();
  const CoordinateFrameReader cfr = cf.data(devc);
  PhaseSpaceSynthesis poly_ps = tsm.exportPhaseSpaceSynthesis(index_key, 0.0, 518693,
                                                              t_nml.getIntValue("-gpos_bits"));
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(index_key);
  poly_ps.upload();
  poly_ag.upload();
  PsSynthesisWriter white_psw = poly_ps.data(CoordinateCycle::WHITE, devc);
  PsSynthesisWriter black_psw = poly_ps.data(CoordinateCycle::BLACK, devc);
  const PsSynthesisReader white_psr(white_psw);
  const PsSynthesisReader black_psr(black_psw);

  // Create an array of random number generators sufficient for all threads on the GPU to have
  // their own.
  const CoreKlManager launcher(gpu, poly_ag);
  const int nscramble_block = gpu.getSMPCount() * 4;
  Hybrid<ullint2> rng_states(nscramble_block * small_block_size);
  initXoroshiro128pArray(&rng_states, t_nml.getIntValue("-ig_seed"), 25, gpu);
  ullint2* rng_states_ptr = rng_states.data(devc);
  
  // Measure the update time so that it can be subtracted from the particle migration timings.
  const std::string sys_base = getBaseName(tsm.getTopologyFile(system_idx));
  std::string sys_key, sys_ext;
  splitPath(sys_base, &sys_key, &sys_ext);
  if (use_dual_cg) {
    sys_key += " (Dual Grid)";
  }
  const int scramble_time  = timer->addCategory(sys_key + " Scramble");
  const int migration_time = timer->addCategory(sys_key + " Migration");
  timer->assignTime(0);
  for (int i = 0; i < iter; i++) {
    if (i & 0x1) {
      kScrambleCoordinates<<<nscramble_block, small_block_size>>>(white_psw, cfr, rng_states_ptr);
    }
    else {
      kScrambleCoordinates<<<nscramble_block, small_block_size>>>(black_psw, cfr, rng_states_ptr);
    }
  }
  cudaDeviceSynchronize();
  timer->assignTime(scramble_time);  
  const PrecisionModel prec = (std::type_index(typeid(T)).hash_code() == double_type_index) ?
                              PrecisionModel::DOUBLE : PrecisionModel::SINGLE;
  
  // Branch the cases of one or two neighbor lists, as in other programs
  if (use_dual_cg) {
    CellGrid<T, Tacc, T, T4> cg_qq(&poly_ps, poly_ag, (0.5 * elec_cutoff), cell_padding, 4,
                                   NonbondedTheme::ELECTROSTATIC);
    CellGridWriter<T, Tacc, T, T4> white_cg_qqw = cg_qq.data(CoordinateCycle::WHITE, devc);
    CellGridWriter<T, Tacc, T, T4> black_cg_qqw = cg_qq.data(CoordinateCycle::BLACK, devc);
    const CellOriginsReader white_corg_qq = cg_qq.getRulers(CoordinateCycle::WHITE, devc);
    const CellOriginsReader black_corg_qq = cg_qq.getRulers(CoordinateCycle::BLACK, devc);
    CellGrid<T, Tacc, T, T4> cg_lj(&poly_ps, poly_ag, (0.5 * vdw_cutoff), cell_padding, 4,
                                   NonbondedTheme::VAN_DER_WAALS);
    CellGridWriter<T, Tacc, T, T4> white_cg_ljw = cg_lj.data(CoordinateCycle::WHITE, devc);
    CellGridWriter<T, Tacc, T, T4> black_cg_ljw = cg_lj.data(CoordinateCycle::BLACK, devc);
    const CellOriginsReader white_corg_lj = cg_lj.getRulers(CoordinateCycle::WHITE, devc);
    const CellOriginsReader black_corg_lj = cg_lj.getRulers(CoordinateCycle::BLACK, devc);
    cg_qq.upload();
    cg_lj.upload();
    timer->assignTime(0);
    const int chain_count = white_cg_qqw.total_chain_count + white_cg_ljw.total_chain_count;
    const int2 bt_i  = launcher.getMigrationKernelDims(prec, NeighborListKind::DUAL, 1,
                                                       white_psw.gpos_bits, chain_count);
    const int2 bt_ii = launcher.getMigrationKernelDims(prec, NeighborListKind::DUAL, 2,
                                                       white_psw.gpos_bits, chain_count);
    for (int i = 0; i < iter; i++) {
      if (i & 0x1) {
        kScrambleCoordinates<<<nscramble_block, small_block_size>>>(white_psw, cfr,
                                                                    rng_states_ptr);
        launchMigration(&black_cg_qqw, &black_cg_ljw, black_corg_qq, black_corg_lj, black_psr,
                        bt_i, bt_ii);
      }
      else {
        kScrambleCoordinates<<<nscramble_block, small_block_size>>>(black_psw, cfr,
                                                                    rng_states_ptr);
        launchMigration(&white_cg_qqw, &white_cg_ljw, white_corg_qq, white_corg_lj, white_psr,
                        bt_i, bt_ii);
      }
    }
    cudaDeviceSynchronize();
    timer->assignTime(migration_time);
  }
  else {
    CellGrid<T, Tacc, T, T4> cg(&poly_ps, poly_ag, (0.5 * vdw_cutoff), cell_padding, 4,
                                NonbondedTheme::ALL);
    CellGridWriter<T, Tacc, T, T4> white_cgw = cg.data(CoordinateCycle::WHITE, devc);
    CellGridWriter<T, Tacc, T, T4> black_cgw = cg.data(CoordinateCycle::BLACK, devc);
    const CellOriginsReader white_corg = cg.getRulers(CoordinateCycle::WHITE, devc);
    const CellOriginsReader black_corg = cg.getRulers(CoordinateCycle::BLACK, devc);
    cg.upload();
    timer->assignTime(0);
    const int2 bt_i  = launcher.getMigrationKernelDims(prec, NeighborListKind::MONO, 1,
                                                       white_psw.gpos_bits,
                                                       white_cgw.total_chain_count);
    const int2 bt_ii = launcher.getMigrationKernelDims(prec, NeighborListKind::MONO, 2,
                                                       white_psw.gpos_bits,
                                                       white_cgw.total_chain_count);

    // CHECK
    printf("Launch  I on %3d %3d (gpos_bits = %2d total chains = %4d)\n", bt_i.x, bt_i.y,
           white_psw.gpos_bits, white_cgw.total_chain_count);
    printf("Launch II on %3d %3d\n", bt_ii.x, bt_ii.y);
    // END CHECK
    
    for (int i = 0; i < iter; i++) {
      if (i & 0x1) {
        kScrambleCoordinates<<<nscramble_block, small_block_size>>>(white_psw, cfr,
                                                                    rng_states_ptr);
        launchMigration(&black_cgw, black_corg, black_psr, bt_i, bt_ii);
      }
      else {
        kScrambleCoordinates<<<nscramble_block, small_block_size>>>(black_psw, cfr,
                                                                    rng_states_ptr);
        launchMigration(&white_cgw, white_corg, white_psr, bt_i, bt_ii);
      }
    }
    cudaDeviceSynchronize();
    timer->assignTime(migration_time);
  }
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Engage the testing environment and the GPU
  StopWatch timer;
  HpcConfig gpu_config(ExceptionResponse::WARN);
  std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> engage_gpu(1);

  // Get additional command-line inputs
  CommandLineParser clip("test_periodic_kernels", "A program for testing kernels associated with "
                         "calculations in periodic boundary conditions and molecular dynamics.",
                         { "-timings" });
  clip.addStandardAmberInputs({ "-p", "-c", "-o", "-ig_seed" });
  clip.addStandardBenchmarkingInputs({ "-iter", "-trials", "-cutoff", "-elec_cutoff",
                                       "-vdw_cutoff", "-pad", "-replicas" });
  NamelistEmulator *t_nml = clip.getNamelistPointer();
  t_nml->addKeyword("-pert", NamelistType::REAL, std::to_string(0.1));
  t_nml->addHelp("-pert", "The Gaussian width by which to perturb atomic positions, relative to "
                 "their initial states, between neighbor list migrations.");
  t_nml->addKeyword("-gpos_bits", NamelistType::INTEGER, std::to_string(36));
  t_nml->addHelp("-gpos_bits", "The number of btis after the decimal with which to store particle "
                 "positions");
  t_nml->addKeyword("-dual_cg", NamelistType::BOOLEAN);
  t_nml->addHelp("-dual_cg", "Specify dual cell grids to enforce separate electrostatic "
                 "and van-der Waals neighbor list grids.  Specifying distinct electrostatic and "
                 "van-der Waals cutoffs will also trigger two neighbor lists to be built.");
  TestEnvironment oe(argc, argv, &clip, TmpdirStatus::NOT_REQUIRED, ExceptionResponse::SILENT);
  clip.parseUserInput(argc, argv);
  
  // Create a series of test systems
  const std::vector<std::string> all_top = t_nml->getAllStringValues("-p");
  const std::vector<std::string> all_crd = t_nml->getAllStringValues("-c");
  TestSystemManager tsm(all_top, all_crd);
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    testSystemMigration<float, int, float4>(tsm, i, *t_nml, &timer, gpu);
  }
  
  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());     
}
