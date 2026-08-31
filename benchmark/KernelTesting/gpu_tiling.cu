#include <cuda.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <nvml.h>
#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/core_kernel_manager.h"
#include "../../src/Accelerator/ptx_macros.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/bspline.h"
#include "../../src/Math/summation.h"
#include "../../src/Math/hpc_summation.cuh"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Random/random.h"
#include "../../src/Random/hpc_random.h"
#include "../../src/MolecularMechanics/mm_controls.h"
#include "../../src/Namelists/command_line_parser.h"
#include "../../src/Namelists/namelist_emulator.h"
#include "../../src/Namelists/nml_files.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/hpc_valence_potential.h"
#include "../../src/Potential/hpc_nonbonded_potential.h"
#include "../../src/Potential/local_exclusionmask.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Structure/local_arrangement.h"
#include "../../src/Synthesis/implicit_solvent_workspace.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/systemcache.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/static_mask_synthesis.h"
#include "../../src/Synthesis/synthesis_abstracts.h"
#include "../../src/Synthesis/nonbonded_workunit.h"
#include "../../src/Synthesis/valence_workunit.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/unit_test_enumerators.h"

using stormm::card::HpcConfig;
using stormm::card::Hybrid;
using stormm::card::HybridFormat;
using stormm::card::HybridTargetLevel;
using stormm::card::GpuDetails;
#ifndef STORMM_USE_HPC
using stormm::data_types::int4;
#endif
using stormm::data_types::int95_t;
using stormm::data_types::isFloatingPointScalarType;
using stormm::data_types::llint;
using stormm::data_types::ullint;
using stormm::data_types::ullint2;
using stormm::data_types::ullint4;
using stormm::energy::LocalExclusionMask;
using stormm::energy::LocalExclusionMaskReader;
using stormm::energy::testExclusion;
using stormm::errors::rtWarn;
using stormm::review::stormmSplash;
using stormm::testing::TestSystemManager;
using namespace stormm::constants;
using namespace stormm::diskutil;
using namespace stormm::errors;
using namespace stormm::energy;
using namespace stormm::mm;
using namespace stormm::namelist;
using namespace stormm::numerics;
using namespace stormm::parse;
using namespace stormm::restraints;
using namespace stormm::review;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

#include "../../src/Numerics/accumulation.cui"
#include "../../src/Potential/evaluate_localmask.cui"

//-------------------------------------------------------------------------------------------------
// GPU Kernel to compute interactions within a tile
//-------------------------------------------------------------------------------------------------
__global__ void __launch_bounds__(256, 5) 
                computeTileInteractions (PsSynthesisWriter poly_psw,
                                         const SyNonbondedKit<float, float2> poly_nbk,
                                         const int4 *t_assgn, const int num_tiles,
                                         const LocalExclusionMaskReader lemr) {
  __shared__ volatile int total_warps;

  // Get the total number of warps (blocks * threads)
  if(threadIdx.x == 0){
    total_warps = (blockDim.x >> warp_bits) * gridDim.x;
  }
  __syncthreads();
  
  // Calculate the thread's warp index within the block, its lane index within the warp
  const int warp_idx = (blockDim.x >> warp_bits) * blockIdx.x + (threadIdx.x >> warp_bits);
  const int lane_idx = (threadIdx.x & warp_bits_mask_int);

  // Iterate over the tiles, each block handles multiple tiles
  for(int i = warp_idx; i < num_tiles; i += total_warps) {
    int4 instr = t_assgn[i];  // Instruction to fetch atom counts for current tile
    float recv_x, recv_y, recv_z, recv_q;
    int recv_topo_idx, recv_lj_idx;
    
    if (lane_idx < instr.z) {
      if (poly_psw.gpos_bits > globalpos_scale_nonoverflow_bits) {
        recv_x = int95ToDouble(poly_psw.xcrd[instr.x + lane_idx],
                               poly_psw.xcrd_ovrf[instr.x + lane_idx]) * poly_psw.inv_gpos_scale;
        recv_y = int95ToDouble(poly_psw.ycrd[instr.x + lane_idx],
                               poly_psw.ycrd_ovrf[instr.x + lane_idx]) * poly_psw.inv_gpos_scale;
        recv_z = int95ToDouble(poly_psw.zcrd[instr.x + lane_idx],
                               poly_psw.zcrd_ovrf[instr.x + lane_idx]) * poly_psw.inv_gpos_scale;
      }
      else {
        recv_x = (float)poly_psw.xcrd[instr.x + lane_idx] * poly_psw.inv_gpos_scale;
        recv_y = (float)poly_psw.ycrd[instr.x + lane_idx] * poly_psw.inv_gpos_scale;
        recv_z = (float)poly_psw.zcrd[instr.x + lane_idx] * poly_psw.inv_gpos_scale;
      }
      recv_topo_idx = instr.x + lane_idx;
      recv_lj_idx = poly_nbk.lj_idx[instr.x + lane_idx];
      recv_q = poly_nbk.charge[instr.x + lane_idx] * poly_nbk.coulomb;
    }
    else {
      recv_x = (float)1.0;
      recv_y = (float)1.0;
      recv_z = (float)(-lane_idx);
      recv_topo_idx = -1;
      recv_lj_idx = -1;
      recv_q = (float)0.0;
    }
    
    float send_x, send_y, send_z, send_q;
    int send_topo_idx, send_lj_idx;
    ullint send_prof;
    if (lane_idx < instr.w) {
      if (poly_psw.gpos_bits > globalpos_scale_nonoverflow_bits) {
        send_x = int95ToDouble(poly_psw.xcrd[instr.y + lane_idx],
                               poly_psw.xcrd_ovrf[instr.y + lane_idx]) * poly_psw.inv_gpos_scale;
        send_y = int95ToDouble(poly_psw.ycrd[instr.y + lane_idx],
                               poly_psw.ycrd_ovrf[instr.y + lane_idx]) * poly_psw.inv_gpos_scale;
        send_z = int95ToDouble(poly_psw.zcrd[instr.y + lane_idx],
                               poly_psw.zcrd_ovrf[instr.y + lane_idx]) * poly_psw.inv_gpos_scale;
      }
      else {
        send_x = (float)poly_psw.xcrd[instr.y + lane_idx] * poly_psw.inv_gpos_scale;
        send_y = (float)poly_psw.ycrd[instr.y + lane_idx] * poly_psw.inv_gpos_scale;
        send_z = (float)poly_psw.zcrd[instr.y + lane_idx] * poly_psw.inv_gpos_scale;    
      }
      send_topo_idx = instr.y + lane_idx;
      send_lj_idx = poly_nbk.lj_idx[instr.y + lane_idx];
      send_q = poly_nbk.charge[instr.y + lane_idx] * poly_nbk.coulomb;
      send_prof = lemr.profiles[lemr.prof_idx[instr.y + lane_idx]];
    }
    else {
      send_x = (float)1.0;
      send_y = (float)1.0;
      send_z = (float)(-lane_idx);
      send_topo_idx = -1;
      send_lj_idx = -1;
      send_q = (float)0.0;
      send_prof = 0x1fffffffffffffffU;
    }

    float send_acc_fx = 0.0f;
    float send_acc_fy = 0.0f;
    float send_acc_fz = 0.0f;
    float recv_acc_fx = 0.0f;
    float recv_acc_fy = 0.0f;
    float recv_acc_fz = 0.0f;

    for (int j = 0; j < warp_size_int; j++) {
      // Compute interactions between send and recv atoms
      // Calculate distance between recv and send atoms
      const float dx = send_x - recv_x;
      const float dy = send_y - recv_y; 
      const float dz = send_z - recv_z;
      const float r2 = dx*dx + dy*dy + dz*dz;
      const float r = sqrtf(r2);

      // Only compute forces if both atoms are real and not excluded
      if (recv_topo_idx >= 0 && send_topo_idx >= 0 &&
          (! devcEvaluateLocalMask(recv_topo_idx, send_topo_idx, send_prof, lemr.aux_masks))) {

        // Coulomb force calculation
        const float inv_r = 1.0f / r;
        const float inv_r2 = inv_r*inv_r;
        float force = -(recv_q * send_q) * inv_r * inv_r2;

        // Get LJ parameters and calculate LJ forces
        int ljparm = (poly_nbk.n_lj_types[0] * recv_lj_idx + send_lj_idx) +
                     poly_nbk.ljabc_offsets[0];
        const float2 ljab = poly_nbk.ljab_coeff[ljparm];
        const float inv_r4 = inv_r2 * inv_r2;
        const float inv_r6 = inv_r4 * inv_r2;
        const float inv_r12 = inv_r6 * inv_r6;
        force += (12.0f * ljab.x * inv_r12 - 6.0f * ljab.y * inv_r6) * inv_r2;

        // Accumulate forces
        const float fx = force * dx;
        const float fy = force * dy;
        const float fz = force * dz;

        // Add contributions to sending atom force accumulators (in registers)
        send_acc_fx += fx;
        send_acc_fy += fy;
        send_acc_fz += fz;
        recv_acc_fx -= fx;
        recv_acc_fy -= fy;
        recv_acc_fz -= fz;

        // Shift the send atoms 1 to the left (lane 0 shifts to max lane (31))
        const int next_lane = ((lane_idx + 1) & warp_bits_mask_int);
        send_x = SHFL(send_x, next_lane);
        send_y = SHFL(send_y, next_lane);
        send_z = SHFL(send_z, next_lane);
        send_q = SHFL(send_q, next_lane);
        send_acc_fx = SHFL(send_acc_fx, next_lane);
        send_acc_fy = SHFL(send_acc_fy, next_lane);
        send_acc_fz = SHFL(send_acc_fz, next_lane);
        send_topo_idx = SHFL(send_topo_idx, next_lane);
        send_lj_idx = SHFL(send_lj_idx, next_lane);
      }
    }

    // Add the register accumulators to global arrays
    atomicSplit(send_acc_fx * poly_psw.frc_scale, instr.y + lane_idx, poly_psw.xfrc,
                poly_psw.xfrc_ovrf);
    atomicSplit(send_acc_fy * poly_psw.frc_scale, instr.y + lane_idx, poly_psw.yfrc,
                poly_psw.yfrc_ovrf);
    atomicSplit(send_acc_fz * poly_psw.frc_scale, instr.y + lane_idx, poly_psw.zfrc,
                poly_psw.zfrc_ovrf);
    atomicSplit(recv_acc_fx * poly_psw.frc_scale, instr.x + lane_idx, poly_psw.xfrc,
                poly_psw.xfrc_ovrf);
    atomicSplit(recv_acc_fy * poly_psw.frc_scale, instr.x + lane_idx, poly_psw.yfrc,
                poly_psw.yfrc_ovrf);
    atomicSplit(recv_acc_fz * poly_psw.frc_scale, instr.x + lane_idx, poly_psw.zfrc,
                poly_psw.zfrc_ovrf);
  }
}

//-------------------------------------------------------------------------------------------------
// C++ Function to compute interactions within a system
//-------------------------------------------------------------------------------------------------
void computeTileInteractions (PsSynthesisWriter *poly_psw, 
                              const SyNonbondedKit<float, float2> &poly_nbk,
                              const LocalExclusionMaskReader &lemr) {
  for(int i = 0; i < poly_psw->system_count; ++i) {
    const int atom_count = poly_psw->atom_counts[i];
    const int atom_start = poly_psw->atom_starts[i];
    const int ljabc_start = poly_nbk.ljabc_offsets[i];
    const int nlj_types = poly_nbk.n_lj_types[i];
    
    for(int j = atom_start; j < atom_start + atom_count; ++j){
      const double atom_jx = hostInt95ToDouble(poly_psw->xcrd[j], poly_psw->xcrd_ovrf[j]);
      const double atom_jy = hostInt95ToDouble(poly_psw->ycrd[j], poly_psw->ycrd_ovrf[j]);
      const double atom_jz = hostInt95ToDouble(poly_psw->zcrd[j], poly_psw->zcrd_ovrf[j]);
      const double atom_jq = poly_nbk.charge[j] * poly_nbk.coulomb;

      const int ljidx_j = poly_nbk.lj_idx[j];

      for(int k = atom_start; k < j; ++k){
        double fmag, dx, dy, dz;
        if(testExclusion(lemr, j, k)) {
          fmag = 0.0;
          dx = 0.0;
          dy = 0.0;
          dz = 0.0;
        }
        else {
          const double atom_kx = hostInt95ToDouble(poly_psw->xcrd[k], poly_psw->xcrd_ovrf[k]);
          const double atom_ky = hostInt95ToDouble(poly_psw->ycrd[k], poly_psw->ycrd_ovrf[k]);
          const double atom_kz = hostInt95ToDouble(poly_psw->zcrd[k], poly_psw->zcrd_ovrf[k]);
          const double atom_kq = poly_nbk.charge[k] * poly_nbk.coulomb;

          const int ljidx_k = poly_nbk.lj_idx[k];
          const int ljparm_jk = nlj_types * ljidx_j + ljidx_k + ljabc_start;
          const float2 ljab_jk = poly_nbk.ljab_coeff[ljparm_jk];
          dx = atom_kx - atom_jx;
          dy = atom_ky - atom_jy;
          dz = atom_kz - atom_jz;
          double distance_sq = dx*dx + dy*dy + dz*dz;
          double distance = sqrt(distance_sq);
          double fmag = -(atom_jq * atom_kq) / (distance * distance_sq);
          const double inv_distance = 1.0 / distance;
          const double inv_distance_sq = inv_distance * inv_distance;
          const double inv_distance_4 = inv_distance_sq * inv_distance_sq;
          const double inv_distance_6 = inv_distance_sq * inv_distance_4;
          const double inv_distance_8 = inv_distance_4  * inv_distance_4;
          fmag += ((6.0 * ljab_jk.y) - (12.0 * ljab_jk.x * inv_distance_6)) * inv_distance_8;
        }
        fmag *= poly_psw->frc_scale;
        const double force_x = (fmag * dx);
        const double force_y = (fmag * dy);
        const double force_z = (fmag * dz);
        
        // Update the forces on atom j and atom k
        const int95_t ifx = hostDoubleToInt95(force_x);
        const int95_t ify = hostDoubleToInt95(force_y);
        const int95_t ifz = hostDoubleToInt95(force_z);
        hostSplitFPSum(&poly_psw->xfrc[j], &poly_psw->xfrc_ovrf[j], ifx);
        hostSplitFPSum(&poly_psw->yfrc[j], &poly_psw->yfrc_ovrf[j], ify);
        hostSplitFPSum(&poly_psw->zfrc[j], &poly_psw->zfrc_ovrf[j], ifz);
        hostSplitFPSubtract(&poly_psw->xfrc[k], &poly_psw->xfrc_ovrf[k], ifx);
        hostSplitFPSubtract(&poly_psw->yfrc[k], &poly_psw->yfrc_ovrf[k], ify);
        hostSplitFPSubtract(&poly_psw->zfrc[k], &poly_psw->zfrc_ovrf[k], ifz);
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Main function
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {
  
  // Engage the testing environment and the GPU
  StopWatch timer;
  HpcConfig gpu_config(ExceptionResponse::WARN);
  std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> engage_gpu(1);
  
  // Time categories
  const int time_load  = timer.addCategory("Load input files");
  const int time_build = timer.addCategory("Class object setup");
  const int time_instr_load = timer.addCategory("Setting up tile interactions");
  const int time_cpp_tiling = timer.addCategory("C++ Tiling Interactions");
  const int time_gpu_basic_tiling = timer.addCategory("GPU Tiling Interactions (Basic)");
  const int time_gpu_alt_tiling = timer.addCategory("GPU Tiling (Alt. sending and receiving)");
  
  // Program Specific Command Line Inputs
  CommandLineParser clip("gpu_tiling", "A benchmarking program to test GPU tile strategies.",
                         {"-timings"});
  clip.addStandardAmberInputs("-p", "-c", "-ig_seed");
  clip.addStandardBenchmarkingInputs("-fp_bits", "-iter");
  NamelistEmulator *t_nml = clip.getNamelistPointer();

  // Initialize the testing environment
  TestEnvironment oe(argc, argv, &clip, TmpdirStatus::NOT_REQUIRED, ExceptionResponse::SILENT);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Read command line instructions
  clip.parseUserInput(argc, argv);
  const std::string inpcrd_file = t_nml->getStringValue("-c");
  const std::string topology_file = t_nml->getStringValue("-p");
  const int ig_seed = t_nml->getIntValue("-ig_seed");
  const int fp_bits = t_nml->getIntValue("-fp_bits");
  const int iterations = t_nml->getIntValue("-iter");
  
  // Check that a system has been provided and load the basic topology and input coordinates
  TestSystemManager tsm(std::vector<std::string>(1, topology_file),
                        std::vector<std::string>(1, inpcrd_file), ExceptionResponse::DIE);
  timer.assignTime(time_load);

  // Create basic objects and load system data
  PhaseSpaceSynthesis poly_ps = tsm.exportPhaseSpaceSynthesis(std::vector<int>(1, 0),
                                                              0.0, ig_seed, 40, 44, fp_bits);

  int total_tile_count = 0;
  for(int i = 0; i < poly_ps.getSystemCount(); ++i){
    PhaseSpace ps = poly_ps.exportSystem(i);
    PhaseSpaceWriter psw = ps.data();
    if (psw.unit_cell != UnitCellType::NONE) {
      imageCoordinates<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.natom, psw.umat, psw.invu,
                                       psw.unit_cell, ImagingMethod::PRIMARY_UNIT_CELL);
    }
    poly_ps.importSystem(ps, i);
    int atom_counts = psw.natom;
    int tiles = (atom_counts + 31) / 32;
    total_tile_count += (tiles * (tiles + 1)) / 2;
  }

  // x: Atom Start for the receiving atoms
  // y: Atom Start for the sending atoms
  // z: Number of receiving atoms
  // w: Number of sending atoms
  Hybrid<int4> tile_assignments(total_tile_count, "tile_assignments");

  // Stage initial data and tile interaction instructions
  int counter = 0;
  PsSynthesisWriter poly_psw = poly_ps.data();
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {

    // Range over all tiles in the current system
    int atom_start = poly_psw.atom_starts[i];
    int atom_count = poly_psw.atom_counts[i];
    int tiles = (atom_count + 31) / 32;

    // Assign tile assignments for each tile in the system
    for (int j = 0; j < tiles; j++) {
      int start_recv = atom_start + (j * 32);  // Starting index for receiving atoms
      int end_recv = std::min(atom_start + (j + 1) * 32, atom_start + atom_count);  
      int num_recv = end_recv - start_recv;  // Number of receiving atoms

      for (int k = j; k < tiles; ++k) {
        int start_send = atom_start + (k * 32);  // Starting index for sending atoms
        int end_send = std::min(atom_start + (k + 1) * 32, atom_start + atom_count);  
        int num_send = end_send - start_send;  // Number of sending atoms

        // Fill the tile assignments
        tile_assignments.putHost({ start_recv, start_send, num_recv, num_send }, counter);
        counter++;
      }
    }
  }
  timer.assignTime(time_build);
  PsSynthesisWriter host_psw = poly_ps.data();
  Xoshiro256ppGenerator xrs(ig_seed);
  for (int i = 1; i < host_psw.system_count; i++) {
    const size_t aoff = host_psw.atom_starts[i];
    addRandomNoise(&xrs, &host_psw.xcrd[aoff], &host_psw.xcrd_ovrf[aoff], &host_psw.ycrd[aoff],
                   &host_psw.ycrd_ovrf[aoff], &host_psw.zcrd[aoff], &host_psw.zcrd_ovrf[aoff],
                   host_psw.atom_counts[i], 0.01, host_psw.gpos_scale);
  }
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(std::vector<int>(1, 0));
  
  LocalExclusionMask lem(poly_ag);
  SyNonbondedKit<float, float2> host_nbk = poly_ag.getSinglePrecisionNonbondedKit();

  // Run through the C++ function and collect data first (Ground Truth)
  // Create copies for C++ ground truth comparison
  computeTileInteractions(&host_psw, host_nbk, lem.data());

  // Upload basic objects to the GPU
  poly_ps.upload();
  poly_ag.upload();
  tile_assignments.upload();
  lem.upload();
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  PsSynthesisWriter devc_psw = poly_ps.data(devc);
  SyNonbondedKit<float, float2> poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(devc);
  const int4 *t_assgn = tile_assignments.data(devc);
  for (int i = 0; i < iterations; ++i) {
    poly_ps.initializeForces(gpu, devc);
    computeTileInteractions<<<5 * gpu.getSMPCount(), 256>>>(devc_psw, poly_nbk, t_assgn, 
                                                            total_tile_count, lem.data());
  }
  cudaDeviceSynchronize();
  timer.assignTime(time_gpu_basic_tiling);
  poly_ps.download();
  poly_ag.download();
  CoordinateFrame cpu_frc(poly_ps.getAtomCount(0), poly_ps.getUnitCellType(),
                          HybridFormat::HOST_MOUNTED);
  CoordinateFrame gpu_frc(poly_ps.getAtomCount(0), poly_ps.getUnitCellType(),
                          HybridFormat::HOST_MOUNTED);
  const CoordinateFrameReader cpu_cfr = cpu_frc.data();
  const CoordinateFrameReader gpu_cfr = gpu_frc.data();
  std::vector<double> cpu_chk(cpu_cfr.natom, 0.0);
  std::vector<double> gpu_chk(cpu_cfr.natom, 0.0);

  coordCopy(&cpu_frc, poly_ps, 0, TrajectoryKind::FORCES, HybridTargetLevel::HOST, 
            HybridTargetLevel::HOST);
  coordCopy(&gpu_frc, poly_ps, 0, TrajectoryKind::FORCES, HybridTargetLevel::HOST,
            HybridTargetLevel::DEVICE);
  
  for (int j = 0; j < cpu_cfr.natom; j++) {
    cpu_chk[j] = cpu_cfr.xcrd[j];
    gpu_chk[j] = gpu_cfr.xcrd[j]; 
  }

  check(gpu_chk, RelationalOperator::EQUAL, Approx(cpu_chk).margin(1.0e-2),
        "Forces on particles calculated using the GPU kernel running in " + 
        getEnumerationName(PrecisionModel::SINGLE) + "-precision do not agree with CPU forces " +
        "calculated in " + getEnumerationName(PrecisionModel::DOUBLE) + "-precision.");

  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
}
