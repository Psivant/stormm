// -*-c++-*-
#include <cuda.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <nvml.h>
#include <string>
#include <vector>
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/core_kernel_manager.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/Constants/scaling.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/MolecularMechanics/mm_controls.h"
#include "../../src/Namelists/command_line_parser.h"
#include "../../src/Namelists/nml_files.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Potential/cacheresource.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/hpc_valence_potential.h"
#include "../../src/Potential/hpc_nonbonded_potential.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Reporting/error_format.h"
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
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/stopwatch.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::diskutil;
using namespace stormm::errors;
using namespace stormm::energy;
using namespace stormm::stmath;
using namespace stormm::mm;
using namespace stormm::namelist;
using namespace stormm::numerics;
using namespace stormm::parse;
using namespace stormm::restraints;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Constant expressions to guide testing
//-------------------------------------------------------------------------------------------------
constexpr int null_kernel_repeats = 2500;
constexpr int args_kernel_repeats = 2500;
constexpr int gmem_kernel_repeats = 2500;

//-------------------------------------------------------------------------------------------------
// Get a SystemCache object containing all topologies and coordinates in a pair of directories.
//
// Arguments:
//   topol_path:  A series of strings that will be joined into the topology directory name
//   coord_path:  A series of strings that will be joined into the coordinate directory name
//   oe:          Contains critical shell variables such as the $STORMM source path where the
//                named directories are expected to reside
//-------------------------------------------------------------------------------------------------
SystemCache directorySweep(const std::vector<std::string> &topol_path,
                           const std::vector<std::string> &coord_path, const TestEnvironment &oe) {
  
  // Collect coordinates and topologies
  const char osc = osSeparator();
  std::string buffer("&files\n  -p ");
  buffer += oe.getStormmSourcePath() + osc + "benchmark";
  for (size_t i = 0; i < topol_path.size(); i++) {
    buffer += osc + topol_path[i];
  }
  buffer += "\n  -c ";
  buffer += oe.getStormmSourcePath() + osc + "benchmark";
  for (size_t i = 0; i < topol_path.size(); i++) {
    buffer += osc + coord_path[i];
  }
  buffer += "\n&end\n";
  const TextFile tf(buffer, TextOrigin::RAM);
  int start_line = 0;
  FilesControls fcon(tf, &start_line);
  return SystemCache(fcon, ExceptionResponse::SILENT, MapRotatableGroups::NO);
}

//-------------------------------------------------------------------------------------------------
// Replicate a single topology and coordinate system many times, then run kernels to obtain
// timings.
//
// Arguments:
//   ag_vec:      The topologies to replicate
//   ps_vec:      The coordinates to replicate
//   nrep:        The number of replicas to make
//   mmctrl:      Step counter and progress counters for all work units
//   gpu:         Details of the GPU to use in calculations
//   timer:       Object to record the timings
//   prec:        The precision model in which to perform calculations
//   eval_frc:    Flag to indicate that forces should be evaluated
//   eval_nrg:    Flag to indicate that energies should be evaluated
//   acc_meth:    The accumulation method to use in each kernel
//   purpose:     Indicate whether to accumulate forces or (instead) to move atoms
//   batch_name:  Name of the batch, for reporting purposes
//-------------------------------------------------------------------------------------------------
void replicaProcessing(const std::vector<AtomGraph*> &ag_vec,
                       const std::vector<PhaseSpace> &ps_vec, const int nrep,
                       MolecularMechanicsControls *mmctrl, const GpuDetails &gpu,
                       StopWatch *timer, const PrecisionModel prec, const EvaluateForce eval_frc,
                       const EvaluateEnergy eval_nrg, const AccumulationMethod acc_meth,
                       const VwuGoal purpose, const std::string &batch_name,
                       const int iter = 100) {
  std::vector<int> ag_idx = tileVector(incrementingSeries<int>(0, ag_vec.size()), nrep);
  AtomGraphSynthesis poly_ag(ag_vec, ag_idx, ExceptionResponse::SILENT, gpu);
  PhaseSpaceSynthesis poly_ps(ps_vec, ag_vec, ag_idx);
  StaticExclusionMaskSynthesis poly_se;
  switch (poly_ag.getUnitCellType()) {
  case UnitCellType::NONE:
    poly_se = StaticExclusionMaskSynthesis(ag_vec, ag_idx);
    poly_ag.loadNonbondedWorkUnits(poly_se);
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    break;
  }
  SeMaskSynthesisReader poly_ser = poly_se.data();
  CoreKlManager launcher(gpu, poly_ag);
  ScoreCard sc(poly_ag.getSystemCount(), 1, 32);
  const int2 valence_lp = launcher.getValenceKernelDims(prec, eval_frc, eval_nrg, acc_meth,
                                                        purpose, ClashResponse::NONE);
  int2 nonbond_lp, gbr_lp, gbd_lp;
  switch (poly_ag.getUnitCellType()) {
  case UnitCellType::NONE:
    nonbond_lp = launcher.getNonbondedKernelDims(prec, NbwuKind::TILE_GROUPS, eval_frc, eval_nrg,
                                                 acc_meth, ImplicitSolventModel::NONE,
                                                 ClashResponse::NONE);
    gbr_lp = launcher.getBornRadiiKernelDims(prec, NbwuKind::TILE_GROUPS, acc_meth,
                                             ImplicitSolventModel::NONE);
    gbd_lp = launcher.getBornDerivativeKernelDims(prec, NbwuKind::TILE_GROUPS, acc_meth,
                                                  ImplicitSolventModel::NONE);
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    nonbond_lp = { 1, 1 };
    break;
  }
  CacheResource valence_tb_space(valence_lp.x, maximum_valence_work_unit_atoms);
  CacheResource nonbond_tb_space(nonbond_lp.x, small_block_max_atoms);
  mmctrl->primeWorkUnitCounters(launcher, eval_frc, eval_nrg, ClashResponse::NONE,
                                VwuGoal::ACCUMULATE, prec, prec, poly_ag);
  Thermostat tstat(poly_ag, ThermostatKind::NONE);
  ImplicitSolventWorkspace isw(poly_ag.getSystemAtomOffsets(), poly_ag.getSystemAtomCounts(),
                               prec);

  // Upload the critical components
  poly_ag.upload();
  poly_se.upload();
  poly_ps.upload();

  // Some common variables for either branch
  const ValenceKernelSize kwidth = poly_ag.getValenceThreadBlockSize();
  const std::string valk_name = valenceKernelKey(prec, eval_frc, eval_nrg, acc_meth, purpose,
                                                 ClashResponse::NONE, kwidth);
  std::string nnbk_name;
  switch (poly_ag.getUnitCellType()) {
  case UnitCellType::NONE:
    nnbk_name = nonbondedKernelKey(prec, NbwuKind::TILE_GROUPS, eval_frc, eval_nrg, acc_meth,
                                   ImplicitSolventModel::NONE, ClashResponse::NONE);
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    break;
  }
  const int sys_val_time = timer->addCategory(batch_name + " on " + valk_name + " (" +
                                              std::to_string(nrep) + ")");
  const int sys_nb_time = timer->addCategory(batch_name + " on " + nnbk_name + " (" +
                                             std::to_string(nrep) + ")");
  
  // Obtain abstracts outside the inner loop, in case this is a significant contributor to the
  // run time.  Forces will only be initialized once, and thereafter calculated repeatedly to test
  // only the run time of the one kernel.
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(devc_tier);
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag.getDoublePrecisionNonbondedKit(devc_tier);
      const SyRestraintKit<double,
                           double2,
                           double4> poly_rk = poly_ag.getDoublePrecisionRestraintKit(devc_tier);
      const SyAtomUpdateKit<double,
                            double2,
                            double4> poly_auk = poly_ag.getDoublePrecisionAtomUpdateKit(devc_tier);
      MMControlKit<double> ctrl = mmctrl->dpData(devc_tier);
      ThermostatWriter<double> tstw = tstat.dpData(devc_tier);
      ISWorkspaceKit<double> iswk = isw.dpData(devc_tier);
      PsSynthesisWriter poly_psw = poly_ps.data(devc_tier);
      ScoreCardWriter scw = sc.data(devc_tier);
      CacheResourceKit<double> gmem_rval = valence_tb_space.dpData(devc_tier);
      CacheResourceKit<double> gmem_rnnb = nonbond_tb_space.dpData(devc_tier);
      poly_ps.initializeForces(gpu, devc_tier);
      timer->assignTime(0);

      // Test the valence kernel
      for (int i = 0; i < iter; i++) {
        ctrl.step += 1;
        launchValence(poly_vk, poly_rk, &ctrl, &poly_psw, poly_auk, &tstw, &scw, &gmem_rval,
                      eval_frc, eval_nrg, purpose, valence_lp);
      }
      cudaDeviceSynchronize();
      timer->assignTime(sys_val_time);

      // Test the non-bonded kernels
      switch (poly_ag.getUnitCellType()) {
      case UnitCellType::NONE:
        poly_ps.initializeForces(gpu, devc_tier);
        timer->assignTime(0);
        for (int i = 0; i < iter; i++) {
          ctrl.step += 1;
          launchNonbonded(NbwuKind::TILE_GROUPS, poly_nbk, poly_ser, &ctrl, &poly_psw, &tstw, &scw,
                          &gmem_rnnb, &iswk, eval_frc, eval_nrg, nonbond_lp, gbr_lp, gbd_lp);
        }
        cudaDeviceSynchronize();
        timer->assignTime(sys_nb_time);
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        break;
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(devc_tier);
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(devc_tier);
      const SyRestraintKit<float,
                           float2,
                           float4> poly_rk = poly_ag.getSinglePrecisionRestraintKit(devc_tier);
      const SyAtomUpdateKit<float,
                            float2,
                            float4> poly_auk = poly_ag.getSinglePrecisionAtomUpdateKit(devc_tier);
      MMControlKit<float> ctrl = mmctrl->spData(devc_tier);
      ThermostatWriter<float> tstw = tstat.spData(devc_tier);
      ISWorkspaceKit<float> iswk = isw.spData(devc_tier);
      PsSynthesisWriter poly_psw = poly_ps.data(devc_tier);
      ScoreCardWriter scw = sc.data(devc_tier);
      CacheResourceKit<float> gmem_rval = valence_tb_space.spData(devc_tier);
      CacheResourceKit<float> gmem_rnnb = nonbond_tb_space.spData(devc_tier);
      poly_ps.initializeForces(gpu, devc_tier);      
      timer->assignTime(0);

      // Test the valence kernel
      for (int i = 0; i < iter; i++) {
        ctrl.step += 1;
        launchValence(poly_vk, poly_rk, &ctrl, &poly_psw, poly_auk, &tstw, &scw, &gmem_rval,
                      eval_frc, eval_nrg, purpose, acc_meth, valence_lp);
      }
      cudaDeviceSynchronize();
      timer->assignTime(sys_val_time);
      
      // Test the non-bonded kernels
      switch (poly_ag.getUnitCellType()) {
      case UnitCellType::NONE:
        poly_ps.initializeForces(gpu, devc_tier);      
        timer->assignTime(0);
        for (int i = 0; i < iter; i++) {
          ctrl.step += 1;
          launchNonbonded(NbwuKind::TILE_GROUPS, poly_nbk, poly_ser, &ctrl, &poly_psw, &tstw, &scw,
                          &gmem_rnnb, &iswk, eval_frc, eval_nrg, acc_meth, nonbond_lp, gbr_lp,
                          gbd_lp);
        }
        cudaDeviceSynchronize();
        timer->assignTime(sys_nb_time);
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        break;
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
// Run a batch of molecules based on a particular folder within the Topologies and Coordinates
// subdirectories of the STORMM benchmarking suite.
//
// Arguments:
//   batch_name:  Name of the folder
//   gpu:         Details of the GPU
//   oe:          Contains the name of the STORMM source directory, and other environment variables
//   timer:       Object to track the timings
//   iter:        The number of iterations for which to run each kernel
//-------------------------------------------------------------------------------------------------
void runBatch(const std::string &batch_name, const GpuDetails &gpu, const TestEnvironment &oe,
              StopWatch *timer, const int iter = 100) {
  
  // Read the molecules in this batch
  const std::vector<std::string> topols = { "Topologies", batch_name, ".*_ff1.*SB.top" };
  const std::vector<std::string> coords = { "Coordinates", batch_name, ".*_ff1.*SB.inpcrd"};
  SystemCache sc = directorySweep(topols, coords, oe);
  MolecularMechanicsControls mmctrl;

  // Loop over the dipeptides one at a time, make syntheses of each of them individually, and
  // test kernels.
  const int mol_count = sc.getSystemCount();
  const std::vector<int> batch_multiplier = { 1, 2, 3, 4, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50,
                                              55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 125,
                                              150, 175, 200, 225, 250, 275, 300, 325, 350, 375,
                                              400, 425, 450, 475, 500, 625, 750, 875, 1000, 1125,
                                              1250, 1375, 1500, 1625, 1750, 1875, 2000, 2250, 2500,
                                              2750, 3000, 3250, 3500, 3750, 4000, 4250, 4500, 4750,
                                              5000 };
  std::vector<AtomGraph*> ag_vec(mol_count);
  std::vector<PhaseSpace> ps_vec;
  ps_vec.reserve(mol_count);
  for (int i = 0; i < mol_count; i++) {
    ag_vec[i] = sc.getSystemTopologyPointer(i);
    ps_vec.push_back(sc.getCoordinates(i));
  }
  for (size_t i = 0; i < batch_multiplier.size(); i++) {
    const int ncopy = batch_multiplier[i];
    replicaProcessing(ag_vec, ps_vec, ncopy, &mmctrl, gpu, timer, PrecisionModel::SINGLE,
                      EvaluateForce::YES, EvaluateEnergy::NO, AccumulationMethod::SPLIT,
                      VwuGoal::ACCUMULATE, batch_name, iter);
    replicaProcessing(ag_vec, ps_vec, ncopy, &mmctrl, gpu, timer, PrecisionModel::SINGLE,
                      EvaluateForce::YES, EvaluateEnergy::YES, AccumulationMethod::SPLIT,
                      VwuGoal::ACCUMULATE, batch_name, iter);
    replicaProcessing(ag_vec, ps_vec, ncopy, &mmctrl, gpu, timer, PrecisionModel::SINGLE,
                      EvaluateForce::YES, EvaluateEnergy::NO, AccumulationMethod::WHOLE,
                      VwuGoal::ACCUMULATE, batch_name, iter);
    replicaProcessing(ag_vec, ps_vec, ncopy, &mmctrl, gpu, timer, PrecisionModel::SINGLE,
                      EvaluateForce::YES, EvaluateEnergy::YES, AccumulationMethod::WHOLE,
                      VwuGoal::ACCUMULATE, batch_name, iter);

    // Only do double-precision calculations for low replica numbers--this can be strenuous on
    // many architectures, particularly in the non-bonded kernel.
    if (batch_multiplier[i] < 3) {
      replicaProcessing(ag_vec, ps_vec, ncopy, &mmctrl, gpu, timer, PrecisionModel::DOUBLE,
                        EvaluateForce::YES, EvaluateEnergy::NO, AccumulationMethod::SPLIT,
                        VwuGoal::ACCUMULATE, batch_name, iter);
      replicaProcessing(ag_vec, ps_vec, ncopy, &mmctrl, gpu, timer, PrecisionModel::DOUBLE,
                        EvaluateForce::YES, EvaluateEnergy::YES, AccumulationMethod::SPLIT,
                        VwuGoal::ACCUMULATE, batch_name, iter);
    }
  }

  // Test many replicas of each small molecule system
  const int nimg = gpu.getSMPCount() * 64;
  for (size_t i = 0; i < mol_count; i++) {
    std::vector<AtomGraph*> ss_agv(1, sc.getSystemTopologyPointer(i));
    std::vector<PhaseSpace> ss_psv(1, sc.getCoordinates(i));
    const std::string ss_name = getBaseName(ss_agv[0]->getFileName()) + " " +
                                std::to_string(ss_agv[0]->getAtomCount());
    replicaProcessing(ss_agv, ss_psv, nimg, &mmctrl, gpu, timer, PrecisionModel::SINGLE,
                      EvaluateForce::YES, EvaluateEnergy::NO, AccumulationMethod::SPLIT,
                      VwuGoal::ACCUMULATE, ss_name, iter);
    replicaProcessing(ss_agv, ss_psv, nimg, &mmctrl, gpu, timer, PrecisionModel::SINGLE,
                      EvaluateForce::YES, EvaluateEnergy::YES, AccumulationMethod::SPLIT,
                      VwuGoal::ACCUMULATE, ss_name, iter);
    replicaProcessing(ss_agv, ss_psv, nimg, &mmctrl, gpu, timer, PrecisionModel::SINGLE,
                      EvaluateForce::YES, EvaluateEnergy::NO, AccumulationMethod::WHOLE,
                      VwuGoal::ACCUMULATE, ss_name, iter);
    replicaProcessing(ss_agv, ss_psv, nimg, &mmctrl, gpu, timer, PrecisionModel::SINGLE,
                      EvaluateForce::YES, EvaluateEnergy::YES, AccumulationMethod::WHOLE,
                      VwuGoal::ACCUMULATE, ss_name, iter);
  }
}

//-------------------------------------------------------------------------------------------------
// A kernel that does nothing.  Points its finger but there's no one around.
//-------------------------------------------------------------------------------------------------
__global__ void kNothing() {
}

//-------------------------------------------------------------------------------------------------
// Test a null kernel, launching blocks of a specified thread count and grid size.
//
// Arguments:
//   block_size:           Thread count of the blocks to launch
//   block_count_per_smp:  The number of blocks to launch per streaming multiprocessor
//   gpu:                  Details of the GPU to use in calculations
//   timer:                Object to collect timings data
//-------------------------------------------------------------------------------------------------
void testNullKernel(const int block_size, const int block_count_per_smp, const GpuDetails &gpu,
                    StopWatch *timer) {
  const int nblocks = block_count_per_smp * gpu.getSMPCount();
  const int timings_sect = timer->addCategory("Null Kernel Launch, " + std::to_string(block_size) +
                                              " x " + std::to_string(block_count_per_smp));
  timer->assignTime(0);
  for (int trial = 0; trial < 8; trial++) {
    for (int i = 0; i < null_kernel_repeats; i++) {
      kNothing<<<nblocks, block_size>>>();
    }
    cudaDeviceSynchronize();
    timer->assignTime(timings_sect);
  }
}

//-------------------------------------------------------------------------------------------------
// A kernel that takes arguments but still does nothing.
//-------------------------------------------------------------------------------------------------
__global__ void kTakeArguments(const ValenceKit<float> vk, const NonbondedKit<float> nbk,
                               const VirtualSiteKit<float> vsk, const ConstraintKit<float> cnk) {
}

//-------------------------------------------------------------------------------------------------
// Test a kernel that accepts arguments (a little more than 1kB worth, if the sizes of various
// topology abstracts do not change significantly).
//
// Arguments:
//   block_size:           Thread count of the blocks to launch
//   block_count_per_smp:  The number of blocks to launch per streaming multiprocessor
//   gpu:                  Details of the GPU to use in calculations
//   timer:                Object to collect timings data
//-------------------------------------------------------------------------------------------------
void testArgumentLoadedKernel(const int block_size, const int block_count_per_smp,
                              const GpuDetails &gpu, StopWatch *timer) {
  const int nblocks = block_count_per_smp * gpu.getSMPCount();
  const int timings_sect = timer->addCategory("Args Kernel Launch, " + std::to_string(block_size) +
                                              " x " + std::to_string(block_count_per_smp));
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  AtomGraph ag;
  ValenceKit<float> ag_vk = ag.getSinglePrecisionValenceKit(devc_tier);
  NonbondedKit<float> ag_nbk = ag.getSinglePrecisionNonbondedKit(devc_tier);
  VirtualSiteKit<float> ag_vsk = ag.getSinglePrecisionVirtualSiteKit(devc_tier);
  ConstraintKit<float> ag_cnk = ag.getSinglePrecisionConstraintKit(devc_tier);
  timer->assignTime(0);
  for (int trial = 0; trial < 8; trial++) {
    for (int i = 0; i < args_kernel_repeats; i++) {
      kTakeArguments<<<nblocks, block_size>>>(ag_vk, ag_nbk, ag_vsk, ag_cnk);
    }
    cudaDeviceSynchronize();
    timer->assignTime(timings_sect);
  }
}

//-------------------------------------------------------------------------------------------------
// A kernel that runs a lap through global memory.
//-------------------------------------------------------------------------------------------------
__global__ void kGmemAccess(int* read_from, int* write_to, int threads_active) {
  if (threadIdx.x < threads_active) {
    const int pos = (blockIdx.x * blockDim.x) + threadIdx.x;
    const int var = read_from[pos];
    write_to[pos] = var * 2;
  }
}

//-------------------------------------------------------------------------------------------------
// Test a kernel that takes no arguments and merely reads, then writes elements of global memory.
// Both reads and writes are cold operations on non-cached addresses.
//
// Arguments:
//   block_size:           Thread count of the blocks to launch
//   block_count_per_smp:  The number of blocks to launch per streaming multiprocessor
//   gpu:                  Details of the GPU to use in calculations
//   timer:                Object to collect timings data
//-------------------------------------------------------------------------------------------------
void testGmemKernel(const int block_size, const int block_count_per_smp, const GpuDetails &gpu,
                    StopWatch *timer) {
  const int nthreads = gpu.getSMPCount() * block_size * block_count_per_smp;
  const int nblocks = block_count_per_smp * gpu.getSMPCount();
  Hybrid<int> read_from(nthreads);
  Hybrid<int> write_to(nthreads);
  int* read_ptr = read_from.data();
  int neg_fac = 1;
  for (int i = 0; i < nthreads; i++) {
    read_ptr[i] = neg_fac * i;
    neg_fac *= -1;
  }
  read_from.upload();
  int* read_devc_ptr = read_from.data(HybridTargetLevel::DEVICE);
  int* write_devc_ptr = write_to.data(HybridTargetLevel::DEVICE);
  const int none_timings = timer->addCategory("GMEM Kernel launch, zero warps active");
  const int half_timings = timer->addCategory("GMEM Kernel launch, half warps active");
  const int full_timings = timer->addCategory("GMEM Kernel launch, full warps active");
  timer->assignTime(0);
  for (int trial = 0; trial < 8; trial++) {
    for (int i = 0; i < gmem_kernel_repeats; i++) {
      kGmemAccess<<<nblocks, block_size>>>(read_devc_ptr, write_devc_ptr, 0);
    }
    cudaDeviceSynchronize();
    timer->assignTime(none_timings);
    for (int i = 0; i < gmem_kernel_repeats; i++) {
      kGmemAccess<<<nblocks, block_size>>>(read_devc_ptr, write_devc_ptr, block_size / 2);
    }
    cudaDeviceSynchronize();
    timer->assignTime(half_timings);
    for (int i = 0; i < gmem_kernel_repeats; i++) {
      kGmemAccess<<<nblocks, block_size>>>(read_devc_ptr, write_devc_ptr, block_size);
    }
    cudaDeviceSynchronize();
    timer->assignTime(full_timings);
  }
}

//-------------------------------------------------------------------------------------------------
// Test kernel timings on a command-line specified topology and coordinate set.
//
// Arguments:
//   top_name:  Name of the topology (if this is actually a coordinate file, it will be taken as
//              such)
//   crd_name:  Name of the coordinate file (if this is actually a topology, it will be taken as
//              such)
//   iter:      The number of iterations for which to run various kernels
//   replicas:  The number of times to repliate the system within a synthesis used in computations
//-------------------------------------------------------------------------------------------------
void runAdditionalTest(const std::string &top_name, const std::string &crd_name,
                       const GpuDetails &gpu, StopWatch *timer, const int iter = 100,
                       const int replicas = 100) {

  // Take in the one system and check its properties
  AtomGraph ag;
  std::vector<PhaseSpace> psv;
  try {
    ag = AtomGraph(top_name);
    psv.emplace_back(crd_name);
  }
  catch (std::runtime_error) {
    try {
      ag = AtomGraph(crd_name);
      psv.emplace_back(top_name);
    }
    catch (std::runtime_error) {
      rtErr("The files " + top_name + " and " + crd_name + " do not correspond to a valid pair of "
            "topology and coordinate files.", "runPeriodicTest");
    }
  }
  if (ag.getAtomCount() != psv[0].getAtomCount()) {
    rtErr("The topology and coordinate set have different numbers of atoms (" +
          std::to_string(ag.getAtomCount()) + " in the topology, " +
          std::to_string(psv[0].getAtomCount()) + " in the coordinate set).", "runPeriodicTest");
  }
  std::vector<AtomGraph*> agv(1, &ag);
  MolecularMechanicsControls mmctrl;
  const std::string test_name = getBaseName(ag.getFileName());
  replicaProcessing(agv, psv, replicas, &mmctrl, gpu, timer, PrecisionModel::SINGLE,
                    EvaluateForce::YES, EvaluateEnergy::NO, AccumulationMethod::SPLIT,
                    VwuGoal::ACCUMULATE, test_name, iter);
  replicaProcessing(agv, psv, replicas, &mmctrl, gpu, timer, PrecisionModel::SINGLE,
                    EvaluateForce::YES, EvaluateEnergy::YES, AccumulationMethod::SPLIT,
                    VwuGoal::ACCUMULATE, test_name, iter);
  replicaProcessing(agv, psv, replicas, &mmctrl, gpu, timer, PrecisionModel::SINGLE,
                    EvaluateForce::YES, EvaluateEnergy::NO, AccumulationMethod::WHOLE,
                    VwuGoal::ACCUMULATE, test_name, iter);
  replicaProcessing(agv, psv, replicas, &mmctrl, gpu, timer, PrecisionModel::SINGLE,
                    EvaluateForce::YES, EvaluateEnergy::YES, AccumulationMethod::WHOLE,
                    VwuGoal::ACCUMULATE, test_name, iter);
  replicaProcessing(agv, psv, replicas, &mmctrl, gpu, timer, PrecisionModel::DOUBLE,
                    EvaluateForce::YES, EvaluateEnergy::NO, AccumulationMethod::SPLIT,
                    VwuGoal::ACCUMULATE, test_name, iter);
  replicaProcessing(agv, psv, replicas, &mmctrl, gpu, timer, PrecisionModel::DOUBLE,
                    EvaluateForce::YES, EvaluateEnergy::YES, AccumulationMethod::SPLIT,
                    VwuGoal::ACCUMULATE, test_name, iter);
}
//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  StopWatch timer;
  HpcConfig gpu_config(ExceptionResponse::WARN);
  std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> engage_gpu(1);

  // Prepare program-specific command line inputs
  CommandLineParser clip("test_nonperiodic_kernels", "A program for timing general kernel "
                         "performance in non-periodic systems.", { "-timings" });
  clip.addStandardAmberInputs("-p", "-c");
  clip.addStandardBenchmarkingInputs("-iter", "-replicas");
  NamelistEmulator *t_nml = clip.getNamelistPointer();
  t_nml->addKeyword("-klaunch", NamelistType::BOOLEAN);
  t_nml->addHelp("-klaunch", "Test the baseline kernel launch latency on the GPU.");

  // Initialize a test environment, for wall time tracking and possible testing.
  TestEnvironment oe(argc, argv, &clip, TmpdirStatus::NOT_REQUIRED, ExceptionResponse::SILENT);

  // Parse command-line information
  clip.parseUserInput(argc, argv);
  const int iter = t_nml->getIntValue("-iter");
  const int replicas = t_nml->getIntValue("-replicas");
  const bool test_kernel_launch = t_nml->getBoolValue("-klaunch");
  const std::string extra_top = t_nml->getStringValue("-p");
  const std::string extra_crd = t_nml->getStringValue("-c");
  
  // Check user input
  if (iter < 0 || iter > 100000) {
    rtErr("The number of kernel iterations must be a positive integer between 1 and 100000.",
          "main");
  }
  if (replicas <= 0 || replicas >= 65536) {
    rtErr("The number of system replicas must be a positive integer less than 65536.", "main");
  }
  
  // Test some basic kernels to examine the launch latency effects of different characteristics.
  if (test_kernel_launch) {
    testNullKernel(tiny_block_size, large_block_size / tiny_block_size, gpu, &timer);
    testNullKernel(small_block_size, large_block_size / small_block_size, gpu, &timer);
    testNullKernel(medium_block_size, large_block_size / medium_block_size, gpu, &timer);
    testNullKernel(large_block_size, 1, gpu, &timer);
    testArgumentLoadedKernel(tiny_block_size, large_block_size / tiny_block_size, gpu, &timer);
    testArgumentLoadedKernel(small_block_size, large_block_size / small_block_size, gpu, &timer);
    testArgumentLoadedKernel(medium_block_size, large_block_size / medium_block_size, gpu, &timer);
    testArgumentLoadedKernel(large_block_size, 1, gpu, &timer);
    testGmemKernel(tiny_block_size, large_block_size / tiny_block_size, gpu, &timer);
    testGmemKernel(small_block_size, large_block_size / small_block_size, gpu, &timer);
    testGmemKernel(medium_block_size, large_block_size / medium_block_size, gpu, &timer);
    testGmemKernel(large_block_size, 1, gpu, &timer);
  }

  // Run different classes of molecules.  This will stress-test the code as well as provide
  // performance curves with different sizes of molecules.
  runBatch("Dipeptides", gpu, oe, &timer, iter);
  runBatch("Tripeptides", gpu, oe, &timer, iter);
  runBatch("Tetrapeptides", gpu, oe, &timer, iter);
  if (extra_top.size() > 0) {
    runAdditionalTest(extra_top, extra_crd, gpu, &timer, iter, replicas);
  }
  
  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
}
