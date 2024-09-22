#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Accelerator/core_kernel_manager.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/Math/math_enumerators.h"
#include "../../src/MolecularMechanics/mm_controls.h"
#include "../../src/Namelists/command_line_parser.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/cellgrid.h"
#include "../../src/Potential/energy_enumerators.h"
#ifdef STORMM_USE_HPC
#  include "../../src/Potential/hpc_pme_potential.h"
#endif
#include "../../src/Potential/local_exclusionmask.h"
#include "../../src/Potential/pme_potential.h"
#include "../../src/Potential/pme_util.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/tile_manager.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/hpc_phasespace_synthesis.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinate_copy.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/unit_test_enumerators.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::namelist;
using namespace stormm::parse;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;

//-------------------------------------------------------------------------------------------------
// Determine the number and sizes of atom batches needed to span a given total number of atoms.
//
// Arguments:
//   total_atoms:  The total number of atoms to cover
//   max_batch:    The maximum batch size to use in covering all atoms
//   min_batch:    The minimum batch size to use in covering all atoms
//-------------------------------------------------------------------------------------------------
std::vector<int> batchList(const int total_atoms, const int max_batch, const int min_batch) {
  int ncovered = 0;
  int bsize = max_batch;
  int bmask = min_batch - 1;
  std::vector<int> result;
  const int twice_min_batch = min_batch * 2;
  while (ncovered < total_atoms) {
    if (total_atoms - ncovered >= bsize || bsize == min_batch ||
        (bsize == twice_min_batch && total_atoms - ncovered > (bsize >> 1))) {
      result.push_back(bsize);
      ncovered += bsize;
    }
    else if (bsize >= min_batch) {
      bsize >>= 1;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Analyze the theoretical occupancy of the tower-plate kernel.
//
// Arguments:
//   cg:       The cell grid to analyze
//   sys_idx:  Index of the system of interest within the cell grid
//   cutoff:   The particle-particle cutoff determining valid interactions
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void reportTowerPlateOccupancy(const CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> &cg, const int sys_idx,
                               const double cutoff) {

  // Construct the tower-plate arrangement
  const std::vector<int> tw_rel_a = {  0,  0,  0,  0,  0 };
  const std::vector<int> pl_rel_a = { -2, -1,  0,  1,  2, -2, -1,  0,  1,  2, -2, -1 };
  const std::vector<int> tw_rel_b = {  0,  0,  0,  0,  0 };
  const std::vector<int> pl_rel_b = { -2, -2, -2, -2, -2, -1, -1, -1, -1, -1,  0,  0 };
  const std::vector<int> tw_rel_c = { -2, -1,  0,  1,  2 };
  const std::vector<int> pl_rel_c = {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 };
  const CellGridReader<Tcoord, Tacc, Tcalc, Tcoord4> cgr = cg.data();
  const ullint sys_dims = cgr.system_cell_grids[sys_idx];
  const int na_cell = ((sys_dims >> 28) & 0xfff);
  const int nb_cell = ((sys_dims >> 40) & 0xfff);
  const int nc_cell = ((sys_dims >> 52) & 0xfff);
  const int cell_start = (sys_dims & 0xfffffff);
  const int nabc_cells = na_cell * nb_cell * nc_cell;
  std::vector<int> tower_cells(5), plate_cells(12), tower_offsets(5), plate_offsets(12);
  std::vector<int> tower_prefix(6), plate_prefix(13);
  double tp_total_work = 0.0;
  double tp_ideal_work = 0.0;
  double tp_successful = 0.0;
  const int xfrm_offset = sys_idx * roundUp(9, warp_size_int);
  std::vector<double> cg_invu(9);
  for (int i = 0; i < 9; i++) {
    cg_invu[i] = cgr.system_cell_invu[xfrm_offset + i];
  }
  const double cut_sq = cutoff * cutoff;
  std::vector<int> populations(256, 0);
  for (int i = 0; i < na_cell; i++) {
    for (int j = 0; j < nb_cell; j++) {
      for (int k = 0; k < nc_cell; k++) {

        // Determine the cell identities
        for (int m = 0; m < 5; m++) {
          int c_a = i + tw_rel_a[m];
          int c_b = j + tw_rel_b[m];
          int c_c = k + tw_rel_c[m];
          c_a += ((c_a < 0) - (c_a >= na_cell)) * na_cell;
          c_b += ((c_b < 0) - (c_b >= nb_cell)) * nb_cell;
          c_c += ((c_c < 0) - (c_c >= nc_cell)) * nc_cell;
          tower_cells[m] = (((c_c * nb_cell) + c_b) * na_cell) + c_a;
          const uint2 cldat = cgr.cell_limits[tower_cells[m]];
          tower_offsets[m] = cldat.x;
          tower_prefix[m] = (cldat.y >> 16);
        }
        tower_prefix[5] = 0;
        prefixSumInPlace<int>(&tower_prefix, PrefixSumType::EXCLUSIVE);
        populations[tower_prefix[3] - tower_prefix[2]] += 1;
        for (int m = 0; m < 12; m++) {
          int c_a = i + pl_rel_a[m];
          int c_b = j + pl_rel_b[m];
          int c_c = k + pl_rel_c[m];
          c_a += ((c_a < 0) - (c_a >= na_cell)) * na_cell;
          c_b += ((c_b < 0) - (c_b >= nb_cell)) * nb_cell;
          c_c += ((c_c < 0) - (c_c >= nc_cell)) * nc_cell;
          plate_cells[m] = (((c_c * nb_cell) + c_b) * na_cell) + c_a;
          const uint2 cldat = cgr.cell_limits[plate_cells[m]];
          plate_offsets[m] = cldat.x;
          plate_prefix[m] = (cldat.y >> 16);
        }
        plate_prefix[12] = 0;
        prefixSumInPlace<int>(&plate_prefix, PrefixSumType::EXCLUSIVE);
        tp_ideal_work += static_cast<double>(tower_prefix[5] * plate_prefix[12]);

        // Determine the number and sizes of tower sets for tower-plate interactions
        const std::vector<int> tower_sets = batchList(tower_prefix[5], 32, 8);
        const std::vector<int> plate_sets = batchList(plate_prefix[12], 32, 8);
        const std::vector<int> centr_sets = batchList(tower_prefix[3] - tower_prefix[2], 32, 8);
        const std::vector<int> lower_sets = batchList(tower_prefix[2], 32, 8);
        const double tower_atoms_tiled = sum<int>(tower_sets);
        const double plate_atoms_tiled = sum<int>(plate_sets);
        const double centr_atoms_tiled = sum<int>(centr_sets);
        const double lower_atoms_tiled = sum<int>(lower_sets);
        tp_total_work += (tower_atoms_tiled * plate_atoms_tiled) +
                          (centr_atoms_tiled * ((centr_atoms_tiled / 2) + lower_atoms_tiled));

        // Stage the tower and plate atoms
        std::vector<double> tower_x(tower_prefix[5]);
        std::vector<double> tower_y(tower_prefix[5]);
        std::vector<double> tower_z(tower_prefix[5]);
        std::vector<double> plate_x(plate_prefix[12]);
        std::vector<double> plate_y(plate_prefix[12]);
        std::vector<double> plate_z(plate_prefix[12]);
        size_t npt = 0;
        for (int m = 0; m < 5; m++) {
          const uint plim = tower_offsets[m] + tower_prefix[m + 1] - tower_prefix[m];
          const double stack_mult = m - 2;
          const double x_del = stack_mult * cg_invu[6];
          const double y_del = stack_mult * cg_invu[7];
          const double z_del = stack_mult * cg_invu[8];
          for (uint pos = tower_offsets[m]; pos < plim; pos++) {
            const Tcoord4 atom_img = cgr.image[pos];
            tower_x[npt] = atom_img.x + x_del;
            tower_y[npt] = atom_img.y + y_del;
            tower_z[npt] = atom_img.z + z_del;
            npt++;
          }
        }
        npt = 0;
        for (int m = 0; m < 12; m++) {
          const uint plim = plate_offsets[m] + plate_prefix[m + 1] - plate_prefix[m];
          const double stack_mult_x = static_cast<double>(m - ((m / 5) * 5) - 2);
          const double stack_mult_y = static_cast<double>((m / 5) - 2);
          const double x_del = (stack_mult_x * cg_invu[0]) + (stack_mult_y * cg_invu[3]);
          const double y_del =                               (stack_mult_y * cg_invu[4]);
          for (uint pos = plate_offsets[m]; pos < plim; pos++) {
            const Tcoord4 atom_img = cgr.image[pos];
            plate_x[npt] = atom_img.x + x_del;
            plate_y[npt] = atom_img.y + y_del;
            plate_z[npt] = atom_img.z;
            npt++;
          }
        }

        // Calculate the total number of interactions that are in range
        int n_success = 0;
        for (int m = 0; m < tower_prefix[5]; m++) {
          for (int  n = 0; n < plate_prefix[12]; n++) {
            const double dx = plate_x[n] - tower_x[m];
            const double dy = plate_y[n] - tower_y[m];
            const double dz = plate_z[n] - tower_z[m];
            const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
            if (r2 < cut_sq) {
              n_success++;
            }
          }
        }
        for (int m = tower_prefix[2]; m < tower_prefix[3]; m++) {
          for (int n = tower_prefix[2]; n < m; n++) {
            const double dx = tower_x[n] - tower_x[m];
            const double dy = tower_y[n] - tower_y[m];
            const double dz = tower_z[n] - tower_z[m];
            const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
            if (r2 < cut_sq) {
              n_success++;
            }
          }
          for (int n = tower_prefix[0]; n < tower_prefix[2]; n++) {
            const double dx = tower_x[n] - tower_x[m];
            const double dy = tower_y[n] - tower_y[m];
            const double dz = tower_z[n] - tower_z[m];
            const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
            if (r2 < cut_sq) {
              n_success++;
            }
          }
        }
        tp_successful += static_cast<double>(n_success);

        // Calculate the number of trivial exclusions based on distance to line and distance to
        // plane.  The first plane of interest is defined by the unit cell A and C vectors, the
        // second by the unit cell B and C vectors.  The AB plane (the Cartesian XY plane) is also
        // of interest and most simple to compute in the general case (just consider the particle's
        // Cartesian Z coordinate), but as a culling criterion for the tower itself it could be
        // useful.  The tower atoms would need to be re-ordered so that the small percentage of
        // tower atoms that fall into this category of things that could never interact with
        // anything in the plate are more likely to be concentrated into a single batch.
        // Otherwise, the distances to six critical lines will define whether a particle is within
        // range of anything in the central column.
        
        
      }
    }
  }
  printf("  Cutoff %9.4lf : Real Atom Content     %9.4lf\n"
         "                     Successful Real Pairs %9.4lf\n"
         "                     Success in All Pairs  %9.4lf\n", cutoff,
         (tp_ideal_work / tp_total_work) * 100.0, (tp_successful / tp_ideal_work) * 100.0,
         (tp_successful / tp_total_work) * 100.0);
  int pop_max = 255;
  while (pop_max >= 0 && populations[pop_max] == 0) {
    pop_max--;
  }
  printf("  Cell Population  Count\n  ---------------  -----\n");
  for (int i = 0; i < pop_max; i++) {
    printf("        %3d      %5d\n", i, populations[i]);
  }
  printf("\n");
}

//-------------------------------------------------------------------------------------------------
// Analyze the theoretical occupancy of a kernel based on the half-shell method with a trivial
// rejection criterion.
//
// Arguments:
//   
//-------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------
// Check the forces computed by the GPU kernel against CPU results.
//
// Arguments:
//   poly_ps:    The synthesis of coordinates, including forces on all particles
//   eval_frc:   Indicate whether forces were evaluated as part of the run
//   timer:      Time tracking object, to absorb the time running the prior CPU calculation as well
//               as the test to compare the results
//   time_test:  The integer code for the test time tracking
//-------------------------------------------------------------------------------------------------
void checkPairBenchmarkForces(const PhaseSpaceSynthesis &poly_ps, const EvaluateForce eval_frc,
                              StopWatch *timer, const int time_test) {
#ifdef STORMM_USE_HPC
  const UnitCellType system_uc = poly_ps.getUnitCellType();
  CoordinateFrame cpu_frc(poly_ps.getAtomCount(0), system_uc, HybridFormat::HOST_MOUNTED);
  CoordinateFrame gpu_frc(poly_ps.getAtomCount(0), system_uc, HybridFormat::HOST_MOUNTED);
  const CoordinateFrameReader cpu_cfr = cpu_frc.data();
  const CoordinateFrameReader gpu_cfr = gpu_frc.data();
  std::vector<double> cpu_chk(cpu_cfr.natom, 0.0);
  std::vector<double> gpu_chk(cpu_cfr.natom, 0.0);
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    switch (eval_frc) {
    case EvaluateForce::YES:
      coordCopy(&cpu_frc, poly_ps, i, TrajectoryKind::FORCES, HybridTargetLevel::HOST,
                HybridTargetLevel::HOST);
      coordCopy(&gpu_frc, poly_ps, i, TrajectoryKind::FORCES, HybridTargetLevel::HOST,
                HybridTargetLevel::DEVICE);
      for (int j = 0; j < cpu_cfr.natom; j++) {
        cpu_chk[j] = cpu_cfr.xcrd[j];
        gpu_chk[j] = gpu_cfr.xcrd[j];
      }
      check(gpu_chk, RelationalOperator::EQUAL, Approx(cpu_chk).margin(1.0e-2), "Forces on "
            "particles in system " + std::to_string(i) + " along the Cartesian X axis calculated "
            "using the GPU kernel running in " + getEnumerationName(PrecisionModel::SINGLE) +
            "-precision do not agree with CPU forces calculated in " +
            getEnumerationName(PrecisionModel::DOUBLE) + "-precision.");
      for (int j = 0; j < cpu_cfr.natom; j++) {
        cpu_chk[j] = cpu_cfr.ycrd[j];
        gpu_chk[j] = gpu_cfr.ycrd[j];
      }
      check(gpu_chk, RelationalOperator::EQUAL, Approx(cpu_chk).margin(1.0e-2), "Forces on "
            "particles in system " + std::to_string(i) + " along the Cartesian Y axis calculated "
            "using the GPU kernel running in " + getEnumerationName(PrecisionModel::SINGLE) +
            "-precision do not agree with CPU forces calculated in " +
            getEnumerationName(PrecisionModel::DOUBLE) + "-precision.");
      for (int j = 0; j < cpu_cfr.natom; j++) {
        cpu_chk[j] = cpu_cfr.zcrd[j];
        gpu_chk[j] = gpu_cfr.zcrd[j];
      }
      check(gpu_chk, RelationalOperator::EQUAL, Approx(cpu_chk).margin(1.0e-2), "Forces on "
            "particles in system " + std::to_string(i) + " along the Cartesian Z axis calculated "
            "using the GPU kernel running in " + getEnumerationName(PrecisionModel::SINGLE) +
            "-precision do not agree with CPU forces calculated in " +
            getEnumerationName(PrecisionModel::DOUBLE) + "-precision.");
      break;
    case EvaluateForce::NO:
      break;
    }
  }
#endif
  timer->assignTime(time_test);
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main (const int argc, const char* argv[]) {

  // Baseline variables
  StopWatch timer;
#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> force_gpu_to_engage(1);
#endif
  
  // Lay out time categories for program profiling
  const int time_load  = timer.addCategory("Load input files");
  const int time_build = timer.addCategory("Class object setup");
  const int time_test  = timer.addCategory("CPU testing");
  const int time_init  = timer.addCategory("Initialization");
  const int time_pairs = timer.addCategory("Pairs kernel with initialization");
  
  // Take in command line inputs
  CommandLineParser clip("pair_interactions", "A benchmarking program for measuring the rate at "
                         "which various GPU kernels can compute all pairwise particle-particle "
                         "interactions within one or more systems.", { "-timings" });
  clip.addStandardAmberInputs("-c", "-p", "-ig_seed");
  clip.addStandardBenchmarkingInputs({ "-iter", "-trials", "-replicas", "-cutoff", "-elec_cutoff",
                                       "-vdw_cutoff", "-pad" });

  // Custom inputs for this benchmarking program
  NamelistEmulator* t_nml = clip.getNamelistPointer();
  t_nml->addKeyword("-dual_cg", NamelistType::BOOLEAN);
  t_nml->addHelp("-dual_cg", "Force the use of dual neighbor lists, if the electrostatic and "
                 "van-der Waals cutoffs are not already distinct.");
  t_nml->addKeyword("-skip_cpu_check", NamelistType::BOOLEAN);
  t_nml->addHelp("-skip_cpu_check", "Skip a CPU_based check on the forces computed by the GPU.");
  t_nml->addKeyword("-occupancy", NamelistType::BOOLEAN);
  t_nml->addHelp("-occupancy", "Compute the occupancy of GPU warps during the pairs calculation.");
  t_nml->addKeyword("-eval_nrg", NamelistType::BOOLEAN);
  t_nml->addHelp("-eval_nrg", "Request that the overall system non-bonded energy components be "
                 "evaluated, in addition to the forces on particles.");
  t_nml->addKeyword("-omit_frc", NamelistType::BOOLEAN);
  t_nml->addHelp("-omit_frc", "Request that the forces on particles be omitted from the GPU "
                 "calculation.  This will also omit the CPU check on forces, and compel an "
                 "evaluation of the energy so as to have at least one quantity for the GPU to "
                 "compute.");
  t_nml->addKeyword("-fp_bits", NamelistType::INTEGER, std::to_string(24));
  t_nml->addHelp("-fp_bits", "The number of fixed precision bits after the point (values in "
                 "kcal/mol-A) with which to accumulate forces on all particles.");
  t_nml->addKeyword("-warp_mult", NamelistType::INTEGER, std::to_string(1));
  t_nml->addHelp("-warp_mult", "The number of warps to devote to processing each neighbor list "
                 "cell's assigned pair interactions.");

  // Initialize the testing environment such that it cooperates with this program's own
  // CommandLineParser to read user input.
  TestEnvironment oe(argc, argv, &clip, TmpdirStatus::NOT_REQUIRED, ExceptionResponse::SILENT);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Read command line instructions.
  clip.parseUserInput(argc, argv);
  const std::string inpcrd_file = t_nml->getStringValue("-c");
  const std::string topology_file = t_nml->getStringValue("-p");
  bool use_dual_cg = t_nml->getBoolValue("-dual_cg");
  bool skip_cpu_check = t_nml->getBoolValue("-skip_cpu_check");
  bool est_occupancy = t_nml->getBoolValue("-occupancy");
  int n_trials = t_nml->getIntValue("-trials");
  int n_repeats = t_nml->getIntValue("-iter");
  int n_replicas = t_nml->getIntValue("-replicas");
  int fp_bits = t_nml->getIntValue("-fp_bits");
  int warp_mult = t_nml->getIntValue("-warp_mult");
  int ig_seed = t_nml->getIntValue("-ig_seed");
  double elec_cutoff = t_nml->getRealValue("-elec_cutoff");
  double vdw_cutoff = t_nml->getRealValue("-vdw_cutoff");
  double cutoff_pad = t_nml->getRealValue("-pad");
  EvaluateEnergy eval_nrg = (t_nml->getBoolValue("-eval_nrg")) ? EvaluateEnergy::YES :
                                                                 EvaluateEnergy::NO;
  EvaluateForce eval_frc;
  if (t_nml->getBoolValue("-omit_frc")) {

    // Force the evaluation of energy if forces are to be omitted.
    eval_nrg = EvaluateEnergy::YES;
    eval_frc = EvaluateForce::NO;
  }
  else {
    eval_frc = EvaluateForce::YES;
  }

  // Input checks
  if (n_replicas <= 0) {
    rtErr("A replica count of " + std::to_string(n_replicas) + " is invalid.\n",
          "pair_interactions");
  }
  
  // A Hybrid object was created to engage the GPU.  Absorb any bootup time into "miscellaneous."
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
  }
  
  // Check that a system has been provided.  Load the basic topology and input coordinates.
#ifdef STORMM_USE_HPC
  TestSystemManager tsm(std::vector<std::string>(1, topology_file),
                        std::vector<std::string>(1, inpcrd_file), ExceptionResponse::DIE);
  timer.assignTime(time_load);
  
  // Stage the input parameters: cutoffs for each interaction.
  MolecularMechanicsControls mmctrl(0.01, 1, 1, warp_mult, elec_cutoff, vdw_cutoff);

  // Create the basic objects
  PhaseSpaceSynthesis poly_ps = tsm.exportPhaseSpaceSynthesis(std::vector<int>(n_replicas, 0));
  PsSynthesisWriter host_psw = poly_ps.data();
  Xoshiro256ppGenerator xrs(ig_seed);
  for (int i = 1; i < host_psw.system_count; i++) {
    const size_t aoff = host_psw.atom_starts[i];
    addRandomNoise(&xrs, &host_psw.xcrd[aoff], &host_psw.xcrd_ovrf[aoff], &host_psw.ycrd[aoff],
                   &host_psw.ycrd_ovrf[aoff], &host_psw.zcrd[aoff], &host_psw.zcrd_ovrf[aoff],
                   host_psw.atom_counts[i], 0.01, host_psw.gpos_scale_f);
  }
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(std::vector<int>(n_replicas, 0));
  ScoreCard sc(poly_ps.getSystemCount(), 1, 32);
  use_dual_cg = (use_dual_cg || (fabs(elec_cutoff - vdw_cutoff) > 1.0e-6));
  if (use_dual_cg == false) {
    elec_cutoff = vdw_cutoff;
  }
  const NeighborListKind layout = (use_dual_cg) ? NeighborListKind::DUAL : NeighborListKind::MONO;
  CoreKlManager launcher(gpu, poly_ag);
  PPITable direct_space_table(NonbondedTheme::ELECTROSTATIC, BasisFunctions::MIXED_FRACTIONS,
                              TableIndexing::SQUARED_ARG, elec_cutoff);
  const double ew_coeff = direct_space_table.getEwaldCoefficient();
  timer.assignTime(time_build);
  
  // Upload the basic objects
  poly_ps.upload();
  poly_ag.upload();
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps.data(devc);
  SyNonbondedKit<float, float2> poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(devc);
  MMControlKit<float> ctrl = mmctrl.spData(devc);
  PPIKit<float, float4> nrg_tab = direct_space_table.spData(devc);
  ScoreCardWriter scw = sc.data(devc);
  timer.assignTime(time_build);

  // Create the neighbor list cell grids and run calculations
  if (use_dual_cg) {
    CellGrid<float, int, float, float4> cg_qq(&poly_ps, &poly_ag, 0.5 * elec_cutoff, cutoff_pad, 4,
                                              NonbondedTheme::ELECTROSTATIC);
    CellGrid<float, int, float, float4> cg_lj(&poly_ps, &poly_ag, 0.5 * vdw_cutoff, cutoff_pad, 4,
                                              NonbondedTheme::VAN_DER_WAALS);
    const TinyBoxPresence has_tiny_box = (cg_qq.getTinyBoxPresence() == TinyBoxPresence::YES ||
                                          cg_lj.getTinyBoxPresence() == TinyBoxPresence::YES) ?
                                         TinyBoxPresence::YES : TinyBoxPresence::NO;
    const int2 tp_bt = launcher.getPMEPairsKernelDims(PrecisionModel::SINGLE,
                                                      PrecisionModel::SINGLE,
                                                      NeighborListKind::DUAL, has_tiny_box,
                                                      eval_frc, eval_nrg, ClashResponse::NONE);
    TileManager tlmn(tp_bt);
    TilePlan tlpn = tlmn.data(devc);
    mmctrl.primeWorkUnitCounters(launcher, eval_frc, eval_nrg, ClashResponse::NONE,
                                 VwuGoal::MOVE_PARTICLES, PrecisionModel::SINGLE,
                                 PrecisionModel::SINGLE, QMapMethod::ACC_SHARED,
                                 PrecisionModel::SINGLE, float_type_index, 5, layout, has_tiny_box,
                                 poly_ag);
    mmctrl.upload();
    cg_qq.upload();
    cg_lj.upload();

    // Make the appropriate exclusion masks
    LocalExclusionMask lem(poly_ag, NonbondedTheme::ALL);
    lem.upload();
    const LocalExclusionMaskReader lemr = lem.data(devc);
    
    // Take abstracts and run the kernel repeatedly.
    CellGridWriter<float, int, float, float4> cg_qqw = cg_qq.data(devc);
    CellGridWriter<float, int, float, float4> cg_ljw = cg_lj.data(devc);
    const PsSynthesisBorders sysbrd = cg_qq.getUnitCellTransforms(devc);
    timer.assignTime(time_build);
    
    // Run the initialization of forces multiple times to gain a sense of the overhead paid in the
    // test.
    for (int i = 0; i < n_trials; i++) {
      for (int j = 0; j < n_repeats; j++) {
        cg_qq.initializeForces(devc, gpu);
        cg_lj.initializeForces(devc, gpu);
      }
      cudaDeviceSynchronize();
      timer.assignTime(time_init);
    }

    // Run the pair interactions evaluation and continue clearing the buffers
    for (int i = 0; i < n_trials; i++) {
      for (int j = 0; j < n_repeats; j++) {
        cg_qq.initializeForces(devc, gpu);
        cg_lj.initializeForces(devc, gpu);
        switch (has_tiny_box) {
        case TinyBoxPresence::NO:
          launchPMEPairs(poly_nbk, lemr, nrg_tab, &cg_qqw, &cg_ljw, &tlpn, &scw, &ctrl, eval_frc,
                         eval_nrg, tp_bt, 0.0, 0.0);
          break;
        case TinyBoxPresence::YES:
          launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cg_qqw, &cg_ljw, &tlpn, &scw, &ctrl,
                         eval_frc, eval_nrg, tp_bt, 0.0, 0.0);
          break;
        }          
        ctrl.step +=1;
      }
      cudaDeviceSynchronize();
      timer.assignTime(time_pairs);
    }
    poly_ps.initializeForces(gpu, devc);
    cg_qq.contributeForces(devc, gpu);
    cg_lj.contributeForces(devc, gpu);

    // Estimate the thread occupancy of the GPU during the calculation
    if (est_occupancy) {
      for (int i = 0; i < poly_ps.getSystemCount(); i++) {
        reportTowerPlateOccupancy(cg_qq, i, elec_cutoff);
        reportTowerPlateOccupancy(cg_lj, i, vdw_cutoff);
      }
    }
  }
  else {
    CellGrid<float, int, float, float4> cg(&poly_ps, &poly_ag, 0.5 * vdw_cutoff, cutoff_pad, 4,
                                           NonbondedTheme::ALL);
    const TinyBoxPresence has_tiny_box = cg.getTinyBoxPresence();
    const int2 tp_bt = launcher.getPMEPairsKernelDims(PrecisionModel::SINGLE,
                                                      PrecisionModel::SINGLE,
                                                      NeighborListKind::DUAL,
                                                      cg.getTinyBoxPresence(), eval_frc, eval_nrg,
                                                      ClashResponse::NONE);
    TileManager tlmn(tp_bt);
    TilePlan tlpn = tlmn.data(devc);
    mmctrl.primeWorkUnitCounters(launcher, eval_frc, eval_nrg, ClashResponse::NONE,
                                 VwuGoal::MOVE_PARTICLES, PrecisionModel::SINGLE,
                                 PrecisionModel::SINGLE, QMapMethod::ACC_SHARED,
                                 PrecisionModel::SINGLE, float_type_index, 5, layout,
                                 cg.getTinyBoxPresence(), poly_ag);
    mmctrl.upload();
    cg.upload();

    // Make the appropriate exclusion masks
    LocalExclusionMask lem(poly_ag, NonbondedTheme::ALL);
    lem.upload();
    const LocalExclusionMaskReader lemr = lem.data(devc);
    
    // Take abstracts and run the kernel repeatedly.
    CellGridWriter<float, int, float, float4> cgw = cg.data(devc);
    const PsSynthesisBorders sysbrd = cg.getUnitCellTransforms(devc);
    timer.assignTime(time_build);

    // Test the basic force initialization
    for (int i = 0; i < n_trials; i++) {
      for (int j = 0; j < n_repeats; j++) {
        cg.initializeForces(devc, gpu);
      }
      cudaDeviceSynchronize();
      timer.assignTime(time_init);
    }

    // Run the pair interactions evaluation and continue clearing the buffers
    for (int i = 0; i < n_trials; i++) {
      for (int j = 0; j < n_repeats; j++) {
        cg.initializeForces(devc, gpu);
        switch (cg.getTinyBoxPresence()) {
        case TinyBoxPresence::NO:
          launchPMEPairs(poly_nbk, lemr, nrg_tab, &cgw, &tlpn, &scw, &ctrl, eval_frc, eval_nrg,
                         tp_bt, 0.0, 0.0);
          break;
        case TinyBoxPresence::YES:
          launchPMEPairs(poly_nbk, lemr, nrg_tab, sysbrd, &cgw, &tlpn, &scw, &ctrl, eval_frc,
                         eval_nrg, tp_bt, 0.0, 0.0);
          break;
        }
        ctrl.step +=1;
      }
      cudaDeviceSynchronize();
      timer.assignTime(time_pairs);
    }
    poly_ps.initializeForces(gpu, devc);
    cg.contributeForces(devc, gpu);

    // Estimate the thread occupancy of the GPU during the calculation
    if (est_occupancy) {
      for (int i = 0; i < poly_ps.getSystemCount(); i++) {
        reportTowerPlateOccupancy(cg, i, vdw_cutoff);
      }
    }
  }

  // Run a separate check on the validity of the forces
  if (skip_cpu_check == false) {
    for (int i = 0; i < poly_ps.getSystemCount(); i++) {
      PhaseSpace ps_i(poly_ps.getAtomCount(i), poly_ps.getUnitCellType(),
                      HybridFormat::HOST_ONLY);
      coordCopy(&ps_i, poly_ps, i);
      ps_i.initializeForces();
      const AtomGraph *ag_i = poly_ps.getSystemTopologyPointer(i);
      const LocalExclusionMask lem_i(ag_i);
      evaluateParticleParticleEnergy(&ps_i, ag_i, lem_i, PrecisionModel::SINGLE, elec_cutoff,
                                     vdw_cutoff, ew_coeff);
      coordCopy(&poly_ps, i, ps_i);
    }
    checkPairBenchmarkForces(poly_ps, eval_frc, &timer, time_test);
  }
#else // STORMM_USE_HPC
  rtWarn("This benchmarking program requires GPU support to run.", "pair_interactions");
#endif // STORMM_USE_HPC

  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return 0;
}
