// -*-c++-*-
#include "../../src/Accelerator/core_kernel_manager.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/DataTypes/mixed_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/hpc_reduction.h"
#include "../../src/Math/reduction_abstracts.h"
#include "../../src/Math/reduction_bridge.h"
#include "../../src/Math/reduction_enumerators.h"
#include "../../src/Math/rounding.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/MolecularMechanics/hpc_minimization.h"
#include "../../src/MolecularMechanics/line_minimization.h"
#include "../../src/MolecularMechanics/mm_controls.h"
#include "../../src/MolecularMechanics/mm_evaluation.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/cacheresource.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/hpc_nonbonded_potential.h"
#include "../../src/Potential/hpc_valence_potential.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Structure/virtual_site_handling.h"
#include "../../src/Structure/hpc_virtual_site_handling.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/static_mask_synthesis.h"
#include "../../src/Synthesis/synthesis_abstracts.h"
#include "../../src/Synthesis/nonbonded_workunit.h"
#include "../../src/Synthesis/valence_workunit.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinate_copy.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/thermostat.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::stmath;
using namespace stormm::mm;
using namespace stormm::parse;
using namespace stormm::random;
using namespace stormm::structure;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Check that all forces on equivalent systems are equal, to within a reasonable tolerance.  It is
// assumed that systems described by the same topology are identical--no restraints, or at least
// no unique restraints, are present.
//
// Arguments:
//   poly_ps:  Compilation of coordinates and forces
//   poly_ag:  Compilation of topologies for all systems
//   err_msg:  Error to display if forces between equivalent systems do not match up
//   c_limit:  Tolerance for coordinates to disagree before being deemed inconsistent
//   f_limit:  Tolerance for forces to disagree before being deemed inconsistent
//   step_no:  Number of the minimization step
//-------------------------------------------------------------------------------------------------
bool checkConsistency(const PhaseSpaceSynthesis &poly_ps, const AtomGraphSynthesis &poly_ag,
                      const std::string &err_msg, const double c_limit, const double f_limit,
                      const int step_no) {
  int n_coord_mismatch = 0;
  int n_force_mismatch = 0;
  int n_syscrd_mismatch = 0;
  int n_sysfrc_mismatch = 0;
  PsSynthesisReader poly_psr = poly_ps.data();
  const std::vector<int> ag_indices = poly_ag.getTopologyIndices();
  std::vector<bool> covered(poly_psr.system_count, false);
  std::vector<bool> syscrd_ok(poly_psr.system_count, true);
  std::vector<bool> sysfrc_ok(poly_psr.system_count, true);
  for (int i = 0; i < poly_psr.system_count; i++) {
    if (covered[i]) {
      continue;
    }
    covered[i] = true;
    const int iag_no = ag_indices[i];
    for (int j = i + 1; j < poly_psr.system_count; j++) {
      if (covered[j] == false && ag_indices[j] == iag_no) {
        covered[j] = true;
        int itrack = poly_psr.atom_starts[i];
        const int jtrack_lim = poly_psr.atom_starts[j] + poly_psr.atom_counts[j];
        for (int jtrack = poly_psr.atom_starts[j]; jtrack < jtrack_lim; jtrack++) {
          if (fabs(poly_psr.xcrd[itrack] - poly_psr.xcrd[jtrack]) > f_limit ||
              fabs(poly_psr.ycrd[itrack] - poly_psr.ycrd[jtrack]) > f_limit ||
              fabs(poly_psr.zcrd[itrack] - poly_psr.zcrd[jtrack]) > f_limit) {
            n_coord_mismatch++;
            if (syscrd_ok[j]) {
              n_syscrd_mismatch++;
              syscrd_ok[j] = false;
            }
          }
          if (fabs(poly_psr.xfrc[itrack] - poly_psr.xfrc[jtrack]) > f_limit ||
              fabs(poly_psr.yfrc[itrack] - poly_psr.yfrc[jtrack]) > f_limit ||
              fabs(poly_psr.zfrc[itrack] - poly_psr.zfrc[jtrack]) > f_limit) {
            n_force_mismatch++;
            if (sysfrc_ok[j]) {
              n_sysfrc_mismatch++;
              sysfrc_ok[j] = false;
            }
          }
          itrack++;
        }
      }
    }
  }
  if (n_coord_mismatch > 0) {
    rtWarn("A total of " + std::to_string(n_coord_mismatch) + " atoms' coordinates were "
           "inconsistent among the first instance of a system and subsequent replicas with the "
           "same topology.  In all, " + std::to_string(n_syscrd_mismatch) + " systems displayed "
           "errors.  Checked at: " + err_msg + ", step " + std::to_string(step_no) + ".");
  }
  if (n_force_mismatch > 0) {
    rtWarn("A total of " + std::to_string(n_force_mismatch) + " atoms' forces were inconsistent "
           "among the first instance of a system and subsequent replicas with the same topology.  "
           "In all, " + std::to_string(n_sysfrc_mismatch) + " systems displayed errors.  Checked "
           "at: " + err_msg + ", step " + std::to_string(step_no) + ".");
  }
  return (n_coord_mismatch == 0 && n_force_mismatch == 0);
}

//-------------------------------------------------------------------------------------------------
// Set forces and coordinates for all systems to be exactly the same.
//
// Arguments:
//   poly_ps:  Compilation of coordinates and forces
//   poly_ag:  Compilation of topologies for all systems
//-------------------------------------------------------------------------------------------------
void mandateEquality(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag) {
  PsSynthesisWriter poly_psr = poly_ps->data();
  const std::vector<int> ag_indices = poly_ag.getTopologyIndices();
  std::vector<bool> covered(poly_psr.system_count, false);
  for (int i = 0; i < poly_psr.system_count; i++) {
    if (covered[i]) {
      continue;
    }
    covered[i] = true;
    const int iag_no = ag_indices[i];
    for (int j = i + 1; j < poly_psr.system_count; j++) {
      if (ag_indices[j] == iag_no) {
        covered[j] = true;
        int itrack = poly_psr.atom_starts[i];
        const int jtrack_lim = poly_psr.atom_starts[j] + poly_psr.atom_counts[j];
        for (int jtrack = poly_psr.atom_starts[j]; jtrack < jtrack_lim; jtrack++) {
          poly_psr.xfrc[jtrack] = poly_psr.xfrc[itrack];
          poly_psr.yfrc[jtrack] = poly_psr.yfrc[itrack];
          poly_psr.zfrc[jtrack] = poly_psr.zfrc[itrack];
          itrack++;
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Check that the components of the PhaseSpaceSynthesis relevant to conjugate gradient line
// minimizations are consistent with respect to identical systems.
//
// Arguments:
//   poly_ps:      Compilation of coordinates and forces
//   line_record:  Workspace for performing line minimizations
//   mol_id_vec:   List of molecule IDs found in the synthesis poly_ps
//   step_number:  Number of the minimization step
//   stage:        Stage of teh Conjugate Gradient minimization cycle
//-------------------------------------------------------------------------------------------------
void checkLineMinimizationComponents(PhaseSpaceSynthesis *poly_ps, LineMinimization *line_record,
                                     const std::vector<int> &mol_id_vec, const int step_number,
                                     const std::string &stage) {
  poly_ps->download();
  line_record->download();
  const PsSynthesisWriter poly_psw = poly_ps->data();
  const LinMinWriter lmw = line_record->data();

  // Check the components, system-by-system
  const int nsys = mol_id_vec.size();
  int n_unique_sys = maxValue(mol_id_vec) + 1;
  std::vector<int> sys_copies(n_unique_sys, 0);
  std::vector<std::vector<int>> sys_map(n_unique_sys);
  for (int sysid = 0; sysid < nsys; sysid++) {
    sys_copies[mol_id_vec[sysid]] += 1;
  }
  for (int i = 0; i < n_unique_sys; i++) {
    sys_map[i].reserve(sys_copies[i]);
  }
  for (int sysid = 0; sysid < nsys; sysid++) {
    sys_map[mol_id_vec[sysid]].push_back(sysid);
  }
  for (int i = 0; i < n_unique_sys; i++) {

    // Check individual atoms
    const int natom = poly_psw.atom_counts[sys_map[i][0]];
    std::vector<double> x_crd(sys_copies[i], 0.0);
    std::vector<double> y_crd(sys_copies[i], 0.0);
    std::vector<double> z_crd(sys_copies[i], 0.0);
    for (int atomcon = 0; atomcon < natom; atomcon++) {
      for (int j = 0; j < sys_copies[i]; j++) {
        const int sysid = sys_map[i][j];
        const int atom_llim = poly_psw.atom_starts[sysid];
        const int readidx = atom_llim + atomcon;
        x_crd[j] = poly_psw.xcrd[readidx] + (max_llint_accumulation * poly_psw.xcrd_ovrf[readidx]);
        y_crd[j] = poly_psw.ycrd[readidx] + (max_llint_accumulation * poly_psw.ycrd_ovrf[readidx]);
        z_crd[j] = poly_psw.zcrd[readidx] + (max_llint_accumulation * poly_psw.zcrd_ovrf[readidx]);
        x_crd[j] *= poly_psw.inv_gpos_scale;
        y_crd[j] *= poly_psw.inv_gpos_scale;
        z_crd[j] *= poly_psw.inv_gpos_scale;
      }
      addScalarToVector(&x_crd, -(mean(x_crd)));
      addScalarToVector(&y_crd, -(mean(y_crd)));
      addScalarToVector(&z_crd, -(mean(z_crd)));
      const double dx_max = maxAbsValue(x_crd);
      const double dy_max = maxAbsValue(y_crd);
      const double dz_max = maxAbsValue(z_crd);
      if (dx_max > small) {
        printf("A maximum deviation of %12.8lf is observed in X coordinates of atom %4d, system "
               "%4zu.  Step: %4d  Stage: %s\n", dx_max, atomcon, locateValue(x_crd, dx_max),
               step_number, stage.c_str());
        exit(1);
      }
      if (dy_max > small) {
        printf("A maximum deviation of %12.8lf is observed in Y coordinates of atom %4d, system "
               "%4zu.  Step: %4d  Stage: %s\n", dy_max, atomcon, locateValue(y_crd, dy_max),
               step_number, stage.c_str());
        exit(1);
      }
      if (dz_max > small) {
        printf("A maximum deviation of %12.8lf is observed in Z coordinates of atom %4d, system "
               "%4zu.  Step: %4d  Stage: %s\n", dz_max, atomcon, locateValue(z_crd, dz_max),
               step_number, stage.c_str());
        exit(1);
      }
    }

    // Check the line minimization data
    std::vector<double> l_move(sys_copies[i], 0.0);
    std::vector<double> s_move(sys_copies[i], 0.0);
    std::vector<double> mfac_a(sys_copies[i], 0.0);
    std::vector<double> mfac_b(sys_copies[i], 0.0);
    std::vector<double> mfac_c(sys_copies[i], 0.0);
    std::vector<double> nrg_a(sys_copies[i], 0.0);
    std::vector<double> nrg_b(sys_copies[i], 0.0);
    std::vector<double> nrg_c(sys_copies[i], 0.0);
    std::vector<double> nrg_d(sys_copies[i], 0.0);
    for (int j = 0; j < sys_copies[i]; j++) {
      const int sysid = sys_map[i][j];
      l_move[j] = lmw.l_move[sysid];
      s_move[j] = lmw.s_move[sysid];
      mfac_a[j] = lmw.mfac_a[sysid];
      mfac_b[j] = lmw.mfac_b[sysid];
      mfac_c[j] = lmw.mfac_c[sysid];
      nrg_a[j] = lmw.nrg_a[sysid];
      nrg_b[j] = lmw.nrg_b[sysid];
      nrg_c[j] = lmw.nrg_c[sysid];
      nrg_d[j] = lmw.nrg_d[sysid];
    }
    addScalarToVector(&l_move, -mean(l_move));
    addScalarToVector(&s_move, -mean(s_move));
    addScalarToVector(&mfac_a, -mean(mfac_a));
    addScalarToVector(&mfac_b, -mean(mfac_b));
    addScalarToVector(&mfac_c, -mean(mfac_c));
    addScalarToVector(&nrg_a, -mean(nrg_a));
    addScalarToVector(&nrg_b, -mean(nrg_b));
    addScalarToVector(&nrg_c, -mean(nrg_c));
    addScalarToVector(&nrg_d, -mean(nrg_d));
    const double dlmove_max = maxAbsValue(l_move);
    const double dsmove_max = maxAbsValue(s_move);
    const double dmfaca_max = maxAbsValue(mfac_a);
    const double dmfacb_max = maxAbsValue(mfac_b);
    const double dmfacc_max = maxAbsValue(mfac_c);
    const double dnrga_max = maxAbsValue(nrg_a);
    const double dnrgb_max = maxAbsValue(nrg_b);
    const double dnrgc_max = maxAbsValue(nrg_c);
    const double dnrgd_max = maxAbsValue(nrg_d);
    if (dlmove_max > small) {
      printf("A maximum deviation of %12.8lf is observed in the L-move.  Step: %4d  Stage: %s\n",
             dlmove_max, step_number, stage.c_str());
      exit(1);
    }
    if (dsmove_max > small) {
      printf("A maximum deviation of %12.8lf is observed in the S-move.  Step: %4d  Stage: %s\n",
             dsmove_max, step_number, stage.c_str());
      exit(1);
    }
    if (dmfaca_max > small) {
      printf("A maximum deviation of %12.8lf is observed in M-fac (A).  Step: %4d  Stage: %s\n",
             dmfaca_max, step_number, stage.c_str());
      exit(1);
    }
    if (dmfacb_max > small) {
      printf("A maximum deviation of %12.8lf is observed in M-fac (B).  Step: %4d  Stage: %s\n",
             dmfacb_max, step_number, stage.c_str());
      exit(1);
    }
    if (dmfacc_max > small) {
      printf("A maximum deviation of %12.8lf is observed in M-fac (C).  Step: %4d  Stage: %s\n",
             dmfacc_max, step_number, stage.c_str());
      exit(1);
    }
    if (dnrga_max > small) {
      printf("A maximum deviation of %12.8lf is observed in energy (A).  Step: %4d  Stage: %s\n",
             dnrga_max, step_number, stage.c_str());
      printf("Nrg A = [\n");
      int k = 0;
      for (size_t j = 0; j < nrg_a.size(); j++) {
        printf("  %10.4lf", nrg_a[j]);
        k++;
        if (k == 16) {
          k = 0;
          printf("\n");
        }
      }
      if (k > 0) {
        printf("\n");
      }
      printf("];\n");
      exit(1);
    }
    if (dnrgb_max > small) {
      printf("A maximum deviation of %12.8lf is observed in energy (B).  Step: %4d  Stage: %s\n",
             dnrgb_max, step_number, stage.c_str());
      exit(1);
    }
    if (dnrgc_max > small) {
      printf("A maximum deviation of %12.8lf is observed in energy (C).  Step: %4d  Stage: %s\n",
             dnrgc_max, step_number, stage.c_str());
      exit(1);
    }
    if (dnrgd_max > small) {
      printf("A maximum deviation of %12.8lf is observed in energy (D).  Step: %4d  Stage: %s\n",
             dnrgd_max, step_number, stage.c_str());
      exit(1);
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Perform energy minimization with meticulous checks to ensure that the process is consistent and
// reproducible.
//
// Arguments:
//   ag_ptr_vec:          Vector of pointers to unique topologies for systems in ps_vec
//   ps_vec:              Coordinates for all systems, to be replicated in the resulting synthesis
//   mol_id_vec:          Indication of how to replicate various structures described in ps_vec
//   gpu:                 Details of the GPU available
//   do_tests:            Indicate whether tests are possible to run
//   oe:                  Contains the name of the STORMM source path from shell variables as well
//                        as information on whether to write snapshot files or do the comparisons
//   psnap:               Instructions as to whether to begin printing a new snapshot file or
//                        append to an existing one, if snapshots are to be written
//   snap_name:           Name of the snapshot file for final energies of systems
//   prec:                Precision model for arithmetic and fixed-precision representations
//   gpos_bits:           Fixed-precision bits after the decimal in the positional representation
//   frc_bits:            Fixed-precision bits after the decimal in force accumulation
//   maxcyc:              Maximum number of minimization cycles
//   enforce_same_track:  Flag to have the systems explicitly set to keep on the same track if
//                        they start to diverge very slightly
//   check_mm:            Flag to have molecular mechanics energies and forces of the final states
//                        checked
//   frc_tol:             Tolerance for force comparisons
//   nrg_tol:             Tolerance for energy comparisons
//   test_name:           Name given to the collection of tests 
//   timer:               Optional profiling tool, if of interest for a brief performance check
//-------------------------------------------------------------------------------------------------
void metaMinimization(const std::vector<AtomGraph*> &ag_ptr_vec,
                      const std::vector<PhaseSpace> &ps_vec, const std::vector<int> &mol_id_vec,
                      const GpuDetails &gpu, const TestPriority do_tests,
                      const TestEnvironment &oe,
                      const PrintSituation psnap = PrintSituation::OVERWRITE, 
                      const std::string &snap_name = std::string(""),
                      const std::string &var_name = std::string(""),
                      const PrecisionModel prec = PrecisionModel::DOUBLE, const int gpos_bits = 40,
                      const int frc_bits = 40, const int maxcyc = 500,
                      const bool enforce_same_track = true, const bool check_mm = true,
                      const double frc_tol = 1.0e-6, const double nrg_tol = 1.0e-6,
                      const std::string &test_name = std::string(""), StopWatch *timer = nullptr,
                      const ImplicitSolventModel ism_choice = ImplicitSolventModel::NONE) {
  AtomGraphSynthesis poly_ag(ag_ptr_vec, mol_id_vec, ExceptionResponse::SILENT, gpu, timer);
  const NeckGeneralizedBornTable ngb_tables;
  if (ism_choice != ImplicitSolventModel::NONE) {
    const AtomicRadiusSet radii = (ism_choice == ImplicitSolventModel::NECK_GB_II) ?
                                  AtomicRadiusSet::MBONDI3 : AtomicRadiusSet::BONDI;
    poly_ag.setImplicitSolventModel(ism_choice, ngb_tables, radii);
  }
  StaticExclusionMaskSynthesis poly_se(poly_ag.getUniqueTopologies(), mol_id_vec);
  poly_ag.loadNonbondedWorkUnits(poly_se);
  PhaseSpaceSynthesis poly_ps(ps_vec, ag_ptr_vec, mol_id_vec, gpos_bits, 24, 40, frc_bits);
    
  // Create the minimization instructions
  MinimizeControls mincon;
  mincon.setTotalCycles(maxcyc);
  mincon.setSteepestDescentCycles(maxcyc / 10);
  mincon.setClashDampingCycles(0);
    
  // Create separate molecular mechanics control objects based on the minimization operations
  // for each of the ways that the valence kernel gets subdivided. 
  MolecularMechanicsControls mmctrl_fe(mincon);
  MolecularMechanicsControls mmctrl_xe(mincon);
  
  // Track energies in the systems
  ScoreCard sc(mol_id_vec.size(), mincon.getTotalCycles(), 32);
  
  // Obtain kernel launch parameters for the workload
  CoreKlManager launcher(gpu, poly_ag);
      
  // Lay out GPU cache resources
  const int2 vale_fe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::YES,
                                                        EvaluateEnergy::YES,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE, ClashResponse::NONE);
  const int2 vale_xe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::NO,
                                                        EvaluateEnergy::YES,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE, ClashResponse::NONE);
  const int2 vste_mv_lp = launcher.getVirtualSiteKernelDims(prec, VirtualSiteActivity::PLACEMENT);
  const int2 vste_xm_lp = launcher.getVirtualSiteKernelDims(prec,
                                                            VirtualSiteActivity::TRANSMIT_FORCES);
  const NbwuKind nb_work_type = poly_ag.getNonbondedWorkType();
  const ImplicitSolventModel ism_type = poly_ag.getImplicitSolventModel();
  const int2 nonb_lp = launcher.getNonbondedKernelDims(prec, nb_work_type, EvaluateForce::YES,
                                                       EvaluateEnergy::YES,
                                                       AccumulationMethod::SPLIT, ism_type,
                                                       ClashResponse::NONE);
  const int2 gbr_lp = launcher.getBornRadiiKernelDims(prec, nb_work_type,
                                                      AccumulationMethod::SPLIT, ism_type);
  const int2 gbd_lp = launcher.getBornDerivativeKernelDims(prec, nb_work_type,
                                                           AccumulationMethod::SPLIT, ism_type);
  const int2 redu_lp = launcher.getReductionKernelDims(prec, ReductionGoal::CONJUGATE_GRADIENT,
                                                       ReductionStage::ALL_REDUCE);
  CacheResource valence_fe_tb_reserve(vale_fe_lp.x, maximum_valence_work_unit_atoms);
  CacheResource valence_xe_tb_reserve(vale_xe_lp.x, maximum_valence_work_unit_atoms);
  CacheResource nonbond_tb_reserve(nonb_lp.x, small_block_max_atoms);
  ImplicitSolventWorkspace ism_space(poly_ag.getSystemAtomOffsets(),
                                     poly_ag.getSystemAtomCounts(), prec);
  Thermostat heat_bath(poly_ag, ThermostatKind::NONE);
    
  // Upload the synthesis and prime the pumps
  poly_ag.upload();
  poly_ps.upload();
  poly_se.upload();
  mmctrl_fe.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                                  VwuGoal::ACCUMULATE, prec, prec, poly_ag);
  mmctrl_xe.primeWorkUnitCounters(launcher, EvaluateForce::NO, EvaluateEnergy::YES,
                                  VwuGoal::ACCUMULATE, prec, prec, poly_ag);

  // Obtain the appropriate abstracts.  Both precision modes are covered, even though in this test
  // only one is used.
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  const SyValenceKit<double> d_poly_vk = poly_ag.getDoublePrecisionValenceKit(devc);
  const SyValenceKit<float>  f_poly_vk = poly_ag.getSinglePrecisionValenceKit(devc);
  const SyAtomUpdateKit<double, double2, double4> d_poly_auk =
    poly_ag.getDoublePrecisionAtomUpdateKit(devc);
  const SyAtomUpdateKit<float, float2, float4> f_poly_auk =
    poly_ag.getSinglePrecisionAtomUpdateKit(devc);
  const SyNonbondedKit<double, double2> d_poly_nbk = poly_ag.getDoublePrecisionNonbondedKit(devc);
  const SyNonbondedKit<float, float2>  f_poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(devc);
  const SeMaskSynthesisReader poly_ser = poly_se.data(devc);
  const SyRestraintKit<double,
                       double2, double4> d_poly_rk = poly_ag.getDoublePrecisionRestraintKit(devc);
  const SyRestraintKit<float,
                       float2, float4> f_poly_rk = poly_ag.getSinglePrecisionRestraintKit(devc);
  const bool virtual_sites_present = (poly_ag.getVirtualSiteCount() > 0);
  MMControlKit<double> d_ctrl_fe = mmctrl_fe.dpData(devc);
  MMControlKit<double> d_ctrl_xe = mmctrl_xe.dpData(devc);
  MMControlKit<float>  f_ctrl_fe = mmctrl_fe.spData(devc);
  MMControlKit<float>  f_ctrl_xe = mmctrl_xe.spData(devc);
  const CoordinateCycle alt_orientation = getNextCyclePosition(poly_ps.getCyclePosition());
  PsSynthesisWriter poly_psw = poly_ps.data(devc);
  PsSynthesisWriter poly_psw_alt = poly_ps.data(alt_orientation, devc);
  ScoreCardWriter scw = sc.data(devc);
  CacheResourceKit<double> d_vale_fe_tbk = valence_fe_tb_reserve.dpData(devc);
  CacheResourceKit<double> d_vale_xe_tbk = valence_xe_tb_reserve.dpData(devc);
  CacheResourceKit<double> d_nonb_tbk = nonbond_tb_reserve.dpData(devc);
  CacheResourceKit<float> f_vale_fe_tbk = valence_fe_tb_reserve.spData(devc);
  CacheResourceKit<float> f_vale_xe_tbk = valence_xe_tb_reserve.spData(devc);
  CacheResourceKit<float> f_nonb_tbk = nonbond_tb_reserve.spData(devc);
  ISWorkspaceKit<double> d_iswk = ism_space.dpData(devc);
  ISWorkspaceKit<float>  f_iswk = ism_space.spData(devc);
  ReductionKit poly_redk(poly_ag, devc);
  ReductionBridge poly_rbg(poly_ag.getReductionWorkUnitCount());
  ConjGradSubstrate cgsbs(&poly_ps, &poly_rbg, devc);
  LineMinimization line_record(poly_ag.getSystemCount(), mmctrl_fe.getInitialMinimizationStep());
  LinMinWriter lmw = line_record.data(devc);
  ThermostatWriter d_tstw = heat_bath.dpData(devc);
  ThermostatWriter f_tstw = heat_bath.spData(devc);
  
  // Run minimizations
  const int meta_timings = (timer == nullptr) ? 0 : timer->addCategory(test_name);
  const int pert_timings = (timer == nullptr) ? 0 : timer->addCategory("Structure perturbation");
  const int chek_timings = (timer == nullptr) ? 0 : timer->addCategory("Check MM against CPU");
  if (timer != nullptr) {
    timer->assignTime(0);
  }
  const int n_mm_sample = roundUp(mincon.getTotalCycles(), 32) / 32;
  const TrajectoryKind tforce = TrajectoryKind::FORCES; 
  std::vector<double> cpu_total_e(n_mm_sample, 0.0);
  std::vector<double> gpu_total_e(n_mm_sample, 0.0);
  std::vector<double> force_mue(n_mm_sample, 0.0);
  int consistency_failures = 0;
  for (int i = 0; i < mincon.getTotalCycles(); i++) {
    
    // First stage of the cycle: compute forces and obtain the conjugate gradient move.
    poly_ps.initializeForces(gpu, devc);
    ism_space.initialize(devc, CoordinateCycle::WHITE, gpu);
    sc.initialize(devc, gpu);
    switch (prec) {
    case PrecisionModel::DOUBLE:
      if (virtual_sites_present) {
        launchVirtualSitePlacement(&poly_psw_alt, &d_vale_xe_tbk, d_poly_vk, d_poly_auk,
                                   vste_mv_lp);
      }
      launchNonbonded(nb_work_type, d_poly_nbk, poly_ser, &d_ctrl_fe, &poly_psw, &d_tstw, &scw,
                      &d_nonb_tbk, &d_iswk, EvaluateForce::YES, EvaluateEnergy::YES, nonb_lp,
                      gbr_lp, gbd_lp);
      launchValence(d_poly_vk, d_poly_rk, &d_ctrl_fe, &poly_psw, d_poly_auk, &d_tstw, &scw,
                    &d_vale_fe_tbk, EvaluateForce::YES, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                    vale_fe_lp);
      if (virtual_sites_present) {
        launchTransmitVSiteForces(&poly_psw, &d_vale_xe_tbk, d_poly_vk, d_poly_auk, vste_xm_lp);
      }
      d_ctrl_fe.step += 1;
      break;
    case PrecisionModel::SINGLE:
      if (virtual_sites_present) {
        launchVirtualSitePlacement(&poly_psw_alt, &f_vale_xe_tbk, f_poly_vk, f_poly_auk,
                                   vste_mv_lp);
      }
      launchNonbonded(nb_work_type, f_poly_nbk, poly_ser, &f_ctrl_fe, &poly_psw, &f_tstw, &scw,
                      &f_nonb_tbk, &f_iswk, EvaluateForce::YES, EvaluateEnergy::YES,
                      AccumulationMethod::SPLIT, nonb_lp, gbr_lp, gbd_lp);
      launchValence(f_poly_vk, f_poly_rk, &f_ctrl_fe, &poly_psw, f_poly_auk, &f_tstw, &scw,
                    &f_vale_fe_tbk, EvaluateForce::YES, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                    AccumulationMethod::SPLIT, vale_fe_lp);
      if (virtual_sites_present) {
        launchTransmitVSiteForces(&poly_psw, &f_vale_xe_tbk, f_poly_vk, f_poly_auk, vste_xm_lp);
      }
      f_ctrl_fe.step += 1;
      break;
    }
    
    // Check the forces computed for a couple of systems.  This is somewhat redundant, but serves
    // as a sanity check in case other aspects of the energy minimization show problems.
    if (check_mm && (i & 0x1f) == 0) {
      cudaDeviceSynchronize();
      timer->assignTime(meta_timings);
      const int jlim = (3 * i) + 1;
      for (int j = 3 * i; j < jlim; j++) {
        const int jmod = j % poly_ag.getSystemCount();
        PhaseSpace chkj_ps = poly_ps.exportSystem(jmod, devc);
        const std::vector<double> gpu_frc = chkj_ps.getInterlacedCoordinates(tforce);
        chkj_ps.initializeForces();
        ScoreCard tmp_sc(1, 1, 32);
        const AtomGraph *jag_ptr = poly_ag.getSystemTopologyPointer(jmod);
        StaticExclusionMask chkj_se(jag_ptr);
        if (jag_ptr->getImplicitSolventModel() == ImplicitSolventModel::NONE) {
          evalNonbValeMM(&chkj_ps, &tmp_sc, jag_ptr, chkj_se, EvaluateForce::YES, 0);
        }
        else {
          const RestraintApparatus *jra_ptr = poly_ag.getSystemRestraintPointer(jmod);
          evalRestrainedMMGB(&chkj_ps, &tmp_sc, jag_ptr, ngb_tables, chkj_se, jra_ptr,
                             EvaluateForce::YES, 0);
        }
        transmitVirtualSiteForces(&chkj_ps, jag_ptr);
        const std::vector<double> cpu_frc = chkj_ps.getInterlacedCoordinates(tforce);
        cpu_total_e[i / 32] = tmp_sc.reportTotalEnergy(0);
        gpu_total_e[i / 32] = sc.reportTotalEnergy(jmod, devc);
        force_mue[i / 32] = meanUnsignedError(cpu_frc, gpu_frc);
      }
      timer->assignTime(chek_timings);
    }

    // Download and check the forces for each system to verify consistency.  If the forces are
    // consistent enough, set them to be exactly consistent, and do the same with the coordinates,
    // to avoid miniscule roundoff errors that could otherwise creep in over hundreds of steps.
    if (enforce_same_track) {
      poly_ps.download();
      if (checkConsistency(poly_ps, poly_ag, "Force computation", 5.0e-7, 5.0e-7, i)) {
        mandateEquality(&poly_ps, poly_ag);
      }
      else {
        consistency_failures++;
      }
      poly_ps.upload();
    }

    // Perform the conjugate gradient transformation
    if (i == 0) {
      poly_ps.primeConjugateGradientCalculation(gpu, devc);
    }
    switch (prec) {
    case PrecisionModel::DOUBLE:
      launchConjugateGradient(poly_redk, &cgsbs, &d_ctrl_fe, redu_lp);
      break;
    case PrecisionModel::SINGLE:
      launchConjugateGradient(poly_redk, &cgsbs, &f_ctrl_fe, redu_lp);
      break;
    }
    
    // Download and check the conjugate gradient transformation.
    if (enforce_same_track) {
      poly_ps.download();
      if (checkConsistency(poly_ps, poly_ag, "Conjugate gradient transformation",
                           5.0e-7, 5.0e-7, i)) {
        mandateEquality(&poly_ps, poly_ag);
      }
      else {
        consistency_failures++;
      }
      poly_ps.upload();
    }

    // Second stage of the cycle: advance once along the line and recompute the energy.
    launchLineAdvance(prec, &poly_psw, poly_redk, scw, &lmw, 0, redu_lp);
    ism_space.initialize(devc, CoordinateCycle::WHITE, gpu);
    sc.initialize(devc, gpu);
    switch (prec) {
    case PrecisionModel::DOUBLE:
      if (virtual_sites_present) {
        launchVirtualSitePlacement(&poly_psw_alt, &d_vale_xe_tbk, d_poly_vk, d_poly_auk,
                                   vste_mv_lp);
      }
      launchNonbonded(nb_work_type, d_poly_nbk, poly_ser, &d_ctrl_xe, &poly_psw, &d_tstw, &scw,
                      &d_nonb_tbk, &d_iswk, EvaluateForce::NO, EvaluateEnergy::YES, nonb_lp,
                      gbr_lp, gbd_lp);
      launchValence(d_poly_vk, d_poly_rk, &d_ctrl_xe, &poly_psw, d_poly_auk, &d_tstw, &scw,
                    &d_vale_xe_tbk, EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                    vale_xe_lp);
      if (virtual_sites_present) {
        launchTransmitVSiteForces(&poly_psw, &d_vale_xe_tbk, d_poly_vk, d_poly_auk, vste_xm_lp);
      }
      d_ctrl_xe.step += 1;
      break;
    case PrecisionModel::SINGLE:
      if (virtual_sites_present) {
        launchVirtualSitePlacement(&poly_psw_alt, &f_vale_xe_tbk, f_poly_vk, f_poly_auk,
                                   vste_mv_lp);
      }
      launchNonbonded(nb_work_type, f_poly_nbk, poly_ser, &f_ctrl_xe, &poly_psw, &f_tstw, &scw,
                      &f_nonb_tbk, &f_iswk, EvaluateForce::NO, EvaluateEnergy::YES,
                      AccumulationMethod::SPLIT, nonb_lp, gbr_lp, gbd_lp);
      launchValence(f_poly_vk, f_poly_rk, &f_ctrl_xe, &poly_psw, f_poly_auk, &f_tstw, &scw,
                    &f_vale_xe_tbk, EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                    AccumulationMethod::SPLIT, vale_xe_lp);
      if (virtual_sites_present) {
        launchTransmitVSiteForces(&poly_psw, &f_vale_xe_tbk, f_poly_vk, f_poly_auk, vste_xm_lp);
      }
      f_ctrl_xe.step += 1;
      break;
    }

    // Download and check the particle advancement.
    if (enforce_same_track) {
      poly_ps.download();
      if (checkConsistency(poly_ps, poly_ag, "Particle advance I", 1.0e5, 5.0e-7, i)) {
        mandateEquality(&poly_ps, poly_ag);
      }
      else {
        consistency_failures++;
      }
      poly_ps.upload();
    }

    // Third stage of the cycle: advance once more along the line and recompute the energy.
    launchLineAdvance(prec, &poly_psw, poly_redk, scw, &lmw, 1, redu_lp);
    ism_space.initialize(devc, CoordinateCycle::WHITE, gpu);
    sc.initialize(devc, gpu);
    switch (prec) {
    case PrecisionModel::DOUBLE:
      if (virtual_sites_present) {
        launchVirtualSitePlacement(&poly_psw_alt, &d_vale_xe_tbk, d_poly_vk, d_poly_auk,
                                   vste_mv_lp);
      }
      launchNonbonded(nb_work_type, d_poly_nbk, poly_ser, &d_ctrl_xe, &poly_psw, &d_tstw, &scw,
                      &d_nonb_tbk, &d_iswk, EvaluateForce::NO, EvaluateEnergy::YES, nonb_lp,
                      gbr_lp, gbd_lp);
      launchValence(d_poly_vk, d_poly_rk, &d_ctrl_xe, &poly_psw, d_poly_auk, &d_tstw, &scw,
                    &d_vale_xe_tbk, EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                    vale_xe_lp);
      if (virtual_sites_present) {
        launchTransmitVSiteForces(&poly_psw, &d_vale_xe_tbk, d_poly_vk, d_poly_auk, vste_xm_lp);
      }
      d_ctrl_xe.step += 1;
      break;
    case PrecisionModel::SINGLE:
      if (virtual_sites_present) {
        launchVirtualSitePlacement(&poly_psw_alt, &f_vale_xe_tbk, f_poly_vk, f_poly_auk,
                                   vste_mv_lp);
      }
      launchNonbonded(nb_work_type, f_poly_nbk, poly_ser, &f_ctrl_xe, &poly_psw, &f_tstw, &scw,
                      &f_nonb_tbk, &f_iswk, EvaluateForce::NO, EvaluateEnergy::YES,
                      AccumulationMethod::SPLIT, nonb_lp, gbr_lp, gbd_lp);
      launchValence(f_poly_vk, f_poly_rk, &f_ctrl_xe, &poly_psw, f_poly_auk, &f_tstw, &scw,
                    &f_vale_xe_tbk, EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                    AccumulationMethod::SPLIT, vale_xe_lp);
      if (virtual_sites_present) {
        launchTransmitVSiteForces(&poly_psw, &f_vale_xe_tbk, f_poly_vk, f_poly_auk, vste_xm_lp);
      }
      f_ctrl_xe.step += 1;
      break;
    }

    // Download and check the particle advancement.
    if (enforce_same_track) {
      poly_ps.download();
      if (checkConsistency(poly_ps, poly_ag, "Particle advance II", 1.0e5, 5.0e-7, i)) {
        mandateEquality(&poly_ps, poly_ag);
      }
      else {
        consistency_failures++;
      }
      poly_ps.upload();
    }
    
    // Final stage of the cycle: advance a final time along the line and recompute the energy.
    launchLineAdvance(prec, &poly_psw, poly_redk, scw, &lmw, 2, redu_lp);
    ism_space.initialize(devc, CoordinateCycle::WHITE, gpu);
    sc.initialize(devc, gpu);
    switch (prec) {
    case PrecisionModel::DOUBLE:
      if (virtual_sites_present) {
        launchVirtualSitePlacement(&poly_psw_alt, &d_vale_xe_tbk, d_poly_vk, d_poly_auk,
                                   vste_mv_lp);
      }
      launchNonbonded(nb_work_type, d_poly_nbk, poly_ser, &d_ctrl_xe, &poly_psw, &d_tstw, &scw,
                      &d_nonb_tbk, &d_iswk, EvaluateForce::NO, EvaluateEnergy::YES, nonb_lp,
                      gbr_lp, gbd_lp);
      launchValence(d_poly_vk, d_poly_rk, &d_ctrl_xe, &poly_psw, d_poly_auk, &d_tstw, &scw,
                    &d_vale_xe_tbk, EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                    vale_xe_lp);
      if (virtual_sites_present) {
        launchTransmitVSiteForces(&poly_psw, &d_vale_xe_tbk, d_poly_vk, d_poly_auk, vste_xm_lp);
      }
      d_ctrl_xe.step += 1;
      break;
    case PrecisionModel::SINGLE:
      if (virtual_sites_present) {
        launchVirtualSitePlacement(&poly_psw_alt, &f_vale_xe_tbk, f_poly_vk, f_poly_auk,
                                   vste_mv_lp);
      }
      launchNonbonded(nb_work_type, f_poly_nbk, poly_ser, &f_ctrl_xe, &poly_psw, &f_tstw, &scw,
                      &f_nonb_tbk, &f_iswk, EvaluateForce::NO, EvaluateEnergy::YES,
                      AccumulationMethod::SPLIT, nonb_lp, gbr_lp, gbd_lp);
      launchValence(f_poly_vk, f_poly_rk, &f_ctrl_xe, &poly_psw, f_poly_auk, &f_tstw, &scw,
                    &f_vale_xe_tbk, EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                    AccumulationMethod::SPLIT, vale_xe_lp);
      if (virtual_sites_present) {
        launchTransmitVSiteForces(&poly_psw, &f_vale_xe_tbk, f_poly_vk, f_poly_auk, vste_xm_lp);
      }
      f_ctrl_xe.step += 1;
      break;
    }

    // Download and check the particle advancement.
    if (enforce_same_track) {
      poly_ps.download();
      if (checkConsistency(poly_ps, poly_ag, "Particle advance III", 1.0e5, 5.0e-7, i)) {
        mandateEquality(&poly_ps, poly_ag);
      } 
      else {
        consistency_failures++;
      }
      poly_ps.upload();
    }
    
    // Fit a cubic polynomial to guess the best overall advancement, and place the system there.
    launchLineAdvance(prec, &poly_psw, poly_redk, scw, &lmw, 3, redu_lp);
    switch (prec) {
    case PrecisionModel::DOUBLE:
      if (virtual_sites_present) {
        launchVirtualSitePlacement(&poly_psw_alt, &d_vale_xe_tbk, d_poly_vk, d_poly_auk,
                                   vste_mv_lp);
      }
      break;
    case PrecisionModel::SINGLE:
      if (virtual_sites_present) {
        launchVirtualSitePlacement(&poly_psw_alt, &f_vale_xe_tbk, f_poly_vk, f_poly_auk,
                                   vste_mv_lp);
      }
      break;
    }
    
    // Download and check the particle advancement.
    if (enforce_same_track) {
      poly_ps.download();
      if (checkConsistency(poly_ps, poly_ag, "Particle advance IV", 1.0e5, 5.0e-7, i)) {
        mandateEquality(&poly_ps, poly_ag);
      }
      else {
        consistency_failures++;
      }
      poly_ps.upload();
    }
  }
  cudaDeviceSynchronize();
  if (timer != nullptr) {
    timer->assignTime(meta_timings);
  }
  if (enforce_same_track) {
    check(consistency_failures == 0, "Consistency failures occurred during the minimizations.  "
          "The results in identical systems diverged after approximately " +
          std::to_string((consistency_failures + 5) / 6) + " steps.", do_tests);
  }
  if (check_mm) {
    const int ave_atom_count = poly_ag.getAtomCount() / poly_ag.getSystemCount();
    check(cpu_total_e, RelationalOperator::EQUAL, Approx(gpu_total_e).margin(nrg_tol),
          "Energies of relaxed structures did not agree.  Average structure size: " +
          std::to_string(ave_atom_count) + ".  Precision level: " + getEnumerationName(prec) +
          ".  Test name: " + test_name + ".", do_tests);
    check(force_mue, RelationalOperator::EQUAL,
          Approx(std::vector<double>(n_mm_sample, 0.0)).margin(frc_tol), "Snapshots of forces "
          "taken during energy minimization on the GPU do not agree with their CPU-derived "
          "counterparts.  Average structure size: " + std::to_string(ave_atom_count) +
          ".  Precision level: " + getEnumerationName(prec) + ".  Test name: " + test_name +
          ".", do_tests);
  }
  
  // Verify that the energies for all systems meet the expected values
  if (check_mm) {
    const bool snap_exists = (getDrivePathType(snap_name) == DrivePathType::FILE);
    if (snap_exists == false) {
      rtWarn("The snapshot file " + snap_name + " could not be found.  Check the ${STORMM_SOURCE} "
             "environment variable, currently set to " + oe.getStormmSourcePath() + ", for "
             "validity.  Subsequent tests will be skipped.", "test_hpc_minimization");
    }

    // Single-precision minimizations can sometimes end up at different points, depending on the
    // card.  While there seems to be only a handful of different results one can get, it is hard
    // to check in that level of detail.  The results always seem to be self-consistent for the
    // same card, however: the same system goes to the same state, regardless of which SMP it
    // passed through.
    const TestPriority do_snps = (snap_exists && do_tests == TestPriority::CRITICAL) ?
                                 TestPriority::NON_CRITICAL : TestPriority::ABORT;
    const std::vector<double> final_e = sc.reportTotalEnergies(devc);
    const std::string test_var = var_name + ((prec == PrecisionModel::DOUBLE) ? "d" : "f");
    snapshot(snap_name, polyNumericVector(final_e), test_var, 1.0e-5, "Final energies of "
             "energy-minimized structures did not reach their expected values.  Test: " +
             test_name + ".  Precision model: " + getEnumerationName(prec) + ".",
             oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL, psnap, do_snps);
    const int n_unique_systems = maxValue(mol_id_vec) + 1;
    const int nsystems = mol_id_vec.size();
    bool consistent = true;
    std::vector<CombineIDp> failures;
    for (int i = 0; i < n_unique_systems; i++) {
      int ncopies = 0;
      for (int j = 0; j < nsystems; j++) {
        ncopies += (mol_id_vec[j] == i);
      }
      if (ncopies > 0) {
        std::vector<double> copy_e(ncopies);
        ncopies = 0;
        for (int j = 0; j < nsystems; j++) {
          if (mol_id_vec[j] == i) {
            copy_e[ncopies] = final_e[j];
            ncopies++;
          }
        }
        const double nrg_spread = variance(copy_e.data(), ncopies,
                                           VarianceMethod::COEFFICIENT_OF_VARIATION);
        consistent = (consistent && fabs(nrg_spread) < 2.0e-6);
        if (fabs(nrg_spread) >= 2.0e-6) {
          failures.push_back({ i, nrg_spread });
        }
      }
    }
    std::string err_msg("Divergent systems include:");
    if (consistent == false) {
      for (size_t i = 0; i < failures.size(); i++) {
        const int natom = ag_ptr_vec[failures[i].x]->getAtomCount();
        err_msg += " " + std::to_string(failures[i].x) + " (" + std::to_string(natom) +
                   " atoms, " + realToString(failures[i].y, 11, 4, NumberFormat::SCIENTIFIC) +
                   " kcal/mol)";
        if (i < failures.size() - 1LLU) {
          err_msg += ", ";
        }
      }
    }
    check(consistent, "The energies of systems originating at identical coordinates do not "
          "end up at the same values.  " + err_msg, do_tests);
  }

  // Perturb the structures and test the encapsulated energy minimization protocol.  Make a copy
  // as a test against Heisenbugs.
  if (timer != nullptr) {
    timer->assignTime(0);
  }
  Xoshiro256ppGenerator xrs(38175335);
  poly_ps.download();
  PsSynthesisWriter host_psw = poly_ps.data();
  for (int i = 0; i < host_psw.system_count; i++) {
    const int jlim = host_psw.atom_starts[i] + host_psw.atom_counts[i];
    for (int j = host_psw.atom_starts[i]; j < jlim; j++) {
      host_psw.xcrd[j] += 0.1 * (0.5 - xrs.uniformRandomNumber());
    }
  }
  poly_ps.upload();
  PhaseSpaceSynthesis poly_ps_cpy(poly_ps);
  poly_ps_cpy.upload();
  
  if (timer != nullptr) {
    cudaDeviceSynchronize();
    timer->assignTime(pert_timings);
  }
  mincon.setDiagnosticPrintFrequency(mincon.getTotalCycles() / 10);  
  ScoreCard e_refine = launchMinimization(poly_ag, poly_se, &poly_ps, mincon, gpu, prec, 32, timer,
                                          test_name + " (II)");
  cudaDeviceSynchronize();
  e_refine.download();
  e_refine.computePotentialEnergy();
  e_refine.computeTotalEnergy();
  std::vector<std::vector<double>> e_hist_a(poly_ps.getSystemCount());
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    e_hist_a[i] = e_refine.reportHistory(i, HybridTargetLevel::HOST);
  }
  e_refine = launchMinimization(poly_ag, poly_se, &poly_ps_cpy, mincon, gpu, prec, 32, timer,
                                test_name + " (III)");
  cudaDeviceSynchronize();
  e_refine.download();
  e_refine.computePotentialEnergy();
  e_refine.computeTotalEnergy();
  std::vector<std::vector<double>> e_hist_b(poly_ps_cpy.getSystemCount());
  std::vector<int> mismatch(poly_ps_cpy.getSystemCount(), 0);
  const int npts = e_refine.getSampleSize();
  for (int i = 0; i < poly_ps_cpy.getSystemCount(); i++) {
    e_hist_b[i] = e_refine.reportHistory(i, HybridTargetLevel::HOST);
    for (int j = 0; j < npts; j++) {
      if (fabs(e_hist_b[i][j] - e_hist_a[i][j]) > stormm::constants::small) {
        mismatch[i] += 1;
      }
    }
  }
  check(mismatch, RelationalOperator::EQUAL, std::vector<int>(poly_ps.getSystemCount(), 0),
        "Mismatches were detected in back-to-back minimizations of the same batch of perturbed "
        "starting coordinates.  Test: " + test_name + ".  Precision model: " +
        getEnumerationName(prec) + ".", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Read a series of topologies and coordinate files, then process them into a tiled array of
// structures for energy minimization.
//
// Arguments:
//   top_names:  Names of topology files to seek out and read
//   crd_names:  Names of coordinate files to seek out and read
//   tile_list:  Indices of structures to add to the synthesis
//   n_tiles:    The number of times to repeat the tile list when making the synthesis
//   test_name:  Name given to this group of tests
//   oe:         Contains the name of the STORMM source path from shell variables
//   gpu:        Details of the GPU in use
//   test_name:  Name given to this test 
//   psnap:      Instructions as to whether to begin printing a new snapshot file or append to an
//               existing one, if snapshots are to be written
//   snap_name:  Name of the snapshot file for final energies of systems
//   timer:      Time tracking object for optional performance analysis
//-------------------------------------------------------------------------------------------------
void testCompilation(const std::vector<std::string> &top_names,
                     const std::vector<std::string> &crd_names, const std::vector<int> tile_list,
                     const int n_tiles, const double frc_tol, const double nrg_tol,
                     const TestEnvironment &oe, const GpuDetails &gpu,
                     const std::string &test_name, const PrintSituation psnap,
                     const std::string &snap_name, const std::string &var_name, StopWatch *timer,
                     const ImplicitSolventModel ism_choice = ImplicitSolventModel::NONE) {
  const int mol_count = top_names.size();
  if (crd_names.size() != top_names.size()) {
    rtErr("A total of " + std::to_string(top_names.size()) + " topologies and " +
          std::to_string(crd_names.size()) + " coordinate files were provided.  The counts must "
          "match.", "test_hpc_minimization", "testCompilation");
  }
  bool files_exist = true;
  for (int i = 0; i < mol_count; i++) {
    files_exist = (getDrivePathType(top_names[i]) == DrivePathType::FILE &&
                   getDrivePathType(crd_names[i]) == DrivePathType::FILE && files_exist);
  }
  std::vector<AtomGraph> mol_ag;
  std::vector<AtomGraph*> mol_ag_ptr;
  std::vector<PhaseSpace> mol_ps;
  if (files_exist) {
    mol_ag.reserve(mol_count);
    mol_ps.reserve(mol_count);
    mol_ag_ptr.resize(mol_count);
    for (int i = 0; i < mol_count; i++) {
      mol_ag.emplace_back(top_names[i], ExceptionResponse::SILENT);
      mol_ps.emplace_back(crd_names[i]);
      mol_ag_ptr[i] = &mol_ag[i];
    }
  }
  else {
    mol_ag.resize(mol_count);
    mol_ag.resize(mol_count);
    rtWarn("Topology and coordinate files for a number of small molecules and dipeptides were not "
           "found.  Check the ${STORMM_SOURCE} environment variable, currently set to " +
           oe.getStormmSourcePath() + ", for validity.  Subsequent tests will be skipped.",
           "test_hpc_minimization");
  }
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  const int tlen = tile_list.size();
  const int total_mol = tlen * n_tiles;
  std::vector<int> mol_id(total_mol);
  std::vector<int> d_nrg_target(total_mol);
  std::vector<int> f_nrg_target(total_mol);
  
  // The test name determines the content of the target energy vector.  Codify the test name.
  for (int i = 0; i < n_tiles; i++) {
    for (int j = 0; j < tlen; j++) {
      mol_id[(tlen * i) + j] = tile_list[j];
    }
  }
  const PrintSituation x_psnap = (psnap == PrintSituation::OVERWRITE) ? PrintSituation::OVERWRITE :
                                                                        PrintSituation::APPEND;
  metaMinimization(mol_ag_ptr, mol_ps, mol_id, gpu, do_tests, oe, x_psnap, snap_name, var_name,
                   PrecisionModel::DOUBLE, 40, 40, 100, false, true, frc_tol, nrg_tol,
                   test_name + " (fp64)", timer, ism_choice);
  metaMinimization(mol_ag_ptr, mol_ps, mol_id, gpu, do_tests, oe, PrintSituation::APPEND,
                   snap_name, var_name, PrecisionModel::SINGLE, 28, 24, 500, false, true,
                   10.0 * frc_tol, 10.0 * nrg_tol, test_name + " (fp32)", timer, ism_choice);
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {
  TestEnvironment oe(argc, argv);
  StopWatch timer;
  HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  
  // Kernel __shared__ memory configuration
  reductionKernelSetup();
  minimizationKernelSetup();
  
  // Section 1
  section("Minimize a collection of drug molecules and dipeptides");

  // Read small molecules and compile them into a synthesis
  const char osc = osSeparator();
  const std::string base_top_name  = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name  = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string alad_top_name = base_top_name + osc + "ala_dipeptide.top";
  const std::string alad_crd_name = base_crd_name + osc + "ala_dipeptide.inpcrd";
  const std::string brbz_top_name = base_top_name + osc + "bromobenzene_iso.top";
  const std::string brbz_crd_name = base_crd_name + osc + "bromobenzene_iso.inpcrd";
  const std::string lig1_top_name = base_top_name + osc + "stereo_L1.top";
  const std::string lig1_crd_name = base_crd_name + osc + "stereo_L1.inpcrd";
  const std::string lig2_top_name = base_top_name + osc + "symmetry_L1.top";
  const std::string lig2_crd_name = base_crd_name + osc + "symmetry_L1.inpcrd";
  const std::string trpi_top_name = base_top_name + osc + "trpcage.top";
  const std::string trpi_crd_name = base_crd_name + osc + "trpcage.inpcrd";
  const std::string dhfr_top_name = base_top_name + osc + "dhfr_cmap.top";
  const std::string dhfr_crd_name = base_crd_name + osc + "dhfr_cmap.inpcrd";
  const std::string l1vs_top_name = base_top_name + osc + "stereo_L1_vs.top";
  const std::string l1vs_crd_name = base_crd_name + osc + "stereo_L1_vs.inpcrd";
  const std::string l2vs_top_name = base_top_name + osc + "symmetry_L1_vs.top";
  const std::string l2vs_crd_name = base_crd_name + osc + "symmetry_L1_vs.inpcrd";
  const std::string l3vs_top_name = base_top_name + osc + "bromobenzene_vs_iso.top";
  const std::string l3vs_crd_name = base_crd_name + osc + "bromobenzene_vs_iso.inpcrd";
  const std::string l4vs_top_name = base_top_name + osc + "drug_example_vs_iso.top";
  const std::string l4vs_crd_name = base_crd_name + osc + "drug_example_vs_iso.inpcrd";
  const std::vector<std::string> lig_top = { alad_top_name, brbz_top_name, lig1_top_name,
                                             lig2_top_name };
  const std::vector<std::string> lig_crd = { alad_crd_name, brbz_crd_name, lig1_crd_name,
                                             lig2_crd_name };
  const std::vector<std::string> pro_top = { trpi_top_name, dhfr_top_name };
  const std::vector<std::string> pro_crd = { trpi_crd_name, dhfr_crd_name };
  const std::vector<std::string> tvs_top = { l1vs_top_name, l2vs_top_name, trpi_top_name,
                                             l4vs_top_name, brbz_top_name, l3vs_top_name,
                                             lig1_top_name, lig2_top_name };
  const std::vector<std::string> tvs_crd = { l1vs_crd_name, l2vs_crd_name, trpi_crd_name,
                                             l4vs_crd_name, brbz_crd_name, l3vs_crd_name,
                                             lig1_crd_name, lig2_crd_name };  
  const std::string snap_name = oe.getStormmSourcePath() + osc + "test" + osc +
                                "MolecularMechanics" + osc + "min_energy.m";
  // Run small molecule tests
  testCompilation(lig_top, lig_crd, { 0, 1, 2, 3, 0, 1, 2, 3, 0, 3, 1, 2, 2, 1, 3, 0 },
                  256, 1.0e-5, 1.0e-5, oe, gpu, "Small molecules", PrintSituation::OVERWRITE,
                  snap_name, "small_mol_", &timer);

  // Run tests on small proteins
  testCompilation(pro_top, pro_crd, { 0, 1, 0, 1, 1, 1, 0, 0 }, 3, 1.0e-5, 1.0e-3, oe, gpu,
                  "Folded proteins", PrintSituation::APPEND, snap_name, "folded_pro_", &timer);

  // Run tests on small proteins
  testCompilation(pro_top, pro_crd, { 0, 0, 0, 0, 0, 0, 0, 0 }, 8, 1.0e-5, 6.0e-5, oe, gpu,
                  "Trp-cage only", PrintSituation::APPEND, snap_name, "trp_cage_", &timer);

  // Run tests with various GB models
  testCompilation(pro_top, pro_crd, { 0, 0, 0, 0, 0, 0, 0, 0 }, 8, 1.0e-5, 6.0e-5, oe, gpu,
                  "Trp-cage + HCT GB", PrintSituation::APPEND, snap_name, "trp_cage_gb_hct",
                  &timer, ImplicitSolventModel::HCT_GB);

  // Run tests with various GB models
  testCompilation(pro_top, pro_crd, { 0, 0, 0, 0, 0, 0, 0, 0 }, 8, 1.0e-5, 6.0e-5, oe, gpu,
                  "Trp-cage + Neck II GB", PrintSituation::APPEND, snap_name, "trp_cage_gb_nk2",
                  &timer, ImplicitSolventModel::NECK_GB_II);

  // Run tests with various GB models
  testCompilation(pro_top, pro_crd, { 0, 0, 0, 0, 0, 0, 0, 0 }, 8, 1.0e-5, 6.0e-5, oe, gpu,
                  "Trp-cage + OBC GB", PrintSituation::APPEND, snap_name, "trp_cage_gb_obc",
                  &timer, ImplicitSolventModel::OBC_GB);

  // Run tests on larger proteins
  testCompilation(pro_top, pro_crd, { 1, 1, 1, 1, 1, 1, 1, 1 }, 1, 1.0e-5, 6.0e-3, oe, gpu,
                  "DHFR only", PrintSituation::APPEND, snap_name, "dhfr_", &timer);

  // Run tests on systems with virtual sites
  testCompilation(tvs_top, tvs_crd, { 0, 1, 2, 3, 4, 5, 6, 7, 4, 2, 3, 7, 5, 0, 6, 1 }, 4, 1.0e-5,
                  6.0e-3, oe, gpu, "Molecules including some virtual sites",
                  PrintSituation::APPEND, snap_name, "vste_", &timer);  

  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
  return countGlobalTestFailures();
}
