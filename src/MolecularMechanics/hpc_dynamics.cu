// -*-c++-*-
#include "copyright.h"
#include "Constants/behavior.h"
#include "FileManagement/file_enumerators.h"
#include "Math/math_enumerators.h"
#include "Math/rounding.h"
#include "Math/series_ops.h"
#include "Parsing/parsing_enumerators.h"
#include "Potential/hpc_nonbonded_potential.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/trim.h"
#include "Trajectory/hpc_trim.h"
#include "hpc_dynamics.h"
#include "hpc_kinetic.h"

namespace stormm {
namespace mm {

using diskutil::PrintSituation;
using energy::NbwuKind;
using energy::NonbondedTheme;
using energy::launchNonbonded;
using energy::StateVariable;
using energy::ValenceKernelSize;
using parse::NumberFormat;
using stmath::BasisFunctions;
using stmath::FFTMode;
using stmath::incrementingSeries;
using stmath::TableIndexing;
using synthesis::createMaskSynthesis;
using synthesis::ISWorkspaceKit;
using synthesis::SeMaskSynthesisReader;
using topology::UnitCellType;
using trajectory::CoordinateFileKind;
  
//-------------------------------------------------------------------------------------------------
void launchDynamics(const PrecisionModel valence_prec, const PrecisionModel nonbond_prec,
                    const AtomGraphSynthesis &poly_ag, const StaticExclusionMaskSynthesis &poly_se,
                    Thermostat *tst, PhaseSpaceSynthesis *poly_ps, MotionSweeper *mos,
                    const DynamicsControls &dyncon, const ReportControls &repcon,
                    MolecularMechanicsControls *mmctrl_fe, MolecularMechanicsControls *mmctrl_fx,
                    ScoreCard *sc, CacheResource *vale_fe_cache, CacheResource *vale_fx_cache,
                    CacheResource *nonb_fe_cache, CacheResource *nonb_fx_cache,
                    ImplicitSolventWorkspace *ism_space, const AccumulationMethod acc_meth,
                    const SystemCache &sysc, const SynthesisCacheMap &syscmap,
                    const GpuDetails &gpu, const CoreKlManager &launcher, StopWatch *timer,
                    ProgressBar *progress_bar, const std::string &task_name) {

  // Detect the presence of intervention activities: debugging, analysis, perhaps modified force
  // calculations.
  const bool interventions_active = dyna_tk.isActive();
  
  // Extract critical information from the objects.  Begin with a handful of convenient constants.
  const int nstep = dyncon.getStepCount();
  const int ntpr = dyncon.getDiagnosticPrintFrequency();
  const int ntwx = dyncon.getTrajectoryPrintFrequency();
  const int nscm = dyncon.getCenterOfMassMotionPurgeFrequency();
  const NbwuKind nb_work_type = poly_ag.getNonbondedWorkType();
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  const BrokenAsciiCode ascii_recovery = repcon.getAsciiSalvageStyle();

  // Create a place to store coordinates in real-valued format prior to trajectory output.  This
  // is a small amount of memory compared to the allocations for the fixed-precision coordinate
  // synthesis (PhaseSpaceSynthesis) and much smaller than the topology synthesis.
  Condensate staging_zone(poly_ps, nonbond_prec, gpu);
  const int nsys = poly_ps->getSystemCount();
  Hybrid<int2> system_pairs(nsys, "synth_pair_list");
  int2* sysp_ptr = system_pairs.data();
  for (int i = 0; i < nsys; i++) {
    sysp_ptr[i].x = i;
    sysp_ptr[i].y = i;
  }

  // Extract topology abstracts.  Obtain both single- and double-precision variants, as many steps
  // of dynamics will be completed with one or the other.  While one of each will go unused, the
  // cost is trivial.
  const SyNonbondedKit<double, double2> nbk_d = poly_ag.getDoublePrecisionNonbondedKit(devc_tier);
  const SyNonbondedKit<float, float2> nbk_f = poly_ag.getSinglePrecisionNonbondedKit(devc_tier);
  const SyValenceKit<double> vk_d = poly_ag.getDoublePrecisionValenceKit(devc_tier);
  const SyValenceKit<float> vk_f = poly_ag.getSinglePrecisionValenceKit(devc_tier);
  const SyRestraintKit<double,
                       double2,
                       double4_16a> rk_d = poly_ag.getDoublePrecisionRestraintKit(devc_tier);
  const SyRestraintKit<float,
                       float2, float4> rk_f = poly_ag.getSinglePrecisionRestraintKit(devc_tier);
  const SyAtomUpdateKit<double,
                        double2,
                        double4_16a> auk_d = poly_ag.getDoublePrecisionAtomUpdateKit(devc_tier);
  const SyAtomUpdateKit<float,
                        float2,
                        float4> auk_f = poly_ag.getSinglePrecisionAtomUpdateKit(devc_tier);
  const SeMaskSynthesisReader poly_ser = poly_se.data(devc_tier);
  
  // Extract abstracts of the coordinate synthesis oriented towards each point in its time cycle.
  const CoordinateCycle curr_cyc_pos = poly_ps->getCyclePosition();
  const CoordinateCycle next_cyc_pos = getNextCyclePosition(curr_cyc_pos);
  PsSynthesisWriter prm_psw = poly_ps->data(curr_cyc_pos, devc_tier);
  PsSynthesisWriter alt_psw = poly_ps->data(next_cyc_pos, devc_tier);

  // Prepare to remove net center of mass motion with the MotionSweeper at either point in its own
  // time cycle.  While the coordinate synthesis will alternate its current time cycle point with
  // every step, the MotionSweeper will alternate with the frequency of momentum removal (nscm in
  // the user input).
  MotionSweepWriter prm_mosw = mos->data(curr_cyc_pos, devc_tier);
  MotionSweepWriter alt_mosw = mos->data(next_cyc_pos, devc_tier);

  // Obtain abstracts for tracking energies and progress counters.  These are again taken in both
  // single- and double-precision to serve any configuration of the kernels.
  MMControlKit<double> ctrl_fe_d = mmctrl_fe->dpData(devc_tier);
  MMControlKit<double> ctrl_fx_d = mmctrl_fx->dpData(devc_tier);
  MMControlKit<float> ctrl_fe_f = mmctrl_fe->spData(devc_tier);
  MMControlKit<float> ctrl_fx_f = mmctrl_fx->spData(devc_tier);
  ScoreCardWriter scw = sc->data(devc_tier);

  // Obtain abstracts for thread block cache resources.
  CacheResourceKit<double> vale_fe_res_d = vale_fe_cache->dpData(devc_tier);
  CacheResourceKit<double> vale_fx_res_d = vale_fx_cache->dpData(devc_tier);
  CacheResourceKit<double> nonb_fe_res_d = nonb_fe_cache->dpData(devc_tier);
  CacheResourceKit<double> nonb_fx_res_d = nonb_fx_cache->dpData(devc_tier);
  CacheResourceKit<float> vale_fe_res_f = vale_fe_cache->spData(devc_tier);
  CacheResourceKit<float> vale_fx_res_f = vale_fx_cache->spData(devc_tier);
  CacheResourceKit<float> nonb_fe_res_f = nonb_fe_cache->spData(devc_tier);
  CacheResourceKit<float> nonb_fx_res_f = nonb_fx_cache->spData(devc_tier);

  // Obtain abstracts for the implicit solvent workspace.  Its coordinate cycle will be slaved to
  // the coordinate synthesis.
  const ImplicitSolventModel ism_type = poly_ag.getImplicitSolventModel();
  while (ism_space->getCyclePosition() != curr_cyc_pos) {
    ism_space->updateCyclePosition();
  }
  ISWorkspaceKit<double> prm_iswk_d = ism_space->dpData(curr_cyc_pos, devc_tier);
  ISWorkspaceKit<float> prm_iswk_f = ism_space->spData(curr_cyc_pos, devc_tier);
  ISWorkspaceKit<double> alt_iswk_d = ism_space->dpData(next_cyc_pos, devc_tier);
  ISWorkspaceKit<float> alt_iswk_f = ism_space->spData(next_cyc_pos, devc_tier);
  ism_space->initialize(devc_tier, curr_cyc_pos, gpu);

  // Obtain abstracts for the thermostat.
  ThermostatWriter<double> tstw_d = tst->dpData(devc_tier);
  ThermostatWriter<float> tstw_f = tst->spData(devc_tier);
  
  // Get launch parameters for each kernel, the "abstracts" of the kernel manager.  These can be
  // obtained for the specific kernels at the proper precision levels.
  const int2 nonb_bt_fe = launcher.getNonbondedKernelDims(nonbond_prec, nb_work_type,
                                                          EvaluateForce::YES, EvaluateEnergy::YES,
                                                          AccumulationMethod::SPLIT, ism_type,
                                                          ClashResponse::NONE);
  const int2 nonb_bt_fx = launcher.getNonbondedKernelDims(nonbond_prec, nb_work_type,
                                                          EvaluateForce::YES, EvaluateEnergy::NO,
                                                          AccumulationMethod::SPLIT, ism_type,
                                                          ClashResponse::NONE);
  const int2 gbr_bt = launcher.getBornRadiiKernelDims(nonbond_prec, nb_work_type,
                                                      AccumulationMethod::SPLIT, ism_type);
  const int2 gbd_bt = launcher.getBornDerivativeKernelDims(nonbond_prec, nb_work_type,
                                                           AccumulationMethod::SPLIT, ism_type);
  const int2 vale_bt_fe = launcher.getValenceKernelDims(valence_prec, EvaluateForce::YES,
                                                        EvaluateEnergy::YES,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::MOVE_PARTICLES,
                                                        ClashResponse::NONE);
  const int2 vale_bt_fx = launcher.getValenceKernelDims(valence_prec, EvaluateForce::YES,
                                                        EvaluateEnergy::NO,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::MOVE_PARTICLES,
                                                        ClashResponse::NONE);

  // For steps when interventions or checkpointing are needed, the integration process must proceed
  // in stages and paused when complete forces or velocities (at various stages) are available.
  std::vector<IntegrationStage> traj_integration_stages(1, IntegrationStage::VELOCITY_ADVANCE);
  size_t traj_write_stage_idx;
  switch (dyncon.constrainGeometry()) {
  case ApplyConstraints::YES:
    traj_integration_stages.push_back(IntegrationStage::VELOCITY_CONSTRAINT);
    traj_write_stage_idx = 1;
    break;
  case ApplyConstraints::NO:
    traj_write_stage_idx = 0;
    break;
  }
  traj_integration_stages.push_back(IntegrationStage::CALC_KINETIC);
  traj_integration_stages.push_back(IntegrationStage::POSITION_ADVANCE);
  switch (dyncon.constrainGeometry()) {
  case ApplyConstraints::YES:
    traj_integration_stages.push_back(IntegrationStage::GEOMETRY_CONSTRAINT);
    break;
  case ApplyConstraints::NO:
    break;
  }
  const size_t total_traj_intg_stages = traj_integration_stages.size();
  const std::vector<int2> traj_intg_stage_lp(total_traj_intg_stages, vale_bt_fe);

  // Unroll the width of the integration thread blocks, as will be done for the valence
  // interactions kernel in the general launchValence call.
  ValenceKernelSize intg_kwidth;
  if (vale_bt_fe.y > 256) {
    intg_kwidth = ValenceKernelSize::XL;
  }
  else if (vale_bt_fe.y > 128) {
    intg_kwidth = ValenceKernelSize::LG;
  }
  else if (vale_bt_fe.y > 64) {
    intg_kwidth = ValenceKernelSize::MD;
  }
  else {
    intg_kwidth = ValenceKernelSize::SM;
  }
  
  // If there is a valid progress bar, set the total number of reporting cycles.
  const bool show_bar = (progress_bar != nullptr);
  if (show_bar) {
    progress_bar->setTitle(task_name);
    if (ntpr > 0) {
      progress_bar->setCycleCount(nstep / ntpr);
    }
    else {
      progress_bar->setCycleCount(1);
    }
    progress_bar->reset();
  }
  
  // Loop over all steps
  const int traj_freq = dyncon.getTrajectoryPrintFrequency();
  for (int step_idx = 0; step_idx < nstep; step_idx++) {
    PsSynthesisWriter *crd_ptr;
    ISWorkspaceKit<double> *iswk_dptr;
    ISWorkspaceKit<float> *iswk_fptr;
    if (step_idx & 0x1) {
      crd_ptr = &alt_psw;
      iswk_dptr = &alt_iswk_d;
      iswk_fptr = &alt_iswk_f;
    }
    else {
      crd_ptr = &prm_psw;
      iswk_dptr = &prm_iswk_d;
      iswk_fptr = &prm_iswk_f;
    }
    const bool on_energy_step = (step_idx % ntpr == 0);
    const bool on_trajectory_step = (traj_freq > 0 && (step_idx + 1) % traj_freq == 0);
    const bool intervene = (interventions_active && dyna_tk.isActiveOnStep(step_idx));
    
    // Remove any motion of the center of mass, if requested.
    if (nscm > 0 && step_idx > 0 && step_idx % nscm == 0) {
      removeMomentum(poly_ps, poly_ag, mos, gpu);
      mos->updateCyclePosition();
    }
    
    // Initialize energy accumulators if needed.
    if (on_energy_step) {
      sc->initialize(devc_tier, gpu);
    }

    // Prepare to use piecewise kernels or a fused kernel to finish force computations and
    // move particles.
    VwuGoal vale_objective;
    if (on_trajectory_step || intervene) {
      vale_objective = VwuGoal::ACCUMULATE;
    }
    else {
      vale_objective = VwuGoal::MOVE_PARTICLES;
    }
    
    // Perform the non-bonded calculation
    switch (nonbond_prec) {
    case PrecisionModel::DOUBLE:
      if (on_energy_step) {
        launchNonbonded(nb_work_type, nbk_d, poly_ser, &ctrl_fe_d, crd_ptr, &tstw_d, &scw,
                        &nonb_fe_res_d, iswk_dptr, EvaluateForce::YES, EvaluateEnergy::YES,
                        nonb_bt_fe, gbr_bt, gbd_bt, 0.0, 0.0);
      }
      else {
        launchNonbonded(nb_work_type, nbk_d, poly_ser, &ctrl_fx_d, crd_ptr, &tstw_d, &scw,
                        &nonb_fx_res_d, iswk_dptr, EvaluateForce::YES, EvaluateEnergy::NO,
                        nonb_bt_fx, gbr_bt, gbd_bt, 0.0, 0.0);
      }
      break;
    case PrecisionModel::SINGLE:
      if (on_energy_step) {
        launchNonbonded(nb_work_type, nbk_f, poly_ser, &ctrl_fe_f, crd_ptr, &tstw_f, &scw,
                        &nonb_fe_res_f, iswk_fptr, EvaluateForce::YES, EvaluateEnergy::YES,
                        AccumulationMethod::SPLIT, nonb_bt_fe, gbr_bt, gbd_bt, 0.0, 0.0);
      }
      else {
        launchNonbonded(nb_work_type, nbk_f, poly_ser, &ctrl_fx_f, crd_ptr, &tstw_f, &scw,
                        &nonb_fx_res_f, iswk_fptr, EvaluateForce::YES, EvaluateEnergy::NO,
                        AccumulationMethod::SPLIT, nonb_bt_fx, gbr_bt, gbd_bt, 0.0, 0.0);
      }
      break;
    }

    // Perform the valence calculation and atom update.  If interventions or trajectory writing are
    // needed on this step, the VwuGoal parameter (vale_objective) will be modified to have forces
    // written and save other steps for subsequent kernels.
    switch (valence_prec) {
    case PrecisionModel::DOUBLE:
      if (on_energy_step) {
        launchValence(vk_d, rk_d, &ctrl_fe_d, crd_ptr, auk_d, &tstw_d, &scw, &vale_fe_res_d,
                      EvaluateForce::YES, EvaluateEnergy::YES, vale_objective, vale_bt_fe, 0.0,
                      0.0);
      }
      else {
        launchValence(vk_d, rk_d, &ctrl_fx_d, crd_ptr, auk_d, &tstw_d, &scw, &vale_fx_res_d,
                      EvaluateForce::YES, EvaluateEnergy::NO, vale_objective, vale_bt_fx, 0.0,
                      0.0);
      }
      break;
    case PrecisionModel::SINGLE:
      if (on_energy_step) {
        launchValence(vk_f, rk_f, &ctrl_fe_f, crd_ptr, auk_f, &tstw_f, &scw, &vale_fe_res_f,
                      EvaluateForce::YES, EvaluateEnergy::YES, vale_objective,
                      AccumulationMethod::SPLIT, vale_bt_fe, 0.0, 0.0);
      }
      else {
        launchValence(vk_f, rk_f, &ctrl_fx_f, crd_ptr, auk_f, &tstw_f, &scw, &vale_fx_res_f,
                      EvaluateForce::YES, EvaluateEnergy::NO, vale_objective,
                      AccumulationMethod::SPLIT, vale_bt_fx, 0.0, 0.0);
      }
      break;
    }
    if (intervene) {
      dyna_tk.execute(step_idx, IntegrationStage::CALC_FORCES);
    }

    // Finish the dynamics step if necessary, writing trajectory components or performing
    // interventions.
    if (intervene || on_trajectory_step) {
      for (size_t i = 0; i < total_traj_intg_stages; i++) {
        switch (valence_prec) {
        case PrecisionModel::DOUBLE:
          if (on_energy_step) {
            launchIntegrationProcess(crd_ptr, &vale_fe_res_d, &ctrl_fe_d, &scw, vk_d, auk_d,
                                     tstw_d, traj_intg_stage_lp[i], traj_integration_stages[i]);
          }
          else {
            launchIntegrationProcess(crd_ptr, &vale_fe_res_d, &ctrl_fx_d, vk_d, auk_d, tstw_d,
                                     traj_intg_stage_lp[i], traj_integration_stages[i]);
          }
          break;
        case PrecisionModel::SINGLE:
          if (on_energy_step) {
            launchIntegrationProcess(crd_ptr, &vale_fe_res_f, &ctrl_fe_f, &scw, vk_f, auk_f,
                                     tstw_f, traj_intg_stage_lp[i], AccumulationMethod::SPLIT,
                                     intg_kwidth, traj_integration_stages[i]);
          }
          else {

            // Energy calculations must not be performed on steps where the energy is not
            // requested, even if the comprehensive list of piecewise integration procedures
            // contains a calculation of the kinetic energy.
            if (traj_integration_stages[i] != IntegrationStage::CALC_KINETIC) {
              launchIntegrationProcess(crd_ptr, &vale_fe_res_f, &ctrl_fx_f, vk_f, auk_f, tstw_f,
                                       traj_intg_stage_lp[i], AccumulationMethod::SPLIT,
                                       intg_kwidth, traj_integration_stages[i]);
            }
          }
          break;
        }         
        if (on_trajectory_step && i == traj_write_stage_idx) {
          writeSynthesisFrames(*poly_ps, &staging_zone, system_pairs, syscmap,
                               CoordinateFileKind::AMBER_CRD,
                               static_cast<double>(step_idx) * tstw_d.dt,
                               TrajectoryKind::POSITIONS, HybridTargetLevel::DEVICE, gpu,
                               ascii_recovery);
        }
        if (intervene) {
          dyna_tk.execute(step_idx, traj_integration_stages[i]);
        }
      }

      // Clear forces in the alternate buffer so that fresh accumulators will be ready for the next
      // cycle.  This is done by the non-bonded routine (vacuum conditions) or by the Generalized
      // Born radii computation kernel, but the piecewise velocity updates make use of these
      // accumulators to store the Langevin impulses on particles as well as other total force
      // accumulations.
      poly_ps->initializeForces(getNextCyclePosition(poly_ps->getCyclePosition()), gpu,
                                HybridTargetLevel::DEVICE);
    }

    // Log energies if requested.  Refresh work unit progress counters.
    if (on_energy_step) {
      ctrl_fe_d.step += 1;
      ctrl_fe_f.step += 1;
      mmctrl_fe->incrementStep();
      launchTemperatureComputation(poly_ag, &scw, tstw_d.cnst_geom, gpu);
      sc->commit(devc_tier, gpu);
      sc->incrementSampleCount();
      sc->setLastTimeStep(tstw_d.step, HybridTargetLevel::DEVICE);
      if (show_bar) {
        progress_bar->update();
      }
    }
    else {
      ctrl_fx_d.step += 1;
      ctrl_fx_f.step += 1;
      mmctrl_fx->incrementStep();
    }
    
    // Increment the cycle positions.  The thermostat is the official keeper of the time step.
    // Molecular mechanics progress trackers mmctrl_f{e,x} must increment only to ensure that they
    // refresh their counters when appropriate.
    poly_ps->updateCyclePosition();
    tst->incrementStep();
    tstw_d.step += 1;
    tstw_f.step += 1;
  }
  sc->computePotentialEnergy(HybridTargetLevel::DEVICE, gpu);
  sc->computeTotalEnergy(HybridTargetLevel::DEVICE, gpu);
  if (show_bar)	{
    progress_bar->finalizeTerminalOutput();
  }
}

//-------------------------------------------------------------------------------------------------
ScoreCard launchDynamics(const AtomGraphSynthesis &poly_ag,
                         const StaticExclusionMaskSynthesis &poly_se, Thermostat *tst,
                         PhaseSpaceSynthesis *poly_ps, const DynamicsControls &dyncon,
                         const ReportControls &repcon, const SystemCache &sysc,
                         const SynthesisCacheMap &syscmap, const GpuDetails &gpu,
                         const PrecisionModel valence_prec, const PrecisionModel nonbond_prec,
                         const int energy_bits, StopWatch *timer, ProgressBar *progress_bar,
                         const std::string &task_name) {
  
  // The thermostat will have been initialized before being submitted to this function.  Create
  // the energy tracking object.
  const int ntpr   = dyncon.getDiagnosticPrintFrequency();
  const int nframe = (roundUp(dyncon.getStepCount(), ntpr) / ntpr) + 1;
  ScoreCard result(poly_ps->getSystemCount(), nframe, energy_bits);
  MotionSweeper mos(poly_ps);
  mos.uploadAll();

  // Create the kernel launcher based on the GPU and workload.  Use launch parameters to make
  // appropriate allocations of memory resources for thread blocks and molecular mechanics
  // progress trackers.
  const CoreKlManager launcher(gpu, poly_ag);
  const int2 vale_fe_lp = launcher.getValenceKernelDims(valence_prec, EvaluateForce::YES,
                                                        EvaluateEnergy::YES,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::MOVE_PARTICLES,
                                                        ClashResponse::NONE);
  const int2 vale_fx_lp = launcher.getValenceKernelDims(valence_prec, EvaluateForce::YES,
                                                        EvaluateEnergy::NO,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::MOVE_PARTICLES,
                                                        ClashResponse::NONE);
  const NbwuKind nb_work_type = poly_ag.getNonbondedWorkType();
  const ImplicitSolventModel ism_type = poly_ag.getImplicitSolventModel();
  const int2 nonb_fe_lp = launcher.getNonbondedKernelDims(nonbond_prec, nb_work_type,
                                                          EvaluateForce::YES, EvaluateEnergy::YES,
                                                          AccumulationMethod::SPLIT, ism_type,
                                                          ClashResponse::NONE);
  const int2 nonb_fx_lp = launcher.getNonbondedKernelDims(nonbond_prec, nb_work_type,
                                                          EvaluateForce::YES, EvaluateEnergy::NO,
                                                          AccumulationMethod::SPLIT, ism_type,
                                                          ClashResponse::NONE);
  MolecularMechanicsControls mmctrl_fe(dyncon);
  MolecularMechanicsControls mmctrl_fx(dyncon);
  switch (poly_ag.getUnitCellType()) {
  case UnitCellType::NONE:
    mmctrl_fe.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                                    ClashResponse::NONE, VwuGoal::MOVE_PARTICLES, valence_prec,
                                    nonbond_prec, poly_ag);
    mmctrl_fx.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::NO,
                                    ClashResponse::NONE, VwuGoal::MOVE_PARTICLES, valence_prec,
                                    nonbond_prec, poly_ag);
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    break;
  }
  CacheResource vale_fe_cache(vale_fe_lp.x, maximum_valence_work_unit_atoms);
  CacheResource vale_fx_cache(vale_fx_lp.x, maximum_valence_work_unit_atoms);
  CacheResource nonb_fe_cache(nonb_fe_lp.x, small_block_max_atoms);
  CacheResource nonb_fx_cache(nonb_fx_lp.x, small_block_max_atoms);
  
  // Create an implicit solvent workspace for this system.
  ImplicitSolventWorkspace ism_space(poly_ag.getSystemAtomOffsets(),
                                     poly_ag.getSystemAtomCounts(), nonbond_prec);
  launchDynamics(valence_prec, nonbond_prec, poly_ag, poly_se, tst, poly_ps, &mos, dyncon, repcon,
                 &mmctrl_fe, &mmctrl_fx, &result, &vale_fe_cache, &vale_fx_cache, &nonb_fe_cache,
                 &nonb_fx_cache, &ism_space, AccumulationMethod::SPLIT, sysc, syscmap, gpu,
                 launcher, timer, progress_bar, task_name);
  return result;
}

//-------------------------------------------------------------------------------------------------
void launchDynamics(const AtomGraphSynthesis &poly_ag, const LocalExclusionMask &lem,
                    const PPITable &nrg_tab, PhaseSpaceSynthesis *poly_ps,
                    Condensate *staging_zone, CellGrid<double, llint, double, double4_16a> *cg,
                    PMIGrid *pmig, ConvolutionManager *cvol, MotionSweeper *mos, Thermostat *tst,
                    ScoreCard *sc, MolecularMechanicsControls *mmctrl_fe,
                    MolecularMechanicsControls *mmctrl_fx, CacheResource *vale_fe_cache,
                    CacheResource *vale_fx_cache, TileManager *tlmn_fe, TileManager *tlmn_fx,
                    const DynamicsControls &dyncon, const PPPMControls &pmecon,
                    const ReportControls &repcon, const SystemCache &sysc,
                    const SynthesisCacheMap &syscmap, const GpuDetails &gpu,
                    const CoreKlManager &launcher, const PrecisionModel valence_prec,
                    StopWatch *timer, ProgressBar *progress_bar, const std::string &task_name) {
  launchDynamics<double, llint,
                 double, double4_16a>(poly_ps, staging_zone, cg, pmig, cvol, mos, tst, sc,
                                      mmctrl_fe, mmctrl_fx, vale_fe_cache, vale_fx_cache, tlmn_fe,
                                      tlmn_fx, poly_ag, lem, nrg_tab, dyncon, pmecon, repcon, sysc,
                                      syscmap, gpu, launcher, valence_prec, timer, progress_bar,
                                      task_name);
}

//-------------------------------------------------------------------------------------------------
void launchDynamics(const AtomGraphSynthesis &poly_ag, const LocalExclusionMask &lem,
                    const PPITable &nrg_tab, PhaseSpaceSynthesis *poly_ps,
                    Condensate *staging_zone, CellGrid<double, llint, float, double4_16a> *cg,
                    PMIGrid *pmig, ConvolutionManager *cvol, MotionSweeper *mos, Thermostat *tst,
                    ScoreCard *sc, MolecularMechanicsControls *mmctrl_fe,
                    MolecularMechanicsControls *mmctrl_fx, CacheResource *vale_fe_cache,
                    CacheResource *vale_fx_cache, TileManager *tlmn_fe, TileManager *tlmn_fx,
                    const DynamicsControls &dyncon, const PPPMControls &pmecon,
                    const ReportControls &repcon, const SystemCache &sysc,
                    const SynthesisCacheMap &syscmap, const GpuDetails &gpu,
                    const CoreKlManager &launcher, const PrecisionModel valence_prec,
                    StopWatch *timer, ProgressBar *progress_bar, const std::string &task_name) {
  launchDynamics<double, llint,
                 float, double4_16a>(poly_ps, staging_zone, cg, pmig, cvol, mos, tst, sc,
                                     mmctrl_fe, mmctrl_fx, vale_fe_cache, vale_fx_cache, tlmn_fe,
                                     tlmn_fx, poly_ag, lem, nrg_tab, dyncon, pmecon, repcon, sysc,
                                     syscmap, gpu, launcher, valence_prec, timer, progress_bar,
                                     task_name);
}

//-------------------------------------------------------------------------------------------------
void launchDynamics(const AtomGraphSynthesis &poly_ag, const LocalExclusionMask &lem,
                    const PPITable &nrg_tab, PhaseSpaceSynthesis *poly_ps,
                    Condensate *staging_zone, CellGrid<float, int, double, float4> *cg,
                    PMIGrid *pmig, ConvolutionManager *cvol, MotionSweeper *mos, Thermostat *tst,
                    ScoreCard *sc, MolecularMechanicsControls *mmctrl_fe,
                    MolecularMechanicsControls *mmctrl_fx, CacheResource *vale_fe_cache,
                    CacheResource *vale_fx_cache, TileManager *tlmn_fe, TileManager *tlmn_fx,
                    const DynamicsControls &dyncon, const PPPMControls &pmecon,
                    const ReportControls &repcon, const SystemCache &sysc,
                    const SynthesisCacheMap &syscmap, const GpuDetails &gpu,
                    const CoreKlManager &launcher, const PrecisionModel valence_prec,
                    StopWatch *timer, ProgressBar *progress_bar, const std::string &task_name) {
  launchDynamics<float, int, double, float4>(poly_ps, staging_zone, cg, pmig, cvol, mos, tst, sc,
                                             mmctrl_fe, mmctrl_fx, vale_fe_cache, vale_fx_cache,
                                             tlmn_fe, tlmn_fx, poly_ag, lem, nrg_tab, dyncon,
                                             pmecon, repcon, sysc, syscmap, gpu, launcher,
                                             valence_prec, timer, progress_bar, task_name);
}

//-------------------------------------------------------------------------------------------------
void launchDynamics(const AtomGraphSynthesis &poly_ag, const LocalExclusionMask &lem,
                    const PPITable &nrg_tab, PhaseSpaceSynthesis *poly_ps,
                    Condensate *staging_zone, CellGrid<float, int, float, float4> *cg,
                    PMIGrid *pmig, ConvolutionManager *cvol, MotionSweeper *mos, Thermostat *tst,
                    ScoreCard *sc, MolecularMechanicsControls *mmctrl_fe,
                    MolecularMechanicsControls *mmctrl_fx, CacheResource *vale_fe_cache,
                    CacheResource *vale_fx_cache, TileManager *tlmn_fe, TileManager *tlmn_fx,
                    const DynamicsControls &dyncon, const PPPMControls &pmecon,
                    const ReportControls &repcon, const SystemCache &sysc,
                    const SynthesisCacheMap &syscmap, const GpuDetails &gpu,
                    const CoreKlManager &launcher, const PrecisionModel valence_prec,
                    StopWatch *timer, ProgressBar *progress_bar, const std::string &task_name) {
  launchDynamics<float, int, float, float4>(poly_ps, staging_zone, cg, pmig, cvol, mos, tst, sc,
                                            mmctrl_fe, mmctrl_fx, vale_fe_cache, vale_fx_cache,
                                            tlmn_fe, tlmn_fx, poly_ag, lem, nrg_tab, dyncon,
                                            pmecon, repcon, sysc, syscmap, gpu, launcher,
                                            valence_prec, timer, progress_bar, task_name);
}

//-------------------------------------------------------------------------------------------------
ScoreCard launchDynamics(const AtomGraphSynthesis &poly_ag, PhaseSpaceSynthesis *poly_ps,
                         const DynamicsControls &dyncon, const PPPMControls &pmecon,
                         const PrecisionControls &preccon, const ReportControls &repcon,
                         const SystemCache &sysc, const SynthesisCacheMap &syscmap,
                         const GpuDetails &gpu, StopWatch *timer, ProgressBar *progress_bar,
                         const std::string &task_name) {
  const PrecisionModel nonbond_prec = preccon.getNonbondedMethod();
  const PrecisionModel valence_prec = preccon.getValenceMethod();

  // Create a place to store coordinates in real-valued format prior to trajectory output.  This
  // is a small amount of memory compared to the allocations for the fixed-precision coordinate
  // synthesis (PhaseSpaceSynthesis) and much smaller than the topology synthesis.
  Condensate staging_zone(poly_ps, nonbond_prec, gpu);
  
  // This overloaded variant of launchDynamics will create the cell grid neighbor list, thermostat,
  // and all other resources needed by the coordinate and topology syntheses based on the minimal
  // user input data.  While it is not as efficient as the more differentiated variants, which
  // include pre-allocated resources, if called many times, this variant will be preferred by
  // developers who want to minimize complexity at a high level.
  MolecularMechanicsControls mmctrl_fe(dyncon, pmecon);
  MolecularMechanicsControls mmctrl_fx(dyncon, pmecon);
  if (fabs(mmctrl_fe.getLongestCutoff() - mmctrl_fx.getLongestCutoff()) > constants::small) {
    rtErr("The force-only and force+energy molecular mechanics control objects disagree in terms "
          "of the longest cutoff (" +
          realToString(mmctrl_fx.getLongestCutoff(), 9, 5, NumberFormat::STANDARD_REAL) + ", " +
          realToString(mmctrl_fe.getLongestCutoff(), 9, 5, NumberFormat::STANDARD_REAL) + ").",
          "launchDynamics");
  }

  // Having created each mmctrl object with both dyncon and pmecon, either of which might contain
  // user input as to the cutoff, it is now the "source of truth." Both variables will contain the
  // same results for the single, unified cutoff that the simulation will rely upon, so take
  // mmctrl_fe as representative.
  LocalExclusionMask lem(poly_ag);
  int log_tab_bits;
  switch (preccon.getNonbondedMethod()) {
  case PrecisionModel::DOUBLE:
    log_tab_bits = 6;
    break;
  case PrecisionModel::SINGLE:
    log_tab_bits = 5;
    break;
  }
  PPITable nrg_tab(NonbondedTheme::ELECTROSTATIC, BasisFunctions::MIXED_FRACTIONS,
                   TableIndexing::SQUARED_ARG, mmctrl_fe.getElectrostaticCutoff(), 0.0,
                   pmecon.getDirectSumTolerance(), log_tab_bits);
  Thermostat tst(poly_ag, dyncon, sysc, syscmap.getCacheOrigins());
  MotionSweeper mos(poly_ps, preccon.getMomentumConservationBits(), preccon.getCenterOfMassBits());
  const CoreKlManager launcher(gpu, poly_ag);
  switch (poly_ag.getUnitCellType()) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    mmctrl_fe.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                                    ClashResponse::NONE, VwuGoal::MOVE_PARTICLES, valence_prec,
                                    nonbond_prec, poly_ag);
    mmctrl_fx.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::NO,
                                    ClashResponse::NONE, VwuGoal::MOVE_PARTICLES, valence_prec,
                                    nonbond_prec, poly_ag);
    break;
  }
  const int2 vale_fe_lp = launcher.getValenceKernelDims(valence_prec, EvaluateForce::YES,
                                                        EvaluateEnergy::YES,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::MOVE_PARTICLES,
                                                        ClashResponse::NONE);
  const int2 vale_fx_lp = launcher.getValenceKernelDims(valence_prec, EvaluateForce::YES,
                                                        EvaluateEnergy::NO,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::MOVE_PARTICLES,
                                                        ClashResponse::NONE);
  CacheResource vale_fe_cache(vale_fe_lp.x, maximum_valence_work_unit_atoms);
  CacheResource vale_fx_cache(vale_fx_lp.x, maximum_valence_work_unit_atoms);
  
  // Upload the new supporting objects.  The coordinate and topology syntheses are assumed to have
  // been staged on the GPU prior to calling this function.  The workflow then merges with a
  // templated function called on the HPC-compiled side.  The templated function will call various
  // non-templated kernels for valence and nonbonded interactions as well as particle migration.
  lem.upload();
  nrg_tab.upload();
  tst.uploadPartitions();
  mos.uploadAll();
  const int ntpr   = dyncon.getDiagnosticPrintFrequency();
  const int nframe = (roundUp(dyncon.getStepCount(), ntpr) / ntpr) + 1;
  ScoreCard result(poly_ps->getSystemCount(), nframe, preccon.getEnergyScalingBits());
  switch (nonbond_prec) {
  case PrecisionModel::DOUBLE:
    {
      CellGrid<double, llint,
               double, double4_16a> cg(poly_ps, poly_ag, mmctrl_fe.getLongestCutoff(), 0.02,
                                       pmecon.getMeshSubdivisions(), NonbondedTheme::ALL);
      std::vector<CellGrid<double, llint, double, double4_16a>> cg_workspace;
      if (dyna_tk.doNeighborListChecks()) {
        cg_workspace.reserve(1);
        cg_workspace.emplace_back(dyna_tk.getWorkspacePointer(), poly_ag,
                                  mmctrl_fe.getLongestCutoff(), 0.02,
                                  pmecon.getMeshSubdivisions(), NonbondedTheme::ALL);
        dyna_tk.setNeighborList(&cg);
        dyna_tk.setNLWorkspace(&cg_workspace[0]);
      }
      cg.checkViability();
      mmctrl_fe.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                                      ClashResponse::NONE, VwuGoal::MOVE_PARTICLES, valence_prec,
                                      nonbond_prec, QMapMethod::ACC_SHARED, nonbond_prec,
                                      double_type_index, pmecon.getInterpolationOrder(),
                                      NeighborListKind::MONO, cg.getTinyBoxPresence(), poly_ag);
      mmctrl_fx.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::NO,
                                      ClashResponse::NONE, VwuGoal::MOVE_PARTICLES, valence_prec,
                                      nonbond_prec, QMapMethod::ACC_SHARED, nonbond_prec,
                                      double_type_index, pmecon.getInterpolationOrder(),
                                      NeighborListKind::MONO, cg.getTinyBoxPresence(), poly_ag);
      const int2 pair_fe_lp = launcher.getPMEPairsKernelDims(nonbond_prec, nonbond_prec,
                                                             NeighborListKind::MONO,
                                                             cg.getTinyBoxPresence(),
                                                             EvaluateForce::YES,
                                                             EvaluateEnergy::YES,
                                                             ClashResponse::NONE);
      const int2 pair_fx_lp = launcher.getPMEPairsKernelDims(nonbond_prec, nonbond_prec,
                                                             NeighborListKind::MONO,
                                                             cg.getTinyBoxPresence(),
                                                             EvaluateForce::YES,
                                                             EvaluateEnergy::NO,
                                                             ClashResponse::NONE);
      TileManager tlmn_fe(pair_fe_lp);
      TileManager tlmn_fx(pair_fx_lp);
      PMIGrid pmig(&cg, NonbondedTheme::ELECTROSTATIC, pmecon.getInterpolationOrder(),
                   nonbond_prec, FFTMode::OUT_OF_PLACE, preccon.getChargeMeshScalingBits(),
                   preccon.getChargeMeshScalingBits(), gpu, pmecon.getDensityMappingMethod());
      ConvolutionManager cvol(pmig, pmecon.getEwaldCoefficient(), gpu);
      cg.upload();
      pmig.upload();
      cvol.upload();
      mmctrl_fe.upload();
      mmctrl_fx.upload();
      tlmn_fe.upload();
      tlmn_fx.upload();
      launchDynamics<double, llint,
                     double, double4_16a>(poly_ps, &staging_zone, &cg, &pmig, &cvol, &mos, &tst,
                                          &result, &mmctrl_fe, &mmctrl_fx, &vale_fe_cache,
                                          &vale_fx_cache, &tlmn_fe, &tlmn_fx, poly_ag, lem,
                                          nrg_tab, dyncon, pmecon, repcon, sysc, syscmap, gpu,
                                          launcher, valence_prec, timer, progress_bar, task_name);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      CellGrid<float, int, float, float4> cg(poly_ps, poly_ag, mmctrl_fe.getLongestCutoff(), 0.02,
                                             pmecon.getMeshSubdivisions(), NonbondedTheme::ALL);
      std::vector<CellGrid<float, int, float, float4>> cg_workspace;
      if (dyna_tk.doNeighborListChecks()) {
        cg_workspace.reserve(1);
        cg_workspace.emplace_back(dyna_tk.getWorkspacePointer(), poly_ag,
                                  mmctrl_fe.getLongestCutoff(), 0.02,
                                  pmecon.getMeshSubdivisions(), NonbondedTheme::ALL);
        dyna_tk.setNeighborList(&cg);
        dyna_tk.setNLWorkspace(&cg_workspace[0]);
      }
      cg.checkViability();
      mmctrl_fe.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                                      ClashResponse::NONE, VwuGoal::MOVE_PARTICLES, valence_prec,
                                      nonbond_prec, QMapMethod::ACC_SHARED, nonbond_prec,
                                      double_type_index, pmecon.getInterpolationOrder(),
                                      NeighborListKind::MONO, cg.getTinyBoxPresence(), poly_ag);
      mmctrl_fx.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::NO,
                                      ClashResponse::NONE, VwuGoal::MOVE_PARTICLES, valence_prec,
                                      nonbond_prec, QMapMethod::ACC_SHARED, nonbond_prec,
                                      double_type_index, pmecon.getInterpolationOrder(),
                                      NeighborListKind::MONO, cg.getTinyBoxPresence(), poly_ag);
      const int2 pair_fe_lp = launcher.getPMEPairsKernelDims(nonbond_prec, nonbond_prec,
                                                             NeighborListKind::MONO,
                                                             cg.getTinyBoxPresence(),
                                                             EvaluateForce::YES,
                                                             EvaluateEnergy::YES,
                                                             ClashResponse::NONE);
      const int2 pair_fx_lp = launcher.getPMEPairsKernelDims(nonbond_prec, nonbond_prec,
                                                             NeighborListKind::MONO,
                                                             cg.getTinyBoxPresence(),
                                                             EvaluateForce::YES,
                                                             EvaluateEnergy::NO,
                                                             ClashResponse::NONE);
      TileManager tlmn_fe(pair_fe_lp);
      TileManager tlmn_fx(pair_fx_lp);
      PMIGrid pmig(&cg, NonbondedTheme::ELECTROSTATIC, pmecon.getInterpolationOrder(),
                   nonbond_prec, FFTMode::OUT_OF_PLACE, preccon.getChargeMeshScalingBits(),
                   preccon.getChargeMeshScalingBits(), gpu, pmecon.getDensityMappingMethod());
      ConvolutionManager cvol(pmig, pmecon.getEwaldCoefficient(), gpu);
      cg.upload();
      pmig.upload();
      cvol.upload();
      mmctrl_fe.upload();
      mmctrl_fx.upload();
      tlmn_fe.upload();
      tlmn_fx.upload();
      launchDynamics<float, int,
                     float, float4>(poly_ps, &staging_zone, &cg, &pmig, &cvol, &mos, &tst, &result,
                                    &mmctrl_fe, &mmctrl_fx, &vale_fe_cache, &vale_fx_cache,
                                    &tlmn_fe, &tlmn_fx, poly_ag, lem, nrg_tab, dyncon, pmecon,
                                    repcon, sysc, syscmap, gpu, launcher, valence_prec, timer,
                                    progress_bar, task_name);
    }
    break;
  }
  return result;
}

} // namespace mm
} // namespace stormm
