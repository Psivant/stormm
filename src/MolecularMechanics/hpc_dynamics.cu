// -*-c++-*-
#include "copyright.h"
#include "FileManagement/file_enumerators.h"
#include "Math/rounding.h"
#include "Math/series_ops.h"
#include "Potential/hpc_nonbonded_potential.h"
#include "Potential/hpc_valence_potential.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/trajectory_enumerators.h"
#include "Trajectory/trim.h"
#include "hpc_dynamics.h"
#include "hpc_kinetic.h"

namespace stormm {
namespace mm {

using diskutil::PrintSituation;
using energy::CacheResourceKit;
using energy::ClashResponse;
using energy::EvaluateEnergy;
using energy::EvaluateForce;
using energy:: NbwuKind;
using energy::launchNonbonded;
using energy::ScoreCardWriter;
using energy::StateVariable;
using synthesis::createMaskSynthesis;
using synthesis::ISWorkspaceKit;
using synthesis::PsSynthesisWriter;
using synthesis::SyAtomUpdateKit;
using synthesis::SyNonbondedKit;
using synthesis::SyRestraintKit;
using synthesis::SyValenceKit;
using synthesis::SeMaskSynthesisReader;
using synthesis::VwuGoal;
using stmath::roundUp;
using stmath::incrementingSeries;
using topology::UnitCellType;
using trajectory::CoordinateCycle;
using trajectory::CoordinateFileKind;
using trajectory::getNextCyclePosition;
using trajectory::MotionSweepWriter;
using trajectory::removeMomentum;
using trajectory::ThermostatWriter;
  
//-------------------------------------------------------------------------------------------------
void launchDynamics(const PrecisionModel valence_prec, const PrecisionModel nonbond_prec,
                    const AtomGraphSynthesis &poly_ag, const StaticExclusionMaskSynthesis &poly_se,
                    Thermostat *tst, PhaseSpaceSynthesis *poly_ps, MotionSweeper *mos,
                    const DynamicsControls &dyncon, MolecularMechanicsControls *mmctrl_fe,
                    MolecularMechanicsControls *mmctrl_fx, ScoreCard *sc,
                    CacheResource *vale_fe_cache, CacheResource *vale_fx_cache,
                    CacheResource *nonb_fe_cache, CacheResource *nonb_fx_cache,
                    ImplicitSolventWorkspace *ism_space, const AccumulationMethod acc_meth,
                    const SystemCache &sysc, const SynthesisCacheMap &syscmap,
                    const GpuDetails &gpu, const CoreKlManager &launcher, StopWatch *timer,
                    const std::string &task_name) {

  // Extract critical information from the objects.  Begin with a handful of convenient constants.
  const int nstep = dyncon.getStepCount();
  const int ntpr = dyncon.getDiagnosticPrintFrequency();
  const int ntwx = dyncon.getTrajectoryPrintFrequency();
  const int nscm = dyncon.getCenterOfMassMotionPurgeFrequency();
  const NbwuKind nb_work_type = poly_ag.getNonbondedWorkType();
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;

  // Extract topology abstracts.  Obtain both single- and double-precision variants, as many steps
  // of dynamics will be completed with one or the other.  While one of each will go unused, the
  // cost is trivial.
  const SyNonbondedKit<double, double2> nbk_d = poly_ag.getDoublePrecisionNonbondedKit(devc_tier);
  const SyNonbondedKit<float, float2> nbk_f = poly_ag.getSinglePrecisionNonbondedKit(devc_tier);
  const SyValenceKit<double> vk_d = poly_ag.getDoublePrecisionValenceKit(devc_tier);
  const SyValenceKit<float> vk_f = poly_ag.getSinglePrecisionValenceKit(devc_tier);
  const SyRestraintKit<double,
                       double2, double4> rk_d = poly_ag.getDoublePrecisionRestraintKit(devc_tier);
  const SyRestraintKit<float,
                       float2, float4> rk_f = poly_ag.getSinglePrecisionRestraintKit(devc_tier);
  const SyAtomUpdateKit<double,
                        double2,
                        double4> auk_d = poly_ag.getDoublePrecisionAtomUpdateKit(devc_tier);
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

    // Remove any motion of the center of mass, if requested.
    if (nscm > 0 && step_idx > 0 && step_idx % nscm == 0) {
      removeMomentum(poly_ps, poly_ag, mos, gpu);
      mos->updateCyclePosition();
    }
    
    
    // Initialize energy accumulators if needed
    if (on_energy_step) {
      sc->initialize(devc_tier, gpu);
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

    // Perform the valence calculation and atom update
    switch (valence_prec) {
    case PrecisionModel::DOUBLE:
      if (on_energy_step) {
        launchValence(vk_d, rk_d, &ctrl_fe_d, crd_ptr, auk_d, &tstw_d, &scw, &vale_fe_res_d,
                      EvaluateForce::YES, EvaluateEnergy::YES, VwuGoal::MOVE_PARTICLES, vale_bt_fe,
                      0.0, 0.0);
      }
      else {
        launchValence(vk_d, rk_d, &ctrl_fx_d, crd_ptr, auk_d, &tstw_d, &scw, &vale_fx_res_d,
                      EvaluateForce::YES, EvaluateEnergy::NO, VwuGoal::MOVE_PARTICLES, vale_bt_fx,
                      0.0, 0.0);
      }
      break;
    case PrecisionModel::SINGLE:
      if (on_energy_step) {
        launchValence(vk_f, rk_f, &ctrl_fe_f, crd_ptr, auk_f, &tstw_f, &scw, &vale_fe_res_f,
                      EvaluateForce::YES, EvaluateEnergy::YES, VwuGoal::MOVE_PARTICLES,
                      AccumulationMethod::SPLIT, vale_bt_fe, 0.0, 0.0);
      }
      else {
        launchValence(vk_f, rk_f, &ctrl_fx_f, crd_ptr, auk_f, &tstw_f, &scw, &vale_fx_res_f,
                      EvaluateForce::YES, EvaluateEnergy::NO, VwuGoal::MOVE_PARTICLES,
                      AccumulationMethod::SPLIT, vale_bt_fx, 0.0, 0.0);
      }
      break;
    }

    // Log the trajectory frame if requested.  This will skip the initial coordinates but catch
    // the final frame, if the step count is a multiple of the trajectory printing frequency.
    if (traj_freq > 0 && (step_idx + 1) % traj_freq == 0) {
      poly_ps->download();
      const double current_time = static_cast<double>(step_idx) * tstw_d.dt;
      std::vector<int> system_indices(1);
      for (int i = 0; i < poly_ps->getSystemCount(); i++) {
        system_indices[0] = i;
        const int sysc_idx = syscmap.getSystemCacheIndex(i);
        const std::string& traj_name = sysc.getTrajectoryName(sysc_idx);
        poly_ps->printTrajectory(system_indices, traj_name, current_time,
                                 CoordinateFileKind::AMBER_CRD, PrintSituation::APPEND);
      }
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
}

//-------------------------------------------------------------------------------------------------
ScoreCard launchDynamics(const AtomGraphSynthesis &poly_ag,
                         const StaticExclusionMaskSynthesis &poly_se, Thermostat *tst,
                         PhaseSpaceSynthesis *poly_ps, const DynamicsControls &dyncon,
                         const SystemCache &sysc, const SynthesisCacheMap &syscmap,
                         const GpuDetails &gpu, const PrecisionModel valence_prec,
                         const PrecisionModel nonbond_prec, const int energy_bits,
                         StopWatch *timer, const std::string &task_name) {

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
  launchDynamics(valence_prec, nonbond_prec, poly_ag, poly_se, tst, poly_ps, &mos, dyncon,
                 &mmctrl_fe, &mmctrl_fx, &result, &vale_fe_cache, &vale_fx_cache, &nonb_fe_cache,
                 &nonb_fx_cache, &ism_space, AccumulationMethod::SPLIT, sysc, syscmap, gpu,
                 launcher, timer, task_name);
  return result;
}

} // namespace mm
} // namespace stormm
