// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace mm {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tnb_calc, typename Tcoord4>
void launchDynamics(PhaseSpaceSynthesis *poly_ps, Condensate *staging_zone,
                    CellGrid<Tcoord, Tacc, Tnb_calc, Tcoord4> *cg, PMIGrid *pmig,
                    ConvolutionManager *cvol, MotionSweeper *mos, Thermostat *tst, ScoreCard *sc,
                    MolecularMechanicsControls *mmctrl_fe, MolecularMechanicsControls *mmctrl_fx,
                    CacheResource *vale_fe_cache, CacheResource *vale_fx_cache,
                    TileManager *tlmn_fe, TileManager *tlmn_fx, const AtomGraphSynthesis &poly_ag,
                    const LocalExclusionMask &lem, const PPITable &nrg_tab,
                    const DynamicsControls &dyncon, const PPPMControls &pmecon,
                    const ReportControls &repcon, const SystemCache &sysc,
                    const SynthesisCacheMap &syscmap, const GpuDetails &gpu,
                    const CoreKlManager &launcher, const PrecisionModel valence_prec,
                    StopWatch *timer, ProgressBar *progress_bar, const std::string &task_name) {
  const int ntpr   = dyncon.getDiagnosticPrintFrequency();
  const int ntwx   = dyncon.getTrajectoryPrintFrequency();
  const int nframe = (roundUp(dyncon.getStepCount(), ntpr) / ntpr) + 1;

  
  // Construct a vector of the pair corrrespondence between the staging zone and the original
  // coordinate synthesis.
  const int nsys = poly_ps->getSystemCount();
  Hybrid<int2> system_pairs(nsys, "synth_pair_list");
  int2* sysp_ptr = system_pairs.data();
  for (int i = 0; i < nsys; i++) {
    sysp_ptr[i].x = i;
    sysp_ptr[i].y = i;
  }
  
  // Extract critical run and diagnostic parameters from the user input.
  const int nstep = dyncon.getStepCount();
  const int traj_freq = dyncon.getTrajectoryPrintFrequency();
  const int diag_freq = dyncon.getDiagnosticPrintFrequency();
  const int cmpg_freq = dyncon.getCenterOfMassMotionPurgeFrequency();
  const BrokenAsciiCode ascii_recovery = repcon.getAsciiSalvageStyle();

  // Recall the non-bonded precision calculation model.  Various user choices in the &precision
  // input namelist have been codified in the coordinate synthesis representation or the neighbor
  // list object.  The precision of valence calculations is conveyed by an explicit input
  // parameter.
  const size_t ct_tmat = std::type_index(typeid(Tnb_calc)).hash_code();
  const bool use_overflow = pmig->useOverflowAccumulation();
  const bool has_tiny_box = (cg->getTinyBoxPresence() == TinyBoxPresence::YES);
  const PrecisionModel nonbond_prec = (ct_tmat == double_type_index) ? PrecisionModel::DOUBLE :
                                                                       PrecisionModel::SINGLE;
  const PrecisionModel nblist_coord_prec = (ct_tmat == double_type_index) ?
                                           PrecisionModel::DOUBLE : PrecisionModel::SINGLE;
  
  // Produce abstracts of the relevant objects, at each point in the coordinate time cycle.
  const HybridTargetLevel devc = HybridTargetLevel::DEVICE;
  ScoreCardWriter scw = sc->data(devc);
  const CoordinateCycle poly_ps_next_stage = getNextCyclePosition(poly_ps->getCyclePosition());
  PsSynthesisWriter poly_psw = poly_ps->data(devc);
  PsSynthesisWriter poly_psw_alt = poly_ps->data(poly_ps_next_stage, devc);
  PsSynthesisReader poly_psr(poly_psw);
  PsSynthesisReader poly_psr_alt(poly_psw_alt);
  PsSynthesisBorders pssb = poly_ps->borders(devc);
  PsSynthesisBorders pssb_alt = poly_ps->borders(poly_ps_next_stage, devc);
  const LocalExclusionMaskReader lemr = lem.data(devc);
  const PPIKit<double, double4_16a> dppi_direct = nrg_tab.dpData(devc);
  const PPIKit<float, float4> fppi_direct = nrg_tab.spData(devc);
  TilePlan tlpn_fe = tlmn_fe->data(devc);
  TilePlan tlpn_fx = tlmn_fx->data(devc);
  MMControlKit<double> d_ctrl_fe = mmctrl_fe->dpData(devc);
  MMControlKit<double> d_ctrl_fx = mmctrl_fx->dpData(devc);
  MMControlKit<float> f_ctrl_fe = mmctrl_fe->spData(devc);
  MMControlKit<float> f_ctrl_fx = mmctrl_fx->spData(devc);
  const CoordinateCycle cg_next_stage = getNextCyclePosition(cg->getCyclePosition());
  CellGridWriter<void, void, void, void> cgv = cg->templateFreeData(devc);
  CellGridWriter<void, void, void, void> cgv_alt = cg->templateFreeData(cg_next_stage, devc);
  CellGridWriter<double, llint,
                 double, double4_16a> dd_cgw = restoreType<double, llint,
                                                           double, double4_16a>(cgv);
  CellGridWriter<double, llint,
                 double, double4_16a> dd_cgw_alt = restoreType<double, llint,
                                                               double, double4_16a>(cgv_alt);
  CellGridWriter<float, int, float, float4> ff_cgw = restoreType<float, int, float, float4>(cgv);
  CellGridWriter<float, int, float, float4> ff_cgw_alt = restoreType<float, int,
                                                                     float, float4>(cgv_alt);
  CellGridReader<double, llint, double, double4_16a> dd_cgr(dd_cgw);
  CellGridReader<double, llint, double, double4_16a> dd_cgr_alt(dd_cgw_alt);
  CellGridReader<float, int, float, float4> ff_cgr(ff_cgw);
  CellGridReader<float, int, float, float4> ff_cgr_alt(ff_cgw_alt);
  CellGridWriter<void, void, void, void> v_cgw = cg->templateFreeData(devc);
  CellGridWriter<void, void, void, void> v_cgw_alt = cg->templateFreeData(cg_next_stage, devc);
  CellGridReader<void, void, void, void> v_cgr(v_cgw);
  CellGridReader<void, void, void, void> v_cgr_alt(v_cgw_alt);
  CellOriginsReader corg = cg->getRulers(devc);
  CellOriginsReader corg_alt = cg->getRulers(cg_next_stage, devc);
  ConvolutionWriter<double, double2> d_cvolw = cvol->dpData();
  ConvolutionWriter<float, float2> f_cvolw = cvol->spData();
  PMIGridAccumulator pm_acc = pmig->fpData(devc, ExceptionResponse::SILENT);
  PMIGridWriter pm_wrt = pmig->data(devc);
  const PMIGridReader pm_rdr(pm_wrt);
  const SyNonbondedKit<double,
                       double2> dpoly_nbk = poly_ag.getDoublePrecisionNonbondedKit(devc);
  const SyNonbondedKit<float,
                       float2> fpoly_nbk = poly_ag.getSinglePrecisionNonbondedKit(devc);
  ThermostatWriter<double> d_tstw = tst->dpData(devc);
  ThermostatWriter<float> f_tstw = tst->spData(devc);
  const SyValenceKit<double> dpoly_vk = poly_ag.getDoublePrecisionValenceKit(devc);
  const SyValenceKit<float> fpoly_vk = poly_ag.getSinglePrecisionValenceKit(devc);
  const SyAtomUpdateKit<double,
                        double2,
                        double4_16a> dpoly_auk = poly_ag.getDoublePrecisionAtomUpdateKit(devc);
  const SyAtomUpdateKit<float,
                        float2,
                        float4> fpoly_auk = poly_ag.getSinglePrecisionAtomUpdateKit(devc);
  const SyRestraintKit<double,
                       double2,
                       double4_16a> dpoly_rk = poly_ag.getDoublePrecisionRestraintKit(devc);
  const SyRestraintKit<float,
                       float2,
                       float4> fpoly_rk = poly_ag.getSinglePrecisionRestraintKit(devc);
  CacheResourceKit<double> d_vale_fe_res = vale_fe_cache->dpData(devc);
  CacheResourceKit<double> d_vale_fx_res = vale_fx_cache->dpData(devc);
  CacheResourceKit<float> f_vale_fe_res = vale_fe_cache->spData(devc);
  CacheResourceKit<float> f_vale_fx_res = vale_fx_cache->spData(devc);
  
  // Obtain launch parameters for various kernels.
  const int2 pair_fx_lp = launcher.getPMEPairsKernelDims(nblist_coord_prec, nonbond_prec,
                                                         NeighborListKind::MONO,
                                                         cg->getTinyBoxPresence(),
                                                         EvaluateForce::YES, EvaluateEnergy::NO,
                                                         ClashResponse::NONE);
  const int2 pair_fe_lp = launcher.getPMEPairsKernelDims(nblist_coord_prec, nonbond_prec,
                                                         NeighborListKind::MONO,
                                                         cg->getTinyBoxPresence(),
                                                         EvaluateForce::YES, EvaluateEnergy::YES,
                                                         ClashResponse::NONE);
  const int2 qmap_lp = launcher.getDensityMappingKernelDims(pmecon.getDensityMappingMethod(),
                                                            nonbond_prec, nonbond_prec,
                                                            use_overflow, ct_tmat,
                                                            pmecon.getInterpolationOrder());
  const int2 fintrp_lp = launcher.getForceGatheringKernelDims(QMapMethod::GENERAL_PURPOSE,
                                                              nonbond_prec, ct_tmat,
                                                              pmecon.getInterpolationOrder());
  const int2 vale_fx_lp = launcher.getValenceKernelDims(valence_prec, EvaluateForce::YES,
                                                        EvaluateEnergy::NO,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::MOVE_PARTICLES,
                                                        ClashResponse::NONE);
  const int2 vale_fe_lp = launcher.getValenceKernelDims(valence_prec, EvaluateForce::YES,
                                                        EvaluateEnergy::YES,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::MOVE_PARTICLES,
                                                        ClashResponse::NONE);
  const int2 intg_vvi_lp =
    launcher.getIntegrationKernelDims(valence_prec, AccumulationMethod::SPLIT,
                                      IntegrationStage::VELOCITY_ADVANCE);
  const int2 intg_vc_lp =
    launcher.getIntegrationKernelDims(valence_prec, AccumulationMethod::SPLIT,
                                      IntegrationStage::VELOCITY_CONSTRAINT);
  const int2 intg_ke_lp =
    launcher.getIntegrationKernelDims(valence_prec, AccumulationMethod::SPLIT,
                                      IntegrationStage::CALC_KINETIC);
  const int2 intg_vvii_lp =
    launcher.getIntegrationKernelDims(valence_prec, AccumulationMethod::SPLIT,
                                      IntegrationStage::POSITION_ADVANCE);
  const int2 intg_gc_lp =
    launcher.getIntegrationKernelDims(valence_prec, AccumulationMethod::SPLIT,
                                      IntegrationStage::GEOMETRY_CONSTRAINT);
  const int2 migr_one_lp = launcher.getMigrationKernelDims(nblist_coord_prec,
                                                           NeighborListKind::MONO, 1,
                                                           poly_ps->getGlobalPositionBits(),
                                                           cg->getTotalChainCount());
  const int2 migr_two_lp = launcher.getMigrationKernelDims(nblist_coord_prec,
                                                           NeighborListKind::MONO, 2,
                                                           poly_ps->getGlobalPositionBits(),
                                                           cg->getTotalChainCount());
  
  // For steps when interventions or checkpointing are needed, the integration process must proceed
  // in stages and paused when the complete forces or velocities at various stages are available.
  std::vector<IntegrationStage> traj_integration_stages(1, IntegrationStage::VELOCITY_ADVANCE);
  std::vector<int2> traj_intg_stage_lp(1, intg_vvi_lp);
  size_t traj_write_stage_idx;
  switch (dyncon.constrainGeometry()) {
  case ApplyConstraints::YES:
    traj_integration_stages.push_back(IntegrationStage::VELOCITY_CONSTRAINT);
    traj_intg_stage_lp.push_back(intg_vc_lp);
    traj_write_stage_idx = 1;
    break;
  case ApplyConstraints::NO:
    traj_write_stage_idx = 0;
    break;
  }
  traj_integration_stages.push_back(IntegrationStage::POSITION_ADVANCE);
  traj_intg_stage_lp.push_back(intg_vvii_lp);
  switch (dyncon.constrainGeometry()) {
  case ApplyConstraints::YES:
    traj_integration_stages.push_back(IntegrationStage::GEOMETRY_CONSTRAINT);
    traj_intg_stage_lp.push_back(intg_gc_lp);
    break;
  case ApplyConstraints::NO:
    break;
  }
  const size_t total_traj_intg_stages = traj_integration_stages.size();

  // Prepare the progress bar, if there is a valid object provided.
  const bool show_bar = (progress_bar != nullptr);
  if (show_bar) {
    if (ntpr > 0) {
      progress_bar->setCycleCount(nstep / ntpr);
    }
    else {
      progress_bar->setCycleCount(1);
    }
    progress_bar->reset();
  }

  // Load objects to prepare for interventions.  Anomaly reporting is expected to have been set
  // prior to calling one of the low-level functions, to avoid having to carry debugging namelist
  // control information all the way into a production routine.
  if (dyna_tk.isActive()) {
    dyna_tk.setCoordinateSynthesis(poly_ps);
    dyna_tk.setTopologySynthesis(const_cast<AtomGraphSynthesis*>(poly_ag.getSelfPointer()));
    dyna_tk.setGpu(gpu);
    dyna_tk.setNeighborList<Tcoord, Tacc, Tnb_calc, Tcoord4>(cg);
    dyna_tk.setEnergyTracking(sc);
  }
  
  // Loop over all time steps.  Various combinations of precisions for non-bonded calculations and
  // particle position in the neighbor list are enumerated
  const int nscm = dyncon.getCenterOfMassMotionPurgeFrequency();
  for (int step_idx = 0; step_idx < nstep; step_idx++) {
    const bool on_energy_step = (ntpr > 0 && step_idx % ntpr == 0);
    const bool on_trajectory_step = ((ntwx > 0 && step_idx % ntwx == 0) || step_idx == nstep - 1);
    const bool intervene = dyna_tk.isActiveOnStep(step_idx);
    
    // Begin by removing momentum from the system.  Most center of mass motion will be nullified
    // during the normal PME workflow, by accumulating non-conservative forces coming off of the
    // particle-mesh interaction grid and removing that net momentum at every step.
    if (nscm > 0 && step_idx > 0 && step_idx % nscm == 0) {
      removeMomentum(poly_ps, poly_ag, mos, gpu);
      mos->updateCyclePosition();
    }

    // Initialize energy accumulators if needed.  Set pointers to the appropriate molecular
    // mechanics control objects for either precision model.  On energy-yielding steps, a pointer
    // to the energy tracking object's device-facing abstract will be included.  On all other
    // steps, a null pointer will be supplied to forego the energy computation by choosing a
    // separate kernel.
    MMControlKit<double> *d_mmctrl_ptr;
    MMControlKit<float> *f_mmctrl_ptr;
    CacheResourceKit<double> *d_vale_res_ptr;
    CacheResourceKit<float> *f_vale_res_ptr;
    TilePlan *tile_plan_ptr;
    EvaluateEnergy eval_nrg;
    int2 pair_lp, vale_lp;
    const VwuGoal vale_objective = (on_trajectory_step || intervene) ? VwuGoal::ACCUMULATE :
                                                                       VwuGoal::MOVE_PARTICLES;
    ScoreCardWriter *scw_ptr;
    if (on_energy_step) {
      sc->initialize(devc, gpu);
      d_mmctrl_ptr = &d_ctrl_fe;
      f_mmctrl_ptr = &f_ctrl_fe;
      d_vale_res_ptr = &d_vale_fe_res,
      f_vale_res_ptr = &f_vale_fe_res,
      eval_nrg = EvaluateEnergy::YES;
      pair_lp = pair_fe_lp;
      vale_lp = vale_fe_lp;
      tile_plan_ptr = &tlpn_fe;
      scw_ptr = &scw;
    }
    else {
      d_mmctrl_ptr = &d_ctrl_fx;
      f_mmctrl_ptr = &f_ctrl_fx;
      d_vale_res_ptr = &d_vale_fx_res,
      f_vale_res_ptr = &f_vale_fx_res,
      eval_nrg = EvaluateEnergy::NO;
      pair_lp = pair_fx_lp;
      vale_lp = vale_fx_lp;
      tile_plan_ptr = &tlpn_fx;
      scw_ptr = nullptr;
    }

    // Set pointers to the appropriate point in the time cycle for each coordinate object.
    // Pointers for any and all available precision models must be set in order to overcome
    // compilation issues with the templated neighbor list (CellGrid) object.
    PsSynthesisWriter *crdw_ptr;
    PsSynthesisReader *crdr_ptr;
    PsSynthesisBorders *borders_ptr;
    CellGridWriter<double, llint, double, double4_16a> *dd_cgw_ptr;
    CellGridWriter<float, int, float, float4> *ff_cgw_ptr;
    CellGridReader<double, llint, double, double4_16a> *dd_cgr_ptr;
    CellGridReader<float, int, float, float4> *ff_cgr_ptr;
    CellGridWriter<void, void, void, void> *v_cgw_ptr;
    CellGridReader<void, void, void, void> *v_cgr_ptr;
    CellOriginsReader *corg_ptr;
    if (step_idx & 0x1) {
      crdw_ptr = &poly_psw_alt;
      crdr_ptr = &poly_psr_alt;
      borders_ptr = &pssb_alt;
      dd_cgw_ptr = &dd_cgw_alt;
      ff_cgw_ptr = &ff_cgw_alt;
      dd_cgr_ptr = &dd_cgr_alt;
      ff_cgr_ptr = &ff_cgr_alt;
      v_cgw_ptr = &v_cgw_alt;
      v_cgr_ptr = &v_cgr_alt;
      corg_ptr = &corg_alt;
    }
    else {
      crdw_ptr = &poly_psw;
      crdr_ptr = &poly_psr;
      borders_ptr = &pssb;
      dd_cgw_ptr = &dd_cgw;
      ff_cgw_ptr = &ff_cgw;
      dd_cgr_ptr = &dd_cgr;
      ff_cgr_ptr = &ff_cgr;
      v_cgw_ptr = &v_cgw;
      v_cgr_ptr = &v_cgr;
      corg_ptr = &corg;
    }

    // Compute the non-bonded particle-particle interactions.  All data is assumed to reside on the
    // GPU as of the time this function is called, so kernel launches can proceed without further
    // preparation.
    poly_ps->initializeForces(gpu, devc);
    switch (nonbond_prec) {
    case PrecisionModel::DOUBLE:
      if (has_tiny_box) {
        launchPMEPairs(dpoly_nbk, lemr, dppi_direct, *borders_ptr, dd_cgw_ptr, tile_plan_ptr, &scw,
                       d_mmctrl_ptr, EvaluateForce::YES, eval_nrg, pair_lp, 0.0, 0.0);
      }
      else {
        launchPMEPairs(dpoly_nbk, lemr, dppi_direct, dd_cgw_ptr, tile_plan_ptr, &scw, d_mmctrl_ptr,
                       EvaluateForce::YES, eval_nrg, pair_lp, 0.0, 0.0);
      }
      switch (pmecon.getDensityMappingMethod()) {
      case QMapMethod::ACC_SHARED:
      case QMapMethod::AUTOMATIC:
        launchShrAccDensityKernel(&pm_wrt, use_overflow, d_mmctrl_ptr, *v_cgr_ptr,
                                  double_type_index, dpoly_nbk, qmap_lp);
        break;
      case QMapMethod::GENERAL_PURPOSE:
        launchPMIGridInitialization(&pm_acc, gpu);
        launchGenPrpDensityKernel(&pm_acc, *v_cgr_ptr, double_type_index, dpoly_nbk, qmap_lp);
        launchPMIGridRealConversion(&pm_wrt, pm_acc, gpu);
        break;
      }
      applyConvolution(&d_cvolw, *borders_ptr, pm_rdr, gpu, &scw);
      launchGenForceGatheringKernel(v_cgw_ptr, pm_rdr, ct_tmat, *borders_ptr, dpoly_nbk,
                                    fintrp_lp);
      break;
    case PrecisionModel::SINGLE:
      if (has_tiny_box) {
        launchPMEPairs(fpoly_nbk, lemr, fppi_direct, *borders_ptr, ff_cgw_ptr, tile_plan_ptr, &scw,
                       f_mmctrl_ptr, EvaluateForce::YES, eval_nrg, pair_lp, 0.0, 0.0);
      }
      else {
        launchPMEPairs(fpoly_nbk, lemr, fppi_direct, ff_cgw_ptr, tile_plan_ptr, &scw, f_mmctrl_ptr,
                       EvaluateForce::YES, eval_nrg, pair_lp, 0.0, 0.0);
      }
      switch (pmecon.getDensityMappingMethod()) {
      case QMapMethod::ACC_SHARED:
      case QMapMethod::AUTOMATIC:
        launchShrAccDensityKernel(&pm_wrt, use_overflow, f_mmctrl_ptr, *v_cgr_ptr,
                                  float_type_index, fpoly_nbk, qmap_lp);
        break;
      case QMapMethod::GENERAL_PURPOSE:
        launchPMIGridInitialization(&pm_acc, gpu);
        launchGenPrpDensityKernel(&pm_acc, *v_cgr_ptr, float_type_index, fpoly_nbk, qmap_lp);
        launchPMIGridRealConversion(&pm_wrt, pm_acc, gpu);
        break;
      }
      applyConvolution(&f_cvolw, *borders_ptr, pm_rdr, gpu, &scw);
      launchGenForceGatheringKernel(v_cgw_ptr, pm_rdr, ct_tmat, *borders_ptr, fpoly_nbk,
                                    fintrp_lp);
      break;
    }

    // Compute the valence interactions.  Whether this routine also subsumes the particle movement
    // will be determined by the "valence objective" (vale_objective).  On trajectory-yielding or
    // checkpoint steps, forces will be calculated and stored so that particle movement can be
    // driven by further individual kernel launches.
    switch (valence_prec) {
    case PrecisionModel::DOUBLE:
      switch (nonbond_prec) {
      case PrecisionModel::DOUBLE:
        launchValence(dpoly_vk, dpoly_rk, *dd_cgr_ptr, d_mmctrl_ptr, crdw_ptr, dpoly_auk, &d_tstw,
                      &scw, d_vale_res_ptr, EvaluateForce::YES, eval_nrg, vale_objective, vale_lp,
                      0.0, 0.0);
        break;
      case PrecisionModel::SINGLE:
        launchValence(dpoly_vk, dpoly_rk, *ff_cgr_ptr, d_mmctrl_ptr, crdw_ptr, dpoly_auk, &d_tstw,
                      &scw, d_vale_res_ptr, EvaluateForce::YES, eval_nrg, vale_objective, vale_lp,
                      0.0, 0.0);
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      switch (nonbond_prec) {
      case PrecisionModel::DOUBLE:
        launchValence(fpoly_vk, fpoly_rk, *dd_cgr_ptr, f_mmctrl_ptr, crdw_ptr, fpoly_auk, &f_tstw,
                      &scw, f_vale_res_ptr, EvaluateForce::YES, eval_nrg, vale_objective,
                      AccumulationMethod::SPLIT, vale_lp, 0.0, 0.0);
        break;
      case PrecisionModel::SINGLE:
        launchValence(fpoly_vk, fpoly_rk, *ff_cgr_ptr, f_mmctrl_ptr, crdw_ptr, fpoly_auk, &f_tstw,
                      &scw, f_vale_res_ptr, EvaluateForce::YES, eval_nrg, vale_objective,
                      AccumulationMethod::SPLIT, vale_lp, 0.0, 0.0);
        break;
      }
      break;
    }

    // Additional forces may be added now, or checks on the existing forces may be performed.
    if (intervene) {
      dyna_tk.execute(d_tstw.step, IntegrationStage::CALC_FORCES);
    }

    // Perform the particle movement and coordinate printing if needed.
    if (on_trajectory_step || intervene) {

      // On an trajectory-yielding step, forces will have been accumulated but the particles will
      // not have been moved.  Complete the time step with the work below.
      switch (valence_prec) {
      case PrecisionModel::DOUBLE:
        switch (nonbond_prec) {
        case PrecisionModel::DOUBLE:
          for (size_t i = 0; i < total_traj_intg_stages; i++) {
            launchIntegrationProcess(crdw_ptr, d_vale_res_ptr, d_mmctrl_ptr, scw_ptr, dd_cgw_ptr,
                                     dpoly_vk, dpoly_auk, d_tstw, traj_intg_stage_lp[i],
                                     traj_integration_stages[i]);
            if (intervene) {
              dyna_tk.execute(d_tstw.step, traj_integration_stages[i]);
            }

            // Print the trajectory (which may include velocities, such as in the case of a
            // checkpoint step) after constraining velocities.
            if (i == traj_write_stage_idx) {
              writeSynthesisFrames(*poly_ps, staging_zone, system_pairs, syscmap,
                                   CoordinateFileKind::AMBER_CRD,
                                   static_cast<double>(step_idx) * d_tstw.dt,
                                   TrajectoryKind::POSITIONS, devc, gpu, ascii_recovery);
              if (on_energy_step) {

                // If there is no trajectory to write but there is energy to compute, the energy
                // will be computed as part of the fused kernel call above.  When writing
                // trajectory, the energy calculation point coincides with the trajectory
                // transcription point.
                launchIntegrationProcess(crdw_ptr, d_vale_res_ptr, d_mmctrl_ptr, scw_ptr,
                                         dd_cgw_ptr, dpoly_vk, dpoly_auk, d_tstw, intg_ke_lp,
                                         IntegrationStage::CALC_KINETIC);
                if (intervene) {
                  dyna_tk.execute(d_tstw.step, IntegrationStage::CALC_KINETIC);
                }
              }
            }
          }
          break;
        case PrecisionModel::SINGLE:
          for (size_t i = 0; i < total_traj_intg_stages; i++) {
            launchIntegrationProcess(crdw_ptr, d_vale_res_ptr, d_mmctrl_ptr, scw_ptr, ff_cgw_ptr,
                                     dpoly_vk, dpoly_auk, d_tstw, traj_intg_stage_lp[i],
                                     traj_integration_stages[i]);
            if (intervene) {
              dyna_tk.execute(d_tstw.step, traj_integration_stages[i]);
            }
            if (i == traj_write_stage_idx) {
              writeSynthesisFrames(*poly_ps, staging_zone, system_pairs, syscmap,
                                   CoordinateFileKind::AMBER_CRD,
                                   static_cast<double>(step_idx) * d_tstw.dt,
                                   TrajectoryKind::POSITIONS, devc, gpu, ascii_recovery);
              if (on_energy_step) {
                launchIntegrationProcess(crdw_ptr, d_vale_res_ptr, d_mmctrl_ptr, scw_ptr,
                                         ff_cgw_ptr, dpoly_vk, dpoly_auk, d_tstw,
                                         intg_ke_lp, IntegrationStage::CALC_KINETIC);
                if (intervene) {
                  dyna_tk.execute(d_tstw.step, IntegrationStage::CALC_KINETIC);
                }
              }
            }
          }
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (nonbond_prec) {
        case PrecisionModel::DOUBLE:
          for (size_t i = 0; i < total_traj_intg_stages; i++) {
            launchIntegrationProcess(crdw_ptr, f_vale_res_ptr, f_mmctrl_ptr, scw_ptr, dd_cgw_ptr,
                                     fpoly_vk, fpoly_auk, f_tstw, traj_intg_stage_lp[i],
                                     AccumulationMethod::SPLIT, traj_integration_stages[i]);
            if (intervene) {
              dyna_tk.execute(f_tstw.step, traj_integration_stages[i]);
            }
            if (i == traj_write_stage_idx) {
              writeSynthesisFrames(*poly_ps, staging_zone, system_pairs, syscmap,
                                   CoordinateFileKind::AMBER_CRD,
                                   static_cast<double>(step_idx) * d_tstw.dt,
                                   TrajectoryKind::POSITIONS, devc, gpu, ascii_recovery);
              if (on_energy_step) {
                launchIntegrationProcess(crdw_ptr, f_vale_res_ptr, f_mmctrl_ptr, scw_ptr,
                                         dd_cgw_ptr, fpoly_vk, fpoly_auk, f_tstw,
                                         intg_ke_lp, AccumulationMethod::SPLIT,
                                         IntegrationStage::CALC_KINETIC);
                if (intervene) {
                  dyna_tk.execute(f_tstw.step, IntegrationStage::CALC_KINETIC);
                }
              }
            }
          }
          break;
        case PrecisionModel::SINGLE:
          for (size_t i = 0; i < total_traj_intg_stages; i++) {
            launchIntegrationProcess(crdw_ptr, f_vale_res_ptr, f_mmctrl_ptr, scw_ptr, ff_cgw_ptr,
                                     fpoly_vk, fpoly_auk, f_tstw, traj_intg_stage_lp[i],
                                     AccumulationMethod::SPLIT, traj_integration_stages[i]);
            if (intervene) {
              dyna_tk.execute(f_tstw.step, traj_integration_stages[i]);
            }
            if (i == traj_write_stage_idx) {
              writeSynthesisFrames(*poly_ps, staging_zone, system_pairs, syscmap,
                                   CoordinateFileKind::AMBER_CRD,
                                   static_cast<double>(step_idx) * d_tstw.dt,
                                   TrajectoryKind::POSITIONS, devc, gpu, ascii_recovery);
              if (on_energy_step) {
                launchIntegrationProcess(crdw_ptr, f_vale_res_ptr, f_mmctrl_ptr, scw_ptr,
                                         ff_cgw_ptr, fpoly_vk, fpoly_auk, f_tstw,
                                         intg_ke_lp, AccumulationMethod::SPLIT,
                                         IntegrationStage::CALC_KINETIC);
                if (intervene) {
                  dyna_tk.execute(f_tstw.step, IntegrationStage::CALC_KINETIC);
                }
              }
            }
          }
          break;
        }
        break;
      }
    }

    // Migrate particles within the neighbor list
    switch (nonbond_prec) {
    case PrecisionModel::DOUBLE:
      launchMigration(dd_cgw_ptr, *corg_ptr, *crdr_ptr, migr_one_lp, migr_two_lp);
      break;
    case PrecisionModel::SINGLE:
      launchMigration(ff_cgw_ptr, *corg_ptr, *crdr_ptr, migr_one_lp, migr_two_lp);
      break;
    }

    // Update the neighbor list's coordinate cycle and advance step counters.  Catalog energy
    // results from this step, if requested.
    if (on_energy_step) {
      d_ctrl_fe.step += 1;
      f_ctrl_fe.step += 1;
      mmctrl_fe->incrementStep();
      sc->commit(devc, gpu);
      sc->incrementSampleCount();

      // The ScoreCard's last time step can be updated with either the double- or single-precision
      // thermostat abstract.  Both are kept up-to-date with the dynamics progress.
      sc->setLastTimeStep(d_tstw.step, devc);
      if (show_bar) {
        progress_bar->update();
      }
    }
    else {
      d_ctrl_fx.step += 1;
      f_ctrl_fx.step += 1;
      mmctrl_fx->incrementStep();
    }
    poly_ps->updateCyclePosition();
    cg->updateCyclePosition();
    tst->incrementStep();
    d_tstw.step += 1;
    f_tstw.step += 1;
  }
  sc->computePotentialEnergy(devc, gpu);
  sc->computeTotalEnergy(devc, gpu);
  if (show_bar) {
    progress_bar->finalizeTerminalOutput();
  }
}

} // namespace mm
} // namespace stormm
