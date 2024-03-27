// -*-c++-*-
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/hpc_reduction.h"
#include "Math/reduction_abstracts.h"
#include "Math/reduction_enumerators.h"
#include "Math/rounding.h"
#include "Potential/energy_enumerators.h"
#include "Potential/hpc_nonbonded_potential.h"
#include "Potential/hpc_valence_potential.h"
#include "Structure/hpc_virtual_site_handling.h"
#include "Structure/structure_enumerators.h"
#include "Synthesis/nonbonded_workunit.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Synthesis/valence_workunit.h"
#include "Trajectory/thermostat.h"
#include "Trajectory/trajectory_enumerators.h"
#include "hpc_minimization.h"

namespace stormm {
namespace mm {

using constants::verytiny;
using energy::CacheResourceKit;
using energy::ClashResponse;
using energy::EvaluateEnergy;
using energy::EvaluateForce;
using energy::StateVariable;
using stmath::ConjGradSubstrate;
using stmath::RdwuAbstractMap;
using stmath::ReductionGoal;
using stmath::ReductionStage;
using stmath::roundUp;
using stmath::rdwu_abstract_length;
using numerics::chooseAccumulationMethod;
using structure::launchVirtualSitePlacement;
using structure::launchTransmitVSiteForces;
using structure::VirtualSiteActivity;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::NbwuKind;
using synthesis::SeMaskSynthesisReader;
using synthesis::small_block_max_atoms;
using synthesis::SyAtomUpdateKit;
using synthesis::SyNonbondedKit;
using synthesis::SyRestraintKit;
using synthesis::SyValenceKit;
using synthesis::VwuGoal;
using trajectory::CoordinateCycle;
using trajectory::getNextCyclePosition;
using trajectory::Thermostat;
using trajectory::ThermostatKind;
using trajectory::ThermostatWriter;

#include "../Numerics/accumulation.cui"

// Conjugate gradient particle advancement: both single- and double-precision kernels do
// calculations in double-precision, but the double-precision form of the kernel respects the
// extended fixed-precision format.
#define TCALC double
#define FABS_FUNC fabs
#define SQRT_FUNC sqrt
#define LLCONV_FUNC __double2ll_rn
#  define TCALC_IS_DOUBLE
#    define KERNEL_NAME kdLineAdvance
#      include "line_movement.cui"
#    undef KERNEL_NAME
#  undef TCALC_IS_DOUBLE
#  define KERNEL_NAME kfLineAdvance
#    include "line_movement.cui"
#  undef KERNEL_NAME
#undef FABS_FUNC
#undef SQRT_FUNC
#undef LLCONV_FUNC
#undef TCALC

//-------------------------------------------------------------------------------------------------
extern void minimizationKernelSetup() {
  const cudaSharedMemConfig sms_eight = cudaSharedMemBankSizeEightByte;
  if (cudaFuncSetSharedMemConfig(kdLineAdvance, sms_eight) != cudaSuccess) {
    rtErr("Error setting kdLineAdvance __shared__ memory bank size to eight bytes.",
          "minimizationKernelSetup");
  }
  if (cudaFuncSetSharedMemConfig(kfLineAdvance, sms_eight) != cudaSuccess) {
    rtErr("Error setting kdLineAdvance __shared__ memory bank size to eight bytes.",
          "minimizationKernelSetup");
  }
}

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes queryMinimizationKernelRequirements(const PrecisionModel prec) {
  cudaFuncAttributes result;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    if (cudaFuncGetAttributes(&result, kdLineAdvance) != cudaSuccess) {
      rtErr("Error obtaining attributes for kernel kdLineAdvance.",
            "queryMinimizationKernelRequirements");
    }
    break;
  case PrecisionModel::SINGLE:
    if (cudaFuncGetAttributes(&result, kfLineAdvance) != cudaSuccess) {
      rtErr("Error obtaining attributes for kernel kfLineAdvance.",
            "queryMinimizationKernelRequirements");
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
extern void launchLineAdvance(const PrecisionModel prec, PsSynthesisWriter *poly_psw,
                              const ReductionKit &redk, const ScoreCardWriter &scw,
                              LinMinWriter *lmw, const int move_number, const int2 redu_lp) {
  switch (prec) {
  case PrecisionModel::DOUBLE:
    kdLineAdvance<<<redu_lp.x, redu_lp.y>>>(*poly_psw, redk, scw, *lmw, move_number);
    break;
  case PrecisionModel::SINGLE:
    kfLineAdvance<<<redu_lp.x, redu_lp.y>>>(*poly_psw, redk, scw, *lmw, move_number);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchLineAdvance(const PrecisionModel prec, PhaseSpaceSynthesis *poly_ps,
                              const AtomGraphSynthesis &poly_ag, ScoreCard *sc,
                              LineMinimization *line_record, const int move_number,
                              const CoreKlManager &launcher) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  const ReductionKit redk(poly_ag, tier);
  ScoreCardWriter scw = sc->data();
  LinMinWriter lmw = line_record->data();
  const int2 redu_lp = launcher.getReductionKernelDims(prec, ReductionGoal::CONJUGATE_GRADIENT,
                                                       ReductionStage::ALL_REDUCE);
  launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, move_number, redu_lp);
}
  
//-------------------------------------------------------------------------------------------------
extern void launchMinimization(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                               const StaticExclusionMaskSynthesis &poly_se,
                               PhaseSpaceSynthesis *poly_ps, const MinimizeControls &mincon,
                               MolecularMechanicsControls *mmctrl_fe,
                               MolecularMechanicsControls *mmctrl_xe,
                               MolecularMechanicsControls *mmctrl_cdfe,
                               MolecularMechanicsControls *mmctrl_cdxe, ScoreCard *sc,
                               CacheResource *vale_fe_cache, CacheResource *vale_xe_cache,
                               CacheResource *vale_cdfe_cache, CacheResource *vale_cdxe_cache,
                               CacheResource *nonb_cache, CacheResource *nonb_cd_cache,
                               ImplicitSolventWorkspace *ism_space, ReductionBridge *rbg,
                               LineMinimization *line_record, const AccumulationMethod acc_meth,
                               const GpuDetails &gpu, const CoreKlManager &launcher,
                               StopWatch *timer, const std::string &task_name) {
  
  // Obtain abstracts of critical objects, re-using them throughout the inner loop.  Some abstracts
  // are needed by both branches.
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(devc_tier);
  const CoordinateCycle alt_orientation = getNextCyclePosition(poly_ps->getCyclePosition());
  PsSynthesisWriter poly_psw_alt = poly_ps->data(alt_orientation, devc_tier);
  const SeMaskSynthesisReader poly_ser = poly_se.data(devc_tier);
  if (sc->getSystemCount() != poly_ag.getSystemCount()) {
    rtErr("The energy tracking object is not prepared for the number of systems contained in the "
          "topology synthesis (" + std::to_string(sc->getSystemCount()) + " vs. " +
          std::to_string(poly_ag.getSystemCount()) + ").", "launchMinimization");
  }

  // Make sure that the energy tracking object holds space for storing all energy components in
  // all snapshots that might be generated over the course of the energy minimization, plus the
  // final energies of all minimized snapshots.
  const int ntpr = mincon.getDiagnosticPrintFrequency();
  const int total_nrg_snapshots = roundUp(mincon.getTotalCycles(), ntpr) + 1;
  if (sc->getSampleCapacity() != total_nrg_snapshots) {
    sc->reserve(total_nrg_snapshots);
  }
  const int2 vale_fe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::YES,
                                                        EvaluateEnergy::YES,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE, ClashResponse::NONE);
  const int2 vale_xe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::NO,
                                                        EvaluateEnergy::YES,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE, ClashResponse::NONE);
  const int2 vale_cdfe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::YES,
                                                          EvaluateEnergy::YES,
                                                          AccumulationMethod::SPLIT,
                                                          VwuGoal::ACCUMULATE,
                                                          ClashResponse::FORGIVE);
  const int2 vale_cdxe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::NO,
                                                          EvaluateEnergy::YES,
                                                          AccumulationMethod::SPLIT,
                                                          VwuGoal::ACCUMULATE,
                                                          ClashResponse::FORGIVE);
  const int2 vste_mv_lp = launcher.getVirtualSiteKernelDims(prec, VirtualSiteActivity::PLACEMENT);
  const int2 vste_xm_lp = launcher.getVirtualSiteKernelDims(prec,
                                                            VirtualSiteActivity::TRANSMIT_FORCES);
  const NbwuKind nb_work_type = poly_ag.getNonbondedWorkType();
  const ImplicitSolventModel ism_type = poly_ag.getImplicitSolventModel();
  const int2 nonb_lp = launcher.getNonbondedKernelDims(prec, nb_work_type, EvaluateForce::YES,
                                                       EvaluateEnergy::YES,
                                                       AccumulationMethod::SPLIT, ism_type,
                                                       ClashResponse::NONE);
  const int2 nonb_cdlp = launcher.getNonbondedKernelDims(prec, nb_work_type, EvaluateForce::YES,
                                                         EvaluateEnergy::YES,
                                                         AccumulationMethod::SPLIT, ism_type,
                                                         ClashResponse::FORGIVE);
  const int2 gbr_lp = launcher.getBornRadiiKernelDims(prec, nb_work_type,
                                                      AccumulationMethod::SPLIT, ism_type);
  const int2 gbd_lp = launcher.getBornDerivativeKernelDims(prec, nb_work_type,
                                                           AccumulationMethod::SPLIT, ism_type);
  const int2 redu_lp = launcher.getReductionKernelDims(prec, ReductionGoal::CONJUGATE_GRADIENT,
                                                       ReductionStage::ALL_REDUCE);
  line_record->primeMoveLengths(mincon.getInitialStep());
  LinMinWriter lmw = line_record->data(devc_tier);
  const ReductionKit redk(poly_ag, devc_tier);
  ConjGradSubstrate cgsbs(poly_psw, rbg, devc_tier);

  // Conjugate gradient minimization does not use a thermostat in any respect, but such an object
  // must be created in order to get its abstract as a placeholder.
  Thermostat heat_bath(poly_ag, ThermostatKind::NONE);

  // Test whether there are virtual sites in the synthesis
  const bool virtual_sites_present = (poly_ag.getVirtualSiteCount() > 0);
  
  // Progress through minimization cycles will not be measured with the step counter in the
  // molecular mechanics control object--that will increment at four times the rate in the case
  // of conjugate gradient line minimizations.
  const int total_steps = mincon.getTotalCycles();
  const int clash_steps = mincon.getClashDampingCycles();
  const double clash_floor = mincon.getAbsoluteClashDistance();
  const double clash_ratio = mincon.getVdwClashRatio();
  mmctrl_fe->primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                                   VwuGoal::ACCUMULATE, prec, prec, poly_ag);
  mmctrl_xe->primeWorkUnitCounters(launcher, EvaluateForce::NO, EvaluateEnergy::YES,
                                   VwuGoal::ACCUMULATE, prec, prec, poly_ag);
  mmctrl_cdfe->primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                                     ClashResponse::FORGIVE, VwuGoal::ACCUMULATE, prec, prec,
                                     poly_ag);
  mmctrl_cdxe->primeWorkUnitCounters(launcher, EvaluateForce::NO, EvaluateEnergy::YES,
                                     ClashResponse::FORGIVE, VwuGoal::ACCUMULATE, prec, prec,
                                     poly_ag);
  poly_ps->primeConjugateGradientCalculation(gpu, devc_tier);
  int min_timings;
  if (timer != nullptr) {
    min_timings = timer->addCategory(task_name);
    timer->assignTime(0);
  }
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(devc_tier);
      const SyAtomUpdateKit<double, double2, double4> poly_auk =
        poly_ag.getDoublePrecisionAtomUpdateKit(devc_tier);
      const SyNonbondedKit<double, double2> poly_nbk =
        poly_ag.getDoublePrecisionNonbondedKit(devc_tier);
      const SyRestraintKit<double, double2, double4> poly_rk =
        poly_ag.getDoublePrecisionRestraintKit(devc_tier);
      CacheResourceKit<double> vale_fe_tbr = vale_fe_cache->dpData(devc_tier);
      CacheResourceKit<double> vale_xe_tbr = vale_xe_cache->dpData(devc_tier);
      CacheResourceKit<double> nonb_tbr = nonb_cache->dpData(devc_tier);
      CacheResourceKit<double> vale_cdfe_tbr = vale_cdfe_cache->dpData(devc_tier);
      CacheResourceKit<double> vale_cdxe_tbr = vale_cdxe_cache->dpData(devc_tier);
      CacheResourceKit<double> nonb_cd_tbr = nonb_cd_cache->dpData(devc_tier);
      ISWorkspaceKit<double> iswk = ism_space->dpData(devc_tier);
      MMControlKit<double> ctrl_fe = mmctrl_fe->dpData(devc_tier);
      MMControlKit<double> ctrl_xe = mmctrl_xe->dpData(devc_tier);
      MMControlKit<double> ctrl_cdfe = mmctrl_cdfe->dpData(devc_tier);
      MMControlKit<double> ctrl_cdxe = mmctrl_cdxe->dpData(devc_tier);
      ThermostatWriter<double> tstw = heat_bath.dpData(devc_tier);
      if (virtual_sites_present) {
        launchVirtualSitePlacement(&poly_psw_alt, &vale_xe_tbr, poly_vk, poly_auk, vste_mv_lp);
      }
      for (int i = 0; i < total_steps; i++) {
        
        // First stage of the cycle: compute forces and obtain the conjugate gradient move.
        poly_ps->initializeForces(gpu, devc_tier);
        ism_space->initialize(devc_tier, CoordinateCycle::WHITE, gpu);
        sc->initialize(devc_tier, gpu);
        ScoreCardWriter scw = sc->data(devc_tier);
        double clash_progress, tmp_clash_floor, tmp_clash_ratio;
        if (i < clash_steps) {
          clash_progress = static_cast<double>(clash_steps - i) / static_cast<double>(clash_steps);
          tmp_clash_floor = clash_floor * clash_progress;
          tmp_clash_ratio = clash_ratio * clash_progress;
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_cdfe, &poly_psw, &tstw, &scw,
                          &nonb_cd_tbr, &iswk, EvaluateForce::YES, EvaluateEnergy::YES, nonb_cdlp,
                          gbr_lp, gbd_lp, tmp_clash_floor, tmp_clash_ratio);
          launchValence(poly_vk, poly_rk, &ctrl_cdfe, &poly_psw, poly_auk, &tstw, &scw,
                        &vale_cdfe_tbr, EvaluateForce::YES, EvaluateEnergy::YES,
                        VwuGoal::ACCUMULATE, vale_cdfe_lp, tmp_clash_floor, tmp_clash_ratio);
          ctrl_cdfe.step += 1;
        }
        else {
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_fe, &poly_psw, &tstw, &scw,
                          &nonb_tbr, &iswk, EvaluateForce::YES, EvaluateEnergy::YES, nonb_lp,
                          gbr_lp, gbd_lp);
          launchValence(poly_vk, poly_rk, &ctrl_fe, &poly_psw, poly_auk, &tstw, &scw, &vale_fe_tbr,
                        EvaluateForce::YES, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, vale_fe_lp);
        }
        if (virtual_sites_present) {
          launchTransmitVSiteForces(&poly_psw, &vale_fe_tbr, poly_vk, poly_auk, vste_xm_lp);
        }
        if (i % ntpr == 0) {

          // It is not necessary to re-initialize the ScoreCard abstract immediately, as the
          // number of sampled steps has no bearing on how it will accept instantaneous results
          // from the following energy evaluations.
          sc->commit(devc_tier, gpu);
          sc->incrementSampleCount();
          sc->setLastTimeStep(i);
        }
        ctrl_fe.step += 1;
        launchConjugateGradient(redk, &cgsbs, &ctrl_fe, redu_lp);

        // Second stage of the cycle: advance once along the line and recompute the energy.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 0, redu_lp);
        if (virtual_sites_present) {
          launchVirtualSitePlacement(&poly_psw_alt, &vale_xe_tbr, poly_vk, poly_auk, vste_mv_lp);
        }
        sc->initialize(devc_tier, gpu);
        ism_space->initialize(devc_tier, CoordinateCycle::WHITE, gpu);
        if (i < clash_steps) {
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_cdxe, &poly_psw, &tstw, &scw,
                          &nonb_cd_tbr, &iswk, EvaluateForce::NO, EvaluateEnergy::YES, nonb_cdlp,
                          gbr_lp, gbd_lp, tmp_clash_floor, tmp_clash_ratio);
          launchValence(poly_vk, poly_rk, &ctrl_cdxe, &poly_psw, poly_auk, &tstw, &scw,
                        &vale_cdxe_tbr, EvaluateForce::NO, EvaluateEnergy::YES,
                        VwuGoal::ACCUMULATE, vale_cdxe_lp, tmp_clash_floor, tmp_clash_ratio);
          ctrl_cdxe.step += 1;
        }
        else {
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &tstw, &scw,
                          &nonb_tbr, &iswk, EvaluateForce::NO, EvaluateEnergy::YES, nonb_lp,
                          gbr_lp, gbd_lp);
          launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, poly_auk, &tstw, &scw, &vale_xe_tbr,
                        EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, vale_xe_lp);
        }
        if (virtual_sites_present) {
          launchTransmitVSiteForces(&poly_psw, &vale_xe_tbr, poly_vk, poly_auk, vste_xm_lp);
        }
        ctrl_xe.step += 1;

        // Third stage of the cycle: advance once more along the line and recompute the energy.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 1, redu_lp);
        if (virtual_sites_present) {
          launchVirtualSitePlacement(&poly_psw_alt, &vale_xe_tbr, poly_vk, poly_auk, vste_mv_lp);
        }
        sc->initialize(devc_tier, gpu);
        ism_space->initialize(devc_tier, CoordinateCycle::WHITE, gpu);
        if (i < clash_steps) {
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_cdxe, &poly_psw, &tstw, &scw,
                          &nonb_cd_tbr, &iswk, EvaluateForce::NO, EvaluateEnergy::YES, nonb_cdlp,
                          gbr_lp, gbd_lp, tmp_clash_floor, tmp_clash_ratio);
          launchValence(poly_vk, poly_rk, &ctrl_cdxe, &poly_psw, poly_auk, &tstw, &scw,
                        &vale_cdxe_tbr, EvaluateForce::NO, EvaluateEnergy::YES,
                        VwuGoal::ACCUMULATE, vale_cdxe_lp, tmp_clash_floor, tmp_clash_ratio);
          ctrl_cdxe.step += 1;
        }
        else {
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &tstw, &scw,
                          &nonb_tbr, &iswk, EvaluateForce::NO, EvaluateEnergy::YES, nonb_lp,
                          gbr_lp, gbd_lp);
          launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, poly_auk, &tstw, &scw, &vale_xe_tbr,
                        EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, vale_xe_lp);
        }
        if (virtual_sites_present) {
          launchTransmitVSiteForces(&poly_psw, &vale_xe_tbr, poly_vk, poly_auk, vste_xm_lp);
        }
        ctrl_xe.step += 1;

        // Final stage of the cycle: advance a final time along the line and recompute the energy.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 2, redu_lp);
        if (virtual_sites_present) {
          launchVirtualSitePlacement(&poly_psw_alt, &vale_xe_tbr, poly_vk, poly_auk, vste_mv_lp);
        }
        sc->initialize(devc_tier, gpu);
        ism_space->initialize(devc_tier, CoordinateCycle::WHITE, gpu);
        if (i < clash_steps) {
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_cdxe, &poly_psw, &tstw, &scw,
                          &nonb_cd_tbr, &iswk, EvaluateForce::NO, EvaluateEnergy::YES, nonb_cdlp,
                          gbr_lp, gbd_lp, tmp_clash_floor, tmp_clash_ratio);
          launchValence(poly_vk, poly_rk, &ctrl_cdxe, &poly_psw, poly_auk, &tstw, &scw,
                        &vale_cdxe_tbr, EvaluateForce::NO, EvaluateEnergy::YES,
                        VwuGoal::ACCUMULATE, vale_cdxe_lp, tmp_clash_floor, tmp_clash_ratio);
          ctrl_cdxe.step += 1;
        }
        else {
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &tstw, &scw,
                          &nonb_tbr, &iswk, EvaluateForce::NO, EvaluateEnergy::YES, nonb_lp,
                          gbr_lp, gbd_lp);
          launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, poly_auk, &tstw, &scw, &vale_xe_tbr,
                        EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, vale_xe_lp);
        }
        if (virtual_sites_present) {
          launchTransmitVSiteForces(&poly_psw, &vale_xe_tbr, poly_vk, poly_auk, vste_xm_lp);
        }
        ctrl_xe.step += 1;

        // Fit a cubic polynomial to guess the best overall advancement.  Place the system there.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 3, redu_lp);
        if (virtual_sites_present) {
          launchVirtualSitePlacement(&poly_psw_alt, &vale_xe_tbr, poly_vk, poly_auk, vste_mv_lp);
        }
      }

      // One additional energy calculation to get the final energy
      sc->initialize(devc_tier, gpu);
      ScoreCardWriter scw_final = sc->data(devc_tier);
      ism_space->initialize(devc_tier, CoordinateCycle::WHITE, gpu);
      launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &tstw,
                      &scw_final, &nonb_tbr, &iswk, EvaluateForce::NO, EvaluateEnergy::YES,
                      nonb_lp, gbr_lp, gbd_lp);
      launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, poly_auk, &tstw, &scw_final,
                    &vale_xe_tbr, EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                    vale_xe_lp);
      ctrl_xe.step += 1;
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(devc_tier);
      const SyAtomUpdateKit<float, float2, float4> poly_auk =
        poly_ag.getSinglePrecisionAtomUpdateKit(devc_tier);
      const SyNonbondedKit<float, float2> poly_nbk =
        poly_ag.getSinglePrecisionNonbondedKit(devc_tier);
      const SyRestraintKit<float, float2, float4> poly_rk =
        poly_ag.getSinglePrecisionRestraintKit(devc_tier);
      CacheResourceKit<float> vale_fe_tbr = vale_fe_cache->spData(devc_tier);
      CacheResourceKit<float> vale_xe_tbr = vale_xe_cache->spData(devc_tier);
      CacheResourceKit<float> vale_cdfe_tbr = vale_cdfe_cache->spData(devc_tier);
      CacheResourceKit<float> vale_cdxe_tbr = vale_cdxe_cache->spData(devc_tier);
      CacheResourceKit<float> nonb_tbr = nonb_cache->spData(devc_tier);
      CacheResourceKit<float> nonb_cd_tbr = nonb_cd_cache->spData(devc_tier);
      ISWorkspaceKit<float> iswk = ism_space->spData(devc_tier);
      MMControlKit<float> ctrl_fe = mmctrl_fe->spData(devc_tier);
      MMControlKit<float> ctrl_xe = mmctrl_xe->spData(devc_tier);
      MMControlKit<float> ctrl_cdfe = mmctrl_fe->spData(devc_tier);
      MMControlKit<float> ctrl_cdxe = mmctrl_xe->spData(devc_tier);
      ThermostatWriter<float> tstw = heat_bath.spData(devc_tier);
      if (virtual_sites_present) {
        launchVirtualSitePlacement(&poly_psw_alt, &vale_xe_tbr, poly_vk, poly_auk, vste_mv_lp);
      }
      for (int i = 0; i < total_steps; i++) {

        // First stage of the cycle: compute forces and obtain the conjugate gradient move.
        poly_ps->initializeForces(gpu, devc_tier);
        ism_space->initialize(devc_tier, CoordinateCycle::WHITE, gpu);
        sc->initialize(devc_tier, gpu);
        ScoreCardWriter scw = sc->data(devc_tier);
        double clash_progress, tmp_clash_floor, tmp_clash_ratio;
        if (i < clash_steps) {
          clash_progress = static_cast<double>(clash_steps - i) / static_cast<double>(clash_steps);
          tmp_clash_floor = clash_floor * clash_progress;
          tmp_clash_ratio = clash_ratio * clash_progress;
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_cdfe, &poly_psw, &tstw, &scw,
                          &nonb_cd_tbr, &iswk, EvaluateForce::YES, EvaluateEnergy::YES, acc_meth,
                          nonb_cdlp, gbr_lp, gbd_lp, tmp_clash_floor, tmp_clash_ratio);
          launchValence(poly_vk, poly_rk, &ctrl_cdfe, &poly_psw, poly_auk, &tstw, &scw,
                        &vale_cdfe_tbr, EvaluateForce::YES, EvaluateEnergy::YES,
                        VwuGoal::ACCUMULATE, acc_meth, vale_cdfe_lp, tmp_clash_floor,
                        tmp_clash_ratio);
          ctrl_cdfe.step += 1;
        }
        else {
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_fe, &poly_psw, &tstw, &scw,
                          &nonb_tbr, &iswk, EvaluateForce::YES, EvaluateEnergy::YES, acc_meth,
                          nonb_lp, gbr_lp, gbd_lp);
          launchValence(poly_vk, poly_rk, &ctrl_fe, &poly_psw, poly_auk, &tstw, &scw, &vale_fe_tbr,
                        EvaluateForce::YES, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, acc_meth,
                        vale_fe_lp);
        }
        if (virtual_sites_present) {
          launchTransmitVSiteForces(&poly_psw, &vale_xe_tbr, poly_vk, poly_auk, vste_xm_lp);
        }
        if (i % ntpr == 0) {

          // It is not necessary to re-initialize the ScoreCard abstract immediately, as the
          // number of sampled steps has no bearing on how it will accept instantaneous results
          // from the following energy evaluations.
          sc->commit(devc_tier, gpu);
          sc->incrementSampleCount();
          sc->setLastTimeStep(i);
        }
        ctrl_fe.step += 1;
        launchConjugateGradient(redk, &cgsbs, &ctrl_fe, redu_lp);

        // Second stage of the cycle: advance once along the line and recompute the energy.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 0, redu_lp);
        if (virtual_sites_present) {
          launchVirtualSitePlacement(&poly_psw_alt, &vale_xe_tbr, poly_vk, poly_auk, vste_mv_lp);
        }
        sc->initialize(devc_tier, gpu);
        ism_space->initialize(devc_tier, CoordinateCycle::WHITE, gpu);
        if (i < clash_steps) {
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_cdxe, &poly_psw, &tstw, &scw,
                          &nonb_cd_tbr, &iswk, EvaluateForce::NO, EvaluateEnergy::YES, acc_meth,
                          nonb_cdlp, gbr_lp, gbd_lp, tmp_clash_floor, tmp_clash_ratio);
          launchValence(poly_vk, poly_rk, &ctrl_cdxe, &poly_psw, poly_auk, &tstw, &scw,
                        &vale_cdxe_tbr, EvaluateForce::NO, EvaluateEnergy::YES,
                        VwuGoal::ACCUMULATE, acc_meth, vale_cdxe_lp, tmp_clash_floor,
                        tmp_clash_ratio);
          ctrl_cdxe.step += 1;
        }
        else {
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &tstw, &scw,
                          &nonb_tbr, &iswk, EvaluateForce::NO, EvaluateEnergy::YES, acc_meth,
                          nonb_lp, gbr_lp, gbd_lp);
          launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, poly_auk, &tstw, &scw, &vale_xe_tbr,
                        EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, acc_meth,
                        vale_xe_lp);
        }
        if (virtual_sites_present) {
          launchTransmitVSiteForces(&poly_psw, &vale_xe_tbr, poly_vk, poly_auk, vste_xm_lp);
        }
        ctrl_xe.step += 1;

        // Third stage of the cycle: advance once more along the line and recompute the energy.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 1, redu_lp);
        if (virtual_sites_present) {
          launchVirtualSitePlacement(&poly_psw_alt, &vale_xe_tbr, poly_vk, poly_auk, vste_mv_lp);
        }
        sc->initialize(devc_tier, gpu);
        ism_space->initialize(devc_tier, CoordinateCycle::WHITE, gpu);
        if (i < clash_steps) {
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_cdxe, &poly_psw, &tstw, &scw,
                          &nonb_cd_tbr, &iswk, EvaluateForce::NO, EvaluateEnergy::YES, acc_meth,
                          nonb_cdlp, gbr_lp, gbd_lp, tmp_clash_floor, tmp_clash_ratio);
          launchValence(poly_vk, poly_rk, &ctrl_cdxe, &poly_psw, poly_auk, &tstw, &scw,
                        &vale_cdxe_tbr, EvaluateForce::NO, EvaluateEnergy::YES,
                        VwuGoal::ACCUMULATE, acc_meth, vale_cdxe_lp, tmp_clash_floor,
                        tmp_clash_ratio);
          ctrl_cdxe.step += 1;
        }
        else {
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &tstw, &scw,
                          &nonb_tbr, &iswk, EvaluateForce::NO, EvaluateEnergy::YES, acc_meth,
                          nonb_lp, gbr_lp, gbd_lp);
          launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, poly_auk, &tstw, &scw, &vale_xe_tbr,
                        EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, acc_meth,
                        vale_xe_lp);
        }
        if (virtual_sites_present) {
          launchTransmitVSiteForces(&poly_psw, &vale_xe_tbr, poly_vk, poly_auk, vste_xm_lp);
        }
        ctrl_xe.step += 1;

        // Final stage of the cycle: advance a final time along the line and recompute the energy.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 2, redu_lp);
        if (virtual_sites_present) {
          launchVirtualSitePlacement(&poly_psw_alt, &vale_xe_tbr, poly_vk, poly_auk, vste_mv_lp);
        }
        sc->initialize(devc_tier, gpu);
        ism_space->initialize(devc_tier, CoordinateCycle::WHITE, gpu);
        if (i < clash_steps) {
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_cdxe, &poly_psw, &tstw, &scw,
                          &nonb_cd_tbr, &iswk, EvaluateForce::NO, EvaluateEnergy::YES, acc_meth,
                          nonb_cdlp, gbr_lp, gbd_lp, tmp_clash_floor, tmp_clash_ratio);
          launchValence(poly_vk, poly_rk, &ctrl_cdxe, &poly_psw, poly_auk, &tstw, &scw,
                        &vale_cdxe_tbr, EvaluateForce::NO, EvaluateEnergy::YES,
                        VwuGoal::ACCUMULATE, acc_meth, vale_cdxe_lp, tmp_clash_floor,
                        tmp_clash_ratio);
          ctrl_cdxe.step += 1;
        }
        else {
          launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &tstw, &scw,
                          &nonb_tbr, &iswk, EvaluateForce::NO, EvaluateEnergy::YES, acc_meth,
                          nonb_lp, gbr_lp, gbd_lp);
          launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, poly_auk, &tstw, &scw, &vale_xe_tbr,
                        EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE, acc_meth,
                        vale_xe_lp);
        }
        if (virtual_sites_present) {
          launchTransmitVSiteForces(&poly_psw, &vale_xe_tbr, poly_vk, poly_auk, vste_xm_lp);
        }
        ctrl_xe.step += 1;

        // Fit a cubic polynomial to guess the best overall advancement.  Place the system there.
        launchLineAdvance(prec, &poly_psw, redk, scw, &lmw, 3, redu_lp);
        if (virtual_sites_present) {
          launchVirtualSitePlacement(&poly_psw_alt, &vale_xe_tbr, poly_vk, poly_auk, vste_mv_lp);
        }
      }

      // One additional energy calculation to get the final energy
      sc->initialize(devc_tier, gpu);
      ScoreCardWriter scw_final = sc->data(devc_tier);
      ism_space->initialize(devc_tier, CoordinateCycle::WHITE, gpu);
      launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl_xe, &poly_psw, &tstw,
                      &scw_final, &nonb_tbr, &iswk, EvaluateForce::NO, EvaluateEnergy::YES,
                      acc_meth, nonb_lp, gbr_lp, gbd_lp);
      launchValence(poly_vk, poly_rk, &ctrl_xe, &poly_psw, poly_auk, &tstw, &scw_final,
                    &vale_xe_tbr, EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                    acc_meth, vale_xe_lp);
      ctrl_xe.step += 1;
    }
    break;
  }

  // Advance the energy tracking history counter to log the final energy results
  sc->commit(devc_tier, gpu);
  sc->incrementSampleCount();
  sc->setLastTimeStep(total_steps);
  if (timer != nullptr) {
    cudaDeviceSynchronize();
    timer->assignTime(min_timings);
  }
}

//-------------------------------------------------------------------------------------------------
extern ScoreCard launchMinimization(const AtomGraphSynthesis &poly_ag,
                                    const StaticExclusionMaskSynthesis &poly_se,
                                    PhaseSpaceSynthesis *poly_ps, const MinimizeControls &mincon,
                                    const GpuDetails &gpu, const PrecisionModel prec,
                                    const int energy_accumulation_bits, StopWatch *timer,
                                    const std::string &task_name) {

  // Prepare to track the energies of the structures as they undergo geometry optimization.
  const int ntpr   = mincon.getDiagnosticPrintFrequency();
  const int nframe = (roundUp(mincon.getTotalCycles(), ntpr) / ntpr) + 1;
  ScoreCard result(poly_ps->getSystemCount(), nframe, energy_accumulation_bits);

  // Map out all kernels.  Only a few are needed but this is not a lot of work.
  const CoreKlManager launcher(gpu, poly_ag);
  const int2 vale_fe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::YES,
                                                        EvaluateEnergy::YES,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE, ClashResponse::NONE);
  const int2 vale_xe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::NO,
                                                        EvaluateEnergy::YES,
                                                        AccumulationMethod::SPLIT,
                                                        VwuGoal::ACCUMULATE, ClashResponse::NONE);
  const int2 vale_cdfe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::YES,
                                                          EvaluateEnergy::YES,
                                                          AccumulationMethod::SPLIT,
                                                          VwuGoal::ACCUMULATE,
                                                          ClashResponse::FORGIVE);
  const int2 vale_cdxe_lp = launcher.getValenceKernelDims(prec, EvaluateForce::NO,
                                                          EvaluateEnergy::YES,
                                                          AccumulationMethod::SPLIT,
                                                          VwuGoal::ACCUMULATE,
                                                          ClashResponse::FORGIVE);
  const NbwuKind nb_work_type = poly_ag.getNonbondedWorkType();
  const ImplicitSolventModel ism_type = poly_ag.getImplicitSolventModel();
  const int2 nonb_lp = launcher.getNonbondedKernelDims(prec, nb_work_type, EvaluateForce::YES,
                                                       EvaluateEnergy::YES,
                                                       AccumulationMethod::SPLIT, ism_type,
                                                       ClashResponse::NONE);
  const int2 nonb_cdlp = launcher.getNonbondedKernelDims(prec, nb_work_type, EvaluateForce::YES,
                                                         EvaluateEnergy::YES,
                                                         AccumulationMethod::SPLIT, ism_type,
                                                         ClashResponse::FORGIVE);
  const int2 redu_lp = launcher.getReductionKernelDims(prec, ReductionGoal::CONJUGATE_GRADIENT,
                                                       ReductionStage::ALL_REDUCE);

  // Prepare progress tracking objects.
  MolecularMechanicsControls mmctrl_fe(mincon);
  MolecularMechanicsControls mmctrl_xe(mincon);
  MolecularMechanicsControls mmctrl_cdfe(mincon);
  MolecularMechanicsControls mmctrl_cdxe(mincon);

  // Prepare cache space for each kernel
  CacheResource vale_fe_cache(vale_fe_lp.x, maximum_valence_work_unit_atoms);
  CacheResource vale_xe_cache(vale_xe_lp.x, maximum_valence_work_unit_atoms);
  CacheResource vale_cdfe_cache(vale_cdfe_lp.x, maximum_valence_work_unit_atoms);
  CacheResource vale_cdxe_cache(vale_cdxe_lp.x, maximum_valence_work_unit_atoms);
  CacheResource nonb_cache(nonb_lp.x, small_block_max_atoms);
  CacheResource nonb_cd_cache(nonb_cdlp.x, small_block_max_atoms);
  ImplicitSolventWorkspace ism_space(poly_ag.getSystemAtomOffsets(),
                                     poly_ag.getSystemAtomCounts(), prec);
  ReductionBridge poly_rbg(poly_ag.getReductionWorkUnitCount());
  LineMinimization line_record(poly_ag.getSystemCount());
  launchMinimization(prec, poly_ag, poly_se, poly_ps, mincon, &mmctrl_fe, &mmctrl_xe, &mmctrl_cdfe,
                     &mmctrl_cdxe, &result, &vale_fe_cache, &vale_xe_cache, &vale_cdfe_cache,
                     &vale_cdxe_cache, &nonb_cache, &nonb_cd_cache, &ism_space, &poly_rbg,
                     &line_record, chooseAccumulationMethod(poly_ps->getForceAccumulationBits()),
                     gpu, launcher, timer, task_name);
  return result;
}

//-------------------------------------------------------------------------------------------------
extern ScoreCard launchMinimization(AtomGraphSynthesis *poly_ag, PhaseSpaceSynthesis *poly_ps,
                                    const MinimizeControls &mincon, const GpuDetails &gpu,
                                    const PrecisionModel prec, const int energy_accumulation_bits,
                                    StopWatch *timer, const std::string &task_name) {
  switch (poly_ag->getNonbondedWorkType()) {
  case NbwuKind::TILE_GROUPS:
  case NbwuKind::SUPERTILES:
    {
      const StaticExclusionMaskSynthesis poly_se(poly_ag->getUniqueTopologies(),
                                                 poly_ag->getTopologyIndices());
      poly_ag->loadNonbondedWorkUnits(poly_se);
      return launchMinimization(*poly_ag, poly_se, poly_ps, mincon, gpu, prec,
                                energy_accumulation_bits, timer, task_name);
    }
    break;
  case NbwuKind::HONEYCOMB:
    break;
  }
  __builtin_unreachable();
}

} // namespace mm
} // namespace stormm
