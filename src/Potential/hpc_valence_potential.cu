// -*-c++-*-
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Accelerator/gpu_details.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/rounding.h"
#include "Numerics/numeric_enumerators.h"
#include "Numerics/split_fixed_precision.h"
#include "Potential/cellgrid.h"
#include "Potential/energy_enumerators.h"
#include "Synthesis/valence_workunit.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/trajectory_enumerators.h"
#include "hpc_valence_potential.h"

namespace stormm {
namespace energy {

using card::GpuDetails;
using card::CoreKlManager;
using constants::PrecisionModel;
using constants::large_block_size;
using constants::medium_block_size;
using constants::small_block_size;
using constants::twice_warp_bits_mask_int;
using constants::twice_warp_size_int;
using constants::warp_size_int;
using constants::warp_bits;
using constants::warp_bits_mask_int;
using numerics::chooseAccumulationMethod;
using numerics::AccumulationMethod;
using numerics::getEnumerationName;
using stmath::roundUp;
using symbols::asymptotic_to_one_f;
using symbols::asymptotic_to_one_lf;
using symbols::boltzmann_constant_f;
using symbols::gafs_to_kcal_f;
using symbols::inverse_one_minus_asymptote_f;
using symbols::inverse_one_minus_asymptote_lf;
using symbols::inverse_twopi_f;
using symbols::kcal_to_gafs_f;
using symbols::near_to_one_f;
using symbols::near_to_one_lf;
using symbols::pi;
using symbols::pi_f;
using symbols::twopi;
using symbols::twopi_f;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::half_valence_work_unit_atoms;
using synthesis::quarter_valence_work_unit_atoms;
using synthesis::eighth_valence_work_unit_atoms;
using synthesis::VwuAbstractMap;
using synthesis::VwuGoal;
using synthesis::vwu_abstract_length;
using trajectory::ThermostatKind;
using trajectory::ThermostatPartition;
using topology::TorsionKind;
using topology::VirtualSiteKind;

//-------------------------------------------------------------------------------------------------
#ifdef STORMM_USE_CUDA
extern cudaFuncAttributes queryValenceKernelRequirements(const PrecisionModel prec,
                                                         const EvaluateForce eval_frc,
                                                         const EvaluateEnergy eval_nrg,
                                                         const AccumulationMethod acc_meth,
                                                         const VwuGoal purpose,
                                                         const ClashResponse collision_handling,
                                                         const ValenceKernelSize kwidth,
                                                         const bool pme_compatible,
                                                         const PrecisionModel neighbor_prec,
                                                         const NeighborListKind neighbor_layout) {

  // Various details of the calculation lead into different branches and CUDA units.
  switch (prec) {
  case PrecisionModel::DOUBLE:
    if (pme_compatible) {
      switch (neighbor_layout) {
      case NeighborListKind::DUAL:
        return queryValenceKernelRequirementsDPMEDual(eval_nrg, collision_handling, neighbor_prec);
      case NeighborListKind::MONO:
        return queryValenceKernelRequirementsDPME(eval_nrg, collision_handling, neighbor_prec);
      }        
    }
    else {
      return queryValenceKernelRequirementsD(eval_frc, eval_nrg, purpose, collision_handling);
    }
    break;
  case PrecisionModel::SINGLE:
    if (pme_compatible) {
      switch (neighbor_layout) {
      case NeighborListKind::DUAL:
        return queryValenceKernelRequirementsFPMEDual(eval_nrg, acc_meth, collision_handling,
                                                      neighbor_prec);
      case NeighborListKind::MONO:
        return queryValenceKernelRequirementsFPME(eval_nrg, acc_meth, collision_handling,
                                                  neighbor_prec);
      }
    }
    else {
      switch (kwidth) {
      case ValenceKernelSize::XL:
        return queryValenceKernelRequirementsXL(eval_frc, eval_nrg, acc_meth, purpose,
                                                collision_handling);
      case ValenceKernelSize::LG:
        return queryValenceKernelRequirementsLG(eval_frc, eval_nrg, acc_meth, purpose,
                                                collision_handling);
      case ValenceKernelSize::MD:
        return queryValenceKernelRequirementsMD(eval_frc, eval_nrg, acc_meth, purpose,
                                                collision_handling);
      case ValenceKernelSize::SM:
        return queryValenceKernelRequirementsSM(eval_frc, eval_nrg, acc_meth, purpose,
                                                collision_handling);
      }
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes queryPmeValenceKernelRequirements(const PrecisionModel prec,
                                                            const PrecisionModel neighbor_prec,
                                                            const EvaluateForce eval_frc,
                                                            const EvaluateEnergy eval_nrg,
                                                            const AccumulationMethod acc_meth,
                                                            const VwuGoal purpose,
                                                            const ClashResponse clash_handling) {
  return queryValenceKernelRequirements(prec, eval_frc, eval_nrg, acc_meth, purpose,
                                        clash_handling, ValenceKernelSize::XL, true,
                                        neighbor_prec);
}
#endif

//-------------------------------------------------------------------------------------------------
extern void launchValence(const SyValenceKit<float> &poly_vk,
                          const SyRestraintKit<float, float2, float4> &poly_rk,
                          MMControlKit<float> *ctrl, PsSynthesisWriter *poly_psw,
	                  const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                          ThermostatWriter<float> *tstw, ScoreCardWriter *scw,
                          CacheResourceKit<float> *gmem_r, const EvaluateForce eval_force,
                          const EvaluateEnergy eval_energy, const VwuGoal purpose,
                          const AccumulationMethod force_sum, const int2 bt,
	                  const float clash_distance, const float clash_ratio) {
  AccumulationMethod refined_force_sum;
  switch (force_sum) {
  case AccumulationMethod::SPLIT:
  case AccumulationMethod::WHOLE:
    refined_force_sum = force_sum;
    break;
  case AccumulationMethod::AUTOMATIC:
    refined_force_sum = chooseAccumulationMethod(poly_psw->frc_bits);
    break;
  }

  // Use the launch grid to determine what size the kernel is.  This requires that all "XL" kernels
  // have > 256 threads per block, "LG" kernels have > 128 threads per block, and "MD" kernels
  // have > 64 threads per block.  All "SM" kernels launch with 64 threads per block.
  if (bt.y > 256) {
    launchValenceXL(poly_vk, poly_rk, ctrl, poly_psw, poly_auk, tstw, scw, gmem_r, eval_force,
                    eval_energy, purpose, refined_force_sum, bt, clash_distance, clash_ratio);
  }
  else if (bt.y > 128) {
    launchValenceLG(poly_vk, poly_rk, ctrl, poly_psw, poly_auk, tstw, scw, gmem_r, eval_force,
                    eval_energy, purpose, refined_force_sum, bt, clash_distance, clash_ratio);
  }
  else if (bt.y > 64) {
    launchValenceMD(poly_vk, poly_rk, ctrl, poly_psw, poly_auk, tstw, scw, gmem_r, eval_force,
                    eval_energy, purpose, refined_force_sum, bt, clash_distance, clash_ratio);
  }
  else {
    launchValenceSM(poly_vk, poly_rk, ctrl, poly_psw, poly_auk, tstw, scw, gmem_r, eval_force,
                    eval_energy, purpose, refined_force_sum, bt, clash_distance, clash_ratio);
  }
}
  
//-------------------------------------------------------------------------------------------------
extern void launchValence(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                          MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                          Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                          const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                          const VwuGoal purpose, const AccumulationMethod force_sum,
                          const CoreKlManager &launcher, const double clash_distance,
                          const double clash_ratio) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  const ClashResponse mitigate_clash = (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) ?
                                       ClashResponse::FORGIVE : ClashResponse::NONE;
  const int2 bt = launcher.getValenceKernelDims(prec, eval_force, eval_energy,
                                                AccumulationMethod::SPLIT, purpose,
                                                mitigate_clash);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(tier);
      const SyRestraintKit<double, double2, double4_16a> poly_rk =
        poly_ag.getDoublePrecisionRestraintKit(tier);
      const SyAtomUpdateKit<double, double2, double4_16a> poly_auk =
        poly_ag.getDoublePrecisionAtomUpdateKit(tier);
      MMControlKit<double> ctrl = mmctrl->dpData(tier);
      ThermostatWriter tstw = heat_bath->dpData(tier);
      CacheResourceKit<double> gmem_r = tb_space->dpData(tier);
      launchValence(poly_vk, poly_rk, &ctrl, &poly_psw, poly_auk, &tstw, &scw, &gmem_r, eval_force,
                    eval_energy, purpose, bt, clash_distance, clash_ratio);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(tier);
      const SyRestraintKit<float, float2, float4> poly_rk =
        poly_ag.getSinglePrecisionRestraintKit(tier);
      const SyAtomUpdateKit<float, float2, float4> poly_auk =
        poly_ag.getSinglePrecisionAtomUpdateKit(tier);
      MMControlKit<float> ctrl = mmctrl->spData(tier);
      ThermostatWriter tstw = heat_bath->spData(tier);
      CacheResourceKit<float> gmem_r = tb_space->spData(tier);
      launchValence(poly_vk, poly_rk, &ctrl, &poly_psw, poly_auk, &tstw, &scw, &gmem_r, eval_force,
                    eval_energy, purpose, force_sum, bt, clash_distance, clash_ratio);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchValence(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                          MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                          Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                          const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                          const VwuGoal purpose, const CoreKlManager &launcher,
                          const double clash_distance, const double clash_ratio) {
  if (prec == PrecisionModel::DOUBLE || poly_ps->getForceAccumulationBits() <= 24) {
    launchValence(prec, poly_ag, mmctrl, poly_ps, heat_bath, sc, tb_space, eval_force, eval_energy,
                  purpose, AccumulationMethod::SPLIT, launcher, clash_distance, clash_ratio);
  }
  else {
    launchValence(prec, poly_ag, mmctrl, poly_ps, heat_bath, sc, tb_space, eval_force, eval_energy,
                  purpose, AccumulationMethod::WHOLE, launcher, clash_distance, clash_ratio);
  }
}

//-------------------------------------------------------------------------------------------------
void launchValence(PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                   const CellGrid<double, llint, double, double4_16a> &cg,
                   MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                   Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                   const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                   const VwuGoal purpose, const AccumulationMethod force_sum,
                   const CoreKlManager &launcher, const double clash_distance,
                   const double clash_ratio) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  const ClashResponse mitigate_clash = (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) ?
                                       ClashResponse::FORGIVE : ClashResponse::NONE;
  const int2 bt = launcher.getValenceKernelDims(prec, eval_force, eval_energy,
                                                AccumulationMethod::SPLIT, purpose,
                                                mitigate_clash);
  const CellGridReader<double, llint, double, double4_16a> cgr = cg.data(tier);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(tier);
      const SyRestraintKit<double, double2, double4_16a> poly_rk =
        poly_ag.getDoublePrecisionRestraintKit(tier);
      const SyAtomUpdateKit<double, double2, double4_16a> poly_auk =
        poly_ag.getDoublePrecisionAtomUpdateKit(tier);
      MMControlKit<double> ctrl = mmctrl->dpData(tier);
      ThermostatWriter tstw = heat_bath->dpData(tier);
      CacheResourceKit<double> gmem_r = tb_space->dpData(tier);
      launchValence(poly_vk, poly_rk, cgr, &ctrl, &poly_psw, poly_auk, &tstw, &scw, &gmem_r,
                    eval_force, eval_energy, purpose, bt, clash_distance, clash_ratio);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(tier);
      const SyRestraintKit<float, float2, float4> poly_rk =
        poly_ag.getSinglePrecisionRestraintKit(tier);
      const SyAtomUpdateKit<float, float2, float4> poly_auk =
        poly_ag.getSinglePrecisionAtomUpdateKit(tier);
      MMControlKit<float> ctrl = mmctrl->spData(tier);
      ThermostatWriter tstw = heat_bath->spData(tier);
      CacheResourceKit<float> gmem_r = tb_space->spData(tier);
      launchValence(poly_vk, poly_rk, cgr, &ctrl, &poly_psw, poly_auk, &tstw, &scw, &gmem_r,
                    eval_force, eval_energy, purpose, force_sum, bt, clash_distance, clash_ratio);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchValence(PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                   const CellGrid<double, llint, double, double4_16a> &cg,
                   MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                   Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                   const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                   const VwuGoal purpose, const CoreKlManager &launcher,
                   const double clash_distance, const double clash_ratio) {
  if (prec == PrecisionModel::DOUBLE || poly_ps->getForceAccumulationBits() <= 24) {
    launchValence(prec, poly_ag, cg, mmctrl, poly_ps, heat_bath, sc, tb_space, eval_force,
                  eval_energy, purpose, AccumulationMethod::SPLIT, launcher, clash_distance,
                  clash_ratio);
  }
  else {
    launchValence(prec, poly_ag, cg, mmctrl, poly_ps, heat_bath, sc, tb_space, eval_force,
                  eval_energy, purpose, AccumulationMethod::WHOLE, launcher, clash_distance,
                  clash_ratio);
  }
}

//-------------------------------------------------------------------------------------------------
void launchValence(PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                   const CellGrid<double, llint, double, double4_16a> &cg_qq,
                   const CellGrid<double, llint, double, double4_16a> &cg_lj,
                   MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                   Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                   const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                   const VwuGoal purpose, const AccumulationMethod force_sum,
                   const CoreKlManager &launcher, const double clash_distance,
                   const double clash_ratio) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  const ClashResponse mitigate_clash = (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) ?
                                       ClashResponse::FORGIVE : ClashResponse::NONE;
  const int2 bt = launcher.getValenceKernelDims(prec, eval_force, eval_energy,
                                                AccumulationMethod::SPLIT, purpose,
                                                mitigate_clash);
  const CellGridReader<double, llint, double, double4_16a> cgr_qq = cg_qq.data(tier);
  const CellGridReader<double, llint, double, double4_16a> cgr_lj = cg_lj.data(tier);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(tier);
      const SyRestraintKit<double, double2, double4_16a> poly_rk =
        poly_ag.getDoublePrecisionRestraintKit(tier);
      const SyAtomUpdateKit<double, double2, double4_16a> poly_auk =
        poly_ag.getDoublePrecisionAtomUpdateKit(tier);
      MMControlKit<double> ctrl = mmctrl->dpData(tier);
      ThermostatWriter tstw = heat_bath->dpData(tier);
      CacheResourceKit<double> gmem_r = tb_space->dpData(tier);
      launchValence(poly_vk, poly_rk, cgr_qq, cgr_lj, &ctrl, &poly_psw, poly_auk, &tstw, &scw,
                    &gmem_r, eval_force, eval_energy, purpose, bt, clash_distance, clash_ratio);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(tier);
      const SyRestraintKit<float, float2, float4> poly_rk =
        poly_ag.getSinglePrecisionRestraintKit(tier);
      const SyAtomUpdateKit<float, float2, float4> poly_auk =
        poly_ag.getSinglePrecisionAtomUpdateKit(tier);
      MMControlKit<float> ctrl = mmctrl->spData(tier);
      ThermostatWriter tstw = heat_bath->spData(tier);
      CacheResourceKit<float> gmem_r = tb_space->spData(tier);
      launchValence(poly_vk, poly_rk, cgr_qq, cgr_lj, &ctrl, &poly_psw, poly_auk, &tstw, &scw,
                    &gmem_r, eval_force, eval_energy, purpose, force_sum, bt, clash_distance,
                    clash_ratio);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchValence(PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                   const CellGrid<double, llint, double, double4_16a> &cg_qq,
                   const CellGrid<double, llint, double, double4_16a> &cg_lj,
                   MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                   Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                   const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                   const VwuGoal purpose, const CoreKlManager &launcher,
                   const double clash_distance, const double clash_ratio) {
  if (prec == PrecisionModel::DOUBLE || poly_ps->getForceAccumulationBits() <= 24) {
    launchValence(prec, poly_ag, cg_qq, cg_lj,  mmctrl, poly_ps, heat_bath, sc, tb_space,
                  eval_force, eval_energy, purpose, AccumulationMethod::SPLIT, launcher,
                  clash_distance, clash_ratio);
  }
  else {
    launchValence(prec, poly_ag, cg_qq, cg_lj, mmctrl, poly_ps, heat_bath, sc, tb_space,
                  eval_force, eval_energy, purpose, AccumulationMethod::WHOLE, launcher,
                  clash_distance, clash_ratio);
  }
}

//-------------------------------------------------------------------------------------------------
void launchValence(PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                   const CellGrid<float, int, float, float4> &cg,
                   MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                   Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                   const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                   const VwuGoal purpose, const AccumulationMethod force_sum,
                   const CoreKlManager &launcher, const double clash_distance,
                   const double clash_ratio) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  const ClashResponse mitigate_clash = (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) ?
                                       ClashResponse::FORGIVE : ClashResponse::NONE;
  const int2 bt = launcher.getValenceKernelDims(prec, eval_force, eval_energy,
                                                AccumulationMethod::SPLIT, purpose,
                                                mitigate_clash);
  const CellGridReader<float, int, float, float4> cgr = cg.data(tier);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(tier);
      const SyRestraintKit<double, double2, double4_16a> poly_rk =
        poly_ag.getDoublePrecisionRestraintKit(tier);
      const SyAtomUpdateKit<double, double2, double4_16a> poly_auk =
        poly_ag.getDoublePrecisionAtomUpdateKit(tier);
      MMControlKit<double> ctrl = mmctrl->dpData(tier);
      ThermostatWriter tstw = heat_bath->dpData(tier);
      CacheResourceKit<double> gmem_r = tb_space->dpData(tier);
      launchValence(poly_vk, poly_rk, cgr, &ctrl, &poly_psw, poly_auk, &tstw, &scw, &gmem_r,
                    eval_force, eval_energy, purpose, bt, clash_distance, clash_ratio);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(tier);
      const SyRestraintKit<float, float2, float4> poly_rk =
        poly_ag.getSinglePrecisionRestraintKit(tier);
      const SyAtomUpdateKit<float, float2, float4> poly_auk =
        poly_ag.getSinglePrecisionAtomUpdateKit(tier);
      MMControlKit<float> ctrl = mmctrl->spData(tier);
      ThermostatWriter tstw = heat_bath->spData(tier);
      CacheResourceKit<float> gmem_r = tb_space->spData(tier);
      launchValence(poly_vk, poly_rk, cgr, &ctrl, &poly_psw, poly_auk, &tstw, &scw, &gmem_r,
                    eval_force, eval_energy, purpose, force_sum, bt, clash_distance, clash_ratio);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchValence(PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                   const CellGrid<float, int, float, float4> &cg,
                   MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                   Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                   const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                   const VwuGoal purpose, const CoreKlManager &launcher,
                   const double clash_distance, const double clash_ratio) {
  if (prec == PrecisionModel::DOUBLE || poly_ps->getForceAccumulationBits() <= 24) {
    launchValence(prec, poly_ag, cg, mmctrl, poly_ps, heat_bath, sc, tb_space, eval_force,
                  eval_energy, purpose, AccumulationMethod::SPLIT, launcher, clash_distance,
                  clash_ratio);
  }
  else {
    launchValence(prec, poly_ag, cg, mmctrl, poly_ps, heat_bath, sc, tb_space, eval_force,
                  eval_energy, purpose, AccumulationMethod::WHOLE, launcher, clash_distance,
                  clash_ratio);
  }
}

//-------------------------------------------------------------------------------------------------
void launchValence(PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                   const CellGrid<float, int, float, float4> &cg_qq,
                   const CellGrid<float, int, float, float4> &cg_lj,
                   MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                   Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                   const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                   const VwuGoal purpose, const AccumulationMethod force_sum,
                   const CoreKlManager &launcher, const double clash_distance,
                   const double clash_ratio) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  const ClashResponse mitigate_clash = (clash_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) ?
                                       ClashResponse::FORGIVE : ClashResponse::NONE;
  const int2 bt = launcher.getValenceKernelDims(prec, eval_force, eval_energy,
                                                AccumulationMethod::SPLIT, purpose,
                                                mitigate_clash);
  const CellGridReader<float, int, float, float4> cgr_qq = cg_qq.data(tier);
  const CellGridReader<float, int, float, float4> cgr_lj = cg_lj.data(tier);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(tier);
      const SyRestraintKit<double, double2, double4_16a> poly_rk =
        poly_ag.getDoublePrecisionRestraintKit(tier);
      const SyAtomUpdateKit<double, double2, double4_16a> poly_auk =
        poly_ag.getDoublePrecisionAtomUpdateKit(tier);
      MMControlKit<double> ctrl = mmctrl->dpData(tier);
      ThermostatWriter<double> tstw = heat_bath->dpData(tier);
      CacheResourceKit<double> gmem_r = tb_space->dpData(tier);
      launchValence(poly_vk, poly_rk, cgr_qq, cgr_lj, &ctrl, &poly_psw, poly_auk, &tstw, &scw,
                    &gmem_r, eval_force, eval_energy, purpose, bt, clash_distance, clash_ratio);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(tier);
      const SyRestraintKit<float, float2, float4> poly_rk =
        poly_ag.getSinglePrecisionRestraintKit(tier);
      const SyAtomUpdateKit<float, float2, float4> poly_auk =
        poly_ag.getSinglePrecisionAtomUpdateKit(tier);
      MMControlKit<float> ctrl = mmctrl->spData(tier);
      ThermostatWriter<float> tstw = heat_bath->spData(tier);
      CacheResourceKit<float> gmem_r = tb_space->spData(tier);
      launchValence(poly_vk, poly_rk, cgr_qq, cgr_lj, &ctrl, &poly_psw, poly_auk, &tstw, &scw,
                    &gmem_r, eval_force, eval_energy, purpose, force_sum, bt, clash_distance,
                    clash_ratio);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchValence(PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                   const CellGrid<float, int, float, float4> &cg_qq,
                   const CellGrid<float, int, float, float4> &cg_lj,
                   MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                   Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                   const EvaluateForce eval_force, const EvaluateEnergy eval_energy,
                   const VwuGoal purpose, const CoreKlManager &launcher,
                   const double clash_distance, const double clash_ratio) {
  if (prec == PrecisionModel::DOUBLE || poly_ps->getForceAccumulationBits() <= 24) {
    launchValence(prec, poly_ag, cg_qq, cg_lj,  mmctrl, poly_ps, heat_bath, sc, tb_space,
                  eval_force, eval_energy, purpose, AccumulationMethod::SPLIT, launcher,
                  clash_distance, clash_ratio);
  }
  else {
    launchValence(prec, poly_ag, cg_qq, cg_lj, mmctrl, poly_ps, heat_bath, sc, tb_space,
                  eval_force, eval_energy, purpose, AccumulationMethod::WHOLE, launcher,
                  clash_distance, clash_ratio);
  }
}

} // namespace energy
} // namespace stormm
