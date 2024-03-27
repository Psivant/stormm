#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "Numerics/split_fixed_precision.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Trajectory/trajectory_enumerators.h"
#include "mm_controls.h"

namespace stormm {
namespace mm {

using card::HybridKind;
using stmath::ReductionGoal;
using numerics::AccumulationMethod;
using topology::UnitCellType;
using trajectory::IntegrationStage;
  
//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls::MolecularMechanicsControls(const double time_step_in,
                                                       const double rattle_tol_in,
                                                       const double initial_step_in,
                                                       const int sd_cycles_in,
                                                       const int max_cycles_in) :
  step_number{0}, sd_cycles{sd_cycles_in}, max_cycles{max_cycles_in}, time_step{time_step_in},
  rattle_tol{rattle_tol_in}, initial_step{initial_step_in},
  vwu_progress{HybridKind::POINTER, "mm_vwu_counters"},
  velocity_update_progress{HybridKind::POINTER, "mm_vupt_counters"},
  velocity_constraint_progress{HybridKind::POINTER, "mm_vcns_counters"},
  position_update_progress{HybridKind::POINTER, "mm_pupt_counters"},
  geometry_constraint_progress{HybridKind::POINTER, "mm_gcns_counters"},
  nbwu_progress{HybridKind::POINTER, "mm_nbwu_counters"},
  pmewu_progress{HybridKind::POINTER, "mm_pmewu_counters"},
  gbrwu_progress{HybridKind::POINTER, "mm_gbrwu_counters"},
  gbdwu_progress{HybridKind::POINTER, "mm_gbdwu_counters"},
  gather_wu_progress{HybridKind::POINTER, "mm_gtwu_counters"},
  scatter_wu_progress{HybridKind::POINTER, "mm_scwu_counters"},
  all_reduce_wu_progress{HybridKind::POINTER, "mm_rdwu_counters"},
  int_data{24 * warp_size_int, "work_unit_prog_data"}
{
  vwu_progress.setPointer(&int_data,                                  0, 2 * warp_size_int);
  velocity_update_progress.setPointer(&int_data,      2 * warp_size_int, 2 * warp_size_int);
  velocity_constraint_progress.setPointer(&int_data,  4 * warp_size_int, 2 * warp_size_int);
  position_update_progress.setPointer(&int_data,      6 * warp_size_int, 2 * warp_size_int);
  geometry_constraint_progress.setPointer(&int_data,  8 * warp_size_int, 2 * warp_size_int);
  nbwu_progress.setPointer(&int_data,                10 * warp_size_int, 2 * warp_size_int);
  pmewu_progress.setPointer(&int_data,               12 * warp_size_int, 2 * warp_size_int);
  gbrwu_progress.setPointer(&int_data,               14 * warp_size_int, 2 * warp_size_int);
  gbdwu_progress.setPointer(&int_data,               16 * warp_size_int, 2 * warp_size_int);
  gather_wu_progress.setPointer(&int_data,           18 * warp_size_int, 2 * warp_size_int);
  scatter_wu_progress.setPointer(&int_data,          20 * warp_size_int, 2 * warp_size_int);
  all_reduce_wu_progress.setPointer(&int_data,       22 * warp_size_int, 2 * warp_size_int);
}

//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls::MolecularMechanicsControls(const DynamicsControls &user_input) :
    MolecularMechanicsControls(user_input.getTimeStep(), user_input.getRattleTolerance(),
                               default_minimize_dx0, default_minimize_ncyc,
                               user_input.getStepCount())
{}

//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls::MolecularMechanicsControls(const MinimizeControls &user_input) :
    MolecularMechanicsControls(default_dynamics_time_step, default_rattle_tolerance,
                               user_input.getInitialStep(), user_input.getSteepestDescentCycles(),
                               user_input.getTotalCycles())
{}

//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls::
MolecularMechanicsControls(const MolecularMechanicsControls &original) :
  step_number{original.step_number},
  sd_cycles{original.sd_cycles},
  max_cycles{original.max_cycles},
  time_step{original.time_step},
  rattle_tol{original.rattle_tol},
  initial_step{original.initial_step},
  vwu_progress{original.vwu_progress},
  velocity_update_progress{original.velocity_update_progress},
  velocity_constraint_progress{original.velocity_constraint_progress},
  position_update_progress{original.position_update_progress},
  geometry_constraint_progress{original.geometry_constraint_progress},
  nbwu_progress{original.nbwu_progress},
  pmewu_progress{original.pmewu_progress},
  gbrwu_progress{original.gbrwu_progress},
  gbdwu_progress{original.gbdwu_progress},
  gather_wu_progress{original.gather_wu_progress},
  scatter_wu_progress{original.scatter_wu_progress},
  all_reduce_wu_progress{original.all_reduce_wu_progress},
  int_data{original.int_data}
{
  rebasePointers();
}

//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls&
MolecularMechanicsControls::operator=(const MolecularMechanicsControls &other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  step_number = other.step_number;
  sd_cycles = other.sd_cycles;
  max_cycles = other.max_cycles;
  time_step = other.time_step;
  rattle_tol = other.rattle_tol;
  initial_step = other.initial_step;
  vwu_progress = other.vwu_progress;
  velocity_update_progress = other.velocity_update_progress;
  velocity_constraint_progress = other.velocity_constraint_progress;
  position_update_progress = other.position_update_progress;
  geometry_constraint_progress = other.geometry_constraint_progress;
  nbwu_progress = other.nbwu_progress;
  pmewu_progress = other.pmewu_progress;
  gbrwu_progress = other.gbrwu_progress;
  gbdwu_progress = other.gbdwu_progress;
  gather_wu_progress = other.gather_wu_progress;
  scatter_wu_progress = other.scatter_wu_progress;
  all_reduce_wu_progress = other.all_reduce_wu_progress;
  int_data = other.int_data;

  // Repair pointers and return the result
  rebasePointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls::MolecularMechanicsControls(MolecularMechanicsControls &&original) :
  step_number{original.step_number},
  sd_cycles{original.sd_cycles},
  max_cycles{original.max_cycles},
  time_step{original.time_step},
  rattle_tol{original.rattle_tol},
  initial_step{original.initial_step},
  vwu_progress{std::move(original.vwu_progress)},
  velocity_update_progress{std::move(original.velocity_update_progress)},
  velocity_constraint_progress{std::move(original.velocity_constraint_progress)},
  position_update_progress{std::move(original.position_update_progress)},
  geometry_constraint_progress{std::move(original.geometry_constraint_progress)},
  nbwu_progress{std::move(original.nbwu_progress)},
  pmewu_progress{std::move(original.pmewu_progress)},
  gbrwu_progress{std::move(original.gbrwu_progress)},
  gbdwu_progress{std::move(original.gbdwu_progress)},
  gather_wu_progress{std::move(original.gather_wu_progress)},
  scatter_wu_progress{std::move(original.scatter_wu_progress)},
  all_reduce_wu_progress{std::move(original.all_reduce_wu_progress)},
  int_data{std::move(original.int_data)}
{}

//-------------------------------------------------------------------------------------------------
MolecularMechanicsControls&
MolecularMechanicsControls::operator=(MolecularMechanicsControls &&other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  step_number = other.step_number;
  sd_cycles = other.sd_cycles;
  max_cycles = other.max_cycles;
  time_step = other.time_step;
  rattle_tol = other.rattle_tol;
  initial_step = other.initial_step;
  vwu_progress = std::move(other.vwu_progress);
  velocity_update_progress = std::move(other.velocity_update_progress);
  velocity_constraint_progress = std::move(other.velocity_constraint_progress);
  position_update_progress = std::move(other.position_update_progress);
  geometry_constraint_progress = std::move(other.geometry_constraint_progress);
  nbwu_progress = std::move(other.nbwu_progress);
  pmewu_progress = std::move(other.pmewu_progress);
  gbrwu_progress = std::move(other.gbrwu_progress);
  gbdwu_progress = std::move(other.gbdwu_progress);
  gather_wu_progress = std::move(other.gather_wu_progress);
  scatter_wu_progress = std::move(other.scatter_wu_progress);
  all_reduce_wu_progress = std::move(other.all_reduce_wu_progress);
  int_data = std::move(other.int_data);
  return *this;
}

//-------------------------------------------------------------------------------------------------
int MolecularMechanicsControls::getStepNumber() const {
  return step_number;
}

//-------------------------------------------------------------------------------------------------
int MolecularMechanicsControls::getSteepestDescentCycles() const {
  return sd_cycles;
}

//-------------------------------------------------------------------------------------------------
int MolecularMechanicsControls::getTotalCycles() const {
  return max_cycles;
}

//-------------------------------------------------------------------------------------------------
double MolecularMechanicsControls::getTimeStep() const {
  return time_step;
}

//-------------------------------------------------------------------------------------------------
double MolecularMechanicsControls::getRattleTolerance() const {
  return rattle_tol;
}

//-------------------------------------------------------------------------------------------------
double MolecularMechanicsControls::getInitialMinimizationStep() const {
  return initial_step;
}

//-------------------------------------------------------------------------------------------------
int MolecularMechanicsControls::getValenceWorkUnitProgress(const int counter_index,
                                                           const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return vwu_progress.readHost(counter_index);    
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return vwu_progress.readDevice(counter_index);
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int MolecularMechanicsControls::getNonbondedWorkUnitProgress(const int counter_index,
                                                             const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return nbwu_progress.readHost(counter_index);    
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return nbwu_progress.readDevice(counter_index);
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int MolecularMechanicsControls::getPmeWorkUnitProgress(const int counter_index,
                                                       const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    return pmewu_progress.readHost(counter_index);    
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return pmewu_progress.readDevice(counter_index);
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int MolecularMechanicsControls::getReductionWorkUnitProgress(const int counter_index,
                                                             const ReductionStage process,
                                                             const HybridTargetLevel tier) const {
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (process) {
    case ReductionStage::GATHER:
      return gather_wu_progress.readHost(counter_index);
    case ReductionStage::SCATTER:
    case ReductionStage::RESCALE:
      return scatter_wu_progress.readHost(counter_index);
    case ReductionStage::ALL_REDUCE:
      return all_reduce_wu_progress.readHost(counter_index);
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (process) {
    case ReductionStage::GATHER:
      return gather_wu_progress.readDevice(counter_index);
    case ReductionStage::SCATTER:
    case ReductionStage::RESCALE:
      return scatter_wu_progress.readDevice(counter_index);
    case ReductionStage::ALL_REDUCE:
      return all_reduce_wu_progress.readDevice(counter_index);
    }
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
MMControlKit<double> MolecularMechanicsControls::dpData(const HybridTargetLevel tier) {
  return MMControlKit<double>(step_number, sd_cycles, max_cycles, time_step, rattle_tol,
                              initial_step, vwu_progress.data(tier),
                              velocity_update_progress.data(tier),
                              velocity_constraint_progress.data(tier),
                              position_update_progress.data(tier),
                              geometry_constraint_progress.data(tier), nbwu_progress.data(tier),
                              pmewu_progress.data(tier), gbrwu_progress.data(tier),
                              gbdwu_progress.data(tier), gather_wu_progress.data(tier),
                              scatter_wu_progress.data(tier), all_reduce_wu_progress.data(tier));
}

//-------------------------------------------------------------------------------------------------
MMControlKit<float> MolecularMechanicsControls::spData(const HybridTargetLevel tier) {
  return MMControlKit<float>(step_number, sd_cycles, max_cycles, time_step, rattle_tol,
                             initial_step, vwu_progress.data(tier),
                             velocity_update_progress.data(tier),
                             velocity_constraint_progress.data(tier),
                             position_update_progress.data(tier),
                             geometry_constraint_progress.data(tier), nbwu_progress.data(tier),
                             pmewu_progress.data(tier), gbrwu_progress.data(tier),
                             gbdwu_progress.data(tier), gather_wu_progress.data(tier),
                             scatter_wu_progress.data(tier), all_reduce_wu_progress.data(tier));
}

//-------------------------------------------------------------------------------------------------
void MolecularMechanicsControls::primeWorkUnitCounters(const CoreKlManager &launcher,
                                                       const EvaluateForce eval_frc,
                                                       const EvaluateEnergy eval_nrg,
                                                       const ClashResponse softcore,
                                                       const VwuGoal purpose,
                                                       const PrecisionModel valence_prec,
                                                       const PrecisionModel nonbond_prec,
                                                       const QMapMethod qspread_approach,
                                                       const PrecisionModel acc_prec,
                                                       const size_t image_coord_type,
                                                       const int qspread_order,
                                                       const AtomGraphSynthesis &poly_ag) {
  const GpuDetails wgpu = launcher.getGpu();
  const ImplicitSolventModel igb = poly_ag.getImplicitSolventModel();
  const int arch_major = wgpu.getArchMajor();
  const int arch_minor = wgpu.getArchMinor();
  const int smp_count = wgpu.getSMPCount();

  // The numbers of blocks that will be launched in each grid are critical for priming the work
  // unit progress counters.  As such, they (the grid launch size) should be consistent across
  // different variants of each kernel for a given precision level, even if the exact numbers of
  // threads per block have to vary based on what each block of that kernel variant is required to
  // do.  Here, it should suffice to query the launch parameters of just one of the blocks.
  const int2 vwu_lp = launcher.getValenceKernelDims(valence_prec, eval_frc, eval_nrg,
                                                    AccumulationMethod::SPLIT, purpose, softcore);
#if 0
  const int2 vupt_lp = launcher.getIntegrationKernelDims(valence_prec, AccumulationMethod::SPLIT,
                                                         IntegrationStage::VELOCITY_UPDATE);
#endif
  const int2 vcns_lp = launcher.getIntegrationKernelDims(valence_prec, AccumulationMethod::SPLIT,
                                                         IntegrationStage::VELOCITY_CONSTRAINT);
#if 0
  const int2 pupt_lp = launcher.getIntegrationKernelDims(valence_prec, AccumulationMethod::SPLIT,
                                                         IntegrationStage::POSITION_UPDATE);
#endif
  const int2 gcns_lp = launcher.getIntegrationKernelDims(valence_prec, AccumulationMethod::SPLIT,
                                                         IntegrationStage::GEOMETRY_CONSTRAINT);
  switch (poly_ag.getUnitCellType()) {
  case UnitCellType::NONE:
    {
      const int2 gbrwu_lp = launcher.getBornRadiiKernelDims(nonbond_prec,
                                                            poly_ag.getNonbondedWorkType(),
                                                            AccumulationMethod::SPLIT, igb);
      const int2 gbdwu_lp = launcher.getBornDerivativeKernelDims(nonbond_prec,
                                                                 poly_ag.getNonbondedWorkType(),
                                                                 AccumulationMethod::SPLIT, igb);
      const int2 nbwu_lp = launcher.getNonbondedKernelDims(nonbond_prec,
                                                           poly_ag.getNonbondedWorkType(),
                                                           eval_frc, eval_nrg,
                                                           AccumulationMethod::SPLIT, igb,
                                                           softcore);
      for (int i = 0; i < twice_warp_size_int; i++) {
        nbwu_progress.putHost(nbwu_lp.x, i);
        gbrwu_progress.putHost(gbrwu_lp.x, i);
        gbdwu_progress.putHost(gbdwu_lp.x, i);
      }
    }
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    {
      // The launch grid for the density mapping kernel will be identical for short- and long-form
      // fixed precision accumulations.
      const int2 pmewu_lp = launcher.getDensityMappingKernelDims(qspread_approach, nonbond_prec,
                                                                 acc_prec, true, image_coord_type,
                                                                 qspread_order);
      for (int i = 0; i < twice_warp_size_int; i++) {
        pmewu_progress.putHost(pmewu_lp.x, i);
      }
    }
    break;
  }
  const int2 rdwu_lp = launcher.getReductionKernelDims(valence_prec,
                                                       ReductionGoal::CONJUGATE_GRADIENT,
                                                       ReductionStage::ALL_REDUCE);
  const int vwu_block_count  = vwu_lp.x;
#if 0
  const int vupt_block_count = vupt_lp.x;
#endif
  const int vcns_block_count = vcns_lp.x;
#if 0
  const int pupt_block_count = pupt_lp.x;
#endif
  const int gcns_block_count = gcns_lp.x;
  const int gtwu_block_count = rdwu_lp.x;
  const int scwu_block_count = rdwu_lp.x;
  const int rdwu_block_count = rdwu_lp.x;
  for (int i = 0; i < twice_warp_size_int; i++) {
    vwu_progress.putHost(vwu_block_count, i);
    velocity_constraint_progress.putHost(vcns_block_count, i);
    geometry_constraint_progress.putHost(gcns_block_count, i);
    gather_wu_progress.putHost(gtwu_block_count, i);
    scatter_wu_progress.putHost(scwu_block_count, i);
    all_reduce_wu_progress.putHost(rdwu_block_count, i);
  }
#ifdef STORMM_USE_HPC
  upload();
#endif
}

//-------------------------------------------------------------------------------------------------
void MolecularMechanicsControls::primeWorkUnitCounters(const CoreKlManager &launcher,
                                                       const EvaluateForce eval_frc,
                                                       const EvaluateEnergy eval_nrg,
                                                       const ClashResponse softcore,
                                                       const VwuGoal purpose,
                                                       const PrecisionModel valence_prec,
                                                       const PrecisionModel nonbond_prec,
                                                       const AtomGraphSynthesis &poly_ag) {
  primeWorkUnitCounters(launcher, eval_frc, eval_nrg, softcore, purpose, valence_prec,
                        nonbond_prec, QMapMethod::GENERAL_PURPOSE, valence_prec, int_type_index,
                        4, poly_ag);
}

//-------------------------------------------------------------------------------------------------
void MolecularMechanicsControls::primeWorkUnitCounters(const CoreKlManager &launcher,
                                                       const EvaluateForce eval_frc,
                                                       const EvaluateEnergy eval_nrg,
                                                       const VwuGoal purpose,
                                                       const PrecisionModel valence_prec,
                                                       const PrecisionModel nonbond_prec,
                                                       const AtomGraphSynthesis &poly_ag) {
  primeWorkUnitCounters(launcher, eval_frc, eval_nrg, ClashResponse::NONE, purpose, valence_prec,
                        nonbond_prec, QMapMethod::GENERAL_PURPOSE, valence_prec, int_type_index,
                        4, poly_ag);
}

//-------------------------------------------------------------------------------------------------
void MolecularMechanicsControls::primeWorkUnitCounters(const CoreKlManager &launcher,
                                                       const EvaluateForce eval_frc,
                                                       const EvaluateEnergy eval_nrg,
                                                       const VwuGoal purpose,
                                                       const PrecisionModel general_prec,
                                                       const AtomGraphSynthesis &poly_ag) {
  primeWorkUnitCounters(launcher, eval_frc, eval_nrg, ClashResponse::NONE, purpose, general_prec,
                        general_prec, QMapMethod::GENERAL_PURPOSE, general_prec, int_type_index,
                        4, poly_ag);
}

//-------------------------------------------------------------------------------------------------
void MolecularMechanicsControls::incrementStep() {
  step_number += 1;
}
  
#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void MolecularMechanicsControls::upload() {
  int_data.upload();
}

//-------------------------------------------------------------------------------------------------
void MolecularMechanicsControls::download() {
  int_data.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void MolecularMechanicsControls::rebasePointers() {
  vwu_progress.swapTarget(&int_data);
  velocity_update_progress.swapTarget(&int_data);
  velocity_constraint_progress.swapTarget(&int_data);
  position_update_progress.swapTarget(&int_data);
  geometry_constraint_progress.swapTarget(&int_data);
  nbwu_progress.swapTarget(&int_data);
  pmewu_progress.swapTarget(&int_data);
  gbrwu_progress.swapTarget(&int_data);
  gbdwu_progress.swapTarget(&int_data);
  gather_wu_progress.swapTarget(&int_data);
  scatter_wu_progress.swapTarget(&int_data);
  all_reduce_wu_progress.swapTarget(&int_data);
}

} // namespace mm
} // namespace stormm
