// -*-c++-*-
#ifndef STORMM_VALENCE_POTENTIAL_CUH
#define STORMM_VALENCE_POTENTIAL_CUH

#ifdef STORMM_USE_CUDA
#  include <cuda_runtime.h>
#endif
#include "copyright.h"
#include "Accelerator/core_kernel_manager.h"
#include "Accelerator/gpu_details.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/stormm_vector_types.h"
#include "MolecularMechanics/mm_controls.h"
#include "Potential/energy_enumerators.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Trajectory/thermostat.h"
#include "cacheresource.h"
#include "energy_enumerators.h"
#include "scorecard.h"

namespace stormm {
namespace energy {

using card::GpuDetails;
using card::CoreKlManager;
using constants::PrecisionModel;
using mm::MMControlKit;
using mm::MolecularMechanicsControls;
using numerics::AccumulationMethod;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using synthesis::SyAtomUpdateKit;
using synthesis::SyRestraintKit;
using synthesis::SyValenceKit;
using synthesis::VwuGoal;
using trajectory::Thermostat;
using trajectory::ThermostatWriter;

/// \brief Test various subdivisions of a thread block to see whether the workload can be better
///        distributed.  This is perilous if the kernel allocates significant static __shared__
///        memory resources per block, but can be advantageous if the systems are small and there
///        are many to compute.  This assumes that the kernel is compiled with a minimum of one
///        block per multiprocessor in its launch bounds, and the highest possible thread count
///        for that one block.  The function returns a tuple containing the best identified
///        block multiplier and the best corresponding block size, i.e. "split a kernel of 896
///        maximum threads into 7 blocks of 128."
///
/// \param max_threads  The maximum number of threads for which the kernel is compiled
/// \param smp_count    The number of streaming multiprocessors on the device
/// \param vwu_size     Maximum size of each valence work unit
/// \param vwu_count    The total number of valence work units that the kernel will process
int2 testValenceKernelSubdivision(const int max_threads, const int smp_count, const int vwu_size,
                                  const int vwu_count);

/// \brief Obtain information on launch bounds and block-specific requirements for a selected
///        version of the valence interactions kernel.
///
/// \param prec                Precision model for the selected kernel
/// \param eval_frc            Select whether to evaluate forces
/// \param eval_nrg            Select whether to evaluate energies
/// \param acc_meth            Accumulation method for forces
/// \param purpose             Indicate whether the kernel shall deposit its accumulated forces
///                            back into global arrays or use them immediately to move particles
/// \param collision_handling  Indication of whether to forgive clashes
/// \param kwidth              Indication of the desired thread block size in the kernel
cudaFuncAttributes queryValenceKernelRequirements(PrecisionModel prec, EvaluateForce eval_frc,
                                                  EvaluateEnergy eval_nrg,
                                                  AccumulationMethod acc_meth, VwuGoal purpose,
                                                  ClashResponse collision_handling,
                                                  ValenceKernelSize kwidth);

/// \brief Evaluate valence work units and move atoms.
///
/// Overloaded:
///   - Perform work in single or double precision
///   - Compute forces, energy, or both
///   - Move particles, or instead accumulate forces, energies, or both
///
/// \param poly_vk      Valence parameters based on consensus tables from a topology synthesis
/// \param poly_auk     Abstract of the topology critical for atom movement
/// \param poly_ag      Compiled topologies of all systems
/// \param poly_psw     Abstract for coordinates and forces of all systems
/// \param poly_ps      Coordinates, velocities, and forces of all systems
/// \param tstw         Abstract of the thermostat
/// \param heat_bath    Thermostat for keeping each system at a particular temperature
/// \param gmem_r       Exclusive space in global memory arrays reserved for each thread block, to
/// \param scw          Abstract for the energy tracking on all systems
/// \param sc           Energy tracking for all systems (this must be pre-sized to accommodate the
///                     entire synthesis)
/// \param ctrl         Abstract for molecular mechanics progress counters and run bounds
/// \param mmctrl       Progress counters and run bounds
///                     be brought into free L1 cache
/// \param tb_space     Cache resources for the kernel launch
/// \param eval_force   Whether to have forces evaluated
/// \param eval_energy  Whether to have energy components evaluated
/// \param purpose      Whether to move particles or accumulate forces and energies
/// \param bt           Block and thread counts for the kernel launch
/// \param launcher     General repository for launch parameters of all kernels
/// \{
void launchValence(const SyValenceKit<double> &poly_vk,
                   const SyRestraintKit<double, double2, double4> &poly_rk,
                   MMControlKit<double> *ctrl, PsSynthesisWriter *poly_psw,
                   const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                   ThermostatWriter<double> *tstw, ScoreCardWriter *scw,
                   CacheResourceKit<double> *gmem_r, EvaluateForce eval_force,
                   EvaluateEnergy eval_energy, VwuGoal purpose, const int2 bt,
                   double clash_minimum_distance = 0.0, double clash_ratio = 0.0);

void launchValence(const SyValenceKit<float> &poly_vk,
                   const SyRestraintKit<float, float2, float4> &poly_rk,
                   MMControlKit<float> *ctrl, PsSynthesisWriter *poly_psw,
                   const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                   ThermostatWriter<float> *tstw, ScoreCardWriter *scw,
                   CacheResourceKit<float> *gmem_r, EvaluateForce eval_force,
                   EvaluateEnergy eval_energy, VwuGoal purpose, AccumulationMethod force_sum,
                   const int2 bt, float clash_minimum_distance = 0.0, float clash_ratio = 0.0);
  
void launchValence(PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                   MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                   Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                   EvaluateForce eval_force, EvaluateEnergy eval_energy, VwuGoal purpose,
                   AccumulationMethod force_sum, const CoreKlManager &launcher,
                   double clash_minimum_distance = 0.0, double clash_ratio = 0.0);

void launchValence(PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                   MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                   Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                   EvaluateForce eval_force, EvaluateEnergy eval_energy, VwuGoal purpose,
                   const CoreKlManager &launcher, double clash_minimum_distance = 0.0,
                   double clash_ratio = 0.0);
/// \}
  
} // namespace energy
} // namespace stormm

#endif
