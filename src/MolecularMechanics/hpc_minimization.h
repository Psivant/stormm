// -*-c++-*-
#ifndef STORMM_HPC_MINIMIZATION_H
#define STORMM_HPC_MINIMIZATION_H

#include "copyright.h"
#include "Accelerator/core_kernel_manager.h"
#include "Accelerator/gpu_details.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "MolecularMechanics/line_minimization.h"
#include "MolecularMechanics/mm_controls.h"
#include "Namelists/nml_minimize.h"
#include "Numerics/split_fixed_precision.h"
#include "Potential/cacheresource.h"
#include "Potential/scorecard.h"
#include "Math/reduction_bridge.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/implicit_solvent_workspace.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/static_mask_synthesis.h"
#include "UnitTesting/stopwatch.h"

namespace stormm {
namespace mm {

using card::GpuDetails;
using card::CoreKlManager;
using constants::PrecisionModel;
using energy::CacheResource;
using energy::ScoreCard;
using energy::ScoreCardWriter;
using stmath::ReductionBridge;
using stmath::ReductionKit;
using namelist::MinimizeControls;
using numerics::default_energy_scale_bits;
using numerics::AccumulationMethod;
using synthesis::AtomGraphSynthesis;
using synthesis::ImplicitSolventWorkspace;
using synthesis::ISWorkspaceKit;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using synthesis::StaticExclusionMaskSynthesis;
using testing::StopWatch;
  
/// \brief Set the __shared__ memory access size for the conjugate gradient particle advancement
///        kernels.
void minimizationKernelSetup();

/// \brief Get the launch parameters for conjugate gradient particle advancement kernels.
///
/// \param prec  The precision model of the coordinate object (the particle advancement is always
///              performed in double-precision real numbers)
cudaFuncAttributes queryMinimizationKernelRequirements(const PrecisionModel prec);

/// \brief Launch the line minimization kernel in a standalone capacity.  This function expedites
///        debugging and provides abstraction of the line movement function for other applications.
///
/// \param prec         The precision model in use by the coordinate set (indicates whether the
///                     extended fixed-precision number format is in use)
/// \param poly_psw     Coordinates for each system and forces acting on each system
/// \param redk         Reduction work units to inform the limits of each thread block's activity
/// \param scw          Read-only access to the current energy of each system
/// \param lmw          Writeable abstract for a line minimization manager
/// \param move_number  The number of the move--values of 0 to 3 will lead to different behavior
/// \param redu_lp      Launch parameters for the kernel (an abstract of the CoreKlManager object)
/// \{
void launchLineAdvance(PrecisionModel prec, PsSynthesisWriter *poly_psw, const ReductionKit &redk,
                       const ScoreCardWriter &scw, LinMinWriter *lmw, int move_number,
                       const int2 redu_lp);

void launchLineAdvance(PrecisionModel prec, PhaseSpaceSynthesis *poly_ps,
                       const AtomGraphSynthesis &poly_ag, ScoreCard *sc,
                       LineMinimization *line_record, int move_number,
                       const CoreKlManager &launcher);
/// \}

/// \brief Run energy minimizations of all structures in a synthesis based on user input or other
///        data contained in a molecular mechanics control object.
///
/// Overloaded:
///   - Accept minimal arguments and refine the objects internally to construct all necessary
///     abstracts and scratch spaces for the minimization (this will return the energy series
///     generated by the energy minimizations)
///   - Accept an intermediate level of prepared objects (e.g. the non-bonded exclusions, in
///     addition to basic coordinates, topologies, and run controls), then prepare the rest
///     internally (this will also return the energy series generated by the minimizations)
///   - Accept all abstracts and work spaces to run minimization with no temporary allocations
///     (the energy minimizations will accumulate the results in the prepared energy tracking
///     object and leave it modified when the function completes)
///
/// \param prec             The precision to use in general arithmetic.  Also activates an extended
///                         mode of fixed-precision data representation in the coordinates object.
/// \param poly_ag          Compilation of all systems' topologies (a "synthesis" of systems)
/// \param poly_se          Synthesis of exclusion masks, tailored to the topology synthesis
/// \param poly_ps          Compilation of coordinates, forces, and other useful arrays for
///                         minimizing each structure
/// \param mmctrl_fe        Control data for the minimization run, plus progress counters for GPU
///                         kernel management in force + energy computations
/// \param mmctrl_xe        Control data for the minimization run, plus progress counters for GPU
///                         kernel management in energy-only computations
/// \param mmctrl_cdfe      Control data for the minimization run, plus progress counters for GPU
///                         kernel management in force + energy + clash dampiing computations
/// \param mmctrl_cdxe      Control data for the minimization run, plus progress counters for GPU
///                         kernel management in energy + clash damping computations
/// \param sc               Energy tracking object allocated to cover all systems in the synthesis
/// \param vale_fe_cache    Pre-allocated scratch space for valence kernel thread blocks in force +
///                         energy computations
/// \param vale_xe_cache    Pre-allocated scratch space for valence kernel thread blocks in force +
///                         energy computations
/// \param vale_cdfe_cache  Pre-allocated scratch space for valence kernel thread blocks in force +
///                         energy computations with clash damping
/// \param vale_cdxe_cache  Pre-allocated scratch space for valence kernel thread blocks in force +
///                         energy computations with clash damping
/// \param nonb_cache       Pre-allocated scratch space for non-bonded kernel thread blocks
/// \param nonb_cd_cache    Pre-allocated scratch space for non-bonded kernel thread blocks with
///                         clash damping
/// \param ism_space        Pre-allocated space to store intermediate computations for Born radii
///                         and their derivatives
/// \param rbg              Pre-allocated storage space for reduction work units
/// \param line_record      Stores line minimization progress, energies and move lengths
/// \param acc_meth         Choice of force accumulation method
/// \param gpu              Details of the GPU in use
/// \param launcher         Holds parameters for all kernels, will be accessed for the correct
///                         launch configurations of the kernels relevant to the workload at hand
/// \{
void launchMinimization(PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
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
                        LineMinimization *line_record, AccumulationMethod acc_meth,
                        const GpuDetails &gpu, const CoreKlManager &launcher,
                        StopWatch *timer = nullptr,
                        const std::string &task_name = std::string(""));

ScoreCard launchMinimization(const AtomGraphSynthesis &poly_ag,
                             const StaticExclusionMaskSynthesis &poly_se,
                             PhaseSpaceSynthesis *poly_ps, const MinimizeControls &mincon,
                             const GpuDetails &gpu, PrecisionModel prec = PrecisionModel::SINGLE,
                             int energy_accumulation_bits = default_energy_scale_bits, 
                             StopWatch *timer = nullptr,
                             const std::string &task_name = std::string(""));

ScoreCard launchMinimization(const AtomGraphSynthesis &poly_ag, PhaseSpaceSynthesis *poly_ps,
                             const MinimizeControls &mincon, const GpuDetails &gpu,
                             PrecisionModel prec = PrecisionModel::SINGLE,
                             int energy_accumulation_bits = default_energy_scale_bits, 
                             StopWatch *timer = nullptr,
                             const std::string &task_name = std::string(""));
/// \}
  
} // namespace mm
} // namespace stormm

#endif
