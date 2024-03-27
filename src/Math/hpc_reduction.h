// -*-c++-*-
#ifndef STORMM_HPC_REDUCTION_H
#define STORMM_HPC_REDUCTION_H

#include "copyright.h"
#include "Accelerator/core_kernel_manager.h"
#include "MolecularMechanics/mm_controls.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "reduction_abstracts.h"
#include "reduction_bridge.h"
#include "reduction_enumerators.h"

namespace stormm {
namespace stmath {

using card::CoreKlManager;
using constants::PrecisionModel;
using mm::MolecularMechanicsControls;
using mm::MMControlKit;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;

/// \brief Set the __shared__ memory configuration of various reduction kernels to accommodate
///        eight-byte words.
void reductionKernelSetup();
  
/// \brief Obtain the kernel function attributes for one of the reduction kernels.
///
/// \param prec     The precision level at which coordinates and forces are stored--accumulations
///                 by various reduction kernels are always handled in double precision real
///                 numbers, but specifying SINGLE or DOUBLE here indicates whether the kernel
///                 will access overflow arrays for the extended fixed-precision data format.
/// \param purpose  The procedure, e.g. conjugate gradient move vector determination, to be
///                 accomplished by reductions across all systems
/// \param process  Strategy for performing the reduction, one of two stages or all together
cudaFuncAttributes queryReductionKernelRequirements(const PrecisionModel prec,
                                                    const ReductionGoal purpose,
                                                    const ReductionStage process);

/// \brief Launch the appropriate kernel to do a conjugate transformation of the computed forces
///        in preparation for one line step of conjugate gradient energy minimization.
///
/// \param prec      Precision model used to handle the arithmetic for all interactions.  Also
///                  informs whether to make use of extended fixed-precision representations in
///                  the coordinate synthesis object.
/// \param poly_ag   Synthesis of all systems' topologies
/// \param poly_ps   Synthesis of all systems' coordinates, forces, and other data arrays
/// \param rbg       Pre-allocated scratch space used by multi-staged reduction work units
/// \param mmctrl    Molecular mechanics control data and progress counters
/// \param launcher  Holds launch parameters for all kernels, and can be queried for launch
///                  parameters of the conjugate gradient transformation at hand
/// \param redk      Abstract for the reduction work units extracted from the topology synthesis
/// \param cgsbs     Recasts some arrays of the coordinate synthesis and ties them together with
///                  allocations in the reduction bridge
/// \param ctrl      Molecular mechanics control abstract, with pointers to progress counter arrays
/// \param bt        Launch parameters for the conjugate gradient transformation kernel
/// \{
void launchConjugateGradient(const ReductionKit &redk, ConjGradSubstrate *cgsbs,
                             MMControlKit<double> *ctrl, const int2 bt);

void launchConjugateGradient(const ReductionKit &redk, ConjGradSubstrate *cgsbs,
                             MMControlKit<float> *ctrl, const int2 bt);

void launchConjugateGradient(PrecisionModel prec, const AtomGraphSynthesis poly_ag,
                             PhaseSpaceSynthesis *poly_ps, ReductionBridge *rbg,
                             MolecularMechanicsControls *mmctrl, const CoreKlManager &launcher);
/// \}

} // namespace synthesis
} // namespace stormm

#endif
