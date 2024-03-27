// -*-c++-*-
#ifndef STORMM_HPC_RMSD_H
#define STORMM_HPC_RMSD_H

#ifdef STORMM_USE_CUDA
#include <cuda_runtime.h>
#endif
#include "copyright.h"
#include "Accelerator/core_kernel_manager.h"
#include "Accelerator/hybrid.h"
#include "Analysis/comparison_guide.h"
#include "Constants/behavior.h"
#include "Synthesis/phasespace_synthesis.h"
#include "rmsd_plan.h"

namespace stormm {
namespace structure {

using analysis::ComparisonGuide;
using card::Hybrid;
using card::CoreKlManager;
using constants::PrecisionModel;
using synthesis::PhaseSpaceSynthesis;

/// \brief Obtain the kernel requirements for one of the RMSD calculation kernels.
///
/// \param prec   The type of floating point numbers in which the kernel shall work
/// \param order  The order of the calculation: it can be RMSD to reference calculations, or
///               matrix RMSD calculations spanning all pairs of systems with the same topology
cudaFuncAttributes queryRMSDKernelRequirements(PrecisionModel prec, RMSDTask order);

/// \brief Compute positional RMSDs using HPC resources.
///
/// \param cg                Instructions for carrying out structure-to-structure comparisons and
///                          computing necessary RMSD results
/// \param rplan             Plan appropriate to the coordinate synthesis encompassing system
///                          sizes, counts, atomic masses, and possible structural symmetries
/// \param poly_ps           Coordinates of many systems with various numbers of replicas, and maps
///                          illustrating which systems share the same topologies.  Replicas of
///                          systems sharing the same topologies will be compared in terms of RMSD.
/// \param reference_frames  Indices within the list of systems kept by poly_ps for reference
///                          systems to be used in "all to one" RMSD calculations spanning systems
///                          using each unique topology.
/// \param result            Pre-allocated results array with contiguous space for systems making
///                          use of each topology.
/// \param launcher          Source of kernel launch parameters
/// \{
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &poly_ps,
          const Hybrid<int> &reference_frames, Hybrid<double> *result,
          const CoreKlManager &launcher);

void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &poly_ps,
          const Hybrid<int> &reference_frames, Hybrid<float> *result,
          const CoreKlManager &launcher);

void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &poly_ps,
          Hybrid<double> *result, const CoreKlManager &launcher);

void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &poly_ps,
          Hybrid<float> *result, const CoreKlManager &launcher);
/// \}

} // namespace structure
} // namespace stormm

#endif
