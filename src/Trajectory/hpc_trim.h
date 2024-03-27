// -*-c++-*-
#ifndef STORMM_HPC_TRIM_H
#define STORMM_HPC_TRIM_H

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Synthesis/phasespace_synthesis.h"
#include "motion_sweeper.h"

namespace stormm {
namespace trajectory {

using card::GpuDetails;
using synthesis::PsSynthesisReader;
using synthesis::PsSynthesisWriter;
using synthesis::SyAtomUpdateKit;

/// \brief Launch the GPU kernel to accumulate motion of the center of mass for all systems.
///
/// \param mosw      Contains accumulators for center of mass, momentum, and moment of inertia
/// \param poly_auk  Contains masses of all particles in all systems
/// \param poly_psr  Contains coordinates and velocities of all particles in all systems
/// \param gpu       Specifications of the GPU that will perform the calculations.
void launchAccCenterOfMassMotion(MotionSweepWriter *mosw,
                                 const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                                 const PsSynthesisReader &poly_psr, const GpuDetails &gpu);

/// \brief Launch the kernel to recenter and remove center of mass motion in all systems.
///
/// \param poly_psw  Contains coordinates and velocities of all particles in all systems
/// \param mosr      Contains accumulators for center of mass, momentum, and moment of inertia
/// \param gpu       Specifications of the GPU that will perform the calculations.
void launchRemoveCenterOfMassMotion(PsSynthesisWriter *poly_psw, const MotionSweepReader &mosr,
                                    const GpuDetails &gpu);

/// \brief Launch the kernel to accumulate angular momentum and the inertial tensor.  While only
///        valid for non-periodic systems with no boundary conditions, this function relies on the
///        calling function to evaluate the unit cell type of the coordinate synthesis.
///
/// \param mosw      Contains accumulators for center of mass, momentum, and moment of inertia
/// \param poly_auk  Contains masses of all particles in all systems
/// \param poly_psr  Contains coordinates and velocities of all particles in all systems
void launchAccAngularMomentum(MotionSweepWriter *mosw,
                              const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                              const PsSynthesisReader &poly_psr, const GpuDetails &gpu);

/// \brief Launch the kernel to remove the net rotational velocity from each system.
///
/// \param poly_psw  Contains coordinates and velocities of all particles in all systems
/// \param mosr      Contains accumulators for center of mass, momentum, and moment of inertia
/// \param gpu       Specifications of the GPU that will perform the calculations.
void launchRemoveAngularMomentum(PsSynthesisWriter *poly_psw, const MotionSweepReader &mosr,
                                 const GpuDetails &gpu);
  
} // namespace trajectory
} // namespace stormm

#endif
