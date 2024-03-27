// -*-c++-*-
#ifndef STORMM_HPC_INTEGRATION_H
#define STORMM_HPC_INTEGRATION_H

#ifdef STORMM_USE_CUDA
#  include <cuda_runtime.h>
#endif
#include "copyright.h"
#include "Accelerator/core_kernel_manager.h"
#include "Constants/behavior.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "MolecularMechanics/mm_controls.h"
#include "Numerics/numeric_enumerators.h"
#include "Potential/cacheresource.h"
#include "Potential/energy_enumerators.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "thermostat.h"
#include "trajectory_enumerators.h"

namespace stormm {
namespace trajectory {

using card::CoreKlManager;
using constants::PrecisionModel;
using energy::CacheResource;
using energy::CacheResourceKit;
using energy::ValenceKernelSize;
using mm::MolecularMechanicsControls;
using mm::MMControlKit;
using numerics::AccumulationMethod;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using synthesis::SyAtomUpdateKit;
using synthesis::SyValenceKit;
  
/// \brief Query the launch parameters and other attributes of one of the stand-alone integration
///        kernels.
///
/// \param prec      The precision with which arithmetic operations will be carried out.  This also
///                  affects the extent of coordinate representations processed by the kernel.
/// \param acc_meth  The accumulation method for forces in related valence kernels.  Each kernel's
///                  core lines of code are designed to be absracted as an included file inside the
///                  main valence kernels with metching kernel sizes, accumulation methods, and
///                  calculation precision.  The stand-alone kernels may utilize the force arrays
///                  of the valence kernel to buffer critical data, and must therefore emulate the
///                  methods of the valence kernels.
/// \param kwidth    Conveys the size of each work unit that the kernel must be prepared to handle
/// \param process   The process to carrry out, e.g. constrain velocities
cudaFuncAttributes queryIntegrationKernelRequirements(PrecisionModel prec,
                                                      AccumulationMethod acc_meth,
                                                      ValenceKernelSize kwidth,
                                                      IntegrationStage process);

/// \brief Launch a stand-alone kernel to carry out one of several processes in the integration
///        time step.
///
/// Overloaded:
///   - Provide the original coordinate synthesis, topology synthesis, thermostat and other objects
///   - Provide abstracts for all objects
///
/// \param poly_psw   Abstract of the coordinate synthesis, taken at the correct stage of the
///                   time cycle
/// \param poly_ps    Coordinates, velocities, and forces for all systems
/// \param tb_resw    Abstract of the thread block space, prepared based on the precision model of
///                   calculations to be carried out
/// \param tb_res     Private L1 cache resources allocated for all thread blocks
/// \param tstw       Thermostat abstract containing random numbers for the thermostat and also
///                   information about constraint convergence and the time step.  This could be
///                   read-only, but the need to alter data within the thermostat duing the force
///                   calculation, to 
/// \param tst        Thermostat abstract 
/// \param poly_vk    Valence abstract for all systems, including valence work units to guide
///                   constraint work and delegate the movement of atoms amongst different thread
///                   blocks
/// \param poly_auk   Atom update information, including masses and geometric constraint parameters
/// \param poly_ag    Topologies for all systems, able to dispense abstracts for valence work
///                   units, constraints, and atom updates
/// \param lp         Launch parameters for the appropriate GPU kernel (in essence, an abstract of
///                   the launcher object)
/// \param launcher   Holds launch parameters for all major kernels involved in MD and molecular
///                   energy minimization
/// \param prec       The precision with which arithmetic operations will be carried out.
/// \param process    The part of the integration workflow to implement
/// \{
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *tb_resw,
                              MMControlKit<double> *ctrl, const SyValenceKit<double> &poly_vk,
                              const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                              const ThermostatWriter<double> &tstw, const int2 lp,
                              IntegrationStage process);

void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *tb_resw,
                              MMControlKit<float> *ctrl, const SyValenceKit<float> &poly_vk,
                              const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                              const ThermostatWriter<float> &tstw, const int2 lp,
                              AccumulationMethod acc_meth, ValenceKernelSize kwidth,
                              IntegrationStage process);

void launchIntegrationProcess(PhaseSpaceSynthesis *poly_ps, CacheResource *tb_res, Thermostat *tst,
                              MolecularMechanicsControls *mmctrl,
                              const AtomGraphSynthesis &poly_ag, const CoreKlManager &launcher,
                              PrecisionModel prec, AccumulationMethod acc_meth,
                              IntegrationStage process);
/// \}

} // namespace trajectory
} // namespace stormm

#endif
