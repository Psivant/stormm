// -*-c++-*-
#ifndef STORMM_HPC_DEBUG_H
#define STORMM_HPC_DEBUG_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Debug/watcher.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Trajectory/trajectory_enumerators.h"

namespace stormm {
namespace debug {

using card::GpuDetails;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using trajectory::CoordinateCycle;
using trajectory::IntegrationStage;
  
/// \brief Check for large forces in periodic dynamics, whether in the neighbor list or in the
///        fixed-precision coordinate synthesis.
///
/// Overloaded:
///   - Provide a neighbor list object or coordinate synthesis containing coordinates and forces.
///   - Provide the original objects, or abstracts taken at the correct stage of the time cycle
///     if applicable.  This route is taken when engaging the function as an intervention to the
///     main library molecular dynamics protocols.
///   - Indicate whether to take the abstracts at the BLACK or WHITE stage of the coordinate cycle.
///     This route enables developers to check for pollution or initialization failures, but the
///     alternate force arrays are not a typical place to store relevant data.
///
/// \param anom         The event tracking object ("anomalies")
/// \param bugw         Abstract of the event tracking object, with data pointers to memory on the
///                     GPU device ("bug writer", if anomalies indicate bugs)
/// \param poly_ps      The coordinate synthesis
/// \param poly_psr     Abstract of the coordinate synthesis, with data pointers to memory on the
///                     GPU device
/// \param orientation  Indicate whether to take the abstract of poly_ps for the WHITE or BLACK
///                     stage of the coordinate cycle
/// \param step         Step number of the simulation, which will be recorded alongside any
///                     detected events
/// \param when         The point in the MD cycle at which the check is executed, which will be
///                     recorded alongside any detected events
/// \param gpu          Specifications of the GPU that holds the data
/// \{
void checkAtomicForces(WatcherWriter *bugw, const PsSynthesisReader &poly_psr, int step,
                       IntegrationStage when, const GpuDetails &gpu);

void checkAtomicForces(Watcher *anom, const PhaseSpaceSynthesis &poly_ps, int step,
                       IntegrationStage when, const GpuDetails &gpu);

void checkAtomicForces(Watcher *anom, const PhaseSpaceSynthesis &poly_ps,
                       CoordinateCycle orientation, int step, IntegrationStage when,
                       const GpuDetails &gpu);
/// \}

/// \brief Check for high speeds in periodic dynamics, whether in the neighbor list or in the
///        fixed-precision coordinate synthesis.  Descriptions of input variables follow from
///        checkAtomicForces(), above.  There are no overloads involving a neighbor list, as only
///        the coordinate synthesis contains particle velocities.
/// \{
void checkAtomicSpeeds(WatcherWriter *bugw, const PsSynthesisReader &poly_psw, int step,
                       IntegrationStage when, const GpuDetails &gpu);

void checkAtomicSpeeds(Watcher *anom, const PhaseSpaceSynthesis &poly_ps, int step,
                       IntegrationStage when, const GpuDetails &gpu);

void checkAtomicSpeeds(Watcher *anom, const PhaseSpaceSynthesis &poly_ps,
                       CoordinateCycle orientation, int step, IntegrationStage when,
                       const GpuDetails &gpu);
/// \}

/// \brief A general function to gather up debugging activities.  This will reference the
///        global instance of the DynamicsIntervention class (dyna_tk, see dynamics_intervention.h
///        in the src/MolecularMechanics folder) to execute selected checks.
void checkPeriodicDynamics(int step_number, IntegrationStage when);

} // namespace debug
} // namespace stormm

#endif
