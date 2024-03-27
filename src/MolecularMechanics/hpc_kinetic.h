// -*-c++-*-
#ifndef STORMM_HPC_KINETIC_H
#define STORMM_HPC_KINETIC_H

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Potential/scorecard.h"
#include "Trajectory/thermostat.h"
#include "Synthesis/atomgraph_synthesis.h"

namespace stormm {
namespace mm {

using card::GpuDetails;
using energy::ScoreCardWriter;
using synthesis::AtomGraphSynthesis;

/// \brief Evaluate the temperautre for a synthesis of systems.  It is expected that kinetic
///        energies for all systems have already been computed.
///
/// \param poly_ag          The synthesis of topologies.  No abstract will be taken from this,
///                         rather a sole pointer to the appropriate array of degrees of freedom
///                         will be extracted.  One check to ensure that the system counts are
///                         consistent between the synthesis and the energy tracking object will
///                         be performed.
/// \param scw              Writeable abstract of the energy tracking object, taken with pointers
///                         set to memory on the GPU device
/// \param use_constraints  Indicate whether constraints are in effect, which has bearing on the
///                         number of degrees of freedom to consider
/// \param gpu              Details of the GPU that will carry out the very basic calculation.  The
///                         kernel is so simple and short that no optimization or kernel manager is
///                         needed.
void launchTemperatureComputation(const AtomGraphSynthesis &poly_ag, ScoreCardWriter *scw,
                                  const bool use_constraints, const GpuDetails &gpu);

} // namespace mm
} // namespace stormm

#endif
