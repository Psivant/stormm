// -*-c++-*-
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/symbol_values.h"
#include "Potential/energy_enumerators.h"
#include "hpc_kinetic.h"

namespace stormm {
namespace mm {

using card::Hybrid;
using card::HybridTargetLevel;
using energy::StateVariable;
using symbols::boltzmann_constant;

/// \brief Simple kernel for translating total kinetic energies in each system into temperatures.
///
/// \param dof  The number of degrees of freedom in each system, taken with consideration to
///             whether constraints apply
/// \param scw  Writeable abstract of the energy tracking object.  Kinetic energies are expected to
///             already have been computed.
__global__ __launch_bounds__(large_block_size, 1)
void kComputeTemperature(const int* dof, ScoreCardWriter scw) {
  int pos = threadIdx.x + (blockIdx.x * blockDim.x);
  while (pos < scw.system_count) {
    const size_t read_idx = (scw.data_stride * pos) + (int)(StateVariable::KINETIC);
    const size_t write_idx = (scw.data_stride * pos) + (int)(StateVariable::TEMPERATURE_ALL);
    const llint ike = scw.instantaneous_accumulators[read_idx];
    const double tmpr = 2.0 * (double)(ike) / ((double)(dof[pos]) * boltzmann_constant);
    scw.instantaneous_accumulators[write_idx] = __double2ll_rn(tmpr);
    pos += (blockDim.x * gridDim.x);
  }
}

//-------------------------------------------------------------------------------------------------
void launchTemperatureComputation(const AtomGraphSynthesis &poly_ag, ScoreCardWriter *scw,
                                  const bool use_constraints, const GpuDetails &gpu) {
  const Hybrid<int>& system_dof = (use_constraints) ? poly_ag.getConstrainedDegreesOfFreedom() :
                                                      poly_ag.getDegreesOfFreedom();
  const int* dof_ptr = system_dof.data(HybridTargetLevel::DEVICE);
  kComputeTemperature<<<gpu.getSMPCount(), large_block_size>>>(dof_ptr, *scw);
}

} // namespace mm
} // namespace stormm
