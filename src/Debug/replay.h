// -*-c++-*-
#ifndef STORMM_REPLAY_H
#define STORMM_REPLAY_H

#include <vector>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Synthesis/phasespace_synthesis.h"

namespace stormm {
namespace debug {

using card::GpuDetails;
using card::Hybrid;
using synthesis::PhaseSpaceSynthesis;
  
/// \brief Store a series of frames for each system in a synthesis, convering positions,
///        velocities, and forces using a collection of PhaseSpaceSynthesis objects.
class Replay {
public:

  /// \brief Constructors require an existing coordinate synthesis and a depth to which history is
  ///        to be recorded.
  ///
  /// \param current_step  The step of the simulation at which the object is first created
  Replay(const PhaseSpaceSynthesis *poly_ps_in, int depth_in, int current_step = 0);
  
  /// \brief Get the depth of the history contained in the object.
  int getDepth() const;

  /// \brief Get the number of unique systems contained in the synthesis served by the Replay
  ///        object.
  int getSystemCount() const;

  /// \brief Get a pointer to the original coordinate synthesis, which will contain the current
  ///        state of all systems in the calculation.
  const PhaseSpaceSynthesis* getCoordinateSynthesisPointer() const;

  /// \brief Get a pointer to the snapshot of the synthesis recorded for a particular step number,
  ///        relative to the current step number, if such a record is yet (or still) available.
  ///        From the standpoint of the object, the current step number is calculated according
  ///        to the original step index fed to the constructor plus any updates that occur over
  ///        the object's lifetime.
  ///
  /// \param relative_index  The relative index of the step number of interest, in reference to the
  ///                        current state of the synthesis of systems.  This number is expected to
  ///                        be less than or equal to zero, and not larger than the depth of the
  ///                        object.
  const PhaseSpaceSynthesis* getSnapshot(int relative_index = 0) const;

  /// \brief Copy the current state of the synthesis into the object's next available holdings
  ///        slot.
  ///
  /// \param gpu  Details of the GPU performing the calculations (if present, which will also be
  ///             where the device memory is stored)
  void takeSnapshot(const GpuDetails &gpu);
  
private:

  // A handful of integers controls the traffic and contents of the object.
  int depth;                ///< The total number of snapshots of the coordinate synthesis to
                            ///<   allocate in the holdings array
  int next_holdings_index;  ///< Index of the next holdings array element to copy into

  /// The inverse of steps_recorded, a series of indices into the holdings array tracking which
  /// sequence of elements one should access in order to read the state of the coordinate synthesis
  /// at step k, k - 1, k - 2, ...
  std::vector<int> backtrace_order;
  
  /// A series of coordinate syntheses, equivalent to the first, which will store the history of
  /// all particle movements for the preceding few steps.  Pointers to objects in each element of
  /// this array will be incremented so as to keep a continuous, up-to-date record of each move
  /// at any depth without copying more than the contents of one coordinate synthesis per step.
  std::vector<PhaseSpaceSynthesis> holdings;

  /// A list of all systems in the synthesis, prepared as a correspondence map to feed to one of
  /// the coordinate copying overloads.
  Hybrid<int2> system_pairs;
  
  /// A pointer to the original coordinate synthesis.
  PhaseSpaceSynthesis *poly_ps_ptr;

  /// Increment the value of the next holdings index.  The index will tick up so long as there is
  /// space in the object, but upon reaching the available depth will cycle back to zero.
  void incrementNextHoldingsIndex();
};

} // namespace debug
} // namespace stormm

#endif
