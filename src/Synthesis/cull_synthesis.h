// -*-c++-*-
#ifndef STORMM_CULL_SYNTHESIS_H
#define STORMM_CULL_SYNTHESIS_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "atomgraph_synthesis.h"
#include "phasespace_synthesis.h"
#include "systemcache.h"
#include "synthesis_cache_map.h"

namespace stormm {
namespace synthesis {

using card::GpuDetails;
using synthesis::SynthesisCacheMap;
  
/// \brief Cull a topology synthesis to serve a coordinate synthesis based on a list of individual
///        systems which have met certain criteria, whether sanity checks or a selected subset.
///        The function will create and return a new topology synthesis.
///
/// Overloaded:
///   - Provide the topology synthesis by const pointer
///   - Provide the topology synthesis by const reference
///  
/// \param original   The original topology synthesis.  The original object will not be destroyed
///                   or modified, but could be overwritten by the result.
/// \param survivors  A list of system indices which will be carried over into the result.  This
///                   vector may contain duplicate entries, and even make the "culled" synthesis
///                   larger overall than the original.
/// \param gpu        GPU that will be utilized in calculations involving the topology synthesis
/// \{
AtomGraphSynthesis cullSynthesis(const AtomGraphSynthesis *original,
                                 const std::vector<int> &survivors,
                                 const GpuDetails &gpu = null_gpu,
                                 ExceptionResponse policy = ExceptionResponse::WARN);

AtomGraphSynthesis cullSynthesis(const AtomGraphSynthesis &original,
                                 const std::vector<int> &survivors,
                                 const GpuDetails &gpu = null_gpu,
                                 ExceptionResponse policy = ExceptionResponse::WARN);
/// \}

/// \brief Cull a coordinate synthesis based on a mask of system indices.  There is no overloaded
///        variant for a Condensate object, as a new Condensate can be made and tied to the culled
///        PhaseSpaceSynthesis.
///
/// Overloaded:
///   - Provide the coordinate synthesis by const pointer
///   - Provide the coordinate synthesis by const reference
///
/// \param original   The original coordinate synthesis.  The original object will not be destroyed
///                   or modified, but could be overwritten by the result.
/// \param survivors  A list of system indices which will be carried over into the result.  This
///                   vector may contain duplicate entries, and even make the "culled" synthesis
///                   larger overall than the original.
/// \{
PhaseSpaceSynthesis cullSynthesis(const PhaseSpaceSynthesis *original,
                                  const std::vector<int> &survivors, const GpuDetails &gpu);
  
PhaseSpaceSynthesis cullSynthesis(const PhaseSpaceSynthesis &original,
                                  const std::vector<int> &survivors, const GpuDetails &gpu);
/// \}

/// \brief Build a new cache map to reflect how an updated synthesis corresponds to structures in
///        the underlying cache.
///
/// \param original   The original cache map.  This will contain a pointer ot the underlying system
///                   cache obtained from user input.
/// \param sruvivors  A list of system indices which will be carried over into the result.  This
///                   will be compared to the mapping present in the original object in order to
///                   derive the new mapping.
SynthesisCacheMap remapSynthesis(const SynthesisCacheMap &original,
                                 const std::vector<int> &survivors);

} // namespace synthesis
} // namespace stormm

#endif
