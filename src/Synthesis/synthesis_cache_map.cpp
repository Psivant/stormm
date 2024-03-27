#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "Math/series_ops.h"
#include "synthesis_cache_map.h"

namespace stormm {
namespace synthesis {

using card::HybridKind;
using stmath::indexingArray;

//-------------------------------------------------------------------------------------------------
SynthesisMapReader::SynthesisMapReader(const int ncache_in, const int nsynth_in,
                                       const int nlabel_in, const int ntopol_in,
                                       const int* cache_origins_in, const int* topology_origins_in,
                                       const int* label_origins_in, const int* csystem_proj_in,
                                       const int* csystem_bounds_in, const int* clabel_proj_in,
                                       const int* clabel_bounds_in, const int* ctopol_proj_in,
                                       const int* ctopol_bounds_in) :
    ncache{ncache_in}, nsynth{nsynth_in}, nlabel{nlabel_in}, ntopol{ntopol_in},
    cache_origins{cache_origins_in}, topology_origins{topology_origins_in},
    label_origins{label_origins_in}, csystem_proj{csystem_proj_in},
    csystem_bounds{csystem_bounds_in}, clabel_proj{clabel_proj_in},
    clabel_bounds{clabel_bounds_in}, ctopol_proj{ctopol_proj_in}, ctopol_bounds{ctopol_bounds_in}
{}

//-------------------------------------------------------------------------------------------------
SynthesisCacheMap::SynthesisCacheMap() :
    cache_system_count{0}, synthesis_system_count{0}, cache_label_count{0},
    cache_topology_count{0},
    cache_origins{HybridKind::POINTER, "syncache_orig"},
    topology_origins{HybridKind::POINTER, "syncache_orig"},
    label_origins{HybridKind::POINTER, "syncache_orig"},
    sys_projections{HybridKind::POINTER, "syncache_sysproj"},
    sys_projection_bounds{HybridKind::POINTER, "syncache_sp_bounds"},
    label_projections{HybridKind::POINTER, "syncache_labelproj"},
    label_projection_bounds{HybridKind::POINTER, "syncache_lbl_bounds"},
    topology_projections{HybridKind::POINTER, "syncache_agproj"},
    topology_projection_bounds{HybridKind::POINTER, "syncache_ag_bounds"},
    int_data{HybridKind::ARRAY, "syncache_data"},
    sc_ptr{nullptr}, poly_ag_ptr{nullptr}, poly_ps_ptr{nullptr}
{}
    
//-------------------------------------------------------------------------------------------------
SynthesisCacheMap::SynthesisCacheMap(const std::vector<int> &cache_origins_in,
                                     const SystemCache *sc_in,
                                     const AtomGraphSynthesis *poly_ag_in,
                                     const PhaseSpaceSynthesis *poly_ps_in) :
    SynthesisCacheMap()
{
  setCache(cache_origins_in, sc_in);
  if (poly_ag_in != nullptr) {
    setSynthesis(poly_ag_in);
  }
  if (poly_ps_in != nullptr) {
    setSynthesis(poly_ps_in);
  }
  validateCorrespondence();
}

//-------------------------------------------------------------------------------------------------
SynthesisCacheMap::SynthesisCacheMap(const std::vector<int> &cache_origins_in,
                                     const SystemCache &sc_in,
                                     const AtomGraphSynthesis &poly_ag_in,
                                     const PhaseSpaceSynthesis &poly_ps_in) :
    SynthesisCacheMap(cache_origins_in, sc_in.getSelfPointer(), poly_ag_in.getSelfPointer(),
                      poly_ps_in.getSelfPointer())
{}

//-------------------------------------------------------------------------------------------------
SynthesisCacheMap::SynthesisCacheMap(const std::vector<int> &cache_origins_in,
                                     const SystemCache &sc_in,
                                     const PhaseSpaceSynthesis &poly_ps_in) :
    SynthesisCacheMap(cache_origins_in, sc_in.getSelfPointer(), nullptr,
                      poly_ps_in.getSelfPointer())
{}

//-------------------------------------------------------------------------------------------------
SynthesisCacheMap::SynthesisCacheMap(const std::vector<int> &cache_origins_in,
                                     const SystemCache &sc_in,
                                     const AtomGraphSynthesis &poly_ag_in) :
    SynthesisCacheMap(cache_origins_in, sc_in.getSelfPointer(), poly_ag_in.getSelfPointer(),
                      nullptr)
{}

//-------------------------------------------------------------------------------------------------
SynthesisCacheMap::SynthesisCacheMap(const SynthesisCacheMap &original) :
    cache_system_count{original.cache_system_count},
    synthesis_system_count{original.synthesis_system_count},
    cache_label_count{original.cache_label_count},
    cache_topology_count{original.cache_topology_count},
    cache_origins{original.cache_origins},
    topology_origins{original.topology_origins},
    label_origins{original.label_origins},
    sys_projections{original.sys_projections},
    sys_projection_bounds{original.sys_projection_bounds},
    label_projections{original.label_projections},
    label_projection_bounds{original.label_projection_bounds},
    topology_projections{original.topology_projections},
    topology_projection_bounds{original.topology_projection_bounds},
    int_data{original.int_data},
    sc_ptr{original.sc_ptr},
    poly_ag_ptr{original.poly_ag_ptr},
    poly_ps_ptr{original.poly_ps_ptr}
{
  rebasePointers();
}

//-------------------------------------------------------------------------------------------------
SynthesisCacheMap::SynthesisCacheMap(SynthesisCacheMap &&original) :
    cache_system_count{original.cache_system_count},
    synthesis_system_count{original.synthesis_system_count},
    cache_label_count{original.cache_label_count},
    cache_topology_count{original.cache_topology_count},
    cache_origins{std::move(original.cache_origins)},
    topology_origins{std::move(original.topology_origins)},
    label_origins{std::move(original.label_origins)},
    sys_projections{std::move(original.sys_projections)},
    sys_projection_bounds{std::move(original.sys_projection_bounds)},
    label_projections{std::move(original.label_projections)},
    label_projection_bounds{std::move(original.label_projection_bounds)},
    topology_projections{std::move(original.topology_projections)},
    topology_projection_bounds{std::move(original.topology_projection_bounds)},
    int_data{std::move(original.int_data)},
    sc_ptr{original.sc_ptr},
    poly_ag_ptr{original.poly_ag_ptr},
    poly_ps_ptr{original.poly_ps_ptr}
{}

//-------------------------------------------------------------------------------------------------
SynthesisCacheMap& SynthesisCacheMap::operator=(const SynthesisCacheMap &other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }

  // Copy all members
  cache_system_count = other.cache_system_count;
  synthesis_system_count = other.synthesis_system_count;
  cache_label_count = other.cache_label_count;
  cache_topology_count = other.cache_topology_count;
  cache_origins = other.cache_origins;
  topology_origins = other.topology_origins;
  label_origins = other.label_origins;
  sys_projections = other.sys_projections;
  sys_projection_bounds = other.sys_projection_bounds;
  label_projections = other.label_projections;
  label_projection_bounds = other.label_projection_bounds;
  topology_projections = other.topology_projections;
  topology_projection_bounds = other.topology_projection_bounds;
  int_data = other.int_data;
  sc_ptr = other.sc_ptr;
  poly_ag_ptr = other.poly_ag_ptr;
  poly_ps_ptr = other.poly_ps_ptr;

  // Repair pointers and return
  rebasePointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
SynthesisCacheMap& SynthesisCacheMap::operator=(SynthesisCacheMap &&other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }

  // Copy all members
  cache_system_count = std::move(other.cache_system_count);
  synthesis_system_count = std::move(other.synthesis_system_count);
  cache_label_count = std::move(other.cache_label_count);
  cache_topology_count = std::move(other.cache_topology_count);
  cache_origins = std::move(other.cache_origins);
  topology_origins = std::move(other.topology_origins);
  label_origins = std::move(other.label_origins);
  sys_projections = std::move(other.sys_projections);
  sys_projection_bounds = std::move(other.sys_projection_bounds);
  label_projections = std::move(other.label_projections);
  label_projection_bounds = std::move(other.label_projection_bounds);
  topology_projections = std::move(other.topology_projections);
  topology_projection_bounds = std::move(other.topology_projection_bounds);
  int_data = std::move(other.int_data);
  sc_ptr = other.sc_ptr;
  poly_ag_ptr = other.poly_ag_ptr;
  poly_ps_ptr = other.poly_ps_ptr;
  return *this;
}

//-------------------------------------------------------------------------------------------------
int SynthesisCacheMap::getCacheSystemCount() const {
  return cache_system_count;
}

//-------------------------------------------------------------------------------------------------
int SynthesisCacheMap::getSynthesisSystemCount() const {
  return synthesis_system_count;
}

//-------------------------------------------------------------------------------------------------
int SynthesisCacheMap::getCacheLabelCount() const {
  return cache_label_count;
}

//-------------------------------------------------------------------------------------------------
int SynthesisCacheMap::getCacheTopologyCount() const {
  return cache_topology_count;
}

//-------------------------------------------------------------------------------------------------
int SynthesisCacheMap::getSynthesisTopologyCount() const {
  if (poly_ps_ptr != nullptr) {
    return poly_ps_ptr->getUniqueTopologyCount();
  }
  else if (poly_ag_ptr != nullptr) {
    return poly_ag_ptr->getUniqueTopologyCount();
  }
  else {
    rtErr("Neither a topology synthesis nor a coordinate synthesis is attached.",
          "SynthesisCacheMap", "getSynthesisTopologyCount");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SynthesisCacheMap::getSourceGroup(const int query_index) const {
  std::vector<int> result;
  if (query_index >= 0 && query_index >= cache_system_count) {
    const size_t llim = sys_projection_bounds.readHost(query_index);
    const size_t hlim = sys_projection_bounds.readHost(query_index + 1);
    const int* sys_ptr = sys_projections.data();
    for (size_t i = llim; i < hlim; i++) {
      result[i - llim] = sys_ptr[i];
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SynthesisCacheMap::getLabelGroup(const std::string &query_label) const {

  // Match the label to the cache
  const int cache_idx = sc_ptr->getLabelCacheIndex(query_label);
  std::vector<int> result;
  if (cache_idx < sc_ptr->getLabelCount()) {
    const size_t llim = label_projection_bounds.readHost(cache_idx);
    const size_t hlim = label_projection_bounds.readHost(cache_idx + 1);
    result.resize(hlim - llim);
    const int* lbl_ptr = label_projections.data();
    for (size_t i = llim; i < hlim; i++) {
      result[i - llim] = lbl_ptr[i];
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SynthesisCacheMap::getTopologyGroup(const AtomGraph *query_ag) const {

  // Match the topology to the cache
  const int cache_idx = sc_ptr->getTopologyCacheIndex(query_ag);
  std::vector<int> result;
  if (cache_idx < sc_ptr->getTopologyCount()) {
    const size_t llim = topology_projection_bounds.readHost(cache_idx);
    const size_t hlim = topology_projection_bounds.readHost(cache_idx + 1);
    result.resize(hlim - llim);
    const int* top_ptr = topology_projections.data();
    for (size_t i = llim; i < hlim; i++) {
      result[i - llim] = top_ptr[i];
    }
  }
  return result;  
}

//-------------------------------------------------------------------------------------------------
std::vector<int> SynthesisCacheMap::getTopologyGroup(const AtomGraph &query_ag) const {
  return getTopologyGroup(query_ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
int SynthesisCacheMap::getSystemCacheIndex(const int synthesis_index) const {
  return cache_origins.readHost(synthesis_index);
}

//-------------------------------------------------------------------------------------------------
int SynthesisCacheMap::getTopologyCacheIndex(const int synthesis_index) const {
  return topology_origins.readHost(synthesis_index);
}

//-------------------------------------------------------------------------------------------------
int SynthesisCacheMap::getLabelCacheIndex(const int synthesis_index) const {
  return label_origins.readHost(synthesis_index);
}

//-------------------------------------------------------------------------------------------------
int SynthesisCacheMap::getPartitionCount(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return sc_ptr->getSystemCount();
  case SystemGrouping::TOPOLOGY:
    if (poly_ag_ptr != nullptr) {
      return poly_ag_ptr->getUniqueTopologyCount();
    }
    if (poly_ps_ptr != nullptr) {
      return poly_ps_ptr->getUniqueTopologyCount();
    }
    break;
  case SystemGrouping::LABEL:
    return sc_ptr->getLabelCount();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int SynthesisCacheMap::getTotalProjection(int query_index, SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    if (query_index >= 0 && query_index < cache_system_count) {
      return sys_projection_bounds.readHost(query_index + 1) -
             sys_projection_bounds.readHost(query_index);
    }
    break;
  case SystemGrouping::TOPOLOGY:
    if (query_index >= 0 && query_index < cache_topology_count) {
      return topology_projection_bounds.readHost(query_index + 1) -
             topology_projection_bounds.readHost(query_index);
    }
    break;
  case SystemGrouping::LABEL:
    if (query_index >= 0 && query_index < cache_label_count) {
      return label_projection_bounds.readHost(query_index + 1) -
             label_projection_bounds.readHost(query_index);
    }
    break;
  }

  // If no quantity has been found thus far, return zero
  return 0;
}

//-------------------------------------------------------------------------------------------------
int SynthesisCacheMap::getTotalProjection(const AtomGraph *query_ag) const {
  const int top_idx = sc_ptr->getTopologyCacheIndex(query_ag);
  if (top_idx < cache_topology_count) {
    return topology_projection_bounds.readHost(top_idx + 1) -
           topology_projection_bounds.readHost(top_idx);
  }
  else {
    return 0;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int SynthesisCacheMap::getTotalProjection(const AtomGraph &query_ag) const {
  return getTotalProjection(query_ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
int SynthesisCacheMap::getTotalProjection(const std::string &query_label) const {
  const int lbl_idx = sc_ptr->getLabelCacheIndex(query_label);
  if (lbl_idx < cache_label_count) {
    return label_projection_bounds.readHost(lbl_idx + 1) -
           label_projection_bounds.readHost(lbl_idx);
  }
  else {
    return 0;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const SystemCache* SynthesisCacheMap::getCachePointer() const {
  return sc_ptr;
}

//-------------------------------------------------------------------------------------------------
const AtomGraphSynthesis* SynthesisCacheMap::getTopologySynthesisPointer() const {
  return poly_ag_ptr;
}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceSynthesis* SynthesisCacheMap::getCoordinateSynthesisPointer() const {
  return poly_ps_ptr;
}

//-------------------------------------------------------------------------------------------------
const SynthesisCacheMap* SynthesisCacheMap::getSelfPointer() const {
  return this;
}

//-------------------------------------------------------------------------------------------------
const SynthesisMapReader SynthesisCacheMap::data(const HybridTargetLevel tier) const {
  return SynthesisMapReader(cache_system_count, synthesis_system_count, cache_label_count,
                            cache_topology_count, cache_origins.data(tier),
                            topology_origins.data(tier), label_origins.data(tier),
                            sys_projections.data(tier), sys_projection_bounds.data(tier),
                            label_projections.data(tier), label_projection_bounds.data(tier),
                            topology_projections.data(tier),
                            topology_projection_bounds.data(tier));
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void SynthesisCacheMap::upload() {
  int_data.upload();
}

//-------------------------------------------------------------------------------------------------
void SynthesisCacheMap::download() {
  int_data.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void SynthesisCacheMap::setCache(const std::vector<int> &cache_origins_in,
                                 const SystemCache *sc_in) {
  sc_ptr = const_cast<SystemCache*>(sc_in);

  // Check that each system index in the cache origins is valid.
  cache_system_count = sc_ptr->getSystemCount();
  cache_topology_count = sc_ptr->getTopologyCount();
  cache_label_count = sc_ptr->getLabelCount();
  synthesis_system_count = cache_origins_in.size();
  for (int i = 0; i < synthesis_system_count; i++) {
    if (cache_origins_in[i] < 0 || cache_origins_in[i] >= cache_system_count) {
      rtErr("Cache system index " + std::to_string(cache_origins_in[i]) + " is invalid for a "
            "cache of " + std::to_string(cache_system_count) + " systems.", "SynthesisCacheMap",
            "setCache");
    }
  }

  // Compute the necessary size of Hybrid data storage, then allocate
  const size_t padded_chcsys_count = roundUp<size_t>(cache_system_count + 1, warp_size_zu);
  const size_t padded_synsys_count = roundUp<size_t>(synthesis_system_count, warp_size_zu);
  size_t n_int = 6LLU * padded_synsys_count;
  n_int += 3LLU * padded_chcsys_count;
  int_data.resize(n_int);
  
  // Construct the correspondence arrays
  std::vector<int> tmp_sys_projections(synthesis_system_count);
  std::vector<int> tmp_sys_projection_bounds(cache_system_count + 1, 0);
  std::vector<int> tmp_label_projections(synthesis_system_count);
  std::vector<int> tmp_label_projection_bounds(cache_system_count + 1, 0);
  std::vector<int> tmp_topology_projections(synthesis_system_count);
  std::vector<int> tmp_topology_projection_bounds(cache_system_count + 1, 0);
  std::vector<int> cached_label_indices(synthesis_system_count);
  std::vector<int> cached_topology_indices(synthesis_system_count);
  for (int i = 0; i < synthesis_system_count; i++) {
    cached_label_indices[i] = sc_ptr->getSystemLabelIndex(cache_origins_in[i]);
    cached_topology_indices[i] = sc_ptr->getSystemTopologyIndex(cache_origins_in[i]);
  }
  indexingArray(cache_origins_in, &tmp_sys_projections, &tmp_sys_projection_bounds);
  indexingArray(cached_label_indices, &tmp_label_projections, &tmp_label_projection_bounds);
  indexingArray(cached_topology_indices, &tmp_topology_projections,
                &tmp_topology_projection_bounds);
  
  // Load the index maps
  size_t ic = 0;
  ic = cache_origins.putHost(&int_data, cache_origins_in, ic, warp_size_zu);
  ic = topology_origins.putHost(&int_data, cached_topology_indices, ic, warp_size_zu);
  ic = label_origins.putHost(&int_data, cached_label_indices, ic, warp_size_zu);
  ic = sys_projections.putHost(&int_data, tmp_sys_projections, ic, warp_size_zu);
  ic = sys_projection_bounds.putHost(&int_data, tmp_sys_projection_bounds, ic, warp_size_zu);
  ic = topology_projections.putHost(&int_data, tmp_topology_projections, ic, warp_size_zu);
  ic = topology_projection_bounds.putHost(&int_data, tmp_topology_projection_bounds, ic,
                                          warp_size_zu);
  ic = label_projections.putHost(&int_data, tmp_label_projections, ic, warp_size_zu);
  ic = label_projection_bounds.putHost(&int_data, tmp_label_projection_bounds, ic, warp_size_zu);
}

//-------------------------------------------------------------------------------------------------
void SynthesisCacheMap::setCache(const std::vector<int> &cache_origins_in,
                                 const SystemCache &sc_in) {
  setCache(cache_origins_in, sc_in.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
void SynthesisCacheMap::setSynthesis(const AtomGraphSynthesis *poly_ag_in) {
  poly_ag_ptr = const_cast<AtomGraphSynthesis*>(poly_ag_in);
}

//-------------------------------------------------------------------------------------------------
void SynthesisCacheMap::setSynthesis(const AtomGraphSynthesis &poly_ag_in) {
  setSynthesis(poly_ag_in.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
void SynthesisCacheMap::setSynthesis(const PhaseSpaceSynthesis *poly_ps_in) {
  poly_ps_ptr = const_cast<PhaseSpaceSynthesis*>(poly_ps_in);
}

//-------------------------------------------------------------------------------------------------
void SynthesisCacheMap::setSynthesis(const PhaseSpaceSynthesis &poly_ps_in) {
  setSynthesis(poly_ps_in.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
void SynthesisCacheMap::rebasePointers() {
  cache_origins.swapTarget(&int_data);
  topology_origins.swapTarget(&int_data);
  label_origins.swapTarget(&int_data);
  sys_projections.swapTarget(&int_data);
  sys_projection_bounds.swapTarget(&int_data);
  label_projections.swapTarget(&int_data);
  label_projection_bounds.swapTarget(&int_data);
  topology_projections.swapTarget(&int_data);
  topology_projection_bounds.swapTarget(&int_data);
}

//-------------------------------------------------------------------------------------------------
void SynthesisCacheMap::validateCorrespondence() const {

  // Check that the cache origins of each system match those of the system in either synthesis.
  if (sc_ptr == nullptr) {
    rtErr("A systems cache must be referenced.", "SynthesisCacheMap", "validateCorrespondence");
  }
  if (poly_ps_ptr == nullptr && poly_ag_ptr == nullptr) {
    rtErr("A coordinate or topology synthesis must be referenced.", "SynthesisCacheMap",
          "validateCorrespondence");
  }
  const int* cache_orig_ptr = cache_origins.data();
  if (poly_ps_ptr != nullptr) {
    if (poly_ps_ptr->getSystemCount() != synthesis_system_count) {
      rtErr("The coordinate synthesis contains " + std::to_string(poly_ps_ptr->getSystemCount()) +
            " systems, whereas the map understands the origins of " +
            std::to_string(synthesis_system_count) + ".", "SynthesisCacheMap",
            "validateCorrespondence");
    }
    for (int i = 0; i < synthesis_system_count; i++) {
      const AtomGraph* cache_top = sc_ptr->getSystemTopologyPointer(cache_orig_ptr[i]);
      const AtomGraph* synth_top = poly_ps_ptr->getSystemTopologyPointer(i);
      if (cache_top != synth_top && cache_top->getFileName() != synth_top->getFileName()) {
        rtErr("Coordinate synthesis system " + std::to_string(i) + ", mapped to cache system " +
              std::to_string(cache_orig_ptr[i]) + ", references a topology originating in file " +
              synth_top->getFileName() + " while cache references a topology originating in "
              "file " + cache_top->getFileName() + ".", "SynthesisCacheMap",
              "validateCorrespondence");
      }
    }
  }
  if (poly_ag_ptr != nullptr) {
    if (poly_ag_ptr->getSystemCount() != synthesis_system_count) {
      rtErr("The topology synthesis contains " + std::to_string(poly_ps_ptr->getSystemCount()) +
            " systems, whereas the map understands the origins of " +
            std::to_string(synthesis_system_count) + ".", "SynthesisCacheMap",
            "validateCorrespondence");
    }
    for (int i = 0; i < synthesis_system_count; i++) {
      const AtomGraph* cache_top = sc_ptr->getSystemTopologyPointer(cache_orig_ptr[i]);
      const AtomGraph* synth_top = poly_ag_ptr->getSystemTopologyPointer(i);
      if (cache_top != synth_top && cache_top->getFileName() != synth_top->getFileName()) {
        rtErr("Topology synthesis system " + std::to_string(i) + ", mapped to cache system " +
              std::to_string(cache_orig_ptr[i]) + ", references a topology originating in file " +
              synth_top->getFileName() + " while cache references a topology originating in "
              "file " + cache_top->getFileName() + ".", "SynthesisCacheMap",
              "validateCorrespondence");
      }
    }
  }
}

} // namespace synthesis
} // namespace stormm
