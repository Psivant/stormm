// -*-c++-*-
#ifndef STORMM_SYNTHESIS_SYSTEM_CACHE_MAP
#define STORMM_SYNTHESIS_SYSTEM_CACHE_MAP

#include <string>
#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "phasespace_synthesis.h"
#include "atomgraph_synthesis.h"
#include "synthesis_enumerators.h"
#include "systemcache.h"

namespace stormm {
namespace synthesis {

using card::Hybrid;
using card::HybridTargetLevel;
using card::GpuDetails;
  
/// \brief A set of pointers to critical arrays in the map below.
struct SynthesisMapReader {

  /// \brief As with other abstracts, the constructor takes a list of input arguments for all
  ///        member variables.
  SynthesisMapReader(int ncache_in, int nsynth_in, int nlabel_in, int ntopol_in,
                     const int* cache_origins_in, const int* topology_origins_in,
                     const int* label_origins_in, const int* csystem_proj_in,
                     const int* csystem_bounds_in, const int* clabel_proj_in,
                     const int* clabel_bounds_in, const int* ctopol_proj_in,
                     const int* ctopol_bounds_in);

  /// \brief The copy and move constructors take their default forms, but copy and move assignment
  ///        operators are implicitly deleted due to the const-ness of all member variables.
  /// \{
  SynthesisMapReader(const SynthesisMapReader &original) = default;
  SynthesisMapReader(SynthesisMapReader &&original) = default;
  /// \}

  const int ncache;             ///< The number of unique systems in the cache (adding one to this
                                ///<   indicates the length of the csystem_bounds array)
  const int nsynth;             ///< The number of systems in the synthesis, and the length of
                                ///<   several arrays below
  const int nlabel;             ///< The number of unique labels in the cache (adding one to this
                                ///<   indicates the length of the clabel_bounds array)
  const int ntopol;             ///< The number of unique topologies in the cache (adding one to
                                ///<   this indicates the length of the ctopol_bounds array)
  const int* cache_origins;     ///< Cache system origin of each system in the synthesis
  const int* topology_origins;  ///< Indicies of topologies in the cache guiding each system in the
                                ///<   synthesis (the cache may hold a superset of all topologies
                                ///<   used by the synthesis)
  const int* label_origins;     ///< Cache label origin of each system in the synthesis
  const int* csystem_proj;      ///< Consecutive lists of indices for systems in the synthesis
                                ///<   based on each system in the cache
  const int* csystem_bounds;    ///< Bounds array for each cache system's projections into the
                                ///<   synthesis
  const int* clabel_proj;       ///< Consecutive lists of indices for systems in the synthesis
                                ///<   grouped under each label
  const int* clabel_bounds;     ///< Bounds array for each cache label group's projections into the
                                ///<   synthesis
  const int* ctopol_proj;       ///< Consecutive lists of indices for systems in the synthesis
                                ///<   making use of each topology in the cache
  const int* ctopol_bounds;     ///< Bounds array for each of the cached topologies' projections
                                ///<   into the synthesis
};
  
/// \brief Encode a map between the systems of a snythesis and those of a SystemCache.  The
///        synthesis is expected to replicate the systems present in the cache, perhaps in a
///        different order, perhaps to varying degrees for each system.
class SynthesisCacheMap {
public:

  /// \brief The constructor can work with a coordinate or topology synthesis, but in either case
  ///        requires a pairwise map between the synthesis and the cache.  If a synthesis is
  ///        provided, it will be used to check the correspondence.
  /// \{
  SynthesisCacheMap();
  
  SynthesisCacheMap(const std::vector<int> &cache_origins_in, const SystemCache *sc_in,
                    const AtomGraphSynthesis *poly_ag_in = nullptr,
                    const PhaseSpaceSynthesis *poly_ps_in = nullptr);

  SynthesisCacheMap(const std::vector<int> &cache_origins_in, const SystemCache &sc_in,
                    const AtomGraphSynthesis &poly_ag_in, const PhaseSpaceSynthesis &poly_ps_in);

  SynthesisCacheMap(const std::vector<int> &cache_origins_in, const SystemCache &sc_in,
                    const AtomGraphSynthesis &poly_ag_in);

  SynthesisCacheMap(const std::vector<int> &cache_origins_in, const SystemCache &sc_in,
                    const PhaseSpaceSynthesis &poly_ps_in);
  /// \}

  /// \brief The copy and move constructors, as well as assignment operators, must be explicitly
  ///        defined to handle repair of POINTER-kind Hybrid objects.
  ///
  /// \param original  The object to copy or move
  /// \param other     Another object to copy or move into the present one
  /// \{
  SynthesisCacheMap(const SynthesisCacheMap &original);
  SynthesisCacheMap(SynthesisCacheMap &&original);
  SynthesisCacheMap& operator=(const SynthesisCacheMap &original);
  SynthesisCacheMap& operator=(SynthesisCacheMap &&original);
  /// \}

  /// \brief Get the number of systems in the associated cache.
  int getCacheSystemCount() const;

  /// \brief Get the number of systems in the associated synthesis.
  int getSynthesisSystemCount() const;

  /// \brief Get the number of labels in the associated cache.
  int getCacheLabelCount() const;

  /// \brief Get the number of unique topologies in the associated cache.
  int getCacheTopologyCount() const;

  /// \brief Get the number of (unique) topologies in the associated synthesis.  This will check
  ///        the coordinate synthesis pointer, then the topology synthesis pointer, and return the
  ///        number of unique topologies based on the first valid pointer it encounters.
  int getSynthesisTopologyCount() const;

  /// \brief Get a list of all system indices in the synthesis derived from a particular -sys
  ///        keyword entry (a single system within the cache).
  ///
  /// \param query_index  The index of the system of interest within the systems cache
  std::vector<int> getSourceGroup(int query_index) const;
  
  /// \brief Get a list of all system indices in the synthesis derived from systems in the cache
  ///        matching a particular label.
  ///
  /// \param query_label  The label of interest
  std::vector<int> getLabelGroup(const std::string &query_label) const;

  /// \brief Get a list of all system indices in the synthesis derived from systems in the cache
  ///        matching a particular topology.
  ///
  /// Overloaded:
  ///   - Provide the topology by const pointer
  ///   - Provide the topology by const reference
  ///
  /// \param query_ag  The topology of interest
  /// \{
  std::vector<int> getTopologyGroup(const AtomGraph *query_ag) const;  
  std::vector<int> getTopologyGroup(const AtomGraph &query_ag) const;  
  /// \}

  /// \brief Get the system cache index of a system from the synthesis.
  ///
  /// \param synthesis_index  Index of the system of interest within the synthesis
  int getSystemCacheIndex(int synthesis_index) const;

  /// \brief Get the index of a topology in the cache guiding a system in the synthesis.  This
  ///        might not be the same as the unique topology index for the system referenced in the
  ///        synthesis itself, if the synthesis uses a subset of the cache's topologies.
  ///
  /// \param synthesis_index  Index of the system of interest within the synthesis
  int getTopologyCacheIndex(int synthesis_index) const;

  /// \brief Get the system label index of a system from the synthesis.
  ///
  /// \param synthesis_index  Index of the system of interest within the synthesis
  int getLabelCacheIndex(int synthesis_index) const;

  /// \brief Get the number of partitions that the synthesis should be divided into under any of
  ///        the methods for grouping systems.
  ///
  /// \param organization  The chosen method of grouping systems
  int getPartitionCount(SystemGrouping organization) const;
  
  /// \brief Get the number of systems in the synthesis associated with a particular system, label,
  ///        or topology in the cache.  If the index is invalid, this function will return zero.
  ///
  /// Overloaded:
  ///   - Indicate the system, label, or topology by its index in the cache.  The following
  ///     enumerator will determine how the index is interpreted.
  ///   - Indicate the topology by pointer or reference
  ///   - Provide the label string
  ///
  /// \param query_index     Index of the system, topology, or label within the cache
  /// \param topology_index  Index of the topology within the cache
  /// \param label_index     Index of the label within the cache
  /// \param query_label     Label string to seek out within the cache
  /// \param query_ag        Topology pointer to seek out within the cache
  /// \param oganization     The manner in which systems are to be grouped (the setting of this
  ///                        enumeration will determine the interpretation of preceding integer
  ///                        indices and only certain settings will be compatible with preceding
  ///                        label strings or topology pointers)
  /// \{
  int getTotalProjection(int query_index, SystemGrouping organization) const;
  int getTotalProjection(const AtomGraph *query_ag) const;
  int getTotalProjection(const AtomGraph &query_ag) const;
  int getTotalProjection(const std::string &query_label) const;
  /// \}

  /// \brief Get a const pointer to the system cache referenced by this map.
  const SystemCache* getCachePointer() const;

  /// \brief Get a const pointer to the coordinate synthesis referenced by this map.
  const AtomGraphSynthesis* getTopologySynthesisPointer() const;

  /// \brief Get a const pointer to the coordinate synthesis referenced by this map.
  const PhaseSpaceSynthesis* getCoordinateSynthesisPointer() const;

  /// \brief Get a const pointer to the object itself (useful when the object has been passed by
  ///        const reference).
  const SynthesisCacheMap* getSelfPointer() const;
  
  /// \brief Get an abstract of this map.
  const SynthesisMapReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

#ifdef STORMM_USE_HPC
  /// \brief Upload all data to the GPU device
  void upload();

  /// \brief Download all data from the GPU device
  void download();
#endif
  
  /// \brief Set the systems cache along with a list of correspondences.
  ///
  /// Overloaded:
  ///   - Provide a pointer to the cache
  ///   - Provide a reference to the cache
  ///
  /// \param cache_origins_in  System origins in the cache for each system in some synthesis
  /// \param sc_in             The cache to reference
  /// \{
  void setCache(const std::vector<int> &cache_origins_in, const SystemCache *sc_in);
  void setCache(const std::vector<int> &cache_origins_in, const SystemCache &sc_in);
  /// \}

  /// \brief Set the PhaseSpaceSynthesis or AtomGraphSynthesis pointers.  This will be useful if
  ///        the map is used to retrieve lists of conformations or coordinate sets from the
  ///        synthesis.
  ///
  /// Overloaded:
  ///   - Provide the synthesis by const pointer
  ///   - Provide the synthesis by const reference
  ///
  /// \param poly_ag_in  The topology synthesis of interest
  /// \param poly_ps_in  The coordinate synthesis of interest
  /// \{
  void setSynthesis(const AtomGraphSynthesis *poly_ag_in);
  void setSynthesis(const AtomGraphSynthesis &poly_ag_in);
  void setSynthesis(const PhaseSpaceSynthesis *poly_ps_in);
  void setSynthesis(const PhaseSpaceSynthesis &poly_ps_in);
  /// \}

private:
  int cache_system_count;                  ///< Number of systems in the cache
  int synthesis_system_count;              ///< Number of systems in the synthesis, and the
                                           ///<   dimension of many of the arrays below
  int cache_label_count;                   ///< Number of unique labels in the systems cache
  int cache_topology_count;                ///< Number of unique topologies in the systems cache
  Hybrid<int> cache_origins;               ///< Indices of systems in the cache that provided the
                                           ///<   basis for each system in the synthesis
  Hybrid<int> topology_origins;            ///< Indices of topologies in the cache that guide
                                           ///<   each system in the synthesis (the contents of
                                           ///<   this array may differ from those of the
                                           ///<   topology index array for the referenced
                                           ///<   synthesis, if the referenced synthesis does not
                                           ///<   use all topologies from the cache)
  Hybrid<int> label_origins;               ///< Indices of labels in the cache that provided the
                                           ///<   basis for each system in the synthesis
  Hybrid<int> sys_projections;             ///< Lists of systems in the synthesis originating with
                                           ///<   each system in the cache
  Hybrid<int> sys_projection_bounds;       ///< Bounds array for sys_projections
  Hybrid<int> label_projections;           ///< Lists of systems in the synthesis grouped under
                                           ///<   each label in the cache
  Hybrid<int> label_projection_bounds;     ///< Bounds array for label_projections
  Hybrid<int> topology_projections;        ///< Lists of systems in the synthesis governed by each
                                           ///<   topology from the cache, given in the order in
                                           ///<   which topologies appear in the cache.  This might
                                           ///<   not be the same as the apping of systems to
                                           ///<   unique topologies presented in the accompanying
                                           ///<   PhaseSpaceSynthesis.
  Hybrid<int> topology_projection_bounds;  ///< Bounds array for topology_projections
  Hybrid<int> int_data;                    ///< An ARRAY-kind Hybrid object targeted by the
                                           ///<   POINTER-kind Hybrids above

  // Pointers to the objects being linked
  SystemCache *sc_ptr;               ///< The systems cache
  AtomGraphSynthesis *poly_ag_ptr;   ///< The topology synthesis
  PhaseSpaceSynthesis *poly_ps_ptr;  ///< The coordinate synthesis

  /// \brief Repair pointers after copy or copy assignment.
  void rebasePointers();

  /// \brief Validate the synthesis and cache system assignments.
  void validateCorrespondence() const;
};
  
} // namespace synthesis
} // namespace stormm

#endif
