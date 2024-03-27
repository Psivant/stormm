// -*-c++-*-
#ifndef STORMM_COMPARISON_GUIDE_H
#define STORMM_COMPARISON_GUIDE_H

#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/systemcache.h"
#include "Synthesis/synthesis_cache_map.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinate_series.h"

namespace stormm {
namespace analysis {

using card::GpuDetails;
using card::Hybrid;
using card::HybridKind;
using card::HybridTargetLevel;
using synthesis::Condensate;
using synthesis::PhaseSpaceSynthesis;
using synthesis::StructureSource;
using synthesis::SynthesisCacheMap;
using synthesis::SystemGrouping;
using topology::AtomGraph;
using trajectory::CoordinateSeries;
  
/// \brief Read-only abstract for the AnalysisGuide below.
struct CompGuideKit {

  /// \brief The constructor takes a straight list of all relevant constants and pointers.
  CompGuideKit(int system_count_in, int source_count_in, int topology_count_in, int label_count_in,
               int natr_insr_src_in, int natr_insr_top_in, int natr_insr_lbl_in,
               int nata_insr_src_in, int nata_insr_top_in, int nata_insr_lbl_in,
               const int4* atr_member_src_in, const int4* atr_member_top_in,
               const int4* atr_member_lbl_in, const int* atr_group_src_in,
               const int* atr_group_top_in, const int* atr_group_lbl_in,
               const int4* ata_member_src_in, const int4* ata_member_top_in,
               const int4* ata_member_lbl_in, const int* ata_group_src_in,
               const int* ata_group_top_in, const int* ata_group_lbl_in, const int* src_groups_in,
               const int* top_groups_in, const int* lbl_groups_in, const int* src_grp_ag_in,
               const int* top_grp_ag_in, const int* lbl_grp_ag_in, const int* src_group_bounds_in,
               const int* top_group_bounds_in, const int* lbl_group_bounds_in,
               const size_t* se_pair_src_offsets_in, const size_t* se_pair_top_offsets_in,
               const size_t* se_pair_lbl_offsets_in, const size_t* nr_pair_src_offsets_in,
               const size_t* nr_pair_top_offsets_in, const size_t* nr_pair_lbl_offsets_in);

  /// \brief The default copy and move constructors can apply, but the copy and move assignment
  ///        operators will be implicitly deleted due to the presence of const members.
  /// \{
  CompGuideKit(const CompGuideKit &original) = default;
  CompGuideKit(CompGuideKit &&original) = default;
  /// \}

  const int system_count;    ///< The number of systems in whatever underlying collection
  const int source_count;    ///< The number of cache sources found in associated user input
  const int topology_count;  ///< The number of unique topologies needed to describe all systems
                             ///<   (if no user data or systems cache is provided, each unique
                             ///<   topology will become its own cache source and label group)
  const int label_count;     ///< The number of unique labels obtained from user input (or inferred
                             ///<   inferred from deficient input)
  const int natr_insr_src;   ///< Number of block instructions spanning the coordinate synthesis in
                             ///<   all-to-one analyses for groups of systems originating in the
                             ///<   same systems cache source
  const int natr_insr_top;   ///< Number of block instructions for all-to-one analyses on all
                             ///<   groups of systems making use of the same topology
  const int natr_insr_lbl;   ///< Number of block instructions for all-to-one analyses on all
                             ///<   groups of systems making use of a specific cache label
  const int nata_insr_src;   ///< Number of block instructions spanning the coordinate synthesis in
                             ///<   all-to-all analyses for groups of systems originating in the
                             ///<   same systems cache source
  const int nata_insr_top;   ///< Number of block instructions for all-to-all analyses on groups of
                             ///<   systems making use of the same topology
  const int nata_insr_lbl;   ///< Number of block instructions for all-to-one analyses on all
                             ///<   groups of systems making use of a specific cache label

  // The following arrays provide instructions for operations on various types of structure
  // groupings.  Each instruction should be thought of as a five-integer tuple, with the int4 and
  // associated int being read one after another in whatever code is following this guide. 
  const int4* atr_member_src;  ///< Block instructions for all-to-one (All To Reference) analyses
                               ///<   concerning groups originating in each systems cache source
  const int4* atr_member_top;  ///< Block instructions for all-to-one (All To Reference) analyses
                               ///<   concerning groups of systems using the same topology
  const int4* atr_member_lbl;  ///< Block instructions for all-to-one (All To Reference) analyses
                               ///<   concerning groups of systems using the same cache label
  const int* atr_group_src;    ///< Cache system group indices to put the above instructions in
                               ///<   context
  const int* atr_group_top;    ///< Topology group indices to put the above instructions in
                               ///<   context
  const int* atr_group_lbl;    ///< Cache label group indices to put the above instructions in
                               ///<   context
  const int4* ata_member_src;  ///< Block instructions for _A_ll _T_o _A_ll analyses concerning
                               ///<   groups originating in each systems cache source
  const int4* ata_member_top;  ///< Block instructions for _A_ll _T_o _A_ll analyses concerning
                               ///<   groups of systems using the same topology
  const int4* ata_member_lbl;  ///< Block instructions for _A_ll _T_o _A_ll analyses concerning
                               ///<   groups groups of systems using the same cache label
  const int* ata_group_src;    ///< Cache system group indices to put the above instructions in
                               ///<   context
  const int* ata_group_top;    ///< Topology group indices to put the above instructions in
                               ///<   context
  const int* ata_group_lbl;    ///< Cache label group indices to put the above instructions in
                               ///<   context

  // The following arrays provide context for the above instructions.  The instructions themselves
  // contain information about a group and specific structures within that group, but the actual
  // contents of the group are opaque without an understanding of which structures within the
  // collection (particularly if it is a coordinate synthesis) make up the group.  In various
  // contexts, some of these arrays may not hold meaningful data.  A CoordinateSeries will only
  // have a list of incrementing integers 0, 1, 2, ..., N-1 for for all of the (...)_groups
  // variables, and a list of zeros for all of the (...)_ag variables.  If a synthesis was used to
  // construct the parent object with no cache map, then src_(...) and lbl_(...) array will
  // match the top_(...) array (every unique topology will be posited to have come from one cache
  // source with its own unique label).
  const int* src_groups;      ///< Lists of structures derived from the same cache source
  const int* top_groups;      ///< Lists of structures making use of the same topology
  const int* lbl_groups;      ///< Lists of structures falling under the same cache label group
  const int* src_grp_ag;      ///< Toplogies for each structure indexed in src_groups
  const int* top_grp_ag;      ///< Toplogies for each structure indexed in top_groups
  const int* lbl_grp_ag;      ///< Toplogies for each structure indexed in lbl_groups
  const int* src_grp_bounds;  ///< Bounds array for src_groups and src_grp_ag
  const int* top_grp_bounds;  ///< Bounds array for top_groups and top_grp_ag
  const int* lbl_grp_bounds;  ///< Bounds array for lbl_groups and lbl_grp_ag

  // While the results for all-to-one comparisons are expected to be placed in arrays of length
  // equal to the total number of systems in the synthesis / series collection (compartmentalized
  // by the various (...)_group_bounds arrays above for each respective partitioning of
  // structures), all-to-all comparisons require substantially more data to store the results and
  // the bounds between results for each grouping require somewhat more calculation.  Furthermore,
  // there might be a question as to whether the comparison is symmetric and can thus assume that
  // the result of comparing structure i to structure j in group g is the same as the comparison
  // of structure j to structure i.  The arrays below provide the boundaries for such results.
  // Total allocations are indicated by member variables in the parent object, accessible through
  // getter functions on the host CPU.
  const size_t* se_pair_src_offsets;  ///< Symmetrically equivalent pair result offsets for cache
                                      ///<   source structure groupings
  const size_t* se_pair_top_offsets;  ///< Symmetrically equivalent pair result offsets for unique
                                      ///<   topology structure groupings
  const size_t* se_pair_lbl_offsets;  ///< Symmetrically equivalent pair result offsets for cache
                                      ///<   label structure groupings
  const size_t* nr_pair_src_offsets;  ///< Non-reflexive pair result offsets for cache source
                                      ///<   structure groupings
  const size_t* nr_pair_top_offsets;  ///< Non-reflexive pair result offsets for unique topology
                                      ///<   structure groupings
  const size_t* nr_pair_lbl_offsets;  ///< Non-reflexive pair result offsets for cache label
                                      ///<   structure groupings  
};
  
/// \brief A means for computing and storing work instructions to guide comparisons of structures
///        in one of STORMM's multi-system collections (e.g. a PhaseSpaceSynthesis or
///        CoordinateSeries).  Different ways to build the object may or may not permit different
///        groupings of the structures in the collection.
class ComparisonGuide {
public:

  /// \brief Construction can take place using a coordinate synthesis, with or without an
  ///        associated cache map.  It can also take a coordinate series and treat it as a sort
  ///        of synthesis based entirely on a single topology and cache source.
  /// \{
  ComparisonGuide();

  ComparisonGuide(const PhaseSpaceSynthesis *poly_ps_in, const GpuDetails &gpu = null_gpu);

  ComparisonGuide(const PhaseSpaceSynthesis &poly_ps_in, const GpuDetails &gpu = null_gpu);

  ComparisonGuide(const PhaseSpaceSynthesis *poly_ps_in, const SynthesisCacheMap *scmap_in,
                  const GpuDetails &gpu = null_gpu);

  ComparisonGuide(const PhaseSpaceSynthesis &poly_ps_in, const SynthesisCacheMap &scmap_in,
                  const GpuDetails &gpu = null_gpu);
  
  ComparisonGuide(const Condensate *cdns_in, const GpuDetails &gpu = null_gpu);

  ComparisonGuide(const Condensate &cdns_in, const GpuDetails &gpu = null_gpu);

  ComparisonGuide(const Condensate *cdns_in, const AtomGraph *ag_in,
                  const GpuDetails &gpu = null_gpu);

  ComparisonGuide(const Condensate &cdns_in, const AtomGraph &ag_in,
                  const GpuDetails &gpu = null_gpu);

  ComparisonGuide(const Condensate *cdns_in, const SynthesisCacheMap *scmap_in,
                  const GpuDetails &gpu = null_gpu);
  
  ComparisonGuide(const Condensate &cdns_in, const SynthesisCacheMap &scmap_in,
                  const GpuDetails &gpu = null_gpu);

  template <typename T>
  ComparisonGuide(const CoordinateSeries<T> *cs_in, const AtomGraph *ag_in,
                  const GpuDetails &gpu = null_gpu);

  template <typename T>
  ComparisonGuide(const CoordinateSeries<T> &cs_in, const AtomGraph &ag_in,
                  const GpuDetails &gpu = null_gpu);
  /// \}

  /// \brief The copy and move constructors as well as copy and move assignment operators require
  ///        explicit instantiations.
  ///
  /// \param original  An existing object to copy or move
  /// \param other     An existing object to fill the right hand side of an assignment statement
  /// \{
  ComparisonGuide(const ComparisonGuide &original);
  ComparisonGuide(ComparisonGuide &&original);
  ComparisonGuide& operator=(const ComparisonGuide &original);
  ComparisonGuide& operator=(ComparisonGuide &&original);
  /// \}
  
  /// \brief Get the total number of systems involved in the comparisons outlined herein.
  int getSystemCount() const;

  /// \brief Get the number of unique topologies referenced by comparisons outlined in this object.
  int getTopologyCount() const;

  /// \brief Get the basis of the collection of structures.
  StructureSource getBasis() const;

  /// \brief Get the number of system groupings based on any of the available breakdowns.
  ///
  /// \param organization  The manner in which structures of the collection are to be grouped
  int getPartitionCount(SystemGrouping organization) const;

  /// \brief Get the array size to allocate for all-to-one calculations across the collection of
  ///        systems.
  size_t getAllToReferenceOutputSize() const;

  /// \brief Get the array size to allocate for symmetry-equivalent all-to-all calculations across
  ///        the collection of systems.
  ///
  /// \param organization  The manner in which structures of the collection are to be grouped
  size_t getSymmetryEquivalentPairOutputSize(SystemGrouping organization) const;
  
  /// \brief Get the array size to allocate for non-reflexive all-to-all calculations across the
  ///        collection of systems.
  ///
  /// \param organization  The manner in which structures of the collection are to be grouped
  size_t getNonReflexivePairOutputSize(SystemGrouping organization) const;
  
  /// \brief Get the number of all-to-one instructions for a particular subdivision of structures
  ///        in the collection.
  ///
  /// \param organization  The manner in which structures of the collection are to be grouped
  int getATRInstructionCount(SystemGrouping organization) const;

  /// \brief Get the number of all-to-all instructions for a particular subdivision of structures
  ///        in the collection.
  ///
  /// \param organization  The manner in which structures of the collection are to be grouped
  int getATAInstructionCount(SystemGrouping organization) const;

  /// \brief Get information pertaining to the structure members of each group involved in every
  ///        all-to-one comparison instruction for a given subdivision of the structures in the
  ///        collection.  This performs a deep copy of part of the instructions for use in
  ///        CPU-bound applications, an alternative to getting the information from the abstract.
  ///
  /// \param organization  The manner in which structures of the collection are to be grouped
  std::vector<int4> getATRInstructionMembers(SystemGrouping organization) const;

  /// \brief Get information pertaining to the groups involved in each all-to-one comparison
  ///        instruction.  The kth element of the array returned by this function combines with the
  ///        kth element of the four-tuple array returned by getATRInstructionMembers() to complete
  ///        the specification of systems to compare.
  ///
  /// \param organization  The manner in which structures of the collection are to be grouped
  std::vector<int> getATRInstructionGroups(SystemGrouping organization) const;

  /// \brief Get information pertaining to the structure members of each group involved in every
  ///        all-to-all comparison instruction for a given subdivision of the structures in the
  ///        collection.  This performs a deep copy of part of the instructions for use in
  ///        CPU-bound applications, an alternative to getting the information from the abstract.
  ///
  /// \param organization  The manner in which structures of the collection are to be grouped
  std::vector<int4> getATAInstructionMembers(SystemGrouping organization) const;

  /// \brief Get information pertaining to the groups involved in each all-to-all comparison
  ///        instruction.  The kth element of the array returned by this function combines with the
  ///        kth element of the four-tuple array returned by getATAInstructionMembers() to complete
  ///        the specification of systems to compare.
  ///
  /// \param organization  The manner in which structures of the collection are to be grouped
  std::vector<int> getATAInstructionGroups(SystemGrouping organization) const;

  /// \brief Get the synthesis (or series) frame indices of structures involved in each group.
  ///
  /// Overloaded:
  ///   - Get the structure indices of a specific partition for the given grouping method
  ///   - Get the structure indices of all partitions within a given grouping method, concatenated
  ///
  /// \param partition     Index of the partition for the specific grouping method
  /// \param organization  The manner in which structures of the collection are to be grouped
  /// \{
  std::vector<int> getGroupStructureIndices(int partition, SystemGrouping organization) const;
  std::vector<int> getGroupStructureIndices(SystemGrouping organization) const;
  /// \}

  /// \brief Get the indices of topologies governing the structures in each group.  The topologies
  ///        are assumed to be held in some auxiliary array, such as the PhaseSpaceSynthesis member
  ///        variable unique_topologies.  A series references only a single topology and if that is
  ///        the basis of the ComparsonGuide then the output of this function will be an array of
  ///        index zero.
  ///
  /// Overloaded:
  ///   - Get the topology indices of a specific partition for the given grouping method
  ///   - Get the topology indices of all partitions within a given grouping method, concatenated
  ///
  /// \param partition     Index of the partition for the specific grouping method
  /// \param organization  The manner in which structures of the collection are to be grouped
  /// \{
  std::vector<int> getGroupTopologyIndices(int partition, SystemGrouping organization) const;
  std::vector<int> getGroupTopologyIndices(SystemGrouping organization) const;
  /// \}
  
  /// \brief Get the limits of each group in the associated arrays of structure and topology
  ///        indices.
  ///
  /// \param organization  The manner in which structures of the collection are to be grouped
  std::vector<int> getGroupBounds(SystemGrouping organization) const;  

  /// \brief Get the number of structures occupying a specific, numbered partition of one of the
  ///        system groupings, e.g. the number of systems governed by topoogy index 5.
  ///
  /// \param index         The 
  /// \param organization  The manner in which structures of the collection are to be grouped
  
  int getFrameCount(int index, SystemGrouping organizaiton) const;
  
  /// \brief Get the offset for recording results of all-to-one comparisons within a given
  ///        partition of one structure grouping.
  ///
  /// \param index         Index of the partition within the group, e.g. the result offset for all
  ///                      structures governed by the 3rd topology compared to one member of that
  ///                      same partition
  /// \param organization  The manner in which structures of the collection are to be grouped
  size_t getAllToOneResultOffset(int index, SystemGrouping organization) const;

  /// \brief Get the offset for recording results of symmetric all-to-all comparisons within a
  ///        given partition of one structure grouping.
  ///
  /// \param index         Index of the partition within the group, e.g. the result offset for all
  ///                      structures governed by the 3rd topology compared to one member of that
  ///                      same partition
  /// \param organization  The manner in which structures of the collection are to be grouped
  size_t getSymmetryEquivalentResultOffset(int index, SystemGrouping organization) const;

  /// \brief Get the offset for recording results of non-symmetric (non-reflexive) all-to-all
  ///        comparisons within a given partition of one structure grouping.
  ///
  /// \param index         Index of the partition within the group, e.g. the result offset for all
  ///                      structures governed by the 3rd topology compared to one member of that
  ///                      same partition
  /// \param organization  The manner in which structures of the collection are to be grouped
  size_t getNonReflexiveResultOffset(int index, SystemGrouping organization) const;

  /// \brief Get a const pointer to the original synthesis.  This will error out if the object is
  ///        not based on a PhaseSpaceSynthesis.
  const PhaseSpaceSynthesis* getPhaseSpaceSynthesisPointer() const;

  /// \brief Get a const pointer to the original synthesis.  This will error out if the object is
  ///        not based on a Condensate with its own coordinate data.
  const Condensate* getCondensatePointer() const;

  /// \brief Get a const pointer to the map between the synthesis that this object is based upon
  ///        and an underlying systems cache used in construction of the synthesis.  This will
  ///        error out if the object is not basd on a synthesis with an associated map.
  const SynthesisCacheMap* getSynthesisMapPointer() const;

  /// \brief Get a const pointer to the original coordinate series.  The developer must supply the
  ///        template type for the CoordinateSeries in order to re-interpret the arbitrarily stored
  ///        <int> type for the CoordinateSeries pointer.
  template <typename T> const CoordinateSeries<T>* getSeriesPointer() const;

  /// \brief Get the data type of the CoordinateSeries upon which this object is based.
  size_t getCoordinateSeriesTypeID() const;

  /// \brief Get one of the topologies referenced by this guide.
  ///
  /// Overloaded:
  ///   - Get a topology from the list of unique objects based solely on its index in that list
  ///   - Get the topology governing the jth member of the ith partition of the systems according
  ///     to some system grouping method
  ///
  /// \param index         Index of the topology within the unique list
  /// \param organization  The method of grouping systems in the synthesis (or series, but grouping
  ///                      systems in a series would have to be done under carefully arranged
  ///                      circumstances in order to have any meaning)
  /// \param partition     A partition within the grouping method of choice
  /// \param member        A member of the partition to access (most members of the same partition
  ///                      will share the same topology, but for label groups, different members
  ///                      can be governed by different topologies)
  /// \{
  const AtomGraph* getTopologyPointer(int index) const;
  
  const AtomGraph* getTopologyPointer(SystemGrouping organization, int partition,
                                      int member) const;
  /// \}
  
  /// \brief Get the abstract.
  ///
  /// \param tier  Retrieve pointers on the CPU host or GPU device
  CompGuideKit data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

#ifdef STORMM_USE_HPC
  /// \brief Upload data to the GPU device.
  void upload();

  /// \brief Download data to the CPU host.
  void download();
#endif

private:
  int system_count;                ///< The number of systems in whatever underlying collection
  int source_count;                ///< The number of cache sources found in associated user input
  int topology_count;              ///< The number of unique topologies needed to describe all
                                   ///<   systems (if no user data or systems cache is provided,
                                   ///<   each unique topology will become its own cache source and
                                   ///<   label group)
  int label_count;                 ///< The number of unique labels obtained from user input (or
                                   ///<   inferred from deficient input)
  StructureSource basis;           ///< Indicates the nature and maybe the data type of the origin
                                   ///<   of the structures to be compared
  size_t csptr_data_type;          ///< The original data type of the CoordinateSeries pointer,
                                   ///<   stored here for internal reference.  The pointer itself
                                   ///<   (see cs_ptr, below) is re-cast to an arbitrary <int>
                                   ///<   type to prevent the ComparisonGuide object itself from
                                   ///<   taking on a template type requirement.
  int atr_instruction_count_src;   ///< The number of instructions needed to span the coordinate
                                   ///<   synthesis in all-to-reference (ATR) analyses for systems
                                   ///<   grouped under the same cache source.
  int atr_instruction_count_top;   ///< The number of instructions needed to span the coordinate
                                   ///<   synthesis in all-to-reference (ATR) analyses for systems
                                   ///<   grouped under the same unique topology.
  int atr_instruction_count_lbl;   ///< The number of instructions needed to span the coordinate
                                   ///<   synthesis in all-to-reference (ATR) analyses for systems
                                   ///<   grouped under the same cache label.
  int ata_instruction_count_src;   ///< The number of instructions needed to span the coordinate
                                   ///<   synthesis in all-to-all (ATA) analyses for systems
                                   ///<   grouped under the same cache source.
  int ata_instruction_count_top;   ///< The number of instructions needed to span the coordinate
                                   ///<   synthesis in all-to-all (ATA) analyses for systems
                                   ///<   grouped under the same unique topology.
  int ata_instruction_count_lbl;   ///< The number of instructions needed to span the coordinate
                                   ///<   synthesis in all-to-all (ATA) analyses for systems
                                   ///<   grouped under the same cache label.

  /// Block instructions for all-to-reference calculations with various methods of grouping the
  /// systems.  Extensions of each array name correspond to those seen in the SynthesisCacheMap
  /// object (src = source / -sys { ... } keyword entry, top = topology, lbl = cache label).  The
  /// reference structures corresponding to each group will be supplied by an auxiliary array of
  /// integers giving the structures' absolute indices within the synthesis.  Each instruction
  /// gives the index of one structure in the synthesis in its "x" member, the index of the
  /// reference structure in the auxiliary array in its "y" member, provides the unique topology
  /// index of all systems to be accessed (also the RMSD plan index) in the "z" member, and the
  /// number of the structure within its respective group in the "w" member.
  /// \{
  Hybrid<int4> atr_instruction_members_src;
  Hybrid<int4> atr_instruction_members_top;
  Hybrid<int4> atr_instruction_members_lbl;
  /// \}

  /// The instructions above provide offsets within specific system groupings and the indices of
  /// topologies to work from, but the group index is missing.  The following arrays provide group
  /// indices to put each of the instructions above in context.
  /// \{
  Hybrid<int> atr_instruction_groups_src;
  Hybrid<int> atr_instruction_groups_top;
  Hybrid<int> atr_instruction_groups_lbl;
  /// \}

  /// Block instructions for all-to-all calculations with various methods of grouping the systems.
  /// Extensions of each array name correspond to those seen in the SynthesisCacheMap object
  /// (src = source / -sys { ... } keyword entry, top = topology, lbl = cache label).  Each
  /// instruction contains the group indices of two structures in the "x" and "y" members,
  /// which shall be compared according to the topology of index given in the "z" member.  The
  /// low and high 16 bits of the "w" member provide the quantity of additional systems starting
  /// from the "x" and "y" group indices, respectively, that are to participate in the same
  /// comparison.
  /// \{
  Hybrid<int4> ata_instruction_members_src;
  Hybrid<int4> ata_instruction_members_top;
  Hybrid<int4> ata_instruction_members_lbl;
  /// \}

  /// The instructions above provide offsets within specific system groupings and the indices of
  /// topologies to work from, but the group index is missing.  The following arrays provide group
  /// indices to put each of the instructions above in context.
  /// \{
  Hybrid<int> ata_instruction_groups_src;
  Hybrid<int> ata_instruction_groups_top;
  Hybrid<int> ata_instruction_groups_lbl;
  /// \}

  /// Arrays to hold concatenated lists of indices for systems in eahc group within the collection
  /// of systems (coordinate synthesis or coordinate series) to which this ComparisonGuide is
  /// attached.
  /// \{
  Hybrid<int> cache_source_groups;
  Hybrid<int> topology_groups;
  Hybrid<int> cache_label_groups;
  /// \}

  /// Arrays of topology indices associated with each structure index given in one of the above
  /// three arrays for systems grouped by cache source, the underlying topology, or cache label.
  /// These topology indices refer to the collection of topologies associated with whatever
  /// coordinate synthesis this ComparisonGuide references.  If the guide instead references a
  /// coordinate series, all data in these arrays will read zero to reference the single topology
  /// descriptive of every structure in the series.  While it would be possible to compress the
  /// information by realizing that all structures in the same cache source group will be
  /// associated with the same topology and that all structures in a topology group will reference
  /// the topology with an index identical to the group itself, it is helpful to have cognate
  /// arrays in each category which allow the same methods of access.
  /// \{
  Hybrid<int> cache_source_group_topologies;
  Hybrid<int> topology_group_topologies;
  Hybrid<int> cache_label_group_topologies;
  /// \}

  /// Bounds arrays for the various (...)_groups and (...)_group_topologies arrays above.  These
  /// arrays indicate the various sizes of each group of systems.
  /// \{
  Hybrid<int> cache_source_group_bounds;
  Hybrid<int> topology_group_bounds;
  Hybrid<int> cache_label_group_bounds;
  /// \}

  /// An ARRAY-kind Hybrid object targeted by various <int>-type POINTER-kind Hybrids above
  Hybrid<int> int_data;

  /// The total sizes to allocate for symmetrically equivalent (se) and non-reflexive (nr)
  /// all-to-all comparisons.  Arrays below contain offsets for each grouping of systems
  /// (src = cache source, top = topology, lbl = cache label group) within these allocated sizes.
  /// \{
  size_t sepair_result_alloc_src;
  size_t sepair_result_alloc_top;
  size_t sepair_result_alloc_lbl;
  size_t nrpair_result_alloc_src;
  size_t nrpair_result_alloc_top;
  size_t nrpair_result_alloc_lbl;
  /// \}
  
  /// Arrays to hold the offsets for all-to-all results comparing structure i to structure j in
  /// group g, for all j < i.  To the gth offset, code must add ((i - 1) * i / 2) + j to obtain
  /// the array index into which some result shall be stored, whatever the analysis and whatever
  /// the data type of the result.  This will make space for all symmetrically equivalent (se)
  /// pairs of systems in each method of grouping the structures.  Results for comparisons within
  /// different groups are padded by the warp size, for consistency with other data storage
  /// approaches.  Arrays for all-to-one results follow from the bounds arrays for each grouping
  /// (e.g. the (...)_bounds arrays above).
  /// \{
  Hybrid<size_t> sepair_result_offsets_src;
  Hybrid<size_t> sepair_result_offsets_top;
  Hybrid<size_t> sepair_result_offsets_lbl;
  /// \}

  /// Arrays to hold offsets for all-to-all results comparing structure i to structure j in group
  /// group g, for all i and j (both i and j span all structures in group g).  These reflexive,
  /// square (sq) group bounds must be used when the result of comparing structures i -> j might
  /// differ from the result of comparing structures j -> i.
  /// \{
  Hybrid<size_t> nrpair_result_offsets_src;
  Hybrid<size_t> nrpair_result_offsets_top;
  Hybrid<size_t> nrpair_result_offsets_lbl;  
  /// \}

  /// An ARRAY-kind Hybrid object targeted by various <size_t>-type POINTER-kind Hybrids above
  Hybrid<size_t> sizet_data;
  
  // The following are pointers to objects containing coordinates for the comparisons guided by
  // this object.
  PhaseSpaceSynthesis *pps_ptr;   ///< A coordinate synthesis
  Condensate *cdns_ptr;           ///< A Condensate object can hold either a coordinate synthesis
                                  ///<   or series, and must be checked to properly set the basis
                                  ///<   member variable.
  SynthesisCacheMap *scmap_ptr;   ///< A map relating the synthesis to a separate collection of
                                  ///<   systems built from user input
  CoordinateSeries<int> *cs_ptr;  ///< A coordinate series to hold the coordinates.  This is recast
                                  ///<   to <int> in similar fashion to what is done for the
                                  ///<   Condensate object, to negate the need for templating.  The
                                  ///<   original type is recorded in the csptr_data_type member
                                  ///<   variable.

  /// A list of the relevant topologies, in an order valid to the indexing in the various
  /// (...)_group_topologies member variable arrays above.
  std::vector<AtomGraph*> ag_pointers;

  /// \brief Allocate date for system indexing.  This can also serve to repair POINTER-kind Hybrid
  ///        arrays in copy and move operations.
  void allocateSystemIndexing();
  
  /// \brief Record system indexing based on the coordinate source and its possible partitions.
  ///        This routine regularizes the various inputs and lays the foundation for a coherent
  ///        plan to compare structures from the whole set.
  void setSystemIndexing();
  
  /// \brief Generate a set of work units for a particular grouping of the systems in the
  ///        synthesis (or, in a special case, the series).
  ///
  /// \param system_list            Concatenated lists of system indices (referring to the order of
  ///                               each structure's appearance in the associated collection)
  ///                               describing each group (e.g. all structures derived from the
  ///                               same cache label)
  /// \param topology_index_list    List of indices for the topologies governing each structure
  ///                               under the present method of grouping.  These topology indices
  ///                               refer to an array containing all of the unique topologies
  ///                               needed to describe structures in the associated coordinate
  ///                               synthesis (a coordinate series is governed by one topology).
  /// \param bounds_list            Bounds array for system_list, and also implicity for
  ///                               topology_index_list
  /// \param partitions             The number of groupings into which all systems of the synthesis
  ///                               (or series) are divided.  A series can only be divided into one
  ///                               meaningful group.  The trusted length of the bounds_list member
  ///                               variable is partitions + 1.
  /// \param atr_insr_members       Systems (indexed by position in a group) and topologies
  ///                               needed for all-to-one instructions
  /// \param atr_insr_groups        Groups pertaining to each all-to-one instruction.  This array
  ///                               of integers completes the information put forth in the
  ///                               associated array containing information about (group) members.
  /// \param ata_insr_members       Systems (indexed by position in a group) and topologies
  ///                               needed for all-to-all instructions
  /// \param ata_insr_groups        Groups pertaining to each all-to-all instruction.  This array
  ///                               of integers completes the information put forth in the
  ///                               associated array containing information about (group) members.
  /// \param atr_instruction_count  Total number of all-to-one instructions (set and returned)
  /// \param ata_instruction_count  Total number of all-to-all instructions (set and returned)
  /// \param gpu                    Specifications of the GPU available
  void generateWorkUnits(const int* system_list, const int* topology_index_list,
                         const int* bounds_list, int partitions, Hybrid<int4> *atr_instructions,
                         Hybrid<int> *atr_instruction_groups, Hybrid<int4> *ata_instructions,
                         Hybrid<int> *ata_instruction_groups, int *atr_instruction_count,
                         int *ata_instruction_count, const GpuDetails &gpu);

  /// \brief Call the generateWorkUnits() member function to interpret the system indexing and
  ///        produce work units for each grouping of structures.
  ///
  /// \param gpu  Specifications of the GPU available
  void setWorkUnits(const GpuDetails &gpu);

  /// \brief Validate the plan index.
  ///
  /// Overloaded:
  ///   - Provide the partition index and manner of grouping systems
  ///   - Provide the partition index, member within the parition, and manner of grouping systems
  ///
  /// \param index         Index of the plan that will be requested
  /// \param caller        Name of the calling function
  /// \param organization  The manner of partitioning systems within the synthesis
  /// \{
  void validatePartitionIndex(int index, int member, SystemGrouping organization,
                              const char* caller = nullptr) const;

  void validatePartitionIndex(int index, SystemGrouping organization,
                              const char* caller = nullptr) const;
  /// \}
};

} // namespace analysis
} // namespace stormm

#include "comparison_guide.tpp"

#endif
