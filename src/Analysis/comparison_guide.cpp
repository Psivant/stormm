#include "copyright.h"
#include "Math/rounding.h"
#include "comparison_guide.h"

namespace stormm {
namespace analysis {

using stmath::roundUp;
using synthesis::CondensateReader;
using synthesis::SynthesisMapReader;
using synthesis::PsSynthesisReader;
  
//-------------------------------------------------------------------------------------------------
CompGuideKit::CompGuideKit(const int system_count_in, const int source_count_in,
                           const int topology_count_in, const int label_count_in,
                           const int natr_insr_src_in, const int natr_insr_top_in,
                           const int natr_insr_lbl_in, const int nata_insr_src_in,
                           const int nata_insr_top_in, const int nata_insr_lbl_in,
                           const int4* atr_member_src_in, const int4* atr_member_top_in,
                           const int4* atr_member_lbl_in, const int* atr_group_src_in,
                           const int* atr_group_top_in, const int* atr_group_lbl_in,
                           const int4* ata_member_src_in, const int4* ata_member_top_in,
                           const int4* ata_member_lbl_in, const int* ata_group_src_in,
                           const int* ata_group_top_in, const int* ata_group_lbl_in,
                           const int* src_groups_in, const int* top_groups_in,
                           const int* lbl_groups_in, const int* src_grp_ag_in,
                           const int* top_grp_ag_in, const int* lbl_grp_ag_in,
                           const int* src_grp_bounds_in, const int* top_grp_bounds_in,
                           const int* lbl_grp_bounds_in, const size_t* se_pair_src_offsets_in,
                           const size_t* se_pair_top_offsets_in,
                           const size_t* se_pair_lbl_offsets_in,
                           const size_t* nr_pair_src_offsets_in,
                           const size_t* nr_pair_top_offsets_in,
                           const size_t* nr_pair_lbl_offsets_in) :
    system_count{system_count_in}, source_count{source_count_in},
    topology_count{topology_count_in}, label_count{label_count_in},
    natr_insr_src{natr_insr_src_in}, natr_insr_top{natr_insr_top_in},
    natr_insr_lbl{natr_insr_lbl_in}, nata_insr_src{nata_insr_src_in},
    nata_insr_top{nata_insr_top_in}, nata_insr_lbl{nata_insr_lbl_in},
    atr_member_src{atr_member_src_in}, atr_member_top{atr_member_top_in},
    atr_member_lbl{atr_member_lbl_in}, atr_group_src{atr_group_src_in},
    atr_group_top{atr_group_top_in}, atr_group_lbl{atr_group_lbl_in},
    ata_member_src{atr_member_src_in}, ata_member_top{atr_member_top_in},
    ata_member_lbl{atr_member_lbl_in}, ata_group_src{atr_group_src_in},
    ata_group_top{atr_group_top_in}, ata_group_lbl{atr_group_lbl_in},
    src_groups{src_groups_in}, top_groups{top_groups_in}, lbl_groups{lbl_groups_in},
    src_grp_ag{src_grp_ag_in}, top_grp_ag{top_grp_ag_in}, lbl_grp_ag{lbl_grp_ag_in},
    src_grp_bounds{src_grp_bounds_in}, top_grp_bounds{top_grp_bounds_in},
    lbl_grp_bounds{lbl_grp_bounds_in}, se_pair_src_offsets{se_pair_src_offsets_in},
    se_pair_top_offsets{se_pair_top_offsets_in}, se_pair_lbl_offsets{se_pair_lbl_offsets_in},
    nr_pair_src_offsets{nr_pair_src_offsets_in}, nr_pair_top_offsets{nr_pair_top_offsets_in},
    nr_pair_lbl_offsets{nr_pair_lbl_offsets_in}
{}

//-------------------------------------------------------------------------------------------------
ComparisonGuide::ComparisonGuide() :
    system_count{0}, basis{StructureSource::NONE}, csptr_data_type{int_type_index},
    atr_instruction_count_src{0}, atr_instruction_count_top{0}, atr_instruction_count_lbl{0},
    ata_instruction_count_src{0}, ata_instruction_count_top{0}, ata_instruction_count_lbl{0},
    atr_instruction_members_src{HybridKind::ARRAY, "cguide_atr_mmbr_src"},
    atr_instruction_members_top{HybridKind::ARRAY, "cguide_atr_mmbr_top"},
    atr_instruction_members_lbl{HybridKind::ARRAY, "cguide_atr_mmbr_lbl"},
    atr_instruction_groups_src{HybridKind::ARRAY, "cguide_atr_grp_src"},
    atr_instruction_groups_top{HybridKind::ARRAY, "cguide_atr_grp_top"},
    atr_instruction_groups_lbl{HybridKind::ARRAY, "cguide_atr_grp_lbl"},
    ata_instruction_members_src{HybridKind::ARRAY, "cguide_ata_mmbr_src"},
    ata_instruction_members_top{HybridKind::ARRAY, "cguide_ata_mmbr_top"},
    ata_instruction_members_lbl{HybridKind::ARRAY, "cguide_ata_mmbr_lbl"},
    ata_instruction_groups_src{HybridKind::ARRAY, "cguide_ata_grp_src"},
    ata_instruction_groups_top{HybridKind::ARRAY, "cguide_ata_grp_top"},
    ata_instruction_groups_lbl{HybridKind::ARRAY, "cguide_ata_grp_lbl"},
    cache_source_groups{HybridKind::POINTER, "cguide_src_idx"},
    topology_groups{HybridKind::POINTER, "cguide_top_idx"},
    cache_label_groups{HybridKind::POINTER, "cguide_lbl_idx"},
    cache_source_group_topologies{HybridKind::POINTER, "cguide_src_top_idx"},
    topology_group_topologies{HybridKind::POINTER, "cguide_top_top_idx"},
    cache_label_group_topologies{HybridKind::POINTER, "cguide_lbl_top_idx"},
    cache_source_group_bounds{HybridKind::POINTER, "cguide_src_bounds"},
    topology_group_bounds{HybridKind::POINTER, "cguide_top_bounds"},
    cache_label_group_bounds{HybridKind::POINTER, "cguide_lbl_bounds"},
    int_data{HybridKind::ARRAY, "cguide_ints"},
    sepair_result_alloc_src{0}, sepair_result_alloc_top{0}, sepair_result_alloc_lbl{0},
    nrpair_result_alloc_src{0}, nrpair_result_alloc_top{0}, nrpair_result_alloc_lbl{0},
    sepair_result_offsets_src{HybridKind::POINTER, "cguide_seres_src"},
    sepair_result_offsets_top{HybridKind::POINTER, "cguide_seres_top"},
    sepair_result_offsets_lbl{HybridKind::POINTER, "cguide_seres_lbl"},
    nrpair_result_offsets_src{HybridKind::POINTER, "cguide_nrres_src"},
    nrpair_result_offsets_top{HybridKind::POINTER, "cguide_nrres_top"},
    nrpair_result_offsets_lbl{HybridKind::POINTER, "cguide_nrres_lbl"},
    sizet_data{HybridKind::ARRAY, "cguide_sizets"},
    pps_ptr{nullptr}, cdns_ptr{nullptr}, scmap_ptr{nullptr}, cs_ptr{nullptr},
    ag_pointers{}
{}

//-------------------------------------------------------------------------------------------------
ComparisonGuide::ComparisonGuide(const PhaseSpaceSynthesis *poly_ps_in, const GpuDetails &gpu) :
    ComparisonGuide()
{
  pps_ptr = const_cast<PhaseSpaceSynthesis*>(poly_ps_in);
  basis = StructureSource::SYNTHESIS;
  allocateSystemIndexing();
  setSystemIndexing();
  setWorkUnits(gpu);
}

//-------------------------------------------------------------------------------------------------
ComparisonGuide::ComparisonGuide(const PhaseSpaceSynthesis &poly_ps_in, const GpuDetails &gpu) :
    ComparisonGuide(poly_ps_in.getSelfPointer(), gpu)
{}

//-------------------------------------------------------------------------------------------------
ComparisonGuide::ComparisonGuide(const PhaseSpaceSynthesis *poly_ps_in,
                                 const SynthesisCacheMap *scmap_in, const GpuDetails &gpu) :
    ComparisonGuide()
{
  pps_ptr = const_cast<PhaseSpaceSynthesis*>(poly_ps_in);
  scmap_ptr = const_cast<SynthesisCacheMap*>(scmap_in);
  if (scmap_ptr->getCoordinateSynthesisPointer() != pps_ptr) {
    rtErr("The synthesis cache map must reference a synthesis identical to the one supplied to "
          "the guide.", "ComparisonGuide");
  }
  basis = StructureSource::SYNTHESIS;
  allocateSystemIndexing();
  setSystemIndexing();
  setWorkUnits(gpu);
}

//-------------------------------------------------------------------------------------------------
ComparisonGuide::ComparisonGuide(const PhaseSpaceSynthesis &poly_ps_in,
                                 const SynthesisCacheMap &scmap_in, const GpuDetails &gpu) :
    ComparisonGuide(poly_ps_in.getSelfPointer(), scmap_in.getSelfPointer(), gpu)
{}

//-------------------------------------------------------------------------------------------------
ComparisonGuide::ComparisonGuide(const Condensate *cdns_in, const GpuDetails &gpu) :
    ComparisonGuide()
{
  cdns_ptr = const_cast<Condensate*>(cdns_in);
  basis = cdns_ptr->getBasis();
  if (basis != StructureSource::SYNTHESIS) {
    rtErr("A condensate must be based on a PhaseSpaceSynthesis in order to construct a "
          "ComparisonGuide without also supplying a topology governing the systems in it.",
          "ComparisonGuide");
  }
  allocateSystemIndexing();
  setSystemIndexing();
  setWorkUnits(gpu);
}
  
//-------------------------------------------------------------------------------------------------
ComparisonGuide::ComparisonGuide(const Condensate &cdns_in, const GpuDetails &gpu) :
    ComparisonGuide(cdns_in.getSelfPointer(), gpu)
{}

//-------------------------------------------------------------------------------------------------
ComparisonGuide::ComparisonGuide(const Condensate *cdns_in, const AtomGraph *ag_in,
                                 const GpuDetails &gpu) :
    ComparisonGuide()
{
  cdns_ptr = const_cast<Condensate*>(cdns_in);
  basis = cdns_ptr->getBasis();
  ag_pointers.resize(1);
  ag_pointers[0] = const_cast<AtomGraph*>(ag_in);

  // The Condensate must be based on a series if supplying an auxiliary topology, which will be
  // assumed to describe the systems in the series.
  if (basis != StructureSource::SERIES) {
    rtErr("A ComparisonGuide based on a Condensate with an auxiliary topology implies that the "
          "Condensate is itself based on a coordinate series, but the basis of the Condesate is " +
          getEnumerationName(cdns_in->getBasis()) + ".", "ComparisonGuide");
  }
  const CondensateReader cdnsr = cdns_in->data();
  if (cdnsr.system_count == 0) {
    rtErr("No coordinates are present in the Condensate.", "ComparisonGuide");
  }
  if (ag_pointers[0]->getAtomCount() != cdnsr.atom_counts[0]) {
    rtErr("The number of atoms in the governing topology (" +
          std::to_string(ag_pointers[0]->getAtomCount()) + ") must equal the number of atoms in "
          "each frame of the series (" + std::to_string(cdnsr.atom_counts[0]) + ").",
          "ComparisonGuide");
  }

  // If the Condensate source is based on a CoordinateSeries and does NOT allocate memory for its
  // own coordinates, there is a CoordinateSeries object available which contains all of the
  // original data.  Transfer the details of that series to this object and mark it as based on the
  // series, not the Condensate.  This pass-through behavior only applies if the Condensate is
  // provided with no additional information in the form of a cache map.
  if (cdns_in->ownsCoordinates() == false) {
    cs_ptr = const_cast<CoordinateSeries<int>*>(cdns_in->getSeriesPointer<int>());
    csptr_data_type = cdns_in->getCoordinateSeriesTypeID();
    basis = StructureSource::SERIES;
    cdns_ptr = nullptr;
  }
  allocateSystemIndexing();
  setSystemIndexing();
  setWorkUnits(gpu);
}

//-------------------------------------------------------------------------------------------------
ComparisonGuide::ComparisonGuide(const Condensate &cdns_in, const AtomGraph &ag_in,
                                 const GpuDetails &gpu) :
    ComparisonGuide(cdns_in.getSelfPointer(), ag_in.getSelfPointer(), gpu)
{}

//-------------------------------------------------------------------------------------------------
ComparisonGuide::ComparisonGuide(const Condensate *cdns_in, const SynthesisCacheMap *scmap_in,
                                 const GpuDetails &gpu) :
    ComparisonGuide()
{
  cdns_ptr = const_cast<Condensate*>(cdns_in);
  scmap_ptr = const_cast<SynthesisCacheMap*>(scmap_in);
  if (scmap_ptr->getCoordinateSynthesisPointer() != cdns_ptr->getSynthesisPointer()) {
    rtErr("The synthesis cache map must reference a synthesis identical to the one referenced by "
          "the condensate object.", "ComparisonGuide");
  }
  basis = StructureSource::SYNTHESIS;
  allocateSystemIndexing();
  setSystemIndexing();
  setWorkUnits(gpu);
}

//-------------------------------------------------------------------------------------------------
ComparisonGuide::ComparisonGuide(const Condensate &cdns_in, const SynthesisCacheMap &scmap_in,
                                 const GpuDetails &gpu) :
    ComparisonGuide(cdns_in.getSelfPointer(), scmap_in.getSelfPointer(), gpu)
{}

//-------------------------------------------------------------------------------------------------
ComparisonGuide::ComparisonGuide(const ComparisonGuide &original) :
    system_count{original.system_count},
    source_count{original.source_count},
    topology_count{original.topology_count},
    label_count{original.label_count},
    basis{original.basis},
    csptr_data_type{original.csptr_data_type},
    atr_instruction_count_src{original.atr_instruction_count_src},
    atr_instruction_count_top{original.atr_instruction_count_top},
    atr_instruction_count_lbl{original.atr_instruction_count_lbl},
    ata_instruction_count_src{original.ata_instruction_count_src},
    ata_instruction_count_top{original.ata_instruction_count_top},
    ata_instruction_count_lbl{original.ata_instruction_count_lbl},
    atr_instruction_members_src{original.atr_instruction_members_src},
    atr_instruction_members_top{original.atr_instruction_members_top},
    atr_instruction_members_lbl{original.atr_instruction_members_lbl},
    atr_instruction_groups_src{original.atr_instruction_groups_src},
    atr_instruction_groups_top{original.atr_instruction_groups_top},
    atr_instruction_groups_lbl{original.atr_instruction_groups_lbl},
    ata_instruction_members_src{original.ata_instruction_members_src},
    ata_instruction_members_top{original.ata_instruction_members_top},
    ata_instruction_members_lbl{original.ata_instruction_members_lbl},
    ata_instruction_groups_src{original.ata_instruction_groups_src},
    ata_instruction_groups_top{original.ata_instruction_groups_top},
    ata_instruction_groups_lbl{original.ata_instruction_groups_lbl},
    cache_source_groups{original.cache_source_groups},
    topology_groups{original.topology_groups},
    cache_label_groups{original.cache_label_groups},
    cache_source_group_topologies{original.cache_source_group_topologies},
    topology_group_topologies{original.topology_group_topologies},
    cache_label_group_topologies{original.cache_label_group_topologies},
    cache_source_group_bounds{original.cache_source_group_bounds},
    topology_group_bounds{original.topology_group_bounds},
    cache_label_group_bounds{original.cache_label_group_bounds},
    int_data{original.int_data},
    sepair_result_alloc_src{original.sepair_result_alloc_src},
    sepair_result_alloc_top{original.sepair_result_alloc_top},
    sepair_result_alloc_lbl{original.sepair_result_alloc_lbl},
    nrpair_result_alloc_src{original.nrpair_result_alloc_src},
    nrpair_result_alloc_top{original.nrpair_result_alloc_top},
    nrpair_result_alloc_lbl{original.nrpair_result_alloc_lbl},
    sepair_result_offsets_src{original.sepair_result_offsets_src},
    sepair_result_offsets_top{original.sepair_result_offsets_top},
    sepair_result_offsets_lbl{original.sepair_result_offsets_lbl},
    nrpair_result_offsets_src{original.nrpair_result_offsets_src},
    nrpair_result_offsets_top{original.nrpair_result_offsets_top},
    nrpair_result_offsets_lbl{original.nrpair_result_offsets_lbl},
    sizet_data{original.sizet_data},
    pps_ptr{original.pps_ptr},
    cdns_ptr{original.cdns_ptr},
    scmap_ptr{original.scmap_ptr},
    cs_ptr{original.cs_ptr}
{
  // Repair POINTER-kind Hybrid objects
  allocateSystemIndexing();
}

//-------------------------------------------------------------------------------------------------
ComparisonGuide::ComparisonGuide(ComparisonGuide &&original) :
    system_count{original.system_count},
    basis{original.basis},
    source_count{original.source_count},
    topology_count{original.topology_count},
    label_count{original.label_count},
    csptr_data_type{original.csptr_data_type},
    atr_instruction_count_src{original.atr_instruction_count_src},
    atr_instruction_count_top{original.atr_instruction_count_top},
    atr_instruction_count_lbl{original.atr_instruction_count_lbl},
    ata_instruction_count_src{original.ata_instruction_count_src},
    ata_instruction_count_top{original.ata_instruction_count_top},
    ata_instruction_count_lbl{original.ata_instruction_count_lbl},
    atr_instruction_members_src{std::move(original.atr_instruction_members_src)},
    atr_instruction_members_top{std::move(original.atr_instruction_members_top)},
    atr_instruction_members_lbl{std::move(original.atr_instruction_members_lbl)},
    atr_instruction_groups_src{std::move(original.atr_instruction_groups_src)},
    atr_instruction_groups_top{std::move(original.atr_instruction_groups_top)},
    atr_instruction_groups_lbl{std::move(original.atr_instruction_groups_lbl)},
    ata_instruction_members_src{std::move(original.ata_instruction_members_src)},
    ata_instruction_members_top{std::move(original.ata_instruction_members_top)},
    ata_instruction_members_lbl{std::move(original.ata_instruction_members_lbl)},
    ata_instruction_groups_src{std::move(original.ata_instruction_groups_src)},
    ata_instruction_groups_top{std::move(original.ata_instruction_groups_top)},
    ata_instruction_groups_lbl{std::move(original.ata_instruction_groups_lbl)},
    cache_source_groups{std::move(original.cache_source_groups)},
    topology_groups{std::move(original.topology_groups)},
    cache_label_groups{std::move(original.cache_label_groups)},
    cache_source_group_topologies{std::move(original.cache_source_group_topologies)},
    topology_group_topologies{std::move(original.topology_group_topologies)},
    cache_label_group_topologies{std::move(original.cache_label_group_topologies)},
    cache_source_group_bounds{std::move(original.cache_source_group_bounds)},
    topology_group_bounds{std::move(original.topology_group_bounds)},
    cache_label_group_bounds{std::move(original.cache_label_group_bounds)},
    int_data{std::move(original.int_data)},
    sepair_result_alloc_src{original.sepair_result_alloc_src},
    sepair_result_alloc_top{original.sepair_result_alloc_top},
    sepair_result_alloc_lbl{original.sepair_result_alloc_lbl},
    nrpair_result_alloc_src{original.nrpair_result_alloc_src},
    nrpair_result_alloc_top{original.nrpair_result_alloc_top},
    nrpair_result_alloc_lbl{original.nrpair_result_alloc_lbl},
    sepair_result_offsets_src{std::move(original.sepair_result_offsets_src)},
    sepair_result_offsets_top{std::move(original.sepair_result_offsets_top)},
    sepair_result_offsets_lbl{std::move(original.sepair_result_offsets_lbl)},
    nrpair_result_offsets_src{std::move(original.nrpair_result_offsets_src)},
    nrpair_result_offsets_top{std::move(original.nrpair_result_offsets_top)},
    nrpair_result_offsets_lbl{std::move(original.nrpair_result_offsets_lbl)},
    sizet_data{std::move(original.sizet_data)},
    pps_ptr{original.pps_ptr},
    cdns_ptr{original.cdns_ptr},
    scmap_ptr{original.scmap_ptr},
    cs_ptr{original.cs_ptr}
{}

//-------------------------------------------------------------------------------------------------
ComparisonGuide& ComparisonGuide::operator=(const ComparisonGuide &other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // Copy all elements
  system_count = other.system_count;
  source_count = other.source_count;
  topology_count = other.topology_count;
  label_count = other.label_count;
  basis = other.basis;
  csptr_data_type = other.csptr_data_type;
  atr_instruction_count_src = other.atr_instruction_count_src;
  atr_instruction_count_top = other.atr_instruction_count_top;
  atr_instruction_count_lbl = other.atr_instruction_count_lbl;
  ata_instruction_count_src = other.ata_instruction_count_src;
  ata_instruction_count_top = other.ata_instruction_count_top;
  ata_instruction_count_lbl = other.ata_instruction_count_lbl;
  atr_instruction_members_src = other.atr_instruction_members_src;
  atr_instruction_members_top = other.atr_instruction_members_top;
  atr_instruction_members_lbl = other.atr_instruction_members_lbl;
  atr_instruction_groups_src = other.atr_instruction_groups_src;
  atr_instruction_groups_top = other.atr_instruction_groups_top;
  atr_instruction_groups_lbl = other.atr_instruction_groups_lbl;
  ata_instruction_members_src = other.ata_instruction_members_src;
  ata_instruction_members_top = other.ata_instruction_members_top;
  ata_instruction_members_lbl = other.ata_instruction_members_lbl;
  ata_instruction_groups_src = other.ata_instruction_groups_src;
  ata_instruction_groups_top = other.ata_instruction_groups_top;
  ata_instruction_groups_lbl = other.ata_instruction_groups_lbl;
  cache_source_groups = other.cache_source_groups;
  topology_groups = other.topology_groups;
  cache_label_groups = other.cache_label_groups;
  cache_source_group_topologies = other.cache_source_group_topologies;
  topology_group_topologies = other.topology_group_topologies;
  cache_label_group_topologies = other.cache_label_group_topologies;
  cache_source_group_bounds = other.cache_source_group_bounds;
  topology_group_bounds = other.topology_group_bounds;
  cache_label_group_bounds = other.cache_label_group_bounds;
  int_data = other.int_data;
  sepair_result_alloc_src = other.sepair_result_alloc_src;
  sepair_result_alloc_top = other.sepair_result_alloc_top;
  sepair_result_alloc_lbl = other.sepair_result_alloc_lbl;
  nrpair_result_alloc_src = other.nrpair_result_alloc_src;
  nrpair_result_alloc_top = other.nrpair_result_alloc_top;
  nrpair_result_alloc_lbl = other.nrpair_result_alloc_lbl;
  sepair_result_offsets_src = other.sepair_result_offsets_src;
  sepair_result_offsets_top = other.sepair_result_offsets_top;
  sepair_result_offsets_lbl = other.sepair_result_offsets_lbl;
  nrpair_result_offsets_src = other.nrpair_result_offsets_src;
  nrpair_result_offsets_top = other.nrpair_result_offsets_top;
  nrpair_result_offsets_lbl = other.nrpair_result_offsets_lbl;
  sizet_data = other.sizet_data;
  pps_ptr = other.pps_ptr;
  cdns_ptr = other.cdns_ptr;
  scmap_ptr = other.scmap_ptr;
  cs_ptr = other.cs_ptr;

  // Repair POINTER-kind Hybrid objects and return the result
  allocateSystemIndexing();
  return *this;
}

//-------------------------------------------------------------------------------------------------
ComparisonGuide& ComparisonGuide::operator=(ComparisonGuide &&other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }

  // Move arrays or copy various scalar / pointer elements
  system_count = other.system_count;
  source_count = other.source_count;
  topology_count = other.topology_count;
  label_count = other.label_count;
  basis = other.basis;
  csptr_data_type = other.csptr_data_type;
  atr_instruction_count_src = other.atr_instruction_count_src;
  atr_instruction_count_top = other.atr_instruction_count_top;
  atr_instruction_count_lbl = other.atr_instruction_count_lbl;
  ata_instruction_count_src = other.ata_instruction_count_src;
  ata_instruction_count_top = other.ata_instruction_count_top;
  ata_instruction_count_lbl = other.ata_instruction_count_lbl;
  atr_instruction_members_src= std::move(other.atr_instruction_members_src);
  atr_instruction_members_top= std::move(other.atr_instruction_members_top);
  atr_instruction_members_lbl= std::move(other.atr_instruction_members_lbl);
  atr_instruction_groups_src= std::move(other.atr_instruction_groups_src);
  atr_instruction_groups_top= std::move(other.atr_instruction_groups_top);
  atr_instruction_groups_lbl= std::move(other.atr_instruction_groups_lbl);
  ata_instruction_members_src= std::move(other.ata_instruction_members_src);
  ata_instruction_members_top= std::move(other.ata_instruction_members_top);
  ata_instruction_members_lbl= std::move(other.ata_instruction_members_lbl);
  ata_instruction_groups_src= std::move(other.ata_instruction_groups_src);
  ata_instruction_groups_top= std::move(other.ata_instruction_groups_top);
  ata_instruction_groups_lbl= std::move(other.ata_instruction_groups_lbl);
  cache_source_groups= std::move(other.cache_source_groups);
  topology_groups= std::move(other.topology_groups);
  cache_label_groups= std::move(other.cache_label_groups);
  cache_source_group_topologies= std::move(other.cache_source_group_topologies);
  topology_group_topologies= std::move(other.topology_group_topologies);
  cache_label_group_topologies= std::move(other.cache_label_group_topologies);
  cache_source_group_bounds= std::move(other.cache_source_group_bounds);
  topology_group_bounds= std::move(other.topology_group_bounds);
  cache_label_group_bounds= std::move(other.cache_label_group_bounds);
  int_data= std::move(other.int_data);
  sepair_result_alloc_src = other.sepair_result_alloc_src;
  sepair_result_alloc_top = other.sepair_result_alloc_top;
  sepair_result_alloc_lbl = other.sepair_result_alloc_lbl;
  nrpair_result_alloc_src = other.nrpair_result_alloc_src;
  nrpair_result_alloc_top = other.nrpair_result_alloc_top;
  nrpair_result_alloc_lbl = other.nrpair_result_alloc_lbl;
  sepair_result_offsets_src= std::move(other.sepair_result_offsets_src);
  sepair_result_offsets_top= std::move(other.sepair_result_offsets_top);
  sepair_result_offsets_lbl= std::move(other.sepair_result_offsets_lbl);
  nrpair_result_offsets_src= std::move(other.nrpair_result_offsets_src);
  nrpair_result_offsets_top= std::move(other.nrpair_result_offsets_top);
  nrpair_result_offsets_lbl= std::move(other.nrpair_result_offsets_lbl);
  sizet_data= std::move(other.sizet_data);
  pps_ptr = other.pps_ptr;
  cdns_ptr = other.cdns_ptr;
  scmap_ptr = other.scmap_ptr;
  cs_ptr = other.cs_ptr;
  return *this;
}

//-------------------------------------------------------------------------------------------------
int ComparisonGuide::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
int ComparisonGuide::getTopologyCount() const {
  return ag_pointers.size();
}

//-------------------------------------------------------------------------------------------------
StructureSource ComparisonGuide::getBasis() const {
  return basis;
}

//-------------------------------------------------------------------------------------------------
int ComparisonGuide::getPartitionCount(const SystemGrouping organization) const {
  switch (basis) {
  case StructureSource::SYNTHESIS:
    if (scmap_ptr != nullptr) {
      return scmap_ptr->getPartitionCount(organization);
    }
    else {
      return (pps_ptr != nullptr) ? pps_ptr->getUniqueTopologyCount() :
                                    cdns_ptr->getSynthesisPointer()->getUniqueTopologyCount();
    }
    break;
  case StructureSource::SERIES:
    return 1;
  case StructureSource::NONE:
    return 0;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
size_t ComparisonGuide::getAllToReferenceOutputSize() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
size_t
ComparisonGuide::getSymmetryEquivalentPairOutputSize(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return sepair_result_alloc_src;
  case SystemGrouping::TOPOLOGY:
    return sepair_result_alloc_top;
  case SystemGrouping::LABEL:
    return sepair_result_alloc_lbl;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
size_t ComparisonGuide::getNonReflexivePairOutputSize(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return nrpair_result_alloc_src;
  case SystemGrouping::TOPOLOGY:
    return nrpair_result_alloc_top;
  case SystemGrouping::LABEL:
    return nrpair_result_alloc_lbl;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int ComparisonGuide::getATRInstructionCount(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return atr_instruction_count_src;
  case SystemGrouping::TOPOLOGY:
    return atr_instruction_count_top;
  case SystemGrouping::LABEL:
    return atr_instruction_count_lbl;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int ComparisonGuide::getATAInstructionCount(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return ata_instruction_count_src;
  case SystemGrouping::TOPOLOGY:
    return ata_instruction_count_top;
  case SystemGrouping::LABEL:
    return ata_instruction_count_lbl;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int4>
ComparisonGuide::getATRInstructionMembers(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return atr_instruction_members_src.readHost();
  case SystemGrouping::TOPOLOGY:
    return atr_instruction_members_top.readHost();
  case SystemGrouping::LABEL:
    return atr_instruction_members_lbl.readHost();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ComparisonGuide::getATRInstructionGroups(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return atr_instruction_groups_src.readHost();
  case SystemGrouping::TOPOLOGY:
    return atr_instruction_groups_top.readHost();
  case SystemGrouping::LABEL:
    return atr_instruction_groups_lbl.readHost();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int4>
ComparisonGuide::getATAInstructionMembers(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return ata_instruction_members_src.readHost();
  case SystemGrouping::TOPOLOGY:
    return ata_instruction_members_top.readHost();
  case SystemGrouping::LABEL:
    return ata_instruction_members_lbl.readHost();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ComparisonGuide::getATAInstructionGroups(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return ata_instruction_groups_src.readHost();
  case SystemGrouping::TOPOLOGY:
    return ata_instruction_groups_top.readHost();
  case SystemGrouping::LABEL:
    return ata_instruction_groups_lbl.readHost();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ComparisonGuide::getGroupStructureIndices(const int partition,
                                          const SystemGrouping organization) const {
  validatePartitionIndex(partition, organization, "getGroupStructureIndices");
  switch (organization) {
  case SystemGrouping::SOURCE:
    {
      const int llim = cache_source_group_bounds.readHost(partition);
      const int hlim = cache_source_group_bounds.readHost(partition + 1);
      return cache_source_groups.readHost(llim, hlim - llim);
    }
    break;
  case SystemGrouping::TOPOLOGY:
    {
      const int llim = topology_group_bounds.readHost(partition);
      const int hlim = topology_group_bounds.readHost(partition + 1);
      return topology_groups.readHost(llim, hlim - llim);
    }
    break;
  case SystemGrouping::LABEL:
    {
      const int llim = cache_label_group_bounds.readHost(partition);
      const int hlim = cache_label_group_bounds.readHost(partition + 1);
      return cache_label_groups.readHost(llim, hlim - llim);
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ComparisonGuide::getGroupStructureIndices(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return cache_source_groups.readHost();
  case SystemGrouping::TOPOLOGY:
    return topology_groups.readHost();
  case SystemGrouping::LABEL:
    return cache_label_groups.readHost();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ComparisonGuide::getGroupTopologyIndices(const int partition,
                                          const SystemGrouping organization) const {
  validatePartitionIndex(partition, organization, "getGroupStructureIndices");
  switch (organization) {
  case SystemGrouping::SOURCE:
    {
      const int llim = cache_source_group_bounds.readHost(partition);
      const int hlim = cache_source_group_bounds.readHost(partition + 1);
      return cache_source_group_topologies.readHost(llim, hlim - llim);
    }
    break;
  case SystemGrouping::TOPOLOGY:
    {
      const int llim = topology_group_bounds.readHost(partition);
      const int hlim = topology_group_bounds.readHost(partition + 1);
      return topology_group_topologies.readHost(llim, hlim - llim);
    }
    break;
  case SystemGrouping::LABEL:
    {
      const int llim = cache_label_group_bounds.readHost(partition);
      const int hlim = cache_label_group_bounds.readHost(partition + 1);
      return cache_label_group_topologies.readHost(llim, hlim - llim);
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
ComparisonGuide::getGroupTopologyIndices(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return cache_source_group_topologies.readHost();
  case SystemGrouping::TOPOLOGY:
    return topology_group_topologies.readHost();
  case SystemGrouping::LABEL:
    return cache_label_group_topologies.readHost();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ComparisonGuide::getGroupBounds(const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return cache_source_group_bounds.readHost();
  case SystemGrouping::TOPOLOGY:
    return topology_group_bounds.readHost();
  case SystemGrouping::LABEL:
    return cache_label_group_bounds.readHost();
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int ComparisonGuide::getFrameCount(const int index, const SystemGrouping organization) const {
  switch (organization) {
  case SystemGrouping::SOURCE:
    return cache_source_group_bounds.readHost(index + 1) -
           cache_source_group_bounds.readHost(index);
  case SystemGrouping::TOPOLOGY:
    return topology_group_bounds.readHost(index + 1) - topology_group_bounds.readHost(index);
  case SystemGrouping::LABEL:
    return cache_label_group_bounds.readHost(index + 1) - cache_label_group_bounds.readHost(index);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
size_t ComparisonGuide::getAllToOneResultOffset(const int index,
                                                const SystemGrouping organization) const {
  validatePartitionIndex(index, organization, "getAllToOneResultOffset");
  switch (organization) {
  case SystemGrouping::SOURCE:
    return cache_source_group_bounds.readHost(index);
  case SystemGrouping::TOPOLOGY:
    return topology_group_bounds.readHost(index);
  case SystemGrouping::LABEL:
    return cache_label_group_bounds.readHost(index);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
size_t
ComparisonGuide::getSymmetryEquivalentResultOffset(const int index,
                                                   const SystemGrouping organization) const {
  validatePartitionIndex(index, organization, "getSymmetryEquivalentResultOffset");
  switch (organization) {
  case SystemGrouping::SOURCE:
    return sepair_result_offsets_src.readHost(index);
  case SystemGrouping::TOPOLOGY:
    return sepair_result_offsets_top.readHost(index);
  case SystemGrouping::LABEL:
    return sepair_result_offsets_lbl.readHost(index);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
size_t ComparisonGuide::getNonReflexiveResultOffset(const int index,
                                                    const SystemGrouping organization) const {
  validatePartitionIndex(index, organization, "getNonReflexiveResultOffset");
  switch (organization) {
  case SystemGrouping::SOURCE:
    return nrpair_result_offsets_src.readHost(index);
  case SystemGrouping::TOPOLOGY:
    return nrpair_result_offsets_top.readHost(index);
  case SystemGrouping::LABEL:
    return nrpair_result_offsets_lbl.readHost(index);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceSynthesis* ComparisonGuide::getPhaseSpaceSynthesisPointer() const {
  if (pps_ptr == nullptr) {
    rtErr("No such object is referenced.", "ComparisonGuide", "getPhaseSpaceSynthesisPointer");
  }
  return pps_ptr;
}

//-------------------------------------------------------------------------------------------------
const Condensate* ComparisonGuide::getCondensatePointer() const {
  if (cdns_ptr == nullptr) {
    rtErr("No such object is referenced.", "ComparisonGuide", "getCondensatePointer");
  }
  return cdns_ptr;
}

//-------------------------------------------------------------------------------------------------
const SynthesisCacheMap* ComparisonGuide::getSynthesisMapPointer() const {
  if (scmap_ptr == nullptr) {
    rtErr("No such object is referenced.", "ComparisonGuide", "getSynthesisMapPointer");
  }
  return scmap_ptr;
}

//-------------------------------------------------------------------------------------------------
size_t ComparisonGuide::getCoordinateSeriesTypeID() const {
  return csptr_data_type;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* ComparisonGuide::getTopologyPointer(const int index) const {
  if (index < 0 || index >= ag_pointers.size()) {
    rtErr("Index " + std::to_string(index) + " is invalid for a collection of " +
          std::to_string(ag_pointers.size()) + " unique topologies.", "ComparisonGuide",
          "getTopologyPointer");
  }
  return ag_pointers[index];
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* ComparisonGuide::getTopologyPointer(const SystemGrouping organization,
                                                     const int partition, const int member) const {
  validatePartitionIndex(partition, member, organization, "getTopologyPointer");
  int idx;
  switch (organization) {
  case SystemGrouping::SOURCE:
    idx = cache_source_group_bounds.readHost(partition) + member;
    return ag_pointers[cache_source_group_topologies.readHost(idx)];
  case SystemGrouping::TOPOLOGY:
    idx = topology_group_bounds.readHost(partition) + member;
    return ag_pointers[topology_group_topologies.readHost(idx)];
  case SystemGrouping::LABEL:
    idx = cache_label_group_bounds.readHost(partition) + member;
    return ag_pointers[cache_label_group_topologies.readHost(idx)];
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
CompGuideKit ComparisonGuide::data(const HybridTargetLevel tier) const {
  return CompGuideKit(system_count, source_count, topology_count, label_count,
                      atr_instruction_count_src, atr_instruction_count_top,
                      atr_instruction_count_lbl, ata_instruction_count_src,
                      ata_instruction_count_top, ata_instruction_count_lbl,
                      atr_instruction_members_src.data(tier),
                      atr_instruction_members_top.data(tier),
                      atr_instruction_members_lbl.data(tier),
                      atr_instruction_groups_src.data(tier),
                      atr_instruction_groups_top.data(tier),
                      atr_instruction_groups_lbl.data(tier),
                      ata_instruction_members_src.data(tier),
                      ata_instruction_members_top.data(tier),
                      ata_instruction_members_lbl.data(tier),
                      ata_instruction_groups_src.data(tier),
                      ata_instruction_groups_top.data(tier),
                      ata_instruction_groups_lbl.data(tier), cache_source_groups.data(tier),
                      topology_groups.data(tier), cache_label_groups.data(tier),
                      cache_source_group_topologies.data(tier),
                      topology_group_topologies.data(tier),
                      cache_label_group_topologies.data(tier),
                      cache_source_group_bounds.data(tier), topology_group_bounds.data(tier),
                      cache_label_group_bounds.data(tier), sepair_result_offsets_src.data(tier),
                      sepair_result_offsets_top.data(tier), sepair_result_offsets_lbl.data(tier),
                      nrpair_result_offsets_src.data(tier), nrpair_result_offsets_top.data(tier),
                      nrpair_result_offsets_lbl.data(tier));
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void ComparisonGuide::upload() {
  atr_instruction_members_src.upload();
  atr_instruction_members_top.upload();
  atr_instruction_members_lbl.upload();
  atr_instruction_groups_src.upload();
  atr_instruction_groups_top.upload();
  atr_instruction_groups_lbl.upload();
  ata_instruction_members_src.upload();
  ata_instruction_members_top.upload();
  ata_instruction_members_lbl.upload();
  ata_instruction_groups_src.upload();
  ata_instruction_groups_top.upload();
  ata_instruction_groups_lbl.upload();
  int_data.upload();
  sizet_data.upload();
}

//-------------------------------------------------------------------------------------------------
void ComparisonGuide::download() {
  atr_instruction_members_src.download();
  atr_instruction_members_top.download();
  atr_instruction_members_lbl.download();
  atr_instruction_groups_src.download();
  atr_instruction_groups_top.download();
  atr_instruction_groups_lbl.download();
  ata_instruction_members_src.download();
  ata_instruction_members_top.download();
  ata_instruction_members_lbl.download();
  ata_instruction_groups_src.download();
  ata_instruction_groups_top.download();
  ata_instruction_groups_lbl.download();
  int_data.download();
  sizet_data.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void ComparisonGuide::allocateSystemIndexing() {
  
  // Determine the number of systems
  if (pps_ptr != nullptr) {
    system_count = pps_ptr->getSystemCount();
  }
  else if (cdns_ptr != nullptr) {
    system_count = cdns_ptr->getSystemCount();
  }
  else if (cs_ptr != nullptr) {
    system_count = cs_ptr->getFrameCount();
  }

  // Determine the number of sources, topologies, and labels
  if (scmap_ptr != nullptr) {
    const SynthesisMapReader scmapr = scmap_ptr->data();
    source_count = scmapr.ncache;
    topology_count = scmapr.ntopol;
    label_count = scmapr.nlabel;
  }
  else {
    switch (basis) {
    case StructureSource::SYNTHESIS:
      if (pps_ptr != nullptr) {
        topology_count = pps_ptr->getUniqueTopologyCount();
      }
      else if (cdns_ptr != nullptr) {
        topology_count = cdns_ptr->getSynthesisPointer()->getUniqueTopologyCount();
      }
      break;
    case StructureSource::SERIES:
      topology_count = 1;
      break;
    case StructureSource::NONE:
      topology_count = 0;
      break;
    }
    source_count = topology_count;
    label_count = topology_count;
  }
  
  // Determine the number of structures in each method of grouping
  const size_t padded_nsys = roundUp<size_t>(system_count, warp_size_zu);
  size_t nsrc, ntop, nlbl;
  switch (basis) {
  case StructureSource::SYNTHESIS:
    if (scmap_ptr != nullptr) {
      const SynthesisMapReader scmapr = scmap_ptr->data();
      nsrc = scmapr.ncache;
      ntop = scmapr.ntopol;
      nlbl = scmapr.nlabel;
    }
    else {
      const PsSynthesisReader poly_psw = (pps_ptr != nullptr) ?
                                         pps_ptr->data() : cdns_ptr->getSynthesisPointer()->data();
      ntop = poly_psw.unique_topology_count;
      nsrc = ntop;
      nlbl = ntop;
    }
    break;
  case StructureSource::SERIES:
    {
      ntop = 1;
      nsrc = 1;
      nlbl = 1;
    }
    break;
  case StructureSource::NONE:
    break;
  }

  // Allocate (resize) ARRAY-kind Hybrid storage and set POINTER-kind Hybrid objects for system
  // indices.
  const size_t padded_nsrc_x = roundUp<size_t>(nsrc + 1, warp_size_zu);
  const size_t padded_ntop_x = roundUp<size_t>(ntop + 1, warp_size_zu);
  const size_t padded_nlbl_x = roundUp<size_t>(nlbl + 1, warp_size_zu);
  const size_t padded_nsrc = roundUp(nsrc, warp_size_zu);
  const size_t padded_ntop = roundUp(ntop, warp_size_zu);
  const size_t padded_nlbl = roundUp(nlbl, warp_size_zu);
  int_data.resize((6LLU * padded_nsys) + padded_nsrc_x + padded_ntop_x + padded_nlbl_x);
  cache_source_groups.setPointer(&int_data,                             0, system_count);
  topology_groups.setPointer(&int_data,                       padded_nsys, system_count);
  cache_label_groups.setPointer(&int_data,             2LLU * padded_nsys, system_count);
  cache_source_group_topologies.setPointer(&int_data,  3LLU * padded_nsys, system_count);
  topology_group_topologies.setPointer(&int_data,      4LLU * padded_nsys, system_count);
  cache_label_group_topologies.setPointer(&int_data,   5LLU * padded_nsys, system_count);
  size_t ic = 6LLU * padded_nsys;
  cache_source_group_bounds.setPointer(&int_data, ic, nsrc + 1);
  ic += padded_nsrc_x;
  topology_group_bounds.setPointer(&int_data,     ic, ntop + 1);
  ic += padded_ntop_x;
  cache_label_group_bounds.setPointer(&int_data,  ic, nlbl + 1);

  // Allocate arrays for results bounds.  The exact bounds will be calculated in the following
  // section, but the space requirements are already known.
  sizet_data.resize(2LLU * (padded_nsrc + padded_ntop + padded_nlbl));
  ic = 0;
  sepair_result_offsets_src.setPointer(&sizet_data, ic, nsrc);
  ic += padded_nsrc;
  sepair_result_offsets_top.setPointer(&sizet_data, ic, ntop);
  ic += padded_ntop;
  sepair_result_offsets_lbl.setPointer(&sizet_data, ic, nlbl);
  ic += padded_nlbl;
  nrpair_result_offsets_src.setPointer(&sizet_data, ic, nsrc);
  ic += padded_nsrc;
  nrpair_result_offsets_top.setPointer(&sizet_data, ic, ntop);
  ic += padded_ntop;
  nrpair_result_offsets_lbl.setPointer(&sizet_data, ic, nlbl);
}

//-------------------------------------------------------------------------------------------------
void ComparisonGuide::setSystemIndexing() {

  // Fill out the structure index arrays.
  int* src_grp_ptr = cache_source_groups.data();
  int* top_grp_ptr = topology_groups.data();
  int* lbl_grp_ptr = cache_label_groups.data();
  int* src_top_ptr = cache_source_group_topologies.data();
  int* top_top_ptr = topology_group_topologies.data();
  int* lbl_top_ptr = cache_label_group_topologies.data();
  int* src_grp_bounds_ptr = cache_source_group_bounds.data();
  int* top_grp_bounds_ptr = topology_group_bounds.data();
  int* lbl_grp_bounds_ptr = cache_label_group_bounds.data();  
  switch (basis) {
  case StructureSource::SYNTHESIS:
    {
      // The guide is based on a synthesis.  Is there a cache map available?
      const PsSynthesisReader poly_psw = (pps_ptr != nullptr) ?
                                         pps_ptr->data() : cdns_ptr->getSynthesisPointer()->data();
      if (scmap_ptr != nullptr) {
        const SynthesisMapReader scmapr = scmap_ptr->data();
        for (int i = 0; i < poly_psw.system_count; i++) {
          src_grp_ptr[i] = scmapr.csystem_proj[i];
          top_grp_ptr[i] = scmapr.ctopol_proj[i];
          lbl_grp_ptr[i] = scmapr.clabel_proj[i];
        }
        for (int i = 0; i < scmapr.ncache + 1; i++) {
          src_grp_bounds_ptr[i] = scmapr.csystem_bounds[i];
        }
        for (int i = 0; i < scmapr.ntopol + 1; i++) {
          top_grp_bounds_ptr[i] = scmapr.ctopol_bounds[i];
        }
        for (int i = 0; i < scmapr.nlabel + 1; i++) {
          lbl_grp_bounds_ptr[i] = scmapr.clabel_bounds[i];
        }
        ag_pointers = scmap_ptr->getCachePointer()->getTopologyPointer();
        for (int i = 0; i < poly_psw.system_count; i++) {
          src_top_ptr[i] = scmapr.topology_origins[src_grp_ptr[i]];
          top_top_ptr[i] = scmapr.topology_origins[top_grp_ptr[i]];
          lbl_top_ptr[i] = scmapr.topology_origins[lbl_grp_ptr[i]];
        }
      }
      else {

        // In the absence of any map to a user-generated cache of systems, assume that each unique
        // topology in use by the synthesis originates in a single entry and that every entry has
        // its own label group.
        for (int i = 0; i < poly_psw.system_count; i++) {
          src_grp_ptr[i] = poly_psw.common_ag_list[i];
          top_grp_ptr[i] = poly_psw.common_ag_list[i];
          lbl_grp_ptr[i] = poly_psw.common_ag_list[i];
        }
        for (int i = 0; i < poly_psw.unique_topology_count + 1; i++) {
          src_grp_bounds_ptr[i] = poly_psw.common_ag_bounds[i];
          top_grp_bounds_ptr[i] = poly_psw.common_ag_bounds[i];
          lbl_grp_bounds_ptr[i] = poly_psw.common_ag_bounds[i];
        }
        ag_pointers.resize(poly_psw.unique_topology_count);
        const std::vector<int> un_ex = pps_ptr->getUniqueTopologyExampleIndices();
        for (int i = 0; i < poly_psw.unique_topology_count; i++) {
          ag_pointers[i] = const_cast<AtomGraph*>(pps_ptr->getSystemTopologyPointer(un_ex[i]));
        }
        for (int i = 0; i < poly_psw.system_count; i++) {
          src_top_ptr[i] = poly_psw.unique_ag_idx[src_grp_ptr[i]];
          top_top_ptr[i] = poly_psw.unique_ag_idx[top_grp_ptr[i]];
          lbl_top_ptr[i] = poly_psw.unique_ag_idx[lbl_grp_ptr[i]];
        }
      }
    }
    break;
  case StructureSource::SERIES:

    // The guide is based on a series.
    for (int i = 0; i < system_count; i++) {
      src_grp_ptr[i] = 0;
      top_grp_ptr[i] = 0;
      lbl_grp_ptr[i] = 0;
      src_top_ptr[i] = 0;
      top_top_ptr[i] = 0;
      lbl_top_ptr[i] = 0;      
    }
    src_grp_bounds_ptr[0] = 0;
    top_grp_bounds_ptr[1] = system_count;
    lbl_grp_bounds_ptr[0] = 0;
    src_grp_bounds_ptr[1] = system_count;
    top_grp_bounds_ptr[0] = 0;
    lbl_grp_bounds_ptr[1] = system_count;

    // If the guide is based on a series, the topology pointers array will have been resized in
    // the constructor with its first and only element set to the sole topology.
    break;
  case StructureSource::NONE:
    break;
  }

  // Fill out the results offset arrays.  Use the array sizes determined in a preceding call to
  // allocateSystemIndexing() to avoid repeating the deduction based on the coordinate source.
  const int nsrc = cache_source_group_bounds.size() - 1;
  const int ntop = topology_group_bounds.size() - 1;
  const int nlbl = cache_label_group_bounds.size() - 1;
  const int* src_lims = cache_source_group_bounds.data();
  const int* top_lims = topology_group_bounds.data();
  const int* lbl_lims = cache_label_group_bounds.data();
  size_t* se_src_ptr = sepair_result_offsets_src.data();
  size_t* se_top_ptr = sepair_result_offsets_top.data();
  size_t* se_lbl_ptr = sepair_result_offsets_lbl.data();
  size_t* nr_src_ptr = nrpair_result_offsets_src.data();
  size_t* nr_top_ptr = nrpair_result_offsets_top.data();
  size_t* nr_lbl_ptr = nrpair_result_offsets_lbl.data();
  size_t se_src_result_offset = 0;
  size_t se_top_result_offset = 0;
  size_t se_lbl_result_offset = 0;
  size_t nr_src_result_offset = 0;
  size_t nr_top_result_offset = 0;
  size_t nr_lbl_result_offset = 0;
  for (int i = 0; i < nsrc; i++) {
    const size_t ni_src = src_lims[i + 1] - src_lims[i];
    const size_t ni_top = top_lims[i + 1] - top_lims[i];
    const size_t ni_lbl = lbl_lims[i + 1] - lbl_lims[i];
    se_src_ptr[i] = se_src_result_offset;
    se_top_ptr[i] = se_top_result_offset;
    se_lbl_ptr[i] = se_lbl_result_offset;
    nr_src_ptr[i] = nr_src_result_offset;
    nr_top_ptr[i] = nr_top_result_offset;
    nr_lbl_ptr[i] = nr_lbl_result_offset;
    se_src_result_offset += roundUp<size_t>((ni_src - 1LLU) * ni_src / 2LLU, warp_size_zu);
    se_top_result_offset += roundUp<size_t>((ni_top - 1LLU) * ni_top / 2LLU, warp_size_zu);
    se_lbl_result_offset += roundUp<size_t>((ni_lbl - 1LLU) * ni_lbl / 2LLU, warp_size_zu);
    nr_src_result_offset += roundUp(ni_src * ni_src, warp_size_zu);
    nr_top_result_offset += roundUp(ni_top * ni_top, warp_size_zu);
    nr_lbl_result_offset += roundUp(ni_lbl * ni_lbl, warp_size_zu);
  }

  // Set the pair result allocation sizes.
  sepair_result_alloc_src = se_src_result_offset;
  sepair_result_alloc_top = se_top_result_offset;
  sepair_result_alloc_lbl = se_lbl_result_offset;
  nrpair_result_alloc_src = nr_src_result_offset;
  nrpair_result_alloc_top = nr_top_result_offset;
  nrpair_result_alloc_lbl = nr_lbl_result_offset;
}
  
//-------------------------------------------------------------------------------------------------
void ComparisonGuide::generateWorkUnits(const int* system_list, const int* topology_index_list,
                                        const int* bounds_list, const int partitions,
                                        Hybrid<int4> *atr_insr_members,
                                        Hybrid<int> *atr_insr_groups,
                                        Hybrid<int4> *ata_insr_members,
                                        Hybrid<int> *ata_insr_groups,
                                        int *atr_instruction_count, int *ata_instruction_count,
                                        const GpuDetails &gpu) {
  const int system_count = bounds_list[partitions];

  // Compute the number of warps on the GPU and the number of bits per unsigned int.
  const llint gpu_nwarp = (gpu == null_gpu) ?
                          large_block_size / warp_size_int :
                          gpu.getSMPCount() * (gpu.getMaxThreadsPerBlock() / warp_size_int);
  const int uint_bits = sizeof(uint) * 8;
  const int half_int_bits = uint_bits / 2;
  const llint np_atr = system_count;
  llint np_ata = 0LL;
  for (int i = 0; i < partitions; i++) {
    const llint ni = bounds_list[i + 1] - bounds_list[i];
    np_ata += ni * ni;
  }
  llint atr_per_warp = ceil(static_cast<double>(np_atr) / static_cast<double>(gpu_nwarp));
  atr_per_warp = std::min(static_cast<llint>(uint_bits), atr_per_warp);
  llint ata_per_warp = ceil(static_cast<double>(np_ata) / static_cast<double>(gpu_nwarp));
  const int ata_stride = std::min(static_cast<llint>(half_int_bits),
                                  static_cast<llint>(ceil(sqrt(ata_per_warp))));;
  ata_per_warp = ata_stride * ata_stride;

  // Lay out a temporary array for all-to-one instructions.
  std::vector<int4> tmp_atr_insr_members;
  std::vector<int> tmp_atr_insr_groups;
  tmp_atr_insr_members.reserve((np_atr / atr_per_warp) + partitions);
  tmp_atr_insr_groups.reserve(tmp_atr_insr_members.size());
  for (int i = 0; i < partitions; i++) {
    const int group_lower_bound = bounds_list[i];
    const int group_upper_bound = bounds_list[i + 1];
    const int top_idx = topology_index_list[group_lower_bound];
    for (int j = group_lower_bound; j < group_upper_bound; j += atr_per_warp) {
      const int nrep = std::min(static_cast<int>(atr_per_warp), group_upper_bound - j);
      tmp_atr_insr_members.push_back({ j - group_lower_bound, i, top_idx, nrep });
      tmp_atr_insr_groups.push_back(i);
    }
  }
  *atr_instruction_count = tmp_atr_insr_members.size();

  // Lay out a temporary array for all-to-all instructions.
  std::vector<int4> tmp_ata_insr_members;
  std::vector<int> tmp_ata_insr_groups;
  size_t ata_instruction_guess = 0;
  for (int i = 0; i < partitions; i++) {
    const size_t nstride = (bounds_list[i + 1] - bounds_list[i] + ata_stride - 1) / ata_stride;
    ata_instruction_guess += nstride * (nstride + 1) / 2;
  }
  tmp_ata_insr_members.reserve(ata_instruction_guess);
  tmp_ata_insr_groups.reserve(ata_instruction_guess);
  for (int i = 0; i < partitions; i++) {
    const int group_lower_bound = bounds_list[i];
    const int group_upper_bound = bounds_list[i + 1];
    const int top_idx = topology_index_list[group_lower_bound];
    for (int j = group_lower_bound; j < group_upper_bound; j += ata_stride) {
      const int njrep = (j + ata_stride < group_upper_bound) ? ata_stride : group_upper_bound - j;
      for (int k = group_lower_bound; k <= j; k += ata_stride) {
	const int nkrep = (k + ata_stride < group_upper_bound) ? ata_stride :
                                                                 group_upper_bound - k;
	tmp_ata_insr_members.push_back({ j - group_lower_bound, k - group_lower_bound, top_idx,
                                         (nkrep << half_int_bits) | njrep });
	tmp_ata_insr_groups.push_back(i);
      }
    }
  }
  *ata_instruction_count = tmp_ata_insr_members.size();

  // Load the temporary arrays into the object.
  atr_insr_members->resize(*atr_instruction_count);
  ata_insr_members->resize(*ata_instruction_count);
  atr_insr_members->putHost(tmp_atr_insr_members);
  ata_insr_members->putHost(tmp_ata_insr_members);
  atr_insr_groups->resize(*atr_instruction_count);
  ata_insr_groups->resize(*ata_instruction_count);
  atr_insr_groups->putHost(tmp_atr_insr_groups);
  ata_insr_groups->putHost(tmp_ata_insr_groups);
}

//-------------------------------------------------------------------------------------------------
void ComparisonGuide::setWorkUnits(const GpuDetails &gpu) {
  generateWorkUnits(cache_source_groups.data(), cache_source_group_topologies.data(),
                    cache_source_group_bounds.data(), cache_source_group_bounds.size() - 1,
                    &atr_instruction_members_src, &atr_instruction_groups_src,
                    &ata_instruction_members_src, &ata_instruction_groups_src,
                    &atr_instruction_count_src, &ata_instruction_count_src, gpu);
  generateWorkUnits(topology_groups.data(), topology_group_topologies.data(),
                    topology_group_bounds.data(), topology_group_bounds.size() - 1,
                    &atr_instruction_members_top, &atr_instruction_groups_top,
                    &ata_instruction_members_top, &ata_instruction_groups_top,
                    &atr_instruction_count_top, &ata_instruction_count_top, gpu);
  generateWorkUnits(cache_label_groups.data(), cache_label_group_topologies.data(),
                    cache_label_group_bounds.data(), cache_label_group_bounds.size() - 1,
                    &atr_instruction_members_lbl, &atr_instruction_groups_lbl,
                    &ata_instruction_members_lbl, &ata_instruction_groups_lbl,
                    &atr_instruction_count_lbl, &ata_instruction_count_lbl, gpu);
}

//-------------------------------------------------------------------------------------------------
void ComparisonGuide::validatePartitionIndex(const int index, const int member,
                                             const SystemGrouping organization,
                                             const char* caller) const {
  if (index < 0 || index > getPartitionCount(organization)) {
    rtErr(getEnumerationName(organization) + " index " + std::to_string(index) + " is invalid "
          "for a systems cache / synthesis with " +
          std::to_string(getPartitionCount(organization)) + " distinct partitions.",
          "ComparisonGuide", caller);
  }
  int ngroup;
  switch (organization) {
  case SystemGrouping::SOURCE:
    ngroup = cache_source_group_bounds.readHost(index + 1) -
             cache_source_group_bounds.readHost(index);
    break;
  case SystemGrouping::TOPOLOGY:
    ngroup = topology_group_bounds.readHost(index + 1) - topology_group_bounds.readHost(index);
    break;
  case SystemGrouping::LABEL:
    ngroup = cache_label_group_bounds.readHost(index + 1) -
             cache_label_group_bounds.readHost(index);
    break;
  }
  if (member >= ngroup && member > 0) {
    rtErr("Group member " + std::to_string(member) + " is invalid for partition " +
          std::to_string(index) + " when dividing systems by " + getEnumerationName(organization) +
          ".", "ComparisonGuide", caller);
  }
}

//-------------------------------------------------------------------------------------------------
void ComparisonGuide::validatePartitionIndex(const int index, const SystemGrouping organization,
                                             const char* caller) const {
  validatePartitionIndex(index, 0, organization, caller);
}

} // namespace analysis
} // namespace stormm
