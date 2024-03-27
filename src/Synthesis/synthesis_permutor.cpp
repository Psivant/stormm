#include <limits.h>
#include "copyright.h"
#include "Chemistry/chemistry_enumerators.h"
#include "Constants/symbol_values.h"
#include "FileManagement/file_listing.h"
#include "Math/clustering.h"
#include "Math/cluster_manager.h"
#include "Math/summation.h"
#include "Math/vector_ops.h"
#include "Potential/static_exclusionmask.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/topology_util.h"
#include "Structure/structure_enumerators.h"
#include "Trajectory/coordinate_copy.h"
#include "synthesis_permutor.h"

namespace stormm {
namespace synthesis {

using card::HybridKind;
using chemistry::MapRotatableGroups;
using chemistry::permutationsAreLinked;
using diskutil::getBaseName;
using energy::StaticExclusionMask;
using namelist::default_conf_cis_trans_samples;
using namelist::default_conf_cis_trans_set_zero;
using namelist::default_conf_cis_trans_set_one;
using namelist::default_conf_cis_trans_snap;
using namelist::default_conf_rotation_samples;
using namelist::default_conf_rotation_set_zero;
using namelist::default_conf_rotation_set_one;
using namelist::default_conf_rotation_set_two;
using namelist::default_conf_rotation_snap;
using numerics::globalpos_scale_nonoverflow_bits;
using stmath::addScalarToVector;
using stmath::ClusterManager;
using stmath::kMeans;
using stmath::minValue;
using stmath::prefixSumInPlace;
using stmath::PrefixSumType;
using stmath::accumulateBitmask;
using stmath::readBitFromMask;
using stmath::sum;
using structure::BoundaryCondition;
using structure::imageValue;
using structure::ImagingMethod;
using symbols::pi;
using topology::matchTopology;
using topology::NonbondedKit;
using trajectory::coordCopy;

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor() :
    system_count{0}, permutor_map_count{0},
    rotatable_bond_samples{default_conf_rotation_samples},
    cis_trans_bond_samples{default_conf_cis_trans_samples},
    features_objects_accessible{false},
    rotation_setting_snap{stod(std::string(default_conf_rotation_snap)) * pi / 180.0}, 
    cis_trans_setting_snap{stod(std::string(default_conf_cis_trans_snap)) * pi / 180.0},
    total_rotatable_groups{0}, total_cis_trans_groups{0}, total_invertible_groups{0},
    permutor_map_indices{HybridKind::POINTER, "synp_perm_map_idx"},
    permutor_elements{HybridKind::POINTER, "synp_perm_elements"},
    permutor_element_bounds{HybridKind::POINTER, "synp_perm_ele_bnds"},
    system_settings{HybridKind::POINTER, "synp_sys_sett"},
    system_settings_limits{HybridKind::POINTER, "synp_sys_sett_lim"},
    synthesis_data{HybridKind::ARRAY, "synp_synth_data"},
    rotatable_group_atoms{HybridKind::POINTER, "synp_rot_grp_atoms"},
    rotatable_group_bounds{HybridKind::POINTER, "synp_rot_grp_bnds"},
    permutor_rotatable_group_bounds{HybridKind::POINTER, "synp_prot_grp_bnds"},
    cis_trans_group_atoms{HybridKind::POINTER, "synp_ctx_grp_atoms"},
    cis_trans_group_bounds{HybridKind::POINTER, "synp_ctx_grp_bnds"},
    permutor_cis_trans_group_bounds{HybridKind::POINTER, "synp_pctx_grp_bnds"},
    invertible_group_atoms{HybridKind::POINTER, "synp_inv_grp_atoms"},
    invertible_group_bounds{HybridKind::POINTER, "synp_inv_grp_bnds"},
    permutor_invertible_group_bounds{HybridKind::POINTER, "synp_pinv_grp_bnds"},
    invertible_atom_centers{HybridKind::POINTER, "synp_chiral_centers"},
    invertible_group_protocols{HybridKind::POINTER, "synp_chiral_protocols"},
    rotatable_bond_markers{HybridKind::POINTER, "synp_rot_bond_mark"},
    cis_trans_bond_markers{HybridKind::POINTER, "synp_ctx_bond_mark"},
    chiral_markers{HybridKind::POINTER, "synp_chiral_mark"},
    rotatable_bond_settings{HybridKind::ARRAY, "synp_rot_bond_sett"},
    cis_trans_bond_settings{HybridKind::ARRAY, "synp_ctx_bond_sett"},
    sp_rotatable_bond_settings{HybridKind::ARRAY, "synp_rot_bond_sett_sp"},
    sp_cis_trans_bond_settings{HybridKind::ARRAY, "synp_ctx_bond_sett_sp"},
    rotatable_bond_settings_bounds{HybridKind::POINTER, "synp_rot_bond_sbnd"},
    cis_trans_bond_settings_bounds{HybridKind::POINTER, "synp_ctx_bond_sbnd"},
    chiral_settings_bounds{HybridKind::POINTER, "synp_chiral_sbnd"},
    group_data{HybridKind::ARRAY, "synp_group_data"},
    marker_data{HybridKind::ARRAY, "synp_marker_data"},
    chiral_settings{HybridKind::ARRAY, "synp_chiral_sett"},
    light_sampling_foci{HybridKind::POINTER, "synp_light_foci"},
    heavy_sampling_foci{HybridKind::POINTER, "synp_heavy_foci"},
    light_sampling_foci_bounds{HybridKind::POINTER, "synp_light_foci_bounds"},
    heavy_sampling_foci_bounds{HybridKind::POINTER, "synp_heavy_foci_bounds"},
    variables_data{HybridKind::ARRAY, "synp_var_data"},
    general_rotation_settings{}, general_cis_trans_settings{}, rotatable_groups{},
    cis_trans_groups{}, invertible_groups{}, state_trackers{},
    minimal_sampling_replicas{}, light_sampling_replicas{}, heavy_sampling_replicas{},
    exhaustive_sampling_replicas{}, topologies{}, features{}, poly_ps_ptr{nullptr}
{
  // Load the default settings for rotatable bonds
  general_rotation_settings.resize(rotatable_bond_samples);
  general_rotation_settings[0] = stod(std::string(default_conf_rotation_set_zero)) * pi / 180.0;
  general_rotation_settings[1] = stod(std::string(default_conf_rotation_set_one)) * pi / 180.0;
  general_rotation_settings[2] = stod(std::string(default_conf_rotation_set_two)) * pi / 180.0;
  general_cis_trans_settings.resize(cis_trans_bond_samples);
  general_cis_trans_settings[0] = stod(std::string(default_conf_cis_trans_set_zero)) * pi / 180.0;
  general_cis_trans_settings[1] = stod(std::string(default_conf_cis_trans_set_one)) * pi / 180.0;
}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const AtomGraphSynthesis *poly_ag, StopWatch *timer) :
    SynthesisPermutor()
{
  system_count = poly_ag->getSystemCount();
  permutor_map_count = poly_ag->getUniqueTopologyCount();
  topologies = poly_ag->getUniqueTopologies();
  const std::vector<ChemicalFeatures> all_chemfe = temporaryFeatures(timer);
  const std::vector<ChemicalFeatures*> all_chemfe_ptr = temporaryFeaturesPointers(all_chemfe);  
  fillPermutorMaps(all_chemfe_ptr, timer);
  setVariableRanges(general_rotation_settings, general_cis_trans_settings, true, true);
}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const AtomGraphSynthesis &poly_ag, StopWatch *timer) :
    SynthesisPermutor(poly_ag.getSelfPointer(), timer)
{}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const AtomGraphSynthesis *poly_ag,
                                     const ConformerControls &confcon, StopWatch *timer) :
    SynthesisPermutor()
{
  system_count = poly_ag->getSystemCount();
  permutor_map_count = poly_ag->getUniqueTopologyCount();
  topologies = poly_ag->getUniqueTopologies();
  const std::vector<ChemicalFeatures> all_chemfe = temporaryFeatures(timer);
  const std::vector<ChemicalFeatures*> all_chemfe_ptr = temporaryFeaturesPointers(all_chemfe);  
  fillPermutorMaps(all_chemfe_ptr, timer);
  impartControlData(confcon);
}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const AtomGraphSynthesis &poly_ag,
                                     const ConformerControls &confcon, StopWatch *timer) :
    SynthesisPermutor(poly_ag.getSelfPointer(), confcon, timer)
{}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const PhaseSpaceSynthesis *poly_ps, StopWatch *timer) :
    SynthesisPermutor()
{
  system_count = poly_ps->getSystemCount();
  permutor_map_count = poly_ps->getUniqueTopologyCount();
  topologies = poly_ps->getUniqueTopologies();
  const std::vector<ChemicalFeatures> all_chemfe = temporaryFeatures(timer);
  const std::vector<ChemicalFeatures*> all_chemfe_ptr = temporaryFeaturesPointers(all_chemfe);
  fillPermutorMaps(all_chemfe_ptr, timer);
  setVariableRanges(general_rotation_settings, general_cis_trans_settings, true, true);
  applySynthesis(poly_ps, VariableTorsionAdjustment::DO_NOT_CHANGE);
}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const PhaseSpaceSynthesis &poly_ps, StopWatch *timer) :
    SynthesisPermutor(poly_ps.getSelfPointer(), timer)
{}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const PhaseSpaceSynthesis *poly_ps,
                                     const ConformerControls &confcon, StopWatch *timer) :
    SynthesisPermutor()
{
  system_count = poly_ps->getSystemCount();
  permutor_map_count = poly_ps->getUniqueTopologyCount();
  topologies = poly_ps->getUniqueTopologies();
  const std::vector<ChemicalFeatures> all_chemfe = temporaryFeatures(timer);
  const std::vector<ChemicalFeatures*> all_chemfe_ptr = temporaryFeaturesPointers(all_chemfe);
  fillPermutorMaps(all_chemfe_ptr, timer);
  impartControlData(confcon);
  applySynthesis(poly_ps, confcon.getTorsionAdjustmentProtocol());
}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const PhaseSpaceSynthesis &poly_ps,
                                     const ConformerControls &confcon, StopWatch *timer) :
    SynthesisPermutor(poly_ps.getSelfPointer(), confcon, timer)
{}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const std::vector<ChemicalFeatures*> &chemfe_in,
                                     const bool retain_pointers, StopWatch *timer) :
    SynthesisPermutor()
{
  permutor_map_count = chemfe_in.size();
  features_objects_accessible = retain_pointers;
  topologies.resize(permutor_map_count);
  for (int i = 0; i < permutor_map_count; i++) {
    topologies[i] = const_cast<AtomGraph*>(chemfe_in[i]->getTopologyPointer());
  }
  fillPermutorMaps(chemfe_in, timer);
  setVariableRanges(general_rotation_settings, general_cis_trans_settings, true, true);
}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const std::vector<ChemicalFeatures> &chemfe_in,
                                     const bool retain_pointers, StopWatch *timer) :
  SynthesisPermutor(temporaryFeaturesPointers(chemfe_in), retain_pointers, timer)
{}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const std::vector<ChemicalFeatures*> &chemfe_in,
                                     const PhaseSpaceSynthesis *poly_ps,
                                     const ConformerControls &confcon, const bool retain_pointers,
                                     StopWatch *timer) :
    SynthesisPermutor()
{
  permutor_map_count = chemfe_in.size();
  features_objects_accessible = retain_pointers;
  topologies.resize(permutor_map_count);
  for (int i = 0; i < permutor_map_count; i++) {
    topologies[i] = const_cast<AtomGraph*>(chemfe_in[i]->getTopologyPointer());
  }
  fillPermutorMaps(chemfe_in, timer);
  impartControlData(confcon);
  applySynthesis(poly_ps, confcon.getTorsionAdjustmentProtocol());
}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const std::vector<ChemicalFeatures*> &chemfe_in,
                                     const PhaseSpaceSynthesis &poly_ps,
                                     const ConformerControls &confcon, const bool retain_pointers,
                                     StopWatch *timer) :
    SynthesisPermutor(chemfe_in, poly_ps.getSelfPointer(), confcon, retain_pointers, timer)
{}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const std::vector<ChemicalFeatures> &chemfe_in,
                                     const PhaseSpaceSynthesis &poly_ps,
                                     const ConformerControls &confcon, const bool retain_pointers,
                                     StopWatch *timer) :
    SynthesisPermutor(temporaryFeaturesPointers(chemfe_in), poly_ps.getSelfPointer(), confcon,
                      retain_pointers, timer)
{}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(const SynthesisPermutor &original) :
  system_count{original.system_count}, permutor_map_count{original.permutor_map_count},
  rotatable_bond_samples{original.rotatable_bond_samples},
  cis_trans_bond_samples{original.cis_trans_bond_samples},
  features_objects_accessible{original.features_objects_accessible},
  rotation_setting_snap{original.rotation_setting_snap},
  cis_trans_setting_snap{original.cis_trans_setting_snap},
  total_rotatable_groups{original.total_rotatable_groups},
  total_cis_trans_groups{original.total_cis_trans_groups},
  total_invertible_groups{original.total_invertible_groups},
  permutor_map_indices{original.permutor_map_indices},
  permutor_elements{original.permutor_elements},
  permutor_element_bounds{original.permutor_element_bounds},
  system_settings{original.system_settings},
  system_settings_limits{original.system_settings_limits},
  synthesis_data{original.synthesis_data},
  rotatable_group_atoms{original.rotatable_group_atoms},
  rotatable_group_bounds{original.rotatable_group_bounds},
  permutor_rotatable_group_bounds{original.permutor_rotatable_group_bounds},
  cis_trans_group_atoms{original.cis_trans_group_atoms},
  cis_trans_group_bounds{original.cis_trans_group_bounds},
  permutor_cis_trans_group_bounds{original.permutor_cis_trans_group_bounds},
  invertible_group_atoms{original.invertible_group_atoms},
  invertible_group_bounds{original.invertible_group_bounds},
  permutor_invertible_group_bounds{original.permutor_invertible_group_bounds},
  invertible_atom_centers{original.invertible_atom_centers},
  invertible_group_protocols{original.invertible_group_protocols},
  rotatable_bond_markers{original.rotatable_bond_markers},
  cis_trans_bond_markers{original.cis_trans_bond_markers},
  chiral_markers{original.chiral_markers},
  rotatable_bond_settings{original.rotatable_bond_settings},
  cis_trans_bond_settings{original.cis_trans_bond_settings},
  sp_rotatable_bond_settings{original.sp_rotatable_bond_settings},
  sp_cis_trans_bond_settings{original.sp_cis_trans_bond_settings},
  rotatable_bond_settings_bounds{original.rotatable_bond_settings_bounds},
  cis_trans_bond_settings_bounds{original.cis_trans_bond_settings_bounds},
  chiral_settings_bounds{original.chiral_settings_bounds},
  group_data{original.group_data},
  marker_data{original.marker_data},
  chiral_settings{original.chiral_settings},
  light_sampling_foci{original.light_sampling_foci},
  heavy_sampling_foci{original.heavy_sampling_foci},
  light_sampling_foci_bounds{original.light_sampling_foci_bounds},
  heavy_sampling_foci_bounds{original.heavy_sampling_foci_bounds},
  variables_data{original.variables_data},
  general_rotation_settings{original.general_rotation_settings},
  general_cis_trans_settings{original.general_cis_trans_settings},
  rotatable_groups{original.rotatable_groups},
  cis_trans_groups{original.cis_trans_groups},
  invertible_groups{original.invertible_groups},
  state_trackers{original.state_trackers},
  minimal_sampling_replicas{original.minimal_sampling_replicas},
  light_sampling_replicas{original.light_sampling_replicas},
  heavy_sampling_replicas{original.heavy_sampling_replicas},
  exhaustive_sampling_replicas{original.exhaustive_sampling_replicas},
  topologies{original.topologies},
  features{original.features}
{
  // Rebase the POINTER-kind Hybrid objects 
  rebasePointers();
}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor::SynthesisPermutor(SynthesisPermutor &&original) :
  system_count{original.system_count}, permutor_map_count{original.permutor_map_count},
  rotatable_bond_samples{original.rotatable_bond_samples},
  cis_trans_bond_samples{original.cis_trans_bond_samples},
  features_objects_accessible{original.features_objects_accessible},
  rotation_setting_snap{original.rotation_setting_snap},
  cis_trans_setting_snap{original.cis_trans_setting_snap},
  total_rotatable_groups{original.total_rotatable_groups},
  total_cis_trans_groups{original.total_cis_trans_groups},
  total_invertible_groups{original.total_invertible_groups},
  permutor_map_indices{std::move(original.permutor_map_indices)},
  permutor_elements{std::move(original.permutor_elements)},
  permutor_element_bounds{std::move(original.permutor_element_bounds)},
  system_settings{std::move(original.system_settings)},
  system_settings_limits{std::move(original.system_settings_limits)},
  synthesis_data{std::move(original.synthesis_data)},
  rotatable_group_atoms{std::move(original.rotatable_group_atoms)},
  rotatable_group_bounds{std::move(original.rotatable_group_bounds)},
  permutor_rotatable_group_bounds{std::move(original.permutor_rotatable_group_bounds)},
  cis_trans_group_atoms{std::move(original.cis_trans_group_atoms)},
  cis_trans_group_bounds{std::move(original.cis_trans_group_bounds)},
  permutor_cis_trans_group_bounds{std::move(original.permutor_cis_trans_group_bounds)},
  invertible_group_atoms{std::move(original.invertible_group_atoms)},
  invertible_group_bounds{std::move(original.invertible_group_bounds)},
  permutor_invertible_group_bounds{std::move(original.permutor_invertible_group_bounds)},
  invertible_atom_centers{std::move(original.invertible_atom_centers)},
  invertible_group_protocols{std::move(original.invertible_group_protocols)},
  rotatable_bond_markers{std::move(original.rotatable_bond_markers)},
  cis_trans_bond_markers{std::move(original.cis_trans_bond_markers)},
  chiral_markers{std::move(original.chiral_markers)},
  rotatable_bond_settings{std::move(original.rotatable_bond_settings)},
  cis_trans_bond_settings{std::move(original.cis_trans_bond_settings)},
  sp_rotatable_bond_settings{std::move(original.sp_rotatable_bond_settings)},
  sp_cis_trans_bond_settings{std::move(original.sp_cis_trans_bond_settings)},
  rotatable_bond_settings_bounds{std::move(original.rotatable_bond_settings_bounds)},
  cis_trans_bond_settings_bounds{std::move(original.cis_trans_bond_settings_bounds)},
  chiral_settings_bounds{std::move(original.chiral_settings_bounds)},
  group_data{std::move(original.group_data)},
  marker_data{std::move(original.marker_data)},
  chiral_settings{std::move(original.chiral_settings)},
  light_sampling_foci{std::move(original.light_sampling_foci)},
  heavy_sampling_foci{std::move(original.heavy_sampling_foci)},
  light_sampling_foci_bounds{std::move(original.light_sampling_foci_bounds)},
  heavy_sampling_foci_bounds{std::move(original.heavy_sampling_foci_bounds)},
  variables_data{std::move(original.variables_data)},
  general_rotation_settings{std::move(original.general_rotation_settings)},
  general_cis_trans_settings{std::move(original.general_cis_trans_settings)},
  rotatable_groups{std::move(original.rotatable_groups)},
  cis_trans_groups{std::move(original.cis_trans_groups)},
  invertible_groups{std::move(original.invertible_groups)},
  state_trackers{std::move(original.state_trackers)},
  minimal_sampling_replicas{std::move(original.minimal_sampling_replicas)},
  light_sampling_replicas{std::move(original.light_sampling_replicas)},
  heavy_sampling_replicas{std::move(original.heavy_sampling_replicas)},
  exhaustive_sampling_replicas{std::move(original.exhaustive_sampling_replicas)},
  topologies{std::move(original.topologies)},
  features{std::move(original.features)}
{}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor& SynthesisPermutor::operator=(const SynthesisPermutor &other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }

  // Copy all member variables
  system_count = other.system_count;
  permutor_map_count = other.permutor_map_count;
  rotatable_bond_samples = other.rotatable_bond_samples;
  cis_trans_bond_samples = other.cis_trans_bond_samples;
  features_objects_accessible = other.features_objects_accessible;
  rotation_setting_snap = other.rotation_setting_snap;
  cis_trans_setting_snap = other.cis_trans_setting_snap;
  total_rotatable_groups = other.total_rotatable_groups;
  total_cis_trans_groups = other.total_cis_trans_groups;
  total_invertible_groups = other.total_invertible_groups;
  permutor_map_indices = other.permutor_map_indices;
  permutor_elements = other.permutor_elements;
  permutor_element_bounds = other.permutor_element_bounds;
  system_settings = other.system_settings;
  system_settings_limits = other.system_settings_limits;
  synthesis_data = other.synthesis_data;
  rotatable_group_atoms = other.rotatable_group_atoms;
  rotatable_group_bounds = other.rotatable_group_bounds;
  permutor_rotatable_group_bounds = other.permutor_rotatable_group_bounds;
  cis_trans_group_atoms = other.cis_trans_group_atoms;
  cis_trans_group_bounds = other.cis_trans_group_bounds;
  permutor_cis_trans_group_bounds = other.permutor_cis_trans_group_bounds;
  invertible_group_atoms = other.invertible_group_atoms;
  invertible_group_bounds = other.invertible_group_bounds;
  permutor_invertible_group_bounds = other.permutor_invertible_group_bounds;
  invertible_atom_centers = other.invertible_atom_centers;
  invertible_group_protocols = other.invertible_group_protocols;
  rotatable_bond_markers = other.rotatable_bond_markers;
  cis_trans_bond_markers = other.cis_trans_bond_markers;
  chiral_markers = other.chiral_markers;
  rotatable_bond_settings = other.rotatable_bond_settings;
  cis_trans_bond_settings = other.cis_trans_bond_settings;
  sp_rotatable_bond_settings = other.sp_rotatable_bond_settings;
  sp_cis_trans_bond_settings = other.sp_cis_trans_bond_settings;
  rotatable_bond_settings_bounds = other.rotatable_bond_settings_bounds;
  cis_trans_bond_settings_bounds = other.cis_trans_bond_settings_bounds;
  chiral_settings_bounds = other.chiral_settings_bounds;
  group_data = other.group_data;
  marker_data = other.marker_data;
  chiral_settings = other.chiral_settings;
  light_sampling_foci = other.light_sampling_foci;
  heavy_sampling_foci = other.heavy_sampling_foci;
  light_sampling_foci_bounds = other.light_sampling_foci_bounds;
  heavy_sampling_foci_bounds = other.heavy_sampling_foci_bounds;
  variables_data = other.variables_data;
  general_rotation_settings = other.general_rotation_settings;
  general_cis_trans_settings = other.general_cis_trans_settings;
  rotatable_groups = other.rotatable_groups;
  cis_trans_groups = other.cis_trans_groups;
  invertible_groups = other.invertible_groups;
  state_trackers = other.state_trackers;
  minimal_sampling_replicas = other.minimal_sampling_replicas;
  light_sampling_replicas = other.light_sampling_replicas;
  heavy_sampling_replicas = other.heavy_sampling_replicas;
  exhaustive_sampling_replicas = other.exhaustive_sampling_replicas;
  topologies = other.topologies;
  features = other.features;

  // Rebase pointers and return
  rebasePointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
SynthesisPermutor& SynthesisPermutor::operator=(SynthesisPermutor &&other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }

  // Copy all member variables
  system_count = other.system_count;
  permutor_map_count = other.permutor_map_count;
  rotatable_bond_samples = other.rotatable_bond_samples;
  cis_trans_bond_samples = other.cis_trans_bond_samples;
  features_objects_accessible = other.features_objects_accessible;
  rotation_setting_snap = other.rotation_setting_snap;
  cis_trans_setting_snap = other.cis_trans_setting_snap;
  total_rotatable_groups = other.total_rotatable_groups;
  total_cis_trans_groups = other.total_cis_trans_groups;
  total_invertible_groups = other.total_invertible_groups;
  permutor_map_indices = std::move(other.permutor_map_indices);
  permutor_elements = std::move(other.permutor_elements);
  permutor_element_bounds = std::move(other.permutor_element_bounds);
  system_settings = std::move(other.system_settings);
  system_settings_limits = std::move(other.system_settings_limits);
  synthesis_data = std::move(other.synthesis_data);
  rotatable_group_atoms = std::move(other.rotatable_group_atoms);
  rotatable_group_bounds = std::move(other.rotatable_group_bounds);
  permutor_rotatable_group_bounds = std::move(other.permutor_rotatable_group_bounds);
  cis_trans_group_atoms = std::move(other.cis_trans_group_atoms);
  cis_trans_group_bounds = std::move(other.cis_trans_group_bounds);
  permutor_cis_trans_group_bounds = std::move(other.permutor_cis_trans_group_bounds);
  invertible_group_atoms = std::move(other.invertible_group_atoms);
  invertible_group_bounds = std::move(other.invertible_group_bounds);
  permutor_invertible_group_bounds = std::move(other.permutor_invertible_group_bounds);
  invertible_atom_centers = std::move(other.invertible_atom_centers);
  invertible_group_protocols = std::move(other.invertible_group_protocols);
  rotatable_bond_markers = std::move(other.rotatable_bond_markers);
  cis_trans_bond_markers = std::move(other.cis_trans_bond_markers);
  chiral_markers = std::move(other.chiral_markers);
  rotatable_bond_settings = std::move(other.rotatable_bond_settings);
  cis_trans_bond_settings = std::move(other.cis_trans_bond_settings);
  sp_rotatable_bond_settings = std::move(other.sp_rotatable_bond_settings);
  sp_cis_trans_bond_settings = std::move(other.sp_cis_trans_bond_settings);
  rotatable_bond_settings_bounds = std::move(other.rotatable_bond_settings_bounds);
  cis_trans_bond_settings_bounds = std::move(other.cis_trans_bond_settings_bounds);
  chiral_settings_bounds = std::move(other.chiral_settings_bounds);
  group_data = std::move(other.group_data);
  marker_data = std::move(other.marker_data);
  chiral_settings = std::move(other.chiral_settings);
  light_sampling_foci = std::move(other.light_sampling_foci);
  heavy_sampling_foci = std::move(other.heavy_sampling_foci);
  light_sampling_foci_bounds = std::move(other.light_sampling_foci_bounds);
  heavy_sampling_foci_bounds = std::move(other.heavy_sampling_foci_bounds);
  variables_data = std::move(other.variables_data);
  general_rotation_settings = std::move(other.general_rotation_settings);
  general_cis_trans_settings = std::move(other.general_cis_trans_settings);
  rotatable_groups = std::move(other.rotatable_groups);
  cis_trans_groups = std::move(other.cis_trans_groups);
  invertible_groups = std::move(other.invertible_groups);
  state_trackers = std::move(other.state_trackers);
  minimal_sampling_replicas = std::move(other.minimal_sampling_replicas);
  light_sampling_replicas = std::move(other.light_sampling_replicas);
  heavy_sampling_replicas = std::move(other.heavy_sampling_replicas);
  exhaustive_sampling_replicas = std::move(other.exhaustive_sampling_replicas);
  topologies = std::move(other.topologies);
  features = std::move(other.features);

  // As with other move assignment operations, rebasing is not needed.  Simply return the result.
  return *this;
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getPermutorSetCount() const {
  return permutor_map_count;
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getPermutorMapIndex(const AtomGraph *query_ag) const {
  const int map_idx = matchTopology(query_ag, topologies);
  if (map_idx == permutor_map_count) {
    rtErr("No topology originating in file " + query_ag->getFileName() + " was found.",
          "SynthesisPermutor", "getPermutorMapIndex");
  }
  return map_idx;
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getPermutorMapIndex(const AtomGraph &query_ag) const {
  return getPermutorMapIndex(query_ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getPermutorMapIndex(const PhaseSpaceSynthesis *poly_ps,
                                           const int system_index) const {
  return getPermutorMapIndex(poly_ps->getSystemTopologyPointer(system_index));
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getPermutorMapIndex(const PhaseSpaceSynthesis &poly_ps,
                                           const int system_index) const {
  return getPermutorMapIndex(poly_ps.getSelfPointer(), system_index);
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getRotatableBondSampleCount() const {
  return rotatable_bond_samples;
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getCisTransBondSampleCount() const {
  return cis_trans_bond_samples;
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getRotatableBondCount(const int permutor_map_index) const {
  validateMapIndex(permutor_map_index);
  return permutor_rotatable_group_bounds.readHost(permutor_map_index + 1) -
         permutor_rotatable_group_bounds.readHost(permutor_map_index);
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getRotatableBondCount() const {
  return total_rotatable_groups;
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getCisTransBondCount(const int permutor_map_index) const {
  validateMapIndex(permutor_map_index);
  return permutor_cis_trans_group_bounds.readHost(permutor_map_index + 1) -
         permutor_cis_trans_group_bounds.readHost(permutor_map_index);
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getCisTransBondCount() const {
  return total_cis_trans_groups;
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getChiralCenterCount(const int permutor_map_index) const {
  validateMapIndex(permutor_map_index);
  return permutor_invertible_group_bounds.readHost(permutor_map_index + 1) -
         permutor_invertible_group_bounds.readHost(permutor_map_index);
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getChiralCenterCount() const {
  return total_invertible_groups;
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& SynthesisPermutor::getElementSampleCounts(int system_index) const {
  validateSystemIndex(system_index);
  const int map_idx = permutor_map_indices.readHost(system_index);
  return state_trackers[map_idx].getStateLimits();
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>&
SynthesisPermutor::getElementSampleCounts(const AtomGraph *query_ag) const {
  for (int i = 0; i < permutor_map_count; i++) {
    if (topologies[i] == query_ag) {
      return state_trackers[i].getStateLimits();
    }
  }
  rtErr("No topology match was found for an object originating in file " +
        getBaseName(query_ag->getFileName()) + ".", "SynthesisPermutor", "getElementSampleCounts");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>&
SynthesisPermutor::getElementSampleCounts(const AtomGraph &query_ag) const {
  return getElementSampleCounts(query_ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getVariableElementCount(const int system_index) const {
  validateSystemIndex(system_index);
  const int map_index = permutor_map_indices.readHost(system_index);
  return (permutor_rotatable_group_bounds.readHost(map_index + 1) -
          permutor_rotatable_group_bounds.readHost(map_index)) +
         (permutor_cis_trans_group_bounds.readHost(map_index + 1) -
          permutor_cis_trans_group_bounds.readHost(map_index)) +
         (permutor_invertible_group_bounds.readHost(map_index + 1) -
          permutor_invertible_group_bounds.readHost(map_index));
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getElementSampleCount(const ConformationEdit element_kind,
                                             const int index) const {

  // Some blank entries in tuples will reference the -1 element.  Return one in these cases.
  if (index < 0) {
    return 1;
  }
  switch (element_kind) {
  case ConformationEdit::BOND_ROTATION:
    return rotatable_bond_settings_bounds.readHost(index + 1) -
           rotatable_bond_settings_bounds.readHost(index);
  case ConformationEdit::CIS_TRANS_FLIP:
    return cis_trans_bond_settings_bounds.readHost(index + 1) -
           cis_trans_bond_settings_bounds.readHost(index);
  case ConformationEdit::CHIRAL_INVERSION:
    return chiral_settings_bounds.readHost(index + 1) - chiral_settings_bounds.readHost(index);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getElementSampleCount(const CoupledEdit ce) const {
  return getElementSampleCount(ce.edit, ce.index);
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getSystemRotatableBondCount(const int system_index) const {
  validateSystemIndex(system_index);
  const int map_index = permutor_map_indices.readHost(system_index);
  return (permutor_rotatable_group_bounds.readHost(map_index + 1) -
          permutor_rotatable_group_bounds.readHost(map_index));
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getSystemCisTransBondCount(const int system_index) const {
  validateSystemIndex(system_index);
  const int map_index = permutor_map_indices.readHost(system_index);
  return (permutor_cis_trans_group_bounds.readHost(map_index + 1) -
          permutor_cis_trans_group_bounds.readHost(map_index));
}

//-------------------------------------------------------------------------------------------------
int SynthesisPermutor::getSystemChiralCenterCount(const int system_index) const {
  validateSystemIndex(system_index);
  const int map_index = permutor_map_indices.readHost(system_index);
  return (permutor_invertible_group_bounds.readHost(map_index + 1) -
          permutor_invertible_group_bounds.readHost(map_index));
}

//-------------------------------------------------------------------------------------------------
const IsomerPlan& SynthesisPermutor::getRotatableGroup(const int system_index,
                                                       const int group_index) const {
  if (group_index < 0 || group_index >= getRotatableBondCount(system_index)) {
    const int tmp_map_index = permutor_map_indices.readHost(system_index);
    rtErr("Rotatable bond index " + std::to_string(group_index) + " is invalid for the system "
          "based on " + getBaseName(topologies[tmp_map_index]->getFileName()) + " with " +
          std::to_string(getRotatableBondCount(tmp_map_index)) + " rotatable bonds.",
          "SynthesisPermutor", "getRotatableGroup");
  }
  const int map_index = permutor_map_indices.readHost(system_index);
  return rotatable_groups[permutor_rotatable_group_bounds.readHost(map_index) + group_index];
}

//-------------------------------------------------------------------------------------------------
const IsomerPlan& SynthesisPermutor::getCisTransGroup(const int system_index,
                                                      const int group_index) const {
  if (group_index < 0 || group_index >= getCisTransBondCount(system_index)) {
    const int tmp_map_index = permutor_map_indices.readHost(system_index);
    rtErr("Cis-trans isomeric bond index " + std::to_string(group_index) + " is invalid for the "
          "system based on " + getBaseName(topologies[tmp_map_index]->getFileName()) + " with " +
          std::to_string(getCisTransBondCount(tmp_map_index)) + " cis-trans isomeric bonds.",
          "SynthesisPermutor", "getCisTransGroup");
  }
  const int map_index = permutor_map_indices.readHost(system_index);
  return cis_trans_groups[permutor_cis_trans_group_bounds.readHost(map_index) + group_index];
}

//-------------------------------------------------------------------------------------------------
const IsomerPlan& SynthesisPermutor::getInvertibleGroup(const int system_index,
                                                        const int group_index) const {
  if (group_index < 0 || group_index >= getChiralCenterCount(system_index)) {
    const int tmp_map_index = permutor_map_indices.readHost(system_index);
    rtErr("Chiral center number " + std::to_string(group_index) + " is invalid for the "
          "system based on " + getBaseName(topologies[tmp_map_index]->getFileName()) + " with " +
          std::to_string(getChiralCenterCount(tmp_map_index)) + " chiral centers.",
          "SynthesisPermutor", "getInvertibleGroup");
  }
  const int map_index = permutor_map_indices.readHost(system_index);
  return invertible_groups[permutor_invertible_group_bounds.readHost(map_index) + group_index];
}

//-------------------------------------------------------------------------------------------------
const TickCounter<double>& SynthesisPermutor::getStateTracker(const int system_index) const {
  validateSystemIndex(system_index);
  const int map_index = permutor_map_indices.readHost(system_index);
  return state_trackers[map_index];
}

//-------------------------------------------------------------------------------------------------
const TickCounter<double>& SynthesisPermutor::getStateTracker(const AtomGraph *ag) const {
  const int map_idx = matchTopology(ag, topologies);
  if (map_idx == permutor_map_count) {
    rtErr("No topology originating in file " + ag->getFileName() + " was found.",
          "SynthesisPermutor", "getStateTracker");
  }
  return state_trackers[map_idx];
}

//-------------------------------------------------------------------------------------------------
const TickCounter<double>& SynthesisPermutor::getStateTracker(const AtomGraph &ag) const {
  return getStateTracker(ag.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
llint SynthesisPermutor::getReplicaCount(const int query_index,
                                         const SamplingIntensity effort) const {
  validateSystemIndex(query_index);
  const int map_index = permutor_map_indices.readHost(query_index);
  switch (effort) {
  case SamplingIntensity::MINIMAL:
    return minimal_sampling_replicas[map_index];
  case SamplingIntensity::LIGHT:
    return light_sampling_replicas[map_index];
  case SamplingIntensity::HEAVY:
    return heavy_sampling_replicas[map_index];
  case SamplingIntensity::EXHAUSTIVE:
    return exhaustive_sampling_replicas[map_index];
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
llint SynthesisPermutor::getReplicaCount(const AtomGraph *query_ag,
                                         const SamplingIntensity effort) const {
  int map_idx = getPermutorMapIndex(query_ag);
  if (map_idx == permutor_map_count) {
    rtErr("A topology originating in file " + getBaseName(query_ag->getFileName()) + " was not "
          "matched to any in the synthesis.", "SynthesisPermutor", "getReplicaCount");
  }
  switch (effort) {
  case SamplingIntensity::MINIMAL:
    return minimal_sampling_replicas[map_idx];
  case SamplingIntensity::LIGHT:
    return light_sampling_replicas[map_idx];
  case SamplingIntensity::HEAVY:
    return heavy_sampling_replicas[map_idx];
  case SamplingIntensity::EXHAUSTIVE:
    return exhaustive_sampling_replicas[map_idx];
  }
  __builtin_unreachable();  
}

//-------------------------------------------------------------------------------------------------
llint SynthesisPermutor::getReplicaCount(const AtomGraph &query_ag,
                                         const SamplingIntensity effort) const {
  return getReplicaCount(query_ag.getSelfPointer(), effort);
}

//-------------------------------------------------------------------------------------------------
std::vector<llint> SynthesisPermutor::getReplicaCount(const SamplingIntensity effort) const {
  switch (effort) {
  case SamplingIntensity::MINIMAL:
    return std::vector<llint>(minimal_sampling_replicas.begin(), minimal_sampling_replicas.end());
  case SamplingIntensity::LIGHT:
    return std::vector<llint>(light_sampling_replicas.begin(), light_sampling_replicas.end());
  case SamplingIntensity::HEAVY:
    return std::vector<llint>(heavy_sampling_replicas.begin(), heavy_sampling_replicas.end());
  case SamplingIntensity::EXHAUSTIVE:
    return std::vector<llint>(exhaustive_sampling_replicas.begin(),
                              exhaustive_sampling_replicas.end());
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const ChemicalFeatures* SynthesisPermutor::getChemicalFeaturesPointer(const int map_index) const {
  if (features_objects_accessible) {
    rtErr("The underlying chemical features of each map may not be accessible any longer.",
          "SynthesisPointer", "getChemicalFeaturesPointer");
  }
  validateMapIndex(map_index);
  return features[map_index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<ChemicalFeatures*> SynthesisPermutor::getChemicalFeaturesPointer() const {
  if (features_objects_accessible) {
    rtErr("The underlying chemical features of each map may not be accessible any longer.",
          "SynthesisPointer", "getChemicalFeaturesPointer");
  }
  return features;
}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceSynthesis* SynthesisPermutor::getSynthesisPointer() const {
  return poly_ps_ptr;
}

//-------------------------------------------------------------------------------------------------
const SyPermutorKit<double> SynthesisPermutor::dpData(const HybridTargetLevel tier) const {
  return SyPermutorKit<double>(system_count, permutor_map_count, permutor_map_indices.data(tier),
                               permutor_elements.data(tier), permutor_element_bounds.data(tier),
                               system_settings.data(tier), system_settings_limits.data(tier),
                               rotatable_group_atoms.data(tier), rotatable_group_bounds.data(tier),
                               permutor_rotatable_group_bounds.data(tier),
                               cis_trans_group_atoms.data(tier), cis_trans_group_bounds.data(tier),
                               permutor_cis_trans_group_bounds.data(tier),
                               invertible_group_atoms.data(tier),
                               invertible_group_bounds.data(tier),
                               permutor_invertible_group_bounds.data(tier),
                               invertible_atom_centers.data(tier),
                               invertible_group_protocols.data(tier),
                               rotatable_bond_markers.data(tier),
                               cis_trans_bond_markers.data(tier), chiral_markers.data(tier),
                               rotatable_bond_settings.data(tier),
                               cis_trans_bond_settings.data(tier), 
                               rotatable_bond_settings_bounds.data(tier),
                               cis_trans_bond_settings_bounds.data(tier),
                               chiral_settings.data(tier), chiral_settings_bounds.data(tier));
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::upload() {
  permutor_elements.upload();
  synthesis_data.upload();
  rotatable_bond_settings.upload();
  cis_trans_bond_settings.upload();
  sp_rotatable_bond_settings.upload();
  sp_cis_trans_bond_settings.upload();
  group_data.upload();
  marker_data.upload();
  chiral_settings.upload();
  variables_data.upload();
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::download() {
  permutor_elements.download();
  synthesis_data.download();
  rotatable_bond_settings.download();
  cis_trans_bond_settings.download();
  sp_rotatable_bond_settings.download();
  sp_cis_trans_bond_settings.download();
  group_data.download();
  marker_data.download();
  chiral_settings.download();
  variables_data.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void
SynthesisPermutor::defineRotatableBondSettings(const int map_index,
                                               const std::vector<std::vector<double>> &settings) {
  alterMapSettings(map_index, &rotatable_bond_settings, rotatable_bond_settings_bounds,
                   permutor_rotatable_group_bounds, settings);
}

//-------------------------------------------------------------------------------------------------
void
SynthesisPermutor::defineCisTransBondSettings(const int map_index,
                                              const std::vector<std::vector<double>> &settings) {
  alterMapSettings(map_index, &cis_trans_bond_settings, cis_trans_bond_settings_bounds,
                   permutor_cis_trans_group_bounds, settings);
}

//-------------------------------------------------------------------------------------------------
void
SynthesisPermutor::defineChiralCenterSettings(const int map_index,
                                              const std::vector<std::vector<int>> &settings) {
  alterMapSettings(map_index, &chiral_settings, chiral_settings_bounds,
                   permutor_invertible_group_bounds, settings);
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::applySynthesis(const PhaseSpaceSynthesis *poly_ps_in,
                                       const VariableTorsionAdjustment adj) {
  poly_ps_ptr = const_cast<PhaseSpaceSynthesis*>(poly_ps_in);

  // Check that the synthesis is compatible with the available maps.
  const int ntop = poly_ps_ptr->getUniqueTopologyCount();
  if (ntop != permutor_map_count) {
    rtErr("A synthesis with " + std::to_string(ntop) + " unique topologies is incompatible with a "
          "collection of " + std::to_string(permutor_map_count) + " permutor maps.",
          "SynthesisPermutor", "applySynthesis");
  }
  const std::vector<AtomGraph*>& synth_tops = poly_ps_ptr->getUniqueTopologies();
  std::vector<int> topology_correspondence(ntop, 0);
  for (int i = 0; i < ntop; i++) {
    if (synth_tops[i] == topologies[i]) {
      topology_correspondence[i] = i;
    }
    else {
      topology_correspondence[i] = matchTopology(synth_tops[i], topologies);
    }
  }
  for (int i = 0; i < ntop; i++) {
    if (topology_correspondence[i] >= ntop || topology_correspondence[i] < 0) {
      rtErr("The synthesis is based, in part, on a topology originating in file " +
            getBaseName(synth_tops[i]->getFileName()) + ", which is nowhere to be found among the "
            "topologies used to create the permutor maps.", "SynthesisPermutor", "applySynthesis");
    }
  }

  // Allocate synthesis data and lay out maps of what systems correspond to each permutor map.
  PsSynthesisReader poly_psr = poly_ps_in->data();
  system_count = poly_psr.system_count;
  std::vector<int> tmp_permutor_map_indices(poly_psr.system_count);
  int total_var = 0;
  for (int i = 0; i < ntop; i++) {
    const int perm_map_idx = topology_correspondence[i];
    for (int j = poly_psr.common_ag_bounds[i]; j < poly_psr.common_ag_bounds[i + 1]; j++) {
      tmp_permutor_map_indices[poly_psr.common_ag_list[j]] = perm_map_idx;
    }
    const int ni_rotg = permutor_rotatable_group_bounds.readHost(perm_map_idx + 1) -
                        permutor_rotatable_group_bounds.readHost(perm_map_idx);
    const int ni_ctxg = permutor_cis_trans_group_bounds.readHost(perm_map_idx + 1) -
                        permutor_cis_trans_group_bounds.readHost(perm_map_idx);
    const int ni_invg = permutor_invertible_group_bounds.readHost(perm_map_idx + 1) -
                        permutor_invertible_group_bounds.readHost(perm_map_idx);
    total_var += (poly_psr.common_ag_bounds[i + 1] - poly_psr.common_ag_bounds[i]) *
                 (ni_rotg + ni_ctxg + ni_invg);
  }
  std::vector<int> tmp_system_settings_limits;
  tmp_system_settings_limits.reserve(total_var);
  for (int i = 0; i < poly_psr.system_count; i++) {
    const int perm_map_idx = topology_correspondence[poly_psr.unique_ag_idx[i]];
    const int rotg_llim = permutor_rotatable_group_bounds.readHost(perm_map_idx);
    const int rotg_hlim = permutor_rotatable_group_bounds.readHost(perm_map_idx + 1);
    const int ctxg_llim = permutor_cis_trans_group_bounds.readHost(perm_map_idx);
    const int ctxg_hlim = permutor_cis_trans_group_bounds.readHost(perm_map_idx + 1);
    const int invg_llim = permutor_invertible_group_bounds.readHost(perm_map_idx);
    const int invg_hlim = permutor_invertible_group_bounds.readHost(perm_map_idx + 1);
    const std::vector<int>& perm_map_ranges = state_trackers[perm_map_idx].getStateLimits();
    tmp_system_settings_limits.insert(tmp_system_settings_limits.end(), perm_map_ranges.begin(),
                                      perm_map_ranges.end());
  }
  std::vector<int> tmp_system_settings(total_var, 0);
  const size_t ns_data_size = roundUp(tmp_permutor_map_indices.size(), warp_size_zu) +
                              (2LLU * roundUp(tmp_system_settings.size(), warp_size_zu));
  synthesis_data.resize(ns_data_size);
  size_t ic = 0;
  ic = permutor_map_indices.putHost(&synthesis_data, tmp_permutor_map_indices, ic, warp_size_zu);
  ic = system_settings.putHost(&synthesis_data, tmp_system_settings, ic, warp_size_zu);
  ic = system_settings_limits.putHost(&synthesis_data, tmp_system_settings_limits, ic,
                                      warp_size_zu);
  
  // Find the highest number of systems in the synthesis and pre-allocate a vector to hold all of
  // the values for one particular rotatable bond.
  std::vector<double> observed_values;
  switch (adj) {
  case VariableTorsionAdjustment::ADJUST_NEARBY_VALUES:
  case VariableTorsionAdjustment::RESTRICT_TO_NEARBY_VALUES:
  case VariableTorsionAdjustment::CLUSTER_AND_APPLY_VALUES:
    {
      int max_reps = 0;
      for (int i = 0; i < poly_psr.unique_topology_count; i++) {
        max_reps = std::max(max_reps,
                            poly_psr.common_ag_bounds[i + 1] - poly_psr.common_ag_bounds[i]);
      }
      observed_values.resize(max_reps);
    }
    break;
  case VariableTorsionAdjustment::DO_NOT_CHANGE:

    // Return immediately.  Subsequent operations do not apply if the settings will not be altered
    // to conform to observations from the applied coordinate synthesis.
    return;
  }
  
  // Loop over all rotatable bonds, adjusting the available settings as required
  const int max_settings = std::max(std::max(rotatable_bond_samples, cis_trans_bond_samples), 2);
  ClusterManager<double, double> clsmgr(observed_values.size(), max_settings);
  int* rbsb_ptr = rotatable_bond_settings_bounds.data();
  double* rbs_ptr = rotatable_bond_settings.data();
  std::vector<bool> centroid_used(max_settings);
  for (int i = 0; i < ntop; i++) {
    const int perm_map_idx = topology_correspondence[i];
    const int nsharing = poly_psr.common_ag_bounds[i + 1] - poly_psr.common_ag_bounds[i];
    const int rbg_llim = permutor_rotatable_group_bounds.readHost(perm_map_idx);
    const int rbg_hlim = permutor_rotatable_group_bounds.readHost(perm_map_idx + 1);
    const int k_llim = poly_psr.common_ag_bounds[i];
    const int k_hlim = poly_psr.common_ag_bounds[i + 1];
    for (int j = rbg_llim; j < rbg_hlim; j++) {

      // Prepate the vector of torsion angle values and cluster it for the relevant variable.
      switch (adj) {
      case VariableTorsionAdjustment::ADJUST_NEARBY_VALUES:
      case VariableTorsionAdjustment::RESTRICT_TO_NEARBY_VALUES:
      case VariableTorsionAdjustment::CLUSTER_AND_APPLY_VALUES:
        {
          const int4 marks = rotatable_bond_markers.readHost(j);
          for (int k = k_llim; k < k_hlim; k++) {
            observed_values[k - k_llim] = dihedralAngle<double>(marks.x, marks.y, marks.z, marks.w,
                                                                poly_psr,
                                                                poly_psr.common_ag_list[k]);
          }

          // Cluster results as required
          const int nsettings = rbsb_ptr[j + 1] - rbsb_ptr[j];
          kMeans(&clsmgr, observed_values.data(), k_hlim - k_llim, nsettings, 1.0,
                 BoundaryCondition::PERIODIC, symbols::twopi);
        }
        break;
      case VariableTorsionAdjustment::DO_NOT_CHANGE:
        break;
      }

      // Apply the clustering in the prescribed manner.
      switch (adj) {
      case VariableTorsionAdjustment::ADJUST_NEARBY_VALUES:
      case VariableTorsionAdjustment::RESTRICT_TO_NEARBY_VALUES:
        {
          // Loop over all settings and, if they lie within a tolerance of one of the centroids,
          // adjust the setting to meet the centroid.
          const int nsettings = rbsb_ptr[j + 1] - rbsb_ptr[j];
          for (int k = 0; k < nsettings; k++) {
            centroid_used[k] = false;
          }
          for (int k = rbsb_ptr[j]; k < rbsb_ptr[j + 1]; k++) {
            double nearest_centroid_dist = symbols::twopi;
            int nearest_centroid = -1;
            for (int m = 0; m < nsettings; m++) {
              if (centroid_used[m] == false) {
                const double centroid_dist = imageValue(clsmgr.getClusterCentroid(m) - rbs_ptr[k],
                                                        symbols::twopi,
                                                        ImagingMethod::MINIMUM_IMAGE);
                if (centroid_dist < nearest_centroid_dist) {
                  nearest_centroid = m;
                  nearest_centroid_dist = centroid_dist;
                }
              }
            }
            if (nearest_centroid >= 0 &&
                (nearest_centroid_dist < rotation_setting_snap ||
                 adj == VariableTorsionAdjustment::RESTRICT_TO_NEARBY_VALUES)) {
              centroid_used[nearest_centroid] = true;
              rbs_ptr[k] = clsmgr.getClusterCentroid(nearest_centroid);
            }
          }
        }
        break;
      case VariableTorsionAdjustment::CLUSTER_AND_APPLY_VALUES:
        for (int k = rbsb_ptr[j]; k < rbsb_ptr[j + 1]; k++) {
          rbs_ptr[k] = clsmgr.getClusterCentroid(k - rbsb_ptr[j]);
        }
        break;
      case VariableTorsionAdjustment::DO_NOT_CHANGE:
        break;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::applySynthesis(const PhaseSpaceSynthesis &poly_ps_in,
                                       const VariableTorsionAdjustment adj) {
  applySynthesis(poly_ps_in.getSelfPointer(), adj);
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::setVariableRanges(const std::vector<double> &general_rot_settings_in,
                                          const std::vector<double> &general_ctx_settings_in,
                                          const bool sample_cis_trans_states,
                                          const bool sample_chiral_states) {

  // Obtain general settings for rotatable and cis-trans isomeric bonds from user input.
  general_rotation_settings = general_rot_settings_in;
  general_cis_trans_settings = general_ctx_settings_in;
  rotatable_bond_samples = general_rotation_settings.size();
  cis_trans_bond_samples = general_cis_trans_settings.size();

  // Create the state trackers for each permutor map.
  state_trackers.reserve(permutor_map_count);
  const int* prgb_ptr = permutor_rotatable_group_bounds.data();
  const int* pcgb_ptr = permutor_cis_trans_group_bounds.data();
  const int* pigb_ptr = permutor_invertible_group_bounds.data();
  std::vector<double> tmp_rotatable_bond_settings, tmp_cis_trans_bond_settings;
  std::vector<int> tmp_chiral_settings;
  const int n_rotg = permutor_rotatable_group_bounds.readHost(permutor_map_count);
  const int n_ctxg = permutor_cis_trans_group_bounds.readHost(permutor_map_count);
  const int n_invg = permutor_invertible_group_bounds.readHost(permutor_map_count);
  std::vector<int> tmp_rotatable_bond_settings_bounds(n_rotg + 1, 0);
  std::vector<int> tmp_cis_trans_bond_settings_bounds(n_ctxg + 1, 0);
  std::vector<int> tmp_chiral_settings_bounds(n_invg + 1, 0);
  tmp_rotatable_bond_settings.reserve(rotatable_bond_samples * n_rotg);
  tmp_cis_trans_bond_settings.reserve(cis_trans_bond_samples * n_ctxg);
  
  // Lay out each group's settings one by one, to allow for some groups to have different numbers
  // of settings based on user input.
  int rotg_counter = 0;
  int ctxg_counter = 0;
  int invg_counter = 0;
  int rotg_value_counter = 0;
  int ctxg_value_counter = 0;
  int invg_value_counter = 0;
  for (int i = 0; i < permutor_map_count; i++) {
    const int ni_rotg = prgb_ptr[i + 1] - prgb_ptr[i];
    const int ni_ctxg = pcgb_ptr[i + 1] - pcgb_ptr[i];
    const int ni_invg = pigb_ptr[i + 1] - pigb_ptr[i];
    const int ni_var = ni_rotg + ni_ctxg + ni_invg;
    std::vector<std::vector<double>> tracker_settings(ni_var, std::vector<double>());
    int var_counter = 0;
    for (int j = 0; j < ni_rotg; j++) {
      tmp_rotatable_bond_settings_bounds[rotg_counter] = rotatable_bond_samples;
      rotg_counter++;
      for (int k = 0; k < rotatable_bond_samples; k++) {
        tmp_rotatable_bond_settings.push_back(general_rotation_settings[k]);
      }
      tracker_settings[var_counter] = general_rotation_settings;
      var_counter++;
    }
    for (int j = 0; j < ni_ctxg; j++) {
      tmp_cis_trans_bond_settings_bounds[ctxg_counter] = cis_trans_bond_samples;
      ctxg_counter++;
      for (int k = 0; k < cis_trans_bond_samples; k++) {
        tmp_cis_trans_bond_settings.push_back(general_cis_trans_settings[k]);
      }
      tracker_settings[var_counter] = general_cis_trans_settings;
      var_counter++;
    }
    for (int j = 0; j < ni_invg; j++) {
      switch (invertible_groups[invg_counter].getChiralPlan()) {
      case ChiralInversionProtocol::ROTATE:
      case ChiralInversionProtocol::REFLECT:
        tmp_chiral_settings_bounds[invg_counter] = 2;
        tmp_chiral_settings.push_back(static_cast<int>(ChiralOrientation::RECTUS));
        tmp_chiral_settings.push_back(static_cast<int>(ChiralOrientation::SINISTER));
        tracker_settings[var_counter] = { 0.0, 1.0 };
        break;
      case ChiralInversionProtocol::DO_NOT_INVERT:
        tmp_chiral_settings_bounds[invg_counter] = 1;
        tmp_chiral_settings.push_back(static_cast<int>(ChiralOrientation::NONE));
        tracker_settings[var_counter] = { 0.0 };
        break;
      }
      invg_counter++;
      var_counter++;
    }
    state_trackers.emplace_back(tracker_settings);
  }
  prefixSumInPlace(&tmp_rotatable_bond_settings_bounds, PrefixSumType::EXCLUSIVE,
                   "sum the temporary array of rotatable bond settings bounds");
  prefixSumInPlace(&tmp_cis_trans_bond_settings_bounds, PrefixSumType::EXCLUSIVE,
                   "sum the temporary array of cis-trans isomeric bond settings bounds");
  prefixSumInPlace(&tmp_chiral_settings_bounds, PrefixSumType::EXCLUSIVE,
                   "sum the temporary array of chiral center settings bounds");
  
  // Given the completed permutation counters, compute the sampling needed for minimal, light,
  // heavy, and exhaustive sampling.  While computing these samplings, store the coupled variables
  // for rotatable bonds, cis-trans isomeric bonds, and invertible chiral centers.
  minimal_sampling_replicas.resize(permutor_map_count);
  light_sampling_replicas.resize(permutor_map_count);
  heavy_sampling_replicas.resize(permutor_map_count);
  exhaustive_sampling_replicas.resize(permutor_map_count);
  for (int i = 0; i < permutor_map_count; i++) {
    const int ni_var = state_trackers[i].getVariableCount();
    const int ni_rotg = prgb_ptr[i + 1] - prgb_ptr[i];
    const int ni_ctxg = pcgb_ptr[i + 1] - pcgb_ptr[i];
    const int ni_invg = pigb_ptr[i + 1] - pigb_ptr[i];
    std::vector<IsomerPlan> all_mutables;
    all_mutables.reserve(ni_var);
    for (int j = prgb_ptr[i]; j < prgb_ptr[i + 1]; j++) {
      all_mutables.push_back(rotatable_groups[j]);
    }
    for (int j = pcgb_ptr[i]; j < pcgb_ptr[i + 1]; j++) {
      all_mutables.push_back(cis_trans_groups[j]);
    }
    for (int j = pigb_ptr[i]; j < pigb_ptr[i + 1]; j++) {
      all_mutables.push_back(invertible_groups[j]);
    }
    const NonbondedKit<double> nbk = topologies[i]->getDoublePrecisionNonbondedKit();
    std::vector<bool> var_linked(ni_var * ni_var, true);
    for (int j = 1; j < ni_var; j++) {
      for (int k = 0; k < j; k++) {
        const bool jk_link = permutationsAreLinked(all_mutables, j, k, nbk);
        var_linked[(j * ni_var) + k] = jk_link;
        var_linked[(k * ni_var) + j] = jk_link;
      }
    }
    int nminml = 0;
    int nlight = 0;
    int nheavy = 0;
    for (int j = 0; j < ni_var; j++) {
      const int nsj = state_trackers[i].getStateLimits(j);
      nminml += (nsj > 1) * nsj;
      bool one_linked = false;
      for (int k = 0; k < j; k++) {
        const int nsk = state_trackers[i].getStateLimits(k);
        if (var_linked[(j * ni_var) + k]) {
          if (nsj > 1 || nsk > 1) {
            nlight += nsj * nsk;
            one_linked = true;
          }
          bool two_linked = false;
          for (int m = 0; m < k; m++) {
            const int nsm = state_trackers[i].getStateLimits(m);
            if ((nsj > 1 || nsk > 1 || nsm > 1) && var_linked[(k * ni_var) + m]) {
              nheavy += nsj * nsk * nsm;
              two_linked = true;
            }
          }
          if (two_linked == false && (nsj > 1 || nsk > 1)) {
            nheavy += nsj * nsk;
          }
        }
      }
      if (one_linked == false) {
        nlight += (nsj > 1) * nsj;
        nheavy += (nsj > 1) * nsj;
      }
    }
    minimal_sampling_replicas[i] = nminml + (nminml == 0);
    light_sampling_replicas[i] = nlight + (nlight == 0);
    heavy_sampling_replicas[i] = nheavy +  (nheavy == 0);
    
    // Cap the number of replicas for exhaustive sampling at approximately 2^60, a huge number
    // that would only be calculable with supercomputing resources and impractical to store, ever.
    // As with other orders of sampling, ensure that there is at least one conformation allocated
    // for the exhaustive sampling case, even if there are no variable elements in the system.
    if (state_trackers[i].getLogPermutations() < 60.0) {
      exhaustive_sampling_replicas[i] = state_trackers[i].getExactPermutationCount();
    }
    else {
      exhaustive_sampling_replicas[i] = 1.153e18;
    }
    exhaustive_sampling_replicas[i] += (exhaustive_sampling_replicas[i] == 0LL);
  }
  
  // Load the appropriate Hybrid objects
  rotatable_bond_settings.resize(tmp_rotatable_bond_settings.size());
  rotatable_bond_settings.putHost(tmp_rotatable_bond_settings);
  rotatable_bond_settings_bounds.putHost(tmp_rotatable_bond_settings_bounds);
  cis_trans_bond_settings.resize(tmp_cis_trans_bond_settings.size());
  cis_trans_bond_settings.putHost(tmp_cis_trans_bond_settings);
  cis_trans_bond_settings_bounds.putHost(tmp_cis_trans_bond_settings_bounds);
  const std::vector<float> tmp_sp_rotatable_bond_settings(tmp_rotatable_bond_settings.begin(),
                                                          tmp_rotatable_bond_settings.end());
  const std::vector<float> tmp_sp_cis_trans_bond_settings(tmp_cis_trans_bond_settings.begin(),
                                                          tmp_cis_trans_bond_settings.end());
  sp_rotatable_bond_settings.resize(tmp_sp_rotatable_bond_settings.size());
  sp_rotatable_bond_settings.putHost(tmp_sp_rotatable_bond_settings);
  sp_cis_trans_bond_settings.resize(tmp_sp_cis_trans_bond_settings.size());
  sp_cis_trans_bond_settings.putHost(tmp_sp_cis_trans_bond_settings);
  chiral_settings.resize(tmp_chiral_settings.size());
  chiral_settings.putHost(tmp_chiral_settings);
  chiral_settings_bounds.putHost(tmp_chiral_settings_bounds);
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::impartControlData(const ConformerControls &confcon) {
  setVariableRanges(confcon.getRotationSampleValues(), confcon.getCisTransSampleValues(),
                    confcon.sampleCisTrans(), confcon.sampleChirality());
  rotation_setting_snap = confcon.getRotatableBondSnapThreshold();
  cis_trans_setting_snap = confcon.getCisTransBondSnapThreshold();
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::permuteSystem(PhaseSpaceSynthesis *psynth, const int system_index,
                                      const int map_index, const CoupledEdit ce,
                                      const int setting_index, const PrecisionModel prec) const {

  // Check the validity of the request.
  if (system_index < 0 || system_index >= psynth->getSystemCount()) {
    rtErr("System index " + std::to_string(system_index) + " is out of bounds for a collection "
          "of " + std::to_string(psynth->getSystemCount()) + " systems.", "SynthesisPermutor",
          "permuteSystem");
  }
  const int natom = topologies[map_index]->getAtomCount();
  const UnitCellType uc = topologies[map_index]->getUnitCellType();
  CoordinateSeries<double> dcs(natom, 0, uc);
  CoordinateSeries<float> fcs(natom, 0, uc);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    dcs.resize(1);
    break;
  case PrecisionModel::SINGLE:
    fcs.resize(1);
    break;
  }
  switch (prec) {
  case PrecisionModel::DOUBLE:
    coordCopy<double>(&dcs, 0, *psynth, system_index);
    permuteSystem<double>(&dcs, map_index, ce, setting_index);
    coordCopy<double>(psynth, system_index, dcs, 0);
    break;
  case PrecisionModel::SINGLE:
    coordCopy<float>(&fcs, 0, *psynth, system_index);
    permuteSystem<float>(&fcs, map_index, ce, setting_index);
    coordCopy<float>(psynth, system_index, fcs, 0);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::permuteSystem(const int system_index, const CoupledEdit ce,
                                      const int setting_index, const PrecisionModel prec) {
  const int map_index = permutor_map_indices.readHost(system_index);
  permuteSystem(poly_ps_ptr, system_index, map_index, ce, setting_index, prec);
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::randomizeSystem(PhaseSpaceSynthesis *psynth, const int system_index,
                                        const int map_index, Xoshiro256ppGenerator *xrs,
                                        const PrecisionModel prec, const CoupledEdit mask_a,
                                        const CoupledEdit mask_b, const CoupledEdit mask_c) const {

  // Check the validity of the request.
  if (system_index < 0 || system_index >= psynth->getSystemCount()) {
    rtErr("System index " + std::to_string(system_index) + " is out of bounds for a collection "
          "of " + std::to_string(psynth->getSystemCount()) + " systems.", "SynthesisPermutor",
          "permuteSystem");
  }
  const int natom = topologies[map_index]->getAtomCount();
  const UnitCellType uc = topologies[map_index]->getUnitCellType();
  CoordinateSeries<double> dcs(natom, 0, uc);
  CoordinateSeries<float> fcs(natom, 0, uc);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    dcs.resize(1);
    break;
  case PrecisionModel::SINGLE:
    fcs.resize(1);
    break;
  }
  switch (prec) {
  case PrecisionModel::DOUBLE:
    coordCopy<double>(&dcs, 0, *psynth, system_index);
    randomizeSystem<double>(&dcs, map_index, xrs, mask_a, mask_b, mask_c);
    coordCopy<double>(psynth, system_index, dcs, 0);
    break;
  case PrecisionModel::SINGLE:
    coordCopy<float>(&fcs, 0, *psynth, system_index);
    randomizeSystem<float>(&fcs, map_index, xrs, mask_a, mask_b, mask_c);
    coordCopy<float>(psynth, system_index, fcs, 0);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::randomizeSystem(const int system_index, Xoshiro256ppGenerator *xrs,
                                        const PrecisionModel prec, const CoupledEdit mask_a,
                                        const CoupledEdit mask_b, const CoupledEdit mask_c) {
  randomizeSystem(poly_ps_ptr, system_index, permutor_map_indices.readHost(system_index), xrs,
                  prec, mask_a, mask_b, mask_c);
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis
SynthesisPermutor::buildSynthesis(const SamplingIntensity effort, Xoshiro256ppGenerator *xrs,
                                  std::vector<int> *correspondence, const int system_limit,
                                  const int globalpos_scale_bits, const int localpos_scale_bits,
                                  const int velocity_scale_bits, const int force_scale_bits,
                                  const PrecisionModel prec, ClashReport *clrep,
                                  const int iteration_limit, const int max_clashes) const {
  size_t sys_requirement = 0;
  const int* perm_idx_ptr = permutor_map_indices.data();
  std::vector<SamplingIntensity> per_system_effort(system_count, effort);
  std::vector<int> per_system_replicas(system_count + 1, 0);
  llint sum_reps;
  do {
    for (int i = 0; i < system_count; i++) {
      const int perm_map_idx = perm_idx_ptr[i];
      switch (per_system_effort[i]) {
      case SamplingIntensity::MINIMAL:
        per_system_replicas[i] = minimal_sampling_replicas[perm_map_idx];
        break;
      case SamplingIntensity::LIGHT:
        per_system_replicas[i] = light_sampling_replicas[perm_map_idx];
        break;
      case SamplingIntensity::HEAVY:
        per_system_replicas[i] = heavy_sampling_replicas[perm_map_idx];
        break;
      case SamplingIntensity::EXHAUSTIVE:
        if (exhaustive_sampling_replicas[i] >= static_cast<llint>(INT_MAX)) {
          per_system_effort[i] = SamplingIntensity::HEAVY;
          per_system_replicas[i] = heavy_sampling_replicas[perm_map_idx];
        }
        else {
          per_system_replicas[i] = exhaustive_sampling_replicas[perm_map_idx];
        }
        break;
      }
    }
    sum_reps = sum<llint>(per_system_replicas);

    // If the system limit has been exceeded, find the steepest way to reduce the counts.
    if (sum_reps > system_limit) {
      int best_reduction = -1;
      int best_spot, best_map_idx;
      for (int i = 0; i < system_count; i++) {
        const int perm_map_idx = perm_idx_ptr[i];
        int reduction;
        switch (per_system_effort[i]) {
        case SamplingIntensity::MINIMAL:
          reduction = 0;
          break;
        case SamplingIntensity::LIGHT:
          reduction = light_sampling_replicas[perm_map_idx] -
                      minimal_sampling_replicas[perm_map_idx];
          break;
        case SamplingIntensity::HEAVY:
          reduction = heavy_sampling_replicas[perm_map_idx] -
                      light_sampling_replicas[perm_map_idx];
          break;
        case SamplingIntensity::EXHAUSTIVE:
          reduction = exhaustive_sampling_replicas[perm_map_idx] -
                      heavy_sampling_replicas[perm_map_idx];
          break;
        }
        if (reduction > best_reduction) {
          best_reduction = reduction;
          best_spot = i;
          best_map_idx = perm_map_idx;
        }
      }
      if (best_reduction == -1) {
        rtErr("The number of systems in the resulting synthesis cannot be reduced below the "
              "target of " + std::to_string(system_limit) + ".", "SynthesisPermutor",
              "buildSynthesis");
      }
      switch (per_system_effort[best_spot]) {
      case SamplingIntensity::MINIMAL:
        break;
      case SamplingIntensity::LIGHT:
        per_system_effort[best_spot] = SamplingIntensity::MINIMAL;
        per_system_replicas[best_spot] = minimal_sampling_replicas[best_map_idx];
        break;
      case SamplingIntensity::HEAVY:
        per_system_effort[best_spot] = SamplingIntensity::LIGHT;
        per_system_replicas[best_spot] = light_sampling_replicas[best_map_idx];
        break;
      case SamplingIntensity::EXHAUSTIVE:
        per_system_effort[best_spot] = SamplingIntensity::HEAVY;
        per_system_replicas[best_spot] = heavy_sampling_replicas[best_map_idx];
        break;
      }
    }
  } while (sum_reps > system_limit);

  // The design entails a bounds array on the replicas for various systems of the current synthesis
  // in the "x" members, and indicators of what degree of sampling to perform in the "y" members.
  prefixSumInPlace(&per_system_replicas, PrefixSumType::EXCLUSIVE);
  Hybrid<int2> design(system_count + 1, "new_synth_map");
  int2* design_ptr = design.data();
  for (int i = 0; i < system_count; i++) {
    design_ptr[i].x = per_system_replicas[i];
    design_ptr[i].y = static_cast<int>(per_system_effort[i]);
  }
  design_ptr[system_count].x = per_system_replicas[system_count];
  design_ptr[system_count].y = static_cast<int>(SamplingIntensity::MINIMAL);
  Hybrid<uint> placement_success((per_system_replicas[system_count] + warp_bits_mask_int) /
                                 warp_size_int, "success_mask");
  uint* plc_success_ptr = placement_success.data();

#ifdef STORMM_USE_HPC
  // Upload the design so that the GPU can perform placements if possible
  design.upload();
#endif
  
  // Create the resulting synthesis
  const std::vector<AtomGraph*>& orig_tops = poly_ps_ptr->getUniqueTopologies();
  std::vector<PhaseSpace> orig_coords;
  orig_coords.reserve(system_count);
  std::vector<int> coord_index_key(sum_reps), topol_index_key(sum_reps);
  int result_tracker = 0;
  for (int i = 0; i < system_count; i++) {
    orig_coords.push_back(poly_ps_ptr->exportSystem(i));
    const int orig_topol_idx = poly_ps_ptr->getUniqueTopologyIndex(i);
    for (int j = per_system_replicas[i]; j < per_system_replicas[i + 1]; j++) {
      coord_index_key[result_tracker] = i;
      topol_index_key[result_tracker] = orig_topol_idx;
      result_tracker++;
    }
  }
  PhaseSpaceSynthesis trial(orig_coords, coord_index_key, orig_tops, topol_index_key,
                            globalpos_scale_bits, localpos_scale_bits, velocity_scale_bits,
                            force_scale_bits);
  
  // The resulting synthesis is now filled with copies of the synthesis referenced by this
  // SynthesisPermutor object.  Perform the permutations on the new synthesis.
  const bool use_single_prec = (globalpos_scale_bits <= globalpos_scale_nonoverflow_bits);
//#ifdef STORMM_USE_HPC

  // Count the number of successful placements made on the GPU
  //int n_success = 0;
//#else
  int n_success = 0;
  for (int i = 0; i < system_count; i++) {
    const int imap_index = permutor_map_indices.readHost(i);
    if (clrep != nullptr) {
      clrep->setTopologyPointer(topologies[imap_index]);
    }
    const StaticExclusionMask imask(topologies[imap_index]);
    
    // Allocate space in the proper format CoordinateSeries<T> object for a floating-point
    // representation of the system that can help minimize costly fixed-precision conversions and
    // improve the replication of the process that will take place on the GPU.
    const int natom = topologies[imap_index]->getAtomCount();
    const UnitCellType uc = topologies[imap_index]->getUnitCellType();
    CoordinateSeries<double> dcs(natom, 0, uc);
    CoordinateSeries<float> fcs(natom, 0, uc);
    switch (prec) {
    case PrecisionModel::DOUBLE:
      dcs.resize(1);
      break;
    case PrecisionModel::SINGLE:
      fcs.resize(1);
      break;
    }
    int pos = design_ptr[i].x;
    switch (static_cast<SamplingIntensity>(design_ptr[i].y)) {
    case SamplingIntensity::MINIMAL:
      {
        // Loop over all mutable elements of the system.
        const int rotg_llim = permutor_rotatable_group_bounds.readHost(imap_index);
        const int rotg_hlim = permutor_rotatable_group_bounds.readHost(imap_index + 1);
        for (int j = rotg_llim; j < rotg_hlim; j++) {
          const int k_llim = rotatable_bond_settings_bounds.readHost(j);
          const int k_hlim = rotatable_bond_settings_bounds.readHost(j + 1);
          const CoupledEdit focus = { ConformationEdit::BOND_ROTATION, j };
          if (k_hlim > k_llim + 1) {
            for (int k = k_llim; k < k_hlim; k++) {
              if (pos >= design_ptr[i + 1].x) {
                rtErr("The position counter for system " + std::to_string(i) + ", with " +
                      std::to_string(k_hlim - k_llim) + " rotatable bond settings, has exceeded "
                      "the number of replicas (" + std::to_string(pos) + " >= " +
                      std::to_string(design_ptr[i + 1].x) + ").", "SynthesisPermutor",
                      "buildSynthesis");
              }
              switch (prec) {
              case PrecisionModel::DOUBLE:
                coordCopy<double>(&dcs, 0, trial, pos);
                randomizeSystem<double>(&dcs, imap_index, xrs, focus);
                permuteSystem<double>(&dcs, imap_index, focus, k - k_llim);
                if (resolveClashes<double>(&dcs, imap_index, xrs, imask, clrep, iteration_limit,
                                           max_clashes, focus)) {
                  accumulateBitmask(plc_success_ptr, pos);
                  n_success++;
                  coordCopy<double>(&trial, pos, dcs, 0);
                }
                break;
              case PrecisionModel::SINGLE:
                coordCopy<float>(&fcs, 0, trial, pos);
                randomizeSystem<float>(&fcs, imap_index, xrs, focus);
                permuteSystem<float>(&fcs, imap_index, focus, k - k_llim);
                if (resolveClashes<float>(&fcs, imap_index, xrs, imask, clrep, iteration_limit,
                                          max_clashes, focus)) {
                  accumulateBitmask(plc_success_ptr, pos);
                  n_success++;
                  coordCopy<float>(&trial, pos, fcs, 0);
                }
                break;
              }
              pos++;
            }
          }
        }
        const int ctxg_llim = permutor_cis_trans_group_bounds.readHost(imap_index);
        const int ctxg_hlim = permutor_cis_trans_group_bounds.readHost(imap_index + 1);
        for (int j = ctxg_llim; j < ctxg_hlim; j++) {
          const int k_llim = cis_trans_bond_settings_bounds.readHost(j);
          const int k_hlim = cis_trans_bond_settings_bounds.readHost(j + 1);
          const CoupledEdit focus = { ConformationEdit::CIS_TRANS_FLIP, j };
          if (k_hlim > k_llim + 1) {
            for (int k = k_llim; k < k_hlim; k++) {
              if (pos >= design_ptr[i + 1].x) {
                rtErr("The position counter for system " + std::to_string(i) + ", with " +
                      std::to_string(k_hlim - k_llim) + " cis-trans isomeric bond settings, has "
                      "exceeded the number of replicas (" + std::to_string(pos) + " >= " +
                      std::to_string(design_ptr[i + 1].x) + ").", "SynthesisPermutor",
                      "buildSynthesis");
              }
              switch (prec) {
              case PrecisionModel::DOUBLE:
                coordCopy<double>(&dcs, 0, trial, pos);
                randomizeSystem<double>(&dcs, imap_index, xrs, focus);
                permuteSystem<double>(&dcs, imap_index, focus, k - k_llim);
                if (resolveClashes<double>(&dcs, imap_index, xrs, imask, clrep, iteration_limit,
                                           max_clashes, focus)) {
                  accumulateBitmask(plc_success_ptr, pos);
                  n_success++;
                  coordCopy<double>(&trial, pos, dcs, 0);
                }
                break;
              case PrecisionModel::SINGLE:
                coordCopy<float>(&fcs, 0, trial, pos);
                randomizeSystem<float>(&fcs, imap_index, xrs, focus);
                permuteSystem<float>(&fcs, imap_index, focus, k - k_llim);
                if (resolveClashes<float>(&fcs, imap_index, xrs, imask, clrep, iteration_limit,
                                          max_clashes, focus)) {
                  accumulateBitmask(plc_success_ptr, pos);
                  n_success++;
                  coordCopy<float>(&trial, pos, fcs, 0);
                }
                break;
              }
              pos++;
            }
          }
        }
        const int invg_llim = permutor_invertible_group_bounds.readHost(imap_index);
        const int invg_hlim = permutor_invertible_group_bounds.readHost(imap_index + 1);
        for (int j = invg_llim; j < invg_hlim; j++) {
          const int k_llim = chiral_settings_bounds.readHost(j);
          const int k_hlim = chiral_settings_bounds.readHost(j + 1);
          const CoupledEdit focus = { ConformationEdit::CHIRAL_INVERSION, j };
          if (k_hlim > k_llim + 1) {
            for (int k = k_llim; k < k_hlim; k++) {
              if (pos >= design_ptr[i + 1].x) {
                rtErr("The position counter for system " + std::to_string(i) + ", with " +
                      std::to_string(k_hlim - k_llim) + " chiral settings, has exceeded the "
                      "number of replicas (" + std::to_string(pos) + " >= " +
                      std::to_string(design_ptr[i + 1].x) + ").", "SynthesisPermutor",
                      "buildSynthesis");
              }
              switch (prec) {
              case PrecisionModel::DOUBLE:
                coordCopy<double>(&dcs, 0, trial, pos);
                randomizeSystem<double>(&dcs, imap_index, xrs, focus);
                permuteSystem<double>(&dcs, imap_index, focus, k - k_llim);
                if (resolveClashes<double>(&dcs, imap_index, xrs, imask, clrep, iteration_limit,
                                           max_clashes, focus)) {
                  accumulateBitmask(plc_success_ptr, pos);
                  n_success++;
                  coordCopy<double>(&trial, pos, dcs, 0);
                }
                break;
              case PrecisionModel::SINGLE:
                coordCopy<float>(&fcs, 0, trial, pos);
                randomizeSystem<float>(&fcs, imap_index, xrs, focus);
                permuteSystem<float>(&fcs, imap_index, focus, k - k_llim);
                if (resolveClashes<float>(&fcs, imap_index, xrs, imask, clrep, iteration_limit,
                                          max_clashes, focus)) {
                  accumulateBitmask(plc_success_ptr, pos);
                  n_success++;
                  coordCopy<float>(&trial, pos, fcs, 0);
                }
                break;
              }
              pos++;
            }
          }
        }
      }
      break;
    case SamplingIntensity::LIGHT:
      {
        const int j_llim = light_sampling_foci_bounds.readHost(imap_index);
        const int j_hlim = light_sampling_foci_bounds.readHost(imap_index + 1);
        for (int j = j_llim; j < j_hlim; j += 2) {

          // Get the two variables and their types
          const CoupledEdit vara_raw(light_sampling_foci.readHost(j));
          const CoupledEdit varb_raw(light_sampling_foci.readHost(j + 1));

          // Obtain the number of settings available to each variable.  Skip if both variables
          // have only one available setting.
          const int vara_settings = getElementSampleCount(vara_raw);
          const int varb_settings = getElementSampleCount(varb_raw);
          if (vara_settings * varb_settings == 1) {
            continue;
          }
          for (int k = 0; k < vara_settings; k++) {
            for (int m = 0; m < varb_settings; m++) {
              if (pos >= design_ptr[i + 1].x) {
                rtErr("The position counter for system " + std::to_string(i) + ", with " +
                      std::to_string(j_hlim - j_llim) + " pairs of coupled mutable elements, "
                      "has exceeded the number of replicas with pair " + std::to_string(j / 2) +
                      "(" + std::to_string(k) + ", " + std::to_string(m) + ").",
                      "SynthesisPermutor", "buildSynthesis");
              }
              switch (prec) {
              case PrecisionModel::DOUBLE:
                coordCopy<double>(&dcs, 0, trial, pos);
                randomizeSystem<double>(&dcs, imap_index, xrs, vara_raw, varb_raw);
                permuteSystem<double>(&dcs, imap_index, vara_raw, k);
                if (varb_raw.index >= 0) {
                  permuteSystem<double>(&dcs, imap_index, varb_raw, m);
                }
                if (resolveClashes<double>(&dcs, imap_index, xrs, imask, clrep, iteration_limit,
                                           max_clashes, vara_raw, varb_raw)) {
                  accumulateBitmask(plc_success_ptr, pos);
                  n_success++;
                  coordCopy<double>(&trial, pos, dcs, 0);
                }
                break;
              case PrecisionModel::SINGLE:
                coordCopy<float>(&fcs, 0, trial, pos);
                randomizeSystem<float>(&fcs, imap_index, xrs, vara_raw, varb_raw);
                permuteSystem<float>(&fcs, imap_index, vara_raw, k);
                if (varb_raw.index >= 0) {
                  permuteSystem<float>(&fcs, imap_index, varb_raw, m);
                }
                if (resolveClashes<float>(&fcs, imap_index, xrs, imask, clrep, iteration_limit,
                                          max_clashes, vara_raw, varb_raw)) {
                  accumulateBitmask(plc_success_ptr, pos);
                  n_success++;
                  coordCopy<float>(&trial, pos, fcs, 0);
                }
                break;
              }
              pos++;
            }
          }
        }
      }
      break;
    case SamplingIntensity::HEAVY:
      {
        const int j_llim = heavy_sampling_foci_bounds.readHost(imap_index);
        const int j_hlim = heavy_sampling_foci_bounds.readHost(imap_index + 1);
        for (int j = j_llim; j < j_hlim; j += 3) {
        
          // Get the three variables and their types
          const CoupledEdit vara_raw(heavy_sampling_foci.readHost(j));
          const CoupledEdit varb_raw(heavy_sampling_foci.readHost(j + 1));
          const CoupledEdit varc_raw(heavy_sampling_foci.readHost(j + 2));

          // Obtain the number of settings available to each variable.  Skip if all three variables
          // have only one available setting.
          const int vara_settings = getElementSampleCount(vara_raw);
          const int varb_settings = getElementSampleCount(varb_raw);
          const int varc_settings = getElementSampleCount(varc_raw);
          if (vara_settings * varb_settings * varc_settings == 1) {
            continue;
          }
          for (int k = 0; k < vara_settings; k++) {
            for (int m = 0; m < varb_settings; m++) {
              for (int n = 0; n < varc_settings; n++) {
                if (pos >= design_ptr[i + 1].x) {
                  rtErr("The position counter for system " + std::to_string(i) + ", with " +
                        std::to_string(j_hlim - j_llim) + " triplets of coupled mutable "
                        "elements, has exceeded the number of replicas with triplet " +
                        std::to_string(j / 3) + "(" + std::to_string(k) + ", " +
                        std::to_string(m) + ", " + std::to_string(n) + ").", "SynthesisPermutor",
                        "buildSynthesis");
                }
                switch (prec) {
                case PrecisionModel::DOUBLE:
                  coordCopy<double>(&dcs, 0, trial, pos);
                  randomizeSystem<double>(&dcs, imap_index, xrs, vara_raw, varb_raw, varc_raw);
                  permuteSystem<double>(&dcs, imap_index, vara_raw, k);
                  if (varb_raw.index >= 0) {
                    permuteSystem<double>(&dcs, imap_index, varb_raw, m);
                  }
                  if (varc_raw.index >= 0) {
                    permuteSystem<double>(&dcs, imap_index, varc_raw, n);
                  }
                  if (resolveClashes<double>(&dcs, imap_index, xrs, imask, clrep, iteration_limit,
                                             max_clashes, vara_raw, varb_raw, varc_raw)) {
                    accumulateBitmask(plc_success_ptr, pos);
                    n_success++;
                    coordCopy<double>(&trial, pos, dcs, 0);
                  }
                  break;
                case PrecisionModel::SINGLE:
                  coordCopy<float>(&fcs, 0, trial, pos);
                  randomizeSystem<float>(&fcs, imap_index, xrs, vara_raw, varb_raw, varc_raw);
                  permuteSystem<float>(&fcs, imap_index, vara_raw, k);
                  if (varb_raw.index >= 0) {
                    permuteSystem<float>(&fcs, imap_index, varb_raw, m);
                  }
                  if (varc_raw.index >= 0) {
                    permuteSystem<float>(&fcs, imap_index, varc_raw, n);
                  }
                  if (resolveClashes<float>(&fcs, imap_index, xrs, imask, clrep, iteration_limit,
                                            max_clashes, vara_raw, varb_raw, varc_raw)) {
                    accumulateBitmask(plc_success_ptr, pos);
                    n_success++;
                    coordCopy<float>(&trial, pos, fcs, 0);
                  }
                  break;
                }
                pos++;
              }
            }
          }
        }
      }
      break;
    case SamplingIntensity::EXHAUSTIVE:
      {
        // The CPU can make use of the TickCounter object created for this map.  The GPU will need
        // an array of integers for the system, combined with limits on those settings, as is
        // available in the system_settings and system_settings_limits arrays.
        TickCounter itrack = state_trackers[imap_index];
        const int ni_var = itrack.getVariableCount();
        const int rotg_llim = permutor_rotatable_group_bounds.readHost(imap_index);
        const int ni_rotg = permutor_rotatable_group_bounds.readHost(imap_index + 1) - rotg_llim;
        const int ctxg_llim = permutor_cis_trans_group_bounds.readHost(imap_index);
        const int ni_ctxg = permutor_cis_trans_group_bounds.readHost(imap_index + 1) - ctxg_llim;
        const int ni_rotg_ctxg = ni_rotg + ni_ctxg;
        const int invg_llim = permutor_invertible_group_bounds.readHost(imap_index);
        std::vector<int> zero_state(ni_var, 0);
        const std::vector<int>& itrack_watch = itrack.getSettings();
        itrack.set(zero_state);
        do {
          if (pos >= design_ptr[i + 1].x) {
            rtErr("The position counter for system " + std::to_string(i) + ", with exhaustive "
                  "sampling of all mutable elements, has exceeded the number of replicas (" +
                  std::to_string(pos) + " >= " + std::to_string(design_ptr[i + 1].x) + ").",
                  "SynthesisPermutor", "buildSynthesis");
          }
          switch (prec) {
          case PrecisionModel::DOUBLE:
            coordCopy<double>(&dcs, 0, trial, pos);
            break;
          case PrecisionModel::SINGLE:
            coordCopy<float>(&fcs, 0, trial, pos);
            break;
          }
          for (int j = 0; j < ni_var; j++) {
            CoupledEdit vara_raw;
            if (j < ni_rotg) {
              vara_raw = { ConformationEdit::BOND_ROTATION, rotg_llim + j };
            }
            else if (j < ni_rotg_ctxg) {
              vara_raw = { ConformationEdit::CIS_TRANS_FLIP, ctxg_llim + j - ni_rotg };
            }
            else {
              vara_raw = { ConformationEdit::CHIRAL_INVERSION, invg_llim + j - ni_rotg_ctxg };
            }
            switch (prec) {
            case PrecisionModel::DOUBLE:
              permuteSystem<double>(&dcs, imap_index, vara_raw, itrack_watch[j]);
              break;
            case PrecisionModel::SINGLE:
              permuteSystem<float>(&fcs, imap_index, vara_raw, itrack_watch[j]);
              break;
            }
          }
          switch (prec) {
          case PrecisionModel::DOUBLE:
            if (clrep == nullptr || detectClash<double, double>(&dcs, 0,
                                                                clrep->getTopologyPointer(),
                                                                &imask, clrep) == false) {
              accumulateBitmask(plc_success_ptr, pos);
              n_success++;
              coordCopy<double>(&trial, pos, dcs, 0);
            }
            break;
          case PrecisionModel::SINGLE:
            if (clrep == nullptr || detectClash<float, float>(&fcs, 0, clrep->getTopologyPointer(),
                                                              &imask, clrep) == false) {
              accumulateBitmask(plc_success_ptr, pos);
              n_success++;
              coordCopy<float>(&trial, pos, fcs, 0);
            }
            break;
          }
          pos++;
          itrack.advance();
        } while (sum<int>(itrack_watch) > 0);
      }
      break;
    }
  }
//#endif
  
  // Report the correspondence, if requested.
  if (correspondence != nullptr) {
    correspondence->resize(n_success);
    int* corr_data = correspondence->data();
    int cpos = 0;
    for (int i = 0; i < system_count; i++) {
      for (int j = design_ptr[i].x; j < design_ptr[i + 1].x; j++) {
        if (readBitFromMask(plc_success_ptr, j)) {
          corr_data[cpos] = i;
          cpos++;
        }
      }
    }
  }
  
  // Determine whether the trial synthesis contains any invalid structures, or if it can be
  // returned immediately as the result.  If necessary, filter out the invalid structures.
  if (n_success == design_ptr[system_count].x) {
    return trial;
  }
  else {
    std::vector<int> revised_coord_index_key(n_success);
    std::vector<int> revised_topol_index_key(n_success);
    int k = 0;
    for (int i = 0; i < system_count; i++) {
      const int orig_topol_idx = poly_ps_ptr->getUniqueTopologyIndex(i);
      for (int j = design_ptr[i].x; j < design_ptr[i + 1].x; j++) {
        if (readBitFromMask(plc_success_ptr, j)) {
          revised_coord_index_key[k] = i;
          revised_topol_index_key[k] = orig_topol_idx;
          k++;
        }
      }
    }
    PhaseSpaceSynthesis result(orig_coords, revised_coord_index_key, orig_tops,
                               revised_topol_index_key, globalpos_scale_bits, localpos_scale_bits,
                               velocity_scale_bits, force_scale_bits);
    int cpos = 0;
    for (int i = 0; i < system_count; i++) {
      for (int j = design_ptr[i].x; j < design_ptr[i + 1].x; j++) {
        if (readBitFromMask(plc_success_ptr, j)) {
          coordCopy(&result, cpos, trial, j);
          cpos++;
        }
      }
    }
    return result;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis
SynthesisPermutor::buildSynthesis(const SamplingIntensity effort, Xoshiro256ppGenerator *xrs,
                                  std::vector<int> *correspondence, const PrecisionModel prec,
                                  const int system_limit, ClashReport *clrep,
                                  const int iteration_limit, const int max_clashes) const {
  return buildSynthesis(effort, xrs, correspondence, system_limit, default_globalpos_scale_bits,
                        default_localpos_scale_bits, default_velocity_scale_bits,
                        default_force_scale_bits, prec, clrep, iteration_limit, max_clashes);
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis
SynthesisPermutor::buildSynthesis(const SamplingIntensity effort, Xoshiro256ppGenerator *xrs,
                                  const int system_limit, const int globalpos_scale_bits,
                                  const int localpos_scale_bits, const int velocity_scale_bits,
                                  const int force_scale_bits, const PrecisionModel prec,
                                  ClashReport *clrep, const int iteration_limit,
                                  const int max_clashes) const {
  return buildSynthesis(effort, xrs, nullptr, system_limit, globalpos_scale_bits,
                        localpos_scale_bits, velocity_scale_bits, force_scale_bits, prec, clrep,
                        iteration_limit, max_clashes);
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis
SynthesisPermutor::buildSynthesis(const SamplingIntensity effort, Xoshiro256ppGenerator *xrs,
                                  const PrecisionModel prec, const int system_limit, 
                                  ClashReport *clrep, const int iteration_limit,
                                  const int max_clashes) const {
  return buildSynthesis(effort, xrs, nullptr, system_limit, default_globalpos_scale_bits,
                        default_localpos_scale_bits, default_velocity_scale_bits,
                        default_force_scale_bits, prec, clrep, iteration_limit, max_clashes);
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::rebasePointers() {
  permutor_map_indices.swapTarget(&synthesis_data);
  system_settings.swapTarget(&synthesis_data);
  system_settings_limits.swapTarget(&synthesis_data);
  rotatable_group_atoms.swapTarget(&group_data);
  rotatable_group_bounds.swapTarget(&group_data);
  permutor_rotatable_group_bounds.swapTarget(&group_data);  
  cis_trans_group_atoms.swapTarget(&group_data);
  cis_trans_group_bounds.swapTarget(&group_data);
  permutor_cis_trans_group_bounds.swapTarget(&group_data);  
  invertible_group_atoms.swapTarget(&group_data);
  invertible_group_bounds.swapTarget(&group_data);
  permutor_invertible_group_bounds.swapTarget(&group_data);
  invertible_atom_centers.swapTarget(&group_data);
  invertible_group_protocols.swapTarget(&group_data);
  permutor_element_bounds.swapTarget(&group_data);
  light_sampling_foci_bounds.swapTarget(&group_data);
  heavy_sampling_foci_bounds.swapTarget(&group_data);
  permutor_elements.swapTarget(&variables_data);
  light_sampling_foci.swapTarget(&variables_data);
  heavy_sampling_foci.swapTarget(&variables_data);
  rotatable_bond_markers.swapTarget(&marker_data);
  cis_trans_bond_markers.swapTarget(&marker_data);
  chiral_markers.swapTarget(&marker_data);
  rotatable_bond_settings_bounds.swapTarget(&group_data);
  cis_trans_bond_settings_bounds.swapTarget(&group_data);
  chiral_settings_bounds.swapTarget(&group_data);
}

//------------------------------------------------------------------------------------------------
std::vector<ChemicalFeatures> SynthesisPermutor::temporaryFeatures(StopWatch *timer) const {
  std::vector<ChemicalFeatures> result(permutor_map_count);
  result.reserve(topologies.size());
  for (int i = 0; i < permutor_map_count; i++) {
    result.emplace_back(ChemicalFeatures(topologies[i], MapRotatableGroups::YES, 300.0, timer));
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const std::vector<ChemicalFeatures*>
SynthesisPermutor::temporaryFeaturesPointers(const std::vector<ChemicalFeatures> &feat_in) const {
  std::vector<ChemicalFeatures*> result(permutor_map_count);
  for (int i = 0; i < permutor_map_count; i++) {

    // Even passing the vector by reference, pointers to the individual elements will be valid.
    // A pointer to the vector itself would not.
    result[i] = const_cast<ChemicalFeatures*>(&feat_in[i]);
  }  
  return result;
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::fillPermutorMaps(const std::vector<ChemicalFeatures*> &chemfe_in,
                                         StopWatch *timer) {
  
  // Check that, for each set of features, rotatable bond groups have been mapped.  If not, create
  // a new ChemicalFeatures object which does map the rotating groups (the input features are taken
  // as const).  
  features.resize(permutor_map_count);
  std::vector<ChemicalFeatures> fallback_features;
  int n_unmapped = 0;
  for (int i = 0; i < permutor_map_count; i++) {
    n_unmapped += (chemfe_in[i]->rotatableGroupsMapped() == false);
  }
  fallback_features.reserve(n_unmapped);
  n_unmapped = 0;
  for (int i = 0; i < permutor_map_count; i++) {
    if (chemfe_in[i]->rotatableGroupsMapped() == false) {
      fallback_features.emplace_back(chemfe_in[i]->getTopologyPointer(), MapRotatableGroups::YES,
                                     300.0, timer);
      features[i] = &fallback_features[n_unmapped];
      n_unmapped++;
    }
    else {
      features[i] = chemfe_in[i];
    }
  }

  // If any features had to be recalculated in order to get the maps, the temporary objects will
  // not persist and the feature pointers as a whole should be considered unreliable.
  if (n_unmapped > 0) {
    features_objects_accessible = false;
  }

  // Obtain the isomerization groups
  std::vector<int> tmp_permutor_rotatable_group_bounds(permutor_map_count + 1);
  std::vector<int> tmp_permutor_cis_trans_group_bounds(permutor_map_count + 1);
  std::vector<int> tmp_permutor_invertible_group_bounds(permutor_map_count + 1);
  int n_rotg = 0;
  int n_ctxg = 0;
  int n_invg = 0;
  for (int i = 0; i < permutor_map_count; i++) {
    tmp_permutor_rotatable_group_bounds[i] = n_rotg;
    tmp_permutor_cis_trans_group_bounds[i] = n_ctxg;
    tmp_permutor_invertible_group_bounds[i] = n_invg;
    n_rotg += features[i]->getRotatableBondCount();
    n_ctxg += features[i]->getCisTransBondCount();
    n_invg += features[i]->getChiralCenterCount();
  }
  tmp_permutor_rotatable_group_bounds[permutor_map_count] = n_rotg;
  tmp_permutor_cis_trans_group_bounds[permutor_map_count] = n_ctxg;
  tmp_permutor_invertible_group_bounds[permutor_map_count] = n_invg;
  rotatable_groups.reserve(n_rotg);
  cis_trans_groups.reserve(n_ctxg);
  invertible_groups.reserve(n_invg);

  // Lay out the arrays of rotatable and cis-trans isomeric bond settings bounds.  The number of
  // variable elements in each system is known at this time--the exact settings and number of
  // settings for each of them are not, but will be filled in the setVariableRanges() function.
  const std::vector<int> tmp_rotatable_bond_settings_bounds(n_rotg + 1, 0);
  const std::vector<int> tmp_cis_trans_bond_settings_bounds(n_ctxg + 1, 0);
  std::vector<int> tmp_chiral_settings_bounds(n_invg + 1, 0);
  std::vector<int4> tmp_rotatable_bond_markers(n_rotg);
  std::vector<int4> tmp_cis_trans_bond_markers(n_ctxg);
  std::vector<int4> tmp_chiral_markers(n_invg);
  std::vector<int> tmp_invertible_atom_centers(n_invg);
  std::vector<int> tmp_invertible_group_protocols(n_invg);
  n_rotg = 0;
  n_ctxg = 0;
  n_invg = 0;
  int chiral_states = 0;
  for (int i = 0; i < permutor_map_count; i++) {
    const int ni_rotg = features[i]->getRotatableBondCount();
    const int ni_ctxg = features[i]->getCisTransBondCount();
    const int ni_invg = features[i]->getChiralCenterCount();
    const std::vector<IsomerPlan> rotg_vec = features[i]->getRotatableBondGroups();
    for (int j = 0; j < ni_rotg; j++) {
      rotatable_groups.push_back(rotg_vec[j]);
      tmp_rotatable_bond_markers[n_rotg] = { rotg_vec[j].getRootHandle(),
                                             rotg_vec[j].getRootAtom(), rotg_vec[j].getPivotAtom(),
                                             rotg_vec[j].getPivotHandle() };
      n_rotg++;
    }
    const std::vector<IsomerPlan> ctxg_vec = features[i]->getCisTransIsomerizationGroups();
    for (int j = 0; j < ni_ctxg; j++) {
      cis_trans_groups.push_back(ctxg_vec[j]);
      tmp_cis_trans_bond_markers[n_ctxg] = { ctxg_vec[j].getRootHandle(),
                                             ctxg_vec[j].getRootAtom(), ctxg_vec[j].getPivotAtom(),
                                             ctxg_vec[j].getPivotHandle() };
      n_ctxg++;
    }
    const std::vector<IsomerPlan> invg_vec = features[i]->getChiralInversionGroups();
    const std::vector<int4> chiral_base_atoms = features[i]->getChiralArmBaseAtoms();
    const std::vector<int> chiral_center_atoms = features[i]->getChiralCenters();
    for (int j = 0; j < ni_invg; j++) {
      tmp_invertible_atom_centers[n_invg] = chiral_center_atoms[j];
      tmp_invertible_group_protocols[n_invg] = static_cast<int>(invg_vec[j].getChiralPlan());
      invertible_groups.push_back(invg_vec[j]);
      tmp_chiral_markers[n_invg] = chiral_base_atoms[j];
      switch (invg_vec[j].getChiralPlan()) {
      case ChiralInversionProtocol::ROTATE:
      case ChiralInversionProtocol::REFLECT:
        chiral_states += 2;
        break;
      case ChiralInversionProtocol::DO_NOT_INVERT:
        chiral_states += 1;
        break;
      }
      n_invg++;
      tmp_chiral_settings_bounds[n_invg] = chiral_states;
    }
  }
  total_rotatable_groups = n_rotg;
  total_cis_trans_groups = n_ctxg;
  total_invertible_groups = n_invg;

  // Count the numbers of atoms in each distinct system's rotatable bond groups.  Track to ensure
  // that the group mappings do not exceed a maximum of two billion (2^31) atoms, and raise an
  // exception if they do.  The int-format mapping vectors are limited, but any problem set that
  // hits these limits is almost certain to run up against other barriers like the memory
  // available to the card.
  bool problem = false;
  std::vector<int> tmp_rotatable_group_bounds(n_rotg + 1);
  std::vector<int> tmp_cis_trans_group_bounds(n_ctxg + 1);
  std::vector<int> tmp_invertible_group_bounds(n_invg + 1);
  int rotg_atom_count = 0;
  for (int i = 0; i < tmp_permutor_rotatable_group_bounds[permutor_map_count]; i++) {
    tmp_rotatable_group_bounds[i] = rotg_atom_count;
    rotg_atom_count += rotatable_groups[i].getMovingAtomCount() + 2;
    problem = (problem || rotg_atom_count < 0);
  }
  int ctxg_atom_count = 0;
  for (int i = 0; i < tmp_permutor_cis_trans_group_bounds[permutor_map_count]; i++) {
    tmp_cis_trans_group_bounds[i] = ctxg_atom_count;
    ctxg_atom_count += cis_trans_groups[i].getMovingAtomCount() + 2;
    problem = (problem || ctxg_atom_count < 0);
  }
  int invg_atom_count = 0;
  for (int i = 0; i < tmp_permutor_invertible_group_bounds[permutor_map_count]; i++) {
    tmp_invertible_group_bounds[i] = invg_atom_count;
    switch (invertible_groups[i].getChiralPlan()) {
    case ChiralInversionProtocol::ROTATE:
    case ChiralInversionProtocol::REFLECT:
      invg_atom_count += invertible_groups[i].getMovingAtomCount() + 2;
      break;
    case ChiralInversionProtocol::DO_NOT_INVERT:
      invg_atom_count += 2;
      break;
    }
    problem = (problem || invg_atom_count < 0);
  }
  tmp_rotatable_group_bounds[n_rotg] = rotg_atom_count;
  tmp_cis_trans_group_bounds[n_ctxg] = ctxg_atom_count;
  tmp_invertible_group_bounds[n_invg] = invg_atom_count;
  if (problem) {
    rtErr("Rotatable, cis-trans, or invertible atom groups specify too many atoms for available "
          "indexing.  Reduce the number of systems if possible.", "SynthesisPermutor",
          "FillPermutorMaps");
  }
  
  // Allocate temporary storage arrays for the moving atom groups, then fill them.  Compile arrays
  // of each map's elements.  Deduce the couplings between variable elements in each structure.
  // This can be done irrespective of the exact number of settings that each element will have,
  // which is yet to be determined.  The foci arrays will be concatenated with the
  // tmp_permutor_elements array into the same ARRAY-kind Hybrid object from grouped upload and
  // download.
  std::vector<int> tmp_rotatable_group_atoms(rotg_atom_count);
  std::vector<int> tmp_cis_trans_group_atoms(ctxg_atom_count);
  std::vector<int> tmp_invertible_group_atoms(invg_atom_count);
  rotg_atom_count = 0;
  ctxg_atom_count = 0;
  invg_atom_count = 0;
  std::vector<int2> tmp_permutor_elements(n_rotg + n_ctxg + n_invg);
  std::vector<int2> tmp_light_sampling_foci, tmp_heavy_sampling_foci;
  std::vector<int> tmp_permutor_element_bounds(permutor_map_count + 1);
  std::vector<int> tmp_light_sampling_foci_bounds(permutor_map_count + 1);
  std::vector<int> tmp_heavy_sampling_foci_bounds(permutor_map_count + 1);
  int elem_idx = 0;
  int rotg_idx = 0;
  int ctxg_idx = 0;
  int invg_idx = 0;
  for (int i = 0; i < permutor_map_count; i++) {

    // Set the bounds entries
    tmp_permutor_element_bounds[i] = elem_idx;
    tmp_light_sampling_foci_bounds[i] = tmp_light_sampling_foci.size();
    tmp_heavy_sampling_foci_bounds[i] = tmp_heavy_sampling_foci.size();
    
    // Rotatable bond groups
    const std::vector<IsomerPlan> rotg_vec = features[i]->getRotatableBondGroups();
    const int ni_rotg = rotg_vec.size();
    for (int j = 0; j < ni_rotg; j++) {
      tmp_rotatable_group_atoms[rotg_atom_count] = rotg_vec[j].getRootAtom();
      rotg_atom_count++;
      tmp_rotatable_group_atoms[rotg_atom_count] = rotg_vec[j].getPivotAtom();
      rotg_atom_count++;
      const std::vector<int>& mv_atoms = rotg_vec[j].getMovingAtoms();
      const int nij_atoms = rotg_vec[j].getMovingAtomCount();
      for (int k = 0; k < nij_atoms; k++) {
        tmp_rotatable_group_atoms[rotg_atom_count] = mv_atoms[k];
        rotg_atom_count++;
      }

      // Handle the permutor elements
      tmp_permutor_elements[elem_idx] = { static_cast<int>(ConformationEdit::BOND_ROTATION),
                                          rotg_idx };
      elem_idx++;
      rotg_idx++;
    }

    // Cis-trans isomeric bond groups
    const std::vector<IsomerPlan> ctxg_vec = features[i]->getCisTransIsomerizationGroups();
    const int ni_ctxg = ctxg_vec.size();
    for (int j = 0; j < ni_ctxg; j++) {
      tmp_cis_trans_group_atoms[ctxg_atom_count] = ctxg_vec[j].getRootAtom();
      ctxg_atom_count++;
      tmp_cis_trans_group_atoms[ctxg_atom_count] = ctxg_vec[j].getPivotAtom();
      ctxg_atom_count++;
      const std::vector<int>& mv_atoms = ctxg_vec[j].getMovingAtoms();
      const int nij_atoms = ctxg_vec[j].getMovingAtomCount();
      for (int k = 0; k < nij_atoms; k++) {
        tmp_cis_trans_group_atoms[ctxg_atom_count] = mv_atoms[k];
        ctxg_atom_count++;
      }

      // Handle the permutor elements
      tmp_permutor_elements[elem_idx] = { static_cast<int>(ConformationEdit::CIS_TRANS_FLIP),
                                          ctxg_idx };
      elem_idx++;
      ctxg_idx++;
    }

    // Invertible chiral centers
    const std::vector<IsomerPlan> invg_vec = features[i]->getChiralInversionGroups();    
    const int ni_invg = invg_vec.size();
    for (int j = 0; j < ni_invg; j++) {
      tmp_invertible_group_atoms[invg_atom_count] = invg_vec[j].getRootAtom();
      invg_atom_count++;
      tmp_invertible_group_atoms[invg_atom_count] = static_cast<int>(invg_vec[j].getPivotAtom());
      invg_atom_count++;
      switch (invg_vec[j].getChiralPlan()) {
      case ChiralInversionProtocol::ROTATE:
      case ChiralInversionProtocol::REFLECT:
        {
          const std::vector<int>& mv_atoms = invg_vec[j].getMovingAtoms();
          const int nij_atoms = invg_vec[j].getMovingAtomCount();
          for (int k = 0; k < nij_atoms; k++) {
            tmp_invertible_group_atoms[invg_atom_count] = mv_atoms[k];
            invg_atom_count++;
          }
        }
        break;
      case ChiralInversionProtocol::DO_NOT_INVERT:
        break;
      }

      // Handle the permutor elements
      tmp_permutor_elements[elem_idx] = { static_cast<int>(ConformationEdit::CHIRAL_INVERSION),
                                          invg_idx };
      elem_idx++;
      invg_idx++;
    }

    // Compose a single array of all mutable elements (this process will be repeated in the
    // setVariableRanges() function and can be expensive for all of its memory transfer, but there
    // is not a clean way to have the settings in pplace prior to this step).
    const int ni_var = ni_rotg + ni_ctxg + ni_invg;
    std::vector<IsomerPlan> all_mutables;
    all_mutables.reserve(ni_var);
    all_mutables.insert(all_mutables.end(), rotg_vec.begin(), rotg_vec.end());
    all_mutables.insert(all_mutables.end(), ctxg_vec.begin(), ctxg_vec.end());
    all_mutables.insert(all_mutables.end(), invg_vec.begin(), invg_vec.end());
    const NonbondedKit<double> nbk = topologies[i]->getDoublePrecisionNonbondedKit();
    std::vector<bool> var_linked(ni_var * ni_var, true);
    for (int j = 1; j < ni_var; j++) {
      for (int k = 0; k < j; k++) {
        const bool jk_link = permutationsAreLinked(all_mutables, j, k, nbk);
        var_linked[(j * ni_var) + k] = jk_link;
        var_linked[(k * ni_var) + j] = jk_link;
      }
    }

    // When mapping out the coupled variables, some couplings may be identified that actually have
    // only one possible permutation as the variables each have only one state.  Those cases will
    // be skipped in the processing stage.  Without knowing the enumerated possible settings for
    // each variable a-priori, all connected variables must be taken into account.
    std::vector<int2> var_codes(ni_var);
    for (int j = 0; j < ni_rotg; j++) {
      var_codes[j] = { static_cast<int>(ConformationEdit::BOND_ROTATION),
                       rotg_idx - ni_rotg + j };
    }
    for (int j = ni_rotg; j < ni_rotg + ni_ctxg; j++) {
      var_codes[j] = { static_cast<int>(ConformationEdit::CIS_TRANS_FLIP),
                       ctxg_idx - ni_ctxg + j - ni_rotg };
    }
    for (int j = ni_rotg + ni_ctxg; j < ni_var; j++) {
      var_codes[j] = { static_cast<int>(ConformationEdit::CHIRAL_INVERSION),
                       invg_idx - ni_invg + j - (ni_rotg + ni_ctxg) };
    }
    const int2 blank_code = { static_cast<int>(ConformationEdit::BOND_ROTATION), -1 };
    for (int j = 0; j < ni_var; j++) {
      bool one_linked = false;
      for (int k = 0; k < j; k++) {
        if (var_linked[(j * ni_var) + k]) {
          one_linked = true;
          tmp_light_sampling_foci.push_back(var_codes[j]);
          tmp_light_sampling_foci.push_back(var_codes[k]);
          bool two_linked = false;
          for (int m = 0; m < k; m++) {
            if (var_linked[(k * ni_var) + m]) {
              two_linked = true;
              tmp_heavy_sampling_foci.push_back(var_codes[j]);
              tmp_heavy_sampling_foci.push_back(var_codes[k]);
              tmp_heavy_sampling_foci.push_back(var_codes[m]);
            }
          }
          if (two_linked == false) {

            // Add all combinations of { j, k, -1 } to the list of heavy sampling tasks.
            tmp_heavy_sampling_foci.push_back(var_codes[j]);
            tmp_heavy_sampling_foci.push_back(var_codes[k]);
            tmp_heavy_sampling_foci.push_back(blank_code);
          }
        }
      }
      if (one_linked == false) { 

        // Add { j, -1 } to the list of light sampling tasks and { j, -1, -1 } to the list of
        // heavy sampling tasks.
        tmp_light_sampling_foci.push_back(var_codes[j]);
        tmp_light_sampling_foci.push_back(blank_code);
        tmp_heavy_sampling_foci.push_back(var_codes[j]);
        tmp_heavy_sampling_foci.push_back(blank_code);
        tmp_heavy_sampling_foci.push_back(blank_code);
      }
    }
  }

  // Add the final bounds to the permutor elements as well as foci arrays
  tmp_permutor_element_bounds[permutor_map_count] = elem_idx;
  tmp_light_sampling_foci_bounds[permutor_map_count] = tmp_light_sampling_foci.size();
  tmp_heavy_sampling_foci_bounds[permutor_map_count] = tmp_heavy_sampling_foci.size();

  // Allocate Hybrid data and load it with the staging arrays
  const size_t gd_length = roundUp(tmp_rotatable_group_atoms.size(), warp_size_zu) +
                           roundUp(tmp_rotatable_group_bounds.size(), warp_size_zu) +
                           roundUp(tmp_permutor_rotatable_group_bounds.size(), warp_size_zu) +
                           roundUp(tmp_cis_trans_group_atoms.size(), warp_size_zu) +
                           roundUp(tmp_cis_trans_group_bounds.size(), warp_size_zu) +
                           roundUp(tmp_permutor_cis_trans_group_bounds.size(), warp_size_zu) +
                           roundUp(tmp_invertible_group_atoms.size(), warp_size_zu) +
                           roundUp(tmp_invertible_group_bounds.size(), warp_size_zu) +
                           roundUp(tmp_permutor_invertible_group_bounds.size(), warp_size_zu) +
                           roundUp(tmp_invertible_atom_centers.size(), warp_size_zu) +
                           roundUp(tmp_invertible_group_protocols.size(), warp_size_zu) +
                           roundUp(tmp_rotatable_bond_settings_bounds.size(), warp_size_zu) +
                           roundUp(tmp_cis_trans_bond_settings_bounds.size(), warp_size_zu) +
                           roundUp(tmp_chiral_settings_bounds.size(), warp_size_zu) +
                           roundUp(tmp_permutor_element_bounds.size(), warp_size_zu) +
                           roundUp(tmp_light_sampling_foci_bounds.size(), warp_size_zu) +
                           roundUp(tmp_light_sampling_foci_bounds.size(), warp_size_zu);
  group_data.resize(gd_length);
  const size_t mk_length = roundUp(tmp_rotatable_bond_markers.size(), warp_size_zu) +
                           roundUp(tmp_cis_trans_bond_markers.size(), warp_size_zu) +
                           roundUp(tmp_chiral_markers.size(), warp_size_zu);
  marker_data.resize(mk_length);
  const size_t pe_length = roundUp(tmp_permutor_elements.size(), warp_size_zu) +
                           roundUp(tmp_light_sampling_foci.size(), warp_size_zu) +
                           roundUp(tmp_heavy_sampling_foci.size(), warp_size_zu);
  variables_data.resize(pe_length);
  size_t ic = 0;
  ic = rotatable_group_atoms.putHost(&group_data, tmp_rotatable_group_atoms, ic, warp_size_zu);
  ic = rotatable_group_bounds.putHost(&group_data, tmp_rotatable_group_bounds, ic, warp_size_zu);
  ic = permutor_rotatable_group_bounds.putHost(&group_data, tmp_permutor_rotatable_group_bounds,
                                               ic, warp_size_zu);
  ic = cis_trans_group_atoms.putHost(&group_data, tmp_cis_trans_group_atoms, ic, warp_size_zu);
  ic = cis_trans_group_bounds.putHost(&group_data, tmp_cis_trans_group_bounds, ic, warp_size_zu);
  ic = permutor_cis_trans_group_bounds.putHost(&group_data, tmp_permutor_cis_trans_group_bounds,
                                               ic, warp_size_zu);
  ic = invertible_group_atoms.putHost(&group_data, tmp_invertible_group_atoms, ic, warp_size_zu);
  ic = invertible_group_bounds.putHost(&group_data, tmp_invertible_group_bounds, ic, warp_size_zu);
  ic = permutor_invertible_group_bounds.putHost(&group_data, tmp_permutor_invertible_group_bounds,
                                                ic, warp_size_zu);
  ic = invertible_atom_centers.putHost(&group_data, tmp_invertible_atom_centers, ic, warp_size_zu);
  ic = invertible_group_protocols.putHost(&group_data, tmp_invertible_group_protocols, ic,
                                          warp_size_zu);
  ic = rotatable_bond_settings_bounds.putHost(&group_data, tmp_rotatable_bond_settings_bounds,
                                              ic, warp_size_zu);
  ic = cis_trans_bond_settings_bounds.putHost(&group_data, tmp_cis_trans_bond_settings_bounds,
                                              ic, warp_size_zu);
  ic = chiral_settings_bounds.putHost(&group_data, tmp_chiral_settings_bounds, ic, warp_size_zu);
  ic = permutor_element_bounds.putHost(&group_data, tmp_permutor_element_bounds, ic, warp_size_zu);
  ic = light_sampling_foci_bounds.putHost(&group_data, tmp_light_sampling_foci_bounds, ic,
                                          warp_size_zu);
  ic = heavy_sampling_foci_bounds.putHost(&group_data, tmp_heavy_sampling_foci_bounds, ic,
                                          warp_size_zu);
  ic = 0;
  ic = rotatable_bond_markers.putHost(&marker_data, tmp_rotatable_bond_markers, ic, warp_size_zu);
  ic = cis_trans_bond_markers.putHost(&marker_data, tmp_cis_trans_bond_markers, ic, warp_size_zu);
  ic = chiral_markers.putHost(&marker_data, tmp_chiral_markers, ic, warp_size_zu);
  ic = 0;
  ic = permutor_elements.putHost(&variables_data, tmp_permutor_elements, ic, warp_size_zu);
  ic = light_sampling_foci.putHost(&variables_data, tmp_light_sampling_foci, ic, warp_size_zu);
  ic = heavy_sampling_foci.putHost(&variables_data, tmp_heavy_sampling_foci, ic, warp_size_zu);
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::validateSystemIndex(const int system_index) const {
  if (system_index < 0 || system_index >= system_count) {
    rtErr("System index " + std::to_string(system_index) + " is out of bounds for a collection "
          "of " + std::to_string(system_count) + " systems.", "SynthesisPermutor",
          "validateSystemIndex");
  }
}

//-------------------------------------------------------------------------------------------------
void SynthesisPermutor::validateMapIndex(const int map_index) const {
  if (map_index < 0 || map_index >= permutor_map_count) {
    rtErr("Permutor map index " + std::to_string(map_index) + " is out of bounds for a collection "
          "of " + std::to_string(permutor_map_count) + " maps.", "SynthesisPermutor",
          "validateMapIndex");
  }
}

} // namespace synthesis
} // namespace stormm
