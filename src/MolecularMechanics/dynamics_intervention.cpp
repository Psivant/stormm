#include "copyright.h"
#include "Reporting/error_format.h"
#include "Structure/structure_enumerators.h"
#include "Topology/atomgraph_enumerators.h"
#include "dynamics_intervention.h"

namespace stormm {
namespace mm {

using card::HybridKind;
using errors::rtErr;
using structure::ApplyConstraints;
using topology::UnitCellType;
  
//-------------------------------------------------------------------------------------------------
DynamicsIntervention::DynamicsIntervention() :
    active_interventions{false}, investigate_large_forces{false}, check_neighbor_list{false},
    relevant_interval_count{0}, relevant_intervals{}, gpu{null_gpu},

    // Arrays of abstracts, along with descriptors that help accessors to choose the correct
    // abstract for computing on the CPU or GPU.  Additional descriptors may be relevant to
    // choosing the proper abstract given the state of the underlying coordiante object. 
    topology_abstract_tiers{}, dpoly_vks{}, fpoly_vks{}, dpoly_nbks{}, fpoly_nbks{}, dpoly_rks{},
    fpoly_rks{}, dpoly_auks{}, fpoly_auks{},
    coordinate_abstract_tiers{}, coordinate_cycle_stages{}, poly_psws{}, poly_psrs{},
    condensate_abstract_tiers{}, cdnsws{},

    // Details of the simulation neighbor list
    ngbr_list_tcoord{0}, ngbr_list_tacc{0}, ngbr_list_tcalc{0}, ngbr_list_tcoord4{0},
    nl_layout{NeighborListKind::MONO}, neighbor_list_abstract_tiers{},
    neighbor_list_cycle_stages{}, neighbor_list_themes{}, voided_cgws{},

    // More arrays of abstracts
    implicit_solvent_abstract_tiers{}, implicit_solvent_cycle_stages{}, d_iswks{}, f_iswks{},
    exclusion_mask_abstract_tiers{}, lemrs{}, semrs{}, utracking_abstract_tiers{}, scws{},
    anomaly_abstract_tiers{}, bugws{}, hbond_abstract_tiers{}, hbws{},

    // Pointers to the underlying resources attached to this intervention object
    poly_ps_ptr{nullptr}, cdns_ptr{nullptr}, cg_ptr{},
    poly_ag_ptr{nullptr}, isw_ptr{nullptr}, lem_ptr{nullptr}, sc_ptr{nullptr}, anom_ptr{nullptr},

    // Arrays of pointers to sets of analyses attached to this object
    hbond_trk_ptr{},

    // Action queues
    force_calc_actions{}, velocity_adv_actions{}, velocity_cnst_actions{},
    kinetic_calc_actions{}, position_adv_actions{}, geometry_cnst_actions{},

    // Intervals at which to carry out each action in the near-eponymous arrays above
    force_calc_action_intv{}, velocity_adv_action_intv{}, velocity_cnst_action_intv{},
    kinetic_calc_action_intv{}, position_adv_action_intv{}, geometry_cnst_action_intv{},

    // Indices of the analysis or other attached resource used by near-eponymous actions above
    force_calc_action_index{}, velocity_adv_action_index{}, velocity_cnst_action_index{},
    kinetic_calc_action_index{}, position_adv_action_index{}, geometry_cnst_action_index{},

    // Compute resource upon which to carry out each action in the near-eponymous arrays above
    force_calc_action_tier{}, velocity_adv_action_tier{}, velocity_cnst_action_tier{},
    kinetic_calc_action_tier{}, position_adv_action_tier{}, geometry_cnst_action_tier{},

    // Workspaces for extra coordinate and force sets will be allocated as needed
    workspaces{}, nl_workspaces{}
{}

//-------------------------------------------------------------------------------------------------
bool DynamicsIntervention::isActive() const {
  return active_interventions;
}

//-------------------------------------------------------------------------------------------------
bool DynamicsIntervention::isActiveOnStep(const int step_number) const {
  bool matched = false;
  for (int i = 0; i < relevant_interval_count; i++) {
    matched = (matched || (step_number % relevant_intervals[i]) == 0);
  }
  return matched;
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& DynamicsIntervention::getActiveIntervals() const {
  return relevant_intervals;
}

//-------------------------------------------------------------------------------------------------
const GpuDetails& DynamicsIntervention::getGpuDetails() const {
  return gpu;
}

//-------------------------------------------------------------------------------------------------
int DynamicsIntervention::getSMPCount() const {
  return gpu.getSMPCount();
}

//-------------------------------------------------------------------------------------------------
bool DynamicsIntervention::hasPhaseSpaceSynthesis() const {
  return (poly_ps_ptr != nullptr);
}

//-------------------------------------------------------------------------------------------------
bool DynamicsIntervention::hasCondensate() const {
  return (cdns_ptr != nullptr);
}

//-------------------------------------------------------------------------------------------------
bool DynamicsIntervention::hasTopologySynthesis() const {
  return (poly_ag_ptr != nullptr);
}

//-------------------------------------------------------------------------------------------------
bool DynamicsIntervention::hasImplicitSolventWorkspace() const {
  return (isw_ptr != nullptr);
}

//-------------------------------------------------------------------------------------------------
bool DynamicsIntervention::hasExclusionMasks() const {
  if (poly_ag_ptr != nullptr) {
    switch (poly_ag_ptr->getUnitCellType()) {
    case UnitCellType::NONE:
      return (se_ptr != nullptr);
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      return (lem_ptr != nullptr);
    }
  }
  else {
    return (lem_ptr != nullptr || se_ptr != nullptr);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool DynamicsIntervention::hasEnergyTracking() const {
  return (sc_ptr != nullptr);
}

//-------------------------------------------------------------------------------------------------
bool DynamicsIntervention::hasAnomalyReporting() const {
  return (anom_ptr != nullptr);
}

//-------------------------------------------------------------------------------------------------
bool DynamicsIntervention::doNeighborListChecks() const {
  return check_neighbor_list;
}

//-------------------------------------------------------------------------------------------------
PsSynthesisWriter DynamicsIntervention::getCoordinateData(const HybridTargetLevel tier) {
  const size_t n_options = coordinate_abstract_tiers.size();
  const CoordinateCycle curr_cc = poly_ps_ptr->getCyclePosition();
  for (size_t i = 0; i < n_options; i++) {
    if (coordinate_abstract_tiers[i] == tier && coordinate_cycle_stages[i] == curr_cc) {
      return poly_psws[i];
    }
  }
  rtErr("No PhaseSpaceSynthesis abstract was available for the main coordinates on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getCoordinateData");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const PsSynthesisReader
DynamicsIntervention::getCoordinateData(const HybridTargetLevel tier) const {
  return this->getReadOnlyCoordinateData(tier);
}

//-------------------------------------------------------------------------------------------------
const PsSynthesisReader
DynamicsIntervention::getReadOnlyCoordinateData(const HybridTargetLevel tier) const {
  const size_t n_options = coordinate_abstract_tiers.size();
  const CoordinateCycle curr_cc = poly_ps_ptr->getCyclePosition();
  for (size_t i = 0; i < n_options; i++) {
    if (coordinate_abstract_tiers[i] == tier && coordinate_cycle_stages[i] == curr_cc) {
      return poly_psrs[i];
    }
  }
  rtErr("No PhaseSpaceSynthesis abstract was available for the main coordinates on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getReadOnlyCoordinateData");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
CondensateWriter DynamicsIntervention::getCondensateData(const HybridTargetLevel tier) {
  const size_t n_options = condensate_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (condensate_abstract_tiers[i] == tier) {
      return cdnsws[i];
    }
  }
  rtErr("No Condensate abstract was available for the auxiliary coordinates on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getCondensateData");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const CondensateReader
DynamicsIntervention::getCondensateData(const HybridTargetLevel tier) const {
  const size_t n_options = condensate_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (condensate_abstract_tiers[i] == tier) {
      return CondensateReader(cdnsws[i]);
    }
  }
  rtErr("No Condensate abstract was available for the auxiliary coordinates on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getCondensateData");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const SyValenceKit<double>
DynamicsIntervention::getDoublePrecisionValenceKit(const HybridTargetLevel tier) const {
  const size_t n_options = topology_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (topology_abstract_tiers[i] == tier) {
      return dpoly_vks[i];
    }
  }
  rtErr("No valence abstract was available for the topology synthesis on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getDoublePrecisionValenceKit");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const SyValenceKit<float>
DynamicsIntervention::getSinglePrecisionValenceKit(const HybridTargetLevel tier) const {
  const size_t n_options = topology_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (topology_abstract_tiers[i] == tier) {
      return fpoly_vks[i];
    }
  }
  rtErr("No valence abstract was available for the topology synthesis on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getSinglePrecisionValenceKit");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const SyNonbondedKit<double, double2>
DynamicsIntervention::getDoublePrecisionNonbondedKit(const HybridTargetLevel tier) const {
  const size_t n_options = topology_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (topology_abstract_tiers[i] == tier) {
      return dpoly_nbks[i];
    }
  }
  rtErr("No valence abstract was available for the topology synthesis on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getDoublePrecisionNonbondedKit");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const SyNonbondedKit<float, float2>
DynamicsIntervention::getSinglePrecisionNonbondedKit(const HybridTargetLevel tier) const {
  const size_t n_options = topology_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (topology_abstract_tiers[i] == tier) {
      return fpoly_nbks[i];
    }
  }
  rtErr("No valence abstract was available for the topology synthesis on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getSinglePrecisionNonbondedKit");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const SyRestraintKit<double, double2, double4_16a>
DynamicsIntervention::getDoublePrecisionRestraintKit(const HybridTargetLevel tier) const {
  const size_t n_options = topology_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (topology_abstract_tiers[i] == tier) {
      return dpoly_rks[i];
    }
  }
  rtErr("No valence abstract was available for the topology synthesis on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getDoublePrecisionRestraintKit");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const SyRestraintKit<float, float2, float4>
DynamicsIntervention::getSinglePrecisionRestraintKit(const HybridTargetLevel tier) const {
  const size_t n_options = topology_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (topology_abstract_tiers[i] == tier) {
      return fpoly_rks[i];
    }
  }
  rtErr("No valence abstract was available for the topology synthesis on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getSinglePrecisionRestraintKit");
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
const SyAtomUpdateKit<double, double2, double4_16a>
DynamicsIntervention::getDoublePrecisionAtomUpdateKit(const HybridTargetLevel tier) const {
  const size_t n_options = topology_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (topology_abstract_tiers[i] == tier) {
      return dpoly_auks[i];
    }
  }
  rtErr("No valence abstract was available for the topology synthesis on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getDoublePrecisionAtomUpdateKit");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const SyAtomUpdateKit<float, float2, float4>
DynamicsIntervention::getSinglePrecisionAtomUpdateKit(const HybridTargetLevel tier) const {
  const size_t n_options = topology_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (topology_abstract_tiers[i] == tier) {
      return fpoly_auks[i];
    }
  }
  rtErr("No valence abstract was available for the topology synthesis on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getSinglePrecisionAtomUpdateKit");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
ISWorkspaceKit<double>
DynamicsIntervention::getDoublePrecisionISWorkspace(const CoordinateCycle orientation,
                                                    const HybridTargetLevel tier) {
  const size_t n_options = implicit_solvent_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (implicit_solvent_cycle_stages[i] == orientation &&
        implicit_solvent_abstract_tiers[i] == tier) {
      return d_iswks[i];
    }
  }
  rtErr("No implicit solvent workspace kit was available on the " + getEnumerationName(tier) +
        " at stage " + getEnumerationName(orientation) + ".", "DynamicsIntervention",
        "getDoublePrecisionISWorkspace");
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
ISWorkspaceKit<float>
DynamicsIntervention::getSinglePrecisionISWorkspace(const CoordinateCycle orientation,
                                                    const HybridTargetLevel tier) {
  const size_t n_options = implicit_solvent_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (implicit_solvent_cycle_stages[i] == orientation &&
        implicit_solvent_abstract_tiers[i] == tier) {
      return f_iswks[i];
    }
  }
  rtErr("No implicit solvent workspace kit was available on the " + getEnumerationName(tier) +
        " at stage " + getEnumerationName(orientation) + ".", "DynamicsIntervention",
        "getSinglePrecisionISWorkspace");
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
const LocalExclusionMaskReader
DynamicsIntervention::getLocalExclusionMaskData(const HybridTargetLevel tier) const {
  const size_t n_options = exclusion_mask_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (exclusion_mask_abstract_tiers[i] == tier) {
      return lemrs[i];
    }
  }
  rtErr("No exclusion mask abstract was available for the synthesis of systems on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getLocalExclusionMaskData");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const SeMaskSynthesisReader
DynamicsIntervention::getStaticExclusionMaskData(const HybridTargetLevel tier) const {
  const size_t n_options = exclusion_mask_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (exclusion_mask_abstract_tiers[i] == tier) {
      return semrs[i];
    }
  }
  rtErr("No exclusion mask abstract was available for the synthesis of systems on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getStaticExclusionMaskData");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
size_t DynamicsIntervention::getNeighborListCoordinateType() const {
  return ngbr_list_tcoord;
}

//-------------------------------------------------------------------------------------------------
size_t DynamicsIntervention::getNeighborListAccumulationType() const {
  return ngbr_list_tacc;
}

//-------------------------------------------------------------------------------------------------
size_t DynamicsIntervention::getNeighborListCalculationType() const {
  return ngbr_list_tcalc;
}

//-------------------------------------------------------------------------------------------------
size_t DynamicsIntervention::getNeighborListTupleType() const {
  return ngbr_list_tcoord4;
}

//-------------------------------------------------------------------------------------------------
NeighborListKind DynamicsIntervention::getNeighborListLayout() const {
  return nl_layout;
}

//-------------------------------------------------------------------------------------------------
ScoreCardWriter DynamicsIntervention::getEnergyTrackingData(const HybridTargetLevel tier) {
  const size_t n_options = utracking_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (utracking_abstract_tiers[i] == tier) {
      return scws[i];
    }
  }
  rtErr("No energy tracking abstract was available for the calculation on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getEnergyTrackingData");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
WatcherWriter DynamicsIntervention::getAnomalyReportingData(const HybridTargetLevel tier) {
  const size_t n_options = anomaly_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    if (anomaly_abstract_tiers[i] == tier) {
      return bugws[i];
    }
  }
  rtErr("No anomaly reporting abstract was available for the calculation on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getAnomalyReportingData");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
HBondWriter DynamicsIntervention::getHydrogenBondAnalysisData(const int index,
                                                              const HybridTargetLevel tier) {
  validateAnalysisIndex<HydrogenBondAnalysis>(index, hbond_trk_ptr,
                                              "getHydrogenBondAnalysisData");
  for (int i = 2 * index; i < 2 * (index + 1); i++) {
    if (hbond_abstract_tiers[i] == tier) {
      return hbws[i];
    }
  }
  rtErr("No hydrogen bond analysis abstract was available for the calculation on the " +
        getEnumerationName(tier) + ".", "DynamicsIntervention", "getHydrogenBondAnalysisData");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceSynthesis* DynamicsIntervention::getPhaseSpaceSynthesisPointer() const {
  if (poly_ps_ptr == nullptr) {
    rtErr("No coordinate synthesis has yet been attached.", "DynamicsIntervention",
          "getPhaseSpaceSynthesisPointer");
  }
  return poly_ps_ptr;
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis* DynamicsIntervention::getPhaseSpaceSynthesisPointer() {
  if (poly_ps_ptr == nullptr) {
    rtErr("No coordinate synthesis has yet been attached.", "DynamicsIntervention",
          "getPhaseSpaceSynthesisPointer");
  }
  return poly_ps_ptr;
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis* DynamicsIntervention::getWorkspacePointer(const size_t index) {
  if (workspaces.size() <= index) {
    if (poly_ps_ptr != nullptr) {
      const int nsys = poly_ps_ptr->getSystemCount();
      if (workspaces.size() == 0) {
        rtErr("No workspace has been allocated to hold an extra copy of the coordinates, "
              "velocities, and forces associated with the synthesis of " + std::to_string(nsys) +
              " systems.", "DynamicsIntervention", "getWorkspacePointer");
      }
      else {
        rtErr("Workspace index " + std::to_string(index) + " is invalid for a set of " +
              std::to_string(workspaces.size()) + " extra copies of the coordinates, velocities, "
              "and forces associated with the synthesis of " + std::to_string(nsys) +
              " systems.", "DynamicsIntervention", "getWorkspacePointer");
      }
    }
    else {
      if (workspaces.size() == 0) {
        rtErr("No workspace has been allocated to hold an extra copy of the coordinates, "
	      "velocities, and forces.  Furthermore, no coordinate synthesis of the original "
              "problem is attached.", "DynamicsIntervention", "getWorkspacePointer");
      }
      else {
        rtErr("Workspace index " + std::to_string(index) + " is invalid for a set of " +
              std::to_string(workspaces.size()) + " extra copies of the coordinates, velocities, "
              "and forces.  No synthesis of systems is, in fact, attached.",
              "DynamicsIntervention", "getWorkspacePointer");
      }
    }
  }
  return workspaces[index];
}

//-------------------------------------------------------------------------------------------------
Condensate* DynamicsIntervention::getCondensatePointer() {
  if (cdns_ptr == nullptr) {
    rtErr("No coordinate synthesis has yet been attached.", "DynamicsIntervention",
          "getCondensatePointer");
  }
  return cdns_ptr;
}

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis* DynamicsIntervention::getTopologySynthesisPointer() {
  if (poly_ag_ptr == nullptr) {
    rtErr("No topology synthesis has yet been attached.", "DynamicsIntervention",
          "getAtomGraphSynthesisPointer");
  }
  return poly_ag_ptr;
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventWorkspace* DynamicsIntervention::getImplicitSolventWorkspacePointer() {
  if (isw_ptr == nullptr) {
    rtErr("No implicit solvent workspace has yet been attached.", "DynamicsIntervention",
          "getImplicitSolventWorkspacePointer");
  }
  return isw_ptr;
}
  
//-------------------------------------------------------------------------------------------------
const LocalExclusionMask* DynamicsIntervention::getLocalExclusionMaskPointer() const {
  if (lem_ptr == nullptr) {
    rtErr("No set of local exclusion masks has yet been attached.", "DynamicsIntervention",
          "getExclusionMaskPointer");
  }
  return lem_ptr;
};

//-------------------------------------------------------------------------------------------------
const StaticExclusionMaskSynthesis* DynamicsIntervention::getStaticExclusionMaskPointer() const {
  if (se_ptr == nullptr) {
    rtErr("No set of static exclusion masks has yet been attached.", "DynamicsIntervention",
          "getStaticExclusionMaskPointer");
  }
  return se_ptr;
};

//-------------------------------------------------------------------------------------------------
ScoreCard* DynamicsIntervention::getEnergyTrackingPointer() {
  if (sc_ptr == nullptr) {
    rtErr("No energy tracking object has yet been assigned to record data for interventions.",
          "DynamicsIntervention", "getEnergyTracking");
  }
  return sc_ptr;
}

//-------------------------------------------------------------------------------------------------
const ScoreCard* DynamicsIntervention::getEnergyTrackingPointer() const {
  if (sc_ptr == nullptr) {
    rtErr("No energy tracking object has yet been assigned to record data for interventions.",
          "DynamicsIntervention", "getEnergyTracking");
  }
  return sc_ptr;
}

//-------------------------------------------------------------------------------------------------
Watcher* DynamicsIntervention::getAnomalyReportingPointer() {
  if (anom_ptr == nullptr) {
    rtErr("No anomaly reporting object has yet been attached.", "DynamicsIntervention",
          "getAnomalyReporting");
  }
  return anom_ptr;
}

//-------------------------------------------------------------------------------------------------
const Watcher* DynamicsIntervention::getAnomalyReportingPointer() const {
  if (anom_ptr == nullptr) {
    rtErr("No anomaly reporting object has yet been attached.", "DynamicsIntervention",
          "getAnomalyReporting");
  }
  return anom_ptr;
}

//-------------------------------------------------------------------------------------------------
int DynamicsIntervention::getHydrogenBondAnalysisCount() const {
  return hbond_trk_ptr.size();
}

//-------------------------------------------------------------------------------------------------
HydrogenBondAnalysis* DynamicsIntervention::getHydrogenBondAnalysisPointer(const int index) {
  validateAnalysisIndex<HydrogenBondAnalysis>(index, hbond_trk_ptr,
                                              "getHydrogenBondAnalysisPointer");
  return hbond_trk_ptr[index];
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::uploadDebugging() {
  if (this->hasAnomalyReporting()) {
    anom_ptr->upload();
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::downloadDebugging() {
  if (this->hasAnomalyReporting()) {
    anom_ptr->download();
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::uploadAnalysis() {

  // Hydrogen bond analyses
  for (size_t i = 0; i < hbond_trk_ptr.size(); i++) {
    hbond_trk_ptr[i]->upload();
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::downloadAnalysis() {

  // Hydrogen bond analyses
  for (size_t i = 0; i < hbond_trk_ptr.size(); i++) {
    hbond_trk_ptr[i]->download();
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::uploadAnalysisSetup() {

  // Hydrogen bond analyses
  for (size_t i = 0; i < hbond_trk_ptr.size(); i++) {
    hbond_trk_ptr[i]->uploadDefinitions();
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::downloadAnalysisSetup() {

  // Hydrogen bond analyses
  for (size_t i = 0; i < hbond_trk_ptr.size(); i++) {
    hbond_trk_ptr[i]->downloadDefinitions();
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::uploadAnalysisData() {

  // Hydrogen bond analyses
  for (size_t i = 0; i < hbond_trk_ptr.size(); i++) {
    hbond_trk_ptr[i]->uploadAccumulators();
  }
}
  
//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::downloadAnalysisData() {

  // Hydrogen bond analyses
  for (size_t i = 0; i < hbond_trk_ptr.size(); i++) {
    hbond_trk_ptr[i]->downloadAccumulators();
  }
}
#endif

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::setGpu(const GpuDetails &gpu_in) {
  gpu = gpu_in;
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::setCoordinateSynthesis(PhaseSpaceSynthesis *poly_ps_in) {

  // Clear the existing buffers, if there is any allocation.  The abstracts contain const-qualified
  // members, which cannot be moved once they are constructed.
  poly_psws.clear();
  poly_psrs.clear();

  // The coordinate synthesis may take different formats, implying the presence of data on
  // different components of the machine.
  std::vector<HybridTargetLevel> available_tiers;
#ifdef STORMM_USE_HPC
  switch (poly_ps_in->getFormat()) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
    available_tiers = { HybridTargetLevel::HOST, HybridTargetLevel::DEVICE };
    break;
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
    available_tiers = { HybridTargetLevel::HOST };
    break;
  case HybridFormat::DEVICE_ONLY:
    available_tiers = { HybridTargetLevel::DEVICE };
    break;
  }
#else
  available_tiers = { HybridTargetLevel::HOST };
#endif
  const size_t n_options = available_tiers.size();
  poly_psws.reserve(2 * n_options);
  poly_psrs.reserve(2 * n_options);
  coordinate_abstract_tiers.reserve(2 * n_options);
  const std::vector<CoordinateCycle> all_stages = {
    poly_ps_in->getCyclePosition(),
    getNextCyclePosition(poly_ps_in->getCyclePosition())
  };
  for (size_t i = 0; i < n_options; i++) {
    const HybridTargetLevel itier = available_tiers[i];
    for (size_t j = 0; j < all_stages.size(); j++) {
      coordinate_abstract_tiers.push_back(itier);
      coordinate_cycle_stages.push_back(all_stages[j]);
      PsSynthesisWriter t_psw = poly_ps_in->data(all_stages[j], itier);
      poly_psws.push_back(t_psw);
      poly_psrs.emplace_back(t_psw);
    }
  }
  poly_ps_ptr = poly_ps_in;
  
  // Check against any pre-established topology synthesis and other resources
  const int nsys = poly_ps_ptr->getSystemCount();
  std::vector<int> natom_v(nsys);
  for (int i = 0; i < nsys; i++) {
    natom_v[i] = poly_ps_ptr->getAtomCount(i);
  }
  validateContents(natom_v, AttachmentKind::PHASE_SPACE_SYNTHESIS, "setCoordinateSynthesis");
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::setCoordinateSynthesis(Condensate *cdns_in) {

  // The implementation parallels that for the PhaseSpaceSynthesis, except that with a Condensate
  // there is no conept of a coordinate time cycle.
  cdnsws.clear();
#ifdef STORMM_USE_HPC
  switch (cdns_in->getFormat()) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
    condensate_abstract_tiers = { HybridTargetLevel::HOST, HybridTargetLevel::DEVICE };
    break;
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
    condensate_abstract_tiers = { HybridTargetLevel::HOST };
    break;
  case HybridFormat::DEVICE_ONLY:
    condensate_abstract_tiers = { HybridTargetLevel::DEVICE };
    break;
  }
#else
  condensate_abstract_tiers = { HybridTargetLevel::HOST };
#endif
  const size_t n_options = condensate_abstract_tiers.size();
  cdnsws.reserve(n_options);
  for (size_t i = 0; i < n_options; i++) {
    cdnsws.push_back(cdns_in->data(condensate_abstract_tiers[i]));
  }
  cdns_ptr = cdns_in;
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::setTopologySynthesis(AtomGraphSynthesis *poly_ag_in) {

  // Clear the existing buffers, if there is any allocation.  The abstracts contain const-qualified
  // members, which will not play well with dynamic array resizing / growth.
  dpoly_vks.clear();
  fpoly_vks.clear();
  dpoly_nbks.clear();
  fpoly_nbks.clear();
  dpoly_rks.clear();
  fpoly_rks.clear();
  dpoly_auks.clear();
  fpoly_auks.clear();

  // Topology abstracts always come with memory at all available tiers.  This may change in a
  // future release.
#ifdef STORMM_USE_HPC
  topology_abstract_tiers = { HybridTargetLevel::HOST, HybridTargetLevel::DEVICE };
#else
  topology_abstract_tiers = { HybridTargetLevel::HOST };
#endif
  const size_t n_options = topology_abstract_tiers.size();
  dpoly_vks.reserve(n_options);
  fpoly_vks.reserve(n_options);
  dpoly_nbks.reserve(n_options);
  fpoly_nbks.reserve(n_options);
  dpoly_rks.reserve(n_options);
  fpoly_rks.reserve(n_options);
  dpoly_auks.reserve(n_options);
  fpoly_auks.reserve(n_options);
  for (size_t i = 0; i < n_options; i++) {
    dpoly_vks.push_back(poly_ag_in->getDoublePrecisionValenceKit(topology_abstract_tiers[i]));
    fpoly_vks.push_back(poly_ag_in->getSinglePrecisionValenceKit(topology_abstract_tiers[i]));
    dpoly_nbks.push_back(poly_ag_in->getDoublePrecisionNonbondedKit(topology_abstract_tiers[i]));
    fpoly_nbks.push_back(poly_ag_in->getSinglePrecisionNonbondedKit(topology_abstract_tiers[i]));
    dpoly_rks.push_back(poly_ag_in->getDoublePrecisionRestraintKit(topology_abstract_tiers[i]));
    fpoly_rks.push_back(poly_ag_in->getSinglePrecisionRestraintKit(topology_abstract_tiers[i]));
    dpoly_auks.push_back(poly_ag_in->getDoublePrecisionAtomUpdateKit(topology_abstract_tiers[i]));
    fpoly_auks.push_back(poly_ag_in->getSinglePrecisionAtomUpdateKit(topology_abstract_tiers[i]));
  }
  poly_ag_ptr = poly_ag_in;

  // Check that the topology synthesis plausibly matches the coordinate synthesis and any
  // accessories.  Reciprocal checks are included when any of these other resources are attached
  // to the intervention object.
  const int nsys = poly_ag_ptr->getSystemCount();
  std::vector<int> natom_v(nsys);
  for (int i = 0; i < nsys; i++) {
    natom_v[i] = poly_ag_ptr->getAtomCount(i);
  }
  validateContents(natom_v, AttachmentKind::ATOM_GRAPH_SYNTHESIS, "setTopologySynthesis");
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::setImplicitSolventWorkspace(ImplicitSolventWorkspace *isw) {

  // Clear any existing buffers
  d_iswks.clear();
  f_iswks.clear();
  implicit_solvent_cycle_stages.clear();
  implicit_solvent_abstract_tiers.clear();
  
  // Implicit solvent workspaces will always be prepared in the premium memory format.
#ifdef STORMM_USE_HPC
  const std::vector<HybridTargetLevel> all_tiers = { HybridTargetLevel::HOST,
                                                     HybridTargetLevel::DEVICE };
#else
  const std::vector<HybridTargetLevel> all_tiers = { HybridTargetLevel::HOST };
#endif
  const size_t n_options = all_tiers.size();
  d_iswks.reserve(2 * n_options);
  f_iswks.reserve(2 * n_options);
  implicit_solvent_cycle_stages.reserve(2 * n_options);
  implicit_solvent_abstract_tiers.reserve(2 * n_options);
  for (size_t i = 0; i < n_options; i++) {
    const CoordinateCycle current_stage = isw->getCyclePosition();
    const CoordinateCycle next_stage = getNextCyclePosition(current_stage);
    d_iswks.push_back(isw->dpData(current_stage, all_tiers[i]));
    d_iswks.push_back(isw->dpData(next_stage, all_tiers[i]));
    f_iswks.push_back(isw->spData(current_stage, all_tiers[i]));
    f_iswks.push_back(isw->spData(next_stage, all_tiers[i]));
    implicit_solvent_cycle_stages.push_back(current_stage);
    implicit_solvent_cycle_stages.push_back(next_stage);
    implicit_solvent_abstract_tiers.push_back(all_tiers[i]);
    implicit_solvent_abstract_tiers.push_back(all_tiers[i]);
  }

  // The implicit solvent workspace has only an indicator as to its overall size, no system
  // boundaries.  Check this against any attached synthesis objects.
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::setExclusionMasks(const LocalExclusionMask &lem) {
  
  // Clear any existing buffer
  lemrs.clear();

  // Like the topology they are based upon, exclusion masks are always prepared with the premium
  // memory format, availabel on both the CPU host and GPU device with optimized transfers
  // between each resource.
#ifdef STORMM_USE_HPC
  exclusion_mask_abstract_tiers = { HybridTargetLevel::HOST, HybridTargetLevel::DEVICE };
#else
  exclusion_mask_abstract_tiers = { HybridTargetLevel::HOST };
#endif
  const size_t n_options = exclusion_mask_abstract_tiers.size();
  lemrs.reserve(n_options);
  for (size_t i = 0; i < n_options; i++) {
    lemrs.push_back(lem.data(exclusion_mask_abstract_tiers[i]));
  }

  // Store a pointer to the original object.
  lem_ptr = const_cast<LocalExclusionMask*>(lem.getSelfPointer());

  // Check that the set of exclusion masks corresponds to the topology synthesis loaded into the
  // object.
  const AtomGraphSynthesis *tpag = lem.getTopologySynthesisPointer();
  if (poly_ag_ptr != nullptr && tpag != poly_ag_ptr) {
    rtErr("A request was made to load a LocalExclusionMask object based on a topology synthesis "
          "which is inconsistent with the topology synthesis upon which this object is based.  "
          "The LocalExclusionMask serves " + std::to_string(tpag->getSystemCount()) +
          " systems, and begins with a system based on topology " +
          getBaseName(tpag->getSystemTopologyPointer(0)->getFileName()) + ".  In contrast, this "
          "object's underlying topology synthesis contains " +
          std::to_string(poly_ag_ptr->getSystemCount()) + " systems, beginning with one based on "
          "topology " + getBaseName(poly_ag_ptr->getSystemTopologyPointer(0)->getFileName()) + ".",
          "DynamicsIntervention", "setExclusionMasks");
  }
  validateContents(tpag->getSystemCount(), "setExclusionMasks");
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::setExclusionMasks(const StaticExclusionMaskSynthesis &se) {
  
  // Clear any existing buffer
  semrs.clear();

  // Both types of exclusion masks make use of the same array to mark their memory tiers.  For any
  // given simulation, only one abstract will be relevant.
#ifdef STORMM_USE_HPC
  exclusion_mask_abstract_tiers = { HybridTargetLevel::HOST, HybridTargetLevel::DEVICE };
#else
  exclusion_mask_abstract_tiers = { HybridTargetLevel::HOST };
#endif
  const size_t n_options = exclusion_mask_abstract_tiers.size();
  semrs.reserve(n_options);
  for (size_t i = 0; i < n_options; i++) {
    semrs.push_back(se.data(exclusion_mask_abstract_tiers[i]));
  }

  // Store a pointer to the original object.
  se_ptr = const_cast<StaticExclusionMaskSynthesis*>(se.getSelfPointer());

  // Check that the set of exclusion masks corresponds to the topology synthesis loaded into the
  // object.
  const int nsys = se.getSystemCount();
  if ((poly_ag_ptr != nullptr && nsys != poly_ag_ptr->getSystemCount()) ||
      (poly_ps_ptr != nullptr && nsys != poly_ag_ptr->getSystemCount())) {
    const int expected_nsys = (poly_ag_ptr != nullptr) ? poly_ag_ptr->getSystemCount() :
                                                         poly_ps_ptr->getSystemCount();
    rtErr("A request was made to load a StaticExclusionMaskSynthesis object based on a collection "
          "of systems which cannot be consistent with the topology synthesis upon which this "
          "object is based.  The StaticExclusionMask serves " +
          std::to_string(nsys) + " systems, whereas this object's underlying synthesis contains " +
          std::to_string(expected_nsys) + " systems.", "DynamicsIntervention",
          "setExclusionMasks");
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::setEnergyTracking(ScoreCard *sc_in) {

  // Clear any existing buffer
  scws.clear();

  // ScoreCard objects are created with the premium memory layout for rapid communication to and
  // from the GPU device.
#ifdef STORMM_USE_HPC
  utracking_abstract_tiers = { HybridTargetLevel::HOST, HybridTargetLevel::DEVICE };
#else
  utracking_abstract_tiers = { HybridTargetLevel::HOST };
#endif
  const size_t n_options = utracking_abstract_tiers.size();
  scws.reserve(n_options);
  for (size_t i = 0; i < n_options; i++) {
    scws.push_back(sc_in->data(utracking_abstract_tiers[i]));
  }
  sc_ptr = sc_in;
}
  
//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::setAnomalyReporting(Watcher *anom_in) {
  active_interventions = true;

  // Clear any existing buffer
  bugws.clear();
  
  // Watcher objects are created with the premium memory layout, for maximum efficiency in
  // transferring information from the GPU device to the CPU host, which then decides how to
  // proceed when an anomaly is detected.
#ifdef STORMM_USE_HPC
  anomaly_abstract_tiers = { HybridTargetLevel::HOST, HybridTargetLevel::DEVICE };
#else
  anomaly_abstract_tiers = { HybridTargetLevel::HOST };
#endif
  const size_t n_options = anomaly_abstract_tiers.size();
  bugws.reserve(n_options);
  for (size_t i = 0; i < n_options; i++) {
    bugws.push_back(anom_in->data(anomaly_abstract_tiers[i]));
  }
  anom_ptr = anom_in;
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::setAnalysis(HydrogenBondAnalysis *attachment,
                                       const DynamicsControls &dyncon,
                                       const HybridTargetLevel tier) {
  active_interventions = true;
  hbond_trk_ptr.push_back(attachment);
  hbws.clear();
  const int n_hba = hbond_trk_ptr.size();
#ifdef STORMM_USE_HPC
  hbws.reserve(2 * n_hba);
  hbond_abstract_tiers.reserve(2 * n_hba);
  for (int i = 0; i < n_hba; i++) {
    hbond_abstract_tiers.push_back(HybridTargetLevel::HOST);
    hbond_abstract_tiers.push_back(HybridTargetLevel::DEVICE);
  }
#else
  hbws.reserve(n_hba);
  hbond_abstract_tiers = std::vector<HybridTargetLevel>(n_hba, HybridTargetLevel::HOST);
#endif
  for (int i = 0; i < n_hba; i++) {
#ifdef STORMM_USE_HPC
    for (int j = 0; j < 2; j++) {
      hbws.push_back(hbond_trk_ptr[i]->data(hbond_abstract_tiers[(2 * i) + j]));
    }
#else
    hbws.push_back(hbond_trk_ptr[i]->data(hbond_abstract_tiers[i]));
#endif
  }

  // Add the evaluator to the appropriate queue.  Hydrogen bond evaluation will occur after
  // coordinates have been updated, which may include application of geometric constraints.
  IntegrationStage when;
  switch (dyncon.constrainGeometry()) {
  case ApplyConstraints::YES:
    when = IntegrationStage::GEOMETRY_CONSTRAINT;
    break;
  case ApplyConstraints::NO:
    when = IntegrationStage::POSITION_ADVANCE;
    break;
  }
  addAction(evalHydrogenBondAnalysis, when, attachment->getEvaluationInterval(), n_hba - 1, tier);
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::setDebugging(Watcher *anom, const DebugControls &dbgcon,
                                        const ApplyConstraints enforce_constraints,
                                        PhaseSpaceSynthesis *workspace_in,
                                        const HybridTargetLevel tier) {

  // Transfer the user directives from the input file to this object in the form of queued actions
  // in the private arrays of function pointers.
  active_interventions = true;
  setAnomalyReporting(anom);
  switch (workspace_in->getFormat()) {
  case HybridFormat::HOST_ONLY:
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::HOST_MOUNTED:
    break;
  case HybridFormat::UNIFIED:
    rtErr("A workspace with " + getEnumerationName(workspace_in->getFormat()) + " is not accepted "
          "for debugging due to the way the Page Migration Engine performs constant "
          "synchronization between memory on the CPU host and GPU device.  Debugging is intended "
          "to keep scratch work on the host.", "DynamicsIntervention", "setDebugging");
#endif
  }
  int work_idx = (workspaces.size() > 0) ? workspaces.size() - 1 : 0;
  if (workspace_in != nullptr) {
    workspaces.push_back(workspace_in);
    work_idx = workspaces.size() - 1;
  }

  // Checks on the neighbor list will be queued here, but require the attachment of a workspace
  // or workspaces (for CellGrid objects) tailored to the neighbor list in use.  That workspace
  // may not exist until the neighbor lists themselves are constructed, and so will not be
  // required as a formal argument to this function.
  if (dbgcon.checkNeighborListComp()) {
    check_neighbor_list = true;
    switch (poly_ps_ptr->getUnitCellType()) {
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      addAction(execNLCompositionDebug, IntegrationStage::CALC_FORCES, 1, 0, tier);
      break;
    case UnitCellType::NONE:
      rtErr("Neighbor lists do not exist for systems without periodic boundary conditions.",
            "DynamicsIntervention", "setDebugging");
    }
  }
  if (dbgcon.checkForces()) {
    switch (poly_ps_ptr->getUnitCellType()) {
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      addAction(execNLForceDebug, IntegrationStage::CALC_FORCES, dbgcon.getInspectionInterval(),
                0, tier);
      break;
    case UnitCellType::NONE:
      break;
    }
    addAction(execForceDebug, IntegrationStage::CALC_FORCES, dbgcon.getInspectionInterval(),
              work_idx, tier);
  }
  if (dbgcon.reportLargeForces()) {
    investigate_large_forces = true;
    addAction(evalForceAnomalies, IntegrationStage::CALC_FORCES, 1, work_idx, tier);
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::setDebugging(Watcher *anom, const DebugControls &dbgcon,
                                        const DynamicsControls &dyncon,
                                        PhaseSpaceSynthesis *workspace_in,
                                        const HybridTargetLevel tier) {
  setDebugging(anom, dbgcon, dyncon.constrainGeometry(), workspace_in, tier);
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::addAction(DynIntvFuncPtr action, const IntegrationStage process,
                                     const int interval, const int attachment_index,
                                     const HybridTargetLevel tier) {
  switch (process) {
  case IntegrationStage::CALC_FORCES:
    force_calc_actions.push_back(action);
    force_calc_action_intv.push_back(interval);
    force_calc_action_index.push_back(attachment_index);
    force_calc_action_tier.push_back(tier);
    break;
  case IntegrationStage::VELOCITY_ADVANCE:
    velocity_adv_actions.push_back(action);
    velocity_adv_action_intv.push_back(interval);
    velocity_adv_action_index.push_back(attachment_index);
    velocity_adv_action_tier.push_back(tier);
    break;
  case IntegrationStage::VELOCITY_CONSTRAINT:
    velocity_cnst_actions.push_back(action);
    velocity_cnst_action_intv.push_back(interval);
    velocity_cnst_action_index.push_back(attachment_index);
    velocity_cnst_action_tier.push_back(tier);
    break;
  case IntegrationStage::CALC_KINETIC:
    kinetic_calc_actions.push_back(action);
    kinetic_calc_action_intv.push_back(interval);
    kinetic_calc_action_index.push_back(attachment_index);
    kinetic_calc_action_tier.push_back(tier);
    break;
  case IntegrationStage::POSITION_ADVANCE:
    position_adv_actions.push_back(action);
    position_adv_action_intv.push_back(interval);
    position_adv_action_index.push_back(attachment_index);
    position_adv_action_tier.push_back(tier);
    break;
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    geometry_cnst_actions.push_back(action);
    geometry_cnst_action_intv.push_back(interval);
    geometry_cnst_action_index.push_back(attachment_index);
    geometry_cnst_action_tier.push_back(tier);
    break;
  }
  setActiveInterval(interval);
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::addAction(DynIntvFuncPtr action, const IntegrationStage process,
                                     const HybridTargetLevel tier) {
  addAction(action, process, 1, 0, tier);
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::execute(const int step, const IntegrationStage process) {
  switch (process) {
  case IntegrationStage::CALC_FORCES:
    {
      const size_t n_actions = force_calc_actions.size();
      for (size_t i = 0; i < n_actions; i++) {
        if (step % force_calc_action_intv[i] == 0) {
          force_calc_actions[i](step, force_calc_action_index[i], force_calc_action_tier[i]);
        }
      }
    }
    break;
  case IntegrationStage::VELOCITY_ADVANCE:
    {
      const size_t n_actions = velocity_adv_actions.size();
      for (size_t i = 0; i < n_actions; i++) {
        if (step % velocity_adv_action_intv[i] == 0) {
          velocity_adv_actions[i](step, velocity_adv_action_index[i], velocity_adv_action_tier[i]);
        }
      }
    }
    break;
  case IntegrationStage::VELOCITY_CONSTRAINT:
    {
      const size_t n_actions = velocity_cnst_actions.size();
      for (size_t i = 0; i < n_actions; i++) {
        if (step % velocity_cnst_action_intv[i] == 0) {
          velocity_cnst_actions[i](step, velocity_cnst_action_index[i],
                                   velocity_cnst_action_tier[i]);
        }
      }
    }
    break;
  case IntegrationStage::CALC_KINETIC:
    {
      const size_t n_actions = kinetic_calc_actions.size();
      for (size_t i = 0; i < n_actions; i++) {
        if (step % kinetic_calc_action_intv[i] == 0) {
          kinetic_calc_actions[i](step, kinetic_calc_action_index[i], kinetic_calc_action_tier[i]);
        }
      }
    }
    break;
  case IntegrationStage::POSITION_ADVANCE:
    {
      const size_t n_actions = position_adv_actions.size();
      for (size_t i = 0; i < n_actions; i++) {
        if (step % position_adv_action_intv[i] == 0) {
          position_adv_actions[i](step, position_adv_action_index[i], position_adv_action_tier[i]);
        }
      }
    }
    break;
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    {
      const size_t n_actions = geometry_cnst_actions.size();
      for (size_t i = 0; i < n_actions; i++) {
        if (step % geometry_cnst_action_intv[i] == 0) {
          geometry_cnst_actions[i](step, geometry_cnst_action_index[i],
                                   geometry_cnst_action_tier[i]);
        }
      }
    }
    break;
  }
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::upload() {

}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::download() {

}
#endif
  
//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::validateContents(int system_count, const char* caller) const {

  // Check that the energy tracking is prepared to accommodate the synthesis of systems
  if (poly_ag_ptr != nullptr && system_count != poly_ag_ptr->getSystemCount()) {
    rtErr("The new attachment is set up for " + std::to_string(system_count) + " systems, while "
          "the attached topology synthesis contains " +
          std::to_string(poly_ag_ptr->getSystemCount()) + " systems.", "DynamicsIntervention",
          caller);
  }
  if (poly_ps_ptr != nullptr && system_count != poly_ps_ptr->getSystemCount()) {
    rtErr("The new attachment is set up for " + std::to_string(system_count) + " systems, while "
          "the attached coordinate synthesis (PhaseSpaceSynthesis) contains " +
          std::to_string(poly_ps_ptr->getSystemCount()) + " systems.", "DynamicsIntervention",
          caller);
  }
  if (cdns_ptr != nullptr && system_count != cdns_ptr->getSystemCount()) {
    rtErr("The new attachment is set up for " + std::to_string(system_count) + " systems, while "
          "the attached coordinate synthesis (Condensate) contains " +
          std::to_string(cdns_ptr->getSystemCount()) + " systems.", "DynamicsIntervention",
          caller);
  }
  if (lem_ptr != nullptr) {
    const int lem_nsys = lem_ptr->getTopologySynthesisPointer()->getSystemCount();
    if (lem_nsys != system_count) {
      rtErr("The new attachment is set up for " + std::to_string(system_count) + " systems, while "
            "the attached (local) exclusion masks apply to " + std::to_string(lem_nsys) +
            " systems.", "DynamicsIntervention", caller);
    }
  }
  if (se_ptr != nullptr) {
    const int se_nsys = se_ptr->getSystemCount();
    if (se_nsys != system_count) {
      rtErr("The new attachment is set up for " + std::to_string(system_count) + " systems, while "
            "the attached (static) exclusion masks apply to " + std::to_string(se_nsys) +
            " systems.", "DynamicsIntervention", caller);
    }
  }
  if (sc_ptr != nullptr && system_count != sc_ptr->getSystemCount()) {
    rtErr("The new attachment is set up for " + std::to_string(system_count) + " systems, while "
          "the attached energy tracking object accommodates " +
          std::to_string(sc_ptr->getSystemCount()) + " systems.", "DynamicsIntervention", caller);
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::validateContents(const int padded_atom_count, const AttachmentKind atch,
                                            const char* caller) const {

  // Check for existing coordinate or topology syntheses (these attachments will check each other
  // for exact system sizes, but not overall atom counts across all systems).
  switch (atch) {
  case AttachmentKind::IMPLICIT_SOLVENT_WSPC:
    if (poly_ps_ptr != nullptr && poly_ps_ptr->getPaddedAtomCount() != padded_atom_count) {
      rtErr("The padded atom counts of the attached PhaseSpaceSynthesis object and the "
            "ImplicitSolventWorkspace do not match (" +
            std::to_string(poly_ps_ptr->getPaddedAtomCount()) + " vs " +
            std::to_string(isw_ptr->getPaddedAtomCount()) + ").", "DynamicsIntervention",
            "setImplicitSolventWorkspace");
    }
    if (cdns_ptr != nullptr && cdns_ptr->getPaddedAtomCount() != padded_atom_count) {
      rtErr("The padded atom counts of the attached Condensate object and the "
            "ImplicitSolventWorkspace do not match (" +
            std::to_string(cdns_ptr->getPaddedAtomCount()) + " vs " +
            std::to_string(isw_ptr->getPaddedAtomCount()) + ").", "DynamicsIntervention",
            "setImplicitSolventWorkspace");
    }
    if (poly_ag_ptr != nullptr && poly_ag_ptr->getPaddedAtomCount() != padded_atom_count) {
      rtErr("The padded atom counts of the attached AtomGraphSynthesis object and the "
            "ImplicitSolventWorkspace do not match (" +
            std::to_string(poly_ag_ptr->getPaddedAtomCount()) + " vs " +
            std::to_string(isw_ptr->getPaddedAtomCount()) + ").", "DynamicsIntervention",
            "setImplicitSolventWorkspace");
    }
    break;
  case AttachmentKind::ATOM_GRAPH_SYNTHESIS:
  case AttachmentKind::PHASE_SPACE_SYNTHESIS:
  case AttachmentKind::CONDENSATE:
  case AttachmentKind::LOCAL_EXCLUSIONMASK:
  case AttachmentKind::SCORE_CARD:
  case AttachmentKind::WATCHER:
    break;
  }

  // Check for an existing implicit solvent workspace.
  switch (atch) {
  case AttachmentKind::ATOM_GRAPH_SYNTHESIS:
  case AttachmentKind::PHASE_SPACE_SYNTHESIS:
  case AttachmentKind::CONDENSATE:
    if (isw_ptr != nullptr && padded_atom_count != isw_ptr->getPaddedAtomCount()) {
      rtErr("The padded atom counts of an attached ImplicitSolventWorkspace object and an " +
            getEnumerationName(atch) + " object do not agree (" +
            std::to_string(isw_ptr->getPaddedAtomCount()) + " vs " +
            std::to_string(padded_atom_count) + ").", "DynamicsIntervention",
            "setImplicitSolventWorkspace");
    }
    break;
  case AttachmentKind::IMPLICIT_SOLVENT_WSPC:
  case AttachmentKind::LOCAL_EXCLUSIONMASK:
  case AttachmentKind::SCORE_CARD:
  case AttachmentKind::WATCHER:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::validateContents(const std::vector<int> &system_sizes,
                                            const AttachmentKind atch, const char* caller) const {
  const int nsys = system_sizes.size();
  validateContents(nsys, caller);

  // Check for an existing PhaseSpaceSynthesis attachment
  switch (atch) {
  case AttachmentKind::ATOM_GRAPH_SYNTHESIS:
  case AttachmentKind::CONDENSATE:
    if (poly_ps_ptr != nullptr) {
      for (int i = 0; i < nsys; i++) {
        if (poly_ps_ptr->getAtomCount(i) != system_sizes[i]) {
          rtErr("The attached coordinate synthesis holds " +
                std::to_string(poly_ps_ptr->getAtomCount(i)) + " atoms in system index " +
                std::to_string(i) + ", but the new resource holds " +
                std::to_string(system_sizes[i]) + ".", "DynamicsIntervention", caller);
        }
      }
    }
    break;
  case AttachmentKind::PHASE_SPACE_SYNTHESIS:
  case AttachmentKind::IMPLICIT_SOLVENT_WSPC:
  case AttachmentKind::LOCAL_EXCLUSIONMASK:
  case AttachmentKind::SCORE_CARD:
  case AttachmentKind::WATCHER:
    return;
  }

  // Check for an existing AtomGraphSynthesis attachment
  switch (atch) {
  case AttachmentKind::ATOM_GRAPH_SYNTHESIS:
    break;
  case AttachmentKind::PHASE_SPACE_SYNTHESIS:
  case AttachmentKind::CONDENSATE:
    if (poly_ag_ptr != nullptr) {
      for (int i = 0; i < nsys; i++) {
        if (poly_ag_ptr->getAtomCount(i) != system_sizes[i]) {
          rtErr("The attached topology synthesis holds " +
                std::to_string(poly_ag_ptr->getAtomCount(i)) + " atoms in system index " +
                std::to_string(i) + ", but the new resource holds " +
                std::to_string(system_sizes[i]) + ".", "DynamicsIntervention", caller);
        }
      }
    }
    break;
  case AttachmentKind::IMPLICIT_SOLVENT_WSPC:
  case AttachmentKind::LOCAL_EXCLUSIONMASK:
  case AttachmentKind::SCORE_CARD:
  case AttachmentKind::WATCHER:
    return;
  }

  // Check for an existing Condensate attachment
  switch (atch) {
  case AttachmentKind::ATOM_GRAPH_SYNTHESIS:
  case AttachmentKind::PHASE_SPACE_SYNTHESIS:
    if (cdns_ptr != nullptr) {
      for (int i = 0; i < nsys; i++) {
        if (cdns_ptr->getAtomCount(i) != system_sizes[i]) {
          rtErr("The attached coordinate synthesis holds " +
                std::to_string(cdns_ptr->getAtomCount(i)) + " atoms in system index " +
                std::to_string(i) + ", but the new resource holds " +
                std::to_string(system_sizes[i]) + ".", "DynamicsIntervention", caller);
        }
      }
    }
    break;
  case AttachmentKind::CONDENSATE:
  case AttachmentKind::IMPLICIT_SOLVENT_WSPC:
  case AttachmentKind::LOCAL_EXCLUSIONMASK:
  case AttachmentKind::SCORE_CARD:
  case AttachmentKind::WATCHER:
    return;
  }
}
 
//-------------------------------------------------------------------------------------------------
void DynamicsIntervention::setActiveInterval(const int next_intv) {

  // Treat a call for "zero" intervals as a command to execute on every step.
  if (next_intv == 0) {
    relevant_intervals.resize(1);
    relevant_intervals[0] = 1;
    relevant_interval_count = relevant_intervals.size();
    return;
  }

  // Catch spurious input
  if (next_intv < 0) {
    rtErr("An active interval of " + std::to_string(next_intv) + " is invalid.",
          "DynamicsIntervention", "setActiveInterval");
  }
  
  // Is the interval already covered?
  bool is_covered = false;
  const size_t n = relevant_intervals.size();
  for (size_t i = 0; i < n; i++) {
    is_covered = (is_covered || (next_intv % relevant_intervals[i] == 0));
  }
  if (is_covered == false) {
    relevant_intervals.push_back(next_intv);

    // A new interval has been added.  This interval may, in turn, be a factor of other
    // intervals, making them irrelevant.
    for (size_t i = 0; i < relevant_intervals.size(); i++) {
      while (i < relevant_intervals.size() && relevant_intervals[i] != next_intv &&
             (relevant_intervals[i] % next_intv == 0)) {
        const size_t n_i = relevant_intervals.size() - 1;
        for (size_t j = i; j < n_i; j++) {
          relevant_intervals[j] = relevant_intervals[j + 1];
        }
        relevant_intervals.resize(n_i);
      }
    }
  }

  // Revert to making every step "relevant" to interventions if a very large number of different
  // active frequencies are present.  Any subsequent calls to setActiveInterval() will no longer
  // add to the list.
  if (relevant_intervals.size() >= max_relevant_intervals) {
    relevant_intervals.resize(1);
    relevant_intervals[0] = 1;
  }
  relevant_interval_count = relevant_intervals.size();
}

//-------------------------------------------------------------------------------------------------
DynamicsIntervention dyna_tk;
  
} // namespace mm
} // namespace stormm
