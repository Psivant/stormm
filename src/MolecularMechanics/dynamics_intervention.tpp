// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace mm {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void DynamicsIntervention::validateNeighborListTypes(const char* caller) const {
  if (std::type_index(typeid(Tcoord)).hash_code() != ngbr_list_tcoord) {
    rtErr("The data type of the attached CellGrid object's coordinate transforms is '" +
          getStormmScalarTypeName(ngbr_list_tcoord) + "', while the data type in the restored, "
          "templated abstract is '" + getStormmScalarTypeName<Tcoord>() + "'.",
          "DynamicsIntervention", "getNeighborListData");
  }
  if (std::type_index(typeid(Tacc)).hash_code() != ngbr_list_tacc) {
    rtErr("The data type of the attached CellGrid object's accumulators is '" +
          getStormmScalarTypeName(ngbr_list_tacc) + "', while the data type in the restored, "
          "templated abstract is '" + getStormmScalarTypeName<Tacc>() + "'.",
          "DynamicsIntervention", "getNeighborListData");
  }
  if (std::type_index(typeid(Tcalc)).hash_code() != ngbr_list_tcalc) {
    rtErr("The data type used in the attached CellGrid object's calculations is '" +
          getStormmScalarTypeName(ngbr_list_tcalc) + "', while the data type in the restored, "
          "templated abstract is '" + getStormmScalarTypeName<Tcalc>() + "'.",
          "DynamicsIntervention", "getNeighborListData");
  }
  if (std::type_index(typeid(Tcoord4)).hash_code() != ngbr_list_tcoord4) {
    rtErr("The data type of the attached CellGrid object's coordinate / property tuples is '" +
          getHpcVectorTypeName(ngbr_list_tcoord4) + "', while the data type in the restored, "
          "templated abstract is '" + getHpcVectorTypeName<Tcoord4>() + "'.",
          "DynamicsIntervention", "getNeighborListData");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
CellGridWriter<Tcoord, Tacc, Tcalc, Tcoord4>
DynamicsIntervention::getNeighborListData(const NonbondedTheme theme,
                                          const CoordinateCycle orientation,
                                          const HybridTargetLevel tier) {

  // Check the data types
  validateNeighborListTypes<Tcoord, Tacc, Tcalc, Tcoord4>("getNeighborListData");
  const size_t n_options = voided_cgws.size();
  for (size_t i = 0; i < n_options; i++) {
    if (orientation == neighbor_list_cycle_stages[i] && theme == neighbor_list_themes[i] &&
        tier == neighbor_list_abstract_tiers[i]) {
      return restoreType<Tcoord, Tacc, Tcalc, Tcoord4>(voided_cgws[i]);
    }
  }
  rtErr("No neighbor list abstract was available on the " + getEnumerationName(tier) + " with "
        "theme " + getEnumerationName(theme) + " at the " + getEnumerationName(orientation) +
        " stage of the coordinate cycle.", "DynamicsIntervention", "getNeighborListData");
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
CellGridWriter<Tcoord, Tacc, Tcalc, Tcoord4>
DynamicsIntervention::getNeighborListData(const CoordinateCycle orientation,
                                          const HybridTargetLevel tier) {

  // Check that there is no ambiguity in the neighbor list theme
  const size_t n_options = voided_cgws.size();
  if (n_options == 0) {
    rtErr("No neighbor list (CellGrid) object has been attached.", "DynamicsIntervention",
          "getNeighborListData");
  }
  const NonbondedTheme base_theme = voided_cgws[0].theme;
  for (size_t i = 1; i < n_options; i++) {
    if (voided_cgws[i].theme != base_theme) {
      rtErr("Multiple themes were found within the attached neighbor list abstracts (" +
            getEnumerationName(base_theme) + " and " + getEnumerationName(voided_cgws[i].theme) +
            ").  Specify one using the other accessor function variant in order to obtain a "
            "specific abstract.", "DynamicsIntervention", "getNeighborListData");
    }
  }
  return getNeighborListData<Tcoord, Tacc, Tcalc, Tcoord4>(base_theme, orientation, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>*
DynamicsIntervention::getNeighborListPointer(const NonbondedTheme theme) {
  return getNLAttachment<Tcoord, Tacc, Tcalc, Tcoord4>(theme, &cg_ptr, "getNLWorkspacePointer");
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
const CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>*
DynamicsIntervention::getNeighborListPointer(const NonbondedTheme theme) const {
  return getNLAttachment<Tcoord, Tacc, Tcalc, Tcoord4>(theme, cg_ptr, "getNLWorkspacePointer");
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>*
DynamicsIntervention::getNLWorkspacePointer(const NonbondedTheme theme) {
  return getNLAttachment<Tcoord, Tacc, Tcalc, Tcoord4>(theme, &nl_workspaces,
                                                       "getNLWorkspacePointer");
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
const CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>*
DynamicsIntervention::getNLWorkspacePointer(const NonbondedTheme theme) const {
  return getNLAttachment<Tcoord, Tacc, Tcalc, Tcoord4>(theme, nl_workspaces,
                                                       "getNLWorkspacePointer");
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void DynamicsIntervention::setNeighborList(CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> *cg_a,
                                           CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> *cg_b) {

  // Check that multiple neighbor lists, if provided, are aligned in their cycle
  if (cg_b != nullptr) {
    if (cg_a->getCyclePosition() != cg_b->getCyclePosition()) {
      rtErr("Two neighbor lists were provided, but they differ in cycle positions (" +
            getEnumerationName(cg_a->getCyclePosition()) + " vs " +
            getEnumerationName(cg_a->getCyclePosition()) + ").", "DynamicsIntervention",
            "setNeighborListAbstracts");
    }
    nl_layout = NeighborListKind::DUAL;
  }
  else {
    nl_layout = NeighborListKind::MONO;
  }
  
  // Clear any existing abstracts
  neighbor_list_abstract_tiers.clear();
  neighbor_list_cycle_stages.clear();
  voided_cgws.clear();

  // Account for all memory tiers.  The CellGrid object allocates memory with the premium format,
  // expedited for transfers between the CPU host and GPU device.
#ifdef STORMM_USE_HPC
  const std::vector<HybridTargetLevel> all_tiers = { HybridTargetLevel::HOST,
                                                     HybridTargetLevel::DEVICE };
  const size_t n_abstr = 2 * (1 + (cg_b != nullptr)) * all_tiers.size();
#else
  const std::vector<HybridTargetLevel> all_tiers = { HybridTargetLevel::HOST };
  const size_t n_abstr = (1 + (cg_b != nullptr)) * all_tiers.size();
#endif
  const std::vector<CoordinateCycle> all_stages = {
    cg_a->getCyclePosition(),
    getNextCyclePosition(cg_a->getCyclePosition())
  };
  neighbor_list_abstract_tiers.reserve(n_abstr);
  neighbor_list_cycle_stages.reserve(n_abstr);
  neighbor_list_themes.reserve(n_abstr);
  voided_cgws.reserve(n_abstr);
  const size_t n_options = neighbor_list_abstract_tiers.size();
  for (size_t i = 0; i < n_options; i++) {
    const HybridTargetLevel itier = neighbor_list_abstract_tiers[i]; 
    for (size_t j = 0; j < all_stages.size(); j++) {
      neighbor_list_abstract_tiers.push_back(itier);
      neighbor_list_cycle_stages.push_back(all_stages[j]);
      neighbor_list_themes.push_back(cg_a->getTheme());
      voided_cgws.push_back(cg_a->templateFreeData(all_stages[j], itier));
      if (cg_b != nullptr) {
        neighbor_list_abstract_tiers.push_back(itier);
        neighbor_list_cycle_stages.push_back(all_stages[j]);
        voided_cgws.push_back(cg_b->templateFreeData(all_stages[j], itier));
        neighbor_list_themes.push_back(cg_b->getTheme());
      }
    }
  }

  // Record the relevant data types
  ngbr_list_tcoord = std::type_index(typeid(Tcoord)).hash_code();
  ngbr_list_tacc = std::type_index(typeid(Tacc)).hash_code();
  ngbr_list_tcalc = std::type_index(typeid(Tcalc)).hash_code();
  ngbr_list_tcoord4 = std::type_index(typeid(Tcoord4)).hash_code();

  // Set pointers to the original object or objects
  cg_ptr.resize(2, nullptr);
  if (cg_b == nullptr) {
    switch (cg_a->getTheme()) {
    case NonbondedTheme::ELECTROSTATIC:
    case NonbondedTheme::VAN_DER_WAALS:
      rtErr("If a single neighbor list is to be attached, it must address all non-bonded "
            "interactions.  A neighbor list for " + getEnumerationName(cg_a->getTheme()) +
            " alone is invalid.", "DynamicsIntervention", "setNeighborList");
    case NonbondedTheme::ALL:
      cg_ptr[0] = reinterpret_cast<CellGrid<float, int, float, float4>*>(cg_a);
      break;
    }
  }
  else {
    if (cg_a->getTheme() == cg_b->getTheme()) {
      rtErr("If two neighbor lists are to be attached, they must address different non-bonded "
            "interactions. (Current attachments are both for " +
            getEnumerationName(cg_a->getTheme()) + ".)", "DynamicsIntervention",
            "setNeighborList");
    }
    bool one_uses_all = false;
    switch (cg_a->getTheme()) {
    case NonbondedTheme::ELECTROSTATIC:
      cg_ptr[0] = reinterpret_cast<CellGrid<float, int, float, float4>*>(cg_a);
      cg_ptr[1] = reinterpret_cast<CellGrid<float, int, float, float4>*>(cg_b);
      break;
    case NonbondedTheme::VAN_DER_WAALS:
      cg_ptr[1] = reinterpret_cast<CellGrid<float, int, float, float4>*>(cg_a);
      cg_ptr[0] = reinterpret_cast<CellGrid<float, int, float, float4>*>(cg_b);
      break;
    case NonbondedTheme::ALL:
      one_uses_all = true;
      break;
    }
    switch (cg_b->getTheme()) {
    case NonbondedTheme::ELECTROSTATIC:
    case NonbondedTheme::VAN_DER_WAALS:
    case NonbondedTheme::ALL:
      one_uses_all = true;
      break;
    }
    if (one_uses_all) {
      rtErr("If two neighbor lists are to be attached, one must address the electostatic "
            "interactions while the other addresses van-der Waals interactions, exclusively.  A "
            "neighbor list addressing " + getEnumerationName(cg_a->getTheme()) + " interactions "
            "is invalid.", "DynamicsIntervention", "setNeighborList");
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void DynamicsIntervention::setNLWorkspace(CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> *cg_a,
                                          CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> *cg_b) {
  if (cg_ptr[0] == nullptr) {
    rtErr("No neighbor list is yet attached to provide a point of reference for comparisons.",
          "DynamicsIntervention", "setNLWorkspace");
  }
  if (cg_a == nullptr) {
    rtErr("At least one valid neighbor list is required.", "DynamicsIntervention",
          "setNLWorkspace");
  }

  // Check that the attachment configuration matches that of the original neighbor list(s)
  validateNeighborListTypes<Tcoord, Tacc, Tcalc, Tcoord4>("setNLWorkspace");
  
  // Check the neighbor list layout to ensure that the attachments are valid.  The neighbor list
  // workspaces will be attached in the order of their themes: electrostatic (or all) non-bonded
  // interactions, then van-der Waals interactions.  If the neighbor list layout is DUAL, each
  // pair of neighbor lists will be counted as one workspace.  Otherwise, each attached neighbor
  // list in the array of pointers is a contained workspace.
  nl_workspaces.clear();
  switch (nl_layout) {
  case NeighborListKind::MONO:
    if (cg_b != nullptr) {
      rtErr("A second neighbor list is invalid for a " + getEnumerationName(nl_layout) +
            " configuration.", "DynamicsIntervention", "setNLWorkspace");
    }
    if (cg_a->getTheme() != NonbondedTheme::ALL) {
      rtErr("A layout with a single neighbor list requires that the attached neighbor list have "
            "a commensurate theme (" + getEnumerationName(NonbondedTheme::ALL) + ").",
            "DynamicsIntervention", "setNLWorkspace");
    }
    nl_workspaces.push_back(reinterpret_cast<CellGrid<float, int, float, float4>*>(cg_a));;    
    break;
  case NeighborListKind::DUAL:
    if (cg_b == nullptr) {
      rtErr("A second neighbor list is required for a " + getEnumerationName(nl_layout) +
            " configuration.", "DynamicsIntervention", "setNLWorkspace");
    }
    nl_workspaces.resize(2);
    bool one_is_all = false;
    switch (cg_a->getTheme()) {
    case NonbondedTheme::ELECTROSTATIC:
      nl_workspaces[0] = reinterpret_cast<CellGrid<float, int, float, float4>*>(cg_a);
      break;
    case NonbondedTheme::VAN_DER_WAALS:
      nl_workspaces[1] = reinterpret_cast<CellGrid<float, int, float, float4>*>(cg_a);
      break;
    case NonbondedTheme::ALL:
      one_is_all = true;
      break;
    }
    switch (cg_b->getTheme()) {
    case NonbondedTheme::ELECTROSTATIC:
      nl_workspaces[0] = reinterpret_cast<CellGrid<float, int, float, float4>*>(cg_b);
      break;
    case NonbondedTheme::VAN_DER_WAALS:
      nl_workspaces[1] = reinterpret_cast<CellGrid<float, int, float, float4>*>(cg_b);
      break;
    case NonbondedTheme::ALL:
      one_is_all = true;
      break;
    }
    if (one_is_all) {
      rtErr("If the neighbor list layout involves " + getEnumerationName(nl_layout) +
            ", grids, themes " + getEnumerationName(cg_a->getTheme()) + " and " +
            getEnumerationName(cg_b->getTheme()) + " are invalid.", "DynamicsIntervention",
            "setNLWorkspace");
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void DynamicsIntervention::validateAnalysisIndex(const int index, const std::vector<T*> &assays,
                                                 const char* caller) {
  if (index < 0 || index >= assays.size()) {
    rtErr("Index " + std::to_string(index) + " is invalid for a set of " +
          std::to_string(assays.size()) + " analyses.", "DynamicsIntervention", caller);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void
DynamicsIntervention::validateNLAttachment(const NonbondedTheme theme,
                                           const std::vector<CellGrid<float,
                                                                      int,
                                                                      float,
                                                                      float4>*> &attachments,
                                           const char* caller) const {

  // Check the neighbor list typing
  validateNeighborListTypes<Tcoord, Tacc, Tcalc, Tcoord4>(caller);
  
  // Check the neighbor list layout and requested theme for consistency
  switch (nl_layout) {
  case NeighborListKind::DUAL:
    {
      if (attachments.size() < 2) {
        const std::string bridge = (attachments.size() == 0) ? std::string(" were ") :
                                                                std::string(" was ");
        rtErr("Two neighbor lists with different themes are expected to be attached.  Instead, " +
              std::to_string(attachments.size()) + bridge + " found.", "DynamicsIntervention",
              caller);
      }
      switch (theme) {
      case NonbondedTheme::ELECTROSTATIC:
        if (attachments[0]->getTheme() != theme) {
          rtErr("The first attached neighbor list workspace is expected to have theme " +
                getEnumerationName(theme) + ".", "DynamicsIntervention", caller);
        }
        break;
      case NonbondedTheme::VAN_DER_WAALS:
        if (attachments[0]->getTheme() != theme) {
          rtErr("The second attached neighbor list workspace is expected to have theme " +
                getEnumerationName(theme) + ".", "DynamicsIntervention", caller);
        }
        break;
      case NonbondedTheme::ALL:
        rtErr("A neighbor list handling " + getEnumerationName(theme) + "particles is "
              "incompatible with a " + getEnumerationName(nl_layout) + " layout.",
              "DynamicsIntervention", caller);
      }
    }
    break;
  case NeighborListKind::MONO:
    if (attachments.size() == 0) {
      rtErr("No neighbor list has yet been attached.", "DynamicsIntervention", caller);
    }
    if (theme != NonbondedTheme::ALL) {
      rtErr("A non-bonded theme of " + getEnumerationName(theme) + " cannot be requested of a "
            "neighbor list workspace for " + getEnumerationName(nl_layout) + " layout.",
            "DynamicsIntervention", caller);
    }
    if (attachments[0]->getTheme() != NonbondedTheme::ALL) {
      rtErr("The attached neighbor list workspace for " + getEnumerationName(nl_layout) +
            " layout has theme " + getEnumerationName(attachments[0]->getTheme()) + ".",
            "DynamicsIntervention", caller);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>*
DynamicsIntervention::getNLAttachment(const NonbondedTheme theme,
                                      std::vector<CellGrid<float, int,
                                                           float, float4>*> *attachments,
                                      const char* caller) {
  validateNLAttachment<Tcoord, Tacc, Tcalc, Tcoord4>(theme, *attachments, caller);
  switch (nl_layout) {
  case NeighborListKind::DUAL:
    switch (theme) {
    case NonbondedTheme::ELECTROSTATIC:
      return reinterpret_cast<CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>*>(attachments->at(0));
    case NonbondedTheme::VAN_DER_WAALS:
      return reinterpret_cast<CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>*>(attachments->at(1));
    case NonbondedTheme::ALL:
      break;
    }
    break;
  case NeighborListKind::MONO:
    return reinterpret_cast<CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>*>(attachments->at(0));
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
const CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>*
DynamicsIntervention::getNLAttachment(const NonbondedTheme theme,
                                      const std::vector<CellGrid<float, int,
                                                                 float, float4>*> &attachments,
                                      const char* caller) const {
  validateNLAttachment<Tcoord, Tacc, Tcalc, Tcoord4>(theme, attachments, caller);
  switch (nl_layout) {
  case NeighborListKind::DUAL:
    switch (theme) {
    case NonbondedTheme::ELECTROSTATIC:
      return reinterpret_cast<CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>*>(attachments[0]);
    case NonbondedTheme::VAN_DER_WAALS:
      return reinterpret_cast<CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>*>(attachments[1]);
    case NonbondedTheme::ALL:
      break;
    }
    break;
  case NeighborListKind::MONO:
    return reinterpret_cast<CellGrid<Tcoord, Tacc, Tcalc, Tcoord4>*>(attachments[0]);
  }
  __builtin_unreachable();
}

} // namespace mm
} // namespace stormm
