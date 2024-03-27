#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

using constants::CaseSensitivity;
using parse::strcmpCased;


//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ImagingMethod input) {
  switch (input) {
  case ImagingMethod::PRIMARY_UNIT_CELL:
    return std::string("PRIMARY_UNIT_CELL");
  case ImagingMethod::MINIMUM_IMAGE:
    return std::string("MINIMUM_IMAGE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RMSDMethod input) {
  switch (input) {
  case RMSDMethod::ALIGN_MASS:
    return std::string("ALIGN_MASS");
  case RMSDMethod::ALIGN_GEOM:
    return std::string("ALIGN_GEOM");
  case RMSDMethod::NO_ALIGN_MASS:
    return std::string("NO_ALIGN_MASS");
  case RMSDMethod::NO_ALIGN_GEOM:
    return std::string("NO_ALIGN_GEOM");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RMSDAlignmentProtocol input) {
  switch (input) {
  case RMSDAlignmentProtocol::BUILD_CORE:
    return std::string("BUILD_CORE");
  case RMSDAlignmentProtocol::ALIGN_CORE:
    return std::string("ALIGN_CORE");
  case RMSDAlignmentProtocol::ALIGN_ALL:
    return std::string("ALIGN_ALL");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RMSDTask input) {
  switch (input) {
  case RMSDTask::REFERENCE:
    return std::string("REFERENCE");
  case RMSDTask::MATRIX:
    return std::string("MATRIX");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const VirtualSiteActivity input) {
  switch (input) {
  case VirtualSiteActivity::PLACEMENT:
    return std::string("PLACEMENT");
  case VirtualSiteActivity::TRANSMIT_FORCES:
    return std::string("TRANSMIT_FORCES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const GridDetail input) {
  switch (input) {
  case GridDetail::OCCLUSION:
    return std::string("OCCLUSION");
  case GridDetail::OCCLUSION_FIELD:
    return std::string("OCCLUSION_FIELD");
  case GridDetail::NONBONDED_FIELD:
    return std::string("NONBONDED_FIELD");
  case GridDetail::NONBONDED_ATOMIC:
    return std::string("NONBONDED_ATOMIC");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MappingActivity input) {
  switch (input) {
  case MappingActivity::PARTICLE_TO_MESH:
    return std::string("PARTICLE_TO_MESH");
  case MappingActivity::MESH_TO_PARTICLE:
    return std::string("MESH_TO_PARTICLE");
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ClashKind input) {
  switch (input) {
  case ClashKind::VAN_DER_WAALS:
    return std::string("VAN_DER_WAALS");
  case ClashKind::PURE_DISTANCE:
    return std::string("PURE_DISTANCE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const OffMeshProtocol input) {
  switch (input) {
  case OffMeshProtocol::DIE:
    return std::string("DIE");
  case OffMeshProtocol::EXTRAPOLATE:
    return std::string("EXTRAPOLATE");
  case OffMeshProtocol::ZERO:
    return std::string("ZERO");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const SamplingIntensity input) {
  switch (input) {
  case SamplingIntensity::MINIMAL:
    return std::string("MINIMAL");
  case SamplingIntensity::LIGHT:
    return std::string("LIGHT");
  case SamplingIntensity::HEAVY:
    return std::string("HEAVY");
  case SamplingIntensity::EXHAUSTIVE:
    return std::string("EXHAUSTIVE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const BoundaryCondition input) {
  switch (input) {
  case BoundaryCondition::ISOLATED:
    return std::string("ISOLATED");
  case BoundaryCondition::PERIODIC:
    return std::string("PERIODIC");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MeshPosition input) {
  switch (input) {
  case MeshPosition::MOLECULE:
    return std::string("MOLECULE");
  case MeshPosition::ORIGIN:
    return std::string("ORIGIN");
  case MeshPosition::ARBITRARY:
    return std::string("ARBITRARY");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RattleMethod input) {
  switch (input) {
  case RattleMethod::SEQUENTIAL:
    return std::string("SEQUENTIAL");
  case RattleMethod::CENTER_SUM:
    return std::string("CENTER_SUM");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ApplyConstraints input) {
  switch (input) {
  case ApplyConstraints::YES:
    return std::string("YES");
  case ApplyConstraints::NO:
    return std::string("NO");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
RMSDMethod translateRMSDMethod(const std::string &input) {
  if (strcmpCased(input, std::string("align_mass"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("align_mass_weighted"), CaseSensitivity::NO)) {
    return RMSDMethod::ALIGN_MASS;
  }
  else if (strcmpCased(input, std::string("align_geom"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("align_geometry"), CaseSensitivity::NO)) {
    return RMSDMethod::ALIGN_GEOM;
  }
  else if (strcmpCased(input, std::string("noalign_mass"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("noalign_mass_weighted"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("no_align_mass"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("no_align_mass_weighted"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("nonalign_mass"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("nonalign_mass_weighted"), CaseSensitivity::NO)) {
    return RMSDMethod::NO_ALIGN_MASS;
  }
  else if (strcmpCased(input, std::string("noalign_geom"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("noalign_geometry"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("no_align_geom"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("no_align_geometry"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("nonalign_geom"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("nonalign_geometry"), CaseSensitivity::NO)) {
    return RMSDMethod::NO_ALIGN_GEOM;
  }
  else {
    rtErr("\"" + input + "\" is not a valid RMSD calculation method.", "translateRMSDMethod");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
GridDetail translateGridDetail(const std::string &input) {
  if (strcmpCased(input, std::string("occlusion"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("spacefill"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("space_fill"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("spacefilling"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("space_filling"), CaseSensitivity::NO)) {
    return GridDetail::OCCLUSION;
  }
  else if (strcmpCased(input, std::string("nonbonded"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("non_bonded"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("nonbonded_field"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("non_bonded_field"), CaseSensitivity::NO)) {
    return GridDetail::NONBONDED_FIELD;
  }
  else if (strcmpCased(input, std::string("nonbonded_atomic"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("non_bonded_atomic"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("atomic"), CaseSensitivity::NO)) {
    return GridDetail::NONBONDED_ATOMIC;
  }
  else {
    rtErr("\"" + input + "\" is not a valid level of grid detail.", "translateGridDetail");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
SamplingIntensity translateSamplingIntensity(const std::string &input) {
  if (strcmpCased(input, std::string("minimal"), CaseSensitivity::NO)) {
    return SamplingIntensity::MINIMAL;
  }
  else if (strcmpCased(input, std::string("light"), CaseSensitivity::NO)) {
    return SamplingIntensity::LIGHT;
  }
  else if (strcmpCased(input, std::string("heavy"), CaseSensitivity::NO)) {
    return SamplingIntensity::HEAVY;
  }
  else if (strcmpCased(input, std::string("exhaustive"), CaseSensitivity::NO)) {
    return SamplingIntensity::EXHAUSTIVE;
  }
  else {
    rtErr("\"" + input + "\" is not a valid sampling intensity.", "translateSamplingIntensity");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
BoundaryCondition translateBoundaryCondition(const std::string &input) {
  if (strcmpCased(input, std::string("isolated"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("nonperiodic"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("non_periodic"), CaseSensitivity::NO)) {
    return BoundaryCondition::ISOLATED;
  }
  else if (strcmpCased(input, std::string("periodic"), CaseSensitivity::NO)) {
    return BoundaryCondition::PERIODIC;
  }
  else {
    rtErr("\"" + input + "\" is not a valid boundary condition.", "translateBoundaryCondition");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
MeshPosition translateMeshPosition(const std::string &input) {
  if (strcmpCased(input, std::string("molecule"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("structure"), CaseSensitivity::NO)) {
    return MeshPosition::MOLECULE;
  }
  else if (strcmpCased(input, std::string("origin"), CaseSensitivity::NO)) {
    return MeshPosition::ORIGIN;
  }
  else if (strcmpCased(input, std::string("arbitrary"), CaseSensitivity::NO)) {
    return MeshPosition::ARBITRARY;
  }
  else {
    rtErr("\"" + input + "\" is not a valid mesh alignment.", "translateMeshPosition");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
RattleMethod translateRattleMethod(const std::string &input) {
  if (strcmpCased(input, std::string("sequential"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("traditional"), CaseSensitivity::NO)) {
    return RattleMethod::SEQUENTIAL;
  }
  else if (strcmpCased(input, std::string("center_sum"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("bias_free"), CaseSensitivity::NO)) {
    return RattleMethod::CENTER_SUM;
  }
  else {
    rtErr("\"" + input + "\" is not a valid RATTLE approach.", "translateRattleMethod");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
ApplyConstraints translateApplyConstraints(const std::string &input) {
  if (strcmpCased(input, std::string("off"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("false"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("no"), CaseSensitivity::NO)) {
    return ApplyConstraints::NO;
  }
  else if (strcmpCased(input, std::string("on"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("true"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("yes"), CaseSensitivity::NO)) {
    return ApplyConstraints::YES;
  }
  else {
    rtErr("\"" + input + "\" is not a valid constraint instruction.", "translateApplyConstraints");
  }
  __builtin_unreachable();
}

} // namespace structure
} // namespace stormm
