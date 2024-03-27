#include "copyright.h"
#include "Topology/atomgraph_abstracts.h"
#include "Parsing/polynumeric.h"
#include "Parsing/parse.h"
#include "clash_detection.h"

namespace stormm {
namespace structure {

using parse::NumberFormat;
using parse::char4ToString;
using parse::realToString;
using topology::ChemicalDetailsKit;

//-------------------------------------------------------------------------------------------------
ClashReport::ClashReport(const double clash_distance_in, const double clash_ratio_in,
                         const AtomGraph *ag_pointer_in) :
    clash_count{0}, clash_distance{clash_distance_in}, clash_ratio{clash_ratio_in},
    pairs{}, distances{}, kinds{}, ag_pointer{const_cast<AtomGraph*>(ag_pointer_in)}
{}

//-------------------------------------------------------------------------------------------------
ClashReport::ClashReport(const AtomGraph *ag_pointer_in) :
    ClashReport(default_minimize_clash_r0, default_minimize_clash_ratio, ag_pointer_in)
{}

//-------------------------------------------------------------------------------------------------
int ClashReport::getClashCount() const {
  return clash_count;
}

//-------------------------------------------------------------------------------------------------
double ClashReport::getMinimumDistance() const {
  return clash_distance;
}

//-------------------------------------------------------------------------------------------------
double ClashReport::getMinimumSigmaRatio() const {
  return clash_ratio;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* ClashReport::getTopologyPointer() const {
  return ag_pointer;
}

//-------------------------------------------------------------------------------------------------
int2 ClashReport::getClashingPair(const int index) const {
  validateClashIndex(index);
  return pairs[index];
}

//-------------------------------------------------------------------------------------------------
double ClashReport::getClashDistance(const int index) const {
  validateClashIndex(index);
  return distances[index];
}

//-------------------------------------------------------------------------------------------------
ClashKind ClashReport::getClashKind(const int index) const {
  validateClashIndex(index);
  return kinds[index];
}

//-------------------------------------------------------------------------------------------------
std::string ClashReport::getClashDescription(const int clash_index) const {
  validateClashIndex(clash_index);
  if (ag_pointer == nullptr) {
    rtErr("No topology is referenced.", "ClashReport", "getClashDescription");
  }
  const ChemicalDetailsKit cdk = ag_pointer->getChemicalDetailsKit();
  const int atom_i = pairs[clash_index].x;
  const int atom_j = pairs[clash_index].y;
  const int res_i = ag_pointer->getResidueIndex(atom_i);
  const int res_j = ag_pointer->getResidueIndex(atom_j);
  const char4 name_i = cdk.atom_names[atom_i];
  const char4 name_j = cdk.atom_names[atom_j];
  const char4 rname_i = cdk.res_names[atom_i];
  const char4 rname_j = cdk.res_names[atom_j];
  const std::string result = "Atom " + std::to_string(cdk.atom_numbers[atom_i]) + " (" +
                             char4ToString(name_i) + " in " + char4ToString(rname_i) + " " +
                             std::to_string(cdk.res_numbers[atom_i]) + ") clashes with atom " +
                             std::to_string(cdk.atom_numbers[atom_j]) + " (" +
                             char4ToString(name_j) + " in " + char4ToString(rname_j) + " " +
                             std::to_string(cdk.res_numbers[atom_j]) + ") at a distance of " +
                             realToString(distances[clash_index], 7, 4,
                                          NumberFormat::STANDARD_REAL) + ".";
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string ClashReport::getClashDescription(const int clash_index,
                                             const CoordinateFrameReader &crd_object) const {
  const int atom_i = pairs[clash_index].x;
  const int atom_j = pairs[clash_index].y;
  const std::string result = getClashDescription(clash_index) +
                             atomPairCoordinates(crd_object.xcrd[atom_i], crd_object.ycrd[atom_i],
                                                 crd_object.zcrd[atom_i], crd_object.xcrd[atom_j],
                                                 crd_object.ycrd[atom_j], crd_object.zcrd[atom_j]);
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string ClashReport::getClashDescription(const int clash_index,
                                             const CoordinateFrame &crd_object) const {
  return getClashDescription(clash_index, crd_object.data());
}

//-------------------------------------------------------------------------------------------------
std::string ClashReport::getClashDescription(const int clash_index,
                                             const CoordinateFrame *crd_object) const {
  return getClashDescription(clash_index, crd_object->data());
}

//-------------------------------------------------------------------------------------------------
std::string ClashReport::getClashDescription(const int clash_index,
                                             const PhaseSpaceReader &crd_object) const {
  const int atom_i = pairs[clash_index].x;
  const int atom_j = pairs[clash_index].y;
  const std::string result = getClashDescription(clash_index) +
                             atomPairCoordinates(crd_object.xcrd[atom_i], crd_object.ycrd[atom_i],
                                                 crd_object.zcrd[atom_i], crd_object.xcrd[atom_j],
                                                 crd_object.ycrd[atom_j], crd_object.zcrd[atom_j]);
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string ClashReport::getClashDescription(const int clash_index,
                                             const PhaseSpace &crd_object) const {
  return getClashDescription(clash_index, crd_object.data());
}

//-------------------------------------------------------------------------------------------------
std::string ClashReport::getClashDescription(const int clash_index,
                                             const PhaseSpace *crd_object) const {
  return getClashDescription(clash_index, crd_object->data());
}

//-------------------------------------------------------------------------------------------------
std::string ClashReport::getClashDescription(const int clash_index,
                                             const PsSynthesisReader &crd_object,
                                             const int system_index) const {
  validateSystemIndex(system_index, crd_object.system_count);
  const int system_start = crd_object.atom_starts[system_index];
  const int atom_i = pairs[clash_index].x + system_start;
  const int atom_j = pairs[clash_index].y + system_start;
  const double x_i = hostInt95ToDouble(crd_object.xcrd[atom_i],
                                       crd_object.xcrd_ovrf[atom_i]) * crd_object.inv_gpos_scale;
  const double y_i = hostInt95ToDouble(crd_object.ycrd[atom_i],
                                       crd_object.ycrd_ovrf[atom_i]) * crd_object.inv_gpos_scale;
  const double z_i = hostInt95ToDouble(crd_object.zcrd[atom_i],
                                       crd_object.zcrd_ovrf[atom_i]) * crd_object.inv_gpos_scale;
  const double x_j = hostInt95ToDouble(crd_object.xcrd[atom_j],
                                       crd_object.xcrd_ovrf[atom_j]) * crd_object.inv_gpos_scale;
  const double y_j = hostInt95ToDouble(crd_object.ycrd[atom_j],
                                       crd_object.ycrd_ovrf[atom_j]) * crd_object.inv_gpos_scale;
  const double z_j = hostInt95ToDouble(crd_object.zcrd[atom_j],
                                       crd_object.zcrd_ovrf[atom_j]) * crd_object.inv_gpos_scale;
  const std::string result = getClashDescription(clash_index) +
                             atomPairCoordinates(x_i, y_i, z_i, x_j, y_j, z_j);
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string ClashReport::getClashDescription(const int clash_index,
                                             const PhaseSpaceSynthesis &crd_object,
                                             const int system_index) const {
  return getClashDescription(clash_index, crd_object.data(), system_index);
}

//-------------------------------------------------------------------------------------------------
std::string ClashReport::getClashDescription(const int clash_index,
                                             const PhaseSpaceSynthesis *crd_object,
                                             const int system_index) const {
  return getClashDescription(clash_index, crd_object->data(), system_index);
}

//-------------------------------------------------------------------------------------------------
std::string ClashReport::getClashDescription(const int clash_index,
                                             const CondensateReader &crd_object,
                                             const int system_index) const {
  validateSystemIndex(system_index, crd_object.system_count);
  const size_t system_start = crd_object.atom_starts[system_index];
  const size_t atom_i = static_cast<size_t>(pairs[clash_index].x) + system_start;
  const size_t atom_j = static_cast<size_t>(pairs[clash_index].y) + system_start;
  double x_i, y_i, z_i, x_j, y_j, z_j;
  switch (crd_object.mode) {
  case PrecisionModel::DOUBLE:
    x_i = crd_object.xcrd[atom_i];
    y_i = crd_object.ycrd[atom_i];
    z_i = crd_object.zcrd[atom_i];
    break;
  case PrecisionModel::SINGLE:
    x_i = crd_object.xcrd_sp[atom_i];
    y_i = crd_object.ycrd_sp[atom_i];
    z_i = crd_object.zcrd_sp[atom_i];
    break;
  }
  const std::string result = getClashDescription(clash_index) +
                             atomPairCoordinates(x_i, y_i, z_i, x_j, y_j, z_j);
  return result;
}

//-------------------------------------------------------------------------------------------------
void ClashReport::addClash(const int atom_i, const int atom_j, const double distance_in) {

  // Validate the atom indices
  if (ag_pointer == nullptr) {
    rtErr("No topology is referenced.", "ClashReport", "addClash");
  }
  const int natom = ag_pointer->getAtomCount();
  if (atom_i >= natom || atom_j >= natom) {
    rtErr("Atom indices " + std::to_string(atom_i) + " and " + std::to_string(atom_j) +
          " are invalid for a topology with " + std::to_string(natom) + " atoms.", "ClashReport",
          "addClash");
  }
  pairs.push_back({ atom_i, atom_j });
  distances.push_back(distance_in);
  ClashKind kind = (distance_in < clash_distance) ? ClashKind::PURE_DISTANCE :
                                                    ClashKind::VAN_DER_WAALS;
  kinds.push_back(kind);
  clash_count += 1;
}

//-------------------------------------------------------------------------------------------------
void ClashReport::clear() {
  clash_count = 0;
  pairs.resize(0);
  distances.resize(0);
  kinds.resize(0);
}

//-------------------------------------------------------------------------------------------------
void ClashReport::setTopologyPointer(const AtomGraph *ag_pointer_in) {
  ag_pointer = const_cast<AtomGraph*>(ag_pointer_in);
}

//-------------------------------------------------------------------------------------------------
void ClashReport::setTopologyPointer(const AtomGraph &ag_in) {
  setTopologyPointer(ag_in.getSelfPointer());
}

//-------------------------------------------------------------------------------------------------
void ClashReport::setMinimumDistance(const double clash_distance_in) {
  clash_distance = clash_distance_in;
}

//-------------------------------------------------------------------------------------------------
void ClashReport::setMinimumSigmaRatio(const double clash_ratio_in) {
  clash_ratio = clash_ratio_in;
}

//-------------------------------------------------------------------------------------------------
void ClashReport::validateClashIndex(const int index) const {
  if (index < 0 || index >= clash_count) {
    rtErr("Index " + std::to_string(index) + " is invalid for a clash report with " +
          std::to_string(clash_count) + " clashes.", "ClashReport", "getClashingPair");
  }
}

//-------------------------------------------------------------------------------------------------
void ClashReport::validateSystemIndex(const int system_index, const int system_count) const {
  if (system_index >= system_count) {
    rtErr("System index " + std::to_string(system_index) + " is invalid for a synthesis of " +
          std::to_string(system_count) + " systems.", "ClashReport", "getClashDescription");
  }
}

//-------------------------------------------------------------------------------------------------
std::string ClashReport::atomPairCoordinates(const double x_i, const double y_i, const double z_i,
                                             const double x_j, const double y_j,
                                             const double z_j) const {
  std::string result = "\n      XYZ [ " + realToString(x_i, 9, 4, NumberFormat::STANDARD_REAL) +
                       " " + realToString(y_i, 9, 4, NumberFormat::STANDARD_REAL) + " " +
                       realToString(z_i, 9, 4, NumberFormat::STANDARD_REAL) + "] - [ " +
                       realToString(x_j, 9, 4, NumberFormat::STANDARD_REAL) + " " +
                       realToString(y_j, 9, 4, NumberFormat::STANDARD_REAL) + " " +
                       realToString(z_j, 9, 4, NumberFormat::STANDARD_REAL) + "];\n";
  return result;
}

//-------------------------------------------------------------------------------------------------
bool detectClash(const CoordinateFrameReader &cfr, const ValenceKit<double> &vk,
                 const NonbondedKit<double> &nbk, const StaticExclusionMask *mask,
                 const double elec_limit, const double vdw_ratio, ClashReport *summary) {
  return detectClash<double, double>(cfr.xcrd, cfr.ycrd, cfr.zcrd, vk, nbk, mask, elec_limit,
                                     vdw_ratio, 1.0, summary);
}

//-------------------------------------------------------------------------------------------------
bool detectClash(const CoordinateFrame &cf, const AtomGraph &ag, const StaticExclusionMask &mask,
                 const double elec_limit, const double vdw_ratio, ClashReport *summary) {
  return detectClash(cf.data(), ag.getDoublePrecisionValenceKit(),
                     ag.getDoublePrecisionNonbondedKit(), mask.getSelfPointer(), elec_limit,
                     vdw_ratio, summary);
}

//-------------------------------------------------------------------------------------------------
bool detectClash(const CoordinateFrame *cf, const AtomGraph *ag, const StaticExclusionMask *mask,
                 const double elec_limit, const double vdw_ratio, ClashReport *summary) {
  return detectClash(cf->data(), ag->getDoublePrecisionValenceKit(),
                     ag->getDoublePrecisionNonbondedKit(), mask, elec_limit, vdw_ratio, summary);
}

//-------------------------------------------------------------------------------------------------
bool detectClash(const CoordinateFrame &cf, const AtomGraph &ag, const StaticExclusionMask &mask,
                 ClashReport *summary) {
  return detectClash(cf.data(), ag.getDoublePrecisionValenceKit(),
                     ag.getDoublePrecisionNonbondedKit(), mask.getSelfPointer(),
                     summary->getMinimumDistance(), summary->getMinimumSigmaRatio(), summary);
}

//-------------------------------------------------------------------------------------------------
bool detectClash(const CoordinateFrame *cf, const AtomGraph *ag, const StaticExclusionMask *mask,
                 ClashReport *summary) {
  return detectClash(cf->data(), ag->getDoublePrecisionValenceKit(),
                     ag->getDoublePrecisionNonbondedKit(), mask, summary->getMinimumDistance(),
                     summary->getMinimumSigmaRatio(), summary);
}

//-------------------------------------------------------------------------------------------------
bool detectClash(const PhaseSpaceReader &psr, const ValenceKit<double> &vk,
                 const NonbondedKit<double> &nbk, const StaticExclusionMask *mask,
                 const double elec_limit, const double vdw_ratio, ClashReport *summary) {
  return detectClash<double, double>(psr.xcrd, psr.ycrd, psr.zcrd, vk, nbk, mask, elec_limit,
                                     vdw_ratio, 1.0, summary);
}

//-------------------------------------------------------------------------------------------------
bool detectClash(const PhaseSpace &ps, const AtomGraph &ag, const StaticExclusionMask &mask,
                 const double elec_limit, const double vdw_ratio, ClashReport *summary) {
  return detectClash(ps.data(), ag.getDoublePrecisionValenceKit(),
                     ag.getDoublePrecisionNonbondedKit(), mask.getSelfPointer(), elec_limit,
                     vdw_ratio, summary);
}

//-------------------------------------------------------------------------------------------------
bool detectClash(const PhaseSpace *ps, const AtomGraph *ag, const StaticExclusionMask *mask,
                 const double elec_limit, const double vdw_ratio, ClashReport *summary) {
  return detectClash(ps->data(), ag->getDoublePrecisionValenceKit(),
                     ag->getDoublePrecisionNonbondedKit(), mask, elec_limit, vdw_ratio, summary);
}

//-------------------------------------------------------------------------------------------------
bool detectClash(const PhaseSpace &ps, const AtomGraph &ag, const StaticExclusionMask &mask,
                 ClashReport *summary) {
  return detectClash(ps.data(), ag.getDoublePrecisionValenceKit(),
                     ag.getDoublePrecisionNonbondedKit(), mask.getSelfPointer(),
                     summary->getMinimumDistance(), summary->getMinimumSigmaRatio(), summary);
}

//-------------------------------------------------------------------------------------------------
bool detectClash(const PhaseSpace *ps, const AtomGraph *ag, const StaticExclusionMask *mask,
                 ClashReport *summary) {
  return detectClash(ps->data(), ag->getDoublePrecisionValenceKit(),
                     ag->getDoublePrecisionNonbondedKit(), mask,
                     summary->getMinimumDistance(), summary->getMinimumSigmaRatio(), summary);
}

} // namespace structure
} // namespace stormm
