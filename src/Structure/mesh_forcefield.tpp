// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename T> MeshFFKit<T>::MeshFFKit() :
    ljrule{VdwCombiningRule::GEOMETRIC}, coulomb{amber_ancient_bioq},
    coulomb_f{static_cast<float>(amber_ancient_bioq)}, clash_ratio{0.0},
    clash_ratio_f{static_cast<float>(0.0)}, clash_distance{0.0},
    clash_distance_f{static_cast<float>(0.0)}, probe_lja{nullptr}, probe_ljb{nullptr},
    probe_ljsig{nullptr}, softcore_lja{nullptr}, softcore_ljb{nullptr}, softcore_ljc{nullptr},
    softcore_ljd{nullptr}, softcore_lje{nullptr}, softcore_ljf{nullptr}, softcore_qq{nullptr}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshFFKit<T>::MeshFFKit(const VdwCombiningRule ljrule_in, const double coulomb_in,
                        const double clash_ratio_in, const double clash_distance_in,
                        const T* probe_lja_in, const T* probe_ljb_in, const T* probe_ljsig_in,
                        const T* softcore_lja_in, const T* softcore_ljb_in,
                        const T* softcore_ljc_in, const T* softcore_ljd_in,
                        const T* softcore_lje_in, const T* softcore_ljf_in,
                        const T* softcore_qq_in) :
    ljrule{ljrule_in}, coulomb{coulomb_in}, coulomb_f{static_cast<float>(coulomb_in)},
    clash_ratio{clash_ratio_in}, clash_ratio_f{static_cast<float>(clash_ratio_in)},
    clash_distance{clash_distance_in}, clash_distance_f{static_cast<float>(clash_distance_in)},
    probe_lja{probe_lja_in}, probe_ljb{probe_ljb_in}, probe_ljsig{probe_ljsig_in},
    softcore_lja{softcore_lja_in}, softcore_ljb{softcore_ljb_in}, softcore_ljc{softcore_ljc_in},
    softcore_ljd{softcore_ljd_in}, softcore_lje{softcore_lje_in}, softcore_ljf{softcore_ljf_in},
    softcore_qq{softcore_qq_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshForceField<T>::MeshForceField(const VdwCombiningRule mixing_protocol_in,
                                  const double coulomb_constant_in, const double clash_ratio_in,
                                  const double clash_distance_in, const int lj_type_count_in) :
    lj_type_count{0},
    mixing_protocol{mixing_protocol_in}, coulomb_constant{coulomb_constant_in},
    clash_ratio{clash_ratio_in}, clash_distance{clash_distance_in},
    probe_lja{HybridKind::POINTER, "mesh_probe_lja"},
    probe_ljb{HybridKind::POINTER, "mesh_probe_ljb"},
    probe_lj_sigma{HybridKind::POINTER, "mesh_probe_ljsigma"},
    probe_softcore_lja{HybridKind::POINTER, "mesh_softcore_lja"},
    probe_softcore_ljb{HybridKind::POINTER, "mesh_softcore_ljb"},
    probe_softcore_ljc{HybridKind::POINTER, "mesh_softcore_ljc"},
    probe_softcore_ljd{HybridKind::POINTER, "mesh_softcore_ljd"},
    probe_softcore_lje{HybridKind::POINTER, "mesh_softcore_lje"},
    probe_softcore_ljf{HybridKind::POINTER, "mesh_softcore_ljf"},
    elec_softcore{HybridKind::POINTER, "mesh_softcore_qq"},
    softcore_data{HybridKind::ARRAY, "mesh_softcore_data"},
    ref_probe_lja{HybridKind::POINTER, "mesh_rprobe_lja"},
    ref_probe_ljb{HybridKind::POINTER, "mesh_rprobe_ljb"},
    ref_probe_lj_sigma{HybridKind::POINTER, "mesh_rprobe_ljsigma"},
    ref_probe_softcore_lja{HybridKind::POINTER, "mesh_rsc_lja"},
    ref_probe_softcore_ljb{HybridKind::POINTER, "mesh_rsc_ljb"},
    ref_probe_softcore_ljc{HybridKind::POINTER, "mesh_rsc_ljc"},
    ref_probe_softcore_ljd{HybridKind::POINTER, "mesh_rsc_ljd"},
    ref_probe_softcore_lje{HybridKind::POINTER, "mesh_rsc_lje"},
    ref_probe_softcore_ljf{HybridKind::POINTER, "mesh_rsc_ljf"},
    ref_elec_softcore{HybridKind::POINTER, "mesh_rsc_qq"},
    ref_softcore_data{HybridKind::ARRAY, "mesh_rsc_data"}
{
  allocate(lj_type_count_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshForceField<T>::MeshForceField(const VdwCombiningRule mixing_protocol_in,
                                  const double clash_ratio_in, const double clash_distance_in,
                                  const AtomGraph *ag) :
    MeshForceField(mixing_protocol_in, ag->getCoulombConstant(), clash_ratio_in, clash_distance_in,
                   ag->getLJTypeCount())
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshForceField<T>::MeshForceField(const VdwCombiningRule mixing_protocol_in,
                                  const double clash_ratio_in, const double clash_distance_in,
                                  const AtomGraph &ag) :
    MeshForceField(mixing_protocol_in, ag.getCoulombConstant(), clash_ratio_in, clash_distance_in,
                   ag.getLJTypeCount())
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshForceField<T>::MeshForceField(const VdwCombiningRule mixing_protocol_in,
                                  const double coulomb_constant_in, const double clash_ratio_in,
                                  const double clash_distance_in,
                                  const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon) :
    MeshForceField(mixing_protocol_in, coulomb_constant_in, clash_ratio_in, clash_distance_in,
                   probe_sigma.size())
{
  if (probe_sigma.size() != probe_epsilon.size()) {
    rtErr("Consistent inputs must be provided for Lennard-Jones sigma (" +
          std::to_string(probe_sigma.size()) + " parameters) and epsilon (" +
          std::to_string(probe_epsilon.size()) + " parameters).", "MeshForceField");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> MeshForceField<T>::MeshForceField(const MeshForceField<T> &original) :
    lj_type_count{original.lj_type_count},
    mixing_protocol{original.mixing_protocol},
    coulomb_constant{original.coulomb_constant},
    clash_ratio{original.clash_ratio},
    clash_distance{original.clash_distance},
    probe_lja{original.probe_lja},
    probe_ljb{original.probe_ljb},
    probe_lj_sigma{original.probe_lj_sigma},
    probe_softcore_lja{original.probe_softcore_lja},
    probe_softcore_ljb{original.probe_softcore_ljb},
    probe_softcore_ljc{original.probe_softcore_ljc},
    probe_softcore_ljd{original.probe_softcore_ljd},
    probe_softcore_lje{original.probe_softcore_lje},
    probe_softcore_ljf{original.probe_softcore_ljf},
    elec_softcore{original.elec_softcore},
    softcore_data{original.softcore_data},
    ref_probe_lja{original.ref_probe_lja},
    ref_probe_ljb{original.ref_probe_ljb},
    ref_probe_lj_sigma{original.ref_probe_lj_sigma},
    ref_probe_softcore_lja{original.ref_probe_softcore_lja},
    ref_probe_softcore_ljb{original.ref_probe_softcore_ljb},
    ref_probe_softcore_ljc{original.ref_probe_softcore_ljc},
    ref_probe_softcore_ljd{original.ref_probe_softcore_ljd},
    ref_probe_softcore_lje{original.ref_probe_softcore_lje},
    ref_probe_softcore_ljf{original.ref_probe_softcore_ljf},
    ref_elec_softcore{original.ref_elec_softcore},
    ref_softcore_data{original.ref_softcore_data}
{
  // Repair pointers
  rebasePointers();
}

//-------------------------------------------------------------------------------------------------
template <typename T> MeshForceField<T>::MeshForceField(MeshForceField<T> &&original) :
    lj_type_count{original.lj_type_count},
    mixing_protocol{original.mixing_protocol},
    coulomb_constant{original.coulomb_constant},
    clash_ratio{original.clash_ratio},
    clash_distance{original.clash_distance},
    probe_lja{std::move(original.probe_lja)},
    probe_ljb{std::move(original.probe_ljb)},
    probe_lj_sigma{std::move(original.probe_lj_sigma)},
    probe_softcore_lja{std::move(original.probe_softcore_lja)},
    probe_softcore_ljb{std::move(original.probe_softcore_ljb)},
    probe_softcore_ljc{std::move(original.probe_softcore_ljc)},
    probe_softcore_ljd{std::move(original.probe_softcore_ljd)},
    probe_softcore_lje{std::move(original.probe_softcore_lje)},
    probe_softcore_ljf{std::move(original.probe_softcore_ljf)},
    elec_softcore{std::move(original.elec_softcore)},
    softcore_data{std::move(original.softcore_data)},
    ref_probe_lja{std::move(original.ref_probe_lja)},
    ref_probe_ljb{std::move(original.ref_probe_ljb)},
    ref_probe_lj_sigma{std::move(original.ref_probe_lj_sigma)},
    ref_probe_softcore_lja{std::move(original.ref_probe_softcore_lja)},
    ref_probe_softcore_ljb{std::move(original.ref_probe_softcore_ljb)},
    ref_probe_softcore_ljc{std::move(original.ref_probe_softcore_ljc)},
    ref_probe_softcore_ljd{std::move(original.ref_probe_softcore_ljd)},
    ref_probe_softcore_lje{std::move(original.ref_probe_softcore_lje)},
    ref_probe_softcore_ljf{std::move(original.ref_probe_softcore_ljf)},
    ref_elec_softcore{std::move(original.ref_elec_softcore)},
    ref_softcore_data{std::move(original.ref_softcore_data)}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshForceField<T>& MeshForceField<T>::operator=(const MeshForceField<T> &other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  lj_type_count = other.lj_type_count;
  mixing_protocol = other.mixing_protocol;
  coulomb_constant = other.coulomb_constant;
  clash_ratio = other.clash_ratio;
  clash_distance = other.clash_distance;
  probe_lja = other.probe_lja;
  probe_ljb = other.probe_ljb;
  probe_lj_sigma = other.probe_lj_sigma;
  probe_softcore_lja = other.probe_softcore_lja;
  probe_softcore_ljb = other.probe_softcore_ljb;
  probe_softcore_ljc = other.probe_softcore_ljc;
  probe_softcore_ljd = other.probe_softcore_ljd;
  probe_softcore_lje = other.probe_softcore_lje;
  probe_softcore_ljf = other.probe_softcore_ljf;
  elec_softcore = other.elec_softcore;
  softcore_data = other.softcore_data;
  ref_probe_lja = other.ref_probe_lja;
  ref_probe_ljb = other.ref_probe_ljb;
  ref_probe_lj_sigma = other.ref_probe_lj_sigma;
  ref_probe_softcore_lja = other.ref_probe_softcore_lja;
  ref_probe_softcore_ljb = other.ref_probe_softcore_ljb;
  ref_probe_softcore_ljc = other.ref_probe_softcore_ljc;
  ref_probe_softcore_ljd = other.ref_probe_softcore_ljd;
  ref_probe_softcore_lje = other.ref_probe_softcore_lje;
  ref_probe_softcore_ljf = other.ref_probe_softcore_ljf;
  ref_elec_softcore = other.ref_elec_softcore;
  ref_softcore_data = other.ref_softcore_data;

  // Repair pointers
  rebasePointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshForceField<T>& MeshForceField<T>::operator=(MeshForceField<T> &&other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  lj_type_count = other.lj_type_count;
  mixing_protocol = other.mixing_protocol;
  coulomb_constant = other.coulomb_constant;
  clash_ratio = other.clash_ratio;
  clash_distance = other.clash_distance;
  probe_lja = std::move(other.probe_lja);
  probe_ljb = std::move(other.probe_ljb);
  probe_lj_sigma = std::move(other.probe_lj_sigma);
  probe_softcore_lja = std::move(other.probe_softcore_lja);
  probe_softcore_ljb = std::move(other.probe_softcore_ljb);
  probe_softcore_ljc = std::move(other.probe_softcore_ljc);
  probe_softcore_ljd = std::move(other.probe_softcore_ljd);
  probe_softcore_lje = std::move(other.probe_softcore_lje);
  probe_softcore_ljf = std::move(other.probe_softcore_ljf);
  elec_softcore = std::move(other.elec_softcore);
  softcore_data = std::move(other.softcore_data);
  ref_probe_lja = std::move(other.ref_probe_lja);
  ref_probe_ljb = std::move(other.ref_probe_ljb);
  ref_probe_lj_sigma = std::move(other.ref_probe_lj_sigma);
  ref_probe_softcore_lja = std::move(other.ref_probe_softcore_lja);
  ref_probe_softcore_ljb = std::move(other.ref_probe_softcore_ljb);
  ref_probe_softcore_ljc = std::move(other.ref_probe_softcore_ljc);
  ref_probe_softcore_ljd = std::move(other.ref_probe_softcore_ljd);
  ref_probe_softcore_lje = std::move(other.ref_probe_softcore_lje);
  ref_probe_softcore_ljf = std::move(other.ref_probe_softcore_ljf);
  ref_elec_softcore = std::move(other.ref_elec_softcore);
  ref_softcore_data = std::move(other.ref_softcore_data);

  // No pointer repair is needed, as everything was moved.
  return *this;
}

//-------------------------------------------------------------------------------------------------
template <typename T> VdwCombiningRule MeshForceField<T>::getCombiningRule() const {
  return mixing_protocol;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double MeshForceField<T>::getClashRatio() const {
  return clash_ratio;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double MeshForceField<T>::getClashDistance() const {
  return clash_distance;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const MeshFFKit<T> MeshForceField<T>::data(const HybridTargetLevel tier) const {
  return MeshFFKit<T>(mixing_protocol, coulomb_constant, clash_ratio, clash_distance,
                      probe_lja.data(tier), probe_ljb.data(tier), probe_lj_sigma.data(tier),
                      probe_softcore_lja.data(tier), probe_softcore_ljb.data(tier),
                      probe_softcore_ljc.data(tier), probe_softcore_ljd.data(tier),
                      probe_softcore_lje.data(tier), probe_softcore_ljf.data(tier),
                      elec_softcore.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const MeshFFKit<double> MeshForceField<T>::referenceData(const HybridTargetLevel tier) const {
  return MeshFFKit<double>(mixing_protocol, coulomb_constant, clash_ratio, clash_distance,
                           ref_probe_lja.data(tier), ref_probe_ljb.data(tier),
                           ref_probe_lj_sigma.data(tier), ref_probe_softcore_lja.data(tier),
                           ref_probe_softcore_ljb.data(tier), ref_probe_softcore_ljc.data(tier),
                           ref_probe_softcore_ljd.data(tier), ref_probe_softcore_lje.data(tier),
                           ref_probe_softcore_ljf.data(tier), ref_elec_softcore.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const MeshFFKit<void> MeshForceField<T>::templateFreeData(const HybridTargetLevel tier) const {
  return MeshFFKit<void>(mixing_protocol, coulomb_constant, clash_ratio, clash_distance,
                         reinterpret_cast<const void*>(probe_lja.data(tier)),
                         reinterpret_cast<const void*>(probe_ljb.data(tier)),
                         reinterpret_cast<const void*>(probe_lj_sigma.data(tier)),
                         reinterpret_cast<const void*>(probe_softcore_lja.data(tier)),
                         reinterpret_cast<const void*>(probe_softcore_ljb.data(tier)),
                         reinterpret_cast<const void*>(probe_softcore_ljc.data(tier)),
                         reinterpret_cast<const void*>(probe_softcore_ljd.data(tier)),
                         reinterpret_cast<const void*>(probe_softcore_lje.data(tier)),
                         reinterpret_cast<const void*>(probe_softcore_ljf.data(tier)),
                         reinterpret_cast<const void*>(elec_softcore.data(tier)));
}

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
//-------------------------------------------------------------------------------------------------
template <typename T>
const MeshFFKit<T> MeshForceField<T>::deviceViewToHostData() const {
  return MeshFFKit<T>(mixing_protocol, coulomb_constant, clash_ratio, clash_distance,
                      probe_lja.getDeviceValidHostPointer(), probe_ljb.getDeviceValidHostPointer(),
                      probe_lj_sigma.getDeviceValidHostPointer(),
                      probe_softcore_lja.getDeviceValidHostPointer(),
                      probe_softcore_ljb.getDeviceValidHostPointer(),
                      probe_softcore_ljc.getDeviceValidHostPointer(),
                      probe_softcore_ljd.getDeviceValidHostPointer(),
                      probe_softcore_lje.getDeviceValidHostPointer(),
                      probe_softcore_ljf.getDeviceValidHostPointer(),
                      elec_softcore.getDeviceValidHostPointer());
}

//-------------------------------------------------------------------------------------------------
template <typename T> const MeshFFKit<double>
MeshForceField<T>::deviceViewToReferenceHostData() const {
  return MeshFFKit<double>(mixing_protocol, coulomb_constant, clash_ratio, clash_distance,
                           ref_probe_lja.getDeviceValidHostPointer(),
                           ref_probe_ljb.getDeviceValidHostPointer(),
                           ref_probe_lj_sigma.getDeviceValidHostPointer(),
                           ref_probe_softcore_lja.getDeviceValidHostPointer(),
                           ref_probe_softcore_ljb.getDeviceValidHostPointer(),
                           ref_probe_softcore_ljc.getDeviceValidHostPointer(),
                           ref_probe_softcore_ljd.getDeviceValidHostPointer(),
                           ref_probe_softcore_lje.getDeviceValidHostPointer(),
                           ref_probe_softcore_ljf.getDeviceValidHostPointer(),
                           ref_elec_softcore.getDeviceValidHostPointer());
}

//-------------------------------------------------------------------------------------------------
template <typename T> const MeshFFKit<void>
MeshForceField<T>::deviceViewToTemplateFreeHostData() const {
  return MeshFFKit<void>(mixing_protocol, coulomb_constant, clash_ratio, clash_distance,
                         reinterpret_cast<void*>(probe_lja.getDeviceValidHostPointer()),
                         reinterpret_cast<void*>(probe_ljb.getDeviceValidHostPointer()),
                         reinterpret_cast<void*>(probe_lj_sigma.getDeviceValidHostPointer()),
                         reinterpret_cast<void*>(probe_softcore_lja.getDeviceValidHostPointer()),
                         reinterpret_cast<void*>(probe_softcore_ljb.getDeviceValidHostPointer()),
                         reinterpret_cast<void*>(probe_softcore_ljc.getDeviceValidHostPointer()),
                         reinterpret_cast<void*>(probe_softcore_ljd.getDeviceValidHostPointer()),
                         reinterpret_cast<void*>(probe_softcore_lje.getDeviceValidHostPointer()),
                         reinterpret_cast<void*>(probe_softcore_ljf.getDeviceValidHostPointer()),
                         reinterpret_cast<void*>(elec_softcore.getDeviceValidHostPointer()));
}
#  endif

//-------------------------------------------------------------------------------------------------
template <typename T> void MeshForceField<T>::upload() {
  softcore_data.upload();
  ref_softcore_data.upload();
}

//-------------------------------------------------------------------------------------------------
template <typename T> void MeshForceField<T>::download() {
  softcore_data.download();
  ref_softcore_data.download();
}
#endif

//-------------------------------------------------------------------------------------------------
template <typename T>
void MeshForceField<T>::setCombiningRule(const VdwCombiningRule mixing_protocol_in) {
  mixing_protocol = mixing_protocol_in;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MeshForceField<T>::setClashDistance(const double clash_distance_in) {
  clash_distance = clash_distance_in;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MeshForceField<T>::setClashRatio(const double clash_ratio_in) {
  clash_ratio = clash_ratio_in;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MeshForceField<T>::setElecSoftcoreParameter(const double coef, const int pos) {
  elec_softcore.putHost(coef, pos);
  ref_elec_softcore.putHost(coef, pos);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MeshForceField<T>::setElecSoftcoreParameters(const Interpolant stencil_kind) {
  double4 abcd_coef = { 0.0, 0.0, 0.0, 0.0 };
  double2 ef_coef = { 0.0, 0.0 };
  if (clash_distance > constants::tiny) {
    const double inv_cd  = 1.0 / clash_distance;
    const double inv_cd2 = inv_cd * inv_cd;
    const double inv_cd6 = inv_cd2 * inv_cd2 * inv_cd2;
    const double qq_switch   =  inv_cd;
    const double dqq_switch  = -inv_cd2;
    const double d2qq_switch =  2.0 * inv_cd2 * inv_cd;
    const double d3qq_switch = -6.0 * inv_cd2 * inv_cd2;
    switch (stencil_kind) {
    case Interpolant::SMOOTHNESS:
      quinticSoftCore(&abcd_coef, &ef_coef, clash_distance, qq_switch, dqq_switch, d2qq_switch,
                      d3qq_switch, 0.0);
      break;
    case Interpolant::FUNCTION_VALUE:
      cubicSoftCore(&abcd_coef, clash_distance, qq_switch, dqq_switch, 0.0);
      break;
    }
  }
  setElecSoftcoreParameter(abcd_coef.x, 0);
  setElecSoftcoreParameter(abcd_coef.y, 1);
  setElecSoftcoreParameter(abcd_coef.z, 2);
  setElecSoftcoreParameter(abcd_coef.w, 3);
  setElecSoftcoreParameter(ef_coef.x, 4);
  setElecSoftcoreParameter(ef_coef.y, 5);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MeshForceField<T>::setLJCoefficients(const int index, const double lja_in,
                                          const double ljb_in) {
  validateTypeIndex(index);
  probe_lja.putHost(lja_in, index);
  probe_ljb.putHost(ljb_in, index);
  ref_probe_lja.putHost(lja_in, index);
  ref_probe_ljb.putHost(ljb_in, index);
  probe_lj_sigma.putHost(sqrt(cbrt(lja_in / ljb_in)), index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MeshForceField<T>::setLJCoefficients(const AtomGraph *ag,
                                          const std::vector<double> &probe_radius,
                                          const std::vector<double> &probe_well_depth) {

  // Set local template-type pointers
  T* probe_lja_ptr = probe_lja.data();
  T* probe_ljb_ptr = probe_ljb.data();
  T* probe_lj_sigma_ptr = probe_lj_sigma.data();

  // Set local double-precision pointers
  double* ref_probe_lja_ptr = ref_probe_lja.data();
  double* ref_probe_ljb_ptr = ref_probe_ljb.data();
  double* ref_probe_lj_sigma_ptr = ref_probe_lj_sigma.data();
  
  // If there is not the right amount of space to store the Lennard-Jones types listed in the
  // topology, this is an error.  Reallocate, but store the electrostatic softcore potential in
  // the process as this would be blown away.
  if (ag->getLJTypeCount() != lj_type_count) {
    reallocate(ag->getLJTypeCount());
  }
  const std::vector<double> ag_sigma = ag->getLennardJonesSigma<double>();
  const std::vector<double> ag_eps = ag->getLennardJonesEpsilon<double>();    
  for (int i = 0; i < lj_type_count; i++) {
    double t_sig6, t_eps;
    switch (mixing_protocol) {
    case VdwCombiningRule::GEOMETRIC:

      // The sigma parameter will be the square root of the product of the probe's given radius
      // with the topology's parameter for the ith atom index.
      t_sig6 = sqrt(ag_sigma[i] * probe_radius[i]);
      t_eps  = sqrt(ag_eps[i] * probe_well_depth[i]);
      break;
    case VdwCombiningRule::LORENTZ_BERTHELOT:

      // The sigma parameter will be the average of the probe's given radius and the topology's
      // parameter for the ith atom index.
      t_sig6 = 0.5 * (ag_sigma[i] + probe_radius[i]);
      t_eps  = sqrt(ag_eps[i] * probe_well_depth[i]);
      break;
    case VdwCombiningRule::NBFIX:
      t_sig6 = probe_radius[i];
      t_eps = probe_well_depth[i];
      break;
    }
    probe_lj_sigma_ptr[i] = t_sig6;
    ref_probe_lj_sigma_ptr[i] = t_sig6;
    t_sig6 = t_sig6 * t_sig6;
    t_sig6 = t_sig6 * t_sig6 * t_sig6;
    const double lja = 4.0 * t_eps * t_sig6 * t_sig6;
    const double ljb = 4.0 * t_eps * t_sig6;
    probe_lja_ptr[i] = lja;
    probe_ljb_ptr[i] = ljb;
    ref_probe_lja_ptr[i] = lja;
    ref_probe_ljb_ptr[i] = ljb;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MeshForceField<T>::setLJCoefficients(const AtomGraph &ag,
                                          const std::vector<double> &probe_radius,
                                          const std::vector<double> &probe_well_depth) {
  setLJCoefficients(ag.getSelfPointer(), probe_radius, probe_well_depth);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MeshForceField<T>::setLJCoefficients(const AtomGraph *ag, const double probe_radius,
                                          const double probe_well_depth) {
  const size_t nt = ag->getLJTypeCount();
  setLJCoefficients(ag, std::vector<double>(nt, probe_radius),
                    std::vector<double>(nt, probe_well_depth));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MeshForceField<T>::setLJCoefficients(const AtomGraph &ag, const double probe_radius,
                                          const double probe_well_depth) {
  const size_t nt = ag.getLJTypeCount();
  setLJCoefficients(ag.getSelfPointer(), std::vector<double>(nt, probe_radius),
                    std::vector<double>(nt, probe_well_depth));
}


//-------------------------------------------------------------------------------------------------
template <typename T>
void MeshForceField<T>::setLJSoftcoreParameter(const int index, const double coef, const int pos) {
  validateTypeIndex(index);
  switch (pos) {
  case 0:
    probe_softcore_lja.putHost(coef, index);
    ref_probe_softcore_lja.putHost(coef, index);
    break;
  case 1:
    probe_softcore_ljb.putHost(coef, index);
    ref_probe_softcore_ljb.putHost(coef, index);
    break;
  case 2:
    probe_softcore_ljc.putHost(coef, index);
    ref_probe_softcore_ljc.putHost(coef, index);
    break;
  case 3:
    probe_softcore_ljd.putHost(coef, index);
    ref_probe_softcore_ljd.putHost(coef, index);
    break;
  case 4:
    probe_softcore_lje.putHost(coef, index);
    ref_probe_softcore_lje.putHost(coef, index);
    break;
  case 5:
    probe_softcore_ljf.putHost(coef, index);
    ref_probe_softcore_ljf.putHost(coef, index);
    break;
  default:
    rtErr("Lennard-Jones softcore potential functions take at most six coefficients (quintic "
          "order polynomials).", "MeshForceField", "setLJSoftcoreParameter");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MeshForceField<T>::setLJSoftcoreParameters(Interpolant stencil_kind) {

  // Set local template-type pointers
  T* probe_lja_ptr = probe_lja.data();
  T* probe_ljb_ptr = probe_ljb.data();
  T* probe_lj_sigma_ptr = probe_lj_sigma.data();

  // Set local double-precision pointers
  double* ref_probe_lja_ptr = ref_probe_lja.data();
  double* ref_probe_ljb_ptr = ref_probe_ljb.data();
  for (int i = 0; i < lj_type_count; i++) {
    const double lja = ref_probe_lja_ptr[i];
    const double ljb = ref_probe_ljb_ptr[i];

    // Compute the softcore Lennard-Jones polynomial.  For an interpolant based on mixed
    // partial derivatives at each corner, this is a quintic polynomial.
    double4 abcd_coef = { 0.0, 0.0, 0.0, 0.0 };
    double2 ef_coef = { 0.0, 0.0 };
    const double rswitch = clash_ratio * probe_lj_sigma_ptr[i];
    if (rswitch > constants::tiny) {
      const double inv_rs  = 1.0 / rswitch;
      const double inv_rs2 = inv_rs * inv_rs;
      const double inv_rs3 = inv_rs2 * inv_rs;
      const double inv_rs6 = inv_rs3 * inv_rs3;
      const double lj_switch   = ((lja * inv_rs6) - ljb) * inv_rs6;
      const double dlj_switch  = ((  -12.0 * lja * inv_rs6) + (  6.0 * ljb)) * inv_rs6 * inv_rs;
      const double d2lj_switch = ((  156.0 * lja * inv_rs6) + (-42.0 * ljb)) * inv_rs6 * inv_rs2;
      const double d3lj_switch = ((-2184.0 * lja * inv_rs6) + (336.0 * ljb)) * inv_rs6 * inv_rs3;
      switch (stencil_kind) {
      case Interpolant::SMOOTHNESS:
        quinticSoftCore(&abcd_coef, &ef_coef, rswitch, lj_switch, dlj_switch, d2lj_switch,
                        d3lj_switch, 0.0);
        break;
      case Interpolant::FUNCTION_VALUE:
        cubicSoftCore(&abcd_coef, rswitch, lj_switch, dlj_switch, 0.0);
        break;
      }
    }
    setLJSoftcoreParameter(i, abcd_coef.x, 0);
    setLJSoftcoreParameter(i, abcd_coef.y, 1);
    setLJSoftcoreParameter(i, abcd_coef.z, 2);
    setLJSoftcoreParameter(i, abcd_coef.w, 3);
    setLJSoftcoreParameter(i, ef_coef.x, 4);
    setLJSoftcoreParameter(i, ef_coef.y, 5);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void MeshForceField<T>::reallocate(const int lj_type_count_in) {
  std::vector<double> elec_sc_tmp = ref_elec_softcore.readHost();
  allocate(lj_type_count_in);
  ref_elec_softcore.putHost(elec_sc_tmp);
  elec_softcore.putHost(std::vector<T>(elec_sc_tmp.begin(), elec_sc_tmp.end()));
}

//-------------------------------------------------------------------------------------------------
template <typename T> void MeshForceField<T>::allocate(int lj_type_count_in) {

  // Record the number of Lennard-Jones types
  lj_type_count = lj_type_count_in;
  const int padded_nljt = roundUp(lj_type_count, warp_size_int);
  const int padded_elec = roundUp(6, warp_size_int);

  // Allocate and set the template-type parameter arrays
  softcore_data.resize((9 * padded_nljt) + padded_elec);
  probe_lja.setPointer(&softcore_data,                        0, lj_type_count);
  probe_ljb.setPointer(&softcore_data,              padded_nljt, lj_type_count);
  probe_lj_sigma.setPointer(&softcore_data,     2 * padded_nljt, lj_type_count);
  probe_softcore_lja.setPointer(&softcore_data, 3 * padded_nljt, lj_type_count);
  probe_softcore_ljb.setPointer(&softcore_data, 4 * padded_nljt, lj_type_count);
  probe_softcore_ljc.setPointer(&softcore_data, 5 * padded_nljt, lj_type_count);
  probe_softcore_ljd.setPointer(&softcore_data, 6 * padded_nljt, lj_type_count);
  probe_softcore_lje.setPointer(&softcore_data, 7 * padded_nljt, lj_type_count);
  probe_softcore_ljf.setPointer(&softcore_data, 8 * padded_nljt, lj_type_count);
  elec_softcore.setPointer(&softcore_data, 9 * padded_nljt, 6);

  // Allocate and set the double-precision parameter arrays
  ref_softcore_data.resize((9 * padded_nljt) + padded_elec);
  ref_probe_lja.setPointer(&ref_softcore_data,                        0, lj_type_count);
  ref_probe_ljb.setPointer(&ref_softcore_data,              padded_nljt, lj_type_count);
  ref_probe_lj_sigma.setPointer(&ref_softcore_data,     2 * padded_nljt, lj_type_count);
  ref_probe_softcore_lja.setPointer(&ref_softcore_data, 3 * padded_nljt, lj_type_count);
  ref_probe_softcore_ljb.setPointer(&ref_softcore_data, 4 * padded_nljt, lj_type_count);
  ref_probe_softcore_ljc.setPointer(&ref_softcore_data, 5 * padded_nljt, lj_type_count);
  ref_probe_softcore_ljd.setPointer(&ref_softcore_data, 6 * padded_nljt, lj_type_count);
  ref_probe_softcore_lje.setPointer(&ref_softcore_data, 7 * padded_nljt, lj_type_count);
  ref_probe_softcore_ljf.setPointer(&ref_softcore_data, 8 * padded_nljt, lj_type_count);
  ref_elec_softcore.setPointer(&ref_softcore_data, 9 * padded_nljt, 6);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void MeshForceField<T>::rebasePointers() {
  probe_lja.swapTarget(&softcore_data);
  probe_ljb.swapTarget(&softcore_data);
  probe_lj_sigma.swapTarget(&softcore_data);
  probe_softcore_lja.swapTarget(&softcore_data);
  probe_softcore_ljb.swapTarget(&softcore_data);
  probe_softcore_ljc.swapTarget(&softcore_data);
  probe_softcore_ljd.swapTarget(&softcore_data);
  probe_softcore_lje.swapTarget(&softcore_data);
  probe_softcore_ljf.swapTarget(&softcore_data);
  elec_softcore.swapTarget(&softcore_data);
  ref_probe_lja.swapTarget(&ref_softcore_data);
  ref_probe_ljb.swapTarget(&ref_softcore_data);
  ref_probe_lj_sigma.swapTarget(&ref_softcore_data);
  ref_probe_softcore_lja.swapTarget(&ref_softcore_data);
  ref_probe_softcore_ljb.swapTarget(&ref_softcore_data);
  ref_probe_softcore_ljc.swapTarget(&ref_softcore_data);
  ref_probe_softcore_ljd.swapTarget(&ref_softcore_data);
  ref_probe_softcore_lje.swapTarget(&ref_softcore_data);
  ref_probe_softcore_ljf.swapTarget(&ref_softcore_data);
  ref_elec_softcore.swapTarget(&ref_softcore_data);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void MeshForceField<T>::validateTypeIndex(int index) {
  if (index >= lj_type_count) {
    rtErr("Type index " + std::to_string(index) + " is invalid for a table of " +
          std::to_string(lj_type_count) + " Lennard-Jones types.", "MeshForceField",
          "setLJSoftcoreParameter");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> MeshFFKit<T> restoreType(const MeshFFKit<void> *rasa) {
  return MeshFFKit<T>(rasa->ljrule, rasa->clash_ratio, reinterpret_cast<const T*>(rasa->probe_lja),
                      reinterpret_cast<const T*>(rasa->probe_ljb),
                      reinterpret_cast<const T*>(rasa->probe_ljsig),
                      reinterpret_cast<const T*>(rasa->softcore_lja),
                      reinterpret_cast<const T*>(rasa->softcore_ljb),
                      reinterpret_cast<const T*>(rasa->softcore_ljc),
                      reinterpret_cast<const T*>(rasa->softcore_ljd),
                      reinterpret_cast<const T*>(rasa->softcore_lje),
                      reinterpret_cast<const T*>(rasa->softcore_ljf), rasa->clash_distance,
                      rasa->coulomb, reinterpret_cast<const T*>(rasa->softcore_qq));
}

//-------------------------------------------------------------------------------------------------
template <typename T> MeshFFKit<T> restoreType(const MeshFFKit<void> &rasa) {
  return MeshFFKit<T>(rasa.ljrule, rasa.clash_ratio, reinterpret_cast<const T*>(rasa.probe_lja),
                      reinterpret_cast<const T*>(rasa.probe_ljb),
                      reinterpret_cast<const T*>(rasa.probe_ljsig),
                      reinterpret_cast<const T*>(rasa.softcore_lja),
                      reinterpret_cast<const T*>(rasa.softcore_ljb),
                      reinterpret_cast<const T*>(rasa.softcore_ljc),
                      reinterpret_cast<const T*>(rasa.softcore_ljd),
                      reinterpret_cast<const T*>(rasa.softcore_lje),
                      reinterpret_cast<const T*>(rasa.softcore_ljf), rasa.clash_distance,
                      rasa.coulomb, reinterpret_cast<const T*>(rasa.softcore_qq));
}

} // namespace structure
} // namespace stormm
