// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace topology {

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> AtomGraph::getPartialCharge() const {
  return getPartialCharge<T>(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> AtomGraph::getPartialCharge(const int low_index, const int high_index) const {
  atomValidityCheck(low_index, high_index, atom_count, "AtomGraph", "getPartialCharge");
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    std::vector<float> tmpq = sp_atomic_charges.readHost(low_index, high_index - low_index);
    return std::vector<T>(tmpq.begin(), tmpq.end());
  }
  else if (ct == double_type_index) {
    std::vector<double> tmpq = atomic_charges.readHost(low_index, high_index - low_index);
    return std::vector<T>(tmpq.begin(), tmpq.end());
  }
  else {
    rtErr("Invalid request for partial charges in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "AtomGraph", "getPartialCharge");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T AtomGraph::getPartialCharge(const int index) const {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    return sp_atomic_charges.readHost(index);
  }
  else if (ct == double_type_index) {
    return atomic_charges.readHost(index);
  }
  else {
    rtErr("Invalid request for partial charge in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "AtomGraph", "getPartialCharge");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> AtomGraph::getAtomicMass(const MassForm rep) const {
  return getAtomicMass<T>(0, atom_count, rep);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> AtomGraph::getAtomicMass(const int low_index, const int high_index,
                                        const MassForm rep) const {
  atomValidityCheck(low_index, high_index, atom_count, "AtomGraph", "getAtomicMass");
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const int* znum_ptr = atomic_numbers.data();
  if (ct == float_type_index) {
    std::vector<float> tmpm;
    switch (rep) {
    case MassForm::ORDINARY:
      tmpm = sp_atomic_masses.readHost(low_index, high_index - low_index);
      break;
    case MassForm::INVERSE:
      tmpm = sp_inverse_atomic_masses.readHost(low_index, high_index - low_index);
      break;
    }
    for (int i = low_index; i < high_index; i++) {
      tmpm[i] *= (znum_ptr[i] > 0);
    }
    return std::vector<T>(tmpm.begin(), tmpm.end());
  }
  else if (ct == double_type_index) {
    std::vector<double> tmpm;
    switch (rep) {
    case MassForm::ORDINARY:
      tmpm = atomic_masses.readHost(low_index, high_index - low_index);
      break;
    case MassForm::INVERSE:
      tmpm = inverse_atomic_masses.readHost(low_index, high_index - low_index);
      break;
    }
    for (int i = low_index; i < high_index; i++) {
      tmpm[i] *= (znum_ptr[i] > 0);
    }
    return std::vector<T>(tmpm.begin(), tmpm.end());
  }
  else {
    rtErr("Invalid request for atomic masses in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "AtomGraph", "getAtomicMass");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T AtomGraph::getAtomicMass(const int index, const MassForm rep) const {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    switch (rep) {
    case MassForm::ORDINARY:
      return sp_atomic_masses.readHost(index) * (atomic_numbers.readHost(index) > 0);
    case MassForm::INVERSE:
      return sp_inverse_atomic_masses.readHost(index) * (atomic_numbers.readHost(index) > 0);
    }
  }
  else if (ct == double_type_index) {
    switch (rep) {
    case MassForm::ORDINARY:
      return atomic_masses.readHost(index) * (atomic_numbers.readHost(index) > 0);
    case MassForm::INVERSE:
      return inverse_atomic_masses.readHost(index) * (atomic_numbers.readHost(index) > 0);
    }
  }
  else {
    rtErr("Invalid request for atomic mass in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "AtomGraph", "getAtomicMass");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> UreyBradleyTerm<T> AtomGraph::getUreyBradleyTerm(const int index) const {
  if (index < 0 || index > urey_bradley_term_count) {
    rtErr("Topology \"" + title + "\" contains " + std::to_string(urey_bradley_term_count) +
          " Urey-Bradley terms.  Index " + std::to_string(index) + " is invalid.\n", "AtomGraph",
          "getUreyBradleyTerm");
  }
  UreyBradleyTerm<T> result;
  result.i_atom    = urey_bradley_i_atoms.readHost(index);
  result.k_atom    = urey_bradley_k_atoms.readHost(index);
  result.param_idx = urey_bradley_parameter_indices.readHost(index);

  // The only allowed types for the potential parameters are float and double
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    result.keq = sp_urey_bradley_stiffnesses.readHost(result.param_idx);
    result.leq = sp_urey_bradley_equilibria.readHost(result.param_idx);
  }
  else if (ct == double_type_index) {
    result.keq = urey_bradley_stiffnesses.readHost(result.param_idx);
    result.leq = urey_bradley_equilibria.readHost(result.param_idx);
  }
  else {
    rtErr("Invalid request for Urey-Bradley parameters in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "AtomGraph", "getUreyBradleyTerm");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> CharmmImprTerm<T> AtomGraph::getCharmmImprTerm(const int index) const {
  if (index < 0 || index > charmm_impr_term_count) {
    rtErr("Topology \"" + title + "\" contains " + std::to_string(charmm_impr_term_count) +
          " CHARMM improper terms.  Index " + std::to_string(index) + " is invalid.\n",
          "AtomGraph", "getCharmmImprTerm");
  }
  CharmmImprTerm<T> result;
  result.i_atom    = charmm_impr_i_atoms.readHost(index);
  result.j_atom    = charmm_impr_j_atoms.readHost(index);
  result.k_atom    = charmm_impr_k_atoms.readHost(index);
  result.l_atom    = charmm_impr_l_atoms.readHost(index);
  result.param_idx = charmm_impr_parameter_indices.readHost(index);

  // The only allowed types for the potential parameters are float and double
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    result.keq    = sp_charmm_impr_stiffnesses.readHost(result.param_idx);
    result.phi_eq = sp_charmm_impr_phase_angles.readHost(result.param_idx);
  }
  else if (ct == double_type_index) {
    result.keq    = charmm_impr_stiffnesses.readHost(result.param_idx);
    result.phi_eq = charmm_impr_phase_angles.readHost(result.param_idx);
  }
  else {
    rtErr("Invalid request for CHARMM improper dihedral parameters in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "AtomGraph", "getCharmmImprTerm");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> CmapTerm<T> AtomGraph::getCmapTerm(const int index) const {
  if (index < 0 || index > cmap_term_count) {
    rtErr("Topology \"" + title + "\" contains " + std::to_string(cmap_term_count) +
          " CMAP terms.  Index " + std::to_string(index) + " is invalid.\n", "AtomGraph",
          "getCmapTerm");
  }
  CmapTerm<T> result;
  result.i_atom   = cmap_i_atoms.readHost(index);
  result.j_atom   = cmap_j_atoms.readHost(index);
  result.k_atom   = cmap_k_atoms.readHost(index);
  result.l_atom   = cmap_l_atoms.readHost(index);
  result.m_atom   = cmap_m_atoms.readHost(index);
  result.surf_idx = cmap_surface_indices.readHost(index);
  result.surf_dim = cmap_surface_dimensions.readHost(result.surf_idx);

  // The only allowed types for the potential parameters are float and double
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    result.surf = (T*)&sp_cmap_surfaces.data()[cmap_surface_bounds.readHost(result.surf_idx)];
  }
  else if (ct == double_type_index) {
    result.surf = (T*)&cmap_surfaces.data()[cmap_surface_bounds.readHost(result.surf_idx)];
  }
  else {
    rtErr("Invalid request for CMAP parameters in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "AtomGraph", "getCmapTerm");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> BondTerm<T> AtomGraph::getBondTerm(const int index) const {
  if (index < 0 || index > bond_term_count) {
    rtErr("Topology \"" + title + "\" contains " + std::to_string(bond_term_count) + " bond "
          "stretching terms.  Index " + std::to_string(index) + " is invalid.\n", "AtomGraph",
          "getBondTerm");
  }
  BondTerm<T> result;
  result.i_atom    = bond_i_atoms.readHost(index);
  result.j_atom    = bond_j_atoms.readHost(index);
  result.param_idx = bond_parameter_indices.readHost(index);

  // The only allowed types for the potential parameters are float and double
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    result.keq = sp_bond_stiffnesses.readHost(result.param_idx);
    result.leq = sp_bond_equilibria.readHost(result.param_idx);
  }
  else if (ct == double_type_index) {
    result.keq = bond_stiffnesses.readHost(result.param_idx);
    result.leq = bond_equilibria.readHost(result.param_idx);
  }
  else {
    rtErr("Invalid request for bond parameters in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "AtomGraph", "getBondTerm");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> AngleTerm<T> AtomGraph::getAngleTerm(const int index) const {
  if (index < 0 || index > angl_term_count) {
    rtErr("Topology \"" + title + "\" contains " + std::to_string(angl_term_count) + " angle "
          "bending terms.  Index " + std::to_string(index) + " is invalid.\n", "AtomGraph",
          "getAngleTerm");
  }
  AngleTerm<T> result;
  result.i_atom    = angl_i_atoms.readHost(index);
  result.j_atom    = angl_j_atoms.readHost(index);
  result.k_atom    = angl_j_atoms.readHost(index);
  result.param_idx = angl_parameter_indices.readHost(index);

  // The only allowed types for the potential parameters are float and double
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    result.keq = sp_angl_stiffnesses.readHost(result.param_idx);
    result.theta_eq = sp_angl_equilibria.readHost(result.param_idx);
  }
  else if (ct == double_type_index) {
    result.keq = angl_stiffnesses.readHost(result.param_idx);
    result.theta_eq = angl_equilibria.readHost(result.param_idx);
  }
  else {
    rtErr("Invalid request for bond angle parameters in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "AtomGraph", "getAngleTerm");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> DihedralTerm<T> AtomGraph::getDihedralTerm(const int index) const {
  if (index < 0 || index > dihe_term_count) {
    rtErr("Topology \"" + title + "\" contains " + std::to_string(angl_term_count) + " dihedral "
          "torsion terms.  Index " + std::to_string(index) + " is invalid.\n", "AtomGraph",
          "getDihedralTerm");
  }
  DihedralTerm<T> result;
  result.i_atom    = dihe_i_atoms.readHost(index);
  result.j_atom    = dihe_j_atoms.readHost(index);
  result.k_atom    = dihe_j_atoms.readHost(index);
  result.l_atom    = dihe_l_atoms.readHost(index);
  result.param_idx = dihe_parameter_indices.readHost(index);

  // The only allowed types for the potential parameters are float and double
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const int attn14_index = dihe14_parameter_indices.readHost(index);
  if (ct == float_type_index) {
    result.amplitude   = sp_dihe_amplitudes.readHost(result.param_idx);
    result.phase       = sp_dihe_phase_angles.readHost(result.param_idx);
    result.periodicity = sp_dihe_periodicities.readHost(result.param_idx);
    result.elec_screen = sp_attn14_elec_factors.readHost(attn14_index);
    result.vdw_screen  = sp_attn14_vdw_factors.readHost(attn14_index);
  }
  else if (ct == double_type_index) {
    result.amplitude   = dihe_amplitudes.readHost(result.param_idx);
    result.phase       = dihe_phase_angles.readHost(result.param_idx);
    result.periodicity = dihe_periodicities.readHost(result.param_idx);
    result.elec_screen = attn14_elec_factors.readHost(attn14_index);
    result.vdw_screen  = attn14_vdw_factors.readHost(attn14_index);
  }
  else {
    rtErr("Invalid request for dihedral parameters in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "AtomGraph", "getDihedralTerm");
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> T AtomGraph::getVirtualSiteFrameDimension(const int index, const int ndim) {
  bool problem = false;
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    switch (ndim) {
    case 1:
      return sp_virtual_site_frame_dim1.readHost(index);
    case 2:
      return sp_virtual_site_frame_dim2.readHost(index);
    case 3:
      return sp_virtual_site_frame_dim3.readHost(index);
    default: problem = true;
      break;
    }
  }
  else if (ct == double_type_index) {
    switch (ndim) {
    case 1:
      return virtual_site_frame_dim1.readHost(index);
    case 2:
      return virtual_site_frame_dim2.readHost(index);
    case 3:
      return virtual_site_frame_dim3.readHost(index);
    default: problem = true;
      break;
    }
  }
  else {
    rtErr("Invalid request for virtual site frame dimensions in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "AtomGraph", "getUreyBradleyTerm");
  }
  if (problem) {
    rtErr("Virtual sites can contain up to three frame dimensions.  Frame dimension " +
          std::to_string(ndim) + " is invalid.", "getVirtualSiteFrameDimension", "AtomGraph");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
T AtomGraph::getChargeParameter(const int index) const {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    return sp_atomic_charges.readHost(index);
  }
  else if (ct == double_type_index) {
    return atomic_charges.readHost(index);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
T AtomGraph::getLennardJonesSigma(const int index_a, const int index_b) const {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    return sp_lj_sigma_values.readHost((index_a * lj_type_count) + index_b);
  }
  else if (ct == double_type_index) {
    return lj_sigma_values.readHost((index_a * lj_type_count) + index_b);
  }
  else {
    rtErr("Lennard-Jones sigma values may be specified as single- or double-precision real "
          "numbers.", "getLennardJonesSigma", "AtomGraph");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T AtomGraph::getLennardJonesSigma(const int index_a) const {
  return getLennardJonesSigma<T>(index_a, index_a);
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> AtomGraph::getLennardJonesSigma() const {
  std::vector<T> result(lj_type_count);
  for (int i = 0; i < lj_type_count; i++) {
    result[i] = getLennardJonesSigma<T>(i, i);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
T AtomGraph::getLennardJonesEpsilon(const int index_a, const int index_b) const {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    const double lj_a = sp_lj_a_values.readHost((index_a * lj_type_count) + index_b);
    const double lj_b = sp_lj_b_values.readHost((index_a * lj_type_count) + index_b);
    return (lj_a < constants::tiny) ? 0.0 : static_cast<T>(0.25 * lj_b * lj_b / lj_a);
  }
  else if (ct == double_type_index) {
    const double lj_a = lj_a_values.readHost((index_a * lj_type_count) + index_b);
    const double lj_b = lj_b_values.readHost((index_a * lj_type_count) + index_b);
    return (lj_a < constants::tiny) ? 0.0 : static_cast<T>(0.25 * lj_b * lj_b / lj_a);
  }
  else {
    rtErr("Lennard-Jones epsilon values may be specified as single- or double-precision real "
          "numbers.", "getLennardJonesEpsilon", "AtomGraph");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T AtomGraph::getLennardJonesEpsilon(const int index_a) const {
  return getLennardJonesEpsilon<T>(index_a, index_a);
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> AtomGraph::getLennardJonesEpsilon() const {
  std::vector<T> result(lj_type_count);
  for (int i = 0; i < lj_type_count; i++) {
    result[i] = getLennardJonesEpsilon<T>(i, i);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> AtomGraph::getAtomPBRadius() const {
  return getAtomPBRadius<T>(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> AtomGraph::getAtomPBRadius(const int low_index, const int high_index) const {
  atomValidityCheck(low_index, high_index, atom_count, "AtomGraph", "getAtomPBRadius");
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    std::vector<float> tmpr = sp_atomic_pb_radii.readHost(low_index, high_index - low_index);
    return std::vector<T>(tmpr.begin(), tmpr.end());
  }
  else if (ct == double_type_index) {
    std::vector<double> tmpr = atomic_pb_radii.readHost(low_index, high_index - low_index);
    return std::vector<T>(tmpr.begin(), tmpr.end());
  }
  else {
    rtErr("Invalid request for atomic PB Radii in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "AtomGraph", "getAtomPBRadius");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T AtomGraph::getAtomPBRadius(const int index) const {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    return sp_atomic_pb_radii.readHost(index);
  }
  else if (ct == double_type_index) {
    return atomic_pb_radii.readHost(index);
  }
  else {
    rtErr("Invalid request for atomic PB Radius in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "AtomGraph", "getAtomPBRadius");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> AtomGraph::getGBScreeningFactor() const {
  return getGBScreeningFactor<T>(0, atom_count);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> AtomGraph::getGBScreeningFactor(const int low_index, const int high_index) const {
  atomValidityCheck(low_index, high_index, atom_count, "AtomGraph", "getGBScreeningFactor");
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    std::vector<float> tmpf = sp_gb_screening_factors.readHost(low_index, high_index - low_index);
    return std::vector<T>(tmpf.begin(), tmpf.end());
  }
  else if (ct == double_type_index) {
    std::vector<double> tmpf = gb_screening_factors.readHost(low_index, high_index - low_index);
    return std::vector<T>(tmpf.begin(), tmpf.end());
  }
  else {
    rtErr("Invalid request for atomic GB screening factors in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "AtomGraph",
          "getGBScreeningFactor");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> T AtomGraph::getGBScreeningFactor(const int index) const {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == float_type_index) {
    return sp_gb_screening_factors.readHost(index);
  }
  else if (ct == double_type_index) {
    return gb_screening_factors.readHost(index);
  }
  else {
    rtErr("Invalid request for atomic GB screening factor in format " +
          std::string(std::type_index(typeid(T)).name()) + ".", "AtomGraph",
          "getGBScreeningFactor");
  }
  __builtin_unreachable();
}

} // namespace topology
} // namespace stormm
