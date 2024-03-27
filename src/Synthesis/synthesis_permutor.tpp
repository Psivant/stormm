// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
template <typename T>
SyPermutorKit<T>::SyPermutorKit(int system_count_in, int perm_map_count_in,
                                const int* perm_map_idx_in, const int2* perm_elements_in,
                                const int* perm_element_bounds_in, const int* system_settings_in,
                                const int* system_settings_limits_in, const int* rot_grp_atoms_in,
                                const int* rot_grp_bounds_in, const int* prm_rot_grp_bounds_in,
                                const int* ctx_grp_atoms_in, const int* ctx_grp_bounds_in,
                                const int* prm_ctx_grp_bounds_in, const int* inv_grp_atoms_in,
                                const int* inv_grp_bounds_in, const int* prm_inv_grp_bounds_in,
                                const int* chiral_atoms_in, const int* chiral_protocols_in,
                                const int4* rot_bond_markers_in, const int4* ctx_bond_markers_in,
                                const int4* chiral_markers_in, const T* rot_bond_settings_in,
                                const T* ctx_bond_settings_in,
                                const int* rot_bond_settings_bounds_in,
                                const int* ctx_bond_settings_bounds_in,
                                const int* chiral_settings_in,
                                const int* chiral_settings_bounds_in) :
    system_count{system_count_in}, perm_map_count{perm_map_count_in},
    perm_map_idx{perm_map_idx_in}, perm_elements{perm_elements_in},
    perm_element_bounds{perm_element_bounds_in}, system_settings{system_settings_in},
    system_settings_limits{system_settings_limits_in}, rot_grp_atoms{rot_grp_atoms_in},
    rot_grp_bounds{rot_grp_bounds_in}, prm_rot_grp_bounds{prm_rot_grp_bounds_in},
    ctx_grp_atoms{ctx_grp_atoms_in}, ctx_grp_bounds{ctx_grp_bounds_in},
    prm_ctx_grp_bounds{prm_ctx_grp_bounds_in}, inv_grp_atoms{inv_grp_atoms_in},
    inv_grp_bounds{inv_grp_bounds_in}, prm_inv_grp_bounds{prm_inv_grp_bounds_in},
    chiral_atoms{chiral_atoms_in}, chiral_protocols{chiral_protocols_in},
    rot_bond_markers{rot_bond_markers_in}, ctx_bond_markers{ctx_bond_markers_in},
    chiral_markers{chiral_markers_in}, rot_bond_settings{rot_bond_settings_in},
    ctx_bond_settings{ctx_bond_settings_in}, rot_bond_settings_bounds{rot_bond_settings_bounds_in},
    ctx_bond_settings_bounds{ctx_bond_settings_bounds_in},
    chiral_settings_bounds{chiral_settings_bounds_in}, chiral_settings{chiral_settings_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
void SynthesisPermutor::alterMapSettings(const int map_index, Hybrid<T> *variable_settings,
                                         const Hybrid<int> &variable_settings_bounds,
                                         const Hybrid<int> &permutor_map_bounds,
                                         const std::vector<std::vector<T>> &settings) {
  validateMapIndex(map_index);

  // Check that the array of settings conforms to the map's current shape.
  const int map_llim = permutor_map_bounds.readHost(map_index);
  const int nvar = permutor_map_bounds.readHost(map_index + 1) - map_llim;
  if (static_cast<int>(settings.size()) != nvar) {
    rtErr("The number of variables in the map (" + std::to_string(nvar) + ") does not agree with "
          "input settings for " + std::to_string(settings.size()) + " variables.",
          "SynthesisPermutor", "defineCisTransBondSettings");
  }
  for (int i = 0; i < nvar; i++) {
    const int llim = variable_settings_bounds.readHost(map_llim + i);
    const int hlim = variable_settings_bounds.readHost(map_llim + i + 1);
    if (static_cast<int>(settings[i].size()) != hlim - llim) {
      rtErr("The number of states available to rotatable bond " + std::to_string(i) + " of map " +
            std::to_string(map_index) + ", " + std::to_string(hlim - llim) + ", does not match "
            "the input of " + std::to_string(settings[i].size()) + " states.", "SynthesisPermutor",
            "defineRotatableBondSettings");
    }
    T* vs_ptr = variable_settings->data();
    for (int j = llim; j < hlim; j++) {
      vs_ptr[j] = settings[i][j - llim];
    }    
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void SynthesisPermutor::permuteSystem(CoordinateSeries<T> *cs, const int map_index,
                                      const CoupledEdit ce, const int setting_index) const {
  const bool tcalc_is_double = (std::type_index(typeid(T)).hash_code() == double_type_index);
  switch (ce.edit) {
  case ConformationEdit::BOND_ROTATION:
    if (ce.index < 0 || ce.index >= total_rotatable_groups) {
      rtErr("Rotatable bond index " + std::to_string(ce.index) + " is invalid for a series of " +
            std::to_string(total_rotatable_groups) + " spanning " +
            std::to_string(permutor_map_count) + " topologies.", "SynthesisPermutor",
            "permuteSystem");
    }
    else {

      // Get the atoms at either end of the bond.
      const int* rotg_atom_ptr = rotatable_group_atoms.data();
      const int g_offset = rotatable_group_bounds.readHost(ce.index);
      const int g_limit = rotatable_group_bounds.readHost(ce.index + 1);
      const int root_atom  = rotg_atom_ptr[g_offset];
      const int pivot_atom = rotg_atom_ptr[g_offset + 1];
      const int4 atom_defs = rotatable_bond_markers.readHost(ce.index);
      const int setting_offset = rotatable_bond_settings_bounds.readHost(ce.index);
      const T current_angle = dihedralAngle<T, T>(atom_defs.x, atom_defs.y, atom_defs.z,
                                                  atom_defs.w, cs, 0);
      const T target_angle = rotatable_bond_settings.readHost(setting_offset + setting_index);
      const T delta = target_angle - current_angle;
      if (fabs(delta) >= 1.0e-5) {
        rotateAboutBond<T, T>(cs->data(), 0, root_atom, pivot_atom, &rotg_atom_ptr[g_offset + 2],
                              g_limit - g_offset - 2, delta);
      }
    }
    break;
  case ConformationEdit::CIS_TRANS_FLIP:
    if (ce.index < 0 || ce.index >= total_cis_trans_groups) {
      rtErr("Cis-trans isomeric bond index " + std::to_string(ce.index) + " is invalid for a "
            "series of " + std::to_string(total_cis_trans_groups) + " spanning " +
            std::to_string(permutor_map_count) + " topologies.", "SynthesisPermutor",
            "permuteSystem");
    }
    else {

      // Get the atoms at either end of the bond.
      const int* ctxg_atom_ptr = cis_trans_group_atoms.data();
      const int g_offset = cis_trans_group_bounds.readHost(ce.index);
      const int g_limit = cis_trans_group_bounds.readHost(ce.index + 1);
      const int root_atom  = ctxg_atom_ptr[g_offset];
      const int pivot_atom = ctxg_atom_ptr[g_offset + 1];
      const int4 atom_defs = cis_trans_bond_markers.readHost(ce.index);
      const int setting_offset = cis_trans_bond_settings_bounds.readHost(ce.index);
      const T target_angle = cis_trans_bond_settings.readHost(setting_offset + setting_index);
      const T current_angle = dihedralAngle<T, T>(atom_defs.x, atom_defs.y, atom_defs.z,
                                                  atom_defs.w, cs, 0);
      const T delta = target_angle - current_angle;
      if (fabs(delta) >= 1.0e-5) {
        rotateAboutBond<T, T>(cs->data(), 0, root_atom, pivot_atom, &ctxg_atom_ptr[g_offset + 2],
                              g_limit - g_offset - 2, delta);
      }
    }
    break;
  case ConformationEdit::CHIRAL_INVERSION:
    if (ce.index < 0 || ce.index >= total_invertible_groups) {
      rtErr("Chiral center index " + std::to_string(ce.index) + " is invalid for a "
            "series of " + std::to_string(total_invertible_groups) + " spanning " +
            std::to_string(permutor_map_count) + " topologies.", "SynthesisPermutor",
            "permuteSystem");
    }

    // Check that the system is NOT in the state that is required before flipping it.
    const int* invac_ptr = invertible_atom_centers.data();
    const int* invgp_ptr = invertible_group_protocols.data();
    const ChiralInversionProtocol cen_prot =
      static_cast<ChiralInversionProtocol>(invgp_ptr[ce.index]);
    const int setting_offset = chiral_settings_bounds.readHost(ce.index);
    const int cen_target = chiral_settings.readHost(setting_offset + setting_index);
    const int4 cen_branches = chiral_markers.readHost(ce.index);
    const int cen_chiral_state = getChiralOrientation<T>(cs, 0, invac_ptr[ce.index],
                                                         cen_branches.x, cen_branches.y,
                                                         cen_branches.z, cen_branches.w);
    if (cen_prot != ChiralInversionProtocol::DO_NOT_INVERT && cen_target != cen_chiral_state) {

      // The whole array of chiral centers becomes relevant due to the recursive nature of the
      // CPU function for inverting chirality if the protocol is to "REFLECT" (and only one center
      // per molecule may be subject to reflection).  If a center undergoes reflection, the
      // chirality of every center in the molecule will flip.  Therefore, after a reflection
      // procedure, every other chiral center (excluding those marked "DO_NOT_INVERT") must be
      // inverted using a rotational procedure.
      const int llim = permutor_invertible_group_bounds.readHost(map_index);
      const int hlim = permutor_invertible_group_bounds.readHost(map_index + 1);
      std::vector<int> tmp_chiral_atoms(hlim - llim);
      std::vector<ChiralInversionProtocol> tmp_chiral_protocols(hlim - llim);
      for (int i = llim; i < hlim; i++) {
        tmp_chiral_atoms[i - llim] = invac_ptr[i];
        tmp_chiral_protocols[i - llim] = static_cast<ChiralInversionProtocol>(invgp_ptr[i]);
      }

      // Rebuild the isomerization plans based on the synthesis permutor's tables.  This is not as
      // efficient as possible, but the cost of copying the plans is still minor compared to the
      // cost of manipulating the coordinates, and this process will help to verify the permutor
      // maps in C++ code.
      std::vector<IsomerPlan> tmp_chiral_plans;
      tmp_chiral_plans.reserve(hlim - llim);
      const AtomGraph *ag_ptr = topologies[map_index];
      for (int i = llim; i < hlim; i++) {
        const int g_llim = invertible_group_bounds.readHost(i);
        const int g_hlim = invertible_group_bounds.readHost(i + 1);
        const std::vector<int> moving_atoms = invertible_group_atoms.readHost(g_llim + 2,
                                                                              g_hlim - g_llim - 2);
        tmp_chiral_plans.emplace_back(ConformationEdit::CHIRAL_INVERSION,
                                      static_cast<ChiralInversionProtocol>(invgp_ptr[i]),
                                      invertible_group_atoms.readHost(g_llim),
                                      invertible_group_atoms.readHost(g_llim + 1), moving_atoms,
                                      ag_ptr);
      }
      flipChiralCenter<T, T>(cs, 0, ce.index - llim, tmp_chiral_atoms, tmp_chiral_protocols,
                             tmp_chiral_plans);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void SynthesisPermutor::randomizeSystem(CoordinateSeries<T> *cs, const int map_index,
                                        Xoshiro256ppGenerator *xrs, const CoupledEdit mask_a,
                                        const CoupledEdit mask_b, const CoupledEdit mask_c) const {

  // Loop over rotatable bonds.  Skip the operation if the feature is one of up to three masked
  // features (the masks' indices will be set to -1 if they are to be disregarded), or if there is
  // only one valid choice for the settings in that feature.
  const int rotg_llim = permutor_rotatable_group_bounds.readHost(map_index);
  const int rotg_hlim = permutor_rotatable_group_bounds.readHost(map_index + 1);
  for (int i = rotg_llim; i < rotg_hlim; i++) {
    if ((mask_a.edit == ConformationEdit::BOND_ROTATION && mask_a.index == i) ||
        (mask_b.edit == ConformationEdit::BOND_ROTATION && mask_b.index == i) ||
        (mask_c.edit == ConformationEdit::BOND_ROTATION && mask_c.index == i)) {
      continue;
    }
    const int nchoice = rotatable_bond_settings_bounds.readHost(i + 1) -
                        rotatable_bond_settings_bounds.readHost(i);
    if (nchoice == 1) {
      continue;
    }
    const int selection = xrs->uniformRandomNumber() * static_cast<double>(nchoice);
    permuteSystem<T>(cs, map_index, { ConformationEdit::BOND_ROTATION, i },
                     selection);
  }

  // Loop over cis-trans flips
  const int ctxg_llim = permutor_cis_trans_group_bounds.readHost(map_index);
  const int ctxg_hlim = permutor_cis_trans_group_bounds.readHost(map_index + 1);
  for (int i = ctxg_llim; i < ctxg_hlim; i++) {
    if ((mask_a.edit == ConformationEdit::CIS_TRANS_FLIP && mask_a.index == i) ||
        (mask_b.edit == ConformationEdit::CIS_TRANS_FLIP && mask_b.index == i) ||
        (mask_c.edit == ConformationEdit::CIS_TRANS_FLIP && mask_c.index == i)) {
      continue;
    }
    const int nchoice = cis_trans_bond_settings_bounds.readHost(i + 1) -
                        cis_trans_bond_settings_bounds.readHost(i);
    if (nchoice == 1) {
      continue;
    }
    const int selection = xrs->uniformRandomNumber() * static_cast<double>(nchoice);
    permuteSystem<T>(cs, map_index, { ConformationEdit::CIS_TRANS_FLIP, i },
                     selection);
  }

  // Loop over chiral inversions
  const int invg_llim = permutor_invertible_group_bounds.readHost(map_index);
  const int invg_hlim = permutor_invertible_group_bounds.readHost(map_index + 1);
  for (int i = invg_llim; i < invg_hlim; i++) {
    if ((mask_a.edit == ConformationEdit::CHIRAL_INVERSION && mask_a.index == i) ||
        (mask_b.edit == ConformationEdit::CHIRAL_INVERSION && mask_b.index == i) ||
        (mask_c.edit == ConformationEdit::CHIRAL_INVERSION && mask_c.index == i)) {
      continue;
    }
    const int nchoice = chiral_settings_bounds.readHost(i + 1) -
                        chiral_settings_bounds.readHost(i);
    if (nchoice == 1) {
      continue;
    }
    const int selection = xrs->uniformRandomNumber() * static_cast<double>(nchoice);
    permuteSystem<T>(cs, map_index, { ConformationEdit::CHIRAL_INVERSION, i },
                     selection);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
bool SynthesisPermutor::resolveClashes(CoordinateSeries<T> *cs, const int map_index,
                                       Xoshiro256ppGenerator *xrs, const StaticExclusionMask &excl,
                                       ClashReport *summary, const int iteration_limit,
                                       const int max_clashes, const CoupledEdit mask_a,
                                       const CoupledEdit mask_b, const CoupledEdit mask_c) const {
  if (summary != nullptr) {
    bool clashing;
    int iter = 0;
    do {
      detectClash<T, T>(cs, 0, summary->getTopologyPointer(), &excl, summary->getMinimumDistance(),
                        summary->getMinimumSigmaRatio(), summary);
      clashing = (summary->getClashCount() > max_clashes);
      if (clashing) {
        randomizeSystem<T>(cs, map_index, xrs, mask_a, mask_b, mask_c);
      }
      iter++;
    } while (clashing && iter < iteration_limit);
    return (clashing == false);
  }
  else {
    return true;
  }
  __builtin_unreachable();
}

} // namespace synthesis
} // namespace stormm
