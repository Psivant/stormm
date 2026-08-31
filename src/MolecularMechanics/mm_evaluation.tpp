// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace mm {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
void evalValeMM(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                const double* invu, const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                Tforce* zfrc, ScoreCard *sc, const ValenceKit<Tcalc> &vk,
                const NonbondedKit<Tcalc> &nbk, const EvaluateForce eval_force,
                const int system_index, const Tcalc inv_gpos_factor, const Tcalc force_factor,
                const Tcalc clash_distance, const Tcalc clash_ratio) {
  evaluateBondTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc,
                                           zfrc, sc, eval_force, system_index, inv_gpos_factor,
                                           force_factor);
  evaluateAngleTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                            yfrc, zfrc, sc, eval_force, system_index,
                                            inv_gpos_factor, force_factor);
  evaluateDihedralTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                               yfrc, zfrc, sc, eval_force, system_index,
                                               inv_gpos_factor, force_factor);
  evaluateUreyBradleyTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell,
                                                  xfrc, yfrc, zfrc, sc, eval_force, system_index,
                                                  inv_gpos_factor, force_factor);
  evaluateCharmmImproperTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell,
                                                     xfrc, yfrc, zfrc, sc, eval_force,
                                                     system_index, inv_gpos_factor, force_factor);
  evaluateCmapTerms<Tcoord, Tforce, Tcalc>(vk, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc,
                                           zfrc, sc, eval_force, system_index, inv_gpos_factor,
                                           force_factor);
  evaluateAttenuated14Terms<Tcoord, Tforce, Tcalc>(vk, nbk, xcrd, ycrd, zcrd, umat, invu,
                                                   unit_cell, xfrc, yfrc, zfrc, sc, eval_force,
                                                   eval_force, system_index, inv_gpos_factor,
                                                   force_factor, clash_distance, clash_ratio);
}

//-----------------------------------------------------------------------------------------------
template <typename Tc, typename Tc2, typename Tc4>
void evalValeMM(PsSynthesisWriter *poly_psw, ScoreCard* sc, const SyValenceKit<Tc> &poly_vk,
                const SyAtomUpdateKit<Tc, Tc2, Tc4> &poly_auk,
                const EvaluateForce eval_force, const VwuTask activity,
                const double clash_distance, const double clash_ratio, const int step,
                const SyRestraintKit<Tc, Tc2, Tc4> &poly_rk) {

  // Allocate space to mock the GPU's __shared__ memory buffers
  std::vector<int2> vwu_map(vwu_abstract_length);
  std::vector<llint> sh_xfrc(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_yfrc(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_zfrc(maximum_valence_work_unit_atoms);
  std::vector<int> sh_xfrc_ovrf(maximum_valence_work_unit_atoms);
  std::vector<int> sh_yfrc_ovrf(maximum_valence_work_unit_atoms);
  std::vector<int> sh_zfrc_ovrf(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_xcrd(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_ycrd(maximum_valence_work_unit_atoms);
  std::vector<llint> sh_zcrd(maximum_valence_work_unit_atoms);
  std::vector<int> sh_xcrd_ovrf(maximum_valence_work_unit_atoms);
  std::vector<int> sh_ycrd_ovrf(maximum_valence_work_unit_atoms);
  std::vector<int> sh_zcrd_ovrf(maximum_valence_work_unit_atoms);
  std::vector<Tc> sh_charges(maximum_valence_work_unit_atoms);
  std::vector<int> sh_lj_idx(maximum_valence_work_unit_atoms);

  // Set pointers with regard to the detail present in the synthesis
  llint* xfrc_ptr = sh_xfrc.data();
  llint* yfrc_ptr = sh_yfrc.data();
  llint* zfrc_ptr = sh_zfrc.data();
  int *xfrc_ovrf_ptr, *yfrc_ovrf_ptr, *zfrc_ovrf_ptr;
  if (poly_psw->frc_bits > force_scale_nonoverflow_bits) {
    xfrc_ovrf_ptr = sh_xfrc_ovrf.data();
    yfrc_ovrf_ptr = sh_yfrc_ovrf.data();
    zfrc_ovrf_ptr = sh_zfrc_ovrf.data();
  }
  else {
    xfrc_ovrf_ptr = nullptr;
    yfrc_ovrf_ptr = nullptr;
    zfrc_ovrf_ptr = nullptr;
  }
  llint* xcrd_ptr = sh_xcrd.data();
  llint* ycrd_ptr = sh_ycrd.data();
  llint* zcrd_ptr = sh_zcrd.data();
  int *xcrd_ovrf_ptr, *ycrd_ovrf_ptr, *zcrd_ovrf_ptr;
  if (poly_psw->gpos_bits > globalpos_scale_nonoverflow_bits) {
    xcrd_ovrf_ptr = sh_xcrd_ovrf.data();
    ycrd_ovrf_ptr = sh_ycrd_ovrf.data();
    zcrd_ovrf_ptr = sh_zcrd_ovrf.data();
  }
  else {
    xcrd_ovrf_ptr = nullptr;
    ycrd_ovrf_ptr = nullptr;
    zcrd_ovrf_ptr = nullptr;
  }
  
  // Determine the transform stride and extract the energy scaling factor.
  const int xfrm_stride = roundUp(9, warp_size_int);
  const Tc nrg_scale_factor = sc->getEnergyScalingFactor<Tc>();
  
  // Loop over all valence work units.  Check that any given interaction is scheduled for
  // evaluation.
  for (int vwu_idx = 0; vwu_idx < poly_vk.nvwu; vwu_idx++) {

    // Initialize energy accumulators
    llint bond_acc = 0LL;
    llint angl_acc = 0LL;
    llint dihe_acc = 0LL;
    llint impr_acc = 0LL;
    llint ubrd_acc = 0LL;
    llint cimp_acc = 0LL;
    llint cmap_acc = 0LL;
    llint qq14_acc = 0LL;
    llint lj14_acc = 0LL;
    llint rest_acc = 0LL;

    // Extract the valence work unit's abstract.
    for (int j = 0; j < vwu_abstract_length; j++) {
      vwu_map[j] = poly_vk.vwu_abstracts[(vwu_idx * vwu_abstract_length) + j];
    }

    // Determine the system index.
    const int sys_idx = vwu_map[static_cast<size_t>(VwuAbstractMap::SYSTEM_ID)].x;

    // Import the atoms of the work unit.
    const int2 import_limits = vwu_map[static_cast<size_t>(VwuAbstractMap::IMPORT)];
    for (int j = import_limits.x; j < import_limits.y; j++) {
      const size_t synth_atom = poly_vk.vwu_imports[j];
      const size_t local_atom = j - import_limits.x;
      sh_charges[local_atom] = poly_vk.charges[synth_atom];
      sh_lj_idx[local_atom] = poly_vk.lj_idx[synth_atom];
      sh_xfrc[local_atom] = 0LL;
      sh_yfrc[local_atom] = 0LL;
      sh_zfrc[local_atom] = 0LL;
      if (poly_psw->frc_bits > force_scale_nonoverflow_bits) {
        sh_xfrc_ovrf[local_atom] = 0;
        sh_yfrc_ovrf[local_atom] = 0;
        sh_zfrc_ovrf[local_atom] = 0;
      }
      sh_xcrd[local_atom] = poly_psw->xcrd[synth_atom];
      sh_ycrd[local_atom] = poly_psw->ycrd[synth_atom];
      sh_zcrd[local_atom] = poly_psw->zcrd[synth_atom];
      if (poly_psw->gpos_bits > globalpos_scale_nonoverflow_bits) {
        sh_xcrd_ovrf[local_atom] = poly_psw->xcrd_ovrf[synth_atom];
        sh_ycrd_ovrf[local_atom] = poly_psw->ycrd_ovrf[synth_atom];
        sh_zcrd_ovrf[local_atom] = poly_psw->zcrd_ovrf[synth_atom];
      }
    }
    double* umat = &poly_psw->umat[xfrm_stride * sys_idx];
    double* invu = &poly_psw->invu[xfrm_stride * sys_idx];

    // Perform bond energy calculation
    if (activity == VwuTask::BOND || activity == VwuTask::UBRD || activity == VwuTask::ALL_TASKS) {
      const int2 cbnd_limits = vwu_map[static_cast<size_t>(VwuAbstractMap::CBND)];
      const int2 nrg_limits = vwu_map[static_cast<int>(VwuAbstractMap::CBND_NRG)];
      for (int pos = cbnd_limits.x; pos < cbnd_limits.y; pos++) {
        const uint2 t_insr = poly_vk.cbnd_insr[pos];
        const bool is_urey_bradley = ((t_insr.x >> 20) & 0x1);
        const size_t i_atom = (t_insr.x & 0x3ff);
        const size_t j_atom = ((t_insr.x >> 10) & 0x3ff);
        if ((activity == VwuTask::BOND && is_urey_bradley) ||
            (activity == VwuTask::UBRD && (! is_urey_bradley))) {
          continue;
        }
        const Tc keq = (is_urey_bradley) ? poly_vk.ubrd_keq[t_insr.y] : poly_vk.bond_keq[t_insr.y];
        const Tc leq = (is_urey_bradley) ? poly_vk.ubrd_leq[t_insr.y] :
                                           std::abs(poly_vk.bond_leq[t_insr.y]);
        const Tc du = evalHarmonicStretch<llint,
                                          llint, Tc>(i_atom, j_atom, keq, leq, xcrd_ptr, ycrd_ptr,
                                                     zcrd_ptr, umat, invu, poly_psw->unit_cell,
                                                     xfrc_ptr, yfrc_ptr, zfrc_ptr, eval_force,
                                                     poly_psw->inv_gpos_scale, poly_psw->frc_scale,
                                                     xcrd_ovrf_ptr, ycrd_ovrf_ptr, zcrd_ovrf_ptr,
                                                     xfrc_ovrf_ptr, yfrc_ovrf_ptr, zfrc_ovrf_ptr);
        if (readBitFromMask(&poly_vk.cbnd_acc[nrg_limits.x], pos - cbnd_limits.x) == 1) {
          if (is_urey_bradley) {
            ubrd_acc += llround(du * nrg_scale_factor);
          }
          else {
            bond_acc += llround(du * nrg_scale_factor);
          }
        }
      }
    }

    // Perform angle energy calculation
    if (activity == VwuTask::ANGL || activity == VwuTask::ALL_TASKS) {
      const int2 angl_limits = vwu_map[static_cast<size_t>(VwuAbstractMap::ANGL)];
      const int2 nrg_limits = vwu_map[static_cast<int>(VwuAbstractMap::ANGL_NRG)];
      for (int pos = angl_limits.x; pos < angl_limits.y; pos++) {
        const uint2 t_insr = poly_vk.angl_insr[pos];
        const size_t i_atom = (t_insr.x & 0x3ff);
        const size_t j_atom = ((t_insr.x >> 10) & 0x3ff);
        const size_t k_atom = ((t_insr.x >> 20) & 0x3ff);
        const Tc du = evalHarmonicBend<llint,
                                       llint, Tc>(i_atom, j_atom, k_atom,
                                                  poly_vk.angl_keq[t_insr.y],
                                                  poly_vk.angl_theta[t_insr.y], xcrd_ptr, ycrd_ptr,
                                                  zcrd_ptr, umat, invu, poly_psw->unit_cell,
                                                  xfrc_ptr, yfrc_ptr, zfrc_ptr, eval_force,
                                                  poly_psw->inv_gpos_scale, poly_psw->frc_scale,
                                                  xcrd_ovrf_ptr, ycrd_ovrf_ptr, zcrd_ovrf_ptr,
                                                  xfrc_ovrf_ptr, yfrc_ovrf_ptr, zfrc_ovrf_ptr);
        if (readBitFromMask(&poly_vk.angl_acc[nrg_limits.x], pos - angl_limits.x) == 1) {
          angl_acc += llround(du * nrg_scale_factor);
        }
      }
    }

    // Perform dihedral energy calculation
    if (activity == VwuTask::DIHE || activity == VwuTask::CIMP || activity == VwuTask::INFR14 ||
        activity == VwuTask::ALL_TASKS) {
      const int2 cdhe_limits = vwu_map[static_cast<size_t>(VwuAbstractMap::CDHE)];
      const int2 nrg_limits = vwu_map[static_cast<int>(VwuAbstractMap::CDHE_NRG)];
      const int ljabc_offset = poly_vk.ljabc_offsets[sys_idx];
      const int nljt = poly_vk.n_lj_types[sys_idx];
      const Tc value_one = 1.0;
      for (int pos = cdhe_limits.x; pos < cdhe_limits.y; pos++) {
        const uint2 t_insr = poly_vk.cdhe_insr[pos];
        const size_t i_atom = (t_insr.x & 0x3ff);
        const size_t l_atom = (t_insr.y & 0x3ff);
        if (activity == VwuTask::INFR14 || activity == VwuTask::ALL_TASKS) {
          const int attn_idx = ((t_insr.y >> 10) & 0x1f);
          if (attn_idx > 0) {
            const Vec2<Tc> uc =
              evalAttenuated14Pair<llint,
                                   llint, Tc>(i_atom, l_atom, attn_idx, poly_vk.coulomb,
                                              sh_charges.data(), sh_lj_idx.data(),
                                              poly_vk.attn14_elec, poly_vk.attn14_vdw,
                                              poly_vk.lja_14_coeff, poly_vk.ljb_14_coeff,
                                              poly_vk.lj_14_sigma, ljabc_offset, nljt, xcrd_ptr,
                                              ycrd_ptr, zcrd_ptr, umat, invu, poly_psw->unit_cell,
                                              xfrc_ptr, yfrc_ptr, zfrc_ptr, eval_force, eval_force,
                                              poly_psw->inv_gpos_scale, poly_psw->frc_scale,
                                              clash_distance, clash_ratio, xcrd_ovrf_ptr,
                                              ycrd_ovrf_ptr, zcrd_ovrf_ptr, xfrc_ovrf_ptr,
                                              yfrc_ovrf_ptr, zfrc_ovrf_ptr);
            if (readBitFromMask(&poly_vk.cdhe_acc[nrg_limits.x], pos - cdhe_limits.x) == 1) {
              qq14_acc += llround(uc.x * nrg_scale_factor);
              lj14_acc += llround(uc.y * nrg_scale_factor);
            }
          }
          if (activity == VwuTask::INFR14) {
            continue;
          }
        }

        // Skip CHARMM improper interactions if only standard dihedrals are of interest, and
        // skip standard dihedrals if only CHARMM impropers are of interest.
        const bool is_charmm_improper = ((t_insr.x >> 30) & 0x1);
        if ((activity == VwuTask::DIHE && is_charmm_improper) ||
            (activity == VwuTask::CIMP && (! is_charmm_improper))) {
          continue;
        }
        const size_t j_atom = ((t_insr.x >> 10) & 0x3ff);
        const size_t k_atom = ((t_insr.x >> 20) & 0x3ff);
        const llint param_idx = ((t_insr.y >> 16) & 0xffff);
        const bool second_term = ((t_insr.y >> 15) & 0x1);
        const TorsionKind kind = (t_insr.x >> 31) ? TorsionKind::IMPROPER : TorsionKind::PROPER;
        DihedralStyle dihe_style;
        Tc ampl, phi, freq;
        if (is_charmm_improper) {
          dihe_style = DihedralStyle::HARMONIC;
          ampl = poly_vk.cimp_keq[param_idx];
          phi  = poly_vk.cimp_phi[param_idx];
          freq = value_one;
        }
        else {
          dihe_style = DihedralStyle::COSINE;
          ampl = poly_vk.dihe_amp[param_idx];
          phi  = poly_vk.dihe_phi[param_idx];
          freq = poly_vk.dihe_freq[param_idx];
        }
        Tc du = evalDihedralTwist<llint,
                                  llint, Tc>(i_atom, j_atom, k_atom, l_atom, ampl, phi, freq,
                                             dihe_style, xcrd_ptr, ycrd_ptr, zcrd_ptr, umat, invu,
                                             poly_psw->unit_cell, xfrc_ptr, yfrc_ptr, zfrc_ptr,
                                             eval_force, poly_psw->inv_gpos_scale,
                                             poly_psw->frc_scale, xcrd_ovrf_ptr, ycrd_ovrf_ptr,
                                             zcrd_ovrf_ptr, xfrc_ovrf_ptr, yfrc_ovrf_ptr,
                                             zfrc_ovrf_ptr);
        if (second_term) {

          // Any second term will, by construction, be a cosine-based dihedral.  Computing it in
          // this way foregoes the efficiency that could be gained by combining both terms into a
          // single twist applied to the same four atoms.  GPU kernels will exploit the
          // opportunity.  
          const uint xinsr = poly_vk.cdhe_ovrt_insr[pos];
          const int param_ii_idx = (xinsr & 0xfffff);
          du += evalDihedralTwist<llint,
                                  llint, Tc>(i_atom, j_atom, k_atom, l_atom,
                                             poly_vk.dihe_amp[param_ii_idx],
                                             poly_vk.dihe_phi[param_ii_idx],
                                             poly_vk.dihe_freq[param_ii_idx], dihe_style, xcrd_ptr,
                                             ycrd_ptr, zcrd_ptr, umat, invu, poly_psw->unit_cell,
                                             xfrc_ptr, yfrc_ptr, zfrc_ptr, eval_force,
                                             poly_psw->inv_gpos_scale, poly_psw->frc_scale,
                                             xcrd_ovrf_ptr, ycrd_ovrf_ptr, zcrd_ovrf_ptr,
                                             xfrc_ovrf_ptr, yfrc_ovrf_ptr, zfrc_ovrf_ptr);
        }
        if (readBitFromMask(&poly_vk.cdhe_acc[nrg_limits.x], pos - cdhe_limits.x) == 1) {
          if (is_charmm_improper) {
            cimp_acc += llround(du * nrg_scale_factor);
          }
          else {
            if (kind == TorsionKind::PROPER) {
              dihe_acc += llround(du * nrg_scale_factor);
            }
            else {
              impr_acc += llround(du * nrg_scale_factor);
            }
          }
        }
      }
    }

    // Evaluate CMAP interactions
    if (activity == VwuTask::CMAP || activity == VwuTask::ALL_TASKS) {
      const int2 cmap_limits = vwu_map[static_cast<int>(VwuAbstractMap::CMAP)];
      const int2 nrg_limits = vwu_map[static_cast<int>(VwuAbstractMap::CMAP_NRG)];
      for (int pos = cmap_limits.x; pos < cmap_limits.y; pos++) {
        const uint2 t_insr = poly_vk.cmap_insr[pos];
        const int i_atom = (t_insr.x & 0x3ff);
        const int j_atom = ((t_insr.x >> 10) & 0x3ff);
        const int k_atom = ((t_insr.x >> 20) & 0x3ff);
        const int l_atom = (t_insr.y & 0x3ff);
        const int m_atom = ((t_insr.y >> 10) & 0x3ff);
        const int surf_idx = (t_insr.y >> 20);
        const double du = evalCmap<llint,
                                   llint, Tc>(poly_vk.cmap_patches, poly_vk.cmap_patch_bounds,
                                              surf_idx, poly_vk.cmap_dim[surf_idx], i_atom, j_atom,
                                              k_atom, l_atom, m_atom, xcrd_ptr, ycrd_ptr, zcrd_ptr,
                                              umat, invu, poly_psw->unit_cell, xfrc_ptr,
                                              yfrc_ptr, zfrc_ptr, eval_force,
                                              poly_psw->inv_gpos_scale, poly_psw->frc_scale,
                                              xcrd_ovrf_ptr, ycrd_ovrf_ptr, zcrd_ovrf_ptr,
                                              xfrc_ovrf_ptr, yfrc_ovrf_ptr, zfrc_ovrf_ptr);
        if (readBitFromMask(&poly_vk.cmap_acc[nrg_limits.x], pos - cmap_limits.x) == 1) {
          cmap_acc += llround(du * nrg_scale_factor);
        }
      }
    }

    // Evaluate any remaining 1:4 interactions
    if (activity == VwuTask::INFR14 || activity == VwuTask::ALL_TASKS) {
      const int2 infr14_limits = vwu_map[static_cast<int>(VwuAbstractMap::INFR14)];
      const int2 nrg_limits = vwu_map[static_cast<int>(VwuAbstractMap::INFR14_NRG)];
      const int ljabc_offset = poly_vk.ljabc_offsets[sys_idx];
      const int nljt = poly_vk.n_lj_types[sys_idx];
      for (int pos = infr14_limits.x; pos < infr14_limits.y; pos++) {
        const uint t_insr = poly_vk.infr14_insr[pos];
        const int i_atom = (t_insr & 0x3ff);
        const int l_atom = ((t_insr >> 10) & 0x3ff);
        const int attn_idx = (t_insr >> 20);
        if (attn_idx == 0) {
          continue;
        }
        const Vec2<Tc> uc = evalAttenuated14Pair<llint,
                                                 llint, Tc>(i_atom, l_atom, attn_idx,
                                                            poly_vk.coulomb, sh_charges.data(),
                                                            sh_lj_idx.data(), poly_vk.attn14_elec,
                                                            poly_vk.attn14_vdw,
                                                            poly_vk.lja_14_coeff,
                                                            poly_vk.ljb_14_coeff,
                                                            poly_vk.lj_14_sigma, ljabc_offset,
                                                            nljt, xcrd_ptr, ycrd_ptr, zcrd_ptr,
                                                            umat, invu, poly_psw->unit_cell,
                                                            xfrc_ptr, yfrc_ptr, zfrc_ptr,
                                                            eval_force, eval_force,
                                                            poly_psw->inv_gpos_scale,
                                                            poly_psw->frc_scale, clash_distance,
                                                            clash_ratio, xcrd_ovrf_ptr,
                                                            ycrd_ovrf_ptr, zcrd_ovrf_ptr,
                                                            xfrc_ovrf_ptr, yfrc_ovrf_ptr,
                                                            zfrc_ovrf_ptr);
        if (readBitFromMask(&poly_vk.infr14_acc[nrg_limits.x], pos - infr14_limits.x) == 1) {
          qq14_acc += llround(uc.x * nrg_scale_factor);
          lj14_acc += llround(uc.y * nrg_scale_factor);
        }
      }
    }

    // The following restraint terms are only evaluated if the restraint abstract is nontrivial.
    if (poly_rk.rposn_step_bounds != nullptr) {
      
      // Evaluate positional restraints
      if (activity == VwuTask::RPOSN || activity == VwuTask::ALL_TASKS) {
        const int2 rposn_limits = vwu_map[static_cast<int>(VwuAbstractMap::RPOSN)];
        const int2 nrg_limits = vwu_map[static_cast<int>(VwuAbstractMap::RPOSN_NRG)];
        for (int pos = rposn_limits.x; pos < rposn_limits.y; pos++) {
          const uint2 t_insr = poly_rk.rposn_insr[pos];
          const int p_atom = (t_insr.x & 0x3ff);
          const int kr_param_idx = ((t_insr.x >> 10) & 0x1fffff);
          const int xyz_param_idx = t_insr.y;
          const Tc du = evalPosnRestraint<llint,
                                          llint, Tc>(p_atom, step,
                                                     poly_rk.rposn_step_bounds[kr_param_idx].x,
                                                     poly_rk.rposn_step_bounds[kr_param_idx].y,
                                                     poly_rk.rposn_init_xy[xyz_param_idx],
                                                     poly_rk.rposn_finl_xy[xyz_param_idx],
                                                     poly_rk.rposn_init_z[xyz_param_idx],
                                                     poly_rk.rposn_finl_z[xyz_param_idx],
                                                     poly_rk.rposn_init_k[kr_param_idx],
                                                     poly_rk.rposn_finl_k[kr_param_idx],
                                                     poly_rk.rposn_init_r[kr_param_idx],
                                                     poly_rk.rposn_finl_r[kr_param_idx],
                                                     xcrd_ptr, ycrd_ptr, zcrd_ptr, umat, invu,
                                                     poly_psw->unit_cell, xfrc_ptr, yfrc_ptr,
                                                     zfrc_ptr, eval_force,
                                                     poly_psw->inv_gpos_scale, poly_psw->frc_scale,
                                                     xcrd_ovrf_ptr, ycrd_ovrf_ptr, zcrd_ovrf_ptr,
                                                     xfrc_ovrf_ptr, yfrc_ovrf_ptr, zfrc_ovrf_ptr);
          if (readBitFromMask(&poly_rk.rposn_acc[nrg_limits.x], pos - rposn_limits.x) == 1) {
            rest_acc += llround(du * nrg_scale_factor);
          }
        }
      }

      // Evaluate distance restraints
      if (activity == VwuTask::RBOND || activity == VwuTask::ALL_TASKS) {
        const int2 rbond_limits = vwu_map[static_cast<int>(VwuAbstractMap::RBOND)];
        const int2 nrg_limits = vwu_map[static_cast<int>(VwuAbstractMap::RBOND_NRG)];
        for (int pos = rbond_limits.x; pos < rbond_limits.y; pos++) {
          const uint2 t_insr = poly_rk.rbond_insr[pos];
          const int i_atom = (t_insr.x & 0x3ff);
          const int j_atom = ((t_insr.x >> 10) & 0x3ff);
          const int param_idx = t_insr.y;
          const Tc du = evalBondRestraint<llint,
                                          llint,
                                          Tc, Tc2, Tc4>(i_atom, j_atom, step,
                                                        poly_rk.rbond_step_bounds[param_idx].x,
                                                        poly_rk.rbond_step_bounds[param_idx].y,
                                                        poly_rk.rbond_init_k[param_idx],
                                                        poly_rk.rbond_finl_k[param_idx],
                                                        poly_rk.rbond_init_r[param_idx],
                                                        poly_rk.rbond_finl_r[param_idx], xcrd_ptr,
                                                        ycrd_ptr, zcrd_ptr, umat, invu,
                                                        poly_psw->unit_cell, xfrc_ptr, yfrc_ptr,
                                                        zfrc_ptr, eval_force,
                                                        poly_psw->inv_gpos_scale,
                                                        poly_psw->frc_scale, xcrd_ovrf_ptr,
                                                        ycrd_ovrf_ptr, zcrd_ovrf_ptr,
                                                        xfrc_ovrf_ptr, yfrc_ovrf_ptr,
                                                        zfrc_ovrf_ptr);
          if (readBitFromMask(&poly_rk.rbond_acc[nrg_limits.x], pos - rbond_limits.x) == 1) {
            rest_acc += llround(du * nrg_scale_factor);
          }
        }
      }

      // Evaluate angle restraints
      if (activity == VwuTask::RANGL || activity == VwuTask::ALL_TASKS) {
        const int2 rangl_limits = vwu_map[static_cast<int>(VwuAbstractMap::RANGL)];
        const int2 nrg_limits = vwu_map[static_cast<int>(VwuAbstractMap::RANGL_NRG)];
        for (int pos = rangl_limits.x; pos < rangl_limits.y; pos++) {
          const uint2 t_insr = poly_rk.rangl_insr[pos];
          const int i_atom = (t_insr.x & 0x3ff);
          const int j_atom = ((t_insr.x >> 10) & 0x3ff);
          const int k_atom = ((t_insr.x >> 20) & 0x3ff);
          const int param_idx = t_insr.y;
          const Tc du = evalAnglRestraint<llint,
                                          llint,
                                          Tc, Tc2, Tc4>(i_atom, j_atom, k_atom, step,
                                                        poly_rk.rangl_step_bounds[param_idx].x,
                                                        poly_rk.rangl_step_bounds[param_idx].y,
                                                        poly_rk.rangl_init_k[param_idx],
                                                        poly_rk.rangl_finl_k[param_idx],
                                                        poly_rk.rangl_init_r[param_idx],
                                                        poly_rk.rangl_finl_r[param_idx], xcrd_ptr,
                                                        ycrd_ptr, zcrd_ptr, umat, invu,
                                                        poly_psw->unit_cell, xfrc_ptr, yfrc_ptr,
                                                        zfrc_ptr, eval_force,
                                                        poly_psw->inv_gpos_scale,
                                                        poly_psw->frc_scale, xcrd_ovrf_ptr,
                                                        ycrd_ovrf_ptr, zcrd_ovrf_ptr,
                                                        xfrc_ovrf_ptr, yfrc_ovrf_ptr,
                                                        zfrc_ovrf_ptr);
          if (readBitFromMask(&poly_rk.rangl_acc[nrg_limits.x], pos - rangl_limits.x) == 1) {
            rest_acc += llround(du * nrg_scale_factor);
          }
        }
      }

      // Evaluate dihedral restraints
      if (activity == VwuTask::RDIHE || activity == VwuTask::ALL_TASKS) {
        const int2 rdihe_limits = vwu_map[static_cast<int>(VwuAbstractMap::RDIHE)];
        const int2 nrg_limits = vwu_map[static_cast<int>(VwuAbstractMap::RDIHE_NRG)];
        for (int pos = rdihe_limits.x; pos < rdihe_limits.y; pos++) {
          const uint2 t_insr = poly_rk.rdihe_insr[pos];
          const int i_atom = (t_insr.x & 0x3ff);
          const int j_atom = ((t_insr.x >> 10) & 0x3ff);
          const int k_atom = ((t_insr.x >> 20) & 0x3ff);
          const int l_atom = (t_insr.y & 0x3ff);
          const int param_idx = ((t_insr.y >> 10) & 0x3fffff);
          const Tc du = evalDiheRestraint<llint,
                                          llint,
                                          Tc, Tc2, Tc4>(i_atom, j_atom, k_atom, l_atom, step,
                                                        poly_rk.rdihe_step_bounds[param_idx].x,
                                                        poly_rk.rdihe_step_bounds[param_idx].y,
                                                        poly_rk.rdihe_init_k[param_idx],
                                                        poly_rk.rdihe_finl_k[param_idx],
                                                        poly_rk.rdihe_init_r[param_idx],
                                                        poly_rk.rdihe_finl_r[param_idx],
                                                        xcrd_ptr, ycrd_ptr, zcrd_ptr, umat, invu,
                                                        poly_psw->unit_cell, xfrc_ptr, yfrc_ptr,
                                                        zfrc_ptr, eval_force,
                                                        poly_psw->inv_gpos_scale,
                                                        poly_psw->frc_scale, xcrd_ovrf_ptr,
                                                        ycrd_ovrf_ptr, zcrd_ovrf_ptr,
                                                        xfrc_ovrf_ptr, yfrc_ovrf_ptr,
                                                        zfrc_ovrf_ptr);
          if (readBitFromMask(&poly_rk.rdihe_acc[nrg_limits.x], pos - rdihe_limits.x) == 1) {
            rest_acc += llround(du * nrg_scale_factor);
          }
        }
      }
    }

    // Write results back to the main force arrays.  In serial processing mode all forces can be
    // contributed back without risk of race conditions.
    switch (eval_force) {
    case EvaluateForce::NO:
      break;
    case EvaluateForce::YES:
      {
        const int2 manip_limits = vwu_map[static_cast<size_t>(VwuAbstractMap::MANIPULATE)];
        for (int j = import_limits.x; j < import_limits.y; j++) {
          const size_t synth_atom = poly_vk.vwu_imports[j];
          const size_t local_atom = j - import_limits.x;
          const int manip_segment = (local_atom >> 5) + manip_limits.x;
          const int manip_bitpos = (local_atom & 0x1f);
          const uint2 manip_mask = poly_auk.vwu_manip[manip_segment];
          if ((manip_mask.y >> manip_bitpos) & 0x1) {
            if (poly_psw->frc_bits > force_scale_nonoverflow_bits) {
              const int95_t nfx = { sh_xfrc[local_atom], sh_xfrc_ovrf[local_atom] };
              const int95_t nfy = { sh_yfrc[local_atom], sh_yfrc_ovrf[local_atom] };
              const int95_t nfz = { sh_zfrc[local_atom], sh_zfrc_ovrf[local_atom] };
              hostSplitFPSum(&poly_psw->xfrc[synth_atom], &poly_psw->xfrc_ovrf[synth_atom], nfx);
              hostSplitFPSum(&poly_psw->yfrc[synth_atom], &poly_psw->yfrc_ovrf[synth_atom], nfy);
              hostSplitFPSum(&poly_psw->zfrc[synth_atom], &poly_psw->zfrc_ovrf[synth_atom], nfz);
            }
            else {
              poly_psw->xfrc[synth_atom] += sh_xfrc[local_atom];
              poly_psw->yfrc[synth_atom] += sh_yfrc[local_atom];
              poly_psw->zfrc[synth_atom] += sh_zfrc[local_atom];
            }
          }
        }
      }
      break;
    }

    // Contribute results to the instantaneous states.
    commitVwuEnergies(bond_acc, angl_acc, dihe_acc, impr_acc, ubrd_acc, cimp_acc, cmap_acc,
                      qq14_acc, lj14_acc, rest_acc, sys_idx, activity, sc);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
void evalValeRestMM(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                    const double* invu, const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                    Tforce* zfrc, ScoreCard *sc, const ValenceKit<Tcalc> &vk,
                    const NonbondedKit<Tcalc> &nbk, const RestraintKit<Tcalc, Tcalc2, Tcalc4> &rar,
                    const EvaluateForce eval_force, const int step, const int system_index,
                    const Tcalc inv_gpos_factor, const Tcalc force_factor,
                    const Tcalc clash_distance, const Tcalc clash_ratio) {
  evalValeMM<Tcoord, Tforce, Tcalc>(xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                                    vk, nbk, eval_force, system_index, inv_gpos_factor,
                                    force_factor, clash_distance, clash_ratio);
  evaluateRestraints<Tcoord, Tforce,
                     Tcalc, Tcalc2, Tcalc4>(rar, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                            yfrc, zfrc, sc, eval_force, system_index, step,
                                            inv_gpos_factor, force_factor);
}

//-----------------------------------------------------------------------------------------------
template <typename Tc, typename Tc2, typename Tc4>
void evalValeRestMM(PsSynthesisWriter *poly_psw, ScoreCard* sc, const SyValenceKit<Tc> &poly_vk,
                    const SyRestraintKit<Tc, Tc2, Tc4> &poly_rk,
                    const SyAtomUpdateKit<Tc, Tc2, Tc4> &poly_auk,
                    const EvaluateForce eval_force, const VwuTask activity, const int step,
                    const Tc clash_distance, const Tc  clash_ratio) {
  evalValeMM<Tc, Tc2, Tc4>(poly_psw, sc, poly_vk, poly_auk, eval_force, activity,
                           clash_distance, clash_ratio, step, poly_rk);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
void evalNonbValeMM(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                    const double* invu, const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                    Tforce* zfrc, ScoreCard *sc, const ValenceKit<Tcalc> &vk,
                    const NonbondedKit<Tcalc> &nbk, const StaticExclusionMaskReader &ser,
                    const EvaluateForce eval_force, const int system_index,
                    const Tcalc inv_gpos_factor, const Tcalc force_factor,
                    const Tcalc clash_distance, const Tcalc clash_ratio) {
  evaluateNonbondedEnergy<Tcoord, Tforce, Tcalc>(nbk, ser, xcrd, ycrd, zcrd, umat, invu, unit_cell,
                                                 xfrc, yfrc, zfrc, sc, eval_force, eval_force,
                                                 system_index, inv_gpos_factor, force_factor,
                                                 clash_distance, clash_ratio);
  evalValeMM<Tcoord, Tforce, Tcalc>(xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                                    vk, nbk, eval_force, system_index, inv_gpos_factor,
                                    force_factor, clash_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
void evalNonbValeRestMM(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                        const double* umat, const double* invu, const UnitCellType unit_cell,
                        Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, ScoreCard *sc,
                        const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                        const StaticExclusionMaskReader &ser,
                        const RestraintKit<Tcalc, Tcalc2, Tcalc4> &rar,
                        const EvaluateForce eval_force, const int system_index, const int step,
                        const Tcalc inv_gpos_factor, const Tcalc force_factor,
                        const Tcalc clash_distance, const Tcalc clash_ratio) {
  evaluateNonbondedEnergy<Tcoord, Tforce, Tcalc>(nbk, ser, xcrd, ycrd, zcrd, umat, invu, unit_cell,
                                                 xfrc, yfrc, zfrc, sc, eval_force, eval_force,
                                                 system_index, inv_gpos_factor, force_factor,
                                                 clash_distance, clash_ratio);
  evalValeMM<Tcoord, Tforce, Tcalc>(xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc,
                                    vk, nbk, eval_force, system_index, inv_gpos_factor,
                                    force_factor, clash_distance, clash_ratio);
  evaluateRestraints<Tcoord, Tforce,
                     Tcalc, Tcalc2, Tcalc4>(rar, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                            yfrc, zfrc, sc, eval_force, system_index, step,
                                            inv_gpos_factor, force_factor);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
void evalRestrainedMMGB(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                        const double* umat, const double* invu, const UnitCellType unit_cell,
                        Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, ScoreCard *sc,
                        const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                        const StaticExclusionMaskReader &ser, const ImplicitSolventKit<Tcalc> &isk,
                        const NeckGeneralizedBornKit<Tcalc> &neck_gbk,
                        Tforce* effective_gb_radii, Tforce *psi, Tforce *sumdeijda,
                        const RestraintKit<Tcalc, Tcalc2, Tcalc4> &rar,
                        const EvaluateForce eval_force, const int system_index, const int step,
                        const Tcalc inv_gpos_factor, const Tcalc force_factor,
                        const Tcalc clash_distance, const Tcalc clash_ratio) {
  evalNonbValeRestMM(xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, sc, vk, nbk, ser,
                     rar, eval_force, system_index, step, inv_gpos_factor, force_factor,
                     clash_distance, clash_ratio);
  evaluateGeneralizedBornEnergy<Tcoord, Tforce, Tcalc>(nbk, ser, isk, neck_gbk, xcrd, ycrd, zcrd,
                                                       xfrc, yfrc, zfrc, effective_gb_radii, psi,
                                                       sumdeijda, sc, eval_force, system_index,
                                                       inv_gpos_factor, force_factor);
}

} // namespace mm
} // namespace stormm
