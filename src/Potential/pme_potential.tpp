// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
double2 cellToCellInteractions(PhaseSpaceWriter *psw, const std::vector<int> &cell_list,
                               const std::vector<int> &cell_bounds, const int ncell_a,
                               const int ncell_b, const int ncell_c, const int cell_idx_a,
                               const int cell_idx_b, const int cell_idx_c, const int i_pair,
                               const int j_pair, const int k_pair, const NonbondedKit<Tcalc> &nbk,
                               const LocalExclusionMaskReader &lemr, const Tcalc elec_cutoff,
                               const Tcalc vdw_cutoff, const Tcalc qqew_coeff,
                               const Tcalc ljew_cutoff, const VdwSumMethod vdw_sum,
                               const EvaluateForce eval_frc, const NonbondedTheme theme) {

  // Initialize the energy (always returned) and determine the nature of the calculation.
  const bool tcalc_is_double = std::type_index(typeid(Tcalc)).hash_code();
  const Tcalc value_zero = 0.0;
  const Tcalc value_one  = 1.0;
  const Tcalc qqcut_sq = elec_cutoff * elec_cutoff;
  const Tcalc ljcut_sq = vdw_cutoff * vdw_cutoff;
  const Tcalc sqrt_pi = sqrt(symbols::pi);
  const Tcalc qq_bfac = 2.0 * qqew_coeff / sqrt_pi;
  double2 result = { 0.0, 0.0 };
  const int slab_ab = ncell_a * ncell_b;
  int ipca = cell_idx_a + i_pair;
  int ipcb = cell_idx_b + j_pair;
  int ipcc = cell_idx_c + k_pair;
  ipca += ((ipca < 0) - (ipca >= ncell_a)) * ncell_a;
  ipcb += ((ipcb < 0) - (ipcb >= ncell_b)) * ncell_b;
  ipcc += ((ipcc < 0) - (ipcc >= ncell_c)) * ncell_c;
  const int cell_idx = (((cell_idx_c * ncell_b) + cell_idx_b) * ncell_a) + cell_idx_a;
  const int pair_idx = (((ipcc * ncell_b) + ipcb) * ncell_a) + ipca;
  const int illim = cell_bounds[cell_idx];
  const int ihlim = cell_bounds[cell_idx + 1];
  const int jmax = cell_bounds[pair_idx + 1] - cell_bounds[pair_idx];
  std::vector<Tcalc> dx(jmax), dy(jmax), dz(jmax);
  for (int i = illim; i < ihlim; i++) {
    const int iatom = cell_list[i];
    const int jllim = cell_bounds[pair_idx];
    const int jhlim = (pair_idx == cell_idx) ? i : cell_bounds[pair_idx + 1];
    const int nj = jhlim - jllim;
    for (int j = 0; j < nj; j++) {
      const int jatom = cell_list[jllim + j];
      dx[j] = psw->xcrd[jatom] - psw->xcrd[iatom];
      dy[j] = psw->ycrd[jatom] - psw->ycrd[iatom];
      dz[j] = psw->zcrd[jatom] - psw->zcrd[iatom];
    }
    const Tcalc scaled_qi = nbk.coulomb_constant * nbk.charge[iatom];
    imageCoordinates<Tcalc, Tcalc>(&dx, &dy, &dz, psw->umat, psw->invu, psw->unit_cell,
                                   ImagingMethod::MINIMUM_IMAGE, 1.0);
    for (int j = 0; j < nj; j++) {
      const int jatom = cell_list[jllim + j];
      const Tcalc r2 = (dx[j] * dx[j]) + (dy[j] * dy[j]) + (dz[j] * dz[j]);
      const Tcalc invr2 = value_one / r2;
      const Tcalc invr = (tcalc_is_double) ? sqrt(invr2) : sqrtf(invr2);
      Tcalc fmag = value_zero;
      if (testExclusion(lemr, iatom, jatom)) {

        // Subtract the primary image term for electrostatic and / or van-der Waals potentials.
        // For van-der Waals potentials involving any sort of cutoff, it is assumed that any such
        // exclusions will take place within the cutoff where the full potential is in effect.
        switch (theme) {
        case NonbondedTheme::ELECTROSTATIC:
        case NonbondedTheme::ALL:
          {
            const Tcalc kqq = scaled_qi * nbk.charge[jatom];
            const Tcalc u_quant = erfc(qqew_coeff / invr) * invr;
            result.x += kqq * (u_quant - invr);
            if (eval_frc == EvaluateForce::YES) {
              const Tcalc exp_quant = qq_bfac * exp(-qqew_coeff * qqew_coeff * r2);
              fmag -= kqq * (((exp_quant + u_quant) * invr) - invr2) * invr;
            }
          }
          break;
        case NonbondedTheme::VAN_DER_WAALS:
          break;
        }
        switch (theme) {
        case NonbondedTheme::VAN_DER_WAALS:
        case NonbondedTheme::ALL:
          switch (vdw_sum) {
          case VdwSumMethod::PME:
            {
              const int atyp_i = nbk.lj_idx[iatom];
              const int atyp_j = nbk.lj_idx[jatom];
              const size_t atyp_ij = (nbk.n_lj_types * atyp_i) + atyp_j;
              const Tcalc invr4 = invr2 * invr2;
              const Tcalc invr6 = invr4 * invr2;
              result.y += nbk.ljb_coeff[atyp_ij] * invr6;
              if (eval_frc == EvaluateForce::YES) {
                fmag += ((nbk.lja_coeff[atyp_ij] * invr6) - nbk.ljb_coeff[atyp_ij]) * invr4 *
                        invr4;
              }
            }
            break;
          case VdwSumMethod::CUTOFF:
          case VdwSumMethod::SMOOTH:
            break;
          }
          break;
        case NonbondedTheme::ELECTROSTATIC:
          break;
        }
      }
      else {

        // Evaluate the electrostatic component of the potential
        switch (theme) {
        case NonbondedTheme::ELECTROSTATIC:
        case NonbondedTheme::ALL:
          if (r2 < qqcut_sq) {
            Tcalc u_quant, exp_quant;
            const Tcalc kqq = scaled_qi * nbk.charge[jatom];
            if (tcalc_is_double) {
              const Tcalc r = sqrt(r2);
              u_quant = erfc(qqew_coeff * r) * invr;
              result.x += kqq * u_quant;
              if (eval_frc == EvaluateForce::YES) {
                exp_quant = qq_bfac * exp(-qqew_coeff * qqew_coeff * r2);
                fmag = -kqq * (exp_quant + u_quant) * invr2;
              }
            }
            else {
              const Tcalc r = sqrtf(r2);
              u_quant = erfcf(qqew_coeff * r) * invr;
              result.x += kqq * u_quant;
              if (eval_frc == EvaluateForce::YES) {
                exp_quant = qq_bfac * expf(-qqew_coeff * qqew_coeff * r2);
                fmag = -kqq * (exp_quant + u_quant) * invr2;
              }
            }
          }
          else {
            fmag = value_zero;
          }
          break;
        case NonbondedTheme::VAN_DER_WAALS:
          fmag = value_zero;
          break;
        }

        // Evaluate the van-der Waals component of the potential
        switch (theme) {
        case NonbondedTheme::VAN_DER_WAALS:
        case NonbondedTheme::ALL:
          if (r2 < ljcut_sq) {
            const int atyp_i = nbk.lj_idx[iatom];
            const int atyp_j = nbk.lj_idx[jatom];
            const size_t atyp_ij = (nbk.n_lj_types * atyp_i) + atyp_j;
            switch (vdw_sum) {
            case VdwSumMethod::CUTOFF:
              {
                const Tcalc invr4 = invr2 * invr2;
                const Tcalc invr6 = invr4 * invr2;
                result.y += ((nbk.lja_coeff[atyp_ij] * invr6) - nbk.ljb_coeff[atyp_ij]) * invr6;
                if (eval_frc == EvaluateForce::YES) {
                  const Tcalc invr8 = invr4 * invr4;
                  fmag -= ((12.0 * nbk.lja_coeff[atyp_ij] * invr6) -
                           (6.0 * nbk.ljb_coeff[atyp_ij])) * invr8;
                }
              }
              break;
            case VdwSumMethod::SMOOTH:
            case VdwSumMethod::PME:
              break;
            }
          }
          break;
        case NonbondedTheme::ELECTROSTATIC:
          break;
        }
      }
      
      // Distribute the accumulated force across the pair of interacting particles
      if (eval_frc == EvaluateForce::YES && std::abs(fmag) > value_zero) {
        const Tcalc fx = fmag * dx[j];
        const Tcalc fy = fmag * dy[j];
        const Tcalc fz = fmag * dz[j];
        psw->xfrc[iatom] += fx;
        psw->yfrc[iatom] += fy;
        psw->zfrc[iatom] += fz;
        psw->xfrc[jatom] -= fx;
        psw->yfrc[jatom] -= fy;
        psw->zfrc[jatom] -= fz;
      }
    }
  }

  // Return the accumulated energy as a tuple of double-precision values
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
double2 evaluateParticleParticleEnergy(PhaseSpaceWriter *psw, const NonbondedKit<Tcalc> &nbk,
                                       const LocalExclusionMaskReader &lemr,
                                       const Tcalc elec_cutoff, const Tcalc vdw_cutoff,
                                       const Tcalc qqew_coeff, const Tcalc ljew_coeff,
                                       const VdwSumMethod vdw_sum, const EvaluateForce eval_frc,
                                       const NonbondedTheme theme) {
  
  // Compute the Ewald coefficient, if this has not already been done.
  const double actual_qqew_coeff = (qqew_coeff < 1.0e-6) ?
                                   ewaldCoefficient(elec_cutoff, default_dsum_tol) : qqew_coeff;
  const double actual_ljew_coeff = (ljew_coeff < 1.0e-6) ?
                                   ewaldCoefficient(vdw_cutoff, default_dsum_tol) : ljew_coeff;
  
  // Construct a neighbor list internally, using a method separate from what would be used by the
  // more complex cell grid.
  const double max_cutoff = std::max(elec_cutoff, vdw_cutoff);
  std::vector<int> cell_residency(psw->natom);
  double norm_abl, norm_acl, norm_bcl;
  hessianNormalWidths<double>(psw->invu, &norm_bcl, &norm_acl, &norm_abl);
  const int ncell_a = floor(norm_bcl / max_cutoff);
  const int ncell_b = floor(norm_acl / max_cutoff);
  const int ncell_c = floor(norm_abl / max_cutoff);
  const double dncell_a = ncell_a;
  const double dncell_b = ncell_b;
  const double dncell_c = ncell_c;
  const double fcell_a = 1.0 / static_cast<double>(ncell_a);
  const double fcell_b = 1.0 / static_cast<double>(ncell_b);
  const double fcell_c = 1.0 / static_cast<double>(ncell_c);
  const int natom = psw->natom;
  const double* xcrd = psw->xcrd;
  const double* ycrd = psw->ycrd;
  const double* zcrd = psw->zcrd;
  const double* umat = psw->umat;
  for (int i = 0; i < natom; i++) {
    double fx = (umat[0] * xcrd[i]) + (umat[3] * ycrd[i]) + (umat[6] * zcrd[i]);
    double fy = (umat[1] * xcrd[i]) + (umat[4] * ycrd[i]) + (umat[7] * zcrd[i]);
    double fz = (umat[2] * xcrd[i]) + (umat[5] * ycrd[i]) + (umat[8] * zcrd[i]);
    fx -= floor(fx);
    fy -= floor(fy);
    fz -= floor(fz);
    const int cidx_a  = (fx * dncell_a);
    const int cidx_b  = (fy * dncell_b);
    const int cidx_c  = (fz * dncell_c);
    cell_residency[i] = (((cidx_c * ncell_b) + cidx_b) * ncell_a) + cidx_a;
  }
  const int ncell = ncell_a * ncell_b * ncell_c;
  std::vector<int> cell_bounds(ncell + 1);
  std::vector<int> cell_list(natom);
  indexingArray(cell_residency, &cell_list, &cell_bounds);
  double2 result = { 0.0, 0.0 };
  const int slab_ab = ncell_a * ncell_b;
  for (int cidx = 0; cidx < ncell; cidx++) {
    const int cidx_c = cidx / (slab_ab);
    const int cidx_b = (cidx - (cidx_c * slab_ab)) / ncell_a;
    const int cidx_a = cidx - (cidx_c * slab_ab) - (cidx_b * ncell_a);
    
    // Loop over all cells in the half-bowl, computing non-bonded interactions according to the
    // interaction theme.
    const int imin = (ncell_a > 2) ? -1 : -(cidx_a == 0);
    const int jmin = (ncell_b > 2) ? -1 : -(cidx_b == 0);
    const int kmin = (ncell_c > 2) ? -1 : -(cidx_c == 0);
    const int imax = (ncell_a > 2) ? 1 : 1 - (cidx_a == 0);
    const int jmax = (ncell_b > 2) ? 1 : 1 - (cidx_b == 0);
    for (int i = imin; i <= imax; i++) {
      for (int j = jmin; j <= jmax; j++) {
        for (int k = kmin; k <= 0; k++) {
          if (k < 0 || j < 0 || (j == 0 && i <= 0)) {
            const double2 tmp_e = cellToCellInteractions<Tcalc>(psw, cell_list, cell_bounds,
                                                                ncell_a, ncell_b, ncell_c,
                                                                cidx_a, cidx_b, cidx_c, i, j, k,
                                                                nbk, lemr, elec_cutoff, vdw_cutoff,
                                                                actual_qqew_coeff,
                                                                actual_ljew_coeff, vdw_sum,
                                                                eval_frc, theme);
            result.x += tmp_e.x;
            result.y += tmp_e.y;
          }
        }
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcalc2, typename Tcalc4>
double2 basicTileInteractions(const std::vector<Tcalc> &a_xpos, const std::vector<Tcalc> &a_ypos,
                              const std::vector<Tcalc> &a_zpos, const std::vector<Tcalc> &b_xpos,
                              const std::vector<Tcalc> &b_ypos, const std::vector<Tcalc> &b_zpos,
                              const std::vector<Tcalc> &scl_aq, const std::vector<Tcalc> &bq,
                              const std::vector<int> &ofs_aljidx, const std::vector<int> &bljidx,
                              const std::vector<int> &top_aidx, const std::vector<int> &top_bidx,
                              const std::vector<uint> &img_aidx, const std::vector<uint> &img_bidx,
                              const int system_index, const int na, const int nb,
                              const bool self_interaction, const bool small_box,
                              CellGridWriter<Tcoord, Tacc, Tcalc, Tcalc4> *cgw,
                              const PsSynthesisReader &poly_psr,
                              const SyNonbondedKit<Tcalc, Tcalc2> &poly_nbk,
                              const LocalExclusionMaskReader &lemr, const Tcalc elec_cutsq,
                              const Tcalc vdw_cutsq, const Tcalc qqew_coeff,
                              const Tcalc ljew_coeff, const VdwSumMethod vdw_sum,
                              const EvaluateForce eval_frc, const NonbondedTheme theme) {
  
  // Some useful constants
  const bool tcoord_is_real = isFloatingPointScalarType<Tcoord>();
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const Tcalc value_zero = 0.0;
  const Tcalc value_one  = 1.0;
  const Tcalc qq_bfac = 2.0 * qqew_coeff / sqrt(symbols::pi);
  const Tcalc frc_scale = poly_psr.frc_scale;
  double2 result = { 0.0, 0.0 };
  const int xfrm_offset = system_index * roundUp(9, warp_size_int);
  for (int i = 0; i < na; i++) {
    const int jlim = (self_interaction) ? i : nb;
    for (int j = 0; j < jlim; j++) {
      Tcalc dx = b_xpos[j] - a_xpos[i];
      Tcalc dy = b_ypos[j] - a_ypos[i];
      Tcalc dz = b_zpos[j] - a_zpos[i];
      if (small_box) {

        // When the box is less than five cells in any particular direction, it becomes important
        // to re-image all interactions.
        imageCoordinates<Tcalc, Tcalc>(&dx, &dy, &dz, &poly_psr.umat[xfrm_offset],
                                       &poly_psr.invu[xfrm_offset], poly_psr.unit_cell,
                                       ImagingMethod::MINIMUM_IMAGE, 1.0);
      }
      const Tcalc r2 = (dx * dx) + (dy * dy) + (dz * dz);
      const Tcalc invr2 = value_one / r2;
      const Tcalc invr = (tcalc_is_double) ? sqrt(invr2) : sqrtf(invr2);
      Tcalc fmag = value_zero;
      if (testExclusion(lemr, top_aidx[i], top_bidx[j])) {

        // Compute the electrostatic excluded interaction
        switch (theme) {
        case NonbondedTheme::ELECTROSTATIC:
        case NonbondedTheme::ALL:
          {
            const Tcalc kqq = scl_aq[i] * bq[j];
            const Tcalc u_quant = erfc(qqew_coeff / invr) * invr;
            result.x += kqq * (u_quant - invr);
            if (eval_frc == EvaluateForce::YES) {
              const Tcalc exp_quant = qq_bfac * exp(-qqew_coeff * qqew_coeff * r2);
              fmag -= kqq * (((exp_quant + u_quant) * invr) - invr2) * invr;
            }
          }
          break;
        case NonbondedTheme::VAN_DER_WAALS:
          break;
        }
        switch (theme) {
        case NonbondedTheme::VAN_DER_WAALS:
        case NonbondedTheme::ALL:

          // Compute the excluded van-der Waals interaction
          switch (vdw_sum) {
          case VdwSumMethod::PME:
            {
              const size_t atyp_ij = ofs_aljidx[i] + bljidx[j];
              const Tcalc invr4 = invr2 * invr2;
              const Tcalc invr6 = invr4 * invr2;
              const Tcalc2 tlj_ab = poly_nbk.ljab_coeff[atyp_ij];
              result.y += tlj_ab.y * invr6;
              if (eval_frc == EvaluateForce::YES) {
                fmag += ((12.0 * (tlj_ab.x * invr6)) - (6.0 * tlj_ab.y)) * invr4 * invr4;
              }
            }
            break;
          case VdwSumMethod::CUTOFF:
          case VdwSumMethod::SMOOTH:
            break;
          }
          break;
        case NonbondedTheme::ELECTROSTATIC:
          break;
        }
      }
      else {

        // Compute the non-excluded electrostatic potential
        switch (theme) {
        case NonbondedTheme::ELECTROSTATIC:
        case NonbondedTheme::ALL:
          if (r2 < elec_cutsq) {
            Tcalc u_quant, exp_quant;
            const Tcalc kqq = scl_aq[i] * bq[j];
            if (tcalc_is_double) {
              u_quant = erfc(qqew_coeff / invr) * invr;
              result.x += kqq * u_quant;
              if (eval_frc == EvaluateForce::YES) {
                exp_quant = qq_bfac * exp(-qqew_coeff * qqew_coeff * r2);
                fmag -= kqq * (exp_quant + u_quant) * invr2;
              }
            }
            else {
              u_quant = erfcf(qqew_coeff / invr) * invr;
              result.x += kqq * u_quant;
              if (eval_frc == EvaluateForce::YES) {
                exp_quant = qq_bfac * expf(-qqew_coeff * qqew_coeff * r2);
                fmag -= kqq * (exp_quant + u_quant) * invr2;
              }
            }
          }
          break;
        case NonbondedTheme::VAN_DER_WAALS:
          break;
        }

        // Compute the non-excluded van-der Waals potential
        switch (theme) {
        case NonbondedTheme::VAN_DER_WAALS:
        case NonbondedTheme::ALL:
          if (r2 < vdw_cutsq) {
            const size_t atyp_ij = ofs_aljidx[i] + bljidx[j];
            switch (vdw_sum) {
            case VdwSumMethod::CUTOFF:
              {
                const Tcalc invr4 = invr2 * invr2;
                const Tcalc invr6 = invr4 * invr2;
                const Tcalc2 tlj_ab = poly_nbk.ljab_coeff[atyp_ij];
                result.y += ((tlj_ab.x * invr6) - tlj_ab.y) * invr6;
                if (eval_frc == EvaluateForce::YES) {
                  fmag -= ((12.0 * tlj_ab.x * invr6) - (6.0 * tlj_ab.y)) * invr4 * invr4;
                }
              }
              break;
            case VdwSumMethod::SMOOTH:
            case VdwSumMethod::PME:
              break;
            }
          }
          break;
        case NonbondedTheme::ELECTROSTATIC:
          break;
        }
      }

      // Distribute the accumulated force across the pair of interacting particles
      const size_t a_atom_idx = img_aidx[i];
      const size_t b_atom_idx = img_bidx[j];
      if (eval_frc == EvaluateForce::YES && std::abs(fmag) > value_zero) {
        fmag *= frc_scale;
        const Tcalc fx = fmag * dx;
        const Tcalc fy = fmag * dy;
        const Tcalc fz = fmag * dz;
        hostSplitAccumulation(fx, &cgw->xfrc[a_atom_idx], &cgw->xfrc_ovrf[a_atom_idx]);
        hostSplitAccumulation(fy, &cgw->yfrc[a_atom_idx], &cgw->yfrc_ovrf[a_atom_idx]);
        hostSplitAccumulation(fz, &cgw->zfrc[a_atom_idx], &cgw->zfrc_ovrf[a_atom_idx]);
        hostSplitAccumulation(-fx, &cgw->xfrc[b_atom_idx], &cgw->xfrc_ovrf[b_atom_idx]);
        hostSplitAccumulation(-fy, &cgw->yfrc[b_atom_idx], &cgw->yfrc_ovrf[b_atom_idx]);
        hostSplitAccumulation(-fz, &cgw->zfrc[b_atom_idx], &cgw->zfrc_ovrf[b_atom_idx]);
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcalc2, typename Tcoord4>
double2 towerPlatePairInteractions(CellGridWriter<Tcoord, Tacc, Tcalc, Tcoord4> *cgw,
                                   const PsSynthesisReader &poly_psr,
                                   const int system_idx, const int cell_idx,
                                   const SyNonbondedKit<Tcalc, Tcalc2> &poly_nbk,
                                   const LocalExclusionMaskReader &lemr,
                                   const Tcalc elec_cutoff, const Tcalc vdw_cutoff,
                                   const Tcalc qqew_coeff, const Tcalc ljew_coeff,
                                   const VdwSumMethod vdw_sum, const EvaluateForce eval_frc,
                                   const NonbondedTheme theme) {

  // Some useful constants
  const bool tcoord_is_real = isFloatingPointScalarType<Tcoord>();
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const uint2 two_zeros = { 0U, 0U };
  const Tcalc value_zero = 0.0;
  const Tcalc value_one  = 1.0;
  const Tcalc elec_cutsq = elec_cutoff * elec_cutoff;
  const Tcalc vdw_cutsq = vdw_cutoff * vdw_cutoff;
  const Tcalc qq_bfac = 2.0 * qqew_coeff / sqrt(symbols::pi);
  const int sys_ljidx_offset = poly_nbk.ljabc_offsets[system_idx];
  const int nsys_lj_types = poly_nbk.n_lj_types[system_idx];

  // The tower is organized along the unit cell C axis, to maximize the coalescence of data along
  // the A axis for the larger number of cells in the plate.  Develop a list of the atom limits
  // for the tower and plate.  The final task is to calculate the interactions of the plate's
  // central cell (the third of five in the tower) with themselves and the first two batches of
  // atoms in the tower.
  const ullint footprint = cgw->system_cell_grids[system_idx];
  const int cell_start = (footprint & 0xfffffffLLU);
  const int ncell_a = ((footprint >> 28) & 0xfffLLU);
  const int ncell_b = ((footprint >> 40) & 0xfffLLU);
  const int ncell_c = (footprint >> 52);
  const bool small_box = (ncell_a < 5 || ncell_b < 5 || ncell_c < 5);
  const int cell_count = ncell_a * ncell_b * ncell_c;
  const int cell_cidx = cell_idx / (ncell_a * ncell_b);
  const int cell_bidx = (cell_idx - (cell_cidx * ncell_a * ncell_b)) / ncell_a;
  const int cell_aidx = cell_idx - (((cell_cidx * ncell_b) + cell_bidx) * ncell_a);
  std::vector<uint2> tower_atoms(5, two_zeros), plate_atoms(12, two_zeros);
  const int tower_kmax = 1 + (ncell_c > 4);
  for (int i = -2; i <= tower_kmax; i++) {
    int tower_cidx = cell_cidx + i;
    tower_cidx += ((tower_cidx < 0) - (tower_cidx >= ncell_c)) * ncell_c;
    const int celli_idx = cell_start + (((tower_cidx * ncell_b) + cell_bidx) * ncell_a) +
                          cell_aidx;
    tower_atoms[i + 2] = cgw->cell_limits[celli_idx];
  }
  int plate_imin, plate_imax;
  if (ncell_a > 4) {
    plate_imin = -2;
    plate_imax = 2;
  }
  else {
    if (cell_aidx < 2) {
      plate_imin = -2;
      plate_imax = 1;
    }
    else {
      plate_imin = -1;
      plate_imax = 2;
    }
  }
  const int plate_jmin = -2 + (ncell_b == 4 && cell_bidx >= 2);
  const int plate_jmax = 0;
  for (int i = plate_imin; i <= plate_imax; i++) {
    for (int j = plate_jmin; j <= plate_jmax; j++) {
      if (i < 0 || j < 0) {
        int plate_aidx = cell_aidx + i;
        plate_aidx += ((plate_aidx < 0) - (plate_aidx >= ncell_a)) * ncell_a;
        int plate_bidx = cell_bidx + j;
        plate_bidx += (plate_bidx < 0) * ncell_b;
        const int cellij_idx = cell_start + (((cell_cidx * ncell_b) + plate_bidx) * ncell_a) +
                               plate_aidx;
        plate_atoms[(5 * (j + 2)) + i + 2] = cgw->cell_limits[cellij_idx];
      }
    }
  }
  
  // Compute the total numbers of atoms in the tower and in the plate.
  std::vector<int> tower_prefix(6, 0);
  std::vector<int> plate_prefix(13, 0);
  uint ntower = 0;
  for (size_t i = 0; i < 5; i++) {
    ntower += (tower_atoms[i].y >> 16);
    tower_prefix[i + 1] = ntower;
  }
  uint nplate = 0;
  for (size_t i = 0; i < 12; i++) {
    nplate += (plate_atoms[i].y >> 16);
    plate_prefix[i + 1] = nplate;
  }
  
  // Loop over all interactions between the tower and plate, respecting exclusions.
  double2 result = { 0.0, 0.0 };
  int tw_base = 0;
  std::vector<Tcalc> tw_xpos(warp_size_int), tw_ypos(warp_size_int), tw_zpos(warp_size_int);
  std::vector<Tcalc> pl_xpos(warp_size_int), pl_ypos(warp_size_int), pl_zpos(warp_size_int);
  std::vector<uint> tower_img_idx(warp_size_int), plate_img_idx(warp_size_int);
  std::vector<int> tower_topl_idx(warp_size_int), plate_topl_idx(warp_size_int);
  std::vector<Tcalc> tower_atom_q(warp_size_int), plate_atom_q(warp_size_int);
  std::vector<int> tower_atom_ljidx(warp_size_int), plate_atom_ljidx(warp_size_int);
  const Tcalc tw_xshft = cgw->system_cell_invu[(system_idx * warp_size_int) + 6];
  const Tcalc tw_yshft = cgw->system_cell_invu[(system_idx * warp_size_int) + 7];
  const Tcalc tw_zshft = cgw->system_cell_invu[(system_idx * warp_size_int) + 8];

  // Compute interactions between the tower and the plate.
  while (tw_base < ntower) {

    // Take in the tower atoms's coordinates and arrange them as appropriate relative to the
    // central cell.  Read their topological indices.  This is done as a batch, as it will be on
    // the GPU, here for the purpose of re-using the cell homes found for the plate atoms.
    const int tw_batch = (ntower - tw_base > warp_size_int) ? warp_size_int : ntower - tw_base;
    for (int i = 0; i < tw_batch; i++) {
      const int tower_cell = findBin(tower_prefix, tw_base + i);
      const int tower_cell_pos = tw_base + i - tower_prefix[tower_cell];
      tower_topl_idx[i] = cgw->nonimg_atom_idx[tower_atoms[tower_cell].x + tower_cell_pos];
      tower_img_idx[i] = tower_atoms[tower_cell].x + tower_cell_pos;
      const Tcoord4 tw_xyzq = cgw->image[tower_img_idx[i]];
      if (tcoord_is_real) {
        tw_xpos[i] = tw_xyzq.x;
        tw_ypos[i] = tw_xyzq.y;
        tw_zpos[i] = tw_xyzq.z;
      }
      else {
        tw_xpos[i] = static_cast<Tcalc>(tw_xyzq.x) * cgw->inv_lpos_scale;
        tw_ypos[i] = static_cast<Tcalc>(tw_xyzq.y) * cgw->inv_lpos_scale;
        tw_zpos[i] = static_cast<Tcalc>(tw_xyzq.z) * cgw->inv_lpos_scale;
      }
      const Tcalc dtw_cell = tower_cell - 2;
      tw_xpos[i] += dtw_cell * tw_xshft;
      tw_ypos[i] += dtw_cell * tw_yshft;
      tw_zpos[i] += dtw_cell * tw_zshft;
      switch (theme) {
      case NonbondedTheme::ELECTROSTATIC:
      case NonbondedTheme::ALL:
        tower_atom_q[i] = sourceMagnitude(NonbondedTheme::ELECTROSTATIC, cgw->theme, tw_xyzq.w,
                                          tcoord_is_real, system_idx, poly_nbk) * poly_nbk.coulomb;
        break;
      case NonbondedTheme::VAN_DER_WAALS:
        break;
      }
      switch (theme) {
      case NonbondedTheme::VAN_DER_WAALS:
      case NonbondedTheme::ALL:
        tower_atom_ljidx[i] = (sourceIndex(NonbondedTheme::VAN_DER_WAALS, cgw->theme, tw_xyzq.w,
                                           tcoord_is_real) * nsys_lj_types) + sys_ljidx_offset;
        break;
      case NonbondedTheme::ELECTROSTATIC:
        break;
      }
    }
    
    // To apply a "trivial rejection" test, based on whether atoms in the plate might be anywhere
    // within range of some atom in the tower group, does not appear to be worthwile--it could
    // save perhaps 15-20% of the particles at what could be an even greater cost in terms of
    // additional evaluations and threads waiting, all the while increasing the complexity of the
    // code.
    int pl_base = 0;
    const Tcalc pli_xshft = cgw->system_cell_invu[(system_idx * warp_size_int) + 0];
    const Tcalc pli_yshft = cgw->system_cell_invu[(system_idx * warp_size_int) + 1];
    const Tcalc pli_zshft = cgw->system_cell_invu[(system_idx * warp_size_int) + 2];
    const Tcalc plj_xshft = cgw->system_cell_invu[(system_idx * warp_size_int) + 3];
    const Tcalc plj_yshft = cgw->system_cell_invu[(system_idx * warp_size_int) + 4];
    const Tcalc plj_zshft = cgw->system_cell_invu[(system_idx * warp_size_int) + 5];
    while (pl_base < nplate) {
      const int pl_batch = (nplate - pl_base > warp_size_int) ? warp_size_int : nplate - pl_base;
      for (int i = 0; i < pl_batch; i++) {
        const int plate_cell = findBin(plate_prefix, pl_base + i);
        const int plate_cell_pos = pl_base + i - plate_prefix[plate_cell];
        plate_topl_idx[i] = cgw->nonimg_atom_idx[plate_atoms[plate_cell].x + plate_cell_pos];
        plate_img_idx[i] = plate_atoms[plate_cell].x + plate_cell_pos;
        const Tcoord4 pl_xyzq = cgw->image[plate_img_idx[i]];
        if (tcoord_is_real) {
          pl_xpos[i] = pl_xyzq.x;
          pl_ypos[i] = pl_xyzq.y;
          pl_zpos[i] = pl_xyzq.z;
        }
        else {
          pl_xpos[i] = static_cast<Tcalc>(pl_xyzq.x) * cgw->inv_lpos_scale;
          pl_ypos[i] = static_cast<Tcalc>(pl_xyzq.y) * cgw->inv_lpos_scale;
          pl_zpos[i] = static_cast<Tcalc>(pl_xyzq.z) * cgw->inv_lpos_scale;
        }
        const int plj = plate_cell / 5;
        const int pli = plate_cell - (5 * plj);
        const Tcalc dipl_cell = pli - 2;
        const Tcalc djpl_cell = plj - 2;
        pl_xpos[i] += (dipl_cell * pli_xshft) + (djpl_cell * plj_xshft);
        pl_ypos[i] += (dipl_cell * pli_yshft) + (djpl_cell * plj_yshft);
        pl_zpos[i] += (dipl_cell * pli_zshft) + (djpl_cell * plj_zshft);
        switch (theme) {
        case NonbondedTheme::ELECTROSTATIC:
        case NonbondedTheme::ALL:
          plate_atom_q[i] = sourceMagnitude(NonbondedTheme::ELECTROSTATIC, cgw->theme, pl_xyzq.w,
                                            tcoord_is_real, system_idx, poly_nbk);
          break;
        case NonbondedTheme::VAN_DER_WAALS:
          break;
        }
        switch (theme) {
        case NonbondedTheme::VAN_DER_WAALS:
        case NonbondedTheme::ALL:
          plate_atom_ljidx[i] = sourceIndex(NonbondedTheme::VAN_DER_WAALS, cgw->theme, pl_xyzq.w,
                                            tcoord_is_real);
          break;
        case NonbondedTheme::ELECTROSTATIC:
          break;
        }
      }

      // Loop over tiles of the tower and plate.  These will never be on-diagonal tiles.
      const double2 tmp_e = basicTileInteractions(tw_xpos, tw_ypos, tw_zpos, pl_xpos, pl_ypos,
                                                  pl_zpos, tower_atom_q, plate_atom_q,
                                                  tower_atom_ljidx, plate_atom_ljidx,
                                                  tower_topl_idx, plate_topl_idx, tower_img_idx,
                                                  plate_img_idx, system_idx, tw_batch, pl_batch,
                                                  false, small_box, cgw, poly_psr, poly_nbk, lemr,
                                                  elec_cutsq, vdw_cutsq, qqew_coeff, ljew_coeff,
                                                  vdw_sum, eval_frc, theme);
      result.x += tmp_e.x;
      result.y += tmp_e.y;
      pl_base += pl_batch;
    }
    tw_base += tw_batch;
  }

  // Compute interactions within the home cell, always the third cell of the tower.
  std::vector<Tcalc> cr_xpos(warp_size_int), cr_ypos(warp_size_int), cr_zpos(warp_size_int);
  std::vector<Tcalc> xc_xpos(warp_size_int), xc_ypos(warp_size_int), xc_zpos(warp_size_int);
  std::vector<uint> core_img_idx(warp_size_int), xcore_img_idx(warp_size_int);
  std::vector<int> core_topl_idx(warp_size_int), xcore_topl_idx(warp_size_int);
  std::vector<Tcalc> core_atom_q(warp_size_int), xcore_atom_q(warp_size_int);
  std::vector<int> core_atom_ljidx(warp_size_int), xcore_atom_ljidx(warp_size_int);
  const int ncore = (tower_atoms[2].y >> 16);
  int cr_base = 0;
  while (cr_base < ncore) {
    const int cr_batch = (ncore - cr_base > warp_size_int) ? warp_size_int : ncore - cr_base;
    for (int i = 0; i < cr_batch; i++) {
      core_img_idx[i] = tower_atoms[2].x + cr_base + i;
      core_topl_idx[i] = cgw->nonimg_atom_idx[core_img_idx[i]];
      xcore_img_idx[i] = core_img_idx[i];
      xcore_topl_idx[i] = core_topl_idx[i];
      const Tcoord4 cr_xyzq = cgw->image[core_img_idx[i]];
      if (tcoord_is_real) {
        cr_xpos[i] = cr_xyzq.x;
        cr_ypos[i] = cr_xyzq.y;
        cr_zpos[i] = cr_xyzq.z;
      }
      else {
        cr_xpos[i] = static_cast<Tcalc>(cr_xyzq.x) * cgw->inv_lpos_scale;
        cr_ypos[i] = static_cast<Tcalc>(cr_xyzq.y) * cgw->inv_lpos_scale;
        cr_zpos[i] = static_cast<Tcalc>(cr_xyzq.z) * cgw->inv_lpos_scale;
      }
      xc_xpos[i] = cr_xpos[i];
      xc_ypos[i] = cr_ypos[i];
      xc_zpos[i] = cr_zpos[i];
      switch (theme) {
      case NonbondedTheme::ELECTROSTATIC:
      case NonbondedTheme::ALL:
        xcore_atom_q[i] = sourceMagnitude(NonbondedTheme::ELECTROSTATIC, cgw->theme, cr_xyzq.w,
                                          tcoord_is_real, system_idx, poly_nbk);
        core_atom_q[i] = xcore_atom_q[i] * poly_nbk.coulomb;
        break;
      case NonbondedTheme::VAN_DER_WAALS:
        break;
      }
      switch (theme) {
      case NonbondedTheme::VAN_DER_WAALS:
      case NonbondedTheme::ALL:
        xcore_atom_ljidx[i] = sourceIndex(NonbondedTheme::VAN_DER_WAALS, cgw->theme, cr_xyzq.w,
                                          tcoord_is_real);
        core_atom_ljidx[i] = (xcore_atom_ljidx[i] * nsys_lj_types) + sys_ljidx_offset;
        break;
      case NonbondedTheme::ELECTROSTATIC:
        break;
      }
    }

    // The first tile in this line will be self-interacting.
    double2 tmp_e = basicTileInteractions(cr_xpos, cr_ypos, cr_zpos, xc_xpos, xc_ypos, xc_zpos,
                                          core_atom_q, xcore_atom_q, core_atom_ljidx,
                                          xcore_atom_ljidx, core_topl_idx, xcore_topl_idx,
                                          core_img_idx, xcore_img_idx, system_idx, cr_batch,
                                          cr_batch, true, small_box, cgw, poly_psr, poly_nbk, lemr,
                                          elec_cutsq, vdw_cutsq, qqew_coeff, ljew_coeff, vdw_sum,
                                          eval_frc, theme);
    result.x += tmp_e.x;
    result.y += tmp_e.y;

    // Subsequent tiles will not be self-interacting.
    int xc_base = cr_batch;
    while (xc_base < ncore) {
      const int xc_batch = (ncore - xc_base > warp_size_int) ? warp_size_int : ncore - xc_base;
      for (int i = 0; i < xc_batch; i++) {
        xcore_img_idx[i] = tower_atoms[2].x + xc_base + i;
        xcore_topl_idx[i] = cgw->nonimg_atom_idx[xcore_img_idx[i]];
        const Tcoord4 xc_xyzq = cgw->image[xcore_img_idx[i]];
        if (tcoord_is_real) {
          xc_xpos[i] = xc_xyzq.x;
          xc_ypos[i] = xc_xyzq.y;
          xc_zpos[i] = xc_xyzq.z;
        }
        else {
          xc_xpos[i] = static_cast<Tcalc>(xc_xyzq.x) * cgw->inv_lpos_scale;
          xc_ypos[i] = static_cast<Tcalc>(xc_xyzq.y) * cgw->inv_lpos_scale;
          xc_zpos[i] = static_cast<Tcalc>(xc_xyzq.z) * cgw->inv_lpos_scale;
        }
        switch (theme) {
        case NonbondedTheme::ELECTROSTATIC:
        case NonbondedTheme::ALL:
          xcore_atom_q[i] = sourceMagnitude(NonbondedTheme::ELECTROSTATIC, cgw->theme, xc_xyzq.w,
                                            tcoord_is_real, system_idx, poly_nbk);
          break;
        case NonbondedTheme::VAN_DER_WAALS:
          break;
        }
        switch (theme) {
        case NonbondedTheme::VAN_DER_WAALS:
        case NonbondedTheme::ALL:
          xcore_atom_ljidx[i] = sourceIndex(NonbondedTheme::VAN_DER_WAALS, cgw->theme, xc_xyzq.w,
                                            tcoord_is_real);
          break;
        case NonbondedTheme::ELECTROSTATIC:
          break;
        }
      }
      tmp_e = basicTileInteractions(cr_xpos, cr_ypos, cr_zpos, xc_xpos, xc_ypos, xc_zpos,
                                    core_atom_q, xcore_atom_q, core_atom_ljidx, xcore_atom_ljidx,
                                    core_topl_idx, xcore_topl_idx, core_img_idx, xcore_img_idx,
                                    system_idx, cr_batch, xc_batch, false, small_box, cgw,
                                    poly_psr, poly_nbk, lemr, elec_cutsq, vdw_cutsq, qqew_coeff,
                                    ljew_coeff, vdw_sum, eval_frc, theme);
      result.x += tmp_e.x;
      result.y += tmp_e.y;
      xc_base += xc_batch;
    }
    cr_base += cr_batch;
  }

  // Compute interactions between the home cell and the lower half of the tower.
  std::vector<Tcalc> ft_xpos(warp_size_int), ft_ypos(warp_size_int), ft_zpos(warp_size_int);
  std::vector<uint> foot_img_idx(warp_size_int);
  std::vector<int> foot_topl_idx(warp_size_int);
  std::vector<Tcalc> foot_atom_q(warp_size_int);
  std::vector<int> foot_atom_ljidx(warp_size_int);
  int nfoot = (tower_atoms[0].y >> 16) + (tower_atoms[1].y >> 16);
  cr_base = 0;
  while (cr_base < ncore) {
    const int cr_batch = (ncore - cr_base > warp_size_int) ? warp_size_int : ncore - cr_base;
    for (int i = 0; i < cr_batch; i++) {
      core_img_idx[i] = tower_atoms[2].x + cr_base + i;
      core_topl_idx[i] = cgw->nonimg_atom_idx[core_img_idx[i]];
      const Tcoord4 cr_xyzq = cgw->image[core_img_idx[i]];
      if (tcoord_is_real) {
        cr_xpos[i] = cr_xyzq.x;
        cr_ypos[i] = cr_xyzq.y;
        cr_zpos[i] = cr_xyzq.z;
      }
      else {
        cr_xpos[i] = static_cast<Tcalc>(cr_xyzq.x) * cgw->inv_lpos_scale;
        cr_ypos[i] = static_cast<Tcalc>(cr_xyzq.y) * cgw->inv_lpos_scale;
        cr_zpos[i] = static_cast<Tcalc>(cr_xyzq.z) * cgw->inv_lpos_scale;
      }
      switch (theme) {
      case NonbondedTheme::ELECTROSTATIC:
      case NonbondedTheme::ALL:
        core_atom_q[i] = sourceMagnitude(NonbondedTheme::ELECTROSTATIC, cgw->theme, cr_xyzq.w,
                                         tcoord_is_real, system_idx, poly_nbk) * poly_nbk.coulomb;
        break;
      case NonbondedTheme::VAN_DER_WAALS:
        break;
      }
      switch (theme) {
      case NonbondedTheme::VAN_DER_WAALS:
      case NonbondedTheme::ALL:
        core_atom_ljidx[i] = (sourceIndex(NonbondedTheme::VAN_DER_WAALS, cgw->theme, cr_xyzq.w,
                                          tcoord_is_real) * nsys_lj_types) +
                             sys_ljidx_offset;
        break;
      case NonbondedTheme::ELECTROSTATIC:
        break;
      }
    }
    
    // No tiles will be self-interacting.
    int ft_base = (ncell_c > 4 || cell_cidx >= 2) ? tower_prefix[0] : tower_prefix[1];
    while (ft_base < nfoot) {
      const int ft_batch = (nfoot - ft_base > warp_size_int) ? warp_size_int : nfoot - ft_base;
      for (int i = 0; i < ft_batch; i++) {
        const int foot_cell = findBin(tower_prefix, ft_base + i);
        const int foot_cell_pos = ft_base + i - tower_prefix[foot_cell];
        foot_img_idx[i] = tower_atoms[foot_cell].x + foot_cell_pos;
        foot_topl_idx[i] = cgw->nonimg_atom_idx[foot_img_idx[i]];
        const Tcoord4 ft_xyzq = cgw->image[foot_img_idx[i]];
        if (tcoord_is_real) {
          ft_xpos[i] = ft_xyzq.x;
          ft_ypos[i] = ft_xyzq.y;
          ft_zpos[i] = ft_xyzq.z;
        }
        else {
          ft_xpos[i] = static_cast<Tcalc>(ft_xyzq.x) * cgw->inv_lpos_scale;
          ft_ypos[i] = static_cast<Tcalc>(ft_xyzq.y) * cgw->inv_lpos_scale;
          ft_zpos[i] = static_cast<Tcalc>(ft_xyzq.z) * cgw->inv_lpos_scale;
        }
        const Tcalc dft_cell = foot_cell - 2;
        ft_xpos[i] += dft_cell * tw_xshft;
        ft_ypos[i] += dft_cell * tw_yshft;
        ft_zpos[i] += dft_cell * tw_zshft;
        switch (theme) {
        case NonbondedTheme::ELECTROSTATIC:
        case NonbondedTheme::ALL:
          foot_atom_q[i] = sourceMagnitude(NonbondedTheme::ELECTROSTATIC, cgw->theme, ft_xyzq.w,
                                           tcoord_is_real, system_idx, poly_nbk);
          break;
        case NonbondedTheme::VAN_DER_WAALS:
          break;
        }
        switch (theme) {
        case NonbondedTheme::VAN_DER_WAALS:
        case NonbondedTheme::ALL:
          foot_atom_ljidx[i] = sourceIndex(NonbondedTheme::VAN_DER_WAALS, cgw->theme, ft_xyzq.w,
                                           tcoord_is_real);
          break;
        case NonbondedTheme::ELECTROSTATIC:
          break;
        }
      }
      const double2 tmp_e = basicTileInteractions(cr_xpos, cr_ypos, cr_zpos, ft_xpos, ft_ypos,
                                                  ft_zpos, core_atom_q, foot_atom_q,
                                                  core_atom_ljidx, foot_atom_ljidx, core_topl_idx,
                                                  foot_topl_idx, core_img_idx, foot_img_idx,
                                                  system_idx, cr_batch, ft_batch, false, small_box,
                                                  cgw, poly_psr, poly_nbk, lemr, elec_cutsq,
                                                  vdw_cutsq, qqew_coeff, ljew_coeff, vdw_sum,
                                                  eval_frc, theme);

      result.x += tmp_e.x;
      result.y += tmp_e.y;
      ft_base += ft_batch;
    }
    cr_base += cr_batch;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcalc2, typename Tcoord4>
void evaluateParticleParticleEnergy(const PhaseSpaceSynthesis &poly_ps,
                                    CellGridWriter<void, void, void, void> *cgw_v,
                                    ScoreCardWriter *scw, const PsSynthesisReader &poly_psr,
                                    const SyNonbondedKit<Tcalc, Tcalc2> &poly_nbk,
                                    const LocalExclusionMaskReader &lemr,
                                    const Tcalc elec_cutoff, const Tcalc vdw_cutoff,
                                    const Tcalc qqew_coeff, const Tcalc ljew_coeff,
                                    const VdwSumMethod vdw_sum, const EvaluateForce eval_frc,
                                    const NonbondedTheme theme) {

  // Restore the type of the CellGrid abstract.
  CellGridWriter<Tcoord, Tacc, Tcalc, Tcoord4> cgw = restoreType<Tcoord, Tacc,
                                                                 Tcalc, Tcoord4>(cgw_v);

  // Check that the requested type of calculation is feasible with the supplied CellGrid.
  switch (cgw.theme) {
  case NonbondedTheme::ELECTROSTATIC:
  case NonbondedTheme::VAN_DER_WAALS:
    if (theme != cgw.theme) {
      rtErr("Calculation of " + getEnumerationName(theme) + " is not possible with a CellGrid "
            "object containing " + getEnumerationName(cgw.theme) + " particle information",
            "evaluateParticleParticleEnergy");
    }
    break;
  case NonbondedTheme::ALL:
    break;
  }
  
  // Loop over all systems in the associated synthesis.  For each system, loop over the cells in
  // the system and, based on the size of the cell grid, determine the neighbors to import and
  // calculate pair interactions.  Each cell is at leat half the non-bonded cutoff between any two
  // opposing faces, implying that interactions must be computed among 63 pairs of cells (including
  // the self-interaction of each cell with itself).
  for (int sys_idx = 0; sys_idx < cgw.system_count; sys_idx++) {
    double2 result = { 0.0, 0.0 };
    const ullint footprint = cgw.system_cell_grids[sys_idx];
    const int ncell_a = ((footprint >> 28) & 0xfff);
    const int ncell_b = ((footprint >> 40) & 0xfff);
    const int ncell_c = ((footprint >> 52) & 0xfff);
    const int cell_count = ncell_a * ncell_b * ncell_c;
    for (int i = 0; i < cell_count; i++) {
      const double2 tmp_e = towerPlatePairInteractions(&cgw, poly_ps.data(), sys_idx, i, poly_nbk,
                                                       lemr, elec_cutoff, vdw_cutoff, qqew_coeff,
                                                       ljew_coeff, vdw_sum, eval_frc, theme);
      result.x += tmp_e.x;
      result.y += tmp_e.y;
    }
    bool add_elec = true;
    bool add_vdw = true;
    switch (theme) {
    case NonbondedTheme::ELECTROSTATIC:
      add_vdw = false;
      break;
    case NonbondedTheme::VAN_DER_WAALS:
      add_elec = false;
      break;
    case NonbondedTheme::ALL:
      break;
    }
    if (add_elec) {
      add(scw, StateVariable::ELECTROSTATIC, llround(result.x * scw->nrg_scale_lf), sys_idx);
    }
    if (add_vdw) {
      add(scw, StateVariable::VDW, llround(result.y * scw->nrg_scale_lf), sys_idx);
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void evaluateParticleParticleEnergy(CellGrid<Tcoord, Tacc, Tcalc, Tcoord4> *cg, ScoreCard *sc,
                                    const LocalExclusionMask &lema, const Tcalc elec_cutoff,
                                    const Tcalc vdw_cutoff, const Tcalc qqew_coeff,
                                    const Tcalc ljew_coeff, const VdwSumMethod vdw_sum,
                                    const EvaluateForce eval_frc, const NonbondedTheme theme) {
  const PhaseSpaceSynthesis *poly_ps = cg->getCoordinateSynthesisPointer();
  const AtomGraphSynthesis *poly_ag = cg->getTopologySynthesisPointer();
  CellGridWriter<void, void, void, void> cgw_v = cg->templateFreeData();
  ScoreCardWriter scw = sc->data();
  if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
    evaluateParticleParticleEnergy<Tcoord,
                                   Tacc,
                                   double,
                                   double2,
                                   Tcoord4>(*poly_ps, &cgw_v, &scw, poly_ps->data(),
                                            poly_ag->getDoublePrecisionNonbondedKit(), lema.data(),
                                            elec_cutoff, vdw_cutoff, qqew_coeff, ljew_coeff,
                                            vdw_sum, eval_frc, theme);
  }
  else {
    evaluateParticleParticleEnergy<Tcoord,
                                   Tacc,
                                   float,
                                   float2,
                                   Tcoord4>(*poly_ps, &cgw_v, &scw, poly_ps->data(),
                                            poly_ag->getSinglePrecisionNonbondedKit(), lema.data(),
                                            elec_cutoff, vdw_cutoff, qqew_coeff, ljew_coeff,
                                            vdw_sum, eval_frc, theme);
  }
}

} // namespace energy
} // namespace stormm
