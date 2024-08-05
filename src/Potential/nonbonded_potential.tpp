// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
double2 evaluateNonbondedEnergy(const NonbondedKit<Tcalc> nbk, const StaticExclusionMaskReader ser,
                                const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                                const double* umat, const double* invu,
                                const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                                Tforce* zfrc, ScoreCard *ecard,
                                const EvaluateForce eval_elec_force,
                                const EvaluateForce eval_vdw_force, const int system_index,
                                const Tcalc inv_gpos_factor, const Tcalc force_factor,
                                const Tcalc clash_distance, const Tcalc clash_ratio) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const bool tforce_is_sgnint = isSignedIntegralScalarType<Tforce>();
  const Tcalc value_zero  = 0.0;
  const Tcalc value_nil   = 1.0e-6;
  const Tcalc value_one   = 1.0;
  const Tcalc value_three = 3.0;
  const Tcalc value_four  = 4.0;
  
  // Initialize the energy result as two separate accumulators for each of a pair of quantities
  double ele_energy = 0.0;
  double vdw_energy = 0.0;
  llint ele_acc = 0LL;
  llint vdw_acc = 0LL;
  const Tcalc nrg_scale_factor = ecard->getEnergyScalingFactor<Tcalc>();
  const bool clash_check = (clash_distance > value_nil || clash_ratio > value_nil);
  
  // Allocate arrays for tile coordinates and accumulated forces, akin to GPU protocols.
  std::vector<Tcalc> cachi_xcrd(tile_length), cachi_ycrd(tile_length), cachi_zcrd(tile_length);
  std::vector<Tcalc> cachj_xcrd(tile_length), cachj_ycrd(tile_length), cachj_zcrd(tile_length);
  std::vector<Tcalc> cachi_q(tile_length), cachj_q(tile_length);
  std::vector<int> cachi_ljidx(tile_length), cachj_ljidx(tile_length);
  std::vector<Tcalc> cachi_xfrc(tile_length), cachi_yfrc(tile_length), cachi_zfrc(tile_length);
  std::vector<Tcalc> cachj_xfrc(tile_length), cachj_yfrc(tile_length), cachj_zfrc(tile_length);
  
  // Perform nested loops over all supertiles and all tiles within them
  for (int sti = 0; sti < ser.supertile_stride_count; sti++) {
    for (int stj = 0; stj <= sti; stj++) {

      // Tile dimensions and locations
      const int stni_atoms = std::min(ser.natom - (supertile_length * sti), supertile_length);
      const int stnj_atoms = std::min(ser.natom - (supertile_length * stj), supertile_length);
      const int ni_tiles = (stni_atoms + tile_length - 1) / tile_length;
      const int nj_tiles = (stnj_atoms + tile_length - 1) / tile_length;

      // Access the supertile's map index: if zero, there are no exclusions to worry about
      const int stij_map_index = ser.supertile_map_idx[(stj * ser.supertile_stride_count) + sti];
      const int diag_supertile = (sti == stj);

      // The outer loops can proceed until the branch about exclusions
      for (int ti = 0; ti < ni_tiles; ti++) {
        const int ni_atoms = std::min(stni_atoms - (ti * tile_length), tile_length);
        const int tjlim = (nj_tiles * (1 - diag_supertile)) + (diag_supertile * ti);
        for (int tj = 0; tj <= tjlim; tj++) {
          const int nj_atoms = std::min(stnj_atoms - (tj * tile_length), tile_length);
          const int diag_tile = diag_supertile * (ti == tj);

          // Pre-cache the atom positions to avoid conversions in the inner loops.  Pre-cache
          // properties and local force accumulators to mimic GPU activity.
          for (int i = 0; i < ni_atoms; i++) {
            const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
            if (tcoord_is_sgnint) {
              cachi_xcrd[i] = static_cast<Tcalc>(xcrd[atom_i]) * inv_gpos_factor;
              cachi_ycrd[i] = static_cast<Tcalc>(ycrd[atom_i]) * inv_gpos_factor;
              cachi_zcrd[i] = static_cast<Tcalc>(zcrd[atom_i]) * inv_gpos_factor;
            }
            else {
              cachi_xcrd[i] = xcrd[atom_i];
              cachi_ycrd[i] = ycrd[atom_i];
              cachi_zcrd[i] = zcrd[atom_i];
            }
            cachi_xfrc[i] = value_zero;
            cachi_yfrc[i] = value_zero;
            cachi_zfrc[i] = value_zero;
            cachi_q[i] = nbk.coulomb_constant * nbk.charge[atom_i];
            cachi_ljidx[i] = nbk.lj_idx[atom_i];
          }
          for (int j = 0; j < nj_atoms; j++) {
            const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
            if (tcoord_is_sgnint) {
              cachj_xcrd[j] = static_cast<Tcalc>(xcrd[atom_j]) * inv_gpos_factor;
              cachj_ycrd[j] = static_cast<Tcalc>(ycrd[atom_j]) * inv_gpos_factor;
              cachj_zcrd[j] = static_cast<Tcalc>(zcrd[atom_j]) * inv_gpos_factor;
            }
            else {
              cachj_xcrd[j] = xcrd[atom_j];
              cachj_ycrd[j] = ycrd[atom_j];
              cachj_zcrd[j] = zcrd[atom_j];
            }
            cachj_xfrc[j] = value_zero;
            cachj_yfrc[j] = value_zero;
            cachj_zfrc[j] = value_zero;
            cachj_q[j] = nbk.charge[atom_j];
            cachj_ljidx[j] = nbk.lj_idx[atom_j];              
          }

          // Branch for different types of tiles
          if (stij_map_index == 0) {

            // Branch for clash detection and mitigation
            if (clash_check) {

              // Exclusion-free loops with clash checking
              for (int i = 0; i < ni_atoms; i++) {
                const Tcalc qi = cachi_q[i];
                const int ljt_i = cachi_ljidx[i];
                const int jlim = (nj_atoms * (1 - diag_tile)) + (diag_tile * i);
                const Tcalc atomi_x = cachi_xcrd[i];
                const Tcalc atomi_y = cachi_ycrd[i];
                const Tcalc atomi_z = cachi_zcrd[i];
                Tcalc ele_contrib = value_zero;
                Tcalc vdw_contrib = value_zero;
                for (int j = 0; j < jlim; j++) {
                  const Tcalc dx = cachj_xcrd[j] - atomi_x;
                  const Tcalc dy = cachj_ycrd[j] - atomi_y;
                  const Tcalc dz = cachj_zcrd[j] - atomi_z;
                  const Tcalc r = (tcalc_is_double) ? sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                                                      sqrtf((dx * dx) + (dy * dy) + (dz * dz));
                  Tcalc fmag = value_zero;
                  const Tcalc qiqj = qi * cachj_q[j];
                  if (eval_elec_force == EvaluateForce::YES) {
                    quadraticCoreElectrostatics<Tcalc>(r, clash_distance, qiqj, &ele_contrib,
                                                       &fmag);
                  }
                  else {
                    quadraticCoreElectrostatics<Tcalc>(r, clash_distance, qiqj, &ele_contrib,
                                                       nullptr);
                  }
                  const int ljt_j = cachj_ljidx[j];
                  const Tcalc lja = nbk.lja_coeff[(ljt_j * nbk.n_lj_types) + ljt_i];
                  const Tcalc ljb = nbk.ljb_coeff[(ljt_j * nbk.n_lj_types) + ljt_i];
                  if (eval_vdw_force == EvaluateForce::YES) {
                    quarticCoreLennardJones<Tcalc>(r, clash_ratio, lja, ljb,
                                                   nbk.lj_sigma[(ljt_j * nbk.n_lj_types) + ljt_i],
                                                   &vdw_contrib, &fmag);
                  }
                  else {
                    quarticCoreLennardJones<Tcalc>(r, clash_ratio, lja, ljb,
                                                   nbk.lj_sigma[(ljt_j * nbk.n_lj_types) + ljt_i],
                                                   &vdw_contrib, nullptr);
                  }

                  // Evaluate the force, if requested
                  if (eval_elec_force == EvaluateForce::YES ||
                      eval_vdw_force == EvaluateForce::YES) {
                    const Tcalc fmag_dx = fmag * dx;
                    const Tcalc fmag_dy = fmag * dy;
                    const Tcalc fmag_dz = fmag * dz;
                    cachi_xfrc[i] += fmag_dx;
                    cachi_yfrc[i] += fmag_dy;
                    cachi_zfrc[i] += fmag_dz;
                    cachj_xfrc[j] -= fmag_dx;
                    cachj_yfrc[j] -= fmag_dy;
                    cachj_zfrc[j] -= fmag_dz;
                  }
                }
            
                // Contribute what would be thread-accumulated energies to the totals
                ele_energy += ele_contrib;
                ele_acc += llround(ele_contrib * nrg_scale_factor);
                vdw_energy += vdw_contrib;
                vdw_acc += llround(vdw_contrib * nrg_scale_factor);
              }
            }
            else {

              // Exclusion-free loops without clash checking
              for (int i = 0; i < ni_atoms; i++) {
                const Tcalc qi = cachi_q[i];
                const int ljt_i = cachi_ljidx[i];
                const int jlim = (nj_atoms * (1 - diag_tile)) + (diag_tile * i);
                const Tcalc atomi_x = cachi_xcrd[i];
                const Tcalc atomi_y = cachi_ycrd[i];
                const Tcalc atomi_z = cachi_zcrd[i];
                Tcalc ele_contrib = value_zero;
                Tcalc vdw_contrib = value_zero;
                for (int j = 0; j < jlim; j++) {
                  const Tcalc dx = cachj_xcrd[j] - atomi_x;
                  const Tcalc dy = cachj_ycrd[j] - atomi_y;
                  const Tcalc dz = cachj_zcrd[j] - atomi_z;
                  Tcalc invr2, invr;
                  if (tcalc_is_double) {
                    invr2 = 1.0 / ((dx * dx) + (dy * dy) + (dz * dz));
                    invr = sqrt(invr2);
                  }
                  else {
                    invr2 = value_one / ((dx * dx) + (dy * dy) + (dz * dz));
                    invr = sqrtf(invr2);
                  }
                  const Tcalc invr4 = invr2 * invr2;
                  const Tcalc qiqj = qi * cachj_q[j];
                  ele_contrib += qiqj * invr;
                  const int ljt_j = cachj_ljidx[j];
                  const Tcalc lja = nbk.lja_coeff[(ljt_j * nbk.n_lj_types) + ljt_i];
                  const Tcalc ljb = nbk.ljb_coeff[(ljt_j * nbk.n_lj_types) + ljt_i];
                  vdw_contrib += ((lja * invr4 * invr4) - (ljb * invr2)) * invr4;

                  // Evaluate the force, if requested
                  if (eval_elec_force == EvaluateForce::YES ||
                      eval_vdw_force == EvaluateForce::YES) {
                    Tcalc fmag = (eval_elec_force == EvaluateForce::YES) ? -(qiqj * invr * invr2) :
                                                                           value_zero;
                    if (eval_vdw_force == EvaluateForce::YES) {
                      if (tcalc_is_double) {
                        fmag += ((6.0 * ljb) - (12.0 * lja * invr4 * invr2)) * invr4 * invr4;
                      }
                      else {
                        fmag += ((6.0f * ljb) - (12.0f * lja * invr4 * invr2)) * invr4 * invr4;
                      }
                    }
                    const Tcalc fmag_dx = fmag * dx;
                    const Tcalc fmag_dy = fmag * dy;
                    const Tcalc fmag_dz = fmag * dz;
                    cachi_xfrc[i] += fmag_dx;
                    cachi_yfrc[i] += fmag_dy;
                    cachi_zfrc[i] += fmag_dz;
                    cachj_xfrc[j] -= fmag_dx;
                    cachj_yfrc[j] -= fmag_dy;
                    cachj_zfrc[j] -= fmag_dz;
                  }
                }
            
                // Contribute what would be thread-accumulated energies to the totals
                ele_energy += ele_contrib;
                ele_acc += llround(ele_contrib * nrg_scale_factor);
                vdw_energy += vdw_contrib;
                vdw_acc += llround(vdw_contrib * nrg_scale_factor);
              }
            }
          }
          else {

            // Get the tile's mask to check exclusions with each interaction
            const int tij_map_index = ser.tile_map_idx[stij_map_index +
                                                       (tj * tile_lengths_per_supertile) + ti];
            if (clash_check) {

              // Perform clash checking in addition to exclusion checking
              for (int i = 0; i < ni_atoms; i++) {
                const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
                const uint mask_i = ser.mask_data[tij_map_index + i];
                const Tcalc atomi_x = cachi_xcrd[i];
                const Tcalc atomi_y = cachi_ycrd[i];
                const Tcalc atomi_z = cachi_zcrd[i];
                const Tcalc qi = cachi_q[i];
                const int ljt_i = cachi_ljidx[i];
                const int jlim = (nj_atoms * (1 - diag_tile)) + (diag_tile * i);
                Tcalc ele_contrib = value_zero;
                Tcalc vdw_contrib = value_zero;
                for (int j = 0; j < jlim; j++) {
                  if ((mask_i >> j) & 0x1) {
                    continue;
                  }
                  const Tcalc dx = cachj_xcrd[j] - atomi_x;
                  const Tcalc dy = cachj_ycrd[j] - atomi_y;
                  const Tcalc dz = cachj_zcrd[j] - atomi_z;
                  const Tcalc r = (tcalc_is_double) ? sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                                                      sqrtf((dx * dx) + (dy * dy) + (dz * dz));
                  Tcalc fmag = value_zero;
                  const Tcalc qiqj = qi * cachj_q[j];
                  if (eval_elec_force == EvaluateForce::YES) {
                    quadraticCoreElectrostatics<Tcalc>(r, clash_distance, qiqj, &ele_contrib,
                                                       &fmag);
                  }
                  else {
                    quadraticCoreElectrostatics<Tcalc>(r, clash_distance, qiqj, &ele_contrib,
                                                       nullptr);
                  }
                  const int ljt_j = cachj_ljidx[j];
                  const Tcalc lja = nbk.lja_coeff[(ljt_j * nbk.n_lj_types) + ljt_i];
                  const Tcalc ljb = nbk.ljb_coeff[(ljt_j * nbk.n_lj_types) + ljt_i];
                  if (eval_vdw_force == EvaluateForce::YES) {
                    quarticCoreLennardJones<Tcalc>(r, clash_ratio, lja, ljb,
                                                   nbk.lj_sigma[(ljt_j * nbk.n_lj_types) + ljt_i],
                                                   &vdw_contrib, &fmag);
                  }
                  else {
                    quarticCoreLennardJones<Tcalc>(r, clash_ratio, lja, ljb,
                                                   nbk.lj_sigma[(ljt_j * nbk.n_lj_types) + ljt_i],
                                                   &vdw_contrib, nullptr);
                  }

                  // Evaluate the force, if requested
                  if (eval_elec_force == EvaluateForce::YES ||
                      eval_vdw_force == EvaluateForce::YES) {
                    const Tcalc fmag_dx = fmag * dx;
                    const Tcalc fmag_dy = fmag * dy;
                    const Tcalc fmag_dz = fmag * dz;
                    cachi_xfrc[i] += fmag_dx;
                    cachi_yfrc[i] += fmag_dy;
                    cachi_zfrc[i] += fmag_dz;
                    cachj_xfrc[j] -= fmag_dx;
                    cachj_yfrc[j] -= fmag_dy;
                    cachj_zfrc[j] -= fmag_dz;
                  }                
                }

                // Contribute what would be thread-accumulated energies to the totals
                ele_energy += ele_contrib;
                ele_acc += llround(ele_contrib * nrg_scale_factor);
                vdw_energy += vdw_contrib;
                vdw_acc += llround(vdw_contrib * nrg_scale_factor);
              }
            }
            else {

              // Proceed without clash checking
              for (int i = 0; i < ni_atoms; i++) {
                const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
                const uint mask_i = ser.mask_data[tij_map_index + i];
                const Tcalc atomi_x = cachi_xcrd[i];
                const Tcalc atomi_y = cachi_ycrd[i];
                const Tcalc atomi_z = cachi_zcrd[i];
                const Tcalc qi = cachi_q[i];
                const int ljt_i = cachi_ljidx[i];
                const int jlim = (nj_atoms * (1 - diag_tile)) + (diag_tile * i);
                Tcalc ele_contrib = 0.0;
                Tcalc vdw_contrib = 0.0;
                for (int j = 0; j < jlim; j++) {
                  if ((mask_i >> j) & 0x1) {
                    continue;
                  }
                  const Tcalc dx = cachj_xcrd[j] - atomi_x;
                  const Tcalc dy = cachj_ycrd[j] - atomi_y;
                  const Tcalc dz = cachj_zcrd[j] - atomi_z;
                  Tcalc invr2, invr;
                  if (tcalc_is_double) {
                    invr2 = 1.0 / ((dx * dx) + (dy * dy) + (dz * dz));
                    invr = sqrt(invr2);
                  }
                  else {
                    invr2 = value_one / ((dx * dx) + (dy * dy) + (dz * dz));
                    invr = sqrtf(invr2);
                  }
                  const Tcalc invr4 = invr2 * invr2;
                  const Tcalc qiqj = qi * cachj_q[j];
                  ele_contrib += qiqj * invr;
                  const int ljt_j = cachj_ljidx[j];
                  const Tcalc lja = nbk.lja_coeff[(ljt_j * nbk.n_lj_types) + ljt_i];
                  const Tcalc ljb = nbk.ljb_coeff[(ljt_j * nbk.n_lj_types) + ljt_i];
                  vdw_contrib += ((lja * invr4 * invr4) - (ljb * invr2)) * invr4;

                  // Evaluate the force, if requested
                  if (eval_elec_force == EvaluateForce::YES ||
                      eval_vdw_force == EvaluateForce::YES) {
                    Tcalc fmag = (eval_elec_force == EvaluateForce::YES) ? -(qiqj * invr * invr2) :
                                                                           value_zero;
                    if (eval_vdw_force == EvaluateForce::YES) {
                      if (tcalc_is_double) {
                        fmag += ((6.0 * ljb) - (12.0 * lja * invr4 * invr2)) * invr4 * invr4;
                      }
                      else {
                        fmag += ((6.0f * ljb) - (12.0f * lja * invr4 * invr2)) * invr4 * invr4;
                      }
                    }
                    const Tcalc fmag_dx = fmag * dx;
                    const Tcalc fmag_dy = fmag * dy;
                    const Tcalc fmag_dz = fmag * dz;
                    cachi_xfrc[i] += fmag_dx;
                    cachi_yfrc[i] += fmag_dy;
                    cachi_zfrc[i] += fmag_dz;
                    cachj_xfrc[j] -= fmag_dx;
                    cachj_yfrc[j] -= fmag_dy;
                    cachj_zfrc[j] -= fmag_dz;
                  }                
                }

                // Contribute what would be thread-accumulated energies to the totals
                ele_energy += ele_contrib;
                ele_acc += llround(ele_contrib * nrg_scale_factor);
                vdw_energy += vdw_contrib;
                vdw_acc += llround(vdw_contrib * nrg_scale_factor);
              }
            }
          }

          // Contribute cached forces back to global arrays.
          if (eval_elec_force == EvaluateForce::YES || eval_vdw_force == EvaluateForce::YES) {
            for (int i = 0; i < ni_atoms; i++) {
              const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
              if (tforce_is_sgnint) {
                xfrc[atom_i] += llround(cachi_xfrc[i] * force_factor);              
                yfrc[atom_i] += llround(cachi_yfrc[i] * force_factor);
                zfrc[atom_i] += llround(cachi_zfrc[i] * force_factor);
              }
              else {
                xfrc[atom_i] += cachi_xfrc[i];
                yfrc[atom_i] += cachi_yfrc[i];
                zfrc[atom_i] += cachi_zfrc[i];
              }
            }
            for (int j = 0; j < nj_atoms; j++) {
              const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
              if (tforce_is_sgnint) {
                xfrc[atom_j] += llround(cachj_xfrc[j] * force_factor);              
                yfrc[atom_j] += llround(cachj_yfrc[j] * force_factor);
                zfrc[atom_j] += llround(cachj_zfrc[j] * force_factor);
              }
              else {
                xfrc[atom_j] += cachj_xfrc[j];
                yfrc[atom_j] += cachj_yfrc[j];
                zfrc[atom_j] += cachj_zfrc[j];
              }
            }
          }
        }
      }
    }
  }

  // Contribute results
  ecard->contribute(StateVariable::ELECTROSTATIC, ele_acc, system_index);
  ecard->contribute(StateVariable::VDW, vdw_acc, system_index);

  // Return the double-precision energy sums, if of interest
  return { ele_energy, vdw_energy };
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double2 evaluateNonbondedEnergy(const NonbondedKit<Tcalc> &nbk,
                                const StaticExclusionMaskReader &ser,
                                const CoordinateSeriesReader<Tcoord> csr, ScoreCard *ecard,
                                const int system_index, const Tcalc clash_distance,
                                const Tcalc clash_ratio) {
  const size_t atom_os = static_cast<size_t>(system_index) *
                         roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t xfrm_os = static_cast<size_t>(system_index) * roundUp<size_t>(9, warp_size_zu);
  return evaluateNonbondedEnergy<Tcoord, Tcoord, Tcalc>(nbk, ser, &csr.xcrd[atom_os],
                                                        &csr.ycrd[atom_os], &csr.zcrd[atom_os],
                                                        &csr.umat[xfrm_os], &csr.invu[xfrm_os],
                                                        csr.unit_cell, nullptr, nullptr, nullptr,
                                                        ecard, EvaluateForce::NO,
                                                        EvaluateForce::NO, system_index,
                                                        csr.inv_gpos_scale, clash_distance,
                                                        clash_ratio);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
double evaluateGeneralizedBornEnergy(const NonbondedKit<Tcalc> nbk,
                                     const StaticExclusionMaskReader ser,
                                     const ImplicitSolventKit<Tcalc> isk,
                                     const NeckGeneralizedBornKit<Tcalc> ngb_kit,
                                     const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                                     Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                                     Tforce *effective_gb_radii, Tforce *psi, Tforce *sumdeijda,
                                     ScoreCard *ecard, const EvaluateForce eval_force,
                                     const int system_index, const Tcalc inv_gpos_factor,
                                     const Tcalc force_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const bool tforce_is_sgnint = isSignedIntegralScalarType<Tforce>();
  const Tcalc v_zero = 0.0;
  const Tcalc v_one  = 1.0;
  const Tcalc v_two  = 2.0;
  const Tcalc v_thre = 3.0;
  const Tcalc v_four = 4.0;
  const Tcalc v_half = 0.5;
  const Tcalc v_qrtr = 0.25;
  const Tcalc v_pthr = 0.3;
  const Tcalc v_opei = 1.8;

  // Complete the implicit solvent model and initialize the energy result
  const ImplicitSolventRecipe<Tcalc> isr(isk, ngb_kit);
  double egb_energy = 0.0;
  llint egb_acc = 0LL;
  const Tcalc nrg_scale_factor = ecard->getEnergyScalingFactor<Tcalc>();

  // Initialize psi
  for (int i = 0; i < nbk.natom; i++) {
    psi[i] = 0.0;
  }
  
  // Return zero if there is no implicit solvent model.  This is a reference calculation and
  // intended to keep clean code.
  if (isr.igb == ImplicitSolventModel::NONE) {
    ecard->contribute(StateVariable::GENERALIZED_BORN, egb_acc, system_index);
    return egb_energy;
  }

  // Shorthand for some constants defined in Constants/generalized_born.h
  const Tcalc gta  = (tcalc_is_double) ? gb_taylor_a_lf  : gb_taylor_a_f;
  const Tcalc gtb  = (tcalc_is_double) ? gb_taylor_b_lf  : gb_taylor_b_f;
  const Tcalc gtc  = (tcalc_is_double) ? gb_taylor_c_lf  : gb_taylor_c_f;
  const Tcalc gtd  = (tcalc_is_double) ? gb_taylor_d_lf  : gb_taylor_d_f;
  const Tcalc gtdd = (tcalc_is_double) ? gb_taylor_dd_lf : gb_taylor_dd_f;
  const Tcalc gte  = (tcalc_is_double) ? gb_taylor_e_lf  : gb_taylor_e_f;
  const Tcalc gtf  = (tcalc_is_double) ? gb_taylor_f_lf  : gb_taylor_f_f;
  const Tcalc gtg  = (tcalc_is_double) ? gb_taylor_g_lf  : gb_taylor_g_f;
  const Tcalc gth  = (tcalc_is_double) ? gb_taylor_h_lf  : gb_taylor_h_f;
  const Tcalc gthh = (tcalc_is_double) ? gb_taylor_hh_lf : gb_taylor_hh_f;

  // Perform nested loops over all supertiles, and all tiles within them, to compute the GB radii
  std::vector<Tcalc> cachi_xcrd(tile_length), cachi_ycrd(tile_length), cachi_zcrd(tile_length);
  std::vector<Tcalc> cachj_xcrd(tile_length), cachj_ycrd(tile_length), cachj_zcrd(tile_length);
  std::vector<Tcalc> cachi_xfrc(tile_length), cachi_yfrc(tile_length), cachi_zfrc(tile_length);
  std::vector<Tcalc> cachj_xfrc(tile_length), cachj_yfrc(tile_length), cachj_zfrc(tile_length);
  std::vector<Tcalc> cachi_radii(tile_length), cachj_radii(tile_length);
  std::vector<Tcalc> cachi_screen(tile_length), cachj_screen(tile_length);
  std::vector<Tcalc> cachi_psi(tile_length), cachj_psi(tile_length);
  std::vector<int> cachi_neck_idx(tile_length), cachj_neck_idx(tile_length);
  for (int sti = 0; sti < ser.supertile_stride_count; sti++) {
    for (int stj = 0; stj <= sti; stj++) {

      // Tile dimensions and locations
      const int stni_atoms = std::min(ser.natom - (supertile_length * sti), supertile_length);
      const int stnj_atoms = std::min(ser.natom - (supertile_length * stj), supertile_length);
      const int ni_tiles = (stni_atoms + tile_length - 1) / tile_length;
      const int nj_tiles = (stnj_atoms + tile_length - 1) / tile_length;

      // Access the supertile's map index: if zero, there are no exclusions to worry about
      const int diag_supertile = (sti == stj);

      // The outer loops can proceed until the branch about exclusions
      for (int ti = 0; ti < ni_tiles; ti++) {
        const int ni_atoms = std::min(stni_atoms - (ti * tile_length), tile_length);
        const int tjlim = (nj_tiles * (1 - diag_supertile)) + (diag_supertile * ti);
        for (int tj = 0; tj <= tjlim; tj++) {
          const int nj_atoms = std::min(stnj_atoms - (tj * tile_length), tile_length);
          const int diag_tile = diag_supertile * (ti == tj);

          // Pre-cache the atom positions to avoid conversions in the inner loops.  Pre-cache
          // properties and local force accumulators to mimic GPU activity.
          for (int i = 0; i < ni_atoms; i++) {
            const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
            if (tcoord_is_sgnint) {
              cachi_xcrd[i] = static_cast<Tcalc>(xcrd[atom_i]) * inv_gpos_factor;
              cachi_ycrd[i] = static_cast<Tcalc>(ycrd[atom_i]) * inv_gpos_factor;
              cachi_zcrd[i] = static_cast<Tcalc>(zcrd[atom_i]) * inv_gpos_factor;
            }
            else {
              cachi_xcrd[i] = xcrd[atom_i];
              cachi_ycrd[i] = ycrd[atom_i];
              cachi_zcrd[i] = zcrd[atom_i];
            }
            cachi_radii[i]  = isr.pb_radii[atom_i] - isr.gb_offset;
            cachi_screen[i] = isr.gb_screen[atom_i];
            cachi_psi[i]    = 0.0;
            if (isr.igb == ImplicitSolventModel::NECK_GB ||
                isr.igb == ImplicitSolventModel::NECK_GB_II) {
              cachi_neck_idx[i] = isr.neck_gb_idx[atom_i];
            }
          }
          for (int j = 0; j < nj_atoms; j++) {
            const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
            if (tcoord_is_sgnint) {
              cachj_xcrd[j] = static_cast<Tcalc>(xcrd[atom_j]) * inv_gpos_factor;
              cachj_ycrd[j] = static_cast<Tcalc>(ycrd[atom_j]) * inv_gpos_factor;
              cachj_zcrd[j] = static_cast<Tcalc>(zcrd[atom_j]) * inv_gpos_factor;
            }
            else {
              cachj_xcrd[j] = xcrd[atom_j];
              cachj_ycrd[j] = ycrd[atom_j];
              cachj_zcrd[j] = zcrd[atom_j];
            }
            cachj_radii[j]  = isr.pb_radii[atom_j] - isr.gb_offset;
            cachj_screen[j] = isr.gb_screen[atom_j];
            cachj_psi[j]    = 0.0;
            if (isr.igb == ImplicitSolventModel::NECK_GB ||
                isr.igb == ImplicitSolventModel::NECK_GB_II) {
              cachj_neck_idx[j] = isr.neck_gb_idx[atom_j];
            }
          }

          // Loop over all atoms in the tile
          for (int i = 0; i < ni_atoms; i++) {
            const int jlim = (nj_atoms * (1 - diag_tile)) + (diag_tile * i);
            const Tcalc atomi_radius = cachi_radii[i];
            const Tcalc atomi_inv_radius = v_one / atomi_radius;
            for (int j = 0; j < jlim; j++) {
              const Tcalc dx   = cachi_xcrd[i] - cachj_xcrd[j];
              const Tcalc dy   = cachi_ycrd[i] - cachj_ycrd[j];
              const Tcalc dz   = cachi_zcrd[i] - cachj_zcrd[j];
              const Tcalc r2   = (dx * dx) + (dy * dy) + (dz * dz);
              const Tcalc r = (tcalc_is_double) ? sqrt(r2) : sqrtf(r2);
              const Tcalc invr = v_one / r;
              const Tcalc atomj_radius = cachj_radii[j];
              const Tcalc atomj_inv_radius = v_one / atomj_radius;

              // First computation: atom I -> atom J
              const Tcalc sj = cachj_screen[j] * atomj_radius;
              const Tcalc sj2 = sj * sj;
              if (r > v_four * sj) {
                const Tcalc invr2 = invr * invr;
                const Tcalc tmpsd = sj2 * invr2;
                const Tcalc dumbo = gta +
                                    tmpsd * (gtb + tmpsd * (gtc + tmpsd * (gtd + tmpsd * gtdd)));
                cachi_psi[i] -= sj * tmpsd * invr2 * dumbo;
              }
              else if (r > atomi_radius + sj) {
                if (tcalc_is_double) {
                  cachi_psi[i] -= v_half * ((sj / (r2 - sj2)) +
                                            (v_half * invr * log((r - sj) / (r + sj))));
                }
                else {
                  cachi_psi[i] -= v_half * ((sj / (r2 - sj2)) +
                                            (v_half * invr * logf((r - sj) / (r + sj))));
                }
              }
              else if (r > fabs(atomi_radius - sj)) {
                const Tcalc theta = v_half * atomi_inv_radius * invr *
                                    (r2 + (atomi_radius * atomi_radius) - sj2);
                const Tcalc uij   = v_one / (r + sj);
                if (tcalc_is_double) {
                  cachi_psi[i] -= v_qrtr * ((atomi_inv_radius * (v_two - theta)) -
                                            uij + (invr * log(atomi_radius * uij)));
                }
                else {
                  cachi_psi[i] -= v_qrtr * ((atomi_inv_radius * (v_two - theta)) -
                                            uij + (invr * logf(atomi_radius * uij)));
                }
              }
              else if (atomi_radius < sj) {
                if (tcalc_is_double) {
                  cachi_psi[i] -= v_half * ((sj / (r2 - sj2)) + (v_two * atomi_inv_radius) +
                                            (v_half * invr * log((sj - r) / (sj + r))));
                }
                else {
                  cachi_psi[i] -= v_half * ((sj / (r2 - sj2)) + (v_two * atomi_inv_radius) +
                                            (v_half * invr * logf((sj - r) / (sj + r))));
                }
              }

              // Second computation: atom J -> atom I
              const Tcalc si = cachi_screen[i] * atomi_radius;
              const Tcalc si2 = si * si;
              if (r > v_four * si) {
                const Tcalc invr2  = invr * invr;
                const Tcalc tmpsd  = si2 * invr2;
                const Tcalc dumbo  = gta +
                                     tmpsd * (gtb + tmpsd * (gtc + tmpsd * (gtd + tmpsd * gtdd)));
                cachj_psi[j] -= si * tmpsd * invr2 * dumbo;
              }
              else if (r > atomj_radius + si) {
                if (tcalc_is_double) {
                  cachj_psi[j] -= v_half * ((si / (r2 - si2)) +
                                            (v_half * invr * log((r - si) / (r + si))));
                }
                else {
                  cachj_psi[j] -= v_half * ((si / (r2 - si2)) +
                                            (v_half * invr * logf((r - si) / (r + si))));
                }
              }
              else if (r > fabs(atomj_radius - si)) {
                const Tcalc theta = v_half * atomj_inv_radius * invr *
                                    (r2 + (atomj_radius * atomj_radius) - si2);
                const Tcalc uij   = v_one / (r + si);
                if (tcalc_is_double) {
                  cachj_psi[j] -= v_qrtr * (atomj_inv_radius * (v_two - theta) - uij +
                                            invr * log(atomj_radius * uij));
                }
                else {
                  cachj_psi[j] -= v_qrtr * (atomj_inv_radius * (v_two - theta) - uij +
                                            invr * logf(atomj_radius * uij));
                }
              }
              else if (atomj_radius < si) {
                if (tcalc_is_double) {
                  cachj_psi[j] -= v_half * ((si / (r2 - si2)) + (v_two * atomj_inv_radius) +
                                            (v_half * invr * log((si - r) / (si + r))));
                }
                else {
                  cachj_psi[j] -= v_half * ((si / (r2 - si2)) + (v_two * atomj_inv_radius) +
                                            (v_half * invr * logf((si - r) / (si + r))));
                }
              }

              // Neck GB contribution
              if ((isr.igb == ImplicitSolventModel::NECK_GB ||
                   isr.igb == ImplicitSolventModel::NECK_GB_II) &&
                  r < cachi_radii[i] + cachj_radii[j] + (v_two * isr.gb_offset) + isr.gb_neckcut) {

                // First computation: atom I -> atom J
                const int ij_table_idx = (isr.table_size * cachj_neck_idx[j]) +
                                         cachi_neck_idx[i];
                Tcalc mdist  = r - isr.neck_max_sep[ij_table_idx];
                Tcalc mdist2 = mdist * mdist;
                Tcalc mdist6 = mdist2 * mdist2 * mdist2;
                const Tcalc ij_neck = isr.neck_max_val[ij_table_idx] /
                                              (v_one + mdist2 + (v_pthr * mdist6));
                cachi_psi[i] -= isr.gb_neckscale * ij_neck;

                // Second computation: atom J -> atom I
                const int ji_table_idx = (isr.table_size * cachi_neck_idx[i]) +
                                         cachj_neck_idx[j];
                mdist  = r - isr.neck_max_sep[ji_table_idx];
                mdist2 = mdist * mdist;
                mdist6 = mdist2 * mdist2 * mdist2;
                const Tcalc ji_neck = isr.neck_max_val[ji_table_idx] /
                                      (v_one + mdist2 + (v_pthr * mdist6));
                cachj_psi[j] -= isr.gb_neckscale * ji_neck;
              }
            }
          }

          // Contribute the local psi accumulators
          for (int i = 0; i < ni_atoms; i++) {
            const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
            if (tforce_is_sgnint) {
              psi[atom_i] += llround(cachi_psi[i] * force_factor);
            }
            else {
              psi[atom_i] += cachi_psi[i];
            }
          }
          for (int j = 0; j < nj_atoms; j++) {
            const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
            if (tforce_is_sgnint) {
              psi[atom_j] += llround(cachj_psi[j] * force_factor);
            }
            else {
              psi[atom_j] += cachj_psi[j];
            }
          }
        }
      }
    }
  }

  // Make a second pass to finalize the effective GB radii
  for (int i = 0; i < nbk.natom; i++) {
    switch (isr.igb) {
    case ImplicitSolventModel::HCT_GB:
      {
        // Original (Hawkins-Craemer-Truhlar) effective radii
        const Tcalc atomi_inv_radius = v_one / (isr.pb_radii[i] - isr.gb_offset);
        const Tcalc psi_i = psi[i];
        Tcalc egbi;
        if (tforce_is_sgnint) {
          egbi = v_one / (atomi_inv_radius + (psi_i / force_factor));
        }
        else {
          egbi = v_one / (atomi_inv_radius + psi_i);
        }
        if (egbi < 0.0) {
          egbi = 30.0;
        }
        effective_gb_radii[i] = (tforce_is_sgnint) ? llround(egbi * force_factor) : egbi;
      }
      break;
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      {
        // "GBAO" formulas
        const Tcalc atomi_radius = isr.pb_radii[i] - isr.gb_offset;
        const Tcalc atomi_inv_radius = v_one / atomi_radius;
        const Tcalc psi_i = psi[i];        
        const Tcalc fipsi = (tforce_is_sgnint) ? (psi_i / force_factor) * (-atomi_radius) :
                                                 psi_i * (-atomi_radius);
        Tcalc egbi;
        if (tcalc_is_double) {
          egbi = v_one / (atomi_inv_radius -
                          tanh((isr.gb_alpha[i] - (isr.gb_beta[i] * fipsi) +
                                (isr.gb_gamma[i] * fipsi * fipsi)) * fipsi) / 
                          isr.pb_radii[i]);
        }
        else {
          egbi = v_one / (atomi_inv_radius -
                          tanhf((isr.gb_alpha[i] - (isr.gb_beta[i] * fipsi) +
                                 (isr.gb_gamma[i] * fipsi * fipsi)) * fipsi) / 
                          isr.pb_radii[i]);
        }
        effective_gb_radii[i] = (tforce_is_sgnint) ? llround(egbi * force_factor) : egbi;
      }
      break;
    case ImplicitSolventModel::NONE:
      break;
    }
  }
  
  // Compute inherent Generalized Born energies and initialize an array for solvent forces
  for (int i = 0; i < nbk.natom; i++) {
    const Tcalc atomi_q = nbk.charge[i];
    const Tcalc egbi = effective_gb_radii[i];
    const Tcalc atomi_radius = (tforce_is_sgnint) ? egbi / force_factor : egbi;
    const Tcalc expmkf = (tcalc_is_double) ?
                         exp(-default_gb_kscale * isr.kappa * atomi_radius) / isr.dielectric :
                         expf(-default_gb_kscale * isr.kappa * atomi_radius) / isr.dielectric;
    const Tcalc dielfac = v_one - expmkf;
    const Tcalc atmq2h = v_half * atomi_q * atomi_q * nbk.coulomb_constant;
    const Tcalc atmqd2h = atmq2h * dielfac;
    const Tcalc contrib = -atmqd2h / atomi_radius;
    egb_energy += contrib;
    egb_acc += static_cast<llint>(contrib * nrg_scale_factor);
    if (eval_force == EvaluateForce::YES) {
      if (tforce_is_sgnint) {
        sumdeijda[i] = llround((atmqd2h - (default_gb_kscale * isr.kappa *
                                           atmq2h * expmkf * atomi_radius)) * force_factor);
      }
      else {
        sumdeijda[i] = atmqd2h - (default_gb_kscale * isr.kappa * atmq2h * expmkf * atomi_radius);
      }
    }
  }
  
  // Due to the lack of exclusions, the Generalized Born reference calculation is a much simpler
  // pair of nested loops over all atoms without self-interactions or double-counting.  However,
  // the tiling continues to add a degree of complexity.
  for (int sti = 0; sti < ser.supertile_stride_count; sti++) {
    for (int stj = 0; stj <= sti; stj++) {

      // Tile dimensions and locations
      const int stni_atoms = std::min(ser.natom - (supertile_length * sti), supertile_length);
      const int stnj_atoms = std::min(ser.natom - (supertile_length * stj), supertile_length);
      const int ni_tiles = (stni_atoms + tile_length - 1) / tile_length;
      const int nj_tiles = (stnj_atoms + tile_length - 1) / tile_length;

      // Access the supertile's map index: if zero, there are no exclusions to worry about
      const int diag_supertile = (sti == stj);

      // The outer loops can proceed until the branch about exclusions
      for (int ti = 0; ti < ni_tiles; ti++) {
        const int ni_atoms = std::min(stni_atoms - (ti * tile_length), tile_length);
        const int tjlim = (nj_tiles * (1 - diag_supertile)) + (diag_supertile * ti);
        for (int tj = 0; tj <= tjlim; tj++) {
          const int nj_atoms = std::min(stnj_atoms - (tj * tile_length), tile_length);
          const int diag_tile = diag_supertile * (ti == tj);

          // Pre-cache the atom positions to avoid conversions in the inner loops.  Pre-cache
          // properties and initialize local force accumulators to mimic GPU activity.
          for (int i = 0; i < ni_atoms; i++) {
            const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
            if (tcoord_is_sgnint) {
              cachi_xcrd[i] = static_cast<Tcalc>(xcrd[atom_i]) * inv_gpos_factor;
              cachi_ycrd[i] = static_cast<Tcalc>(ycrd[atom_i]) * inv_gpos_factor;
              cachi_zcrd[i] = static_cast<Tcalc>(zcrd[atom_i]) * inv_gpos_factor;
            }
            else {
              cachi_xcrd[i] = xcrd[atom_i];
              cachi_ycrd[i] = ycrd[atom_i];
              cachi_zcrd[i] = zcrd[atom_i];
            }
            cachi_xfrc[i] = 0.0;
            cachi_yfrc[i] = 0.0;
            cachi_zfrc[i] = 0.0;

            // Re-use the cachi_radii array for the effective, no longer the offset intrinsic,
            // Born radii.  Re-use the cachi_screen array for atomic partial charges.  Re-use
            // the cachi_psi array to accumulate sumdeijda.
            cachi_radii[i]  = (tforce_is_sgnint) ?
                              static_cast<Tcalc>(effective_gb_radii[atom_i]) / force_factor :
                              effective_gb_radii[atom_i];
            cachi_screen[i] = nbk.charge[atom_i];
            cachi_psi[i] = 0.0;
          }
          for (int j = 0; j < nj_atoms; j++) {
            const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
            if (tcoord_is_sgnint) {
              cachj_xcrd[j] = static_cast<Tcalc>(xcrd[atom_j]) * inv_gpos_factor;
              cachj_ycrd[j] = static_cast<Tcalc>(ycrd[atom_j]) * inv_gpos_factor;
              cachj_zcrd[j] = static_cast<Tcalc>(zcrd[atom_j]) * inv_gpos_factor;
            }
            else {
              cachj_xcrd[j] = xcrd[atom_j];
              cachj_ycrd[j] = ycrd[atom_j];
              cachj_zcrd[j] = zcrd[atom_j];
            }
            cachj_xfrc[j] = 0.0;
            cachj_yfrc[j] = 0.0;
            cachj_zfrc[j] = 0.0;

            // See above for the re-use of these local detail caches and accumulators.
            cachj_radii[j]  = (tforce_is_sgnint) ?
                              static_cast<Tcalc>(effective_gb_radii[atom_j]) / force_factor :
                              effective_gb_radii[atom_j];
            cachj_screen[j] = nbk.charge[atom_j];
            cachj_psi[j] = 0.0;
          }

          // Perform a nested loop over all atoms in the tile to compute the energy
          // and the bulk of the force.
          for (int i = 0; i < ni_atoms; i++) {
            const Tcalc atomi_x = cachi_xcrd[i];
            const Tcalc atomi_y = cachi_ycrd[i];
            const Tcalc atomi_z = cachi_zcrd[i];
            const Tcalc atomi_q = nbk.coulomb_constant * cachi_screen[i];
            const Tcalc atomi_radius = cachi_radii[i];
            const int jlim = (nj_atoms * (1 - diag_tile)) + (diag_tile * i);
            Tcalc contrib = 0.0;
            for (int j = 0; j < jlim; j++) {
              const Tcalc dx = cachj_xcrd[j] - atomi_x;
              const Tcalc dy = cachj_ycrd[j] - atomi_y;
              const Tcalc dz = cachj_zcrd[j] - atomi_z;
              const Tcalc r2 = (dx * dx) + (dy * dy) + (dz * dz);
              const Tcalc qiqj = atomi_q * cachj_screen[j];
              const Tcalc ij_born_radius = atomi_radius * cachj_radii[j];
              Tcalc efac, fgbi, fgbk, expmkf;
              if (tcalc_is_double) {
                efac = exp(-r2 / (v_four * ij_born_radius));
                fgbi = v_one / sqrt(r2 + (ij_born_radius * efac));
                fgbk = -isr.kappa * default_gb_kscale / fgbi;
                expmkf = exp(fgbk) / isr.dielectric;
              }
              else {
                efac = expf(-r2 / (v_four * ij_born_radius));
                fgbi = v_one / sqrtf(r2 + (ij_born_radius * efac));
                fgbk = -isr.kappa * default_gb_kscale / fgbi;
                expmkf = expf(fgbk) / isr.dielectric;
              }
              const Tcalc dielfac = v_one - expmkf;
              contrib -= qiqj * dielfac * fgbi;
              if (eval_force == EvaluateForce::YES) {
                const Tcalc temp4 = fgbi * fgbi * fgbi;
                const Tcalc temp6 = qiqj * temp4 * (dielfac + (fgbk * expmkf));
                const Tcalc fmag = temp6 * (v_one - (v_qrtr * efac));
                const Tcalc temp5 = v_half * efac * temp6 * (ij_born_radius + (v_qrtr * r2));
                cachi_psi[i] += atomi_radius * temp5;
                cachj_psi[j] += cachj_radii[j] * temp5;
                cachi_xfrc[i] += fmag * dx;
                cachi_yfrc[i] += fmag * dy;
                cachi_zfrc[i] += fmag * dz;
                cachj_xfrc[j] -= fmag * dx;
                cachj_yfrc[j] -= fmag * dy;
                cachj_zfrc[j] -= fmag * dz;
              }
            }

            // Accumulate the energy contribution for what would be one thread's work in the
            // tile on the GPU.
            egb_acc += llround(contrib * nrg_scale_factor);
            egb_energy += contrib;
          }

          // Contribute local force accumulators back to global arrays, if appropriate
          if (eval_force == EvaluateForce::YES) {
            for (int i = 0; i < ni_atoms; i++) {
              const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
              if (tforce_is_sgnint) {
                xfrc[atom_i] += llround(cachi_xfrc[i] * force_factor);
                yfrc[atom_i] += llround(cachi_yfrc[i] * force_factor);
                zfrc[atom_i] += llround(cachi_zfrc[i] * force_factor);
                sumdeijda[atom_i] += llround(cachi_psi[i] * force_factor);
              }
              else {
                xfrc[atom_i] += cachi_xfrc[i];
                yfrc[atom_i] += cachi_yfrc[i];
                zfrc[atom_i] += cachi_zfrc[i];
                sumdeijda[atom_i] += cachi_psi[i];
              }
            }
            for (int j = 0; j < nj_atoms; j++) {
              const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
              if (tforce_is_sgnint) {
                xfrc[atom_j] += llround(cachj_xfrc[j] * force_factor);
                yfrc[atom_j] += llround(cachj_yfrc[j] * force_factor);
                zfrc[atom_j] += llround(cachj_zfrc[j] * force_factor);
                sumdeijda[atom_j] += llround(cachj_psi[j] * force_factor);
              }
              else {
                xfrc[atom_j] += cachj_xfrc[j];
                yfrc[atom_j] += cachj_yfrc[j];
                zfrc[atom_j] += cachj_zfrc[j];
                sumdeijda[atom_j] += cachj_psi[j];
              }
            }
          }
        }
      }
    }
  }

  // Contribute results
  ecard->contribute(StateVariable::GENERALIZED_BORN, egb_acc, system_index);

  // Return immediately if forces are not required, obviating the need for further indentation
  if (eval_force == EvaluateForce::NO) {
    return egb_energy;
  }

  // A third pair of nested loops over all atoms is needed to fold in derivatives of the
  // effective Born radii to the forces on each atom.  Begin by updating the energy derivative
  // factors (sumdeijda) for each atom, then roll into the nested loops.
  switch (isr.igb) {
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::NONE:
    break;
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    for (int i = 0; i < nbk.natom; i++) {
      const Tcalc atomi_radius = isr.pb_radii[i] - isr.gb_offset;
      const Tcalc psi_i = psi[i];
      const Tcalc fipsi = (tforce_is_sgnint) ? (psi_i / force_factor) * (-atomi_radius) :
                                               psi_i * (-atomi_radius);
      Tcalc thi;
      if (tcalc_is_double) {
        thi = tanh((isr.gb_alpha[i] -
                    (isr.gb_beta[i] - (isr.gb_gamma[i] * fipsi)) * fipsi) * fipsi);
      }
      else {
        thi = tanhf((isr.gb_alpha[i] -
                     (isr.gb_beta[i] - (isr.gb_gamma[i] * fipsi)) * fipsi) * fipsi);          
      }
      const Tcalc sdi_current = (tforce_is_sgnint) ?
                                static_cast<Tcalc>(sumdeijda[i]) / force_factor :
                                sumdeijda[i];
      const Tcalc sdi_multiplier = (isr.gb_alpha[i] -
                                    ((v_two * isr.gb_beta[i]) -
                                     (v_thre * isr.gb_gamma[i] * fipsi)) * fipsi) *
                                   (v_one - thi * thi) * atomi_radius / isr.pb_radii[i];
      sumdeijda[i] = (tforce_is_sgnint) ? llround(sdi_current * sdi_multiplier * force_factor) :
                                          sdi_current * sdi_multiplier;
    }
    break;
  }  
  for (int sti = 0; sti < ser.supertile_stride_count; sti++) {
    for (int stj = 0; stj <= sti; stj++) {

      // Tile dimensions and locations
      const int stni_atoms = std::min(ser.natom - (supertile_length * sti), supertile_length);
      const int stnj_atoms = std::min(ser.natom - (supertile_length * stj), supertile_length);
      const int ni_tiles = (stni_atoms + tile_length - 1) / tile_length;
      const int nj_tiles = (stnj_atoms + tile_length - 1) / tile_length;

      // Access the supertile's map index: if zero, there are no exclusions to worry about
      const int diag_supertile = (sti == stj);

      // The outer loops can proceed until the branch about exclusions
      for (int ti = 0; ti < ni_tiles; ti++) {
        const int ni_atoms = std::min(stni_atoms - (ti * tile_length), tile_length);
        const int tjlim = (nj_tiles * (1 - diag_supertile)) + (diag_supertile * ti);
        for (int tj = 0; tj <= tjlim; tj++) {
          const int nj_atoms = std::min(stnj_atoms - (tj * tile_length), tile_length);
          const int diag_tile = diag_supertile * (ti == tj);

          // Pre-cache the atom positions to avoid conversions in the inner loops.  Pre-cache
          // properties and initialize local force accumulators to mimic GPU activity.
          for (int i = 0; i < ni_atoms; i++) {
            const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
            if (tcoord_is_sgnint) {
              cachi_xcrd[i] = static_cast<Tcalc>(xcrd[atom_i]) * inv_gpos_factor;
              cachi_ycrd[i] = static_cast<Tcalc>(ycrd[atom_i]) * inv_gpos_factor;
              cachi_zcrd[i] = static_cast<Tcalc>(zcrd[atom_i]) * inv_gpos_factor;
            }
            else {
              cachi_xcrd[i] = xcrd[atom_i];
              cachi_ycrd[i] = ycrd[atom_i];
              cachi_zcrd[i] = zcrd[atom_i];
            }
            cachi_xfrc[i] = 0.0;
            cachi_yfrc[i] = 0.0;
            cachi_zfrc[i] = 0.0;

            // The cachi_radii array returns to its original role caching the offset intrinsic
            // Born radii and the cachi_screen array caches the GB screening factors.  Neck GB
            // factors return to their cache home.  However, the cachi_psi array handles
            // sum_deijda once more.
            cachi_radii[i]  = isr.pb_radii[atom_i] - isr.gb_offset;
            cachi_screen[i] = isr.gb_screen[atom_i];
            if (isr.igb == ImplicitSolventModel::NECK_GB ||
                isr.igb == ImplicitSolventModel::NECK_GB_II) {
              cachi_neck_idx[i] = isr.neck_gb_idx[atom_i];
            }
            cachi_psi[i] = (tforce_is_sgnint) ?
                           static_cast<Tcalc>(sumdeijda[atom_i]) / force_factor :
                           sumdeijda[atom_i];
          }
          for (int j = 0; j < nj_atoms; j++) {
            const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
            if (tcoord_is_sgnint) {
              cachj_xcrd[j] = static_cast<Tcalc>(xcrd[atom_j]) * inv_gpos_factor;
              cachj_ycrd[j] = static_cast<Tcalc>(ycrd[atom_j]) * inv_gpos_factor;
              cachj_zcrd[j] = static_cast<Tcalc>(zcrd[atom_j]) * inv_gpos_factor;
            }
            else {
              cachj_xcrd[j] = xcrd[atom_j];
              cachj_ycrd[j] = ycrd[atom_j];
              cachj_zcrd[j] = zcrd[atom_j];
            }
            cachj_xfrc[j] = 0.0;
            cachj_yfrc[j] = 0.0;
            cachj_zfrc[j] = 0.0;

            // See above for the re-use of these local detail caches and accumulators.
            cachj_radii[j]  = isr.pb_radii[atom_j] - isr.gb_offset;
            cachj_screen[j] = isr.gb_screen[atom_j];
            if (isr.igb == ImplicitSolventModel::NECK_GB ||
                isr.igb == ImplicitSolventModel::NECK_GB_II) {
              cachj_neck_idx[j] = isr.neck_gb_idx[atom_j];
            }
            cachj_psi[j] = (tforce_is_sgnint) ?
                           static_cast<Tcalc>(sumdeijda[atom_j]) / force_factor :
                           sumdeijda[atom_j];
          }
          for (int i = 0; i < ni_atoms; i++) {
            const Tcalc atomi_x = cachi_xcrd[i];
            const Tcalc atomi_y = cachi_ycrd[i];
            const Tcalc atomi_z = cachi_zcrd[i];
            const Tcalc atomi_radius = cachi_radii[i];
            const Tcalc atomi_inv_radius = v_one / atomi_radius;
            const int jlim = (nj_atoms * (1 - diag_tile)) + (diag_tile * i);
            for (int j = 0; j < jlim; j++) {
              const Tcalc dx = cachj_xcrd[j] - atomi_x;
              const Tcalc dy = cachj_ycrd[j] - atomi_y;
              const Tcalc dz = cachj_zcrd[j] - atomi_z;
              const Tcalc atomj_radius = cachj_radii[j];
              const Tcalc atomj_inv_radius = v_one / atomj_radius;
              const Tcalc r2 = (dx * dx) + (dy * dy) + (dz * dz);
              const Tcalc invr = (tcalc_is_double) ? v_one / sqrt(r2) : v_one / sqrtf(r2);
              const Tcalc invr2 = invr * invr;
              const Tcalc r = r2 * invr;

              // First computation: atom I -> atom J
              const Tcalc sj = cachj_screen[j] * atomj_radius;
              const Tcalc sj2 = sj * sj;
              Tcalc datmpi, datmpj;
              if (r > v_four * sj) {
                const Tcalc tmpsd  = sj2 * invr2;
                const Tcalc dumbo  = gte +
                                     tmpsd * (gtf + tmpsd * (gtg + tmpsd * (gth + tmpsd * gthh)));
                datmpi = tmpsd * sj * invr2 * invr2 * dumbo;
              }
              else if (r > atomi_radius + sj) {
                const Tcalc temp1  = v_one / (r2 - sj2);
                if (tcalc_is_double) {
                  datmpi = (temp1 * sj * (-v_half * invr2 + temp1)) +
                           (v_qrtr * invr * invr2 * log((r - sj) / (r + sj)));
                }
                else {
                  datmpi = (temp1 * sj * (-v_half * invr2 + temp1)) +
                           (v_qrtr * invr * invr2 * logf((r - sj) / (r + sj)));
                }
              }
              else if (r > fabs(atomi_radius - sj)) {
                const Tcalc temp1  = v_one / (r + sj);
                const Tcalc invr3  = invr2 * invr;
                if (tcalc_is_double) {
                  datmpi = -v_qrtr * ((-v_half * (r2 - atomi_radius * atomi_radius + sj2) *
                                       invr3 * atomi_inv_radius * atomi_inv_radius) +
                                      (invr * temp1 * (temp1 - invr)) -
                                      (invr3 * log(atomi_radius * temp1)));
                }
                else {
                  datmpi = -v_qrtr * ((-v_half * (r2 - atomi_radius * atomi_radius + sj2) *
                                       invr3 * atomi_inv_radius * atomi_inv_radius) +
                                      (invr * temp1 * (temp1 - invr)) -
                                      (invr3 * logf(atomi_radius * temp1)));
                }
              }
              else if (atomi_radius < sj) {
                const Tcalc temp1  = v_one / (r2 - sj2);
                if (tcalc_is_double) {
                  datmpi = -v_half * ((sj * invr2 * temp1) - (v_two * sj * temp1 * temp1) -
                                      (v_half * invr2 * invr * log((sj - r) / (sj + r))));
                }
                else {
                  datmpi = -v_half * ((sj * invr2 * temp1) - (v_two * sj * temp1 * temp1) -
                                      (v_half * invr2 * invr * logf((sj - r) / (sj + r))));
                }
              }
              else {
                datmpi = v_zero;
              }

              // Second computation: atom J -> atom I
              const Tcalc si = cachi_screen[i] * atomi_radius;
              const Tcalc si2 = si * si;
              if (r > v_four * si) {
                const Tcalc tmpsd  = si2 * invr2;
                const Tcalc dumbo  = gte +
                                     tmpsd * (gtf + tmpsd * (gtg + tmpsd * (gth + tmpsd * gthh)));
                datmpj = tmpsd * si * invr2 * invr2 * dumbo;
              }
              else if (r > atomj_radius + si) {
                const Tcalc temp1  = v_one / (r2 - si2);
                if (tcalc_is_double) {
                  datmpj = (temp1 * si * (-v_half * invr2 + temp1)) +
                           (v_qrtr * invr * invr2 * log((r - si) / (r + si)));
                }
                else {
                  datmpj = (temp1 * si * (-v_half * invr2 + temp1)) +
                           (v_qrtr * invr * invr2 * logf((r - si) / (r + si)));
                }
              }
              else if (r > fabs(atomj_radius - si)) {
                const Tcalc temp1 = v_one / (r + si);
                const Tcalc invr3 = invr2 * invr;
                if (tcalc_is_double) {
                  datmpj = -v_qrtr * ((-v_half * (r2 - atomj_radius * atomj_radius + si2) *
                                       invr3 * atomj_inv_radius * atomj_inv_radius) +
                                      (invr * temp1 * (temp1 - invr)) -
                                      (invr3 * log(atomj_radius * temp1)));
                }
                else {
                  datmpj = -v_qrtr * ((-v_half * (r2 - atomj_radius * atomj_radius + si2) *
                                       invr3 * atomj_inv_radius * atomj_inv_radius) +
                                      (invr * temp1 * (temp1 - invr)) -
                                      (invr3 * logf(atomj_radius * temp1)));
                }
              }
              else if (atomj_radius < si) {
                const Tcalc temp1  = v_one / (r2 - si2);
                if (tcalc_is_double) {
                  datmpj = -v_half * ((si * invr2 * temp1) - (v_two * si * temp1 * temp1) -
                                      (v_half * invr2 * invr * log((si - r) / (si + r))));
                }
                else {
                  datmpj = -v_half * ((si * invr2 * temp1) - (v_two * si * temp1 * temp1) -
                                      (v_half * invr2 * invr * logf((si - r) / (si + r))));
                }
              }
              else {
                datmpj = v_zero;
              }

              // Neck GB contributions
              if ((isr.igb == ImplicitSolventModel::NECK_GB ||
                   isr.igb == ImplicitSolventModel::NECK_GB_II) &&
                  r < cachi_radii[i] + cachj_radii[j] + (v_two * isr.gb_offset) + isr.gb_neckcut) {

                // First computation: atom I -> atom J
                const int ij_table_idx = (isr.table_size * cachj_neck_idx[j]) +
                                         cachi_neck_idx[i];
                Tcalc mdist = r - isr.neck_max_sep[ij_table_idx];
                Tcalc mdist2 = mdist * mdist;
                Tcalc mdist6 = mdist2 * mdist2 * mdist2;
                Tcalc temp1 = v_one + mdist2 + (v_pthr * mdist6);
                temp1 = temp1 * temp1 * r;
                datmpi += (((v_two * mdist) + (v_opei * mdist2 * mdist2 * mdist)) *
                           isr.neck_max_val[ij_table_idx] * isr.gb_neckscale) / temp1;

                // Second computation: atom J -> atom I
                const int ji_table_idx = (isr.table_size * cachi_neck_idx[i]) +
                                         cachj_neck_idx[j];
                mdist = r - isr.neck_max_sep[ji_table_idx];
                mdist2 = mdist * mdist;
                mdist6 = mdist2 * mdist2 * mdist2;
                temp1 = v_one + mdist2 + (v_pthr * mdist6);
                temp1 = temp1 * temp1 * r;
                datmpj += (((v_two * mdist) + (v_opei * mdist2 * mdist2 * mdist)) *
                           isr.neck_max_val[ji_table_idx] * isr.gb_neckscale) / temp1;
              }

              // Contribute the derivatives to the force arrays
              const Tcalc fmag = (datmpi * cachi_psi[i]) + (datmpj * cachj_psi[j]);
              cachi_xfrc[i] -= fmag * dx;
              cachi_yfrc[i] -= fmag * dy;
              cachi_zfrc[i] -= fmag * dz;
              cachj_xfrc[j] += fmag * dx;
              cachj_yfrc[j] += fmag * dy;
              cachj_zfrc[j] += fmag * dz;
            }
          }

          // Contribute local force accumulators back to the global arrays
          for (int i = 0; i < ni_atoms; i++) {
            const int atom_i = i + (ti * tile_length) + (sti * supertile_length);
            if (tforce_is_sgnint) {
              xfrc[atom_i] += llround(cachi_xfrc[i] * force_factor);
              yfrc[atom_i] += llround(cachi_yfrc[i] * force_factor);
              zfrc[atom_i] += llround(cachi_zfrc[i] * force_factor);
            }
            else {
              xfrc[atom_i] += cachi_xfrc[i];
              yfrc[atom_i] += cachi_yfrc[i];
              zfrc[atom_i] += cachi_zfrc[i];
            }
          }
          for (int j = 0; j < nj_atoms; j++) {
            const int atom_j = j + (tj * tile_length) + (stj * supertile_length);
            if (tforce_is_sgnint) {
              xfrc[atom_j] += llround(cachj_xfrc[j] * force_factor);
              yfrc[atom_j] += llround(cachj_yfrc[j] * force_factor);
              zfrc[atom_j] += llround(cachj_zfrc[j] * force_factor);
            }
            else {
              xfrc[atom_j] += cachj_xfrc[j];
              yfrc[atom_j] += cachj_yfrc[j];
              zfrc[atom_j] += cachj_zfrc[j];
            }
          }
        }
      }
    }
  }

  // Return the double-precision energy sum, if of interest
  return egb_energy;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double evaluateGeneralizedBornEnergy(const NonbondedKit<Tcalc> nbk,
                                     const StaticExclusionMaskReader ser,
                                     const ImplicitSolventKit<Tcalc> isk,
                                     const NeckGeneralizedBornKit<Tcalc> ngb_kit,
                                     const CoordinateSeriesReader<Tcoord> csr, ScoreCard *ecard,
                                     const int system_index, const int force_scale_bits) {
  const size_t atom_os = static_cast<size_t>(system_index) *
                         roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t xfrm_os = static_cast<size_t>(system_index) * roundUp<size_t>(9, warp_size_zu);

  // Compute the force scaling factor based on the requested bit count.  While this routine will
  // not compute forces per se, the scaling factor is still relevant to the accumulators for Born
  // radii.
  const Tcalc force_scale = (isSignedIntegralScalarType<Tcoord>()) ?
                            pow(2.0, force_scale_bits) : 1.0;
  std::vector<Tcoord> effective_gb_radii(csr.natom);
  std::vector<Tcoord> psi(csr.natom);
  std::vector<Tcoord> sumdeijda(csr.natom);
  return evaluateGeneralizedBornEnergy<Tcoord,
                                       Tcoord, Tcalc>(nbk, ser, isk, ngb_kit,
                                                      &csr.xcrd[atom_os], &csr.ycrd[atom_os],
                                                      &csr.zcrd[atom_os], nullptr,
                                                      nullptr, nullptr, effective_gb_radii.data(),
                                                      psi.data(), sumdeijda.data(), ecard,
                                                      EvaluateForce::NO, system_index,
                                                      csr.inv_gpos_scale, force_scale);
}

} // namespace energy
} // namespace stormm
