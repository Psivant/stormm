// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void rattlePositions(Tcoord* xdev, Tcoord* ydev, Tcoord* zdev, Tcoord* xvel_dev, Tcoord* yvel_dev,
                     Tcoord* zvel_dev, const Tcoord* xref, const Tcoord* yref, const Tcoord* zref,
                     const ConstraintKit<Tcalc> &cnk, const Tcalc dt, const Tcalc tol,
                     const int max_iter, const RattleMethod style, const Tcalc gpos_scale_factor,
                     const Tcalc vel_scale_factor) {
  const bool tcoord_is_real  = (isSignedIntegralScalarType<Tcoord>() == false);
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const Tcoord zero = 0.0;
  const Tcalc inv_gpos_scale_factor = static_cast<Tcalc>(1.0) / gpos_scale_factor;
  const Tcalc vel_scale_dt = (1.0 / dt) * vel_scale_factor * inv_gpos_scale_factor;
  std::vector<Tcoord> xdev_init, ydev_init, zdev_init;
  for (int i = 0; i < cnk.ngroup; i++) {
    const int atom_llim = cnk.group_bounds[i];
    const int atom_hlim = cnk.group_bounds[i + 1];
    const int nbnd = atom_hlim - atom_llim - 1;
    const int catom_idx = cnk.group_list[atom_llim];
    const int param_idx = cnk.group_param_idx[i];
    const int param_llim = cnk.group_param_bounds[param_idx];
    const Tcalc catom_inv_mass = cnk.group_inv_masses[param_llim];

    // Record the original positions of all particles.  This is critical for the post-hoc velocity
    // update accounting for the move made by the constraints.
    xdev_init.resize(nbnd + 1);
    ydev_init.resize(nbnd + 1);
    zdev_init.resize(nbnd + 1);
    for (int j = 0; j <= nbnd; j++) {
      const size_t jatom = cnk.group_list[atom_llim + j];
      xdev_init[j] = xdev[jatom];
      ydev_init[j] = ydev[jatom];
      zdev_init[j] = zdev[jatom];
    }

    // Iterate to converge the geometry towards correct bond lengths.
    int iter = 0;
    bool done = false;
    while ((! done) && iter < max_iter) {
      Tcoord cmove_x = zero;
      Tcoord cmove_y = zero;
      Tcoord cmove_z = zero;
      done = true;
      for (int j = 0; j < nbnd; j++) {
        const int hatom_idx = cnk.group_list[atom_llim + j + 1];

        // Compute the displacement after the most recent update (this is expected to be a state
        // for which construction is in progress)
        Tcalc dx, dy, dz;
        if (tcoord_is_real) {
          dx = xdev[hatom_idx] - xdev[catom_idx];
          dy = ydev[hatom_idx] - ydev[catom_idx];
          dz = zdev[hatom_idx] - zdev[catom_idx];
        }
        else {
          dx = static_cast<Tcalc>(xdev[hatom_idx] - xdev[catom_idx]) * inv_gpos_scale_factor;
          dy = static_cast<Tcalc>(ydev[hatom_idx] - ydev[catom_idx]) * inv_gpos_scale_factor;
          dz = static_cast<Tcalc>(zdev[hatom_idx] - zdev[catom_idx]) * inv_gpos_scale_factor;
        }
        const Tcalc r2 = (dx * dx) + (dy * dy) + (dz * dz);
        const Tcalc l2_target = cnk.group_sq_lengths[param_llim + j + 1];
        const Tcalc delta = l2_target - r2;
        if (fabs(delta) > tol) {
          done = false;
          Tcalc dx_ref, dy_ref, dz_ref;
          if (tcoord_is_real) {
            dx_ref = xref[hatom_idx] - xref[catom_idx];
            dy_ref = yref[hatom_idx] - yref[catom_idx];
            dz_ref = zref[hatom_idx] - zref[catom_idx];
          }
          else {
            dx_ref = static_cast<Tcalc>(xref[hatom_idx] - xref[catom_idx]) * inv_gpos_scale_factor;
            dy_ref = static_cast<Tcalc>(yref[hatom_idx] - yref[catom_idx]) * inv_gpos_scale_factor;
            dz_ref = static_cast<Tcalc>(zref[hatom_idx] - zref[catom_idx]) * inv_gpos_scale_factor;
          }
          const Tcalc dot_disp = (dx * dx_ref) + (dy * dy_ref) + (dz * dz_ref);
          const Tcalc hatom_inv_mass = cnk.group_inv_masses[param_llim + j + 1];
          const Tcalc term = (tcalc_is_double) ?
                             0.6  * delta / (dot_disp * (catom_inv_mass + hatom_inv_mass)) :
                             0.6f * delta / (dot_disp * (catom_inv_mass + hatom_inv_mass));
          if (tcoord_is_real) {
            const Tcalc ca_term = term * catom_inv_mass;
            const Tcalc ha_term = term * hatom_inv_mass;
            switch (style) {
            case RattleMethod::SEQUENTIAL:
              xdev[catom_idx] -= static_cast<Tcoord>(dx_ref * ca_term);
              ydev[catom_idx] -= static_cast<Tcoord>(dy_ref * ca_term);
              zdev[catom_idx] -= static_cast<Tcoord>(dz_ref * ca_term);
              xdev[hatom_idx] += static_cast<Tcoord>(dx_ref * ha_term);
              ydev[hatom_idx] += static_cast<Tcoord>(dy_ref * ha_term);
              zdev[hatom_idx] += static_cast<Tcoord>(dz_ref * ha_term);
              break;
            case RattleMethod::CENTER_SUM:
              cmove_x -= static_cast<Tcoord>(dx_ref * ca_term);
              cmove_y -= static_cast<Tcoord>(dy_ref * ca_term);
              cmove_z -= static_cast<Tcoord>(dz_ref * ca_term);
              xdev[hatom_idx] += static_cast<Tcoord>(dx_ref * ha_term);
              ydev[hatom_idx] += static_cast<Tcoord>(dy_ref * ha_term);
              zdev[hatom_idx] += static_cast<Tcoord>(dz_ref * ha_term);
              break;
            }
          }
          else {
            const Tcalc ca_term = term * catom_inv_mass * gpos_scale_factor;
            const Tcalc ha_term = term * hatom_inv_mass * gpos_scale_factor;
            switch (style) {
            case RattleMethod::SEQUENTIAL:
              xdev[catom_idx] -= llround(dx_ref * ca_term);
              ydev[catom_idx] -= llround(dy_ref * ca_term);
              zdev[catom_idx] -= llround(dz_ref * ca_term);
              xdev[hatom_idx] += llround(dx_ref * ha_term);
              ydev[hatom_idx] += llround(dy_ref * ha_term);
              zdev[hatom_idx] += llround(dz_ref * ha_term);
              break;
            case RattleMethod::CENTER_SUM:
              cmove_x -= llround(dx_ref * ca_term);
              cmove_y -= llround(dy_ref * ca_term);
              cmove_z -= llround(dz_ref * ca_term);
              xdev[hatom_idx] += llround(dx_ref * ha_term);
              ydev[hatom_idx] += llround(dy_ref * ha_term);
              zdev[hatom_idx] += llround(dz_ref * ha_term);
              break;
            }
          }
        }
      }

      // If the central atom's motion was accrued over all motions of the distal atoms, adjust its
      // position now.
      switch (style) {
      case RattleMethod::SEQUENTIAL:
        break;
      case RattleMethod::CENTER_SUM:
        xdev[catom_idx] += cmove_x;
        ydev[catom_idx] += cmove_y;
        zdev[catom_idx] += cmove_z;
        break;
      }
      iter++;
    }
    if (iter == max_iter && (! done)) {
      rtWarn("Maximum iteration count (" + std::to_string(iter) + ") was reached attempting to "
             "RATTLE a bond group centered at atom " + std::to_string(catom_idx) + " (" +
             std::to_string(nbnd) + " bonds in all).", "rattlePositions");
    }

    // Apply the post-hoc velocity correction.  The particles moved by the amount of the correction
    // in the space of one time step.
    for (int j = 0; j <= nbnd; j++) {
      const size_t jatom = cnk.group_list[atom_llim + j];
      Tcoord dxvel, dyvel, dzvel;
      const Tcalc xmove = xdev[jatom] - xdev_init[j];
      const Tcalc ymove = ydev[jatom] - ydev_init[j];
      const Tcalc zmove = zdev[jatom] - zdev_init[j];
      if (tcoord_is_real) {
        dxvel = xmove * vel_scale_dt;
        dyvel = ymove * vel_scale_dt;
        dzvel = zmove * vel_scale_dt;
      }
      else {
        dxvel = llround(xmove * vel_scale_dt);
        dyvel = llround(ymove * vel_scale_dt);
        dzvel = llround(zmove * vel_scale_dt);
      }
      xvel_dev[jatom] += dxvel;
      yvel_dev[jatom] += dyvel;
      zvel_dev[jatom] += dzvel;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void rattlePositions(PhaseSpaceWriter *psw, const ConstraintKit<Tcalc> &cnk, const Tcalc dt,
                     const Tcalc tol, const int max_iter, const RattleMethod style) {
  rattlePositions<double, Tcalc>(psw->xalt, psw->yalt, psw->zalt, psw->vxalt, psw->vyalt,
                                 psw->vzalt, psw->xcrd, psw->ycrd, psw->zcrd, cnk, dt, tol,
                                 max_iter, style);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2, typename T4>
void rattlePositions(PsSynthesisWriter *poly_psw, const SyValenceKit<T> &poly_vk,
                     const SyAtomUpdateKit<T, T2, T4> &poly_auk, const T dt, const T tol,
                     const int max_iter) {
  const bool tcalc_is_double = (std::type_index(typeid(T)).hash_code() == double_type_index);
  const T zero = 0.0;

  // Allocate arrays to hold the cached data  
  std::vector<int> imported_atom_ids(maximum_valence_work_unit_atoms);
  std::vector<llint> imported_xcrd(maximum_valence_work_unit_atoms);
  std::vector<llint> imported_ycrd(maximum_valence_work_unit_atoms);
  std::vector<llint> imported_zcrd(maximum_valence_work_unit_atoms);
  std::vector<int> imported_xcrd_ovrf, imported_ycrd_ovrf, imported_zcrd_ovrf;
  if (tcalc_is_double) {
    imported_xcrd_ovrf.resize(maximum_valence_work_unit_atoms);
    imported_ycrd_ovrf.resize(maximum_valence_work_unit_atoms);
    imported_zcrd_ovrf.resize(maximum_valence_work_unit_atoms);
  }

  // Loop over all valence work units, which contain the constraint instructions
  for (int vwu_idx = 0; vwu_idx < poly_vk.nvwu; vwu_idx++) {

    // Import atoms as the work unit would on the GPU.
    const int2 import_limits = poly_vk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                     static_cast<int>(VwuAbstractMap::IMPORT)];
    for (int i = import_limits.x; i < import_limits.y; i++) {
      const int local_idx  = i - import_limits.x;
      const size_t global_idx = poly_vk.vwu_imports[i];
      imported_atom_ids[local_idx] = global_idx;
      imported_xcrd[local_idx] = poly_psw->xalt[global_idx];
      imported_ycrd[local_idx] = poly_psw->yalt[global_idx];
      imported_zcrd[local_idx] = poly_psw->zalt[global_idx];
      if (tcalc_is_double) {
        imported_xcrd_ovrf[local_idx] = poly_psw->xalt_ovrf[global_idx];
        imported_ycrd_ovrf[local_idx] = poly_psw->yalt_ovrf[global_idx];
        imported_zcrd_ovrf[local_idx] = poly_psw->zalt_ovrf[global_idx];
      }
    }
    
    // Loop over instruction batches, honoring warp coalescence
    const int2 insr_limits = poly_vk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +  
                                                   static_cast<int>(VwuAbstractMap::CGROUP)];
    std::vector<uint2> warp_insr(warp_size_int);
    
    // Cycle and process a full warp's worth of instructions, as if performing warp-synchronous
    // programming.  The reference coordinates will not change, and in the GPU kernel can be
    // pre-loaded into registers just as the atoms' inverse masses.  The sum of the inverse
    // masses is pre-computed as a parameter to help avoid some catastrophic cancellation.  On
    // the GPU, all of this information can be read directly from the reference abstract, without
    // the need for a cached intermediary, as it will only need to be read once: the registers
    // are the cache space.
    std::vector<T> dx_ref(warp_size_int), dy_ref(warp_size_int), dz_ref(warp_size_int);
    std::vector<T> ca_invmass(warp_size_int), ph_invmass(warp_size_int), l2_target(warp_size_int);
    std::vector<T> combined_invmass(warp_size_int), ca_move_x(warp_size_int);
    std::vector<T> ca_move_y(warp_size_int), ca_move_z(warp_size_int);  
    std::vector<int95_t> ica_move_x(warp_size_int), ica_move_y(warp_size_int);
    std::vector<int95_t> ica_move_z(warp_size_int);
    for (int i = insr_limits.x; i < insr_limits.y; i += warp_size_int) {
      const int batch_size = std::min(warp_size_int, insr_limits.y - i);
      for (int lane = 0; lane < batch_size; lane++) {
        const uint2 tinsr = poly_auk.cnst_insr[i + lane];
        const int central_atom = (tinsr.x & 0x3ff);
        const int peripheral_atom = ((tinsr.x >> 10) & 0x3ff);
        const size_t ca_global_idx = imported_atom_ids[central_atom];
        const size_t ph_global_idx = imported_atom_ids[peripheral_atom];

        // Check that the constraint is valid
        if (central_atom == 0 && peripheral_atom == 0) {
          continue;
        }
        
        // Compute the reference displacement
        if (tcalc_is_double) {
          const int95_t idx_ref = hostInt95Sum(poly_psw->xcrd[ph_global_idx],
                                               poly_psw->xcrd_ovrf[ph_global_idx],
                                               -(poly_psw->xcrd[ca_global_idx]),
                                               -(poly_psw->xcrd_ovrf[ca_global_idx]));
          const int95_t idy_ref = hostInt95Sum(poly_psw->ycrd[ph_global_idx],
                                               poly_psw->ycrd_ovrf[ph_global_idx],
                                               -(poly_psw->ycrd[ca_global_idx]),
                                               -(poly_psw->ycrd_ovrf[ca_global_idx]));
          const int95_t idz_ref = hostInt95Sum(poly_psw->zcrd[ph_global_idx],
                                               poly_psw->zcrd_ovrf[ph_global_idx],
                                               -(poly_psw->zcrd[ca_global_idx]),
                                               -(poly_psw->zcrd_ovrf[ca_global_idx]));
          dx_ref[lane] = hostInt95ToDouble(idx_ref) * poly_psw->inv_gpos_scale_f;
          dy_ref[lane] = hostInt95ToDouble(idy_ref) * poly_psw->inv_gpos_scale_f;
          dz_ref[lane] = hostInt95ToDouble(idz_ref) * poly_psw->inv_gpos_scale_f;
        }
        else {
          const llint idx_ref = poly_psw->xcrd[ph_global_idx] - poly_psw->xcrd[ca_global_idx];
          const llint idy_ref = poly_psw->ycrd[ph_global_idx] - poly_psw->ycrd[ca_global_idx];
          const llint idz_ref = poly_psw->zcrd[ph_global_idx] - poly_psw->zcrd[ca_global_idx];
          dx_ref[lane] = static_cast<T>(idx_ref) * poly_psw->inv_gpos_scale_f;
          dy_ref[lane] = static_cast<T>(idy_ref) * poly_psw->inv_gpos_scale_f;
          dz_ref[lane] = static_cast<T>(idz_ref) * poly_psw->inv_gpos_scale_f;
        }
        ca_invmass[lane]       = poly_auk.inv_masses[ca_global_idx];
        ph_invmass[lane]       = poly_auk.inv_masses[ph_global_idx];
        l2_target[lane]        = poly_auk.cnst_grp_params[tinsr.y].x;
        combined_invmass[lane] = poly_auk.cnst_grp_params[tinsr.y].y;
      }
      bool converged;
      int iter = 0;
      do {
        converged = true;
        for (int lane = 0; lane < batch_size; lane++) {
          const uint2 tinsr = poly_auk.cnst_insr[i + lane];
          const int central_atom = (tinsr.x & 0x3ff);
          const int peripheral_atom = ((tinsr.x >> 10) & 0x3ff);

          // Check that the constraint is valid
          if (central_atom == 0 && peripheral_atom == 0) {
            continue;
          }

          // Regardless of the calculation model, the coordinates are taken in based on the
          // representation in the original synthesis.
          T dx, dy, dz;
          if (tcalc_is_double) {
            const int95_t idx = hostInt95Sum(imported_xcrd[peripheral_atom],
                                             imported_xcrd_ovrf[peripheral_atom],
                                             -imported_xcrd[central_atom],
                                             -imported_xcrd_ovrf[central_atom]);
            const int95_t idy = hostInt95Sum(imported_ycrd[peripheral_atom],
                                             imported_ycrd_ovrf[peripheral_atom],
                                             -imported_ycrd[central_atom],
                                             -imported_ycrd_ovrf[central_atom]);
            const int95_t idz = hostInt95Sum(imported_zcrd[peripheral_atom],
                                             imported_zcrd_ovrf[peripheral_atom],
                                             -imported_zcrd[central_atom],
                                             -imported_zcrd_ovrf[central_atom]);
            dx = hostInt95ToDouble(idx) * poly_psw->inv_gpos_scale_f;
            dy = hostInt95ToDouble(idy) * poly_psw->inv_gpos_scale_f;
            dz = hostInt95ToDouble(idz) * poly_psw->inv_gpos_scale_f;
          } 
          else {
            dx = static_cast<T>(imported_xcrd[peripheral_atom] - imported_xcrd[central_atom]) *
                 poly_psw->inv_gpos_scale_f;
            dy = static_cast<T>(imported_ycrd[peripheral_atom] - imported_ycrd[central_atom]) *
                 poly_psw->inv_gpos_scale_f;
            dz = static_cast<T>(imported_zcrd[peripheral_atom] - imported_zcrd[central_atom]) *
                 poly_psw->inv_gpos_scale_f;
          }

          // Proceed with the RATTLE iteration
          const T r2 = (dx * dx) + (dy * dy) + (dz * dz);
          const T delta = (l2_target[lane] - r2);
          if (std::abs(delta) > tol) {
            converged = false;
            const T dot_disp = (dx * dx_ref[lane]) + (dy * dy_ref[lane]) + (dz * dz_ref[lane]);
            const T term = (tcalc_is_double) ?
                           0.6  * delta / (dot_disp * combined_invmass[lane]) :
                           0.6f * delta / (dot_disp * combined_invmass[lane]);
            ca_move_x[lane] = -dx_ref[lane] * term * ca_invmass[lane];
            ca_move_y[lane] = -dy_ref[lane] * term * ca_invmass[lane];
            ca_move_z[lane] = -dz_ref[lane] * term * ca_invmass[lane];
            const T ph_move_x = dx_ref[lane] * term * ph_invmass[lane];
            const T ph_move_y = dy_ref[lane] * term * ph_invmass[lane];
            const T ph_move_z = dz_ref[lane] * term * ph_invmass[lane];

            // Hub-and-spoke constraints assume that the peripheral atoms are connected only to
            // the central atom, not to each other.  This permits the peripheral atom's position
            // to be updated now.  In terms of 32-bit floating point numbers, nudging the atom in
            // its fixed-point representation is a more precise approach than adding the move to
            // the real number representation and then converting the entirety back to a
            // fixed-point representation.  By the definitions of the constraint groups, this can
            // be done without atomic operations.
            if (tcalc_is_double) {
              int95_t iph_move_x = hostDoubleToInt95(ph_move_x * poly_psw->gpos_scale_f);
              int95_t iph_move_y = hostDoubleToInt95(ph_move_y * poly_psw->gpos_scale_f);
              int95_t iph_move_z = hostDoubleToInt95(ph_move_z * poly_psw->gpos_scale_f);
              iph_move_x = hostSplitFPSum(iph_move_x, imported_xcrd[peripheral_atom],
                                          imported_xcrd_ovrf[peripheral_atom]);
              iph_move_y = hostSplitFPSum(iph_move_y, imported_ycrd[peripheral_atom],
                                          imported_ycrd_ovrf[peripheral_atom]);
              iph_move_z = hostSplitFPSum(iph_move_z, imported_zcrd[peripheral_atom],
                                          imported_zcrd_ovrf[peripheral_atom]);
              imported_xcrd[peripheral_atom] = iph_move_x.x;
              imported_ycrd[peripheral_atom] = iph_move_y.x;
              imported_zcrd[peripheral_atom] = iph_move_z.x;
              imported_xcrd_ovrf[peripheral_atom] = iph_move_x.y;
              imported_ycrd_ovrf[peripheral_atom] = iph_move_y.y;
              imported_zcrd_ovrf[peripheral_atom] = iph_move_z.y;
            }
            else {
              imported_xcrd[peripheral_atom] += llround(ph_move_x * poly_psw->gpos_scale_f);
              imported_ycrd[peripheral_atom] += llround(ph_move_y * poly_psw->gpos_scale_f);
              imported_zcrd[peripheral_atom] += llround(ph_move_z * poly_psw->gpos_scale_f);
            }
          }
          else {
            ca_move_x[lane] = zero;
            ca_move_y[lane] = zero;
            ca_move_z[lane] = zero;
          }
        }

        // If any constraint failed to converge, have all groups update the positions of the
        // central atom.  On the GPU, the thread with the lowest lane index handling any constraint
        // in a group will read the central atom moves of other threads handling constraints in the
        // same group.  After assembling the central atom moves, the first thread handling a
        // constraint in any given group will update the central atom's position in the imported
        // coordinates array.  On the subsequent pass, all threads will read the updated position.
        if (converged == false) {
          for (int lane = 0; lane < batch_size; lane++) {
            const uint2 tinsr = poly_auk.cnst_insr[i + lane];

            // Check that the constraint is valid.
            if ((tinsr.x & 0xfffff) > 0) {
              if (tcalc_is_double) {
                ica_move_x[lane] = hostDoubleToInt95(ca_move_x[lane] * poly_psw->gpos_scale_f);
                ica_move_y[lane] = hostDoubleToInt95(ca_move_y[lane] * poly_psw->gpos_scale_f);
                ica_move_z[lane] = hostDoubleToInt95(ca_move_z[lane] * poly_psw->gpos_scale_f);
              }
              else {
                ica_move_x[lane].x = llround(ca_move_x[lane] * poly_psw->gpos_scale_f);
                ica_move_y[lane].x = llround(ca_move_y[lane] * poly_psw->gpos_scale_f);
                ica_move_z[lane].x = llround(ca_move_z[lane] * poly_psw->gpos_scale_f);
                ica_move_x[lane].y = 0;
                ica_move_y[lane].y = 0;
                ica_move_z[lane].y = 0;
              }
            }
            else {
              ica_move_x[lane] = { 0LL, 0 };
              ica_move_y[lane] = { 0LL, 0 };
              ica_move_z[lane] = { 0LL, 0 };
            }
          }
          for (int lane = 0; lane < batch_size; lane++) {
            const uint2 tinsr = poly_auk.cnst_insr[i + lane];

            // Check that the constraint is valid.
            if ((tinsr.x & 0xfffff) == 0) {
              continue;
            }

            // Apply the moves on the central atom.
            const int cg_base_lane = ((tinsr.x >> 20) & 0xff);
            const int cg_total_lanes = ((tinsr.x >> 28) & 0xf);
            const int central_atom = (tinsr.x & 0x3ff);
            if (lane == cg_base_lane) {
              for (int j = 1; j < cg_total_lanes; j++) {
                if (tcalc_is_double) {
                  ica_move_x[lane] = hostSplitFPSum(ica_move_x[lane], ica_move_x[lane + j]);
                  ica_move_y[lane] = hostSplitFPSum(ica_move_y[lane], ica_move_y[lane + j]);
                  ica_move_z[lane] = hostSplitFPSum(ica_move_z[lane], ica_move_z[lane + j]);
                }
                else {
                  ica_move_x[lane].x += ica_move_x[lane + j].x;
                  ica_move_y[lane].x += ica_move_y[lane + j].x;
                  ica_move_z[lane].x += ica_move_z[lane + j].x;
                }
              }
              if (tcalc_is_double) {
                const int95_t ix_update = hostSplitFPSum(ica_move_x[lane],
                                                         imported_xcrd[central_atom],
                                                         imported_xcrd_ovrf[central_atom]);
                const int95_t iy_update = hostSplitFPSum(ica_move_y[lane],
                                                         imported_ycrd[central_atom],
                                                         imported_ycrd_ovrf[central_atom]);
                const int95_t iz_update = hostSplitFPSum(ica_move_z[lane],
                                                         imported_zcrd[central_atom],
                                                         imported_zcrd_ovrf[central_atom]);
                imported_xcrd[central_atom] = ix_update.x;
                imported_ycrd[central_atom] = iy_update.x;
                imported_zcrd[central_atom] = iz_update.x;
                imported_xcrd_ovrf[central_atom] = ix_update.y;
                imported_ycrd_ovrf[central_atom] = iy_update.y;
                imported_zcrd_ovrf[central_atom] = iz_update.y;
              }
              else {
                imported_xcrd[central_atom] += ica_move_x[lane].x;
                imported_ycrd[central_atom] += ica_move_y[lane].x;
                imported_zcrd[central_atom] += ica_move_z[lane].x;
              }
            }
          }
        }
        iter++;
      } while (converged == false && iter < max_iter);
    }

    // The global developing positions arrays (the "next" stage of the time cycle, which appears
    // as various "alt" arrays in the abstract) hold the positions prior to applying constraints.
    // Use these global arrays to determine how far each particle moved, then update the velocities
    // before overwriting the global developing positions arrays.  On the GPU, it will be necessary
    // to stash the newly updated, unconstrained coordinates in the force arrays prior to beginning
    // constraint iterations.  The forces are no longer needed by the time positional constraints
    // are called.  As soon as positional constraints are in place, new forces can be computed.
    const int nbits = sizeof(uint) * 8;
    const int2 mask_limits = poly_vk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                   static_cast<int>(VwuAbstractMap::MANIPULATE)];
    for (int i = import_limits.x; i < import_limits.y; i++) {
      const int mask_elem = ((i - import_limits.x) / nbits);
      const int mask_pos = (i - import_limits.x) - (mask_elem * nbits);
      if ((poly_auk.vwu_manip[mask_limits.x + mask_elem].y >> mask_pos) & 0x1) {
        const int local_idx  = i - import_limits.x;
        const size_t global_idx = poly_vk.vwu_imports[i];
        imported_atom_ids[local_idx] = global_idx;

        // Determine the implicit velocity delta.
        double mv_x, mv_y, mv_z;
        if (tcalc_is_double) {
          const int95_t imv_x = hostInt95Sum(imported_xcrd[local_idx],
                                             imported_xcrd_ovrf[local_idx],
                                             -poly_psw->xalt[global_idx],
                                             -poly_psw->xalt_ovrf[global_idx]);
          const int95_t imv_y = hostInt95Sum(imported_ycrd[local_idx],
                                             imported_ycrd_ovrf[local_idx],
                                             -poly_psw->yalt[global_idx],
                                             -poly_psw->yalt_ovrf[global_idx]);
          const int95_t imv_z = hostInt95Sum(imported_zcrd[local_idx],
                                             imported_zcrd_ovrf[local_idx],
                                             -poly_psw->zalt[global_idx],
                                             -poly_psw->zalt_ovrf[global_idx]);
          const T mv_x = hostInt95ToDouble(imv_x) * poly_psw->inv_gpos_scale_f;
          const T mv_y = hostInt95ToDouble(imv_y) * poly_psw->inv_gpos_scale_f;
          const T mv_z = hostInt95ToDouble(imv_z) * poly_psw->inv_gpos_scale_f;
          const int95_t vel_xadj = hostDoubleToInt95(mv_x / dt * poly_psw->vel_scale_f);
          const int95_t vel_yadj = hostDoubleToInt95(mv_y / dt * poly_psw->vel_scale_f);
          const int95_t vel_zadj = hostDoubleToInt95(mv_z / dt * poly_psw->vel_scale_f);
          const int95_t n_vx = hostSplitFPSum(vel_xadj, poly_psw->vxalt[global_idx],
                                              poly_psw->vxalt_ovrf[global_idx]);
          const int95_t n_vy = hostSplitFPSum(vel_yadj, poly_psw->vyalt[global_idx],
                                              poly_psw->vyalt_ovrf[global_idx]);
          const int95_t n_vz = hostSplitFPSum(vel_zadj, poly_psw->vzalt[global_idx],
                                              poly_psw->vzalt_ovrf[global_idx]);
          poly_psw->vxalt[global_idx] = n_vx.x;
          poly_psw->vyalt[global_idx] = n_vy.x;
          poly_psw->vzalt[global_idx] = n_vz.x;
          poly_psw->vxalt_ovrf[global_idx] = n_vx.y;
          poly_psw->vyalt_ovrf[global_idx] = n_vy.y;
          poly_psw->vzalt_ovrf[global_idx] = n_vz.y;
        }
        else {
          const llint imv_x = imported_xcrd[local_idx] - poly_psw->xalt[global_idx];
          const llint imv_y = imported_ycrd[local_idx] - poly_psw->yalt[global_idx];
          const llint imv_z = imported_zcrd[local_idx] - poly_psw->zalt[global_idx];
          const T mv_x = static_cast<T>(imv_x) * poly_psw->inv_gpos_scale_f;
          const T mv_y = static_cast<T>(imv_y) * poly_psw->inv_gpos_scale_f;
          const T mv_z = static_cast<T>(imv_z) * poly_psw->inv_gpos_scale_f;
          const llint vel_xadj = llround(mv_x / dt * poly_psw->vel_scale_f);
          const llint vel_yadj = llround(mv_y / dt * poly_psw->vel_scale_f);
          const llint vel_zadj = llround(mv_z / dt * poly_psw->vel_scale_f);
          poly_psw->vxalt[global_idx] += vel_xadj;
          poly_psw->vyalt[global_idx] += vel_yadj;
          poly_psw->vzalt[global_idx] += vel_zadj;
        }

        // Overwrite the global arrays
        poly_psw->xalt[global_idx] = imported_xcrd[local_idx];
        poly_psw->yalt[global_idx] = imported_ycrd[local_idx];
        poly_psw->zalt[global_idx] = imported_zcrd[local_idx];
        if (tcalc_is_double) {
          poly_psw->xalt_ovrf[global_idx] = imported_xcrd_ovrf[local_idx];
          poly_psw->yalt_ovrf[global_idx] = imported_ycrd_ovrf[local_idx];
          poly_psw->zalt_ovrf[global_idx] = imported_zcrd_ovrf[local_idx];
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void rattleVelocities(Tcoord* xvel_dev, Tcoord* yvel_dev, Tcoord *zvel_dev, const Tcoord* xcrd_ref,
                      const Tcoord* ycrd_ref, const Tcoord* zcrd_ref,
                      const ConstraintKit<Tcalc> &cnk, const Tcalc dt, const Tcalc tol,
                      const int max_iter, const RattleMethod style, const Tcalc gpos_scale_factor,
                      const Tcalc vel_scale_factor) {
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const bool tcoord_is_real = (isSignedIntegralScalarType<Tcoord>() == false);
  const Tcalc inv_gpos_scale_factor = static_cast<Tcalc>(1.0) / gpos_scale_factor;
  const Tcalc inv_vel_scale_factor = static_cast<Tcalc>(1.0) / vel_scale_factor;
  const Tcalc rtoldt = tol / dt;
  const Tcoord zero = 0.0;
  for (int i = 0; i < cnk.ngroup; i++) {
    const int atom_llim = cnk.group_bounds[i];
    const int atom_hlim = cnk.group_bounds[i + 1];
    const int nbnd = atom_hlim - atom_llim - 1;
    const int catom_idx = cnk.group_list[atom_llim];
    const int param_idx = cnk.group_param_idx[i];
    const int param_llim = cnk.group_param_bounds[param_idx];
    const Tcalc catom_inv_mass = cnk.group_inv_masses[param_llim];
    int iter = 0;
    bool done = false;
    while ((! done) && iter < max_iter) {
      Tcoord caccl_x = zero;
      Tcoord caccl_y = zero;
      Tcoord caccl_z = zero;
      done = true;
      for (int j = 0; j < nbnd; j++) {
        const int hatom_idx = cnk.group_list[atom_llim + j + 1];

        // Compute the current displacement
        Tcalc dx, dy, dz, dvx, dvy, dvz;
        if (tcoord_is_real) {
          dx = xcrd_ref[hatom_idx] - xcrd_ref[catom_idx];
          dy = ycrd_ref[hatom_idx] - ycrd_ref[catom_idx];
          dz = zcrd_ref[hatom_idx] - zcrd_ref[catom_idx];
          dvx = xvel_dev[hatom_idx] - xvel_dev[catom_idx];
          dvy = yvel_dev[hatom_idx] - yvel_dev[catom_idx];
          dvz = zvel_dev[hatom_idx] - zvel_dev[catom_idx];
        }
        else {
          dx = static_cast<Tcalc>(xcrd_ref[hatom_idx] - xcrd_ref[catom_idx]) *
               inv_gpos_scale_factor;
          dy = static_cast<Tcalc>(ycrd_ref[hatom_idx] - ycrd_ref[catom_idx]) *
               inv_gpos_scale_factor;
          dz = static_cast<Tcalc>(zcrd_ref[hatom_idx] - zcrd_ref[catom_idx]) *
               inv_gpos_scale_factor;
          dvx = static_cast<Tcalc>(xvel_dev[hatom_idx] - xvel_dev[catom_idx]) *
                inv_vel_scale_factor;
          dvy = static_cast<Tcalc>(yvel_dev[hatom_idx] - yvel_dev[catom_idx]) *
                inv_vel_scale_factor;
          dvz = static_cast<Tcalc>(zvel_dev[hatom_idx] - zvel_dev[catom_idx]) *
                inv_vel_scale_factor;
        }
        const Tcalc dot_rv = (dx * dvx) + (dy * dvy) + (dz * dvz);
        const Tcalc l2_target = cnk.group_sq_lengths[param_llim + j + 1];
        const Tcalc hatom_inv_mass = cnk.group_inv_masses[param_llim + j + 1];
        const Tcalc term = (tcalc_is_double) ?
                           -1.2  * dot_rv / ((catom_inv_mass + hatom_inv_mass) * l2_target) :
                           -1.2f * dot_rv / ((catom_inv_mass + hatom_inv_mass) * l2_target);
        if (fabs(term) > rtoldt) {
          done = false;
          if (tcoord_is_real) {
            const Tcalc ca_term = term * catom_inv_mass;
            const Tcalc ha_term = term * hatom_inv_mass;
            switch (style) {
            case RattleMethod::SEQUENTIAL:
              xvel_dev[catom_idx] -= dx * ca_term;
              yvel_dev[catom_idx] -= dy * ca_term;
              zvel_dev[catom_idx] -= dz * ca_term;
              xvel_dev[hatom_idx] += dx * ha_term;
              yvel_dev[hatom_idx] += dy * ha_term;
              zvel_dev[hatom_idx] += dz * ha_term;
              break;
            case RattleMethod::CENTER_SUM:
              caccl_x -= dx * ca_term;
              caccl_y -= dy * ca_term;
              caccl_z -= dz * ca_term;
              xvel_dev[hatom_idx] += dx * ha_term;
              yvel_dev[hatom_idx] += dy * ha_term;
              zvel_dev[hatom_idx] += dz * ha_term;
              break;
            }
          }
          else {
            const Tcalc ca_term = term * catom_inv_mass * vel_scale_factor;
            const Tcalc ha_term = term * hatom_inv_mass * vel_scale_factor;
            switch (style) {
            case RattleMethod::SEQUENTIAL:
              xvel_dev[catom_idx] -= llround(dx * ca_term);
              yvel_dev[catom_idx] -= llround(dy * ca_term);
              zvel_dev[catom_idx] -= llround(dz * ca_term);
              xvel_dev[hatom_idx] += llround(dx * ha_term);
              yvel_dev[hatom_idx] += llround(dy * ha_term);
              zvel_dev[hatom_idx] += llround(dz * ha_term);
              break;
            case RattleMethod::CENTER_SUM:
              caccl_x -= llround(dx * ca_term);
              caccl_y -= llround(dy * ca_term);
              caccl_z -= llround(dz * ca_term);
              xvel_dev[hatom_idx] += llround(dx * ha_term);
              yvel_dev[hatom_idx] += llround(dy * ha_term);
              zvel_dev[hatom_idx] += llround(dz * ha_term);
              break;
            }
          }
        }
      }  

      // If the central atom's motion was accrued over all motions of the distal atoms, adjust its
      // position now.
      switch (style) {
      case RattleMethod::SEQUENTIAL:
        break;
      case RattleMethod::CENTER_SUM:
        xvel_dev[catom_idx] += caccl_x;
        yvel_dev[catom_idx] += caccl_y;
        zvel_dev[catom_idx] += caccl_z;
        break;
      }
      iter++;
    }
    if (iter == max_iter && (! done)) {
      rtWarn("Maximum iteration count (" + std::to_string(iter) + ") was reached attempting to "
             "RATTLE a bond group centered at atom " + std::to_string(catom_idx) + " (" +
             std::to_string(nbnd) + " bonds in all).", "rattleVelocities");
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void rattleVelocities(PhaseSpaceWriter *psw, const ConstraintKit<Tcalc> &cnk, const Tcalc dt,
                      const Tcalc tol, const int max_iter, const RattleMethod style) {
  rattleVelocities<double, Tcalc>(psw->vxalt, psw->vyalt, psw->vzalt, psw->xcrd, psw->ycrd,
                                  psw->zcrd, cnk, dt, tol, max_iter, style);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T2, typename T4>
void rattleVelocities(PsSynthesisWriter *poly_psw, const SyValenceKit<T> &poly_vk,
                      const SyAtomUpdateKit<T, T2, T4> &poly_auk, const T dt, const T tol,
                      const int max_iter) {
  const bool tcalc_is_double = (std::type_index(typeid(T)).hash_code() == double_type_index);
  const T rtoldt = tol * poly_psw->vel_scale_f / dt;

  // Allocate arrays to hold the cached data  
  std::vector<int> imported_atom_ids(maximum_valence_work_unit_atoms);
  std::vector<llint> imported_xcrd(maximum_valence_work_unit_atoms);
  std::vector<llint> imported_ycrd(maximum_valence_work_unit_atoms);
  std::vector<llint> imported_zcrd(maximum_valence_work_unit_atoms);
  std::vector<int> imported_xcrd_ovrf, imported_ycrd_ovrf, imported_zcrd_ovrf;
  if (tcalc_is_double) {
    imported_xcrd_ovrf.resize(maximum_valence_work_unit_atoms);
    imported_ycrd_ovrf.resize(maximum_valence_work_unit_atoms);
    imported_zcrd_ovrf.resize(maximum_valence_work_unit_atoms);
  }

  // Loop over all valence work units, which contain the constraint instructions
  for (int vwu_idx = 0; vwu_idx < poly_vk.nvwu; vwu_idx++) {
    
    // Import atoms as the work unit would on the GPU.
    const int2 import_limits = poly_vk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                                     static_cast<int>(VwuAbstractMap::IMPORT)];
    for (int i = import_limits.x; i < import_limits.y; i++) {
      const int local_idx  = i - import_limits.x;
      const size_t global_idx = poly_vk.vwu_imports[i];
      imported_atom_ids[local_idx] = global_idx;
      imported_xcrd[local_idx] = poly_psw->xcrd[global_idx];
      imported_ycrd[local_idx] = poly_psw->ycrd[global_idx];
      imported_zcrd[local_idx] = poly_psw->zcrd[global_idx];
      if (tcalc_is_double) {
        imported_xcrd_ovrf[local_idx] = poly_psw->xcrd_ovrf[global_idx];
        imported_ycrd_ovrf[local_idx] = poly_psw->ycrd_ovrf[global_idx];
        imported_zcrd_ovrf[local_idx] = poly_psw->zcrd_ovrf[global_idx];
      }
    }

    // Loop over instruction batches, honoring warp coalescence
    const int2 insr_limits = poly_vk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +  
                                                   static_cast<int>(VwuAbstractMap::CGROUP)];
    std::vector<uint2> warp_insr(warp_size_int);
    std::vector<T> dx_ref(warp_size_int), dy_ref(warp_size_int), dz_ref(warp_size_int);
    std::vector<llint> ca_vx(warp_size_int), ca_vy(warp_size_int), ca_vz(warp_size_int);
    std::vector<llint> ph_vx(warp_size_int), ph_vy(warp_size_int), ph_vz(warp_size_int);
    std::vector<int> ca_vx_ovrf(warp_size_int), ph_vx_ovrf(warp_size_int);
    std::vector<int> ca_vy_ovrf(warp_size_int), ph_vy_ovrf(warp_size_int);
    std::vector<int> ca_vz_ovrf(warp_size_int), ph_vz_ovrf(warp_size_int);
    std::vector<T> ca_invmass(warp_size_int), ph_invmass(warp_size_int);
    std::vector<T> l2_target(warp_size_int), combined_invmass(warp_size_int);
    for (int i = insr_limits.x; i < insr_limits.y; i += warp_size_int) {
      const int batch_size = std::min(warp_size_int, insr_limits.y - i);
      
      // Compute the reference displacement that each thread in the GPU warp would be working with.
      // This can be obtained from the pre-cached coordinates, which will not have changed
      // positions (the forces have been computed, but only velocities have been updated by half a
      // step).  Unlike the single-system routine, which is designed to run on a single CPU thread,
      // the reference displacement (which will not change as the velocities are being constrained)
      // can be computed outside of the iteration loop, as each thread is assigned to a particular
      // bond and there is thus no nested loop within the iterations.  Similarly, the peripheral
      // atom masses and constraint length targets can be taken into "registers" outside of the
      // iteration loop.  Here, the particle velocities are placed in similar arrays as if they
      // were in registers, but it may be necessary to instead swap the forces into the L1 thread
      // block space and use those __shared__ memory arrays to store the velocities as warps
      // iterate to convergence.
      for (int lane = 0; lane < batch_size; lane++) {
        const uint2 tinsr = poly_auk.cnst_insr[i + lane];
        const int central_atom = (tinsr.x & 0x3ff);
        const int peripheral_atom = ((tinsr.x >> 10) & 0x3ff);

        // Check that the constraint is valid
        if (central_atom == 0 && peripheral_atom == 0) {
          continue;
        }

        // Load the displacement, atom velocities, and constraint parameters
        const size_t ca_gbl_idx = imported_atom_ids[central_atom];
        const size_t ph_gbl_idx = imported_atom_ids[peripheral_atom];
        if (tcalc_is_double) {
          const int95_t idx_ref = hostInt95Sum(imported_xcrd[peripheral_atom],
                                               imported_xcrd_ovrf[peripheral_atom],
                                               -imported_xcrd[central_atom],
                                               -imported_xcrd_ovrf[central_atom]);
          const int95_t idy_ref = hostInt95Sum(imported_ycrd[peripheral_atom],
                                               imported_ycrd_ovrf[peripheral_atom],
                                               -imported_ycrd[central_atom],
                                               -imported_ycrd_ovrf[central_atom]);
          const int95_t idz_ref = hostInt95Sum(imported_zcrd[peripheral_atom],
                                               imported_zcrd_ovrf[peripheral_atom],
                                               -imported_zcrd[central_atom],
                                               -imported_zcrd_ovrf[central_atom]);
          dx_ref[lane] = hostInt95ToDouble(idx_ref) * poly_psw->inv_gpos_scale_f;
          dy_ref[lane] = hostInt95ToDouble(idy_ref) * poly_psw->inv_gpos_scale_f;
          dz_ref[lane] = hostInt95ToDouble(idz_ref) * poly_psw->inv_gpos_scale_f;
          ca_vx[lane] = poly_psw->vxalt[ca_gbl_idx];
          ca_vy[lane] = poly_psw->vyalt[ca_gbl_idx];
          ca_vz[lane] = poly_psw->vzalt[ca_gbl_idx];
          ca_vx_ovrf[lane] = poly_psw->vxalt_ovrf[ca_gbl_idx];
          ca_vy_ovrf[lane] = poly_psw->vyalt_ovrf[ca_gbl_idx];
          ca_vz_ovrf[lane] = poly_psw->vzalt_ovrf[ca_gbl_idx];
          ph_vx[lane] = poly_psw->vxalt[ph_gbl_idx];
          ph_vy[lane] = poly_psw->vyalt[ph_gbl_idx];
          ph_vz[lane] = poly_psw->vzalt[ph_gbl_idx];
          ph_vx_ovrf[lane] = poly_psw->vxalt_ovrf[ph_gbl_idx];
          ph_vy_ovrf[lane] = poly_psw->vyalt_ovrf[ph_gbl_idx];
          ph_vz_ovrf[lane] = poly_psw->vzalt_ovrf[ph_gbl_idx];
        }
        else {
          const llint idx_ref = imported_xcrd[peripheral_atom] - imported_xcrd[central_atom];
          const llint idy_ref = imported_ycrd[peripheral_atom] - imported_ycrd[central_atom];
          const llint idz_ref = imported_zcrd[peripheral_atom] - imported_zcrd[central_atom];
          dx_ref[lane] = static_cast<T>(idx_ref) * poly_psw->inv_gpos_scale_f;
          dy_ref[lane] = static_cast<T>(idy_ref) * poly_psw->inv_gpos_scale_f;
          dz_ref[lane] = static_cast<T>(idz_ref) * poly_psw->inv_gpos_scale_f;
          ca_vx[lane] = poly_psw->vxalt[ca_gbl_idx];
          ca_vy[lane] = poly_psw->vyalt[ca_gbl_idx];
          ca_vz[lane] = poly_psw->vzalt[ca_gbl_idx];
          ph_vx[lane] = poly_psw->vxalt[ph_gbl_idx];
          ph_vy[lane] = poly_psw->vyalt[ph_gbl_idx];
          ph_vz[lane] = poly_psw->vzalt[ph_gbl_idx];
        }
        ca_invmass[lane]       = poly_auk.inv_masses[ca_gbl_idx];
        ph_invmass[lane]       = poly_auk.inv_masses[ph_gbl_idx];
        l2_target[lane]        = poly_auk.cnst_grp_params[tinsr.y].x;
        combined_invmass[lane] = poly_auk.cnst_grp_params[tinsr.y].y;
      }
      bool converged = false;
      int iter = 0;
      while (converged == false && iter < max_iter) {
        converged = true;
        for (int lane = 0; lane < batch_size; lane++) {
          const uint2 tinsr = poly_auk.cnst_insr[i + lane];
          const int central_atom = (tinsr.x & 0x3ff);
          const int peripheral_atom = ((tinsr.x >> 10) & 0x3ff);

          // Compute the current velocity differential
          T dvx, dvy, dvz;
          if (tcalc_is_double) {
            dvx = hostInt95ToDouble(hostInt95Sum( ph_vx[lane],  ph_vx_ovrf[lane],
                                                 -ca_vx[lane], -ca_vx_ovrf[lane])) *
                  poly_psw->inv_vel_scale_f;
            dvy = hostInt95ToDouble(hostInt95Sum( ph_vy[lane],  ph_vy_ovrf[lane],
                                                 -ca_vy[lane], -ca_vy_ovrf[lane])) *
                  poly_psw->inv_vel_scale_f;
            dvz = hostInt95ToDouble(hostInt95Sum( ph_vz[lane], ph_vz_ovrf[lane],
                                                 -ca_vz[lane], -ca_vz_ovrf[lane])) *
                  poly_psw->inv_vel_scale_f;
          }
          else {
            dvx = static_cast<T>(ph_vx[lane] - ca_vx[lane]) * poly_psw->inv_vel_scale_f;
            dvy = static_cast<T>(ph_vy[lane] - ca_vy[lane]) * poly_psw->inv_vel_scale_f;
            dvz = static_cast<T>(ph_vz[lane] - ca_vz[lane]) * poly_psw->inv_vel_scale_f;
          }
          
          // Compute the dot product of the reference displacement and the velocity differential.
          // These vectors should be nearly orthogonal, at least to within a tolerance when scaled
          // by the combined inverse masses and the squared distance target.
          const T dot_rv = (dx_ref[lane] * dvx) + (dy_ref[lane] * dvy) + (dz_ref[lane] * dvz);
          const T term = (tcalc_is_double) ?
                         -1.2  * dot_rv * poly_psw->vel_scale_f /
                         (combined_invmass[lane] * l2_target[lane]) :
                         -1.2f * dot_rv * poly_psw->vel_scale_f /
                         (combined_invmass[lane] * l2_target[lane]);
          if (fabs(term) > rtoldt) {
            converged = false;
            if (tcalc_is_double) {
              const int95_t phvx_updt = hostDoubleToInt95(dx_ref[lane] * term * ph_invmass[lane]);
              const int95_t phvy_updt = hostDoubleToInt95(dy_ref[lane] * term * ph_invmass[lane]);
              const int95_t phvz_updt = hostDoubleToInt95(dz_ref[lane] * term * ph_invmass[lane]);
              const int95_t nph_vx = hostSplitFPSum(phvx_updt, ph_vx[lane], ph_vx_ovrf[lane]);
              const int95_t nph_vy = hostSplitFPSum(phvy_updt, ph_vy[lane], ph_vy_ovrf[lane]);
              const int95_t nph_vz = hostSplitFPSum(phvz_updt, ph_vz[lane], ph_vz_ovrf[lane]);
              ph_vx[lane] = nph_vx.x;
              ph_vy[lane] = nph_vy.x;
              ph_vz[lane] = nph_vz.x;
              ph_vx_ovrf[lane] = nph_vx.y;
              ph_vy_ovrf[lane] = nph_vy.y;
              ph_vz_ovrf[lane] = nph_vz.y;
              const int95_t cavx_updt = hostDoubleToInt95(dx_ref[lane] * term * ca_invmass[lane]);
              const int95_t cavy_updt = hostDoubleToInt95(dy_ref[lane] * term * ca_invmass[lane]);
              const int95_t cavz_updt = hostDoubleToInt95(dz_ref[lane] * term * ca_invmass[lane]);

              // On the CPU, only the leader thread for each group will update its central atom
              // velocity.  Other threads will replace their information about the central atom's
              // velocity with the result from the iteration.  This is in preparation for a second
              // pass, mimicking the __shfl() operation, in which the leader thread will accumulate
              // other threads' contributions to the central atom velocity.
              const int leader_lane = ((tinsr.x >> 20) & 0xff);
              if (lane == leader_lane) {
                const int95_t nca_vx = hostInt95Sum(-cavx_updt.x, -cavx_updt.y, ca_vx[lane],
                                                    ca_vx_ovrf[lane]);
                const int95_t nca_vy = hostInt95Sum(-cavy_updt.x, -cavy_updt.y, ca_vy[lane],
                                                    ca_vy_ovrf[lane]);
                const int95_t nca_vz = hostInt95Sum(-cavz_updt.x, -cavz_updt.y, ca_vz[lane],
                                                    ca_vz_ovrf[lane]);
                ca_vx[lane] = nca_vx.x;
                ca_vy[lane] = nca_vy.x;
                ca_vz[lane] = nca_vz.x;
                ca_vx_ovrf[lane] = nca_vx.y;
                ca_vy_ovrf[lane] = nca_vy.y;
                ca_vz_ovrf[lane] = nca_vz.y;
              }
              else {
                ca_vx[lane] = -cavx_updt.x;
                ca_vy[lane] = -cavy_updt.x;
                ca_vz[lane] = -cavz_updt.x;
                ca_vx_ovrf[lane] = -cavx_updt.y;
                ca_vy_ovrf[lane] = -cavy_updt.y;
                ca_vz_ovrf[lane] = -cavz_updt.y;
              }
            }
            else {
              ph_vx[lane] += llround(dx_ref[lane] * term * ph_invmass[lane]);
              ph_vy[lane] += llround(dy_ref[lane] * term * ph_invmass[lane]);
              ph_vz[lane] += llround(dz_ref[lane] * term * ph_invmass[lane]);
              const llint cavx_update = llround(dx_ref[lane] * term * ca_invmass[lane]);
              const llint cavy_update = llround(dy_ref[lane] * term * ca_invmass[lane]);
              const llint cavz_update = llround(dz_ref[lane] * term * ca_invmass[lane]);
              const int leader_lane = ((tinsr.x >> 20) & 0xff);
              if (lane == leader_lane) {
                ca_vx[lane] -= cavx_update;
                ca_vy[lane] -= cavy_update;
                ca_vz[lane] -= cavz_update;
              }
              else {
                ca_vx[lane] = -cavx_update;
                ca_vy[lane] = -cavy_update;
                ca_vz[lane] = -cavz_update;
              }
            }
          }
          else {
            const int leader_lane = ((tinsr.x >> 20) & 0xff);
            if (lane != leader_lane) {
              if (tcalc_is_double) {
                ca_vx[lane] = 0LL;
                ca_vy[lane] = 0LL;
                ca_vz[lane] = 0LL;
                ca_vx_ovrf[lane] = 0;
                ca_vy_ovrf[lane] = 0;
                ca_vz_ovrf[lane] = 0;
              }
              else {
                ca_vx[lane] = 0LL;
                ca_vy[lane] = 0LL;
                ca_vz[lane] = 0LL;
              }
            }
          }
        }

        // With the batch done, the CPU code has emulated the work of a synchronized warp.  Begin
        // reducing the moves on central atoms.
        for (int lane = 0; lane < batch_size; lane++) {
          const uint2 tinsr = poly_auk.cnst_insr[i + lane];
          const int central_atom = (tinsr.x & 0x3ff);
          const int leader_lane = ((tinsr.x >> 20) & 0xff);
          if (lane == leader_lane) {
            const int n_partners = (tinsr.x >> 28);
            if (tcalc_is_double) {
              for (int j = 1; j < n_partners; j++) {
                const int95_t nca_vx = hostInt95Sum(ca_vx[lane], ca_vx_ovrf[lane], ca_vx[lane + j],
                                                    ca_vx_ovrf[lane + j]);
                const int95_t nca_vy = hostInt95Sum(ca_vy[lane], ca_vy_ovrf[lane], ca_vy[lane + j],
                                                    ca_vy_ovrf[lane + j]);
                const int95_t nca_vz = hostInt95Sum(ca_vz[lane], ca_vz_ovrf[lane], ca_vz[lane + j],
                                                    ca_vz_ovrf[lane + j]);
                ca_vx[lane] = nca_vx.x;
                ca_vy[lane] = nca_vy.x;
                ca_vz[lane] = nca_vz.x;
                ca_vx_ovrf[lane] = nca_vx.y;
                ca_vy_ovrf[lane] = nca_vy.y;
                ca_vz_ovrf[lane] = nca_vz.y;
              }
            }
            else {
              for (int j = 1; j < n_partners; j++) {
                ca_vx[lane] += ca_vx[lane + j];
                ca_vy[lane] += ca_vy[lane + j];
                ca_vz[lane] += ca_vz[lane + j];
              }
            }
          }
        }

        // Broadcast the moves on central atoms to all threads assigned to each constraint group.
        for (int lane = 0; lane < batch_size; lane++) {
          const uint2 tinsr = poly_auk.cnst_insr[i + lane];
          const int leader_lane = ((tinsr.x >> 20) & 0xff);
          if (lane != leader_lane) {
            if (tcalc_is_double) {
              ca_vx[lane] = ca_vx[leader_lane];
              ca_vy[lane] = ca_vy[leader_lane];
              ca_vz[lane] = ca_vz[leader_lane];
              ca_vx_ovrf[lane] = ca_vx_ovrf[leader_lane];
              ca_vy_ovrf[lane] = ca_vy_ovrf[leader_lane];
              ca_vz_ovrf[lane] = ca_vz_ovrf[leader_lane];
            }
            else {
              ca_vx[lane] = ca_vx[leader_lane];
              ca_vy[lane] = ca_vy[leader_lane];
              ca_vz[lane] = ca_vz[leader_lane];
            }
          }
        }

        // Increment the iteration count
        iter++;
      }

      // Check for non-converged groups
      if (converged == false && iter == max_iter) {
        rtErr("A constraint group did not converge when iterated in a coordinate synthesis.",
              "rattleVelocities");
      }

      // Write the constrained velocities back to main memory.
      for (int lane = 0; lane < batch_size; lane++) {
        const uint2 tinsr = poly_auk.cnst_insr[i + lane];
        const int central_atom = (tinsr.x & 0x3ff);
        const int peripheral_atom = ((tinsr.x >> 10) & 0x3ff);

        // Check that the constraint is valid
        if (central_atom == 0 && peripheral_atom == 0) {
          continue;
        }

        // All threads managing any instruction will have a peripheral atom to cover.  However, it
        // must be ascertained that the work unit has authority to move that atom.
        const int nbits = sizeof(uint) * 8;
        const int2 mask_limits =
          poly_vk.vwu_abstracts[(vwu_idx * vwu_abstract_length) +
                                static_cast<int>(VwuAbstractMap::MANIPULATE)];
        const int ph_updt_elem = peripheral_atom / nbits;
        const int ph_updt_bit = peripheral_atom - (ph_updt_elem * nbits);
        if ((poly_auk.vwu_manip[mask_limits.x + ph_updt_elem].y >> ph_updt_bit) & 0x1) {
          const size_t ph_gbl_idx = imported_atom_ids[peripheral_atom];
          poly_psw->vxalt[ph_gbl_idx] = ph_vx[lane];
          poly_psw->vyalt[ph_gbl_idx] = ph_vy[lane];
          poly_psw->vzalt[ph_gbl_idx] = ph_vz[lane];            
          if (tcalc_is_double) {
            poly_psw->vxalt_ovrf[ph_gbl_idx] = ph_vx_ovrf[lane];
            poly_psw->vyalt_ovrf[ph_gbl_idx] = ph_vy_ovrf[lane];
            poly_psw->vzalt_ovrf[ph_gbl_idx] = ph_vz_ovrf[lane];            
          }
        }
        
        // Leader lanes will write the velocities of their central atom.
        const int leader_lane = ((tinsr.x >> 20) & 0xff);
        if (lane == leader_lane) {

          // Check that the work unit has authority to move the central atom.
          const int ca_updt_elem = central_atom / nbits;
          const int ca_updt_bit  = central_atom - (ca_updt_elem * nbits);
          if ((poly_auk.vwu_manip[mask_limits.x + ca_updt_elem].y >> ca_updt_bit) & 0x1) {
            const size_t ca_gbl_idx = imported_atom_ids[central_atom];
            poly_psw->vxalt[ca_gbl_idx] = ca_vx[lane];
            poly_psw->vyalt[ca_gbl_idx] = ca_vy[lane];
            poly_psw->vzalt[ca_gbl_idx] = ca_vz[lane];
            if (tcalc_is_double) {
              poly_psw->vxalt_ovrf[ca_gbl_idx] = ca_vx_ovrf[lane];
              poly_psw->vyalt_ovrf[ca_gbl_idx] = ca_vy_ovrf[lane];
              poly_psw->vzalt_ovrf[ca_gbl_idx] = ca_vz_ovrf[lane];
            }
          }
        }
      }
    }
  }
}

} // namespace structure
} // namespace stormm
