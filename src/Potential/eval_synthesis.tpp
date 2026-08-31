// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

using namespace generalized_born_defaults;
using topology::ImplicitSolventModel;

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc2>
void evalSyNonbondedTileGroups(const SyNonbondedKit<Tcalc, Tcalc2> synbk,
                               const SeMaskSynthesisReader syse,
                               PsSynthesisWriter *psyw, ScoreCard *ecard,
                               const NonbondedTask task, const EvaluateForce eval_elec_force,
                               const EvaluateForce eval_vdw_force,
                               const Tcalc clash_minimum_distance, const Tcalc clash_ratio,
                               ISWorkspaceKit<Tcalc> *iswk) {

  // Guard against Generalized Born computations when no GB model is in effect, or computations
  // without a valid ImplicitSolventWorkspace when a GB model is in effect.
  const bool gb_engaged = (synbk.igb != ImplicitSolventModel::NONE);
  switch (synbk.igb) {
  case ImplicitSolventModel::NONE:
    switch (task) {
    case NonbondedTask::GB_RADII:
    case NonbondedTask::GB_RADII_DERIVATIVES:
      return;
    case NonbondedTask::PARTICLE_PARTICLE:
      break;
    }
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    if (iswk == nullptr) {
      rtErr("An ImplicitSolventWorkspace abstract must be provided for Generalized Born "
            "calculations involving task " + getEnumerationName(task) + ".",
            "evalSyNonbondedTileGroups");
    }
    break;
  }

  // Critical indexing and offsets for each work unit
  std::vector<int> sh_nbwu_abstract(tile_groups_wu_abstract_length);
  std::vector<int> sh_system_indices(small_block_max_imports);
  std::vector<int> sh_n_lj_types(small_block_max_imports);
  std::vector<int> sh_ljabc_offsets(small_block_max_imports);

  // Pre-computed information for rapidly manipulating the particles of any one tile
  std::vector<Tcalc> sh_tile_xcog(small_block_max_imports);
  std::vector<Tcalc> sh_tile_ycog(small_block_max_imports);
  std::vector<Tcalc> sh_tile_zcog(small_block_max_imports);
  std::vector<Tcalc> sh_tile_tpts(small_block_max_imports);

  // L1-cached coordinates of particles.  These will be accessed repeatedly as the list of tile
  // instructions gets processed.
  std::vector<llint> lc_xcrd(small_block_max_atoms);
  std::vector<llint> lc_ycrd(small_block_max_atoms);
  std::vector<llint> lc_zcrd(small_block_max_atoms);
  std::vector<int> lc_xcrd_overflow(small_block_max_atoms);
  std::vector<int> lc_ycrd_overflow(small_block_max_atoms);
  std::vector<int> lc_zcrd_overflow(small_block_max_atoms);
  std::vector<Tcalc> lc_charge(small_block_max_atoms);
  std::vector<int>   lc_lj_idx(small_block_max_atoms);
  std::vector<Tcalc> lc_gb_radius(small_block_max_atoms);
  std::vector<Tcalc> lc_gb_screen(small_block_max_atoms);
  std::vector<int>   lc_neck_idx(small_block_max_atoms);

  // Local force accumulators, stored in __shared__ on the GPU to do atomic operations to L1.
  std::vector<llint> sh_xfrc(small_block_max_atoms);
  std::vector<llint> sh_yfrc(small_block_max_atoms);
  std::vector<llint> sh_zfrc(small_block_max_atoms);
  std::vector<int> sh_xfrc_overflow(small_block_max_atoms);
  std::vector<int> sh_yfrc_overflow(small_block_max_atoms);
  std::vector<int> sh_zfrc_overflow(small_block_max_atoms);
  std::vector<llint> sh_psi(small_block_max_atoms);
  std::vector<int> sh_psi_overflow(small_block_max_atoms);
  std::vector<Tcalc> sh_sum_deijda(small_block_max_atoms);
  std::vector<int> sh_sum_deijda_overflow(small_block_max_atoms);

  // Arrays for mocking the register content of various threads
  std::vector<Tcalc> reg_xcrd(tile_length * 2), reg_xfrc(tile_length * 2);
  std::vector<Tcalc> reg_ycrd(tile_length * 2), reg_yfrc(tile_length * 2);
  std::vector<Tcalc> reg_zcrd(tile_length * 2), reg_zfrc(tile_length * 2);
  std::vector<uint>  reg_excl(tile_length * 2);
  std::vector<int>   reg_lj_idx(tile_length * 2), reg_neck_idx(tile_length * 2);
  std::vector<Tcalc> reg_charge(tile_length * 2), reg_sum_deijda(tile_length * 2);
  std::vector<Tcalc> reg_radius(tile_length * 2), reg_screen(tile_length * 2);

  // Constants must be cast in the proper calculation precision
  const Tcalc v_one  = 1.0;
  const Tcalc v_two = 2.0;
  const Tcalc v_thre = 3.0;
  const Tcalc v_four = 4.0;
  const Tcalc v_half = 0.5;
  const Tcalc v_qrtr = 0.25;
  const Tcalc v_pthr = 0.3;
  const Tcalc v_opei = 1.8;
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc sqrt_coul = (tcalc_is_double) ? sqrt(synbk.coulomb) : sqrtf(synbk.coulomb);
  const Tcalc nrg_scale_factor = ecard->getEnergyScalingFactor<Tcalc>();
  const bool do_elec_force = (eval_elec_force == EvaluateForce::YES);
  const bool do_vdw_force  = (eval_vdw_force == EvaluateForce::YES);
  const bool do_either_force = (do_elec_force || do_vdw_force);
  const bool do_neck = (synbk.igb == ImplicitSolventModel::NECK_GB ||
                        synbk.igb == ImplicitSolventModel::NECK_GB_II);
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
  
  // Initialize the appropriate energy terms in all systems
  switch (task) {
  case NonbondedTask::PARTICLE_PARTICLE:
    for (int i = 0; i < psyw->system_count; i++) {
      ecard->initialize(StateVariable::ELECTROSTATIC, i);
      ecard->initialize(StateVariable::VDW, i);
      if (gb_engaged) {
        ecard->initialize(StateVariable::GENERALIZED_BORN, i);
      }
    }
    break;
  case NonbondedTask::GB_RADII:
  case NonbondedTask::GB_RADII_DERIVATIVES:
    break;
  }
  
  // Loop over all non-bonded work units within the topology synthesis.  Each holds within it a
  // list of atom imports and a series of tiles.  Within each tile, all atoms pertain to the same
  // system, but within each work unit, depending on the boundary conditions and the non-bonded
  // list type, different tiles may correspond to different systems.
  for (int nbwu_idx = 0; nbwu_idx < synbk.nnbwu; nbwu_idx++) {

    // Import the abstract.
    for (int pos = 0; pos < tile_groups_wu_abstract_length; pos++) {
      sh_nbwu_abstract[pos] = synbk.nbwu_abstracts[(tile_groups_wu_abstract_length * nbwu_idx) +
                                                   pos];
    }
    const int ntile_sides = sh_nbwu_abstract[0];
    for (int pos = 0; pos < ntile_sides; pos++) {
      const int system_idx   = sh_nbwu_abstract[pos + 28];
      sh_system_indices[pos] = system_idx;
      sh_n_lj_types[pos]     = synbk.n_lj_types[system_idx];
      sh_ljabc_offsets[pos]  = synbk.ljabc_offsets[system_idx];
    }
    const int tile_insr_start = sh_nbwu_abstract[small_block_max_imports + 6];
    const int tile_insr_end   = sh_nbwu_abstract[small_block_max_imports + 7];
    const uint init_mask = sh_nbwu_abstract[(2 * small_block_max_imports) + 8];
    
    // Import atoms into the appropriate arrays.  Prepare to compute the center of geometry for
    // all imported atoms.
    for (int pos = 0; pos < ntile_sides; pos++) {
      const int atom_start_idx = sh_nbwu_abstract[pos + 1];
      const int system_idx     = sh_system_indices[pos];
      const int tside_count = hostGetTileSideAtomCount(sh_nbwu_abstract, pos);
      
      // Pre-compute the centers of geometry for each batch of tile_length atoms, storing the
      // results (totals plus weights) in floating-point format.  When it comes time to do actual
      // tiles, combine the results for the abscissa and ordinate atoms, divide by the combined
      // weight, and use that number to shift the tile atoms to the optimal, precision-preserving
      // locations in the fixed-precision format prior to converting to floating point numbers.
      Tcalc x_cog = 0.0;
      Tcalc y_cog = 0.0;
      Tcalc z_cog = 0.0;
      Tcalc t_pts = 0.0;
      llint gb_self_acc = 0LL;
      for (int i = 0; i < tside_count; i++) {
        const size_t localpos = (tile_length * pos) + i;
        const size_t synthpos = atom_start_idx + i;
        lc_xcrd[localpos]    = psyw->xcrd[synthpos];
        lc_ycrd[localpos]    = psyw->ycrd[synthpos];
        lc_zcrd[localpos]    = psyw->zcrd[synthpos];
        if (psyw->gpos_bits > globalpos_scale_nonoverflow_bits) {
          lc_xcrd_overflow[localpos] = psyw->xcrd_ovrf[synthpos];
          lc_ycrd_overflow[localpos] = psyw->ycrd_ovrf[synthpos];
          lc_zcrd_overflow[localpos] = psyw->zcrd_ovrf[synthpos];
        }
        sh_xfrc[localpos]    = 0LL;
        sh_yfrc[localpos]    = 0LL;
        sh_zfrc[localpos]    = 0LL;
        if (psyw->frc_bits > force_scale_nonoverflow_bits) {
          sh_xfrc_overflow[localpos] = 0;
          sh_yfrc_overflow[localpos] = 0;
          sh_zfrc_overflow[localpos] = 0;
        }
        
        // Pre-scale all charges by the square root of Coulomb's constant so that it will carry
        // through in all subsequent electrostatic calculations.  On the GPU, the latency involved
        // in this step will likely hide all of the cost of the extra multiplication.
        lc_charge[localpos]  = synbk.charge[synthpos] * sqrt_coul;

        // The Lennard-Jones indices recorded here are specific to each system.  For each tile,
        // it is still critical to know the relevant system's total number of Lennard-Jones types
        // as well as its table offset.  On the GPU, these pieces of information will be imported
        // into __shared__ memory for convenient access.
        lc_lj_idx[localpos]  = synbk.lj_idx[synthpos];

        // For Generalized Born computations, baseline atomic radii and screening factors must be
        // loaded.
        if (task == NonbondedTask::GB_RADII || task == NonbondedTask::GB_RADII_DERIVATIVES) {
          lc_gb_radius[localpos]  = synbk.pb_radii[synthpos] - synbk.gb_offset;
          lc_gb_screen[localpos]  = synbk.gb_screen[synthpos];
          sh_psi[localpos]        = 0LL;
          if (tcalc_is_double) {
            sh_psi_overflow[localpos] = 0;
          }
          if (do_neck) {
            lc_neck_idx[localpos] = synbk.neck_gb_idx[synthpos];
          }
        }
        else if (task == NonbondedTask::PARTICLE_PARTICLE && gb_engaged) {
          const Tcalc pb_radius = synbk.pb_radii[synthpos];
          const Tcalc gb_radius = pb_radius - synbk.gb_offset;
          const Tcalc inv_gb_radius = v_one / gb_radius;
          Tcalc psival;
          if (tcalc_is_double) {
            psival = hostInt95ToDouble(iswk->psi[synthpos], iswk->psi_ovrf[synthpos]) *
                     iswk->inv_fp_scale;
          }
          else {
            psival = static_cast<Tcalc>(iswk->psi[synthpos]) * iswk->inv_fp_scale;
          }
          Tcalc egbi = 0.0;
          switch (synbk.igb) {
          case ImplicitSolventModel::HCT_GB:
            egbi = v_one / (inv_gb_radius + psival);
            if (egbi < 0.0) {
              egbi = 30.0;
            }
            break;
          case ImplicitSolventModel::OBC_GB:
          case ImplicitSolventModel::OBC_GB_II:
          case ImplicitSolventModel::NECK_GB:
          case ImplicitSolventModel::NECK_GB_II:
            {
              const Tcalc fipsi = psival * (-gb_radius);
              egbi = v_one / (inv_gb_radius -
                              ((tcalc_is_double) ?
                               tanh((synbk.gb_alpha[synthpos] -
                                     (synbk.gb_beta[synthpos] -
                                      (synbk.gb_gamma[synthpos] * fipsi)) * fipsi) * fipsi) :
                               tanhf((synbk.gb_alpha[synthpos] -
                                      (synbk.gb_beta[synthpos] -
                                       (synbk.gb_gamma[synthpos] * fipsi)) * fipsi) * fipsi)) /
                               pb_radius);
            }
            break;
          case ImplicitSolventModel::NONE:
            break;
          }
          lc_gb_radius[localpos] = egbi;
          sh_sum_deijda[localpos] = 0.0;
          if (tcalc_is_double) {
            sh_sum_deijda_overflow[localpos] = 0;
          }
          if ((init_mask >> pos) & 0x1) {
            const Tcalc atomq = lc_charge[localpos];
            Tcalc expmkf;
            if (tcalc_is_double) {
              expmkf = exp(-default_gb_kscale * synbk.kappa * egbi) / synbk.dielectric;
            }
            else {
              expmkf = expf(-default_gb_kscale_f * synbk.kappa * egbi) / synbk.dielectric;
            }
            const Tcalc dielfac = v_one - expmkf;
            const Tcalc atmq2h = v_half * atomq * atomq;
            const Tcalc atmqd2h = atmq2h * dielfac;
            const Tcalc gb_self_nrg = -atmqd2h / egbi;
            gb_self_acc += llround(gb_self_nrg * nrg_scale_factor);
            if (do_either_force) {
              const Tcalc sdi = atmqd2h - ((tcalc_is_double) ? default_gb_kscale :
                                             default_gb_kscale_f) * synbk.kappa * atmq2h *
                                expmkf * egbi;
              sh_sum_deijda[localpos] = sdi;
            }
          }
        }

        // Center of geometry computation--coordinates are not scaled to real units at this stage
        if (psyw->gpos_bits > globalpos_scale_nonoverflow_bits) {
          x_cog += hostInt95ToDouble(lc_xcrd[localpos], lc_xcrd_overflow[localpos]);
          y_cog += hostInt95ToDouble(lc_ycrd[localpos], lc_ycrd_overflow[localpos]);
          z_cog += hostInt95ToDouble(lc_zcrd[localpos], lc_zcrd_overflow[localpos]);
        }
        else {
          x_cog += static_cast<Tcalc>(lc_xcrd[localpos]);
          y_cog += static_cast<Tcalc>(lc_ycrd[localpos]);
          z_cog += static_cast<Tcalc>(lc_zcrd[localpos]);
        }
        t_pts += v_one;
      }

      // Accumulate GB self energy contributions if needed
      if (task == NonbondedTask::PARTICLE_PARTICLE && gb_engaged) {
        ecard->add(StateVariable::GENERALIZED_BORN, gb_self_acc, system_idx);
      }
      
      // As on the GPU, mock atoms get properties that will deliver neutered effects when
      // interacting with one another or with real particles.  Serial summation of the center of
      // geometry obviates the need for zeroing contributions from these mock atoms.
      for (int i = tside_count; i < tile_length; i++) {
        const size_t localpos = (tile_length * pos) + i;
        const double xdum = static_cast<double>(1024 * i) * psyw->gpos_scale;
        const double ydum = static_cast<double>(1152 * i) * psyw->gpos_scale;
        const double zdum = static_cast<double>(1280 * i) * psyw->gpos_scale;
        if (psyw->gpos_bits > globalpos_scale_nonoverflow_bits) {
          const int95_t ixdum = hostDoubleToInt95(xdum);
          const int95_t iydum = hostDoubleToInt95(ydum);
          const int95_t izdum = hostDoubleToInt95(zdum);
          lc_xcrd[localpos] = ixdum.x;
          lc_ycrd[localpos] = iydum.x;
          lc_zcrd[localpos] = izdum.x;
          lc_xcrd_overflow[localpos] = ixdum.y;
          lc_ycrd_overflow[localpos] = iydum.y;
          lc_zcrd_overflow[localpos] = izdum.y;
        }
        else {
          lc_xcrd[localpos] = static_cast<llint>(xdum);
          lc_ycrd[localpos] = static_cast<llint>(ydum);
          lc_zcrd[localpos] = static_cast<llint>(zdum);
        }
        lc_gb_radius[localpos]  = v_one;
        sh_sum_deijda[localpos] = 0.0;
      }

      // Log the center of coordinates.
      sh_tile_xcog[pos] = x_cog;
      sh_tile_ycog[pos] = y_cog;
      sh_tile_zcog[pos] = z_cog;
      sh_tile_tpts[pos] = t_pts;
    }

    // Load the sum_deijda accumulators and prepare baseline GB radii for derivative force
    // calculations, following the protocol of gbderivative_tilegroups.cui.
    if (task == NonbondedTask::GB_RADII_DERIVATIVES) {
      for (int pos = 0; pos < ntile_sides; pos++) {
        const int atom_start_idx = sh_nbwu_abstract[pos + 1];
        const int tside_count = hostGetTileSideAtomCount(sh_nbwu_abstract, pos);
        for (int i = 0; i < tside_count; i++) {
          const size_t localpos = (tile_length * pos) + i;
          const size_t synthpos = atom_start_idx + i;
          const Tcalc gb_radius = lc_gb_radius[localpos];
          Tcalc sdi_current;
          if (tcalc_is_double) {
            sdi_current = hostInt95ToDouble(iswk->sum_deijda[synthpos],
                                            iswk->sum_deijda_ovrf[synthpos]) *
                          iswk->inv_fp_scale;
          }
          else {
            sdi_current = static_cast<Tcalc>(iswk->sum_deijda[synthpos]) * iswk->inv_fp_scale;
          }
          switch (synbk.igb) {
          case ImplicitSolventModel::NONE:
          case ImplicitSolventModel::HCT_GB:
            sh_sum_deijda[localpos] = sdi_current;
            break;
          case ImplicitSolventModel::OBC_GB:
          case ImplicitSolventModel::OBC_GB_II:
          case ImplicitSolventModel::NECK_GB:
          case ImplicitSolventModel::NECK_GB_II:
            {
              Tcalc psival;
              if (tcalc_is_double) {
                psival = hostInt95ToDouble(iswk->psi[synthpos], iswk->psi_ovrf[synthpos]) *
                         iswk->inv_fp_scale;
              }
              else {
                psival = static_cast<Tcalc>(iswk->psi[synthpos]) * iswk->inv_fp_scale;
              }
              const Tcalc fipsi = psival * (-gb_radius);
              const Tcalc thi = (tcalc_is_double) ?
                                tanh((synbk.gb_alpha[synthpos] -
                                      (synbk.gb_beta[synthpos] -
                                       (synbk.gb_gamma[synthpos] * fipsi)) * fipsi) * fipsi) :
                                tanhf((synbk.gb_alpha[synthpos] -
                                       (synbk.gb_beta[synthpos] -
                                        (synbk.gb_gamma[synthpos] * fipsi)) * fipsi) * fipsi);
              const Tcalc sdi_multiplier = (synbk.gb_alpha[synthpos] -
                                            ((v_two * synbk.gb_beta[synthpos]) -
                                             (v_thre * synbk.gb_gamma[synthpos] * fipsi)) *
                                            fipsi) * (v_one - thi * thi) * gb_radius /
                                           synbk.pb_radii[synthpos];
              sh_sum_deijda[localpos] = sdi_current * sdi_multiplier;
            }
            break;
          }
        }
      }
    }
    
    // Loop over tile instructions
    for (int pos = tile_insr_start; pos < tile_insr_end; pos++) {
      uint2 tinsr = synbk.nbwu_insr[pos];
      for (int i = 0; i < 2 * tile_length; i++) {
        reg_excl[i] = syse.mask_data[tinsr.y + i];
      }
      const int local_absc_start = (tinsr.x & 0xffff);
      const int local_ordi_start = ((tinsr.x >> 16) & 0xffff);
      const int absc_import_idx = local_absc_start >> tile_length_bits;
      const int ordi_import_idx = local_ordi_start >> tile_length_bits;
      const bool on_diagonal = (absc_import_idx == ordi_import_idx);
      
      // The system index is needed in order to know where to accumulate the resulting energy
      const int system_idx = sh_system_indices[absc_import_idx];
      
      // On the GPU, the atomic coordinates stored in signed long long integers will be converted
      // to floating point numbers at the beginning of the tile calculation, but to preserve as
      // much information as possible that conversion will also involve centering those atoms on
      // the tile's center of geometry.
      const Tcalc inv_tile_pts = v_one /
                                 (sh_tile_tpts[absc_import_idx] + sh_tile_tpts[ordi_import_idx]);
      const Tcalc tx_cog = (sh_tile_xcog[absc_import_idx] + sh_tile_xcog[ordi_import_idx]) *
                           inv_tile_pts;
      const Tcalc ty_cog = (sh_tile_ycog[absc_import_idx] + sh_tile_ycog[ordi_import_idx]) *
                           inv_tile_pts;
      const Tcalc tz_cog = (sh_tile_zcog[absc_import_idx] + sh_tile_zcog[ordi_import_idx]) *
                           inv_tile_pts;
      if (psyw->gpos_bits > globalpos_scale_nonoverflow_bits) {
        const int95_t x_center = hostDoubleToInt95(tx_cog);
        const int95_t y_center = hostDoubleToInt95(ty_cog);
        const int95_t z_center = hostDoubleToInt95(tz_cog);

        // Based on knowledge of the proper centering, re-import the coordinates after making the
        // shift in precision-preserving integer arithmetic.
        for (int i = 0; i < tile_length; i++) {
          const size_t ilabsc = i + local_absc_start;
          const size_t ilordi = i + local_ordi_start;
          const size_t iplust = i + tile_length;
          const int95_t ix_absc = hostInt95Subtract(lc_xcrd[ilabsc], lc_xcrd_overflow[ilabsc],
                                                    x_center.x, x_center.y);
          const int95_t iy_absc = hostInt95Subtract(lc_ycrd[ilabsc], lc_ycrd_overflow[ilabsc],
                                                    y_center.x, y_center.y);
          const int95_t iz_absc = hostInt95Subtract(lc_zcrd[ilabsc], lc_zcrd_overflow[ilabsc],
                                                    z_center.x, z_center.y);
          reg_xcrd[i]   = hostInt95ToDouble(ix_absc) * psyw->inv_gpos_scale;
          reg_ycrd[i]   = hostInt95ToDouble(iy_absc) * psyw->inv_gpos_scale;
          reg_zcrd[i]   = hostInt95ToDouble(iz_absc) * psyw->inv_gpos_scale;
          const int95_t ix_ordi = hostInt95Subtract(lc_xcrd[ilordi], lc_xcrd_overflow[ilordi],
                                                    x_center.x, x_center.y);
          const int95_t iy_ordi = hostInt95Subtract(lc_ycrd[ilordi], lc_ycrd_overflow[ilordi],
                                                    y_center.x, y_center.y);
          const int95_t iz_ordi = hostInt95Subtract(lc_zcrd[ilordi], lc_zcrd_overflow[ilordi],
                                                    z_center.x, z_center.y);
          reg_xcrd[iplust] = hostInt95ToDouble(ix_ordi) * psyw->inv_gpos_scale;
          reg_ycrd[iplust] = hostInt95ToDouble(iy_ordi) * psyw->inv_gpos_scale;
          reg_zcrd[iplust] = hostInt95ToDouble(iz_ordi) * psyw->inv_gpos_scale;
        }
      }
      else {
        const llint x_center = static_cast<llint>(tx_cog);
        const llint y_center = static_cast<llint>(ty_cog);
        const llint z_center = static_cast<llint>(tz_cog);

        // Based on knowledge of the proper centering, re-import the coordinates after making the
        // shift in precision-preserving integer arithemtic.
        for (int i = 0; i < tile_length; i++) {
          const size_t ilabsc = i + local_absc_start;
          const size_t ilordi = i + local_ordi_start;
          const size_t iplust = i + tile_length;
          reg_xcrd[i]   = static_cast<Tcalc>(lc_xcrd[ilabsc] - x_center) * psyw->inv_gpos_scale;
          reg_ycrd[i]   = static_cast<Tcalc>(lc_ycrd[ilabsc] - y_center) * psyw->inv_gpos_scale;
          reg_zcrd[i]   = static_cast<Tcalc>(lc_zcrd[ilabsc] - z_center) * psyw->inv_gpos_scale;
          reg_xcrd[iplust] = static_cast<Tcalc>(lc_xcrd[ilordi] - x_center) * psyw->inv_gpos_scale;
          reg_ycrd[iplust] = static_cast<Tcalc>(lc_ycrd[ilordi] - y_center) * psyw->inv_gpos_scale;
          reg_zcrd[iplust] = static_cast<Tcalc>(lc_zcrd[ilordi] - z_center) * psyw->inv_gpos_scale;
        }
      }

      // Load parameter indices and real-valued accumulators into mock register arrays, regardless
      // of whether overflow buffers are active in the coordinate representation.
      for (int i = 0; i < tile_length; i++) {
        const size_t ilabsc = i + local_absc_start;
        const size_t ilordi = i + local_ordi_start;
        const size_t iplust = i + tile_length;
        reg_lj_idx[i] = lc_lj_idx[ilabsc];
        reg_charge[i] = lc_charge[ilabsc];
        reg_lj_idx[iplust] = lc_lj_idx[ilordi];
        reg_charge[iplust] = lc_charge[ilordi];
        switch (synbk.igb) {
        case ImplicitSolventModel::NONE:
          reg_radius[i] = 1.0;
          reg_screen[i] = 1.0;
          reg_radius[iplust] = 1.0;
          reg_screen[iplust] = 1.0;
          break;
        case ImplicitSolventModel::HCT_GB:
        case ImplicitSolventModel::OBC_GB:
        case ImplicitSolventModel::OBC_GB_II:
        case ImplicitSolventModel::NECK_GB:
        case ImplicitSolventModel::NECK_GB_II:
          reg_radius[i]       = lc_gb_radius[ilabsc];
          reg_screen[i]       = lc_gb_screen[ilabsc];
          reg_radius[iplust]  = lc_gb_radius[ilordi];
          reg_screen[iplust]  = lc_gb_screen[ilordi];

          // The following cases need only be initialized if some GB-related process is in effect.
          if (do_neck) {
            reg_neck_idx[i]      = lc_neck_idx[ilabsc];
            reg_neck_idx[iplust] = lc_neck_idx[ilordi];
          }
          if (task == NonbondedTask::GB_RADII_DERIVATIVES) {
            reg_sum_deijda[i]      = sh_sum_deijda[ilabsc];
            reg_sum_deijda[iplust] = sh_sum_deijda[ilordi];
          }
          break;
        }
      }

      // Initialize forces, if appropriate
      switch (task) {
      case NonbondedTask::PARTICLE_PARTICLE:
      case NonbondedTask::GB_RADII_DERIVATIVES:
        for (int i = 0; i < 2 * tile_length; i++) {
          reg_xfrc[i] = 0.0;
          reg_yfrc[i] = 0.0;
          reg_zfrc[i] = 0.0;
        }
        break;
      case NonbondedTask::GB_RADII:
        break;
      }

      // Scan over the entire tile--if the tile runs past the end of the system's atoms, there
      // will be a solid mask of excluded interactions.
      switch (task) {
      case NonbondedTask::PARTICLE_PARTICLE:
        {
          const int nljt = sh_n_lj_types[absc_import_idx];
          const int lj_offset = sh_ljabc_offsets[absc_import_idx];
          Tcalc elec_nrg = 0.0;
          Tcalc vdw_nrg  = 0.0;
          Tcalc gb_nrg   = 0.0;
          for (int i = 0; i < tile_length; i++) {
            const uint i_mask = reg_excl[i];
            const Tcalc xi = reg_xcrd[i];
            const Tcalc yi = reg_ycrd[i];
            const Tcalc zi = reg_zcrd[i];
            const Tcalc qi = reg_charge[i];
            const int ilj_idx = (nljt * reg_lj_idx[i]) + lj_offset;
            const Tcalc atomi_radius = reg_radius[i];

            // For particle-particle electrostatic and van-der Waals interactions, all
            // interactions of on-diagonal tiles are marked excluded for j >= i.  Generalized Born
            // interactions are never excluded by the mask, but on-diagonal tiles receive a factor
            // of one half and self-interactions are omitted.
            for (int j = tile_length; j < 2 * tile_length; j++) {
              const bool pair_excluded = ((i_mask >> (j - tile_length)) & 0x1);
              const Tcalc dx       = reg_xcrd[j] - xi;
              const Tcalc dy       = reg_ycrd[j] - yi;
              const Tcalc dz       = reg_zcrd[j] - zi;
              const Tcalc r2       = (dx * dx) + (dy * dy) + (dz * dz);
              const Tcalc dr       = (tcalc_is_double) ? sqrt(r2) : sqrtf(r2);
              const Tcalc invr     = v_one / dr;
              const Tcalc invr2    = invr * invr;
              const Tcalc invr4    = invr2 * invr2;
              const Tcalc qqij     = reg_charge[j] * qi;
              const int   ij_ljidx = reg_lj_idx[j] + ilj_idx;

              // Generalized Born pair terms (always evaluated when a GB model is active)
              const Tcalc atomij_radius = atomi_radius * reg_radius[j];
              Tcalc efac    = 0.0;
              Tcalc fgbi    = 0.0;
              Tcalc fgbk    = 0.0;
              Tcalc expmkf  = 0.0;
              Tcalc dielfac = 0.0;
              if (gb_engaged) {
                if (tcalc_is_double) {
                  efac = exp(-r2 / (v_four * atomij_radius));
                  fgbi = v_one / sqrt(r2 + (atomij_radius * efac));
                  fgbk = -synbk.kappa * default_gb_kscale / fgbi;
                  expmkf = exp(fgbk) / synbk.dielectric;
                }
                else {
                  efac = expf(-r2 / (v_four * atomij_radius));
                  fgbi = v_one / sqrtf(r2 + (atomij_radius * efac));
                  fgbk = -synbk.kappa * default_gb_kscale_f / fgbi;
                  expmkf = expf(fgbk) / synbk.dielectric;
                }
                dielfac = v_one - expmkf;
                Tcalc gb_pair_nrg = -qqij * dielfac * fgbi;
                if (on_diagonal) {
                  if (i != j - tile_length) {
                    gb_nrg += v_half * gb_pair_nrg;
                  }
                }
                else {
                  gb_nrg += gb_pair_nrg;
                }
              }

              // Log the vacuum electrostatic and van-der Waals energies when not excluded.
              Tcalc2 ljab;
              if (! pair_excluded) {
                ljab = synbk.ljab_coeff[ij_ljidx];
                elec_nrg += qqij * invr;
                vdw_nrg  += ((ljab.x * invr4 * invr2) - ljab.y) * invr4 * invr2;
              }

              // Compute the forces and contribute them to accumulators.
              if (do_either_force) {
                Tcalc fmag = 0.0;
                if (gb_engaged) {
                  const Tcalc temp4 = fgbi * fgbi * fgbi;
                  const Tcalc temp6 = qqij * temp4 * (dielfac + (fgbk * expmkf));
                  fmag = temp6 * (v_one - (v_qrtr * efac));
                  Tcalc temp5 = v_half * efac * temp6 * (atomij_radius + (v_qrtr * r2));
                  if (on_diagonal) {
                    fmag *= v_half;
                    temp5 *= v_half;
                    if (i == j - tile_length) {
                      fmag = 0.0;
                      temp5 = 0.0;
                    }
                  }
                  const size_t ilabsc = i + local_absc_start;
                  const size_t ilordi = (j - tile_length) + local_ordi_start;
                  sh_sum_deijda[ilabsc] += temp5 * atomi_radius;
                  sh_sum_deijda[ilordi] += temp5 * reg_radius[j];
                }
                if (! pair_excluded) {
                  if (do_elec_force) {
                    fmag -= qqij * invr * invr2;
                  }
                  if (do_vdw_force) {
                    if (tcalc_is_double) {
                      fmag += ((6.0 * ljab.y) - (12.0 * ljab.x * invr2 * invr4)) * invr4 * invr4;
                    }
                    else {
                      fmag += ((6.0f * ljab.y) - (12.0f * ljab.x * invr2 * invr4)) * invr4 * invr4;
                    }
                  }
                }
                const Tcalc fmag_dx = fmag * dx;
                const Tcalc fmag_dy = fmag * dy;
                const Tcalc fmag_dz = fmag * dz;
                reg_xfrc[i] += fmag_dx;
                reg_yfrc[i] += fmag_dy;
                reg_zfrc[i] += fmag_dz;
                reg_xfrc[j] -= fmag_dx;
                reg_yfrc[j] -= fmag_dy;
                reg_zfrc[j] -= fmag_dz;
              }
            }
          }

          // There is no need to test whether each work unit is responsible for accumulating the
          // energy computed in a tile.  However, each tile will need to contribute its result to
          // the energy accumulator for a particular system, due to the fact that work units can
          // contain tiles from different systems.
          const llint elec_acc = llround(elec_nrg * nrg_scale_factor);
          ecard->add(StateVariable::ELECTROSTATIC, elec_acc, system_idx);
          const llint vdw_acc  = llround(vdw_nrg * nrg_scale_factor);
          ecard->add(StateVariable::VDW, vdw_acc, system_idx);
          if (gb_engaged) {
            const llint gb_acc = llround(gb_nrg * nrg_scale_factor);
            ecard->add(StateVariable::GENERALIZED_BORN, gb_acc, system_idx);
          }
        }
        break;
      case NonbondedTask::GB_RADII:
        for (int i = 0; i < tile_length; i++) {
          const uint i_mask = reg_excl[i];
          const Tcalc xi = reg_xcrd[i];
          const Tcalc yi = reg_ycrd[i];
          const Tcalc zi = reg_zcrd[i];
          const Tcalc atomi_radius = reg_radius[i];
          const Tcalc atomi_inv_radius = v_one / atomi_radius;
          for (int j = tile_length; j < 2 * tile_length; j++) {
            const Tcalc dx = reg_xcrd[j] - xi;
            const Tcalc dy = reg_ycrd[j] - yi;
            const Tcalc dz = reg_zcrd[j] - zi;
            const Tcalc r2 = (dx * dx) + (dy * dy) + (dz * dz);
            const Tcalc r  = (tcalc_is_double) ? sqrt(r2) : sqrtf(r2);
            const Tcalc invr = v_one / r;
            const Tcalc atomj_radius = reg_radius[j];
            const Tcalc atomj_inv_radius = v_one / atomj_radius;

            // First computation: atom I -> atom J
            const Tcalc sj = reg_screen[j] * atomj_radius;
            const Tcalc sj2 = sj * sj;
            Tcalc t_psi = 0.0;
            if (r > v_four * sj) {
              const Tcalc invr2 = invr * invr;
              const Tcalc tmpsd = sj2 * invr2;
              const Tcalc dumbo = gta + tmpsd * (gtb + tmpsd * (gtc + tmpsd * (gtd + tmpsd *
                                                                               gtdd)));
              t_psi -= sj * tmpsd * invr2 * dumbo;
            }
            else if (r > atomi_radius + sj) {
              if (tcalc_is_double) {
                t_psi -= v_half * ((sj / (r2 - sj2)) +
                                   (v_half * invr * log((r - sj) / (r + sj))));
              }
              else {
                t_psi -= v_half * ((sj / (r2 - sj2)) +
                                   (v_half * invr * logf((r - sj) / (r + sj))));
              }
            }
            else if (r > fabs(atomi_radius - sj)) {
              const Tcalc theta = v_half * atomi_inv_radius * invr *
                                  (r2 + (atomi_radius * atomi_radius) - sj2);
              const Tcalc uij   = v_one / (r + sj);
              if (tcalc_is_double) {
                t_psi -= v_qrtr * ((atomi_inv_radius * (v_two - theta)) -
                                   uij + (invr * log(atomi_radius * uij)));
              }
              else {
                t_psi -= v_qrtr * ((atomi_inv_radius * (v_two - theta)) -
                                   uij + (invr * logf(atomi_radius * uij)));
              }
            }
            else if (atomi_radius < sj) {
              if (tcalc_is_double) {
                t_psi -= v_half * ((sj / (r2 - sj2)) + (v_two * atomi_inv_radius) +
                                   (v_half * invr * log((sj - r) / (sj + r))));
              }
              else {
                t_psi -= v_half * ((sj / (r2 - sj2)) + (v_two * atomi_inv_radius) +
                                   (v_half * invr * logf((sj - r) / (sj + r))));
              }
            }

            // Second computation: atom J -> atom I
            const Tcalc si = reg_screen[i] * atomi_radius;
            const Tcalc si2 = si * si;
            Tcalc o_psi = 0.0;
            if (r > v_four * si) {
              const Tcalc invr2  = invr * invr;
              const Tcalc tmpsd  = si2 * invr2;
              const Tcalc dumbo  = gta + tmpsd * (gtb + tmpsd * (gtc + tmpsd * (gtd + tmpsd *
                                                                                gtdd)));
              o_psi -= si * tmpsd * invr2 * dumbo;
            }
            else if (r > atomj_radius + si) {
              if (tcalc_is_double) {
                o_psi -= v_half * ((si / (r2 - si2)) +
                                   (v_half * invr * log((r - si) / (r + si))));
              }
              else {
                o_psi -= v_half * ((si / (r2 - si2)) +
                                   (v_half * invr * logf((r - si) / (r + si))));
              }
            }
            else if (r > fabs(atomj_radius - si)) {
              const Tcalc theta = v_half * atomj_inv_radius * invr *
                                  (r2 + (atomj_radius * atomj_radius) - si2);
              const Tcalc uij   = v_one / (r + si);
              if (tcalc_is_double) {
                o_psi -= v_qrtr * (atomj_inv_radius * (v_two - theta) - uij +
                                   invr * log(atomj_radius * uij));
              }
              else {
                o_psi -= v_qrtr * (atomj_inv_radius * (v_two - theta) - uij +
                                   invr * logf(atomj_radius * uij));
              }
            }
            else if (atomj_radius < si) {
              if (tcalc_is_double) {
                o_psi -= v_half * ((si / (r2 - si2)) + (v_two * atomj_inv_radius) +
                                   (v_half * invr * log((si - r) / (si + r))));
              }
              else {
                o_psi -= v_half * ((si / (r2 - si2)) + (v_two * atomj_inv_radius) +
                                   (v_half * invr * logf((si - r) / (si + r))));
              }
            }

            // Neck GB contributions
            if (do_neck &&
                r < atomi_radius + atomj_radius + (v_two * synbk.gb_offset) + synbk.gb_neckcut &&
                reg_neck_idx[i] >= 0 && reg_neck_idx[j] >= 0) {
              const int ij_table_idx = (synbk.neck_table_size * reg_neck_idx[j]) +
                                       reg_neck_idx[i];
              Tcalc mdist  = r - synbk.neck_limits[ij_table_idx].x;
              Tcalc mdist2 = mdist * mdist;
              Tcalc mdist6 = mdist2 * mdist2 * mdist2;
              t_psi -= synbk.gb_neckscale * synbk.neck_limits[ij_table_idx].y /
                       (v_one + mdist2 + (v_pthr * mdist6));
              const int ji_table_idx = (synbk.neck_table_size * reg_neck_idx[i]) +
                                       reg_neck_idx[j];
              mdist  = r - synbk.neck_limits[ji_table_idx].x;
              mdist2 = mdist * mdist;
              mdist6 = mdist2 * mdist2 * mdist2;
              o_psi -= synbk.gb_neckscale * synbk.neck_limits[ji_table_idx].y /
                       (v_one + mdist2 + (v_pthr * mdist6));
            }

            if (on_diagonal) {
              t_psi *= v_half;
              o_psi *= v_half;
              if (i == j - tile_length) {
                t_psi = 0.0;
                o_psi = 0.0;
              }
            }

            const size_t ilabsc = i + local_absc_start;
            const size_t ilordi = j - tile_length + local_ordi_start;
            if (tcalc_is_double) {
              const int95_t npsi_i = hostInt95Sum(sh_psi[ilabsc], sh_psi_overflow[ilabsc],
                                                  t_psi * iswk->fp_scale);
              sh_psi[ilabsc] = npsi_i.x;
              sh_psi_overflow[ilabsc] = npsi_i.y;
              const int95_t npsi_j = hostInt95Sum(sh_psi[ilordi], sh_psi_overflow[ilordi],
                                                  o_psi * iswk->fp_scale);
              sh_psi[ilordi] = npsi_j.x;
              sh_psi_overflow[ilordi] = npsi_j.y;
            }
            else {
              sh_psi[ilabsc] += llround(t_psi * iswk->fp_scale);
              sh_psi[ilordi] += llround(o_psi * iswk->fp_scale);
            }
          }
        }
        break;
      case NonbondedTask::GB_RADII_DERIVATIVES:
        {
          for (int i = 0; i < tile_length; i++) {
            const uint i_mask = reg_excl[i];
            const Tcalc xi = reg_xcrd[i];
            const Tcalc yi = reg_ycrd[i];
            const Tcalc zi = reg_zcrd[i];
            const Tcalc atomi_radius = reg_radius[i];
            const Tcalc atomi_inv_radius = v_one / atomi_radius;
            for (int j = tile_length; j < 2 * tile_length; j++) {
              const Tcalc dx = reg_xcrd[j] - xi;
              const Tcalc dy = reg_ycrd[j] - yi;
              const Tcalc dz = reg_zcrd[j] - zi;
              const Tcalc atomj_radius = reg_radius[j];
              const Tcalc atomj_inv_radius = v_one / atomj_radius;
              const Tcalc r2 = (dx * dx) + (dy * dy) + (dz * dz);
              const Tcalc invr = (tcalc_is_double) ? v_one / sqrt(r2) :
                                                     v_one / sqrtf(r2);
              const Tcalc invr2 = invr * invr;
              const Tcalc r = r2 * invr;

              // First computation: atom I -> atom J
              const Tcalc sj = reg_screen[j] * atomj_radius;
              const Tcalc sj2 = sj * sj;
              Tcalc datmpi, datmpj;
              if (r > v_four * sj) {
                const Tcalc tmpsd  = sj2 * invr2;
                const Tcalc dumbo  = gte + tmpsd * (gtf + tmpsd * (gtg + tmpsd * (gth + tmpsd *
                                                                                  gthh)));
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
                datmpi = 0.0;
              }

              // Second computation: atom J -> atom I
              const Tcalc si = reg_screen[i] * atomi_radius;
              const Tcalc si2 = si * si;
              if (r > v_four * si) {
                const Tcalc tmpsd  = si2 * invr2;
                const Tcalc dumbo  = gte + tmpsd * (gtf + tmpsd * (gtg + tmpsd * (gth + tmpsd *
                                                                                  gthh)));
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
                datmpj = 0.0;
              }

              // Neck GB contributions
              if (do_neck &&
                  r < atomi_radius + atomj_radius + (v_two * synbk.gb_offset) + synbk.gb_neckcut &&
                  reg_neck_idx[i] >= 0 && reg_neck_idx[j] >= 0) {
                const int ij_table_idx = (synbk.neck_table_size * reg_neck_idx[j]) +
                                         reg_neck_idx[i];
                Tcalc mdist = r - synbk.neck_limits[ij_table_idx].x;
                Tcalc mdist2 = mdist * mdist;
                Tcalc mdist6 = mdist2 * mdist2 * mdist2;
                Tcalc temp1 = v_one + mdist2 + (v_pthr * mdist6);
                temp1 = temp1 * temp1 * r;
                datmpi += (((v_two * mdist) + (v_opei * mdist2 * mdist2 * mdist)) *
                           synbk.neck_limits[ij_table_idx].y * synbk.gb_neckscale) / temp1;
                const int ji_table_idx = (synbk.neck_table_size * reg_neck_idx[i]) +
                                         reg_neck_idx[j];
                mdist = r - synbk.neck_limits[ji_table_idx].x;
                mdist2 = mdist * mdist;
                mdist6 = mdist2 * mdist2 * mdist2;
                temp1 = v_one + mdist2 + (v_pthr * mdist6);
                temp1 = temp1 * temp1 * r;
                datmpj += (((v_two * mdist) + (v_opei * mdist2 * mdist2 * mdist)) *
                           synbk.neck_limits[ji_table_idx].y * synbk.gb_neckscale) / temp1;
              }
              Tcalc fmag = (datmpi * reg_sum_deijda[i]) + (datmpj * reg_sum_deijda[j]);
              if (on_diagonal) {
                fmag *= v_half;
                if (i == j - tile_length) {
                  fmag = 0.0;
                }
              }
              reg_xfrc[i] -= fmag * dx;
              reg_yfrc[i] -= fmag * dy;
              reg_zfrc[i] -= fmag * dz;
              reg_xfrc[j] += fmag * dx;
              reg_yfrc[j] += fmag * dy;
              reg_zfrc[j] += fmag * dz;
            }
          }
        }
        break;
      }
      
      // Store forces computed in the tile in local accumulators.  This mirrors what will happen on
      // the GPU as forces stored in floating point numbers in registers will be converted back to
      // fixed-precision and contributed to __shared__ arrays at the end of the cycle.
      switch (task) {
      case NonbondedTask::PARTICLE_PARTICLE:
      case NonbondedTask::GB_RADII_DERIVATIVES:
        if (do_either_force || task == NonbondedTask::GB_RADII_DERIVATIVES) {
          for (int i = 0; i < tile_length; i++) {
            const size_t ilabsc = i + local_absc_start;
            const size_t ilordi = i + local_ordi_start;
            const size_t iplust = i + tile_length;
            if (psyw->frc_bits > force_scale_nonoverflow_bits) {
              const int95_t nfx = hostInt95Sum(sh_xfrc[ilabsc], sh_xfrc_overflow[ilabsc],
                                               reg_xfrc[i] * psyw->frc_scale);
              const int95_t nfy = hostInt95Sum(sh_yfrc[ilabsc], sh_yfrc_overflow[ilabsc],
                                               reg_yfrc[i] * psyw->frc_scale);
              const int95_t nfz = hostInt95Sum(sh_zfrc[ilabsc], sh_zfrc_overflow[ilabsc],
                                               reg_zfrc[i] * psyw->frc_scale);
              sh_xfrc[ilabsc] = nfx.x;
              sh_yfrc[ilabsc] = nfy.x;
              sh_zfrc[ilabsc] = nfz.x;
              sh_xfrc_overflow[ilabsc] = nfx.y;
              sh_yfrc_overflow[ilabsc] = nfy.y;
              sh_zfrc_overflow[ilabsc] = nfz.y;
              const int95_t pfx = hostInt95Sum(sh_xfrc[ilordi], sh_xfrc_overflow[ilordi],
                                               reg_xfrc[iplust] * psyw->frc_scale);
              const int95_t pfy = hostInt95Sum(sh_yfrc[ilordi], sh_yfrc_overflow[ilordi],
                                               reg_yfrc[iplust] * psyw->frc_scale);
              const int95_t pfz = hostInt95Sum(sh_zfrc[ilordi], sh_zfrc_overflow[ilordi],
                                               reg_zfrc[iplust] * psyw->frc_scale);
              sh_xfrc[ilordi] = pfx.x;
              sh_yfrc[ilordi] = pfy.x;
              sh_zfrc[ilordi] = pfz.x;
              sh_xfrc_overflow[ilordi] = pfx.y;
              sh_yfrc_overflow[ilordi] = pfy.y;
              sh_zfrc_overflow[ilordi] = pfz.y;
            }
            else {
              sh_xfrc[ilabsc] += llround(reg_xfrc[i] * psyw->frc_scale);
              sh_yfrc[ilabsc] += llround(reg_yfrc[i] * psyw->frc_scale);
              sh_zfrc[ilabsc] += llround(reg_zfrc[i] * psyw->frc_scale);
              sh_xfrc[ilordi] += llround(reg_xfrc[iplust] * psyw->frc_scale);
              sh_yfrc[ilordi] += llround(reg_yfrc[iplust] * psyw->frc_scale);
              sh_zfrc[ilordi] += llround(reg_zfrc[iplust] * psyw->frc_scale);
            }
          }
        }
        break;
      case NonbondedTask::GB_RADII:
        break;
      }
    }

    // Contribute local force accumulators back to global, or local psi accumulators in the case
    // of GB radii computations.
    switch (task) {
    case NonbondedTask::PARTICLE_PARTICLE:
    case NonbondedTask::GB_RADII_DERIVATIVES:
      for (int pos = 0; pos < ntile_sides; pos++) {
        const int atom_start_idx = sh_nbwu_abstract[pos + 1];
        const int tside_count = hostGetTileSideAtomCount(sh_nbwu_abstract, pos);
        for (int i = 0; i < tside_count; i++) {
          const size_t localpos = (tile_length * pos) + i;
          const size_t synthpos = atom_start_idx + i;
          if (psyw->frc_bits > force_scale_nonoverflow_bits) {
            const int95_t nfx = hostInt95Sum(psyw->xfrc[synthpos], psyw->xfrc_ovrf[synthpos],
                                             sh_xfrc[localpos], sh_xfrc_overflow[localpos]);
            const int95_t nfy = hostInt95Sum(psyw->yfrc[synthpos], psyw->yfrc_ovrf[synthpos],
                                             sh_yfrc[localpos], sh_yfrc_overflow[localpos]);
            const int95_t nfz = hostInt95Sum(psyw->zfrc[synthpos], psyw->zfrc_ovrf[synthpos],
                                             sh_zfrc[localpos], sh_zfrc_overflow[localpos]);
            psyw->xfrc[synthpos] = nfx.x;
            psyw->yfrc[synthpos] = nfy.x;
            psyw->zfrc[synthpos] = nfz.x;
            psyw->xfrc_ovrf[synthpos] = nfx.y;
            psyw->yfrc_ovrf[synthpos] = nfy.y;
            psyw->zfrc_ovrf[synthpos] = nfz.y;
          }
          else {
            psyw->xfrc[synthpos] += sh_xfrc[localpos];
            psyw->yfrc[synthpos] += sh_yfrc[localpos];
            psyw->zfrc[synthpos] += sh_zfrc[localpos];
          }
        }
      }
      break;
    case NonbondedTask::GB_RADII:
      for (int pos = 0; pos < ntile_sides; pos++) {
        const int atom_start_idx = sh_nbwu_abstract[pos + 1];
        const int tside_count = hostGetTileSideAtomCount(sh_nbwu_abstract, pos);
        for (int i = 0; i < tside_count; i++) {
          const size_t localpos = (tile_length * pos) + i;
          const size_t synthpos = atom_start_idx + i;
          if (tcalc_is_double) {
            const int95_t npsi = hostInt95Sum(iswk->psi[synthpos], iswk->psi_ovrf[synthpos],
                                              sh_psi[localpos], sh_psi_overflow[localpos]);
            iswk->psi[synthpos] = npsi.x;
            iswk->psi_ovrf[synthpos] = npsi.y;
          }
          else {
            iswk->psi[synthpos] += sh_psi[localpos];
          }
        }
      }
      break;
    }

    // Commit Generalized Born radii derivative accumulators from pair interactions
    if (task == NonbondedTask::PARTICLE_PARTICLE && gb_engaged && do_either_force) {
      for (int pos = 0; pos < ntile_sides; pos++) {
        const int atom_start_idx = sh_nbwu_abstract[pos + 1];
        const int tside_count = hostGetTileSideAtomCount(sh_nbwu_abstract, pos);
        for (int i = 0; i < tside_count; i++) {
          const size_t localpos = (tile_length * pos) + i;
          const size_t synthpos = atom_start_idx + i;
          if (tcalc_is_double) {
            const int95_t nsdi = hostInt95Sum(iswk->sum_deijda[synthpos],
                                              iswk->sum_deijda_ovrf[synthpos],
                                              llround(sh_sum_deijda[localpos] * iswk->fp_scale));
            iswk->sum_deijda[synthpos] = nsdi.x;
            iswk->sum_deijda_ovrf[synthpos] = nsdi.y;
          }
          else {
            iswk->sum_deijda[synthpos] += llround(sh_sum_deijda[localpos] * iswk->fp_scale);
          }
        }
      }
    }
  }
}

} // namespace energy
} // namespace stormm
