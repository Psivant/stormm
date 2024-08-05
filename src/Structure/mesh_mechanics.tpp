// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc, typename Tcoord, typename Tobs>
double interpolate(const BackgroundMeshReader<Tdata> &bgmr, const MeshFFKit<Tcalc> &mnbk,
                   const Tcalc* prop_a, const Tcalc* prop_b, const int* prop_idx,
                   const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                   const int* xcrd_ovrf, const int* ycrd_ovrf, const int* zcrd_ovrf,
                   const int natom, Tobs* vnrg, Tcoord* xfrc, Tcoord* yfrc, Tcoord* zfrc,
                   int* xfrc_ovrf, int* yfrc_ovrf, int* zfrc_ovrf, const OffMeshProtocol policy,
                   const int system_index, const int crd_scaling_bits,
                   const int frc_scaling_bits, const int nrg_scaling_bits) {

  // Constants for the calculation precision
  const Tcalc value_zero        = 0.0;
  const Tcalc value_one         = 1.0;
  const Tcalc value_sixteen     = 4.0;
  const Tcoord coord_value_zero = 0.0;
  const size_t sixty_four   = 64;
  const Tcalc dfactors[4] = { 0.0, 1.0, 2.0, 3.0 };

  // Check whether the grid works in a format that would restrict the fixed-precision math to just
  // the first 64 bits of the representation.
  const size_t data_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (data_ct == double_type_index);  
  const bool coordinate_overflow_active = (xcrd_ovrf != nullptr && ycrd_ovrf != nullptr &&
                                           zcrd_ovrf != nullptr);
  const bool compute_forces = (xfrc != nullptr && yfrc != nullptr && zfrc != nullptr);
  const bool force_overflow_active = (xfrc_ovrf != nullptr && yfrc_ovrf != nullptr &&
                                      zfrc_ovrf != nullptr);
  const bool coord_is_integral = (crd_scaling_bits > 0);
  const Tcalc coordinate_scale = pow(2.0, crd_scaling_bits);
  const Tcalc coordinate_inv_scale = static_cast<Tcalc>(1.0) / coordinate_scale;
  const bool force_is_integral = (frc_scaling_bits > 0);
  const bool obs_is_integral = (isFloatingPointScalarType<Tobs>() == false);
  const Tcalc force_scale = pow(2.0, frc_scaling_bits);
  const Tcalc nrg_scale_factor = pow(2.0, nrg_scaling_bits);
  const Tcalc mesh_orig_x = hostInt95ToDouble(bgmr.dims.orig_x) * bgmr.dims.inv_scale;
  const Tcalc mesh_orig_y = hostInt95ToDouble(bgmr.dims.orig_y) * bgmr.dims.inv_scale;
  const Tcalc mesh_orig_z = hostInt95ToDouble(bgmr.dims.orig_z) * bgmr.dims.inv_scale;
  const std::vector<Tcalc> mesh_umat = { static_cast<Tcalc>(bgmr.dims.umat[0]),
                                         static_cast<Tcalc>(bgmr.dims.umat[1]),
                                         static_cast<Tcalc>(bgmr.dims.umat[2]),
                                         static_cast<Tcalc>(bgmr.dims.umat[3]),
                                         static_cast<Tcalc>(bgmr.dims.umat[4]),
                                         static_cast<Tcalc>(bgmr.dims.umat[5]),
                                         static_cast<Tcalc>(bgmr.dims.umat[6]),
                                         static_cast<Tcalc>(bgmr.dims.umat[7]),
                                         static_cast<Tcalc>(bgmr.dims.umat[8]) };
  const Tcalc real_na = bgmr.dims.na;
  const Tcalc real_nb = bgmr.dims.nb;
  const Tcalc real_nc = bgmr.dims.nc;
  
  // Loop over all particles.  Initialize the energy contribution.
  llint nrg_acc = 0LL;
  for (int pos = 0; pos < natom; pos++) {
    switch (bgmr.field) {
    case NonbondedPotential::CLASH:
      {
        const int pidx = (pos >> 5);
        const int prem = pos - (pidx << 5);
        if ((static_cast<uint>(prop_idx[pidx]) >> prem) == 0) {
          continue;
        }
      }
      break;
    case NonbondedPotential::ELECTROSTATIC:
      break;
    case NonbondedPotential::VAN_DER_WAALS:
      {
        // Determine whether the self-interaction parameters of the atom in the molecule being
        // moved around the mesh match those of the probe.  The subtle differences between an
        // AtomGraph and AtomGraphSynthesis, in this respect, present a challenge for abstracting
        // this function to both objects.  For Lennard-Jones parameters, prop_a contains the
        // Lennard-Jones sigma parameters and prop_b contains the Lennard-Jones B coefficients,
        // from which the epsilon can be readily extracted.
        const size_t pidx = prop_idx[pos];
        const Tcalc p_sigma = sqrt(cbrt(prop_a[pidx] / prop_b[pidx]));
        const Tcalc p_sigma3 = p_sigma * p_sigma * p_sigma;
        const Tcalc p_epsilon = 0.25 * prop_b[pidx] / (p_sigma3 * p_sigma3);
        if (fabs(bgmr.probe_radius - p_sigma) > 1.0e-5 ||
            fabs(bgmr.well_depth - p_epsilon) > 1.0e-5) {
          continue;
        }
      }
      break;
    }
    
    // Determine the grid bin using a transformation into mesh space followed by a local check and,
    // if necessary, a correction.  The precision of the mesh coefficients will determine the
    // precision of the fixed-precision coordinate representation, following the standard set forth
    // in validateScalingBits() above.
    int95_t ixcrd, iycrd, izcrd;
    if (coord_is_integral) {

      // The coordinates may not have the same precision as the positions expected by the mesh.
      // Harmonize the representations.
      if (crd_scaling_bits == bgmr.dims.scale_bits) {
        ixcrd.x = xcrd[pos];
        iycrd.x = ycrd[pos];
        izcrd.x = zcrd[pos];
        if (coordinate_overflow_active) {
          ixcrd.y = xcrd_ovrf[pos];
          iycrd.y = ycrd_ovrf[pos];
          izcrd.y = zcrd_ovrf[pos];       
        }
      }
      else {
        int95_t orep_xcrd, orep_ycrd, orep_zcrd;
        if (coordinate_overflow_active) {
          orep_xcrd = { static_cast<llint>(xcrd[pos]), xcrd_ovrf[pos] };
          orep_ycrd = { static_cast<llint>(ycrd[pos]), ycrd_ovrf[pos] };
          orep_zcrd = { static_cast<llint>(zcrd[pos]), zcrd_ovrf[pos] };
        }
        else {
          orep_xcrd = { static_cast<llint>(xcrd[pos]), 0 };
          orep_ycrd = { static_cast<llint>(ycrd[pos]), 0 };
          orep_zcrd = { static_cast<llint>(zcrd[pos]), 0 };
        }
        ixcrd = hostChangeFPBits(orep_xcrd, crd_scaling_bits, bgmr.dims.scale_bits);
        iycrd = hostChangeFPBits(orep_ycrd, crd_scaling_bits, bgmr.dims.scale_bits);
        izcrd = hostChangeFPBits(orep_zcrd, crd_scaling_bits, bgmr.dims.scale_bits);
      }
    }
    else {
      ixcrd = hostDoubleToInt95(xcrd[pos] * bgmr.dims.scale);
      iycrd = hostDoubleToInt95(ycrd[pos] * bgmr.dims.scale);
      izcrd = hostDoubleToInt95(zcrd[pos] * bgmr.dims.scale);
    }

    // Obtain the coordinates of the atom, relative to the mesh origin, in the precision of the
    // transformation matrices.  Estimate the appropriate mesh element, then test by computing
    // the location of the atom relative to the origin of this element.
    const int95_t ipt_relx = hostSplitFPSubtract(ixcrd, bgmr.dims.orig_x.x, bgmr.dims.orig_x.y);
    const int95_t ipt_rely = hostSplitFPSubtract(iycrd, bgmr.dims.orig_y.x, bgmr.dims.orig_y.y);
    const int95_t ipt_relz = hostSplitFPSubtract(izcrd, bgmr.dims.orig_z.x, bgmr.dims.orig_z.y);
    const Tcalc pt_relx = hostInt95ToDouble(ipt_relx) * bgmr.dims.inv_scale;
    const Tcalc pt_rely = hostInt95ToDouble(ipt_rely) * bgmr.dims.inv_scale;
    const Tcalc pt_relz = hostInt95ToDouble(ipt_relz) * bgmr.dims.inv_scale;
    Tcalc pt_grid_a = (mesh_umat[0] * pt_relx) + (mesh_umat[3] * pt_rely) +
                      (mesh_umat[6] * pt_relz);
    Tcalc pt_grid_b = (mesh_umat[4] * pt_rely) + (mesh_umat[7] * pt_relz);
    Tcalc pt_grid_c = (mesh_umat[8] * pt_relz);
    switch (bgmr.dims.bounds) {
    case BoundaryCondition::ISOLATED:
      switch (policy) {
      case OffMeshProtocol::DIE:
      case OffMeshProtocol::EXTRAPOLATE:
        break;
      case OffMeshProtocol::ZERO:

        // No contribution to the energy takes place.  Set applicable forces to zero.
        if (compute_forces) {
          xfrc[pos] = coord_value_zero;
          yfrc[pos] = coord_value_zero;
          zfrc[pos] = coord_value_zero;
          if (force_overflow_active) {
            xfrc_ovrf[pos] = 0;
            yfrc_ovrf[pos] = 0;
            zfrc_ovrf[pos] = 0;
          }
        }
      }
      break;
    case BoundaryCondition::PERIODIC:
      if (tcalc_is_double) {
        pt_grid_a -= floor(pt_grid_a / real_na) * real_na;
        pt_grid_b -= floor(pt_grid_b / real_nb) * real_nb;
        pt_grid_c -= floor(pt_grid_c / real_nc) * real_nc;
      }
      else {
        pt_grid_a -= floorf(pt_grid_a / real_na) * real_na;
        pt_grid_b -= floorf(pt_grid_b / real_nb) * real_nb;
        pt_grid_c -= floorf(pt_grid_c / real_nc) * real_nc;
      }
      break;
    }
    Tcalc fl_grid_a, fl_grid_b, fl_grid_c;
    if (tcalc_is_double) {
      fl_grid_a = floor(pt_grid_a);
      fl_grid_b = floor(pt_grid_b);
      fl_grid_c = floor(pt_grid_c);
    }
    else {
      fl_grid_a = floorf(pt_grid_a);
      fl_grid_b = floorf(pt_grid_b);
      fl_grid_c = floorf(pt_grid_c);
    }
    int cell_a = fl_grid_a;
    int cell_b = fl_grid_b;
    int cell_c = fl_grid_c;
    Tcalc frac_a = pt_grid_a - fl_grid_a;
    Tcalc frac_b = pt_grid_b - fl_grid_b;
    Tcalc frac_c = pt_grid_c - fl_grid_c;
    
    // Check for off-lattice points
    switch (bgmr.dims.bounds) {
    case BoundaryCondition::ISOLATED:
      switch (policy) {
      case OffMeshProtocol::DIE:
        if (cell_a < 0 || cell_a >= bgmr.dims.na || cell_b < 0 || cell_b >= bgmr.dims.nb ||
            cell_c < 0 || cell_c >= bgmr.dims.nc) {
          rtErr("Point [ " + realToString(hostInt95ToDouble(ixcrd) * bgmr.dims.inv_scale) + ", " +
                realToString(hostInt95ToDouble(iycrd) * bgmr.dims.inv_scale) + ", " +
                realToString(hostInt95ToDouble(izcrd) * bgmr.dims.inv_scale) + "] is off the mesh "
                "lattice.", "interpolate");
        }
        break;
      case OffMeshProtocol::EXTRAPOLATE:

        // Move the mesh element index to one that is in bounds, and adjust the fractional
        // coordinates as appropriate.
        if (cell_a < 0) {
          frac_a += static_cast<Tcalc>(cell_a);
          cell_a = 0;
        }
        else if (cell_a >= bgmr.dims.na - 1) {
          frac_a += static_cast<Tcalc>(cell_a - bgmr.dims.na + 1);
          cell_a = bgmr.dims.na - 1;
        }
        if (cell_b < 0) {
          frac_b += static_cast<Tcalc>(cell_b);
          cell_b = 0;
        }
        else if (cell_b >= bgmr.dims.nb - 1) {
          frac_b += static_cast<Tcalc>(cell_b - bgmr.dims.nb + 1);
          cell_b = bgmr.dims.nb - 1;
        }
        if (cell_c < 0) {
          frac_c += static_cast<Tcalc>(cell_c);
          cell_c = 0;
        }
        else if (cell_c >= bgmr.dims.nc - 1) {
          frac_c += static_cast<Tcalc>(cell_c - bgmr.dims.nc + 1);
          cell_c = bgmr.dims.nc - 1;
        }
        break;
      case OffMeshProtocol::ZERO:

        // No contribution to the energy takes place.  Set applicable forces to zero and move on.
        if (cell_a < 0 || cell_a >= bgmr.dims.na || cell_b < 0 || cell_b >= bgmr.dims.nb ||
            cell_c < 0 || cell_c >= bgmr.dims.nc) {
          if (compute_forces) {
            xfrc[pos] = coord_value_zero;
            yfrc[pos] = coord_value_zero;
            zfrc[pos] = coord_value_zero;
            if (force_overflow_active) {
              xfrc_ovrf[pos] = 0;
              yfrc_ovrf[pos] = 0;
              zfrc_ovrf[pos] = 0;
            }
          }
        }
        continue;
      }
      break;
    case BoundaryCondition::PERIODIC:
      break;
    }

    // Perform the interpolations and return the vector of results
    const size_t element_idx = (((cell_c * bgmr.dims.nb) + cell_b) * bgmr.dims.na) + cell_a;
    switch (bgmr.kind) {
    case GridDetail::OCCLUSION:
      {
        int cubelet_a_idx = frac_a * value_sixteen;
        int cubelet_b_idx = frac_b * value_sixteen;
        int cubelet_c_idx = frac_c * value_sixteen;
        const int cube_a_idx = cubelet_a_idx / 4;
        const int cube_b_idx = cubelet_b_idx / 4;
        const int cube_c_idx = cubelet_c_idx / 4;
        cubelet_a_idx -= cube_a_idx * 4;
        cubelet_b_idx -= cube_b_idx * 4;
        cubelet_c_idx -= cube_c_idx * 4;
        const size_t cube_idx = (((cube_c_idx * 4) + cube_b_idx) * 4) + cube_a_idx;
        const int cubelet_idx = (((cubelet_c_idx * 4) + cubelet_b_idx) * 4) + cubelet_a_idx;
        const ullint tmp_coef = bgmr.coeffs[(element_idx * sixty_four) + cube_idx];
        if ((tmp_coef >> cubelet_idx) & 0x1LLU) {

          // Add to the energy.  There are no forces for an occlusion mesh.
          nrg_acc += llround(bgmr.occ_cost * nrg_scale_factor);
        }
      }
      break;
    case GridDetail::OCCLUSION_FIELD:
    case GridDetail::NONBONDED_FIELD:
    case GridDetail::NONBONDED_ATOMIC:
      {
        const bool a_is_zero = (std::abs(frac_a) < constants::tiny);
        const bool b_is_zero = (std::abs(frac_b) < constants::tiny);
        const bool c_is_zero = (std::abs(frac_c) < constants::tiny);
        Tcalc tmp_sum = value_zero;

        // Force evaluation is triggered by the presence of an array to place the forces in.
        // If this is not present, only energies will be evaluated.  As with other molecular
        // mechanics computations in CPU host code, energy will always be evaluated.
        if (compute_forces) {
          Tcalc tmp_da  = value_zero;
          Tcalc tmp_db  = value_zero;
          Tcalc tmp_dc  = value_zero;
          Tcalc av = value_one;
          Tcalc d_av = (a_is_zero) ? value_zero : value_one / frac_a;
          for (int i = 0; i < 4; i++) {
            Tcalc bv = value_one;
            Tcalc d_bv = (b_is_zero) ? value_zero : value_one / frac_b;
            for (int j = 0; j < 4; j++) {
              Tcalc cv = value_one;
              Tcalc d_cv = (c_is_zero) ? value_zero : value_one / frac_c;
              for (int k = 0; k < 4; k++) {
                const size_t coef_idx = (((k * 4) + j) * 4) + i;
                const Tcalc tcoef = bgmr.coeffs[(sixty_four * element_idx) + coef_idx];
                tmp_sum += tcoef * av * bv * cv;
                tmp_da  += tcoef * dfactors[i] * d_av * bv * cv;
                tmp_db  += tcoef * dfactors[j] * av * d_bv * cv;
                tmp_dc  += tcoef * dfactors[k] * av * bv * d_cv;
                cv *= frac_c;
                if (c_is_zero) {
                  d_cv = (k == 0);
                }
                else {
                  d_cv *= frac_c;
                }
              }
              bv *= frac_b;
              if (b_is_zero) {
                d_bv = (j == 0);
              }
              else {
                d_bv *= frac_b;
              }
            }
            av *= frac_a;
            if (a_is_zero) {
              d_av = (i == 0);
            }
            else {
              d_av *= frac_a;
            }
          }

          // Convert the forces to appropriate units.  The first derivatives, continuous throughout
          // the tricubic mesh, must be interpreted through the chain rules for converting the
          // coordinate system in the a, b, and c axes of the mesh cell back to Cartesian
          // coordinates.
          switch (bgmr.field) {
          case NonbondedPotential::CLASH:
          case NonbondedPotential::VAN_DER_WAALS:

            // For a probe of specific Lennard-Jones parameters, the mesh itself contains the
            // exact value of the interaction.
            break;
          case NonbondedPotential::ELECTROSTATIC:

            // The charge is contained in the first property array.  The Coulomb constant must be
            // applied, and is conveyed by the mesh's non-bonded kit.
            tmp_sum *= prop_a[pos];
            break;
          }
          const Tcalc tmp_fx = (mesh_umat[0] * tmp_da);
          const Tcalc tmp_fy = (mesh_umat[3] * tmp_da) + (mesh_umat[4] * tmp_db);
          const Tcalc tmp_fz = (mesh_umat[6] * tmp_da) + (mesh_umat[7] * tmp_db) +
                               (mesh_umat[8] * tmp_dc);
          if (force_overflow_active) {
            const int95_t ixfrc = hostDoubleToInt95(tmp_da);
            const int95_t iyfrc = hostDoubleToInt95(tmp_db);
            const int95_t izfrc = hostDoubleToInt95(tmp_dc);
            const int95_t ipxfrc = hostSplitFPSum(ixfrc, xfrc[pos], xfrc_ovrf[pos]);
            const int95_t ipyfrc = hostSplitFPSum(iyfrc, yfrc[pos], yfrc_ovrf[pos]);
            const int95_t ipzfrc = hostSplitFPSum(izfrc, zfrc[pos], zfrc_ovrf[pos]);
            xfrc[pos] = ipxfrc.x;
            yfrc[pos] = ipyfrc.x;
            zfrc[pos] = ipzfrc.x;
            xfrc_ovrf[pos] = ipxfrc.y;
            yfrc_ovrf[pos] = ipyfrc.y;
            zfrc_ovrf[pos] = ipzfrc.y;
          }
          else {
            if (force_is_integral) {
              xfrc[pos] += llround(tmp_fx * force_scale);
              yfrc[pos] += llround(tmp_fy * force_scale);
              zfrc[pos] += llround(tmp_fz * force_scale);
            }
            else {
              xfrc[pos] += tmp_fx;
              yfrc[pos] += tmp_fy;
              zfrc[pos] += tmp_fz;
            }
          }
        }
        else {
          Tcalc av = value_one;
          for (int i = 0; i < 4; i++) {
            Tcalc bv = value_one;
            for (int j = 0; j < 4; j++) {
              Tcalc cv = value_one;
              for (int k = 0; k < 4; k++) {
                const size_t coef_idx = (((k * 4) + j) * 4) + i;
                const Tcalc tcoef = bgmr.coeffs[(sixty_four * element_idx) + coef_idx];
                tmp_sum += tcoef * av * bv * cv;
                cv *= frac_c;
              }
              bv *= frac_b;
            }
            av *= frac_a;
          }
        }
        
        // The function value requires no chain rules to interpret.  Accumulate the energy.
        const llint itmp_sum = llround(tmp_sum * nrg_scale_factor);
        nrg_acc += itmp_sum;
        if (vnrg != nullptr) {
          if (obs_is_integral) {
            vnrg[pos] = itmp_sum;
          }
          else {
            vnrg[pos] = tmp_sum;
          }
        }
      }
      break;
    }
  }

  // Return the accumulated energy in units of kcal/mol
  return (static_cast<double>(nrg_acc) / nrg_scale_factor);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void interpolate(const BackgroundMeshReader<Tcalc> &bgmr, const MeshFFKit<Tcalc> &mnbk,
                 const NonbondedKit<Tcalc> &nbk, PhaseSpaceWriter *psw, ScoreCard *sc,
                 const OffMeshProtocol policy, const int system_index) {
  double nrg_acc;
  switch (bgmr.field) {
  case NonbondedPotential::CLASH:
    break;
  case NonbondedPotential::ELECTROSTATIC:
    nrg_acc = interpolate<Tcalc, double>(bgmr, mnbk, psw->xcrd, psw->ycrd, psw->zcrd, nullptr,
                                         nullptr, nullptr, nbk.charge, nullptr, nullptr,
                                         psw->natom, nullptr, psw->xfrc, psw->yfrc, psw->zfrc,
                                         nullptr, nullptr, nullptr, policy, system_index, 0, 0,
                                         sc->getEnergyScaleBits());
    break;
  case NonbondedPotential::VAN_DER_WAALS:
    nrg_acc = interpolate<Tcalc, double>(bgmr, mnbk, psw->xcrd, psw->ycrd, psw->zcrd, nullptr,
                                         nullptr, nullptr, nbk.lja_coeff, nbk.ljb_coeff,
                                         nbk.lj_idx, psw->natom, nullptr, psw->xfrc, psw->yfrc,
                                         psw->zfrc, nullptr, nullptr, nullptr, policy,
                                         system_index, 0, 0, sc->getEnergyScaleBits());
    break;
  }

  // Store the accumulated energy in the appropriate component of the tracking object.
  switch (bgmr.field) {
  case NonbondedPotential::CLASH:
    sc->contribute(StateVariable::VDW, nrg_acc, system_index);
    break;
  case NonbondedPotential::ELECTROSTATIC:
    sc->contribute(StateVariable::ELECTROSTATIC, nrg_acc, system_index);
    break;
  case NonbondedPotential::VAN_DER_WAALS:
    sc->contribute(StateVariable::VDW, nrg_acc, system_index);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcoord>
double footprint(const BackgroundMeshReader<ullint> &bgmr, const Tcoord* xcrd, const Tcoord* ycrd,
                 const Tcoord* zcrd, const int* xcrd_ovrf, const int* ycrd_ovrf,
                 const int* zcrd_ovrf, const uint* prop_idx, const int natom, const Tcoord* vclash,
                 const OffMeshProtocol policy, const int system_index,
                 const int crd_scaling_bits) {

  // Create a mock MeshFFKit as a placeholder
  const MeshFFKit<Tcalc> mnbk;
  return interpolate<ullint, Tcalc, Tcoord>(bgmr, mnbk, xcrd, ycrd, zcrd, xcrd_ovrf, ycrd_ovrf,
                                            zcrd_ovrf, nullptr, nullptr,
                                            reinterpret_cast<const uint*>(prop_idx), natom, vclash,
                                            nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                            policy, system_index, crd_scaling_bits, 0, 0);
}

} // namespace structure
} // namespace stormm
