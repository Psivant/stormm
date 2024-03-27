// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
Tcalc evalHarmonicStretch(const int i_atom, const int j_atom, const Tcalc stiffness,
                          const Tcalc equilibrium, const Tcoord* xcrd, const Tcoord* ycrd,
                          const Tcoord* zcrd, const double* umat, const double* invu,
                          const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                          const EvaluateForce eval_force, const Tcalc inv_gpos_factor,
                          const Tcalc force_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  Tcalc dx, dy, dz;
  if (isSignedIntegralScalarType<Tcoord>()) {
    dx = static_cast<Tcalc>(xcrd[j_atom] - xcrd[i_atom]) * inv_gpos_factor;
    dy = static_cast<Tcalc>(ycrd[j_atom] - ycrd[i_atom]) * inv_gpos_factor;
    dz = static_cast<Tcalc>(zcrd[j_atom] - zcrd[i_atom]) * inv_gpos_factor;
  }
  else {
    dx = xcrd[j_atom] - xcrd[i_atom];
    dy = ycrd[j_atom] - ycrd[i_atom];
    dz = zcrd[j_atom] - zcrd[i_atom];
  }
  imageCoordinates<Tcalc, Tcalc>(&dx, &dy, &dz, umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);
  const Tcalc dr = (tcalc_ct == double_type_index) ? sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                                                     sqrtf((dx * dx) + (dy * dy) + (dz * dz));
  const Tcalc dl = dr - equilibrium;
  
  // Compute forces
  if (eval_force == EvaluateForce::YES) {
    const Tcalc fmag = 2.0 * stiffness * dl / dr;
    const Tcalc fmag_dx = fmag * dx;
    const Tcalc fmag_dy = fmag * dy;
    const Tcalc fmag_dz = fmag * dz;
    if (isSignedIntegralScalarType<Tforce>()) {
      const Tforce ifmag_dx = llround(fmag_dx * force_factor);
      const Tforce ifmag_dy = llround(fmag_dy * force_factor);
      const Tforce ifmag_dz = llround(fmag_dz * force_factor);
      xfrc[i_atom] += ifmag_dx;
      yfrc[i_atom] += ifmag_dy;
      zfrc[i_atom] += ifmag_dz;
      xfrc[j_atom] -= ifmag_dx;
      yfrc[j_atom] -= ifmag_dy;
      zfrc[j_atom] -= ifmag_dz;
    }
    else {
      xfrc[i_atom] += fmag_dx;
      yfrc[i_atom] += fmag_dy;
      zfrc[i_atom] += fmag_dz;
      xfrc[j_atom] -= fmag_dx;
      yfrc[j_atom] -= fmag_dy;
      zfrc[j_atom] -= fmag_dz;
    }
  }
  return stiffness * dl * dl;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
double evaluateBondTerms(const ValenceKit<Tcalc> vk, const Tcoord* xcrd, const Tcoord* ycrd,
                         const Tcoord* zcrd, const double* umat, const double* invu,
                         const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                         ScoreCard *ecard, const EvaluateForce eval_force, const int system_index,
                         const Tcalc inv_gpos_factor, const Tcalc force_factor) {

  // Use two accumulators: one, a standard double-precision accumulator and the other a
  // fixed-precision long long integer with discretization at the global energy scaling factor
  // (see Constants/fixed_precision.h).
  double bond_energy = 0.0;
  llint bond_acc = 0LL;
  const Tcalc nrg_scale_factor = ecard->getEnergyScalingFactor<Tcalc>();

  // Accumulate the results (energy in both precision models)
  for (int pos = 0; pos < vk.nbond; pos++) {
    const int param_idx = vk.bond_param_idx[pos];
    const double du =
      evalHarmonicStretch<Tcoord, Tforce, Tcalc>(vk.bond_i_atoms[pos], vk.bond_j_atoms[pos],
                                                 vk.bond_keq[param_idx],
                                                 fabs(vk.bond_leq[param_idx]),
                                                 xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                                 yfrc, zfrc, eval_force, inv_gpos_factor,
                                                 force_factor);
    bond_energy += du;
    bond_acc += llround(du * nrg_scale_factor);
  }
  
  // Contribute results
  ecard->contribute(StateVariable::BOND, bond_acc, system_index);
  
  // Return the double-precision bond energy sum, if of interest
  return bond_energy;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double evaluateBondTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesReader<Tcoord> csr,
                         ScoreCard *ecard, const int system_index, const int force_scale_bits) {
  const size_t atom_os = static_cast<size_t>(system_index) *
                         roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t xfrm_os = static_cast<size_t>(system_index) * roundUp<size_t>(9, warp_size_zu);

  // Compute the force scaling factor based on the requested bit count.  While this routine will
  // not compute forces per se, the scaling factor is still relevant to the accumulators for Born
  // radii.
  const Tcalc force_scale = (isSignedIntegralScalarType<Tcoord>()) ?
                            pow(2.0, force_scale_bits) : 1.0;
  return evaluateBondTerms<Tcoord, Tcoord, Tcalc>(vk, &csr.xcrd[atom_os], &csr.ycrd[atom_os],
                                                  &csr.zcrd[atom_os], &csr.umat[xfrm_os],
                                                  &csr.invu[xfrm_os], csr.unit_cell, nullptr,
                                                  nullptr, nullptr, ecard, EvaluateForce::NO,
                                                  system_index, csr.inv_gpos_scale, force_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double evaluateBondTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesWriter<Tcoord> csw,
                         ScoreCard *ecard, const int system_index, const int force_scale_bits) {
  return evaluateBondTerms<Tcoord, Tcalc>(vk, CoordinateSeriesReader(csw), ecard, system_index,
                                          force_scale_bits);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
Tcalc evalHarmonicBend(const int i_atom, const int j_atom, const int k_atom,
                       const Tcalc stiffness, const Tcalc equilibrium, const Tcoord* xcrd,
                       const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                       const double* invu, const UnitCellType unit_cell, Tforce* xfrc,
                       Tforce* yfrc, Tforce* zfrc, const EvaluateForce eval_force,
                       const Tcalc inv_gpos_factor, const Tcalc force_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc value_one = 1.0;
  
  // Compute displacements
  Tcalc ba[3], bc[3];
  if (isSignedIntegralScalarType<Tcoord>()) {
    ba[0] = static_cast<Tcalc>(xcrd[i_atom] - xcrd[j_atom]) * inv_gpos_factor;
    ba[1] = static_cast<Tcalc>(ycrd[i_atom] - ycrd[j_atom]) * inv_gpos_factor;
    ba[2] = static_cast<Tcalc>(zcrd[i_atom] - zcrd[j_atom]) * inv_gpos_factor;
    bc[0] = static_cast<Tcalc>(xcrd[k_atom] - xcrd[j_atom]) * inv_gpos_factor;
    bc[1] = static_cast<Tcalc>(ycrd[k_atom] - ycrd[j_atom]) * inv_gpos_factor;
    bc[2] = static_cast<Tcalc>(zcrd[k_atom] - zcrd[j_atom]) * inv_gpos_factor;
  }
  else {
    ba[0] = xcrd[i_atom] - xcrd[j_atom];
    ba[1] = ycrd[i_atom] - ycrd[j_atom];
    ba[2] = zcrd[i_atom] - zcrd[j_atom];
    bc[0] = xcrd[k_atom] - xcrd[j_atom];
    bc[1] = ycrd[k_atom] - ycrd[j_atom];
    bc[2] = zcrd[k_atom] - zcrd[j_atom];
  }
  imageCoordinates<Tcalc, Tcalc>(&ba[0], &ba[1], &ba[2], umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates<Tcalc, Tcalc>(&bc[0], &bc[1], &bc[2], umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);

  // On to the angle force computation
  const Tcalc mgba = (ba[0] * ba[0]) + (ba[1] * ba[1]) + (ba[2] * ba[2]);
  const Tcalc mgbc = (bc[0] * bc[0]) + (bc[1] * bc[1]) + (bc[2] * bc[2]);
  const Tcalc invbabc = (tcalc_is_double) ? 1.0 / sqrt(mgba * mgbc) :
                                            value_one / sqrtf(mgba * mgbc);
  Tcalc costheta = ((ba[0] * bc[0]) + (ba[1] * bc[1]) + (ba[2] * bc[2])) * invbabc;
  costheta = (costheta < -value_one) ? -value_one : (costheta > value_one) ? value_one : costheta;
  const Tcalc theta = (tcalc_is_double) ? acos(costheta) : acosf(costheta);
  const Tcalc dtheta = theta - equilibrium;

  // Compute forces
  if (eval_force == EvaluateForce::YES) {
    const Tcalc dA = (tcalc_is_double) ?
                     -2.0 * stiffness * dtheta / sqrt(1.0 - (costheta * costheta)) :
                     -2.0f * stiffness * dtheta / sqrtf(value_one - (costheta * costheta));
    const Tcalc sqba = dA / mgba;
    const Tcalc sqbc = dA / mgbc;
    const Tcalc mbabc = dA * invbabc;
    if (isSignedIntegralScalarType<Tforce>()) {
      Tforce iadf[3], icdf[3];
      for (int i = 0; i < 3; i++) {
        iadf[i] = llround(((bc[i] * mbabc) - (costheta * ba[i] * sqba)) * force_factor);
        icdf[i] = llround(((ba[i] * mbabc) - (costheta * bc[i] * sqbc)) * force_factor);
      }      
      xfrc[i_atom] -= iadf[0];
      yfrc[i_atom] -= iadf[1];
      zfrc[i_atom] -= iadf[2];
      xfrc[j_atom] += iadf[0] + icdf[0];
      yfrc[j_atom] += iadf[1] + icdf[1];
      zfrc[j_atom] += iadf[2] + icdf[2];
      xfrc[k_atom] -= icdf[0];
      yfrc[k_atom] -= icdf[1];
      zfrc[k_atom] -= icdf[2];
    }
    else {
      Tcalc adf[3], cdf[3];
      for (int i = 0; i < 3; i++) {
        adf[i] = (bc[i] * mbabc) - (costheta * ba[i] * sqba);
        cdf[i] = (ba[i] * mbabc) - (costheta * bc[i] * sqbc);
      }
      xfrc[i_atom] -= adf[0];
      yfrc[i_atom] -= adf[1];
      zfrc[i_atom] -= adf[2];
      xfrc[j_atom] += adf[0] + cdf[0];
      yfrc[j_atom] += adf[1] + cdf[1];
      zfrc[j_atom] += adf[2] + cdf[2];
      xfrc[k_atom] -= cdf[0];
      yfrc[k_atom] -= cdf[1];
      zfrc[k_atom] -= cdf[2];
    }
  }
  return stiffness * dtheta * dtheta;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
double evaluateAngleTerms(const ValenceKit<Tcalc> vk, const Tcoord* xcrd, const Tcoord* ycrd,
                          const Tcoord* zcrd, const double* umat, const double* invu,
                          const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                          ScoreCard *ecard, const EvaluateForce eval_force, const int system_index,
                          const Tcalc inv_gpos_factor, const Tcalc force_factor) {

  // As in the bond interaction computation above, use two accumulators.  The full double-precision
  // result will be returned if there is interest in comparing to the fixed-precision accumulation.
  // Also as above, the virial is not computed in full double precision; that quantity may affect
  // the way the system changes volume, but because its accumulation is spread over all of the
  // system's interactions and it in turn affects the motion of all particles collectively, it is
  // not of as much interest as the error in the fixed-precision net force on a single particle.
  double angl_energy = 0.0;
  llint angl_acc = 0LL;
  const Tcalc nrg_scale_factor = ecard->getEnergyScalingFactor<Tcalc>();

  // Accumulate results by looping over all angle bending terms.
  for (int pos = 0; pos < vk.nangl; pos++) {
    const int param_idx = vk.angl_param_idx[pos];
    const double du =
      evalHarmonicBend<Tcoord, Tforce, Tcalc>(vk.angl_i_atoms[pos], vk.angl_j_atoms[pos],
                                              vk.angl_k_atoms[pos], vk.angl_keq[param_idx],
                                              vk.angl_theta[param_idx], xcrd, ycrd, zcrd, umat,
                                              invu, unit_cell, xfrc, yfrc, zfrc, eval_force,
                                              inv_gpos_factor, force_factor);
    angl_energy += du;
    angl_acc += llround(du * nrg_scale_factor);
  }

  // Contribute results
  ecard->contribute(StateVariable::ANGLE, angl_acc, system_index);

  // Return the double-precision angle energy sum, if of interest
  return angl_energy;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double evaluateAngleTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesReader<Tcoord> csr,
                          ScoreCard *ecard, const int system_index, const int force_scale_bits) {
  const size_t atom_os = static_cast<size_t>(system_index) *
                         roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t xfrm_os = static_cast<size_t>(system_index) * roundUp<size_t>(9, warp_size_zu);
  const Tcalc force_scale = (isSignedIntegralScalarType<Tcoord>()) ?
                            pow(2.0, force_scale_bits) : 1.0;
  return evaluateAngleTerms<Tcoord, Tcoord, Tcalc>(vk, &csr.xcrd[atom_os], &csr.ycrd[atom_os],
                                                   &csr.zcrd[atom_os], &csr.umat[xfrm_os],
                                                   &csr.invu[xfrm_os], csr.unit_cell, nullptr,
                                                   nullptr, nullptr, ecard, EvaluateForce::NO,
                                                   system_index, csr.inv_gpos_scale, force_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double evaluateAngleTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesWriter<Tcoord> csw,
                          ScoreCard *ecard, const int system_index, const int force_scale_bits) {
  return evaluateAngleTerms<Tcoord, Tcalc>(vk, CoordinateSeriesReader(csw), ecard, system_index,
                                           force_scale_bits);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
Tcalc evalDihedralTwist(const int i_atom, const int j_atom, const int k_atom, const int l_atom,
                        const Tcalc amplitude, const Tcalc phase_angle, const Tcalc frequency,
                        const DihedralStyle kind, const Tcoord* xcrd, const Tcoord* ycrd,
                        const Tcoord* zcrd, const double* umat, const double* invu,
                        const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                        const EvaluateForce eval_force, const Tcalc inv_gpos_factor,
                        const Tcalc force_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc value_one = 1.0;

  // Compute displacements
  Tcalc ab[3], bc[3], cd[3], crabbc[3], crbccd[3], scr[3];
  if (isSignedIntegralScalarType<Tcoord>()) {
    ab[0] = static_cast<Tcalc>(xcrd[j_atom] - xcrd[i_atom]) * inv_gpos_factor;
    ab[1] = static_cast<Tcalc>(ycrd[j_atom] - ycrd[i_atom]) * inv_gpos_factor;
    ab[2] = static_cast<Tcalc>(zcrd[j_atom] - zcrd[i_atom]) * inv_gpos_factor;
    bc[0] = static_cast<Tcalc>(xcrd[k_atom] - xcrd[j_atom]) * inv_gpos_factor;
    bc[1] = static_cast<Tcalc>(ycrd[k_atom] - ycrd[j_atom]) * inv_gpos_factor;
    bc[2] = static_cast<Tcalc>(zcrd[k_atom] - zcrd[j_atom]) * inv_gpos_factor;
    cd[0] = static_cast<Tcalc>(xcrd[l_atom] - xcrd[k_atom]) * inv_gpos_factor;
    cd[1] = static_cast<Tcalc>(ycrd[l_atom] - ycrd[k_atom]) * inv_gpos_factor;
    cd[2] = static_cast<Tcalc>(zcrd[l_atom] - zcrd[k_atom]) * inv_gpos_factor;
  }
  else {
    ab[0] = xcrd[j_atom] - xcrd[i_atom];
    ab[1] = ycrd[j_atom] - ycrd[i_atom];
    ab[2] = zcrd[j_atom] - zcrd[i_atom];
    bc[0] = xcrd[k_atom] - xcrd[j_atom];
    bc[1] = ycrd[k_atom] - ycrd[j_atom];
    bc[2] = zcrd[k_atom] - zcrd[j_atom];
    cd[0] = xcrd[l_atom] - xcrd[k_atom];
    cd[1] = ycrd[l_atom] - ycrd[k_atom];
    cd[2] = zcrd[l_atom] - zcrd[k_atom];
  }
  imageCoordinates<Tcalc, Tcalc>(&ab[0], &ab[1], &ab[2], umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates<Tcalc, Tcalc>(&bc[0], &bc[1], &bc[2], umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates<Tcalc, Tcalc>(&cd[0], &cd[1], &cd[2], umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);

  // Compute cross products and then the angle between the planes
  crossProduct(ab, bc, crabbc);
  crossProduct(bc, cd, crbccd);
  Tcalc costheta = (crabbc[0] * crbccd[0]) + (crabbc[1] * crbccd[1]) + (crabbc[2] * crbccd[2]);
  if (tcalc_is_double) {
    costheta /= sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                     (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  }
  else {
    costheta /= sqrtf((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                      (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  }
  crossProduct(crabbc, crbccd, scr);
  costheta = (costheta < -value_one) ? -value_one : (costheta > value_one) ? value_one : costheta;
  const Tcalc theta = angleVerification(costheta, crabbc, crbccd, bc, scr);
  Tcalc sangle;
  switch (kind) {
  case DihedralStyle::COSINE:
    sangle = (frequency * theta) - phase_angle;
    break;
  case DihedralStyle::HARMONIC:
    sangle = theta - phase_angle;
    break;
  }

  // Compute forces, if requested
  if (eval_force == EvaluateForce::YES) {
    Tcalc fr;
    switch (kind) {
    case DihedralStyle::COSINE:
      fr = (tcalc_is_double) ? amplitude * frequency * sin(sangle) :
                               amplitude * frequency * sinf(sangle);
      break;
    case DihedralStyle::HARMONIC:
      fr = (tcalc_is_double) ? -2.0 * amplitude * sangle : -2.0f * amplitude * sangle;
      break;
    }
    Tcalc mgab, mgbc, mgcd;
    if (tcalc_is_double) {
      mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
    }
    else {
      mgab = sqrtf(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      mgbc = sqrtf(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      mgcd = sqrtf(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
    }
    const Tcalc invab = value_one / mgab;
    const Tcalc invbc = value_one / mgbc;
    const Tcalc invcd = value_one / mgcd;
    const Tcalc cosb = -(ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2]) * invab * invbc;
    const Tcalc cosc = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2]) * invbc * invcd;
    Tcalc isinb2, isinc2;
    if (tcalc_is_double) {
      isinb2 = (cosb * cosb < asymptotic_to_one_lf) ?
               fr / (1.0 - (cosb * cosb)) : fr * inverse_one_minus_asymptote_lf;
      isinc2 = (cosc * cosc < asymptotic_to_one_lf) ?
               fr / (1.0 - (cosc * cosc)) : fr * inverse_one_minus_asymptote_lf;
    }
    else {
      isinb2 = (cosb * cosb < asymptotic_to_one_f) ?
               fr / (value_one - (cosb * cosb)) : fr * inverse_one_minus_asymptote_f;
      isinc2 = (cosc * cosc < asymptotic_to_one_f) ?
               fr / (value_one - (cosc * cosc)) : fr * inverse_one_minus_asymptote_f;
    }
    const Tcalc invabc = invab * invbc;
    const Tcalc invbcd = invbc * invcd;
    for (int i = 0; i < 3; i++) {
      crabbc[i] *= invabc;
      crbccd[i] *= invbcd;
    }

    // Transform the rotational derivatives to cartesian coordinates
    const Tcalc fa = -invab * isinb2;
    const Tcalc fb1 = (mgbc - (mgab * cosb)) * invabc * isinb2;
    const Tcalc fb2 = cosc * invbc * isinc2;
    const Tcalc fc1 = (mgbc - (mgcd * cosc)) * invbcd * isinc2;
    const Tcalc fc2 = cosb * invbc * isinb2;
    const Tcalc fd = -invcd * isinc2;
    if (isSignedIntegralScalarType<Tforce>()) {
      const Tforce ifrc_ix = llround(crabbc[0] * fa * force_factor);
      const Tforce ifrc_jx = llround(((fb1 * crabbc[0]) - (fb2 * crbccd[0])) * force_factor);
      const Tforce ifrc_lx = llround(-fd * crbccd[0] * force_factor);
      xfrc[i_atom] += ifrc_ix;
      xfrc[j_atom] += ifrc_jx;
      xfrc[k_atom] -= ifrc_ix + ifrc_jx + ifrc_lx;
      xfrc[l_atom] += ifrc_lx;
      const Tforce ifrc_iy = llround(crabbc[1] * fa * force_factor);
      const Tforce ifrc_jy = llround(((fb1 * crabbc[1]) - (fb2 * crbccd[1])) * force_factor);
      const Tforce ifrc_ly = llround(-fd * crbccd[1] * force_factor);
      yfrc[i_atom] += ifrc_iy;
      yfrc[j_atom] += ifrc_jy;
      yfrc[k_atom] -= ifrc_iy + ifrc_jy + ifrc_ly;
      yfrc[l_atom] += ifrc_ly;
      const Tforce ifrc_iz = llround(crabbc[2] * fa * force_factor);
      const Tforce ifrc_jz = llround(((fb1 * crabbc[2]) - (fb2 * crbccd[2])) * force_factor);
      const Tforce ifrc_lz = llround(-fd * crbccd[2] * force_factor);
      zfrc[i_atom] += ifrc_iz;
      zfrc[j_atom] += ifrc_jz;
      zfrc[k_atom] -= ifrc_iz + ifrc_jz + ifrc_lz;
      zfrc[l_atom] += ifrc_lz;
    }
    else {
      const Tforce frc_ix = crabbc[0] * fa;
      const Tforce frc_jx = (fb1 * crabbc[0]) - (fb2 * crbccd[0]);
      const Tforce frc_lx = -fd * crbccd[0];
      xfrc[i_atom] += frc_ix;
      xfrc[j_atom] += frc_jx;
      xfrc[k_atom] -= frc_ix + frc_jx + frc_lx;
      xfrc[l_atom] += frc_lx;
      const Tforce frc_iy = crabbc[1] * fa;
      const Tforce frc_jy = (fb1 * crabbc[1]) - (fb2 * crbccd[1]);
      const Tforce frc_ly = -fd * crbccd[1];
      yfrc[i_atom] += frc_iy;
      yfrc[j_atom] += frc_jy;
      yfrc[k_atom] -= frc_iy + frc_jy + frc_ly;
      yfrc[l_atom] += frc_ly;
      const Tforce frc_iz = crabbc[2] * fa;
      const Tforce frc_jz = (fb1 * crabbc[2]) - (fb2 * crbccd[2]);
      const Tforce frc_lz = -fd * crbccd[2];
      zfrc[i_atom] += frc_iz;
      zfrc[j_atom] += frc_jz;
      zfrc[k_atom] -= frc_iz + frc_jz + frc_lz;
      zfrc[l_atom] += frc_lz;
    }
  }
  switch (kind) {
  case DihedralStyle::COSINE:
    return (tcalc_is_double) ? amplitude * (1.0 + cos(sangle)) :
                               amplitude * (value_one + cosf(sangle));
  case DihedralStyle::HARMONIC:
    return amplitude * sangle * sangle;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
double2 evaluateDihedralTerms(const ValenceKit<Tcalc> vk, const Tcoord* xcrd, const Tcoord* ycrd,
                              const Tcoord* zcrd, const double* umat, const double* invu,
                              const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                              Tforce* zfrc, ScoreCard *ecard, const EvaluateForce eval_force,
                              const int system_index, const Tcalc inv_gpos_factor,
                              const Tcalc force_factor) {
  double2 dihe_energy = { 0.0, 0.0 };
  llint proper_acc   = 0LL;
  llint improper_acc = 0LL;
  const Tcalc nrg_scale_factor = ecard->getEnergyScalingFactor<Tcalc>();

  // Accumulate results by looping over all dihedral terms.
  for (int pos = 0; pos < vk.ndihe; pos++) {
    const int param_idx = vk.dihe_param_idx[pos];
    const double du =
      evalDihedralTwist<Tcoord, Tforce, Tcalc>(vk.dihe_i_atoms[pos], vk.dihe_j_atoms[pos],
                                               vk.dihe_k_atoms[pos], vk.dihe_l_atoms[pos],
                                               vk.dihe_amp[param_idx], vk.dihe_phi[param_idx],
                                               vk.dihe_freq[param_idx], DihedralStyle::COSINE,
                                               xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                               yfrc, zfrc, eval_force, inv_gpos_factor,
                                               force_factor);

    // Contribute the result to the correct pile: proper or improper
    const TorsionKind kind = static_cast<TorsionKind>(vk.dihe_modifiers[pos].w);
    switch (kind) {
    case TorsionKind::PROPER:
    case TorsionKind::PROPER_NO_14:
      dihe_energy.x += du;
      proper_acc += llround(du * nrg_scale_factor);
      break;
    case TorsionKind::IMPROPER:
    case TorsionKind::IMPROPER_NO_14:
      dihe_energy.y += du;
      improper_acc += llround(du * nrg_scale_factor);
      break;
    }
  }

  // Contribute results
  ecard->contribute(StateVariable::PROPER_DIHEDRAL, proper_acc, system_index);
  ecard->contribute(StateVariable::IMPROPER_DIHEDRAL, improper_acc, system_index);

  // Return the double-precision dihedral energy sum, if of interest
  return dihe_energy;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double2 evaluateDihedralTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesReader<Tcoord> csr,
                              ScoreCard *ecard, const int system_index,
                              const int force_scale_bits) {
  const size_t atom_os = static_cast<size_t>(system_index) *
                         roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t xfrm_os = static_cast<size_t>(system_index) * roundUp<size_t>(9, warp_size_zu);
  const Tcalc force_scale = (isSignedIntegralScalarType<Tcoord>()) ?
                            pow(2.0, force_scale_bits) : 1.0;
  return evaluateDihedralTerms<Tcoord, Tcoord, Tcalc>(vk, &csr.xcrd[atom_os], &csr.ycrd[atom_os],
                                                      &csr.zcrd[atom_os], &csr.umat[xfrm_os],
                                                      &csr.invu[xfrm_os], csr.unit_cell, nullptr,
                                                      nullptr, nullptr, ecard, EvaluateForce::NO,
                                                      system_index, csr.inv_gpos_scale,
                                                      force_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double2 evaluateDihedralTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesWriter<Tcoord> csw,
                              ScoreCard *ecard, const int system_index,
                              const int force_scale_bits) {
  return evaluateDihedralTerms<Tcoord, Tcalc>(vk, CoordinateSeriesReader(csw), ecard, system_index,
                                              force_scale_bits);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
double evaluateUreyBradleyTerms(const ValenceKit<Tcalc> vk, const Tcoord* xcrd,
                                const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                                const double* invu, const UnitCellType unit_cell, Tforce* xfrc,
                                Tforce* yfrc, Tforce* zfrc, ScoreCard *ecard,
                                const EvaluateForce eval_force, const int system_index,
                                const Tcalc inv_gpos_factor, const Tcalc force_factor) {

  // Accumulate the results in two numerical precision models by looping over all terms
  double ubrd_energy = 0.0;
  llint ubrd_acc = 0LL;
  const Tcalc nrg_scale_factor = ecard->getEnergyScalingFactor<Tcalc>();

  // Accumulate the results (energy in both precision models)
  for (int pos = 0; pos < vk.nubrd; pos++) {
    const int param_idx = vk.ubrd_param_idx[pos];
    const double du =
      evalHarmonicStretch<Tcoord, Tforce, Tcalc>(vk.ubrd_i_atoms[pos], vk.ubrd_k_atoms[pos],
                                                 vk.ubrd_keq[param_idx], vk.ubrd_leq[param_idx],
                                                 xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc,
                                                 yfrc, zfrc, eval_force, inv_gpos_factor,
                                                 force_factor);
    ubrd_energy += du;
    ubrd_acc += llround(du * nrg_scale_factor);
  }

  // Contribute results
  ecard->contribute(StateVariable::UREY_BRADLEY, ubrd_acc, system_index);

  // Return the double-precision Urey-Bradley energy sum, if of interest
  return ubrd_energy;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double evaluateUreyBradleyTerms(const ValenceKit<Tcalc> vk,
                                const CoordinateSeriesReader<Tcoord> csr, ScoreCard *ecard,
                                const int system_index, const int force_scale_bits) {
  const size_t atom_os = static_cast<size_t>(system_index) *
                         roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t xfrm_os = static_cast<size_t>(system_index) * roundUp<size_t>(9, warp_size_zu);
  const Tcalc force_scale = (isSignedIntegralScalarType<Tcoord>()) ?
                            pow(2.0, force_scale_bits) : 1.0;
  return evaluateUreyBradleyTerms<Tcoord,
                                  Tcoord, Tcalc>(vk, &csr.xcrd[atom_os], &csr.ycrd[atom_os],
                                                 &csr.zcrd[atom_os], &csr.umat[xfrm_os],
                                                 &csr.invu[xfrm_os], csr.unit_cell, nullptr,
                                                 nullptr, nullptr, ecard, EvaluateForce::NO,
                                                 system_index, csr.inv_gpos_scale, force_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double evaluateUreyBradleyTerms(const ValenceKit<Tcalc> vk,
                                const CoordinateSeriesWriter<Tcoord> csw, ScoreCard *ecard,
                                const int system_index, const int force_scale_bits) {
  return evaluateUreyBradleyTerms<Tcoord, Tcalc>(vk, CoordinateSeriesReader(csw), ecard,
                                                 system_index, force_scale_bits);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
double evaluateCharmmImproperTerms(const ValenceKit<Tcalc> vk, const Tcoord* xcrd,
                                   const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                                   const double* invu, const UnitCellType unit_cell, Tforce* xfrc,
                                   Tforce* yfrc, Tforce* zfrc, ScoreCard *ecard,
                                   const EvaluateForce eval_force, const int system_index,
                                   const Tcalc inv_gpos_factor, const Tcalc force_factor) {
  double cimp_energy = 0.0;
  llint cimp_acc = 0LL;
  const Tcalc nrg_scale_factor = ecard->getEnergyScalingFactor<Tcalc>();

  // Accumulate results by looping over all CHARMM improper terms.
  for (int pos = 0; pos < vk.ncimp; pos++) {

    // Get parameters for an angle between atoms i, j, and k
    const int param_idx = vk.cimp_param_idx[pos];
    const double du =
      evalDihedralTwist<Tcoord, Tforce, Tcalc>(vk.cimp_i_atoms[pos], vk.cimp_j_atoms[pos],
                                               vk.cimp_k_atoms[pos], vk.cimp_l_atoms[pos],
                                               vk.cimp_keq[param_idx], vk.cimp_phi[param_idx],
                                               1.0, DihedralStyle::HARMONIC, xcrd, ycrd, zcrd,
                                               umat, invu, unit_cell, xfrc, yfrc, zfrc,
                                               eval_force, inv_gpos_factor, force_factor);
    cimp_energy += du;
    cimp_acc += llround(du * nrg_scale_factor);
  }

  // Contribute results
  ecard->contribute(StateVariable::CHARMM_IMPROPER, cimp_acc, system_index);

  // Return the double-precision CHARMM improper energy sum, if of interest
  return cimp_energy;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double evaluateCharmmImproperTerms(const ValenceKit<Tcalc> vk,
                                   const CoordinateSeriesReader<Tcoord> csr, ScoreCard *ecard,
                                   const int system_index, const int force_scale_bits) {
  const size_t atom_os = static_cast<size_t>(system_index) *
                         roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t xfrm_os = static_cast<size_t>(system_index) * roundUp<size_t>(9, warp_size_zu);
  const Tcalc force_scale = (isSignedIntegralScalarType<Tcoord>()) ?
                            pow(2.0, force_scale_bits) : 1.0;
  return evaluateCharmmImproperTerms<Tcoord,
                                     Tcoord, Tcalc>(vk, &csr.xcrd[atom_os], &csr.ycrd[atom_os],
                                                    &csr.zcrd[atom_os], &csr.umat[xfrm_os],
                                                    &csr.invu[xfrm_os], csr.unit_cell, nullptr,
                                                    nullptr, nullptr, ecard, EvaluateForce::NO,
                                                    system_index, csr.inv_gpos_scale, force_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double evaluateCharmmImproperTerms(const ValenceKit<Tcalc> vk,
                                   const CoordinateSeriesWriter<Tcoord> csw, ScoreCard *ecard,
                                   const int system_index, const int force_scale_bits) {
  return evaluateCharmmImproperTerms<Tcoord, Tcalc>(vk, CoordinateSeriesReader(csw), ecard,
                                                    system_index, force_scale_bits);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
Tcalc evalCmap(const Tcalc* cmap_patches, const int* cmap_patch_bounds, const int surf_idx,
               const int surf_dim, const int i_atom, const int j_atom, const int k_atom,
               const int l_atom, const int m_atom, const Tcoord* xcrd, const Tcoord* ycrd,
               const Tcoord* zcrd, const double* umat, const double* invu,
               const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
               const EvaluateForce eval_force, const Tcalc inv_gpos_factor,
               const Tcalc force_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc value_one = 1.0;

  // Compute displacements
  std::vector<Tcalc> acoef(16);
  Tcalc ab[3], bc[3], cd[3], de[3], crabbc[3], crbccd[3], crcdde[3], scr_phi[3], scr_psi[3];
  Tcalc phi_progression[4], psi_progression[4], acoef_psi[4];
  for (int i = 0; i < 4; i++) {
    acoef_psi[i] = 0.0;
  }
  if (isSignedIntegralScalarType<Tcoord>()) {
    ab[0] = static_cast<Tcalc>(xcrd[j_atom] - xcrd[i_atom]) * inv_gpos_factor;
    ab[1] = static_cast<Tcalc>(ycrd[j_atom] - ycrd[i_atom]) * inv_gpos_factor;
    ab[2] = static_cast<Tcalc>(zcrd[j_atom] - zcrd[i_atom]) * inv_gpos_factor;
    bc[0] = static_cast<Tcalc>(xcrd[k_atom] - xcrd[j_atom]) * inv_gpos_factor;
    bc[1] = static_cast<Tcalc>(ycrd[k_atom] - ycrd[j_atom]) * inv_gpos_factor;
    bc[2] = static_cast<Tcalc>(zcrd[k_atom] - zcrd[j_atom]) * inv_gpos_factor;
    cd[0] = static_cast<Tcalc>(xcrd[l_atom] - xcrd[k_atom]) * inv_gpos_factor;
    cd[1] = static_cast<Tcalc>(ycrd[l_atom] - ycrd[k_atom]) * inv_gpos_factor;
    cd[2] = static_cast<Tcalc>(zcrd[l_atom] - zcrd[k_atom]) * inv_gpos_factor;
    de[0] = static_cast<Tcalc>(xcrd[m_atom] - xcrd[l_atom]) * inv_gpos_factor;
    de[1] = static_cast<Tcalc>(ycrd[m_atom] - ycrd[l_atom]) * inv_gpos_factor;
    de[2] = static_cast<Tcalc>(zcrd[m_atom] - zcrd[l_atom]) * inv_gpos_factor;
  }
  else {
    ab[0] = xcrd[j_atom] - xcrd[i_atom];
    ab[1] = ycrd[j_atom] - ycrd[i_atom];
    ab[2] = zcrd[j_atom] - zcrd[i_atom];
    bc[0] = xcrd[k_atom] - xcrd[j_atom];
    bc[1] = ycrd[k_atom] - ycrd[j_atom];
    bc[2] = zcrd[k_atom] - zcrd[j_atom];
    cd[0] = xcrd[l_atom] - xcrd[k_atom];
    cd[1] = ycrd[l_atom] - ycrd[k_atom];
    cd[2] = zcrd[l_atom] - zcrd[k_atom];
    de[0] = xcrd[m_atom] - xcrd[l_atom];
    de[1] = ycrd[m_atom] - ycrd[l_atom];
    de[2] = zcrd[m_atom] - zcrd[l_atom];
  }
  imageCoordinates<Tcalc, Tcalc>(&ab[0], &ab[1], &ab[2], umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates<Tcalc, Tcalc>(&bc[0], &bc[1], &bc[2], umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates<Tcalc, Tcalc>(&cd[0], &cd[1], &cd[2], umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates<Tcalc, Tcalc>(&de[0], &de[1], &de[2], umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);
    
  // Compute the first dihedral
  crossProduct(ab, bc, crabbc);
  crossProduct(bc, cd, crbccd);
  Tcalc cos_phi = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
  if (tcalc_is_double) {
    cos_phi /= sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                    (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  }
  else {
    cos_phi /= sqrtf((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                     (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  }
  crossProduct(crabbc, crbccd, scr_phi);
  cos_phi = (cos_phi < -value_one) ? -value_one : (cos_phi > value_one) ? value_one : cos_phi;
  Tcalc phi = angleVerification(cos_phi, crabbc, crbccd, bc, scr_phi);
  if (tcalc_is_double) {
    phi += pi;
    phi = (phi < 0.0) ? phi + twopi : (phi >= twopi) ? phi - twopi : phi;
  }
  else {
    phi += pi_f;
    phi = (phi < 0.0f) ? phi + twopi_f : (phi >= twopi_f) ? phi - twopi_f : phi;
  }

  // Compute the second dihedral
  crossProduct(cd, de, crcdde);
  Tcalc cos_psi = crbccd[0]*crcdde[0] + crbccd[1]*crcdde[1] + crbccd[2]*crcdde[2];
  if (tcalc_is_double) {
    cos_psi /= sqrt((crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]) *
                    (crcdde[0]*crcdde[0] + crcdde[1]*crcdde[1] + crcdde[2]*crcdde[2]));
  }
  else {
    cos_psi /= sqrtf((crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]) *
                     (crcdde[0]*crcdde[0] + crcdde[1]*crcdde[1] + crcdde[2]*crcdde[2]));
  }
  crossProduct(crbccd, crcdde, scr_psi);
  cos_psi = (cos_psi < -value_one) ? -value_one : (cos_psi > value_one) ? value_one : cos_psi;
  Tcalc psi = angleVerification(cos_psi, crbccd, crcdde, cd, scr_psi);
  if (tcalc_is_double) {
    psi += pi;
    psi = (psi < 0.0) ? psi + twopi : (psi >= twopi) ? psi - twopi : psi;
  }
  else {
    psi += pi_f;
    psi = (psi < 0.0f) ? psi + twopi_f : (psi >= twopi_f) ? psi - twopi_f : psi;
  }
  
  // Compute the patch index (idx_phi, idx_psi) of the CMAP
  const Tcalc dsurf_dim = static_cast<Tcalc>(surf_dim);
  Tcalc phi_grid, psi_grid;
  if (tcalc_is_double) {
    phi_grid = phi * dsurf_dim * inverse_twopi;
    psi_grid = psi * dsurf_dim * inverse_twopi;
  }
  else {
    phi_grid = phi * dsurf_dim * inverse_twopi_f;
    psi_grid = psi * dsurf_dim * inverse_twopi_f;
  }
  const int idx_phi = phi_grid;
  const int idx_psi = psi_grid;
  const Tcalc phifrac = phi_grid - idx_phi;
  const Tcalc psifrac = psi_grid - idx_psi;

  // Draw in the matrix of spline values and derivatives
  const int patch_idx = cmap_patch_bounds[surf_idx] + (((idx_psi * surf_dim) + idx_phi) * 16);
  for (int i = 0; i < 16; i++) {
    acoef[i] = cmap_patches[patch_idx + i];
  }

  // Perform the matrix multiplications to obtain the bicubic spline coefficients
  phi_progression[0] = value_one;
  psi_progression[0] = value_one;
  for (int i = 1; i < 4; i++) {
    phi_progression[i] = phi_progression[i - 1] * phifrac;
    psi_progression[i] = psi_progression[i - 1] * psifrac;
  }
  matrixVectorMultiply<Tcalc>(acoef.data(), psi_progression, acoef_psi, 4, 4, 1.0, 1.0, 0.0);
  const Tcalc contrib = (phi_progression[0] * acoef_psi[0]) + (phi_progression[1] * acoef_psi[1]) +
                        (phi_progression[2] * acoef_psi[2]) + (phi_progression[3] * acoef_psi[3]);

  // Compute forces, if requested
  if (eval_force == EvaluateForce::YES) {

    // The derivatives along phi and psi follow from the energy expression, evaluation of the
    // guiding matrix equation for the bicubic spline potential:
    //
    //                                                    [ a00 a01 a02 a03     [ 1
    // Energy = [ 1 + phifrac + phifrac^2 + phifrac^3 ] *   a10 a11 a12 a13  *    psifrac
    //                                                      a20 a21 a22 a23       psifrac^2
    //                                                      a30 a31 a32 a33 ]     psifrac^3 ]
    Tcalc dphi, dpsi, mgab, mgbc, mgcd, mgde;
    if (tcalc_is_double) {
      dphi  = (((3.0 * acoef[15] * phifrac) + (2.0 * acoef[14])) * phifrac) + acoef[13];
      dphi *= psifrac;
      dphi +=       (((3.0 * acoef[11] * phifrac) + (2.0 * acoef[10])) * phifrac) + acoef[ 9];
      dphi *= psifrac;
      dphi +=       (((3.0 * acoef[ 7] * phifrac) + (2.0 * acoef[ 6])) * phifrac) + acoef[ 5];
      dphi *= psifrac;
      dphi +=       (((3.0 * acoef[ 3] * phifrac) + (2.0 * acoef[ 2])) * phifrac) + acoef[ 1];
      dpsi  = (((3.0 * acoef[15] * psifrac) + (2.0 * acoef[11])) * psifrac) + acoef[ 7];
      dpsi *= phifrac;
      dpsi +=       (((3.0 * acoef[14] * psifrac) + (2.0 * acoef[10])) * psifrac) + acoef[ 6];
      dpsi *= phifrac;
      dpsi +=       (((3.0 * acoef[13] * psifrac) + (2.0 * acoef[ 9])) * psifrac) + acoef[ 5];
      dpsi *= phifrac;
      dpsi +=       (((3.0 * acoef[12] * psifrac) + (2.0 * acoef[ 8])) * psifrac) + acoef[ 4];
      dphi *= dsurf_dim / twopi;
      dpsi *= dsurf_dim / twopi;
      mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
      mgde = sqrt(de[0]*de[0] + de[1]*de[1] + de[2]*de[2]);
    }
    else {
      dphi  = (((3.0f * acoef[15] * phifrac) + (2.0f * acoef[14])) * phifrac) + acoef[13];
      dphi *= psifrac;
      dphi +=       (((3.0f * acoef[11] * phifrac) + (2.0f * acoef[10])) * phifrac) + acoef[ 9];
      dphi *= psifrac;
      dphi +=       (((3.0f * acoef[ 7] * phifrac) + (2.0f * acoef[ 6])) * phifrac) + acoef[ 5];
      dphi *= psifrac;
      dphi +=       (((3.0f * acoef[ 3] * phifrac) + (2.0f * acoef[ 2])) * phifrac) + acoef[ 1];
      dpsi  = (((3.0f * acoef[15] * psifrac) + (2.0f * acoef[11])) * psifrac) + acoef[ 7];
      dpsi *= phifrac;
      dpsi +=       (((3.0f * acoef[14] * psifrac) + (2.0f * acoef[10])) * psifrac) + acoef[ 6];
      dpsi *= phifrac;
      dpsi +=       (((3.0f * acoef[13] * psifrac) + (2.0f * acoef[ 9])) * psifrac) + acoef[ 5];
      dpsi *= phifrac;
      dpsi +=       (((3.0f * acoef[12] * psifrac) + (2.0f * acoef[ 8])) * psifrac) + acoef[ 4];
      dphi *= dsurf_dim / twopi;
      dpsi *= dsurf_dim / twopi;
      mgab = sqrtf(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      mgbc = sqrtf(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      mgcd = sqrtf(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
      mgde = sqrtf(de[0]*de[0] + de[1]*de[1] + de[2]*de[2]);
    }

    // With the derivative in hand, evaluate the transformation of coordinates for either the
    // phi or psi dihedrals.
    const Tcalc invab = value_one / mgab;
    const Tcalc invbc = value_one / mgbc;
    const Tcalc invcd = value_one / mgcd;
    const Tcalc invde = value_one / mgde;
    const Tcalc invabc = invab * invbc;
    const Tcalc invbcd = invbc * invcd;
    const Tcalc invcde = invcd * invde;
    for (int i = 0; i < 3; i++) {
      crabbc[i] *= invabc;
      crbccd[i] *= invbcd;
      crcdde[i] *= invcde;
    }

    // Feed the gradient, negative of the derivative, into the functions below
    dphi *= -value_one;
    dpsi *= -value_one;

    // Phi accumulation: transform the rotational derivatives to cartesian coordinates
    const Tcalc phi_cosb = -(ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2]) * invab * invbc;
    const Tcalc phi_cosc = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2]) * invbc * invcd;
    Tcalc phi_isinb2, phi_isinc2;
    if (tcalc_is_double) {
      phi_isinb2 = (phi_cosb * phi_cosb < asymptotic_to_one_lf) ?
                   dphi / (1.0 - (phi_cosb * phi_cosb)) : dphi * inverse_one_minus_asymptote_lf;
      phi_isinc2 = (phi_cosc * phi_cosc < asymptotic_to_one_lf) ?
                   dphi / (1.0 - (phi_cosc * phi_cosc)) : dphi * inverse_one_minus_asymptote_lf;
    }
    else {
      phi_isinb2 = (phi_cosb * phi_cosb < asymptotic_to_one_f) ?
                   dphi / (value_one - (phi_cosb * phi_cosb)) :
                   dphi * inverse_one_minus_asymptote_f;
      phi_isinc2 = (phi_cosc * phi_cosc < asymptotic_to_one_f) ?
                   dphi / (value_one - (phi_cosc * phi_cosc)) :
                   dphi * inverse_one_minus_asymptote_f;      
    }
    const Tcalc phi_fa = -invab * phi_isinb2;
    const Tcalc phi_fb1 = (mgbc - (mgab * phi_cosb)) * invabc * phi_isinb2;
    const Tcalc phi_fb2 = phi_cosc * invbc * phi_isinc2;
    const Tcalc phi_fc1 = (mgbc - (mgcd * phi_cosc)) * invbcd * phi_isinc2;
    const Tcalc phi_fc2 = phi_cosb * invbc * phi_isinb2;
    const Tcalc phi_fd = -invcd * phi_isinc2;

    // Psi accumulation: transform the rotational derivatives to cartesian coordinates
    const Tcalc psi_cosb = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2]) * invbc * invcd;
    const Tcalc psi_cosc = -(cd[0]*de[0] + cd[1]*de[1] + cd[2]*de[2]) * invcd * invde;
    Tcalc psi_isinb2, psi_isinc2;
    if (tcalc_is_double) {
      psi_isinb2 = (psi_cosb * psi_cosb < asymptotic_to_one_lf) ?
                   dpsi / (1.0 - (psi_cosb * psi_cosb)) : dpsi * inverse_one_minus_asymptote_lf;
      psi_isinc2 = (psi_cosc * psi_cosc < asymptotic_to_one_lf) ?
                   dpsi / (1.0 - (psi_cosc * psi_cosc)) : dpsi * inverse_one_minus_asymptote_lf;
    }
    else {
      psi_isinb2 = (psi_cosb * psi_cosb < asymptotic_to_one_f) ?
                   dpsi / (value_one - (psi_cosb * psi_cosb)) :
                   dpsi * inverse_one_minus_asymptote_f;
      psi_isinc2 = (psi_cosc * psi_cosc < asymptotic_to_one_f) ?
                   dpsi / (value_one - (psi_cosc * psi_cosc)) :
                   dpsi * inverse_one_minus_asymptote_f;
    }
    const Tcalc psi_fa = -invbc * psi_isinb2;
    const Tcalc psi_fb1 = (mgcd - (mgbc * psi_cosb)) * invbcd * psi_isinb2;
    const Tcalc psi_fb2 = psi_cosc * invcd * psi_isinc2;
    const Tcalc psi_fc1 = (mgcd - (mgde * psi_cosc)) * invcde * psi_isinc2;
    const Tcalc psi_fc2 = psi_cosb * invcd * psi_isinb2;
    const Tcalc psi_fd = -invde * psi_isinc2;

    if (isSignedIntegralScalarType<Tforce>()) {

      // Accumulate the phi dihedral forces
      Tforce ifrc_ix = llround(crabbc[0] * phi_fa * force_factor);
      Tforce ifrc_jx = llround(((phi_fb1 * crabbc[0]) - (phi_fb2 * crbccd[0])) * force_factor);
      Tforce ifrc_lx = llround(-phi_fd * crbccd[0] * force_factor);
      xfrc[i_atom] += ifrc_ix;
      xfrc[j_atom] += ifrc_jx;
      xfrc[k_atom] -= ifrc_ix + ifrc_jx + ifrc_lx;
      xfrc[l_atom] += ifrc_lx;
      Tforce ifrc_iy = llround(crabbc[1] * phi_fa * force_factor);
      Tforce ifrc_jy = llround(((phi_fb1 * crabbc[1]) - (phi_fb2 * crbccd[1])) * force_factor);
      Tforce ifrc_ly = llround(-phi_fd * crbccd[1] * force_factor);
      yfrc[i_atom] += ifrc_iy;
      yfrc[j_atom] += ifrc_jy;
      yfrc[k_atom] -= ifrc_iy + ifrc_jy + ifrc_ly;
      yfrc[l_atom] += ifrc_ly;
      Tforce ifrc_iz = llround(crabbc[2] * phi_fa * force_factor);
      Tforce ifrc_jz = llround(((phi_fb1 * crabbc[2]) - (phi_fb2 * crbccd[2])) * force_factor);
      Tforce ifrc_lz = llround(-phi_fd * crbccd[2] * force_factor);
      zfrc[i_atom] += ifrc_iz;
      zfrc[j_atom] += ifrc_jz;
      zfrc[k_atom] -= ifrc_iz + ifrc_jz + ifrc_lz;
      zfrc[l_atom] += ifrc_lz;

      // Accumulate the psi dihedral forces
      ifrc_jx = llround(crbccd[0] * psi_fa * force_factor);
      Tforce ifrc_kx = llround(((psi_fb1 * crbccd[0]) - (psi_fb2 * crcdde[0])) * force_factor);
      Tforce ifrc_mx = llround(-psi_fd * crcdde[0] * force_factor);
      xfrc[j_atom] += ifrc_jx;
      xfrc[k_atom] += ifrc_kx;
      xfrc[l_atom] -= ifrc_jx + ifrc_kx + ifrc_mx;
      xfrc[m_atom] += ifrc_mx;
      ifrc_jy = llround(crbccd[1] * psi_fa * force_factor);
      Tforce ifrc_ky = llround(((psi_fb1 * crbccd[1]) - (psi_fb2 * crcdde[1])) * force_factor);
      Tforce ifrc_my = llround(-psi_fd * crcdde[1] * force_factor);
      yfrc[j_atom] += ifrc_jy;
      yfrc[k_atom] += ifrc_ky;
      yfrc[l_atom] -= ifrc_jy + ifrc_ky + ifrc_my;
      yfrc[m_atom] += ifrc_my;
      ifrc_jz = llround(crbccd[2] * psi_fa * force_factor);
      Tforce ifrc_kz = llround(((psi_fb1 * crbccd[2]) - (psi_fb2 * crcdde[2])) * force_factor);
      Tforce ifrc_mz = llround(-psi_fd * crcdde[2] * force_factor);
      zfrc[j_atom] += ifrc_jz;
      zfrc[k_atom] += ifrc_kz;
      zfrc[l_atom] -= ifrc_jz + ifrc_kz + ifrc_mz;
      zfrc[m_atom] += ifrc_mz;
    }
    else {
      
      // Accumulate the phi dihedral forces
      Tforce frc_ix = crabbc[0] * phi_fa;
      Tforce frc_jx = (phi_fb1 * crabbc[0]) - (phi_fb2 * crbccd[0]);
      Tforce frc_lx = -phi_fd * crbccd[0];
      xfrc[i_atom] += frc_ix;
      xfrc[j_atom] += frc_jx;
      xfrc[k_atom] -= frc_ix + frc_jx + frc_lx;
      xfrc[l_atom] += frc_lx;
      Tforce frc_iy = crabbc[1] * phi_fa;
      Tforce frc_jy = (phi_fb1 * crabbc[1]) - (phi_fb2 * crbccd[1]);
      Tforce frc_ly = -phi_fd * crbccd[1];
      yfrc[i_atom] += frc_iy;
      yfrc[j_atom] += frc_jy;
      yfrc[k_atom] -= frc_iy + frc_jy + frc_ly;
      yfrc[l_atom] += frc_ly;
      Tforce frc_iz = crabbc[2] * phi_fa;
      Tforce frc_jz = (phi_fb1 * crabbc[2]) - (phi_fb2 * crbccd[2]);
      Tforce frc_lz = -phi_fd * crbccd[2];
      zfrc[i_atom] += frc_iz;
      zfrc[j_atom] += frc_jz;
      zfrc[k_atom] -= frc_iz + frc_jz + frc_lz;
      zfrc[l_atom] += frc_lz;

      // Accumulate the psi dihedral forces
      frc_jx = crbccd[0] * psi_fa;
      Tforce frc_kx = (psi_fb1 * crbccd[0]) - (psi_fb2 * crcdde[0]);
      Tforce frc_mx = -psi_fd * crcdde[0];
      xfrc[j_atom] += frc_jx;
      xfrc[k_atom] += frc_kx;
      xfrc[l_atom] -= frc_jx + frc_kx + frc_mx;
      xfrc[m_atom] += frc_mx;
      frc_jy = crbccd[1] * psi_fa;
      Tforce frc_ky = (psi_fb1 * crbccd[1]) - (psi_fb2 * crcdde[1]);
      Tforce frc_my = -psi_fd * crcdde[1];
      yfrc[j_atom] += frc_jy;
      yfrc[k_atom] += frc_ky;
      yfrc[l_atom] -= frc_jy + frc_ky + frc_my;
      yfrc[m_atom] += frc_my;
      frc_jz = crbccd[2] * psi_fa;
      Tforce frc_kz = (psi_fb1 * crbccd[2]) - (psi_fb2 * crcdde[2]);
      Tforce frc_mz = -psi_fd * crcdde[2];
      zfrc[j_atom] += frc_jz;
      zfrc[k_atom] += frc_kz;
      zfrc[l_atom] -= frc_jz + frc_kz + frc_mz;
      zfrc[m_atom] += frc_mz;
    }
  }

  // Return the double-precision energy result.  It will be converted to a fixed-precision
  // quantity later.
  return contrib;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
double evaluateCmapTerms(const ValenceKit<Tcalc> vk, const Tcoord* xcrd, const Tcoord* ycrd,
                         const Tcoord* zcrd, const double* umat, const double* invu,
                         const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                         ScoreCard *ecard, const EvaluateForce eval_force,
                         const int system_index, const Tcalc inv_gpos_factor,
                         const Tcalc force_factor) {
  double cmap_energy = 0.0;
  llint cmap_acc = 0LL;
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();

  // Accumulate results by looping over all CMAP (here, meaning coupled dihedral) terms.
  for (int pos = 0; pos < vk.ncmap; pos++) {

    // Get parameters for an angle between atoms i, j, and k
    const int i_atom = vk.cmap_i_atoms[pos];
    const int j_atom = vk.cmap_j_atoms[pos];
    const int k_atom = vk.cmap_k_atoms[pos];
    const int l_atom = vk.cmap_l_atoms[pos];
    const int m_atom = vk.cmap_m_atoms[pos];
    const int surf_idx = vk.cmap_surf_idx[pos];
    const int surf_dim = vk.cmap_dim[surf_idx];
    const double contrib =
      evalCmap<Tcoord, Tforce, Tcalc>(vk.cmap_patches, vk.cmap_patch_bounds, vk.cmap_surf_idx[pos],
                                      vk.cmap_dim[surf_idx], i_atom, j_atom, k_atom, l_atom,
                                      m_atom, xcrd, ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc,
                                      zfrc, eval_force, inv_gpos_factor, force_factor);
    cmap_energy += contrib;
    cmap_acc += llround(contrib * nrg_scale_factor);
  }

  // Contribute results
  ecard->contribute(StateVariable::CMAP, cmap_acc, system_index);

  // Return the double-precision CMAP energy sum, if of interest
  return cmap_energy;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double evaluateCmapTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesReader<Tcoord> csr,
                         ScoreCard *ecard, const int system_index, const int force_scale_bits) {
  const size_t atom_os = static_cast<size_t>(system_index) *
                         roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t xfrm_os = static_cast<size_t>(system_index) * roundUp<size_t>(9, warp_size_zu);
  const Tcalc force_scale = (isSignedIntegralScalarType<Tcoord>()) ?
                            pow(2.0, force_scale_bits) : 1.0;
  return evaluateCmapTerms<Tcoord, Tcoord, Tcalc>(vk, &csr.xcrd[atom_os], &csr.ycrd[atom_os],
                                                  &csr.zcrd[atom_os], &csr.umat[xfrm_os],
                                                  &csr.invu[xfrm_os], csr.unit_cell, nullptr,
                                                  nullptr, nullptr, ecard, EvaluateForce::NO,
                                                  system_index, csr.inv_gpos_scale, force_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double evaluateCmapTerms(const ValenceKit<Tcalc> vk, const CoordinateSeriesWriter<Tcoord> csw,
                         ScoreCard *ecard, const int system_index, const int force_scale_bits) {
  return evaluateCmapTerms<Tcoord, Tcalc>(vk, CoordinateSeriesReader(csw), ecard, system_index,
                                          force_scale_bits);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
Vec2<Tcalc> evaluateAttenuated14Pair(const int i_atom, const int l_atom, const int attn_idx,
                                     const Tcalc coulomb_constant, const Tcalc* charges,
                                     const int* lj_param_idx, const Tcalc* attn14_elec_factors,
                                     const Tcalc* attn14_vdw_factors, const Tcalc* lja_14_coeff,
                                     const Tcalc* ljb_14_coeff, const Tcalc* lj_14_sigma,
                                     const int ljtab_offset, const int n_lj_types,
                                     const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                                     const double* umat, const double* invu,
                                     const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                                     Tforce* zfrc, const EvaluateForce eval_elec_force,
                                     const EvaluateForce eval_vdw_force,
                                     const Tcalc inv_gpos_factor, const Tcalc force_factor,
                                     const Tcalc clash_distance, const Tcalc clash_ratio) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc value_zero  = 0.0;
  const Tcalc value_nil   = 1.0e-6;
  const Tcalc value_half  = 1.0e-6;
  const Tcalc value_one   = 1.0;
  const Tcalc value_two   = 2.0;
  Tcalc dx, dy, dz;
  if (isSignedIntegralScalarType<Tcoord>()) {
    dx = static_cast<Tcalc>(xcrd[l_atom] - xcrd[i_atom]) * inv_gpos_factor;
    dy = static_cast<Tcalc>(ycrd[l_atom] - ycrd[i_atom]) * inv_gpos_factor;
    dz = static_cast<Tcalc>(zcrd[l_atom] - zcrd[i_atom]) * inv_gpos_factor;
  }
  else {
    dx = xcrd[l_atom] - xcrd[i_atom];
    dy = ycrd[l_atom] - ycrd[i_atom];
    dz = zcrd[l_atom] - zcrd[i_atom];
  }
  imageCoordinates<Tcalc, Tcalc>(&dx, &dy, &dz, umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);
  const int ilj_t = lj_param_idx[i_atom];
  const int jlj_t = lj_param_idx[l_atom] + ljtab_offset;
  const Tcalc vdw_scale = attn14_vdw_factors[attn_idx];
  const size_t ij_ljidx = (ilj_t * n_lj_types) + jlj_t;
  const Tcalc lja = lja_14_coeff[ij_ljidx] / vdw_scale;
  const Tcalc ljb = ljb_14_coeff[ij_ljidx] / vdw_scale;
  const Tcalc ele_scale = attn14_elec_factors[attn_idx];
  const Tcalc qiqj = (coulomb_constant * charges[i_atom] * charges[l_atom]) / ele_scale;
  const bool mitigate_clashes = (clash_distance > value_nil || clash_ratio > value_nil);
  Tcalc fmag = value_zero;
  Tcalc ele_contrib = value_zero;
  Tcalc vdw_contrib = value_zero;
  if (mitigate_clashes) {
    const Tcalc r = (tcalc_is_double) ? sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                                        sqrtf((dx * dx) + (dy * dy) + (dz * dz));
    if (eval_elec_force == EvaluateForce::YES) {
      quadraticCoreElectrostatics<Tcalc>(r, clash_distance, qiqj, &ele_contrib, &fmag);
    }
    else {
      quadraticCoreElectrostatics<Tcalc>(r, clash_distance, qiqj, &ele_contrib, nullptr);
    }
    if (eval_vdw_force == EvaluateForce::YES) {
      quarticCoreLennardJones<Tcalc>(r, clash_ratio, lja, ljb, lj_14_sigma[ij_ljidx],
                                     &vdw_contrib, &fmag);
    }
    else {
      quarticCoreLennardJones<Tcalc>(r, clash_ratio, lja, ljb, lj_14_sigma[ij_ljidx],
                                     &vdw_contrib, nullptr);
    }
  }
  else {

    // Standard evaluation of the Lennard-Jones and Coulomb potentials
    const Tcalc invr2 = value_one / ((dx * dx) + (dy * dy) + (dz * dz));
    const Tcalc invr = (tcalc_is_double) ? sqrt(invr2) : sqrtf(invr2);
    const Tcalc invr4 = invr2 * invr2;
    ele_contrib = qiqj * invr;
    vdw_contrib = (lja * invr4 * invr4 * invr4) - (ljb * invr4 * invr2);
    if (eval_elec_force == EvaluateForce::YES || eval_vdw_force == EvaluateForce::YES) {
      fmag = (eval_elec_force == EvaluateForce::YES) ? -(qiqj * invr * invr2) : 0.0;
      if (eval_vdw_force == EvaluateForce::YES) {
        if (tcalc_is_double) {
          fmag += ((6.0 * ljb) - (12.0 * lja * invr4 * invr2)) * invr4 * invr4;
        }
        else {
          fmag += ((6.0f * ljb) - (12.0f * lja * invr4 * invr2)) * invr4 * invr4;
        }
      }
    }
  }

  // Accumulate forces, if their calculation was requested
  if (eval_elec_force == EvaluateForce::YES || eval_vdw_force == EvaluateForce::YES) {
    if (isSignedIntegralScalarType<Tforce>()) {
      const Tforce ifmag_dx = llround(fmag * dx * force_factor);
      const Tforce ifmag_dy = llround(fmag * dy * force_factor);
      const Tforce ifmag_dz = llround(fmag * dz * force_factor);
      xfrc[i_atom] += ifmag_dx;
      yfrc[i_atom] += ifmag_dy;
      zfrc[i_atom] += ifmag_dz;
      xfrc[l_atom] -= ifmag_dx;
      yfrc[l_atom] -= ifmag_dy;
      zfrc[l_atom] -= ifmag_dz;
    }
    else {
      const Tcalc fmag_dx = fmag * dx;
      const Tcalc fmag_dy = fmag * dy;
      const Tcalc fmag_dz = fmag * dz;
      xfrc[i_atom] += fmag_dx;
      yfrc[i_atom] += fmag_dy;
      zfrc[i_atom] += fmag_dz;
      xfrc[l_atom] -= fmag_dx;
      yfrc[l_atom] -= fmag_dy;
      zfrc[l_atom] -= fmag_dz;
    }
  }
  return Vec2<Tcalc>(ele_contrib, vdw_contrib);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
Vec2<Tcalc> evaluateAttenuated14Pair(const int i_atom, const int l_atom, const int attn_idx,
                                     const Tcalc coulomb_constant, const Tcalc* charges,
                                     const int* lj_param_idx, const Tcalc* attn14_elec_factors,
                                     const Tcalc* attn14_vdw_factors, const Tcalc* lja_14_coeff,
                                     const Tcalc* ljb_14_coeff, const Tcalc* lj_14_sigma,
                                     const int n_lj_types, const Tcoord* xcrd, const Tcoord* ycrd,
                                     const Tcoord* zcrd, const double* umat, const double* invu,
                                     const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                                     Tforce* zfrc, const EvaluateForce eval_elec_force,
                                     const EvaluateForce eval_vdw_force,
                                     const Tcalc inv_gpos_factor, const Tcalc force_factor,
                                     const Tcalc clash_distance, const Tcalc clash_ratio) {
  return evaluateAttenuated14Pair(i_atom, l_atom, attn_idx, coulomb_constant, charges,
                                  lj_param_idx, attn14_elec_factors, attn14_vdw_factors,
                                  lja_14_coeff, ljb_14_coeff, lj_14_sigma, 0, n_lj_types, xcrd,
                                  ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc,
                                  eval_elec_force, eval_vdw_force, inv_gpos_factor, force_factor,
                                  clash_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
double2 evaluateAttenuated14Terms(const ValenceKit<Tcalc> vk, const NonbondedKit<Tcalc> nbk,
                                  const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                                  const double* umat, const double* invu,
                                  const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc,
                                  Tforce* zfrc, ScoreCard *ecard,
                                  const EvaluateForce eval_elec_force,
                                  const EvaluateForce eval_vdw_force, const int system_index,
                                  const Tcalc inv_gpos_factor, const Tcalc force_factor,
                                  const Tcalc clash_distance, const Tcalc clash_ratio) {
  double vdw_energy = 0.0;
  llint vdw_acc = 0LL;
  double ele_energy = 0.0;
  llint ele_acc = 0LL;
  const Tcalc nrg_scale_factor = ecard->getEnergyScalingFactor<Tcalc>();

  // Accumulate results in both electrostatic and van-der Waals (here, Lennard-Jones) energies by
  // looping over all 1:4 interactions.  The two 1:4 "non-bonded" energy components will be
  // combined into a tuple to return the single-threaded double-precision results.  Not only are
  // there two energies to accumulate, there are two possible sources of 1:4 interactions.  Start
  // with dihedrals that control 1:4 pairs between their I and L atoms.
  for (int pos = 0; pos < vk.ndihe; pos++) {

    // The zero index points to a pair interaction with zero electrostatic and zero Lennard-Jones
    // strength (full attenuation, a complete exclusion).  This is most often the case in dihedrals
    // that are part of a larger cosine series affecting the same atoms (only one of the dihedrals
    // will carry the responsibility of computing the I and L atom 1:4 interaction).  Another case
    // in which it can happen is dihedrals on either side of a six-membered ring.
    const int attn_idx = vk.dihe14_param_idx[pos];
    if (attn_idx == 0) {
      continue;
    }
    const Vec2<Tcalc> uc =
      evaluateAttenuated14Pair<Tcoord, Tforce, Tcalc>(vk.dihe_i_atoms[pos], vk.dihe_l_atoms[pos],
                                                      attn_idx, nbk.coulomb_constant, nbk.charge,
                                                      nbk.lj_idx, vk.attn14_elec, vk.attn14_vdw,
                                                      nbk.lja_14_coeff, nbk.ljb_14_coeff,
                                                      nbk.lj_14_sigma, nbk.n_lj_types, xcrd, ycrd,
                                                      zcrd, umat, invu, unit_cell, xfrc, yfrc,
                                                      zfrc, eval_elec_force, eval_vdw_force,
                                                      inv_gpos_factor, force_factor,
                                                      clash_distance, clash_ratio);
    ele_energy += uc.x;
    vdw_energy += uc.y;
    ele_acc += llround(uc.x * nrg_scale_factor);
    vdw_acc += llround(uc.y * nrg_scale_factor);
  }

  // Evaluate additional, inferred 1:4 attenuated interactions.  These occur between virtual sites
  // V and other atoms or virtual sites that are 1:4 to the parent atoms of V.
  for (int pos = 0; pos < vk.ninfr14; pos++) {
    const int attn_idx = vk.infr14_param_idx[pos];
    if (attn_idx == 0) {
      continue;
    }
    const Vec2<double> uc =
      evaluateAttenuated14Pair<Tcoord,
                               Tforce, Tcalc>(vk.infr14_i_atoms[pos], vk.infr14_l_atoms[pos],
                                              attn_idx, nbk.coulomb_constant, nbk.charge,
                                              nbk.lj_idx, vk.attn14_elec, vk.attn14_vdw,
                                              nbk.lja_14_coeff, nbk.ljb_14_coeff, nbk.lj_14_sigma,
                                              nbk.n_lj_types, xcrd, ycrd, zcrd, umat, invu,
                                              unit_cell, xfrc, yfrc, zfrc, eval_elec_force,
                                              eval_vdw_force, inv_gpos_factor, force_factor,
                                              clash_distance, clash_ratio);
    ele_energy += uc.x;
    vdw_energy += uc.y;
    ele_acc += llround(uc.x * nrg_scale_factor);
    vdw_acc += llround(uc.y * nrg_scale_factor);
  }

  // Contribute results
  ecard->contribute(StateVariable::ELEC_ONE_FOUR, ele_acc, system_index);
  ecard->contribute(StateVariable::VDW_ONE_FOUR, vdw_acc, system_index);

  // Return the double-precision energy sums, if of interest
  return { ele_energy, vdw_energy };
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double2 evaluateAttenuated14Terms(const ValenceKit<Tcalc> vk,
                                  const CoordinateSeriesReader<Tcoord> csr, ScoreCard *ecard,
                                  const int system_index, const int force_scale_bits,
                                  const Tcalc clash_distance, const Tcalc clash_ratio) {
  const size_t atom_os = static_cast<size_t>(system_index) *
                         roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t xfrm_os = static_cast<size_t>(system_index) * roundUp<size_t>(9, warp_size_zu);
  const Tcalc force_scale = (isSignedIntegralScalarType<Tcoord>()) ?
                            pow(2.0, force_scale_bits) : 1.0;
  return evaluateAttenuated14Terms<Tcoord,
                                   Tcoord, Tcalc>(vk, &csr.xcrd[atom_os], &csr.ycrd[atom_os],
                                                  &csr.zcrd[atom_os], &csr.umat[xfrm_os],
                                                  &csr.invu[xfrm_os], csr.unit_cell, nullptr,
                                                  nullptr, nullptr, ecard, EvaluateForce::NO,
                                                  system_index, csr.inv_gpos_scale, force_scale,
                                                  clash_distance, clash_ratio);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
double2 evaluateAttenuated14Terms(const ValenceKit<Tcalc> vk,
                                  const CoordinateSeriesWriter<Tcoord> csw, ScoreCard *ecard,
                                  const int system_index, const int force_scale_bits,
                                  const Tcalc clash_distance, const Tcalc clash_ratio) {
  return evaluateAttenuated14Terms<Tcoord, Tcalc>(vk, CoordinateSeriesReader(csw), ecard,
                                                  system_index, force_scale_bits, clash_distance,
                                                  clash_ratio);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
Tcalc evalPosnRestraint(const int p_atom, const int step_number, const int init_step,
                        const int finl_step, const Tcalc2 init_xy, const Tcalc2 finl_xy,
                        const Tcalc init_z, const Tcalc finl_z, const Tcalc2 init_keq,
                        const Tcalc2 finl_keq, const Tcalc4 init_r, const Tcalc4 finl_r,
                        const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                        const double* umat, const double* invu, const UnitCellType unit_cell,
                        Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, const EvaluateForce eval_force,
                        const Tcalc inv_gpos_factor, const Tcalc force_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Vec2<Tcalc> mixwt = computeRestraintMixture<Tcalc>(step_number, init_step, finl_step);
  Tcalc dx, dy, dz;
  if (isSignedIntegralScalarType<Tcoord>()) {
    dx = (static_cast<Tcalc>(xcrd[p_atom]) * inv_gpos_factor) -
         ((mixwt.x * init_xy.x) + (mixwt.y * finl_xy.x));
    dy = (static_cast<Tcalc>(ycrd[p_atom]) * inv_gpos_factor) -
         ((mixwt.x * init_xy.y) + (mixwt.y * finl_xy.y));
    dz = (static_cast<Tcalc>(zcrd[p_atom]) * inv_gpos_factor) -
         ((mixwt.x * init_z) + (mixwt.y * finl_z));
  }
  else {
    dx = xcrd[p_atom] - ((mixwt.x * init_xy.x) + (mixwt.y * finl_xy.x));
    dy = ycrd[p_atom] - ((mixwt.x * init_xy.y) + (mixwt.y * finl_xy.y));
    dz = zcrd[p_atom] - ((mixwt.x * init_z) + (mixwt.y * finl_z));
  }
  imageCoordinates<Tcalc, Tcalc>(&dx, &dy, &dz, umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);
  const Tcalc dr = (tcalc_is_double) ? sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                                       sqrtf((dx * dx) + (dy * dy) + (dz * dz));
  const Vec3<Tcalc> rst_eval = restraintDelta(Vec2<Tcalc>(init_keq.x, init_keq.y),
                                              Vec2<Tcalc>(finl_keq.x, finl_keq.y),
                                              Vec4<Tcalc>(init_r.x, init_r.y, init_r.z, init_r.w),
                                              Vec4<Tcalc>(finl_r.x, finl_r.y, finl_r.z, finl_r.w),
                                              mixwt, dr);
  
  // Compute forces
  if (eval_force == EvaluateForce::YES) {
    if (dr < constants::tiny) {

      // The case of positional restraints has a wrinkle when particles are already at their
      // exact target locations.  The force will likely be zero anyway, but it's possible to
      // define a positional restraint that forces a particle to be some finite distance away
      // from the target point, which would imply a non-zero force when the particle is at
      // the target location.  The proper direction of that force is an arbitrary thing in
      // such a case, so subdivide it among all three dimensions.
      const Tcalc fmag = (tcalc_is_double) ? 2.0  * rst_eval.x * rst_eval.y / sqrt(3.0) :
                                             2.0f * rst_eval.x * rst_eval.y / sqrtf(3.0f);
      if (isSignedIntegralScalarType<Tforce>()) {
        const Tforce ifmag = llround(fmag * force_factor);
        xfrc[p_atom] -= ifmag;
        yfrc[p_atom] -= ifmag;
        zfrc[p_atom] -= ifmag;        
      }
      else {
        xfrc[p_atom] -= fmag;
        yfrc[p_atom] -= fmag;
        zfrc[p_atom] -= fmag;
      }
    }
    else {
      const Tcalc fmag = (tcalc_is_double) ? 2.0  * rst_eval.x * rst_eval.y / dr :
                                             2.0f * rst_eval.x * rst_eval.y / dr;
      if (isSignedIntegralScalarType<Tforce>()) {
        xfrc[p_atom] -= llround(fmag * dx * force_factor);
        yfrc[p_atom] -= llround(fmag * dy * force_factor);
        zfrc[p_atom] -= llround(fmag * dz * force_factor);
      }
      else {
        xfrc[p_atom] -= fmag * dx;
        yfrc[p_atom] -= fmag * dy;
        zfrc[p_atom] -= fmag * dz;
      }
    }
  }
  return rst_eval.z;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
Tcalc evalBondRestraint(const int i_atom, const int j_atom, const int step_number,
                        const int init_step, const int finl_step, const Tcalc2 init_keq,
                        const Tcalc2 finl_keq, const Tcalc4 init_r, const Tcalc4 finl_r,
                        const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                        const double* umat, const double* invu, const UnitCellType unit_cell,
                        Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, const EvaluateForce eval_force,
                        const Tcalc inv_gpos_factor, const Tcalc force_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Vec2<Tcalc> mixwt = computeRestraintMixture<Tcalc>(step_number, init_step, finl_step);
  Tcalc dx, dy, dz;
  if (isSignedIntegralScalarType<Tcoord>()) {
    dx = static_cast<Tcalc>(xcrd[j_atom] - xcrd[i_atom]) * inv_gpos_factor;
    dy = static_cast<Tcalc>(ycrd[j_atom] - ycrd[i_atom]) * inv_gpos_factor;
    dz = static_cast<Tcalc>(zcrd[j_atom] - zcrd[i_atom]) * inv_gpos_factor;
  }
  else {
    dx = xcrd[j_atom] - xcrd[i_atom];
    dy = ycrd[j_atom] - ycrd[i_atom];
    dz = zcrd[j_atom] - zcrd[i_atom];
  }
  imageCoordinates<Tcalc, Tcalc>(&dx, &dy, &dz, umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);
  const Tcalc dr = (tcalc_is_double) ? sqrt((dx * dx) + (dy * dy) + (dz * dz)) :
                                       sqrtf((dx * dx) + (dy * dy) + (dz * dz));
  const Vec3<Tcalc> rst_eval = restraintDelta(Vec2<Tcalc>(init_keq.x, init_keq.y),
                                              Vec2<Tcalc>(finl_keq.x, finl_keq.y),
                                              Vec4<Tcalc>(init_r.x, init_r.y, init_r.z, init_r.w),
                                              Vec4<Tcalc>(finl_r.x, finl_r.y, finl_r.z, finl_r.w),
                                              mixwt, dr);
  if (eval_force == EvaluateForce::YES) {
    const Tcalc fmag = (tcalc_is_double) ? 2.0  * rst_eval.x * rst_eval.y / dr :
                                           2.0f * rst_eval.x * rst_eval.y / dr;
    if (isSignedIntegralScalarType<Tforce>()) {
      const Tforce ifmag_dx = llround(fmag * dx * force_factor);
      const Tforce ifmag_dy = llround(fmag * dy * force_factor);
      const Tforce ifmag_dz = llround(fmag * dz * force_factor);
      xfrc[i_atom] += ifmag_dx;
      yfrc[i_atom] += ifmag_dy;
      zfrc[i_atom] += ifmag_dz;
      xfrc[j_atom] -= ifmag_dx;
      yfrc[j_atom] -= ifmag_dy;
      zfrc[j_atom] -= ifmag_dz;
    }
    else {
      const Tcalc fmag_dx = fmag * dx;
      const Tcalc fmag_dy = fmag * dy;
      const Tcalc fmag_dz = fmag * dz;
      xfrc[i_atom] += fmag_dx;
      yfrc[i_atom] += fmag_dy;
      zfrc[i_atom] += fmag_dz;
      xfrc[j_atom] -= fmag_dx;
      yfrc[j_atom] -= fmag_dy;
      zfrc[j_atom] -= fmag_dz;
    }
  }
  return rst_eval.z;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
Tcalc evalAnglRestraint(const int i_atom, const int j_atom, const int k_atom,
                        const int step_number, const int init_step, const int finl_step,
                        const Tcalc2 init_keq, const Tcalc2 finl_keq, const Tcalc4 init_r,
                        const Tcalc4 finl_r, const Tcoord* xcrd, const Tcoord* ycrd,
                        const Tcoord* zcrd, const double* umat, const double* invu,
                        const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                        const EvaluateForce eval_force, const Tcalc inv_gpos_factor,
                        const Tcalc force_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc value_one = 1.0;
  Tcalc ba[3], bc[3];
  if (isSignedIntegralScalarType<Tcoord>()) {
    ba[0] = static_cast<Tcalc>(xcrd[i_atom] - xcrd[j_atom]) *inv_gpos_factor;
    ba[1] = static_cast<Tcalc>(ycrd[i_atom] - ycrd[j_atom]) *inv_gpos_factor;
    ba[2] = static_cast<Tcalc>(zcrd[i_atom] - zcrd[j_atom]) *inv_gpos_factor;
    bc[0] = static_cast<Tcalc>(xcrd[k_atom] - xcrd[j_atom]) *inv_gpos_factor;
    bc[1] = static_cast<Tcalc>(ycrd[k_atom] - ycrd[j_atom]) *inv_gpos_factor;
    bc[2] = static_cast<Tcalc>(zcrd[k_atom] - zcrd[j_atom]) *inv_gpos_factor;
  }
  else {
    ba[0] = xcrd[i_atom] - xcrd[j_atom];
    ba[1] = ycrd[i_atom] - ycrd[j_atom];
    ba[2] = zcrd[i_atom] - zcrd[j_atom];
    bc[0] = xcrd[k_atom] - xcrd[j_atom];
    bc[1] = ycrd[k_atom] - ycrd[j_atom];
    bc[2] = zcrd[k_atom] - zcrd[j_atom];
  }
  imageCoordinates<Tcalc, Tcalc>(&ba[0], &ba[1], &ba[2], umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates<Tcalc, Tcalc>(&bc[0], &bc[1], &bc[2], umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);

  // On to the angle force computation
  const Tcalc mgba = (ba[0] * ba[0]) + (ba[1] * ba[1]) + (ba[2] * ba[2]);
  const Tcalc mgbc = (bc[0] * bc[0]) + (bc[1] * bc[1]) + (bc[2] * bc[2]);
  const Tcalc invbabc = (tcalc_is_double) ? value_one / sqrt(mgba * mgbc) :
                                            value_one / sqrtf(mgba * mgbc);
  Tcalc costheta = ((ba[0] * bc[0]) + (ba[1] * bc[1]) + (ba[2] * bc[2])) * invbabc;
  costheta = (costheta < -value_one) ? -value_one : (costheta > value_one) ? value_one : costheta;
  const Tcalc theta = (tcalc_is_double) ? acos(costheta) : acosf(costheta);
  const Vec2<Tcalc> mixwt = computeRestraintMixture<Tcalc>(step_number, init_step, finl_step);
  const Vec3<Tcalc> rst_eval = restraintDelta(Vec2<Tcalc>(init_keq.x, init_keq.y),
                                              Vec2<Tcalc>(finl_keq.x, finl_keq.y),
                                              Vec4<Tcalc>(init_r.x, init_r.y, init_r.z, init_r.w),
                                              Vec4<Tcalc>(finl_r.x, finl_r.y, finl_r.z, finl_r.w),
                                              mixwt, theta);

  // Compute forces
  if (eval_force == EvaluateForce::YES) {
    const Tcalc dA = (tcalc_is_double) ?
                     -2.0  * rst_eval.x * rst_eval.y / sqrt(1.0 - (costheta * costheta)) :
                     -2.0f * rst_eval.x * rst_eval.y / sqrtf(value_one - (costheta * costheta));
    const Tcalc sqba = dA / mgba;
    const Tcalc sqbc = dA / mgbc;
    const Tcalc mbabc = dA * invbabc;
    if (isSignedIntegralScalarType<Tforce>()) {
      Tforce iadf[3], icdf[3];
      for (int i = 0; i < 3; i++) {
        iadf[i] = llround(((bc[i] * mbabc) - (costheta * ba[i] * sqba)) * force_factor);
        icdf[i] = llround(((ba[i] * mbabc) - (costheta * bc[i] * sqbc)) * force_factor);
      }
      xfrc[i_atom] -= iadf[0];
      yfrc[i_atom] -= iadf[1];
      zfrc[i_atom] -= iadf[2];
      xfrc[j_atom] += iadf[0] + icdf[0];
      yfrc[j_atom] += iadf[1] + icdf[1];
      zfrc[j_atom] += iadf[2] + icdf[2];
      xfrc[k_atom] -= icdf[0];
      yfrc[k_atom] -= icdf[1];
      zfrc[k_atom] -= icdf[2];
    }
    else {
      Tcalc adf[3], cdf[3];
      for (int i = 0; i < 3; i++) {
        adf[i] = (bc[i] * mbabc) - (costheta * ba[i] * sqba);
        cdf[i] = (ba[i] * mbabc) - (costheta * bc[i] * sqbc);
      }
      xfrc[i_atom] -= adf[0];
      yfrc[i_atom] -= adf[1];
      zfrc[i_atom] -= adf[2];
      xfrc[j_atom] += adf[0] + cdf[0];
      yfrc[j_atom] += adf[1] + cdf[1];
      zfrc[j_atom] += adf[2] + cdf[2];
      xfrc[k_atom] -= cdf[0];
      yfrc[k_atom] -= cdf[1];
      zfrc[k_atom] -= cdf[2];
    }
  }
  return rst_eval.z;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
Tcalc evalDiheRestraint(const int i_atom, const int j_atom, const int k_atom, const int l_atom,
                        const int step_number, const int init_step, const int finl_step,
                        const Tcalc2 init_keq, const Tcalc2 finl_keq, const Tcalc4 init_r,
                        const Tcalc4 finl_r, const Tcoord* xcrd, const Tcoord* ycrd,
                        const Tcoord* zcrd, const double* umat, const double* invu,
                        const UnitCellType unit_cell, Tforce* xfrc, Tforce* yfrc, Tforce* zfrc,
                        const EvaluateForce eval_force, const Tcalc inv_gpos_factor,
                        const Tcalc force_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const Tcalc value_one = 1.0;
  Tcalc ab[3], bc[3], cd[3], crabbc[3], crbccd[3], scr[3];
  if (isSignedIntegralScalarType<Tcoord>()) {
    ab[0] = static_cast<Tcalc>(xcrd[j_atom] - xcrd[i_atom]) * inv_gpos_factor;
    ab[1] = static_cast<Tcalc>(ycrd[j_atom] - ycrd[i_atom]) * inv_gpos_factor;
    ab[2] = static_cast<Tcalc>(zcrd[j_atom] - zcrd[i_atom]) * inv_gpos_factor;
    bc[0] = static_cast<Tcalc>(xcrd[k_atom] - xcrd[j_atom]) * inv_gpos_factor;
    bc[1] = static_cast<Tcalc>(ycrd[k_atom] - ycrd[j_atom]) * inv_gpos_factor;
    bc[2] = static_cast<Tcalc>(zcrd[k_atom] - zcrd[j_atom]) * inv_gpos_factor;
    cd[0] = static_cast<Tcalc>(xcrd[l_atom] - xcrd[k_atom]) * inv_gpos_factor;
    cd[1] = static_cast<Tcalc>(ycrd[l_atom] - ycrd[k_atom]) * inv_gpos_factor;
    cd[2] = static_cast<Tcalc>(zcrd[l_atom] - zcrd[k_atom]) * inv_gpos_factor;
  }
  else {
    ab[0] = xcrd[j_atom] - xcrd[i_atom];
    ab[1] = ycrd[j_atom] - ycrd[i_atom];
    ab[2] = zcrd[j_atom] - zcrd[i_atom];
    bc[0] = xcrd[k_atom] - xcrd[j_atom];
    bc[1] = ycrd[k_atom] - ycrd[j_atom];
    bc[2] = zcrd[k_atom] - zcrd[j_atom];
    cd[0] = xcrd[l_atom] - xcrd[k_atom];
    cd[1] = ycrd[l_atom] - ycrd[k_atom];
    cd[2] = zcrd[l_atom] - zcrd[k_atom];
  }
  imageCoordinates<Tcalc, Tcalc>(&ab[0], &ab[1], &ab[2], umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates<Tcalc, Tcalc>(&bc[0], &bc[1], &bc[2], umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);
  imageCoordinates<Tcalc, Tcalc>(&cd[0], &cd[1], &cd[2], umat, invu, unit_cell,
                                 ImagingMethod::MINIMUM_IMAGE);

  // Compute cross products and then the angle between the planes
  crossProduct(ab, bc, crabbc);
  crossProduct(bc, cd, crbccd);
  Tcalc costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
  if (tcalc_is_double) {
    costheta /= sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                     (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  }
  else {
    costheta /= sqrtf((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                      (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
  }
  crossProduct(crabbc, crbccd, scr);
  costheta = (costheta < -value_one) ? -value_one : (costheta > value_one) ? value_one : costheta;
  Tcalc theta;
  if (tcalc_is_double) {
    theta = (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0) ? acos(costheta) : -acos(costheta);
  }
  else {
    theta = (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0f) ?  acosf(costheta) :
                                                                  -acosf(costheta);
  }
  const Vec2<Tcalc> mixwt = computeRestraintMixture<Tcalc>(step_number, init_step, finl_step);

  // As part of the setup, the restraint has been arranged so that r1, r2, r3, and r4 are
  // monotonically increasing and span at most two pi radians.  The center of this arrangement
  // may not be at zero, but will be within the range [-pi, pi).  Image the angle to align with
  // the center of the restraint displacements r2 and r3.
  Tcalc midpoint, midpoint_delta;
  if (tcalc_is_double) {
    midpoint = 0.5 * (mixwt.x * (init_r.y + init_r.z) + mixwt.y * (finl_r.y + finl_r.z));
    midpoint_delta = imageValue(theta - midpoint, twopi, ImagingMethod::MINIMUM_IMAGE);
  }
  else {
    midpoint = 0.5f * (mixwt.x * (init_r.y + init_r.z) + mixwt.y * (finl_r.y + finl_r.z));
    midpoint_delta = imageValue(theta - midpoint, twopi_f, ImagingMethod::MINIMUM_IMAGE);
  }
  theta += midpoint_delta - (theta - midpoint);
  const Vec3<Tcalc> rst_eval = restraintDelta(Vec2<Tcalc>(init_keq.x, init_keq.y),
                                              Vec2<Tcalc>(finl_keq.x, finl_keq.y),
                                              Vec4<Tcalc>(init_r.x, init_r.y, init_r.z, init_r.w),
                                              Vec4<Tcalc>(finl_r.x, finl_r.y, finl_r.z, finl_r.w),
                                              mixwt, theta);
  
  // Compute forces
  if (eval_force == EvaluateForce::YES) {
    const Tcalc fr = -2.0 * rst_eval.x * rst_eval.y;
    Tcalc mgab, mgbc, mgcd;
    if (tcalc_is_double) {
      mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
    }
    else {
      mgab = sqrtf(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
      mgbc = sqrtf(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
      mgcd = sqrtf(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
    }
    const Tcalc invab = value_one / mgab;
    const Tcalc invbc = value_one / mgbc;
    const Tcalc invcd = value_one / mgcd;
    const Tcalc cosb = -(ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2]) * invab * invbc;
    const Tcalc cosc = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2]) * invbc * invcd;
    Tcalc isinb2, isinc2;
    if (tcalc_is_double) {
      isinb2 = (cosb * cosb < asymptotic_to_one_lf) ?
               fr / (1.0 - (cosb * cosb)) : fr * inverse_one_minus_asymptote_lf;
      isinc2 = (cosc * cosc < asymptotic_to_one_lf) ?
               fr / (1.0 - (cosc * cosc)) : fr * inverse_one_minus_asymptote_lf;
    }
    else {
      isinb2 = (cosb * cosb < asymptotic_to_one_f) ?
               fr / (value_one - (cosb * cosb)) : fr * inverse_one_minus_asymptote_f;
      isinc2 = (cosc * cosc < asymptotic_to_one_f) ?
               fr / (value_one - (cosc * cosc)) : fr * inverse_one_minus_asymptote_f;
    }
    const Tcalc invabc = invab * invbc;
    const Tcalc invbcd = invbc * invcd;
    for (int i = 0; i < 3; i++) {
      crabbc[i] *= invabc;
      crbccd[i] *= invbcd;
    }

    // Transform the rotational derivatives to Cartesian coordinates
    const Tcalc fa = -invab * isinb2;
    const Tcalc fb1 = (mgbc - (mgab * cosb)) * invabc * isinb2;
    const Tcalc fb2 = cosc * invbc * isinc2;
    const Tcalc fc1 = (mgbc - (mgcd * cosc)) * invbcd * isinc2;
    const Tcalc fc2 = cosb * invbc * isinb2;
    const Tcalc fd = -invcd * isinc2;
    if (isSignedIntegralScalarType<Tforce>()) {
      const Tforce ifrc_ix = llround(crabbc[0] * fa * force_factor);
      const Tforce ifrc_jx = llround(((fb1 * crabbc[0]) - (fb2 * crbccd[0])) * force_factor);
      const Tforce ifrc_lx = llround(-fd * crbccd[0] * force_factor);
      xfrc[i_atom] += ifrc_ix;
      xfrc[j_atom] += ifrc_jx;
      xfrc[k_atom] -= ifrc_ix + ifrc_jx + ifrc_lx;
      xfrc[l_atom] += ifrc_lx;
      const Tforce ifrc_iy = llround(crabbc[1] * fa * force_factor);
      const Tforce ifrc_jy = llround(((fb1 * crabbc[1]) - (fb2 * crbccd[1])) * force_factor);
      const Tforce ifrc_ly = llround(-fd * crbccd[1] * force_factor);
      yfrc[i_atom] += ifrc_iy;
      yfrc[j_atom] += ifrc_jy;
      yfrc[k_atom] -= ifrc_iy + ifrc_jy + ifrc_ly;
      yfrc[l_atom] += ifrc_ly;
      const Tforce ifrc_iz = llround(crabbc[2] * fa * force_factor);
      const Tforce ifrc_jz = llround(((fb1 * crabbc[2]) - (fb2 * crbccd[2])) * force_factor);
      const Tforce ifrc_lz = llround(-fd * crbccd[2] * force_factor);
      zfrc[i_atom] += ifrc_iz;
      zfrc[j_atom] += ifrc_jz;
      zfrc[k_atom] -= ifrc_iz + ifrc_jz + ifrc_lz;
      zfrc[l_atom] += ifrc_lz;
    }
    else {
      const Tforce frc_ix = crabbc[0] * fa;
      const Tforce frc_jx = (fb1 * crabbc[0]) - (fb2 * crbccd[0]);
      const Tforce frc_lx = -fd * crbccd[0];
      xfrc[i_atom] += frc_ix;
      xfrc[j_atom] += frc_jx;
      xfrc[k_atom] -= frc_ix + frc_jx + frc_lx;
      xfrc[l_atom] += frc_lx;
      const Tforce frc_iy = crabbc[1] * fa;
      const Tforce frc_jy = (fb1 * crabbc[1]) - (fb2 * crbccd[1]);
      const Tforce frc_ly = -fd * crbccd[1];
      yfrc[i_atom] += frc_iy;
      yfrc[j_atom] += frc_jy;
      yfrc[k_atom] -= frc_iy + frc_jy + frc_ly;
      yfrc[l_atom] += frc_ly;
      const Tforce frc_iz = crabbc[2] * fa;
      const Tforce frc_jz = (fb1 * crabbc[2]) - (fb2 * crbccd[2]);
      const Tforce frc_lz = -fd * crbccd[2];
      zfrc[i_atom] += frc_iz;
      zfrc[j_atom] += frc_jz;
      zfrc[k_atom] -= frc_iz + frc_jz + frc_lz;
      zfrc[l_atom] += frc_lz;
    }
  }
  return rst_eval.z;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
double evaluateRestraints(const RestraintKit<Tcalc, Tcalc2, Tcalc4> rar, const Tcoord* xcrd,
                          const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                          const double* invu, const UnitCellType unit_cell, Tforce* xfrc,
                          Tforce* yfrc, Tforce* zfrc, ScoreCard *ecard,
                          const EvaluateForce eval_force, const int system_index,
                          const int step_number, const Tcalc inv_gpos_factor,
                          const Tcalc force_factor) {
  double rest_energy = 0.0;
  llint rest_acc = 0LL;
  const Tcalc nrg_scale_factor = ecard->getEnergyScalingFactor<Tcalc>();

  // Accumulate results by looping over all restraint terms
  for (int i = 0; i < rar.nposn; i++) {
    const double contrib =
      evalPosnRestraint(rar.rposn_atoms[i], step_number, rar.rposn_init_step[i],
                        rar.rposn_finl_step[i], rar.rposn_init_xy[i], rar.rposn_finl_xy[i],
                        rar.rposn_init_z[i], rar.rposn_finl_z[i], rar.rposn_init_keq[i],
                        rar.rposn_finl_keq[i], rar.rposn_init_r[i], rar.rposn_finl_r[i], xcrd,
                        ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc, eval_force,
                        inv_gpos_factor, force_factor);
    rest_energy += contrib;
    rest_acc += llround(contrib * nrg_scale_factor);
  }
  for (int pos = 0; pos < rar.nbond; pos++) {
    const double contrib =
      evalBondRestraint<Tcoord, Tforce,
                        Tcalc, Tcalc2, Tcalc4>(rar.rbond_i_atoms[pos], rar.rbond_j_atoms[pos],
                                               step_number, rar.rbond_init_step[pos],
                                               rar.rbond_finl_step[pos], rar.rbond_init_keq[pos],
                                               rar.rbond_finl_keq[pos], rar.rbond_init_r[pos],
                                               rar.rbond_finl_r[pos], xcrd, ycrd, zcrd, umat, invu,
                                               unit_cell, xfrc, yfrc, zfrc, eval_force,
                                               inv_gpos_factor, force_factor);
    rest_energy += contrib;
    rest_acc += llround(contrib * nrg_scale_factor);
  }
  for (int pos = 0; pos < rar.nangl; pos++) {
    const double contrib =
      evalAnglRestraint<Tcoord, Tforce,
                        Tcalc, Tcalc2, Tcalc4>(rar.rangl_i_atoms[pos], rar.rangl_j_atoms[pos],
                                               rar.rangl_k_atoms[pos], step_number,
                                               rar.rangl_init_step[pos], rar.rangl_finl_step[pos],
                                               rar.rangl_init_keq[pos], rar.rangl_finl_keq[pos],
                                               rar.rangl_init_r[pos], rar.rangl_finl_r[pos], xcrd,
                                               ycrd, zcrd, umat, invu, unit_cell, xfrc, yfrc, zfrc,
                                               eval_force, inv_gpos_factor, force_factor);
    rest_energy += contrib;
    rest_acc += llround(contrib * nrg_scale_factor);
  }
  for (int pos = 0; pos < rar.ndihe; pos++) {
    const double contrib =
      evalDiheRestraint<Tcoord, Tforce,
                        Tcalc, Tcalc2, Tcalc4>(rar.rdihe_i_atoms[pos], rar.rdihe_j_atoms[pos],
                                               rar.rdihe_k_atoms[pos], rar.rdihe_l_atoms[pos],
                                               step_number, rar.rdihe_init_step[pos],
                                               rar.rdihe_finl_step[pos], rar.rdihe_init_keq[pos],
                                               rar.rdihe_finl_keq[pos], rar.rdihe_init_r[pos],
                                               rar.rdihe_finl_r[pos], xcrd, ycrd, zcrd, umat, invu,
                                               unit_cell, xfrc, yfrc, zfrc, eval_force,
                                               inv_gpos_factor, force_factor);
    rest_energy += contrib;
    rest_acc += llround(contrib * nrg_scale_factor);
  }

  // Contribute results
  ecard->contribute(StateVariable::RESTRAINT, rest_acc, system_index);

  // Return the double-precision energy sum, if of interest
  return rest_energy;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc, typename Tcalc2, typename Tcalc4>
double evaluateRestraints(const RestraintKit<Tcalc, Tcalc2, Tcalc4> rar,
                          const CoordinateSeriesReader<Tcoord> csr, ScoreCard *ecard,
                          const int system_index, const int step_number,
                          const int force_scale_bits) {
  const size_t atom_os = static_cast<size_t>(system_index) *
                         roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t xfrm_os = static_cast<size_t>(system_index) * roundUp<size_t>(9, warp_size_zu);
  const Tcalc force_scale = (isSignedIntegralScalarType<Tcoord>()) ?
                            pow(2.0, force_scale_bits) : 1.0;
  return evaluateRestraints<Tcoord, Tcoord,
                            Tcalc, Tcalc2, Tcalc4>(rar, csr.xcrd, csr.ycrd, csr.zcrd, csr.umat,
                                                   csr.invu, csr.unit_cell, csr.xfrc, csr.yfrc,
                                                   csr.zfrc, ecard, EvaluateForce::NO,
                                                   system_index, step_number);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc, typename Tcalc2, typename Tcalc4>
double evaluateRestraints(const RestraintKit<Tcalc, Tcalc2, Tcalc4> rar,
                          const CoordinateSeriesWriter<Tcoord> csw, ScoreCard *ecard,
                          const int system_index, const int step_number,
                          const int force_scale_bits) {
  return evaluateRestraints<Tcoord, Tcalc, Tcalc2, Tcalc4>(rar, CoordinateSeriesReader(csw), ecard,
                                                           system_index, step_number,
                                                           force_scale_bits);
}

} // namespace energy
} // namespace stormm
