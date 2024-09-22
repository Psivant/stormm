// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace mm {

//-------------------------------------------------------------------------------------------------
template <typename Tforce, typename Tcalc>
void computeGradientMove(Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, Tforce* xprv_move,
                         Tforce* yprv_move, Tforce* zprv_move, Tcalc* x_cg_temp, Tcalc* y_cg_temp,
                         Tcalc* z_cg_temp, const int natom, const int step, const int sd_steps,
                         const Tcalc force_factor) {
  const Tcalc value_one = 1.0;
  const Tcalc inv_force_factor = 1.0 / force_factor;
  const bool tforce_is_sgnint = isSignedIntegralScalarType<Tforce>();
  
  // Apply the conjugate gradient protocol
  if (step >= sd_steps) {

    // Both gg and dgg will be accumulated in double-precision on the GPU, as the order of
    // operations can be enforced.
    double gg = 0.0;
    double dgg = 0.0;
    if (tforce_is_sgnint) {
      for (int i = 0; i < natom; i++) {
        const Tcalc xprv_tmv = static_cast<Tcalc>(xprv_move[i]) * inv_force_factor;
        const Tcalc yprv_tmv = static_cast<Tcalc>(yprv_move[i]) * inv_force_factor;
        const Tcalc zprv_tmv = static_cast<Tcalc>(zprv_move[i]) * inv_force_factor;
        gg += xprv_tmv * xprv_tmv;
        gg += yprv_tmv * yprv_tmv;
        gg += zprv_tmv * zprv_tmv;
        const Tcalc xfrc_t  = static_cast<Tcalc>(xfrc[i]) * inv_force_factor;
        const Tcalc yfrc_t  = static_cast<Tcalc>(yfrc[i]) * inv_force_factor;
        const Tcalc zfrc_t  = static_cast<Tcalc>(zfrc[i]) * inv_force_factor;
        const Tcalc xdiff_t = static_cast<Tcalc>(xfrc[i] - xprv_move[i]) * inv_force_factor;
        const Tcalc ydiff_t = static_cast<Tcalc>(yfrc[i] - yprv_move[i]) * inv_force_factor;
        const Tcalc zdiff_t = static_cast<Tcalc>(zfrc[i] - zprv_move[i]) * inv_force_factor;
        dgg += xdiff_t * xfrc_t;
        dgg += ydiff_t * yfrc_t;
        dgg += zdiff_t * zfrc_t;
      }
    }
    else {
      for (int i = 0; i < natom; i++) {
        gg += xprv_move[i] * xprv_move[i];
        gg += yprv_move[i] * yprv_move[i];
        gg += zprv_move[i] * zprv_move[i];
        dgg += (xfrc[i] - xprv_move[i]) * xfrc[i];
        dgg += (yfrc[i] - yprv_move[i]) * yfrc[i];
        dgg += (zfrc[i] - zprv_move[i]) * zfrc[i];
      }
    }
    const double gam = (step == sd_steps) ? 0.0 : dgg / gg;
    if (tforce_is_sgnint) {
      for (int i = 0; i < natom; i++) {
        xprv_move[i] = xfrc[i];
        yprv_move[i] = yfrc[i];
        zprv_move[i] = zfrc[i];
        const Tcalc xprv_tmv = static_cast<Tcalc>(xprv_move[i]) * inv_force_factor;
        const Tcalc yprv_tmv = static_cast<Tcalc>(yprv_move[i]) * inv_force_factor;
        const Tcalc zprv_tmv = static_cast<Tcalc>(zprv_move[i]) * inv_force_factor;
        x_cg_temp[i] = xprv_tmv + (gam * x_cg_temp[i]);
        y_cg_temp[i] = yprv_tmv + (gam * y_cg_temp[i]);
        z_cg_temp[i] = zprv_tmv + (gam * z_cg_temp[i]);
        xfrc[i] = llround(x_cg_temp[i] * force_factor);
        yfrc[i] = llround(y_cg_temp[i] * force_factor);
        zfrc[i] = llround(z_cg_temp[i] * force_factor);
      }
    }
    else {
      for (int i = 0; i < natom; i++) {
        xprv_move[i] = xfrc[i];
        yprv_move[i] = yfrc[i];
        zprv_move[i] = zfrc[i];
        x_cg_temp[i] = xprv_move[i] + (gam * x_cg_temp[i]);
        y_cg_temp[i] = yprv_move[i] + (gam * y_cg_temp[i]);
        z_cg_temp[i] = zprv_move[i] + (gam * z_cg_temp[i]);
        xfrc[i] = x_cg_temp[i];
        yfrc[i] = y_cg_temp[i];
        zfrc[i] = z_cg_temp[i];
      }
    }
  }

  // Normalize the force vector to get the move.  This will be a double-precision accumulation on
  // the GPU as well, as the number is sensitive and the order of operations can be enforced.
  double msum = 0.0;
  if (tforce_is_sgnint) {
    for (int i = 0; i < natom; i++) {
      const Tcalc fx = static_cast<Tcalc>(xfrc[i]) * inv_force_factor;
      const Tcalc fy = static_cast<Tcalc>(yfrc[i]) * inv_force_factor;
      const Tcalc fz = static_cast<Tcalc>(zfrc[i]) * inv_force_factor;
      msum += fx * fx;
      msum += fy * fy;
      msum += fz * fz;
    }
  }
  else {
    for (int i = 0; i < natom; i++) {
      const Tcalc fx = xfrc[i];
      const Tcalc fy = yfrc[i];
      const Tcalc fz = zfrc[i];
      msum += fx * fx;
      msum += fy * fy;
      msum += fz * fz;
    }
  }
  msum = 1.0 / sqrt(msum);
  if (tforce_is_sgnint) {
    for (int i = 0; i < natom; i++) {
      xfrc[i] = llround(static_cast<double>(xfrc[i]) * msum);
      yfrc[i] = llround(static_cast<double>(yfrc[i]) * msum);
      zfrc[i] = llround(static_cast<double>(zfrc[i]) * msum);
    }
  }
  else {
    for (int i = 0; i < natom; i++) {
      xfrc[i] *= msum;
      yfrc[i] *= msum;
      zfrc[i] *= msum;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
void moveParticles(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const Tforce* xmove,
                   const Tforce* ymove, const Tforce* zmove, const double* umat,
                   const double* invu, const UnitCellType unit_cell,
                   const VirtualSiteKit<Tcalc> &vsk, const int natom, const Tcalc dist,
                   const Tcalc gpos_factor, const Tcalc force_factor) {
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const bool tforce_is_sgnint = isSignedIntegralScalarType<Tforce>();
  if (tcoord_is_sgnint) {
    if (tforce_is_sgnint) {
      const Tcalc gfdist = dist * gpos_factor / force_factor;
      for (int i = 0; i < natom; i++) {
        xcrd[i] += llround(static_cast<Tcalc>(xmove[i]) * gfdist);
        ycrd[i] += llround(static_cast<Tcalc>(ymove[i]) * gfdist);
        zcrd[i] += llround(static_cast<Tcalc>(zmove[i]) * gfdist);
      }
    }
    else {
      const Tcalc gdist = dist * gpos_factor;
      for (int i = 0; i < natom; i++) {
        xcrd[i] += llround(xmove[i] * gdist);
        ycrd[i] += llround(ymove[i] * gdist);
        zcrd[i] += llround(zmove[i] * gdist);
      }      
    }
  }
  else {
    if (tforce_is_sgnint) {
      const Tcalc fdist = dist / force_factor;
      for (int i = 0; i < natom; i++) {
        xcrd[i] += static_cast<Tcalc>(xmove[i]) * fdist;
        ycrd[i] += static_cast<Tcalc>(ymove[i]) * fdist;
        zcrd[i] += static_cast<Tcalc>(zmove[i]) * fdist;
      }
    }
    else {
      for (int i = 0; i < natom; i++) {
        xcrd[i] += xmove[i] * dist;
        ycrd[i] += ymove[i] * dist;
        zcrd[i] += zmove[i] * dist;
      }
    }
  }
  placeVirtualSites(xcrd, ycrd, zcrd, umat, invu, unit_cell, vsk);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc, typename Tcalc2, typename Tcalc4>
ScoreCard minimize(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, Tforce* xfrc, Tforce* yfrc,
                   Tforce* zfrc, Tforce* xprv_move, Tforce* yprv_move, Tforce* zprv_move,
                   Tcalc* x_cg_temp, Tcalc* y_cg_temp, Tcalc* z_cg_temp,
                   const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                   const ImplicitSolventKit<Tcalc> &isk, const NeckGeneralizedBornKit<Tcalc> &ngbk,
                   const RestraintKit<Tcalc, Tcalc2, Tcalc4> &rar,
                   const VirtualSiteKit<Tcalc> &vsk, const StaticExclusionMaskReader &ser,
                   const MinimizeControls &mincon, const int nrg_scale_bits,
                   const Tcalc gpos_factor, const Tcalc force_factor) {

  // Pre-allocate the Generalized Born buffers, if necessary
  std::vector<double> effective_gb_radii, psi, sumdeijda;
  switch (isk.igb) {
  case ImplicitSolventModel::NONE:
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    effective_gb_radii.resize(vk.natom);
    psi.resize(vk.natom);
    sumdeijda.resize(vk.natom);
    break;
  }
  
  // Loop for the requested number of cycles
  const int cd_steps = mincon.getClashDampingCycles();
  const Tcalc clash_distance = mincon.getAbsoluteClashDistance();
  const Tcalc clash_ratio = mincon.getVdwClashRatio();
  ScoreCard sc(1, mincon.getTotalCycles() + 1, nrg_scale_bits), sc_temp(1);
  const Tcalc inv_gpos_factor = static_cast<Tcalc>(1.0) / gpos_factor;
  Tcalc move_scale = mincon.getInitialStep() * gpos_factor;
  Tcalc evec[4], mvec[4], abcd_coefs[4], d_abcd[3], amat[16], inva[16];
  for (int i = 0; i < vk.natom; i++) {
    xprv_move[i] = 0.0;
    yprv_move[i] = 0.0;
    zprv_move[i] = 0.0;
    x_cg_temp[i] = 0.0;
    y_cg_temp[i] = 0.0;
    z_cg_temp[i] = 0.0;
  }
  for (int step = 0; step < mincon.getTotalCycles(); step++) {

    // Compute the forces on all particles
    for (int i = 0; i < vk.natom; i++) {
      xfrc[i] = 0.0;
      yfrc[i] = 0.0;
      zfrc[i] = 0.0;
    }
    if (step < cd_steps) {
      const Tcalc cd_progress = static_cast<Tcalc>(cd_steps - step) / static_cast<Tcalc>(cd_steps);
      const Tcalc tmp_clash_distance = clash_distance * cd_progress;
      const Tcalc tmp_clash_ratio    = clash_ratio * cd_progress;
      switch (isk.igb) {
      case ImplicitSolventModel::NONE:
        evalNonbValeRestMM<Tcoord, Tforce,
                           Tcalc, Tcalc2, Tcalc4>(xcrd, ycrd, zcrd, nullptr, nullptr,
                                                  UnitCellType::NONE, xfrc, yfrc, zfrc, &sc, vk,
                                                  nbk, ser, rar, EvaluateForce::YES, 0, step,
                                                  inv_gpos_factor, force_factor,
                                                  tmp_clash_distance, tmp_clash_ratio);
        break;
      case ImplicitSolventModel::HCT_GB:
      case ImplicitSolventModel::OBC_GB:
      case ImplicitSolventModel::OBC_GB_II:
      case ImplicitSolventModel::NECK_GB:
      case ImplicitSolventModel::NECK_GB_II:
        evalRestrainedMMGB<Tcoord, Tforce,
                           Tcalc, Tcalc2, Tcalc4>(xcrd, ycrd, zcrd, nullptr, nullptr,
                                                  UnitCellType::NONE, xfrc, yfrc, zfrc, &sc, vk,
                                                  nbk, ser, isk, ngbk, effective_gb_radii.data(),
                                                  psi.data(), sumdeijda.data(), rar,
                                                  EvaluateForce::YES, 0, step, inv_gpos_factor,
                                                  force_factor, tmp_clash_distance,
                                                  tmp_clash_ratio);
        break;
      }
    }
    else {
      switch (isk.igb) {
      case ImplicitSolventModel::NONE:
        evalNonbValeRestMM<Tcoord, Tforce,
                           Tcalc, Tcalc2, Tcalc4>(xcrd, ycrd, zcrd, nullptr, nullptr,
                                                  UnitCellType::NONE, xfrc, yfrc, zfrc, &sc, vk,
                                                  nbk, ser, rar, EvaluateForce::YES, 0, step,
                                                  inv_gpos_factor, force_factor);
        break;
      case ImplicitSolventModel::HCT_GB:
      case ImplicitSolventModel::OBC_GB:
      case ImplicitSolventModel::OBC_GB_II:
      case ImplicitSolventModel::NECK_GB:
      case ImplicitSolventModel::NECK_GB_II:
        evalRestrainedMMGB<Tcoord, Tforce,
                           Tcalc, Tcalc2, Tcalc4>(xcrd, ycrd, zcrd, nullptr, nullptr,
                                                  UnitCellType::NONE, xfrc, yfrc, zfrc, &sc, vk,
                                                  nbk, ser, isk, ngbk, effective_gb_radii.data(),
                                                  psi.data(), sumdeijda.data(), rar,
                                                  EvaluateForce::YES, 0, step, inv_gpos_factor,
                                                  force_factor);
        break;
      }
    }
    transmitVirtualSiteForces<Tcalc, Tcalc>(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, nullptr, nullptr,
                                            UnitCellType::NONE, vsk);
    if (step % mincon.getDiagnosticPrintFrequency() == 0) {
      sc.commit(StateVariable::ALL_STATES);
      sc.incrementSampleCount();
      sc.setLastTimeStep(step);
    }
    evec[0] = sc.reportTotalEnergy();
    mvec[0] = 0.0;

    // Generate the move
    computeGradientMove<Tforce, Tcalc>(xfrc, yfrc, zfrc, xprv_move, yprv_move, zprv_move,
                                       x_cg_temp, y_cg_temp, z_cg_temp, nbk.natom, step,
                                       mincon.getSteepestDescentCycles(), force_factor);
    
    // Implement the move three times, compute the energy, and arrive at a minimum value.
    double move_scale_factor = 1.0;
    for (int i = 0; i < 3; i++) {
      sc_temp.initialize();
      moveParticles<Tcoord, Tforce, Tcalc>(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, nullptr, nullptr,
                                           UnitCellType::NONE, vsk, nbk.natom,
                                           move_scale * move_scale_factor, force_factor);
      if (step < cd_steps) {
        const Tcalc cd_progress = static_cast<Tcalc>(cd_steps - step) /
                                  static_cast<Tcalc>(cd_steps);
        const Tcalc tmp_clash_distance = clash_distance * cd_progress;
        const Tcalc tmp_clash_ratio    = clash_ratio * cd_progress;
        switch (isk.igb) {
        case ImplicitSolventModel::NONE:
          evalNonbValeRestMM<Tcoord, Tforce,
                             Tcalc, Tcalc2, Tcalc4>(xcrd, ycrd, zcrd, nullptr, nullptr,
                                                    UnitCellType::NONE, xfrc, yfrc, zfrc, &sc_temp,
                                                    vk, nbk, ser, rar, EvaluateForce::NO, 0, step,
                                                    inv_gpos_factor, force_factor,
                                                    tmp_clash_distance, tmp_clash_ratio);
          break;
        case ImplicitSolventModel::HCT_GB:
        case ImplicitSolventModel::OBC_GB:
        case ImplicitSolventModel::OBC_GB_II:
        case ImplicitSolventModel::NECK_GB:
        case ImplicitSolventModel::NECK_GB_II:
          evalRestrainedMMGB<Tcoord, Tforce,
                             Tcalc, Tcalc2, Tcalc4>(xcrd, ycrd, zcrd, nullptr, nullptr,
                                                    UnitCellType::NONE, xfrc, yfrc, zfrc, &sc_temp,
                                                    vk, nbk, ser, isk, ngbk,
                                                    effective_gb_radii.data(), psi.data(),
                                                    sumdeijda.data(), rar, EvaluateForce::NO, 0,
                                                    step, inv_gpos_factor, force_factor,
                                                    tmp_clash_distance, tmp_clash_ratio);
          break;
        }
      }
      else {
        switch (isk.igb) {
        case ImplicitSolventModel::NONE:
          evalNonbValeRestMM<Tcoord, Tforce,
                             Tcalc, Tcalc2, Tcalc4>(xcrd, ycrd, zcrd, nullptr, nullptr,
                                                    UnitCellType::NONE, xfrc, yfrc, zfrc, &sc_temp,
                                                    vk, nbk, ser, rar, EvaluateForce::NO, 0, step,
                                                    inv_gpos_factor, force_factor);
          break;
        case ImplicitSolventModel::HCT_GB:
        case ImplicitSolventModel::OBC_GB:
        case ImplicitSolventModel::OBC_GB_II:
        case ImplicitSolventModel::NECK_GB:
        case ImplicitSolventModel::NECK_GB_II:
          evalRestrainedMMGB<Tcoord, Tforce,
                             Tcalc, Tcalc2, Tcalc4>(xcrd, ycrd, zcrd, nullptr, nullptr,
                                                    UnitCellType::NONE, xfrc, yfrc, zfrc, &sc_temp,
                                                    vk, nbk, ser, isk, ngbk,
                                                    effective_gb_radii.data(), psi.data(),
                                                    sumdeijda.data(), rar, EvaluateForce::NO, 0,
                                                    step, inv_gpos_factor, force_factor);
          break;
        }
      }
      evec[i + 1] = sc_temp.reportTotalEnergy();
      mvec[i + 1] = mvec[i] + move_scale_factor;
      if (evec[i + 1] < evec[0]) {
        const double idecay = 0.01 * static_cast<double>(i);
        move_scale_factor *= 1.05 - idecay;
      }
      else {
        const double decay_factor = 0.75 + (0.05 * static_cast<double>(i));
        move_scale_factor *= decay_factor;
      }
    }

    // Solve the linear system of equations to arrive at a polynomial Ax^3 + Bx^2 + Cx + D = evec
    // to solve the energy surface for x = 0.0 (gets evec[0]), mvec[1] (gets evec[1]), mvec[2],
    // and mvec[3].  The move vector mvec is kept in units of the initial move, to protect the
    // numerics of solving this system of linear equations as the step size becomes very small.
    for (int i = 0; i < 4; i++) {
      const double x = mvec[i];
      amat[i     ] = x * x * x;
      amat[i +  4] = x * x;
      amat[i +  8] = x;
      amat[i + 12] = 1.0;
      abcd_coefs[i] = 0.0;
    }
    invertSquareMatrix(amat, inva, 4);
    matrixVectorMultiply(inva, evec, abcd_coefs, 4, 4);

    // The cubic polynomial will have at most one local minimum.  Find the local minimum, and if
    // it is in the range [mvec[0] = 0.0, mvec[3]] (inclusive of the endpoints), then take it by
    // moving all coordinates that distance from their origin at the outset of the move.
    // Otherwise, compare evec[0] and evec[3] to find the new coordinates (leave them where they
    // are after completing the third move thus far, or reset them to where they started).  The
    // step size is adjusting the whole time.
    d_abcd[0] = 3.0 * abcd_coefs[0];
    d_abcd[1] = 2.0 * abcd_coefs[1];
    d_abcd[2] = abcd_coefs[2];
    const double sqrt_arg = (d_abcd[1] * d_abcd[1]) - (4.0 * d_abcd[0] * d_abcd[2]);
    if (sqrt_arg < 0.0 || (evec[0] < evec[1] && evec[0] < evec[2] && evec[0] < evec[3])) {
      
      // The cubic equation has no minima or maxima, or no move was successful in reducing the
      // energy.  Check the extrema of the range to ascertain the correct move.  If no move was
      // able to reduce the energy, reset the positions of all particles to the origin of the
      // line minimization process.
      if (evec[0] < evec[3]) {
        moveParticles<Tforce, Tcalc>(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, nullptr, nullptr,
                                     UnitCellType::NONE, vsk, vk.natom, -mvec[3] * move_scale,
                                     force_factor);
        move_scale *= 0.8;
      }
      else {
        move_scale *= 1.2;
      }
    }
    else {

      // Finish solving the quadratic formula to obtain both solutions and evaluate the cubic
      // polynomial's second derivative at each point.
      const double ext_i  = (-d_abcd[1] + sqrt(sqrt_arg)) / (2.0 * d_abcd[0]);
      const double ext_ii = (-d_abcd[1] - sqrt(sqrt_arg)) / (2.0 * d_abcd[0]);
      const double min_pos = ((2.0 * d_abcd[0] * ext_i) + d_abcd[1] > 0.0) ? ext_i : ext_ii;
      if (min_pos <= 0.0) {
        if (evec[3] > evec[0]) {
          moveParticles<Tforce, Tcalc>(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, nullptr, nullptr,
                                       UnitCellType::NONE, vsk, vk.natom, -mvec[3] * move_scale,
                                       force_factor);
          move_scale *= 0.8;
        }
        else {
          move_scale *= 1.2;
        }
      }
      else if (min_pos < mvec[3]) {

        // Make a final check that the minimum value of the function is an improvement over
        // the extrema. (This is only the value of the estimate function, not the actual value
        // of the system energy.) Otherwise, move the particles to a point within the range
        // scored by the four data points.
        const double epred = (((((abcd_coefs[0] * min_pos) + abcd_coefs[1]) * min_pos) +
                               abcd_coefs[2]) * min_pos) + abcd_coefs[3];
        if (evec[0] < epred) {
          moveParticles<Tforce, Tcalc>(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, nullptr, nullptr,
                                       UnitCellType::NONE, vsk, vk.natom, -mvec[3] * move_scale,
                                       force_factor);
          move_scale *= 0.8;
        }
        else if (epred < evec[3]) {
          moveParticles<Tforce, Tcalc>(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, nullptr, nullptr,
                                       UnitCellType::NONE, vsk, vk.natom,
                                       (min_pos - mvec[3]) * move_scale, force_factor);
          if (min_pos > 0.6) {
            move_scale *= 1.05;
          }
          else {
            move_scale *= 0.95;
          }
        }
        else {
          move_scale *= 1.05;
        }
      }
      else {
        move_scale *= 1.05;
      }
    }
  }

  // Perform one final energy and force evaluation, without clash checking (use the unmodified
  // potential function)
  switch (isk.igb) {
  case ImplicitSolventModel::NONE:
    evalNonbValeRestMM<Tcoord, Tforce,
                       Tcalc, Tcalc2, Tcalc4>(xcrd, ycrd, zcrd, nullptr, nullptr,
                                              UnitCellType::NONE, xfrc, yfrc, zfrc, &sc, vk, nbk,
                                              ser, rar, EvaluateForce::YES, 0,
                                              mincon.getTotalCycles());
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    evalRestrainedMMGB<Tcoord, Tforce,
                       Tcalc, Tcalc2, Tcalc4>(xcrd, ycrd, zcrd, nullptr, nullptr,
                                              UnitCellType::NONE, xfrc, yfrc, zfrc, &sc, vk, nbk,
                                              ser, isk, ngbk, effective_gb_radii.data(),
                                              psi.data(), sumdeijda.data(), rar,
                                              EvaluateForce::YES, 0, mincon.getTotalCycles());
  }
  transmitVirtualSiteForces<Tcalc, Tcalc>(xcrd, ycrd, zcrd, xfrc, yfrc, zfrc, nullptr, nullptr,
                                          UnitCellType::NONE, vsk);

  // Always commit the final state
  sc.commit(StateVariable::ALL_STATES);
  sc.incrementSampleCount();
  sc.setLastTimeStep(mincon.getTotalCycles());
  return sc;
}

} // namespace mm
} // namespace stormm
