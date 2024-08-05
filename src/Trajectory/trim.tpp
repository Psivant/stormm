// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename T>
bool checkInertialTensor(const std::vector<T> &inrt, const int natom,
                         const ExceptionResponse policy) {
  const T det_inert = std::abs(leibnizDeterminant(inrt));
  if (det_inert  < constants::small) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The rotational moment of a system containing " + std::to_string(natom) + " particles "
            "is zero along one or more axes, as detected by the determinant of the inertial "
            "moment tensor being " + realToString(det_inert, 10, 3, NumberFormat::SCIENTIFIC) +
            ".", "removeMomentum");
    case ExceptionResponse::WARN:
      rtWarn("The rotational moment of a system containing " + std::to_string(natom) +
             " particles is zero along one or more axes, as detected by the determinant of the "
             "inertial moment tensor being " +
             realToString(det_inert, 10, 3, NumberFormat::SCIENTIFIC) + ".  No removal of "
             "rotational moentum will be attempted.", "removeMomentum");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    return false;
  }
  else {
    return true;
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tmass, typename Tcalc>
void removeMomentum(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, int* xcrd_ovrf, int* ycrd_ovrf,
                    int* zcrd_ovrf, Tcoord* xvel, Tcoord* yvel, Tcoord* zvel, int* xvel_ovrf,
                    int* yvel_ovrf, int* zvel_ovrf, const Tmass* masses,
                    const UnitCellType unit_cell, const int natom, const Tcalc gpos_scale,
                    const Tcalc vel_scale, const ExceptionResponse policy) {

  // Avoid dividing by zero.  Return immediately if there are no particles.
  if (natom == 0) {
    return;
  }
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const bool tcoord_is_integral = isSignedIntegralScalarType<Tcoord>();
  const Tcalc zero = 0.0;
  const Tcalc one = 1.0;
  const Tcalc inv_gpos_scale = one / gpos_scale;
  const Tcalc inv_vel_scale = one / vel_scale;
  const Tcalc half_pi = (tcalc_is_double) ? 0.5 * pi : 0.5 * pi_f;
  Tcalc momx  = zero;
  Tcalc momy  = zero;
  Tcalc momz  = zero;
  Tcalc tmass = zero;
  for (int i = 0; i < natom; i++) {
    const Tcalc mss_i = masses[i];
    if (tcoord_is_integral) {
      Tcalc vx, vy, vz;
      if (tcalc_is_double) {

        // Extended coordinates and double-precision calculations are implicitly linked, as in
        // other calculations and GPU kernels involving the PhaseSpaceSynthesis.
        vx = hostInt95ToDouble(xvel[i], xvel_ovrf[i]) * inv_vel_scale;
        vy = hostInt95ToDouble(yvel[i], yvel_ovrf[i]) * inv_vel_scale;
        vz = hostInt95ToDouble(zvel[i], zvel_ovrf[i]) * inv_vel_scale;
      }
      else {
        vx = static_cast<Tcalc>(xvel[i]) * inv_vel_scale;
        vy = static_cast<Tcalc>(yvel[i]) * inv_vel_scale;
        vz = static_cast<Tcalc>(zvel[i]) * inv_vel_scale;
      }
      momx += mss_i * vx;
      momy += mss_i * vy;
      momz += mss_i * vz;
    }
    else {
      momx += mss_i * xvel[i];
      momy += mss_i * yvel[i];
      momz += mss_i * zvel[i];
    }

    // Masses are always expressed as real numbers
    tmass += mss_i;
  }

  // Normalize the system's net velocity by its total mass.  If there is no mass (which should
  // never happen), bail out.
  if (tmass < constants::tiny) {
    return;
  }
  tmass = one / tmass;
  momx *= tmass;
  momy *= tmass;
  momz *= tmass;
  
  // Subtract the net velocity from all particles
  if (tcoord_is_integral) {
    if (tcalc_is_double) {
      const int95_t idvx = hostDoubleToInt95(momx * vel_scale);
      const int95_t idvy = hostDoubleToInt95(momy * vel_scale);
      const int95_t idvz = hostDoubleToInt95(momz * vel_scale);
      for (int i = 0; i < natom; i++) {
        const int95_t vx_update = hostInt95Subtract(xvel[i], xvel_ovrf[i], idvx.x, idvx.y);
        const int95_t vy_update = hostInt95Subtract(yvel[i], yvel_ovrf[i], idvy.x, idvy.y);
        const int95_t vz_update = hostInt95Subtract(zvel[i], zvel_ovrf[i], idvz.x, idvz.y);
        xvel[i] = vx_update.x;
        yvel[i] = vy_update.x;
        zvel[i] = vz_update.x;
        xvel_ovrf[i] = vx_update.y;
        yvel_ovrf[i] = vy_update.y;
        zvel_ovrf[i] = vz_update.y;
      }
    }
    else {
      const llint idvx = llround(momx * vel_scale);
      const llint idvy = llround(momy * vel_scale);
      const llint idvz = llround(momz * vel_scale);
      for (int i = 0; i < natom; i++) {
        xvel[i] -= idvx;
        yvel[i] -= idvy;
        zvel[i] -= idvz;
      }
    }
  }
  else {
    for (int i = 0; i < natom; i++) {
      xvel[i] -= momx;
      yvel[i] -= momy;
      zvel[i] -= momz;
    }
  }

  // Remove angular momentum, if applicable.  Systems in periodic boundary conditions are exempt.
  switch (unit_cell) {
  case UnitCellType::NONE:
    {
      // Compute the center of mass and moment of inertia.
      Tcalc comx = zero;
      Tcalc comy = zero;
      Tcalc comz = zero;
      for (int i = 0; i < natom; i++) {
        const Tcalc mss_i = masses[i];
        if (tcoord_is_integral) {
          if (tcalc_is_double) {
            comx += mss_i * hostInt95ToDouble(xcrd[i], xcrd_ovrf[i]) * inv_gpos_scale;
            comy += mss_i * hostInt95ToDouble(ycrd[i], ycrd_ovrf[i]) * inv_gpos_scale;
            comz += mss_i * hostInt95ToDouble(zcrd[i], zcrd_ovrf[i]) * inv_gpos_scale;
          }
          else {
            comx += mss_i * static_cast<Tcalc>(xcrd[i]) * inv_gpos_scale;
            comy += mss_i * static_cast<Tcalc>(ycrd[i]) * inv_gpos_scale;
            comz += mss_i * static_cast<Tcalc>(zcrd[i]) * inv_gpos_scale;
          }
        }
        else {
          comx += mss_i * xcrd[i];
          comy += mss_i * ycrd[i];
          comz += mss_i * zcrd[i];
        }
      }

      // Re-use the total inverse mass computed earlier
      comx *= tmass;
      comy *= tmass;
      comz *= tmass;

      // Move the system's center of mass to the origin.
      if (tcoord_is_integral) {
        if (tcalc_is_double) {
          const int95_t icomx = hostDoubleToInt95(comx * gpos_scale);
          const int95_t icomy = hostDoubleToInt95(comy * gpos_scale);
          const int95_t icomz = hostDoubleToInt95(comz * gpos_scale);
          for (int i = 0; i < natom; i++) {
            const int95_t ix_disp = hostInt95Subtract(xcrd[i], xcrd_ovrf[i], icomx.x, icomx.y);
            const int95_t iy_disp = hostInt95Subtract(ycrd[i], ycrd_ovrf[i], icomy.x, icomy.y);
            const int95_t iz_disp = hostInt95Subtract(zcrd[i], zcrd_ovrf[i], icomz.x, icomz.y);
            xcrd[i] = ix_disp.x;
            ycrd[i] = iy_disp.x;
            zcrd[i] = iz_disp.x;
            xcrd_ovrf[i] = ix_disp.y;
            ycrd_ovrf[i] = iy_disp.y;
            zcrd_ovrf[i] = iz_disp.y;
          }
        }
        else {
          const llint icomx = llround(comx * gpos_scale);
          const llint icomy = llround(comy * gpos_scale);
          const llint icomz = llround(comz * gpos_scale);
          for (int i = 0; i < natom; i++) {
            xcrd[i] -= icomx;
            ycrd[i] -= icomy;
            zcrd[i] -= icomz;
          }
        }
      }
      else {
        for (int i = 0; i < natom; i++) {
          xcrd[i] -= comx;
          ycrd[i] -= comy;
          zcrd[i] -= comz;
        }
      }
      
      // Compute the angular momentum, given the center of mass.  With net translational momentum
      // removed, the calculation is properly grounded.
      Tcalc r[3], v[3], rcv[3];
      Tcalc ang_mom[3] = { zero, zero, zero };
      std::vector<Tcalc> inertia(9, zero);
      Tcalc xx = 0.0;
      Tcalc xy = 0.0;
      Tcalc xz = 0.0;
      Tcalc yy = 0.0;
      Tcalc yz = 0.0;
      Tcalc zz = 0.0;
      if (tcoord_is_integral) {
        if (tcalc_is_double) {
          const int95_t icomx = hostDoubleToInt95(comx * gpos_scale);
          const int95_t icomy = hostDoubleToInt95(comy * gpos_scale);
          const int95_t icomz = hostDoubleToInt95(comz * gpos_scale);
          for (int i = 0; i < natom; i++) {
            const Tcalc mss_i = masses[i];
            r[0] = hostInt95ToDouble(xcrd[i], xcrd_ovrf[i]) * inv_gpos_scale;
            r[1] = hostInt95ToDouble(ycrd[i], ycrd_ovrf[i]) * inv_gpos_scale;
            r[2] = hostInt95ToDouble(zcrd[i], zcrd_ovrf[i]) * inv_gpos_scale;
            v[0] = hostInt95ToDouble(xvel[i], xvel_ovrf[i]) * inv_vel_scale;
            v[1] = hostInt95ToDouble(yvel[i], yvel_ovrf[i]) * inv_vel_scale;
            v[2] = hostInt95ToDouble(zvel[i], zvel_ovrf[i]) * inv_vel_scale;
            crossProduct<Tcalc>(r, v, rcv);
            ang_mom[0] += mss_i * rcv[0];
            ang_mom[1] += mss_i * rcv[1];
            ang_mom[2] += mss_i * rcv[2];
            xx += mss_i * r[0] * r[0];
            xy += mss_i * r[0] * r[1];
            xz += mss_i * r[0] * r[2];
            yy += mss_i * r[1] * r[1];
            yz += mss_i * r[1] * r[2];
            zz += mss_i * r[2] * r[2];
          }
        }
        else {
          const llint icomx = llround(comx * gpos_scale);
          const llint icomy = llround(comy * gpos_scale);
          const llint icomz = llround(comz * gpos_scale);
          for (int i = 0; i < natom; i++) {
            const Tcalc mss_i = masses[i];
            r[0] = static_cast<Tcalc>(xcrd[i]) * inv_gpos_scale;
            r[1] = static_cast<Tcalc>(ycrd[i]) * inv_gpos_scale;
            r[2] = static_cast<Tcalc>(zcrd[i]) * inv_gpos_scale;
            v[0] = static_cast<Tcalc>(xvel[i]) * inv_vel_scale;
            v[1] = static_cast<Tcalc>(yvel[i]) * inv_vel_scale;
            v[2] = static_cast<Tcalc>(zvel[i]) * inv_vel_scale;
            crossProduct<Tcalc>(r, v, rcv);
            ang_mom[0] += mss_i * rcv[0];
            ang_mom[1] += mss_i * rcv[1];
            ang_mom[2] += mss_i * rcv[2];
            xx += mss_i * r[0] * r[0];
            xy += mss_i * r[0] * r[1];
            xz += mss_i * r[0] * r[2];
            yy += mss_i * r[1] * r[1];
            yz += mss_i * r[1] * r[2];
            zz += mss_i * r[2] * r[2];
          }
        }
      }
      else {
        for (int i = 0; i < natom; i++) {
          const Tcalc mss_i = masses[i];
          r[0] = xcrd[i];
          r[1] = ycrd[i];
          r[2] = zcrd[i];
          v[0] = xvel[i];
          v[1] = yvel[i];
          v[2] = zvel[i];
          crossProduct<Tcalc>(r, v, rcv);
          ang_mom[0] += mss_i * rcv[0];
          ang_mom[1] += mss_i * rcv[1];
          ang_mom[2] += mss_i * rcv[2];
          xx += mss_i * r[0] * r[0];
          xy += mss_i * r[0] * r[1];
          xz += mss_i * r[0] * r[2];
          yy += mss_i * r[1] * r[1];
          yz += mss_i * r[1] * r[2];
          zz += mss_i * r[2] * r[2];
        }
      }
      inertia[0] = yy + zz;
      inertia[1] = -xy;
      inertia[2] = -xz;
      inertia[3] = -xy;
      inertia[4] = xx + zz;
      inertia[5] = -yz;
      inertia[6] = -xz;
      inertia[7] = -yz;
      inertia[8] = xx + yy;

      // The rotational velocity omega may be defined by the angular momentum, the sum of r x v for
      // all particles, divided by the moment of inertia I.  As the moment of inertia is a tensor,
      // invert the matrix in order to calculate I^(-1) L = omega.
      std::vector<Tcalc> inv_inertia(9);
      if (checkInertialTensor(inertia, natom, policy) == false) {
        return;
      }
      
      // This code will only be reached if the intertial tensor is invertible.
      invertSquareMatrix(inertia, &inv_inertia);
      const Tcalc rot_vel[3] = { (inv_inertia[0] * ang_mom[0]) + (inv_inertia[3] * ang_mom[1]) +
                                 (inv_inertia[6] * ang_mom[2]),
                                 (inv_inertia[1] * ang_mom[0]) + (inv_inertia[4] * ang_mom[1]) +
                                 (inv_inertia[7] * ang_mom[2]),
                                 (inv_inertia[2] * ang_mom[0]) + (inv_inertia[5] * ang_mom[1]) +
                                 (inv_inertia[8] * ang_mom[2]) };
      Tcalc vcr[3];      
      if (tcoord_is_integral) {
        if (tcalc_is_double) {

          // In order to clear angular velocity about the Cartesian x, y, or z axes, add the
          // {y, z}, {x, z}, or {x, y} velocities implied by the particle's moment arm about the
          // same two axes and the negative of the net rotational velocity.
          for (int i = 0; i < natom; i++) {
            r[0] = hostInt95ToDouble(xcrd[i], xcrd_ovrf[i]) * inv_gpos_scale;
            r[1] = hostInt95ToDouble(ycrd[i], ycrd_ovrf[i]) * inv_gpos_scale;
            r[2] = hostInt95ToDouble(zcrd[i], zcrd_ovrf[i]) * inv_gpos_scale;
            crossProduct(rot_vel, r, vcr);
            const int95_t xvn = hostInt95Sum(xvel[i], xvel_ovrf[i], -vcr[0] * vel_scale);
            const int95_t yvn = hostInt95Sum(yvel[i], yvel_ovrf[i], -vcr[1] * vel_scale);
            const int95_t zvn = hostInt95Sum(zvel[i], zvel_ovrf[i], -vcr[2] * vel_scale);
            xvel[i] = xvn.x;
            yvel[i] = yvn.x;
            zvel[i] = zvn.x;
            xvel_ovrf[i] = xvn.y;
            yvel_ovrf[i] = yvn.y;
            zvel_ovrf[i] = zvn.y;
          }
        }
        else {
          const llint icomx = llround(comx * gpos_scale);
          const llint icomy = llround(comy * gpos_scale);
          const llint icomz = llround(comz * gpos_scale);
          for (int i = 0; i < natom; i++) {
            r[0] = static_cast<Tcalc>(xcrd[i]) * inv_gpos_scale;
            r[1] = static_cast<Tcalc>(ycrd[i]) * inv_gpos_scale;
            r[2] = static_cast<Tcalc>(zcrd[i]) * inv_gpos_scale;
            crossProduct(rot_vel, r, vcr);
            xvel[i] -= llround(vcr[0] * vel_scale);
            yvel[i] -= llround(vcr[1] * vel_scale);
            zvel[i] -= llround(vcr[2] * vel_scale);
          }
        }
      }
      else {
        for (int i = 0; i < natom; i++) {
          r[0] = xcrd[i];
          r[1] = ycrd[i];
          r[2] = zcrd[i];
          crossProduct(rot_vel, r, vcr);
          xvel[i] -= vcr[0];
          yvel[i] -= vcr[1];
          zvel[i] -= vcr[2];
        }
      }
    }
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    break;
  }
}
  
} // namespace trajectory
} // namespace stormm
