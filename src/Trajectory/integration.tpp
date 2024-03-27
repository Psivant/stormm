// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void velocityVerletVelocityUpdate(const Tcoord* xvel, const Tcoord* yvel, const Tcoord* zvel,
                                  const Tcoord* xfrc, const Tcoord* yfrc, const Tcoord* zfrc,
                                  const int natom, const double* masses, Tcoord* vxalt,
                                  Tcoord* vyalt, Tcoord* vzalt,
                                  const ThermostatReader<Tcalc> &tstr,
                                  const Tcalc vel_scale_factor, const Tcalc frc_scale_factor) {
  const Tcalc zero = 0.0;
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const bool tcoord_is_integral = isSignedIntegralScalarType<Tcoord>();
  Tcalc dt_factor = (tcalc_is_double) ? 0.5  * tstr.dt * kcal_to_gafs :
                                        0.5f * tstr.dt * kcal_to_gafs_f;
  if (tcoord_is_integral) {
    dt_factor *= vel_scale_factor / frc_scale_factor;
  }
  switch (tstr.kind) {
  case ThermostatKind::NONE:
  case ThermostatKind::ANDERSEN:
  case ThermostatKind::BERENDSEN:
    for (int i = 0; i < natom; i++) {
      const Tcalc mss_i = masses[i];
      const Tcalc hmdt = (mss_i > constants::small) ? dt_factor / mss_i : zero;
      if (tcoord_is_integral) {
        const Tcalc xpush = hmdt * static_cast<Tcalc>(xfrc[i]);
        const Tcalc ypush = hmdt * static_cast<Tcalc>(yfrc[i]);
        const Tcalc zpush = hmdt * static_cast<Tcalc>(zfrc[i]);
        vxalt[i] = xvel[i] + static_cast<Tcoord>(llround(xpush));
        vyalt[i] = yvel[i] + static_cast<Tcoord>(llround(ypush));
        vzalt[i] = zvel[i] + static_cast<Tcoord>(llround(zpush));
      }
      else {
        vxalt[i] = xvel[i] + (hmdt * xfrc[i]);
        vyalt[i] = yvel[i] + (hmdt * yfrc[i]);
        vzalt[i] = zvel[i] + (hmdt * zfrc[i]);
      }
    }
    break;
  case ThermostatKind::LANGEVIN:
    switch (tstr.layout) {
    case ThermostatPartition::COMMON:
      {
        Tcalc sdfac;
        if (tcalc_is_double) {
          sdfac = sqrt(tstr.gamma_ln * boltzmann_constant *
                       getAtomTemperatureTarget<Tcalc>(tstr) / dt_factor);
        }
        else {
          sdfac = sqrtf(tstr.gamma_ln * boltzmann_constant_f *
                        getAtomTemperatureTarget<Tcalc>(tstr) / dt_factor);
        }
        if (tcoord_is_integral) {
          sdfac *= frc_scale_factor;
        }
        const size_t rnd_offset = static_cast<size_t>((tstr.step % tstr.depth) * 3) *
                                  static_cast<size_t>(tstr.padded_natom);
        for (int i = 0; i < natom; i++) {
          const Tcalc mss_i = masses[i];
          Tcalc rsd, hmdt;
          if (mss_i > constants::small) {
            rsd = (tcalc_is_double) ? sdfac * sqrt(mss_i) : sdfac * sqrtf(mss_i);
            hmdt = dt_factor / mss_i;
          }
          else {
            rsd = zero;
            hmdt = zero;
          }
          Tcalc xbump, ybump, zbump;
          switch (tstr.rng_mode) {
          case PrecisionModel::DOUBLE:
            xbump = rsd * tstr.cache[rnd_offset +                           i];
            ybump = rsd * tstr.cache[rnd_offset +      tstr.padded_natom  + i];
            zbump = rsd * tstr.cache[rnd_offset + (2 * tstr.padded_natom) + i];
            break;
          case PrecisionModel::SINGLE:
            xbump = rsd * tstr.sp_cache[rnd_offset +                           i];
            ybump = rsd * tstr.sp_cache[rnd_offset +      tstr.padded_natom  + i];
            zbump = rsd * tstr.sp_cache[rnd_offset + (2 * tstr.padded_natom) + i];
            break;
          }
          
          // The forces acting on each particle will not be modified with the Langevin bumps in
          // this case.  The same composite forces will be re-computed for the second half update.
          if (tcoord_is_integral) {
            const Tcoord ixbump = llround(xbump);
            const Tcoord iybump = llround(ybump);
            const Tcoord izbump = llround(zbump);
            const Tcalc fx_comp = xfrc[i] + ixbump;
            const Tcalc fy_comp = yfrc[i] + iybump;
            const Tcalc fz_comp = zfrc[i] + izbump;
            vxalt[i] = (xvel[i] + (hmdt * fx_comp)) * tstr.ln_implicit;
            vyalt[i] = (yvel[i] + (hmdt * fy_comp)) * tstr.ln_implicit;
            vzalt[i] = (zvel[i] + (hmdt * fz_comp)) * tstr.ln_implicit;
          }
          else {
            vxalt[i] = (xvel[i] + (hmdt * (xfrc[i] + xbump))) * tstr.ln_implicit;
            vyalt[i] = (yvel[i] + (hmdt * (yfrc[i] + ybump))) * tstr.ln_implicit;
            vzalt[i] = (zvel[i] + (hmdt * (zfrc[i] + zbump))) * tstr.ln_implicit;
          }
        }
      }
      break;
    case ThermostatPartition::SYSTEMS:
      rtErr("A single system cannot run with a thermostat configured for multiple systems.",
            "velocityVerletVelocityUpdate");
    case ThermostatPartition::ATOMS:
      {
        const size_t rnd_offset = static_cast<size_t>((tstr.step % tstr.depth) * 3) *
                                  static_cast<size_t>(tstr.padded_natom);
        for (int i = 0; i < tstr.npart; i++) {
          Tcalc sdfac;
          if (tcalc_is_double) {
            sdfac = sqrt(tstr.gamma_ln * boltzmann_constant *
                         getPartitionTemperatureTarget<Tcalc>(tstr, i) / dt_factor);
          }
          else {
            sdfac = sqrtf(tstr.gamma_ln * boltzmann_constant_f *
                          getPartitionTemperatureTarget<Tcalc>(tstr, i) / dt_factor);
          }
          if (tcoord_is_integral) {
            sdfac *= frc_scale_factor;
          }
          const int4 prtn = tstr.partition_bounds[i];
          for (int j = prtn.x; j < prtn.y; j++) {
            const Tcalc mss_j = masses[j];
            Tcalc rsd, hmdt;
            if (mss_j > constants::small) {
              rsd = (tcalc_is_double) ? sdfac * sqrt(mss_j) : sdfac * sqrtf(mss_j);
              hmdt = dt_factor / mss_j;
            }
            else {
              rsd = zero;
              hmdt = zero;
            }
            Tcalc xbump, ybump, zbump;
            switch (tstr.rng_mode) {
            case PrecisionModel::DOUBLE:
              xbump = rsd * tstr.cache[rnd_offset +                           j];
              ybump = rsd * tstr.cache[rnd_offset +      tstr.padded_natom  + j];
              zbump = rsd * tstr.cache[rnd_offset + (2 * tstr.padded_natom) + j];
              break;
            case PrecisionModel::SINGLE:
              xbump = rsd * tstr.sp_cache[rnd_offset +                           j];
              ybump = rsd * tstr.sp_cache[rnd_offset +      tstr.padded_natom  + j];
              zbump = rsd * tstr.sp_cache[rnd_offset + (2 * tstr.padded_natom) + j];
              break;
            }
            if (tcoord_is_integral) {
              const Tcoord ixbump = llround(xbump);
              const Tcoord iybump = llround(ybump);
              const Tcoord izbump = llround(zbump);
              const Tcalc fx_comp = xfrc[j] + ixbump;
              const Tcalc fy_comp = yfrc[j] + iybump;
              const Tcalc fz_comp = zfrc[j] + izbump;
              vxalt[j] = (xvel[j] + (hmdt * fx_comp)) * tstr.ln_implicit;
              vyalt[j] = (yvel[j] + (hmdt * fy_comp)) * tstr.ln_implicit;
              vzalt[j] = (zvel[j] + (hmdt * fz_comp)) * tstr.ln_implicit;
            }
            else {
              vxalt[j] = (xvel[j] + (hmdt * (xfrc[j] + xbump))) * tstr.ln_implicit;
              vyalt[j] = (yvel[j] + (hmdt * (yfrc[j] + ybump))) * tstr.ln_implicit;
              vzalt[j] = (zvel[j] + (hmdt * (zfrc[j] + zbump))) * tstr.ln_implicit;
            }
          }
        }
      }
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void velocityVerletCoordinateUpdate(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                                    const Tcoord* xfrc, const Tcoord* yfrc, const Tcoord* zfrc,
                                    const int natom, const double* masses, Tcoord* xalt,
                                    Tcoord* yalt, Tcoord* zalt, Tcoord* vxalt, Tcoord* vyalt,
                                    Tcoord* vzalt, const ThermostatReader<Tcalc> &tstr,
                                    const Tcalc gpos_scale_factor, const Tcalc vel_scale_factor,
                                    const Tcalc frc_scale_factor) {
  const Tcalc zero = 0.0;
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const bool tcoord_is_integral = isSignedIntegralScalarType<Tcoord>();
  Tcalc dt_factor = (tcalc_is_double) ? 0.5  * tstr.dt * kcal_to_gafs :
                                        0.5f * tstr.dt * kcal_to_gafs_f;
  if (tcoord_is_integral) {
    dt_factor *= vel_scale_factor / frc_scale_factor;
  }
  switch (tstr.kind) {
  case ThermostatKind::NONE:
  case ThermostatKind::ANDERSEN:
  case ThermostatKind::BERENDSEN:
    for (int i = 0; i < natom; i++) {
      const Tcalc mss_i = masses[i];
      const Tcalc hmdt = (mss_i > constants::small) ? dt_factor / mss_i : zero;
      if (tcoord_is_integral) {
        const Tcalc xpush = hmdt * static_cast<Tcalc>(xfrc[i]);
        const Tcalc ypush = hmdt * static_cast<Tcalc>(yfrc[i]);
        const Tcalc zpush = hmdt * static_cast<Tcalc>(zfrc[i]);
        vxalt[i] += static_cast<Tcoord>(llround(xpush));
        vyalt[i] += static_cast<Tcoord>(llround(ypush));
        vzalt[i] += static_cast<Tcoord>(llround(zpush));
      }
      else {
        vxalt[i] += hmdt * xfrc[i];
        vyalt[i] += hmdt * yfrc[i];
        vzalt[i] += hmdt * zfrc[i];
      }
    }
    break;
  case ThermostatKind::LANGEVIN:
    switch (tstr.layout) {
    case ThermostatPartition::COMMON:
      {
        Tcalc sdfac;
        if (tcalc_is_double) {
          sdfac = sqrt(tstr.gamma_ln * boltzmann_constant *
                       getAtomTemperatureTarget<Tcalc>(tstr) / dt_factor);
        }
        else {
          sdfac = sqrtf(tstr.gamma_ln * boltzmann_constant_f *
                        getAtomTemperatureTarget<Tcalc>(tstr) / dt_factor);
        }
        if (tcoord_is_integral) {
          sdfac *= frc_scale_factor;
        }
        const size_t rnd_offset = static_cast<size_t>((tstr.step % tstr.depth) * 3) *
                                  static_cast<size_t>(tstr.padded_natom);
        for (int i = 0; i < natom; i++) {
          const Tcalc mss_i = masses[i];
          Tcalc rsd, hmdt;
          if (mss_i > constants::small) {
            rsd = (tcalc_is_double) ? sdfac * sqrt(mss_i) : sdfac * sqrtf(mss_i);
            hmdt = dt_factor / mss_i;
          }
          else {
            rsd = zero;
            hmdt = zero;
          }
          Tcalc xbump, ybump, zbump;
          switch (tstr.rng_mode) {
          case PrecisionModel::DOUBLE:
            xbump = rsd * tstr.cache[rnd_offset +                           i];
            ybump = rsd * tstr.cache[rnd_offset +      tstr.padded_natom  + i];
            zbump = rsd * tstr.cache[rnd_offset + (2 * tstr.padded_natom) + i];
            break;
          case PrecisionModel::SINGLE:
            xbump = rsd * tstr.sp_cache[rnd_offset +                           i];
            ybump = rsd * tstr.sp_cache[rnd_offset +      tstr.padded_natom  + i];
            zbump = rsd * tstr.sp_cache[rnd_offset + (2 * tstr.padded_natom) + i];
            break;
          }
          if (tcoord_is_integral) {
            const Tcoord ixbump = llround(xbump);
            const Tcoord iybump = llround(ybump);
            const Tcoord izbump = llround(zbump);
            const Tcalc fx_comp = xfrc[i] + ixbump;
            const Tcalc fy_comp = yfrc[i] + iybump;
            const Tcalc fz_comp = zfrc[i] + izbump;
            vxalt[i] = (vxalt[i] * tstr.ln_explicit) + (hmdt * fx_comp);
            vyalt[i] = (vyalt[i] * tstr.ln_explicit) + (hmdt * fy_comp);
            vzalt[i] = (vzalt[i] * tstr.ln_explicit) + (hmdt * fz_comp);
          }
          else {
            vxalt[i] = (vxalt[i] * tstr.ln_explicit) + (hmdt * (xfrc[i] + xbump));
            vyalt[i] = (vyalt[i] * tstr.ln_explicit) + (hmdt * (yfrc[i] + ybump));
            vzalt[i] = (vzalt[i] * tstr.ln_explicit) + (hmdt * (zfrc[i] + zbump));
          }
        }
      }
      break;
    case ThermostatPartition::SYSTEMS:
      rtErr("A single system cannot run with a thermostat configured for multiple systems.",
            "velocityVerletVelocityUpdate");
    case ThermostatPartition::ATOMS:
      {
        const size_t rnd_offset = static_cast<size_t>((tstr.step % tstr.depth) * 3) *
                                  static_cast<size_t>(tstr.padded_natom);
        for (int i = 0; i < tstr.npart; i++) {
          Tcalc sdfac;
          if (tcalc_is_double) {
            sdfac = sqrt(tstr.gamma_ln * boltzmann_constant *
                         getPartitionTemperatureTarget<Tcalc>(tstr, i) / dt_factor);
          }
          else {
            sdfac = sqrtf(tstr.gamma_ln * boltzmann_constant_f *
                          getPartitionTemperatureTarget<Tcalc>(tstr, i) / dt_factor);
          }
          if (tcoord_is_integral) {
            sdfac *= frc_scale_factor;
          }
          const int4 prtn = tstr.partition_bounds[i];
          for (int j = prtn.x; j < prtn.y; j++) {
            const Tcalc mss_j = masses[j];
            Tcalc rsd, hmdt;
            if (mss_j > constants::small) {
              rsd = (tcalc_is_double) ? sdfac * sqrt(mss_j) : sdfac * sqrtf(mss_j);
              hmdt = dt_factor / mss_j;
            }
            else {
              rsd = zero;
              hmdt = zero;
            }
            Tcalc xbump, ybump, zbump;
            switch (tstr.rng_mode) {
            case PrecisionModel::DOUBLE:
              xbump = rsd * tstr.cache[rnd_offset +                           j];
              ybump = rsd * tstr.cache[rnd_offset +      tstr.padded_natom  + j];
              zbump = rsd * tstr.cache[rnd_offset + (2 * tstr.padded_natom) + j];
              break;
            case PrecisionModel::SINGLE:
              xbump = rsd * tstr.sp_cache[rnd_offset +                           j];
              ybump = rsd * tstr.sp_cache[rnd_offset +      tstr.padded_natom  + j];
              zbump = rsd * tstr.sp_cache[rnd_offset + (2 * tstr.padded_natom) + j];
              break;
            }
            if (tcoord_is_integral) {
              const Tcoord ixbump = llround(xbump);
              const Tcoord iybump = llround(ybump);
              const Tcoord izbump = llround(zbump);
              const Tcalc fx_comp = xfrc[j] + ixbump;
              const Tcalc fy_comp = yfrc[j] + iybump;
              const Tcalc fz_comp = zfrc[j] + izbump;
              vxalt[j] = (static_cast<Tcalc>(vxalt[j]) * tstr.ln_explicit) + (hmdt * fx_comp);
              vyalt[j] = (static_cast<Tcalc>(vyalt[j]) * tstr.ln_explicit) + (hmdt * fy_comp);
              vzalt[j] = (static_cast<Tcalc>(vzalt[j]) * tstr.ln_explicit) + (hmdt * fz_comp);
            }
            else {
              vxalt[j] = (vxalt[j] * tstr.ln_explicit) + (hmdt * (xfrc[j] + xbump));
              vyalt[j] = (vyalt[j] * tstr.ln_explicit) + (hmdt * (yfrc[j] + ybump));
              vzalt[j] = (vzalt[j] * tstr.ln_explicit) + (hmdt * (zfrc[j] + zbump));
            }
          }
        }
      }
      break;
    }
    break;
  }

  // With the complete velocity update, move the particles forward.
  if (tcoord_is_integral) {
    const Tcalc vcdt_factor = tstr.dt * gpos_scale_factor / vel_scale_factor;
    for (int i = 0; i < natom; i++) {
      xalt[i] = xcrd[i] + (vxalt[i] * vcdt_factor);
      yalt[i] = ycrd[i] + (vyalt[i] * vcdt_factor);
      zalt[i] = zcrd[i] + (vzalt[i] * vcdt_factor);
    }
  }
  else {
    for (int i = 0; i < natom; i++) {
      xalt[i] = xcrd[i] + (vxalt[i] * tstr.dt);
      yalt[i] = ycrd[i] + (vyalt[i] * tstr.dt);
      zalt[i] = zcrd[i] + (vzalt[i] * tstr.dt);
    }
  }
}

} // namespace structure
} // namespace stormm
