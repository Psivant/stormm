// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename T>
ThermostatWriter<T>::ThermostatWriter(const ThermostatKind kind_in,
                                      const ThermostatPartition layout_in, const int natom_in,
                                      const size_t padded_natom_in, const int npart_in,
                                      const int step_in, const int depth_in,
                                      const int init_evolution_in, const int end_evolution_in,
                                      const T init_temperature_in, const T final_temperature_in,
                                      const T dt_in, const bool cnst_geom_in,
                                      const T rattle_tol_in, const int rattle_iter_in,
                                      const int andr_cyc_in, const T gamma_ln_in,
                                      const T ln_implicit_in, const T ln_explicit_in,
                                      const int4* partition_bounds_in,
                                      const T* init_temperatures_in,
                                      const T* final_temperatures_in, ullint2* state_xy_in,
                                      ullint2* state_zw_in, const PrecisionModel rng_mode_in,
                                      double* cache_in, float* sp_cache_in) :
    kind{kind_in}, layout{layout_in}, natom{natom_in}, padded_natom{padded_natom_in},
    npart{npart_in}, step{step_in}, depth{depth_in}, init_evolution{init_evolution_in},
    end_evolution{end_evolution_in}, init_temperature{init_temperature_in},
    final_temperature{final_temperature_in}, dt{dt_in}, cnst_geom{cnst_geom_in},
    rattle_tol{rattle_tol_in}, rattle_iter{rattle_iter_in}, andr_cyc{andr_cyc_in},
    gamma_ln{gamma_ln_in}, ln_implicit{ln_implicit_in}, ln_explicit{ln_explicit_in},
    partition_bounds{partition_bounds_in}, init_temperatures{init_temperatures_in},
    final_temperatures{final_temperatures_in}, state_xy{state_xy_in}, state_zw{state_zw_in},
    rng_mode{rng_mode_in}, cache{cache_in}, sp_cache{sp_cache_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
ThermostatReader<T>::ThermostatReader(const ThermostatKind kind_in,
                                      const ThermostatPartition layout_in, const int natom_in,
                                      const size_t padded_natom_in, const int npart_in,
                                      const int step_in, const int depth_in,
                                      const int init_evolution_in, const int end_evolution_in,
                                      const T init_temperature_in, const T final_temperature_in,
                                      const T dt_in, const bool cnst_geom_in,
                                      const T rattle_tol_in, const int rattle_iter_in,
                                      const int andr_cyc_in, const T gamma_ln_in,
                                      const T ln_implicit_in, const T ln_explicit_in,
                                      const int4* partition_bounds_in,
                                      const T* init_temperatures_in,
                                      const T* final_temperatures_in, const ullint2* state_xy_in,
                                      const ullint2* state_zw_in, const PrecisionModel rng_mode_in,
                                      const double* cache_in, const float* sp_cache_in) :
    kind{kind_in}, layout{layout_in}, natom{natom_in}, padded_natom{padded_natom_in},
    npart{npart_in}, step{step_in}, depth{depth_in}, init_evolution{init_evolution_in},
    end_evolution{end_evolution_in}, init_temperature{init_temperature_in},
    final_temperature{final_temperature_in}, dt{dt_in}, cnst_geom{cnst_geom_in},
    rattle_tol{rattle_tol_in}, rattle_iter{rattle_iter_in}, andr_cyc{andr_cyc_in},
    gamma_ln{gamma_ln_in}, ln_implicit{ln_implicit_in}, ln_explicit{ln_explicit_in},
    partition_bounds{partition_bounds_in}, init_temperatures{init_temperatures_in},
    final_temperatures{final_temperatures_in}, state_xy{state_xy_in}, state_zw{state_zw_in},
    rng_mode{rng_mode_in}, cache{cache_in}, sp_cache{sp_cache_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
ThermostatReader<T>::ThermostatReader(const ThermostatWriter<T> &w) :
    kind{w.kind}, layout{w.layout}, natom{w.natom}, padded_natom{w.padded_natom}, npart{w.npart},
    step{w.step}, depth{w.depth}, init_evolution{w.init_evolution}, end_evolution{w.end_evolution},
    init_temperature{w.init_temperature}, final_temperature{w.final_temperature}, dt{w.dt},
    cnst_geom{w.cnst_geom}, rattle_tol{w.rattle_tol}, rattle_iter{w.rattle_iter},
    andr_cyc{w.andr_cyc}, gamma_ln{w.gamma_ln}, ln_implicit{w.ln_implicit},
    ln_explicit{w.ln_explicit}, init_temperatures{w.init_temperatures},
    final_temperatures{w.final_temperatures}, state_xy{w.state_xy}, state_zw{w.state_zw},
    rng_mode{w.rng_mode}, cache{w.cache}, sp_cache{w.sp_cache}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> Thermostat::getTemperatureSpread(const int step_number) const {
  const size_t nt = partition_count;
  std::vector<T> result(nt);
  switch (bath_layout) {
  case ThermostatPartition::COMMON:
    if (step_number <= initial_evolution_step) {
      result[0] = initial_temperature;
    }
    else if (step_number >= final_evolution_step) {
      result[0] = final_temperature;
    }
    else {
      const double evo = step_number - initial_evolution_step;
      const double prog = evo / static_cast<double>(final_evolution_step - initial_evolution_step);
      result[0] = (prog * (final_temperature - initial_temperature)) + initial_temperature;
    }
    break;
  case ThermostatPartition::SYSTEMS:
  case ThermostatPartition::ATOMS:
    const double* init_tvec = initial_temperatures.data();
    const double* finl_tvec = final_temperatures.data();
    if (step_number <= initial_evolution_step) {
      for (size_t i = 0; i < nt; i++) {
        result[i] = init_tvec[i];
      }
    }
    else if (step_number >= final_evolution_step) {
      for (size_t i = 0; i < nt; i++) {
        result[i] = finl_tvec[i];
      }
    }
    else {
      const double evo = step_number - initial_evolution_step;
      const double prog = evo / static_cast<double>(final_evolution_step - initial_evolution_step);
      for (size_t i = 0; i < nt; i++) {
        result[i] = (prog * (finl_tvec[i] - init_tvec[i])) + init_tvec[i];
      }
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc Thermostat::getAtomTarget(const int atom_index) const {
  if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
    const ThermostatReader<double> tstr = dpData();
    return getAtomTemperatureTarget<double>(tstr, atom_index);
  }
  else {
    const ThermostatReader<float> tstr = spData();
    return getAtomTemperatureTarget<float>(tstr, atom_index);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc Thermostat::getPartitionTarget(const int index) const {
  if (std::type_index(typeid(Tcalc)).hash_code() == double_type_index) {
    const ThermostatReader<double> tstr = dpData();
    return getPartitionTemperatureTarget<double>(tstr, index);
  }
  else {
    const ThermostatReader<float> tstr = spData();
    return getPartitionTemperatureTarget<float>(tstr, index);
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc getAtomTemperatureTarget(const ThermostatReader<Tcalc> &tstr, const int atom_index) {
  switch (tstr.layout) {
  case ThermostatPartition::COMMON:
    if (tstr.step <= tstr.init_evolution) {
      return tstr.init_temperature;
    }
    else if (tstr.step >= tstr.end_evolution) {
      return tstr.final_temperature;
    }
    else {
      const Tcalc evln = static_cast<Tcalc>(tstr.step - tstr.init_evolution) /
                         static_cast<Tcalc>(tstr.end_evolution - tstr.init_evolution);
      return tstr.init_temperature + (evln * (tstr.final_temperature - tstr.init_temperature));
    }
    break;
  case ThermostatPartition::SYSTEMS:
    
    break;
  case ThermostatPartition::ATOMS:
    {
      // The atom index is expected to have been supplied, but the partition index is not known.
      // Find the partition index, then the teraget temperature.  For setting temperatures or
      // Langvein bumps on the CPU host, a means of looping over all atoms in each compartment is
      // taken.  This function is used to query the target temperature for other purposes.
      if (atom_index < 0 || atom_index >= tstr.padded_natom) {
        rtErr("Atom index " + std::to_string(atom_index) + " is invalid for a thermostat "
              "regulating " + std::to_string(tstr.natom) + " atoms in all, even after considering "
              "padding.", "getAtomTemperatureTarget");
      }
      if (tstr.step <= tstr.init_evolution) {
        return tstr.init_temperatures[atom_index];
      }
      else if (tstr.step >= tstr.end_evolution) {
        return tstr.final_temperatures[atom_index];
      }
      else {
        const Tcalc evln = static_cast<Tcalc>(tstr.step - tstr.init_evolution) /
                           static_cast<Tcalc>(tstr.end_evolution - tstr.init_evolution);
        return tstr.init_temperatures[atom_index] + (evln * (tstr.final_temperatures[atom_index] -
                                                             tstr.init_temperatures[atom_index]));
      }
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc getAtomTemperatureTarget(const ThermostatWriter<Tcalc> &tstw, const int atom_index) {
  ThermostatReader tstr(tstw);
  return getAtomTemperatureTarget(tstr, atom_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc getPartitionTemperatureTarget(const ThermostatReader<Tcalc> &tstr, const int index) {
  switch (tstr.layout) {
  case ThermostatPartition::COMMON:
    if (tstr.step <= tstr.init_evolution) {
      return tstr.init_temperature;
    }
    else if (tstr.step >= tstr.end_evolution) {
      return tstr.final_temperature;
    }
    else {
      const Tcalc evln = static_cast<Tcalc>(tstr.step - tstr.init_evolution) /
                         static_cast<Tcalc>(tstr.end_evolution - tstr.init_evolution);
      return tstr.init_temperature + (evln * (tstr.final_temperature - tstr.init_temperature));
    }
    break;
  case ThermostatPartition::SYSTEMS:
  case ThermostatPartition::ATOMS:
    if (index < 0 || index >= tstr.npart) {
      rtErr("Partition index " + std::to_string(index) + " is invalid for a thermostat "
            "regulating " + std::to_string(tstr.npart) + " atoms in all, even after considering "
            "padding.", "getPartitionTemperatureTarget");
    }
    if (tstr.step <= tstr.init_evolution) {
      return tstr.init_temperatures[index];
    }
    else if (tstr.step >= tstr.end_evolution) {
      return tstr.final_temperatures[index];
    }
    else {
      const Tcalc evln = static_cast<Tcalc>(tstr.step - tstr.init_evolution) /
                         static_cast<Tcalc>(tstr.end_evolution - tstr.init_evolution);
      return tstr.init_temperatures[index] + (evln * (tstr.final_temperatures[index] -
                                                      tstr.init_temperatures[index]));
    }
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
Tcalc getPartitionTemperatureTarget(const ThermostatWriter<Tcalc> &tstw, const int index) {
  ThermostatReader tstr(tstw);
  return getPartitionTemperatureTarget(tstr, index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void andersenVelocityReset(Tcoord* xvel, Tcoord* yvel, Tcoord* zvel, int* xvel_ovrf,
                           int* yvel_ovrf, int* zvel_ovrf, const Tcalc* inv_masses,
                           const ThermostatReader<Tcalc> &tstr, const Tcalc vel_scale_factor) {
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  const bool tcoord_is_integral = isSignedIntegralScalarType<Tcoord>();
  const size_t natom_zu = tstr.natom;
  const size_t padded_natom_zu = tstr.padded_natom;
  
  // Determine the layer of the random cache to mine
  const size_t cache_offset = (tstr.andr_cyc > 0) ?
                              static_cast<size_t>((tstr.step % tstr.andr_cyc) * 3LLU) *
                              static_cast<size_t>(tstr.padded_natom) : 0;
  switch (tstr.layout) {
  case ThermostatPartition::COMMON:
    {
      const Tcalc tval = getAtomTemperatureTarget<Tcalc>(tstr);
      Tcalc ebeta;
      if (tcalc_is_double) {
        ebeta = (tval > constants::small) ? sqrt(boltzmann_constant_gafs * tval) : 0.0;
      }
      else {
        ebeta = (tval > constants::small) ? sqrtf(boltzmann_constant_gafs_f * tval) : 0.0;
      }
      for (size_t i = 0; i < natom_zu; i++) {
        const Tcalc mss_i = inv_masses[i];
        const Tcalc efac = (mss_i > constants::small) ? ebeta * sqrt(mss_i) : 0.0;
        Tcalc vx, vy, vz;
        switch (tstr.rng_mode) {
        case PrecisionModel::DOUBLE:
          vx = efac * tstr.cache[cache_offset +                            i];
          vy = efac * tstr.cache[cache_offset +         padded_natom_zu  + i];
          vz = efac * tstr.cache[cache_offset + (2LLU * padded_natom_zu) + i];
          break;
        case PrecisionModel::SINGLE:
          vx = efac * tstr.sp_cache[cache_offset +                            i];
          vy = efac * tstr.sp_cache[cache_offset +         padded_natom_zu  + i];
          vz = efac * tstr.sp_cache[cache_offset + (2LLU * padded_natom_zu) + i];
          break;
        }
        if (tcoord_is_integral) {
          if (tcalc_is_double) {
            const int95_t ivx = hostDoubleToInt95(vx * vel_scale_factor);
            const int95_t ivy = hostDoubleToInt95(vy * vel_scale_factor);
            const int95_t ivz = hostDoubleToInt95(vz * vel_scale_factor);
            xvel[i] = ivx.x;
            yvel[i] = ivy.x;
            zvel[i] = ivz.x;
            xvel_ovrf[i] = ivx.y;
            yvel_ovrf[i] = ivy.y;
            zvel_ovrf[i] = ivz.y;
          }
          else {
            xvel[i] = llround(vx * vel_scale_factor);
            yvel[i] = llround(vy * vel_scale_factor);
            zvel[i] = llround(vz * vel_scale_factor);
          }
        }
        else {
          xvel[i] = vx;
          yvel[i] = vy;
          zvel[i] = vz;
        }
      }
    }
    break;
  case ThermostatPartition::SYSTEMS:
  case ThermostatPartition::ATOMS:
    for (int i = 0; i < tstr.npart; i++) {
      const Tcalc tval = getPartitionTemperatureTarget<Tcalc>(tstr, i);
      Tcalc ebeta;
      if (tcalc_is_double) {
        ebeta = (tval > constants::small) ? sqrt(boltzmann_constant_gafs * tval) : 0.0;
      }
      else {
        ebeta = (tval > constants::small) ? sqrtf(boltzmann_constant_gafs_f * tval) : 0.0;
      }
      const int4 prtn = tstr.partition_bounds[i];
      for (int j = prtn.x; j < prtn.y; j++) {
        const Tcalc mss_j = inv_masses[j];
        const Tcalc efac = (mss_j > constants::small) ? ebeta * sqrt(mss_j) : 0.0;
        Tcalc vx, vy, vz;
        switch (tstr.rng_mode) {
        case PrecisionModel::DOUBLE:
          vx = efac * tstr.cache[cache_offset +                            j];
          vy = efac * tstr.cache[cache_offset +         padded_natom_zu  + j];
          vz = efac * tstr.cache[cache_offset + (2LLU * padded_natom_zu) + j];
          break;
        case PrecisionModel::SINGLE:
          vx = efac * tstr.sp_cache[cache_offset +                            j];
          vy = efac * tstr.sp_cache[cache_offset +         padded_natom_zu  + j];
          vz = efac * tstr.sp_cache[cache_offset + (2LLU * padded_natom_zu) + j];
          break;
        }
        if (tcoord_is_integral) {
          if (tcalc_is_double) {
            const int95_t ivx = hostDoubleToInt95(vx * vel_scale_factor);
            const int95_t ivy = hostDoubleToInt95(vy * vel_scale_factor);
            const int95_t ivz = hostDoubleToInt95(vz * vel_scale_factor);
            xvel[j] = ivx.x;
            yvel[j] = ivy.x;
            zvel[j] = ivz.x;
            xvel_ovrf[j] = ivx.y;
            yvel_ovrf[j] = ivy.y;
            zvel_ovrf[j] = ivz.y;
          }
          else {
            xvel[j] = llround(vx * vel_scale_factor);
            yvel[j] = llround(vy * vel_scale_factor);
            zvel[j] = llround(vz * vel_scale_factor);
          }
        }
        else {
          xvel[j] = vx;
          yvel[j] = vy;
          zvel[j] = vz;
        }
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void andersenVelocityReset(Tcoord* xvel, Tcoord* yvel, Tcoord* zvel, const Tcalc* inv_masses,
                           const ThermostatReader<Tcalc> &tstr, const Tcalc vel_scale_factor) {
  andersenVelocityReset<Tcoord, Tcalc>(xvel, yvel, zvel, nullptr, nullptr, nullptr, inv_masses,
                                       tstr, vel_scale_factor);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc, typename Tcalc2, typename Tcalc4>
void andersenVelocityReset(PsSynthesisWriter *poly_psw,
                           const SyAtomUpdateKit<Tcalc, Tcalc2, Tcalc4> &poly_auk,
                           const ThermostatReader<Tcalc> &tstr) {
  andersenVelocityReset<llint, Tcalc>(poly_psw->xvel, poly_psw->yvel, poly_psw->zvel,
                                      poly_psw->xvel_ovrf, poly_psw->yvel_ovrf,
                                      poly_psw->zvel_ovrf, poly_auk.inv_masses, tstr,
                                      poly_psw->vel_scale);

}

} // namespace trajectory
} // namespace stormm
