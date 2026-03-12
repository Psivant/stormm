"""ctypes-based Python API for selected STORMM dynamics workflows.

Design notes for contributors:
- This module is intentionally thin.  It forwards work to compiled STORMM C++ code.
- All raw handles are opaque pointers created/destroyed by the bridge C API.
- We keep explicit argtypes/restype declarations to make failures deterministic.
"""

from __future__ import annotations

import ctypes
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union


class StormmError(RuntimeError):
  """Raised when the bridge reports an execution or validation error."""

  pass


def _shared_lib_extension() -> str:
  """Map platform name to shared-library suffix."""

  if sys.platform == "darwin":
    return ".dylib"
  if sys.platform.startswith("win"):
    return ".dll"
  return ".so"


def _candidate_bridge_paths() -> List[Path]:
  """Search order for the C bridge shared library.

  STORMM_PYBRIDGE can override the location, otherwise we check the common build path and
  package-local fallbacks.
  """

  ext = _shared_lib_extension()
  repo_root = Path(__file__).resolve().parents[2]
  env_override = os.environ.get("STORMM_PYBRIDGE")

  candidates: List[Path] = []
  if env_override:
    candidates.append(Path(env_override).expanduser())
  candidates.append(repo_root / "build" / f"libstormm_pybridge{ext}")
  candidates.append(repo_root / "build" / f"stormm_pybridge{ext}")
  candidates.append(Path(__file__).resolve().parent / f"libstormm_pybridge{ext}")
  candidates.append(Path(__file__).resolve().parent / f"stormm_pybridge{ext}")
  return candidates


def _load_library() -> ctypes.CDLL:
  """Load the compiled bridge library or fail with actionable diagnostics."""

  checked_paths: List[str] = []
  for candidate in _candidate_bridge_paths():
    checked_paths.append(str(candidate))
    if candidate.exists():
      return ctypes.CDLL(str(candidate))
  joined = "\n".join(checked_paths)
  raise StormmError(
      "Could not locate the STORMM Python bridge library.\n"
      "Build it with: cmake --build build --target stormm_pybridge\n"
      "Searched paths:\n" + joined
  )


_lib = _load_library()

# ------------------------------------------------------------------------------
# C ABI signatures
# ------------------------------------------------------------------------------
# Explicit signatures keep ctypes from making unsafe type assumptions.
_lib.stormm_get_last_error.restype = ctypes.c_char_p

_lib.stormm_state_variable_count.argtypes = []
_lib.stormm_state_variable_count.restype = ctypes.c_int
_lib.stormm_state_variable_name.argtypes = [ctypes.c_int]
_lib.stormm_state_variable_name.restype = ctypes.c_char_p

# AtomGraph constructor / query signatures.
_lib.stormm_atomgraph_create.argtypes = [ctypes.c_char_p]
_lib.stormm_atomgraph_create.restype = ctypes.c_void_p
_lib.stormm_atomgraph_destroy.argtypes = [ctypes.c_void_p]
_lib.stormm_atomgraph_destroy.restype = ctypes.c_int
_lib.stormm_atomgraph_get_atom_count.argtypes = [ctypes.c_void_p]
_lib.stormm_atomgraph_get_atom_count.restype = ctypes.c_int
_lib.stormm_atomgraph_get_residue_count.argtypes = [ctypes.c_void_p]
_lib.stormm_atomgraph_get_residue_count.restype = ctypes.c_int
_lib.stormm_atomgraph_get_total_mass.argtypes = [ctypes.c_void_p]
_lib.stormm_atomgraph_get_total_mass.restype = ctypes.c_double
_lib.stormm_atomgraph_get_atomic_masses.argtypes = [ctypes.c_void_p,
                                                     ctypes.POINTER(ctypes.c_double), ctypes.c_int]
_lib.stormm_atomgraph_get_atomic_masses.restype = ctypes.c_int

# PhaseSpace constructor / query signatures.
_lib.stormm_phasespace_create.argtypes = [ctypes.c_char_p, ctypes.c_void_p]
_lib.stormm_phasespace_create.restype = ctypes.c_void_p
_lib.stormm_phasespace_destroy.argtypes = [ctypes.c_void_p]
_lib.stormm_phasespace_destroy.restype = ctypes.c_int
_lib.stormm_phasespace_get_atom_count.argtypes = [ctypes.c_void_p]
_lib.stormm_phasespace_get_atom_count.restype = ctypes.c_int
_lib.stormm_phasespace_get_interlaced_coordinates.argtypes = [ctypes.c_void_p, ctypes.c_int,
                                                               ctypes.POINTER(ctypes.c_double),
                                                               ctypes.c_int]
_lib.stormm_phasespace_get_interlaced_coordinates.restype = ctypes.c_int
_lib.stormm_phasespace_export.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_double]
_lib.stormm_phasespace_export.restype = ctypes.c_int

# DynamicsControls mutator signatures.
_lib.stormm_dynamics_controls_create.argtypes = []
_lib.stormm_dynamics_controls_create.restype = ctypes.c_void_p
_lib.stormm_dynamics_controls_destroy.argtypes = [ctypes.c_void_p]
_lib.stormm_dynamics_controls_destroy.restype = ctypes.c_int
_lib.stormm_dynamics_controls_get_step_count.argtypes = [ctypes.c_void_p]
_lib.stormm_dynamics_controls_get_step_count.restype = ctypes.c_int
_lib.stormm_dynamics_controls_set_step_count.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.stormm_dynamics_controls_set_step_count.restype = ctypes.c_int
_lib.stormm_dynamics_controls_set_time_step.argtypes = [ctypes.c_void_p, ctypes.c_double]
_lib.stormm_dynamics_controls_set_time_step.restype = ctypes.c_int
_lib.stormm_dynamics_controls_set_diagnostic_frequency.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.stormm_dynamics_controls_set_diagnostic_frequency.restype = ctypes.c_int
_lib.stormm_dynamics_controls_set_trajectory_frequency.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.stormm_dynamics_controls_set_trajectory_frequency.restype = ctypes.c_int
_lib.stormm_dynamics_controls_set_com_purge_frequency.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.stormm_dynamics_controls_set_com_purge_frequency.restype = ctypes.c_int
_lib.stormm_dynamics_controls_set_thermostat_kind.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.stormm_dynamics_controls_set_thermostat_kind.restype = ctypes.c_int
_lib.stormm_dynamics_controls_set_temperature.argtypes = [ctypes.c_void_p, ctypes.c_double,
                                                          ctypes.c_double]
_lib.stormm_dynamics_controls_set_temperature.restype = ctypes.c_int
_lib.stormm_dynamics_controls_set_langevin_frequency.argtypes = [ctypes.c_void_p, ctypes.c_double]
_lib.stormm_dynamics_controls_set_langevin_frequency.restype = ctypes.c_int
_lib.stormm_dynamics_controls_set_andersen_frequency.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.stormm_dynamics_controls_set_andersen_frequency.restype = ctypes.c_int
_lib.stormm_dynamics_controls_set_seed.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.stormm_dynamics_controls_set_seed.restype = ctypes.c_int

# One-system dynamics execution signature.
_lib.stormm_dynamics_run.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int,
                                     ctypes.c_double, ctypes.c_char_p, ctypes.c_char_p,
                                     ctypes.POINTER(ctypes.c_double),
                                     ctypes.POINTER(ctypes.c_double),
                                     ctypes.POINTER(ctypes.c_int),
                                     ctypes.POINTER(ctypes.c_double), ctypes.c_int]
_lib.stormm_dynamics_run.restype = ctypes.c_int

# AtomGraphSynthesis signatures.
_lib.stormm_atomgraph_synthesis_create.argtypes = [ctypes.POINTER(ctypes.c_void_p), ctypes.c_int]
_lib.stormm_atomgraph_synthesis_create.restype = ctypes.c_void_p
_lib.stormm_atomgraph_synthesis_destroy.argtypes = [ctypes.c_void_p]
_lib.stormm_atomgraph_synthesis_destroy.restype = ctypes.c_int
_lib.stormm_atomgraph_synthesis_get_system_count.argtypes = [ctypes.c_void_p]
_lib.stormm_atomgraph_synthesis_get_system_count.restype = ctypes.c_int
_lib.stormm_atomgraph_synthesis_get_unique_topology_count.argtypes = [ctypes.c_void_p]
_lib.stormm_atomgraph_synthesis_get_unique_topology_count.restype = ctypes.c_int
_lib.stormm_atomgraph_synthesis_get_atom_count.argtypes = [ctypes.c_void_p]
_lib.stormm_atomgraph_synthesis_get_atom_count.restype = ctypes.c_int
_lib.stormm_atomgraph_synthesis_get_atom_count_system.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.stormm_atomgraph_synthesis_get_atom_count_system.restype = ctypes.c_int

# PhaseSpaceSynthesis signatures.
_lib.stormm_phasespace_synthesis_create.argtypes = [ctypes.POINTER(ctypes.c_void_p),
                                                    ctypes.POINTER(ctypes.c_void_p), ctypes.c_int]
_lib.stormm_phasespace_synthesis_create.restype = ctypes.c_void_p
_lib.stormm_phasespace_synthesis_destroy.argtypes = [ctypes.c_void_p]
_lib.stormm_phasespace_synthesis_destroy.restype = ctypes.c_int
_lib.stormm_phasespace_synthesis_get_system_count.argtypes = [ctypes.c_void_p]
_lib.stormm_phasespace_synthesis_get_system_count.restype = ctypes.c_int
_lib.stormm_phasespace_synthesis_get_unique_topology_count.argtypes = [ctypes.c_void_p]
_lib.stormm_phasespace_synthesis_get_unique_topology_count.restype = ctypes.c_int
_lib.stormm_phasespace_synthesis_get_padded_atom_count.argtypes = [ctypes.c_void_p]
_lib.stormm_phasespace_synthesis_get_padded_atom_count.restype = ctypes.c_int
_lib.stormm_phasespace_synthesis_get_atom_count.argtypes = [ctypes.c_void_p, ctypes.c_int]
_lib.stormm_phasespace_synthesis_get_atom_count.restype = ctypes.c_int
_lib.stormm_phasespace_synthesis_get_interlaced_coordinates.argtypes = [ctypes.c_void_p,
                                                                         ctypes.c_int, ctypes.c_int,
                                                                         ctypes.POINTER(ctypes.c_double),
                                                                         ctypes.c_int]
_lib.stormm_phasespace_synthesis_get_interlaced_coordinates.restype = ctypes.c_int


def _last_error_message() -> str:
  """Fetch the bridge's thread-local error string."""

  raw = _lib.stormm_get_last_error()
  if raw is None:
    return "Unknown STORMM error."
  return raw.decode("utf-8")


def _check_ptr(value: ctypes.c_void_p, context: str) -> ctypes.c_void_p:
  """Raise a Python exception when a creation call returns null."""

  if not value:
    raise StormmError(f"{context} failed: {_last_error_message()}")
  return ctypes.c_void_p(value)


def _check_rc(rc: int, context: str) -> int:
  """Raise a Python exception when a bridge call returns a negative code."""

  if rc < 0:
    raise StormmError(f"{context} failed: {_last_error_message()}")
  return rc


def _encode_optional_path(path: Optional[Union[str, Path]]) -> Optional[bytes]:
  """Normalize Python paths to absolute UTF-8 bytes for C calls."""

  if path is None:
    return None
  return str(Path(path).expanduser().resolve()).encode("utf-8")


_TRAJECTORY_KIND = {"positions": 0, "velocities": 1, "forces": 2}
_THERMOSTAT_KIND = {"none": 0, "andersen": 1, "langevin": 2, "berendsen": 3}


def _normalize_trajectory_kind(kind: str) -> int:
  """Validate and translate user-facing trajectory kind strings."""

  key = kind.lower().strip()
  if key not in _TRAJECTORY_KIND:
    raise ValueError(f"Unsupported trajectory kind '{kind}'.")
  return _TRAJECTORY_KIND[key]


def _normalize_thermostat_kind(kind: Union[str, int, None]) -> int:
  """Validate and translate thermostat selectors."""

  if kind is None:
    return -1
  if isinstance(kind, int):
    if kind < -1 or kind > 3:
      raise ValueError("Thermostat code must be -1, 0, 1, 2, or 3.")
    return kind
  key = kind.lower().strip()
  if key not in _THERMOSTAT_KIND:
    raise ValueError(f"Unsupported thermostat '{kind}'.")
  return _THERMOSTAT_KIND[key]


def state_variable_names() -> List[str]:
  """Return canonical STORMM state variable names in enum order."""

  count = _check_rc(_lib.stormm_state_variable_count(), "state variable count query")
  names: List[str] = []
  for i in range(count):
    raw = _lib.stormm_state_variable_name(i)
    if raw is None:
      raise StormmError(f"state variable name query failed: {_last_error_message()}")
    names.append(raw.decode("utf-8"))
  return names


class AtomGraph:
  """Python owner for a STORMM AtomGraph handle."""

  def __init__(self, topology_path: Union[str, Path]) -> None:
    path = str(Path(topology_path).expanduser().resolve()).encode("utf-8")
    self._handle = _check_ptr(_lib.stormm_atomgraph_create(path), "AtomGraph creation")

  def close(self) -> None:
    if getattr(self, "_handle", None) and self._handle.value:
      _check_rc(_lib.stormm_atomgraph_destroy(self._handle), "AtomGraph destruction")
      self._handle = ctypes.c_void_p()

  def __del__(self) -> None:
    try:
      self.close()
    except Exception:
      pass

  @property
  def atom_count(self) -> int:
    return _check_rc(_lib.stormm_atomgraph_get_atom_count(self._handle), "AtomGraph atom count")

  @property
  def residue_count(self) -> int:
    return _check_rc(_lib.stormm_atomgraph_get_residue_count(self._handle),
                     "AtomGraph residue count")

  @property
  def total_mass(self) -> float:
    mass = _lib.stormm_atomgraph_get_total_mass(self._handle)
    if mass < 0.0:
      raise StormmError(f"AtomGraph total mass query failed: {_last_error_message()}")
    return float(mass)

  def atomic_masses(self) -> List[float]:
    needed = _check_rc(_lib.stormm_atomgraph_get_atomic_masses(self._handle, None, 0),
                       "Atomic mass size query")
    out = (ctypes.c_double * needed)()
    _check_rc(_lib.stormm_atomgraph_get_atomic_masses(self._handle, out, needed),
              "Atomic mass query")
    return list(out)


class PhaseSpace:
  """Python owner for a STORMM PhaseSpace handle."""

  def __init__(self, coordinate_path: Union[str, Path], atomgraph: AtomGraph) -> None:
    # Hold a topology reference so callers cannot accidentally destroy it first.
    self._atomgraph_ref = atomgraph
    path = str(Path(coordinate_path).expanduser().resolve()).encode("utf-8")
    self._handle = _check_ptr(_lib.stormm_phasespace_create(path, atomgraph._handle),
                              "PhaseSpace creation")

  def close(self) -> None:
    if getattr(self, "_handle", None) and self._handle.value:
      _check_rc(_lib.stormm_phasespace_destroy(self._handle), "PhaseSpace destruction")
      self._handle = ctypes.c_void_p()

  def __del__(self) -> None:
    try:
      self.close()
    except Exception:
      pass

  @property
  def atom_count(self) -> int:
    return _check_rc(_lib.stormm_phasespace_get_atom_count(self._handle), "PhaseSpace atom count")

  def interlaced_coordinates(self, kind: str = "positions") -> List[float]:
    kcode = _normalize_trajectory_kind(kind)
    needed = _check_rc(
        _lib.stormm_phasespace_get_interlaced_coordinates(self._handle, kcode, None, 0),
        "PhaseSpace coordinate size query"
    )
    out = (ctypes.c_double * needed)()
    _check_rc(_lib.stormm_phasespace_get_interlaced_coordinates(self._handle, kcode, out, needed),
              "PhaseSpace coordinate query")
    return list(out)

  def coordinates(self, kind: str = "positions") -> List[Tuple[float, float, float]]:
    xyz = self.interlaced_coordinates(kind=kind)
    return [(xyz[i], xyz[i + 1], xyz[i + 2]) for i in range(0, len(xyz), 3)]

  def export(self, output_file: Union[str, Path], current_time: float = 0.0) -> None:
    out_path = str(Path(output_file).expanduser().resolve()).encode("utf-8")
    _check_rc(_lib.stormm_phasespace_export(self._handle, out_path, float(current_time)),
              "PhaseSpace export")


class DynamicsControls:
  """Python owner for STORMM DynamicsControls."""

  def __init__(self) -> None:
    self._handle = _check_ptr(_lib.stormm_dynamics_controls_create(), "DynamicsControls creation")

  def close(self) -> None:
    if getattr(self, "_handle", None) and self._handle.value:
      _check_rc(_lib.stormm_dynamics_controls_destroy(self._handle),
                "DynamicsControls destruction")
      self._handle = ctypes.c_void_p()

  def __del__(self) -> None:
    try:
      self.close()
    except Exception:
      pass

  @property
  def step_count(self) -> int:
    return _check_rc(_lib.stormm_dynamics_controls_get_step_count(self._handle),
                     "DynamicsControls step count")

  def set_step_count(self, step_count: int) -> "DynamicsControls":
    _check_rc(_lib.stormm_dynamics_controls_set_step_count(self._handle, int(step_count)),
              "set step count")
    return self

  def set_time_step(self, time_step_fs: float) -> "DynamicsControls":
    _check_rc(_lib.stormm_dynamics_controls_set_time_step(self._handle, float(time_step_fs)),
              "set time step")
    return self

  def set_diagnostic_print_frequency(self, frequency: int) -> "DynamicsControls":
    _check_rc(_lib.stormm_dynamics_controls_set_diagnostic_frequency(self._handle, int(frequency)),
              "set diagnostic print frequency")
    return self

  def set_trajectory_print_frequency(self, frequency: int) -> "DynamicsControls":
    _check_rc(_lib.stormm_dynamics_controls_set_trajectory_frequency(self._handle, int(frequency)),
              "set trajectory print frequency")
    return self

  def set_center_of_mass_purge_frequency(self, frequency: int) -> "DynamicsControls":
    _check_rc(_lib.stormm_dynamics_controls_set_com_purge_frequency(self._handle, int(frequency)),
              "set center-of-mass purge frequency")
    return self

  def set_thermostat_kind(self, kind: Union[str, int]) -> "DynamicsControls":
    code = _normalize_thermostat_kind(kind)
    if code < 0:
      raise ValueError("DynamicsControls thermostat kind cannot be None.")
    _check_rc(_lib.stormm_dynamics_controls_set_thermostat_kind(self._handle, code),
              "set thermostat kind")
    return self

  def set_temperature(self, initial_kelvin: float,
                      final_kelvin: Optional[float] = None) -> "DynamicsControls":
    if final_kelvin is None:
      final_kelvin = initial_kelvin
    _check_rc(_lib.stormm_dynamics_controls_set_temperature(self._handle, float(initial_kelvin),
                                                            float(final_kelvin)),
              "set thermostat temperatures")
    return self

  def set_langevin_frequency(self, collision_frequency: float) -> "DynamicsControls":
    _check_rc(_lib.stormm_dynamics_controls_set_langevin_frequency(self._handle,
                                                                   float(collision_frequency)),
              "set Langevin frequency")
    return self

  def set_andersen_frequency(self, reassignment_frequency: int) -> "DynamicsControls":
    _check_rc(_lib.stormm_dynamics_controls_set_andersen_frequency(self._handle,
                                                                   int(reassignment_frequency)),
              "set Andersen frequency")
    return self

  def set_seed(self, random_seed: int) -> "DynamicsControls":
    _check_rc(_lib.stormm_dynamics_controls_set_seed(self._handle, int(random_seed)),
              "set thermostat seed")
    return self


class AtomGraphSynthesis:
  """Python owner for STORMM AtomGraphSynthesis."""

  def __init__(self, atomgraphs: Sequence[AtomGraph]) -> None:
    if not atomgraphs:
      raise ValueError("AtomGraphSynthesis needs at least one AtomGraph.")
    self._atomgraphs = list(atomgraphs)
    handles = (ctypes.c_void_p * len(atomgraphs))(*[ag._handle.value for ag in atomgraphs])
    self._handle = _check_ptr(_lib.stormm_atomgraph_synthesis_create(handles, len(atomgraphs)),
                              "AtomGraphSynthesis creation")

  def close(self) -> None:
    if getattr(self, "_handle", None) and self._handle.value:
      _check_rc(_lib.stormm_atomgraph_synthesis_destroy(self._handle),
                "AtomGraphSynthesis destruction")
      self._handle = ctypes.c_void_p()

  def __del__(self) -> None:
    try:
      self.close()
    except Exception:
      pass

  @property
  def system_count(self) -> int:
    return _check_rc(_lib.stormm_atomgraph_synthesis_get_system_count(self._handle),
                     "AtomGraphSynthesis system count")

  @property
  def unique_topology_count(self) -> int:
    return _check_rc(_lib.stormm_atomgraph_synthesis_get_unique_topology_count(self._handle),
                     "AtomGraphSynthesis unique topology count")

  @property
  def atom_count(self) -> int:
    return _check_rc(_lib.stormm_atomgraph_synthesis_get_atom_count(self._handle),
                     "AtomGraphSynthesis atom count")

  def atom_count_of(self, system_index: int) -> int:
    return _check_rc(_lib.stormm_atomgraph_synthesis_get_atom_count_system(self._handle,
                                                                            int(system_index)),
                     "AtomGraphSynthesis system atom count")


class PhaseSpaceSynthesis:
  """Python owner for STORMM PhaseSpaceSynthesis."""

  def __init__(self, phasespaces: Sequence[PhaseSpace], atomgraphs: Sequence[AtomGraph]) -> None:
    if not phasespaces:
      raise ValueError("PhaseSpaceSynthesis needs at least one system.")
    if len(phasespaces) != len(atomgraphs):
      raise ValueError("PhaseSpace and AtomGraph sequences must have identical lengths.")

    self._phasespaces = list(phasespaces)
    self._atomgraphs = list(atomgraphs)
    ps_handles = (ctypes.c_void_p * len(phasespaces))(*[ps._handle.value for ps in phasespaces])
    ag_handles = (ctypes.c_void_p * len(atomgraphs))(*[ag._handle.value for ag in atomgraphs])
    self._handle = _check_ptr(_lib.stormm_phasespace_synthesis_create(ps_handles, ag_handles,
                                                                       len(phasespaces)),
                              "PhaseSpaceSynthesis creation")

  def close(self) -> None:
    if getattr(self, "_handle", None) and self._handle.value:
      _check_rc(_lib.stormm_phasespace_synthesis_destroy(self._handle),
                "PhaseSpaceSynthesis destruction")
      self._handle = ctypes.c_void_p()

  def __del__(self) -> None:
    try:
      self.close()
    except Exception:
      pass

  @property
  def system_count(self) -> int:
    return _check_rc(_lib.stormm_phasespace_synthesis_get_system_count(self._handle),
                     "PhaseSpaceSynthesis system count")

  @property
  def unique_topology_count(self) -> int:
    return _check_rc(_lib.stormm_phasespace_synthesis_get_unique_topology_count(self._handle),
                     "PhaseSpaceSynthesis unique topology count")

  @property
  def padded_atom_count(self) -> int:
    return _check_rc(_lib.stormm_phasespace_synthesis_get_padded_atom_count(self._handle),
                     "PhaseSpaceSynthesis padded atom count")

  def atom_count_of(self, system_index: int) -> int:
    return _check_rc(_lib.stormm_phasespace_synthesis_get_atom_count(self._handle,
                                                                     int(system_index)),
                     "PhaseSpaceSynthesis system atom count")

  def interlaced_coordinates(self, system_index: int, kind: str = "positions") -> List[float]:
    kcode = _normalize_trajectory_kind(kind)
    needed = _check_rc(_lib.stormm_phasespace_synthesis_get_interlaced_coordinates(
        self._handle, int(system_index), kcode, None, 0), "PhaseSpaceSynthesis coordinate size query")
    out = (ctypes.c_double * needed)()
    _check_rc(_lib.stormm_phasespace_synthesis_get_interlaced_coordinates(
        self._handle, int(system_index), kcode, out, needed), "PhaseSpaceSynthesis coordinate query")
    return list(out)

  def coordinates(self, system_index: int, kind: str = "positions") -> List[Tuple[float, float, float]]:
    xyz = self.interlaced_coordinates(system_index, kind=kind)
    return [(xyz[i], xyz[i + 1], xyz[i + 2]) for i in range(0, len(xyz), 3)]


@dataclass(frozen=True)
class DynamicsResult:
  """Structured return values from `run_dynamics`."""

  potential_energy: float
  total_energy: float
  sample_count: int
  instantaneous_states: Dict[str, float]


def run_dynamics(phasespace: PhaseSpace, atomgraph: AtomGraph, controls: Optional[DynamicsControls] = None,
                 thermostat: Union[str, int, None] = None, temperature_kelvin: float = 298.15,
                 trajectory_file: Optional[Union[str, Path]] = None,
                 restart_file: Optional[Union[str, Path]] = None) -> DynamicsResult:
  """Run STORMM dynamics and return energies/states in Python-native types."""

  state_count = _check_rc(_lib.stormm_state_variable_count(), "state variable count query")
  state_values = (ctypes.c_double * state_count)()
  potential = ctypes.c_double(0.0)
  total = ctypes.c_double(0.0)
  sample_count = ctypes.c_int(0)
  thermostat_code = _normalize_thermostat_kind(thermostat)
  trajectory_bytes = _encode_optional_path(trajectory_file)
  restart_bytes = _encode_optional_path(restart_file)
  controls_handle = controls._handle if controls is not None else None

  rc = _lib.stormm_dynamics_run(phasespace._handle, atomgraph._handle, controls_handle,
                                thermostat_code, float(temperature_kelvin),
                                trajectory_bytes, restart_bytes,
                                ctypes.byref(potential), ctypes.byref(total),
                                ctypes.byref(sample_count), state_values, state_count)
  _check_rc(rc, "Dynamics execution")

  names = state_variable_names()
  states = {name: state_values[i] for i, name in enumerate(names)}
  return DynamicsResult(potential_energy=float(potential.value),
                        total_energy=float(total.value),
                        sample_count=int(sample_count.value),
                        instantaneous_states=states)


__all__ = [
    "StormmError",
    "AtomGraph",
    "PhaseSpace",
    "AtomGraphSynthesis",
    "PhaseSpaceSynthesis",
    "DynamicsControls",
    "DynamicsResult",
    "run_dynamics",
    "state_variable_names",
]
