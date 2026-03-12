#!/usr/bin/env python3
"""End-to-end smoke test for the Python STORMM dynamics bridge.

The script demonstrates the intended user workflow:
1. Create AtomGraph and PhaseSpace objects from topology/coordinate files.
2. Build AtomGraphSynthesis and PhaseSpaceSynthesis from those objects.
3. Configure DynamicsControls and run dynamics.
4. Read all results back in Python for post-processing.
"""

from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
# Allow this script to run directly from the repository without package installation.
sys.path.insert(0, str(REPO_ROOT / "python"))

from stormm import dynamics


def center_of_mass(coords, masses):
  """Compute a Cartesian center of mass in pure Python."""

  total_mass = sum(masses)
  mx = sum(masses[i] * coords[i][0] for i in range(len(coords))) / total_mass
  my = sum(masses[i] * coords[i][1] for i in range(len(coords))) / total_mass
  mz = sum(masses[i] * coords[i][2] for i in range(len(coords))) / total_mass
  return (mx, my, mz)


def main():
  # Pick two compatible systems (both isolated-boundary) to exercise synthesis constructors.
  top_a = REPO_ROOT / "test" / "Topology" / "ala_dipeptide.top"
  crd_a = REPO_ROOT / "test" / "Trajectory" / "ala_dipeptide.inpcrd"
  top_b = REPO_ROOT / "test" / "Topology" / "trpcage.top"
  crd_b = REPO_ROOT / "test" / "Trajectory" / "trpcage.inpcrd"

  # Construct primary single-system objects from STORMM input files.
  ag_a = dynamics.AtomGraph(top_a)
  ps_a = dynamics.PhaseSpace(crd_a, ag_a)
  ag_b = dynamics.AtomGraph(top_b)
  ps_b = dynamics.PhaseSpace(crd_b, ag_b)

  # Construct synthesis objects from existing C++ objects (no manual hardcoding of arrays).
  ag_syn = dynamics.AtomGraphSynthesis([ag_a, ag_b])
  ps_syn = dynamics.PhaseSpaceSynthesis([ps_a, ps_b], [ag_a, ag_b])

  print(f"System A atoms: {ag_a.atom_count}, residues: {ag_a.residue_count}")
  print(f"System B atoms: {ag_b.atom_count}, residues: {ag_b.residue_count}")
  print(f"AtomGraphSynthesis systems: {ag_syn.system_count}, unique topologies: "
        f"{ag_syn.unique_topology_count}, total atoms: {ag_syn.atom_count}")
  print(f"PhaseSpaceSynthesis systems: {ps_syn.system_count}, unique topologies: "
        f"{ps_syn.unique_topology_count}, padded atoms: {ps_syn.padded_atom_count}")

  # Pull coordinates and masses into Python and perform an analysis calculation.
  masses_a = ag_a.atomic_masses()
  coords_a_before = ps_a.coordinates("positions")
  com_before = center_of_mass(coords_a_before, masses_a)
  print(f"COM before dynamics: {com_before}")

  # Configure integration controls in Python through direct wrappers on DynamicsControls setters.
  controls = dynamics.DynamicsControls()
  controls.set_step_count(20)
  controls.set_time_step(1.0)
  controls.set_diagnostic_print_frequency(5)
  controls.set_trajectory_print_frequency(10)
  controls.set_center_of_mass_purge_frequency(1000)
  controls.set_thermostat_kind("langevin")
  controls.set_temperature(300.0)
  controls.set_langevin_frequency(0.001)
  controls.set_seed(20260311)

  # Run CPU dynamics and request trajectory/restart output files for inspection.
  traj_out = Path("/tmp/stormm_python_demo.crd")
  rst_out = Path("/tmp/stormm_python_demo.inpcrd")
  result = dynamics.run_dynamics(
      ps_a,
      ag_a,
      controls=controls,
      thermostat="langevin",
      temperature_kelvin=300.0,
      trajectory_file=traj_out,
      restart_file=rst_out,
  )

  # Re-read coordinates from the updated PhaseSpace handle and compare with initial values.
  coords_a_after = ps_a.coordinates("positions")
  com_after = center_of_mass(coords_a_after, masses_a)
  print(f"COM after dynamics:  {com_after}")

  first_a_before = coords_a_before[0]
  first_a_after = coords_a_after[0]
  print("First atom displacement: "
        f"({first_a_after[0] - first_a_before[0]:.6f}, "
        f"{first_a_after[1] - first_a_before[1]:.6f}, "
        f"{first_a_after[2] - first_a_before[2]:.6f})")

  print(f"Potential energy: {result.potential_energy:.6f} kcal/mol")
  print(f"Total energy:     {result.total_energy:.6f} kcal/mol")
  print(f"Score samples:    {result.sample_count}")
  # Inspect selected state variables returned from the STORMM ScoreCard.
  print("Selected states:")
  for key in ["BOND", "ANGLE", "VDW", "ELECTROSTATIC", "KINETIC", "POTENTIAL_ENERGY",
              "TOTAL_ENERGY"]:
    print(f"  {key:>18s}: {result.instantaneous_states[key]: .6f}")

  # Access coordinates from PhaseSpaceSynthesis to verify synthesis data is queryable in Python.
  syn_crd_b = ps_syn.coordinates(1, "positions")
  print(f"Synthesis system 1 first atom: {syn_crd_b[0]}")
  print(f"Trajectory output: {traj_out}")
  print(f"Restart output:    {rst_out}")


if __name__ == "__main__":
  main()
