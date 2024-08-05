// -*-c++-*-
#include <cuda.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <nvml.h>
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/core_kernel_manager.h"
#include "../../src/Constants/scaling.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/rounding.h"
#include "../../src/MolecularMechanics/mm_evaluation.h"
#include "../../src/Namelists/nml_files.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Potential/cacheresource.h"
#include "../../src/Potential/hpc_valence_potential.h"
#include "../../src/Potential/hpc_nonbonded_potential.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Restraints/restraint_apparatus.h"
#include "../../src/Structure/hpc_virtual_site_handling.h"
#include "../../src/Structure/structure_enumerators.h"
#include "../../src/Structure/virtual_site_handling.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/implicit_solvent_workspace.h"
#include "../../src/Synthesis/synthesis_abstracts.h"
#include "../../src/Synthesis/nonbonded_workunit.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/systemcache.h"
#include "../../src/Synthesis/valence_workunit.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/thermostat.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "assemble_restraints.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::errors;
using namespace stormm::diskutil;
using namespace stormm::energy;
using namespace stormm::stmath;
using namespace stormm::mm;
using namespace stormm::numerics;
using namespace stormm::parse;
using namespace stormm::random;
using namespace stormm::restraints;
using namespace stormm::review;
using namespace stormm::structure;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Compute forces due to valence interactions acting on a series of systems using topology and
// coordinate compilations.  Check the results against the accumulations for individual systems.
// This function assumes that the topology and coordinate compilations have already been uploaded
// to the device.
//
// Arguments:
//   poly_ps:           Coordinates for many systems
//   mmctrl:            Molecular mechanics progress counters
//   valence_tb_space:  Thread block resources, pre-allocated on the GPU
//   nonbond_tb_space:  Thread block resources, pre-allocated on the GPU
//   poly_ag:           Topologies for many systems
//   facc_method:       Force accumulation method
//   prec:              Precision level at which to perform the calculations (may not be
//                      compatible with all force accumulation methods)
//   gpu:               Details of the GPU to use
//   mue_tol:           Tolerance for mean unsigned error in forces 
//   max_error_tol:     Tolerance for maximum unsigned error in forces 
//-------------------------------------------------------------------------------------------------
void checkCompilationForces(PhaseSpaceSynthesis *poly_ps, MolecularMechanicsControls *mmctrl,
                            CacheResource *valence_tb_space, CacheResource *nonbond_tb_space,
                            const AtomGraphSynthesis &poly_ag,
                            const StaticExclusionMaskSynthesis &poly_se,
                            const AccumulationMethod facc_method, const PrecisionModel prec,
                            const GpuDetails &gpu, const CoreKlManager &launcher,
                            const double mue_tol, const double max_error_tol,
                            const TestPriority do_tests,
                            const std::string &side_note = std::string(""),
                            const bool do_valence_tests = true) {

  // Prepare for GPU-based calculations
  const int nsys = poly_ps->getSystemCount();
  std::vector<double> frc_mues(nsys);
  const Approx frc_mue_tolerance = Approx(std::vector<double>(nsys, 0.0)).margin(mue_tol);
  const Approx frc_max_error_tolerance =
    Approx(std::vector<double>(nsys, 0.0)).margin(max_error_tol);
  std::vector<double> frc_max_errors(nsys);
  ScoreCard sc(nsys, 1, 32);
  const TrajectoryKind frcid = TrajectoryKind::FORCES;
  ImplicitSolventWorkspace ism_space(poly_ag.getSystemAtomOffsets(),
                                     poly_ag.getSystemAtomCounts(), prec);
  Thermostat heat_bath(poly_ag, ThermostatKind::NONE);

  // Clear forces.  Perform calculations for valence interactions and restraints.
  poly_ps->initializeForces(gpu, HybridTargetLevel::DEVICE);
  mmctrl->incrementStep();
  launchValence(prec, poly_ag, mmctrl, poly_ps, &heat_bath, &sc, valence_tb_space,
                EvaluateForce::YES, EvaluateEnergy::NO, VwuGoal::ACCUMULATE, facc_method,
                launcher);
  int total_restraints = 0;
  for (int i = 0; i < nsys; i++) {
    PhaseSpace host_result = poly_ps->exportSystem(i, HybridTargetLevel::HOST);
    PhaseSpace devc_result = poly_ps->exportSystem(i, HybridTargetLevel::DEVICE);
    host_result.initializeForces();
    ScoreCard isc(1, 1, 32);
    evalValeRestMM(&host_result, &isc, poly_ag.getSystemTopologyPointer(i),
                   *(poly_ag.getSystemRestraintPointer(i)), EvaluateForce::YES, 0);
    const std::vector<double> devc_frc = devc_result.getInterlacedCoordinates(frcid);
    const std::vector<double> host_frc = host_result.getInterlacedCoordinates(frcid);
    frc_mues[i] = meanUnsignedError(devc_frc, host_frc);
    frc_max_errors[i] = maxAbsoluteDifference(devc_frc, host_frc);
    total_restraints += poly_ag.getSystemRestraintPointer(i)->getTotalRestraintCount();
  }
  const std::string restraint_presence = (total_restraints > 0) ? "with" : "without";
  const std::string end_note = (side_note.size() > 0LLU) ? "  " + side_note : std::string("");
  if (do_valence_tests) {
    check(frc_mues, RelationalOperator::EQUAL, frc_mue_tolerance, "Forces obtained by the "
          "valence interaction kernel, operating on systems " + restraint_presence + " external " +
          "restraints, exceed the tolerance for mean unsigned errors in their vector components.  "
          "Force accumulation method: " + getAccumulationMethodName(facc_method) +
          ".  Precision level in the calculation: " + getEnumerationName(prec) + "." + end_note,
          do_tests);
    check(frc_max_errors, RelationalOperator::EQUAL, frc_max_error_tolerance, "Forces obtained "
          "by the valence interaction kernel, operating on systems " + restraint_presence +
          " external restraints, exceed the maximum allowed errors for forces acting on any one "
          "particle.  Force accumulation method: " + getAccumulationMethodName(facc_method) +
          ".  Precision level in the calculation: " + getEnumerationName(prec) + "." + end_note,
          do_tests);
  }

  // Clear forces again.  Compute non-bonded interactions.
  poly_ps->initializeForces(gpu, HybridTargetLevel::DEVICE);
  mmctrl->incrementStep();
  launchNonbonded(prec, poly_ag, poly_se, mmctrl, poly_ps, &heat_bath, &sc, nonbond_tb_space,
                  &ism_space, EvaluateForce::YES, EvaluateEnergy::NO, facc_method, launcher);
  if (poly_ag.getImplicitSolventModel() != ImplicitSolventModel::NONE) {
    ism_space.download();
  }
  const ISWorkspaceKit<double> ismr = ism_space.dpData();
  const NeckGeneralizedBornTable ngb_tables;
  const int total_atoms = sum<int>(poly_ag.getSystemAtomCounts());
  std::vector<double> cpu_psi_values(total_atoms), gpu_psi_values(total_atoms);
  int gb_counter = 0;
  for (int i = 0; i < nsys; i++) {
    PhaseSpace host_result = poly_ps->exportSystem(i, HybridTargetLevel::HOST);
    PhaseSpace devc_result = poly_ps->exportSystem(i, HybridTargetLevel::DEVICE);
    host_result.initializeForces();
    ScoreCard isc(1, 1, 32);
    const AtomGraph *iag_ptr = poly_ag.getSystemTopologyPointer(i);
    const StaticExclusionMask ise(iag_ptr);
    const double2 tnbe = evaluateNonbondedEnergy(iag_ptr, ise, &host_result, &isc,
                                                 EvaluateForce::YES, EvaluateForce::YES, 0);
    if (poly_ag.getImplicitSolventModel() != ImplicitSolventModel::NONE) {
      PhaseSpaceWriter hostw = host_result.data();
      std::vector<double> eff_gb_radii(hostw.natom, 0.0);
      std::vector<double> psi(hostw.natom, 0.0);
      std::vector<double> sum_deijda(hostw.natom, 0.0);
      const NonbondedKit<double> inbk = iag_ptr->getDoublePrecisionNonbondedKit();
      const ValenceKit<double> ivk = iag_ptr->getDoublePrecisionValenceKit();
      const ImplicitSolventKit<double> iisk = iag_ptr->getDoublePrecisionImplicitSolventKit();
      const ImplicitSolventRecipe<double> isr(iisk, ngb_tables.dpData());
      const double tgbe = evaluateGeneralizedBornEnergy(inbk, ise.data(), iisk,
                                                        ngb_tables.dpData(), hostw.xcrd,
                                                        hostw.ycrd, hostw.zcrd, hostw.xfrc,
                                                        hostw.yfrc, hostw.zfrc,
                                                        eff_gb_radii.data(), psi.data(),
                                                        sum_deijda.data(), &isc,
                                                        EvaluateForce::YES, 0);
      const int sys_ofs = poly_ag.getSystemAtomOffsets().readHost(i);
      for (int j = 0; j < hostw.natom; j++) {
        cpu_psi_values[gb_counter] = psi[j];
        switch (prec) {
        case PrecisionModel::DOUBLE:
          gpu_psi_values[gb_counter] = hostInt95ToDouble(ismr.psi[sys_ofs + j],
                                                         ismr.psi_ovrf[sys_ofs + j]) *
                                       ismr.inv_fp_scale;
          break;
        case PrecisionModel::SINGLE:
          gpu_psi_values[gb_counter] = static_cast<double>(ismr.psi[sys_ofs + j]) *
                                       ismr.inv_fp_scale;
          break;
        }
        gb_counter++;
      }
    }
    std::vector<double> devc_frc = devc_result.getInterlacedCoordinates(frcid);
    std::vector<double> host_frc = host_result.getInterlacedCoordinates(frcid);
    
    // These systems contain some hard clashes, which generate very large forces.  This is good
    // for testing the split force accumulation, but not for discerning values that are truly
    // inaccurate.  Check the forces individually and clean out large values that are within
    // relative tolerances.
    const int natom = iag_ptr->getAtomCount();
    for (int j = 0; j < 3 * natom; j++) {
      if ((fabs(devc_frc[j]) >= 200.0 || fabs(host_frc[j]) >= 200.0) &&
          1.0 - (host_frc[j] / devc_frc[j]) <= max_error_tol) {
        devc_frc[j] = 0.0;
        host_frc[j] = 0.0;
      }
    }
    frc_mues[i] = meanUnsignedError(devc_frc, host_frc);
    frc_max_errors[i] = maxAbsoluteDifference(devc_frc, host_frc);
  }
  if (poly_ag.getImplicitSolventModel() != ImplicitSolventModel::NONE) {
    check(gpu_psi_values, RelationalOperator::EQUAL, Approx(cpu_psi_values).margin(max_error_tol),
          "Values for Psi computed on the GPU do not agree with their CPU counterparts.  "
          "Accumulation method: " + getAccumulationMethodName(facc_method) + ".  Precision level "
          "in the calculation: " + getEnumerationName(prec) + "." + end_note, do_tests);
  }
  check(frc_mues, RelationalOperator::EQUAL, frc_mue_tolerance, "Forces obtained by the "
        "non-bonded interaction kernel, operating on systems " + restraint_presence +
        " external restraints, exceed the tolerance for mean unsigned errors in their vector "
        "components.  Force accumulation method: " + getAccumulationMethodName(facc_method) +
        ".  Precision level in the calculation: " + getEnumerationName(prec) + "." + end_note,
        do_tests);
  check(frc_max_errors, RelationalOperator::EQUAL, frc_max_error_tolerance, "Forces obtained "
        "by the non-bonded interaction kernel, operating on systems " + restraint_presence +
        " external restraints, exceed the maximum allowed errors for forces acting on any one "
        "particle.  Force accumulation method: " + getAccumulationMethodName(facc_method) +
        ".  Precision level in the calculation: " + getEnumerationName(prec) + "." + end_note,
        do_tests);
}

//-------------------------------------------------------------------------------------------------
// Compute energies of a series of systems, using topology and coordinate compilations on the GPU,
// due to valence interactions.  Check the results against the accumulations for individual
// systems.  This function assumes that the topology and coordinate compilations have already been
// uploaded to the device.
//
// Arguments:
//   poly_ps:           Coordinates of all systems
//   mmctrl:            Molecular mechanics control data and progress counters
//   valence_tb_space:  Thread-block specific L1 scratch space allocations (valence work)
//   nonbond_tb_space:  Thread-block specific L1 scratch space allocations (non-bonded work)
//   poly_ag:           Collated topologies and parameters for all systems
//   poly_se:           Static exclusion masks for all systems
//   prec:              The precision level at which to operate
//   gpu:               Details of the GPU chosen for calculations
//   launcher:          Repository of kernel launch parameters
//   bond_tol:          Tolerance for bond energy calculations
//   angl_tol:          Tolerance for angle energy calculations
//   dihe_tol:          Tolerance for dihedral energy calculations
//   impr_tol:          Tolerance for CHARMM improper dihedral energy calculations
//   ubrd_tol:          Tolerance for Urey-Bradley energy calculations
//   cimp_tol:          Tolerance for CHARMM improper energy calculations
//   cmap_tol:          Tolerance for CMAP energy calculations
//   lj14_tol:          Tolerance for Lennard-Jones 1:4 energy calculations
//   qq14_tol:          Tolerance for electrostatic 1:4 energy calculations
//   rstr_tol:          Tolerance for restraint energy calculations
//   do_tests:          Indication that tests should be performed or aborted
//-------------------------------------------------------------------------------------------------
void checkCompilationEnergies(PhaseSpaceSynthesis *poly_ps, MolecularMechanicsControls *mmctrl,
                              CacheResource *valence_tb_space, CacheResource *nonbond_tb_space,
                              const AtomGraphSynthesis &poly_ag,
                              const StaticExclusionMaskSynthesis &poly_se,
                              const PrecisionModel prec, const GpuDetails &gpu,
                              const CoreKlManager &launcher, const double bond_tol,
                              const double angl_tol, const double dihe_tol, const double impr_tol,
                              const double ubrd_tol, const double cimp_tol, const double cmap_tol,
                              const double lj14_tol, const double qq14_tol, const double rstr_tol,
                              const double ljnb_tol, const double qqnb_tol, const double gbnb_tol,
                              const TestPriority do_tests, const bool do_valence_tests = true) {
  const int nsys = poly_ps->getSystemCount();
  ScoreCard sc(nsys, 1, 32);
  Thermostat heat_bath(poly_ag, ThermostatKind::NONE);
  poly_ps->initializeForces(gpu, HybridTargetLevel::DEVICE);
  mmctrl->incrementStep();
  launchValence(prec, poly_ag, mmctrl, poly_ps, &heat_bath, &sc, valence_tb_space,
                EvaluateForce::NO, EvaluateEnergy::YES, VwuGoal::ACCUMULATE,
                AccumulationMethod::SPLIT, launcher);
  ImplicitSolventWorkspace ism_space(poly_ag.getSystemAtomOffsets(),
                                     poly_ag.getSystemAtomCounts(), prec);
  const NeckGeneralizedBornTable ngb_tables;
  launchNonbonded(prec, poly_ag, poly_se, mmctrl, poly_ps, &heat_bath, &sc, nonbond_tb_space,
                  &ism_space, EvaluateForce::NO, EvaluateEnergy::YES, AccumulationMethod::SPLIT,
                  launcher);
  sc.download();
  std::vector<double> cpu_bond(nsys), gpu_bond(nsys), cpu_angl(nsys), gpu_angl(nsys);
  std::vector<double> cpu_dihe(nsys), gpu_dihe(nsys), cpu_impr(nsys), gpu_impr(nsys);
  std::vector<double> cpu_ubrd(nsys), gpu_ubrd(nsys), cpu_cimp(nsys), gpu_cimp(nsys);
  std::vector<double> cpu_cmap(nsys), gpu_cmap(nsys), cpu_qq14(nsys), gpu_qq14(nsys);
  std::vector<double> cpu_lj14(nsys), gpu_lj14(nsys), cpu_rstr(nsys), gpu_rstr(nsys);
  std::vector<double> cpu_ljnb(nsys), gpu_ljnb(nsys), cpu_qqnb(nsys), gpu_qqnb(nsys);
  std::vector<double> cpu_gbnb(nsys), gpu_gbnb(nsys);
  int nrstr = 0;
  for (int i = 0; i < nsys; i++) {
    PhaseSpace devc_result = poly_ps->exportSystem(i, HybridTargetLevel::DEVICE);
    PhaseSpace host_result = poly_ps->exportSystem(i, HybridTargetLevel::HOST);
    host_result.initializeForces();
    ScoreCard isc(1, 1, 32);
    const StaticExclusionMask sysi_se(poly_ag.getSystemTopologyPointer(i));
    if (poly_ag.getImplicitSolventModel() != ImplicitSolventModel::NONE) {
      evalRestrainedMMGB(&host_result, &isc, poly_ag.getSystemTopologyPointer(i),
                         ngb_tables, sysi_se, poly_ag.getSystemRestraintPointer(i),
                         EvaluateForce::NO, 0, 0);
    }
    else {
      evalNonbValeRestMM(&host_result, &isc, poly_ag.getSystemTopologyPointer(i), sysi_se,
                         poly_ag.getSystemRestraintPointer(i), EvaluateForce::NO, 0, 0);
    }
    nrstr += poly_ag.getSystemRestraintPointer(i)->getTotalRestraintCount();
    gpu_bond[i] =  sc.reportInstantaneousStates(StateVariable::BOND, i);
    cpu_bond[i] = isc.reportInstantaneousStates(StateVariable::BOND, 0);
    gpu_angl[i] =  sc.reportInstantaneousStates(StateVariable::ANGLE, i);
    cpu_angl[i] = isc.reportInstantaneousStates(StateVariable::ANGLE, 0);
    gpu_dihe[i] =  sc.reportInstantaneousStates(StateVariable::PROPER_DIHEDRAL, i);
    cpu_dihe[i] = isc.reportInstantaneousStates(StateVariable::PROPER_DIHEDRAL, 0);
    gpu_impr[i] =  sc.reportInstantaneousStates(StateVariable::IMPROPER_DIHEDRAL, i);
    cpu_impr[i] = isc.reportInstantaneousStates(StateVariable::IMPROPER_DIHEDRAL, 0);
    gpu_ubrd[i] =  sc.reportInstantaneousStates(StateVariable::UREY_BRADLEY, i);
    cpu_ubrd[i] = isc.reportInstantaneousStates(StateVariable::UREY_BRADLEY, 0);
    gpu_cimp[i] =  sc.reportInstantaneousStates(StateVariable::CHARMM_IMPROPER, i);
    cpu_cimp[i] = isc.reportInstantaneousStates(StateVariable::CHARMM_IMPROPER, 0);
    gpu_cmap[i] =  sc.reportInstantaneousStates(StateVariable::CMAP, i);
    cpu_cmap[i] = isc.reportInstantaneousStates(StateVariable::CMAP, 0);
    gpu_qq14[i] =  sc.reportInstantaneousStates(StateVariable::ELEC_ONE_FOUR, i);
    cpu_qq14[i] = isc.reportInstantaneousStates(StateVariable::ELEC_ONE_FOUR, 0);
    gpu_lj14[i] =  sc.reportInstantaneousStates(StateVariable::VDW_ONE_FOUR, i);
    cpu_lj14[i] = isc.reportInstantaneousStates(StateVariable::VDW_ONE_FOUR, 0);
    gpu_rstr[i] =  sc.reportInstantaneousStates(StateVariable::RESTRAINT, i);
    cpu_rstr[i] = isc.reportInstantaneousStates(StateVariable::RESTRAINT, 0);
    gpu_ljnb[i] =  sc.reportInstantaneousStates(StateVariable::VDW, i);
    cpu_ljnb[i] = isc.reportInstantaneousStates(StateVariable::VDW, 0);
    gpu_qqnb[i] =  sc.reportInstantaneousStates(StateVariable::ELECTROSTATIC, i);
    cpu_qqnb[i] = isc.reportInstantaneousStates(StateVariable::ELECTROSTATIC, 0);
    if (poly_ag.getImplicitSolventModel() != ImplicitSolventModel::NONE) {
      gpu_gbnb[i] =  sc.reportInstantaneousStates(StateVariable::ELECTROSTATIC, i);
      cpu_gbnb[i] = isc.reportInstantaneousStates(StateVariable::ELECTROSTATIC, 0);
    }
  }
  if (do_valence_tests) {
    check(gpu_bond, RelationalOperator::EQUAL, Approx(cpu_bond).margin(bond_tol), "Bond energies "
          "computed on the CPU and GPU do not agree.  Precision level in the calculation: " +
          getEnumerationName(prec) + ".", do_tests);
    check(gpu_angl, RelationalOperator::EQUAL, Approx(cpu_angl).margin(angl_tol), "Angle energies "
          "computed on the CPU and GPU do not agree.  Precision level in the calculation: " +
          getEnumerationName(prec) + ".", do_tests);
    check(gpu_dihe, RelationalOperator::EQUAL, Approx(cpu_dihe).margin(dihe_tol), "Proper "
          "dihedral energies computed on the CPU and GPU do not agree.  Precision level in the "
          "calculation: " + getEnumerationName(prec) + ".", do_tests);
    check(gpu_impr, RelationalOperator::EQUAL, Approx(cpu_impr).margin(impr_tol), "Improper "
          "dihedral energies computed on the CPU and GPU do not agree.  Precision level in the "
          "calculation: " + getEnumerationName(prec) + ".", do_tests);
    check(gpu_ubrd, RelationalOperator::EQUAL, Approx(cpu_ubrd).margin(ubrd_tol), "Urey-Bradley "
          "energies computed on the CPU and GPU do not agree.  Precision level in the "
          "calculation: " + getEnumerationName(prec) + ".", do_tests);
    check(gpu_cimp, RelationalOperator::EQUAL, Approx(cpu_cimp).margin(cimp_tol), "CHARMM "
          "improper dihedral energies computed on the CPU and GPU do not agree.  Precision level "
          "in the calculation: " + getEnumerationName(prec) + ".", do_tests);
    check(gpu_cmap, RelationalOperator::EQUAL, Approx(cpu_cmap).margin(cmap_tol), "CMAP "
          "energies computed on the CPU and GPU do not agree.  Precision level in the "
          "calculation: " + getEnumerationName(prec) + ".", do_tests);
    check(gpu_qq14, RelationalOperator::EQUAL, Approx(cpu_qq14).margin(qq14_tol), "Electrostatic "
          "1:4 energies computed on the CPU and GPU do not agree.  Precision level in the "
          "calculation: " + getEnumerationName(prec) + ".", do_tests);
    check(gpu_lj14, RelationalOperator::EQUAL, Approx(cpu_lj14).margin(lj14_tol), "Lennard-Jones "
          "1:4 energies computed on the CPU and GPU do not agree.  Precision level in the "
          "calculation: " + getEnumerationName(prec) + ".", do_tests);
    if (nrstr > 0) {
      check(gpu_rstr, RelationalOperator::EQUAL, Approx(cpu_rstr).margin(rstr_tol), "Restraint "
            "energies computed on the CPU and GPU do not agree.  Precision level in the "
            "calculation: " + getEnumerationName(prec) + ".", do_tests);    
    }
  }
  check(gpu_qqnb, RelationalOperator::EQUAL, Approx(cpu_qqnb).margin(qqnb_tol), "Electrostatic "
        "non-bonded energies computed on the CPU and GPU do not agree.  Precision level in the "
        "calculation: " + getEnumerationName(prec) + ".", do_tests);
  check(gpu_ljnb, RelationalOperator::EQUAL, Approx(cpu_ljnb).margin(ljnb_tol), "Lennard-Jones "
        "non-bonded energies computed on the CPU and GPU do not agree.  Precision level in the "
        "calculation: " + getEnumerationName(prec) + ".", do_tests);  
  if (poly_ag.getImplicitSolventModel() != ImplicitSolventModel::NONE) {
    check(gpu_gbnb, RelationalOperator::EQUAL, Approx(cpu_gbnb).margin(gbnb_tol), "Generalized "
          "Born non-bonded energies computed on the CPU and GPU do not agree.  Precision level in "
          "the calculation: " + getEnumerationName(prec) + ".", do_tests);
  }
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  StopWatch timer;

  // Section 1
  section("Topology compilation and staging");
  
  // Section 2
  section("Coordinate compilation and staging");

  // Section 3
  section("Vacuum energy and force calculations");

  // Section 4
  section("Generalized Born energy and force calculations");

  // Section 5
  section("Virtual site placement");

  // Get the GPU specs.  Set of parameters for the work units and launch grids.  Preparing a cache
  // resource for the largest possible number of valence work unit atoms obviates the need to
  // worry about resizing it case by case, and does not contribute to cache pollution (individual
  // blocks' buffers are already padded by the warp size, hence each takes a set number of cache
  // lines to read).
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  int nblocks = gpu.getSMPCount();
  int nthreads = gpu.getMaxThreadsPerBlock();
  if (gpu.getArchMajor() == 6 && gpu.getArchMinor() == 1) {
    nblocks *= 2;
    nthreads /= 2;
  }
  
  // Collect coordinates and topologies
  const char osc = osSeparator();
  std::string buffer("&files\n  -p ");
  buffer += oe.getStormmSourcePath() + osc + "test" + osc + "Namelists" + osc + "topol" + osc +
            ".*.top\n  -c ";
  buffer += oe.getStormmSourcePath() + osc + "test" + osc + "Namelists" + osc + "coord" + osc +
            ".*.inpcrd\n&end\n";
  const TextFile tf(buffer, TextOrigin::RAM);
  int start_line = 0;
  FilesControls fcon(tf, &start_line);
  const SystemCache sysc(fcon, ExceptionResponse::SILENT, MapRotatableGroups::NO,
                         PrintSituation::OPEN_NEW, &timer);

  // Form the syntheses of topologies and coordinates
  section(1);
  const int nsys = sysc.getSystemCount();
  const TestPriority do_tests = (nsys > 0) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (nsys == 0) {
    rtWarn("No system topology and coordinate pairs were detected.  Subsequent tests will be "
           "skipped.", "test_hpc_synthesis");
  }
  std::vector<int> topology_indices(nsys, 0);
  for (int i = 0; i < nsys; i++) {
    topology_indices[i] = i;
  }
  AtomGraphSynthesis poly_ag(sysc.getSystemTopologyPointer(), topology_indices,
                             ExceptionResponse::WARN, gpu, &timer);
  StaticExclusionMaskSynthesis poly_se(poly_ag.getUniqueTopologies(),
                                       poly_ag.getTopologyIndices());
  poly_ag.loadNonbondedWorkUnits(poly_se);
  PhaseSpaceSynthesis poly_ps(sysc.getCoordinates(), sysc.getSystemTopologyPointer());
  PhaseSpaceSynthesis poly_ps_dbl(sysc.getCoordinates(), sysc.getSystemTopologyPointer(),
                                  36, 24, 36, 40);
  PhaseSpaceSynthesis poly_ps_sdbl(sysc.getCoordinates(), sysc.getSystemTopologyPointer(),
                                   72, 24, 36, 72);
  CoreKlManager launcher(gpu, poly_ag);
  check(poly_ag.getSystemCount(), RelationalOperator::EQUAL, poly_ps.getSystemCount(),
        "PhaseSpaceSynthesis and AtomGraphSynthesis objects formed from the same SystemCache have "
        "different numbers of systems inside of them.", do_tests);
  
  // Upload the compiled systems and check the results
  poly_ag.upload();
  poly_se.upload();
  poly_ps.upload();
  poly_ps_dbl.upload();
  poly_ps_sdbl.upload();
  std::vector<double> gpu_charges = poly_ag.getPartialCharges<double>(HybridTargetLevel::DEVICE);
  int padded_atom_count = 0;
  for (int i = 0; i < nsys; i++) {
    if (i == nsys - 1) {
      padded_atom_count += poly_ag.getSystemTopologyPointer(i)->getAtomCount();
    }
    else {
      padded_atom_count += roundUp(poly_ag.getSystemTopologyPointer(i)->getAtomCount(),
                                   warp_size_int);
    }
  }
  std::vector<double> rbt_charges(padded_atom_count, 0.0);
  int atom_offset = 0;
  for (int i = 0; i < nsys; i++) {
    const AtomGraph *iag_ptr = poly_ag.getSystemTopologyPointer(i);
    std::vector<double> ichg = iag_ptr->getPartialCharge<double>();
    const int natom = iag_ptr->getAtomCount();
    for (int j = 0; j < natom; j++) {
      rbt_charges[atom_offset + j] = ichg[j];
    }
    const int padded_natom = roundUp(natom, warp_size_int);
    if (i < nsys - 1) {
      for (int j = natom; j < padded_natom; j++) {
        gpu_charges[atom_offset + j] = 0.0;
      }
    }
    atom_offset += padded_natom;
  }
  check(gpu_charges, RelationalOperator::EQUAL, rbt_charges, "Charges pulled from the GPU in an "
        "AtomGraphSynthesis object do not meet expectations.", do_tests);

  // Check that different coordinate precisions translate to similar structures
  section(2);
  std::vector<int> sd_mismatch(nsys, 0);
  std::vector<int> sld_mismatch(nsys, 0);
  for (int i = 0; i < nsys; i++) {
    const CoordinateFrame s_frm  = poly_ps.exportCoordinates(i);
    const CoordinateFrame d_frm  = poly_ps_dbl.exportCoordinates(i);
    const CoordinateFrame ld_frm = poly_ps_sdbl.exportCoordinates(i);
    const CoordinateFrameReader s_frm_r  = s_frm.data();
    const CoordinateFrameReader d_frm_r  = d_frm.data();
    const CoordinateFrameReader ld_frm_r = ld_frm.data();
    const int natom = s_frm_r.natom;
    sd_mismatch[i] += (maxAbsoluteDifference(s_frm_r.xcrd, d_frm_r.xcrd, natom) > 1.0e-6);
    sd_mismatch[i] += (maxAbsoluteDifference(s_frm_r.ycrd, d_frm_r.ycrd, natom) > 1.0e-6);
    sd_mismatch[i] += (maxAbsoluteDifference(s_frm_r.zcrd, d_frm_r.zcrd, natom) > 1.0e-6);
    sld_mismatch[i] += (maxAbsoluteDifference(s_frm_r.xcrd, ld_frm_r.xcrd, natom) > 1.0e-6);
    sld_mismatch[i] += (maxAbsoluteDifference(s_frm_r.ycrd, ld_frm_r.ycrd, natom) > 1.0e-6);
    sld_mismatch[i] += (maxAbsoluteDifference(s_frm_r.zcrd, ld_frm_r.zcrd, natom) > 1.0e-6);
  }
  check(sd_mismatch, RelationalOperator::EQUAL, std::vector<int>(nsys, 0), "Coordinates of "
        "various sytems do not align when shifting from " +
        std::to_string(poly_ps.getGlobalPositionBits()) + " to " +
        std::to_string(poly_ps_dbl.getGlobalPositionBits()) + ".");
  check(sd_mismatch, RelationalOperator::EQUAL, std::vector<int>(nsys, 0), "Coordinates of "
        "various sytems do not align when shifting from " +
        std::to_string(poly_ps.getGlobalPositionBits()) + " to " +
        std::to_string(poly_ps_sdbl.getGlobalPositionBits()) + ".");

  // Allocate resources for various kernels.  The cache resources are allocated to accommodate
  // any foreseeable launch configuration.
  section(3);
  CacheResource valence_tb_space(8 * nblocks, maximum_valence_work_unit_atoms);
  CacheResource nonbond_tb_space(5 * nblocks, small_block_max_atoms);
  MolecularMechanicsControls mmctrl;
  mmctrl.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                               VwuGoal::ACCUMULATE, PrecisionModel::DOUBLE, PrecisionModel::DOUBLE,
                               poly_ag);
  ScoreCard sc(nsys, 1, 32);
  
  // Launch the valence evaluation kernel for small systems with only bonds, angles, dihedrals,
  // and 1:4 attenuated interactions.
  checkCompilationForces(&poly_ps_dbl, &mmctrl, &valence_tb_space, &nonbond_tb_space, poly_ag,
                         poly_se, AccumulationMethod::SPLIT, PrecisionModel::DOUBLE, gpu,
                         launcher, 3.5e-6, 2.0e-6, do_tests);
  checkCompilationEnergies(&poly_ps, &mmctrl, &valence_tb_space, &nonbond_tb_space, poly_ag,
                           poly_se, PrecisionModel::DOUBLE, gpu, launcher, 1.0e-6, 1.0e-6, 1.0e-6,
                           1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6,
                           1.0e-6, 1.0e-6, do_tests);

  // Check that super-high precision forms of the coordinates and force accumulation are OK in
  // double-precision mode.
  checkCompilationForces(&poly_ps_sdbl, &mmctrl, &valence_tb_space, &nonbond_tb_space,
                         poly_ag, poly_se, AccumulationMethod::SPLIT,
                         PrecisionModel::DOUBLE, gpu, launcher, 3.5e-7, 5.0e-7, do_tests);
  
  // Reconfigure the launch coordination for single-precision calculations
  mmctrl.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                               VwuGoal::ACCUMULATE, PrecisionModel::SINGLE, PrecisionModel::SINGLE,
                               poly_ag);
  checkCompilationForces(&poly_ps_dbl, &mmctrl, &valence_tb_space, &nonbond_tb_space, poly_ag,
                         poly_se, AccumulationMethod::SPLIT, PrecisionModel::SINGLE, gpu,
                         launcher, 3.5e-5, 2.0e-4, do_tests);
  checkCompilationForces(&poly_ps, &mmctrl, &valence_tb_space, &nonbond_tb_space, poly_ag, poly_se,
                         AccumulationMethod::WHOLE, PrecisionModel::SINGLE, gpu, launcher,
                         3.5e-5, 2.0e-4, do_tests);
  checkCompilationEnergies(&poly_ps, &mmctrl, &valence_tb_space, &nonbond_tb_space, poly_ag,
                           poly_se, PrecisionModel::SINGLE, gpu, launcher, 1.5e-5, 1.5e-5, 5.0e-6,
                           1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 6.0e-6, 2.2e-5, 1.0e-6,
                           1.0e-1, 3.5e-5, 3.5e-5, do_tests);

  // Create a set of larger systems, now involving CMAPs and other CHARMM force field terms
  const std::string topology_base = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string trpi_top_name = topology_base + osc + "trpcage.top";
  const std::string dhfr_top_name = topology_base + osc + "dhfr_cmap.top";
  const std::string alad_top_name = topology_base + osc + "ala_dipeptide.top";
  const std::string coordinate_base = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string trpi_crd_name = coordinate_base + osc + "trpcage.inpcrd";
  const std::string dhfr_crd_name = coordinate_base + osc + "dhfr_cmap.inpcrd";
  const std::string alad_crd_name = coordinate_base + osc + "ala_dipeptide.inpcrd";  
  const bool files_exist = (getDrivePathType(trpi_top_name) == DrivePathType::FILE &&
                            getDrivePathType(dhfr_top_name) == DrivePathType::FILE &&
                            getDrivePathType(alad_top_name) == DrivePathType::FILE &&
                            getDrivePathType(trpi_crd_name) == DrivePathType::FILE &&
                            getDrivePathType(dhfr_crd_name) == DrivePathType::FILE &&
                            getDrivePathType(alad_crd_name) == DrivePathType::FILE);
  AtomGraph trpi_ag, dhfr_ag, alad_ag;
  PhaseSpace trpi_ps, dhfr_ps, alad_ps;
  if (files_exist) {
    trpi_ag.buildFromPrmtop(trpi_top_name, ExceptionResponse::SILENT);
    dhfr_ag.buildFromPrmtop(dhfr_top_name, ExceptionResponse::SILENT);
    alad_ag.buildFromPrmtop(alad_top_name, ExceptionResponse::SILENT);
    trpi_ps.buildFromFile(trpi_crd_name);
    dhfr_ps.buildFromFile(dhfr_crd_name);
    alad_ps.buildFromFile(alad_crd_name);
  }
  else {
    rtWarn("Files for several systems in implicit solvent were not found.  Check the "
           "${STORMM_SOURCE} environment variable for validity.  Subsequent tests will be "
           "skipped.");
  }

  // Read some larger topologies, with CHARMM CMAP and other force field terms
  const std::vector<AtomGraph*> bigger_tops = { &trpi_ag, &dhfr_ag, &alad_ag };
  const std::vector<PhaseSpace> bigger_crds = { trpi_ps, dhfr_ps, alad_ps };
  PhaseSpaceSynthesis big_poly_ps(bigger_crds, bigger_tops);
  PhaseSpaceSynthesis big_poly_ps_dbl(bigger_crds, bigger_tops, 36, 24, 36, 40);
  PhaseSpaceSynthesis big_poly_ps_sdbl(bigger_crds, bigger_tops, 72, 24, 36, 72);
  const std::vector<int> big_top_indices = { 0, 1, 2 };
  AtomGraphSynthesis big_poly_ag(bigger_tops, big_top_indices, ExceptionResponse::SILENT, gpu,
                                 &timer);
  StaticExclusionMaskSynthesis big_poly_se(big_poly_ag.getUniqueTopologies(),
                                           big_poly_ag.getTopologyIndices());
  big_poly_ag.loadNonbondedWorkUnits(big_poly_se);
  CoreKlManager big_launcher(gpu, big_poly_ag);
  big_poly_ag.upload();
  big_poly_se.upload();
  big_poly_ps.upload();
  big_poly_ps_dbl.upload();
  big_poly_ps_sdbl.upload();

  // Reconfigure the work unit progress counters and launch double-precision calculations
  mmctrl.primeWorkUnitCounters(big_launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                               VwuGoal::ACCUMULATE, PrecisionModel::DOUBLE, PrecisionModel::DOUBLE,
                               big_poly_ag);
  checkCompilationForces(&big_poly_ps_dbl, &mmctrl, &valence_tb_space, &nonbond_tb_space,
                         big_poly_ag, big_poly_se, AccumulationMethod::SPLIT,
                         PrecisionModel::DOUBLE, gpu, big_launcher, 3.5e-6, 2.5e-5, do_tests);
  checkCompilationEnergies(&big_poly_ps, &mmctrl, &valence_tb_space, &nonbond_tb_space,
                           big_poly_ag, big_poly_se, PrecisionModel::DOUBLE, gpu, big_launcher,
                           1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 6.0e-6, 1.0e-6, 1.0e-6,
                           1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, do_tests);
  checkCompilationForces(&big_poly_ps_sdbl, &mmctrl, &valence_tb_space, &nonbond_tb_space,
                         big_poly_ag, big_poly_se, AccumulationMethod::SPLIT,
                         PrecisionModel::DOUBLE, gpu, big_launcher, 3.5e-6, 2.5e-5, do_tests);

  // Reconfigure again and launch single-precision calculations on the "big" systems
  mmctrl.primeWorkUnitCounters(big_launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                               VwuGoal::ACCUMULATE, PrecisionModel::SINGLE, PrecisionModel::SINGLE,
                               big_poly_ag);
  checkCompilationForces(&big_poly_ps, &mmctrl, &valence_tb_space, &nonbond_tb_space, big_poly_ag,
                         big_poly_se, AccumulationMethod::SPLIT, PrecisionModel::SINGLE,
                         gpu, big_launcher, 7.5e-5, 3.0e-3, do_tests);
  checkCompilationForces(&big_poly_ps, &mmctrl, &valence_tb_space, &nonbond_tb_space, big_poly_ag,
                         big_poly_se, AccumulationMethod::WHOLE, PrecisionModel::SINGLE,
                         gpu, big_launcher, 7.5e-5, 3.0e-3, do_tests);
  checkCompilationEnergies(&big_poly_ps, &mmctrl, &valence_tb_space, &nonbond_tb_space,
                           big_poly_ag, big_poly_se, PrecisionModel::SINGLE, gpu, big_launcher,
                           1.5e-4, 6.5e-5, 9.0e-5, 1.5e-5, 6.0e-5, 3.0e-5, 1.0e-5, 7.5e-5, 2.2e-4,
                           1.0e-6, 7.5e-4, 7.5e-3, 7.5e-3, do_tests);
  
  // Tweak the original synthesis to incorporate a Generalized Born model.
  section(4);
  const NeckGeneralizedBornTable ngb_tables;
  const std::vector<ImplicitSolventModel> amber_isms = {
    ImplicitSolventModel::OBC_GB, ImplicitSolventModel::HCT_GB, ImplicitSolventModel::OBC_GB_II,
    ImplicitSolventModel::NECK_GB, ImplicitSolventModel::NECK_GB_II };
  const std::vector<AtomicRadiusSet> appropriate_ism_radii = {
    AtomicRadiusSet::AMBER6, AtomicRadiusSet::MBONDI2, AtomicRadiusSet::MBONDI,
    AtomicRadiusSet::BONDI, AtomicRadiusSet::MBONDI3 };

  // The compiled topology synthesis must have some implicit solvent model set in order to trigger
  // the initialization of the relevant counters.  Otherwise, one initialization for the whole
  // series is fine.
  poly_ag.setImplicitSolventModel(amber_isms[0], ngb_tables, appropriate_ism_radii[0]);
  mmctrl.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                               VwuGoal::ACCUMULATE, PrecisionModel::DOUBLE, PrecisionModel::DOUBLE,
                               poly_ag); 
  for (size_t i = 0; i < amber_isms.size(); i++) {
    poly_ag.setImplicitSolventModel(amber_isms[i], ngb_tables, appropriate_ism_radii[i]);
    poly_ag.upload();
    const std::string side_note = "GB model: " + getEnumerationName(amber_isms[i]) +
                                  " (" + getEnumerationName(appropriate_ism_radii[i]) +
                                  " radii).";
    checkCompilationForces(&poly_ps_dbl, &mmctrl, &valence_tb_space, &nonbond_tb_space, poly_ag,
                           poly_se, AccumulationMethod::SPLIT, PrecisionModel::DOUBLE, gpu,
                           launcher, 3.5e-6, 2.0e-6, do_tests, side_note, false);
    checkCompilationEnergies(&poly_ps_dbl, &mmctrl, &valence_tb_space, &nonbond_tb_space, poly_ag,
                             poly_se, PrecisionModel::DOUBLE, gpu, launcher, 1.0e-6, 1.0e-6,
                             1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6,
                             1.0e-6, 1.0e-6, 1.0e-6, do_tests, false);
    checkCompilationForces(&poly_ps_sdbl, &mmctrl, &valence_tb_space, &nonbond_tb_space, poly_ag,
                           poly_se, AccumulationMethod::SPLIT, PrecisionModel::DOUBLE, gpu,
                           launcher, 3.5e-6, 2.0e-6, do_tests, side_note, false);
    checkCompilationEnergies(&poly_ps_sdbl, &mmctrl, &valence_tb_space, &nonbond_tb_space, poly_ag,
                             poly_se, PrecisionModel::DOUBLE, gpu, launcher, 1.0e-6, 1.0e-6,
                             1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6,
                             1.0e-6, 1.0e-6, 1.0e-6, do_tests, false);
  }
  mmctrl.primeWorkUnitCounters(launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                               VwuGoal::ACCUMULATE, PrecisionModel::SINGLE, PrecisionModel::SINGLE,
                               poly_ag);
  for (size_t i = 0; i < amber_isms.size(); i++) {
    poly_ag.setImplicitSolventModel(amber_isms[i], ngb_tables, appropriate_ism_radii[i]);
    poly_ag.upload();
    const std::string side_note = "GB model: " + getEnumerationName(amber_isms[i]) +
                                  " (" + getEnumerationName(appropriate_ism_radii[i]) +
                                  " radii).";
    checkCompilationForces(&poly_ps, &mmctrl, &valence_tb_space, &nonbond_tb_space, poly_ag,
                           poly_se, AccumulationMethod::SPLIT, PrecisionModel::SINGLE, gpu,
                           launcher, 9.0e-6, 2.5e-4, do_tests, side_note, false);
    checkCompilationEnergies(&poly_ps, &mmctrl, &valence_tb_space, &nonbond_tb_space, poly_ag,
                             poly_se, PrecisionModel::SINGLE, gpu, launcher, 1.5e-5, 1.5e-5,
                             5.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 6.0e-6, 2.2e-5, 1.0e-6,
                             1.0e-1, 3.5e-5, 3.5e-5, do_tests, false);
  }

  // Read some topologies with virtual sites.  First, test the forces that appear to act on the
  // virtual sites.  Add restraints to these ligands.
  section(5);
  const std::string brbz_top_name = topology_base + osc + "bromobenzene_vs_iso.top";
  const std::string lig1_top_name = topology_base + osc + "stereo_L1_vs.top";
  const std::string lig2_top_name = topology_base + osc + "symmetry_L1_vs.top";
  const std::string lig3_top_name = topology_base + osc + "drug_example_vs_iso.top";
  const std::string lig4_top_name = topology_base + osc + "drug_example_iso.top";
  const std::string brbz_crd_name = coordinate_base + osc + "bromobenzene_vs_iso.inpcrd";
  const std::string lig1_crd_name = coordinate_base + osc + "stereo_L1_vs.inpcrd";
  const std::string lig2_crd_name = coordinate_base + osc + "symmetry_L1_vs.inpcrd";
  const std::string lig3_crd_name = coordinate_base + osc + "drug_example_vs_iso.inpcrd";
  const std::string lig4_crd_name = coordinate_base + osc + "drug_example_iso.inpcrd";
  const bool ligands_exist = (getDrivePathType(brbz_top_name) == DrivePathType::FILE &&
                              getDrivePathType(lig1_top_name) == DrivePathType::FILE &&
                              getDrivePathType(lig2_top_name) == DrivePathType::FILE &&
                              getDrivePathType(lig3_top_name) == DrivePathType::FILE &&
                              getDrivePathType(lig4_top_name) == DrivePathType::FILE &&
                              getDrivePathType(brbz_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(lig1_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(lig2_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(lig3_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(lig4_crd_name) == DrivePathType::FILE);
  AtomGraph brbz_ag, lig1_ag, lig2_ag, lig3_ag, lig4_ag;
  PhaseSpace brbz_ps, lig1_ps, lig2_ps, lig3_ps, lig4_ps;
  if (ligands_exist) {
    brbz_ag.buildFromPrmtop(brbz_top_name);
    lig1_ag.buildFromPrmtop(lig1_top_name);
    lig2_ag.buildFromPrmtop(lig2_top_name);
    lig3_ag.buildFromPrmtop(lig3_top_name);
    lig4_ag.buildFromPrmtop(lig4_top_name);
    brbz_ps.buildFromFile(brbz_crd_name);
    lig1_ps.buildFromFile(lig1_crd_name);
    lig2_ps.buildFromFile(lig2_crd_name);
    lig3_ps.buildFromFile(lig3_crd_name);
    lig4_ps.buildFromFile(lig4_crd_name);
  }
  RestraintApparatus brbz_ra = assembleRestraints(&brbz_ag, brbz_ps);
  RestraintApparatus lig1_ra = assembleRestraints(&lig1_ag, lig1_ps);
  RestraintApparatus lig2_ra = assembleRestraints(&lig2_ag, lig2_ps);
  RestraintApparatus lig3_ra = assembleRestraints(&lig3_ag, lig3_ps);
  RestraintApparatus lig4_ra = assembleRestraints(&lig4_ag, lig3_ps);
  const std::vector<AtomGraph*> ligand_ag_list = { &brbz_ag, &lig1_ag, &lig2_ag, &lig3_ag,
                                                   &lig4_ag };
  const std::vector<PhaseSpace> ligand_ps_list = {  brbz_ps,  lig1_ps,  lig2_ps,  lig3_ps,
                                                    lig4_ps };
  const std::vector<RestraintApparatus*> ligand_ra_list = { &brbz_ra, &lig1_ra, &lig2_ra,
                                                            &lig3_ra, &lig4_ra };
  const std::vector<int> ligand_minitile = { 0, 1, 2, 3, 4, 0, 1, 4, 3, 2, 2, 3, 2, 1, 0, 4 };
  const int lm_size = ligand_minitile.size();
  const int n_ligand_tile_reps = 7;
  const int nligands = n_ligand_tile_reps * lm_size;
  std::vector<int> ligand_tiling(112, 0);
  for (int i = 0; i < nligands; i += lm_size) {
    for (int j = 0; j < lm_size; j++) {
      ligand_tiling[i + j] = ligand_minitile[j];
    }
  }
  PhaseSpaceSynthesis ligand_poly_ps(ligand_ps_list, ligand_ag_list, ligand_tiling);
  PhaseSpaceSynthesis ligand_poly_ps_dbl(ligand_ps_list, ligand_ag_list, ligand_tiling, 40, 24,
                                         36, 44);
  PhaseSpaceSynthesis ligand_poly_ps_sdbl(ligand_ps_list, ligand_ag_list, ligand_tiling, 72, 24,
                                          36, 72);
  AtomGraphSynthesis ligand_poly_ag(ligand_ag_list, ligand_ra_list, ligand_tiling, ligand_tiling,
                                    ExceptionResponse::WARN, gpu, &timer);
  StaticExclusionMaskSynthesis ligand_poly_se(ligand_poly_ag.getUniqueTopologies(),
                                              ligand_poly_ag.getTopologyIndices());
  ligand_poly_ag.loadNonbondedWorkUnits(ligand_poly_se);
  CoreKlManager ligand_launcher(gpu, ligand_poly_ag);
  ligand_poly_ag.upload();
  ligand_poly_se.upload();
  ligand_poly_ps.upload();
  ligand_poly_ps_dbl.upload();
  ligand_poly_ps_sdbl.upload();

  // Reconfigure the progress counters for the large array of ligand systems and launch
  // double-precision calculations.
  mmctrl.primeWorkUnitCounters(ligand_launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                               VwuGoal::ACCUMULATE, PrecisionModel::DOUBLE, PrecisionModel::DOUBLE,
                               ligand_poly_ag);
  checkCompilationForces(&ligand_poly_ps_dbl, &mmctrl, &valence_tb_space, &nonbond_tb_space,
                         ligand_poly_ag, ligand_poly_se, AccumulationMethod::SPLIT,
                         PrecisionModel::DOUBLE, gpu, ligand_launcher, 3.5e-6, 2.0e-6, do_tests);
  checkCompilationEnergies(&ligand_poly_ps, &mmctrl, &valence_tb_space, &nonbond_tb_space,
                           ligand_poly_ag, ligand_poly_se, PrecisionModel::DOUBLE, gpu,
                           ligand_launcher, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6,
                           6.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, 1.0e-6, do_tests);

  // Remarkably, the "super-double" representation with ultra-high precision coordinates cannot
  // get forces any closer than 5.0e-7 kcal/mol-A, due to something that happens in the last bit
  // of one of the improper dihedral terms.  The geometry is planar, but depending on whether the
  // final bit is 0 or 1, the angle can change by a much more significant amount, and not even a
  // guard against the arccos instability will change the outcome.
  checkCompilationForces(&ligand_poly_ps_sdbl, &mmctrl, &valence_tb_space, &nonbond_tb_space,
                         ligand_poly_ag, ligand_poly_se, AccumulationMethod::SPLIT,
                         PrecisionModel::DOUBLE, gpu, ligand_launcher, 3.5e-8, 5.0e-7, do_tests);

  // Reconfigure one final time and perform single-precision calculations on the ligands
  mmctrl.primeWorkUnitCounters(ligand_launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                               VwuGoal::ACCUMULATE, PrecisionModel::SINGLE, PrecisionModel::SINGLE,
                               ligand_poly_ag);
  checkCompilationForces(&ligand_poly_ps, &mmctrl, &valence_tb_space, &nonbond_tb_space,
                         ligand_poly_ag, ligand_poly_se, AccumulationMethod::SPLIT,
                         PrecisionModel::SINGLE, gpu, ligand_launcher, 7.5e-5, 3.0e-3, do_tests);
  checkCompilationForces(&ligand_poly_ps, &mmctrl, &valence_tb_space, &nonbond_tb_space,
                         ligand_poly_ag, ligand_poly_se, AccumulationMethod::WHOLE,
                         PrecisionModel::SINGLE, gpu, ligand_launcher, 7.5e-5, 3.0e-3, do_tests);
  checkCompilationEnergies(&ligand_poly_ps, &mmctrl, &valence_tb_space, &nonbond_tb_space,
                           ligand_poly_ag, ligand_poly_se, PrecisionModel::SINGLE, gpu,
                           ligand_launcher, 1.5e-4, 2.2e-5, 9.0e-5, 1.5e-5, 6.0e-5, 3.0e-5,
                           6.0e-6, 7.5e-5, 2.2e-4, 1.0e-5, 2.2e-4, 3.0e-5, 3.0e-5, do_tests);

  // Scramble and replace virtual sites on the GPU
  Xoroshiro128pGenerator xrs(618608377);
  for (int i = 0; i < nligands; i++) {
    const AtomGraph *iag_ptr = ligand_poly_ag.getSystemTopologyPointer(i);
    PsSynthesisWriter ligand_psw = ligand_poly_ps.data();
    PsSynthesisWriter ligand_dbl_psw = ligand_poly_ps_dbl.data();
    PsSynthesisWriter ligand_sdbl_psw = ligand_poly_ps_sdbl.data();
    const ChemicalDetailsKit icdk = iag_ptr->getChemicalDetailsKit();
    for (int j = 0; j < icdk.natom; j++) {
      if (icdk.z_numbers[j] == 0 || true) {
        const double pert_fac = (icdk.z_numbers[j] == 0) ? 1.0 : 0.05;
        const int synth_idx = ligand_psw.atom_starts[i] + j;
        hostSplitAccumulation(xrs.gaussianRandomNumber() * ligand_psw.gpos_scale * pert_fac,
                              &ligand_psw.xcrd[synth_idx], &ligand_psw.xcrd_ovrf[synth_idx]);
        hostSplitAccumulation(xrs.gaussianRandomNumber() * ligand_psw.gpos_scale * pert_fac,
                              &ligand_psw.ycrd[synth_idx], &ligand_psw.ycrd_ovrf[synth_idx]);
        hostSplitAccumulation(xrs.gaussianRandomNumber() * ligand_psw.gpos_scale * pert_fac,
                              &ligand_psw.zcrd[synth_idx], &ligand_psw.zcrd_ovrf[synth_idx]);
        hostSplitAccumulation(xrs.gaussianRandomNumber() * ligand_dbl_psw.gpos_scale * pert_fac,
                              &ligand_dbl_psw.xcrd[synth_idx],
                              &ligand_dbl_psw.xcrd_ovrf[synth_idx]);
        hostSplitAccumulation(xrs.gaussianRandomNumber() * ligand_dbl_psw.gpos_scale * pert_fac,
                              &ligand_dbl_psw.ycrd[synth_idx],
                              &ligand_dbl_psw.ycrd_ovrf[synth_idx]);
        hostSplitAccumulation(xrs.gaussianRandomNumber() * ligand_dbl_psw.gpos_scale * pert_fac,
                              &ligand_dbl_psw.zcrd[synth_idx],
                              &ligand_dbl_psw.zcrd_ovrf[synth_idx]);
        hostSplitAccumulation(xrs.gaussianRandomNumber() * ligand_sdbl_psw.gpos_scale * pert_fac,
                              &ligand_sdbl_psw.xcrd[synth_idx],
                              &ligand_sdbl_psw.xcrd_ovrf[synth_idx]);
        hostSplitAccumulation(xrs.gaussianRandomNumber() * ligand_sdbl_psw.gpos_scale * pert_fac,
                              &ligand_sdbl_psw.ycrd[synth_idx],
                              &ligand_sdbl_psw.ycrd_ovrf[synth_idx]);
        hostSplitAccumulation(xrs.gaussianRandomNumber() * ligand_sdbl_psw.gpos_scale * pert_fac,
                              &ligand_sdbl_psw.zcrd[synth_idx],
                              &ligand_sdbl_psw.zcrd_ovrf[synth_idx]);
      }
    }
  }
  ligand_poly_ps.upload();
  ligand_poly_ps_dbl.upload();
  ligand_poly_ps_sdbl.upload();

  // Copy the host-based scrambled virtual site positions and replace them using vetted
  // single-system methods.  Replace the virtual sites on the GPU.  Check the results against
  // CPU-based computations.
  std::vector<CoordinateFrame> cpulig_cf_vec, cpulig_dbl_cf_vec, cpulig_sdbl_cf_vec;
  cpulig_cf_vec.reserve(nligands);
  cpulig_dbl_cf_vec.reserve(nligands);
  cpulig_sdbl_cf_vec.reserve(nligands);
  for (int i = 0; i < nligands; i++) {
    const AtomGraph *iag_ptr = ligand_poly_ag.getSystemTopologyPointer(i);
    cpulig_cf_vec.emplace_back(ligand_poly_ps.exportCoordinates(i));
    cpulig_dbl_cf_vec.emplace_back(ligand_poly_ps_dbl.exportCoordinates(i));
    cpulig_sdbl_cf_vec.emplace_back(ligand_poly_ps_sdbl.exportCoordinates(i));
    placeVirtualSites(&cpulig_cf_vec[i], iag_ptr);
    placeVirtualSites(&cpulig_dbl_cf_vec[i], iag_ptr);
    placeVirtualSites(&cpulig_sdbl_cf_vec[i], iag_ptr);
  }

  // By convention, the standalone virtual site placement kernels act on frame positions stored in
  // the developing positions arrays of the synthesis.  In order to mock this effect, the
  // coordinate syntheses must be advanced one cycle.  This will be undone after the following
  // tests.
  const VirtualSiteActivity v_activity = VirtualSiteActivity::PLACEMENT;
  ligand_poly_ps.updateCyclePosition();
  ligand_poly_ps_dbl.updateCyclePosition();
  ligand_poly_ps_sdbl.updateCyclePosition();
  launchVirtualSiteHandling(PrecisionModel::SINGLE, v_activity, &ligand_poly_ps, &valence_tb_space,
                           ligand_poly_ag, ligand_launcher);
  launchVirtualSiteHandling(PrecisionModel::DOUBLE, v_activity, &ligand_poly_ps_dbl,
                            &valence_tb_space, ligand_poly_ag, ligand_launcher);
  launchVirtualSiteHandling(PrecisionModel::DOUBLE, v_activity, &ligand_poly_ps_sdbl,
                            &valence_tb_space, ligand_poly_ag, ligand_launcher);
  ligand_poly_ps.updateCyclePosition();
  ligand_poly_ps_dbl.updateCyclePosition();
  ligand_poly_ps_sdbl.updateCyclePosition();
  std::vector<CoordinateFrame> gpulig_cf_vec, gpulig_dbl_cf_vec, gpulig_sdbl_cf_vec;
  gpulig_cf_vec.reserve(nligands);
  gpulig_dbl_cf_vec.reserve(nligands);
  gpulig_sdbl_cf_vec.reserve(nligands);
  int nvs_atom = 0;
  for (int i = 0; i < nligands; i++) {

    // It's not as efficient to extract coordinates from the GPU level in this way, but the
    // systems are not large and this leaves the CPU coordinates untouched.
    const TrajectoryKind tcoord = TrajectoryKind::POSITIONS;
    const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
    gpulig_cf_vec.emplace_back(ligand_poly_ps.exportCoordinates(i, tcoord, devc_tier));
    gpulig_dbl_cf_vec.emplace_back(ligand_poly_ps_dbl.exportCoordinates(i, tcoord, devc_tier));
    gpulig_sdbl_cf_vec.emplace_back(ligand_poly_ps_sdbl.exportCoordinates(i, tcoord, devc_tier));
    const AtomGraph *iag_ptr = ligand_poly_ag.getSystemTopologyPointer(i);
    nvs_atom += iag_ptr->getVirtualSiteCount();
  }
  std::vector<double> cpu_positions(3 * nvs_atom), gpu_positions(3 * nvs_atom);
  std::vector<double> cpu_dbl_positions(3 * nvs_atom), gpu_dbl_positions(3 * nvs_atom);
  std::vector<double> cpu_sdbl_positions(3 * nvs_atom), gpu_sdbl_positions(3 * nvs_atom);
  int vscon = 0;
  for (int i = 0; i < nligands; i++) {
    const AtomGraph *iag_ptr = ligand_poly_ag.getSystemTopologyPointer(i);
    CoordinateFrameWriter cpu_cfr      = cpulig_cf_vec[i].data();
    CoordinateFrameWriter cpu_dbl_cfr  = cpulig_dbl_cf_vec[i].data();
    CoordinateFrameWriter cpu_sdbl_cfr = cpulig_sdbl_cf_vec[i].data();
    CoordinateFrameWriter gpu_cfr      = gpulig_cf_vec[i].data();
    CoordinateFrameWriter gpu_dbl_cfr  = gpulig_dbl_cf_vec[i].data();
    CoordinateFrameWriter gpu_sdbl_cfr = gpulig_sdbl_cf_vec[i].data();
    const ChemicalDetailsKit icdk = iag_ptr->getChemicalDetailsKit();
    for (int j = 0; j < icdk.natom; j++) {
      if (icdk.z_numbers[j] == 0) {
        cpu_positions[ 3 * vscon     ] = cpu_cfr.xcrd[j];
        cpu_positions[(3 * vscon) + 1] = cpu_cfr.ycrd[j];
        cpu_positions[(3 * vscon) + 2] = cpu_cfr.zcrd[j];
        gpu_positions[ 3 * vscon     ] = gpu_cfr.xcrd[j];
        gpu_positions[(3 * vscon) + 1] = gpu_cfr.ycrd[j];
        gpu_positions[(3 * vscon) + 2] = gpu_cfr.zcrd[j];
        cpu_dbl_positions[ 3 * vscon     ] = cpu_dbl_cfr.xcrd[j];
        cpu_dbl_positions[(3 * vscon) + 1] = cpu_dbl_cfr.ycrd[j];
        cpu_dbl_positions[(3 * vscon) + 2] = cpu_dbl_cfr.zcrd[j];
        gpu_dbl_positions[ 3 * vscon     ] = gpu_dbl_cfr.xcrd[j];
        gpu_dbl_positions[(3 * vscon) + 1] = gpu_dbl_cfr.ycrd[j];
        gpu_dbl_positions[(3 * vscon) + 2] = gpu_dbl_cfr.zcrd[j];
        cpu_sdbl_positions[ 3 * vscon     ] = cpu_sdbl_cfr.xcrd[j];
        cpu_sdbl_positions[(3 * vscon) + 1] = cpu_sdbl_cfr.ycrd[j];
        cpu_sdbl_positions[(3 * vscon) + 2] = cpu_sdbl_cfr.zcrd[j];
        gpu_sdbl_positions[ 3 * vscon     ] = gpu_sdbl_cfr.xcrd[j];
        gpu_sdbl_positions[(3 * vscon) + 1] = gpu_sdbl_cfr.ycrd[j];
        gpu_sdbl_positions[(3 * vscon) + 2] = gpu_sdbl_cfr.zcrd[j];
      }
    }
  }
  check(cpu_positions, RelationalOperator::EQUAL, gpu_positions, "Positions of virtual sites "
        "placed by the GPU kernel do not agree with those placed by the CPU function.  "
        "Precision level: " + getEnumerationName(PrecisionModel::SINGLE) + ".  Bits after "
        "the decimal: " + std::to_string(ligand_poly_ps.getGlobalPositionBits()) + ".", do_tests);
  check(cpu_dbl_positions, RelationalOperator::EQUAL, gpu_dbl_positions, "Positions of virtual "
        "sites placed by the GPU kernel do not agree with those placed by the CPU function.  "
        "Precision level: " + getEnumerationName(PrecisionModel::DOUBLE) + ".  Bits after "
        "the decimal: " + std::to_string(ligand_poly_ps_dbl.getGlobalPositionBits()) + ".",
        do_tests);
  check(cpu_sdbl_positions, RelationalOperator::EQUAL, gpu_sdbl_positions, "Positions of virtual "
        "sites placed by the GPU kernel do not agree with those placed by the CPU function.  "
        "Precision level: " + getEnumerationName(PrecisionModel::DOUBLE) + ".  Bits after the "
        "decimal: " + std::to_string(ligand_poly_ps_sdbl.getGlobalPositionBits()) + ".", do_tests);
  
  // Now, download the updated positions to the CPU.  Recompute forces with the new positions
  // and compare with CPU-derived forces on those new positions (checking the virtual site
  // positions by exporting CoordinateFrame objects and then work in them did not change the
  // virtual site positions on the CPU memory of the original PhaseSpaceSynthesis objects).
  ligand_poly_ps.upload();
  ligand_poly_ps_dbl.upload();
  ligand_poly_ps_sdbl.upload();
  mmctrl.primeWorkUnitCounters(ligand_launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                               VwuGoal::ACCUMULATE, PrecisionModel::SINGLE, PrecisionModel::SINGLE,
                               ligand_poly_ag);
  checkCompilationForces(&ligand_poly_ps, &mmctrl, &valence_tb_space, &nonbond_tb_space,
                         ligand_poly_ag, ligand_poly_se, AccumulationMethod::SPLIT,
                         PrecisionModel::SINGLE, gpu, ligand_launcher, 7.5e-5, 3.0e-3, do_tests,
                         "(Following virtual site replacement.)");
  mmctrl.primeWorkUnitCounters(ligand_launcher, EvaluateForce::YES, EvaluateEnergy::YES,
                               VwuGoal::ACCUMULATE, PrecisionModel::DOUBLE, PrecisionModel::DOUBLE,
                               ligand_poly_ag);
  checkCompilationForces(&ligand_poly_ps_dbl, &mmctrl, &valence_tb_space, &nonbond_tb_space,
                         ligand_poly_ag, ligand_poly_se, AccumulationMethod::SPLIT,
                         PrecisionModel::DOUBLE, gpu, ligand_launcher, 3.5e-6, 2.0e-6, do_tests,
                         "(Following virtual site replacement.)");
  checkCompilationForces(&ligand_poly_ps_sdbl, &mmctrl, &valence_tb_space, &nonbond_tb_space,
                         ligand_poly_ag, ligand_poly_se, AccumulationMethod::SPLIT,
                         PrecisionModel::DOUBLE, gpu, ligand_launcher, 3.5e-6, 2.0e-6, do_tests,
                         "(Following virtual site replacement.)");

  // Download the forces computed on all particles (prior to virtual site force transmission)
  ligand_poly_ps.download(TrajectoryKind::FORCES);
  ligand_poly_ps_dbl.download(TrajectoryKind::FORCES);
  ligand_poly_ps_sdbl.download(TrajectoryKind::FORCES);

  // Transmit forces from the virtual sites back to their frame atoms.  Attempt similar processes
  // on the CPU, using PhaseSpace objects to process one system at a time.
  const VirtualSiteActivity t_activity = VirtualSiteActivity::TRANSMIT_FORCES;
  launchVirtualSiteHandling(PrecisionModel::SINGLE, t_activity, &ligand_poly_ps, &valence_tb_space,
                           ligand_poly_ag, ligand_launcher);
  launchVirtualSiteHandling(PrecisionModel::DOUBLE, t_activity, &ligand_poly_ps_dbl,
                            &valence_tb_space, ligand_poly_ag, ligand_launcher);
  launchVirtualSiteHandling(PrecisionModel::DOUBLE, t_activity, &ligand_poly_ps_sdbl,
                            &valence_tb_space, ligand_poly_ag, ligand_launcher);

  // Check the force transmission kernel with the equivalent CPU routines
  std::vector<PhaseSpace> cpulig_ps_vec, cpulig_dbl_ps_vec, cpulig_sdbl_ps_vec;
  std::vector<PhaseSpace> gpulig_ps_vec, gpulig_dbl_ps_vec, gpulig_sdbl_ps_vec;
  cpulig_ps_vec.reserve(nligands);
  cpulig_dbl_ps_vec.reserve(nligands);
  cpulig_sdbl_ps_vec.reserve(nligands);
  gpulig_ps_vec.reserve(nligands);
  gpulig_dbl_ps_vec.reserve(nligands);
  gpulig_sdbl_ps_vec.reserve(nligands);
  std::vector<double> cpu_forces, cpu_dbl_forces, cpu_sdbl_forces;
  std::vector<double> gpu_forces, gpu_dbl_forces, gpu_sdbl_forces;
  for (int i = 0; i < nligands; i++) {

    // Obtain coordinates and un-transmitted forces on the CPU
    cpulig_ps_vec.emplace_back(ligand_poly_ps.exportSystem(i));
    cpulig_dbl_ps_vec.emplace_back(ligand_poly_ps_dbl.exportSystem(i));
    cpulig_sdbl_ps_vec.emplace_back(ligand_poly_ps_sdbl.exportSystem(i));

    // Obtain coordinates and transmitted forces from the GPU (virtual site forces should be zero)
    const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
    gpulig_ps_vec.emplace_back(ligand_poly_ps.exportSystem(i, devc_tier));
    gpulig_dbl_ps_vec.emplace_back(ligand_poly_ps_dbl.exportSystem(i, devc_tier));
    gpulig_sdbl_ps_vec.emplace_back(ligand_poly_ps_sdbl.exportSystem(i, devc_tier));
    
    // Transmit forces, via CPU routines, in the respective individual system coordinate objects
    const AtomGraph *iag_ptr = ligand_poly_ag.getSystemTopologyPointer(i);
    transmitVirtualSiteForces(&cpulig_ps_vec[i], *iag_ptr);
    transmitVirtualSiteForces(&cpulig_dbl_ps_vec[i], *iag_ptr);
    transmitVirtualSiteForces(&cpulig_sdbl_ps_vec[i], *iag_ptr);

    // Test forces from each process for consistency
    const TrajectoryKind tforce = TrajectoryKind::FORCES;
    std::vector<double> cpu_ilfrc      = cpulig_ps_vec[i].getInterlacedCoordinates(tforce);
    std::vector<double> cpu_dbl_ilfrc  = cpulig_dbl_ps_vec[i].getInterlacedCoordinates(tforce);
    std::vector<double> cpu_sdbl_ilfrc = cpulig_sdbl_ps_vec[i].getInterlacedCoordinates(tforce);
    std::vector<double> gpu_ilfrc      = gpulig_ps_vec[i].getInterlacedCoordinates(tforce);
    std::vector<double> gpu_dbl_ilfrc  = gpulig_dbl_ps_vec[i].getInterlacedCoordinates(tforce);
    std::vector<double> gpu_sdbl_ilfrc = gpulig_sdbl_ps_vec[i].getInterlacedCoordinates(tforce);
    cpu_forces.insert(cpu_forces.end(), cpu_ilfrc.begin(), cpu_ilfrc.end());
    cpu_dbl_forces.insert(cpu_dbl_forces.end(), cpu_dbl_ilfrc.begin(), cpu_dbl_ilfrc.end());
    cpu_sdbl_forces.insert(cpu_sdbl_forces.end(), cpu_sdbl_ilfrc.begin(), cpu_sdbl_ilfrc.end());
    gpu_forces.insert(gpu_forces.end(), gpu_ilfrc.begin(), gpu_ilfrc.end());
    gpu_dbl_forces.insert(gpu_dbl_forces.end(), gpu_dbl_ilfrc.begin(), gpu_dbl_ilfrc.end());
    gpu_sdbl_forces.insert(gpu_sdbl_forces.end(), gpu_sdbl_ilfrc.begin(), gpu_sdbl_ilfrc.end());
  }
  check(cpu_forces, RelationalOperator::EQUAL, Approx(gpu_forces).margin(1.0e-5),
        "Force transmission from virtual sites produces different results when accomplished by "
        "C++ functions and the standalone GPU kernel.  Precision model: " +
        getEnumerationName(PrecisionModel::SINGLE) + ".");
  check(cpu_dbl_forces, RelationalOperator::EQUAL, Approx(gpu_dbl_forces).margin(1.0e-5),
        "Force transmission from virtual sites produces different results when accomplished by "
        "C++ functions and the standalone GPU kernel.  Precision model: " +
        getEnumerationName(PrecisionModel::DOUBLE) + ".  Bits after the decimal (coordinates / "
        "forces): " + std::to_string(ligand_poly_ps_dbl.getGlobalPositionBits()) + " / " +
        std::to_string(ligand_poly_ps_dbl.getForceAccumulationBits()) + ".");
  check(cpu_sdbl_forces, RelationalOperator::EQUAL, Approx(gpu_sdbl_forces).margin(1.0e-5),
        "Force transmission from virtual sites produces different results when accomplished by "
        "C++ functions and the standalone GPU kernel.  Precision model: " +
        getEnumerationName(PrecisionModel::DOUBLE) + ".  Bits after the decimal (coordinates / "
        "forces): " + std::to_string(ligand_poly_ps_sdbl.getGlobalPositionBits()) + " / " +
        std::to_string(ligand_poly_ps_sdbl.getForceAccumulationBits()) + ".");

  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
