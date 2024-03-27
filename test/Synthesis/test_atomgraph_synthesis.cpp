#include <vector>
#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/generalized_born.h"
#include "../../src/Constants/scaling.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/eval_synthesis.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Restraints/bounded_restraint.h"
#include "../../src/Restraints/restraint_apparatus.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/static_mask_synthesis.h"
#include "../../src/Synthesis/synthesis_abstracts.h"
#include "../../src/Synthesis/synthesis_enumerators.h"
#include "../../src/Synthesis/valence_workunit.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/unit_test.h"
#include "assemble_restraints.h"

using stormm::constants::ExceptionResponse;
using stormm::constants::verytiny;
#ifndef STORMM_USE_HPC
using stormm::data_types::double2;
using stormm::data_types::double3;
using stormm::data_types::double4;
using stormm::data_types::float2;
#endif
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::errors::rtWarn;
using stormm::random::Xoroshiro128pGenerator;
using stormm::restraints::BoundedRestraint;
using stormm::restraints::RestraintApparatus;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using namespace stormm::energy;
using namespace stormm::generalized_born_defaults;
using namespace stormm::synthesis;
using namespace stormm::topology;
using namespace stormm::trajectory;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Simple enumerator to call for non-bonded computations based on a topology synthesis.
//-------------------------------------------------------------------------------------------------
enum class EvaluateNonbonded {
  YES, NO
};

//-------------------------------------------------------------------------------------------------
// Compare the forces from a single system and one member of a multi-system representation.
// Return the maximum absolute difference in Cartesian force components on any atom.
//
// Arguments:
//   ps:     A single system's coordinates, velocities, and forces
//   psy:    The multi-system object
//   sysid:  System ID within the multi-system object
//-------------------------------------------------------------------------------------------------
double getForceDeviation(const PhaseSpace &ps, const PhaseSpaceSynthesis *psy, int sysid) {
  const CoordinateFrame sys_frc = psy->exportCoordinates(sysid, TrajectoryKind::FORCES);
  const std::vector<double> frc_a = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> frc_b = sys_frc.getInterlacedCoordinates();
  return maxAbsoluteDifference(frc_a, frc_b);
}

//-------------------------------------------------------------------------------------------------
// Evaluate basic force field terms and CHARMM force field extensions evaluated with a PhaseSpace-
// and AtomGraph-Synthesis representation of many systems.  Check energies and forces reported
// against a simpler calculation of the corresponding quantity.
//
// Arguments:
//   poly_ag:   The synthesis of topologies
//   poly_ps:   The synthesis of coordinate / velocity / force representations
//   do_tests:  Indication that tests are possible
//-------------------------------------------------------------------------------------------------
void checkSynthesis(const AtomGraphSynthesis &poly_ag, const StaticExclusionMaskSynthesis &syse,
                    PhaseSpaceSynthesis *poly_ps, const TestPriority do_tests,
                    const EvaluateNonbonded do_nonbonded = EvaluateNonbonded::NO) {
  
  // Get the valence abstract and prepare for energy calculations
  SyValenceKit<double> syvk = poly_ag.getDoublePrecisionValenceKit();
  SyAtomUpdateKit<double, double2, double4> syauk = poly_ag.getDoublePrecisionAtomUpdateKit();
  SyRestraintKit<double, double2, double4> syrk = poly_ag.getDoublePrecisionRestraintKit();
  ScoreCard sc(poly_ps->getSystemCount(), 1, 32);

  // Bonds
  poly_ps->initializeForces();
  evalSyValenceEnergy<double, double2, double4>(syvk, syauk, syrk, poly_ps->data(), &sc,
                                                EvaluateForce::YES, VwuTask::BOND,
                                                VwuGoal::ACCUMULATE, 0);
  const int nsys = poly_ps->getSystemCount();
  Approx error_limits(std::vector<double>(nsys, 0.0), ComparisonType::ABSOLUTE, verytiny);
  std::vector<double> bond_nrg, bond_nrg_answer, bond_frc_deviations;
  ScoreCard tmp_sc(1, 1, 32);    
  for (int i = 0; i < nsys; i++) {
    bond_nrg.push_back(sc.reportInstantaneousStates(StateVariable::BOND, i));
    PhaseSpace psi = poly_ps->exportSystem(i);
    psi.initializeForces();
    bond_nrg_answer.push_back(evaluateBondTerms(poly_ag.getSystemTopologyPointer(i), &psi, &tmp_sc,
                                                EvaluateForce::YES));
    bond_frc_deviations.push_back(getForceDeviation(psi, poly_ps, i));
  }
  error_limits.setValues(std::vector<double>(nsys, 1.0e-5));
  check(bond_nrg, RelationalOperator::EQUAL, Approx(bond_nrg_answer).margin(5.5e-7), "Bond "
        "energies computed using the synthesis methods are inconsistent with those computed using "
        "a simpler approach.", do_tests);
  check(bond_frc_deviations, RelationalOperator::LESS_THAN, error_limits, "Forces due to bond "
        "interactions are inconsistent with those computed using a simpler approach.", do_tests);

  // Typical harmonic angles
  poly_ps->initializeForces();
  evalSyValenceEnergy<double, double2, double4>(syvk, syauk, syrk, poly_ps->data(), &sc,
                                                EvaluateForce::YES, VwuTask::ANGL,
                                                VwuGoal::ACCUMULATE, 0);
  std::vector<double> angl_nrg, angl_nrg_answer, angl_frc_deviations;
  for (int i = 0; i < nsys; i++) {
    angl_nrg.push_back(sc.reportInstantaneousStates(StateVariable::ANGLE, i));
    PhaseSpace psi = poly_ps->exportSystem(i);
    psi.initializeForces();
    angl_nrg_answer.push_back(evaluateAngleTerms(poly_ag.getSystemTopologyPointer(i), &psi,
                                                 &tmp_sc, EvaluateForce::YES));
    angl_frc_deviations.push_back(getForceDeviation(psi, poly_ps, i));
  }
  error_limits.setValues(std::vector<double>(nsys, 5.0e-6));
  check(angl_frc_deviations, RelationalOperator::LESS_THAN, error_limits, "Forces due to harmonic "
        "angle interactions are inconsistent with those computed using a simpler approach.",
        do_tests);
  check(angl_nrg, RelationalOperator::EQUAL, Approx(angl_nrg_answer).margin(3.1e-7),
        "Harmonic angle energies computed using the synthesis methods are inconsistent with those "
        "computed using a simpler approach.", do_tests);

  // Cosine-based dihedrals
  poly_ps->initializeForces();
  evalSyValenceEnergy<double, double2, double4>(syvk, syauk, syrk, poly_ps->data(), &sc,
                                                EvaluateForce::YES, VwuTask::DIHE,
                                                VwuGoal::ACCUMULATE, 0);
  std::vector<double> dihe_nrg, impr_nrg, dihe_nrg_answer, impr_nrg_answer, dihe_frc_deviations;
  for (int i = 0; i < nsys; i++) {
    dihe_nrg.push_back(sc.reportInstantaneousStates(StateVariable::PROPER_DIHEDRAL, i));
    impr_nrg.push_back(sc.reportInstantaneousStates(StateVariable::IMPROPER_DIHEDRAL, i));
    PhaseSpace psi = poly_ps->exportSystem(i);
    psi.initializeForces();
    const double2 du = evaluateDihedralTerms(poly_ag.getSystemTopologyPointer(i), &psi, &tmp_sc,
                                             EvaluateForce::YES);
    dihe_nrg_answer.push_back(du.x);
    impr_nrg_answer.push_back(du.y);
    dihe_frc_deviations.push_back(getForceDeviation(psi, poly_ps, i));
  }
  error_limits.setValues(std::vector<double>(nsys, 1.5e-6));
  check(dihe_frc_deviations, RelationalOperator::LESS_THAN, error_limits, "Forces due to "
        "cosine-based dihedral interactions are inconsistent with those computed using a simpler "
        "approach.", do_tests);
  check(dihe_nrg, RelationalOperator::EQUAL, Approx(dihe_nrg_answer).margin(1.0e-7),
        "Cosine-based dihedral energies computed using the synthesis methods are inconsistent "
        "with those computed using a simpler approach.", do_tests);
  check(impr_nrg, RelationalOperator::EQUAL, Approx(impr_nrg_answer).margin(4.0e-9),
        "Cosine-based improper dihedral energies computed using the synthesis methods are "
        "inconsistent with those computed using a simpler approach.", do_tests);

  // General 1:4 interactions
  poly_ps->initializeForces();
  evalSyValenceEnergy<double, double2, double4>(syvk, syauk, syrk, poly_ps->data(), &sc,
                                                EvaluateForce::YES, VwuTask::INFR14,
                                                VwuGoal::ACCUMULATE, 0);
  std::vector<double> qq14_nrg, lj14_nrg, qq14_nrg_answer, lj14_nrg_answer, attn14_frc_deviations;
  for (int i = 0; i < nsys; i++) {
    qq14_nrg.push_back(sc.reportInstantaneousStates(StateVariable::ELEC_ONE_FOUR, i));
    lj14_nrg.push_back(sc.reportInstantaneousStates(StateVariable::VDW_ONE_FOUR, i));
    PhaseSpace psi = poly_ps->exportSystem(i);
    psi.initializeForces();
    const double2 du = evaluateAttenuated14Terms(poly_ag.getSystemTopologyPointer(i), &psi,
                                                 &tmp_sc, EvaluateForce::YES, EvaluateForce::YES);
    qq14_nrg_answer.push_back(du.x);
    lj14_nrg_answer.push_back(du.y);
    attn14_frc_deviations.push_back(getForceDeviation(psi, poly_ps, i));
  }
  error_limits.setValues(std::vector<double>(nsys, 4.5e-6));
  check(attn14_frc_deviations, RelationalOperator::LESS_THAN, error_limits, "Forces due to "
        "attenuated 1:4 interactions are inconsistent with those computed using a simpler "
        "approach.", do_tests);
  check(qq14_nrg, RelationalOperator::EQUAL, Approx(qq14_nrg_answer).margin(1.5e-7),
        "Attenuated 1:4 electrostatic energies computed using the synthesis methods are "
        "inconsistent with those computed using a simpler approach.", do_tests);
  check(lj14_nrg, RelationalOperator::EQUAL, Approx(lj14_nrg_answer).margin(2.5e-7),
        "Attenuated 1:4 van-der Waals energies computed using the synthesis methods are "
        "inconsistent with those computed using a simpler approach.", do_tests);

  // Urey-Bradley interactions
  poly_ps->initializeForces();
  evalSyValenceEnergy<double, double2, double4>(syvk, syauk, syrk, poly_ps->data(), &sc,
                                                EvaluateForce::YES, VwuTask::UBRD,
                                                VwuGoal::ACCUMULATE, 0);
  std::vector<double> ubrd_nrg, ubrd_nrg_answer, ubrd_frc_deviations;
  for (int i = 0; i < nsys; i++) {
    ubrd_nrg.push_back(sc.reportInstantaneousStates(StateVariable::UREY_BRADLEY, i));
    PhaseSpace psi = poly_ps->exportSystem(i);
    psi.initializeForces();
    ubrd_nrg_answer.push_back(evaluateUreyBradleyTerms(poly_ag.getSystemTopologyPointer(i), &psi,
                                                       &tmp_sc, EvaluateForce::YES));
    ubrd_frc_deviations.push_back(getForceDeviation(psi, poly_ps, i));
  }
  error_limits.setValues(std::vector<double>(nsys, 4.5e-6));
  check(ubrd_frc_deviations, RelationalOperator::LESS_THAN, error_limits, "Forces due to "
        "Urey-Bradley interactions are inconsistent with those computed using a simpler approach.",
        do_tests);
  check(ubrd_nrg, RelationalOperator::EQUAL, Approx(ubrd_nrg_answer).margin(1.5e-7),
        "Urey-Bradley interaction energies computed using the synthesis methods are inconsistent "
        "with those computed using a simpler approach.", do_tests);
  
  // CHARMM improper dihedral interactions
  poly_ps->initializeForces();
  evalSyValenceEnergy<double, double2, double4>(syvk, syauk, syrk, poly_ps->data(), &sc,
                                                EvaluateForce::YES, VwuTask::CIMP,
                                                VwuGoal::ACCUMULATE, 0);
  std::vector<double> cimp_nrg, cimp_nrg_answer, cimp_frc_deviations;
  for (int i = 0; i < nsys; i++) {
    cimp_nrg.push_back(sc.reportInstantaneousStates(StateVariable::CHARMM_IMPROPER, i));
    PhaseSpace psi = poly_ps->exportSystem(i);
    psi.initializeForces();
    cimp_nrg_answer.push_back(evaluateCharmmImproperTerms(poly_ag.getSystemTopologyPointer(i),
                                                          &psi, &tmp_sc, EvaluateForce::YES));
    cimp_frc_deviations.push_back(getForceDeviation(psi, poly_ps, i));
  }
  error_limits.setValues(std::vector<double>(nsys, 4.5e-6));
  check(cimp_frc_deviations, RelationalOperator::LESS_THAN, error_limits, "Forces due to "
        "CHARMM improper interactions are inconsistent with those computed using a simpler "
        "approach.", do_tests);
  check(cimp_nrg, RelationalOperator::EQUAL, Approx(cimp_nrg_answer).margin(1.5e-7),
        "CHARMM improper interaction energies computed using the synthesis methods are "
        "inconsistent with those computed using a simpler approach.", do_tests);
  
  // CMAP interactions
  poly_ps->initializeForces();
  evalSyValenceEnergy<double, double2, double4>(syvk, syauk, syrk, poly_ps->data(), &sc,
                                                EvaluateForce::YES, VwuTask::CMAP,
                                                VwuGoal::ACCUMULATE, 0);
  std::vector<double> cmap_nrg, cmap_nrg_answer, cmap_frc_deviations;
  for (int i = 0; i < nsys; i++) {
    cmap_nrg.push_back(sc.reportInstantaneousStates(StateVariable::CMAP, i));
    PhaseSpace psi = poly_ps->exportSystem(i);
    psi.initializeForces();
    cmap_nrg_answer.push_back(evaluateCmapTerms(poly_ag.getSystemTopologyPointer(i), &psi, &tmp_sc,
                                                EvaluateForce::YES));
    cmap_frc_deviations.push_back(getForceDeviation(psi, poly_ps, i));
  }
  error_limits.setValues(std::vector<double>(nsys, 1.0e-6));
  check(cmap_frc_deviations, RelationalOperator::LESS_THAN, error_limits, "Forces due to CMAP "
        "energy surface contributions are inconsistent with those computed using a simpler "
        "approach.", do_tests);
  check(cmap_nrg, RelationalOperator::EQUAL, Approx(cmap_nrg_answer).margin(3.1e-7),
        "CMAP energies computed using the synthesis methods are inconsistent with those "
        "computed using a simpler approach.", do_tests);

  // Various restraints
  poly_ps->initializeForces();
  evalSyValenceEnergy<double, double2, double4>(syvk, syauk, syrk, poly_ps->data(), &sc,
                                                EvaluateForce::YES, VwuTask::RPOSN,
                                                VwuGoal::ACCUMULATE, 0);
  evalSyValenceEnergy<double, double2, double4>(syvk, syauk, syrk, poly_ps->data(), &sc,
                                                EvaluateForce::YES, VwuTask::RBOND,
                                                VwuGoal::ACCUMULATE, 0);
  evalSyValenceEnergy<double, double2, double4>(syvk, syauk, syrk, poly_ps->data(), &sc,
                                                EvaluateForce::YES, VwuTask::RANGL,
                                                VwuGoal::ACCUMULATE, 0);
  evalSyValenceEnergy<double, double2, double4>(syvk, syauk, syrk, poly_ps->data(), &sc,
                                                EvaluateForce::YES, VwuTask::RDIHE,
                                                VwuGoal::ACCUMULATE, 0);
  std::vector<double> rstr_nrg, rstr_nrg_answer, rstr_frc_deviations;
  for (int i = 0; i < nsys; i++) {
    rstr_nrg.push_back(sc.reportInstantaneousStates(StateVariable::RESTRAINT, i));
    PhaseSpace psi = poly_ps->exportSystem(i);
    psi.initializeForces();
    rstr_nrg_answer.push_back(evaluateRestraints(poly_ag.getSystemRestraintPointer(i), &psi,
                                                 &tmp_sc, EvaluateForce::YES));
    rstr_frc_deviations.push_back(getForceDeviation(psi, poly_ps, i));
  }
  error_limits.setValues(std::vector<double>(nsys, 1.0e-6));
  check(rstr_frc_deviations, RelationalOperator::LESS_THAN, error_limits, "Forces due to "
        "restraint contributions are inconsistent with those computed using a simpler "
        "approach.", do_tests);
  check(rstr_nrg, RelationalOperator::EQUAL, Approx(rstr_nrg_answer).margin(3.1e-7),
        "Restraint energy penalties computed using the synthesis methods are inconsistent with "
        "those computed using a simpler approach.", do_tests);

  // Non-bonded interactions
  if (do_nonbonded == EvaluateNonbonded::YES) {
    poly_ps->initializeForces();
    evalSyNonbondedEnergy<double>(poly_ag, syse, poly_ps, &sc, EvaluateForce::YES,
                                  EvaluateForce::YES);
  }
}

//-------------------------------------------------------------------------------------------------
// Check that the atomic partial charges obtained by looking up charge parameter indices and then
// reading from the array of unique partial charges match those obtained by reading directly from
// the array of atomic partial charges.  Accessing the index and then the parameter may seem like
// an extra step, but in some contexts the charge and Lennard-Jones parameter indices can be
// encoded in the same 32-bit word, saving memory bandwidth and conserving L1 cache.
//
// Arguments:
//   poly_ag:   The topology synthesis to check
//   do_tests:  Indicate whether the tests are possible absed on previous file reading
//-------------------------------------------------------------------------------------------------
void inspectChargeIndexing(const AtomGraphSynthesis &poly_ag, const TestPriority do_tests) {
  const SyNonbondedKit<double, double2> synbk = poly_ag.getDoublePrecisionNonbondedKit();
  int total_atoms = 0;
  for (int i = 0; i < synbk.nsys; i++) {
    total_atoms += synbk.atom_counts[i];
  }
  std::vector<double> indexed_q(total_atoms), stored_q(total_atoms);
  int pos = 0;
  for (int i = 0; i < synbk.nsys; i++) {
    const int jllim = synbk.atom_offsets[i];
    const int jhlim = jllim + synbk.atom_counts[i];
    for (int j = jllim; j < jhlim; j++) {
      indexed_q[pos] = synbk.q_params[synbk.q_idx[j]];
      stored_q[pos]  = synbk.charge[j];
    }
  }
  check(indexed_q, RelationalOperator::EQUAL, stored_q, "Charges obtained by lookup in the "
        "parameter array do not agree with those obtained by direct reference to the unrolled "
        "array.", do_tests);
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Obtain environment variables or command-line input, if available
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  StopWatch timer;

  // Section 1
  section("Test AtomGraphSynthesis layout");

  // Section 2
  section("Test molecular mechanics potential calculations");

  // Section 3
  section("Apply implicit solvent models");
  
  // Section 4
  section("Traps for bad input");
  
  // Create some vectors of random numbers, then upload them and test what happens when perturbing
  // atomic coordinates by these numbers.
  const char osc = osSeparator();
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string tip3p_top_name = base_top_name + osc + "tip3p.top";
  const std::string tip4p_top_name = base_top_name + osc + "tip4p.top";
  const std::string trpcage_top_name = base_top_name + osc + "trpcage_in_water.top";
  const std::string trpcage_nbfix_top_name = base_top_name + osc + "trpcage_in_water_nbfix.top";
  const std::string ubiquitin_top_name = base_top_name + osc + "ubiquitin.top";
  const std::string drug_top_name = base_top_name + osc + "drug_example_vs.top";
  const std::string brbz_top_name = base_top_name + osc + "bromobenzene_vs.top";
  const bool files_exist = (getDrivePathType(tip3p_top_name) == DrivePathType::FILE &&
                            getDrivePathType(tip4p_top_name) == DrivePathType::FILE &&
                            getDrivePathType(trpcage_top_name) == DrivePathType::FILE &&
                            getDrivePathType(trpcage_nbfix_top_name) == DrivePathType::FILE &&
                            getDrivePathType(ubiquitin_top_name) == DrivePathType::FILE &&
                            getDrivePathType(drug_top_name) == DrivePathType::FILE &&
                            getDrivePathType(brbz_top_name) == DrivePathType::FILE);
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  AtomGraph tip3p_ag, tip4p_ag, trpcage_ag, trpcage2_ag, trpcage3_ag, nbfix_ag, ubiquitin_ag;
  AtomGraph drug_ag, brbz_ag;
  if (files_exist) {
    tip3p_ag.buildFromPrmtop(tip3p_top_name);
    tip4p_ag.buildFromPrmtop(tip4p_top_name);
    trpcage_ag.buildFromPrmtop(trpcage_top_name);
    trpcage2_ag.buildFromPrmtop(trpcage_top_name);
    trpcage3_ag.buildFromPrmtop(trpcage_top_name);
    nbfix_ag.buildFromPrmtop(trpcage_nbfix_top_name);
    ubiquitin_ag.buildFromPrmtop(ubiquitin_top_name);
    drug_ag.buildFromPrmtop(drug_top_name);
    brbz_ag.buildFromPrmtop(brbz_top_name);
  }
  else {
    rtWarn("The topology files for the TIP3P and TIP4P water boxes as well as two versions of the "
           "solvated Trp-cage miniprotein, ubiquitin, and two drug molecules must be available in "
           "${STORMM_SOURCE}/test/ subdirectories Topology and Trajectory, respectively.  Check "
           "the $STORMM_SOURCE environment variable to make sure that it is set properly.  A "
           "number of tests will be skipped.", "test_atomgraph_synthesis");
  }

  // Set one of the Trp-cage systems to have a different topology source name.  This will trick
  // the synthesis generator into thinking that this is a unique system, whereas the third copy
  // is not, leaving six (not seven) total topologies in the synthesis.
  trpcage2_ag.setSource(trpcage_top_name + "_2");

  // Create the synthesis
  const std::vector<AtomGraph*> all_tops = { &tip3p_ag, &tip4p_ag, &trpcage_ag, &trpcage2_ag,
                                             &trpcage3_ag, &nbfix_ag, &ubiquitin_ag, &drug_ag,
                                             &brbz_ag};
  const std::vector<int> system_ids = { 0, 1, 2, 3, 4, 3, 3, 5, 2, 1, 1, 3, 6, 7, 8 };
  AtomGraphSynthesis poly_ag(all_tops, system_ids, ExceptionResponse::SILENT,
                             null_gpu, &timer);
  const StaticExclusionMaskSynthesis poly_se(poly_ag.getUniqueTopologies(),
                                             poly_ag.getTopologyIndices());
  poly_ag.loadNonbondedWorkUnits(poly_se);

  // Create a copy of the synthesis
  AtomGraphSynthesis poly_ag_copy(poly_ag);
  const StaticExclusionMaskSynthesis poly_se_copy(poly_ag_copy.getUniqueTopologies(),
                                                  poly_ag_copy.getTopologyIndices());
  poly_ag_copy.loadNonbondedWorkUnits(poly_se_copy);

  // Check various descriptors
  section(1);
  check(poly_ag.getAtomCount(), RelationalOperator::EQUAL, 50570, "The topology synthesis does "
        "not contain the expected number of atoms.", do_tests);
  check(poly_ag.getVirtualSiteCount(), RelationalOperator::EQUAL, 1264, "The topology synthesis "
        "does not contain the expected number of virtual sites.", do_tests);
  std::vector<int> valence_term_counts;
  if (files_exist) {
    valence_term_counts = { poly_ag.getBondTermCount(), poly_ag.getAngleTermCount(),
                            poly_ag.getDihedralTermCount(), poly_ag.getUreyBradleyTermCount(),
                            poly_ag.getCharmmImproperTermCount(), poly_ag.getCmapTermCount() };
  }
  else {
    valence_term_counts = { 0, 0, 0, 0, 0, 0 };
  }
  const std::vector<int> valence_term_answer = { 50608, 6823, 16765, 0, 0, 71 };
  check(valence_term_counts, RelationalOperator::EQUAL, valence_term_answer, "The topology "
        "synthesis contains incorrect numbers of some valence terms.", do_tests);

  // Get the coordinates for all structures
  section(2);
  PhaseSpace tip3p_ps, tip4p_ps, trpcage_ps, trpcage2_ps, trpcage3_ps, nbfix_ps, ubiquitin_ps;
  PhaseSpace drug_ps, brbz_ps;
  const std::string base_crd_name  = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string tip3p_crd_name     = base_crd_name + osc + "tip3p.inpcrd";
  const std::string tip4p_crd_name     = base_crd_name + osc + "tip4p.inpcrd";
  const std::string trpcage_crd_name   = base_crd_name + osc + "trpcage_in_water.inpcrd";
  const std::string ubiquitin_crd_name = base_crd_name + osc + "ubiquitin.inpcrd";
  const std::string drug_crd_name      = base_crd_name + osc + "drug_example_vs.inpcrd";
  const std::string brbz_crd_name      = base_crd_name + osc + "bromobenzene_vs.inpcrd";
  const bool coords_exist = (getDrivePathType(tip3p_crd_name) == DrivePathType::FILE &&
                             getDrivePathType(tip4p_crd_name) == DrivePathType::FILE &&
                             getDrivePathType(trpcage_crd_name) == DrivePathType::FILE &&
                             getDrivePathType(ubiquitin_crd_name) == DrivePathType::FILE &&
                             getDrivePathType(drug_crd_name) == DrivePathType::FILE &&
                             getDrivePathType(brbz_crd_name) == DrivePathType::FILE);
  if (coords_exist) {
    tip3p_ps.buildFromFile(tip3p_crd_name);
    tip4p_ps.buildFromFile(tip4p_crd_name);
    trpcage_ps.buildFromFile(trpcage_crd_name);
    ubiquitin_ps.buildFromFile(ubiquitin_crd_name);
    drug_ps.buildFromFile(drug_crd_name);
    brbz_ps.buildFromFile(brbz_crd_name);
  }
  else {
    rtWarn("Coordinates for the periodic systems needed to accompany the first AtomGraphSynthesis "
           "were not found.  Check the installation and the ${STORMM_SOURCE} environment "
           "variable.  Subsequent tests will be skipped.\n");
  }
  const TestPriority do_per_eval = (coords_exist && files_exist) ? TestPriority::CRITICAL :
                                                                   TestPriority::ABORT;
  
  const std::string base_pept_top_name = oe.getStormmSourcePath() + osc + "test" + osc +
                                         "Namelist" + osc + "topol";
  const std::string base_pept_crd_name = oe.getStormmSourcePath() + osc + "test" + osc +
                                         "Namelist" + osc + "coord";
  const std::string brbz_vs_top_name = base_top_name + osc + "bromobenzene_vs.top";
  const std::string brbz_vs_crd_name = base_crd_name + osc + "bromobenzene_vs.inpcrd";
  const std::vector<AtomGraph*> ag_list = poly_ag.getSystemTopologyPointer();
  std::vector<PhaseSpace> ps_list;
  ps_list.reserve(poly_ag.getSystemCount());
  for (int i = 0; i < poly_ag.getSystemCount(); i++) {
    const AtomGraph *ag_ptr = poly_ag.getSystemTopologyPointer(i);
    if (ag_ptr == &tip3p_ag) {
      ps_list.push_back(tip3p_ps);
    }
    else if (ag_ptr == &tip4p_ag) {
      ps_list.push_back(tip4p_ps);
    }
    else if (ag_ptr == &trpcage_ag || ag_ptr == &trpcage2_ag || ag_ptr == &trpcage3_ag ||
             ag_ptr == &nbfix_ag) {
      ps_list.push_back(trpcage_ps);
    }
    else if (ag_ptr == &ubiquitin_ag) {
      ps_list.push_back(ubiquitin_ps);
    }
    else if (ag_ptr == &drug_ag) {
      ps_list.push_back(drug_ps);
    }
    else if (ag_ptr == &brbz_ag) {
      ps_list.push_back(brbz_ps);
    }
  }
  PhaseSpaceSynthesis poly_ps(ps_list, ag_list);
  checkSynthesis(poly_ag, poly_se, &poly_ps, do_tests, EvaluateNonbonded::NO);
  checkSynthesis(poly_ag_copy, poly_se_copy, &poly_ps, do_tests, EvaluateNonbonded::NO);

  // Check the charge indexing.
  inspectChargeIndexing(poly_ag_copy, do_per_eval);

  // Prepare some more systems
  const std::string tiso_top_name = base_top_name + osc + "trpcage.top";
  const std::string tiso_crd_name = base_crd_name + osc + "trpcage.inpcrd";
  const std::string brbi_top_name = base_top_name + osc + "bromobenzene_vs_iso.top";
  const std::string brbi_crd_name = base_crd_name + osc + "bromobenzene_vs_iso.inpcrd";
  const std::string lig1_top_name = base_top_name + osc + "stereo_L1_vs.top";
  const std::string lig1_crd_name = base_crd_name + osc + "stereo_L1_vs.inpcrd";
  const std::string lig2_top_name = base_top_name + osc + "symmetry_L1_vs.top";
  const std::string lig2_crd_name = base_crd_name + osc + "symmetry_L1_vs.inpcrd";
  const std::string dhfr_top_name = base_top_name + osc + "dhfr_cmap.top";
  const std::string dhfr_crd_name = base_crd_name + osc + "dhfr_cmap.inpcrd";
  const bool new_exist = (getDrivePathType(tiso_top_name) == DrivePathType::FILE &&
                          getDrivePathType(brbi_top_name) == DrivePathType::FILE &&
                          getDrivePathType(lig1_top_name) == DrivePathType::FILE &&
                          getDrivePathType(lig2_top_name) == DrivePathType::FILE &&
                          getDrivePathType(dhfr_top_name) == DrivePathType::FILE &&
                          getDrivePathType(tiso_crd_name) == DrivePathType::FILE &&
                          getDrivePathType(brbi_crd_name) == DrivePathType::FILE &&
                          getDrivePathType(lig1_crd_name) == DrivePathType::FILE &&
                          getDrivePathType(lig2_crd_name) == DrivePathType::FILE &&
                          getDrivePathType(dhfr_crd_name) == DrivePathType::FILE);
  AtomGraph tiso_ag, brbi_ag, lig1_ag, lig2_ag, dhfr_ag;
  PhaseSpace tiso_ps, brbi_ps, lig1_ps, lig2_ps, dhfr_ps;
  if (new_exist) {
    tiso_ag.buildFromPrmtop(tiso_top_name, ExceptionResponse::SILENT);
    brbi_ag.buildFromPrmtop(brbi_top_name, ExceptionResponse::SILENT);
    lig1_ag.buildFromPrmtop(lig1_top_name, ExceptionResponse::SILENT);
    lig2_ag.buildFromPrmtop(lig2_top_name, ExceptionResponse::SILENT);
    dhfr_ag.buildFromPrmtop(dhfr_top_name, ExceptionResponse::SILENT);
    tiso_ps.buildFromFile(tiso_crd_name);
    brbi_ps.buildFromFile(brbi_crd_name);
    lig1_ps.buildFromFile(lig1_crd_name);
    lig2_ps.buildFromFile(lig2_crd_name);
    dhfr_ps.buildFromFile(dhfr_crd_name);
  }
  else {
    rtWarn("Files corresponding to various systems in isolated boundary conditions were not "
           "found.  Check the ${STORMM_SOURCE} environment variable.  The necessary directories "
           "are the same as for other files needed by this test program.  Subsequent tests will "
           "be skipped.", "test_atomgraph_synthesis");
  }
  const TestPriority do_new_tests = (new_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;

  // Create some restraints and apply them, then check the synthesis implementation
  RestraintApparatus tiso_ra = assembleRestraints(&tiso_ag, tiso_ps);
  RestraintApparatus brbi_ra = assembleRestraints(&brbi_ag, brbi_ps);
  RestraintApparatus lig1_ra = assembleRestraints(&lig1_ag, lig1_ps);
  RestraintApparatus lig2_ra = assembleRestraints(&lig2_ag, lig2_ps);
  RestraintApparatus dhfr_ra = assembleRestraints(&dhfr_ag, dhfr_ps);
  std::vector<AtomGraph*> agn_list = { &tiso_ag, &brbi_ag, &lig1_ag, &lig2_ag, &dhfr_ag };
  std::vector<RestraintApparatus*> rsn_list = { &tiso_ra, &brbi_ra, &lig1_ra, &lig2_ra, &dhfr_ra };
  std::vector<PhaseSpace> psn_list = { tiso_ps, brbi_ps, lig1_ps, lig2_ps, dhfr_ps };
  std::vector<int> system_list = { 0, 1, 2, 3, 4, 4, 4, 4, 4, 4 };
  AtomGraphSynthesis poly_agn_rst(agn_list, rsn_list, system_list, system_list,
                                  ExceptionResponse::SILENT, null_gpu, &timer);
  const StaticExclusionMaskSynthesis poly_sen(poly_agn_rst.getUniqueTopologies(),
                                              poly_agn_rst.getTopologyIndices());
  poly_agn_rst.loadNonbondedWorkUnits(poly_sen);
  PhaseSpaceSynthesis poly_psn(psn_list, agn_list, system_list);
  checkSynthesis(poly_agn_rst, poly_sen, &poly_psn, do_new_tests, EvaluateNonbonded::YES);

  // Create a copy of the new topology synthesis (including restraints), and check whether
  // an exclusion mask object created for the first topology synthesis (which should be equivalent)
  // works in conjunction with the new synthesis.
  AtomGraphSynthesis poly_agn_rst_copy(poly_agn_rst);
  poly_agn_rst_copy.loadNonbondedWorkUnits(poly_sen);
  checkSynthesis(poly_agn_rst_copy, poly_sen, &poly_psn, do_new_tests, EvaluateNonbonded::YES);

  // Check the charge indexing.
  inspectChargeIndexing(poly_agn_rst, do_new_tests);
  
  // Apply implicit solvent models to the synthesis
  NeckGeneralizedBornTable ngb_tab;
  poly_agn_rst.setImplicitSolventModel(ImplicitSolventModel::NECK_GB_II, ngb_tab,
                                       AtomicRadiusSet::MBONDI3);
  const int nsys = poly_agn_rst.getSystemCount();
  std::vector<int> radius_mistakes(nsys, 0);
  std::vector<int> sp_radius_mistakes(nsys, 0);
  std::vector<int> abg_mistakes(nsys, 0);
  std::vector<int> sp_abg_mistakes(nsys, 0);
  const SyNonbondedKit<double, double2> ism_nbk = poly_agn_rst.getDoublePrecisionNonbondedKit();
  const SyNonbondedKit<float, float2> ism_nbk_sp = poly_agn_rst.getSinglePrecisionNonbondedKit();
  for (int i = 0; i < nsys; i++) {
    const AtomGraph *iag_ptr = poly_agn_rst.getSystemTopologyPointer(i);
    const ImplicitSolventKit<double> isk = iag_ptr->getDoublePrecisionImplicitSolventKit();
    const ImplicitSolventKit<float> isk_sp = iag_ptr->getSinglePrecisionImplicitSolventKit();
    const int natom = iag_ptr->getAtomCount();
    const int aoffs = ism_nbk.atom_offsets[i];
    for (int j = 0; j < natom; j++) {
      if (fabs(isk.pb_radii[j] - ism_nbk.pb_radii[aoffs + j]) > stormm::constants::tiny) {
        radius_mistakes[i] += 1;
      }
      if (fabs(isk_sp.pb_radii[j] - ism_nbk_sp.pb_radii[aoffs + j]) > stormm::constants::tiny) {
        sp_radius_mistakes[i] += 1;
      }
      if (fabs(isk.gb_alpha[j] - ism_nbk.gb_alpha[aoffs + j]) > stormm::constants::tiny ||
          fabs(isk.gb_beta[j]  - ism_nbk.gb_beta[aoffs + j])  > stormm::constants::tiny ||
          fabs(isk.gb_gamma[j] - ism_nbk.gb_gamma[aoffs + j]) > stormm::constants::tiny) {
        abg_mistakes[i] += 1;
      }
      if (fabs(isk_sp.gb_alpha[j] - ism_nbk_sp.gb_alpha[aoffs + j]) > stormm::constants::tiny ||
          fabs(isk_sp.gb_beta[j]  - ism_nbk_sp.gb_beta[aoffs + j])  > stormm::constants::tiny ||
          fabs(isk_sp.gb_gamma[j] - ism_nbk_sp.gb_gamma[aoffs + j]) > stormm::constants::tiny) {
        sp_abg_mistakes[i] += 1;
      }
    }
  }
  check(radius_mistakes, RelationalOperator::EQUAL, std::vector<int>(nsys, 0), "Radii entered "
        "into the synthesis when setting all systems to MBondi3 disagree with the underlying "
        "topologies.  Precision setting: " + getEnumerationName(PrecisionModel::DOUBLE) + ".",
        do_new_tests);
  check(sp_radius_mistakes, RelationalOperator::EQUAL, std::vector<int>(nsys, 0), "Radii entered "
        "into the synthesis when setting all systems to MBondi3 disagree with the underlying "
        "topologies.  Precision setting: " + getEnumerationName(PrecisionModel::SINGLE) + ".",
        do_new_tests);
  check(abg_mistakes, RelationalOperator::EQUAL, std::vector<int>(nsys, 0), "Alpha, Beta, and "
        "Gamma atomic parameters entered into the synthesis when setting all systems to MBondi3 "
        "disagree with the underlying topologies.  Precision setting: " +
        getEnumerationName(PrecisionModel::DOUBLE) + ".", do_new_tests);
  check(sp_abg_mistakes, RelationalOperator::EQUAL, std::vector<int>(nsys, 0), "Alpha, Beta, and "
        "Gamma atomic parameters entered into the synthesis when setting all systems to MBondi3 "
        "disagree with the underlying topologies.  Precision setting: " +
        getEnumerationName(PrecisionModel::SINGLE) + ".", do_new_tests);

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
