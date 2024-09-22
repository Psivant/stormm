#include "copyright.h"
#include "../../src/Chemistry/chemical_features.h"
#include "../../src/Constants/behavior.h"
#include "../../src/FileManagement/file_util.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/MolecularMechanics/minimization.h"
#include "../../src/MolecularMechanics/mm_evaluation.h"
#include "../../src/Namelists/nml_minimize.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Potential/static_exclusionmask.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Restraints/bounded_restraint.h"
#include "../../src/Restraints/restraint_builder.h"
#include "../../src/Restraints/restraint_apparatus.h"
#include "../../src/Structure/virtual_site_handling.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/file_snapshot.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::chemistry;
using namespace stormm::diskutil;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::mm;
using namespace stormm::namelist;
using namespace stormm::random;
using namespace stormm::restraints;
using namespace stormm::review;
using namespace stormm::structure;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Test that each atom of the molecule is in a local minimum with respect to displacement in the
// x, y, and z directions.
//
// Arguments:
//   ps:             Coordinates, box binformation, and forces
//   ag:             System topology
//   ra:             Restraint collection acting on the system, as a supplement to the topology
//   se:             Permanent map of the non-bonded exclusions for all pairs of atoms
//   final_step_no:  Number of the final minimization step (in the event that time- or
//                   step-dependent restraints are in effect)
//   grad_tol:       Tolerance for the gradient at any one atom, in kcal/mol-A
//   do_test:        Indicator of whether the test is feasible
//-------------------------------------------------------------------------------------------------
void testLocalMinimum(PhaseSpace *ps, const AtomGraph &ag, const RestraintApparatus &ra,
                      const StaticExclusionMask &se, const int final_step_no,
                      const double grad_tol, const TestPriority do_test) {
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  PhaseSpaceWriter psw = ps->data();
  std::vector<int> local_minimum(3 * cdk.natom, 1);
  const double test_displacement = 0.0008;
  ScoreCard lsc(1, 16, 32);
  evalNonbValeRestMM(ps, &lsc, ag, se, ra, EvaluateForce::NO, 0, final_step_no);
  const double e0 = lsc.reportTotalEnergy();
  for (int i = 0; i < cdk.natom; i++) {

    // Skip extra points
    if (cdk.z_numbers[i] == 0) {
      continue;
    }

    // Perturb the system along the X direction
    psw.xcrd[i] += test_displacement;
    placeVirtualSites(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, vsk);
    lsc.initialize();
    evalNonbValeRestMM(ps, &lsc, ag, se, ra, EvaluateForce::NO, 0, final_step_no);
    const double epx = lsc.reportTotalEnergy();
    psw.xcrd[i] -= 2.0 * test_displacement;
    placeVirtualSites(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, vsk);    
    lsc.initialize();
    evalNonbValeRestMM(ps, &lsc, ag, se, ra, EvaluateForce::NO, 0, final_step_no);
    const double enx = lsc.reportTotalEnergy();
    psw.xcrd[i] += test_displacement;

    // Perturb the system along the Y direction
    psw.ycrd[i] += test_displacement;
    placeVirtualSites(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, vsk);    
    lsc.initialize();
    evalNonbValeRestMM(ps, &lsc, ag, se, ra, EvaluateForce::NO, 0, final_step_no);
    const double epy = lsc.reportTotalEnergy();
    psw.ycrd[i] -= 2.0 * test_displacement;
    placeVirtualSites(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, vsk);    
    lsc.initialize();
    evalNonbValeRestMM(ps, &lsc, ag, se, ra, EvaluateForce::NO, 0, final_step_no);
    const double eny = lsc.reportTotalEnergy();
    psw.ycrd[i] += test_displacement;

    // Perturb the system along the Z direction
    psw.zcrd[i] += test_displacement;
    placeVirtualSites(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, vsk);    
    lsc.initialize();
    evalNonbValeRestMM(ps, &lsc, ag, se, ra, EvaluateForce::NO, 0, final_step_no);
    const double epz = lsc.reportTotalEnergy();
    psw.zcrd[i] -= 2.0 * test_displacement;
    placeVirtualSites(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu, psw.unit_cell, vsk);    
    lsc.initialize();
    evalNonbValeRestMM(ps, &lsc, ag, se, ra, EvaluateForce::NO, 0, final_step_no);
    const double enz = lsc.reportTotalEnergy();
    psw.zcrd[i] += test_displacement;

    const double xgrad = (epx - enx) / (2.0 * test_displacement);
    const double ygrad = (epy - eny) / (2.0 * test_displacement);
    const double zgrad = (epz - enz) / (2.0 * test_displacement);
    
    // Check that the original position is lower in energy
    local_minimum[(3 * i)    ] = (fabs(xgrad) < grad_tol);
    local_minimum[(3 * i) + 1] = (fabs(ygrad) < grad_tol);
    local_minimum[(3 * i) + 2] = (fabs(zgrad) < grad_tol);
  }
  check(local_minimum, RelationalOperator::EQUAL, std::vector<int>(3 * cdk.natom, 1),
        "The system described by topology " + getBaseName(ag.getFileName()) + " was not left in a "
        "local minimum by energy minimization.", do_test);
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv, TmpdirStatus::REQUIRED);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  StopWatch timer;

  // Section 1
  section("Minimize a dipeptide");

  // Read files
  const char osc = osSeparator();
  const std::string base_top_name  = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name  = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string alad_top_name = base_top_name + osc + "ala_dipeptide.top";
  const std::string alad_crd_name = base_crd_name + osc + "ala_dipeptide.inpcrd";
  const std::string brbz_top_name = base_top_name + osc + "bromobenzene.top";
  const std::string brbz_crd_name = base_crd_name + osc + "bromobenzene.inpcrd";
  const std::string bbvs_top_name = base_top_name + osc + "bromobenzene_vs.top";
  const std::string bbvs_crd_name = base_crd_name + osc + "bromobenzene_vs.inpcrd";
  const std::vector<std::string> all_top = { alad_top_name, brbz_top_name, bbvs_top_name };
  const std::vector<std::string> all_crd = { alad_crd_name, brbz_crd_name, bbvs_crd_name };
  const int system_count = all_top.size();
  bool files_exist = true;
  for (int i = 0; i < system_count; i++) {
    files_exist = (getDrivePathType(all_top[i]) == DrivePathType::FILE &&
                   getDrivePathType(all_crd[i]) == DrivePathType::FILE && files_exist);
  }
  if (files_exist == false) {
    rtWarn("Topology and coordinate files must be available for tests to proceed.  Check the "
           "${STORMM_SOURCE} environment variable, currently set to " + oe.getStormmSourcePath() +
           ", to verify that it is valid.  Subsequent tests will be skipped.",
           "test_minimization");
  }
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  std::vector<AtomGraph> all_ag;
  std::vector<PhaseSpace> all_ps;
  std::vector<ChemicalFeatures> all_chemfe;
  std::vector<RestraintApparatus> all_ra;
  std::vector<StaticExclusionMask> all_se;
  all_ag.reserve(system_count);
  all_ps.reserve(system_count);
  all_chemfe.reserve(system_count);
  all_ra.reserve(system_count);
  all_se.reserve(system_count);
  if (files_exist) {
    for (int i = 0; i < system_count; i++) {
      all_ag.emplace_back(all_top[i], ExceptionResponse::SILENT);
      all_ps.emplace_back(all_crd[i]);
      all_chemfe.emplace_back(&all_ag[i], all_ps[i], MapRotatableGroups::YES, 300.0, &timer);
      const std::vector<BoundedRestraint> brs = applyHydrogenBondPreventors(&all_ag[i],
                                                                            all_chemfe[i], 16.0); 
      all_ra.emplace_back(brs, &all_ag[i]);
      all_se.emplace_back(&all_ag[i]);
    }
  }
  else {
    all_ag.resize(system_count);
    all_ps.resize(system_count);
    all_chemfe.resize(system_count);
    all_ra.resize(system_count);
    all_se.resize(system_count);
  }
  
  // Try the dipeptide--this systems contains CMAPs in addition to basic Amber force field terms
  section(1);
  MinimizeControls mincon;
  mincon.setTotalCycles(1000);
  mincon.setClashDampingCycles(0);
  Xoroshiro128pGenerator xrs(53018479);
  if (files_exist) {
    for (size_t i = 0; i < all_ps.size(); i++) {
      PhaseSpaceWriter psw = all_ps[i].data();
      addRandomNoise(&xrs, psw.xcrd, psw.ycrd, psw.zcrd, psw.natom, 0.1, 1.0);
    }
    timer.assignTime(0);
    const int alad_timings = timer.addCategory("Minimize Ala dipeptide");
    const ScoreCard alad_ene = minimize(&all_ps[0], all_ag[0], all_ra[0], all_se[0], mincon);
    timer.assignTime(alad_timings);
    const int brbz_timings = timer.addCategory("Minimize Bromobenzene");
    const ScoreCard brbz_ene = minimize(&all_ps[1], all_ag[1], all_ra[1], all_se[1], mincon);
    timer.assignTime(brbz_timings);
    const int bbvs_timings = timer.addCategory("Minimize Bromobenzene-VS");
    const ScoreCard bbvs_ene = minimize(&all_ps[2], all_ag[2], all_ra[2], all_se[2], mincon);
    timer.assignTime(bbvs_timings);
  }
  for (int i = 0; i < system_count; i++) {
    testLocalMinimum(&all_ps[i], all_ag[i], all_ra[i], all_se[i], mincon.getTotalCycles(),
                     1.0e-2, do_tests);
  }
  
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
