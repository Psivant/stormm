#include <algorithm>
#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Chemistry/atommask.h"
#include "../../src/Chemistry/chemical_features.h"
#include "../../src/Chemistry/chemistry_enumerators.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/mixed_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/reduction.h"
#include "../../src/Math/reduction_abstracts.h"
#include "../../src/Math/reduction_workunit.h"
#include "../../src/Math/rounding.h"
#include "../../src/Math/sorting_enumerators.h"
#include "../../src/Math/summation.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Namelists/nml_files.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/static_exclusionmask.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Potential/eval_valence_workunit.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Restraints/restraint_apparatus.h"
#include "../../src/Restraints/restraint_builder.h"
#include "../../src/Structure/local_arrangement.h"
#include "../../src/Synthesis/brickwork.h"
#include "../../src/Synthesis/nonbonded_workunit.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/static_mask_synthesis.h"
#include "../../src/Synthesis/synthesis_abstracts.h"
#include "../../src/Synthesis/synthesis_cache_map.h"
#include "../../src/Synthesis/systemcache.h"
#include "../../src/Synthesis/valence_workunit.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/trajectory_enumerators.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::chemistry::AtomMask;
using stormm::chemistry::ChemicalFeatures;
using stormm::chemistry::MapRotatableGroups;
using stormm::constants::ExceptionResponse;
using stormm::constants::verytiny;
using stormm::constants::tiny;
using stormm::constants::warp_size_int;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getBaseName;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::openOutputFile;
using stormm::diskutil::osSeparator;
using stormm::errors::rtWarn;
using stormm::namelist::FilesControls;
using stormm::parse::findStringInVector;
using stormm::parse::TextFile;
using stormm::parse::TextOrigin;
using stormm::random::Xoroshiro128pGenerator;
using stormm::random::addRandomNoise;
using stormm::restraints::applyHydrogenBondPreventors;
using stormm::restraints::applyPositionalRestraints;
using stormm::restraints::BoundedRestraint;
using stormm::restraints::RestraintApparatus;
using stormm::restraints::RestraintKit;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using stormm::structure::distance;
using stormm::structure::angle;
using stormm::structure::dihedralAngle;
using namespace stormm::card;
using namespace stormm::data_types;
using namespace stormm::energy;
using namespace stormm::stmath;
using namespace stormm::numerics;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Check the coverage of a simple task form a ValenceWorkUnit.
//
// Arguments:
//   accd:              Accumulator directives: coverage is only counted if the accumulator bit is
//                      set to 1, which means energy is accumulated by the work unit
//   taskid:            List of simple task topological index numbers
//   coverage:          Coverage array (accumulated and returned)
//   range_problem_in:  Current state of the array range indexing problem detector
//-------------------------------------------------------------------------------------------------
bool checkSimpleTaskCoverage(const std::vector<uint> &accd, const std::vector<int> &taskid,
                             std::vector<int> *coverage, const bool range_problem_in) {
  bool range_problem = range_problem_in;
  const int n_items = coverage->size();
  int* cov_ptr = coverage->data();
  for (size_t j = 0; j < taskid.size(); j++) {
    if (readBitFromMask(accd, j) == 0) {
      continue;
    }
    if (taskid[j] >= 0 && taskid[j] < n_items) {
      cov_ptr[taskid[j]] += 1;
    }
    else {
      range_problem = true;
    }
  }
  return range_problem;
}

//-------------------------------------------------------------------------------------------------
// Check that the naive distance between two particles is the distance computed upon re-imaging.
// Return TRUE is this is so, FALSE otherwise.
//
// Arguments:
//   i:     The first particle
//   j:     The second particle
//   cfr:   Coordinates of the particles, plus box information
//-------------------------------------------------------------------------------------------------
bool checkNaiveDistance(const int i, const int j, const CoordinateFrameReader &cfr) {
  const double reim_dist  = distance(i, j, cfr);
  const double naive_dist = distance<double, double>(i, j, cfr.xcrd, cfr.ycrd, cfr.zcrd, nullptr,
                                                     nullptr, UnitCellType::NONE);
  return (fabs(reim_dist - naive_dist) < stormm::constants::tiny);
}

//-------------------------------------------------------------------------------------------------
// Run a series of tests using valence work units.
//
// Arguments:
//   top_name:  Name of the topology to use
//   crd_name:  Name of the coordinate file to use
//   oe:        Operating environment information (for error reporting)
//   my_prng:   Random number generator (modified by use inside this function)
//-------------------------------------------------------------------------------------------------
void runValenceWorkUnitTests(const std::string &top_name, const std::string &crd_name,
                             const TestEnvironment &oe, Xoroshiro128pGenerator *my_prng) {
  const bool files_exist = (getDrivePathType(top_name) == DrivePathType::FILE &&
                            getDrivePathType(crd_name) == DrivePathType::FILE);
  if (files_exist == false) {
    rtWarn("The topology and input coordinates for a critical system appear to be missing.  Check "
           "the ${STORMM_SOURCE} variable (currently " + oe.getStormmSourcePath() + ") to make "
           "sure that " + top_name + " and " + crd_name + " valid paths.", "test_synthesis");
  }
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  AtomGraph ag  = (files_exist) ? AtomGraph(top_name, ExceptionResponse::SILENT) :
                                  AtomGraph();
  PhaseSpace ps = (files_exist) ? PhaseSpace(crd_name) : PhaseSpace();
  const CoordinateFrameReader cfr(ps);

  // Create a set of restraints on the structure, starting from backbone positional restraints
  // and adding things that randomly keep various two-, three-, and four-point meansurements near
  // their observed values.
  std::vector<BoundedRestraint> mol_rstr = applyPositionalRestraints(&ag, cfr, "@N,CA,C,O", 1.4);
  for (size_t i = 0; i < mol_rstr.size(); i++) {
    const double3 trgt = mol_rstr[i].getTargetSite();
    mol_rstr[i].setTargetSite({ trgt.x + 0.5 - my_prng->uniformRandomNumber(),
                                 trgt.y + 0.5 - my_prng->uniformRandomNumber(),
                                 trgt.z + 0.5 - my_prng->uniformRandomNumber() });
  }
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const int solute_extent = ag.getLastSoluteAtom();
  const int min_iinc = std::max(solute_extent / 8, 4);
  for (int i = min_iinc; i < solute_extent; i += min_iinc) {
    for (int j = 0; j < i; j += min_iinc) {
      const double orig_distance = distance(i, j, ps) + 0.5 - my_prng->uniformRandomNumber();
      if (checkNaiveDistance(i, j, cfr) == false) {
        continue;
      }
      mol_rstr.emplace_back(i, j, &ag, 1.0, 1.4, 0.0, orig_distance - 0.25, orig_distance + 0.25,
                            orig_distance + 100.0);
      if (checkNaiveDistance(i, i - 2, cfr) && checkNaiveDistance(j, i - 2, cfr)) {
        const double orig_angle = angle(i, j, i - 2, ps) +
                                  (0.1 * (0.5 - my_prng->uniformRandomNumber()));
        if (orig_angle > 0.1 && orig_angle < 3.0) {
          mol_rstr.emplace_back(i, j, i - 2, &ag, 2.3, 0.9, 0.0, orig_angle - 0.05,
                                orig_angle + 0.02, stormm::symbols::pi - stormm::constants::tiny);
        }
      }
      if (checkNaiveDistance(i, i - 1, cfr) && checkNaiveDistance(j, i - 1, cfr) &&
          checkNaiveDistance(j, j + 1, cfr)) {
        const double orig_dihedral = dihedralAngle(i, i - 1, j, j + 1, ps) +
                                     (0.8 * (0.5 - my_prng->uniformRandomNumber()));
        mol_rstr.emplace_back(i, i - 1, j, j + 1, &ag, 2.3, 0.9, orig_dihedral - 1.0,
                              orig_dihedral - 0.03, orig_dihedral + 0.03, orig_dihedral + 1.0);
      }
    }
  }

  // Create the restraint apparatus, then the valence work units.
  RestraintApparatus ra(mol_rstr);
  ValenceDelegator vdel(&ag, &ra);
  const std::vector<ValenceWorkUnit> all_vwu = buildValenceWorkUnits(&vdel);
  const int n_vwu = all_vwu.size();

  // Check the coverage in valence work units with an independent tally of each atom and term
  std::vector<int> bond_coverage(vk.nbond, 0);
  std::vector<int> angl_coverage(vk.nangl, 0);
  std::vector<int> dihe_coverage(vk.ndihe, 0);
  std::vector<int> ubrd_coverage(vk.nubrd, 0);
  std::vector<int> cimp_coverage(vk.ncimp, 0);
  std::vector<int> cmap_coverage(vk.ncmap, 0);
  std::vector<int> infr_coverage(vk.ninfr14, 0);
  bool bond_range_problem = false;
  bool angl_range_problem = false;
  bool dihe_range_problem = false;
  bool ubrd_range_problem = false;
  bool cimp_range_problem = false;
  bool cmap_range_problem = false;
  bool infr_range_problem = false;
  for (int i = 0; i < n_vwu; i++) {
    const std::vector<uint2> cbnd_insr = all_vwu[i].getCompositeBondInstructions();
    const std::vector<uint> cbnd_accd  = all_vwu[i].getAccumulationFlags(VwuTask::CBND);
    const std::vector<int> cbnd_taskid = all_vwu[i].getCompositeBondTaskList();
    for (size_t j = 0; j < cbnd_insr.size(); j++) {
      if (readBitFromMask(cbnd_accd, j) == 0) {
        continue;
      }
      if ((cbnd_insr[j].x >> 20) & 0x1) {
        if (cbnd_taskid[j] >= 0 && cbnd_taskid[j] < vk.nubrd) {
          ubrd_coverage[cbnd_taskid[j]] += 1;
        }
        else {
          ubrd_range_problem = true;
        }
      }
      else {
        if (cbnd_taskid[j] >= 0 && cbnd_taskid[j] < vk.nbond) {
          bond_coverage[cbnd_taskid[j]] += 1;
        }
        else {
          bond_range_problem = true;
        }
      }
    }
    angl_range_problem = checkSimpleTaskCoverage(all_vwu[i].getAccumulationFlags(VwuTask::ANGL),
                                                 all_vwu[i].getSimpleTaskList(VwuTask::ANGL),
                                                 &angl_coverage, angl_range_problem);
    const std::vector<uint3> cdhe_insr = all_vwu[i].getCompositeDihedralInstructions();
    const std::vector<uint> cdhe_accd  = all_vwu[i].getAccumulationFlags(VwuTask::CDHE);
    const std::vector<int2> cdhe_taskid = all_vwu[i].getCompositeDihedralTaskList();
    for (size_t j = 0; j < cdhe_insr.size(); j++) {
      if (readBitFromMask(cdhe_accd, j) == 0) {
        continue;
      }
      if ((cdhe_insr[j].x >> 30) & 0x1) {
        if (cdhe_taskid[j].x >= 0 && cdhe_taskid[j].x < vk.ncimp) {
          cimp_coverage[cdhe_taskid[j].x] += 1;
        }
        else {
          cimp_range_problem = true;
        }
      }
      else {
        if (cdhe_taskid[j].x >= 0 && cdhe_taskid[j].x < vk.ndihe) {
          dihe_coverage[cdhe_taskid[j].x] += 1;
        }
        else {
          dihe_range_problem = true;          
        }
        if (cdhe_taskid[j].y >= 0) {
          if (cdhe_taskid[j].y < vk.ndihe) {
            dihe_coverage[cdhe_taskid[j].y] += 1;
          }
          else {
            dihe_range_problem = true;          
          }
        }
      }
    }
    cmap_range_problem = checkSimpleTaskCoverage(all_vwu[i].getAccumulationFlags(VwuTask::CMAP),
                                                 all_vwu[i].getSimpleTaskList(VwuTask::CMAP),
                                                 &cmap_coverage, cmap_range_problem);
    infr_range_problem = checkSimpleTaskCoverage(all_vwu[i].getAccumulationFlags(VwuTask::INFR14),
                                                 all_vwu[i].getSimpleTaskList(VwuTask::INFR14),
                                                 &infr_coverage, infr_range_problem);
  }
  const std::vector<int> bond_coverage_answer(vk.nbond, 1);
  const std::vector<int> angl_coverage_answer(vk.nangl, 1);
  const std::vector<int> dihe_coverage_answer(vk.ndihe, 1);
  const std::vector<int> ubrd_coverage_answer(vk.nubrd, 1);
  const std::vector<int> cimp_coverage_answer(vk.ncimp, 1);
  const std::vector<int> cmap_coverage_answer(vk.ncmap, 1);
  const std::vector<int> infr_coverage_answer(vk.ninfr14, 1);
  if (vk.nbond > 0) {
    check(bond_range_problem == false, "Composite bond instructions reference a bad topology bond "
          "index in valence work units for topology " + ag.getFileName() + ".", do_tests);
    check(bond_coverage, RelationalOperator::EQUAL, bond_coverage_answer, "Bond accumulation is "
          "incorrect in valence work units for topology " + ag.getFileName() + ".", do_tests);
  }
  if (vk.nangl > 0) {
    check(angl_range_problem == false, "Composite angle instructions reference a bad topology "
          "index in valence work units for topology " + ag.getFileName() + ".", do_tests);
    check(angl_coverage, RelationalOperator::EQUAL, angl_coverage_answer, "Angle accumulation is "
          "incorrect in valence work units for topology " + ag.getFileName() + ".", do_tests);
  }
  if (vk.ndihe > 0) {
    check(dihe_range_problem == false, "Composite dihedral instructions reference a bad topology "
          "dihedral index in valence work units for topology " + ag.getFileName() + ".", do_tests);
    check(dihe_coverage, RelationalOperator::EQUAL, dihe_coverage_answer, "Dihedral accumulation "
          "is incorrect in valence work units for topology " + ag.getFileName() + ".", do_tests);
  }
  if (vk.nubrd > 0) {
    check(ubrd_range_problem == false, "Composite bond instructions reference a bad topology "
          "Urey-Bradley index in valence work units for topology " + ag.getFileName() + ".",
          do_tests);
    check(ubrd_coverage, RelationalOperator::EQUAL, ubrd_coverage_answer, "Urey-Bradley "
          "accumulation is incorrect in valence work units for topology " + ag.getFileName() + ".",
          do_tests);
  }
  if (vk.ncimp > 0) {
    check(cimp_range_problem == false, "Composite dihedral instructions reference a bad CHARMM "
          "improper topology index in valence work units for topology " + ag.getFileName() + ".",
          do_tests);
    check(cimp_coverage, RelationalOperator::EQUAL, cimp_coverage_answer, "CHARMM improper "
          "accumulation is incorrect in valence work units for topology " + ag.getFileName() + ".",
          do_tests);
  }
  if (vk.ncmap > 0) {
    check(cmap_range_problem == false, "CMAP instructions reference a bad topology index in "
          "valence work units for topology " + ag.getFileName() + ".", do_tests);
    check(cmap_coverage, RelationalOperator::EQUAL, cmap_coverage_answer, "CMAP surface term "
          "accumulation is incorrect in valence work units for topology " + ag.getFileName() + ".",
          do_tests);
  }
  if (vk.ninfr14 > 0) {
    check(infr_range_problem == false, "Inferred 1:4 instructions reference a bad topology index "
          "in valence work units for topology " + ag.getFileName() + ".", do_tests);
    check(infr_coverage, RelationalOperator::EQUAL, infr_coverage_answer, "Inferred 1:4 "
          "non-bonded interaction accumulation is incorrect in valence work units for topology " +
          ag.getFileName() + ".", do_tests);
  }

  // Check each component of the energies the work units evaluate
  ScoreCard sc(2);
  PhaseSpace ps_vwu(ps);
  evaluateBondTerms(&ag, &ps, &sc, EvaluateForce::YES, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::BOND);
  const std::vector<double> bond_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> bond_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  ps.initializeForces();
  ps_vwu.initializeForces();
  evaluateAngleTerms(&ag, &ps, &sc, EvaluateForce::YES, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::ANGL);
  const std::vector<double> angl_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> angl_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  ps.initializeForces();
  ps_vwu.initializeForces();
  evaluateDihedralTerms(&ag, &ps, &sc, EvaluateForce::YES, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::DIHE);
  const std::vector<double> dihe_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> dihe_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  ps.initializeForces();
  ps_vwu.initializeForces();
  evaluateUreyBradleyTerms(&ag, &ps, &sc, EvaluateForce::YES, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::UBRD);
  const std::vector<double> ubrd_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> ubrd_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  ps.initializeForces();
  ps_vwu.initializeForces();
  evaluateCharmmImproperTerms(&ag, &ps, &sc, EvaluateForce::YES, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::CIMP);
  const std::vector<double> cimp_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> cimp_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  ps.initializeForces();
  ps_vwu.initializeForces();
  evaluateCmapTerms(&ag, &ps, &sc, EvaluateForce::YES, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::CMAP);
  const std::vector<double> cmap_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> cmap_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  ps.initializeForces();
  ps_vwu.initializeForces();
  evaluateAttenuated14Terms(&ag, &ps, &sc, EvaluateForce::YES, EvaluateForce::YES, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::INFR14);
  const std::vector<double> attn_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> attn_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  ps.initializeForces();
  ps_vwu.initializeForces();
  
  // Restraint energies need to be plucked from the results immediately, as the reference
  // calculations only evaluate restraints all at once whereas the ValenceWorkUnit evaluation
  // implements them one type at a time.
  evaluateRestraints(&ra, &ps, &sc, EvaluateForce::YES, 0);
  const double rest_e = sc.reportInstantaneousStates(StateVariable::RESTRAINT, 0);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::RPOSN);
  const double rposn_e = sc.reportInstantaneousStates(StateVariable::RESTRAINT, 1);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::RBOND);
  const double rbond_e = sc.reportInstantaneousStates(StateVariable::RESTRAINT, 1);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::RANGL);
  const double rangl_e = sc.reportInstantaneousStates(StateVariable::RESTRAINT, 1);
  evalValenceWorkUnits(&ag, &ps_vwu, &ra, &sc, 1, all_vwu, EvaluateForce::YES, VwuTask::RDIHE);
  const double rdihe_e = sc.reportInstantaneousStates(StateVariable::RESTRAINT, 1);
  const std::vector<double> rstr_frc = ps_vwu.getInterlacedCoordinates(TrajectoryKind::FORCES);
  const std::vector<double> rstr_frc_ref = ps.getInterlacedCoordinates(TrajectoryKind::FORCES);
  
  // Read off other energy terms
  const std::vector<double> bond_e = sc.reportInstantaneousStates(StateVariable::BOND);
  const std::vector<double> angl_e = sc.reportInstantaneousStates(StateVariable::ANGLE);
  const std::vector<double> dihe_e = sc.reportInstantaneousStates(StateVariable::PROPER_DIHEDRAL);
  const std::vector<double> impr_e =
    sc.reportInstantaneousStates(StateVariable::IMPROPER_DIHEDRAL);
  const std::vector<double> ubrd_e = sc.reportInstantaneousStates(StateVariable::UREY_BRADLEY);
  const std::vector<double> cimp_e = sc.reportInstantaneousStates(StateVariable::CHARMM_IMPROPER);
  const std::vector<double> cmap_e = sc.reportInstantaneousStates(StateVariable::CMAP);
  const std::vector<double> qq14_e =
    sc.reportInstantaneousStates(StateVariable::ELEC_ONE_FOUR);
  const std::vector<double> lj14_e =
    sc.reportInstantaneousStates(StateVariable::VDW_ONE_FOUR);
  
  // Check that energy terms computed with ValenceWorkUnits match the reference calculations
  check(bond_e[0], RelationalOperator::EQUAL, Approx(bond_e[1]).margin(1.0e-6), "Bond energies "
        "computed with ValenceWorkUnits do not agree with the reference calculations (topology " +
        top_name + ").", do_tests);
  check(angl_e[0], RelationalOperator::EQUAL, Approx(angl_e[1]).margin(1.0e-6), "Angle energies "
        "computed with ValenceWorkUnits do not agree with the reference calculations (topology " +
        top_name + ").", do_tests);
  check(dihe_e[0], RelationalOperator::EQUAL, Approx(dihe_e[1]).margin(1.0e-6), "Cosine-based "
        "proper dihedral energies computed with ValenceWorkUnits do not agree with the reference "
        "calculations (topology " + top_name + ").", do_tests);
  check(impr_e[0], RelationalOperator::EQUAL, Approx(impr_e[1]).margin(1.0e-6), "Cosine-based "
        "improper dihedral energies computed with ValenceWorkUnits do not agree with the "
        "reference calculations (topology " + top_name + ").", do_tests);
  check(ubrd_e[0], RelationalOperator::EQUAL, Approx(ubrd_e[1]).margin(1.0e-6), "Urey-Bradley "
        "energies computed with ValenceWorkUnits do not agree with the reference calculations "
        "(topology " + top_name + ").", do_tests);
  check(cimp_e[0], RelationalOperator::EQUAL, Approx(cimp_e[1]).margin(1.0e-6), "CHARMM improper "
        "energies computed with ValenceWorkUnits do not agree with the reference calculations "
        "(topology " + top_name + ").", do_tests);
  check(cmap_e[0], RelationalOperator::EQUAL, Approx(cmap_e[1]).margin(1.0e-6), "CMAP energies "
        "computed with ValenceWorkUnits do not agree with the reference calculations (topology " +
        top_name + ").", do_tests);
  check(qq14_e[0], RelationalOperator::EQUAL, Approx(qq14_e[1]).margin(1.0e-6), "Electrostatic "
        "1:4 interactions computed with ValenceWorkUnits do not agree with the reference "
        "calculations (topology " + top_name + ").", do_tests);
  check(lj14_e[0], RelationalOperator::EQUAL, Approx(lj14_e[1]).margin(1.0e-6), "Electrostatic "
        "1:4 interactions computed with ValenceWorkUnits do not agree with the reference "
        "calculations (topology " + top_name + ").", do_tests);
  check(rest_e, RelationalOperator::EQUAL, Approx(rposn_e + rbond_e + rangl_e + rdihe_e, 1.0e-6),
        "Restraint penalty energies computed with ValenceWorkUnits do not agree with the "
        "reference calculations (topology " + top_name + ").", do_tests);

  // Check that forces computed with ValenceWorkUnits match the reference calculations
  check(bond_frc, RelationalOperator::EQUAL, Approx(bond_frc_ref).margin(1.0e-6), "Forces due to "
        "bonds computed with ValenceWorkUnits do not agree with the reference calculations "
        "(topology " + top_name + ").", do_tests);
  check(angl_frc, RelationalOperator::EQUAL, Approx(angl_frc_ref).margin(1.0e-6), "Forces due to "
        "angles computed with ValenceWorkUnits do not agree with the reference calculations "
        "(topology " + top_name + ").", do_tests);
  check(dihe_frc, RelationalOperator::EQUAL, Approx(dihe_frc_ref).margin(1.0e-6), "Forces due to "
        "proper and improper (cosine-based) dihedral terms computed with ValenceWorkUnits do not "
        "agree with the reference calculations (topology " + top_name + ").", do_tests);
  check(attn_frc, RelationalOperator::EQUAL, Approx(attn_frc_ref).margin(1.0e-6), "Forces due to "
        "attenuated 1:4 interactions computed with ValenceWorkUnits do not agree with the "
        "reference calculations (topology " + top_name + ").", do_tests);
  check(ubrd_frc, RelationalOperator::EQUAL, Approx(ubrd_frc_ref).margin(1.0e-6), "Forces due to "
        "Urey-Bradley terms computed with ValenceWorkUnits do not agree with the reference "
        "calculations (topology " + top_name + ").", do_tests);
  check(cimp_frc, RelationalOperator::EQUAL, Approx(cimp_frc_ref).margin(1.0e-6), "Forces due to "
        "CHARMM improper dihedrals computed with ValenceWorkUnits do not agree with the reference "
        "calculations (topology " + top_name + ").", do_tests);
  check(cmap_frc, RelationalOperator::EQUAL, Approx(cmap_frc_ref).margin(1.0e-6), "Forces due to "
        "CMAP potential surfaces computed with ValenceWorkUnits do not agree with the reference "
        "calculations (topology " + top_name + ").", do_tests);
  check(rstr_frc, RelationalOperator::EQUAL, Approx(rstr_frc_ref).margin(1.0e-6), "Forces due to "
        "restraint potentials computed with ValenceWorkUnits do not agree with the reference "
        "calculations (topology " + top_name + ").", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Accumulate contributing atoms base on a series of force terms.
//
// Arguments:
//   item_count:   The number of energy terms of some type, the trusted length of {i,j,...}_atoms
//   result:       Array of all atoms that influence any one atom (accumulated and returned)
//   i_atoms:      The first list of atoms associated with whatever energy term
//   j_atoms:      The second list of atoms associated with whatever energy term
//   k_atoms:      The third list of atoms associated with whatever energy term
//   l_atoms:      The fourth list of atoms associated with whatever energy term
//   m_atoms:      The fifth list of atoms associated with whatever energy term
//-------------------------------------------------------------------------------------------------
void accumulateContributingAtoms(const int item_count, std::vector<std::vector<int>> *result,
                                 const int* i_atoms, const int* j_atoms = nullptr,
                                 const int* k_atoms = nullptr, const int* l_atoms = nullptr,
                                 const int* m_atoms = nullptr) {
  std::vector<int>* result_data = result->data();
  for (int pos = 0; pos < item_count; pos++) {
    std::vector<int> atom_list;
    atom_list.push_back(i_atoms[pos]);
    int order = 1;
    if (j_atoms != nullptr) {
      atom_list.push_back(j_atoms[pos]);
      order++;
    }
    if (k_atoms != nullptr) {
      atom_list.push_back(k_atoms[pos]);
      order++;
    }
    if (l_atoms != nullptr) {
      atom_list.push_back(l_atoms[pos]);
      order++;
    }
    if (m_atoms != nullptr) {
      atom_list.push_back(m_atoms[pos]);
      order++;
    }
    for (int i = 0; i < order - 1; i++) {
      for (int j = i + 1; j < order; j++) {
        result_data[atom_list[i]].push_back(atom_list[j]);
        result_data[atom_list[j]].push_back(atom_list[i]);
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Determine a list of all atoms that contribute to forces on any particular atom through valence
// terms, restraints, or virtual sites.  This is done by looping over all topology terms and
// restraints while extending lists pertaining to each atom of the topology.
//
// Arguments:
//   ag:      The system topology
//   ra:      The system restraint apparatus, a sort of supplemental energy surface
//-------------------------------------------------------------------------------------------------
std::vector<std::vector<int>> getAtomForceContributors(const AtomGraph &ag,
                                                       const RestraintApparatus &ra) {
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  const RestraintKit<double, double2, double4> rar = ra.dpData();

  // Initialize the result.  Every atom list includes the atom itself.
  std::vector<std::vector<int>> result(vk.natom, std::vector<int>());
  for (int i = 0; i < vk.natom; i++) {
    result[i].push_back(i);
  }
  accumulateContributingAtoms(vk.nbond, &result, vk.bond_i_atoms, vk.bond_j_atoms);
  accumulateContributingAtoms(vk.nangl, &result, vk.angl_i_atoms, vk.angl_j_atoms,
                              vk.angl_k_atoms);
  accumulateContributingAtoms(vk.ndihe, &result, vk.dihe_i_atoms, vk.dihe_j_atoms, vk.dihe_k_atoms,
                              vk.dihe_l_atoms);
  accumulateContributingAtoms(vk.nubrd, &result, vk.ubrd_i_atoms, vk.ubrd_k_atoms);
  accumulateContributingAtoms(vk.ncimp, &result, vk.cimp_i_atoms, vk.cimp_j_atoms, vk.cimp_k_atoms,
                              vk.cimp_l_atoms);
  accumulateContributingAtoms(vk.ncmap, &result, vk.cmap_i_atoms, vk.cmap_j_atoms, vk.cmap_k_atoms,
                              vk.cmap_l_atoms, vk.cmap_m_atoms);
  accumulateContributingAtoms(vk.ninfr14, &result, vk.infr14_i_atoms, vk.infr14_l_atoms);
  accumulateContributingAtoms(rar.nposn, &result, rar.rposn_atoms);
  accumulateContributingAtoms(rar.nbond, &result, rar.rbond_i_atoms, rar.rbond_j_atoms);
  accumulateContributingAtoms(rar.nangl, &result, rar.rangl_i_atoms, rar.rangl_j_atoms,
                              rar.rangl_k_atoms);
  accumulateContributingAtoms(rar.ndihe, &result, rar.rdihe_i_atoms, rar.rdihe_j_atoms,
                              rar.rdihe_k_atoms, rar.rdihe_l_atoms);

  // All atoms of each virtual site frame should be included with every atom of the frame
  for (int pos = 0; pos < vsk.nsite; pos++) {
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[pos])) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
      accumulateContributingAtoms(1, &result, &vsk.frame1_idx[pos], &vsk.frame2_idx[pos]);
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      accumulateContributingAtoms(1, &result, &vsk.frame1_idx[pos], &vsk.frame2_idx[pos],
                                  &vsk.frame3_idx[pos]);
      break;
    case VirtualSiteKind::FIXED_4:
      accumulateContributingAtoms(1, &result, &vsk.frame1_idx[pos], &vsk.frame2_idx[pos],
                                  &vsk.frame3_idx[pos], &vsk.frame4_idx[pos]);
      break;
    case VirtualSiteKind::NONE:
      break;
    }
  }

  // Finally, all atoms that contribute to forces on any virtual site should be included with
  // every frame atom.
  for (int pos = 0; pos < vsk.nsite; pos++) {
    const int vatom = vsk.vs_atoms[pos];
    result[vsk.frame1_idx[pos]].insert(result[vsk.frame1_idx[pos]].end(),
                                       result[vatom].begin(), result[vatom].end());
    result[vsk.frame2_idx[pos]].insert(result[vsk.frame2_idx[pos]].end(),
                                       result[vatom].begin(), result[vatom].end());
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[pos])) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      result[vsk.frame3_idx[pos]].insert(result[vsk.frame3_idx[pos]].end(),
                                         result[vatom].begin(), result[vatom].end());
      break;
    case VirtualSiteKind::FIXED_4:
      result[vsk.frame3_idx[pos]].insert(result[vsk.frame3_idx[pos]].end(),
                                         result[vatom].begin(), result[vatom].end());
      result[vsk.frame4_idx[pos]].insert(result[vsk.frame4_idx[pos]].end(),
                                         result[vatom].begin(), result[vatom].end());
      break;
    case VirtualSiteKind::NONE:
      break;
    }
  }

  // Reduce each vector to unique values and return the result
  for (int i = 0; i < vk.natom; i++) {
    reduceUniqueValues(&result[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Check that the non-bonded work units for a system cover all tiles along the diagonal and below
// it, no more or less.
//
// Arguments:
//   nbwu_vec:  Vector of non-bonded work units
//   ag:        Topology for the system of interest (this contains the atom count and, if needed,
//              the file name for error reporting purposes)
//   do_tests:  Indication that tests should be performed at all
//-------------------------------------------------------------------------------------------------
void checkNonbondedWorkUnitCoverage(const std::vector<NonbondedWorkUnit> &nbwu_vec,
                                    const AtomGraph &ag, const TestPriority do_tests) {
  const int natom = ag.getAtomCount();
  const int ntile = (natom + tile_length - 1) / tile_length;
  std::vector<bool> covered(ntile * ntile, false);
  const size_t nnbwu = nbwu_vec.size();
  int nduplicate = 0;
  for (size_t wu_idx = 0; wu_idx < nnbwu; wu_idx++) {
    const std::vector<int> abstract = nbwu_vec[wu_idx].getAbstract();
    std::vector<int> imports(abstract[0]);
    for (int i = 0; i < abstract[0]; i++) {
      imports[i] = abstract[i + 1] / tile_length;
    }
    const std::vector<uint2> tinsr = nbwu_vec[wu_idx].getTileInstructions();

    // Translate the instruction back to the absolute tile numbers in the system (the absolute
    // tile number of any atom is the absolute topological atom index number divided by the tile
    // length).
    const size_t insr_count = tinsr.size();
    for (size_t i = 0; i < insr_count; i++) {
      const int absc_side = imports[(tinsr[i].x & 0xffff) / tile_length];
      const int ordi_side = imports[((tinsr[i].x >> 16) & 0xffff) / tile_length];
      const int ao_idx = (absc_side * ntile) + ordi_side;
      const int oa_idx = (ordi_side * ntile) + absc_side;
      if (covered[ao_idx] || covered[oa_idx]) {
        nduplicate++;
      }
      covered[ao_idx] = true;
      covered[oa_idx] = true;
    }
  }
  int nmissing = 0;
  for (int i = 0; i < ntile; i++) {
    for (int j = 0; j <= i; j++) {
      nmissing += (covered[(j * ntile) + i] == false);
    }
  }
  check(nduplicate, RelationalOperator::EQUAL, 0, "Duplicate non-bonded work unit tiles were "
        "discovered in the topology for " + getBaseName(ag.getFileName()) + ", containing " +
        std::to_string(natom) + " atoms.  Some interactions would be double-counted.", do_tests);
  check(nmissing, RelationalOperator::EQUAL, 0, "A total of " + std::to_string(nmissing) +
        " non-bonded tiles were missing in the work units for the topology described in " +
        getBaseName(ag.getFileName()) + ", containing " + std::to_string(natom) + " atoms.",
        do_tests);
}

//-------------------------------------------------------------------------------------------------
// Run a mock reduction by construction work units and then checking the results against a direct,
// single-threaded approach over the same series of "systems."
//
// Arguments:
//   atom_counts:   Sizes of each system in the collection / mock synthesis
//   atom_offsets:  Starting indices of each system in the concatenated list
//   prop_x:        Vector containing individual "atom" readouts, the first of up to three
//                  properties to reduce across each system
//   prop_y:        [Optional] Vector containing the second property to reduce across each system
//   prop_z:        [Optional] Vector containing the third property to reduce across each system
//-------------------------------------------------------------------------------------------------
void testReduction(const std::vector<int> &atom_counts, const std::vector<int> &atom_offsets,
                   const std::vector<double> &prop_x, const std::vector<double> &prop_y = {},
                   const std::vector<double> &prop_z = {}) {
  const bool have_y = (prop_y.size() == prop_x.size());
  const bool have_z = (prop_z.size() == prop_x.size());
  
  // The mock GPU will have 84 streaming multiprocessors.
  GpuDetails gpu;
  gpu.setSMPCount(84);
  const int nprop = 1 + (prop_y.size() == prop_x.size()) + (prop_z.size() == prop_x.size());
  const std::vector<ReductionWorkUnit> rwu_vec = buildReductionWorkUnits(atom_offsets, atom_counts,
                                                                         gpu, nprop);
  const size_t nsys = atom_counts.size();
  std::vector<double> chksum_x(nsys, 0.0), chksum_y(nsys, 0.0), chksum_z(nsys, 0.0);
  int max_stash = 0;
  const size_t nrwu = rwu_vec.size();
  for (size_t i = 0; i < nrwu; i++) {
    max_stash = std::max(max_stash, rwu_vec[i].getResultIndex() + 1);
  }
  std::vector<double> rwusum_x(max_stash, 0.0), rwusum_y(max_stash, 0.0), rwusum_z(max_stash, 0.0);
  for (size_t i = 0; i < nsys; i++) {
    chksum_x[i] = sum<double>(&prop_x[atom_offsets[i]], atom_counts[i]);
    if (have_y) {
      chksum_y[i] = sum<double>(&prop_y[atom_offsets[i]], atom_counts[i]);
    }
    if (have_z) {
      chksum_z[i] = sum<double>(&prop_z[atom_offsets[i]], atom_counts[i]);
    }
  }
  for (size_t i = 0; i < nrwu; i++) {
    const int min_wu_idx = rwu_vec[i].getAtomStart();
    const int max_wu_idx = rwu_vec[i].getAtomEnd();
    const int wu_res_idx = rwu_vec[i].getResultIndex();
    rwusum_x[wu_res_idx] = sum<double>(&prop_x[min_wu_idx], max_wu_idx - min_wu_idx);
    if (have_y) {
      rwusum_y[wu_res_idx] = sum<double>(&prop_y[min_wu_idx], max_wu_idx - min_wu_idx);
    }
    if (have_z) {
      rwusum_z[wu_res_idx] = sum<double>(&prop_z[min_wu_idx], max_wu_idx - min_wu_idx);
    }
  }
  std::vector<int2> dep_chk(nsys);
  std::vector<bool> dep_found(nsys, false);
  std::vector<double> rwu_result_x(nsys), rwu_result_y(nsys), rwu_result_z(nsys);
  bool dep_mismatch = false;
  for (size_t i = 0; i < nrwu; i++) {
    const int dep_start = rwu_vec[i].getDependencyStart();
    const int dep_end   = rwu_vec[i].getDependencyEnd();
    const int sys_idx   = rwu_vec[i].getSystemIndex();
    if (dep_found[sys_idx]) {
      dep_mismatch = (dep_mismatch ||
                      dep_start != dep_chk[sys_idx].x || dep_end != dep_chk[sys_idx].y);
    }
    else {
      dep_found[sys_idx] = true;
      dep_chk[sys_idx].x = dep_start;
      dep_chk[sys_idx].y = dep_end;
    }
  }
  for (size_t i = 0; i < nsys; i++) {
    rwu_result_x[i] = sum<double>(&rwusum_x[dep_chk[i].x], dep_chk[i].y - dep_chk[i].x);
    if (have_y) {
      rwu_result_y[i] = sum<double>(&rwusum_y[dep_chk[i].x], dep_chk[i].y - dep_chk[i].x);
    }
    if (have_z) {
      rwu_result_z[i] = sum<double>(&rwusum_z[dep_chk[i].x], dep_chk[i].y - dep_chk[i].x);
    }
  }
  check(dep_mismatch == false, "Dependencies for various reduction work units are inconsistent "
        "across work units serving the same system in a collection with between " +
        std::to_string(minValue(atom_counts)) + " and " + std::to_string(maxValue(atom_counts)) +
        " items per system.");
  check(rwu_result_x, RelationalOperator::EQUAL, Approx(chksum_x).margin(tiny), "Reduction work "
        "units did not properly guide the summation of X values in a collection with between " +
        std::to_string(minValue(atom_counts)) + " and " + std::to_string(maxValue(atom_counts)) +
        " items per system.");
  if (have_y) {
    check(rwu_result_y, RelationalOperator::EQUAL, Approx(chksum_y).margin(tiny), "Reduction work "
          "units did not properly guide the summation of Y values in a collection with between " +
          std::to_string(minValue(atom_counts)) + " and " + std::to_string(maxValue(atom_counts)) +
          " items per system.");
  }
  if (have_z) {
    check(rwu_result_z, RelationalOperator::EQUAL, Approx(chksum_z).margin(tiny), "Reduction work "
          "units did not properly guide the summation of Z values in a collection with between " +
          std::to_string(minValue(atom_counts)) + " and " + std::to_string(maxValue(atom_counts)) +
          " items per system.");
  }

  // Make a more detailed mockup of the reduction work units, closer in implementation to what
  // will happen on the GPU.
  std::vector<int> rwu_abstracts(nrwu * rdwu_abstract_length);
  RdwuPerSystem rps = RdwuPerSystem::ONE;
  for (int i = 0; i < nrwu; i++) {
    const std::vector<int> tr_abs = rwu_vec[i].getAbstract();
    for (int j = 0; j < rdwu_abstract_length; j++) {
      rwu_abstracts[(rdwu_abstract_length * i) + j] = tr_abs[j];
    }
    if (rwu_vec[i].getDependencyEnd() - rwu_vec[i].getDependencyStart() > 1) {
      rps = RdwuPerSystem::MULTIPLE;
    }
  }
  std::vector<double> tmp_gathered_x(nrwu), tmp_gathered_y(nrwu), tmp_gathered_z(nrwu);
  std::vector<double> copy_prop_x(prop_x);
  std::vector<double> copy_prop_y(prop_y);
  std::vector<double> copy_prop_z(prop_z);
  ReductionKit redk(nrwu, rps, rwu_abstracts.data(), atom_counts.data());
  const double* y_ptr = (prop_y.size() == prop_x.size()) ? copy_prop_y.data() : nullptr;
  const double* z_ptr = (prop_z.size() == prop_x.size()) ? copy_prop_z.data() : nullptr;
  GenericRdSubstrate rsbs(prop_x.data(), y_ptr, z_ptr, tmp_gathered_x.data(),
                          tmp_gathered_y.data(), tmp_gathered_z.data(), copy_prop_x.data(),
                          copy_prop_y.data(), copy_prop_z.data());
  evalReduction(&rsbs, redk, ReductionStage::ALL_REDUCE, ReductionGoal::NORMALIZE);
  std::vector<double> norm_prop_x(prop_x), norm_prop_y(prop_y), norm_prop_z(prop_z);
  double* norm_xptr = norm_prop_x.data();
  double* norm_yptr = norm_prop_y.data();
  double* norm_zptr = norm_prop_z.data();
  for (int i = 0; i < nsys; i++) {
    const double xmag = magnitude(&norm_xptr[atom_offsets[i]], atom_counts[i]);
    const double ymag = (have_y) ? magnitude(&norm_yptr[atom_offsets[i]], atom_counts[i]) : 0.0;
    const double zmag = (have_z) ? magnitude(&norm_zptr[atom_offsets[i]], atom_counts[i]) : 0.0;
    const double norm_factor = 1.0 / sqrt((xmag * xmag) + (ymag * ymag) + (zmag * zmag));
    elementwiseMultiply(&norm_xptr[atom_offsets[i]], atom_counts[i], norm_factor);
    if (have_y) {
      elementwiseMultiply(&norm_yptr[atom_offsets[i]], atom_counts[i], norm_factor);
    }
    if (have_z) {
      elementwiseMultiply(&norm_zptr[atom_offsets[i]], atom_counts[i], norm_factor);
    }
  }
  check(norm_prop_x, RelationalOperator::EQUAL, Approx(copy_prop_x).margin(tiny),
        "Normalization by reduction work units was not accomplished correctly in the first of " +
        std::to_string(nprop) + " arrays.");
  if (have_y) {
    check(norm_prop_y, RelationalOperator::EQUAL, Approx(copy_prop_y).margin(tiny),
          "Normalization by reduction work units was not accomplished correctly in the second "
          "of " + std::to_string(nprop) + " arrays.");
  }
  if (have_z) {
    check(norm_prop_z, RelationalOperator::EQUAL, Approx(copy_prop_z).margin(tiny),
          "Normalization by reduction work units was not accomplished correctly in the third "
          "of " + std::to_string(nprop) + " arrays.");
  }

  // Reset the arrays and try reduction for centering on zero.
  std::vector<double> cent_prop_x(prop_x), cent_prop_y(prop_y), cent_prop_z(prop_z);
  for (int i = 0; i < nsys; i++) {
    const int jmin = atom_offsets[i];
    const int jmax = jmin + atom_counts[i];
    for (int j = jmin; j < jmax; j++) {
      copy_prop_x[j] = prop_x[j];
      if (have_y) {
        copy_prop_y[j] = prop_y[j];
      }
      if (have_z) {
        copy_prop_z[j] = prop_z[j];
      }
    }
  }
  evalReduction(&rsbs, redk, ReductionStage::ALL_REDUCE, ReductionGoal::CENTER_ON_ZERO);
  double* cent_xptr = cent_prop_x.data();
  double* cent_yptr = cent_prop_y.data();
  double* cent_zptr = cent_prop_z.data();
  for (int i = 0; i < nsys; i++) {
    const double xave = mean(&cent_xptr[atom_offsets[i]], atom_counts[i]);
    addScalarToVector(&cent_xptr[atom_offsets[i]], atom_counts[i], -xave);
    if (have_y) {
      const double yave = mean(&cent_yptr[atom_offsets[i]], atom_counts[i]);
      addScalarToVector(&cent_yptr[atom_offsets[i]], atom_counts[i], -yave);
    }
    if (have_z) {
      const double zave = mean(&cent_zptr[atom_offsets[i]], atom_counts[i]);
      addScalarToVector(&cent_zptr[atom_offsets[i]], atom_counts[i], -zave);
    }
  }
  check(cent_prop_x, RelationalOperator::EQUAL, Approx(copy_prop_x).margin(tiny),
        "Zero-centering by reduction work units was not accomplished correctly in the first of " +
        std::to_string(nprop) + " arrays.");
  if (have_y) {
    check(cent_prop_y, RelationalOperator::EQUAL, Approx(copy_prop_y).margin(tiny),
          "Zero-centering by reduction work units was not accomplished correctly in the second "
          "of " + std::to_string(nprop) + " arrays.");
  }
  if (have_z) {
    check(cent_prop_z, RelationalOperator::EQUAL, Approx(copy_prop_z).margin(tiny),
          "Zero-centering by reduction work units was not accomplished correctly in the third "
          "of " + std::to_string(nprop) + " arrays.");
  }
}

//-------------------------------------------------------------------------------------------------
// Verify that a given set of volume partitions exactly covers all volume of each periodic system.
//-------------------------------------------------------------------------------------------------
void verifyBrickCoverage(const Brickwork &bw) {
  const int nsys = bw.getSystemCount();
  const int nbricks = bw.getBrickCount();
  std::vector<std::vector<int>> coverage(nsys);
  for (int i = 0; i < nsys; i++) {
    const int3 idims = bw.getSystemDimensions(i);
    coverage[i] = std::vector<int>(idims.x * idims.y * idims.z, 0);
  }
  for (int bci = 0; bci < nbricks; bci++) {
    const int sysid = bw.getSystemMembership(bci);
    const int3 br_orig = bw.getBrickOrigin(bci);
    const int3 br_len  = bw.getBrickLengths(bci);
    for (int k = 0; k < br_len.z; k++) {
      for (int j = 0; j < br_len.y; j++) {
        const int jk_idx = (((br_orig.z + k) * br_len.y) + br_orig.y + j) * br_len.x;
        for (int i = 0; i < br_len.x; i++) {
          coverage[sysid][jk_idx + br_orig.x + i] += 1;
        }
      }
    }
  }
  std::vector<double> mean_coverage(nsys);
  std::vector<double> max_coverage(nsys);
  std::vector<double> min_coverage(nsys);
  for (int i = 0; i < nsys; i++) {
    mean_coverage[i] = mean(coverage[i]);
  }
  check(mean_coverage, RelationalOperator::EQUAL, std::vector<double>(nsys, 1.0), "The mean "
        "coverage of volume partitions is incorrect.");
}

//-------------------------------------------------------------------------------------------------
// Test the way valence work units are sized.
//-------------------------------------------------------------------------------------------------
void testVwuSizingAlgorithm() {

  // Try mock ligands with small hypothetical GPUs
  Xoroshiro128pGenerator xrs(88014239);
  std::vector<int> vwu_cap(18);
  std::vector<ValenceKernelSize> kwidths(18);
  std::vector<int> atom_counts(32);
  uniformRand(&xrs, &atom_counts, 48.0, 1.0);
  for (int i = 0; i < 32; i++) {
    atom_counts[i] += 16;
  }
  int iplc = 0;
  for (int i = 32; i > 0; i /= 2) {
    vwu_cap[iplc] = calculateValenceWorkUnitSize(atom_counts, i, &kwidths[iplc]);
    iplc++;
  }
  atom_counts.resize(768);
  uniformRand(&xrs, &atom_counts, 48.0, 1.0);
  for (int i = 0; i < 768; i++) {
    atom_counts[i] += 16;
  }
  for (int i = 32; i > 0; i /= 2) {
    vwu_cap[iplc] = calculateValenceWorkUnitSize(atom_counts, i, &kwidths[iplc]);
    iplc++;
  }
  atom_counts.resize(9216);
  uniformRand(&xrs, &atom_counts, 48.0, 1.0);
  for (int i = 0; i < 9216; i++) {
    atom_counts[i] += 16;
  }
  for (int i = 32; i > 0; i /= 2) {
    vwu_cap[iplc] = calculateValenceWorkUnitSize(atom_counts, i, &kwidths[iplc]);
    iplc++;
  }
  const std::vector<int> vwu_cap_ans = { 51, 68, 68, 68, 68, 68,
                                         51, 68, 68, 68, 68, 68,
                                         68, 68, 68, 68, 68, 68 };
  check(vwu_cap, RelationalOperator::EQUAL, vwu_cap_ans, "The calcualted valence work unit sizes "
        "did not meet expectations.");
  std::vector<int> kwidths_int(iplc);
  for (int i = 0; i < iplc; i++) {
    kwidths_int[i] = static_cast<int>(kwidths[i]);
  }
  std::vector<int> kwidths_ans(kwidths.size(), static_cast<int>(ValenceKernelSize::SM));
  kwidths_ans[0] = static_cast<int>(ValenceKernelSize::XL);
  kwidths_ans[1] = static_cast<int>(ValenceKernelSize::XL);
  kwidths_ans[2] = static_cast<int>(ValenceKernelSize::LG);
  kwidths_ans[3] = static_cast<int>(ValenceKernelSize::MD);
  check(kwidths_int, RelationalOperator::EQUAL, kwidths_ans, "The selected valence thread block "
        "sizes did not meet expectations.");

  // Try protein-sized systems with an approximation of an Ampere-era NVIDIA GPU
  const std::vector<int> prot_atom_counts = { 2000, 4000, 6000, 8000, 10000, 20000, 30000, 40000 };
  const int nprot_trials = prot_atom_counts.size();
  std::vector<int> prot_vwu_cap(nprot_trials), prot_kwidths_int(nprot_trials);
  std::vector<ValenceKernelSize> prot_kwidths(nprot_trials);
  for (int i = 0; i < nprot_trials; i++) {
    const std::vector<int> tmp_atom_counts(1, prot_atom_counts[i]);
    prot_vwu_cap[i] = calculateValenceWorkUnitSize(tmp_atom_counts, 68, &prot_kwidths[i]);
  }
  const std::vector<int> prot_vwu_cap_ans = { 34, 51, 68, 85, 102, 170, 238, 323 };
  check(prot_vwu_cap, RelationalOperator::EQUAL, prot_vwu_cap_ans, "The valence work unit sizes "
        "selected for a series of protein-sized systems do not meet expectations.");
  for (int i = 0; i < nprot_trials; i++) {
    prot_kwidths_int[i] = static_cast<int>(prot_kwidths[i]);
  }
  const std::vector<int> prot_kwidths_ans(nprot_trials, static_cast<int>(ValenceKernelSize::XL));
  check(prot_kwidths_int, RelationalOperator::EQUAL, prot_kwidths_ans, "The valence thread block "
        "sizes selected for a series of protein-sized systems do not meet expectations.");
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

  // Section 1
  section("SystemCache construction");

  // Section 2
  section("PhaseSpaceSynthesis layout");

  // Section 3
  section("Valence work unit construction");

  // Section 4
  section("Static exclusion mask compilation");

  // Section 5
  section("Non-bonded work unit construction");

  // Section 6
  section("Reduction work unit construction");

  // Section 7
  section("Periodic unit cell subdivision");

  // Test the construction of a systems cache (one possible precursor to a synthesis of
  // coordinates and / or topologies)
  section(1);
  const char osc = osSeparator();
  const TestPriority test_sysc = (oe.getTemporaryDirectoryAccess()) ? TestPriority::CRITICAL :
                                                                      TestPriority::ABORT;
  const std::string mdin_name = oe.getTemporaryDirectoryPath() + osc + "cachetest.in";
  const std::string testdir =  oe.getStormmSourcePath() + osc + "test" + osc;
  if (oe.getTemporaryDirectoryAccess()) {
    std::ofstream foutp = openOutputFile(mdin_name, PrintSituation::OVERWRITE, "Prepare a brief "
                                         "input file for testing the SystemCache machinery.");
    std::string buffer("&files\n  -p ");
    buffer += testdir + "Namelists" + osc + "topol" + osc + ".*.top\n  -c ";
    buffer += testdir + "Namelists" + osc + "coord" + osc + ".*.inpcrd\n&end\n";
    foutp.write(buffer.c_str(), buffer.size());
    foutp.close();
    oe.logFileCreated(mdin_name);
  }
  else {
    rtWarn("The temporary directory " + oe.getTemporaryDirectoryPath() + " is not writeable.  "
           "This will lead to tests either failing or being impossible to run.", "test_synthesis");
  }
  const bool mdin_exists = (getDrivePathType(mdin_name) == DrivePathType::FILE);
  const TextFile tf = (mdin_exists) ? TextFile(mdin_name) : TextFile();
  int start_line = 0;
  FilesControls fcon(tf, &start_line);
  SystemCache sysc(fcon, ExceptionResponse::SILENT);
  const int nsys = sysc.getSystemCount();
  check(nsys, RelationalOperator::EQUAL, 16, "The number of systems detected with regular "
        "expression searching did not meet expectations.  This may indicate a problem with the "
        "${STORMM_SOURCE} environment variable that did not show up in test program setup.",
        test_sysc);
  const std::vector<AtomGraph*> ag_ptr = sysc.getSystemTopologyPointer();
  const std::vector<PhaseSpace*> ps_ptr = sysc.getCoordinatePointer();
  std::vector<int> ag_atoms(nsys);
  std::vector<int> ps_atoms(nsys);
  for (int i = 0; i < nsys; i++) {
    ag_atoms[i] = ag_ptr[i]->getAtomCount();
    ps_atoms[i] = ps_ptr[i]->getAtomCount();
  }
  check(ag_atoms, RelationalOperator::EQUAL, ps_atoms, "The numbers of atoms found with pointers "
        "matched to topologies and coordinate sets of the same SystemCache object disagree.",
        test_sysc);
  std::vector<int> ag_ref_atoms(nsys);
  std::vector<int> ps_ref_atoms(nsys);
  for (int i = 0; i < nsys; i++) {
    const AtomGraph  &ag_ref = sysc.getSystemTopology(i);
    const PhaseSpace &ps_ref = sysc.getCoordinates(i);
    ag_ref_atoms[i] = ag_ref.getAtomCount();
    ps_ref_atoms[i] = ps_ref.getAtomCount();
  }
  check(ag_ref_atoms, RelationalOperator::EQUAL, ps_ref_atoms, "The numbers of atoms found with "
        "references to matching topologies and coordinate sets of the same SystemCache object "
        "disagree.", test_sysc);
  check(ag_ref_atoms, RelationalOperator::EQUAL, ps_atoms, "The numbers of atoms found with "
        "pointers and references matched to topologies and coordinate sets of the same "
        "SystemCache object disagree.", test_sysc);
  const std::vector<std::string> base_aa = { "arg", "trp", "gly", "gly_arg", "gly_trp", "gly_phe",
                                             "tyr", "gly_pro", "gly_tyr", "phe", "gly", "gly_gly",
                                             "pro", "lys", "gly_ala", "gly_lys", "ala" };
  const int nbase_aa = base_aa.size();
  std::vector<std::string> rst_outputs, trj_outputs;
  rst_outputs.reserve(nbase_aa);
  trj_outputs.reserve(nbase_aa);
  for (int i = 0; i < nbase_aa; i++) {
    rst_outputs.push_back("md_" + base_aa[i] + ".rst");
    trj_outputs.push_back("md_" + base_aa[i] + ".crd");
  }
  bool rst_outside_list = false;
  bool trj_outside_list = false;
  for (int i = 0; i < nsys; i++) {
    rst_outside_list = (rst_outside_list ||
                        (findStringInVector(rst_outputs,
                                            sysc.getCheckpointName(i)) == nbase_aa));
    trj_outside_list = (trj_outside_list ||
                        (findStringInVector(trj_outputs,
                                            sysc.getTrajectoryName(i)) == nbase_aa));
  }
  check(rst_outside_list == false, "One of the SystemCache's checkpoint file names was outside "
        "the expected list.  The name extension feature, for differentiating names of multiple "
        "systems under the same label, may be malfunctioning.", test_sysc);
  check(trj_outside_list == false, "One of the SystemCache's trajectory file names was outside "
        "the expected list.  The name extension feature, for differentiating names of multiple "
        "systems under the same label, may be malfunctioning.", test_sysc);

  // Make a new systems cache and verify that multiple frames under the same label are properly
  // handled.
  std::string multi_label_deck("&files\n");
  multi_label_deck += "  -sys { -c " + testdir + "MoleculeFormat" + osc + "sulfonamide_rots.sdf "
                      "-p " + testdir + "Topology" + osc + "sulfonamide.top -label sulfonamide "
                      "frame_end -1 }";
  multi_label_deck += "  -sys { -c " + testdir + "Trajectory" + osc + "stereo_L1_vs.inpcrd -p " +
                      testdir + "Topology" + osc + "stereo_L1_vs.top -label stereo -n 16 }"
                      "\n&end\n";
  const TextFile multi_label_tf(multi_label_deck, TextOrigin::RAM);
  start_line = 0;
  const FilesControls multi_label_fcon(multi_label_tf, &start_line);
  SystemCache multi_label_sysc(multi_label_fcon, ExceptionResponse::SILENT);
  check(multi_label_sysc.getCheckpointName(10), RelationalOperator::EQUAL, "md_0_10.rst",
        "The SystemCache object does not produced the expected non-colliding name for a frame of "
        "an SD file.", test_sysc);
  check(multi_label_sysc.getCheckpointName(38), RelationalOperator::EQUAL, "md_1_2.rst",
        "The SystemCache object does not produced the expected non-colliding name for a replica "
        "of an coordinate system first read from an Amber ASCII restart file.", test_sysc);

  std::string cat_output_deck("&files\n");
  cat_output_deck += "  -sys { -c " + testdir + "MoleculeFormat" + osc + "sulfonamide_rots.sdf "
                     "-p " + testdir + "Topology" + osc + "sulfonamide.top -label sulfonamide "
                     "frame_end -1 -r min_sulf.sdf r_kind SDF }";
  cat_output_deck += "  -sys { -c " + testdir + "Trajectory" + osc + "stereo_L1_vs.inpcrd -p " +
                     testdir + "Topology" + osc + "stereo_L1_vs.top -label stereo -n 16 "
                     "-x min_stereo.crd r_kind AMBER_CRD }\n&end\n";
  const TextFile cat_output_tf(cat_output_deck, TextOrigin::RAM);
  start_line = 0;
  const FilesControls cat_output_fcon(cat_output_tf, &start_line);
  SystemCache cat_output_sysc(cat_output_fcon, ExceptionResponse::SILENT);
  std::vector<int> prnt_situ_ans(52, static_cast<int>(PrintSituation::APPEND));
  std::vector<int> prnt_situ(52);
  prnt_situ_ans[0] = static_cast<int>(PrintSituation::OPEN_NEW);
  prnt_situ_ans[36] = prnt_situ_ans[0];
  const CoordinateFileRole chk_purpose = CoordinateFileRole::CHECKPOINT;
  for (int i = 0; i < cat_output_sysc.getSystemCount(); i++) {
    prnt_situ[i] = static_cast<int>(cat_output_sysc.getPrintingProtocol(chk_purpose, i));
  }
  check(prnt_situ, RelationalOperator::EQUAL, prnt_situ_ans, "The printing protocols for multiple "
        "systems grouped under the same label feeding into a common output file able to accept "
        "multiple frames do not meet expectations.");

  // Create some topologies and coordinate sets.
  Xoroshiro128pGenerator my_prng(oe.getRandomSeed());
  const std::string base_crd_name = testdir + "Trajectory";
  const std::string base_top_name = testdir + "Topology";
  const std::string tip3p_crd_name = base_crd_name + osc + "tip3p.inpcrd";
  const std::string tip3p_top_name = base_top_name + osc + "tip3p.top";
  const std::string tip4p_crd_name = base_crd_name + osc + "tip4p.inpcrd";
  const std::string tip4p_top_name = base_top_name + osc + "tip4p.top";
  const std::string trpcage_crd_name = base_crd_name + osc + "trpcage_in_water.inpcrd";
  const std::string trpcage_top_name = base_top_name + osc + "trpcage_in_water.top";
  const bool files_exist = (getDrivePathType(tip3p_crd_name) == DrivePathType::FILE &&
                            getDrivePathType(tip3p_top_name) == DrivePathType::FILE &&
                            getDrivePathType(tip4p_crd_name) == DrivePathType::FILE &&
                            getDrivePathType(tip4p_top_name) == DrivePathType::FILE &&
                            getDrivePathType(trpcage_crd_name) == DrivePathType::FILE &&
                            getDrivePathType(trpcage_top_name) == DrivePathType::FILE);
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  PhaseSpace tip3p_ps, tip4p_ps, trpcage_ps, tip3p_ps_2, tip4p_ps_2, trpcage_ps_2;
  AtomGraph tip3p_ag, tip4p_ag, trpcage_ag;
  if (files_exist) {
    tip3p_ag.buildFromPrmtop(tip3p_top_name);
    tip3p_ps.buildFromFile(tip3p_crd_name, CoordinateFileKind::AMBER_INPCRD);
    tip4p_ag.buildFromPrmtop(tip4p_top_name);
    tip4p_ps.buildFromFile(tip4p_crd_name, CoordinateFileKind::AMBER_INPCRD);
    trpcage_ag.buildFromPrmtop(trpcage_top_name);
    trpcage_ps.buildFromFile(trpcage_crd_name, CoordinateFileKind::AMBER_INPCRD);
  }
  else {
    rtWarn("The topology and coordinate files for the TIP3P and TIP4P water boxes as well as the "
           "Trp-cage miniprotein in water must be available in ${STORMM_SOURCE}/test/ "
           "subdirectories Topology and Trajectory, respectively.  Check the $STORMM_SOURCE "
           "environment variable to make sure that it is set properly.  A number of tests will "
           "be skipped.", "test_phase_space_synthesis");
  }
  
  // Make a PhaseSpaceSynthesis the meticulous way
  section(2);
  std::vector<PhaseSpace>  psv = { tip3p_ps, tip4p_ps, trpcage_ps, tip3p_ps, tip4p_ps,
                                   tip3p_ps, tip3p_ps, tip4p_ps, trpcage_ps };
  const std::vector<AtomGraph*> agv = { &tip3p_ag, &tip4p_ag, &trpcage_ag, &tip3p_ag, &tip4p_ag,
                                        &tip3p_ag, &tip3p_ag, &tip4p_ag, &trpcage_ag };
  for (size_t i = 3; i < psv.size(); i++) {
    PhaseSpaceWriter pswi = psv[i].data();
    for (int j = 0; j < pswi.natom; j++) {
      pswi.xcrd[j] += 0.02 * my_prng.gaussianRandomNumber();
      pswi.ycrd[j] += 0.02 * my_prng.gaussianRandomNumber();
      pswi.zcrd[j] += 0.02 * my_prng.gaussianRandomNumber();
    }
  }
  PhaseSpaceSynthesis psynth(psv, agv);

  // Try extracting a system from it
  PhaseSpace tip3p_ps_copy(tip3p_ps.getAtomCount(), tip3p_ps.getUnitCellType());
  psynth.extractSystem(&tip3p_ps_copy, 3);
  PhaseSpaceWriter tip3p_orig_writer = tip3p_ps.data();
  PhaseSpaceWriter tip3p_muta_writer = psv[3].data();
  PhaseSpaceWriter tip3p_copy_writer = tip3p_ps_copy.data();
  std::vector<double> y_orig(tip3p_orig_writer.natom);
  std::vector<double> y_muta(tip3p_orig_writer.natom);
  std::vector<double> y_copy(tip3p_orig_writer.natom);
  for (int i = 0; i < tip3p_orig_writer.natom; i++) {
    y_orig[i] = tip3p_orig_writer.ycrd[i];
    y_muta[i] = tip3p_muta_writer.ycrd[i];
    y_copy[i] = tip3p_copy_writer.ycrd[i];
  }
  check(y_muta, RelationalOperator::EQUAL, y_copy, "The PhaseSpaceSynthesis object was not able "
        "to return the correct image of one of its systems.", do_tests);
  const std::vector<double> scaling_answer = {
    default_globalpos_scale_lf, default_localpos_scale_lf, default_velocity_scale_lf,
    default_force_scale_lf, default_inverse_globalpos_scale_lf, default_inverse_localpos_scale_lf,
    default_inverse_velocity_scale_lf, default_inverse_force_scale_lf };
  const std::vector<int> scale_bits_answer = {
    default_globalpos_scale_bits, default_localpos_scale_bits, default_velocity_scale_bits,
    default_force_scale_bits };
  PsSynthesisWriter psynth_w = psynth.data();
  const std::vector<double> scaling_result = { psynth_w.gpos_scale, psynth_w.lpos_scale,
                                               psynth_w.vel_scale, psynth_w.frc_scale,
                                               psynth_w.inv_gpos_scale, psynth_w.inv_lpos_scale,
                                               psynth_w.inv_vel_scale, psynth_w.inv_frc_scale };
  const std::vector<float> scaling_result_f = { psynth_w.gpos_scale_f, psynth_w.lpos_scale_f,
                                                psynth_w.vel_scale_f, psynth_w.frc_scale_f,
                                                psynth_w.inv_gpos_scale_f,
                                                psynth_w.inv_lpos_scale_f,
                                                psynth_w.inv_vel_scale_f,
                                                psynth_w.inv_frc_scale_f };
  const std::vector<int> scale_bits_result = { psynth_w.gpos_bits, psynth_w.lpos_bits,
                                               psynth_w.vel_bits, psynth_w.frc_bits };
  check(scaling_result, RelationalOperator::EQUAL,
        Approx(scaling_answer, ComparisonType::RELATIVE, stormm::constants::verytiny),
        "Double-precision scaling constants found in the PhaseSpaceSynthesis object's writer do "
        "not meet expectations.", do_tests);
  check(scaling_result_f, RelationalOperator::EQUAL,
        Approx(scaling_answer, ComparisonType::RELATIVE, stormm::constants::verytiny),
        "Single-precision scaling constants found in the PhaseSpaceSynthesis object's writer do "
        "not meet expectations.", do_tests);
  check(scale_bits_result, RelationalOperator::EQUAL, scale_bits_answer, "Fixed-precision bit "
        "counts found in the PhaseSpaceSynthesis object's writer do not meet expectations.",
        do_tests);

  // Check that the copy assignment works
  PhaseSpaceSynthesis psynth_copy = psynth;
  std::vector<double> orig_xyz, copy_xyz;
  for (int i = 0; i < psynth.getSystemCount(); i++) {
    const CoordinateFrame icf = psynth.exportCoordinates(i);
    const std::vector<double> icf_xyz = icf.getInterlacedCoordinates();
    orig_xyz.insert(orig_xyz.end(), icf_xyz.begin(), icf_xyz.end());
    const CoordinateFrame ccf = psynth_copy.exportCoordinates(i);
    const std::vector<double> ccf_xyz = ccf.getInterlacedCoordinates();
    copy_xyz.insert(copy_xyz.end(), ccf_xyz.begin(), ccf_xyz.end());
  }
  check(orig_xyz, RelationalOperator::EQUAL, copy_xyz, "Atom positions were not properly "
        "replicated by a deep copy of the PhaseSpaceSynthesis.", do_tests);
  PsSynthesisWriter copy_psw = psynth_copy.data();
  Xoroshiro128pGenerator xrs_dc;
  const size_t np_atom_count = orig_xyz.size();
  for (int i = 0; i < copy_psw.system_count; i++) {
    addRandomNoise(&xrs_dc, copy_psw.xcrd, copy_psw.xcrd_ovrf, np_atom_count, 1.0,
                   copy_psw.gpos_scale);
    addRandomNoise(&xrs_dc, copy_psw.ycrd, copy_psw.ycrd_ovrf, np_atom_count, 1.0,
                   copy_psw.gpos_scale);
    addRandomNoise(&xrs_dc, copy_psw.zcrd, copy_psw.zcrd_ovrf, np_atom_count, 1.0,
                   copy_psw.gpos_scale);
  }
  std::vector<double> orig_xyz2, copy_xyz2;
  for (int i = 0; i < psynth.getSystemCount(); i++) {
    const CoordinateFrame icf = psynth.exportCoordinates(i);
    const std::vector<double> icf_xyz = icf.getInterlacedCoordinates();
    orig_xyz2.insert(orig_xyz2.end(), icf_xyz.begin(), icf_xyz.end());
    const CoordinateFrame ccf = psynth_copy.exportCoordinates(i);
    const std::vector<double> ccf_xyz = ccf.getInterlacedCoordinates();
    copy_xyz2.insert(copy_xyz2.end(), ccf_xyz.begin(), ccf_xyz.end());
  }
  check(orig_xyz, RelationalOperator::EQUAL, orig_xyz2, "Atom positions of the original "
        "PhaseSpaceSynthesis were corrupted by modifying the deep copy.", do_tests);
  const double noise_introduced = meanUnsignedError(orig_xyz, copy_xyz2);
  check(noise_introduced, RelationalOperator::EQUAL, Approx(9.10463).margin(3.0e-5),
        "Modifying the deep copy was ineffective.");

  // Create a SystemCache containing all of the systems of the first synthesis, plus meta data.
  std::string mock_sysc_deck("&files\n");
  mock_sysc_deck += "  -sys { -c " + tip3p_crd_name + " -p " + tip3p_top_name + 
                    " -label water -x tip3p_output.crd -r tip3p_output.rst }";
  mock_sysc_deck += "  -sys { -c " + tip4p_crd_name + " -p " + tip4p_top_name +
                    " -label water -x tip4p_output.crd -r tip4p_output.rst }";
  mock_sysc_deck += "  -sys { -c " + trpcage_crd_name + " -p " + trpcage_top_name +
                    " -label protein  -x trpcage_output.crd -r trpcage_output.rst }";
  mock_sysc_deck += "  -sys { -c " + tip3p_crd_name + " -p " + tip3p_top_name + 
                    " -label steam -x steam_output.crd -r steam_output.rst }";
  mock_sysc_deck += "  -sys { -c " + tip4p_crd_name + " -p " + tip4p_top_name +
                    " -label steam -x steam_output.crd -r steam_output.rst }";
  mock_sysc_deck += "&end\n";
  const TextFile mock_sysc_tf(mock_sysc_deck, TextOrigin::RAM);
  start_line = 0;
  const FilesControls mock_sysc_fcon(mock_sysc_tf, &start_line);
  SystemCache mock_sysc(mock_sysc_fcon, ExceptionResponse::SILENT);
  const std::vector<int> synth_index = { 0, 1, 2, 3, 4, 0, 3, 1, 2 };
  SynthesisCacheMap mock_map(synth_index, mock_sysc, psynth);
  CHECK_THROWS_SOFT(SynthesisCacheMap scmap(tileVector(synth_index, 2), mock_sysc, psynth),
                    "A SynthesisCacheMap was constructed with an insufficient number of "
                    "cache-to-synthesis indices.", do_tests);
  const std::vector<int> tip3p_derivations   = mock_map.getTopologyGroup(&tip3p_ag);
  const std::vector<int> tip4p_derivations   = mock_map.getTopologyGroup(&tip4p_ag);
  const std::vector<int> steam_derivations   = mock_map.getLabelGroup("steam");
  const std::vector<int> protein_derivations = mock_map.getLabelGroup("protein");
  const std::vector<int> trpcage_derivations = mock_map.getTopologyGroup(&trpcage_ag);
  const std::vector<int> tip3p_derivations_ans = { 0, 3, 5, 6 };
  const std::vector<int> tip4p_derivations_ans = { 1, 4, 7 };
  const std::vector<int> steam_derivations_ans = { 3, 4, 6 };
  check(tip3p_derivations, RelationalOperator::EQUAL, tip3p_derivations_ans, "Systems derived "
        "from the TIP3P topology were not correctly identified by the map.", do_tests);
  check(tip4p_derivations, RelationalOperator::EQUAL, tip4p_derivations_ans, "Systems derived "
        "from the TIP3P topology were not correctly identified by the map.", do_tests);
  check(steam_derivations, RelationalOperator::EQUAL, steam_derivations_ans, "Systems associated "
        "with the 'steam' label were not correctly identified by the map.", do_tests);
  check(protein_derivations, RelationalOperator::EQUAL, trpcage_derivations, "Systems associated "
        "with the 'protein' label should match those associated with the Trp-cage topology in the "
        "map.", do_tests);
  
  // Make a second coordinate synthesis, this time with different fixed-precision settings
  PhaseSpaceSynthesis psynth2(psv, agv, 24, 25, 40, 28);
  PsSynthesisWriter psynth_w2 = psynth2.data();
  const std::vector<double> scaling_answer2 = { pow(2.0, 24), pow(2.0, 25), pow(2.0, 40),
                                                pow(2.0, 28), 1.0 / pow(2.0, 24),
                                                1.0 / pow(2.0, 25), 1.0 / pow(2.0, 40),
                                                1.0 / pow(2.0, 28) };
  const std::vector<float> scaling_answer2_f(scaling_answer2.begin(), scaling_answer2.end());
  const std::vector<int> scale_bits_answer2 = { 24, 25, 40, 28 };
  const std::vector<double> scaling_result2 = { psynth_w2.gpos_scale, psynth_w2.lpos_scale,
                                                psynth_w2.vel_scale, psynth_w2.frc_scale,
                                                psynth_w2.inv_gpos_scale, psynth_w2.inv_lpos_scale,
                                                psynth_w2.inv_vel_scale, psynth_w2.inv_frc_scale };
  const std::vector<float> scaling_result2_f = { psynth_w2.gpos_scale_f, psynth_w2.lpos_scale_f,
                                                 psynth_w2.vel_scale_f, psynth_w2.frc_scale_f,
                                                 psynth_w2.inv_gpos_scale_f,
                                                 psynth_w2.inv_lpos_scale_f,
                                                 psynth_w2.inv_vel_scale_f,
                                                 psynth_w2.inv_frc_scale_f };
  const std::vector<int> scale_bits_result2 = { psynth_w2.gpos_bits, psynth_w2.lpos_bits,
                                                psynth_w2.vel_bits, psynth_w2.frc_bits };
  check(scaling_result2, RelationalOperator::EQUAL,
        Approx(scaling_answer2, ComparisonType::RELATIVE, stormm::constants::verytiny),
        "Double-precision scaling constants found in the PhaseSpaceSynthesis object's writer do "
        "not meet expectations.", do_tests);
  check(scaling_result2_f, RelationalOperator::EQUAL,
        Approx(scaling_answer2, ComparisonType::RELATIVE, stormm::constants::verytiny),
        "Single-precision scaling constants found in the PhaseSpaceSynthesis object's writer do "
        "not meet expectations.", do_tests);
  check(scale_bits_result2, RelationalOperator::EQUAL, scale_bits_answer2, "Fixed-precision bit "
        "counts found in the PhaseSpaceSynthesis object's writer do not meet expectations.");
  psynth2.extractSystem(&tip3p_ps_copy, 3);
  for (int i = 0; i < tip3p_orig_writer.natom; i++) {
    y_orig[i] = tip3p_orig_writer.ycrd[i];
    y_muta[i] = tip3p_muta_writer.ycrd[i];
    y_copy[i] = tip3p_copy_writer.ycrd[i];
  }
  check(y_muta, RelationalOperator::EQUAL, Approx(y_copy).margin(1.0e-6),
        "The PhaseSpaceSynthesis object returns an incorrect image of one of its systems, even "
        "after compensating for a lower global position resolution", do_tests);
  check(y_muta, RelationalOperator::NOT_EQUAL, Approx(y_copy).margin(1.0e-8),
        "The PhaseSpaceSynthesis object returns an image of one of its systems with higher "
        "fidelity to the original than expected.", do_tests);

  // Try extracting the unique topologies and unique topology indices from a PhaseSpaceSynthesis
  const std::vector<AtomGraph*> psv_tops = psynth.getUniqueTopologies();
  const std::vector<int> psv_top_idx = psynth.getUniqueTopologyExampleIndices();
  const std::vector<int> psv_top_idx_ans = { 0, 1, 2 };
  const bool top_ordered = (psv_tops[0] == &tip3p_ag && psv_tops[1] == &tip4p_ag &&
                            psv_tops[2] == &trpcage_ag);
  check(top_ordered, "The order of unique topologies produced by a PhaseSpaceSynthesis object is "
        "not as expected.", do_tests);
  check(psv_top_idx, RelationalOperator::EQUAL, psv_top_idx_ans, "The sequence of unique topology "
        "indices produced by a PhaseSpaceSynthesis object is not as expected.", do_tests);

  // Check the high-precision modes of a phase-space synthesis: are the extended coordinate,
  // velocity, and force arrays reporting the correct results?
  PhaseSpaceSynthesis highres(psv, agv, 54, 24, 60, 55);
  PsSynthesisWriter highres_w = highres.data();
  std::vector<double> pos_mues(highres_w.system_count);
  std::vector<double> vel_mues(highres_w.system_count);
  std::vector<double> frc_mues(highres_w.system_count);
  for (int i = 0; i < highres_w.system_count; i++) {
    const int jlim = highres_w.atom_starts[i] + highres_w.atom_counts[i];
    PhaseSpace sysi_ps = highres.exportSystem(i);
    PhaseSpaceWriter sysi_psw = sysi_ps.data();
    for (int j = highres_w.atom_starts[i]; j < jlim; j++) {
      const size_t jlocal = j - highres_w.atom_starts[i];
      const double tx_frc = 1024.0 * (0.5 - my_prng.uniformRandomNumber());
      const double ty_frc = 1024.0 * (0.5 - my_prng.uniformRandomNumber());
      const double tz_frc = 1024.0 * (0.5 - my_prng.uniformRandomNumber());
      sysi_psw.xfrc[jlocal] = tx_frc;
      sysi_psw.yfrc[jlocal] = ty_frc;
      sysi_psw.zfrc[jlocal] = tz_frc;
      hostDoubleToInt95(tx_frc * highres_w.frc_scale, &highres_w.xfrc[j], &highres_w.xfrc_ovrf[j]);
      hostDoubleToInt95(ty_frc * highres_w.frc_scale, &highres_w.yfrc[j], &highres_w.yfrc_ovrf[j]);
      hostDoubleToInt95(tz_frc * highres_w.frc_scale, &highres_w.zfrc[j], &highres_w.zfrc_ovrf[j]);
      const double tx_vel = 64.0 * (0.5 - my_prng.uniformRandomNumber());
      const double ty_vel = 64.0 * (0.5 - my_prng.uniformRandomNumber());
      const double tz_vel = 64.0 * (0.5 - my_prng.uniformRandomNumber());
      sysi_psw.xvel[jlocal] = tx_vel;
      sysi_psw.yvel[jlocal] = ty_vel;
      sysi_psw.zvel[jlocal] = tz_vel;
      hostDoubleToInt95(tx_vel * highres_w.vel_scale, &highres_w.xvel[j], &highres_w.xvel_ovrf[j]);
      hostDoubleToInt95(ty_vel * highres_w.vel_scale, &highres_w.yvel[j], &highres_w.yvel_ovrf[j]);
      hostDoubleToInt95(tz_vel * highres_w.vel_scale, &highres_w.zvel[j], &highres_w.zvel_ovrf[j]);
    }
    const PhaseSpace sysi_rb = highres.exportSystem(i);
    const TrajectoryKind tjpos = TrajectoryKind::POSITIONS;
    const TrajectoryKind tjvel = TrajectoryKind::VELOCITIES;
    const TrajectoryKind tjfrc = TrajectoryKind::FORCES;
    const std::vector<double> ps_coords = sysi_ps.getInterlacedCoordinates(tjpos);
    const std::vector<double> ps_velocs = sysi_ps.getInterlacedCoordinates(tjvel);
    const std::vector<double> ps_forces = sysi_ps.getInterlacedCoordinates(tjfrc);
    const std::vector<double> rb_coords = sysi_rb.getInterlacedCoordinates(tjpos);
    const std::vector<double> rb_velocs = sysi_rb.getInterlacedCoordinates(tjvel);
    const std::vector<double> rb_forces = sysi_rb.getInterlacedCoordinates(tjfrc);
    pos_mues[i] = meanUnsignedError(ps_coords, rb_coords);
    vel_mues[i] = meanUnsignedError(ps_velocs, rb_velocs);
    frc_mues[i] = meanUnsignedError(ps_forces, rb_forces);
  }
  const Approx dead_on = Approx(std::vector<double>(highres_w.system_count, 0.0)).margin(verytiny);
  check(pos_mues, RelationalOperator::EQUAL, dead_on, "Positions recorded in a high-precision "
        "PhaseSpaceSynthesis object do not produce the correct values when a PhaseSpace object is "
        "exported.", do_tests);
  check(vel_mues, RelationalOperator::EQUAL, dead_on, "Velocities recorded in a high-precision "
        "PhaseSpaceSynthesis object do not produce the correct values when a PhaseSpace object is "
        "exported.", do_tests);
  check(frc_mues, RelationalOperator::EQUAL, dead_on, "Forces recorded in a high-precision "
        "PhaseSpaceSynthesis object do not produce the correct values when a PhaseSpace object is "
        "exported.", do_tests);

  // Prepare valence work units for the array of topologies
  section(3);
  std::vector<RestraintApparatus> ra_vec;
  for (int i = 0; i < sysc.getTopologyCount(); i++) {
    const int example_system_idx = sysc.getSystemExampleIndex(i);
    const AtomGraph *ag_i = sysc.getTopologyPointer(i);
    const PhaseSpace &ps_i = sysc.getCoordinates(example_system_idx);
    const CoordinateFrameReader cfr_i(ps_i);
    const ChemicalFeatures chemfe_i(ag_i, ps_i, MapRotatableGroups::YES);
    const AtomMask bkbn_i(":* & @CA,N,C,O", ag_i, chemfe_i, cfr_i);
    ra_vec.emplace_back(applyHydrogenBondPreventors(ag_i, chemfe_i, 64.0, 3.1));
    ra_vec[i].addRestraints(applyPositionalRestraints(ag_i, cfr_i, bkbn_i, 16.0));
  }
  bool force_partner_counts_match = true;
  bool force_partners_match = true;
  for (int i = 0; i < sysc.getTopologyCount(); i++) {
    ValenceDelegator vdel(sysc.getTopologyPointer(i), &ra_vec[i]);
    const std::vector<ValenceWorkUnit> vwu_i = buildValenceWorkUnits(&vdel);
    const AtomGraph &ag_i = sysc.getTopology(i);
    const std::vector<std::vector<int>> fcontrib = getAtomForceContributors(ag_i, ra_vec[i]);
    const ValenceKit<double> vk = ag_i.getDoublePrecisionValenceKit();
    const VirtualSiteKit<double> vsk = ag_i.getDoublePrecisionVirtualSiteKit();
    for (int j = 0; j < vk.natom; j++) {
      const std::vector<int> vwu_deps = vdel.findForcePartners(j);
      force_partner_counts_match = (force_partner_counts_match &&
                                    vwu_deps.size() == fcontrib[j].size());
      const std::vector<ValueWithCounter<int>> missing_dependencies =
        findUnmatchedValues(vwu_deps, fcontrib[j], UniqueValueHandling::CONFIRM_ALL_COPIES);
      force_partners_match = (missing_dependencies.size() == 0LLU);
    }
  }
  check(force_partner_counts_match, "The counts of force-relevant partners in ValenceWorkUnit "
        "objects made for a series of amino acid dipeptides do not agree with numbers computed "
        "through an alternative method.", do_tests);
  check(force_partners_match, "Lists of force-relevant partners in ValenceWorkUnit objects made "
        "for a series of amino acid dipeptides do not agree with those assembled through an "
        "alternative method.", do_tests);

  // Check the valence work unit sizing for various GPU sizes and workloads.
  testVwuSizingAlgorithm();
  
  // Run diagnotics of the valence work units on a simple system, one with only a single work unit
  runValenceWorkUnitTests(sysc.getSystemTopology(0).getFileName(),
                          sysc.getCoordinates(0).getFileName(), oe, &my_prng);

  // Read a larger topology that will be forced to split its contents among several work units
  const std::string dhfr_crd_name = base_crd_name + osc + "dhfr_cmap.inpcrd";
  const std::string dhfr_top_name = base_top_name + osc + "dhfr_cmap.top";
  runValenceWorkUnitTests(dhfr_top_name, dhfr_crd_name, oe, &my_prng);

  // Read a solvated topology that will be forced to split its contents among work units, and
  // include water molecules with SETTLE constraint groups.
  const std::string ubiq_crd_name = base_crd_name + osc + "ubiquitin.inpcrd";
  const std::string ubiq_top_name = base_top_name + osc + "ubiquitin.top";
  runValenceWorkUnitTests(ubiq_top_name, ubiq_crd_name, oe, &my_prng);

  // Read additional systems with implicit solvent, then construct exclusion masks and make a
  // compilation of those masks for a synthesis of the topologies.
  const std::string trpi_crd_name = base_crd_name + osc + "trpcage.inpcrd";
  const std::string trpi_top_name = base_top_name + osc + "trpcage.top";
  const std::string lig1_crd_name = base_crd_name + osc + "stereo_L1.inpcrd";
  const std::string lig1_top_name = base_top_name + osc + "stereo_L1.top";
  const std::string lig2_crd_name = base_crd_name + osc + "stereo_L1_vs.inpcrd";
  const std::string lig2_top_name = base_top_name + osc + "stereo_L1_vs.top";
  const bool semk_exist = (getDrivePathType(trpi_crd_name) == DrivePathType::FILE &&
                           getDrivePathType(trpi_top_name) == DrivePathType::FILE &&
                           getDrivePathType(dhfr_crd_name) == DrivePathType::FILE &&
                           getDrivePathType(dhfr_top_name) == DrivePathType::FILE &&
                           getDrivePathType(lig1_crd_name) == DrivePathType::FILE &&
                           getDrivePathType(lig1_top_name) == DrivePathType::FILE &&
                           getDrivePathType(lig2_crd_name) == DrivePathType::FILE &&
                           getDrivePathType(lig2_top_name) == DrivePathType::FILE);
  AtomGraph trpi_ag, lig1_ag, lig2_ag, dhfr_ag;
  PhaseSpace trpi_ps, lig1_ps, lig2_ps, dhfr_ps;
  if (semk_exist) {
    trpi_ag.buildFromPrmtop(trpi_top_name, ExceptionResponse::SILENT);
    dhfr_ag.buildFromPrmtop(dhfr_top_name, ExceptionResponse::SILENT);
    lig1_ag.buildFromPrmtop(lig1_top_name, ExceptionResponse::SILENT);
    lig2_ag.buildFromPrmtop(lig2_top_name, ExceptionResponse::SILENT);
    trpi_ps.buildFromFile(trpi_crd_name);
    dhfr_ps.buildFromFile(dhfr_crd_name);
    lig1_ps.buildFromFile(lig1_crd_name);
    lig2_ps.buildFromFile(lig2_crd_name);
  }
  else {
    rtWarn("Files for additional systems in isolated boundary conditions were not found.  Check "
           "the ${STORMM_SOURCE} environment variable.  Subsequent tests involving these files "
           "will be skipped.", "test_synthesis");
  }
  const TestPriority do_semk_tests = (semk_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  StaticExclusionMask dhfr_se(&dhfr_ag);
  StaticExclusionMask trpi_se(&trpi_ag);
  StaticExclusionMask lig1_se(&lig1_ag);
  StaticExclusionMask lig2_se(&lig2_ag);
  const std::vector<StaticExclusionMask*> semask_list = { &dhfr_se, &trpi_se, &lig1_se, &lig2_se };
  const std::vector<int> system_indices = { 0, 1, 2, 3, 0, 1, 1, 3, 3, 2, 1 };
  const StaticExclusionMaskSynthesis poly_se(semask_list, system_indices);
  const SeMaskSynthesisReader poly_ser = poly_se.data();
  section(4);
  std::vector<int> poly_ser_natom(system_indices.size()), ag_natom(system_indices.size());
  std::vector<int> mask_mismatches(system_indices.size(), 0);
  for (size_t sysid = 0; sysid < system_indices.size(); sysid++) {
    const StaticExclusionMaskReader ser = semask_list[system_indices[sysid]]->data();
    const int natom = poly_ser.atom_counts[sysid];
    poly_ser_natom[sysid] = natom;
    ag_natom[sysid] = semask_list[system_indices[sysid]]->getTopologyPointer()->getAtomCount();
    const int nst = (poly_ser.atom_counts[sysid] + supertile_length - 1) / supertile_length;
    const int synth_st_offset = poly_ser.supertile_map_bounds[sysid];

    // Loop over all pairs (double-counting interactions) so as to test the exclusion masks from
    // both directions.
    for (int i = 0; i < natom; i++) {
      const int sti  = i / supertile_length;
      const int ti   = (i - (sti * supertile_length)) / tile_length;
      const int posi = i - (sti * supertile_length) - (ti * tile_length);
      for (int j = 0; j < natom; j++) {
        const int stj  = j / supertile_length;
        const int tj   = (j - (stj * supertile_length)) / tile_length;
        const int posj = j - (stj * supertile_length) - (tj * tile_length);

        // Get the single-system tile mask
        const int onesys_stij_map_idx = ser.supertile_map_idx[(stj * ser.supertile_stride_count) +
                                                              sti];
        const int onesys_tij_map_idx = ser.tile_map_idx[onesys_stij_map_idx +
                                                        (tj * tile_lengths_per_supertile) + ti];
        const uint onesys_mask_i = ser.mask_data[onesys_tij_map_idx + posi];

        // Get the synthesis tile mask
        const int synth_stij_map_idx = poly_ser.supertile_map_idx[(stj * nst) + sti +
                                                                  synth_st_offset];
        const int synth_tij_map_idx = poly_ser.tile_map_idx[synth_stij_map_idx + ti +
                                                            (tj * tile_lengths_per_supertile)];
        const uint synth_mask_i = poly_ser.mask_data[synth_tij_map_idx + posi];

        // Compare the two exclusion masks
        if (((onesys_mask_i >> posj) & 0x1) != ((synth_mask_i >> posj) & 0x1)) {
          mask_mismatches[sysid] += 1;
        }
      }
    }
  }
  check(poly_ser_natom, RelationalOperator::EQUAL, ag_natom, "The systems' numbers of atoms "
        "stored in the compilation of static exclusion masks do not match the numbers of atoms "
        "referenced from the original topologies.", do_semk_tests);
  check(mask_mismatches, RelationalOperator::EQUAL, std::vector<int>(system_indices.size(), 0),
        "Exclusions referenced from the compilation of systems do not match those obtained from "
        "the individual systems' static exclusion masks.", do_semk_tests);
  section(5);
  const std::vector<NonbondedWorkUnit> dhfr_nbv = buildNonbondedWorkUnits(dhfr_se);
  const std::vector<NonbondedWorkUnit> trpi_nbv = buildNonbondedWorkUnits(trpi_se);
  const std::vector<NonbondedWorkUnit> lig1_nbv = buildNonbondedWorkUnits(lig1_se);
  const std::vector<NonbondedWorkUnit> lig2_nbv = buildNonbondedWorkUnits(lig2_se);
  std::vector<int> nbt_counts(4, 0);
  for (size_t i = 0; i < dhfr_nbv.size(); i++) {
    nbt_counts[0] += dhfr_nbv[i].getTileCount();
  }
  for (size_t i = 0; i < trpi_nbv.size(); i++) {
    nbt_counts[1] += trpi_nbv[i].getTileCount();
  }
  for (size_t i = 0; i < lig1_nbv.size(); i++) {
    nbt_counts[2] += lig1_nbv[i].getTileCount();
  }
  for (size_t i = 0; i < lig2_nbv.size(); i++) {
    nbt_counts[3] += lig2_nbv[i].getTileCount();
  }
  const std::vector<int> nbt_counts_answer = { 12246, 190, 15, 21 };
  const std::vector<size_t> nbwu_counts = { dhfr_nbv.size(), trpi_nbv.size(), lig1_nbv.size(),
                                            lig2_nbv.size() };
  const std::vector<size_t> nbwu_counts_answer = { 1531LLU, 24LLU, 2LLU, 3LLU };
  check(nbwu_counts, RelationalOperator::EQUAL, nbwu_counts_answer, "The number of non-bonded "
        "work units obtained for four systems in isolated boundary conditions do not meet "
        "expectations.", do_semk_tests);
  check(nbt_counts, RelationalOperator::EQUAL, nbt_counts_answer, "The number of non-bonded tiles "
        "obtained for four systems in isolated boundary conditions do not meet expectations.",
        do_semk_tests);
  const std::vector<NonbondedWorkUnit> poly_nbv = buildNonbondedWorkUnits(poly_se);
  int poly_nbt_count = 0;
  for (size_t i = 0; i < poly_nbv.size(); i++) {
    poly_nbt_count += poly_nbv[i].getTileCount();
  }
  check(poly_nbt_count, RelationalOperator::EQUAL, 25345, "The number of non-bonded tiles "
        "obtained for the compiled set of systems does not meet expectations.", do_semk_tests);
  checkNonbondedWorkUnitCoverage(dhfr_nbv, dhfr_ag, do_semk_tests);
  checkNonbondedWorkUnitCoverage(trpi_nbv, trpi_ag, do_semk_tests);
  checkNonbondedWorkUnitCoverage(lig1_nbv, lig1_ag, do_semk_tests);
  checkNonbondedWorkUnitCoverage(lig2_nbv, lig2_ag, do_semk_tests);

  // Test reduction work unit construction using a series of dummy vectors as test data
  section(6);
  std::vector<int> small_system_counts(100), small_system_offsets(100);
  int running_total = 0;
  for (int i = 0; i < 100; i++) {
    small_system_counts[i]  = 48.0 * (1.5 - my_prng.uniformRandomNumber());
    small_system_offsets[i] = running_total;
    running_total += roundUp(small_system_counts[i], warp_size_int);
  }
  std::vector<double> small_system_prop_x(running_total);
  std::vector<double> small_system_prop_y(running_total);
  std::vector<double> small_system_prop_z(running_total);
  for (int i = 0; i < running_total; i++) {
    small_system_prop_x[i] = my_prng.uniformRandomNumber();
    small_system_prop_y[i] = my_prng.uniformRandomNumber();
    small_system_prop_z[i] = my_prng.uniformRandomNumber();
  }
  testReduction(small_system_counts, small_system_offsets, small_system_prop_x);
  testReduction(small_system_counts, small_system_offsets, small_system_prop_x,
                small_system_prop_y);
  testReduction(small_system_counts, small_system_offsets, small_system_prop_x,
                small_system_prop_y, small_system_prop_z);
  std::vector<int> large_system_counts(100), large_system_offsets(100);
  running_total = 0;
  for (int i = 0; i < 100; i++) {
    large_system_counts[i]  = 480.0 * (1.5 - my_prng.uniformRandomNumber());
    large_system_offsets[i] = running_total;
    running_total += roundUp(large_system_counts[i], warp_size_int);
  }
  std::vector<double> large_system_prop_x(running_total);
  std::vector<double> large_system_prop_y(running_total);
  std::vector<double> large_system_prop_z(running_total);
  for (int i = 0; i < running_total; i++) {
    large_system_prop_x[i] = my_prng.uniformRandomNumber();
    large_system_prop_y[i] = my_prng.uniformRandomNumber();
    large_system_prop_z[i] = my_prng.uniformRandomNumber();
  }
  testReduction(large_system_counts, large_system_offsets, large_system_prop_x);
  testReduction(large_system_counts, large_system_offsets, large_system_prop_x,
                large_system_prop_y);
  testReduction(large_system_counts, large_system_offsets, large_system_prop_x,
                large_system_prop_y, large_system_prop_z);

  // Test subdivisions of periodic systems
  section(7);
  std::vector<int3> cell_dims(1);
  cell_dims[0] = { 10, 12, 10 };
  Brickwork lego(cell_dims, 5, 16, 1, 0);
  lego.subdivide();
  check(lego.getBrickCount(), RelationalOperator::EQUAL, 54, "The subdivision of a single grid "
        "did not proceed as expected.");
  lego.subdivide(68);
  check(lego.getBrickCount(), RelationalOperator::EQUAL, 68, "The subdivision of a single grid "
        "did not proceed as expected when presented with a larger target number of bricks.");
  verifyBrickCoverage(lego);
  cell_dims.push_back({ 8, 8, 4 });
  cell_dims.push_back({ 18, 17, 14 });
  cell_dims.push_back({ 9, 9, 6 });
  cell_dims.push_back({ 50, 2, 15 });
  Brickwork lego_ii(cell_dims, 5, 16, 1, 0);
  lego_ii.subdivide(68);
  verifyBrickCoverage(lego_ii);
  CHECK_THROWS(Brickwork bad_bw(cell_dims, 8, 12, 1, 0, 4, 15), "A Brickwork with an inadequate "
               "maximum non-halo volume for its cross-sectional area was created.");
  Brickwork lego_iii(cell_dims, 5, 16, 2, 0);
  lego_iii.subdivide(68);
  verifyBrickCoverage(lego_iii);

  // Summary evaluation
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
