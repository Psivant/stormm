#include "copyright.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/ForceField/forcefield_element.h"
#include "../../src/ForceField/forcefield_enumerators.h"
#include "../../src/MolecularMechanics/minimization.h"
#include "../../src/MolecularMechanics/mm_evaluation.h"
#include "../../src/Namelists/nml_minimize.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/static_exclusionmask.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Topology/atomgraph_stage.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::diskutil;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::mm;
using namespace stormm::modeling;
using namespace stormm::namelist;
using namespace stormm::parse;
using namespace stormm::review;
using namespace stormm::topology;
using namespace stormm::trajectory;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Initialize the test environment
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1: Test bond modifications
  section("Test bond parameter modifications");

  // Section 2: Test angle modifications
  section("Test angle parameter modifications");

  // Section 3: Test dihedral modifications
  section("Test dihedral parameter modifications");

  // Section 4: Test topology staging
  section("Test topology staging");
  
  // Load a series of topologies, then modify them a bit at a time
  const char osc = osSeparator();
  const std::string all_base_name   = oe.getStormmSourcePath() + osc + "test";
  const std::string topology_home   = "Topology";
  const std::string trajectory_home = "Trajectory";
  const std::string chemistry_home  = "Chemistry";
  const std::vector<std::string> top_names = {
    topology_home + osc + "stereo_L1", topology_home + osc + "stereo_L1_vs",
    topology_home + osc + "symmetry_L1", topology_home + osc + "symmetry_L1_vs",
    topology_home + osc + "bromobenzene", topology_home + osc + "bromobenzene_vs",
    chemistry_home + osc + "lig1_c8h8", chemistry_home + osc + "lig2_c8h8",
    chemistry_home + osc + "lig3_1fsj", chemistry_home + osc + "morphine_like"
  };
  const std::vector<std::string> crd_names = {
    trajectory_home + osc + "stereo_L1", trajectory_home + osc + "stereo_L1_vs",
    trajectory_home + osc + "symmetry_L1", trajectory_home + osc + "symmetry_L1_vs",
    trajectory_home + osc + "bromobenzene", trajectory_home + osc + "bromobenzene_vs",
    chemistry_home + osc + "lig1_c8h8", chemistry_home + osc + "lig2_c8h8",
    chemistry_home + osc + "lig3_1fsj", chemistry_home + osc + "morphine_like",
  };
  TestSystemManager tsm(all_base_name, "top", top_names, all_base_name, "inpcrd", crd_names);
  const TestPriority do_tests = tsm.getTestingStatus();
  std::vector<AtomGraph> ag_list;
  std::vector<PhaseSpace> ps_list;
  std::vector<StaticExclusionMask> se_list;
  const int system_count = tsm.getSystemCount();
  ag_list.reserve(system_count);
  ps_list.reserve(system_count);
  se_list.reserve(system_count);
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    ag_list.push_back(tsm.exportAtomGraph(i));
    ps_list.push_back(tsm.exportPhaseSpace(i));
    se_list.emplace_back(StaticExclusionMask(ag_list[i]));
  }

  // Evaluate energies with the canonical topologies
  ScoreCard sc_orig(system_count);
  MinimizeControls mincon;
  mincon.setSteepestDescentCycles(25);
  mincon.setTotalCycles(50);
  for (size_t i = 0; i < system_count; i++) {
    RestraintApparatus ra(&ag_list[i]);
    ScoreCard min_sc = minimize(&ps_list[i], &ag_list[i], se_list[i], mincon, 30);
    evalNonbValeMM(&ps_list[i], &sc_orig, ag_list[i], se_list[i], EvaluateForce::NO, i);
  }
  const std::vector<double> bond_e0 = sc_orig.reportInstantaneousStates(StateVariable::BOND);
  const std::vector<double> angl_e0 = sc_orig.reportInstantaneousStates(StateVariable::ANGLE);
  const std::vector<double> dihe_e0 =
    sc_orig.reportInstantaneousStates(StateVariable::PROPER_DIHEDRAL);
  const std::vector<double> impr_e0 =
    sc_orig.reportInstantaneousStates(StateVariable::IMPROPER_DIHEDRAL);
  ForceFieldElement ca_ha_bond(ParameterKind::BOND, stringToChar4("ca"), stringToChar4("ha"));
  ForceFieldElement bl_bm_bond(ParameterKind::BOND, stringToChar4("bl"), stringToChar4("bm"));
  ca_ha_bond.setStiffness(100.0);
  ca_ha_bond.setEquilibrium(1.5);
  bl_bm_bond.setStiffness(120.0);
  bl_bm_bond.setEquilibrium(1.7);
  std::vector<AtomGraph> ag_mods;
  std::vector<PhaseSpace> ps_mods;
  for (size_t i = 0; i < system_count; i++) {
    ag_mods.emplace_back(ag_list[i]);
    ps_mods.emplace_back(ps_list[i]);
  }
  ScoreCard sc_mods(system_count);
  for (size_t i = 0; i < system_count; i++) {
    ca_ha_bond.apply(&ag_mods[i]);
    bl_bm_bond.apply(&ag_mods[i]);
    RestraintApparatus ra(&ag_mods[i]);
    ScoreCard min_sc = minimize(&ps_mods[i], &ag_mods[i], se_list[i], mincon, 30);
    evalNonbValeMM(&ps_mods[i], &sc_mods, ag_mods[i], se_list[i], EvaluateForce::NO, i);
  }
  const std::vector<double> bond_em = sc_mods.reportInstantaneousStates(StateVariable::BOND);
  const std::vector<double> angl_em= sc_mods.reportInstantaneousStates(StateVariable::ANGLE);
  const std::vector<double> dihe_em =
    sc_mods.reportInstantaneousStates(StateVariable::PROPER_DIHEDRAL);
  const std::vector<double> impr_em =
    sc_mods.reportInstantaneousStates(StateVariable::IMPROPER_DIHEDRAL);

  // Try deconstructing a topology and putting it back together
  section(4);
  AtomGraphStage ags(ag_list[0]);
  AtomGraph nchiral_sys = ags.exportTopology();
  check(nchiral_sys.getAtomCount(), RelationalOperator::EQUAL, ag_list[0].getAtomCount(),
        "The AtomGraphStage does not convey the atom count correctly when initialized based on an "
        "existing topology.", do_tests);
  AtomGraphStage ags_vs(ag_list[5]);
  AtomGraph nbromo_vs_sys = ags_vs.exportTopology();
  check(nbromo_vs_sys.getAtomCount(), RelationalOperator::EQUAL, ag_list[5].getAtomCount(),
        "The AtomGraphStage does not convey the atom count correctly when initialized based on an "
        "existing topology with virtual sites.", do_tests);
  AtomGraphStage ion_web(100, { 0, 10, 20, 30, 40, 50, 100 });
  AtomGraph nion_sys = ion_web.exportTopology();
  check(nion_sys.getAtomCount(), RelationalOperator::EQUAL, 100, "The AtomGraphStage does not "
        "convey the atom count correctly when initialized based on a number of particles.");
  
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
