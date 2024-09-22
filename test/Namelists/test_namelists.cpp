#include <regex>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "copyright.h"
#include "../../src/Chemistry/chemical_features.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/ForceField/forcefield_enumerators.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Math/statistics.h"
#include "../../src/MoleculeFormat/molecule_format_enumerators.h"
#include "../../src/Namelists/input.h"
#include "../../src/Namelists/nml_conformer.h"
#include "../../src/Namelists/nml_dynamics.h"
#include "../../src/Namelists/nml_ffmorph.h"
#include "../../src/Namelists/nml_files.h"
#include "../../src/Namelists/nml_mesh.h"
#include "../../src/Namelists/nml_minimize.h"
#include "../../src/Namelists/nml_random.h"
#include "../../src/Namelists/nml_receptor.h"
#include "../../src/Namelists/nml_remd.h"
#include "../../src/Namelists/nml_report.h"
#include "../../src/Namelists/nml_restraint.h"
#include "../../src/Namelists/nml_solvent.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/reporting_enumerators.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::chemistry::ChemicalFeatures;
using stormm::constants::ExceptionResponse;
using stormm::constants::tiny;
using stormm::diskutil::osSeparator;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::energy::EvaluateForce;
using stormm::energy::evaluateRestraints;
using stormm::energy::getEnumerationName;
using stormm::energy::ScoreCard;
using stormm::energy::StateVariable;
using stormm::errors::rtWarn;
using stormm::stmath::mean;
using stormm::stmath::variance;
using stormm::stmath::VarianceMethod;
using stormm::modeling::ForceFieldElement;
using stormm::modeling::ParameterKind;
using stormm::parse::separateText;
using stormm::parse::strcmpCased;
using stormm::parse::TextOrigin;
using stormm::review::OutputScope;
using stormm::review::OutputSyntax;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using stormm::restraints::RestraintApparatus;
using stormm::structure::DataRequestKind;
using stormm::topology::AtomGraph;
using stormm::topology::AtomicRadiusSet;
using stormm::topology::ImplicitSolventModel;
using stormm::trajectory::CoordinateFrame;
using namespace stormm::namelist;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Test what should be bad namelist input and make sure that it throws.
//
// Arguments:
//   nml_name:       Name of the namelist to test (i.e. &files)
//   content:        Content of the erroneous namelist.  This function add the trailing &end card.
//   error_message:  Error message to display if the content does NOT throw an exception.  This
//                   function adds "in the (namelist name) namelist" to the end.
//-------------------------------------------------------------------------------------------------
void testBadNamelist(const std::string &nml_name, const std::string &content,
                     const std::string &error_message) {
  const std::string bad_idea = std::string("&") + nml_name + " " + content + " &end";
  const TextFile bad_input(bad_idea, TextOrigin::RAM);
  int start_line = 0;
  bool found_nml;
  
  // Redirect stdout and record the number of test failures thus far.  This is to suppress Alert
  // messages that occur leading up to an ultimate namelist failure, which is supposed to happen
  // and thus should not call for any attention.
  FILE fp_old = *stdout;
  *stdout = *fopen("/dev/null", "w");
  const int initial_failures = gbl_test_results.getFailureCount();
  const std::string updated_error = error_message + " in the &" + nml_name + " namelist.";
  
  // Run the test, which may send messages to stdout
  if (strcmpCased(nml_name, "minimize")) {
    CHECK_THROWS(MinimizeControls t_mincon(bad_input, &start_line, &found_nml), updated_error);
  }
  else if (strcmpCased(nml_name, "dynamics")) {
    CHECK_THROWS(DynamicsControls t_dyncon(bad_input, &start_line, &found_nml), updated_error);
  }
  else if (strcmpCased(nml_name, "random")) {
    CHECK_THROWS(RandomControls t_rngcon(bad_input, &start_line, &found_nml), updated_error);
  }
  else if (strcmpCased(nml_name, "remd")) {
  	CHECK_THROWS(RemdControls t_remcon(bad_input, &start_line, &found_nml), updated_error);
  }
  else if (strcmpCased(nml_name, "solvent")) {
    CHECK_THROWS(SolventControls t_watcon(bad_input, &start_line, &found_nml), updated_error);
  }
  else if (strcmpCased(nml_name, "ffmorph")) {
    CHECK_THROWS(FFMorphControls t_ffmcon(bad_input, &start_line, &found_nml), updated_error);
  }
  else if (strcmpCased(nml_name, "report")) {
    CHECK_THROWS(ReportControls t_repcon(bad_input, &start_line, &found_nml), updated_error);
  }
  else if (strcmpCased(nml_name, "conformer")) {
    CHECK_THROWS(ConformerControls t_confcon(bad_input, &start_line, &found_nml), updated_error);
  }
  else if (strcmpCased(nml_name, "receptor")) {
    CHECK_THROWS(ReceptorControls t_repcon(bad_input, &start_line, &found_nml), updated_error);
  }
  else if (strcmpCased(nml_name, "mesh")) {
    CHECK_THROWS(MeshControls t_meshcon(bad_input, &start_line, &found_nml), updated_error);
  }
  else {
    rtErr("The namelist &" + nml_name + " does not pair with any known case.", "test_namelists");
  }

  // Reset stdout and reprint the error message if there was a new failure.
  *stdout = fp_old;
  if (gbl_test_results.getFailureCount() > initial_failures) {
    const std::string parsed_msg = terminalFormat(error_message + " in the &" + nml_name +
                                                  " namelist.", "", "", 14, 0, 14);
    printf("Check FAILED: %s\n", parsed_msg.c_str());
  }
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1
  section("Test the &files namelist");

  // Section 2
  section("Test the &minimize namelist");

  // Section 3
  section("Test the &random namelist");

  // Section 4
  section("Test the &solvent namelist");

  // Section 5
  section("Test the &ffmorph namelist");

  // Section 6
  section("Test the &restraint namelist");

  // Section 7
  section("The the &conformer namelist");
  
  // Section 8
  section("Test the &report namelist");

  // Section 9
  section("Test the &mesh namelist");

  // Section 10
  section("Test the &receptor namelist");

  // Section 11
  section("Test the &dynamics namelist");

  // Section 12
  section("Test the REMD Namelist");
  
  // The files namelist is perhaps the most complex due to its interchangeable defaults, and
  // will be critical to the operation of any STORMM app
  section(1);
  const char osc = osSeparator();
  const std::string input_base = oe.getStormmSourcePath() + osc + "test" + osc + "Namelists";
  const std::string main_file = input_base + osc + "testrun.in";
  const bool input_exists = (getDrivePathType(main_file) == DrivePathType::FILE);
  const TestPriority do_tests = (input_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;

  // Read the input text from disk and fill in the $STORMM_SOURCE paths to ensure that it works on
  // any system, regardless of the installation and build paths.
  const TextFile tf_disk = (input_exists) ? TextFile(main_file) : TextFile();
  std::string input_str = "";
  for (int i = 0; i < tf_disk.getLineCount(); i++) {
    std::string tfd_line = tf_disk.getLineAsString(i);
    tfd_line = std::regex_replace(tfd_line, std::regex(std::string("test/Namelists/topol/")),
                                  input_base + osc + "topol" + osc);
    tfd_line = std::regex_replace(tfd_line, std::regex(std::string("test/Namelists/coord/")),
                                  input_base + osc + "coord" + osc);
    input_str += tfd_line + "\n";
  }
  const TextFile tf(input_str, TextOrigin::RAM);
  int start_line = 0;
  FilesControls filcon(tf, &start_line);
  const int nfreetop = filcon.getFreeTopologyCount();
  const int nfreecrd = filcon.getFreeCoordinatesCount();
  check(nfreetop, RelationalOperator::EQUAL, 6, "The &files namelist did not record the expected "
        "number of free topologies.", do_tests);
  check(nfreecrd, RelationalOperator::EQUAL, 6, "The &files namelist did not record the expected "
        "number of free input coordinate sets.", do_tests);
  check(filcon.getSystemDefinitionCount(), RelationalOperator::EQUAL, 2, "The &files namelist did "
        "not record the expected number of defined systems.", do_tests);
  std::vector<AtomGraph> ags;
  ags.reserve(nfreetop);
  std::vector<int> top_atom_counts(nfreetop);
  for (int i = 0; i < nfreetop; i++) {
    ags.push_back(AtomGraph(filcon.getFreeTopologyName(i), ExceptionResponse::SILENT));
    top_atom_counts[i] = ags[i].getAtomCount();
  }
  check(sum<int>(top_atom_counts), RelationalOperator::EQUAL, 213, "The free topologies named in "
        "the &files namelist do not contain the expected, combined number of atoms.", do_tests);
  std::vector<CoordinateFrame> cfs;
  ags.reserve(nfreecrd);
  std::vector<int> crd_atom_counts(nfreecrd);
  for (int i = 0; i < nfreecrd; i++) {
    cfs.push_back(CoordinateFrame(filcon.getFreeCoordinateName(i), CoordinateFileKind::UNKNOWN,
                                  0));
    crd_atom_counts[i] = cfs[i].getAtomCount();
  }
  check(sum<int>(crd_atom_counts), RelationalOperator::EQUAL, 213, "The free topologies named in "
        "the &files namelist do not contain the expected, combined number of atoms.", do_tests);

  // The minimization namelist contains integer and real-valued numbers
  std::string bad;
  section(2);
  start_line = 0;
  bool found_nml;
  MinimizeControls mincon(tf, &start_line, &found_nml);
  check(mincon.getTotalCycles(), RelationalOperator::EQUAL, 1000, "The &minimize namelist did "
        "not convey the correct total number of cycles.", do_tests);
  check(mincon.getSteepestDescentCycles(), RelationalOperator::EQUAL, 100, "The &minimize "
        "namelist did not convey the correct number of steepest descent cycles.", do_tests);
  check(mincon.getElectrostaticCutoff(), RelationalOperator::EQUAL, Approx(9.0).margin(tiny),
        "The &minimize namelist did not convey the correct electrostatic cutoff.", do_tests);
  check(mincon.getLennardJonesCutoff(), RelationalOperator::EQUAL, Approx(9.0).margin(tiny),
        "The &minimize namelist did not convey the correct van der-Waals cutoff.", do_tests);
  check(mincon.getInitialStep(), RelationalOperator::EQUAL, Approx(0.02).margin(tiny),
        "The &minimize namelist did not convey the correct initial step size.", do_tests);
  check(mincon.getConvergenceTarget(), RelationalOperator::EQUAL, Approx(0.00008).margin(tiny),
        "The &minimize namelist did not convey the correct convergence criterion.", do_tests);
  testBadNamelist("minimize", "maxcyc = -1", "A negative cycle count was accepted");
  testBadNamelist("minimize", "ncyc = 50, maxcyc = 49", "The number of steepest descent cycles "
                  "was allowed to exceed the number of all cycles");
  testBadNamelist("minimize", "ncyc = 5, maxcyc = 9, dx0 = 0.0", "An initial step size of zero "
                  "was accepted");
  testBadNamelist("minimize", "drms = 0.00000000001", "An exceedingly small convergence step "
                  "was permitted");
  
  // The random number control namelist provides a number of useful options
  section(3);
  start_line = 0;
  RandomControls rngcon(tf, &start_line, &found_nml);
  check(rngcon.getRandomSeed(), RelationalOperator::EQUAL, 67108863, "The &random namelist did "
        "not convey the correct random seed.", do_tests);
  check(rngcon.getStreamCount(), RelationalOperator::EQUAL, 128, "The &random namelist did "
        "not convey the correct stream count.", do_tests);
  check(rngcon.getProductionStride(), RelationalOperator::EQUAL, 32, "The &random namelist did "
        "not convey the correct production stride.", do_tests);
  check(rngcon.getWarmupCycleCount(), RelationalOperator::EQUAL, 63, "The &random namelist did "
        "not convey the correct warmup cycle count.", do_tests);
  testBadNamelist("random", "igseed=0", "A random seed of zero was accepted");
  testBadNamelist("random", "igstreams = -5", "A negative number of random streams was accepted");
  testBadNamelist("random", "igstreams = 0", "A count of zero random streams was accepted");
  testBadNamelist("random", "igstride = -2", "A negative random production stride was accepted");

  // The solvent namelist provides a means of specifying the solvent model
  section(4);
  start_line = 0;
  SolventControls watcon(tf, &start_line, &found_nml);
  check(watcon.getBornRadiiCutoff(), RelationalOperator::EQUAL, Approx(25.0).margin(tiny),
        "The &solvent namelist reports the wrong Born radius cutoff (rgbmax).", do_tests);
  check(watcon.getInternalDielectric(), RelationalOperator::EQUAL, Approx(1.5).margin(tiny),
        "The &solvent namelist reports the wrong internal dieletric (intdiel).", do_tests);
  check(watcon.getExternalDielectric(), RelationalOperator::EQUAL, Approx(77.9).margin(tiny),
        "The &solvent namelist reports the wrong external dieletric (extdiel).", do_tests);
  check(watcon.getPBRadiiSet() == AtomicRadiusSet::MBONDI2, "The &solvent namelist reports the "
        "wrong PB radii set (pbradii)", do_tests);
  check(watcon.getImplicitSolventModel() == ImplicitSolventModel::OBC_GB_II, "The &solvent "
        "namelist reports the wrong implicit solvent model (igb)", do_tests);
  testBadNamelist("solvent", "pbradii = MBONDI4", "An invalid radius set was accepted");
  testBadNamelist("solvent", "pbradii = MBONDI2, pbradii = MBONDI", "Duplicate specification of "
                  "the radius set keyword was allowed");
  testBadNamelist("solvent", "rgbmax = -16.0", "A negative Born radius cutoff was allowed");
  testBadNamelist("solvent", "intdiel = 0.001", "A ridiculously small internal dielectric was "
                  "allowed");
  testBadNamelist("solvent", "extdiel = 0.005", "A ridiculously small external dielectric was "
                  "allowed");
  
  // The force field morphing namelist allows parameter modifications at run time
  section(5);
  start_line = 0;
  FFMorphControls ffmcon(tf, &start_line, &found_nml);
  check(ffmcon.getEditCount(ParameterKind::BOND), RelationalOperator::EQUAL, 3, "The &ffmorph "
        "namelist does not convey the expected number of bond parameter edits.", do_tests);
  check(ffmcon.getEditCount(ParameterKind::ANGLE), RelationalOperator::EQUAL, 2, "The &ffmorph "
        "namelist does not convey the expected number of angle parameter edits.", do_tests);
  check(ffmcon.getEditCount(ParameterKind::DIHEDRAL), RelationalOperator::EQUAL, 3, "The &ffmorph "
        "namelist does not convey the expected number of dihedral parameter edits.", do_tests);
  check(ffmcon.getEditCount(ParameterKind::UREY_BRADLEY), RelationalOperator::EQUAL, 1,
        "The &ffmorph namelist does not convey the expected number of Urey-bradley interaction "
        "edits.", do_tests);
  check(ffmcon.getModelEdit(ParameterKind::BOND, 0).testStiffnessModification(), "The stiffness "
        "parameter of a bond parameter is not properly marked for modification.", do_tests);
  check(ffmcon.getModelEdit(ParameterKind::BOND, 1).testEquilibriumModification() == false,
        "The equilibrium parameter of a bond parameter is not properly marked for modification.",
        do_tests);
  check(ffmcon.getModelEdit(ParameterKind::ANGLE, 1).testStiffnessModification() == false,
        "The stiffness constant of an angle parameter is not properly marked for modification.",
        do_tests);
  check(ffmcon.getModelEdit(ParameterKind::DIHEDRAL, 1).testAmplitudeModification(),
        "The amplitude of a dihedral parameter is not properly marked for modification.",
        do_tests);
  check(ffmcon.getModelEdit(ParameterKind::DIHEDRAL, 1).getAmplitude(), RelationalOperator::EQUAL,
        Approx(2.07).margin(tiny), "The amplitude of a dihedral parameter is not properly marked "
        "for modification.", do_tests);
  check(ffmcon.getModelEdit(ParameterKind::DIHEDRAL, 0).getPhaseAngle(), RelationalOperator::EQUAL,
        Approx(0.24 * stormm::symbols::pi / 180.0).margin(tiny), "The amplitude of a dihedral "
        "parameter is not properly marked for modification.", do_tests);
  testBadNamelist("ffmorph", "bond { -ti CX -tj CT }", "An incomplete bond keyword entry was "
                  "accepted");
  testBadNamelist("ffmorph", "bond { -ti CX -tj CT -blah }", "A bond keyword entry with a "
                  "spurious subkey was accepted");
  testBadNamelist("ffmorph", "angle { -ti CX -tj CT -k 57.0 }", "An angle keyword entry with "
                  "too few atom types was accepted");
  testBadNamelist("ffmorph", "angle { -ti CX -tk CT -theta0 107.8 }", "An angle keyword entry "
                  "with too few atom types was accepted");
  testBadNamelist("ffmorph", "angle { -tj CT -tk CT -theta0 107.8 }", "An angle keyword "
                  "entry with too few atom types was accepted");
  testBadNamelist("ffmorph", "angle { -ti CX -tj CT -tk CT }", "An angle keyword "
                  "entry with no modifications was accepted");
  testBadNamelist("ffmorph", "dihedral { -ti CG -tk CV -tl CN -n 2 -amp 0.97 }", "A dihedral "
                  "keyword entry with too few atom types was accepted");
  testBadNamelist("ffmorph", "dihedral { -ti CG -tk CV -tj OS -tl CN -tl H1 -n 2 -amp 0.97 }",
                  "A dihedral keyword entry with repeated atom types was accepted");

  // The restraint namelist can build individual restraints as well as ensembles of them
  section(6);
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  TestSystemManager tsm(base_top_name, "top",
                        { "stereo_L1_vs", "symmetry_L1_vs", "drug_example_vs_iso",
                          "bromobenzene_vs_iso", "med_1", "med_2", "med_3", "med_4", "med_5" },
                        base_crd_name, "inpcrd",
                        { "stereo_L1_vs", "symmetry_L1_vs", "drug_example_vs_iso",
                          "bromobenzene_vs_iso", "med_1", "med_2", "med_3", "med_4", "med_5" });
  const std::string rst_nml_a("&restraint\n  ensemble heavy_dihedrals\n  penalty 50.0, fbhw 0.0\n"
                              "&end");
  const TextFile nml_tf(rst_nml_a, TextOrigin::RAM);
  start_line = 0;
  RestraintControls rst_ctrl_a(nml_tf, &start_line, nullptr);
  std::vector<int> rst_count(tsm.getSystemCount(), 0);
  std::vector<RestraintApparatus> ra_vec;
  std::vector<double> ra_mean_kvals(tsm.getSystemCount()), ra_std_kvals(tsm.getSystemCount());
  ra_vec.reserve(tsm.getSystemCount());
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    const CoordinateFrame cf_i = tsm.exportCoordinateFrame(i);
    const ChemicalFeatures chemfe_i(tsm.getTopologyPointer(i), cf_i);
    const std::vector<BoundedRestraint> r_i = rst_ctrl_a.getRestraint(tsm.getTopologyPointer(i),
                                                                      chemfe_i, cf_i);
    rst_count[i] = r_i.size();
    ra_vec.emplace_back(r_i, tsm.getTopologyPointer(i));
    std::vector<double>k_vals(r_i.size());
    for (size_t j = 0; j < r_i.size(); j++) {
      k_vals[j] = r_i[j].getStiffness().x;
    }
    ra_mean_kvals[i] = mean(k_vals);
    ra_std_kvals[i] = variance(k_vals, VarianceMethod::STANDARD_DEVIATION);
  }
  const std::vector<int> rst_count_ans = { 88, 36, 36, 8, 47, 124, 51, 90, 81 };
  check(rst_count, RelationalOperator::EQUAL, rst_count_ans, "The number of heavy-atom dihedral "
        "restraints found in each system do not meet expectations.", tsm.getTestingStatus());
  check(ra_mean_kvals, RelationalOperator::EQUAL, std::vector<double>(9, 50.0), "The restraint "
        "stiffnesses were not applied correctly.", tsm.getTestingStatus());
  check(ra_std_kvals, RelationalOperator::EQUAL, std::vector<double>(9, 0.0), "The restraint "
        "stiffnesses were not applied correctly.  All stiffnesses should be identical.",
        tsm.getTestingStatus());
  std::vector<double> rstr_e(tsm.getSystemCount());
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    ScoreCard sc(1, 1, 32);
    PhaseSpace ps_i = tsm.exportPhaseSpace(i);
    evaluateRestraints(ra_vec[i], &ps_i, &sc, EvaluateForce::YES);
    rstr_e[i] = sc.reportTotalEnergy();
  }
  check(rstr_e, RelationalOperator::EQUAL, std::vector<double>(9, 0.0), "Holding restraints "
        "applied to a variety of systems show non-zero initial energies.  This indicates that the "
        "restraint targets, and possibly other aspects of the restraints, do not match the "
        "structures to which they were assigned.", tsm.getTestingStatus());

  // The conformer namelist contains keywords with multiple default values.
  section(7);
  start_line = 0;
  ConformerControls confcon(tf, &start_line, &found_nml);
  check(confcon.getRotationSampleValues().size(), RelationalOperator::EQUAL, 3, "The number of "
        "rotation sampling defaults in the &conformer namelist is incorrect.");
  const std::string conf_nml_a("&conformer\n  trial_limit 10, final_states 6,\n  "
                               "core_mask { atoms \"@N,CA,C & !(:ACE,NME)\" },\n  "
                               "rotation_sample 65.2, rotation_sample 181.8,\n  "
                               "rotation_sample -68.7, rotation_sample 58.3,\n  "
                               "rotation_sample -179.6,\n&end");
  const TextFile confinpa_tf(conf_nml_a, TextOrigin::RAM);
  start_line = 0;
  ConformerControls confcon_a(confinpa_tf, &start_line, nullptr);
  check(confcon_a.getRotationSampleValues().size(), RelationalOperator::EQUAL, 5, "The number of "
        "rotation sampling values in the &conformer namelist is incorrect once more than the "
        "default number have been added.");
  const std::string conf_nml_b("&conformer\n  trial_limit 10, final_states 6,\n  "
                               "core_mask { atoms \"@N,CA,C & !(:ACE,NME)\" },\n  "
                               "rotation_sample 65.0, rotation_sample 178.9,\n&end");
  const TextFile confinpb_tf(conf_nml_b, TextOrigin::RAM);
  start_line = 0;
  ConformerControls confcon_b(confinpb_tf, &start_line, nullptr);
  check(confcon_b.getRotationSampleValues().size(), RelationalOperator::EQUAL, 2, "The number of "
        "rotation sampling values in the &conformer namelist is incorrect when fewer than the "
        "default number have been included.");
  std::vector<double> all_rotation_settings;
  const std::vector<double> rsamp   = confcon.getRotationSampleValues();
  const std::vector<double> rsamp_a = confcon_a.getRotationSampleValues();
  const std::vector<double> rsamp_b = confcon_b.getRotationSampleValues();
  all_rotation_settings.insert(all_rotation_settings.end(), rsamp.begin(), rsamp.end());
  all_rotation_settings.insert(all_rotation_settings.end(), rsamp_a.begin(), rsamp_a.end());
  all_rotation_settings.insert(all_rotation_settings.end(), rsamp_b.begin(), rsamp_b.end());
  std::vector<double> all_rotation_settings_ans = { 60.0, -180.0, -60.0, 65.2, -178.2, -68.7,
                                                    58.3, -179.6, 65.0, 178.9 };
  for (size_t i = 0; i < all_rotation_settings_ans.size(); i++) {
    all_rotation_settings_ans[i] *= stormm::symbols::pi / 180.0;
    const double idiff = all_rotation_settings_ans[i] - all_rotation_settings[i];

    // Correct for possible platform-specific rounding issues that should not change the validity
    // of the rotation angles.
    if (fabs(idiff - stormm::symbols::twopi) < 1.0e-6) {
      all_rotation_settings_ans[i] -= stormm::symbols::twopi;
    }
    else if (fabs(idiff + stormm::symbols::twopi) < 1.0e-6) {
      all_rotation_settings_ans[i] += stormm::symbols::twopi;
    }
  }
  check(all_rotation_settings, RelationalOperator::EQUAL, all_rotation_settings_ans,
        "The rotation angles found for a series of namelists do not meet expectations.", do_tests);
  testBadNamelist("conformer", "rotation_sample_count 5, rotation_sample 150.7", "A &conformer "
                  "namelist specifying five rotation samples and one specific value was "
                  "accepted");
  testBadNamelist("conformer", "cis_trans_sample_count 5, cis_trans_sample 150.7", "A &conformer "
                  "namelist specifying five rotation samples and one specific value was "
                  "accepted");
  testBadNamelist("conformer", "cis_trans_sample_count 4, cis_trans_sample 150.7, "
                  "cis_trans_sample 1.6, cis_trans_sample -0.5, cis_trans_sample 178.2, "
                  "cis_trans_sample 0.2", "A &conformer namelist specifying four rotation samples "
                  "and five specific values was accepted");
  const std::string conf_nml_c("&conformer\n  trial_limit 10, final_states 6,\n  "
                               "core_mask { atoms \"@N,CA,C & !(:ACE,NME)\" },\n  "
                               "rotation_sample 179.0, rotation_sample 179.5,\n  "
                               "rotation_sample 60.7, rotation_sample -59.3,\n  "
                               "cis_trans_sample 181.4\n&end\n");
  const TextFile confinpc_tf(conf_nml_c, TextOrigin::RAM);
  start_line = 0;
  ConformerControls confcon_c(confinpc_tf, &start_line, nullptr);
  check(confcon_c.getCisTransSampleValues(), RelationalOperator::EQUAL,
        std::vector<double>(1, -178.6 * stormm::symbols::pi / 180.0), "The number of cis-trans "
        "sampling values in the &conformer namelist is incorrect when fewer than the default "
        "number have been included.");
  check(confcon_c.getCoreAtomMask(), RelationalOperator::EQUAL, "@N,CA,C & !(:ACE,NME)",
        "The core atom mask of a &conformer namelist is not transcribed correctly.");
  
  // The report namelist needs to be able to convey information for an SD file, among other
  // features of the output.
  section(8);
  const std::string rep_nml_a("&report\n  sdf_item { -title BOND_E -label sulfonamide -energy "
                              "bond }\n  sdf_item { -title ANGLE_E -label sulfonamide -energy "
                              "HarmonicAngle }\n  energy bond\n  energy angle\n  scope full  "
                              "syntax matplotlib &end\n");
  const TextFile repinpa_tf(rep_nml_a, TextOrigin::RAM);
  start_line = 0;
  ReportControls repcon_a(repinpa_tf, &start_line, nullptr);
  check(repcon_a.getReportedQuantityCount(), RelationalOperator::EQUAL, 3, "The number of "
        "reported quantities in a &report namelist does not meet expectations.");
  check(repcon_a.getReportedQuantityCount() >= 3 &&
        repcon_a.getReportedQuantities()[1] == StateVariable::ANGLE &&
        repcon_a.getReportedQuantities()[2] == StateVariable::UREY_BRADLEY, "Details of the "
        "reported quantities in a &report namelist do not meet expectations.");
  check(repcon_a.getSDFileDataRequestCount(), RelationalOperator::EQUAL, 2, "The number of data "
        "requests in a &report namelist is incorrect.");
  check(repcon_a.getOutputSyntax() == OutputSyntax::MATPLOTLIB, "The ouptut format conveyed by a "
        "&report namelist was not correct.");
  check(repcon_a.getOutputScope() == OutputScope::FULL, "The ouptut scope conveyed by a &report "
        "namelist was not correct.");
  testBadNamelist("report", "sdf_item { -title SomeTitle -energy ANGLE }", "A composite energy "
                  "term would be printed to a SD data item");
  testBadNamelist("report", "sdf_item { -energy HarmonicAngle }", "A data item with no title was "
                  "accepted");
  testBadNamelist("report", "sdf_item { -title PrintAngle -parameter HarmonicAngle -typeI CT "
                  "-typeJ CN }", "A data item printing angle force field parameters was "
                  "accepted without the correct number of atom types");
  testBadNamelist("report", "sdf_item { -title PrintBond -parameter bond -typeI CT -typeJ CN "
                  "-typeK CB }", "A data item printing bond force field parameters was accepted "
                  "without the correct number of atom types");
  testBadNamelist("report", "sdf_item { -title PrintAngle -parameter anglep -typeI CT -typeJ CN "
                  "-typeK CB }", "A data item with a nonsensical parameter name was accepted for "
                  "SD file output");
  testBadNamelist("report", "sdf_item { -title PrintDihedral -parameter dihedral -typeI CT "
                  "-typeJ CN -typeK CB -typeL OT }", "A data item printing dihedral force field "
                  "parameters was accepted, but this is a composite of proper, improper, and "
                  "CHARMM improper energy terms unsuitable for a single SD file data item");
  testBadNamelist("report", "sdf_item { -title PrintDihedral -parameter HarmonicAngle -typeI CT "
                  "-typeJ CN -typeK CB -message blah }", "A data item printing angle force field "
                  "parameters was accepted, but there is an extra messaged tacked on");
  testBadNamelist("report", "outlier_sigmas -9.1", "An invalid outlier criterion was accepted");
  
  // Test the automatic correction of an SD file item when conflicting directives are supplied
  const std::string rep_nml_b("&report\n  sdf_item { -title BOND_E -label sulfonamide -energy "
                              "bond }\n  sdf_item { -title ANGLE_E -label sulfonamide -energy "
                              "HarmonicAngle -parameter bond -typeI CA -typeJ NV }\n  energy "
                              "bond\n  energy angle\n&end\n");
  const TextFile repinpb_tf(rep_nml_b, TextOrigin::RAM);
  start_line = 0;
  ReportControls repcon_b(repinpb_tf, &start_line, nullptr, ExceptionResponse::SILENT);
  check(repcon_b.getSDFileDataRequest(1).getKind() == DataRequestKind::STATE_VARIABLE,
        "The subject of an SD file data item request is not correctly conveyed by a &report "
        "namelist.");
  const std::string rep_nml_c("&report\n  sdf_item { -title BOND_E -label sulfonamide -energy "
                              "bond }\n  sdf_item { -title ANGLE_E -label sulfonamide -energy "
                              "HarmonicAngle }\n  energy bond\n  energy angle\n  scope full  "
                              "syntax matplotlib\n  state volume, state temp, state pressure\n"
                              "&end\n");
  const TextFile repinpc_tf(rep_nml_c, TextOrigin::RAM);
  start_line = 0;
  ReportControls repcon_c(repinpc_tf, &start_line, nullptr, ExceptionResponse::SILENT);
  check(repcon_c.getReportedQuantityCount(), RelationalOperator::EQUAL, 6, "The number of "
        "reported quantities, including requested states, in a &report namelist does not meet "
        "expectations.");
  check(repcon_c.getReportedQuantityCount() >= 6 &&
        repcon_c.getReportedQuantities()[4] == StateVariable::VOLUME &&
        repcon_c.getReportedQuantities()[5] == StateVariable::TEMPERATURE_ALL, "Details of the "
        "reported quantities in a &report namelist do not meet expectations.");

  // The receptor namelist must come with a label group for the receptor structures.  Whether that
  // label group actually exists in a separate &files namelist will be tested at runtime.
  section(10);
  const std::string dock_nml_a("&receptor\n  label_group = \"proteins\",\n"
                               "  mesh_position = \"arbitrary\"\n&end\n");
  const TextFile dock_tf_a(dock_nml_a, TextOrigin::RAM);
  start_line = 0;
  ReceptorControls dock_a(dock_tf_a, &start_line, nullptr, ExceptionResponse::SILENT);
  check(dock_a.getLabelGroup(), RelationalOperator::EQUAL, std::string("proteins"), "The label "
        "group for receptor structures was not conveyed correctly by the &receptor namelist.");
  testBadNamelist("receptor", "label_group = unprotected, boundary = plox",
                  "Input was accepted with an unrecognized boundary condition");
  testBadNamelist("receptor", "label_group = unprotected, mesh_position = phlox",
                  "Input was accepted with an unrecognized mesh alignment");
  testBadNamelist("receptor", "label_group = unprotected, potential = blox",
                  "Input was accepted with an unrecognized potential form");
  
  // The mesh namelist groups parameters for the actual mesh, abstracting this functionality for
  // multiple situations in which a mesh might be wanted.
  section(9);
  const std::string mesh_nml_a("&mesh\n  mesh_dim = 64, mesh_dim_b = 48,\n  mesh_spacing = 0.8, "
                               "mesh_spacing_a = 0.75, mesh_alpha = 94.0\n  mesh_origin_x = 17.5, "
                               "mesh_origin_y = 19.4, mesh_origin_z = 14.2\n&end\n");
  const TextFile mesh_tf_a(mesh_nml_a, TextOrigin::RAM);
  start_line = 0;
  MeshControls mesh_a(mesh_tf_a, &start_line, nullptr, ExceptionResponse::SILENT);
  const std::vector<double> mesh_spacings = { mesh_a.getSpacing(UnitCellAxis::A),
                                              mesh_a.getSpacing(UnitCellAxis::B),
                                              mesh_a.getSpacing(UnitCellAxis::C) };
  const std::vector<double> mesh_spacings_ans = { 0.75, 0.8, 0.8 };
  check(mesh_spacings, RelationalOperator::EQUAL, mesh_spacings_ans, "Mesh spacings were not "
        "conveyed as expected by a &mesh namelist.\n");
  const std::vector<double> mesh_angles = { mesh_a.getAlpha(), mesh_a.getBeta(),
                                            mesh_a.getGamma() };
  const std::vector<double> mesh_angles_ans = { 94.0 * stormm::symbols::pi / 180.0,
                                                0.5 * stormm::symbols::pi,
                                                0.5 * stormm::symbols::pi };
  check(mesh_angles, RelationalOperator::EQUAL, mesh_angles_ans, "Mesh angles were not conveyed "
        "as expected by a &mesh namelist.");
  const std::vector<int> mesh_dims = { mesh_a.getAxisElementCount(UnitCellAxis::A),
                                       mesh_a.getAxisElementCount(UnitCellAxis::B),
                                       mesh_a.getAxisElementCount(UnitCellAxis::C) };
  const std::vector<int> mesh_dims_ans = { 64, 48, 64 };
  check(mesh_dims, RelationalOperator::EQUAL, mesh_dims_ans, "Mesh dimensions were not conveyed "
        "as expected by a &mesh namelist.");
  const std::string mesh_nml_b("&mesh\n  mesh_dim = 64, mesh_spacing = 1.0,\n  "
                               "mesh_alpha = 15.0, mesh_beta = 16.1, mesh_gamma = 107.8\n&end\n");
  const TextFile mesh_tf_b(mesh_nml_b, TextOrigin::RAM);
  start_line = 0;
  CHECK_THROWS(MeshControls mesh_b(mesh_tf_b, &start_line, nullptr, ExceptionResponse::DIE),
               "A nonsensical mesh element was cleared for mesh generation.");
  start_line = 0;
  MeshControls mesh_bx(mesh_tf_b, &start_line, nullptr, ExceptionResponse::SILENT);
  const std::vector<double> mesh_adj_dims = { mesh_bx.getAlpha(), mesh_bx.getBeta(),
                                              mesh_bx.getGamma() };
  const std::vector<double> mesh_adj_dims_ans = { 0.877027949, 0.887203219, 1.735450688 };
  check(mesh_adj_dims, RelationalOperator::EQUAL, mesh_adj_dims_ans, "The adjusted angles values "
        "for an impossible mesh element design do not meet expectations.");

  // Test dynamics namelist variables.
  section(11);
  const std::string dynamics_nml_a("&dynamics\n  nstlim = 57, ntpr = 19, ntwx = 19, nscm = 3,\n"
                                   "  dt = 1.5, rigid_geom = off, tol = 5.4e-7, "
                                   "rattle_iter = 45,\n  rattle_style center_sum, ntt = 3, "
                                   "tevo_start = 5, tevo_end = 8\n  vrand = 6, gamma_ln = 0.004, "
                                   "tcache_depth = 6,\n  thermostat_seed = 21858302, "
                                   "tcache_config double\n&end\n");
  const TextFile dyna_tf_a(dynamics_nml_a, TextOrigin::RAM);
  start_line = 0;
  DynamicsControls dyna_a(dyna_tf_a, &start_line, nullptr, ExceptionResponse::SILENT);
  const int step_count_a = dyna_a.getStepCount();
  const int ntpr_a = dyna_a.getDiagnosticPrintFrequency();
  const int ntwx_a = dyna_a.getTrajectoryPrintFrequency();
  const int nscm_a = dyna_a.getCenterOfMassMotionPurgeFrequency();
  const double dt_a = dyna_a.getTimeStep();
  const ApplyConstraints cnst_a = dyna_a.constrainGeometry();
  const double tol_a = dyna_a.getRattleTolerance();
  const int cnst_iter_a = dyna_a.getRattleIterations();
  const RattleMethod cnst_meth_a = dyna_a.getCpuRattleMethod();
  const ThermostatKind thrm_a = dyna_a.getThermostatKind();
  const int tevo_ia = dyna_a.getThermostatEvolutionStart();
  const int tevo_fa = dyna_a.getThermostatEvolutionEnd();
  const int depth_a = dyna_a.getThermostatCacheDepth();
  const int seed_a = dyna_a.getThermostatSeed();
  const PrecisionModel cconfig_a = dyna_a.getThermostatCacheConfig();
  check(step_count_a, RelationalOperator::EQUAL, 57, "The total number of dynamics steps recorded "
        "from a &dynamics namelist does not meet expectations.");
  check(ntpr_a, RelationalOperator::EQUAL, 19, "The diagnostic print frequency recorded from a "
        "&dynamics namelist does not meet expectations.");
  check(ntwx_a, RelationalOperator::EQUAL, 19, "The trajectory print frequency recorded from a "
        "&dynamics namelist does not meet expectations.");
  check(nscm_a, RelationalOperator::EQUAL, 3, "The momentum purge frequency recorded from a "
        "&dynamics namelist does not meet expectations.");
  check(dt_a, RelationalOperator::EQUAL, 1.5, "The time step recorded from a &dynamics namelist "
        "does not meet expectations.");
  check(cnst_a == ApplyConstraints::NO, "The constraint directive recorded from a &dynamics "
        "namelist does not meet expectations.");
  check(tol_a, RelationalOperator::EQUAL, Approx(5.4e-7, ComparisonType::RELATIVE, 2.0e-2),
        "The constraint tolerance recorded from a &dynamics namelist does not meet expectations.");
  check(cnst_iter_a, RelationalOperator::EQUAL, 45, "The constraint maximum iterations setting "
        "recorded from a &dynamics namelist does not meet expectations.");
  check(cnst_meth_a == RattleMethod::CENTER_SUM, "The CPU-based constraint iteration method "
        "recorded from a &dynamics namelist does not meet expectations.");
  check(thrm_a == ThermostatKind::LANGEVIN, "The type of thermostat read from a &dynamics "
        "namelist does not meet expectations.");
  check(tevo_ia, RelationalOperator::EQUAL, 5, "The beginning of temperature evolution read from "
        "a &dynamics namelist does not meet expectations.");
  check(tevo_fa, RelationalOperator::EQUAL, 8, "The conclusion of temperature evolution read from "
        "a &dynamics namelist does not meet expectations.");
  check(depth_a, RelationalOperator::EQUAL, 6, "The random number cache depth recorded from a "
        "&dynamics namelist does not meet expectations.");
  check(seed_a, RelationalOperator::EQUAL, 21858302, "The random number seed recorded from a "
        "&dynamics namelist does not meet expectations.");
  check(cconfig_a == PrecisionModel::DOUBLE, "The random number cache configuration recorded from "
        "a &dynamics namelist does not meet expectations.");
  testBadNamelist("dynamics", "nstlim = -1", "Input was accepted with a nonsensical number of "
                  "steps");
  testBadNamelist("dynamics", "ntpr = -1000", "Input was accepted with a nonsensical diagnostic "
                  "reporting frequency");
  testBadNamelist("dynamics", "tcache_config = \"medium\", ntpr = 50", "Input was accepted with "
                  "an invalid random number cache configuration");
  
  // Testing the REMD Namelist
  section(12);
  const std::string remd_nml_a("&remd\n  total_swaps = 10000, remd_type Temperature,\n"
                               "  freq_swaps = 100, swap_store Successful,\n"
                               "  temp_distribution = \"Van Der Spoel\",\n"
                               "  exchange_probability = 0.2, tolerance = 0.0001,\n"
                               "  max_replicas = 1000, low_temperature = 293.7,\n"
                               "  high_temperature = 393.7\n&end\n");
  const TextFile remd_tf_a(remd_nml_a, TextOrigin::RAM);
  start_line = 0;
  RemdControls remd_a(remd_tf_a, &start_line, nullptr, ExceptionResponse::SILENT);
  
  const int total_steps = remd_a.getTotalSwapCount();
  const std::string remd_type = remd_a.getRemdType();
  const int freq_swaps = remd_a.getFrequencyOfSwaps();
  const std::string swap_storage = remd_a.getSwapStore();
  const std::string temp_dist = remd_a.getTemperatureDistributionMethod();
  const double exchange_probability = remd_a.getExchangeProbability();
  const double tolerance = remd_a.getTolerance();
  const int max_replicas = remd_a.getMaxReplicas();
  const double low_temperature = remd_a.getInitialTemperature();
  const double high_temperature = remd_a.getEquilibriumTemperature();
  check(total_steps, RelationalOperator::EQUAL, 10000, "The total number of swaps recorded from "
        "the &remd namelist do not meet expectations.");
  check(remd_type, RelationalOperator::EQUAL, "Temperature", "The type of REMD recorded from the "
        "&remd namelist do not meet expectations.");
  check(freq_swaps, RelationalOperator::EQUAL, 100, "The number of frequency of swaps recorded "
        "from the &remd namelist does not meet expectations.");
  check(swap_storage, RelationalOperator::EQUAL, "Successful", "The Swap Storage method recorded "
        "from the &remd namelist does not meet expectations.");
  check(temp_dist, RelationalOperator::EQUAL, "Van Der Spoel", "The temperature distribution "
        "algorithm recorded from the &remd namelist does not meet expectations.");
  check(exchange_probability, RelationalOperator::EQUAL, 0.2, "The Exchange Probability recorded "
        "from the &remd namelist does not meet expectations.");
  check(tolerance, RelationalOperator::EQUAL, 0.0001, "The tolerance recorded from the &remd "
        "namelist does not meet expectations.");
  check(max_replicas, RelationalOperator::EQUAL, 1000, "The maximum number of replicas recorded "
        "from the &remd namelist does not meet expectations.");
  check(low_temperature, RelationalOperator::EQUAL, 293.7, "The low temperatures recorded "
        "from the &remd namelist does not meet expectations.");
  check(high_temperature, RelationalOperator::EQUAL, 393.7, "The high temperatures recorded "
        "from the &remd namelist does not meet expectations.");
  testBadNamelist("remd", "exchange_probability = 2.1", "Input was accepted with a nonsensical "
                  "exchange probability.");
  testBadNamelist("remd", "max_replicas = -100", "Input was accepted with a nonsensical max "
                  "replica count.");
  testBadNamelist("remd", "freq_swaps = -100", "Input was accepted with a nonsensical frequency "
                  "of steps to attempt before attempting to do a swap.");
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
