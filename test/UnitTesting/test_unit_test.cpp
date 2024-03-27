#include <string>
#include "copyright.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_util.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/series_ops.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Parsing/parsing_enumerators.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/coordinate_series.h"
#include "../../src/UnitTesting/dissect_textfile.h"
#include "../../src/UnitTesting/file_snapshot.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/Random/random.h"

using stormm::constants::warp_size_int;
#ifndef STORMM_USE_HPC
using stormm::data_types::int4;
#endif
using stormm::errors::rtWarn;
using stormm::stmath::incrementingSeries;
using stormm::parse::polyNumericVector;
using stormm::parse::TextFile;
using stormm::parse::TextOrigin;
using stormm::random::Ran2Generator;
using stormm::random::Xoshiro256ppGenerator;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using namespace stormm::diskutil;
using namespace stormm::trajectory;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv, TmpdirStatus::REQUIRED);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1
  section("Structs and methods");
  
  // Section 2
  section("Snapshotting");

  // Section 3
  section("Test system manager");
  
  // Section 4
  section("Timings");

  // Section 5
  section("Text file parsing");
  
  // Perform basic checks of the Approx object
  StopWatch section_timer("Trial Timer");
  section_timer.addCategory(unitTestSectionName(1));
  section(1);
  Approx gray_number(9.6, ComparisonType::ABSOLUTE, 0.11);
  check(gray_number.test(9.7), "Approx object fails to perform a real-to-real scalar comparison.");
  Approx gray_vector(std::vector<double>{9.8, 8.1, 8.0, 4.5, 7.3}, ComparisonType::ABSOLUTE,
                     0.011);
  check(gray_vector.test(9.7) == false, "A scalar-to-vector approximate comparison was judged "
        "successful.");
  check(gray_vector.test(std::vector<double>{9.79, 8.11, 7.99, 4.49, 7.31}), "Approx object fails "
        "to pass a correct vector-to-vector comparison.");
  check(gray_vector.test(std::vector<double>{9.79, 8.08, 7.99, 4.49, 7.31}) == false,
        "Approx object fails to reject an unacceptable real vector-to-vector comparison.");
  check(gray_vector.test(std::vector<double>{9.8, 8.1, 8.0}) == false, "Approx object performs a "
        "comparison on vectors of different lengths and returns success.");
  check(gray_number, RelationalOperator::GREATER_THAN, 9.4, "Approximate comparison fails to "
        "process a scalar greater-than inequality.");
  check((gray_number > 9.8) == false, "Approximate comparison fails to properly declare a "
        "scalar greater-than inequality to be false.");
  check(4.5, RelationalOperator::LESS_THAN, gray_number, "Approximate comparison fails to "
        "process a scalar less-than inequality.");
  check((9.8 < gray_number) == false, "Approximate comparison fails to properly declare a "
        "scalar less-than inequality to be false.");
  check(4, RelationalOperator::GREATER_THAN_OR_EQUAL, Approx(4), "Approximate comparison fails "
        "to identify a true greater-than-or-equal scalar comparison.");
  check(6, RelationalOperator::GREATER_THAN_OR_EQUAL, Approx(4), "Approximate comparison fails "
        "to identify a true greater-than-or-equal scalar comparison.");
  check(4, RelationalOperator::LESS_THAN_OR_EQUAL, Approx(4), "Approximate comparison fails "
        "to identify a true less-than-or-equal scalar comparison.");
  check(6, RelationalOperator::LESS_THAN_OR_EQUAL, Approx(8), "Approximate comparison fails "
        "to identify a true less-than-or-equal scalar comparison.");
  check(std::vector<double>{9.79, 8.08, 7.99, 4.49, 7.31}, RelationalOperator::LE,
        gray_vector, "Approximate less-than-or-equal comparisons should give the benefit of the "
        "doubt, and let tolerances work in favor of any inequality comparisons.  This was not the "
        "case for a comparison of two vectors.");
  check((std::vector<double>{11.79, 8.08, 7.99, 4.49, 7.31} < gray_vector) == false,
        "Approximate comparison fails to properly identify a less-than comparison of two vectors "
        "as false.  Any element violating the inequality should render the entire comparison "
        "false.");
  check(std::vector<double>{11.79, 80.08, 71.99, 40.49, 17.31}, RelationalOperator::GREATER_THAN,
        gray_vector, "Approximate comparison fails to process a greater-than comparison of two "
        "vectors.");
  check(6, RelationalOperator::GREATER_THAN_OR_EQUAL, Approx(4), "Approximate comparison fails "
        "to identify a true greater-than-or-equal scalar comparison.");
  check(4, RelationalOperator::LESS_THAN_OR_EQUAL, Approx(4), "Approximate comparison fails "
        "to identify a true less-than-or-equal scalar comparison.");
  check(6, RelationalOperator::LESS_THAN_OR_EQUAL, Approx(8), "Approximate comparison fails "
        "to identify a true less-than-or-equal scalar comparison.");
  check(908153, RelationalOperator::GE, 908153, "Inferred approximate comparison is too harsh in "
        "comparing two medium-sized integer values by greater-than-or-equal inequality.");
  check(908153, RelationalOperator::LE, 908153, "Inferred approximate comparison is too harsh in "
        "comparing two medium-sized integer values by less-than-or-equal inequality.");
  check(-90153, RelationalOperator::LT, 90153, "Inferred approximate comparison misses an obvious "
        "comparison of a significant number with its own negative value.");
  const std::vector<double> short_vector = {11.79, 80.08, 71.99, 40.49};
  check((short_vector > gray_vector) == false, "Approximate comparison processes a greater-than "
        "comparison of two vectors with different lengths and returns success.");
  check((short_vector < gray_vector) == false, "Approximate comparison processes a less-than "
        "comparison of two vectors with two different lengths and returns success.");
  check((short_vector >= gray_vector) == false, "Approximate comparison processes a "
        "greater-than-or-equal comparison of two vectors with different lengths, returning "
        "success.");
  check((short_vector <= gray_vector) == false, "Approximate comparison processes a "
        "less-than-or-equal comparison of two vectors, returning success despite them being of "
        "different lengths.");
  section_timer.assignTime(1);

  // Examine file snapshotting by first recording data, then re-reading it.
  section_timer.addCategory(unitTestSectionName(2));
  section(2);
  const int n_pts = 100;
  Ran2Generator prng(71277);
  std::vector<double> x(n_pts, 0.0);
  std::vector<double> y(n_pts, 0.0);
  std::vector<double> z(n_pts, 0.0);
  for (int i = 0; i < n_pts; i++) {
    x[i] = prng.gaussianRandomNumber();
    y[i] = prng.gaussianRandomNumber();
    z[i] = prng.gaussianRandomNumber();
  }
  TestPriority snp_tests = (oe.getTemporaryDirectoryAccess()) ? TestPriority::CRITICAL :
                                                                TestPriority::ABORT;
  if (snp_tests == TestPriority::ABORT) {
    rtWarn("The temporary directory " + oe.getTemporaryDirectoryPath() + " is not writeable and "
           "will therefore not support some subsequent tests.  Make sure that the $STORMM_TMPDIR "
           "environment variable is set to a directory where you have write permissions.",
           "test_unit_test");
  }
  const std::string snp_file = oe.getTemporaryDirectoryPath() + osSeparator() + "xyz_randoms.m";
  snapshot(snp_file, polyNumericVector(x), "x_randoms", 1.0e-4, "Failed to record a shapshot file "
           "when explicitly told to do so.", SnapshotOperation::SNAPSHOT, 1.0e-8,
           NumberFormat::STANDARD_REAL, PrintSituation::OPEN_NEW, snp_tests);
  snapshot(snp_file, polyNumericVector(y), "y_randoms", 1.0e-4, "Failed to record a shapshot file "
           "when explicitly told to do so.", SnapshotOperation::SNAPSHOT, 1.0e-8,
           NumberFormat::STANDARD_REAL, PrintSituation::APPEND, snp_tests);
  snapshot(snp_file, polyNumericVector(z), "z_randoms", 1.0e-4, "Failed to record a shapshot file "
           "when explicitly told to do so.", SnapshotOperation::SNAPSHOT, 1.0e-8,
           NumberFormat::STANDARD_REAL, PrintSituation::APPEND, snp_tests);
  snapshot(snp_file, polyNumericVector(x), "x_randoms", 1.0e-4, "Failed to read x_randoms from " +
           snp_file + ".", SnapshotOperation::COMPARE, 1.0e-8, NumberFormat::STANDARD_REAL,
           PrintSituation::APPEND, snp_tests);
  snapshot(snp_file, polyNumericVector(y), "y_randoms", 1.0e-4, "Failed to read y_randoms from " +
           snp_file + ".", SnapshotOperation::COMPARE, 1.0e-8, NumberFormat::STANDARD_REAL,
           PrintSituation::APPEND, snp_tests);
  snapshot(snp_file, polyNumericVector(z), "z_randoms", 1.0e-4, "Failed to read z_randoms from " +
           snp_file + ".", SnapshotOperation::COMPARE, 1.0e-8, NumberFormat::STANDARD_REAL,
           PrintSituation::APPEND, snp_tests);
  CHECK_THROWS_SOFT(snapshot(snp_file, polyNumericVector(x), "x_randoms", 1.0e-4,
                             "This should not work.", SnapshotOperation::SNAPSHOT, 1.0e-8,
                             NumberFormat::STANDARD_REAL, PrintSituation::OPEN_NEW),
                    "Snapshotting was able to open a new file " + snp_file +
                    " despite it already existing.", snp_tests);
  const std::string snp_txt_file = oe.getTemporaryDirectoryPath() + osSeparator() + "verbiage.txt";
  const std::string some_text("# Text written to a snapshot file can take any form.\n"
                              "The only restriction is that it not contain special markers\n"
                              "beginning with '|>>>' and followed by a label name or 'End'.\n");
  snapshot(snp_txt_file, some_text, "ktxt", "A text block was not properly conveyed to or from a "
           "snapshot file.", SnapshotOperation::SNAPSHOT, PrintSituation::OPEN_NEW, snp_tests);
  snapshot(snp_txt_file, some_text, "ktxt", "A text block was not properly conveyed to or from a "
           "snapshot file.", SnapshotOperation::COMPARE, PrintSituation::OPEN_NEW, snp_tests);
  if (snp_tests == TestPriority::CRITICAL) {
    oe.logFileCreated(snp_file);
    oe.logFileCreated(snp_txt_file);
  }
  section_timer.assignTime(2);

  // Test the test system manager
  section(3);
  const char osc = osSeparator();
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string top_name_ext = "top";
  const std::string crd_name_ext = "inpcrd";
  const std::vector<std::string> sysnames = { "symmetry_L1", "stereo_L1", "med_1", "med_4" };
  TestSystemManager tsm(base_top_name, top_name_ext, sysnames, base_crd_name, crd_name_ext,
                        sysnames);
  std::vector<std::string> all_top_names, all_crd_names;
  all_top_names.reserve(sysnames.size());
  all_crd_names.reserve(sysnames.size());
  for (size_t i = 0; i < sysnames.size(); i++) {
    all_top_names.push_back(base_top_name + osc + sysnames[i] + "." + top_name_ext);
    all_crd_names.push_back(base_crd_name + osc + sysnames[i] + "." + crd_name_ext);
  }
  bool systems_exist = true;
  for (size_t i = 0; i < sysnames.size(); i++) {
    systems_exist = (systems_exist &&
                     getDrivePathType(all_top_names[i]) == DrivePathType::FILE &&
                     getDrivePathType(all_crd_names[i]) == DrivePathType::FILE);
  }
  const TestPriority do_tests = (systems_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  check(tsm.getTestingStatus() == do_tests, "A hand-coded search of the files composing the "
        "TestSystemManager object does not produce the same testing status.");
  check(tsm.getTopologyFile(1), RelationalOperator::EQUAL, all_top_names[1],
        "The TestSystemManager does not report the expected name for one of its topology files.",
        tsm.getTestingStatus());
  check(tsm.getCoordinateFile(1), RelationalOperator::EQUAL, all_crd_names[1],
        "The TestSystemManager does not report the expected name for one of its topology files.",
        tsm.getTestingStatus());
  check(tsm.getSystemCount(), RelationalOperator::EQUAL, 4, "The number of systems in the "
        "TestSystemManager is incorrect.", tsm.getTestingStatus());
  const CoordinateFrame cf_one = tsm.exportCoordinateFrame(1);
  check(cf_one.getAtomCount(), RelationalOperator::EQUAL, 78, "The TestSystemManager does not "
        "export CoordinateFrame objects correctly.", tsm.getTestingStatus());
  const PhaseSpace ps_two = tsm.exportPhaseSpace(2);
  check(ps_two.getAtomCount(), RelationalOperator::EQUAL, 44, "The TestSystemManager does not "
        "export PhaseSpace objects correctly.", tsm.getTestingStatus());
  const PhaseSpaceSynthesis all_coords = tsm.exportPhaseSpaceSynthesis({ 0, 1, 2, 3 }, 0.0,
                                                                       1, 24);
  const CoordinateFrame cf_two = all_coords.exportCoordinates(2);
  const std::vector<double> cf_two_xyz = cf_two.getInterlacedCoordinates();
  const std::vector<double> ps_two_xyz = ps_two.getInterlacedCoordinates();
  const double roundoff_err = meanUnsignedError(cf_two_xyz, ps_two_xyz);
  check(roundoff_err, RelationalOperator::LESS_THAN, 5.0e-8, "The TestSystemManager's synthesis "
        "products contain unexpected roundoff error, or perhaps the wrong order of systems.",
        tsm.getTestingStatus());
  const double perturbation_sigma = 0.1;
  const int igseed = 915087;
  const CoordinateSeries<float> cs_three = tsm.exportCoordinateSeries<float>(3, 4,
                                                                             perturbation_sigma,
                                                                             igseed);
  Xoshiro256ppGenerator main_xrs(igseed);
  CoordinateSeries<double> cs_three_main(tsm.exportCoordinateFrame(3), 4);
  CoordinateSeriesWriter<double> csr = cs_three_main.data();
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < csr.natom; j++) {
      const size_t atom_ij = (i * roundUp(csr.natom, warp_size_int)) + j;
      csr.xcrd[atom_ij] += 0.1 * main_xrs.gaussianRandomNumber();
      csr.ycrd[atom_ij] += 0.1 * main_xrs.gaussianRandomNumber();
      csr.zcrd[atom_ij] += 0.1 * main_xrs.gaussianRandomNumber();
    }
  }
  const CoordinateFrame cs_outcome      = cs_three.exportFrame(1);
  const CoordinateFrame cs_outcome_main = cs_three_main.exportFrame(1);
  const std::vector<double> frm_one_xyz      = cs_outcome.getInterlacedCoordinates();
  const std::vector<double> frm_one_xyz_main = cs_outcome_main.getInterlacedCoordinates();
  const double noise_roundoff = meanUnsignedError(frm_one_xyz, frm_one_xyz_main);
  check(noise_roundoff, RelationalOperator::LESS_THAN, 1.0e-7, "The TestSystemManager does not "
        "introduce the expected noise into a CoordinateSeries when requested.",
        tsm.getTestingStatus());
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis({ 0, 1, 2, 3, 0, 2, 3, 1 });
  check(poly_ag.getSystemCount(), RelationalOperator::EQUAL, 8, "The atomgraph synthesis created "
        "by the TestSystemManager object does not contain the expected number of systems.",
        tsm.getTestingStatus());
  check(poly_ag.getUniqueTopologyCount(), RelationalOperator::EQUAL, 4, "The atomgraph synthesis "
        "created by the TestSystemManager object does not contain the expected number of unique "
        "topologies.", tsm.getTestingStatus());
  check(poly_ag.getSystemTopologyPointer(6)->getAtomCount(), RelationalOperator::EQUAL,
        tsm.getTopologyPointer(3)->getAtomCount(), "The topologies present in the synthesis "
        "created by the TestSystemManager do not meet expectations.", tsm.getTestingStatus());

  // Try making a system manager with missing components.   The tests below rely on the original,
  // complete tsm variable to execute as this manager contains all of the components that *should*
  // be present.
  std::vector<std::string> bad_sysnames = sysnames;
  bad_sysnames.push_back("this_does_not_exist.nope");
  TestSystemManager bad_tsm(base_top_name, top_name_ext, bad_sysnames, base_crd_name, crd_name_ext,
                            bad_sysnames, ExceptionResponse::SILENT);
  check(bad_tsm.getTestingStatus() == TestPriority::ABORT, "A TestSystemManager with missing "
        "components returns TestPriority::" + getEnumerationName(bad_tsm.getTestingStatus()) +
        " when it should signal " + getEnumerationName(TestPriority::ABORT) + ".",
        tsm.getTestingStatus());
  check(bad_tsm.getTestingStatus(3) == TestPriority::CRITICAL, "A TestSystemManager with missing "
        "components returns " + getEnumerationName(bad_tsm.getTestingStatus(3)) + " for a "
        "specific system, but should return " + getEnumerationName(TestPriority::CRITICAL) + ".",
        tsm.getTestingStatus());
  TestSystemManager ign_tsm(base_top_name, top_name_ext, bad_sysnames, base_crd_name, crd_name_ext,
                            bad_sysnames, ExceptionResponse::SILENT, TestPriority::NON_CRITICAL);
  check(ign_tsm.getTestingStatus(4) == TestPriority::NON_CRITICAL, "A TestSystemManager with "
        "missing components returns " + getEnumerationName(bad_tsm.getTestingStatus(4)) +
        "for tests on a missing system, but should return " +
        getEnumerationName(TestPriority::NON_CRITICAL) + ".", tsm.getTestingStatus());

  // Test synthesis production from the TestSystemManager
  std::vector<std::string> small_mols = { "symmetry_C1", "bromobenzene_iso", "bromobenzene_vs",
                                          "drug_example_dry", "drug_example", "trpcage",
                                          "trpcage_in_water", "ubiquitin", "med_1",
                                          "symmetry_C2_in_water", "med_3",
                                          "symmetry_C3_in_water", "dna" };
  TestSystemManager smol_tsm(base_top_name, top_name_ext, small_mols, base_crd_name, crd_name_ext,
                             small_mols, ExceptionResponse::SILENT);
  const std::vector<UnitCellType> pbc_unit_cells = { UnitCellType::ORTHORHOMBIC,
                                                     UnitCellType::TRICLINIC };
  AtomGraphSynthesis all_pbc = smol_tsm.exportAtomGraphSynthesis(pbc_unit_cells);
  AtomGraphSynthesis all_iso = smol_tsm.exportAtomGraphSynthesis(UnitCellType::NONE);
  check(all_pbc.getSystemCount(), RelationalOperator::EQUAL, 7, "The number of systems with "
        "periodic boundary conditions extracted from a TestSystemManager does not meet "
        "expectations.", smol_tsm.getTestingStatus());
  check(all_iso.getSystemCount(), RelationalOperator::EQUAL, 6, "The number of systems with "
        "isolated boundary conditions extracted from a TestSystemManager does not meet "
        "expectations.", smol_tsm.getTestingStatus());
  PhaseSpaceSynthesis all_pbc_sys = smol_tsm.exportPhaseSpaceSynthesis(pbc_unit_cells);
  PhaseSpaceSynthesis all_iso_sys = smol_tsm.exportPhaseSpaceSynthesis(UnitCellType::NONE);
  check(all_pbc.getSystemCount(), RelationalOperator::EQUAL, all_pbc_sys.getSystemCount(),
        "The coordinate and topology syntheses produced by a TestSystemManager have dissimilar "
        "numbers of systems in periodic boundary conditions.", smol_tsm.getTestingStatus());
  check(all_iso.getSystemCount(), RelationalOperator::EQUAL, all_iso_sys.getSystemCount(),
        "The coordinate and topology syntheses produced by a TestSystemManager have dissimilar "
        "numbers of systems in isolated boundary conditions.", smol_tsm.getTestingStatus());
  std::vector<int> pbc_atom_counts_ag(all_pbc.getSystemCount());
  std::vector<int> pbc_atom_counts_ps(all_pbc.getSystemCount());
  for (size_t i = 0; i < all_pbc.getSystemCount(); i++) {
    pbc_atom_counts_ag[i] = (all_pbc.getSystemTopologyPointer(i)->getAtomCount());
    pbc_atom_counts_ps[i] = (all_pbc_sys.getSystemTopologyPointer(i)->getAtomCount());
  }
  check(pbc_atom_counts_ag, RelationalOperator::EQUAL, pbc_atom_counts_ps, "The number of atoms "
        "in systems with periodic boundary conditions differ in the topology and coordinate "
        "syntheses.", smol_tsm.getTestingStatus());
  CHECK_THROWS_SOFT(AtomGraphSynthesis bad_poly_ag =
                    smol_tsm.exportAtomGraphSynthesis(incrementingSeries(0, 13)), "A topology "
                    "synthesis containing systems with mismatched boundary conditions was created "
                    "by the TestSystemManager.",  smol_tsm.getTestingStatus());

  // Test the benchmarking 
  section_timer.addCategory(unitTestSectionName(3) + ", Part A");
  section(4);
  double t = 0.0;
  for (int i = 0; i < 7500 * n_pts; i++) {
    t += prng.gaussianRandomNumber();
  }
  section_timer.assignTime(3);
  section_timer.addCategory(unitTestSectionName(3) + ", Part B");
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 7500 * n_pts; j++) {
      t += prng.gaussianRandomNumber();
    }
    section_timer.assignTime(4);
  }
  section_timer.addCategory(unitTestSectionName(3) + ", Part C");
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 15000 * n_pts; j++) {
      t += prng.gaussianRandomNumber();
    }
    section_timer.assignTime(5);
  }
  check(section_timer.getCategorySamples(4), RelationalOperator::EQUAL, 4, "An incorrect number "
        "of samples is recorded for " + section_timer.getCategoryName(4) + ".");
  check(section_timer.getCategorySamples(3), RelationalOperator::EQUAL, 1, "An incorrect number "
        "of samples is recorded for " + section_timer.getCategoryName(3) + ".");
  check(section_timer.getCategoryMinimumTime(3), RelationalOperator::EQUAL,
        section_timer.getCategoryMaximumTime(3), "The minimum and maximum intervals recorded for "
        "a timing category with only one entry should be identical.");
  check(section_timer.getCategoryAverageInterval(5) /
        section_timer.getCategoryAverageInterval(4), RelationalOperator::EQUAL,
        Approx(2.0).margin(0.1), "Timings should roughly double for producing twice the quantity "
        "of random numbers.", TestPriority::NON_CRITICAL);

  // Test the finer points of text file parsing
  section(5);
  const std::string test_tfstr("Two gluggish dworps gleeped bretly in a thiggy glipment.  They "
                               "really gleeped, and it was bretly.");
  const TextFile test_tf(test_tfstr, TextOrigin::RAM);
  const std::vector<char> text_separators = { ',', '.' };
  const std::vector<char> text_delimiters = { '[', ']', '(', ')', '\n' };
  const int n_two = countInstances(test_tf, "Two", text_separators, text_delimiters);
  const int n_dwo = countInstances(test_tf, "dworps", text_separators, text_delimiters);
  const int n_glp = countInstances(test_tf, "glipment", text_separators, text_delimiters);
  const int n_gbr = countInstances(test_tf, "gleeped bretly", text_separators, text_delimiters);
  const std::vector<std::string> dual_words = { "gleeped", "bretly" };
  const int n_dwd = countInstances(test_tf, dual_words, text_separators, text_delimiters);
  const std::vector<std::string> comma_words = { "gleeped", ",", "and" };
  const int n_cma = countInstances(test_tf, comma_words, text_separators, text_delimiters);
  check(n_two, RelationalOperator::EQUAL, 1, "The word 'Two' was not found the expected number of "
        "times in a short paragraph.");
  check(n_dwo, RelationalOperator::EQUAL, 1, "The word 'dworps' was not found the expected number "
        "of times in a short paragraph.");
  check(n_glp, RelationalOperator::EQUAL, 1, "The word 'glipment' was not found the expected "
        "number of times in a short paragraph.");
  check(n_gbr, RelationalOperator::EQUAL, 1, "The words 'gleeped bretly' were not found the "
        "expected number of times in a short paragraph.");
  check(n_dwd, RelationalOperator::EQUAL, 1, "The sequence of words 'gleeped bretly' was not "
        "found the expected number of times in a short paragraph.");
  check(n_cma, RelationalOperator::EQUAL, 1, "The sequence of words 'gleeped, and' was not "
        "found the expected number of times in a short paragraph.");
  const std::string poem_tfstr("Mary (had a) little lamb,\n  Its fleece was white as snow\nAnd "
                               "everywhere that Mary went\nThe lamb was sure to go.");
  const TextFile poem_tf(poem_tfstr, TextOrigin::RAM);
  const int n_mary = countInstances(poem_tf, "Mary", text_separators, text_delimiters);
  const int n_have = countInstances(poem_tf, { "Mary", "had", "a", "little" }, text_separators,
                                    text_delimiters);
  const int n_lnbr = countInstances(poem_tf, { "lamb", ",", "Its", "fleece" }, text_separators,
                                    text_delimiters);
  check(n_mary, RelationalOperator::EQUAL, 2, "The number of instances of 'Mary' counted in a "
        "short poem is incorrect.");
  check(n_have, RelationalOperator::EQUAL, 1, "The number of instances of 'Mary had a little' "
        "counted in a short poem is incorrect.");
  check(n_lnbr, RelationalOperator::EQUAL, 1, "The   number of instances of 'lamb, Its fleece' "
        "counted in a short poem is incorrect.");
  const std::string tech_tfstr("Some text can be hidden behind comments\n"
                               "To hide it from being accessible to grep.\n"
                               "Do it // like this.  Some text here is not greppable.\n"
                               "See?  Some text is not greppable.\n");
  const TextFile tech_tf(tech_tfstr, TextOrigin::RAM);
  const std::vector<std::string> some_text_vec = { "Some", "text" };
  const int n_smtx = countInstances(tech_tf, some_text_vec, text_separators, text_delimiters);
  const int n_smtx_comm = countInstances(tech_tf, some_text_vec, text_separators,
                                         text_delimiters, {}, { TextGuard("//") });
  check(n_smtx, RelationalOperator::EQUAL, 3, "The words 'Some text' do not appear the expected "
        "number of times in a short paragraph with technical markup.");
  check(n_smtx_comm, RelationalOperator::EQUAL, 2, "The words 'Some text' are not protected "
        "from analysis by comment symbols as they should be.");
  
  // Print results
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
