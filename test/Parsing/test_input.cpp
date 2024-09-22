#include "copyright.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Namelists/input.h"
#include "../../src/Namelists/user_settings.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/parsing_enumerators.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::constants::ExceptionResponse;
using stormm::diskutil::osSeparator;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::errors::rtWarn;
using stormm::parse::separateText;
using stormm::parse::TextOrigin;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using namespace stormm::namelist;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// This function will wrap the creation of the &lisa namelist.  In practice, most namelists should
// be created inside a separate function like this which fills out the keywords and then
// immediately reads the data starting from a point in a file as indicated in the arguments.  Even
// additional copies of the same namelist can, in principle, be constructed and then adjust with
// user-specified inputs in this way.  Declaring the NamelistEmulator results of such namelist
// reader functions to be const variables when they arrive back in the calling function implicitly
// prevents accidental editing of user input once it has been read, and also allows a program to
// read multiple instances of a given namelist.
//
// Arguments:
//   tf:          The text file containing user input, committed to RAM
//   start_line:  First line of the text file at whcih to begin searching for relevant data
//-------------------------------------------------------------------------------------------------
NamelistEmulator lisaInput(const TextFile &tf, int *start_line) {

  NamelistEmulator t_nml("lisa", CaseSensitivity::YES, ExceptionResponse::DIE, "An amicable and "
                        "level-headed perpetual second grader.");
  t_nml.addKeyword(NamelistElement("FoodPrefs", NamelistType::STRING, "vegetarian"));
  t_nml.addKeyword(NamelistElement("Instrument", NamelistType::STRING, "saxophone"));
  t_nml.addKeyword(NamelistElement("Grade", NamelistType::INTEGER, "2"));
  t_nml.addKeyword(NamelistElement("Age", NamelistType::INTEGER, "7"));
  t_nml.addKeyword(NamelistElement("GenesisTub", { "Source", "Product", "Population", "Height" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::INTEGER, NamelistType::REAL },
                                   { "baby molar", "Lutherans", "619", "1.1" },
                                   DefaultIsObligatory::YES, InputRepeats::YES,
                                   "Lisa's terrarium becomes the site of an amazing abiogenesis "
                                   "experiment.",
                                   { "The object Lisa places into the terrarium that brings about "
                                     "the miracle of life", "The sort of life she creates", "The "
                                     "number of living organisms arising from the experiment",
                                     "The average height of creatures produced in the "
                                     "terrarium" }));

  // This function will automatically read as far as the end and not wrap back ot the start
  // of the file.  Other reader functions may include additional arguments to control the portion
  // of the file that gets searched.
  *start_line = readNamelist(tf, &t_nml, *start_line, WrapTextSearch::NO, tf.getLineCount());

  return t_nml;
}

//-------------------------------------------------------------------------------------------------
// Create a generic namelist with simple keywords corresponding to each type.
//-------------------------------------------------------------------------------------------------
NamelistEmulator genericNamelist() {
  NamelistEmulator t_nml("dtopic", CaseSensitivity::YES, ExceptionResponse::DIE, "A minimal "
                         "testing namelist");
  t_nml.addKeyword("kw_real", NamelistType::REAL);
  t_nml.addKeyword("kw_string", NamelistType::STRING);
  t_nml.addKeyword("kw_integer", NamelistType::INTEGER);
  return t_nml;
}
//-------------------------------------------------------------------------------------------------
// The UserSettings class
//-------------------------------------------------------------------------------------------------
UserSettings spaceSettings(const std::string &remdynamics) {

  // Testing the UserSettings Class
  const char* argv[] = {"", "-i", remdynamics.c_str()};
  UserSettings ui(3, argv, AppName::DYNAMICS);
  check(ui.getRemdPresence(), "The UserInput object created from " + remdynamics + " did not "
        "contain a &remd namelist.");
  check(ui.getDynamicsPresence(), "The UserInput object created from " + remdynamics + " did not "
        "contain a &dynamics namelist.");
  if (ui.getRemdPresence() && ui.getDynamicsPresence()) {
    return ui;
  }
  else {
    rtErr("The input file was not read correctly.", "spaceSettings");
  }
  __builtin_unreachable();
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
  section("Namelist creation and default initialization");

  // Section 2
  section("TextFile parsing into input namelists");

  // Section 3
  section("Extracting information from namelists");

  // Section 4
  section("Testing the UserSettings Class for multiple namelists");

  // Try making a NamelistEmulator with a number of NamelistElement objects in it
  section(1);
  NamelistEmulator nml_a("monteburns", CaseSensitivity::YES);
  nml_a.addHelp("The world's worst boss, or perhaps just a typical aristocratic tycoon with too "
                "much inherited wealth and not enough worldly experience.  Either way, his "
                "comical inability to understand the common man is emblematic of the 0.001%.");
  nml_a.addKeyword(NamelistElement("xdim", NamelistType::REAL, "4.0"));
  nml_a.addKeyword(NamelistElement("ydim", NamelistType::REAL, "4.0"));
  nml_a.addKeyword(NamelistElement("zdim", NamelistType::REAL, "8.0"));
  nml_a.addHelp("xdim", "The width of Monte Burns's desk.");
  nml_a.addHelp("ydim", "The vertical height of Marge Simpson's painting of Mr. Burns.");
  nml_a.addHelp("zdim", "The height of Monte Burns's recycling factory, in feet.");
  nml_a.addKeyword(NamelistElement("nbox", NamelistType::INTEGER, "2"));
  nml_a.addKeyword(NamelistElement("boxcolor", NamelistType::STRING, "Red",
                                     DefaultIsObligatory::YES, InputRepeats::YES));
  nml_a.addKeyword(NamelistElement("boxmodel", NamelistType::STRING, "Studebaker",
                                     DefaultIsObligatory::NO, InputRepeats::NO, "The type of "
                                     "vintage automobile that Homer and Monte Burns were driving "
                                     "when Chief Wiggums saw them."));
  check(nml_a.getKeywordCount(), RelationalOperator::EQUAL, 6, "Namelist " + nml_a.getTitle() +
        " reports an incorrect number of keywords.");
  check(nml_a.getRealValue("ydim"), RelationalOperator::EQUAL, Approx(4.0), "Namelist " +
        nml_a.getTitle() + " carries the wrong default value for ydim");
  check(nml_a.getIntValue("nbox"), RelationalOperator::EQUAL, 2, "Namelist " + nml_a.getTitle() +
        " records the wrong default INTEGER value for nbox.");
  check(nml_a.getStringValue("boxmodel") == "Studebaker", "Namelist " + nml_a.getTitle() +
        " records the wrong boxmodel STRING value.");
  CHECK_THROWS(const int test_int = nml_a.getIntValue("xdim"), "An INTEGER value was extracted "
               "from a REAL namelist keyword.");
  CHECK_THROWS(const std::string test_str = nml_a.getStringValue("zdim"), "A STRING value was "
               "extracted from a REAL namelist keyword.");
  NamelistEmulator nml_b("wiggums", CaseSensitivity::NO, ExceptionResponse::DIE, "A stereotypical "
                         "cop, so much so that one wonders if his character will go the way of "
                         "Apu.");
  nml_b.addKeyword(NamelistElement("rank", NamelistType::STRING, "Chief"));
  nml_b.addKeyword(NamelistElement("reply", NamelistType::STRING, "More of a burgundy..."));
  nml_b.addKeyword(NamelistElement("weight", NamelistType::REAL, "220.05"));
  nml_b.addKeyword({ NamelistElement("parking", NamelistType::INTEGER, "45"),
                     NamelistElement("badge", NamelistType::REAL, "54.57") });
  nml_b.addKeyword(NamelistElement("deputies", NamelistType::INTEGER));
  nml_b.addKeyword(NamelistElement("show", { "producer", "animator", "vocalist", "seasons" },
                                     { NamelistType::STRING, NamelistType::STRING,
                                       NamelistType::STRING, NamelistType::INTEGER },
                                     { "Matt Groenig", "David Silverman", "Hank Azaria", "33" }));
  nml_b.addHelp("show", "animator", "The chief animation specialist.");
  CHECK_THROWS(nml_b.addKeyword(NamelistElement("show", NamelistType::STRUCT, "randomly")),
               "A STRUCT namelist element was added with no member variables, just a string "
               "default value.");
  CHECK_THROWS(nml_b.addHelp("show", "janitor", "No such person on the payroll."), "Namelist "
               "addHelp() attempts to add a help message to a keyword that does not exist.");
  nml_b.addKeyword(NamelistElement("broadcast", { "seasons", "network", "neilsen" },
                                     { NamelistType::INTEGER, NamelistType::STRING,
                                       NamelistType::REAL }, { "33", "Fox", "4.7" }));
  CHECK_THROWS(nml_b.addKeyword(NamelistElement("height", NamelistType::REAL, "5' 11\"")),
               "A REAL keyword was added to a namelist with a STRING value default.");
  CHECK_THROWS(nml_b.addKeyword(NamelistElement("show", NamelistType::REAL, "706.24")),
               "A keyword was added to a namelist with the same name as an existing keyword.");

  // Try reading a namelist from an example file
  section(2);
  const std::string namelist_file = oe.getStormmSourcePath() + osSeparator() + "test" +
    osSeparator() + "Parsing" + osSeparator() + "simpsons.txt";
  const bool file_exists = (getDrivePathType(namelist_file) == DrivePathType::FILE);
  const TestPriority file_check = (file_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  const TextFile tf = (file_exists) ? TextFile(namelist_file) : TextFile();
  int start_line = 0;
  if (file_exists) {
    start_line = readNamelist(tf, &nml_a, start_line, WrapTextSearch::NO, tf.getLineCount());  
    start_line = readNamelist(tf, &nml_b, start_line, WrapTextSearch::NO, tf.getLineCount());  
  }
  else {
    rtWarn("The file " + namelist_file + " was not found and is critical to many subsequent "
           "tests.  Make sure that the $STORMM_SOURCE environment variable is set properly, to "
           "the source tree where src/ and test/ subdirectories can be found.  Subsequent tests "
           "that rely on this file will be skipped.", "test_input");
  }
  std::vector<std::string> help_words_b = separateText(nml_b.getHelp());
  check(help_words_b[1], RelationalOperator::EQUAL, "stereotypical",
        "User documentation produced with the getHelp() function does not meet expectations.",
        file_check);
  check(nml_a.getHelp("xdim"), RelationalOperator::EQUAL, "The width of Monte Burns's desk.",
        "User documentation produced with the getHelp() function for a specific keyword does not "
        "meet expectations.");
  check(nml_b.getHelp("show", "animator"), RelationalOperator::EQUAL,
        "The chief animation specialist.", "User documentation produced with the getHelp() "
        "function for a keyword / subkey pair does not meet expectations.");

  // A function encapsulates the creation and parsing of the &lisa namelist (see above)
  start_line = 34;
  NamelistEmulator nml_c = lisaInput(tf, &start_line);

  // Also try reading other &lisa namelists from the same file that have deliberate errors
  start_line = 0;
  CHECK_THROWS_SOFT(NamelistEmulator nml_c1 = lisaInput(tf, &start_line), "Reading namelist &lisa "
                    "did not throw an exception as expected when receiving an unknown keyword.",
                    file_check);
  start_line = 27;
  CHECK_THROWS_SOFT(NamelistEmulator nml_c2 = lisaInput(tf, &start_line), "Reading namelist &lisa "
                    "did not throw an exception as expected when reading a real valued number in "
                    "place of an integer.", file_check);
  start_line = 48;
  CHECK_THROWS_SOFT(NamelistEmulator nml_c3 = lisaInput(tf, &start_line), "Namelist &lisa was "
                    "read with two specifications of the same sub-key inside data associated with "
                    "a STRUCT keyword.", file_check);

  // Additional checks on what can be added to a namelist.  The object nml_c was not declared
  // const so that these modifications could be attempted.  In practice, it's possible, but
  // probably not the best way to code.
  section(1);
  CHECK_THROWS(nml_c.addKeyword(NamelistElement("Poindextrose", NamelistType::REAL, "15.4",
                                                  DefaultIsObligatory::YES, InputRepeats::NO)),
               "A namelist element was added with an obligatory default and a non-extensible "
               "value list."); 
  CHECK_THROWS(nml_c.addKeyword(NamelistElement("XX", NamelistType::REAL, "this_is_no-number")),
               "A NamelistElement for a REAL keyword with a non-numerical default value was "
               "successfully added to namelist \"" + nml_c.getTitle() + "\".");
  CHECK_THROWS(nml_c.addKeyword(NamelistElement("YY", NamelistType::INTEGER, "1093bt")),
               "A NamelistElement for an INTEGER keyword with a non-numerical default value was "
               "successfully added to namelist \"" + nml_c.getTitle() + "\".");
  CHECK_THROWS(nml_c.getHelp("XX"), "A badly formed namelist keyword was added to namelist \"" +
               nml_c.getTitle() + "\".");
  CHECK_THROWS_SOFT(int iz = readNamelist(tf, &nml_c, 48, WrapTextSearch::NO, tf.getLineCount()),
                    "A namelist with existing data was put through another input reading cycle.",
                    file_check);

  // Interpret the namelist data
  section(3);
  check(nml_a.getIntValue("nbox"), RelationalOperator::EQUAL, 57, "Namelist " + nml_a.getTitle() +
        " produces the wrong value for its nbox keyword.", file_check);
  const std::vector<std::string> car_colors = nml_a.getAllStringValues("boxcolor");
  check(car_colors.size(), RelationalOperator::EQUAL, 4, "Retrieval of a multi-entry STRING "
        "keyword's data failed.", file_check);
  check(car_colors[3], RelationalOperator::EQUAL, "And a 'third' value", "Entries of a STRING "
        "keyword's data were not retrieved as expected.", file_check);
  const std::vector<double> gt_heights       = nml_c.getAllRealValues("GenesisTub", "Height");
  const std::vector<int> gt_pops             = nml_c.getAllIntValues("GenesisTub", "Population");
  const std::vector<std::string> gt_products = nml_c.getAllStringValues("GenesisTub", "Product");
  const std::vector<std::string> gt_sources  = nml_c.getAllStringValues("GenesisTub", "Source");
  CHECK_THROWS_SOFT(const std::vector<std::string> gtb =
                    nml_c.getAllStringValues("GenesisTub", "Height"), "A request for STRING "
                    "values from a namelist STRUCT object was allowed to access the data as if it "
                    "were REALs.", file_check);
  CHECK_THROWS_SOFT(const std::vector<int> gtb2 = nml_c.getAllIntValues("GenesisTub", "blah"),
                    "A request for non-existent INTEGER values from a namelist STRUCT object was "
                    "allowed to go through.", file_check);
  check(gt_heights.size(), RelationalOperator::EQUAL, 4, "A series of real values was not "
        "properly returned for a sub-key of a namelist STRUCT keword.", file_check);
  check(gt_heights.size(), RelationalOperator::EQUAL, 4, "A series of real values was not "
        "properly returned for a sub-key of a namelist STRUCT keword.", file_check);
  check(gt_pops[0], RelationalOperator::EQUAL, 619, "The default value for an INTEGER sub-key "
        "within a STRUCT keyword with an obligatory default entry was not properly maintained.");
  const double gth2 = (file_exists) ? gt_heights[2] : 0.0;
  check(gth2, RelationalOperator::EQUAL, 1.1, "The default value for an INTEGER sub-key "
        "within a STRUCT keyword was not properly applied to an entry lacking a specification "
        "for this value.", file_check);
  const std::string gts3 = (file_exists) ? gt_sources[3] : "NONE";
  check(gts3, RelationalOperator::EQUAL, "starfish", "A STRING value for a sub-key of a repeated "
        "STRUCT keyword did not carry through correctly.", file_check);
  const double gth3 = (file_exists) ? gt_heights[3] : 0.0;
  check(gth3, RelationalOperator::EQUAL, 4.20, "A REAL value for a sub-key of a repeated STRUCT "
        "keyword did not carry through correctly.", file_check);

  // Check that missing input variables do no change existing values external to the namelist.
  const std::string nml_d_input_a("&dtopic\n&end\n");
  const TextFile nml_d_tf_a(nml_d_input_a, TextOrigin::RAM);
  NamelistEmulator nml_d = genericNamelist();
  start_line = 0;
  start_line = readNamelist(nml_d_tf_a, &nml_d, start_line, WrapTextSearch::NO,
                            nml_d_tf_a.getLineCount());
  double my_real = -0.6;
  std::string my_string("not_your_average_string");
  int my_int = 9;
  nml_d.assignVariable(&my_real, "kw_real");
  nml_d.assignVariable(&my_string, "kw_string");
  nml_d.assignVariable(&my_int, "kw_integer");
  check(my_real, RelationalOperator::EQUAL, -0.6, "A real-valued variable was altered by reading "
        "from an unspecified namelist keyword.");
  check(my_string, RelationalOperator::EQUAL, "not_your_average_string", "A string variable was "
        "altered by reading from an unspecified namelist keyword.");
  check(my_int, RelationalOperator::EQUAL, 9, "An integer variable was altered by reading from an "
        "unspecified namelist keyword.");

  // Check that the NamelistEmulator object is able to assign values as required.
  const std::string nml_d_input_b("&dtopic\n  kw_real 5.7, kw_string \"some_string\" "
                                  "kw_integer 8,\n&end\n");
  const TextFile nml_d_tf_b(nml_d_input_b, TextOrigin::RAM);
  nml_d = genericNamelist();
  start_line = 0;
  start_line = readNamelist(nml_d_tf_b, &nml_d, start_line, WrapTextSearch::NO,
                            nml_d_tf_b.getLineCount());
  nml_d.assignVariable(&my_real, "kw_real");
  nml_d.assignVariable(&my_string, "kw_string");
  nml_d.assignVariable(&my_int, "kw_integer");
  std::vector<int> trip_int = { -1, -1, -1 };
  nml_d.assignVariable(&trip_int[0], &trip_int[1], &trip_int[2], "kw_integer");
  check(my_real, RelationalOperator::EQUAL, 5.7, "A real-valued variable was not assigned "
        "correctly.");
  check(my_string, RelationalOperator::EQUAL, "some_string", "A string variable was not assigned "
        "correctly.");
  check(my_int, RelationalOperator::EQUAL, 8, "An integer variable was not assigned correctly.");
  check(trip_int, RelationalOperator::EQUAL, std::vector<int>(3, 8), "Integers were not assigned "
        "as expected in triplicate by a generic &namelist object.");
  
  // Testing the UserSettings Class
  section(4);
  const std::string remdynamics = oe.getStormmSourcePath() + osSeparator() + "test" +
    osSeparator() + "Parsing" + osSeparator() + "remdynamics.txt";
  UserSettings ui = spaceSettings(remdynamics);
  check(ui.getDynamicsPresence(), RelationalOperator::EQUAL, true, "Dynamics Namelist was not"
        " detected successfully.");
  check(ui.getRemdPresence(), RelationalOperator::EQUAL, true, "REMD Namelist was not detected "
        "successfully.");
  check(ui.getFilesPresence(), RelationalOperator::EQUAL, false, "Files Namelist was detected, but"
        " not actually present.");
  DynamicsControls dyna_a = ui.getDynamicsNamelistInfo();
  RemdControls remd_a = ui.getRemdNamelistInfo();
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

  // Summary evaluation
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
