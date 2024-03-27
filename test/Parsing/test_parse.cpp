#include "copyright.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Parsing/ascii_numbers.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Topology/amber_prmtop_util.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::data_types::ullint;
#ifndef STORMM_USE_HPC
using stormm::data_types::uchar4;
#endif
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::errors::rtWarn;
using stormm::stmath::sum;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using stormm::review::protectText;
using stormm::topology::amberPrmtopData;

using namespace stormm::parse;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1 (this illustrates "forward declaration" of test sections)
  section("TextFile manipulation and abstraction");

  // Section 2
  section("TextGuard operation, putting bounds on quoted or commented text within a character "
          "stream");

  // Section 3
  section("Test number format verification, extraction, and display");

  // Section 4
  section("Test character conversion and ASCII text table correspondence");

  // Section 5
  section("Test formatted column data reading");

  // Section 6
  section("Test string comparisons with wildcards");

  // Section 7
  section("Test string manipulation");
  
  // Test the TextFile object and its features
  section(1);
  const std::string something_path = oe.getStormmSourcePath() + osSeparator() + "test" +
                                     osSeparator() + "Parsing" + osSeparator() + "something.txt";
  const bool sfile_exists = (getDrivePathType(something_path) == DrivePathType::FILE);
  const TextFile tf = (sfile_exists) ? TextFile(something_path) : TextFile();
  if (sfile_exists == false) {
    rtWarn("Text file " + something_path + " was not found.  Make sure that the $STORMM_SOURCE "
           "environment variable is set to the root of the source tree, where src/ and test/ "
           "subdirectories can be found.  Many subsequent tests will be skipped until this file "
           "is available.", "test_parse");
  }
  const TestPriority sfile_check = (sfile_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  check(tf.getFileName(), RelationalOperator::EQUAL, oe.getStormmSourcePath() + osSeparator() +
        "test" + osSeparator() + "Parsing" + osSeparator() + "something.txt", "Text file name was "
        "not properly identified.", sfile_check);
  check(tf.getLineCount(), RelationalOperator::EQUAL, 40, "TextFile object records an incorrect "
        "number of lines.", sfile_check);
  const int llim5 = (sfile_exists) ? tf.getLineLimits(5) : 0;
  check(llim5, RelationalOperator::EQUAL, 150, "TextFile line limits are incorrect.", sfile_check);
  const char gt54 = (sfile_exists) ? tf.getChar(54) : ' ';
  const char gt56 = (sfile_exists) ? tf.getChar(56) : ' ';
  check(gt54 == 'e' && gt56 == 't', "TextFile contains incorrect text.", sfile_check);
  const TextFileReader tfr = tf.data();
  check(tfr.line_count, RelationalOperator::EQUAL, tf.getLineCount(), "TextFileReader nested "
        "struct does not produce a line count in agreement with the original struct.",
        sfile_check);
  const bool comp15 = (sfile_exists) ? (tfr.line_limits[15] == tf.getLineLimits(15)) : false;
  check(comp15, "TextFileReader does not report line limits in agreement with the original "
        "struct.", sfile_check);
  if (oe.getTemporaryDirectoryAccess() == false) {
    rtWarn("Write access to the temporary directory, " + oe.getTemporaryDirectoryPath() + ", is "
           "unavailable.  Subsequent tests for writing a TextFile object will be skipped.",
           "test_parse");
  }
  const TestPriority do_write_test = oe.getTemporaryDirectoryAccess() ? TestPriority::CRITICAL :
                                                                        TestPriority::ABORT;
  const std::string text_to_obj("There is\n  some text here.\n\"What about a quote?\""
                                "\n\n\nblah.\n");
  const TextFile aether(text_to_obj, TextOrigin::RAM);
  check(aether.getLineCount(), RelationalOperator::EQUAL, 6, "The number of lines in a TextFile "
        "object built internally from a string is incorrect.");
  const std::vector<int> aether_lims = { aether.getLineLimits(0), aether.getLineLimits(1),
                                         aether.getLineLimits(2), aether.getLineLimits(3),
                                         aether.getLineLimits(4), aether.getLineLimits(5),
                                         aether.getLineLimits(6) };
  const std::vector<int> aether_lims_ans = { 0, 8, 25, 46, 46, 46, 51 };
  check(aether_lims, RelationalOperator::EQUAL, aether_lims_ans, "The line demarcations of a "
        "TextFile object built internally from a string are not correct.");
  CHECK_THROWS(aether.write(), "An attempt was made to write a TextFile object with no inherent "
               "name or surrogate name.");
  const std::string aether_fi = oe.getTemporaryDirectoryPath() + osSeparator() + "aether.txt";
  TextFile aether_read;
  std::vector<int> aether_read_lims(7, 0);
  if (oe.getTemporaryDirectoryAccess()) {
    aether.write(aether_fi, PrintSituation::OPEN_NEW);
    oe.logFileCreated(aether_fi);
    aether_read = TextFile(aether_fi);
    for (int i = 0; i < 7; i++) {
      aether_read_lims[i] = aether_read.getLineLimits(i);
    }
  }
  check(aether_read_lims, RelationalOperator::EQUAL, aether_lims_ans, "The line demarcations of a "
        "TextFile object built internally, written to disk, and read back are not correct.",
        do_write_test);
  
  // Test the comment and quotation masking
  section(2);
  const std::vector<TextGuard> comments = { TextGuard("//"),
                                            TextGuard("/*", "*/", LineSpan::MULTIPLE) };
  const std::vector<TextGuard> quotations = { TextGuard("'", "'"),
                                              TextGuard("\"", "\"", LineSpan::MULTIPLE) };
  const std::vector<bool> commented_text = markGuardedText(tfr, comments, quotations);
  const std::vector<bool> quoted_text = markGuardedText(tfr, quotations, comments);
  check(tfr.line_limits[tfr.line_count], RelationalOperator::EQUAL, commented_text.size(),
        "Size of comments mask vector does not match the reference text file.", sfile_check);
  check(tfr.line_limits[tfr.line_count], RelationalOperator::EQUAL, commented_text.size(),
        "Size of quotations mask vector does not match the reference text file.", sfile_check);
  check(sum<int>(commented_text), RelationalOperator::EQUAL, 683, "An incorrect number of "
        "characters were marked as being within comments in \"" + tf.getFileName() + "\".",
        sfile_check);
  check(sum<int>(quoted_text), RelationalOperator::EQUAL, 269, "An incorrect number of "
        "characters were marked as being within quotations in \"" + tf.getFileName() + "\".",
        sfile_check);

  // Test scope evaluation
  const std::string exprsn("90 + (7 * 8) + ((t - v) * g + 8) / x");
  const std::vector<TextGuard> scope_characters = { TextGuard("(", ")"), TextGuard("{", "}"),
                                                    TextGuard("[", "]") };
  const std::vector<int> ex_ops = resolveScopes(exprsn, scope_characters);
  const std::vector<int> expected_scp = { 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 2,
                                          2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0 };
  check(ex_ops, RelationalOperator::EQUAL, expected_scp, "Scope resolution did not produce the "
        "correct result.", sfile_check);

  // Test some of the exceptions thrown by text guarding
  const std::vector<TextGuard> braces = { TextGuard("{", "}", LineSpan::MULTIPLE) };
  const std::vector<TextGuard> brackets = { TextGuard("[", "]", LineSpan::SINGLE) };
  const std::vector<bool> braced_text = markGuardedText(tfr, braces, comments + quotations);
  check(sum<int>(braced_text), RelationalOperator::EQUAL, 48, "An incorrect number of "
        "characters were marked as being within braces in \"" + tf.getFileName() + "\".",
        sfile_check);
  CHECK_THROWS(const TextGuard bad_txg_a("", "|"), "Initializing a TextGuard object without a "
               "left-hand, initiating delimiter did not throw an error.");
  CHECK_THROWS(const TextGuard bad_txg_b("|", "", LineSpan::MULTIPLE), "Initializing a TextGuard "
               "object with multiple line span capability and no right-hand, terminating "
               "delimiter did not throw an error.");
  CHECK_THROWS_SOFT(const std::vector<bool> bracketed_text =
                    markGuardedText(tfr, brackets, comments + quotations), "Failure to terminate "
                    "a single-line text guard sequence (\"[\", \"]\") in file \"" +
                    tf.getFileName() + "\" did not throw an error.", sfile_check);

  // Test number interpretation
  section(3);
  double dbl_real = 0.000201;
  check(realToString(dbl_real, 4), RelationalOperator::EQUAL, "0.0002", "realToString() reports "
        "an incorrect representation.");
  check(realToString(dbl_real, 6), RelationalOperator::EQUAL, "0.000201", "realToString() reports "
        "an incorrect representation.");
  check(realToString(dbl_real, 10, 6), RelationalOperator::EQUAL, "  0.000201", "realToString() "
        "reports an incorrect representation.");
  check(realToString(dbl_real, 10, 6, NumberFormat::STANDARD_REAL,
                     NumberPrintStyle::LEADING_ZEROS), RelationalOperator::EQUAL, "000.000201",
        "realToString() reports an incorrect representation.");
  check(realDecimalPlaces(dbl_real, 4), RelationalOperator::EQUAL, 4, "Incorrect number of "
        "decimal places reported by realDecimalPlaces for " + realToString(dbl_real, 6) + ".");
  CHECK_THROWS(const std::string test_dbl_str_a =
               realToString(dbl_real, 15, 4, NumberFormat::INTEGER), "An incorrect numerical "
               "format was passed to realToString() without being trapped.");
  CHECK_THROWS(const std::string test_dbl_str_b = realToString(dbl_real, 65), "An obscenely long "
               "real number format was permitted.");
  CHECK_THROWS(const std::string test_dbl_str_c = realToString(dbl_real, 16, 78), "An obscene "
               "number of decimal places was permitted in a real number.");
  CHECK_THROWS(const std::string test_dbl_str_d = realToString(dbl_real, 16, 15), "The number of "
               "decimal places was permitted to exceed the overall format in a real number.");
  CHECK_THROWS(const std::string test_dbl_str_e =
               realToString(dbl_real, 16, 12, NumberFormat::SCIENTIFIC), "The number of "
               "decimal places was permitted to exceed the overall (scientific notation) format "
               "in a real number.");
  dbl_real = 686.4;
  check(realDecimalPlaces(dbl_real), RelationalOperator::EQUAL, 1, "Incorrect number of "
        "decimal places reported by realDecimalPlaces() for " + realToString(dbl_real, 5) + ".");
  check(realToString(dbl_real, 12, 5, NumberFormat::SCIENTIFIC), RelationalOperator::EQUAL,
        " 6.86400e+02", "realToString() reports an incorrect representation.");
  check(realToString(dbl_real, 12, 5, NumberFormat::SCIENTIFIC, NumberPrintStyle::LEADING_ZEROS),
        RelationalOperator::EQUAL, "06.86400e+02", "realToString() reports an incorrect "
        "representation.");
  check(realToString(-dbl_real, 12, 5, NumberFormat::SCIENTIFIC, NumberPrintStyle::LEADING_ZEROS),
        RelationalOperator::EQUAL, "-6.86400e+02", "realToString() reports an incorrect "
        "representation.");
  check(realToString(-2.0 * dbl_real, 14, 5, NumberFormat::SCIENTIFIC,
                     NumberPrintStyle::LEADING_ZEROS), RelationalOperator::EQUAL, "-001.37280e+03",
        "realToString() reports an incorrect representation.");
  check(minimalRealFormat(2.0e-9, 1.0e-6), RelationalOperator::EQUAL, "2.0e-9", "A real number "
        "was not represented its most compact form.");
  check(minimalRealFormat(2.1245e-01, 1.0e-6), RelationalOperator::EQUAL, "0.21245", "A real "
        "number was not represented its most compact form.");
  check(minimalRealFormat(-7.06e+08, 1.0e-6), RelationalOperator::EQUAL, "-7.06e8", "A real "
        "number was not represented its most compact form.");
  check(minimalRealFormat(-1.071593e+08, 1.0e-6), RelationalOperator::EQUAL, "-107159300",
        "A real number was not represented its most compact form.");

  // Test character conversion
  section(4);
  check(uppercase('a') == 'A', "Conversion of 'a' to 'A' failed.");
  check(uppercase('g') == 'G', "Conversion of 'g' to 'G' failed.");
  check(uppercase('z') == 'Z', "Conversion of 'z' to 'Z' failed.");
  check(uppercase('#') == '#', "Uppercase conversion of '#' changed the character.");
  check(lowercase('A') == 'a', "Conversion of 'A' to 'a' failed.");
  check(lowercase('G') == 'g', "Conversion of 'G' to 'g' failed.");
  check(lowercase('Z') == 'z', "Conversion of 'Z' to 'a' failed.");
  check(lowercase('<') == '<', "Lowercase conversion of '<' changed the character.");
  check(lowercase(std::string("@!&%AJOfjoTP782bi__<")), RelationalOperator::EQUAL,
        "@!&%ajofjotp782bi__<", "Lowercase conversion of a mixed C++ string failed.");
  check(uppercase(std::string("ggz_a@!&%AJOfjoTP782bi!*x")), RelationalOperator::EQUAL,
        "GGZ_A@!&%AJOFJOTP782BI!*X", "Uppercase conversion of a mixed C++ string failed.");
  check(lowercase("@!&%AJOfjoTP782bi__<"), RelationalOperator::EQUAL,
        "@!&%ajofjotp782bi__<", "Lowercase conversion of a mixed, const C-style string failed.");
  check(uppercase("ttyl78_0bama"), RelationalOperator::EQUAL, "TTYL78_0BAMA",
        "Lowercase conversion of a mixed, const C-style string failed.");
  char test_str[32];
  snprintf(test_str, 32, "upaiwof89xxv");
  uppercase(test_str);
  check(std::string(test_str), RelationalOperator::EQUAL, "UPAIWOF89XXV",
        "Uppercase conversion of a mixed C-style string failed.");
  snprintf(test_str, 32, "upaiwof89xxv");
  uppercase(test_str, 4);
  check(std::string(test_str), RelationalOperator::EQUAL, "UPAIwof89xxv",
        "Uppercase partial conversion of a mixed C-style string failed.");
  snprintf(test_str, 32, "VDNnL@)*cCB");
  lowercase(test_str);
  check(std::string(test_str), RelationalOperator::EQUAL, "vdnnl@)*ccb",
        "Lowercase conversion of a mixed C-style string failed.");
  snprintf(test_str, 32, "VDNnL@)*cCB");
  lowercase(test_str, 4);
  check(std::string(test_str), RelationalOperator::EQUAL, "vdnnL@)*cCB",
        "Lowercase partial conversion of a mixed C-style string failed.");
  const uchar4 rgb_code_i = { 255, 255, 255, 0 };
  const std::string str_rgb_i = rgbHexCode(rgb_code_i);
  check(str_rgb_i, RelationalOperator::EQUAL, "ffffff", "Conversion of an eight-bit color triplet "
        "does not meet expectations.");
  const uchar4 rgb_code_ii = {   5, 68, 170, 0 };
  const std::string str_rgb_ii = rgbHexCode(rgb_code_ii);
  check(str_rgb_ii, RelationalOperator::EQUAL, "0544aa", "Conversion of an eight-bit color "
        "triplet does not meet expectations.");
  
  // Test formatted column data reading: starting with integers
  section(5);
  const std::string spsh_path = oe.getStormmSourcePath() + osSeparator() + "test" + osSeparator() +
                                "Parsing" + osSeparator() + "formatted_numbers.txt";
  const bool pfile_exists = (getDrivePathType(spsh_path) == DrivePathType::FILE);
  if (pfile_exists == false) {
    rtWarn("Data file " + spsh_path + " was not found.  Make sure that the $STORMM_SOURCE "
           "environment variable is set to the root of the source tree, where src/ and test/ "
           "subdirectories can be found.  Many subsequent tests will be skipped until this file "
           "is available.", "test_parse");
  }
  const TextFile spreadsheet = (pfile_exists) ? TextFile(spsh_path) : TextFile();
  const TestPriority pfile_check = (pfile_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  PolyNumeric one_pn;
  one_pn.i = 0;
  const std::vector<PolyNumeric> blank_pn = std::vector<PolyNumeric>(1, one_pn);
  const std::vector<PolyNumeric> apn_ints = (pfile_exists) ?
    amberPrmtopData(spreadsheet, 1, NumberFormat::INTEGER, 4, 4, 14) : blank_pn;
  const std::vector<PolyNumeric> bpn_ints = (pfile_exists) ?
    amberPrmtopData(spreadsheet, 1, NumberFormat::INTEGER, 4, 4, 12, 14) : blank_pn;
  const std::vector<int> astd_ints = intFromPolyNumeric(apn_ints);
  const std::vector<int> bstd_ints = intFromPolyNumeric(bpn_ints);
  check(astd_ints == Approx(bstd_ints), "Integer vectors read from " + tf.getFileName() +
        " should be identical, irrespective of the soft limit on reading by amberPrmtopData().",
        pfile_check);
  check(astd_ints[5], RelationalOperator::EQUAL, 90, "Formatter integers read by "
        "amberPrmtopData() from " + tf.getFileName() + " do not appear to be correct.",
        pfile_check);
  const std::vector<PolyNumeric> pn_dbls = (pfile_exists) ?
    amberPrmtopData(spreadsheet, 8, NumberFormat::STANDARD_REAL, 5, 14, 23) : blank_pn;
  CHECK_THROWS_SOFT(const std::vector<PolyNumeric> spn_dbls =
                    amberPrmtopData(spreadsheet, 8, NumberFormat::INTEGER, 5, 14, 23), "Formatted "
                    "integer reading by amberPrmtopData() from " + tf.getFileName() + " passes "
                    "despite the numbers having general format.", pfile_check);

  // Test formatted data reading with general-format reals
  const double pn16d = (pfile_exists) ? pn_dbls[16].d : 0.0;
  check(pn16d, RelationalOperator::EQUAL, Approx(20.29945545).margin(1.0e-9), "General format "
        "real numbers read by amberPrmtopData() from " + tf.getFileName() + " are incorrect.",
        pfile_check);
  CHECK_THROWS_SOFT(const std::vector<PolyNumeric> pn_dbls_ii =
                    amberPrmtopData(spreadsheet, 8, NumberFormat::STANDARD_REAL, 5, 14, 24),
                    "Formatted real number reading by amberPrmtopData() from " + tf.getFileName() +
                    " passes despite there not being enough numbers in the file.", pfile_check);
  CHECK_THROWS_SOFT(const std::vector<PolyNumeric> pn_dbls_iii =
                    amberPrmtopData(spreadsheet, 8, NumberFormat::STANDARD_REAL, 5, 14, 10, 18),
                    "Formatted real number reading by amberPrmtopData() from " + tf.getFileName() +
                    " passes despite there being too many numbers in the file.", pfile_check);
  CHECK_THROWS_SOFT(const std::vector<PolyNumeric> pn_dbls_iv =
                    amberPrmtopData(spreadsheet, 8, NumberFormat::STANDARD_REAL, 4, 14, 23),
                    "Formatted real number reading by amberPrmtopData() from " + tf.getFileName() +
                    " passes despite there being too many numbers per line in the file.",
                    pfile_check);
  CHECK_THROWS_SOFT(const std::vector<PolyNumeric> pn_dbls_iv =
                    amberPrmtopData(spreadsheet, 8, NumberFormat::STANDARD_REAL, 6, 14, 23),
                    "Formatted real number reading by amberPrmtopData() from " +
                    spreadsheet.getFileName() + " passes despite there being insufficient numbers "
                    "per line in the file.", pfile_check);

  // Test scientific notation reading
  const std::vector<PolyNumeric> spn_dbls = (pfile_exists) ?
    amberPrmtopData(spreadsheet, 15, NumberFormat::SCIENTIFIC, 5, 14, 17) : blank_pn;
  std::vector<double> std_dbls = doubleFromPolyNumeric(spn_dbls);
  const double stdb3 = (pfile_exists) ? std_dbls[3] : 0.0;
  check(stdb3, RelationalOperator::EQUAL, -49.123, "Formatted scientific notation number "
        "reading by amberPrmtopData() from " + spreadsheet.getFileName() + " fails.", pfile_check);
  check(sum<double>(std_dbls), RelationalOperator::EQUAL, Approx(-249.7637).margin(1.0e-8),
        "Formatted scientific notation number reading by amberPrmtopData() from " +
        spreadsheet.getFileName() + " fails in the sum.", pfile_check);
  CHECK_THROWS_SOFT(const std::vector<PolyNumeric> spn_dbls2 =
                    amberPrmtopData(spreadsheet, 8, NumberFormat::SCIENTIFIC, 5, 14, 23),
                    "Formatted number reading in scientific notation by amberPrmtopData() from " +
                    spreadsheet.getFileName() + " passes despite the numbers having general "
                    "format.", pfile_check);

  // Try char4 data and operator overloading
  const std::vector<PolyNumeric> pn_chrs = (pfile_exists) ?
    amberPrmtopData(spreadsheet, 23, NumberFormat::CHAR4, 20, 4, 98) : blank_pn;
  std::string pn_example = (pfile_exists) ? char4ToString(pn_chrs[47].c4) : std::string("");
  check(pn_example, RelationalOperator::EQUAL, "CD1 ", "Vectorized char4 data read by "
        "amberPrmtopData() from " + spreadsheet.getFileName() + " does not appear to be correct.",
        pfile_check);
  pn_example = (pfile_exists) ? pn_example + pn_chrs[48].c4 : pn_example;
  check(pn_example, RelationalOperator::EQUAL, "CD1 HD1 ", "Operator += overloading produces "
        "incorrect results.", pfile_check);

  // Try reading single numbers from formatted data
  check(0.578, RelationalOperator::EQUAL, Approx(readRealValue("670.578924", 2, 5)), "Real value "
        "reading from a character string fails.");
  check(789, RelationalOperator::EQUAL, Approx(readRealValue("670.578924", 5, 3)), "Integer value "
        "reading from a character string fails.");
  const std::string sv_number("91855.293");
  check(55.2, RelationalOperator::EQUAL, Approx(readRealValue(sv_number, 3, 4)), "Real value "
        "reading from a Standard Template Library string fails.");
  check(185, RelationalOperator::EQUAL, Approx(readIntegerValue(sv_number, 1, 3)), "Integer value "
        "reading from a Standard Template Library string fails.");
  const std::string tf_base("Blah blah\n   0.5 456105\nThis is a text line.\n");
  const TextFile sv_tf(tf_base, TextOrigin::RAM);
  check(456, RelationalOperator::EQUAL, Approx(readRealValue(sv_tf, 1, 7, 3)), "Integer value "
        "interpretation from a TextFile object fails.");
  check(0.5, RelationalOperator::EQUAL, Approx(readRealValue(sv_tf, 1, 2, 4)), "Integer value "
        "interpretation from a TextFile object fails.");
  
  // Test custom string comparison functions
  section(6);
  const std::string tw_a("ThisWord");
  const std::string tw_b("thisword");
  const std::string tw_c("ThisWord T7J");
  const std::string tw_d("thisword t7j");
  check(strncmpCased(tw_a, tw_b, CaseSensitivity::YES) == false, "Case-sensitive "
        "string comparison fails to differentiate \"" + tw_a + "\" and \"" + tw_b + "\".");
  check(strncmpCased(tw_a, tw_b, CaseSensitivity::NO), "Case-insensitive "
        "string comparison still differentiates \"" + tw_a + "\" and \"" + tw_b + "\".");
  check(strncmpCased(tw_c, tw_d, CaseSensitivity::NO), "Case-insensitive "
        "string comparison still differentiates \"" + tw_c + "\" and \"" + tw_d + "\".");
  check(strcmpCased(tw_a, tw_b, CaseSensitivity::YES) == false, "Case-sensitive "
        "string comparison fails to differentiate \"" + tw_a + "\" and \"" + tw_b + "\".");
  check(strcmpCased(tw_a, tw_b, CaseSensitivity::NO), "Case-insensitive "
        "string comparison still differentiates \"" + tw_a + "\" and \"" + tw_b + "\".");
  check(strcmpCased(tw_c, tw_d, CaseSensitivity::NO), "Case-insensitive "
        "string comparison still differentiates \"" + tw_c + "\" and \"" + tw_d + "\".");
  check(strcmpCased(tw_a, tw_c, CaseSensitivity::YES) == false, "Case-sensitive unsized "
        "string comparison fails to differentiate \"" + tw_a + "\" and \"" + tw_d + "\".");
  check(strcmpCased(tw_a, tw_d, CaseSensitivity::NO) == false, "Case-insensitive unsized "
        "string comparison fails to differentiate \"" + tw_a + "\" and \"" + tw_d + "\".");
  const std::string tw_regexp_a("Th*W*?d");
  const std::vector<WildCardKind> wildcards_a = { WildCardKind::NONE, WildCardKind::NONE,
                                                  WildCardKind::FREE_STRETCH, WildCardKind::NONE,
                                                  WildCardKind::FREE_STRETCH,
                                                  WildCardKind::FREE_CHARACTER,
                                                  WildCardKind::NONE };
  check(strcmpWildCard(tw_a, tw_regexp_a, wildcards_a), "Case-sensitive string matching with "
        "wildcards fails to equate \"" + tw_a + "\" with \"" + tw_regexp_a + "\"."); 
  check(strcmpWildCard(tw_b, tw_regexp_a, wildcards_a) == false, "Case-sensitive string matching "
        "with wildcards equates \"" + tw_b + "\" with \"" + tw_regexp_a + "\".");
  const std::string tw_e("  This is a very long string with white space at either end.  ");
  const std::string tw_regexp_b(" * ");
  const std::vector<WildCardKind> wildcards_b = { WildCardKind::NONE, WildCardKind::FREE_STRETCH,
                                                  WildCardKind::NONE };
  const std::string tw_regexp_c("T*.");
  check(strcmpWildCard(tw_e, tw_regexp_c, wildcards_b), "Case-sensitive string matching "
        "with wildcards fails to equate \"" + tw_e + "\" with \"" + tw_regexp_c + "\".");
  const std::string tw_regexp_d("T*k");
  check(strcmpWildCard(tw_e, tw_regexp_d, wildcards_b) == false, "Case-sensitive string matching "
        "with wildcards incorrectly equates \"" + tw_e + "\" with \"" + tw_regexp_d + "\".");

  // Test how numbers and strings can be formatted with special considerations to tabulated,
  // terminal output.
  section(7);
  std::vector<std::string> many_params(10);
  many_params[0] = "Times Radio";
  many_params[1] = "October 15, 2024";
  many_params[2] = "Londoners in Shock at Outcome of Keating Summit";
  many_params[3] = "newsy";
  many_params[4] = "people";
  many_params[5] = "hearing";
  many_params[6] = "story";
  many_params[7] = "common";
  many_params[8] = "money";
  many_params[9] = "Downing";
  const int general_length = justifyStrings(&many_params, JustifyText::RIGHT, 10, 2, 4);
  check(general_length, RelationalOperator::EQUAL, 11, "The consensus length of a series of "
        "strings does not meet expectations after right-hand justification.");
  check(many_params[3].size(), RelationalOperator::EQUAL, 11, "A set of strings was not set to "
        "the consensus lengths as expected.");
  check(many_params[2].size(), RelationalOperator::EQUAL, 47, "A very long string should not "
        "have been altered by the formatStrings() function.");
  check(many_params[9], RelationalOperator::EQUAL, "    Downing", "Strings were not pre-pended "
        "with white space as expected after right-hand justification.");
  const std::string all_blank = "          ";
  check(removeLeadingWhiteSpace(all_blank).size(), RelationalOperator::EQUAL, 0, "A completely "
        "blank string was not reduced to zero length by removal of leading white space as "
        "expected.");
  check(removeTailingWhiteSpace(all_blank).size(), RelationalOperator::EQUAL, 0, "A completely "
        "blank string was not reduced to zero length by removal of tailing white space as "
        "expected.");
  check(addTailingWhiteSpace(tw_a, 12), RelationalOperator::EQUAL, "ThisWord    ", "Tailing "
        "whitespace was not added in the expected manner to a simple string.");
  const int lgen_length = justifyStrings(&many_params, JustifyText::LEFT, 10, 2, 4);
  check(lgen_length, RelationalOperator::EQUAL, 11, "The consensus length of a series of "
        "strings does not meet expectations after left-hand justification.");
  check(many_params[4], RelationalOperator::EQUAL, "people     ", "Strings were not appended with "
        "the expected white space after left-hand justification.");

  // Test the format determination of various number series.
  const std::vector<double> column_num_a = { 5.68, 67.45, 53.22, 0.005, 9.997 };
  const std::vector<double> column_num_b = { 67.5, 32.47846, 77.5, -99.499 };
  const std::vector<double> column_num_c = { 5.68, 3.45, 1.22, -0.5, 9.997 };
  const std::vector<double> column_num_d = { 5, 16, 7, 9, -2, 44 };
  const int align_fmt_a = findAlignmentWidth<double>(column_num_a, 3);
  const int align_fmt_b = findAlignmentWidth<double>(column_num_b, 1);
  const int align_fmt_c = findAlignmentWidth<double>(column_num_c, 3);
  const int align_fmt_d = findAlignmentWidth<double>(column_num_d, 0);
  const int align_fmt_e = findAlignmentWidth<double>(column_num_b, 0);
  check(align_fmt_a, RelationalOperator::EQUAL, 6, "The column format width of a series of "
        "numbers does not meet expectations.");
  check(align_fmt_b, RelationalOperator::EQUAL, 5, "The column format width of a series of "
        "numbers does not meet expectations when a negative number is included.");
  check(align_fmt_c, RelationalOperator::EQUAL, 6, "The column format width of a series of "
        "numbers does not meet expectations when one of the numbers is a negative fraction.");
  check(align_fmt_d, RelationalOperator::EQUAL, 2, "The column format width of a series of "
        "numbers does not meet expectations with integer data when the number of decimal places "
        "is zero.");
  check(align_fmt_e, RelationalOperator::EQUAL, 3, "The column format width of a series of "
        "numbers does not meet expectations with real data when the number of decimal places is "
        "zero.");
  
  // Check text protection for output purposes.
  const std::string my_paragraph = "This is a paragraph containing many words, including some "
                                   "which are particularly long such as "
                                   "\"antidisestablishmentarianism.\" Printing this much text may "
                                   "require multiple line breaks.  In fact, the text itself can "
                                   "contain line breaks:\n\nAll of the text is to be protected by "
                                   "a hash symbol.";
  const std::string protected_paragraph = protectText(my_paragraph, '#', 28);
  const std::string protected_ans = "# This is a paragraph\n# containing many words,\n"
                                    "# including some which are\n# particularly long such as\n"
                                    "# \"antidisestablishmentarianism.\"\n"
                                    "# Printing this much text\n# may require multiple line\n"
                                    "# breaks.  In fact, the text\n# itself can contain line\n"
                                    "# breaks:\n#\n# All of the text is to be\n"
                                    "# protected by a hash\n# symbol.\n";
  check(protected_paragraph, RelationalOperator::EQUAL, protected_ans, "A paragraph was not "
        "protected for output in the manner expected.");

  // Test conversion of base-10 integers to alphabetical base-26 strings.
  const std::vector<ullint> number_range = { 18, 29, 383, 19836732 };
  std::string string_range;
  for (int i = 0; i < 4; i++) {
    string_range += alphabetNumber(number_range[i]);
    if (i < 3) {
      string_range += "_";
    }
  }
  const std::string string_range_ans("s_ad_nt_aqjpgg");
  check(string_range, RelationalOperator::EQUAL, string_range_ans, "Positive integers were not "
        "converted to alphanumeric strings as expected.");
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
