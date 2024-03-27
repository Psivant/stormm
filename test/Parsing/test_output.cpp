#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/FileManagement/file_enumerators.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Parsing/parsing_enumerators.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/ordered_list.h"
#include "../../src/Reporting/reporting_enumerators.h"
#include "../../src/Reporting/report_table.h"
#include "../../src/Reporting/section_contents.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/unit_test_enumerators.h"

using namespace stormm::diskutil;
using namespace stormm::errors;
using namespace stormm::parse;
using namespace stormm::random;
using namespace stormm::review;
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

  // Section 1
  section("Ordered list creation");

  // Section 2
  section("Table creation");

  // Section 3
  section("Assembling output sections");

  // Section 4
  section("Assembling a complete output file");
  
  // Create an ordered list and test its output
  section(1);
  OrderedList list_a(ListEnumeration::BULLET);
  list_a.addItem("The first item in the list");
  list_a.addItem("The second item in the list");
  list_a.addItem("The third item in the list");
  check(list_a.getItemCount(), RelationalOperator::EQUAL, 3, "The number of items in an ordered "
        "list is incorrect.");
  check(list_a.getItem(2), RelationalOperator::EQUAL,
        std::string("  * The third item in the list\n"), "An ordered list does not produce the "
        "expected string when queried for one of its items.");
  list_a.addItemBefore("The first item in the list", "The second");
  check(list_a.getItem(1), RelationalOperator::EQUAL,
        std::string("  * The first item in the list\n"), "An ordered list does not produce the "
        "expected string after using the addItemBefore() method.");
  list_a.addNestedItem("A nested item beneath the second item, which calls itself the first.", 1);
  check(list_a.getNestedItemCount(), RelationalOperator::EQUAL, 1, "The presence of a nested "
        "item is not registered in an ordered list.");
  check(list_a.getNestedItem(1, 0).substr(0, 16), RelationalOperator::EQUAL, "    - A nested i",
        "An ordered list does not produce the expected string when queried for one of its nested "
        "items.");
  list_a.addItemAfter("Item B", 0);
  check(list_a.getNestedItem(2, 0).substr(0, 16), RelationalOperator::EQUAL, "    - A nested i",
        "An ordered list does not produce the expected string when queried for one of its nested "
        "items, following insertion of a new main list item.");
  CHECK_THROWS(list_a.addItemAfter("Item C", 5), "A list of five items allowed a new item to be "
               "added after the sixth.");
  CHECK_THROWS(list_a.getNestedItem(3, 0), "A nested item was produced for a main list item that "
               "has no nested items.");
  CHECK_THROWS(list_a.getNestedItem(1, 0), "A nested item was produced for a main list item that "
               "no longer has any nested items.");
  list_a.setNestedBullet(':');
  check(list_a.getNestedItem(2, 0).substr(0, 16), RelationalOperator::EQUAL, "    : A nested i",
        "An ordered list does properly update its nested bullet item.");
  list_a.setStyle(ListEnumeration::NUMBERED);
  check(list_a.getItem(2), RelationalOperator::EQUAL,
        std::string("  3) The first item in the list\n"), "An ordered list does not produce the "
        "expected string after changing the enumeration style.");
  list_a.setNestedStyle(ListEnumeration::ALPHABETIC);
  list_a.addNestedItem("Nested i", 1);
  list_a.addNestedItem("Nested ii", 1);
  list_a.addNestedItem("Nested iii", 1);
  list_a.addNestedItem("Nested iv", 1);
  list_a.addNestedItem("Nested v", 1);
  list_a.addNestedItem("Bi-nested i", 4);
  list_a.addNestedItem("Bi-nested ii", 4);
  list_a.addNestedItem("Tri-nested x", 0);
  check(list_a.getNestedItemCount(), RelationalOperator::EQUAL, 9, "After numerous additions, the "
        "list no longer as the expected number of nested items.");
  const std::vector<int> nested_item_counts_ans = { 1, 5, 1, 0, 2 };
  std::vector<int> nested_item_counts(list_a.getItemCount());
  for (int i = 0; i < list_a.getItemCount(); i++) {
    nested_item_counts[i] = list_a.getNestedItemCount(i);
  }
  check(nested_item_counts, RelationalOperator::EQUAL, nested_item_counts_ans, "The numbers of "
        "nested items distributed over the list do not meet expectations.");
  const char osc = osSeparator();
  const std::string snp_file = oe.getStormmSourcePath() + osc + "test" + osc + "Parsing" + osc +
                               "output_snapshot.txt";
  const bool snpfi_exists = (getDrivePathType(snp_file) == DrivePathType::FILE);
  const TestPriority do_snp_tests = (snpfi_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (snpfi_exists == false && oe.takeSnapshot() == SnapshotOperation::COMPARE) {
    rtWarn("The snapshot file " + snp_file + " was not found.  Check the STORMM source path "
           "environemnt variable, currently " + oe.getStormmSourcePath() + ", to ensure its "
           "validity and that is contains the test" + osc + "Parsing" + osc + " directory.  Tests "
           "depending on this file will be skipped.", "test_output");
  }
  const std::string list_a_str = list_a.printList(80);
  snapshot(snp_file, list_a_str, "list_basic_output", "A list compiled from previous tests does "
           "print in the expected form.", oe.takeSnapshot(), PrintSituation::OVERWRITE,
           do_snp_tests);
  const TextFile list_a_tf = list_a.printList(OutputSyntax::STANDALONE, 80);
  snapshot(snp_file, list_a_tf, "list_prot_output", "A list compiled from previous tests does not "
           "print in the expected form when protected.", oe.takeSnapshot(), PrintSituation::APPEND,
           do_snp_tests);
  const std::string block_text("This blocked text will run way over forty characters, for the "
                               "purpose of testing how the STORMM libraries can protect a stream "
                               "of text, with or without initial indentation.  Typing more text "
                               "affords larger numbers of lines, and thus more chances to see a "
                               "word wrap differently (based on the initial indentation).  This "
                               "effect is especially pronounced if there is a spread of long and "
                               "short words.");
  const std::string frmt_block_text = indentText(block_text, 0, 40);
  snapshot(snp_file, frmt_block_text, "frmt_block", "A block of plain formatted text does not "
           "print as expected.", oe.takeSnapshot(), PrintSituation::APPEND, do_snp_tests);
  const std::string prot_block_text = protectText(block_text, '|', 40);
  snapshot(snp_file, prot_block_text, "prot_block", "A block of comment-protected text does not "
           "print as expected.", oe.takeSnapshot(), PrintSituation::APPEND, do_snp_tests);
  const std::string idnt_block_text = indentText(block_text, 2, 38);
  snapshot(snp_file, idnt_block_text, "idnt_block", "A block of indented, formatted text does not "
           "print as expected.", oe.takeSnapshot(), PrintSituation::APPEND, do_snp_tests);
  const std::string prid_block_text = protectText(idnt_block_text, '%', 40);
  snapshot(snp_file, prid_block_text, "prid_block", "A block of comment-protected text with "
           "pre-existing indentation does not print as expected.", oe.takeSnapshot(),
           PrintSituation::APPEND, do_snp_tests);
  CHECK_THROWS(const std::string& bad_str = toRoman(maximum_roman_numeral + 1), "An integer "
               "outside the allowed range was converted to a Roman numeral.");
  CHECK_THROWS(const std::string& bad_str = toRoman(-5), "A negative number was converted to a "
               "Roman numeral.");
  check(toRoman(1), RelationalOperator::EQUAL, "i", "Roman numeral conversion produces an "
        "unexpected result.");
  check(toRoman(24), RelationalOperator::EQUAL, "xxiv", "Roman numeral conversion produces an "
        "unexpected result.");
  check(toRoman(maximum_roman_numeral + 4, ExceptionResponse::SILENT), RelationalOperator::EQUAL,
        " ", "A Roman numeral outside the allowed range should return white space in "
        "fault-tolerant conditions.");
  
  // Create a table of results and test its output.
  section(2);
  Xoshiro256ppGenerator xrs;
  const std::vector<double> dbl_uni_data = uniformRand(&xrs, 20, 1.0);
  const std::vector<std::string> dbl_cols = { "Column A", "Segment B", "Partition C", "Group D" };
  const std::vector<int> dbl_dp = { 2, 3, 4, 3 };
  const ReportTable rtab_a(dbl_uni_data, dbl_cols, dbl_dp, "uni_dbl");
  const std::string rtab_a_str = rtab_a.printTable(OutputSyntax::MATRIX_PKG, 100);
  snapshot(snp_file, rtab_a_str, "rand_uni_table", "A table of double-precision random numbers "
           "uniformly selected on the interval [0, 1) was not printed to string as expected.",
           oe.takeSnapshot(), PrintSituation::APPEND, do_snp_tests);
  const std::vector<double> dbl_gss_data = gaussianRand(&xrs, 20, 1.0);
  const ReportTable rtab_b(dbl_gss_data, dbl_cols, dbl_dp, "gss_dbl");
  const std::string rtab_b_str = rtab_b.printTable(OutputSyntax::MATRIX_PKG, 100);
  snapshot(snp_file, rtab_b_str, "rand_gss_table", "A table of double-precision random numbers "
           "selected on the normal distribution was not printed to string as expected.",
           oe.takeSnapshot(), PrintSituation::APPEND, do_snp_tests);

  // Create additional tables 
  std::vector<int> int_uni_data(20);
  for (int i = 0; i < 20; i++) {
    int_uni_data[i] = static_cast<int>(50.0 * (xrs.uniformRandomNumber() - 0.5));
  }
  const std::vector<std::string> int_cols = { "Part A", "Part B", "Sanskrit Domesticity", "D" };
  const ReportTable rtab_c(int_uni_data, int_cols, "uni_int");
  const std::string rtab_c_str = rtab_c.printTable(OutputSyntax::MATRIX_PKG, 100);
  snapshot(snp_file, rtab_c_str, "rand_int_table", "A table of random integers selected over the "
           "selected on the normal distribution was not printed to string as expected.",
           oe.takeSnapshot(), PrintSituation::APPEND, do_snp_tests);
  const TextFile rtab_c_tf(rtab_c_str, TextOrigin::RAM);
  snapshot(snp_file, rtab_c_tf, "rand_int_table", "A table of random integers converted to a "
           "TextFile was not printed as expected.", oe.takeSnapshot(), PrintSituation::APPEND,
           do_snp_tests);

  // Assemble the results into the contents of a mock output section and print.
  section(3);
  SectionContents scon;
  scon.setTitle("Mock section with an excessively long title to check the wrapping behavior of "
                "the header and its associated section number");
  scon.reserve(SectionComponent::NARRATIVE, 4);
  scon.reserve(SectionComponent::LIST, 1);
  scon.reserve(SectionComponent::TABLE, 1);
  const std::string explanation_a("The SectionContents object organizes the results of multiple "
                                  "narrative paragraphs, lists, and tables into a coherent, "
                                  "formatted result that is amenable to a plotting program of the "
                                  "user's choice.  Many SectionContents objects can be strung "
                                  "together in an array and then printed as a coherent output "
                                  "file, but they can also be printed individually.");
  scon.addNarration(explanation_a);
  check(scon.getComponentCount(), RelationalOperator::EQUAL, 1, "A section with one narrative "
        "paragraph reports the wrong number of components.");
  const std::string explanation_b("In general, sections begin with a double horizontal rule "
                                  "('==...==' preceded by the appropriate comment character) "
                                  "while subsections begin with a single horizontal rule "
                                  "('--...--', again preceded by the appropriate comment "
                                  "character).  Different section components are spaced by a "
                                  "single blank line, but no blank space is placed between the "
                                  "leading horizontal rule and the first component.");
  scon.addNarration(explanation_b);
  check(scon.getComponentCount(SectionComponent::NARRATIVE), RelationalOperator::EQUAL, 2,
        "A SectionContents object does not report the expected number of narrative components.");
  OrderedList list_b(ListEnumeration::ROMAN, ListEnumeration::NUMBERED);
  list_b.addItem("Lists can be effective ways to communicate the meaning of output.  The "
                 "wrapping of list items is as important as the wrapping of other text, and must "
                 "also respect the maximum marker sizes of the list items.");
  list_b.addNestedItem("Nesting provides further organization.");
  list_b.addNestedItem("Only one layer of nesting is available.");
  list_b.addItem("A paragraph of narration will typically lead into a list.");
  list_b.addItem("Lists in the output are typically indented to set them apart.");
  list_b.addNestedItem("Nested list items are further indented.");
  list_b.addNestedItem("Nested items follow their own numbering scheme.");
  list_b.addItem("The entirety of a list will be protected within a comment block.");
  scon.addList(list_b);
  scon.addNarration("Narration can also precede a table in the output.");
  const std::vector<std::string> exmp_cols = { "Column A", "Column B", "Partition C", "Segment D",
                                               "Column E", "Group F", "Part G", "H", "I", "J" };
  const std::vector<double> exmp_data = gaussianRand(&xrs, 50, 1.0);
  const std::vector<int> exmp_dp = { 2, 3, 4, 2, 7, 8, 7, 6, 6, 6 };
  const ReportTable exmp_rtab(exmp_data, exmp_cols, exmp_dp, "pseudo_rngs");
  scon.addTable(exmp_rtab);
  scon.addNarration("The tables will be printed within the allowed, overall format width of the "
                    "output file, which itself will be user-specifiable.  If columns do not fit "
                    "in the space allowed, they will be printed further down.  Different columns "
                    "can be printed with different numbers of decimal places, if desired.  Column "
                    "headings will be broken into separate words and stacked on multiple lines, "
                    "as is reasonable given the width of data elements, to keep the table compact "
                    "and tidy.");
  const std::vector<OutputSyntax> all_formats = { OutputSyntax::MATPLOTLIB,
                                                  OutputSyntax::MATRIX_PKG,
                                                  OutputSyntax::STANDALONE };
  const std::vector<std::string> fmt_codes = { "mplb", "mpkg", "stdl" };
  for (size_t i = 0; i < all_formats.size(); i++) {
    snapshot(snp_file, scon.sectionAsString(all_formats[i]), "section_contents_" + fmt_codes[i],
             "The SectionContents object does not produce the output expected for " +
             getEnumerationName(all_formats[i]) + ".", oe.takeSnapshot(), PrintSituation::APPEND,
             do_snp_tests);
  }

  // Combine multiple sections into a single output file
  section(4);
  SectionContents sca, scb, scb_i, scb_ii, scc;
  sca.setTitle("Primary");
  sca.addNarration("The first section, it stands alone.");
  OrderedList list_c(ListEnumeration::ROMAN, ListEnumeration::NUMBERED);
  list_c.addItem("Item one.");
  list_c.addNestedItem("Sub-item the first.");
  list_c.addNestedItem("Sub-item the second.");
  list_c.addItem("Item two.");
  list_c.addItem("Item three.");
  sca.addList(list_c);
  sca.addNarration("Closing statement for the first section.");
  scb.setTitle("Secondary, Main");
  scb.addNarration("This section will have sub-sections.  Sub-sections do not get indented "
                   "further than their parent sections but they do have unique borders.");
  scb_i.setTitle("Secondary, First");
  scb_i.designateSubsection();
  scb_i.addNarration("This section contains a small table of randoms.");
  scb_i.addTable(rtab_a);
  scb_i.addNarration("And, perhaps, something to close it out.");
  scb_ii.setTitle("Secondary, Second");
  scb_ii.designateSubsection();
  scb_ii.addNarration("Does anything more need to be said?");
  scc.setTitle("Tertiary");
  scc.addNarration("Any components to this section should be printed as per usual, with no "
                   "differences in style from prior sections.");
  const std::vector<SectionContents> all_sections = { sca, scb, scb_i, scb_ii, scc };
  const std::string all_sectstr = printAllSections(all_sections, OutputSyntax::MATRIX_PKG);
  snapshot(snp_file, all_sectstr, "multi_section", "A composition of multiple sections, with "
           "splash and watermark, does not print as expected.", oe.takeSnapshot(),
           PrintSituation::APPEND, do_snp_tests);
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
