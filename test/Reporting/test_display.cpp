#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Reporting/present_field.h"
#include "../../src/Reporting/render_options.h"
#include "../../src/Structure/background_mesh.h"
#include "../../src/Structure/mesh_parameters.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/UnitTesting/dissect_textfile.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::structure;
using namespace stormm::testing;
using namespace stormm::topology;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main (const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  StopWatch timer;
  const int file_read_wtm  = timer.addCategory("File Reading");
  const int mesh_build_wtm = timer.addCategory("Mesh Generation");
  const int file_write_wtm = timer.addCategory("File Writing");
  
  // Section 1
  section("Test molecule display");

  // Create a receptor, a ligand, a field, and a scene with them all
  const char osc = osSeparator();
  const std::string topology_base = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string coordinate_base = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::vector<std::string> molecule_names = { "tamavidin", "biotin" };
  TestSystemManager tsm(topology_base, "top", molecule_names, coordinate_base, "inpcrd",
                        molecule_names);
  const CoordinateFrame tama_cf = tsm.exportCoordinateFrame(0);
  const CoordinateFrame biot_cf = tsm.exportCoordinateFrame(1);
  const CoordinateFrameReader biot_cfr = biot_cf.data();
  AtomGraph tama_ag = tsm.exportAtomGraph(0);
  tama_ag.modifyAtomMobility(0, tama_ag.getAtomCount(), MobilitySetting::OFF);
  if (oe.getDisplayTimingsOrder()) timer.assignTime(file_read_wtm);
  const MeshParameters mps(9, 9, 9, minValue(biot_cfr.xcrd, biot_cfr.natom) - 4.0,
                           minValue(biot_cfr.ycrd, biot_cfr.natom) - 4.0,
                           minValue(biot_cfr.zcrd, biot_cfr.natom) - 4.0, 2.0, 40,
                           Interpolant::FUNCTION_VALUE);
  const BackgroundMesh<double> bgm(GridDetail::NONBONDED_FIELD, NonbondedPotential::VAN_DER_WAALS,
                                   tama_ag, tama_cf, 3.1, 0.15,
                                   VdwCombiningRule::LORENTZ_BERTHELOT, mps, {}, {}, 0.0, 1.0,
                                   PrecisionModel::DOUBLE);
  if (oe.getDisplayTimingsOrder()) timer.assignTime(mesh_build_wtm);
  const MeshParameters mps_report(37, 37, 37, minValue(biot_cfr.xcrd, biot_cfr.natom) - 4.0,
                                  minValue(biot_cfr.ycrd, biot_cfr.natom) - 4.0,
                                  minValue(biot_cfr.zcrd, biot_cfr.natom) - 4.0, 0.5 - 1.0e-8, 40,
                                  Interpolant::FUNCTION_VALUE);
  RenderOptions ropt;
  ropt.setPotentialFieldBorderDisplay(true);
  const double4 iso_a_color = { 1.0, 0.4, 0.0, 1.0 };
  const uchar4  iso_b_color = {  64,   0, 192, 255 };
  ropt.addIsosurface(-1.5, iso_a_color, SurfaceRender::WIRE);
  ropt.addIsosurface(-2.0, iso_b_color, SurfaceRender::SOLID);
  const std::string scene_file = oe.getTemporaryDirectoryPath() + osc + "test_mesh.m";
  printToFile(bgm, mps_report, scene_file, PrintSituation::OPEN_NEW, GridFileSyntax::MATRIX_PKG,
              ropt);
  if (oe.getDisplayTimingsOrder()) timer.assignTime(file_write_wtm);
  oe.logFileCreated(scene_file);
  const std::string scene_data_file = getSceneDataFileName(scene_file);
  oe.logFileCreated(scene_data_file);

  // The above output entails a significant output file that could potentially change in minute
  // ways which are irrelevant to its overall function.  However, certain features of the output
  // can be expected to remain constant in order for it to work as expected.  Use parsing methods
  // in the unit testing library to read the output as a text file and verify that it has the
  // correct features.
  const TextFile scene_display_tf(scene_file);
  const TextFile scene_data_tf(scene_data_file);
  const double2 numbers_per_line = meanWordsPerLine(scene_data_tf);
  check(numbers_per_line.x, RelationalOperator::EQUAL, 37.0, "The number of values placed per "
        "line in the potential field otput for a " +
        getEnumerationName(GridFileSyntax::MATRIX_PKG) + " data file is inconsistent with the "
        "expected result.", tsm.getTestingStatus());
  check(numbers_per_line.y, RelationalOperator::EQUAL, 0.0, "The number of words per line in the "
        "potential field output of a " + getEnumerationName(GridFileSyntax::MATRIX_PKG) +
        " data file must be consistent in order for it to be read correctly.",
        tsm.getTestingStatus());
  const std::vector<std::string> reshaping = { "dock_field", "=", "reshape", "(", "dock_field",
                                               ",", "na", ",", "nb", ",", "nc", ")" };
  const std::vector<char> matrix_pkg_separators = {',', '(', ')', '*', ';', '[', ']', '\n' };
  const int n_reshaping = countInstances(scene_display_tf, reshaping, matrix_pkg_separators);
  check(n_reshaping, RelationalOperator::EQUAL, 1, "The reshaping command for displaying a "
        "docking field was not found as expected.");
  const std::vector<std::string> border_path = { "border_path", "=", "(", "invU", "*",
                                                 "border_path'", ")", "'" };
  const int n_border_path = countInstances(scene_display_tf, border_path, matrix_pkg_separators);
  check(n_border_path, RelationalOperator::EQUAL, 1, "The border path generation command was not "
        "found as expected.", tsm.getTestingStatus());
  const TextFile plot3_grep = grep(scene_display_tf, "(plot3)(.*)");
  check(plot3_grep.getLineCount(), RelationalOperator::EQUAL, 7, "The number of plot3() commands "
        "found in matrix package output script did notmeet expectations.", tsm.getTestingStatus());
  
  if (oe.getDisplayTimingsOrder()) timer.assignTime(0);
  
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
