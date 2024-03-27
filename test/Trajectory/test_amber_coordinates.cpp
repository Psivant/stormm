#include "copyright.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/coordinate_series.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/trajectory_enumerators.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::constants::small;
using stormm::data_types::llint;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::diskutil::PrintSituation;
using stormm::errors::rtWarn;
using stormm::stmath::computeBoxTransform;
using stormm::stmath::elementwiseDivide;
using stormm::stmath::elementwiseMultiply;
using stormm::stmath::extractBoxDimensions;
using stormm::stmath::mean;
using stormm::parse::NumberFormat;
using stormm::parse::polyNumericVector;
using stormm::random::Xoroshiro128pGenerator;
using stormm::random::Ran2Generator;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using stormm::symbols::pi;
using stormm::topology::UnitCellType;
using namespace stormm::testing;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// A function that creates a PhaseSpace object filled with random coordinates normally distributed
// about the origin.
//
// Arguments:
//   particle_count:  The number of particles to create
//   sigma:           Width of the normal distribution
//   igseed:          Random number generator seed
//-------------------------------------------------------------------------------------------------
PhaseSpace createRandomParticles(const int particle_count, const double sigma = 5.0,
                                 const int igseed = -1) {
  const int actual_seed = (igseed == -1) ? particle_count : igseed;
  Xoroshiro128pGenerator xrs_rng(actual_seed);
  PhaseSpace result(particle_count);
  PhaseSpaceWriter psw = result.data();
  for (int i = 0; i < particle_count; i++) {
    psw.xcrd[i] = sigma * xrs_rng.gaussianRandomNumber();
    psw.ycrd[i] = sigma * xrs_rng.gaussianRandomNumber();
    psw.zcrd[i] = sigma * xrs_rng.gaussianRandomNumber();
  }
  return result;
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
  section("Coordinate reading");

  // Section 2
  section("Coordinate manipulation");

  // Section 3
  section("Coordinate object integration");

  // Section 4
  section("Coordinate object C++ handling");

  // Section 5
  section("Coordinate series handling");
  
  // If one or more of the coordinate files is missing, this will cause tests to be skipped.
  // Do not report immediately, as one missing file likely means other missing files.  Wait
  // until the end and then alert the user that some tests were aborted.
  bool missing_files = false;

  // Test (delayed) file reading for an input coordinate set
  section(1);
  const char osc = osSeparator();
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string tip3p_crd_name = base_crd_name + osc + "tip3p.inpcrd";
  const bool tip3p_exists = (getDrivePathType(tip3p_crd_name) == DrivePathType::FILE);
  const TestPriority do_tip3p = (tip3p_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  PhaseSpace tip3p;
  if (tip3p_exists) {
    tip3p.buildFromFile(tip3p_crd_name, CoordinateFileKind::AMBER_INPCRD);    
  }
  else {
    missing_files = true;
  }
  PhaseSpaceWriter tip3pwr = tip3p.data();
  check(tip3pwr.natom, RelationalOperator::EQUAL, 768, "The TIP3P Amber input coordinates "
        "file " + tip3p.getFileName() + " contains the wrong number of atoms.", do_tip3p);
  const CoordinateFileKind t3p_kind = (tip3p_exists) ? detectCoordinateFileKind(tip3p_crd_name) :
                                                       CoordinateFileKind::UNKNOWN;
  check(t3p_kind == CoordinateFileKind::AMBER_INPCRD, "Coordinate file type detection did not "
        "succeed for an Amber ASCII input coordinate file.  The file was marked " +
        getEnumerationName(t3p_kind) + ".", do_tip3p);
  double oh_mean = 0.0;
  if (tip3p_exists) {
    for (int i = 0; i < tip3pwr.natom; i += 3) {
      const double oh1x = tip3pwr.xcrd[i + 1] - tip3pwr.xcrd[i];
      const double oh1y = tip3pwr.ycrd[i + 1] - tip3pwr.ycrd[i];
      const double oh1z = tip3pwr.zcrd[i + 1] - tip3pwr.zcrd[i];
      const double oh2x = tip3pwr.xcrd[i + 2] - tip3pwr.xcrd[i];
      const double oh2y = tip3pwr.ycrd[i + 2] - tip3pwr.ycrd[i];
      const double oh2z = tip3pwr.zcrd[i + 2] - tip3pwr.zcrd[i];
      oh_mean += sqrt((oh1x * oh1x) + (oh1y * oh1y) + (oh1z * oh1z));
      oh_mean += sqrt((oh2x * oh2x) + (oh2y * oh2y) + (oh2z * oh2z));
    }
    oh_mean /= 2.0 * static_cast<double>(tip3pwr.natom / 3);
  }
  check(oh_mean, RelationalOperator::EQUAL, Approx(0.9572).margin(1.0e-3), "O-H bond lengths in "
        "TIP3P water deviate too wildly from the prescribed 0.9572 Angstroms.", do_tip3p);
  check(tip3pwr.ycrd[32], RelationalOperator::EQUAL, Approx(7.9019706).margin(1.0e-7),
        "Coordinates read from an Amber .inpcrd file are incorrect.", do_tip3p);
  std::vector<double> gdims(6);
  for (int i = 0; i < 3; i++) {
    gdims[i] = 23.3662499;
    gdims[i + 3] = 109.4712190 * pi / 180.0;
  }
  std::vector<double> box_xfrm(9);
  std::vector<double> inv_xfrm(9);
  computeBoxTransform(gdims, &box_xfrm, &inv_xfrm);
  const std::vector<double> t3box = (tip3p_exists) ? tip3p.getBoxSpaceTransform() :
                                                     std::vector<double>(9, 0.0);
  check(t3box, RelationalOperator::EQUAL, Approx(box_xfrm).margin(1.0e-7), "The box "
        "transformation was not computed properly for " + tip3p.getFileName() + ".", do_tip3p);
  const std::vector<double> t3inv = (tip3p_exists) ? tip3p.getInverseTransform() :
                                                     std::vector<double>(9, 0.0);
  check(t3inv, RelationalOperator::EQUAL, Approx(inv_xfrm).margin(1.0e-7), "The inverse "
        "transformation was not computed properly for " + tip3p.getFileName() + ".", do_tip3p);

  // Test the copy constructor.  Stretch the original object after copying to ensure independence
  // of the copy.
  section(2);
  const PhaseSpace tip3p_copy = tip3p;
  const PhaseSpaceReader tip3pcr = tip3p_copy.data();
  oh_mean = 0.0;
  double orig_oh_mean = 0.0;
  if (tip3p_exists) {
    for (int i = 0; i < tip3pcr.natom; i++) {
      tip3pwr.xcrd[i] *= 1.5;
      tip3pwr.ycrd[i] *= 1.5;
      tip3pwr.zcrd[i] *= 1.5;
    }
    for (int i = 0; i < tip3pcr.natom; i += 3) {
      const double oh1x = tip3pcr.xcrd[i + 1] - tip3pcr.xcrd[i];
      const double oh1y = tip3pcr.ycrd[i + 1] - tip3pcr.ycrd[i];
      const double oh1z = tip3pcr.zcrd[i + 1] - tip3pcr.zcrd[i];
      const double oh2x = tip3pcr.xcrd[i + 2] - tip3pcr.xcrd[i];
      const double oh2y = tip3pcr.ycrd[i + 2] - tip3pcr.ycrd[i];
      const double oh2z = tip3pcr.zcrd[i + 2] - tip3pcr.zcrd[i];
      oh_mean += sqrt((oh1x * oh1x) + (oh1y * oh1y) + (oh1z * oh1z));
      oh_mean += sqrt((oh2x * oh2x) + (oh2y * oh2y) + (oh2z * oh2z));
    }
    oh_mean /= 2.0 * static_cast<double>(tip3pcr.natom / 3);
    for (int i = 0; i < tip3pwr.natom; i += 3) {
      const double oh1x = tip3pwr.xcrd[i + 1] - tip3pwr.xcrd[i];
      const double oh1y = tip3pwr.ycrd[i + 1] - tip3pwr.ycrd[i];
      const double oh1z = tip3pwr.zcrd[i + 1] - tip3pwr.zcrd[i];
      const double oh2x = tip3pwr.xcrd[i + 2] - tip3pwr.xcrd[i];
      const double oh2y = tip3pwr.ycrd[i + 2] - tip3pwr.ycrd[i];
      const double oh2z = tip3pwr.zcrd[i + 2] - tip3pwr.zcrd[i];
      orig_oh_mean += sqrt((oh1x * oh1x) + (oh1y * oh1y) + (oh1z * oh1z));
      orig_oh_mean += sqrt((oh2x * oh2x) + (oh2y * oh2y) + (oh2z * oh2z));
    }
    orig_oh_mean /= 2.0 * static_cast<double>(tip3pwr.natom / 3);
  }
  check(oh_mean, RelationalOperator::EQUAL, Approx(0.9572).margin(1.0e-3), "The PhaseSpace "
        "object's copy constructor does not create an independent copy.", do_tip3p);
  check(orig_oh_mean, RelationalOperator::EQUAL, Approx(1.4358).margin(1.0e-3), "Manipulating a "
        "PhaseSpace object by scaling its coordinates does not produce the expected effect.",
        do_tip3p);
  const std::vector<double> cp3inv = (tip3p_exists) ? tip3p_copy.getInverseTransform() :
                                                      std::vector<double>(9, 0.0);
  check(cp3inv, RelationalOperator::EQUAL, Approx(inv_xfrm).margin(1.0e-7), "The inverse "
        "transformation was not computed properly for " + tip3p.getFileName() + ".", do_tip3p);
  check(tip3p_copy.getUnitCellType() == UnitCellType::TRICLINIC, "The unit cell type of " +
        tip3p.getFileName() + " is incorrectly interpreted.", do_tip3p);

  // Reset the original tip3p PhaseSpace object
  if (tip3p_exists) {
    for (int i = 0; i < tip3pcr.natom; i++) {
      tip3pwr.xcrd[i] /= 1.5;
      tip3pwr.ycrd[i] /= 1.5;
      tip3pwr.zcrd[i] /= 1.5;
    }
  }

  // Try reading some additional coordinate files: trpcage in isolated boundary conditions
  section(1);
  const std::string trpcage_crd_name = base_crd_name + osc + "trpcage.inpcrd";
  const bool trpcage_exists = (getDrivePathType(trpcage_crd_name) == DrivePathType::FILE);
  const TestPriority do_trpcage = (trpcage_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  PhaseSpace trpcage;
  if (trpcage_exists) {
    trpcage.buildFromFile(trpcage_crd_name, CoordinateFileKind::AMBER_INPCRD);
  }
  else {
    missing_files = true;
  }
  check(trpcage.getAtomCount(), RelationalOperator::EQUAL, 304, "The number of atoms in the "
        "Trp-cage (isolated boundary conditions) system is incorrect.", do_trpcage);
  check(trpcage.getUnitCellType() == UnitCellType::NONE, "Isolated boundary conditions are not "
        "properly interpreted in the Trp-cage system.", do_trpcage);
  PhaseSpaceWriter trpcage_wr = trpcage.data();
  check(trpcage_wr.zcrd[51], RelationalOperator::EQUAL, Approx(4.233).margin(1.0e-7),
        "Coordinates of the Trp-cage system (isolated boundary conditions) are incorect.",
        do_trpcage);
  const std::string trpcage_iso_snapshot = base_crd_name + osc + "trpcage_iso_sample.m";
  const std::string tip5p_vel_snapshot = base_crd_name + osc + "tip5p_velocity_sample.m";
  const bool snps_exist = (getDrivePathType(trpcage_iso_snapshot) == DrivePathType::FILE &&
                           getDrivePathType(tip5p_vel_snapshot) == DrivePathType::FILE);
  const TestPriority snap_check = (snps_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  snapshot(trpcage_iso_snapshot, polyNumericVector(trpcage.getInterlacedCoordinates(20, 100)),
           "trp_part_crd", 1.0e-4, "Reconstructed partial coordinates of the Trp-cage system are "
           "incorrect.", oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL,
           PrintSituation::OVERWRITE, snap_check);

  // Read TIP5P in a periodic, cubic unit cell.  It has velocities.
  const std::string tip5p_crd_name = base_crd_name + osc + "tip5p.rst";
  const bool tip5p_exists = (getDrivePathType(trpcage_crd_name) == DrivePathType::FILE);
  const TestPriority do_tip5p = (tip5p_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  PhaseSpace tip5p;
  if (tip5p_exists) {
    tip5p.buildFromFile(tip5p_crd_name, CoordinateFileKind::AMBER_ASCII_RST);
  }
  else {
    missing_files = true;
  }
  check(tip5p.getAtomCount(), RelationalOperator::EQUAL, 1080, "The number of atoms in the TIP5P "
        "system is incorrect.", do_tip5p);
  const CoordinateFileKind t5p_kind = (tip5p_exists) ? detectCoordinateFileKind(tip5p_crd_name) :
                                                       CoordinateFileKind::UNKNOWN;
  check(t5p_kind == CoordinateFileKind::AMBER_ASCII_RST, "Coordinate file type detection did not "
        "succeed for an Amber ASCII restart file.  The file was marked " +
        getEnumerationName(t5p_kind) + ".", do_tip5p);
  PhaseSpaceWriter tip5pwr = tip5p.data();
  check(tip5pwr.unit_cell == UnitCellType::ORTHORHOMBIC, "An orthorhombic unit cell is not "
        "correctly interpreted in the TIP5P system.", do_tip5p);
  snapshot(tip5p_vel_snapshot,
           polyNumericVector(tip5p.getInterlacedCoordinates(20, 100, TrajectoryKind::VELOCITIES)),
           "t5p_part_vel", 1.0e-4, "Reconstructed partial velocities of the TIP5P system are "
           "incorrect.", oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL,
           PrintSituation::OVERWRITE, snap_check);
  
  // Read ubiquitin, solvated in a very small octahedral cell (there is an odd number of atoms)
  const std::string ubiquitin_crd_name = base_crd_name + osc + "ubiquitin.inpcrd";
  const bool ubiquitin_exists = (getDrivePathType(trpcage_crd_name) == DrivePathType::FILE);
  const TestPriority do_ubiquitin = (ubiquitin_exists) ? TestPriority::CRITICAL :
                                                         TestPriority::ABORT;
  PhaseSpace ubiquitin;
  if (ubiquitin_exists) {
    ubiquitin.buildFromFile(ubiquitin_crd_name, CoordinateFileKind::AMBER_INPCRD);
  }
  else {
    missing_files = true;
  }
  check(ubiquitin.getAtomCount(), RelationalOperator::EQUAL, 3105, "The number of atoms in the "
        "ubiquitin system is incorrect.");
  double ub_lx, ub_ly, ub_lz, ub_alpha, ub_beta, ub_gamma;
  extractBoxDimensions(&ub_lx, &ub_ly, &ub_lz, &ub_alpha, &ub_beta, &ub_gamma,
                       ubiquitin.getInverseTransform().data());
  std::vector<double> box_dimensions = { ub_lx, ub_ly, ub_lz, ub_alpha, ub_beta, ub_gamma };
  std::vector<double> box_target(6);
  for (int i = 0; i < 3; i++) {
    box_target[i    ] =  36.2646985;
    box_target[i + 3] = 109.4712190 * pi / 180.0;
  }
  check(box_dimensions, RelationalOperator::EQUAL, Approx(box_target).margin(small),
        "Ubiquitin is solvated in a truncated, octahedral box, but this was not found by "
        "extracting the box dimensions from the inverse transformation matrix.", do_ubiquitin);

  // Read Trp-cage as a CoordinateFrame object
  CoordinateFrame trpcage_coord_frame;
  if (trpcage_exists) {
    trpcage_coord_frame.buildFromFile(trpcage_crd_name, CoordinateFileKind::AMBER_INPCRD);
  }
  check(trpcage_coord_frame.getInterlacedCoordinates(), RelationalOperator::EQUAL,
        trpcage.getInterlacedCoordinates(), "Coordinates read by the CoordinateFrame and "
        "PhaseSpace objects do not agree.", do_trpcage);
  check(trpcage_coord_frame.getUnitCellType() == UnitCellType::NONE, "The Trp-cage system read as "
        "a CoordinateFrame object was not correctly found to have isolated boundary conditions.",
        do_trpcage);

  // Try making a CoordinateFrame based on an existing PhaseSpace object
  CoordinateFrame tip5p_coord_frame(&tip5p);
  CoordinateFrameWriter tip5p_coord_frame_wr = tip5p_coord_frame.data();
  check(tip5p_coord_frame.getAtomCount(), RelationalOperator::EQUAL, tip5p.getAtomCount(),
        "Making an ARRAY-kind CoordinateFrame based on a PhaseSpace object obtains the wrong "
        "number of atoms.", do_tip5p);
  check(tip5p_coord_frame.getUnitCellType() == tip5p.getUnitCellType(), "Making an ARRAY-kind "
        "CoordinateFrame based on a PhaseSpace object obtains the wrong unit cell type.",
        do_tip5p);
 
  // Try making a CoordinateFrameWriter object based on an existing PhaseSpace object
  section(3);
  CoordinateFrameWriter tip5p_access(&tip5p);
  check(tip5p_access.natom, RelationalOperator::EQUAL, tip5p.getAtomCount(), "Making a "
        "CoordinateFrameWriter based on a PhaseSpace object obtains the wrong number of atoms.",
        do_tip5p);
  check(tip5p_access.unit_cell == tip5p.getUnitCellType(), "Making a CoordinateFrameWriter "
        "based on a PhaseSpace object obtains the wrong unit cell type.", do_tip5p);
  const std::vector<double> tip5p_orig_crd = tip5p.getInterlacedCoordinates();
  if (tip5p_exists) {
    for (int i = 0; i < tip5p_access.natom; i++) {
      tip5p_access.xcrd[i] *= 1.7;
      tip5p_access.ycrd[i] *= 0.8;
      tip5p_access.zcrd[i] *= 2.9;
    }
  }
  check(tip5p_coord_frame.getInterlacedCoordinates(), RelationalOperator::EQUAL,
        tip5p_orig_crd, "The ARRAY-kind CoordinateFrame object does not keep an independent set "
        "of coordinates after being created from a PhaseSpace object.", do_tip5p);
  box_target[0] = 18.9668845;
  box_target[1] = 18.7716607;
  box_target[2] = 18.7928724;
  for (int i = 3; i < 6; i++) {
    box_target[i] = 0.5 * pi;
  }
  check(tip5p_coord_frame.getBoxDimensions(), RelationalOperator::EQUAL,
        Approx(box_target).margin(1.0e-7), "Box dimensions stored in a CoordinateFrame object "
        "built from an input file do not meet expectations.", do_tip5p);
  std::vector<double> access_dims(6);
  for (int i = 0; i < 6; i++) {
    access_dims[i] = tip5p_access.boxdim[i];
  }
  check(access_dims, RelationalOperator::EQUAL, Approx(box_target).margin(1.0e-7),
        "Box dimensions obtained from a CoordinateFrameWriter object targeting a PhaseSpace "
        "object do not meet expectations.", do_tip5p);

  // Try reading an entire trajectory of Trp-cage
  const std::vector<int> frames = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  const std::string trpcage_traj_name = base_crd_name + osc + "trpcage.crd";
  const bool traj_exists = (getDrivePathType(trpcage_traj_name) == DrivePathType::FILE);
  if (traj_exists == false) {
    rtWarn("A trajectory of Trp-cage coordinates was not found.  Check ${STORMM_SOURCE} to ensure "
           "that " + trpcage_traj_name + " becomes a valid path.", "test_amber_coordinates");
  }
  const TestPriority do_traj_tests = (traj_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  std::vector<CoordinateFrame> multi_trp = getSelectedFrames(trpcage_traj_name,
                                                             trpcage.getAtomCount(),
                                                             UnitCellType::NONE, frames);
  const std::vector<double> p11_x_answer = {  -3.406,   7.255,   8.085,   6.384,  -0.730,
                                              -1.376,   2.180,  -0.051,   4.558,   4.740 };
  const std::vector<double> p11_y_answer = {   7.351,   5.837,   0.704,   2.707,  -8.253,
                                              -2.722,  -8.488,  -8.305,  -3.565,   7.792 };
  const std::vector<double> p11_z_answer = {   4.455,  -0.182,   1.200,   0.033,  -3.852,
                                               8.617,   0.455,  -3.302,  -4.791,  -0.505 };
  std::vector<double> p11_x(10), p11_y(10), p11_z(10);
  const std::vector<double> mean_x_answer = { -0.2738026316,   0.0413717105,  -0.4480756579,
                                              -0.9473421053,   0.6269210526,   0.4081710526,
                                               0.1894046053,  -0.0134210526,   0.7516776316,
                                               1.3314769737 };
  const std::vector<double> mean_y_answer = { -0.0609835526,  -0.0383782895,   0.6164375000,
                                              -0.7236546053,  -0.5499868421,   0.7145296053,
                                               1.1919638158,   0.2099901316,   0.0292993421,
                                               0.0207138158 };
  const std::vector<double> mean_z_answer = { -0.4562664474,   0.1306907895,   0.3557631579,
                                               0.3126217105,  -0.0445131579,  -0.1440822368,
                                               0.1110000000,  -0.5185756579,   1.3198717105,
                                              -0.8598453947 };
  std::vector<double> mean_x(10), mean_y(10), mean_z(10);
  for (int i = 0; i < 10; i++) {
    CoordinateFrameWriter cfw = multi_trp[i].data();
    p11_x[i] = cfw.xcrd[10];
    p11_y[i] = cfw.ycrd[10];
    p11_z[i] = cfw.zcrd[10];
    mean_x[i] = mean(cfw.xcrd, cfw.natom);
    mean_y[i] = mean(cfw.ycrd, cfw.natom);
    mean_z[i] = mean(cfw.zcrd, cfw.natom);
  }
  check(p11_x, RelationalOperator::EQUAL, p11_x_answer, "Cartesian X coordinates of a series of "
        "trajectory frames taken in isolated boundary conditions do not meet expectations.",
        do_traj_tests);
  check(p11_y, RelationalOperator::EQUAL, p11_y_answer, "Cartesian Y coordinates of a series of "
        "trajectory frames taken in isolated boundary conditions do not meet expectations.",
        do_traj_tests);
  check(p11_z, RelationalOperator::EQUAL, p11_z_answer, "Cartesian Z coordinates of a series of "
        "trajectory frames taken in isolated boundary conditions do not meet expectations.",
        do_traj_tests);
  check(mean_x, RelationalOperator::EQUAL, mean_x_answer, "Cartesian X centers of geometry for a "
        "series of trajectory frames taken in isolated boundary conditions do not meet "
        "expectations.", do_traj_tests);
  check(mean_y, RelationalOperator::EQUAL, mean_y_answer, "Cartesian Y centers of geometry for a "
        "series of trajectory frames taken in isolated boundary conditions do not meet "
        "expectations.", do_traj_tests);
  check(mean_z, RelationalOperator::EQUAL, mean_z_answer, "Cartesian Z centers of geometry for a "
        "series of trajectory frames taken in isolated boundary conditions do not meet "
        "expectations.", do_traj_tests);
  
  // Report missing files
  if (missing_files) {
    rtWarn("One or more coordinate and topology files required by this test program were not "
           "found, causing some tests to be skipped.  Make sure that the $STORMM_SOURCE "
           "environment variable is set to the root directory with src/ and test/ subdirectories, "
           "then re-run the tests to ensure that the libraries are working as intended.",
           "test_phasespace");
  }
  if (snps_exist == false && oe.takeSnapshot() != SnapshotOperation::SNAPSHOT) {
    rtWarn("One or more of the snapshot files required by this test program were not found, "
           "causing some tests to be skipped.  Make sure that the $STORMM_SOURCE environment "
           "variable is set to the root directory with src/ and test/ subdirectories.  If this "
           "error was encountered without also triggering a warning about the actual topology and "
           "coordinate files, the installation may be corrupted.", "test_phasespace");
  }

  // Try using the PhaseSpace object's copy, copy assignment, and move constructors
  section(4);
  std::vector<PhaseSpace> psvec;
  psvec.reserve(4);
  TestPriority do_aggregate_tests;
  if (ubiquitin_exists && tip5p_exists && tip3p_exists && trpcage_exists) {
    psvec.push_back(PhaseSpace(ubiquitin_crd_name, CoordinateFileKind::AMBER_INPCRD));
    psvec.push_back(PhaseSpace(tip5p_crd_name, CoordinateFileKind::AMBER_INPCRD));
    psvec.push_back(PhaseSpace(tip3p_crd_name, CoordinateFileKind::AMBER_INPCRD));
    psvec.push_back(PhaseSpace(trpcage_crd_name, CoordinateFileKind::AMBER_INPCRD));
    do_aggregate_tests = TestPriority::CRITICAL;
  }
  else {
    psvec.push_back(PhaseSpace());
    psvec.push_back(PhaseSpace());
    psvec.push_back(PhaseSpace());
    psvec.push_back(PhaseSpace());
    do_aggregate_tests = TestPriority::ABORT;
  }
  std::vector<double> fifth_atoms(12);
  const std::vector<double> fifth_atom_answer = {  16.1156366,  16.6966148,  28.8721505,
                                                    5.7666007,   6.4934467,   8.8166586,
                                                    6.6229624,   9.9389213,  16.9374141,
                                                   -8.6080000,   3.1350000,  -1.6180000 };
  for (int i = 0; i < 4; i++) {
    std::vector<double> xyz = psvec[i].getInterlacedCoordinates(4, 5);
    for (int j = 0; j < 3; j++) {
      fifth_atoms[(3 * i) + j] = xyz[j];
    }
  }
  check(fifth_atoms, RelationalOperator::EQUAL, Approx(fifth_atom_answer).margin(1.0e-6),
        "Coordinates of atoms from a Standard Template Library vector of PhaseSpace objects "
        "created with the push_back method are incorrect.");
  psvec.resize(0);
  for (int i = 0; i < 4; i++) {
    psvec.push_back(createRandomParticles(400 - (4 * i)));
  }
  const std::vector<double> set_sums_ans = {   37.2396341144,  -68.2296975604,   26.1908040068,
                                              -92.2686284816, -121.8449964662,  -16.1506079272,
                                               29.2925057928,   29.3206610555, -139.9787246260,
                                               -0.8906798012, -180.4192460200,  -25.0029609269 };
  std::vector<double> set_sums(12);
  for (int i = 0; i < 4; i++) {
    double xsum = 0.0;
    double ysum = 0.0;
    double zsum = 0.0;
    if (i < 2) {

      // Test the copy assignment constructor
      PhaseSpace pscopy;
      pscopy = psvec[i];
      const PhaseSpaceWriter psw = pscopy.data();
      for (int j = 0; j < psw.natom; j++) {
        xsum += psw.xcrd[j];
        ysum += psw.ycrd[j];
        zsum += psw.zcrd[j];
      }
    }
    else {

      // Test the copy constructor
      const PhaseSpace pscopy(psvec[i]);
      const PhaseSpaceReader psr = pscopy.data();
      for (int j = 0; j < psr.natom; j++) {
        xsum += psr.xcrd[j];
        ysum += psr.ycrd[j];
        zsum += psr.zcrd[j];
      }      
    }
    set_sums[(3 * i)    ] = xsum;
    set_sums[(3 * i) + 1] = ysum;
    set_sums[(3 * i) + 2] = zsum;
  }
  check(set_sums, RelationalOperator::EQUAL, Approx(set_sums_ans).margin(1.0e-8), "Sums of "
        "random particle coordinates obtained from PhaseSpace objects after std::vector "
        "push_back() and subsequent copy assignment operations do not meet expectations.");

  // Create a CoordinateSeries based on several copies of a coordinate frame
  const std::string stereo_crd_name = base_crd_name + osc + "stereo_L1.inpcrd";
  const std::string stereo_trj_name = base_crd_name + osc + "stereo_L1x.crd";
  const bool stereo_exists = (getDrivePathType(stereo_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(stereo_trj_name) == DrivePathType::FILE);
  const TestPriority do_stereo = (stereo_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (stereo_exists == false) {
    rtWarn("Coordinate files for a highly stereo-isomerized ligand, " + stereo_crd_name + " and " +
           stereo_trj_name + " were not found.  Check the ${STORMM_SOURCE} variable, currently "
           "set to " + oe.getStormmSourcePath() + " to make sure it is valid.  Subsequent tests "
           "will be skipped.", "test_amber_coordinates");
  }
  const CoordinateFrame stro_cf = (stereo_exists) ? CoordinateFrame(stereo_crd_name) :
                                                    CoordinateFrame();
  const CoordinateFrameReader stro_cfr = stro_cf.data();
  CoordinateSeries<double> stro_cs = (stereo_exists) ?
                                     CoordinateSeries<double>(stereo_trj_name,
                                                              stro_cf.getAtomCount()) :
                                     CoordinateSeries<double>();
  std::vector<double> xyz_sum(3, 0.0);
  for (int i = 0; i < stro_cs.getFrameCount(); i++) {
    std::vector<double> frm_crd = stro_cs.getInterlacedCoordinates<double>(i, 0, 3);
    for (int j = 0; j < 3; j++) {
      for (int k = 0; k < 3; k++) {
        xyz_sum[k] += frm_crd[(3 * j) + k];
      }
    }
  }
  const std::vector<double> xyz_sum_ans = { 346.572, -149.484, 2220.577 };
  check(xyz_sum, RelationalOperator::EQUAL, Approx(xyz_sum_ans).margin(1.0e-8), "Sums of "
        "coordinates extracted from a CoordinateSeries object do not meet expectations.",
        do_stereo);
  check(stro_cs.getAtomCount(), RelationalOperator::EQUAL, stro_cf.getAtomCount(),
        "The trajectory read from " + stereo_trj_name + " does not contain the correct number of "
        "atoms.", do_stereo);
  check(stro_cs.getFrameCount(), RelationalOperator::EQUAL, 20, "The trajectory read from " +
        stereo_trj_name + " does not contain the correct number of frames.", do_stereo);
  std::vector<double> fr2_crd = stro_cs.getInterlacedCoordinates<double>(2);
  for (int i = 0; i < 3; i++) {
    xyz_sum[i] = 0.0;
  }
  for (int i = 0; i < stro_cfr.natom; i++) {
    xyz_sum[0] += fr2_crd[(3 * i)    ];
    xyz_sum[1] += fr2_crd[(3 * i) + 1];
    xyz_sum[2] += fr2_crd[(3 * i) + 2];
  }
  const std::vector<double> fr2_sum_ans = { 482.8030, -198.3540, 2850.508 };
  check(xyz_sum, RelationalOperator::EQUAL, Approx(fr2_sum_ans).margin(1.0e-8), "Sums of "
        "coordinates extracted from one frame of a CoordinateSeries object do not meet "
        "expectations.", do_stereo);
  stro_cs.pushBack(stro_cf);
  check(stro_cs.getFrameCount(), RelationalOperator::EQUAL, 21, "The pushBack() member function "
        "of the CoordinateSeries object does not extend the series as expected.", do_stereo);
  check(stro_cf.getInterlacedCoordinates(), RelationalOperator::EQUAL,
        stro_cs.getInterlacedCoordinates<double>(20), "The pushBack() member function of the "
        "CoordinateSeries object does not extend the series as expected.", do_stereo);
  if (stereo_exists) {
    stro_cs.importFromFile(stereo_trj_name);
  }
  check(stro_cs.getFrameCount(), RelationalOperator::EQUAL, 41, "Importing a second time from a "
        "trajectory file did not extend a CoordinateSeries object in the expected manner.",
        do_stereo);
  check(stro_cs.getInterlacedCoordinates<double>(4), RelationalOperator::EQUAL,
        stro_cs.getInterlacedCoordinates<double>(25), "Correspondence between frames read from "
        "the same trajectory file is not maintained in a CoordinateSeries object.", do_stereo);
  for (int i = 0; i < 5; i++) {
    stro_cs.pushBack(stro_cf);
  }
  stro_cs.shrinkToFit();
  check(stro_cs.getInterlacedCoordinates<double>(20), RelationalOperator::EQUAL,
        stro_cs.getInterlacedCoordinates<double>(43), "Frames that should match exactly no longer "
        "do after applying the CoordinateSeries object's shirnkToFit method.", do_stereo);
  CoordinateSeries<int> istro_cs = changeCoordinateSeriesType<double, int>(stro_cs, 14);
  
  const CoordinateFrame frame17 = istro_cs.exportFrame(17);
  const std::vector<double> stro_drep = stro_cs.getInterlacedCoordinates<double>(17);
  const std::vector<double> stro_di14rep = frame17.getInterlacedCoordinates();
  check(maxAbsoluteDifference(stro_di14rep, stro_drep), RelationalOperator::LESS_THAN, 1.0e-4,
        "A 14 bit fixed-precision representation of the coordinate series does not reproduce the "
        "original, double-precision data with the expected accuracy.", do_stereo);
  CoordinateSeries<int> tip5p_cs(tip5p, 15, 16);
  const std::vector<double> tip5p_drep = tip5p.getInterlacedCoordinates();
  const CoordinateFrame frame9 = tip5p_cs.exportFrame(9);
  const std::vector<double> tip5p_di16rep = frame9.getInterlacedCoordinates();
  check(maxAbsoluteDifference(tip5p_di16rep, tip5p_drep), RelationalOperator::LESS_THAN, 1.5e-5,
        "A 16 bit fixed-precision representation of the coordinate series does not reproduce the "
        "original, double-precision data with the expected accuracy.", do_tip5p);
  check(tip5p_cs.getUnitCellType() == tip5p.getUnitCellType(), "A CoordinateSeries object does "
        "not adopt the original PhaseSpace object's unit cell type (" +
        getEnumerationName(tip5p.getUnitCellType()) + ").", do_tip5p);
  CHECK_THROWS_SOFT(tip5p_cs.pushBack(trpcage), "A CoordinateSeries accepts new coordinates from "
                    "a system with the wrong particle number.", do_tip5p);
  CHECK_THROWS_SOFT(tip5p_cs.importFromFile(tip3p_crd_name), "A CoordinateSeries accepts new "
                    "coordinates from an Amber input coordinates file with the wrong particle "
                    "number.", do_tip3p);

  // Check coordinate series construction from alternative template types
  const CoordinateSeries<llint> tip5p_cs_64(tip5p_cs, 24);
  const CoordinateFrame frame9_int64 = tip5p_cs_64.exportFrame(9);
  check(tip5p_cs_64.getInterlacedCoordinates<double>(9), RelationalOperator::EQUAL,
        Approx(tip5p_cs.getInterlacedCoordinates<double>(9)).margin(stormm::constants::tiny),
        "A CoordinateSeries constructed from another of a lower fixed precision representation "
        "does not preserve the correct information.", do_tip5p);
  Xoroshiro128pGenerator xrs_crd(1984232);
  CoordinateSeries<double> tip5p_cs_dbl(tip5p_cs);
  CoordinateSeriesWriter<double> tip5p_csw = tip5p_cs_dbl.data();
  const size_t padded_natom_zu = roundUp<size_t>(tip5p_csw.natom, stormm::constants::warp_size_zu);
  for (size_t i = 0; i < tip5p_csw.nframe; i++) {
    for (size_t j = 0; j < tip5p_csw.natom; j++) {
      const size_t atom_idx = (padded_natom_zu * i) + j;
      tip5p_csw.xcrd[atom_idx] += 0.05 * xrs_crd.gaussianRandomNumber();
      tip5p_csw.ycrd[atom_idx] += 0.05 * xrs_crd.gaussianRandomNumber();
      tip5p_csw.zcrd[atom_idx] += 0.05 * xrs_crd.gaussianRandomNumber();
    }
  }
  CoordinateSeries<int> tip5p_cs_rng(tip5p_cs_dbl, 16);
  std::vector<int> perturbation(3 * tip5p.getAtomCount());
  const std::vector<int> orig_tip5p_crd = tip5p_cs.getInterlacedCoordinates<int>(5);
  const std::vector<int> pert_tip5p_crd = tip5p_cs_rng.getInterlacedCoordinates<int>(5);
  const std::vector<int> fp32_tip5p_crd = tip5p_cs_dbl.getInterlacedCoordinates<int>(5, 16);
  for (int i = 0; i < perturbation.size(); i++) {
    perturbation[i] = pert_tip5p_crd[i] - orig_tip5p_crd[i];
  }
  snapshot(tip5p_vel_snapshot, polyNumericVector(perturbation), "t5p_pert_xyz", 1.0e-4,
           "Perturbations to TIP5P water positions do not meet expectations.", oe.takeSnapshot(),
           1.0e-8, NumberFormat::INTEGER, PrintSituation::APPEND, snap_check);
  check(pert_tip5p_crd, RelationalOperator::EQUAL, fp32_tip5p_crd, "Extracting coordinates from a "
        "CoordinateSeries as fixed precision integers does not yield the same result as "
        "converting the series to a different object in the same fixed precision format.",
        do_tip5p);
  const std::vector<int> highres_orig_tip5p_crd = tip5p_cs.getInterlacedCoordinates<int>(5, 18);
  std::vector<int> mult_tip5p_crd(orig_tip5p_crd);
  elementwiseMultiply(&mult_tip5p_crd, 4);
  check(highres_orig_tip5p_crd, RelationalOperator::EQUAL, mult_tip5p_crd, "Conversion of data in "
        "a CoordinateSeries with fixed-precision representation does not go to a deeper "
        "fixed-precision representation with the expected result.", do_tip5p);
  const std::vector<int> lowres_orig_tip5p_crd = tip5p_cs.getInterlacedCoordinates<int>(5, 14);
  elementwiseDivide(&mult_tip5p_crd, 16);
  check(lowres_orig_tip5p_crd, RelationalOperator::EQUAL, Approx(mult_tip5p_crd).margin(0.01),
        "Conversion of data in a CoordinateSeries with fixed-precision representation does not go "
        "to a shallower fixed-precision representation with the expected result.", do_tip5p);

  // Summary evaluation
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
