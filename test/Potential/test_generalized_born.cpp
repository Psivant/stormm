#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/static_exclusionmask.h"
#include "../../src/Potential/nonbonded_potential.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::llint;
using stormm::llint_type_index;
#ifndef STORMM_USE_HPC
using stormm::double2;
using stormm::int2;
using stormm::int3;
#endif
using stormm::constants::ExceptionResponse;
using stormm::data_types::getStormmScalarTypeName;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::energy::StateVariable;
using stormm::energy::StaticExclusionMask;
using stormm::errors::rtWarn;
using stormm::parse::char4ToString;
using stormm::parse::NumberFormat;
using stormm::parse::polyNumericVector;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using stormm::symbols::amber_ancient_bioq;
using stormm::topology::AtomGraph;
using stormm::topology::AtomicRadiusSet;
using stormm::topology::ImplicitSolventModel;
using stormm::topology::ImplicitSolventKit;
using stormm::topology::NonbondedKit;
using stormm::trajectory::CoordinateFileKind;
using stormm::trajectory::PhaseSpace;
using stormm::trajectory::PhaseSpaceWriter;
using stormm::trajectory::TrajectoryKind;
using namespace stormm::energy;
using namespace stormm::testing;
using namespace stormm::generalized_born_defaults;

//-------------------------------------------------------------------------------------------------
// Compute forces on selected atoms using a finite difference scheme.
//
// Arguments:
//   ag:               System topology
//   se:               Static exclusion mask (the exclusions apply only to the non-bonded
//                     electrostatic and van-der Waals interactions, but this object also holds
//                     information on the tile sizes to help organize the loop structure)
//   ps:               System coordinates and forces
//   ngb_tab:          Tables of constants for "Neck" GB
//   sample_interval:  Sampling interval for computing finite-difference forces (the calculation
//                     scales as 4 * (N / sampling interval) * (1/2) N^2, so for larger systems not
//                     every atom's forces can be evaluated in a cost-effective test!)
//-------------------------------------------------------------------------------------------------
std::vector<double> forceByFiniteDifference(const AtomGraph &ag, const StaticExclusionMask &se,
                                            PhaseSpace *ps,
                                            const NeckGeneralizedBornTable &ngb_tab,
                                            const int sample_interval) {
  ScoreCard sc(1, 1, 32);
  PhaseSpaceWriter psw = ps->data();
  const EvaluateForce evfrc = EvaluateForce::YES;

  // Perturb atomic positions by an amount that the double-precision coordinate arrays will be
  // able to represent exactly.
  const double disc = 0.000003814697265625;

  // Get the baseline energy, then do perturbations in X, Y, and Z for selected atoms
  std::vector<double> result;
  for (int i = 0; i < psw.natom; i += sample_interval) {
    const double midpt = evaluateGeneralizedBornEnergy(ag, se, ngb_tab, ps, &sc, evfrc, 0);
    psw.xcrd[i] += disc;
    const double plusx = evaluateGeneralizedBornEnergy(ag, se, ngb_tab, ps, &sc, evfrc, 0);
    psw.xcrd[i] -= disc;
    psw.ycrd[i] += disc;
    const double plusy = evaluateGeneralizedBornEnergy(ag, se, ngb_tab, ps, &sc, evfrc, 0);
    psw.ycrd[i] -= disc;
    psw.zcrd[i] += disc;
    const double plusz = evaluateGeneralizedBornEnergy(ag, se, ngb_tab, ps, &sc, evfrc, 0);
    psw.zcrd[i] -= disc;

    // Gradient is rise over run, but force is the negative of gradient: if the particle moves
    // positively in X and the energy rises, then the force will be to push the particle backwards,
    // in the negative X direction, to move back down the energy gradient.
    result.push_back(-(plusx - midpt) / disc);
    result.push_back(-(plusy - midpt) / disc);
    result.push_back(-(plusz - midpt) / disc);
  }

  // Wipe the force arrays clean before returning the outcome
  ps->initializeForces();
  return result;
}

//-------------------------------------------------------------------------------------------------
// Compute the non-bonded forces on a system using the precision model specified by a particular
// parameter set, coordinate, and force representation.
//
// Arguments:
//   nbk:         Non-bonded parameters in a particular precision
//   semask:      Non-bonded exclusions list
//   ps:          Coordinates of the system in double precision
//   gb_ref_frc:  Reference Generalized Born forces on all particles
//   tol:         Tolerance to apply to each force component on every particle when comparing to
//                the reference
//   do_tests:    Indicator of whether to pursue tests
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
void testNBPrecisionModel(const NonbondedKit<Tcalc> nbk, const ImplicitSolventKit<Tcalc> isk,
                          const NeckGeneralizedBornKit<Tcalc> ngb_kit,
                          const StaticExclusionMask &semask, const PhaseSpace &ps,
                          const std::vector<double> &gb_ref_frc, const double tol,
                          const TestPriority do_tests) {
  const size_t tcoord_ct = std::type_index(typeid(Tcoord)).hash_code();
  const size_t tforce_ct = std::type_index(typeid(Tforce)).hash_code();
  const CoordinateSeries<Tcoord> crd(ps, 1, 26 * (tcoord_ct == llint_type_index));
  CoordinateSeries<Tforce> frc(ps, 1, 23 * (tforce_ct == llint_type_index));
  CoordinateSeries<Tforce> xtra(ps, 1, 23 * (tforce_ct == llint_type_index));
  CoordinateSeriesReader crdr = crd.data();
  CoordinateSeriesWriter frcw = frc.data();
  CoordinateSeriesWriter xtrw = xtra.data();
  ScoreCard sc(1, 16, 32);
  const Tforce zero = 0.0;
  for (int i = 0; i < frcw.natom; i++) {
    frcw.xcrd[i] = zero;
    frcw.ycrd[i] = zero;
    frcw.zcrd[i] = zero;
  }
  evaluateGeneralizedBornEnergy<Tcoord, Tforce, Tcalc>(nbk, semask.data(), isk, ngb_kit, crdr.xcrd,
                                                       crdr.ycrd, crdr.zcrd, frcw.xcrd, frcw.ycrd,
                                                       frcw.zcrd, xtrw.xcrd, xtrw.ycrd, xtrw.zcrd,
                                                       &sc, EvaluateForce::YES, 0,
                                                       crdr.inv_gpos_scale, frcw.gpos_scale);
  const CoordinateFrame gb_frame = frc.exportFrame(0);
  const std::vector<double> gb_result = gb_frame.getInterlacedCoordinates();
  check(gb_result, RelationalOperator::EQUAL, Approx(gb_ref_frc).margin(tol),
        "Generalized Born forces do not agree with the reference when computed in " +
        getStormmScalarTypeName<Tcalc>() + " with " + getStormmScalarTypeName<Tcoord>() +
        " coordinates and " + getStormmScalarTypeName<Tforce>() + " force accumulation.",
        do_tests);
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
  StopWatch timer;
  const int input_timings =  timer.addCategory("Input file parsing");
  const int gb_timings    =  timer.addCategory("Generalized born");
  
  // Section 1
  section("ImplicitSolventRecipe abstract");

  // Section 2
  section("HCT GB");

  // Section 3
  section("OBC GB");

  // Section 4
  section("Neck GB");

  // Section 5
  section("Coordinate representations");
  
  // Locate topologies and coordinate files
  const char osc = osSeparator();
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string base_ptl_name = oe.getStormmSourcePath() + osc + "test" + osc + "Potential";
  const std::string trpi_top_name = base_top_name + osc + "trpcage.top";
  const std::string trpi_crd_name = base_crd_name + osc + "trpcage.inpcrd";
  const std::string dhfr_top_name = base_top_name + osc + "dhfr_cmap.top";
  const std::string dhfr_crd_name = base_crd_name + osc + "dhfr_cmap.inpcrd";
  const std::string alad_top_name = base_top_name + osc + "ala_dipeptide.top";
  const std::string alad_crd_name = base_crd_name + osc + "ala_dipeptide.inpcrd";
  const bool systems_exist = (getDrivePathType(trpi_top_name) == DrivePathType::FILE &&
                              getDrivePathType(trpi_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(dhfr_top_name) == DrivePathType::FILE &&
                              getDrivePathType(dhfr_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(alad_top_name) == DrivePathType::FILE &&
                              getDrivePathType(alad_crd_name) == DrivePathType::FILE);
  const TestPriority do_tests = (systems_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (systems_exist == false) {
    rtWarn("Files for the Trp-cage miniprotein, DHFR globular protein, and alanine dipeptide were "
           "not found.  These files should be found in the ${STORMM_SOURCE}/test/Topology and "
           "${STORMM_SOURCE}/test/Trajectory directories.  Check the $STORMM_SOURCE environment "
           "variable.  A number of tests will be skipped.", "test_generalized_born");
  }
  const std::string trpi_snapshot(base_ptl_name + osc + "trpcage_gb_forces.m");
  const std::string dhfr_snapshot(base_ptl_name + osc + "dhfr_gb_forces.m");
  const std::string alad_snapshot(base_ptl_name + osc + "ala_dipeptide_gb_forces.m");
  const bool snaps_exist = (getDrivePathType(trpi_snapshot) == DrivePathType::FILE &&
                            getDrivePathType(dhfr_snapshot) == DrivePathType::FILE &&
                            getDrivePathType(alad_snapshot) == DrivePathType::FILE);
  const TestPriority do_snaps = (snaps_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (snaps_exist == false && oe.takeSnapshot() != SnapshotOperation::SNAPSHOT) {
    rtWarn("Snapshot files for the Trp-cage miniprotein, DHFR globular protein, and alanine "
           "dipeptide were not found.  These should be found in the "
           "${STORMM_SOURCE}/test/Potential directory.  Some tests will be skipped.",
           "test_generalized_born");
  }

  // Read topologies and coordinates.  Suppress charge discretization and normalization to
  // maintain the best possible agreement with results from Amber's sander program.
  AtomGraph trpi_ag, dhfr_ag, alad_ag;
  PhaseSpace trpi_ps, dhfr_ps, alad_ps;
  const double charge_discretization_inc = 9.31322574615478515625E-10;
  const double charge_rounding_tol = 1.0e-10;
  const double amber_default_qq14 = 1.2;
  const double amber_default_lj14 = 2.0;
  const double charmm_default_qq14 = 1.0;
  const double charmm_default_lj14 = 1.0;
  const double sander_dielectric = 78.5;
  timer.assignTime(0);
  if (systems_exist) {
    trpi_ag.buildFromPrmtop(trpi_top_name, ExceptionResponse::SILENT, amber_ancient_bioq,
                            amber_default_qq14, amber_default_lj14, charge_rounding_tol,
                            charge_discretization_inc);
    trpi_ps.buildFromFile(trpi_crd_name, CoordinateFileKind::AMBER_INPCRD);
    dhfr_ag.buildFromPrmtop(dhfr_top_name, ExceptionResponse::SILENT, amber_ancient_bioq,
                            charmm_default_qq14, charmm_default_lj14, charge_rounding_tol,
                            charge_discretization_inc);
    dhfr_ps.buildFromFile(dhfr_crd_name, CoordinateFileKind::AMBER_INPCRD);
    alad_ag.buildFromPrmtop(alad_top_name, ExceptionResponse::SILENT, amber_ancient_bioq,
                            amber_default_qq14, amber_default_lj14, charge_rounding_tol,
                            charge_discretization_inc);
    alad_ps.buildFromFile(alad_crd_name, CoordinateFileKind::AMBER_INPCRD);
  }
  timer.assignTime(input_timings);
  const StaticExclusionMask trpi_se = (systems_exist) ? StaticExclusionMask(&trpi_ag) :
                                                        StaticExclusionMask();
  const StaticExclusionMask dhfr_se = (systems_exist) ? StaticExclusionMask(&dhfr_ag) :
                                                        StaticExclusionMask();
  const StaticExclusionMask alad_se = (systems_exist) ? StaticExclusionMask(&alad_ag) :
                                                        StaticExclusionMask();

  // Apply an implicit solvent model to each topology and check the result
  section(1);
  std::vector<AtomGraph*> all_topologies = { &trpi_ag, &dhfr_ag, &alad_ag };
  std::vector<PhaseSpace*> all_coords = { &trpi_ps, &dhfr_ps, &alad_ps };
  std::vector<StaticExclusionMask> all_semasks = { trpi_se, dhfr_se, alad_se };
  const size_t system_count = all_topologies.size();
  for (size_t i = 0; i < system_count; i++) {
    all_topologies[i]->setImplicitSolventModel(ImplicitSolventModel::HCT_GB, sander_dielectric);
  }
  ImplicitSolventKit<double> trpi_isk = all_topologies[0]->getDoublePrecisionImplicitSolventKit();
  ImplicitSolventKit<double> dhfr_isk = all_topologies[1]->getDoublePrecisionImplicitSolventKit();
  check(trpi_isk.dielectric, RelationalOperator::EQUAL, sander_dielectric,
        "The implicit solvent dielectric was not properly set.", do_tests);
  check(trpi_isk.saltcon, RelationalOperator::EQUAL, 0.0, "The implicit solvent salt "
        "concentration was not properly set.", do_tests);
  check(dhfr_isk.igb == ImplicitSolventModel::HCT_GB, "The Hawkins / Cramer / Truhlar GB model "
        "was not imparted to systems as expected.", do_tests);
  check(sum<double>(trpi_isk.pb_radii, trpi_isk.natom), RelationalOperator::EQUAL,
        Approx(444.95).margin(1.0e-5), "PB radii for the Trp-cage system were not conveyed to the "
        "ImplicitSolventKit abstract as expected.", do_tests);
  check(sum<double>(dhfr_isk.gb_screen, dhfr_isk.natom), RelationalOperator::EQUAL,
        Approx(1998.81).margin(1.0e-5), "GB screening factors for the DHFR system were not "
        "conveyed to the ImplicitSolventKit abstract as expected.", do_tests);

  // Compute Generalized Born energy and forces.  The DHFR system departs most significantly (0.05
  // kcal/mol) from the GB energies of the equivalent calculation in Amber's sander facility, due
  // to the relatively large number of charge increments (75) that the default topology settings
  // distribute over its atoms.  Other systems, which distribute 6 (Trp-cage) and 3 (alanine
  // dipeptide) charge increments of ~1.19e-7 across their most charged atoms, agree to within
  // 0.0002 kcal/mol of the sander results.  All systems can be brought into near perfect agreement
  // with the sander calculations by disabling charge discretization and normalization.
  const NeckGeneralizedBornTable ngb_tab;
  ScoreCard all_systems_sc(3, 1, 32);
  const int trpi_idx = 0;
  const int dhfr_idx = 1;
  const int alad_idx = 2;
  const EvaluateForce evfrc = EvaluateForce::YES;
  const std::vector<double> hct_gb_answer   = { -313.0368696, -2480.8487515, -14.7193347 };
  const std::vector<double> obc_gb_answer   = { -320.6017219, -2504.3368641, -16.0329696 };
  const std::vector<double> obc2_gb_answer  = { -298.2900710, -2401.7427251, -14.9957762 };
  const std::vector<double> neck_gb_answer  = { -311.6397166, -2503.8809649, -14.8329479 };
  const std::vector<double> neck2_gb_answer = { -328.6694773, -2686.2812806, -14.8989435 };

  // Check energy and force computations for each successive GB model
  section(2);
  std::vector<double> gb_energy(system_count);
  std::vector<std::vector<double>> gb_forces(system_count);
  const TrajectoryKind tkind = TrajectoryKind::FORCES;
  for (size_t i = 0; i < system_count; i++) {
    all_topologies[i]->setImplicitSolventModel(ImplicitSolventModel::HCT_GB, 78.5);
    all_coords[i]->initializeForces();
    timer.assignTime(0);
    gb_energy[i] = evaluateGeneralizedBornEnergy(*(all_topologies[i]), all_semasks[i], ngb_tab,
                                                 all_coords[i], &all_systems_sc, evfrc, i);
    timer.assignTime(gb_timings);
    gb_forces[i] = all_coords[i]->getInterlacedCoordinates(tkind);
  }
  check(all_systems_sc.reportInstantaneousStates(StateVariable::GENERALIZED_BORN),
        RelationalOperator::EQUAL, Approx(hct_gb_answer).margin(1.0e-4), "Hawkins / Cramer / "
        "Truhlar Generalized Born energies do not meet expectations.", do_tests);
  snapshot(trpi_snapshot, polyNumericVector(gb_forces[0]), "trpcage_hct_gb_frc",
           NumberFormat::STANDARD_REAL, "Forces due to Hawkins / Cramer / Truhlar Generalized "
           "Born interactions in the Trp-cage system do not meet expectations.", oe.takeSnapshot(),
           1.0e-6, 1.0e-12, PrintSituation::OVERWRITE, do_snaps);
  snapshot(dhfr_snapshot, polyNumericVector(gb_forces[1]), "dhfr_hct_gb_frc",
           NumberFormat::STANDARD_REAL, "Forces due to Hawkins / Cramer / Truhlar Generalized "
           "Born interactions in the DHFR system do not meet expectations.", oe.takeSnapshot(),
           1.0e-6, 1.0e-12, PrintSituation::OVERWRITE, do_snaps);
  snapshot(alad_snapshot, polyNumericVector(gb_forces[2]), "ala_hct_gb_frc",
           NumberFormat::STANDARD_REAL, "Forces due to Hawkins / Cramer / Truhlar Generalized "
           "Born interactions in the alanine dipeptide system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::OVERWRITE, do_snaps);  
  section(3);
  for (size_t i = 0; i < system_count; i++) {
    all_topologies[i]->setImplicitSolventModel(ImplicitSolventModel::OBC_GB, 78.5);
    all_coords[i]->initializeForces();
    timer.assignTime(0);
    gb_energy[i] = evaluateGeneralizedBornEnergy(*(all_topologies[i]), all_semasks[i], ngb_tab,
                                                 all_coords[i], &all_systems_sc, evfrc, i);
    timer.assignTime(gb_timings);
    gb_forces[i] = all_coords[i]->getInterlacedCoordinates(tkind);
  }
  check(all_systems_sc.reportInstantaneousStates(StateVariable::GENERALIZED_BORN),
        RelationalOperator::EQUAL, Approx(obc_gb_answer).margin(1.0e-4), "Onufriev / Bashford / "
        "Case (model I) Generalized Born energies do not meet expectations.", do_tests);
  snapshot(trpi_snapshot, polyNumericVector(gb_forces[0]), "trpcage_obc_gb_frc",
           NumberFormat::STANDARD_REAL, "Forces due to Onufriev / Bashford / Case (model I) "
           "Generalized Born interactions in the Trp-cage system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::APPEND, do_snaps);
  snapshot(dhfr_snapshot, polyNumericVector(gb_forces[1]), "dhfr_obc_gb_frc",
           NumberFormat::STANDARD_REAL, "Forces due to Onufriev / Bashford / Case (model I) "
           "Generalized Born interactions in the DHFR system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::APPEND, do_snaps);
  snapshot(alad_snapshot, polyNumericVector(gb_forces[2]), "ala_obc_gb_frc",
           NumberFormat::STANDARD_REAL, "Forces due to Onufriev / Bashford / Case (model I) "
           "Generalized Born interactions in the alanine dipeptide system do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::APPEND, do_snaps);
  for (size_t i = 0; i < system_count; i++) {
    all_topologies[i]->setImplicitSolventModel(ImplicitSolventModel::OBC_GB_II, 78.5);
    all_coords[i]->initializeForces();
    timer.assignTime(0);
    gb_energy[i] = evaluateGeneralizedBornEnergy(*(all_topologies[i]), all_semasks[i], ngb_tab,
                                                 all_coords[i], &all_systems_sc, evfrc, i);
    timer.assignTime(gb_timings);
    gb_forces[i] = all_coords[i]->getInterlacedCoordinates(tkind);
  }
  check(all_systems_sc.reportInstantaneousStates(StateVariable::GENERALIZED_BORN),
        RelationalOperator::EQUAL, Approx(obc2_gb_answer).margin(1.0e-4), "Onufriev / Bashford / "
        "Case (model II) Generalized Born energies do not meet expectations.", do_tests);
  snapshot(trpi_snapshot, polyNumericVector(gb_forces[0]), "trpcage_obc2_gb_frc",
           NumberFormat::STANDARD_REAL, "Forces due to Onufriev / Bashford / Case (model II) "
           "Generalized Born interactions in the Trp-cage system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::APPEND, do_snaps);
  snapshot(dhfr_snapshot, polyNumericVector(gb_forces[1]), "dhfr_obc2_gb_frc",
           NumberFormat::STANDARD_REAL, "Forces due to Onufriev / Bashford / Case (model II) "
           "Generalized Born interactions in the DHFR system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::APPEND, do_snaps);
  snapshot(alad_snapshot, polyNumericVector(gb_forces[2]), "ala_obc2_gb_frc",
           NumberFormat::STANDARD_REAL, "Forces due to Onufriev / Bashford / Case (model II) "
           "Generalized Born interactions in the alanine dipeptide system do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::APPEND, do_snaps);

  // Some topologies are not acceptable for "neck" GB calculations (they must have Bondi radii).
  // This is checked in the topology and radii are reassigned (based on bondi PB radii) if
  // deficiencies are found.  Note that this permanently modifies the PB radii in the topology,
  // which cannot be undone by simply reverting to a different implicit solvent model that is
  // agnostic to the radii.  Do this to all topologies except alanine dipeptide, which already has
  // acceptable (albeit not Bondi) radii for the "neck" GB calculations.
  section(4);
  for (size_t i = 0; i < system_count; i++) {
    const AtomicRadiusSet rset_apply = (all_topologies[i]->getFileName() == alad_top_name) ?
                                       AtomicRadiusSet::NONE : AtomicRadiusSet::BONDI;
    all_topologies[i]->setImplicitSolventModel(ImplicitSolventModel::NECK_GB, 78.5, 0.0,
                                               rset_apply);
    all_coords[i]->initializeForces();
    timer.assignTime(0);
    gb_energy[i] = evaluateGeneralizedBornEnergy(*(all_topologies[i]), all_semasks[i], ngb_tab,
                                                 all_coords[i], &all_systems_sc, evfrc, i);
    timer.assignTime(gb_timings);
    gb_forces[i] = all_coords[i]->getInterlacedCoordinates(tkind);
  }
  check(all_systems_sc.reportInstantaneousStates(StateVariable::GENERALIZED_BORN),
        RelationalOperator::EQUAL, Approx(neck_gb_answer).margin(1.0e-4), "Mongan \"Neck\" "
        "(model I) Generalized Born energies do not meet expectations.", do_tests);
  snapshot(trpi_snapshot, polyNumericVector(gb_forces[0]), "trpcage_neck_gb_frc",
           NumberFormat::STANDARD_REAL, "Forces due to Mongan \"Neck\" (model I) Generalized Born "
           "interactions in the Trp-cage system do not meet expectations.", oe.takeSnapshot(),
           1.0e-6, 1.0e-12, PrintSituation::APPEND, do_snaps);
  snapshot(dhfr_snapshot, polyNumericVector(gb_forces[1]), "dhfr_neck_gb_frc",
           NumberFormat::STANDARD_REAL, "Forces due to Mongan \"Neck\" (model I) Generalized Born "
           "interactions in the DHFR system do not meet expectations.", oe.takeSnapshot(), 1.0e-6,
           1.0e-12, PrintSituation::APPEND, do_snaps);
  snapshot(alad_snapshot, polyNumericVector(gb_forces[2]), "ala_neck_gb_frc",
           NumberFormat::STANDARD_REAL, "Forces due to Mongan \"Neck\" (model I) Generalized Born "
           "interactions in the alanine dipeptide system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::APPEND, do_snaps);
  for (size_t i = 0; i < system_count; i++) {
    const AtomicRadiusSet rset_apply = (all_topologies[i]->getFileName() == alad_top_name) ?
                                       AtomicRadiusSet::NONE : AtomicRadiusSet::MBONDI3;
    all_topologies[i]->setImplicitSolventModel(ImplicitSolventModel::NECK_GB_II, 78.5, 0.0,
                                               rset_apply);
    all_coords[i]->initializeForces();
    timer.assignTime(0);
    gb_energy[i] = evaluateGeneralizedBornEnergy(*(all_topologies[i]), all_semasks[i], ngb_tab,
                                                 all_coords[i], &all_systems_sc, evfrc, i);
    timer.assignTime(gb_timings);
    gb_forces[i] = all_coords[i]->getInterlacedCoordinates(tkind);
  }
  check(all_systems_sc.reportInstantaneousStates(StateVariable::GENERALIZED_BORN),
        RelationalOperator::EQUAL, Approx(neck2_gb_answer).margin(1.0e-4), "Mongan \"Neck\" "
        "(model II) Generalized Born energies do not meet expectations.", do_tests);
  snapshot(trpi_snapshot, polyNumericVector(gb_forces[0]), "trpcage_neck2_gb_frc",
           NumberFormat::STANDARD_REAL, "Forces due to Mongan \"Neck\" (model II) Generalized "
           "Born interactions in the Trp-cage system do not meet expectations.", oe.takeSnapshot(),
           1.0e-6, 1.0e-12, PrintSituation::APPEND, do_snaps);
  snapshot(dhfr_snapshot, polyNumericVector(gb_forces[1]), "dhfr_neck2_gb_frc",
           NumberFormat::STANDARD_REAL, "Forces due to Mongan \"Neck\" (model II) Generalized "
           "Born interactions in the DHFR system do not meet expectations.", oe.takeSnapshot(),
           1.0e-6, 1.0e-12, PrintSituation::APPEND, do_snaps);
  snapshot(alad_snapshot, polyNumericVector(gb_forces[2]), "ala_neck2_gb_frc",
           NumberFormat::STANDARD_REAL, "Forces due to Mongan \"Neck\" (model II) Generalized "
           "Born interactions in the alanine dipeptide system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::APPEND, do_snaps);

  // Attempt some finite difference force calculations to check the self-consistency of the
  // Generalized Born energies and forces.
  const std::vector<ImplicitSolventModel> conditions = { ImplicitSolventModel::HCT_GB,
                                                         ImplicitSolventModel::OBC_GB,
                                                         ImplicitSolventModel::OBC_GB_II,
                                                         ImplicitSolventModel::NECK_GB,
                                                         ImplicitSolventModel::NECK_GB_II };
  const std::vector<AtomicRadiusSet> shapes = { AtomicRadiusSet::MBONDI,
                                                AtomicRadiusSet::MBONDI2,
                                                AtomicRadiusSet::AMBER6,
                                                AtomicRadiusSet::BONDI,
                                                AtomicRadiusSet::MBONDI3 };
  for (size_t i = 0; i < shapes.size(); i++) {
    switch(conditions[i]) {
    case ImplicitSolventModel::HCT_GB:
      section(2);
      break;
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
      section(3);
      break;
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      section(4);
      break;
    case ImplicitSolventModel::NONE:
      break;
    }
    
    // Set conditions and re-evaluate Generalized Born forces in each system.  Skip DHFR as it
    // is expensive, less stable for comparing finite difference results, and not telling us much
    // at this stage.
    for (size_t j = 0; j < system_count; j++) {
      if (all_topologies[j]->getFileName() == dhfr_top_name) {
        continue;
      }
      all_topologies[j]->setImplicitSolventModel(conditions[i], 80.0, 0.0, shapes[i]);
      all_coords[j]->initializeForces();
      timer.assignTime(0);
      gb_energy[j] = evaluateGeneralizedBornEnergy(*(all_topologies[j]), all_semasks[j], ngb_tab,
                                                   all_coords[j], &all_systems_sc, evfrc, j);
      timer.assignTime(gb_timings);
      gb_forces[j] = all_coords[j]->getInterlacedCoordinates(tkind);
    }

    // Prepare samples of the forces in each system
    const int trpi_natom = trpi_ag.getAtomCount();
    const int dhfr_natom = dhfr_ag.getAtomCount();
    const int alad_natom = alad_ag.getAtomCount();
    std::vector<double> trpi_fsample;
    for (int j = 0; j < trpi_natom; j += 60) {
      trpi_fsample.push_back(gb_forces[0][(3 * j)    ]);
      trpi_fsample.push_back(gb_forces[0][(3 * j) + 1]);
      trpi_fsample.push_back(gb_forces[0][(3 * j) + 2]);
    }
    std::vector<double> alad_fsample;
    for (int j = 0; j < alad_natom; j += 1) {
      alad_fsample.push_back(gb_forces[2][(3 * j)    ]);
      alad_fsample.push_back(gb_forces[2][(3 * j) + 1]);
      alad_fsample.push_back(gb_forces[2][(3 * j) + 2]);
    }
    const std::vector<double> trpi_fd = forceByFiniteDifference(trpi_ag, trpi_se, &trpi_ps,
                                                                ngb_tab, 60);
    const std::vector<double> alad_fd = forceByFiniteDifference(alad_ag, alad_se, &alad_ps,
                                                                ngb_tab, 1);
    check(trpi_fsample, RelationalOperator::EQUAL, Approx(trpi_fd).margin(1.0e-4),
          "Forces computed in the Trp-cage system with " +
          getEnumerationName(conditions[i]) + " and " +
          getEnumerationName(shapes[i]) + " do not agree with their finite-difference "
          "approximations.", do_tests);
    check(alad_fsample, RelationalOperator::EQUAL, Approx(alad_fd).margin(1.0e-4),
          "Forces computed in the alanine dipeptide system with " +
          getEnumerationName(conditions[i]) + " and " +
          getEnumerationName(shapes[i]) + " do not agree with their finite-difference "
          "approximations.", do_tests);
  }

  // Read a Trp-cage trajectory and evaluate the results in double, single, and fixed-precision
  // representations.
  const std::string trpcage_traj = base_crd_name + osc + "trpcage.crd";
  const bool traj_exists = (getDrivePathType(trpcage_traj) == DrivePathType::FILE &&
                            systems_exist);
  if (getDrivePathType(trpcage_traj) != DrivePathType::FILE) {
    rtWarn("A trajectory of Trp-cage conformations (isolated boundary conditions) was not found "
           "in " + trpcage_traj + ".  Check the ${STORMM_SOURCE} environment variable for "
           "validity.  Subsequent tests will be skipped.", "test_generalized_born");
  }
  const TestPriority do_traj_tests = (traj_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  const CoordinateSeries<double> trpcage_csd = traj_exists ?
                                               CoordinateSeries<double>(trpcage_traj,
                                                                        trpi_ag.getAtomCount()) :
                                               CoordinateSeries<double>();
  const CoordinateSeries<float> trpcage_csf(trpcage_csd);
  const CoordinateSeries<llint> trpcage_csi(trpcage_csd, 26);
  timer.assignTime(0);
  const int nframe = trpcage_csd.getFrameCount();
  std::vector<double> gbe_dd(nframe), gbe_fd(nframe), gbe_id(nframe);
  std::vector<double> gbe_df(nframe), gbe_ff(nframe), gbe_if(nframe);
  ScoreCard traj_sc(20, 1, 32);
  const NonbondedKit<double> trpi_nbk_d = trpi_ag.getDoublePrecisionNonbondedKit();
  const ImplicitSolventKit<double> trpi_isk_d = trpi_ag.getDoublePrecisionImplicitSolventKit();
  const NeckGeneralizedBornKit<double> ngb_kit_d = ngb_tab.dpData();
  for (int i = 0; i < nframe; i++) {
    gbe_dd[i] = evaluateGeneralizedBornEnergy<double, double>(trpi_nbk_d, trpi_se.data(),
                                                              trpi_isk_d, ngb_kit_d,
                                                              trpcage_csd.data(), &traj_sc, i);
    gbe_fd[i] = evaluateGeneralizedBornEnergy<float, double>(trpi_nbk_d, trpi_se.data(),
                                                             trpi_isk_d, ngb_kit_d,
                                                             trpcage_csf.data(), &traj_sc, i);
    gbe_id[i] = evaluateGeneralizedBornEnergy<llint, double>(trpi_nbk_d, trpi_se.data(),
                                                             trpi_isk_d, ngb_kit_d,
                                                             trpcage_csi.data(), &traj_sc, i, 28);
  }
  const NonbondedKit<float> trpi_nbk_f = trpi_ag.getSinglePrecisionNonbondedKit();
  const ImplicitSolventKit<float> trpi_isk_f = trpi_ag.getSinglePrecisionImplicitSolventKit();
  const NeckGeneralizedBornKit<float> ngb_kit_f = ngb_tab.spData();
  for (int i = 0; i < nframe; i++) {
    gbe_df[i] = evaluateGeneralizedBornEnergy<double, float>(trpi_nbk_f, trpi_se.data(),
                                                             trpi_isk_f, ngb_kit_f,
                                                             trpcage_csd.data(), &traj_sc, i);
    gbe_ff[i] = evaluateGeneralizedBornEnergy<float, float>(trpi_nbk_f, trpi_se.data(),
                                                            trpi_isk_f, ngb_kit_f,
                                                            trpcage_csf.data(), &traj_sc, i);
    gbe_if[i] = evaluateGeneralizedBornEnergy<llint, float>(trpi_nbk_f, trpi_se.data(),
                                                            trpi_isk_f, ngb_kit_f,
                                                            trpcage_csi.data(), &traj_sc, i, 22);
  }
  section(5);
  check(meanUnsignedError(gbe_dd, gbe_fd), RelationalOperator::LESS_THAN, 8.0e-6, "Representing "
        "coordinates in single-precision (calculations in double-precision) exacts a greater than "
        "expected toll on the accuracy of (Mongan Neck II) Generalized Born calculations in the "
        "Trp-cage system.", do_traj_tests);
  check(meanUnsignedError(gbe_dd, gbe_id), RelationalOperator::LESS_THAN, 8.0e-7, "Representing "
        "coordinates in fixed-precision (calculations in double-precision) exacts a greater than "
        "expected toll on the accuracy of (Mongan Neck II) Generalized Born calculations in the "
        "Trp-cage system.", do_traj_tests);
  check(meanUnsignedError(gbe_dd, gbe_df), RelationalOperator::LESS_THAN, 6.5e-5, "Performing "
        "calculations in single-precision exacts a greater than expected toll on the accuracy of "
        "(Mongan Neck II) Generalized Born calculations in the Trp-cage system.", do_traj_tests);
  check(meanUnsignedError(gbe_dd, gbe_ff), RelationalOperator::LESS_THAN, 7.5e-5, "Performing "
        "calculations in single-precision (with coordinates also represented in single-precision) "
        "exacts a greater than expected toll on the accuracy of (Mongan Neck II) Generalized Born "
        "calculations in the Trp-cage system.", do_traj_tests);
  check(meanUnsignedError(gbe_dd, gbe_ff), RelationalOperator::LESS_THAN, 8.5e-5, "Performing "
        "calculations in single-precision (with coordinates represented in fixed-precision) "
        "exacts a greater than expected toll on the accuracy of (Mongan Neck II) Generalized Born "
        "calculations in the Trp-cage system.", do_traj_tests);
  trpi_ps.initializeForces();
  ScoreCard trpi_sc(1, 1, 32);
  evaluateGeneralizedBornEnergy(trpi_ag, trpi_se, ngb_tab, &trpi_ps, &trpi_sc, EvaluateForce::YES);
  const std::vector<double> trpi_ref_gb_forces = trpi_ps.getInterlacedCoordinates(tkind);
  testNBPrecisionModel<double, float, double>(trpi_nbk_d, trpi_isk_d, ngb_kit_d, trpi_se, trpi_ps,
                                              trpi_ref_gb_forces, 7.0e-6, do_tests);
  testNBPrecisionModel<double, float, float>(trpi_nbk_f, trpi_isk_f, ngb_kit_f, trpi_se, trpi_ps,
                                             trpi_ref_gb_forces, 3.0e-5, do_tests);
  testNBPrecisionModel<double, llint, double>(trpi_nbk_d, trpi_isk_d, ngb_kit_d, trpi_se, trpi_ps,
                                              trpi_ref_gb_forces, 1.0e-5, do_tests);
  testNBPrecisionModel<double, llint, float>(trpi_nbk_f, trpi_isk_f, ngb_kit_f, trpi_se, trpi_ps,
                                             trpi_ref_gb_forces, 2.3e-5, do_tests);
  testNBPrecisionModel<float, double, double>(trpi_nbk_d, trpi_isk_d, ngb_kit_d, trpi_se, trpi_ps,
                                              trpi_ref_gb_forces, 2.0e-5, do_tests);
  testNBPrecisionModel<float, double, float>(trpi_nbk_f, trpi_isk_f, ngb_kit_f, trpi_se, trpi_ps,
                                             trpi_ref_gb_forces, 2.4e-5, do_tests);
  testNBPrecisionModel<float, float, double>(trpi_nbk_d, trpi_isk_d, ngb_kit_d, trpi_se, trpi_ps,
                                             trpi_ref_gb_forces, 1.7e-5, do_tests);
  testNBPrecisionModel<float, float, float>(trpi_nbk_f, trpi_isk_f, ngb_kit_f, trpi_se, trpi_ps,
                                            trpi_ref_gb_forces, 3.0e-5, do_tests);
  testNBPrecisionModel<float, llint, double>(trpi_nbk_d, trpi_isk_d, ngb_kit_d, trpi_se, trpi_ps,
                                             trpi_ref_gb_forces, 2.0e-5, do_tests);
  testNBPrecisionModel<float, llint, float>(trpi_nbk_f, trpi_isk_f, ngb_kit_f, trpi_se, trpi_ps,
                                            trpi_ref_gb_forces, 2.2e-5, do_tests);
  testNBPrecisionModel<llint, double, double>(trpi_nbk_d, trpi_isk_d, ngb_kit_d, trpi_se, trpi_ps,
                                              trpi_ref_gb_forces, 6.0e-7, do_tests);
  testNBPrecisionModel<llint, double, float>(trpi_nbk_f, trpi_isk_f, ngb_kit_f, trpi_se, trpi_ps,
                                             trpi_ref_gb_forces, 1.8e-5, do_tests);
  testNBPrecisionModel<llint, float, double>(trpi_nbk_d, trpi_isk_d, ngb_kit_d, trpi_se, trpi_ps,
                                             trpi_ref_gb_forces, 1.0e-5, do_tests);
  testNBPrecisionModel<llint, float, float>(trpi_nbk_f, trpi_isk_f, ngb_kit_f, trpi_se, trpi_ps,
                                            trpi_ref_gb_forces, 3.2e-5, do_tests);
  testNBPrecisionModel<llint, llint, double>(trpi_nbk_d, trpi_isk_d, ngb_kit_d, trpi_se, trpi_ps,
                                             trpi_ref_gb_forces, 1.0e-5, do_tests);
  testNBPrecisionModel<llint, llint, float>(trpi_nbk_f, trpi_isk_f, ngb_kit_f, trpi_se, trpi_ps,
                                            trpi_ref_gb_forces, 1.8e-5, do_tests);
  
  // Print results
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
