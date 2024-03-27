#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/coordinate_series.h"
#include "../../src/UnitTesting/unit_test.h"

#ifndef STORMM_USE_HPC
using stormm::double2;
#endif
using stormm::llint;
using stormm::constants::ExceptionResponse;
using stormm::data_types::getStormmScalarTypeName;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::errors::rtWarn;
using stormm::parse::NumberFormat;
using stormm::parse::polyNumericVector;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using stormm::topology::AtomGraph;
using stormm::trajectory::CoordinateFileKind;
using stormm::trajectory::CoordinateSeries;
using stormm::trajectory::PhaseSpace;
using stormm::trajectory::TrajectoryKind;
using namespace stormm::energy;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Evaluate forces for a particular combination of coordinate, force, and calculation precision.
//
// Arguments:
//   vk:                 Valence parameters for the system
//   csr:                Original coordinates of the system
//   force_accumulator:  Set of arrays in which to accumulate forces (this is a repurposing of the
//                       CoordinateSeries object)
//   sys_title:          Title of the system (for error reporting purposes)
//   bond_ref:           Reference forces due to harmonic bond terms
//   angl_ref:           Reference forces due to harmonic angles
//   dihe_ref:           Reference forces due to cosine-based dihedrals
//   ubrd_ref:           Reference forces due to Urey-Bradley harmonic angles
//   cimp_ref:           Reference forces due to CHARMM improper dihedrals
//   cmap_ref:           Reference CMAP forces
//   tol:                Tolerance for judging success.  The amount will be scaled depending on the
//                       particular test.
//   do_tests:           Indicator that tests are even possible, passed down from a check in the
//                       main program
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
void evalAlternateForces(const ValenceKit<Tcalc> vk, const CoordinateSeriesReader<Tcoord> csr,
                         CoordinateSeries<Tforce> *force_accumulator, const std::string &sys_title,
                         const std::vector<double> &bond_ref, const std::vector<double> &angl_ref,
                         const std::vector<double> &dihe_ref, const std::vector<double> &ubrd_ref,
                         const std::vector<double> &cimp_ref, const std::vector<double> &cmap_ref, 
                         const double tol, const TestPriority do_tests) {
  CoordinateSeriesWriter<Tforce> faw = force_accumulator->data();
  const Tforce zero = 0.0;
  for (int i = 0; i < faw.natom; i++) {
    faw.xcrd[i] = zero;
    faw.ycrd[i] = zero;
    faw.zcrd[i] = zero;
  }
  ScoreCard alt_sc(1, 16, 32);
  evaluateBondTerms<Tcoord, Tforce, Tcalc>(vk, csr.xcrd, csr.ycrd, csr.zcrd, csr.umat, csr.invu,
                                           csr.unit_cell, faw.xcrd, faw.ycrd, faw.zcrd, &alt_sc,
                                           EvaluateForce::YES, 0, csr.inv_gpos_scale,
                                           faw.gpos_scale);
  const CoordinateFrame bond_frame = force_accumulator->exportFrame(0);
  const std::vector<double> bond_frc = bond_frame.getInterlacedCoordinates();
  const RelationalOperator rel_eq = RelationalOperator::EQUAL;
  check(bond_frc, rel_eq, Approx(bond_ref).margin(tol), "Forces due to harmonic bond "
        "stretching interactions in the " + sys_title + " system do not meet expectations when "
        "recomputed in " + getStormmScalarTypeName<Tcalc>() + ", with accumulation in " +
        getStormmScalarTypeName<Tforce>() + ", based on " + getStormmScalarTypeName<Tcoord>() +
        " coordinates.", do_tests);
  for (int i = 0; i < faw.natom; i++) {
    faw.xcrd[i] = zero;
    faw.ycrd[i] = zero;
    faw.zcrd[i] = zero;
  }
  evaluateAngleTerms<Tcoord, Tforce, Tcalc>(vk, csr.xcrd, csr.ycrd, csr.zcrd, csr.umat, csr.invu,
                                            csr.unit_cell, faw.xcrd, faw.ycrd, faw.zcrd, &alt_sc,
                                            EvaluateForce::YES, 0, csr.inv_gpos_scale,
                                            faw.gpos_scale);
  const CoordinateFrame angl_frame = force_accumulator->exportFrame(0);
  const std::vector<double> angl_frc = angl_frame.getInterlacedCoordinates();
  check(angl_frc, rel_eq, Approx(angl_ref).margin(0.5 * tol), "Forces due to harmonic angle "
        "bending interactions in the " + sys_title + " system do not meet expectations when "
        "recomputed in " + getStormmScalarTypeName<Tcalc>() + ", with accumulation in " +
        getStormmScalarTypeName<Tforce>() + ", based on " + getStormmScalarTypeName<Tcoord>() +
        " coordinates.", do_tests);
  for (int i = 0; i < faw.natom; i++) {
    faw.xcrd[i] = zero;
    faw.ycrd[i] = zero;
    faw.zcrd[i] = zero;
  }
  evaluateDihedralTerms<Tcoord, Tforce, Tcalc>(vk, csr.xcrd, csr.ycrd, csr.zcrd, csr.umat,
                                               csr.invu, csr.unit_cell, faw.xcrd, faw.ycrd,
                                               faw.zcrd, &alt_sc, EvaluateForce::YES, 0,
                                               csr.inv_gpos_scale, faw.gpos_scale);
  const CoordinateFrame dihe_frame = force_accumulator->exportFrame(0);
  const std::vector<double> dihe_frc = dihe_frame.getInterlacedCoordinates();
  check(dihe_frc, rel_eq, Approx(dihe_ref).margin(tol), "Forces due to cosine-based "
        "dihedral interactions in the " + sys_title + " system do not meet expectations when "
        "recomputed in " + getStormmScalarTypeName<Tcalc>() + ", with accumulation in " +
        getStormmScalarTypeName<Tforce>() + ", based on " + getStormmScalarTypeName<Tcoord>() +
        " coordinates.", do_tests);
  if (vk.nubrd > 0) {
    for (int i = 0; i < faw.natom; i++) {
      faw.xcrd[i] = zero;
      faw.ycrd[i] = zero;
      faw.zcrd[i] = zero;
    }
    evaluateUreyBradleyTerms<Tcoord, Tforce, Tcalc>(vk, csr.xcrd, csr.ycrd, csr.zcrd, csr.umat,
                                                    csr.invu, csr.unit_cell, faw.xcrd, faw.ycrd,
                                                    faw.zcrd, &alt_sc, EvaluateForce::YES, 0,
                                                    csr.inv_gpos_scale, faw.gpos_scale);
    const CoordinateFrame ubrd_frame = force_accumulator->exportFrame(0);
    const std::vector<double> ubrd_frc = ubrd_frame.getInterlacedCoordinates();
    check(ubrd_frc, rel_eq, Approx(ubrd_ref).margin(0.5 * tol), "Forces due to Urey-Bradley angle "
          "stretching interactions in the " + sys_title + " system do not meet expectations when "
          "recomputed in " + getStormmScalarTypeName<Tcalc>() + ", with accumulation in " +
          getStormmScalarTypeName<Tforce>() + ", based on " + getStormmScalarTypeName<Tcoord>() +
          " coordinates.", do_tests);
  }
  if (vk.ncimp > 0) {
    for (int i = 0; i < faw.natom; i++) {
      faw.xcrd[i] = zero;
      faw.ycrd[i] = zero;
      faw.zcrd[i] = zero;
    }
    evaluateCharmmImproperTerms<Tcoord, Tforce, Tcalc>(vk, csr.xcrd, csr.ycrd, csr.zcrd, csr.umat,
                                                       csr.invu, csr.unit_cell, faw.xcrd, faw.ycrd,
                                                       faw.zcrd, &alt_sc, EvaluateForce::YES, 0,
                                                       csr.inv_gpos_scale, faw.gpos_scale);
    const CoordinateFrame cimp_frame = force_accumulator->exportFrame(0);
    const std::vector<double> cimp_frc = cimp_frame.getInterlacedCoordinates();
    check(cimp_frc, rel_eq, Approx(cimp_ref).margin(tol), "Forces due to CHARMM harmonic improper "
          "dihedral interactions in the " + sys_title + " system do not meet expectations when "
          "recomputed in " + getStormmScalarTypeName<Tcalc>() + ", with accumulation in " +
          getStormmScalarTypeName<Tforce>() + ", based on " + getStormmScalarTypeName<Tcoord>() +
          " coordinates.", do_tests);
  }
  if (vk.ncmap > 0) {
    for (int i = 0; i < faw.natom; i++) {
      faw.xcrd[i] = zero;
      faw.ycrd[i] = zero;
      faw.zcrd[i] = zero;
    }
    evaluateCmapTerms<Tcoord, Tforce, Tcalc>(vk, csr.xcrd, csr.ycrd, csr.zcrd, csr.umat, csr.invu,
                                             csr.unit_cell, faw.xcrd, faw.ycrd, faw.zcrd, &alt_sc,
                                             EvaluateForce::YES, 0, csr.inv_gpos_scale,
                                             faw.gpos_scale);
    const CoordinateFrame cmap_frame = force_accumulator->exportFrame(0);
    const std::vector<double> cmap_frc = cmap_frame.getInterlacedCoordinates();
    check(cmap_frc, rel_eq, Approx(cmap_ref).margin(tol), "Forces due to CMAP 2D spline surface "
          "interactions in the " + sys_title + " system do not meet expectations when recomputed "
          "in " + getStormmScalarTypeName<Tcalc>() + ", with accumulation in " +
          getStormmScalarTypeName<Tforce>() + ", based on " + getStormmScalarTypeName<Tcoord>() +
          " coordinates.", do_tests);
  }
}

//-------------------------------------------------------------------------------------------------
// Evaluate the various valence energies and forces for a system.
//
// Arguments:
//   ag:              System topology
//   cf:              Original coordinates of the system in the highest-precision model
//   sys_title:       System title (for error reporting purposes)
//   valence_energy:  Vector of target valence energies computed with the highest-precision model
//   bond_frc:        Reference forces due to harmonic bond terms
//   angl_frc:        Reference forces due to harmonic angles
//   dihe_frc:        Reference forces due to cosine-based dihedrals
//   ubrd_frc:        Reference forces due to Urey-Bradley harmonic angles
//   cimp_frc:        Reference forces due to CHARMM improper dihedrals
//   cmap_frc:        Reference CMAP forces
//   ef_tol:          Tolerance for energy calculations when single-precision floating point
//                    operations are expected to be the primary source of error
//   ei_tol:          Tolerance for energy calculations when fixed-precision coordinates are
//                    expected to be the primary source of error
//   frcf_tol:        Tolerance for force calculations when single-precision floating point
//                    operations are expected to be the primary source of error
//   frci_tol:        Tolerance for force calculations when fixed-precision coordinates or force
//                    accumulators are expected to be the primary source of error
//-------------------------------------------------------------------------------------------------
void evalAlternateCoords(const AtomGraph &ag, const CoordinateFrame &cf,
                         const std::string &sys_title, const std::vector<double> &valence_energy,
                         const std::vector<double> &bond_frc, const std::vector<double> &angl_frc,
                         const std::vector<double> &dihe_frc, const std::vector<double> &ubrd_frc,
                         const std::vector<double> &cimp_frc, const std::vector<double> &cmap_frc,
                         const double ef_tol, const double ei_tol, const double frcf_tol,
                         const double frci_tol, const TestPriority do_tests)
{
  const CoordinateSeries<double> csd(cf, 1);
  const CoordinateSeries<float> csf(cf, 1);
  const CoordinateSeries<llint> csi(cf, 1, 26);
  const ValenceKit<double> vkd = ag.getDoublePrecisionValenceKit();  
  const ValenceKit<float> vkf = ag.getSinglePrecisionValenceKit();  
  ScoreCard alt_sc(1, 16, 32);
  const double2 dihe_energy_fd =
    evaluateDihedralTerms<float, double>(vkd, csf.data(), &alt_sc, 0);
  const std::vector<double> valence_energy_fd = {
    evaluateBondTerms<float, double>(vkd, csf.data(), &alt_sc, 0),
    evaluateAngleTerms<float, double>(vkd, csf.data(), &alt_sc, 0),
    dihe_energy_fd.x, dihe_energy_fd.y,
    evaluateUreyBradleyTerms<float, double>(vkd, csf.data(), &alt_sc, 0),
    evaluateCharmmImproperTerms<float, double>(vkd, csf.data(), &alt_sc, 0),
    evaluateCmapTerms<float, double>(vkd, csf.data(), &alt_sc, 0)
  };
  const double2 dihe_energy_id =
    evaluateDihedralTerms<llint, double>(vkd, csi.data(), &alt_sc, 0);
  const std::vector<double> valence_energy_id = {
    evaluateBondTerms<llint, double>(vkd, csi.data(), &alt_sc, 0),
    evaluateAngleTerms<llint, double>(vkd, csi.data(), &alt_sc, 0),
    dihe_energy_id.x, dihe_energy_id.y,
    evaluateUreyBradleyTerms<llint, double>(vkd, csi.data(), &alt_sc, 0),
    evaluateCharmmImproperTerms<llint, double>(vkd, csi.data(), &alt_sc, 0),
    evaluateCmapTerms<llint, double>(vkd, csi.data(), &alt_sc, 0)
  };
  section(3);
  const RelationalOperator rel_eq = RelationalOperator::EQUAL;
  check(valence_energy_fd, rel_eq, Approx(valence_energy).margin(ef_tol), sys_title +
        " valence energies do not meet expectations when calculated in double-precision using a "
        "single-precision coordinate representation.", do_tests);
  check(valence_energy_id, rel_eq, Approx(valence_energy).margin(ei_tol), sys_title +
        " valence energies do not meet expectations when calculated in double-precision using a "
        "fixed-precision coordinate representation.", do_tests);
  CoordinateSeries<double> fv_d(csf);
  CoordinateSeries<float> fv_f(csf);
  CoordinateSeries<llint> fv_i(csf, 23);
  const CoordinateSeriesReader<double> csdr = csd.data();
  const CoordinateSeriesReader<float> csfr  = csf.data();
  const CoordinateSeriesReader<llint> csir  = csi.data();
  CoordinateSeriesWriter<double> fv_dw = fv_d.data();
  CoordinateSeriesWriter<float> fv_fw  = fv_f.data();
  CoordinateSeriesWriter<llint> fv_iw  = fv_i.data();
  section(4);
  evalAlternateForces(vkd, csdr, &fv_d, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc, stormm::constants::tiny, do_tests);
  evalAlternateForces(vkd, csfr, &fv_d, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc, frcf_tol, do_tests);
  evalAlternateForces(vkd, csir, &fv_d, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc, frci_tol, do_tests);
  evalAlternateForces(vkd, csdr, &fv_f, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc, frcf_tol, do_tests);
  evalAlternateForces(vkd, csfr, &fv_f, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc, frcf_tol, do_tests);
  evalAlternateForces(vkd, csir, &fv_f, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc, frcf_tol, do_tests);
  evalAlternateForces(vkd, csdr, &fv_i, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc, frcf_tol, do_tests);
  evalAlternateForces(vkd, csfr, &fv_i, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc, frcf_tol, do_tests);
  evalAlternateForces(vkd, csir, &fv_i, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc, frci_tol, do_tests);
  evalAlternateForces(vkf, csdr, &fv_d, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc,       frcf_tol, do_tests);
  evalAlternateForces(vkf, csfr, &fv_d, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc, 2.0 * frcf_tol, do_tests);
  evalAlternateForces(vkf, csir, &fv_d, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc,       frcf_tol, do_tests);
  evalAlternateForces(vkf, csdr, &fv_f, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc, 2.0 * frcf_tol, do_tests);
  evalAlternateForces(vkf, csfr, &fv_f, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc, 3.0 * frcf_tol, do_tests);
  evalAlternateForces(vkf, csir, &fv_f, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc, 2.0 * frcf_tol, do_tests);
  evalAlternateForces(vkf, csdr, &fv_i, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc,       frcf_tol, do_tests);
  evalAlternateForces(vkf, csfr, &fv_i, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc, 2.0 * frcf_tol, do_tests);
  evalAlternateForces(vkf, csir, &fv_i, sys_title, bond_frc, angl_frc, dihe_frc, ubrd_frc,
                      cimp_frc, cmap_frc,       frcf_tol, do_tests);
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
  section("Valence energy evaluation");

  // Section 2
  section("Valence force evaluation");

  // Section 3
  section("Energy with alternate representations");

  // Section 4
  section("Forces with alternate representations");

  // Locate topologies and coordinate files
  const char osc = osSeparator();
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string base_ptl_name = oe.getStormmSourcePath() + osc + "test" + osc + "Potential";
  const std::string trpcage_top_name = base_top_name + osc + "trpcage.top";
  const std::string trpcage_crd_name = base_crd_name + osc + "trpcage.inpcrd";
  const std::string dhfr_top_name = base_top_name + osc + "dhfr_cmap.top";
  const std::string dhfr_crd_name = base_crd_name + osc + "dhfr_cmap.inpcrd";
  const std::string alad_top_name = base_top_name + osc + "ala_dipeptide.top";
  const std::string alad_crd_name = base_crd_name + osc + "ala_dipeptide.inpcrd";
  const bool systems_exist = (getDrivePathType(trpcage_top_name) == DrivePathType::FILE &&
                              getDrivePathType(trpcage_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(dhfr_top_name) == DrivePathType::FILE &&
                              getDrivePathType(dhfr_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(alad_top_name) == DrivePathType::FILE &&
                              getDrivePathType(alad_crd_name) == DrivePathType::FILE);
  const TestPriority do_tests = (systems_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (systems_exist == false) {
    rtWarn("Files for the Trp-cage miniprotein, the DHFR globular protein (with CHARMM potential "
           "details), and alanine dipeptide (with ff19SB) were not found.  These files should be "
           "found in the ${STORMM_SOURCE}/test/Topology and ${STORMM_SOURCE}/test/Trajectory "
           "directories.  Check the $STORMM_SOURCE environment variable.  A number of tests will "
           "be skipped.", "test_valence_evaluation");
  }
  AtomGraph trpcage_ag, dhfr_ag, alad_ag;
  PhaseSpace trpcage_ps, dhfr_ps, alad_ps;
  if (systems_exist) {
    trpcage_ag.buildFromPrmtop(trpcage_top_name, ExceptionResponse::SILENT);
    trpcage_ps.buildFromFile(trpcage_crd_name, CoordinateFileKind::AMBER_INPCRD);
    dhfr_ag.buildFromPrmtop(dhfr_top_name, ExceptionResponse::SILENT);
    dhfr_ps.buildFromFile(dhfr_crd_name, CoordinateFileKind::AMBER_INPCRD);
    alad_ag.buildFromPrmtop(alad_top_name, ExceptionResponse::SILENT);
    alad_ps.buildFromFile(alad_crd_name, CoordinateFileKind::AMBER_INPCRD);
  }
  ScoreCard all_systems_sc(3);
  const int trpcage_idx = 0;
  const int dhfr_idx = 1;
  const int alad_idx = 2;
  const TrajectoryKind tkind = TrajectoryKind::FORCES;
  trpcage_ps.initializeForces();
  const double trpcage_bond_e = evaluateBondTerms(trpcage_ag, &trpcage_ps, &all_systems_sc,
                                                  EvaluateForce::YES, trpcage_idx);
  const std::vector<double> trpcage_bond_frc = trpcage_ps.getInterlacedCoordinates(tkind);
  trpcage_ps.initializeForces();
  const double trpcage_angl_e = evaluateAngleTerms(trpcage_ag, &trpcage_ps, &all_systems_sc,
                                                   EvaluateForce::YES, trpcage_idx);
  const std::vector<double> trpcage_angl_frc = trpcage_ps.getInterlacedCoordinates(tkind);
  trpcage_ps.initializeForces();
  const double2 trpcage_dihe_e = evaluateDihedralTerms(trpcage_ag, &trpcage_ps, &all_systems_sc,
                                                       EvaluateForce::YES, trpcage_idx);
  const std::vector<double> trpcage_dihe_frc = trpcage_ps.getInterlacedCoordinates(tkind);
  trpcage_ps.initializeForces();
  const double trpcage_ubrd_e = evaluateUreyBradleyTerms(trpcage_ag, &trpcage_ps, &all_systems_sc,
                                                         EvaluateForce::YES, trpcage_idx);
  const double trpcage_cimp_e = evaluateCharmmImproperTerms(trpcage_ag, &trpcage_ps,
                                                            &all_systems_sc, EvaluateForce::YES,
                                                            trpcage_idx);
  const double trpcage_cmap_e = evaluateCmapTerms(trpcage_ag, &trpcage_ps, &all_systems_sc,
                                                  EvaluateForce::YES, trpcage_idx);
  dhfr_ps.initializeForces();
  const double dhfr_bond_e = evaluateBondTerms(dhfr_ag, &dhfr_ps, &all_systems_sc,
                                               EvaluateForce::YES, dhfr_idx);
  const std::vector<double> dhfr_bond_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  dhfr_ps.initializeForces();
  const double dhfr_angl_e = evaluateAngleTerms(dhfr_ag, &dhfr_ps, &all_systems_sc,
                                                EvaluateForce::YES, dhfr_idx);
  const std::vector<double> dhfr_angl_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  dhfr_ps.initializeForces();
  const double2 dhfr_dihe_e = evaluateDihedralTerms(dhfr_ag, &dhfr_ps, &all_systems_sc,
                                                    EvaluateForce::YES, dhfr_idx);
  const std::vector<double> dhfr_dihe_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  dhfr_ps.initializeForces();
  const double dhfr_ubrd_e = evaluateUreyBradleyTerms(dhfr_ag, &dhfr_ps, &all_systems_sc,
                                                      EvaluateForce::YES, dhfr_idx);
  const std::vector<double> dhfr_ubrd_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  dhfr_ps.initializeForces();
  const double dhfr_cimp_e = evaluateCharmmImproperTerms(dhfr_ag, &dhfr_ps, &all_systems_sc,
                                                         EvaluateForce::YES, dhfr_idx);
  const std::vector<double> dhfr_cimp_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  dhfr_ps.initializeForces();
  const double dhfr_cmap_e = evaluateCmapTerms(dhfr_ag, &dhfr_ps, &all_systems_sc,
                                               EvaluateForce::YES, dhfr_idx);
  const std::vector<double> dhfr_cmap_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  alad_ps.initializeForces();
  const double alad_bond_e = evaluateBondTerms(alad_ag, &alad_ps, &all_systems_sc,
                                               EvaluateForce::YES, alad_idx);
  const std::vector<double> alad_bond_frc = alad_ps.getInterlacedCoordinates(tkind);
  alad_ps.initializeForces();
  const double alad_angl_e = evaluateAngleTerms(alad_ag, &alad_ps, &all_systems_sc,
                                                EvaluateForce::YES, alad_idx);
  const std::vector<double> alad_angl_frc = alad_ps.getInterlacedCoordinates(tkind);
  alad_ps.initializeForces();
  const double2 alad_dihe_e = evaluateDihedralTerms(alad_ag, &alad_ps, &all_systems_sc,
                                                    EvaluateForce::YES, alad_idx);
  const std::vector<double> alad_dihe_frc = alad_ps.getInterlacedCoordinates(tkind);
  alad_ps.initializeForces();
  const double alad_ubrd_e = evaluateUreyBradleyTerms(alad_ag, &alad_ps, &all_systems_sc,
                                                      EvaluateForce::YES, alad_idx);
  const double alad_cimp_e = evaluateCharmmImproperTerms(alad_ag, &alad_ps, &all_systems_sc,
                                                         EvaluateForce::YES, alad_idx);
  const double alad_cmap_e = evaluateCmapTerms(alad_ag, &alad_ps, &all_systems_sc,
                                               EvaluateForce::YES, alad_idx);
  const std::vector<double> alad_cmap_frc = alad_ps.getInterlacedCoordinates(tkind);

  // Collect pertinent energy results
  const std::vector<double> trpcage_acc = all_systems_sc.reportInstantaneousStates(trpcage_idx);
  const std::vector<double> dhfr_acc    = all_systems_sc.reportInstantaneousStates(dhfr_idx);
  const std::vector<double> alad_acc    = all_systems_sc.reportInstantaneousStates(alad_idx);

  const std::vector<double> trpcage_valence_energy_answer =
    {   12.38440829,   69.48736586, -142.41607189,    0.77975874,    0.00000000,    0.00000000,
         0.00000000 };
  const std::vector<double> dhfr_valence_energy_answer =
    {  147.47163346,  439.91217087,  754.04779195,    0.00000000,   31.87732134,   18.58590078,
       -85.10501247 };
  const std::vector<double> alad_valence_energy_answer =
    {    0.38340254,    0.54479669,    2.42761153,    0.00708592,    0.00000000,    0.00000000,
        -0.38861139 };
  const std::vector<double> trpcage_valence_energy_result =
    { trpcage_bond_e, trpcage_angl_e, trpcage_dihe_e.x, trpcage_dihe_e.y, trpcage_ubrd_e,
      trpcage_cimp_e, trpcage_cmap_e };
  const std::vector<double> dhfr_valence_energy_result =
    { dhfr_bond_e, dhfr_angl_e, dhfr_dihe_e.x, dhfr_dihe_e.y, dhfr_ubrd_e, dhfr_cimp_e,
      dhfr_cmap_e };
  const std::vector<double> alad_valence_energy_result =
    { alad_bond_e, alad_angl_e, alad_dihe_e.x, alad_dihe_e.y, alad_ubrd_e, alad_cimp_e,
      alad_cmap_e };
  section(1);
  const RelationalOperator rel_eq = RelationalOperator::EQUAL;
  check(trpcage_valence_energy_result, rel_eq, trpcage_valence_energy_answer,
        "Trp-cage valence energies do not meet expectations.", do_tests);
  check(dhfr_valence_energy_result, rel_eq, dhfr_valence_energy_answer,
        "DHFR valence energies do not meet expectations.", do_tests);
  check(alad_valence_energy_result, rel_eq, alad_valence_energy_answer,
        "Alanine dipeptide valence energies do not meet expectations.", do_tests);

  // Re-compute valence energies with a CoordinateFrame abstract (no force computations)
  ScoreCard secondary_sc(3);
  const CoordinateFrame trpcage_cf(trpcage_ps);
  const CoordinateFrame dhfr_cf(dhfr_ps);
  const CoordinateFrame alad_cf(alad_ps);
  const double2 trpcage_dihe_energy_ii = evaluateDihedralTerms(trpcage_ag, trpcage_cf,
                                                               &secondary_sc, trpcage_idx);
  const std::vector<double> trpcage_valence_energy_ii = {
    evaluateBondTerms(trpcage_ag, trpcage_cf, &secondary_sc, trpcage_idx),
    evaluateAngleTerms(trpcage_ag, trpcage_cf, &secondary_sc, trpcage_idx),
    trpcage_dihe_energy_ii.x, trpcage_dihe_energy_ii.y,
    evaluateUreyBradleyTerms(trpcage_ag, trpcage_cf, &secondary_sc, trpcage_idx),
    evaluateCharmmImproperTerms(trpcage_ag, trpcage_cf, &secondary_sc, trpcage_idx),
    evaluateCmapTerms(trpcage_ag, trpcage_cf, &secondary_sc, trpcage_idx)
  };
  const double2 dhfr_dihe_energy_ii = evaluateDihedralTerms(dhfr_ag, dhfr_cf, &secondary_sc,
                                                            dhfr_idx);
  const std::vector<double> dhfr_valence_energy_ii = {
    evaluateBondTerms(dhfr_ag, dhfr_cf, &secondary_sc, dhfr_idx),
    evaluateAngleTerms(dhfr_ag, dhfr_cf, &secondary_sc, dhfr_idx),
    dhfr_dihe_energy_ii.x, dhfr_dihe_energy_ii.y,
    evaluateUreyBradleyTerms(dhfr_ag, dhfr_cf, &secondary_sc, dhfr_idx),
    evaluateCharmmImproperTerms(dhfr_ag, dhfr_cf, &secondary_sc, dhfr_idx),
    evaluateCmapTerms(dhfr_ag, dhfr_cf, &secondary_sc, dhfr_idx)
  };
  const double2 alad_dihe_energy_ii = evaluateDihedralTerms(alad_ag, alad_cf, &secondary_sc,
                                                            alad_idx);
  const std::vector<double> alad_valence_energy_ii = {
    evaluateBondTerms(alad_ag, alad_cf, &secondary_sc, alad_idx),
    evaluateAngleTerms(alad_ag, alad_cf, &secondary_sc, alad_idx),
    alad_dihe_energy_ii.x, alad_dihe_energy_ii.y,
    evaluateUreyBradleyTerms(alad_ag, alad_cf, &secondary_sc, alad_idx),
    evaluateCharmmImproperTerms(alad_ag, alad_cf, &secondary_sc, alad_idx),
    evaluateCmapTerms(alad_ag, alad_cf, &secondary_sc, alad_idx)
  };
  check(trpcage_valence_energy_ii, rel_eq, trpcage_valence_energy_answer,
        "Trp-cage valence energies do not meet expectations when computed with a CoordinateFrame "
        "abstract.", do_tests);
  check(dhfr_valence_energy_ii, rel_eq, dhfr_valence_energy_answer,
        "DHFR valence energies do not meet expectations when computed with a CoordinateFrame "
        "abstract.", do_tests);
  check(alad_valence_energy_ii, rel_eq, alad_valence_energy_answer,
        "Alanine dipeptide valence energies do not meet expectations when computed with a "
        "CoordinateFrame abstract.", do_tests);
  
  // Check for the existence of snapshot files
  const std::string trpcage_snapshot(base_ptl_name + osc + "trpcage_details.m");
  const std::string dhfr_snapshot(base_ptl_name + osc + "dhfr_details.m");
  const std::string alad_snapshot(base_ptl_name + osc + "ala_dipeptide_details.m");
  const bool snps_exist = (getDrivePathType(trpcage_snapshot) == DrivePathType::FILE &&
                           getDrivePathType(dhfr_snapshot) == DrivePathType::FILE &&
                           getDrivePathType(alad_snapshot) == DrivePathType::FILE);
  const TestPriority snap_check = (snps_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (snps_exist == false && oe.takeSnapshot() != SnapshotOperation::SNAPSHOT) {
    rtWarn("Snapshot files " + alad_snapshot + " were not found.  These files contain reference "
           "forces for checking the valence energy derivative calculations.  Check that the "
           "${STORMM_SOURCE} environment variable is set properly so that these snapshot files "
           "may be found in ${STORMM_SOURCE}/test/Potential/.  Subsequent tests will be skipped "
           "until these reference files are available.", "test_valence_evaluation");
  }
  section(2);
  snapshot(trpcage_snapshot, polyNumericVector(trpcage_bond_frc), "trpcage_bond",
           NumberFormat::SCIENTIFIC, "Forces due to harmonic bond stretching interactions in the "
           "Trp-cage system do not meet expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12,
           PrintSituation::OVERWRITE, snap_check);
  snapshot(trpcage_snapshot, polyNumericVector(trpcage_angl_frc), "trpcage_angl",
           NumberFormat::SCIENTIFIC, "Forces due to harmonic angle bending interactions in the "
           "Trp-cage system do not meet expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12,
           PrintSituation::APPEND, snap_check);
  snapshot(trpcage_snapshot, polyNumericVector(trpcage_dihe_frc), "trpcage_dihe",
           NumberFormat::SCIENTIFIC, "Forces due to cosine-based dihedral (proper and improper) "
           "interactions in the Trp-cage system do not meet expectations.", oe.takeSnapshot(),
           1.0e-8, 1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_bond_frc), "dhfr_bond", NumberFormat::SCIENTIFIC,
           "Forces due to harmonic bond stretching interactions in the DHFR system do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::OVERWRITE,
           snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_angl_frc), "dhfr_angl", NumberFormat::SCIENTIFIC,
           "Forces due to harmonic angle bending interactions in the DHFR system do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::APPEND,
           snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_dihe_frc), "dhfr_dihe", NumberFormat::SCIENTIFIC,
           "Forces due to cosine-based dihedral interactions in the DHFR system do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::APPEND,
           snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_ubrd_frc), "dhfr_ubrd", NumberFormat::SCIENTIFIC,
           "Forces due to Urey-Bradley interactions in the DHFR system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_cimp_frc), "dhfr_cimp", NumberFormat::SCIENTIFIC,
           "Forces due to CHARMM improper terms in the DHFR system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_cmap_frc), "dhfr_cmap", NumberFormat::SCIENTIFIC,
           "Forces due to CMAP terms in the DHFR system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(alad_snapshot, polyNumericVector(alad_bond_frc), "ala_dipeptide_bond_forces",
           NumberFormat::SCIENTIFIC, "Forces due to harmonic bond stretching interactions in the "
           "Ala dipeptide system do not meet expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12,
           PrintSituation::OVERWRITE, snap_check);
  snapshot(alad_snapshot, polyNumericVector(alad_angl_frc), "ala_dipeptide_angl_forces",
           NumberFormat::SCIENTIFIC, "Forces due to harmonic angle bending interactions in the "
           "Ala dipeptide system do not meet expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12,
           PrintSituation::APPEND, snap_check);
  snapshot(alad_snapshot, polyNumericVector(alad_dihe_frc), "ala_dipeptide_dihe_forces",
           NumberFormat::SCIENTIFIC, "Forces due to cosine-based dihedral (proper and improper) "
           "interactions in the Ala dipeptide system do not meet expectations.", oe.takeSnapshot(),
           1.0e-8, 1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(alad_snapshot, polyNumericVector(alad_cmap_frc), "ala_dipeptide_cmap_forces",
           NumberFormat::SCIENTIFIC, "Forces due to ff19SB CMAP terms in the Ala dipeptide system "
           "do not meet expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::APPEND,
           snap_check);

  // Make alternate coordinate representations and recompute the energies
  evalAlternateCoords(trpcage_ag, trpcage_cf, "Trp-cage", trpcage_valence_energy_ii,
                      trpcage_bond_frc, trpcage_angl_frc, trpcage_dihe_frc, trpcage_bond_frc,
                      trpcage_bond_frc, trpcage_bond_frc, 2.0e-4, 1.0e-5, 9.0e-4, 2.5e-5,
                      do_tests);
  evalAlternateCoords(dhfr_ag, dhfr_cf, "DHFR", dhfr_valence_energy_ii, dhfr_bond_frc,
                      dhfr_angl_frc, dhfr_dihe_frc, dhfr_ubrd_frc, dhfr_cimp_frc, dhfr_cmap_frc,
                      5.0e-4, 1.0e-5, 7.0e-3, 2.5e-5, do_tests);
  evalAlternateCoords(alad_ag, alad_cf, "Alanine dipeptide", alad_valence_energy_ii, alad_bond_frc,
                      alad_angl_frc, alad_dihe_frc, alad_bond_frc, alad_bond_frc, alad_cmap_frc,
                      1.0e-4, 1.0e-5, 7.0e-4, 1.5e-5, do_tests);

  // Print results
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
