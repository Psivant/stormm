#include "copyright.h"
#include "../../src/Chemistry/atom_equivalence.h"
#include "../../src/Chemistry/chemical_features.h"
#include "../../src/Chemistry/chemistry_enumerators.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/clash_detection.h"
#include "../../src/Structure/global_manipulation.h"
#include "../../src/Structure/isomerization.h"
#include "../../src/Structure/local_arrangement.h"
#include "../../src/Structure/rmsd.h"
#include "../../src/Structure/rmsd_plan.h"
#include "../../src/Structure/structure_enumerators.h"
#include "../../src/Synthesis/condensate.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/synthesis_enumerators.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/coordinate_copy.h"
#include "../../src/Trajectory/coordinate_series.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::constants::tiny;
using stormm::constants::warp_size_int;
using stormm::constants::warp_size_zu;
#ifndef STORMM_USE_HPC
using stormm::data_types::int4;
#endif
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getBaseName;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::energy::evaluateBondTerms;
using stormm::energy::evaluateAngleTerms;
using stormm::energy::ScoreCard;
using stormm::errors::rtWarn;
using stormm::parse::char4ToString;
using stormm::random::Xoshiro256ppGenerator;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using stormm::stmath::computeBoxTransform;
using stormm::stmath::meanUnsignedError;
using stormm::symbols::twopi;
using stormm::symbols::pi;
using namespace stormm::chemistry;
using namespace stormm::structure;
using namespace stormm::synthesis;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Check rotatable bonds throughout a structure.  If the topology supplied pertains to a system
// with many molecules, only the first will be analyzed.  Options are provided to take snapshots
// of detailed particle positions.
//
// Arguments:
//   ag:            System topology
//   ps:            Coordinates of the system
//   chemfe:        Chemical features of the system
//   oe:            System environment variables (for snapshotting behavior)
//   snp_var_name:  The variable name under which to store or access snapshot data
//   expectation:   The manner in which to approach the snapshot file, if a variable name has been
//                  supplied to direct the flow of data into it (when updating snapshot files)
//-------------------------------------------------------------------------------------------------
void checkRotationalSampling(const AtomGraph &ag, const PhaseSpace &ps,
                             const ChemicalFeatures &chemfe, const TestEnvironment &oe,
                             const TestPriority do_tests,
                             const std::string &snp_var_name = std::string(""),
                             const PrintSituation expectation = PrintSituation::APPEND) {
  ScoreCard sc(1);
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  const CoordinateFrameReader cfr(ps);
  const double orig_bond_e = evaluateBondTerms(vk, cfr, &sc, 0);
  const double orig_angl_e = evaluateAngleTerms(vk, cfr, &sc, 0);
  section(1);
  
  // Rotate about various bonds.  This will generate all sorts of clashes, but bond and angle
  // energies should be unaffected.
  const std::vector<IsomerPlan> rt_grp = chemfe.getRotatableBondGroups();
  const int nrt = rt_grp.size();
  PhaseSpace rotation_copy(ps);
  PhaseSpaceWriter psw = rotation_copy.data();
  const PhaseSpaceReader psr = ps.data();
  std::vector<double> rot_crd(3 * (cdk.mol_limits[1] - cdk.mol_limits[0]) * nrt);
  std::vector<double> repos_dev(nrt), ubond_dev(nrt), uangl_dev(nrt);
  int rcpos = 0;
  for (int i = 0; i < nrt; i++) {
    rotateAboutBond(&rotation_copy, rt_grp[i].getRootAtom(), rt_grp[i].getPivotAtom(),
                    rt_grp[i].getMovingAtoms(), 2.0 * stormm::symbols::pi / 3.0);

    // Record the positions of atoms in the first molecule
    for (int j = cdk.mol_limits[0]; j < cdk.mol_limits[1]; j++) {
      rot_crd[rcpos] = psw.xcrd[cdk.mol_contents[j]];
      rcpos++;
      rot_crd[rcpos] = psw.ycrd[cdk.mol_contents[j]];
      rcpos++;
      rot_crd[rcpos] = psw.zcrd[cdk.mol_contents[j]];
      rcpos++;
    }

    // Record the bond and angle energies
    const double new_bond_e = evaluateBondTerms(vk, CoordinateFrameReader(rotation_copy), &sc, 0);
    const double new_angl_e = evaluateAngleTerms(vk, CoordinateFrameReader(rotation_copy), &sc, 0);
    ubond_dev[i] = fabs(new_bond_e - orig_bond_e);
    uangl_dev[i] = fabs(new_angl_e - orig_angl_e);

    // Reverse the rotation
    rotateAboutBond(&rotation_copy, rt_grp[i].getRootAtom(), rt_grp[i].getPivotAtom(),
                    rt_grp[i].getMovingAtoms(), -2.0 * stormm::symbols::pi / 3.0);

    // Check that the molecule was returned to its original state
    repos_dev[i] = rmsd(ps, rotation_copy, ag, RMSDMethod::ALIGN_GEOM, cdk.mol_limits[0],
                        cdk.mol_limits[1]);
  }
  const char osc = osSeparator();
  const std::string base_iso_path = oe.getStormmSourcePath() + osc + "test" + osc + "Structure";
  const std::string rcrd_snapshot = base_iso_path + osc + "rotated_coords.m";
  const bool snps_exist = (getDrivePathType(rcrd_snapshot) == DrivePathType::FILE);
  if (snps_exist == false && oe.takeSnapshot() == SnapshotOperation::COMPARE) {
    rtWarn("The snapshot file " + rcrd_snapshot + " was not found.  Check the ${STORMM_SOURCE} "
           "environment variable, currently set to " + oe.getStormmSourcePath() +
           ", for validity.  Subsequent tests will be skipped.", "test_isomerization");
  }
  const TestPriority do_snps = (snps_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (snp_var_name.size() > 0) {
    snapshot(rcrd_snapshot, polyNumericVector(rot_crd), snp_var_name, 1.0e-6, "Coordinates "
             "obtained from rotation of selected bonds do not meet expectations.",
             oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL, expectation, do_snps);
  }
  const Approx target(std::vector<double>(nrt, 0.0), ComparisonType::ABSOLUTE, 1.0e-6); 
  check(repos_dev, RelationalOperator::EQUAL, target, "Reversing the rotational operations of "
        "selected bonds does not return the molecule to its original state.", do_tests);
  check(ubond_dev, RelationalOperator::EQUAL, target, "Bond energies are changed by rotation "
        "about a bond.", do_tests);
  check(uangl_dev, RelationalOperator::EQUAL, target, "Angle energies are changed by rotation "
        "about a bond.", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Check chiral centers throughout a structure.  If the topology supplied pertains to a system
// with many molecules, only the first will be analyzed.  Options are provided to take snapshots
// of detailed particle positions resulting from inverting the various chiral centers.
//
// Arguments:
//   ag:            System topology
//   ps:            Coordinates of the system
//   chemfe:        Chemical features of the system
//   oe:            System environment variables (for snapshotting behavior)
//   snp_var_name:  The variable name under which to store or access snapshot data
//   expectation:   The manner in which to approach the snapshot file, if a variable name has been
//                  supplied to direct the flow of data into it (when updating snapshot files)
//-------------------------------------------------------------------------------------------------
void checkChiralSampling(const AtomGraph &ag, const PhaseSpace &ps,
                         const ChemicalFeatures &chemfe, const TestEnvironment &oe,
                         const TestPriority do_tests,
                         const std::string &snp_var_name = std::string(""),
                         const PrintSituation expectation = PrintSituation::APPEND) {
  ScoreCard sc(1);
  const std::vector<IsomerPlan> inv_grp = chemfe.getChiralInversionGroups();
  const std::vector<ChiralInversionProtocol> protocols = chemfe.getChiralInversionMethods();
  const std::vector<int> centers = chemfe.getChiralCenters();
  const int nchiral = protocols.size();
  if (nchiral == 0) {
    return;
  }
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  const CoordinateFrameReader cfr(ps);
  const double orig_bond_e = evaluateBondTerms(vk, cfr, &sc, 0);
  PhaseSpace inversion_copy(ps);
  PhaseSpaceWriter psw = inversion_copy.data();
  const PhaseSpaceReader psr = ps.data();
  int invcpos = 0;
  std::vector<double> inverted_crd(3 * nchiral * (cdk.mol_limits[1] - cdk.mol_limits[0]));
  std::vector<double> repos_dev(nchiral), ubond_dev(nchiral);
  const std::vector<int4> cbase = chemfe.getChiralArmBaseAtoms();
  for (int i = 0; i < nchiral; i++) {
    const int orig_chiral_orientation = getChiralOrientation(inversion_copy, centers[i],
                                                             cbase[i].x, cbase[i].y, cbase[i].z,
                                                             cbase[i].w);
    flipChiralCenter(&inversion_copy, i, centers, protocols, inv_grp);
    const int new_chiral_orientation = getChiralOrientation(inversion_copy, centers[i], cbase[i].x,
                                                            cbase[i].y, cbase[i].z, cbase[i].w);
    if (chemfe.getChiralInversionMethods(i) != ChiralInversionProtocol::DO_NOT_INVERT) {
      check(orig_chiral_orientation * new_chiral_orientation == -1, "Chiral center " +
            std::to_string(i) + ", between atoms " + char4ToString(cdk.atom_names[cbase[i].x]) +
            ", " + char4ToString(cdk.atom_names[cbase[i].y]) + ", " +
            char4ToString(cdk.atom_names[cbase[i].z]) + ", and " +
            char4ToString(cdk.atom_names[cbase[i].w]) + " surrounding atom " +
            char4ToString(cdk.atom_names[centers[i]]) + " in the topology based on " +
            getBaseName(ag.getFileName()) + ", was not inverted as it should have "
            "been.", do_tests);
    }
    
    // Record the positions of atoms in the first molecule
    for (int j = cdk.mol_limits[0]; j < cdk.mol_limits[1]; j++) {
      inverted_crd[invcpos] = psw.xcrd[cdk.mol_contents[j]];
      invcpos++;
      inverted_crd[invcpos] = psw.ycrd[cdk.mol_contents[j]];
      invcpos++;
      inverted_crd[invcpos] = psw.zcrd[cdk.mol_contents[j]];
      invcpos++;
    }
    
    // Record the bond and angle energies
    const double new_bond_e = evaluateBondTerms(vk, CoordinateFrameReader(inversion_copy), &sc, 0);
    ubond_dev[i] = fabs(new_bond_e - orig_bond_e);
    
    // Reverse the inversion to check that the molecule recovers its initial state
    flipChiralCenter(&inversion_copy, i, centers, protocols, inv_grp);
    repos_dev[i] = rmsd(ps, inversion_copy, ag, RMSDMethod::ALIGN_GEOM, cdk.mol_limits[0],
                        cdk.mol_limits[1]);
  }
  const char osc = osSeparator();
  const std::string base_iso_path = oe.getStormmSourcePath() + osc + "test" + osc + "Structure";
  const std::string invcrd_snapshot = base_iso_path + osc + "inverted_coords.m";
  const bool snps_exist = (getDrivePathType(invcrd_snapshot) == DrivePathType::FILE);
  if (snps_exist == false && oe.takeSnapshot() == SnapshotOperation::COMPARE) {
    rtWarn("The snapshot file " + invcrd_snapshot + " was not found.  Check the ${STORMM_SOURCE} "
           "environment variable, currently set to " + oe.getStormmSourcePath() +
           ", for validity.  Subsequent tests will be skipped.", "test_isomerization");
  }
  const TestPriority do_snps = (snps_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (snp_var_name.size() > 0) {
    snapshot(invcrd_snapshot, polyNumericVector(inverted_crd), snp_var_name, 1.0e-6, "Coordinates "
             "obtained from inversion of selected chiral centers do not meet expectations.",
             oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL, expectation, do_snps);
  }
  const Approx target(std::vector<double>(nchiral, 0.0), ComparisonType::ABSOLUTE, 1.0e-6);
  check(repos_dev, RelationalOperator::EQUAL, target, "Reversing the inversion operations of "
        "selected chiral centers does not return the molecule to its original state.", do_tests);
  check(ubond_dev, RelationalOperator::EQUAL, target, "Bond energies are changed by inversion "
        "of a chiral center.", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Run a series of tests comparing various structural manipulations across different coordinate
// objects.
//
// Arguments:
//   cfv:         An array of coordinate sets
//   psv:         Array of PhaseSpace objects, velocities and forces being irrelevant
//   hr_poly_ps:  High-resolution synthesis, a collection of systems matching cfv and psv
//   lr_poly_ps:  Low-resolution synthesis, similar to hr_poly_ps
//   cdns_hr:     Condensate based on one of the syntheses
//   tsm:         Cache of test systems, provided to produce a signal about the testing status
//   effect:      Descriptor of the global coordinate change operation
//   lp_mue_tol:  Mean unsigned error tolerance for manipulations in low-precision syntheses
//   lp_abs_tol:  Absolute tolerance for manipulations in low-precision synthesis representations
//-------------------------------------------------------------------------------------------------
void testStructuralCorrespondence(const std::vector<CoordinateFrame> &cfv,
                                  const std::vector<PhaseSpace> &psv,
                                  const PhaseSpaceSynthesis &hr_poly_ps,
                                  const PhaseSpaceSynthesis &lr_poly_ps,
                                  const Condensate &cdns_hr, const TestSystemManager &tsm,
                                  const std::string &effect, const double lp_mue_tol = 1.0e-6,
                                  const double lp_abs_tol = 1.0e-5) {
  const int ntrials = cfv.size();
  std::vector<double> cf_ps_mue(ntrials);
  std::vector<double> cf_ps_abs(ntrials);
  std::vector<double> cf_hsyn_mue(ntrials);
  std::vector<double> cf_hsyn_abs(ntrials);
  std::vector<double> cf_lsyn_mue(ntrials);
  std::vector<double> cf_lsyn_abs(ntrials);
  std::vector<double> cf_hcdns_mue(ntrials);
  std::vector<double> cf_hcdns_abs(ntrials);
  for (int i = 0; i < ntrials; i++) {
    const std::vector<double> cf_xyz = cfv[i].getInterlacedCoordinates();
    const std::vector<double> ps_xyz = cfv[i].getInterlacedCoordinates();
    const std::vector<double> hpoly_ps_xyz = hr_poly_ps.getInterlacedCoordinates(i);
    const std::vector<double> lpoly_ps_xyz = lr_poly_ps.getInterlacedCoordinates(i);
    const std::vector<double> hcdns_xyz = cdns_hr.getInterlacedCoordinates(i);
    cf_ps_mue[i] = meanUnsignedError(cf_xyz, ps_xyz);
    cf_ps_abs[i] = maxAbsoluteDifference(cf_xyz, ps_xyz);
    cf_hsyn_mue[i] = meanUnsignedError(cf_xyz, hpoly_ps_xyz);
    cf_hsyn_abs[i] = maxAbsoluteDifference(cf_xyz, hpoly_ps_xyz);
    cf_lsyn_mue[i] = meanUnsignedError(cf_xyz, lpoly_ps_xyz);
    cf_lsyn_abs[i] = maxAbsoluteDifference(cf_xyz, lpoly_ps_xyz);
    cf_hcdns_mue[i] = meanUnsignedError(cf_xyz, hcdns_xyz);
    cf_hcdns_abs[i] = maxAbsoluteDifference(cf_xyz, hcdns_xyz);
  }
  const Approx zero_target(std::vector<double>(ntrials, 0.0));
  check(cf_ps_mue, RelationalOperator::EQUAL, zero_target.margin(tiny), "The mean unsigned "
        "error due to " + effect + " in the CoordinateFrame and PhaseSpace objects is "
        "unexpectedly large.", tsm.getTestingStatus());
  check(cf_ps_abs, RelationalOperator::EQUAL, zero_target.margin(1.0e-8), "The maximum "
        "deviation due to " + effect + " in the CoordinateFrame and PhaseSpace objects is "
        "unexpectedly large.", tsm.getTestingStatus());
  check(cf_hsyn_mue, RelationalOperator::EQUAL, zero_target.margin(tiny), "The mean unsigned "
        "error due to " + effect + " in the CoordinateFrame and high-precision "
        "PhaseSpaceSynthesis objects is unexpectedly large.", tsm.getTestingStatus());
  check(cf_hsyn_abs, RelationalOperator::EQUAL, zero_target.margin(1.0e-8), "The maximum "
        "deviation due to " + effect + " in the CoordinateFrame and high-precision "
        "PhaseSpaceSynthesis objects is unexpectedly large.", tsm.getTestingStatus());
  check(cf_lsyn_mue, RelationalOperator::EQUAL, zero_target.margin(lp_mue_tol), "The mean "
        "unsigned error due to " + effect + " in the CoordinateFrame and low-precision "
        "PhaseSpaceSynthesis objects is unexpectedly large.", tsm.getTestingStatus());
  check(cf_lsyn_abs, RelationalOperator::EQUAL, zero_target.margin(lp_abs_tol), "The maximum "
        "deviation due to " + effect + " in the CoordinateFrame and low-precision "
        "PhaseSpaceSynthesis objects is unexpectedly large.", tsm.getTestingStatus());
  check(cf_hcdns_mue, RelationalOperator::EQUAL, zero_target.margin(tiny), "The mean unsigned "
        "error due to " + effect + " in the CoordinateFrame and high-precision Condensate "
        "objects is unexpectedly large.", tsm.getTestingStatus());
  check(cf_hcdns_abs, RelationalOperator::EQUAL, zero_target.margin(1.0e-8), "The maximum "
        "deviation due to " + effect + " in the CoordinateFrame and high-precision Condensate "
        "objects is unexpectedly large.", tsm.getTestingStatus());
}

//-------------------------------------------------------------------------------------------------
// Check rotational, translational, and chiral sampling operations using a synthesis of structures.
// In the functions above, operations on individual frames are tested in typical double-precision
// arithmetic.  Here, the various operations are performed on each coordinate representation and
// compared to one another.
//
// Arguments:
//   tsm:         Manager object for multiple test systems
//   igseed:      Random number generator seed value
//   ntrials:     The number of different systems to attempt perturbing, and to place in each
//                synthesis object
//   nrot_perm:   Number of bond rotations to perform on each structure with such features
//   nctx_perm:   Number of cis-trans flips to perform on each structure with such features
//   ninv_perm:   Number of chiral inversions to perform on each structure with such features
//   lp_mue_tol:  Mean unsigned error tolerance for manipulations in low-precision syntheses
//   lp_abs_tol:  Absolute tolerance for manipulations in low-precision synthesis representations
//-------------------------------------------------------------------------------------------------
void checkStructureManipulation(const TestSystemManager &tsm, const int igseed,
                                const int ntrials = 8, const int nrot_perm = 16,
                                const int nctx_perm = 4, const int ninv_perm = 1,
                                const double lp_mue_tol = 1.0e-6,
                                const double lp_abs_tol = 1.0e-5) {

  // Boot up a random number generator.
  Xoshiro256ppGenerator xrs(igseed);
  
  // Extract random systems and make single-system, single-structure coordinate objects.
  std::vector<int> sequence(ntrials);
  for (int i = 0; i < ntrials; i++) {
    sequence[i] = xrs.uniformRandomNumber() * tsm.getSystemCount();
  }
  std::vector<CoordinateFrame> cfv;
  std::vector<PhaseSpace> psv;
  std::vector<AtomGraph*> agv(ntrials);
  std::vector<ChemicalFeatures> chemv;
  cfv.reserve(ntrials);
  psv.reserve(ntrials);
  chemv.reserve(ntrials);
  for (int i = 0; i < ntrials; i++) {
    cfv.push_back(tsm.exportCoordinateFrame(sequence[i]));
    psv.push_back(tsm.exportPhaseSpace(sequence[i]));
    agv[i] = const_cast<AtomGraph*>(tsm.getTopologyPointer(sequence[i]));
    chemv.emplace_back(agv[i], cfv[i], MapRotatableGroups::YES);
  }
  std::vector<CoordinateFrame> w_cfv = cfv;
  std::vector<PhaseSpace> w_psv = psv;
  
  // Prepare a coordinate series for a random system and produce a condensate based on that.
  const int series_no = xrs.uniformRandomNumber() * tsm.getSystemCount();
  CoordinateSeries<float> cs(tsm.exportCoordinateFrame(series_no), 4);
  Condensate cdns_cs(cs, PrecisionModel::SINGLE);
  const AtomGraph& ag_for_cs = tsm.getTopologyReference(series_no);
  
  // Prepare a high-resolution synthesis and a condensate based on that.
  PhaseSpaceSynthesis  hr_poly_ps = tsm.exportPhaseSpaceSynthesis(sequence, 0.0, 10983, 48);
  PhaseSpaceSynthesis whr_poly_ps = tsm.exportPhaseSpaceSynthesis(sequence, 0.0, 10983, 48);
  Condensate cdns_hr(hr_poly_ps, PrecisionModel::DOUBLE);
  Condensate cdns_whr(whr_poly_ps, PrecisionModel::DOUBLE);

  // Prepare a low-resolution synthesis.
  PhaseSpaceSynthesis  lr_poly_ps = tsm.exportPhaseSpaceSynthesis(sequence, 0.0, 10983, 24);
  PhaseSpaceSynthesis wlr_poly_ps = tsm.exportPhaseSpaceSynthesis(sequence, 0.0, 10983, 24);

  // Select a sequence of bond rotations, cis-trans flips, and chiral inversions for each
  // structure.  Apply them to each individual frame, or a randomly selected frame from each
  // each synthesis object, then compare the results.
  std::vector<int> rotg_counts(ntrials);
  std::vector<int> ctxg_counts(ntrials);
  std::vector<int> invg_counts(ntrials);

  // Perform bond rotations on each structure.
  for (int i = 0; i < ntrials; i++) {
    const int n_rotg = chemv[i].getRotatableBondCount();
    std::vector<int> rot_groups(nrot_perm);
    std::vector<double> rot_changes(nrot_perm);
    rotg_counts[i] = n_rotg;
    if (n_rotg > 0) {
      for (int j = 0; j < nrot_perm; j++) {
        rot_groups[j] = xrs.uniformRandomNumber() * static_cast<double>(n_rotg);
        rot_changes[j] = xrs.gaussianRandomNumber() * twopi;
      }
      const std::vector<IsomerPlan> rotg_moves = chemv[i].getRotatableBondGroups();
      for (int j = 0; j < nrot_perm; j++) {
        const size_t rgj = rot_groups[j];
        rotateAboutBond(&w_cfv[i], rotg_moves[rgj].getRootAtom(), rotg_moves[rgj].getPivotAtom(),
                        rotg_moves[rgj].getMovingAtoms(), rot_changes[j]);
        rotateAboutBond(&w_psv[i], rotg_moves[rgj].getRootAtom(), rotg_moves[rgj].getPivotAtom(),
                        rotg_moves[rgj].getMovingAtoms(), rot_changes[j]);
        rotateAboutBond<double>(&whr_poly_ps, i, rotg_moves[rgj].getRootAtom(),
                                rotg_moves[rgj].getPivotAtom(), rotg_moves[rgj].getMovingAtoms(),
                                rot_changes[j]);
        rotateAboutBond<float>(&wlr_poly_ps, i, rotg_moves[rgj].getRootAtom(),
                               rotg_moves[rgj].getPivotAtom(), rotg_moves[rgj].getMovingAtoms(),
                               rot_changes[j]);
        rotateAboutBond<double>(&cdns_whr, i, rotg_moves[rgj].getRootAtom(),
                                rotg_moves[rgj].getPivotAtom(), rotg_moves[rgj].getMovingAtoms(),
                                rot_changes[j]);
      }
    }
  }
  if (sum<int>(rotg_counts) > 0) {
    testStructuralCorrespondence(w_cfv, w_psv, whr_poly_ps, wlr_poly_ps, cdns_whr, tsm,
                                 "bond rotation", lp_mue_tol, lp_abs_tol);
  }

  // Reset the writeable coordinate objects to prevent errors in the previous tests from cascading
  // into subsequent tests.
  for (int i = 0; i < ntrials; i++) {
    coordCopy(&w_cfv[i], cfv[i]);
    coordCopy(&w_psv[i], psv[i]);
    coordCopy(&whr_poly_ps, i, hr_poly_ps, i);
    coordCopy(&wlr_poly_ps, i, lr_poly_ps, i);
    coordCopy(&cdns_whr, i, cdns_hr, i);
  }

  // Perform cis-transisomerizations, with slight perturbations to prevent two flips from returning
  // a structure to its starting point.
  for (int i = 0; i < ntrials; i++) {
    const int n_ctxg = chemv[i].getCisTransBondCount();
    std::vector<int> ctx_groups(nctx_perm);
    std::vector<double> ctx_changes(nctx_perm);  
    ctxg_counts[i] = n_ctxg;
    if (n_ctxg > 0) {
      for (int j = 0; j < nctx_perm; j++) {
        ctx_groups[j] = xrs.uniformRandomNumber() * static_cast<double>(n_ctxg);
        ctx_changes[j] = (round(xrs.uniformRandomNumber() * 3.0) +
                          (xrs.uniformRandomNumber() * 0.1)) * pi;
      }
      const std::vector<IsomerPlan> ctxg_moves = chemv[i].getCisTransIsomerizationGroups();
      for (int j = 0; j < nctx_perm; j++) {
        const size_t cgj = ctx_groups[j];
        rotateAboutBond(&w_cfv[i], ctxg_moves[cgj].getRootAtom(), ctxg_moves[cgj].getPivotAtom(),
                        ctxg_moves[cgj].getMovingAtoms(), ctx_changes[j]);
        rotateAboutBond(&w_psv[i], ctxg_moves[cgj].getRootAtom(), ctxg_moves[cgj].getPivotAtom(),
                        ctxg_moves[cgj].getMovingAtoms(), ctx_changes[j]);
        rotateAboutBond<double>(&whr_poly_ps, i, ctxg_moves[cgj].getRootAtom(),
                                ctxg_moves[cgj].getPivotAtom(), ctxg_moves[cgj].getMovingAtoms(),
                                ctx_changes[j]);
        rotateAboutBond<float>(&wlr_poly_ps, i, ctxg_moves[cgj].getRootAtom(),
                               ctxg_moves[cgj].getPivotAtom(), ctxg_moves[cgj].getMovingAtoms(),
                               ctx_changes[j]);
        rotateAboutBond<double>(&cdns_whr, i, ctxg_moves[cgj].getRootAtom(),
                                ctxg_moves[cgj].getPivotAtom(), ctxg_moves[cgj].getMovingAtoms(),
                                ctx_changes[j]);        
      }
    }
  }
  if (sum<int>(ctxg_counts) > 0) {
    testStructuralCorrespondence(w_cfv, w_psv, whr_poly_ps, wlr_poly_ps, cdns_whr, tsm,
                                 "cis-trans isomerization", lp_mue_tol, lp_abs_tol);
  }

  // Reset the writeable coordinate objects to prevent errors in the previous tests from cascading
  // into subsequent tests.
  for (int i = 0; i < ntrials; i++) {
    coordCopy(&w_cfv[i], cfv[i]);
    coordCopy(&w_psv[i], psv[i]);
    coordCopy(&whr_poly_ps, i, hr_poly_ps, i);
    coordCopy(&wlr_poly_ps, i, lr_poly_ps, i);
    coordCopy(&cdns_whr, i, cdns_hr, i);
  }

  // Perform chiral inversions.
  for (int i = 0; i < ntrials; i++) {
    const int n_invg = chemv[i].getChiralCenterCount();
    std::vector<int> inv_groups(ninv_perm);
    std::vector<ChiralOrientation> inv_changes(ninv_perm);
    invg_counts[i] = n_invg;
    if (n_invg > 0) {
      for (int j = 0; j < ninv_perm; j++) {
        inv_groups[j] = xrs.uniformRandomNumber() * static_cast<double>(n_invg);
        switch (chemv[i].getChiralInversionMethods(inv_groups[j])) {
        case ChiralInversionProtocol::ROTATE:
        case ChiralInversionProtocol::REFLECT:
          inv_changes[j] = (xrs.uniformRandomNumber() < 0.5) ? ChiralOrientation::RECTUS :
                                                               ChiralOrientation::SINISTER;
          break;
        case ChiralInversionProtocol::DO_NOT_INVERT:

          // Use the non-assigned value here to denote the lack of an instruction for changing
          // the chirality of a group that indeed has chirality but is not supposed to change.
          inv_changes[j] = ChiralOrientation::NONE;
          break;
        }
      }
      const std::vector<IsomerPlan> invg_moves = chemv[i].getChiralInversionGroups();
      const std::vector<int> invg_centers = chemv[i].getChiralCenters();
      const std::vector<ChiralInversionProtocol> invg_protocols =
        chemv[i].getChiralInversionMethods();
      for (int j = 0; j < ninv_perm; j++) {
        flipChiralCenter(&w_cfv[i], inv_groups[j], invg_centers, invg_protocols, invg_moves);
        flipChiralCenter(&w_psv[i], inv_groups[j], invg_centers, invg_protocols, invg_moves);
        flipChiralCenter<double>(&whr_poly_ps, i, inv_groups[j], invg_centers, invg_protocols,
                                 invg_moves);
        flipChiralCenter<float>(&wlr_poly_ps, i, inv_groups[j], invg_centers, invg_protocols,
                                invg_moves);
        flipChiralCenter<double>(&cdns_whr, i, inv_groups[j], invg_centers, invg_protocols,
                                 invg_moves);
      }
    }
  }
  if (sum<int>(invg_counts) > 0) {
    testStructuralCorrespondence(w_cfv, w_psv, whr_poly_ps, wlr_poly_ps, cdns_whr, tsm,
                                 "chiral inversion", lp_mue_tol, lp_abs_tol);
  }
}

//-------------------------------------------------------------------------------------------------
// Test the RMSD calculation guide on a variety of structures.  This will takes lists of system
// topologies and coordinates, replicate each to varying degrees, and apply a small perturbation
// plus rotational and translational motion to make RMSD calculations and alignments meaningful.
//
// Arguments:
//   ag_list:              List of system topologies describing each structure in ps_list
//   ps_list:              List of initial structures (each will be replicated, then randomly
//                         perturbed, rotated, and translated)
//   approach:             The RMSD calculation method to use
//   do_tests:             Pre-determined indicator of whether unit tests are possible
//   frame_counts:         The number of replicas to apply to each system
//   default_frame_count:  The number of replicas to apply to any remaining systems if the
//                         frame_counts vector runs out of information
//-------------------------------------------------------------------------------------------------
void testRMSDGuide(const std::vector<AtomGraph> &ag_list, const std::vector<PhaseSpace> &ps_list,
                   const RMSDMethod approach, const TestPriority do_tests,
                   const std::vector<int> &frame_counts = { 16, 14, 8, 7, 19 },
                   const int default_frame_count = 16) {
  std::vector<CoordinateSeries<double>> structure_pile;
  const SystemGrouping top_grp = SystemGrouping::TOPOLOGY;
  for (size_t item = 0; item < ag_list.size(); item++) {
    const CoordinateFrame item_cf(ps_list[item]);
    const AtomEquivalence item_eq(ag_list[item], item_cf);
    RMSDPlan item_rplan(item_eq, approach);
    const int nframe = (item < frame_counts.size()) ? frame_counts[item] : default_frame_count;
    CoordinateSeries<double> item_series(ps_list[item], nframe);
    CoordinateSeriesWriter<double> item_seriesw = item_series.data();
    Xoshiro256ppGenerator xrs(71842203 + item);
    for (int i = 1; i < item_seriesw.nframe; i++) {
      const int frame_offset = roundUp(item_seriesw.natom, warp_size_int) * i;
      for (int j = 0; j < item_seriesw.natom; j++) {
        const int ij_pos = frame_offset + j;
        item_seriesw.xcrd[ij_pos] += 0.25 * (0.5 - xrs.uniformRandomNumber());
        item_seriesw.ycrd[ij_pos] += 0.25 * (0.5 - xrs.uniformRandomNumber());
        item_seriesw.zcrd[ij_pos] += 0.25 * (0.5 - xrs.uniformRandomNumber());
      }
    }
    ComparisonGuide item_cg(item_series, ag_list[item]);
    Hybrid<double> item_rmsds(item_cg.getAllToReferenceOutputSize());
    Hybrid<double> item_pair_rmsds(item_cg.getSymmetryEquivalentPairOutputSize(top_grp));
    rmsd(item_cg, item_rplan, item_cf, item_series, &item_rmsds);
    rmsd(item_cg, item_rplan, item_series, &item_pair_rmsds);

    // Check basic features of the RMSDPlan object
    if (item == 0 && ps_list.size() == 1) {
      RMSDPlanReader<double> item_rpr = item_rplan.dpData();
      check(item_rpr.plan_count, RelationalOperator::EQUAL, 1, "Only one plan should be present "
            "in the RMSD guide tailored for a single system.", do_tests);
    }

    // Remember the structures
    structure_pile.push_back(item_series);
  }

  // Shuffle the structures (but do not re-order the list of structures instantiating any given
  // topology) and create a synthesis of them.
  std::vector<PhaseSpace> shuffled_ps;
  std::vector<AtomGraph*> shuffled_ag;
  int ns = 0;
  for (size_t i = 0; i < ag_list.size(); i++) {
    ns += structure_pile[i].getFrameCount();
  }
  shuffled_ps.reserve(ns);
  shuffled_ag.reserve(ns);
  int added = 1;
  int add_frame = 0;
  while (added > 0) {
    added = 0;
    for (size_t i = 0; i < ag_list.size(); i++) {
      if (add_frame < structure_pile[i].getFrameCount()) {
        const CoordinateFrame cf = structure_pile[i].exportFrame(add_frame);
        const CoordinateFrameReader cfr = cf.data();
        shuffled_ps.emplace_back(cfr.natom, cfr.unit_cell);
        shuffled_ps.back().fill(cfr.xcrd, cfr.ycrd, cfr.zcrd, TrajectoryKind::POSITIONS,
                                CoordinateCycle::WHITE, 0, cfr.boxdim);
        shuffled_ag.push_back(const_cast<AtomGraph*>(ag_list[i].getSelfPointer()));
        added++;
      }
    }
    add_frame++;
  }
  PhaseSpaceSynthesis poly_ps(shuffled_ps, shuffled_ag, 40);
  const RMSDPlan poly_rplan(poly_ps, approach);
  const ComparisonGuide poly_cg(poly_ps);
  Hybrid<double> poly_rmsds(poly_cg.getAllToReferenceOutputSize());
  Hybrid<double> poly_pair_rmsds(poly_cg.getSymmetryEquivalentPairOutputSize(top_grp));
  Hybrid<int> sys_example_idx(poly_ps.getUniqueTopologyExampleIndices());
  rmsd(poly_cg, poly_rplan, poly_ps, sys_example_idx, &poly_rmsds);
  rmsd(poly_cg, poly_rplan, poly_ps, &poly_pair_rmsds);
  std::vector<double> replica_rmsds, replica_pair_rmsds;
  for (size_t i = 0; i < ag_list.size(); i++) {
    const ChemicalDetailsKit icdk = ag_list[i].getChemicalDetailsKit();

    // Recompute the reference RMSDs
    replica_rmsds.push_back(0.0);
    const int jlim = structure_pile[i].getFrameCount();
    const CoordinateFrame iref = structure_pile[i].exportFrame(0);
    const CoordinateFrameReader irefr = iref.data();
    for (int j = 1; j < jlim; j++) {
      const CoordinateFrame cf_ij = structure_pile[i].exportFrame(j);
      replica_rmsds.push_back(rmsd(irefr, cf_ij.data(), icdk, approach));
    }

    // Recompute the RMSD matrices
    for (size_t j = 1; j < jlim; j++) {
      const CoordinateFrame cf_ij = structure_pile[i].exportFrame(j);
      const CoordinateFrameReader cf_ijr = cf_ij.data();
      for (size_t k = 0; k < j; k++) {
        const CoordinateFrame cf_ik = structure_pile[i].exportFrame(k);
        replica_pair_rmsds.push_back(rmsd(cf_ijr, cf_ik.data(), icdk, approach));
      }
    }
    const size_t padded_triangle = roundUp(replica_pair_rmsds.size(), warp_size_zu);
    for (size_t j = replica_pair_rmsds.size(); j < padded_triangle; j++) {
      replica_pair_rmsds.push_back(0.0);
    }
  }
  check(poly_rmsds.readHost(), RelationalOperator::EQUAL, Approx(replica_rmsds).margin(tiny),
        "RMSD calculations to a reference structure performed across a synthesis of structures "
        "did not match the list of values computed individually.", do_tests);
  check(poly_pair_rmsds.readHost(), RelationalOperator::EQUAL,
        Approx(replica_pair_rmsds).margin(tiny), "RMSD matrix calculations to a reference "
        "structure performed across a synthesis of structures did not match the list of values "
        "computed individually.", do_tests);
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
  section("Rotation about a bond");

  // Section 2
  section("Selected chiral inversions");

  // Section 3
  section("RMSD calculations");

  // Section 4
  section("Symmetry-reduced RMSD");

  // Section 5
  section("Clash detection");
  
  // Get a handful of realistic systems
  const char osc = osSeparator();
  const std::string base_top_path = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_path = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string drug_top_path = base_top_path + osc + "drug_example.top";
  const std::string drug_crd_path = base_crd_path + osc + "drug_example.inpcrd";
  const std::string lig1_top_path = base_top_path + osc + "stereo_L1.top";
  const std::string lig1_crd_path = base_crd_path + osc + "stereo_L1.inpcrd";
  const std::string lig2_top_path = base_top_path + osc + "symmetry_L1.top";
  const std::string lig2_crd_path = base_crd_path + osc + "symmetry_L1.inpcrd";
  const std::string trpc_top_path = base_top_path + osc + "trpcage.top";
  const std::string trpc_crd_path = base_crd_path + osc + "trpcage.inpcrd";
  const bool files_exist = (getDrivePathType(drug_top_path) == DrivePathType::FILE &&
                            getDrivePathType(drug_crd_path) == DrivePathType::FILE &&
                            getDrivePathType(lig1_top_path) == DrivePathType::FILE &&
                            getDrivePathType(lig1_crd_path) == DrivePathType::FILE &&
                            getDrivePathType(lig2_top_path) == DrivePathType::FILE &&
                            getDrivePathType(lig2_crd_path) == DrivePathType::FILE &&
                            getDrivePathType(trpc_top_path) == DrivePathType::FILE &&
                            getDrivePathType(trpc_crd_path) == DrivePathType::FILE);
  if (files_exist == false) {
    rtWarn("Files for various drug molecules and a miniprotein were not found.  Check the "
           "$STORMM_SOURCE environment variable to ensure that " + drug_top_path + " and " +
           drug_crd_path + " become valid paths.  Some tests will be skipped",
           "test_local_arrangement");
  }
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  AtomGraph drug_ag = (files_exist) ? AtomGraph(drug_top_path) : AtomGraph();
  PhaseSpace drug_ps = (files_exist) ? PhaseSpace(drug_crd_path,
                                                  CoordinateFileKind::AMBER_INPCRD) :
                                       PhaseSpace();
  const ChemicalFeatures drug_feat = (files_exist) ? ChemicalFeatures(&drug_ag, drug_ps,
                                                                      MapRotatableGroups::YES) :
                                                     ChemicalFeatures();
  AtomGraph trpc_ag = (files_exist) ? AtomGraph(trpc_top_path) : AtomGraph();
  PhaseSpace trpc_ps = (files_exist) ? PhaseSpace(trpc_crd_path,
                                                  CoordinateFileKind::AMBER_INPCRD) :
                                       PhaseSpace();
  const ChemicalFeatures trpc_feat = (files_exist) ? ChemicalFeatures(&trpc_ag, trpc_ps,
                                                                      MapRotatableGroups::YES) :
                                                     ChemicalFeatures();
  AtomGraph lig1_ag = (files_exist) ? AtomGraph(lig1_top_path) : AtomGraph();
  PhaseSpace lig1_ps = (files_exist) ? PhaseSpace(lig1_crd_path,
                                                  CoordinateFileKind::AMBER_INPCRD) :
                                       PhaseSpace();
  const ChemicalFeatures lig1_feat = (files_exist) ? ChemicalFeatures(&lig1_ag, lig1_ps,
                                                                      MapRotatableGroups::YES) :
                                                     ChemicalFeatures();
  AtomGraph lig2_ag = (files_exist) ? AtomGraph(lig2_top_path) : AtomGraph();
  PhaseSpace lig2_ps = (files_exist) ? PhaseSpace(lig2_crd_path,
                                                  CoordinateFileKind::AMBER_INPCRD) :
                                       PhaseSpace();
  const ChemicalFeatures lig2_feat = (files_exist) ? ChemicalFeatures(&lig2_ag, lig2_ps,
                                                                      MapRotatableGroups::YES) :
                                                     ChemicalFeatures();
  
  // Rotate bonds within each system
  checkRotationalSampling(drug_ag, drug_ps, drug_feat, oe, do_tests, "drug_rot_iso",
                          PrintSituation::OVERWRITE);
  checkRotationalSampling(trpc_ag, trpc_ps, trpc_feat, oe, do_tests);
  checkRotationalSampling(lig1_ag, lig1_ps, lig1_feat, oe, do_tests, "lig1_rot_iso");
  checkRotationalSampling(lig2_ag, lig2_ps, lig2_feat, oe, do_tests, "lig2_rot_iso");

  // Isomerize chiral centers
  checkChiralSampling(drug_ag, drug_ps, drug_feat, oe, do_tests, "drug_chir_iso",
                      PrintSituation::OVERWRITE);
  checkChiralSampling(trpc_ag, trpc_ps, trpc_feat, oe, do_tests);
  checkChiralSampling(lig1_ag, lig1_ps, lig1_feat, oe, do_tests, "lig1_chir_iso");
  checkChiralSampling(lig2_ag, lig2_ps, lig2_feat, oe, do_tests, "lig2_chir_iso");

  // Check global structural manipulations on all coordinate objects.
  const std::string pro_top_path = oe.getStormmSourcePath() + osc + "test" + osc + "Namelists" +
                                   osc + "topol";
  const std::string pro_crd_path = oe.getStormmSourcePath() + osc + "test" + osc + "Namelists" +
                                   osc + "coord";
  const std::vector<std::string> pro_names = { "gly_phe", "tyr", "gly_lys", "ala", "gly_pro",
                                               "arg", "gly_arg" };
  TestSystemManager pro_tsm(pro_top_path, "top", pro_names, pro_crd_path, "inpcrd", pro_names);
  checkStructureManipulation(pro_tsm, 39583203, 12);
  const std::vector<std::string> drug_names = { "med_1", "symmetry_L1_vs", "drug_example_iso",
                                                "med_2", "stereo_L1", "stereo_L1_vs", "med_4",
                                                "med_3", "bromobenzene_iso", "symmetry_C2",
                                                "symmetry_C3", "symmetry_C6" };
  TestSystemManager drug_tsm(base_top_path, "top", drug_names, base_crd_path, "inpcrd",
                             drug_names);
  checkStructureManipulation(drug_tsm, 39583203, 17, 16, 4, 1, 3.0e-5, 1.2e-4);

  // Test RMSD computations on simple structures
  section(3);
  CoordinateFrame starfish_a(7);
  CoordinateFrameWriter strfa = starfish_a.data();
  strfa.xcrd[1] =  1.0;
  strfa.ycrd[2] =  1.0;
  strfa.xcrd[3] = -1.0;
  strfa.ycrd[4] = -1.0;
  strfa.zcrd[5] =  1.0;
  strfa.zcrd[6] = -1.0;
  CoordinateFrame starfish_b(starfish_a);
  CoordinateFrameWriter strfb = starfish_b.data();
  rotateCoordinates(strfb.xcrd, strfb.ycrd, strfb.zcrd, 0.0, 0.0, stormm::symbols::pi / 4.0, 0,
                    strfb.natom);
  double rms_no_align = rmsd<double, double>(strfa.xcrd, strfa.ycrd, strfa.zcrd, strfb.xcrd,
                                             strfb.ycrd, strfb.zcrd, nullptr,
                                             RMSDMethod::NO_ALIGN_GEOM, 0, strfa.natom);
  double rms_align = rmsd<double, double>(strfa.xcrd, strfa.ycrd, strfa.zcrd, strfb.xcrd,
                                          strfb.ycrd, strfb.zcrd, nullptr, RMSDMethod::ALIGN_GEOM,
                                          0, strfa.natom);
  check(rms_no_align, RelationalOperator::EQUAL, Approx(0.578562967).margin(1.0e-8),
        "Positional (non-aligned) RMSD computed for coordinates pre-shifted to their respective "
        "centers of mass does not produce the expected result.");
  check(rms_align, RelationalOperator::EQUAL, Approx(0.0).margin(1.0e-8), "Quaternion-aligned, "
        "positional RMSD computed for coordinates pre-shifted to their respective centers of mass "
        "does not produce the expected result.");
  translateCoordinates(strfb.xcrd, strfb.ycrd, strfb.zcrd, 5.0, 4.8, 9.7, 0, strfb.natom);
  rms_no_align = rmsd<double, double>(strfa.xcrd, strfa.ycrd, strfa.zcrd, strfb.xcrd, strfb.ycrd,
                                      strfb.zcrd, nullptr, RMSDMethod::NO_ALIGN_GEOM, 0,
                                      strfa.natom);
  rms_align = rmsd<double, double>(strfa.xcrd, strfa.ycrd, strfa.zcrd, strfb.xcrd, strfb.ycrd,
                                   strfb.zcrd, nullptr, RMSDMethod::ALIGN_GEOM, 0, strfa.natom);
  check(rms_no_align, RelationalOperator::EQUAL, Approx(11.935859211).margin(1.0e-8),
        "Positional (non-aligned) RMSD computed for coordinates differing in their respective "
        "centers of mass does not produce the expected result.");
  check(rms_align, RelationalOperator::EQUAL, Approx(0.0).margin(1.0e-8), "Quaternion-aligned, "
        "positional RMSD computed for coordinates differing in their respective centers of mass "
        "does not produce the expected result.");
  rotateCoordinates(strfb.xcrd, strfb.ycrd, strfb.zcrd, 0.1, -0.3, -0.25, 0, strfb.natom);
  rms_align = rmsd<double, double>(strfa.xcrd, strfa.ycrd, strfa.zcrd, strfb.xcrd, strfb.ycrd,
                                   strfb.zcrd, nullptr, RMSDMethod::ALIGN_GEOM, 0, strfa.natom);
  check(rms_align, RelationalOperator::EQUAL, Approx(0.0).margin(1.0e-8), "Quaternion-aligned, "
        "positional RMSD computed for coordinates differing in their respective centers of mass, "
        "one rotated a second time for more frustration, does not produce the expected result.");

  // Test an RMSD calculation guide
  section(4);
  testRMSDGuide({ lig2_ag }, { lig2_ps }, RMSDMethod::ALIGN_MASS, do_tests);
  testRMSDGuide({ lig1_ag, lig2_ag, trpc_ag, }, { lig1_ps, lig2_ps, trpc_ps },
                RMSDMethod::ALIGN_MASS, do_tests);

  // Test clash detection
  section(5);
  const std::string example_top_path = oe.getStormmSourcePath() + osc + "test" + osc +
                                       "Namelists" + osc + "topol";
  const std::string example_crd_path = oe.getStormmSourcePath() + osc + "test" + osc +
                                       "Namelists" + osc + "coord";
  const std::vector<std::string> sysnames = { "arg", "gly_phe", "lys" };
  TestSystemManager clashers(example_top_path, "top", sysnames, example_crd_path, "inpcrd",
                             sysnames);
  std::vector<ClashReport> impacts;
  std::vector<int> busted;
  const std::vector<int> busted_ans = { 1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0 };
  bool all_ok = true;
  Xoshiro256ppGenerator xrs(15923);
  for (int rpos = 0; rpos < 36; rpos++) {
    const double alpha = twopi * (xrs.uniformRandomNumber() - 0.5);
    const double beta  = twopi * (xrs.uniformRandomNumber() - 0.5);
    const double gamma = twopi * (xrs.uniformRandomNumber() - 0.5);
    for (size_t i = 0; i < sysnames.size(); i++) {
      const std::string clash_crd_name = oe.getStormmSourcePath() + osc + "test" + osc +
                                         "Structure" + osc + "clash_" + sysnames[i] + ".crd";
      const bool clsh_file_exists = (getDrivePathType(clash_crd_name) == DrivePathType::FILE);
      const TestPriority do_clshi = (clsh_file_exists) ? TestPriority::CRITICAL :
                                                         TestPriority::ABORT;
      const AtomGraph iag = clashers.getTopologyReference(i);
      const StaticExclusionMask imask(iag);
      CoordinateSeries<double> ics(iag.getAtomCount());
      if (clsh_file_exists) {
        ics.importFromFile(clash_crd_name);
      }
      const int nfrm = ics.getFrameCount();
      for (int j = 0; j < nfrm; j++) {
        impacts.emplace_back();
        busted.push_back(detectClash(ics, j, iag, imask, default_minimize_clash_r0,
                                     default_minimize_clash_ratio, &impacts[impacts.size() - 1]));

        // Rotate coordinates in preparation for another round.
        rotateCoordinates(&ics, j, iag, alpha, beta, gamma);
      }
    }
    if (rpos == 0) {
      check(busted, RelationalOperator::EQUAL, busted_ans, "Clashes are not detected where "
            "expected in various coordiante series.", clashers.getTestingStatus());
    }
    else {
      for (size_t i = 0; i < busted_ans.size(); i++) {
        all_ok = (all_ok && busted[i] == busted_ans[i]);        
      }
    }
    busted.resize(0);
    impacts.resize(0);
  }
  check(all_ok, "Clashes are not detected after certain rotations of systems that clearly contain "
        "conflicting atoms.", clashers.getTestingStatus());
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
