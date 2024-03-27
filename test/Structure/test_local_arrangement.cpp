#include "copyright.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Math/series_ops.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Namelists/nml_dynamics.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/parsing_enumerators.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Potential/nonbonded_potential.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/local_arrangement.h"
#include "../../src/Structure/rattle.h"
#include "../../src/Structure/structure_enumerators.h"
#include "../../src/Structure/structure_ops.h"
#include "../../src/Structure/virtual_site_handling.h"
#include "../../src/Synthesis/condensate.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/synthesis_enumerators.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Topology/atomgraph_analysis.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/coordinate_copy.h"
#include "../../src/Trajectory/coordinate_series.h"
#include "../../src/Trajectory/thermostat.h"
#include "../../src/Trajectory/trajectory_enumerators.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"
#ifdef STORMM_USE_HPC
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Trajectory/hpc_integration.h"
#endif

#ifdef STORMM_USE_HPC
using stormm::card::HpcConfig;
#endif
using stormm::card::GpuDetails;
using stormm::constants::small;
using stormm::constants::tiny;
using stormm::data_types::double_type_index;
#ifndef STORMM_USE_HPC
using stormm::data_types::double2;
using stormm::data_types::double3;
using stormm::data_types::double4;
using stormm::data_types::float2;
using stormm::data_types::float3;
using stormm::data_types::float4;
using stormm::data_types::int2;
using stormm::data_types::uint2;
#endif
using stormm::data_types::llint;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getBaseName;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::energy::EvaluateForce;
using stormm::energy::evaluateNonbondedEnergy;
using stormm::energy::ScoreCard;
using stormm::energy::StaticExclusionMask;
using stormm::errors::rtWarn;
using stormm::namelist::DynamicsControls;
using stormm::parse::char4ToString;
using stormm::parse::NumberFormat;
using stormm::random::Xoroshiro128pGenerator;
using stormm::random::Xoshiro256ppGenerator;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using stormm::stmath::computeBoxTransform;
using stormm::symbols::pi;
using stormm::topology::AtomGraph;
using stormm::topology::ConstraintKit;
using stormm::topology::listVirtualSiteFrameTypes;
using stormm::topology::UnitCellType;
using stormm::topology::ValenceKit;
using stormm::topology::VirtualSiteKind;
using namespace stormm::stmath;
using namespace stormm::synthesis;
using namespace stormm::structure;
using namespace stormm::testing;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Scramble the positions of virtual sites in a system.  Scatter frame atoms between box images
// but do not alter their re-imaged relative positions.
//
// Arguments:
//   ps:   Coordinates and box dimensions of the system
//   ag:   System topology (identifies virtual sites)
//   xsr:  Random number generator
//-------------------------------------------------------------------------------------------------
void scrambleSystemCoordinates(PhaseSpace *ps, const AtomGraph &ag, Xoshiro256ppGenerator *xsr) {
  PhaseSpaceWriter psw = ps->data();
  if (psw.unit_cell == UnitCellType::NONE) {
    return;
  }

  // Scramble and recover the virtual site locations
  for (int i = 0; i < ag.getAtomCount(); i++) {
    if (ag.getAtomicNumber(i) == 0) {
      psw.xcrd[i] += xsr->gaussianRandomNumber();
      psw.ycrd[i] += xsr->gaussianRandomNumber();
      psw.zcrd[i] += xsr->gaussianRandomNumber();
    }
    else {
      const double indx = floor(4.0 * (xsr->uniformRandomNumber() - 0.5));
      const double indy = floor(4.0 * (xsr->uniformRandomNumber() - 0.5));
      const double indz = floor(4.0 * (xsr->uniformRandomNumber() - 0.5));
      psw.xcrd[i] += (psw.invu[0] * indx) + (psw.invu[3] * indy) + (psw.invu[6] * indz);
      psw.ycrd[i] += (psw.invu[1] * indx) + (psw.invu[4] * indy) + (psw.invu[7] * indz);
      psw.zcrd[i] += (psw.invu[2] * indx) + (psw.invu[5] * indy) + (psw.invu[8] * indz);
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Center a system on its first atom and re-image it according to the minimum image convention.
//
// Arguments:
//   ps:   Coordinates and box dimensions of the system
//-------------------------------------------------------------------------------------------------
void centerAndReimageSystem(PhaseSpace *ps) {
  PhaseSpaceWriter psw = ps->data();
  const double dx = psw.xcrd[0];
  const double dy = psw.ycrd[0];
  const double dz = psw.zcrd[0];
  for (int i = 0; i < psw.natom; i++) {
    psw.xcrd[i] -= dx;
    psw.ycrd[i] -= dy;
    psw.zcrd[i] -= dz;
  }
  imageCoordinates(ps, ImagingMethod::MINIMUM_IMAGE);
}

//-------------------------------------------------------------------------------------------------
// Check the virtual site placement of a system using code independent of the virtual site
// positioning functions.  Each virtual site gets up to three double-precision results, testing
// its characteristics against expected results.  Tests for the structure will be run internally.
//
// Arguments:
//   ps:        Coordinates of the system
//   ag:        System topology
//   answers:   Packed array of expected results
//   do_tests:  Indicator of whether it is worth pursuing the tests
//-------------------------------------------------------------------------------------------------
void checkVirtualSiteMetrics(const PhaseSpace &ps, const AtomGraph &ag,
                             const std::vector<double3> &answers, const TestPriority do_tests) {
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  const PhaseSpaceReader psr = ps.data();
  std::vector<double3> result(vsk.nsite, { 0.0, 0.0, 0.0 });
  for (int i = 0; i < vsk.nsite; i++) {
    const int vsite_atom = vsk.vs_atoms[i];
    const int parent_atom = vsk.frame1_idx[i];
    const int frame2_atom = vsk.frame2_idx[i];
    const int param_idx = vsk.vs_param_idx[i];
    const std::vector<double> p_vs = { psr.xcrd[vsite_atom] - psr.xcrd[parent_atom],
                                       psr.ycrd[vsite_atom] - psr.ycrd[parent_atom],
                                       psr.zcrd[vsite_atom] - psr.zcrd[parent_atom] };
    const std::vector<double> p_f2 = { psr.xcrd[frame2_atom] - psr.xcrd[parent_atom],
                                       psr.ycrd[frame2_atom] - psr.ycrd[parent_atom],
                                       psr.zcrd[frame2_atom] - psr.zcrd[parent_atom] };
    const double p_vs_distance = sqrt((p_vs[0] * p_vs[0]) + (p_vs[1] * p_vs[1]) +
                                      (p_vs[2] * p_vs[2]));
    const double p_f2_distance = sqrt((p_f2[0] * p_f2[0]) + (p_f2[1] * p_f2[1]) +
                                      (p_f2[2] * p_f2[2]));
    int frame3_atom, frame4_atom;
    double p_f3_distance, p_f4_distance;
    std::vector<double> p_f3(3), p_f4(3);
    const VirtualSiteKind kind = static_cast<VirtualSiteKind>(vsk.vs_types[param_idx]);

    // Get details of the the third frame atom
    switch (kind) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
    case VirtualSiteKind::FIXED_4:
      frame3_atom = vsk.frame3_idx[i];
      p_f3 = { psr.xcrd[frame3_atom] - psr.xcrd[parent_atom],
               psr.ycrd[frame3_atom] - psr.ycrd[parent_atom],
               psr.zcrd[frame3_atom] - psr.zcrd[parent_atom] };
      break;
    case VirtualSiteKind::NONE:
      break;
    }

    // Get details of the the fourth frame atom
    switch (kind) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      break;
    case VirtualSiteKind::FIXED_4:
      frame4_atom = vsk.frame4_idx[i];
      p_f4 = { psr.xcrd[frame4_atom] - psr.xcrd[parent_atom],
               psr.ycrd[frame4_atom] - psr.ycrd[parent_atom],
               psr.zcrd[frame4_atom] - psr.zcrd[parent_atom] };
      break;
    case VirtualSiteKind::NONE:
      break;
    }

    // Evaluate each frame
    switch (kind) {
    case VirtualSiteKind::FLEX_2:
      {
        // Check the ratio of distances between the virtual site and its parent atom
        result[i].x = p_vs_distance / p_f2_distance;

        // Check co-linearity and fill in the second slot
        result[i].y = dot(p_f2, p_vs) / (magnitude(p_f2) * magnitude(p_vs));
      }
      break;
    case VirtualSiteKind::FIXED_2:
      {
        // Check the distance between the virtual site and its parent atom
        result[i].x = p_vs_distance;

        // Check co-linearity and fill in the second slot
        result[i].y = dot(p_f2, p_vs) / (magnitude(p_f2) * magnitude(p_vs));
      }
      break;
    case VirtualSiteKind::FLEX_3:
      {
        // Check the distance between the virtual site and its parent atom, comparing it to the
        // expected distance.
        std::vector<double> overall_displacement(3);
        for (int j = 0; j < 3; j++) {
          overall_displacement[j] = (p_f2[j] * vsk.dim1[i]) + (p_f3[j] * vsk.dim2[i]);
        }
        result[i].x = magnitude(p_vs) / magnitude(overall_displacement);

        // Check the co-planarity of the virtual site with its frame and fill in the second slot.
        // This should be zero.
        result[i].y = pointPlaneDistance(p_f2, p_f3, p_vs);
      }
      break;
    case VirtualSiteKind::FIXED_3:
      {
        // Check the distance between the virtual site and its parent atom.
        result[i].x = magnitude(p_vs);

        // Check the co-planarity of the virtual site with its frame and fill in the second slot.
        // This should be zero.
        result[i].y = pointPlaneDistance(p_f2, p_f3, p_vs);
      }
      break;
    case VirtualSiteKind::FAD_3:
      {
        // Check the distance between the virtual site and its parent atom.
        result[i].x = magnitude(p_vs);

        // Check the co-planarity of the virtual site with its frame and fill in the second slot.
        // This should be zero.
        result[i].y = pointPlaneDistance(p_f2, p_f3, p_vs);

        // Check the angle between the virtual site and the frame atoms.
        const std::vector<double> f2_f3 = { psr.xcrd[frame3_atom] - psr.xcrd[frame2_atom],
                                            psr.ycrd[frame3_atom] - psr.ycrd[frame2_atom],
                                            psr.zcrd[frame3_atom] - psr.zcrd[frame2_atom] };
        const std::vector<double> f2_p = { psr.xcrd[parent_atom] - psr.xcrd[frame2_atom],
                                           psr.ycrd[parent_atom] - psr.ycrd[frame2_atom],
                                           psr.zcrd[parent_atom] - psr.zcrd[frame2_atom] };
        const std::vector<double> f23_t_f2p = perpendicularComponent(f2_f3, f2_p);
        const std::vector<double> pvs_t_f2p = perpendicularComponent(p_vs, f2_p);
        result[i].z = angleBetweenVectors(p_f2, p_vs) + angleBetweenVectors(f23_t_f2p, pvs_t_f2p);
      }
      break;
    case VirtualSiteKind::OUT_3:
      {
        // Check the distance between the virtual site and its parent atom.
        result[i].x = magnitude(p_vs);

        // Check the distance between this point and the plane, comparing it to the distance
        // prescribed by the ratio of the cross product between the vectors defining the plane.
        std::vector<double> f2_x_f3(3);
        crossProduct(p_f2, p_f3, &f2_x_f3);
        result[i].y = pointPlaneDistance(p_f2, p_f3, p_vs) - (magnitude(f2_x_f3) * vsk.dim1[i] *
                                                              vsk.dim2[i] * vsk.dim3[i]);

        // Compute the distance of the virtual site from the axis perpendicular to the plane
        // running through the parent atom.
        std::vector<double> vs_t = perpendicularComponent(p_vs, f2_x_f3);
        result[i].z = magnitude(vs_t);
      }
      break;
    case VirtualSiteKind::FIXED_4:
      {
        // Check the distance between the virtual site and its parent atom.
        result[i].x = magnitude(p_vs);
      }
      break;
    case VirtualSiteKind::NONE:
      break;
    }
  }

  // Compile a list of all frame types found in this system.
  const std::string frame_type_list = listVirtualSiteFrameTypes(ag);
  
  // Compare the first result
  std::vector<double> rvec(vsk.nsite), avec(vsk.nsite);
  for (int i = 0; i < vsk.nsite; i++) {
    rvec[i] = result[i].x;
    avec[i] = answers[i].x;
  }
  check(rvec, RelationalOperator::EQUAL, Approx(avec).margin(small), "Metric 1, measuring "
        "relative and absolute distances between virtual particles and their parent atoms, fails "
        "in a system described by topology " + ag.getFileName() + ".  Virtual site types present "
        "in this system: " + frame_type_list + ".", do_tests);
  for (int i = 0; i < vsk.nsite; i++) {
    rvec[i] = result[i].y;
    avec[i] = answers[i].y;
  }
  check(rvec, RelationalOperator::EQUAL, Approx(avec).margin(small), "Metric 2, measuring dot "
        "products and expected distances to test colinearity or point-plane distances, fails in a "
        "system described by topology " + ag.getFileName() + ".  Virtual site types present in "
        "this system: " + frame_type_list + ".", do_tests);
  for (int i = 0; i < vsk.nsite; i++) {
    rvec[i] = result[i].z;
    avec[i] = answers[i].z;
  }
  check(rvec, RelationalOperator::EQUAL, Approx(avec).margin(1.0e-7), "Metric 3, measuring "
        "miscellaneous quantities like the distance between the virtual site and a normal vector "
        "(Out-3) or the angle made by a virtual site with its nearby frame atoms (FAD-3), fails "
        "in a system described by topology " + ag.getFileName() + ".  Virtual site types present "
        "in this system: " + frame_type_list + ".", do_tests);        
}

//-------------------------------------------------------------------------------------------------
// Check the virtual site force transfer, first with analytic comparisons of the torque on the
// frame atoms and then by finite difference computations.  Electrostatic forces exerted by all
// particles on each other, in the minimum image convention as all of the systems have been placed
// near the origin, will provide forces to transmit from the virtual sites to their frame atoms.
//
// Arguments:
//   ps:        Coordinates and forces for the system
//   ag:        System topology
//   do_tests:  Indicator of whether to attempt tests
//-------------------------------------------------------------------------------------------------
void checkVirtualSiteForceXfer(PhaseSpace *ps, const AtomGraph *ag, const TestPriority do_tests) {

  // The virtual sites should already be in position as of the time this function is called, but
  // place them again rather than resting on that assumption.
  placeVirtualSites(ps, ag);
  const StaticExclusionMask se(ag);
  ScoreCard sc(1);
  double2 base_qlj;
  std::string vs_type_list;
  switch (do_tests) {
  case TestPriority::CRITICAL:
  case TestPriority::NON_CRITICAL:
    base_qlj = evaluateNonbondedEnergy(*ag, se, ps, &sc, EvaluateForce::YES,
                                       EvaluateForce::NO, 0);
    vs_type_list = listVirtualSiteFrameTypes(ag);
    break;
  case TestPriority::ABORT:
    break;
  }
  
  // Compute torque on the molecule, transmit forces to frame atoms, recalculate torque
  const double3 pre_xfer_torque = molecularTorque(ag, ps, 0);
  transmitVirtualSiteForces(ps, ag);
  const double3 pst_xfer_torque = molecularTorque(ag, ps, 0);
  const std::vector<double> t_before = { pre_xfer_torque.x, pre_xfer_torque.y, pre_xfer_torque.z };
  const std::vector<double> t_after  = { pst_xfer_torque.x, pst_xfer_torque.y, pst_xfer_torque.z };
  check(t_before, RelationalOperator::EQUAL, t_after, "Net torque on a molecule bearing virtual "
        "sites of type(s) " + vs_type_list + " is not conserved when transmitting the forces from "
        "virtual sites to frame atoms.", do_tests);

  // Move each of the virtual sites' frame atoms slightly, reposition the sites, recalculate the
  // energy, and compare the analytic forces.
  const std::vector<double> analytic_frc = ps->getInterlacedCoordinates(TrajectoryKind::FORCES);
  const VirtualSiteKit<double> vsk = ag->getDoublePrecisionVirtualSiteKit();
  const int natom = ag->getAtomCount();
  std::vector<bool> is_frame_atom(natom);
  const int site_limit = (ag->getMoleculeCount() > 1 && natom > 200) ? 1 : vsk.nsite;
  for (int i = 0; i < site_limit; i++) {
    is_frame_atom[vsk.frame1_idx[i]] = true;
    is_frame_atom[vsk.frame2_idx[i]] = true;
    switch (static_cast<VirtualSiteKind>(vsk.vs_types[i])) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      is_frame_atom[vsk.frame3_idx[i]] = true;
      break;
    case VirtualSiteKind::FIXED_4:
      is_frame_atom[vsk.frame3_idx[i]] = true;
      is_frame_atom[vsk.frame4_idx[i]] = true;
      break;
    case VirtualSiteKind::NONE:
      break;
    }
  }
  const double fd_perturbation = 0.000001;
  PhaseSpaceWriter psw = ps->data();
  std::vector<double> finite_difference_frc(3 * natom, 0.0);
  std::vector<double> analytic_frc_sample, finite_difference_frc_sample;
  for (int i = 0; i < natom; i++) {
    if (is_frame_atom[i] == false) {
      continue;
    }
    psw.xcrd[i] += fd_perturbation;
    placeVirtualSites(ps, ag);
    const double2 xpos_qlj = evaluateNonbondedEnergy(*ag, se, ps, &sc);
    psw.xcrd[i] -= fd_perturbation;
    psw.ycrd[i] += fd_perturbation;
    placeVirtualSites(ps, ag);
    const double2 ypos_qlj = evaluateNonbondedEnergy(*ag, se, ps, &sc);
    psw.ycrd[i] -= fd_perturbation;
    psw.zcrd[i] += fd_perturbation;
    placeVirtualSites(ps, ag);
    const double2 zpos_qlj = evaluateNonbondedEnergy(*ag, se, ps, &sc);
    psw.zcrd[i] -= fd_perturbation;
    finite_difference_frc[ 3 * i     ] = (base_qlj.x - xpos_qlj.x) / fd_perturbation;
    finite_difference_frc[(3 * i) + 1] = (base_qlj.x - ypos_qlj.x) / fd_perturbation;
    finite_difference_frc[(3 * i) + 2] = (base_qlj.x - zpos_qlj.x) / fd_perturbation;
    analytic_frc_sample.push_back(analytic_frc[ 3 * i     ]);
    analytic_frc_sample.push_back(analytic_frc[(3 * i) + 1]);
    analytic_frc_sample.push_back(analytic_frc[(3 * i) + 2]);
    finite_difference_frc_sample.push_back(analytic_frc[ 3 * i     ]);
    finite_difference_frc_sample.push_back(analytic_frc[(3 * i) + 1]);
    finite_difference_frc_sample.push_back(analytic_frc[(3 * i) + 2]);
  }
  check(analytic_frc_sample, RelationalOperator::EQUAL,
        Approx(finite_difference_frc_sample).margin(1.0e-4), "Analytic forces and those computed "
        "by a finite difference scheme do not agree in a system bearing virtual sites of type(s)" +
        vs_type_list + ".", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Check the local geometry calculations for all types of coordinate objects.
//
// Arguments:
//   tsm:        A collection of test systems, including topologies and coordinates, along with
//               information verifying their availability
//   uc_choice:  Choices for the unit cell types to test as a group
//-------------------------------------------------------------------------------------------------
void checkGeometrics(const TestSystemManager &tsm, const std::vector<UnitCellType> &uc_choice) {
  const std::vector<int> uc_examples = tsm.getQualifyingSystems(uc_choice);
  const int n_examples = uc_examples.size();
  std::vector<CoordinateFrame> tsm_cf;
  std::vector<PhaseSpace> tsm_ps;
  std::vector<CoordinateSeries<llint>> tsm_cs;
  std::vector<int> atom_counts(n_examples);
  tsm_cf.reserve(n_examples);
  tsm_ps.reserve(n_examples);
  tsm_cs.reserve(n_examples);
  for (int i = 0; i < n_examples; i++) {
    tsm_cf.push_back(tsm.exportCoordinateFrame(uc_examples[i]));
    tsm_ps.push_back(tsm.exportPhaseSpace(uc_examples[i]));
    tsm_cs.push_back(tsm.exportCoordinateSeries<llint>(uc_examples[i], 4, 0.0, 918404, 24));
    atom_counts[i] = tsm_cf.back().getAtomCount();
  }
  PhaseSpaceSynthesis tsm_poly_ps = tsm.exportPhaseSpaceSynthesis(uc_examples, 0.0, 97, 40);
  Condensate tsm_cdns(tsm_poly_ps, PrecisionModel::DOUBLE);
  std::vector<double> test_distances(n_examples * 5), test_distances_ans(n_examples * 5);
  std::vector<double> test_angles(n_examples * 5), test_angles_ans(n_examples * 5);
  std::vector<double> test_torsions(n_examples * 5), test_torsions_ans(n_examples * 5);
  for (int i = 0; i < n_examples; i++) {
    const int quarter_atom = atom_counts[i] / 4;
    const int halfway_atom = atom_counts[i] / 2;
    const int thr_qrt_atom = (3 * atom_counts[i]) / 4;
    test_distances[(n_examples * i)    ] = distance(0, halfway_atom, tsm_cf[i]);
    test_distances[(n_examples * i) + 1] = distance(0, halfway_atom, tsm_ps[i]);
    test_distances[(n_examples * i) + 2] = distance<llint, double>(0, halfway_atom,
                                                                   tsm_cs[i], 3);
    test_distances[(n_examples * i) + 3] = distance<double>(0, halfway_atom, tsm_poly_ps, i);
    test_distances[(n_examples * i) + 4] = distance<double>(0, halfway_atom, tsm_cdns, i);
    for (int j = 0; j < 5; j++) {
      test_distances_ans[(n_examples * i) + j] = test_distances[(n_examples * i) + j];
    }
    test_angles[(n_examples * i)    ] = angle(0, quarter_atom, halfway_atom, tsm_cf[i]);
    test_angles[(n_examples * i) + 1] = angle(0, quarter_atom, halfway_atom, tsm_ps[i]);
    test_angles[(n_examples * i) + 2] = angle<llint, double>(0, quarter_atom, halfway_atom,
                                                             tsm_cs[i], 3);
    test_angles[(n_examples * i) + 3] = angle<double>(0, quarter_atom, halfway_atom, tsm_poly_ps,
                                                      i);
    test_angles[(n_examples * i) + 4] = angle<double>(0, quarter_atom, halfway_atom, tsm_cdns, i);
    for (int j = 0; j < 5; j++) {
      test_angles_ans[(n_examples * i) + j] = test_angles[(n_examples * i) + j];
    }
    test_torsions[(n_examples * i)    ] = dihedralAngle(0, quarter_atom, halfway_atom,
                                                        thr_qrt_atom, tsm_cf[i]);
    test_torsions[(n_examples * i) + 1] = dihedralAngle(0, quarter_atom, halfway_atom,
                                                        thr_qrt_atom, tsm_ps[i]);
    test_torsions[(n_examples * i) + 2] = dihedralAngle<llint, double>(0, quarter_atom,
                                                                       halfway_atom, thr_qrt_atom,
                                                                       tsm_cs[i], 3);
    test_torsions[(n_examples * i) + 3] = dihedralAngle<double>(0, quarter_atom, halfway_atom,
                                                                thr_qrt_atom, tsm_poly_ps, i);
    test_torsions[(n_examples * i) + 4] = dihedralAngle<double>(0, quarter_atom, halfway_atom,
                                                                thr_qrt_atom, tsm_cdns, i);
    for (int j = 0; j < 5; j++) {
      test_torsions_ans[(n_examples * i) + j] = test_torsions[(n_examples * i) + j];
    }    
  }
  const std::string sys_desc = (uc_choice[0] == UnitCellType::NONE) ? "isolated" : "periodic";
  check(test_distances, RelationalOperator::EQUAL, Approx(test_distances_ans).margin(1.0e-5),
        "Distance measurements made in a collection of systems with " + sys_desc + " boundary "
        "conditions do not meet expectations.", tsm.getTestingStatus(uc_examples));
  check(test_angles, RelationalOperator::EQUAL, Approx(test_angles_ans).margin(1.0e-5),
        "Angle measurements made in a collection of systems with " + sys_desc + " boundary "
        "conditions do not meet expectations.", tsm.getTestingStatus(uc_examples));
  check(test_torsions, RelationalOperator::EQUAL, Approx(test_torsions_ans).margin(1.0e-5),
        "Dihedral angle measurements made in a collection of systems with " + sys_desc +
        " boundary conditions do not meet expectations.", tsm.getTestingStatus(uc_examples));
}

//-------------------------------------------------------------------------------------------------
// Extract the bond lengths of various atom pairs and fill a new array to use in later analyses.
//
// Arguments:
//   vk:        The valence interactions of the topology of interest
//   pairs:     The list of atom pairs
//   do_tests:  Indicate whether tests are possible
//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<double> listBondLengths(const ValenceKit<T> &vk, const std::vector<int2> &pairs,
                                    const TestPriority do_tests) {
  std::vector<double> result;
  const int nb = pairs.size();
  result.reserve(nb);
  bool all_found = true;
  for (int pos = 0; pos < nb; pos++) {
    const int atom_i = pairs[pos].x;
    const int atom_j = pairs[pos].y;
    bool found = false;
    for (int k = vk.bond_asgn_bounds[atom_i]; k < vk.bond_asgn_bounds[atom_i + 1]; k++) {
      if (vk.bond_asgn_atoms[k] == atom_j) {
        result.push_back(vk.bond_leq[vk.bond_asgn_index[k]]);
        found = true;
      }
    }
    all_found = (all_found && found);
  }
  check(all_found, "Not all bonds stipulated in the constraints were found in the topology.",
        do_tests);
  return result;
}

//-------------------------------------------------------------------------------------------------
// Check the lengths of constraints bonds against the stated values.  This routine will also check
// that the constraints are properly arranged in a singular topology.
//
// Arguments:
//   ps:        Coordinates of the system of interest.  The RATTLE'd coordinates are expected to be
//              in the BLACK cycle position.
//   ag:        Topology of the system of interest
//   prec:      The precision model in which RATTLE was performed, used to select the bond
//              equilibria
//   style:     The manner in which RATTLE groups were evaluated (for error reporting purposes)
//   tol:       The criterion for a successful test of the applied constraints, not necessarily the
//              RATTLE tolerance but close it to at the least
//   do_tests:  Indicate whether tests are possible
//-------------------------------------------------------------------------------------------------
void testBondLengths(const PhaseSpace &ps, const AtomGraph &ag, const PrecisionModel prec,
                     const RattleMethod style, const double tol, const TestPriority do_tests) {

  // Get the target bond lengths
  std::vector<double> constrained_bond_lengths;
  std::vector<int2> constrained_atom_pairs;
  ConstraintKit<double> cnk = ag.getDoublePrecisionConstraintKit();
  int nb = cnk.group_bounds[cnk.ngroup] - cnk.ngroup;
  constrained_bond_lengths.reserve(nb);
  constrained_atom_pairs.reserve(nb);
  for (int i = 0; i < cnk.ngroup; i++) {
    const int atom_i = cnk.group_list[cnk.group_bounds[i]];
    for (int j = cnk.group_bounds[i] + 1; j < cnk.group_bounds[i + 1]; j++) {
      constrained_atom_pairs.push_back({ atom_i, cnk.group_list[j] });
    }
  }
  check(nb, RelationalOperator::EQUAL, constrained_atom_pairs.size(), "The actual number of "
        "constrained bonds is not what would be expected based on the number of groups and the "
        "group sizes.");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    constrained_bond_lengths = listBondLengths<double>(ag.getDoublePrecisionValenceKit(),
                                                       constrained_atom_pairs, do_tests);
    break;
  case PrecisionModel::SINGLE:
    constrained_bond_lengths = listBondLengths<float>(ag.getSinglePrecisionValenceKit(),
                                                      constrained_atom_pairs, do_tests);
    break;
  }

  // Compute the actual distances between the atoms
  const PhaseSpaceReader psr = ps.data(CoordinateCycle::BLACK);
  std::vector<double> rattled_distances(nb);
  for (int pos = 0; pos < nb; pos++) {
    const size_t atom_i = constrained_atom_pairs[pos].x;
    const size_t atom_j = constrained_atom_pairs[pos].y;
    const double dx = psr.xcrd[atom_j] - psr.xcrd[atom_i];
    const double dy = psr.ycrd[atom_j] - psr.ycrd[atom_i];
    const double dz = psr.zcrd[atom_j] - psr.zcrd[atom_i];
    rattled_distances[pos] = sqrt((dx * dx) + (dy * dy) + (dz * dz));
  }
  check(rattled_distances, RelationalOperator::EQUAL, Approx(constrained_bond_lengths).margin(tol),
        "Bond length constraints did not achieve the stipulated tolerance.  Precision model: " +
        getEnumerationName(prec) + ".  Rattle method: " + getEnumerationName(style) + ".",
        do_tests);
}

//-------------------------------------------------------------------------------------------------
// Check the momentum of constrained groups before and after application of RATTLE velocity
// constraints.  Descriptions of input arguments follow from testBondLengths(), above, with a
// significant modification:
//
// Arguments:
//   ps:       Coordinates and velocities of the system of interest.  The RATTLE'd coordinates and
//             velocities expected to be in the BLACK cycle position while the un-RATTLE'd
//             coordinates and velocities are expected to be in the WHITE cycle position.
//-------------------------------------------------------------------------------------------------
void testGroupMomenta(const PhaseSpace &ps, const AtomGraph &ag, const PrecisionModel prec,
                      const RattleMethod style, const double tol, const TestPriority do_tests) {

  ConstraintKit<double> cnk = ag.getDoublePrecisionConstraintKit();
  const PhaseSpaceReader psr = ps.data();
  std::vector<double> masses(psr.natom);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    for (int i = 0; i < psr.natom; i++) {
      masses[i] = ag.getAtomicMass<double>(i);
    }
    break;
  case PrecisionModel::SINGLE:
    for (int i = 0; i < psr.natom; i++) {
      masses[i] = ag.getAtomicMass<float>(i);
    }
    break;
  }

  // Get the original and the RATTLE-adjusted momenta.
  std::vector<double> orig_xmnt(cnk.ngroup, 0.0), rattled_xmnt(cnk.ngroup, 0.0);
  std::vector<double> orig_ymnt(cnk.ngroup, 0.0), rattled_ymnt(cnk.ngroup, 0.0);
  std::vector<double> orig_zmnt(cnk.ngroup, 0.0), rattled_zmnt(cnk.ngroup, 0.0);
  std::vector<double> orig_vel, rattled_vel;
  orig_vel.reserve(3 * cnk.group_bounds[cnk.ngroup]);
  rattled_vel.reserve(3 * cnk.group_bounds[cnk.ngroup]);
  for (int i = 0; i < cnk.ngroup; i++) {
    for (int j = cnk.group_bounds[i]; j < cnk.group_bounds[i + 1]; j++) {
      const size_t atom_ij = cnk.group_list[j];
      const double mass_ij = masses[atom_ij];
      orig_xmnt[i] += psr.xvel[atom_ij] * mass_ij;
      orig_ymnt[i] += psr.yvel[atom_ij] * mass_ij;
      orig_zmnt[i] += psr.zvel[atom_ij] * mass_ij;
      rattled_xmnt[i] += psr.vxalt[atom_ij] * mass_ij;
      rattled_ymnt[i] += psr.vyalt[atom_ij] * mass_ij;
      rattled_zmnt[i] += psr.vzalt[atom_ij] * mass_ij;
      orig_vel.push_back(psr.xvel[atom_ij]);
      orig_vel.push_back(psr.yvel[atom_ij]);
      orig_vel.push_back(psr.zvel[atom_ij]);
      rattled_vel.push_back(psr.vxalt[atom_ij]);
      rattled_vel.push_back(psr.vyalt[atom_ij]);
      rattled_vel.push_back(psr.vzalt[atom_ij]);
    }
  }
  check(orig_xmnt, RelationalOperator::EQUAL, rattled_xmnt, "Application of RATTLE velocity "
        "constraints (" + getEnumerationName(prec) + " precision, " + getEnumerationName(style) +
        " method) did not conserve momentum of the connected groups along the X axis.", do_tests);
  check(orig_ymnt, RelationalOperator::EQUAL, rattled_ymnt, "Application of RATTLE velocity "
        "constraints (" + getEnumerationName(prec) + " precision, " + getEnumerationName(style) +
        " method) did not conserve momentum of the connected groups along the Y axis.", do_tests);
  check(orig_zmnt, RelationalOperator::EQUAL, rattled_zmnt, "Application of RATTLE velocity "
        "constraints (" + getEnumerationName(prec) + " precision, " + getEnumerationName(style) +
        " method) did not conserve momentum of the connected groups along the Z axis.", do_tests);

  // Check the test's validity: the RATTLE'd velocities should differ from the original velocities.
  check(rmsError(orig_vel, rattled_vel), RelationalOperator::GREATER_THAN, 0.01, "The original "
        "and final (rattled) velocities did not register enough difference to constitute a valid "
        "test of the method.");
}

//-------------------------------------------------------------------------------------------------
// Check the operations of RATTLE constraints.
//
// Arguments:
//   tsm:     The collection of test systems
//   prec:    The precision model in which RATTLE was performed, used to select the bond
//            equilibria
//   style:   The manner in which RATTLE groups were evaluated (for error reporting purposes)
//   xsr:     Source of random numbers for adding noise to coordinates and velocities  
//-------------------------------------------------------------------------------------------------
void checkRattleFunctionality(const TestSystemManager &tsm, const PrecisionModel prec,
                              const RattleMethod style, Xoshiro256ppGenerator *xsr) {

  // Count the number of systems with RATTLE hub-and-spoke groups
  const int nsys = tsm.getSystemCount();
  int ntop_with_cnst = 0;
  int ntop_with_sett = 0;
  std::vector<bool> has_cnst(nsys, false);
  std::vector<bool> has_sett(nsys, false);
  for (int i = 0; i < nsys; i++) {
    const AtomGraph& ag = tsm.getTopologyReference(i);
    if (ag.getConstraintGroupCount() > 0) {
      has_cnst[i] = true;
      ntop_with_cnst++;
    }
    if (ag.getSettleGroupCount() > 0) {
      has_sett[i] = true;
      ntop_with_sett++;
    }
  }

  // Extract systems one by one, perturb the relevant atoms, and apply RATTLE.
  for (int i = 0; i < nsys; i++) {
    if (has_cnst[i] == false) {
      continue;
    }

    // Create the coordinates and add random velocities
    PhaseSpace ps = tsm.exportPhaseSpace(i);
    PhaseSpaceWriter psw = ps.data(CoordinateCycle::WHITE);
    PhaseSpaceWriter psw_alt = ps.data(CoordinateCycle::BLACK);
    addRandomNoise(xsr, psw.xvel, psw.yvel, psw.zvel, psw.natom, 1.0);
    CoordinateFrame cf(psw.natom);

    // There is currently no way to copy one coordinate cycle point's contents to the other in
    // a PhaseSpace object.  Use an intermediary object.
    coordCopy(&cf, ps, TrajectoryKind::POSITIONS, CoordinateCycle::WHITE);
    coordCopy(&ps, TrajectoryKind::POSITIONS, CoordinateCycle::BLACK, cf);

    // Perturb the positions of atoms.
    const AtomGraph &ag = tsm.getTopologyReference(i);
    const ConstraintKit cnk = ag.getDoublePrecisionConstraintKit();
    for (int i = 0; i < cnk.ngroup; i++) {
      const int llim = cnk.group_bounds[i];
      const int hlim = cnk.group_bounds[i + 1];
      for (int j = llim; j < hlim; j++) {
        const int jatom = cnk.group_list[j];
        psw_alt.xcrd[jatom] += xsr->gaussianRandomNumber() * 0.05;
        psw_alt.ycrd[jatom] += xsr->gaussianRandomNumber() * 0.05;
        psw_alt.zcrd[jatom] += xsr->gaussianRandomNumber() * 0.05;
      }
    }

    // Restore the correct bond lengths with RATTLE.
    rattlePositions(&ps, ag, prec, 1.0, 1.0e-6, 100, style);
    testBondLengths(ps, ag, prec, style, 2.0e-6, tsm.getTestingStatus(i));
    
    // Further perturb the random velocities, as if performing a force calculation.
    addRandomNoise(xsr, psw_alt.xvel, psw_alt.yvel, psw_alt.zvel, psw.natom, 1.0);
    coordCopy(&cf, ps, TrajectoryKind::VELOCITIES, CoordinateCycle::BLACK);
    coordCopy(&ps, TrajectoryKind::VELOCITIES, CoordinateCycle::WHITE, cf);

    // Restore the correct velocities
    rattleVelocities(&ps, ag, PrecisionModel::DOUBLE, 1.0, 1.0e-6, 100, style);
    testGroupMomenta(ps, ag, prec, style, 2.0e-6, tsm.getTestingStatus(i));
  }
}

//-------------------------------------------------------------------------------------------------
// Encapsulate the constraint instruction tests for a particular precision model.
//
// Arguments:
//   poly_ag:   The synthesis of topologies (included here so that its underlying individual
//              topologies may be queried)
//   syvk:      Valence parameters and valence work units from the topology synthesis
//   poly_auk:  Atom update information from the topology synthesis
//   do_tests:  Indicate whether tests are possible
//-------------------------------------------------------------------------------------------------
template <typename T, typename T2, typename T4>
void synthesisConstraintInsrScan(const AtomGraphSynthesis &poly_ag, const SyValenceKit<T> syvk,
                                 const SyAtomUpdateKit<T, T2, T4> &poly_auk,
                                 const TestPriority do_tests) {

  // Begin by laying out the arrays of all constrained bonds.
  std::vector<int> central_atom_id, peripheral_atom_id;
  std::vector<double> cnst_sum_atom_invmass, cnst_length;
  for (int sysid = 0; sysid < poly_ag.getSystemCount(); sysid++) {
    const AtomGraph* ag_ptr = poly_ag.getSystemTopologyPointer(sysid);
    const ConstraintKit<double> cnk = ag_ptr->getDoublePrecisionConstraintKit();
    const int system_offset = poly_ag.getAtomOffset(sysid);
    for (int i = 0; i < cnk.ngroup; i++) {
      const int atom_llim = cnk.group_bounds[i];
      const int atom_hlim = cnk.group_bounds[i + 1];
      const int param_idx = cnk.group_param_idx[i];
      const int param_llim = cnk.group_param_bounds[param_idx];
      for (int j = atom_llim + 1; j < atom_hlim; j++) { 
        central_atom_id.push_back(cnk.group_list[atom_llim] + system_offset);
        peripheral_atom_id.push_back(cnk.group_list[j] + system_offset);
        const double sq_len = cnk.group_sq_lengths[param_llim + j - atom_llim];
        const double inv_ms = cnk.group_inv_masses[param_llim + j - atom_llim];
        cnst_sum_atom_invmass.push_back(cnk.group_inv_masses[param_llim] + inv_ms);
        cnst_length.push_back(sq_len);
      }
    }
  }
  const int constraint_count = central_atom_id.size();
  std::vector<int> instruction_covered(constraint_count, 0);
  std::vector<int> atom_is_central(constraint_count);
  std::vector<int> atom_is_central_bounds(poly_ag.getPaddedAtomCount() + 1);
  std::vector<int> atom_is_peripheral(constraint_count);
  std::vector<int> atom_is_peripheral_bounds(poly_ag.getPaddedAtomCount() + 1);
  indexingArray(central_atom_id, &atom_is_central, &atom_is_central_bounds);
  indexingArray(peripheral_atom_id, &atom_is_peripheral, &atom_is_peripheral_bounds);

  // Scan all constraint instructions
  std::string param_mistakes;
  int n_param_errors = 0;
  for (int vwuidx = 0; vwuidx < syvk.nvwu; vwuidx++) {

    // Create a cache of the atom imports, much as would be done on the GPU.
    const int2 atom_import_limits = syvk.vwu_abstracts[(vwu_abstract_length * vwuidx) +
                                                       static_cast<int>(VwuAbstractMap::IMPORT)];
    std::vector<int> atom_import_ids(atom_import_limits.y - atom_import_limits.x);
    for (int i = atom_import_limits.x; i < atom_import_limits.y; i++) {
      atom_import_ids[i - atom_import_limits.x] = syvk.vwu_imports[i];
    }

    // Loop over all constraint group instructions, again as would be done on the GPU.
    const int2 cnst_insr_limits = syvk.vwu_abstracts[(vwu_abstract_length * vwuidx) +
                                                     static_cast<int>(VwuAbstractMap::CGROUP)];
    for (int i = cnst_insr_limits.x; i < cnst_insr_limits.y; i++) {

      // Determine which atoms this constraint applies to.
      const uint2 tinsr = poly_auk.cnst_insr[i];
      const int central_atom = atom_import_ids[tinsr.x & 0x3ff];
      const int peripheral_atom = atom_import_ids[(tinsr.x >> 10) & 0x3ff];
      const int param_idx = tinsr.y;

      // Search for the bond constraint among the lists compiled for individual topologies
      const int jlim = atom_is_central_bounds[central_atom + 1];
      for (int j = atom_is_central_bounds[central_atom]; j < jlim; j++) {
        const int jcnst = atom_is_central[j];
        if (central_atom_id[jcnst] != central_atom) {
          rtErr("Atom index " + std::to_string(central_atom) + " is expected to be the central "
                "atom of constraint index " + std::to_string(jcnst) + ", but the central atom "
                "is " + std::to_string(central_atom_id[jcnst]) + " instead.",
                "checkSynthesisConstraintInstructions");
        }
        if (peripheral_atom_id[jcnst] == peripheral_atom) {
          instruction_covered[jcnst] += 1;
          if (fabs(poly_auk.cnst_grp_params[param_idx].x - cnst_length[jcnst]) > 1.0e-6 ||
              fabs(poly_auk.cnst_grp_params[param_idx].y -
                   cnst_sum_atom_invmass[jcnst]) > 1.0e-6) {
            if (n_param_errors < 8) {
              const int sysid = syvk.vwu_abstracts[(vwu_abstract_length * vwuidx) +
                                                   static_cast<int>(VwuAbstractMap::SYSTEM_ID)].x;
              const int top_ca_idx = central_atom - poly_ag.getAtomOffset(sysid);
              const int top_pl_idx = peripheral_atom - poly_ag.getAtomOffset(sysid);
              const AtomGraph* ag_ptr = poly_ag.getSystemTopologyPointer(sysid);
              if (n_param_errors > 0) {
                param_mistakes += ", ";
              }
              param_mistakes += char4ToString(ag_ptr->getAtomName(top_ca_idx)) + " - " +
                                char4ToString(ag_ptr->getAtomName(top_pl_idx)) +
                                "(synthesis indices " + std::to_string(central_atom) + ", " +
                std::to_string(peripheral_atom) + " : expected " +
                realToString(cnst_length[jcnst], 9, 6, NumberFormat::STANDARD_REAL) + " A and " +
                realToString(cnst_sum_atom_invmass[jcnst], 9, 6, NumberFormat::STANDARD_REAL) +
                " mol/g, found " +
                realToString(poly_auk.cnst_grp_params[param_idx].x, 9, 6,
                             NumberFormat::STANDARD_REAL) + " A and " +
                realToString(poly_auk.cnst_grp_params[param_idx].y, 9, 6,
                             NumberFormat::STANDARD_REAL) + " mol/g)";
              n_param_errors++;
            }
          }
        }
      }
    }
  }
  
  // Check the coverage.
  check(instruction_covered, RelationalOperator::GREATER_THAN,
        std::vector<int>(constraint_count, 0), "Constraint instructions of the synthesis did not "
        "cover all constraints found in the individual topologies.", do_tests);

  // Check the parameters.
  check(n_param_errors, RelationalOperator::EQUAL, 0, "Parameters for the constrained bonds were "
        "not rendered as expected.  Examples of erroneous parameter pairs include " +
        param_mistakes + ".", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Check that the AtomGraphSynthesis constructs constraint instructions that fulfill all of the
// constraints in the underlying topologies.
//
// Arguments:
//   tsm:        The collection of test systems, whcih will be parsed for periodic as well as
//               non-periodic systems to make topology syntheses
//   uc_choice:  The type(s) of boundary conditions to seek in choosing systems for the synthesis
//   prec:       The precision model in which to test constraint instructions
//-------------------------------------------------------------------------------------------------
void checkSynthesisConstraintInstructions(const TestSystemManager &tsm,
                                          const std::vector<UnitCellType> &uc_choice,
                                          const PrecisionModel prec) {

  // Create a synthesis of non-periodic systems
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(uc_choice);
  
  // Run through instructions and check off constraints as they are found.  Each thread will get
  // an instruction and then constrain a bond with it.
  switch (prec) {
  case PrecisionModel::DOUBLE:
    synthesisConstraintInsrScan(poly_ag, poly_ag.getDoublePrecisionValenceKit(),
                                poly_ag.getDoublePrecisionAtomUpdateKit(), tsm.getTestingStatus());
      
    break;
  case PrecisionModel::SINGLE:
    synthesisConstraintInsrScan(poly_ag, poly_ag.getSinglePrecisionValenceKit(),
                                poly_ag.getSinglePrecisionAtomUpdateKit(), tsm.getTestingStatus());
    break;
  }
}

//-------------------------------------------------------------------------------------------------
// Check that positional constraints have been correctly applied to a structure within a coordinate
// synthesis.  In question is whether the protocol for fixed-precision coordinates (which, on the
// CPU, is also written to emulate actions that the GPU will take) is able to reproduce the results
// of the CPU double-precision methods.
//
// Arguments:
//   
//-------------------------------------------------------------------------------------------------
template <typename T>
void checkSystemPositionConstraints(const PhaseSpaceWriter &ref_psw,
                                    const PhaseSpaceReader &test_psw, const ConstraintKit<T> &cnk,
                                    const T tol, const std::string &topology_name,
                                    const TestPriority do_tests) {
  std::vector<double> lengths_found;
  std::vector<double> target_lengths;
  for (int i = 0; i < cnk.ngroup; i++) {
    const int llim = cnk.group_bounds[i];
    const int hlim = cnk.group_bounds[i + 1];
    const int cen_atom = cnk.group_list[llim];
    const int param_idx = cnk.group_param_idx[i];
    const int parm_llim = cnk.group_param_bounds[param_idx];
    for (int j = llim + 1; j < hlim; j++) {
      const int dst_atom = cnk.group_list[j];
      const double dx = ref_psw.xcrd[dst_atom] - ref_psw.xcrd[cen_atom];
      const double dy = ref_psw.ycrd[dst_atom] - ref_psw.ycrd[cen_atom];
      const double dz = ref_psw.zcrd[dst_atom] - ref_psw.zcrd[cen_atom];
      const double r  = sqrt((dx * dx) + (dy * dy) + (dz * dz));
      const double r_target = sqrt(cnk.group_sq_lengths[parm_llim + j - llim]);
      lengths_found.push_back(r);
      target_lengths.push_back(r_target);
    }
  }
  check(lengths_found, RelationalOperator::EQUAL, Approx(target_lengths).margin(tol),
        "Bond lengths in a system taken from " + topology_name + " were not constrained to the "
        "expected lengths after positional constraints were applied to the PhaseSpace "
        "representation.", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Check that the results of iterative RATTLE velocity constraints achieved the proper result.
//
// Arguments:
//   poly_ps:   The synthesis of coordinates
//   poly_ag:   The synthesis of topologies
//   prec:      Precision model for the calculations
//   tol:       The tolerance by which RATTLE was applied
//   chk_tol:   The consistency demanded between fixed-precision and floating point constrained
//              coordinate results
//   do_tests:  Indicate whether tests are possible
//-------------------------------------------------------------------------------------------------
template <typename T, typename T2, typename T4>
void checkSynthesisRattleV(const PhaseSpaceSynthesis &poly_ps, const PhaseSpaceSynthesis &raw,
                           const SyValenceKit<T> &syvk, const SyAtomUpdateKit<T, T2, T4> &poly_auk,
                           const T dt, const T tol, const T chk_tol, const HybridTargetLevel tier,
                           const TestPriority do_tests) {

  // Recover the precision model.  This is the price of being able to template code surrounding the
  // topology synthesis.
  const PrecisionModel prec = (std::type_index(typeid(T)).hash_code() == double_type_index) ?
                              PrecisionModel::DOUBLE : PrecisionModel::SINGLE;
  
  // Extract systems from the "raw" coordinates, before constraints were applied.  Use the
  // constraint procedures for single systems as an independent check on the protocols encoded in
  // the topology synthesis and operating on the fixed-precision coordinates.
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    PhaseSpace ref_ps = raw.exportSystem(i, tier);
    const AtomGraph *ag_ptr = raw.getSystemTopologyPointer(i);
    rattleVelocities(&ref_ps, ag_ptr, prec, dt, tol, 30, RattleMethod::CENTER_SUM);
    
    // Compare the alternate approach to processing the raw coordinates to that obtained with the
    // synthesis-based method.  
    PhaseSpaceWriter ref_psw = ref_ps.data();
    const PhaseSpace test_ps = poly_ps.exportSystem(i, tier);

    // Compare the constrained velocities from the double-precision reference calculation to those
    // obtained from the synthesis method, mixing single- or double-precision arithmetic with a
    // fixed-precision coordinate (position and velocity) representation.
    const std::vector<double> ref_xyz =
      ref_ps.getInterlacedCoordinates(CoordinateCycle::BLACK, TrajectoryKind::VELOCITIES);
    const std::vector<double> test_xyz =
      test_ps.getInterlacedCoordinates(CoordinateCycle::BLACK, TrajectoryKind::VELOCITIES);
    check(ref_xyz, RelationalOperator::EQUAL, Approx(test_xyz).margin(chk_tol), "Constrained "
          "velocities from " + getEnumerationName(prec) + "-precision arithmetic operating on "
          "coordinates with " + std::to_string(poly_ps.getGlobalPositionBits()) + " bits of "
          "positional precision and " + std::to_string(poly_ps.getVelocityBits()) + " bits of "
          "velocity precision do not agree with a double-precision standard.  The system topology "
          "is " + getBaseName(ag_ptr->getFileName()) + ".", do_tests);
  }
}

//-------------------------------------------------------------------------------------------------
// Check that the results of iterative RATTLE position constraints achieved the proper result.
// Descriptions of input parameters follow from checkSynthesisRattleV() above.
//-------------------------------------------------------------------------------------------------
template <typename T, typename T2, typename T4>
void checkSynthesisRattleC(const PhaseSpaceSynthesis &poly_ps, const PhaseSpaceSynthesis &raw,
                           const SyValenceKit<T> &syvk, const SyAtomUpdateKit<T, T2, T4> &poly_auk,
                           const T dt, const T tol, const T chk_tol, const HybridTargetLevel tier,
                           const TestPriority do_tests) {

  // Repeat the process above to extract systems from each coordinate synthesis, apply positional
  // RATTLE constraints to the "raw" (unconstrained) coordinates, and check that the results match.
  const PrecisionModel prec = (std::type_index(typeid(T)).hash_code() == double_type_index) ?
                              PrecisionModel::DOUBLE : PrecisionModel::SINGLE;
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    PhaseSpace ref_ps = raw.exportSystem(i, tier);
    const AtomGraph *ag_ptr = raw.getSystemTopologyPointer(i);
    rattlePositions(&ref_ps, ag_ptr, prec, dt, tol, 30, RattleMethod::CENTER_SUM);
    PhaseSpaceWriter ref_psw = ref_ps.data();
    const PhaseSpace test_ps = poly_ps.exportSystem(i, tier);
    
    // Check positions against the CPU results.
    const std::vector<double> ref_xyz =
      ref_ps.getInterlacedCoordinates(CoordinateCycle::BLACK, TrajectoryKind::POSITIONS);
    const std::vector<double> test_xyz =
      test_ps.getInterlacedCoordinates(CoordinateCycle::BLACK, TrajectoryKind::POSITIONS);
    check(ref_xyz, RelationalOperator::EQUAL, Approx(test_xyz).margin(chk_tol), "Constrained "
          "positions from " + getEnumerationName(prec) + "-precision arithmetic operating on "
          "coordinates with " + std::to_string(poly_ps.getGlobalPositionBits()) + " bits of "
          "positional precision and " + std::to_string(poly_ps.getVelocityBits()) + " bits of "
          "velocity precision do not agree with a double-precision standard.  The system topology "
          "is " + getBaseName(ag_ptr->getFileName()) + ".", do_tests);

    // Test velocities to ensure that the correction to account for positional adjustments was
    // successful.
    const std::vector<double> ref_vxyz =
      ref_ps.getInterlacedCoordinates(CoordinateCycle::BLACK, TrajectoryKind::VELOCITIES);
    const std::vector<double> test_vxyz =
      test_ps.getInterlacedCoordinates(CoordinateCycle::BLACK, TrajectoryKind::VELOCITIES);
    check(ref_vxyz, RelationalOperator::EQUAL, Approx(test_vxyz).margin(chk_tol),
          "Constraint-corrected velocities from " + getEnumerationName(prec) + "-precision "
          "arithmetic operating on coordinates with " +
          std::to_string(poly_ps.getGlobalPositionBits()) + " bits of positional precision and " +
          std::to_string(poly_ps.getVelocityBits()) + " bits of velocity precision do not agree "
          "with a double-precision standard.  The system topology is " +
          getBaseName(ag_ptr->getFileName()) + ".", do_tests);
  }
}

//-------------------------------------------------------------------------------------------------
// Implement RATTLE velocity and position constraints for a coordinate synthesis.  
//
// Arguments:
//   tsm:        The collection of test systems, whcih will be parsed for periodic as well as
//               non-periodic systems to make topology syntheses
//   xsr:        Sourcce of random numbers
//   uc_choice:  The type(s) of boundary conditions to seek in choosing systems for the synthesis
//   prec:       The precision model in which to test constraint instructions
//-------------------------------------------------------------------------------------------------
void checkRattle(const TestSystemManager &tsm, Xoshiro256ppGenerator *xsr,
                 const std::vector<UnitCellType> &uc_choice, const PrecisionModel prec,
                 const GpuDetails &gpu = null_gpu) {

  // Create the synthesis
  const size_t nuc = uc_choice.size();
  std::vector<int> index_key;
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    const AtomGraph* iag_ptr = tsm.getTopologyPointer(i);
    for (size_t j = 0; j < nuc; j++) {
      if (iag_ptr->getUnitCellType() == uc_choice[j]) {
        index_key.push_back(i);
      }
    }
  }
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(index_key);
#ifdef STORMM_USE_HPC
  const std::vector<AtomGraph*>& unique_tops = poly_ag.getUniqueTopologies();
  std::vector<StaticExclusionMask> unique_se_masks;
  unique_se_masks.reserve(unique_tops.size());
  for (size_t i = 0; i < unique_tops.size(); i++) {
    unique_se_masks.push_back(StaticExclusionMask(unique_tops[i]));
  }
  StaticExclusionMaskSynthesis poly_se(unique_se_masks, poly_ag.getTopologyIndices());
  poly_ag.loadNonbondedWorkUnits(poly_se, InitializationTask::GENERAL_DYNAMICS);
  const CoreKlManager launcher(gpu, poly_ag);
  poly_ag.upload();
  poly_se.upload();
#endif
  int gbits, vbits;
  double rtol, chk_tol;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    gbits = 60;
    vbits = 66;
    rtol = 1.0e-8;
    chk_tol= 1.2e-6;
    break;
  case PrecisionModel::SINGLE:
    gbits = 28;
    vbits = 40;
    rtol = 1.0e-6;
    chk_tol= 1.5e-6;
    break;
  }
  PhaseSpaceSynthesis poly_ps = tsm.exportPhaseSpaceSynthesis(index_key, 0.0, 5748424, gbits,
                                                              vbits);
  DynamicsControls dyncon;
  dyncon.setThermostatGroup(273.15, 273.15);
  dyncon.setGeometricConstraints();
  std::vector<CoordinateFrame> force_bumps;
  const int nsys = poly_ps.getSystemCount();
  force_bumps.reserve(nsys);
  for (int sysid = 0; sysid < nsys; sysid++) {

    // Copy the system's coordinates.  Set the bond lengths exactly for reference.  Copy those
    // positions back.  Then, scramble the positions and set those as the "developing" positions.
    CoordinateFrame cf = poly_ps.exportCoordinates(sysid);
    CoordinateFrameWriter cfw = cf.data();
    const AtomGraph *ag_ptr = poly_ps.getSystemTopologyPointer(sysid);
    const ConstraintKit<double> cnk = ag_ptr->getDoublePrecisionConstraintKit();
    for (int i = 0; i < cnk.ngroup; i++) {
      const int llim = cnk.group_bounds[i];
      const int hlim = cnk.group_bounds[i + 1];
      const int cen_atom = cnk.group_list[llim];
      const int param_idx = cnk.group_param_idx[i];
      const int parm_llim = cnk.group_param_bounds[param_idx];
      for (int j = llim + 1; j < hlim; j++) {
        const int dst_atom = cnk.group_list[j];
        const double dx = cfw.xcrd[dst_atom] - cfw.xcrd[cen_atom];
        const double dy = cfw.ycrd[dst_atom] - cfw.ycrd[cen_atom];
        const double dz = cfw.zcrd[dst_atom] - cfw.zcrd[cen_atom];
        const double r  = sqrt((dx * dx) + (dy * dy) + (dz * dz));
        const double r_target = sqrt(cnk.group_sq_lengths[parm_llim + j - llim]);
        const double normfac = r_target / r;
        const double norm_dx = dx * normfac;
        const double norm_dy = dy * normfac;
        const double norm_dz = dz * normfac;
        cfw.xcrd[dst_atom] = cfw.xcrd[cen_atom] + norm_dx;
        cfw.ycrd[dst_atom] = cfw.ycrd[cen_atom] + norm_dy;
        cfw.zcrd[dst_atom] = cfw.zcrd[cen_atom] + norm_dz;
      }
    }
    coordCopy(&poly_ps, sysid, TrajectoryKind::POSITIONS, CoordinateCycle::WHITE, cf);

    // Posit a set of velocities for the system.  Use the CoordinateFrame object as a convenient
    // shuttle for the information and necessary type conversions.
    PhaseSpace ps = poly_ps.exportSystem(sysid);
    Thermostat tst(ag_ptr, ThermostatKind::ANDERSEN, 273.15);
    const ThermostatReader<double> tstr = tst.dpData();
    velocityKickStart(&ps, ag_ptr, &tst, dyncon, EnforceExactTemperature::YES);

    // The initialized velocities are in a PhaseSpace object.  Transfer the velocities to the
    // coordinate synthesis.
    coordCopy(&cf, ps, TrajectoryKind::VELOCITIES, CoordinateCycle::WHITE);
    coordCopy(&poly_ps, sysid, TrajectoryKind::VELOCITIES, CoordinateCycle::WHITE, cf);

    // Perturb the velocities, as if by applying a force for half a time step.  Remember the
    // 'forces' (deltas in the velocity of each particle) so that they can applied again later.
    // Store the perturbed velocities in the "developing" BLACK time cycle stage of the
    // coordinate synthesis.
    force_bumps.push_back(cf);
    addRandomNoise(xsr, cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.natom, 0.05);
    CoordinateFrameWriter fb_cfw = force_bumps[sysid].data();
    for (int i = 0; i < cfw.natom; i++) {
      fb_cfw.xcrd[i] = cfw.xcrd[i] - fb_cfw.xcrd[i];
      fb_cfw.ycrd[i] = cfw.ycrd[i] - fb_cfw.ycrd[i];
      fb_cfw.zcrd[i] = cfw.zcrd[i] - fb_cfw.zcrd[i];
    }
    coordCopy(&poly_ps, sysid, TrajectoryKind::VELOCITIES, CoordinateCycle::BLACK, cf);

    // Perturb the coordinates to produce an effect as if they have a history.  This is not a
    // function of the velocities or the velocity perturbation, which will be applied to the
    // coordinates later in the analysis.  Store the perturbed coordinates in the BLACK
    // time cycle of the coordinate synthesis.  Until the actual move, this stage is the history.
    coordCopy(&cf, ps, TrajectoryKind::POSITIONS, CoordinateCycle::WHITE);
    addRandomNoise(xsr, cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.natom, 0.05);
    coordCopy(&poly_ps, sysid, TrajectoryKind::POSITIONS, CoordinateCycle::BLACK, cf);
  }

  // Upload coordinates in the synthesis to the GPU.  This is, in effect, making a copy while the
  // CPU processes coordinates on the host memory.  However, the GPU coordinates will be processed
  // with specific kernels to check for agreement with CPU results.
#ifdef STORMM_USE_HPC
  poly_ps.upload();
  const int2 intg_lp = launcher.getIntegrationKernelDims(prec, AccumulationMethod::SPLIT,
                                                         IntegrationStage::VELOCITY_CONSTRAINT);
  Thermostat velcns_tst(poly_ag, ThermostatKind::NONE, 298.15, PrecisionModel::SINGLE, 41939158,
                        gpu);
  CacheResource intg_cache(intg_lp.x, poly_ag.getValenceWorkUnitSize());
  MolecularMechanicsControls mmctrl_ie(dyncon);
  mmctrl_ie.primeWorkUnitCounters(launcher, EvaluateForce::NO, EvaluateEnergy::YES,
                                  VwuGoal::MOVE_PARTICLES, prec, poly_ag);
  launchIntegrationProcess(&poly_ps, &intg_cache, &velcns_tst, &mmctrl_ie, poly_ag, launcher, prec,
                           AccumulationMethod::SPLIT, IntegrationStage::VELOCITY_CONSTRAINT);
  const std::vector<HybridTargetLevel> tiers = { HybridTargetLevel::HOST,
                                                 HybridTargetLevel::DEVICE };
#else
  const std::vector<HybridTargetLevel> tiers = { HybridTargetLevel::HOST };
#endif
  
  // Make a copy of the coordinate synthesis as it has been constructed: reference positions (in
  // the WHITE stage of the time cycle, which it holds as current) set such that the bond lengths
  // are all correct, random reference velocities are based on some state at 273 Kelvin, and there
  // are some sort of evolving velocities.  If running in CPU mode, check the work of the CPU
  // emulator only.  If running in GPU mode, check both the CPU emulator and the real GPU kernel.
  PhaseSpaceSynthesis raw = poly_ps;
  rattleVelocities(&poly_ps, poly_ag, prec, 1.0, rtol);
  for (size_t i = 0; i < tiers.size(); i++) {
    switch (prec) {
    case PrecisionModel::DOUBLE:
      checkSynthesisRattleV<double,
                            double2, double4>(poly_ps, raw, poly_ag.getDoublePrecisionValenceKit(),
                                              poly_ag.getDoublePrecisionAtomUpdateKit(), 1.0, rtol,
                                              chk_tol, tiers[i], tsm.getTestingStatus(index_key));
      break;
    case PrecisionModel::SINGLE:
      checkSynthesisRattleV<float,
                            float2, float4>(poly_ps, raw, poly_ag.getSinglePrecisionValenceKit(),
                                            poly_ag.getSinglePrecisionAtomUpdateKit(), 1.0, rtol,
                                            chk_tol, tiers[i], tsm.getTestingStatus(index_key));
      break;
    }
  }
  
  // Check that the coordinate cycle for the synthesis remains "WHITE."
  check(getEnumerationName(poly_ps.getCyclePosition()), RelationalOperator::EQUAL,
        getEnumerationName(CoordinateCycle::WHITE), "The coordinate cycle of a synthesis was "
        "not as expected after the velocity constraints protocol.  This may indicate a change "
        "earlier in the workflow.");
  
  // Further advance the velocities according to the prior-computed force bumps, then advance
  // positions.  This will replace the positions in the BLACK time point such that they are
  // now "under construction."
  for (int sysid = 0; sysid < nsys; sysid++) {
    CoordinateFrame vel = poly_ps.exportCoordinates(sysid, CoordinateCycle::BLACK,
                                                    TrajectoryKind::VELOCITIES);
    CoordinateFrameWriter velw = vel.data();
    CoordinateFrameWriter fbw = force_bumps[sysid].data();
    for (int i = 0; i < velw.natom; i++) {
      velw.xcrd[i] += fbw.xcrd[i];
      velw.ycrd[i] += fbw.ycrd[i];
      velw.zcrd[i] += fbw.zcrd[i];
    }
    coordCopy(&poly_ps, sysid, TrajectoryKind::VELOCITIES, CoordinateCycle::BLACK, vel);
    CoordinateFrame pos = poly_ps.exportCoordinates(sysid, CoordinateCycle::WHITE,
                                                    TrajectoryKind::POSITIONS);
    CoordinateFrameWriter posw = pos.data();

    // Advance assuming a 1fs time step.
    for (int i = 0; i < posw.natom; i++) {
      posw.xcrd[i] += velw.xcrd[i];
      posw.ycrd[i] += velw.ycrd[i];
      posw.zcrd[i] += velw.zcrd[i];
    }
    coordCopy(&poly_ps, sysid, TrajectoryKind::POSITIONS, CoordinateCycle::BLACK, pos);
  }
#ifdef STORMM_USE_HPC
  poly_ps.upload();
#endif

  // Update the copy (unconstrained reference) of the coordinate synthesis.
  Hybrid<int2> pairs(nsys);
  for (int i = 0; i < nsys; i++) {
    pairs.putHost({ i, i }, i);
  }
  coordCopy(&raw, poly_ps, pairs);
#ifdef STORMM_USE_HPC
  coordCopy(&raw, poly_ps, pairs, HybridTargetLevel::DEVICE, HybridTargetLevel::DEVICE);
#endif
  // Constrain the positions and adjust the velocities as appropriate.
  rattlePositions(&poly_ps, poly_ag, prec, 1.0, rtol);
#ifdef STORMM_USE_HPC
  launchIntegrationProcess(&poly_ps, &intg_cache, &velcns_tst, &mmctrl_ie, poly_ag, launcher, prec,
                           AccumulationMethod::SPLIT, IntegrationStage::GEOMETRY_CONSTRAINT);
#endif
  for (size_t i = 0; i < tiers.size(); i++) {
    switch (prec) {
    case PrecisionModel::DOUBLE:
      checkSynthesisRattleC<double,
                            double2, double4>(poly_ps, raw, poly_ag.getDoublePrecisionValenceKit(),
                                              poly_ag.getDoublePrecisionAtomUpdateKit(), 1.0, rtol,
                                              chk_tol, tiers[i], tsm.getTestingStatus(index_key));
      break;
    case PrecisionModel::SINGLE:
      checkSynthesisRattleC<float,
                            float2, float4>(poly_ps, raw, poly_ag.getSinglePrecisionValenceKit(),
                                            poly_ag.getSinglePrecisionAtomUpdateKit(), 1.0, rtol,
                                            chk_tol, tiers[i], tsm.getTestingStatus(index_key));
      break;
    }
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

  // Prep the GPU
#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  const Hybrid<int> array_to_trigger_gpu_mapping(1);
#else
  const GpuDetails gpu = null_gpu;
#endif

  // Section 1
  section("Basic checks on re-imaging calculations");

  // Section 2
  section("Check distance, angle, and dihedral calculations");

  // Section 3
  section("Virtual site placement");

  // Section 4
  section("Virtual site force transmission");

  // Section 5
  section("Constraint application");
  
  // Get a realistic system
  const char osc = osSeparator();
  const std::string base_top_path = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_path = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string drug_top_path = base_top_path + osc + "drug_example.top";
  const std::string drug_crd_path = base_crd_path + osc + "drug_example.inpcrd";
  const bool files_exist = (getDrivePathType(drug_top_path) == DrivePathType::FILE &&
                            getDrivePathType(drug_crd_path) == DrivePathType::FILE);
  AtomGraph drug_ag = (files_exist) ? AtomGraph(drug_top_path) : AtomGraph();
  CoordinateFrame drug_cf = (files_exist) ? CoordinateFrame(drug_crd_path,
                                                            CoordinateFileKind::AMBER_INPCRD) :
                                            CoordinateFrame();
  if (files_exist == false) {
    rtWarn("Files for a drug molecule, in water and inside a periodic box, were not found.  Check "
           "the $STORMM_SOURCE environment variable to ensure that " + drug_top_path + " and " +
           drug_crd_path + " become valid paths.  Some tests will be skipped",
           "test_local_arrangement");
  }
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;

  // Image a single number (needed for NMR restraints with periodicity, for example)
  section(1);
  const std::vector<double> point_samples = {  0.75,  0.25,  0.35, -0.65, -1.87,  2.33 };
  const double trial_range = 0.5;
  const std::vector<double> primu_answer  = {  0.25,  0.25,  0.35,  0.35,  0.13,  0.33 };
  const std::vector<double> minim_answer  = { -0.25, -0.25, -0.15, -0.15,  0.13, -0.17 };
  std::vector<double> primu_result(point_samples.size()), minim_result(point_samples.size());
  for (size_t i = 0; i < point_samples.size(); i++) {
    primu_result[i] = imageValue(point_samples[i], trial_range, ImagingMethod::PRIMARY_UNIT_CELL);
    minim_result[i] = imageValue(point_samples[i], trial_range, ImagingMethod::MINIMUM_IMAGE);
  }
  check(primu_result, RelationalOperator::EQUAL, Approx(primu_answer).margin(tiny), "Value "
        "imaging into a [0, range) interval does not work as expected.");
  check(minim_result, RelationalOperator::EQUAL, Approx(minim_answer).margin(tiny), "Value "
        "imaging into a [-0.5 * range, 0.5 * range) interval does not work as expected.");
  
  // Create a fake rectilinear system
  const std::vector<double> rectilinear_box = { 15.0, 24.0, 18.0, 0.5 * pi, 0.5  * pi, 0.5 * pi };
  const std::vector<double> rectilinear_x_crd = {  5.1,  0.7,  7.9, 27.9 };
  const std::vector<double> rectilinear_y_crd = {  3.8, -3.4, 13.5, -9.4 };
  const std::vector<double> rectilinear_z_crd = { -4.9,  1.6, -3.1, -38.2 };
  std::vector<double> rectilinear_umat(9), rectilinear_invu(9);
  computeBoxTransform(rectilinear_box, &rectilinear_umat, &rectilinear_invu);
  double x1 = rectilinear_x_crd[0];
  double y1 = rectilinear_y_crd[0];
  double z1 = rectilinear_z_crd[0];
  imageCoordinates<double, double>(&x1, &y1, &z1, rectilinear_umat.data(), rectilinear_invu.data(),
                                   UnitCellType::ORTHORHOMBIC, ImagingMethod::MINIMUM_IMAGE);
  check(std::vector<double>({ x1, y1, z1}), RelationalOperator::EQUAL,
        std::vector<double>({ 5.1, 3.8, -4.9 }), "Re-imaging a single point in a rectilinear box "
        "by the minimum image convention fails.");
  imageCoordinates<double, double>(&x1, &y1, &z1, rectilinear_umat.data(), rectilinear_invu.data(),
                                   UnitCellType::ORTHORHOMBIC, ImagingMethod::PRIMARY_UNIT_CELL);
  check(std::vector<double>({ x1, y1, z1}), RelationalOperator::EQUAL,
        std::vector<double>({ 5.1, 3.8, 13.1 }), "Re-imaging a single point in a rectilinear box "
        "into the primary unit cell (first octant) fails.");
  std::vector<double> rect_x_copy(rectilinear_x_crd);
  std::vector<double> rect_y_copy(rectilinear_y_crd);
  std::vector<double> rect_z_copy(rectilinear_z_crd);
  imageCoordinates<double, double>(&rect_x_copy, &rect_y_copy, &rect_z_copy,
                                   rectilinear_umat.data(), rectilinear_invu.data(),
                                   UnitCellType::ORTHORHOMBIC, ImagingMethod::MINIMUM_IMAGE);
  check(rect_x_copy, RelationalOperator::EQUAL, std::vector<double>({ 5.1, 0.7, -7.1, -2.1 }),
        "Re-imaging rectilinear Cartesian X coordinates by the minimum image convention fails.");  
  check(rect_y_copy, RelationalOperator::EQUAL, std::vector<double>({ 3.8, -3.4, -10.5, -9.4 }),
        "Re-imaging rectilinear Cartesian Y coordinates by the minimum image convention fails.");  
  check(rect_z_copy, RelationalOperator::EQUAL, std::vector<double>({ -4.9,  1.6, -3.1, -2.2 }),
        "Re-imaging rectilinear Cartesian Z coordinates by the minimum image convention fails.");

  // Create a dense spread of points that will need lots of re-imaging in the rectilinear box
  const int npts = 500;
  std::vector<double> dense_x_crd(npts), dense_y_crd(npts), dense_z_crd(npts);
  Xoshiro256ppGenerator xsr(918733245);
  for (int i = 0; i < npts; i++) {
    dense_x_crd[i] = 1000.0 * (xsr.uniformRandomNumber() - 0.5);
    dense_y_crd[i] = 1000.0 * (xsr.uniformRandomNumber() - 0.5);
    dense_z_crd[i] = 1000.0 * (xsr.uniformRandomNumber() - 0.5);
  }
  std::vector<double> dense_x_copy(dense_x_crd);
  std::vector<double> dense_y_copy(dense_y_crd);
  std::vector<double> dense_z_copy(dense_z_crd);
  imageCoordinates<double, double>(&dense_x_copy, &dense_y_copy, &dense_z_copy,
                                   rectilinear_umat.data(), rectilinear_invu.data(),
                                   UnitCellType::ORTHORHOMBIC, ImagingMethod::MINIMUM_IMAGE);
  std::vector<double> box_x_disp(npts), box_y_disp(npts), box_z_disp(npts);
  for (int i = 0; i < npts; i++) {
    box_x_disp[i] = (dense_x_copy[i] - dense_x_crd[i]) * rectilinear_umat[0];
    box_y_disp[i] = (dense_y_copy[i] - dense_y_crd[i]) * rectilinear_umat[4];
    box_z_disp[i] = (dense_z_copy[i] - dense_z_crd[i]) * rectilinear_umat[8];
    box_x_disp[i] -= round(box_x_disp[i]);
    box_y_disp[i] -= round(box_y_disp[i]);
    box_z_disp[i] -= round(box_z_disp[i]);
  }
  check(box_x_disp, RelationalOperator::EQUAL, Approx(std::vector<double>(npts, 0.0)).margin(tiny),
        "Re-imaging by the minimum image convention moved particles by a non-integral number of "
        "box lengths in the X dimension.");
  check(box_y_disp, RelationalOperator::EQUAL, Approx(std::vector<double>(npts, 0.0)).margin(tiny),
        "Re-imaging by the minimum image convention moved particles by a non-integral number of "
        "box lengths in the Y dimension.");
  check(box_z_disp, RelationalOperator::EQUAL, Approx(std::vector<double>(npts, 0.0)).margin(tiny),
        "Re-imaging by the minimum image convention moved particles by a non-integral number of "
        "box lengths in the Z dimension.");
  std::vector<double> min_rect_crd = { minValue(dense_x_copy), minValue(dense_y_copy),
                                       minValue(dense_z_copy) };
  std::vector<double> max_rect_crd = { maxValue(dense_x_copy), maxValue(dense_y_copy),
                                       maxValue(dense_z_copy) };
  std::vector<double> min_rect_ans = { -0.5 * rectilinear_box[0], -0.5 * rectilinear_box[1],
                                       -0.5 * rectilinear_box[2] };
  std::vector<double> max_rect_ans = { 0.5 * rectilinear_box[0], 0.5 * rectilinear_box[1],
                                       0.5 * rectilinear_box[2] };
  check(min_rect_crd, RelationalOperator::GE, Approx(min_rect_ans).margin(tiny), "Re-imaging by "
        "the minimum image convention puts some particles outside the minimum expected range.");
  check(max_rect_crd, RelationalOperator::LT, Approx(max_rect_ans).margin(tiny), "Re-imaging by "
        "the minimum image convention puts some particles outside the maximum expected range.");
  imageCoordinates<double, double>(&dense_x_copy, &dense_y_copy, &dense_z_copy,
                                   rectilinear_umat.data(), rectilinear_invu.data(),
                                   UnitCellType::ORTHORHOMBIC, ImagingMethod::PRIMARY_UNIT_CELL);
  for (int i = 0; i < npts; i++) {
    box_x_disp[i] = (dense_x_copy[i] - dense_x_crd[i]) * rectilinear_umat[0];
    box_y_disp[i] = (dense_y_copy[i] - dense_y_crd[i]) * rectilinear_umat[4];
    box_z_disp[i] = (dense_z_copy[i] - dense_z_crd[i]) * rectilinear_umat[8];
    box_x_disp[i] -= round(box_x_disp[i]);
    box_y_disp[i] -= round(box_y_disp[i]);
    box_z_disp[i] -= round(box_z_disp[i]);
  }
  check(box_x_disp, RelationalOperator::EQUAL, Approx(std::vector<double>(npts, 0.0)).margin(tiny),
        "Re-imaging into the primary unit cell moved particles by a non-integral number of box "
        "lengths in the X dimension.");
  check(box_y_disp, RelationalOperator::EQUAL, Approx(std::vector<double>(npts, 0.0)).margin(tiny),
        "Re-imaging into the primary unit cell moved particles by a non-integral number of box "
        "lengths in the Y dimension.");
  check(box_z_disp, RelationalOperator::EQUAL, Approx(std::vector<double>(npts, 0.0)).margin(tiny),
        "Re-imaging into the primary unit cell moved particles by a non-integral number of box "
        "lengths in the Z dimension.");
  min_rect_crd = { minValue(dense_x_copy), minValue(dense_y_copy), minValue(dense_z_copy) };
  max_rect_crd = { maxValue(dense_x_copy), maxValue(dense_y_copy), maxValue(dense_z_copy) };
  min_rect_ans = { 0.0, 0.0, 0.0 };
  max_rect_ans = { rectilinear_box[0], rectilinear_box[1], rectilinear_box[2] };
  check(min_rect_crd, RelationalOperator::GE, Approx(min_rect_ans).margin(tiny), "Re-imaging into "
        "the primary unit cell puts some particles outside the minimum expected range.");
  check(max_rect_crd, RelationalOperator::LT, Approx(max_rect_ans).margin(tiny), "Re-imaging into "
        "the primary unit cell puts some particles outside the maximum expected range.");
  
  // Create and test a fake triclinic system
  const std::vector<double> triclinic_box = { 16.0, 17.0, 15.0, 0.6 * pi, 0.53 * pi, 0.55 * pi };
  std::vector<double> triclinic_umat(9), triclinic_invu(9);
  computeBoxTransform(triclinic_box, &triclinic_umat, &triclinic_invu);
  for (int i = 0; i < npts; i++) {
    dense_x_copy[i] = dense_x_crd[i];
    dense_y_copy[i] = dense_y_crd[i];
    dense_z_copy[i] = dense_z_crd[i];
  }
  imageCoordinates<double, double>(&dense_x_copy, &dense_y_copy, &dense_z_copy,
                                   triclinic_umat.data(), triclinic_invu.data(),
                                   UnitCellType::TRICLINIC, ImagingMethod::MINIMUM_IMAGE);
  for (int i = 0; i < npts; i++) {
    const double dx = dense_x_copy[i] - dense_x_crd[i];
    const double dy = dense_y_copy[i] - dense_y_crd[i];
    const double dz = dense_z_copy[i] - dense_z_crd[i];
    box_x_disp[i] = (triclinic_umat[0] * dx) + (triclinic_umat[3] * dy) + (triclinic_umat[6] * dz);
    box_y_disp[i] = (triclinic_umat[1] * dx) + (triclinic_umat[4] * dy) + (triclinic_umat[7] * dz);
    box_z_disp[i] = (triclinic_umat[2] * dx) + (triclinic_umat[5] * dy) + (triclinic_umat[8] * dz);
    box_x_disp[i] -= round(box_x_disp[i]);
    box_y_disp[i] -= round(box_y_disp[i]);
    box_z_disp[i] -= round(box_z_disp[i]);
  }
  std::vector<double> non_integral_disp = { maxAbsValue(box_x_disp), maxAbsValue(box_y_disp),
                                            maxAbsValue(box_z_disp) };
  check(non_integral_disp, RelationalOperator::EQUAL,
        Approx(std::vector<double>(3, 0.0)).margin(tiny), "Re-imaging in a triclinic unit cell by "
        "the minimum image convention moves particles by non-integral box lengths.");
  std::vector<double> frac_tric_x(npts), frac_tric_y(npts), frac_tric_z(npts);
  std::vector<double> tmp_crd(3);
  std::vector<double> tmp_frac(3, 0.0);
  for (int i = 0; i < npts; i++) {
    tmp_crd[0] = dense_x_copy[i];
    tmp_crd[1] = dense_y_copy[i];
    tmp_crd[2] = dense_z_copy[i];
    matrixVectorMultiply(triclinic_umat.data(), tmp_crd.data(), tmp_frac.data(), 3, 3, 1.0, 1.0,
                         0.0);
    frac_tric_x[i] = tmp_frac[0];
    frac_tric_y[i] = tmp_frac[1];
    frac_tric_z[i] = tmp_frac[2];
  }
  std::vector<double> min_tric_frac = { minValue(frac_tric_x), minValue(frac_tric_y),
                                        minValue(frac_tric_z) };
  std::vector<double> max_tric_frac = { maxValue(frac_tric_x), maxValue(frac_tric_y),
                                        maxValue(frac_tric_z) };
  check(min_tric_frac, RelationalOperator::GE, std::vector<double>({ -0.5, -0.5, -0.5 }),
        "Re-imaging in a triclinic unit cell by the minimum image convention puts some particles "
        "outside of the expected range.");
  check(max_tric_frac, RelationalOperator::LT, std::vector<double>({ 0.5, 0.5, 0.5 }),
        "Re-imaging in a triclinic unit cell by the minimum image convention puts some particles "
        "outside of the expected range.");
  for (int i = 0; i < npts; i++) {
    dense_x_copy[i] = dense_x_crd[i];
    dense_y_copy[i] = dense_y_crd[i];
    dense_z_copy[i] = dense_z_crd[i];
  }
  imageCoordinates<double, double>(&dense_x_copy, &dense_y_copy, &dense_z_copy,
                                   triclinic_umat.data(), triclinic_invu.data(),
                                   UnitCellType::TRICLINIC, ImagingMethod::PRIMARY_UNIT_CELL);
  for (int i = 0; i < npts; i++) {
    const double dx = dense_x_copy[i] - dense_x_crd[i];
    const double dy = dense_y_copy[i] - dense_y_crd[i];
    const double dz = dense_z_copy[i] - dense_z_crd[i];
    box_x_disp[i] = (triclinic_umat[0] * dx) + (triclinic_umat[3] * dy) + (triclinic_umat[6] * dz);
    box_y_disp[i] = (triclinic_umat[1] * dx) + (triclinic_umat[4] * dy) + (triclinic_umat[7] * dz);
    box_z_disp[i] = (triclinic_umat[2] * dx) + (triclinic_umat[5] * dy) + (triclinic_umat[8] * dz);
    box_x_disp[i] -= round(box_x_disp[i]);
    box_y_disp[i] -= round(box_y_disp[i]);
    box_z_disp[i] -= round(box_z_disp[i]);
  }
  non_integral_disp = { maxAbsValue(box_x_disp), maxAbsValue(box_y_disp),
                        maxAbsValue(box_z_disp) };
  check(non_integral_disp, RelationalOperator::EQUAL,
        Approx(std::vector<double>(3, 0.0)).margin(tiny), "Re-imaging to the primary unit cell in "
        "triclinic system moves particles by non-integral box lengths.");
  for (int i = 0; i < npts; i++) {
    tmp_crd[0] = dense_x_copy[i];
    tmp_crd[1] = dense_y_copy[i];
    tmp_crd[2] = dense_z_copy[i];
    matrixVectorMultiply(triclinic_umat.data(), tmp_crd.data(), tmp_frac.data(), 3, 3, 1.0, 1.0,
                         0.0);
    frac_tric_x[i] = tmp_frac[0];
    frac_tric_y[i] = tmp_frac[1];
    frac_tric_z[i] = tmp_frac[2];
  }
  min_tric_frac[0] = minValue(frac_tric_x);
  min_tric_frac[1] = minValue(frac_tric_y);
  min_tric_frac[2] = minValue(frac_tric_z);
  max_tric_frac[0] = maxValue(frac_tric_x);
  max_tric_frac[1] = maxValue(frac_tric_y);
  max_tric_frac[2] = maxValue(frac_tric_z);
  check(min_tric_frac, RelationalOperator::GE, std::vector<double>({ 0.0, 0.0, 0.0 }),
        "Re-imaging in a triclinic unit cell by the minimum image convention puts some particles "
        "outside of the expected range.");
  check(max_tric_frac, RelationalOperator::LT, std::vector<double>({ 1.0, 1.0, 1.0 }),
        "Re-imaging in a triclinic unit cell by the minimum image convention puts some particles "
        "outside of the expected range.");  

  // Check internal coordinate computations
  section(2);
  check(distance(0, 1, drug_cf), RelationalOperator::EQUAL, Approx(1.3656697990).margin(small),
        "The distance between two atoms is not computed correctly.");
  CoordinateFrameWriter cfw = drug_cf.data();
  cfw.xcrd[9] +=       cfw.boxdim[0];
  cfw.ycrd[9] += 2.0 * cfw.boxdim[1];
  cfw.zcrd[9] -= 3.0 * cfw.boxdim[2];
  check(distance(9, 10, drug_cf), RelationalOperator::EQUAL, Approx(1.5113186295).margin(small),
        "The distance between two atoms is not computed correctly after one of them has been "
        "pushed into another unit cell image.");
  check(angle(13, 46, 47, drug_cf), RelationalOperator::EQUAL, Approx(1.9124059543).margin(small),
        "The angle made by three atoms in a drug molecule is not computed correctly.");
  cfw.xcrd[22] -= 4.0 * cfw.boxdim[0];
  cfw.ycrd[22] +=       cfw.boxdim[1];
  cfw.zcrd[22] -= 2.0 * cfw.boxdim[2];
  check(angle(40, 22, 50, drug_cf), RelationalOperator::EQUAL, Approx(1.8791980388).margin(small),
        "The angle made by three atoms in a drug molecule is not computed correctly.");
  check(dihedralAngle(9, 10, 11, 12, drug_cf), RelationalOperator::EQUAL,
        Approx(-3.1069347104).margin(small), "The dihedral made by four atoms in a drug molecule "
        "is not computed correctly.");
  check(dihedralAngle(20, 9, 10, 11, drug_cf), RelationalOperator::EQUAL,
        Approx(-1.6233704738).margin(small), "The dihedral made by four atoms in a drug molecule "
        "is not computed correctly.");
  const std::vector<std::string> collection = { "med_1", "med_3", "symmetry_C1", "tip3p",
                                                "drug_example_dry", "ubiquitin", "stereo_L1_vs" };
  TestSystemManager tsm(base_top_path, "top", collection, base_crd_path, "inpcrd", collection);
  std::vector<CoordinateFrame> tsm_cf;
  tsm_cf.reserve(tsm.getSystemCount());
  std::vector<PhaseSpace> tsm_ps;
  tsm_ps.reserve(tsm.getSystemCount());
  std::vector<CoordinateSeries<llint>> tsm_cs;
  tsm_cs.reserve(tsm.getSystemCount());
  std::vector<int> tsm_atom_counts(tsm.getSystemCount());
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    tsm_cf.push_back(tsm.exportCoordinateFrame(i));
    tsm_ps.push_back(tsm.exportPhaseSpace(i));
    tsm_cs.push_back(tsm.exportCoordinateSeries<llint>(i, 4, 0.0, 918404, 24));
    tsm_atom_counts[i] = tsm_cf.back().getAtomCount();
  }
  const std::vector<int> iso_examples = tsm.getQualifyingSystems(UnitCellType::NONE);
  const std::vector<int> pbc_examples = tsm.getQualifyingSystems({ UnitCellType::ORTHORHOMBIC,
                                                                   UnitCellType::TRICLINIC });
  const int n_iso_sys = tsm.getSystemCount(UnitCellType::NONE);
  const int n_pbc_sys = tsm.getSystemCount({ UnitCellType::ORTHORHOMBIC,
                                             UnitCellType::TRICLINIC });
  check(n_iso_sys + n_pbc_sys, RelationalOperator::EQUAL, tsm.getSystemCount(), "The number of "
        "systems with isolated boundary conditions plus the number with periodic boundary "
        "conditions does not match the total number of systems in the TestSystemManager.",
        tsm.getTestingStatus());
  checkGeometrics(tsm, { UnitCellType::NONE });
  checkGeometrics(tsm, { UnitCellType::ORTHORHOMBIC, UnitCellType::TRICLINIC });
  
  // Check the placement of virtual sites, frame type by frame type.  Scramble and recover the
  // virtual site locations, then place the entire system back in the primary unit cell and
  // check the geometries.
  section(3);
  const std::string brbz_top_path = base_top_path + osc + "bromobenzene_vs.top";
  const std::string brbz_crd_path = base_crd_path + osc + "bromobenzene_vs.inpcrd";
  const std::string stro_top_path = base_top_path + osc + "stereo_L1_vs.top";
  const std::string stro_crd_path = base_crd_path + osc + "stereo_L1_vs.inpcrd";
  const std::string symm_top_path = base_top_path + osc + "symmetry_L1_vs.top";
  const std::string symm_crd_path = base_crd_path + osc + "symmetry_L1_vs.inpcrd";
  const std::string dgvs_top_path = base_top_path + osc + "drug_example_vs.top";
  const std::string dgvs_crd_path = base_crd_path + osc + "drug_example_vs.inpcrd";
  const bool vsfi_exist = (getDrivePathType(brbz_top_path) == DrivePathType::FILE &&
                           getDrivePathType(brbz_crd_path) == DrivePathType::FILE &&
                           getDrivePathType(stro_top_path) == DrivePathType::FILE &&
                           getDrivePathType(stro_crd_path) == DrivePathType::FILE &&
                           getDrivePathType(symm_top_path) == DrivePathType::FILE &&
                           getDrivePathType(symm_crd_path) == DrivePathType::FILE &&
                           getDrivePathType(symm_top_path) == DrivePathType::FILE &&
                           getDrivePathType(symm_crd_path) == DrivePathType::FILE);
  const TestPriority do_vs_tests = (vsfi_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  AtomGraph brbz_ag = (vsfi_exist) ? AtomGraph(brbz_top_path) : AtomGraph();
  PhaseSpace brbz_ps = (vsfi_exist) ? PhaseSpace(brbz_crd_path, CoordinateFileKind::AMBER_INPCRD) :
                                      PhaseSpace();
  const int n_brbz_vs = brbz_ag.getVirtualSiteCount();
  std::vector<int> brbz_frame_type_answer(n_brbz_vs);
  brbz_frame_type_answer[0] = static_cast<int>(VirtualSiteKind::FIXED_2);
  for (int i = 1; i < n_brbz_vs; i++) {
    brbz_frame_type_answer[i] = static_cast<int>(VirtualSiteKind::OUT_3);
  }
  std::vector<int> brbz_detected_frame_types(n_brbz_vs);
  for (int i = 0; i < n_brbz_vs; i++) {
    brbz_detected_frame_types[i] = static_cast<int>(brbz_ag.getVirtualSiteFrameType(i));
  }
  check(brbz_frame_type_answer, RelationalOperator::EQUAL, brbz_detected_frame_types,
        "The bromobenzene system, containing a mixture of FIXED_2 and OUT_3 virtual sites, did "
        "not correctly report its virtual site content.", do_vs_tests);  
  PhaseSpaceWriter brbz_psw = brbz_ps.data();
  scrambleSystemCoordinates(&brbz_ps, brbz_ag, &xsr);
  placeVirtualSites(&brbz_ps, brbz_ag);
  centerAndReimageSystem(&brbz_ps);
  const std::vector<double3> brbz_answers = { {  0.800000000, -1.000000000,  0.000000000 },
                                              {  0.434257916,  1.000000000,  0.000000000 },
                                              {  0.434257916,  1.000000000,  0.000000000 },
                                              {  0.429810021,  1.000000000,  0.000000000 },
                                              {  0.429810021,  1.000000000,  0.000000000 },
                                              {  0.441946790,  1.000000000,  0.000000000 },
                                              {  0.441946790,  1.000000000,  0.000000000 },
                                              {  0.440347769,  1.000000000,  0.000000000 },
                                              {  0.440347769,  1.000000000,  0.000000000 },
                                              {  0.438562255,  1.000000000,  0.000000000 },
                                              {  0.438562255,  1.000000000,  0.000000000 },
                                              {  0.421362406,  1.000000000,  0.000000000 },
                                              {  0.421362406,  1.000000000,  0.000000000 } };
  checkVirtualSiteMetrics(brbz_ps, brbz_ag, brbz_answers, do_vs_tests);
  AtomGraph stro_ag = (vsfi_exist) ? AtomGraph(stro_top_path) : AtomGraph();
  PhaseSpace stro_ps = (vsfi_exist) ? PhaseSpace(stro_crd_path, CoordinateFileKind::AMBER_INPCRD) :
                                      PhaseSpace();
  PhaseSpaceWriter stro_psw = stro_ps.data();
  scrambleSystemCoordinates(&stro_ps, stro_ag, &xsr);
  placeVirtualSites(&stro_ps, stro_ag);
  centerAndReimageSystem(&stro_ps);
  const std::vector<double3> stro_answers = { {  1.000000000,  0.000000000,  0.000000000 },
                                              {  0.450000000,  1.000000000,  0.000000000 },
                                              {  0.300000000, -1.000000000,  0.000000000 },
                                              {  0.250000000,  0.000000000,  0.000000000 } };
  checkVirtualSiteMetrics(stro_ps, stro_ag, stro_answers, do_vs_tests);
  AtomGraph symm_ag = (vsfi_exist) ? AtomGraph(symm_top_path) : AtomGraph();
  PhaseSpace symm_ps = (vsfi_exist) ? PhaseSpace(symm_crd_path, CoordinateFileKind::AMBER_INPCRD) :
                                      PhaseSpace();
  PhaseSpaceWriter symm_psw = symm_ps.data();
  scrambleSystemCoordinates(&symm_ps, symm_ag, &xsr);
  placeVirtualSites(&symm_ps, symm_ag);
  centerAndReimageSystem(&symm_ps);
  const std::vector<double3> symm_answers = { {  0.500000000,  0.000000000,  2.094395102 },
                                              {  0.500000000,  0.000000000,  2.094395102 },
                                              {  0.500000000,  0.000000000,  2.094395102 },
                                              {  0.500000000,  0.000000000,  2.094395102 },
                                              {  0.500000000,  0.000000000,  2.094395102 },
                                              {  0.500000000,  0.000000000,  2.094395102 },
                                              {  0.450000000,  0.000000000,  0.000000000 },
                                              {  0.450000000,  0.000000000,  0.000000000 },
                                              {  0.450000000,  0.000000000,  0.000000000 } };
  checkVirtualSiteMetrics(symm_ps, symm_ag, symm_answers, do_vs_tests);
  AtomGraph dgvs_ag = (vsfi_exist) ? AtomGraph(dgvs_top_path) : AtomGraph();
  PhaseSpace dgvs_ps = (vsfi_exist) ? PhaseSpace(dgvs_crd_path, CoordinateFileKind::AMBER_INPCRD) :
                                      PhaseSpace();
  PhaseSpaceWriter dgvs_psw = dgvs_ps.data();
  scrambleSystemCoordinates(&dgvs_ps, dgvs_ag, &xsr);
  placeVirtualSites(&dgvs_ps, dgvs_ag);
  centerAndReimageSystem(&dgvs_ps);
  const std::vector<double3> dgvs_answers = { {  0.333300000,  0.000000000,  0.00000000 },
                                              {  0.333300000,  0.000000000,  0.00000000 },
                                              {  0.420000000,  0.000000000,  0.00000000 } };
  checkVirtualSiteMetrics(dgvs_ps, dgvs_ag, dgvs_answers, do_vs_tests);
  
  // Check force transmission from virtual sites
  section(4);
  checkVirtualSiteForceXfer(&brbz_ps, &brbz_ag, do_vs_tests);
  checkVirtualSiteForceXfer(&stro_ps, &stro_ag, do_vs_tests);
  checkVirtualSiteForceXfer(&symm_ps, &symm_ag, do_vs_tests);
  checkVirtualSiteForceXfer(&dgvs_ps, &dgvs_ag, do_vs_tests);

  // Check restraint applications through RATTLE and SETTLE
  section(5);
  checkRattleFunctionality(tsm, PrecisionModel::DOUBLE, RattleMethod::SEQUENTIAL, &xsr);
  checkRattleFunctionality(tsm, PrecisionModel::DOUBLE, RattleMethod::CENTER_SUM, &xsr);
  checkRattleFunctionality(tsm, PrecisionModel::SINGLE, RattleMethod::SEQUENTIAL, &xsr);
  checkRattleFunctionality(tsm, PrecisionModel::SINGLE, RattleMethod::CENTER_SUM, &xsr);
  checkSynthesisConstraintInstructions(tsm, { UnitCellType::NONE }, PrecisionModel::DOUBLE);
  checkSynthesisConstraintInstructions(tsm, { UnitCellType::NONE }, PrecisionModel::SINGLE);
  const std::vector<UnitCellType> pbct = { UnitCellType::ORTHORHOMBIC, UnitCellType::TRICLINIC };
  checkSynthesisConstraintInstructions(tsm, pbct, PrecisionModel::DOUBLE);
  checkSynthesisConstraintInstructions(tsm, pbct, PrecisionModel::SINGLE);
  checkRattle(tsm, &xsr, std::vector<UnitCellType>(1, UnitCellType::NONE),
              PrecisionModel::DOUBLE, gpu);
  checkRattle(tsm, &xsr, std::vector<UnitCellType>(1, UnitCellType::NONE),
              PrecisionModel::SINGLE, gpu);

  // Summary evaluation
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
