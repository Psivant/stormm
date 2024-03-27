#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/Chemistry/chemical_features.h"
#include "../../src/Chemistry/chemistry_enumerators.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/statistics.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Potential/layered_potential_metrics.h"
#include "../../src/Potential/layered_potential.h"
#include "../../src/Potential/pme_util.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Potential/soft_core_potentials.h"
#include "../../src/Potential/static_exclusionmask.h"
#include "../../src/Potential/nonbonded_potential.h"
#include "../../src/Potential/valence_potential.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/isomerization.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/coordinate_copy.h"
#include "../../src/Trajectory/coordinate_series.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/trajectory_enumerators.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

#ifndef STORMM_USE_HPC
using stormm::double2;
using stormm::double4;
using stormm::int2;
using stormm::int3;
#endif
using stormm::chemistry::ChemicalFeatures;
using stormm::chemistry::IsomerPlan;
using stormm::chemistry::MapRotatableGroups;
using stormm::constants::ExceptionResponse;
using stormm::constants::tiny;
using stormm::data_types::llint;
using stormm::data_types::llint_type_index;
using stormm::data_types::getStormmScalarTypeName;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::energy::ewaldCoefficient;
using stormm::energy::elecPMEDirectSpace;
using stormm::errors::rtWarn;
using stormm::parse::char4ToString;
using stormm::parse::NumberFormat;
using stormm::parse::polyNumericVector;
using stormm::random::Xoshiro256ppGenerator;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using stormm::stmath::mean;
using stormm::stmath::variance;
using stormm::stmath::VarianceMethod;
using stormm::structure::rotateAboutBond;
using stormm::topology::AtomGraph;
using stormm::topology::NonbondedKit;
using namespace stormm::energy;
using namespace stormm::testing;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Check that the sending atom exclusion masks match the receiving atom masks for each tile of
// a StaticExclusionMask object.
//
// Arguments:
//   se:  Non-bonded all-to-all exclusion mask for the system
//-------------------------------------------------------------------------------------------------
int2 checkReflexiveMarks(const StaticExclusionMask &se) {
  int2 result = { 0, 0 };
  const int natom = se.getAtomCount();
  const int nsptile = (natom + 255) / 256;
  for (int sti = 0; sti < nsptile; sti++) {
    const int ni_tile = (sti < nsptile - 1) ? 16 : (natom - (sti * 256)) / 16;
    for (int stj = 0; stj < nsptile; stj++) {
      const int nj_tile = (stj < nsptile - 1) ? 16 : (natom - (stj * 256)) / 16;
      for (int ti = 0; ti < ni_tile; ti++) {
        for (int tj = 0; tj < nj_tile; tj++) {

          // Get the mask for this tile
          std::vector<uint> cmask = se.getTileExclusions(sti, stj, ti, tj);

          // Loop over the first 16 atoms, assemble the corresponding unsigned ints for the
          // second 16, and compare.
          for (int i = 0; i < 16; i++) {
            uint rec_mask = 0x0;
            for (int j = 0; j < 16; j++) {
              rec_mask |= (((cmask[16 + j] >> i) & 0x1) << i);
            }
            if (cmask[i] != rec_mask) {
              for (int j = 0; j < 16; j++) {
                const bool send_excl = ((cmask[i] >> j) & 0x1);
                const bool recv_excl = ((cmask[16 + j] >> i) & 0x1);
                result.x += (send_excl && (recv_excl == false));
                result.y += ((send_excl == false) && recv_excl);
              }
            }
          }
        }
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Check all exclusions from the original topology against those marked in a StaticExclusionMask
// object.
//
// Arguments:
//   se:  Non-bonded all-to-all exclusion mask for the system
//-------------------------------------------------------------------------------------------------
int3 checkMarkedExclusions(const StaticExclusionMask &se) {
  const int natom = se.getTopologyPointer()->getAtomCount();
  const NonbondedKit<double> nbk = se.getTopologyPointer()->getDoublePrecisionNonbondedKit();
  int3 n_errors = { 0, 0, 0 };
  for (int i = 0; i < natom; i++) {
    for (int j = 0; j <= i; j++) {
      const bool is_excl = se.testExclusion(i, j);
      if (j == i && is_excl == false) {
        n_errors.x += 1;
      }
      bool topol_excl = false;
      for (int k = nbk.nb11_bounds[i]; k < nbk.nb11_bounds[i + 1]; k++) {
        topol_excl = (topol_excl || j == nbk.nb11x[k]);
      }
      for (int k = nbk.nb12_bounds[i]; k < nbk.nb12_bounds[i + 1]; k++) {
        topol_excl = (topol_excl || j == nbk.nb12x[k]);
      }
      for (int k = nbk.nb13_bounds[i]; k < nbk.nb13_bounds[i + 1]; k++) {
        topol_excl = (topol_excl || j == nbk.nb13x[k]);
      }
      for (int k = nbk.nb14_bounds[i]; k < nbk.nb14_bounds[i + 1]; k++) {
        topol_excl = (topol_excl || j == nbk.nb14x[k]);
      }
      n_errors.y += (topol_excl && (is_excl == false));
      n_errors.z += ((topol_excl == false) && is_excl && j != i);
    }
  }
  return n_errors;
}

//-------------------------------------------------------------------------------------------------
// Compute the non-bonded forces on a system using the precision model specified by a particular
// parameter set, coordinate, and force representation.
//
// Arguments:
//   nbk:           Non-bonded parameters in a particular precision
//   semask:        Non-bonded exclusions list
//   ps:            Coordinates of the system in double precision
//   elec_ref_frc:  Reference electrostatic forces on all particles
//   vdw_ref_frc:   Reference van-der Waals forces on all particles
//   tol:           Tolerance to apply to each force component on every particle when comparing to
//                  the reference
//   do_tests:      Indicator of whether to pursue tests
//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tforce, typename Tcalc>
void testNBPrecisionModel(const NonbondedKit<Tcalc>nbk, const StaticExclusionMask &semask,
                          const PhaseSpace &ps, const std::vector<double> &elec_ref_frc,
                          const std::vector<double> &vdw_ref_frc, const double tol,
                          const TestPriority do_tests) {
  const size_t tcoord_ct = std::type_index(typeid(Tcoord)).hash_code();
  const size_t tforce_ct = std::type_index(typeid(Tforce)).hash_code();
  const CoordinateSeries<Tcoord> crd(ps, 1, 26 * (tcoord_ct == llint_type_index));
  CoordinateSeries<Tforce> frc(ps, 1, 23 * (tforce_ct == llint_type_index));
  CoordinateSeriesReader crdr = crd.data();
  CoordinateSeriesWriter frcw = frc.data();
  ScoreCard sc(1, 16, 32);
  const Tforce zero = 0.0;
  for (int i = 0; i < frcw.natom; i++) {
    frcw.xcrd[i] = zero;
    frcw.ycrd[i] = zero;
    frcw.zcrd[i] = zero;
  }
  evaluateNonbondedEnergy<Tcoord, Tforce, Tcalc>(nbk, semask.data(), crdr.xcrd, crdr.ycrd,
                                                 crdr.zcrd, crdr.umat, crdr.invu, crdr.unit_cell,
                                                 frcw.xcrd, frcw.ycrd, frcw.zcrd, &sc,
                                                 EvaluateForce::YES, EvaluateForce::NO, 0,
                                                 crdr.inv_gpos_scale, frcw.gpos_scale);
  const CoordinateFrame elec_frame = frc.exportFrame(0);
  const std::vector<double> elec_result = elec_frame.getInterlacedCoordinates();
  for (int i = 0; i < frcw.natom; i++) {
    frcw.xcrd[i] = zero;
    frcw.ycrd[i] = zero;
    frcw.zcrd[i] = zero;
  }
  evaluateNonbondedEnergy<Tcoord, Tforce, Tcalc>(nbk, semask.data(), crdr.xcrd, crdr.ycrd,
                                                 crdr.zcrd, crdr.umat, crdr.invu, crdr.unit_cell,
                                                 frcw.xcrd, frcw.ycrd, frcw.zcrd, &sc,
                                                 EvaluateForce::NO, EvaluateForce::YES, 0,
                                                 crdr.inv_gpos_scale, frcw.gpos_scale);
  const CoordinateFrame vdw_frame = frc.exportFrame(0);
  const std::vector<double> vdw_result = vdw_frame.getInterlacedCoordinates();
  check(elec_result, RelationalOperator::EQUAL, Approx(elec_ref_frc).margin(tol),
        "Electrostatic forces do not agree with the reference when computed in " +
        getStormmScalarTypeName<Tcalc>() + " with " + getStormmScalarTypeName<Tcoord>() +
        " coordinates and " + getStormmScalarTypeName<Tforce>() + " force accumulation.",
        do_tests);
  check(vdw_result, RelationalOperator::EQUAL, Approx(vdw_ref_frc).margin(tol),
        "van-der Waals forces do not agree with the reference when computed in " +
        getStormmScalarTypeName<Tcalc>() + " with " + getStormmScalarTypeName<Tcoord>() +
        " coordinates and " + getStormmScalarTypeName<Tforce>() + " force accumulation.",
        do_tests);
}

//-------------------------------------------------------------------------------------------------
// Test the softcore Coulomb potential.
//
// Arguments:
//   qiqj:            Product of the two atoms' charges (with any implicit attenuation effects)
//   clash_distance:  Minimum distance between the two particles before the softcore potential
//                    kicks in
//   rinc:            [Optional] The sampling size for the test.  The quadratic soft core Coulomb
//                    function wll be tested over a range from rinc to twice the clash distance.
//   order:           The order of softcore potential to use.  Valid inputs are 1 (linear),
//                    2 (quadratic), or 3 (cubic).
//-------------------------------------------------------------------------------------------------
void testSoftCoreCoulomb(const double qiqj, const double clash_distance,
                         const double rinc = 0.001, const int order = 2) {
  const int npts = ceil(2.0 * clash_distance / rinc);
  std::vector<double> analytic_coulomb_e(npts);
  std::vector<double> analytic_coulomb_f(npts);
  std::vector<double> computed_e(npts);
  std::vector<double> computed_f(npts);
  std::vector<double> finite_difference_f(npts);
  bool e_on_track = true;
  bool f_on_track = true;
  bool e_has_kink = false;
  bool f_has_kink = false;
  for (int i = 0; i < npts; i++) {
    const double r = static_cast<double>(i + 1) * rinc;
    analytic_coulomb_e[i] = qiqj / r;
    analytic_coulomb_f[i] = -qiqj / (r * r * r);
    quadraticCoreElectrostatics<double>(r, clash_distance, qiqj, &computed_e[i], &computed_f[i]);
    if (r >= clash_distance) {
      e_on_track = (e_on_track && fabs(computed_e[i] - analytic_coulomb_e[i]) < tiny);
      f_on_track = (f_on_track && fabs(computed_f[i] - analytic_coulomb_f[i]) < tiny);
    }
    if (i > 1) {
      e_has_kink = (e_has_kink || fabs(computed_e[i] - computed_e[i - 1]) >
                                  10.0 * fabs(computed_e[i - 1] - computed_e[i - 2]));
      f_has_kink = (f_has_kink || fabs(computed_f[i] - computed_f[i - 1]) >
                                  10.0 * fabs(computed_f[i - 1] - computed_f[i - 2]));
    }
    double p_nrg = 0.0;
    double n_nrg = 0.0;
    if (order == 1) {
      linearCoreElectrostatics<double>(r + 5.0e-8, clash_distance, qiqj, &p_nrg, nullptr);
      linearCoreElectrostatics<double>(r - 5.0e-8, clash_distance, qiqj, &n_nrg, nullptr);
    }
    else if (order == 2) {
      quadraticCoreElectrostatics<double>(r + 5.0e-8, clash_distance, qiqj, &p_nrg, nullptr);
      quadraticCoreElectrostatics<double>(r - 5.0e-8, clash_distance, qiqj, &n_nrg, nullptr);
    }
    else if (order == 3) {
      cubicCoreElectrostatics<double>(r + 5.0e-8, clash_distance, qiqj, &p_nrg, nullptr);
      cubicCoreElectrostatics<double>(r - 5.0e-8, clash_distance, qiqj, &n_nrg, nullptr);
    }
    else {
      rtErr("Invalid softcore potential order " + std::to_string(order) + " for electrostatic "
            "modification.", "testSoftCoreElectrostatics");
    }

    // The extra division by r is necessary as the force magnitudes emerging from the
    // quadraticCoreElectrostatics() function are divided by the distance to prepare for computing
    // force components.
    finite_difference_f[i] = (p_nrg - n_nrg) / (1.0e-7) / r;
  }
  check(e_on_track, "Energies computed for the non-softcore range of the Coulomb function do not "
        "meet expectations.");
  check(f_on_track, "Forces computed for the non-softcore range of the Coulomb function do not "
        "meet expectations.");
  check(e_has_kink == false, "The quadratic soft core Coulomb function evaluates with a kink in "
        "its energy.");
  check(f_has_kink == false, "The quadratic soft core Coulomb function evaluates with a kink in "
        "its force.");
  check(computed_f, RelationalOperator::EQUAL, Approx(finite_difference_f).margin(1.0e-5),
        "Soft-core Coulomb forces do not agree with a finite difference approximation.");
}

//-------------------------------------------------------------------------------------------------
// Test the softcore Lennard-Jones potential.
//
// Arguments:
//   lja:          The Lennard-Jones A parameter for the atom pair
//   ljb:          The Lennard-Jones B parameter for the atom pair
//   clash_ratio:  Minimum ratio of the distance between the two particles and their pariwise
//                 LennardJones sigma, below which the softcore potential engages
//   rinc:         [Optional] The sampling size for the test.  Quartic soft core Lennard-Jones
//                 evaluations wll be tested over a range from rinc to twice the clash distance.
//   order:        The order of softcore potential to use.  Valid inputs are 3 (cubic),
//                 4 (quartic).
//-------------------------------------------------------------------------------------------------
void testSoftCoreLennardJones(const double lja, const double ljb, const double clash_ratio,
                              const double rinc = 0.001, const int order = 4) {
  const double sigma = (ljb > 1.0e-6) ? sqrt(cbrt(lja / ljb)) : 0.0;
  const int npts = ceil(2.0 * sigma / rinc);
  std::vector<double> analytic_lj_e(npts);
  std::vector<double> analytic_lj_f(npts);
  std::vector<double> computed_e(npts);
  std::vector<double> computed_f(npts);
  std::vector<double> finite_difference_f(npts);
  bool e_on_track = true;
  bool f_on_track = true;
  bool e_has_kink = false;
  bool f_has_kink = false;
  for (int i = 0; i < npts; i++) {
    const double r = static_cast<double>(i + 1) * rinc;
    const double invr = 1.0 / r;
    const double invr2 = invr * invr;
    const double invr4 = invr2 * invr2;
    analytic_lj_e[i] = ((lja * invr2 * invr4) - ljb) * invr2 * invr4;
    analytic_lj_f[i] = ((6.0 * ljb) - (12.0 * lja * invr2 * invr4)) * invr4 * invr4;
    if (order == 3) {
      cubicCoreLennardJones<double>(r, clash_ratio, lja, ljb, sigma, &computed_e[i],
                                    &computed_f[i]);
    }
    else if (order == 4) {
      quarticCoreLennardJones<double>(r, clash_ratio, lja, ljb, sigma, &computed_e[i],
                                      &computed_f[i]);
    }
    else {
      rtErr("Invalid softcore potential order " + std::to_string(order) + " for Lennard-Jones "
            "modification.", "testSoftCoreLennardJones");
    }
    if (r >= clash_ratio * sigma) {
      e_on_track = (e_on_track && fabs(computed_e[i] - analytic_lj_e[i]) < tiny);
      f_on_track = (f_on_track && fabs(computed_f[i] - analytic_lj_f[i]) < tiny);
    }
    if (i > 1) {
      e_has_kink = (e_has_kink || fabs(computed_e[i] - computed_e[i - 1]) >
                                  20.0 * fabs(computed_e[i - 1] - computed_e[i - 2]));
      f_has_kink = (f_has_kink || fabs(computed_f[i] - computed_f[i - 1]) >
                                  20.0 * fabs(computed_f[i - 1] - computed_f[i - 2]));
    }
    double p_nrg = 0.0;
    double n_nrg = 0.0;
    if (order == 3) {
      cubicCoreLennardJones<double>(r + 5.0e-7, clash_ratio, lja, ljb, sigma, &p_nrg, nullptr);
      cubicCoreLennardJones<double>(r - 5.0e-7, clash_ratio, lja, ljb, sigma, &n_nrg, nullptr);
    }
    else if (order == 4) {
      quarticCoreLennardJones<double>(r + 5.0e-7, clash_ratio, lja, ljb, sigma, &p_nrg, nullptr);
      quarticCoreLennardJones<double>(r - 5.0e-7, clash_ratio, lja, ljb, sigma, &n_nrg, nullptr);
    }
    
    // The extra division by r is necessary as the force magnitudes emerging from the
    // quadraticCoreElectrostatics() function are divided by the distance to prepare for computing
    // force components.
    finite_difference_f[i] = (p_nrg - n_nrg) / (1.0e-6) / r;
  }
  check(e_on_track, "Energies computed for the non-softcore range of the Lennard-Jones function "
        "do not meet expectations.");
  check(f_on_track, "Forces computed for the non-softcore range of the Lennard-Jones function do "
        "not meet expectations.");
  check(e_has_kink == false, "The quartic softcore Lennard-Jones function evaluates with a kink "
        "in its energy.");
  check(f_has_kink == false, "The quartic softcore Lennard-Jones function evaluates with a kink "
        "in its force.");
  check(computed_f, RelationalOperator::EQUAL, Approx(finite_difference_f).margin(3.0e-4),
        "Soft-core Lennard-Jones forces do not agree with a finite difference approximation.");
}

//-------------------------------------------------------------------------------------------------
// Test the cubic-, quartic-, and quintic softcore potentials with a typical Lennard-Jones
// function and a standard protein force field.  This will provide a lot of tests as to whether,
// for a reasonable range of applications, the functions could contain artificial minima.
//
// Arguments:
//   ag:            Topology for a protein system with numerous Lennard-Jones pair interactions
//   sigma_factor:  Factor of the sigma parameter at which the handoff between the standard
//                  potential and the softcore functions occurs
//   slope_zero:    Target first derivative for the potentials as the inter-particle distance goes
//                  to zero
//   zero_ratio3:   Expected average ratio of the Lennard Jones well depth to the maximum of the
//                  cubic softcore function at r = 0
//   zero_ratio4:   Expected average ratio of the Lennard Jones well depth to the maximum of the
//                  quartic softcore function at r = 0
//   zero_ratio5:   Expected average ratio of the Lennard Jones well depth to the maximum of the
//                  quintic softcore function at r = 0
//   do_tests:      Indicate whether the tests can be run
//-------------------------------------------------------------------------------------------------
void testSplineSoftCore(const AtomGraph &ag, const double sigma_factor, const double slope_zero,
                        const double zero_ratio3, const double zero_ratio4,
                        const double zero_ratio5, const TestPriority do_tests) {
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  bool cubic_slope_always_negative = true;
  bool quartic_slope_always_negative = true;
  bool quintic_slope_always_negative = true;
  std::vector<double> all_ij3_zratio, all_ij4_zratio, all_ij5_zratio;
  for (int i = 0; i < nbk.n_lj_types; i++) {
    for (int j = 0; j <= i; j++) {
      const size_t ij_ljidx = (j * nbk.n_lj_types) + i;
      const double ij_sigma = nbk.lj_sigma[ij_ljidx];
      const double ij_well_min = pow(2.0, 1.0 / 6.0) * ij_sigma;
      const double rswitch = ij_sigma * sigma_factor;
      const double lja = nbk.lja_coeff[ij_ljidx];
      const double ljb = nbk.ljb_coeff[ij_ljidx];
      const double ij_epsilon = (ij_sigma > 1.0e-6) ? 0.25 * ljb / pow(ij_sigma, 6.0) : 0.0;
      const double inv_rs = 1.0 / rswitch;
      const double inv_rs2 = inv_rs * inv_rs;
      const double inv_rs3 = inv_rs2 * inv_rs;
      const double inv_rs4 = inv_rs2 * inv_rs2;
      const double inv_rs6 = inv_rs3 * inv_rs3;
      const double inv_rs8 = inv_rs4 * inv_rs4;
      const double f_rswitch = ((lja * inv_rs6) - ljb) * inv_rs6;
      const double df_rswitch  = ((  -12.0 * lja * inv_rs6) + (  6.0 * ljb)) * inv_rs6 * inv_rs;
      const double d2f_rswitch = ((  156.0 * lja * inv_rs6) - ( 42.0 * ljb)) * inv_rs6 * inv_rs2;
      const double d3f_rswitch = ((-2184.0 * lja * inv_rs6) + (336.0 * ljb)) * inv_rs6 * inv_rs3;
      
      // Test the cubic softcore potential
      double4 abcd_coefs;
      cubicSoftCore(&abcd_coefs, rswitch, f_rswitch, df_rswitch, slope_zero);
      const double cbfd_a = 3.0 * abcd_coefs.x;
      const double cbfd_b = 2.0 * abcd_coefs.y;
      const double cbfd_c = abcd_coefs.z;

      // Test the minimum and maximum points for nonzero sigma.  Test whether the function's 
      // first derivative crosses zero using the quadratic formula and test any qualifying points.
      std::vector<double> cubic_dd_roots;
      cubic_dd_roots.push_back(0.0);
      cubic_dd_roots.push_back(ij_sigma);
      double discriminant = (cbfd_b * cbfd_b) - (4.0 * cbfd_a * cbfd_c);
      if (discriminant > 0.0) {
        const double zero_p = (-cbfd_b + sqrt(discriminant)) / (2.0 * cbfd_a);
        const double zero_n = (-cbfd_b - sqrt(discriminant)) / (2.0 * cbfd_a);
        if (zero_p > 0.0 && zero_p < ij_sigma) {
          cubic_dd_roots.push_back(zero_p);
        }
        if (zero_n > 0.0 && zero_n < ij_sigma) {
          cubic_dd_roots.push_back(zero_n);
        }        
      }
      for (size_t i = 0; i < cubic_dd_roots.size(); i++) {
        const double r = cubic_dd_roots[i];
        const double cbfd_val = (((cbfd_a * r) + cbfd_b) * r) + cbfd_c;
        if (cbfd_val > 1.0e-6) {
          cubic_slope_always_negative = false;
        }
      }
      
      // Add the ratio of the cubic function's maximum height (which, if the previous tests have
      // passed, occurs at zero) to the original Lennard-Jones well depth.
      if (ij_sigma > 1.0e-6) {
        all_ij3_zratio.push_back(abcd_coefs.w / ij_epsilon);
      }
      
      // Test the quartic softcore potential
      double e_coef;
      quarticSoftCore(&abcd_coefs, &e_coef, rswitch, f_rswitch, df_rswitch, d2f_rswitch,
                      slope_zero);
      const double qrfd_a = 4.0 * abcd_coefs.x;
      const double qrfd_b = 3.0 * abcd_coefs.y;
      const double qrfd_c = 2.0 * abcd_coefs.z;
      const double qrfd_d = abcd_coefs.w;
      const double qrf2d_a = 12.0 * abcd_coefs.x;
      const double qrf2d_b =  6.0 * abcd_coefs.y;
      const double qrf2d_c =  2.0 * abcd_coefs.z;

      // Test the minimum and maximum points for nonzero sigma.  Find the roots of the second
      // derivative using the quadratic formula and test any qualifying points.
      std::vector<double> quartic_dd_roots;
      quartic_dd_roots.push_back(0.0);
      quartic_dd_roots.push_back(ij_sigma);
      discriminant = (qrf2d_b * qrf2d_b) - (4.0 * qrf2d_a * qrf2d_c);
      if (discriminant > 0.0) {
        const double zero_p = (-qrf2d_b + sqrt(discriminant)) / (2.0 * qrf2d_a);
        const double zero_n = (-qrf2d_b - sqrt(discriminant)) / (2.0 * qrf2d_a);
        if (zero_p > 0.0 && zero_p < ij_sigma) {
          quartic_dd_roots.push_back(zero_p);
        }
        if (zero_n > 0.0 && zero_n < ij_sigma) {
          quartic_dd_roots.push_back(zero_n);
        }
      }
      for (size_t i = 0; i < quartic_dd_roots.size(); i++) {
        const double r = quartic_dd_roots[i];
        const double qrfd_val = (((((qrfd_a * r) + qrfd_b) * r) + qrfd_c) * r) + qrfd_d;
        if (qrfd_val > 1.0e-6) {
          quartic_slope_always_negative = false;
        }
      }

      // Add the ratio of the quartic function's maximum height (which, if the previous tests have
      // passed, occurs at zero) to the original Lennard-Jones well depth.
      if (ij_sigma > 1.0e-6) {
        all_ij4_zratio.push_back(e_coef / ij_epsilon);
      }
      
      // Test the quintic softcore potential
      double2 ef_coefs;
      quinticSoftCore(&abcd_coefs, &ef_coefs, rswitch, f_rswitch, df_rswitch, d2f_rswitch,
                      d3f_rswitch, slope_zero);
      const double qnfd_a = 5.0 * abcd_coefs.x;
      const double qnfd_b = 4.0 * abcd_coefs.y;
      const double qnfd_c = 3.0 * abcd_coefs.z;
      const double qnfd_d = 2.0 * abcd_coefs.w;
      const double qnfd_e = ef_coefs.x;
      const double qnf2d_a = 20.0 * abcd_coefs.x;
      const double qnf2d_b = 12.0 * abcd_coefs.y;
      const double qnf2d_c =  6.0 * abcd_coefs.z;
      const double qnf2d_d =  2.0 * abcd_coefs.w;

      // Iterate over the limits of the softcore function to see if there are any places that
      // the second derivative hits zero, indicating an extremum in the first derivative which
      // would then need to be checked.
      std::vector<double> quintic_dd_roots;
      quintic_dd_roots.push_back(0.0);
      quintic_dd_roots.push_back(ij_sigma);
      const double rdisc = 0.0001;
      const double rlim = std::min(ij_sigma, rswitch - rdisc);
      for (double r = 0.0; r < rlim; r += rdisc) {
        const double rp = r + rdisc;
        const double ddf_low  = (((((qnf2d_a * r) + qnf2d_b) * r) + qnf2d_c) * r) + qnf2d_d;
        const double ddf_high = (((((qnf2d_a * rp) + qnf2d_b) * rp) + qnf2d_c) * rp) + qnf2d_d;
        if (ddf_low == 0.0) {
          quintic_dd_roots.push_back(r);
        }
        else if (ddf_high == 0.0) {
          quintic_dd_roots.push_back(rp);
        }
        else if ((ddf_low < 0.0 && ddf_high > 0.0) || (ddf_low > 0.0 && ddf_high < 0.0)) {

          // Solve the linear approximation for the root
          const double rise = ddf_high - ddf_low;
          const double local_slope = rise / rdisc;
          const double dr = -ddf_low / local_slope;
          quintic_dd_roots.push_back(r + dr);
        }
      }
      for (size_t i = 0; i < quintic_dd_roots.size(); i++) {
        const double r = quintic_dd_roots[i];
        const double qnfd_val = (((((((qnfd_a * r) + qnfd_b) * r) + qnfd_c) * r) + qnfd_d) * r) +
                                qnfd_e;
        if (qnfd_val > 1.0e-6) {
          quintic_slope_always_negative = false;
        }
      }

      // Add the ratio of the quintic function's maximum height (which, if the previous tests have
      // passed, occurs at zero) to the original Lennard-Jones well depth.
      if (ij_sigma > 1.0e-6) {
        all_ij5_zratio.push_back(ef_coefs.y / ij_epsilon);
      }
    }
  }
  check(cubic_slope_always_negative, "The cubic softcore potentials do not always deliver a "
        "negative slope when taking over from a Lennard-Jones function.", do_tests);
  check(quartic_slope_always_negative, "The quartic softcore potentials do not always deliver a "
        "negative slope when taking over from a Lennard-Jones function.", do_tests);
  check(quintic_slope_always_negative, "The quintic softcore potentials do not always deliver a "
        "negative slope when taking over from a Lennard-Jones function.", do_tests);
  switch (do_tests) {
  case TestPriority::CRITICAL:
    break;
  case TestPriority::NON_CRITICAL:
  case TestPriority::ABORT:
    if (all_ij3_zratio.size() == 0) {
      all_ij3_zratio = std::vector(3, 0.0);
    }
    if (all_ij4_zratio.size() == 0) {
      all_ij4_zratio = std::vector(3, 0.0);
    }
    if (all_ij5_zratio.size() == 0) {
      all_ij5_zratio = std::vector(3, 0.0);
    }
    break;
  }
  check(mean(all_ij3_zratio), RelationalOperator::EQUAL, zero_ratio3, "The ratio of the maximum "
        "function height to the Lennard-Jones well depth in a cubic soft core potential does "
        "not meet expectations.", do_tests);
  check(mean(all_ij4_zratio), RelationalOperator::EQUAL, zero_ratio4, "The ratio of the maximum "
        "function height to the Lennard-Jones well depth in a quartic soft core potential does "
        "not meet expectations.", do_tests);
  check(mean(all_ij5_zratio), RelationalOperator::EQUAL, zero_ratio5, "The ratio of the maximum "
        "function height to the Lennard-Jones well depth in a quintic soft core potential does "
        "not meet expectations.", do_tests);
  check(variance(all_ij3_zratio, VarianceMethod::STANDARD_DEVIATION), RelationalOperator::EQUAL,
        Approx(0.0).margin(1.0e-5), "The ratio of the maximum function height to the "
        "Lennard-Jones minimum is not entirely consistent for rswitch = " +
        realToString(sigma_factor, 9, 4, NumberFormat::STANDARD_REAL) + " * sigma in a cubic "
        "soft core potential.", do_tests);
  check(variance(all_ij4_zratio, VarianceMethod::STANDARD_DEVIATION), RelationalOperator::EQUAL,
        Approx(0.0).margin(1.0e-5), "The ratio of the maximum function height to the "
        "Lennard-Jones minimum is not entirely consistent for rswitch = " +
        realToString(sigma_factor, 9, 4, NumberFormat::STANDARD_REAL) + " * sigma in a quartic "
        "soft core potential.", do_tests);
  check(variance(all_ij5_zratio, VarianceMethod::STANDARD_DEVIATION), RelationalOperator::EQUAL,
        Approx(0.0).margin(1.0e-5), "The ratio of the maximum function height to the "
        "Lennard-Jones minimum is not entirely consistent for rswitch = " +
        realToString(sigma_factor, 9, 4, NumberFormat::STANDARD_REAL) + " * sigma in a quintic "
        "soft core potential.", do_tests);
}

//-------------------------------------------------------------------------------------------------
// Run tests of the benchmark PME direct space interaction functions.
//-------------------------------------------------------------------------------------------------
void runPMEDirectSpaceTests() {
  const double ew_coeff = ewaldCoefficient(12.0, 1.0e-7);
  const double kcoul = stormm::symbols::charmm_gromacs_bioq;

  // Run direct space one-dimensional tests
  const int npts = 128;
  std::vector<double> du_eval(npts), du_fd(npts);
  std::vector<double> d2u_eval(npts), d2u_fd(npts), d3u_eval(npts), d3u_fd(npts);
  const double fd_dr = pow(2.0, -18.0);
  double r = 0.9;
  for (int i = 0; i < npts; i++) {
    const double r_p = r + fd_dr;
    const double r_n = r - fd_dr;

    // First derivative
    du_eval[i] = elecPMEDirectSpace(ew_coeff, kcoul, r, r * r, 1);
    const double u_p = kcoul * erfc(ew_coeff * r_p) / r_p;
    const double u_n = kcoul * erfc(ew_coeff * r_n) / r_n;
    du_fd[i] = (u_p - u_n) / (2.0 * fd_dr);

    // Second derivative
    d2u_eval[i] = elecPMEDirectSpace(ew_coeff, kcoul, r, r * r, 2);
    const double du_p = elecPMEDirectSpace(ew_coeff, kcoul, r_p, r_p * r_p, 1);
    const double du_n = elecPMEDirectSpace(ew_coeff, kcoul, r_n, r_n * r_n, 1);
    d2u_fd[i] = (du_p - du_n) / (2.0 * fd_dr);

    // Third derivative
    d3u_eval[i] = elecPMEDirectSpace(ew_coeff, kcoul, r, r * r, 3);
    const double d2u_p = elecPMEDirectSpace(ew_coeff, kcoul, r_p, r_p * r_p, 2);
    const double d2u_n = elecPMEDirectSpace(ew_coeff, kcoul, r_n, r_n * r_n, 2);
    d3u_fd[i] = (d2u_p - d2u_n) / (2.0 * fd_dr);

    // Increment the range variable
    r += 0.1;
  }
  check(du_eval, RelationalOperator::EQUAL, du_fd, "The first derivative of the electrostatic PME "
        "direct-space interaction was not computed as expected.");
  check(d2u_eval, RelationalOperator::EQUAL, d2u_fd, "The second derivative of the electrostatic "
        "PME direct-space interaction was not computed as expected.");
  check(d3u_eval, RelationalOperator::EQUAL, d3u_fd, "The third derivative of the electrostatic "
        "PME direct-space interaction was not computed as expected.");
}

//-------------------------------------------------------------------------------------------------
// Run tests of the layered potential decomposition in the Coulomb potential.
//
// Arguments:
//   oe:                Testing environment details, including the path to the STORMM ${tmpdir}
//   u_type:            The type of potential to decompose
//   prt_prt_cutoff:    Cutoff between particle-particle interactions
//   first_pass_write:  Indicate how to open the snapshot file if it needs to be revised
//-------------------------------------------------------------------------------------------------
void runLayeredCoulombTests(const TestEnvironment &oe, const DecomposablePotential u_type,
                            const double prt_prt_cutoff = 5.4,
                            const PrintSituation first_pass_write = PrintSituation::APPEND) {
  LayeredPotentialMetrics lpm;
  lpm.setForm(u_type);
  lpm.setBoundaryCondition(BoundaryCondition::ISOLATED);
  lpm.setCutoff(prt_prt_cutoff);
  lpm.setRangeCompounding(2.0);
  lpm.setMaximumRange(100.0);
  lpm.setEwaldCoefficient(ewaldCoefficient(4.0 * prt_prt_cutoff, 1.0e-10));

  // Perform general checks on the management class for only one of the potential types--these are
  // agnostic to the form.
  if (u_type == DecomposablePotential::ELECTROSTATIC) {
    check(lpm.getCutoff(0), RelationalOperator::EQUAL, 5.4, "The LayeredPotentialMetrics object "
          "did not return the expected value of the cutoff.");
    check(lpm.getCutoff(1), RelationalOperator::EQUAL, 10.8, "The LayeredPotentialMetrics object "
          "did not return the expected value of the cutoff.");
    check(lpm.getCutoff(4), RelationalOperator::EQUAL, 86.4, "The LayeredPotentialMetrics object "
          "did not return the expected value of the cutoff.");
    const std::vector<double> etest = { 1.25, 1.75, 2.25 };
    lpm.setExponentFactor(etest);
    std::vector<double> rpt_efac(9);
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        rpt_efac[(3 * j) + i] = lpm.getExponentFactor(i, j);        
      }
    }
    const std::vector<double> chk_efac = { 0.000, 0.000, 0.000,
                                           1.250, 1.750, 2.250,
                                           0.625, 0.875, 1.125 };
    check(rpt_efac, RelationalOperator::EQUAL, chk_efac, "The exponential factors for a layereed "
          "potential are not computed as expected.");
  }
  LayeredPotential<double, double4> lnrg(lpm);
  const double4 layer_a_smc = lnrg.getSmoothingCoefficients(0);
  const double4 layer_b_smc = lnrg.getSmoothingCoefficients(1);
  const double4 layer_d_smc = lnrg.getSmoothingCoefficients(3);
  const std::vector<double> rpt_coef = { layer_a_smc.x, layer_a_smc.y, layer_a_smc.z,
                                         layer_a_smc.w, layer_b_smc.x, layer_b_smc.y,
                                         layer_b_smc.z, layer_b_smc.w, layer_d_smc.x,
                                         layer_d_smc.y, layer_d_smc.z, layer_d_smc.w };
  std::vector<double> chk_coef;
  switch (u_type) {
  case DecomposablePotential::ELECTROSTATIC:
    chk_coef = {    0.00000000,    0.00000000,    0.00000000,    0.00000000,
                 -102.48056614,  112.29930802,  -35.47125849,   87.14366796,
                  -25.62014153,   28.07482701,   -8.86781462,   21.78591699 };
    break;
  case DecomposablePotential::DISPERSION:
    chk_coef = {    0.00000000,    0.00000000,    0.00000000,    0.00000000,
                    0.00084762,   -0.00093321,    0.00029850,   -0.00025324,
                    0.00000021,   -0.00000023,    0.00000007,   -0.00000006 };
    break;
  case DecomposablePotential::ELEC_PME_DIRECT:
    chk_coef = {    0.00000000,    0.00000000,    0.00000000,    0.00000000,
                  -69.52872751,   71.38184892,  -21.75725140,   28.17895772,
                   -0.00003361,    0.00004208,   -0.00001487,    0.00000643 };
    break;
  }
  check(rpt_coef, RelationalOperator::EQUAL, chk_coef, "The fitted coefficients for various "
        "layers of the " + getEnumerationName(u_type) + " potential decomposition do not meet "
        "expectations.");
  const int npts = 500;
  std::vector<double> base_layer(npts), first_layer(npts), second_layer(npts), third_layer(npts);
  double rp = 0.5;
  for (int i = 0; i < npts; i++) {
    base_layer[i] = lnrg.getAnalyticValue(0, rp);
    first_layer[i] = lnrg.getAnalyticValue(1, rp);
    second_layer[i] = lnrg.getAnalyticValue(2, rp);
    third_layer[i] = lnrg.getAnalyticValue(3, rp);
    rp += 0.1;
  }
  const char osc = osSeparator();
  const std::string snp_name = oe.getStormmSourcePath() + osc + "test" + osc + "Potential" + osc +
                               "potential_decomp.m";
  const bool snp_exists = (getDrivePathType(snp_name) == DrivePathType::FILE);
  const TestPriority do_snps = (snp_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (snp_exists == false && oe.takeSnapshot() != SnapshotOperation::SNAPSHOT) {
    rtWarn("The snapshot file " + snp_name + " could not be found.  Check the STORMM source "
           "path for validity.", "runHAILCoulombicTests");
  }
  std::string var_ext;
  switch (u_type) {
  case DecomposablePotential::ELECTROSTATIC:
    var_ext = "_elec";
    break;
  case DecomposablePotential::DISPERSION:
    var_ext = "_disp";
    break;
  case DecomposablePotential::ELEC_PME_DIRECT:
    var_ext = "_epme";
    break;
  }
  snapshot(snp_name, polyNumericVector(base_layer), "base" + var_ext , 1.0e-6, "Values of the "
           "decomposed " + getEnumerationName(u_type) + " particle-particle interaction do not "
           "meet expectations.", oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL,
           first_pass_write, do_snps);
  snapshot(snp_name, polyNumericVector(first_layer), "first" + var_ext, 1.0e-6, "Values of the "
           "decomposed " + getEnumerationName(u_type) + " first mesh-based potential do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL,
           PrintSituation::APPEND, do_snps);
  snapshot(snp_name, polyNumericVector(second_layer), "second" + var_ext, 1.0e-6, "Values of the "
           "decomposed " + getEnumerationName(u_type) + " second mesh-based potential do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL,
           PrintSituation::APPEND, do_snps);
  snapshot(snp_name, polyNumericVector(third_layer), "third" + var_ext, 1.0e-6, "Values of the "
           "decomposed " + getEnumerationName(u_type) + " third mesh-based potential do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-8, NumberFormat::STANDARD_REAL,
           PrintSituation::APPEND, do_snps);

  // Check derivatives of the potential by successive finite difference approximations.
  const std::vector<std::string> deriv_names = { "first", "second", "third" };
  for (int dlev = 1; dlev < 4; dlev++) {
    std::vector<double> base_dlyr(npts), first_dlyr(npts), second_dlyr(npts), third_dlyr(npts);
    std::vector<double> base_fd(npts), first_fd(npts), second_fd(npts), third_fd(npts);
    const double delta = pow(2.0, -20.0);
    rp = 0.5;
    for (int i = 0; i < npts; i++) {
      base_dlyr[i] = lnrg.getAnalyticDerivative(0, rp, rp * rp, dlev);
      first_dlyr[i] = lnrg.getAnalyticDerivative(1, rp, rp * rp, dlev);
      second_dlyr[i] = lnrg.getAnalyticDerivative(2, rp, rp * rp, dlev);
      third_dlyr[i] = lnrg.getAnalyticDerivative(3, rp, rp * rp, dlev);
      const double rp_p = rp + delta;
      const double rp_n = rp - delta;
      const double base_p = lnrg.getAnalyticDerivative(0, rp_p, rp_p * rp_p, dlev - 1);
      const double base_n = lnrg.getAnalyticDerivative(0, rp_n, rp_n * rp_n, dlev - 1);
      base_fd[i] = (base_p - base_n) / (2.0 * delta);
      const double first_p = lnrg.getAnalyticDerivative(1, rp_p, rp_p * rp_p, dlev - 1);
      const double first_n = lnrg.getAnalyticDerivative(1, rp_n, rp_n * rp_n, dlev - 1);
      first_fd[i] = (first_p - first_n) / (2.0 * delta);
      const double second_p = lnrg.getAnalyticDerivative(2, rp_p, rp_p * rp_p, dlev - 1);
      const double second_n = lnrg.getAnalyticDerivative(2, rp_n, rp_n * rp_n, dlev - 1);
      second_fd[i] = (second_p - second_n) / (2.0 * delta);
      const double third_p = lnrg.getAnalyticDerivative(3, rp_p, rp_p * rp_p, dlev - 1);
      const double third_n = lnrg.getAnalyticDerivative(3, rp_n, rp_n * rp_n, dlev - 1);
      third_fd[i] = (third_p - third_n) / (2.0 * delta);
    }
    double dthr_tol;
    switch (u_type) {
    case DecomposablePotential::ELECTROSTATIC:
      dthr_tol = 3.0e-6;
      break;
    case DecomposablePotential::DISPERSION:
      dthr_tol = 1.2e-5;
      break;
    case DecomposablePotential::ELEC_PME_DIRECT:
      dthr_tol = 3.0e-6;
      break;
    }
    check(base_dlyr, RelationalOperator::EQUAL, Approx(base_fd).margin(dthr_tol), "The " +
          deriv_names[dlev - 1] + " derivative of the decomposed " + getEnumerationName(u_type) +
          " potential does not agree with a finite difference approximation in the "
          "particle-particle layer.");
    check(first_dlyr, RelationalOperator::EQUAL, first_fd, "The " + deriv_names[dlev - 1] +
          " derivative of the decomposed " + getEnumerationName(u_type) + " potential does not "
          "agree with a finite difference approximation in the first mesh-based layer.");
    check(second_dlyr, RelationalOperator::EQUAL, second_fd, "The " + deriv_names[dlev - 1] +
          " derivative of the decomposed " + getEnumerationName(u_type) + " potential does not "
          "agree with a finite difference approximation in the second mesh-based layer.");
    check(third_dlyr, RelationalOperator::EQUAL, third_fd, "The " + deriv_names[dlev - 1] +
          " derivative of the decomposed " + getEnumerationName(u_type) + " potential does not "
          "agree with a finite difference approximation in the third mesh-based layer.");
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
  StopWatch timer;
  const int input_timings  = timer.addCategory("Input file parsing");
  const int sem_timings    = timer.addCategory("Static exclusion mask");
  const int attn14_timings = timer.addCategory("Compute 1:4 electrostatics");
  const int nb_timings     = timer.addCategory("Compute non-bonded interactions");
  
  // Section 1
  section("Tests of the StaticExclusionMask object");

  // Section 2
  section("Non-bonded 1:4 energy evaluation");

  // Section 3
  section("Non-bonded energy evaluation");

  // Section 4
  section("Non-bonded all-to-all force evaluation");

  // Section 5
  section("Fixed-precision energy accumulator checks");

  // Section 6
  section("Single-precision coordinate performance");

  // Section 7
  section("Fixed-precision coordinate performance");

  // Section 8
  section("Energy tracking object validation");

  // Section 9
  section("Soft core potential evaluation");

  // Section 10
  section("PME direct space interactions");

  // Section 11
  section("Layered potential decomposition");

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
  const std::string trpw_top_name = base_top_name + osc + "trpcage_in_water.top";
  const std::string trpw_crd_name = base_crd_name + osc + "trpcage_in_water.inpcrd";
  const std::string trpp_top_name = base_top_name + osc + "trpcage_no_z.top";
  const bool systems_exist = (getDrivePathType(trpi_top_name) == DrivePathType::FILE &&
                              getDrivePathType(trpi_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(dhfr_top_name) == DrivePathType::FILE &&
                              getDrivePathType(dhfr_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(alad_top_name) == DrivePathType::FILE &&
                              getDrivePathType(alad_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(trpw_top_name) == DrivePathType::FILE &&
                              getDrivePathType(trpw_crd_name) == DrivePathType::FILE &&
                              getDrivePathType(trpw_top_name) == DrivePathType::FILE);
  const TestPriority do_tests = (systems_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (systems_exist == false) {
    rtWarn("Files for the Trp-cage miniprotein (with and without water), the DHFR globular "
           "protein (with CHARMM potential details), and alanine dipeptide (with ff19SB) were not "
           "found.  These files should be found in the ${STORMM_SOURCE}/test/Topology and "
           "${STORMM_SOURCE}/test/Trajectory directories.  Check the $STORMM_SOURCE environment "
           "variable.  A number of tests will be skipped.", "test_nonbonded_evaluation");
  }

  // Read topologies and coordinates
  AtomGraph trpi_ag, dhfr_ag, alad_ag, trpw_ag, trpp_ag;
  PhaseSpace trpi_ps, dhfr_ps, alad_ps, trpw_ps;
  if (systems_exist) {
    trpi_ag.buildFromPrmtop(trpi_top_name, ExceptionResponse::SILENT);
    trpi_ps.buildFromFile(trpi_crd_name, CoordinateFileKind::AMBER_INPCRD);
    dhfr_ag.buildFromPrmtop(dhfr_top_name, ExceptionResponse::SILENT);
    dhfr_ps.buildFromFile(dhfr_crd_name, CoordinateFileKind::AMBER_INPCRD);
    alad_ag.buildFromPrmtop(alad_top_name, ExceptionResponse::SILENT);
    alad_ps.buildFromFile(alad_crd_name, CoordinateFileKind::AMBER_INPCRD);
    trpw_ag.buildFromPrmtop(trpw_top_name, ExceptionResponse::SILENT);
    trpw_ps.buildFromFile(trpw_crd_name, CoordinateFileKind::AMBER_INPCRD);
    trpp_ag.buildFromPrmtop(trpw_top_name, ExceptionResponse::SILENT);
  }
  timer.assignTime(input_timings);
  
  // Prepare exclusion masks for each system (the constructor will work even if the
  // system has zero atoms)
  const StaticExclusionMask trpi_semask(&trpi_ag);
  const StaticExclusionMask dhfr_semask(&dhfr_ag);
  const StaticExclusionMask alad_semask(&alad_ag);
  const StaticExclusionMask trpw_semask(&trpw_ag);
  timer.assignTime(sem_timings);

  // Check exclusions for three systems against the original topologies
  section(1);
  const std::vector<const StaticExclusionMask*> my_masks = { &trpi_semask, &alad_semask,
                                                             &trpw_semask };
  const std::vector<std::string> my_systems = { "Trp-cage (unsolvated)", "alanine dipeptide",
                                                "Trp-cage (solvated)" };
  for (size_t i = 0; i < my_masks.size(); i++) {
    const int3 n_errors = checkMarkedExclusions(*(my_masks[i]));
    check(n_errors.x == 0, "The StaticExclusionMask object recorded " +
          std::to_string(n_errors.x) + " exclusions along the diagonal in the " + my_systems[i] +
          " system.  These should be handled implicitly by any non-bonded loop.", do_tests);
    check(n_errors.y == 0, "The StaticExclusionMask object failed to record " +
          std::to_string(n_errors.y) + " exclusions marked in the topology for the " +
          my_systems[i] + " system.", do_tests);
    check(n_errors.z == 0, "The StaticExclusionMask object recorded " +
          std::to_string(n_errors.z) + " exclusions not in the original topology for the " +
          my_systems[i] + " system.", do_tests);
  }

  // Check the reflexive nature of the exclusion mask tiles
  for (size_t i = 0; i < my_masks.size(); i++) {
    const int2 n_errors = checkReflexiveMarks(*(my_masks[i]));
    check(n_errors.x == 0, "The StaticExclusionMask object recorded " +
          std::to_string(n_errors.x) + " exclusions in the sending atoms but not in the receiving "
          "atoms for the " + my_systems[i] + " system.", do_tests);
    check(n_errors.y == 0, "The StaticExclusionMask object recorded " +
          std::to_string(n_errors.x) + " exclusions in the receiving atoms but not in the sending "
          "atoms for the " + my_systems[i] + " system.", do_tests);
  }

  // Accumulate energies and forces
  ScoreCard all_systems_sc(5), secondary_sc(5);
  const int trpi_idx = 0;
  const int dhfr_idx = 1;
  const int alad_idx = 2;
  const int trpw_idx = 3;
  const int trpp_idx = 4;
  TrajectoryKind tkind = TrajectoryKind::FORCES;

  // Check for the existence of snapshot files
  const std::string trpi_snapshot(base_ptl_name + osc + "trpcage_nb_details.m");
  const std::string dhfr_snapshot(base_ptl_name + osc + "dhfr_nb_details.m");
  const std::string alad_snapshot(base_ptl_name + osc + "ala_dipeptide_nb_details.m");
  const bool snps_exist = (getDrivePathType(trpi_snapshot) == DrivePathType::FILE &&
                           getDrivePathType(dhfr_snapshot) == DrivePathType::FILE &&
                           getDrivePathType(alad_snapshot) == DrivePathType::FILE);
  const TestPriority snap_check = (snps_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (snps_exist == false && oe.takeSnapshot() != SnapshotOperation::SNAPSHOT) {
    rtWarn("Snapshot files " + trpi_snapshot + ", " + dhfr_snapshot + ", and " + alad_snapshot +
           " were not found.  These files contain reference forces for checking the valence "
           "energy derivative calculations.  Check that the ${STORMM_SOURCE} environment variable "
           "is set properly so that these snapshot files may be found in "
           "${STORMM_SOURCE}/test/Potential/.  Subsequent tests will be skipped until these "
           "reference files are available.", "test_nonbonded_evaluation");
  }

  // Compute 1:4 electrostatic and Lennard-Jones forces on the isolated Trp-cage system
  section(2);
  trpi_ps.initializeForces();
  timer.assignTime(0);
  const double2 trpi_14_e = evaluateAttenuated14Terms(trpi_ag, &trpi_ps, &all_systems_sc,
                                                      EvaluateForce::YES, EvaluateForce::NO,
                                                      trpi_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> trpi_14_elec_frc = trpi_ps.getInterlacedCoordinates(tkind);
  trpi_ps.initializeForces();
  timer.assignTime(0);
  evaluateAttenuated14Terms(trpi_ag, &trpi_ps, &secondary_sc, EvaluateForce::NO,
                            EvaluateForce::YES, trpi_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> trpi_14_vdw_frc = trpi_ps.getInterlacedCoordinates(tkind);
  check(trpi_14_e.x, RelationalOperator::EQUAL, Approx(1458.0998129).margin(1.0e-6),
        "Electrostatic 1:4 energy for Trp-cage (AMBER ff99SB force field) was not computed "
        "correctly.", do_tests);
  check(trpi_14_e.y, RelationalOperator::EQUAL, Approx(62.6481811).margin(1.0e-6), "van-der Waals "
        "1:4 energy for Trp-cage (AMBER ff99SB force field) was not computed correctly.",
        do_tests);
  snapshot(trpi_snapshot, polyNumericVector(trpi_14_elec_frc), "trpcage_qq14_forces",
           NumberFormat::SCIENTIFIC, "Forces due to electrostatic 1:4 non-bonded interactions in "
           "the Trp-cage (isolated boundary conditions) system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::OVERWRITE, snap_check);
  snapshot(trpi_snapshot, polyNumericVector(trpi_14_vdw_frc), "trpcage_lj14_forces",
           NumberFormat::SCIENTIFIC, "Forces due to Lennard-Jones 1:4 non-bonded interactions in "
           "the Trp-cage (isolated boundary conditions) system do not meet expectations.",
           oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::APPEND, snap_check);
  timer.assignTime(0);
  const double2 trpi_14_e_ii = evaluateAttenuated14Terms(trpi_ag, CoordinateFrame(trpi_ps),
                                                         &secondary_sc, trpi_idx);
  timer.assignTime(attn14_timings);
  check(trpi_14_e_ii.x, RelationalOperator::EQUAL, Approx(1458.0998129).margin(1.0e-6),
        "Electrostatic 1:4 energy for Trp-cage (AMBER ff99SB force field) was not computed "
        "correctly when evaluating energy only with a coordinate frame object.", do_tests);
  check(trpi_14_e_ii.y, RelationalOperator::EQUAL, Approx(62.6481811).margin(1.0e-6),
        "van-der Waals 1:4 energy for Trp-cage (AMBER ff99SB force field) was not computed "
        "correctly when evaluating energy only with a coordinate frame object.", do_tests);
  
  // Compute 1:4 electrostatic and Lennard-Jones forces on the DHFR system
  dhfr_ps.initializeForces();
  timer.assignTime(0);
  const double2 dhfr_14_e = evaluateAttenuated14Terms(dhfr_ag, &dhfr_ps, &all_systems_sc,
                                                      EvaluateForce::YES, EvaluateForce::NO,
                                                      dhfr_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> dhfr_14_elec_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  dhfr_ps.initializeForces();
  timer.assignTime(0);
  evaluateAttenuated14Terms(dhfr_ag, &dhfr_ps, &secondary_sc, EvaluateForce::NO,
                            EvaluateForce::YES, dhfr_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> dhfr_14_vdw_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  check(dhfr_14_e.x, RelationalOperator::EQUAL, Approx(6507.3376751).margin(1.0e-6),
        "Electrostatic 1:4 energy for DHFR (CHARMM force field) was not computed correctly.",
        do_tests);
  check(dhfr_14_e.y, RelationalOperator::EQUAL, Approx(367.0925927).margin(1.0e-6),
        "van-der Waals 1:4 energy for DHFR (CHARMM force field, with special 1:4 LJ coefficients) "
        "was not computed correctly.", do_tests);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_14_elec_frc), "dhfr_qq14_forces",
           NumberFormat::SCIENTIFIC, "Forces due to electrostatic 1:4 non-bonded interactions in "
           "the DHFR system do not meet expectations.", oe.takeSnapshot(), 1.0e-6, 1.0e-12,
           PrintSituation::OVERWRITE, snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_14_vdw_frc), "dhfr_lj14_forces",
           NumberFormat::SCIENTIFIC, "Forces due to Lennard-Jones 1:4 non-bonded interactions in "
           "the DHFR system do not meet expectations.", oe.takeSnapshot(), 1.0e-6, 1.0e-12,
           PrintSituation::APPEND, snap_check);

  // Compute 1:4 electrostatic and Lennard-Jones forces on the alanine dipeptide system
  alad_ps.initializeForces();
  timer.assignTime(0);
  const double2 alad_14_e = evaluateAttenuated14Terms(alad_ag, &alad_ps, &all_systems_sc,
                                                      EvaluateForce::YES, EvaluateForce::NO,
                                                      alad_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> alad_14_elec_frc = alad_ps.getInterlacedCoordinates(tkind);
  alad_ps.initializeForces();
  timer.assignTime(0);
  evaluateAttenuated14Terms(alad_ag, &alad_ps, &secondary_sc, EvaluateForce::NO,
                            EvaluateForce::YES, alad_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> alad_14_vdw_frc = alad_ps.getInterlacedCoordinates(tkind);
  check(alad_14_e.x, RelationalOperator::EQUAL, Approx(46.8072425).margin(1.0e-6), "Electrostatic "
        "1:4 energy for alanine dipeptide (ff19SB force field) was not computed correctly.",
        do_tests);
  check(alad_14_e.y, RelationalOperator::EQUAL, Approx(3.1589264).margin(1.0e-6), "van-der Waals "
        "1:4 energy for alanine dipeptide (ff19SB force field) was not computed correctly.",
        do_tests);
  snapshot(alad_snapshot, polyNumericVector(alad_14_elec_frc), "ala_dipeptide_qq14_forces",
           NumberFormat::SCIENTIFIC, "Forces due to electrostatic 1:4 non-bonded interactions in "
           "the Ala dipeptide system do not meet expectations.", oe.takeSnapshot(), 1.0e-6,
           1.0e-12, PrintSituation::OVERWRITE, snap_check);
  snapshot(alad_snapshot, polyNumericVector(alad_14_vdw_frc), "ala_dipeptide_lj14_forces",
           NumberFormat::SCIENTIFIC, "Forces due to Lennard-Jones 1:4 non-bonded interactions in "
           "the Ala dipeptide system do not meet expectations.", oe.takeSnapshot(), 1.0e-6,
           1.0e-12, PrintSituation::APPEND, snap_check);

  // Compute 1:4 electrostatic and Lennard-Jones forces on the solvated Trp-cage system
  trpw_ps.initializeForces();
  timer.assignTime(0);
  const double2 trpw_14_e = evaluateAttenuated14Terms(trpw_ag, &trpw_ps, &all_systems_sc,
                                                      EvaluateForce::YES, EvaluateForce::NO,
                                                      trpw_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> trpw_14_elec_frc = trpw_ps.getInterlacedCoordinates(tkind);
  trpw_ps.initializeForces();
  timer.assignTime(0);
  evaluateAttenuated14Terms(trpw_ag, &trpw_ps, &secondary_sc, EvaluateForce::YES,
                            EvaluateForce::NO, trpw_idx);
  timer.assignTime(attn14_timings);
  const std::vector<double> trpw_14_vdw_frc = trpw_ps.getInterlacedCoordinates(tkind);
  check(trpw_14_e.x, RelationalOperator::EQUAL, Approx(1458.0996855).margin(1.0e-6),
        "Electrostatic 1:4 energy for solvated Trp-cage (AMBER ff99SB force field) was not "
        "computed correctly.", do_tests);
  check(trpw_14_e.y, RelationalOperator::EQUAL, Approx(62.6481797).margin(1.0e-6), "van-der Waals "
        "1:4 energy for solvated Trp-cage (AMBER ff99SB force field) was not computed correctly.",
        do_tests);
  timer.assignTime(0);
  const double2 trpw_14_e_ii = evaluateAttenuated14Terms(trpw_ag, CoordinateFrame(trpw_ps),
                                                         &all_systems_sc, trpw_idx);
  timer.assignTime(attn14_timings);
  check(trpw_14_e_ii.x, RelationalOperator::EQUAL, Approx(1458.0996855).margin(1.0e-6),
        "Electrostatic 1:4 energy for solvated Trp-cage (AMBER ff99SB force field) was not "
        "computed correctly when evaluating energy only with a CoordinateFrame abstract.",
        do_tests);
  check(trpw_14_e_ii.y, RelationalOperator::EQUAL, Approx(62.6481797).margin(1.0e-6),
        "van-der Waals 1:4 energy for solvated Trp-cage (AMBER ff99SB force field) was not "
        "computed correctly when evaluating energy only with a CoordinateFrame abstract.",
        do_tests);

  // Compute 1:4 electrostatic energies on the Trp-cage system with no explicit scaling factors in
  // the topology
  trpi_ps.initializeForces();
  timer.assignTime(0);
  const double2 trpp_14_e = evaluateAttenuated14Terms(trpp_ag, &trpi_ps, &all_systems_sc,
                                                      EvaluateForce::YES, EvaluateForce::NO,
                                                      trpp_idx);
  timer.assignTime(attn14_timings);
  trpi_ps.initializeForces();
  check(trpp_14_e.x, RelationalOperator::EQUAL, Approx(1458.0996843).margin(1.0e-6),
        "Electrostatic 1:4 energy for Trp-cage (AMBER ff99SB force field) was not computed "
        "correctly.", do_tests);
  check(trpp_14_e.y, RelationalOperator::EQUAL, Approx(62.6481811).margin(1.0e-6), "van-der Waals "
        "1:4 energy for Trp-cage (AMBER ff99SB force field) was not computed correctly.",
        do_tests);

  // Compute non-bonded energy for the isolated Trp-cage system
  section(3);
  trpi_ps.initializeForces();
  timer.assignTime(0);
  const double2 trpi_nonb_e = evaluateNonbondedEnergy(trpi_ag, trpi_semask, &trpi_ps,
                                                      &all_systems_sc, EvaluateForce::YES,
                                                      EvaluateForce::NO, trpi_idx);
  timer.assignTime(nb_timings);
  const std::vector<double> trpi_elec_frc = trpi_ps.getInterlacedCoordinates(tkind);
  trpi_ps.initializeForces();
  timer.assignTime(0);
  evaluateNonbondedEnergy(trpi_ag, trpi_semask, &trpi_ps, &secondary_sc, EvaluateForce::NO,
                          EvaluateForce::YES, trpi_idx);
  timer.assignTime(nb_timings);
  const std::vector<double> trpi_vdw_frc = trpi_ps.getInterlacedCoordinates(tkind);
  check(trpi_nonb_e.x, RelationalOperator::EQUAL, Approx(-1816.1478854).margin(1.0e-6),
        "Electrostatic energy for isolated Trp-cage (AMBER ff99SB force field) was not computed "
        "correctly.", do_tests);
  check(trpi_nonb_e.y, RelationalOperator::EQUAL, Approx(-119.3500215).margin(1.0e-6),
        "van-der Waals energy for isolated Trp-cage (AMBER ff99SB force field) was not computed "
        "correctly.", do_tests);

  // Compute non-bonded energy for the DHFR system
  dhfr_ps.initializeForces();
  timer.assignTime(0);
  const double2 dhfr_nonb_e = evaluateNonbondedEnergy(dhfr_ag, dhfr_semask, &dhfr_ps,
                                                      &all_systems_sc, EvaluateForce::YES,
                                                      EvaluateForce::NO, dhfr_idx);
  timer.assignTime(nb_timings);
  const std::vector<double> dhfr_elec_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  dhfr_ps.initializeForces();
  timer.assignTime(0);
  evaluateNonbondedEnergy(dhfr_ag, dhfr_semask, &dhfr_ps, &secondary_sc, EvaluateForce::NO,
                          EvaluateForce::YES, dhfr_idx);
  timer.assignTime(nb_timings);
  const std::vector<double> dhfr_vdw_frc = dhfr_ps.getInterlacedCoordinates(tkind);
  check(dhfr_nonb_e.x, RelationalOperator::EQUAL, Approx(-10036.4655324).margin(1.0e-6),
        "Electrostatic energy for DHFR (CHARMM force field) was not computed correctly.",
        do_tests);
  check(dhfr_nonb_e.y, RelationalOperator::EQUAL, Approx(-1009.1558078).margin(1.0e-6),
        "van-der Waals energy for DHFR (CHARMM force field) was not computed correctly.",
        do_tests);
  timer.assignTime(0);
  const double2 dhfr_nonb_e_ii = evaluateNonbondedEnergy(dhfr_ag, dhfr_semask,
                                                         CoordinateFrame(dhfr_ps),
                                                         &secondary_sc, dhfr_idx);
  timer.assignTime(nb_timings);
  check(dhfr_nonb_e_ii.x, RelationalOperator::EQUAL, Approx(-10036.4655324).margin(1.0e-6),
        "Electrostatic energy for DHFR (CHARMM force field) was not computed correctly when "
        "evaluating energy only with a CoordinateFrame abstract.", do_tests);
  check(dhfr_nonb_e_ii.y, RelationalOperator::EQUAL, Approx(-1009.1558078).margin(1.0e-6),
        "van-der Waals energy for DHFR (CHARMM force field) was not computed correctly when "
        "evaluating energy only with a CoordinateFrame abstract.", do_tests);

  // Compute non-bonded energy for the alanine dipeptide system
  alad_ps.initializeForces();
  timer.assignTime(0);
  const double2 alad_nonb_e = evaluateNonbondedEnergy(alad_ag, alad_semask, &alad_ps,
                                                      &all_systems_sc, EvaluateForce::YES,
                                                      EvaluateForce::NO, alad_idx);
  timer.assignTime(nb_timings);
  const std::vector<double> alad_elec_frc = alad_ps.getInterlacedCoordinates(tkind);
  alad_ps.initializeForces();
  timer.assignTime(0);
  evaluateNonbondedEnergy(alad_ag, alad_semask, &alad_ps, &secondary_sc, EvaluateForce::NO,
                          EvaluateForce::YES, alad_idx);
  timer.assignTime(nb_timings);
  const std::vector<double> alad_vdw_frc = alad_ps.getInterlacedCoordinates(tkind);
  check(alad_nonb_e.x, RelationalOperator::EQUAL, Approx(-78.8463310).margin(1.0e-6),
        "Electrostatic energy for alanine dipeptide (ff19SB force field) was not computed "
        "correctly.", do_tests);
  check(alad_nonb_e.y, RelationalOperator::EQUAL, Approx(-1.2452480).margin(1.0e-6),
        "van-der Waals energy for alanine dipeptide (ff19SB force field) was not computed "
        "correctly.", do_tests);
  timer.assignTime(0);
  const double2 alad_nonb_e_ii = evaluateNonbondedEnergy(alad_ag, alad_semask,
                                                         CoordinateFrame(alad_ps),
                                                         &secondary_sc, alad_idx);
  timer.assignTime(nb_timings);
  check(alad_nonb_e_ii.x, RelationalOperator::EQUAL, Approx(-78.8463310).margin(1.0e-6),
        "Electrostatic energy for alanine dipeptide (ff19SB force field) was not computed "
        "correctly when evaluating energy only with a CoordinateFrame abstract.", do_tests);
  check(alad_nonb_e_ii.y, RelationalOperator::EQUAL, Approx(-1.2452480).margin(1.0e-6),
        "van-der Waals energy for alanine dipeptide (ff19SB force field) was not computed "
        "correctly when evaluating energy only with a CoordinateFrame abstract.", do_tests);  

  // Check forces for the above three systems (skip the Trp-cage system with no Z numbers or
  // explicit scaling factors (old form of Amber prmtop), as well as the Trp-cage system in water
  // (it's really a periodic box, useful for computing a meaningful 1:4 non-bonded energy sum).
  section(4);
  snapshot(trpi_snapshot, polyNumericVector(trpi_elec_frc), "trpcage_qq_forces",
           NumberFormat::SCIENTIFIC, "Forces due to electrostatic non-bonded (exclusive of 1:4) "
           "interactions in the Trp-cage (isolated boundary conditions) system do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::APPEND,
           snap_check);
  snapshot(trpi_snapshot, polyNumericVector(trpi_vdw_frc), "trpcage_lj_forces",
           NumberFormat::SCIENTIFIC, "Forces due to Lennard-Jones non-bonded (exclusive of 1:4) "
           "interactions in the Trp-cage (isolated boundary conditions) system do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-6, 1.0e-12, PrintSituation::APPEND,
           snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_elec_frc), "dhfr_qq_forces",
           NumberFormat::SCIENTIFIC, "Forces due to electrostatic non-bonded (exclusive of 1:4) "
           "interactions in the DHFR system do not meet expectations.", oe.takeSnapshot(), 1.0e-6,
           1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(dhfr_snapshot, polyNumericVector(dhfr_vdw_frc), "dhfr_lj_forces",
           NumberFormat::SCIENTIFIC, "Forces due to Lennard-Jones non-bonded (exclusive of 1:4) "
           "interactions in the DHFR system do not meet expectations.", oe.takeSnapshot(), 1.0e-6,
           1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(alad_snapshot, polyNumericVector(alad_elec_frc), "ala_dipeptide_qq_forces",
           NumberFormat::SCIENTIFIC, "Forces due to electrostatic non-bonded (exclusive of 1:4) "
           "interactions in the Ala dipeptide system do not meet expectations.", oe.takeSnapshot(),
           1.0e-6, 1.0e-12, PrintSituation::APPEND, snap_check);
  snapshot(alad_snapshot, polyNumericVector(alad_vdw_frc), "ala_dipeptide_lj_forces",
           NumberFormat::SCIENTIFIC, "Forces due to Lennard-Jones non-bonded (exclusive of 1:4) "
           "interactions in the Ala dipeptide system do not meet expectations.", oe.takeSnapshot(),
           1.0e-6, 1.0e-12, PrintSituation::APPEND, snap_check);

  // Check the relevant, accumulated energies against the double-precision standard
  section(5);
  const std::vector<double> trpi_acc = all_systems_sc.reportInstantaneousStates(trpi_idx);
  const std::vector<double> dhfr_acc = all_systems_sc.reportInstantaneousStates(dhfr_idx);
  const std::vector<double> alad_acc = all_systems_sc.reportInstantaneousStates(alad_idx);
  const std::vector<double> trpi_nbe_answer = {   1458.0998129,     62.6481811,
                                                 -1816.1478854,   -119.3500215 };
  const std::vector<double> dhfr_nbe_answer = {   6507.3376751,    367.0925927,
                                                -10036.4655324,  -1009.1558078 };
  const std::vector<double> alad_nbe_answer = {     46.8072425,      3.1589264,
                                                   -78.8463310,     -1.2452480 };
  const int qq14_idx = static_cast<size_t>(StateVariable::ELEC_ONE_FOUR);
  const int lj14_idx = static_cast<size_t>(StateVariable::VDW_ONE_FOUR);
  const int qqnb_idx = static_cast<size_t>(StateVariable::ELECTROSTATIC);
  const int ljnb_idx = static_cast<size_t>(StateVariable::VDW);
  const std::vector<double> trpi_acc_result = { trpi_acc[qq14_idx], trpi_acc[lj14_idx],
                                                trpi_acc[qqnb_idx], trpi_acc[ljnb_idx] };
  const std::vector<double> dhfr_acc_result = { dhfr_acc[qq14_idx], dhfr_acc[lj14_idx],
                                                dhfr_acc[qqnb_idx], dhfr_acc[ljnb_idx] };
  const std::vector<double> alad_acc_result = { alad_acc[qq14_idx], alad_acc[lj14_idx],
                                                alad_acc[qqnb_idx], alad_acc[ljnb_idx] };
  check(trpi_acc_result, RelationalOperator::EQUAL, Approx(trpi_nbe_answer).margin(2.0e-5),
        "Trp-cage (isolated boundary conditions) non-bonded energy accumulators do not reproduce "
        "the double-precision targets to within tolerances.", do_tests);
  check(dhfr_acc_result, RelationalOperator::EQUAL, Approx(dhfr_nbe_answer).margin(1.0e-3),
        "DHFR non-bonded energy accumulators do not reproduce the double-precision targets to "
        "within tolerances.", do_tests);
  check(alad_acc_result, RelationalOperator::EQUAL, Approx(alad_nbe_answer).margin(1.0e-5),
        "Alanine dipeptide non-bonded energy accumulators do not reproduce the double-precision "
        "targets to within tolerances.", do_tests);
  timer.assignTime(0);

  // Read a Trp-cage trajectory and evaluate the results in double, single, and fixed-precision
  // representations.
  const std::string trpcage_traj = base_crd_name + osc + "trpcage.crd";
  const bool traj_exists = (getDrivePathType(trpcage_traj) == DrivePathType::FILE &&
                            systems_exist);
  if (getDrivePathType(trpcage_traj) != DrivePathType::FILE) {
    rtWarn("A trajectory of Trp-cage conformations (isolated boundary conditions) was not found "
           "in " + trpcage_traj	+ ".  Check the ${STORMM_SOURCE} environment variable for "
           "validity.  Subsequent tests will be skipped.", "test_generalized_born");
  }
  const TestPriority do_traj_tests = (traj_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  const CoordinateSeries<double> trpcage_csd = (traj_exists) ?
                                               CoordinateSeries<double>(trpcage_traj,
                                                                        trpi_ag.getAtomCount()) :
                                               CoordinateSeries<double>();
  const CoordinateSeries<float> trpcage_csf(trpcage_csd);
  const CoordinateSeries<llint> trpcage_csi(trpcage_csd, 26);
  timer.assignTime(0);
  const int nframe = trpcage_csd.getFrameCount();
  std::vector<double2> traj_nbe_dd(nframe), traj_nbe_fd(nframe), traj_nbe_id(nframe);
  std::vector<double2> traj_nbe_df(nframe), traj_nbe_ff(nframe), traj_nbe_if(nframe);
  ScoreCard traj_sc(20, 1, 32);
  const NonbondedKit<double> trpi_nbk_d = trpi_ag.getDoublePrecisionNonbondedKit();
  for (int i = 0; i < nframe; i++) {
    traj_nbe_dd[i] = evaluateNonbondedEnergy<double, double>(trpi_nbk_d, trpi_semask.data(),
                                                             trpcage_csd.data(), &traj_sc, i);
    traj_nbe_fd[i] = evaluateNonbondedEnergy<float, double>(trpi_nbk_d, trpi_semask.data(),
                                                            trpcage_csf.data(), &traj_sc, i);
    traj_nbe_id[i] = evaluateNonbondedEnergy<llint, double>(trpi_nbk_d, trpi_semask.data(),
                                                            trpcage_csi.data(), &traj_sc, i);
  }
  const NonbondedKit<float> trpi_nbk_f = trpi_ag.getSinglePrecisionNonbondedKit();
  for (int i = 0; i < nframe; i++) {
    traj_nbe_df[i] = evaluateNonbondedEnergy<double, float>(trpi_nbk_f, trpi_semask.data(),
                                                            trpcage_csd.data(), &traj_sc, i);
    traj_nbe_ff[i] = evaluateNonbondedEnergy<float, float>(trpi_nbk_f, trpi_semask.data(),
                                                           trpcage_csf.data(), &traj_sc, i);
    traj_nbe_if[i] = evaluateNonbondedEnergy<llint, float>(trpi_nbk_f, trpi_semask.data(),
                                                          trpcage_csi.data(), &traj_sc, i);
  }
  double mue_elec_fd = 0.0;
  double mue_elec_id = 0.0;
  double mue_elec_df = 0.0;
  double mue_elec_ff = 0.0;
  double mue_elec_if = 0.0;
  double mue_vdw_fd = 0.0;
  double mue_vdw_id = 0.0;
  double mue_vdw_df = 0.0;
  double mue_vdw_ff = 0.0;
  double mue_vdw_if = 0.0;
  for (int i = 0; i < nframe; i++) {
    mue_elec_fd += fabs(traj_nbe_fd[i].x - traj_nbe_dd[i].x);
    mue_elec_id += fabs(traj_nbe_id[i].x - traj_nbe_dd[i].x);
    mue_elec_df += fabs(traj_nbe_df[i].x - traj_nbe_dd[i].x);
    mue_elec_ff += fabs(traj_nbe_ff[i].x - traj_nbe_dd[i].x);
    mue_elec_if += fabs(traj_nbe_if[i].x - traj_nbe_dd[i].x);
    mue_vdw_fd += fabs(traj_nbe_fd[i].y - traj_nbe_dd[i].y);
    mue_vdw_id += fabs(traj_nbe_id[i].y - traj_nbe_dd[i].y);
    mue_vdw_df += fabs(traj_nbe_df[i].y - traj_nbe_dd[i].y);
    mue_vdw_ff += fabs(traj_nbe_ff[i].y - traj_nbe_dd[i].y);
    mue_vdw_if += fabs(traj_nbe_if[i].y - traj_nbe_dd[i].y);
  }
  const double dnframe = nframe;
  mue_elec_fd /= dnframe;
  mue_elec_id /= dnframe;
  mue_elec_df /= dnframe;
  mue_elec_ff /= dnframe;
  mue_elec_if /= dnframe;
  mue_vdw_fd /= dnframe;
  mue_vdw_id /= dnframe;
  mue_vdw_df /= dnframe;
  mue_vdw_ff /= dnframe;
  mue_vdw_if /= dnframe;
  section(6);
  check(mue_elec_fd, RelationalOperator::LESS_THAN, 2.0e-5, "Mean unsigned error for "
        "electrostatic energies of Trp-cage frames is higher than expected when representing "
        "coordinates in single precision and calculating in double precision.", do_traj_tests);
  check(mue_vdw_fd, RelationalOperator::LESS_THAN, 5.0e-6, "Mean unsigned error for "
        "van-der Waals energies of Trp-cage frames is higher than expected when representing "
        "coordinates in single precision and calculating in double precision.", do_traj_tests);
  check(mue_elec_df, RelationalOperator::LESS_THAN, 1.0e-4, "Mean unsigned error for "
        "electrostatic energies of Trp-cage frames is higher than expected when representing "
        "coordinates in double precision and calculating in single precision.", do_traj_tests);
  check(mue_vdw_df, RelationalOperator::LESS_THAN, 1.0e-5, "Mean unsigned error for "
        "van-der Waals energies of Trp-cage frames is higher than expected when representing "
        "coordinates in double precision and calculating in single precision.", do_traj_tests);
  check(mue_elec_ff, RelationalOperator::LESS_THAN, 1.0e-4, "Mean unsigned error for "
        "electrostatic energies of Trp-cage frames is higher than expected when representing "
        "coordinates in single precision and calculating in single precision.", do_traj_tests);
  check(mue_vdw_ff, RelationalOperator::LESS_THAN, 1.0e-5, "Mean unsigned error for "
        "van-der Waals energies of Trp-cage frames is higher than expected when representing "
        "coordinates in single precision and calculating in single precision.", do_traj_tests);
  section(7);
  check(mue_elec_id, RelationalOperator::LESS_THAN, 2.0e-6, "Mean unsigned error for "
        "electrostatic energies of Trp-cage frames is higher than expected when representing "
        "coordinates in fixed precision and calculating in double precision.", do_traj_tests);
  check(mue_vdw_id, RelationalOperator::LESS_THAN, 1.0e-6, "Mean unsigned error for "
        "van-der Waals energies of Trp-cage frames is higher than expected when representing "
        "coordinates in fixed precision and calculating in double precision.", do_traj_tests);
  check(mue_elec_if, RelationalOperator::LESS_THAN, 1.0e-4, "Mean unsigned error for "
        "electrostatic energies of Trp-cage frames is higher than expected when representing "
        "coordinates in fixed precision and calculating in single precision.", do_traj_tests);
  check(mue_vdw_if, RelationalOperator::LESS_THAN, 1.0e-5, "Mean unsigned error for "
        "van-der Waals energies of Trp-cage frames is higher than expected when representing "
        "coordinates in fixed precision and calculating in single precision.", do_traj_tests);
  section(6);
  testNBPrecisionModel<double, float, double>(trpi_nbk_d, trpi_semask, trpi_ps, trpi_elec_frc,
                                              trpi_vdw_frc, 5.0e-6, do_tests);
  testNBPrecisionModel<double, float, float>(trpi_nbk_f, trpi_semask, trpi_ps, trpi_elec_frc,
                                             trpi_vdw_frc, 5.0e-5, do_tests);
  testNBPrecisionModel<double, llint, double>(trpi_nbk_d, trpi_semask, trpi_ps, trpi_elec_frc,
                                              trpi_vdw_frc, 7.0e-7, do_tests);
  testNBPrecisionModel<double, llint, float>(trpi_nbk_f, trpi_semask, trpi_ps, trpi_elec_frc,
                                             trpi_vdw_frc, 5.0e-5, do_tests);
  testNBPrecisionModel<float, double, double>(trpi_nbk_d, trpi_semask, trpi_ps, trpi_elec_frc,
                                              trpi_vdw_frc, 5.0e-5, do_tests);
  testNBPrecisionModel<float, double, float>(trpi_nbk_f, trpi_semask, trpi_ps, trpi_elec_frc,
                                             trpi_vdw_frc, 5.0e-5, do_tests);
  testNBPrecisionModel<float, float, double>(trpi_nbk_d, trpi_semask, trpi_ps, trpi_elec_frc,
                                             trpi_vdw_frc, 5.0e-5, do_tests);
  testNBPrecisionModel<float, float, float>(trpi_nbk_f, trpi_semask, trpi_ps, trpi_elec_frc,
                                            trpi_vdw_frc, 5.0e-5, do_tests);
  testNBPrecisionModel<float, llint, double>(trpi_nbk_d, trpi_semask, trpi_ps, trpi_elec_frc,
                                             trpi_vdw_frc, 5.0e-5, do_tests);
  testNBPrecisionModel<float, llint, float>(trpi_nbk_f, trpi_semask, trpi_ps, trpi_elec_frc,
                                            trpi_vdw_frc, 5.0e-5, do_tests);
  testNBPrecisionModel<llint, double, double>(trpi_nbk_d, trpi_semask, trpi_ps, trpi_elec_frc,
                                              trpi_vdw_frc, 1.8e-6, do_tests);
  testNBPrecisionModel<llint, double, float>(trpi_nbk_f, trpi_semask, trpi_ps, trpi_elec_frc,
                                             trpi_vdw_frc, 5.0e-5, do_tests);
  testNBPrecisionModel<llint, float, double>(trpi_nbk_d, trpi_semask, trpi_ps, trpi_elec_frc,
                                             trpi_vdw_frc, 5.0e-6, do_tests);
  testNBPrecisionModel<llint, float, float>(trpi_nbk_f, trpi_semask, trpi_ps, trpi_elec_frc,
                                            trpi_vdw_frc, 5.0e-5, do_tests);
  testNBPrecisionModel<llint, llint, double>(trpi_nbk_d, trpi_semask, trpi_ps, trpi_elec_frc,
                                             trpi_vdw_frc, 5.0e-6, do_tests);
  testNBPrecisionModel<llint, llint, float>(trpi_nbk_f, trpi_semask, trpi_ps, trpi_elec_frc,
                                            trpi_vdw_frc, 5.0e-5, do_tests);

  // Check the energy tracking object (placed in this test program for convenience)
  section(8);
  Xoshiro256ppGenerator xrs;
  int test_no = 0;
  const int nvar = static_cast<int>(StateVariable::ALL_STATES);
  const int nsys = 17;
  std::vector<int> tsc_sample_sizes(10);
  std::vector<int> tsc_sample_size_ans(10);
  std::vector<double> tsc_ave_dev(10 * nsys * nvar), tsc_ave_dev_ans(10 * nsys * nvar);
  std::vector<double> tsc_std_dev(10 * nsys * nvar), tsc_std_dev_ans(10 * nsys * nvar);
  bool lateral_survey = true;
  bool sc_importing_tested = false;
  for (int npts = 3; npts <= 48; npts += 5) {
    ScoreCard trial_sc(nsys, npts, 32);
    ScoreCardWriter trial_scw = trial_sc.data();
    std::vector<std::vector<double>> fake_bond_e(nsys, std::vector<double>(npts));
    std::vector<std::vector<double>> fake_angl_e(nsys, std::vector<double>(npts));
    std::vector<std::vector<double>> fake_dihe_e(nsys, std::vector<double>(npts));
    std::vector<std::vector<double>> fake_impr_e(nsys, std::vector<double>(npts));
    std::vector<std::vector<double>> fake_ubrd_e(nsys, std::vector<double>(npts));
    std::vector<std::vector<double>> fake_cimp_e(nsys, std::vector<double>(npts));
    std::vector<std::vector<double>> fake_cmap_e(nsys, std::vector<double>(npts));
    std::vector<std::vector<double>> fake_qq14_e(nsys, std::vector<double>(npts));
    std::vector<std::vector<double>> fake_lj14_e(nsys, std::vector<double>(npts));
    std::vector<std::vector<double>> fake_qqnb_e(nsys, std::vector<double>(npts));
    std::vector<std::vector<double>> fake_ljnb_e(nsys, std::vector<double>(npts));
    std::vector<std::vector<double>> fake_gbrn_e(nsys, std::vector<double>(npts));
    std::vector<std::vector<double>> fake_rstr_e(nsys, std::vector<double>(npts));
    std::vector<std::vector<double>> fake_kine_e(nsys, std::vector<double>(npts));
    for (int j = 0; j < npts; j++) {
      for (int i = 0; i < nsys; i++) {
        fake_bond_e[i][j] =   20.0 * xrs.uniformRandomNumber();
        fake_angl_e[i][j] =   40.0 * xrs.uniformRandomNumber();
        fake_dihe_e[i][j] =  125.0 * (0.50 - xrs.uniformRandomNumber());
        fake_impr_e[i][j] =   25.0 * xrs.uniformRandomNumber();
        fake_ubrd_e[i][j] =   35.0 * xrs.uniformRandomNumber();
        fake_cimp_e[i][j] =   25.0 * xrs.uniformRandomNumber();
        fake_cmap_e[i][j] =   15.0 * (0.50 - xrs.uniformRandomNumber());
        fake_qq14_e[i][j] =  109.5 * (0.25 - xrs.uniformRandomNumber());
        fake_lj14_e[i][j] =   75.5 * (0.05 - xrs.uniformRandomNumber());
        fake_qqnb_e[i][j] =  250.0 * (0.50 - xrs.uniformRandomNumber());
        fake_ljnb_e[i][j] =  150.0 * (0.05 - xrs.uniformRandomNumber());
        fake_gbrn_e[i][j] = -(0.90 + (0.10 * xrs.uniformRandomNumber())) * fake_qqnb_e[i][j];
        fake_rstr_e[i][j] =   15.0 * xrs.uniformRandomNumber();
        fake_kine_e[i][j] =  255.0 * (0.5 + xrs.uniformRandomNumber());
        const llint ll_bond_e = static_cast<llint>(fake_bond_e[i][j] * trial_scw.nrg_scale_lf);
        const llint ll_angl_e = static_cast<llint>(fake_angl_e[i][j] * trial_scw.nrg_scale_lf);
        const llint ll_dihe_e = static_cast<llint>(fake_dihe_e[i][j] * trial_scw.nrg_scale_lf);
        const llint ll_impr_e = static_cast<llint>(fake_impr_e[i][j] * trial_scw.nrg_scale_lf);
        const llint ll_ubrd_e = static_cast<llint>(fake_ubrd_e[i][j] * trial_scw.nrg_scale_lf);
        const llint ll_cimp_e = static_cast<llint>(fake_cimp_e[i][j] * trial_scw.nrg_scale_lf);
        const llint ll_cmap_e = static_cast<llint>(fake_cmap_e[i][j] * trial_scw.nrg_scale_lf);
        const llint ll_qq14_e = static_cast<llint>(fake_qq14_e[i][j] * trial_scw.nrg_scale_lf);
        const llint ll_lj14_e = static_cast<llint>(fake_lj14_e[i][j] * trial_scw.nrg_scale_lf);
        const llint ll_qqnb_e = static_cast<llint>(fake_qqnb_e[i][j] * trial_scw.nrg_scale_lf);
        const llint ll_ljnb_e = static_cast<llint>(fake_ljnb_e[i][j] * trial_scw.nrg_scale_lf);
        const llint ll_gbrn_e = static_cast<llint>(fake_gbrn_e[i][j] * trial_scw.nrg_scale_lf);
        const llint ll_rstr_e = static_cast<llint>(fake_rstr_e[i][j] * trial_scw.nrg_scale_lf);
        const llint ll_kine_e = static_cast<llint>(fake_kine_e[i][j] * trial_scw.nrg_scale_lf);
        trial_sc.contribute(StateVariable::BOND, ll_bond_e, i);
        trial_sc.contribute(StateVariable::ANGLE, ll_angl_e, i);
        trial_sc.contribute(StateVariable::PROPER_DIHEDRAL, ll_dihe_e, i);
        trial_sc.contribute(StateVariable::IMPROPER_DIHEDRAL, ll_impr_e, i);
        trial_sc.contribute(StateVariable::UREY_BRADLEY, ll_ubrd_e, i);
        trial_sc.contribute(StateVariable::CHARMM_IMPROPER, ll_cimp_e, i);
        trial_sc.contribute(StateVariable::CMAP, ll_cmap_e, i);
        trial_sc.contribute(StateVariable::ELEC_ONE_FOUR, ll_qq14_e, i);
        trial_sc.contribute(StateVariable::VDW_ONE_FOUR, ll_lj14_e, i);
        trial_sc.contribute(StateVariable::ELECTROSTATIC, ll_qqnb_e, i);
        trial_sc.contribute(StateVariable::VDW, ll_ljnb_e, i);
        trial_sc.contribute(StateVariable::GENERALIZED_BORN, ll_gbrn_e, i);
        trial_sc.contribute(StateVariable::RESTRAINT, ll_rstr_e, i);
        trial_sc.contribute(StateVariable::KINETIC, ll_kine_e, i);
      }
      trial_sc.incrementSampleCount();
    }
    tsc_sample_sizes[test_no] = trial_sc.getSampleSize();
    tsc_sample_size_ans[test_no] = npts;
    for (int i = 0; i < nvar; i++) {
      const StateVariable istate = static_cast<StateVariable>(i);

      // Loop over the systems to get longitudinal measurements of mean and standard deviations
      // for each potential energy quantity.
      for (int j = 0; j < nsys; j++) {
        const int tsc_idx = (test_no * nvar * nsys) + (i * nsys) + j;
        tsc_ave_dev[tsc_idx] = trial_sc.reportAverageStates(istate, j);
        tsc_std_dev[tsc_idx] = trial_sc.reportVarianceOfStates(istate, j);
        switch (istate) {
        case StateVariable::BOND:
          tsc_ave_dev_ans[tsc_idx] = mean(fake_bond_e[j]);
          tsc_std_dev_ans[tsc_idx] = variance(fake_bond_e[j], VarianceMethod::STANDARD_DEVIATION);
          break;
        case StateVariable::ANGLE:
          tsc_ave_dev_ans[tsc_idx] = mean(fake_angl_e[j]);
          tsc_std_dev_ans[tsc_idx] = variance(fake_angl_e[j], VarianceMethod::STANDARD_DEVIATION);
          break;
        case StateVariable::PROPER_DIHEDRAL:
          tsc_ave_dev_ans[tsc_idx] = mean(fake_dihe_e[j]);
          tsc_std_dev_ans[tsc_idx] = variance(fake_dihe_e[j], VarianceMethod::STANDARD_DEVIATION);
          break;
        case StateVariable::IMPROPER_DIHEDRAL:
          tsc_ave_dev_ans[tsc_idx] = mean(fake_impr_e[j]);
          tsc_std_dev_ans[tsc_idx] = variance(fake_impr_e[j], VarianceMethod::STANDARD_DEVIATION);
          break;
        case StateVariable::UREY_BRADLEY:
          tsc_ave_dev_ans[tsc_idx] = mean(fake_ubrd_e[j]);
          tsc_std_dev_ans[tsc_idx] = variance(fake_ubrd_e[j], VarianceMethod::STANDARD_DEVIATION);
          break;
        case StateVariable::CHARMM_IMPROPER:
          tsc_ave_dev_ans[tsc_idx] = mean(fake_cimp_e[j]);
          tsc_std_dev_ans[tsc_idx] = variance(fake_cimp_e[j], VarianceMethod::STANDARD_DEVIATION);
          break;
        case StateVariable::CMAP:
          tsc_ave_dev_ans[tsc_idx] = mean(fake_cmap_e[j]);
          tsc_std_dev_ans[tsc_idx] = variance(fake_cmap_e[j], VarianceMethod::STANDARD_DEVIATION);
          break;
        case StateVariable::ELECTROSTATIC:
          tsc_ave_dev_ans[tsc_idx] = mean(fake_qqnb_e[j]);
          tsc_std_dev_ans[tsc_idx] = variance(fake_qqnb_e[j], VarianceMethod::STANDARD_DEVIATION);
          break;
        case StateVariable::VDW:
          tsc_ave_dev_ans[tsc_idx] = mean(fake_ljnb_e[j]);
          tsc_std_dev_ans[tsc_idx] = variance(fake_ljnb_e[j], VarianceMethod::STANDARD_DEVIATION);
          break;
        case StateVariable::ELEC_ONE_FOUR:
          tsc_ave_dev_ans[tsc_idx] = mean(fake_qq14_e[j]);
          tsc_std_dev_ans[tsc_idx] = variance(fake_qq14_e[j], VarianceMethod::STANDARD_DEVIATION);
          break;
        case StateVariable::VDW_ONE_FOUR:
          tsc_ave_dev_ans[tsc_idx] = mean(fake_lj14_e[j]);
          tsc_std_dev_ans[tsc_idx] = variance(fake_lj14_e[j], VarianceMethod::STANDARD_DEVIATION);
          break;
        case StateVariable::GENERALIZED_BORN:
          tsc_ave_dev_ans[tsc_idx] = mean(fake_gbrn_e[j]);
          tsc_std_dev_ans[tsc_idx] = variance(fake_gbrn_e[j], VarianceMethod::STANDARD_DEVIATION);
          break;
        case StateVariable::RESTRAINT:
          tsc_ave_dev_ans[tsc_idx] = mean(fake_rstr_e[j]);
          tsc_std_dev_ans[tsc_idx] = variance(fake_rstr_e[j], VarianceMethod::STANDARD_DEVIATION);
          break;
        case StateVariable::KINETIC:
          tsc_ave_dev_ans[tsc_idx] = mean(fake_kine_e[j]);
          tsc_std_dev_ans[tsc_idx] = variance(fake_kine_e[j], VarianceMethod::STANDARD_DEVIATION);
          break;
        case StateVariable::PRESSURE:
        case StateVariable::VIRIAL_11:
        case StateVariable::VIRIAL_12:
        case StateVariable::VIRIAL_22:
        case StateVariable::VIRIAL_13:
        case StateVariable::VIRIAL_23:
        case StateVariable::VIRIAL_33:
        case StateVariable::VOLUME:
        case StateVariable::TEMPERATURE_ALL:
        case StateVariable::TEMPERATURE_PROTEIN:
        case StateVariable::TEMPERATURE_LIGAND:
        case StateVariable::TEMPERATURE_SOLVENT:
        case StateVariable::DU_DLAMBDA:
        case StateVariable::POTENTIAL_ENERGY:
        case StateVariable::TOTAL_ENERGY:
        case StateVariable::ALL_STATES:
          tsc_ave_dev_ans[tsc_idx] = 0.0;
          tsc_std_dev_ans[tsc_idx] = 0.0;
          break;
        }
      }

      // Take lateral measurements over all "systems" at each time point, recording the mean and
      // standard deviations.
      if (i < 14) {
        const std::vector<double2> ehist = trial_sc.reportHistory(istate);
        std::vector<double2> ehist_chk(npts);
        for (int j = 0; j < npts; j++) {
          std::vector<double> survey(nsys);
          for (int k = 0; k < nsys; k++) {
            switch (istate) {
            case StateVariable::BOND:
              survey[k] = fake_bond_e[k][j];
              break;
            case StateVariable::ANGLE:
              survey[k] = fake_angl_e[k][j];
              break;
            case StateVariable::PROPER_DIHEDRAL:
              survey[k] = fake_dihe_e[k][j];
              break;
            case StateVariable::IMPROPER_DIHEDRAL:
              survey[k] = fake_impr_e[k][j];
              break;
            case StateVariable::UREY_BRADLEY:
              survey[k] = fake_ubrd_e[k][j];
              break;
            case StateVariable::CHARMM_IMPROPER:
              survey[k] = fake_cimp_e[k][j];
              break;
            case StateVariable::CMAP:
              survey[k] = fake_cmap_e[k][j];
              break;
            case StateVariable::ELECTROSTATIC:
              survey[k] = fake_qqnb_e[k][j];
              break;
            case StateVariable::VDW:
              survey[k] = fake_ljnb_e[k][j];
              break;
            case StateVariable::ELEC_ONE_FOUR:
              survey[k] = fake_qq14_e[k][j];
              break;
            case StateVariable::VDW_ONE_FOUR:
              survey[k] = fake_lj14_e[k][j];
              break;
            case StateVariable::GENERALIZED_BORN:
              survey[k] = fake_gbrn_e[k][j];
              break;
            case StateVariable::RESTRAINT:
              survey[k] = fake_rstr_e[k][j];
              break;
            case StateVariable::KINETIC:
              survey[k] = fake_kine_e[k][j];
              break;
            case StateVariable::PRESSURE:
            case StateVariable::VIRIAL_11:
            case StateVariable::VIRIAL_12:
            case StateVariable::VIRIAL_22:
            case StateVariable::VIRIAL_13:
            case StateVariable::VIRIAL_23:
            case StateVariable::VIRIAL_33:
            case StateVariable::VOLUME:
            case StateVariable::TEMPERATURE_ALL:
            case StateVariable::TEMPERATURE_PROTEIN:
            case StateVariable::TEMPERATURE_LIGAND:
            case StateVariable::TEMPERATURE_SOLVENT:
            case StateVariable::DU_DLAMBDA:
            case StateVariable::POTENTIAL_ENERGY:
            case StateVariable::TOTAL_ENERGY:
            case StateVariable::ALL_STATES:
              break;
            }
          }
          ehist_chk[j].x = mean(survey);
          ehist_chk[j].y = variance(survey, VarianceMethod::STANDARD_DEVIATION);
        }
        for (int j = 0; j < npts; j++) {
          lateral_survey = (lateral_survey && fabs(ehist[j].x - ehist_chk[j].x) < 1.0e-6 &&
                            fabs(ehist[j].y - ehist_chk[j].y) < 1.0e-6);
        }
      }
    }
    test_no++;

    // Try assembling the pieces of the original ScoreCard into a large one
    if (sc_importing_tested == false) {
      sc_importing_tested = true;
      ScoreCard multi_sc(4 * nsys, 5, 24);
      const std::vector<double> trial_tote = trial_sc.reportTotalEnergies();
      std::vector<double> tote_target(4 * nsys);
      for (int i = 0; i < nsys; i++) {
        multi_sc.import(trial_sc, i, i);
        multi_sc.import(trial_sc, (2 * nsys) - (i + 1), i);
        multi_sc.import(trial_sc, (2 * nsys) + i, i);
        multi_sc.import(trial_sc, (3 * nsys) + i, i);
        tote_target[i] = trial_tote[i];
        tote_target[(2 * nsys) - (i + 1)] = trial_tote[i];
        tote_target[(2 * nsys) + i] = trial_tote[i];
        tote_target[(3 * nsys) + i] = trial_tote[i];
      }
      CHECK_THROWS(multi_sc.import(trial_sc, (4 * nsys) + 8, 0), "Importation of information into "
                   "a non-existent entry of the destination was allowed.");
      CHECK_THROWS(multi_sc.import(trial_sc, 0, (3 * nsys) + 2), "Importation of information from "
                   "a non-existent entry of the source material was allowed.");
      const std::vector<double> multi_tote = multi_sc.reportTotalEnergies();
      check(multi_tote, RelationalOperator::EQUAL, tote_target, "Energies imported into a new "
            "ScoreCard do not align with a similar series accumulated from the original source.");
    }
  }
  check(tsc_sample_sizes, RelationalOperator::EQUAL, tsc_sample_size_ans, "Sample sizes for "
        "various time series are not reported correctly by the ScoreCard object.");
  check(tsc_ave_dev, RelationalOperator::EQUAL, tsc_ave_dev_ans, "Mean values for "
        "various time series are not reported correctly by the ScoreCard object.");
  check(tsc_std_dev, RelationalOperator::EQUAL, tsc_std_dev_ans, "Standard deviations for "
        "various time series are not reported correctly by the ScoreCard object.");
  check(lateral_survey, "A survey of the lateral mean values and standard deviations among "
        "familiar potential energy terms, taken point by point throughout the history, failed to "
        "confirm the results reported by the tracking object.");
  
  // Check the non-bonded softcore potentials by iteratively feeding them different values of
  // the displacement and checking the results against analytic numbers.
  section(9);
  for (int order = 1; order <= 3; order++) {
    testSoftCoreCoulomb( 0.54 *  0.31, 0.9, order);
    testSoftCoreCoulomb( 0.47 * -0.24, 0.9, order);
    testSoftCoreCoulomb(-0.17 *  0.84, 0.5, order);
  }
  for (int order = 3; order <= 4; order++) {
    testSoftCoreLennardJones( 79716.0, 109.35, 0.6, order);
    testSoftCoreLennardJones(110272.8, 124.25, 0.5, order);
  }
  testSplineSoftCore(trpi_ag, 1.0, 0.0, 8.0, 50.0, 206.0, do_tests);
  
  // Compute the potential and forces on some configurations that have clashes.
  const std::string base_nml_topl = oe.getStormmSourcePath() + osc + "test" + osc + "Namelists" +
                                    osc + "topol";
  const std::string base_nml_icrd = oe.getStormmSourcePath() + osc + "test" + osc + "Namelists" +
                                    osc + "coord";
  TestSystemManager tsm(base_nml_topl, "top", { "gly_phe", "gly_trp", "gly_lys" }, base_nml_icrd,
                        "inpcrd", { "gly_phe", "gly_trp", "gly_lys" });
  std::vector<PhaseSpace> clashing_structures;
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    const AtomGraph* iag_ptr = tsm.getTopologyPointer(i);
    const PhaseSpace ips = tsm.exportPhaseSpace(i);
    const PhaseSpaceReader ipsr = ips.data();
    CoordinateFrame icf(ips);
    CoordinateFrameWriter icfw = icf.data();
    const ChemicalFeatures ichemfe(iag_ptr, icf, MapRotatableGroups::YES);
    const std::vector<IsomerPlan> irgroups = ichemfe.getRotatableBondGroups();
    const int nrot = ichemfe.getRotatableBondCount();
    for (int j = 0; j < nrot; j++) {
      coordCopy(&icfw, ipsr, TrajectoryKind::POSITIONS);
      rotateAboutBond(icfw.xcrd, icfw.ycrd, icfw.zcrd, irgroups[j].getRootAtom(),
                      irgroups[j].getPivotAtom(), irgroups[j].getMovingAtoms().data(),
                      irgroups[j].getMovingAtoms().size(),
                      xrs.uniformRandomNumber() * stormm::symbols::twopi);
    }
  }

  // Check the non-bonded PME direct space interaction and its derivatives
  section(10);
  runPMEDirectSpaceTests();

  // Check the potential decomposition
  section(11);
  runLayeredCoulombTests(oe, DecomposablePotential::ELECTROSTATIC, 5.4, PrintSituation::OVERWRITE);
  runLayeredCoulombTests(oe, DecomposablePotential::DISPERSION, 5.4);
  runLayeredCoulombTests(oe, DecomposablePotential::ELEC_PME_DIRECT);
  
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
