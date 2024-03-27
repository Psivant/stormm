#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/mesh_kernel_manager.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Math/math_enumerators.h"
#include "../../src/Math/statistical_enumerators.h"
#include "../../src/Math/tricubic_cell.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/scorecard.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/present_field.h"
#include "../../src/Reporting/reporting_enumerators.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/background_mesh.h"
#include "../../src/Structure/mesh_mechanics.h"
#include "../../src/Structure/mesh_parameters.h"
#include "../../src/Structure/structure_enumerators.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::constants::getEnumerationName;
using stormm::constants::PrecisionModel;
using stormm::data_types::ullint;
#ifndef STORMM_USE_HPC
using stormm::data_types::double3;
using stormm::data_types::double4;
#endif
using stormm::energy::ScoreCard;
using stormm::parse::minimalRealFormat;
using stormm::review::GridFileSyntax;
using stormm::review::printToFile;
using stormm::stmath::computeBoxTransform;
using stormm::stmath::extractBoxDimensions;
using stormm::stmath::getEnumerationName;
using stormm::stmath::maxAbsValue;
using stormm::stmath::maxValue;
using stormm::stmath::minValue;
using stormm::stmath::mean;
using stormm::stmath::readBitFromMask;
using stormm::stmath::sum;
using stormm::stmath::TricubicCell;
using stormm::stmath::TricubicStencil;
using stormm::stmath::variance;
using stormm::stmath::VarianceMethod;
using stormm::numerics::hostInt95ToDouble;
using stormm::random::Xoshiro256ppGenerator;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using namespace stormm::card;
using namespace stormm::diskutil;
using namespace stormm::structure;
using namespace stormm::testing;
using namespace stormm::topology;

// A tolerance on the fraction of failed points
constexpr int ok_point_ratio = 10000;

//-------------------------------------------------------------------------------------------------
// Produce a message to describe the geometry of a mesh element in the event of a test failure.
//
// Arguments:
//   invu:          Transformation matrix taking fractional coordinates on the mesh into Cartesian
//                  space
//   coeff_type:    Type of the coefficients stored on the mesh
//   build_type:    Floating-point data type used in mesh calculations
//   stencil_type:  The stencil used to construct the interpolant
//-------------------------------------------------------------------------------------------------
template <typename T> std::string meshDescriptor(const T* invu, const std::string &coeff_type,
                                                 const std::string &build_type,
                                                 const std::string &stencil_type) {
  T box_a, box_b, box_c, box_alpha, box_beta, box_gamma;
  const T degree_conv = 180.0 / stormm::symbols::pi;
  extractBoxDimensions(&box_a, &box_b, &box_c, &box_alpha, &box_beta, &box_gamma, invu);
  const NumberFormat nmb_fmt = NumberFormat::STANDARD_REAL;
  const std::string mesh_err("Mesh element dimensions: [ " + realToString(box_a, 9, 4, nmb_fmt) +
                             ", " + realToString(box_b, 9, 4, nmb_fmt) + ", " +
                             realToString(box_c, 9, 4, nmb_fmt) + ", " +
                             realToString(box_alpha * degree_conv, 9, 4, nmb_fmt) + ", " +
                             realToString(box_beta  * degree_conv, 9, 4, nmb_fmt) + ", " +
                             realToString(box_gamma * degree_conv, 9, 4, nmb_fmt) +
                             " ].  Coefficient data type: " + coeff_type + ".  Build precision: " +
                             build_type + ".  Stencil type: " + stencil_type + ".");
  return mesh_err;
}

//-------------------------------------------------------------------------------------------------
// Create a TricubicCell object with a specified origin and dimensons based on the electrostatic
// potential due to a molecule's rigid atoms.  This routine will return a tuple containing the
// Cartesian X, Y, and Z components of the electrostatic force in the "x", "y", and "z" members,
// plus the electrostatic energy in the "w" member.
//
// Arguments:
//   ag:        Topology of the molecule, including atom mobilities and atomic partial charges
//   cf:        Coordinates of the molecule
//   origin:    Cartesian coordinates for the chosen origin of the TricubicCell object
//   bounds:    3 x 3 matrix specifying the inverse transformation matrix for the TricubicCell
//   test_pt:   Point at which to test the electrostatic potential and force due to the molecule
//-------------------------------------------------------------------------------------------------
double4 checkTricubicElectrostatics(const AtomGraph &ag, const CoordinateFrame &cf,
                                    const double3 origin, const std::vector<double> &bounds,
                                    const double3 test_pt) {

  // Determine the coordinates for the box boundaries
  std::vector<double> xbnd(8), ybnd(8), zbnd(8);
  for (int i = 0; i < 2; i++) {
    const double di = i;
    for (int j = 0; j < 2; j++) {
      const double dj = j;
      for (int k = 0; k < 2; k++) {
        const int pos = (((k * 2) + j) * 2) + i;
        const double dk = k;
        xbnd[pos] = (bounds[0] * di) + (bounds[3] * dj) + (bounds[6] * dk) + origin.x;
        ybnd[pos] =                    (bounds[4] * dj) + (bounds[7] * dk) + origin.y;
        zbnd[pos] =                                       (bounds[8] * dk) + origin.z;
      }
    }
  }
  
  // Accumulate the potentials at each boundary point
  std::vector<double> u(8, 0.0), du_dx(8, 0.0), du_dy(8, 0.0), du_dz(8, 0.0), du_dxx(8, 0.0);
  std::vector<double> du_dxy(8, 0.0), du_dxz(8, 0.0), du_dyy(8, 0.0), du_dyz(8, 0.0);
  std::vector<double> du_dxxx(8, 0.0), du_dxxy(8, 0.0), du_dxxz(8, 0.0), du_dxyy(8, 0.0);
  std::vector<double> du_dxyz(8, 0.0);
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  const CoordinateFrameReader cfr = cf.data();
  const std::vector<bool> mobile_atoms = ag.getAtomMobility();
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      for (int k = 0; k < 2; k++) {
        const int pt_idx = (((k * 2) + j) * 2) + i;
        const double xloc = xbnd[pt_idx];
        const double yloc = ybnd[pt_idx];
        const double zloc = zbnd[pt_idx];
        for (int m = 0; m < nbk.natom; m++) {
          if (mobile_atoms[m]) {
            continue;
          }
          const double dx = xloc - cfr.xcrd[m];
          const double dy = yloc - cfr.ycrd[m];
          const double dz = zloc - cfr.zcrd[m];
          const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
          const double r = sqrt(r2);
          const double invr = 1.0 / r;
          const double invr2 = invr * invr;
          const double scaled_q = nbk.charge[m] * nbk.coulomb_constant;
          const double4 qq = {        scaled_q * invr,               -scaled_q * invr2,
                                2.0 * scaled_q * invr2 * invr, -6.0 * scaled_q * invr2 * invr2 };
          u[pt_idx] += qq.x;
          du_dx[pt_idx] += radialFirstDerivative<double>(qq.y, dx, r);
          du_dy[pt_idx] += radialFirstDerivative<double>(qq.y, dy, r);
          du_dz[pt_idx] += radialFirstDerivative<double>(qq.y, dz, r);
          du_dxx[pt_idx] += radialSecondDerivative<double>(qq.y, qq.z, dx, r);
          du_dxy[pt_idx] += radialSecondDerivative<double>(qq.y, qq.z, dx, dy, r, r2);
          du_dxz[pt_idx] += radialSecondDerivative<double>(qq.y, qq.z, dx, dz, r, r2);
          du_dyy[pt_idx] += radialSecondDerivative<double>(qq.y, qq.z, dy, r);
          du_dyz[pt_idx] += radialSecondDerivative<double>(qq.y, qq.z, dy, dz, r, r2);
          du_dxxx[pt_idx] += radialThirdDerivative<double>(qq.y, qq.z, qq.w, dx, r, r2);
          du_dxxy[pt_idx] += radialThirdDerivative<double>(qq.y, qq.z, qq.w, dx, dy, r, r2);
          du_dxxz[pt_idx] += radialThirdDerivative<double>(qq.y, qq.z, qq.w, dx, dz, r, r2);
          du_dxyy[pt_idx] += radialThirdDerivative<double>(qq.y, qq.z, qq.w, dy, dx, r, r2);
          du_dxyz[pt_idx] += radialThirdDerivative<double>(qq.y, qq.z, qq.w, dx, dy, dz, r, r2);
        }
      }
    }
  }

  // Form the TricubicCell object
  const TricubicStencil tc_weights;
  const std::vector<double> tc_bounds = { origin.x, origin.y, origin.z, bounds[0], bounds[1],
                                          bounds[2], bounds[3], bounds[4], bounds[5], bounds[6],
                                          bounds[7], bounds[8] };
  TricubicCell<double> tcc(tc_weights, tc_bounds, u, du_dx, du_dy, du_dz, du_dxx, du_dxy, du_dxz,
                           du_dyy, du_dyz, du_dxxx, du_dxxy, du_dxxz, du_dxyy, du_dxyz);
  const double intrp_e = tcc.evaluate(test_pt.x, test_pt.y, test_pt.z);
  const double3 intrp_f = tcc.derivative<double3>(test_pt.x, test_pt.y, test_pt.z);
  return { intrp_f.x, intrp_f.y, intrp_f.z, intrp_e };
}

//-------------------------------------------------------------------------------------------------
// Check the non-bonded potentials and forces at the exact mesh points against an independent
// calculation.
//
// Arguments:
//   ele_mesh:  Mesh containing the electrostatic potential due to the molecule
//   vdw_mesh:  Mesh containing the Lennard-Jones potential due to the molecule (must reference
//              the same molecule in elec_mesh)
//   xrs:       Random number generator, used to create scattered points
//   tsm:       Collection of test systems, tracking whether tests are feasible
//-------------------------------------------------------------------------------------------------
template <typename T>
void checkNonbondedPotentials(const BackgroundMesh<T> &ele_mesh, const BackgroundMesh<T> &vdw_mesh,
                              Xoshiro256ppGenerator *xrs, const TestSystemManager &tsm) {

  // Check compatibility
  const BackgroundMeshReader ele_r = ele_mesh.data();
  const BackgroundMeshReader vdw_r = vdw_mesh.data();
  const MeshFoundation& ele_bss = ele_mesh.getMolecularBasis();
  const MeshFoundation& vdw_bss = vdw_mesh.getMolecularBasis();
  if (ele_bss.getCoordinatePointer() != vdw_bss.getCoordinatePointer() ||
      ele_bss.getTopologyPointer() != vdw_bss.getTopologyPointer() ||
      ele_r.dims.na != vdw_r.dims.na || ele_r.dims.nb != vdw_r.dims.nb ||
      ele_r.dims.nc != vdw_r.dims.nc) {
    rtErr("Meshes must reference the same coordinates and topologies, and have the same sizes.\n",
          "checkNonbondedPotentials");
  }
  
  // Compute the coordinates of each mesh point
  const int na_plus_one = ele_r.dims.na + 1;
  const int nb_plus_one = ele_r.dims.nb + 1;
  const int nc_plus_one = ele_r.dims.nc + 1;
  std::vector<double> gpt_ax(na_plus_one), gpt_ay(na_plus_one), gpt_az(na_plus_one);
  std::vector<double> gpt_bx(nb_plus_one), gpt_by(nb_plus_one), gpt_bz(nb_plus_one);
  std::vector<double> gpt_cx(nc_plus_one), gpt_cy(nc_plus_one), gpt_cz(nc_plus_one);
  hostInt95ToDouble(gpt_ax.data(), ele_r.rlrs.avec_abs_x, ele_r.rlrs.avec_abs_x_ovrf, na_plus_one,
                    ele_r.dims.inv_scale);
  hostInt95ToDouble(gpt_ay.data(), ele_r.rlrs.avec_abs_y, ele_r.rlrs.avec_abs_y_ovrf, na_plus_one,
                    ele_r.dims.inv_scale);
  hostInt95ToDouble(gpt_az.data(), ele_r.rlrs.avec_abs_z, ele_r.rlrs.avec_abs_z_ovrf, na_plus_one,
                    ele_r.dims.inv_scale);
  hostInt95ToDouble(gpt_bx.data(), ele_r.rlrs.bvec_x, ele_r.rlrs.bvec_x_ovrf, nb_plus_one,
                    ele_r.dims.inv_scale);
  hostInt95ToDouble(gpt_by.data(), ele_r.rlrs.bvec_y, ele_r.rlrs.bvec_y_ovrf, nb_plus_one,
                    ele_r.dims.inv_scale);
  hostInt95ToDouble(gpt_bz.data(), ele_r.rlrs.bvec_z, ele_r.rlrs.bvec_z_ovrf, nb_plus_one,
                    ele_r.dims.inv_scale);
  hostInt95ToDouble(gpt_cx.data(), ele_r.rlrs.cvec_x, ele_r.rlrs.cvec_x_ovrf, nc_plus_one,
                    ele_r.dims.inv_scale);
  hostInt95ToDouble(gpt_cy.data(), ele_r.rlrs.cvec_y, ele_r.rlrs.cvec_y_ovrf, nc_plus_one,
                    ele_r.dims.inv_scale);
  hostInt95ToDouble(gpt_cz.data(), ele_r.rlrs.cvec_z, ele_r.rlrs.cvec_z_ovrf, nc_plus_one,
                    ele_r.dims.inv_scale);

  // Compute potentials at each grid point using an independent method
  const AtomGraph *ag = ele_bss.getTopologyPointer();
  std::vector<double> lj_check, lj_check_ans, qq_check, qq_check_ans;
  const std::vector<double> sig = ag->getLennardJonesSigma<double>();
  const std::vector<double> eps = ag->getLennardJonesEpsilon<double>();
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const MeshFFKit<T> ele_mnbk = ele_mesh.getNonbondedKit();
  const MeshFFKit<T> vdw_mnbk = vdw_mesh.getNonbondedKit();
  const MeshFFKit<double> ele_mnbk_d = ele_mesh.getReferenceNonbondedKit();
  const MeshFFKit<double> vdw_mnbk_d = vdw_mesh.getReferenceNonbondedKit();
  std::vector<double> lja(nbk.n_lj_types), ljb(nbk.n_lj_types), mesh_ljsig(nbk.n_lj_types);
  std::vector<double> lja_watch(nbk.n_lj_types), ljb_watch(nbk.n_lj_types);
  std::vector<double> ljsig_watch(nbk.n_lj_types);
  for (int i = 0; i < nbk.n_lj_types; i++) {
    double sigfac;
    switch (vdw_mesh.getCombiningRule()) {
    case VdwCombiningRule::GEOMETRIC:
      sigfac = sqrt(sig[i] * vdw_r.probe_radius);
      break;
    case VdwCombiningRule::LORENTZ_BERTHELOT:
      sigfac = 0.5 * (sig[i] + vdw_r.probe_radius);
      break;
    case VdwCombiningRule::NBFIX:
      break;
    }
    mesh_ljsig[i] = sigfac;
    sigfac *= sigfac * sigfac;
    sigfac *= sigfac;
    const double epsfac = 4.0 * sqrt(eps[i] * vdw_r.well_depth);
    lja[i] = epsfac * sigfac * sigfac;
    ljb[i] = epsfac * sigfac;
    lja_watch[i] = vdw_mnbk_d.probe_lja[i];
    ljb_watch[i] = vdw_mnbk_d.probe_ljb[i];
    ljsig_watch[i] = vdw_mnbk_d.probe_ljsig[i];
  }
  check(lja_watch, RelationalOperator::EQUAL, lja, "Probe Lennard-Jones A parameters in the mesh "
        "non-bonded kit do not agree with independent calculations.", tsm.getTestingStatus());
  check(ljb_watch, RelationalOperator::EQUAL, ljb, "Probe Lennard-Jones B parameters in the mesh "
        "non-bonded kit do not agree with independent calculations.", tsm.getTestingStatus());
  check(ljsig_watch, RelationalOperator::EQUAL, mesh_ljsig, "Probe Lennard-Jones pairwise sigma "
        "parameters in the mesh non-bonded kit do not agree with independent calculations.",
        tsm.getTestingStatus());
  const std::vector<bool> mobile_atoms = ag->getAtomMobility();
  const CoordinateFrameReader cfr = ele_bss.getCoordinatePointer()->data();
  for (int i = 0; i < ele_r.dims.na; i++) {
    for (int j = 0; j < ele_r.dims.nb; j++) {
      for (int k = 0; k < ele_r.dims.nc; k++) {
        double u_vdw_raw  = 0.0;
        double u_vdw_soft = 0.0;
        double u_ele_soft = 0.0;
        for (int m = 0; m < cfr.natom; m++) {
          if (mobile_atoms[m]) {
            continue;
          }
          const double gpt_x = gpt_ax[i] + gpt_bx[j] + gpt_cx[k];
          const double gpt_y = gpt_ay[i] + gpt_by[j] + gpt_cy[k];
          const double gpt_z = gpt_az[i] + gpt_bz[j] + gpt_cz[k];
          const double dx = gpt_x - cfr.xcrd[m];
          const double dy = gpt_y - cfr.ycrd[m];
          const double dz = gpt_z - cfr.zcrd[m];
          const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
          const double r = sqrt(r2);
          const double invr  = 1.0 / r;
          const double invr2 = invr * invr;
          const double invr6 = invr2 * invr2 * invr2;
          const int m_idx = nbk.lj_idx[m];
          const double m_q = nbk.charge[m];
          u_vdw_raw += ((lja[m_idx] * invr6) - ljb[m_idx]) * invr6;
          if (r < vdw_mnbk.clash_ratio * vdw_mnbk.probe_ljsig[m_idx]) {
            switch (vdw_r.dims.stencil_kind) {
            case Interpolant::SMOOTHNESS:
              u_vdw_soft += (((((((((vdw_mnbk.softcore_lja[m_idx] * r) +
                                    vdw_mnbk.softcore_ljb[m_idx]) * r) +
                                  vdw_mnbk.softcore_ljc[m_idx]) * r) +
                                vdw_mnbk.softcore_ljd[m_idx]) * r) +
                              vdw_mnbk.softcore_lje[m_idx]) * r) + vdw_mnbk.softcore_ljf[m_idx];
              break;
            case Interpolant::FUNCTION_VALUE:
              u_vdw_soft += (((((vdw_mnbk.softcore_lja[m_idx] * r) +
                                vdw_mnbk.softcore_ljb[m_idx]) * r) +
                              vdw_mnbk.softcore_ljc[m_idx]) * r) + vdw_mnbk.softcore_ljd[m_idx];
              break;
            }
          }
          else {
            u_vdw_soft += ((lja[m_idx] * invr6) - ljb[m_idx]) * invr6;
          }
          if (r < ele_mnbk.clash_distance) {
            switch (ele_r.dims.stencil_kind) {
            case Interpolant::SMOOTHNESS:
              u_ele_soft += m_q * ((((((((((ele_mnbk.softcore_qq[0] * r) +
                                           ele_mnbk.softcore_qq[1]) * r) +
                                         ele_mnbk.softcore_qq[2]) * r) +
                                       ele_mnbk.softcore_qq[3]) * r) +
                                     ele_mnbk.softcore_qq[4]) * r) + ele_mnbk.softcore_qq[5]);
              break;
            case Interpolant::FUNCTION_VALUE:
              u_ele_soft += (((((ele_mnbk.softcore_qq[0] * r) + ele_mnbk.softcore_qq[1]) * r) +
                              ele_mnbk.softcore_qq[2]) * r) + ele_mnbk.softcore_qq[3];
              break;
            }
          }
          else {
            u_ele_soft += m_q * invr;
          }
        }
        if (u_vdw_raw < 64.0) {
          lj_check_ans.push_back(u_vdw_soft);
          qq_check_ans.push_back(u_ele_soft * nbk.coulomb_constant);
          const size_t mesh_coef_idx = 64 * ((((k * ele_r.dims.nb) + j) * ele_r.dims.na) + i);
          lj_check.push_back(vdw_r.coeffs[mesh_coef_idx]);
          qq_check.push_back(ele_r.coeffs[mesh_coef_idx]);
        }
      }
    }
  }
  const std::string mesh_msg = meshDescriptor(ele_r.dims.invu, getStormmScalarTypeName<T>(),
                                              getEnumerationName(vdw_mesh.getBuildPrecision()),
                                              getEnumerationName(vdw_r.dims.stencil_kind));
  double ele_etol, vdw_etol;
  if (sizeof(T) == 4) {
    switch (ele_mesh.getBuildPrecision()) {
    case PrecisionModel::DOUBLE:
      ele_etol = 1.0e-5;
      break;
    case PrecisionModel::SINGLE:
      ele_etol = 9.0e-5;
      break;
    }
    switch (vdw_mesh.getBuildPrecision()) {
    case PrecisionModel::DOUBLE:
      vdw_etol = 3.5e-5;
      break;
    case PrecisionModel::SINGLE:
      vdw_etol = 5.0e-4;
      break;
    }
  }
  else {
    switch (ele_mesh.getBuildPrecision()) {
    case PrecisionModel::DOUBLE:
      ele_etol = 1.0e-6;
      break;
    case PrecisionModel::SINGLE:
      ele_etol = 7.5e-5;
      break;
    }
    switch (vdw_mesh.getBuildPrecision()) {
    case PrecisionModel::DOUBLE:
      vdw_etol = 1.0e-6;
      break;
    case PrecisionModel::SINGLE:
      vdw_etol = 5.0e-4;
      break;
    }
  }
  check(lj_check, RelationalOperator::EQUAL, Approx(lj_check_ans).margin(vdw_etol),
        "Lennard-Jones energies computed for mesh grid points do not match the constant "
        "coefficients for the corresponding mesh elements.  Mesh data type: " +
        getStormmScalarTypeName<T>() + ".  Precision for mesh construction: " +
        getEnumerationName(vdw_mesh.getBuildPrecision()) + ".  Clash ratio: " +
        minimalRealFormat(vdw_mnbk.clash_ratio, 1.0e-4) + ".  " + mesh_msg,
        tsm.getTestingStatus());
  check(qq_check, RelationalOperator::EQUAL, Approx(qq_check_ans).margin(ele_etol),
        "Electrostatic energies computed for mesh grid points do not match the constant "
        "coefficients for the corresponding mesh elements.  Mesh data type: " +
        getStormmScalarTypeName<T>() + ".  Precision for mesh construction: " +
        getEnumerationName(ele_mesh.getBuildPrecision()) + ".  Clash distance: " +
        minimalRealFormat(ele_mnbk.clash_ratio, 1.0e-4) + ".  " + mesh_msg,
        tsm.getTestingStatus());

  // Choose a series of points and test energies
  const double d_na = ele_r.dims.na;
  const double d_nb = ele_r.dims.nb;
  const double d_nc = ele_r.dims.nc;
  const double mesh_orig_x = hostInt95ToDouble(ele_r.dims.orig_x) * ele_r.dims.inv_scale;
  const double mesh_orig_y = hostInt95ToDouble(ele_r.dims.orig_y) * ele_r.dims.inv_scale;
  const double mesh_orig_z = hostInt95ToDouble(ele_r.dims.orig_z) * ele_r.dims.inv_scale;
  const int npts = 64;
  std::vector<double> xdump(npts), ydump(npts), zdump(npts), uvdw_dump(npts), uele_dump(npts);
  std::vector<double> fx_vdw_dump(npts, 0.0), fy_vdw_dump(npts, 0.0), fz_vdw_dump(npts, 0.0);
  std::vector<double> fx_ele_dump(npts, 0.0), fy_ele_dump(npts, 0.0), fz_ele_dump(npts, 0.0);
  for (int i = 0; i < npts; i++) {
    bool high_vdw;
    int pt_aidx, pt_bidx, pt_cidx;
    double ptx, pty, ptz, u_vdw, u_elec, fx_vdw, fy_vdw, fz_vdw, fx_ele, fy_ele, fz_ele;
    do {
      const double a_loc = xrs->uniformRandomNumber() * d_na;
      const double b_loc = xrs->uniformRandomNumber() * d_nb;
      const double c_loc = xrs->uniformRandomNumber() * d_nc;
      ptx = (a_loc * ele_r.dims.invu[0]) + (b_loc * ele_r.dims.invu[3]) +
            (c_loc * ele_r.dims.invu[6]);
      pty = (b_loc * ele_r.dims.invu[4]) + (c_loc * ele_r.dims.invu[7]);
      ptz = (c_loc * ele_r.dims.invu[8]);
      ptx += mesh_orig_x;
      pty += mesh_orig_y;
      ptz += mesh_orig_z;
      high_vdw = false;
      u_vdw  = 0.0;
      u_elec = 0.0;
      fx_ele = 0.0;
      fy_ele = 0.0;
      fz_ele = 0.0;
      fx_vdw = 0.0;
      fy_vdw = 0.0;
      fz_vdw = 0.0;
      for (int j = 0; j < cfr.natom; j++) {
        if (mobile_atoms[j]) {
          continue;
        }
        const double dx = ptx - cfr.xcrd[j];
        const double dy = pty - cfr.ycrd[j];
        const double dz = ptz - cfr.zcrd[j];
        const double r = sqrt((dx * dx) + (dy * dy) + (dz * dz));
        const int j_idx = nbk.lj_idx[j];
        if (r < 2.0 * sig[j_idx]) {
          high_vdw = true;
          break;
        }
        const double invr = 1.0 / r;
        const double invr2 = invr * invr;
        const double invr6 = invr2 * invr2 * invr2;
        u_vdw  += ((lja[j_idx] * invr6) - ljb[j_idx]) * invr6;
        u_elec += nbk.charge[j] * invr;
        const double fmag_ele = -nbk.charge[j] * invr2 * invr;
        fx_ele += fmag_ele * dx;
        fy_ele += fmag_ele * dy;
        fz_ele += fmag_ele * dz;
        const double fmag_vdw = ((6.0 * ljb[j_idx]) - (12.0 * lja[j_idx] * invr6)) * invr6 * invr2;
        fx_vdw += fmag_vdw * dx;
        fy_vdw += fmag_vdw * dy;
        fz_vdw += fmag_vdw * dz;
      }
    } while (high_vdw == true);

    // Collect the points and their exact non-bonded energies to compare to interpolated values
    xdump[i] = ptx;
    ydump[i] = pty;
    zdump[i] = ptz;
    uvdw_dump[i] = u_vdw;
    uele_dump[i] = u_elec * nbk.coulomb_constant;
    fx_vdw_dump[i] = fx_vdw;
    fy_vdw_dump[i] = fy_vdw;
    fz_vdw_dump[i] = fz_vdw;
    fx_ele_dump[i] = fx_ele * nbk.coulomb_constant;
    fy_ele_dump[i] = fy_ele * nbk.coulomb_constant;
    fz_ele_dump[i] = fz_ele * nbk.coulomb_constant;
  }
  const std::vector<T> test_charges(npts, 1.0);
  const std::vector<T> test_lja(1, 4.0 * vdw_r.well_depth * pow(vdw_r.probe_radius, 12.0));
  const std::vector<T> test_ljb(1, 4.0 * vdw_r.well_depth * pow(vdw_r.probe_radius, 6.0));
  const std::vector<int> test_lj_indices(npts, 0);
  ScoreCard sc(1);
  std::vector<T> uvec_ele_save(npts, 0.0), uvec_vdw_save(npts, 0.0);
  std::vector<double> fx_ele_extracted(npts, 0.0), fx_vdw_extracted(npts, 0.0);
  std::vector<double> fy_ele_extracted(npts, 0.0), fy_vdw_extracted(npts, 0.0);
  std::vector<double> fz_ele_extracted(npts, 0.0), fz_vdw_extracted(npts, 0.0);
  const double uele_extracted =
    interpolate<T, T, double, T>(ele_r, ele_mnbk, test_charges.data(), nullptr, nullptr,
                                 xdump.data(), ydump.data(), zdump.data(), nullptr, nullptr,
                                 nullptr, npts, uvec_ele_save.data(), fx_ele_extracted.data(),
                                 fy_ele_extracted.data(), fz_ele_extracted.data(), nullptr,
                                 nullptr, nullptr, OffMeshProtocol::DIE, 0, 0, 0,
                                 sc.getEnergyScaleBits());
  const double uvdw_extracted = 
    interpolate<T, T, double, T>(vdw_r, vdw_mnbk, test_lja.data(), test_ljb.data(), 
                                 test_lj_indices.data(), xdump.data(), ydump.data(), zdump.data(),
                                 nullptr, nullptr, nullptr, npts, uvec_vdw_save.data(),
                                 fx_vdw_extracted.data(), fy_vdw_extracted.data(),
                                 fz_vdw_extracted.data(), nullptr, nullptr, nullptr,
                                 OffMeshProtocol::DIE, 0, 0, 0, sc.getEnergyScaleBits());
  const std::vector<double> mesh_invu = { ele_r.dims.invu[0], ele_r.dims.invu[1],
                                          ele_r.dims.invu[2], ele_r.dims.invu[3],
                                          ele_r.dims.invu[4], ele_r.dims.invu[5],
                                          ele_r.dims.invu[6], ele_r.dims.invu[7],
                                          ele_r.dims.invu[8] };
  std::vector<double> tric_chk_u(npts), tric_chk_fx(npts), tric_chk_fy(npts), tric_chk_fz(npts);
  for (int i = 0; i < npts; i++) {

    // Check the results of the mesh with a minimal calculation on a TricubicCell object
    // representing just one of its elements.  Use the electrostatics case.
    const double rel_x = xdump[i] - mesh_orig_x;
    const double rel_y = ydump[i] - mesh_orig_y;
    const double rel_z = zdump[i] - mesh_orig_z;
    const int exmp_cell_a = (rel_x * ele_r.dims.umat[0]) + (rel_y * ele_r.dims.umat[3]) +
                            (rel_z * ele_r.dims.umat[6]);
    const int exmp_cell_b = (rel_y * ele_r.dims.umat[4]) + (rel_z * ele_r.dims.umat[7]);
    const int exmp_cell_c = (rel_z * ele_r.dims.umat[8]);
    const double3 exmp_origin = { (static_cast<double>(exmp_cell_a) * ele_r.dims.invu[0]) +
                                  (static_cast<double>(exmp_cell_b) * ele_r.dims.invu[3]) +
                                  (static_cast<double>(exmp_cell_c) * ele_r.dims.invu[6]) +
                                  mesh_orig_x,
                                  (static_cast<double>(exmp_cell_b) * ele_r.dims.invu[4]) +
                                  (static_cast<double>(exmp_cell_c) * ele_r.dims.invu[7]) +
                                  mesh_orig_y,
                                  (static_cast<double>(exmp_cell_c) * ele_r.dims.invu[8]) +
                                  mesh_orig_z };
    const double3 test_pt = { xdump[i], ydump[i], zdump[i] };
    const double4 utcc = checkTricubicElectrostatics(*ag, *(ele_mesh.getCoordinatePointer()),
                                                     exmp_origin, mesh_invu, test_pt);
    tric_chk_u[i]  = utcc.w;
    tric_chk_fx[i] = utcc.x;
    tric_chk_fy[i] = utcc.y;
    tric_chk_fz[i] = utcc.z;
  }
  check(sum<double>(uvec_ele_save), RelationalOperator::EQUAL,
        Approx(uele_extracted).margin(1.0e-6), "The sum of observations on an electrostatic mesh "
        "does not agree with the quantity reported by the overall interpolation process.",
        tsm.getTestingStatus());
  check(sum<double>(uvec_vdw_save), RelationalOperator::EQUAL,
        Approx(uvdw_extracted).margin(1.e-6), "The sum of observations on a Lennard-Jones mesh "
        "does not agree with the quantity reported by the overall interpolation process.",
        tsm.getTestingStatus());
  check(uvdw_extracted, RelationalOperator::EQUAL, Approx(sum<double>(uvdw_dump)).margin(1.5e-3),
        "Van-der Waals (Lennard-Jones) energies interpolated for a TIP3P-like solvent probe do "
        "not meet expectations.  " + mesh_msg, tsm.getTestingStatus());
  check(uele_extracted, RelationalOperator::EQUAL, Approx(sum<double>(uele_dump)).margin(1.0e-3),
        "Electrostatic energies interpolated from a mesh do not meet expectations.  " + mesh_msg,
        tsm.getTestingStatus());
  check(sum<double>(tric_chk_u), RelationalOperator::EQUAL,
        Approx(sum<double>(uele_dump)).margin(6.0e-3), "Electrostatic energies interpolated in "
        "the tricubic model cell do not meet expectations.  " + mesh_msg, tsm.getTestingStatus());
  check(tric_chk_fx, RelationalOperator::EQUAL, Approx(fx_ele_extracted).margin(1.0e-3),
        "Electrostatic X-axis forces interpolated in the tricubic model cell do not meet "
        "expectations.  " + mesh_msg, tsm.getTestingStatus());
  check(tric_chk_fy, RelationalOperator::EQUAL, Approx(fy_ele_extracted).margin(1.0e-3),
        "Electrostatic Y-axis forces interpolated in the tricubic model cell do not meet "
        "expectations.  " + mesh_msg, tsm.getTestingStatus());
  check(tric_chk_fz, RelationalOperator::EQUAL, Approx(fz_ele_extracted).margin(1.0e-3),
        "Electrostatic Z-axis forces interpolated in the tricubic model cell do not meet "
        "expectations.  " + mesh_msg, tsm.getTestingStatus());
  check(fx_ele_dump, RelationalOperator::EQUAL, Approx(fx_ele_extracted).margin(5.0e-3),
        "Electrostatic X-axis forces interpolated from a mesh do not meet expectations.  " +
        mesh_msg, tsm.getTestingStatus());
  check(fy_ele_dump, RelationalOperator::EQUAL, Approx(fy_ele_extracted).margin(5.0e-3),
        "Electrostatic Y-axis forces interpolated from a mesh do not meet expectations.  " +
        mesh_msg, tsm.getTestingStatus());
  check(fz_ele_dump, RelationalOperator::EQUAL, Approx(fz_ele_extracted).margin(5.0e-3),
        "Electrostatic Z-axis forces interpolated from a mesh do not meet expectations.  " +
        mesh_msg, tsm.getTestingStatus());
  check(fx_vdw_dump, RelationalOperator::EQUAL, Approx(fx_vdw_extracted).margin(7.0e-4),
        "Van-der Waals X-axis forces interpolated from a mesh do not meet expectations.  " +
        mesh_msg, tsm.getTestingStatus());
  check(fy_vdw_dump, RelationalOperator::EQUAL, Approx(fy_vdw_extracted).margin(7.0e-4),
        "Van-der Waals Y-axis forces interpolated from a mesh do not meet expectations.  " +
        mesh_msg, tsm.getTestingStatus());
  check(fz_vdw_dump, RelationalOperator::EQUAL, Approx(fz_vdw_extracted).margin(7.0e-4),
        "Van-der Waals Z-axis forces interpolated from a mesh do not meet expectations.  " +
        mesh_msg, tsm.getTestingStatus());
}

//-------------------------------------------------------------------------------------------------
// Check the bits in an occlusion map using independent mapping methods.
//
// Arguments:
//   lattice:  The mesh of bitstrings indicating whether space is occluded or not
//-------------------------------------------------------------------------------------------------
void checkBitSettings(const BackgroundMesh<ullint> &lattice, const std::string &platform,
                      const TestPriority do_tests) {
  const BackgroundMeshReader<ullint> lattice_r = lattice.data();
  const MeshParamKit mps = lattice_r.dims;
  const MeshFoundation& lattice_bss = lattice.getMolecularBasis();
  const NonbondedKit<double> nbk =
    lattice_bss.getTopologyPointer()->getDoublePrecisionNonbondedKit();
  const CoordinateFrame *cf = lattice_bss.getCoordinatePointer();
  const CoordinateFrameReader cfr = cf->data();
  const double dorig_x = hostInt95ToDouble(mps.orig_x) * mps.inv_scale;
  const double dorig_y = hostInt95ToDouble(mps.orig_y) * mps.inv_scale;
  const double dorig_z = hostInt95ToDouble(mps.orig_z) * mps.inv_scale;
  const double probe_rad = lattice.getProbeRadius();

  // Get a vector of the mobile atoms.
  const std::vector<uint> frozen_mask = lattice_bss.getFrozenAtomMask();
  std::vector<bool> frozen_atoms(nbk.natom, false);
  for (int i = 0; i < nbk.natom; i++) {
    if (readBitFromMask(frozen_mask, i) == 1) {
      frozen_atoms[i] = true;
    }
  }

  // Access the colored elements of the mesh and confirm that they all lie within the proper
  // radius of some atom.
  std::vector<double> color_radii(nbk.natom);
  for (int pos = 0; pos < nbk.natom; pos++) {
    color_radii[pos] = (0.5 * nbk.lj_sigma[(nbk.n_lj_types + 1) * nbk.lj_idx[pos]]) + probe_rad;
  }
  std::vector<int> short_list;
  short_list.reserve(16);
  std::vector<double> best_failed_ranges;
  int n_valid_pts = 0;
  int n_atomless = 0;
  for (int i = 0; i < mps.na; i++) {
    for (int j = 0; j < mps.nb; j++) {
      for (int k = 0; k < mps.nc; k++) {

        // Make a short list of atoms that are within striking distance of this element.
        const double di = i;
        const double dj = j;
        const double dk = k;
        const double dpi = di + 0.5;
        const double dpj = dj + 0.5;
        const double dpk = dk + 0.5;
        const double elem_midx = dorig_x + 
                                 (dpi * mps.invu[0]) + (dpj * mps.invu[3]) + (dpk * mps.invu[6]);
        const double elem_midy = dorig_y +
                                 (dpi * mps.invu[1]) + (dpj * mps.invu[4]) + (dpk * mps.invu[7]);
        const double elem_midz = dorig_z +
                                 (dpi * mps.invu[2]) + (dpj * mps.invu[5]) + (dpk * mps.invu[8]);
        short_list.resize(0);
        for (int pos = 0; pos < nbk.natom; pos++) {
          const double dx = elem_midx - cfr.xcrd[pos];
          const double dy = elem_midy - cfr.ycrd[pos];
          const double dz = elem_midz - cfr.zcrd[pos];
          if (frozen_atoms[pos] &&
              sqrt((dx * dx) + (dy * dy) + (dz * dz)) < color_radii[pos] + (0.5 * mps.max_span)) {
            short_list.push_back(pos);
          }
        }
        const int nshort = short_list.size();

        // Loop over the 4096 mesh points in this element.
        const int coeff_idx = 64 * ((((k * mps.nb) + j) * mps.na) + i);
        for (int m = 0; m < 64; m++) {
          if (lattice_r.coeffs[coeff_idx + m] == 0LLU) {
            continue;
          }
          const int cube_k = m / 16;
          const int cube_j = (m - (cube_k * 16)) / 4;
          const int cube_i = m - (cube_k * 16) - (cube_j * 4);
          const double dcbi = di + (static_cast<double>(cube_i) * 0.25);
          const double dcbj = dj + (static_cast<double>(cube_j) * 0.25);
          const double dcbk = dk + (static_cast<double>(cube_k) * 0.25);
          const ullint lrcm = lattice_r.coeffs[coeff_idx + m];
          for (int n = 0; n < 64; n++) {
            if ((lrcm >> n) & 0x1) {

              // This mesh point is colored.  Find at least one atom that might be near enough.
              const int cubelet_k = n / 16;
              const int cubelet_j = (n - (cubelet_k * 16)) / 4;
              const int cubelet_i = n - (cubelet_k * 16) - (cubelet_j * 4); 
              const double dcbli = dcbi + ((static_cast<double>(cubelet_i) + 0.5) * 0.0625);
              const double dcblj = dcbj + ((static_cast<double>(cubelet_j) + 0.5) * 0.0625);
              const double dcblk = dcbk + ((static_cast<double>(cubelet_k) + 0.5) * 0.0625);
              const double pt_x = dorig_x + (dcbli * mps.invu[0]) + (dcblj * mps.invu[3]) +
                                  (dcblk * mps.invu[6]);
              const double pt_y = dorig_y + (dcbli * mps.invu[1]) + (dcblj * mps.invu[4]) +
                                  (dcblk * mps.invu[7]);
              const double pt_z = dorig_z + (dcbli * mps.invu[2]) + (dcblj * mps.invu[5]) +
                                  (dcblk * mps.invu[8]);
              bool in_range = false;
              double best_range = 1.0e30;
              for (int pos = 0; pos < nshort; pos++) {
                const int pidx = short_list[pos];
                const double dx = pt_x - cfr.xcrd[pidx];
                const double dy = pt_y - cfr.ycrd[pidx];
                const double dz = pt_z - cfr.zcrd[pidx];
                const double t_range = sqrt((dx * dx) + (dy * dy) + (dz * dz));
                if (t_range < color_radii[pidx]) {
                  in_range = true;
                  break;
                }
                else {
                  best_range = std::min(best_range, t_range - color_radii[pidx]);
                }
              }
              if (in_range) {
                n_valid_pts++;
              }
              else if (nshort > 0) {
                best_failed_ranges.push_back(best_range);
              }
              else {
                n_atomless++;
              }
            }
          }
        }
      }
    }
  }
  std::string max_failed_range, mean_failed_range, stdev_failed_range;
  double max_failed_rng = 0.0;
  double mean_failed_rng = 0.0;
  double stdev_failed_rng = 0.0;
  if (best_failed_ranges.size() > 0) {
    max_failed_rng = maxValue(best_failed_ranges);
    mean_failed_rng = mean(best_failed_ranges);
    stdev_failed_rng = variance(best_failed_ranges, VarianceMethod::STANDARD_DEVIATION);
    max_failed_range = realToString(max_failed_rng, 9, 4, NumberFormat::STANDARD_REAL);
    mean_failed_range = realToString(mean_failed_rng, 9, 4, NumberFormat::STANDARD_REAL);
    stdev_failed_range = realToString(stdev_failed_rng, 9, 4, NumberFormat::STANDARD_REAL);
  }

  // In single-precision computations on the GPU, roundoff error may create the illusion of errors
  // in the test.  Tolerate a very small degree of failed ranges.
  if (best_failed_ranges.size() < n_valid_pts / ok_point_ratio && stdev_failed_rng < 1.0e-4) {
    best_failed_ranges.resize(0);
  }
  check(best_failed_ranges.size() == 0, "A total of " + std::to_string(best_failed_ranges.size()) +
        " points on the mesh were not found in be in range of any atoms.  The maximum discrepancy "
        "in ranges was " + max_failed_range + ", with mean value " + mean_failed_range + " +/- " +
        stdev_failed_range + ".  There are " + std::to_string(n_valid_pts) + " valid points "
        "colored on the mesh.  There are " + std::to_string(n_atomless) + " points colored that "
        "could not possibly be in range of any atoms.  Topology involved: " +
        getBaseName(lattice_bss.getTopologyPointer()->getFileName()) + ".  The mesh was computed "
        "by the " + platform + ".", do_tests);

  // Loop over atoms and recolor a separate array for comparison.  Check that this procedure does
  // not mark any points that were unmarked in the original mesh.
  std::vector<bool> refmask(mps.na * mps.nb * mps.nc * 64 * 64, false);
  std::vector<double> unmarked;
  for (int pos = 0; pos < nbk.natom; pos++) {
    if (frozen_atoms[pos] == false) {
      continue;
    }
    const double atom_x = cfr.xcrd[pos];
    const double atom_y = cfr.ycrd[pos];
    const double atom_z = cfr.zcrd[pos];
    const double sq_radius = color_radii[pos] * color_radii[pos];
    const double dx = atom_x - dorig_x;
    const double dy = atom_y - dorig_y;
    const double dz = atom_z - dorig_z;
    const int gcenx = floor(16.0 * ((mps.umat[0] * dx) + (mps.umat[3] * dy) + (mps.umat[6] * dz)));
    const int gceny = floor(16.0 * ((mps.umat[1] * dx) + (mps.umat[4] * dy) + (mps.umat[7] * dz)));
    const int gcenz = floor(16.0 * ((mps.umat[2] * dx) + (mps.umat[5] * dy) + (mps.umat[8] * dz)));
    const int gwide_x = ceil(16.0 * color_radii[pos] / mps.widths[0]);
    const int gwide_y = ceil(16.0 * color_radii[pos] / mps.widths[1]);
    const int gwide_z = ceil(16.0 * color_radii[pos] / mps.widths[2]);
    const int gmin_x = std::max(0, gcenx - gwide_x);
    const int gmin_y = std::max(0, gceny - gwide_y);
    const int gmin_z = std::max(0, gcenz - gwide_z);
    const int gmax_x = std::min(gcenx + gwide_x + 1, 16 * mps.na);
    const int gmax_y = std::min(gceny + gwide_y + 1, 16 * mps.nb);
    const int gmax_z = std::min(gcenz + gwide_z + 1, 16 * mps.nc);
    for (int i = gmin_x; i < gmax_x; i++) {
      const int element_i = i / 16;
      const int cube_i = (i - (element_i * 16)) / 4;
      const int cubelet_i = i - (element_i * 16) - (cube_i * 4);
      const double gu_x = (static_cast<double>(i) + 0.5) * 0.0625;
      for (int j = gmin_y; j < gmax_y; j++) {
        const int element_j = j / 16;
        const int cube_j = (j - (element_j * 16)) / 4;
        const int cubelet_j = j - (element_j * 16) - (cube_j * 4);
        const double gu_y = (static_cast<double>(j) + 0.5) * 0.0625;
        for (int k = gmin_z; k < gmax_z; k++) {
          const double gu_z = (static_cast<double>(k) + 0.5) * 0.0625;
          const double pt_x = dorig_x + (mps.invu[0] * gu_x) + (mps.invu[3] * gu_y) +
                              (mps.invu[6] * gu_z);
          const double pt_y = dorig_y + (mps.invu[1] * gu_x) + (mps.invu[4] * gu_y) +
                              (mps.invu[7] * gu_z);
          const double pt_z = dorig_z + (mps.invu[2] * gu_x) + (mps.invu[5] * gu_y) +
                              (mps.invu[8] * gu_z);
          const double dx = pt_x - atom_x;
          const double dy = pt_y - atom_y;
          const double dz = pt_z - atom_z;
          if ((dx * dx) + (dy * dy) + (dz * dz) < sq_radius) {
            const int element_k = k / 16;
            const int cube_k = (k - (element_k * 16)) / 4;
            const int cubelet_k = k - (element_k * 16) - (cube_k * 4);
            const int element_idx = (((element_k * mps.nb) + element_j) * mps.na) + element_i;
            const int cube_idx = (((cube_k * 4) + cube_j) * 4) + cube_i;
            const int cubelet_idx = (((cubelet_k * 4) + cubelet_j) * 4) + cubelet_i;
            refmask[(((element_idx * 64) + cube_idx) * 64) + cubelet_idx] = true;
            
            // Check against the original mesh
            const ullint t_coeff = lattice_r.coeffs[(element_idx * 64) + cube_idx];
            if (((t_coeff >> cubelet_idx) & 0x1) == 0) {
              unmarked.push_back(sqrt((dx * dx) + (dy * dy) + (dz * dz)));
            }
          }
        }
      }
    }
  }
  
  // Check that there are no points marked by the above procedure that were not marked in the
  // original mesh.
  std::string min_unmarked_range, max_unmarked_range, mean_unmarked_range, stdev_unmarked_range;
  double stdev_unmarked_rng = 0.0;
  if (unmarked.size() > 0) {
    min_unmarked_range = realToString(minValue(unmarked), 9, 4, NumberFormat::STANDARD_REAL);
    max_unmarked_range = realToString(maxValue(unmarked), 9, 4, NumberFormat::STANDARD_REAL);
    mean_unmarked_range = realToString(mean(unmarked), 9, 4, NumberFormat::STANDARD_REAL);
    stdev_unmarked_rng = variance(unmarked, VarianceMethod::STANDARD_DEVIATION);
    stdev_unmarked_range = realToString(stdev_unmarked_rng, 9, 4, NumberFormat::STANDARD_REAL);
  }
  if (mean_unmarked_range.size() < n_valid_pts / ok_point_ratio &&
      stdev_unmarked_rng < 1.0e-4) {
    unmarked.resize(0);
  }
  check(unmarked.size() == 0, "The independent checking procedure found " +
        std::to_string(unmarked.size()) + " points that were not marked by the mesh generation "
        "for topology " + getBaseName(lattice_bss.getTopologyPointer()->getFileName()) + ".  The "
        "points have distances of " + min_unmarked_range + " to " + max_unmarked_range +
        " Angstroms from relevant atoms with a mean distance of " + mean_unmarked_range + " +/- " +
        stdev_unmarked_range + " Angstroms.  The mesh was computed by the " + platform + ".",
        do_tests);

  // Check that the original mesh does not contain marked points that the independent procedure
  // did not find.
  int n_overmarked = 0;
  int n_undermarked = 0;
  int n_fence = 0;
  std::vector<double> overmarked_ranges;
  std::vector<double> undermarked_ranges;
  for (int i = 0; i < mps.na; i++) {
    for (int j = 0; j < mps.nb; j++) {
      for (int k = 0; k < mps.nc; k++) {
        const int element_idx = (((k * mps.nb) + j) * mps.na) + i;
        const int element_offset = 64 * element_idx;
        for (int m = 0; m < 64; m++) {
          if (lattice_r.coeffs[element_offset + m] == 0LLU) {
            continue;
          }
          const ullint t_coeff = lattice_r.coeffs[element_offset + m];
          for (int n = 0; n < 64; n++) {
            if (((t_coeff >> n) & 0x1) && refmask[(((element_idx * 64) + m) * 64) + n] == false) {
              
              // Should the point have been marked, in fact?
              const int cube_k = m / 16;
              const int cube_j = (m - (16 * cube_k)) / 4;
              const int cube_i = (m & 3);
              const int cubelet_k = n / 16;
              const int cubelet_j = (n - (16 * cubelet_k)) / 4;
              const int cubelet_i = (n & 3);
              const double ginc_a = i + (0.25 * static_cast<double>(cube_i)) +
                                    (0.0625 * static_cast<double>(cubelet_i));
              const double ginc_b = j + (0.25 * static_cast<double>(cube_j)) +
                                    (0.0625 * static_cast<double>(cubelet_j));
              const double ginc_c = k + (0.25 * static_cast<double>(cube_k)) +
                                    (0.0625 * static_cast<double>(cubelet_k));
              const double pt_x = dorig_x + (ginc_a * mps.invu[0]) + (ginc_b * mps.invu[3]) +
                                  (ginc_c * mps.invu[6]);
              const double pt_y = dorig_y + (ginc_a * mps.invu[1]) + (ginc_b * mps.invu[4]) +
                                  (ginc_c * mps.invu[7]);
              const double pt_z = dorig_z + (ginc_a * mps.invu[2]) + (ginc_b * mps.invu[5]) +
                                  (ginc_c * mps.invu[8]);
              double min_discrepancy = 0.0;
              for (int pos = 0; pos < nbk.natom; pos++) {
                if (frozen_atoms[pos] == false) {
                  continue;
                }
                const double dx = pt_x - cfr.xcrd[pos];
                const double dy = pt_y - cfr.ycrd[pos];
                const double dz = pt_z - cfr.zcrd[pos];
                const double r = sqrt((dx * dx) + (dy * dy) + (dz * dz));
                min_discrepancy = (pos == 0) ? r - color_radii[pos] :
                                               std::min(min_discrepancy, r - color_radii[pos]);
              }
              if (min_discrepancy > 0.0) {
                n_overmarked++;
                overmarked_ranges.push_back(min_discrepancy);
              }
              else if (min_discrepancy < 0.0) {
                n_undermarked++;
                undermarked_ranges.push_back(min_discrepancy);
              }
              else {
                n_fence++;
              }
            }
          }
        }
      }
    }
  }
  std::string overmarked_max, overmarked_mean, overmarked_stdev;
  double overmarked_stdev_val = 0.0;
  if (overmarked_ranges.size() > 0) {
    overmarked_max = realToString(maxValue(overmarked_ranges), 9, 4, NumberFormat::STANDARD_REAL);
    overmarked_mean = realToString(mean(overmarked_ranges), 9, 4, NumberFormat::STANDARD_REAL);
    overmarked_stdev_val = variance(overmarked_ranges, VarianceMethod::STANDARD_DEVIATION);
    overmarked_stdev = realToString(overmarked_stdev_val, 9, 4, NumberFormat::STANDARD_REAL);
  }
  std::string undermarked_max, undermarked_mean, undermarked_stdev;
  double undermarked_stdev_val = 0.0;
  if (undermarked_ranges.size() > 0) {
    undermarked_max = realToString(maxAbsValue(undermarked_ranges), 9, 4,
                                   NumberFormat::STANDARD_REAL);
    undermarked_mean = realToString(fabs(mean(undermarked_ranges)), 9, 4,
                                    NumberFormat::STANDARD_REAL);
    undermarked_stdev_val = variance(undermarked_ranges, VarianceMethod::STANDARD_DEVIATION);
    undermarked_stdev = realToString(undermarked_stdev_val, 9, 4, NumberFormat::STANDARD_REAL);
  }
  if (n_overmarked < n_valid_pts / ok_point_ratio && overmarked_stdev_val < 1.0e-4) {
    n_overmarked = 0;
  }
  check(n_overmarked == 0, "A total of " + std::to_string(n_overmarked) + " points were marked as "
        "occluded in the mesh for topology " +
        getBaseName(lattice_bss.getTopologyPointer()->getFileName()) + ", but not found in the "
        "independent procedure.  The over-marked points extend as far as " + overmarked_max +
        " Angstroms from any atoms' valid occlusive boundaries, with a mean value of " +
        overmarked_mean + " +/- " + overmarked_stdev + " Angstroms.  The mesh was computed by "
        "the " + platform + ".", do_tests);
  if (n_undermarked < n_valid_pts / ok_point_ratio && undermarked_stdev_val < 1.0e-4) {
    n_undermarked = 0;
  }
  check(n_undermarked == 0, "A total of " + std::to_string(n_undermarked) + " points were marked "
        "as occluded in the mesh for topology " +
        getBaseName(lattice_bss.getTopologyPointer()->getFileName()) + ", but apparently missed "
        "by the independent procedure.  The test itself may need to be debugged.  The "
        "under-marked points extend as far as " + undermarked_max + " Angstroms inside of the "
        "atoms' valid occlusive boundaries, with a mean value of " + undermarked_mean + " +/- " +
        undermarked_stdev + " Angstroms.  The mesh was computed by the " + platform + ".",
        do_tests);
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  StopWatch timer("Mesh Generation Testing");
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1
  section("Establish an occlusion mesh object");
  
  // Section 2
  section("Establish a non-bonded field mesh object");

  // Read the topologies for systems of various sizes and complexities.  Set their atoms to be
  // mobile.
  char osc = osSeparator();
  const std::string crd_base_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::string top_base_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::vector<std::string> molecule_names = { "trpcage", "sodium", "chloride", "nacl",
                                                    "drug_example_dry", "tamavidin" };
  TestSystemManager tsm(top_base_name, "top", molecule_names, crd_base_name, "inpcrd",
                        molecule_names);
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    AtomGraph *iag_ptr = tsm.getTopologyPointer(i);
    iag_ptr->modifyAtomMobility(0, iag_ptr->getAtomCount(), MobilitySetting::OFF);
#ifdef STORMM_USE_HPC
    iag_ptr->upload();
#endif
  }

  // Add timings categories
  timer.addCategory("Trial, Trp-Cage Occlusion Mesh (CPU)");
  timer.addCategory("Occlusion Mesh Checking");
  timer.addCategory("Trial, Trp-Cage Electrostatic Mesh (CPU)");
  timer.addCategory("Trial, Trp-Cage van-der Waals Mesh (CPU)");
  timer.addCategory("Trial, Na+ Occlusion Mesh (CPU)");
  timer.addCategory("Trial, Cl- Occlusion Mesh (CPU)");
  timer.addCategory("Trial, NaCl Occlusion Mesh (CPU)");
  timer.addCategory("Trial, Na+ Tric. Occ. Mesh (CPU)");
  timer.addCategory("Trial, NaCl Tric. Occ. Mesh (CPU)");
  timer.addCategory("Trial, Tamavidin Elec. Mesh (CPU)");
#ifdef STORMM_USE_HPC
  timer.addCategory("Trial, Trp-Cage Occlusion Mesh (GPU)");
  timer.addCategory("Trial, Na+ Occlusion Mesh (GPU)");
  timer.addCategory("Trial, Cl- Occlusion Mesh (GPU)");
  timer.addCategory("Trial, NaCl Occlusion Mesh (GPU)");
  timer.addCategory("Trial, NaCl Tric. Occ. Mesh (GPU)");
  timer.addCategory("Trial, Trp-Cage Tric. Occ. Mesh (GPU)");
  timer.addCategory("Trial, Tamavidin Elec. Mesh (GPU)");
#endif
  // Make an exclusion mesh out of the most complex system.
  section(1);
  AtomGraph *trpi_ag = tsm.getTopologyPointer(0);
  trpi_ag->modifyAtomMobility(0, trpi_ag->getAtomCount(), MobilitySetting::OFF);
  const std::vector<bool> trpi_mobile = trpi_ag->getAtomMobility();
  CoordinateFrame trpi_cf = tsm.exportCoordinateFrame(0);
  timer.assignTime(0);
  BackgroundMesh<ullint> bgm_a(GridDetail::OCCLUSION, NonbondedPotential::CLASH);
  timer.assignTime("Trial, Trp-Cage Occlusion Mesh (CPU)");
  MeshFoundation *bgm_a_bss = bgm_a.getMolecularBasis();
  bgm_a_bss->setTopology(tsm.getTopologyPointer(0));
  bgm_a_bss->setCoordinates(&trpi_cf);
  bgm_a.setMeshParameters(10.0, 2.5);
  bgm_a.setProbeRadius(1.4);
  bgm_a.computeField();
  BackgroundMeshWriter<ullint> bgmw_a = bgm_a.data();
  timer.assignTime(0);
  checkBitSettings(bgm_a, "CPU", tsm.getTestingStatus());
  timer.assignTime("Occlusion Mesh Checking");
  const std::vector<double> trpi_mesh_dims = { static_cast<double>(bgmw_a.dims.na),
                                               static_cast<double>(bgmw_a.dims.nb),
                                               static_cast<double>(bgmw_a.dims.nc),
                                               bgmw_a.dims.invu[0], bgmw_a.dims.invu[1],
                                               bgmw_a.dims.invu[2], bgmw_a.dims.invu[3],
                                               bgmw_a.dims.invu[4], bgmw_a.dims.invu[5],
                                               bgmw_a.dims.invu[6], bgmw_a.dims.invu[7],
                                               bgmw_a.dims.invu[8] };
  const std::vector<double> trpi_mesh_dims_ans = { 18.0, 18.0, 15.0,  2.5,  0.0,  0.0,
                                                    0.0,  2.5,  0.0,  0.0,  0.0,  2.5 };
  check(trpi_mesh_dims, RelationalOperator::EQUAL, trpi_mesh_dims_ans, "The mesh dimensions "
        "computed for Trp-cage do not meet expectations.", tsm.getTestingStatus());
  timer.assignTime(0);
  
  // Make some non-bonded field meshes out of the most complex system.  Use no softcore potential
  // for the first tests, then implement a significant softcore potential and test whether its
  // values are properly computed.
  section(2);
  BackgroundMesh<double> trpi_ortho_ele(GridDetail::NONBONDED_FIELD,
                                        NonbondedPotential::ELECTROSTATIC,
                                        tsm.getTopologyPointer(0), &trpi_cf, 2.5, 1.5);
  std::vector<double> umat_diamond(9), invu_diamond(9), dims_diamond(6);
  for (int i = 0; i < 3; i++) {
    dims_diamond[i    ] = 1.5;
    dims_diamond[i + 3] = stormm::symbols::tetrahedral_angle;
  }
  computeBoxTransform(dims_diamond, &umat_diamond, &invu_diamond);
  BackgroundMesh<double> trpi_tricl_ele(GridDetail::NONBONDED_FIELD,
                                        NonbondedPotential::ELECTROSTATIC,
                                        tsm.getTopologyPointer(0), &trpi_cf, 2.5, invu_diamond,
                                        40, 1.0, PrecisionModel::DOUBLE);
  timer.assignTime("Trial, Trp-Cage Electrostatic Mesh (CPU)");
  BackgroundMesh<double> trpi_ortho_vdw(GridDetail::NONBONDED_FIELD,
                                        NonbondedPotential::VAN_DER_WAALS, 3.15061, 0.1521,
                                        VdwCombiningRule::LORENTZ_BERTHELOT,
                                        tsm.getTopologyPointer(0), &trpi_cf, 2.5, 1.5, 40, {}, {},
                                        0.0, 0.8, PrecisionModel::DOUBLE);
  BackgroundMesh<double> trpi_tricl_vdw(GridDetail::NONBONDED_FIELD,
                                        NonbondedPotential::VAN_DER_WAALS, 3.15061, 0.1521,
                                        VdwCombiningRule::LORENTZ_BERTHELOT,
                                        tsm.getTopologyPointer(0), &trpi_cf, 2.5, invu_diamond, 40,
                                        {}, {}, 0.0, 0.8, PrecisionModel::DOUBLE);
  timer.assignTime("Trial, Trp-Cage van-der Waals Mesh (CPU)");
  timer.assignTime("Trial, Trp-Cage Electrostatic Mesh (CPU)");

  // Check the non-bonded energies at mesh points
  Xoshiro256ppGenerator xrs;
  checkNonbondedPotentials<double>(trpi_ortho_ele, trpi_ortho_vdw, &xrs, tsm);
  checkNonbondedPotentials<double>(trpi_tricl_ele, trpi_tricl_vdw, &xrs, tsm);

#ifdef STORMM_USE_HPC  
  HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  const MeshKlManager mesh_launcher(gpu);
#endif
  
  // Create Lennard-Jones meshes with various combining rules, data types, and interpolants.  Test
  // them in conjunction with new electrostatic meshes.
  CoordinateFrame drug_mol = tsm.exportCoordinateFrame(4);
  const CoordinateFrameReader drug_molr = drug_mol.data();
  std::vector<PrecisionModel> prec_mdl = { PrecisionModel::DOUBLE, PrecisionModel::SINGLE };
  std::vector<Interpolant> fitting = { Interpolant::SMOOTHNESS, Interpolant::FUNCTION_VALUE };
  for (size_t i = 0; i < prec_mdl.size(); i++) {
    for (size_t j = 0; j < fitting.size(); j++) {
      MeshParameters drug_mps(10, 10, 11, drug_molr.xcrd[drug_molr.natom / 2],
                              drug_molr.ycrd[drug_molr.natom / 2],
                              drug_molr.zcrd[drug_molr.natom / 2], invu_diamond, 40, fitting[j]);
      BackgroundMesh<float> drug_ljgeom_cpu(GridDetail::NONBONDED_FIELD,
                                            NonbondedPotential::VAN_DER_WAALS,
                                            tsm.getTopologyPointer(4), &drug_mol, 1.44, 0.16,
                                            VdwCombiningRule::GEOMETRIC, drug_mps, {}, {}, 0.0,
                                            1.0, prec_mdl[i]);
      BackgroundMesh<float> drug_elec_cpu(GridDetail::NONBONDED_FIELD,
                                          NonbondedPotential::ELECTROSTATIC,
                                          tsm.getTopologyPointer(4), &drug_mol, 1.4, 0.0,
                                          VdwCombiningRule::GEOMETRIC, drug_mps, {}, {}, 1.0, 0.0,
                                          prec_mdl[i]);
      checkNonbondedPotentials<float>(drug_elec_cpu, drug_ljgeom_cpu, &xrs, tsm);
#ifdef STORMM_USE_HPC
      BackgroundMesh<float> drug_ljgeom_gpu(GridDetail::NONBONDED_FIELD,
                                            NonbondedPotential::VAN_DER_WAALS,
                                            tsm.getTopologyPointer(4), &drug_mol, 1.44, 0.16,
                                            VdwCombiningRule::GEOMETRIC, drug_mps, {}, {}, 0.0,
                                            1.0, prec_mdl[i], mesh_launcher,
                                            HybridTargetLevel::HOST);
      BackgroundMesh<float> drug_elec_gpu(GridDetail::NONBONDED_FIELD,
                                          NonbondedPotential::ELECTROSTATIC,
                                          tsm.getTopologyPointer(4), &drug_mol, 1.4, 0.0,
                                          VdwCombiningRule::GEOMETRIC, drug_mps, {}, {}, 1.0, 0.0,
                                          prec_mdl[i], mesh_launcher, HybridTargetLevel::HOST);
      checkNonbondedPotentials<float>(drug_elec_gpu, drug_ljgeom_gpu, &xrs, tsm);
#endif
    }
  }
  
  // Assemble various occlusion meshes
  CoordinateFrame sodium_cf = tsm.exportCoordinateFrame(1);
  CHECK_THROWS_SOFT(BackgroundMesh<ullint> bad_news(GridDetail::OCCLUSION, 3.0,
                                                    VdwCombiningRule::LORENTZ_BERTHELOT,
                                                    tsm.getTopologyPointer(3), &sodium_cf,
                                                    { -5.0, -5.0, -5.0, 6.25, 6.25, 6.25 }, 1.125),
                    "A mesh was constructed with inconsistent numbers of atoms in the topology "
                    "and coordinate object.", tsm.getTestingStatus());
  
  // Create additional occlusion meshes with orthogonal grid axes
  CoordinateFrame chloride_cf = tsm.exportCoordinateFrame(2);
  CoordinateFrame nacl_cf = tsm.exportCoordinateFrame(3);
  timer.assignTime(0);
  BackgroundMesh<ullint> sodium_occm_cpu(GridDetail::OCCLUSION, 3.0,
                                         VdwCombiningRule::LORENTZ_BERTHELOT,
                                         tsm.getTopologyPointer(1), &sodium_cf,
                                         { -5.0, -5.0, -5.0, 6.25, 6.25, 6.25 }, 1.125);
  timer.assignTime("Trial, Na+ Occlusion Mesh (CPU)");
  BackgroundMesh<ullint> chloride_occm_cpu(GridDetail::OCCLUSION, 3.0,
                                           VdwCombiningRule::LORENTZ_BERTHELOT,
                                           tsm.getTopologyPointer(2), &chloride_cf,
                                           { -6.0, -6.0, -6.0, 7.5, 8.5, 9.1 }, 2.25);
  timer.assignTime("Trial, Cl- Occlusion Mesh (CPU)");
  BackgroundMesh<ullint> nacl_occm_cpu(GridDetail::OCCLUSION, 3.0,
                                       VdwCombiningRule::LORENTZ_BERTHELOT,
                                       tsm.getTopologyPointer(3), &nacl_cf,
                                       { -6.2, -6.0, -6.1, 7.5, 7.7, 8.4 }, 2.25);
  timer.assignTime("Trial, NaCl Occlusion Mesh (CPU)");
  checkBitSettings(sodium_occm_cpu, "CPU", tsm.getTestingStatus());
  checkBitSettings(chloride_occm_cpu, "CPU", tsm.getTestingStatus());
  checkBitSettings(nacl_occm_cpu, "CPU", tsm.getTestingStatus());
  timer.assignTime("Occlusion Mesh Checking");

  // Create occlusion meshes with non-orthogonal grid axes
  const double tetra = 109.4712206 * stormm::symbols::pi / 180.0;
  const std::vector<double> tric_box = { 1.2, 1.2, 1.2, tetra, tetra, tetra };
  std::vector<double> tric_umat(9), tric_invu(9);
  computeBoxTransform(tric_box, &tric_umat, &tric_invu);
  const MeshParameters triclinic_dims(24, 32, 36, -15.0, -15.0, -15.0, tric_invu, 36);
  BackgroundMesh<ullint> na_tric_occm_cpu(GridDetail::OCCLUSION, NonbondedPotential::CLASH,
                                            tsm.getTopologyPointer(1), &sodium_cf, 1.4, 0.0,
                                            VdwCombiningRule::LORENTZ_BERTHELOT, triclinic_dims);
  timer.assignTime("Trial, Na+ Tric. Occ. Mesh (CPU)");
  BackgroundMesh<ullint> nacl_tric_occm_cpu(GridDetail::OCCLUSION, NonbondedPotential::CLASH,
                                            tsm.getTopologyPointer(3), &nacl_cf, 1.4, 0.0,
                                            VdwCombiningRule::LORENTZ_BERTHELOT, triclinic_dims);
  timer.assignTime("Trial, NaCl Tric. Occ. Mesh (CPU)");
  checkBitSettings(na_tric_occm_cpu, "CPU", tsm.getTestingStatus());
  checkBitSettings(nacl_tric_occm_cpu, "CPU", tsm.getTestingStatus());
  timer.assignTime("Occlusion Mesh Checking");

  // Engage occlusion mesh GPU tests
#ifdef STORMM_USE_HPC
  trpi_cf.upload();
  sodium_cf.upload();
  chloride_cf.upload();
  nacl_cf.upload();

  // To avoid confusion in the timings, do a "dry run" to engage the first GPU kernel and toss
  // the result into "miscellaneous"
  if (oe.getDisplayTimingsOrder()) {
    BackgroundMesh<ullint> scratch_mesh(GridDetail::OCCLUSION, 3.0,
                                        VdwCombiningRule::LORENTZ_BERTHELOT,
                                        tsm.getTopologyPointer(1), &sodium_cf,
                                        { -5.0, -5.0, -5.0, 6.25, 6.25, 6.25 }, 1.125, 40, {},
                                        PrecisionModel::SINGLE, mesh_launcher,
                                        HybridTargetLevel::DEVICE);
  }
  timer.assignTime(0);
  BackgroundMesh<ullint> sodium_occm_gpu(GridDetail::OCCLUSION, 3.0,
                                         VdwCombiningRule::LORENTZ_BERTHELOT,
                                         tsm.getTopologyPointer(1), &sodium_cf,
                                         { -5.0, -5.0, -5.0, 6.25, 6.25, 6.25 }, 1.125, 40, {},
                                         PrecisionModel::SINGLE, mesh_launcher,
                                         HybridTargetLevel::DEVICE);
  timer.assignTime("Trial, Na+ Occlusion Mesh (GPU)");
  BackgroundMesh<ullint> chloride_occm_gpu(GridDetail::OCCLUSION, 3.0,
                                           VdwCombiningRule::LORENTZ_BERTHELOT,
                                           tsm.getTopologyPointer(2), &chloride_cf,
                                           { -6.0, -6.0, -6.0, 7.5, 8.5, 9.1 }, 2.25, 40, {},
                                           PrecisionModel::SINGLE, mesh_launcher,
                                           HybridTargetLevel::DEVICE);
  timer.assignTime("Trial, Cl- Occlusion Mesh (GPU)");
  BackgroundMesh<ullint> nacl_occm_gpu(GridDetail::OCCLUSION, 3.0,
                                       VdwCombiningRule::LORENTZ_BERTHELOT,
                                       tsm.getTopologyPointer(3), &nacl_cf,
                                       { -6.2, -6.0, -6.1, 7.5, 7.7, 8.4 }, 2.25, 40, {},
                                       PrecisionModel::SINGLE, mesh_launcher,
                                       HybridTargetLevel::DEVICE);
  timer.assignTime("Trial, NaCl Occlusion Mesh (GPU)");
  BackgroundMesh<ullint> trpi_occm_gpu(GridDetail::OCCLUSION, 3.0,
                                       VdwCombiningRule::LORENTZ_BERTHELOT,
                                       tsm.getTopologyPointer(0), &trpi_cf,
                                       { -20.8, -20.3, -20.5, 25.1, 25.7, 26.0 }, 3.0, 40, {},
                                       PrecisionModel::SINGLE, mesh_launcher,
                                       HybridTargetLevel::DEVICE);
  timer.assignTime("Trial, Trp-Cage Occlusion Mesh (GPU)");
  BackgroundMesh<ullint> nacl_tric_occm_gpu(GridDetail::OCCLUSION, NonbondedPotential::CLASH,
                                            tsm.getTopologyPointer(3), &nacl_cf, 1.4, 0.0,
                                            VdwCombiningRule::LORENTZ_BERTHELOT, triclinic_dims,
                                            {}, {}, 0.0, 0.0, PrecisionModel::SINGLE,
                                            mesh_launcher, HybridTargetLevel::DEVICE);
  timer.assignTime("Trial, NaCl Tric. Occ. Mesh (GPU)");
  BackgroundMesh<ullint> trpi_tric_occm_gpu(GridDetail::OCCLUSION, NonbondedPotential::CLASH,
                                            tsm.getTopologyPointer(0), &trpi_cf, 1.4, 0.0,
                                            VdwCombiningRule::LORENTZ_BERTHELOT, triclinic_dims,
                                            {}, {}, 0.0, 0.0, PrecisionModel::DOUBLE,
                                            mesh_launcher, HybridTargetLevel::DEVICE);
  timer.assignTime("Trial, Trp-Cage Tric. Occ. Mesh (GPU)");
  checkBitSettings(sodium_occm_gpu, "GPU", tsm.getTestingStatus());
  checkBitSettings(chloride_occm_gpu, "GPU", tsm.getTestingStatus());
  checkBitSettings(nacl_occm_gpu, "GPU", tsm.getTestingStatus());
  checkBitSettings(trpi_occm_gpu, "GPU", tsm.getTestingStatus());
  checkBitSettings(nacl_tric_occm_gpu, "GPU", tsm.getTestingStatus());
  checkBitSettings(trpi_tric_occm_gpu, "GPU", tsm.getTestingStatus());
  timer.assignTime("Occlusion Mesh Checking");

  // If timings were requested, perform some additional tests to compare CPU and GPU times
  if (oe.getDisplayTimingsOrder()) {
    const int trpi_cpu_mesh_timings = timer.addCategory("Trp-Cage Occlusion Mesh, CPU");
    const int trpi_gpu_mesh_timings = timer.addCategory("Trp-Cage Occlusion Mesh, GPU");
    timer.assignTime(0);
    BackgroundMesh<ullint> trpi_occm_cpu2(GridDetail::OCCLUSION, 3.0,
                                          VdwCombiningRule::LORENTZ_BERTHELOT,
                                          tsm.getTopologyPointer(0), &trpi_cf,
                                          { -20.0, -20.0, -20.0, 25.0, 25.0, 25.0 }, 1.5, 40);
    timer.assignTime(trpi_cpu_mesh_timings);
    BackgroundMesh<ullint> trpi_occm_gpu2(GridDetail::OCCLUSION, 3.0,
                                          VdwCombiningRule::LORENTZ_BERTHELOT,
                                          tsm.getTopologyPointer(0), &trpi_cf,
                                          { -20.0, -20.0, -20.0, 25.0, 25.0, 25.0 }, 1.5, 40, {},
                                          PrecisionModel::SINGLE, mesh_launcher,
                                          HybridTargetLevel::DEVICE);
    timer.assignTime(trpi_gpu_mesh_timings);
    CoordinateFrame tama_cf = tsm.exportCoordinateFrame(5);
    tama_cf.upload();
    const CoordinateFrameReader tama_cfr = tama_cf.data();
    const AtomGraph tama_ag = tsm.exportAtomGraph(5);
    MeshParameters field_mps(36, 36, 36, minValue(tama_cfr.xcrd, tama_cfr.natom) - 4.0,
                             minValue(tama_cfr.ycrd, tama_cfr.natom) - 4.0,
                             minValue(tama_cfr.zcrd, tama_cfr.natom) - 4.0, 1.0, 40,
                             Interpolant::SMOOTHNESS);
    timer.assignTime(0);
    const BackgroundMesh<float> tama_bgm_cpu(GridDetail::NONBONDED_FIELD,
                                             NonbondedPotential::ELECTROSTATIC,
                                             tama_ag, tama_cf, 3.1, 0.15,
                                             VdwCombiningRule::LORENTZ_BERTHELOT, field_mps, {},
                                             {}, 1.0, 1.0, PrecisionModel::SINGLE);
    timer.assignTime("Trial, Tamavidin Elec. Mesh (CPU)");
    const BackgroundMesh<float> tama_bgm_gpu(GridDetail::NONBONDED_FIELD,
                                             NonbondedPotential::ELECTROSTATIC,
                                             tama_ag, tama_cf, 3.1, 0.15,
                                             VdwCombiningRule::LORENTZ_BERTHELOT, field_mps, {},
                                             {}, 1.0, 1.0, PrecisionModel::SINGLE, mesh_launcher,
                                             HybridTargetLevel::DEVICE);
    timer.assignTime("Trial, Tamavidin Elec. Mesh (GPU)");
  }
#endif

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
