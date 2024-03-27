#include <climits>
#include "copyright.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "Constants/fixed_precision.h"
#include "Math/matrix_ops.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "mesh_parameters.h"

namespace stormm {
namespace structure {

using numerics::min_localpos_scale_bits;
using parse::realToString;
using parse::NumberFormat;
using stmath::extractBoxDimensions;
using stmath::hessianNormalWidths;
using stmath::invertSquareMatrix;

//-------------------------------------------------------------------------------------------------
MeshParamKit::MeshParamKit(const int na_in, const int nb_in, const int nc_in,
                           const int95_t orig_x_in, const int95_t orig_y_in,
                           const int95_t orig_z_in, const double scale_in,
                           const double inv_scale_in, const int scale_bits_in,
                           const double* umat_in, const double* invu_in,
                           const double* full_umat_in, const double* full_invu_in,
                           const double* widths_in, const int95_t* fp_invu_in,
                           const double max_span_in, const BoundaryCondition bounds_in,
                           const Interpolant stencil_kind_in, const UnitCellType unit_cell_in) :
    na{na_in}, nb{nb_in}, nc{nc_in}, orig_x{orig_x_in}, orig_y{orig_y_in}, orig_z{orig_z_in},
    scale{scale_in}, inv_scale{inv_scale_in}, scale_f{static_cast<float>(scale_in)},
    inv_scale_f{static_cast<float>(inv_scale_in)}, scale_bits{scale_bits_in},
    umat{ umat_in[0], umat_in[1], umat_in[2], umat_in[3], umat_in[4], umat_in[5], umat_in[6],
          umat_in[7], umat_in[8] },
    invu{ invu_in[0], invu_in[1], invu_in[2], invu_in[3], invu_in[4], invu_in[5], invu_in[6],
          invu_in[7], invu_in[8] },
    full_umat{ full_umat_in[0], full_umat_in[1], full_umat_in[2], full_umat_in[3], full_umat_in[4],
               full_umat_in[5], full_umat_in[6], full_umat_in[7], full_umat_in[8] },
    full_invu{ full_invu_in[0], full_invu_in[1], full_invu_in[2], full_invu_in[3], full_invu_in[4],
               full_invu_in[5], full_invu_in[6], full_invu_in[7], full_invu_in[8] },
    widths{ widths_in[0], widths_in[1], widths_in[2] },
    fp_invu{ fp_invu_in[0], fp_invu_in[1], fp_invu_in[2], fp_invu_in[3], fp_invu_in[4],
             fp_invu_in[5], fp_invu_in[6], fp_invu_in[7], fp_invu_in[8] },
    max_span{max_span_in}, bounds{bounds_in}, stencil_kind{stencil_kind_in},
    unit_cell{unit_cell_in}
{}

//-------------------------------------------------------------------------------------------------
MeshParameters::MeshParameters(const int na_in, const int nb_in, const int nc_in,
                               const double origin_x_in, const double origin_y_in,
                               const double origin_z_in,
                               const std::vector<double> &element_vectors,
                               const int scale_bits_in, const Interpolant stencil_kind_in) :
    na{na_in}, nb{nb_in}, nc{nc_in}, origin_x{0LL, 0}, origin_y{0LL, 0}, origin_z{0LL, 0},
    scale_bits{scale_bits_in}, scale_factor{pow(2.0, scale_bits)},
    inverse_scale_factor{1.0 / scale_factor}, unit_cell{UnitCellType::NONE},
    boundary{BoundaryCondition::ISOLATED}, element_umat{},
    stencil_kind{stencil_kind_in},
    sp_element_umat{}, element_invu{}, sp_element_invu{}, full_umat{}, full_invu{}, widths{},
    fp_element_invu{}
{
  origin_x = hostDoubleToInt95(origin_x_in * scale_factor);
  origin_y = hostDoubleToInt95(origin_y_in * scale_factor);
  origin_z = hostDoubleToInt95(origin_z_in * scale_factor);
  defineElement(element_vectors);
  validateMeshDimensions();
  validateFixedPrecisionBits();
  maximum_span = maximumSpan();
}

//-------------------------------------------------------------------------------------------------
MeshParameters::MeshParameters() :
    MeshParameters(1, 1, 1, 0.0, 0.0, 0.0, { 1.0, 1.0, 1.0 })
{}

//-------------------------------------------------------------------------------------------------
MeshParameters::MeshParameters(const int na_in, const int nb_in, const int nc_in,
                               const double origin_x_in, const double origin_y_in,
                               const double origin_z_in, const double element_x,
                               const double element_y, const double element_z,
                               const int scale_bits_in, const Interpolant stencil_kind_in) :
    MeshParameters(na_in, nb_in, nc_in, origin_x_in, origin_y_in, origin_z_in,
                   { element_x, element_y, element_z }, scale_bits_in, stencil_kind_in)
{}

//-------------------------------------------------------------------------------------------------
MeshParameters::MeshParameters(const int na_in, const int nb_in, const int nc_in,
                               const double origin_x_in, const double origin_y_in,
                               const double origin_z_in, const double element_width,
                               const int scale_bits_in, const Interpolant stencil_kind_in) :
    MeshParameters(na_in, nb_in, nc_in, origin_x_in, origin_y_in, origin_z_in,
                   { element_width, element_width, element_width }, scale_bits_in, stencil_kind_in)
{}

//-------------------------------------------------------------------------------------------------
int MeshParameters::getAxisElementCount(const UnitCellAxis dim) const {
  switch (dim) {
  case UnitCellAxis::A:
    return na;
  case UnitCellAxis::B:
    return nb;
  case UnitCellAxis::C:
    return nc;
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
int MeshParameters::getAxisElementCount(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return na;
  case CartesianDimension::Y:
    return nb;
  case CartesianDimension::Z:
    return nc;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double MeshParameters::getMeshOrigin(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return hostInt95ToDouble(origin_x) * inverse_scale_factor;
  case CartesianDimension::Y:
    return hostInt95ToDouble(origin_y) * inverse_scale_factor;
  case CartesianDimension::Z:
    return hostInt95ToDouble(origin_z) * inverse_scale_factor;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int95_t> MeshParameters::getMeshOriginAsFP() const {
  std::vector<int95_t> result(3);
  result[0] = origin_x;
  result[1] = origin_y;
  result[2] = origin_z;
  return result;
}

//-------------------------------------------------------------------------------------------------
int95_t MeshParameters::getMeshOriginAsFP(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return origin_x;
  case CartesianDimension::Y:
    return origin_y;
  case CartesianDimension::Z:
    return origin_z;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
UnitCellType MeshParameters::getMeshCellType() const {
  return unit_cell;
}

//-------------------------------------------------------------------------------------------------
BoundaryCondition MeshParameters::getBoundaryConditions() const {
  return boundary;
}

//-------------------------------------------------------------------------------------------------
Interpolant MeshParameters::getStencilKind() const {
  return stencil_kind;
}

//-------------------------------------------------------------------------------------------------
std::vector<int95_t> MeshParameters::getMeshElementVectorAsFP(const UnitCellAxis dim) const {
  std::vector<int95_t> result(3);
  const int icol = static_cast<int>(dim);
  for (int i = 0; i < 3; i++) {
    result[i] = fp_element_invu[i + (3 * icol)];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int95_t> MeshParameters::getMeshInverseTransformAsFP() const {
  std::vector<int95_t> result(9);
  for (int i = 0; i < 9; i++) {
    result[i] = fp_element_invu[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int MeshParameters::getScalingBits() const {
  return scale_bits;
}

//-------------------------------------------------------------------------------------------------
double MeshParameters::getScalingFactor() const {
  return scale_factor;
}

//-------------------------------------------------------------------------------------------------
double MeshParameters::getInverseScalingFactor() const {
  return inverse_scale_factor;
}

//-------------------------------------------------------------------------------------------------
std::vector<int95_t> MeshParameters::getAxisCoordinates(const UnitCellAxis mesh_axis,
                                                        const CartesianDimension cart_axis) const {
  std::vector<int95_t> result;
  const int cdim = static_cast<int>(cart_axis);
  switch (mesh_axis) {
  case UnitCellAxis::A:
    result.resize(na + 1);
    result[0] = origin_x;
    for (int i = 0; i < na; i++) {
      result[i + 1] = hostSplitFPSum(result[i], fp_element_invu[cdim]);
    }
    break;
  case UnitCellAxis::B:
    result.resize(nb + 1);
    result[0] = origin_y;
    for (int i = 0; i < na; i++) {
      result[i + 1] = hostSplitFPSum(result[i], fp_element_invu[3 + cdim]);
    }
    break;
  case UnitCellAxis::C:
    result.resize(nc + 1);
    result[0] = origin_z;
    for (int i = 0; i < nc; i++) {
      result[i + 1] = hostSplitFPSum(result[i], fp_element_invu[6 + cdim]);
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double MeshParameters::getMaximumSpan() const {
  return maximum_span;
}

//-------------------------------------------------------------------------------------------------
std::string MeshParameters::printDimensions() const {
  const NumberFormat nfmt = NumberFormat::STANDARD_REAL;
  double xlen, ylen, zlen, alpha, beta, gamma;
  extractBoxDimensions(&xlen, &ylen, &zlen, &alpha, &beta, &gamma, element_invu);
  const double pifac = 180.0 / symbols::pi;
  const std::string result("Mesh dimensions: [ lengths " + realToString(xlen, 8, 5, nfmt) +
                           " x " + realToString(ylen, 8, 5, nfmt) + " x " +
                           realToString(zlen, 8, 5, nfmt) + ", angles " +
                           realToString(alpha * pifac, 8, 5, nfmt) + " x " +
                           realToString(beta * pifac, 8, 5, nfmt) + " x " +
                           realToString(gamma * pifac, 8, 5, nfmt) + "].");
  return result;
}

//-------------------------------------------------------------------------------------------------
MeshParamKit MeshParameters::data() const {
  return MeshParamKit(na, nb, nc, origin_x, origin_y, origin_z, scale_factor, inverse_scale_factor,
                      scale_bits, element_umat, element_invu, full_umat, full_invu, widths,
                      fp_element_invu, maximum_span, boundary, stencil_kind, unit_cell);
}

//-------------------------------------------------------------------------------------------------
const MeshParameters* MeshParameters::getSelfPointer() const  {
  return this;
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setMeshDimension(const int n_in, const UnitCellAxis mesh_axis) {
  if (n_in < 0) {
    rtErr("A mesh dimension of " + std::to_string(n_in) + " along the " +
          getEnumerationName(mesh_axis) + " is invalid.", "MeshParameters", "setMeshDimension");
  }
  switch (mesh_axis) {
  case UnitCellAxis::A:
    na = n_in;
  case UnitCellAxis::B:
    nb = n_in;
  case UnitCellAxis::C:
    nc = n_in;
  }
  validateMeshDimensions();
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setMeshDimension(const std::vector<int> &n_in) {
  if (n_in.size() != 3LLU) {
    rtErr("A vector of three elements is required (" + std::to_string(n_in.size()) + " provided).",
          "MeshParameters", "setMeshDimension");
  }
  na = n_in[0];
  nb = n_in[1];
  nc = n_in[2];
  validateMeshDimensions();  
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setOrigin(const double v, CartesianDimension cart_axis) {
  setOrigin(hostDoubleToInt95(v), cart_axis);
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setOrigin(const int95_t v, CartesianDimension cart_axis) {
  switch (cart_axis) {
  case CartesianDimension::X:
    origin_x = v;
    break;
  case CartesianDimension::Y:
    origin_y = v;
    break;
  case CartesianDimension::Z:
    origin_z = v;
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setOrigin(const std::vector<double> &v) {
  if (v.size() != 3LLU) {
    rtErr("A vector of three elements is required (" + std::to_string(v.size()) + " provided).",
          "MeshParameters", "setOrigin");
  }
  std::vector<int95_t> vfp = { hostDoubleToInt95(v[0]), hostDoubleToInt95(v[1]),
                               hostDoubleToInt95(v[2]) };
  setOrigin(vfp);
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setOrigin(const std::vector<int95_t> &v) {
  if (v.size() != 3LLU) {
    rtErr("A vector of three elements is required (" + std::to_string(v.size()) + " provided).",
          "MeshParameters", "setOrigin");
  }
  origin_x = v[0];
  origin_y = v[1];
  origin_z = v[2];
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setScalingBits(const int scale_bits_in) {
  scale_bits = scale_bits_in;
  validateFixedPrecisionBits();
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setBoundaryCondition(const BoundaryCondition boundary_in) {
  boundary = boundary_in;
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::setStencilKind(const Interpolant stencil_kind_in) {
  stencil_kind = stencil_kind_in;
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::defineElement(const std::vector<double> &element_vectors) {
  if (element_vectors.size() == 9LLU) {
    for (int i = 0; i < 9; i++) {
      element_invu[i] = element_vectors[i];
      sp_element_invu[i] = element_vectors[i];
      fp_element_invu[i] = hostDoubleToInt95(element_vectors[i] * scale_factor);
    }

    // A copy of the box vectors must be created so that the inverse transformation matrix, which
    // is given, does not itself get inverted while computing the transformation into box space.
    std::vector<double> copy_invu(9);
    for (int i = 0; i < 9; i++) {
      copy_invu[i] = element_invu[i];
    }
    invertSquareMatrix(copy_invu.data(), element_umat, 3);
    for (int i = 0; i < 9; i++) {
      sp_element_umat[i] = element_umat[i];
    }
  }
  else if (element_vectors.size() == 3LLU) {
    for (int i = 0; i < 9; i++) {
      element_umat[i] = 0.0;
      sp_element_umat[i] = 0.0;
      element_invu[i] = 0.0;
      sp_element_invu[i] = 0.0;
      fp_element_invu[i] = { 0LL, 0 };
    }
    for (int i = 0; i < 3; i++) {
      element_umat[4 * i] = 1.0 / element_vectors[i];
      sp_element_umat[4 * i] = element_umat[4 * i];
      element_invu[4 * i] = element_vectors[i];
      sp_element_invu[4 * i] = element_invu[4 * i];
      fp_element_invu[4 * i] = hostDoubleToInt95(element_vectors[i] * scale_factor);
    }
  }
  else {
    rtErr("The mesh element is defined by a 3x3 matrix.  A total of " +
          std::to_string(element_vectors.size()) + " elements were provided.", "MeshParameters");
  }

  // Define the unit cell for the entire mesh
  const double nvals[3] = { static_cast<double>(na), static_cast<double>(nb),
                            static_cast<double>(nc) };
  for (int i = 0; i < 3; i++) {
    for (int j = i; j < 3; j++) {
      full_invu[(j * 3) + i] = element_invu[(j * 3) + i] * nvals[j];
      full_umat[(j * 3) + i] = element_umat[(j * 3) + i] / nvals[i];
    }
  }
  
  // The determineUnitCellTypeByShape() function will define what would be very small unit cells,
  // in particular 1 x 1 x 1 Angstrom, as "NONE" type.  This is because there must be a unit cell
  // defined for some purposes, and such a unit cell would be an impossibly small simulation.
  // However, for mesh elements, the default "NONE" type unit cell dimensions are quite common.
  // Use a basic off-diagonal check instead.
  if (fabs(element_invu[1]) > constants::tiny || fabs(element_invu[2]) > constants::tiny ||
      fabs(element_invu[3]) > constants::tiny || fabs(element_invu[5]) > constants::tiny ||
      fabs(element_invu[6]) > constants::tiny || fabs(element_invu[7]) > constants::tiny) {
    unit_cell = UnitCellType::TRICLINIC;
  }
  else {
    unit_cell = UnitCellType::ORTHORHOMBIC;
  }

  // Use the Hessian Normal form to compute the number of mesh elements to search in each direction
  // and color all of the necessary elements around each atom.
  hessianNormalWidths(element_invu, &widths[0], &widths[1], &widths[2]);
  for (int i = 0; i < 3; i++) {
    sp_widths[i] = widths[i];
  }
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::validateMeshDimensions() const {
  const llint total_elements = static_cast<llint>(na) * static_cast<llint>(nb) *
                               static_cast<llint>(nc);
  if (total_elements > INT_MAX) {
    rtErr("The total number of elements on the mesh cannot exceed " + std::to_string(INT_MAX) +
          " (currently " + std::to_string(total_elements) + ").", "MeshParameters",
          "validateMeshDimensions");
  }
}

//-------------------------------------------------------------------------------------------------
void MeshParameters::validateFixedPrecisionBits() const {
  if (scale_bits < min_localpos_scale_bits) {
    rtErr("A minimum of " + std::to_string(min_localpos_scale_bits) + " must be stored after the "
          "decimal in fixed-precision mesh representations (" + std::to_string(scale_bits) +
          " specified).", "MeshParameters", "validateFixedPrecisionBits");
  }
  if (scale_bits > max_mesh_definition_bits) {
    rtErr("A maximum of " + std::to_string(max_mesh_definition_bits) + " is permitted for "
          "defining the locations of mesh vertex elements (" + std::to_string(scale_bits) +
          " specified).", "MeshParameters", "validateFixedPrecisionBits");
  }
}

//-------------------------------------------------------------------------------------------------
double MeshParameters::maximumSpan() const {
  double max_span = 0.0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {

      // Take the (i, j) vertex on the lower face parallel to the XY plane and compare it to the
      // caddy-corner vertex on the upper face parallel to the XY plane.
      const int ci = 1 - i;
      const int cj = 1 - j;
      const double di  = i;
      const double dci = ci;
      const double dj  = j;
      const double dcj = cj;
      const double x_ij = (di * element_invu[0]) + (dj * element_invu[3]);
      const double y_ij = (di * element_invu[1]) + (dj * element_invu[4]);
      const double z_ij = (di * element_invu[2]) + (dj * element_invu[5]);
      const double x_cij = (dci * element_invu[0]) + (dcj * element_invu[3]) + element_invu[6];
      const double y_cij = (dci * element_invu[1]) + (dcj * element_invu[4]) + element_invu[7];
      const double z_cij = (dci * element_invu[2]) + (dcj * element_invu[5]) + element_invu[8];
      const double dx = x_cij - x_ij;
      const double dy = y_cij - y_ij;
      const double dz = z_cij - z_ij;
      max_span = std::max(max_span, (dx * dx) + (dy * dy) + (dz * dz));
    }
  }
  return sqrt(max_span);
}

//-------------------------------------------------------------------------------------------------
MeshParameters getMeasurements(const AtomGraph *ag, const CoordinateFrame *cf,
                               const double padding, const double spacing,
                               const int scale_bits_in) {
  return getMeasurements(ag, cf, padding, std::vector<double>(3, spacing), scale_bits_in);
}

//-------------------------------------------------------------------------------------------------
MeshParameters getMeasurements(const AtomGraph *ag, const CoordinateFrame *cf,
                               const std::vector<double> &mesh_bounds, const double spacing,
                               const int scale_bits_in) {
  return getMeasurements(mesh_bounds, std::vector<double>(3, spacing), scale_bits_in);
}

//-------------------------------------------------------------------------------------------------
MeshParameters getMeasurements(const AtomGraph *ag, const CoordinateFrame *cf,
                               const double padding, const std::vector<double> &spacing,
                               const int scale_bits_in) {
  if (cf == nullptr || ag == nullptr) {
    rtErr("When creating a mesh based on a zone surrounding a system of interest, the topology "
          "and coordinates must be defined.  Supply the mesh boundaries explicitly to avoid "
          "relying on pre-defined coordinates, or submit the system to the mesh prior to calling "
          "this function.", "getMeasurements");
  }
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const CoordinateFrameReader cfr = cf->data();
  double xmin, ymin, zmin, xmax, ymax, zmax;
  bool points_unset = true;
  const int lj_idx_offset = nbk.n_lj_types + 1;
  if (spacing.size() == 3) {

    // Build an orthorhombic mesh.
    for (int pos = 0; pos < nbk.natom; pos++) {
      if (ag->getAtomMobility(pos)) {
        continue;
      }
      const size_t plj_idx = lj_idx_offset * nbk.lj_idx[pos];
      const double atom_radius = 0.5 * pow(nbk.lja_coeff[plj_idx] / nbk.ljb_coeff[plj_idx],
                                           1.0 / 6.0);
      if (points_unset) {
        xmin = cfr.xcrd[pos] - atom_radius;
        xmax = cfr.xcrd[pos] + atom_radius;
        ymin = cfr.ycrd[pos] - atom_radius;
        ymax = cfr.ycrd[pos] + atom_radius;
        zmin = cfr.zcrd[pos] - atom_radius;
        zmax = cfr.zcrd[pos] + atom_radius;
        points_unset = false;
      }
      else {
        xmin = std::min(xmin, cfr.xcrd[pos] - atom_radius);
        xmax = std::max(xmax, cfr.xcrd[pos] + atom_radius);
        ymin = std::min(ymin, cfr.ycrd[pos] - atom_radius);
        ymax = std::max(ymax, cfr.ycrd[pos] + atom_radius);
        zmin = std::min(zmin, cfr.zcrd[pos] - atom_radius);
        zmax = std::max(zmax, cfr.zcrd[pos] + atom_radius);
      }
    }
    xmin -= padding;
    xmax += padding;
    ymin -= padding;
    ymax += padding;
    zmin -= padding;
    zmax += padding;
    const std::vector<double> limits = { xmin, ymin, zmin, xmax, ymax, zmax };
    return getMeasurements(limits, spacing, scale_bits_in);
  }
  else if (spacing.size() == 9) {

    // Determine the required dimensions and appropriate origin of a triclinic mesh to enclose the
    // system by the requested padding buffer.
    const std::vector<double> widths = hessianNormalWidths<double>(spacing);
    std::vector<double> umat(9);
    invertSquareMatrix(spacing.data(), umat.data(), 3);
    for (int pos = 0; pos < nbk.natom; pos++) {
      if (ag->getAtomMobility(pos)) {
        continue;
      }
      const size_t plj_idx = lj_idx_offset * nbk.lj_idx[pos];
      const double atom_radius = 0.5 * pow(nbk.lja_coeff[plj_idx] / nbk.ljb_coeff[plj_idx],
                                           1.0 / 6.0);
      const double radfx = atom_radius / widths[0];
      const double radfy = atom_radius / widths[1];
      const double radfz = atom_radius / widths[2];
      const double xfrac = (umat[0] * cfr.xcrd[pos]) + (umat[3] * cfr.ycrd[pos]) +
                           (umat[6] * cfr.zcrd[pos]);
      const double yfrac = (umat[1] * cfr.xcrd[pos]) + (umat[4] * cfr.ycrd[pos]) +
                           (umat[7] * cfr.zcrd[pos]);
      const double zfrac = (umat[2] * cfr.xcrd[pos]) + (umat[5] * cfr.ycrd[pos]) +
                           (umat[8] * cfr.zcrd[pos]);
      if (points_unset) {
        xmin = xfrac - radfx;
        xmax = xfrac + radfx;
        ymin = yfrac - radfy;
        ymax = yfrac + radfy;
        zmin = zfrac - radfz;
        zmax = zfrac + radfz;
        points_unset = false;
      }
      else {
        xmin = std::min(xmin, xfrac - radfx);
        xmax = std::max(xmax, xfrac + radfx);
        ymin = std::min(ymin, yfrac - radfy);
        ymax = std::max(ymax, yfrac + radfy);
        zmin = std::min(zmin, zfrac - radfz);
        zmax = std::max(zmax, zfrac + radfz);
      }
    }
    const double padfx = padding / widths[0];
    const double padfy = padding / widths[1];
    const double padfz = padding / widths[2];
    xmin -= padfx;
    xmax += padfx;
    ymin -= padfy;
    ymax += padfy;
    zmin -= padfz;
    zmax += padfz;
    const double xdim = ceil(xmax - xmin);
    const double ydim = ceil(ymax - ymin);
    const double zdim = ceil(zmax - zmin);
    const double x_overhang = 0.5 * (xdim - (xmax - xmin));
    const double y_overhang = 0.5 * (ydim - (ymax - ymin));
    const double z_overhang = 0.5 * (zdim - (zmax - zmin));
    xmin -= x_overhang;
    ymin -= y_overhang;
    zmin -= z_overhang;
    const double xorig = (spacing[0] * xmin) + (spacing[3] * ymin) + (spacing[6] * zmin);
    const double yorig = (spacing[1] * xmin) + (spacing[4] * ymin) + (spacing[7] * zmin);
    const double zorig = (spacing[2] * xmin) + (spacing[5] * ymin) + (spacing[8] * zmin);
    const std::vector<double> limits = { xorig, yorig, zorig, xdim, ydim, zdim };
    return getMeasurements(limits, spacing, scale_bits_in);
  }
  else {
    rtErr("The mesh element dimensions must contain three orthogonal lengths or a full matrix "
          "with nine elements.", "getMeasurements");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
MeshParameters getMeasurements(const std::vector<double> &mesh_bounds,
                               const std::vector<double> &spacing, const int scale_bits_in) {
  if (mesh_bounds.size() != 6) {
    rtErr("An array of six elements, the minimum X, Y, and Z Cartesian coordinates followed by "
          "the maximum coordinates, is required.  " + std::to_string(mesh_bounds.size()) +
          " elements were provided.", "getMeasurements");
  }
  if (spacing.size() != 3 && spacing.size() != 9) {
    rtErr("An array of three elements, the length, width, and height (Cartesian X, Y, and Z "
          "dimensions) of a rectilinear mesh element, or nine elements defining the bounding "
          "vectors of a triclininc element, is required.  " + std::to_string(spacing.size()) +
          " elements were provided.", "getMeasurements");
  }
  if (spacing.size() == 3) {

    // In the orthorhombic case, assume that the mesh bounds represent minimum Cartesian X, Y,
    // and Z coordiantes followed by maximum X, Y, and Z coordinates.
    const std::vector<double> mesh_limits = {
      std::min(mesh_bounds[0], mesh_bounds[3]), std::min(mesh_bounds[1], mesh_bounds[4]),
      std::min(mesh_bounds[2], mesh_bounds[5]), std::max(mesh_bounds[0], mesh_bounds[3]),
      std::max(mesh_bounds[1], mesh_bounds[4]), std::max(mesh_bounds[2], mesh_bounds[5]) };
    const int pna = ceil((mesh_limits[3] - mesh_limits[0]) / spacing[0]);
    const int pnb = ceil((mesh_limits[4] - mesh_limits[1]) / spacing[1]);
    const int pnc = ceil((mesh_limits[5] - mesh_limits[2]) / spacing[2]);
    const double dnx = static_cast<double>(pna) * spacing[0];
    const double dny = static_cast<double>(pnb) * spacing[1];
    const double dnz = static_cast<double>(pnc) * spacing[2];
    const double overshoot_x = 0.5 * (dnx - (mesh_limits[3] - mesh_limits[0]));
    const double overshoot_y = 0.5 * (dny - (mesh_limits[4] - mesh_limits[1]));
    const double overshoot_z = 0.5 * (dnz - (mesh_limits[5] - mesh_limits[2]));
    return MeshParameters(pna, pnb, pnc, mesh_limits[0] - overshoot_x,
                          mesh_limits[1] - overshoot_y, mesh_limits[2] - overshoot_z, spacing,
                          scale_bits_in);
  }
  else {

    // In the non-orthorhombic case, assume that the mesh bounds represent the minimum Cartesian
    // X, Y, and Z coordinates (the origin) followed by mesh element counts along the a, b, and c
    // unit cell axes.
    const int pna = round(mesh_bounds[3]);
    const int pnb = round(mesh_bounds[4]);
    const int pnc = round(mesh_bounds[5]);
    if (fabs(static_cast<double>(pna) - mesh_bounds[3]) > 1.0e-4 ||
        fabs(static_cast<double>(pnb) - mesh_bounds[4]) > 1.0e-4 ||
        fabs(static_cast<double>(pnc) - mesh_bounds[5]) > 1.0e-4) {
      rtErr("A non-orthorhombic cell must have integral values for the number of elements along "
            "its a, b, and c axes, as expressed in the final three elements of the bounds array.  "
            "The values " + realToString(mesh_bounds[3], 9, 4, NumberFormat::STANDARD_REAL) +
            ", " + realToString(mesh_bounds[4], 9, 4, NumberFormat::STANDARD_REAL) + ", " +
            realToString(mesh_bounds[5], 9, 4, NumberFormat::STANDARD_REAL) + " are invalid.",
            "getMeasurements");
    }
    return MeshParameters(pna, pnb, pnc, mesh_bounds[0], mesh_bounds[1], mesh_bounds[2], spacing,
                          scale_bits_in);
  }
  __builtin_unreachable();
}

} // namespace structure
} // namespace stormm
