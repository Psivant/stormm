#include "copyright.h"
#include "Numerics/split_fixed_precision.h"
#include "Parsing/parsing_enumerators.h"
#include "Reporting/error_format.h"
#include "mesh_rulers.h"

namespace stormm {
namespace structure {

using card::HybridKind;
using numerics::fixedPrecisionGrid;
using numerics::hostInt95Sum;
using numerics::hostSplitFPSum;
using numerics::hostInt95ToDouble;
using parse::NumberFormat;

//-------------------------------------------------------------------------------------------------
MeshRulerKit::MeshRulerKit(const llint* avec_x_in, const llint* avec_y_in, const llint* avec_z_in,
                           const llint* bvec_x_in, const llint* bvec_y_in, const llint* bvec_z_in,
                           const llint* cvec_x_in, const llint* cvec_y_in, const llint* cvec_z_in,
                           const llint* avec_abs_x_in, const llint* avec_abs_y_in,
                           const llint* avec_abs_z_in, const int* avec_x_ovrf_in,
                           const int* avec_y_ovrf_in, const int* avec_z_ovrf_in,
                           const int* bvec_x_ovrf_in, const int* bvec_y_ovrf_in,
                           const int* bvec_z_ovrf_in, const int* cvec_x_ovrf_in,
                           const int* cvec_y_ovrf_in, const int* cvec_z_ovrf_in,
                           const int* avec_abs_x_ovrf_in, const int* avec_abs_y_ovrf_in,
                           const int* avec_abs_z_ovrf_in) :
  avec_x{avec_x_in}, avec_y{avec_y_in}, avec_z{avec_z_in}, bvec_x{bvec_x_in}, bvec_y{bvec_y_in},
  bvec_z{bvec_z_in}, cvec_x{cvec_x_in}, cvec_y{cvec_y_in}, cvec_z{cvec_z_in},
  avec_abs_x{avec_abs_x_in}, avec_abs_y{avec_abs_y_in}, avec_abs_z{avec_abs_z_in},
  avec_x_ovrf{avec_x_ovrf_in}, avec_y_ovrf{avec_y_ovrf_in}, avec_z_ovrf{avec_z_ovrf_in},
  bvec_x_ovrf{bvec_x_ovrf_in}, bvec_y_ovrf{bvec_y_ovrf_in}, bvec_z_ovrf{bvec_z_ovrf_in},
  cvec_x_ovrf{cvec_x_ovrf_in}, cvec_y_ovrf{cvec_y_ovrf_in}, cvec_z_ovrf{cvec_z_ovrf_in},
  avec_abs_x_ovrf{avec_abs_x_ovrf_in}, avec_abs_y_ovrf{avec_abs_y_ovrf_in},
  avec_abs_z_ovrf{avec_abs_z_ovrf_in}
{}
  
//-------------------------------------------------------------------------------------------------
MeshRulers::MeshRulers(const MeshParameters &mps) :
    a_line_x{HybridKind::POINTER, "mesh_avector_x"},
    a_line_y{HybridKind::POINTER, "mesh_avector_y"},
    a_line_z{HybridKind::POINTER, "mesh_avector_z"},
    b_line_x{HybridKind::POINTER, "mesh_bvector_x"},
    b_line_y{HybridKind::POINTER, "mesh_bvector_y"},
    b_line_z{HybridKind::POINTER, "mesh_bvector_z"},
    c_line_x{HybridKind::POINTER, "mesh_cvector_x"},
    c_line_y{HybridKind::POINTER, "mesh_cvector_y"},
    c_line_z{HybridKind::POINTER, "mesh_cvector_z"},
    a_abs_line_x{HybridKind::POINTER, "mesh_avector_abs_x"},
    a_abs_line_y{HybridKind::POINTER, "mesh_avector_abs_y"},
    a_abs_line_z{HybridKind::POINTER, "mesh_avector_abs_z"},
    a_line_x_overflow{HybridKind::POINTER, "mesh_avec_x_ovrf"},
    a_line_y_overflow{HybridKind::POINTER, "mesh_avec_y_ovrf"},
    a_line_z_overflow{HybridKind::POINTER, "mesh_avec_z_ovrf"},
    b_line_x_overflow{HybridKind::POINTER, "mesh_bvec_x_ovrf"},
    b_line_y_overflow{HybridKind::POINTER, "mesh_bvec_y_ovrf"},
    b_line_z_overflow{HybridKind::POINTER, "mesh_bvec_z_ovrf"},
    c_line_x_overflow{HybridKind::POINTER, "mesh_cvec_x_ovrf"},
    c_line_y_overflow{HybridKind::POINTER, "mesh_cvec_y_ovrf"},
    c_line_z_overflow{HybridKind::POINTER, "mesh_cvec_z_ovrf"},
    a_abs_line_x_overflow{HybridKind::POINTER, "mesh_avec_abs_x_ovrf"},
    a_abs_line_y_overflow{HybridKind::POINTER, "mesh_avec_abs_y_ovrf"},
    a_abs_line_z_overflow{HybridKind::POINTER, "mesh_avec_abs_z_ovrf"},
    int_data{HybridKind::ARRAY, "mesh_int_data"},
    llint_data{HybridKind::ARRAY, "mesh_llint_data"},
    mps_pointer{const_cast<MeshParameters*>(mps.getSelfPointer())}
{
  // Allocate the proper amount of memory and set POINTER-kind Hybrid objects.
  const MeshParamKit mpsk = mps.data();
  const int padded_na = roundUp(mpsk.na + 1, warp_size_int);
  const int padded_nb = roundUp(mpsk.nb + 1, warp_size_int);
  const int padded_nc = roundUp(mpsk.nc + 1, warp_size_int);
  llint_data.resize(3 * ((2 * padded_na) + padded_nb + padded_nc));
  int_data.resize(3 * ((2 * padded_na) + padded_nb + padded_nc));
  a_line_x.setPointer(&llint_data,             0, mpsk.na + 1);
  a_line_y.setPointer(&llint_data,     padded_na, mpsk.na + 1);
  a_line_z.setPointer(&llint_data, 2 * padded_na, mpsk.na + 1);
  a_line_x_overflow.setPointer(&int_data,             0, mpsk.na + 1);
  a_line_y_overflow.setPointer(&int_data,     padded_na, mpsk.na + 1);
  a_line_z_overflow.setPointer(&int_data, 2 * padded_na, mpsk.na + 1);
  int thus_far = 3 * padded_na;
  b_line_x.setPointer(&llint_data,                   thus_far, mpsk.nb + 1);
  b_line_y.setPointer(&llint_data,       padded_nb + thus_far, mpsk.nb + 1);
  b_line_z.setPointer(&llint_data, (2 * padded_nb) + thus_far, mpsk.nb + 1);
  b_line_x_overflow.setPointer(&int_data,                   thus_far, mpsk.nb + 1);
  b_line_y_overflow.setPointer(&int_data,       padded_nb + thus_far, mpsk.nb + 1);
  b_line_z_overflow.setPointer(&int_data, (2 * padded_nb) + thus_far, mpsk.nb + 1);
  thus_far += 3 * padded_nb;
  c_line_x.setPointer(&llint_data,                   thus_far, mpsk.nc + 1);
  c_line_y.setPointer(&llint_data,       padded_nc + thus_far, mpsk.nc + 1);
  c_line_z.setPointer(&llint_data, (2 * padded_nc) + thus_far, mpsk.nc + 1);
  c_line_x_overflow.setPointer(&int_data,                   thus_far, mpsk.nc + 1);
  c_line_y_overflow.setPointer(&int_data,       padded_nc + thus_far, mpsk.nc + 1);
  c_line_z_overflow.setPointer(&int_data, (2 * padded_nc) + thus_far, mpsk.nc + 1);
  thus_far += 3 * padded_nc;
  a_abs_line_x.setPointer(&llint_data,                   thus_far, mpsk.na + 1);
  a_abs_line_y.setPointer(&llint_data,       padded_na + thus_far, mpsk.na + 1);
  a_abs_line_z.setPointer(&llint_data, (2 * padded_na) + thus_far, mpsk.na + 1);
  a_abs_line_x_overflow.setPointer(&int_data,                   thus_far, mpsk.na + 1);
  a_abs_line_y_overflow.setPointer(&int_data,       padded_na + thus_far, mpsk.na + 1);
  a_abs_line_z_overflow.setPointer(&int_data, (2 * padded_na) + thus_far, mpsk.na + 1);
  
  // Map out the axes of the mesh, drawing three lines radiating from the mesh origin along its
  // "a", "b", and "c" axes and marking the Cartesian coordinates of tick marks for the origins of
  // each grid point in lines [ i_a, 0, 0 ], [ 0, j_b, 0 ], or [ 0, 0, k_c ] on these axes.  All
  // locations are given relative to the origin, such that adding the coordinates of the origin to
  // the coordinates of point i_a from the "a" axis line, the coordinates of point j_b from the
  // "b" axis line, and point k_c from the "c" axis line will give the location of the origin for
  // grid element { i_a, j_b, k_c }.  A fourth axis, the "a" axis points with the origin
  // coordinates "baked in", is provided for convenience.  Given that that the majority of
  // transactions with the mesh will involve a whole warp computing the 64-term tricubic polynomial
  // for one particle after just one thread computes the proper element and broadcasts the relative
  // displacments within that element, this savings may be significant.
  fixedPrecisionGrid(&a_line_x, &a_line_x_overflow, { 0LL, 0 }, mpsk.fp_invu[0]);
  fixedPrecisionGrid(&a_line_y, &a_line_y_overflow, { 0LL, 0 }, mpsk.fp_invu[1]);
  fixedPrecisionGrid(&a_line_z, &a_line_z_overflow, { 0LL, 0 }, mpsk.fp_invu[2]);
  fixedPrecisionGrid(&b_line_x, &b_line_x_overflow, { 0LL, 0 }, mpsk.fp_invu[3]);
  fixedPrecisionGrid(&b_line_y, &b_line_y_overflow, { 0LL, 0 }, mpsk.fp_invu[4]);
  fixedPrecisionGrid(&b_line_z, &b_line_z_overflow, { 0LL, 0 }, mpsk.fp_invu[5]);
  fixedPrecisionGrid(&c_line_x, &c_line_x_overflow, { 0LL, 0 }, mpsk.fp_invu[6]);
  fixedPrecisionGrid(&c_line_y, &c_line_y_overflow, { 0LL, 0 }, mpsk.fp_invu[7]);
  fixedPrecisionGrid(&c_line_z, &c_line_z_overflow, { 0LL, 0 }, mpsk.fp_invu[8]);
  fixedPrecisionGrid(&a_abs_line_x, &a_abs_line_x_overflow, mpsk.orig_x, mpsk.fp_invu[0]);
  fixedPrecisionGrid(&a_abs_line_y, &a_abs_line_y_overflow, mpsk.orig_y, mpsk.fp_invu[1]);
  fixedPrecisionGrid(&a_abs_line_z, &a_abs_line_z_overflow, mpsk.orig_z, mpsk.fp_invu[2]);
}

//-------------------------------------------------------------------------------------------------
MeshRulers::MeshRulers(const MeshRulers &original) :
  a_line_x{original.a_line_x},
  a_line_y{original.a_line_y},
  a_line_z{original.a_line_z},
  b_line_x{original.b_line_x},
  b_line_y{original.b_line_y},
  b_line_z{original.b_line_z},
  c_line_x{original.c_line_x},
  c_line_y{original.c_line_y},
  c_line_z{original.c_line_z},
  a_abs_line_x{original.a_abs_line_x},
  a_abs_line_y{original.a_abs_line_y},
  a_abs_line_z{original.a_abs_line_z},
  a_line_x_overflow{original.a_line_x_overflow},
  a_line_y_overflow{original.a_line_y_overflow},
  a_line_z_overflow{original.a_line_z_overflow},
  b_line_x_overflow{original.b_line_x_overflow},
  b_line_y_overflow{original.b_line_y_overflow},
  b_line_z_overflow{original.b_line_z_overflow},
  c_line_x_overflow{original.c_line_x_overflow},
  c_line_y_overflow{original.c_line_y_overflow},
  c_line_z_overflow{original.c_line_z_overflow},
  a_abs_line_x_overflow{original.a_abs_line_x_overflow},
  a_abs_line_y_overflow{original.a_abs_line_y_overflow},
  a_abs_line_z_overflow{original.a_abs_line_z_overflow},
  int_data{original.int_data},
  llint_data{original.llint_data}
{
  // Repair the object's pointers
  rebasePointers();
}

//-------------------------------------------------------------------------------------------------
MeshRulers::MeshRulers(MeshRulers &&original) :
  a_line_x{std::move(original.a_line_x)},
  a_line_y{std::move(original.a_line_y)},
  a_line_z{std::move(original.a_line_z)},
  b_line_x{std::move(original.b_line_x)},
  b_line_y{std::move(original.b_line_y)},
  b_line_z{std::move(original.b_line_z)},
  c_line_x{std::move(original.c_line_x)},
  c_line_y{std::move(original.c_line_y)},
  c_line_z{std::move(original.c_line_z)},
  a_abs_line_x{std::move(original.a_abs_line_x)},
  a_abs_line_y{std::move(original.a_abs_line_y)},
  a_abs_line_z{std::move(original.a_abs_line_z)},
  a_line_x_overflow{std::move(original.a_line_x_overflow)},
  a_line_y_overflow{std::move(original.a_line_y_overflow)},
  a_line_z_overflow{std::move(original.a_line_z_overflow)},
  b_line_x_overflow{std::move(original.b_line_x_overflow)},
  b_line_y_overflow{std::move(original.b_line_y_overflow)},
  b_line_z_overflow{std::move(original.b_line_z_overflow)},
  c_line_x_overflow{std::move(original.c_line_x_overflow)},
  c_line_y_overflow{std::move(original.c_line_y_overflow)},
  c_line_z_overflow{std::move(original.c_line_z_overflow)},
  a_abs_line_x_overflow{std::move(original.a_abs_line_x_overflow)},
  a_abs_line_y_overflow{std::move(original.a_abs_line_y_overflow)},
  a_abs_line_z_overflow{std::move(original.a_abs_line_z_overflow)},
  int_data{std::move(original.int_data)},
  llint_data{std::move(original.llint_data)}
{}

//-------------------------------------------------------------------------------------------------
MeshRulers& MeshRulers::operator=(const MeshRulers &other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  a_line_x = other.a_line_x;
  a_line_y = other.a_line_y;
  a_line_z = other.a_line_z;
  b_line_x = other.b_line_x;
  b_line_y = other.b_line_y;
  b_line_z = other.b_line_z;
  c_line_x = other.c_line_x;
  c_line_y = other.c_line_y;
  c_line_z = other.c_line_z;
  a_abs_line_x = other.a_abs_line_x;
  a_abs_line_y = other.a_abs_line_y;
  a_abs_line_z = other.a_abs_line_z;
  a_line_x_overflow = other.a_line_x_overflow;
  a_line_y_overflow = other.a_line_y_overflow;
  a_line_z_overflow = other.a_line_z_overflow;
  b_line_x_overflow = other.b_line_x_overflow;
  b_line_y_overflow = other.b_line_y_overflow;
  b_line_z_overflow = other.b_line_z_overflow;
  c_line_x_overflow = other.c_line_x_overflow;
  c_line_y_overflow = other.c_line_y_overflow;
  c_line_z_overflow = other.c_line_z_overflow;
  a_abs_line_x_overflow = other.a_abs_line_x_overflow;
  a_abs_line_y_overflow = other.a_abs_line_y_overflow;
  a_abs_line_z_overflow = other.a_abs_line_z_overflow;
  int_data = other.int_data;
  llint_data = other.llint_data;

  // Repair the object's pointers
  rebasePointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
MeshRulers& MeshRulers::operator=(MeshRulers &&other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  a_line_x = std::move(other.a_line_x);
  a_line_y = std::move(other.a_line_y);
  a_line_z = std::move(other.a_line_z);
  b_line_x = std::move(other.b_line_x);
  b_line_y = std::move(other.b_line_y);
  b_line_z = std::move(other.b_line_z);
  c_line_x = std::move(other.c_line_x);
  c_line_y = std::move(other.c_line_y);
  c_line_z = std::move(other.c_line_z);
  a_abs_line_x = std::move(other.a_abs_line_x);
  a_abs_line_y = std::move(other.a_abs_line_y);
  a_abs_line_z = std::move(other.a_abs_line_z);
  a_line_x_overflow = std::move(other.a_line_x_overflow);
  a_line_y_overflow = std::move(other.a_line_y_overflow);
  a_line_z_overflow = std::move(other.a_line_z_overflow);
  b_line_x_overflow = std::move(other.b_line_x_overflow);
  b_line_y_overflow = std::move(other.b_line_y_overflow);
  b_line_z_overflow = std::move(other.b_line_z_overflow);
  c_line_x_overflow = std::move(other.c_line_x_overflow);
  c_line_y_overflow = std::move(other.c_line_y_overflow);
  c_line_z_overflow = std::move(other.c_line_z_overflow);
  a_abs_line_x_overflow = std::move(other.a_abs_line_x_overflow);
  a_abs_line_y_overflow = std::move(other.a_abs_line_y_overflow);
  a_abs_line_z_overflow = std::move(other.a_abs_line_z_overflow);
  int_data = std::move(other.int_data);
  llint_data = std::move(other.llint_data);

  // No pointer repair is needed as everything was moved.
  return *this;
}

//-------------------------------------------------------------------------------------------------
double3 MeshRulers::getMeshOrigin() const {
  const double inv_scale = mps_pointer->getInverseScalingFactor();
  return { hostInt95ToDouble(a_abs_line_x.readHost(0), a_abs_line_x_overflow.readHost(0)) *
           inv_scale,
           hostInt95ToDouble(a_abs_line_y.readHost(0), a_abs_line_y_overflow.readHost(0)) *
           inv_scale,
           hostInt95ToDouble(a_abs_line_z.readHost(0), a_abs_line_z_overflow.readHost(0)) *
           inv_scale };
}

//-------------------------------------------------------------------------------------------------
double3 MeshRulers::getRealLocation(const double3 mesh_loc) const {
  double image_x = mesh_loc.x;
  double image_y = mesh_loc.y;
  double image_z = mesh_loc.z;
  const MeshParamKit mpsk = mps_pointer->data();
  const double d_na = mpsk.na;
  const double d_nb = mpsk.nb;
  const double d_nc = mpsk.nc;
  switch (mpsk.bounds) {
  case BoundaryCondition::ISOLATED:
    if (image_x < 0.0 || image_x >= d_na || image_y < 0.0 || image_y >= d_nb ||
        image_z < 0.0 || image_z >= d_nc) {
      rtErr("A point at " + realToString(image_x, 9, 4, NumberFormat::STANDARD_REAL) + ", " +
            realToString(image_y, 9, 4, NumberFormat::STANDARD_REAL) + ", " +
            realToString(image_z, 9, 4, NumberFormat::STANDARD_REAL) + " is off the mesh ("
            "dimensions " + std::to_string(mpsk.na) + ", " + std::to_string(mpsk.nb) + ", " +
            std::to_string(mpsk.nc) + ").", "MeshRulers", "getRealLocation");
    }
    break;
  case BoundaryCondition::PERIODIC:
    image_x -= floor(image_x / d_na) * d_na;
    image_y -= floor(image_y / d_nb) * d_nb;
    image_z -= floor(image_z / d_nc) * d_nc;
    break;
  }
  int mesh_a = floor(image_x);
  int mesh_b = floor(image_y);
  int mesh_c = floor(image_z);
  const double mesh_da = image_x - static_cast<double>(mesh_a);
  const double mesh_db = image_y - static_cast<double>(mesh_b);
  const double mesh_dc = image_z - static_cast<double>(mesh_c);
  int95_t ielem_x = hostInt95Sum(a_abs_line_x.readHost(mesh_a),
                                 a_abs_line_x_overflow.readHost(mesh_a),
                                 b_line_x.readHost(mesh_b), b_line_x_overflow.readHost(mesh_b));
  int95_t ielem_y = hostInt95Sum(a_abs_line_y.readHost(mesh_a),
                                 a_abs_line_y_overflow.readHost(mesh_a),
                                 b_line_y.readHost(mesh_b), b_line_y_overflow.readHost(mesh_b));
  int95_t ielem_z = hostInt95Sum(a_abs_line_z.readHost(mesh_a),
                                 a_abs_line_z_overflow.readHost(mesh_a),
                                 b_line_z.readHost(mesh_b), b_line_z_overflow.readHost(mesh_b));
  ielem_x = hostSplitFPSum(ielem_x, c_line_x.readHost(mesh_c), c_line_x_overflow.readHost(mesh_c));
  ielem_y = hostSplitFPSum(ielem_y, c_line_y.readHost(mesh_c), c_line_y_overflow.readHost(mesh_c));
  ielem_z = hostSplitFPSum(ielem_z, c_line_z.readHost(mesh_c), c_line_z_overflow.readHost(mesh_c));
  const double inv_scale = mpsk.inv_scale;
  const double dx = (hostInt95ToDouble(ielem_x) * mpsk.inv_scale) +
                    (mpsk.invu[0] * mesh_da) + (mpsk.invu[3] * mesh_db) + (mpsk.invu[6] * mesh_dc);
  const double dy = (hostInt95ToDouble(ielem_y) * mpsk.inv_scale) +
                                               (mpsk.invu[4] * mesh_db) + (mpsk.invu[7] * mesh_dc);
  const double dz = (hostInt95ToDouble(ielem_z) * mpsk.inv_scale) + (mpsk.invu[8] * mesh_dc);
  return { dx, dy, dz };
}

//-------------------------------------------------------------------------------------------------
double3 MeshRulers::getMeshLocation(const double3 real_loc) const {
  double image_x = real_loc.x;
  double image_y = real_loc.y;
  double image_z = real_loc.z;
  const MeshParamKit mpsk = mps_pointer->data();
  const double d_na = mpsk.na;
  const double d_nb = mpsk.nb;
  const double d_nc = mpsk.nc;
  double ma = (mpsk.umat[0] * image_x) + (mpsk.umat[3] * image_y) + (mpsk.umat[6] * image_z);
  double mb =                            (mpsk.umat[4] * image_y) + (mpsk.umat[7] * image_z);
  double mc =                                                       (mpsk.umat[8] * image_z);
  switch (mpsk.bounds) {
  case BoundaryCondition::ISOLATED:
    if (ma < 0.0 || ma >= d_na || mb < 0.0 || mb > d_nb || mc < 0.0 || mc > d_nc) {
      rtErr("A point at " + realToString(image_x, 9, 4, NumberFormat::STANDARD_REAL) + ", " +
            realToString(image_y, 9, 4, NumberFormat::STANDARD_REAL) + ", " +
            realToString(image_z, 9, 4, NumberFormat::STANDARD_REAL) + " is off the mesh ("
            "dimensions " + std::to_string(mpsk.na) + ", " + std::to_string(mpsk.nb) + ", " +
            std::to_string(mpsk.nc) + ").", "MeshRulers", "getMeshLocation");
    }
    break;
  case BoundaryCondition::PERIODIC:
    ma -= floor(ma / d_na) * d_na;
    mb -= floor(mb / d_nb) * d_nb;
    mc -= floor(mc / d_nc) * d_nc;
    break;
  }
  return { ma, mb, mc };
}

//-------------------------------------------------------------------------------------------------
const MeshRulerKit MeshRulers::data(const HybridTargetLevel tier) const {
  return MeshRulerKit(a_line_x.data(tier), a_line_y.data(tier), a_line_z.data(tier),
                      b_line_x.data(tier), b_line_y.data(tier), b_line_z.data(tier),
                      c_line_x.data(tier), c_line_y.data(tier), c_line_z.data(tier),
                      a_abs_line_x.data(tier), a_abs_line_y.data(tier), a_abs_line_z.data(tier),
                      a_line_x_overflow.data(tier), a_line_y_overflow.data(tier),
                      a_line_z_overflow.data(tier), b_line_x_overflow.data(tier),
                      b_line_y_overflow.data(tier), b_line_z_overflow.data(tier),
                      c_line_x_overflow.data(tier), c_line_y_overflow.data(tier),
                      c_line_z_overflow.data(tier), a_abs_line_x_overflow.data(tier),
                      a_abs_line_y_overflow.data(tier), a_abs_line_z_overflow.data(tier));
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void MeshRulers::upload() {
  int_data.upload();
  llint_data.upload();
}

//-------------------------------------------------------------------------------------------------
void MeshRulers::download() {
  int_data.download();
  llint_data.download();
}

#  ifdef STORMM_USE_CUDA
//-------------------------------------------------------------------------------------------------
const MeshRulerKit MeshRulers::deviceViewToHostData() const {
  return MeshRulerKit(a_line_x.getDeviceValidHostPointer(), a_line_y.getDeviceValidHostPointer(),
                      a_line_z.getDeviceValidHostPointer(), b_line_x.getDeviceValidHostPointer(),
                      b_line_y.getDeviceValidHostPointer(), b_line_z.getDeviceValidHostPointer(),
                      c_line_x.getDeviceValidHostPointer(), c_line_y.getDeviceValidHostPointer(),
                      c_line_z.getDeviceValidHostPointer(),
                      a_abs_line_x.getDeviceValidHostPointer(),
                      a_abs_line_y.getDeviceValidHostPointer(),
                      a_abs_line_z.getDeviceValidHostPointer(),
                      a_line_x_overflow.getDeviceValidHostPointer(),
                      a_line_y_overflow.getDeviceValidHostPointer(),
                      a_line_z_overflow.getDeviceValidHostPointer(),
                      b_line_x_overflow.getDeviceValidHostPointer(),
                      b_line_y_overflow.getDeviceValidHostPointer(),
                      b_line_z_overflow.getDeviceValidHostPointer(),
                      c_line_x_overflow.getDeviceValidHostPointer(),
                      c_line_y_overflow.getDeviceValidHostPointer(),
                      c_line_z_overflow.getDeviceValidHostPointer(),
                      a_abs_line_x_overflow.getDeviceValidHostPointer(),
                      a_abs_line_y_overflow.getDeviceValidHostPointer(),
                      a_abs_line_z_overflow.getDeviceValidHostPointer());
}
#  endif
#endif
//-------------------------------------------------------------------------------------------------
void MeshRulers::rebasePointers() {
  a_line_x.swapTarget(&llint_data);
  a_line_y.swapTarget(&llint_data);
  a_line_z.swapTarget(&llint_data);
  b_line_x.swapTarget(&llint_data);
  b_line_y.swapTarget(&llint_data);
  b_line_z.swapTarget(&llint_data);
  c_line_x.swapTarget(&llint_data);
  c_line_y.swapTarget(&llint_data);
  c_line_z.swapTarget(&llint_data);
  a_abs_line_x.swapTarget(&llint_data);
  a_abs_line_y.swapTarget(&llint_data);
  a_abs_line_z.swapTarget(&llint_data);
  a_line_x_overflow.swapTarget(&int_data);
  a_line_y_overflow.swapTarget(&int_data);
  a_line_z_overflow.swapTarget(&int_data);
  b_line_x_overflow.swapTarget(&int_data);
  b_line_y_overflow.swapTarget(&int_data);
  b_line_z_overflow.swapTarget(&int_data);
  c_line_x_overflow.swapTarget(&int_data);
  c_line_y_overflow.swapTarget(&int_data);
  c_line_z_overflow.swapTarget(&int_data);
  a_abs_line_x_overflow.swapTarget(&int_data);
  a_abs_line_y_overflow.swapTarget(&int_data);
  a_abs_line_z_overflow.swapTarget(&int_data);
}
  
} // namespace structure
} // namespace stormm
