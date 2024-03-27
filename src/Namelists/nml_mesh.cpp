#include "copyright.h"
#include "Constants/symbol_values.h"
#include "Math/matrix_ops.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Reporting/error_format.h"
#include "Structure/mesh_parameters.h"
#include "input.h"
#include "nml_mesh.h"

namespace stormm {
namespace namelist {

using constants::getEnumerationName;
using energy::translateNonbondedPotential;
using energy::getEnumerationName;
using parse::minimalRealFormat;
using parse::NumberFormat;
using stmath::computeBoxTransform;
using structure::default_mesh_scaling_bits;
using structure::getEnumerationName;
using structure::translateBoundaryCondition;
using structure::translateGridDetail;
using structure::translateMeshPosition;
using structure::translateRMSDMethod;

//-------------------------------------------------------------------------------------------------
MeshControls::MeshControls(const ExceptionResponse policy_in) :
  policy{policy_in},
  mesh_points_a{default_mesh_grid_dim},
  mesh_points_b{default_mesh_grid_dim},
  mesh_points_c{default_mesh_grid_dim},
  mesh_spacing_a{default_mesh_grid_spacing},
  mesh_spacing_b{default_mesh_grid_spacing},
  mesh_spacing_c{default_mesh_grid_spacing},
  mesh_alpha{default_mesh_grid_angle},
  mesh_beta{default_mesh_grid_angle},
  mesh_gamma{default_mesh_grid_angle},
  mesh_origin_x{0.0}, mesh_origin_y{0.0}, mesh_origin_z{0.0},
  buffer_width_a{0.0}, buffer_width_b{0.0}, buffer_width_c{0.0},
  kind{std::string(default_mesh_grid_detail)},
  potential{std::string(default_mesh_grid_potential)},
  boundaries{std::string(default_mesh_grid_boundary)},
  mesh_scaling_bits{default_mesh_scaling_bits},
  clash_distance{default_mesh_elec_damping_range},
  clash_ratio{default_mesh_vdw_damping_ratio},
  nml_transcript{"mesh"}
{}

//-------------------------------------------------------------------------------------------------
MeshControls::MeshControls(const TextFile &tf, int *start_line, bool *found_nml,
                           const ExceptionResponse policy_in, const WrapTextSearch wrap) :
  MeshControls(policy_in)
{
  const NamelistEmulator t_nml = meshInput(tf, start_line, found_nml, policy, wrap);
  nml_transcript = t_nml;

  // Mesh dimensions are set to default values during the object's construction, but can be
  // replaced with other values if the associated keywords are specified in the user input.
  t_nml.assignVariable(&mesh_points_a, &mesh_points_b, &mesh_points_c, "mesh_dim");
  t_nml.assignVariable(&mesh_points_a, "mesh_dim_a");
  t_nml.assignVariable(&mesh_points_b, "mesh_dim_b");
  t_nml.assignVariable(&mesh_points_c, "mesh_dim_c");
  t_nml.assignVariable(&mesh_spacing_a, &mesh_spacing_b, &mesh_spacing_c, "mesh_spacing");
  t_nml.assignVariable(&mesh_spacing_a, "mesh_spacing_a");
  t_nml.assignVariable(&mesh_spacing_b, "mesh_spacing_b");
  t_nml.assignVariable(&mesh_spacing_c, "mesh_spacing_c");
  t_nml.assignVariable(&mesh_alpha, symbols::pi / 180.0, "mesh_alpha");
  t_nml.assignVariable(&mesh_beta,  symbols::pi / 180.0, "mesh_beta");
  t_nml.assignVariable(&mesh_gamma, symbols::pi / 180.0, "mesh_gamma");
  t_nml.assignVariable(&mesh_origin_x, "mesh_origin_x");
  t_nml.assignVariable(&mesh_origin_y, "mesh_origin_y");
  t_nml.assignVariable(&mesh_origin_z, "mesh_origin_z");
  t_nml.assignVariable(&buffer_width_a, &buffer_width_b, &buffer_width_c, "buffer_width");
  t_nml.assignVariable(&buffer_width_a, "buffer_width_a");
  t_nml.assignVariable(&buffer_width_b, "buffer_width_b");
  t_nml.assignVariable(&buffer_width_c, "buffer_width_c");
  setPotential(t_nml.getStringValue("potential"));
  setBoundaries(t_nml.getStringValue("boundary"));
  mesh_scaling_bits = t_nml.getIntValue("scaling_bits");
  validateMeshElement();
}

//-------------------------------------------------------------------------------------------------
int MeshControls::getAxisElementCount(const UnitCellAxis dim) const {
  switch (dim) {
  case UnitCellAxis::A:
    return mesh_points_a;
  case UnitCellAxis::B:
    return mesh_points_b;
  case UnitCellAxis::C:
    return mesh_points_c;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int MeshControls::getAxisElementCount(const CartesianDimension dim) const {
  return getAxisElementCount(static_cast<UnitCellAxis>(dim));
}

//-------------------------------------------------------------------------------------------------
double MeshControls::getSpacing(const UnitCellAxis dim) const {
  switch (dim) {
  case UnitCellAxis::A:
    return mesh_spacing_a;
  case UnitCellAxis::B:
    return mesh_spacing_b;
  case UnitCellAxis::C:
    return mesh_spacing_c;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double MeshControls::getSpacing(const CartesianDimension dim) const {
  return getSpacing(static_cast<UnitCellAxis>(dim));
}

//-------------------------------------------------------------------------------------------------
double MeshControls::getAlpha() const {
  return mesh_alpha;
}

//-------------------------------------------------------------------------------------------------
double MeshControls::getBeta() const {
  return mesh_beta;
}

//-------------------------------------------------------------------------------------------------
double MeshControls::getGamma() const {
  return mesh_gamma;
}

//-------------------------------------------------------------------------------------------------
double MeshControls::getOrigin(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return mesh_origin_x;
  case CartesianDimension::Y:
    return mesh_origin_y;
  case CartesianDimension::Z:
    return mesh_origin_z;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double MeshControls::getBufferWidth(const UnitCellAxis dim) const {
  switch (dim) {
  case UnitCellAxis::A:
    return buffer_width_a;
  case UnitCellAxis::B:
    return buffer_width_b;
  case UnitCellAxis::C:
    return buffer_width_c;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double MeshControls::getBufferWidth(const CartesianDimension dim) const {
  return getBufferWidth(static_cast<UnitCellAxis>(dim));
}

//-------------------------------------------------------------------------------------------------
GridDetail MeshControls::getDetail() const {
  return translateGridDetail(kind);
}

//-------------------------------------------------------------------------------------------------
NonbondedPotential MeshControls::getPotential() const {
  return translateNonbondedPotential(potential);
}

//-------------------------------------------------------------------------------------------------
BoundaryCondition MeshControls::getBoundaries() const {
  return translateBoundaryCondition(boundaries);
}

//-------------------------------------------------------------------------------------------------
int MeshControls::getScalingBits() const {
  return mesh_scaling_bits;
}

//-------------------------------------------------------------------------------------------------
double MeshControls::getElecClashDistance() const {
  return clash_distance;
}

//-------------------------------------------------------------------------------------------------
double MeshControls::getVdwClashRatio() const {
  return clash_ratio;
}

//-------------------------------------------------------------------------------------------------
const NamelistEmulator& MeshControls::getTranscript() const {
  return nml_transcript;
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setElementCount(const int mesh_points_in, const UnitCellAxis dim) {
  if (mesh_points_in <= 0) {
    rtErr("Invalid mesh dimension " + std::to_string(mesh_points_in) + " along the " +
          getEnumerationName(dim) + " axis.", "MeshControls", "setElementCount");
  }
  switch (dim) {
  case UnitCellAxis::A:
    mesh_points_a = mesh_points_in;
    break;
  case UnitCellAxis::B:
    mesh_points_b = mesh_points_in;
    break;
  case UnitCellAxis::C:
    mesh_points_c = mesh_points_in;
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setElementCount(const int mesh_points_in,
                                           const CartesianDimension dim) {
  setElementCount(mesh_points_in, static_cast<UnitCellAxis>(dim));
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setSpacing(const double mesh_spacing_in, const UnitCellAxis dim) {
  if (mesh_spacing_in <= 0.0) {
    rtErr("Invalid mesh spacing " + minimalRealFormat(mesh_spacing_in, 1.0e-4) + ".",
          "MeshControls", "setSpacing");
  }
  switch (dim) {
  case UnitCellAxis::A:
    mesh_spacing_a = mesh_spacing_in;
    break;
  case UnitCellAxis::B:
    mesh_spacing_b = mesh_spacing_in;
    break;
  case UnitCellAxis::C:
    mesh_spacing_c = mesh_spacing_in;
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setSpacing(const double mesh_spacing_in, const CartesianDimension dim) {
  setSpacing(mesh_spacing_in, static_cast<UnitCellAxis>(dim));
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setAlphaAngle(const double alpha_in) {
  mesh_alpha = alpha_in;
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setBetaAngle(const double beta_in) {
  mesh_beta = beta_in;
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setGammaAngle(const double gamma_in) {
  mesh_gamma = gamma_in;
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setOrigin(const double mesh_origin_in, const CartesianDimension dim) {
  switch (dim) {
  case CartesianDimension::X:
    mesh_origin_x = mesh_origin_in;
    break;
  case CartesianDimension::Y:
    mesh_origin_y = mesh_origin_in;
    break;
  case CartesianDimension::Z:
    mesh_origin_z = mesh_origin_in;
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setBufferWidth(const double buffer_width_in) {
  buffer_width_a = buffer_width_in;
  buffer_width_b = buffer_width_in;
  buffer_width_c = buffer_width_in;
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setBufferWidth(const double buffer_width_in, const UnitCellAxis dim) {
  switch (dim) {
  case UnitCellAxis::A:
    buffer_width_a = buffer_width_in;
  case UnitCellAxis::B:
    buffer_width_b = buffer_width_in;
  case UnitCellAxis::C:
    buffer_width_c = buffer_width_in;
  }
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setBufferWidth(const double buffer_width_in, const CartesianDimension dim) {
  setBufferWidth(buffer_width_in, static_cast<UnitCellAxis>(dim));
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setDetail(const std::string &kind_in) {
  kind = kind_in;
  try {
    const GridDetail interp = translateGridDetail(kind);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("An invalid mesh type " + kind + " was specified.", "MeshControls", "setDetail");
    case ExceptionResponse::WARN:
      rtWarn("An invalid mesh type " + kind + " was specified.  The default of " +
             std::string(default_mesh_grid_detail) + " will be restored.", "MeshControls",
             "setDetail");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    kind = std::string(default_mesh_grid_detail);
  }
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setDetail(const GridDetail kind_in) {
  kind = getEnumerationName(kind_in);
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setPotential(const std::string &potential_in) {
  potential = potential_in;
  try {
    const NonbondedPotential interp = translateNonbondedPotential(potential);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("An invalid potential form " + potential + " was specified.", "MeshControls",
            "setPotential");
    case ExceptionResponse::WARN:
      rtWarn("An invalid potential form " + potential + " was specified.  The default of " +
             std::string(default_mesh_grid_potential) + " will be restored.",
             "MeshControls", "setPotential");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    potential = std::string(default_mesh_grid_potential);
  }
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setPotential(const NonbondedPotential potential_in) {
  potential = getEnumerationName(potential_in);
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setBoundaries(const std::string &boundaries_in) {
  boundaries = boundaries_in;
  try {
    const BoundaryCondition interp = translateBoundaryCondition(boundaries);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("An invalid boundary condition " + boundaries + " was specified.", "MeshControls",
            "setBoundaries");
    case ExceptionResponse::WARN:
      rtWarn("An invalid boundary condition " + boundaries + " was specified.  The default of " +
             std::string(default_mesh_grid_boundary) + " will be restored.",
             "MeshControls", "setBoundaries");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    boundaries = std::string(default_mesh_grid_boundary);
  }
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setBoundaries(const BoundaryCondition boundaries_in) {
  boundaries = getEnumerationName(boundaries_in);
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setScalingBits(const int mesh_scaling_bits_in) {
  mesh_scaling_bits = mesh_scaling_bits_in;
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setElecClashDistance(const double clash_distance_in) {
  clash_distance = clash_distance_in;
}

//-------------------------------------------------------------------------------------------------
void MeshControls::setVdwClashRatio(const double clash_ratio_in) {
  clash_ratio = clash_ratio_in;
}

//-------------------------------------------------------------------------------------------------
void MeshControls::validateMeshElement() {

  // Try computing the transformation matrix into box space
  bool problem = false;
  std::vector<double> umat_test(9), invu_test(9);
  try {
    computeBoxTransform(mesh_spacing_a, mesh_spacing_b, mesh_spacing_c, mesh_alpha, mesh_beta,
                        mesh_gamma, umat_test.data(), invu_test.data());
  }
  catch (std::runtime_error) {
    problem = true;
  }

  // Check that the transformation matrices are valid.
  for (int i = 0; i < 9; i++) {
    problem = (problem || (std::isnan(umat_test[i]) || std::isinf(umat_test[i]) ||
                           std::isnan(invu_test[i]) || std::isinf(invu_test[i])));
  }
  if (problem) {
    const double rad_factor = 180.0 / symbols::pi;
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The mesh dimensions [ " +
            realToString(mesh_spacing_a, 9, 4, NumberFormat::STANDARD_REAL) + ", " +
            realToString(mesh_spacing_b, 9, 4, NumberFormat::STANDARD_REAL) + ", " +
            realToString(mesh_spacing_c, 9, 4, NumberFormat::STANDARD_REAL) +
            " (Angstroms) ] with angles [ " +
            realToString(mesh_alpha * rad_factor, 9, 4, NumberFormat::STANDARD_REAL) + ", " +
            realToString(mesh_beta * rad_factor, 9, 4, NumberFormat::STANDARD_REAL) + ", " +
            realToString(mesh_gamma * rad_factor, 9, 4, NumberFormat::STANDARD_REAL) +
            " (degrees) ] produce an error when attempting to compute transformation matrices.",
            "MeshControls", "validateMeshElement");
    case ExceptionResponse::WARN:
      rtWarn("The mesh dimensions [ " +
             realToString(mesh_spacing_a, 9, 4, NumberFormat::STANDARD_REAL) + ", " +
             realToString(mesh_spacing_b, 9, 4, NumberFormat::STANDARD_REAL) + ", " +
             realToString(mesh_spacing_c, 9, 4, NumberFormat::STANDARD_REAL) +
             " (Angstroms) ] with angles [ " +
             realToString(mesh_alpha * rad_factor, 9, 4, NumberFormat::STANDARD_REAL) + ", " +
             realToString(mesh_beta * rad_factor, 9, 4, NumberFormat::STANDARD_REAL) + ", " +
             realToString(mesh_gamma * rad_factor, 9, 4, NumberFormat::STANDARD_REAL) +
             " (degrees) ] produce an error when attempting to compute transformation matrices.  "
             "Any negative box edge lengths will be converted to their absolute values, and zero "
             "or extremely small box edge lengths will be converted to the default value of " +
             minimalRealFormat(default_mesh_grid_spacing, 1.0e-4) + " Angstroms.  Box angles will "
             "be iteratively moved back towards their default value of " +
             minimalRealFormat(default_mesh_grid_angle * rad_factor, 1.0e-4) + " degrees, until "
             "a viable matrix is produced.  Check the output for the outcome.  These inputs "
             "should be reviewed and corrected.", "MeshControls", "validateMeshElement");
      break;
    case ExceptionResponse::SILENT:
      break;
    }

    // Adjust box lengths
    mesh_spacing_a = fabs(mesh_spacing_a);
    mesh_spacing_b = fabs(mesh_spacing_b);
    mesh_spacing_c = fabs(mesh_spacing_c);
    if (mesh_spacing_a < 1.0e-6) {
      mesh_spacing_a =  default_mesh_grid_spacing;
    }
    if (mesh_spacing_b < 1.0e-6) {
      mesh_spacing_b =  default_mesh_grid_spacing;
    }
    if (mesh_spacing_c < 1.0e-6) {
      mesh_spacing_c =  default_mesh_grid_spacing;
    }

    // Adjust angles
    const double alpha_incr = 0.01 * (default_mesh_grid_angle - mesh_alpha);
    const double beta_incr  = 0.01 * (default_mesh_grid_angle - mesh_beta);
    const double gamma_incr = 0.01 * (default_mesh_grid_angle - mesh_gamma);
    int niter = 0;
    while (problem) {

      // The box length adjustment may have corrected the problem.  Test before incrementing any
      // of the angles.
      problem = false;
      try {
        computeBoxTransform(mesh_spacing_a, mesh_spacing_b, mesh_spacing_c, mesh_alpha, mesh_beta,
                        mesh_gamma, umat_test.data(), invu_test.data());
      }
      catch (std::runtime_error) {
        problem = true;
      }
      for (int i = 0; i < 9; i++) {
         problem = (problem || (std::isnan(umat_test[i]) || std::isinf(umat_test[i]) ||
                                std::isnan(invu_test[i]) || std::isinf(invu_test[i])));
      }
      if (problem) {
        mesh_alpha += alpha_incr;
        mesh_beta  += beta_incr;
        mesh_gamma += gamma_incr;
      }
      niter++;
      if (niter > 100) {
        rtErr("The mesh element could not be validated with the input dimensions.\n",
              "MeshControls", "validateMeshElement");
      }
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
NamelistEmulator meshInput(const TextFile &tf, int *start_line, bool *found,
                               const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("mesh", CaseSensitivity::AUTOMATIC, policy, "Collects details of a "
                         "mesh density field for docking and conformer selection in STORMM.");

  // Keyword: mesh_dim{_a, _b, _c}
  t_nml.addKeyword("mesh_dim_a", NamelistType::INTEGER);
  t_nml.addHelp("mesh_dim_a", "The number of mesh points to space at regular intervals along the "
                "Cartesian X axis (mesh a vector).  Negative values in this keyword imply that "
                "the dimension shall be computed at run time based on the system at hand.");
  t_nml.addKeyword("mesh_dim_b", NamelistType::INTEGER);
  t_nml.addHelp("mesh_dim_b", "The number of mesh points to space at regular intervals along the "
                "mesh b vector.  Negative values in this keyword imply that the dimension shall "
                "be computed at run time based on the system at hand.");
  t_nml.addKeyword("mesh_dim_c", NamelistType::INTEGER);
  t_nml.addHelp("mesh_dim_c", "The number of mesh points to space at regular intervals along the "
                "mesh c vector.  Negative values in this keyword imply that the dimension shall "
                "be computed at run time based on the system at hand.");
  t_nml.addKeyword("mesh_dim", NamelistType::INTEGER);
  t_nml.addHelp("mesh_dim", "General mesh dimension to apply along all three mesh axes.  "
                "Specifying a dimension for a particular axis will override this general "
                "directive.");

  // Keyword: mesh_spacing{_a, _b, _c}
  t_nml.addKeyword("mesh_spacing_a", NamelistType::REAL);
  t_nml.addHelp("mesh_spacing_a", "The regular spacing between mesh points along the mesh a "
                "vector, units of Angstroms.");
  t_nml.addKeyword("mesh_spacing_b", NamelistType::REAL);
  t_nml.addHelp("mesh_spacing_b", "The regular spacing between mesh points along the mesh b "
                "vector, units of Angstroms.");
  t_nml.addKeyword("mesh_spacing_c", NamelistType::REAL);
  t_nml.addHelp("mesh_spacing_c", "The regular spacing between mesh points along the mesh c "
                "vector, units of Angstroms.");
  t_nml.addKeyword("mesh_spacing", NamelistType::REAL);
  t_nml.addHelp("mesh_spacing", "General mesh spacing to apply along all three mesh axes.  "
                "Specifying the spacing for a particular axis will override this general "
                "directive.");

  // Keyword: mesh_alpha
  t_nml.addKeyword("mesh_alpha", NamelistType::REAL);
  t_nml.addHelp("mesh_alpha", "Angle between the mesh b and c vectors, in units of radians.");

  // Keyword: mesh_beta
  t_nml.addKeyword("mesh_beta", NamelistType::REAL);
  t_nml.addHelp("mesh_beta", "Angle between the mesh a and c vectors, in units of radians.");

  // Keyword: mesh_gamma
  t_nml.addKeyword("mesh_gamma", NamelistType::REAL);
  t_nml.addHelp("mesh_gamma", "Angle between the mesh a and b vectors, in units of radians.");

  // Keyword: mesh_origin{_x, _y, _z}
  t_nml.addKeyword("mesh_origin_x", NamelistType::REAL);
  t_nml.addHelp("mesh_origin_x", "Origin of the mesh along the Cartesian X axis, in units of "
                "Angstroms.");
  t_nml.addKeyword("mesh_origin_y", NamelistType::REAL);
  t_nml.addHelp("mesh_origin_y", "Origin of the mesh along the Cartesian Y axis, in units of "
                "Angstroms.");
  t_nml.addKeyword("mesh_origin_z", NamelistType::REAL);
  t_nml.addHelp("mesh_origin_z", "Origin of the mesh along the Cartesian Z axis, in units of "
                "Angstroms.");

  // Keyword: buffer_width{_a, _b, _c}
  t_nml.addKeyword("buffer_width_a", NamelistType::REAL);
  t_nml.addHelp("buffer_width_a", "The minimum distance between van-der Waals spheres of the "
                "mesh and the mesh boundary faces parallel to the b x c plane.");
  t_nml.addKeyword("buffer_width_b", NamelistType::REAL);
  t_nml.addHelp("buffer_width_b", "The minimum distance between van-der Waals spheres of the "
                "mesh and the mesh boundary faces parallel to the a x c plane.");
  t_nml.addKeyword("buffer_width_c", NamelistType::REAL);
  t_nml.addHelp("buffer_width_c", "The minimum distance between van-der Waals spheres of the "
                "mesh and the mesh boundary faces parallel to the a x b plane.");
  t_nml.addKeyword("buffer_width", NamelistType::REAL);
  t_nml.addHelp("buffer_width", "The minimum distance between van-der Waals spheres of the "
                "mesh and any mesh boundary faces.");
  
  // Keyword: potential
  t_nml.addKeyword("potential", NamelistType::STRING,
                   std::string(default_mesh_grid_potential));
  t_nml.addHelp("potential", "The type of non-bonded mesh potential to map onto the mesh (the "
                "mesh will express this potential in units of kcal/mol).");

  // Keyword: boundary
  t_nml.addKeyword("boundary", NamelistType::STRING,
                   std::string(default_mesh_grid_boundary));
  t_nml.addHelp("boundary", "The type of boundary conditions in which the mesh exists.  Applying "
                "periodic boundary conditions here implies that the mesh molecule will be "
                "treated as periodic in the given mesh dimensions.  Options include ISOLATED, "
                "NONPERIODIC, or NON_PERIODIC (isolated boundary conditions), and PERIODIC "
                "(periodic boundary conditions).");
  
  // Keyword: scaling_bits
  t_nml.addKeyword("scaling_bits", NamelistType::INTEGER,
                   std::to_string(default_mesh_scaling_bits));
  t_nml.addHelp("scaling_bits", "The number of bits after the decimal to use in setting mesh "
                "vertex positions as well as computing potential quantities.");

  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  return t_nml;
}
  
} // namespace namelist
} // namespace stormm
