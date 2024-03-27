// -*-c++-*-
#ifndef STORMM_PUREMESH_H
#define STORMM_PUREMESH_H

#include <algorithm>
#include <string>
#include <vector>
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/mesh_kernel_manager.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/matrix_ops.h"
#include "Math/one_dimensional_splines.h"
#include "Math/radial_derivatives.h"
#include "Math/rounding.h"
#include "Math/tricubic_cell.h"
#include "Math/vector_ops.h"
#include "Namelists/nml_mesh.h"
#include "Numerics/split_fixed_precision.h"
#include "Parsing/polynumeric.h"
#include "Potential/energy_enumerators.h"
#include "Reporting/reporting_enumerators.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_copy.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "local_arrangement.h"
#include "mesh_forcefield.h"
#include "mesh_foundation.h"
#include "mesh_parameters.h"
#include "mesh_rulers.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

using card::GpuDetails;
using card::Hybrid;
using card::HybridKind;
using card::HybridTargetLevel;
using card::MeshKlManager;
using constants::CartesianDimension;
using constants::PrecisionModel;
using constants::UnitCellAxis;
using data_types::getStormmScalarTypeName;
using data_types::isScalarType;
using data_types::isFloatingPointScalarType;
using data_types::isSignedIntegralScalarType;
using energy::NonbondedPotential;
using energy::cubicSoftCore;
using energy::quinticSoftCore;
using energy::VdwCombiningRule;
using parse::NumberFormat;
using namelist::default_mesh_density_averaging_order;
using namelist::default_mesh_elec_damping_range;
using namelist::default_mesh_vdw_damping_ratio;
using numerics::default_globalpos_scale_bits;
using numerics::fixedPrecisionGrid;
using review::GridFileSyntax;
using review::getEnumerationName;
using stmath::addScalarToVector;
using stmath::convertData;
using stmath::elementwiseMultiply;
using stmath::evaluateCubicSpline;
using stmath::evaluateQuinticSpline;
using stmath::fvStencilCoordinates;
using stmath::hessianNormalWidths;
using stmath::invertSquareMatrix;
using stmath::maxAbsoluteDifference;
using stmath::radialFirstDerivative;
using stmath::radialSecondDerivative;
using stmath::radialThirdDerivative;
using stmath::roundUp;
using stmath::TricubicCell;
using stmath::TricubicStencil;
using topology::AtomGraph;
using topology::NonbondedKit;
using trajectory::coordCopy;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
using trajectory::PhaseSpace;

/// \brief The templated, writeable abstract of a BackgroundMesh object.
template <typename Tdata> struct BackgroundMeshWriter {

  /// \brief The constructor takes arguments for all member variables.
  BackgroundMeshWriter(const MeshParamKit &measurements, GridDetail kind_in,
                       NonbondedPotential field_in, const MeshRulerKit &rulers,
                       Tdata* coeffs_in, double coeff_scale_in, double probe_radius_in,
                       double well_depth_in, double occ_cost_in, const MeshBasicsKit &mbss_in);

  /// \brief The default copy and move constructors will be valid for this object.  Const members
  ///        negate the use of default copy and move assignment operators.
  ///
  /// \param original  The object to copy or move
  /// \{
  BackgroundMeshWriter(const BackgroundMeshWriter<Tdata> &original) = default;
  BackgroundMeshWriter(BackgroundMeshWriter<Tdata> &&original) = default;
  /// \}

  const MeshParamKit dims;         ///< Dimensions of the mesh.  These are pre-established.
  const GridDetail kind;           ///< The type of mesh, also pre-established.
  const NonbondedPotential field;  ///< The field described by the mesh, also pre-established.
  const MeshRulerKit rlrs;         ///< Coordinate series for tick marks along each mesh axis

  // The content of the mesh is what is actually writeable.
  Tdata* coeffs;              ///< Coefficients for all mesh elements.  In an OCCLUSION mesh, these
                              ///<   are the bit-packed masks for each cubelet, 64 cubelets making
                              ///<   one element.  In a NONBONDED_FIELD or NONBONDED_ATOMIC mesh,
                              ///<   the coefficients are tricubic splines for each element.
  const double coeff_scale;   ///< Scaling factor to take mesh coefficients from internal units
                              ///<   (kcal/mol, kcal/mol-A, kcal/mol-A2, and kcal/mol-A3, all based
                              ///<   on the same factor) into a fixed-precision representation.
                              ///<   this conversion is only performed when accumulating meshes.
                              ///<   Meshes used in actual energy calculations will be assumed to
                              ///<   be of real scalar types (float, double).
  const float coeff_scale_f;  ///< Single-precision variant of coeff_scale

  // The following data is critical to producing the mesh, and can be necessary for interpreting
  // the potentials read from it.
  const double probe_radius;    ///< The probe radius to use when mapping a Lennard-Jones mesh with
                                ///<   geometric combining rules, or when mapping a clash potential
                                ///<   to the mesh
  const double well_depth;      ///< The Lennard-Jones well depth for a test particle used to map a
                                ///<   Lennard-Jones potential of any mixing rules to the mesh.
  const double occ_cost;        ///< Energetic penalty of an occlusion interaction
  const MeshBasicsKit mbss;     ///< Collection of pointers for essential elements about the
                                ///<   molecular system underlying the mesh, including a neighbor
                                ///<   list for each cell (if relevant) and a mask of frozen atoms
};

/// \brief The templated, writeable abstract of a BackgroundMesh object.
template <typename Tdata> struct BackgroundMeshReader {

  /// \brief The constructor takes arguments for all member variables.
  BackgroundMeshReader(const MeshParamKit &dims_in, GridDetail kind_in,
                       NonbondedPotential field_in, const MeshRulerKit &rulers,
                       const Tdata* coeffs_in, double coeff_scale_in, double probe_radius_in,
                       double well_depth_in, double occ_cost_in, const MeshBasicsKit &mbss_in);

  /// \brief The default copy and move constructors will be valid for this object.  Const members
  ///        negate the use of default copy and move assignment operators.
  ///
  /// \param original  The object to copy or move
  /// \{
  BackgroundMeshReader(const BackgroundMeshReader<Tdata> &original) = default;
  BackgroundMeshReader(BackgroundMeshReader<Tdata> &&original) = default;
  /// \}

  const MeshParamKit dims;         ///< Dimensions of the mesh.  These are pre-established.
  const GridDetail kind;           ///< The type of mesh, also pre-established.
  const NonbondedPotential field;  ///< The field described by the mesh, also pre-established.
  const MeshRulerKit rlrs;         ///< Coordinate series for tick marks along each mesh axis 
  const Tdata* coeffs;             ///< Coefficients for all mesh elements.  In an OCCLUSION mesh,
                                   ///<   these are the bit-packed masks for each cubelet, 64
                                   ///<   cubelets making one element.  In a NONBONDED_FIELD or
                                   ///<   NONBONDED_ATOMIC mesh, the coefficients are tricubic
                                   ///<   splines for each element.
  const double coeff_scale;        ///< Scaling factor to take mesh coefficients from internal
                                   ///<   units (kcal/mol, kcal/mol-A, kcal/mol-A2, and
                                   ///<   kcal/mol-A3, all based on the same factor) into a
                                   ///<   fixed-precision representation.  This conversion is only
                                   ///<   performed when accumulating meshes.  Meshes used in
                                   ///<   actual energy calculations will be assumed to be of real
                                   ///<   scalar types (float, double).
  const float coeff_scale_f;       ///< Single-precision variant of coeff_scale
  
  // The following data can be necessary for interpreting potentials read from the mesh.
  const double probe_radius;    ///< The probe radius to use when mapping a Lennard-Jones mesh with
                                ///<   geometric combining rules, or when mapping a clash potential
                                ///<   to the mesh
  const double well_depth;      ///< The Lennard-Jones well depth for a test particle used to map a
                                ///<   Lennard-Jones potential of any mixing rules to the mesh.
  const double occ_cost;        ///< Energetic penalty of an occlusion interaction
  const MeshBasicsKit mbss;     ///< Collection of pointers for essential elements about the
                                ///<   molecular system underlying the mesh, including a neighbor
                                ///<   list for each cell (if relevant) and a mask of frozen atoms
};
  
/// \brief A workspace for constructing a pure potential mesh based on the frozen atoms of a
///        large molecule.  If the large molecule has nonrigid components, they must be excluded
///        from contributing to the grid.  In addition, any atoms up to 1:4 (connected by three
///        bonds or less) must also be excluded from the grid-based potential.  Computations on
///        these atoms will not be accurate off the grid, but since they are frozen the
///        consequences are mitigated.
template <typename T> class BackgroundMesh {
public:

  /// \brief The constructor takes all dimension parameters plus an indication of what type of
  ///        potential, the molecular system, and what mask of atoms is to be mapped.  Variants
  ///        include different ways to define the limits of the mesh.  If a GPU is available, it
  ///        will be used to compute the mesh.
  ///
  /// Overloaded:
  ///   - Various combinations of parameters specific to relevant mesh kinds
  ///   - Templated constructors for including a coordinate series rather than a single structure
  ///     in occlusion field and specific nonbonded field meshes
  ///
  /// \param ag                 System topology (markings in its mobile_atoms array will be used
  ///                           to determine which atoms to map)
  /// \param cf                 Cartesian coordinates of all particles
  /// \param buffer             Breadth around the molecule drawing an orthhombic region in which
  ///                           to map the mesh
  /// \param mesh_bounds        Boundaries of the mesh, a six-element vector describing the lower
  ///                           and upper Cartesian X, Y, and Z limits of an orthorhombic region
  ///                           in which to define a rectilinear mesh.
  /// \param spacing            Grid spacings for the mesh cells.  Provide either a single number
  ///                           or three dimensions for the Cartesian X, Y, and Z widths of
  ///                           rectilinear elements.
  /// \param measurements_in    A full description of the mesh parameters.  This provides a means
  ///                           to define non-orthorhombic meshes.
  /// \param clash_distance_in  The absolute distance at which electrostatic interactions will be
  ///                           deemed in conflict and a softcore potential will take over
  /// \param clash_ratio_in     The ratio of the pairwise van-der Waals sigma parameter at which a
  ///                           van-der Waals interaction will be declared in conflict and the
  ///                           softcore potential form will take over
  /// \param scale_bits         Number of bits after the decimal when taking coordinates in units
  ///                           of Angstroms into the fixed-precision format of the mesh framework
  /// \param averaging_order    The width of the discrete point patch from an occlusion mesh
  ///                           measured for a single frame that will be averaged to compute a
  ///                           contribution to a continuous field density
  ///                           Angstroms into the fixed-precision format of the mesh framework
  /// \param launcher           Manager for launching mesh-based kernels (dispenses launch grid
  ///                           parameters)
  /// \param prec               Precision model to use in GPU-based mesh construction computations
  /// \param availability       Indicate whether the particle coordinate data on the CPU or GPU
  ///                           should be trusted
  /// \{
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in,
                 const MeshParameters &measurements_in = MeshParameters(),
                 const PrecisionModel build_precision_in = PrecisionModel::SINGLE);

  // General meshes with all available details  
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in,
                 const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                 double probe_radius_in, double well_depth_in, VdwCombiningRule mixing_protocol_in,
                 const MeshParameters &measurements_in = MeshParameters(),
                 const std::vector<double> &probe_sigma = {},
                 const std::vector<double> &probe_epsilon = {},
                 double clash_distance_in = default_mesh_elec_damping_range,
                 double clash_ratio_in = default_mesh_vdw_damping_ratio,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in,
                 const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                 double probe_radius_in, double well_depth_in, VdwCombiningRule mixing_protocol_in,
                 const MeshParameters &measurements_in = MeshParameters(),
                 const std::vector<double> &probe_sigma = {},
                 const std::vector<double> &probe_epsilon = {},
                 double clash_distance_in = default_mesh_elec_damping_range,
                 double clash_ratio_in = default_mesh_vdw_damping_ratio,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);
  
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, double buffer, double spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 const std::vector<double> &probe_epsilon = {},
                 double clash_distance_in = default_mesh_elec_damping_range,
                 double clash_ratio_in = default_mesh_vdw_damping_ratio,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph &ag_in,
                 const CoordinateFrame &cf_in, double buffer, double spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 const std::vector<double> &probe_epsilon = {},
                 double clash_distance_in = default_mesh_elec_damping_range,
                 double clash_ratio_in = default_mesh_vdw_damping_ratio,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);
  
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, double buffer, const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 const std::vector<double> &probe_epsilon = {},
                 double clash_distance_in = default_mesh_elec_damping_range,
                 double clash_ratio_in = default_mesh_vdw_damping_ratio,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph &ag_in,
                 const CoordinateFrame &cf_in, double buffer, const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 const std::vector<double> &probe_epsilon = {},
                 double clash_distance_in = default_mesh_elec_damping_range,
                 double clash_ratio_in = default_mesh_vdw_damping_ratio,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);
  
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, const std::vector<double> &mesh_bounds,
                 double spacing, int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 const std::vector<double> &probe_epsilon = {},
                 double clash_distance_in = default_mesh_elec_damping_range,
                 double clash_ratio_in = default_mesh_vdw_damping_ratio,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph &ag_in,
                 const CoordinateFrame &cf_in, const std::vector<double> &mesh_bounds,
                 double spacing, int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 const std::vector<double> &probe_epsilon = {},
                 double clash_distance_in = default_mesh_elec_damping_range,
                 double clash_ratio_in = default_mesh_vdw_damping_ratio,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, const std::vector<double> &mesh_bounds,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 const std::vector<double> &probe_epsilon = {},
                 double clash_distance_in = default_mesh_elec_damping_range,
                 double clash_ratio_in = default_mesh_vdw_damping_ratio,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph &ag_in,
                 const CoordinateFrame &cf_in, const std::vector<double> &mesh_bounds,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 const std::vector<double> &probe_epsilon = {},
                 double clash_distance_in = default_mesh_elec_damping_range,
                 double clash_ratio_in = default_mesh_vdw_damping_ratio,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  // Electrostatics-based meshes
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, double buffer, double spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 double clash_distance_in = default_mesh_elec_damping_range,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph &ag_in,
                 const CoordinateFrame &cf_in, double buffer, double spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 double clash_distance_in = default_mesh_elec_damping_range,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, double buffer,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 double clash_distance_in = default_mesh_elec_damping_range,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph &ag_in,
                 const CoordinateFrame &cf_in, double buffer,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 double clash_distance_in = default_mesh_elec_damping_range,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, const std::vector<double> &mesh_bounds,
                 double spacing, int scale_bits_in = default_globalpos_scale_bits,
                 double clash_distance_in = default_mesh_elec_damping_range,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph &ag_in,
                 const CoordinateFrame &cf_in, const std::vector<double> &mesh_bounds,
                 double spacing, int scale_bits_in = default_globalpos_scale_bits,
                 double clash_distance_in = default_mesh_elec_damping_range,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph *ag_in,
                 const CoordinateFrame *cf_in, const std::vector<double> &mesh_bounds,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 double clash_distance_in = default_mesh_elec_damping_range,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph &ag_in,
                 const CoordinateFrame &cf_in, const std::vector<double> &mesh_bounds,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 double clash_distance_in = default_mesh_elec_damping_range,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  // Clash-based meshes
  BackgroundMesh(GridDetail kind_in, double probe_radius_in, VdwCombiningRule mixing_protocol_in,
                 const AtomGraph *ag_in, const CoordinateFrame *cf_in, double buffer,
                 double spacing, int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, double probe_radius_in, VdwCombiningRule mixing_protocol_in,
                 const AtomGraph &ag_in, const CoordinateFrame &cf_in, double buffer,
                 double spacing, int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, double probe_radius_in, VdwCombiningRule mixing_protocol_in,
                 const AtomGraph *ag_in, const CoordinateFrame *cf_in, double buffer,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, double probe_radius_in, VdwCombiningRule mixing_protocol_in,
                 const AtomGraph &ag_in, const CoordinateFrame &cf_in, double buffer,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, double probe_radius_in, VdwCombiningRule mixing_protocol_in,
                 const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                 const std::vector<double> &mesh_bounds, double spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, double probe_radius_in, VdwCombiningRule mixing_protocol_in,
                 const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                 const std::vector<double> &mesh_bounds, double spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, double probe_radius_in, VdwCombiningRule mixing_protocol_in,
                 const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                 const std::vector<double> &mesh_bounds, const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  BackgroundMesh(GridDetail kind_in, double probe_radius_in, VdwCombiningRule mixing_protocol_in,
                 const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                 const std::vector<double> &mesh_bounds, const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 const std::vector<double> &probe_sigma = {},
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  // General meshes with all available details
  template <typename Tcoord>
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in,
                 const AtomGraph &ag_in, const CoordinateSeries<Tcoord> &cs_in,
                 double probe_radius_in, double well_depth_in, VdwCombiningRule mixing_protocol_in,
                 const MeshParameters &measurements_in = MeshParameters(),
                 int averaging_order = default_mesh_density_averaging_order,
                 const std::vector<double> &probe_sigma = {},
                 const std::vector<double> &probe_epsilon = {},
                 double clash_distance_in = default_mesh_elec_damping_range,
                 double clash_ratio_in = default_mesh_vdw_damping_ratio,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  template <typename Tcoord>
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph &ag_in,
                 const CoordinateSeries<Tcoord> &cs_in, double buffer, double spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 int averaging_order = default_mesh_density_averaging_order,
                 const std::vector<double> &probe_sigma = {},
                 const std::vector<double> &probe_epsilon = {},
                 double clash_distance_in = default_mesh_elec_damping_range,
                 double clash_ratio_in = default_mesh_vdw_damping_ratio,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  template <typename Tcoord>
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph &ag_in,
                 const CoordinateSeries<Tcoord> &cs_in, double buffer,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 int averaging_order = default_mesh_density_averaging_order,
                 const std::vector<double> &probe_sigma = {},
                 const std::vector<double> &probe_epsilon = {},
                 double clash_distance_in = default_mesh_elec_damping_range,
                 double clash_ratio_in = default_mesh_vdw_damping_ratio,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  template <typename Tcoord>
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph &ag_in,
                 const CoordinateSeries<Tcoord> &cs_in, const std::vector<double> &mesh_bounds,
                 double spacing, int scale_bits_in = default_globalpos_scale_bits,
                 int averaging_order = default_mesh_density_averaging_order,
                 const std::vector<double> &probe_sigma = {},
                 const std::vector<double> &probe_epsilon = {},
                 double clash_distance_in = default_mesh_elec_damping_range,
                 double clash_ratio_in = default_mesh_vdw_damping_ratio,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  template <typename Tcoord>
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, double probe_radius_in,
                 double well_depth_in, VdwCombiningRule mixing_protocol_in, const AtomGraph &ag_in,
                 const CoordinateSeries<Tcoord> &cs_in, const std::vector<double> &mesh_bounds,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 int averaging_order = default_mesh_density_averaging_order,
                 const std::vector<double> &probe_sigma = {},
                 const std::vector<double> &probe_epsilon = {},
                 double clash_distance_in = default_mesh_elec_damping_range,
                 double clash_ratio_in = default_mesh_vdw_damping_ratio,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  // Electrostatics-based meshes
  template <typename Tcoord>
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph &ag_in,
                 const CoordinateSeries<Tcoord> &cs_in, double buffer, double spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 double clash_distance_in = default_mesh_elec_damping_range,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  template <typename Tcoord>
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph &ag_in,
                 const CoordinateSeries<Tcoord> &cs_in, double buffer,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 double clash_distance_in = default_mesh_elec_damping_range,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  template <typename Tcoord>
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph &ag_in,
                 const CoordinateSeries<Tcoord> &cs_in, const std::vector<double> &mesh_bounds,
                 double spacing, int scale_bits_in = default_globalpos_scale_bits,
                 double clash_distance_in = default_mesh_elec_damping_range,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  template <typename Tcoord>
  BackgroundMesh(GridDetail kind_in, NonbondedPotential field_in, const AtomGraph &ag_in,
                 const CoordinateSeries<Tcoord> &cs_in, const std::vector<double> &mesh_bounds,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 double clash_distance_in = default_mesh_elec_damping_range,
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  // Clash-based meshes
  template <typename Tcoord>
  BackgroundMesh(GridDetail kind_in, double probe_radius_in, const AtomGraph &ag_in,
                 const CoordinateSeries<Tcoord> &cs_in, double buffer, double spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 int averaging_order = default_mesh_density_averaging_order,
                 const std::vector<double> &probe_sigma = {},
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  template <typename Tcoord>
  BackgroundMesh(GridDetail kind_in, double probe_radius_in, const AtomGraph &ag_in,
                 const CoordinateSeries<Tcoord> &cs_in, double buffer,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 int averaging_order = default_mesh_density_averaging_order,
                 const std::vector<double> &probe_sigma = {},
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  template <typename Tcoord>
  BackgroundMesh(GridDetail kind_in, double probe_radius_in, const AtomGraph &ag_in,
                 const CoordinateSeries<Tcoord> &cs_in, const std::vector<double> &mesh_bounds,
                 double spacing, int scale_bits_in = default_globalpos_scale_bits,
                 int averaging_order = default_mesh_density_averaging_order,
                 const std::vector<double> &probe_sigma = {},
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);

  template <typename Tcoord>
  BackgroundMesh(GridDetail kind_in, double probe_radius_in, const AtomGraph &ag_in,
                 const CoordinateSeries<Tcoord> &cs_in, const std::vector<double> &mesh_bounds,
                 const std::vector<double> &spacing,
                 int scale_bits_in = default_globalpos_scale_bits,
                 int averaging_order = default_mesh_density_averaging_order,
                 const std::vector<double> &probe_sigma = {},
                 PrecisionModel prec = PrecisionModel::SINGLE,
                 const MeshKlManager &launcher = MeshKlManager(),
                 HybridTargetLevel availability = HybridTargetLevel::HOST);
  /// \}

  /// \brief Copy and move constructors as well as assignment operators can be set to their
  ///        defaults because the object is composed of others, which may contain POINTER-kind
  ///        Hybrid objects but have their own copy and move assignment operators.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Assignments are made based on this pre-existing object.
  /// \{
  BackgroundMesh(const BackgroundMesh<T> &original) = default;
  BackgroundMesh(BackgroundMesh<T> &&original) = default;
  BackgroundMesh& operator=(const BackgroundMesh<T> &original) = default;
  BackgroundMesh& operator=(BackgroundMesh<T> &&original) = default;
  /// \}

  /// \brief Get an object describing the mesh dimensions
  const MeshParameters& getDimensions() const;
  
  /// \brief Get the collection of molecular details underlying the mesh, including the topology
  ///        and associated coordinates, plus any neighbor list.
  ///
  /// Overloaded:
  ///   - Return a const reference for a const object
  ///   - Return a mutable pointer for a non-const object
  /// \{
  const MeshFoundation& getMolecularBasis() const;
  MeshFoundation* getMolecularBasis();
  /// \}

  /// \brief Get a const pointer to the topology responsible for creating this mesh.
  const AtomGraph* getTopologyPointer() const;

  /// \brief Get a const pointer to the coordinates responsible for creating this mesh.
  const CoordinateFrame* getCoordinatePointer() const;

  /// \brief Get a const reference to the array of structures (each stored as a CoordinateFrame)
  ///        reponsible for creating this mesh.
  template <typename Tcoord>
  const CoordinateSeries<Tcoord>* getEnsemblePointer() const;

  /// \brief Get the codified data type of the ensemble.
  size_t getEnsembleTypeCode() const;
  
  /// \brief Get the type of mesh.
  GridDetail getMeshKind() const;

  /// \brief Get the non-bonded potential expressed on the mesh.  This will return an error if the
  ///        mesh is not associated with a non-bonded potential.
  NonbondedPotential getNonbondedPotential() const;

  /// \brief Get the probe radius.  Valid for all types of meshes.
  double getProbeRadius() const;

  /// \brief Get the non-bonded probe well depth.  This is valid only for non-bonded fields, in
  ///        particular those associated with a van-der Waals potential.
  double getWellDepth() const;

  /// \brief Get the absolute distance at which softcore electrostatic interactions take over.
  ///        This is provided for utility and accesses the nested nonbonded_model MeshForceField
  ///        class member variable.
  double getClashDistance() const;

  /// \brief Get the ratio of the van-der Waals (Lennard-Jones) sigma parameter for interacting
  ///        pairs of particles at which softcore van-der Waals interactions take over.  This is
  ///        provided for utility and accesses the nested nonbonded_model MeshForceField
  ///        class member variable.
  double getClashRatio() const;

  /// \brief Get the penalty associated with a collision on an occlusion mesh.
  double getOcclusionPenalty() const;
  
  /// \brief Get the mixing protocol used for van-der Waals (as expressed by the Lennard-Jones
  ///        model) interactions.
  VdwCombiningRule getCombiningRule() const;

  /// \brief Get the precisionmodel under which the mesh was calculated.
  PrecisionModel getBuildPrecision() const;
  
  /// \brief Get the scaling factor for mesh coefficients in fixed-precision representations.
  double getCoefficientScalingFactor() const;

  /// \brief Get the number of bits after the decimal in fixed-precision mesh coefficient
  ///        representations.
  int getCoefficientScalingBits() const;

  /// \brief Get an abstract of the mesh.
  ///
  /// Overloaded:
  ///   - Get the reader for a const object
  ///   - Get the writer for a mutable object
  ///
  /// \param tier  Obtain the data at the level of the CPU host or GPU device
  /// \{
  const BackgroundMeshReader<T> data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  BackgroundMeshWriter<T> data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Get an abstract of the mesh with any templating removed.
  ///
  /// Overloaded:
  ///   - Get the reader for a const object
  ///   - Get the writer for a mutable object
  ///
  /// \param tier  Obtain the data at the level of the CPU host or GPU device
  /// \{
  const BackgroundMeshReader<void>
  templateFreeData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  BackgroundMeshWriter<void>
  templateFreeData(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  /// \brief Get an abstract of the mesh with any templating removed.
  ///
  /// Overloaded:
  ///   - Get the reader for a const object
  ///   - Get the writer for a mutable object
  ///
  /// \{
  BackgroundMeshWriter<void> deviceViewToTemplateFreeHostData();
  const BackgroundMeshReader<void> deviceViewToTemplateFreeHostData() const;
  /// \}
#  endif
#endif

  /// \brief Get the abstract of the mesh in the precision of the mesh's coefficient data.
  ///
  /// \param tier  Obtain the data at the level of the CPU host or GPU device
  const MeshFFKit<T>
  getNonbondedKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the abstract of the mesh in double precision.
  ///
  /// \param tier  Obtain the data at the level of the CPU host or GPU device
  const MeshFFKit<double>
  getReferenceNonbondedKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a template-free form of the of the mesh nonbonded abstract in the mesh coefficient
  ///        data type, useful for passing between the C++- and HPC-compiled code objects.
  ///
  /// \param tier  Obtain the data at the level of the CPU host or GPU device
  const MeshFFKit<void>
  templateFreeNonbondedKit(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  /// \brief Get the template-free form of the critical non-bonded parameters in their "reference"
  ///        double-precision form, visible on the device but mounted on the host.
  const MeshFFKit<double> deviceViewToReferenceNonbondedKit() const;

  /// \brief Get the template-free form of the critical non-bonded parameters, visible on the
  ///        device but mounted on the host.
  const MeshFFKit<void> deviceViewToTemplateFreeNonbondedKit() const;
#  endif
  /// \brief Upload all data to the device layer.
  void upload();

  /// \brief Download all data from the device layer.
  void download();

  /// \brief Upload the "frame" data--including rulers for the a, b, and c mesh axes as well as
  ///        the bitmask of relevant atoms from the referenced topology, to the mesh.
  void uploadFrame();

  /// \brief Download the "frame" data--including rulers for the a, b, and c mesh axes as well as
  ///        the bitmask of relevant atoms from the referenced topology, to the mesh.
  void downloadFrame();

  /// \brief Upload the non-bonded constants determined for softcore potentials on this mesh.
  void uploadForceField();

  /// \brief Download the non-bonded constants determined for softcore potentials on this mesh.
  void downloadForceField();
#endif
  
  /// \brief Set the dimensions of the mesh.  Variants of this member function that do not include
  ///        a coordinate pointer or topology may assume that one has already been set.  Variants
  ///        that do accept a coordinate pointer will set that as the mesh's coordinate set.
  ///
  /// Overloaded:
  ///   - Take coordinate and topology pointers with a single parameter to define the mesh limits.
  ///   - Take a single parameter to define the mesh limits around a set of coordinates defined by
  ///     a particular topology, and assume that the mesh already has these things.
  ///   - Take the minimum and maximum Cartesian X, Y, and Z limits of a rectilinear mesh (this
  ///     will not require a system to implement)
  ///   - Take a single parameter to define isotropic mesh elements
  ///   - Take three or nine parameters to define anisotropic rectilinear or triclinic mesh
  ///     elements
  ///   - Indicate the scaling bits to be used in fixed-precision arithmetic when determining
  ///     positions on the mesh, or leave that parameter unchanged
  ///
  /// \param ag_in          New topology to build the mesh around (this will be incorporated into
  ///                       the object)
  /// \param cf_in          New coordinate set to build the mesh around (this will be incorporated
  ///                       into the object)
  /// \param padding        Region around the molecule of interest to spread the mesh
  /// \param mesh_bounds    Minimum and maximum Cartesian limits of the mesh
  /// \param spacing        Width (and length, and height, or bounding vectors) of mesh elements
  /// \param scale_bits_in  The number of bits after the decimal to use in fixed-precision
  ///                       computations of particle positions on the grid.  The default value of
  ///                       -100 is nonsensical and will be interpreted as "leave the current
  ///                       setting as it is." Limiting the number of function overloads in this
  ///                       way is a step towards controlling code bloat due to this large,
  ///                       templated class.
  /// \{
  void setMeshParameters(const AtomGraph *ag_in, const CoordinateFrame *cf_in, double padding,
                         double spacing, int scale_bits_in = -100);

  void setMeshParameters(const AtomGraph *ag_in, const CoordinateFrame *cf_in, double padding,
                         const std::vector<double> &spacing, int scale_bits_in = -100);

  void setMeshParameters(double padding, double spacing, int scale_bits_in = -100);

  void setMeshParameters(double padding, const std::vector<double> &spacing,
                         int scale_bits_in = -100);

  void setMeshParameters(const std::vector<double> &mesh_bounds, double spacing,
                         int scale_bits_in = -100);

  void setMeshParameters(const std::vector<double> &mesh_bounds,
                         const std::vector<double> &spacing, int scale_bits_in = -100);
  /// \}

  /// \brief Set the probe radius, meaning either the Lennard-Jones potential sigma radius or
  ///        the clash probe radius, depending on the mesh type.  This includes a validity check.
  ///
  /// \param probe_radius_in  The probe radius to set
  void setProbeRadius(double probe_radius_in);

  /// \brief Set the Lennard-Jones well depth for the probe that will generate the potential.  This
  ///        includes a validity check.
  ///
  /// \param well_depth_in  The well depth to set
  void setWellDepth(double well_depth_in);

  /// \brief Set the penalty associated with striking occupied volume on an occlusion mesh.  This
  ///        is intended to be modified afeter making the mesh, if at all, as the actual value has
  ///        no bearing on the way the bitmask is constructed.
  ///
  /// \param occlusion_penalty_in  The occlusion penalty to set
  void setOcclusionPenalty(double occlusion_penalty_in);
  
  /// \brief Check the combining rule that will be used to make the probe interact with the
  ///        receptor on any mesh (with the exception of an electrostatic field).  With geometric
  ///        combining rules, it is possible to tailor a single mesh for all particles that might
  ///        interact with the receptor.  However, with Lorentz-Berthelot rules or any case of
  ///        non-conformant pair rules, new grids are required for a rigorous description of each
  ///        particle type that might interact with the mesh potential.
  ///
  /// \param mixing_protocol_in  The method to use
  /// \param probe_sigma         Probe pairwise sigma parameters with each Lennard-Jones type in
  ///                            the topology (the topology is assumed to be set by the time this
  ///                            function is called).  For NBFIX combining rules, this array must
  ///                            have entries for each Lennard-Jones type.  For other combining
  ///                            rules, this array should have no contents.
  /// \param probe_epsilon       Probe pairwise epsilon parameters (see probe_sigma for
  ///                            stipulations about this array)
  void validateCombiningRule(VdwCombiningRule mixing_protocol_in,
                             const std::vector<double> &probe_sigma,
                             const std::vector<double> &probe_epsilon);

  /// \brief Set the number of bits after the decimal to be used in fixed-precision representations
  ///        of coordinates.
  ///
  /// \param scaling_bits_in  The desired bit count
  void setCoefficientScalingBits(int scaling_bits_in);

  /// \brief Compute the appropriate field for the mesh.  This is called automatically by the
  ///        constructor if enough information is provided.
  ///
  /// \param launcher      Manager for mesh kernels, dispenses launch parameters
  /// \param prec          Precision model to use in GPU-based mesh construction computations
  /// \param availability  Level at which the data is expected to be available
  void computeField(const MeshKlManager &launcher = MeshKlManager(),
                    PrecisionModel prec = PrecisionModel::SINGLE,
                    HybridTargetLevel availability = HybridTargetLevel::HOST,
                    const std::vector<double> &probe_sigma = {},
                    const std::vector<double> &probe_epsilon = {});

private:

  /// Measurements for the mesh, including various transformation matrices and fixed-precision
  /// representations.
  MeshParameters measurements;

  /// The detail level of the mesh--this carries implications for the allowed data types.  An
  /// "OCCLUSION" mesh must be of type ullint (64-bit unsigned integer, to encode 4 x 4 x 4
  /// cubelets of bitwise information for "occluded" or "non-occluded." Such meshes can be very
  /// dense, as every mesh element constitutes 4 x 4 x 4 cubelets, or a total of 512 bytes of
  /// information, as much as 20k mesh points on each side.  In contrast, a "NONBONDED_FIELD" type
  /// grid must have data type float or double to encode the appropriate non-bonded potential or
  /// influence field (each Lennard-Jones type and Generalized Born parameter combination will get
  /// its own field).  A "NONBONDED_ATOMIC" type mesh stores field data in floating-point scalar
  /// data as well as local neighbor lists in an additional array of integer data, with a bounds
  /// array.
  GridDetail kind;

  /// The type of non-bonded potential to map to the mesh.  This is irrelevant for "OCCLUSION"-kind
  /// meshes (the effective potential is a clash based on the Lennard-Jones sigma radii and a probe
  /// width).
  NonbondedPotential field;

  /// The probe radius to use in computing clashes for an OCCLUSION-type mesh.
  double probe_radius;

  /// Lennard-Jones well depth of the probe particle's self interaction.  If the mesh is of type
  /// NONBONDED_FIELD or NONBONDED_ATOMIC, the probe_radius will be taken to mean the particle's
  /// self-interacting sigma parameter, and this is its epsilon.  For electrostatic fields, neither
  /// this number nor probe_radius has any relevance.
  double well_depth;

  /// The penalty associated with placing an atom inside an occupied volume element when sampling
  /// an occlusion mesh.
  double occlusion_penalty;
  
  /// The precision used to build the mesh.  Meshes can have float-type coefficient data after
  /// being constructed in double-precision math, and vice-versa (although the latter case would
  /// be most useful as an experiment to capture the effects of low-precision calculations).
  PrecisionModel build_precision;

  /// Cartesian coordinates of the three mesh axes, containing tick marks for each element.
  MeshRulers tick_marks;

  /// The non-bonded model components laying out how the mesh's underlying topology and the
  /// probe particle will interact.
  MeshForceField<T> nonbonded_model;
  
  /// Coefficients for tricubic spline functions spanning each grid element.  Coefficients for
  /// element {Ea, Eb, Ec} traversing the "a", "b", and "c" axes, respectively, are held
  /// consecutively in a 64-element series with order {i, j, k} = (0, 0, 0), (0, 0, 1), ...,
  /// (0, 0, 3), (1, 0, 0), (1, 0, 1), ..., (1, 0, 3), ..., (3, 3, 3).  The 64-element series for
  /// element {Ea, Eb, Ec} will be found at 64 * ((((Ec * nb) + Eb) * na) + Ea), where {na, nb, nc}
  /// are the mesh dimensions along each axis.  This matches the order of elements found in the
  /// output of the TricubicCell() object (see Math/tricubic_cell.h).  In an occlusion mesh, the
  /// 64 numbers are 64-bit unsigned integers, each a bitstring for a 4 x 4 x 4 subdivision of the
  /// 4 x 4 x 4 divisions of the space within the element.  Division {Di, Dj, Dk} is found in array
  /// index ((((Dk * 4) + Dj) * 4) + Di), and the subdivisions of the 64 bits follow a similar
  /// pattern for optimal localization, access, and coloring.
  Hybrid<T> coefficients;

  // Mesh coefficients may come in fixed-precision format.  These coefficients mirror contents of
  // other templated classes to track the conversion.
  int coefficient_scale_bits;        ///< Number of bits after the decimal in mesh coefficients,
                                     ///<   the base-2 logarithm of coefficient_scale 
  double coefficient_scale;          ///< Scaling factor for fixed-precision representations of
                                     ///<   mesh coefficients.
  double inverse_coefficient_scale;  ///< Scaling factor to take fixed-precision coefficient
                                     ///<   representations back into internal units of kcal/mol,
                                     ///<   kcal/mol-A, etc.  It should be noted that the one
                                     ///<   scaling factor applies to potentials, derivatives,
                                     ///<   and all mixed partial derivatives, despite these
                                     ///<   quantities all having different units.

  /// The molecular basis of the mesh
  MeshFoundation basis;
  
  /// \brief Allocate memory for the mesh coefficients, axis coordinates, and frozen atoms mask.
  ///        This does not allocate space for the concatenated mesh element neighbor lists or their
  ///        bounds array does it allocate space for the non-bonded parameter arrays.
  void allocate();

  /// \brief Given the underlying topology and critical information about the mesh probe (sometimes
  ///        referred to as the "solvent probe"), compute the necessary pairwise parameters for
  ///        interactions with each atom type.  Given the switching distances, compute softcore
  ///        polynomial coefficients.  Given the stencil for interpolation (for which this routine
  ///        must be called after the measurements member variable is set, as is done in the
  ///        constructor), compute the appropriate softcore polynomials.
  ///
  /// \param mixing_protocol_in  The method to use in combining Lennard-Jones parameters (this will
  ///                            be validated based on the field type to be projected onto the
  ///                            mesh).  This is passed on to the nonbonded_model member variable.
  /// \param clash_ratio_in      Ratio of the Lennard-Jones sigma parameter at which the softcore
  ///                            handoff occurs (passed on to the nonbonded_model member variable)
  /// \param clash_distance_in   Absolute distance at which electrostatic interactions transition
  ///                            to the softcore potential
  /// \param probe_sigma         Probe pairwise sigma parameters with each Lennard-Jones type in
  ///                            the topology (the topology is assumed to be set by the time this
  ///                            function is called).  For NBFIX combining rules, this array must
  ///                            have entries for each Lennard-Jones type.  For other combining
  ///                            rules, this array should have no contents.
  /// \param probe_epsilon       Probe pairwise epsilon parameters (see probe_sigma for
  ///                            stipulations about this array)
  void setNonbondedModel(VdwCombiningRule mixing_protocol_in, double clash_ratio_in,
                         double clash_distance_in, const std::vector<double> &probe_sigma,
                         const std::vector<double> &probe_epsilon);

  /// \brief Certain types of meshes are restricted to certain data types.  This function will
  ///        ensure the correct relationships.
  void validateMeshKind() const;

  /// \brief Validate the fixed-precision scaling based on the mesh's data type.  Single-precision
  ///        nonbonded field meshes and clash maps must not overflow 64-bit positional
  ///        representations.
  void validateScalingBits() const;

  /// \brief Color the occlusion mesh based on a particular set of coordinates.
  ///
  /// \param launcher      Manager for launching mesh-based kernels (dispenses launch grid
  ///                      parameters)
  /// \param prec          Precision model to use in GPU-based mesh construction computations
  /// \param availability  Level at which the data is expected to be available
  /// \param probe_sigma   A list of pairwise sigma parameters for the solvent probe interacting
  ///                      with each Lennard-Jones atom type in the underlying topology.  This must
  ///                      be provided if and only if the van-der Waals mixing rule is set to
  ///                      NBFIX, and in such a case it must address every Lennard-Jones atom type.
  void colorOcclusionMesh(const MeshKlManager &launcher,
                          PrecisionModel prec = PrecisionModel::SINGLE,
                          HybridTargetLevel availability = HybridTargetLevel::HOST,
                          const std::vector<double> &probe_sigma = {});

  /// \brief Fill in an occlusion field based on a one-to-one transformation of the occupancy seen
  ///        for a series of structures.  The structures must be aligned, by means outside of the
  ///        mesh object itself, prior to using them to assemble the mesh.  Descriptions of input
  ///        parameters for this function follow from colorOcclusionMesh() above.
  void mapOccupancyField(const MeshKlManager &launcher,
                         PrecisionModel prec = PrecisionModel::SINGLE,
                         HybridTargetLevel availability = HybridTargetLevel::HOST,
                         const std::vector<double> &probe_sigma = {});
  
  /// \brief Common function for various NONBONDED_FIELD meshes making use of a series of grids
  ///        for the function values plus partial and mixed partial derivatives.  This works for
  ///        CPU-based code.  For an M x N x P mesh, annd grids measure (M + 1) x (N + 1) x (P + 1)
  ///        points.  The grids are populated with Cartesian X, Y, and Z derivatives and will feed
  ///        into functions for computing spline coefficients for any shape of mesh.
  ///
  /// Overloaded:
  ///   - Provide only the grids of partial derivatives needed for an orthorhombic mesh
  ///   - Provide the grids of partial drivatives needed for a general, triclinic mesh
  ///
  /// \param u_grid       Grid of function values
  /// \param dudx_grid    Grid of Cartesian X partial derivatives
  /// \param dudy_grid    Grid of Cartesian Y partial derivatives
  /// \param dudz_grid    Grid of Cartesian Z partial derivatives
  /// \param dudxx_grid   Grid of Cartesian X second derivatives, for non-orthorhombic meshes
  /// \param dudxy_grid   Grid of Cartesian X and Y double derivatives
  /// \param dudxz_grid   Grid of Cartesian X and Z double derivatives
  /// \param dudyy_grid   Grid of Cartesian Y second derivatives, for non-orthorhombic meshes
  /// \param dudyz_grid   Grid of Cartesian Y and Z double derivatives
  /// \param dudxxx_grid  Grid of X triple derivatives, for non-orthorhombic meshes
  /// \param dudxxy_grid  Grid of X, X, Y partial derivatives, for non-orthorhombic meshes
  /// \param dudxxy_grid  Grid of X, X, Z partial derivatives, for non-orthorhombic meshes
  /// \param dudxxy_grid  Grid of X, Y, Y partial derivatives, for non-orthorhombic meshes
  /// \param dudxyz_grid  Grid of X, Y, Z triple partial derivatives
  /// \{
  void computeCoefficients(const std::vector<double> &u_grid,
                           const std::vector<double> &dudx_grid,
                           const std::vector<double> &dudy_grid,
                           const std::vector<double> &dudz_grid,
                           const std::vector<double> &dudxx_grid,
                           const std::vector<double> &dudxy_grid,
                           const std::vector<double> &dudxz_grid,
                           const std::vector<double> &dudyy_grid,
                           const std::vector<double> &dudyz_grid,
                           const std::vector<double> &dudxxx_grid,
                           const std::vector<double> &dudxxy_grid,
                           const std::vector<double> &dudxxz_grid,
                           const std::vector<double> &dudxyy_grid,
                           const std::vector<double> &dudxyz_grid,
                           const TricubicStencil &tc_weights);

  void computeCoefficients(const std::vector<double> &u_grid,
                           const std::vector<double> &dudx_grid,
                           const std::vector<double> &dudy_grid,
                           const std::vector<double> &dudz_grid,
                           const std::vector<double> &dudxy_grid,
                           const std::vector<double> &dudxz_grid,
                           const std::vector<double> &dudyz_grid,
                           const std::vector<double> &dudxyz_grid,
                           const TricubicStencil &tc_weights);
  /// \}

  /// \brief The CPU-based algorithm for coloring the non-bonded field mesh.  Additional templating
  ///        is provided to enable the use of single- and double-precision numerics.
  ///
  /// \param nbk           Non-bonded parameters from the underlying topology
  /// \param sigma_table   A table of sigma radii by which the probe interacts with the atoms of
  ///                      the mesh-mapped receptor
  /// \param eps_table     A table of well depth values by which the probe interacts with the
  ///                      atoms of the mesh-mapped receptor
  template <typename Tcalc, typename Tcalc2, typename Tcalc3, typename Tcalc4>
  void colorNonbondedField(const NonbondedKit<Tcalc> &nbk, const TricubicStencil &tc_weights,
                           const std::vector<double> &sigma_table = {},
                           const std::vector<double> &eps_table = {});
  
  /// \brief Map the non-bonded potential due to the receptor's rigid atoms, without any neighbor
  ///        lists (which implies element-specific mesh exclusions).  Descriptions of input
  ///        arguments follow from colorNonbondedField() above, in addition to:
  ///
  /// \param launcher      Manager for launching mesh-based kernels (dispenses launch grid
  ///                      parameters)
  /// \param prec          Precision model to use in GPU-based mesh construction computations
  /// \param availability  Level at which the data is expected to be available
  void mapPureNonbondedPotential(const MeshKlManager &launcher,
                                 PrecisionModel prec = PrecisionModel::SINGLE,
                                 HybridTargetLevel availability = HybridTargetLevel::HOST,
                                 const std::vector<double> &sigma_table = {},
                                 const std::vector<double> &eps_table = {});
};

/// \brief Restore the true data type of a void-casted BackgroundMesh abstract.  This will be
///        compiled by both C++ and HPC compilers that build executable objects dependent on the
///        BackgroundMesh class, such that executable objects built by either compiler will have
///        the means to remove and restore templated character to abstracts of BackgroundMesh.
///
/// Overloaded:
///   - Produce a read-only templated abstract based on a read-only void-casted abstract
///   - Produce a writeable, templated abstract based on a writeable void-casted abstract
///   - Overloads in other libraries produce outputs to match their own input types
///
/// \param rasa  Abstract of the original mesh object, free of templated pointers (inspired by the
///              phrase "tabula rasa", a blank slate)
/// \{
template <typename T> BackgroundMeshWriter<T> restoreType(BackgroundMeshWriter<void> *rasa);
template <typename T> BackgroundMeshWriter<T> restoreType(const BackgroundMeshWriter<void> &rasa);
template <typename T> BackgroundMeshReader<T> restoreType(const BackgroundMeshReader<void> *rasa);
template <typename T> BackgroundMeshReader<T> restoreType(const BackgroundMeshReader<void> &rasa);
/// \}

/// \brief Check two meshes for compatibility in combining operations.
///
/// \param target        The mesh being accumulated.  This must be a long long int-type mesh in
///                      order to enable fixed-precision accumulation.
/// \param contribution  The mesh to contribute.  Dimensions will be checked against those of the
///                      accumulating target.
template <typename Ttarget, typename Tcontrib>
void checkMeshCompatibility(const BackgroundMeshWriter<Ttarget> &target,
                            const BackgroundMeshReader<Tcontrib> &contrib);

/// \brief Evaluate an electrostatic potential which switches to a cubic softcore function at a
///        given distance.
///
/// \param abcd     Coefficients for the cubic softcore polynomial potential (unscaled by Coulomb's
///                 constant or the atom charge)
/// \param r        Distance from the atom
/// \param rswitch  Distance below which the Coulomb potential hands off to the cubic polynomial
/// \param q        Partial charge on the atom
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evalElecCSC(const Tcalc4 abcd, Tcalc r, Tcalc rswitch, Tcalc q);

/// \brief Evaluate an electrostatic potential which switches to a cubic softcore function at a
///        given distance.  Do not evaluate the derivative.  Descriptions of input parameters
///        follow from evalElecCSC() above.
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc evalElecCSCND(const Tcalc4 abcd, Tcalc r, Tcalc rswitch, Tcalc q);

/// \brief Evaluate a Lennard-Jones potential which switches to a cubic softcore function at a
///        given distance.
///
/// \param abcd     Coefficients for the cubic softcore polynomial potential
/// \param r        Distance from the atom
/// \param rswitch  Distance below which the Lennard-Jones potential hands off to the cubic
///                 polynomial
/// \param lja      Pairwise Lennard-Jones A coefficient
/// \param ljb      Pairwise Lennard-Jones B coefficient
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc2 evalLennardJonesCSC(const Tcalc4 abcd, Tcalc r, Tcalc rswitch, Tcalc lja, Tcalc ljb);

/// \brief Evaluate a Lennard-Jones potential which switches to a cubic softcore function at a
///        given distance.  Do not evaluate the derivative.  Descriptions of input parameters
///        follow from evalLennardJonesCSC() above.
template <typename Tcalc4, typename Tcalc2, typename Tcalc>
Tcalc evalLennardJonesCSCND(const Tcalc4 abcd, Tcalc r, Tcalc rswitch, Tcalc lja, Tcalc ljb);
  
#ifdef STORMM_USE_HPC
/// \brief Launch a kernel to color the occlusion mesh using the GPU.  This is, in practicality, a
///        private function of the BackgroundMesh class, but the class templating presents a
///        problem cross over to the HPC layer if a combination of an ordinary C++ compiler and
///        the HPC compiler are used.  Therefore, it appears as a free function with the template
///        types stripped from the mesh abstract and supplied instead as codes.
///
/// \param bgmw          Collection of pointers and array lengths to critical components of the
///                      background mesh, including the probe radius.  The data type of the mesh
///                      is required to be unsigned long long int.
/// \param ag            Pointer to the topology of interest (this is provided as a pointer to the
///                      original object due to the way that the proper abstract will be obtained
///                      based on the prec parameter, and the availability of the topology
///                      describing the coordinate set as a pointer in the BackgroundMesh object
///                      that produced the bgmw input parameter)
/// \param cfr           A set of atomic positions (includes the number of atoms, the trusted
///                      length of the array of radii)
/// \param radii         A collection of atomic radii
/// \param launcher      Kernel launch manager specializing in the mesh-based kernels
void launchColorOcclusionMesh(BackgroundMeshWriter<void> *bgmw, const AtomGraph *ag,
                              const CoordinateFrameReader &cfr, PrecisionModel prec,
                              const MeshKlManager &launcher);
#endif

/// \brief Accumulate a mesh in fixed precision, based on another mesh created from a single
///        structure snapshot.
///
/// Overloaded:
///   - Accept abstracts of the relevant meshes
///   - Accept the objects and, optionally, information about any available GPU
///
/// \param target        The mesh being accumulated.  This must be a long long int-type mesh in
///                      order to enable fixed-precision accumulation.
/// \param contribution  The mesh to contribute.  Dimensions will be checked against those of the
///                      accumulating target.
/// \param gpu           Details of the GPU available
/// \param availability  Information about the location of relevant mesh coefficient data.  This
///                      applies to both meshes and the accumulation will not, by itself, perform
///                      uploading or downloading of either mesh.
/// \{
void accumulateOcclusionMesh(BackgroundMeshWriter<llint> *target,
                             const BackgroundMeshReader<llint> &contribution);

void accumulateOcclusionMesh(BackgroundMesh<llint> *target,
                             const BackgroundMesh<llint> &contribution,
                             const GpuDetails &gpu = null_gpu,
                             HybridTargetLevel availability = HybridTargetLevel::HOST);
/// \}

#ifdef STORMM_USE_HPC
/// \brief Launch a kernel to perform accumulation of an occlusion mesh into an occlusion field.
///        Descriptions of input parameters follow from accumulateOcclusionMesh() above.
void launchAccOcclusionMesh(BackgroundMeshWriter<llint> *target,
                            const BackgroundMeshReader<llint> &contribution,
                            const GpuDetails &gpu);
#endif

/// \brief Accumulate a continuous field potential on a mesh.  Overloading and descriptions of
///        input parameters follow from accumulateOcclusionMesh() above, with the addition of:
///
/// \param target        The mesh that will accumulate a contribution
/// \param contribution  The mesh to be contributed
/// \param ct_contrib    Codified data type of coefficients in the contribution mesh
/// \param gpu           Specifications of the available GPU
/// \param availability  Indicate whether the coefficients data is available on the CPU host or
///                      GPU device
/// \{
template <typename T>
void accumulateNonbondedFieldMesh(BackgroundMeshWriter<llint> *target,
                                  const BackgroundMeshReader<void> &contribution,
                                  size_t ct_contrib, const GpuDetails &gpu);

template <typename T>
void accumulateNonbondedFieldMesh(BackgroundMesh<llint> *target,
                                  const BackgroundMesh<T> &contribution,
                                  const GpuDetails &gpu, HybridTargetLevel availability);
/// \}

#ifdef STORMM_USE_HPC
/// \brief Launch a templated kernel to produce a non-bonded field mesh using the avaiable GPU.
///
/// \param target        The mesh being accumulated.  Accumulation will occur in long-long integer
///                      format and then be converted to the appropriate coefficient type (must be
///                      double or float).
/// \param ct_targ       Codified data type of the mesh being constructed
/// \param mnbk          Softcore parameters for assembling the mesh
/// \param prec          Precision model in which to compute the mesh-based potentials
/// \param stn_xfrm      Transformation matrix for taking observations in each mesh element (the
///                      stencil) into tricubic polynomial coefficients
/// \param ag            Pointer to the target mesh's underlying topology
/// \param cfr           Read-only abstract for the underlying coordinates.  The coordinates must
///                      be stored on the same memory tier (CPU host or GPU device) as the topology
///                      non-bonded data.
/// \param launcher      Object containing launch parameters for the appropriate kernel
/// \param availability  Indicate whether the coefficients data is available on the CPU host or
///                      GPU device
void launchColorNonbondedFieldMesh(BackgroundMeshWriter<void> *target, size_t ct_targ,
                                   const MeshFFKit<double> &mnbk, PrecisionModel prec,
                                   const double* stn_xfrm, const AtomGraph *ag,
                                   const CoordinateFrameReader &cfr,
                                   const MeshKlManager &launcher,
                                   const HybridTargetLevel availability = HybridTargetLevel::HOST);

/// \brief Launch a templated kernel to perform accumulation on a mesh holding a non-bonded field.
///        Descriptions of parameters follow from launchColorNonbondedFieldMesh() above, in
///        addition to:
///
/// \param ct_contrib    Codified data type for the mesh being contributed
/// \param gpu           Details of the GPU available
/// \param contribution  The mesh to contribute.  Dimensions will be checked against those of the
///                      accumulating target.
void launchAccNonbondedFieldMesh(BackgroundMeshWriter<llint> *target,
                                 const BackgroundMeshReader<void> &contribution,
                                 const size_t ct_contrib, const GpuDetails &gpu);
#endif

/// \brief Compute boundary values for a patch of a non-bonded field in isolated boundary
///        conditions using natural boundary stipulations.
///
/// \param u_patch  The patch of potential values, modified and returned
/// \param dim_a    Length of the patch along the mesh A axis
/// \param dim_b    Length of the patch along the mesh B axis
/// \param dim_c    Length of the patch along the mesh C axis
/// \param start_a  The patch point at which to start editing along the mesh (and patch) A axis
/// \param start_b  The patch point at which to start editing along the mesh (and patch) B axis
/// \param start_c  The patch point at which to start editing along the mesh (and patch) C axis
/// \param end_a    The number of patch points along the mesh A axis
/// \param end_b    The number of patch points along the mesh B axis
/// \param end_c    The number of patch points along the mesh C axis
/// \param inc_a    Increment by which to step along the mesh (and patch) A axis
/// \param inc_b    Increment by which to step along the mesh (and patch) B axis
/// \param inc_c    Increment by which to step along the mesh (and patch) C axis
template <typename T>
void inferPatchBoundary(std::vector<T> *u_patch, int dim_a, int dim_b, int dim_c, int start_a,
                        int start_b, int start_c, int end_a, int end_b, int end_c, int inc_a,
                        int inc_b, int inc_c);

/// \brief Transfer the results of a batched finite difference calculation to the mesh.  This
///        function takes into account the appropriate offset to apply in setting each element's
///        derivatives and will copy the results to up to eight mesh elements.
///
/// \param patch            Array of derivatives to place on the mesh
/// \param patch_dim_a      Size of the patch array along the A axis
/// \param patch_dim_b      Size of the patch array along the B axis
/// \param patch_dim_c      Size of the patch array along the C axis
/// \param occfieldw        Writeable mesh abstract with dimensions and the coefficients array
/// \param a_index          Mesh index, along the A axis, at which to begin transcribing results
///                         from the patch
/// \param b_index          Mesh index, along the B axis, at which to begin transcribing results
///                         from the patch
/// \param c_index          Mesh index, along the C axis, at which to begin transcribing results
///                         from the patch
/// \param readout_a_start  Patch index, along the mesh A axis, at which to begin reading output to
///                         place on the mesh
/// \param readout_b_start  Patch index, along the mesh B axis, at which to begin reading output to
///                         place on the mesh
/// \param readout_c_start  Patch index, along the mesh C axis, at which to begin reading output to
///                         place on the mesh
/// \param element_offset   The offset within each element at which the derivatives encoded in the
///                         patch can be applied to elements of the mesh.  As the mesh accumulates
///                         values and derivatives, each element has 64 slots for data.  Values of
///                         the function f(a,b,c) go in slots [0-7], df/dx in slots [8-15],
///                         followed by df/dy, df/dz
template <typename T>
void transferPatchResult(const std::vector<T> &patch, int patch_dim_a, int patch_dim_b,
                         int patch_dim_c, BackgroundMeshWriter<T> *occfieldw, int a_index,
                         int b_index, int c_index, int readout_a_start, int readout_b_start,
                         int readout_c_start, int element_offset);

/// \brief Compute the derivative of a non-bonded field by numerical methods (quadratic curve
///        fitting based on the central point and its neighbors to either side).  The neighbors are
///        obtained with considerations to boundary conditions, and natural boundary conditions are
///        applied at the extremes of isolated meshes.
///
/// \param occfieldw      Writeable abstract for the mesh, containing coefficients and dimensions
/// \param a_index        Mesh index along the A axis where the patch has its origin
/// \param b_index        Mesh index along the B axis where the patch has its origin
/// \param c_index        Mesh index along the C axis where the patch has its origin
/// \param max_occlusion  The maximum value of the function, for clamping derivatives.  Set to a
///                       very large number or a negative value to avoid clamping derivatives.
template <typename T>
void patchNumericDeriv(BackgroundMeshWriter<T> *occfieldw, int a_index, int b_index, int c_index,
                       T max_occlusion);

/// \brief Compute derivatives for an occlusion mesh based on a complete set of potential values.
///        Rather than mirroring the approach for CMAP potentials, this function will automatically
///        set the derivatives to zero everywhere the occlusion is 0% or 100%.  Elsewhere, the
///        derivatives will be set based on quadratic interpolation of the neighboring points, or
///        linear interpolation (natural boundary conditions) if the point is on the mesh's
///        isolated boundary.
///
/// \param occfield      The occlusion mesh field with its potentials solved
/// \param gpu           Details of the gpu to use in computing the derivatives
/// \param availability  Indicate whether the potentials are available in the GPU device or CPU
///                      host memory
template <typename T>
void occlusionFieldDerivatives(BackgroundMesh<T> *occfield, const GpuDetails &gpu,
                               HybridTargetLevel availability);

#ifdef STORMM_USE_HPC
/// \brief Launch the templated kernel to compute derivatives on an occlusion field mesh, given a
///        void-casted abstract and the code for the original type.
///
/// \param bgmw     Writeable abstract for the background mesh, with its coefficients pointer cast
///                 to void to remove templated character
/// \param ct_mesh  Codified data type of the coefficients in the background mesh
void launchOccFieldDerivativeCalc(BackgroundMeshWriter<void> *bgmw, size_t ct_mesh);
#endif
  
} // namespace structure
} // namespace stormm

#include "background_mesh.tpp"
#include "background_mesh_freefunc.tpp"

#endif
