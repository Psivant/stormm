// -*-c++-*-
#include "copyright.h"
#include "DataTypes/common_types.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tdata>
BackgroundMeshWriter<Tdata>::BackgroundMeshWriter(const MeshParamKit &dims_in, GridDetail kind_in,
                                                  NonbondedPotential field_in,
                                                  const MeshRulerKit &rulers, Tdata* coeffs_in,
                                                  double coeff_scale_in, double probe_radius_in,
                                                  double well_depth_in, double occ_cost_in,
                                                  const MeshBasicsKit &mbss_in) :
    dims{dims_in}, kind{kind_in}, field{field_in}, rlrs{rulers}, coeffs{coeffs_in},
    coeff_scale{coeff_scale_in}, coeff_scale_f{static_cast<float>(coeff_scale_in)},
    probe_radius{probe_radius_in}, well_depth{well_depth_in}, occ_cost{occ_cost_in}, mbss{mbss_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename Tdata>
BackgroundMeshReader<Tdata>::BackgroundMeshReader(const MeshParamKit &dims_in, GridDetail kind_in,
                                                  NonbondedPotential field_in,
                                                  const MeshRulerKit &rulers,
                                                  const Tdata* coeffs_in, double coeff_scale_in,
                                                  double probe_radius_in, double well_depth_in,
                                                  double occ_cost_in,
                                                  const MeshBasicsKit &mbss_in) :
    dims{dims_in}, kind{kind_in}, field{field_in}, rlrs{rulers}, coeffs{coeffs_in},
    coeff_scale{coeff_scale_in}, coeff_scale_f{static_cast<float>(coeff_scale_in)},
    probe_radius{probe_radius_in}, well_depth{well_depth_in}, occ_cost{occ_cost_in}, mbss{mbss_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const MeshParameters &measurements_in,
                                  const PrecisionModel build_precision_in) :
    measurements{measurements_in}, kind{kind_in}, field{field_in}, probe_radius{0.0},
    well_depth{0.0}, occlusion_penalty{1.0},
    build_precision{build_precision_in},
    tick_marks{measurements_in},
    nonbonded_model{},
    coefficients{HybridKind::ARRAY, "mesh_tricubic_coef"},
    coefficient_scale_bits{40},
    coefficient_scale{pow(2.0, coefficient_scale_bits)},
    inverse_coefficient_scale{1.0 / coefficient_scale},
    basis{}
{
  validateMeshKind();
  allocate();
}
  
//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const MeshParameters &measurements_in,
                                  const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon,
                                  const double clash_distance_in, const double clash_ratio_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
  BackgroundMesh(kind_in, field_in, measurements_in, prec)
{
  // Set the system
  basis = MeshFoundation(cf_in, ag_in);
  
  // Validate the mesh dimensions, then allocate memory
  validateScalingBits();
  allocate();

  // Map the potential and atomic near-neighbor interactions
  setProbeRadius(probe_radius_in);
  setWellDepth(well_depth_in);
  setNonbondedModel(mixing_protocol_in, clash_ratio_in, clash_distance_in, probe_sigma,
                    probe_epsilon);
  computeField(launcher, prec, availability, probe_sigma, probe_epsilon);
  switch (kind) {
  case GridDetail::OCCLUSION:
  case GridDetail::OCCLUSION_FIELD:
  case GridDetail::NONBONDED_FIELD:
    break;
  case GridDetail::NONBONDED_ATOMIC:
    basis.computeNeighborLists(measurements, tick_marks, launcher, prec, availability);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const MeshParameters &measurements_in,
                                  const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon,
                                  const double clash_distance_in, const double clash_ratio_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher, 
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in.getSelfPointer(), cf_in.getSelfPointer(),
                   probe_radius_in, well_depth_in, mixing_protocol_in, measurements_in,
                   probe_sigma, probe_epsilon, clash_distance_in, clash_ratio_in, prec, launcher,
                   availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const double buffer, const double spacing,
                                  const int scale_bits_in, const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon,
                                  const double clash_distance_in, const double clash_ratio_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in, cf_in, probe_radius_in, well_depth_in,
                   mixing_protocol_in, getMeasurements(ag_in, cf_in, buffer,
                                                       std::vector<double>(3, spacing),
                                                       scale_bits_in), probe_sigma, probe_epsilon,
                   clash_distance_in, clash_ratio_in, prec, launcher, availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                                  const double buffer, const double spacing,
                                  const int scale_bits_in, const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon,
                                  const double clash_distance_in, const double clash_ratio_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, probe_radius_in, well_depth_in, mixing_protocol_in,
                   ag_in.getSelfPointer(), cf_in.getSelfPointer(), buffer, spacing, scale_bits_in,
                   probe_sigma, probe_epsilon, clash_distance_in, clash_ratio_in, prec, launcher,
                   availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const double buffer, const std::vector<double> &spacing,
                                  const int scale_bits_in, const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon,
                                  const double clash_distance_in, const double clash_ratio_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in, cf_in, probe_radius_in, well_depth_in,
                   mixing_protocol_in, getMeasurements(ag_in, cf_in, buffer, spacing,
                                                       scale_bits_in), probe_sigma, probe_epsilon,
                   clash_distance_in, clash_ratio_in, prec, launcher, availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                                  const double buffer, const std::vector<double> &spacing,
                                  const int scale_bits_in, const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon,
                                  const double clash_distance_in, const double clash_ratio_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, probe_radius_in, well_depth_in, mixing_protocol_in,
                   ag_in.getSelfPointer(), cf_in.getSelfPointer(), buffer, spacing, scale_bits_in,
                   probe_sigma, probe_epsilon, clash_distance_in, clash_ratio_in, prec, launcher,
                   availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const std::vector<double> &mesh_bounds, const double spacing,
                                  const int scale_bits_in, const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon,
                                  const double clash_distance_in, const double clash_ratio_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in, cf_in, probe_radius_in, well_depth_in,
                   mixing_protocol_in, getMeasurements(mesh_bounds,
                                                       std::vector<double>(3, spacing),
                                                       scale_bits_in), probe_sigma, probe_epsilon,
                   clash_distance_in, clash_ratio_in, prec, launcher, availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                                  const std::vector<double> &mesh_bounds, const double spacing,
                                  const int scale_bits_in, const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon,
                                  const double clash_distance_in, const double clash_ratio_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, probe_radius_in, well_depth_in, mixing_protocol_in,
                   ag_in.getSelfPointer(), cf_in.getSelfPointer(), mesh_bounds, spacing,
                   scale_bits_in, probe_sigma, probe_epsilon, clash_distance_in, clash_ratio_in,
                   launcher, availability, prec)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon,
                                  const double clash_distance_in, const double clash_ratio_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in, cf_in, probe_radius_in, well_depth_in,
                   mixing_protocol_in, getMeasurements(mesh_bounds, spacing, scale_bits_in),
                   probe_sigma, probe_epsilon, clash_distance_in, clash_ratio_in, prec, launcher,
                   availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon,
                                  const double clash_distance_in, const double clash_ratio_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, probe_radius_in, well_depth_in, mixing_protocol_in,
                   ag_in.getSelfPointer(), cf_in.getSelfPointer(), mesh_bounds, spacing,
                   scale_bits_in, probe_sigma, probe_epsilon, clash_distance_in, clash_ratio_in,
                   prec, launcher, availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const double buffer, const double spacing,
                                  const int scale_bits_in, const double clash_distance_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
  BackgroundMesh(kind_in, field_in, ag_in, cf_in, 0.0, 0.0, VdwCombiningRule::LORENTZ_BERTHELOT,
                 getMeasurements(ag_in, cf_in, buffer, std::vector<double>(3, spacing),
                                 scale_bits_in), {}, {}, clash_distance_in,
                 default_mesh_vdw_damping_ratio, prec, launcher, availability)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                                  const double buffer, const double spacing,
                                  const int scale_bits_in, const double clash_distance_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in,  ag_in.getSelfPointer(), cf_in.getSelfPointer(), buffer,
                   spacing, scale_bits_in, clash_distance_in, prec, launcher, availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const double buffer, const std::vector<double> &spacing,
                                  const int scale_bits_in, const double clash_distance_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in, cf_in, 0.0, 0.0, VdwCombiningRule::LORENTZ_BERTHELOT,
                   getMeasurements(ag_in, cf_in, buffer, spacing, scale_bits_in), {}, {},
                   clash_distance_in, default_mesh_vdw_damping_ratio, prec, launcher,
                   availability)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                                  const double buffer, const std::vector<double> &spacing,
                                  const int scale_bits_in, const double clash_distance_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in.getSelfPointer(), cf_in.getSelfPointer(), buffer,
                   spacing, scale_bits_in, clash_distance_in, prec, launcher, availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const std::vector<double> &mesh_bounds, const double spacing,
                                  const int scale_bits_in, const double clash_distance_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in, cf_in, 0.0, 0.0, VdwCombiningRule::LORENTZ_BERTHELOT,
                   getMeasurements(mesh_bounds, std::vector<double>(3, spacing), scale_bits_in),
                   {}, {}, clash_distance_in, default_mesh_vdw_damping_ratio, prec, launcher,
                   availability)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                                  const std::vector<double> &mesh_bounds, const double spacing,
                                  const int scale_bits_in, const double clash_distance_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in.getSelfPointer(), cf_in.getSelfPointer(), mesh_bounds,
                   spacing, scale_bits_in, clash_distance_in, prec, launcher, availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const double clash_distance_in, const PrecisionModel prec,
                                  const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in, cf_in, 0.0, 0.0, VdwCombiningRule::LORENTZ_BERTHELOT,
                   getMeasurements(mesh_bounds, spacing, scale_bits_in), {}, {}, clash_distance_in,
                   default_mesh_vdw_damping_ratio, prec, launcher, availability)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const double clash_distance_in, const PrecisionModel prec,
                                  const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in.getSelfPointer(), cf_in.getSelfPointer(), mesh_bounds,
                   spacing, scale_bits_in, clash_distance_in, prec, launcher, availability)
{}
  
//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const double buffer, const double spacing,
                                  const int scale_bits_in, const std::vector<double> &probe_sigma,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, NonbondedPotential::CLASH, ag_in, cf_in, probe_radius_in, 0.0,
                   mixing_protocol_in,
                   getMeasurements(ag_in, cf_in, buffer, std::vector<double>(3, spacing),
                                   scale_bits_in), probe_sigma, {},
                   default_mesh_elec_damping_range, default_mesh_vdw_damping_ratio, prec, launcher,
                   availability)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                                  const double buffer, const double spacing,
                                  const int scale_bits_in, const std::vector<double> &probe_sigma,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, probe_radius_in, mixing_protocol_in, ag_in.getSelfPointer(),
                   cf_in.getSelfPointer(), buffer, spacing, scale_bits_in, probe_sigma, {},
                   prec, launcher, availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const double buffer, const std::vector<double> &spacing,
                                  const int scale_bits_in, const std::vector<double> &probe_sigma,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, NonbondedPotential::CLASH, ag_in, cf_in, probe_radius_in, 0.0,
                   mixing_protocol_in,
                   getMeasurements(ag_in, cf_in, buffer, spacing, scale_bits_in), probe_sigma, {},
                   default_mesh_elec_damping_range, default_mesh_vdw_damping_ratio, prec, launcher,
                   availability)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                                  const double buffer, const std::vector<double> &spacing,
                                  const int scale_bits_in, const std::vector<double> &probe_sigma,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, probe_radius_in, mixing_protocol_in, ag_in.getSelfPointer(),
                   cf_in.getSelfPointer(), buffer, spacing, scale_bits_in, probe_sigma, {},
                   prec, launcher, availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const std::vector<double> &mesh_bounds,
                                  const double spacing, const int scale_bits_in,
                                  const std::vector<double> &probe_sigma,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, NonbondedPotential::CLASH, ag_in, cf_in, probe_radius_in, 0.0,
                   mixing_protocol_in,
                   getMeasurements(mesh_bounds, std::vector<double>(3, spacing), scale_bits_in),
                   probe_sigma, {}, default_mesh_elec_damping_range,
                   default_mesh_vdw_damping_ratio, prec, launcher, availability)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                                  const std::vector<double> &mesh_bounds,
                                  const double spacing, const int scale_bits_in,
                                  const std::vector<double> &probe_sigma,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, probe_radius_in, mixing_protocol_in, ag_in.getSelfPointer(),
                   cf_in.getSelfPointer(), mesh_bounds, spacing, scale_bits_in, probe_sigma,
                   prec, launcher, availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const std::vector<double> &probe_sigma,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, NonbondedPotential::CLASH, ag_in, cf_in, probe_radius_in, 0.0,
                   mixing_protocol_in, getMeasurements(mesh_bounds, spacing, scale_bits_in),
                   probe_sigma, {}, default_mesh_elec_damping_range,
                   default_mesh_vdw_damping_ratio, prec, launcher, availability)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph &ag_in, const CoordinateFrame &cf_in,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const std::vector<double> &probe_sigma,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, probe_radius_in, mixing_protocol_in, ag_in.getSelfPointer(),
                   cf_in.getSelfPointer(), mesh_bounds, spacing, scale_bits_in, probe_sigma, prec,
                   launcher, availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tcoord>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph &ag_in, const CoordinateSeries<Tcoord> &cs_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const MeshParameters &measurements_in, const int averaging_order,
                                  const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon,
                                  const double clash_distance_in, const double clash_ratio_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
  BackgroundMesh(kind_in, field_in, measurements_in, prec)
{
  // Set the system
  basis = MeshFoundation(cs_in, ag_in);

  // Devise new mesh parameters based on the frame averaging at hand.  Occlusion meshes will
  // require different parameters, as their bitwise densities are averaged in a
  MeshParameters snapshot_measurements = measurements_in;
  const MeshParamKit measr = measurements_in.data();
  double a_vector_rescale, b_vector_rescale, c_vector_rescale;
  switch (kind_in) {
  case GridDetail::NONBONDED_FIELD:
    a_vector_rescale = 1.0;
    b_vector_rescale = 1.0;
    c_vector_rescale = 1.0;
    break;
  case GridDetail::OCCLUSION_FIELD:
    {
      // The "top-hat" function density computation will take place on a mesh of size and
      // discretization related to the continuous field mesh, element for element.  If the
      // continuous density mesh derives from the series of top-hat meshes with an averaging order
      // of 4, then each sub-element (1/64th of a full element) of the top-hat mesh corresponds to
      // an element of the continuous mesh, and the top-hat mesh element will be scaled up 16x.
      // At the other end of the spectrum, if the continuous density mesh derives from the series
      // of top-hat meshes with an averaging order of 16, the elements of the top-hat occlusion
      // meshes will have the same dimensions as the continuous density field mesh.  In all cases,
      // the top-hat mesh origin will be shifted by half of the element width of the continuous
      // mesh along all three lattice vectors.  If the mesh spans a periodic system, the capacity
      // to accept lower sampling orders than 16 will depend on the continuous mesh element counts
      // having enough factors of two along each dimension.
      switch (averaging_order) {
      case 4:

        // Check that the dimensions of the continuous mesh are divisible by four.  If not, check
        // for a factor of 2, and if that is not found, do not interpolate along that dimension.
        a_vector_rescale = ((measr.na % 4) == 0) ? 4.0 : ((measr.na % 2) == 0) ? 2.0 : 1.0;
        b_vector_rescale = ((measr.nb % 4) == 0) ? 4.0 : ((measr.nb % 2) == 0) ? 2.0 : 1.0;
        c_vector_rescale = ((measr.nc % 4) == 0) ? 4.0 : ((measr.nc % 2) == 0) ? 2.0 : 1.0;
        break;
      case 8:

        // Check that the dimensions of the continuous mesh are divisible by two.  If not, do not
        // interpolate along that dimension.
        a_vector_rescale = ((measr.na % 2) == 0) ? 2.0 : 1.0;
        b_vector_rescale = ((measr.nb % 2) == 0) ? 2.0 : 1.0;
        c_vector_rescale = ((measr.nc % 2) == 0) ? 2.0 : 1.0;
        break;
      case 32:
        a_vector_rescale = 0.5;
        b_vector_rescale = 0.5;
        c_vector_rescale = 0.5;
        break;
      case 16:
      default:
        a_vector_rescale = 1.0;
        b_vector_rescale = 1.0;
        c_vector_rescale = 1.0;
        break;
      }

      // Redefine the mesh for this snapshot
      const int na_sm = snapshot_measurements.getAxisElementCount(UnitCellAxis::A);
      const int nb_sm = snapshot_measurements.getAxisElementCount(UnitCellAxis::B);
      const int nc_sm = snapshot_measurements.getAxisElementCount(UnitCellAxis::C);
      std::vector<double> invu = snapshot_measurements.getMeshInverseTransform<double>();
      int period_factor;
      switch (measurements_in.getBoundaryConditions()) {
      case BoundaryCondition::PERIODIC:
        period_factor = 0;
        break;
      case BoundaryCondition::ISOLATED:
        period_factor = 1;
        break;
      }
      const int na_update = round(static_cast<double>(na_sm) / a_vector_rescale) + period_factor;
      const int nb_update = round(static_cast<double>(nb_sm) / b_vector_rescale) + period_factor;
      const int nc_update = round(static_cast<double>(nc_sm) / c_vector_rescale) + period_factor;
      for (int i = 0; i < 3; i++) {
        invu[i    ] *= a_vector_rescale;
        invu[i + 3] *= b_vector_rescale;
        invu[i + 6] *= c_vector_rescale;
      }
      double orig_update_x = snapshot_measurements.getMeshOrigin(CartesianDimension::X);
      double orig_update_y = snapshot_measurements.getMeshOrigin(CartesianDimension::Y);
      double orig_update_z = snapshot_measurements.getMeshOrigin(CartesianDimension::Z);
      orig_update_x -= (0.5 / a_vector_rescale) * (invu[0] + invu[3] + invu[6]);
      orig_update_y -= (0.5 / b_vector_rescale) * (invu[1] + invu[4] + invu[7]);
      orig_update_z -= (0.5 / c_vector_rescale) * (invu[2] + invu[5] + invu[8]);
      snapshot_measurements = MeshParameters(na_update, nb_update, nc_update, orig_update_x,
                                             orig_update_y, orig_update_z, invu,
                                             snapshot_measurements.getScalingBits());
    }
    break;
  case GridDetail::OCCLUSION:
  case GridDetail::NONBONDED_ATOMIC:
    break;
  }

  // Formulate a mesh with fixed-precision values in which to perform the accumulation.  Even for a
  // continuous field computed on a single thread, accumulation in fixed precision is important to
  // make the final mesh-based potential independent of the order in which snapshots are presented.
  // This is also an opportunity to check whether any continuous field mesh might exceed the
  // available format.
  switch (kind_in) {
  case GridDetail::NONBONDED_FIELD:
    switch (field_in) {
    case NonbondedPotential::ELECTROSTATIC:
    case NonbondedPotential::VAN_DER_WAALS:
    case NonbondedPotential::CLASH:
      break;
    }
    break;
  case GridDetail::OCCLUSION_FIELD:
  case GridDetail::OCCLUSION:
  case GridDetail::NONBONDED_ATOMIC:
    break;
  }
  BackgroundMesh<llint> proto(kind_in, field_in, measurements_in);
  proto.setCoefficientScalingBits(40);
  
  // Loop over each frame of the series, create the appropriate mesh for the accumulated
  // potential, and sum the coefficients of the actual mesh.
  const int nframe = cs_in.getFrameCount();
  CoordinateFrame cf_hold(cs_in.getAtomCount());
  const GpuDetails& gpu = launcher.getGpu();
#ifdef STORMM_USE_HPC
  const HybridTargetLevel target_tier = (gpu == null_gpu) ? HybridTargetLevel::HOST :
                                                            HybridTargetLevel::DEVICE;
#else
  const HybridTargetLevel target_tier = HybridTargetLevel::HOST;
#endif
  BackgroundMeshWriter<llint> protow = proto.data(target_tier);
  for (int i = 0; i < nframe; i++) {
    switch (kind_in) {
    case GridDetail::OCCLUSION:
    case GridDetail::NONBONDED_ATOMIC:
      rtErr("A mesh of type " + getEnumerationName(kind_in) + " must be created using a single "
            "structure provided through a CoordinateFrame object, not a series of " +
            std::to_string(nframe) + " frames.", "BackgroundMesh");
    case GridDetail::OCCLUSION_FIELD:
    case GridDetail::NONBONDED_FIELD:
      {
        // Copy coordinates of the ith structure from the tier at which they are available tier
        // in the series to the appropriate tier of the holding frame.

        // The coordinate copy operation will try to raise the availability of the coordinates to
        // the DEVICE tier if a GPU is available.
        coordCopy<Tcoord>(&cf_hold, cs_in, i, target_tier, availability, gpu);
        if (kind_in == GridDetail::OCCLUSION_FIELD) {
          const BackgroundMesh<llint> t_mesh(GridDetail::OCCLUSION, NonbondedPotential::CLASH,
                                             basis.getTopologyPointer(), &cf_hold, probe_radius_in,
                                             well_depth_in, mixing_protocol_in,
                                             snapshot_measurements, probe_sigma, probe_epsilon,
                                             clash_distance_in, clash_ratio_in, prec, launcher,
                                             target_tier);
          accumulateOcclusionMesh(&proto, t_mesh, gpu, target_tier);
        }
        else {
          const BackgroundMesh<T> t_mesh(kind_in, field_in, basis.getTopologyPointer(), cf_hold,
                                         probe_radius_in, well_depth_in, mixing_protocol_in,
                                         snapshot_measurements, probe_sigma, probe_epsilon, prec,
                                         launcher, target_tier);
          const size_t ct_self = std::type_index(typeid(T)).hash_code();
          accumulateNonbondedFieldMesh(&proto, t_mesh, gpu, target_tier);
        }
      }
      break;
    }
  }
  
  // Validate the mesh dimensions, then allocate memory
  validateScalingBits();
  allocate();

  // Map the potential and atomic near-neighbor interactions for the ensemble-averaged mesh at
  // hand, as will have been done for each of the contributing meshes.
  setProbeRadius(probe_radius_in);
  setWellDepth(well_depth_in);
  setNonbondedModel(mixing_protocol_in, clash_ratio_in, clash_distance_in, probe_sigma,
                    probe_epsilon);
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tcoord>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph &ag_in, const CoordinateSeries<Tcoord> &cs_in,
                                  const double buffer, const double spacing,
                                  const int scale_bits_in, const int averaging_order,
                                  const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon,
                                  const double clash_distance_in, const double clash_ratio_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in, cs_in, probe_radius_in, well_depth_in,
                   mixing_protocol_in,
                   getMeasurements(ag_in, cs_in, buffer, std::vector<double>(3, spacing),
                                   scale_bits_in), averaging_order, probe_sigma, probe_epsilon,
                   clash_distance_in, clash_ratio_in, prec, launcher, availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tcoord>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph &ag_in, const CoordinateSeries<Tcoord> &cs_in,
                                  const double buffer, const std::vector<double> &spacing,
                                  const int scale_bits_in, const int averaging_order,
                                  const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon,
                                  const double clash_distance_in, const double clash_ratio_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in, cs_in, probe_radius_in, well_depth_in,
                   mixing_protocol_in, getMeasurements(ag_in, cs_in, buffer, spacing,
                                                       scale_bits_in), averaging_order,
                   probe_sigma, probe_epsilon, clash_distance_in, clash_ratio_in, prec, launcher,
                   availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tcoord>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph &ag_in, const CoordinateSeries<Tcoord> &cs_in,
                                  const std::vector<double> &mesh_bounds, const double spacing,
                                  const int scale_bits_in, const int averaging_order,
                                  const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon,
                                  const double clash_distance_in, const double clash_ratio_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in, cs_in, probe_radius_in, well_depth_in,
                   mixing_protocol_in,
                   getMeasurements(mesh_bounds, std::vector<double>(3, spacing), scale_bits_in),
                   averaging_order, probe_sigma, probe_epsilon, clash_distance_in, clash_ratio_in,
                   prec, launcher, availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tcoord>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const double probe_radius_in, const double well_depth_in,
                                  const VdwCombiningRule mixing_protocol_in,
                                  const AtomGraph &ag_in, const CoordinateSeries<Tcoord> &cs_in,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const int averaging_order,
                                  const std::vector<double> &probe_sigma,
                                  const std::vector<double> &probe_epsilon,
                                  const double clash_distance_in, const double clash_ratio_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in, cs_in, probe_radius_in, well_depth_in,
                   mixing_protocol_in,
                   getMeasurements(mesh_bounds, spacing, scale_bits_in), averaging_order,
                   probe_sigma, probe_epsilon, clash_distance_in, clash_ratio_in, prec, launcher,
                   availability)
{}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tcoord>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph &ag_in, const CoordinateSeries<Tcoord> &cs_in,
                                  const double buffer, const double spacing,
                                  const int scale_bits_in, const double clash_distance_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in,  ag_in, cs_in, 0.0, 0.0, VdwCombiningRule::LORENTZ_BERTHELOT,
                   getMeasurements(ag_in, cs_in, buffer, std::vector<double>(3, spacing),
                                   scale_bits_in), 1, {}, {}, clash_distance_in, 1.0, prec,
                   launcher, availability)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tcoord>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph &ag_in, const CoordinateSeries<Tcoord> &cs_in,
                                  const double buffer, const std::vector<double> &spacing,
                                  const int scale_bits_in, const double clash_distance_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
  BackgroundMesh(kind_in, field_in, ag_in, cs_in, 0.0, 0.0, VdwCombiningRule::LORENTZ_BERTHELOT,
                 getMeasurements(ag_in, cs_in, buffer, spacing, scale_bits_in), 1, {}, {},
                 clash_distance_in, 1.0, prec, launcher, availability)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tcoord>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph &ag_in, const CoordinateSeries<Tcoord> &cs_in,
                                  const std::vector<double> &mesh_bounds, const double spacing,
                                  const int scale_bits_in, const double clash_distance_in,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in, cs_in, 0.0, 0.0, VdwCombiningRule::LORENTZ_BERTHELOT,
                   getMeasurements(mesh_bounds, std::vector<double>(3, spacing), scale_bits_in),
                   1, {}, {}, clash_distance_in, 1.0, prec, launcher, availability)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tcoord>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const NonbondedPotential field_in,
                                  const AtomGraph &ag_in, const CoordinateSeries<Tcoord> &cs_in,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const double clash_distance_in, const PrecisionModel prec,
                                  const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
    BackgroundMesh(kind_in, field_in, ag_in, cs_in, 0.0, 0.0, VdwCombiningRule::LORENTZ_BERTHELOT,
                   getMeasurements(mesh_bounds, spacing, scale_bits_in), 1, {}, {},
                   clash_distance_in, 1.0, prec, launcher, availability)
{
  // This overload provides a way to get electrostatic potential meshes without specifying an
  // irrelevant van-der Waals type or clash probe radius.
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tcoord>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const AtomGraph &ag_in, const CoordinateSeries<Tcoord> &cs_in,
                                  const double buffer, const double spacing,
                                  const int scale_bits_in, const int averaging_order,
                                  const std::vector<double> &probe_sigma,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, ag_in, cs_in, probe_radius_in, 0.0,
                 VdwCombiningRule::LORENTZ_BERTHELOT,
                 getMeasurements(ag_in, cs_in, buffer, std::vector<double>(3, spacing),
                                 scale_bits_in), averaging_order, probe_sigma, {}, 1.0, 1.0, prec,
                 launcher, availability)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tcoord>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const AtomGraph &ag_in, const CoordinateSeries<Tcoord> &cs_in,
                                  const double buffer, const std::vector<double> &spacing,
                                  const int scale_bits_in, const int averaging_order,
                                  const std::vector<double> &probe_sigma,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, ag_in, cs_in, probe_radius_in, 0.0,
                 VdwCombiningRule::LORENTZ_BERTHELOT,
                 getMeasurements(ag_in, cs_in, buffer, spacing, scale_bits_in), averaging_order,
                 probe_sigma, {}, 1.0, 1.0, prec, launcher, availability)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tcoord>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const AtomGraph &ag_in, const CoordinateSeries<Tcoord> &cs_in,
                                  const std::vector<double> &mesh_bounds,
                                  const double spacing, const int scale_bits_in,
                                  const int averaging_order,
                                  const std::vector<double> &probe_sigma,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, ag_in, cs_in, probe_radius_in, 0.0,
                 VdwCombiningRule::LORENTZ_BERTHELOT,
                 getMeasurements(mesh_bounds, std::vector<double>(3, spacing), scale_bits_in),
                 averaging_order, probe_sigma, {}, 1.0, 1.0, prec, launcher, availability)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tcoord>
BackgroundMesh<T>::BackgroundMesh(const GridDetail kind_in, const double probe_radius_in,
                                  const AtomGraph &ag_in, const CoordinateSeries<Tcoord> &cs_in,
                                  const std::vector<double> &mesh_bounds,
                                  const std::vector<double> &spacing, const int scale_bits_in,
                                  const int averaging_order,
                                  const std::vector<double> &probe_sigma,
                                  const PrecisionModel prec, const MeshKlManager &launcher,
                                  const HybridTargetLevel availability) :
  BackgroundMesh(kind_in, NonbondedPotential::CLASH, ag_in, cs_in, probe_radius_in, 0.0,
                 VdwCombiningRule::LORENTZ_BERTHELOT,
                 getMeasurements(mesh_bounds, spacing, scale_bits_in), averaging_order,
                 probe_sigma, {}, 1.0, 1.0, prec, launcher, availability)
{
  // This overload provides a way to get an occlusion mask without specifying an irrelevant van-der
  // Waals atom type or the type of non-bonded field.
}

//-------------------------------------------------------------------------------------------------
template <typename T> const MeshParameters& BackgroundMesh<T>::getDimensions() const {
  return measurements;
}

//-------------------------------------------------------------------------------------------------
template <typename T> const MeshFoundation& BackgroundMesh<T>::getMolecularBasis() const {
  return basis;
}

//-------------------------------------------------------------------------------------------------
template <typename T> MeshFoundation* BackgroundMesh<T>::getMolecularBasis() {
  return &basis;
}

//-------------------------------------------------------------------------------------------------
template <typename T> const AtomGraph* BackgroundMesh<T>::getTopologyPointer() const {
  return basis.getTopologyPointer();
}

//-------------------------------------------------------------------------------------------------
template <typename T> const CoordinateFrame* BackgroundMesh<T>::getCoordinatePointer() const {
  return basis.getCoordinatePointer();
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tcoord>
const CoordinateSeries<Tcoord>* BackgroundMesh<T>::getEnsemblePointer() const {
  return basis.getEnsemblePointer<Tcoord>();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
size_t BackgroundMesh<T>::getEnsembleTypeCode() const {
  return basis.getEnsembleTypeCode();
}

//-------------------------------------------------------------------------------------------------
template <typename T> GridDetail BackgroundMesh<T>::getMeshKind() const {
  return kind;
}

//-------------------------------------------------------------------------------------------------
template <typename T> NonbondedPotential BackgroundMesh<T>::getNonbondedPotential() const {
  switch (kind) {
  case GridDetail::OCCLUSION:
  case GridDetail::OCCLUSION_FIELD:
    rtErr("This mesh is associated with occlusion (type " + getEnumerationName(kind) + "), not a "
          "non-bonded potential.", "BackgroundMesh", "getNonbondedPotential");
  case GridDetail::NONBONDED_FIELD:
  case GridDetail::NONBONDED_ATOMIC:
    break;
  }
  return field;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double BackgroundMesh<T>::getProbeRadius() const {
  return probe_radius;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double BackgroundMesh<T>::getWellDepth() const {
  switch (kind) {
  case GridDetail::OCCLUSION:
  case GridDetail::OCCLUSION_FIELD:
    rtErr("This mesh is associated with occusion (type " + getEnumerationName(kind) + "), not a "
          "non-bonded potential.", "BackgroundMesh", "getWellDepth");
  case GridDetail::NONBONDED_FIELD:
  case GridDetail::NONBONDED_ATOMIC:
    switch (field) {
    case NonbondedPotential::ELECTROSTATIC:
    case NonbondedPotential::CLASH:
      rtErr("Meshes associated with a " + getEnumerationName(field) + " non-bonded potential "
            "do not have a valid well depth to report.", "BackgroundMesh", "getWellDepth");
    case NonbondedPotential::VAN_DER_WAALS:
      break;
    }
    break;
  }
  return well_depth;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double BackgroundMesh<T>::getClashDistance() const {
  return nonbonded_model.getClashDistance();
}

//-------------------------------------------------------------------------------------------------
template <typename T> double BackgroundMesh<T>::getClashRatio() const {
  return nonbonded_model.getClashRatio();
}

//-------------------------------------------------------------------------------------------------
template <typename T> double BackgroundMesh<T>::getOcclusionPenalty() const {
  return occlusion_penalty;
}

//-------------------------------------------------------------------------------------------------
template <typename T> VdwCombiningRule BackgroundMesh<T>::getCombiningRule() const {
  switch (kind) {
  case GridDetail::OCCLUSION:
  case GridDetail::OCCLUSION_FIELD:
    rtErr("This mesh is associated with occusion (type " + getEnumerationName(kind) + "), not a "
          "non-bonded potential.", "BackgroundMesh", "getCombiningRule");
  case GridDetail::NONBONDED_FIELD:
  case GridDetail::NONBONDED_ATOMIC:
    switch (field) {
    case NonbondedPotential::ELECTROSTATIC:
    case NonbondedPotential::CLASH:
      rtErr("Meshes associated with a " + getEnumerationName(field) + " non-bonded potential "
            "do not have a valid mixing protocol to report.", "BackgroundMesh",
            "getCombiningRule");
    case NonbondedPotential::VAN_DER_WAALS:
      break;
    }
    break;
  }
  return nonbonded_model.getCombiningRule();
}

//-------------------------------------------------------------------------------------------------
template <typename T> PrecisionModel BackgroundMesh<T>::getBuildPrecision() const {
  return build_precision;
}

//-------------------------------------------------------------------------------------------------
template <typename T> double BackgroundMesh<T>::getCoefficientScalingFactor() const {
  return coefficient_scale;
}

//-------------------------------------------------------------------------------------------------
template <typename T> int BackgroundMesh<T>::getCoefficientScalingBits() const {
  return coefficient_scale_bits;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const BackgroundMeshReader<T> BackgroundMesh<T>::data(const HybridTargetLevel tier) const {
  return BackgroundMeshReader<T>(measurements.data(), kind, field, tick_marks.data(tier),
                                 coefficients.data(tier), coefficient_scale, probe_radius,
                                 well_depth, occlusion_penalty, basis.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMeshWriter<T> BackgroundMesh<T>::data(const HybridTargetLevel tier) {
  return BackgroundMeshWriter<T>(measurements.data(), kind, field, tick_marks.data(tier),
                                 coefficients.data(tier), coefficient_scale, probe_radius,
                                 well_depth, occlusion_penalty, basis.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T> const BackgroundMeshReader<void>
BackgroundMesh<T>::templateFreeData(const HybridTargetLevel tier) const {
  return BackgroundMeshReader<void>(measurements.data(), kind, field, tick_marks.data(tier),
                                    reinterpret_cast<void*>(coefficients.data(tier)),
                                    coefficient_scale, probe_radius, well_depth, occlusion_penalty,
                                    basis.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMeshWriter<void> BackgroundMesh<T>::templateFreeData(const HybridTargetLevel tier) {
  return BackgroundMeshWriter<void>(measurements.data(), kind, field, tick_marks.data(tier),
                                    reinterpret_cast<void*>(coefficients.data(tier)),
                                    coefficient_scale, probe_radius, well_depth, occlusion_penalty,
                                    basis.data(tier));
}

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMeshWriter<void> BackgroundMesh<T>::deviceViewToTemplateFreeHostData() {
  void* coeff_ptr = reinterpret_cast<void*>(coefficients.getDeviceValidHostPointer());
  return BackgroundMeshWriter<void>(measurements.data(), kind, field,
                                    tick_marks.deviceViewToHostData(), coeff_ptr,
                                    coefficient_scale, probe_radius, well_depth,
                                    basis.deviceViewToHostData());
}

//-------------------------------------------------------------------------------------------------
template <typename T> const BackgroundMeshReader<void>
BackgroundMesh<T>::deviceViewToTemplateFreeHostData() const {
  const void* coeff_ptr = reinterpret_cast<const void*>(coefficients.getDeviceValidHostPointer());
  return BackgroundMeshReader<void>(measurements.data(), kind, field,
                                    tick_marks.deviceViewToHostData(), coeff_ptr,
                                    coefficient_scale, probe_radius, well_depth,
                                    basis.deviceViewToHostData());
}
#  endif
#endif

//-------------------------------------------------------------------------------------------------
template <typename T> const MeshFFKit<T>
BackgroundMesh<T>::getNonbondedKit(const HybridTargetLevel tier) const {
  return nonbonded_model.data(tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T> const MeshFFKit<double>
BackgroundMesh<T>::getReferenceNonbondedKit(const HybridTargetLevel tier) const {
  return nonbonded_model.referenceData(tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T> const MeshFFKit<void>
BackgroundMesh<T>::templateFreeNonbondedKit(const HybridTargetLevel tier) const {
  return nonbonded_model.templateFreeData(tier);
}

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
//-------------------------------------------------------------------------------------------------
template <typename T> const MeshFFKit<double>
BackgroundMesh<T>::deviceViewToReferenceNonbondedKit() const {
  return nonbonded_model.deviceViewToReferenceHostData();
}

//-------------------------------------------------------------------------------------------------
template <typename T> const MeshFFKit<void>
BackgroundMesh<T>::deviceViewToTemplateFreeNonbondedKit() const {
  return nonbonded_model.deviceViewToTemplateFreeHostData();
}
#  endif

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::upload() {
  coefficients.upload();
  tick_marks.upload();
  nonbonded_model.upload();
  basis.upload();
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::download() {
  coefficients.download();
  tick_marks.download();
  nonbonded_model.download();
  basis.download();
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::uploadFrame() {
  basis.upload();
  tick_marks.upload();
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::downloadFrame() {
  basis.download();
  tick_marks.download();
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::uploadForceField() {
  nonbonded_model.upload();
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::downloadForceField() {
  nonbonded_model.download();
}
#endif

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                          const double padding, const double spacing,
                                          const int scale_bits_in) {
  setMeshParameters(ag_in, cf_in, padding, std::vector<double>(3, spacing), scale_bits_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const AtomGraph *ag_in, const CoordinateFrame *cf_in,
                                          const double padding,
                                          const std::vector<double> &spacing,
                                          const int scale_bits_in) {
  if (ag_in != nullptr) {
    basis.setTopology(ag_in);
  }
  if (cf_in != nullptr) {
    basis.setCoordinates(cf_in);
  }
  const int actual_scale_bits = (scale_bits_in == -100) ?
                                measurements.getScalingBits() : scale_bits_in;
  measurements = getMeasurements(basis.getTopologyPointer(), basis.getCoordinatePointer(), padding,
                                 spacing, actual_scale_bits);
  validateScalingBits();
  allocate();
  tick_marks = MeshRulers(measurements);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const double padding, const double spacing,
                                          const int scale_bits_in) {
  setMeshParameters(basis.getTopologyPointer(), basis.getCoordinatePointer(), padding,
                    std::vector<double>(3, spacing), scale_bits_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const double padding, const std::vector<double> &spacing,
                                          const int scale_bits_in) {
  setMeshParameters(basis.getTopologyPointer(), basis.getCoordinatePointer(), padding, spacing,
                    scale_bits_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const std::vector<double> &mesh_bounds,
                                          const double spacing, const int scale_bits_in) {
  setMeshParameters(mesh_bounds, std::vector<double>(3, spacing), scale_bits_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setMeshParameters(const std::vector<double> &mesh_bounds,
                                          const std::vector<double> &spacing,
                                          const int scale_bits_in) {
  const int actual_scale_bits = (scale_bits_in == -100) ?
                                measurements.getScalingBits() : scale_bits_in;
  measurements = getMeasurements(mesh_bounds, spacing, actual_scale_bits);
  validateScalingBits();
  allocate();
  tick_marks = MeshRulers(measurements);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setProbeRadius(const double probe_radius_in) {
  if (probe_radius_in < 0.0) {
    rtErr("A probe radius of " + realToString(probe_radius, 7, 4, NumberFormat::STANDARD_REAL) +
          " Angstroms is invalid.", "BackgroundMesh", "setProbeRadius");
  }
  probe_radius = probe_radius_in;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setWellDepth(const double well_depth_in) {
  if (well_depth_in < 0.0) {
    rtErr("A negative well depth of " +
          realToString(probe_radius, 7, 4, NumberFormat::STANDARD_REAL) + " kcal/mol is invalid.  "
          "Use positive numbers to define the depth of a potential energy minimum.",
          "BackgroundMesh", "setWellDepth");
  }
  well_depth = well_depth_in;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setOcclusionPenalty(const double occlusion_penalty_in) {
  occlusion_penalty = occlusion_penalty_in;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::validateCombiningRule(const VdwCombiningRule mixing_protocol_in,
                                              const std::vector<double> &probe_sigma,
                                              const std::vector<double> &probe_epsilon) {

  // Check that the combining rule fits with any provided Lennard-Jones parameter arrays specific
  // to the probe used to draw the mesh.
  switch (mixing_protocol_in) {
  case VdwCombiningRule::GEOMETRIC:
  case VdwCombiningRule::LORENTZ_BERTHELOT:
    if (probe_sigma.size() > 0 || probe_epsilon.size() > 0) {
      rtErr("Pair-specific Lennard-Jones parameters between the probe and atoms of the topology "
            "are only valid if the combining rule is set to " +
            getEnumerationName(VdwCombiningRule::NBFIX) + ".", "BackgroundMesh",
            "validateCombiningRule");
    }
    break;
  case VdwCombiningRule::NBFIX:
    if (probe_sigma.size() != basis.getTopologyPointer()->getLJTypeCount() ||
        probe_sigma.size() != probe_epsilon.size()) {
      rtErr("For " + getEnumerationName(VdwCombiningRule::NBFIX) + " combining rules, all "
            "pairwise parameters between atoms in the topology and the probe must be provided.",
            "BackgroundMesh", "validateCombiningRule");
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setCoefficientScalingBits(const int scaling_bits_in) {

  // Perform checks on the validity of the requested fixed-precision representation.
  if (isFloatingPointScalarType<T>() && scaling_bits_in != 0) {
    rtErr("Mesh coefficients expressed in a floating-point representation cannot have a nonzero "
          "fixed-precision bit count on its coefficients (" + std::to_string(scaling_bits_in) +
          " requested).", "BackgroundMesh", "setCoefficientScalingBits");
  }
  if (scaling_bits_in < 0 || scaling_bits_in > 40) {
    rtErr("A fixed-precision representation of a mesh-based potential cannot have more than 40 "
          "bits after the decimal in its coefficients (" + std::to_string(scaling_bits_in) +
          " requested).", "BackgroundMesh", "setCoefficientScalingBits");
  }
  coefficient_scale_bits = scaling_bits_in;
  coefficient_scale = pow(2.0, scaling_bits_in);
  inverse_coefficient_scale = 1.0 / coefficient_scale;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::computeField(const MeshKlManager &launcher, const PrecisionModel prec,
                                     const HybridTargetLevel availability,
                                     const std::vector<double> &probe_sigma,
                                     const std::vector<double> &probe_epsilon) {
  
  // Loop over all atoms and apply the potential
  switch (kind) {
  case GridDetail::OCCLUSION:
    colorOcclusionMesh(launcher, prec, availability, probe_sigma);
    break;
  case GridDetail::OCCLUSION_FIELD:
    mapOccupancyField(launcher, prec, availability, probe_sigma);
    break;
  case GridDetail::NONBONDED_FIELD:
    mapPureNonbondedPotential(launcher, prec, availability, probe_sigma, probe_epsilon);
    break;
  case GridDetail::NONBONDED_ATOMIC:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::allocate() {
  MeshParamKit mps = measurements.data();
  coefficients.resize(64LLU * static_cast<size_t>(mps.na * mps.nb * mps.nc));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::setNonbondedModel(const VdwCombiningRule mixing_protocol_in,
                                          const double clash_ratio_in,
                                          const double clash_distance_in,
                                          const std::vector<double> &probe_sigma,
                                          const std::vector<double> &probe_epsilon) {

  // Check that the combining rule fits with any provided Lennard-Jones parameter arrays specific
  // to the probe used to draw the mesh.
  const AtomGraph *ag_ptr = basis.getTopologyPointer();
  switch (mixing_protocol_in) {
  case VdwCombiningRule::GEOMETRIC:
  case VdwCombiningRule::LORENTZ_BERTHELOT:
    if (probe_sigma.size() > 0 || probe_epsilon.size() > 0) {
      rtErr("Pair-specific Lennard-Jones parameters between the probe and atoms of the topology "
            "are only valid if the combining rule is set to " +
            getEnumerationName(VdwCombiningRule::NBFIX) + ".", "BackgroundMesh",
            "validateCombiningRule");
    }
    break;
  case VdwCombiningRule::NBFIX:
    if (probe_sigma.size() != ag_ptr->getLJTypeCount() ||
        probe_sigma.size() != probe_epsilon.size()) {
      rtErr("For " + getEnumerationName(VdwCombiningRule::NBFIX) + " combining rules, all "
            "pairwise parameters between atoms in the topology and the probe must be provided.",
            "BackgroundMesh", "validateCombiningRule");
    }
    break;
  }

  // Allocate the non-bonded model and set its overall descriptors
  nonbonded_model = MeshForceField<T>(mixing_protocol_in, clash_ratio_in, clash_distance_in,
                                      ag_ptr);

  // Return immediately if there are no allocations to perform.  Run checks on the data arrays with
  // respect to any van-der Waals potential.  The non-bonded parameters of the mesh come together
  // one piece at a time, owing to the composition and the point at which it is convenient to
  // handle each aspect of the construction.  As of the time this function is called, the
  // non-bonded parameter arrays should have been allocated and set in terms of the Coulomb
  // constant, but not yet populated with coefficients
  // to serve any 
  switch (field) {
  case NonbondedPotential::CLASH:
    if (probe_sigma.size() > 0) {
      nonbonded_model.setLJCoefficients(ag_ptr, probe_sigma, probe_epsilon);
    }
    else {
      nonbonded_model.setLJCoefficients(ag_ptr, probe_radius, well_depth);
    }
    break;
  case NonbondedPotential::VAN_DER_WAALS:

    // The check above will ensure that probe_sigma and probe_epsilon are substantial if and only
    // if the combining rule is NBFIX.
    if (probe_sigma.size() > 0) {
      nonbonded_model.setLJCoefficients(ag_ptr, probe_sigma, probe_epsilon);
    }
    else {
      nonbonded_model.setLJCoefficients(ag_ptr, probe_radius, well_depth);
    }
    nonbonded_model.setLJSoftcoreParameters(measurements.getStencilKind());
    break;
  case NonbondedPotential::ELECTROSTATIC:
    nonbonded_model.setElecSoftcoreParameters(measurements.getStencilKind());
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::validateMeshKind() const {
  switch (kind) {
  case GridDetail::OCCLUSION:    
    if (std::type_index(typeid(T)).hash_code() != ullint_type_index) {
      if (isScalarType<T>()) {
        rtErr("An occlusion mask requires " + getStormmScalarTypeName<ullint>() + " data type, "
              "not " + getStormmScalarTypeName<T>() + ".", "BackgroundMesh", "validateMeshKind");
      }
      else {
        rtErr("An occlusion mask requires " + getStormmScalarTypeName<ullint>() + " data type.",
              "BackgroundMesh", "validateMeshKind");
      }
    }
    break;
  case GridDetail::OCCLUSION_FIELD:
  case GridDetail::NONBONDED_FIELD:
  case GridDetail::NONBONDED_ATOMIC:
    if (isFloatingPointScalarType<T>() == false) {
      rtErr("An non-bonded potential field requires " + getStormmScalarTypeName<float>() + " or " +
            getStormmScalarTypeName<double>() + " data type, not " + getStormmScalarTypeName<T>() +
            ".", "BackgroundMesh", "validateMeshKind");
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void BackgroundMesh<T>::validateScalingBits() const {
  const size_t ct = std::type_index(typeid(T)).hash_code();  
  if (ct == ullint_type_index || ct == float_type_index) {
    const std::vector<double> mesh_origin = measurements.getMeshOrigin<double>();
    std::vector<double> side_a = measurements.getMeshElementVector<double>(UnitCellAxis::A);
    std::vector<double> side_b = measurements.getMeshElementVector<double>(UnitCellAxis::B);
    std::vector<double> side_c = measurements.getMeshElementVector<double>(UnitCellAxis::C);
    elementwiseMultiply<double>(&side_a, measurements.getAxisElementCount(UnitCellAxis::A));
    elementwiseMultiply<double>(&side_b, measurements.getAxisElementCount(UnitCellAxis::B));
    elementwiseMultiply<double>(&side_c, measurements.getAxisElementCount(UnitCellAxis::C));
    double max_log_bound = 0.0;
    for (int i = 0; i < 2; i++) {
      const double dix = static_cast<double>(i) * side_a[0];
      const double diy = static_cast<double>(i) * side_a[1];
      const double diz = static_cast<double>(i) * side_a[2];
      for (int j = 0; j < 2; j++) {
        const double djx = static_cast<double>(j) * side_b[0];
        const double djy = static_cast<double>(j) * side_b[1];
        const double djz = static_cast<double>(j) * side_b[2];
        for (int k = 0; k < 2; k++) {
          const double dkx = static_cast<double>(k) * side_c[0];
          const double dky = static_cast<double>(k) * side_c[1];
          const double dkz = static_cast<double>(k) * side_c[2];
          const std::vector<double> mesh_corner = { mesh_origin[0] + dix + djx + dkx,
                                                    mesh_origin[1] + diy + djy + dky,
                                                    mesh_origin[2] + diz + djz + dkz };
          for (int m = 0; m < 3; m++) {
            if (fabs(mesh_corner[m]) > 0.0) {
              max_log_bound = std::max(max_log_bound, log2(fabs(mesh_corner[m])));
            }
          }
        }
      }
    }
    if (max_log_bound + measurements.getScalingBits() >= 63) {
      rtErr("Occlusion meshes, and non-bonded field meshes using single-precision coefficients, "
            "must keep all elements in the range +/- 32768.0.", "BackgroundMesh",
            "validateScalingBits");
    }
    if (measurements.getScalingBits() > mesh_nonoverflow_bits) {
      rtErr("Occlusion meshes, and non-bonded field meshes using single-precision coefficients, "
            "are assumed to not need positional representations in the extended fixed-precision "
            "model (the overflow bits).  A maximum of " + std::to_string(mesh_nonoverflow_bits) +
            " bits can be used in the represenation (this will rival or exceed the precision of "
            " 64-bit floating point representations for most of space).  A setting of " +
            std::to_string(measurements.getScalingBits()) + " is not acceptable.",
            "BackgroundMesh", "validateScalingBits");
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::colorOcclusionMesh(const MeshKlManager &launcher,
                                           const PrecisionModel prec,
                                           const HybridTargetLevel availability,
                                           const std::vector<double> &probe_sigma) {

  // Check that the coordinates and topology have consistent numbers of atoms.  The coordinates
  // must be represented by a single-structure CoordinateFrame object.
  switch (basis.getReadyFrameCount()) {
  case 0:
    rtErr("An occlusion mesh requires a topology and a structure to have been specified.",
          "BackgroundMesh", "colorOcclusionMesh");
  case 1:
    break;
  default:
    rtErr("An occlusion mesh cannot be formed based on a series of frames.  The OCCLUSION_FIELD "
          "setting is appropriate with such inputs.", "BackgroundMesh", "colorOcclusionMesh");
  }
  const AtomGraph *ag_ptr = basis.getTopologyPointer();
  const CoordinateFrame *cf_ptr = basis.getCoordinatePointer();
  
  // Use the HPC kernel to color the mesh if a GPU is available
#ifdef STORMM_USE_HPC
  if (launcher.getGpu() != null_gpu) {

    // Upload the rulers and frozen atom content
    tick_marks.upload();
    basis.upload();

    // Launch the kernel
    BackgroundMeshWriter<void> bgmw = templateFreeData(HybridTargetLevel::DEVICE);
    switch (availability) {
    case HybridTargetLevel::HOST:
      {
        const CoordinateFrameReader cfr = cf_ptr->deviceViewToHostData();
        launchColorOcclusionMesh(&bgmw, ag_ptr, cfr, prec, launcher);
      }
      break;
    case HybridTargetLevel::DEVICE:
      {
        const CoordinateFrameReader cfr = cf_ptr->data(HybridTargetLevel::DEVICE);
        launchColorOcclusionMesh(&bgmw, ag_ptr, cfr, prec, launcher);
      }
      break;
    }

    // Download the coefficients
    coefficients.download();
    return;
  }
#endif
  
  // Color the mesh on the CPU
  const NonbondedKit<double> nbk = ag_ptr->getDoublePrecisionNonbondedKit();
  const std::vector<bool> mobile_atom = ag_ptr->getAtomMobility();
  const MeshParamKit mps = measurements.data();
  const CoordinateFrameReader cfr = cf_ptr->data();
  
  // Initialize the mesh in CPU memory.
  const int n_elem = mps.na * mps.nb * mps.nc;
  T* coeff_ptr = coefficients.data();
  for (int pos = 0; pos < n_elem; pos++) {
    coeff_ptr[pos] = 0LLU;
  }

  // Replicate the mesh's grid vectors.  Subtract the origin coordinates from the "b" and "c"
  // vectors so that a[i] + b[j] + c[k] = the Cartesian coordinates of any grid point (i,j,k).
  const MeshRulerKit rlrs = tick_marks.data();
  std::vector<double> avx(mps.na + 1), avy(mps.na + 1), avz(mps.na + 1);
  std::vector<double> bvx(mps.nb + 1), bvy(mps.nb + 1), bvz(mps.nb + 1);
  std::vector<double> cvx(mps.nc + 1), cvy(mps.nc + 1), cvz(mps.nc + 1);
  hostInt95ToDouble(avx.data(), avy.data(), avz.data(), rlrs.avec_abs_x, rlrs.avec_abs_x_ovrf,
                    rlrs.avec_abs_y, rlrs.avec_abs_y_ovrf, rlrs.avec_abs_z, rlrs.avec_abs_z_ovrf,
                    mps.na + 1, mps.inv_scale);
  hostInt95ToDouble(bvx.data(), bvy.data(), bvz.data(), rlrs.bvec_x, rlrs.bvec_x_ovrf, rlrs.bvec_y,
                    rlrs.bvec_y_ovrf, rlrs.bvec_z, rlrs.bvec_z_ovrf, mps.nb + 1, mps.inv_scale);
  hostInt95ToDouble(cvx.data(), cvy.data(), cvz.data(), rlrs.cvec_x, rlrs.cvec_x_ovrf, rlrs.cvec_y,
                    rlrs.cvec_y_ovrf, rlrs.cvec_z, rlrs.cvec_z_ovrf, mps.nc + 1, mps.inv_scale);
  const double dorig_x = hostInt95ToDouble(mps.orig_x) * mps.inv_scale;
  const double dorig_y = hostInt95ToDouble(mps.orig_y) * mps.inv_scale;
  const double dorig_z = hostInt95ToDouble(mps.orig_z) * mps.inv_scale;

  // Prepare a buffer to hold subgrid results
  std::vector<ullint> cube_buffer(64, 0LLU);

  // Loop over all atoms in the system
  const int lj_idx_offset = nbk.n_lj_types + 1;
  for (int pos = 0; pos < nbk.natom; pos++) {
    if (mobile_atom[pos]) {
      continue;
    }
    const size_t plj_idx = lj_idx_offset * nbk.lj_idx[pos];
    const double color_radius = (0.5 * nbk.lj_sigma[plj_idx]) + probe_radius;
    const double color_radius_sq = color_radius * color_radius;
    const double atom_x = cfr.xcrd[pos];
    const double atom_y = cfr.ycrd[pos];
    const double atom_z = cfr.zcrd[pos];
    const double atom_dx = atom_x - dorig_x;
    const double atom_dy = atom_y - dorig_y;
    const double atom_dz = atom_z - dorig_z;
    const int icenx = floor((mps.umat[0] * atom_dx) + (mps.umat[3] * atom_dy) +
                            (mps.umat[6] * atom_dz));
    const int iceny = floor((mps.umat[1] * atom_dx) + (mps.umat[4] * atom_dy) +
                            (mps.umat[7] * atom_dz));
    const int icenz = floor((mps.umat[2] * atom_dx) + (mps.umat[5] * atom_dy) +
                            (mps.umat[8] * atom_dz));
    const int pad_a = ceil(color_radius / mps.widths[0]);
    const int pad_b = ceil(color_radius / mps.widths[1]);
    const int pad_c = ceil(color_radius / mps.widths[2]);
    const int ixmin = std::max(icenx - pad_a, 0);
    const int iymin = std::max(iceny - pad_b, 0);
    const int izmin = std::max(icenz - pad_c, 0);
    const int ixmax = std::min(icenx + pad_a + 1, mps.na);
    const int iymax = std::min(iceny + pad_b + 1, mps.nb);
    const int izmax = std::min(icenz + pad_c + 1, mps.nc);
    for (int i = ixmin; i < ixmax; i++) {
      for (int j = iymin; j < iymax; j++) {
        for (int k = izmin; k < izmax; k++) {
          
          // Color the buffer for this atom and this element
          const double base_x = avx[i] + bvx[j] + cvx[k];
          const double base_y = avy[i] + bvy[j] + cvy[k];
          const double base_z = avz[i] + bvz[j] + cvz[k];
          for (size_t m = 0LLU; m < 64LLU; m++) {
            cube_buffer[m] = 0LLU;
          }
          int grid_i = 0;
          for (double di = 0.03125; di < 1.0; di += 0.0625) {
            int grid_j = 0;
            for (double dj = 0.03125; dj < 1.0; dj += 0.0625) {
              int grid_k = 0;
              for (double dk = 0.03125; dk < 1.0; dk += 0.0625) {
                const double xpt = base_x +
                                   (di * mps.invu[0]) + (dj * mps.invu[3]) + (dk * mps.invu[6]);
                const double ypt = base_y +
                                   (di * mps.invu[1]) + (dj * mps.invu[4]) + (dk * mps.invu[7]);
                const double zpt = base_z +
                                   (di * mps.invu[2]) + (dj * mps.invu[5]) + (dk * mps.invu[8]);
                const double dx = xpt - atom_x;
                const double dy = ypt - atom_y;
                const double dz = zpt - atom_z;
                if ((dx * dx) + (dy * dy) + (dz * dz) < color_radius_sq) {
                  const int cubelet_i = grid_i / 4;
                  const int cubelet_j = grid_j / 4;
                  const int cubelet_k = grid_k / 4;
                  const int cubelet_idx = (((cubelet_k * 4) + cubelet_j) * 4) + cubelet_i;
                  const int bit_i = grid_i - (4 * cubelet_i);
                  const int bit_j = grid_j - (4 * cubelet_j);
                  const int bit_k = grid_k - (4 * cubelet_k);
                  const int bit_idx = (((bit_k * 4) + bit_j) * 4) + bit_i;
                  cube_buffer[cubelet_idx] |= (0x1LLU << bit_idx);
                }
                grid_k++;
              }
              grid_j++;
            }
            grid_i++;
          }

          // Accumulate the atom's mapping onto the grid
          const size_t coef_base_idx = 64LLU *
                                       static_cast<size_t>((((k * mps.nb) + j) * mps.na) + i);
          for (size_t m = 0LLU; m < 64LLU; m++) {

            // This is necessary due to the templated nature of the BackgroundMesh object.  In
            // practice, the only way to get here is to have the coefficients array already be of
            // type ullint.
            ullint tcoef = coeff_ptr[coef_base_idx + m];
            tcoef |= cube_buffer[m];
            coeff_ptr[coef_base_idx + m] = tcoef;
          }
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::mapOccupancyField(const MeshKlManager &launcher, const PrecisionModel prec,
                                          const HybridTargetLevel availability,
                                          const std::vector<double> &probe_sigma) {

  // Check that the coordinates and topology have consistent numbers of atoms.  The coordinates
  // must be represented by a single-structure CoordinateFrame object.
  if (basis.getReadyFrameCount() == 0) {
    rtErr("An occlusion field requires a topology and a coordinate series to have been specified.",
          "BackgroundMesh", "mapOccupancyField");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::computeCoefficients(const std::vector<double> &u_grid,
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
                                            const TricubicStencil &tc_weights) {
  
  // Get the measurements abstract locally--better than passing a reference of the same size
  // and then having to de-reference it many times.
  const MeshParamKit mps = measurements.data();

  // Lay out the mesh cell bounds
  std::vector<double> tc_bounds;
  const UnitCellType unit_cell = measurements.getMeshCellType();
  switch (unit_cell) {
  case UnitCellType::NONE:
    break;
  case UnitCellType::ORTHORHOMBIC:
    tc_bounds.resize(6);
    tc_bounds[3] = hostInt95ToDouble(mps.fp_invu[0]) * mps.inv_scale;
    tc_bounds[4] = hostInt95ToDouble(mps.fp_invu[4]) * mps.inv_scale;
    tc_bounds[5] = hostInt95ToDouble(mps.fp_invu[8]) * mps.inv_scale;
    break;
  case UnitCellType::TRICLINIC:
    tc_bounds.resize(12);
    for (int i = 0; i < 9; i++) {
      tc_bounds[i + 3] = hostInt95ToDouble(mps.fp_invu[i]) * mps.inv_scale;
    }
    break;
  }

  // Set constants and pointers to mesh lines for more rapid access
  size_t bounds_offset;
  switch (measurements.getBoundaryConditions()) {
  case BoundaryCondition::ISOLATED:
    bounds_offset = 1;
    break;
  case BoundaryCondition::PERIODIC:
    bounds_offset = 0;
    break;
  }
  const size_t ngrid_a = mps.na + bounds_offset;
  const size_t ngrid_b = mps.nb + bounds_offset;
  const size_t ngrid_c = mps.nc + bounds_offset;
  const MeshRulerKit rlrs = tick_marks.data();
  T* coeff_ptr = coefficients.data();
  
  // Compute coefficeints
  std::vector<double> u_elem(8), dudx_elem(8), dudy_elem(8), dudz_elem(8);
  std::vector<double> dudxx_elem(8), dudxy_elem(8), dudxz_elem(8), dudyy_elem(8), dudyz_elem(8);
  std::vector<double> dudxxx_elem(8), dudxxy_elem(8), dudxxz_elem(8), dudxyy_elem(8);
  std::vector<double> dudxyz_elem(8);
  for (int i = 0; i < mps.na; i++) {
    const int95_t mesh_ax = { rlrs.avec_abs_x[i], rlrs.avec_abs_x_ovrf[i] };
    const int95_t mesh_ay = { rlrs.avec_abs_y[i], rlrs.avec_abs_y_ovrf[i] };
    const int95_t mesh_az = { rlrs.avec_abs_z[i], rlrs.avec_abs_z_ovrf[i] };
    for (int j = 0; j < mps.nb; j++) {
      const int95_t mesh_abx = hostSplitFPSum(mesh_ax, rlrs.bvec_x[j], rlrs.bvec_x_ovrf[j]);
      const int95_t mesh_aby = hostSplitFPSum(mesh_ay, rlrs.bvec_y[j], rlrs.bvec_y_ovrf[j]);
      const int95_t mesh_abz = hostSplitFPSum(mesh_az, rlrs.bvec_z[j], rlrs.bvec_z_ovrf[j]);
      for (int k = 0; k < mps.nc; k++) {
        const int95_t mesh_abcx = hostSplitFPSum(mesh_abx, rlrs.cvec_x[k], rlrs.cvec_x_ovrf[k]);
        const int95_t mesh_abcy = hostSplitFPSum(mesh_aby, rlrs.cvec_y[k], rlrs.cvec_y_ovrf[k]);
        const int95_t mesh_abcz = hostSplitFPSum(mesh_abz, rlrs.cvec_z[k], rlrs.cvec_z_ovrf[k]);
        const size_t nijk_base = 8LLU * static_cast<size_t>((((k * mps.nb) + j) * mps.na) + i);

        // Compose the input vectors
        switch (measurements.getBoundaryConditions()) {
        case BoundaryCondition::ISOLATED:
          for (int ci = 0; ci < 2; ci++) {
            const size_t ici = i + ci;
            for (int cj = 0; cj < 2; cj++) {
              const size_t jcj = j + cj;
              for (int ck = 0; ck < 2; ck++) {
                const size_t nv = (((ck * 2) + cj) * 2) + ci;
                const size_t kck = k + ck;
                const size_t nijk = (((kck * ngrid_b) + jcj) * ngrid_a) + ici;
                u_elem[nv]    = u_grid[nijk];
                dudx_elem[nv] = dudx_grid[nijk];
                dudy_elem[nv] = dudy_grid[nijk];
                dudz_elem[nv] = dudz_grid[nijk];
                switch (measurements.getStencilKind()) {
                case Interpolant::SMOOTHNESS:
                  dudxy_elem[nv]  = dudxy_grid[nijk];
                  dudxz_elem[nv]  = dudxz_grid[nijk];
                  dudyz_elem[nv]  = dudyz_grid[nijk];
                  dudxyz_elem[nv] = dudxyz_grid[nijk];
                  if (unit_cell == UnitCellType::TRICLINIC) {
                    dudxx_elem[nv]  = dudxx_grid[nijk];
                    dudyy_elem[nv]  = dudyy_grid[nijk];
                    dudxxx_elem[nv] = dudxxx_grid[nijk];
                    dudxxy_elem[nv] = dudxxy_grid[nijk];
                    dudxxz_elem[nv] = dudxxz_grid[nijk];
                    dudxyy_elem[nv] = dudxyy_grid[nijk];
                  }
                  break;
                case Interpolant::FUNCTION_VALUE:
                  dudxy_elem[nv]  = dudxy_grid[nijk_base + nv];
                  dudxz_elem[nv]  = dudxz_grid[nijk_base + nv];
                  dudyz_elem[nv]  = dudyz_grid[nijk_base + nv];
                  dudxyz_elem[nv] = dudxyz_grid[nijk_base + nv];
                  break;
                }
              }
            }
          }
          break;
        case BoundaryCondition::PERIODIC:

          // The input grids are expected to be the same dimensions as the mesh, but it order to
          // get the outer boundaries of points at any of the upper mesh limits, the indices
          // referenced in each of the computed grids must wrap.
          for (int ci = 0; ci < 2; ci++) {
            size_t ici = i + ci;
            ici -= (ici == ngrid_a) * ngrid_a;
            for (int cj = 0; cj < 2; cj++) {
              size_t jcj = j + cj;
              jcj -= (jcj == ngrid_b) * ngrid_b;
              for (int ck = 0; ck < 2; ck++) {
                const size_t nv = (((ck * 2) + cj) * 2) + ci;
                size_t kck = k + ck;
                kck -= (kck == ngrid_c) * ngrid_c;
                const size_t nijk = (((kck * ngrid_b) + jcj) * ngrid_a) + ici;
                u_elem[nv] = u_grid[nijk];
                dudx_elem[nv] = dudx_grid[nijk];
                dudy_elem[nv] = dudy_grid[nijk];
                dudz_elem[nv] = dudz_grid[nijk];
                switch (measurements.getStencilKind()) {
                case Interpolant::SMOOTHNESS:
                  dudxy_elem[nv] = dudxy_grid[nijk];
                  dudxz_elem[nv] = dudxz_grid[nijk];
                  dudyz_elem[nv] = dudyz_grid[nijk];
                  dudxyz_elem[nv] = dudxyz_grid[nijk];
                  if (unit_cell == UnitCellType::TRICLINIC) {
                    dudxx_elem[nv] = dudxx_grid[nijk];
                    dudyy_elem[nv] = dudyy_grid[nijk];
                    dudxxx_elem[nv] = dudxxx_grid[nijk];
                    dudxxy_elem[nv] = dudxxy_grid[nijk];
                    dudxxz_elem[nv] = dudxxz_grid[nijk];
                    dudxyy_elem[nv] = dudxyy_grid[nijk];
                  }
                  break;
                case Interpolant::FUNCTION_VALUE:
                  dudxy_elem[nv]  = dudxy_grid[nijk_base + nv];
                  dudxz_elem[nv]  = dudxz_grid[nijk_base + nv];
                  dudyz_elem[nv]  = dudyz_grid[nijk_base + nv];
                  dudxyz_elem[nv] = dudxyz_grid[nijk_base + nv];                  
                  break;
                }
              }
            }
          }
          break;
        }

        // Compose the mesh element.  Complete the bounds array by adding the origin, then
        // compute the tricubic coefficients.
        tc_bounds[0] = hostInt95ToDouble(mesh_abcx) * mps.inv_scale;
        tc_bounds[1] = hostInt95ToDouble(mesh_abcy) * mps.inv_scale;
        tc_bounds[2] = hostInt95ToDouble(mesh_abcz) * mps.inv_scale;
        std::vector<double> tc_coeffs;
        switch (measurements.getStencilKind()) {
        case Interpolant::SMOOTHNESS:
          switch (unit_cell) {
          case UnitCellType::NONE:
            break;
          case UnitCellType::ORTHORHOMBIC:
            {
              const TricubicCell<double> tc_elem(tc_weights, tc_bounds, u_elem, dudx_elem,
                                                 dudy_elem, dudz_elem, dudxy_elem, dudxz_elem,
                                                 dudyz_elem, dudxyz_elem);
              tc_coeffs = tc_elem.getCoefficients();
            }
            break;
          case UnitCellType::TRICLINIC:
            {
              const TricubicCell<double> tc_elem(tc_weights, tc_bounds, u_elem, dudx_elem,
                                                 dudy_elem, dudz_elem, dudxx_elem, dudxy_elem,
                                                 dudxz_elem, dudyy_elem, dudyz_elem, dudxxx_elem,
                                                 dudxxy_elem, dudxxz_elem, dudxyy_elem,
                                                 dudxyz_elem);
              tc_coeffs = tc_elem.getCoefficients();
            }
            break;
          }
          break;
        case Interpolant::FUNCTION_VALUE:
          {
            const TricubicCell<double> tc_elem(tc_weights, tc_bounds, u_elem, dudx_elem, dudy_elem,
                                               dudz_elem, dudxy_elem, dudxz_elem, dudyz_elem,
                                               dudxyz_elem);
            tc_coeffs = tc_elem.getCoefficients();
          }            
          break;
        }
        const size_t tc_offset = static_cast<size_t>((((k * mps.nb) + j) * mps.na) + i) * 64LLU;
        for (size_t cpos = 0; cpos < 64; cpos++) {
          coeff_ptr[tc_offset + cpos] = tc_coeffs[cpos];
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::computeCoefficients(const std::vector<double> &u_grid,
                                            const std::vector<double> &dudx_grid,
                                            const std::vector<double> &dudy_grid,
                                            const std::vector<double> &dudz_grid,
                                            const std::vector<double> &dudxy_grid,
                                            const std::vector<double> &dudxz_grid,
                                            const std::vector<double> &dudyz_grid,
                                            const std::vector<double> &dudxyz_grid,
                                            const TricubicStencil &tc_weights) {
  if (measurements.getMeshCellType() == UnitCellType::TRICLINIC) {
    rtErr("Derivatives for d2/dx2, d2/dy2, d3/dx3, d3/dx2y, d3/dx2z, and d3/dxy2 must be provided "
          "for a non-orthorhombic mesh.", "BackgroundMesh", "computeCoefficients");
  }
  computeCoefficients(u_grid, dudx_grid, dudy_grid, dudz_grid, {}, dudxy_grid, dudxz_grid, {},
                      dudyz_grid, {}, {}, {}, {}, dudxyz_grid, tc_weights);
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tcalc, typename Tcalc2, typename Tcalc3, typename Tcalc4>
void BackgroundMesh<T>::colorNonbondedField(const NonbondedKit<Tcalc> &nbk,
                                            const TricubicStencil &tc_weights,
                                            const std::vector<double> &sigma_table,
                                            const std::vector<double> &eps_table) {
  const std::vector<bool> mobile_atom = basis.getTopologyPointer()->getAtomMobility();
  const MeshParamKit mps = measurements.data();
  const CoordinateFrameReader cfr = basis.getCoordinatePointer()->data();
  const Tcalc value_one   = 1.0;
  const Tcalc value_two   = 2.0;
  const Tcalc value_six   = 6.0;
  const Tcalc value_n12   = -12.0;
  const Tcalc value_n42   = -42.0;
  const Tcalc value_156   = 156.0;
  const Tcalc value_336   = 336.0;
  const Tcalc value_n2184 = -2184.0;
  const bool tcalc_is_double = std::type_index(typeid(Tcalc)).hash_code();
  const T coeff_scl = coefficient_scale;
  
  // Initialize the mesh in CPU memory.
  const int n_elem = mps.na * mps.nb * mps.nc;
  T* coeff_ptr = coefficients.data();
  const T value_zero = 0.0;
  for (int pos = 0; pos < n_elem; pos++) {
    coeff_ptr[pos] = value_zero;
  }

  // Set pointers to mesh lines for more rapid access
  const MeshRulerKit rlrs = tick_marks.data();
  const MeshFFKit<double> mnbk = nonbonded_model.referenceData();
  
  // Create a series of eight three-dimensional grids that will yield the mesh coefficients for
  // each element.  This allocates a small amount of additional memory, relative to the mesh that
  // will be produced, but reduces the pre-computations for getting those elements by a factor
  // approaching eight.  Since this is all being done on the CPU, use std::vector<double> objects
  // and compute the coefficients in double-precision (64-bit), even if they will be used in
  // 32-bit representations.  Here, "grid" will refer to the three-dimensional series of regular
  // points on which the electrostatic potential U, dU/dx, dU/dy, ..., d2U/dxdy, ..., d3U/dxdydz
  // are computed, while "mesh" will refer to the collection of coefficients for each corresponding
  // element.
  const UnitCellType unit_cell = measurements.getMeshCellType();
  UnitCellType imaging_cell;
  size_t bounds_offset;
  switch (measurements.getBoundaryConditions()) {
  case BoundaryCondition::ISOLATED:
    imaging_cell = UnitCellType::NONE;
    bounds_offset = 1;
    break;
  case BoundaryCondition::PERIODIC:
    imaging_cell = unit_cell;
    bounds_offset = 0;
    break;
  }
  const size_t ngrid_a = mps.na + bounds_offset;
  const size_t ngrid_b = mps.nb + bounds_offset;
  const size_t ngrid_c = mps.nc + bounds_offset;
  const size_t ngabc = ngrid_a * ngrid_b * ngrid_c;
  const size_t nmabc = static_cast<size_t>(mps.na) * static_cast<size_t>(mps.nb) *
                       static_cast<size_t>(mps.nc) * 8LLU;
  std::vector<llint> u_grid(ngabc, 0LL);
  std::vector<llint> dudx_grid(ngabc, 0LL), dudy_grid(ngabc, 0LL), dudz_grid(ngabc, 0LL);
  std::vector<llint> dudxy_grid, dudxz_grid, dudyz_grid, dudxyz_grid;
  std::vector<llint> dudxx_grid, dudyy_grid, dudxxx_grid, dudxxy_grid, dudxxz_grid, dudxyy_grid;
  switch (measurements.getStencilKind()) {
  case Interpolant::SMOOTHNESS:
    dudxy_grid.resize(ngabc, 0LL);
    dudxz_grid.resize(ngabc, 0LL);
    dudyz_grid.resize(ngabc, 0LL);
    dudxyz_grid.resize(ngabc, 0LL);
    if (unit_cell == UnitCellType::TRICLINIC) {
      dudxx_grid.resize(ngabc, 0LL);
      dudyy_grid.resize(ngabc, 0LL);
      dudxxx_grid.resize(ngabc, 0LL);
      dudxxy_grid.resize(ngabc, 0LL);
      dudxxz_grid.resize(ngabc, 0LL);
      dudxyy_grid.resize(ngabc, 0LL);
    }
    break;
  case Interpolant::FUNCTION_VALUE:

    // For interpolating with a determinant weighted towards function values or additional
    // first derivatives, the additional second derivatives from the trace of the tensor as well
    // as other third partial derivatives are not needed.  However, some information for fitting
    // the interpolant coefficients will become specific to each mesh element, which raises some
    // of the memory requirements in this scheme by a factor approaching eight.
    dudxy_grid.resize(nmabc, 0LL);
    dudxz_grid.resize(nmabc, 0LL);
    dudyz_grid.resize(nmabc, 0LL);
    dudxyz_grid.resize(nmabc, 0LL);
    break;
  }

  // Prepare for the periodic case by making transformation matrices spanning the entire mesh.
  const double real_na = mps.na;
  const double real_nb = mps.nb;
  const double real_nc = mps.nc;
  const std::vector<double> region_invu = { real_na * mps.invu[0], real_na * mps.invu[1],
                                            real_na * mps.invu[2],
                                            real_nb * mps.invu[3], real_nb * mps.invu[4],
                                            real_nb * mps.invu[5],
                                            real_nc * mps.invu[6], real_nc * mps.invu[7],
                                            real_nc * mps.invu[8] };
  std::vector<double> region_umat(9);
  invertSquareMatrix(region_invu, &region_umat);
  
  // If necessary, pre-compute the coordinates of additional points for the mesh element stencil.
  std::vector<int95_t> xy_box_x(8), xy_box_y(8), xy_box_z(8);
  std::vector<int95_t> xz_box_x(8), xz_box_y(8), xz_box_z(8);
  std::vector<int95_t> yz_box_x(8), yz_box_y(8), yz_box_z(8);
  std::vector<int95_t> cn_box_x(8), cn_box_y(8), cn_box_z(8);
  switch (measurements.getStencilKind()) {
  case Interpolant::SMOOTHNESS:
    break;
  case Interpolant::FUNCTION_VALUE:
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        for (int k = 0; k < 2; k++) {
          const size_t pt_idx = (((k * 2) + j) * 2) + i;
          double dx, dy, dz;
          fvStencilCoordinates<double>(0.0, 0.0, 0.0, UnitCellAxis::C, i, j, k, mps.invu,
                                       &dx, &dy, &dz);
          xy_box_x[pt_idx] = hostDoubleToInt95(dx * mps.scale);
          xy_box_y[pt_idx] = hostDoubleToInt95(dy * mps.scale);
          xy_box_z[pt_idx] = hostDoubleToInt95(dz * mps.scale);
          fvStencilCoordinates<double>(0.0, 0.0, 0.0, UnitCellAxis::B, i, j, k, mps.invu,
                                       &dx, &dy, &dz);
          xz_box_x[pt_idx] = hostDoubleToInt95(dx * mps.scale);
          xz_box_y[pt_idx] = hostDoubleToInt95(dy * mps.scale);
          xz_box_z[pt_idx] = hostDoubleToInt95(dz * mps.scale);
          fvStencilCoordinates<double>(0.0, 0.0, 0.0, UnitCellAxis::A, i, j, k, mps.invu,
                                       &dx, &dy, &dz);
          yz_box_x[pt_idx] = hostDoubleToInt95(dx * mps.scale);
          yz_box_y[pt_idx] = hostDoubleToInt95(dy * mps.scale);
          yz_box_z[pt_idx] = hostDoubleToInt95(dz * mps.scale);
          fvStencilCoordinates<double>(0.0, 0.0, 0.0, i, j, k, mps.invu, &dx, &dy, &dz);
          cn_box_x[pt_idx] = hostDoubleToInt95(dx * mps.scale);
          cn_box_y[pt_idx] = hostDoubleToInt95(dy * mps.scale);
          cn_box_z[pt_idx] = hostDoubleToInt95(dz * mps.scale);
        }
      }
    }
    break;
  }
  
  // Retrieve coefficients for the general softcore electrostatic function.  It uses a slope of
  // zero as the inter-particle distance approaches zero.
  Tcalc4 elec_softcore_abcd;
  Tcalc2 elec_softcore_ef;
  switch (field) {
  case NonbondedPotential::CLASH:
  case NonbondedPotential::VAN_DER_WAALS:
    break;
  case NonbondedPotential::ELECTROSTATIC:
    elec_softcore_abcd = { static_cast<Tcalc>(mnbk.softcore_qq[0]),
                           static_cast<Tcalc>(mnbk.softcore_qq[1]),
                           static_cast<Tcalc>(mnbk.softcore_qq[2]),
                           static_cast<Tcalc>(mnbk.softcore_qq[3]) };
    elec_softcore_ef = { static_cast<Tcalc>(mnbk.softcore_qq[4]),
                         static_cast<Tcalc>(mnbk.softcore_qq[5]) };
    break;
  }
  for (int pos = 0; pos < nbk.natom; pos++) {
    if (mobile_atom[pos]) {
      continue;
    }
    
    // Compute the particle position in the mesh's native fixed-precision format
    const int95_t atom_x = hostDoubleToInt95(cfr.xcrd[pos] * mps.scale);
    const int95_t atom_y = hostDoubleToInt95(cfr.ycrd[pos] * mps.scale);
    const int95_t atom_z = hostDoubleToInt95(cfr.zcrd[pos] * mps.scale);

    // Get the particle's charge.  The mesh coefficients will be computed in kcal/mol-A^n,
    // n = { 0, 1, 2, 3 } depending on the degree of the derivative.
    Tcalc atom_q, rlj_lim, pair_lja, pair_ljb;
    Tcalc4 lj_softcore_abcd;
    Tcalc2 lj_softcore_ef;
    switch (field) {
    case NonbondedPotential::CLASH:
      {
        const size_t atom_lj_idx = nbk.lj_idx[pos];
        rlj_lim  = mnbk.probe_ljsig[atom_lj_idx] * mnbk.clash_ratio;
      }
      break;
    case NonbondedPotential::VAN_DER_WAALS:
      {
        const size_t atom_lj_idx = nbk.lj_idx[pos];
        rlj_lim  = mnbk.probe_ljsig[atom_lj_idx] * mnbk.clash_ratio;
        pair_lja = mnbk.probe_lja[atom_lj_idx];
        pair_ljb = mnbk.probe_ljb[atom_lj_idx];
        lj_softcore_abcd = { static_cast<Tcalc>(mnbk.softcore_lja[atom_lj_idx]),
                             static_cast<Tcalc>(mnbk.softcore_ljb[atom_lj_idx]),
                             static_cast<Tcalc>(mnbk.softcore_ljc[atom_lj_idx]),
                             static_cast<Tcalc>(mnbk.softcore_ljd[atom_lj_idx]) };
        lj_softcore_ef = { static_cast<Tcalc>(mnbk.softcore_lje[atom_lj_idx]),
                           static_cast<Tcalc>(mnbk.softcore_ljf[atom_lj_idx]) };
      }
      break;
    case NonbondedPotential::ELECTROSTATIC:
      atom_q = nbk.charge[pos] * nbk.coulomb_constant;
      break;
    }
    
    // Loop over all grid points
    for (size_t i = 0; i < ngrid_a; i++) {
      const int95_t mesh_ax = { rlrs.avec_abs_x[i], rlrs.avec_abs_x_ovrf[i] };
      const int95_t mesh_ay = { rlrs.avec_abs_y[i], rlrs.avec_abs_y_ovrf[i] };
      const int95_t mesh_az = { rlrs.avec_abs_z[i], rlrs.avec_abs_z_ovrf[i] };
      for (size_t j = 0; j < ngrid_b; j++) {
        const int95_t mesh_abx = hostSplitFPSum(mesh_ax, rlrs.bvec_x[j], rlrs.bvec_x_ovrf[j]);
        const int95_t mesh_aby = hostSplitFPSum(mesh_ay, rlrs.bvec_y[j], rlrs.bvec_y_ovrf[j]);
        const int95_t mesh_abz = hostSplitFPSum(mesh_az, rlrs.bvec_z[j], rlrs.bvec_z_ovrf[j]);
        for (size_t k = 0; k < ngrid_c; k++) {
          const int95_t mesh_abcx = hostSplitFPSum(mesh_abx, rlrs.cvec_x[k], rlrs.cvec_x_ovrf[k]);
          const int95_t mesh_abcy = hostSplitFPSum(mesh_aby, rlrs.cvec_y[k], rlrs.cvec_y_ovrf[k]);
          const int95_t mesh_abcz = hostSplitFPSum(mesh_abz, rlrs.cvec_z[k], rlrs.cvec_z_ovrf[k]);
          
          // Compute the displacements using the mesh's fixed-precision representation, then
          // immediately convert to double for real-valued computations.
          const int95_t fp_disp_x = hostSplitFPSubtract(mesh_abcx, atom_x.x, atom_x.y);
          const int95_t fp_disp_y = hostSplitFPSubtract(mesh_abcy, atom_y.x, atom_y.y);
          const int95_t fp_disp_z = hostSplitFPSubtract(mesh_abcz, atom_z.x, atom_z.y);
          Tcalc disp_x = hostInt95ToDouble(fp_disp_x) * mps.inv_scale;
          Tcalc disp_y = hostInt95ToDouble(fp_disp_y) * mps.inv_scale;
          Tcalc disp_z = hostInt95ToDouble(fp_disp_z) * mps.inv_scale;
          switch (mps.bounds) {
          case BoundaryCondition::ISOLATED:
            break;
          case BoundaryCondition::PERIODIC:
            imageCoordinates<Tcalc, Tcalc>(&disp_x, &disp_y, &disp_z, region_umat.data(),
                                           region_invu.data(), unit_cell,
                                           ImagingMethod::MINIMUM_IMAGE);
            break;
          }
          const Tcalc r2 = (disp_x * disp_x) + (disp_y * disp_y) + (disp_z * disp_z);
          const Tcalc r  = (tcalc_is_double) ? sqrt(r2) : sqrtf(r2);
          const Tcalc invr = value_one / r;
          const Tcalc invr2 = value_one / r2;

          // Compute one of a variety of potentials.  Fill out the chosen stencil in order to fit
          // mesh coefficients for the element.
          Tcalc u, du, d2u, d3u, du_dx, du_dy, du_dz, du_dxy, du_dxz, du_dyz, du_dxx, du_dyy;
          Tcalc du_dxxx, du_dxxy, du_dxxz, du_dxyy, du_dxyz;
          Tcalc blk_xy[8], blk_xz[8], blk_yz[8], blk_cn[8];
          switch (measurements.getStencilKind()) {
          case Interpolant::SMOOTHNESS:
            {
              Tcalc4 u_pkg;
              switch (field) {
              case NonbondedPotential::ELECTROSTATIC:
                if (r < mnbk.clash_distance) {
                  const Tcalc3 tmp_pkg =
                    evaluateQuinticSpline<Tcalc4, Tcalc3, Tcalc2, Tcalc>(elec_softcore_abcd,
                                                                         elec_softcore_ef, r);
                  u_pkg.x = tmp_pkg.x;
                  u_pkg.y = tmp_pkg.y;
                  u_pkg.z = tmp_pkg.z;

                  // While the quintic spline does not have a continuous thrid derivative, the
                  // quintic softcore function is designed to have this feature.  The softcore
                  // polynomial functions are making use of the spline mechanics due to broader
                  // similarities, but are not tabulated splines themselves.
                  if (tcalc_is_double) {
                    u_pkg.w = (((60.0 * elec_softcore_abcd.x * r) +
                                (24.0 * elec_softcore_abcd.y)) * r) + (6.0 * elec_softcore_abcd.z);
                  }
                  else {
                    u_pkg.w = (((60.0f * elec_softcore_abcd.x * r) +
                                (24.0f * elec_softcore_abcd.y)) * r) +
                              (6.0f * elec_softcore_abcd.z);
                  }
                }
                else {
                  u_pkg.x = invr;
                  u_pkg.y = -invr2;
                  u_pkg.z = value_two * invr * invr2;
                  u_pkg.w = -value_six * invr2 * invr2;
                }
                u_pkg.x *= atom_q;
                u_pkg.y *= atom_q;
                u_pkg.z *= atom_q;
                u_pkg.w *= atom_q;
                break;
              case NonbondedPotential::VAN_DER_WAALS:
                if (r < rlj_lim) {
                  const Tcalc3 tmp_pkg =
                    evaluateQuinticSpline<Tcalc4, Tcalc3, Tcalc2, Tcalc>(lj_softcore_abcd,
                                                                         lj_softcore_ef, r);
                  u_pkg.x = tmp_pkg.x;
                  u_pkg.y = tmp_pkg.y;
                  u_pkg.z = tmp_pkg.z;
                  if (tcalc_is_double) {
                    u_pkg.w = (((60.0 * lj_softcore_abcd.x * r) +
                                (24.0 * lj_softcore_abcd.y)) * r) + (6.0 * lj_softcore_abcd.z);
                  }
                  else {
                    u_pkg.w = (((60.0f * lj_softcore_abcd.x * r) +
                                (24.0f * lj_softcore_abcd.y)) * r) + (6.0f * lj_softcore_abcd.z);
                  }
                }
                else {
                  const Tcalc invr3 = invr2 * invr;
                  const Tcalc invr6 = invr3 * invr3;
                  u_pkg.x = ((pair_lja * invr6) - pair_ljb) * invr6;
                  u_pkg.y = ((value_n12   * pair_lja * invr6) + (value_six * pair_ljb)) * invr6 *
                            invr;
                  u_pkg.z = ((value_156   * pair_lja * invr6) + (value_n42 * pair_ljb)) * invr6 *
                            invr2;
                  u_pkg.w = ((value_n2184 * pair_lja * invr6) + (value_336 * pair_ljb)) * invr6 *
                            invr3;
                }
                break;
              case NonbondedPotential::CLASH:
                break;
              }

              // Repackage the function value and first three derivatives, all of which are
              // continuous in the quintic spline.
              u = u_pkg.x;
              du = u_pkg.y;
              d2u = u_pkg.z;
              d3u = u_pkg.w;
            }
            break;
          case Interpolant::FUNCTION_VALUE:
            {
              Tcalc2 u_pkg;
              switch (field) {
              case NonbondedPotential::ELECTROSTATIC:
                u_pkg = evalElecCSC<Tcalc4, Tcalc2, Tcalc>(elec_softcore_abcd, r,
                                                           mnbk.clash_distance, atom_q);
                break;
              case NonbondedPotential::VAN_DER_WAALS:
                u_pkg = evalLennardJonesCSC<Tcalc4, Tcalc2, Tcalc>(lj_softcore_abcd, r, rlj_lim,
                                                                   pair_lja, pair_ljb);
                break;
              case NonbondedPotential::CLASH:
                break;
              }
              
              // Repackage the function value and its first derivative for processing partial
              // derivatives at the mesh vertex (a corner of up to eight mesh elements).
              u = u_pkg.x;
              du = u_pkg.y;
              
              // In addition to the corners, the function value must be evaluated at various other
              // points.  Store the relevant distances, which are all that matter when only
              // potentials are needed, in the buffers that will eventually hold those energies.
              // This work is only needed if the i, j, and k indices are less than the mesh's true
              // na, nb, and nc dimensions, even on a mesh with ISOLATED boundary conditions.
              if (i < mps.na && j < mps.nb && k < mps.nc) {
                for (size_t m = 0; m < 8; m++) {
                  Tcalc4 rp = distance<Tcalc4, Tcalc>(atom_x, atom_y, atom_z,
                                                      hostSplitFPSum(mesh_abcx, xy_box_x[m]),
                                                      hostSplitFPSum(mesh_abcy, xy_box_y[m]),
                                                      hostSplitFPSum(mesh_abcz, xy_box_z[m]),
                                                      region_umat.data(), region_invu.data(),
                                                      imaging_cell, mps.scale);
                  blk_xy[m] = rp.x;
                  rp = distance<Tcalc4, Tcalc>(atom_x, atom_y, atom_z,
                                               hostSplitFPSum(mesh_abcx, xz_box_x[m]),
                                               hostSplitFPSum(mesh_abcy, xz_box_y[m]),
                                               hostSplitFPSum(mesh_abcz, xz_box_z[m]),
                                               region_umat.data(), region_invu.data(),
                                               imaging_cell, mps.scale);
                  blk_xz[m] = rp.x;
                  rp = distance<Tcalc4, Tcalc>(atom_x, atom_y, atom_z,
                                               hostSplitFPSum(mesh_abcx, yz_box_x[m]),
                                               hostSplitFPSum(mesh_abcy, yz_box_y[m]),
                                               hostSplitFPSum(mesh_abcz, yz_box_z[m]),
                                               region_umat.data(), region_invu.data(),
                                               imaging_cell, mps.scale);
                  blk_yz[m] = rp.x;
                  rp = distance<Tcalc4, Tcalc>(atom_x, atom_y, atom_z,
                                               hostSplitFPSum(mesh_abcx, cn_box_x[m]),
                                               hostSplitFPSum(mesh_abcy, cn_box_y[m]),
                                               hostSplitFPSum(mesh_abcz, cn_box_z[m]),
                                               region_umat.data(), region_invu.data(),
                                               imaging_cell, mps.scale);
                  blk_cn[m] = rp.x;
                }
                switch (field) {
                case NonbondedPotential::ELECTROSTATIC:
                  for (size_t m = 0; m < 8; m++) {
                    blk_xy[m] = evalElecCSCND<Tcalc4, Tcalc2, Tcalc>(elec_softcore_abcd, blk_xy[m],
                                                                     mnbk.clash_distance, atom_q);
                    blk_xz[m] = evalElecCSCND<Tcalc4, Tcalc2, Tcalc>(elec_softcore_abcd, blk_xz[m],
                                                                     mnbk.clash_distance, atom_q);
                    blk_yz[m] = evalElecCSCND<Tcalc4, Tcalc2, Tcalc>(elec_softcore_abcd, blk_yz[m],
                                                                     mnbk.clash_distance, atom_q);
                    blk_cn[m] = evalElecCSCND<Tcalc4, Tcalc2, Tcalc>(elec_softcore_abcd, blk_cn[m],
                                                                     mnbk.clash_distance, atom_q);
                  }
                  break;
                case NonbondedPotential::VAN_DER_WAALS:
                  for (size_t m = 0; m < 8; m++) {
                    blk_xy[m] = evalLennardJonesCSCND<Tcalc4, Tcalc2, Tcalc>(lj_softcore_abcd,
                                                                             blk_xy[m], rlj_lim,
                                                                             pair_lja, pair_ljb);
                    blk_xz[m] = evalLennardJonesCSCND<Tcalc4, Tcalc2, Tcalc>(lj_softcore_abcd,
                                                                             blk_xz[m], rlj_lim,
                                                                             pair_lja, pair_ljb);
                    blk_yz[m] = evalLennardJonesCSCND<Tcalc4, Tcalc2, Tcalc>(lj_softcore_abcd,
                                                                             blk_yz[m], rlj_lim,
                                                                             pair_lja, pair_ljb);
                    blk_cn[m] = evalLennardJonesCSCND<Tcalc4, Tcalc2, Tcalc>(lj_softcore_abcd,
                                                                             blk_cn[m], rlj_lim,
                                                                             pair_lja, pair_ljb);
                  }
                  break;
                case NonbondedPotential::CLASH:
                  break;
                }
              }
            }
            break;
          }

          // Accumulate the function derivatives.  First derivatives along each Cartesian axis at
          // the mesh element corners are obligatory and calculated by all stencils.
          du_dx = radialFirstDerivative<Tcalc>(du, disp_x, r);
          du_dy = radialFirstDerivative<Tcalc>(du, disp_y, r);
          du_dz = radialFirstDerivative<Tcalc>(du, disp_z, r);

          // Accumulate additional derivatives or function values at other points inside the mesh
          // element, as necessary to determine a solution to the interpolant coefficients.
          switch (measurements.getStencilKind()) {
          case Interpolant::SMOOTHNESS:
            du_dxy = radialSecondDerivative<Tcalc>(du, d2u, disp_x, disp_y, r, r2);
            du_dxz = radialSecondDerivative<Tcalc>(du, d2u, disp_x, disp_z, r, r2);
            du_dyz = radialSecondDerivative<Tcalc>(du, d2u, disp_y, disp_z, r, r2);
            du_dxyz = radialThirdDerivative<Tcalc>(du, d2u, d3u, disp_x, disp_y, disp_z, r, r2);
            if (unit_cell == UnitCellType::TRICLINIC) {
              du_dxx = radialSecondDerivative<Tcalc>(du, d2u, disp_x, r);
              du_dyy = radialSecondDerivative<Tcalc>(du, d2u, disp_y, r);
              du_dxxx = radialThirdDerivative<Tcalc>(du, d2u, d3u, disp_x, r, r2);
              du_dxxy = radialThirdDerivative<Tcalc>(du, d2u, d3u, disp_x, disp_y, r, r2);
              du_dxxz = radialThirdDerivative<Tcalc>(du, d2u, d3u, disp_x, disp_z, r, r2);
              du_dxyy = radialThirdDerivative<Tcalc>(du, d2u, d3u, disp_y, disp_x, r, r2);
            }
            break;
          case Interpolant::FUNCTION_VALUE:
            break;
          }

          // Log the results.  All interpolants make use of function values and first Cartesian
          // derivatives at the mesh element corners.
          const size_t nijk = (((k * ngrid_b) + j) * ngrid_a) + i;
          if (tcalc_is_double) {
            u_grid[nijk] += llround(u * coeff_scl);
            dudx_grid[nijk] += llround(du_dx * coeff_scl);
            dudy_grid[nijk] += llround(du_dy * coeff_scl);
            dudz_grid[nijk] += llround(du_dz * coeff_scl);
          }
          else {
            u_grid[nijk] += llroundf(u * coeff_scl);
            dudx_grid[nijk] += llroundf(du_dx * coeff_scl);
            dudy_grid[nijk] += llroundf(du_dy * coeff_scl);
            dudz_grid[nijk] += llroundf(du_dz * coeff_scl);
          }

          // Interpolants that only use first derivatives carry the advantage that additional
          // mixed partial derivatives do not need to be calculated.  However, for interpolants
          // fitted using second- and third-order derivatives spanning cells with non-othogonal
          // axes, additional partial derivatives need to be computed and stored.  In contrast,
          // storage requirements for other kernels grow in some arrays 
          switch (measurements.getStencilKind()) {
          case Interpolant::SMOOTHNESS:
            if (tcalc_is_double) {
              dudxy_grid[nijk]  += llround(du_dxy * coeff_scl);
              dudxz_grid[nijk]  += llround(du_dxz * coeff_scl);
              dudyz_grid[nijk]  += llround(du_dyz * coeff_scl);
              dudxyz_grid[nijk] += llround(du_dxyz * coeff_scl);
            }
            else {
              dudxy_grid[nijk]  += llroundf(du_dxy * coeff_scl);
              dudxz_grid[nijk]  += llroundf(du_dxz * coeff_scl);
              dudyz_grid[nijk]  += llroundf(du_dyz * coeff_scl);
              dudxyz_grid[nijk] += llroundf(du_dxyz * coeff_scl);
            }
            if (unit_cell == UnitCellType::TRICLINIC) {
              if (tcalc_is_double) {
                dudxx_grid[nijk]  += llround(du_dxx * coeff_scl);
                dudyy_grid[nijk]  += llround(du_dyy * coeff_scl);
                dudxxx_grid[nijk] += llround(du_dxxx * coeff_scl);
                dudxxy_grid[nijk] += llround(du_dxxy * coeff_scl);
                dudxxz_grid[nijk] += llround(du_dxxz * coeff_scl);
                dudxyy_grid[nijk] += llround(du_dxyy * coeff_scl);
              }
              else {
                dudxx_grid[nijk]  += llroundf(du_dxx * coeff_scl);
                dudyy_grid[nijk]  += llroundf(du_dyy * coeff_scl);
                dudxxx_grid[nijk] += llroundf(du_dxxx * coeff_scl);
                dudxxy_grid[nijk] += llroundf(du_dxxy * coeff_scl);
                dudxxz_grid[nijk] += llroundf(du_dxxz * coeff_scl);
                dudxyy_grid[nijk] += llroundf(du_dxyy * coeff_scl);
              }
            }
            break;
          case Interpolant::FUNCTION_VALUE:
            {
              // Because additional partial derivatives are not needed when building an interpolant
              // based on additional function values, there is no distinction between triclinic and
              // orthorhombic unit cells at this stage.  However, the arrays holding this
              // information, relating to points within each element and off of the mesh lattice,
              // are only valid for i, j, and k less than the mesh's own na, nb, and nc dimensions.
              if (i < mps.na && j < mps.nb && k < mps.nc) {
                const size_t nijk_element = 8LLU *
                                            static_cast<size_t>((((k * mps.nb) + j) * mps.na) + i);
                if (tcalc_is_double) {
                  for (size_t m = 0; m < 8; m++) {
                    dudxy_grid[nijk_element + m]  += llround(blk_xy[m] * coeff_scl);
                    dudxz_grid[nijk_element + m]  += llround(blk_xz[m] * coeff_scl);
                    dudyz_grid[nijk_element + m]  += llround(blk_yz[m] * coeff_scl);
                    dudxyz_grid[nijk_element + m] += llround(blk_cn[m] * coeff_scl);
                  }
                }
                else {
                  for (size_t m = 0; m < 8; m++) {
                    dudxy_grid[nijk_element + m]  += llroundf(blk_xy[m] * coeff_scl);
                    dudxz_grid[nijk_element + m]  += llroundf(blk_xz[m] * coeff_scl);
                    dudyz_grid[nijk_element + m]  += llroundf(blk_yz[m] * coeff_scl);
                    dudxyz_grid[nijk_element + m] += llroundf(blk_cn[m] * coeff_scl);
                  }
                }
              }
            }
            break;
          }
        }
      }
    }
  }
  
  // Mimic the GPU-based protocol by converting the long long integer accumulators to
  // double-precision real numbers prior to multiplying by double-precision transformations to
  // obtain each elements' coefficients, which can then be stored in the mesh's native type.  By
  // converting one array at a time and then resizing the original to zero, memory is conserved.
  std::vector<double> ru_grid = convertData<llint, double>(&u_grid, coeff_scl);
  std::vector<double> rdudx_grid = convertData<llint, double>(&dudx_grid, coeff_scl);
  std::vector<double> rdudy_grid = convertData<llint, double>(&dudy_grid, coeff_scl);
  std::vector<double> rdudz_grid = convertData<llint, double>(&dudz_grid, coeff_scl);
  std::vector<double> rdudxy_grid = convertData<llint, double>(&dudxy_grid, coeff_scl);
  std::vector<double> rdudxz_grid = convertData<llint, double>(&dudxz_grid, coeff_scl);
  std::vector<double> rdudyz_grid = convertData<llint, double>(&dudyz_grid, coeff_scl);
  std::vector<double> rdudxyz_grid = convertData<llint, double>(&dudxyz_grid, coeff_scl);
  std::vector<double> rdudxx_grid, rdudyy_grid;
  std::vector<double> rdudxxx_grid, rdudxxy_grid, rdudxxz_grid, rdudxyy_grid;
  switch (measurements.getStencilKind()) {
  case Interpolant::SMOOTHNESS:
    if (unit_cell == UnitCellType::TRICLINIC) {
      rdudxx_grid = convertData<llint, double>(&dudxx_grid, coeff_scl);
      rdudyy_grid = convertData<llint, double>(&dudyy_grid, coeff_scl);
      rdudxxx_grid = convertData<llint, double>(&dudxxx_grid, coeff_scl);
      rdudxxy_grid = convertData<llint, double>(&dudxxy_grid, coeff_scl);
      rdudxxz_grid = convertData<llint, double>(&dudxxz_grid, coeff_scl);
      rdudxyy_grid = convertData<llint, double>(&dudxyy_grid, coeff_scl);
    }
    break;
  case Interpolant::FUNCTION_VALUE:
    break;
  }
  
  // Convert the computed grid data into a functional mesh
  switch (unit_cell) {
  case UnitCellType::ORTHORHOMBIC:
    computeCoefficients(ru_grid, rdudx_grid, rdudy_grid, rdudz_grid, rdudxy_grid, rdudxz_grid,
                        rdudyz_grid, rdudxyz_grid, tc_weights);
    break;
  case UnitCellType::TRICLINIC:
    computeCoefficients(ru_grid, rdudx_grid, rdudy_grid, rdudz_grid, rdudxx_grid, rdudxy_grid,
                        rdudxz_grid, rdudyy_grid, rdudyz_grid, rdudxxx_grid, rdudxxy_grid,
                        rdudxxz_grid, rdudxyy_grid, rdudxyz_grid, tc_weights);
    break;
  case UnitCellType::NONE:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMesh<T>::mapPureNonbondedPotential(const MeshKlManager &launcher,
                                                  const PrecisionModel prec,
                                                  const HybridTargetLevel availability,
                                                  const std::vector<double> &sigma_table,
                                                  const std::vector<double> &eps_table) {

  // Compute the weights matrix
  TricubicStencil tc_weights(measurements.getStencilKind());
  
  // Use the HPC kernel to color the mesh if a GPU is available
#ifdef STORMM_USE_HPC
  tc_weights.upload();
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  const CoordinateFrameReader cfr = (availability == HybridTargetLevel::HOST) ?
                                    basis.getCoordinatePointer()->deviceViewToHostData() :
                                    basis.getCoordinatePointer()->data(devc_tier);
  if (launcher.getGpu() != null_gpu) {

    // Upload critical components of the mesh object so that they will be available on the GPU.
    // No changes will be made to the coordinate object or the topology.
    uploadFrame();
    uploadForceField();

    // Color the mesh on the GPU.  During evaluations on the mesh, the mesh non-bonded parameters
    // are taken in whatever form the mesh coefficients are found in.  However, for mesh
    // construction, there is a possibility of using double precision arithmetic for coefficients
    // that will ultimately be stored in single precision.  Therefore, always provide
    // double-precision non-bonded mesh parameters, and convert them to the appropriate calculation
    // mode in the kernel.  This mirrors what occurs on the CPU.
    const MeshFFKit<double> mnbk = getReferenceNonbondedKit(devc_tier);
    BackgroundMeshWriter<void> bgmw = templateFreeData(devc_tier);
    launchColorNonbondedFieldMesh(&bgmw, std::type_index(typeid(T)).hash_code(), mnbk, prec,
                                  tc_weights.data(), basis.getTopologyPointer(), cfr, launcher,
                                  availability);

    // Synchronize the computed coefficients on the host and device.
    cudaDeviceSynchronize();
    coefficients.download();
    return;
  }
  else {
#endif
    // Color the mesh on the CPU.  The mesh non-bonded parameters are taken in whatever form the
    // mesh coefficients are found in.  If the calculation mode is SINGLE but the coefficients are
    // double-precision, these coefficients will be converted to single-precision values before
    // being used to compute mesh-based potentials.
    const AtomGraph *ag_ptr = basis.getTopologyPointer();
    switch (prec) {
    case PrecisionModel::DOUBLE:
      {
        const NonbondedKit<double> nbk = ag_ptr->getDoublePrecisionNonbondedKit();
        colorNonbondedField<double, double2, double3, double4>(nbk, tc_weights, sigma_table,
                                                               eps_table);
      }
      break;
    case PrecisionModel::SINGLE:
      {
        const NonbondedKit<float> nbk = ag_ptr->getSinglePrecisionNonbondedKit();
        colorNonbondedField<float, float2, float3, float4>(nbk, tc_weights, sigma_table,
                                                           eps_table);
      }
      break;
    }
#ifdef STORMM_USE_HPC
  }
#endif
}

} // namespace structure
} // namespace stormm
