// -*-c++-*-
#ifndef STORMM_NML_MESH_H
#define STORMM_NML_MESH_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/symbol_values.h"
#include "Parsing/textfile.h"
#include "Potential/energy_enumerators.h"
#include "Structure/structure_enumerators.h"
#include "namelist_element.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using constants::CartesianDimension;
using constants::ExceptionResponse;
using constants::UnitCellAxis;
using energy::NonbondedPotential;
using parse::TextFile;
using parse::WrapTextSearch;
using structure::BoundaryCondition;
using structure::GridDetail;
using structure::MeshPosition;
using structure::RMSDMethod;

/// \brief Default parameters for the &mesh namelist
constexpr int default_mesh_grid_dim = -1;
constexpr int default_mesh_density_averaging_order = 8;
constexpr double default_mesh_grid_angle = 0.5 * symbols::pi;
constexpr double default_mesh_grid_origin = 0.0;
constexpr double default_mesh_grid_spacing = 1.0;
constexpr double default_mesh_elec_damping_range = 0.5;
constexpr double default_mesh_vdw_damping_ratio = 0.8;
constexpr char default_mesh_grid_detail[] = "occlusion";
constexpr char default_mesh_grid_potential[] = "clash";
constexpr char default_mesh_grid_boundary[] = "isolated";

/// \brief Encapsulate the data extracted from a &receptor namelist to define a grid-mapped
///        representation of a rigid macromolecular structure.
///
class MeshControls {
public:
  
  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \param found_nml   Indication of whether the namelist was found in the input file
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &conformer namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  MeshControls(ExceptionResponse policy_in = ExceptionResponse::DIE);

  MeshControls(const TextFile &tf, int *start_line, bool *found_nml,
               ExceptionResponse policy_in = ExceptionResponse::DIE,
               WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  MeshControls(const MeshControls &original) = default;
  MeshControls(MeshControls &&original) = default;
  MeshControls& operator=(const MeshControls &original) = default;
  MeshControls& operator=(MeshControls &&original) = default;
  /// \}

  /// \brief Get the number of grid points in a particular dimension.  This can serve to query a
  ///        particular piece of information without invoking the more complex but comprehensive
  ///        MeshParameters object.
  ///
  /// Overloaded:
  ///   - Technically more correct, specify the unit cell axis
  ///   - Specify the axis by analogy to Cartesian axes (X -> a, Y -> b, Z -> c)
  ///
  /// \param dim  The axis along which to measure the number of mesh points
  /// \{
  int getAxisElementCount(UnitCellAxis dim) const;
  int getAxisElementCount(CartesianDimension dim) const;
  /// \}

  /// \brief Get the mesh spacing along a particular axis.  Overloads and input parameters to this
  ///        function follow from getAxisElementCount() above.
  /// \{
  double getSpacing(UnitCellAxis dim) const;
  double getSpacing(CartesianDimension dim) const;
  /// \}

  /// \brief Get the mesh alpha angle, between the b and c unit cell axes.  The value is returned
  ///        in units of radians.
  double getAlpha() const;

  /// \brief Get the mesh beta angle, between the a and c unit cell axes, in units of radians.
  double getBeta() const;

  /// \brief Get the mesh gamma angle, between the a and b unit cell axes, in units of radians.
  double getGamma() const;
  
  /// \brief Get the origin of the mesh in one Cartesian dimension.
  double getOrigin(CartesianDimension dim) const;

  /// \brief Get buffer widths between mesh boundaries and the receptor structure or structures.
  ///
  /// Overloaded:
  ///   - Technically more correct, specify the unit cell axis
  ///   - Specify the axis by analogy to Cartesian axes (X -> a, Y -> b, Z -> c)
  ///
  /// \param dim  The Cartesian or unit cell axis normal to the face of interest (along which the
  ///             width is relevant)
  /// \{
  double getBufferWidth(UnitCellAxis dim) const; 
  double getBufferWidth(CartesianDimension dim) const;
  /// \}

  /// \brief Get the type of content (potential energy field) that the mesh will represent.  In
  ///        general, programs that use meshes will have specific applications for them, and
  ///        therefore have a concept of what their content should be.  This keyword will often be
  ///        specified to associate its parameters with one of the meshes that a given program
  ///        intends to make, e.g. "the dimensions in this namelist apply to the clash grid,
  ///        whereas the electrostatic non-bonded field uses parameters specified in a separate
  ///        &mesh namelist."
  GridDetail getDetail() const;

  /// \brief Get the type of potential to be expressed on the mesh.  This complements the
  ///        getContent() function and its associated "kind" keyword to specify a particular
  ///        non-bonded potential in situations where programs utilize separate meshes for
  ///        electrostatics and van-der Waals potentials.
  NonbondedPotential getPotential() const;
  
  /// \brief Get the boundary conditions of the mesh
  BoundaryCondition getBoundaries() const;
  
  /// \brief Get the manner in which the mesh is aligned to the rigid molecule it represents.
  MeshPosition getPosition() const;

  /// \brief Get the number of bits after the decimal to be used in composing the mesh-based
  ///        field as well as the positions of its vertices.
  int getScalingBits() const;

  /// \brief Get the switching range for softcore electrostatic potentials.  This will only
  ///        apply to meshes representing an electrostatic potential.
  double getElecClashDistance() const;

  /// \brief Get the switching ratio for softcore van-der Waals potentials.  This will only
  ///        apply to meshes representing a Lennard-Jones (or, perhaps other) van-der Waals
  ///        potential.
  double getVdwClashRatio() const;
  
  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator& getTranscript() const;

  /// \brief Set the number of mesh points along a particular axis.
  ///
  /// Overloaded:
  ///   - Technically more correct, specify the unit cell axis
  ///   - Specify the axis by analogy to Cartesian axes (X -> a, Y -> b, Z -> c)
  ///
  /// \param mesh_points_in  The desired number of mesh points
  /// \param dim             The unit cell axis in question
  /// \{
  void setElementCount(int mesh_points_in, UnitCellAxis dim);
  void setElementCount(int mesh_points_in, CartesianDimension dim);
  /// \}

  /// \brief Set the mesh spacing along a particular axis.  Overloading of this function follows
  ///        from setElementCount() above.
  ///
  /// \param mesh_spacing_in  The desired distance between mesh points, in Angstroms
  /// \param dim              The unit cell axis in question
  /// \{
  void setSpacing(double mesh_spacing_in, UnitCellAxis dim);
  void setSpacing(double mesh_spacing_in, CartesianDimension dim);
  /// \}

  /// \brief Set the mesh alpha angle.
  ///
  /// \param alpha_in  The angle to set, in units of radians.
  void setAlphaAngle(double alpha_in);

  /// \brief Set the mesh beta angle.
  ///
  /// \param beta_in  The angle to set, in units of radians.
  void setBetaAngle(double beta_in);

  /// \brief Set the mesh gamma angle.
  ///
  /// \param gamma_in  The angle to set, in units of radians.
  void setGammaAngle(double gamma_in);

  /// \brief Set the mesh origin along one Cartesian axis.
  ///
  /// \param mesh_origin_in  The desired mesh origin coordinate
  /// \param dim             The Cartesian axis of interest
  void setOrigin(double mesh_origin_in, CartesianDimension dim);

  /// \brief Set the buffer distance between the mesh boundary and the nearest van-der Waals
  ///        sphere.
  ///
  /// Overloaded:
  ///   - Set a single parameter for all three unit cell or Cartesian axes
  ///   - Set separate parameters for secific unit cell or Cartesian axes
  ///
  /// \param buffer_width_in  The buffer width to set, in units of Angstroms
  /// \param dim              The unit cell or Cartesian axis normal to the face of interest
  /// \{
  void setBufferWidth(double buffer_width_in);
  void setBufferWidth(double buffer_width_in, UnitCellAxis dim);
  void setBufferWidth(double buffer_width_in, CartesianDimension dim);
  /// \}

  /// \brief Set the type of mesh that this namelist will be assumed to describe.
  ///
  /// Overloaded:
  ///   - Provide the enumerated value explicitly
  ///   - Translate a string into the appropriate enumeration
  ///
  /// \param kind_in  The type of energetic field
  /// \{
  void setDetail(const std::string &kind_in);
  void setDetail(GridDetail kind_in);
  /// \}
  
  /// \brief Set the potential energy field that the mesh will represent.
  ///
  /// Overloaded:
  ///   - Provide the enumerated value explicitly
  ///   - Translate a string into the appropriate enumeration
  ///
  /// \param potential_in  The type of energetic field
  /// \{
  void setPotential(const std::string &potential_in);
  void setPotential(NonbondedPotential potential_in);
  /// \}

  /// \brief Set the boundary conditions for the mesh.  Overloading in this function follows from
  ///        setPotential() above.
  ///
  /// \param boundaries_in  The chosen boundary conditions
  /// \{
  void setBoundaries(const std::string &boundaries_in);
  void setBoundaries(BoundaryCondition boundaries_in);
  /// \}

  /// \brief Set the number of bits after the decimal to be used in mesh calculations.
  ///
  /// \param scaling_bits_in  The bit count
  void setScalingBits(int mesh_scaling_bits_in);

  /// \brief Set the electrostatic clash distance.
  ///
  /// \param clash_distance_in  The point at which electrostatic interactions switch over to a
  ///                           softcore potential (inverted parabola, with a maximum at r = -1.0)
  void setElecClashDistance(double clash_distance_in);

  /// \brief Set the van-der Waals clash ratio.
  ///
  /// \param clash_ratio_in  The proportion of the van-der Waals (Lennard-Jones) sigma parameter
  ///                        at which interactions for a given pair switch over to a softcore
  ///                        potential (inverted quartic function, with a maximum at r = -1.0)
  void setVdwClashRatio(double clash_ratio_in);
  
private:
  ExceptionResponse policy;  ///< The course to take when encountering bad input
  int mesh_points_a;         ///< Number of mesh points along the mesh (or unit cell) a vector
  int mesh_points_b;         ///< Number of mesh points along the mesh (or unit cell) b vector
  int mesh_points_c;         ///< Number of mesh points along the mesh (or unit cell) c vector
  double mesh_spacing_a;     ///< Regular spacing between mesh points along the unit cell a vector
  double mesh_spacing_b;     ///< Regular spacing between mesh points along the unit cell b vector
  double mesh_spacing_c;     ///< Regular spacing between mesh points along the unit cell c vector
  double mesh_alpha;         ///< Unit cell angle between the b and c unit cell vectors (can also
                             ///<   apply to a paralellepiped mesh in isolated boundary conditions)
  double mesh_beta;          ///< Unit cell angle between the a and c unit cell vectors (can also
                             ///<   apply to a paralellepiped mesh in isolated boundary conditions)
  double mesh_gamma;         ///< Unit cell angle between the a and b unit cell vectors (can also
                             ///<   apply to a paralellepiped mesh in isolated boundary conditions)
  double mesh_origin_x;      ///< Mesh origin along the Cartesian X axis
  double mesh_origin_y;      ///< Mesh origin along the Cartesian Y axis
  double mesh_origin_z;      ///< Mesh origin along the Cartesian Z axis
  double buffer_width_a;     ///< In isolated boundary conditions, this is the minimum distance
                             ///<   between any van-der Waals sphere surface and the nearest (b x c
                             ///<   plane) mesh boundary face.  In periodic boundary conditions,
                             ///<   the same concept applies although the distances technically
                             ///<   refer to the unit cell.
  double buffer_width_b;     ///< Minimum distance between any van-der Waals sphere of the receptor
                             ///<   and the plane of one of the mesh boundary's a x c faces
  double buffer_width_c;     ///< Minimum distance between any van-der Waals sphere of the receptor
                             ///<   and the plane of one of the mesh boundary's a x b faces
  std::string kind;          ///< The type of mesh, e.g. a clash grid, a non-bonded field
  std::string potential;     ///< The potential field that the mesh will display
  std::string boundaries;    ///< Boundary conditions on the mesh
  int mesh_scaling_bits;     ///< Number of bits after the decimal to use when computing the
                             ///<   positions of mesh vertices and the potential expressed on it
  double clash_distance;     ///< The absolute distance at which softcore electrostatic
                             ///<   interactions take over
  double clash_ratio;        ///< The proportion of the van-der Waals sigma parameter for a given
                             ///<   pair of atoms at which softcore van-der Waals interactions
                             ///<   take over

  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;

  /// \brief Test the mesh element's six dimensions (three lengths and three box angles) in order
  ///        to ensure that they produce a valid parallelepiped.
  void validateMeshElement();
};

/// \brief Free function to read the &mesh namelist.
///
/// \param tf          Text of file containing the input deck, read into RAM
/// \param start_line  Line of the input file at which to begin the scan
/// \param found       Indicator that the namelist was found in the input file
/// \param policy      Response to bad inputs
/// \param wrap        Indicate that the search for a &conformer namelist should carry on from the
///                    beginning of an input file if no such namelist is found starting from the
///                    original starting point
NamelistEmulator meshInput(const TextFile &tf, int *start_line, bool *found,
                               ExceptionResponse policy = ExceptionResponse::DIE,
                               WrapTextSearch wrap = WrapTextSearch::NO);

} // namespace namelist
} // namespace stormm

#endif
