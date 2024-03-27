// -*-c++-*-
#ifndef STORMM_BACKGROUND_MESH_COMPOSITION_H
#define STORMM_BACKGROUND_MESH_COMPOSITION_H

#include "copyright.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/combograph_ljmodel.h"
#include "mesh_forcefield.h"
#include "mesh_foundation.h"
#include "mesh_parameters.h"
#include "mesh_rulers.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

using topology::ComboGraphLJModel;
using topology::PairLJInteraction;
  
/// \brief A class to hold the details of multiple BackgroundMesh objects taking the same molecular
///        basis.  It is termed a composition, not a synthesis, as it deals with multiple aspects
///        of a single system (or multiple aspects of many snapshots of the same system), not a
///        series of many systems.
template <typename T>
class BackgroundMeshComposition {
public:

  /// \brief The constructor takes various components (composition in the C++ sense of the word)
  ///        and will immediately allocate.  If the developer wishes to modify any of the
  ///        components post-facto, the meshes will likely need to be re-allocated but this is
  ///        a valid way to construct the composition of meshes and the defaults do not imply a
  ///        large amount of memory to allocate.  All of the constructors will build the nonbonded
  ///        model based on information provided in submitted topologies or other arrays, although
  ///        it remains possible to build the nonbonded_model member variable by direct developer
  ///        intervention after the initial construction.  The calculation each mesh in the
  ///        composition will not happen until the evaluate() member function is called.
  ///
  /// \param purpose_in         Indicate the type of representation to build in a series of meshes
  /// \param probe_radii_in     Radii of occlusion probes to use.  The length of this list will
  ///                           determine the number of meshes produced.
  /// \param ag_other           Topology of the system that will interact with the meshes, given
  ///                           when a nonbonded field is to be created
  /// \param poly_ag_other      The synthesis of topologies that will interact with the meshes,
  ///                           given when a nonbonded field is to be created
  /// \param lj_rule_in         Indicate the manner in which probe radii or Lennard-Jones
  ///                           parameters from supplied topologies will interact with those of
  ///                           atoms in the meshes' underlying topology (as specified in basis_in)
  /// \param clash_forgiveness  Constant ratio scaling all probe : atom interactions when creating
  ///                           occlusion meshes
  /// \param clash_ratio_in     The ratio of sigma parameters at which the Lennard-Jones potential
  ///                           hands off to a softcore polynomial
  /// \param clash_distance_in  The absolute distance at which electrostatic interactions
  ///                           transition to a softcore polynomial
  /// \{
  BackgroundMeshComposition(const MeshFoundation &basis_in, const MeshParameters &measurements_in,
                            const MeshRulers &tick_marks, const GridDetail purpose_in,
                            const std::vector<double> probe_radii, VdwCombiningRule lj_rule_in,
                            double clash_forgiveness, const std::vector<PairLJInteraction> &edits);

  BackgroundMeshComposition(const MeshFoundation &basis_in, const MeshParameters &measurements_in,
                            const MeshRulers &tick_marks, GridDetail purpose_in,
                            const AtomGraph &ag_other, VdwCombiningRule lj_rule_in,
                            double clash_ratio_in, double clash_distance_in,
                            const std::vector<PairLJInteraction> &edits);

  BackgroundMeshComposition(const MeshFoundation &basis_in, const MeshParameters &measurements_in,
                            const MeshRulers &tick_marks, GridDetail purpose_in,
                            const AtomGraphSynthesis &poly_ag_other, VdwCombiningRule lj_rule_in,
                            double clash_ratio_in, double clash_distance_in,
                            const std::vector<PairLJInteraction> &edits);
  /// \}

  /// \brief With no POINTER-kind Hybrid objects or other pointers to repair, the default copy and
  ///        move constructors are valid.  With no const member variables, the default copy and
  ///        move assignment operators are valid.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right hand side of an assignment statement
  /// \{
  BackgroundMeshComposition(const BackgroundMeshComposition &original) = default;
  BackgroundMeshComposition(BackgroundMeshComposition &&original) = default;
  BackgroundMeshComposition& operator=(const BackgroundMeshComposition &original) = default;
  BackgroundMeshComposition& operator=(BackgroundMeshComposition &&original) = default;
  /// \}
  
  /// \brief Get the molecular basis of the meshes.
  ///
  /// Overloaded:
  ///   - Return a const reference for a const object
  ///   - Return a mutable pointer for a non-const object
  /// \{
  const MeshFoundation& getMolecularBasis() const;
  MeshFoundation* getMolecularBasis();
  /// \}

  /// \brief Get the dimensions, boundaries, and other settings for all of the meshes contained
  ///        within this object.
  ///
  /// Overloaded:
  ///   - Return a const reference for a const object
  ///   - Return a mutable pointer for a non-const object
  /// \{
  const MeshParameters& getDimensions() const;
  MeshParameters* getDimensions();
  /// \}
  
  /// \brief Set the dimensions, boundaries, and other details common to all meshes contained
  ///        within this object.  This will automatically refresh the rulers (tick_marks member
  ///        variable).
  ///
  /// Overloaded:
  ///   - Compute mesh dimensions based on the underlying topology and coordinates, given the shape
  ///     of a single mesh element already present in the mesh composition's internally stored
  ///     MeshParameters.  Only the A, B, and C axis element counts will change.
  ///   - Set the mesh dimensions explicitly with a new MeshParameters object.
  ///
  /// \param measurements_in  The measurements to set
  /// \{
  void setMeasurements();
  void setMeasurements(const MeshParameters &measurements_in);
  /// \}

  /// \brief Refresh the rulers for the A, B, and C axes of the underlying meshes.  This must be
  ///        called explicitly if manipulating the measurements via a pointer to the measurements
  ///        member variable.  This will be called implicity if submitting an entire new
  ///        MeshParameters object to reset measurements.  The internal measurements will be
  ///        referenced when setting tick_marks.
  void setRulers();

  /// \brief Set the object's non-bonded model based on a series of topologies, a particular
  ///        combining rule, and a specific interpolation scheme.  If provided, record the specific
  ///        mesh indices referenced by each atom type.  This is also the means to set the type of
  ///        representation: for occlusion meshes, provide the list of specific probe radii (all
  ///        atoms will access the mesh index most appropriate for their Lennard-Jones sigma self
  ///        interactions).  If modeling a nonbonded field, electrostatics will implicitly be
  ///        calculated as part of the model.
  ///
  /// Overloaded:
  ///   - Make provisions for all atom types in an AtomGraph
  ///   - Make provisions for all meshes in an AtomGraphSynthesis
  ///
  /// \param ag
  /// \param poly_ag
  /// \param lj_rule
  /// \param probe_radii
  /// \param well_depths
  void setForceField(const AtomGraphSynthesis &poly_ag);
  
private:

  // These components are applicable to any of the meshes contained within.
  MeshFoundation basis;            ///< Contains a pointer to the topology underlying the mesh, as
                                   ///<   well as a mask of the rigid atoms and a pointer to the
                                   ///<   underlying coordinate object.
  MeshParameters measurements;     ///< The mesh dimensions, including decriptors for boundary
                                   ///<   conditions and the type of interpolation stencil.  These
                                   ///<   are common to all meshes in the composition.
  MeshRulers tick_marks;           ///< Cartesian coordinates of markers along the meshes' A, B,
                                   ///<   and C axes denoting the corners of each mesh element.
                                   ///<   These rulers are common to all meshes.
  MeshForceField nonbonded_model;  ///< Critical non-bonded rules governing the way the underlying
                                   ///<   system is expressed on the mesh, as well as softcore
                                   ///<   potential polynomial coefficients to moderate the effects
                                   ///<   of nearby particles on mesh vertices

  /// The grid detail indicates the purpose of this mesh composition.  Combined with one or more
  /// topologies of systems to be evaluated against each mesh (distinct from the topology of the
  /// underlying system, which gives detail to each mesh), or perhaps an array of probe radii, this
  /// setting determines how many and what sort of meshes to create.
  GridDetail purpose;
  
  /// Most of the following elements have cognates in the singleton BackgroundMesh object, but here
  /// are expressed as arrays to concatenate all meshes in a line.  This counter tracks the total
  /// number of unique meshes expressed for the underlying system.
  int mesh_count;

  /// Index of the electrostatic mesh.  If there is no electrostatic mesh, e.g. in an occlusion or
  /// occlusion field setup, this index is meaningless.
  int electrostatic_mesh_index;
  
  /// The non-bonded potentials expressed on each mesh.  Work units will be designed with mesh
  /// indices for each atom to access listed explicitly.
  std::vector<NonbondedPotential> mesh_contents;

  /// Radii of probes used in each mesh--if the meshes are nonbonded fields, these radii are
  /// Lennard-Jones sigma parameters.  If the meshes indicated occlusion, then these are just radii
  /// which combine with the receptor's Lennard-Jones sigma parameters by the Lorentz-Berthelot
  /// combination rule, and will be incrporated into the nonbonded_model parameter arrays.
  std::vector<double> probe_radii;

  /// Well depth for interactions with each mesh probe, usually combined with Lennard-Jones epsilon
  /// parameters of the atoms in the meshes' underlying topology by geometric or Lorentz-Berthelot
  /// rules.  These values are held for convenience, to indicate the properties of each mesh.  For
  /// an electrostatic mesh the value in this array carries no meaning.
  std::vector<double> well_depths;

  /// Coefficients for all meshes.  Each mesh element contains 64 numbers of the specified data
  /// type, indicating coefficients of a tricubic polynomial or occlusion bitmask, arranged with
  /// the A index varying the fastest, followed by B and then C indices.  See the BackgroundMesh
  /// and TricubicCell ojbects for more information.  Because there are 64 values per mesh element,
  /// each mesh is naturally byte-aligned.  Subsequent meshes in the list are arranged back-to-back
  /// in this array.
  Hybrid<T> coefficients.
};

} // namespace structure
} // namespace stormm

#include "background_mesh_composition.tpp"

#endif
