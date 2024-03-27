// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMeshComposition<T>::
BackgroundMeshComposition(const MeshFoundation &basis_in, const MeshParameters &measurements_in,
                          const MeshRulers &tick_marks, const GridDetail purpose_in,
                          const std::vector<double> probe_radii, const VdwCombiningRule lj_rule_in,
                          const double clash_forgiveness,
                          const std::vector<PairLJInteraction> &edits) :
  basis{basis_in},
  measurements{measurements_in},
  tick_marks{measurements_in},
  nonbonded_model{},
  purpose{purpose_in},
  mesh_count{static_cast<int>(probe_radii.size())},
  electrostatic_mesh_index{-1},
  mesh_contents{NonbondedPotential::CLASH},
  probe_radii{probe_radii_in},
  well_depths{std::vector<double>(probe_radii_in.size(), 0.0)},
  coefficients{HybridKind::ARRAY, "mesh_comp_coefficients"}
{
  // In this construction, the meshes will serve a series of probe radii and indicate whether there
  // is a clash.  All atoms of any system interacting with the mesh will be classified by the
  // nearest available radius to their own Lennard-Jones sigma parameters.  An occlusion mesh will
  // simply answer "yes" or "no" as to whether placing an atom at some position is valid.  An
  // occlusion field will produce an energy relating to whether it is probable that placing an atom
  // of a given radius at some position is viable.
  switch (purpose) {
  case GridDetail::OCCLUSION:
    if (std::type_index(typeid(T)).hash_code() != ullint_type_index) {
      rtErr(getEnumerationName(purpose) + " meshes must operate with a 64-bit unsigned integer "
            "data type.", "BackgroundMeshComposition");
    }
    break;
  case GridDetail::OCCLUSION_FIELD:
    if (isFloatingPointScalarType<T>() == false) {
      rtErr(getEnumerationName(purpose) + " meshes must operate with a floating-point data type.",
            "BackgroundMeshComposition");
    }
    break;
  case GridDetail::NONBONDED_FIELD:
  case GridDetail::NONBONDED_ATOMIC:
    rtErr("The proper constructor for nonbonded fields includes the topology or topologies of "
          "systems that will interact with each mesh as an argument, plus an indication of the "
          "Lennard-Jones combining rules.", "BackgroundMeshComposition");
  }
  allocate();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMeshComposition<T>::
BackgroundMeshComposition(const MeshFoundation &basis_in, const MeshParameters &measurements_in,
                          const MeshRulers &tick_marks, const GridDetail purpose_in,
                          const AtomGraph &ag_other, const VdwCombiningRule lj_rule_in,
                          const double clash_ratio_in, const double clash_distance_in,
                          const std::vector<PairLJInteraction> &edits) :
{
  // In this construction, the meshes will serve the electrostatics of a system and any of its
  // Lennard-Jones atom types.  The non-bonded model will be set based on the supplied topology.
  switch (purpose) {
  case GridDetail::OCCLUSION:
  case GridDetail::OCCLUSION_FIELD:
    rtErr("The proper constructor for occlusion fields includes a list of relevant probe radii, "
          "and possibly a specification of the Lennard-Jones combining rules to determine how "
          "those radii interact with the Lennard-Jones sigma parameters of the underlying "
          "topology.", "BackgroundMeshComposition");
  case GridDetail::NONBONDED_FIELD:
  case GridDetail::NONBONDED_ATOMIC:
    if (isFloatingPointScalarType<T>() == false) {
      rtErr(getEnumerationName(purpose) + " meshes must operate with a floating-point data type.",
            "BackgroundMeshComposition");
    }
    break;
  }
  if (basis_in.getTopologyPointer() != nullptr) {
    const ComboGraphLJModel lj_tables(basis_in.getTopologyPointer(), poly_ag_other, lj_rule_in,
                                      edits);
  }
  allocate();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BackgroundMeshComposition<T>:
BackgroundMeshComposition(const MeshFoundation &basis_in, const MeshParameters &measurements_in,
                          const MeshRulers &tick_marks, const GridDetail purpose_in,
                          const AtomGraphSynthesis &poly_ag_other,
                          const VdwCombiningRule lj_rule_in, const double clash_ratio_in,
                          const double clash_distance_in,
                          const std::vector<PairLJInteraction> &edits) :
{
  // In this construction, the meshes again will serve the electrostatics and a series of
  // Lennard-Jones atom types, but spanning a whole series of topologies.  The non-bonded model
  // will be set based on the supplied topology.
  switch (purpose) {
  case GridDetail::OCCLUSION:
  case GridDetail::OCCLUSION_FIELD:
    rtErr("The proper constructor for occlusion fields includes a list of relevant probe radii, "
          "and possibly a specification of the Lennard-Jones combining rules to determine how "
          "those radii interact with the Lennard-Jones sigma parameters of the underlying "
          "topology.", "BackgroundMeshComposition");
  case GridDetail::NONBONDED_FIELD:
  case GridDetail::NONBONDED_ATOMIC:
    if (isFloatingPointScalarType<T>() == false) {
      rtErr(getEnumerationName(purpose) + " meshes must operate with a floating-point data type.",
            "BackgroundMeshComposition");
    }
    break;
  }

  // Count the number of meshes that will be needed.  The electrostatic mesh will get index zero.
  mesh_count = 1;
  electrostatic_mesh_index = 0;
  if (basis_in.getTopologyPointer() != nullptr) {
    const ComboGraphLJModel lj_tables(basis_in.getTopologyPointer(), poly_ag_other, lj_rule_in,
                                      edits);
  }
  allocate();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BackgroundMeshComposition<T>::allocate() {
  const MeshParamKit mpsk = measurements.data();
  coefficients.resize(static_cast<size_t>(mpsk.na * mpsk.nb * mpsk.nc) * 64LLU *
                      static_cast<size_t>(mesh_count));
}

} // namespace structure
} // namespace stormm
