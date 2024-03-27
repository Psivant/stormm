// -*-c++-*-
#ifndef STORMM_MESH_MECHANICS_H
#define STORMM_MESH_MECHANICS_H

#include "copyright.h"
#include "Math/vector_ops.h"
#include "Numerics/split_fixed_precision.h"
#include "Parsing/parse.h"
#include "Potential/energy_enumerators.h"
#include "Potential/scorecard.h"
#include "background_mesh.h"
#include "mesh_forcefield.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

using energy::ScoreCard;
using energy::StateVariable;
using numerics::hostChangeFPBits;
using parse::realToString;
using stmath::accumulateBitmask;

/// \brief Interpolate the mesh-based forces and potential (or, for an occlusion mesh, simply read
///        the bitmask) for a collection of Cartesian coordinates.  This works for a single mesh.
///        The value of the potential is returned in the "w" member of each tuple, while the
///        "x", "y", and "z" members respectively hold the Cartesian X, Y, and Z force components.
///
/// Overloaded:
///   - Provide coordinates and properties as C-style arrays with a trusted length
///   - Provide coordinates and properties as Standard Template Library vectors
///   - Provide coordinates and properties as Hybrid objects
///
/// \param bgmr              Background mesh abstract with pointers to mesh values and dimensions
/// \param mnbk              Non-bonded softcore parameters and other constants specific to the
///                          mesh.  The data type of this object must match that of other
///                          properties and the requested calculation precision.
/// \param prop_a            Charge or Lennard-Jones multiplier parameters of particles.  For
///                          a charge mesh, provide an array of atomic partial charges in the
///                          internal atomic units (Coulomb's constant will have already been
///                          folded into the mesh).  Provide the Lennard-Jones A coefficients of
///                          the molecule being superimposed onto the mesh otherwise, as these
///                          will be used to determine sigma and epsilon and then compare those
///                          values to the mesh's inherent probe radius and well depth.
/// \param prob_b            Lennard-Jones B coefficients for the self interactions of each atom
///                          in the molecule superimposed on the mesh.  Like prop_a in the case of
///                          Lennard-Jones interactions this array is indexed by prop_idx.
/// \param prop_idx          Array of Lennard-Jones indices for each atom of the molecule to be
///                          superimposed on the mesh
/// \param xcrd              Cartesian X coordinates of particles
/// \param ycrd              Cartesian Y coordinates of particles
/// \param zcrd              Cartesian Z coordinates of particles
/// \param xcrd_ovrf         Overflow bits for Cartesian X coordinates of particles
/// \param ycrd_ovrf         Overflow bits for Cartesian Y coordinates of particles
/// \param zcrd_ovrf         Overflow bits for Cartesian Z coordinates of particles
/// \param natom             The maximum number of atoms to interpolate, and the trusted length of
///                          the coordinate and force arrays
/// \param sc                Energy tracking object, which will hold the sum of all function values
///                          found at interpolated points on the mesh.
/// \param crd_scaling_bits  Number of bits after the decimal in fixed-precision coordinate
///                            representations
/// \param 
/// \{
template <typename Tdata, typename Tcalc, typename Tcoord, typename Tobs>
double interpolate(const BackgroundMeshReader<Tdata> &bgmr, const MeshFFKit<Tcalc> &mnbk,
                   const Tcalc* prop_a, const Tcalc* prop_b, const int* prop_idx,
                   const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                   const int* xcrd_ovrf, const int* ycrd_ovrf, const int* zcrd_ovrf,
                   int natom, Tobs* vnrg, Tcoord* xfrc, Tcoord* yfrc, Tcoord* zfrc,
                   int* xfrc_ovrf, int* yfrc_ovrf, int* zfrc_ovrf,
                   OffMeshProtocol policy = OffMeshProtocol::DIE, int system_index = 0,
                   int crd_scaling_bits = 0, int frc_scaling_bits = 0, int nrg_scaling_bits = 0);
  
template <typename Tcalc, typename Tcoord>
void interpolate(const BackgroundMeshReader<Tcalc> &bgmr, const MeshFFKit<Tcalc> &mnbk,
                 const NonbondedKit<Tcalc> &nbk, PhaseSpaceWriter *psw, ScoreCard *sc,
                 OffMeshProtocol policy = OffMeshProtocol::DIE, int system_index = 0);
/// \}

/// \brief Compute the "footprint," the number of clashes, for a molecule in an occlusion map.
///        Descriptions of parameters follow from interpolate() above, the function which is
///        called internally.
template <typename Tcalc, typename Tcoord>
double footprint(const BackgroundMeshReader<ullint> &bgmr, const Tcoord* xcrd, const Tcoord* ycrd,
                 const Tcoord* zcrd, const int* xcrd_ovrf, const int* ycrd_ovrf,
                 const int* zcrd_ovrf, const uint* prop_idx, int natom, Tcoord* vclash,
                 OffMeshProtocol policy = OffMeshProtocol::DIE, int system_index = 0,
                 int crd_scaling_bits = 0);
  
} // namespace structure
} // namespace stormm

#include "mesh_mechanics.tpp"

#endif
