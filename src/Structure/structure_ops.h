// -*-c++-*-
#ifndef STORMM_STRUCTURE_OPS_H
#define STORMM_STRUCTURE_OPS_H

#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"

namespace stormm {
namespace structure {

using topology::AtomGraph;
using trajectory::CoordinateFrame;
using trajectory::PhaseSpace;
  
/// \brief Compute the location of a molecule's center of mass.
///
/// Overloaded:
///   - Take a topology as a const reference or a const pointer
///   - Take coordinates as a const reference or const pointer to a PhaseSpace or CoordinateFrame
///     object.
///   - Take direct pointers to molecular member atoms, masses, and coordinates, with explicit
///     array reading limits.
///
/// \param ag            System topology
/// \param ps            System coordinates
/// \param cf            System coordinates
/// \param mol_index     Index of the molecule of interest
/// \param mol_start     Beginning index of the molecule of interest in the mol_contents array
/// \param mol_end       Upper bound of atoms for the molecule of interest in the mol_contents
///                      array
/// \param xcrd          Cartesian X coordinates of the molecule
/// \param ycrd          Cartesian Y coordinates of the molecule
/// \param zcrd          Cartesian Z coordinates of the molecule
/// \param masses        Pointer to the atomic masses
/// \param mol_contents  Pointer to the array of atoms (topological indices) contained in the
///                      molecule of interest.  If a nullptr is provided, the order will assumed to
///                      proceed from mol_start to mol_end.
/// \{
double3 centerOfMass(const AtomGraph &ag, const PhaseSpace &ps, int mol_index);

double3 centerOfMass(const AtomGraph *ag, const PhaseSpace *ps, int mol_index);

double3 centerOfMass(const AtomGraph &ag, const CoordinateFrame &cf, int mol_index);

double3 centerOfMass(const AtomGraph *ag, const CoordinateFrame *cf, int mol_index);

double3 centerOfMass(const double* xcrd, const double* ycrd, const double* zcrd,
                     const double* masses, int mol_start = 0, int mol_end = 0,
                     const int* mol_contents = nullptr);
/// \}

/// \brief Posit that a molecular structure is a rigid body free to rotate in space, and compute
///        the torque about its center of mass.  Assume that the molecule is whole, not broken by
///        some sort of imaging artifact.
///
/// Overloaded:
///   - Take a topology as a const reference or a const pointer
///   - Take direct pointers to molecular member atoms, masses, coordinates, and forces, with
///     explicit array reading limits.
///
/// \param ag            System topology
/// \param ps            System coordinates
/// \param mol_index     Index of the molecule of interest
/// \param mol_start     Beginning index of the molecule of interest in the mol_contents array
/// \param mol_end       Upper bound of atoms for the molecule of interest in the mol_contents
///                      array
/// \param xcrd          Cartesian X coordinates of the molecule
/// \param ycrd          Cartesian Y coordinates of the molecule
/// \param zcrd          Cartesian Z coordinates of the molecule
/// \param xfrc          Cartesian X forces acting on each atom in the molecule
/// \param yfrc          Cartesian Y forces acting on each atom in the molecule
/// \param zfrc          Cartesian Z forces acting on each atom in the molecule
/// \param masses        Pointer to the atomic masses
/// \param mol_contents  Pointer to the array of atoms (topological indices) contained in the
///                      molecule of interest
/// \{
double3 molecularTorque(const AtomGraph &ag, const PhaseSpace &ps, int mol_index = 0);

double3 molecularTorque(const AtomGraph *ag, const PhaseSpace *ps, int mol_index = 0);

double3 molecularTorque(const double* xcrd, const double* ycrd, const double* zcrd,
                        const double* xfrc, const double* yfrc, const double* zfrc,
                        const double* masses, const int* mol_contents, int mol_start, int mol_end);
/// \}


} // namespace structure
} // namespace stormm

#endif
