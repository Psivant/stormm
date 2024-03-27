// -*-c++-*-
#ifndef STORMM_STRUCTURE_ENUMERATORS_H
#define STORMM_STRUCTURE_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace structure {
  
/// \brief There are two typical ways of re-imaging particles.
enum class ImagingMethod {
  PRIMARY_UNIT_CELL,  ///< Place coordinates such that, in fractional space, they lie in the first
                      ///<   octant, all components within [0, 1).
  MINIMUM_IMAGE,      ///< Place coordinates (more typically, displacements) such that, in
                      ///<   fractional space, all components fall in the range [ -0.5, 0.5 ).
};

/// \brief List methods for computing (positional) Root Mean Squared Deviation (RMSD) between
///        two sets of coordinates
enum class RMSDMethod {
  ALIGN_MASS,     ///< Align the two structures, centering each by their respective centers of
                  ///<   mass, prior to computing mass-weighted positional RMSD
  ALIGN_GEOM,     ///< Align the two structures, centering each by their respective centers of
                  ///<   geometry, prior to computing positional RMSD with no mass weighting
  NO_ALIGN_MASS,  ///< Do not align the two structures prior to computing positional RMSD with
                  ///<   mass weighting
  NO_ALIGN_GEOM   ///< Do not align the two structures prior to computing positional RMSD without
                  ///<   mass weighting
};

/// \brief List the strategies by which RMSD can be computed.  Depending on the mass of equivalent
///        atom groups, it may be necessary to do a combinatorial search of all symmetry-related
///        atoms until enough of the molecule's mass is placed, then test the minutiae of smaller
///        symmetric groups.  On the other end of the spectrum, there may be no symmetric atom
///        groups, which would allow the RMSD to be computed in a single pass.  If no alignment is
///        necessary, then it is possible to sample any symmetry-related groups in order of their
///        dependence on one another, in any order.
enum class RMSDAlignmentProtocol {
  BUILD_CORE,  ///< Sample the possible arrangements of the largest symmetry-related groups until
               ///<   enough mass has been placed, then proceed to...
  ALIGN_CORE,  ///< Align the molecule's massive core (starting with asymmetric atoms) and
               ///<   afterwards test the various arrangements of each symmetry-related group in
               ///<   the order of their dependencies.
  ALIGN_ALL    ///< With no symmetry-related groups, perform a standard RMSD calculation according
               ///<   to the RMSDMethod indicated.
};

/// \brief There are two orders of RMSD calculation: all to one (reference structure), or all to
///        all (matrix).
enum class RMSDTask {
  REFERENCE,  ///< Compute the RMSD of all structures to a reference structure
  MATRIX      ///< Compute the RMSD of all structures following a particular topology to all others
              ///<   following that same topology.
};

/// \brief Virtual site standalone functiosn and kernels call into two categories.
enum class VirtualSiteActivity {
  PLACEMENT,       ///< Place virtual sites after motion of the underlying frame atoms
  TRANSMIT_FORCES  ///< Transmit forces accumulated on virtual sites to their frame atoms
};

/// \brief Enumerate the various levels of grid detail by which a receptor may be mapped.
enum class GridDetail {
  OCCLUSION,        ///< The grid contains only an occlusion mask, with one bit per grid point
                    ///<   stored in a manner like std::vector<bool>.
  OCCLUSION_FIELD,  ///< The grid stores a real-valued occlusion potential, ranging from zero (no
                    ///<   occlusion) to one (maximum occlusion) and multiplied by some constant to
                    ///<   impose an overall wieght of the occlusive restraining potential.  The
                    ///<   occlusive value at each grid point is a function of the total occupancy
                    ///<   at the grid point averaged over multiple aligned receptor structures,
                    ///<   and the function itself can be tuned to reach unity with a steeper or
                    ///<   steadier slope as the occupany goes to 100%.
  NONBONDED_FIELD,  ///< The grid stores electrostatic, Generalized Born radial derivatives,
                    ///<   Generalized Born pairwise energy, and van-der Waals potentials with
                    ///<   tricibic interpolation.
  NONBONDED_ATOMIC  ///< The grid stores all of the potentials found in the NONBONDED_FIELD case,
                    ///<   but only from atom contributions of rigid components of the receptor
                    ///<   which lie far enough from the bin that tricubic interpolation is a sound
                    ///<   approximation.  Other rigid atoms will be represented in a list of
                    ///<   explicit particles with unique segments for each bin.
};

/// \brief List the various operations that can be associated with a mesh grid.
enum class MappingActivity {
  PARTICLE_TO_MESH,  ///< Construct the mesh based on particles
  MESH_TO_PARTICLE   ///< Interpolate, or read mesh values, to map its information back to
                     ///<   particles
};
  
/// \brief The types of clashes that can occur.
enum class ClashKind {
  VAN_DER_WAALS,  ///< The distance between two particles is less than the minimum van-der Waals
                  ///<   (Lennard-Jones) sigma ratio.
  PURE_DISTANCE   ///< The distance between two particles is less than the minimum absolute
                  ///<   distance allowed before a clash is declared
};

/// \brief The choice of behavior when a point selected for interpolation is off of the mesh that
///        provides the potential.  This can only go into effect when the mesh exists in ISOLATED
///        boundary conditions.
enum class OffMeshProtocol {
  DIE,          ///< This is an error.  Do not proceed further.
  EXTRAPOLATE,  ///< Particles that fall off then mesh should have their potentials and forces
                ///<   extrapolated, rather than interpolated, based on the nearest mesh element.
  ZERO          ///< Particles that fall off the mesh experience zero force and contribute zero
                ///<   energy.
};

/// \brief List the various general levels for system sampling.  The maximum number of trials per
///        system will be capped at a user-specifiable quantity given by the "max_system_trials"
///        keyword, or some value determined on a system-by-system basis determined by one of the
///        enumerations herein.
enum class SamplingIntensity {
  MINIMAL,    ///< Sample each rotamer with a user-specified number of torsion angles across the
              ///<   rotatable bond ("rotation_samples").
  LIGHT,      ///< Sample each pair of nearby rotamers, cis-trans isomerizations, and chiral
              ///<   centers, applying all combinations of torsion angles with the sampling density
              ///<   given by the "rotation_samples" keyword.  For each state, other degrees of
              ///<   freedom will be set randomly.
  HEAVY,      ///< Sample each pair of nearby rotamers, applying all combinations of torsion angles
              ///<   with the sampling density given by the "rotation_samples" keyword up to four
              ///<   times, pending the availability of other degrees of freedom which can be set
              ///<   randomly to generate unique starting conditions.
  EXHAUSTIVE  ///< Sample all possible combinations of rotamers, isomerizable cis-trans bonds, and
              ///<   invertible chiral centers, up to the hard limit set forth by the
              ///<   "max_system_trials" keyword.
};

/// \brief Enumerate the types of boundary conditions for a system or range of values.
enum class BoundaryCondition {
  ISOLATED,  ///< There are no boundaries on the range--the variables are isolated in infinite
             ///<   space with no re-imaging considerations.
  PERIODIC   ///< There is a defined length to the space in which variables reside, and to exceed
             ///<   one limit of the range is to re-enter starting at the opposite limit.
};

/// \brief Enumerate the ways in which a mesh can be positioned relative to a biomolecule.
enum class MeshPosition {
  MOLECULE,  ///< The mesh is positioned such that it encloses a rigid molecule, which it
             ///<   supposedly represents, with equal lengths of space in all directions.
             ///<   Non-orthorhombic meshes will be placed according to the minimum distance of
             ///<   any atom to a mesh face (as computed by Hessian normal form), and if the mesh
             ///<   is not large enough to enclose the molecule completely then the overhang will
             ///<   be equalized over all sides.  No rotation of the molecule or mesh coordinate
             ///<   system will be considered in this placement.
  ORIGIN,    ///< The mesh's origin will be placed at the origin of the coordinate system for the
             ///<   molecule it is tailored to represent.
  ARBITRARY  ///< The mesh's origin is an arbitrary point in space.
};

/// \brief Enumerate different approaches to implementing RATTLE distance constraints.
enum class RattleMethod {
  SEQUENTIAL,   ///< Loop over all bonds in the RATTLE group sequentially, perturbing the central
                ///<   atom with each bond and using that perturbed position as the starting point
                ///<   for the next bond.  This approach depends on the order in which bonds
                ///<   appear.
  CENTER_SUM,   ///< The central atom position is kept fixed as the total movement induced by each
                ///<   bond is accrued.  The distal atoms (hydrogens) move with each constraint
                ///<   computation.
#if 0
  SEQ_UNITARY,  ///< The SEQUENTIAL iteration protocol is applied, but constrained groups with only
                ///<   one bound hydrogen will be solved analytically, always entering the
                ///<   adjustment to set the geometry or particle velocities to their target values
                ///<   within the limits of machine precision.
  CNS_UNITARY   ///< The CENTER_SUM iteration protocol is applied, and constrained groups with only
                ///<   one hydrogen will always enter the adjustment to set the constraint to its
                ///<   target values insofar as the machine precision permits.
#endif
};

/// \brief List the choices for constraint application.
enum class ApplyConstraints {
  YES,  ///< Yes, apply geometric constraints
  NO    ///< No, do not apply geometric constraints
};

/// \brief Return a human-readable string describing each enumerated value.  The enumerator in the
///        argument determines the overload to be used.
///
/// \param input  The enumeration of interest
/// \{
std::string getEnumerationName(ImagingMethod input);
std::string getEnumerationName(RMSDMethod input);
std::string getEnumerationName(RMSDAlignmentProtocol input);
std::string getEnumerationName(RMSDTask input);
std::string getEnumerationName(VirtualSiteActivity input);
std::string getEnumerationName(GridDetail input);
std::string getEnumerationName(MappingActivity input);
std::string getEnumerationName(ClashKind input);
std::string getEnumerationName(OffMeshProtocol input);
std::string getEnumerationName(SamplingIntensity input);
std::string getEnumerationName(BoundaryCondition input);
std::string getEnumerationName(MeshPosition input);
std::string getEnumerationName(RattleMethod input);
std::string getEnumerationName(ApplyConstraints input);
/// \}

/// \brief Interpret string input to recognize an RMSD calculation method (which could include
///        alignment or no alignment).  This enumerator class can also be used to specify a
///        particular method for aligning structures in preparation for other procedures.
///
/// \param input  The human-readable sampling keyword input
RMSDMethod translateRMSDMethod(const std::string &input);
  
/// \brief Interpret string input to recognize a specific type of non-bonded potential for a mesh
///        representation of a receptor or other macromolecule.
///
/// \param input  The human-readable sampling keyword input
GridDetail translateGridDetail(const std::string &input);
  
/// \brief Interpret string input to recognize a specific level of sampling effort.
///
/// \param input  The human-readable sampling keyword input
SamplingIntensity translateSamplingIntensity(const std::string &input);

/// \brief Interpret string input to recognize a specific type of boundary conditions.
///
/// \param input  The human-readable sampling keyword input
BoundaryCondition translateBoundaryCondition(const std::string &input);

/// \brief Interpret string input to recognize a specific type of mesh positioning relative to the
///        molecular structure that the mesh shall represent.
///
/// \param input  The human-readable sampling keyword input
MeshPosition translateMeshPosition(const std::string &input);

/// \brief Interpret string input to recognize a specific approach to implementing RATTLE.
///
/// \param input  The human-readable sampling keyword input
RattleMethod translateRattleMethod(const std::string &input);

/// \brief Interpret a string input to recognize a constraint application directive.
///
/// \param input  The constraint application instruction
ApplyConstraints translateApplyConstraints(const std::string &input);

} // namespace structure
} // namespace stormm

#endif
