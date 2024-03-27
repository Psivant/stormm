// -*-c++-*-
#ifndef	STORMM_RESTRAINT_BUILDER_H
#define	STORMM_RESTRAINT_BUILDER_H

#include <vector>
#include "copyright.h"
#include "Chemistry/atommask.h"
#include "bounded_restraint.h"

namespace stormm {
namespace restraints {

using chemistry::AtomMask;
using topology::AtomGraph;
using trajectory::CoordinateFrameReader;

/// \brief The distance at which to begin applying a harmonic restraint to prevent two atoms from
///        approaching closer to each other.  Units of Angstroms.
constexpr double default_hb_prevention_limit = 3.1;
  
/// \brief Unguarded struct to hold elemens of a flat bottom restraint.  Two such objects can hold
///        the initial and final states of a restraint, as well as the steps at which each state
///        takes effect.
struct FlatBottomPlan {

  /// \brief The default constructor sets stiffnesses to zero and the activation step to zero.
  ///        Additional constructors touch on the most frequently used combinations of parameters.
  FlatBottomPlan();
  FlatBottomPlan(double k_in, double r_in);
  FlatBottomPlan(double k_in, double r2_in, double r3_in);
  FlatBottomPlan(double k2_in, double k3_in, double r2_in, double r3_in);
  FlatBottomPlan(double k2_in, double k3_in, double r1_in, double r2_in, double r3_in,
                 double r4_in, int step_in);
  
  int activation_step;  // Step number at which this state of the restraint will take (full) effect
  double k2;            // Left-hand harmonic stiffness
  double k3;            // Right-hand harmonic stiffness
  double r1;            // Leftmost boundary between linear and harmonic behavior
  double r2;            // Leftmost boundary of the flat-bottom potential
  double r3;            // Rightmost boundary of the flat-bottom potential
  double r4;            // Rightmost boundard between harmonic and linear behavior
};
  
/// \brief Build a restraint to hold two atoms at a particular distance from one another.
///
/// Overloaded:
///   - Supply one set of restraint parameters for a static restraint, or one that fully activates
///     at a specific time step
///   - Supply two parameter sets to make the endpoints of a time-based restraint
///  
/// \param ag       System topology (needed to supply the pointer for the BoundedRestraint object)
/// \param atom_i   The first atom of the angle
/// \param atom_j   The second atom of the angle
/// \param fb_init  Initial (or static) parameters for the restraint
/// \param fb_init  Final parameters for the restraint
/// \{
BoundedRestraint applyDistanceRestraint(const AtomGraph *ag, int atom_i, int atom_j,
                                        const FlatBottomPlan fb_init);
BoundedRestraint applyDistanceRestraint(const AtomGraph *ag, int atom_i, int atom_j,
                                        const FlatBottomPlan fb_init,
                                        const FlatBottomPlan fb_final);
/// \}

/// \brief Build a restraint to hold three atoms in a particular angle arrangement.
///
/// Overloaded:
///   - Supply one set of restraint parameters for a static restraint, or one that fully activates
///     at a specific time step
///   - Supply two parameter sets to make the endpoints of a time-based restraint
///
/// \param ag       System topology (needed to supply the pointer for the BoundedRestraint object)
/// \param atom_i   The first atom of the angle
/// \param atom_j   The second atom of the angle
/// \param atom_k   The third atom of the angle
/// \param fb_init  Initial (or static) parameters for the restraint
/// \param fb_init  Final parameters for the restraint
/// \{
BoundedRestraint applyAngleRestraint(const AtomGraph *ag, int atom_i, int atom_j, int atom_k,
                                     const FlatBottomPlan fb_init);
BoundedRestraint applyAngleRestraint(const AtomGraph *ag, int atom_i, int atom_j, int atom_k,
                                     const FlatBottomPlan fb_init, const FlatBottomPlan fb_final);
/// \}

/// \brief Build a restraint to hold a dihedral towards a particular value.  This can rotate the
///        group about a particular bond, or secure a conformation in a position to which it has
///        already been rotated.
///
/// Overloaded:
///   - Supply one set of restraint parameters for a static restraint, or one that fully activates
///     at a specific time step
///   - Supply two parameter sets to make the endpoints of a time-based restraint
///  
/// \param ag       System topology (needed to supply the pointer for the BoundedRestraint object)
/// \param atom_i   The first atom of the dihedral
/// \param atom_j   The second atom of the dihedral (first atom of the bond about which to rotate)
/// \param atom_k   The third atom of the dihedral (second atom of the bond about which to rotate)
/// \param atom_l   The fourth atom of the dihedral
/// \param fb_init  Initial (or static) parameters for the restraint
/// \param fb_init  Final parameters for the restraint
/// \{
BoundedRestraint applyDihedralRestraint(const AtomGraph *ag, int atom_i, int atom_j, int atom_k,
                                        int atom_l, const FlatBottomPlan fb_init);
BoundedRestraint applyDihedralRestraint(const AtomGraph *ag, int atom_i, int atom_j, int atom_k,
                                        int atom_l, const FlatBottomPlan fb_init,
                                        const FlatBottomPlan fb_final);
/// \}

/// \brief Perform basic checks the validity of the topology and coordinate frame needed by the
///        following restraint builders.
///
/// \param ag      Pointer to the topology for the system of interest
/// \param cframe  Coordinates of the system in its current state
/// \param mask    Atom mask (must match the topology by some basic checks))
void restraintTopologyChecks(const AtomGraph *ag, const CoordinateFrameReader &cframe,
                             const AtomMask &mask);
  
/// \brief Apply positional restraints to a topology based on an AtomMask, with general parameters
///        for the permittivity and stiffness.
///
/// Overloaded:
///   - Provide an atom mask to indicate the atoms to restrain
///   - Provide a string that makes an atom mask to indicate the atoms to restrain
///   - Provide a std::vector of integers to indicate the atoms to restrain
///
/// \param ag                    Pointer to the topology for the system of interest
/// \param ref_cf                Reference coordinates of the system, for positional targets
/// \param mask                  Atom mask (must match the topology by some basic checks)
/// \param masked_atoms          List of masked atom indices, based on the original topology
/// \param displacement_penalty  Harmonic stiffness constant for displacing particles away from
///                              the target locations.  This is k3 in the terminology of Amber NMR
///                              restraints.
/// \param displacement_onset    Point at which the restraint penalty begins for particles moving
///                              away from their target coordinates.  This is r3 in the Amber NMR
///                              restraint terminology and the default is 0.0.
/// \param displacement_plateau  Maximum displacment away from the target location beyond which
///                              the penalty goes linear and the restraining force is constant.
///                              This is r4 in the Amber NMR restraint terminology and the default
///                              behavior is to set this 5.0 Angstroms past displacement_onset.
/// \param proximity_penalty     Atoms that are not at least proximity_onset far from the target
///                              location will be subject to this harmonic penalty as their
///                              proximity to the target point increases (default 0.0, no
///                              proximity penalty is applied).  k2 in the Amber NMR nomenclature.
/// \param proximity_onset       The distance beneath which a penalty is applied for atoms being
///                              too close to the target location.  This is r2 in the Amber NMR
///                              nomenclature.
/// \param proximity_plateau     The distance beneath which the proximity penality increases
///                              linearly rather than quadratically and the penalty force
///                              flatlines.  This is r1 in the Amber NMR nomeclature.
/// \{
std::vector<BoundedRestraint>
applyPositionalRestraints(const AtomGraph *ag, const CoordinateFrameReader &ref_cfr,
                          const AtomMask &mask, double displacement_penalty,
                          double displacement_onset = 0.0, double displacement_plateau = 16.0,
                          double proximity_penalty = 0.0, double proximity_onset = 0.0,
                          double proximity_plateau = 0.0);

std::vector<BoundedRestraint>
applyPositionalRestraints(const AtomGraph *ag, const CoordinateFrameReader &ref_cfr,
                          const std::string &mask, double displacement_penalty,
                          double displacement_onset = 0.0, double displacement_plateau = 16.0,
                          double proximity_penalty = 0.0, double proximity_onset = 0.0,
                          double proximity_plateau = 0.0);

std::vector<BoundedRestraint>
applyPositionalRestraints(const AtomGraph *ag, const CoordinateFrameReader &ref_cfr,
                          const std::vector<int> &masked_atoms, double displacement_penalty,
                          double displacement_onset = 0.0, double displacement_plateau = 16.0,
                          double proximity_penalty = 0.0, double proximity_onset = 0.0,
                          double proximity_plateau = 0.0);

/// \}

/// \brief Apply restraints to inter-atomic distances such that the participating atoms are not
///        involved in a topological bond or bond angle exclusion (separated by at least three
///        bonds) but also within some cutoff of one another.  This permits a user to restrain the
///        structural geometry at a middle length scale, while allowing local arrangments to relax
///        and some degree of long-ranged motion.
///
/// Overloaded:
///   - Provide an atom mask to indicate the atoms to restrain
///   - Provide a string that makes an atom mask to indicate the atoms to restrain
///   - Provide a list of atom indices to restrain
///
/// \param ag                      Pointer to the topology for the system of interest
/// \param cfr                     Coordinates of the system in its current state
/// \param mask                    Mask of atoms to be held in place, relative to one another
/// \param masked_atoms            List of masked atom indices, based on the original topology
/// \param penalty                 The stiffness of the harmonic penalty to apply.  This is
///                                identical for distance, angle, and dihedral restraints.
/// \param flat_bottom_half_width  The flat bottom of the restraint well extends this amount in
///                                either direction from the current, observed value of the
///                                distance, angle, or dihedral coordinate, in the appropriate
///                                unit system.
/// \param cutoff                  The maximum displacement of any two atoms in the original
///                                structure, beyond which no restraints shall be created.
/// \{
std::vector<BoundedRestraint>
applyDistanceRestraints(const AtomGraph *ag, const CoordinateFrameReader &ref_cfr,
                        const AtomMask &mask, double penalty, double flat_bottom_half_width,
                        double cutoff);

std::vector<BoundedRestraint>
applyDistanceRestraints(const AtomGraph *ag, const CoordinateFrameReader &ref_cfr,
                        const std::string &mask, double penalty, double flat_bottom_half_width,
                        double cutoff);

std::vector<BoundedRestraint>
applyDistanceRestraints(const AtomGraph *ag, const CoordinateFrameReader &ref_cfr,
                        const std::vector<int> &masked_atoms, double penalty,
                        double flat_bottom_half_width, double cutoff);
/// \}
  
/// \brief Build restraints needed to maintain elements of the conformation not intended to change
///        their shape.  This will be accomplished by distance, angle, and dihedral restraints
///        between heavy atoms not intended to move.  
///
/// Overloaded:
///   - Provide an atom mask to indicate the atoms to restrain
///   - Provide a string that makes an atom mask to indicate the atoms to restrain
///   - Provide a list of atom indices to restrain
///
/// \param ag                      Pointer to the topology for the system of interest
/// \param cfr                     Coordinates of the system in its current state
/// \param mask                    Mask of atoms to be held in place, relative to one another
/// \param masked_atoms            List of masked atom indices, based on the original topology
/// \param penalty                 The stiffness of the harmonic penalty to apply.  This is
///                                identical for distance, angle, and dihedral restraints.
/// \param flat_bottom_half_width  The flat bottom of the restraint well extends this amount in
///                                either direction from the current, observed value of the
///                                distance, angle, or dihedral coordinate, in the appropriate
///                                unit system.
/// \param harmonic_maximum        The maximum value of the restraint in its harmonic range,
///                                defining the point at which the penalty goes linear.
/// \{
std::vector<BoundedRestraint>
applyHoldingRestraints(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                       const AtomMask &mask, double penalty, double flat_bottom_half_width = 0.0,
                       double harmonic_maximum = 16.0);

std::vector<BoundedRestraint>
applyHoldingRestraints(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                       const std::string &mask, double penalty,
                       double flat_bottom_half_width = 0.0, double harmonic_maximum = 16.0);
/// \{

/// \brief Apply a set of restraints to prevent the formation of hydrogen bonds between donors
///        and potential acceptors on any one molecule within a system.  The restraints operate
///        by activating a harmonic restraint penalty if the donor and acceptor heavy atoms come
///        too close, and otherwise apply no penalty for the atoms moving away from one another.
///        Restraints against hydrogen bonding will only be applied to donors and acceptors that
///        are not already coutned among the non-bonded exclusions (1:4 neighbors or closer).
///
/// Overloaded:
///   - Accept the topology and coordinates (this will implicitly construct the chemical features,
///     which carries a risk of being a lengthy process)
///   - Accept the topology and pre-computed chemical features of the system
///
/// \param ag              System topology (needed to get atomic numebrs and bonding patterns)
/// \param cfr             Coordinates of the system in its current state
/// \param chemfe          Chemical patterns detected within the topology (identifies hydrogen
///                        donor and lone pair-bearing acceptor atoms)
/// \param penalty         Restraint stiffness constant applied when donors and acceptors approach
/// \param approach_point  Distance below which the harmonic penalty engages
/// \{
std::vector<BoundedRestraint>
applyHydrogenBondPreventors(const AtomGraph *ag, const ChemicalFeatures &chemfe,
                            const double penalty,
                            const double approach_point = default_hb_prevention_limit);

std::vector<BoundedRestraint>
applyHydrogenBondPreventors(const AtomGraph *ag, const CoordinateFrameReader &cfr,
                            const double penalty,
                            const double approach_point = default_hb_prevention_limit);
/// \}
  
} // namespace restraints
} // namespace stormm

#endif
