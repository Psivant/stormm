// -*-c++-*-
#ifndef STORMM_ENERGY_ENUMERATORS_H
#define STORMM_ENERGY_ENUMERATORS_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace energy {

/// \brief Enumerate the choices on whether to evaluate the energy... yes or no.  CPU functions
///        are obligated to evaluate the energy, but GPU functions may or may not due to register
///        pressure that it creates.
enum class EvaluateEnergy {
  NO, YES
};

/// \brief Enumerate the choices on whether to evaluate the force... yes or no.
enum class EvaluateForce {
  NO, YES
};

/// \brief Enumerate the choices on whether to evaluate the virial... yes or no.
enum class EvaluateVirial {
  NO, YES
};

/// \brief Enumerate the types of dihedral potentials that can be evaluated in the context of the
///        abstracted routine.
enum class DihedralStyle {
  COSINE, HARMONIC
};

/// \brief Enumerate all state variables that STORMM will track.  This must always increment by one
///        from one enumeration to the next, start at zero, and always end with ALL_STATES.
enum class StateVariable {
  BOND = 0,               ///< Harmonic bond stretching energy
  ANGLE,                  ///< Harmonic angle bending energy
  PROPER_DIHEDRAL,        ///< Proper dihedral energy
  IMPROPER_DIHEDRAL,      ///< Improper (plane-enforcing) dihedral energy
  UREY_BRADLEY,           ///< CHARMM Urey-Bradley angle stretching energy
  CHARMM_IMPROPER,        ///< CHARMM harmonic improper torsion energy
  CMAP,                   ///< Correction map (typically, coupled dihedral-dihedral) energy
  VDW,                    ///< van-der Waals (typically, Lennard Jones) energy from non-bonded
                          ///<   interactions
  VDW_ONE_FOUR,           ///< van-der Waals (typically, Lennard Jones) energy from 1-4 attenuated
                          ///<   interactions
  ELECTROSTATIC,          ///< Electrostatic energy from non-bonded interactions
  ELEC_ONE_FOUR, ///< Electrostatic energy from 1-4 attenuated interactions
  GENERALIZED_BORN,       ///< Generalized Born (implicit solvent) energy
  RESTRAINT,              ///< Energy due to flat-bottom bimodal harmonic potential restraints
  KINETIC,                ///< Energy due to particle motion
  PRESSURE,               ///< System pressure (only computed if virials are accumulated)
  VIRIAL_11,              ///< Virial tensor (1,1) element
  VIRIAL_12,              ///< Virial tensor (1,2) element
  VIRIAL_22,              ///< Virial tensor (2,2) element
  VIRIAL_13,              ///< Virial tensor (1,3) element
  VIRIAL_23,              ///< Virial tensor (2,3) element
  VIRIAL_33,              ///< Virial tensor (3,3) element
  VOLUME,                 ///< Unit cell volume (only relevant to periodic simulations)
  TEMPERATURE_ALL,        ///< Overall system temperature
  TEMPERATURE_PROTEIN,    ///< Temperature of atoms in the "protein" component of the system (this,
                          ///<   like the other temperature subcategories that follow, is an
                          ///<   arbitrary subset of the atoms defined by a special mask in the
                          ///<   control input)
  TEMPERATURE_LIGAND,     ///< Temperature of atoms in the "ligand" component of the system
  TEMPERATURE_SOLVENT,    ///< Temperature of atoms in the "solvent" component of the system
  DU_DLAMBDA,             ///< Derivative of the mixed potential energy function with respect to
                          ///<   the mixing parameter Lambda (relevant to thermodynamic
                          ///<   integration applications only)
  POTENTIAL_ENERGY,       ///< Sum of all potential energy contributions in the system
  TOTAL_ENERGY,           ///< Sum of all potential and kinetic energy components in the system
  ALL_STATES              ///< This must always be the final entry.  The number of tracked
                          ///<   quantities is equal to the value of this entry (ALL_STATES does
                          ///<   not define its own index in the subsequent tracking arrays).
};

/// \brief Enumerate the specific kinds of non-bonded potentials.  This is typically used in the
///        context of a mesh-based potential for a rigid or semi-rigid molecule.
enum class NonbondedPotential {
  ELECTROSTATIC,  ///< Interactions between charges
  VAN_DER_WAALS,  ///< Dispersion iteractions between particles, likely a Lennard-Jones potential
  CLASH           ///< Zero (no clash) or one (clash) based on the van-der Waals radii and the
                  ///<   width of some probe sphere representing heavy atoms in the ligand
};

/// \brief Enumerate the specific kinds of non-bonded potentials.  This is typically used in the
///        context of a neighbor list, decribing what potential the neighbors are relevant to.
enum class NonbondedTheme {
  ELECTROSTATIC,  ///< Interactions between charges
  VAN_DER_WAALS,  ///< Dispersion iteractions between particles, likely a Lennard-Jones potential
  ALL             ///< All non-bonded interactions are covered, typically both electrostatic and
                  ///<   van der Waals
};

/// \brief Various non-bonded potentials can be decomposed into smoothed components which add to
///        recover the whole but can individually be interpolated from meshes with varying degrees
///        of accuracy.
enum class DecomposablePotential {
  ELECTROSTATIC,   ///< The complete Coulombic potential
  DISPERSION,      ///< Dispersion interactions decaying as the sixth power of the distance between
                   ///<   particles.  These interactions may be modified by a switching function to
                   ///<   apply geometric mixing rules after a short-ranged domain in which
                   ///<   Lorentz-Berthelot or NBFix comining rules dominate
  ELEC_PME_DIRECT  ///< The direct particle-particle interactions found in the electrostatic
                   ///<   (Smooth) Particle-Mesh Ewald method
};
  
/// \brief List the methods for combining Lennard-Jones (or other van-der Waals potential)
///        parameters in pairs of interacting atoms.
enum class VdwCombiningRule {
  LORENTZ_BERTHELOT,  ///< Sigma parameters s_i and s_j are combined by taking (s_i + s_j) / 2.
                      ///<   Epsilon parameters are always combined by taking sqrt(E_i * E_j).
  GEOMETRIC,          ///< Sigma parameters s_i and s_j are combined by taking sqrt(s_i * s_j)
  NBFIX               ///< Pair-specific combinations, with some off-diagonal terms not conforming
                      ///<   to either of the other rules
};

/// \brief Functions and kernels can be configured to dampen the effects of clashes (at a minor
///        expense in computation and registers), or not.
enum class ClashResponse {
  NONE,    ///< Do not attempt to dampen the effects of clashes between particles
  FORGIVE  ///< Forgive clashes by not letting the perceived interparticle distance drop below
           ///<   some multiple of the pairwise Lennard-Jones sigma value, or if that is also zero
           ///<   then some arbitrary minimum value.
};

/// \brief List the different modes in which to extract energies from a ScoreCard object.
enum class EnergySample {
  TIME_SERIES,  ///< Return the entire time series of all values stored for a state variable
  FINAL,        ///< Return the current (most recent, or final) value of a state variable
  TIME_AVERAGE  ///< Return the average of the time series of a state variable
};

/// \brief Different shapes for the scaffold used to build a scaffold for a softcore potential
///        function.  Each of these functions are monotonic.
enum class SplineScaffold {
  LINEAR,           ///< The scaffold is a line of points between the function value at the
                    ///<   switching distance r = rs and r = 0.
  SIGMOIDAL,        ///< The scaffold is a sigmoidal function of the form
                    ///<   1 / (1 + e^(-(8r / rs) + 4)).  This normalizes the sigmoidal behavior
                    ///<   such that, over the interval [0, rs], the S-shaped curve traces the
                    ///<   progress that the classic logistic function 1 / (1 + e^-r) would trace
                    ///<   over the interval [ -4, 4 ].
  QUADRATIC_HILL,   ///< The scaffold is an inverted parabola with the value of the original
                    ///<   function at the switching point rs as well as zero derivative at r = 0.
  QUADRATIC_VALLEY  ///< The scaffold is a parabola with the value of the original function and
                    ///<   zero derivative at the switching distance r = rs and the target value at
                    ///<   r = 0.
};

/// \brief The actions to take with a cell grid, necessary to evaluate finite non-bonded
///        contributions and exclusions in a piecemeal fashion.  These actions are useful for
///        analysis of individual functions regarding the most complex component of a typical MD
///        simulation.
enum class CellGridAction {
  INIT_FORCES,       ///< Initialize forces, setting all accumulators of the current image to zero.
  XFER_FORCES,       ///< Transfer forces from the current image to the associated coordinate
                     ///<   synthesis.
  UPDATE_IMG_COORD,  ///< Update the coordinates, the first step in turning the current image into
                     ///<   the new image.  This will leave the CellGrid in an intermediate state,
                     ///<   without a valid current image until executing processes associated with
                     ///<   UPDATE_IMG_CELLS.
  UPDATE_IMG_CELLS   ///< Complete the new image by re-organizing cells.
};

/// \brief Differentiate between different strategies of mapping (spreading) particle density to
///        the particle-mesh interaction grid.
enum class QMapMethod {
  ACC_SHARED,       ///< Force the use of the __shared__ cache accumulation method in whatever
                    ///<   precision model
  GENERAL_PURPOSE,  ///< Force the use of the general-purpose method in whatever precision model
  AUTOMATIC         ///< Use the fastest kernel available based on internal heuristics and
                    ///<   availability under a given precision model
};

/// \brief Enumerate various strategies, including non-automation, for choosing the particle-mesh
///        interaction grid parameters in &namelist input.
enum class PMIStrategy {
  RECOMMENDED,           ///< Take the recommended settings, which will ensure stability in most
                         ///<   simulations without excessive calculations for excessive accuracy.
                         ///<   These settings are akin to the accuracy of Amber default settings,
                         ///<   perhaps slightly better.
  TIGHT,                 ///< Use settings that will produce a level of accuracy, about one order
                         ///<   of magnitude higher than the recommended settings.
  VERY_TIGHT,            ///< Use settings that will produce an extreme level of accuracy, about
                         ///<   two orders of magnitude higher than the recommended strategy.
  RECOMMENDED_PM_HEAVY,  ///< Take settings that would produce the recommended level of accuracy,
                         ///<   but shifting more work into particle-mesh interactions.
  RECOMMENDED_PP_HEAVY,  ///< Take settings that would produce the recommended level of accuracy,
                         ///<   but shifting more work into particle-particle interactions.
  TIGHT_PM_HEAVY,        ///< Take settings that would produce a tight level of accuracy, but
                         ///<   shifting more work into the particle-mesh interactions.
  TIGHT_PP_HEAVY,        ///< Take settings that would produce a tight level of accuracy, but
                         ///<   shifting more work into the particle-mesh interactions.
  NO_AUTOMATION          ///< Do not use any automated settings.
};

/// \brief Enumerate the available sizes of the valence work unit kernel.
enum class ValenceKernelSize {
  XL = 0,  ///< Launch the kernel with the largest possible block size, up to 512 threads per
           ///<   block.
  LG,      ///< Launch the kernel with a large block size, up to 256 threads per block.
  MD,      ///< Launch the kernel with a medium block size, up to 128 threads per block.
  SM       ///< Launch the kernel with the smallest possible block size of 64 threads per block.
};
  
/// \brief Get a human-readable name for the enumerations detailed above.
///
/// \param input  The enumeration of interest
/// \{
std::string getEnumerationName(EvaluateForce input);
std::string getEnumerationName(EvaluateEnergy input);
std::string getEnumerationName(EvaluateVirial input);
std::string getEnumerationName(DihedralStyle input);
std::string getEnumerationName(StateVariable input);
std::string getEnumerationName(NonbondedPotential input);
std::string getEnumerationName(NonbondedTheme input);
std::string getEnumerationName(DecomposablePotential input);
std::string getEnumerationName(VdwCombiningRule input);
std::string getEnumerationName(ClashResponse input);
std::string getEnumerationName(EnergySample input);
std::string getEnumerationName(SplineScaffold input);
std::string getEnumerationName(CellGridAction input);
std::string getEnumerationName(QMapMethod input);
std::string getEnumerationName(PMIStrategy input);
std::string getEnumerationName(ValenceKernelSize input);
/// \}

/// \brief Produce a potential enumeration based on a selection of recognized input strings.
///
/// \param input  The string describing the potential
NonbondedPotential translateNonbondedPotential(const std::string &input);

/// \brief Produce a potential enumeration based on a selection of recognized input strings.
///
/// \param input  The string describing the potential
NonbondedTheme translateNonbondedTheme(const std::string &input);

/// \brief Produce an energy sampling enumeration based on a selection of recognized input strings.
///
/// \param input  The string describing the sampling
EnergySample translateEnergySample(const std::string &input);
  
/// \brief Produce an enumerated accuracy target for the particle-mesh interaction grid based on a
///        selection of recognized input strings.
///
/// \param input  The string describing the accuracy level
PMIStrategy translatePMIStrategy(const std::string &input);
  
} // namespace energy
} // namespace stormm

#endif
