// -*-c++-*-
#ifndef STORMM_ATOMGRAPH_ENUMERATORS_H
#define STORMM_ATOMGRAPH_ENUMERATORS_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/citation.h"

namespace stormm {
namespace topology {

using constants::CaseSensitivity;
using constants::ExceptionResponse;
using parse::Citation;

/// \brief Enumerate different kinds of topologies that STORMM will read.  Each comes with its own
///        defaults for Coulomb's constant and the default 1:4 screening factors.
enum class TopologyKind {
  AMBER,    ///< Amber format, Coulomb's constant C = 332.05221729, van-der Waals screening 2.0,
            ///<   electrostatic screening 1.2
  CHARMM,   ///< CHARMM PSF format, Coulomb's constant C = 332.063711, van-der Waals screening 1.0,
            ///<   electrostatic screening 1.0
  GROMACS,  ///< GROMACS format, Coulomb's constant C = 332.063711, van-der Waals screening 1.0,
            ///<   electrostatic screening 1.0
  OPENMM    ///< OPENMM format, Coulomb's constant C = 332.063711, van-der Waals screening 1.0,
            ///<   electrostatic screening 1.0
};

/// \brief An enumerator to list the various integer attributes of an Amber topology from its
///        preamble.
enum class TopologyDescriptor {
  ATOM_COUNT = 0, ATOM_TYPE_COUNT, BONDS_WITH_HYDROGEN, BONDS_WITHOUT_HYDROGEN,
  ANGLES_WITH_HYDROGEN, ANGLES_WITHOUT_HYDROGEN, DIHEDRALS_WITH_HYDROGEN,
  DIHEDRALS_WITHOUT_HYDROGEN, NHPARM_UNUSED, ADDLES_CREATED, TOTAL_EXCLUDED_ATOMS,
  RESIDUE_COUNT, NBONA_UNUSED, NTHETA_UNUSED, NPHIA_UNUSED, BOND_TYPE_COUNT,
  ANGLE_TYPE_COUNT, DIHEDRAL_TYPE_COUNT, NATYP_UNUSED, NPHB_UNUSED, PERTURBATION,
  BOND_PERTURBATIONS, ANGLE_PERTURBATIONS, DIHEDRAL_PERTURBATIONS, BONDS_IN_PERTURBED_GROUP,
  ANGLES_IN_PERTURBED_GROUP, DIHEDRALS_IN_PERTURBED_GROUP, BOX_TYPE_INDEX,
  ATOM_COUNT_LARGEST_RESIDUE, CAP, EXTRA_POINT_COUNT, PIMD_SLICE_COUNT, N_VALUES
};

/// \brief Secondary enumerator to list all of the previous enumerators in terms of their sander
///        variable names, which may be more familiar to developers.  While the previous enumerator
///        may be more performant, the frequency of accessing most of these delimiters is low so
///        either method should be effective.
enum class SanderDescriptor {
  NATOM = 0, NTYPES, NBONH, MBONA, NTHETH, MTHETA, NPHIH, MPHIA, NHPARM,
  NPARM, NNB, NRES, NBONA, NTHETA, NPHIA, NUMBND, NUMANG, NPTRA, NATYP, NPHB, IFPERT,
  NBPER, NGPER, NDPER, MBPER, MGPER, MDPER, IFBOX, NMXRS, IFCAP, NUMEXTRA, NCOPY, N_VALUES
};

/// \brief List the basic unit cell types (this gives some more color to sander's IFBOX term)
enum class UnitCellType {
  NONE = 0, ORTHORHOMBIC, TRICLINIC
};

/// \brief Enumerate the ways to modify an atom's mobility in the toplogy
enum class MobilitySetting {
  OFF, ON, TOGGLE
};

/// \brief An enumerator to toggle SHAKE and RATTLE, depending on the integrator at hand.  It's
///        all bond length constraints.
enum class ShakeSetting {
  OFF, ON
};

/// \brief An enumerator to toggle SETTLE for rigid water molecules (or perhaps other things with
///        a central atom connected to two other atoms of identical mass connected by bonds of
///        equal length).
enum class SettleSetting {
  OFF, ON
};

/// \brief An enumerator to toggle free energy perturbations.
enum class PerturbationSetting {
  OFF, ON
};

/// \brief An enumerator to toggle the Amber prmtop solvent cap feature (very old)
enum class SolventCapSetting {
  OFF, ON
};

/// \brief An enumerator to toggle polarization
enum class PolarizationSetting {
  OFF, ON
};

/// \brief Indicate whether a field of a topology is a requirement or if it can be skipped
enum class TopologyRequirement {
  OPTIONAL, ESSENTIAL
};

/// \brief Offer a choice between retrieving masses or inverse masses
enum class MassForm {
  ORDINARY, INVERSE
};

/// \brief Valence term modifiers / attributes.  These will occupy the x, y, z, and (for dihedrals)
///        the w components of a char4 vector.  Each attribute holds one of up to 128 different
///        settings, offering great flexibility and possiblities for making new types of valence
///        interactions.
/// \{
enum class ConstraintStatus : char {
  FREE, CONSTRAINED
};
enum class HydrogenContent : char {
  NO_HYDROGEN, HAS_HYDROGEN
};
enum class ForceFieldFamily : char {
  BASIC, AMBER, CHARMM, OPENMM
};
enum class TorsionKind : char {
  PROPER, PROPER_NO_14, IMPROPER, IMPROPER_NO_14
};
/// \}

/// \brief Enumerate the multitude of available Generalized Born and other imnplicit solvent
///        models that gas-phase, isolated boundary conditions can implement.
enum class ImplicitSolventModel {
  NONE,       ///< No GB model, equivalent to gas-phase Coulomb interactions (igb = 6 in sander)
  HCT_GB,     ///< Hawkins / Cramer / Truhlar Generalized Born (igb = 1 in sander)
  OBC_GB,     ///< Onufriev / Bashford / Case Generalized Born (igb = 2 in sander)
  OBC_GB_II,  ///< Onufriev / Bashford / Case Generalized Born (igb = 5 in sander)
  NECK_GB,    ///< Mongan's neck Generalized Born (igb = 7 in sander)
  NECK_GB_II, ///< Mongan's neck Generalized Born, model II (igb = 8 in sander)
};

/// \brief Enumerate the radii sets that can, in some cases, be mixed and matched with implicit
///        solvent models (some models require a specific set, or at least place limits on the
///        largest and smallest radii one can use).
enum class AtomicRadiusSet {
  NONE, BONDI, AMBER6, MBONDI, MBONDI2, MBONDI3, PARSE
};

/// \brief Enumerate commonly used water models
enum class WaterModel {

  // Default and fallback states
  NONE,           ///< The system has no identifiable water molecules with formula H2O (or D2O)
  UNKNOWN,        ///< The water model could not be identified from anything in this list
  CHIMERA,        ///< The water model appears to be chimeric, probably indicating an error
  MULTIPLE,       ///< There are multiple water models in play
  MULTI_CHIMERA,  ///< There are multiple water models in play, and at least one is chimeric

  // Three-site, rigid water models
  OPC3,           ///< OPC-3 Point, S. Izadi and A. Onufriev, J. Chem. Phys. 145:074501 (2016)
  SPC,            ///< SPC, HJC Berendsen et al., "Intermolecular Forces" pp. 331 (Springer,
                  ///<   Dordrecht, 1981)
  SPC_E,          ///< SPC/E, HJC Berendsen et al., J. Phys. Chem. 91:6269-6271 (1987)
  SPC_EB,         ///< SPC/E-b, K Takemura and A Kitao, J. Phys. Chem. B 116:6279-6287 (2012)
  SPC_HW,         ///< SPC/Heavy Water, JR Grigera, J. Chem. Phys. 114:8064-8067 (2001)
  TIP3P,          ///< TIP3P, WL Jorgensen et al., J. Chem. Phys. 79:926-935 (1983)
  TIP3P_EW,       ///< TIP3P-Ewald (model F), DL Price and CL Brooks, J. Chem. Phys.
                  ///<   121:10096-10103 (2004)
  TIP3P_CHARMM,   ///< CHARMM implementation of TIP3P, with nonzero Hydrogen Lennard-Jones

  // Three-site flexible water models
  SPC_FW,         ///< SPC-Flexible Water, Y Wu et al., J. Chem. Phys. 124:024503 (2006)
  TIP3P_FW,       ///< TIP3P-Flexible Water, Y Wu et al., J. Chem. Phys. 124:024503 (2006)

  // Four-site, rigid water models
  OPC,            ///< OPC, S Izadi et al., J. Phys. Chem. Lett. 5:3863-3871 (2014)
  TIP4P,          ///< TIP4P, WL Jorgensen and JD Madura, Mol. Phys. 56:1381-1392 (1985)
  TIP4P_EW,       ///< TIP4P-Ewald, HW Horn et al., J. Chem. Phys. 120:9665-9678 (2004)
  TIP4P_2005,     ///< TIP4P-2005, JLF Abascal and C Vega, J. Chem. Phys. 123:234505 (2005)
  TIP4P_ICE,      ///< TIP4P/Ice, JLF Abascal et al., J. Chem. Phys. 122:234511 (2005)
  TIP4P_EPS,      ///< TIP4P/Epsilon, R Fuentes-Azcatl and MC Barbosa, Physica A 444:86-94 (2016)
  TIP4P_D,        ///< TIP4P-D, Piana et al., J. Phys. Chem. B 119:5113-5123 (2015)

  // Four-site, flexible water models
  TIP4P_2005F,    ///< TIP4P-2005, MA Gonzalez and JLF Abascal, J. Chem. Phys. 135:224516 (2011)

  // Five-site, rigid water models
  TIP5P,          ///< TIP5P, MW Mahoney and WL Jorgensen, J. Chem. Phys. 112:8910-8922 (2000)
  TIP5P_EW,       ///< TIP5P-Ewald, SW Rick, J. Chem. Phys. 120:6085-6093 (2004)
  TIP5P_2018,     ///< TIP5P-2018, Y Khalak et al., J. Chem. Phys. 149:224507 (2018)
};

/// \brief Enumerate the virtual site types available in Amber.
enum class VirtualSiteKind {
  NONE = 0, ///< No frame type
  FLEX_2,   ///< Flexible two-atom frame: the distance between the virtual site and its parent
            ///<   atom scales with the distance between the parent atom and frame atom 2
  FIXED_2,  ///< Fixed distance, two-atom frame: the distance between the virtual site and its
            ///<   parent atom is fixed regardless of the way the frame stretches
  FLEX_3,   ///< Flexible three-atom frame: the distance between the virtual site and its parent
            ///<   atom scales with the distance between the parent atom and frame atoms 2 and 3
  FIXED_3,  ///< Fixed distance, three-atom frame: the distance between the virtual site and its
            ///<   parent atom is fixed, along a line between the parent atom and a point between
            ///<   frame atoms 2 and 3 that does stretch with the distance between those atoms.
            ///<   Mathematically, this is the equivalent of the FIXED_2 frame expanded to a third
            ///<   dimension.
  FAD_3,    ///< Fixed distance, fixed angle three-atom frame: the virtual site is placed at a
            ///<   a fixed distance from its parent atom, such that the virtual site makes a fixed
            ///<   angle with its parent atom and frame atom 2, the orientation of the angle being
            ///<   by the position of frame atom 3.
  OUT_3,    ///< Out of plane, flexible three-atom frame: the virtual site is placed at a point
            ///<   determined in the manner of FLEX_3, then moved out of plane by some proportion
            ///<   of the cross product of the vectors between the parent atom and frame atoms 2
            ///<   or 3.
  FIXED_4   ///< Fixed distance, four-atom frame: this places the virtual site at a set distance
            ///<   from its parent atom, along a vector determined by a cross product of vectors
            ///<   between frame atoms 2, 3, and 4
};
  
/// \brief Produce a human-readable string corresponding to the enumeration of interest.  Various
///        overloads of this function (in this and other libaries and namespaces) serve different
///        enum class objects.
///
/// \param input  The enumeration of interest
/// \{
std::string getEnumerationName(TopologyKind input); 
std::string getEnumerationName(UnitCellType input); 
std::string getEnumerationName(MobilitySetting input);
std::string getEnumerationName(ShakeSetting input);
std::string getEnumerationName(SettleSetting input);
std::string getEnumerationName(PerturbationSetting input);
std::string getEnumerationName(SolventCapSetting input);
std::string getEnumerationName(PolarizationSetting input);
std::string getEnumerationName(TopologyRequirement input);
std::string getEnumerationName(MassForm input);
std::string getEnumerationName(ConstraintStatus input);
std::string getEnumerationName(HydrogenContent input);
std::string getEnumerationName(ForceFieldFamily input);
std::string getEnumerationName(TorsionKind input);
std::string getEnumerationName(ImplicitSolventModel input);
std::string getEnumerationName(AtomicRadiusSet input);
std::string getEnumerationName(WaterModel input);
std::string getEnumerationName(VirtualSiteKind input);
/// \}

/// \brief Translate the numerical input for the implicit solvent model into one of the recognized
///        models available in STORMM.  This serves as the validator to the implicit solvent model
///        in the SolventControls object.
///
/// \param igb_val  Numerical value of the implicit solvent model, as read from the parent namelist
/// \param policy   Set the response to bad inputs
ImplicitSolventModel
translateImplicitSolventModel(int igb_val, ExceptionResponse policy = ExceptionResponse::DIE);

/// \brief Translate string input into one of the enumerated AtomicRadiusSet values.  This serves
///        as the validator to the radius set in the SolventControls object.
///
/// \param pb_radii_in  Input string describing the PB radii (case-insensitive)
/// \param policy       Set the response to bad inputs
AtomicRadiusSet translateAtomicRadiusSet(const std::string &pb_radii_in,
                                         ExceptionResponse policy = ExceptionResponse::DIE);

/// \brief Produce the most relevant citation for a water model based on its enumeration.
///
/// \param wm  The water model enumerated value
Citation getWaterModelCitation(WaterModel wm);
  
} // namespace topology
} // namespace stormm

#endif
