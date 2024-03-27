// -*-c++-*-
#ifndef STORMM_FORCEFIELD_ENUMERATORS_H
#define STORMM_FORCEFIELD_ENUMERATORS_H

#include "copyright.h"

namespace stormm {
namespace modeling {

/// \brief List the possible force field terms, from valence interactions to van-der Waals
///        approximations to charges and virtual site frames.
enum class ParameterKind {
  BOND,               ///< Harmonic bond interaction with two atom types, a stiffness constant,
                      ///<   and an equilibrium length
  ANGLE,              ///< Harmonic angle interaction with three atom types, a stiffness constant,
                      ///<   and an equilibrium angle
  DIHEDRAL,           ///< Cosine-based dihedral term with four atom types, an amplitude, a
                      ///<   phase angle, and a periodicity
  UREY_BRADLEY,       ///< Urey-Bradley angle interaction with three atom types (even though the
                      ///<   type of the middle atom does not explicitly show up in the resulting
                      ///<   topology), a stiffness parameter, and an equilibrium distance
  CHARMM_IMPROPER,    ///< CHARMM harmonic improper dihedral term with four atom types, a stiffness
                      ///<   constant, and a phase angle (almost certainly pi radians)
  CMAP,               ///< Correction MAP 2D bicubic spline surface with five atom and residue
                      ///<   names plus an N x N grid of surface values (the dimension will be
                      ///<   inferred from the square root of the grid size)
  ATTN_14_SCALE,      ///< Attenuated 1:4 non-bonded interaction with four atom types plus charge
                      ///<   and van-der Waals scaling factors
  CHARGE,             ///< Partial atomic charge with atom and residue names plus a magnitude
  LENNARD_JONES,      ///< Lennard-Jones parameter tuple with an atom type, sigma, epsilon, and
                      ///<   possibly rho
  BUCKINGHAM,         ///< Buckingham parameter tuple with an atom type, sigma, epsilon, and rho
                      ///<   constants
  VIRTUAL_SITE_FRAME, ///< Virtual site frame with atom and residue names and three dimensions
  NONE                ///< No interaction
};

/// \brief Get the string interpretation of a ParameterKind enumeration.
///
/// \param kind  The enumeration of interest
std::string getEnumerationName(ParameterKind kind);
  
} // namespace modeling
} // namespace stormm

#endif
