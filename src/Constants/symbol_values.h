// -*-c++-*-
#ifndef STORMM_SYMBOLS_H
#define STORMM_SYMBOLS_H

#include "copyright.h"

namespace stormm {
namespace symbols {

/// \brief Pi and related quantities, represented in values that can be stored exactly in
///        double-precision floating point numbers according to the IEEE_754 standard.
constexpr double pi = 3.141592653589793115997963468544185161590576171875;
constexpr double twopi = 2.0 * pi;
constexpr double inverse_pi = 1.0 / pi;
constexpr double inverse_twopi = 1.0 / twopi;
constexpr float pi_f = 3.1415927410125732421875f;
constexpr float twopi_f = 6.283185482025146484375f;
constexpr float inverse_pi_f = 0.3183098733425140380859375f;
constexpr float inverse_twopi_f = 0.15915493667125701904296875f;

/// \brief The base of the natural logarithm, e.
constexpr double evalue = 2.718281828459045090795598298427648842334747314453125;
constexpr float evalue_f = evalue;
  
/// The tetrahedral bond angle, in radians to twenty decimal places
constexpr double tetrahedral_angle = 1.91063323624901859610;

/// The "biological charge" conversion, taking atomic units of charge, Coulomb's constant, and
/// Angstroms into kilocalories per mole.  Amber and CHARMM disagree on this value.  This is the
/// energy, in kcal/mol, of a pair of charges with 1 a.u. charge being brought from infinity to
/// within 1 Angstrom of each other.
/// \{
constexpr double amber_ancient_bioq  = 332.0522172900000;
constexpr double charmm_gromacs_bioq = 332.0636974382250;
/// \}

/// \brief Convert energy quantities computed with internal units of A, fs, and g/mol into
///        kcal/mol.
constexpr double gafs_to_kcal  = 1.0 / 0.0004184;
constexpr float gafs_to_kcal_f = gafs_to_kcal;
constexpr double kcal_to_gafs = 0.0004184;
constexpr float kcal_to_gafs_f = kcal_to_gafs;
  
/// \brief Avogadro's number
constexpr double avogadro_number = 6.02214076e+23;
  
/// \brief Boltzmann's constant in kcal/mol-K
constexpr double boltzmann_constant = (1.38064852e-23) / 4184.0 * avogadro_number;
constexpr double boltzmann_constant_gafs = boltzmann_constant * kcal_to_gafs;
constexpr float boltzmann_constant_f = boltzmann_constant;
constexpr float boltzmann_constant_gafs_f = boltzmann_constant_gafs;
  
/// \brief Hartree to kcal/mol conversion
constexpr double hartree_to_kcal = 627.509474;

/// \brief Bohr to Angstrom conversion factor
constexpr double bohr_to_angstrom = 0.529177210903;

/// \brief Angstrom to Bohr conversion factor
constexpr double angstrom_to_bohr = 1.889726124626;

/// \brief Values which approach one from below.  They are used in dihedral and similar
///        computations to detect when a value is nearing 1.0 and might generate a singularity
///        in some denominator, or otherwise become numerically ill-conditioned.  The first is
///        nearly the closest value to 1.0 that can be represented in an IEEE-754 format 32-bit
///        floating point number (the format can go a couple of bits' worth of precision closer,
///        but this is 1 part in about a million).  The second is essential for single-precision
///        dihedral computations to guard against instability in the arccosine function.  A
///        long-float (double-precision form of the same near-to-one limit is provided to
///        safeguard double-precision dihedral computations.
/// \{
constexpr double asymptotic_to_one_lf = 0.99999904632568359375;
constexpr float  asymptotic_to_one_f  = (float)asymptotic_to_one_lf;
constexpr float  near_to_one_f        = 0.99993896484375f;
constexpr float  near_to_one_lf       = 0.999999992549419403076171875;
/// \}

/// \brief A value which captures 1 / (1 - asymptotic_to_one), to put a cap on the value of such
///        fractions 1 / (1 - x) as x -> 1.
/// \{
constexpr double inverse_one_minus_asymptote_lf = 1048576.0;
constexpr float  inverse_one_minus_asymptote_f = (float)1048576.0;
/// \}

} // namespace constants
} // namespace stormm

// Put common symbols in the stormm namespace
namespace stormm {
using symbols::pi;
using symbols::twopi;
using symbols::inverse_pi;
using symbols::inverse_twopi;
using symbols::pi_f;
using symbols::twopi_f;
using symbols::inverse_pi_f;
using symbols::inverse_twopi_f;
using symbols::evalue;
using symbols::evalue_f;
} // namespace stormm

#endif
