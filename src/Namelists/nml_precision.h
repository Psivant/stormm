// -*-c++-*-
#ifndef STORMM_NML_PRECISION_H
#define STORMM_NML_PRECISION_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "Parsing/textfile.h"
#include "input.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using constants::PrecisionModel;
using parse::WrapTextSearch;
using numerics::default_charge_mesh_scale_bits;
using numerics::default_energy_scale_bits;
using numerics::default_force_scale_bits;
using numerics::default_globalpos_scale_bits;
using numerics::default_localpos_scale_bits;
using numerics::default_velocity_scale_bits;

/// \brief While the bit counts for fixed-precision accumulation and coordinate representations
///        are held in the fixed_precision.h header file (see src/Constants/), additional default
///        settings are listed here.
/// \{
constexpr char default_precision_valence_method[]   = "single";
constexpr char default_precision_nonbonded_method[] = "single";
constexpr char default_precision_pme_method[]       = "single";
constexpr double default_precision_constraint_tol   = 1.0e-5;
constexpr double min_precision_constraint_tol       = 5.0e-9;
/// \}

/// \brief Object to encapsulate energy precision control information.  Like other namelist
///        encapsualtors, this object can take input file data as part of its construction, or
///        by a series of setters.  Validation of each piece of data is handled as it appears
///        either in the contructor or via setters.  Getter functions dispense the internal
///        information to any application using STORMM libraries.
class PrecisionControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \param found_nml   Indication of whether the namelist was found in the input file
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &conformer namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  PrecisionControls(ExceptionResponse policy_in = ExceptionResponse::DIE,
                    WrapTextSearch wrap = WrapTextSearch::NO);
  PrecisionControls(const TextFile &tf, int *start_line, bool *found_nml,
                    ExceptionResponse policy_in = ExceptionResponse::DIE,
                    WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  PrecisionControls(const PrecisionControls &original) = default;
  PrecisionControls(PrecisionControls &&original) = default;
  PrecisionControls& operator=(const PrecisionControls &original) = default;
  PrecisionControls& operator=(PrecisionControls &&original) = default;
  /// \}

  /// \brief Get the global position fixed-precision bits.
  int getGlobalPosScalingBits() const;

  /// \brief Get the local position fixed-precision bits.
  int getLocalPosScalingBits() const;

  /// \brief Get the velocity fixed-precision bits.
  int getVelocityScalingBits() const;

  /// \brief Get the force accumulation fixed-precision bits.
  int getForceScalingBits() const;

  /// \brief Get the energy accumulation fixed-precision bits.
  int getEnergyScalingBits() const;

  /// \brief Get the PME charge accumulation fixed-precision bits.
  int getChargeMeshScalingBits() const;

  /// \brief Get the SHAKE / RATTLE convergence tolerance
  double getBondConstraintTolerance() const;
  
  /// \brief Get the precision level for valence computations.
  PrecisionModel getValenceMethod() const;

  /// \brief Get the precision level for short-ranged, non-bonded computations.
  PrecisionModel getNonbondedMethod() const;

  /// \brief Get the precision level for Particle-Mesh Ewald computations.
  PrecisionModel getParticleMeshEwaldMethod() const;

  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator& getTranscript() const;
  
  /// \brief Set the global position fixed-precision bits.
  ///
  /// \param bitval  The number of bits after the decimal (will be checked for validity in light
  ///                of hard-coded limits)
  void setGlobalPosScalingBits(int bitval);

  /// \brief Set the local position fixed-precision bits.
  ///
  /// \param bitval  The number of bits after the decimal (will be checked for validity in light
  ///                of hard-coded limits)
  void setLocalPosScalingBits(int bitval);

  /// \brief Set the fixed-precision velocity scaling bits.
  ///
  /// \param bitval  The number of bits after the decimal (will be checked for validity in light
  ///                of hard-coded limits)
  void setVelocityScalingBits(int bitval);

  /// \brief Set the fixed-precision force accumulation bits.
  ///
  /// \param bitval  The number of bits after the decimal (will be checked for validity in light
  ///                of hard-coded limits)
  void setForceScalingBits(int bitval);
  
  /// \brief Set the fixed-precision energy accumulation bits.
  ///
  /// \param bitval  The number of bits after the decimal (will be checked for validity in light
  ///                of hard-coded limits)
  void setEnergyScalingBits(int bitval);
  
  /// \brief Set the fixed-precision charge mesh accumulation bits.
  ///
  /// \param bitval  The number of bits after the decimal (will be checked for validity in light
  ///                of hard-coded limits and charge mesh filling methods)
  void setChargeMeshScalingBits(int bitval);

  /// \brief Set the bond constraint convergence tolerance.
  ///
  /// param tol  The tolerance to which SHAKE- or RATTLE-constrained bonds must agree with their
  ///            equilibrium lengths
  void setBondConstraintTolerance(double tol);

  /// \brief Set the precision model to use in valence term computations.
  ///
  /// param pmodel  The requested precision model
  void setValenceMethod(PrecisionModel pmodel);

  /// \brief Set the precision model to use in non-bonded short-ranged pair computations.
  ///
  /// param pmodel  The requested precision model
  void setNonbondedMethod(PrecisionModel pmodel);

  /// \brief Set the precision model to use in Particle Mesh Ewald (charge density accumulation,
  ///        Fast Fourier Transforms, convolution, and force interpolation) computations.
  ///
  /// param pmodel  The requested precision model
  void setParticleMeshEwaldMethod(PrecisionModel pmodel);

private:
  ExceptionResponse policy;        ///< Set the behavior when bad inputs are encountered.  DIE =
                                   ///<   abort program, WARN = warn the user, and likely reset to
                                   ///<   the default value if one is available, SILENT = do not
                                   ///<   warn the user, but also likely reset to the default value
                                   ///<   if one is available.
  int globalpos_scale_bits;        ///< Global position scaling bits
  int localpos_scale_bits;         ///< Local position scaling bits
  int velocity_scale_bits;         ///< Particle velocity scaling bits
  int force_scale_bits;            ///< Force accumulation scaling bits
  int energy_scale_bits;           ///< Energy accumulation scaling bits
  int charge_mesh_scale_bits;      ///< Charge mesh accumulation scaling bits
  double bond_constraint_tol;      ///< Tolerance for bond constraints (tol in sander and pmemd)
  PrecisionModel valence_method;   ///< Valence term computation precision model
  PrecisionModel nonbonded_method; ///< Non-bonded short-ranged computation precision model
  PrecisionModel pme_method;       ///< Particle-Mesh Ewald computation precision model

  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;
  
  /// \brief Validate the bond constraint tolerance.  It cannot be too small, lest calculations
  ///        become unstable.
  void validateBondConstraintTol();
};
  
/// \brief Produce a namelist for specifying precision model directives, similar to those found
///        in the &cntrl namelist of sander or pmemd.  This namelist can be included as part of
///        many computations.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will wrap back to the beginning of the TextFile object, if needed,
///                    to find a &minimize namelist) 
/// \param found       Indicator that the namelist was found in the input file
/// \param policy      Reaction to exceptions encountered during namelist reading
/// \param wrap        Indicate that the search for a &precision namelist should carry on from the
///                    beginning of an input file if no such namelist is found starting from the
///                    original starting point
NamelistEmulator precisionInput(const TextFile &tf, int *start_line, bool *found,
                                ExceptionResponse policy = ExceptionResponse::DIE,
                                WrapTextSearch wrap = WrapTextSearch::NO);

} // namespace namelist
} // namespace stormm

#endif
