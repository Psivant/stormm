// -*-c++-*-
#ifndef STORMM_NML_SOLVENT_H
#define STORMM_NML_SOLVENT_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/textfile.h"
#include "Topology/atomgraph_enumerators.h"
#include "input.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using constants::ExceptionResponse;
using parse::WrapTextSearch;
using topology::ImplicitSolventModel;
using topology::AtomicRadiusSet;

/// \brief Default values for the implicit solvent model
/// \{
constexpr int default_solvent_igb = 0;
constexpr double default_solvent_rgbmax = 1000.0;
constexpr double default_solvent_intdiel = 1.0;
constexpr double default_solvent_extdiel = 78.5;
constexpr double default_solvent_saltcon = 0.0;
constexpr char default_solvent_pbradii[] = "none";
/// \}

/// \brief Object to encapsulate the data that can be extracted from a &solvent namelist.  Typical
///        C++ construction.  Having this object, like others accompanying their respective
///        namelists, makes the transition between namelists and custom-formatted control objects
///        much simpler.
class SolventControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM 
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \param found_nml   Indication that the namelist was found
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &solvent namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  SolventControls(ExceptionResponse policy_in = ExceptionResponse::DIE,
                  WrapTextSearch wrap = WrapTextSearch::NO);
  SolventControls(const TextFile &tf, int *start_line, bool *found_nml,
                  ExceptionResponse policy_in = ExceptionResponse::DIE,
                  WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  SolventControls(const SolventControls &original) = default;
  SolventControls(SolventControls &&original) = default;
  SolventControls& operator=(const SolventControls &original) = default;
  SolventControls& operator=(SolventControls &&original) = default;
  /// \}

  /// \brief Get the type of implicit solvent (some flavor of Generalized Born)
  ImplicitSolventModel getImplicitSolventModel() const;

  /// \brief Get the Born radius calculation cutoff
  double getBornRadiiCutoff() const;

  /// \brief Get the internal dielectric constant
  double getInternalDielectric() const;

  /// \brief Get the external dielectric constant
  double getExternalDielectric() const;

  /// \brief Get the salt concentration
  double getSaltConcentration() const;

  /// \brief Get the Poisson-Boltzmann radii set
  AtomicRadiusSet getPBRadiiSet() const;

  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator& getTranscript() const;
  
  /// \brief Set the implicit solvent model
  ///
  /// Overloaded:
  ///   - Take the integer code for the implicit solvent model (and apply validation)
  ///   - Take the enumeration for the implicit solvent model
  ///
  /// \param ism_in  Implicit solvent model of choice
  /// \{
  void setImplicitSolventModel(int ism_in);
  void setImplicitSolventModel(ImplicitSolventModel ism_in);
  /// \}

  /// \brief Set the Born radius calculation cutoff
  ///
  /// \param rgbmax_in  The Born radius calculation cutoff of choice
  void setBornRadiiCutoff(double rgbmax_in);

  /// \brief Set the internal dielectric
  ///
  /// \param idiel_in  The internal dielectric to apply
  void setInternalDielectric(double idiel_in);
  
  /// \brief Set the external dielectric
  ///
  /// \param ediel_in  The external dielectric to apply
  void setExternalDielectric(double ediel_in);

  /// \brief Set the salt concentration, the concentration of 1:1 monovalent ion pairs in the
  ///        implicit solvent.
  ///
  /// \param saltcon_in  The concentration to set
  void setSaltConcentration(double saltcon_in);

  /// \brief Choose the Poisson-Boltzmann radii set
  ///
  /// Overloaded:
  ///   - Use string input (and translate the string, checking for validity)
  ///   - Use enumeration input
  ///
  /// \param pbrad_in  The radii set to select (will be checked for validity)
  /// \{
  void choosePBRadiiSet(const std::string &pbrad_in);
  void choosePBRadiiSet(AtomicRadiusSet pbrad_in);
  /// \}
  
private:
  ExceptionResponse policy;       ///< Set the behavior when bad inputs are encountered.  DIE =
                                  ///<   abort program, WARN = warn the user, and likely reset to
                                  ///<   the default value if one is available, SILENT = do not
                                  ///<   warn the user, but also likely reset to the default value
                                  ///<   if one is available.
  ImplicitSolventModel gb_style;  ///< The type of Generalized Born or other implicit solvent
                                  ///<   model to use
  double born_radii_cutoff;       ///< Cutoff for pairs of atoms to participate in Generalized
                                  ///<   Born radii (equivalent to rgbmax in sander)
  double internal_dielectric;     ///< Set the internal dielectric constant for all molecules (in
                                  ///<   implicit solvent situations), or the dielectric constant
                                  ///<   for all electrostatic interactions in an explicit
                                  ///<   solvent simulation.
  double external_dielectric;     ///< Set the external diectric for continuum solvent (implicit
                                  ///<   solvent applications only)
  double salt_concentration;      ///< Sets the salt concentration for all molecules (applicable
                                  ///<   to implicit solvent situations only)
  AtomicRadiusSet pb_radii;       ///< Defines the Poisson-Boltzmann radii set (and baseline GB
                                  ///<   radii).  Acceptable values in the input include "Bondi",
                                  ///<   "Amber6", "mBondi", "mBond2", "mBondi3", or "none".

  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;
  
  /// \brief Validate the Born radii cutoff.
  void validateBornRadiiCutoff();

  /// \brief Validate the internal dielectric constant.
  void validateInternalDielectric();

  /// \brief Validate the external dielectric constant.
  void validateExternalDielectric();

  /// \brief Validate the salt concentration.
  void validateSaltConcentration();
};

/// \brief Produce a namelist for defining the implicit solvent model, replicating various inputs
///        in the &cntrl namelist of sander or pmemd.  This is a separate namelist from the
///        molecular dynamics input as well as the Particle-Mesh Ewald namelist.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will wrap back to the beginning of the TextFile object, if needed,
///                    to find an &solvent namelist)
/// \param found       Indication that the namelist was found (passed back to the calling function)
/// \param policy      Reaction to exceptions encountered during namelist reading
/// \param wrap        Indicate that the search for a &solvent namelist should carry on from the
///                    beginning of an input file if no such namelist is found starting from the
///                    original starting point
NamelistEmulator solventInput(const TextFile &tf, int *start_line, bool *found,
                              ExceptionResponse policy = ExceptionResponse::DIE,
                              WrapTextSearch wrap = WrapTextSearch::NO);

} // namespace namelist
} // namespace stormm

#endif
