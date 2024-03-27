// -*-c++-*-
#ifndef STORMM_NML_FFMORPH_H
#define STORMM_NML_FFMORPH_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "ForceField/forcefield_element.h"
#include "ForceField/forcefield_enumerators.h"
#include "Parsing/textfile.h"
#include "input.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using modeling::ForceFieldElement;
using modeling::ParameterKind;
using parse::WrapTextSearch;

/// \brief Object to encapsulate force field morphing operations.  Take information from an input
///        file or a series of setters and validate each piece of data as it appears with private
///        member functions.
class FFMorphControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &ffmorph namelist
  /// \param found_nml   Indicator of whether namelist input was found
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for an &ffmorph namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  FFMorphControls(ExceptionResponse policy_in = ExceptionResponse::DIE,
                  WrapTextSearch wrap = WrapTextSearch::NO);
  FFMorphControls(const TextFile &tf, int *start_line, bool *found_nml,
                  ExceptionResponse policy_in = ExceptionResponse::DIE,
                  WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  FFMorphControls(const FFMorphControls &original) = default;
  FFMorphControls(FFMorphControls &&original) = default;
  FFMorphControls& operator=(const FFMorphControls &original) = default;
  FFMorphControls& operator=(FFMorphControls &&original) = default;
  /// \}

  /// \brief Get the number of edits for a particular force field term.
  ///
  /// \param kind  The type of force field parameter edit of interest
  int getEditCount(ParameterKind kind) const;

  /// \brief Get a specific edit from within the namelist's holdings.
  ///
  /// \param kind   The type of force field parameter edit of interest
  /// \param index  Index of the edit from the appropriate list
  ForceFieldElement getModelEdit(ParameterKind kind, int index) const;

  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator& getTranscript() const;
  
private:
  ExceptionResponse policy;       ///< Set the behavior when bad inputs are encountered.  DIE =
                                  ///<   abort program, WARN = warn the user, and likely reset to
                                  ///<   the default value if one is available, SILENT = do not
                                  ///<   warn the user, but also likely reset to the default value
                                  ///<   if one is available.
  
  /// Harmonic bond parameters
  std::vector<ForceFieldElement> harmonic_bonds;

  /// Harmonic angle parameters, three-point angle bending with harmonic penalties based on
  /// deviations of the I-J-K angle from an ideal value
  std::vector<ForceFieldElement> harmonic_angles;

  /// Cosine-based dihedral parameters (proper or improper)
  std::vector<ForceFieldElement> cosine_dihedrals;

  /// Urey-Bradley angle parameters, two-point stretching terms with harmonic penalties based on
  /// deviations of the I-K angle from an ideal separation
  std::vector<ForceFieldElement> urey_bradley_angles;  

  /// CHARMM improper dihedral parameters
  std::vector<ForceFieldElement> charmm_impropers;

  /// CMAP bicubic spline surface terms
  std::vector<ForceFieldElement> cmap_surfaces;

  /// Attenuated 1:4 interaction scaling factor pairs
  std::vector<ForceFieldElement> attn14_scalings;

  /// Charge parameters
  std::vector<ForceFieldElement> charge_properties;

  /// Lennard-Jones parameters, or other van-der Waals approximations
  std::vector<ForceFieldElement> van_der_waals_properties;

  /// Virtual site frames--the non-bonded properties of the virtual site are controlled by
  /// entries in van_der_waals_properties and charge_properties, respectively.
  std::vector<ForceFieldElement> virtual_sites;

  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;
};

/// \brief Produce a namelist for translating user input into lists of ForceFieldElement objects,
///        which can then be applied to topologies to experiment with new parameter settings.
///
/// \param tf          Input text file to scan immediately after the namelist is created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    will wrap back to the beginning of the file in search of a unique
///                    &files namelist)
/// \param found       Indicator that the namelist was detected in the input file (returned)
/// \param policy      Reaction to exceptions encountered during namelist reading
/// \param wrap        Indicate that the search for an &ffmorph namelist should carry on from the
///                    beginning of an input file if no such namelist is found starting from the
///                    original starting point
NamelistEmulator ffmorphInput(const TextFile &tf, int *start_line, bool *found,
                              ExceptionResponse policy = ExceptionResponse::DIE,
                              WrapTextSearch wrap = WrapTextSearch::NO);
  
} // namespace namelist
} // namespace stormm

#endif
