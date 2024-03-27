// -*-c++-*-
#ifndef STORMM_NML_PPPM_H
#define STORMM_NML_PPPM_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/textfile.h"
#include "Potential/cellgrid.h"
#include "Potential/energy_enumerators.h"
#include "input.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using constants::ExceptionResponse;
using energy::default_mesh_ticks;
using energy::maximum_cell_width;
using energy::NonbondedTheme;
using energy::PMIStrategy;
using parse::TextFile;
using parse::WrapTextSearch;

/// \brief Default values for periodic simulations and molecular mechanics calculations 
/// \{
constexpr double max_pp_cutoff = 2.0 * maximum_cell_width;
constexpr NonbondedTheme default_pppm_theme = NonbondedTheme::ELECTROSTATIC;
constexpr PMIStrategy default_pppm_accuracy = PMIStrategy::NO_AUTOMATION;
/// \}

/// \brief Object to encapsulate electrostatic and Lennard-Jones particle-mesh interaction
///        controls.  Like other namelist encapsualtors, this object can take input file data as
///        part of its construction, or by a series of setters.  Validation of each piece of data
///        is handled as it appears either in the contructor or via setters.  Getter functions
///        dispense the internal information to any application using STORMM libraries.
class PPPMControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \param found_nml   Indication that the namelist was found (passed back to calling function)
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &random namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  PPPMControls(ExceptionResponse policy_in = ExceptionResponse::DIE,
               WrapTextSearch wrap = WrapTextSearch::NO);
  PPPMControls(const TextFile &tf, int *start_line, bool *found_nml,
               ExceptionResponse policy_in = ExceptionResponse::DIE,
               WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  PPPMControls(const PPPMControls &original) = default;
  PPPMControls(PPPMControls &&original) = default;
  PPPMControls& operator=(const PPPMControls &original) = default;
  PPPMControls& operator=(PPPMControls &&original) = default;
  /// \}

  /// \brief Get the non-bonded potential theme
  NonbondedTheme getTheme() const;

  /// \brief Get the interpolation order.
  int getInterpolationOrder() const;

  /// \brief Get the Ewald coefficient.
  double getEwaldCoefficient() const;

  /// \brief Get the Gaussian particle density spreading sigma.
  double getGaussianWidth() const;
  
  /// \brief Get the particle-particle interaction cutoff.
  double getCutoff() const;

  /// \brief Get the direct sum tolerance.
  double getDirectSumTolerance() const;

  /// \brief Get the number of particle-mesh interaction grid ticks per spatial decomposition
  ///        cell.
  int getMeshSubdivisions() const;

  /// \brief Get the strategy for choosing particle-mesh interaction grid parameters.  This is not
  ///        passed along to other objects, and is reported for purposes of interrogating the
  ///        contents of an object of this class.
  PMIStrategy getStrategy() const;

  /// \brief Set the non-bonded potential type.
  ///
  /// Overloaded:
  ///   - Set the theme based on an explicit enumeration
  ///   - Set the theme based on a human-intepretable, or user-written, keyword string
  ///
  /// \param theme_in  The theme to apply
  /// \{
  void setTheme(NonbondedTheme theme_in);
  void setTheme(const std::string &theme_in);
  /// \}

  /// \brief Set the interpolation order.
  ///
  /// \param order_in  The interpoation order to set
  void setInterpolationOrder(int order_in);

  /// \brief Set the Ewald coefficient.
  ///
  /// \param ewald_coefficient_in  The value a in, for example, the expression erf(a * r) / r
  void setEwaldCoefficient(double ewald_coefficient_in);

  /// \brief Set the Gaussian particle density spreading sigma.
  void setGaussianWidth(double spread_in);

  /// \brief Set the particle-particle interaction cutoff.
  void setCutoff(double cutoff_in);

  /// \brief Set the direct sum tolerance.
  void setDirectSumTolerance(double dsum_tol_in);

  /// \brief Set the number of particle-interaction grid elements per spatial decomposition cell.
  void setMeshSubdivisions(int mesh_ticks_in);

  /// \brief Set the strategy for determining various parameters of the particle-mesh interaction
  ///        grid.  Overloading follows from other keywords which control an enumeration setting.
  ///
  /// \param strat_in  Rough description of the accuracy level sought in non-bonded calculations
  /// \{
  void setStrategy(PMIStrategy strat_in);
  void setStrategy(const std::string &strat_in);
  /// \}

  /// \brief Apply the current strategy, to the extent possible, for filling in missing parameters
  ///        from the &pppm namelist in the interest of obtaining a particular level of accuracy.
  void applyStrategy();

private:
  ExceptionResponse policy;  ///< Set the behavior when bad inputs are encountered.  DIE = abort
                             ///<   program, WARN = warn the user, and likely reset to the default
                             ///<   value if one is available, SILENT = do not warn the user, but
                             ///<   also likely reset to the default value if one is available.
  NonbondedTheme theme;      ///< Specify whether these PPPM parameters apply to electrostatic
                             ///<   charges or Lennard-Jones sources.
  int order;                 ///< The order of B-spline interpolation used to map particles to the
                             ///<   particle-mesh interaction grid.  This is the "mesh" in many
                             ///<   codes' parlance, but STORMM often refers to it as the
                             ///<   "particle-mesh interaction grid" (see the PMIGrid class) due
                             ///<   to the existence of an alternative "mesh" based on tricubic
                             ///<   interpolation used by its docking applications.
  double ewald_coefficient;  ///< The "Ewald coefficient", a, appearing in the expression
                             ///<   erf(a * r) / r.  This is (1/2) (1 / sigma), where sigma is the
                             ///<   standard deviation of the Gaussian distribution by which the
                             ///<   particle density is smeared in its representation on the mesh.
  double cutoff;             ///< Define the particle-particle interaction cutoff for interactions
                             ///<   of the type defined in theme.  This setting, if provided, is
                             ///<   intended to override settings elsewhere, and is provided
                             ///<   because the particle-particle interaction cutoff is tightly
                             ///<   coupled to the mesh density when seeking a given overall level
                             ///<   of accuracy.  The cutoff will apply only to the type of
                             ///<   interactions indicated by the non-bonded potential theme.
  double dsum_tol;           ///< The "direct sum tolerance", a measure of the ratio of the
                             ///<   difference in the interaction of two spherical Gaussian charges
                             ///<   and that of two point charges when separated by the cutoff.
  double mesh_ticks;         ///< The number of particle-mesh interaction grid points spanning the
                             ///<   spatial decomposition cell in all directions.  This is not
                             ///<   equivalent to nfft in Amber codes, the length of each size of
                             ///<   the particle-mesh interaction grid.  In STORMM, the density /
                             ///<   potential grid length is a multiple of the neighbor list
                             ///<   spatial decomposition grid lengths.  This number mesh_ticks
                             ///<   will be arried over to the program's CellGrid object(s).
  PMIStrategy strat;         ///< Define the strategy for laying out the particle-mesh interaction
                             ///<   grids from a human-comprehensible perspective.  Various inputs
                             ///<   can suggest settings for the number of mesh ticks, cutoff,
                             ///<   and direct sum tolerance as needed to fill in different
                             ///<   settings.
  
  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;

  /// \brief Validate the interpolation order based on minimum requirements for force calculations
  ///        and the upper limits of what the kernels can handle.  Return TRUE if the proposed
  ///        interpolation order is valid, FALSE otherwise.
  ///
  /// \param order_in  The order to validate
  bool validateOrder(const int order_in) const;

  /// \brief Validate the Ewald coefficient (this is also invoked by attempting to set the
  ///        particle density Gaussian width).
  ///
  /// \param ewald_coefficient_in
  bool validateEwaldCoefficient(double ewald_coefficient_in) const;

  /// \brief Validate the particle-particle interaction cutoff.
  ///
  /// \param cutoff_in  The particle-particle interaction cutoff to set
  bool validateCutoff(double cutoff_in) const;

  /// \brief Validate the direct sum tolerance.
  ///
  /// \param tol_in  The direct sum tolerance to validate
  bool validateDirectSumTolerance(double tol_in) const;

  /// \brief Validate the number of particle-mesh interaction grid points spanning each spatial
  ///        decomposition cell.
  ///
  /// \param mesh_ticks_in  The number of mesh ticks spanning each spatial decomposition cell
  bool validateMeshSubdivisions(int mesh_ticks_in) const;
};

/// \brief Produce a namelist for specifying particle-particle / particle-mesh non-bonded
///        interactions.  This can be used for both electrostatics as well as Lennard-Jones
///        potentials, and comprises specifics about the density interpolation on the mesh
///        (referred to in the documentation as "mapping," also called "spreading" in other
///        programs) as well as the Gaussian-based particle smearing.  The concept of PPPM, or
///        "P3M" encompasses a variety of particle smearing functions, but the Gaussians and
///        B-splines used by the form that many call "Particle-Mesh Ewald" (PME) is among the
///        most successful due to the rapid convergence of numerical integration of Gaussians
///        with respect to the grid discretization.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will wrap back to the beginning of the TextFile object, if
///                    necessary, to find a &random namelist)
/// \param found       Indication that the namelist was found (passed back to calling function)
/// \param policy      Reaction to exceptions encountered during namelist reading
/// \param wrap        Indicate that the search for a &random namelist should carry on from the
///                    beginning of an input file if no such namelist is found starting from the
///                    original starting point
NamelistEmulator pppmInput(const TextFile &tf, int *start_line, bool *found,
                           ExceptionResponse policy = ExceptionResponse::DIE,
                           WrapTextSearch wrap = WrapTextSearch::NO);

} // namespace namelist
} // namespace stormm

#endif
