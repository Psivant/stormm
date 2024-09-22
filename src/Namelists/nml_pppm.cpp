#include <cmath>
#include "copyright.h"
#include "DataTypes/common_types.h"
#include "Parsing/parse.h"
#include "Potential/pme_util.h"
#include "Reporting/error_format.h"
#include "namelist_element.h"
#include "nml_pppm.h"

namespace stormm {
namespace namelist {

using energy::default_pme_cutoff;
using energy::default_dsum_tol;
using energy::default_pme_grid_spacing_target;
using energy::default_charge_mapping_order;
using energy::ewaldCoefficient;
using energy::maximum_ewald_coefficient;
using energy::minimum_ewald_coefficient;
using energy::translateNonbondedTheme;
using energy::translatePMIStrategy;
using parse::NumberFormat;
using parse::realToString;

//-------------------------------------------------------------------------------------------------
PPPMControls::PPPMControls(const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    policy{policy_in}, theme{default_pppm_theme}, order{default_charge_mapping_order},
    ewald_coefficient{0.0}, cutoff{default_pme_cutoff}, dsum_tol{default_dsum_tol},
    mesh_ticks{default_mesh_ticks}, strat{PMIStrategy::NO_AUTOMATION}, nml_transcript{"pppm"}
{}

//-------------------------------------------------------------------------------------------------
PPPMControls::PPPMControls(const TextFile &tf, int *start_line, bool *found_nml,
                           const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    PPPMControls(policy_in)
{
  // Create the namelist object for the &pppm namelist
  NamelistEmulator t_nml = pppmInput(tf, start_line, found_nml, policy, wrap);
  nml_transcript = t_nml;

  // Translate input values
  setTheme(t_nml.getStringValue("theme"));
  t_nml.assignVariable(&order, "order");
  t_nml.assignVariable(&ewald_coefficient, "ew_coeff");
  t_nml.assignVariable(&cutoff, "cutoff");
  t_nml.assignVariable(&dsum_tol, "dsum_tol");
  t_nml.assignVariable(&mesh_ticks, "mesk_ticks");
  setStrategy(t_nml.getStringValue("accuracy"));

  // Immediately apply the new strategy.  A PMIStrategy is only valid when the object is created
  // with this constructor.  The original &pppm namelist will be referenced to see what was or
  // was not specified by the user.
  applyStrategy();
}

//-------------------------------------------------------------------------------------------------
NonbondedTheme PPPMControls::getTheme() const {
  return theme;
}

//-------------------------------------------------------------------------------------------------
int PPPMControls::getInterpolationOrder() const {
  return order;
}

//-------------------------------------------------------------------------------------------------
double PPPMControls::getEwaldCoefficient() const {
  return ewald_coefficient;
}

//-------------------------------------------------------------------------------------------------
double PPPMControls::getGaussianWidth() const {
  return 0.5 / ewald_coefficient;
}

//-------------------------------------------------------------------------------------------------
double PPPMControls::getCutoff() const {
  return cutoff;
}

//-------------------------------------------------------------------------------------------------
double PPPMControls::getDirectSumTolerance() const {
  return dsum_tol;
}

//-------------------------------------------------------------------------------------------------
int PPPMControls::getMeshSubdivisions() const {
  return mesh_ticks;
}

//-------------------------------------------------------------------------------------------------
PMIStrategy PPPMControls::getStrategy() const {
  return strat;
}

//-------------------------------------------------------------------------------------------------
void PPPMControls::setTheme(const NonbondedTheme theme_in) {
  theme = theme_in;
}

//-------------------------------------------------------------------------------------------------
void PPPMControls::setTheme(const std::string &theme_in) {
  try {
    theme = translateNonbondedTheme(theme_in);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Invalid non-bonded potential type " + theme_in + ".", "PPPMControls", "setTheme");
    case ExceptionResponse::WARN:
      rtWarn("Invalid non-bonded potential type " + theme_in + ".  The current value of " +
             getEnumerationName(theme) + " will remain in effect.", "PPPMControls", "setTheme");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void PPPMControls::setInterpolationOrder(const int order_in) {
  if (validateOrder(order_in)) {
    order = order_in;
  }
}

//-------------------------------------------------------------------------------------------------
void PPPMControls::setEwaldCoefficient(const double ewald_coefficient_in) {
  if (validateEwaldCoefficient(ewald_coefficient_in)) {
    ewald_coefficient = ewald_coefficient_in;
  }
}

//-------------------------------------------------------------------------------------------------
void PPPMControls::setCutoff(const double cutoff_in) {
  if (validateCutoff(cutoff_in)) {
    cutoff = cutoff_in;
  }
}

//-------------------------------------------------------------------------------------------------
void PPPMControls::setDirectSumTolerance(const double dsum_tol_in) {
  if (validateDirectSumTolerance(dsum_tol_in)) {
    dsum_tol = dsum_tol_in;
  }
}

//-------------------------------------------------------------------------------------------------
void PPPMControls::setMeshSubdivisions(const int mesh_ticks_in) {
  if (validateMeshSubdivisions(mesh_ticks_in)) {
    mesh_ticks = mesh_ticks_in;
  }
}

//-------------------------------------------------------------------------------------------------
void PPPMControls::setStrategy(const PMIStrategy strat_in) {
  strat = strat_in;
}

//-------------------------------------------------------------------------------------------------
void PPPMControls::setStrategy(const std::string &strat_in) {
  try {
    strat = translatePMIStrategy(strat_in);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Invalid particle-mesh accuracy setting " + strat_in + ".", "PPPMControls",
            "setStrategy");
    case ExceptionResponse::WARN:
      rtWarn("Invalid particle-mesh accuracy setting " + strat_in + ".  The current setting of " +
             getEnumerationName(strat) + " will be retained.", "PPPMControls", "setStrategy");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void PPPMControls::applyStrategy() {
  bool mod_mesh_ticks = (nml_transcript.getKeywordStatus("mesh_ticks") == InputStatus::DEFAULT);
  bool mod_order      = (nml_transcript.getKeywordStatus("order") == InputStatus::DEFAULT);
  bool mod_cutoff     = (nml_transcript.getKeywordStatus("cutoff") == InputStatus::DEFAULT);
  
  bool cutoff_problem = false;
  switch (strat) {
  case PMIStrategy::RECOMMENDED:
    if (mod_mesh_ticks) {
      if (mod_order) {
        mesh_ticks = 4;
        order = 5;
      }
      else {
        if (order == 4) {
          mesh_ticks = 5;
        }
        else if (order == 5) {
          mesh_ticks = 4;
        }
        else if (order == 6) {
          mesh_ticks = 3;
        }
        else {
          rtErr("An order of " + std::to_string(order) + " cannot be reconciled with " +
                getEnumerationName(strat) + ".", "PPPMControls", "applyStrategy");
        }
      }
    }
    else {
      if (mod_order) {
        if (mesh_ticks == 3) {
          order = 6;
        }
        else if (mesh_ticks == 4) {
          order = 5;
        }
        else if (mesh_ticks == 5) {
          order = 4;
        }
        else {
          rtErr("A spatial cell subdivision into " + std::to_string(mesh_ticks) + " grid "
                "points cannot be reconciled with " + getEnumerationName(strat) + ".",
                "PPPMControls", "applyStrategy");
        }
      }
      else {
        rtErr("Either the number of spatial cell subdivisions or the interpolation order must "
              "be left unset by the user in order to implement " + getEnumerationName(strat) +
              ".", "PPPMControls", "applyStrategy");
      }
      ewald_coefficient = ewaldCoefficient(cutoff, 1.0e-5);
    }
    break;
  case PMIStrategy::TIGHT:
    if (mod_mesh_ticks) {
      if (mod_order) {
        mesh_ticks = 5;
        order = 5;
      }
      else {
        if (order == 4) {
          mesh_ticks = 6;
        }
        else if (order == 5) {
          mesh_ticks = 5;
        }
        else if (order == 6) {
          mesh_ticks = 4;
        }
        else {
          rtErr("An order of " + std::to_string(order) + " cannot be reconciled with " +
                getEnumerationName(strat) + ".", "PPPMControls", "applyStrategy");
        }
      }
    }
    else {
      if (mod_order) {
        if (mesh_ticks == 4) {
          order = 6;
        }
        else if (mesh_ticks == 5) {
          order = 5;
        }
        else if (mesh_ticks == 6) {
          order = 4;
        }
        else {
          rtErr("A spatial cell subdivision into " + std::to_string(mesh_ticks) + " grid "
                "points cannot be reconciled with " + getEnumerationName(strat) + ".",
                "PPPMControls", "applyStrategy");
        }
      }
      else {
        rtErr("Either the number of spatial cell subdivisions or the interpolation order must "
              "be left unset by the user in order to implement " + getEnumerationName(strat) +
              ".", "PPPMControls", "applyStrategy");
      }
    }
    ewald_coefficient = ewaldCoefficient(cutoff, 1.0e-6);
    break;
  case PMIStrategy::VERY_TIGHT:
    if (mod_mesh_ticks) {
      if (mod_order) {
        mesh_ticks = 5;
        order = 5;
      }
      else {
        if (order == 5) {
          mesh_ticks = 6;
        }
        else if (order == 6) {
          mesh_ticks = 5;
        }
        else {
          rtErr("An order of " + std::to_string(order) + " cannot be reconciled with " +
                getEnumerationName(strat) + ".", "PPPMControls", "applyStrategy");
        }
      }
    }
    else {
      if (mod_order) {
        if (mesh_ticks == 5) {
          order = 6;
        }
        else if (mesh_ticks == 6) {
          order = 5;
        }
        else {
          rtErr("A spatial cell subdivision into " + std::to_string(mesh_ticks) + " grid "
                "points cannot be reconciled with " + getEnumerationName(strat) + ".",
                "PPPMControls", "applyStrategy");
        }
      }
      else {
        rtErr("Either the number of spatial cell subdivisions or the interpolation order must "
              "be left unset by the user in order to implement " + getEnumerationName(strat) +
              ".", "PPPMControls", "applyStrategy");
      }
      ewald_coefficient = ewaldCoefficient(cutoff, 1.0e-7);
    }
    break;
  case PMIStrategy::RECOMMENDED_PM_HEAVY:
    if (mod_cutoff) {
      cutoff = 0.8 * default_pme_cutoff;
    }
    else {
      cutoff_problem = true;
    }
    break;
  case PMIStrategy::RECOMMENDED_PP_HEAVY:
    if (mod_cutoff) {
      cutoff = 1.25 * default_pme_cutoff;
    }
    else {
      cutoff_problem = true;
    }
    break;
  case PMIStrategy::TIGHT_PM_HEAVY:
    if (mod_cutoff) {
      cutoff = 0.8 * default_pme_cutoff;
    }
    else {
      cutoff_problem = true;
    }
    break;
  case PMIStrategy::TIGHT_PP_HEAVY:
    if (mod_cutoff) {
      cutoff = 1.25 * default_pme_cutoff;
    }
    else {
      cutoff_problem = true;
    }
    break;
  case PMIStrategy::NO_AUTOMATION:
    break;
  }
  if (cutoff_problem) {
    rtErr("The cutoff must be modifiable in order to apply " + getEnumerationName(strat) + ".",
          "PPPMControls", "applyStrategy");
  }

  // Apply the recommended, or tight criteria, if appropriate, for a modified cutoff  
  switch (strat) {
  case PMIStrategy::RECOMMENDED:
  case PMIStrategy::TIGHT:
  case PMIStrategy::VERY_TIGHT:
  case PMIStrategy::NO_AUTOMATION:
    break;
  case PMIStrategy::RECOMMENDED_PM_HEAVY:
  case PMIStrategy::RECOMMENDED_PP_HEAVY:
    {
      const PMIStrategy save_strat = strat;
      strat = PMIStrategy::RECOMMENDED;
      applyStrategy();
      strat = save_strat;
    }
    break;
  case PMIStrategy::TIGHT_PM_HEAVY:
  case PMIStrategy::TIGHT_PP_HEAVY:
    {
      const PMIStrategy save_strat = strat;
      strat = PMIStrategy::TIGHT;
      applyStrategy();
      strat = save_strat;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
bool PPPMControls::validateOrder(const int order_in) const {
  bool ok = true;
  if (order_in < 4 || order_in > 8) {
    switch (policy) {
    case ExceptionResponse::DIE:
      if (order_in < 4) {
        rtErr("An interpolation order of " + std::to_string(order_in) + " is too low, and will "
              "make it too difficult to compute reasonable energies and forces.", "PPPMControls",
              "validateOrder");
      }
      else {
        rtErr("An interpolation order of " + std::to_string(order_in) + " is excessively high, "
              "impractical and not covered by the available interpolation methods.",
              "PPPMControls", "validateOrder");
      }
      break;
    case ExceptionResponse::WARN:
      if (order_in < 4) {
        rtWarn("An interpolation order of " + std::to_string(order_in) + " is too low, and will "
               "make it too difficult to compute reasonable energies and forces.  The current "
               "value of " + std::to_string(order_in) + " will be retained.", "PPPMControls",
               "validateOrder");
      }
      else {
        rtWarn("An interpolation order of " + std::to_string(order_in) + " is excessively high, "
               "impractical and not covered by the available interpolation methods.  The current "
               "value of " + std::to_string(order_in) + " will be retained.", "PPPMControls",
               "validateOrder");
      }
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    ok = false;
  }
  return ok;
}

//-------------------------------------------------------------------------------------------------
bool PPPMControls::validateEwaldCoefficient(const double ewald_coefficient_in) const {
  bool ok = true;
  if (ewald_coefficient_in < minimum_ewald_coefficient ||
      ewald_coefficient_in >= maximum_ewald_coefficient) {
    const std::string descriptor = (ewald_coefficient_in < minimum_ewald_coefficient) ?
                                   "wide" : "narrow";
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Particle density spread by a Gaussian function with RMS width " +
            realToString(0.5 / ewald_coefficient_in, 7, 4, NumberFormat::STANDARD_REAL) +
            " Angstroms is too " + descriptor + " for a reasonable particle-mesh interaction.",
            "PPPMControls", "validateEwaldCoefficient");
    case ExceptionResponse::WARN:
      rtWarn("Particle density spread by a Gaussian function with RMS width " +
             realToString(0.5 / ewald_coefficient_in, 7, 4, NumberFormat::STANDARD_REAL) +
             " Angstroms is too " + descriptor + " for a reasonable particle-mesh interaction.  "
             "The current value of " +
             realToString(0.5 / ewald_coefficient, 7, 4, NumberFormat::STANDARD_REAL) + " will be "
             "retained.", "PPPMControls", "validateEwaldCoefficient");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    ok = false;
  }
  return ok;
}

//-------------------------------------------------------------------------------------------------
bool PPPMControls::validateCutoff(const double cutoff_in) const {
  bool ok = true;
  if (cutoff_in >= max_pp_cutoff) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A particle-particle interaction cutoff of " +
            realToString(cutoff_in, 7, 4, NumberFormat::STANDARD_REAL) + " is too large to be "
            "practical.", "PPPMControls", "validateCutoff");
    case ExceptionResponse::WARN:
      rtWarn("A particle-particle interaction cutoff of " +
             realToString(cutoff_in, 7, 4, NumberFormat::STANDARD_REAL) + " is too large to be "
             "practical.  The current value of " +
             realToString(cutoff_in, 7, 4, NumberFormat::STANDARD_REAL) + " will be retained.",
             "PPPMControls", "validateCutoff");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    ok = false;
  }
  return ok;
}

//-------------------------------------------------------------------------------------------------
bool PPPMControls::validateDirectSumTolerance(const double tol_in) const {
  bool ok = true;
  if (tol_in >= 1.0e-4) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("To discard a proportion of " +
            realToString(tol_in, 7, 4, NumberFormat::STANDARD_REAL) + " of the direct, "
            "particle-particle interactions would lead to forces that are too imprecise.",
            "PPPMControls", "validateCutoff");
    case ExceptionResponse::WARN:
      rtWarn("To discard a proportion of " +
             realToString(tol_in, 7, 4, NumberFormat::STANDARD_REAL) + " of the direct, "
             "particle-particle interactions would lead to forces that are too imprecise.  "
             "The current value of " + realToString(tol_in, 7, 4, NumberFormat::STANDARD_REAL) +
             " will be retained.", "PPPMControls", "validateCutoff");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    ok = false;
  }
  if (ewaldCoefficient(max_pp_cutoff, tol_in) < minimum_ewald_coefficient) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("To discard a proportion of " + 
            realToString(tol_in, 7, 4, NumberFormat::STANDARD_REAL) + " of the direct, "
            "particle-particle interactions is impractical, even with the maximum "
            "particle-particle interaction cutoff.", "PPPMControls", "validateCutoff");
    case ExceptionResponse::WARN:
      rtWarn("To discard a proportion of " + 
             realToString(tol_in, 7, 4, NumberFormat::STANDARD_REAL) + " of the direct, "
             "particle-particle interactions is impractical, even with the maximum "
             "particle-particle interaction cutoff.  The current value of " +
             realToString(tol_in, 7, 4, NumberFormat::STANDARD_REAL) + " will be retained.",
             "PPPMControls", "validateCutoff");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    ok = false;      
  }
  return ok;
}

//-------------------------------------------------------------------------------------------------
bool PPPMControls::validateMeshSubdivisions(const int mesh_ticks_in) const {
  bool ok = true;
  if (mesh_ticks_in < 2) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("At least two grid elements of the the particle-mesh interaction grid must lie "
            "within each spatial decomposition cell.  It is not possible to get realistic forces, "
            "even with the maximum interpolation order, otherwise.", "PPPMControls",
            "validateMeshTicks");
    case ExceptionResponse::WARN:
      rtWarn("At least two grid elements of the the particle-mesh interaction grid must lie "
             "within each spatial decomposition cell.  It is not possible to get realistic "
             "forces, even with the maximum interpolation order, otherwise.  The current value "
             "of " + std::to_string(mesh_ticks) + " will be maintained.", "PPPMControls",
            "validateMeshTicks");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    ok = false;
  }
  return ok;
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator pppmInput(const TextFile &tf, int *start_line, bool *found,
                           const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("pppm", CaseSensitivity::AUTOMATIC, policy, "Namelist containing "
                         "parameters for the particle-particle, particle-mesh interactions in a "
                         "simulation with periodic boundary conditions.");
  t_nml.addKeyword("theme", NamelistType::STRING, getEnumerationName(default_pppm_theme));
  t_nml.addHelp("theme", "The non-bonded potential to be treated by this particle-particle, "
                "particle-mesh decomposition.  Default " + getEnumerationName(default_pppm_theme) +
                ".");
  t_nml.addKeyword("order", NamelistType::INTEGER, std::to_string(default_charge_mapping_order));
  t_nml.addHelp("order", "Order of interpolation for particle density mapping onto the "
                "particle-mesh interaction grid.  The default of " +
                std::to_string(default_charge_mapping_order) + " pairs well with the default "
                "number of mesh subdivisions (" + std::to_string(default_mesh_ticks) + ".");
  t_nml.addKeyword("ew_coeff", NamelistType::REAL);
  t_nml.addHelp("ew_coeff", "The 'Ewald coefficient' used in splitting the potential.  This is "
                "'a' in the equation erf(a * r) / (r^n) for the particle-mesh interaction "
                "component of a potential of the form 1 / (r^n).  Units of this quantity are "
                "inverse Angstroms.");
  t_nml.addKeyword("gaussian", NamelistType::REAL);
  t_nml.addHelp("gss_sigma", "The RMS width of the Gaussian used in splitting the potential, in "
                "units of Angstroms.  This is related to the 'Ewald Coefficient' ew_coeff by the "
                "expression ew_coeff = 0.5 / gss_sigma, and splits a potential of the form 1 / "
                "(r^n) into a particle-mesh interaction erf(0.5 * r / gss_sigma) / (r^n) and a "
                "particle-particle interaction erfc(0.5 * r / gss_sigma) / (r^n).");
  t_nml.addKeyword("cutoff", NamelistType::REAL, std::to_string(default_pme_cutoff));
  t_nml.addHelp("cutoff", "The particle-particle interaction cutoff, in units of Angstroms.  At "
                "the edge of this cutoff, the interaction between two point particles and two "
                "spherical Gaussian density distributions will differ by proportion indicated by "
                "the direct sum tolerance (dsum_tol).");
  t_nml.addKeyword("dsum_tol", NamelistType::REAL, std::to_string(default_dsum_tol));
  t_nml.addHelp("dsum_tol", "The proportion by which the interaction of two point particles and "
                "the interaction of two spherical Gaussian distributions (e.g. of charge, of "
                "dispersion sources) will differ at the particle-particle interaction cutoff.  "
                "This indicates the proportion of the 'direct' particle-particle interaction "
                "which will be discarded.");
  t_nml.addKeyword("mesh_ticks", NamelistType::INTEGER, std::to_string(default_mesh_ticks));
  t_nml.addHelp("mesh_ticks", "The number of particle-mesh interaction grid spacings that will "
                "span each spatial decomposition cell in all directions.  Spatial decomposition "
                "cells are parallelepiped volumes, at least half the particle-particle "
                "interaction cutoff in width (see the 'cutoff' keyword in this namelist).  Nested "
                "loops over all particles in one spatial decomposition cell and any other cells "
                "within +/-2 of its position on the cell grid will yield a complete accounting of "
                "all particle-particle interactions.  This coupling of the particle-mesh "
                "interaction grid spacing to another grid sized by the particle-particle "
                "interaction cutoff is convenient for STORMM's algorithms and also enforces a "
                "natural relationship between the Gaussian spread of particle density and the "
                "particle-mesh interaction grid spacing needed to maintain a particular level of "
                "accuracy (balancing errors in particle-particle and particle-mesh interactions) "
                "for a given interpolation order.  Combinations of mesh_ticks = 4 and order = 5 "
                "will pair well with dsum_tol = 5.0e-6 (the STORMM defaults), yielding "
                "non-bonded calculations with accuracies intermediate between Amber and CHARMM "
                "default parameters.  The cutoff itself does not have much significance in P3M--"
                "the point is that interactions be infinite--and in this scheme similar "
                "accuracies will be obtained for a variety of cutoffs, e.g. 7.0A, 10.0A.  "
                "Modifying the cutoff parameter thus alters the number of spatial decomposition "
                "cells and thus the density of the FFT grid, moving effort between "
                "particle-particle and particle-mesh interactions.  The cutoff can then be chosen "
                "to optimize the time cost of the simulation without altering the accuracy of the "
                "forces guiding dynamics.");
  t_nml.addKeyword("accuracy", NamelistType::STRING, getEnumerationName(default_pppm_accuracy));
  t_nml.addHelp("accuracy", "Set the accuracy level for nonbonded calculations.  This can be used "
                "to fill in additional settings, based on incomplete user input for mesh_ticks, "
                "order, and dsum_tol.  A suitable choice of cutoff can also be indicated which "
                "will shift work between particle-particle pair interactions and FFT "
                "computations.  Valid inputs include 'RECOMMENDED / STANDARD / DEFAULT' (these "
                "will alter any unset parameters so as to not change the overall accuracy which "
                "would be produced be default settings), 'TIGHT / HIGH / FINE' (these will alter "
                "any unset parameters to produce accuracies roughly one order of magnitude better "
                "than the default settings), 'VERYTIGHT / VERY_TIGHT / VERYFINE / VERY_FINE / "
                "VERYHIGH / VERY_HIGH' (these will produce accuracies roughly to orders of "
                "magnitude better than the default settings), 'GRIDHEAVY / GRID_HEAVY / PMHEAVY / "
                "PM_HEAVY' (these will produce accuracies commensurate with the default settings "
                "but with a shorter cutoff, more FFT work), 'PAIRHEAVY / PAIR_HEAVY / PPHEAVY / "
                "PP_HEAVY' (these will produce default accuracies with a longer cutoff, more "
                "pairwise work), and also 'TIGHT_(...)' with extensions for pairwise- or "
                "FFT-heavy protocols.  The default value of 'NONE / NO_AUTOMATION' will not "
                "alter any settings aside from user inputs or other existing defaults.");
  return t_nml;
}

} // namespace namelist
} // namespace stormm
