#include "copyright.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "Topology/atomgraph_enumerators.h"
#include "namelist_element.h"
#include "nml_solvent.h"

namespace stormm {
namespace namelist {

using constants::CaseSensitivity;
using parse::strncmpCased;
using topology::translateAtomicRadiusSet;
using topology::translateImplicitSolventModel;
  
//-------------------------------------------------------------------------------------------------
SolventControls::SolventControls(const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    policy{policy_in},
    gb_style{translateImplicitSolventModel(default_solvent_igb)},
    born_radii_cutoff{default_solvent_rgbmax},
    internal_dielectric{default_solvent_intdiel},
    external_dielectric{default_solvent_extdiel},
    salt_concentration{default_solvent_saltcon},
    pb_radii{translateAtomicRadiusSet(std::string(default_solvent_pbradii))},
    nml_transcript{"solvent"}
{}

//-------------------------------------------------------------------------------------------------
SolventControls::SolventControls(const TextFile &tf, int *start_line, bool *found_nml,
                                 const ExceptionResponse policy_in, const WrapTextSearch wrap) :
  SolventControls(policy_in)
{
  NamelistEmulator t_nml = solventInput(tf, start_line, found_nml, policy, wrap);
  nml_transcript = t_nml;
  gb_style            = translateImplicitSolventModel(t_nml.getIntValue("igb"));
  born_radii_cutoff   = t_nml.getRealValue("rgbmax");
  internal_dielectric = t_nml.getRealValue("intdiel");
  external_dielectric = t_nml.getRealValue("extdiel");
  salt_concentration  = t_nml.getRealValue("saltcon");

  // If the radius set is unspecified, there are Generalized Born solvent models that demand a
  // particular set of radii, and others that will do best with a default set.  Make the choice for
  // the user if there is no specification.
  if (t_nml.getKeywordStatus("pbradii") == InputStatus::MISSING) {
    switch (gb_style) {
    case ImplicitSolventModel::NONE:
      pb_radii = AtomicRadiusSet::NONE;
      break;
    case ImplicitSolventModel::HCT_GB:
      pb_radii = AtomicRadiusSet::MBONDI;
      break;
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
      pb_radii = AtomicRadiusSet::MBONDI2;
      break;
    case ImplicitSolventModel::NECK_GB:
      pb_radii = AtomicRadiusSet::BONDI;
      break;
    case ImplicitSolventModel::NECK_GB_II:
      pb_radii = AtomicRadiusSet::MBONDI3;
      break;
    }
  }
  else {
    pb_radii = translateAtomicRadiusSet(t_nml.getStringValue("pbradii"));
  }

  // Checks on the input
  validateBornRadiiCutoff();
  validateInternalDielectric();
  validateExternalDielectric();
  validateSaltConcentration();
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventModel SolventControls::getImplicitSolventModel() const {
  return gb_style;
}

//-------------------------------------------------------------------------------------------------
double SolventControls::getBornRadiiCutoff() const {
  return born_radii_cutoff;
}

//-------------------------------------------------------------------------------------------------
double SolventControls::getInternalDielectric() const {
  return internal_dielectric;
}

//-------------------------------------------------------------------------------------------------
double SolventControls::getExternalDielectric() const {
  return external_dielectric;
}

//-------------------------------------------------------------------------------------------------
double SolventControls::getSaltConcentration() const {
  return salt_concentration;
}

//-------------------------------------------------------------------------------------------------
AtomicRadiusSet SolventControls::getPBRadiiSet() const {
  return pb_radii;
}

//-------------------------------------------------------------------------------------------------
const NamelistEmulator& SolventControls::getTranscript() const {
  return nml_transcript;
}

//-------------------------------------------------------------------------------------------------
void SolventControls::setImplicitSolventModel(const int ism_in) {
  gb_style = translateImplicitSolventModel(ism_in, policy);
}

//-------------------------------------------------------------------------------------------------
void SolventControls::setImplicitSolventModel(const ImplicitSolventModel ism_in) {
  gb_style = ism_in;
}

//-------------------------------------------------------------------------------------------------
void SolventControls::setBornRadiiCutoff(const double rgbmax_in) {
  born_radii_cutoff = rgbmax_in;
  validateBornRadiiCutoff();
}

//-------------------------------------------------------------------------------------------------
void SolventControls::setInternalDielectric(const double idiel_in) {
  internal_dielectric = idiel_in;
  validateInternalDielectric();
}

//-------------------------------------------------------------------------------------------------
void SolventControls::setExternalDielectric(const double ediel_in) {
  external_dielectric = ediel_in;
  validateExternalDielectric();
}

//-------------------------------------------------------------------------------------------------
void SolventControls::setSaltConcentration(double saltcon_in) {
  salt_concentration = saltcon_in;
  validateSaltConcentration();
}

//-------------------------------------------------------------------------------------------------
void SolventControls::choosePBRadiiSet(const std::string &pbrad_in) {
  pb_radii = translateAtomicRadiusSet(pbrad_in, policy);
}

//-------------------------------------------------------------------------------------------------
void SolventControls::choosePBRadiiSet(const AtomicRadiusSet pbrad_in) {
  pb_radii = pbrad_in;
}

//-------------------------------------------------------------------------------------------------
void SolventControls::validateBornRadiiCutoff() {
  if (born_radii_cutoff < 0.0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Generalized Born radii computation cannot proceed with a cutoff of " +
            std::to_string(born_radii_cutoff) + ".", "SolventControls", "validateBornRadii");
    case ExceptionResponse::WARN:
      rtWarn("Generalized Born radii computation cannot proceed with a cutoff of " +
             std::to_string(born_radii_cutoff) + ".  The default value of " +
             std::to_string(default_solvent_rgbmax) + " (no cutoff on Born radius contributions) "
             "will be restored.", "SolventControls", "validateBornRadii");
      born_radii_cutoff = default_solvent_rgbmax;
      break;
    case ExceptionResponse::SILENT:
      born_radii_cutoff = default_solvent_rgbmax;
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void SolventControls::validateInternalDielectric() {
  if (internal_dielectric < 0.01) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The dielectric constant cannot be so small as to amplify electrostatics to huge "
            "proportions.  A value of " + std::to_string(internal_dielectric) + " will likely "
            "break calculations.", "SolventControls", "validateInternalDielectric");
    case ExceptionResponse::WARN:
      rtWarn("The dielectric constant cannot be so small as to amplify electrostatics to huge "
             "proportions.  A value of " + std::to_string(internal_dielectric) + " will likely "
             "break calculations and will therefore be replaced by the default value of " +
             std::to_string(default_solvent_intdiel) + ".", "SolventControls",
             "validateInternalDielectric");
      internal_dielectric = default_solvent_intdiel;
      break;
    case ExceptionResponse::SILENT:
      internal_dielectric = default_solvent_intdiel;
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void SolventControls::validateExternalDielectric() {
  if (external_dielectric < 0.01) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The dielectric constant cannot be so small as to amplify electrostatics to huge "
            "proportions.  A value of " + std::to_string(external_dielectric) + " will likely "
            "break calculations.", "SolventControls", "validateExternalDielectric");
    case ExceptionResponse::WARN:
      rtWarn("The dielectric constant cannot be so small as to amplify electrostatics to huge "
             "proportions.  A value of " + std::to_string(external_dielectric) + " will likely "
             "break calculations and will therefore be replaced by the default value of " +
             std::to_string(default_solvent_extdiel) + ".", "SolventControls",
             "validateExternalDielectric");
      external_dielectric = default_solvent_extdiel;
      break;
    case ExceptionResponse::SILENT:
      external_dielectric = default_solvent_extdiel;
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void SolventControls::validateSaltConcentration() {
  if (salt_concentration < 0.0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A negative salt concentration of " + std::to_string(salt_concentration) +
            " is nonsensical.", "SolventControls", "validateSaltConcentration");
    case ExceptionResponse::WARN:
      rtErr("A negative salt concentration of " + std::to_string(salt_concentration) +
            " is nonsensical and will be replaced by the default value of " +
            std::to_string(default_solvent_saltcon) + ".", "SolventControls",
            "validateSaltConcentration");
      salt_concentration = default_solvent_saltcon;
      break;
    case ExceptionResponse::SILENT:
      salt_concentration = default_solvent_saltcon;
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator solventInput(const TextFile &tf, int *start_line, bool *found,
                              const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("solvent", CaseSensitivity::AUTOMATIC, policy, "Wraps directives needed "
                         "to define the solvent model of a molecular mechanical system.  This is "
                         "not the same as the explicit water model in use, and generally refers "
                         "to implicit solvent methods.");
  t_nml.addKeyword(NamelistElement("igb", NamelistType::INTEGER,
                                   std::to_string(default_solvent_igb)));
  t_nml.addKeyword(NamelistElement("intdiel", NamelistType::REAL,
                                   std::to_string(default_solvent_intdiel)));
  t_nml.addKeyword(NamelistElement("extdiel", NamelistType::REAL,
                                   std::to_string(default_solvent_extdiel)));
  t_nml.addKeyword(NamelistElement("saltcon", NamelistType::REAL,
                                   std::to_string(default_solvent_saltcon)));
  t_nml.addKeyword(NamelistElement("rgbmax", NamelistType::REAL,
                                   std::to_string(default_solvent_rgbmax)));
  t_nml.addKeyword(NamelistElement("pbradii", NamelistType::STRING,
                                   std::string(default_solvent_pbradii)));
  t_nml.addHelp("igb", "Definition of the Generalized Born model to use.  Numerical settings "
                "follow those found in the eponymous keyword of Amber sander's &cntrl namelist: "
                "0 or 6 (no Generalized Born model is used), 1 (Hawkins, Cramer, Truhlar), 2 or 5 "
                "(distinct models contributed by Onufriev, Bashford, and Case), 7 or 8 (distinct "
                "models contributed by Mongan, the \"neck\" GB method).");
  t_nml.addHelp("intdiel", "Internal dielectric constant.  No value other than the default of "
                "1.0 has been well tested.");
  t_nml.addHelp("extdiel", "External dielectric constant for the surrounding continuum solvent.");
  t_nml.addHelp("saltcon", "Concentration of 1-1 monovalent ions in the continuum solvent.");
  t_nml.addHelp("rgbmax", "Maximum distance at which two particles can contribute to each other's "
                "effective Born radii.  Not to be confused with the non-bonded cutoff, which is "
                "the distance past which actual interactions between particles are discarded.");
  t_nml.addHelp("pbradii", "Poisson-Boltzmann radius set (this defines the baseline Generalized "
                "Born radii as well).  Acceptable values include \"Bondi\", \"Amber6\", "
                "\"mBondi\", \"mBond2\", \"mBondi3\", or \"none\".");

  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  
  return t_nml;
}

} // namespace namelist
} // namespace stormm
