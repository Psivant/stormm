#include "copyright.h"
#include "Constants/scaling.h"
#include "Numerics/split_fixed_precision.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "namelist_element.h"
#include "nml_precision.h"

namespace stormm {
namespace namelist {

using constants::CaseSensitivity;
using constants::translatePrecisionModel;
using numerics::checkGlobalPositionBits;
using numerics::checkLocalPositionBits;
using numerics::checkVelocityBits;
using numerics::checkForceBits;
using numerics::checkEnergyBits;
using numerics::checkChargeMeshBits;
using parse::NumberFormat;
using parse::realToString;
  
//-------------------------------------------------------------------------------------------------
PrecisionControls::PrecisionControls(const ExceptionResponse policy_in,
                                     const WrapTextSearch wrap) :
    policy{policy_in},
    globalpos_scale_bits{default_globalpos_scale_bits},
    localpos_scale_bits{default_localpos_scale_bits},
    velocity_scale_bits{default_velocity_scale_bits},
    force_scale_bits{default_force_scale_bits},
    energy_scale_bits{default_energy_scale_bits},
    charge_mesh_scale_bits{default_charge_mesh_scale_bits},
    bond_constraint_tol{default_precision_constraint_tol},
    valence_method{translatePrecisionModel(std::string(default_precision_valence_method))},
    nonbonded_method{translatePrecisionModel(std::string(default_precision_nonbonded_method))},
    pme_method{translatePrecisionModel(std::string(default_precision_pme_method))},
    nml_transcript{"precision"}
{}

//-------------------------------------------------------------------------------------------------
PrecisionControls::PrecisionControls(const TextFile &tf, int *start_line, bool *found_nml,
                                     const ExceptionResponse policy_in,
                                     const WrapTextSearch wrap) :
    PrecisionControls(policy_in)
{
  NamelistEmulator t_nml = precisionInput(tf, start_line, found_nml, policy, wrap);
  nml_transcript = t_nml;
  globalpos_scale_bits = t_nml.getIntValue("globalpos_bits");
  localpos_scale_bits = t_nml.getIntValue("localpos_bits");
  velocity_scale_bits = t_nml.getIntValue("velocity_bits");
  force_scale_bits = t_nml.getIntValue("force_bits");
  energy_scale_bits = t_nml.getIntValue("energy_bits");
  charge_mesh_scale_bits = t_nml.getIntValue("charge_mesh_bits");
  bond_constraint_tol = t_nml.getRealValue("bond_constraint_tol");
  const std::string valence_method_str   = t_nml.getStringValue("valence");
  const std::string nonbonded_method_str = t_nml.getStringValue("nonbonded");
  const std::string pme_method_str = t_nml.getStringValue("pmegrid");
  valence_method = translatePrecisionModel(valence_method_str);
  nonbonded_method = translatePrecisionModel(nonbonded_method_str);
  pme_method = translatePrecisionModel(pme_method_str);

  // Validate input
  checkGlobalPositionBits(globalpos_scale_bits);
  checkLocalPositionBits(localpos_scale_bits);
  checkVelocityBits(velocity_scale_bits);
  checkForceBits(force_scale_bits);
  checkEnergyBits(energy_scale_bits);
  checkChargeMeshBits(charge_mesh_scale_bits, pme_method);
  validateBondConstraintTol();
}

//-------------------------------------------------------------------------------------------------
int PrecisionControls::getGlobalPosScalingBits() const {
  return globalpos_scale_bits;
}

//-------------------------------------------------------------------------------------------------
int PrecisionControls::getLocalPosScalingBits() const {
  return localpos_scale_bits;
}

//-------------------------------------------------------------------------------------------------
int PrecisionControls::getVelocityScalingBits() const {
  return velocity_scale_bits;
}

//-------------------------------------------------------------------------------------------------
int PrecisionControls::getForceScalingBits() const {
  return force_scale_bits;
}

//-------------------------------------------------------------------------------------------------
int PrecisionControls::getEnergyScalingBits() const {
  return energy_scale_bits;
}

//-------------------------------------------------------------------------------------------------
int PrecisionControls::getChargeMeshScalingBits() const {
  return charge_mesh_scale_bits;
}

//-------------------------------------------------------------------------------------------------
double PrecisionControls::getBondConstraintTolerance() const {
  return bond_constraint_tol;
}

//-------------------------------------------------------------------------------------------------
PrecisionModel PrecisionControls::getValenceMethod() const {
  return valence_method;
}

//-------------------------------------------------------------------------------------------------
PrecisionModel PrecisionControls::getNonbondedMethod() const {
  return nonbonded_method;
}

//-------------------------------------------------------------------------------------------------
PrecisionModel PrecisionControls::getParticleMeshEwaldMethod() const {
  return pme_method;
}

//-------------------------------------------------------------------------------------------------
const NamelistEmulator& PrecisionControls::getTranscript() const {
  return nml_transcript;
}

//-------------------------------------------------------------------------------------------------
void PrecisionControls::setGlobalPosScalingBits(const int bitval) {
  globalpos_scale_bits = bitval;
  checkGlobalPositionBits(globalpos_scale_bits);
}

//-------------------------------------------------------------------------------------------------
void PrecisionControls::setLocalPosScalingBits(const int bitval) {
  localpos_scale_bits = bitval;
  checkLocalPositionBits(localpos_scale_bits);
}

//-------------------------------------------------------------------------------------------------
void PrecisionControls::setVelocityScalingBits(const int bitval) {
  velocity_scale_bits = bitval;
  checkVelocityBits(velocity_scale_bits);
}

//-------------------------------------------------------------------------------------------------
void PrecisionControls::setForceScalingBits(const int bitval) {
  force_scale_bits = bitval;
  checkForceBits(force_scale_bits);
}

//-------------------------------------------------------------------------------------------------
void PrecisionControls::setEnergyScalingBits(const int bitval) {
  energy_scale_bits = bitval;
  checkEnergyBits(energy_scale_bits);
}

//-------------------------------------------------------------------------------------------------
void PrecisionControls::setChargeMeshScalingBits(const int bitval) {
  charge_mesh_scale_bits = bitval;
  checkChargeMeshBits(charge_mesh_scale_bits, pme_method);
}

//-------------------------------------------------------------------------------------------------
void PrecisionControls::setBondConstraintTolerance(const double tol) {
  bond_constraint_tol = tol;
  validateBondConstraintTol();
}

//-------------------------------------------------------------------------------------------------
void PrecisionControls::setValenceMethod(const PrecisionModel pmodel) {
  valence_method = pmodel;
}

//-------------------------------------------------------------------------------------------------
void PrecisionControls::setNonbondedMethod(const PrecisionModel pmodel) {
  nonbonded_method = pmodel;
}

//-------------------------------------------------------------------------------------------------
void PrecisionControls::setParticleMeshEwaldMethod(const PrecisionModel pmodel) {
  pme_method = pmodel;
}

//-------------------------------------------------------------------------------------------------
void PrecisionControls::validateBondConstraintTol() {
  if (bond_constraint_tol < min_precision_constraint_tol) {
    rtErr("The tightest tolerance on bond constraints is set to " +
          realToString(min_precision_constraint_tol, 11, 4, NumberFormat::SCIENTIFIC) +
          " to ensure stability in the convergence.  A value of " +
          realToString(bond_constraint_tol, 11, 4, NumberFormat::SCIENTIFIC) + " is too small."
          "PrecisionControls", "validateBondConstraintTol");
  }
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator precisionInput(const TextFile &tf, int *start_line, bool *found,
                                const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("precision", CaseSensitivity::AUTOMATIC, policy, "Wraps directives for "
                         "tuning the precision model and accumulation of molecular mechanics "
                         "calculations.");
  t_nml.addKeyword(NamelistElement("globalpos_bits", NamelistType::INTEGER,
                                   std::to_string(default_globalpos_scale_bits)));
  t_nml.addKeyword(NamelistElement("localpos_bits", NamelistType::INTEGER,
                                   std::to_string(default_localpos_scale_bits)));
  t_nml.addKeyword(NamelistElement("velocity_bits", NamelistType::INTEGER,
                                   std::to_string(default_velocity_scale_bits)));
  t_nml.addKeyword(NamelistElement("force_bits", NamelistType::INTEGER,
                                   std::to_string(default_force_scale_bits)));
  t_nml.addKeyword(NamelistElement("energy_bits", NamelistType::INTEGER,
                                   std::to_string(default_energy_scale_bits)));
  t_nml.addKeyword(NamelistElement("charge_mesh_bits", NamelistType::INTEGER,
                                   std::to_string(default_charge_mesh_scale_bits)));
  t_nml.addKeyword(NamelistElement("bond_constraint_tol", NamelistType::REAL,
                                   realToString(default_precision_constraint_tol, 11, 4,
                                                NumberFormat::SCIENTIFIC)));
  t_nml.addKeyword(NamelistElement("valence", NamelistType::STRING,
                                   std::string(default_precision_valence_method)));
  t_nml.addKeyword(NamelistElement("nonbonded", NamelistType::STRING,
                                   std::string(default_precision_nonbonded_method)));
  t_nml.addKeyword(NamelistElement("pmegrid", NamelistType::STRING,
                                   std::string(default_precision_pme_method)));
  t_nml.addHelp("globalpos_bits", "Number of bits after the decimal in fixed-precision "
                "representations of global particle coordinates.");
  t_nml.addHelp("localpos_bits", "Number of bits after the decimal in fixed-precision "
                "representations of local particle coordinates.");
  t_nml.addHelp("velocity_bits", "Number of bits after the decimal in fixed-precision "
                "representations of particle velocities.");
  t_nml.addHelp("force_bits", "Number of bits after the decimal in fixed-precision "
                "force accumulation.");
  t_nml.addHelp("energy_bits", "Number of bits after the decimal in fixed-precision "
                "energy accumulation.");
  t_nml.addHelp("charge_mesh_bits", "Number of bits after the decimal in fixed-precision "
                "charge density accumulation on the PME grid.  Charges are accumulated in atomic "
                "units, with 32 bit signed integer accumulators if the precision model is SINGLE "
                "or SINGLE_PLUS, or 64 bit signed integer accumulators if the precision model is "
                "DOUBLE.");
  t_nml.addHelp("bond_constraint_tol", "Tolerance to which constrained bonds must meet their "
                "equilibrium lengths.  The criterion is for the deviation of the squared actual "
                "length, in Angstroms, to differ from the squared target length of the bond by "
                "less than this value.");
  t_nml.addHelp("valence", "Precision model to use in valence calculations.  Choices include "
                "'single', 'single_plus', and 'double.'");
  t_nml.addHelp("nonbonded", "Precision model to use in non-bonded short-ranged calculations.  "
                "Choices include 'single', 'single_plus', and 'double.'");
  t_nml.addHelp("pmegrid", "Precision model to use in mesh-based Particle-Mesh Ewald "
                "calculations.  Choices include 'single', 'single_plus', and 'double.'");
  
  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  return t_nml;
}

} // namespace namelist
} // namespace stormm
