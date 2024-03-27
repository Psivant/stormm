#include <cmath>
#include "copyright.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/parse.h"
#include "layered_potential_metrics.h"

namespace stormm {
namespace energy {

using energy::minimum_ewald_coefficient;
using energy::maximum_ewald_coefficient;
using parse::realToString;
using parse::NumberFormat;

//-------------------------------------------------------------------------------------------------
LayeredPotentialMetrics::LayeredPotentialMetrics(const DecomposablePotential form_in,
                                                 const BoundaryCondition boundaries_in,
                                                 const double coulomb_in,
                                                 const VdwCombiningRule mixing_rule_in,
                                                 const std::vector<double> &exponent_c_in,
                                                 const double particle_particle_cutoff,
                                                 const double range_multiplier_in,
                                                 const double range_limit_in,
                                                 const double vdw_transition_midpoint_in,
                                                 const double vdw_transition_intensity_in,
                                                 const double ewald_coefficient_in,
                                                 const std::vector<double> &box_vectors_in) :
    form{form_in}, boundaries{boundaries_in}, layer_count{0},
    cutoffs{std::vector<double>(1, particle_particle_cutoff)},
    coulomb_constant{coulomb_in}, mixing_rule{mixing_rule_in}, exponent_c{},
    range_multiplier{range_multiplier_in}, range_limit{range_limit_in},
    vdw_transition_midpoint{vdw_transition_midpoint_in},
    vdw_transition_intensity{vdw_transition_intensity_in}, ewald_coefficient{ewald_coefficient_in},
    box_vectors{box_vectors_in}
{
  setExponentFactor(exponent_c_in);
  setCutoff(particle_particle_cutoff);
}
  
//-------------------------------------------------------------------------------------------------
DecomposablePotential LayeredPotentialMetrics::getForm() const {
  return form;
}

//-------------------------------------------------------------------------------------------------
BoundaryCondition LayeredPotentialMetrics::getBoundaryCondition() const {
  return boundaries;
}

//-------------------------------------------------------------------------------------------------
int LayeredPotentialMetrics::getLayerCount() const {
  return layer_count;
}

//-------------------------------------------------------------------------------------------------
double LayeredPotentialMetrics::getCutoff(const int layer_index) const {
  validateLayerIndex(layer_index, "getCutoff");
  return cutoffs[layer_index];
}

//-------------------------------------------------------------------------------------------------
double LayeredPotentialMetrics::getCoulombConstant() const {
  return coulomb_constant;
}

//-------------------------------------------------------------------------------------------------
VdwCombiningRule LayeredPotentialMetrics::getMixingRule() const {
  return mixing_rule;
}

//-------------------------------------------------------------------------------------------------
double LayeredPotentialMetrics::getExponentFactor(const int factor_index,
                                                  const int layer_index) const {
  if (factor_index < 0 || factor_index > 2) {
    rtErr("The smoothing function is a series of three exponential terms.  Specify an index in "
          "the range [0, 2].", "LayeredPotentialMetrics", "getExponentFactor");
  }
  validateLayerIndex(layer_index, "getExponentFactor");
  switch (factor_index) {
  case 0:
    return exponent_c[layer_index].x;
  case 1:
    return exponent_c[layer_index].y;
  case 2:
    return exponent_c[layer_index].z;
  default:
    return 0.0;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double LayeredPotentialMetrics::getRangeCompounding() const {
  return range_multiplier;
}

//-------------------------------------------------------------------------------------------------
double LayeredPotentialMetrics::getMaximumRange() const {
  return range_limit;
}

//-------------------------------------------------------------------------------------------------
double LayeredPotentialMetrics::getVdwTransitionMidpoint() const {
  return vdw_transition_midpoint;
}

//-------------------------------------------------------------------------------------------------
double LayeredPotentialMetrics::getVdwTransitionIntensity() const {
  return vdw_transition_intensity;
}

//-------------------------------------------------------------------------------------------------
double LayeredPotentialMetrics::getEwaldCoefficient() const {
  return ewald_coefficient;
}

//-------------------------------------------------------------------------------------------------
const std::vector<double>& LayeredPotentialMetrics::getUnitCellVectors() const {
  return box_vectors;
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::setForm(const DecomposablePotential form_in) {
  form = form_in;
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::setBoundaryCondition(const BoundaryCondition boundary_in) {
  boundaries = boundary_in;
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::setCutoff(const double cutoff_in) {
  if (cutoff_in < 0.0) {
    rtErr("Negative cutoffs are invalid.", "LayeredPotentialMetrics", "setCutoff");
  }
  cutoffs.resize(1);
  cutoffs[0] = cutoff_in;
  computeLayers();
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::setCoulombConstant(const double coulomb_in) {
  coulomb_constant = coulomb_in;
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::setMixingRule(const VdwCombiningRule mixing_rule_in) {
  mixing_rule = mixing_rule_in;
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::setExponentFactor(const int factor_index, const double exp_c_in) {
  if (factor_index < 0 || factor_index > 2) {
    rtErr("The smoothing function is a series of three exponential terms.  Specify an index in "
          "the range [0, 2].", "LayeredPotentialMetrics", "setExponentFactor");
  }
  exponent_c.resize(2);
  switch (factor_index) {
  case 0:
    exponent_c[1].x = exp_c_in;
    break;
  case 1:
    exponent_c[1].y = exp_c_in;
    break;
  case 2:
    exponent_c[1].z = exp_c_in;
    break;
  default:
    break;
  }
  computeLayers();
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::setExponentFactor(const std::vector<double> &exp_c_in) {
  if (exp_c_in.size() != 3) {
    rtErr("The smoothing function is a series of three exponential terms.  Provide them all in a "
          "single vector.", "LayeredPotentialMetrics", "setExponentFactor");
  }
  exponent_c.resize(2);
  exponent_c[0] = { 0.0, 0.0, 0.0 };
  exponent_c[1].x = exp_c_in[0];
  exponent_c[1].y = exp_c_in[1];
  exponent_c[1].z = exp_c_in[2];
  computeLayers();
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::setExponentFactor(const double3 exp_c_in) {
  exponent_c.resize(2);
  exponent_c[0] = { 0.0, 0.0, 0.0 };
  exponent_c[1] = exp_c_in;
  computeLayers();
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::setRangeCompounding(const double range_multiplier_in) {
  if (range_multiplier_in < hail_minimum_range_multiplier - constants::small ||
      range_multiplier_in > hail_maximum_range_multiplier + constants::small) {
    const NumberFormat std_real = NumberFormat::STANDARD_REAL;
    rtErr("The range multiplier (" + realToString(range_multiplier_in, 4, 2, std_real) + ") for "
          "successive layers must be at least " +
          realToString(hail_minimum_range_multiplier, 3, 1, std_real) + "and no greater than " +
          realToString(hail_maximum_range_multiplier, 3, 1, std_real) + ".",
          "LayeredPotentialMetrics", "setRangeCompounding");
  }
  range_multiplier = range_multiplier_in;
  computeLayers();
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::setMaximumRange(const double range_limit_in) {
  if (range_limit_in < hail_minimum_range_limit) {
    const NumberFormat std_real = NumberFormat::STANDARD_REAL;
    rtErr("The range limit (" + realToString(range_limit_in, 5, 2, std_real) + " can be no less "
          "than " + realToString(range_limit_in, 5, 2, std_real) + ".",
          "LayeredPotentialMetrics", "setMaximumRange");
  }

  // The range limit is set if isolated boundary conditions apply.  Otherwise, there is no effect.
  switch (boundaries) {
  case BoundaryCondition::ISOLATED:
    range_limit = range_limit_in;
    break;
  case BoundaryCondition::PERIODIC:
    break;
  }
  computeLayers();
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::setVdwTransitionMidpoint(const double midpoint_in) {
  vdw_transition_midpoint = midpoint_in;
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::setVdwTransitionIntensity(const double intensity_in) {
  vdw_transition_intensity = intensity_in;
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::setEwaldCoefficient(const double ewald_coefficient_in) {
  if (ewald_coefficient_in < minimum_ewald_coefficient ||
      ewald_coefficient_in > maximum_ewald_coefficient) {
    rtErr("An Ewald coefficient of " +
          realToString(ewald_coefficient_in, 9, 5, NumberFormat::STANDARD_REAL) + " suggests an "
          "unreasonable density spreading effect.", "LayeredPotentialMetrics",
          "setEwaldCoefficient");
  }
  ewald_coefficient = ewald_coefficient_in;
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::setUnitCellVectors(const std::vector<double> &box_vectors_in) {
  if (box_vectors_in.size() != 9 || fabs(box_vectors_in[1]) > constants::tiny ||
      fabs(box_vectors_in[2]) > constants::tiny || fabs(box_vectors_in[5]) > constants::tiny) {
    const NumberFormat std_real = NumberFormat::STANDARD_REAL;
    rtErr("The box vectors must be a 3 x 3 matrix and a valid transformation.  Provided: { " +
          realToString(box_vectors_in[0], 9, 4, std_real) + ", " +
          realToString(box_vectors_in[1], 9, 4, std_real) + ", " +
          realToString(box_vectors_in[2], 9, 4, std_real) + " } x { " +
          realToString(box_vectors_in[3], 9, 4, std_real) + ", " +
          realToString(box_vectors_in[4], 9, 4, std_real) + ", " +
          realToString(box_vectors_in[5], 9, 4, std_real) + " } x { " +
          realToString(box_vectors_in[6], 9, 4, std_real) + ", " +
          realToString(box_vectors_in[7], 9, 4, std_real) + ", " +
          realToString(box_vectors_in[8], 9, 4, std_real) + " }.", "LayeredPotentialMetrics",
          "setUnitCellVectors");
  }
  box_vectors = box_vectors_in;

  // The range limit will be set if periodic boundary conditions apply.
  switch (boundaries) {
  case BoundaryCondition::ISOLATED:
    break;
  case BoundaryCondition::PERIODIC:
    range_limit = 0.0;
    for (int i = 0; i < 8; i++) {
      const int ic = i / 4;
      const int ib = (i - (4 * ic)) / 2;
      const int ia = i - (((2 * ic) + ib) * 2);
      const double da = static_cast<double>(ia) * 0.5;
      const double db = static_cast<double>(ib) * 0.5;
      const double dc = static_cast<double>(ic) * 0.5;
      const double dx = (box_vectors[0] * da) + (box_vectors[3] * db) + (box_vectors[6] * dc);
      const double dy =                         (box_vectors[4] * db) + (box_vectors[7] * dc);
      const double dz =                                                 (box_vectors[8] * dc);
      range_limit = std::max(range_limit, sqrt((dx * dx) + (dy * dy) + (dz * dz))); 
    }
    if (range_limit < hail_minimum_range_limit) {
      rtErr("The unit cell appears to be too small, implying a range limit of only " +
            realToString(range_limit, 7, 4, NumberFormat::STANDARD_REAL) + ".",
            "LayeredPotentialMetrics", "setUnitCellVectors");
    }
    break;
  }
  computeLayers();
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::validateLayerIndex(const int layer_index, const char* caller) const {
  if (layer_index < 0 || layer_index >= layer_count) {
    rtErr("Layer index " + std::to_string(layer_index) + " is invalid for an approximation with " +
          std::to_string(layer_count) + " layers.", "LayeredPotentialMetrics", caller);
  }
}

//-------------------------------------------------------------------------------------------------
void LayeredPotentialMetrics::computeLayers() {
  
  // Compute the number of layers given the particle-particle cutoff and the range multiplier
  layer_count = static_cast<int>(ceil(log(range_limit / cutoffs[0]) / log(range_multiplier))) + 1;
  
  // Run the range compounding and compute exponential factors
  cutoffs.resize(layer_count);
  exponent_c.resize(layer_count);
  for (int i = 1; i < layer_count; i++) {
    cutoffs[i] = cutoffs[0] * pow(range_multiplier, i);
    if (i > 1) {
      const double emult = cutoffs[0] / cutoffs[i - 1];
      exponent_c[i].x = exponent_c[1].x * emult;
      exponent_c[i].y = exponent_c[1].y * emult;
      exponent_c[i].z = exponent_c[1].z * emult;
    }
  }
}

} // namespace energy
} // namespace stormm
