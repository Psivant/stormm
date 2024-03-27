// -*-c++-*-
#include "copyright.h"
#include "layered_potential.h"

namespace stormm {
namespace energy {

//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
LayeredPotential<T, T4>::LayeredPotential(const LayeredPotentialMetrics &parameters_in) :
    parameters{parameters_in}
{
  solveLayerCoefficients();
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
const LayeredPotentialMetrics& LayeredPotential<T, T4>::getParameters() const {
  return parameters;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
double4 LayeredPotential<T, T4>::getSmoothingCoefficients(const int layer) const {
  validateLayerIndex(layer, "getSmoothingCoefficients");
  return smoothing_coefficients[layer];
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
T LayeredPotential<T, T4>::getSplinedValue(const int layer, const T r, const T r2) const {

}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
T LayeredPotential<T, T4>::getAnalyticValue(const int layer, const T r, const T r2) const {
  
  // All layers will have zero value beyond the cutoff for which they are applicable.  The
  // parameters object's getCutoff() member function will check the layer index's validity.
  if (r >= parameters.getCutoff(layer)) {
    return 0.0;
  }
  
  // Determine whether the interaction distance lies within the regime at which the potential is
  // smoothed.  This is true if distance is less than the cutoff of the previous layer.  There is
  // no smoothed region for the lowest-level layer.
  const bool smoothed_regime = (layer > 0 && r < parameters.getCutoff(layer - 1));
  const bool tcalc_is_double = (std::type_index(typeid(T)).hash_code() == double_type_index);
  T result;
  const T value_zero = 0.0;
  const T value_one  = 1.0;
  if (smoothed_regime) {
    const double4 smc = smoothing_coefficients[layer];
    const T rh = (layer > 0) ? parameters.getCutoff(layer - 1) : value_zero;
    const T a_e = parameters.getExponentFactor(0, layer);
    const T b_e = parameters.getExponentFactor(1, layer);
    const T c_e = parameters.getExponentFactor(2, layer);
    if (tcalc_is_double) {
      result = (smc.x * exp(a_e * (r - rh))) + (smc.y * exp(b_e * (r - rh))) +
               (smc.z * exp(c_e * (r - rh))) + smc.w;
    }
    else {

      // Perform the calculation in 32-bit floating point numbers with the expf library function.
      // Specify the template type, which will be float4 if this branch is taken, to suppress
      // compiler warnings. 
      const T4 fsmc = { static_cast<T>(smc.x), static_cast<T>(smc.y),
                        static_cast<T>(smc.z), static_cast<T>(smc.w) };
      result = (fsmc.x * expf(a_e * (r - rh))) + (fsmc.y * expf(b_e * (r - rh))) +
               (fsmc.z * expf(c_e * (r - rh))) + fsmc.w;
    }
  }
  else {
    switch (parameters.getForm()) {
    case DecomposablePotential::ELECTROSTATIC:
      result = static_cast<T>(parameters.getCoulombConstant()) / r;
      break;
    case DecomposablePotential::DISPERSION:
      result = (tcalc_is_double) ? -value_one / pow(r, 6.0) : -value_one / powf(r, 6.0f);
      break;
    case DecomposablePotential::ELEC_PME_DIRECT:
      {
        const T use_r2 = (r2 <= 0.0) ? r * r : r2;
        result = elecPMEDirectSpace<T>(parameters.getEwaldCoefficient(),
                                       parameters.getCoulombConstant(), r, use_r2, 0);
      }
      break;
    }
  }

  // Subtract the contributions from further on (higher layers) in the decomposition.  These will
  // all be sampled in their smoothed regimes.
  if (layer < parameters.getLayerCount()) {
    const double4 smc = smoothing_coefficients[layer + 1];
    const T rh = parameters.getCutoff(layer);
    const T a_e = parameters.getExponentFactor(0, layer + 1);
    const T b_e = parameters.getExponentFactor(1, layer + 1);
    const T c_e = parameters.getExponentFactor(2, layer + 1);
    if (tcalc_is_double) {
      result -= (smc.x * exp(a_e * (r - rh))) + (smc.y * exp(b_e * (r - rh))) +
                (smc.z * exp(c_e * (r - rh))) + smc.w;
    }
    else {

      // Perform the calculation in 32-bit floating point numbers with the expf library function
      const T4 fsmc = { static_cast<T>(smc.x), static_cast<T>(smc.y),
                        static_cast<T>(smc.z), static_cast<T>(smc.w) };
      result -= (fsmc.x * expf(a_e * (r - rh))) + (fsmc.y * expf(b_e * (r - rh))) +
                (fsmc.z * expf(c_e * (r - rh))) + fsmc.w;
    }
  }
  
  // The value of a sigmoidal switching function may come into play.  The role of this switching
  // function is the activate the Lennard-Jones interaction by a geometric combining rule as the
  // distance approaches the particle-particle cutoff.  Although all layers of the potential have
  // nonzero character within the particle-particle cutoff (the extent of the lowest layer of the
  // decomposition) they are not subject to the switching function in this regime.  Instead, the
  // complete switching is done within the confines of the particle-particle layer.  This makes
  // interpolation with the other layers simpler, as there is a singlecurve of roughly sigmoidal
  // character (the smoothed Lennard-Jones potential), not two competing forms (the smoothed
  // Lennard-Jones plus the switching function).  Given the need to make the switching for all
  // purposes complete with practically zero derivatives within the space of the particle-particle
  // cutoff, it is reasonable to apply the switching function only to the particle-particle
  // interactions.  In the lowest layer only, the complete geometric Lennard-Jones character will
  // be subtracted off while the model's own Lorentz-Berthelot or NBFix character added back in its
  // place with the complementary switching function.
  if (layer == 0) {
    switch (parameters.getForm()) {
    case DecomposablePotential::ELECTROSTATIC:
    case DecomposablePotential::ELEC_PME_DIRECT:
      break;
    case DecomposablePotential::DISPERSION:
      switch (parameters.getMixingRule()) {
      case VdwCombiningRule::GEOMETRIC:
        break;
      case VdwCombiningRule::LORENTZ_BERTHELOT:
      case VdwCombiningRule::NBFIX:
        if (tcalc_is_double) {
          const double4 sigm = sigmoid(r, parameters.getVdwTransitionMidpoint(),
                                       parameters.getVdwTransitionIntensity());
          result -= sigm.x / pow(r, 6.0);
        }
        else {
          const float4 sigm = sigmoidf(r, parameters.getVdwTransitionMidpoint(),
                                       parameters.getVdwTransitionIntensity());
          result -= sigm.x / powf(r, 6.0f);
        }
        break;
      }
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
T LayeredPotential<T, T4>::getSplinedDerivative(const int layer, const T r, const T r2,
                                                const int order) const {

}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
T LayeredPotential<T, T4>::getAnalyticDerivative(const int layer, const T r, const T r2,
                                                 const int order) const {

  // Redirect immediately to the value calculation for order zero, to avoid performing various
  // input checks twice.
  if (order == 0) {
    return getAnalyticValue(layer, r);
  }
  if (order < 0 || order > 3) {
    rtErr("Derivatves of order " + std::to_string(order) + " are not supported.",
          "LayeredPotential", "getAnalyticDerivative");
  }
  const bool r2_provided = (r2 < 0.0);
  const T use_r2 = (r2_provided) ? r * r : r2;
  
  // All layers will have zero value and zero derivative beyond the cutoff for which they are
  // applicable.
  if (r >= parameters.getCutoff(layer)) {
    return 0.0;
  }
  const bool smoothed_regime = (layer > 0 && r < parameters.getCutoff(layer - 1));
  const bool tcalc_is_double = (std::type_index(typeid(T)).hash_code() == double_type_index);
  T result;
  const T value_zero = 0.0;
  const T value_one  = 1.0;
  const T value_two  = 2.0;
  const T value_six  = 6.0;
  const T value_n42  = -42.0;
  const T value_336  = 336.0;
  if (smoothed_regime) {
    const double4 smc = smoothing_coefficients[layer];
    const T rh = (layer > 0) ? parameters.getCutoff(layer - 1) : value_zero;
    const T a_e = parameters.getExponentFactor(0, layer);
    const T b_e = parameters.getExponentFactor(1, layer);
    const T c_e = parameters.getExponentFactor(2, layer);
    T pa_e, pb_e, pc_e;
    if (order == 1) {
      pa_e = a_e;
      pb_e = b_e;
      pc_e = c_e;
    }
    else if (order == 2) {
      pa_e = a_e * a_e;
      pb_e = b_e * b_e;
      pc_e = c_e * c_e;
    }
    if (tcalc_is_double) {
      if (order == 3) {
        pa_e = pow(a_e, 3.0);
        pb_e = pow(b_e, 3.0);
        pc_e = pow(c_e, 3.0);
      }
      result = (smc.x * pa_e * exp(a_e * (r - rh))) + (smc.y * pb_e * exp(b_e * (r - rh))) +
               (smc.z * pc_e * exp(c_e * (r - rh)));
    }
    else {
      if (order == 3) {
        pa_e = powf(a_e, 3.0f);
        pb_e = powf(b_e, 3.0f);
        pc_e = powf(c_e, 3.0f);
      }
      const T4 fsmc = { static_cast<T>(smc.x), static_cast<T>(smc.y),
                        static_cast<T>(smc.z), static_cast<T>(smc.w) };
      result = (fsmc.x * pa_e * expf(a_e * (r - rh))) + (fsmc.y * pb_e * expf(b_e * (r - rh))) +
               (fsmc.z * pc_e * expf(c_e * (r - rh)));
    }
  }
  else {
    switch (parameters.getForm()) {
    case DecomposablePotential::ELECTROSTATIC:
      {
        const T kcoul = parameters.getCoulombConstant();
        switch (order) {
        case 0:

          // This case was diverted to the function for computing function values above
          break;
        case 1:
          result = -kcoul / use_r2;
          break;
        case 2:
          if (r2_provided) {
            result = value_two * kcoul / (r2 * r);
          }
          else {
            if (tcalc_is_double) {
              result = value_two * kcoul / pow(r, 3.0);
            }
            else {
              result = value_two * kcoul / powf(r, 3.0f);
            }
          }
          break;
        case 3:
          if (r2_provided) {
            result = -value_six * kcoul / (r2 * r2);
          }
          else {
            if (tcalc_is_double) {
              result = -value_six * kcoul / pow(r, 4.0);
            }
            else {
              result = -value_six * kcoul / powf(r, 4.0f);
            }
          }
          break;
        }
      }
      break;
    case DecomposablePotential::DISPERSION:
      switch (order) {
      case 0:
        break;
      case 1:
        result = (tcalc_is_double) ? value_six / pow(r, 7.0) : value_six / powf(r, 7.0f);
        break;
      case 2:
        result = (tcalc_is_double) ? value_n42 / pow(r, 8.0) : value_n42 / powf(r, 8.0f);
        break;
      case 3:
        result = (tcalc_is_double) ? value_336 / pow(r, 9.0) : value_336 / powf(r, 9.0f);
        break;
      }
      break;
    case DecomposablePotential::ELEC_PME_DIRECT:
      {
        const T use_r2 = (r2 <= 0.0) ? r * r : r2;
        result = elecPMEDirectSpace<T>(parameters.getEwaldCoefficient(),
                                       parameters.getCoulombConstant(), r, use_r2, order);
      }
      break;
    }
  }

  // Subtract the contributions from further on (higher layers) in the decomposition.  These will
  // all be sampled in their smoothed regimes.
  if (layer < parameters.getLayerCount()) {
    const double4 smc = smoothing_coefficients[layer + 1];
    const T rh = parameters.getCutoff(layer);
    const T a_e = parameters.getExponentFactor(0, layer + 1);
    const T b_e = parameters.getExponentFactor(1, layer + 1);
    const T c_e = parameters.getExponentFactor(2, layer + 1);
    T pa_e, pb_e, pc_e;
    if (order == 1) {
      pa_e = a_e;
      pb_e = b_e;
      pc_e = c_e;
    }
    else if (order == 2) {
      pa_e = a_e * a_e;
      pb_e = b_e * b_e;
      pc_e = c_e * c_e;
    }
    if (tcalc_is_double) {
      if (order == 3) {
        pa_e = pow(a_e, 3.0);
        pb_e = pow(b_e, 3.0);
        pc_e = pow(c_e, 3.0);
      }
      result -= (smc.x * pa_e * exp(a_e * (r - rh))) + (smc.y * pb_e * exp(b_e * (r - rh))) +
                (smc.z * pc_e * exp(c_e * (r - rh)));
    }
    else {
      if (order == 3) {
        pa_e = powf(a_e, 3.0f);
        pb_e = powf(b_e, 3.0f);
        pc_e = powf(c_e, 3.0f);
      }
      const T4 fsmc = { static_cast<T>(smc.x), static_cast<T>(smc.y),
                        static_cast<T>(smc.z), static_cast<T>(smc.w) };
      result -= (fsmc.x * pa_e * expf(a_e * (r - rh))) + (fsmc.y * pb_e * expf(b_e * (r - rh))) +
                (fsmc.z * pc_e * expf(c_e * (r - rh)));
    }
  }

  // As in the case of the potential itself, handle the effects of the subtracted geometric rule
  // Lennard-Jones interaction on the derivative.
  if (layer == 0) {
    switch (parameters.getForm()) {
    case DecomposablePotential::ELECTROSTATIC:
    case DecomposablePotential::ELEC_PME_DIRECT:
      break;
    case DecomposablePotential::DISPERSION:
      if (tcalc_is_double) {
        const double4 sigm = sigmoid(r, parameters.getVdwTransitionMidpoint(),
                                     parameters.getVdwTransitionIntensity());
        const T f_sw = value_one - sigm.x;
        const T df_sw = -sigm.y;
        const T d2f_sw = -sigm.z;
        const T d3f_sw = -sigm.w;
        const T u   = -value_one / pow(r, 6.0);
        const T du  =  value_six / pow(r, 7.0);
        switch (order) {
        case 1:
          result -= (u * df_sw) + (du * f_sw); 
          break;
        case 2:
          {
            const T d2u = value_n42 / pow(r, 8.0);
            result -= (u * d2f_sw) + (value_two * (du * df_sw)) + (d2u * f_sw);
          }
          break;
        case 3:
          {
            const T d2u =   value_n42 / pow(r, 8.0);
            const T d3u = value_336 / pow(r, 9.0);
            result -= (u * d3f_sw) + (3.0 * (du * d2f_sw)) + (3.0 * (d2u * df_sw)) + (d3u * f_sw);
          }
          break;
        }
      }
      else {
        const float4 sigm = sigmoidf(r, parameters.getVdwTransitionMidpoint(),
                                     parameters.getVdwTransitionIntensity());
        const T f_sw = value_one - sigm.x;
        const T df_sw = -sigm.y;
        const T d2f_sw = -sigm.z;
        const T d3f_sw = -sigm.w;
        const T u   = -value_one / powf(r, 6.0f);
        const T du  =  value_six / powf(r, 7.0f);
        switch (order) {
        case 1:
          result -= (u * df_sw) + (du * f_sw); 
          break;
        case 2:
          {
            const T d2u = value_n42 / powf(r, 8.0f);
            result -= (u * d2f_sw) + (value_two * (du * df_sw)) + (d2u * f_sw);
          }
          break;
        case 3:
          {
            const T d2u = value_n42 / powf(r, 8.0f);
            const T d3u = value_336 / powf(r, 9.0f);
            result -= (u * d3f_sw) + (3.0f * (du * d2f_sw)) + (3.0 * (d2u * df_sw)) + (d3u * f_sw);
          }
          break;
        }
      }
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
const NrgLayerKit<T4> LayeredPotential<T, T4>::data(const int layer_index,
                                                    const HybridTargetLevel tier) const {
  
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
void LayeredPotential<T, T4>::validateLayerIndex(const int layer_index, const char* caller) const {
  if (layer_index < 0 || layer_index >= parameters.getLayerCount()) {
    rtErr("Layer index " + std::to_string(layer_index) + " is invalid for an approximation with " +
          std::to_string(parameters.getLayerCount()) + " layers.", "LayeredPotential", caller);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename T4>
void LayeredPotential<T, T4>::solveLayerCoefficients() {
  const int nlayer = parameters.getLayerCount();
  smoothing_coefficients.resize(nlayer);
  smoothing_coefficients[0] = { 0.0, 0.0, 0.0, 0.0 };
  std::vector<double> fmat(16), tvec(4), xvec(4);
  for (int i = 1; i < nlayer; i++) {

    // Set up the matrix to solve for coefficients of the exponential series.  The smoothed
    // functional form for layer k out of a total of n layers has the form:
    //
    //   F[k](r) = A exp(a_e(k) (r - Rh)) + B exp(b_e(k) (r - Rh)) + C exp(c_e(k) (r - Rh)) + D,
    //             r <= Rh
    //           = (Original potential), r > Rh
    //
    // where Rh is the "handoff point", equal to the cutoff of the previous layer.  Coefficients A,
    // B, C, and D are fitted to ensure that the smoothed form is of class C^3.  The total
    // potential at layer k has the form F[k](r) - F[k+1](r) - F[k+2](r) - ... - F[m](r), where m
    // is the total number of layers.  The coefficients a_e(k), b_e(k), and c_e(k) are derived from
    // three baseline coefficients set as parameters of the model, scaled for the kth layer by the
    // ratio of the cutoff in the base layer to the cutoff in the kth layer.
    const double a_e = parameters.getExponentFactor(0, i);
    const double b_e = parameters.getExponentFactor(1, i);
    const double c_e = parameters.getExponentFactor(2, i);
    fmat[ 0] = 1.0;
    fmat[ 4] = 1.0;
    fmat[ 8] = 1.0;
    fmat[12] = 1.0;
    fmat[ 1] = a_e;
    fmat[ 5] = b_e;
    fmat[ 9] = c_e;
    fmat[13] = 0.0;
    fmat[ 2] = a_e * a_e;
    fmat[ 6] = b_e * b_e;
    fmat[10] = c_e * c_e;
    fmat[14] = 0.0;

    // Use pow() to compute integer powers of numbers rather than mult when accuracy is more
    // important than speed.
    fmat[ 3] = pow(a_e, 3.0);
    fmat[ 7] = pow(b_e, 3.0);
    fmat[11] = pow(c_e, 3.0);
    fmat[15] = 0.0;

    // The target potential takes one of the prescribed forms.
    const double rh = parameters.getCutoff(i - 1);
    switch (parameters.getForm()) {
    case DecomposablePotential::ELECTROSTATIC:
      {
        const double kcoul = parameters.getCoulombConstant();
        tvec[0] = kcoul / rh;
        tvec[1] = -kcoul / (rh * rh);
        tvec[2] = 2.0 * kcoul / pow(rh, 3.0);
        tvec[3] = -6.0 * kcoul / pow(rh, 4.0);
      }
      break;
    case DecomposablePotential::DISPERSION:

      // A strict geometric combining rule is the simplest to formulate in terms of the long-ranged
      // potential.  There is no need to switch it off at short range for generating the mesh-based
      // layers of the potential.
      tvec[0] =   -1.0 / pow(rh, 6.0);
      tvec[1] =    6.0 / pow(rh, 7.0);
      tvec[2] =  -42.0 / pow(rh, 8.0);
      tvec[3] =  336.0 / pow(rh, 9.0);
      break;
    case DecomposablePotential::ELEC_PME_DIRECT:
      for (int j = 0; j < 4; j++) {
        tvec[j] = elecPMEDirectSpace<T>(parameters.getEwaldCoefficient(),
                                        parameters.getCoulombConstant(), rh, rh * rh, j);
      }
      break;
    }
    
    // Invert the matrix and solve for the coefficients at this layer.  The matrix must be
    // refreshed with each layer because the exponent scaling factors diminish as the cutoff gets
    // longer.
    qrSolver(&fmat, &xvec, &tvec);
    smoothing_coefficients[i] = { xvec[0], xvec[1], xvec[2], xvec[3] };
  }
}

} // namespace energy
} // namespace stormm
